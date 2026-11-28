import re,os,glob,sys,numpy as np,pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
S="PATH/TO/WORKING_DIR/"  # set to your working directory (see benchmark/README.md)
ZIP={8:"exhaust_8_with_3_seeds",16:"exhaust_16_with_3_seeds",32:"exhaus_32_with_3_seeds"}
SEEDS={1:"seed_1",42:"seed_42",123:"seed123"}
STR={"cocryst":"co-crystallized","fpocket":"fpocket","p2rank":"p2rank","blind":"blind"}
def refat(pid,R):
    o=[]
    for ln in open(R+f"raw_ligands/{pid}_ligand.pdb"):
        if ln.startswith(("ATOM","HETATM")):
            el=ln[76:78].strip().upper() or re.sub(r"\d","",ln[12:16]).strip().upper()
            if el in("H","D"): continue
            o.append((el,float(ln[30:38]),float(ln[38:46]),float(ln[46:54])))
    return o
def models(p):
    ms=[];cur=None
    for ln in open(p,errors="ignore"):
        if ln.startswith("MODEL"): cur=[]
        elif ln.startswith(("ATOM","HETATM")) and cur is not None:
            el=ln[77:79].strip().upper()
            if el in("H","HD"): continue
            el={"A":"C","OA":"O","NA":"N","SA":"S"}.get(el,el); cur.append((el,float(ln[30:38]),float(ln[38:46]),float(ln[46:54])))
        elif ln.startswith("ENDMDL"): ms.append(cur); cur=None
    return ms
def pdir(R,seed,strat):
    sd=[n for n in os.listdir(R) if n.startswith("seed") and re.sub(r"\D","",n)==str(seed)][0]
    for n in os.listdir(R+sd):
        if n.lower().replace("-","")==STR[strat].replace("-",""): return R+sd+"/"+n+"/"
    raise SystemExit(f"no {strat} in {R+sd}")
rows=[];geo=[]
for ex,z in ZIP.items():
    R=S+f"{z}/{z}/"
    refs={}
    for seed in SEEDS:
        for st in STR:
            D=pdir(R,seed,st)
            for f in sorted(glob.glob(D+"docking/*.pdbqt")):
                pid=os.path.basename(f)[:4]
                if pid not in refs: refs[pid]=refat(pid,R)
                r=refs[pid]; re_=np.array([a[0] for a in r]); rx=np.array([a[1:] for a in r])
                ms=models(f); idn=[];hun=[]
                for m in ms:
                    pe=np.array([a[0] for a in m]); px=np.array([a[1:] for a in m])
                    if len(pe)!=len(re_): idn.append(np.nan);hun.append(np.nan);continue
                    idn.append(float(np.sqrt(((px-rx)**2).sum(1).mean())))
                    C=cdist(px,rx)**2; C[pe[:,None]!=re_[None,:]]=1e6
                    i,j=linear_sum_assignment(C); hun.append(float(np.sqrt(C[i,j].mean())) if C[i,j].max()<1e5 else np.nan)
                rows.append(dict(ex=ex,seed=seed,strategy=st,pdb=pid,n_models=len(ms),id_top1=idn[0],id_best=np.nanmin(idn) if not np.all(np.isnan(idn)) else np.nan,nn_top1=hun[0],nn_best=np.nanmin(hun) if not np.all(np.isnan(hun)) else np.nan))
                if ex==16 and seed==42 or True:
                    g=dict(l.split(" = ") for l in open(D+f"grid/{pid}/grid_config.txt").read().split("\n") if " = " in l)
                    c=np.array([float(g[f"center_{a}"]) for a in "xyz"]); sz=np.array([float(g[f"size_{a}"]) for a in "xyz"])
                    inside=(np.abs(rx-c)<=sz/2).all(1).mean()
                    geo.append(dict(ex=ex,seed=seed,strategy=st,pdb=pid,vol=sz.prod(),sx=sz[0],sy=sz[1],sz=sz[2],offset=float(np.linalg.norm(c-rx.mean(0))),inside=inside*100,full=bool(inside==1)))
    print(ex,"done",flush=True)
pd.DataFrame(rows).to_csv(S+"../final/astex_naive_scored.csv",index=False); pd.DataFrame(geo).to_csv(S+"../final/astex_geometry.csv",index=False)
