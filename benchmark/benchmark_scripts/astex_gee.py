import pandas as pd,numpy as np,statsmodels.api as sm,statsmodels.formula.api as smf,warnings; warnings.filterwarnings("ignore")
S="/tmp/claude-0/-home-claude/fa06f76e-b715-54d8-a3bb-7849315bf476/scratchpad/"
G=pd.read_csv(S+"final/astex_geometry.csv"); G=G[(G.ex==16)&(G.seed==42)]
X=pd.read_csv(S+"final/astex_merged.csv"); X=X[X.ex==16].merge(G[["strategy","pdb","vol","offset"]],on=["strategy","pdb"])
R=S+"astex/exhaust_16_with_3_seeds/exhaust_16_with_3_seeds/"
ha={p:sum(1 for ln in open(R+f"raw_ligands/{p}_ligand.pdb") if ln.startswith(("ATOM","HETATM")) and (ln[76:78].strip().upper() not in("H","D"))) for p in X.pdb.unique()}
X["ha"]=X.pdb.map(ha); X["y"]=(X.sp_top1<2).astype(int); X["lv"]=np.log10(X.vol)
X["strategy"]=pd.Categorical(X.strategy,["cocryst","fpocket","p2rank","blind"])
res={}
for nm,f in [("unadj","y~C(strategy)"),("adj_vol","y~C(strategy)+lv"),("adj_vol_ha","y~C(strategy)+lv+ha")]:
    m=smf.gee(f,"pdb",data=X,family=sm.families.Binomial(),cov_struct=sm.cov_struct.Exchangeable()).fit(); ci=m.conf_int()
    res[nm]={k.split("[T.")[1][:-1]:(np.exp(m.params[k]),np.exp(ci.loc[k,0]),np.exp(ci.loc[k,1]),m.pvalues[k]) for k in m.params.index if "strategy" in k}
    if "lv" in m.params: res[nm]["lv(per log10 unit)"]=(np.exp(m.params["lv"]),np.exp(ci.loc["lv",0]),np.exp(ci.loc["lv",1]),m.pvalues["lv"])
for k,v in res.items():
    print(k); [print("  ",a,[round(float(x),3) for x in b]) for a,b in v.items()]
print(X.groupby("strategy").lv.agg(["mean","std","min","max"]).round(2)); print("corr lv vs blind",np.corrcoef(X.lv,(X.strategy=="blind"))[0,1])
# within-strategy relation of volume to success (cocryst, fpocket)
for st in ["cocryst","fpocket","blind"]:
    g=X[X.strategy==st]; print(st,"corr(log vol,success)=%.2f"%np.corrcoef(g.lv,g.y)[0,1],"corr(ha,vol)=%.2f"%np.corrcoef(g.ha,g.lv)[0,1])
import pickle; pickle.dump(res,open(S+"final/gee_res.pkl","wb"))
