"""Counterparts of Tables 2, 3, 5 and the enclosure analysis at every exhaustiveness (8, 16, 32). Writes astex_extra.json."""
import pandas as pd,numpy as np,itertools,json,warnings,statsmodels.api as sm,statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")
S="/tmp/claude-0/-home-claude/fa06f76e-b715-54d8-a3bb-7849315bf476/scratchpad/final/"
X=pd.read_csv(S+"astex_merged.csv"); G=pd.read_csv(S+"astex_geometry.csv"); G=G[(G.ex==16)&(G.seed==42)][["strategy","pdb","vol","full"]]
ST=["cocryst","fpocket","p2rank","blind"]; rng=np.random.default_rng(1); out={}
def mcn(a,b):
    n10=int(((a==1)&(b==0)).sum()); n01=int(((a==0)&(b==1)).sum()); n=n10+n01
    return 1.0 if n==0 else min(1,2*stats.binom.cdf(min(n10,n01),n,.5))
for ex in (8,16,32):
    D=X[X.ex==ex].copy(); D["y1"]=(D.sp_top1<2).astype(int); D["yb"]=(D.sp_best<2).astype(int); pdbs=sorted(D.pdb.unique()); r={}
    for ep,y,col in (("top1","y1","sp_top1"),("best","yb","sp_best")):
        rows={}
        for s in ST:
            g=D[D.strategy==s]; per=[int(g[g.seed==sd][y].sum()) for sd in (1,42,123)]; rate=np.array(per)/85*100
            frac=g.pivot(index="pdb",columns="seed",values=y).reindex(pdbs).values.mean(1)
            bs=[frac[rng.integers(0,85,85)].mean()*100 for _ in range(10000)]
            rows[s]=dict(per=per,mean=rate.mean(),sd=rate.std(ddof=1),lo=np.percentile(bs,2.5),hi=np.percentile(bs,97.5))
        cons={s:(D[D.strategy==s].pivot(index="pdb",columns="seed",values=y).reindex(pdbs).sum(1)>=2).astype(int).values for s in ST}
        pairs=list(itertools.combinations(ST,2)); p=[mcn(cons[a],cons[b]) for a,b in pairs]; ph=multipletests(p,method="holm")[1]
        # Cochran Q
        M=np.column_stack([cons[s] for s in ST]); k=4; Tj=M.sum(0); Ri=M.sum(1); N=M.sum()
        Q=(k-1)*(k*(Tj**2).sum()-N**2)/(k*N-(Ri**2).sum()); pQ=stats.chi2.sf(Q,k-1)
        # RMSD subset
        piv={s:D[D.strategy==s].pivot(index="pdb",columns="seed",values=col).reindex(pdbs) for s in ST}
        ok=np.all([piv[s].notna().all(1).values for s in ST],axis=0); med={s:piv[s].median(1).values[ok] for s in ST}
        fr=stats.friedmanchisquare(*[med[s] for s in ST]); pw=[stats.wilcoxon(med[a],med[b]).pvalue for a,b in pairs]; pwh=multipletests(pw,method="holm")[1]
        r[ep]=dict(rows=rows,pairs=[f"{a}|{b}" for a,b in pairs],mcn=list(ph),Q=Q,pQ=pQ,n=int(ok.sum()),fr=fr.statistic,pfr=fr.pvalue,wil=list(pwh),
                   med={s:[float(np.median(med[s])),*np.percentile(med[s],[25,75])] for s in ST})
    # enclosure and GEE (top-1, seed mean, unresolved = failure)
    Z=D.merge(G,on=["strategy","pdb"]); Z["lv"]=np.log10(Z.vol); Z["strategy"]=pd.Categorical(Z.strategy,ST)
    enc={}
    for s in ST:
        z=Z[Z.strategy==s]; f=z.groupby("pdb").agg(y=("y1","mean"),full=("full","first"))
        enc[s]=dict(n_in=int(f.full.sum()),n_out=int((~f.full).sum()),y_in=float(f[f.full].y.mean()*100),y_out=float(f[~f.full].y.mean()*100) if (~f.full).any() else None)
    r["enc"]=enc; gee={}
    for nm,fm in (("unadj","y1~C(strategy)"),("adjvol","y1~C(strategy)+lv")):
        m=smf.gee(fm,"pdb",data=Z,family=sm.families.Binomial(),cov_struct=sm.cov_struct.Exchangeable()).fit(); ci=m.conf_int()
        gee[nm]={k.split("[T.")[1][:-1]:[float(np.exp(m.params[k])),float(np.exp(ci.loc[k,0])),float(np.exp(ci.loc[k,1]))] for k in m.params.index if "strategy" in k}
        if "lv" in m.params: gee[nm]["lv"]=[float(np.exp(m.params["lv"])),float(np.exp(ci.loc["lv",0])),float(np.exp(ci.loc["lv",1]))]
    r["gee"]=gee; out[ex]=r
json.dump(out,open(S+"astex_extra.json","w"),default=float,indent=1)
r=out[16]; print({s:(round(v['mean'],1),round(v['lo'],1),round(v['hi'],1)) for s,v in r['top1']['rows'].items()}); print([round(x,4) for x in r['top1']['mcn']], r['top1']['Q'],r['top1']['pQ'],r['top1']['n'],r['top1']['fr']); print(r['enc'],r['gee'])
