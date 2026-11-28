import pandas as pd, numpy as np, json, itertools
S="PATH/TO/WORKING_DIR/"  # set to your working directory (see benchmark/README.md)
E=pd.read_csv(S+"esr1_out/ESR1_all_strategies_scored.csv"); E["target"]="ESR1"
A=pd.concat([E,pd.read_csv(S+"mp_out/MAPK1_PPARG_all_scored.csv")],ignore_index=True)
ST=["cocryst","fpocket","p2rank","blind"]; TG=["ESR1","MAPK1","PPARG"]; MET=["t1","bp","site"]; B=10000
P=A[A.primary_analysis]; C=P[~P.cognate]; K=P[P.cognate]
rng=np.random.default_rng(20261125)
mats={}; idx={}
for t in TG:
    c=C[C.target==t]; R=sorted(c.receptor.unique()); L=sorted(c.ligand.unique())
    arr=np.full((len(ST)*len(MET),len(R),len(L)),np.nan)
    for si,s in enumerate(ST):
        g=c[c.strategy==s]; ri=g.receptor.map({r:i for i,r in enumerate(R)}).values; li=g.ligand.map({l:i for i,l in enumerate(L)}).values
        for mi,m in enumerate(MET): arr[si*3+mi,ri,li]=g[m].astype(float).values
    mats[t]=arr; idx[t]=(rng.integers(0,len(R),(B,len(R))),rng.integers(0,len(L),(B,len(L))))
sums={t:np.zeros((B,12)) for t in TG}; cnts={t:np.zeros((B,12)) for t in TG}
for t in TG:
    arr=mats[t]; RI,LI=idx[t]
    for b in range(B):
        v=arr[:,RI[b][:,None],LI[b][None,:]]; ok=~np.isnan(v); sums[t][b]=np.nansum(v,axis=(1,2)); cnts[t][b]=ok.sum(axis=(1,2))
boot={t:sums[t]/cnts[t] for t in TG}
boot["POOLED_pair"]=sum(sums[t] for t in TG)/sum(cnts[t] for t in TG); boot["POOLED_eq"]=np.mean([boot[t] for t in TG],axis=0)
def ci(x): return np.percentile(x,[2.5,97.5])*100
rows=[]
for t in TG+["POOLED_pair","POOLED_eq"]:
    for si,s in enumerate(ST):
        r=dict(target=t,strategy=s)
        for mi,m in enumerate(MET):
            x=boot[t][:,si*3+mi]; lo,hi=ci(x)
            if t in TG: pt=C[(C.target==t)&(C.strategy==s)][m].mean()*100
            elif t=="POOLED_pair": pt=C[C.strategy==s][m].mean()*100
            else: pt=np.mean([C[(C.target==tt)&(C.strategy==s)][m].mean() for tt in TG])*100
            r.update({m:pt,m+"_lo":lo,m+"_hi":hi})
        if t in TG:
            k=K[(K.target==t)&(K.strategy==s)]; r.update(n_cross=int(C[(C.target==t)&(C.strategy==s)].shape[0]),n_cog=len(k),cog_t1=int(k.t1.sum()),cog_bp=int(k.bp.sum()),cog_site=int(k.site.sum()),med_rmsd=C[(C.target==t)&(C.strategy==s)].top1_rmsd.median())
        rows.append(r)
T=pd.DataFrame(rows); T.to_csv(S+"final/cross_main.csv",index=False)
dr=[]
for t in TG+["POOLED_pair","POOLED_eq"]:
    for a,b in itertools.combinations(ST,2):
        for mi,m in enumerate(MET):
            x=boot[t][:,ST.index(a)*3+mi]-boot[t][:,ST.index(b)*3+mi]; lo,hi=ci(x)
            pt=(T[(T.target==t)&(T.strategy==a)][m].values[0]-T[(T.target==t)&(T.strategy==b)][m].values[0])
            dr.append(dict(target=t,a=a,b=b,metric=m,diff=pt,lo=lo,hi=hi))
pd.DataFrame(dr).to_csv(S+"final/cross_paired.csv",index=False)
# tier, flagged, conditional
tier=C.groupby(["target","tier","strategy"]).agg(n=("t1","size"),t1=("t1","mean"),bp=("bp","mean"),site=("site","mean")).reset_index(); tier.to_csv(S+"final/cross_tier.csv",index=False)
Fl=A[A.noneq_site&~A.cognate]; fl=Fl.groupby(["target","strategy"]).agg(n=("t1","size"),t1=("t1","sum"),bp=("bp","sum"),site=("site","sum")).reset_index(); fl.to_csv(S+"final/cross_flagged.csv",index=False)
cond=[]
for t in TG:
    for s in ST:
        g=C[(C.target==t)&(C.strategy==s)]; cond.append(dict(target=t,strategy=s,n=len(g),n_site=int(g.site.sum()),t1_given_site=g[g.site].t1.mean()*100,t1_not_site=g[~g.site].t1.mean()*100 if (~g.site).any() else np.nan,gap=(g.bp.mean()-g.t1.mean())*100))
pd.DataFrame(cond).to_csv(S+"final/cross_cond.csv",index=False)
pd.set_option("display.width",250); print(T.round(1).to_string()); print(pd.DataFrame(dr).query("metric in ['site','t1'] and target in ['ESR1','MAPK1','PPARG','POOLED_pair']").round(1).to_string())
print(pd.DataFrame(cond).round(1).to_string()); print(fl.to_string())
