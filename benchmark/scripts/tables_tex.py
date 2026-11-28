import pandas as pd, numpy as np
S="PATH/TO/WORKING_DIR/"  # set to your working directory (see benchmark/README.md)
T=pd.read_csv(S+"final/cross_main.csv"); tier=pd.read_csv(S+"final/cross_tier.csv"); fl=pd.read_csv(S+"final/cross_flagged.csv"); cond=pd.read_csv(S+"final/cross_cond.csv"); pr=pd.read_csv(S+"final/cross_paired.csv"); cl=pd.read_csv(S+"final/ref_clash.csv")
A=pd.concat([pd.read_csv(S+"esr1_out/ESR1_all_strategies_scored.csv").assign(target="ESR1"),pd.read_csv(S+"mp_out/MAPK1_PPARG_all_scored.csv")],ignore_index=True)
M=pd.read_csv("benchmark/data/cross_docking/inputs/manifest.csv")
ST=["cocryst","fpocket","p2rank","blind"]; NM=dict(cocryst="Co-crystallized",fpocket="fpocket",p2rank="P2Rank",blind="Blind"); TG=["ESR1","MAPK1","PPARG"]
f1=lambda x:f"{x:.1f}"
def ci(r,m): return f"{r[m]:.1f} ({r[m+'_lo']:.1f}--{r[m+'_hi']:.1f})"
# ---- main-text table
def main_table():
    L=[]
    for t in TG+["POOLED_pair","POOLED_eq"]:
        for s in ST:
            r=T[(T.target==t)&(T.strategy==s)].iloc[0]
            if t in TG:
                lab=t if s=="cocryst" else ""; ck=f"{int(r.cog_t1)}/{int(r.n_cog)}"; cb=f"{int(r.cog_bp)}/{int(r.n_cog)}"
            else:
                lab={"POOLED_pair":"All (pairs)","POOLED_eq":"All (targets)"}[t] if s=="cocryst" else ""
                sub=T[(T.target.isin(TG))&(T.strategy==s)]; ck=f"{int(sub.cog_t1.sum())}/{int(sub.n_cog.sum())}"; cb=f"{int(sub.cog_bp.sum())}/{int(sub.n_cog.sum())}"
            L.append(f"    {lab} & {NM[s]} & {ck} & {cb} & {ci(r,'t1')} & {ci(r,'bp')} & {ci(r,'site')} \\\\")
        L.append("    \\addlinespace")
    return "\n".join(L[:-1])
# ---- SI tables
def si_tier():
    L=[]
    for t in TG:
        for s in ST:
            k=T[(T.target==t)&(T.strategy==s)].iloc[0]
            L.append(f"{t} & {NM[s]} & cognate & {int(k.n_cog)} & {int(k.cog_t1)}/{int(k.n_cog)} & {int(k.cog_bp)}/{int(k.n_cog)} & {int(k.cog_site)}/{int(k.n_cog)} \\\\")
            for ti in ["easy","medium","hard"]:
                g=tier[(tier.target==t)&(tier.strategy==s)&(tier.tier==ti)]
                if len(g)==0: continue
                g=g.iloc[0]; n=int(g.n); L.append(f" &  & {ti} & {n} & {g.t1*n:.0f}/{n} ({g.t1*100:.1f}\\%) & {g.bp*n:.0f}/{n} ({g.bp*100:.1f}\\%) & {g.site*n:.0f}/{n} ({g.site*100:.1f}\\%) \\\\")
            L.append("\\addlinespace")
    return "\n".join(L)
def si_flag():
    L=[]
    for t in ["MAPK1","PPARG"]:
        for s in ST:
            g=fl[(fl.target==t)&(fl.strategy==s)].iloc[0]; L.append(f"{t if s=='cocryst' else ''} & {NM[s]} & {int(g.n)} & {int(g.t1)} & {int(g.bp)} & {int(g.site)} \\\\")
        L.append("\\addlinespace")
    return "\n".join(L)
def si_cond():
    L=[]
    for t in TG:
        for s in ST:
            g=cond[(cond.target==t)&(cond.strategy==s)].iloc[0]; L.append(f"{t if s=='cocryst' else ''} & {NM[s]} & {int(g.n)} & {int(g.n_site)} & {g.t1_given_site:.1f} & {g.gap:.1f} \\\\")
        L.append("\\addlinespace")
    return "\n".join(L)
def si_clash():
    c=A[(A.strategy=="cocryst")&A.primary_analysis&~A.cognate].merge(cl[["pair_id","clash"]],on="pair_id"); L=[]
    for t in TG:
        for flag,nm in [(False,"no reference atom within 2.5 \\AA"),(True,"$\\geq$1 reference atom within 2.5 \\AA")]:
            g=c[(c.target==t)&(c.clash==flag)]; L.append(f"{t if not flag else ''} & {nm} & {len(g)} & {int(g.t1.sum())} ({g.t1.mean()*100:.1f}) & {int(g.bp.sum())} ({g.bp.mean()*100:.1f}) & {int(g.site.sum())} ({g.site.mean()*100:.1f}) \\\\")
        L.append("\\addlinespace")
    return "\n".join(L)
def si_paired():
    L=[]; lab={"ESR1":"ESR1","MAPK1":"MAPK1","PPARG":"PPARG","POOLED_pair":"All (pairs)","POOLED_eq":"All (targets)"}
    pairs=[("cocryst","blind"),("p2rank","blind"),("fpocket","blind"),("cocryst","fpocket"),("cocryst","p2rank"),("p2rank","fpocket")]
    for t in ["ESR1","MAPK1","PPARG","POOLED_pair","POOLED_eq"]:
        for i,(a,b) in enumerate(pairs):
            def g(metric):
                x=pr[(pr.target==t)&(pr.metric==metric)&(((pr.a==a)&(pr.b==b))|((pr.a==b)&(pr.b==a)))].iloc[0]
                sign=1 if x.a==a else -1; d=sign*x["diff"]; lo,hi=(x.lo,x.hi) if sign==1 else (-x.hi,-x.lo); return f"{d:+.1f} ({lo:+.1f} to {hi:+.1f})"
            L.append(f"{lab[t] if i==0 else ''} & {NM[a]} $-$ {NM[b]} & {g('t1')} & {g('site')} \\\\")
        L.append("\\addlinespace")
    return "\n".join(L)
def si_pairs():
    L=[]; fmt=lambda x:"--" if pd.isna(x) else f"{x:.2f}"
    for t in TG:
        mt=M[M.target==t]
        for r in mt.itertuples():
            if r.target=="ESR1" and r.ligand=="6chw": continue
            row=[f"{t} {r.receptor.upper()}/{r.ligand.upper()}", "cog" if r.cognate else r.tier, f"{r.pocket_CA_rmsd:.2f}", "F" if r.noneq_site else ""]
            for s in ST:
                g=A[(A.pair_id==r.pair_id)&(A.strategy==s)]
                row+= [fmt(g.top1_rmsd.values[0]),fmt(g.best_rmsd.values[0])] if len(g) else ["--","--"]
            L.append(" & ".join(row)+" \\\\")
    return "\n".join(L)
def si_lig():
    L=[]
    Lg=M[["target","ligand","ligand_load"]].drop_duplicates()
    from rdkit import Chem, RDLogger; RDLogger.DisableLog("rdApp.*")
    for r in Lg.sort_values(["target","ligand"]).itertuples():
        m=Chem.MolFromMolFile(S+f"cd_pkg/crossdock_out/ligands/{r.target}_{r.ligand}.sdf",removeHs=True,sanitize=False); ha=sum(1 for a in m.GetAtoms() if a.GetAtomicNum()>1)
        ld=r.ligand_load; ch=""; kind="read directly (mol2)"
        if ld!="mol2":
            kind="rebuilt" if ld.startswith("rebuilt") else "repaired"; 
            if "halogen" in ld: kind="repaired + halogen re-attached"
            ch=ld.split("charge=")[1].rstrip(")")
        note="covalent, not docked" if (r.target=="ESR1" and r.ligand=="6chw") else ""
        L.append(f"{r.target} & {r.ligand.upper()} & {kind} & {ch if ch else '--'} & {ha} & {note} \\\\")
    return "\n".join(L)
if __name__=="__main__":
    for f in [main_table,si_tier,si_flag,si_cond,si_clash,si_paired,si_lig]: open(S+f"final/{f.__name__}.tex","w").write(f())
    print(open(S+"final/main_table.tex").read()); print(len(si_pairs().splitlines()))
    open(S+"final/si_pairs.tex","w").write(si_pairs())
