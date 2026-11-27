import pandas as pd,numpy as np
S="/tmp/claude-0/-home-claude/fa06f76e-b715-54d8-a3bb-7849315bf476/scratchpad/"
G=pd.read_csv(S+"final/astex_geometry.csv"); G=G[(G.ex==16)&(G.seed==42)]
st=["cocryst","fpocket","p2rank","blind"]; P={s:G[G.strategy==s].set_index("pdb") for s in st}
L=[]
for pdb in sorted(P["cocryst"].index):
    r=[pdb]+[f"{P[s].loc[pdb,'vol']:,.0f}".replace(",","{,}") for s in st]+[f"{P[s].loc[pdb,'offset']:.1f}" for s in st[1:]]+[f"{P[s].loc[pdb,'inside']:.0f}"+("" if P[s].loc[pdb,'full'] else "$^{\\dagger}$") for s in st[1:]]
    L.append(" & ".join(r)+" \\\\")
open(S+"final/si_s8.tex","w").write("\n".join(L))
M=pd.read_csv(S+"final/astex_merged.csv"); R=M[M.sp_top1.notna()].copy()
NM={"cocryst":"Co-cryst.","fpocket":"fpocket","p2rank":"P2Rank","blind":"Blind"}
L=[]
for (ex,seed),g in R.groupby(["ex","seed"]):
    for s in st:
        h=g[g.strategy==s]; a=int(((h.nn_top1<2)&~(h.sp_top1<2)).sum()); b=int(((h.nn_best<2)&~(h.sp_best<2)).sum()); r_=int((~(h.nn_top1<2)&(h.sp_top1<2)).sum())+int((~(h.nn_best<2)&(h.sp_best<2)).sum())
        L.append(f"{ex} & {seed} & {NM[s]} & {len(h)} & {a} & {b} & {r_} \\\\")
open(S+"final/si_s9a.tex","w").write("\n".join(L))
D=R[((R.nn_top1<2)&~(R.sp_top1<2))|((R.nn_best<2)&~(R.sp_best<2))].sort_values(["pdb","ex","seed","strategy"])
L=[]
for r in D.itertuples():
    t="T" if (r.nn_top1<2 and not r.sp_top1<2) else ""; b="B" if (r.nn_best<2 and not r.sp_best<2) else ""
    L.append(f"{r.pdb} & {r.ex} & {r.seed} & {NM[r.strategy]} & {r.sp_top1:.2f} & {r.nn_top1:.2f} & {r.sp_best:.2f} & {r.nn_best:.2f} & {t}{b} \\\\")
open(S+"final/si_s9b.tex","w").write("\n".join(L)); print(len(L),"discordant rows; summary rows",36*4//4*4)
