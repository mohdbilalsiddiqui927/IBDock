import numpy as np, pandas as pd, pickle, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
for f in fm.findSystemFonts():
    if "texgyrepagella" in f.lower(): fm.fontManager.addfont(f)
plt.rcParams.update({"font.family":["TeX Gyre Pagella","DejaVu Serif"],"font.size":8,"axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":.9,"xtick.major.width":.9,"ytick.major.width":.9,"pdf.fonttype":42,"axes.labelsize":8.5,"legend.fontsize":7.5})
S="/tmp/claude-0/-home-claude/fa06f76e-b715-54d8-a3bb-7849315bf476/scratchpad/final/"
G=pd.read_csv(S+"astex_geometry.csv"); G=G[(G.ex==16)&(G.seed==42)]
M=pd.read_csv(S+"astex_merged.csv"); M=M[M.ex==16].merge(G[["strategy","pdb","vol","full"]],on=["strategy","pdb"])
M["y"]=(M.sp_top1<2)
ST=["cocryst","fpocket","p2rank","blind"]; LAB=["Co-\ncrystal","fpocket","P2Rank","Blind"]
COL=dict(zip(ST,["#2b5aa0","#18512f","#550022","#554400"])); INK="#222222"
fig,axs=plt.subplots(1,3,figsize=(7.4,3.1),gridspec_kw=dict(width_ratios=[1,1.15,1.05]))
rng=np.random.default_rng(1)
ax=axs[0]
for i,s in enumerate(ST):
    g=G[G.strategy==s]; v=np.log10(g.vol)
    ax.scatter(i+rng.uniform(-.28,.28,len(v)),v,s=7,color=COL[s],alpha=.55,lw=0,zorder=2)
    ax.hlines(np.median(v),i-.35,i+.35,color="black",lw=1.4,zorder=4)
ax.set_xticks(range(4)); ax.set_xticklabels(LAB,rotation=0,fontsize=7.2); ax.set_ylabel("Box volume (log$_{10}$ Å$^3$)",fontweight="bold")
ax.set_title("Box volume",fontsize=8.5,pad=6); ax.text(-.2,1.04,"a",transform=ax.transAxes,fontsize=12,fontweight="bold")
ax=axs[1]; x=0
for i,s in enumerate(ST):
    m=M[M.strategy==s]; 
    for flag,hatch,off in [(True,None,-.2),(False,"////",.2)]:
        g=m[m.full==flag]; 
        if len(g)==0: continue
        n=g.pdb.nunique(); v=100*g.y.mean()
        ax.bar(i+off,v,.36,color=COL[s],edgecolor=INK,lw=.7,hatch=hatch,alpha=1 if flag else .35,zorder=2)
        ax.text(i+off,v+2,f"n={n}" if flag else f"n={n}\n0%",ha="center",va="bottom",fontsize=6.3)
ax.set_xticks(range(4)); ax.set_xticklabels(LAB,fontsize=7.2); ax.set_ylim(0,80); ax.set_ylabel("Top-1 success (%, seed mean)",fontweight="bold")
ax.set_title("Success by box enclosure",fontsize=8.5,pad=6); ax.text(-.2,1.04,"b",transform=ax.transAxes,fontsize=12,fontweight="bold")
from matplotlib.patches import Patch
ax.legend([Patch(fc="#888",ec=INK),Patch(fc="#ddd",ec=INK,hatch="////")],["Ligand fully inside box","Ligand not fully inside"],loc="upper right",frameon=False,fontsize=6.6)
ax=axs[2]; res=pickle.load(open(S+"gee_res.pkl","rb"))
names=["fpocket","p2rank","blind"]; y=0; ticks=[];lab=[]
for s in names:
    for key,mk,fill in [("unadj","o",COL[s]),("adj_vol","s","white")]:
        o,lo,hi,_=res[key][s]
        ax.errorbar(o,y,xerr=[[o-lo],[hi-o]],fmt=mk,color=COL[s],mfc=fill,mec=COL[s],ms=4.5,capsize=2,lw=1,zorder=3); y-=1
    y-=.6
ax.axvline(1,color="black",lw=.8,ls="--",zorder=1); ax.set_xscale("log"); ax.set_xlim(.12,12)
ax.set_yticks([-.5,-3.1,-5.7]); ax.set_yticklabels(["fpocket","P2Rank","Blind"],fontsize=7.2); ax.set_xlabel("Odds ratio vs co-crystal (log scale)",fontweight="bold")
ax.set_title("Effect of volume adjustment",fontsize=8.5,pad=6); ax.text(-.2,1.04,"c",transform=ax.transAxes,fontsize=12,fontweight="bold")
from matplotlib.lines import Line2D
ax.legend([Line2D([],[],marker="o",ls="",mfc="#555",mec="#555",ms=4.5),Line2D([],[],marker="s",ls="",mfc="white",mec="#555",ms=4.5)],["Unadjusted","Adjusted for volume"],loc="center right",bbox_to_anchor=(1.0,.62),frameon=False,fontsize=6.6)
fig.tight_layout(w_pad=1.6); fig.savefig(S+"figs/Figure_box_enclosure.pdf",bbox_inches="tight"); fig.savefig(S+"figs/Figure_box_enclosure.png",dpi=300,bbox_inches="tight")
