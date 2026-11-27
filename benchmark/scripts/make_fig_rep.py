import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D
for f in fm.findSystemFonts():
    if "texgyrepagella" in f.lower(): fm.fontManager.addfont(f)
plt.rcParams.update({"font.family":["TeX Gyre Pagella","DejaVu Serif"],"font.size":8,"axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":.9,"pdf.fonttype":42,"axes.labelsize":8.5,"legend.fontsize":7.5})
S="/tmp/claude-0/-home-claude/fa06f76e-b715-54d8-a3bb-7849315bf476/scratchpad/final/"
ST=["cocryst","fpocket","p2rank","blind"]; LAB=["Co-crystal","fpocket","P2Rank","Blind"]; COL=dict(zip(ST,["#2b5aa0","#18512f","#550022","#554400"])); INK="#222222"
top1={"ABL1":[0.78,13.15,12.40,0.80],"KIT":[0.78,12.75,12.82,11.83]}; best={"ABL1":[0.78,11.78,11.66,0.80],"KIT":[0.78,6.75,10.28,11.72]}
score={"ABL1":[-12.31,-10.58,-10.65,-12.32],"KIT":[-12.61,-11.63,-11.43,-9.74],"SRC":[-9.53,-10.30,-10.04,-9.83],"EGFR":[-9.62,-9.60,-9.58,-9.61]}
fig,axs=plt.subplots(1,2,figsize=(7.2,3.0),gridspec_kw=dict(width_ratios=[1,1.15]))
ax=axs[0]
for gi,k in enumerate(["ABL1","KIT"]):
    for si,s in enumerate(ST):
        x=gi*5+si; ax.bar(x,top1[k][si],.85,color=COL[s],edgecolor=INK,lw=.7,zorder=2)
        ax.scatter(x,best[k][si],marker="D",s=16,facecolor="white",edgecolor="black",lw=.8,zorder=5)
    ax.text(gi*5+1.5,-.03,k+(" (original target)" if k=="ABL1" else " (repurposed target)"),ha="center",va="top",fontsize=7.4,fontweight="bold",transform=ax.get_xaxis_transform())
ax.set_yscale("log"); ax.set_ylim(.5,30); ax.axhline(2,color="black",lw=.9,ls="--",zorder=1); ax.text(9.9,2.15,"2 Å",ha="right",va="bottom",fontsize=7)
ax.set_xticks([]); ax.set_xlim(-.8,9.8); ax.set_ylabel("RMSD to crystal imatinib (Å, log)",fontweight="bold"); ax.set_title("Pose recovery",fontsize=8.5,pad=6); ax.text(-.2,1.04,"a",transform=ax.transAxes,fontsize=12,fontweight="bold")
ax=axs[1]
for gi,k in enumerate(["ABL1","KIT","SRC","EGFR"]):
    for si,s in enumerate(ST):
        ax.bar(gi*5+si,-score[k][si],.85,bottom=0,color=COL[s],edgecolor=INK,lw=.7,zorder=2)
    ax.text(gi*5+1.5,-.01,k,ha="center",va="top",fontsize=7.8,fontweight="bold",transform=ax.get_xaxis_transform())
    ax.text(gi*5+1.5,-.11,"IC$_{50}$ "+{"ABL1":"38 nM","KIT":"100 nM","SRC":">100 µM","EGFR":">100 µM"}[k],ha="center",va="top",fontsize=6.5,transform=ax.get_xaxis_transform())
ax.set_ylim(6,13.5); ax.set_xticks([]); ax.set_xlim(-.8,18.8); ax.set_ylabel("Top-1 Vina score (−kcal mol$^{-1}$)",fontweight="bold"); ax.set_title("Docking score (not comparable across receptors)",fontsize=8.5,pad=6); ax.text(-.15,1.04,"b",transform=ax.transAxes,fontsize=12,fontweight="bold")
h=[plt.Rectangle((0,0),1,1,fc=COL[s],ec=INK,lw=.7) for s in ST]+[Line2D([],[],marker="D",ls="",mfc="white",mec="black",ms=4)]
fig.legend(h,LAB+["Best saved pose"],loc="lower center",ncol=5,frameon=False,bbox_to_anchor=(.5,-.04))
fig.tight_layout(w_pad=1.8,rect=(0,.07,1,1)); fig.savefig(S+"figs/Figure_repurposing.pdf",bbox_inches="tight"); fig.savefig(S+"figs/Figure_repurposing.png",dpi=300,bbox_inches="tight")
