import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D
for f in fm.findSystemFonts():
    if "texgyrepagella" in f.lower(): fm.fontManager.addfont(f)
plt.rcParams.update({"font.family":["TeX Gyre Pagella","DejaVu Serif"],"font.size":8,"axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":.9,"xtick.major.width":.9,"ytick.major.width":.9,"pdf.fonttype":42,"axes.labelsize":8.5,"legend.fontsize":7.5})
S="/tmp/claude-0/-home-claude/fa06f76e-b715-54d8-a3bb-7849315bf476/scratchpad/final/"
T=pd.read_csv(S+"cross_main.csv"); Cd=pd.read_csv(S+"cross_cond.csv")
ST=["cocryst","fpocket","p2rank","blind"]; LAB=["Co-crystal","fpocket","P2Rank","Blind"]; TG=["ESR1","MAPK1","PPARG"]
COL=dict(zip(ST,["#2b5aa0","#18512f","#550022","#554400"])); INK="#222222"
fig,axs=plt.subplots(1,3,figsize=(7.4,3.1),gridspec_kw=dict(width_ratios=[1,1,1]))
def grp(ax,fn,ylab,ylim,lab,ci=True,cog=False):
    for gi,t in enumerate(TG):
        for si,s in enumerate(ST):
            x=gi*5+si; r=T[(T.target==t)&(T.strategy==s)].iloc[0]; v,lo,hi=fn(r,t,s)
            ax.bar(x,v,.85,color=COL[s],edgecolor=INK,lw=.7,zorder=2)
            if ci: ax.errorbar(x,v,yerr=[[v-lo],[hi-v]],color="black",capsize=2,lw=.9,zorder=4)
            if cog: ax.scatter(x,100*r.cog_t1/r.n_cog,marker="D",s=14,facecolor="white",edgecolor="black",lw=.8,zorder=6)
        ax.text(gi*5+1.5,-ylim*.07,t,ha="center",va="top",fontsize=8.5,fontweight="bold")
    ax.set_ylim(0,ylim); ax.set_xlim(-.8,13.8); ax.set_xticks([]); ax.set_ylabel(ylab,fontweight="bold"); ax.text(-.2,1.04,lab,transform=ax.transAxes,fontsize=12,fontweight="bold")
grp(axs[0],lambda r,t,s:(r.t1,r.t1_lo,r.t1_hi),"Cross-docking Top-1 success (%)",100,"a",cog=True)
grp(axs[1],lambda r,t,s:(r.site,r.site_lo,r.site_hi),"Site capture (%)",100,"b")
def cond(r,t,s): c=Cd[(Cd.target==t)&(Cd.strategy==s)].iloc[0]; return c.t1_given_site,c.t1_given_site,c.t1_given_site
grp(axs[2],cond,"Top-1 given site (%)",100,"c",ci=False)
for ax,tt in zip(axs,["Pose recovery","Site capture","Pose given site"]): ax.set_title(tt,fontsize=8.5,pad=6)
h=[plt.Rectangle((0,0),1,1,fc=COL[s],ec=INK,lw=.7) for s in ST]+[Line2D([],[],marker="D",ls="",mfc="white",mec="black",ms=4)]
fig.legend(h,LAB+["Cognate (redocking) Top-1"],loc="lower center",ncol=5,frameon=False,bbox_to_anchor=(.5,-.04))
fig.tight_layout(w_pad=1.6,rect=(0,.05,1,1)); fig.savefig(S+"figs/Figure6_crossdocking.pdf",bbox_inches="tight"); fig.savefig(S+"figs/Figure6_crossdocking.png",dpi=600,bbox_inches="tight")
