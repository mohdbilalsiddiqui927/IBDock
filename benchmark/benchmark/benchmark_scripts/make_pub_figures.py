"""Publication figures for IBDock Provenire (main Figures 4-9, SI Figures S1-S4).
Every panel is written as its own file (Figure_4a.svg/.pdf/.png ...) and as a composite (Figure_4.*) built from the same panel code.
No grid lines; Okabe-Ito palette; text kept as text in SVG; file names contain no spaces."""
import numpy as np, pandas as pd, pickle, re, os, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy import stats
import sys, glob
# usage: python make_pub_figures.py [DATA_ROOT] [OUT_DIR]   (defaults: ../Data  and  ../Figures)
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=sys.argv[1] if len(sys.argv)>1 else os.path.join(HERE,"..","Data"); OUT=(sys.argv[2] if len(sys.argv)>2 else os.path.join(HERE,"..","Figures"))+os.sep
def _find(n):
    h=glob.glob(os.path.join(ROOT,"**",n),recursive=True)
    if not h: raise FileNotFoundError(n+" under "+ROOT)
    return h[0]
class _D(str):
    def __add__(self,n): return _find(n)
D=_D("")
os.makedirs(OUT+"Main",exist_ok=True); os.makedirs(OUT+"SI",exist_ok=True)
plt.rcParams.update({"font.family":["TeX Gyre Pagella","Palatino Linotype","DejaVu Serif"],"mathtext.default":"regular","font.size":7.5,"axes.labelsize":7.8,"axes.titlesize":8,"xtick.labelsize":7.2,"ytick.labelsize":7.2,
 "legend.fontsize":6.4,"axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":.8,"xtick.major.width":.8,"ytick.major.width":.8,"xtick.major.size":3,"ytick.major.size":3,
 "axes.grid":False,"svg.fonttype":"none","pdf.fonttype":42,"ps.fonttype":42,"savefig.transparent":False,"figure.facecolor":"white"})
ST=["cocryst","fpocket","p2rank","blind"]; LAB=dict(cocryst="Co-crystal",fpocket="fpocket",p2rank="P2Rank",blind="Blind"); LAB2=["Co-crystal","fpocket","P2Rank","Blind"]
COL=dict(cocryst="#2b5aa0",fpocket="#18512f",p2rank="#550022",blind="#554400"); MK=dict(cocryst="o",fpocket="s",p2rank="^",blind="D"); INK="#1a1a1a"; GREY="#BBBBBB"; MID="#8C8C8C"
M=pd.read_csv(D+"astex_merged.csv"); G=pd.read_csv(D+"astex_geometry.csv"); G16=G[(G.ex==16)&(G.seed==42)]
M16=M[M.ex==16].copy(); M16["y1"]=(M16.sp_top1<2).astype(int); M16["yb"]=(M16.sp_best<2).astype(int); PDB=sorted(M16.pdb.unique())
T2=dict(top1=dict(cocryst=(59.6,49.4,69.4),fpocket=(51.4,40.8,61.6),p2rank=(36.9,27.1,47.5),blind=(29.0,21.6,36.9)),best=dict(cocryst=(72.5,63.1,81.6),fpocket=(63.9,53.7,73.7),p2rank=(44.7,34.5,55.3),blind=(37.3,29.0,45.9)))  # Table 2
def piv(df,col,s): return df[df.strategy==s].pivot(index="pdb",columns="seed",values=col).reindex(PDB)
def letter(ax,l,x=None,y=None):
    if not l: return
    fw,fh=ax.figure.get_size_inches(); p=ax.get_position()
    ax.text(((getattr(ax,'_cell',0)+.06)/fw-p.x0)/p.width,((fh-.24)/fh-p.y0)/p.height,l,transform=ax.transAxes,fontsize=11,fontweight="bold",va="bottom",ha="left")
def xlab_strats(ax,pos=None,rot=0):
    ax.set_xticks(range(4) if pos is None else pos); ax.set_xticklabels(LAB2,rotation=rot)
def strat_legend(ax,loc="upper right",**k): ax.legend([Patch(fc=COL[s],ec=INK,lw=.5) for s in ST],LAB2,frameon=False,loc=loc,**k)
from matplotlib.axes import Axes
_orig_legend=Axes.legend
def _leg(self,handles=None,labels=None,**k):
    wide=getattr(self,"_wide",False); nc=3 if wide else (1 if getattr(self,"_onecol",False) else 2)
    return _orig_legend(self,handles,labels,loc="upper center",bbox_to_anchor=(.5,getattr(self,"_legy",-.2)),ncol=nc,frameon=False,fontsize=6.3,handlelength=1.2,columnspacing=.9,labelspacing=.3,handletextpad=.4)
Axes.legend=_leg
# ---------------- Figure 4 ----------------
def f4a(ax,l=None,legend=True):
    seed_rates={}
    for ep,y in (("top1","y1"),("best","yb")):
        for s in ST: seed_rates[(ep,s)]=[100*piv(M16,y,s)[sd].sum()/85 for sd in (1,42,123)]
    unres={(ep,s):100*piv(M16,col,s).isna().sum().mean()/85 for ep,col in (("top1","sp_top1"),("best","sp_best")) for s in ST}
    w=.36
    for i,s in enumerate(ST):
        for j,ep in enumerate(("top1","best")):
            x=i+(j-.5)*(w+.04); m,lo,hi=T2[ep][s]; un=unres[(ep,s)]
            ax.bar(x,m,w,color=COL[s],edgecolor=INK,lw=.6,alpha=1 if j==1 else .55,zorder=2)
            ax.bar(x,100-un-m,w,bottom=m,color="#E4E4E4",edgecolor=INK,lw=.4,zorder=1)
            ax.bar(x,un,w,bottom=100-un,color="white",edgecolor=MID,lw=.4,hatch="////",zorder=1)
            ax.errorbar(x,m,yerr=[[m-lo],[hi-m]],color=INK,capsize=2,lw=.8,zorder=4)
            ax.scatter([x]*3,seed_rates[(ep,s)],s=7,facecolor="white",edgecolor=INK,lw=.6,zorder=5)
            ax.text(x,min(hi+1.2,99),f"{m:.1f}",ha="center",va="bottom",fontsize=6.3,color=INK,bbox=dict(fc="white",ec="none",pad=.6,alpha=.85),zorder=6)
    ax.set_xlim(-.6,3.6); ax.set_ylim(0,100); ax.set_xticks(range(4)); ax.set_xticklabels(LAB2,fontweight="bold"); ax.tick_params(axis="x",pad=12,length=0)
    for i in range(4):
        for j,t in enumerate(("Top-1","Best")): ax.text(i+(j-.5)*(w+.04),-2,t,ha="center",va="top",fontsize=6.2)
    ax.set_ylabel("Complexes (% of 85)"); letter(ax,l)
    if legend: ax.legend([Patch(fc=MID,ec=INK,lw=.5),Patch(fc="#E4E4E4",ec=INK,lw=.5),Patch(fc="white",ec=MID,hatch="////",lw=.5),Line2D([],[],marker="o",ls="",mfc="white",mec=INK,ms=3.5),Line2D([],[],color=INK,lw=.8)],
        ["Success (RMSD < 2.0 Å)","Resolved, not success","sPyRMSD unresolved","Individual seed","95% bootstrap CI"],loc="upper center",bbox_to_anchor=(.5,-.17),ncol=3,frameon=False,fontsize=6.5,columnspacing=1,handlelength=1.4)
def f4b(ax,l=None,legend=True):
    rng=np.random.default_rng(3); ok=np.all([piv(M16,"sp_top1",s).notna().all(1).values for s in ST],axis=0)
    for i,s in enumerate(ST):
        for j,(col,ep) in enumerate((("sp_top1","Top-1"),("sp_best","Best"))):
            v=piv(M16,col,s).median(1).values[ok]; x=i+(j-.5)*.42
            ax.scatter(x+rng.uniform(-.1,.1,len(v)),v,s=4,color=MID,alpha=.5,lw=0,zorder=1)
            q=np.percentile(v,[10,25,50,75,90]); ax.add_patch(plt.Rectangle((x-.14,q[1]),.28,q[3]-q[1],fc=COL[s],alpha=1 if j==1 else .55,ec=INK,lw=.7,zorder=3))
            ax.plot([x-.14,x+.14],[q[2]]*2,color=INK,lw=1.4,zorder=4); ax.plot([x,x],[q[0],q[1]],color=INK,lw=.8,zorder=3); ax.plot([x,x],[q[3],q[4]],color=INK,lw=.8,zorder=3)
            ax.text(x,.135,ep,ha="center",va="top",fontsize=6.3,transform=ax.get_xaxis_transform()) if False else None
    ax.set_yscale("log"); ax.set_ylim(.12,60); ax.set_yticks([.2,.5,1,2,5,10,20,50]); ax.set_yticklabels(["0.2","0.5","1","2","5","10","20","50"])
    ax.axhline(2,color=INK,lw=.8,ls="--",zorder=2); pass
    ax.set_xlim(-.6,3.6); ax.set_xticks(range(4)); ax.set_xticklabels(LAB2,fontweight="bold"); ax.set_ylabel("sPyRMSD (Å, log)"); letter(ax,l)
    if legend: ax.legend([Patch(fc=MID,ec=INK,alpha=.55,lw=.5),Patch(fc=MID,ec=INK,lw=.5),Line2D([],[],color=INK,lw=.8,ls="--")],["Top-1 pose (left)","Best pose (right)","2.0 Å threshold"])
# ---------------- Figure 5 ----------------
def f5a(ax,l=None,legend=True):
    rng=np.random.default_rng(1)
    for i,s in enumerate(ST):
        v=G16[G16.strategy==s].vol.values; ax.scatter(i+rng.uniform(-.28,.28,len(v)),v,s=7,color=COL[s],alpha=.6,lw=0,zorder=2); ax.hlines(np.median(v),i-.38,i+.38,color=INK,lw=1.6,zorder=4)
        ax.text(i+.42,np.median(v),f"{np.median(v):,.0f}",va="center",ha="left",fontsize=6.3)
    ax.set_yscale("log"); ax.set_ylim(1e3,2e6); ax.set_xlim(-.6,3.95); ax.set_xticks(range(4)); ax.set_xticklabels(LAB2,fontsize=6.8,rotation=30,ha="right",rotation_mode="anchor"); ax.set_ylabel("Box volume (Å³, log)"); letter(ax,l)
    if legend: ax.legend([Line2D([],[],color=INK,lw=1.6)],["Median"],frameon=False,loc="lower right",fontsize=6.5,bbox_to_anchor=(1,.02))
def f5b(ax,l=None,legend=True):
    Z=M16.merge(G16[["strategy","pdb","vol","full"]],on=["strategy","pdb"])
    for i,s in enumerate(ST):
        z=Z[Z.strategy==s]
        for flag,off,al,h in ((True,-.2,1,None),(False,.2,.4,"////")):
            g=z[z.full==flag]; n=g.pdb.nunique()
            if n==0: continue
            v=100*g.y1.mean(); ax.bar(i+off,v,.36,color=COL[s],edgecolor=INK,lw=.6,alpha=al,hatch=h,zorder=2); ax.text(i+off,v+1.5,f"{v:.0f}%\nn={n}",ha="center",va="bottom",fontsize=6.0,linespacing=1.1)
    ax.set_ylim(0,92); ax.set_xlim(-.6,3.6); ax.set_xticks(range(4)); ax.set_xticklabels(LAB2,fontsize=6.8); ax.set_ylabel("Top-1 success (%)"); letter(ax,l)
    if legend: ax.legend([Patch(fc=MID,ec=INK,lw=.5),Patch(fc="#DDDDDD",ec=INK,hatch="////",lw=.5)],["Ligand fully inside box","Ligand not fully inside box"],loc="upper right",bbox_to_anchor=(1,1.02),frameon=False,fontsize=6.5)
def f5c(ax,l=None,legend=True):
    res=pickle.load(open(D+"gee_res.pkl","rb")); y=0; ticks=[]
    for s in ("fpocket","p2rank","blind"):
        y0=y
        for key,mk,fill in (("unadj","o",COL[s]),("adj_vol","s","white")):
            o,lo,hi,_=res[key][s]; ax.errorbar(o,y,xerr=[[o-lo],[hi-o]],fmt=mk,color=COL[s],mfc=fill,mec=COL[s],ms=4.5,capsize=2,lw=1,zorder=3); y-=1
        ticks.append((y0+y+1)/2); y-=.7
    ax.axvline(1,color=INK,lw=.8,ls="--",zorder=1); ax.set_xscale("log"); ax.set_xlim(.12,12); ax.set_xticks([.2,.5,1,2,5,10]); ax.set_xticklabels(["0.2","0.5","1","2","5","10"])
    ax.set_yticks(ticks); ax.set_yticklabels(["fpocket","P2Rank","Blind"]); ax.set_xlabel("Odds ratio vs co-crystal (95% CI)"); letter(ax,l,x=-.28)
    if legend: ax.legend([Line2D([],[],marker="o",ls="",mfc=MID,mec=MID,ms=4.5),Line2D([],[],marker="s",ls="",mfc="white",mec=MID,ms=4.5)],["Unadjusted","Adjusted for log box volume"],loc="upper center",bbox_to_anchor=(.35,-.2),ncol=2,columnspacing=1,frameon=False,fontsize=6.5)
# ---------------- Figure 6 ----------------
def _spread(v,gap):
    o=np.argsort(v); out=np.array(v,float)
    for _ in range(50):
        for a,b in zip(o[:-1],o[1:]):
            if out[b]-out[a]<gap: m=(out[a]+out[b])/2; out[a]=m-gap/2; out[b]=m+gap/2
    return out
def _grid_line(ax,y,l,ylab,legend):
    ends={}
    for s in ST:
        r=[]; per=[]
        for ex in (8,16,32):
            d=M[(M.ex==ex)&(M.strategy==s)]; ps=[100*(d[d.seed==sd][y]<2).sum()/85 for sd in (1,42,123)]; per.append(ps); r.append(np.mean(ps))
        ax.plot([8,16,32],r,color=COL[s],marker=MK[s],ms=4.5,lw=1.3,zorder=3)
        for e,ps in zip((8,16,32),per): ax.scatter([e]*3,ps,s=9,facecolor="white",edgecolor=COL[s],lw=.7,zorder=4)
        ends[s]=r[-1]
    ys=_spread([ends[s] for s in ST],4.2)
    for s,yy in zip(ST,ys): ax.text(33.6,yy,LAB[s],color=INK,va="center",fontsize=6.8)
    ax.set_xlim(6,40); ax.set_xticks([8,16,32]); ax.set_ylim(0,80); ax.set_xlabel("Vina exhaustiveness"); ax.set_ylabel(ylab); letter(ax,l)
    if legend: ax.legend([Line2D([],[],marker="o",ls="",mfc="white",mec=MID,ms=3.5)],["Individual seed"],loc="lower right",frameon=False,fontsize=6.5)
def f6a(ax,l=None,legend=True): _grid_line(ax,"sp_top1",l,"Top-1 success (%)",legend)
def f6b(ax,l=None,legend=True): _grid_line(ax,"sp_best",l,"Best-pose success (%)",legend)
def f6c(ax,l=None,legend=True):
    C=dict(all="#2b5aa0",flip="#a9bbdb",none="#CFCFCF"); y=0; yt=[]; yl=[]
    for si,s in enumerate(ST):
        for ex in (8,16,32):
            d=M[(M.ex==ex)&(M.strategy==s)].pivot(index="pdb",columns="seed",values="sp_top1").reindex(PDB); k=(d<2).sum(1)
            a=int((k==3).sum()); f=int(((k>0)&(k<3)).sum()); n=int((k==0).sum())
            ax.barh(y,a,.8,color=C["all"],edgecolor=INK,lw=.4); ax.barh(y,f,.8,left=a,color=C["flip"],edgecolor=INK,lw=.4); ax.barh(y,n,.8,left=a+f,color=C["none"],edgecolor=INK,lw=.4)
            ax.text(a/2,y,str(a),ha="center",va="center",fontsize=6.2,color="white") if a>=6 else None
            ax.text(a+f+n/2,y,str(n),ha="center",va="center",fontsize=6.2,color=INK) if n>=6 else None
            ax.text(a+f/2,y,str(f),ha="center",va="center",fontsize=6.0,color=INK) if f>=9 else None
            yt.append(y); yl.append(f"ex {ex}"); y+=1
        fw_=ax.figure.get_size_inches()[0]; aw_=ax.get_position().width*fw_; lx=-(ax.get_position().x0*fw_-ax._cell-.14)/aw_
        ax.text(lx,y-2,LAB[s],transform=ax.get_yaxis_transform(),rotation=90,ha="left",va="center",fontweight="bold",fontsize=7.5) ; y+=.7
    ax.set_yticks(yt); ax.set_yticklabels(yl,fontsize=6.8); ax.invert_yaxis(); ax.set_xlim(0,85); ax.set_xlabel("Complexes (n = 85)"); letter(ax,l,x=-.30)
    if legend: ax.legend([Patch(fc=C["all"],ec=INK,lw=.4),Patch(fc=C["flip"],ec=INK,lw=.4),Patch(fc=C["none"],ec=INK,lw=.4)],["Success in all 3 seeds","Outcome flips between seeds","Failure in all 3 seeds"],loc="upper center",bbox_to_anchor=(.45,-.14),ncol=3,frameon=False,fontsize=6.5,columnspacing=1,handlelength=1.2)
# ---------------- Figure 7 ----------------
def _bridge(ax,xc,yc,l,title):
    X=M16.dropna(subset=[xc,yc]); n=len(X)
    for s in ST[::-1]:
        z=X[X.strategy==s]; ax.scatter(z[xc],z[yc],s=6,color=COL[s],alpha=.55,lw=0,zorder=2)
    ax.plot([.1,60],[.1,60],color=INK,lw=.8,zorder=1); ax.axhline(2,color=MID,lw=.7,ls="--"); ax.axvline(2,color=MID,lw=.7,ls="--")
    rho=stats.spearmanr(X[xc],X[yc])[0]; mae=(X[xc]-X[yc]).abs().mean()
    ax.text(.04,.96,f"n = {n}\nSpearman ρ = {rho:.3f}\nMAE = {mae:.3f} Å",transform=ax.transAxes,va="top",ha="left",fontsize=6.8)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(.12,60); ax.set_ylim(.12,60); tk=[.2,.5,1,2,5,10,20,50]
    ax.set_xticks(tk); ax.set_yticks(tk); ax.set_xticklabels([str(t) for t in tk]); ax.set_yticklabels([str(t) for t in tk]); ax.set_xlabel("sPyRMSD (Å, log)"); ax.set_ylabel("Bridge RMSD (Å, log)"); ax.set_aspect("equal"); letter(ax,l,x=-.2); ax.set_title(title,fontsize=7.5,pad=4); print("bridge",xc,n,round(rho,3),round(mae,3))
def f7a(ax,l=None,legend=True):
    _bridge(ax,"sp_top1","br_top1",l,"Top-1 pose")
    if legend: ax.legend([Line2D([],[],marker="o",ls="",color=COL[s],ms=3.5) for s in ST]+[Line2D([],[],color=INK,lw=.8)],LAB2+["Identity"],loc="lower right",frameon=False,fontsize=6.3)
def f7b(ax,l=None,legend=True):
    _bridge(ax,"sp_best","br_best",l,"Best pose")
    if legend: ax.legend([Line2D([],[],marker="o",ls="",color=COL[s],ms=3.5) for s in ST]+[Line2D([],[],color=INK,lw=.8)],LAB2+["Identity"],loc="lower right",frameon=False,fontsize=6.3)
# ---------------- Figure 8 ----------------
T=pd.read_csv(D+"cross_main.csv"); Cd=pd.read_csv(D+"cross_cond.csv"); TG=["ESR1","MAPK1","PPARG"]
def _cross(ax,fn,ylab,l,ci=True,cog=False,legend=True):
    for gi,t in enumerate(TG):
        for si,s in enumerate(ST):
            x=gi*5+si; r=T[(T.target==t)&(T.strategy==s)].iloc[0]; v,lo,hi=fn(r,t,s)
            ax.bar(x,v,.85,color=COL[s],edgecolor=INK,lw=.6,zorder=2)
            if ci: ax.errorbar(x,v,yerr=[[v-lo],[hi-v]],color=INK,capsize=1.8,lw=.8,zorder=4)
            if cog: ax.scatter(x,100*r.cog_t1/r.n_cog,marker="D",s=13,facecolor="white",edgecolor=INK,lw=.8,zorder=6)
    ax.set_xticks([1.5,6.5,11.5]); ax.set_xticklabels(TG,fontweight="bold"); ax.tick_params(axis="x",length=0); ax.set_ylim(0,100); ax.set_xlim(-.8,13.8); ax.set_ylabel(ylab); letter(ax,l)
    if legend:
        h=[Patch(fc=COL[s],ec=INK,lw=.5) for s in ST]+([Line2D([],[],marker="D",ls="",mfc="white",mec=INK,ms=3.8)] if cog else [])
        ax.legend(h,LAB2+(["Cognate (redocking)"] if cog else []),loc="upper right",bbox_to_anchor=(1,.86) if cog else (1,1),frameon=False,fontsize=6.3,handlelength=1.2)
def f8a(ax,l=None,legend=True): _cross(ax,lambda r,t,s:(r.t1,r.t1_lo,r.t1_hi),"Cross-docking Top-1 (%)",l,cog=True,legend=legend)
def f8b(ax,l=None,legend=True): _cross(ax,lambda r,t,s:(r.site,r.site_lo,r.site_hi),"Site capture (%)",l,legend=legend)
def f8c(ax,l=None,legend=True):
    def c(r,t,s): v=Cd[(Cd.target==t)&(Cd.strategy==s)].iloc[0].t1_given_site; return v,v,v
    _cross(ax,c,"Top-1 given site capture (%)",l,ci=False,legend=legend)
# ---------------- Figure 9 ----------------
top1={"ABL1":[0.78,13.15,12.40,0.80],"KIT":[0.78,12.75,12.82,11.83]}; best={"ABL1":[0.78,11.78,11.66,0.80],"KIT":[0.78,6.75,10.28,11.72]}
score={"ABL1":[-12.31,-10.58,-10.65,-12.32],"KIT":[-12.61,-11.63,-11.43,-9.74],"SRC":[-9.53,-10.30,-10.04,-9.83],"EGFR":[-9.62,-9.60,-9.58,-9.61]}
IC={"ABL1":"38 nM","KIT":"100 nM","SRC":">100 µM","EGFR":">100 µM"}
def f9a(ax,l=None,legend=True):
    for gi,k in enumerate(["ABL1","KIT"]):
        for si,s in enumerate(ST):
            x=gi*5+si; ax.bar(x,top1[k][si],.85,color=COL[s],edgecolor=INK,lw=.6,zorder=2); ax.scatter(x,best[k][si],marker="D",s=15,facecolor="white",edgecolor=INK,lw=.8,zorder=5)
    ax.set_yscale("log"); ax.set_ylim(.5,30); ax.set_yticks([.5,1,2,5,10,20]); ax.set_yticklabels(["0.5","1","2","5","10","20"]); ax.axhline(2,color=INK,lw=.8,ls="--",zorder=1); ax.text(9.8,2.15,"2.0 Å",ha="right",va="bottom",fontsize=6.5)
    ax.set_xticks([1.5,6.5]); ax.set_xticklabels(["ABL1 (1IEP)","KIT (1T46)"],fontweight="bold"); ax.tick_params(axis="x",length=0); ax.set_xlim(-.8,9.8); ax.set_ylabel("RMSD to crystal (Å, log)"); letter(ax,l)
    if legend: ax.legend([Patch(fc=COL[s],ec=INK,lw=.5) for s in ST]+[Line2D([],[],marker="D",ls="",mfc="white",mec=INK,ms=3.8)],LAB2+["Best saved pose"],loc="upper left",bbox_to_anchor=(0,1.02),frameon=False,fontsize=6.3,ncol=2,columnspacing=1,handlelength=1.2)
def f9b(ax,l=None,legend=True):
    for gi,k in enumerate(["ABL1","KIT","SRC","EGFR"]):
        for si,s in enumerate(ST): ax.bar(gi*5+si,-score[k][si],.85,color=COL[s],edgecolor=INK,lw=.6,zorder=2)
    ax.set_ylim(0,14.5); ax.set_xlim(-.8,18.8); ax.set_xticks([1.5,6.5,11.5,16.5]); ax.set_xticklabels([f"{k}\nIC$_{{50}}$ {IC[k]}" for k in score],fontsize=7); ax.tick_params(axis="x",length=0)
    ax.set_ylabel("Top-1 Vina score (−kcal/mol)"); letter(ax,l)
    if legend: ax.legend([Patch(fc=COL[s],ec=INK,lw=.5) for s in ST],LAB2,loc="upper right",frameon=False,fontsize=6.3,ncol=4,columnspacing=.9,handlelength=1.1)
# ---------------- SI ----------------
def fS1(ax,l=None,col="sp_top1",title=""):
    for s in ST:
        v=np.sort(piv(M16,col,s).median(1).dropna().values); y=100*np.arange(1,len(v)+1)/85
        ax.step(np.r_[v,60],np.r_[y,y[-1]],where="post",color=COL[s],lw=1.5,zorder=3)
        k=100*(v<2).sum()/85; ax.scatter([2],[k],s=22,color=COL[s],edgecolor=INK,lw=.6,zorder=5,marker=MK[s])
    ax.axvline(2,color=INK,lw=.8,ls="--",zorder=1); ax.set_xscale("log"); ax.set_xlim(.12,60); ax.set_xticks([.2,.5,1,2,5,10,20,50]); ax.set_xticklabels(["0.2","0.5","1","2","5","10","20","50"])
    ax.set_ylim(0,100); ax.set_xlabel("RMSD threshold (Å, log)"); ax.set_ylabel("Complexes below (% of 85)"); letter(ax,l); ax.set_title(title,fontsize=7.5,pad=4)
    ax.legend([Line2D([],[],color=COL[s],lw=1.5,marker=MK[s],ms=4,mec=INK,mew=.5) for s in ST],LAB2,loc="lower right",frameon=False,fontsize=6.5)
def fS1a(ax,l=None,legend=True): fS1(ax,l,"sp_top1","Top-1 pose")
def fS1b(ax,l=None,legend=True): fS1(ax,l,"sp_best","Best pose")
def _consensus():
    out={}
    for s in ST:
        a=piv(M16,"sp_top1",s); b=piv(M16,"sp_best",s); c=piv(M16,"br_top1",s); cb=piv(M16,"br_best",s)
        succ=((a<2).sum(1)>=2); near=(~succ)&((b<2).sum(1)>=2); unr=(a.isna().sum(1)>=2)&~succ
        out[s]=dict(succ=int(succ.sum()),near=int((near&~unr).sum()),unr=int(unr.sum()))
        out[s]["none"]=85-out[s]["succ"]-out[s]["near"]-out[s]["unr"]
        sp=(a.notna().sum(1)>=2); br=(c.notna().sum(1)>=2); out[s].update(both=int((sp&br).sum()),sponly=int((sp&~br).sum()),bronly=int((~sp&br).sum()),neither=int((~sp&~br).sum()))
    return out
CON=_consensus(); print(CON)
def _stack(ax,rows,cats,leg,l,tint=True):
    for i,s in enumerate(ST):
        b=0
        for (k,fcn,h) in cats:
            v=CON[s][k]; fc=fcn(s); ax.bar(i,v,.62,bottom=b,color=fc,edgecolor=INK,lw=.5,hatch=h,zorder=2)
            if v>=3: ax.text(i,b+v/2,str(v),ha="center",va="center",fontsize=6.8,color="white" if (k in("succ","both") and s!="p2rank") else INK,bbox=dict(fc="white",ec="none",pad=1) if k in("unr","neither") else None)
            b+=v
    ax.set_ylim(0,85); ax.set_xlim(-.6,3.6); xlab_strats(ax); ax.set_ylabel("Complexes (n = 85)"); letter(ax,l); ax.legend(*leg,loc="upper center",bbox_to_anchor=(.5,-.12),frameon=False,fontsize=6.5,ncol=2,columnspacing=1,handlelength=1.3)
def tint(c,f=.45): 
    c=np.array(matplotlib.colors.to_rgb(c)); return tuple(c+(1-c)*(1-f))
def fS2(ax,l=None,legend=True):
    cats=[("succ",lambda s:COL[s],None),("near",lambda s:tint(COL[s]),None),("none",lambda s:"#D9D9D9",None),("unr",lambda s:"white","////")]
    h=[Patch(fc=MID,ec=INK,lw=.5),Patch(fc="#DDDDDD",ec=INK,lw=.5),Patch(fc="#D9D9D9",ec=INK,lw=.5),Patch(fc="white",ec=MID,hatch="////",lw=.5)]
    _stack(ax,None,cats,(h,["Top-1 pose < 2.0 Å (strategy colour)","Near-native pose saved, not ranked first (light tint)","No pose < 2.0 Å generated","sPyRMSD unresolved"]),l)
def fS3(ax,l=None,legend=True):
    cats=[("both",lambda s:COL[s],None),("sponly",lambda s:tint(COL[s]),None),("bronly",lambda s:"#BDBDBD",None),("neither",lambda s:"white","////")]
    h=[Patch(fc=MID,ec=INK,lw=.5),Patch(fc="#DDDDDD",ec=INK,lw=.5),Patch(fc="#BDBDBD",ec=INK,lw=.5),Patch(fc="white",ec=MID,hatch="////",lw=.5)]
    _stack(ax,None,cats,(h,["Both routes numerical (strategy colour)","sPyRMSD only (light tint)","Coordinate bridge only","Neither route"]),l)
def _s4data():
    ok=np.all([piv(M16,"sp_top1",s).notna().all(1).values for s in ST],axis=0); sc={s:piv(M16,"score",s).mean(1).values[ok] for s in ST}; return ok,sc
def fS4a(ax,l=None,legend=True):
    from statsmodels.stats.multitest import multipletests
    ok,sc=_s4data(); parts=ax.violinplot([sc[s] for s in ST],positions=range(4),widths=.8,showextrema=False)
    for b,s in zip(parts["bodies"],ST): b.set_facecolor(COL[s]); b.set_edgecolor(INK); b.set_linewidth(.6); b.set_alpha(.8)
    for i,s in enumerate(ST):
        q=np.percentile(sc[s],[25,50,75]); ax.plot([i,i],[q[0],q[2]],color=INK,lw=3,solid_capstyle="butt",zorder=4); ax.scatter(i,q[1],s=14,color="white",edgecolor=INK,zorder=5,lw=.6)
    import itertools
    pr=list(itertools.combinations(range(4),2)); p=[stats.wilcoxon(sc[ST[a]],sc[ST[b]]).pvalue for a,b in pr]; ph=multipletests(p,method="holm")[1]; y=1.0; print("S4a holm",[round(x,4) for x in ph])
    for (a,b),q in zip(pr,ph):
        if q<.05:
            ax.plot([a,a,b,b],[y,y+.35,y+.35,y],color=INK,lw=.7); ax.text((a+b)/2,y+.3,"***" if q<.001 else "**" if q<.01 else "*",ha="center",va="bottom",fontsize=7); y+=1.2
    ax.set_ylim(-14.5,y-.6); ax.set_xlim(-.6,3.6); xlab_strats(ax); ax.set_ylabel("Top-1 Vina score (kcal/mol)"); letter(ax,l)
def fS4b(ax,l=None,legend=True):
    rng=np.random.default_rng(5); ok,sc=_s4data(); rec_all={}
    for i,s in enumerate(ST):
        a=piv(M16,"sp_top1",s); score=piv(M16,"score",s).mean(1); rec=((a.fillna(99)<=2).sum(1)>=2); rec_all[s]=(score,rec)
        okS=pd.Series(ok,index=PDB); sc_=score[okS]; r=rec[okS]
        for j,(flag,off,al) in enumerate(((True,-.2,1),(False,.2,.35))):
            v=sc_[r==flag].values
            if len(v)<3: continue
            pt=ax.violinplot([v],positions=[i+off],widths=.36,showextrema=False)["bodies"][0]; pt.set_facecolor(COL[s] if flag else "#CFCFCF"); pt.set_edgecolor(INK); pt.set_linewidth(.5); pt.set_alpha(.85 if flag else .9)
            q=np.percentile(v,[25,50,75]); ax.plot([i+off]*2,[q[0],q[2]],color=INK,lw=2.4,solid_capstyle="butt",zorder=4); ax.scatter(i+off,q[1],s=10,color="white",edgecolor=INK,lw=.5,zorder=5)
        y=r.values.astype(int); x=-sc_.values
        if y.min()!=y.max():
            auc=stats.mannwhitneyu(x[y==1],x[y==0]).statistic/((y==1).sum()*(y==0).sum()); bs=[]
            for _ in range(5000):
                k=rng.integers(0,len(y),len(y)); yy=y[k]; xx=x[k]
                if yy.min()==yy.max(): continue
                bs.append(stats.mannwhitneyu(xx[yy==1],xx[yy==0]).statistic/((yy==1).sum()*(yy==0).sum()))
            lo,hi=np.percentile(bs,[2.5,97.5]); ax.text(i,-3.0,f"AUC {auc:.2f}\n({lo:.2f}–{hi:.2f})",ha="center",va="bottom",fontsize=6.3); print("AUC",s,round(auc,2),round(lo,2),round(hi,2))
    ax.set_ylim(-14.5,-1.5); ax.set_xlim(-.6,3.6); xlab_strats(ax); ax.set_ylabel("Top-1 Vina score (kcal/mol)"); letter(ax,l)
    ax.legend([Patch(fc=MID,ec=INK,lw=.5),Patch(fc="#CFCFCF",ec=INK,lw=.5)],["Recovered (strategy colour)","Not recovered"],loc="lower right",frameon=False,fontsize=6.5,bbox_to_anchor=(1,-.02))
# ---------------- writers ----------------
def _post(path):
    t=open(path).read()
    t=re.sub(r"font: ([^;\"]*?)'(TeX Gyre Pagella|Palatino Linotype|DejaVu Serif)'","font: \\1'Palatino Linotype', Palatino, 'TeX Gyre Pagella', 'Book Antiqua', serif",t)
    t=re.sub(r"font-family: [^;\"]*;","font-family: 'Palatino Linotype', Palatino, 'TeX Gyre Pagella', 'Book Antiqua', serif;",t); open(path,"w").write(t)
def _save(fig,name,sub):
    for ext in ("svg","pdf"): fig.savefig(f"{OUT}{sub}/{name}.{ext}")
    fig.savefig(f"{OUT}{sub}/{name}.png",dpi=600); _post(f"{OUT}{sub}/{name}.svg"); plt.close(fig)
CW,CH=3.5,3.55
def _axes(fig,x0,B,aw,h,L):
    W,H=fig.get_size_inches(); a=fig.add_axes([(x0+L)/W,B/H,aw/W,h/H]); a._cell=x0; return a
def panel(fn,name,sub,letter_=None,tweak=None):
    fig=plt.figure(figsize=(CW,CH)); ax=_axes(fig,0,1.12,CW-.66-.12,2.11,.66)
    if tweak: tweak(ax)
    fn(ax,letter_,True); _save(fig,name,sub)
def comp(fns,name,sub,letters,tweaks):
    n=len(fns); W=7.2; cw=W/n; fig=plt.figure(figsize=(W,CH)); L=.66 if n==2 else .56
    for i,(f,l) in enumerate(zip(fns,letters)):
        ax=_axes(fig,i*cw,1.12,cw-L-.12,2.11,L)
        if tweaks[i]: tweaks[i](ax)
        f(ax,l,True)
    _save(fig,name,sub)
def t_legy(y,one=False):
    def f(ax): ax._legy=y; ax._onecol=one
    return f
T4a=t_legy(-.23); T5a=t_legy(-.36); T5c=t_legy(-.26); T9b=t_legy(-.31); T6=t_legy(-.17)
panel(f4a,"Figure_4a","Main","a",T4a); panel(f4b,"Figure_4b","Main","b",t_legy(-.17)); comp([f4a,f4b],"Figure_4","Main","ab",[T4a,t_legy(-.17)])
panel(f5a,"Figure_5a","Main","a",T5a); panel(f5b,"Figure_5b","Main","b",t_legy(-.2,True)); panel(f5c,"Figure_5c","Main","c",T5c); comp([f5a,f5b,f5c],"Figure_5","Main","abc",[T5a,t_legy(-.2,True),T5c])
panel(f6a,"Figure_6a","Main","a",T6); panel(f6b,"Figure_6b","Main","b",T6)
fig=plt.figure(figsize=(3.5,4.7)); ax=_axes(fig,0,.95,2.4,3.3,1.0); ax._legy=-.1; ax._onecol=True; f6c(ax,"c",True); _save(fig,"Figure_6c","Main")
W=7.2; H2=CH+4.0; fig=plt.figure(figsize=(W,H2))
for i,(f,l) in enumerate(((f6a,"a"),(f6b,"b"))):
    ax=_axes(fig,i*W/2,H2-CH+1.12,W/2-.66-.12,2.11,.66); ax._legy=-.17; f(ax,l,True)
ax=_axes(fig,0,.7,W-1.2-.15,3.2,1.2); ax._wide=True; ax._legy=-.1; f6c(ax,"c",True); _save(fig,"Figure_6","Main")
panel(f7a,"Figure_7a","Main","a",t_legy(-.17)); panel(f7b,"Figure_7b","Main","b",t_legy(-.17)); comp([f7a,f7b],"Figure_7","Main","ab",[t_legy(-.17)]*2)
panel(f8a,"Figure_8a","Main","a",t_legy(-.17)); panel(f8b,"Figure_8b","Main","b",t_legy(-.17)); panel(f8c,"Figure_8c","Main","c",t_legy(-.17)); comp([f8a,f8b,f8c],"Figure_8","Main","abc",[t_legy(-.17),t_legy(-.17),t_legy(-.17)])
panel(f9a,"Figure_9a","Main","a",t_legy(-.17)); panel(f9b,"Figure_9b","Main","b",T9b); comp([f9a,f9b],"Figure_9","Main","ab",[t_legy(-.17),T9b])
panel(fS1a,"Figure_S1a","SI","a",t_legy(-.2)); panel(fS1b,"Figure_S1b","SI","b",t_legy(-.2)); comp([fS1a,fS1b],"Figure_S1","SI","ab",[t_legy(-.2)]*2)
panel(fS2,"Figure_S2","SI",None,t_legy(-.12,True)); panel(fS3,"Figure_S3","SI",None,t_legy(-.12,True))
panel(fS4a,"Figure_S4a","SI","a",None); panel(fS4b,"Figure_S4b","SI","b",t_legy(-.12)); comp([fS4a,fS4b],"Figure_S4","SI","ab",[None,t_legy(-.12)])
print("done")
