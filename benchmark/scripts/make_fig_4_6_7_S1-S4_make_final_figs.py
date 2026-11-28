import numpy as np, pandas as pd, matplotlib, itertools
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Patch
from scipy import stats
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests
for f in fm.findSystemFonts():
    if "texgyrepagella" in f.lower(): fm.fontManager.addfont(f)
plt.rcParams.update({"font.family":"TeX Gyre Pagella","font.size":8,"axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":.9,
 "xtick.major.width":.9,"ytick.major.width":.9,"xtick.major.size":3,"ytick.major.size":3,"axes.unicode_minus":True,"hatch.linewidth":.9,
 "pdf.fonttype":42,"ps.fonttype":42,"legend.fontsize":7.5,"axes.labelsize":8.5,"xtick.labelsize":8,"ytick.labelsize":8})
S=["Co-crystallized","F-pocket","P2Rank","Blind"]; LAB=["Co-crystal","fpocket","P2Rank","Blind"]
COL=dict(zip(S,["#2b5aa0","#18512f","#550022","#554400"])); TINT=dict(zip(S,["#a0b4d6","#96b0a2","#b38c99","#b3aa88"])); GREY="#c5c9b9"; INK="#222222"
d=pd.read_pickle("s3_new.pkl"); PDBS=sorted(d.pdb.unique()); N=85
c=np.load("out/cubes.npz",allow_pickle=True); agg=pd.read_csv("out/agg.csv").query("endpoint_def=='primary'"); rates=pd.read_csv("out/rates.csv").query("endpoint_def=='primary'")
EX=[8,16,32]; SD=[1,42,123]
def cb(col):
    a=np.full((3,4,N),np.nan)
    for j,sd in enumerate(SD):
        g=d[(d.ex==16)&(d.seed==sd)]
        for k,s in enumerate(S): a[j,k]=g[g.proto==s].set_index("pdb").reindex(PDBS)[col].values
    return a
T1,BP,C1,CB,SC=cb("top1"),cb("best"),cb("cb_top1"),cb("cb_best"),cb("score")
def save(fig,name):
    fig.savefig("figs_final/"+name+".pdf",bbox_inches="tight"); fig.savefig("figs_final/"+name+".png",dpi=600,bbox_inches="tight"); plt.close(fig)
def pl(ax,t,x=-.13,y=1.04): ax.text(x,y,t,transform=ax.transAxes,fontsize=12,fontweight="bold")
def maj(b): return b.sum(0)>=2
# ===== Figure 4 =====
fig,(a1,a2)=plt.subplots(1,2,figsize=(7.2,3.5),gridspec_kw=dict(width_ratios=[1,1.05])); x=0
for k,s in enumerate(S):
    for ep,arr,nm in (("Top-1",T1,"Top-1"),("Best-pose",BP,"Best")):
        a=agg[(agg.endpoint==ep)&(agg.exhaustiveness==16)&(agg.strategy==s)].iloc[0]; un=np.isnan(arr[:,k,:]).mean()*100
        a1.bar(x,a.mean_pct,.9,color=COL[s],edgecolor=INK,lw=.8,zorder=2); a1.bar(x,100-un-a.mean_pct,.9,bottom=a.mean_pct,color=GREY,edgecolor=INK,lw=.8,zorder=2)
        a1.bar(x,un,.9,bottom=100-un,color="white",edgecolor="#555",hatch="//",lw=.8,zorder=2)
        a1.errorbar(x,a.mean_pct,yerr=[[a.mean_pct-a.boot_lo],[a.boot_hi-a.mean_pct]],color="black",capsize=2.5,lw=1.1,zorder=4)
        r=rates[(rates.endpoint==ep)&(rates.exhaustiveness==16)&(rates.strategy==s)].pct.values
        a1.scatter(x+np.linspace(-.22,.22,3),r,s=7,color="white",edgecolor="black",lw=.6,zorder=5)
        a1.text(x,a.mean_pct*.5,f"{a.mean_pct:.1f}",ha="center",va="center",fontsize=7.5,color="white",fontweight="bold",zorder=6); a1.text(x,-4,nm,ha="center",va="top",fontsize=7); x+=1
    a1.text(x-1.5,-14,LAB[k],ha="center",va="top",fontsize=8.5,fontweight="bold"); x+=.6
a1.set_ylim(0,100); a1.set_xlim(-.7,x-.3); a1.set_xticks([]); a1.set_ylabel("% of all 85 complexes",fontweight="bold"); pl(a1,"a",-.2)
M1=np.nanmedian(T1,axis=0); MB=np.nanmedian(BP,axis=0); ok=~np.isnan(T1).any(axis=(0,1))&~np.isnan(BP).any(axis=(0,1))
pos=[];dat=[];cols=[];x=0;rng=np.random.default_rng(3)
for k,s in enumerate(S):
    for arr,tint in ((M1,True),(MB,False)): dat.append(arr[k][ok]); pos.append(x); cols.append(TINT[s] if tint else COL[s]); x+=1
    x+=.8
bp=a2.boxplot(dat,positions=pos,widths=.62,patch_artist=True,showfliers=False,whis=(10,90),medianprops=dict(color="black",lw=1.4),whiskerprops=dict(lw=.9),capprops=dict(lw=.9),boxprops=dict(lw=.8,edgecolor="#333"))
for b,cc in zip(bp["boxes"],cols): b.set_facecolor(cc)
for p,v in zip(pos,dat): a2.scatter(rng.normal(p,.08,len(v)),v,s=3.5,color="black",alpha=.4,lw=0,zorder=1)
a2.axhline(2,color="black",ls="--",lw=1); a2.text((pos[1]+pos[2])/2,2.12,"2.0 Å",fontsize=7,va="bottom",ha="center")
a2.set_yscale("log"); a2.set_yticks([.2,.5,1,2,5,10,20,50]); a2.set_yticklabels(["0.2","0.5","1","2","5","10","20","50"]); a2.minorticks_off()
a2.set_xticks(pos); a2.set_xticklabels(["Top-1","Best"]*4,fontsize=7)
for k in range(4): a2.text((pos[2*k]+pos[2*k+1])/2,-.1,LAB[k],transform=a2.get_xaxis_transform(),ha="center",va="top",fontsize=8.5,fontweight="bold")
a2.set_ylabel("sPyRMSD (Å), log scale",fontweight="bold"); pl(a2,"b",-.17); N4B=int(ok.sum())
fig.legend(handles=[Patch(fc="white",ec="black",label="Success (RMSD < 2.0 Å)"),Patch(fc=GREY,ec="black",label="Resolved, not success"),Patch(fc="white",ec="#555",hatch="//",label="sPyRMSD unresolved")],loc="lower left",ncol=3,frameon=False,bbox_to_anchor=(.01,-.1))
fig.tight_layout(w_pad=2); save(fig,"Figure4_astex_performance")
# ===== Figure 5 =====
fig,axs=plt.subplots(1,2,figsize=(7.2,3.7))
for ax,(A,B_,ttl,lab) in zip(axs,((T1,C1,"Top-1 pose","a"),(BP,CB,"Best pose","b"))):
    xs=[];ys=[]
    for k,s in enumerate(S):
        m=~np.isnan(A[:,k,:])&~np.isnan(B_[:,k,:]); xs+=list(A[:,k,:][m]); ys+=list(B_[:,k,:][m]); ax.scatter(A[:,k,:][m],B_[:,k,:][m],s=14,color=COL[s],edgecolor="black",lw=.3,alpha=.8,label=LAB[k],zorder=3)
    xs,ys=np.array(xs),np.array(ys); rho=stats.spearmanr(xs,ys)[0]; mae=np.abs(xs-ys).mean(); lim=[.15,60]
    ax.plot(lim,lim,color="#444",lw=1,label="Identity (y = x)",zorder=2); ax.axvline(2,color="black",ls="--",lw=.8); ax.axhline(2,color="black",ls="--",lw=.8)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim); tk=[.2,.5,1,2,5,10,20,50]
    for f_ in (ax.set_xticks,ax.set_yticks): f_(tk)
    ax.set_xticklabels([str(t) for t in tk]); ax.set_yticklabels([str(t) for t in tk]); ax.minorticks_off(); ax.set_xlabel("sPyRMSD (Å)",fontweight="bold"); ax.set_ylabel("Coordinate-bridge RMSD (Å)",fontweight="bold")
    ax.text(.04,.96,f"n = {len(xs)}\nρ = {rho:.3f}\nMAE = {mae:.3f} Å",transform=ax.transAxes,va="top",fontsize=7.5); ax.set_title(ttl,fontweight="bold",fontsize=9); pl(ax,lab,-.2)
    if lab=="a": ax.legend(loc="lower right",frameon=False,fontsize=7); print("fig5",lab,len(xs),round(rho,3),round(mae,3))
fig.tight_layout(); save(fig,"Figure5_bridge_agreement")
# ===== Figure 6 new: grid =====
fig=plt.figure(figsize=(7.2,5.6)); gs=fig.add_gridspec(2,4,height_ratios=[1.2,1],hspace=.55,wspace=.25)
for i,(ep,lab) in enumerate((("Top-1","a"),("Best-pose","b"))):
    ax=fig.add_subplot(gs[0,2*i:2*i+2])
    for k,s in enumerate(S):
        a=agg[(agg.endpoint==ep)&(agg.strategy==s)].sort_values("exhaustiveness"); ax.plot(range(3),a.mean_pct,color=COL[s],lw=1.6,zorder=2)
        for j,e in enumerate(EX): ax.scatter(np.full(3,j)+np.linspace(-.07,.07,3),rates[(rates.endpoint==ep)&(rates.strategy==s)&(rates.exhaustiveness==e)].pct.values,s=9,color="white",edgecolor=COL[s],lw=.8,zorder=3)
        ax.scatter(range(3),a.mean_pct,s=26,color=COL[s],edgecolor="black",lw=.7,zorder=4)
        off={"Top-1":[0,0,2.5,-2.5],"Best-pose":[0,0,2.8,-2.8]}[ep][k]; ax.text(2.13,a.mean_pct.values[-1]+off,LAB[k],va="center",fontsize=7.5,fontweight="bold")
    ax.set_xticks(range(3)); ax.set_xticklabels(EX); ax.set_xlim(-.25,3.15); ax.set_ylim(0,80); ax.set_xlabel("Vina exhaustiveness"); ax.set_title(ep,loc="left",fontweight="bold",fontsize=9); ax.yaxis.grid(True,color="#e4e3df",lw=.6); ax.set_axisbelow(True); pl(ax,lab,-.2,1.06)
    if i==0: ax.set_ylabel("Success, RMSD < 2.0 Å (% of 85)",fontweight="bold")
a_=c["top1"]; CS,CF,CN="#2b5aa0","#d9822b","#c5c9b9"; axs=[fig.add_subplot(gs[1,k]) for k in range(4)]
for k,(ax,s) in enumerate(zip(axs,S)):
    for i,e in enumerate(EX):
        t=a_[i,:,k,:].sum(0); n3,n0=(t==3).sum(),(t==0).sum(); nf=85-n3-n0; left=0
        for v,cc in ((n3,CS),(nf,CF),(n0,CN)):
            ax.barh(i,v,left=left,color=cc,edgecolor="white",lw=1,height=.7)
            if v>=8: ax.text(left+v/2,i,str(v),ha="center",va="center",fontsize=7,color="white" if cc!=CN else INK)
            left+=v
    ax.set_yticks(range(3)); ax.set_yticklabels([f"ex {e}" for e in EX] if k==0 else []); ax.invert_yaxis(); ax.set_xlim(0,85); ax.set_xticks([0,40,85]); ax.set_title(LAB[k],fontweight="bold",loc="left",fontsize=8.5)
    if k==0: pl(ax,"c",-.55,1.2)
fig.text(.5,.058,"Complexes (n = 85), Top-1 outcome",ha="center",fontsize=8)
fig.legend(handles=[Patch(color=CS,label="success in all 3 seeds"),Patch(color=CF,label="outcome flips between seeds"),Patch(color=CN,label="failure in all 3 seeds")],loc="lower center",ncol=3,frameon=False,bbox_to_anchor=(.5,-.04))
save(fig,"Figure_exhaustiveness_seed")
# ===== SI S1 =====
fig,axs=plt.subplots(1,2,figsize=(7.2,3.3),sharey=True)
for ax,(arr,ttl,lab) in zip(axs,((T1,"Top-1","a"),(BP,"Best-pose","b"))):
    m=np.nanmedian(arr,axis=0); xg=np.logspace(np.log10(.15),np.log10(80),400)
    for k,s in enumerate(S): ax.plot(xg,[(m[k]<=g).sum()/N*100 for g in xg],color=COL[s],lw=1.6,label=LAB[k]); ax.scatter([2],[(m[k]<=2).sum()/N*100],s=24,color=COL[s],edgecolor="black",zorder=5)
    ax.axvline(2,color="black",ls="--",lw=.8); ax.set_xscale("log"); ax.set_xticks([.2,.5,1,2,5,10,20,50]); ax.set_xticklabels(["0.2","0.5","1","2","5","10","20","50"]); ax.minorticks_off(); ax.set_ylim(0,100)
    ax.set_xlabel(f"{ttl} sPyRMSD (Å)",fontweight="bold"); ax.yaxis.grid(True,color="#e4e3df",lw=.6); ax.set_axisbelow(True); pl(ax,lab,-.2)
    if lab=="a": ax.set_ylabel("% of all 85 complexes below threshold",fontweight="bold"); ax.legend(frameon=False,loc="upper left")
fig.tight_layout(); save(fig,"SI_Figures/FigureS1_ecdf")
# ===== SI S2 + S3 =====
def stack(ax,rows,legend):
    for k,(vals,s) in enumerate(zip(rows,S)):
        b=0
        for v,cc,h,tc in zip(vals,[COL[s],TINT[s],GREY,"white"],[None,None,None,"//"],["white",INK,INK,INK]):
            ax.bar(k,v/N*100,.62,bottom=b,color=cc,edgecolor="#555" if h else INK,hatch=h,lw=.8)
            if v>0: ax.text(k,b+v/N*50,str(v),ha="center",va="center",fontsize=8,color=tc,bbox=dict(fc="white",ec="none",alpha=.8,pad=.5) if h else None)
            b+=v/N*100
    ax.set_xticks(range(4)); ax.set_xticklabels(LAB,fontweight="bold"); ax.set_ylim(0,100); ax.set_ylabel("% of all 85 complexes",fontweight="bold")
rows2=[];rows3=[]
for k,s in enumerate(S):
    un=maj(np.isnan(T1[:,k,:])); t1=maj(np.nan_to_num(T1[:,k,:],nan=99)<=2)&~un; bpm=maj(np.nan_to_num(BP[:,k,:],nan=99)<=2)&~t1&~un
    rows2.append((int(t1.sum()),int(bpm.sum()),N-int(t1.sum())-int(bpm.sum())-int(un.sum()),int(un.sum())))
    sp=maj(~np.isnan(T1[:,k,:])); br=maj(~np.isnan(C1[:,k,:])); rows3.append((int((sp&br).sum()),int((sp&~br).sum()),int((~sp&br).sum()),N-int((sp|br).sum())))
fig,ax=plt.subplots(figsize=(4.2,4.2)); stack(ax,rows2,None)
fig.legend(handles=[Patch(fc="#555",ec="black",label="Top-1 pose < 2.0 Å (strategy color)"),Patch(fc="#dddddd",ec="black",label="Near-native pose generated, not ranked first (light tint)"),Patch(fc=GREY,ec="black",label="No pose < 2.0 Å generated"),Patch(fc="white",ec="#555",hatch="//",label="sPyRMSD unresolved")],loc="lower center",frameon=False,bbox_to_anchor=(.5,-.2))
save(fig,"SI_Figures/FigureS2_failure_modes")
fig,ax=plt.subplots(figsize=(4.2,4.2)); stack(ax,rows3,None)
fig.legend(handles=[Patch(fc="#555",ec="black",label="Both routes numerical (strategy color)"),Patch(fc="#dddddd",ec="black",label="sPyRMSD only (light tint)"),Patch(fc=GREY,ec="black",label="Coordinate bridge only"),Patch(fc="white",ec="#555",hatch="//",label="Neither route")],loc="lower center",frameon=False,bbox_to_anchor=(.5,-.2))
save(fig,"SI_Figures/FigureS3_correspondence_coverage"); print("S2",rows2,"S3",rows3)
# ===== SI S4 =====
okc=~np.isnan(T1).any(axis=(0,1)); sm=np.nanmean(SC,axis=0)
fig,(a1,a2)=plt.subplots(1,2,figsize=(7.2,3.7))
pt=a1.violinplot([sm[k][okc] for k in range(4)],positions=range(4),showextrema=False,widths=.8)
for p_,s in zip(pt["bodies"],S): p_.set_facecolor(COL[s]); p_.set_alpha(.6); p_.set_edgecolor("#333")
for k in range(4): q1,m_,q3=np.percentile(sm[k][okc],[25,50,75]); a1.vlines(k,q1,q3,color="black",lw=3,zorder=3); a1.scatter(k,m_,s=16,color="white",edgecolor="black",zorder=4)
pw=list(itertools.combinations(range(4),2)); ph=multipletests([stats.wilcoxon(sm[i][okc],sm[j][okc]).pvalue for i,j in pw],method="holm")[1]; top=sm[:,okc].max(); h=0
for (i,j),pp in sorted(zip(pw,ph),key=lambda z:z[0][1]-z[0][0]):
    if pp<.05: h+=.9; y=top+.3+h; a1.plot([i,i,j,j],[y-.2,y,y,y-.2],color="black",lw=.8); a1.text((i+j)/2,y+.05,"***" if pp<1e-3 else "**" if pp<1e-2 else "*",ha="center",va="bottom",fontsize=8)
a1.set_xticks(range(4)); a1.set_xticklabels(LAB,fontweight="bold"); a1.set_ylabel("Top-1 Vina score (kcal/mol)",fontweight="bold"); a1.set_title(f"Docking score by strategy (n = {okc.sum()})",fontweight="bold",fontsize=9); pl(a1,"a",-.2)
for k,s in enumerate(S):
    rec=maj(np.nan_to_num(T1[:,k,:],nan=99)<=2)[okc]; v=sm[k][okc]
    for off,msk,fc in ((-.2,rec,COL[s]),(.2,~rec,"#d6d6d6")):
        b=a2.violinplot([v[msk]],positions=[k+off],showextrema=False,widths=.38)["bodies"][0]; b.set_facecolor(fc); b.set_alpha(.65 if fc!="#d6d6d6" else .9); b.set_edgecolor("#555")
        q1,m_,q3=np.percentile(v[msk],[25,50,75]); a2.vlines(k+off,q1,q3,color="black",lw=3,zorder=3); a2.scatter(k+off,m_,s=14,color="white",edgecolor="black",zorder=4)
    auc=roc_auc_score(rec.astype(int),-v); r2=np.random.default_rng(5); bs=[]
    for _ in range(5000):
        ix=r2.integers(0,len(v),len(v))
        if 0<rec[ix].sum()<len(ix): bs.append(roc_auc_score(rec[ix].astype(int),-v[ix]))
    lo,hi=np.percentile(bs,[2.5,97.5]); a2.text(k,-1.7,f"AUC {auc:.2f}\n[{lo:.2f}, {hi:.2f}]",ha="center",fontsize=7.5); print("auc",S[k],round(auc,2),round(lo,2),round(hi,2))
a2.set_xticks(range(4)); a2.set_xticklabels(LAB,fontweight="bold"); a2.set_ylim(-14.8,-1.0); a2.set_ylabel("Top-1 Vina score (kcal/mol)",fontweight="bold"); a2.set_title("Docking score and pose recovery",fontweight="bold",fontsize=9); pl(a2,"b",-.2)
a2.legend(handles=[Patch(fc="#888",ec="#555",label="Recovered (strategy color)"),Patch(fc="#d6d6d6",ec="#555",label="Not recovered")],loc="lower right",frameon=False,fontsize=7)
fig.tight_layout(); save(fig,"SI_Figures/FigureS4_score_and_recovery"); print("n4b",N4B)
