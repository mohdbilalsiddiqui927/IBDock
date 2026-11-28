import pandas as pd, numpy as np
from rdkit import Chem, RDLogger; RDLogger.DisableLog("rdApp.*")
from rdkit.Chem import rdMolDescriptors
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
S="PATH/TO/WORKING_DIR/"  # set to your working directory (see benchmark/README.md)
A=pd.concat([pd.read_csv(S+"esr1_out/ESR1_all_strategies_scored.csv").assign(target="ESR1"),pd.read_csv(S+"mp_out/MAPK1_PPARG_all_scored.csv")],ignore_index=True)
P=A[A.primary_analysis]; C=P[~P.cognate]; pd.set_option("display.width",220)
print("TIER cocryst/fpocket/p2rank/blind  n,T1%"); print(C.groupby(["target","tier","strategy"]).agg(n=("t1","size"),t1=("t1","mean")).unstack(2).round(3).to_string())
print("\nESR1 block"); print(C[C.target=="ESR1"].groupby(["block","strategy"]).agg(n=("t1","size"),t1=("t1","mean")).unstack(1).round(3).to_string())
e=C[(C.target=="ESR1")&(C.strategy!="blind")]; print("pocketRMSD~top1_rmsd rho",round(spearmanr(e.pocket_CA_rmsd,e.top1_rmsd)[0],2))
rawpath={"ESR1":"esr1_raw/fd04b6c8-fpocket/fpocket/prepared_receptors/","MAPK1":"mp_raw/MAPK1/MAPK1/fpocket/prepared_receptors/","PPARG":"mp_raw/PPARG/PPARG/fpocket/prepared_receptors/"}
# ligand props
lp={}
for t in ["ESR1","MAPK1","PPARG"]:
    for l in C[C.target==t].ligand.unique():
        m=Chem.MolFromMolFile(S+f"cd_pkg/crossdock_out/ligands/{t}_{l}.sdf",removeHs=True,sanitize=False)
        try: Chem.SanitizeMol(m); rb=rdMolDescriptors.CalcNumRotatableBonds(m)
        except Exception: rb=np.nan
        lp[(t,l)]=(sum(1 for a in m.GetAtoms() if a.GetAtomicNum()>1),rb)
for t in ["ESR1","MAPK1","PPARG"]:
    pl=C[(C.target==t)&(C.strategy!="blind")].groupby("ligand").t1.mean(); df=pd.DataFrame({"t1":pl}); df["HA"]=[lp[(t,l)][0] for l in df.index]; df["RB"]=[lp[(t,l)][1] for l in df.index]
    print(t,"ligand-level Spearman T1~HA",round(spearmanr(df.t1,df.HA)[0],2),"T1~RB",round(spearmanr(df.t1,df.RB,nan_policy="omit")[0],2),"n lig",len(df),"| T1 range",df.t1.min().round(2),df.t1.max().round(2))
# clash
res=[]
for t in ["ESR1","MAPK1","PPARG"]:
    X=C[(C.target==t)&(C.strategy=="cocryst")]; tr={}
    for r in X.receptor.unique():
        pts=[[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in open(S+rawpath[t]+f"{r}.pdb") if l.startswith(("ATOM","HETATM")) and l[17:20].strip()!="LIG" and not l[12:16].strip().startswith("H") and l[76:78].strip()!="H"]; tr[r]=cKDTree(np.array(pts))
    for p,r in zip(X.pair_id,X.receptor):
        m=Chem.MolFromMolFile(S+f"cd_pkg/crossdock_out/reference_poses/{t}/{p}.sdf",removeHs=False,sanitize=False); xyz=m.GetConformer().GetPositions()[[a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum()>1]]
        res.append((t,p,tr[r].query(xyz)[0].min()))
cl=pd.DataFrame(res,columns=["target","pair_id","mind"]); cl["clash"]=cl.mind<2.5; cl.to_csv(S+"final/ref_clash.csv",index=False)
X=C[C.strategy=="cocryst"].merge(cl[["pair_id","mind","clash"]],on="pair_id")
print("\nCLASH (cross, cocryst box): frac and outcomes"); print(X.groupby(["target","clash"]).agg(n=("t1","size"),t1=("t1","mean"),bp=("bp","mean"),site=("site","mean")).round(3).to_string()); print(X.groupby("target").clash.mean().round(3).to_dict())
K=P[(P.cognate)]; Kc=K[K.strategy=="cocryst"].merge(cl[["pair_id"]],on="pair_id",how="left")
print("\nmedian Top-1 RMSD cross by target (cocryst):",C[C.strategy=="cocryst"].groupby("target").top1_rmsd.median().round(2).to_dict())
print("Cognate counts cocryst:",K[K.strategy=="cocryst"].groupby("target").agg(n=("t1","size"),t1=("t1","sum"),bp=("bp","sum")).to_dict())
print("n pocketRMSD tiers per target:",C[C.strategy=="cocryst"].groupby(["target","tier"]).size().to_dict())
