import numpy as np
from pathlib import Path
from rdkit import Chem, RDLogger; RDLogger.DisableLog("rdApp.*")
from spyrmsd.rmsd import symmrmsd
S=Path("PATH/TO/WORKING_DIR/"  # set to your working directory (see benchmark/README.md))
def models(p):
    ms=[];cur=None;sc=[]
    for ln in p.read_bytes().decode("latin1").splitlines():
        if ln.startswith("MODEL"): cur=[]
        elif ln.startswith("REMARK VINA RESULT"): sc.append(float(ln.split()[3]))
        elif ln.startswith(("ATOM","HETATM")) and cur is not None:
            el=ln[77:79].strip().upper() or ln[12:14].strip().upper()
            if el in("H","HD","HS"): continue
            el={"A":"C","OA":"O","NA":"N","SA":"S"}.get(el,el); cur.append((el,float(ln[30:38]),float(ln[38:46]),float(ln[46:54])))
        elif ln.startswith("ENDMDL"): ms.append(cur); cur=None
    return ms,sc
RAD={6:.76,7:.71,8:.66,16:1.05,9:.57,17:1.02}
def adj(z,xyz):
    d=np.linalg.norm(xyz[:,None]-xyz[None],axis=-1); r=np.array([RAD.get(int(a),.8) for a in z]); A=((d<r[:,None]+r[None]+.45)&(d>.1)).astype(int); np.fill_diagonal(A,0); return A
sym={"C":6,"N":7,"O":8,"S":16,"F":9,"CL":17}
import sys
ST=sys.argv[1]
for p in ["1IEP","1T46"]:
    ref=Chem.MolFromMolFile(str(S/f"rep_pkg/reference_poses_for_scoring/{p}_STI_crystal.sdf"),removeHs=False,sanitize=False)
    hv=[a.GetIdx() for a in ref.GetAtoms() if a.GetAtomicNum()>1]; rc=ref.GetConformer().GetPositions()[hv]; ra=np.array([ref.GetAtomWithIdx(i).GetAtomicNum() for i in hv]); am=Chem.GetAdjacencyMatrix(ref)[np.ix_(hv,hv)]
    ms,sc=models(S/f"rep_raw/{ST}/{ST}/docking/{p}_imatinib.pdbqt")
    rm=[]
    for md in ms:
        z=np.array([sym[e] for e,*_ in md]); xyz=np.array([[x,y,zz] for _,x,y,zz in md])
        rm.append(symmrmsd(rc,xyz,ra,z,am,adj(z,xyz),center=False,minimize=False))
    cd=[np.linalg.norm(np.array([[x,y,zz] for _,x,y,zz in md]).mean(0)-rc.mean(0)) for md in ms]
    print(p,"heavy ref",len(ra),"pose",len(ms[0]),"\n scores",sc,"\n rmsd",np.round(rm,2),"\n centroid d",np.round(cd,1))
