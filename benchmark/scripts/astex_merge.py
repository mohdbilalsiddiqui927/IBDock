import docx,pandas as pd,numpy as np,re
S="PATH/TO/WORKING_DIR/"  # set to your working directory (see benchmark/README.md)
M={"a":(16,42),"b":(16,1),"c":(16,123),"d":(8,42),"e":(8,1),"f":(8,123),"g":(32,42),"h":(32,1),"i":(32,123)}
P={"Co-crystallized":"cocryst","P2Rank":"p2rank","F-pocket":"fpocket","Blind":"blind"}
rows=[]
import glob
for f in glob.glob(S+"astex/s3/*.docx"):
    k=re.search(r"S3([a-i])_",f).group(1); ex,seed=M[k]
    t=docx.Document(f).tables[0]
    for r in list(t.rows)[1:]:
        c=[x.text.strip() for x in r.cells]
        if len(c)<8 or c[1] not in{x for x in []} and not re.match(r"^[0-9A-Z]{4}$",c[1]): continue
        num=lambda v: float(v) if re.match(r"^-?\d+(\.\d+)?$",v) else np.nan
        rows.append(dict(ex=ex,seed=seed,pdb=c[1],strategy=P[c[2]],score=num(c[3]),sp_top1=num(c[4]),sp_best=num(c[5]),br_top1=num(c[6]),br_best=num(c[7]),raw4=c[4]))
D=pd.DataFrame(rows); print(len(D)); 
N=pd.read_csv(S+"final/astex_naive_scored.csv"); X=D.merge(N,on=["ex","seed","strategy","pdb"],how="outer",indicator=True); print(X._merge.value_counts().to_dict())
X.to_csv(S+"final/astex_merged.csv",index=False)
print(D[D.sp_top1.isna()].raw4.value_counts().head())
