"""Cross-docking manifest for LIT-PCBA templates (ESR1 ago+ant pooled, MAPK1, PPARG).
For every ordered pair (receptor A, ligand B) of the same target:
  * B is superposed onto A by sequence alignment + iterative core C-alpha fit (3 A trimming, 5 cycles);
  * the crystal pose of B is rewritten in A's frame (reference pose for RMSD only);
  * the docking INPUT stays ligand B in its own crystal frame (no pose leak).
Rules fixed before docking: tier by pocket C-alpha RMSD (<1.5 easy, 1.5-2.5 medium, >2.5 hard);
pairs whose superposed centroid mismatch > 6 A are flagged non-equivalent site and excluded from the primary analysis.
"""
import os, glob, sys, collections, numpy as np, pandas as pd
from Bio.Align import PairwiseAligner
from Bio.Data.IUPACData import protein_letters_3to1 as t31
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors, rdDetermineBonds
RDLogger.DisableLog("rdApp.*")
OUT = "crossdock_out"; os.makedirs(f"{OUT}/ligands", exist_ok=True); os.makedirs(f"{OUT}/reference_poses", exist_ok=True)
SETS = {"ESR1": [("esr_ago/ESR1_ago", "ago"), ("esr_ant/ESR1_ant", "ant")], "MAPK1": [("mapk1/MAPK1", "")], "PPARG": [("pparg/PPARG", "")]}
CENT_CUT, TIER = 6.0, (1.5, 2.5)

def atoms(path):
    out = []; sec = None
    for ln in open(path):
        if ln.startswith("@<TRIPOS>"): sec = ln.strip(); continue
        if sec == "@<TRIPOS>ATOM" and ln.strip():
            p = ln.split(); out.append((p[1], np.array(p[2:5], float), p[5], p[7] if len(p) > 7 else ""))
    return out

def manual_graph(path):
    """Own mol2 reader (RDKit rejects some SYBYL types such as C.cat): elements + coordinates + connectivity, all single bonds."""
    at, bd, sec = [], [], None
    for ln in open(path):
        if ln.startswith("@<TRIPOS>"): sec = ln.strip(); continue
        if not ln.strip(): continue
        p = ln.split()
        if sec == "@<TRIPOS>ATOM": at.append((p[5].split(".")[0], [float(x) for x in p[2:5]]))
        elif sec == "@<TRIPOS>BOND": bd.append((int(p[1]) - 1, int(p[2]) - 1))
    m = Chem.RWMol(); conf = Chem.Conformer(len(at))
    for i, (el, xyz) in enumerate(at):
        a = Chem.Atom(el.capitalize() if len(el) > 1 else el); a.SetNoImplicit(True); m.AddAtom(a); conf.SetAtomPosition(i, xyz)
    for i, j in bd: m.AddBond(i, j, Chem.BondType.SINGLE)
    m = m.GetMol(); m.AddConformer(conf); return m

HAL_MAX = {"I": 2.5, "Br": 2.2, "Cl": 2.0, "F": 1.6}
def join_halogens(raw):
    """a lone halogen atom detached from the rest of the ligand (missing bond in the mol2) is re-attached to its nearest heavy atom."""
    rw = Chem.RWMol(raw); conf = rw.GetConformer(); added = []
    frags = Chem.GetMolFrags(rw)
    if len(frags) < 2: return raw, added
    for fr in frags:
        if len(fr) == 1 and rw.GetAtomWithIdx(fr[0]).GetSymbol() in HAL_MAX:
            i = fr[0]; el = rw.GetAtomWithIdx(i).GetSymbol(); xi = np.array(list(conf.GetAtomPosition(i)))
            best = min(((np.linalg.norm(xi - np.array(list(conf.GetAtomPosition(a.GetIdx())))), a.GetIdx()) for a in rw.GetAtoms()
                        if a.GetIdx() != i and a.GetAtomicNum() > 1 and a.GetIdx() not in fr))
            if best[0] <= HAL_MAX[el]: rw.AddBond(i, best[1], Chem.BondType.SINGLE); added.append((el, best[1], round(best[0], 2)))
    return rw.GetMol(), added

def load_ligand(path):
    """returns (mol, method). mol2 via RDKit first; otherwise (or if detached fragments) bond orders re-perceived on the mol2 connectivity."""
    m = Chem.MolFromMol2File(path, removeHs=False)
    if m is not None and len(Chem.GetMolFrags(m)) == 1: return m, "mol2"
    raw = Chem.MolFromMol2File(path, removeHs=False, sanitize=False)
    src = "repaired"
    if raw is None or m is not None: raw = manual_graph(path); src = "rebuilt"
    raw, added = join_halogens(raw)
    if added: src += f"+halogen-bond{added}"
    for ch in (0, 1, -1, 2, -2):
        mm = Chem.RWMol(raw)
        for b in mm.GetBonds(): b.SetBondType(Chem.BondType.SINGLE); b.SetIsAromatic(False)
        for a in mm.GetAtoms(): a.SetIsAromatic(False); a.SetFormalCharge(0); a.SetNoImplicit(True); a.SetNumExplicitHs(0)
        mm = mm.GetMol()
        try:
            rdDetermineBonds.DetermineBondOrders(mm, charge=ch, allowChargedFragments=True)
            Chem.SanitizeMol(mm)
            ref = {frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx())) for b in raw.GetBonds()}
            new = {frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx())) for b in mm.GetBonds()}
            if ref == new: return mm, f"{src}(charge={ch})"
        except Exception:
            continue
    return None, "FAILED"

al = PairwiseAligner(); al.mode = "global"; al.match_score = 2; al.mismatch_score = -1; al.open_gap_score = -4; al.extend_gap_score = -0.5
def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0); U, S, Vt = np.linalg.svd((P - Pc).T @ (Q - Qc)); d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T; return R, Qc - R @ Pc

def strip_receptor(src, dst, subst_name):
    """write a copy of a protein mol2 without the atoms of one substructure name (bonds and ids renumbered)."""
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    lines = open(src).read().split("\n"); out = []; sec = None; keep = {}; newid = 0; atom_rows = []; bond_rows = []; head = []; tail = []
    for ln in lines:
        if ln.startswith("@<TRIPOS>"): sec = ln.strip()
        if sec == "@<TRIPOS>ATOM" and not ln.startswith("@") and ln.strip():
            p = ln.split()
            if p[7] == subst_name: continue
            newid += 1; keep[int(p[0])] = newid; p[0] = str(newid); atom_rows.append("  ".join(p)); continue
        if sec == "@<TRIPOS>BOND" and not ln.startswith("@") and ln.strip():
            p = ln.split()
            if int(p[1]) in keep and int(p[2]) in keep: bond_rows.append(p); continue
            continue
        if sec == "@<TRIPOS>SUBSTRUCTURE" and not ln.startswith("@") and ln.strip():
            p = ln.split()
            if p[1] == subst_name or int(p[2]) not in keep: continue
            p[2] = str(keep[int(p[2])]); out.append("  ".join(p)); continue
        out.append(ln)
    res = []; i = 0
    for ln in out:
        res.append(ln)
        if ln.startswith("@<TRIPOS>ATOM"): res += atom_rows
        if ln.startswith("@<TRIPOS>BOND"): res += [f"{k+1} {b[1]} {b[2]} {b[3]}" for k, b in enumerate(bond_rows)] and [f"     {k+1} {keep[int(b[1])]} {keep[int(b[2])]} {b[3]}" for k, b in enumerate(bond_rows)]
    txt = "\n".join(res)
    # fix molecule counts line (atoms bonds substructures)
    import re
    nsub = 0
    if "@<TRIPOS>SUBSTRUCTURE" in res:
        for ln in res[res.index("@<TRIPOS>SUBSTRUCTURE") + 1:]:
            if ln.startswith("@"): break
            if ln.strip(): nsub += 1
    k = res.index("@<TRIPOS>MOLECULE"); cnt = res[k + 2].split(); cnt[0], cnt[1], cnt[2] = str(newid), str(len(bond_rows)), str(nsub); res[k + 2] = " ".join(cnt)
    open(dst, "w").write("\n".join(res))

rows, qc = [], []
flag_cov = {}
for target, folders in SETS.items():
    T = {}
    for folder, tag in folders:
        for f in sorted(glob.glob(folder + "/*_protein.mol2")):
            pid = os.path.basename(f)[:4]
            pa = atoms(f); lp = f.replace("_protein", "_ligand")
            mol, method = load_ligand(lp)
            if mol is None: print("LIGAND FAILED", target, pid); continue
            conf = mol.GetConformer(); X = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
            heavy_idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
            ca = [(s, xyz) for n, xyz, t, s in pa if n == "CA"]
            seq = "".join(t31.get(s[:3].capitalize(), "X") for s, _ in ca)
            stray = [a for a in pa if a[3] == "<0>"]
            prot_heavy = np.array([a[1] for a in pa if not a[2].startswith("H") and a[3] != "<0>"])
            mind = min(np.linalg.norm(prot_heavy - X[i], axis=1).min() for i in heavy_idx)
            frags = len(Chem.GetMolFrags(Chem.RemoveHs(mol)))
            charge = Chem.GetFormalCharge(mol)
            flag_cov[(target, pid)] = bool(mind < 1.9)
            T[pid] = dict(mol=mol, X=X, heavy=heavy_idx, seq=seq, CA=np.array([x for _, x in ca]), tag=tag, folder=folder, method=method)
            hm = Chem.RemoveHs(mol)
            qc.append(dict(target=target, template=pid, state=tag, ligand_load=method, heavy_atoms=hm.GetNumHeavyAtoms(), MW=round(Descriptors.MolWt(hm), 1),
                           rot_bonds=rdMolDescriptors.CalcNumRotatableBonds(hm), net_charge=charge, fragments=frags,
                           min_ligand_protein_heavy_dist=round(float(mind), 2), covalent_suspect=bool(mind < 1.9), n_CA=len(ca),
                           receptor_stray_atoms_removed=len(stray)))
            if stray: strip_receptor(f, f"{OUT}/receptors_cleaned/{target}_{pid}_protein_cleaned.mol2", "<0>")
            sdf = f"{OUT}/ligands/{target}_{pid}.sdf"; w = Chem.SDWriter(sdf); mol.SetProp("_Name", f"{target}_{pid}"); w.write(mol); w.close()
    ids = sorted(T)
    cache = {}
    for ri in ids:          # receptor A
        for lj in ids:      # ligand B (from template lj)
            A, B = T[ri], T[lj]
            if ri == lj:
                R, t = np.eye(3), np.zeros(3); n_al = len(A["seq"]); ident = 1.0; core = n_al; prm = 0.0; cd = 0.0
            else:
                aln = al.align(A["seq"], B["seq"])[0]; ia, ib = [], []
                for (s1, e1), (s2, e2) in zip(*aln.aligned): ia += list(range(s1, e1)); ib += list(range(s2, e2))
                ia, ib = np.array(ia), np.array(ib); ident = float(np.mean([A["seq"][a] == B["seq"][b] for a, b in zip(ia, ib)]))
                Q, P = A["CA"][ia], B["CA"][ib]; keep = np.ones(len(P), bool)
                for _ in range(5):
                    if keep.sum() < 10: break
                    R, t = kabsch(P[keep], Q[keep]); keep = np.linalg.norm(P @ R.T + t - Q, axis=1) < 3.0
                R, t = kabsch(P[keep], Q[keep]); d = np.linalg.norm(P @ R.T + t - Q, axis=1)
                LA = A["X"][A["heavy"]]
                pocket = [k for k, a in enumerate(ia) if np.min(np.linalg.norm(LA - A["CA"][a], axis=1)) <= 8.0]
                prm = float(np.sqrt((d[pocket] ** 2).mean())) if len(pocket) >= 8 else np.nan
                n_al, core = len(ia), int(keep.sum())
                cd = float(np.linalg.norm((B["X"][B["heavy"]] @ R.T + t).mean(0) - LA.mean(0)))
            # reference pose of ligand B in receptor-A frame
            mB = Chem.Mol(B["mol"]); conf = mB.GetConformer()
            for i in range(mB.GetNumAtoms()):
                conf.SetAtomPosition(i, (B["X"][i] @ R.T + t).tolist())
            os.makedirs(f"{OUT}/reference_poses/{target}", exist_ok=True)
            ref = f"reference_poses/{target}/{target}_rec-{ri}__lig-{lj}.sdf"
            w = Chem.SDWriter(f"{OUT}/{ref}"); mB.SetProp("_Name", f"{target}_rec-{ri}_lig-{lj}"); w.write(mB); w.close()
            if ri == lj: block = f"within-{A['tag']}" if A["tag"] else "within"; tier = "cognate"
            else:
                block = (f"within-{A['tag']}" if A["tag"] == B["tag"] else "cross-state") if A["tag"] else "within"
                tier = "easy" if prm < TIER[0] else "medium" if prm < TIER[1] else "hard"
            noneq = bool(cd > CENT_CUT)
            rows.append(dict(pair_id=f"{target}_rec-{ri}__lig-{lj}", target=target, block=block, receptor=ri, receptor_state=A["tag"], ligand=lj, ligand_state=B["tag"],
                             cognate=(ri == lj), aligned_residues=n_al, seq_identity=round(ident, 3), core_CA_fit=core, pocket_CA_rmsd=round(prm, 2) if prm == prm else np.nan,
                             centroid_mismatch_A=round(cd, 2), tier=tier, noneq_site=noneq, primary_analysis=bool((not noneq) and not flag_cov.get((target, lj), False)),
                             receptor_file=f"{A['folder']}/{ri}_protein.mol2", ligand_input_file=f"ligands/{target}_{lj}.sdf",
                             box_ligand_file_for_co_crystallized=f"{A['folder']}/{ri}_ligand.mol2", reference_pose_file=ref))

man = pd.DataFrame(rows); q = pd.DataFrame(qc)
# ligand-level flags into the manifest
flag = q.set_index(["target", "template"])
man["ligand_covalent_suspect"] = [flag.loc[(t, l), "covalent_suspect"] for t, l in zip(man.target, man.ligand)]
man["ligand_load"] = [flag.loc[(t, l), "ligand_load"] for t, l in zip(man.target, man.ligand)]
man["receptor_covalent_suspect"] = [flag.loc[(t, r), "covalent_suspect"] for t, r in zip(man.target, man.receptor)]
man.to_csv(f"{OUT}/manifest.csv", index=False); q.to_csv(f"{OUT}/ligand_qc.csv", index=False)
print("pairs", len(man)); print(man.groupby(["target", "block"]).size())
print(man[~man.cognate].groupby(["target", "tier"]).size().unstack(fill_value=0))
print("noneq_site flagged:", man.groupby("target").noneq_site.sum().to_dict(), "| total", int(man.noneq_site.sum()))
print(q.ligand_load.value_counts().to_dict()); print(q[q.ligand_load != "mol2"][["target", "template", "ligand_load", "net_charge"]])
print("covalent suspects:", q[q.covalent_suspect][["target", "template", "min_ligand_protein_heavy_dist"]].to_dict("records"))
print("fragments>1:", q[q.fragments > 1][["target", "template", "fragments"]].to_dict("records"))
print("NaN pocket rmsd (non-cognate):", int(man[~man.cognate].pocket_CA_rmsd.isna().sum()))
