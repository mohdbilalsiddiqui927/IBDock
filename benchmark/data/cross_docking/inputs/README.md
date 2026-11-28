# IBDock cross-docking set (LIT-PCBA templates): ESR1 (ago+ant pooled), MAPK1, PPARG

1,350 ordered pairs (receptor A, ligand B), 60 templates. Diagonal pairs (A = B) are cognate-redocking controls.
Primary analysis set: 1,270 pairs (`primary_analysis == True`).

## What to dock
For every row with `primary_analysis == True`:
- receptor: `receptor_file` (path inside the original LIT-PCBA folder). For ESR1 `5ufx` use `receptors_cleaned/ESR1_5ufx_protein_cleaned.mol2`
  (the original file contains 7 stray carbon atoms, group `<0>`, that coincide with its own ligand).
- ligand INPUT: `ligand_input_file` (crystal 3D conformer of ligand B in ITS OWN frame, hydrogens as in LIT-PCBA). Do not use the reference pose as input.
- strategies: the same four as the Astex run. For "Co-crystallized" the box comes from `box_ligand_file_for_co_crystallized`
  (receptor A's own ligand), never from ligand B. fpocket, P2Rank and blind do not use any ligand.
- settings: exhaustiveness 16, seed 42 (pilot); then seeds 1 and 123. Same receptor prep as the Astex pipeline.
- scoring: symmetry-corrected heavy-atom RMSD (sPyRMSD) between each docked pose and `reference_pose_file`
  (ligand B's crystal pose superposed into receptor A's frame). Success = RMSD <= 2.0 A; Top-1 and best-pose as in Table S3.
Dockings per (exhaustiveness, seed): 1,270 x 4 = 5,080 (all 1,350 pairs: 5,400).

## Rules fixed BEFORE docking (put these in the Methods)
1. Superposition: pairwise sequence alignment, then iterative C-alpha fit (5 cycles, residues > 3 A from the fit dropped). No residue-number matching.
2. Difficulty tier from `pocket_CA_rmsd` (C-alpha RMSD over aligned residues within 8 A of ligand A): easy < 1.5 A, medium 1.5-2.5 A, hard > 2.5 A. Diagonal = "cognate".
3. `noneq_site`: after superposition, the centroid of ligand B is > 6 A from the centroid of ligand A (different sub-site). Docked but excluded from the primary analysis; report separately. (ESR1 0, MAPK1 14, PPARG 36 pairs.)
4. Covalent ligands are excluded from the primary analysis (Vina is non-covalent): ESR1 `6chw` ligand (C20 bonded to Cys530 SG, 1.62 A) = 30 pairs.
5. Primary headline number = mean over all receptors (all pairs), never best-receptor selection. Report cognate vs cross-docking, block (within-ago, within-ant, cross-state) and tier.

## Ligand repair log (see ligand_qc.csv, column ligand_load)
50 of 60 ligand mol2 files read directly. 10 needed repair because RDKit could not read or kekulize them:
ESR1 5aau, 6chw; MAPK1 4qp3, 4qp4, 4qta, 5ax3, 5buj; PPARG 3r8a, 5y2t, 5z5s. Connectivity is taken from the mol2; bond orders and formal charges were re-perceived
from the 3D structure (net charge chosen by trial, so charge states are INFERRED, check them). 5ax3: the mol2 lacked the C-I bond, re-attached at 2.14 A.
Coordinates are never changed, so the crystal poses are exact.

## Columns of manifest.csv
pair_id, target, block (within / within-ago / within-ant / cross-state), receptor, ligand, cognate, aligned_residues, seq_identity, core_CA_fit,
pocket_CA_rmsd, centroid_mismatch_A, tier, noneq_site, primary_analysis, receptor_file, ligand_input_file, box_ligand_file_for_co_crystallized,
reference_pose_file, ligand_covalent_suspect, ligand_load, receptor_covalent_suspect.

## Not checked
Which PPARG/MAPK1 ligands bind which sub-site; PPARG ligands other than the covalent screen (distance only); protonation states of receptors; that a ligand is non-covalent beyond the 1.9 A distance screen.
Regenerate everything with build_manifest.py (needs the original LIT-PCBA folders).

## Batch files (batches/) - run in this order
Each row = one receptor-ligand job (one-to-one batch). Run every row with the 4 strategies, exhaustiveness 16, seed 42.
Paths: `ESR1_ago/`, `ESR1_ant/`, `MAPK1/`, `PPARG/` are the original LIT-PCBA target folders; `ligands/`, `receptors_cleaned/`, `reference_poses/` are in this package.
`ligand_file` (SDF, all 60 ligands) is always usable; `ligand_file_original_mol2` is filled only for the 50 ligands whose original mol2 was readable.
| file | pairs | dockings (x4 strategies) |
|---|---|---|
| stage1_ESR1_within-ago_within-ant.csv (pilot; 29 cognate controls) | 435 | 1,740 |
| stage2_ESR1_cross-state.csv | 435 | 1,740 |
| stage3_MAPK1.csv | 211 | 844 |
| stage4_PPARG.csv | 189 | 756 |
| stage5_optional_flagged-non-equivalent-site.csv (docked, analysed separately) | 50 | 200 |
Primary stages 1-4: 1,270 pairs, 5,080 dockings per (exhaustiveness, seed). Not docked at all: 30 pairs with the covalent 6chw ligand.
