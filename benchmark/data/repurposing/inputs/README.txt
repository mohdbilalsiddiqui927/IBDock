IMATINIB RETROSPECTIVE REPURPOSING (ABL1 -> KIT; negative controls SRC, EGFR)
All 4 PDBs reduced to chain A (1IEP has chains A,B).
A_cocrystal_complex_pdbs/raw_proteins : protein + own ligand HETATM (1IEP/1T46 STI, 2SRC ANP, 1M17 AQ4) -> Co-crystallised run
B_protein_only_pdbs_fpocket_P2Rank_blind/raw_proteins : ligand removed -> fpocket, P2Rank, Blind runs
raw_ligands_same_for_A_and_B/imatinib.sdf : PubChem CID 5291 graph, re-embedded in 3D (not a crystal pose)
reference_poses_for_scoring/ : crystal STI of 1IEP and 1T46 (do NOT give to IBDock)
Jobs: 4 proteins x 1 ligand x 4 strategies = 16 dockings per (exhaustiveness, seed); 3 x 3 -> 144 dockings.
Settings same as Astex/cross-docking (padding too). Scoring done by me: imatinib vs crystal STI in 1T46 and 1IEP, <= 2 A.
Notes: 2SRC PTR kept as protein residue; PO4 and Cl in 1T46/1IEP removed; 2SRC is autoinhibited and 1M17 is DFG-in (weak-pocket negatives).
