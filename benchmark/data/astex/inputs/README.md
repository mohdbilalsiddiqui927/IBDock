# Astex benchmark inputs (85 complexes)

Everything needed to rerun the Astex benchmark; docked poses and Vina logs are not included (see the Zenodo dataset).

- `raw_proteins/` original PDB files used as input.
- `raw_ligands/` crystallographic ligand of each complex, extracted from the PDB (HETATM).
- `prepared_ligands/` ligand files after IBDock preparation (PDBQT and protonated PDB).
- `prepared_receptors/` receptor PDBQT files after IBDock preparation (waters and non-selected ligands removed).
- `grid/<strategy>/<PDB>/grid_config.txt` Vina search box (centre and size, A) for each strategy:
  `co-crystallized`, `fpocket`, `p2rank`, `blind`. The box does not depend on seed or exhaustiveness
  (identical in all nine runs).
- `grid_source.txt` (method, ligand-selection mode, ligand) exists only for 1S19 and 1OYT, because earlier
  versions of the app did not write it. For the other 83 complexes only `grid_config.txt` was stored.

Known note: for 1JLA and 1MEH the stored co-crystal box centre equals the mean of the ligand together with a
modified residue (CSD A:280 and CSO A:319 respectively), not the ligand alone.
