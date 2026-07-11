# IBDock

[![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20macOS%20%7C%20Linux-lightgrey)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

IBDock is a browser-based GUI for running AutoDock Vina docking jobs without touching the command line. It handles everything in one place — protein and ligand preparation, grid box generation, batch docking across multiple receptor–ligand pairs, 3D pose visualisation, and re-docking RMSD validation — and exports a PDF report when you are done.

It was built because setting up Vina manually is fragmented across too many tools, and because existing GUIs either cap you at one ligand at a time or break silently on real-world IUPAC chemical names.

![IBDock Results tab showing docking scores, affinity heatmap, and per-ligand distributions for six protein-ligand pairs](docs/screenshots/tab4_results.png)

---

## Features

- **Batch docking** — dock any number of ligands against one or more receptors in a single session, running all jobs in parallel via `concurrent.futures`
- **Automated grid box derivation** — five selectable modes covering co-crystallised ligand extraction, P2Rank, fpocket, blind docking, and manual entry; the Auto mode cascades through them automatically
- **Built-in RMSD validation** — symmetric heavy-atom RMSD against crystal reference poses, with Excellent / Pass / Borderline / Fail classification following Warren et al. (2006)
- **Interactive 3D pose viewer** — powered by 3Dmol.js, embedded directly in the browser; no PyMOL or external viewer required
- **PDF report export** — formatted multi-page validation report generated automatically from your results
- **Cross-platform** — runs on Windows, macOS, and Linux
- **No command-line required** — the entire workflow is GUI-driven from a browser tab

---

## Requirements

**Python 3.9 or later.**

```bash
pip install -r requirements.txt
```

Four external tools must be installed separately:

| Tool | Version | Download |
|------|---------|----------|
| MGLTools | 1.5.7 | https://ccsb.scripps.edu/mgltools/downloads/ |
| AutoDock Vina | ≥ 1.2.0 | https://github.com/ccsb-scripps/AutoDock-Vina/releases |
| Open Babel | ≥ 3.1.0 | https://openbabel.org/wiki/Category:Installation |
| P2Rank | 2.4 | https://github.com/rdk/p2rank/releases *(optional)* |
| fpocket | ≥ 4.0 | https://github.com/Discngine/fpocket *(optional)* |

P2Rank and fpocket are only needed if you want automated pocket prediction. On Windows,
both must be run through WSL — IBDock handles this automatically if you prefix the
executable paths with `wsl` in the Settings panel.

---

## Installation

```bash
git clone https://github.com/BilalSiddiqui/ibdock.git
cd ibdock
pip install -r requirements.txt
streamlit run IBDock.py
```

This opens IBDock in your browser at `http://localhost:8501`. On first launch, open
the **⚙ Settings** panel in the sidebar and point IBDock at your MGLTools, Vina, and
Open Babel executables. Paths are saved to `IBDock_config.json` and remembered across
sessions.

---

## Quick Start

The `example/` folder contains a prepared 3OCB system (AKT1 kinase) ready to dock:

```
example/
├── receptor/   3OCB.pdb                  raw PDB file as downloaded from RCSB
├── ligand/     3OCB_ligand.sdf           co-crystallised ligand in SDF format
└── reference/  crystal_pose_3OCB.pdb     crystal pose for RMSD validation
```

1. Upload `3OCB.pdb` in **Protein Prep** — the grid box is detected automatically
2. Upload `3OCB_ligand.sdf` in **Ligand Prep**
3. Click **Run Docking** in the Docking tab
4. Upload `crystal_pose_3OCB.pdb` in **Validation** and compute RMSD

**Expected result:** affinity ≈ −9.2 kcal/mol · RMSD ≈ 0.648 Å (Excellent)

---

## Workflow

IBDock is organised around seven sequential tabs:

| Tab | What it does |
|-----|-------------|
| **Protein Prep** | Upload PDB, remove waters/HETATM, compute grid box, generate PDBQT |
| **Ligand Prep** | Convert SDF / MOL2 / PDB to PDBQT with Gasteiger charges |
| **Docking** | Run all receptor–ligand pairs in parallel via AutoDock Vina |
| **Results** | Sortable results table with affinity, ΔE, ligand efficiency; heatmap; CSV export |
| **Pose Viewer** | Interactive 3D viewer; download receptor / ligand / complex as PDB |
| **Validation** | Symmetric RMSD vs crystal poses; colour-coded classification; PDF export |
| **About & Cite** | Version info, dependency list, formatted citation |

Annotated screenshots of every tab are in [`docs/screenshots/`](docs/screenshots/).

---

## Grid Box Modes

IBDock offers five modes for defining the Vina search space, selectable from the
radio button in the Protein Prep tab:

| Mode | When to use |
|------|-------------|
| **Auto (recommended)** | Default — cascades through co-crystallised ligand → P2Rank → fpocket → blind docking, with a warning at each fallback |
| **Co-crystallised ligand** | Crystal structure with a known drug-like co-crystallised ligand |
| **P2Rank pocket prediction** | Apo structure or when no co-crystallised ligand is present |
| **fpocket pocket prediction** | Geometry-based alternative to P2Rank |
| **Blind docking** | Unknown binding site; covers the full protein extent |

All modes write the final grid parameters to `grid_config.txt` alongside your results
for full reproducibility.

---

## Validation

Re-docking was performed across six structurally diverse protein–ligand complexes
(AutoDock Vina 1.2.3, exhaustiveness = 16, MGLTools 1.5.7, Open Babel 3.1.1).
All six systems passed the 2.0 Å acceptance criterion (Warren et al., 2006).

| PDB | Protein | Family | RMSD (Å) | Result |
|-----|---------|--------|----------|--------|
| 3OCB | AKT1 | Kinase (AGC family) | 0.648 | ✅ Excellent |
| 6C9H | AMPK | Kinase (CAMK family) | 0.948 | ✅ Excellent |
| 6SFO | MAPK14 | Kinase (CMGC family) | 0.646 | ✅ Excellent |
| 1ERR | Oestrogen Receptor α | Nuclear Receptor | 0.811 | ✅ Excellent |
| 4H3X | MMP-9 | Matrix Metalloprotease | 1.408 | ✅ Pass |
| 1HPX | HIV-1 Protease | Aspartic Protease | 1.735 | ✅ Pass |

**Mean RMSD: 1.033 ± 0.445 Å across all six systems. Zero Borderline. Zero Fail.**

Classification: Excellent < 1.0 Å · Pass 1.0–2.0 Å · Borderline 2.0–3.0 Å · Fail ≥ 3.0 Å

---

## Tests

```bash
pip install pytest
pytest Tests/test_ibdock.py -v
```

## Contributing

Bug reports, feature requests, and pull requests are welcome.
Please read [CONTRIBUTING.md](CONTRIBUTING.md) before submitting.

---



Please cite the tools IBDock depends on:

- **AutoDock Vina:** Trott & Olson (2010) *J. Comput. Chem.* 31:455–461 · Eberhardt et al. (2021) *J. Chem. Inf. Model.* 61:3891–3898
- **MGLTools:** Morris et al. (2009) *J. Comput. Chem.* 30:2785–2791
- **Open Babel:** O'Boyle et al. (2011) *J. Cheminform.* 3:33
- **P2Rank:** Krivak & Hoksza (2018) *J. Cheminform.* 10:39
- **fpocket:** Le Guilloux et al. (2009) *BMC Bioinformatics* 10:168

---

## License

MIT — see [LICENSE](LICENSE).
