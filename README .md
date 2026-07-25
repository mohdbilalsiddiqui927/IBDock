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
- **Multi-chain-safe grid centring** — when a structure contains multiple crystallographic copies of the same complex (common in the PDB), ligand coordinates are restricted to a single representative chain rather than pooled across copies, avoiding a grid centre that corresponds to no real binding site
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
git clone https://github.com/Mohdbilalsiddiqui927/ibdock.git
cd ibdock
pip install -r requirements.txt
streamlit run IBDock.py
```

This opens IBDock in your browser at `http://localhost:8501`. On first launch, open
the **⚙ Settings** panel in the sidebar and point IBDock at your MGLTools, Vina, and
Open Babel executables. Paths are saved to `IBDock_config.json` and remembered across
sessions.

**Note (Windows):** if the `python` command does not work in your terminal, use the
Python launcher instead: `py -3 -m pip install -r requirements.txt` and `py -3 -m streamlit run IBDock.py`.

---

## Quick Start

The `example/` folder contains a prepared 1HNN system (human PNMT, bound to the
inhibitor SK&F 29661) ready to dock:

```
example/
├── receptor/   1HNN.pdb          raw PDB file as downloaded from RCSB
├── ligand/     1HNN_ligand.pdb   co-crystallised ligand
└── reference/  1HNN_ligand.pdb   crystal pose for RMSD validation
```

1. Upload `1HNN.pdb` in **Protein Prep** — the grid box is detected automatically
2. Upload `1HNN_ligand.pdb` in **Ligand Prep**
3. Click **Run Docking** in the Docking tab
4. Upload `1HNN_ligand.pdb` (reference copy) in **Validation** and compute RMSD

**Expected result:** RMSD ≈ 1.02 Å (Pass)

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

For structures containing multiple crystallographic copies of the same complex
(e.g. several chains each with their own bound ligand), grid derivation is
automatically restricted to a single representative chain rather than pooling
ligand coordinates across copies — see Features above.

All modes write the final grid parameters to `grid_config.txt` alongside your results
for full reproducibility.

---

## Validation

Re-docking was performed across the full 85-complex Astex Diverse Set (Hartshorn et
al., 2007), using AutoDock Vina 1.2.3, exhaustiveness = 16, MGLTools 1.5.7, Open Babel
3.1.1. Top-scored (Top-1) pose pass rate: **42/85 (49.4%)** at the standard RMSD ≤ 2.0 Å
criterion, consistent with published AutoDock Vina performance on the same benchmark
(Buttenschoen et al., 2024).

A worked example is provided in the `example/` folder using 1HNN (human
phenylethanolamine N-methyltransferase, PNMT), co-crystallised with its inhibitor
1,2,3,4-tetrahydroisoquinoline-7-sulfonamide (SK&F 29661). Re-docking this system with
IBDock using the co-crystallised ligand for automated grid box derivation yielded an
RMSD of 1.020 Å against the crystal pose (Pass, < 2.0 Å threshold).

---

## Tests

```bash
pip install pytest
pytest Tests/test_ibdock.py -v
```

---

## Contributing

Bug reports, feature requests, and pull requests are welcome.
Please read [CONTRIBUTING.md](CONTRIBUTING.md) before submitting.

---

## Citation

If IBDock is useful in your work, please cite:

```bibtex
@article{ibdock2026,
  author = {Siddiqui, Bilal and Tiwari, Sakshi and Khan, Imran A and Abdin, M Z},
  title  = {{IBDock}: Browser-Based Batch Docking and Automated Validation
            with {AutoDock Vina}}
}
```

*(Full journal details will be added upon publication.)*

Please also cite the tools IBDock depends on:

- **AutoDock Vina:** Trott & Olson (2010) *J. Comput. Chem.* 31:455–461 · Eberhardt et al. (2021) *J. Chem. Inf. Model.* 61:3891–3898
- **MGLTools:** Morris et al. (2009) *J. Comput. Chem.* 30:2785–2791
- **Open Babel:** O'Boyle et al. (2011) *J. Cheminform.* 3:33
- **P2Rank:** Krivak & Hoksza (2018) *J. Cheminform.* 10:39
- **fpocket:** Le Guilloux et al. (2009) *BMC Bioinformatics* 10:168

---

## License

MIT — see [LICENSE](LICENSE).

