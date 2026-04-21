# 𝒊-Dock — Batch Molecular Docking with AutoDock Vina

<p align="center">
  <img src="https://img.shields.io/badge/Python-3.8%2B-blue?logo=python" />
  <img src="https://img.shields.io/badge/Streamlit-1.31%2B-FF4B4B?logo=streamlit" />
  <img src="https://img.shields.io/badge/Platform-Windows%20%7C%20macOS%20%7C%20Linux-lightgrey" />
  <img src="https://img.shields.io/badge/License-MIT-green" />
</p>

**𝒊-Dock** is a zero-command-line, browser-based GUI for **batch virtual screening** using [AutoDock Vina](https://vina.scripps.edu/). It covers the full pipeline — protein preparation, ligand preparation, grid box definition, parallel docking, results analysis, and interactive 3D pose viewing — in a single Streamlit application.

---

## Features

| Step | What it does |
|------|-------------|
| **1 · Project Setup** | Create or load a project folder; all files are stored locally |
| **2 · Protein Prep** | Upload PDB files → prepare PDBQT with MGLTools |
| **3 · Ligand Prep** | Upload SDF / MOL2 / PDB → convert & prepare PDBQT with Open Babel + MGLTools |
| **4 · Grid Box** | Manual entry or auto-detect from pocket tools (P2Rank / fpocket) |
| **5 · Docking** | Parallel AutoDock Vina jobs with live progress log |
| **6 · Results** | Ranked table, heatmap, box plots, scatter, per-ligand bar chart, CSV export |
| **7 · Pose Viewer** | Interactive 3D viewer (py3Dmol) — receptor + best ligand pose |
| **8 · About** | Citation references for all underlying tools |

---

## ️ Prerequisites

Install the following **external tools** before running 𝒊-Dock. The app lets you configure their paths in the sidebar.

| Tool | Purpose | Download |
|------|---------|----------|
| **MGLTools 1.5.7** | Protein & ligand preparation | [scripps.edu](https://ccsb.scripps.edu/mgltools/downloads/) |
| **AutoDock Vina 1.2+** | Docking engine | [vina.scripps.edu](https://vina.scripps.edu/downloads/) |
| **Open Babel 3.x** | Ligand format conversion | [openbabel.org](https://openbabel.org/wiki/Get_Open_Babel) |
| **P2Rank** *(optional)* | ML-based pocket prediction | [github.com/rdk/p2rank](https://github.com/rdk/p2rank) |
| **fpocket** *(optional)* | Geometry-based pocket detection | [github.com/Discngine/fpocket](https://github.com/Discngine/fpocket) |

> **Windows users:** P2Rank and fpocket require WSL (Windows Subsystem for Linux).

---

##  Installation

```bash
# 1. Clone the repository
git clone https://github.com/YOUR_USERNAME/i-dock.git
cd i-dock

# 2. Create and activate a virtual environment (recommended)
python -m venv venv
# Windows:
venv\Scripts\activate
# macOS / Linux:
source venv/bin/activate

# 3. Install Python dependencies
pip install -r requirements.txt

# 4. Launch the app
streamlit run i-DOCK.py
```

The app will open at **http://localhost:8501** in your default browser.

---

##  Configuration

On first run, open the **sidebar** and set the paths to your local tool installations:

| Setting | Default (Windows) |
|---------|------------------|
| MGLTools python | `C:\MGLTools-1.5.7\python.exe` |
| prepare_receptor4.py | `C:\MGLTools-1.5.7\AutoDockTools\Utilities24\prepare_receptor4.py` |
| prepare_ligand4.py | `C:\MGLTools-1.5.7\AutoDockTools\Utilities24\prepare_ligand4.py` |
| AutoDock Vina | `C:\vina\vina.exe` |
| Open Babel | `C:\Program Files\OpenBabel-3.1.1\obabel.exe` |

Paths are saved automatically to `<project_dir>/idock_config.json` and reloaded on next launch.

---

## 📁 Project Structure

```
i-dock/
├── i-DOCK.py            ← Main Streamlit application
├── requirements.txt     ← Python dependencies
├── .gitignore
└── README.md

# Created automatically when you run a project:
<your_project>/
├── idock_config.json    ← Saved tool paths
├── prep_receptors/      ← Prepared receptor PDBQT files
├── prep_ligands/        ← Prepared ligand PDBQT files
├── results_vina/        ← Docking output (PDBQT poses + log TXT)
└── ...
```

---

##  Results & Analysis

After docking completes, the **Results** tab shows:

- **Ranked table** — all protein × ligand affinities (kcal/mol), sortable and downloadable as CSV
- **Heatmap** — protein vs ligand affinity matrix
- **Box plot** — affinity distribution per protein
- **Scatter plot** — affinity rank vs score with -7 kcal/mol threshold line
- **Bar chart** — per-ligand comparison across proteins

---

##  Citations

If you use 𝒊-Dock in published research, please cite the underlying tools:

**AutoDock Vina 1.2:**
> Eberhardt J, Santos-Martins D, Tillack AF, Forli S. (2021). AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. *J Chem Inf Model.* 61(8):3891–3898. DOI: 10.1021/acs.jcim.1c00203

**AutoDock Vina (original):**
> Trott O, Olson AJ. (2010). AutoDock Vina: Improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *J Comput Chem.* 31(2):455–461. DOI: 10.1002/jcc.21334

**MGLTools / AutoDockTools:**
> Morris GM et al. (2009). AutoDock4 and AutoDockTools4: Automated docking with selective receptor flexibility. *J Comput Chem.* 30(16):2785–2791. DOI: 10.1002/jcc.21256

**Open Babel:**
> O'Boyle NM et al. (2011). Open Babel: An open chemical toolbox. *J Cheminform.* 3:33. DOI: 10.1186/1758-2946-3-33

**P2Rank:**
> Krivák R, Hoksza D. (2018). P2Rank: machine learning based tool for rapid and accurate prediction of ligand binding sites from protein structure. *J Cheminform.* 10:39. DOI: 10.1186/s13321-018-0285-8

**fpocket:**
> Le Guilloux V, Schmidtke P, Tuffery P. (2009). Fpocket: An open source platform for ligand pocket detection. *BMC Bioinformatics.* 10:168. DOI: 10.1186/1471-2105-10-168

---

## ⚠️ Disclaimer

Docking scores are computational approximations. Results should be interpreted alongside experimental validation and domain expertise. 𝒊-Dock is not intended for clinical decision-making.

---

## 📄 License

This project is released under the **MIT License**. See [LICENSE](LICENSE) for details.
AutoDock Vina, MGLTools, Open Babel, P2Rank, and fpocket are subject to their respective licences.
