# Changelog

All notable changes to IBDock are documented here.  
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [1.0.0] — Initial Release

### Added
- Full AutoDock Vina batch docking pipeline as a Streamlit GUI
- Project-based file management with per-project config persistence
- **Tab 1 — Project Setup:** create/load project folder
- **Tab 2 — Protein Prep:** upload PDB → PDBQT via MGLTools
- **Tab 3 — Ligand Prep:** upload SDF/MOL2/PDB → PDBQT via Open Babel + MGLTools
- **Tab 4 — Grid Box:** manual entry + P2Rank / fpocket auto-detection (WSL)
- **Tab 5 — Docking:** parallel Vina jobs with live progress, log streaming, elapsed/ETA timer
- **Tab 6 — Results:** ranked table, heatmap, box plot, scatter plot, bar chart, CSV export
- **Tab 7 — Pose Viewer:** interactive 3D viewer (py3Dmol) for receptor + ligand poses
- **Tab 8 — About:** tool credits and citation references
- Cross-platform defaults for Windows, macOS, and Linux
- Config auto-save to `IBdock_config.json` in project directory
