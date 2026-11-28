# Figures index
Every panel is a separate file (SVG with live text, PDF, 600 dpi PNG), named as in the manuscript (panels a, b, c; composites are assembled in the manuscript source). File names have no spaces. Figures 1-3 and 8 are artwork (PNG + SVG).
Style: Palatino (TeX Gyre Pagella in the files; the SVG asks for Palatino Linotype first); colours Co-crystal #2b5aa0, fpocket #18512f, P2Rank #550022, Blind #554400; no grid lines.
Main figures are in `figures/Main/`, SI figures in `figures/SI/`.

| Figure | Files | Content |
|---|---|---|
| 1 | Figure1_Workflow | IBDock workflow architecture (artwork) |
| 2 | Figure2_final | Search-space selection and grid provenance (artwork) |
| 3 | Figure3_final | Reference-pose validation architecture (artwork) |
| 4 | Figure_4a, 4b | Astex success by strategy; per-complex sPyRMSD |
| 5 | Figure_5a, 5b, 5c | Box volume; success by ligand enclosure; GEE odds ratios |
| 6 | Figure_6a, 6b, 6c | Exhaustiveness/seed sensitivity; seed consistency |
| 7 | Figure_7a, 7b | Bridge vs sPyRMSD (Top-1, best) |
| 8 | Figure_8 | Aging case study (6SFO pose, 2D map, cross-protein matrix, ranking); data in `benchmark/data/aging_case_study/` (artwork) |
| 9 | Figure_9a, 9b, 9c | Cross-docking Top-1; site capture; Top-1 given site |
| 10 | Figure_10a, 10b | Imatinib RMSD; Vina score |
| S1 | Figure_S1a, S1b | ECDF of sPyRMSD |
| S2 | Figure_S2 | Outcome breakdown per strategy |
| S3 | Figure_S3 | Correspondence-route coverage |
| S4 | Figure_S4a, S4b | Score distribution; score vs recovery (ROC-AUC) |

Reproduce Figures 4-7, 9, 10 and S1-S4:
`python benchmark/scripts/make_pub_figures.py benchmark/data figures`
(needs numpy, pandas, scipy, statsmodels, matplotlib; reads `astex_merged.csv`, `astex_geometry.csv`, `gee_res.pkl`, `cross_main.csv`, `cross_cond.csv`; the Figure 10 values are embedded and equal SI Table S11). Figure 8 panel values are SI Tables S5, S5b and S5c.
