# Figures index
| No. | Content | Files | Script | Data used |
|---|---|---|---|---|
| Figure 1 | IBDock workflow architecture | figures/Main/Figure_1 | artwork (no script) | n/a |
| Figure 2 | Search-space selection and grid provenance | Figure_2 | artwork | n/a |
| Figure 3 | Reference-pose validation architecture | Figure_3 | artwork | n/a |
| Figure 4 | Astex pose recovery across the four strategies (exh. 16, 3 seeds) | Figure_4 | benchmark/scripts/make_fig_4_6_7_S1-S4_make_final_figs.py | benchmark/data/Astex_analysis/astex_merged.csv |
| Figure 5 | Box geometry and Top-1 success (enclosure) | Figure_5 | benchmark/scripts/make_fig_box.py | benchmark/data/Astex_analysis/astex_geometry.csv, astex_merged.csv |
| Figure 6 | Sensitivity to exhaustiveness and seed | Figure_6 | benchmark/scripts/make_fig_4_6_7_S1-S4_make_final_figs.py | astex_merged.csv (all 36 conditions) |
| Figure 7 | Coordinate-bridge vs sPyRMSD agreement | Figure_7 | benchmark/scripts/make_fig_4_6_7_S1-S4_make_final_figs.py | astex_merged.csv (sp_*, br_* columns) |
| Figure 8 | Cross-docking by target and strategy | Figure_8 | benchmark/scripts/make_fig_cross.py (after cross_final.py) | benchmark/data/CrossDocking/*.csv |
| Figure 9 | Imatinib repurposing, four kinases x four strategies | Figure_9 | benchmark/scripts/make_fig_rep.py | benchmark/data/Repurposing/ |
| Figure S1 | ECDF of Top-1 and best-pose sPyRMSD | figures/SI/Figure_S1 | make_fig_4_6_7_S1-S4_make_final_figs.py | astex_merged.csv |
| Figure S2 | Outcome breakdown per strategy | Figure_S2 | same | astex_merged.csv |
| Figure S3 | Correspondence-route coverage | Figure_S3 | same | astex_merged.csv |
| Figure S4 | Top-1 score distribution and recovery | Figure_S4 | same | astex_merged.csv |

Caveat on reproducibility: `make_fig_4_6_7_S1-S4_make_final_figs.py` was written against intermediate files (`s3_new.pkl`, `out/agg.csv`, `out/rates.csv`, `out/cubes.npz`) built from `benchmark/data/Astex_analysis/astex_merged.csv` and the Table S3 .docx files. Those intermediates are not in this package, so that script is provided as the record of the plotting code and styling, and needs those inputs to be rebuilt before it will run. Scripts for Figures 5, 8 and 9 read the CSVs here after editing the path constants at the top.
