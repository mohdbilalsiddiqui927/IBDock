# Scripts
Astex: `astex_naive.py` (naive element-matched Hungarian baseline + box geometry) -> `astex_merge.py` -> `astex_gee.py` (logistic GEE) -> `astex_si_tex.py` (SI tables).
Cross-docking: `cross_final.py` (scoring, bootstrap), `text_stats.py`, `tables_tex.py`, `make_fig_cross.py`.
Repurposing: `rep_score.py` / `rep_score2.py <strategy>` (sPyRMSD of imatinib poses vs crystal), `make_fig_rep.py`.
Figures: `make_fig_box.py`, `make_fig_4_6_7_S1-S4_make_final_figs.py` (see FIGURES_INDEX.md caveat).
`build_*.py` generate the LaTeX from templates; they are provenance only.
Only `make_pub_figures.py` is portable (`python make_pub_figures.py <data_root> <out_dir>`). The other scripts are the authors' working scripts: set the working-directory placeholder at the top (`S`/`WORKDIR`) before running.
Python 3.11+, numpy, pandas, scipy 1.18.1, statsmodels 0.15.0, rdkit, spyrmsd, matplotlib.
