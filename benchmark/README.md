# Astex Diverse Set Benchmark

This folder reproduces the statistical analysis reported in the manuscript
for the 85-complex Astex Diverse Set redocking benchmark (Tables 2-4 and
the associated in-text statistics).

Docking itself was run directly through the IBDock Provenire application
(batch execution across all 85 complexes x 4 search-space strategies); no
separate docking script is included here. This folder covers the
**statistical analysis of the resulting per-complex dataset** only.

## Contents

```
benchmark/
├── README.md                      # this file
├── astex_diverse_set_85x4.csv     # the 340-row per-complex x strategy dataset
└── scripts/
    ├── statistical_analysis.py    # reproduces Tables 2, 3, 4 and in-text stats
    └── results/                   # output of statistical_analysis.py (generated)
```

## Dataset schema

`astex_diverse_set_85x4.csv` has one row per complex x strategy
(85 complexes x 4 strategies = 340 rows):

| Column | Type | Description |
|---|---|---|
| `PDB_ID` | str | e.g. `1G9V` |
| `Protocol` | str | `Co-crystallized`, `P2Rank`, `fpocket`, or `Blind` |
| `Score` | float | Vina mode-1 affinity (kcal/mol) |
| `sPy_Top1` | float | sPyRMSD Top-1 pose RMSD (Å); blank if unresolved |
| `sPy_Best` | float | sPyRMSD best-pose RMSD (Å); blank if unresolved |
| `Bridge_Top1` | float | coordinate-bridge Top-1 RMSD (Å); blank if unresolved |
| `Bridge_Best` | float | coordinate-bridge best-pose RMSD (Å); blank if unresolved |
| `sPy_status` | str | `Numerical` or `Unresolved` |
| `Bridge_status` | str | `Numerical` or `Unresolved` |

Run `python scripts/statistical_analysis.py --write-template` to generate
an empty CSV with exactly these headers.

## Reproducing the results

```bash
pip install pandas numpy scipy statsmodels
cd benchmark
python scripts/statistical_analysis.py astex_diverse_set_85x4.csv --outdir scripts/results
```

This prints, and writes as CSV to `scripts/results/`:

- **Table 2** — Top-1 / best-pose success per strategy, resolved-case and
  all-85 views, with Wilson 95% confidence intervals
- **Cochran's Q** (binary endpoint, df = 3) for Top-1 and best-pose
- **Table 3** — Holm-adjusted pairwise comparisons: exact two-sided McNemar
  (binary, n = 85) and paired Wilcoxon signed-rank (continuous RMSD,
  complexes resolved under all four strategies)
- **Friedman's Q** on the continuous RMSD subset
- **Table 4** — coordinate-bridge applicability and success rates per
  strategy
- **Bridge vs. sPyRMSD agreement** — Spearman ρ, mean/median absolute
  difference, over all strategy-complex pairs where both methods returned
  a numerical value

Success criterion throughout: RMSD < 2.0 Å. For the primary (all-85)
endpoint, unresolved sPyRMSD correspondence is coded as non-success,
matching the manuscript's conservative primary endpoint; resolved-case
rates are reported alongside as a sensitivity view, exactly as described
in the manuscript's Statistical Analysis section.
