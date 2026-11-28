# Benchmark data and analysis

Everything behind the manuscript's numbers. Conventions: strategies are `cocryst` (co-crystallized ligand), `fpocket`, `p2rank`, `blind`; distances in angstrom; scores in kcal/mol; AutoDock Vina 1.2.7; success = sPyRMSD < 2.0 A against the crystal pose; unresolved correspondence counts as non-success in the primary endpoint.

## Design (what each dataset is)
| Dataset | Design | Settings |
|---|---|---|
| Astex Diverse Set | 85 complexes x 4 strategies | exhaustiveness 8, 16, 32 x seeds 1, 42, 123 (3,060 dockings); exhaustiveness 16 is the reference slice in the main tables |
| Cross-docking | ESR1, MAPK1, PPARG (LIT-PCBA templates), 1,320 receptor-ligand pairs x 4 strategies | exhaustiveness 16, seed 42 only |
| Repurposing | imatinib into ABL1 (1IEP), KIT (1T46), SRC (2SRC), EGFR (1M17) x 4 strategies | exhaustiveness 16, seed 42 only |
| Aging case study | 7 aging-associated complexes x 4 strategies (28 conditions); plus 7 ligands x 7 receptors with P2Rank boxes (49 dockings) | exhaustiveness 16, seed 42 only |

## data/astex/
- `raw_runs/astex_ex{8|16|32}_seed{1|42|123}_{strategy}_docking_results.csv` (36 files): raw IBDock output, one row per complex. Columns: Protein, Ligand, Best Affinity, Mode 2, delta E Mode1->2, Ligand Efficiency, Heavy Atoms, Modes Saved, Log File. Docked PDBQT poses are in the Zenodo archive, not here.
- `analysis/astex_merged.csv`: one row per complex x condition (3,060 rows). `ex`, `seed`, `pdb`, `strategy`, `score` (Top-1 Vina score); `sp_top1`/`sp_best` sPyRMSD of the Top-1/best saved pose (empty = unresolved); `br_top1`/`br_best` coordinate-bridge RMSD; `raw4` raw-coordinate value as recorded in the Table S3 files; `n_models` saved poses; `id_*` identity-mapping baseline (rejected, kept for the record); `nn_top1`/`nn_best` naive element-matched optimal-assignment RMSD; `_merge` join status.
- `analysis/astex_naive_scored.csv`: the naive-assignment scores alone.
- `analysis/astex_geometry.csv`: search box per complex x condition: `vol` (A^3), `sx, sy, sz` edge lengths, `offset` (box centre to reference-ligand centroid), `inside` (% of reference heavy atoms inside the box), `full` (all inside). Boxes are identical across seeds and exhaustiveness levels.
- `analysis/astex_extra.json`: results behind SI Tables S13-S15 (all three exhaustiveness levels).
- `analysis/Table_S3a-i_*.docx`: per-run sPyRMSD and bridge tables from which `astex_merged.csv` was built.

## data/cross_docking/
- `inputs/`: pair manifest (`manifest.csv`), ligand preparation record (`ligand_qc.csv`), receptor, ligand, and reference-pose files, and the batch lists. `inputs/README.md` describes the design as planned; the final analysis uses exhaustiveness 16, seed 42 only.
- `scored/`: one row per pair x strategy in `ALL3_primary_scored.csv`, `ESR1_all_strategies_scored.csv`, `MAPK1_PPARG_all_scored.csv` (`top1_rmsd`, `best_rmsd`, `t1`/`bp` success flags, `site` site capture, `cognate`, `tier`, `pocket_CA_rmsd`, `centroid_mismatch_A`); summaries `cross_main.csv` (rates with two-way cluster-bootstrap CIs), `cross_tier.csv`, `cross_cond.csv` (success given/not given site capture), `cross_paired.csv` (paired differences), `cross_flagged.csv`, `ref_clash.csv` (minimum receptor-reference distance, clash flag), `ESR1_summary_table.csv`.

## data/repurposing/
`inputs/` (receptor PDBs with and without ligand, imatinib structure, crystal reference poses; `inputs/README.txt`), `raw_docking_outputs/<strategy>/` (prepared files, grids, Vina logs, docked poses), `activity_table_VERIFIED.csv` (imatinib activities from ChEMBL, compound CHEMBL941). Scored values are in SI Table S11.

## data/aging_case_study/
See `data/aging_case_study/README.md`. Vina logs (`*/results/`), poses (`*/docking/`), boxes (`*/grid/`), prepared and raw inputs. `p2rank/results/docking_results.csv` holds the 49 cross-protein scores (SI Tables S5b, S5c); the 28 self-docking conditions are in SI Table S5.

## scripts/
Run order: `astex_naive.py` -> `astex_merge.py` -> `astex_gee.py` -> `astex_extra.py` -> `astex_extra_tex.py` -> `astex_si_tex.py` (Astex); `cross_final.py` -> `text_stats.py` -> `tables_tex.py` -> `make_fig_cross.py` (cross-docking); `rep_score.py`/`rep_score2.py <strategy>` -> `make_fig_rep.py` (repurposing); `make_fig_box.py` (Figure 5). `build_*.py` assemble the LaTeX from templates (provenance only). `make_fig_4_6_7_S1-S4_make_final_figs.py` is the plotting code for Figures 4, 6, 7 and S1-S4; it needs intermediate files (`s3_new.pkl`, `out/*.csv`) built from `astex_merged.csv` and is provided as the record of the plotting code.
Install: `pip install -r requirements-analysis.txt`. `make_pub_figures.py` takes the data and output folders as arguments and runs as is. The other scripts contain a working-directory placeholder (`WORKDIR`, `S`) at the top that points to the authors' working layout; edit it before running.

## Known limitations of this package
- Docked PDBQT poses for the Astex runs are not included (size); they are in the Zenodo archive.
- The cross-docking P2Rank condition used the fixed 22 A fallback box (see manuscript Methods).
- SI Table S11 values for repurposing are not stored as a separate scored CSV.
