#!/usr/bin/env python3
"""
statistical_analysis.py
========================
Reproduces the statistical analysis reported in:

    Siddiqui, B.; Tiwari, S.; Khan, I. A.; Abdin, M. Z. "IBDock Provenire:
    Search-Space Capture and Pose Validation for Molecular Docking."
    J. Chem. Inf. Model. (2026), in review.

Given the per-complex, per-strategy Astex Diverse Set benchmark dataset
(same schema as the manuscript's Supporting Information Table S2), this
script reproduces:

    * Table 2  - Astex Redocking Results by Search-Space Strategy
                 (Top-1 / best-pose success, resolved-case and all-85 views)
    * Table 3  - Holm-Adjusted Paired Comparisons across Search-Space
                 Strategies (exact two-sided McNemar on the binary endpoint;
                 paired Wilcoxon signed-rank on the continuous RMSD endpoint)
    * Table 4  - Coordinate-Bridge Performance, Applicability, and Agreement
                 with sPyRMSD (Spearman rho, mean/median absolute difference)
    * In-text statistics - Cochran's Q (binary, df=3) and Friedman's Q
                 (continuous RMSD, resolved-under-all-four subset)

Expected input CSV
-------------------
One row per complex x strategy (340 rows for the full Astex Diverse Set):

    PDB_ID          str    e.g. "1G9V"
    Protocol        str    one of: "Co-crystallized", "P2Rank", "fpocket", "Blind"
                           (case-insensitive; "F-pocket" is also accepted)
    Score           float  Vina mode-1 affinity (kcal/mol)          [not used
                           in the statistics below, carried through only]
    sPy_Top1        float  sPyRMSD Top-1 pose RMSD (Angstrom), blank/NaN if
                           sPyRMSD could not establish a correspondence
    sPy_Best        float  sPyRMSD best-pose RMSD (Angstrom), blank/NaN if
                           unresolved
    Bridge_Top1     float  coordinate-bridge Top-1 RMSD (Angstrom), blank/NaN
                           if the bridge could not establish a correspondence
    Bridge_Best     float  coordinate-bridge best-pose RMSD (Angstrom),
                           blank/NaN if unresolved
    sPy_status      str    "Numerical" or "Unresolved"
    Bridge_status   str    "Numerical" or "Unresolved"

A template header with these exact column names is written to
``template_astex_dataset.csv`` if you run this script with --write-template.

Success criterion (a priori, per the manuscript): RMSD < 2.0 Angstrom.
For the primary (all-85) endpoint, unresolved sPyRMSD correspondence is
coded as non-success (S = 0); resolved-case rates are reported alongside
as a sensitivity view.

Usage
-----
    python statistical_analysis.py astex_diverse_set_85x4.csv --outdir results/

Dependencies: pandas, numpy, scipy, statsmodels
    pip install pandas numpy scipy statsmodels
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import friedmanchisquare, spearmanr, wilcoxon
from statsmodels.stats.contingency_tables import mcnemar
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportion_confint

# ----------------------------------------------------------------------- #
# Configuration
# ----------------------------------------------------------------------- #

SUCCESS_THRESHOLD_A = 2.0  # Angstrom

# Canonical strategy order used throughout (matches manuscript Tables 2-4)
STRATEGY_ORDER = ["Co-crystallized", "fpocket", "P2Rank", "Blind"]

# Accepts common spelling variants found across the SI dataset / user exports
_STRATEGY_ALIASES = {
    "co-crystallized": "Co-crystallized",
    "co-crystallized ligand": "Co-crystallized",
    "cocrystallized": "Co-crystallized",
    "co-crystal": "Co-crystallized",
    "fpocket": "fpocket",
    "f-pocket": "fpocket",
    "p2rank": "P2Rank",
    "blind": "Blind",
    "blind docking": "Blind",
}

# The six pairwise comparisons reported in Table 3, in the manuscript's order
PAIRWISE_COMPARISONS = [
    ("Co-crystallized", "fpocket"),
    ("Co-crystallized", "P2Rank"),
    ("Co-crystallized", "Blind"),
    ("fpocket", "P2Rank"),
    ("fpocket", "Blind"),
    ("P2Rank", "Blind"),
]


# ----------------------------------------------------------------------- #
# Data loading
# ----------------------------------------------------------------------- #

def write_template(path: Path) -> None:
    cols = [
        "PDB_ID", "Protocol", "Score", "sPy_Top1", "sPy_Best",
        "Bridge_Top1", "Bridge_Best", "sPy_status", "Bridge_status",
    ]
    pd.DataFrame(columns=cols).to_csv(path, index=False)
    print(f"Wrote empty template with correct headers to: {path}")


def load_dataset(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)

    required = {
        "PDB_ID", "Protocol", "sPy_Top1", "sPy_Best",
        "Bridge_Top1", "Bridge_Best", "sPy_status", "Bridge_status",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Input CSV is missing required column(s): {sorted(missing)}\n"
            "Run with --write-template to see the expected schema."
        )

    df["Protocol"] = (
        df["Protocol"].astype(str).str.strip().str.lower().map(_STRATEGY_ALIASES)
    )
    bad = df["Protocol"].isna()
    if bad.any():
        raise ValueError(
            "Unrecognized Protocol value(s) in rows: "
            f"{df.index[bad].tolist()}. Expected one of "
            f"{sorted(set(_STRATEGY_ALIASES.values()))}."
        )

    for col in ("sPy_Top1", "sPy_Best", "Bridge_Top1", "Bridge_Best"):
        df[col] = pd.to_numeric(df[col], errors="coerce")

    for col in ("sPy_status", "Bridge_status"):
        df[col] = df[col].astype(str).str.strip().str.capitalize()
        bad = ~df[col].isin(["Numerical", "Unresolved"])
        if bad.any():
            raise ValueError(
                f"Column '{col}' must be 'Numerical' or 'Unresolved'; "
                f"found unexpected value(s) in rows {df.index[bad].tolist()}."
            )

    n_complexes = df["PDB_ID"].nunique()
    n_rows = len(df)
    print(f"Loaded {n_rows} rows covering {n_complexes} unique complexes.")
    if n_rows != n_complexes * len(STRATEGY_ORDER):
        print(
            f"WARNING: expected {n_complexes * len(STRATEGY_ORDER)} rows "
            f"({n_complexes} complexes x {len(STRATEGY_ORDER)} strategies) "
            f"but found {n_rows}. Check for missing/duplicate rows before "
            "trusting the results below."
        )
    return df


# ----------------------------------------------------------------------- #
# Binary success matrices (complexes x strategies), unresolved -> 0
# ----------------------------------------------------------------------- #

def binary_matrix(df: pd.DataFrame, rmsd_col: str, status_col: str) -> pd.DataFrame:
    """Wide matrix: rows = PDB_ID, columns = strategy, values = 0/1 success.

    Unresolved correspondence is coded as non-success (0), matching the
    manuscript's conservative all-85 primary endpoint.
    """
    success = (df[rmsd_col] < SUCCESS_THRESHOLD_A) & (df[status_col] == "Numerical")
    wide = (
        df.assign(success=success.astype(int))
        .pivot(index="PDB_ID", columns="Protocol", values="success")
        .reindex(columns=STRATEGY_ORDER)
    )
    return wide


# ----------------------------------------------------------------------- #
# Table 2: per-strategy resolved/all-85 success rates with Wilson CI
# ----------------------------------------------------------------------- #

def table2_strategy_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    n_all = df["PDB_ID"].nunique()
    for strategy in STRATEGY_ORDER:
        sub = df[df["Protocol"] == strategy]
        resolved = (sub["sPy_status"] == "Numerical").sum()

        top1_success_all = int(
            ((sub["sPy_Top1"] < SUCCESS_THRESHOLD_A) & (sub["sPy_status"] == "Numerical")).sum()
        )
        best_success_all = int(
            ((sub["sPy_Best"] < SUCCESS_THRESHOLD_A) & (sub["sPy_status"] == "Numerical")).sum()
        )

        resolved_sub = sub[sub["sPy_status"] == "Numerical"]
        top1_success_resolved = int((resolved_sub["sPy_Top1"] < SUCCESS_THRESHOLD_A).sum())
        best_success_resolved = int((resolved_sub["sPy_Best"] < SUCCESS_THRESHOLD_A).sum())

        ci_top1_lo, ci_top1_hi = proportion_confint(
            top1_success_all, n_all, alpha=0.05, method="wilson"
        )
        ci_best_lo, ci_best_hi = proportion_confint(
            best_success_all, n_all, alpha=0.05, method="wilson"
        )

        rows.append({
            "Search space": strategy,
            "sPyRMSD resolved": f"{resolved}/{n_all}",
            "Top-1 resolved cases": f"{top1_success_resolved}/{resolved} "
                                     f"({100 * top1_success_resolved / resolved:.1f}%)" if resolved else "n/a",
            "Top-1 all 85": f"{top1_success_all}/{n_all} ({100 * top1_success_all / n_all:.1f}%)",
            "Top-1 95% CI": f"{100 * ci_top1_lo:.1f}-{100 * ci_top1_hi:.1f}%",
            "Best-pose resolved cases": f"{best_success_resolved}/{resolved} "
                                         f"({100 * best_success_resolved / resolved:.1f}%)" if resolved else "n/a",
            "Best-pose all 85": f"{best_success_all}/{n_all} ({100 * best_success_all / n_all:.1f}%)",
            "Best-pose 95% CI": f"{100 * ci_best_lo:.1f}-{100 * ci_best_hi:.1f}%",
        })
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------- #
# Cochran's Q (binary endpoint, all 85, df = k-1 = 3)
# ----------------------------------------------------------------------- #

def cochrans_q(binary_wide: pd.DataFrame) -> tuple[float, int, float]:
    """Manual Cochran's Q implementation (no extra dependency needed).

    x: (n_subjects x k_conditions) 0/1 matrix.
    Returns (Q, df, p_value).
    """
    from scipy.stats import chi2

    x = binary_wide.to_numpy(dtype=float)
    n, k = x.shape
    col_sums = x.sum(axis=0)          # successes per strategy
    row_sums = x.sum(axis=1)          # successes per complex

    numerator = k * (k - 1) * np.sum((col_sums - col_sums.mean()) ** 2)
    denominator = k * row_sums.sum() - np.sum(row_sums ** 2)
    if denominator == 0:
        return float("nan"), k - 1, float("nan")

    Q = numerator / denominator
    df = k - 1
    p = chi2.sf(Q, df)
    return Q, df, p


# ----------------------------------------------------------------------- #
# Table 3: Holm-adjusted pairwise McNemar (binary) and Wilcoxon (RMSD)
# ----------------------------------------------------------------------- #

def pairwise_mcnemar(binary_wide: pd.DataFrame) -> pd.DataFrame:
    raw_p = []
    labels = []
    for a, b in PAIRWISE_COMPARISONS:
        xa = binary_wide[a].to_numpy()
        xb = binary_wide[b].to_numpy()
        # 2x2 contingency table of (a, b) paired outcomes
        table = np.array([
            [((xa == 1) & (xb == 1)).sum(), ((xa == 1) & (xb == 0)).sum()],
            [((xa == 0) & (xb == 1)).sum(), ((xa == 0) & (xb == 0)).sum()],
        ])
        result = mcnemar(table, exact=True, correction=False)
        raw_p.append(result.pvalue)
        labels.append(f"{a} vs {b}")

    _, holm_p, _, _ = multipletests(raw_p, method="holm")
    return pd.DataFrame({"Comparison": labels, "raw_p": raw_p, "holm_p": holm_p})


def pairwise_wilcoxon(resolved_wide: pd.DataFrame) -> pd.DataFrame:
    """resolved_wide: complexes resolved under ALL four strategies only."""
    raw_p = []
    labels = []
    for a, b in PAIRWISE_COMPARISONS:
        stat = wilcoxon(resolved_wide[a], resolved_wide[b])
        raw_p.append(stat.pvalue)
        labels.append(f"{a} vs {b}")

    _, holm_p, _, _ = multipletests(raw_p, method="holm")
    return pd.DataFrame({"Comparison": labels, "raw_p": raw_p, "holm_p": holm_p})


# ----------------------------------------------------------------------- #
# Friedman's test on continuous RMSD (complexes resolved under all 4)
# ----------------------------------------------------------------------- #

def resolved_under_all_four(df: pd.DataFrame, rmsd_col: str, status_col: str) -> pd.DataFrame:
    wide_rmsd = df.pivot(index="PDB_ID", columns="Protocol", values=rmsd_col).reindex(
        columns=STRATEGY_ORDER
    )
    wide_status = df.pivot(index="PDB_ID", columns="Protocol", values=status_col).reindex(
        columns=STRATEGY_ORDER
    )
    mask = (wide_status == "Numerical").all(axis=1)
    return wide_rmsd[mask]


def friedman_test(resolved_wide: pd.DataFrame) -> tuple[float, float]:
    stat = friedmanchisquare(*[resolved_wide[s] for s in STRATEGY_ORDER])
    return stat.statistic, stat.pvalue


# ----------------------------------------------------------------------- #
# Table 4: coordinate-bridge performance and agreement with sPyRMSD
# ----------------------------------------------------------------------- #

def table4_bridge_summary(df: pd.DataFrame) -> pd.DataFrame:
    n_all = df["PDB_ID"].nunique()
    rows = []
    for strategy in STRATEGY_ORDER:
        sub = df[df["Protocol"] == strategy]
        bridge_numerical = (sub["Bridge_status"] == "Numerical").sum()
        bridge_unresolved = n_all - bridge_numerical

        top1_all = int(
            ((sub["Bridge_Top1"] < SUCCESS_THRESHOLD_A) & (sub["Bridge_status"] == "Numerical")).sum()
        )
        best_all = int(
            ((sub["Bridge_Best"] < SUCCESS_THRESHOLD_A) & (sub["Bridge_status"] == "Numerical")).sum()
        )
        resolved_sub = sub[sub["Bridge_status"] == "Numerical"]
        top1_resolved = int((resolved_sub["Bridge_Top1"] < SUCCESS_THRESHOLD_A).sum())
        best_resolved = int((resolved_sub["Bridge_Best"] < SUCCESS_THRESHOLD_A).sum())

        both_numerical = (
            (sub["sPy_status"] == "Numerical") & (sub["Bridge_status"] == "Numerical")
        ).sum()

        rows.append({
            "Metric / Strategy": strategy,
            "Bridge numerical": f"{bridge_numerical}/{n_all}",
            "Bridge unresolved": f"{bridge_unresolved}/{n_all}",
            "Top-1 success, all 85": f"{top1_all}/{n_all} ({100 * top1_all / n_all:.1f}%)",
            "Best-pose success, all 85": f"{best_all}/{n_all} ({100 * best_all / n_all:.1f}%)",
            "Top-1 success, bridge-resolved": f"{top1_resolved}/{bridge_numerical} "
                                               f"({100 * top1_resolved / bridge_numerical:.1f}%)" if bridge_numerical else "n/a",
            "Best-pose success, bridge-resolved": f"{best_resolved}/{bridge_numerical} "
                                                   f"({100 * best_resolved / bridge_numerical:.1f}%)" if bridge_numerical else "n/a",
            "Both methods numerical": int(both_numerical),
        })
    return pd.DataFrame(rows)


def bridge_agreement(df: pd.DataFrame, which: str) -> dict:
    """which: 'Top1' or 'Best'. Returns Spearman rho, mean/median abs diff, n."""
    spy_col, bridge_col = f"sPy_{which}", f"Bridge_{which}"
    both = df[(df["sPy_status"] == "Numerical") & (df["Bridge_status"] == "Numerical")]
    n = len(both)
    if n == 0:
        return {"n": 0, "rho": float("nan"), "mean_abs_diff": float("nan"), "median_abs_diff": float("nan")}

    rho, p = spearmanr(both[spy_col], both[bridge_col])
    abs_diff = (both[spy_col] - both[bridge_col]).abs()
    return {
        "n": n,
        "rho": rho,
        "p_value": p,
        "mean_abs_diff": abs_diff.mean(),
        "median_abs_diff": abs_diff.median(),
    }


# ----------------------------------------------------------------------- #
# Main
# ----------------------------------------------------------------------- #

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv_path", nargs="?", help="Path to the 85x4 Astex benchmark dataset CSV")
    parser.add_argument("--outdir", default="results", help="Directory to write output CSVs/report (default: results/)")
    parser.add_argument("--write-template", action="store_true", help="Write an empty template CSV with the expected columns and exit")
    args = parser.parse_args()

    if args.write_template:
        write_template(Path("template_astex_dataset.csv"))
        return

    if not args.csv_path:
        parser.error("csv_path is required unless --write-template is given")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_dataset(Path(args.csv_path))

    # --- Table 2 -------------------------------------------------------
    t2 = table2_strategy_summary(df)
    t2.to_csv(outdir / "table2_strategy_summary.csv", index=False)
    print("\n=== Table 2: Astex Redocking Results by Search-Space Strategy ===")
    print(t2.to_string(index=False))

    # --- Cochran's Q (binary, all 85) -----------------------------------
    top1_wide = binary_matrix(df, "sPy_Top1", "sPy_status")
    best_wide = binary_matrix(df, "sPy_Best", "sPy_status")

    q_top1, df_top1, p_top1 = cochrans_q(top1_wide)
    q_best, df_best, p_best = cochrans_q(best_wide)
    print(f"\nCochran's Q (Top-1):     Q = {q_top1:.3f}, df = {df_top1}, P = {p_top1:.3e}")
    print(f"Cochran's Q (Best-pose): Q = {q_best:.3f}, df = {df_best}, P = {p_best:.3e}")

    # --- Table 3a: pairwise McNemar + Holm -------------------------------
    mcnemar_top1 = pairwise_mcnemar(top1_wide)
    mcnemar_best = pairwise_mcnemar(best_wide)

    # --- Continuous RMSD subset (resolved under all four strategies) ----
    resolved_top1 = resolved_under_all_four(df, "sPy_Top1", "sPy_status")
    resolved_best = resolved_under_all_four(df, "sPy_Best", "sPy_status")
    print(f"\nComplexes resolved under all four strategies: n = {len(resolved_top1)}")

    friedman_q_top1, friedman_p_top1 = friedman_test(resolved_top1)
    friedman_q_best, friedman_p_best = friedman_test(resolved_best)
    print(f"Friedman's Q (Top-1 RMSD):     Q = {friedman_q_top1:.3f}, P = {friedman_p_top1:.3e}")
    print(f"Friedman's Q (Best-pose RMSD): Q = {friedman_q_best:.3f}, P = {friedman_p_best:.3e}")

    wilcoxon_top1 = pairwise_wilcoxon(resolved_top1)
    wilcoxon_best = pairwise_wilcoxon(resolved_best)

    t3 = pd.DataFrame({
        "Comparison": mcnemar_top1["Comparison"],
        "McNemar Top-1 (Holm p)": mcnemar_top1["holm_p"],
        "McNemar Best-pose (Holm p)": mcnemar_best["holm_p"],
        "Wilcoxon Top-1 (Holm p)": wilcoxon_top1["holm_p"],
        "Wilcoxon Best-pose (Holm p)": wilcoxon_best["holm_p"],
    })
    t3.to_csv(outdir / "table3_pairwise_comparisons.csv", index=False)
    print("\n=== Table 3: Holm-Adjusted Paired Comparisons ===")
    print(t3.to_string(index=False))

    # --- Table 4: coordinate bridge ------------------------------------
    t4 = table4_bridge_summary(df)
    t4.to_csv(outdir / "table4_bridge_summary.csv", index=False)
    print("\n=== Table 4: Coordinate-Bridge Performance ===")
    print(t4.to_string(index=False))

    agree_top1 = bridge_agreement(df, "Top1")
    agree_best = bridge_agreement(df, "Best")
    print("\n=== Bridge vs sPyRMSD Agreement ===")
    print(
        f"Top-1:      n = {agree_top1['n']}, rho = {agree_top1['rho']:.3f}, "
        f"mean |diff| = {agree_top1['mean_abs_diff']:.3f} A, "
        f"median |diff| = {agree_top1['median_abs_diff']:.3f} A"
    )
    print(
        f"Best-pose:  n = {agree_best['n']}, rho = {agree_best['rho']:.3f}, "
        f"mean |diff| = {agree_best['mean_abs_diff']:.3f} A, "
        f"median |diff| = {agree_best['median_abs_diff']:.3f} A"
    )

    pd.DataFrame([
        {"endpoint": "Top-1", **agree_top1},
        {"endpoint": "Best-pose", **agree_best},
    ]).to_csv(outdir / "bridge_agreement.csv", index=False)

    print(f"\nAll output tables written to: {outdir.resolve()}")


if __name__ == "__main__":
    sys.exit(main())
