#!/usr/bin/env python3
"""Statistical analysis of the IBDock Astex Diverse Set benchmark.

Reproduces every statistic reported in the manuscript and Supporting
Information from the per-complex benchmark table (data/astex_benchmark_S2.csv,
identical to Table S2):

  * Table 2   pose recovery (all-85 and resolved-case rates, Wilson 95% CIs)
  * Table 3   Cochran's Q + exact McNemar (Holm), Friedman + Wilcoxon (Holm)
  * Table 4   coordinate-bridge applicability, success and agreement
  * Table S4  binary success / unresolved inputs, final statistical summary
  * Docking-score analysis behind Figures S4a/b (Holm Wilcoxon, AUC + bootstrap)

Usage
    python analyze_astex.py                       # writes ./results
    python analyze_astex.py --check               # also verify the headline numbers
    python analyze_astex.py --data my.csv --out out_dir

Randomness: only the bootstrap CI of the AUC uses random numbers (fixed seed 0).

SPDX-License-Identifier: MIT
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import scipy
from scipy import stats

from ibdock_common import (LONG_LABELS, N_COMPLEXES, STRATEGIES, THRESHOLD,
                           Benchmark, cochran_q, load_benchmark, mcnemar_exact,
                           pairwise, score_auc, wilson_ci)

POSES = ("top1", "best")
# Column order of Table S4 in the Supporting Information
S4_ORDER = ("Co-crystallized", "P2Rank", "F-pocket", "Blind")


def analyze(b: Benchmark) -> dict:
    R: dict = {"threshold_A": THRESHOLD, "n_complexes": N_COMPLEXES}

    # ---- Table 2: pose recovery -------------------------------------------------
    t2 = {}
    for s in STRATEGIES:
        res = int(b.resolved(s).sum())
        row = {"resolved": res, "unresolved": N_COMPLEXES - res}
        for i, pose in enumerate(POSES):
            k = int(b.success(s, "spy", i).sum())
            lo, hi = wilson_ci(k, N_COMPLEXES)
            row[pose] = {"successes": k, "pct_all85": 100 * k / N_COMPLEXES,
                         "ci95_pct": [100 * lo, 100 * hi], "pct_resolved": 100 * k / res}
        t2[s] = row
    R["table2"] = t2

    # ---- Table 3 (binary): Cochran's Q and exact McNemar with Holm --------------
    R["binary"] = {}
    for i, pose in enumerate(POSES):
        B = np.column_stack([b.success(s, "spy", i) for s in STRATEGIES]).astype(int)
        q, p = cochran_q(B)
        R["binary"][pose] = {
            "cochran_q": q, "df": len(STRATEGIES) - 1, "p": p,
            "mcnemar_holm": pairwise({s: b.success(s, "spy", i) for s in STRATEGIES}, mcnemar_exact)}

    # ---- Table 3 (continuous): common-resolved subset ---------------------------
    common = b.common_resolved()
    R["common_subset"] = {"n": int(common.sum()),
                          "excluded": [b.ids[j] for j in np.where(~common)[0]]}
    R["continuous"] = {}
    for i, pose in enumerate(POSES):
        vals = {s: b.spy[s][common, i] for s in STRATEGIES}
        q, p = stats.friedmanchisquare(*[vals[s] for s in STRATEGIES])
        R["continuous"][pose] = {
            "median_q1_q3": {s: [float(x) for x in np.percentile(vals[s], [50, 25, 75])]
                             for s in STRATEGIES},
            "friedman_q": float(q), "friedman_p": float(p),
            "wilcoxon_holm": pairwise(vals, lambda x, y: stats.wilcoxon(x, y).pvalue)}

    # ---- Table 4: coordinate bridge --------------------------------------------
    t4 = {}
    for s in STRATEGIES:
        nb = int(b.resolved(s, "bridge").sum())
        both = int((b.resolved(s) & b.resolved(s, "bridge")).sum())
        row = {"bridge_numerical": nb, "bridge_unresolved": N_COMPLEXES - nb, "both_numerical": both}
        for i, pose in enumerate(POSES):
            k = int(b.success(s, "bridge", i).sum())
            row[pose] = {"successes": k, "pct_all85": 100 * k / N_COMPLEXES,
                         "pct_bridge_resolved": 100 * k / nb}
        t4[s] = row
    R["table4"] = t4

    # ---- agreement between the two routes ---------------------------------------
    R["agreement"] = {}
    for i, pose in enumerate(POSES):
        x, y, tags = [], [], []
        for s in STRATEGIES:
            m = b.resolved(s) & b.resolved(s, "bridge")
            x.extend(b.spy[s][m, i]); y.extend(b.bridge[s][m, i])
            tags.extend((b.ids[j], s) for j in np.where(m)[0])
        x, y = np.array(x), np.array(y)
        d = np.abs(x - y)
        R["agreement"][pose] = {
            "n_pairs": int(len(x)), "spearman_rho": float(stats.spearmanr(x, y)[0]),
            "mean_abs_diff_A": float(d.mean()), "median_abs_diff_A": float(np.median(d)),
            "max_abs_diff_A": float(d.max()), "within_0.05_A_pct": float(100 * np.mean(d <= 0.05)),
            "classification_disagreements": [
                {"complex": t[0], "strategy": t[1], "spyrmsd": float(a), "bridge": float(c)}
                for t, a, c in zip(tags, x, y) if (a < THRESHOLD) != (c < THRESHOLD)]}

    # ---- unresolved-case bookkeeping (Table S1 header numbers) -------------------
    unres = {s: {b.ids[j] for j in np.where(~b.resolved(s))[0]} for s in STRATEGIES}
    shared = set.intersection(*unres.values())
    bridge_un = {s: {b.ids[j] for j in np.where(~b.resolved(s, "bridge"))[0]} for s in STRATEGIES}
    bu = set.union(*bridge_un.values())
    R["unresolved"] = {
        "spyrmsd_counts": {s: len(unres[s]) for s in STRATEGIES},
        "spyrmsd_shared_by_all": sorted(shared),
        "spyrmsd_extra_per_strategy": {s: sorted(unres[s] - shared) for s in STRATEGIES},
        "bridge_unresolved": sorted(bu),
        "bridge_unresolved_identical_across_strategies": all(v == bu for v in bridge_un.values()),
        "unique_complexes_either_route": len(set.union(*unres.values()) | bu),
        "unresolved_by_both_routes": sorted(set.union(*unres.values()) & bu)}

    # ---- docking-score analysis (Figures S4a/b) ---------------------------------
    R["scores"] = {
        "mean_kcal_mol": {s: float(b.score[s].mean()) for s in STRATEGIES},
        "wilcoxon_holm": pairwise(b.score, lambda x, y: stats.wilcoxon(x, y).pvalue),
        "auc": {s: score_auc(b, s) for s in STRATEGIES}}
    return R


def s4_table(b: Benchmark) -> list:
    """Rows of Table S4: 1 = RMSD < 2.0 A, 0 = failure or unresolved."""
    cols = {(s, i): b.success(s, "spy", i).astype(int) for s in S4_ORDER for i in (0, 1)}
    return [[b.ids[j]] + [int(cols[(s, i)][j]) for s in S4_ORDER for i in (0, 1)]
            for j in range(N_COMPLEXES)]


def sci(p: float) -> str:
    """P value as in Table 3: four significant figures, scientific below 1e-3."""
    return f"{p:.4g}" if p >= 1e-3 else f"{p:.3e}"


def _csv(path: Path):
    return open(path, "w", newline="", encoding="utf-8")


def write_outputs(R: dict, b: Benchmark, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "statistics_summary.json", "w", encoding="utf-8") as f:
        json.dump(R, f, indent=2)
    with _csv(out / "tableS4_binary_inputs.csv") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["pdb_id"] + [f"{s}_{p}" for s in S4_ORDER for p in ("top1", "best")])
        w.writerows(s4_table(b))
    with _csv(out / "table2_pose_recovery.csv") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["strategy", "sPyRMSD_resolved", "top1_resolved", "top1_resolved_pct", "top1_all85",
                    "top1_all85_pct", "top1_ci95", "best_resolved", "best_resolved_pct", "best_all85",
                    "best_all85_pct", "best_ci95"])
        for s in STRATEGIES:
            r = R["table2"][s]; t, bp = r["top1"], r["best"]
            w.writerow([LONG_LABELS[s], f'{r["resolved"]}/85', f'{t["successes"]}/{r["resolved"]}',
                        f'{t["pct_resolved"]:.1f}', f'{t["successes"]}/85', f'{t["pct_all85"]:.1f}',
                        f'{t["ci95_pct"][0]:.1f}-{t["ci95_pct"][1]:.1f}', f'{bp["successes"]}/{r["resolved"]}',
                        f'{bp["pct_resolved"]:.1f}', f'{bp["successes"]}/85', f'{bp["pct_all85"]:.1f}',
                        f'{bp["ci95_pct"][0]:.1f}-{bp["ci95_pct"][1]:.1f}'])
    with _csv(out / "table3_pairwise_holm.csv") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["comparison", "mcnemar_top1", "mcnemar_best", "wilcoxon_top1", "wilcoxon_best"])
        for k in R["binary"]["top1"]["mcnemar_holm"]:
            w.writerow([k.replace("|", " vs "), sci(R["binary"]["top1"]["mcnemar_holm"][k]),
                        sci(R["binary"]["best"]["mcnemar_holm"][k]),
                        sci(R["continuous"]["top1"]["wilcoxon_holm"][k]),
                        sci(R["continuous"]["best"]["wilcoxon_holm"][k])])
    with _csv(out / "table4_bridge.csv") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["metric"] + [LONG_LABELS[s] for s in STRATEGIES])
        t4 = R["table4"]
        w.writerow(["bridge_numerical"] + [f'{t4[s]["bridge_numerical"]}/85' for s in STRATEGIES])
        w.writerow(["bridge_unresolved"] + [f'{t4[s]["bridge_unresolved"]}/85' for s in STRATEGIES])
        for pose in POSES:
            w.writerow([f"{pose}_success_all85"] + [f'{t4[s][pose]["successes"]}/85 ({t4[s][pose]["pct_all85"]:.1f}%)' for s in STRATEGIES])
        for pose in POSES:
            w.writerow([f"{pose}_success_bridge_resolved"] + [f'{t4[s][pose]["successes"]}/{t4[s]["bridge_numerical"]} ({t4[s][pose]["pct_bridge_resolved"]:.1f}%)' for s in STRATEGIES])
        w.writerow(["both_methods_numerical"] + [f'{t4[s]["both_numerical"]}/85' for s in STRATEGIES])

    L = []
    add = L.append
    add("IBDock Astex benchmark - statistical summary")
    add(f"scipy {scipy.__version__}, numpy {np.__version__}; success = RMSD < {THRESHOLD} A; unresolved = failure\n")
    for pose in POSES:
        add(f"[{pose}] all-85 successes: " + ", ".join(f'{LONG_LABELS[s]} {R["table2"][s][pose]["successes"]}/85 ({R["table2"][s][pose]["pct_all85"]:.1f}%)' for s in STRATEGIES))
        bn = R["binary"][pose]
        add(f"[{pose}] Cochran Q = {bn['cochran_q']:.3f} (df {bn['df']}), P = {sci(bn['p'])}")
        add(f"[{pose}] McNemar-Holm: " + "; ".join(f"{k.replace('|',' vs ')} {sci(v)}" for k, v in bn["mcnemar_holm"].items()))
    add(f"\ncommon-resolved subset: n = {R['common_subset']['n']} (excluded: {', '.join(R['common_subset']['excluded'])})")
    for pose in POSES:
        c = R["continuous"][pose]
        add(f"[{pose}] median (Q1-Q3) A: " + "; ".join(f"{LONG_LABELS[s]} {v[0]:.2f} ({v[1]:.2f}-{v[2]:.2f})" for s, v in c["median_q1_q3"].items()))
        add(f"[{pose}] Friedman Q = {c['friedman_q']:.3f}, P = {sci(c['friedman_p'])}")
        add(f"[{pose}] Wilcoxon-Holm: " + "; ".join(f"{k.replace('|',' vs ')} {sci(v)}" for k, v in c["wilcoxon_holm"].items()))
    add("")
    for pose in POSES:
        a = R["agreement"][pose]
        add(f"[{pose}] bridge vs sPyRMSD: n = {a['n_pairs']}, rho = {a['spearman_rho']:.3f}, mean |d| = {a['mean_abs_diff_A']:.3f} A, "
            f"median |d| = {a['median_abs_diff_A']:.2f} A, max |d| = {a['max_abs_diff_A']:.2f} A, disagreements = {len(a['classification_disagreements'])} "
            f"({', '.join(sorted({d['complex'] for d in a['classification_disagreements']}))})")
    u = R["unresolved"]
    add(f"\nunresolved sPyRMSD: {u['spyrmsd_counts']}; shared by all four = {len(u['spyrmsd_shared_by_all'])}; extras = {u['spyrmsd_extra_per_strategy']}")
    add(f"bridge unresolved: {u['bridge_unresolved']} (identical across strategies: {u['bridge_unresolved_identical_across_strategies']}); "
        f"unique complexes unresolved by either route = {u['unique_complexes_either_route']}; by both = {u['unresolved_by_both_routes']}")
    add("\ndocking scores: mean kcal/mol " + ", ".join(f"{LONG_LABELS[s]} {v:.2f}" for s, v in R["scores"]["mean_kcal_mol"].items()))
    add("score Wilcoxon-Holm: " + "; ".join(f"{k.replace('|',' vs ')} {sci(v)}" for k, v in R["scores"]["wilcoxon_holm"].items()))
    add("score AUC [95% CI]: " + "; ".join(f"{LONG_LABELS[s]} {v['auc']:.2f} [{v['lo']:.2f}, {v['hi']:.2f}] (n {v['n_recovered']}/{v['n_failed']})" for s, v in R["scores"]["auc"].items()))
    (out / "statistics_summary.txt").write_text("\n".join(L) + "\n", encoding="utf-8")


# Regression guard: headline values that the manuscript/SI must show (benchmark
# build with the Meeko ligand-preparation fallback). Fails loudly if the input
# table or the analysis code changes.
EXPECTED = {
    "top1_successes": {"Co-crystallized": 51, "F-pocket": 44, "P2Rank": 31, "Blind": 20},
    "best_successes": {"Co-crystallized": 61, "F-pocket": 54, "P2Rank": 36, "Blind": 27},
    "resolved": {"Co-crystallized": 70, "F-pocket": 69, "P2Rank": 69, "Blind": 69},
    "common_subset_n": 69, "paired_n": 265,
    "bridge_numerical_each": 81, "unique_unresolved_either_route": 19,
    "disagreements": {"top1": 4, "best": 4},
    "rho": {"top1": 0.948, "best": 0.944}, "mae": {"top1": 0.115, "best": 0.097},
}


def check(R: dict, b: Benchmark) -> list:
    fails = []

    def eq(name, got, want, tol=0.0):
        bad = abs(got - want) > tol if isinstance(want, float) else got != want
        if bad:
            fails.append(f"{name}: got {got}, expected {want}")

    for s in STRATEGIES:
        eq(f"top1 {s}", R["table2"][s]["top1"]["successes"], EXPECTED["top1_successes"][s])
        eq(f"best {s}", R["table2"][s]["best"]["successes"], EXPECTED["best_successes"][s])
        eq(f"resolved {s}", R["table2"][s]["resolved"], EXPECTED["resolved"][s])
        eq(f"bridge numerical {s}", R["table4"][s]["bridge_numerical"], EXPECTED["bridge_numerical_each"])
    eq("common subset", R["common_subset"]["n"], EXPECTED["common_subset_n"])
    eq("paired n", R["agreement"]["top1"]["n_pairs"], EXPECTED["paired_n"])
    eq("unique unresolved", R["unresolved"]["unique_complexes_either_route"], EXPECTED["unique_unresolved_either_route"])
    for pose in POSES:
        a = R["agreement"][pose]
        eq(f"disagreements {pose}", len(a["classification_disagreements"]), EXPECTED["disagreements"][pose])
        eq(f"rho {pose}", round(a["spearman_rho"], 3), EXPECTED["rho"][pose], 5e-4)
        eq(f"MAE {pose}", round(a["mean_abs_diff_A"], 3), EXPECTED["mae"][pose], 5e-4)
    col = [sum(r[1 + 2 * k + i] for r in s4_table(b)) for k in range(4) for i in (0, 1)]
    eq("Table S4 column sums", col, [51, 61, 31, 36, 44, 54, 20, 27])
    return fails


def main(argv=None) -> int:
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data", type=Path, default=here / "data" / "astex_benchmark_S2.csv")
    ap.add_argument("--out", type=Path, default=Path("results"))
    ap.add_argument("--check", action="store_true", help="verify headline numbers; exit 1 on mismatch")
    args = ap.parse_args(argv)
    bench = load_benchmark(args.data)
    R = analyze(bench)
    write_outputs(R, bench, args.out)
    print(f"wrote results to {args.out.resolve()}")
    print((args.out / "statistics_summary.txt").read_text(encoding="utf-8"))
    if args.check:
        fails = check(R, bench)
        print("CHECK:", "PASS (all headline numbers reproduced)" if not fails else "FAIL\n  " + "\n  ".join(fails))
        return 1 if fails else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
