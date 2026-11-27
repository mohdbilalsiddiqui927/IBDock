"""Shared utilities for the IBDock Astex Diverse Set benchmark analysis.

Provides data loading with validation, the statistical helpers used for
Tables 2-4 and S4, and the plotting palette shared by the figure script.

All success calls use the strict criterion  RMSD < 2.0 A  (the same rule as
the IBDock application); unresolved correspondences are never counted as
successes, which gives the conservative "all-85" endpoint of the manuscript.

SPDX-License-Identifier: MIT
"""
from __future__ import annotations

import csv
import itertools
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import stats

THRESHOLD = 2.0          # Angstrom, strict "<" (matches the IBDock application)
N_COMPLEXES = 85
# Order used in Table 2 and in every figure (performance-ranked in the manuscript)
STRATEGIES = ("Co-crystallized", "F-pocket", "P2Rank", "Blind")
LABELS = {"Co-crystallized": "Co-crystal", "F-pocket": "fpocket",
          "P2Rank": "P2Rank", "Blind": "Blind"}
LONG_LABELS = {"Co-crystallized": "Co-crystallized ligand", "F-pocket": "fpocket",
               "P2Rank": "P2Rank", "Blind": "Blind docking"}

# Strategy colours (hex). MID is used for lines/points, DARK for filled bars/violins;
# Figure 4b pairs MID (Top-1) with DARK (best pose).
PALETTE_MID = {"Co-crystallized": "#3771C8", "F-pocket": "#217844",
               "P2Rank": "#800033", "Blind": "#806600"}
PALETTE_DARK = {"Co-crystallized": "#2C5AA0", "F-pocket": "#16502D",
                "P2Rank": "#550022", "Blind": "#554400"}
NEUTRAL = "#C4C8B7"      # "resolved, not success" / "no pose generated"


@dataclass
class Benchmark:
    """Per-complex benchmark table; arrays are ordered by complex number.

    spy / bridge have shape (85, 2) = (Top-1, best pose); NaN = unresolved.
    """
    ids: list
    score: dict
    spy: dict
    bridge: dict

    def success(self, strategy: str, route: str = "spy", pose: int = 0) -> np.ndarray:
        """Boolean success vector (RMSD < THRESHOLD); NaN counts as failure."""
        arr = (self.spy if route == "spy" else self.bridge)[strategy][:, pose]
        with np.errstate(invalid="ignore"):
            return np.nan_to_num(arr, nan=np.inf) < THRESHOLD

    def resolved(self, strategy: str, route: str = "spy") -> np.ndarray:
        arr = (self.spy if route == "spy" else self.bridge)[strategy][:, 0]
        return ~np.isnan(arr)

    def common_resolved(self) -> np.ndarray:
        """Complexes with a numerical sPyRMSD under all four strategies."""
        return np.logical_and.reduce([self.resolved(s) for s in STRATEGIES])


def _float(x: str) -> float:
    x = x.strip()
    return float("nan") if x == "" else float(x)


def load_benchmark(path) -> Benchmark:
    """Read the tidy benchmark CSV (Table S2) and validate its structure."""
    with open(Path(path), newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if len(rows) != N_COMPLEXES * len(STRATEGIES):
        raise ValueError(f"expected {N_COMPLEXES * len(STRATEGIES)} rows, found {len(rows)}")
    by_pdb = {}
    for r in rows:
        by_pdb.setdefault((int(r["complex_no"]), r["pdb_id"]), {})[r["protocol"]] = r
    if len(by_pdb) != N_COMPLEXES:
        raise ValueError(f"expected {N_COMPLEXES} complexes, found {len(by_pdb)}")
    keys = sorted(by_pdb)
    ids = [k[1] for k in keys]
    score = {s: np.zeros(N_COMPLEXES) for s in STRATEGIES}
    spy = {s: np.full((N_COMPLEXES, 2), np.nan) for s in STRATEGIES}
    bridge = {s: np.full((N_COMPLEXES, 2), np.nan) for s in STRATEGIES}
    for idx, k in enumerate(keys):
        entry = by_pdb[k]
        if set(entry) != set(STRATEGIES):
            raise ValueError(f"{k[1]}: protocols {sorted(entry)} != {sorted(STRATEGIES)}")
        for s in STRATEGIES:
            r = entry[s]
            score[s][idx] = float(r["vina_score_kcal_mol"])
            spy[s][idx] = [_float(r["spyrmsd_top1_A"]), _float(r["spyrmsd_best_A"])]
            bridge[s][idx] = [_float(r["bridge_top1_A"]), _float(r["bridge_best_A"])]
            # status columns must agree with the presence of numerical values
            if (r["spyrmsd_status"] == "Numerical") == np.isnan(spy[s][idx, 0]):
                raise ValueError(f"{k[1]}/{s}: sPyRMSD status inconsistent with values")
            if (r["bridge_status"] == "Numerical") == np.isnan(bridge[s][idx, 0]):
                raise ValueError(f"{k[1]}/{s}: bridge status inconsistent with values")
    return Benchmark(ids, score, spy, bridge)


# ----------------------------------------------------------------- statistics
def wilson_ci(k: int, n: int, z: float = 1.959964) -> tuple:
    """Wilson score interval for a proportion, as fractions."""
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def holm_adjust(pvals) -> list:
    """Holm step-down adjusted P values (same order as input)."""
    pvals = np.asarray(pvals, float)
    order = np.argsort(pvals)
    adj = np.zeros(len(pvals))
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, min(1.0, (len(pvals) - rank) * pvals[idx]))
        adj[idx] = running
    return adj.tolist()


def cochran_q(binary: np.ndarray) -> tuple:
    """Cochran's Q for an (n_subjects, k_conditions) 0/1 matrix -> (Q, P)."""
    b = np.asarray(binary, float)
    k = b.shape[1]
    col, row, total = b.sum(0), b.sum(1), b.sum()
    q = (k - 1) * (k * (col ** 2).sum() - total ** 2) / (k * total - (row ** 2).sum())
    return float(q), float(stats.chi2.sf(q, k - 1))


def mcnemar_exact(a: np.ndarray, b: np.ndarray) -> float:
    """Exact two-sided McNemar P (binomial test on the discordant pairs)."""
    n01 = int(np.sum(a & ~b))
    n10 = int(np.sum(~a & b))
    n = n01 + n10
    return 1.0 if n == 0 else float(stats.binomtest(min(n01, n10), n, 0.5).pvalue)


def pairwise(values: dict, test) -> dict:
    """Holm-adjust `test(x, y)` over all six strategy pairs."""
    pairs = list(itertools.combinations(STRATEGIES, 2))
    raw = [test(values[a], values[b]) for a, b in pairs]
    return {f"{a}|{b}": p for (a, b), p in zip(pairs, holm_adjust(raw))}


def score_auc(bench: Benchmark, strategy: str, n_boot: int = 2000, seed: int = 0) -> dict:
    """AUC for 'recovered poses score better than failed poses' (Top-1, resolved
    cases only) with a percentile bootstrap 95% CI (fixed seed)."""
    ok = bench.resolved(strategy)
    rec_mask = ok & bench.success(strategy)
    rec = bench.score[strategy][rec_mask]
    non = bench.score[strategy][ok & ~rec_mask]

    def auc(x, y):
        return stats.mannwhitneyu(-x, -y)[0] / (len(x) * len(y))

    rng = np.random.default_rng(seed)
    boots = [auc(rng.choice(rec, len(rec)), rng.choice(non, len(non))) for _ in range(n_boot)]
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return {"auc": float(auc(rec, non)), "lo": float(lo), "hi": float(hi),
            "n_recovered": int(len(rec)), "n_failed": int(len(non))}
