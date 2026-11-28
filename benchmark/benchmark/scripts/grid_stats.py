"""IBDock Astex grid statistics: 3 exhaustiveness x 3 seeds x 4 search-space strategies.
Mirrors the manuscript's analysis (Wilson CI, Cochran Q, exact McNemar + Holm, Friedman + Wilcoxon + Holm)
and adds seed-reproducibility and exhaustiveness analyses.
Input : s3_new.pkl (parsed Table S3a-i).  Output: out/*.csv, out/stats.pkl
Success = RMSD <= 2.0 A (no value equals exactly 2.00 in the data, so < and <= coincide).
Primary endpoint: sPyRMSD route, all 85 complexes in the denominator, unresolved = non-success (as in the manuscript).
Sensitivity endpoint: sPyRMSD, else coordinate-bridge RMSD ("combined"); 2BR1 has neither and stays a non-success.
"""
import itertools, os, pickle
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.contingency_tables import cochrans_q, mcnemar
from statsmodels.stats.multitest import multipletests

OUT = "out"; os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20261113)
STRATS = ["Co-crystallized", "F-pocket", "P2Rank", "Blind"]
EXS, SEEDS = [8, 16, 32], [1, 42, 123]
df = pd.read_pickle("s3_new.pkl")
df["comb_t1"] = df.top1.fillna(df.cb_top1)
df["comb_bp"] = df.best.fillna(df.cb_best)
PDBS = sorted(df.pdb.unique()); N = len(PDBS)
assert N == 85 and df.groupby(["ex", "seed"]).size().eq(340).all()

def cube(col, thr=2.0):
    """binary array [ex, seed, strat, complex]; NaN -> 0 (non-success)"""
    a = np.zeros((3, 3, 4, N), dtype=int)
    for i, ex in enumerate(EXS):
        for j, sd in enumerate(SEEDS):
            g = df[(df.ex == ex) & (df.seed == sd)]
            for k, s in enumerate(STRATS):
                h = g[g.proto == s].set_index("pdb").reindex(PDBS)[col]
                a[i, j, k] = (h <= thr).fillna(False).astype(int).values
    return a

def rmsd_cube(col):
    a = np.full((3, 3, 4, N), np.nan)
    for i, ex in enumerate(EXS):
        for j, sd in enumerate(SEEDS):
            g = df[(df.ex == ex) & (df.seed == sd)]
            for k, s in enumerate(STRATS):
                a[i, j, k] = g[g.proto == s].set_index("pdb").reindex(PDBS)[col].values
    return a

ENDPOINTS = {"Top-1": ("top1", "comb_t1"), "Best-pose": ("best", "comb_bp")}
B = {(ep, kind): cube(c) for ep, cols in ENDPOINTS.items() for kind, c in zip(("primary", "combined"), cols)}
R = {(ep, "primary"): rmsd_cube(cols[0]) for ep, cols in ENDPOINTS.items()}

def wilson(k, n, z=1.96):
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * (c - h), 100 * (c + h)

def holm(p):
    return multipletests(p, method="holm")[1]

def mcnemar_exact(x, y):
    n10 = int(((x == 1) & (y == 0)).sum()); n01 = int(((x == 0) & (y == 1)).sum())
    if n10 + n01 == 0: return 1.0, n01, n10
    return mcnemar([[0, n10], [n01, 0]], exact=True).pvalue, n01, n10

def cq(a):  # a: complexes x raters
    if len(np.unique(a)) < 2 or (a.sum(1) % a.shape[1] == 0).all():
        return np.nan, 1.0
    r = cochrans_q(a, return_object=True); return r.statistic, r.pvalue

res = {}
# ---------- 1. per-condition rates (all-85, Wilson) ----------
rows = []
for kind in ("primary", "combined"):
    for ep in ENDPOINTS:
        a = B[(ep, kind)]
        for i, ex in enumerate(EXS):
            for j, sd in enumerate(SEEDS):
                for k, s in enumerate(STRATS):
                    kk = int(a[i, j, k].sum()); lo, hi = wilson(kk, N)
                    rows.append([kind, ep, ex, sd, s, kk, N, 100 * kk / N, lo, hi])
res["rates"] = pd.DataFrame(rows, columns=["endpoint_def", "endpoint", "exhaustiveness", "seed", "strategy", "k", "n", "pct", "wilson_lo", "wilson_hi"])

# ---------- 2. seed-aggregated rates (mean/SD across seeds + cluster bootstrap over complexes) ----------
def boot_mean_ci(x, reps=10000):
    """x: complexes -> per-complex seed-mean success; percentile bootstrap over complexes"""
    idx = rng.integers(0, len(x), (reps, len(x))); m = x[idx].mean(1) * 100
    return np.percentile(m, [2.5, 97.5])
rows = []
for kind in ("primary", "combined"):
    for ep in ENDPOINTS:
        a = B[(ep, kind)]
        for i, ex in enumerate(EXS):
            for k, s in enumerate(STRATS):
                per_seed = a[i, :, k].sum(1) * 100 / N
                cm = a[i, :, k].mean(0)  # per-complex success fraction over seeds
                lo, hi = boot_mean_ci(cm)
                rows.append([kind, ep, ex, s, per_seed.mean(), per_seed.std(ddof=1), per_seed.min(), per_seed.max(), lo, hi])
res["agg"] = pd.DataFrame(rows, columns=["endpoint_def", "endpoint", "exhaustiveness", "strategy", "mean_pct", "sd_pct", "min_pct", "max_pct", "boot_lo", "boot_hi"])

# ---------- 3. strategy effect: per condition (paper's design) ----------
rows = []
pairs = list(itertools.combinations(range(4), 2))
def strat_tests(a_cs, label):  # a_cs: [strat, complex] binary
    out = []
    Q, p = cq(a_cs.T)
    pm = [mcnemar_exact(a_cs[x], a_cs[y]) for x, y in pairs]
    ph = holm([m[0] for m in pm])
    for (x, y), m, h in zip(pairs, pm, ph):
        out.append([label, Q, p, f"{STRATS[x]} vs {STRATS[y]}", m[1], m[2], m[0], h])
    return out
for ep in ENDPOINTS:
    a = B[(ep, "primary")]
    for i, ex in enumerate(EXS):
        for j, sd in enumerate(SEEDS):
            for r in strat_tests(a[i, j], f"ex{ex}_seed{sd}"):
                rows.append([ep, ex, sd] + r[1:])
        cons = (a[i].sum(0) >= 2).astype(int)  # seed-consensus: success in >=2 of 3 seeds
        for r in strat_tests(cons, f"ex{ex}_consensus"):
            rows.append([ep, ex, "consensus"] + r[1:])
res["strategy_tests"] = pd.DataFrame(rows, columns=["endpoint", "exhaustiveness", "seed", "Q", "Q_p", "comparison", "gain_n01", "loss_n10", "McNemar_p", "Holm_p"])

# RMSD route: Friedman + Wilcoxon + Holm on complexes numeric under all 4 strategies (median over seeds)
rows = []
for ep in ENDPOINTS:
    r_ = R[(ep, "primary")]
    for i, ex in enumerate(EXS):
        ok = ~np.isnan(r_[i]).any(axis=(0, 1))  # numeric in all seeds & strategies
        med = np.nanmedian(r_[i][:, :, ok], axis=0)  # strat x complexes
        fr = stats.friedmanchisquare(*med)
        pw = [stats.wilcoxon(med[x], med[y]).pvalue for x, y in pairs]
        ph = holm(pw)
        for (x, y), p_, h in zip(pairs, pw, ph):
            rows.append([ep, ex, int(ok.sum()), fr.statistic, fr.pvalue, f"{STRATS[x]} vs {STRATS[y]}", p_, h])
res["rmsd_tests"] = pd.DataFrame(rows, columns=["endpoint", "exhaustiveness", "n_complexes", "Friedman_Q", "Friedman_p", "comparison", "Wilcoxon_p", "Holm_p"])

rows = []
for ep in ENDPOINTS:
    r_ = R[(ep, "primary")]
    for i, ex in enumerate(EXS):
        ok = ~np.isnan(r_[i]).any(axis=(0, 1))
        for k, s in enumerate(STRATS):
            v = np.nanmedian(r_[i][:, k, :][:, ok], axis=0)
            rows.append([ep, ex, s, int(ok.sum()), np.median(v), *np.percentile(v, [25, 75])])
res["rmsd_summary"] = pd.DataFrame(rows, columns=["endpoint", "exhaustiveness", "strategy", "n", "median", "q1", "q3"])

# ---------- 4. seed reproducibility ----------
rows = []
def fleiss_kappa(a):  # a: complexes x raters binary
    n, m = a.shape; p1 = a.mean(); pe = p1 ** 2 + (1 - p1) ** 2
    Pi = ((a.sum(1) * (a.sum(1) - 1)) + ((m - a.sum(1)) * (m - a.sum(1) - 1))) / (m * (m - 1))
    return (Pi.mean() - pe) / (1 - pe) if pe < 1 else np.nan
for ep in ENDPOINTS:
    a = B[(ep, "primary")]
    for i, ex in enumerate(EXS):
        pv = []
        for k, s in enumerate(STRATS):
            x = a[i, :, k].T; Q, p = cq(x); pv.append(p)
            rows.append([ep, ex, s, Q, p, (x.min(1) == x.max(1)).mean() * 100, int((x.min(1) != x.max(1)).sum()), fleiss_kappa(x)])
        ph = holm(pv)
        for t_, h in zip(range(4), ph): rows[-4 + t_].append(h)
        x = a[i].transpose(2, 1, 0).reshape(-1, 3)  # all 340 cells x seeds
        Q, p = cq(x)
        rows.append([ep, ex, "All strategies (340 cells, descriptive)", Q, p, (x.min(1) == x.max(1)).mean() * 100, int((x.min(1) != x.max(1)).sum()), fleiss_kappa(x), np.nan])
res["seed_repro"] = pd.DataFrame(rows, columns=["endpoint", "exhaustiveness", "strategy", "Q_seeds", "Q_p", "pct_cells_identical", "n_cells_flipping", "fleiss_kappa", "Holm_p_over_strategies"])

# RMSD agreement between seeds (Top-1 numeric cells)
rows = []
r_ = R[("Top-1", "primary")]
for i, ex in enumerate(EXS):
    dif, rho = [], []
    for j1, j2 in itertools.combinations(range(3), 2):
        x = r_[i, j1].ravel(); y = r_[i, j2].ravel(); ok = ~np.isnan(x) & ~np.isnan(y)
        dif.append(np.abs(x[ok] - y[ok])); rho.append(stats.spearmanr(x[ok], y[ok]).statistic)
    d = np.concatenate(dif)
    rows.append([ex, len(d) // 3, np.mean(rho), np.median(d), np.percentile(d, 90), (d > 2).mean() * 100, (d < 0.5).mean() * 100])
res["rmsd_seed_agreement"] = pd.DataFrame(rows, columns=["exhaustiveness", "cells_per_pair", "mean_spearman", "median_abs_diff_A", "p90_abs_diff_A", "pct_pairs_diff_gt_2A", "pct_pairs_diff_lt_0.5A"])

# ---------- 5. exhaustiveness effect ----------
rows = []; lev = [(0, 1), (1, 2), (0, 2)]; names = {(0, 1): "8 -> 16", (1, 2): "16 -> 32", (0, 2): "8 -> 32"}
for ep in ENDPOINTS:
    a = B[(ep, "primary")]
    for k, s in enumerate(STRATS + ["All strategies"]):
        # per-complex seed-mean success fraction (complexes x ex); for All: average over strategies
        if s == "All strategies": m = a.mean(axis=1).mean(axis=1)  # ex x complex
        else: m = a[:, :, k].mean(axis=1)
        for (u, v) in lev:
            d = m[v] - m[u]
            idx = rng.integers(0, N, (10000, N)); bs = d[idx].mean(1) * 100
            lo, hi = np.percentile(bs, [2.5, 97.5])
            nz = d[d != 0]
            wp = stats.wilcoxon(nz).pvalue if len(nz) else 1.0
            rows.append([ep, s, names[(u, v)], d.mean() * 100, lo, hi, wp, np.nan])
tmp = pd.DataFrame(rows, columns=["endpoint", "strategy", "exhaustiveness_step", "delta_pct_points", "boot_lo", "boot_hi", "Wilcoxon_p", "Holm_p"])
# Holm within endpoint over the 4 single strategies x 3 steps (12 tests)
for ep in ENDPOINTS:
    msk = (tmp.endpoint == ep) & (tmp.strategy != "All strategies")
    tmp.loc[msk, "Holm_p"] = holm(tmp.loc[msk, "Wilcoxon_p"].values)
res["exh_effect"] = tmp

# per-seed Cochran Q across exhaustiveness (n=85 per strategy), consensus too
rows = []
for ep in ENDPOINTS:
    a = B[(ep, "primary")]
    for k, s in enumerate(STRATS):
        for j, sd in enumerate(SEEDS):
            Q, p = cq(a[:, j, k].T); rows.append([ep, s, sd, Q, p])
        cons = (a[:, :, k].sum(1) >= 2).astype(int)
        Q, p = cq(cons.T); rows.append([ep, s, "consensus", Q, p])
res["exh_Q"] = pd.DataFrame(rows, columns=["endpoint", "strategy", "seed", "Q_exhaustiveness", "Q_p"])

# interaction: does the 8 -> 32 gain differ between strategies? (Friedman over strategies on per-complex gain)
rows = []
for ep in ENDPOINTS:
    a = B[(ep, "primary")]
    g = [a[2, :, k].mean(0) - a[0, :, k].mean(0) for k in range(4)]
    fr = stats.friedmanchisquare(*g)
    rows.append([ep, fr.statistic, fr.pvalue] + [np.mean(x) * 100 for x in g])
res["interaction"] = pd.DataFrame(rows, columns=["endpoint", "Friedman_Q", "Friedman_p"] + [f"gain_{s}" for s in STRATS])

# ---------- 6. sensitivity: strategy ordering & Q with combined endpoint ----------
rows = []
for ep in ENDPOINTS:
    a = B[(ep, "combined")]
    for i, ex in enumerate(EXS):
        for j, sd in enumerate(SEEDS):
            Q, p = cq(a[i, j]); rows.append([ep, ex, sd, Q, p, " > ".join(STRATS[k] for k in np.argsort(-a[i, j].sum(1), kind="stable"))])
res["combined_strategy_Q"] = pd.DataFrame(rows, columns=["endpoint", "exhaustiveness", "seed", "Q", "Q_p", "order"])

# ---------- 7. resolution map ----------
rs = df[(df.ex == 16) & (df.seed == 42)].copy()
rs["spy"] = rs.top1.notna(); rs["bridge"] = rs.cb_top1.notna()
res["resolution"] = rs.groupby("proto")[["spy", "bridge"]].sum().assign(either=rs.assign(e=rs.spy | rs.bridge).groupby("proto").e.sum(), n=85)

pickle.dump(res, open(f"{OUT}/stats.pkl", "wb"))
for k, v in res.items(): v.to_csv(f"{OUT}/{k}.csv", index=False)
print({k: v.shape for k, v in res.items()})

np.savez(f"{OUT}/cubes.npz", top1=B[("Top-1", "primary")], best=B[("Best-pose", "primary")],
         r_top1=R[("Top-1", "primary")], r_best=R[("Best-pose", "primary")], pdbs=np.array(PDBS))
