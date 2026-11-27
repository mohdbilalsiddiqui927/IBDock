#!/usr/bin/env python3
"""Regenerate the benchmark figures of the IBDock manuscript and Supporting Information.

  Figure 4a  success rate (Top-1 / best pose) with Wilson 95% CI and unresolved cases
  Figure 4b  sPyRMSD distributions, complexes resolved under all four strategies
  Figure 5a/b  coordinate-bridge vs sPyRMSD agreement (Top-1 / best pose)
  Figure S1a/b cumulative RMSD curves (Top-1 / best pose)
  Figure S2  outcome decomposition (recovered / not ranked first / not sampled / unresolved)
  Figure S3  applicability of the two correspondence routes
  Figure S4a/b docking score by strategy; score vs pose recovery (AUC)

Figures 1-3 and 6 are schematics/composites and are not generated here.

Usage
    python make_figures.py                          # all figures -> ./figures (SVG + PNG)
    python make_figures.py --only 4a 5a --dpi 600   # default PNG resolution is 1000 dpi
    python make_figures.py --formats svg

Outputs are deterministic (fixed RNG seeds, no timestamps in the files).
SVG text is kept as editable text; the font family is written as the stack
'Palatino Linotype', 'TeX Gyre Pagella', Palatino, 'Book Antiqua', serif so the
figures render in Palatino on any machine that has one of these fonts.

SPDX-License-Identifier: MIT
"""
from __future__ import annotations

import argparse
import re
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, MultipleLocator, NullLocator
from scipy import stats

from ibdock_common import (LABELS, NEUTRAL, PALETTE_DARK, PALETTE_MID, STRATEGIES,
                           THRESHOLD, N_COMPLEXES, Benchmark, load_benchmark, pairwise,
                           score_auc, wilson_ci)

FONT_CANDIDATES = ["Palatino Linotype", "TeX Gyre Pagella", "Palatino", "Book Antiqua"]
CSS_STACK = "'Palatino Linotype','TeX Gyre Pagella',Palatino,'Book Antiqua',serif"
LOG_TICKS = [0.2, 0.5, 1, 2, 5, 10, 20, 50]
LOG_TICK_LABELS = ["0.2", "0.5", "1", "2", "5", "10", "20", "50"]
MINUS = lambda t: t.replace("-", "\u2212")  # typographic minus for axis labels
ANGSTROM = "\u00c5"


def setup_style() -> str:
    """Choose the first installed Palatino-class font and set matplotlib defaults."""
    installed = {f.name for f in font_manager.fontManager.ttflist}
    family = next((f for f in FONT_CANDIDATES if f in installed), None)
    if family is None:
        family = "DejaVu Serif"
        warnings.warn("No Palatino-class font found; falling back to DejaVu Serif "
                      "(layout will differ slightly). Install TeX Gyre Pagella or Palatino Linotype.")
    plt.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "svg.hashsalt": "ibdock-astex-figures",
        "font.family": family,
        "font.size": 9, "axes.labelsize": 9.6, "xtick.labelsize": 8.6, "ytick.labelsize": 9,
        "legend.fontsize": 8, "axes.linewidth": 0.8, "xtick.major.width": 0.8,
        "ytick.major.width": 0.8, "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.direction": "out", "ytick.direction": "out", "legend.frameon": False,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.facecolor": "white", "savefig.facecolor": "white"})
    return family


def save(fig, outdir: Path, name: str, formats, dpi: int, family: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    if "svg" in formats:
        path = outdir / f"{name}.svg"
        fig.savefig(path, bbox_inches="tight", pad_inches=0.06, metadata={"Date": None})
        text = path.read_text(encoding="utf-8")
        text = re.sub(r"font-family:\s*'%s'" % re.escape(family), "font-family: " + CSS_STACK, text)
        path.write_text(text, encoding="utf-8")
    if "png" in formats:
        fig.savefig(outdir / f"{name}.png", bbox_inches="tight", pad_inches=0.06, dpi=dpi,
                    metadata={"Software": None})
    plt.close(fig)


def tint(color: str, a: float):
    return tuple(np.array(to_rgb(color)) * a + (1 - a))


def panel_letter(ax, letter: str, dx: float, dy: float = 1.055) -> None:
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=14, fontweight="bold")


def log_axis(ax, which: str = "y") -> None:
    axis = ax.yaxis if which == "y" else ax.xaxis
    axis.set_major_locator(FixedLocator(LOG_TICKS)); axis.set_minor_locator(NullLocator())
    (ax.set_yticklabels if which == "y" else ax.set_xticklabels)(LOG_TICK_LABELS)


# --------------------------------------------------------------------- Figure 4
def fig4a(b: Benchmark):
    fig, ax = plt.subplots(figsize=(7.0, 4.1))
    gw = 0.78
    bw = gw / 2 * 0.92
    for g, s in enumerate(STRATEGIES):
        res = int(b.resolved(s).sum())
        for i, pose in enumerate(("Top-1", "Best")):
            k = int(b.success(s, "spy", i).sum())
            x = g - gw / 2 + bw / 2 + i * bw * 1.08
            p = 100 * k / N_COMPLEXES
            lo, hi = wilson_ci(k, N_COMPLEXES)
            ax.bar(x, p, width=bw, color=PALETTE_DARK[s], edgecolor="black", linewidth=0.8, zorder=3)
            ax.bar(x, 100 * (res - k) / N_COMPLEXES, bottom=p, width=bw, color=NEUTRAL,
                   edgecolor="black", linewidth=0.7, zorder=2)
            ax.bar(x, 100 * (N_COMPLEXES - res) / N_COMPLEXES, bottom=100 * res / N_COMPLEXES, width=bw,
                   color="white", hatch="////", edgecolor="#555555", linewidth=0.7, zorder=2)
            ax.errorbar(x, p, yerr=[[p - 100 * lo], [100 * hi - p]], fmt="none", ecolor="black",
                        elinewidth=1, capsize=2.8, capthick=1, zorder=4)
            ax.text(x, 100 * hi + 4.5, f"{p:.1f}", ha="center", fontsize=7.4, zorder=5,
                    bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.9))
            ax.text(x, -7.2, pose, ha="center", fontsize=7.4)
    ax.set_xticks(range(4)); ax.set_xticklabels([LABELS[s] for s in STRATEGIES], fontsize=9.4)
    ax.tick_params(axis="x", pad=18)
    ax.set_ylim(0, 100); ax.yaxis.set_major_locator(MultipleLocator(20))
    ax.set_ylabel("% of all 85 complexes")
    panel_letter(ax, "a", -0.075)
    ax.legend(handles=[Patch(fc="white", ec="black", lw=0.8, label=f"Success (RMSD < {THRESHOLD:.1f} {ANGSTROM})"),
                       Patch(fc=NEUTRAL, ec="black", lw=0.8, label="Resolved, not success"),
                       Patch(fc="white", ec="#555555", hatch="////", lw=0.8, label="sPyRMSD unresolved")],
              loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=3, fontsize=8,
              handlelength=1.4, columnspacing=1.2)
    return fig


def fig4b(b: Benchmark):
    common = b.common_resolved()
    fig, ax = plt.subplots(figsize=(7.0, 4.3))
    rng = np.random.default_rng(7)
    pos, ticks, labels, centers = 0.0, [], [], []
    for s in STRATEGIES:
        start = pos
        for i, lab in enumerate(("Top-1", "Best")):
            d = b.spy[s][common, i]
            bp = ax.boxplot([d], positions=[pos], widths=0.6, patch_artist=True, showfliers=False,
                            medianprops=dict(color="black", lw=1.3), whiskerprops=dict(lw=0.8),
                            capprops=dict(lw=0.8), boxprops=dict(lw=0.8))
            bp["boxes"][0].set_facecolor(PALETTE_MID[s] if i == 0 else PALETTE_DARK[s])
            bp["boxes"][0].set_alpha(0.45 if i == 0 else 0.85)
            ax.scatter(pos + rng.uniform(-0.17, 0.17, len(d)), d, s=10, color="black", alpha=0.4,
                       linewidth=0, zorder=3)
            ticks.append(pos); labels.append(lab); pos += 1
        centers.append((start + pos - 1) / 2); pos += 0.7
    ax.axhline(THRESHOLD, color="black", ls="--", lw=0.9)
    ax.text(1.85, 2.3, f"{THRESHOLD:.1f} {ANGSTROM}", fontsize=8, ha="center", va="bottom",
            bbox=dict(fc="white", ec="none", pad=0.8))
    ax.set_yscale("log"); ax.set_ylim(0.15, 60); log_axis(ax, "y")
    ax.set_ylabel(f"sPyRMSD ({ANGSTROM}), log scale")
    ax.set_xticks(ticks); ax.set_xticklabels(labels, fontsize=8.2)
    for c, s in zip(centers, STRATEGIES):
        ax.text(c, -0.165, LABELS[s], ha="center", fontsize=9.4, transform=ax.get_xaxis_transform())
    panel_letter(ax, "b", -0.075)
    ax.text(0.5, 1.1, f"n = {int(common.sum())} complexes resolved under all four strategies",
            transform=ax.transAxes, ha="center", fontsize=8, style="italic")
    return fig


# --------------------------------------------------------------------- Figure 5
def _paired(b: Benchmark, pose: int):
    xs, ys, who = [], [], []
    for s in STRATEGIES:
        m = b.resolved(s) & b.resolved(s, "bridge")
        xs.extend(b.spy[s][m, pose]); ys.extend(b.bridge[s][m, pose]); who.extend([s] * int(m.sum()))
    return np.array(xs), np.array(ys), np.array(who)


def fig5(b: Benchmark, pose: int):
    xs, ys, who = _paired(b, pose)
    rho = stats.spearmanr(xs, ys)[0]
    mae = np.mean(np.abs(xs - ys))
    fig, ax = plt.subplots(figsize=(5.1, 5.1))
    for s in STRATEGIES:
        m = who == s
        ax.scatter(xs[m], ys[m], s=24, color=PALETTE_MID[s], edgecolor="black", linewidth=0.45,
                   alpha=0.82, label=LABELS[s], zorder=3)
    lo, hi = 0.15, 60
    ax.plot([lo, hi], [lo, hi], color="black", lw=1.1, zorder=2, alpha=0.75, label="Identity (y = x)")
    ax.axvline(THRESHOLD, color="black", ls="--", lw=0.8, zorder=1)
    ax.axhline(THRESHOLD, color="black", ls="--", lw=0.8, zorder=1)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
    log_axis(ax, "x"); log_axis(ax, "y")
    ax.set_xlabel(f"sPyRMSD ({ANGSTROM})"); ax.set_ylabel(f"Coordinate-bridge RMSD ({ANGSTROM})")
    ax.set_title("Top-1 pose" if pose == 0 else "Best pose", fontsize=11, fontweight="bold", pad=10)
    ax.text(0.04, 0.965, f"n = {len(xs)}\n\u03c1 = {rho:.3f}\nMAE = {mae:.3f} {ANGSTROM}",
            transform=ax.transAxes, ha="left", va="top", fontsize=8, linespacing=1.35,
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.8))
    panel_letter(ax, "ab"[pose], -0.21, 1.065)
    ax.legend(loc="lower right", fontsize=7.2, handlelength=1.4, labelspacing=0.35)
    return fig


# ------------------------------------------------------------------ SI figures
def figS1(b: Benchmark, pose: int):
    th = np.linspace(0, 8, 801)
    fig, ax = plt.subplots(figsize=(4.8, 4.0))
    for s in STRATEGIES:
        a = b.spy[s][b.resolved(s), pose]
        ax.plot(th, [100 * np.sum(a < t) / N_COMPLEXES for t in th], color=PALETTE_MID[s], lw=1.9,
                label=LABELS[s], zorder=3)
        ax.scatter([THRESHOLD], [100 * np.sum(a < THRESHOLD) / N_COMPLEXES], s=26, color=PALETTE_MID[s],
                   edgecolor="black", linewidth=0.7, zorder=4)
    ax.axvline(THRESHOLD, color="black", ls="--", lw=0.9, zorder=1)
    ax.set_xlim(0, 8); ax.set_ylim(0, 100); ax.yaxis.set_major_locator(MultipleLocator(20))
    ax.set_xlabel(f"RMSD threshold ({ANGSTROM})")
    ax.set_ylabel("Complexes with RMSD below threshold\n(% of all 85)")
    ax.set_title("Top-1 pose" if pose == 0 else "Best pose", fontsize=11, fontweight="bold", pad=10)
    panel_letter(ax, "ab"[pose], -0.16, 1.07)
    ax.legend(loc="lower right", fontsize=8, handlelength=1.6, labelspacing=0.4)
    return fig


def _stacked(ax, segments, colors, hatches, text_colors, bar_width=0.6):
    bottom = np.zeros(len(STRATEGIES))
    for vals, cols, hatch, tcol in zip(segments, colors, hatches, text_colors):
        for x, (b0, v, c) in enumerate(zip(bottom, vals, cols)):
            h = v * 100 / N_COMPLEXES
            ax.bar(x, h, bottom=b0, color=c, hatch=hatch, edgecolor="#555555" if hatch else "black",
                   linewidth=0.8, width=bar_width)
            if v >= 5:
                ax.text(x, b0 + h / 2, str(v), ha="center", va="center", fontsize=8.2,
                        color=tcol[x] if isinstance(tcol, list) else tcol,
                        bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.9) if hatch else None)
            elif v > 0:
                ax.text(x + bar_width / 2 + 0.05, b0 + h / 2, str(v), ha="left", va="center", fontsize=8)
        bottom += np.array(vals) * 100 / N_COMPLEXES


def _stacked_figure(b, segments, text_colors, legend_labels, figsize, xlim=None):
    fig, ax = plt.subplots(figsize=figsize)
    cols = [[PALETTE_DARK[s] for s in STRATEGIES], [tint(PALETTE_DARK[s], 0.45) for s in STRATEGIES],
            [NEUTRAL] * 4, ["white"] * 4]
    _stacked(ax, segments, cols, [None, None, None, "////"], text_colors)
    ax.set_xticks(range(4)); ax.set_xticklabels([LABELS[s] for s in STRATEGIES])
    ax.set_ylim(0, 100); ax.set_yticks([0, 20, 40, 60, 80, 100]); ax.set_ylabel("% of all 85 complexes")
    if xlim:
        ax.set_xlim(*xlim)
    swatches = [Patch(fc="#555555", ec="black", lw=.8, label=legend_labels[0]),
                Patch(fc="#BDBDBD", ec="black", lw=.8, label=legend_labels[1]),
                Patch(fc=NEUTRAL, ec="black", lw=.8, label=legend_labels[2]),
                Patch(fc="white", ec="#555555", hatch="////", lw=.8, label=legend_labels[3])]
    ax.legend(handles=swatches[::-1], loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=1,
              handlelength=1.4, fontsize=8)
    return fig


def figS2(b: Benchmark):
    seg = [[], [], [], []]
    for s in STRATEGIES:
        res = int(b.resolved(s).sum())
        t1 = int(b.success(s, "spy", 0).sum()); bp = int(b.success(s, "spy", 1).sum())
        for k, v in enumerate([t1, bp - t1, res - bp, N_COMPLEXES - res]):
            seg[k].append(v)
    return _stacked_figure(
        b, seg, [["white"] * 4, ["black"] * 4, "black", "black"],
        [f"Top-1 pose < {THRESHOLD:.1f} {ANGSTROM} (strategy color)",
         "Near-native pose generated, not ranked first (light tint)",
         f"No pose < {THRESHOLD:.1f} {ANGSTROM} generated", "sPyRMSD unresolved"], (4.8, 4.4))


def figS3(b: Benchmark):
    seg = [[], [], [], []]
    for s in STRATEGIES:
        a, c = b.resolved(s), b.resolved(s, "bridge")
        for k, v in enumerate([int((a & c).sum()), int((a & ~c).sum()), int((~a & c).sum()), int((~a & ~c).sum())]):
            seg[k].append(v)
    return _stacked_figure(
        b, seg, [["white"] * 4, "black", "black", "black"],
        ["Both routes numerical (strategy color)", "sPyRMSD only (light tint)",
         "Coordinate bridge only", "Neither route"], (4.8, 4.2), xlim=(-0.6, 3.75))


def _violins(ax, data, positions, colors, width=0.6):
    parts = ax.violinplot(data, positions=positions, widths=width, showmeans=False,
                          showmedians=False, showextrema=False)
    for body, c in zip(parts["bodies"], colors):
        body.set_facecolor(c); body.set_edgecolor("black"); body.set_alpha(0.6); body.set_linewidth(0.7)
    for p, d in zip(positions, data):
        q1, med, q3 = np.percentile(d, [25, 50, 75])
        ax.plot([p, p], [np.min(d), np.max(d)], color="black", lw=0.9, zorder=3)
        ax.plot([p, p], [q1, q3], color="black", lw=3.2, zorder=4, solid_capstyle="butt")
        ax.scatter([p], [med], color="white", edgecolor="black", s=15, zorder=5, linewidth=0.9)


def _score_axis(ax):
    ticks = [-14, -12, -10, -8, -6, -4]
    ax.set_yticks(ticks); ax.set_yticklabels([MINUS(str(t)) for t in ticks])
    ax.set_ylabel("Top-1 Vina score (kcal/mol)")


def figS4a(b: Benchmark):
    holm = pairwise(b.score, lambda x, y: stats.wilcoxon(x, y).pvalue)
    fig, ax = plt.subplots(figsize=(4.8, 4.2))
    _violins(ax, [b.score[s] for s in STRATEGIES], range(4), [PALETTE_DARK[s] for s in STRATEGIES])
    for k, s in enumerate(STRATEGIES):
        ax.text(k, -15.2, MINUS(f"mean {b.score[s].mean():.2f}"), ha="center", fontsize=8)
    for other, y, a in (("P2Rank", -2.6, 2), ("F-pocket", -1.7, 1), ("Co-crystallized", -0.8, 0)):
        p = holm[f"{other}|Blind"]
        ax.plot([a, a, 3, 3], [y - .25, y, y, y - .25], color="black", lw=0.8)
        ax.text((a + 3) / 2, y + .15, "P < 0.001" if p < 1e-3 else f"P = {p:.3f}", ha="center", va="bottom", fontsize=8)
    ax.set_ylim(-16.2, 0.3); _score_axis(ax)
    ax.set_xticks(range(4)); ax.set_xticklabels([LABELS[s] for s in STRATEGIES])
    ax.set_title("Docking score by strategy", fontsize=11, fontweight="bold", pad=10)
    panel_letter(ax, "a", -0.16, 1.07)
    return fig


def figS4b(b: Benchmark):
    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    pos = 0.0
    ticks = []
    for s in STRATEGIES:
        ok = b.resolved(s); rec = ok & b.success(s)
        _violins(ax, [b.score[s][rec], b.score[s][ok & ~rec]], [pos - 0.2, pos + 0.2],
                 [PALETTE_DARK[s], "#B8B8B8"], width=0.34)
        r = score_auc(b, s)
        ax.text(pos, 0.975, f"AUC {r['auc']:.2f}\n[{r['lo']:.2f}, {r['hi']:.2f}]",
                transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=7.8, linespacing=1.2)
        ticks.append(pos); pos += 1.3
    ax.set_xticks(ticks); ax.set_xticklabels([LABELS[s] for s in STRATEGIES])
    ax.set_ylim(-14.8, 0.3); _score_axis(ax)
    panel_letter(ax, "b", -0.16, 1.07)
    ax.set_title("Docking score and pose recovery", fontsize=11, fontweight="bold", pad=10)
    ax.legend(handles=[Patch(fc="#555555", alpha=.6, ec="black", label=f"Recovered, RMSD < {THRESHOLD:.1f} {ANGSTROM} (strategy color)"),
                       Patch(fc="#B8B8B8", alpha=.6, ec="black", label="Not recovered")],
              loc="lower right", fontsize=7.6, frameon=True, facecolor="white", edgecolor="none", framealpha=0.95)
    return fig


FIGURES = {
    "4a": ("Main", "Figure4a_success_rate", lambda b: fig4a(b)),
    "4b": ("Main", "Figure4b_rmsd_distribution", lambda b: fig4b(b)),
    "5a": ("Main", "Figure5a_agreement_top1", lambda b: fig5(b, 0)),
    "5b": ("Main", "Figure5b_agreement_bestpose", lambda b: fig5(b, 1)),
    "S1a": ("SI", "FigureS1a_ecdf_top1", lambda b: figS1(b, 0)),
    "S1b": ("SI", "FigureS1b_ecdf_bestpose", lambda b: figS1(b, 1)),
    "S2": ("SI", "FigureS2_failure_modes", lambda b: figS2(b)),
    "S3": ("SI", "FigureS3_correspondence_coverage", lambda b: figS3(b)),
    "S4a": ("SI", "FigureS4a_score_by_strategy", lambda b: figS4a(b)),
    "S4b": ("SI", "FigureS4b_score_vs_outcome", lambda b: figS4b(b)),
}


def main(argv=None) -> int:
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data", type=Path, default=here / "data" / "astex_benchmark_S2.csv")
    ap.add_argument("--outdir", type=Path, default=Path("figures"))
    ap.add_argument("--formats", nargs="+", default=["svg", "png"], choices=["svg", "png"])
    ap.add_argument("--dpi", type=int, default=1000, help="PNG resolution (default 1000)")
    ap.add_argument("--only", nargs="+", choices=sorted(FIGURES), help="subset of figures")
    args = ap.parse_args(argv)
    bench = load_benchmark(args.data)
    family = setup_style()
    print(f"font: {family}; data: {args.data}")
    for key in (args.only or FIGURES):
        sub, name, build = FIGURES[key]
        save(build(bench), args.outdir / sub, name, args.formats, args.dpi, family)
        print(f"  wrote {sub}/{name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
