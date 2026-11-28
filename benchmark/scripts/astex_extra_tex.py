import json
S="PATH/TO/WORKING_DIR/"  # set to your working directory (see benchmark/README.md)
J=json.load(open(S+"astex_extra.json")); ST=["cocryst","fpocket","p2rank","blind"]; NM=dict(cocryst="Co-crystallized ligand",fpocket="fpocket",p2rank="P2Rank",blind="Blind docking")
def P(p): return "$<$0.001" if p<1e-3 else (f"{p:.3f}" if p>=.01 else f"{p:.4f}")
o=[]
o.append(r"""\clearpage
\section*{Table S13. Astex results at exhaustiveness 8 and 32 (counterpart of main-text Table 2)}
Same definitions as Table 2: success = sPyRMSD $<$ 2.0~\AA; all 85 complexes in the denominator, unresolved correspondence coded as non-success; SD is the sample standard deviation of the three per-seed rates; CI is a 10{,}000-resample percentile bootstrap over complexes (resampling uncertainty, so intervals can differ from Table 2 in the last digit). Per-seed counts are for seeds 1/42/123. Exhaustiveness 16 is in Table 2.
{\footnotesize
\begin{longtable}{@{}rlrrrrrr@{}}
\toprule
Exh. & Search space & \multicolumn{3}{c}{Top-1} & \multicolumn{3}{c}{Best-pose}\\
\cmidrule(lr){3-5}\cmidrule(lr){6-8}
 & & mean $\pm$ SD (\%) & 95\% CI & k/85 & mean $\pm$ SD (\%) & 95\% CI & k/85 \\
\midrule
\endhead
""")
for ex in ("8","32"):
    for s in ST:
        a=J[ex]["top1"]["rows"][s]; b=J[ex]["best"]["rows"][s]
        o.append(f"{ex} & {NM[s]} & {a['mean']:.1f} $\\pm$ {a['sd']:.1f} & {a['lo']:.1f}--{a['hi']:.1f} & {'/'.join(map(str,a['per']))} & {b['mean']:.1f} $\\pm$ {b['sd']:.1f} & {b['lo']:.1f}--{b['hi']:.1f} & {'/'.join(map(str,b['per']))} \\\\")
    o.append("\\midrule" if ex=="8" else "")
o.append("\\bottomrule\n\\end{longtable}}\n")
o.append(r"""
\section*{Table S14. Paired comparisons among strategies at exhaustiveness 8, 16 and 32 (counterpart of main-text Table 3)}
Binary columns: exact two-sided McNemar $P$ on the seed-consensus outcome (success in at least two of three seeds; $n$ = 85), Holm-adjusted across the six comparisons within each endpoint and exhaustiveness. Continuous columns: paired Wilcoxon $P$ on the per-complex median RMSD over seeds for complexes with a numerical sPyRMSD under all four strategies and all three seeds ($n$ in the last row of each block), Holm-adjusted likewise. Cochran's $Q$ and Friedman's $Q$ are the omnibus tests.
{\footnotesize
\begin{longtable}{@{}rlllll@{}}
\toprule
Exh. & Comparison & McNemar Top-1 & McNemar best & Wilcoxon Top-1 & Wilcoxon best \\
\midrule
\endhead
""")
lab=lambda x:{"cocryst":"Co-crystallized","fpocket":"fpocket","p2rank":"P2Rank","blind":"blind"}[x]
for ex in ("8","16","32"):
    t,b=J[ex]["top1"],J[ex]["best"]
    for i,pr in enumerate(t["pairs"]):
        a,c=pr.split("|")
        o.append(f"{ex if i==0 else ''} & {lab(a)} vs {lab(c)} & {P(t['mcn'][i])} & {P(b['mcn'][i])} & {P(t['wil'][i])} & {P(b['wil'][i])} \\\\")
    o.append(f" & \\emph{{Omnibus}} & $Q$ = {t['Q']:.2f} ({P(t['pQ'])}) & $Q$ = {b['Q']:.2f} ({P(b['pQ'])}) & Friedman {t['fr']:.2f} ({P(t['pfr'])}), $n$ = {t['n']} & Friedman {b['fr']:.2f} ({P(b['pfr'])}), $n$ = {b['n']} \\\\\n\\midrule")
o.append("\\end{longtable}}\n")
o.append(r"""
\section*{Table S15. Box enclosure and volume adjustment at each exhaustiveness (counterpart of main-text Tables 4 and 5)}
Boxes were identical across seeds and exhaustiveness, so enclosure is the same at every level; success changes. Top-1 success (seed mean, unresolved counted as failure) is split by whether the box fully encloses the crystal ligand (every heavy atom inside). Odds ratios are from the logistic GEE of Table 5 (exchangeable, clustered on complex, three seeds, post hoc), versus the co-crystallized-ligand strategy, unadjusted and adjusted for log$_{10}$ box volume. For blind docking the adjusted interval is wide because volume and the blind indicator are correlated ($r$ = 0.95).
{\footnotesize
\begin{longtable}{@{}rlrrrrll@{}}
\toprule
Exh. & Strategy & Fully inside & Top-1 inside (\%) & Not inside & Top-1 not inside (\%) & OR unadjusted & OR volume-adjusted \\
\midrule
\endhead
""")
for ex in ("8","16","32"):
    for i,s in enumerate(ST):
        e=J[ex]["enc"][s]; g=J[ex]["gee"]
        ou="1 (ref.)" if s=="cocryst" else "{:.2f} ({:.2f}--{:.2f})".format(*g["unadj"][s]); oa="1 (ref.)" if s=="cocryst" else "{:.2f} ({:.2f}--{:.2f})".format(*g["adjvol"][s])
        yo="--" if e["y_out"] is None else f"{e['y_out']:.1f}"
        o.append(f"{ex if i==0 else ''} & {NM[s]} & {e['n_in']} & {e['y_in']:.1f} & {e['n_out']} & {yo} & {ou} & {oa} \\\\")
    o.append("\\midrule")
o.append("\\bottomrule\n\\end{longtable}}\n")
open(S+"si_extra.tex","w").write("\n".join(o))
for ex in ("8","32"):
    print(ex,{s:round(J[ex]["top1"]["rows"][s]["mean"],1) for s in ST},[round(x,4) for x in J[ex]["top1"]["mcn"]],J[ex]["enc"]["p2rank"],J[ex]["gee"]["adjvol"]["blind"])
