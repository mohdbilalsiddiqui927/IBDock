import sys,pandas as pd,numpy as np,pickle
S="/tmp/claude-0/-home-claude/fa06f76e-b715-54d8-a3bb-7849315bf476/scratchpad/"
p=S+"v7/IBDock_JCIM_v7.tex"; t=open(p).read()
def rep(old,new):
    global t
    if t.count(old)!=1: sys.exit(f"ANCHOR {t.count(old)}: {old[:90]!r}")
    t=t.replace(old,new)
def span(a,b,new,incl=True):
    global t
    i=t.index(a); j=t.index(b,i)+(len(b) if incl else 0); t=t[:i]+new+t[j:]
# ---- tables from data
M=pd.read_csv(S+"final/astex_merged.csv"); R=M[M.sp_top1.notna()]
def row(lab,g):
    n=len(g); nv=int((g.nn_top1<2).sum()); sp=int((g.sp_top1<2).sum()); a=int(((g.nn_top1<2)&~(g.sp_top1<2)).sum()); b=int((~(g.nn_top1<2)&(g.sp_top1<2)).sum())
    return f"    {lab} & {n} & {nv} ({100*nv/n:.1f}) & {sp} ({100*sp/n:.1f}) & {a} & {b} \\\\"
NM={"cocryst":"Co-crystallized","fpocket":"fpocket","p2rank":"P2Rank","blind":"Blind"}
rows=[row(NM[s],R[(R.ex==16)&(R.strategy==s)]) for s in NM]+["    \\midrule",row("All four strategies (exhaustiveness 16)",R[R.ex==16]),row("All 36 conditions",R)]
G=pd.read_csv(S+"final/astex_geometry.csv"); G=G[(G.ex==16)&(G.seed==42)]
# ---- methods
span("To separate the two properties experimentally, a\nvolume-matched control was defined:","exact construction and whether runs are complete}.",
"The contribution of box volume was then assessed with a logistic generalized estimating equation (Statistical Analysis). An experimental volume-matched control, in which the co-crystallized-ligand center is combined with the box volume of each other strategy, was not run; this is stated as a limitation.")
span("generalized estimating equation (exchangeable working correlation,","covariates \\todo{finalize covariate set}.",
"generalized estimating equation (exchangeable working correlation,\nclustered on complex) for Top-1 success at exhaustiveness 16 over the three seeds (255 observations per strategy; unresolved sPyRMSD counted as failure, as in the primary endpoint), with strategy, log$_{10}$ box volume, and ligand heavy-atom count as covariates. The model was fitted post hoc, and rotatable-bond count was not included.")
rep("Naive and sPyRMSD success calls were compared\nwith an exact McNemar test on discordant pairs.","Naive and sPyRMSD success calls were compared by counting discordant poses; no test was applied, because poses of the same complex across seeds are not independent.")
span("also scored with a naive heavy-atom mapping \\todo{define: identity","(file-order) mapping and/or spatial nearest-neighbour assignment}.",
"also scored with a naive mapping that ignores the molecular graph: each docked heavy atom was assigned to a distinct reference heavy atom of the same element by minimizing the total squared distance (optimal one-to-one assignment), and RMSD was computed over that assignment. A file-order mapping was not used because atom order differs between the reference and docked files in 80 of 85 complexes and the docked files do not retain atom names. Every graph-consistent mapping is also a valid assignment, so the naive RMSD cannot exceed the sPyRMSD value.")
# ---- results: box geometry
span("sampling per unit volume, so the gradient in Table~\\ref{tbl:astex}","heavy-atom containment per strategy.}",
"""sampling per unit volume, so the gradient in Table~\\ref{tbl:astex}
cannot be attributed to site information alone. Median box volumes were 5{,}832~\\AA$^3$ (co-crystallized ligand), 18{,}720~\\AA$^3$ (fpocket), 10{,}648~\\AA$^3$ (P2Rank, a fixed 22~\\AA\\ cube), and 462{,}080~\\AA$^3$ (blind), the blind box being 79 times the co-crystallized box. Co-crystallized boxes are centered on the native ligand and contain all of its heavy atoms by construction. The fpocket box was centered a median 2.9~\\AA\\ from the native-ligand centroid and enclosed the whole ligand in 78 of 85 complexes, the P2Rank box 4.3~\\AA\\ (IQR 1.9--21.5) and 54 of 85, and the blind box 13.2~\\AA\\ and all 85. Top-1 success depended on enclosure: it was 56.0\\% (fpocket) and 58.0\\% (P2Rank) when the box enclosed the whole ligand, against 0.0\\% in the 7 fpocket and 31 P2Rank complexes where it did not (seed-mean over exhaustiveness 16, resolved and unresolved complexes counted together). The P2Rank deficit relative to the co-crystallized-ligand box (59.6\\% overall) therefore reflects boxes that missed the site, whereas blind docking enclosed the ligand in every complex and still recovered 29.0\\%, which points to sampling over a much larger volume or pose selection rather than to site loss.""")
rep("    Median box volume (\\AA$^3$) [IQR] & \\tbd & \\tbd & \\tbd & \\tbd \\\\\n    Median centroid offset (\\AA) [IQR] & \\tbd & \\tbd & \\tbd & \\tbd \\\\\n    Ligand atoms inside box, median (\\%) & \\tbd & \\tbd & \\tbd & \\tbd \\\\\n    Complexes with ligand fully inside box & \\tbd & \\tbd & \\tbd & \\tbd \\\\",
"    Median box volume (10$^3$~\\AA$^3$) [IQR] & 5.8 [4.9--7.2] & 18.7 [14.9--24.6] & 10.6 [10.6--10.6] & 462 [309--690] \\\\\n    Median centroid offset (\\AA) [IQR] & 0.0 [0.0--0.0] & 2.9 [1.4--4.3] & 4.3 [1.9--21.5] & 13.2 [9.1--19.3] \\\\\n    Ligand atoms inside box, median (\\%) & 100 & 100 & 100 & 100 \\\\\n    Complexes with ligand fully inside box & 85/85 & 78/85 & 54/85 & 85/85 \\\\")
rep("geometry rows are to be computed from the stored\n  \\fn{grid\\_config.txt} files.","geometry rows were computed from the stored\n  \\fn{grid\\_config.txt} files and the crystal ligand (per-complex values in Table~S8). Boxes were identical across seeds and exhaustiveness. Containment counts a ligand as inside only if every heavy atom lies within the box.")
# GEE table
res=pickle.load(open(S+"final/gee_res.pkl","rb"))
f=lambda v:f"{v[0]:.2f} ({v[1]:.2f}--{v[2]:.2f})"
grows="\n".join(f"    {NM[s]} & {f(res['unadj'][s])} & {f(res['adj_vol'][s])} & {f(res['adj_vol_ha'][s])} \\\\" for s in ["fpocket","p2rank","blind"])
span("    Strategy (vs co-crystallized) & Unadjusted OR (95\\% CI) & Adjusted OR (95\\% CI) & Volume-matched Top-1 \\\\","    Blind   & \\tbd & \\tbd & \\tbd \\\\",
"    Strategy (vs co-crystallized) & Unadjusted OR (95\\% CI) & Adjusted for volume & Adjusted for volume and ligand size \\\\\n    \\midrule\n"+grows)
span("  \\textsuperscript{\\emph{a}} Logistic GEE clustered on complex;","co-crystallized-ligand center combined with that strategy's box volume.",
"  \\textsuperscript{\\emph{a}} Logistic GEE (exchangeable correlation) clustered on complex; Top-1 success (sPyRMSD $<$ 2.0~\\AA, unresolved = failure) at exhaustiveness 16, three seeds, 1{,}020 observations. Odds ratios below 1 mean lower odds than the co-crystallized-ligand strategy. In the volume-adjusted model the odds ratio per tenfold increase in box volume was 0.41 (0.18--0.92). Log volume and the blind-docking indicator were correlated at $r = 0.95$. Post hoc analysis.")
span("\\todo{Interpretation to be written once Table~\\ref{tbl:gee} is filled:","survives.}",
"""The unadjusted odds ratios reproduce the gradient. After adjustment for box volume the fpocket deficit disappeared (odds ratio 1.13) and the blind deficit was no longer distinguishable from no effect, but the interval for blind docking (0.34--6.66) is too wide to exclude either no effect or a large one, because blind docking accounts for nearly all of the variation in box volume ($r = 0.95$). For blind docking the data therefore cannot separate site information from box size. The P2Rank deficit persisted (odds ratio 0.50, 0.33--0.75) at a nearly constant volume, in line with the missed boxes described above. We therefore describe the effect as that of the strategy as implemented, a bundle of site placement and box size, and do not attribute the blind gap to site information alone. An experimental volume-matched control was not run and would be needed to apportion it.""")
rep("the volume-matched control\n(Table~\\ref{tbl:gee}) remains necessary.","the experimental volume-matched control, which was not run, remains necessary (Table~\\ref{tbl:gee} gives the model-based adjustment).")
# naive results
span("Table~\\ref{tbl:naive} compares success calls under the naive mapping and","for symmetric\nligands.}",
"""Table~\\ref{tbl:naive} compares success calls under the naive assignment and under sPyRMSD, pose by pose, for the complexes with a numerical sPyRMSD (70 of 85). Because any graph-consistent mapping is also a valid assignment, the naive RMSD cannot exceed the sPyRMSD value, and the discordant calls were all of one kind: successes under the naive assignment that sPyRMSD did not confirm. At exhaustiveness 16 this applied to 28 of 839 Top-1 poses (3.3\\%; 9, 9, 0, and 10 for the co-crystallized-ligand, fpocket, P2Rank, and blind strategies), and over all 36 conditions to 86 of 2{,}518 Top-1 poses (3.4\\%) and 30 of 2{,}518 best poses (1.2\\%). No pose was a success under sPyRMSD and a failure under the naive assignment. Six complexes accounted for all 28 discordant Top-1 calls at exhaustiveness 16 (1Q41 and 1N1M, 9 each; 1XM6, 6; 1HNN, 2; 1LPZ and 1OPK, 1 each). The effect on any single success rate is therefore small (at most 4.8\\% of resolved poses, for blind docking) but systematic in direction and concentrated in a few ligands. The naive assignment also returns a value for the 15 complexes for which sPyRMSD returns none (128 of 181 such Top-1 poses at exhaustiveness 16 would be called successes), but no graph-based check supports those mappings, which is why they are kept as unresolved.""")
span("  \\caption{Success Calls under Naive Mapping versus sPyRMSD (Top-1)$^{a}$}","  \\textsuperscript{\\emph{a}} Restricted to complexes with a numerical\n  sPyRMSD; \\todo{state denominator}.",
"""  \\caption{Success Calls under Naive Assignment versus sPyRMSD (Top-1)$^{a}$}
  \\label{tbl:naive}
  \\centering
  \\small
  \\begin{tabular}{lccccc}
    \\toprule
    Strategy & Poses & Naive success & sPyRMSD success & Naive only & sPyRMSD only \\\\
    \\midrule
"""+"\n".join(rows)+"""
    \\bottomrule
  \\end{tabular}

  \\raggedright\\footnotesize
  \\textsuperscript{\\emph{a}} Poses are the Top-1 poses of the three seeds pooled (exhaustiveness 16) for the 70 complexes with a numerical sPyRMSD (1R58 under blind docking, seed 123, is unresolved, so $n = 209$). Success counts with percentages in parentheses are of poses, and so differ from the all-85 rates in Table~\\ref{tbl:astex}. Naive only: success under the naive assignment and failure under sPyRMSD; sPyRMSD only: the reverse. Poses of one complex across seeds are not independent, so no test is given.""")
# abstract / conclusion / limitations
rep("Symmetry-aware correspondence changed \\tbd\\ of \\tbd\\ success calls relative to a naive atom mapping.","A naive element-matched atom assignment called 86 of 2{,}518 Top-1 poses (3.4\\%) successes that sPyRMSD did not confirm, and never the reverse. Strategy and box volume are bundled: adjusting for volume removed the fpocket deficit but could not separate the blind-docking deficit from its 79-fold larger box, and the P2Rank deficit mainly reflected boxes that missed the ligand.")
rep("because strategy\nalso sets box volume, \\todo{state how much of this survives volume\nadjustment and the volume-matched control}.","because strategy\nalso sets box volume, a model adjusted for volume removed the fpocket deficit, could not resolve the blind deficit (blind docking carries nearly all the variation in volume), and left the P2Rank deficit, which reflects boxes that missed the ligand in 31 of 85 complexes.")
rep("and \\tbd\\ when it was not.","and without it only for ABL1 under blind docking, not for KIT under any strategy.")
rep("In\naddition, strategy and box volume co-vary by construction, and 15/85","In\naddition, strategy and box volume co-vary by construction, no experimental volume-matched control was run, the volume-adjusted model cannot separate site information from box size for blind docking, and 15/85")
rep("  geometry; naive versus symmetry-aware per-pose calls.","  geometry (per complex); naive versus symmetry-aware calls (summary and discordant poses).")
rep("Tables S8--S9 \\todo{fill from the Astex per-pose data}: box","Tables S8--S9: box")

span("    Heavy-atom count mismatch & \\tbd & \\tbd \\\\","    Other & \\tbd & \\tbd \\\\",
"""    Heavy-atom count mismatch & 0 & none (docked and reference counts equal in all 3{,}060 runs) \\\\
    Sanitization / parsing failure & 0 & none recorded \\\\
    Graph automorphism / bond-perception mismatch & 11 & 1J3J, 1K3U, 1L7F, 1OYT, 1T40, 1U4D, 1UML, 1UOU, 1XOQ, 1YGC, 2BSM \\\\
    Covalent or modified ligand & 0 & none recorded \\\\
    Other or not recorded & 4 & 1GPK, 1HP0 (matching-search limit), 1HQ2 (none recorded), 2BR1 (upstream graph construction) \\\\""")
rep("Cause categories are provisional; assign\n  from Table~S1.","Causes are those documented in Table~S1; for most complexes they are inferred and were not logged by the software.")

rep("  \\begin{tabular}{lcc}\n    \\toprule\n    Cause & Complexes & PDB IDs","  \\begin{tabular}{p{4.6cm}cp{7.2cm}}\n    \\toprule\n    Cause & Complexes & PDB IDs")

rep("strategy effect before and after adjustment for box volume, together\nwith the volume-matched control.","strategy effect before and after adjustment for box volume (the\nexperimental volume-matched control was not run).")

rep("otherwise a predefined 12~\\AA\\ base extent is\nsubstituted in each dimension before padding.","otherwise a predefined 12~\\AA\\ base extent is\nsubstituted in each dimension before padding. In the benchmark runs the pocket-point extent was not used: the fallback (a 22~\\AA\\ cube at $p = 5$~\\AA) was applied in 84 of 85 Astex complexes and for all 60 cross-docking receptors, so the P2Rank condition tests the predicted center and not a predicted extent.")
rep("10{,}648~\\AA$^3$ (P2Rank, a fixed 22~\\AA\\ cube)","10{,}648~\\AA$^3$ (P2Rank, the fixed 22~\\AA\\ fallback cube in 84 of 85 complexes)")
rep("Co-crystallized boxes are centered on the native ligand and contain all of its heavy atoms by construction.","Co-crystallized boxes are centered on the native ligand and contain all of its heavy atoms by construction; their dimensions were reproduced exactly by eq~\\ref{eqn:box} ($p = 5$~\\AA) from the crystal-ligand heavy atoms in 83 of 85 complexes, and the boxes of 1JLA and 1MEH were larger than the equation gives for the crystal ligand alone (cause not traced).")
rep("To separate the two properties, Table~\\ref{tbl:gee} reports","""\\begin{figure*}[htbp]
  \\centering
  \\includegraphics[width=\\textwidth]{Figure_box_enclosure}
  \\caption{Box geometry and Top-1 success in the Astex set (exhaustiveness 16). (a) Box volume per complex (points) and median (bar) for each strategy. (b) Top-1 success (seed mean, sPyRMSD $<$ 2.0~\\AA, unresolved counted as failure) split by whether the box fully encloses the crystal ligand; $n$ is the number of complexes. The 7 fpocket and 31 P2Rank boxes that did not enclose the ligand gave 0\\% success. (c) Odds ratios (95\\% CI, logistic GEE, log scale) relative to the co-crystallized-ligand strategy before and after adjustment for log box volume; the adjusted interval for blind docking spans 0.34--6.66.}
  \\label{fgr:box}
\\end{figure*}

To separate the two properties, Table~\\ref{tbl:gee} reports""")
rep("\\subsection{Positioning Relative to Existing Docking Workflows}","""\\begin{figure*}[htbp]
  \\centering
  \\includegraphics[width=0.9\\textwidth]{Figure_repurposing}
  \\caption{Imatinib docked into four kinase structures with each search-space strategy (exhaustiveness 16, seed 42). (a) RMSD of the Top-1 pose (bars) and the best saved pose (diamonds) to crystal imatinib in ABL1 (1IEP) and KIT (1T46), log scale; the dashed line marks 2.0~\\AA. (b) Top-1 Vina score for all four receptors, with the reported biochemical IC$_{50}$; scores are not comparable across receptors.}
  \\label{fgr:rep}
\\end{figure*}

\\subsection{Positioning Relative to Existing Docking Workflows}""")
rep("Cross-docking used a single seed and exhaustiveness 16, so seed\nsensitivity","The P2Rank box size was the fixed 22~\\AA\\ fallback in 84 of 85 Astex complexes and in every cross-docking receptor, so that condition does not test P2Rank's extent prediction. Cross-docking used a single seed and exhaustiveness 16, so seed\nsensitivity")

rep("Top-1 success depended on enclosure:","Top-1 success depended on enclosure (Figure~\\ref{fgr:box}b):")
rep("The unadjusted odds ratios reproduce the gradient.","The unadjusted odds ratios reproduce the gradient (Figure~\\ref{fgr:box}c).")
rep("with all four search-space strategies (Table~\\ref{tbl:repurp}).","with all four search-space strategies (Table~\\ref{tbl:repurp}, Figure~\\ref{fgr:rep}).")

span("endpoint. \\todo{State the seed and exhaustiveness used for this panel","feasible, repeat with the\nthree seeds.}","endpoint. All 28 conditions were run once, at exhaustiveness 16 and seed 42 (confirmed from the Vina logs), and the panel was not repeated across seeds.")
rep("Spearman correlation and absolute differences. \\todo{deposit\n\\fn{grid\\_stats.py} and its inputs with the release}","Spearman correlation and absolute differences.")
rep("\\todo{Add them to the repository and update the Zenodo archive.}","The analysis scripts for the grid statistics (\\fn{grid\\_stats.py}), cross-docking scoring and bootstrap, box-geometry, naive-assignment, and GEE analyses, and the repurposing scoring are part of the release. \\todo{Add the cross-docking, repurposing, and new analysis files to the repository and update the Zenodo archive before submission.}")
open(p,"w").write(t); print("ok")
