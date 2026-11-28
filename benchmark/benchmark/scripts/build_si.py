import sys
S="PATH/TO/WORKING_DIR/"  # set to your working directory (see benchmark/README.md)
t=open(S+"fp/latex/IBDock_JCIM_SI_v3.tex").read()
f=lambda n:open(S+f"final/{n}.tex").read().rstrip()

a=t.index("\\section*{Table S8."); b=t.index("\\section*{Table S10. Cross-docking}")
s89=r"""\section*{Table S8. Box geometry and site localization per complex}

Box volume, offset between the box center and the native-ligand centroid, and percentage of native-ligand heavy atoms inside the box, computed from the stored \texttt{grid\_config.txt} files and the crystal ligand. Boxes were identical across the nine exhaustiveness--seed runs. Co-crystallized boxes are centered on the native ligand (offset 0 \AA, 100\% inside) and are not repeated. $^{\dagger}$ the box does not enclose every ligand heavy atom. Edge lengths and per-run values are in the accompanying CSV file.

{\scriptsize
\begin{longtable}{@{}lrrrrrrrrrr@{}}
\toprule
 & \multicolumn{4}{c}{Volume (\AA$^3$)} & \multicolumn{3}{c}{Offset (\AA)} & \multicolumn{3}{c}{Inside (\%)} \\
\cmidrule(lr){2-5}\cmidrule(lr){6-8}\cmidrule(lr){9-11}
PDB & Co-cryst. & fpocket & P2Rank & Blind & fpocket & P2Rank & Blind & fpocket & P2Rank & Blind \\
\midrule
\endfirsthead
\toprule
PDB & Co-cryst. & fpocket & P2Rank & Blind & fpocket & P2Rank & Blind & fpocket & P2Rank & Blind \\
\midrule
\endhead
\bottomrule
\endlastfoot
"""+open(S+"final/si_s8.tex").read()+r"""
\end{longtable}}

\section*{Table S9. Naive assignment versus symmetry-aware correspondence}

The naive assignment pairs each docked heavy atom with a distinct reference heavy atom of the same element by minimizing the total squared distance (optimal one-to-one assignment). Restricted to poses with a numerical sPyRMSD. Because any graph-consistent mapping is also a valid assignment, naive RMSD $\leq$ sPyRMSD, and every discordant call is a success under the naive assignment that sPyRMSD did not confirm.

\subsection*{Table S9a. Discordant calls per run}
{\footnotesize
\begin{longtable}{@{}rrlrrrr@{}}
\toprule
Exhaustiveness & Seed & Strategy & Poses & Naive-only Top-1 & Naive-only best pose & sPyRMSD-only (either) \\
\midrule
\endfirsthead
\toprule
Exhaustiveness & Seed & Strategy & Poses & Naive-only Top-1 & Naive-only best pose & sPyRMSD-only (either) \\
\midrule
\endhead
\bottomrule
\endlastfoot
"""+open(S+"final/si_s9a.tex").read()+r"""
\end{longtable}}

\subsection*{Table S9b. Poses with a discordant Top-1 or best-pose call}
{\footnotesize
\begin{longtable}{@{}lrrlrrrrl@{}}
\toprule
PDB & Ex. & Seed & Strategy & sPyRMSD Top-1 & Naive Top-1 & sPyRMSD best & Naive best & Discordant \\
\midrule
\endfirsthead
\toprule
PDB & Ex. & Seed & Strategy & sPyRMSD Top-1 & Naive Top-1 & sPyRMSD best & Naive best & Discordant \\
\midrule
\endhead
\bottomrule
\endlastfoot
"""+open(S+"final/si_s9b.tex").read()+r"""
\end{longtable}}
\noindent{\footnotesize RMSD in \AA. T: Top-1 call discordant; B: best-pose call discordant. Per-pose values for all 3{,}060 runs are in the accompanying CSV file.}

"""
t=t[:a]+s89+t[b:]
i=t.index("\\section*{Table S10. Cross-docking}"); j=t.index("\\begin{figure*}[h]",i)
new=r"""\section*{Table S10. Cross-docking}

All values are for exhaustiveness 16 and seed 42 (the single cross-docking run). Success is sPyRMSD $<$ 2.0~\AA\ against the ligand crystal pose superposed into the receptor frame. Non-cognate pairs in the primary analysis are summarized in the cross-docking table of the main text; the tables below give the breakdowns. ``Site'' is the number or percentage of Top-1 poses with centroid within 4~\AA\ of the reference centroid (descriptive, added after the first target was scored).

\subsection*{Table S10a. Success by pocket-difficulty tier}
\begin{table}[H]
\centering\footnotesize
\begin{tabular}{lllrrrr}
\toprule
Target & Strategy & Tier & $n$ & Top-1 & Best-pose & Site \\
\midrule
"""+f("si_tier")+r"""
\bottomrule
\end{tabular}
\par\smallskip\raggedright\footnotesize Tiers by 8~\AA\ pocket C$\alpha$ RMSD between receptor and the ligand's own receptor: easy $<$ 1.5, medium 1.5--2.5, hard $>$ 2.5~\AA. Cognate rows are counts, non-cognate rows count pairs (primary analysis).
\end{table}

\subsection*{Table S10b. Flagged non-equivalent-site pairs (centroid offset $>$ 6~\AA)}
\begin{table}[H]
\centering\footnotesize
\begin{tabular}{llrrrr}
\toprule
Target & Strategy & Pairs & Top-1 success & Best-pose success & Top-1 at site \\
\midrule
"""+f("si_flag")+r"""
\bottomrule
\end{tabular}
\par\smallskip\raggedright\footnotesize Counts of pairs. ESR1 has no flagged pairs.
\end{table}

\subsection*{Table S10c. Top-1 success conditional on the site being reached}
\begin{table}[H]
\centering\footnotesize
\begin{tabular}{llrrrr}
\toprule
Target & Strategy & Pairs & Top-1 at site (count) & Top-1 success given site (\%) & Best-pose minus Top-1 (points) \\
\midrule
"""+f("si_cond")+r"""
\bottomrule
\end{tabular}
\end{table}

\subsection*{Table S10d. Reference-pose clash and Top-1 success (co-crystallized-ligand box, non-cognate pairs)}
\begin{table}[H]
\centering\footnotesize
\begin{tabular}{llrrrr}
\toprule
Target & Reference pose & $n$ & Top-1 (\%) & Best-pose (\%) & Site (\%) \\
\midrule
"""+f("si_clash")+r"""
\bottomrule
\end{tabular}
\par\smallskip\raggedright\footnotesize Counts, with percentages in parentheses. Clash: any reference heavy atom within 2.5~\AA\ of a receptor heavy atom (descriptive, added post hoc).
\end{table}

\subsection*{Table S10e. Paired differences between strategies}
\begin{table}[H]
\centering\footnotesize
\begin{tabular}{llcc}
\toprule
Target & Comparison & $\Delta$ Top-1 (points, 95\% CI) & $\Delta$ Site capture (points, 95\% CI) \\
\midrule
"""+f("si_paired")+r"""
\bottomrule
\end{tabular}
\par\smallskip\raggedright\footnotesize Non-cognate pairs, primary analysis. Intervals: two-way cluster bootstrap with receptor and ligand indices shared across strategies (10{,}000 resamples, fixed seed). Pooled rows are descriptive.
\end{table}

\subsection*{Table S10f. Per-pair sPyRMSD (\AA), all 1{,}320 docked pairs}
{\scriptsize
\begin{longtable}{@{}lllrrrrrrrrr@{}}
\toprule
Receptor/ligand & Tier & Pocket RMSD & F & \multicolumn{2}{c}{Co-cryst.} & \multicolumn{2}{c}{fpocket} & \multicolumn{2}{c}{P2Rank} & \multicolumn{2}{c}{Blind} \\
\cmidrule(lr){5-6}\cmidrule(lr){7-8}\cmidrule(lr){9-10}\cmidrule(lr){11-12}
& & & & Top-1 & Best & Top-1 & Best & Top-1 & Best & Top-1 & Best \\
\midrule
\endfirsthead
\toprule
Receptor/ligand & Tier & Pocket RMSD & F & Top-1 & Best & Top-1 & Best & Top-1 & Best & Top-1 & Best \\
\midrule
\endhead
\bottomrule
\endlastfoot
"""+f("si_pairs")+r"""
\end{longtable}}
\noindent{\footnotesize F = flagged non-equivalent-site pair (analysed separately); cog = cognate pair; -- = not scored. The covalent ESR1 ligand 6CHW was not docked. The pairs of two PPARG ligands (1ZGY and 5Y2T) were scored with a geometry-perceived reference graph (Methods).}

\section*{Table S11. Repurposing case study}

All four strategies, exhaustiveness 16, seed 42; imatinib (37 heavy atoms) docked into each receptor. sPyRMSD against crystal imatinib (STI) is given for 1IEP and 1T46 only; ``centroid'' is the distance between the centroids of the Top-1 pose and the crystal ligand.

\begin{table}[H]
\centering\footnotesize
\begin{tabular}{llllrrrrr}
\toprule
Receptor & Strategy & Box source; size (\AA) & Score & Top-1 RMSD & Best RMSD & Centroid & Modes \\
\midrule
ABL1 (1IEP) & Co-crystallized & STI A:201; 20 $\times$ 28 $\times$ 24 & $-$12.31 & 0.78 & 0.78 & 0.3 & 9 \\
 & fpocket & 32 $\times$ 42 $\times$ 30 & $-$10.58 & 13.15 & 11.78 & 3.9 & 9 \\
 & P2Rank & 22 $\times$ 22 $\times$ 22 & $-$10.65 & 12.40 & 11.66 & 11.1 & 9 \\
 & Blind & 68 $\times$ 60 $\times$ 80 & $-$12.32 & 0.80 & 0.80 & 0.3 & 7 \\
KIT (1T46) & Co-crystallized & STI A:3; 20 $\times$ 18 $\times$ 28 & $-$12.61 & 0.78 & 0.78 & 0.2 & 8 \\
 & fpocket & 34 $\times$ 24 $\times$ 44 & $-$11.63 & 12.75 & 6.75 & 3.6 & 9 \\
 & P2Rank & 22 $\times$ 22 $\times$ 22 & $-$11.43 & 12.82 & 10.28 & 3.5 & 9 \\
 & Blind & 72 $\times$ 70 $\times$ 70 & $-$9.74 & 11.83 & 11.72 & 10.7 & 9 \\
SRC (2SRC) & Co-crystallized & ANP A:1; 24 $\times$ 20 $\times$ 16 & $-$9.53 & n/a & n/a & n/a & 9 \\
 & fpocket & 36 $\times$ 32 $\times$ 26 & $-$10.30 & n/a & n/a & n/a & 9 \\
 & P2Rank & 22 $\times$ 22 $\times$ 22 & $-$10.04 & n/a & n/a & n/a & 9 \\
 & Blind & 80 $\times$ 86 $\times$ 78 & $-$9.83 & n/a & n/a & n/a & 9 \\
EGFR (1M17) & Co-crystallized & AQ4 A:999; 26 $\times$ 16 $\times$ 18 & $-$9.62 & n/a & n/a & n/a & 9 \\
 & fpocket & 22 $\times$ 22 $\times$ 24 & $-$9.60 & n/a & n/a & n/a & 9 \\
 & P2Rank & 22 $\times$ 22 $\times$ 22 & $-$9.58 & n/a & n/a & n/a & 9 \\
 & Blind & 110 $\times$ 82 $\times$ 68 & $-$9.61 & n/a & n/a & n/a & 9 \\
\bottomrule
\end{tabular}
\par\smallskip\raggedright\footnotesize Scores in kcal~mol$^{-1}$; RMSD and centroid distance in \AA. The pdbqt files retain poses within 3~kcal~mol$^{-1}$ of the best, so the number of modes varies. Box volumes differ among receptors and strategies, so scores are not comparable across proteins.
\end{table}

\subsection*{Table S11b. Imatinib activity record}
\begin{table}[H]
\centering\footnotesize
\begin{tabular}{llp{2.6cm}p{5.2cm}}
\toprule
Kinase (PDB) & IC$_{50}$ & Assay & Source (ChEMBL document)\\
\midrule
ABL1 (1IEP) & 38~nM & biochemical, v-Abl & Zimmermann \textit{et al.} 1997, \textit{Bioorg. Med. Chem. Lett.} (CHEMBL1130562); same value in a 2005 \textit{J. Med. Chem.} report (CHEMBL1140446) \\
KIT (1T46) & 100~nM & biochemical, recombinant wild-type c-Kit & 2005 \textit{J. Med. Chem.} (CHEMBL1141098; CHEMBL1140446) \\
SRC (2SRC) & $>$100~$\mu$M & biochemical & CHEMBL1130562; CHEMBL1140446 (all seven ChEMBL SRC IC$_{50}$ records are ``$>$'' relations) \\
EGFR (1M17) & $>$100~$\mu$M & biochemical & CHEMBL1130562 (an outlier of 0.11~nM, CHEMBL5214926, was not used) \\
\bottomrule
\end{tabular}
\par\smallskip\raggedright\footnotesize Compound CHEMBL941. Values read from a ChEMBL activity export; spread of wild-type biochemical IC$_{50}$: ABL1 19.9--7{,}700~nM ($n = 42$, median 228~nM), KIT 16--10{,}000~nM ($n = 30$, median 290~nM). Cell-based KIT values were 46--97~nM.
\end{table}

\section*{Table S12. Cross-docking ligand preparation record}

\begin{longtable}{@{}llllrl@{}}
\toprule
Target & Ligand & Preparation & Formal charge & Heavy atoms & Note \\
\midrule
\endhead
"""+f("si_lig")+r"""
\bottomrule
\end{longtable}
\noindent{\footnotesize ``Repaired'': bond orders re-perceived from the template; ``rebuilt'': ligand rebuilt from geometry. Formal charge is shown where it was inferred during repair. Ten ligands were repaired: ESR1 5AAU and 6CHW; MAPK1 4QP3, 4QP4, 4QTA, 5AX3, 5BUJ; PPARG 3R8A, 5Y2T, 5Z5S.}

"""
t=t[:i]+new+t[j:]
# cond table header: simplify (avoid unexplained column)

import re as _re
t=t.replace("\\TODO{state the exhaustiveness and seed used for this panel (expected: 16 and 42); if it was run once, say so.}","All 28 conditions were run once, at exhaustiveness 16 and seed 42 (confirmed from the Vina logs).")
t=_re.sub(r"\\TODO\{record the failure mechanism, or state that none was recorded\.\}","No failure mechanism was recorded for this case.",t)
open(S+"v7/IBDock_JCIM_SI_v4.tex","w").write(t); print("ok")
