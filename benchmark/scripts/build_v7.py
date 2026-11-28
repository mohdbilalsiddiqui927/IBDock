import re, shutil, sys
S="PATH/TO/WORKING_DIR/"  # set to your working directory (see benchmark/README.md)
t=open(S+"fp/latex/IBDock_JCIM_v6.tex").read()
def rep(old,new,cnt=1):
    global t
    n=t.count(old)
    if n!=cnt: sys.exit(f"ANCHOR COUNT {n}!={cnt}: {old[:80]!r}")
    t=t.replace(old,new)
def span(start,end,new,incl_end=False):
    global t
    i=t.index(start); j=t.index(end,i)+(len(end) if incl_end else 0); t=t[:i]+new+t[j:]
tab=lambda n:open(S+f"final/{n}.tex").read()

# ---------------- ABSTRACT
span("\\begin{abstract}\n","\\end{abstract}", r"""\begin{abstract}
Docking outcomes depend on how the search space is defined and how poses are validated, yet both choices are rarely recorded in a form that can be audited. We quantify how much the search-space strategy alone changes pose recovery, and we present IBDock \textit{Provenire} (hereafter IBDock), an open-source AutoDock Vina workflow that stores search-space provenance as project state and makes atom correspondence explicit during reference-based validation, retaining cases without a complete correspondence as unresolved. On the 85-complex Astex Diverse Set, four search-space strategies evaluated over three random seeds and three exhaustiveness levels gave seed-mean Top-1 success (RMSD $<$ 2.0~\AA; exhaustiveness 16) of 59.6\% (co-crystallized ligand), 51.4\% (fpocket), 36.9\% (P2Rank), and 29.0\% (blind). Between seeds, the outcome of 1--3 of 85 complexes changed for the localized strategies but 17--26 for blind docking, which also gained 11.8 points in Top-1 success from exhaustiveness 8 to 32. Symmetry-aware correspondence changed \tbd\ of \tbd\ success calls relative to a naive atom mapping. In cross-docking of 1{,}320 receptor--ligand pairs across three LIT-PCBA targets (ESR1, MAPK1, PPARG), Top-1 success with the co-crystallized-ligand box fell from 93\%, 67\%, and 67\% for cognate pairs to 48\%, 6\%, and 5\% for non-cognate pairs. Blind docking gave the lowest Top-1 success in every target, and the co-crystallized-ligand and P2Rank boxes placed the Top-1 pose at the reference site significantly more often than blind docking in each target; once the site was reached, Top-1 success was similar across strategies within a target. In a retrospective repurposing case, imatinib docked into the ABL1 and KIT crystal structures with the co-crystallized-ligand box reproduced the crystal pose to 0.78~\AA\ in both; with the fpocket, P2Rank, and blind boxes the ABL1 pose was recovered only by blind docking (0.80~\AA) and the KIT pose by none (best saved pose 6.75~\AA\ or worse). Code and the complete benchmark dataset are released openly.
""")

# ---------------- INTRO
span("Here we present IBDock \\textit{Provenire} (hereafter IBDock), a","%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\\section{IBDock Workflow and Implementation}", r"""The question addressed here is therefore narrower than whether one platform outperforms another: how much does a workflow choice that is rarely recorded change the answer? We treat search-space definition as an explicit methodological variable and measure its effect, together with the effect of the atom correspondence used to score a pose. IBDock \textit{Provenire} (hereafter IBDock) is a reproducibility-oriented AutoDock Vina workflow built so that both variables are stored and auditable: the provenance of the search space is written to the project, and atom correspondence is explicit in reference-based validation. The focus is not a new docking engine or pocket-prediction model. We do not benchmark other platforms head-to-head, because differences in preparation, engine, and sampling would confound the comparison; instead we quantify the effect of the choices that a single-strategy workflow fixes by design. We use the Astex Diverse Set\cite{Hartshorn2007} to ask four questions. (i) How much does the search-space strategy change pose recovery, separating the Vina-ranked pose (Top-1) from the best recovered pose? (ii) Because a strategy fixes both where the box is placed and how large it is, how much of that effect is site localization and how much is box volume at a fixed sampling effort? (iii) How stable are the results to the random seed and to exhaustiveness, and how many success calls change when symmetry-aware correspondence replaces a naive atom mapping? (iv) Does the effect of the strategy persist outside redocking, in cross-docking on three targets and in a retrospective repurposing case, and where does it stop mattering?

""")

# ---------------- CROSS-DOCKING DESIGN
rep("separately (\\tbd\\ pairs; MAPK1 $n = 14$, PPARG $n = 36$, ESR1 none).","separately (50 pairs; MAPK1 $n = 14$, PPARG $n = 36$, ESR1 none).")
rep("""The ESR1 receptor of 5UFX contained a stray group coincident
with its ligand, which was removed.""","""The ESR1 receptor of 5UFX contained a stray group coincident
with its ligand, which was removed. Poses were scored from the heavy
atoms of the docked pdbqt file against the heavy atoms of the reference.
Because pdbqt atom order differs from the ligand file, the pose graph was
rebuilt from interatomic distances (covalent-radius sum plus 0.45~\\AA);
in two pairs the tolerance was varied (0.20--0.55~\\AA) until the pose
graph matched the ligand graph. For two PPARG ligands (1ZGY and 5Y2T) the
reference bond table lacked one ring bond (26 bonds, where the ligand has
27) and the two graphs could not be matched; the same distance rule was
therefore applied to the reference coordinates (30 pairs per strategy).
Redocking of these two ligands with the co-crystallized-ligand box gave
2.26~\\AA\\ (1ZGY) and 0.74~\\AA\\ (5Y2T).""")
rep("""run with the settings of Table~\\ref{tbl:protocol}: exhaustiveness 16 and
seed 42 for the primary run, then seeds 1 and 123 \\todo{(state whether
completed)}. Success was""","""run with the settings of Table~\\ref{tbl:protocol} at exhaustiveness 16 and
a single seed (42), the reference slice of the Astex grid; sensitivity to
exhaustiveness and seed was characterized on the Astex Diverse Set and was
not repeated in cross-docking. Success was""")
rep("""as the control for each target.

\\paragraph{Uncertainty.}""","""as the control for each target. Two descriptive quantities were added after
the first target (MAPK1) had been scored and showed Top-1 success near the
floor, and neither was a pre-specified endpoint or enters a hypothesis
test: site capture, defined as the centroid of the Top-1 pose lying within
4~\\AA\\ of the centroid of the reference pose; and reference-pose clash,
defined as any reference heavy atom lying within 2.5~\\AA\\ of a receptor
heavy atom.

\\paragraph{Uncertainty.}""")
rep("""Pooling across targets is descriptive only.""","""For cognate pairs (15--29 per target) success counts are given instead of
intervals. Pooling across targets is descriptive only and is reported both
pair-weighted and as the unweighted mean of the three targets.""")

# ---------------- REPURPOSING DESIGN
rep("\\todo{citation: Nagar et al.}","\\cite{Nagar2002}")
rep("""It was docked into each receptor with all four search-space
strategies, exhaustiveness 8, 16, and 32, and seeds 1, 42, and 123
(4 receptors $\\times$ 4 strategies $\\times$ 9 settings = 144 dockings).""","""It was docked into each receptor with all four search-space
strategies at exhaustiveness 16 and seed 42 (4 receptors $\\times$ 4
strategies = 16 dockings). For the co-crystallized-ligand strategy the
box was built from the ligand of each structure (imatinib in 1IEP and
1T46, AMP-PNP in 2SRC, erlotinib in 1M17), so the boxes of the two
negative controls differ in size from the imatinib boxes.""")
rep("""Reported imatinib activities (IC$_{50}$) for each kinase were
tabulated with their sources (Table~\\ref{tbl:repurp}; \\todo{verify each
value against the primary paper or ChEMBL and give the assay}).""","""Reported imatinib activities (IC$_{50}$) for each kinase were
taken from ChEMBL (compound CHEMBL941) and are tabulated with their
sources (Table~\\ref{tbl:repurp}).""")

# ---------------- RESULTS: CROSS-DOCKING
cross_new=r"""\subsection{Cross-Docking}

Cross-docking covered 1{,}320 docked pairs (59 cognate, 1{,}261 non-cognate), of which 1{,}270 (59 cognate, 1{,}211 non-cognate) formed the primary analysis after the 50 non-equivalent-site pairs were set aside. Cognate pairs, the redocking control, recovered the native pose (Top-1, co-crystallized-ligand box) in 27 of 29 pairs (ESR1), 10 of 15 (MAPK1), and 10 of 15 (PPARG), 79.7\% pooled (47 of 59), which is at or above the Astex Top-1 rate of 59.6\% (Table~\ref{tbl:astex}). The preparation, box construction, and scoring chain therefore reproduces redocking in this set.

For non-cognate pairs, Top-1 success fell to 47.9\% (95\% CI 35.9--59.7) in ESR1, 5.6\% (0.6--13.3) in MAPK1, and 5.2\% (0.0--12.7) in PPARG, drops of 45, 61, and 61 points relative to cognate pairs (Table~\ref{tbl:cross}, Figure~\ref{fgr:cross}a). Best-pose success was 70.2\%, 15.8\%, and 14.9\%: in ESR1 the native pose was generated for most pairs and ranked first for about two thirds of those, whereas in MAPK1 and PPARG it was generated for fewer than one pair in six. Success declined with pocket difficulty in ESR1 (co-crystallized-ligand box: 51.8\% easy, 47.9\% medium, 38.5\% hard; Table~S10), but the decline was shallow compared with the cognate-to-cross drop, and pocket C$\alpha$ RMSD was weakly related to Top-1 RMSD (Spearman $\rho = 0.13$). In MAPK1 and PPARG, too few pairs were medium or hard to resolve a gradient (hard: 0 of 7 and 0 of 8). In ESR1, ligand identity mattered more than pocket difference: ligand-level Top-1 ranged from 2\% to 85% and fell with ligand size ($\rho = -0.50$ with heavy-atom count and $-0.43$ with rotatable bonds, 29 ligands), and pairs that crossed the agonist--antagonist boundary were the hardest (39.5--42.1\% for the three localized boxes, against 49.0--51.9\% within agonist and 59.2--63.3\% within antagonist structures). No such relation was seen in MAPK1 or PPARG ($|\rho| \leq 0.12$).

\paragraph{Strategy.} In ESR1, where non-cognate success was high enough to resolve differences, the co-crystallized-ligand, fpocket, and P2Rank boxes were indistinguishable in Top-1 (47.9, 47.4, and 49.2\%; paired differences within $\pm$1.8 points, all intervals including zero), and blind docking was 23--25 points lower (24.5\%; difference to the co-crystallized-ligand box 23.4 points, 95\% CI 14.4--33.1). In MAPK1 and PPARG, Top-1 success sat at the floor for every strategy (4.6--5.6\% and 2.9--5.2\%) and no difference could be resolved, so the Top-1 endpoint preserves the redocking ordering only where it has power. A more sensitive quantity is whether the Top-1 pose reached the reference site (Table~\ref{tbl:cross}, Figure~\ref{fgr:cross}b). In all three targets the order was co-crystallized ligand $\geq$ P2Rank $\geq$ fpocket $>$ blind (ESR1 93.2, 91.8, 83.7, and 47.1\%; MAPK1 65.8, 63.3, 56.6, and 45.9\%; PPARG 49.4, 46.6, 42.0, and 32.8\%). The co-crystallized-ligand and P2Rank boxes reached the site significantly more often than blind docking in each target (paired differences 13.8--46.1 points; every interval excluded zero), whereas fpocket was separated from blind docking only in ESR1. In ESR1, fpocket reached the site less often than the co-crystallized-ligand box (9.5 points, 95\% CI 4.4--15.8) and than P2Rank (8.1 points, 3.6--14.0).

Conditional on the Top-1 pose having reached the site, Top-1 success was similar across strategies within a target: 51.4, 56.7, 53.6, and 52.0\% in ESR1; 8.5, 9.0, 8.9, and 10.0\% in MAPK1; and 10.5, 12.3, 11.1, and 8.8\% in PPARG (Figure~\ref{fgr:cross}c; no pose that missed the site succeeded). The differences between strategies in cross-docking therefore arise mainly from how often the pose reaches the site and not from pose quality once it is there. Best-pose exceeded Top-1 by 21--22 points in ESR1 with the three localized boxes (5 points with blind docking) and by 4--10 points in MAPK1 and PPARG, so pose selection, more than pose generation, limited Top-1 where the native pose was generated at all.

\paragraph{Why MAPK1 and PPARG sit near the floor.} A rigid receptor must accommodate a reference pose that was superposed from another structure. In 94\% (MAPK1), 85\% (PPARG), and 86\% (ESR1) of the non-cognate pairs at least one reference heavy atom lay within 2.5~\AA\ of a receptor heavy atom. Top-1 success with the co-crystallized-ligand box was higher when the reference was clash-free: 77.2\% (88 of 114) against 43.3\% (315 of 727) in ESR1, 19.2\% (5 of 26) against 2.7\% (4 of 148) in PPARG, and 0 of 12 against 6.0\% (11 of 184) in MAPK1 (Table~S10). The clash-free subsets in MAPK1 and PPARG are small, so this explains the gap rather than quantifying it. The flagged non-equivalent-site pairs (14 in MAPK1, 36 in PPARG) were analysed separately: Top-1 success was 1 of 14 in MAPK1 under every strategy and 0--2 of 36 in PPARG, and the Top-1 pose reached the reference site in 2--5 of 14 and 4--7 of 36 (Table~S10).

Cross-docking therefore supports a narrower statement than the redocking benchmark. The strategy effect is present and has one consistent mechanism, namely how often the Top-1 pose reaches the site (blind docking 17--46 points below the co-crystallized-ligand box in the three targets), and once the site is reached the strategies are alike. It does not support a Top-1 strategy ordering in MAPK1 or PPARG, where rigid-receptor docking against superposed reference poses is near the floor for every strategy, including the box that carries the most information. Pooling the targets (pair-weighted Top-1 34.9, 34.5, 35.8, and 18.2\%; unweighted mean of the three targets 19.6, 19.2, 20.0, and 10.7\%) is descriptive, and the pair-weighted figure is dominated by ESR1.

\begin{table*}[htbp]
  \caption{Cross-Docking Success by Target and Search-Space Strategy (Exhaustiveness 16, Seed 42)$^{a}$}
  \label{tbl:cross}
  \centering
  \footnotesize
  \begin{tabular}{llcccccc}
    \toprule
    & & \multicolumn{2}{c}{Cognate (count)} & \multicolumn{3}{c}{Non-cognate, \% (95\% CI)} \\
    \cmidrule(lr){3-4}\cmidrule(lr){5-7}
    Target & Strategy & Top-1 & Best-pose & Top-1 & Best-pose & Site capture$^{b}$ \\
    \midrule
@@MAIN@@
    \bottomrule
  \end{tabular}

  \raggedright\footnotesize
  \textsuperscript{\emph{a}} Primary analysis (flagged non-equivalent-site pairs excluded). Non-cognate pairs: ESR1 841, MAPK1 196, PPARG 174 per strategy; cognate pairs: ESR1 29, MAPK1 15, PPARG 15. Success = sPyRMSD $<$ 2.0~\AA\ against the reference pose superposed into the receptor frame; averages are over all receptors. Confidence intervals: two-way cluster bootstrap (receptors and ligands resampled independently) within target, 10{,}000 resamples; cognate pairs are given as counts. Pooled rows are descriptive: ``pairs'' weights every non-cognate pair equally; ``targets'' is the unweighted mean of the three targets. \textsuperscript{\emph{b}} Top-1 pose centroid within 4~\AA\ of the reference-pose centroid (descriptive, added after the first target was scored).
\end{table*}

\begin{figure*}[htbp]
  \centering
  \includegraphics[width=\textwidth]{Figure6_crossdocking}
  \caption{Cross-docking by target and search-space strategy (exhaustiveness 16, seed 42; non-cognate pairs in the primary analysis). (a) Top-1 success (RMSD $<$ 2.0~\AA) with 95\% two-way cluster-bootstrap intervals; open diamonds mark the Top-1 success of the cognate (redocking) pairs of the same target and strategy. (b) Fraction of Top-1 poses whose centroid lies within 4~\AA\ of the reference-pose centroid (site capture). (c) Top-1 success among the pairs whose Top-1 pose reached the site (proportions without intervals; $n$ in Table~S10).}
  \label{fgr:cross}
\end{figure*}

"""
span("\\subsection{Cross-Docking}\n","\\subsection{Repurposing Case Study}\n",cross_new)

# ---------------- RESULTS: REPURPOSING
rep_new=r"""\subsection{Repurposing Case Study}

Imatinib was docked into the ABL1 (1IEP), KIT (1T46), SRC (2SRC), and EGFR (1M17) structures with all four search-space strategies (Table~\ref{tbl:repurp}). With the co-crystallized-ligand box, the Top-1 pose in KIT, the repurposed target, was 0.78~\AA\ from the crystal pose (Vina score $-12.61$~kcal~mol$^{-1}$), and in ABL1, the original target, 0.78~\AA\ ($-12.31$~kcal~mol$^{-1}$). Because this box is defined from the imatinib molecule in the same structure, this shows that the preparation and scoring chain recovers the pose when the site is given; it does not show that the site would be found.

Without that information the outcome depended on the strategy and the target. Blind docking recovered the ABL1 pose (Top-1 0.80~\AA, $-12.32$~kcal~mol$^{-1}$), whereas fpocket and P2Rank did not (Top-1 13.15 and 12.40~\AA; no saved pose closer than 11.7~\AA). In KIT, none of the three strategies recovered the crystal pose (Top-1 12.75, 12.82, and 11.83~\AA\ for fpocket, P2Rank, and blind; best saved pose 6.75, 10.28, and 11.72~\AA). For fpocket and P2Rank the Top-1 pose in KIT lay within 3.6~\AA\ of the crystal-ligand centroid, so the site was reached but the orientation was not recovered. Imatinib's crystal pose was therefore recovered in 3 of 8 runs for its two known targets, and for the repurposed target KIT only when the site was given. These are single-seed runs; blind docking was the least seed-stable strategy on the Astex set, so the ABL1 blind result should not be read as an advantage of blind docking.

Top-1 scores for the two kinases that imatinib inhibits were lower than those of the two it does not under the co-crystallized-ligand, fpocket, and P2Rank boxes, with margins as small as 0.3~kcal~mol$^{-1}$ (ABL1 $-10.58$ against SRC $-10.30$ with fpocket), and the ordering did not hold under blind docking (KIT $-9.74$ against SRC $-9.83$). The boxes of the negative controls were built from ATP-site ligands of a different size and scores are not comparable across proteins, so these comparisons are descriptive and not tested. The study shows retrospective recovery of a known target and pose and how it depends on the search-space strategy. It does not show that docking predicts repurposing or activity, and SRC and EGFR are weak negative controls (autoinhibited and DFG-in pockets).

\begin{table}[htbp]
  \caption{Repurposing Case Study: Imatinib Docked into Four Kinase Structures (Exhaustiveness 16, Seed 42)$^{a}$}
  \label{tbl:repurp}
  \centering
  \footnotesize\setlength{\tabcolsep}{3pt}
  \begin{tabular}{llccccc}
    \toprule
    Kinase (PDB) & IC$_{50}$ & Strategy & Score & Top-1 RMSD & Best RMSD & Recovered \\
    \midrule
    ABL1 (1IEP) & 38~nM\textsuperscript{\emph{b}} & Co-crystallized & $-$12.31 & 0.78 & 0.78 & Yes \\
     & & fpocket & $-$10.58 & 13.15 & 11.78 & No \\
     & & P2Rank & $-$10.65 & 12.40 & 11.66 & No \\
     & & Blind & $-$12.32 & 0.80 & 0.80 & Yes \\
    \addlinespace
    KIT (1T46) & 100~nM\textsuperscript{\emph{b}} & Co-crystallized & $-$12.61 & 0.78 & 0.78 & Yes \\
     & & fpocket & $-$11.63 & 12.75 & 6.75 & No \\
     & & P2Rank & $-$11.43 & 12.82 & 10.28 & No \\
     & & Blind & $-$9.74 & 11.83 & 11.72 & No \\
    \addlinespace
    SRC (2SRC) & $>$100~$\mu$M\textsuperscript{\emph{b}} & Co-crystallized & $-$9.53 & n/a & n/a & n/a \\
     & & fpocket & $-$10.30 & n/a & n/a & n/a \\
     & & P2Rank & $-$10.04 & n/a & n/a & n/a \\
     & & Blind & $-$9.83 & n/a & n/a & n/a \\
    \addlinespace
    EGFR (1M17) & $>$100~$\mu$M\textsuperscript{\emph{b}} & Co-crystallized & $-$9.62 & n/a & n/a & n/a \\
     & & fpocket & $-$9.60 & n/a & n/a & n/a \\
     & & P2Rank & $-$9.58 & n/a & n/a & n/a \\
     & & Blind & $-$9.61 & n/a & n/a & n/a \\
    \bottomrule
  \end{tabular}

  \raggedright\footnotesize
  \textsuperscript{\emph{a}} Single run per receptor and strategy (seed 42); reference = crystal imatinib in 1IEP and 1T46; recovery = Top-1 sPyRMSD $<$ 2.0~\AA. Score = Top-1 Vina score (kcal~mol$^{-1}$); RMSD in \AA; best RMSD is over the saved poses (within 3~kcal~mol$^{-1}$ of the best). Vina scores are not comparable across receptors. \textsuperscript{\emph{b}} Biochemical IC$_{50}$ as recorded in ChEMBL (compound CHEMBL941): ABL1 (v-Abl) 38~nM and SRC and EGFR $>$100~$\mu$M from \cite{Zimmermann1997}; KIT (recombinant c-Kit) 100~nM from a 2005 \textit{J. Med. Chem.} report (ChEMBL document CHEMBL1141098; \todo{cite after confirming the title}). Reported values depend on the assay: wild-type biochemical IC$_{50}$ values in ChEMBL span 19.9--7{,}700~nM for ABL1 and 16--10{,}000~nM for KIT.
\end{table}

"""
i=t.index("\\subsection{Repurposing Case Study}\n"); j=t.index("\\subsection{Positioning Relative to Existing Docking Workflows}")
t=t[:i]+rep_new+t[j:]

# ---------------- POSITIONING paragraph
rep("""would introduce additional sources of variation.
""","""would introduce additional sources of variation.

What the benchmark can say about these tools is indirect but quantifiable. Table~\\ref{tbl:compare-open} lists tools that offer a single search-space rule (BatchVina-GUI, Anyone Can Dock) and tools that offer several; the evidence here is that the choice among such rules is not neutral. Under one Vina protocol, the rules available when the ligand site is unknown (fpocket, P2Rank, blind) spanned 29.0--51.4\\% Top-1 success on the Astex set at exhaustiveness 16, and in cross-docking the co-crystallized-ligand and P2Rank boxes placed the Top-1 pose at the reference site significantly more often than blind docking in every target. A workflow that fixes one rule therefore fixes part of its answer, and a workflow that does not record which rule was used leaves that part unrecoverable. These figures are operating points under a single protocol and not predictions of how another platform would perform, because preparation and engine differences remain.
""")

# ---------------- LIMITATIONS
rep("""The cross-docking set contains three targets and 60 structures, so the
results are target-specific and are not pooled inferentially; ligand
charge states for ten repaired ligands were inferred; reference poses
depend on a sequence-based superposition, and pairs whose binding site
is not equivalent were analysed separately; covalent ligands were
excluded.""","""Cross-docking used a single seed and exhaustiveness 16, so seed
sensitivity was characterized on the Astex set only; the set has three
targets and 60 structures, so results are reported per target and pooled
only descriptively. Top-1 success in MAPK1 and PPARG was near the floor
for every strategy, so those targets cannot rank strategies by Top-1 pose
recovery; the reference poses are superposed from other structures and
85--94\\% of them clash with the receptor at 2.5~\\AA, which a rigid
receptor cannot accommodate, so these results reflect the rigid-receptor
limit more than the search-space strategy. The site-capture and clash
analyses were added after the first target had been scored and are
descriptive. Ligand charge states for ten repaired ligands were inferred;
two PPARG reference bond tables lacked a ring bond and were scored with a
geometry-perceived graph; reference poses depend on a sequence-based
superposition, and pairs whose binding site is not equivalent were
analysed separately; covalent ligands were excluded.""")
rep("""The repurposing case is retrospective and uses one compound;""","""The repurposing case is retrospective, uses one compound and a single seed, and with
the co-crystallized-ligand box the site is given;""")

# ---------------- CONCLUSIONS
rep("""In cross-docking, \\todo{one sentence: cognate versus non-cognate Top-1 and
whether the ordering held}. In a retrospective case, imatinib's pose was
\\todo{recovered / not recovered} in KIT and ABL1 by \\todo{strategies}.""","""In cross-docking on three targets, Top-1 success with the co-crystallized-ligand box fell from 93\\%, 67\\%, and 67\\% for cognate pairs to 48\\%, 6\\%, and 5\\% for non-cognate pairs; the strategies differed mainly in how often the Top-1 pose reached the reference site, blind docking being the least likely to do so in every target, and a Top-1 ordering could be resolved only in ESR1. In a retrospective case, imatinib's crystal pose was recovered to 0.78~\\AA\\ in KIT and ABL1 when the site was given, and \\tbd\\ when it was not.""")

# ---------------- SI LIST in main
rep("""  \\item Tables S8--S12 \\todo{fill when available}: box geometry; naive versus
  symmetry-aware per-pose calls; cross-docking (pair list, tier, flag,
  per-pair outcome); repurposing (per run, per strategy); ligand
  preparation and repair record for the 60 cross-docking ligands.""","""  \\item Tables S8--S9 \\todo{fill from the Astex per-pose data}: box
  geometry; naive versus symmetry-aware per-pose calls.
  \\item Table S10: cross-docking summaries (tier, flagged pairs, site
  capture, reference-pose clash, paired differences) and per-pair
  outcomes for all 1{,}320 pairs under four strategies.
  \\item Table S11: repurposing, four strategies
  and the imatinib activity record.
  \\item Table S12: ligand preparation and repair record for the 60
  cross-docking ligands.""")
rep("Figures S1--S4 \\todo{regenerate from the grid}:","Figures S1--S4:")
open(S+"v7/IBDock_JCIM_v7.tex","w").write(t.replace("@@MAIN@@",tab("main_table")))
print("main ok")
