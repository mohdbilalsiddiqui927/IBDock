# Contributing to IBDock Provenire

Thanks for your interest in IBDock Provenire. This project accompanies a
peer-reviewed manuscript, so the focus is on correctness and reproducibility
of the docking/validation pipeline described there.

## Reporting issues

If you hit a bug, a crash, or unexpected behavior:

1. Check that you're on the version described in [`requirements.txt`](requirements.txt)
   and the [System Requirements](README.md#system-requirements) table in the README.
2. Open a [GitHub Issue](../../issues) including:
   - Your OS and Python version
   - The exact steps that triggered the problem
   - Any error message or traceback, in full
   - If relevant, the input files (or a minimal example) that reproduce it

## Suggesting features

Open an issue describing the use case before submitting a pull request, so
the change can be discussed first — this avoids wasted effort on features
that don't fit the project's scope (batch AutoDock Vina docking with
search-space provenance and reference-based pose validation).

## Submitting changes

1. Fork the repository and create a branch from `main`.
2. Keep changes focused — one feature or fix per pull request.
3. Test your changes locally before submitting (`streamlit run "IBDock provenire.py"`,
   and if relevant, `python benchmark/scripts/statistical_analysis.py --write-template`
   to confirm the analysis pipeline still runs).
4. Describe what the change does and why in the pull request description.
5. Be responsive to review feedback — this keeps the pipeline's outputs
   consistent with what's reported in the accompanying manuscript.

## Code of conduct

Be respectful and constructive. This is an academic research tool used by a
small team; keep discussion focused on the software and the science.

## Questions

For questions about the methodology behind IBDock Provenire (search-space
strategies, pose validation routes, benchmarking approach), see the
accompanying manuscript, cited in [`CITATION.cff`](CITATION.cff). For
questions about using the software itself, open an issue.
