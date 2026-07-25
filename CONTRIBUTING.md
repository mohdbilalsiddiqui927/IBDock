# Contributing to IBDock

Contributions are welcome — bug reports, fixes, and feature suggestions all help.

## Found a bug?

Open an issue and include:

- Your OS, Python version, and the versions of MGLTools, Vina, and Open Babel
- What you did and what you expected to happen
- The full error message or traceback from the terminal or Streamlit panel
- The input PDB or ligand file if the bug is file-specific (strip any confidential data)

The most common issues are wrong tool paths in `IBDock_config.json`, WSL path
formatting on Windows, and ligand filenames with special characters. Worth checking
those first.

## Want to suggest a feature?

Open an issue describing what you want to do and why. If it is a significant change,
discuss it before writing code — it saves everyone time.

## Submitting a fix or feature

```bash
git clone https://github.com/mohdbilalsiddiqui927/IBDock.git
cd IBDock
pip install -r requirements.txt
pip install pytest
pytest Tests/   # make sure everything passes before you start
```

Keep pull requests small and focused. One thing per PR. Add or update a test in
`Tests/test_ibdock.py` if your change touches a computational function.

A few things that will make review straightforward:

- Keep UI logic (Streamlit calls) and computation logic in separate functions — the computational ones need to be testable without a running Streamlit session
- Use type hints on new functions, consistent with the existing code
- Do not add new pip dependencies without raising it in an issue first

## Questions

Open a GitHub Discussion rather than an issue if you are not sure about something
or want to talk through an approach before writing code.
