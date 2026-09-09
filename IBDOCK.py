# =============================================================================
# IBDock — Batch Molecular Docking with AutoDock Vina
#  A user-friendly Streamlit application for virtual screening
#
#  Requirements:
#      pip install streamlit numpy pandas matplotlib seaborn rdkit spyrmsd
#  External tools (portable bundled paths preferred; configurable in sidebar):
#      • MGLTools 1.5.7  — protein & ligand preparation
#      • Open Babel 3.x  — ligand format conversion
#      • AutoDock Vina   — docking engine
#      • fpocket / P2Rank — pocket detection (bundled native Windows tools)
#
#  Usage:
#      streamlit run IBDock_PORTABLE.py
# =============================================================================

import io
import json
import os
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import streamlit as st

# Optional at import time, required only for the Validation/RMSD module.
# Installation: pip install rdkit spyrmsd
try:
    from rdkit import Chem as _Chem
    from spyrmsd import molecule as _spymol
    from spyrmsd import rmsd as _spyrmsd_rmsd
    _RMSD_DEPS_OK = True
    _RMSD_DEPS_ERROR = ""
except ImportError as _rmsd_import_error:
    _Chem = None
    _spymol = None
    _spyrmsd_rmsd = None
    _RMSD_DEPS_OK = False
    _RMSD_DEPS_ERROR = str(_rmsd_import_error)

# ──────────────────────────────────────────────────────────────────────────────
# CONFIG PERSISTENCE  (saved per-project in project_dir/IBDock_config.json)
# ──────────────────────────────────────────────────────────────────────────────
import platform as _platform
import sys as _sys

_os = _platform.system()

# -----------------------------------------------------------------------------
# PORTABLE TOOL ROOT
# -----------------------------------------------------------------------------
# For the distributed Windows build, external scientific tools live beside the
# application in: IBDock\tools\.  This deliberately uses an ONEDIR build so
# the bundled MGLTools/P2Rank/OpenBabel files remain ordinary files on disk.
if getattr(_sys, "frozen", False):
    APP_DIR = Path(_sys.executable).resolve().parent
else:
    APP_DIR = Path(__file__).resolve().parent

# Resolve bundled scientific tools in both source and PyInstaller ONEDIR builds.
# PyInstaller places bundled data under _internal.
if getattr(_sys, "frozen", False):
    _tool_candidates = [
        APP_DIR / "_internal" / "tools",
        APP_DIR / "tools",
    ]
    BUNDLED_TOOLS = next(
        (p for p in _tool_candidates if p.exists()),
        APP_DIR / "_internal" / "tools",
    )
else:
    BUNDLED_TOOLS = APP_DIR / "tools"


def _first_existing(*paths):
    for p in paths:
        p = Path(p)
        if p.exists():
            return str(p)
    return ""


def _bundled_defaults():
    """Return portable defaults, preferring bundled tools over system installs."""
    mgl = BUNDLED_TOOLS / "MGLTools"
    obabel = _first_existing(
        BUNDLED_TOOLS / "OpenBabel" / "obabel.exe",
        mgl / "OpenBabel-2.3.2" / "obabel.exe",
    )
    vina = _first_existing(BUNDLED_TOOLS / "Vina" / "vina.exe")
    fpocket = _first_existing(
        BUNDLED_TOOLS / "fpocket-master" / "bin" / "fpocket.exe",
        BUNDLED_TOOLS / "fpocket" / "fpocket.exe",
        BUNDLED_TOOLS / "fpocket" / "fpocket.bat",
    )
    p2rank = _first_existing(
        BUNDLED_TOOLS / "P2Rank" / "prank.bat",
        BUNDLED_TOOLS / "P2Rank" / "prank.exe",
    )

    # User's currently validated Windows installation paths are retained as a
    # fallback while the portable tools are being assembled.
    if _os == "Windows":
        mgl_root = Path(r"C:\Program Files (x86)\MGLTools-1.5.7")
        mgl_python = _first_existing(mgl / "python.exe", mgl_root / "python.exe")
        prep_rec = _first_existing(
            mgl / "Lib" / "site-packages" / "AutoDockTools" / "Utilities24" / "prepare_receptor4.py",
            mgl_root / "Lib" / "site-packages" / "AutoDockTools" / "Utilities24" / "prepare_receptor4.py",
        )
        prep_lig = _first_existing(
            mgl / "Lib" / "site-packages" / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py",
            mgl_root / "Lib" / "site-packages" / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py",
        )
        obabel = obabel or _first_existing(
            mgl_root / "OpenBabel-2.3.2" / "obabel.exe",
        )
        

        # Pocket detection is fully native on Windows.
        fpocket = fpocket or _first_existing(
            BUNDLED_TOOLS / "fpocket-master" / "bin" / "fpocket.exe",
            BUNDLED_TOOLS / "fpocket" / "fpocket.exe",
        )
        p2rank = p2rank or str(BUNDLED_TOOLS / "P2Rank" / "prank.bat")

        return dict(
            mgl_python=mgl_python,
            prep_rec=prep_rec,
            prep_lig=prep_lig,
            vina_path=vina,
            obabel_path=obabel,
            fpocket_path=fpocket,
            p2rank_path=p2rank,
        )

    if _os == "Darwin":
        return dict(
            mgl_python="/opt/homebrew/bin/pythonsh",
            prep_rec="/opt/homebrew/share/mgltools/AutoDockTools/Utilities24/prepare_receptor4.py",
            prep_lig="/opt/homebrew/share/mgltools/AutoDockTools/Utilities24/prepare_ligand4.py",
            vina_path="/opt/homebrew/bin/vina",
            obabel_path="/opt/homebrew/bin/obabel",
            fpocket_path="",
            p2rank_path="",
        )

    return dict(
        mgl_python="/usr/bin/pythonsh",
        prep_rec="/usr/share/mgltools/AutoDockTools/Utilities24/prepare_receptor4.py",
        prep_lig="/usr/share/mgltools/AutoDockTools/Utilities24/prepare_ligand4.py",
        vina_path="/usr/bin/vina",
        obabel_path="/usr/bin/obabel",
        fpocket_path="fpocket",
        p2rank_path="prank",
    )


_DEF = _bundled_defaults()

def _load_config(project_dir: Path) -> dict:
    cfg_file = project_dir / "IBDock_config.json"

    if cfg_file.exists():
        try:
            cfg = {**_DEF, **json.loads(cfg_file.read_text())}
        except Exception:
            cfg = dict(_DEF)
    else:
        cfg = dict(_DEF)

    # On Windows, always prefer the bundled native pocket tools.
    if _os == "Windows":
        bundled_p2rank = _first_existing(
            BUNDLED_TOOLS / "P2Rank" / "prank.bat",
            BUNDLED_TOOLS / "P2Rank" / "prank.exe",
            BUNDLED_TOOLS / "P2Rank" / "prank",
        )
        bundled_fpocket = _first_existing(
            BUNDLED_TOOLS / "fpocket-master" / "bin" / "fpocket.exe",
            BUNDLED_TOOLS / "fpocket" / "fpocket.exe",
        )

        if bundled_p2rank:
            cfg["p2rank_path"] = bundled_p2rank

        if bundled_fpocket:
            cfg["fpocket_path"] = bundled_fpocket

    return cfg

def _save_config(project_dir: Path, cfg: dict):
    try:
        (project_dir / "IBDock_config.json").write_text(json.dumps(cfg, indent=2))
    except Exception:
        pass

# ──────────────────────────────────────────────────────────────────────────────
# PAGE CONFIG
# ──────────────────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="IBDock — Molecular Docking",
    page_icon="",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ──────────────────────────────────────────────────────────────────────────────
# GLOBAL STYLE  — JOSS-aligned academic design
# ──────────────────────────────────────────────────────────────────────────────
_CSS = """
<style>
/* ── Fonts ── */
@import url('https://fonts.googleapis.com/css2?family=Source+Serif+4:ital,wght@0,300;0,400;0,600;1,400&family=IBM+Plex+Mono:wght@400;600&family=Inter:wght@400;500;600&display=swap');

/* ── Base ── */
html, body, [class*="css"] {
    font-family: 'Inter', 'Helvetica Neue', Arial, sans-serif;
    font-size: 14px;
    color: #1a1a1a;
}

/* ── Sidebar ── */
[data-testid="stSidebar"] {
    background: #f7f7f5;
    border-right: 1px solid #ddd;
}
[data-testid="stSidebar"] * { color: #2c2c2c !important; }
[data-testid="stSidebar"] input {
    background: #fff !important;
    border: 1px solid #ccc !important;
    color: #1a1a1a !important;
    border-radius: 3px !important;
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.72rem !important;
}
[data-testid="stSidebar"] .stButton button {
    background: #1a4a7a !important;
    color: #ffffff !important;
    border: none !important;
    border-radius: 3px !important;
    font-weight: 500 !important;
    font-size: 0.82rem !important;
    letter-spacing: 0.01em !important;
}
[data-testid="stSidebar"] .stButton button:hover {
    background: #0f3560 !important;
    color: #ffffff !important;
}
[data-testid="stSidebar"] .stButton button p,
[data-testid="stSidebar"] .stButton button span,
[data-testid="stSidebar"] .stButton button div {
    color: #ffffff !important;
}

/* ── Main area ── */
.main .block-container { padding-top: 1rem; max-width: 1180px; }

/* ── Primary buttons ── */
.stButton > button[kind="primary"] {
    background: #1a4a7a !important;
    color: #fff !important;
    border: none !important;
    border-radius: 3px !important;
    font-weight: 500 !important;
    font-size: 0.88rem !important;
    padding: 0.5rem 1.8rem !important;
    letter-spacing: 0.01em !important;
    transition: background 0.15s !important;
}
.stButton > button[kind="primary"]:hover {
    background: #0f3560 !important;
}

/* ── Tabs — compact, academic ── */
.stTabs [data-baseweb="tab-list"] {
    gap: 0;
    border-bottom: 2px solid #ccc;
    background: transparent;
}
.stTabs [data-baseweb="tab"] {
    font-family: 'Inter', sans-serif !important;
    font-weight: 500 !important;
    font-size: 0.80rem !important;
    color: #555 !important;
    padding: 0.45rem 0.9rem !important;
    border-radius: 0 !important;
    background: transparent !important;
    border: none !important;
    text-transform: uppercase;
    letter-spacing: 0.04em;
    white-space: nowrap;
}
.stTabs [aria-selected="true"] {
    color: #1a4a7a !important;
    border-bottom: 2px solid #1a4a7a !important;
    margin-bottom: -2px;
    background: transparent !important;
}
.stTabs [data-baseweb="tab"]:hover {
    color: #1a4a7a !important;
    background: #f0f4f8 !important;
}

/* ── Metrics ── */
[data-testid="stMetric"] {
    background: #fafafa;
    border: 1px solid #e0e0e0;
    border-radius: 3px;
    padding: 0.7rem 1rem;
}
[data-testid="stMetricLabel"] { font-size: 0.72rem !important; color: #666 !important; font-weight: 600 !important; text-transform: uppercase; letter-spacing: 0.06em; }
[data-testid="stMetricValue"] { font-family: 'IBM Plex Mono', monospace !important; font-size: 1.35rem !important; color: #1a1a1a !important; }

/* ── Expander ── */
.streamlit-expanderHeader {
    font-weight: 500 !important;
    color: #1a1a1a !important;
    font-size: 0.88rem !important;
}

/* ── Dataframe ── */
[data-testid="stDataFrame"] { border: 1px solid #ddd; border-radius: 3px; overflow: hidden; }

/* ── Step indicator ── */
.step-bar { display:flex; gap:0; margin-bottom:1.4rem; border-top: 1px solid #ddd; border-bottom: 1px solid #ddd; }
.step-item { flex:1; padding:0.5rem 0.4rem; text-align:center; font-family:'Inter',sans-serif; font-size:0.72rem; font-weight:600; letter-spacing:0.04em; text-transform:uppercase; border-right:1px solid #ddd; cursor:default; }
.step-item:last-child { border-right:none; }
.step-done  { background:#e8f4ee; color:#1a6640; }
.step-active{ background:#e8eef6; color:#1a4a7a; border-bottom: 2px solid #1a4a7a; }
.step-idle  { background:#fafafa; color:#999; }
.step-num   { display:inline-block; width:18px; height:18px; border-radius:50%; font-size:0.65rem; line-height:18px; margin-right:5px; background:#ddd; color:#555; vertical-align:middle; }
.step-active .step-num { background:#1a4a7a; color:#fff; }
.step-done .step-num   { background:#1a6640; color:#fff; }
.step-label { vertical-align:middle; }

/* ── Info cards (landing / guide) ── */
.info-card { background:#fafafa; border:1px solid #e0e0e0; border-radius:3px; padding:1rem 1.2rem; margin-bottom:0.7rem; }
.info-card h4 { margin:0 0 0.35rem; font-size:0.88rem; color:#1a1a1a; font-weight:600; font-family:'Inter',sans-serif; }
.info-card p  { margin:0; font-size:0.82rem; color:#555; line-height:1.6; }

/* ── Preflight checklist ── */
.pf-row { display:flex; align-items:center; gap:10px; padding:0.4rem 0.7rem; border-radius:3px; margin-bottom:5px; font-size:0.83rem; font-weight:500; }
.pf-ok  { background:#f0f8f4; border:1px solid #b8dfc8; color:#1a5c38; }
.pf-warn{ background:#fdf8ed; border:1px solid #e8d49a; color:#7a5800; }
.pf-err { background:#fdf2f2; border:1px solid #e8b8b8; color:#8a1010; }

/* ── Guide accordion ── */
.guide-param { background:#fafafa; border-left:3px solid #1a4a7a; border-radius:0 3px 3px 0; padding:0.55rem 0.85rem; margin-bottom:0.45rem; font-size:0.82rem; }
.guide-param strong { color:#1a1a1a; display:block; margin-bottom:2px; font-family:'Inter',sans-serif; }
.guide-param span   { color:#555; line-height:1.55; }
.guide-tag { display:inline-block; background:#e8eef6; color:#1a4a7a; border:1px solid #b8cce4; border-radius:2px; font-size:0.66rem; font-weight:600; padding:1px 7px; margin-left:6px; vertical-align:middle; text-transform:uppercase; letter-spacing:0.04em; }
.guide-tag.adv { background:#fdf5ec; color:#8a4800; border-color:#e8c898; }

/* ── Section headings ── */
h3 { font-family: 'Source Serif 4', Georgia, serif !important; font-size: 1.25rem !important; font-weight: 600 !important; color: #1a1a1a !important; letter-spacing: -0.01em !important; }
h4 { font-family: 'Inter', sans-serif !important; font-size: 0.95rem !important; font-weight: 600 !important; color: #333 !important; }

/* ── Code / mono ── */
code { font-family: 'IBM Plex Mono', monospace !important; font-size: 0.82em; background: #f0f0ee; padding: 1px 5px; border-radius: 2px; }
</style>
"""

# Inject CSS — try st.html (Streamlit 1.31+), fall back to components
try:
    st.html(_CSS)
except AttributeError:
    import streamlit.components.v1 as _comp_css
    _comp_css.html(_CSS, height=0, scrolling=False)

# ──────────────────────────────────────────────────────────────────────────────
# CONSTANTS
# ──────────────────────────────────────────────────────────────────────────────
METAL_IONS = {"ZN", "MG", "FE", "CA", "MN", "CU", "K", "NA", "NI", "CO",
              "CD", "HG", "PT", "AU", "AG", "BA", "SR", "LI", "RB", "CS"}

# Residue names that are NEVER a drug-like co-crystallised ligand.
# These are crystallographic additives, buffer components, solvents, ions,
# and cryoprotectants that appear as HETATM records in PDB files.
# When computing the grid box centre from co-crystallised ligands, all of
# these are excluded so only the true drug-like ligand(s) are used.
EXCLUDE_RESNAMES = {
    # Water
    "HOH", "WAT", "H2O", "DOD",
    # Ions — monatomic
    "NA",  "MG",  "K",   "CA",  "MN",  "FE",  "CO",  "NI",  "CU",  "ZN",
    "SE",  "BR",  "CD",  "I",   "CS",  "BA",  "HG",  "PT",  "AU",  "AG",
    "LI",  "RB",  "SR",  "PB",  "TL",  "IN",  "YB",  "SM",  "TB",  "GD",
    "CL",  "F",   "IOD",
    # Common polyatomic ions / buffer salts
    "SO4", "PO4", "NO3", "ACT", "ACY", "FMT", "SUC", "TRS", "MES",
    "HEP", "EPE", "BIS", "CIT", "TAR", "OXL", "MLT", "PHO",
    # Cryoprotectants and precipitants
    "GOL", "GLY", "PG4", "PGE", "PG6", "PEG", "PE4", "P6G", "1PE", "2PE",
    "EDO", "EGL", "MPD", "IPA", "EOH", "ACE", "ACN", "DMS", "MSO",
    "DMF", "DMU", "IMD",
    # Reducing agents and their reaction adducts with surface cysteines —
    # NOT ligands; excluding them prevents false grid-centroid shifts.
    "BME", "CME", "DTT", "TCEP", "MRD", "B3P", "BTB", "PGO",
    # Detergents common in membrane protein crystals
    "BOG", "DDM", "OG",  "NG",  "LMT", "LDA",
    # Polyamines / crystallisation additives
    "SPD", "SPM", "PUT",
    # Monosaccharides often added to buffers
    "XYL", "FRU", "GAL", "MAN", "FUC",
}


# ──────────────────────────────────────────────────────────────────────────────
# HELPER UTILITIES
# ──────────────────────────────────────────────────────────────────────────────

def run_cmd(cmd: list):
    """Run a subprocess command, raising RuntimeError on failure."""
    subprocess.run(
        [str(c) for c in cmd],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=False,
    )


def validate_tool(path_str: str, label: str) -> dict:
    """Validate a native external tool."""
    result = {"label": label, "ok": False, "version": None, "detail": ""}

    if not path_str or not path_str.strip():
        result["detail"] = "Path not specified"
        return result

    path = Path(path_str.strip())

    if path.suffix.lower() == ".py":
        if path.exists():
            result["ok"] = True
            result["detail"] = f"Script found ✔  ({path.stat().st_size:,} bytes)"
        else:
            result["detail"] = f"File not found: {path}"
        return result

    if not path.exists():
        result["detail"] = f"File not found: {path}"
        return result

    if path.suffix.lower() not in {".bat", ".cmd"} and not os.access(str(path), os.X_OK):
        result["detail"] = f"File exists but is not executable: {path}"
        return result

    for flag in ["--version", "-version", "--help", ""]:
        try:
            cmd = [str(path)] + ([flag] if flag else [])
            if path.suffix.lower() in {".bat", ".cmd"}:
                cmd = ["cmd.exe", "/c"] + cmd
            proc = subprocess.run(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                text=True, timeout=10, shell=False
            )
            out = (proc.stdout + proc.stderr).strip()
            if out:
                result["version"] = next(
                    (l for l in out.splitlines() if l.strip()), ""
                )[:120]
                break
        except Exception:
            continue

    result["ok"] = True
    result["detail"] = "Executable found and runs ✔"
    return result


# ──────────────────────────────────────────────────────────────────────────────
# STRUCTURE HELPERS
# ──────────────────────────────────────────────────────────────────────────────

def extract_ligand_atoms(lines: list, min_heavy_atoms: int = 7):
    """Return (coords_array, resnames_set) for drug-like co-crystallised ligand atoms.

    Filtering strategy (in order):
    1. Skip any residue name in EXCLUDE_RESNAMES (ions, solvents, buffers,
       cryoprotectants).
    2. Group remaining HETATM records by (residue name, CHAIN) and count
       heavy atoms per group.
    3. Skip any group whose heavy atom count is below min_heavy_atoms
       (default 7). This catches small crystallographic additives not in
       EXCLUDE_RESNAMES such as glycerol fragments, acetate, formate, and
       single-atom ions.
    4. Restrict to a SINGLE chain: when a structure contains multiple
       crystallographic copies of the complex (common — e.g. 4 copies in
       the asymmetric unit), pooling ligand coordinates from every chain
       together produces a centroid that doesn't correspond to any real
       binding site (confirmed root cause of a ~10.7 A grid-centering
       error on PDB 1HVY, which has 4 chains each with its own ligand).
       Instead, pick the single chain with the most qualifying ligand
       atoms and use only that chain's groups — this still allows
       multiple co-located HETATM groups within ONE binding site (e.g. a
       metal cofactor + organic ligand) to be combined together, since
       those legitimately belong to the same site.

    Returns a numpy array of (x, y, z) coords and a set of surviving residue names.
    """
    from collections import defaultdict

    groups = defaultdict(list)   # (resname, chain) -> list of (x, y, z) for heavy atoms
    for line in lines:
        if not line.startswith("HETATM"):
            continue
        rn = line[17:20].strip()
        if rn in EXCLUDE_RESNAMES:
            continue
        atom_name = line[12:16].strip()
        if atom_name.upper().startswith(("H", "D")):
            continue                       # skip hydrogens / deuteriums
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            chain = line[21]
            groups[(rn, chain)].append(xyz)
        except ValueError:
            continue

    # Apply minimum size filter — drop anything too small to be a drug ligand
    qualifying = {key: xyzs for key, xyzs in groups.items() if len(xyzs) >= min_heavy_atoms}
    if not qualifying:
        return np.array([]), set()

    # Restrict to the single chain with the most total qualifying ligand
    # atoms — avoids pooling coordinates from multiple, spatially separate
    # crystallographic copies of the same complex.
    atoms_per_chain = defaultdict(int)
    for (rn, chain), xyzs in qualifying.items():
        atoms_per_chain[chain] += len(xyzs)
    best_chain = max(atoms_per_chain, key=atoms_per_chain.get)

    coords, resnames = [], set()
    for (rn, chain), xyzs in qualifying.items():
        if chain != best_chain:
            continue
        coords.extend(xyzs)
        resnames.add(rn)

    return (np.array(coords) if coords else np.array([])), resnames



def parse_pdbqt_heavy_atoms(text: str) -> list:
    """
    Extract heavy-atom lines from a PDBQT/PDB string.

    Hydrogen detection uses ONLY the atom-name field (cols 12-15).
    We deliberately do NOT check the last-token AD4/Vina type column
    because obabel assigns the type "H" to heavy atoms when it cannot
    perceive bond orders (e.g. raw HETATM lines without CONECT records
    from a crystal PDB), which would incorrectly filter out every heavy
    atom and leave ref=1.

    Standard PDB atom-name rules guarantee that hydrogen atoms always
    have names starting with H (e.g. H, HA, HB1, HD, HS) or with a
    digit followed by H (e.g. 1H, 2HB). Both cases are handled.
    """
    heavy = []
    for line in text.splitlines(keepends=True):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if len(line) < 17:
            continue
        atom_name = line[12:16].strip().upper() if len(line) > 16 else ""
        # Digit-prefixed names like "1H ", "2HB" are also hydrogens
        is_h = atom_name.startswith("H") or (
            len(atom_name) >= 2 and atom_name[0].isdigit() and atom_name[1] == "H"
        )
        if is_h:
            continue
        heavy.append(line)
    return heavy


def extract_protein_atoms(lines: list) -> np.ndarray:
    """Return Nx3 array of ATOM record coordinates."""
    atoms = []
    for line in lines:
        if line.startswith("ATOM"):
            try:
                atoms.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
    return np.array(atoms) if atoms else np.array([])


def clean_protein(lines: list, ligand_resnames: set) -> list:
    """Keep ATOM records and metal HETATMs; drop water and docked ligand atoms."""
    out = []
    for line in lines:
        if line.startswith("ATOM"):
            out.append(line)
        elif line.startswith("HETATM"):
            rn = line[17:20].strip()
            if rn in METAL_IONS:
                out.append(line)
            elif rn not in ligand_resnames and rn not in EXCLUDE_RESNAMES:
                out.append(line)
    out.append("END\n")
    return out


def round_even(v) -> int:
    return int(np.ceil(v / 2) * 2)


def cap_grid(size: list, max_size: int = 120) -> list:
    return [min(int(s), max_size) for s in size]


def fix_pdbqt_atom_names(pdbqt_path):
    """Truncate atom names longer than 4 characters (MGLTools bug)."""
    fixed = []
    with open(pdbqt_path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                name = line[12:16].strip()[:4]
                line = line[:12] + name.ljust(4) + line[16:]
            fixed.append(line)
    with open(pdbqt_path, "w") as fh:
        fh.writelines(fixed)


def strip_receptor_hydrogens(pdbqt_path):
    """
    Remove NONPOLAR hydrogen atoms from a receptor PDBQT, checked by AD4
    atom type (column 77-79), not atom name. MGLTools names many hydrogens
    with digit-prefixed identifiers (e.g. '1HD2', '2HG1') for geminal H's,
    which a name-based "starts with H" check misses entirely — confirmed
    empirically: 100% of H/HD atoms produced by this pipeline use that
    naming convention, so the previous name-based check silently removed
    nothing from any receptor.

    Polar hydrogens (AD4 type 'HD') are explicitly preserved — required for
    Vina's hydrogen-bond donor scoring term.
    """
    kept = []
    with open(pdbqt_path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")) and line[77:79].strip() == "H":
                continue  # drop nonpolar H only; keep HD (polar) and all heavy atoms
            kept.append(line)
    with open(pdbqt_path, "w") as fh:
        fh.writelines(kept)



# ──────────────────────────────────────────────────────────────────────────────
# POCKET DETECTION
# ──────────────────────────────────────────────────────────────────────────────

def _run_native_p2rank(pdb_file, p2rank_exec: str, padding: float = 5.0) -> tuple:
    """Run a native Windows P2Rank launcher and locate its prediction CSVs.

    This is used when a bundled prank.bat/prank.exe is present.  The function
    intentionally searches the generated files instead of assuming one exact
    P2Rank output layout, because P2Rank distributions can differ by version.
    """
    exe = Path(p2rank_exec)
    if not exe.exists():
        raise RuntimeError(f"P2Rank executable not found: {exe}")

    pdb = Path(pdb_file).resolve()
    work_dir = pdb.parent
    cmd = [str(exe), "predict", "-f", str(pdb)]
    kwargs = {"stdout": subprocess.PIPE, "stderr": subprocess.PIPE, "text": True, "cwd": str(work_dir)}
    if exe.suffix.lower() in {".bat", ".cmd"}:
        cmd = ["cmd.exe", "/c"] + cmd

    result = subprocess.run(cmd, **kwargs)
    if result.returncode != 0:
        raise RuntimeError(f"P2Rank failed: {(result.stderr or result.stdout).strip()}")

    stem = pdb.stem
    candidates = list(work_dir.rglob(f"{stem}*predictions.csv"))
    if not candidates:
        # Some distributions write into a fixed output directory under the
        # P2Rank installation. Search beside the executable as a second pass.
        candidates = list(exe.parent.rglob(f"{stem}*predictions.csv"))
    if not candidates:
        raise RuntimeError("P2Rank completed but no predictions.csv file was found.")

    pred_file = max(candidates, key=lambda p: p.stat().st_mtime)
    df = pd.read_csv(pred_file)
    df.columns = df.columns.str.strip()
    if df.empty or not {"center_x", "center_y", "center_z"}.issubset(df.columns):
        raise RuntimeError(f"Invalid P2Rank prediction file: {pred_file}")

    best = df.iloc[0]
    center = np.array([best["center_x"], best["center_y"], best["center_z"]], dtype=float)

    size = None
    points = list(pred_file.parent.glob(f"{stem}*_points.csv"))
    if points:
        try:
            pts_df = pd.read_csv(points[0])
            pts_df.columns = pts_df.columns.str.strip()
            if "pocket" in pts_df.columns:
                pts_df = pts_df[pts_df["pocket"] == 1]
            if len(pts_df) >= 4 and {"x", "y", "z"}.issubset(pts_df.columns):
                xyz = pts_df[["x", "y", "z"]].values.astype(float)
                extent = xyz.max(axis=0) - xyz.min(axis=0)
                size = cap_grid([round_even(v) for v in (extent + 2 * padding)])
        except Exception:
            size = None

    if size is None:
        size = cap_grid([round_even(v) for v in (np.array([12.0, 12.0, 12.0]) + 2 * padding)])

    return center, size


def run_p2rank(pdb_file, p2rank_exec: str, padding: float = 5.0) -> tuple:
    """Run the bundled native Windows P2Rank executable."""
    if not p2rank_exec or not str(p2rank_exec).strip():
        raise RuntimeError("P2Rank path is not configured.")

        raise RuntimeError(
            "Please use the bundled native P2Rank tool."
        )

    return _run_native_p2rank(pdb_file, p2rank_exec, padding)


def run_fpocket(pdb_file, fpocket_exec: str):
    """Run the bundled native Windows fpocket executable."""
    pdb_file = Path(pdb_file)

    if not pdb_file.exists():
        raise RuntimeError(f"fpocket input PDB not found: {pdb_file}")

    fpocket_exec = str(fpocket_exec or "").strip()
    if not fpocket_exec:
        raise RuntimeError("Native fpocket executable is not configured.")

        raise RuntimeError(
            "Please use the bundled native fpocket.exe."
        )

    exe = Path(fpocket_exec)
    if not exe.exists():
        raise RuntimeError(f"fpocket executable not found: {exe}")

    out_dir = pdb_file.parent / f"{pdb_file.stem}_out"
    run_dir = pdb_file.parent
    run_pdb = pdb_file

    # fpocket always creates "<input_stem>_out". On Windows, an older
    # fpocket output can be locked by Explorer, a viewer, antivirus, or
    # another process. Do not fail the entire pocket-detection method just
    # because that previous directory cannot be removed. Instead, run fpocket
    # in a unique temporary working directory using a copy of the PDB.
    if out_dir.exists():
        try:
            shutil.rmtree(out_dir)
        except (PermissionError, OSError):
            unique_tag = f"{os.getpid()}_{int(time.time() * 1000)}"
            run_dir = pdb_file.parent / f".{pdb_file.stem}_fpocket_{unique_tag}"
            run_dir.mkdir(parents=True, exist_ok=False)
            run_pdb = run_dir / pdb_file.name
            shutil.copy2(pdb_file, run_pdb)
            out_dir = run_dir / f"{pdb_file.stem}_out"

    result = subprocess.run(
        [str(exe), "-f", run_pdb.name],
        cwd=str(run_dir),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        errors="replace",
    )

    if result.returncode != 0:
        error_text = (
            result.stderr.strip()
            or result.stdout.strip()
            or "unknown fpocket error"
        )
        raise RuntimeError(f"fpocket failed: {error_text}")

    pocket_file = out_dir / "pockets" / "pocket1_atm.pdb"
    if not pocket_file.exists():
        raise RuntimeError(
            f"fpocket completed but pocket file was not found: {pocket_file}"
        )

    return out_dir


def get_fpocket_center_and_size(pocket_file, padding: float = 5.0) -> tuple:
    """Parse an fpocket pocket atom file and return (center, size).

    fpocket writes pocket1_atm.pdb containing the residue atoms lining the
    cavity.  The geometric extent of those atoms is a real measurement of the
    binding site — far more accurate than a hardcoded 26 Å constant.

    Returns:
        (center, size) — both lists/arrays of 3 Angstrom values.
    """
    coords = []
    with open(pocket_file) as fh:
        for line in fh:
            if line.startswith(("HETATM", "ATOM")):
                try:
                    coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except ValueError:
                    continue
    if not coords:
        raise RuntimeError(f"No atoms in fpocket file: {pocket_file}")
    arr    = np.array(coords)
    center = arr.mean(axis=0)
    extent = arr.max(axis=0) - arr.min(axis=0)
    size   = cap_grid([round_even(v) for v in (extent + 2 * padding)])
    return center, size


# ──────────────────────────────────────────────────────────────────────────────
# VINA DOCKING
# ──────────────────────────────────────────────────────────────────────────────

def run_single_docking(job: tuple):
    """
    Run one AutoDock Vina job.
    job = (receptor, ligand, config_file, vina_exe,
           exhaustiveness, num_modes, energy_range, cores, dock_dir, result_dir, seed)
    Returns ("success"|"failed", protein_name, ligand_name, message).
    """
    rec, lig, config_file, vina_exe, exhaustiveness, num_modes, \
        energy_range, cores_per_job, dock_dir, result_dir, seed = job

    protein_name = Path(rec).stem.replace("_receptor", "")
    ligand_name  = Path(lig).stem
    out_pdbqt    = dock_dir   / f"{protein_name}_{ligand_name}.pdbqt"
    log_txt      = result_dir / f"{protein_name}_{ligand_name}.txt"

    try:
        cmd = [
            str(vina_exe),
            "--receptor",      str(rec),
            "--ligand",        str(lig),
            "--config",        str(config_file),
            "--exhaustiveness",str(exhaustiveness),
            "--num_modes",     str(num_modes),
            "--energy_range",  str(energy_range),
            "--cpu",           str(cores_per_job),
            "--out",           str(out_pdbqt),
        ]
        # Fixed seed => reproducible results run-to-run (Vina otherwise seeds
        # from system entropy, so re-running the same job gives different
        # poses each time — important for a validation study others need to
        # be able to reproduce).
        if seed is not None:
            cmd += ["--seed", str(seed)]

        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, shell=False,
        )
        if result.returncode != 0:
            return ("failed", protein_name, ligand_name, result.stderr.strip())

        log_txt.write_text(result.stdout)

        if not out_pdbqt.exists() or out_pdbqt.stat().st_size == 0:
            return ("failed", protein_name, ligand_name, "Empty output file")

        return ("success", protein_name, ligand_name, "")

    except Exception as exc:
        return ("failed", protein_name, ligand_name, str(exc))



def parse_vina_log(log_file) -> float | None:
    """Extract best affinity (mode 1) from a Vina log file."""
    import re
    try:
        with open(log_file) as fh:
            for line in fh:
                m = re.match(r'^\s*1\s+(-?\d+\.\d+)', line)
                if m:
                    return float(m.group(1))
    except Exception:
        pass
    return None


def parse_vina_all_modes(log_file) -> list:
    """Return list of (mode_number, affinity) from a Vina log."""
    import re
    modes = []
    try:
        with open(log_file) as fh:
            for line in fh:
                m = re.match(r'^\s*(\d+)\s+(-?\d+\.\d+)\s+\S+\s+\S+', line)
                if m:
                    modes.append((int(m.group(1)), float(m.group(2))))
    except Exception:
        pass
    return sorted(modes, key=lambda x: x[0])


# ──────────────────────────────────────────────────────────────────────────────
# POSE EXTRACTION & RMSD
# ──────────────────────────────────────────────────────────────────────────────

def pdbqt_to_pdb_text(pdbqt_path, first_pose_only: bool = True) -> str:
    """
    Convert a PDBQT file to plain PDB text for 3D viewing.
    If first_pose_only, returns only the first MODEL block.
    """
    lines = Path(pdbqt_path).read_text(errors="replace").splitlines(keepends=True)
    has_models = any(l.startswith("MODEL") for l in lines)

    if not has_models:
        return "".join(l for l in lines if l.startswith(("ATOM", "HETATM", "END")))

    out, in_model, model_count = [], False, 0
    for line in lines:
        if line.startswith("MODEL"):
            model_count += 1
            if first_pose_only and model_count > 1:
                break
            in_model = True
            out.append(line)
        elif line.startswith("ENDMDL"):
            out.append(line)
            in_model = False
        elif in_model:
            out.append(line)
    return "".join(out)


def pdbqt_to_sdf_text(pdbqt_path, obabel_path: str, first_pose_only: bool = True, add_hydrogens: bool = True) -> str:
    """
    Convert a PDBQT ligand pose to SDF format using Open Babel.
    SDF preserves explicit bond tables (bond orders, ring perception),
    which py3Dmol needs to render complex ligands correctly.
    Falls back to empty string on failure (caller should use PDB fallback).
    """
    import tempfile
    pdbqt_path = str(pdbqt_path)

    # If we only want the first pose, write a temp PDBQT with just MODEL 1
    if first_pose_only:
        lines = Path(pdbqt_path).read_text(errors="replace").splitlines(keepends=True)
        has_models = any(l.startswith("MODEL") for l in lines)
        if has_models:
            first_block, in_model, done = [], False, False
            for line in lines:
                if done:
                    break
                if line.startswith("MODEL"):
                    in_model = True
                    first_block.append(line)
                elif line.startswith("ENDMDL"):
                    first_block.append(line)
                    done = True
                elif in_model:
                    first_block.append(line)
            tmp_pdbqt = tempfile.NamedTemporaryFile(suffix=".pdbqt", delete=False, mode="w")
            tmp_pdbqt.writelines(first_block)
            tmp_pdbqt.flush()
            tmp_pdbqt.close()
            pdbqt_path = tmp_pdbqt.name

    try:
        tmp_sdf = tempfile.NamedTemporaryFile(suffix=".sdf", delete=False)
        tmp_sdf.close()
        obabel_cmd = [obabel_path, pdbqt_path, "-O", tmp_sdf.name]
        if add_hydrogens:
            obabel_cmd.append("-h")

        result = subprocess.run(
            obabel_cmd,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
        )
        if result.returncode == 0 and Path(tmp_sdf.name).exists():
            sdf_content = Path(tmp_sdf.name).read_text(errors="replace")
            if sdf_content.strip():
                return sdf_content
        return ""
    except Exception:
        return ""
    finally:
        try:
            import os
            os.unlink(tmp_sdf.name)
        except Exception:
            pass
        if first_pose_only:
            try:
                import os
                os.unlink(pdbqt_path)
            except Exception:
                pass


def extract_vina_pose_by_mode(pdbqt_path: str, mode_number: int) -> list:
    """Return ATOM/HETATM lines for a specific Vina mode."""
    lines = Path(pdbqt_path).read_text(errors="replace").splitlines(keepends=True)
    current_model, in_target, coords = 0, False, []
    for line in lines:
        if line.startswith("MODEL"):
            try:
                current_model = int(line.split()[1])
            except (IndexError, ValueError):
                current_model += 1
            in_target = (current_model == mode_number)
        elif line.startswith("ENDMDL"):
            if in_target:
                break
            in_target = False
        elif in_target and line.startswith(("ATOM", "HETATM")):
            coords.append(line)
    if not coords and mode_number == 1:
        coords = [l for l in lines if l.startswith(("ATOM", "HETATM"))]
    return coords


def _count_pdbqt_models(pdbqt_path) -> int:
    """
    Count how many Vina pose (MODEL) blocks a docked PDBQT contains, so the
    best-pose RMSD sweep below knows how many modes to try. Falls back to 1
    for older/single-pose PDBQTs with no MODEL records at all.
    """
    import re
    try:
        text = Path(pdbqt_path).read_text(errors="replace")
    except Exception:
        return 0
    n = len(re.findall(r"^MODEL\s+\d+", text, flags=re.MULTILINE))
    if n > 0:
        return n
    return 1 if ("ATOM" in text or "HETATM" in text) else 0


def _is_hydrogen(line: str) -> bool:
    """
    Return True if a PDBQT/PDB ATOM/HETATM line represents a hydrogen.

    Uses ONLY the atom-name field (cols 12-15). The AD4/Vina type column
    (last token) is NOT used because obabel assigns type "H" to heavy atoms
    when bond perception fails on CONECT-less crystal PDB files, which would
    cause every heavy atom in the reference to be filtered out (ref=1 bug).

    Covers standard naming conventions:
      - Names starting with H:  H, HA, HB, HD, HS, HZ ...
      - Digit-prefixed H names: 1H, 2HB, 3HG ... (old PDB convention)
    """
    if len(line) < 17:
        return False
    atom_name = line[12:16].strip().upper()
    return atom_name.startswith("H") or (
        len(atom_name) >= 2 and atom_name[0].isdigit() and atom_name[1] == "H"
    )


def _coordinate_bridge_rmsd(ref_heavy_lines: list, pose_heavy_lines: list,
                             lig_stem: str, prep_lig_dir,
                             ref_sdf_text: str = "",
                             tolerance: float = 0.02,
                             max_automorphisms: int = 2000):
    """
    Best-case RMSD path: establishes the reference<->docked-pose atom
    correspondence via EXACT coordinate identity through the ligand
    preparation pipeline's own provenance, instead of via bond-graph
    matching or spatial nearest-neighbor guessing.

    Why this can succeed where compute_rmsd_from_sdf()'s graph-based
    stages cannot: AutoDock Vina never reorders, adds, or removes atoms
    in a docked pose PDBQT relative to the PRE-DOCKING prepared ligand
    PDBQT it was given -- it only rotates/translates them. And
    MGLTools' torsion-tree reordering from the (optionally protonated)
    input structure into that prepared PDBQT is itself just a
    reordering of the SAME 3D coordinates, not a recomputation of them.
    So whenever the uploaded reference has the same coordinates as
    whatever originally went into ligand prep for this ligand (true
    whenever the reference is the same co-crystallised structure used
    to build raw_ligands/{lig_stem}, not a separately regenerated one),
    we can chain:
        reference atom -[exact coordinate match]-> prepared PDBQT row
                        -[same row, Vina invariant]-> docked pose row
    with no bond-perception step, and no assumption that the pose is
    spatially close to the reference. That second point is what makes
    this succeed on cases the graph/spatial fallbacks in
    compute_rmsd_from_sdf() cannot: a rotatable substituent that moved
    far from its reference position during docking still has a fixed,
    correct correspondence here, because the correspondence was never
    spatial to begin with -- it was fixed by ligand-prep bookkeeping
    before docking happened.

    SYMMETRY REFINEMENT (second layer): the coordinate bridge fixes
    atom IDENTITY (this is atom 5 in both files), but says nothing
    about molecular SYMMETRY -- if the ligand has a locally symmetric
    group (a phenyl ring, a carboxylate, a tert-butyl), that group can
    flip into a chemically-indistinguishable orientation during docking
    without the pose being "wrong". The bridge alone would score that
    flip as a large per-atom displacement, even though a symmetry-
    equivalent relabelling of the SAME correspondence gives a much
    smaller, chemically honest RMSD. This refinement computes the
    reference ligand's own automorphism group (from ref_sdf_text, a
    cleanly-bonded source -- NOT the pose's PDBQT-derived bonds, which
    is exactly the unreliable input the bridge was built to avoid) and
    minimizes RMSD over that group, on top of the bridge's fixed
    identity mapping. This never changes WHICH atoms correspond, only
    which chemically-equivalent labelling of that correspondence is
    scored. If ref_sdf_text is unavailable, doesn't parse, or its
    heavy-atom order doesn't line up element-for-element with
    ref_heavy_lines, this refinement is silently skipped and the plain
    bridge RMSD is returned unchanged -- it is a strict enhancement,
    never a requirement.

    Returns (rmsd, note) on success. Returns None if the bridge does
    not apply (no prepared-ligand PDBQT retained for this ligand, atom
    counts don't match, or some reference atom has no unambiguous exact
    coordinate match) -- callers should fall back to
    compute_rmsd_from_sdf() in that case; nothing here is a substitute
    for it, only a preferred first attempt.
    """
    prep_pdbqt = Path(prep_lig_dir) / f"{lig_stem}.pdbqt"
    if not prep_pdbqt.exists():
        return None

    _TWO_LETTER_ELEMENTS = {
        "CL", "BR", "NA", "MG", "CA", "ZN", "FE", "MN", "CU", "SE", "SI", "AL",
    }

    def _element_from_name(name: str) -> str:
        """
        Derive an element symbol from a PDB/PDBQT atom-name field (cols
        12-16). Used instead of trusting the trailing columns, because
        PDBQT replaces the standard PDB element column with charge +
        AutoDock atom type (e.g. "0.266  A", "-0.219 NA") -- reading
        that as an element is wrong, not just imprecise. Handles a
        leading digit (old-style "1HB" hydrogen naming) and a short
        list of common two-letter elements; anything else falls back
        to the first alphabetic character, which is exactly what
        _is_hydrogen() already relies on for H-detection.
        """
        s = name.strip().upper()
        if s and s[0].isdigit():
            s = s[1:]
        if len(s) >= 2 and s[:2] in _TWO_LETTER_ELEMENTS:
            return s[:2]
        return s[:1]

    def _parse_heavy(lines):
        out = []
        for line in lines:
            if not line.startswith(("ATOM", "HETATM")) or _is_hydrogen(line):
                continue
            elem = _element_from_name(line[12:16])
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            out.append((elem, x, y, z))
        return out

    try:
        prep_lines = prep_pdbqt.read_text(errors="replace").splitlines(keepends=True)
    except Exception:
        return None

    prep_heavy = _parse_heavy(prep_lines)
    ref_heavy = _parse_heavy(ref_heavy_lines)
    pose_heavy = _parse_heavy(pose_heavy_lines)

    if not prep_heavy or len(prep_heavy) != len(pose_heavy):
        # Not the pose that came from this prepared ligand -- filename
        # matching upstream should prevent this, but never silently
        # mismap if it happens anyway.
        return None

    if len(ref_heavy) != len(prep_heavy):
        return None

    # ref_row -> prep_row via exact (tolerance-gated) coordinate match.
    # Each prep atom is consumed at most once.
    used = set()
    ref_to_prep = {}
    for ri, (relem, rx, ry, rz) in enumerate(ref_heavy):
        best_j, best_d = None, None
        for j, (pelem, px, py, pz) in enumerate(prep_heavy):
            if j in used or pelem != relem:
                continue
            d = ((rx - px) ** 2 + (ry - py) ** 2 + (rz - pz) ** 2) ** 0.5
            if d <= tolerance and (best_d is None or d < best_d):
                best_j, best_d = j, d
        if best_j is None:
            return None
        used.add(best_j)
        ref_to_prep[ri] = best_j

    # Compose ref_row -> prep_row -> pose_row (same row index in the
    # docked pose as in the prepared ligand -- the Vina invariant).
    n = len(ref_to_prep)

    def _rmsd_for_mapping(idx_map):
        """RMSD for a given ref_row -> pose_row correspondence (a dict
        or list indexable by ref row)."""
        total = 0.0
        for ri in range(n):
            relem, rx, ry, rz = ref_heavy[ri]
            pose_row = idx_map[ri]
            pelem, px, py, pz = pose_heavy[pose_row]
            if pelem != relem:
                return None  # sanity check: element must still agree
            total += (rx - px) ** 2 + (ry - py) ** 2 + (rz - pz) ** 2
        return round((total / n) ** 0.5, 4)

    base_rmsd = _rmsd_for_mapping(ref_to_prep)
    if base_rmsd is None:
        return None  # sanity check: element must still agree at this row

    # -----------------------------------------------------------------
    # Symmetry refinement (optional second layer, see docstring). Never
    # changes WHICH atoms correspond -- only searches for a chemically
    # equivalent relabelling of that same correspondence that scores a
    # lower, more honest RMSD when the ligand has genuine local
    # symmetry (e.g. a flipped phenyl ring or carboxylate).
    # -----------------------------------------------------------------
    best_rmsd = base_rmsd
    symmetry_note = ""

    if ref_sdf_text and ref_sdf_text.strip():
        try:
            ref_mol = _Chem.MolFromMolBlock(ref_sdf_text, removeHs=False, sanitize=True)
            if ref_mol is not None:
                h_idx = [a.GetIdx() for a in ref_mol.GetAtoms() if a.GetAtomicNum() == 1]
                if h_idx:
                    rw = _Chem.RWMol(ref_mol)
                    for idx in sorted(h_idx, reverse=True):
                        rw.RemoveAtom(idx)  # descending order preserves relative order of the rest
                    ref_mol = rw.GetMol()

                # Sanity check: heavy-atom element sequence must line up
                # 1:1 with ref_heavy before any automorphism is trusted --
                # if Open Babel's SDF conversion reordered atoms relative
                # to the raw reference lines, automorphism indices would
                # silently refer to the wrong atoms.
                mol_elems = [a.GetSymbol().upper() for a in ref_mol.GetAtoms()]
                line_elems = [e.upper() for (e, _, _, _) in ref_heavy]
                if len(mol_elems) == len(line_elems) and mol_elems == line_elems:
                    matches = ref_mol.GetSubstructMatches(
                        ref_mol, uniquify=False, useChirality=True,
                        maxMatches=max_automorphisms,
                    )
                    n_automorphisms = len(matches)
                    for sigma in matches:
                        if len(sigma) != n:
                            continue
                        # sigma(i) is the ref atom equivalent to atom i;
                        # score that relabelling against the bridge's
                        # already-fixed identity mapping.
                        relabelled_map = [ref_to_prep[sigma[i]] for i in range(n)]
                        cand = _rmsd_for_mapping(relabelled_map)
                        if cand is not None and cand < best_rmsd:
                            best_rmsd = cand
                    if n_automorphisms > 1 and best_rmsd < base_rmsd:
                        symmetry_note = (
                            f" A {n_automorphisms}-fold local symmetry in the "
                            "reference was found and used to refine the score "
                            f"(unrefined bridge RMSD was {base_rmsd} \u00c5)."
                        )
        except Exception:
            pass  # symmetry refinement is optional -- fall through to base_rmsd

    return best_rmsd, (
        "Exact atom correspondence via ligand-preparation provenance "
        "(reference matched to the pre-docking prepared ligand PDBQT by "
        "coordinate identity; Vina preserves atom order from that file "
        "into every docked pose) \u2014 no bond-graph matching or spatial "
        "correspondence guessing was needed or used." + symmetry_note
    )


def pdb_to_sdf_text(pdb_text: str, obabel_path: str) -> str:
    """
    Convert a reference ligand PDB (text) to SDF using Open Babel, the same
    way pdbqt_to_sdf_text() already does for docked poses. This is what
    makes redocking RMSD naming-agnostic: both molecules end up as SDF
    graphs with bond tables, and compute_rmsd_from_sdf() matches them by
    graph structure, not by atom name — which matters because ligands
    prepared from SDF/MOL2 (see the screening pipeline above) get
    Open-Babel-generated PDB-style atom names that will not, in general,
    match the RCSB chemical-component-dictionary names on a reference
    ligand downloaded straight from the PDB.
    """
    import tempfile
    tmp_pdb = tmp_sdf = None
    try:
        tmp_pdb = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode="w")
        tmp_pdb.write(pdb_text)
        tmp_pdb.flush()
        tmp_pdb.close()
        tmp_sdf = tempfile.NamedTemporaryFile(suffix=".sdf", delete=False)
        tmp_sdf.close()
        result = subprocess.run(
            [obabel_path, tmp_pdb.name, "-O", tmp_sdf.name],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
        )
        if result.returncode == 0 and Path(tmp_sdf.name).exists():
            sdf_content = Path(tmp_sdf.name).read_text(errors="replace")
            if sdf_content.strip():
                return sdf_content
        return ""
    except Exception:
        return ""
    finally:
        for tmp in (tmp_pdb, tmp_sdf):
            try:
                if tmp:
                    os.unlink(tmp.name)
            except Exception:
                pass


def extract_vina_pose_by_mode(pdbqt_path: str, mode_number: int) -> list:
    """Return ATOM/HETATM lines for a specific Vina mode."""
    lines = Path(pdbqt_path).read_text(errors="replace").splitlines(keepends=True)
    current_model, in_target, coords = 0, False, []
    for line in lines:
        if line.startswith("MODEL"):
            try:
                current_model = int(line.split()[1])
            except (IndexError, ValueError):
                current_model += 1
            in_target = (current_model == mode_number)
        elif line.startswith("ENDMDL"):
            if in_target:
                break
            in_target = False
        elif in_target and line.startswith(("ATOM", "HETATM")):
            coords.append(line)
    if not coords and mode_number == 1:
        coords = [l for l in lines if l.startswith(("ATOM", "HETATM"))]
    return coords


def compute_rmsd_from_lines(ref_lines: list, pose_lines: list):
    """
    Fast pre-check ONLY: verifies the reference and docked pose plausibly
    have the same number of heavy atoms before we bother invoking Open Babel
    + RDKit. This is NOT the RMSD calculation itself (that's
    compute_rmsd_from_sdf, below) — do not use this function's absence of an
    atom-count error to mean the RMSD is valid.

    Returns (n_ref_heavy, n_pose_heavy, error_string). error_string is None
    if the counts are compatible.
    """
    def _count_heavy(lines):
        n = 0
        for line in lines:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if _is_hydrogen(line):
                continue
            n += 1
        return n

    n_ref, n_pose = _count_heavy(ref_lines), _count_heavy(pose_lines)
    if n_ref == 0:
        return n_ref, n_pose, "Reference has no heavy atoms"
    if n_pose == 0:
        return n_ref, n_pose, "Pose has no heavy atoms"
    if n_ref != n_pose:
        return n_ref, n_pose, f"Atom count mismatch: ref={n_ref}, pose={n_pose}"
    return n_ref, n_pose, None


def compute_rmsd_from_sdf(ref_sdf_text: str, pose_sdf_text: str):
    """
    Two-stage redocking RMSD validation.

    Stage 1: sPyRMSD symmetry-aware graph matching on the reference and
             docked pose exactly as independently parsed.
    Stage 2: a cheap RDKit sanitize/normalize pass on both molecules,
             then sPyRMSD again — this fixes cases where the raw parse
             just wasn't sanitized correctly (valence/aromaticity
             flags), not genuine connectivity differences.
    No RMSD is returned if neither stage produces a valid graph-
    isomorphism correspondence. Cases where bond perception itself
    disagrees between reference and pose (not just sanitization) are
    reported as failures rather than resolved via spatial atom
    correspondence, since a nearest-atom mapping doesn't guarantee the
    atoms are chemically the same one — only graph-verified
    correspondences are used to compute RMSD.

    RMSD is heavy-atom RMSD in the original receptor coordinate frame.
    No translation/rotation minimization is applied.
    """
    if not _RMSD_DEPS_OK:
        return None, (
            "Validation RMSD dependencies are unavailable. "
            f"Import error: {_RMSD_DEPS_ERROR}"
        ), None

    def _strip_h(mol):
        h_idx = [
            a.GetIdx() for a in mol.GetAtoms()
            if a.GetAtomicNum() == 1
        ]
        if not h_idx:
            return mol
        rw = _Chem.RWMol(mol)
        for idx in sorted(h_idx, reverse=True):
            rw.RemoveAtom(idx)
        return rw.GetMol()

    def _extract_pdbqt_model(text, model_number=1):
        """
        Extract a single MODEL from a multi-pose PDBQT file.

        AutoDock Vina writes multiple docking poses as MODEL blocks.
        Validation must compare the reference ligand against one
        individual pose, not the concatenated coordinates of all poses.
        """
        if not text:
            return text

        lines = text.splitlines()

        has_models = any(
            line.startswith("MODEL")
            for line in lines
        )

        if not has_models:
            return text

        selected = []
        inside = False

        for line in lines:
            if line.startswith("MODEL"):
                parts = line.split()

                try:
                    number = int(parts[1])
                except Exception:
                    number = None

                inside = (number == model_number)

                if inside:
                    selected.append(line)

                continue

            if line.startswith("ENDMDL"):
                if inside:
                    selected.append(line)
                inside = False
                continue

            if inside:
                selected.append(line)

        if not selected:
            return None

        return "\n".join(selected) + "\n"


    def _parse(text, sanitize=True):
        # PDBQT docking output may contain multiple Vina poses.
        # Keep only MODEL 1 for validation.
        extracted = _extract_pdbqt_model(text, model_number=1)

        if extracted is None:
            return None

        # MolBlock parsing is retained for SDF/MolBlock input.
        mol = _Chem.MolFromMolBlock(
            extracted,
            removeHs=True,
            sanitize=sanitize
        )

        if mol is None and sanitize:
            mol = _Chem.MolFromMolBlock(
                text,
                removeHs=False,
                sanitize=False
            )
            if mol is not None:
                try:
                    _Chem.SanitizeMol(
                        mol,
                        sanitizeOps=(
                            _Chem.SANITIZE_ALL
                            ^ _Chem.SANITIZE_PROPERTIES
                        )
                    )
                except Exception:
                    pass

        if mol is None:
            return None

        return _strip_h(mol)

    def _spyrmsd(ref_mol, pose_mol):
        r = _spymol.Molecule.from_rdkit(ref_mol)
        q = _spymol.Molecule.from_rdkit(pose_mol)

        r.strip()
        q.strip()

        result = _spyrmsd_rmsd.symmrmsd(
            r.coordinates,
            q.coordinates,
            r.atomicnums,
            q.atomicnums,
            r.adjacency_matrix,
            q.adjacency_matrix,
            center=False,
            minimize=False,
        )

        value = (
            float(result[0])
            if isinstance(result, (list, tuple, np.ndarray))
            else float(result)
        )

        return round(value, 4)

    def _normalize_bonds(mol):
        """
        Rebuild a sanitized copy where possible. This is deliberately
        conservative: connectivity is preserved and no arbitrary bonds
        are created.
        """
        cp = _Chem.Mol(mol)

        try:
            _Chem.SanitizeMol(cp)
        except Exception:
            try:
                _Chem.SanitizeMol(
                    cp,
                    sanitizeOps=(
                        _Chem.SANITIZE_ALL
                        ^ _Chem.SANITIZE_PROPERTIES
                    )
                )
            except Exception:
                pass

        return _strip_h(cp)

    ref = _parse(ref_sdf_text)
    pose = _parse(pose_sdf_text)

    if ref is None:
        return None, "Could not parse reference ligand with RDKit.", None

    if pose is None:
        return None, "Could not parse docked pose with RDKit.", None

    n_ref = ref.GetNumAtoms()
    n_pose = pose.GetNumAtoms()

    if n_ref == 0:
        return None, "Reference ligand contains no heavy atoms.", None

    if n_pose == 0:
        return None, "Docked pose contains no heavy atoms.", None

    if n_ref != n_pose:
        return None, (
            f"Heavy-atom count mismatch: reference={n_ref}, "
            f"pose={n_pose}. Re-docking RMSD requires the same ligand."
        ), None

    debug = []  # human-readable per-stage diagnostics, surfaced on total failure

    # ---------------------------------------------------------------
    # STAGE 1 - primary sPyRMSD symmetry-aware graph matching
    # ---------------------------------------------------------------
    try:
        value = _spyrmsd(ref, pose)

        return (
            value,
            None,
            "sPyRMSD symmetry-corrected graph matching "
            "(primary validation)"
        )

    except Exception as exc1:
        debug.append(f"stage1 (sPyRMSD, raw): {exc1}")

    # ---------------------------------------------------------------
    # STAGE 2 - bond-graph repair + sPyRMSD
    # ---------------------------------------------------------------
    # RDKit sanitize/normalize (cheap, unchanged from before) — fixes
    # cases where the raw parse just wasn't sanitized correctly
    # (valence/aromaticity flags), not genuine connectivity
    # differences. If the reference's and the pose's independently
    # perceived bond graphs are still not isomorphic after this, that
    # is treated as a genuine failure rather than papered over with a
    # spatial (non-graph-verified) atom correspondence.
    # ---------------------------------------------------------------
    try:
        ref_norm = _normalize_bonds(ref)
        pose_norm = _normalize_bonds(pose)

        value = _spyrmsd(ref_norm, pose_norm)

        return (
            value,
            None,
            "sPyRMSD symmetry-corrected graph matching "
            "(after RDKit graph normalization)"
        )

    except Exception as exc2a:
        debug.append(f"stage2a (sPyRMSD, normalized): {exc2a}")

    return None, (
        "No atom correspondence could be established after sPyRMSD "
        "(raw) or RDKit graph normalization. RMSD was not "
        "calculated. Per-stage diagnostics: " + " | ".join(debug)
    ), "NO RMSD"


# ──────────────────────────────────────────────────────────────────────────────
# BINDING-SITE CONTACT ANALYSIS
# (geometric, distance-based residue/interaction mapping for the pose viewer
#  and for figure/caption generation -- same spirit as Schrodinger's ligand
#  interaction diagram, built without any external interaction-analysis
#  dependency. This is a HEURISTIC: distance-only, no donor/acceptor angle
#  geometry, no salt-bridge/pi-cation detection. Treat flagged interactions
#  as candidates for a figure, not as a validated substitute for a tool like
#  PLIP if the manuscript needs a rigorously defined interaction footprint.)
# ──────────────────────────────────────────────────────────────────────────────

# Side-chain atoms capable of acting as H-bond donors / acceptors, by residue.
# Backbone N (donor) and O (acceptor) apply to every residue and are added
# automatically in compute_ligand_contacts().
_AA_DONOR_ATOMS = {
    "SER": {"OG"}, "THR": {"OG1"}, "TYR": {"OH"}, "CYS": {"SG"},
    "HIS": {"ND1", "NE2"}, "TRP": {"NE1"}, "LYS": {"NZ"},
    "ARG": {"NE", "NH1", "NH2"}, "ASN": {"ND2"}, "GLN": {"NE2"},
}
_AA_ACCEPTOR_ATOMS = {
    "SER": {"OG"}, "THR": {"OG1"}, "TYR": {"OH"},
    "ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"},
    "ASN": {"OD1"}, "GLN": {"OE1"},
    "HIS": {"ND1", "NE2"}, "MET": {"SD"}, "CYS": {"SG"},
}
_HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "CYS"}
# Ring atoms used for aromatic-ring centroid (pi-stacking proxy). TRP's
# 6-membered ring only (indole benzo ring) is used for simplicity.
_AROMATIC_RING_ATOMS = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
}
# Side-chain atoms that carry a formal charge at physiological pH, used for
# salt-bridge detection against a formally-charged ligand atom (from RDKit;
# unavailable in the PDB-fallback path, which has no bond/charge info).
# HIS is included as a weak base -- it's often neutral, so a HIS "salt
# bridge" here is a softer call than LYS/ARG and worth a second look.
_AA_BASIC_ATOMS  = {"LYS": {"NZ"}, "ARG": {"NE", "NH1", "NH2"}, "HIS": {"ND1", "NE2"}}
_AA_ACIDIC_ATOMS = {"ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}}

_CONTACT_COLORS = {
    "H-bond": "#e63946",
    "salt-bridge": "#2a9d8f",
    "hydrophobic": "#f4a300",
    "pi-stacking": "#8338ec",
    "contact": "#adb5bd",
}


def _parse_receptor_pdb_residues(pdb_text: str) -> dict:
    """Parse ATOM records into {(chain, resnum, resname): [atoms]}."""
    residues: dict = {}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            atom_name = line[12:16].strip()
            resname   = line[17:20].strip()
            chain     = line[21].strip() or "A"
            resnum    = int(line[22:26])
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            element   = (line[76:78].strip() or atom_name[0]).strip().capitalize()
        except (ValueError, IndexError):
            continue
        if atom_name.startswith("H") or element == "H":
            continue
        residues.setdefault((chain, resnum, resname), []).append(
            {"name": atom_name, "element": element, "x": x, "y": y, "z": z}
        )
    return residues


def _ligand_atoms_from_sdf(sdf_text: str):
    """Return (atoms, aromatic_ring_centroids) parsed via RDKit, or None on
    failure. RDKit gives correct bond orders/aromaticity, so donor/acceptor
    and ring perception is far more reliable than a coordinate-only guess."""
    if not _RMSD_DEPS_OK or not sdf_text.strip():
        return None
    try:
        mol = _Chem.MolFromMolBlock(sdf_text, sanitize=True, removeHs=False)
        if mol is None:
            mol = _Chem.MolFromMolBlock(sdf_text, sanitize=False, removeHs=False)
        if mol is None:
            return None
        conf = mol.GetConformer()
        atoms = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "H":
                continue
            idx = atom.GetIdx()
            pos = conf.GetAtomPosition(idx)
            sym = atom.GetSymbol()
            try:
                n_h = atom.GetTotalNumHs()
            except Exception:
                n_h = 0
            try:
                formal_charge = atom.GetFormalCharge()
            except Exception:
                formal_charge = 0
            atoms.append({
                "symbol": sym, "x": pos.x, "y": pos.y, "z": pos.z,
                "is_donor": sym in ("N", "O") and n_h > 0,
                "is_acceptor": sym in ("N", "O", "F"),
                "formal_charge": formal_charge,
            })
        ring_centroids = []
        try:
            ri = mol.GetRingInfo()
            for ring in ri.AtomRings():
                if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                    pts = np.array([[conf.GetAtomPosition(i).x,
                                      conf.GetAtomPosition(i).y,
                                      conf.GetAtomPosition(i).z] for i in ring])
                    ring_centroids.append(pts.mean(axis=0))
        except Exception:
            pass
        return atoms, ring_centroids
    except Exception:
        return None


def _ligand_atoms_from_pdb_fallback(pdb_text: str):
    """Coordinate-only fallback when SDF/RDKit isn't available. Donor/acceptor
    flags are a crude element-based guess (no bond info to check H count)."""
    atoms = []
    for line in pdb_text.splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        try:
            atom_name = line[12:16].strip()
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            element = (line[76:78].strip() or atom_name[0]).strip().capitalize()
        except (ValueError, IndexError):
            continue
        if element == "H":
            continue
        atoms.append({
            "symbol": element, "x": x, "y": y, "z": z,
            "is_donor": element in ("N", "O"),
            "is_acceptor": element in ("N", "O"),
        })
    return atoms


def compute_ligand_contacts(receptor_pdb_text: str, ligand_text: str, lig_format: str,
                             pocket_cutoff: float = 4.5, hbond_cutoff: float = 3.5,
                             backbone_hbond_cutoff: float = 3.0,
                             pi_cutoff: float = 5.5,
                             salt_bridge_cutoff: float = 4.0) -> list:
    """Distance-based binding-site contact analysis. Returns one dict per
    contact residue, sorted by closest approach:

        {chain, resnum, resname, residue, min_dist, interaction, detail,
         centroid, hbond_pair}

    interaction is one of "salt-bridge", "H-bond", "pi-stacking",
    "hydrophobic", "contact", chosen by whichever category has the closest
    qualifying atom pair (salt-bridge wins ties against a plain H-bond,
    since a charge-complementary polar contact is the more specific label).

    Backbone N/O is present in EVERY residue, so it is by far the biggest
    source of false-positive "H-bond" calls once a generic distance cutoff
    is applied across an entire pocket (a backbone amide that merely
    happens to pass within range of some ligand heteroatom, with no check
    on whether it's actually oriented to donate/accept, looks identical to
    a real one). Side-chain donor/acceptor atoms (SER-OG, ASN-ND2, etc.)
    are a much more specific signal -- a residue's side chain pointing at
    the ligand is comparatively rare -- so they keep the full hbond_cutoff
    while bare backbone N/O is held to a tighter backbone_hbond_cutoff.
    This does not add real donor-H...acceptor angle geometry (PDB files
    generally lack explicit backbone H coordinates), it just stops treating
    "any backbone atom within 3.5 A" as equally strong evidence as a
    genuine side-chain H-bond.

    Salt-bridge detection needs RDKit-derived formal charges on the ligand
    (the SDF path); the PDB-fallback path has no bond/charge info, so this
    category will simply never fire there. HIS is included as a basic
    residue but is often neutral at physiological pH -- treat a HIS
    salt-bridge call as a softer, worth-checking-manually result.

    hbond_pair (if not None) is {"res":[x,y,z], "lig":[x,y,z]} for drawing
    a dashed line in the viewer -- used for both "H-bond" and "salt-bridge".
    """
    residues = _parse_receptor_pdb_residues(receptor_pdb_text)
    if not residues:
        return []

    parsed = _ligand_atoms_from_sdf(ligand_text) if lig_format == "sdf" else None
    if parsed is not None:
        lig_atoms, lig_ring_centroids = parsed
    else:
        lig_atoms, lig_ring_centroids = _ligand_atoms_from_pdb_fallback(ligand_text), []

    if not lig_atoms:
        return []

    lig_xyz = np.array([[a["x"], a["y"], a["z"]] for a in lig_atoms])

    contacts = []
    for (chain, resnum, resname), atoms in residues.items():
        res_xyz = np.array([[a["x"], a["y"], a["z"]] for a in atoms])
        dmat = np.linalg.norm(res_xyz[:, None, :] - lig_xyz[None, :, :], axis=-1)
        min_dist = float(dmat.min())
        if min_dist > pocket_cutoff:
            continue
        i_res, i_lig = np.unravel_index(np.argmin(dmat), dmat.shape)

        # Find the single BEST (minimum-distance) qualifying pair in each
        # candidate category, then classify by whichever category is
        # actually closest -- not by which category was checked first.
        # (Checking H-bond first and breaking on the first cutoff-satisfying
        # pair let a distant, marginal backbone contact outrank an obviously
        # dominant close hydrophobic contact; this compares real distances.)
        best = {"salt-bridge": (np.inf, None, None), "H-bond": (np.inf, None, None),
                "pi-stacking": (np.inf, None, None), "hydrophobic": (np.inf, None, None)}
        # Tie-break order when two categories land on the identical closest
        # pair (a charge-complementary polar contact satisfies both the
        # generic H-bond search and the salt-bridge search at the same
        # distance) -- prefer the more chemically specific label.
        _priority = {"salt-bridge": 0, "H-bond": 1, "pi-stacking": 2, "hydrophobic": 3}

        donor_names    = {"N"} | _AA_DONOR_ATOMS.get(resname, set())
        acceptor_names = {"O"} | _AA_ACCEPTOR_ATOMS.get(resname, set())
        polar_names    = donor_names | acceptor_names
        sidechain_polar_names = polar_names - {"N", "O"}  # side-chain-specific atoms only
        for ai, ra in enumerate(atoms):
            if ra["name"] not in polar_names:
                continue
            # Backbone N/O gets the tighter cutoff; side-chain donor/acceptor
            # atoms (more specific, less prone to coincidental proximity)
            # keep the full hbond_cutoff.
            this_cutoff = hbond_cutoff if ra["name"] in sidechain_polar_names else backbone_hbond_cutoff
            for aj, la in enumerate(lig_atoms):
                if not (la.get("is_donor") or la.get("is_acceptor")):
                    continue
                d = float(np.linalg.norm(res_xyz[ai] - lig_xyz[aj]))
                if d <= this_cutoff and d < best["H-bond"][0]:
                    best["H-bond"] = (d, ai, aj)

        # --- Salt bridge: basic residue atom + anionic ligand atom, or
        # acidic residue atom + cationic ligand atom. Requires RDKit-derived
        # formal charges on the ligand (SDF path); the PDB-fallback path has
        # no charge info, so every ligand atom reads formal_charge=0 and this
        # search simply finds nothing there -- same limitation as elsewhere.
        basic_names, acidic_names = _AA_BASIC_ATOMS.get(resname), _AA_ACIDIC_ATOMS.get(resname)
        if basic_names or acidic_names:
            for ai, ra in enumerate(atoms):
                if basic_names and ra["name"] in basic_names:
                    want_negative_lig = True
                elif acidic_names and ra["name"] in acidic_names:
                    want_negative_lig = False
                else:
                    continue
                for aj, la in enumerate(lig_atoms):
                    fc = la.get("formal_charge", 0)
                    if want_negative_lig and fc >= 0:
                        continue
                    if not want_negative_lig and fc <= 0:
                        continue
                    d = float(np.linalg.norm(res_xyz[ai] - lig_xyz[aj]))
                    if d <= salt_bridge_cutoff and d < best["salt-bridge"][0]:
                        best["salt-bridge"] = (d, ai, aj)

        if resname in _AROMATIC_RING_ATOMS and lig_ring_centroids:
            ring_names = _AROMATIC_RING_ATOMS[resname]
            ring_pts = np.array([[a["x"], a["y"], a["z"]] for a in atoms if a["name"] in ring_names])
            if len(ring_pts) >= 3:
                res_centroid = ring_pts.mean(axis=0)
                for lc in lig_ring_centroids:
                    d = float(np.linalg.norm(res_centroid - lc))
                    if d <= pi_cutoff and d < best["pi-stacking"][0]:
                        best["pi-stacking"] = (d, res_centroid, lc)

        if resname in _HYDROPHOBIC_RES:
            for ai, ra in enumerate(atoms):
                if ra["element"] != "C":
                    continue
                for aj, la in enumerate(lig_atoms):
                    if la["symbol"] != "C":
                        continue
                    d = float(np.linalg.norm(res_xyz[ai] - lig_xyz[aj]))
                    if d <= 4.5 and d < best["hydrophobic"][0]:
                        best["hydrophobic"] = (d, ai, aj)

        # Pick the category with the smallest qualifying distance; ties
        # (identical closest pair found by two searches) break by _priority.
        interaction = "contact"
        detail = f"{atoms[i_res]['name']} \u2194 nearest ligand atom ({min_dist:.2f} \u00c5)"
        hbond_pair = None
        best_cat, (best_d, ref1, ref2) = min(best.items(), key=lambda kv: (kv[1][0], _priority[kv[0]]))
        if best_d < np.inf:
            interaction = best_cat
            min_dist = min(min_dist, best_d)
            if best_cat in ("H-bond", "salt-bridge"):
                tag = "ionic pair" if best_cat == "salt-bridge" else "lig. atom"
                detail = f"{atoms[ref1]['name']} \u2194 {tag} ({best_d:.2f} \u00c5)"
                hbond_pair = {"res": res_xyz[ref1].tolist(), "lig": lig_xyz[ref2].tolist()}
            elif best_cat == "pi-stacking":
                detail = f"ring centroid \u2194 ring centroid ({best_d:.2f} \u00c5)"
            elif best_cat == "hydrophobic":
                detail = f"{atoms[ref1]['name']} \u2194 lig. C ({best_d:.2f} \u00c5)"

        contacts.append({
            "chain": chain, "resnum": resnum, "resname": resname,
            "residue": f"{resname}{resnum}",
            "min_dist": round(min_dist, 2),
            "interaction": interaction,
            "detail": detail,
            "centroid": res_xyz.mean(axis=0).tolist(),
            "hbond_pair": hbond_pair,
        })

    contacts.sort(key=lambda c: c["min_dist"])
    return contacts


def filter_pdb_to_residues(pdb_text: str, keys: set) -> str:
    """Return only the ATOM lines of pdb_text whose (chain, resnum) is in
    keys, preserving original formatting (used to build the 3D highlight
    layer without reconstructing PDB records by hand)."""
    lines = []
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            chain  = line[21].strip() or "A"
            resnum = int(line[22:26])
        except (ValueError, IndexError):
            continue
        if (chain, resnum) in keys:
            lines.append(line)
    return "\n".join(lines)


# ──────────────────────────────────────────────────────────────────────────────
# 3D VIEWER
# ──────────────────────────────────────────────────────────────────────────────

def build_viewer_html(receptor_pdb: str, ligand_pdb: str,
                      protein_name: str, ligand_name: str,
                      affinity: float | None,
                      lig_format: str = "pdb",
                      contacts: list | None = None,
                      contact_pdb: str = "") -> str:
    """Build a self-contained py3Dmol HTML viewer.

    lig_format should be 'sdf' when ligand string was converted via Open Babel
    (preserves bond orders so py3Dmol renders the correct structure).
    Falls back to 'pdb' if SDF conversion was unavailable.

    contacts / contact_pdb (optional): output of compute_ligand_contacts() /
    filter_pdb_to_residues(). When provided, adds a third model containing
    just the binding-site residues, colour-coded and labelled by interaction
    type (H-bond / hydrophobic / pi-stacking / contact), with dashed lines
    for H-bond pairs -- a Schrodinger-style interaction overlay on the pose.
    """
    aff_str  = f"{affinity:.2f} kcal/mol" if affinity is not None else "N/A"
    rec_js      = receptor_pdb.replace("\\","\\\\").replace("'","\\'").replace("\n","\\n")
    lig_js      = ligand_pdb.replace("\\","\\\\").replace("'","\\'").replace("\n","\\n")
    lig_fmt_js  = lig_format  # 'sdf' or 'pdb'

    has_contacts = bool(contacts) and bool(contact_pdb.strip())
    contact_js   = contact_pdb.replace("\\", "\\\\").replace("'", "\\'").replace("\n", "\\n")
    # Strip the raw 'atoms' payload isn't present here (compute_ligand_contacts
    # already omits it); dump only what the viewer needs.
    contact_json = json.dumps([
        {k: c[k] for k in ("chain", "resnum", "resname", "residue",
                            "min_dist", "interaction", "detail",
                            "centroid", "hbond_pair")}
        for c in (contacts or [])
    ])
    contact_colors_json = json.dumps(_CONTACT_COLORS)

    return f"""<!DOCTYPE html><html>
<head>
<meta charset="utf-8">
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono&family=Inter:wght@400;500;600&display=swap" rel="stylesheet">
<script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.6.4/jquery.min.js"></script>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
  * {{ box-sizing:border-box; margin:0; padding:0; }}
  body {{ background:#f7f7f5; font-family:'Inter',sans-serif; color:#1a1a1a; }}
  #viewer {{ width:100%; height:440px; position:relative; border:1px solid #ddd; border-bottom:none; overflow:hidden; }}
  #controls {{
    padding:8px 14px; background:#fafafa;
    border:1px solid #ddd; border-top:none;
    display:flex; gap:10px; flex-wrap:wrap; align-items:center; font-size:12px;
  }}
  #info {{
    padding:7px 14px; background:#f3f3f1;
    border:1px solid #ddd; border-top:none;
    font-size:11px; color:#555;
    display:flex; align-items:center; gap:8px; flex-wrap:wrap;
  }}
  #info b {{ color:#1a4a7a; }}
  .ctrl-label {{ color:#666; font-size:10px; font-weight:600; text-transform:uppercase; letter-spacing:0.06em; white-space:nowrap; }}
  select {{
    background:#fff; color:#1a1a1a;
    border:1px solid #ccc; border-radius:2px;
    padding:3px 7px; font-size:11.5px; font-family:'Inter',sans-serif; cursor:pointer;
  }}
  select:focus {{ outline:none; border-color:#1a4a7a; }}
  button {{
    background:#fff; color:#1a4a7a;
    border:1px solid #ccc; border-radius:2px;
    padding:4px 11px; cursor:pointer; font-size:11.5px; font-weight:500;
    transition:all 0.12s;
  }}
  button:hover {{ background:#e8eef6; border-color:#1a4a7a; }}
  .badge {{
    margin-left:auto; font-family:'IBM Plex Mono',monospace; font-size:10.5px;
    color:#444; border:1px solid #ccc; padding:2px 9px; border-radius:2px;
    background:#fff; white-space:nowrap;
  }}
</style>
</head>
<body>
<div id="viewer"></div>
<div id="controls">
  <div class="ctrl-label">Receptor</div>
  <select id="rec_style">
    <option value="cartoon">Cartoon</option>
    <option value="surface">Surface</option>
    <option value="line">Lines</option>
    <option value="stick">Sticks</option>
    <option value="sphere">Spheres</option>
    <option value="hidden">Hidden</option>
  </select>
  <div class="ctrl-label" style="margin-left:6px">Ligand</div>
  <select id="lig_style">
    <option value="stick">Sticks</option>
    <option value="sphere">Spheres</option>
    <option value="ball_stick">Ball+Stick</option>
    <option value="line">Lines</option>
  </select>
  <div class="ctrl-label" style="margin-left:6px">Colour</div>
  <select id="rec_color">
    <option value="spectrum">Spectrum</option>
    <option value="chain">By Chain</option>
    <option value="ss">Secondary Structure</option>
    <option value="white">White</option>
  </select>
  <div class="ctrl-label" style="margin-left:6px">Background</div>
  <select id="viewer_bg">
    <option value="#1a1f2e">Dark</option>
    <option value="#000000">Black</option>
    <option value="#ffffff">White</option>
    <option value="#f2f2f2">Light Grey</option>
    <option value="#808080">Grey</option>
  </select>
  <button onclick="viewer.zoomTo();">Reset View</button>
  <button onclick="viewer.spin(!spinning); spinning=!spinning; this.textContent=spinning?'Stop':'Spin';">Spin</button>
  <button onclick="downloadPosePNG();">Download PNG</button>
  {'''<label style="display:flex;align-items:center;gap:4px;font-size:11px;color:#666;margin-left:6px;cursor:pointer;">
    <input type="checkbox" id="show_contacts" checked style="cursor:pointer;">
    Binding-site residues
  </label>''' if has_contacts else ''}
  <span class="badge">AutoDock Vina &nbsp;&middot;&nbsp; {aff_str}</span>
</div>
<div id="info">
  <b>Protein:</b> {protein_name} &nbsp;&middot;&nbsp; <b>Ligand:</b> {ligand_name}
  &nbsp;&middot;&nbsp; <b>Best affinity:</b> {aff_str}
  &nbsp;&middot;&nbsp; <span style="color:#999;font-size:10px;">Drag to rotate &middot; Scroll to zoom &middot; Right-drag to translate</span>
  {f'''&nbsp;&middot;&nbsp;
  <span style="color:#e63946;">&#9632;</span> H-bond
  <span style="color:#2a9d8f;">&#9632;</span> Salt bridge
  <span style="color:#f4a300;">&#9632;</span> Hydrophobic
  <span style="color:#8338ec;">&#9632;</span> &pi;-stacking
  <span style="color:#adb5bd;">&#9632;</span> Contact''' if has_contacts else ''}
</div>
<script>
var viewer = $3Dmol.createViewer("viewer", {{backgroundColor:"#1a1f2e"}});
var spinning = false;
viewer.addModel('{rec_js}', 'pdb');
viewer.addModel('{lig_js}', '{lig_fmt_js}');
var hasContacts = {str(has_contacts).lower()};
if (hasContacts) {{
  viewer.addModel('{contact_js}', 'pdb');
}}
var models = viewer.getModelList();
var recModel = models[0];
var ligModel = models[1];

var contactData   = {contact_json};
var contactColors = {contact_colors_json};

function applyContactOverlay(show) {{
  if (!hasContacts) return;
  viewer.setStyle({{model: 2}}, {{}});
  viewer.removeAllLabels();
  if (show) {{
    contactData.forEach(function(c) {{
      var color = contactColors[c.interaction] || contactColors["contact"];
      viewer.setStyle({{model: 2, chain: c.chain, resi: c.resnum}},
                       {{stick: {{color: color, radius: 0.11}}}});
      viewer.addLabel(c.residue, {{
        position: {{x: c.centroid[0], y: c.centroid[1], z: c.centroid[2]}},
        backgroundColor: color, backgroundOpacity: 0.75,
        fontColor: "white", fontSize: 10, borderThickness: 0, inFront: true,
      }});
      if (c.hbond_pair) {{
        viewer.addLine({{
          start: {{x: c.hbond_pair.res[0], y: c.hbond_pair.res[1], z: c.hbond_pair.res[2]}},
          end:   {{x: c.hbond_pair.lig[0], y: c.hbond_pair.lig[1], z: c.hbond_pair.lig[2]}},
          color: color, dashed: true, linewidth: 2,
        }});
      }}
    }});
  }}
  viewer.render();
}}
if (hasContacts) {{
  var contactToggle = document.getElementById("show_contacts");
  if (contactToggle) contactToggle.onchange = function() {{ applyContactOverlay(this.checked); }};
}}

function downloadPosePNG() {{
  try {{
    var uri = viewer.pngURI();
    var a = document.createElement("a");
    a.href = uri;
    a.download = "IBDock_pose_view.png";
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
  }} catch (e) {{
    alert("PNG export failed: " + e);
  }}
}}

document.getElementById("viewer_bg").onchange = function() {{
  viewer.setBackgroundColor(this.value);
  viewer.render();
}};

function applyStyles() {{
  var rs = document.getElementById("rec_style").value;
  var ls = document.getElementById("lig_style").value;
  var rc = document.getElementById("rec_color").value;
  viewer.setStyle({{model: 0}}, {{}});
  viewer.setStyle({{model: 1}}, {{}});
  var colorScheme = rc === "spectrum" ? "spectrum" : rc === "chain" ? "chain" : rc === "ss" ? "ssPyMol" : {{color:"white"}};
  if (rs !== "hidden") {{
    var recOpts = {{}};
    if (rs === "cartoon") recOpts = {{cartoon: {{color: colorScheme}}}};
    else if (rs === "surface") {{ viewer.addSurface($3Dmol.SurfaceType.VDW, {{opacity:0.6, color:"white"}}, {{model:0}}); recOpts = {{cartoon:{{color:colorScheme,opacity:0.3}}}}; }}
    else if (rs === "line") recOpts = {{line: {{colorscheme: colorScheme}}}};
    else if (rs === "stick") recOpts = {{stick: {{colorscheme: colorScheme}}}};
    else if (rs === "sphere") recOpts = {{sphere: {{colorscheme: colorScheme, scale:0.4}}}};
    viewer.setStyle({{model: 0}}, recOpts);
  }}
  if (ls === "stick") viewer.setStyle({{model: 1}}, {{stick: {{colorscheme:"greenCarbon"}}}});
  else if (ls === "sphere") viewer.setStyle({{model: 1}}, {{sphere: {{colorscheme:"greenCarbon"}}}});
  else if (ls === "ball_stick") viewer.setStyle({{model: 1}}, {{stick: {{colorscheme:"greenCarbon"}}, sphere: {{colorscheme:"greenCarbon", scale:0.3}}}});
  else if (ls === "line") viewer.setStyle({{model: 1}}, {{line: {{colorscheme:"greenCarbon"}}}});
  viewer.render();
}}
document.getElementById("rec_style").onchange = applyStyles;
document.getElementById("lig_style").onchange = applyStyles;
document.getElementById("rec_color").onchange = applyStyles;
applyStyles();
applyContactOverlay(true);
viewer.zoomTo({{model: 1}});
viewer.zoom(0.9);
viewer.render();
</script>
</body></html>"""


# ──────────────────────────────────────────────────────────────────────────────
# APP HEADER
# ──────────────────────────────────────────────────────────────────────────────
st.markdown("""
<div style="
    display:flex; align-items:baseline; justify-content:space-between;
    padding: 0.8rem 0 0.9rem;
    border-bottom: 2px solid #1a1a1a;
    margin-bottom: 1.2rem;
">
    <div>
        <span style="font-family:'Source Serif 4',Georgia,serif; font-size:1.7rem; font-weight:600;
                    color:#1a1a1a; letter-spacing:-0.02em; line-height:1;">IBDock</span>
        <span style="font-family:'Inter',sans-serif; font-size:0.78rem;
                    color:#777; margin-left:14px; letter-spacing:0.02em;">
            Batch Molecular Docking &nbsp;&middot;&nbsp; AutoDock Vina &nbsp;&middot;&nbsp; Zero Command Line
        </span>
    </div>
    <span style="font-family:'IBM Plex Mono',monospace; font-size:0.7rem; color:#555;
                 border:1px solid #ccc; padding:2px 10px; border-radius:2px;">
        v1.0.0
    </span>
</div>
""", unsafe_allow_html=True)

# ──────────────────────────────────────────────────────────────────────────────
# SIDEBAR — CONFIGURATION
# ──────────────────────────────────────────────────────────────────────────────
st.sidebar.markdown("""
<div style="padding:0.6rem 0 0.5rem; border-bottom:1px solid #ccc; margin-bottom:0.7rem;">
    <div style="font-family:'Inter',sans-serif; font-size:0.78rem; font-weight:600;
                color:#333; letter-spacing:0.06em; text-transform:uppercase;">Configuration</div>
</div>
""", unsafe_allow_html=True)

st.sidebar.markdown("**Project Directory**")
_proj_input = st.sidebar.text_input(
    "Project path",
    value=st.session_state.get("_last_project_dir", ""),
    placeholder=r"e.g. D:\my_docking_project",
    key="_proj_input_field",
)
project_dir = Path(_proj_input.strip() or ".")
if _proj_input.strip():
    st.session_state["_last_project_dir"] = _proj_input.strip()

# Live project-dir feedback
if _proj_input.strip() and not project_dir.exists():
    st.sidebar.error("❌ Directory does not exist")

# Load per-project config (after we have project_dir candidate)
_cfg = _load_config(project_dir) if project_dir.exists() else dict(_DEF)

st.sidebar.markdown("---")
st.sidebar.markdown("**Tool Paths**")

mgl_python = st.sidebar.text_input(
    "MGLTools pythonsh / python.exe", _cfg["mgl_python"],
    help="Path to the MGLTools Python interpreter.",
)
prep_rec = st.sidebar.text_input("prepare_receptor4.py", _cfg["prep_rec"])
prep_lig = st.sidebar.text_input("prepare_ligand4.py",   _cfg["prep_lig"])
vina_path   = st.sidebar.text_input("AutoDock Vina executable", _cfg["vina_path"])
obabel_path = st.sidebar.text_input("Open Babel (obabel)",      _cfg["obabel_path"])

st.sidebar.markdown("---")
st.sidebar.caption("Portable build: bundled native tools are detected automatically from the application's `tools` folder. No external subsystem is required.")

st.sidebar.markdown("**Pocket Detection** *(optional)*")
fpocket_path = st.sidebar.text_input(
    "fpocket (bundled native tool)", _cfg["fpocket_path"],
    help="Used when pocket detection is selected in Protein Prep.",
)
p2rank_path = st.sidebar.text_input(
    "P2Rank (bundled native tool)", _cfg["p2rank_path"],
    help="P2Rank and fpocket use bundled native Windows tools.",
)

st.sidebar.markdown("---")

# Save config whenever paths change
if project_dir.exists():
    _save_config(project_dir, dict(
        mgl_python=mgl_python, prep_rec=prep_rec, prep_lig=prep_lig,
        vina_path=vina_path, obabel_path=obabel_path,
        fpocket_path=fpocket_path, p2rank_path=p2rank_path,
    ))

if st.sidebar.button("Validate All Tools", use_container_width=True):
    tools = [
        (mgl_python,   "MGLTools python.exe"),
        (prep_rec,     "prepare_receptor4.py"),
        (prep_lig,     "prepare_ligand4.py"),
        (vina_path,    "AutoDock Vina"),
        (obabel_path,  "Open Babel"),
        (fpocket_path, "fpocket"),
        (p2rank_path,  "P2Rank"),
    ]
    ph = st.sidebar.empty()
    ph.info("⏳ Validating…")
    results = [validate_tool(p, l) for p, l in tools]
    ph.empty()

    n_ok = sum(r["ok"] for r in results)
    if n_ok == len(results):
        st.sidebar.success(f"✅ All {n_ok} tools OK")
    else:
        st.sidebar.warning(f"⚠️ {n_ok}/{len(results)} tools OK")

    for r in results:
        icon = "✅" if r["ok"] else "❌"
        with st.sidebar.expander(f"{icon} {r['label']}", expanded=not r["ok"]):
            if r["ok"]:
                st.caption(f"**Status:** {r['detail']}")
                if r["version"]:
                    st.caption(f"**Version:** `{r['version']}`")
            else:
                st.error(r["detail"])



# ──────────────────────────────────────────────────────────────────────────────
# LANDING SCREEN  (shown when no valid project dir is set)
# ──────────────────────────────────────────────────────────────────────────────
if not project_dir.exists() or str(project_dir) == ".":
    st.markdown("""
<div style="max-width:720px; margin:2rem auto 0;">

  <div style="margin-bottom:2rem; border-bottom:1px solid #e0e0e0; padding-bottom:1.2rem;">
    <div style="font-family:'Source Serif 4',Georgia,serif; font-size:1.5rem; font-weight:600;
                color:#1a1a1a; letter-spacing:-0.02em;">Welcome to IBDock</div>
    <div style="font-family:'Inter',sans-serif; font-size:0.83rem; color:#666;
                margin-top:4px;">
        A graphical interface for batch virtual screening with AutoDock Vina.
        Set your project directory in the sidebar to begin.
    </div>
  </div>

  <div style="display:grid; grid-template-columns:1fr 1fr 1fr; gap:0.8rem; margin-bottom:1.8rem;">
    <div class="info-card">
      <h4>1 &mdash; Set Project Folder</h4>
      <p>Enter your project directory path in the sidebar. IBDock will create all subfolders automatically.</p>
    </div>
    <div class="info-card">
      <h4>2 &mdash; Add Input Files</h4>
      <p>Place <code>.pdb</code> proteins in <code>raw_proteins/</code> and ligands (<code>.sdf</code> / <code>.mol2</code>) in <code>raw_ligands/</code>.</p>
    </div>
    <div class="info-card">
      <h4>3 &mdash; Run the Pipeline</h4>
      <p>Work through the tabs in order: Protein Prep &rarr; Ligand Prep &rarr; Docking &rarr; Results.</p>
    </div>
  </div>

  <div style="background:#fafafa; border:1px solid #e0e0e0; border-radius:3px; padding:1.1rem 1.4rem; margin-bottom:1.2rem;">
    <div style="font-weight:600; color:#1a1a1a; font-size:0.85rem; margin-bottom:0.5rem; text-transform:uppercase; letter-spacing:0.05em; font-family:'Inter',sans-serif;">Expected folder structure</div>
    <code style="font-size:0.78rem; color:#333; line-height:1.9; display:block;">
      my_project/<br>
      &nbsp;&nbsp;&nbsp;&nbsp;raw_proteins/&nbsp;&nbsp;&nbsp;&nbsp;&larr; .pdb files<br>
      &nbsp;&nbsp;&nbsp;&nbsp;raw_ligands/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&larr; .sdf / .mol2 / .pdb files
    </code>
    <div style="margin-top:0.5rem; font-size:0.77rem; color:#777; font-family:'Inter',sans-serif;">All other folders (prepared_receptors/, docking/, results/, &hellip;) are created automatically.</div>
  </div>

</div>
    """, unsafe_allow_html=True)
    st.stop()

raw_prot = project_dir / "raw_proteins"
raw_lig  = project_dir / "raw_ligands"

# Do not globally stop the application when input folders are missing.
# Individual workflow tabs handle missing input folders themselves.
raw_prot_exists = raw_prot.exists()
raw_lig_exists = raw_lig.exists()

# Create output directories
grid_dir     = project_dir / "grid"
prep_rec_dir = project_dir / "prepared_receptors"
prep_lig_dir = project_dir / "prepared_ligands"
dock_dir     = project_dir / "docking"
result_dir   = project_dir / "results"
for d in [grid_dir, prep_rec_dir, prep_lig_dir, dock_dir, result_dir]:
    d.mkdir(exist_ok=True)

# ── File counts (used in tab labels + preflight) ───────────────────────────
_n_raw_prot  = len(list(raw_prot.glob("*.pdb")))
_n_raw_lig   = len([f for f in raw_lig.glob("*") if f.suffix.lower() in {".sdf", ".mol2", ".pdb"}])
_n_prep_rec  = len(list(prep_rec_dir.glob("*_receptor.pdbqt")))
_n_prep_lig  = len(list(prep_lig_dir.glob("*.pdbqt")))
_n_docked    = len(list(dock_dir.glob("*.pdbqt")))
_n_results   = len(list(result_dir.glob("*.txt")))


# ──────────────────────────────────────────────────────────────────────────────
# TABS
# ──────────────────────────────────────────────────────────────────────────────
# Tab labels — concise, no emoji clutter
_t1 = f"Protein Prep{f' ({_n_prep_rec})' if _n_prep_rec else f' ({_n_raw_prot} PDB)' if _n_raw_prot else ''}"
_t2 = f"Ligand Prep{f' ({_n_prep_lig})' if _n_prep_lig else f' ({_n_raw_lig})' if _n_raw_lig else ''}"
_t3 = f"Docking{f' ({_n_docked})' if _n_docked else ''}"
_t4 = f"Results{f' ({_n_results})' if _n_results else ''}"

tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs([
    _t1, _t2, _t3, _t4,
    "Pose Viewer",
    "Validation",
    "About & Cite",
])


# ══════════════════════════════════════════════════════════════════════════════
# TAB 1 — PROTEIN PREPARATION
# ══════════════════════════════════════════════════════════════════════════════
with tab1:
    st.markdown("### Protein Preparation")
    st.caption(
        "Cleans PDB files, auto-detects the binding site, computes a grid box, "
        "and converts to PDBQT format using MGLTools."
    )

    with st.expander(" Preparation settings", expanded=True):
        c1, c2 = st.columns(2)
        with c1:
            remove_waters      = st.checkbox("Remove water molecules",                True,  key="prot_waters")
            remove_non_protein = st.checkbox("Remove non-protein HETATM",            True,  key="prot_hetatm")
            remove_lone_pairs  = st.checkbox("Remove lone electron pairs",            True,  key="prot_lp")
            remove_nonpolar_h  = st.checkbox("Remove non-polar hydrogens",            False, key="prot_nph")
        with c2:
            calculate_charges  = st.checkbox("Calculate Gasteiger charges",           True,  key="prot_charges")
            preserve_altloc    = st.checkbox("Keep single altloc conformation",       True,  key="prot_altloc")
            preserve_metals    = st.checkbox("Preserve metal ions (Zn, Mg, Fe, Ca)", True,  key="prot_metals")

    with st.expander(" Grid box settings", expanded=True):
        st.markdown(" Grid centre detection method ")
        grid_centre_method = st.radio(
            "Grid centre detection method",
            [
                "Auto (recommended)",
                "Co-crystallised ligand",
                "P2Rank pocket prediction",
                "fpocket pocket prediction",
                "Blind docking (whole protein)",
            ],
            index=0,
            key="prot_grid_method",
            label_visibility="collapsed",
            help=(
                "**Auto**: tries co-crystallised ligand → P2Rank → fpocket → blind docking.\n\n"
                "**Co-crystallised ligand**: uses HETATM atoms in the PDB as the binding site centre.\n\n"
                "**P2Rank**: machine-learning pocket prediction using the bundled native Windows tool.\n\n"
                "**fpocket**: geometric pocket detection (bundled native Windows tool).\n\n"
                "**Blind docking**: centres on the whole protein — use for unknown binding sites."
            ),
        )
        st.markdown("")
        c1, c2 = st.columns(2)
        with c1:
            padding = st.slider(
                "Grid padding per face (Å)",
                min_value=2.0, max_value=16.0, value=5.0, step=0.5, key="prot_padding",
                help="Added to all geometry-based grids. Box = detected extent + 2× padding. "
                     "8 Å is recommended for virtual screening.",
            )
        with c2:
            protonation = st.radio(
                "Hydrogen treatment",
                ["Add all Hydrogens", "Add missing Hydrogens",
                 "Build bonds + add Hydrogens", "Maintain initial state"],
                index=1, key="prot_protonation",
            )
    st.markdown("")

    if st.button("▶ Prepare Proteins", type="primary", key="run_prot_prep"):
        proteins = list(raw_prot.glob("*.pdb"))
        if not proteins:
            st.error("❌ No .pdb files found in `raw_proteins/`")
            st.stop()

        prepared, failed = [], []
        progress = st.progress(0)
        timer    = st.empty()
        t0       = time.time()

        for i, pdb in enumerate(proteins, 1):
            name = pdb.stem
            # Update progress at START of each item so bar moves immediately
            progress.progress((i - 1) / len(proteins), text=f"Processing {name} ({i}/{len(proteins)})…")

            try:
                lines    = pdb.read_text().splitlines(keepends=True)
                lig_atoms, lig_res = extract_ligand_atoms(lines)

                # ── Grid centre detection ──────────────────────────────────
                _method = grid_centre_method  # user selection from UI
                size    = None  # initialise; set explicitly by pocket methods

                if _method == "Co-crystallised ligand":
                    if lig_atoms.size == 0:
                        st.warning(
                            f"⚠️ **{name}**: No co-crystallised ligand found in PDB. "
                            f"Falling back to blind docking."
                        )
                        atoms = extract_protein_atoms(lines)
                        if atoms.size == 0:
                            raise RuntimeError("No ATOM records found in PDB file")
                        center    = atoms.mean(axis=0)
                        grid_type = "blind (whole protein)"
                    else:
                        center    = lig_atoms.mean(axis=0)
                        grid_type = "co-crystallised ligand"

                elif _method == "P2Rank pocket prediction":
                    try:
                        center, size = run_p2rank(pdb, p2rank_path, padding)
                        grid_type    = "P2Rank pocket prediction"
                    except Exception as p2e:
                        st.warning(f"P2Rank failed for **{name}**: {p2e} — falling back to blind docking.")
                        atoms = extract_protein_atoms(lines)
                        if atoms.size == 0:
                            raise RuntimeError("No ATOM records found in PDB file")
                        center    = atoms.mean(axis=0)
                        size      = None  # will be computed in grid box section
                        grid_type = "blind (whole protein)"

                elif _method == "fpocket pocket prediction":
                    try:
                        pocket_dir  = run_fpocket(pdb, fpocket_path)
                        pocket_file = pocket_dir / "pockets" / "pocket1_atm.pdb"
                        center, size = get_fpocket_center_and_size(pocket_file, padding)
                        grid_type   = "fpocket pocket prediction"
                    except Exception as fpe:
                        st.warning(f"fpocket failed for **{name}**: {fpe} — falling back to blind docking.")
                        atoms = extract_protein_atoms(lines)
                        if atoms.size == 0:
                            raise RuntimeError("No ATOM records found in PDB file")
                        center    = atoms.mean(axis=0)
                        size      = None  # will be computed in grid box section
                        grid_type = "blind (whole protein)"

                elif _method == "Blind docking (whole protein)":
                    atoms = extract_protein_atoms(lines)
                    if atoms.size == 0:
                        raise RuntimeError("No ATOM records found in PDB file")
                    center    = atoms.mean(axis=0)
                    grid_type = "blind (whole protein)"

                else:  # "Auto (recommended)" — cascade
                    if lig_atoms.size > 0:
                        center    = lig_atoms.mean(axis=0)
                        size      = None  # will be computed in grid box section
                        grid_type = "co-crystallised ligand"
                    else:
                        try:
                            center, size = run_p2rank(pdb, p2rank_path, padding)
                            grid_type    = "P2Rank pocket prediction"
                        except Exception as p2e:
                            st.warning(f"P2Rank failed for {name}: {p2e} — trying fpocket…")
                            try:
                                pocket_dir   = run_fpocket(pdb, fpocket_path)
                                pocket_file  = pocket_dir / "pockets" / "pocket1_atm.pdb"
                                center, size = get_fpocket_center_and_size(pocket_file, padding)
                                grid_type    = "fpocket pocket prediction"
                            except Exception as fpe:
                                st.warning(f"fpocket failed for {name}: {fpe} — using blind docking.")
                                atoms = extract_protein_atoms(lines)
                                if atoms.size == 0:
                                    raise RuntimeError("No ATOM records found in PDB file")
                                center    = atoms.mean(axis=0)
                                size      = None  # will be computed in grid box section
                                grid_type = "blind (whole protein)"

                # ── Grid box size ──────────────────────────────────────────
                # P2Rank and fpocket methods already set `size` in their call
                # above (geometry-based).  For other methods, or fallbacks that
                # set size=None, we compute it here from real atom coordinates.
                if size is not None:
                    pass  # already set by pocket detection (real geometry)
                elif "ligand" in grid_type:
                    extent   = lig_atoms.max(axis=0) - lig_atoms.min(axis=0)
                    raw_size = extent + 2 * padding
                    size     = cap_grid([round_even(v) for v in raw_size])
                else:  # blind docking or any unhandled fallback
                    atoms       = extract_protein_atoms(lines)
                    _blind_pad  = max(2 * padding, 16.0)  # min 8 Å per face for blind docking
                    raw_size    = atoms.max(axis=0) - atoms.min(axis=0) + _blind_pad
                    size        = cap_grid([round_even(v) for v in raw_size])

                # ── Validate grid centre is near the protein ───────────────
                _prot_atoms = extract_protein_atoms(lines)
                if _prot_atoms.size > 0:
                    _min_dist = float(np.min(
                        np.linalg.norm(_prot_atoms - center, axis=1)
                    ))
                    if _min_dist > 10.0:
                        st.warning(
                            f"⚠️ **{name}**: grid centre is {_min_dist:.1f} Å from the "
                            f"nearest protein atom — it may be outside the binding site. "
                            f"Check your PDB for multiple HETATM residues or use the "
                            f"Custom Grid Box option in the Docking tab to set it manually."
                        )

                # ── Save grid config ───────────────────────────────────────
                prot_grid_dir = grid_dir / name
                prot_grid_dir.mkdir(exist_ok=True)
                cfg_file = prot_grid_dir / "grid_config.txt"
                cfg_file.write_text(
                    f"center_x = {center[0]:.3f}\n"
                    f"center_y = {center[1]:.3f}\n"
                    f"center_z = {center[2]:.3f}\n\n"
                    f"size_x = {size[0]}\n"
                    f"size_y = {size[1]}\n"
                    f"size_z = {size[2]}\n"
                )

                # ── Clean and convert to PDBQT ─────────────────────────────
                clean_pdb = prep_rec_dir / f"{name}.pdb"
                clean_pdb.write_text("".join(clean_protein(lines, lig_res)))

                remove_flags = []
                if remove_waters:    remove_flags.append("waters")
                if remove_nonpolar_h: remove_flags.append("nphs")
                if remove_lone_pairs: remove_flags.append("lps")
                if preserve_altloc:  remove_flags.append("altlocs")
                U_option = ",".join(remove_flags) or None

                A_option = None
                if protonation == "Add all Hydrogens":           A_option = "hydrogens"
                elif protonation == "Add missing Hydrogens":     A_option = "checkhydrogens"
                elif protonation == "Build bonds + add Hydrogens": A_option = "bonds_hydrogens"

                out_pdbqt = prep_rec_dir / f"{name}_receptor.pdbqt"
                cmd = [mgl_python, prep_rec, "-r", clean_pdb, "-o", out_pdbqt]
                if A_option:  cmd += ["-A", A_option]
                if U_option:  cmd += ["-U", U_option]
                if not calculate_charges: cmd += ["-C"]

                run_cmd(cmd)
                fix_pdbqt_atom_names(out_pdbqt)
                strip_receptor_hydrogens(out_pdbqt)

                st.success(
                    f"✅ **{name}** — Grid: {size[0]}×{size[1]}×{size[2]} Å  |  "
                    f"Centre: ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})  |  "
                    f"Method: {grid_type}"
                )
                prepared.append(name)

            except Exception as exc:
                failed.append(name)
                st.error(f"❌ **{name}** failed")
                with st.expander("Error details"):
                    st.code(str(exc))

            # Update to exact completion after each protein finishes
            progress.progress(i / len(proteins), text=f"Done: {name}")
            elapsed = time.time() - t0
            remaining = elapsed / i * (len(proteins) - i)
            timer.caption(f"⏱ Estimated remaining: {int(remaining)}s")

        progress.progress(1.0, text="Complete!")
        st.success("✅ Protein preparation complete")
        timer.empty()
        st.divider()
        c1, c2 = st.columns(2)
        c1.metric("✅ Prepared", len(prepared))
        c2.metric("❌ Failed",   len(failed))
        if failed:
            st.warning("Failed proteins: " + ", ".join(failed))


# ══════════════════════════════════════════════════════════════════════════════
# TAB 2 — LIGAND PREPARATION
# ══════════════════════════════════════════════════════════════════════════════
with tab2:
    st.markdown("### Ligand Preparation")
    st.caption("Converts SDF / MOL2 / PDB ligands to PDBQT using Open Babel and MGLTools.")

    with st.expander(" Chemistry settings", expanded=True):
        c1, c2 = st.columns(2)
        with c1:
            correct_ph     = st.checkbox("Correct protonation state",   True,  key="lig_ph")
            ph_value       = st.number_input("pH", 0.0, 14.0, 7.4, 0.1, key="lig_phval")
            randomize_pose = st.checkbox("Randomise input pose",        False, key="lig_rand")
            rigid_ligand   = st.checkbox("Treat as rigid (no torsions)", False, key="lig_rigid")
        with c2:
            generate_3d           = st.checkbox("Generate 3D coordinates",      True,  key="lig_3d")
            calc_lig_charges      = st.checkbox("Calculate Gasteiger charges",  True,  key="lig_charges")
            remove_nonpolar_h_lig = st.checkbox("Remove non-polar hydrogens",   True,  key="lig_nph")
            remove_lone_pairs_lig = st.checkbox("Remove lone pairs",            True,  key="lig_lp")
    st.markdown("")

    if st.button("▶ Prepare Ligands", type="primary", key="run_lig_prep"):
        ligands = [f for f in raw_lig.glob("*") if f.suffix.lower() in {".sdf", ".mol2", ".pdb"}]
        if not ligands:
            st.error("❌ No .sdf / .mol2 / .pdb files found in `raw_ligands/`")
            st.stop()

        prepared, failed = [], []
        progress = st.progress(0)
        timer    = st.empty()
        t0       = time.time()

        for i, lig_file in enumerate(ligands, 1):
            progress.progress((i - 1) / len(ligands), text=f"Processing {lig_file.name} ({i}/{len(ligands)})…")
            try:
                pdb_file  = prep_lig_dir / f"{lig_file.stem}.pdb"
                pdbqt_out = prep_lig_dir / f"{lig_file.stem}.pdbqt"

                # Step 1a: Open Babel — pH-based protonation FIRST, as its own
                # call. Combining -p with --gen3d in a single Open Babel
                # invocation silently discards the pH-based protonation
                # (verified directly: acetic acid stays protonated as -COOH
                # at pH 7.4 when -p and --gen3d are combined, regardless of
                # flag order, but correctly deprotonates to -COO- when -p is
                # run alone). Running -p as its own step, then generating 3D
                # coordinates on the ALREADY-protonated structure, preserves
                # the correct protonation state.
                if correct_ph:
                    pdb_protonated = prep_lig_dir / f"{lig_file.stem}_protonated.pdb"
                    cmd_ph = [str(obabel_path), str(lig_file), "-O", str(pdb_protonated),
                              "-p", str(ph_value)]
                    run_cmd(cmd_ph)
                    if not pdb_protonated.exists():
                        raise RuntimeError("Open Babel did not generate a protonated PDB file")
                    step1_input = pdb_protonated
                else:
                    step1_input = lig_file

                # Step 1b: Open Babel — 3D embedding, charges, etc. on top of
                # the (optionally) already-protonated structure.
                # NOTE: --addhydrogens is not a real Open Babel option (it is
                # silently ignored — verified empirically: atom count is
                # identical with or without it) and has been removed. Standard
                # hydrogens are added by --gen3d / -p above and by MGLTools'
                # -A checkhydrogens in Step 2 below.
                cmd = [str(obabel_path), str(step1_input), "-O", str(pdb_file)]
                if generate_3d:      cmd.append("--gen3d")
                if calc_lig_charges: cmd += ["--partialcharge", "gasteiger"]
                if randomize_pose:   cmd.append("--randomize")
                run_cmd(cmd)

                if not pdb_file.exists():
                    raise RuntimeError("Open Babel did not generate a PDB file")

                # Step 2: MGLTools → PDBQT
                rm_flags = []
                if remove_nonpolar_h_lig: rm_flags.append("nphs")
                if remove_lone_pairs_lig: rm_flags.append("lps")
                U_opt = ",".join(rm_flags)

                cmd = [str(mgl_python), str(prep_lig),
                       "-l", pdb_file.name, "-o", pdbqt_out.name,
                       "-A", "checkhydrogens"]
                if U_opt:        cmd += ["-U", U_opt]
                if rigid_ligand: cmd += ["-Z"]

                subprocess.run(
                    cmd, cwd=str(prep_lig_dir),
                    check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                )

                if not pdbqt_out.exists():
                    raise RuntimeError("MGLTools did not generate a PDBQT file")

                st.success(f"✅ **{lig_file.name}**")
                prepared.append(lig_file.name)

            except subprocess.CalledProcessError as exc:
                failed.append(lig_file.name)
                st.error(f"❌ **{lig_file.name}**")
                with st.expander("Error details"):
                    st.code(exc.stderr.decode() if exc.stderr else str(exc))
            except Exception as exc:
                failed.append(lig_file.name)
                st.error(f"❌ **{lig_file.name}**")
                with st.expander("Error details"):
                    st.code(str(exc))

            progress.progress(i / len(ligands), text=f"Done: {lig_file.name}")
            elapsed = time.time() - t0
            remaining = elapsed / i * (len(ligands) - i)
            timer.caption(f"⏱ Estimated remaining: {int(remaining)}s")

        progress.progress(1.0, text="Complete!")
        st.success("✅ Ligand preparation complete")
        timer.empty()
        st.divider()
        c1, c2 = st.columns(2)
        c1.metric("✅ Prepared", len(prepared))
        c2.metric("❌ Failed",   len(failed))
        if failed:
            st.warning(
                "**Tip:** MGLTools struggles with sugars/glycosides. "
                "Try enabling **Treat as rigid** for these compounds."
            )


# ══════════════════════════════════════════════════════════════════════════════
# TAB 3 — DOCKING
# ══════════════════════════════════════════════════════════════════════════════
with tab3:
    st.markdown("### Batch Docking")
    st.caption(
        "Runs AutoDock Vina on all receptor–ligand pairs in parallel. "
        "Grid boxes are read from the configs generated in Protein Prep."
    )

    # ── Pre-flight checklist — always re-read filesystem live ────────────────
    # Re-read counts fresh here (not from cached _n_ variables) so the checklist
    # reflects the actual state after preparation runs complete.
    _pf_raw_prot  = len(list(raw_prot.glob("*.pdb")))
    _pf_raw_lig   = len([f for f in raw_lig.glob("*") if f.suffix.lower() in {".sdf", ".mol2", ".pdb"}])
    _pf_prep_rec  = len(list(prep_rec_dir.glob("*_receptor.pdbqt")))
    _pf_prep_lig  = len(list(prep_lig_dir.glob("*.pdbqt")))
    _pf_all_ready = (_pf_prep_rec > 0 and _pf_prep_lig > 0)

    with st.expander(" Pre-docking Checklist", expanded=not _pf_all_ready):
        _checks = [
            (
                _pf_raw_prot > 0,
                f"Raw proteins: {_pf_raw_prot} PDB file(s) in raw_proteins/",
                "No .pdb files found in raw_proteins/ — add your protein structures."
            ),
            (
                _pf_raw_lig > 0,
                f"Raw ligands: {_pf_raw_lig} file(s) in raw_ligands/",
                "No ligand files (.sdf/.mol2/.pdb) found in raw_ligands/ — add your compounds."
            ),
            (
                _pf_prep_rec > 0,
                f"Proteins prepared: {_pf_prep_rec} receptor PDBQT(s) ready",
                "No prepared receptors found — run Protein Prep first."
            ),
            (
                _pf_prep_lig > 0,
                f"Ligands prepared: {_pf_prep_lig} ligand PDBQT(s) ready",
                "No prepared ligands found — run Ligand Prep first."
            ),
            (
                any((grid_dir / p.stem.replace("_receptor","")).glob("grid_config.txt")
                    for p in prep_rec_dir.glob("*_receptor.pdbqt")) if _pf_prep_rec > 0 else False,
                "Grid configs generated for all prepared receptors",
                "Grid configs missing — re-run Protein Prep to regenerate them."
            ),
            (
                Path(vina_path).exists() if vina_path else False,
                "AutoDock Vina executable found at configured path",
                "Vina not found — check the path in the sidebar."
            ),
        ]
        _all_ok = all(ok for ok, _, _ in _checks)
        if _all_ok:
            st.success("✅ All checks passed — you're ready to dock!")
        else:
            st.warning("⚠️ Fix the issues below before running docking.")
        for ok, msg_ok, msg_fail in _checks:
            cls   = "pf-ok" if ok else "pf-err"
            icon  = "✅" if ok else "❌"
            label = msg_ok if ok else msg_fail
            st.markdown(
                f'<div class="pf-row {cls}">{icon}&nbsp; {label}</div>',
                unsafe_allow_html=True,
            )
        st.markdown("")
    st.markdown("")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Grid / search mode**")
        docking_mode = st.radio(
            "Mode",
            ["Auto (from Protein Prep)", "Custom grid box"],
            key="dock_mode", label_visibility="collapsed",
        )
        if docking_mode == "Custom grid box":
            st.markdown("**Custom grid parameters**")
            cc1, cc2 = st.columns(2)
            with cc1:
                center_x = st.number_input("center_x", value=0.0, key="dock_cx")
                center_y = st.number_input("center_y", value=0.0, key="dock_cy")
                center_z = st.number_input("center_z", value=0.0, key="dock_cz")
            with cc2:
                size_x = st.number_input("size_x", value=20, key="dock_sx")
                size_y = st.number_input("size_y", value=20, key="dock_sy")
                size_z = st.number_input("size_z", value=20, key="dock_sz")

    with col2:
        total_cores = os.cpu_count() or 4
        with st.expander("Vina parameters", expanded=True):
            exhaustiveness = st.number_input(
                "Exhaustiveness", 1, 64, 16, key="vina_exhaust",
                help="Higher = more thorough search. 8 is standard; 16–32 for publication-quality results.",
            )
            num_modes = st.number_input(
                "Binding modes to save", 1, 50, 9, key="vina_modes",
                help="Number of top poses stored per ligand.",
            )
            energy_range = st.number_input(
                "Energy range (kcal/mol)", 1, 20, 3, key="vina_erange",
                help="Maximum energy difference from best pose to include in output.",
            )
            cpu_workers = st.number_input(
                "Parallel jobs", 1, 32, 2, key="vina_workers",
                help="Receptor–ligand pairs run simultaneously.",
            )
            cores_per_job = st.number_input(
                "CPU cores per job", 1, total_cores,
                max(1, total_cores // 2), key="vina_cores",
            )
            use_fixed_seed = st.checkbox(
                "Use fixed random seed (reproducible results)", True, key="vina_use_seed",
                help="Vina's search is stochastic. Without a fixed seed, re-running the "
                     "same job produces different poses each time. Enable this — and "
                     "report the seed value in your Methods section — for reproducible, "
                     "publication-quality results.",
            )
            vina_seed = st.number_input(
                "Random seed", 0, 2**31 - 1, 42, key="vina_seed_val",
                disabled=not use_fixed_seed,
            )
            st.caption(
                f"**{int(cpu_workers)} jobs × {int(cores_per_job)} cores "
                f"= {int(cpu_workers)*int(cores_per_job)} / {total_cores} cores used**"
            )
    st.markdown("")

    # ── Grid Box 3D Visualizer (both Auto and Custom modes) ──────────────────
    st.markdown("---")
    _show_gridbox = st.checkbox(
        "🔲 Preview grid box in 3D",
        value=False,
        key="gridbox_show",
        help="Load an interactive 3D viewer showing the docking search space overlaid on the receptor.",
    )

    _rec_files = sorted(prep_rec_dir.glob("*_receptor.pdbqt")) if _show_gridbox else []

    if _show_gridbox and not _rec_files:
        st.info("⚠️ No prepared receptors found — run Protein Prep first to preview the grid box.")
    elif _show_gridbox:
        # ── Receptor selector (always shown — applies to both modes) ──────
        _rec_names = [f.stem.replace("_receptor", "") for f in _rec_files]
        _auto_caption = (
            "Each receptor has its own auto-computed grid box from Protein Prep."
            if docking_mode == "Auto (from Protein Prep)"
            else "Select which receptor to display the custom grid box on."
        )
        if len(_rec_files) > 1:
            _rec_sel_name = st.selectbox(
                "Receptor to preview",
                _rec_names,
                key="gridbox_rec_sel",
                help=_auto_caption,
            )
        else:
            _rec_sel_name = _rec_names[0]
            st.caption(f"Receptor: **{_rec_sel_name}**")

        _rec_pdbqt_path = prep_rec_dir / f"{_rec_sel_name}_receptor.pdbqt"

        # ── Resolve grid parameters depending on mode ─────────────────────
        _grid_params_ok   = False
        _grid_source_note = ""

        if docking_mode == "Auto (from Protein Prep)":
            _auto_cfg = grid_dir / _rec_sel_name / "grid_config.txt"
            if _auto_cfg.exists():
                _params = {}
                for _line in _auto_cfg.read_text().splitlines():
                    if "=" in _line:
                        _k, _v = _line.split("=", 1)
                        _params[_k.strip()] = float(_v.strip())
                _cx  = _params.get("center_x", 0.0)
                _cy  = _params.get("center_y", 0.0)
                _cz  = _params.get("center_z", 0.0)
                _sx  = _params.get("size_x",  20.0)
                _sy  = _params.get("size_y",  20.0)
                _sz  = _params.get("size_z",  20.0)
                _grid_params_ok   = True
                _grid_source_note = f"Grid from: `{_auto_cfg}`"
            else:
                st.warning(
                    f"⚠️ No grid config found for **{_rec_sel_name}** "
                    f"(`grid/{_rec_sel_name}/grid_config.txt` missing). "
                    "Re-run Protein Prep to generate it."
                )
        else:
            # Custom mode — use the live widget values
            _cx = float(center_x)
            _cy = float(center_y)
            _cz = float(center_z)
            _sx = float(size_x)
            _sy = float(size_y)
            _sz = float(size_z)
            _grid_params_ok   = True
            _grid_source_note = "Grid from: custom input above"

        if _grid_params_ok:
            # ── Style controls ────────────────────────────────────────────
            _gv_c1, _gv_c2, _gv_c3 = st.columns(3)
            with _gv_c1:
                _box_color = st.selectbox(
                    "Box colour", ["cyan", "magenta", "yellow", "lime", "orange", "white"],
                    key="gridbox_color",
                )
            with _gv_c2:
                _rec_style = st.selectbox(
                    "Receptor style", ["cartoon", "surface", "line", "stick"],
                    key="gridbox_rec_style",
                )
            with _gv_c3:
                _box_opacity = st.slider(
                    "Box fill opacity", 0.0, 0.5, 0.15, 0.05,
                    key="gridbox_opacity",
                )

            # ── Show parsed grid values for Auto mode ─────────────────────
            if docking_mode == "Auto (from Protein Prep)":
                _pi_c1, _pi_c2, _pi_c3, _pi_c4, _pi_c5, _pi_c6 = st.columns(6)
                _pi_c1.metric("center_x", f"{_cx:.2f}")
                _pi_c2.metric("center_y", f"{_cy:.2f}")
                _pi_c3.metric("center_z", f"{_cz:.2f}")
                _pi_c4.metric("size_x",   f"{_sx:.0f}")
                _pi_c5.metric("size_y",   f"{_sy:.0f}")
                _pi_c6.metric("size_z",   f"{_sz:.0f}")

            try:
                _rec_text    = _rec_pdbqt_path.read_text(errors="replace")
                _rec_text_js = _rec_text.replace("\\", "\\\\").replace("`", "\\`")

                _style_js = {
                    "cartoon": '{"cartoon": {"color": "spectrum"}}',
                    "surface": '{"surface": {"opacity": 0.7, "color": "white"}}',
                    "line":    '{"line": {"colorscheme": "greenCarbon"}}',
                    "stick":   '{"stick": {"colorscheme": "greenCarbon"}}',
                }[_rec_style]

                _hx  = _sx / 2
                _hy  = _sy / 2
                _hz  = _sz / 2
                _vol = _sx * _sy * _sz

                _corners_js = "[" + ",".join(
                    f"[{_cx+dx:.3f},{_cy+dy:.3f},{_cz+dz:.3f}]"
                    for dx in (-_hx, _hx)
                    for dy in (-_hy, _hy)
                    for dz in (-_hz, _hz)
                ) + "]"

                _mode_label = "Auto" if docking_mode == "Auto (from Protein Prep)" else "Custom"

                _gridbox_html = f"""<!DOCTYPE html>
<html>
<head>
<style>
  body {{ margin:0; background:#1a1a2e; }}
  #viewer {{ width:100%; height:460px; position:relative; }}
  #overlay {{
    position:absolute; top:10px; left:10px; z-index:100;
    background:rgba(0,0,0,0.65); color:#e0e0e0;
    font-family:'IBM Plex Mono',monospace; font-size:11px;
    padding:8px 12px; border-radius:4px; line-height:1.8;
    border:1px solid rgba(255,255,255,0.15); pointer-events:none;
  }}
  #overlay b {{ color:#fff; }}
  #overlay span {{ color:{_box_color}; font-weight:bold; }}
  #mode-badge {{
    position:absolute; top:10px; right:10px; z-index:100;
    background:rgba(26,74,122,0.85); color:#fff;
    font-family:'IBM Plex Mono',monospace; font-size:10px; font-weight:600;
    padding:4px 10px; border-radius:3px; letter-spacing:0.06em;
    text-transform:uppercase; pointer-events:none;
    border:1px solid rgba(255,255,255,0.2);
  }}
</style>
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"></script>
</head>
<body>
<div style="position:relative">
  <div id="viewer"></div>
  <div id="overlay">
    <b>Grid Box</b><br>
    Center:&nbsp;<span>({_cx:.2f}, {_cy:.2f}, {_cz:.2f})</span><br>
    Size:&nbsp;&nbsp;&nbsp;<span>{_sx:.1f} &times; {_sy:.1f} &times; {_sz:.1f} &Aring;</span><br>
    Volume:&nbsp;<span>{_vol:.0f} &Aring;&sup3;</span><br>
    Protein:<span>&nbsp;{_rec_sel_name}</span>
  </div>
  <div id="mode-badge">{_mode_label} grid</div>
</div>
<script>
const viewer = $3Dmol.createViewer(document.getElementById("viewer"), {{backgroundColor:"#1a1a2e"}});
const pdb = `{_rec_text_js}`;
viewer.addModel(pdb, "pdbqt");
viewer.setStyle({{}}, {_style_js});

viewer.addBox({{center:{{x:{_cx},y:{_cy},z:{_cz}}}, dimensions:{{w:{_sx},h:{_sy},d:{_sz}}},
  color:"{_box_color}", opacity:{_box_opacity}, wireframe:false}});
viewer.addBox({{center:{{x:{_cx},y:{_cy},z:{_cz}}}, dimensions:{{w:{_sx},h:{_sy},d:{_sz}}},
  color:"{_box_color}", opacity:0.95, wireframe:true}});

const corners = {_corners_js};
corners.forEach(([x,y,z]) => viewer.addSphere({{
  center:{{x,y,z}}, radius:0.35, color:"{_box_color}", opacity:0.9}}));

viewer.addSphere({{center:{{x:{_cx},y:{_cy},z:{_cz}}},
  radius:0.6, color:"red", opacity:1.0}});

viewer.zoomTo(); viewer.render(); viewer.zoom(0.85);
</script>
</body>
</html>"""

                import streamlit.components.v1 as _grid_comp
                _grid_comp.html(_gridbox_html, height=470, scrolling=False)
                st.caption(
                    f"🔴 Red sphere = box centre  ·  {_grid_source_note}  ·  "
                    "Rotate: left-drag  ·  Zoom: scroll  ·  Pan: right-drag"
                )

            except Exception as _gv_err:
                st.warning(f"⚠️ Grid box preview unavailable: {_gv_err}")

    st.markdown("")

    skip_existing = st.checkbox(
        "⏭ Skip already-docked pairs (resume mode)",
        value=True,
        key="dock_skip_existing",
        help="If a docking output PDBQT already exists and is non-empty, skip that job. "
             "Useful for resuming interrupted runs without re-doing finished jobs.",
    )

    if int(exhaustiveness) < 16:
        st.warning(
            f"⚠️ **Exhaustiveness = {int(exhaustiveness)}** is suitable for quick screening "
            f"but results may not be reproducible between runs. "
            f"Use **≥ 16** for publication-quality results, **≥ 32** for high-confidence poses."
        )

    if st.button("▶ Run Docking", type="primary", key="run_dock"):
        receptors = list(prep_rec_dir.glob("*_receptor.pdbqt"))
        ligands   = list(prep_lig_dir.glob("*.pdbqt"))

        if not receptors:
            st.error("❌ No prepared receptors found. Run Protein Prep first.")
            st.stop()
        if not ligands:
            st.error("❌ No prepared ligands found. Run Ligand Prep first.")
            st.stop()

        # Build job list
        jobs = []
        for rec in receptors:
            pname = rec.stem.replace("_receptor", "")
            if docking_mode == "Custom grid box":
                cfg = dock_dir / f"{pname}_custom.txt"
                cfg.write_text(
                    f"center_x = {center_x}\ncenter_y = {center_y}\ncenter_z = {center_z}\n\n"
                    f"size_x = {size_x}\nsize_y = {size_y}\nsize_z = {size_z}\n"
                )
            else:
                cfg = grid_dir / pname / "grid_config.txt"
            if not cfg.exists():
                st.warning(f"⚠️ No grid config for **{pname}** — skipping.")
                continue
            for lig in ligands:
                # Resume mode: skip if output already exists and is non-empty
                out_pdbqt_check = dock_dir / f"{pname}_{lig.stem}.pdbqt"
                if skip_existing and out_pdbqt_check.exists() and out_pdbqt_check.stat().st_size > 0:
                    continue
                jobs.append((rec, lig, cfg, vina_path,
                             int(exhaustiveness), int(num_modes), int(energy_range),
                             int(cores_per_job), dock_dir, result_dir,
                             int(vina_seed) if use_fixed_seed else None))

        if not jobs:
            st.error("❌ No valid docking jobs. Check protein/ligand files and grid configs.")
            st.stop()

        _total_possible = len(receptors) * len(ligands)
        _skipped_count  = _total_possible - len(jobs)
        _skip_note = f" · {_skipped_count} skipped (already done)" if _skipped_count > 0 else ""
        if not jobs:
            st.success(f"✅ All {_total_possible} jobs already completed (resume mode). Nothing to do.")
            st.stop()
        st.info(f"Running **{len(jobs)}** docking jobs ({len(receptors)} receptors × {len(ligands)} ligands){_skip_note}…")

        progress  = st.progress(0)
        log_area  = st.empty()
        timer     = st.empty()
        t0        = time.time()
        counters  = {"done": 0, "success": 0, "failed": 0}
        log_lines = []

        progress.progress(0, text="Docking in progress…")
        with ThreadPoolExecutor(max_workers=int(cpu_workers)) as pool:
            futures = {pool.submit(run_single_docking, job): job for job in jobs}
            for future in as_completed(futures):
                outcome, pname, lname, msg = future.result()
                counters["done"] += 1
                if outcome == "success":
                    counters["success"] += 1
                    log_lines.append(f"✅  {pname}  +  {lname}")
                else:
                    counters["failed"] += 1
                    log_lines.append(f"❌  {pname}  +  {lname}  — {msg}")

                _pct = counters["done"] / len(jobs)
                _txt = f"{counters['done']}/{len(jobs)} done — ✅ {counters['success']}  ❌ {counters['failed']}"
                progress.progress(_pct, text=_txt)
                log_area.text_area("Live log", "\n".join(log_lines[-30:]),
                                   height=160)
                elapsed   = time.time() - t0
                remaining = elapsed / counters["done"] * (len(jobs) - counters["done"])
                timer.caption(f"⏱ Estimated remaining: {int(remaining)}s")

        progress.progress(1.0)
        st.success("✅ All docking jobs complete")
        timer.empty()
        time.sleep(0.5)  # Allow WebSocket to stabilize before results render
        st.divider()
        c1, c2, c3 = st.columns(3)
        c1.metric("Total jobs",  len(jobs))
        c2.metric("✅ Successful", counters["success"])
        c3.metric("❌ Failed",     counters["failed"])


# ══════════════════════════════════════════════════════════════════════════════
# TAB 4 — RESULTS
# ══════════════════════════════════════════════════════════════════════════════
with tab4:
    st.markdown("### Docking Results")
    st.caption("Scores are parsed directly from Vina log files in `results/`.")

    # Split result file stems into (protein, ligand)
    _known_receptors = sorted(
        {p.stem.replace("_receptor", "") for p in prep_rec_dir.glob("*_receptor.pdbqt")},
        key=len, reverse=True,
    )

    def _split_stem(stem: str):
        for rec in _known_receptors:
            if stem.startswith(rec + "_"):
                return rec, stem[len(rec) + 1:]
            if stem == rec:
                return rec, "unknown"
        parts = stem.split("_", 1)
        return parts[0], (parts[1] if len(parts) > 1 else "unknown")

    rows = []
    for log in result_dir.glob("*.txt"):
        protein, ligand = _split_stem(log.stem)
        all_modes = parse_vina_all_modes(log)
        best_aff  = parse_vina_log(log)

        # Collect per-mode affinities
        mode_affs = {m: a for m, a in all_modes}
        n_modes   = len(all_modes)

        # ΔE between mode 1 and mode 2 (selectivity indicator)
        mode1_aff = mode_affs.get(1, best_aff)
        mode2_aff = mode_affs.get(2, None)
        delta_e   = round(mode2_aff - mode1_aff, 3) if (mode2_aff is not None and mode1_aff is not None) else None

        # Ligand efficiency — count heavy atoms from ONE pose only.
        # The docked output PDBQT contains all 9 Vina modes stacked as MODEL
        # blocks, so reading the whole file gives N_atoms × N_modes (e.g. 297
        # instead of 33). Fix: extract only the first MODEL block.
        # Fallback: use the prepared-ligand PDBQT (single pose, cleanest source).
        dock_pdbqt = dock_dir / f"{protein}_{ligand}.pdbqt"
        prep_pdbqt = prep_lig_dir / f"{ligand}.pdbqt"
        n_heavy = None

        def _count_heavy_first_pose(path: Path) -> int | None:
            """Return heavy-atom count from the first MODEL block (or whole file)."""
            try:
                raw = path.read_text(errors="replace")
                lines = raw.splitlines(keepends=True)
                has_models = any(l.startswith("MODEL") for l in lines)
                if has_models:
                    # Collect only lines inside MODEL 1
                    pose_lines, in_model = [], False
                    for line in lines:
                        if line.startswith("MODEL"):
                            if in_model:          # second MODEL reached → stop
                                break
                            in_model = True
                        elif line.startswith("ENDMDL"):
                            break
                        elif in_model:
                            pose_lines.append(line)
                else:
                    pose_lines = lines
                heavy = parse_pdbqt_heavy_atoms("".join(pose_lines))
                return len(heavy) if heavy else None
            except Exception:
                return None

        # Prefer prepared-ligand PDBQT (single pose) → fall back to docked output
        for _src in (prep_pdbqt, dock_pdbqt):
            if _src.exists():
                n_heavy = _count_heavy_first_pose(_src)
                if n_heavy:
                    break
        le = round(-mode1_aff / n_heavy, 3) if (mode1_aff is not None and n_heavy and n_heavy > 0) else None

        rows.append({
            "Protein":               protein,
            "Ligand":                ligand,
            "Best Affinity (kcal/mol)": best_aff,
            "Mode 2 (kcal/mol)":     mode2_aff,
            "ΔE Mode1→2 (kcal/mol)": delta_e,
            "Ligand Efficiency":     le,
            "Heavy Atoms":           n_heavy,
            "Modes Saved":           n_modes,
            "Log File":              log.name,
        })

    if not rows:
        st.info("No results found yet. Run docking in the ** Docking** tab first.")
    else:
        df = pd.DataFrame(rows).sort_values("Best Affinity (kcal/mol)")

        # ── Summary metrics ────────────────────────────────────────────────
        valid = df["Best Affinity (kcal/mol)"].dropna()
        m1, m2, m3, m4, m5 = st.columns(5)
        m1.metric("Total dockings",    len(df))
        m2.metric("Best (kcal/mol)",   f"{valid.min():.2f}"  if len(valid) else "—")
        m3.metric("Mean (kcal/mol)",   f"{valid.mean():.2f}" if len(valid) else "—")
        m4.metric("Unique ligands",    df["Ligand"].nunique())
        m5.metric("Unique proteins",   df["Protein"].nunique())
        st.markdown("")

        # ── Score table ────────────────────────────────────────────────────
        st.subheader(" Score Table")

        # Column display order
        _display_cols = [
            "Protein", "Ligand",
            "Best Affinity (kcal/mol)", "Mode 2 (kcal/mol)",
            "ΔE Mode1→2 (kcal/mol)", "Ligand Efficiency",
            "Heavy Atoms", "Modes Saved", "Log File",
        ]
        _display_cols = [c for c in _display_cols if c in df.columns]

        # ── Interactive filters ────────────────────────────────────────────
        with st.expander(" Filter & sort", expanded=False):
            fc1, fc2, fc3 = st.columns(3)
            with fc1:
                _prot_opts = ["All proteins"] + sorted(df["Protein"].unique().tolist())
                _prot_sel  = st.selectbox("Protein", _prot_opts, key="res_prot_filter")
            with fc2:
                _aff_min = float(df["Best Affinity (kcal/mol)"].min() if len(valid) else -15)
                _aff_max = float(df["Best Affinity (kcal/mol)"].max() if len(valid) else 0)
                _aff_max = min(_aff_max, 0.0)
                # Guard: slider requires min < max. When all results share the
                # same affinity (e.g. only one docking run), pad by 0.5 kcal/mol.
                if _aff_min >= _aff_max:
                    _aff_min = _aff_max - 0.5
                _aff_rng = st.slider(
                    "Affinity range (kcal/mol)", _aff_min, _aff_max,
                    (_aff_min, _aff_max), step=0.5, key="res_aff_slider",
                )
            with fc3:
                _lig_search = st.text_input("Search ligand name", placeholder="e.g. Imatinib", key="res_lig_search")

        df_view = df.copy()
        if _prot_sel != "All proteins":
            df_view = df_view[df_view["Protein"] == _prot_sel]
        if len(valid) > 0:
            df_view = df_view[
                df_view["Best Affinity (kcal/mol)"].isna() |
                df_view["Best Affinity (kcal/mol)"].between(_aff_rng[0], _aff_rng[1])
            ]
        if _lig_search.strip():
            df_view = df_view[df_view["Ligand"].str.contains(_lig_search.strip(), case=False, na=False)]

        st.caption(f"Showing **{len(df_view)}** of **{len(df)}** results")

        st.dataframe(
            df_view[_display_cols].style.format({
                "Best Affinity (kcal/mol)": lambda v: f"{v:.2f}" if pd.notna(v) else "—",
                "Mode 2 (kcal/mol)":        lambda v: f"{v:.2f}" if pd.notna(v) else "—",
                "ΔE Mode1→2 (kcal/mol)":    lambda v: f"{v:.3f}" if pd.notna(v) else "—",
                "Ligand Efficiency":         lambda v: f"{v:.3f}" if pd.notna(v) else "—",
            }).background_gradient(
                subset=["Best Affinity (kcal/mol)"], cmap="RdYlGn_r", vmin=-12, vmax=-3,
            ),
            use_container_width=True, hide_index=True,
        )
        st.caption(
            "**Best Affinity** = Vina mode 1 score · "
            "**ΔE Mode1→2** = energy gap (larger = more selective pose) · "
            "**Ligand Efficiency** = −(Best Affinity) ÷ Heavy Atoms"
        )

        # ── Download buttons ───────────────────────────────────────────────
        csv_path = result_dir / "docking_results.csv"
        df.to_csv(csv_path, index=False)
        col_dl1, col_dl2, col_note = st.columns([1, 1, 2])
        with col_dl1:
            st.download_button(
                "⬇ Download CSV",
                data=df_view[_display_cols].to_csv(index=False),
                file_name="docking_results.csv",
                mime="text/csv",
                use_container_width=True,
                key="dl_csv_results",
            )
        with col_dl2:
            try:
                import matplotlib.backends.backend_pdf as _pdf_backend
                _pdf_buf = io.BytesIO()
                with _pdf_backend.PdfPages(_pdf_buf) as _pdf:
                    # Page 1: summary table
                    _fig_t, _ax_t = plt.subplots(figsize=(11, max(3, len(df_view) * 0.32 + 1.5)))
                    _ax_t.axis("off")
                    _tbl_cols = [c for c in _display_cols if c != "Log File"]
                    _tbl_data = [
                        [str(df_view.iloc[r][c])[:18] if pd.notna(df_view.iloc[r][c]) else "—"
                         for c in _tbl_cols]
                        for r in range(min(len(df_view), 40))
                    ]
                    _tbl = _ax_t.table(
                        cellText=_tbl_data, colLabels=_tbl_cols,
                        loc="center", cellLoc="center",
                    )
                    _tbl.auto_set_font_size(False)
                    _tbl.set_fontsize(7)
                    _tbl.auto_set_column_width(col=list(range(len(_tbl_cols))))
                    for (row, col), cell in _tbl.get_celld().items():
                        if row == 0:
                            cell.set_facecolor("#0066cc")
                            cell.set_text_props(color="white", fontweight="bold")
                        elif row % 2 == 0:
                            cell.set_facecolor("#f4f8fd")
                    _ax_t.set_title(
                        f"IBDock Docking Report  —  {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}",
                        fontsize=11, fontweight="bold", pad=14, color="#0a2540",
                    )
                    plt.tight_layout()
                    _pdf.savefig(_fig_t, bbox_inches="tight")
                    plt.close(_fig_t)
                    # Page 2: heatmap
                    try:
                        _piv = df_view.pivot(index="Ligand", columns="Protein", values="Best Affinity (kcal/mol)")
                        _fig_h, _ax_h = plt.subplots(figsize=(max(5, len(_piv.columns)*1.6), max(3, len(_piv)*0.7)))
                        sns.heatmap(_piv, annot=True, fmt=".2f", cmap="RdYlGn_r", ax=_ax_h,
                                    linewidths=0.4, linecolor="#e0eaf4",
                                    cbar_kws={"label": "Best Affinity (kcal/mol)"})
                        _ax_h.set_title("Binding Affinity Heatmap (kcal/mol)", fontsize=11, pad=12)
                        plt.tight_layout()
                        _pdf.savefig(_fig_h, bbox_inches="tight")
                        plt.close(_fig_h)
                    except Exception:
                        pass
                    # Page 3: bar chart
                    try:
                        _prots = df_view["Protein"].unique()
                        _fig_b, _axes_b = plt.subplots(
                            1, len(_prots),
                            figsize=(max(5, len(_prots)*3.5), max(3, df_view["Ligand"].nunique()*0.45)),
                            squeeze=False,
                        )
                        for _ax_b, _prot in zip(_axes_b[0], _prots):
                            _sub = df_view[df_view["Protein"]==_prot].sort_values("Best Affinity (kcal/mol)")
                            _clrs = ["#0066cc" if v < -7 else "#1a9e75" if v < -5 else "#ef9f27"
                                     for v in _sub["Best Affinity (kcal/mol)"].fillna(0)]
                            _ax_b.barh(_sub["Ligand"], _sub["Best Affinity (kcal/mol)"],
                                       color=_clrs, edgecolor="white", linewidth=0.3)
                            _ax_b.axvline(-7, color="#e24b4a", linewidth=1, linestyle="--", alpha=0.6)
                            _ax_b.set_xlabel("Best Affinity (kcal/mol)", fontsize=9)
                            _ax_b.set_title(_prot, fontsize=10, fontweight="bold")
                            _ax_b.spines[["top","right"]].set_visible(False)
                        plt.tight_layout()
                        _pdf.savefig(_fig_b, bbox_inches="tight")
                        plt.close(_fig_b)
                    except Exception:
                        pass
                    _d = _pdf.infodict()
                    _d["Title"]   = "IBDock Docking Report"
                    _d["Subject"] = "AutoDock Vina Virtual Screening Results"
                _pdf_buf.seek(0)
                st.download_button(
                    "⬇ Download PDF Report",
                    data=_pdf_buf.getvalue(),
                    file_name=f"IBDock_report_{pd.Timestamp.now().strftime('%Y%m%d_%H%M')}.pdf",
                    mime="application/pdf",
                    use_container_width=True,
                    key="dl_pdf_results",
                )
            except Exception as _pdf_err:
                st.caption(f"PDF unavailable: {_pdf_err}")
        with col_note:
            st.caption(f"Auto-saved to: `{csv_path}`")


        # ── Heatmap ────────────────────────────────────────────────────────
        st.subheader(" Affinity Heatmap")
        _heatmap_placeholder = st.empty()
        try:
            pivot = df.pivot(index="Ligand", columns="Protein", values="Best Affinity (kcal/mol)")
            fig, ax = plt.subplots(
                figsize=(max(5, len(pivot.columns) * 1.6), max(3, len(pivot) * 0.7))
            )
            sns.heatmap(
                pivot, annot=True, fmt=".2f", cmap="RdYlGn_r",
                ax=ax, linewidths=0.4, linecolor="#e0eaf4",
                cbar_kws={"label": "Best Affinity (kcal/mol)"},
            )
            ax.set_title("AutoDock Vina — Binding Affinity (kcal/mol)", fontsize=11, pad=12)
            ax.set_xlabel("")
            ax.set_ylabel("")
            plt.tight_layout()
            _heatmap_placeholder.pyplot(fig)
            plt.close(fig)
        except Exception:
            _heatmap_placeholder.info("Heatmap requires at least one protein and one ligand with valid scores.")

        # ── Bar chart — best affinity per ligand ───────────────────────────
        if len(df) >= 2:
            st.subheader(" Best Affinity Per Ligand")
            _bar_placeholder = st.empty()
            proteins_list = df["Protein"].unique()
            ncols  = min(len(proteins_list), 3)
            fig2, axes = plt.subplots(
                1, len(proteins_list),
                figsize=(max(5, len(proteins_list) * 3.5), max(3, len(df["Ligand"].unique()) * 0.45)),
                squeeze=False,
            )
            for ax, prot in zip(axes[0], proteins_list):
                sub = df[df["Protein"] == prot].sort_values("Best Affinity (kcal/mol)")
                colors = ["#0066cc" if v < -7 else "#1a9e75" if v < -5 else "#ef9f27"
                          for v in sub["Best Affinity (kcal/mol)"].fillna(0)]
                ax.barh(sub["Ligand"], sub["Best Affinity (kcal/mol)"],
                        color=colors, edgecolor="white", linewidth=0.3)
                ax.axvline(-7, color="#e24b4a", linewidth=1, linestyle="--", alpha=0.6, label="-7 kcal/mol")
                ax.set_xlabel("Best Affinity (kcal/mol)", fontsize=9)
                ax.set_title(prot, fontsize=10, fontweight="bold")
                ax.spines[["top", "right"]].set_visible(False)
                ax.legend(fontsize=7)
            plt.tight_layout()
            _bar_placeholder.pyplot(fig2)
            plt.close(fig2)
            st.caption("Blue: strong (< −7) · Green: moderate (−5 to −7) · Orange: weak (> −5)")


# ══════════════════════════════════════════════════════════════════════════════
# TAB 5 — 3D POSE VIEWER
# ══════════════════════════════════════════════════════════════════════════════
import streamlit.components.v1 as _components

with tab5:
    st.markdown("### 3D Pose Viewer")
    st.caption("Interactive 3D visualisation of docked poses.")

    # Build pose index
    _rec_stems = sorted(
        [r.stem.replace("_receptor", "") for r in prep_rec_dir.glob("*_receptor.pdbqt")],
        key=len, reverse=True,
    )

    def _split_pose_stem(stem: str):
        for rs in _rec_stems:
            if stem.startswith(rs + "_") or stem == rs:
                return rs, stem[len(rs):].lstrip("_") or "unknown"
        parts = stem.split("_", 1)
        return parts[0], (parts[1] if len(parts) > 1 else "unknown")

    pose_options = {}
    for f in sorted(dock_dir.glob("*.pdbqt")):
        prot, lig = _split_pose_stem(f.stem)
        pose_options[f"{prot}  +  {lig}"] = (prot, lig, f)

    if not pose_options:
        st.info("No docked poses found. Run docking in the ** Docking** tab first.")
    else:
        pv_c1, pv_c2, pv_c3 = st.columns([3, 1, 1])
        with pv_c1:
            selected = st.selectbox("Select protein + ligand pair", list(pose_options.keys()), key="pv_pair")
        with pv_c2:
            show_all = st.toggle("All poses", False, key="pv_all",
                                 help="Show all 9 Vina binding modes.")
        with pv_c3:
            st.markdown("")  # spacer

        protein_name, ligand_name, pose_file = pose_options[selected]
        receptor_file = prep_rec_dir / f"{protein_name}_receptor.pdbqt"

        # Fuzzy-match receptor if exact name doesn't exist
        if not receptor_file.exists():
            for r in prep_rec_dir.glob("*_receptor.pdbqt"):
                if pose_file.stem.startswith(r.stem.replace("_receptor", "")):
                    receptor_file = r
                    protein_name  = r.stem.replace("_receptor", "")
                    ligand_name   = pose_file.stem[len(protein_name):].lstrip("_")
                    break

        if not receptor_file.exists():
            st.warning(
                f"⚠️ Receptor PDBQT not found for **{pose_file.stem}**. "
                "Re-run Protein Prep to regenerate it."
            )
        else:
            # Affinity lookup
            vina_log  = result_dir / f"{pose_file.stem}.txt"
            affinity  = parse_vina_log(vina_log) if vina_log.exists() else None

            try:
                rec_pdb = pdbqt_to_pdb_text(receptor_file, first_pose_only=False)
                # Convert ligand to SDF for correct bond-order rendering in py3Dmol.
                # PDBQT has no CONECT/bond records; loading it as 'pdb' causes
                # py3Dmol to infer bonds purely from atom distances, which collapses
                # complex ligands into benzene-ring artefacts.
                _obabel = obabel_path  # sidebar-configured Open Babel path
                lig_sdf = ""
                if _obabel:
                    lig_sdf = pdbqt_to_sdf_text(
                        pose_file, _obabel, first_pose_only=not show_all
                    )
                if lig_sdf:
                    lig_pdb    = lig_sdf   # variable name kept for download buttons
                    lig_format = "sdf"
                else:
                    # Fallback: raw PDB text (may render incorrectly for complex ligands)
                    lig_pdb    = pdbqt_to_pdb_text(pose_file, first_pose_only=not show_all)
                    lig_format = "pdb"
            except Exception as exc:
                st.error(f"❌ Could not read pose files: {exc}")
                st.stop()

            # Binding-site contact analysis (residues within 4.5 Å, classified
            # by interaction type). Uses the same rec_pdb/lig_pdb/lig_format
            # already built above, so no extra file conversion is needed.
            try:
                pose_contacts = compute_ligand_contacts(rec_pdb, lig_pdb, lig_format)
                _contact_keys = {(c["chain"], c["resnum"]) for c in pose_contacts}
                pose_contact_pdb = filter_pdb_to_residues(rec_pdb, _contact_keys)
            except Exception as _contact_exc:
                pose_contacts, pose_contact_pdb = [], ""
                st.caption(f"⚠️ Binding-site contact analysis unavailable: {_contact_exc}")

            # Engine badge
            st.markdown(
                '<span style="font-family:\'IBM Plex Mono\',monospace; font-size:0.72rem; color:#444;'
                'border:1px solid #ccc; padding:2px 10px; border-radius:2px;">AutoDock Vina</span>',
                unsafe_allow_html=True,
            )
            st.markdown("")

            # 3D viewer
            viewer_html = build_viewer_html(rec_pdb, lig_pdb, protein_name, ligand_name, affinity,
                                             lig_format, contacts=pose_contacts, contact_pdb=pose_contact_pdb)
            _components.html(viewer_html, height=580, scrolling=False)

            # Stats row
            st.markdown("")
            s1, s2, s3 = st.columns(3)
            s1.metric("Protein",  protein_name)
            s2.metric("Ligand",   ligand_name)
            s3.metric("Best affinity", f"{affinity:.2f} kcal/mol" if affinity else "—")

            # Binding-site residue table (Schrodinger-style interaction summary)
            st.markdown("")
            if pose_contacts:
                st.markdown("#### Binding-site residues")
                st.caption(
                    "Distance-based contact analysis (\u2264 4.5 \u00c5 heavy-atom distance). "
                    "H-bond / salt-bridge / hydrophobic / \u03c0-stacking flags are geometric "
                    "heuristics (no donor/acceptor angle check, no \u03c0-cation detection; "
                    "salt bridges need RDKit-parsed formal charges, so they only resolve when "
                    "the ligand loaded as SDF, and a HIS salt-bridge call is a softer, "
                    "worth-checking-manually result since HIS is often neutral) "
                    "-- treat as a starting point for a figure, not a validated replacement "
                    "for a tool like PLIP if your Methods section needs a rigorously defined footprint."
                )
                _contacts_df = pd.DataFrame([{
                    "Residue": c["residue"],
                    "Chain": c["chain"],
                    "Interaction": c["interaction"],
                    "Min. distance (Å)": c["min_dist"],
                    "Detail": c["detail"],
                } for c in pose_contacts])
                st.dataframe(_contacts_df, use_container_width=True, hide_index=True)
                st.download_button(
                    "⬇ Download binding-site residue table (CSV)",
                    data=_contacts_df.to_csv(index=False),
                    file_name=f"{protein_name}_{ligand_name}_binding_site_residues.csv",
                    mime="text/csv",
                    use_container_width=True,
                )
            else:
                st.caption("No binding-site residues resolved within 4.5 Å (check receptor/ligand parsing above).")

            # Downloads
            st.markdown("")
            d1, d2 = st.columns(2)
            with d1:
                st.download_button(
                    "⬇ Download receptor PDB",
                    data=rec_pdb,
                    file_name=f"{protein_name}_receptor.pdb",
                    mime="text/plain",
                    use_container_width=True,
                )
            with d2:
                st.download_button(
                    "⬇ Download ligand pose PDB",
                    data=lig_pdb,
                    file_name=f"{protein_name}_{ligand_name}_pose.pdb",
                    mime="text/plain",
                    use_container_width=True,
                )

            # Download complex (receptor + ligand together)
            st.markdown("")
            _complex_pdb = rec_pdb.rstrip() + "\nENDMDL\n" + lig_pdb
            st.download_button(
                "⬇ Download full complex PDB  *(receptor + ligand)*",
                data=_complex_pdb,
                file_name=f"{protein_name}_{ligand_name}_complex.pdb",
                mime="text/plain",
                use_container_width=True,
            )

# ══════════════════════════════════════════════════════════════════════════════
# TAB 6 — ABOUT & CITATION
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
# TAB 6 — VALIDATION
# ══════════════════════════════════════════════════════════════════════════════
with tab6:
    st.markdown("### Docking Validation")

    # ── Load results CSV (needed for PDF report) ─────────────────────────────
    results_csv = project_dir / "docking_results.csv"
    df_results  = None
    if results_csv.exists():
        try:
            df_results = pd.read_csv(results_csv)
        except Exception:
            pass

    docked_pdbqts = sorted(dock_dir.glob("*.pdbqt"))


    # ── Re-Docking RMSD vs Crystal / Reference Ligand ────────
    st.markdown("#### Re-Docking RMSD vs Reference Ligand")
    st.caption(
        "Upload a reference ligand PDB (e.g. co-crystallised ligand from the PDB entry) for each protein. "
        "Docked pose 1 is compared against it. **RMSD < 2.0 Å** is the standard pass criterion."
    )

    with st.expander("Important notes on re-docking validity", expanded=False):
        st.markdown(
            "Re-docking RMSD is only valid when the docked molecule is identical to the co-crystallised "
            "reference ligand (same compound, same atom count). If your docked ligands differ from the "
            "native ligand in the reference PDB, atom counts will not match and RMSD is not meaningful. "
            "\n\n**How to obtain a valid reference:** Open the protein entry on RCSB, go to the Ligand tab, "
            "download the co-crystallised ligand as SDF/PDB, then prepare it with Open Babel to match "
            "your docked ligand format. Name the file to include the protein ID — e.g. `6C9H.pdb` or "
            "`6C9H_R734_ref.pdb` — for automatic pairing."
        )

    uploaded_refs = st.file_uploader(
        "Upload reference ligand PDB file(s)",
        type=["pdb"],
        accept_multiple_files=True,
        key="val_ref_upload",
    )

    if uploaded_refs and docked_pdbqts:
        redock_rows = []
        pairs_evaluated = 0

        def _status_band(v):
            if v < 1.0:
                return "✅ Excellent (< 1 Å)"
            elif v < 2.0:
                return "✅ Pass (< 2 Å)"
            elif v < 3.0:
                return "🟡 Borderline (2–3 Å)"
            else:
                return "🔴 Fail (≥ 3 Å)"

        for ref_file in uploaded_refs:
            ref_text  = ref_file.read().decode("utf-8", errors="replace")
            ref_lines = ref_text.splitlines(keepends=True)
            ref_heavy = [l for l in ref_lines
                         if l.startswith(("ATOM", "HETATM")) and not _is_hydrogen(l)]

            if not ref_heavy:
                pairs_evaluated += 1
                redock_rows.append({
                    "Reference File": ref_file.name,
                    "Protein": "—", "Ligand": "—",
                    "Method": "—",
                    "Top-1 RMSD (Å)": None,
                    "Best Pose RMSD (Å)": None,
                    "Best Pose #": None,
                    "Ref Heavy Atoms": 0,
                    "Pose Heavy Atoms": None,
                    "Used as Final": "",
                    "Status": "⚠️ No heavy atoms in reference file",
                })
                continue

            # Match reference filename to docked pairs.
            # Convention: "6C9H_R734_ref.pdb" → matches protein+ligand
            #             "R734_ref.pdb"       → matches ligand name only
            #             "6C9H.pdb"           → matches protein, compare all its ligands
            ref_stem = Path(ref_file.name).stem.lower()
            for suffix in ("_ref", "-ref", "_crystal", "-crystal", "_native", "-native"):
                ref_stem = ref_stem.replace(suffix, "")

            # Convert the reference ligand to SDF once per reference file (bond
            # perception via Open Babel), for graph-based RMSD matching that
            # doesn't depend on atom naming — see compute_rmsd_from_sdf().
            ref_sdf_text = pdb_to_sdf_text(ref_text, obabel_path)

            # Match ONLY docked pairs whose protein ID appears in the reference filename.
            # e.g. "6C9H.pdb" → only 6C9H_*.pdbqt; "6C9H1.pdb" → only 6C9H1_*.pdbqt
            # No fallback to ALL pairs — that creates a combinatorial explosion.
            matched = []
            for dpdbqt in docked_pdbqts:
                if ref_stem in dpdbqt.stem.lower():
                    matched.append(dpdbqt)

            if not matched:
                # Try matching just the protein-name portion of the filename
                _prot_names_rd = [r.stem.replace("_receptor", "")
                                  for r in prep_rec_dir.glob("*_receptor.pdbqt")]
                for pn in sorted(_prot_names_rd, key=len, reverse=True):
                    if pn.lower() in ref_stem or ref_stem.startswith(pn.lower()):
                        for dpdbqt in docked_pdbqts:
                            if dpdbqt.stem.lower().startswith(pn.lower() + "_"):
                                matched.append(dpdbqt)
                        break  # stop after first matching protein

            if not matched:
                pairs_evaluated += 1
                redock_rows.append({
                    "Reference File": ref_file.name,
                    "Protein": "—", "Ligand": "—",
                    "Method": "—",
                    "Top-1 RMSD (Å)": None,
                    "Best Pose RMSD (Å)": None,
                    "Best Pose #": None,
                    "Ref Heavy Atoms": len(ref_heavy),
                    "Pose Heavy Atoms": None,
                    "Used as Final": "",
                    "Status": (
                        f"ℹ️ No matching docked pair found for '{ref_file.name}'. "
                        "Name your reference file as PROTEINID.pdb (e.g. 6C9H.pdb) "
                        "or PROTEINID_LIGANDNAME_ref.pdb (e.g. 6C9H_R734_ref.pdb) to auto-match."
                    ),
                })
                continue

            for dpdbqt in matched:
                pairs_evaluated += 1
                stem = dpdbqt.stem
                protein_names = [r.stem.replace("_receptor", "")
                                 for r in prep_rec_dir.glob("*_receptor.pdbqt")]
                prot, lig = stem, "unknown"
                for rec in sorted(protein_names, key=len, reverse=True):
                    if stem.startswith(rec + "_"):
                        prot = rec
                        lig  = stem[len(rec) + 1:]
                        break

                pose1 = extract_vina_pose_by_mode(str(dpdbqt), 1)
                if not pose1:
                    redock_rows.append({
                        "Reference File": ref_file.name,
                        "Protein": prot, "Ligand": lig,
                        "Method": "—",
                        "Top-1 RMSD (Å)": None,
                        "Best Pose RMSD (Å)": None,
                        "Best Pose #": None,
                        "Ref Heavy Atoms": len(ref_heavy),
                        "Pose Heavy Atoms": None,
                        "Used as Final": "",
                        "Status": "⚠️ Pose 1 missing",
                    })
                    continue

                pose_heavy = [l for l in pose1
                              if l.startswith(("ATOM", "HETATM")) and not _is_hydrogen(l)]

                # Fast pre-check only (line-count based) — catches obviously
                # different molecules before we bother with Open Babel/RDKit.
                # NOTE: passing this check does NOT mean the RMSD below is
                # valid — the real correspondence check happens either via
                # the coordinate bridge below or compute_rmsd_from_sdf()'s
                # graph matching.
                n_ref, n_pose, count_err = compute_rmsd_from_lines(ref_heavy, pose_heavy)

                if count_err:
                    redock_rows.append({
                        "Reference File":        ref_file.name,
                        "Protein":               prot,
                        "Ligand":                lig,
                        "Method":                "—",
                        "Top-1 RMSD (Å)":        None,
                        "Best Pose RMSD (Å)":    None,
                        "Best Pose #":           None,
                        "Ref Heavy Atoms":       n_ref,
                        "Pose Heavy Atoms":      n_pose,
                        "Used as Final":         "",
                        "Status":                f"❌ {count_err}",
                    })
                    continue

                # ---------------------------------------------------------------
                # Compute BOTH validation methods for every complex, rather than
                # short-circuiting on whichever succeeds first. The coordinate
                # bridge (ligand-prep provenance, optionally symmetry-refined)
                # remains the authoritative RMSD whenever it applies, since it
                # carries the strongest correspondence guarantee in the
                # pipeline -- but the graph-based sPyRMSD result is now always
                # computed alongside it too, so every complex gets a genuine,
                # independent second measurement. Each method gets its OWN row
                # in the table (rather than extra columns crammed into one
                # row) so it's clear at a glance which method produced which
                # number; exactly one row per complex is marked "Used as
                # Final" -- that is the value counted in the summary metrics,
                # the bar chart, and the PDF report.
                # ---------------------------------------------------------------
                bridge_result = _coordinate_bridge_rmsd(
                    ref_heavy, pose_heavy, lig, prep_lig_dir,
                    ref_sdf_text=ref_sdf_text,
                )
                bridge_rmsd, bridge_note = (None, "")
                if bridge_result is not None:
                    bridge_rmsd, bridge_note = bridge_result

                graph_rmsd, graph_err, graph_note = (None, None, "")
                if not ref_sdf_text.strip():
                    graph_err = "Could not convert reference ligand to SDF (Open Babel failed)"
                else:
                    pose_sdf_text = pdbqt_to_sdf_text(str(dpdbqt), obabel_path, first_pose_only=True, add_hydrogens=False)
                    if not pose_sdf_text.strip():
                        graph_err = "Could not convert docked pose to SDF (Open Babel failed)"
                    else:
                        graph_rmsd, graph_err, graph_note = compute_rmsd_from_sdf(ref_sdf_text, pose_sdf_text)

                # -----------------------------------------------------------
                # BEST-POSE sweep: in addition to pose 1 (Vina's top-ranked
                # mode, "Top-1"), also try every other returned pose and keep
                # whichever gives the lowest RMSD per method ("Best Pose") --
                # the standard second re-docking metric alongside top-1
                # success rate. Both methods are swept independently, using
                # the exact same per-pose logic as pose 1 above, just looped.
                # -----------------------------------------------------------
                n_poses = _count_pdbqt_models(str(dpdbqt))

                bridge_best_rmsd, bridge_best_pose = bridge_rmsd, (1 if bridge_rmsd is not None else None)
                graph_best_rmsd,  graph_best_pose  = graph_rmsd,  (1 if graph_rmsd  is not None else None)

                if n_poses > 1:
                    for mode_num in range(2, n_poses + 1):
                        pose_n_lines = extract_vina_pose_by_mode(str(dpdbqt), mode_num)
                        if not pose_n_lines:
                            continue
                        pose_n_heavy = [l for l in pose_n_lines
                                        if l.startswith(("ATOM", "HETATM")) and not _is_hydrogen(l)]

                        # Bridge: cheap, pure-Python -- re-run per pose directly.
                        if len(pose_n_heavy) == n_ref:
                            res_n = _coordinate_bridge_rmsd(
                                ref_heavy, pose_n_heavy, lig, prep_lig_dir,
                                ref_sdf_text=ref_sdf_text,
                            )
                            if res_n is not None:
                                v = res_n[0]
                                if bridge_best_rmsd is None or v < bridge_best_rmsd:
                                    bridge_best_rmsd, bridge_best_pose = v, mode_num

                        # Graph: needs a fresh Open Babel conversion per pose
                        # (bond perception is geometry-dependent per pose), so
                        # write this mode's atom lines to a small temp PDBQT
                        # and reuse the same pdbqt_to_sdf_text() + sPyRMSD path
                        # pose 1 already went through above.
                        if ref_sdf_text.strip():
                            import tempfile as _tempfile
                            tmp_pose = _tempfile.NamedTemporaryFile(
                                suffix=".pdbqt", delete=False, mode="w")
                            tmp_pose.write("MODEL 1\n")
                            tmp_pose.writelines(pose_n_lines)
                            tmp_pose.write("ENDMDL\n")
                            tmp_pose.close()
                            try:
                                sdf_n = pdbqt_to_sdf_text(
                                    tmp_pose.name, obabel_path,
                                    first_pose_only=True, add_hydrogens=False,
                                )
                                if sdf_n.strip():
                                    v_n, _, _ = compute_rmsd_from_sdf(ref_sdf_text, sdf_n)
                                    if v_n is not None and (graph_best_rmsd is None or v_n < graph_best_rmsd):
                                        graph_best_rmsd, graph_best_pose = v_n, mode_num
                            finally:
                                try:
                                    os.unlink(tmp_pose.name)
                                except Exception:
                                    pass

                # The bridge is authoritative whenever it applies (strongest
                # guarantee, see _coordinate_bridge_rmsd docstring); the graph
                # match is only ever the FINAL value when the bridge does not
                # apply at all. It is never chosen over a successful bridge
                # result just because the two happen to disagree.
                if bridge_rmsd is not None:
                    final_method = "Coordinate bridge"
                elif graph_rmsd is not None:
                    final_method = "Graph isomorphism (sPyRMSD)"
                else:
                    final_method = None

                common = {
                    "Reference File":   ref_file.name,
                    "Protein":          prot,
                    "Ligand":           lig,
                    "Ref Heavy Atoms":  n_ref,
                    "Pose Heavy Atoms": n_pose,
                }

                # One row for the coordinate bridge, only if it was actually
                # attempted (skipped silently when no prepared-ligand PDBQT
                # was retained for this ligand -- that's "not applicable",
                # not a failure worth a row).
                if bridge_result is not None:
                    b_status = _status_band(bridge_rmsd)
                    if bridge_note:
                        b_status += f"  [{bridge_note}]"
                    redock_rows.append({
                        **common,
                        "Method":             "Coordinate bridge",
                        "Top-1 RMSD (Å)":     bridge_rmsd,
                        "Best Pose RMSD (Å)": bridge_best_rmsd,
                        "Best Pose #":        bridge_best_pose,
                        "Used as Final":      "✓" if final_method == "Coordinate bridge" else "",
                        "Status":             b_status,
                    })

                # One row for the graph-based method, whether it succeeded or
                # failed -- showing why sPyRMSD didn't produce a comparable
                # number is itself useful diagnostic information, not noise.
                if graph_rmsd is not None:
                    g_status = _status_band(graph_rmsd)
                    if graph_note:
                        g_status += f"  [{graph_note}]"
                else:
                    g_status = f"⚠️ {graph_err}" if graph_err else "⚠️ Graph match failed"
                redock_rows.append({
                    **common,
                    "Method":             "Graph isomorphism (sPyRMSD)",
                    "Top-1 RMSD (Å)":     graph_rmsd,
                    "Best Pose RMSD (Å)": graph_best_rmsd,
                    "Best Pose #":        graph_best_pose,
                    "Used as Final":      "✓" if final_method == "Graph isomorphism (sPyRMSD)" else "",
                    "Status":             g_status,
                })

                # Neither method produced a usable value -- nothing more to
                # add; both rows above already carry the failure reasons.

        if redock_rows:
            df_redock = pd.DataFrame(redock_rows)

            # Summary metrics, the bar chart, and the PDF export all count
            # each COMPLEX once -- using only the row marked "Used as Final"
            # -- even though the table itself shows one row per method for
            # transparency. Without this filter, a complex where both
            # methods succeeded would be double-counted.
            df_final = df_redock[df_redock["Used as Final"] == "✓"].copy()

            # Save only the final rows with a real RMSD value for the PDF report
            st.session_state["_pdf_redock_rows"] = [
                r for r in df_final.to_dict("records")
                if r.get("Top-1 RMSD (Å)") is not None
            ]

            valid_top1 = df_final["Top-1 RMSD (Å)"].dropna()
            valid_best = df_final["Best Pose RMSD (Å)"].dropna()

            rd1, rd2, rd3, rd4 = st.columns(4)
            rd1.metric("Pairs evaluated",       pairs_evaluated)
            rd2.metric("✅ Top-1 Pass (< 2 Å)", int((valid_top1 < 2.0).sum()))
            rd3.metric("🔴 Top-1 Fail (≥ 3 Å)", int((valid_top1 >= 3.0).sum()),
                       delta=int(-(valid_top1 >= 3.0).sum()) if (valid_top1 >= 3.0).any() else None,
                       delta_color="inverse")
            rd4.metric("✅ Best-Pose Pass (< 2 Å)", int((valid_best < 2.0).sum()))

            st.dataframe(df_redock, use_container_width=True, hide_index=True)

            if not valid_top1.empty:
                df_rp = df_final.dropna(subset=["Top-1 RMSD (Å)"]).copy()
                df_rp["Label"] = df_rp["Protein"] + "\n" + df_rp["Ligand"]
                rc = [
                    "#1a9e75" if v < 2.0 else "#f0a500" if v < 3.0 else "#d63030"
                    for v in df_rp["Top-1 RMSD (Å)"]
                ]
                fig3, ax3 = plt.subplots(figsize=(8, max(3, len(df_rp) * 0.45)))
                ax3.barh(df_rp["Label"], df_rp["Top-1 RMSD (Å)"],
                         color=rc, alpha=0.82)
                ax3.axvline(2.0, color="#1a9e75", linestyle="--",
                            linewidth=1.2, label="2 Å pass threshold")
                ax3.axvline(3.0, color="#f0a500", linestyle="--",
                            linewidth=1.2, label="3 Å fail threshold")
                ax3.set_xlabel("Top-1 RMSD (Å)", fontsize=9)
                ax3.set_title("Re-Docking RMSD vs Reference Ligand (Top-1 pose)",
                              fontsize=10, fontweight="bold")
                ax3.legend(fontsize=8)
                ax3.invert_yaxis()
                fig3.tight_layout()
                st.pyplot(fig3)
                plt.close(fig3)
    elif uploaded_refs and not docked_pdbqts:
        st.warning("No docked PDBQT files found. Run Docking first.")
    else:
        st.info("Upload one or more reference PDB files above to run re-docking validation.")


    st.divider()

    # ── PDF Report Generator ──────────────────────────────────────────────
    st.markdown("#### Validation Report (PDF)")
    st.caption(
        "Generates a downloadable PDF summarising all validation results: "
        "re-docking RMSD and docking results table."
    )

    if st.button("Generate PDF Report", use_container_width=True, key="gen_pdf"):
        from matplotlib.backends.backend_pdf import PdfPages
        import datetime

        pdf_buf = io.BytesIO()
        with PdfPages(pdf_buf) as pdf:

            # ── Page 1: Title page ────────────────────────────────────────
            fig_t, ax_t = plt.subplots(figsize=(11, 8.5))
            ax_t.axis("off")
            ax_t.text(0.5, 0.72, "IBDock Validation Report",
                      ha="center", va="center", fontsize=26, fontweight="bold",
                      transform=ax_t.transAxes)
            ax_t.text(0.5, 0.62, f"Project: {project_dir.name}",
                      ha="center", va="center", fontsize=14,
                      transform=ax_t.transAxes, color="#444")
            ax_t.text(0.5, 0.54,
                      f"Generated: {datetime.datetime.now().strftime('%Y-%m-%d  %H:%M')}",
                      ha="center", va="center", fontsize=11,
                      transform=ax_t.transAxes, color="#666")
            ax_t.text(0.5, 0.38,
                      "Contents\n"
                      "  1 · Docking Results Summary\n"
                      "  2 · Re-Docking RMSD vs Reference Ligand",
                      ha="center", va="center", fontsize=12,
                      transform=ax_t.transAxes, color="#333",
                      linespacing=1.9)
            fig_t.patch.set_facecolor("#f7f9fc")
            pdf.savefig(fig_t, bbox_inches="tight")
            plt.close(fig_t)

            # ── Page 2: Docking Results Table ─────────────────────────────
            # Reload fresh from disk every time — df_results may be None at button click
            _pdf_csv_path = result_dir / "docking_results.csv"
            _pdf_df = None
            if _pdf_csv_path.exists():
                try:
                    _pdf_df = pd.read_csv(_pdf_csv_path)
                    _pdf_df = _pdf_df.dropna(how="all")
                except Exception:
                    pass
            if _pdf_df is not None and not _pdf_df.empty:
                df_results = _pdf_df  # update for any later use in PDF
            if df_results is not None and not df_results.empty:
                df_show = df_results.copy()
                # Round floats
                for col in df_show.select_dtypes("float").columns:
                    df_show[col] = df_show[col].round(3)
                df_show = df_show.dropna(how="all")

                n_rows = len(df_show)
                fig_h  = max(3.5, 0.38 * n_rows + 1.5)
                fig_r, ax_r = plt.subplots(figsize=(14, fig_h))
                ax_r.axis("off")
                ax_r.set_title("Docking Results Summary",
                               fontsize=14, fontweight="bold", pad=12)

                cols  = list(df_show.columns)
                rows  = df_show.values.tolist()
                tbl   = ax_r.table(cellText=rows, colLabels=cols,
                                   loc="center", cellLoc="center")
                tbl.auto_set_font_size(False)
                tbl.set_fontsize(7.5)
                tbl.auto_set_column_width(list(range(len(cols))))

                # Header styling
                for j in range(len(cols)):
                    tbl[0, j].set_facecolor("#0066cc")
                    tbl[0, j].set_text_props(color="white", fontweight="bold")
                # Alternating rows
                for i in range(1, n_rows + 1):
                    clr = "#f0f5ff" if i % 2 == 0 else "white"
                    for j in range(len(cols)):
                        tbl[i, j].set_facecolor(clr)

                fig_r.tight_layout()
                pdf.savefig(fig_r, bbox_inches="tight")
                plt.close(fig_r)

            # ── Page 3: Re-docking RMSD table (if results exist) ─────────
            # Re-run the re-docking data if available in session state
            _rd_data = st.session_state.get("_pdf_redock_rows", [])
            if _rd_data:
                df_rd = pd.DataFrame(_rd_data)
                n_rd  = len(df_rd)
                fig_rd_h = max(3.5, 0.38 * n_rd + 1.5)
                fig_rd, ax_rd = plt.subplots(figsize=(14, fig_rd_h))
                ax_rd.axis("off")
                ax_rd.set_title("Re-Docking RMSD vs Reference Ligand",
                                fontsize=14, fontweight="bold", pad=12)
                rd_cols = list(df_rd.columns)
                rd_rows = df_rd.round(4).values.tolist()
                tbl_rd  = ax_rd.table(cellText=rd_rows, colLabels=rd_cols,
                                      loc="center", cellLoc="center")
                tbl_rd.auto_set_font_size(False)
                tbl_rd.set_fontsize(7.5)
                tbl_rd.auto_set_column_width(list(range(len(rd_cols))))
                for j in range(len(rd_cols)):
                    tbl_rd[0, j].set_facecolor("#0066cc")
                    tbl_rd[0, j].set_text_props(color="white", fontweight="bold")
                for i in range(1, n_rd + 1):
                    clr = "#f0f5ff" if i % 2 == 0 else "white"
                    for j in range(len(rd_cols)):
                        tbl_rd[i, j].set_facecolor(clr)
                fig_rd.tight_layout()
                pdf.savefig(fig_rd, bbox_inches="tight")
                plt.close(fig_rd)

            # ── PDF metadata ─────────────────────────────────────────────
            d = pdf.infodict()
            d["Title"]   = "IBDock Validation Report"
            d["Author"]  = "IBDock"
            d["Subject"] = "Molecular Docking Validation"

        pdf_buf.seek(0)
        fname = f"IBDock_validation_{datetime.datetime.now().strftime('%Y%m%d_%H%M')}.pdf"
        st.download_button(
            "⬇ Download PDF Report",
            data=pdf_buf.getvalue(),
            file_name=fname,
            mime="application/pdf",
            use_container_width=True,
            key="dl_pdf",
        )
        st.success(f"✅ PDF ready — click above to download **{fname}**")


# ══════════════════════════════════════════════════════════════════════════════


# TAB 7 — ABOUT & CITATION
# ══════════════════════════════════════════════════════════════════════════════
with tab7:
    st.markdown("### About IBDock")
    st.caption("Citation information, tool credits, and version details.")

    c_about, c_cite = st.columns([1, 1], gap="large")

    with c_about:
        st.markdown("""
<div class="info-card">
<h4>What is IBDock?</h4>
<p>
IBDock is a browser-based graphical interface for <b>batch molecular docking</b>
using AutoDock Vina. It integrates the full virtual screening pipeline —
protein preparation, ligand preparation, grid box generation, parallel docking,
and results analysis — into a single zero-command-line application.
</p>
</div>

<div class="info-card">
<h4>Designed for</h4>
<p>
Computational and medicinal chemists and structural biologists who need
to screen compound libraries against one or more protein targets without
writing scripts or managing command-line flags.
</p>
</div>

<div class="info-card">
<h4>Platform</h4>
<p>
Windows &middot; macOS &middot; Linux<br>
The packaged Windows build bundles MGLTools, Vina, Open Babel, P2Rank, and fpocket as native Windows tools.
Pocket detection uses the bundled native P2Rank and fpocket executables.
</p>
</div>
        """, unsafe_allow_html=True)

    with c_cite:
        st.markdown("""
<div class="info-card">
<h4>How to cite</h4>
<p>If you use IBDock in published research, please cite the following tools:</p>
</div>
        """, unsafe_allow_html=True)

    with st.expander("AutoDock Vina", expanded=False):
        st.code(
            "Eberhardt J, Santos-Martins D, Tillack AF, Forli S. (2021). "
            "AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, "
            "and Python Bindings. J Chem Inf Model. 61(8):3891–3898. "
            "DOI: 10.1021/acs.jcim.1c00203",
            language=None,
        )
        st.code(
            "Trott O, Olson AJ. (2010). AutoDock Vina: Improving the speed "
            "and accuracy of docking with a new scoring function, efficient "
            "optimization, and multithreading. J Comput Chem. 31(2):455–461. "
            "DOI: 10.1002/jcc.21334",
            language=None,
        )

    with st.expander("MGLTools / AutoDockTools"):
        st.code(
            "Morris GM, Huey R, Lindstrom W, Sanner MF, Belew RK, Goodsell DS, Olson AJ. (2009). "
            "AutoDock4 and AutoDockTools4: Automated docking with selective receptor flexibility. "
            "J Comput Chem. 30(16):2785–2791. DOI: 10.1002/jcc.21256",
            language=None,
        )

    with st.expander("Open Babel"):
        st.code(
            "O'Boyle NM, Banck M, James CA, Morley C, Vandermeersch T, Hutchison GR. (2011). "
            "Open Babel: An open chemical toolbox. J Cheminform. 3:33. "
            "DOI: 10.1186/1758-2946-3-33",
            language=None,
        )

    with st.expander("P2Rank (pocket prediction)"):
        st.code(
            "Krivák R, Hoksza D. (2018). P2Rank: machine learning based tool "
            "for rapid and accurate prediction of ligand binding sites from "
            "protein structure. J Cheminform. 10:39. DOI: 10.1186/s13321-018-0285-8",
            language=None,
        )

    with st.expander("fpocket (pocket prediction)"):
        st.code(
            "Le Guilloux V, Schmidtke P, Tuffery P. (2009). Fpocket: An open "
            "source platform for ligand pocket detection. BMC Bioinformatics. 10:168. "
            "DOI: 10.1186/1471-2105-10-168",
            language=None,
        )

    st.divider()

    st.markdown("""
<div style="display:grid; grid-template-columns:1fr 1fr 1fr; gap:0.8rem; margin-top:0.5rem;">
  <div class="info-card">
        <h4>Software Stack</h4>
        <p>
          Python 3.x &middot; Streamlit &middot; NumPy &middot; Pandas<br>
          Matplotlib &middot; Seaborn &middot; py3Dmol
        </p>
  </div>
  <div class="info-card">
        <h4>License</h4>
        <p>
         IBDock is open-source software.<br>
          AutoDock Vina, MGLTools, Open Babel, P2Rank, and fpocket
          are subject to their respective licences.
        </p>
  </div>
  <div class="info-card">
        <h4>Disclaimer</h4>
        <p>
          Docking scores are approximations. Results should be interpreted
          alongside experimental validation and expert knowledge.
          Not intended for clinical decision-making.
        </p>
  </div>
</div>
    """, unsafe_allow_html=True)
