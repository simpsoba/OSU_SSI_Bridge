#!/usr/bin/env python3
"""
Goals
-----
Canonical paths and scale constants for lab post-process scripts.

  Durable archive (Shared Drive):
    …/2026-08-21-OSU-SSI-Bridge/opensees_data/
  Junction (repo):
    OSU_SSI_BRIDGE_DATA/  →  that Shared Drive opensees_data folder
  Local workspace (gitignored):
    OSU_SSI_BRIDGE_DATA_LOCAL/opensees_data/   # mirror + mat_extract
    OSU_SSI_BRIDGE_DATA_LOCAL/plots/           # all post-process PNGs
  Git-tracked:
    plot/opensees_data/mat_run_map.json
    plot/opensees_data/TestMatrix_lab_runs.csv

Units
-----
  λ = CYLINDER_LENGTH_SCALE (2.4). Froude time scale = √λ.
  Lab mat Time = real time at **model** scale.
  OpenSees numerical t = **prototype** scale.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

from paths import HERE

# ------------------------------------------------------------
# layout
# ------------------------------------------------------------
REPO = HERE.parent
DRIVE_ROOT = REPO / "OSU_SSI_BRIDGE_DATA"
LOCAL_ROOT = REPO / "OSU_SSI_BRIDGE_DATA_LOCAL"
LOCAL_OPENSEES_DATA = LOCAL_ROOT / "opensees_data"
LOCAL_PLOTS = LOCAL_ROOT / "plots"
LEGACY_PLOTS_ROOT = REPO / "OSU_SSI_PLOTS"
LEGACY_SPACE_ROOT = LOCAL_ROOT / "opensees data"  # old local mirror name

SHARED_DRIVE_ARCHIVE = Path(
    r"G:\Shared drives\Simpson team\Test Data\2026-08-21-OSU-SSI-Bridge"
)
SHARED_DRIVE_OPENSEES_DATA = SHARED_DRIVE_ARCHIVE / "opensees_data"

OPENSEES_DATA_DIR = HERE / "opensees_data"
MAT_RUN_MAP_PATH = OPENSEES_DATA_DIR / "mat_run_map.json"
LAB_RUNS_CSV = OPENSEES_DATA_DIR / "TestMatrix_lab_runs.csv"
LAB_RUNS_CSV_LOCAL = LOCAL_ROOT / "TestMatrix_lab_runs.csv"
MAT_EXTRACT_DIR = LOCAL_OPENSEES_DATA / "mat_extract"
MANIFEST_PATH = LOCAL_OPENSEES_DATA / "backup_manifest.json"
LOCAL_MAT_RUN_MAP_PATH = LOCAL_OPENSEES_DATA / "mat_run_map.json"

# Top-level .mat keys to keep: container `data` plus *OS hybrid channels.
MAT_KEYS = ("data", "tarSigOS", "comSigOS", "meaSigOS", "stateOS")

# Cylinder length scale λ (prototype / model). Froude time scale = √λ.
CYLINDER_LENGTH_SCALE = 2.4
TIME_SCALE_FROUDE = CYLINDER_LENGTH_SCALE**0.5  # ≈ 1.549


# ------------------------------------------------------------
# resolve roots
# ------------------------------------------------------------


def plots_root() -> Path:
    """
    Canonical PNG output dir (LOCAL/plots), with legacy OSU_SSI_PLOTS fallback.

    Args:    none
    Returns: Path (created if needed)
    """
    if LOCAL_PLOTS.is_dir():
        return LOCAL_PLOTS
    if LEGACY_PLOTS_ROOT.is_dir():
        return LEGACY_PLOTS_ROOT
    LOCAL_PLOTS.mkdir(parents=True, exist_ok=True)
    return LOCAL_PLOTS


def lab_runs_csv_path() -> Path:
    """
    Prefer git-tracked CSV; fall back to LOCAL working copy.

    Args:    none
    Returns: Path (may not exist yet)
    """
    if LAB_RUNS_CSV.is_file():
        return LAB_RUNS_CSV
    return LAB_RUNS_CSV_LOCAL


def _has_run_dirs(folder: Path) -> bool:
    """True if folder contains at least one r±… run directory."""
    if not folder.is_dir():
        return False
    for child in folder.iterdir():
        if child.is_dir() and child.name.startswith(("r+", "r-")):
            return True
    return False


def resolve_opensees_data() -> Path | None:
    """
    Folder that holds run dirs.

    Prefer Shared Drive archive, then junction, then LOCAL mirror.
    Also accepts legacy layouts (opensees data/ under junction or shortcut-targets).

    Args:    none
    Returns: Path or None
    """
    if _has_run_dirs(SHARED_DRIVE_OPENSEES_DATA):
        return SHARED_DRIVE_OPENSEES_DATA
    if _has_run_dirs(DRIVE_ROOT):
        return DRIVE_ROOT

    direct = DRIVE_ROOT / "opensees data"
    if _has_run_dirs(direct):
        return direct

    shortcut_root = Path(r"G:\.shortcut-targets-by-id")
    if shortcut_root.is_dir():
        for path in sorted(shortcut_root.glob("*/opensees data")):
            if _has_run_dirs(path):
                return path

    if _has_run_dirs(LOCAL_OPENSEES_DATA):
        return LOCAL_OPENSEES_DATA
    if LOCAL_OPENSEES_DATA.is_dir():
        return LOCAL_OPENSEES_DATA
    return None


def resolve_simulink_dir(root: Path | None = None) -> Path | None:
    """
    Simulink/ under an opensees_data root.

    Args:    root  default resolve_opensees_data()
    Returns: Path or None
    """
    root = root or resolve_opensees_data()
    if root is None:
        return None
    path = root / "Simulink"
    return path if path.is_dir() else None


# ------------------------------------------------------------
# legacy layout migrations (one-shot moves)
# ------------------------------------------------------------


def _merge_tree(src: Path, dst: Path) -> bool:
    """
    Merge src into dst (move missing files/dirs). Does not overwrite.

    Args:    src, dst
    Returns: True if anything moved
    """
    moved = False
    dst.mkdir(parents=True, exist_ok=True)
    for item in src.iterdir():
        target = dst / item.name
        if item.is_dir():
            if target.is_dir():
                moved = _merge_tree(item, target) or moved
            elif not target.exists():
                shutil.move(str(item), str(target))
                moved = True
        elif not target.exists():
            shutil.move(str(item), str(target))
            moved = True
    return moved


def migrate_legacy_local_layout() -> bool:
    """
    Move legacy LOCAL `opensees data/` into `opensees_data/`.

    Args:    none
    Returns: True if a migration happened
    """
    space = LEGACY_SPACE_ROOT
    local = LOCAL_OPENSEES_DATA
    if not space.is_dir():
        return False
    if not local.is_dir():
        shutil.move(str(space), str(local))
        return True
    moved = _merge_tree(space, local)
    try:
        if space.is_dir() and not any(space.iterdir()):
            space.rmdir()
            moved = True
    except OSError:
        pass
    return moved


def migrate_legacy_plots() -> bool:
    """
    Move repo-root OSU_SSI_PLOTS/ into LOCAL/plots/.

    Args:    none
    Returns: True if a migration happened
    """
    if not LEGACY_PLOTS_ROOT.is_dir():
        return False
    if not LOCAL_PLOTS.is_dir():
        shutil.move(str(LEGACY_PLOTS_ROOT), str(LOCAL_PLOTS))
        return True
    moved = _merge_tree(LEGACY_PLOTS_ROOT, LOCAL_PLOTS)
    try:
        if LEGACY_PLOTS_ROOT.is_dir() and not any(LEGACY_PLOTS_ROOT.iterdir()):
            LEGACY_PLOTS_ROOT.rmdir()
            moved = True
    except OSError:
        pass
    return moved


# ------------------------------------------------------------
# mat_run_map + MATLAB helpers
# ------------------------------------------------------------


def load_mat_run_map(path: Path | None = None) -> dict:
    """
    Load mat ↔ run dictionary.

    Args:    path  default MAT_RUN_MAP_PATH
    Returns: dict with mats / runs_without_mat (empty if missing)
    """
    map_path = path or MAT_RUN_MAP_PATH
    if not map_path.is_file():
        return {"mats": {}, "runs_without_mat": []}
    return json.loads(map_path.read_text(encoding="utf-8"))


def save_mat_run_map(data: dict, path: Path | None = None) -> None:
    """
    Write mat_run_map.json (pretty JSON + trailing newline).

    Args:    data, path  default MAT_RUN_MAP_PATH
    Returns: none
    """
    map_path = path or MAT_RUN_MAP_PATH
    map_path.parent.mkdir(parents=True, exist_ok=True)
    map_path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def cell_str(value) -> str:
    """
    Unwrap a MATLAB cell / 0-d ndarray to a plain str.

    Args:    value  scipy.io cell contents
    Returns: str
    """
    import numpy as np

    s = value
    while isinstance(s, np.ndarray) and s.size == 1:
        s = s.flat[0]
    if isinstance(s, bytes):
        return s.decode("utf-8", errors="replace")
    return str(s)


def signal_names(block) -> list[str]:
    """
    signalNames field from an OpenFresco *OS MATLAB struct.

    Args:    block  struct with .signalNames
    Returns: list of str
    """
    return [cell_str(n) for n in block.signalNames.flatten()]


def time_column(names: list[str]) -> int:
    """
    Index of the Time column in a signal-name list.

    Args:    names
    Returns: column index
    Raises:  KeyError if no Time column
    """
    for i, name in enumerate(names):
        if name.strip().lower() == "time":
            return i
    raise KeyError(f"no Time column in {names}")
