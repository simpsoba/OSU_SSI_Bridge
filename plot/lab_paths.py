#!/usr/bin/env python3
"""Paths and helpers for lab dumps (Drive) and the local backup mirror.

  Durable archive (Shared Drive):
    G:\\Shared drives\\Simpson team\\Test Data\\2026-08-21-OSU-SSI-Bridge\\opensees_data\\
  Junction (repo):
    OSU_SSI_BRIDGE_DATA/  →  that Shared Drive opensees_data folder
  Local workspace (gitignored):
    OSU_SSI_BRIDGE_DATA_LOCAL/opensees_data/   # working mirror + mat_extract
    OSU_SSI_BRIDGE_DATA_LOCAL/plots/           # all post-process PNGs
  Git-tracked:
    plot/opensees_data/mat_run_map.json
    plot/opensees_data/TestMatrix_lab_runs.csv
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

from paths import HERE

REPO = HERE.parent
DRIVE_ROOT = REPO / "OSU_SSI_BRIDGE_DATA"
LOCAL_ROOT = REPO / "OSU_SSI_BRIDGE_DATA_LOCAL"
LOCAL_OPENSEES_DATA = LOCAL_ROOT / "opensees_data"
LOCAL_PLOTS = LOCAL_ROOT / "plots"
LEGACY_PLOTS_ROOT = REPO / "OSU_SSI_PLOTS"
LEGACY_SPACE_ROOT = LOCAL_ROOT / "opensees data"  # old local mirror name

# Durable Simpson Shared Drive campaign archive (run dirs + Simulink).
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
# Lab Time in the mats = real time at model scale; OpenSees t = numerical time at prototype.
CYLINDER_LENGTH_SCALE = 2.4
TIME_SCALE_FROUDE = CYLINDER_LENGTH_SCALE**0.5  # ≈ 1.549


def plots_root() -> Path:
    """Canonical plot output dir (LOCAL/plots), with legacy OSU_SSI_PLOTS fallback."""
    if LOCAL_PLOTS.is_dir():
        return LOCAL_PLOTS
    if LEGACY_PLOTS_ROOT.is_dir():
        return LEGACY_PLOTS_ROOT
    LOCAL_PLOTS.mkdir(parents=True, exist_ok=True)
    return LOCAL_PLOTS


def lab_runs_csv_path() -> Path:
    """Prefer git-tracked CSV; fall back to LOCAL working copy."""
    if LAB_RUNS_CSV.is_file():
        return LAB_RUNS_CSV
    return LAB_RUNS_CSV_LOCAL


def _merge_tree(src: Path, dst: Path) -> bool:
    """Merge src into dst (files/dirs). Returns True if anything moved."""
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
    """Move legacy `opensees data/` local tree into `opensees_data/`."""
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
    """Move repo-root OSU_SSI_PLOTS/ into LOCAL/plots/."""
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


def _has_run_dirs(p: Path) -> bool:
    if not p.is_dir():
        return False
    for c in p.iterdir():
        if c.is_dir() and c.name.startswith(("r+", "r-")):
            return True
    return False


def resolve_opensees_data() -> Path | None:
    """Folder that holds run dirs. Prefer Shared Drive archive, then junction, then LOCAL."""
    if _has_run_dirs(SHARED_DRIVE_OPENSEES_DATA):
        return SHARED_DRIVE_OPENSEES_DATA
    if _has_run_dirs(DRIVE_ROOT):
        return DRIVE_ROOT
    # Legacy layouts under junction / My Drive / shortcut
    direct = DRIVE_ROOT / "opensees data"
    if _has_run_dirs(direct):
        return direct
    shortcut_root = Path(r"G:\.shortcut-targets-by-id")
    if shortcut_root.is_dir():
        for p in sorted(shortcut_root.glob("*/opensees data")):
            if _has_run_dirs(p):
                return p
    if _has_run_dirs(LOCAL_OPENSEES_DATA):
        return LOCAL_OPENSEES_DATA
    if LOCAL_OPENSEES_DATA.is_dir():
        return LOCAL_OPENSEES_DATA
    return None


def resolve_simulink_dir(root: Path | None = None) -> Path | None:
    root = root or resolve_opensees_data()
    if root is None:
        return None
    p = root / "Simulink"
    return p if p.is_dir() else None


def load_mat_run_map(path: Path | None = None) -> dict:
    p = path or MAT_RUN_MAP_PATH
    if not p.is_file():
        return {"mats": {}, "runs_without_mat": []}
    return json.loads(p.read_text(encoding="utf-8"))


def save_mat_run_map(data: dict, path: Path | None = None) -> None:
    p = path or MAT_RUN_MAP_PATH
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def cell_str(x) -> str:
    """Unwrap MATLAB cell / nested ndarray to a plain str."""
    import numpy as np

    s = x
    while isinstance(s, np.ndarray) and s.size == 1:
        s = s.flat[0]
    if isinstance(s, bytes):
        return s.decode("utf-8", errors="replace")
    return str(s)


def signal_names(block) -> list[str]:
    return [cell_str(n) for n in block.signalNames.flatten()]


def time_column(names: list[str]) -> int:
    for i, n in enumerate(names):
        if n.strip().lower() == "time":
            return i
    raise KeyError(f"no Time column in {names}")
