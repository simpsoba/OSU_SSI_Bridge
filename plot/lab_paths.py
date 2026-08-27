#!/usr/bin/env python3
"""
Goals
-----
Canonical paths and scale constants for lab post-process scripts.

  Live lab upload (Drive Desktop shortcut; ingest only):
    G:/.shortcut-targets-by-id/…/opensees data/
  Local working mirror (all post-process reads/writes dumps here):
    OSU_SSI_BRIDGE_DATA_LOCAL/opensees_data/
    OSU_SSI_BRIDGE_DATA/  →  junction to that local mirror
  Safety archive (copy from LOCAL after ingest; not used for plotting):
    G:/Shared drives/Simpson team/Test Data/2026-OSU-SSI-Bridge/
      2026-08-19/opensees_data/   # Wed dumps + Simulink/
      2026-08-21/opensees_data/   # Fri dumps + Simulink/
  Plots (always local):
    OSU_SSI_BRIDGE_DATA_LOCAL/plots/
      runs/<dump>/eq/     # PlotEQ* (OpenSees)
      runs/<dump>/os/     # PlotMatOS (mapped mat)
      mats/<stem>/os/     # PlotMatOS (unmapped mat)
      compare/<group>/         # (pairs/ only for pier compares)
      compare/<group>/pairs/   # interim ref vs other (PlotEQComparePairs)
      compare/stateos/         # campaign typeConv3 stacked bars (PlotStateOSBars)
  Git-tracked:
    plot/lab/TestMatrix_lab_runs.csv   # as-run index (Test, dump, mat, …)
    plot/lab/mat_run_map.json          # orphan / pending mats only
    plot/lab/LAB_RUN_MAP.md            # narrative (not a second table)
    plot/lab/1_Monopile_matrix.xlsx    # lab Run Log / schedule

Units
-----
  λ = CYLINDER_LENGTH_SCALE (2.4). Froude time scale = √λ.
  Lab mat Time = real time at **model** scale.
  OpenSees numerical t = **prototype** scale.
  Displacement plots: **mm** (prototype), except deformed-shape x/y axes (m).
"""

from __future__ import annotations

import json
import re
import shutil
from pathlib import Path

from paths import HERE

# ------------------------------------------------------------
# layout
# ------------------------------------------------------------
REPO = HERE.parent
# Working dump root for scripts; junction should target LOCAL_OPENSEES_DATA.
DRIVE_ROOT = REPO / "OSU_SSI_BRIDGE_DATA"
LOCAL_ROOT = REPO / "OSU_SSI_BRIDGE_DATA_LOCAL"
LOCAL_OPENSEES_DATA = LOCAL_ROOT / "opensees_data"
LOCAL_PLOTS = LOCAL_ROOT / "plots"
LEGACY_PLOTS_ROOT = REPO / "OSU_SSI_PLOTS"
LEGACY_SPACE_ROOT = LOCAL_ROOT / "opensees data"  # old local mirror name

# Live flume upload folder (Google Drive Desktop “shortcut-targets”).
LAB_INGEST_SOURCE = Path(
    r"G:\.shortcut-targets-by-id\1s4OLDygoKnrdIonCT7Ht3w9ywrrtK_pz\opensees data"
)

SHARED_DRIVE_ARCHIVE = Path(
    r"G:\Shared drives\Simpson team\Test Data\2026-OSU-SSI-Bridge"
)
# Legacy flat layout (pre date subfolders); still accepted as a fallback.
_SHARED_DRIVE_ARCHIVE_FLAT = SHARED_DRIVE_ARCHIVE / "opensees_data"

LAB_DIR = HERE / "lab"  # git-tracked campaign index (not dump blobs)
MAT_RUN_MAP_PATH = LAB_DIR / "mat_run_map.json"
LAB_RUNS_CSV = LAB_DIR / "TestMatrix_lab_runs.csv"
LAB_RUNS_CSV_LOCAL = LOCAL_ROOT / "TestMatrix_lab_runs.csv"
MAT_EXTRACT_DIR = LOCAL_OPENSEES_DATA / "mat_extract"
MANIFEST_PATH = LOCAL_OPENSEES_DATA / "backup_manifest.json"
LOCAL_MAT_RUN_MAP_PATH = LOCAL_OPENSEES_DATA / "mat_run_map.json"

# Do not push derived LOCAL-only trees into the Shared Drive archive.
ARCHIVE_SKIP_TOP = frozenset({"mat_extract", "plots"})

# Campaign day folder under SHARED_DRIVE_ARCHIVE (YYYY-MM-DD).
_ARCHIVE_DAY_RE = re.compile(r"^(\d{4}-\d{2}-\d{2})$")

# Top-level .mat keys to keep: container `data` plus *OS hybrid channels.
MAT_KEYS = ("data", "tarSigOS", "comSigOS", "meaSigOS", "stateOS")

# Cylinder length scale λ (prototype / model). Froude time scale = √λ.
CYLINDER_LENGTH_SCALE = 2.4
TIME_SCALE_FROUDE = CYLINDER_LENGTH_SCALE**0.5  # ≈ 1.549
M_TO_MM = 1.0e3  # recorder displacement (m) → plot (mm)
# Simulink *OS displacement at model scale (m) → prototype plot (mm).
DISP_M_TO_PROTO_MM = CYLINDER_LENGTH_SCALE * M_TO_MM

# Lab dump folder names: r+01_YYYYMMDD_HHMM_… / r-02_…
_RUN_FOLDER_RE = re.compile(r"^r[+-]?\d+_", re.IGNORECASE)
_OS_PNG_PREFIX = "hist_os_"


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


def run_plot_stem(dump_name: str) -> str:
    """
    plots/runs/<stem>/ folder for a dump: prefer Test ID ``F27``.

    Args:    dump_name  DumpFolder (e.g. r-02_20260821_1601_Storm_Wave)
    Returns: ``F##`` / ``W##`` when mapped, else the dump folder name
    """
    from compare_groups import dump_to_test_id, test_file_slug

    tid = dump_to_test_id().get(dump_name)
    if tid is not None:
        return test_file_slug(tid)
    return dump_name


def run_eq_plots_dir(run_name: str) -> Path:
    """
    OpenSees EQ PNGs for one dump: plots/runs/F##/eq/.

    Args:    run_name  dump folder name (resolved to Test ID when known)
    Returns: Path (not created)
    """
    return plots_root() / "runs" / run_plot_stem(run_name) / "eq"


def run_os_plots_dir(run_name: str) -> Path:
    """
    Simulink *OS PNGs for a mapped dump: plots/runs/F##/os/.

    Args:    run_name  dump folder name
    Returns: Path (not created)
    """
    return plots_root() / "runs" / run_plot_stem(run_name) / "os"


def test_os_plots_dir(test_id: str | int) -> Path:
    """
    Simulink *OS for a Test ID (incl. dry mat-only): plots/runs/F##/os/.

    Args:    test_id  CSV Test cell
    Returns: Path (not created)
    """
    from compare_groups import test_file_slug

    return plots_root() / "runs" / test_file_slug(test_id) / "os"


def mat_os_plots_dir(mat_stem: str) -> Path:
    """
    Simulink *OS PNGs for an unmapped mat: plots/mats/<stem>/os/.

    Args:    mat_stem  .mat file stem (no extension)
    Returns: Path (not created)
    """
    return plots_root() / "mats" / mat_stem / "os"


def compare_plots_dir(group_slug: str) -> Path:
    """
    Overlay PNGs for one physical-model group: plots/compare/<group>/.

    Args:    group_slug  from compare_groups.group_label / dump_to_group
    Returns: Path (not created)
    """
    return plots_root() / "compare" / group_slug


def stateos_compare_plots_dir() -> Path:
    """
    Campaign-wide stateOS compare figures: plots/compare/stateos/.

    Stacked typeConv3 bars (PlotStateOSBars.py), not per-run os/ histories.

    Args:    none
    Returns: Path (not created)
    """
    return plots_root() / "compare" / "stateos"


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
    Working dump folder for all post-process (LOCAL only).

    Never returns Shared Drive or shortcut-targets paths — plotters and
    compare scripts must not read through Drive.

    Args:    none
    Returns: Path or None
    """
    if _has_run_dirs(LOCAL_OPENSEES_DATA):
        return LOCAL_OPENSEES_DATA
    # Junction retargeted to LOCAL: same tree under DRIVE_ROOT.
    if _has_run_dirs(DRIVE_ROOT):
        try:
            if DRIVE_ROOT.resolve() == LOCAL_OPENSEES_DATA.resolve():
                return LOCAL_OPENSEES_DATA
        except OSError:
            pass
    if LOCAL_OPENSEES_DATA.is_dir():
        return LOCAL_OPENSEES_DATA
    return None


def resolve_lab_ingest_source() -> Path | None:
    """
    Live lab upload folder (shortcut-targets). SyncLabBackup reads here only.

    Args:    none
    Returns: Path or None
    """
    if _has_run_dirs(LAB_INGEST_SOURCE):
        return LAB_INGEST_SOURCE
    shortcut_root = Path(r"G:\.shortcut-targets-by-id")
    if shortcut_root.is_dir():
        for path in sorted(shortcut_root.glob("*/opensees data")):
            if _has_run_dirs(path):
                return path
    return None


def archive_day_from_name(name: str) -> str:
    """
    Map a dump folder or mat file name to archive day ``YYYY-MM-DD``.

    Args:    name  e.g. ``0819_GusBridge_…``, ``r-04_20260821_1010_…``,
                   ``0821_GusBridge_….mat``
    Returns: day folder name (default ``2026-08-21`` when unclear)
    """
    s = (name or "").strip()
    if s.startswith("0819") or "20260819" in s:
        return "2026-08-19"
    if s.startswith("0820") or "20260820" in s:
        return "2026-08-20"
    if s.startswith("0821") or "20260821" in s:
        return "2026-08-21"
    m = re.search(r"_(\d{4})(\d{2})(\d{2})_", s)
    if m:
        return f"{m.group(1)}-{m.group(2)}-{m.group(3)}"
    return "2026-08-21"


def archive_opensees_data_for_day(day: str) -> Path:
    """
    ``…/2026-OSU-SSI-Bridge/<day>/opensees_data``.

    Args:    day  ``YYYY-MM-DD``
    Returns: Path (may not exist yet)
    """
    return SHARED_DRIVE_ARCHIVE / day / "opensees_data"


def resolve_archive_opensees_data() -> Path:
    """
    Default Shared Drive archive day root (Fri campaign).

    Prefer date subfolders. Legacy flat ``…/opensees_data`` is only a
    fallback if no day folders exist yet.

    Args:    none
    Returns: archive opensees_data path (may not exist yet)
    """
    fri = archive_opensees_data_for_day("2026-08-21")
    if fri.is_dir() or any(
        (SHARED_DRIVE_ARCHIVE / p).is_dir()
        for p in ("2026-08-19", "2026-08-21")
    ):
        return fri
    if _SHARED_DRIVE_ARCHIVE_FLAT.is_dir():
        return _SHARED_DRIVE_ARCHIVE_FLAT
    return fri


def list_archive_opensees_roots() -> list[Path]:
    """
    All day ``opensees_data`` folders under the campaign archive.

    Args:    none
    Returns: existing day roots (Wed, Fri, …)
    """
    roots: list[Path] = []
    if not SHARED_DRIVE_ARCHIVE.is_dir():
        return roots
    for child in sorted(SHARED_DRIVE_ARCHIVE.iterdir()):
        if child.is_dir() and _ARCHIVE_DAY_RE.match(child.name):
            od = child / "opensees_data"
            if od.is_dir():
                roots.append(od)
    if not roots and _SHARED_DRIVE_ARCHIVE_FLAT.is_dir():
        roots.append(_SHARED_DRIVE_ARCHIVE_FLAT)
    return roots


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


def _move_into(src: Path, dest_dir: Path) -> bool:
    """
    Move src into dest_dir / src.name if the target does not exist.

    Args:    src  file or directory; dest_dir  destination parent
    Returns: True if something was moved
    """
    dest_dir.mkdir(parents=True, exist_ok=True)
    target = dest_dir / src.name
    if target.exists():
        return False
    shutil.move(str(src), str(target))
    return True


def migrate_plots_layout(root: Path | None = None) -> dict[str, int]:
    """
    Reorganize LOCAL/plots into runs/<dump>/{eq,os}, mats/<stem>/os, compare/.

    Idempotent. Moves top-level r±* folders and legacy mat_os/.
    hist_os_*.png under a dump folder go to os/; everything else to eq/.

    Args:    root  default plots_root()
    Returns: counts {runs, mats, os_pngs}
    """
    root = root or plots_root()
    stats = {"runs": 0, "mats": 0, "os_pngs": 0}

    old_mat_os = root / "mat_os"
    if old_mat_os.is_dir():
        for stem_dir in list(old_mat_os.iterdir()):
            if not stem_dir.is_dir():
                continue
            dest = root / "mats" / stem_dir.name / "os"
            dest.mkdir(parents=True, exist_ok=True)
            for item in list(stem_dir.iterdir()):
                if _move_into(item, dest):
                    stats["mats"] += 1
            try:
                if stem_dir.is_dir() and not any(stem_dir.iterdir()):
                    stem_dir.rmdir()
            except OSError:
                pass
        try:
            if old_mat_os.is_dir() and not any(old_mat_os.iterdir()):
                old_mat_os.rmdir()
        except OSError:
            pass

    for child in list(root.iterdir()):
        if not child.is_dir():
            continue
        if child.name in ("runs", "mats", "compare", "mat_os"):
            continue
        if not _RUN_FOLDER_RE.match(child.name):
            continue

        eq_dest = root / "runs" / child.name / "eq"
        os_dest = root / "runs" / child.name / "os"
        eq_dest.mkdir(parents=True, exist_ok=True)

        for item in list(child.iterdir()):
            if item.is_file() and item.name.startswith(_OS_PNG_PREFIX):
                if _move_into(item, os_dest):
                    stats["os_pngs"] += 1
            else:
                if _move_into(item, eq_dest):
                    stats["runs"] += 1
        try:
            if child.is_dir() and not any(child.iterdir()):
                child.rmdir()
        except OSError:
            pass

    return stats


# ------------------------------------------------------------
# lab runs CSV + mat orphan registry
# ------------------------------------------------------------


def load_lab_runs_rows(path: Path | None = None) -> list[dict[str, str]]:
    """
    Read the curated as-run matrix (git-tracked or LOCAL copy).

    Args:    path  default lab_runs_csv_path()
    Returns: list of row dicts (empty if missing)
    """
    import csv

    csv_path = path or lab_runs_csv_path()
    if not csv_path.is_file():
        return []
    with csv_path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def run_to_mat_from_csv(path: Path | None = None) -> dict[str, str]:
    """
    OpenSees dump folder → Simulink .mat file name for paired rows.

    Args:    path  optional CSV override
    Returns: {DumpFolder: MatFile}
    """
    out: dict[str, str] = {}
    for row in load_lab_runs_rows(path):
        dump = (row.get("DumpFolder") or "").strip()
        mat = (row.get("MatFile") or "").strip()
        if dump and mat:
            out[dump] = mat
    return out


def mat_to_run_from_csv(path: Path | None = None) -> dict[str, str]:
    """
    Simulink .mat → OpenSees dump folder for paired rows.

    Args:    path  optional CSV override
    Returns: {MatFile: DumpFolder}
    """
    out: dict[str, str] = {}
    for row in load_lab_runs_rows(path):
        dump = (row.get("DumpFolder") or "").strip()
        mat = (row.get("MatFile") or "").strip()
        if dump and mat:
            out[mat] = dump
    return out


def load_mat_orphans(path: Path | None = None) -> dict:
    """
    Orphan / pending mat registry (slim mat_run_map.json on disk).

    Paired runs live in TestMatrix_lab_runs.csv only.

    Args:    path  default MAT_RUN_MAP_PATH
    Returns: dict with mats_without_dump, pending_*, duplicate_mats
    """
    map_path = path or MAT_RUN_MAP_PATH
    empty = {
        "mats_without_dump": {},
        "pending_mat_upload": [],
        "pending_dump_upload": [],
        "duplicate_mats": {},
    }
    if not map_path.is_file():
        return empty
    data = json.loads(map_path.read_text(encoding="utf-8"))
    # Legacy file: migrate in memory if old "mats" block still present.
    if "mats" in data and "mats_without_dump" not in data:
        orphans = empty.copy()
        for mat_name, info in data.get("mats", {}).items():
            if not info.get("run"):
                orphans["mats_without_dump"][mat_name] = {
                    k: info[k]
                    for k in ("row", "labTrial", "note")
                    if k in info
                }
        orphans["pending_mat_upload"] = list(data.get("pending_mat_upload", []))
        orphans["pending_dump_upload"] = list(data.get("pending_dump_upload", []))
        return orphans
    out = empty.copy()
    out["mats_without_dump"] = dict(data.get("mats_without_dump", {}))
    out["pending_mat_upload"] = list(data.get("pending_mat_upload", []))
    out["pending_dump_upload"] = list(data.get("pending_dump_upload", []))
    out["duplicate_mats"] = dict(data.get("duplicate_mats", {}))
    return out


def _optional_int(text: str) -> int | None:
    s = (text or "").strip()
    if not s:
        return None
    return int(s)


def build_mat_run_catalog(path: Path | None = None) -> dict:
    """
    Combined mat ↔ run view for plotters (CSV pairs + orphan registry).

    Args:    path  optional CSV override for paired rows
    Returns: dict with mats, runs_without_mat, pending_*, duplicate_mats
    """
    orphans = load_mat_orphans()
    mats: dict[str, dict] = {}
    runs_without_mat: list[str] = []

    for row in load_lab_runs_rows(path):
        dump = (row.get("DumpFolder") or "").strip()
        if not dump:
            continue
        mat = (row.get("MatFile") or "").strip()
        if mat:
            mats[mat] = {
                "run": dump,
                "row": (row.get("Run") or "").strip() or None,
                "labTrial": _optional_int(row.get("LabTrial", "")),
                "note": (row.get("Note") or "").strip(),
                "test": (row.get("Test") or "").strip() or None,
            }
        else:
            runs_without_mat.append(dump)

    for mat_name, info in orphans.get("mats_without_dump", {}).items():
        mats[mat_name] = {
            "run": None,
            "row": info.get("row"),
            "labTrial": info.get("labTrial"),
            "note": (info.get("note") or "").strip(),
            "test": None,
        }

    return {
        "mats": mats,
        "runs_without_mat": runs_without_mat,
        "pending_mat_upload": list(orphans.get("pending_mat_upload", [])),
        "pending_dump_upload": list(orphans.get("pending_dump_upload", [])),
        "duplicate_mats": dict(orphans.get("duplicate_mats", {})),
    }


def all_mat_names_for_plot() -> list[str]:
    """
    Mat files to consider for PlotMatOS (CSV + orphans; skip duplicate aliases).

    Args:    none
    Returns: sorted unique .mat names
    """
    catalog = build_mat_run_catalog()
    skip = set(catalog.get("duplicate_mats", {}).keys())
    names = [m for m in catalog.get("mats", {}) if m not in skip]
    return sorted(names)


# ------------------------------------------------------------
# mat_run_map + MATLAB helpers (legacy name: build_mat_run_catalog)
# ------------------------------------------------------------


def load_mat_run_map(path: Path | None = None) -> dict:
    """
    Mat ↔ run catalog for plot scripts (built from CSV + orphan JSON).

    Args:    path  ignored (kept for call-site compatibility)
    Returns: dict with mats / runs_without_mat / pending_* / duplicate_mats
    """
    _ = path
    return build_mat_run_catalog()


def save_mat_run_map(data: dict, path: Path | None = None) -> None:
    """
    Write slim orphan registry (mats_without_dump + pending lists).

    Args:    data  full catalog or orphan-only dict
    Returns: none
    """
    map_path = path or MAT_RUN_MAP_PATH
    if "mats_without_dump" in data:
        payload = {
            "_comment": data.get(
                "_comment",
                "Orphan / pending mats. Paired runs: TestMatrix_lab_runs.csv.",
            ),
            "mats_without_dump": data.get("mats_without_dump", {}),
            "pending_mat_upload": data.get("pending_mat_upload", []),
            "pending_dump_upload": data.get("pending_dump_upload", []),
            "duplicate_mats": data.get("duplicate_mats", {}),
        }
    else:
        orphans: dict = {
            "mats_without_dump": {},
            "pending_mat_upload": list(data.get("pending_mat_upload", [])),
            "pending_dump_upload": list(data.get("pending_dump_upload", [])),
            "duplicate_mats": dict(data.get("duplicate_mats", {})),
        }
        for mat_name, info in data.get("mats", {}).items():
            if not info.get("run"):
                orphans["mats_without_dump"][mat_name] = {
                    k: info[k] for k in ("row", "labTrial", "note") if k in info
                }
        payload = {
            "_comment": "Orphan / pending mats. Paired runs: TestMatrix_lab_runs.csv.",
            **orphans,
        }
    map_path.parent.mkdir(parents=True, exist_ok=True)
    map_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


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
