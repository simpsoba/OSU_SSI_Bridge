#!/usr/bin/env python3
"""
Goals
-----
Entry point for pier-top compare plots by physical-model group.

  python plot/PlotEQCompareRuns.py
  python plot/PlotEQCompareRuns.py <runDir> ...

Writes only pairwise figures (interim ref vs each other dump):

  OSU_SSI_BRIDGE_DATA_LOCAL/plots/compare/<group>/pairs/
    hist_ux_pair_<other>.png
    reference.txt

Same as `python plot/PlotEQComparePairs.py`. Shared dump I/O helpers in this
file are imported by PlotEQComparePairs.

Units (from lab_paths)
----------------------
  λ = 2.4 (cylinder length scale). Froude time scale = √λ.
  OpenSees pier recorder: numerical time t_num (s, prototype), disp (m).
  Mat Time / meaSigOS: lab real time t_lab (s, model), disp (m).

Groups come from TestMatrix_lab_runs.csv via compare_groups.py. Folder
Storm_Wave nicknames on 2026-08-21 are EQ then tsunami, not storm-wave-only.
Interim pairwise ref = fewest D5–95 slowdowns (exclude baseline / single-precision
CuDSS); offline OpenSees (no OpenFresco) is the eventual true reference.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from compare_groups import dump_to_group, group_label, groups_by_dump, legend_labels_for_dumps
from lab_paths import (
    CYLINDER_LENGTH_SCALE,
    MAT_EXTRACT_DIR,
    M_TO_MM,
    TIME_SCALE_FROUDE,
    compare_plots_dir,
    resolve_opensees_data,
    run_to_mat_from_csv,
)
from paths import HERE

# ------------------------------------------------------------
# knobs
# ------------------------------------------------------------
REPO = HERE.parent

DPI = 300
TREC_TOL_S = 1.0  # complete if t_last >= Trec - this (match PlotEQ)
UX_ABS_MAX_MM = 1000.0  # skip dumps that clearly diverged (~100 cm)

# r±NN_YYYYMMDD_HHMM_…
FOLDER_RE = re.compile(r"^r([+-]?\d+)_(\d{8})_(\d{4})_", re.IGNORECASE)

COLORS = (
    "#1565c0",
    "#c45c12",
    "#2e7d32",
    "#6a1b9a",
    "#00838f",
    "#ad1457",
    "#5d4037",
    "#455a64",
    "#0277bd",
    "#ef6c00",
    "#558b2f",
    "#7b1fa2",
)

HELP = """\
usage: python plot/PlotEQCompareRuns.py [runDir ...]

  no args   dumps from TestMatrix_lab_runs.csv (grouped by model knobs)
  runDir    one or more dump folders (still grouped if in lab_runs)
  output    LOCAL/plots/compare/<group>/pairs/hist_ux_pair_*.png
  groups    soilMesh + soilProfile + soilEleType + constitutive +
            expElement + hold + non-default xi  (not solver/np)
  ref       fewest typeConv3->2 in GM D5-95 (see PlotEQComparePairs)
"""


# ------------------------------------------------------------
# 1. FIGURE STYLE
# ------------------------------------------------------------


def apply_paper_style() -> None:
    """Serif paper look shared by all compare PNGs."""
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "cm",
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "legend.fontsize": 7,
            "legend.frameon": False,
            "axes.grid": True,
            "grid.alpha": 0.28,
            "grid.linestyle": ":",
            "grid.linewidth": 0.5,
            "lines.linewidth": 0.95,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.12,
        }
    )


def legend_ncol(n_items: int) -> int:
    """
    How many legend columns keep the row short without crowding.

    Args:    n_items  number of legend entries
    Returns: column count
    """
    if n_items <= 4:
        return n_items
    if n_items <= 10:
        return 4
    if n_items <= 18:
        return 5
    return 6


def place_legend_outside(fig: plt.Figure, ax: plt.Axes, n_items: int) -> None:
    """
    Legend below the axes (axes coords). Relies on bbox_inches=tight at save.

    Args:    fig, ax, n_items
    Returns: none
    """
    if n_items <= 0:
        return
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return
    ncol = legend_ncol(n_items)
    nrows = int(np.ceil(n_items / ncol))
    # one row ≈ 0.055 below the x-label band
    y_anchor = -0.14 - 0.055 * max(0, nrows - 1)
    ax.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, y_anchor),
        ncol=ncol,
        fontsize=6.5,
        frameon=False,
        columnspacing=1.2,
        handlelength=1.45,
        handletextpad=0.35,
        borderaxespad=0.0,
    )


def sym_ylim(ax: plt.Axes, pad: float = 1.05) -> None:
    """
    Center displacement histories on zero: ylim = ±pad · max|y|.

    Args:    ax, pad
    Returns: none (updates ax)
    """
    y_max = 0.0
    for line in ax.get_lines():
        y = np.asarray(line.get_ydata(), dtype=float)
        if y.size and np.isfinite(y).any():
            y_max = max(y_max, float(np.nanmax(np.abs(y))))
    if not np.isfinite(y_max) or y_max <= 0.0:
        y_max = 1.0
    ax.set_ylim(-pad * y_max, pad * y_max)


def dual_proto_model_axes(ax: plt.Axes) -> None:
    """
    Primary (bottom/left) = lab time on prototype scale; secondary = model.

      t_lab,proto = t_lab · √λ ,   u_m = u / λ

    Args:    ax  (data already plotted in prototype units from mat Time)
    Returns: none
    """
    length_scale = CYLINDER_LENGTH_SCALE
    time_scale = TIME_SCALE_FROUDE
    ax.set_xlabel(r"$t_\mathrm{lab}\sqrt{\lambda}$ (s)")
    ax.set_ylabel(r"$\Delta u$ (mm)")
    sec_x = ax.secondary_xaxis(
        "top",
        functions=(lambda t_proto: t_proto / time_scale, lambda t_model: t_model * time_scale),
    )
    sec_x.set_xlabel(r"$t_\mathrm{lab}$ (s)")
    sec_y = ax.secondary_yaxis(
        "right",
        functions=(
            lambda u_proto: u_proto / length_scale,
            lambda u_model: u_model * length_scale,
        ),
    )
    sec_y.set_ylabel(r"$\Delta u/\lambda$ (mm) model scale")


# ------------------------------------------------------------
# 2. DUMP I/O (OpenSees recorders + window_meta)
# ------------------------------------------------------------


def loadtxt_partial(path: Path) -> np.ndarray:
    """
    Load a numeric table; stop at the first bad/truncated line
    (common on incomplete uploads).

    Args:    path
    Returns: array shape (n, ncol), or empty (0, 0)
    """
    try:
        data = np.loadtxt(path)
    except ValueError:
        rows: list[list[float]] = []
        ncol: int | None = None
        with path.open() as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                try:
                    row = [float(x) for x in s.split()]
                except ValueError:
                    break
                if ncol is None:
                    ncol = len(row)
                if len(row) != ncol:
                    break
                rows.append(row)
        data = np.asarray(rows, dtype=float)
    if data.size == 0:
        return np.empty((0, 0))
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def find_pier_top(eq_dir: Path) -> Path | None:
    """
    Pier-top recorder path for a dump folder.

    Args:    eq_dir  OpenSees run dump
    Returns: pier_top_disp.out[.pid] or None
    """
    serial = eq_dir / "pier_top_disp.out"
    if serial.is_file():
        return serial
    shards = sorted(eq_dir.glob("pier_top_disp.out.*"))
    if not shards:
        return None
    for shard in shards:
        if shard.name.endswith(".0"):
            return shard
    return max(shards, key=lambda p: p.stat().st_size)


def read_meta(eq_dir: Path) -> dict[str, str]:
    """
    Parse window_meta.txt[.0] as key → rest-of-line.

    Args:    eq_dir
    Returns: dict (empty if missing)
    """
    meta: dict[str, str] = {}
    for name in ("window_meta.txt.0", "window_meta.txt"):
        path = eq_dir / name
        if not path.is_file():
            continue
        for line in path.read_text().splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            key, _, rest = s.partition(" ")
            meta[key] = rest.strip()
        break
    return meta


def run_duration(eq_dir: Path) -> tuple[float | None, float, bool]:
    """
    How far the pier recorder got vs meta Trec.

    Args:    eq_dir
    Returns: (t_last_s or None, Trec_s, is_complete)
             complete := Trec > 0 and t_last >= Trec - TREC_TOL_S
    """
    meta = read_meta(eq_dir)
    try:
        trec_s = float(meta.get("Trec", 0) or 0)
    except ValueError:
        trec_s = 0.0
    pier = find_pier_top(eq_dir)
    if pier is None:
        return None, trec_s, False
    data = loadtxt_partial(pier)
    if data.size == 0:
        return None, trec_s, False
    t_last_s = float(data[-1, 0])
    complete = trec_s > 0.0 and t_last_s >= trec_s - TREC_TOL_S
    return t_last_s, trec_s, complete


def eq_motion_end_s(eq_dir: Path, trec_s: float) -> float | None:
    """
    Wall-clock end of EQ motion on the OpenSees clock: Trec - freeVibT.

    Args:    eq_dir, trec_s  (s)
    Returns: t_eq_end (s) or None if freeVibT missing / invalid
    """
    meta = read_meta(eq_dir)
    try:
        free_vib_s = float(meta.get("freeVibT", 0) or 0)
    except ValueError:
        return None
    if trec_s > free_vib_s > 0.0:
        return trec_s - free_vib_s
    return None


def short_label(eq_dir: Path, incomplete: bool = False) -> str:
    """
    Legend tag from TestMatrix_lab_runs (fallback: folder name).

    Prefer legend_labels_for_dumps when plotting a group so HH:MM
    disambiguation applies. This helper is for one-off calls.

    Args:    eq_dir, incomplete  (append * when True)
    Returns: label string
    """
    labels = legend_labels_for_dumps(
        [eq_dir.name],
        incomplete={eq_dir.name: incomplete},
    )
    return labels[eq_dir.name]


# ------------------------------------------------------------
# 3. SIMULINK meaSigOS helpers
# ------------------------------------------------------------


def run_to_mat_name() -> dict[str, str]:
    """
    OpenSees dump folder name → Simulink .mat file name.

    Args:    none (reads TestMatrix_lab_runs.csv)
    Returns: {run_folder: mat_name}
    """
    return run_to_mat_from_csv()


def load_mat_mea_feedback(mat_name: str) -> tuple[np.ndarray, np.ndarray] | None:
    """
    Actuator *measurement* (feedback) from SyncLabBackup extract.

    Mats store model-scale Time (s) and meaSigOS primary (m).
    Not the same channel as comSigOS (command) — see PlotMatOS.

    Args:    mat_name  e.g. 0821_GusBridge_rowNeg16_8Core_P1.mat
    Returns: (t_model_s, u_model_m) or None if extract missing
    """
    npz_path = MAT_EXTRACT_DIR / f"{Path(mat_name).stem}.npz"
    if not npz_path.is_file():
        return None
    z = np.load(npz_path, allow_pickle=True)

    # Preferred: SyncLabBackup already pulled Time + first non-Time column
    if "meaSigOS_time" in z and "meaSigOS_primary" in z:
        t_model = np.asarray(z["meaSigOS_time"], dtype=float)
        u_model = np.asarray(z["meaSigOS_primary"], dtype=float)
        return t_model, u_model

    if "meaSigOS_data" not in z:
        return None
    data = np.asarray(z["meaSigOS_data"], dtype=float)
    names = [str(x) for x in z["meaSigOS_signalNames"].tolist()]
    time_col = next(i for i, n in enumerate(names) if n.strip().lower() == "time")
    signal_col = next(i for i in range(data.shape[1]) if i != time_col)
    return data[:, time_col], data[:, signal_col]


def model_disp_to_proto_mm(u_model_m: np.ndarray) -> np.ndarray:
    """
    Model-scale metres → prototype millimetres (relative to first sample).

      u_proto_mm = (u_model - u_model[0]) · λ · 1000

    Args:    u_model_m  (m, model)
    Returns: Δu (mm, prototype)
    """
    return (u_model_m - u_model_m[0]) * CYLINDER_LENGTH_SCALE * M_TO_MM


def model_disp_to_proto_cm(u_model_m: np.ndarray) -> np.ndarray:
    """Deprecated alias for ``model_disp_to_proto_mm`` (returns mm, not cm)."""
    return model_disp_to_proto_mm(u_model_m)


# ------------------------------------------------------------
# 4. DISCOVER / GROUP / WRITE (pairs only)
# ------------------------------------------------------------


def discover_runs() -> list[Path]:
    """
    Every dump under resolve_opensees_data() that has pier_top_disp.

    Args:    none
    Returns: sorted list of dump dirs
    """
    root = resolve_opensees_data()
    if root is None:
        return []
    return sorted(
        p for p in root.iterdir() if p.is_dir() and find_pier_top(p) is not None
    )


def resolve_run_path(dump_name: str, root: Path | None) -> Path | None:
    """
    Find a DumpFolder on Shared Drive / junction, else LOCAL mirror.

    Args:    dump_name, root  (opensees_data or None)
    Returns: path or None
    """
    if root is not None:
        path = root / dump_name
        if path.is_dir():
            return path
    local = REPO / "OSU_SSI_BRIDGE_DATA_LOCAL" / "opensees_data" / dump_name
    return local if local.is_dir() else None


def write_group_plots(
    slug: str,
    runs: list[Path],
    run_mat: dict[str, str],
    label: str,
) -> None:
    """
    Pairwise PNGs for one model-knob group (delegates to PlotEQComparePairs).

    Args:    slug, runs, run_mat, label  (human string for the log)
    Returns: none
    """
    # Lazy import: PlotEQComparePairs imports helpers from this module.
    from PlotEQComparePairs import arias_significant_duration, write_group_pairs

    out_dir = compare_plots_dir(slug) / "pairs"
    print(f"\nPlotEQCompareRuns: group {slug}")
    print(f"  {label}")
    print(f"  {len(runs)} dump(s) -> {out_dir}")
    write_group_pairs(slug, runs, run_mat, label, arias_significant_duration())


def main() -> int:
    if any(a in ("-h", "--help") for a in sys.argv[1:]):
        print(HELP, end="")
        return 0

    root = resolve_opensees_data()
    run_mat = run_to_mat_name()
    dump_group = dump_to_group()
    lab_groups = groups_by_dump()

    # --- CLI: explicit dump folders ---
    if len(sys.argv) > 1:
        buckets: dict[str, list[Path]] = {}
        labels: dict[str, str] = {}
        for arg in sys.argv[1:]:
            path = Path(arg).resolve()
            if not path.is_dir():
                print(f"PlotEQCompareRuns: skip (not a dir) {arg}", file=sys.stderr)
                continue
            slug = dump_group.get(path.name, "cli")
            buckets.setdefault(slug, []).append(path)
            labels.setdefault(slug, slug)
        if not buckets:
            print("PlotEQCompareRuns: no valid run dirs", file=sys.stderr)
            return 1
        for slug, runs in buckets.items():
            write_group_plots(slug, runs, run_mat, labels[slug])
        return 0

    # --- default: curated lab matrix ---
    if not lab_groups:
        print(
            "PlotEQCompareRuns: empty TestMatrix_lab_runs.csv; "
            "falling back to flat discover",
            file=sys.stderr,
        )
        runs = discover_runs()
        if not runs:
            return 1
        write_group_plots("all", runs, run_mat, "all discovered dumps")
        return 0

    n_groups = 0
    for slug, rows in lab_groups.items():
        runs: list[Path] = []
        for row in rows:
            path = resolve_run_path(row["DumpFolder"], root)
            if path is None:
                print(f"PlotEQCompareRuns: missing dump {row['DumpFolder']}")
                continue
            if find_pier_top(path) is None:
                print(f"PlotEQCompareRuns: no pier_top {path.name}")
                continue
            runs.append(path)
        if not runs:
            continue
        write_group_plots(slug, runs, run_mat, group_label(rows[0]))
        n_groups += 1

    print(f"\nPlotEQCompareRuns: {n_groups} group folder(s) under compare/*/pairs/")
    return 0 if n_groups else 1


if __name__ == "__main__":
    raise SystemExit(main())
