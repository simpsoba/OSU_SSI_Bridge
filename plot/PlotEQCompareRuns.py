#!/usr/bin/env python3
"""
Goals
-----
Overlay pier-top Δux for lab dumps that share the same *physical* model
(mesh, soil, element type, constitutive, hybrid hold, non-default ξ).
Solver / np / integrator variants stay in one compare folder so overlays
are fair.

  python plot/PlotEQCompareRuns.py
  python plot/PlotEQCompareRuns.py <runDir> ...

Writes under OSU_SSI_BRIDGE_DATA_LOCAL/plots/compare/<group>/ :

  hist_ux_complete.png          OpenSees pier_top (numerical t), complete only
  hist_ux_all.png               same, complete + incomplete (*)
  hist_ux_*_realtime.png        Simulink meaSigOS feedback (no OpenSees interp)

Units (from lab_paths)
----------------------
  λ = 2.4 (cylinder length scale). Froude time scale = √λ.
  OpenSees pier recorder: prototype time (s), prototype disp (m).
  Mat Time / meaSigOS: model real time (s), model disp (m).
  Realtime plots: primary axes = prototype; secondary = model (_m).

Groups come from TestMatrix_lab_runs.csv via compare_groups.py — not from
folder nicknames alone.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from compare_groups import dump_to_group, group_label, groups_by_dump
from lab_paths import (
    CYLINDER_LENGTH_SCALE,
    MAT_EXTRACT_DIR,
    TIME_SCALE_FROUDE,
    load_mat_run_map,
    plots_root,
    resolve_opensees_data,
)
from paths import HERE

# ------------------------------------------------------------
# knobs
# ------------------------------------------------------------
REPO = HERE.parent
PLOTS_ROOT = plots_root()

DPI = 300
TREC_TOL_S = 1.0  # complete if t_last >= Trec - this (match PlotEQ)
MAX_PLOT_PTS = 40_000  # downsample long DAQ / recorder series for PNG size
UX_ABS_MAX_CM = 100.0  # skip dumps that clearly diverged

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
  output    OSU_SSI_BRIDGE_DATA_LOCAL/plots/compare/<group>/hist_ux_*.png
  groups    soilMesh + soilProfile + soilEleType + constitutive +
            expElement + hold + non-default xi  (not solver/np)
  complete  t_last >= Trec - 1 s
  realtime  meaSigOS (actuator feedback); dual proto/model axes
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
    Primary (bottom/left) = prototype; secondary (top/right) = model.

      t_m = t / √λ ,   u_m = u / λ

    Args:    ax  (data already plotted in prototype units)
    Returns: none
    """
    length_scale = CYLINDER_LENGTH_SCALE
    time_scale = TIME_SCALE_FROUDE
    ax.set_xlabel(r"$t$ (s)")
    ax.set_ylabel(r"$\Delta u$ (cm)")
    sec_x = ax.secondary_xaxis(
        "top",
        functions=(lambda t_proto: t_proto / time_scale, lambda t_model: t_model * time_scale),
    )
    sec_x.set_xlabel(r"$t_\mathrm{m}$ (s)")
    sec_y = ax.secondary_yaxis(
        "right",
        functions=(
            lambda u_proto: u_proto / length_scale,
            lambda u_model: u_model * length_scale,
        ),
    )
    sec_y.set_ylabel(r"$\Delta u_\mathrm{m}$ (cm)")


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
    Compact legend tag: r±NN HH:MM (np=N)[*].

    Args:    eq_dir, incomplete  (append * when True)
    Returns: label string
    """
    meta = read_meta(eq_dir)
    np_str = meta.get("np", "?")
    match = FOLDER_RE.match(eq_dir.name)
    if match:
        hhmm = match.group(3)
        tag = f"r{match.group(1)} {hhmm[:2]}:{hhmm[2:]} (np={np_str})"
    else:
        tag = f"{eq_dir.name} (np={np_str})"
    if incomplete:
        tag += "*"
    return tag


def decimate(
    t: np.ndarray,
    y: np.ndarray,
    n_max: int = MAX_PLOT_PTS,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Thin a time series for plotting (keep endpoints of every stride).

    Args:    t, y  same length; n_max  max points kept
    Returns: (t_plot, y_plot)
    """
    n = len(t)
    if n <= n_max:
        return t, y
    step = int(np.ceil(n / n_max))
    return t[::step], y[::step]


# ------------------------------------------------------------
# 3. SIMULINK meaSigOS (realtime overlays)
# ------------------------------------------------------------


def run_to_mat_name() -> dict[str, str]:
    """
    OpenSees dump folder name → Simulink .mat file name.

    Args:    none (reads mat_run_map.json)
    Returns: {run_folder: mat_name}
    """
    mat_map = load_mat_run_map()
    out: dict[str, str] = {}
    for mat_name, info in mat_map.get("mats", {}).items():
        run = info.get("run")
        if run:
            out[run] = mat_name
    return out


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


def model_disp_to_proto_cm(u_model_m: np.ndarray) -> np.ndarray:
    """
    Model-scale metres → prototype centimetres (relative to first sample).

      u_proto_cm = (u_model - u_model[0]) · λ · 100

    Args:    u_model_m  (m, model)
    Returns: Δu (cm, prototype)
    """
    return (u_model_m - u_model_m[0]) * CYLINDER_LENGTH_SCALE * 100.0


# ------------------------------------------------------------
# 4. PLOT 1 — OpenSees pier_top (numerical prototype clock)
# ------------------------------------------------------------


def plot_compare_opensees(
    runs: list[Path],
    out: Path,
    mark_incomplete: bool,
) -> None:
    """
    Overlay OpenSees pier_top_disp Δux (cm) vs numerical t (s).

    Complete runs: solid. Incomplete: dashed + * in legend; dotted
    vertical at t_last. Grey dashed = EQ motion end (median of group).

    Args:    runs, out, mark_incomplete
    Returns: none (writes PNG)
    """
    out.parent.mkdir(parents=True, exist_ok=True)
    apply_paper_style()
    fig, ax = plt.subplots(figsize=(7.2, 3.6))

    t_eq_ends: list[float] = []
    n_plotted = 0

    for eq_dir in runs:
        pier = find_pier_top(eq_dir)
        if pier is None:
            print(f"PlotEQCompareRuns: skip {eq_dir.name} (no pier_top_disp)")
            continue
        data = loadtxt_partial(pier)
        if data.size == 0 or data.shape[1] < 2:
            print(f"PlotEQCompareRuns: skip {eq_dir.name} (empty pier_top)")
            continue

        t_last_s, trec_s, complete = run_duration(eq_dir)
        t_s = data[:, 0]
        ux_m = data[:, 1] - data[0, 1]
        ux_cm = ux_m * 100.0

        if not np.all(np.isfinite(ux_cm)) or float(np.nanmax(np.abs(ux_cm))) > UX_ABS_MAX_CM:
            print(f"PlotEQCompareRuns: skip {eq_dir.name} (non-physical pier ux)")
            continue

        color = COLORS[n_plotted % len(COLORS)]
        line_style = "-" if complete else "--"
        ax.plot(
            t_s,
            ux_cm,
            color=color,
            lw=1.0 if complete else 0.9,
            ls=line_style,
            label=short_label(eq_dir, incomplete=mark_incomplete and not complete),
        )
        n_plotted += 1

        t_eq = eq_motion_end_s(eq_dir, trec_s)
        if t_eq is not None:
            t_eq_ends.append(t_eq)
        if not complete and t_last_s is not None:
            ax.axvline(float(t_last_s), color=color, lw=0.6, ls=":", alpha=0.5)

    if t_eq_ends:
        ax.axvline(
            float(np.median(t_eq_ends)),
            color="0.45",
            lw=0.8,
            ls="--",
            label="EQ end",
        )

    ax.set_xlabel(r"$t$ (s)")
    ax.set_ylabel(r"$\Delta u_x$ (cm)")
    sym_ylim(ax)
    n_legend = n_plotted + (1 if t_eq_ends else 0)
    place_legend_outside(fig, ax, n_legend)
    fig.savefig(out, dpi=DPI, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    print(f"PlotEQCompareRuns: wrote {out}  ({n_plotted} series)")


# ------------------------------------------------------------
# 5. PLOT 2 — meaSigOS realtime (model mat → prototype axes)
# ------------------------------------------------------------


def plot_compare_realtime(
    runs: list[Path],
    out: Path,
    mark_incomplete: bool,
    run_mat: dict[str, str],
) -> None:
    """
    Overlay Simulink meaSigOS feedback on prototype axes.

    Skips dumps with no mat in mat_run_map (see runs_without_mat).
    Does *not* interpolate onto OpenSees time — DAQ clock only.

    Args:    runs, out, mark_incomplete, run_mat  {dump_name: mat_file}
    Returns: none (writes PNG)
    """
    out.parent.mkdir(parents=True, exist_ok=True)
    apply_paper_style()
    fig, ax = plt.subplots(figsize=(7.2, 3.8))

    t_eq_ends: list[float] = []
    n_plotted = 0

    for eq_dir in runs:
        mat_name = run_mat.get(eq_dir.name)
        if not mat_name:
            print(f"PlotEQCompareRuns: skip realtime {eq_dir.name} (no mat)")
            continue
        pair = load_mat_mea_feedback(mat_name)
        if pair is None:
            print(f"PlotEQCompareRuns: skip realtime {eq_dir.name} (no meaSigOS)")
            continue
        t_model_s, u_model_m = pair
        if t_model_s.size == 0:
            continue

        _, trec_s, complete = run_duration(eq_dir)

        # Scale model mat → prototype for the primary axes
        t_proto_s = t_model_s * TIME_SCALE_FROUDE
        u_proto_cm = model_disp_to_proto_cm(u_model_m)
        t_plot, u_plot = decimate(t_proto_s, u_proto_cm)

        color = COLORS[n_plotted % len(COLORS)]
        line_style = "-" if complete else "--"
        ax.plot(
            t_plot,
            u_plot,
            color=color,
            lw=1.0 if complete else 0.9,
            ls=line_style,
            label=short_label(eq_dir, incomplete=mark_incomplete and not complete),
        )
        n_plotted += 1

        t_eq = eq_motion_end_s(eq_dir, trec_s)
        if t_eq is not None:
            t_eq_ends.append(t_eq)
        if not complete:
            # last lab sample, expressed on the prototype time axis
            ax.axvline(
                float(t_model_s[-1] * TIME_SCALE_FROUDE),
                color=color,
                lw=0.6,
                ls=":",
                alpha=0.5,
            )

    if n_plotted == 0:
        print(f"PlotEQCompareRuns: no series for {out.name}")
        plt.close(fig)
        return

    if t_eq_ends:
        ax.axvline(
            float(np.median(t_eq_ends)),
            color="0.45",
            lw=0.8,
            ls="--",
            label="EQ end",
        )

    dual_proto_model_axes(ax)
    sym_ylim(ax)
    n_legend = n_plotted + (1 if t_eq_ends else 0)
    place_legend_outside(fig, ax, n_legend)
    fig.savefig(out, dpi=DPI, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    print(f"PlotEQCompareRuns: wrote {out}  ({n_plotted} series)")


# ------------------------------------------------------------
# 6. DISCOVER / GROUP / WRITE
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
    Four PNGs for one model-knob group (complete + all × OS + realtime).

    Args:    slug, runs, run_mat, label  (human string for the log)
    Returns: none
    """
    out_dir = PLOTS_ROOT / "compare" / slug
    complete_runs = [r for r in runs if run_duration(r)[2]]
    print(f"\nPlotEQCompareRuns: group {slug}")
    print(f"  {label}")
    print(f"  {len(runs)} dump(s), {len(complete_runs)} complete -> {out_dir}")

    if complete_runs:
        plot_compare_opensees(
            complete_runs,
            out_dir / "hist_ux_complete.png",
            mark_incomplete=False,
        )
        plot_compare_realtime(
            complete_runs,
            out_dir / "hist_ux_complete_realtime.png",
            mark_incomplete=False,
            run_mat=run_mat,
        )

    plot_compare_opensees(
        runs,
        out_dir / "hist_ux_all.png",
        mark_incomplete=True,
    )
    plot_compare_realtime(
        runs,
        out_dir / "hist_ux_all_realtime.png",
        mark_incomplete=True,
        run_mat=run_mat,
    )


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

    print(f"\nPlotEQCompareRuns: {n_groups} group folder(s) under compare/")
    return 0 if n_groups else 1


if __name__ == "__main__":
    raise SystemExit(main())
