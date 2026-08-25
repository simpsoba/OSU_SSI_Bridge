#!/usr/bin/env python3
"""
Goals
-----
Plot Simulink *OS channels on model and prototype clocks, then compare mats
that share the same physical-model group from compare_groups.py.

  python plot/PlotMatOS.py
  python plot/PlotMatOS.py 0821_GusBridge_rowPos1.mat ...

Reads SyncLabBackup extracts from
OSU_SSI_BRIDGE_DATA_LOCAL/opensees_data/mat_extract/ and maps mats to run
folders with plot/opensees_data/mat_run_map.json.

Units
-----
  Mat Time: model-scale lab real time (s).
  Prototype time: Mat Time × √λ, where λ = CYLINDER_LENGTH_SCALE (2.4).
  Mat actuator / *OS displacement: model scale (m).
  Plotted displacement: prototype scale (cm), converted by λ × 100.
  Per-run plots use model time below and prototype time above.
  Compare plots use prototype time.

Writes
------
Per mapped run, under OSU_SSI_BRIDGE_DATA_LOCAL/plots/<run>/:

  hist_os_com.png
  hist_os_state.png
  hist_os_tar_com_mea.png

Unmapped mats go under plots/mat_os/<stem>/. Compare groups go under
plots/compare/<group>/ and write hist_os_com.png and hist_os_state.png.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

from compare_groups import dump_to_group, groups_by_dump
from lab_paths import (
    CYLINDER_LENGTH_SCALE,
    MAT_EXTRACT_DIR,
    TIME_SCALE_FROUDE,
    load_mat_run_map,
    plots_root,
)
from paths import HERE


# ------------------------------------------------------------
# 1. PATHS, UNITS, AND PLOT KNOBS
# ------------------------------------------------------------

REPO = HERE.parent
PLOTS_ROOT = plots_root()

DPI = 300
MAX_PLOT_PTS = 40_000

# Model-scale displacement (m) to prototype displacement (cm).
DISP_M_TO_PROTO_CM = CYLINDER_LENGTH_SCALE * 100.0

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

# OpenFresco typeConv* names used by this campaign.
STATE_LABELS = {
    "typeConv3": "simState",
    "typeConv1": "substep/count",
    "s1": "flag s1",
    "s2": "flag s2",
    "s3": "flag s3",
}

HELP = """\
usage: python plot/PlotMatOS.py [matFile ...]

  no args   every .mat in mat_run_map.json (plus any *.npz extract)
  matFile   one or more .mat names or stems
  needs     SyncLabBackup extracts in OSU_SSI_BRIDGE_DATA_LOCAL/.../mat_extract/
  per-run   hist_os_com.png, hist_os_state.png, hist_os_tar_com_mea.png
  compare   OSU_SSI_BRIDGE_DATA_LOCAL/plots/compare/<group>/hist_os_*.png
"""


# ------------------------------------------------------------
# 2. LABELS AND FIGURE STYLE
# ------------------------------------------------------------


def short_sig(name: str) -> str:
    """
    Keep only the final component of a Simulink signal path.

    Args:    name  signal path
    Returns: short signal name
    """
    signal_name = name.replace("\\", "/").strip()
    if "/" in signal_name:
        signal_name = signal_name.rsplit("/", 1)[-1]
    return signal_name


def state_label(name: str) -> str:
    """
    Give a stateOS signal its campaign-friendly label.

    Args:    name  stateOS signal path
    Returns: display label
    """
    signal_name = short_sig(name)
    return STATE_LABELS.get(signal_name, signal_name)


def apply_paper_style() -> None:
    """
    Apply the serif paper style used by compare PNGs.

    Args:    none
    Returns: none (updates matplotlib rcParams)
    """
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


def legend_ncol(n: int) -> int:
    """
    Choose a compact legend column count.

    Args:    n  number of legend entries
    Returns: legend column count
    """
    if n <= 4:
        return n
    if n <= 10:
        return 4
    if n <= 18:
        return 5
    return 6


def place_legend_outside(ax: plt.Axes, n: int) -> None:
    """
    Place a multirow legend below an axes.

    Args:    ax, n  axes and number of legend entries
    Returns: none (updates ax)
    """
    if n <= 0:
        return
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return

    n_columns = legend_ncol(n)
    n_rows = int(np.ceil(n / n_columns))
    y_anchor = -0.14 - 0.055 * max(0, n_rows - 1)
    ax.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, y_anchor),
        ncol=n_columns,
        fontsize=6.5,
        frameon=False,
        columnspacing=1.2,
        handlelength=1.45,
        handletextpad=0.35,
        borderaxespad=0.0,
    )


def dual_time_xaxis(ax) -> None:
    """
    Add model time below and Froude-scaled prototype time above.

    Args:    ax  axes plotted against model time (s)
    Returns: none (updates ax)
    """
    ax.set_xlabel("t (s) model scale (lab real)")
    prototype_axis = ax.secondary_xaxis(
        "top",
        functions=(
            lambda time_model: time_model * TIME_SCALE_FROUDE,
            lambda time_proto: time_proto / TIME_SCALE_FROUDE,
        ),
    )
    prototype_axis.set_xlabel(
        r"t (s) prototype scale (lab real $\times\sqrt{2.4}$)"
    )


def integer_yaxis(ax) -> None:
    """
    Restrict stateOS major ticks to integer values.

    Args:    ax  stateOS axes
    Returns: none (updates ax)
    """
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, min_n_ticks=2))


# ------------------------------------------------------------
# 3. EXTRACT I/O AND SERIES DISCOVERY
# ------------------------------------------------------------


def decimate(
    t: np.ndarray,
    y: np.ndarray,
    nmax: int = MAX_PLOT_PTS,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Thin a time series to at most nmax plotting points.

    Args:    t  time (s); y  signal values; nmax  maximum points
    Returns: (t_plot, y_plot)
    """
    n_points = len(t)
    if n_points <= nmax:
        return t, y
    step = int(np.ceil(n_points / nmax))
    return t[::step], y[::step]


def load_block(
    z: np.lib.npyio.NpzFile,
    key: str,
) -> tuple[np.ndarray, list[str], np.ndarray] | None:
    """
    Read one tarSigOS, comSigOS, meaSigOS, or stateOS extract block.

    Args:    z  npz extract; key  block prefix
    Returns: (model_time_s, signal_names, data) or None
    """
    time_key = f"{key}_time"
    data_key = f"{key}_data"
    names_key = f"{key}_signalNames"
    if data_key not in z:
        return None

    data = np.asarray(z[data_key], dtype=float)
    names = [str(value) for value in z[names_key].tolist()] if names_key in z else []
    if time_key in z:
        time_model_s = np.asarray(z[time_key], dtype=float)
    else:
        time_column = next(
            index
            for index, name in enumerate(names)
            if name.strip().lower() == "time"
        )
        time_model_s = data[:, time_column]
    return time_model_s, names, data


def non_time_cols(names: list[str], ncol: int) -> list[int]:
    """
    List data columns other than the named Time column.

    Args:    names  column names; ncol  number of data columns
    Returns: non-Time column indices
    """
    time_column = next(
        (
            index
            for index, name in enumerate(names)
            if name.strip().lower() == "time"
        ),
        None,
    )
    return [index for index in range(ncol) if index != time_column]


def out_dir_for_mat(mat_name: str, mmap: dict) -> Path:
    """
    Choose the mapped run plot folder or an unmapped-mat folder.

    Args:    mat_name  mat file name; mmap  mat-to-run mapping
    Returns: output directory
    """
    mat_info = mmap.get("mats", {}).get(mat_name) or mmap.get("mats", {}).get(
        Path(mat_name).name
    )
    run_name = mat_info.get("run") if mat_info else None
    if run_name:
        return PLOTS_ROOT / run_name
    return PLOTS_ROOT / "mat_os" / Path(mat_name).stem


def iter_mat_npz(
    mmap: dict,
    mats: list[str] | None = None,
) -> list[tuple[str, Path]]:
    """
    Find available npz extracts for requested or mapped mats.

    Args:    mmap  mat-to-run mapping; mats  requested mat names or None
    Returns: list of (mat_file_name, npz_path)
    """
    mat_names = mats or sorted(mmap.get("mats", {}).keys())
    if mats is None and MAT_EXTRACT_DIR.is_dir():
        for npz_path in sorted(MAT_EXTRACT_DIR.glob("*.npz")):
            mat_name = f"{npz_path.stem}.mat"
            if mat_name not in mat_names:
                mat_names.append(mat_name)

    pairs: list[tuple[str, Path]] = []
    for mat_name in mat_names:
        npz_path = MAT_EXTRACT_DIR / f"{Path(mat_name).stem}.npz"
        if npz_path.is_file():
            pairs.append((Path(mat_name).name, npz_path))
    return pairs


# ------------------------------------------------------------
# 4. RUN AND COMPARE LABELS
# ------------------------------------------------------------


def compare_label(mat_name: str, mmap: dict) -> str:
    """
    Build a compact r±NN HH:MM (np=N) compare legend label.

    Args:    mat_name  mat file name; mmap  mat-to-run mapping
    Returns: compare legend label
    """
    mat_info = mmap.get("mats", {}).get(Path(mat_name).name) or {}
    run_name = mat_info.get("run")
    if not run_name:
        return Path(mat_name).stem

    folder_match = FOLDER_RE.match(run_name)
    note = mat_info.get("note") or ""
    np_value = ""
    if "np=" in note:
        np_value = note.split("np=", 1)[1].split()[0].rstrip(")")
    if folder_match:
        hhmm = folder_match.group(3)
        label = f"r{folder_match.group(1)} {hhmm[:2]}:{hhmm[2:]}"
        if np_value:
            label += f" (np={np_value})"
        return label
    return run_name


def series_label(mat_name: str, mmap: dict) -> str:
    """
    Build the full run label used in titles and state panels.

    Args:    mat_name  mat file name; mmap  mat-to-run mapping
    Returns: run label or unmapped mat stem
    """
    mat_info = mmap.get("mats", {}).get(Path(mat_name).name) or {}
    run_name = mat_info.get("run")
    note = mat_info.get("note") or ""
    if run_name:
        label = run_name
        if "np=" in note:
            label += f"  ({note})"
        return label
    return Path(mat_name).stem


# ------------------------------------------------------------
# 5. PER-MAT PNGS
# ------------------------------------------------------------


def plot_com(z: np.lib.npyio.NpzFile, out: Path, title: str) -> None:
    """
    Plot the first comSigOS command displacement with dual time axes.

    Args:    z  npz extract; out  PNG path; title  figure title
    Returns: none (writes PNG when comSigOS exists)
    """
    block = load_block(z, "comSigOS")
    if block is None:
        print(f"PlotMatOS: skip comSigOS ({title})")
        return

    time_model_s, names, data = block
    signal_columns = non_time_cols(names, data.shape[1])
    if not signal_columns:
        return

    signal_column = signal_columns[0]
    time_plot, command_plot = decimate(
        time_model_s,
        data[:, signal_column],
    )
    fig, ax = plt.subplots(figsize=(11.0, 4.4), constrained_layout=True)
    ax.plot(
        time_plot,
        command_plot * DISP_M_TO_PROTO_CM,
        color=COLORS[1],
        lw=1.05,
        label=rf"command ($\times\lambda={CYLINDER_LENGTH_SCALE:g}$)",
    )
    dual_time_xaxis(ax)
    ax.set_ylabel(r"command $u$ (cm, prototype)")
    ax.set_title(title)
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, ls=":", alpha=0.45)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(f"PlotMatOS: wrote {out}")


def plot_tar_com_mea(z: np.lib.npyio.NpzFile, out: Path, title: str) -> None:
    """
    Plot first target, command, and measurement *OS displacements together.

    Args:    z  npz extract; out  PNG path; title  figure title
    Returns: none (writes PNG when at least one *OS block exists)
    """
    fig, ax = plt.subplots(figsize=(11.0, 4.6), constrained_layout=True)
    n_plotted = 0

    for block_key, color in zip(("tarSigOS", "comSigOS", "meaSigOS"), COLORS):
        block = load_block(z, block_key)
        if block is None:
            continue
        time_model_s, names, data = block
        signal_columns = non_time_cols(names, data.shape[1])
        if not signal_columns:
            continue

        signal_column = signal_columns[0]
        signal_name = (
            short_sig(names[signal_column])
            if signal_column < len(names)
            else block_key
        )
        time_plot, displacement_plot = decimate(
            time_model_s,
            data[:, signal_column],
        )
        ax.plot(
            time_plot,
            displacement_plot * DISP_M_TO_PROTO_CM,
            color=color,
            lw=1.0,
            label=f"{block_key}: {signal_name}",
        )
        n_plotted += 1

    if n_plotted == 0:
        plt.close(fig)
        print(f"PlotMatOS: skip tar/com/mea ({title})")
        return

    dual_time_xaxis(ax)
    ax.set_ylabel(
        rf"displacement (cm, prototype, $\times\lambda={CYLINDER_LENGTH_SCALE:g}$)"
    )
    ax.set_title(title)
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, ls=":", alpha=0.45)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(f"PlotMatOS: wrote {out}")


def plot_state(z: np.lib.npyio.NpzFile, out: Path, title: str) -> None:
    """
    Plot stateOS simState above its remaining integer flags and counts.

    Args:    z  npz extract; out  PNG path; title  figure title
    Returns: none (writes PNG when stateOS exists)
    """
    block = load_block(z, "stateOS")
    if block is None:
        print(f"PlotMatOS: skip stateOS ({title})")
        return

    time_model_s, names, data = block
    state_columns = non_time_cols(names, data.shape[1])
    if not state_columns:
        return

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(11.0, 6.2),
        sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [2.0, 1.2]},
    )
    state_ax, flags_ax = axes

    sim_state_column = state_columns[0]
    time_plot, state_plot = decimate(
        time_model_s,
        data[:, sim_state_column],
    )
    state_ax.plot(
        time_plot,
        state_plot,
        color=COLORS[0],
        lw=1.05,
        drawstyle="steps-post",
        label=state_label(names[sim_state_column]),
    )
    state_ax.set_ylabel("simState")
    state_ax.set_title(title)
    state_ax.legend(fontsize=8, loc="best")
    state_ax.grid(True, ls=":", alpha=0.45)
    integer_yaxis(state_ax)

    for color_index, state_column in enumerate(state_columns[1:]):
        time_plot, state_plot = decimate(
            time_model_s,
            data[:, state_column],
        )
        flags_ax.plot(
            time_plot,
            state_plot,
            color=COLORS[(color_index + 1) % len(COLORS)],
            lw=0.95,
            drawstyle="steps-post",
            label=state_label(names[state_column]),
        )
    flags_ax.set_ylabel("flags / count")
    flags_ax.legend(fontsize=7.5, loc="best", ncol=2)
    flags_ax.grid(True, ls=":", alpha=0.45)
    integer_yaxis(flags_ax)
    dual_time_xaxis(flags_ax)

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(f"PlotMatOS: wrote {out}")


def process_mat(mat_name: str, mmap: dict) -> None:
    """
    Write all three per-mat PNGs for one extracted mat.

    Args:    mat_name  mat file name; mmap  mat-to-run mapping
    Returns: none
    """
    mat_stem = Path(mat_name).stem
    npz_path = MAT_EXTRACT_DIR / f"{mat_stem}.npz"
    if not npz_path.is_file():
        print(f"PlotMatOS: missing extract {npz_path}", file=sys.stderr)
        return

    extract = np.load(npz_path, allow_pickle=True)
    output_dir = out_dir_for_mat(mat_name, mmap)
    label = series_label(mat_name, mmap)
    plot_com(
        extract,
        output_dir / "hist_os_com.png",
        f"{label} — comSigOS (command $u$, prototype)",
    )
    plot_state(
        extract,
        output_dir / "hist_os_state.png",
        f"{label} — stateOS (OpenFresco; slowdowns)",
    )
    plot_tar_com_mea(
        extract,
        output_dir / "hist_os_tar_com_mea.png",
        f"{label} — tar / com / mea (*OS, prototype)",
    )


# ------------------------------------------------------------
# 6. COMPARE PNGS
# ------------------------------------------------------------


def plot_compare_com(
    pairs: list[tuple[str, Path]],
    mmap: dict,
    out: Path,
) -> None:
    """
    Overlay first comSigOS command histories on prototype time.

    Args:    pairs  (mat name, npz path); mmap  mapping; out  PNG path
    Returns: none (writes PNG when at least one command exists)
    """
    apply_paper_style()
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    n_plotted = 0

    for mat_name, npz_path in pairs:
        extract = np.load(npz_path, allow_pickle=True)
        block = load_block(extract, "comSigOS")
        if block is None:
            continue
        time_model_s, names, data = block
        signal_columns = non_time_cols(names, data.shape[1])
        if not signal_columns:
            continue

        time_proto_s = time_model_s * TIME_SCALE_FROUDE
        time_plot, command_plot = decimate(
            time_proto_s,
            data[:, signal_columns[0]],
        )
        ax.plot(
            time_plot,
            command_plot * DISP_M_TO_PROTO_CM,
            color=COLORS[n_plotted % len(COLORS)],
            lw=1.0,
            label=compare_label(mat_name, mmap),
        )
        n_plotted += 1

    if n_plotted == 0:
        plt.close(fig)
        print("PlotMatOS: no series for compare hist_os_com.png")
        return

    ax.set_xlabel(r"$t$ (s)")
    ax.set_ylabel(r"command $u$ (cm, prototype)")
    place_legend_outside(ax, n_plotted)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=DPI, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    print(f"PlotMatOS: wrote {out}  ({n_plotted} series)")


def plot_compare_state(
    pairs: list[tuple[str, Path]],
    mmap: dict,
    out: Path,
) -> None:
    """
    Stack one prototype-time stateOS simState panel per mat.

    Args:    pairs  (mat name, npz path); mmap  mapping; out  PNG path
    Returns: none (writes PNG when at least one state series exists)
    """
    series: list[tuple[str, np.ndarray, np.ndarray]] = []
    for mat_name, npz_path in pairs:
        extract = np.load(npz_path, allow_pickle=True)
        block = load_block(extract, "stateOS")
        if block is None:
            continue
        time_model_s, names, data = block
        state_columns = non_time_cols(names, data.shape[1])
        if not state_columns:
            continue

        time_proto_s = time_model_s * TIME_SCALE_FROUDE
        time_plot, state_plot = decimate(
            time_proto_s,
            data[:, state_columns[0]],
        )
        series.append(
            (
                series_label(mat_name, mmap),
                time_plot,
                state_plot,
            )
        )

    if not series:
        print("PlotMatOS: no series for compare hist_os_state.png")
        return

    n_series = len(series)
    figure_height = max(6.0, 1.05 * n_series)
    fig, axes = plt.subplots(
        n_series,
        1,
        figsize=(11.0, figure_height),
        sharex=True,
        constrained_layout=True,
    )
    if n_series == 1:
        axes = [axes]

    for series_index, (label, time_plot, state_plot) in enumerate(series):
        ax = axes[series_index]
        color = COLORS[series_index % len(COLORS)]
        ax.plot(
            time_plot,
            state_plot,
            color=color,
            lw=1.0,
            drawstyle="steps-post",
        )
        ax.set_ylabel("simState", fontsize=7)
        ax.grid(True, ls=":", alpha=0.45)
        integer_yaxis(ax)
        # A full legend on every panel would cover the state history.
        ax.text(
            0.01,
            0.92,
            label,
            transform=ax.transAxes,
            fontsize=6.5,
            va="top",
            ha="left",
            color=color,
            bbox={
                "facecolor": "white",
                "alpha": 0.7,
                "edgecolor": "none",
                "pad": 1.5,
            },
        )

    axes[0].set_title("stateOS simState — OpenFresco (slowdowns)")
    axes[-1].set_xlabel(
        r"t (s) lab real $\times\sqrt{2.4}$ (prototype scale)"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(f"PlotMatOS: wrote {out}  ({n_series} panels)")


# ------------------------------------------------------------
# 7. GROUPING AND CLI
# ------------------------------------------------------------


def _mat_group_slug(
    mat_name: str,
    mmap: dict,
    dump_group: dict[str, str],
) -> str:
    """
    Map a mat to its physical-model compare-group slug.

    Args:    mat_name  mat file; mmap  mapping; dump_group  run-to-group map
    Returns: compare-group slug or "unmapped"
    """
    mat_info = mmap.get("mats", {}).get(Path(mat_name).name) or {}
    run_name = mat_info.get("run")
    if run_name and run_name in dump_group:
        return dump_group[run_name]
    return "unmapped"


def write_grouped_compares(
    pairs: list[tuple[str, Path]],
    mmap: dict,
) -> None:
    """
    Write command and state compare PNGs for every represented group.

    Args:    pairs  (mat name, npz path); mmap  mat-to-run mapping
    Returns: none
    """
    dump_group = dump_to_group()
    group_pairs: dict[str, list[tuple[str, Path]]] = {}
    for mat_name, npz_path in pairs:
        group_slug = _mat_group_slug(mat_name, mmap, dump_group)
        group_pairs.setdefault(group_slug, []).append((mat_name, npz_path))

    for group_slug, subgroup_pairs in group_pairs.items():
        output_dir = PLOTS_ROOT / "compare" / group_slug
        print(
            f"PlotMatOS: compare group {group_slug}  "
            f"({len(subgroup_pairs)} mat(s))"
        )
        plot_compare_com(
            subgroup_pairs,
            mmap,
            output_dir / "hist_os_com.png",
        )
        plot_compare_state(
            subgroup_pairs,
            mmap,
            output_dir / "hist_os_state.png",
        )


def main() -> int:
    """
    Parse mat names, write per-mat PNGs, then refresh grouped compares.

    Args:    sys.argv  optional mat file names or stems
    Returns: process exit code
    """
    if any(argument in ("-h", "--help") for argument in sys.argv[1:]):
        print(HELP, end="")
        return 0

    mat_map = load_mat_run_map()
    if len(sys.argv) > 1:
        mat_names = [
            (
                Path(argument).name
                if argument.lower().endswith(".mat")
                else f"{Path(argument).stem}.mat"
            )
            for argument in sys.argv[1:]
        ]
    else:
        mat_names = None

    pairs = iter_mat_npz(mat_map, mat_names)
    if not pairs:
        print("PlotMatOS: no mats to plot", file=sys.stderr)
        return 1

    print(
        f"PlotMatOS: {len(pairs)} mat(s)  "
        f"TIME_SCALE_FROUDE={TIME_SCALE_FROUDE:.4f}"
    )
    for mat_name, _ in pairs:
        process_mat(mat_name, mat_map)
    write_grouped_compares(pairs, mat_map)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
