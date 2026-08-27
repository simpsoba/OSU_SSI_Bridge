#!/usr/bin/env python3
"""
Goals
-----
Plot Simulink *OS channels on model and prototype clocks (per mat / dump).

  python plot/PlotMatOS.py
  python plot/PlotMatOS.py 0821_GusBridge_rowPos1.mat ...

Reads SyncLabBackup extracts from
OSU_SSI_BRIDGE_DATA_LOCAL/opensees_data/mat_extract/ and maps mats to run
folders with TestMatrix_lab_runs.csv (+ orphan mats in mat_run_map.json).

Group-level pier compare lives in PlotEQComparePairs (compare/<group>/pairs/).
This script does not write overlays under compare/.

Units
-----
  Mat Time: model-scale lab real time (s).
  Prototype time: Mat Time × √λ, where λ = CYLINDER_LENGTH_SCALE (2.4).
  Mat actuator / *OS displacement: model scale (m).
  Plotted displacement: prototype scale (mm), converted by λ × 1000.
  Per-run plots use model time below and prototype time above.

Writes
------
Mapped mats → OSU_SSI_BRIDGE_DATA_LOCAL/plots/runs/<dump>/os/ :

  hist_os_com.png
  hist_os_state.png
  hist_os_tar_com_mea.png

Unmapped mats → plots/mats/<stem>/os/.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator
import numpy as np

from lab_paths import (
    CYLINDER_LENGTH_SCALE,
    DISP_M_TO_PROTO_MM,
    MAT_EXTRACT_DIR,
    TIME_SCALE_FROUDE,
    all_mat_names_for_plot,
    build_mat_run_catalog,
    mat_os_plots_dir,
    run_os_plots_dir,
)
from paths import HERE
from gm_duration import arias_significant_duration


# ------------------------------------------------------------
# 1. PATHS, UNITS, AND PLOT KNOBS
# ------------------------------------------------------------

REPO = HERE.parent

DPI = 300

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

# Friendly legend names for stateOS leaves. typeConv3/1 labels are interpretive
# (not from an official OpenFresco map); y-axis uses the raw leaf name.
# Field guide: plot/lab/STATEOS_SIGNALS.md (Seki §2.1).
STATE_LABELS = {
    "typeConv3": "typeConv3",
    "typeConv1": "typeConv1",
    "s1": "flag s1",
    "s2": "flag s2",
    "s3": "flag s3",
}

HELP = """\
usage: python plot/PlotMatOS.py [matFile ...]

  no args   every mapped mat in TestMatrix_lab_runs.csv + orphan registry
  matFile   one or more .mat names or stems
  needs     SyncLabBackup extracts in OSU_SSI_BRIDGE_DATA_LOCAL/.../mat_extract/
  per-run   hist_os_com.png, hist_os_state.png, hist_os_tar_com_mea.png
  note      group pier compares: plot/PlotEQComparePairs.py (pairs/)
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


def dual_time_xaxis(ax) -> None:
    """
    Add model time below and Froude-scaled prototype time above.

    Args:    ax  axes plotted against model time (s)
    Returns: none (updates ax)
    """
    ax.set_xlabel(r"$t_\mathrm{lab}$ (s) model scale")
    prototype_axis = ax.secondary_xaxis(
        "top",
        functions=(
            lambda time_model: time_model * TIME_SCALE_FROUDE,
            lambda time_proto: time_proto / TIME_SCALE_FROUDE,
        ),
    )
    prototype_axis.set_xlabel(
        r"$t_\mathrm{lab}\sqrt{\lambda}$ (s) prototype scale"
    )


def d595_lab_window() -> tuple[float, float] | None:
    """
    D5–95 on the Simulink lab clock (GM prototype ÷ √λ).

    Returns: (t5_lab_s, t95_lab_s) or None
    """
    try:
        duration = arias_significant_duration()
    except (OSError, ValueError) as exc:
        print(f"PlotMatOS: D5-95 unavailable ({exc})")
        return None
    return (
        float(duration.t5_s) / TIME_SCALE_FROUDE,
        float(duration.t95_s) / TIME_SCALE_FROUDE,
    )


def finish_lab_full_zoom(
    ax_f,
    ax_z,
    d595_lab: tuple[float, float] | None,
    *,
    dual_time: bool = True,
    zoom_title: bool = True,
) -> None:
    """Grid + dual time; xlim zoom to lab-mapped D5–95."""
    for ax in (ax_f, ax_z):
        ax.grid(True, ls=":", alpha=0.45)
        if dual_time:
            dual_time_xaxis(ax)
    if d595_lab is not None:
        ax_z.set_xlim(d595_lab[0], d595_lab[1])
        if zoom_title:
            ax_z.set_title(r"D5–95 zoom", fontsize=10, pad=4)


def integer_yaxis(ax, y: np.ndarray | None = None) -> None:
    """
    Restrict stateOS major ticks to integers (OpenFresco flags / counts).

    Args:    ax  stateOS axes; y  optional series used to snap ylim
    Returns: none (updates ax)
    """
    if y is not None and y.size:
        y_min = int(np.floor(float(np.nanmin(y))))
        y_max = int(np.ceil(float(np.nanmax(y))))
        if y_max <= y_min:
            y_max = y_min + 1
        pad = 0.15
        ax.set_ylim(y_min - pad, y_max + pad)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, min_n_ticks=2))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(round(v))}"))


# ------------------------------------------------------------
# 3. EXTRACT I/O AND SERIES DISCOVERY
# ------------------------------------------------------------


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
    Choose runs/<dump>/os/ for a mapped mat, else mats/<stem>/os/.

    Args:    mat_name  mat file name; mmap  build_mat_run_catalog()
    Returns: output directory
    """
    stem = Path(mat_name).name
    mat_info = mmap.get("mats", {}).get(stem)
    run_name = mat_info.get("run") if mat_info else None
    if run_name:
        return run_os_plots_dir(run_name)
    return mat_os_plots_dir(Path(stem).stem)


def iter_mat_npz(
    mmap: dict,
    mats: list[str] | None = None,
) -> list[tuple[str, Path]]:
    """
    Find available npz extracts for requested or mapped mats.

    Args:    mmap  mat-to-run mapping; mats  requested mat names or None
    Returns: list of (mat_file_name, npz_path)
    """
    mat_names: list[str] = list(mats) if mats else all_mat_names_for_plot()
    if mats is None and MAT_EXTRACT_DIR.is_dir():
        skip = set(mmap.get("duplicate_mats", {}))
        for npz_path in sorted(MAT_EXTRACT_DIR.glob("*.npz")):
            mat_name = f"{npz_path.stem}.mat"
            if mat_name not in mat_names and mat_name not in skip:
                mat_names.append(mat_name)

    pairs: list[tuple[str, Path]] = []
    for mat_name in mat_names:
        npz_path = MAT_EXTRACT_DIR / f"{Path(mat_name).stem}.npz"
        if npz_path.is_file():
            pairs.append((Path(mat_name).name, npz_path))
    return pairs


# ------------------------------------------------------------
# 4. RUN LABELS
# ------------------------------------------------------------


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


def plot_com(
    z: np.lib.npyio.NpzFile,
    out: Path,
    title: str,
    d595_lab: tuple[float, float] | None = None,
) -> None:
    """
    Plot the first comSigOS command displacement (full | D5–95 zoom).

    Args:    z  npz extract; out  PNG path; title  figure title
             d595_lab  lab-clock (t5, t95) or None
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
    y = data[:, signal_column] * DISP_M_TO_PROTO_MM
    label = rf"command ($\times\lambda={CYLINDER_LENGTH_SCALE:g}$)"

    if d595_lab is None:
        fig, ax = plt.subplots(figsize=(11.0, 4.4), constrained_layout=True)
        ax.plot(time_model_s, y, color=COLORS[1], lw=1.05, label=label)
        dual_time_xaxis(ax)
        ax.set_ylabel(r"command $u$ (mm, prototype)")
        ax.set_title(title)
        ax.legend(fontsize=8, loc="best")
        ax.grid(True, ls=":", alpha=0.45)
    else:
        fig, (ax_f, ax_z) = plt.subplots(
            1,
            2,
            figsize=(13.2, 4.4),
            sharey=True,
            gridspec_kw={"width_ratios": [1.35, 1.0], "wspace": 0.22},
            constrained_layout=True,
        )
        for ax in (ax_f, ax_z):
            ax.plot(time_model_s, y, color=COLORS[1], lw=1.05, label=label)
        finish_lab_full_zoom(ax_f, ax_z, d595_lab)
        ax_f.set_ylabel(r"command $u$ (mm, prototype)")
        ax_f.set_title(title, fontsize=10)
        ax_f.legend(fontsize=8, loc="best")

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(f"PlotMatOS: wrote {out}")


def plot_tar_com_mea(
    z: np.lib.npyio.NpzFile,
    out: Path,
    title: str,
    d595_lab: tuple[float, float] | None = None,
) -> None:
    """
    Plot first target, command, and measurement *OS (full | D5–95 zoom).

    Args:    z  npz extract; out  PNG path; title  figure title
             d595_lab  lab-clock (t5, t95) or None
    Returns: none (writes PNG when at least one *OS block exists)
    """
    series: list[tuple[np.ndarray, np.ndarray, str, str]] = []
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
        series.append(
            (
                time_model_s,
                data[:, signal_column] * DISP_M_TO_PROTO_MM,
                color,
                f"{block_key}: {signal_name}",
            )
        )

    if not series:
        print(f"PlotMatOS: skip tar/com/mea ({title})")
        return

    ylab = rf"displacement (mm, prototype, $\times\lambda={CYLINDER_LENGTH_SCALE:g}$)"

    def _draw(ax) -> None:
        for t, y, color, lab in series:
            ax.plot(t, y, color=color, lw=1.0, label=lab)

    if d595_lab is None:
        fig, ax = plt.subplots(figsize=(11.0, 4.6), constrained_layout=True)
        _draw(ax)
        dual_time_xaxis(ax)
        ax.set_ylabel(ylab)
        ax.set_title(title)
        ax.legend(fontsize=8, loc="best")
        ax.grid(True, ls=":", alpha=0.45)
    else:
        fig, (ax_f, ax_z) = plt.subplots(
            1,
            2,
            figsize=(13.2, 4.6),
            sharey=True,
            gridspec_kw={"width_ratios": [1.35, 1.0], "wspace": 0.22},
            constrained_layout=True,
        )
        _draw(ax_f)
        _draw(ax_z)
        finish_lab_full_zoom(ax_f, ax_z, d595_lab)
        ax_f.set_ylabel(ylab)
        ax_f.set_title(title, fontsize=10)
        ax_f.legend(fontsize=8, loc="best")

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(f"PlotMatOS: wrote {out}")


def plot_state(
    z: np.lib.npyio.NpzFile,
    out: Path,
    title: str,
    d595_lab: tuple[float, float] | None = None,
) -> None:
    """
    Plot stateOS (full | D5–95 zoom): simState above remaining flags.

    Args:    z  npz extract; out  PNG path; title  figure title
             d595_lab  lab-clock (t5, t95) or None
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

    sim_state_column = state_columns[0]
    type_conv3 = np.rint(data[:, sim_state_column])
    flag_series: list[tuple[np.ndarray, str, str]] = []
    for color_index, state_column in enumerate(state_columns[1:]):
        flag_series.append(
            (
                np.rint(data[:, state_column]),
                COLORS[(color_index + 1) % len(COLORS)],
                state_label(names[state_column]),
            )
        )

    if d595_lab is None:
        fig, axes = plt.subplots(
            2,
            1,
            figsize=(11.0, 6.2),
            sharex=True,
            constrained_layout=True,
            gridspec_kw={"height_ratios": [2.0, 1.2]},
        )
        state_ax, flags_ax = axes
        state_ax.plot(
            time_model_s,
            type_conv3,
            color=COLORS[0],
            lw=1.05,
            drawstyle="steps-post",
            label=state_label(names[sim_state_column]),
        )
        state_ax.set_ylabel(r"typeConv3")
        state_ax.set_title(title)
        state_ax.legend(fontsize=8, loc="best")
        state_ax.grid(True, ls=":", alpha=0.45)
        integer_yaxis(state_ax, type_conv3)
        flag_stack: list[np.ndarray] = []
        for flag_y, color, lab in flag_series:
            flag_stack.append(flag_y)
            flags_ax.plot(
                time_model_s,
                flag_y,
                color=color,
                lw=0.95,
                drawstyle="steps-post",
                label=lab,
            )
        flags_ax.set_ylabel("flags / count")
        flags_ax.legend(fontsize=7.5, loc="best", ncol=2)
        flags_ax.grid(True, ls=":", alpha=0.45)
        if flag_stack:
            integer_yaxis(flags_ax, np.concatenate(flag_stack))
        else:
            integer_yaxis(flags_ax)
        dual_time_xaxis(flags_ax)
    else:
        fig, axes = plt.subplots(
            2,
            2,
            figsize=(13.2, 6.4),
            sharex="col",
            sharey="row",
            constrained_layout=True,
            gridspec_kw={
                "height_ratios": [2.0, 1.2],
                "width_ratios": [1.35, 1.0],
                "wspace": 0.22,
            },
        )
        for col, ax_pair in enumerate((axes[:, 0], axes[:, 1])):
            state_ax, flags_ax = ax_pair
            state_ax.plot(
                time_model_s,
                type_conv3,
                color=COLORS[0],
                lw=1.05,
                drawstyle="steps-post",
                label=state_label(names[sim_state_column]),
            )
            for flag_y, color, lab in flag_series:
                flags_ax.plot(
                    time_model_s,
                    flag_y,
                    color=color,
                    lw=0.95,
                    drawstyle="steps-post",
                    label=lab,
                )
            integer_yaxis(state_ax, type_conv3)
            if flag_series:
                integer_yaxis(
                    flags_ax, np.concatenate([fy for fy, _, _ in flag_series])
                )
            else:
                integer_yaxis(flags_ax)
        finish_lab_full_zoom(axes[0, 0], axes[0, 1], d595_lab, dual_time=False)
        finish_lab_full_zoom(
            axes[1, 0], axes[1, 1], d595_lab, dual_time=True, zoom_title=False
        )
        # sharex=col already synced xlim; re-apply zoom on bottom zoom panel
        if d595_lab is not None:
            axes[0, 1].set_xlim(d595_lab[0], d595_lab[1])
            axes[1, 1].set_xlim(d595_lab[0], d595_lab[1])
        axes[0, 0].set_ylabel(r"typeConv3")
        axes[1, 0].set_ylabel("flags / count")
        axes[0, 0].set_title(title, fontsize=10)
        axes[0, 0].legend(fontsize=8, loc="best")
        axes[1, 0].legend(fontsize=7.5, loc="best", ncol=2)

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
    d595_lab = d595_lab_window()
    plot_com(
        extract,
        output_dir / "hist_os_com.png",
        f"{label} — comSigOS (command $u$, prototype)",
        d595_lab=d595_lab,
    )
    plot_state(
        extract,
        output_dir / "hist_os_state.png",
        f"{label} — stateOS (OpenFresco; slowdowns)",
        d595_lab=d595_lab,
    )
    plot_tar_com_mea(
        extract,
        output_dir / "hist_os_tar_com_mea.png",
        f"{label} — tar / com / mea (*OS, prototype)",
        d595_lab=d595_lab,
    )


# ------------------------------------------------------------
# 6. CLI
# ------------------------------------------------------------


def main() -> int:
    """
    Parse mat names and write per-mat *OS PNGs.

    Args:    sys.argv  optional mat file names or stems
    Returns: process exit code
    """
    if any(argument in ("-h", "--help") for argument in sys.argv[1:]):
        print(HELP, end="")
        return 0

    mat_map = build_mat_run_catalog()
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
