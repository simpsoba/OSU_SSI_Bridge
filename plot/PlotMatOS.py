#!/usr/bin/env python3
"""Plot Simulink *OS channels vs model and prototype time (dual x-axis).

  python plot/PlotMatOS.py
  python plot/PlotMatOS.py 0821_GusBridge_rowPos1.mat ...

Reads extracts under OSU_SSI_BRIDGE_DATA_LOCAL/opensees_data/mat_extract/
(from SyncLabBackup.py). Uses plot/opensees_data/mat_run_map.json for run folders.

Clocks
------
  Lab Time in the mat = real time at **model** scale (DAQ; denser than OpenSees).
  Prototype time = Time × √λ; length scale λ = CYLINDER_LENGTH_SCALE (2.4) in lab_paths.
  Actuator / *OS displacements in mats are model scale → plot as × λ (prototype cm).
  Per-run plots: bottom x = model; top x = prototype. Same samples; not step pairing.
  Compare overlays: x = lab Time × √2.4 (prototype scale), matching hist_ux_*_realtime.

Writes (per mat / run):
  OSU_SSI_BRIDGE_DATA_LOCAL/plots/<run>/hist_os_com.png
  OSU_SSI_BRIDGE_DATA_LOCAL/plots/<run>/hist_os_state.png
  OSU_SSI_BRIDGE_DATA_LOCAL/plots/<run>/hist_os_tar_com_mea.png
  Unmapped mats → …/plots/mat_os/<stem>/

Compare overlays (grouped by model knobs from TestMatrix_lab_runs.csv):
  OSU_SSI_BRIDGE_DATA_LOCAL/plots/compare/<group>/hist_os_com.png
  OSU_SSI_BRIDGE_DATA_LOCAL/plots/compare/<group>/hist_os_state.png
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
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

REPO = HERE.parent
PLOTS_ROOT = plots_root()
DPI = 300
MAX_PLOT_PTS = 40_000
FOLDER_RE = re.compile(r"^r([+-]?\d+)_(\d{8})_(\d{4})_", re.IGNORECASE)
# Model-scale mat disp (m) → prototype cm: × λ × 100
DISP_M_TO_PROTO_CM = CYLINDER_LENGTH_SCALE * 100.0
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

# Friendly names for stateOS columns (OpenFresco typeConv* in this campaign).
STATE_LABELS = {
    "typeConv3": "simState",
    "typeConv1": "substep/count",
    "s1": "flag s1",
    "s2": "flag s2",
    "s3": "flag s3",
}


def short_sig(name: str) -> str:
    s = name.replace("\\", "/").strip()
    if "/" in s:
        s = s.rsplit("/", 1)[-1]
    return s


def state_label(name: str) -> str:
    s = short_sig(name)
    return STATE_LABELS.get(s, s)


def decimate(t: np.ndarray, y: np.ndarray, nmax: int = MAX_PLOT_PTS) -> tuple[np.ndarray, np.ndarray]:
    n = len(t)
    if n <= nmax:
        return t, y
    step = int(np.ceil(n / nmax))
    return t[::step], y[::step]


def apply_paper_style() -> None:
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
    if n <= 4:
        return n
    if n <= 10:
        return 4
    if n <= 18:
        return 5
    return 6


def place_legend_outside(ax: plt.Axes, n: int) -> None:
    if n <= 0:
        return
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return
    ncol = legend_ncol(n)
    nrows = int(np.ceil(n / ncol))
    y = -0.14 - 0.055 * max(0, nrows - 1)
    ax.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, y),
        ncol=ncol,
        fontsize=6.5,
        frameon=False,
        columnspacing=1.2,
        handlelength=1.45,
        handletextpad=0.35,
        borderaxespad=0.0,
    )


def dual_time_xaxis(ax) -> None:
    """Bottom = model time; top = prototype (×√2.4)."""
    ax.set_xlabel("t (s) model scale (lab real)")
    sec = ax.secondary_xaxis(
        "top",
        functions=(
            lambda tm: tm * TIME_SCALE_FROUDE,
            lambda tp: tp / TIME_SCALE_FROUDE,
        ),
    )
    sec.set_xlabel(r"t (s) prototype scale (lab real $\times\sqrt{2.4}$)")


def load_block(z: np.lib.npyio.NpzFile, key: str) -> tuple[np.ndarray, list[str], np.ndarray] | None:
    """Return (time_model, signal_names, data) for tarSigOS / comSigOS / …"""
    t_key = f"{key}_time"
    d_key = f"{key}_data"
    n_key = f"{key}_signalNames"
    if d_key not in z:
        return None
    data = np.asarray(z[d_key], dtype=float)
    names = [str(x) for x in z[n_key].tolist()] if n_key in z else []
    if t_key in z:
        t = np.asarray(z[t_key], dtype=float)
    else:
        ti = next(i for i, n in enumerate(names) if n.strip().lower() == "time")
        t = data[:, ti]
    return t, names, data


def non_time_cols(names: list[str], ncol: int) -> list[int]:
    ti = next((i for i, n in enumerate(names) if n.strip().lower() == "time"), None)
    return [i for i in range(ncol) if i != ti]


def out_dir_for_mat(mat_name: str, mmap: dict) -> Path:
    info = mmap.get("mats", {}).get(mat_name) or mmap.get("mats", {}).get(
        Path(mat_name).name
    )
    run = info.get("run") if info else None
    if run:
        return PLOTS_ROOT / run
    return PLOTS_ROOT / "mat_os" / Path(mat_name).stem


def compare_label(mat_name: str, mmap: dict) -> str:
    """Compact tag for compare legends."""
    info = mmap.get("mats", {}).get(Path(mat_name).name) or {}
    run = info.get("run")
    if not run:
        return Path(mat_name).stem
    m = FOLDER_RE.match(run)
    note = info.get("note") or ""
    np_s = ""
    if "np=" in note:
        np_s = note.split("np=", 1)[1].split()[0].rstrip(")")
    if m:
        hhmm = m.group(3)
        tag = f"r{m.group(1)} {hhmm[:2]}:{hhmm[2:]}"
        if np_s:
            tag += f" (np={np_s})"
        return tag
    return run


def plot_com(z: np.lib.npyio.NpzFile, out: Path, title: str) -> None:
    """Actuator command displacement (comSigOS)."""
    block = load_block(z, "comSigOS")
    if block is None:
        print(f"PlotMatOS: skip comSigOS ({title})")
        return
    t, names, data = block
    cols = non_time_cols(names, data.shape[1])
    if not cols:
        return
    j = cols[0]
    td, yd = decimate(t, data[:, j])
    fig, ax = plt.subplots(figsize=(11.0, 4.4), constrained_layout=True)
    ax.plot(
        td,
        yd * DISP_M_TO_PROTO_CM,
        color=COLORS[1],
        lw=1.05,
        label=rf"comSigOS $\times\lambda$ ($\lambda={CYLINDER_LENGTH_SCALE:g}$)",
    )
    dual_time_xaxis(ax)
    ax.set_ylabel(r"actuator command $u$ (cm, prototype)")
    ax.set_title(title)
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, ls=":", alpha=0.45)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(f"PlotMatOS: wrote {out}")


def plot_tar_com_mea(z: np.lib.npyio.NpzFile, out: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(11.0, 4.6), constrained_layout=True)
    plotted = 0
    for key, color in zip(("tarSigOS", "comSigOS", "meaSigOS"), COLORS):
        block = load_block(z, key)
        if block is None:
            continue
        t, names, data = block
        cols = non_time_cols(names, data.shape[1])
        if not cols:
            continue
        j = cols[0]
        label = short_sig(names[j]) if j < len(names) else key
        td, yd = decimate(t, data[:, j])
        ax.plot(
            td,
            yd * DISP_M_TO_PROTO_CM,
            color=color,
            lw=1.0,
            label=f"{key}: {label}",
        )
        plotted += 1
    if plotted == 0:
        plt.close(fig)
        print(f"PlotMatOS: skip tar/com/mea ({title})")
        return
    dual_time_xaxis(ax)
    ax.set_ylabel(rf"displacement (cm, prototype, $\times\lambda={CYLINDER_LENGTH_SCALE:g}$)")
    ax.set_title(title)
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, ls=":", alpha=0.45)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(f"PlotMatOS: wrote {out}")


def plot_state(z: np.lib.npyio.NpzFile, out: Path, title: str) -> None:
    """OpenFresco stateOS — simState / flags used to spot computational slowdowns."""
    block = load_block(z, "stateOS")
    if block is None:
        print(f"PlotMatOS: skip stateOS ({title})")
        return
    t, names, data = block
    cols = non_time_cols(names, data.shape[1])
    if not cols:
        return
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(11.0, 6.2),
        sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [2.0, 1.2]},
    )
    ax0, ax1 = axes
    # Top: simState (first non-Time col) — clearest slowdown indicator
    j0 = cols[0]
    td, yd = decimate(t, data[:, j0])
    ax0.plot(td, yd, color=COLORS[0], lw=1.05, drawstyle="steps-post", label=state_label(names[j0]))
    ax0.set_ylabel("simState")
    ax0.set_title(title)
    ax0.legend(fontsize=8, loc="best")
    ax0.grid(True, ls=":", alpha=0.45)
    # Bottom: remaining channels
    for k, j in enumerate(cols[1:]):
        td, yd = decimate(t, data[:, j])
        ax1.plot(
            td,
            yd,
            color=COLORS[(k + 1) % len(COLORS)],
            lw=0.95,
            drawstyle="steps-post",
            label=state_label(names[j]),
        )
    ax1.set_ylabel("flags / count")
    ax1.legend(fontsize=7.5, loc="best", ncol=2)
    ax1.grid(True, ls=":", alpha=0.45)
    dual_time_xaxis(ax1)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(f"PlotMatOS: wrote {out}")


def process_mat(mat_name: str, mmap: dict) -> None:
    stem = Path(mat_name).stem
    npz_path = MAT_EXTRACT_DIR / f"{stem}.npz"
    if not npz_path.is_file():
        print(f"PlotMatOS: missing extract {npz_path}", file=sys.stderr)
        return
    z = np.load(npz_path, allow_pickle=True)
    out_dir = out_dir_for_mat(mat_name, mmap)
    tag = series_label(mat_name, mmap)
    plot_com(
        z,
        out_dir / "hist_os_com.png",
        f"{tag} — comSigOS (actuator command, prototype $\\times\\lambda$)",
    )
    plot_state(
        z,
        out_dir / "hist_os_state.png",
        f"{tag} — stateOS (OpenFresco; slowdowns)",
    )
    plot_tar_com_mea(
        z,
        out_dir / "hist_os_tar_com_mea.png",
        f"{tag} — tar / com / mea (*OS, prototype)",
    )


def iter_mat_npz(mmap: dict, mats: list[str] | None = None) -> list[tuple[str, Path]]:
    names = mats or sorted(mmap.get("mats", {}).keys())
    if mats is None and MAT_EXTRACT_DIR.is_dir():
        for p in sorted(MAT_EXTRACT_DIR.glob("*.npz")):
            name = f"{p.stem}.mat"
            if name not in names:
                names.append(name)
    out: list[tuple[str, Path]] = []
    for name in names:
        npz = MAT_EXTRACT_DIR / f"{Path(name).stem}.npz"
        if npz.is_file():
            out.append((Path(name).name, npz))
    return out


def series_label(mat_name: str, mmap: dict) -> str:
    info = mmap.get("mats", {}).get(Path(mat_name).name) or {}
    run = info.get("run")
    note = info.get("note") or ""
    if run:
        tag = run
        if "np=" in note:
            tag += f"  ({note})"
        return tag
    return Path(mat_name).stem


def plot_compare_com(pairs: list[tuple[str, Path]], mmap: dict, out: Path) -> None:
    apply_paper_style()
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    plotted = 0
    for mat_name, npz in pairs:
        z = np.load(npz, allow_pickle=True)
        block = load_block(z, "comSigOS")
        if block is None:
            continue
        t, names, data = block
        cols = non_time_cols(names, data.shape[1])
        if not cols:
            continue
        t_proto = t * TIME_SCALE_FROUDE
        td, yd = decimate(t_proto, data[:, cols[0]])
        ax.plot(
            td,
            yd * DISP_M_TO_PROTO_CM,
            color=COLORS[plotted % len(COLORS)],
            lw=1.0,
            label=compare_label(mat_name, mmap),
        )
        plotted += 1
    if plotted == 0:
        plt.close(fig)
        print("PlotMatOS: no series for compare hist_os_com.png")
        return
    ax.set_xlabel(r"$t$ (s)")
    ax.set_ylabel(r"$u$ (cm)")
    place_legend_outside(ax, plotted)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=DPI, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    print(f"PlotMatOS: wrote {out}  ({plotted} series)")


def plot_compare_state(pairs: list[tuple[str, Path]], mmap: dict, out: Path) -> None:
    """One axes per run (shared x) — stateOS simState for slowdowns."""
    series: list[tuple[str, np.ndarray, np.ndarray]] = []
    for mat_name, npz in pairs:
        z = np.load(npz, allow_pickle=True)
        block = load_block(z, "stateOS")
        if block is None:
            continue
        t, names, data = block
        cols = non_time_cols(names, data.shape[1])
        if not cols:
            continue
        t_proto = t * TIME_SCALE_FROUDE
        td, yd = decimate(t_proto, data[:, cols[0]])
        series.append((series_label(mat_name, mmap), td, yd))
    if not series:
        print("PlotMatOS: no series for compare hist_os_state.png")
        return

    n = len(series)
    fig_h = max(6.0, 1.05 * n)
    fig, axes = plt.subplots(
        n,
        1,
        figsize=(11.0, fig_h),
        sharex=True,
        constrained_layout=True,
    )
    if n == 1:
        axes = [axes]
    for i, (lab, td, yd) in enumerate(series):
        ax = axes[i]
        ax.plot(
            td,
            yd,
            color=COLORS[i % len(COLORS)],
            lw=1.0,
            drawstyle="steps-post",
        )
        ax.set_ylabel("simState", fontsize=7)
        ax.grid(True, ls=":", alpha=0.45)
        # Short run tag inside the panel (full legend would dominate)
        ax.text(
            0.01,
            0.92,
            lab,
            transform=ax.transAxes,
            fontsize=6.5,
            va="top",
            ha="left",
            color=COLORS[i % len(COLORS)],
            bbox={"facecolor": "white", "alpha": 0.7, "edgecolor": "none", "pad": 1.5},
        )
    axes[0].set_title("stateOS simState — OpenFresco (slowdowns)")
    axes[-1].set_xlabel(r"t (s) lab real $\times\sqrt{2.4}$ (prototype scale)")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(f"PlotMatOS: wrote {out}  ({n} panels)")


HELP = """\
usage: python plot/PlotMatOS.py [matFile ...]

  no args   every .mat in mat_run_map.json (plus any *.npz extract)
  matFile   one or more .mat names or stems
  needs     SyncLabBackup extracts in OSU_SSI_BRIDGE_DATA_LOCAL/.../mat_extract/
  per-run   hist_os_com.png, hist_os_state.png, hist_os_tar_com_mea.png
  compare   OSU_SSI_BRIDGE_DATA_LOCAL/plots/compare/<group>/hist_os_*.png
"""


def _mat_group_slug(mat_name: str, mmap: dict, dump_group: dict[str, str]) -> str:
    info = mmap.get("mats", {}).get(Path(mat_name).name) or {}
    run = info.get("run")
    if run and run in dump_group:
        return dump_group[run]
    return "unmapped"


def write_grouped_compares(
    pairs: list[tuple[str, Path]], mmap: dict
) -> None:
    dump_group = dump_to_group()
    buckets: dict[str, list[tuple[str, Path]]] = {}
    for mat_name, npz in pairs:
        slug = _mat_group_slug(mat_name, mmap, dump_group)
        buckets.setdefault(slug, []).append((mat_name, npz))
    for slug, sub in buckets.items():
        out_dir = PLOTS_ROOT / "compare" / slug
        print(f"PlotMatOS: compare group {slug}  ({len(sub)} mat(s))")
        plot_compare_com(sub, mmap, out_dir / "hist_os_com.png")
        plot_compare_state(sub, mmap, out_dir / "hist_os_state.png")


def main() -> int:
    if any(a in ("-h", "--help") for a in sys.argv[1:]):
        print(HELP, end="")
        return 0
    mmap = load_mat_run_map()
    if len(sys.argv) > 1:
        mats = [
            Path(a).name if a.lower().endswith(".mat") else f"{Path(a).stem}.mat"
            for a in sys.argv[1:]
        ]
    else:
        mats = None
    pairs = iter_mat_npz(mmap, mats)
    if not pairs:
        print("PlotMatOS: no mats to plot", file=sys.stderr)
        return 1
    print(f"PlotMatOS: {len(pairs)} mat(s)  TIME_SCALE_FROUDE={TIME_SCALE_FROUDE:.4f}")
    for mat_name, _ in pairs:
        process_mat(mat_name, mmap)
    write_grouped_compares(pairs, mmap)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
