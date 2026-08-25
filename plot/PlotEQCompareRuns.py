#!/usr/bin/env python3
"""Overlay pier-top ux histories, grouped by model knobs (not solver/np).

  python plot/PlotEQCompareRuns.py
  python plot/PlotEQCompareRuns.py <runDir> <runDir> ...

Default: every dump listed in TestMatrix_lab_runs.csv that has pier_top_disp.

Writes under OSU_SSI_BRIDGE_DATA_LOCAL/plots/compare/<group>/ :
  hist_ux_complete.png / hist_ux_all.png              # OpenSees numerical t
  hist_ux_complete_realtime.png / hist_ux_all_realtime.png  # meaSigOS

Groups = same soilMesh, soilProfile, soilEleType, soilConstitutive,
expElementType, holdPierON, and non-default Rayleigh ξ. Solver / np / integrator
variants stay in one folder (fair overlay of the same physical model).

Clocks / scales (lab_paths)
---------------------------
  λ = CYLINDER_LENGTH_SCALE (2.4).
  Lab mat Time = real time at **model** scale.
  Prototype time = Time × √λ; prototype disp = model × λ.
  OpenSees pier recorder = numerical time / disp at **prototype** scale.

Realtime plots use Simulink meaSigOS (actuator feedback) — no OS interp.
Primary axes = prototype; secondary = model.
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

REPO = HERE.parent
PLOTS_ROOT = plots_root()
DPI = 300
TREC_TOL = 1.0  # s; match PlotEQ.truncated_end
MAX_PLOT_PTS = 40_000
# Skip recorder dumps that diverged (short failed runs)
UX_ABS_MAX_CM = 100.0
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
    """Columns for an outside legend: keep rows short without crowding."""
    if n <= 4:
        return n
    if n <= 10:
        return 4
    if n <= 18:
        return 5
    return 6


def place_legend_outside(fig: plt.Figure, ax: plt.Axes, n: int) -> None:
    """Legend below the axes (axes coords); multi-column. Relies on bbox_inches=tight."""
    if n <= 0:
        return
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return
    ncol = legend_ncol(n)
    nrows = int(np.ceil(n / ncol))
    # Negative y in axes fraction: one row ≈ 0.055 below the x-label band
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


def loadtxt_partial(path: Path) -> np.ndarray:
    try:
        a = np.loadtxt(path)
    except ValueError:
        rows: list[list[float]] = []
        ncol: int | None = None
        with path.open() as f:
            for ln in f:
                s = ln.strip()
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
        a = np.asarray(rows, dtype=float)
    if a.size == 0:
        return np.empty((0, 0))
    if a.ndim == 1:
        a = a.reshape(1, -1)
    return a


def find_pier_top(eq: Path) -> Path | None:
    serial = eq / "pier_top_disp.out"
    if serial.is_file():
        return serial
    shards = sorted(eq.glob("pier_top_disp.out.*"))
    if not shards:
        return None
    for p in shards:
        if p.name.endswith(".0"):
            return p
    return max(shards, key=lambda p: p.stat().st_size)


def read_meta(eq: Path) -> dict[str, str]:
    meta: dict[str, str] = {}
    for name in ("window_meta.txt.0", "window_meta.txt"):
        p = eq / name
        if not p.is_file():
            continue
        for ln in p.read_text().splitlines():
            s = ln.strip()
            if not s or s.startswith("#"):
                continue
            k, _, rest = s.partition(" ")
            meta[k] = rest.strip()
        break
    return meta


def run_duration(eq: Path) -> tuple[float | None, float, bool]:
    """(t_last, Trec, is_complete). Complete := t_last >= Trec - tol when Trec>0."""
    meta = read_meta(eq)
    try:
        trec = float(meta.get("Trec", 0) or 0)
    except ValueError:
        trec = 0.0
    pier = find_pier_top(eq)
    if pier is None:
        return None, trec, False
    a = loadtxt_partial(pier)
    if a.size == 0:
        return None, trec, False
    t_last = float(a[-1, 0])
    complete = trec > 0.0 and t_last >= trec - TREC_TOL
    return t_last, trec, complete


def short_label(eq: Path, incomplete: bool = False) -> str:
    """Compact run tag for legends: r±NN HH:MM (np=N)[*]."""
    meta = read_meta(eq)
    np_s = meta.get("np", "?")
    m = FOLDER_RE.match(eq.name)
    if m:
        hhmm = m.group(3)
        tag = f"r{m.group(1)} {hhmm[:2]}:{hhmm[2:]} (np={np_s})"
    else:
        tag = f"{eq.name} (np={np_s})"
    if incomplete:
        tag += "*"
    return tag


def discover_runs() -> list[Path]:
    root = resolve_opensees_data()
    if root is None:
        return []
    out: list[Path] = []
    for p in sorted(root.iterdir()):
        if p.is_dir() and find_pier_top(p) is not None:
            out.append(p)
    return out


def run_to_mat_name() -> dict[str, str]:
    """OpenSees run folder name -> Simulink .mat file name."""
    mmap = load_mat_run_map()
    out: dict[str, str] = {}
    for mat_name, info in mmap.get("mats", {}).items():
        run = info.get("run")
        if run:
            out[run] = mat_name
    return out


def decimate(t: np.ndarray, y: np.ndarray, nmax: int = MAX_PLOT_PTS) -> tuple[np.ndarray, np.ndarray]:
    n = len(t)
    if n <= nmax:
        return t, y
    step = int(np.ceil(n / nmax))
    return t[::step], y[::step]


def load_mat_mea_feedback(mat_name: str) -> tuple[np.ndarray, np.ndarray] | None:
    """Lab Time (s, model) and meaSigOS primary (m, model) from extract."""
    npz = MAT_EXTRACT_DIR / f"{Path(mat_name).stem}.npz"
    if not npz.is_file():
        return None
    z = np.load(npz, allow_pickle=True)
    if "meaSigOS_time" in z and "meaSigOS_primary" in z:
        t = np.asarray(z["meaSigOS_time"], dtype=float)
        u = np.asarray(z["meaSigOS_primary"], dtype=float)
        return t, u
    if "meaSigOS_data" not in z:
        return None
    data = np.asarray(z["meaSigOS_data"], dtype=float)
    names = [str(x) for x in z["meaSigOS_signalNames"].tolist()]
    ti = next(i for i, n in enumerate(names) if n.strip().lower() == "time")
    j = next(i for i in range(data.shape[1]) if i != ti)
    return data[:, ti], data[:, j]


def dual_proto_model_axes(ax) -> None:
    """Primary = prototype; secondary (top/right) = model."""
    lam = CYLINDER_LENGTH_SCALE
    s = TIME_SCALE_FROUDE
    ax.set_xlabel(r"$t$ (s)")
    ax.set_ylabel(r"$\Delta u$ (cm)")
    secx = ax.secondary_xaxis(
        "top",
        functions=(lambda tp: tp / s, lambda tm: tm * s),
    )
    secx.set_xlabel(r"$t_\mathrm{m}$ (s)")
    secy = ax.secondary_yaxis(
        "right",
        functions=(lambda up: up / lam, lambda um: um * lam),
    )
    secy.set_ylabel(r"$\Delta u_\mathrm{m}$ (cm)")


def plot_compare_opensees(
    runs: list[Path],
    out: Path,
    mark_incomplete: bool,
) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    apply_paper_style()
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    t_eq_marks: list[float] = []
    plotted = 0
    for eq in runs:
        pier = find_pier_top(eq)
        if pier is None:
            print(f"PlotEQCompareRuns: skip {eq.name} (no pier_top_disp)")
            continue
        a = loadtxt_partial(pier)
        if a.size == 0 or a.shape[1] < 2:
            print(f"PlotEQCompareRuns: skip {eq.name} (empty pier_top)")
            continue
        t_last, trec, complete = run_duration(eq)
        t = a[:, 0]
        ux = a[:, 1] - a[0, 1]
        if not np.all(np.isfinite(ux)) or float(np.nanmax(np.abs(ux))) * 100.0 > UX_ABS_MAX_CM:
            print(f"PlotEQCompareRuns: skip {eq.name} (non-physical pier ux)")
            continue
        color = COLORS[plotted % len(COLORS)]
        ls = "-" if complete else "--"
        ax.plot(
            t,
            ux * 100.0,
            color=color,
            lw=1.0 if complete else 0.9,
            ls=ls,
            label=short_label(eq, incomplete=mark_incomplete and not complete),
        )
        plotted += 1
        meta = read_meta(eq)
        try:
            fv = float(meta.get("freeVibT", 0) or 0)
            if trec > fv > 0:
                t_eq_marks.append(trec - fv)
        except ValueError:
            pass
        if not complete and t_last is not None:
            ax.axvline(float(t_last), color=color, lw=0.6, ls=":", alpha=0.5)
    if t_eq_marks:
        te = float(np.median(t_eq_marks))
        ax.axvline(te, color="0.45", lw=0.8, ls="--", label="EQ end")
    ax.set_xlabel(r"$t$ (s)")
    ax.set_ylabel(r"$\Delta u_x$ (cm)")
    n_leg = plotted + (1 if t_eq_marks else 0)
    place_legend_outside(fig, ax, n_leg)
    fig.savefig(out, dpi=DPI, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    print(f"PlotEQCompareRuns: wrote {out}  ({plotted} series)")


def plot_compare_realtime(
    runs: list[Path],
    out: Path,
    mark_incomplete: bool,
    run_mat: dict[str, str],
) -> None:
    """meaSigOS actuator feedback vs lab Time; dual proto/model axes. No OS interp."""
    out.parent.mkdir(parents=True, exist_ok=True)
    apply_paper_style()
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    t_eq_marks: list[float] = []
    plotted = 0
    lam = CYLINDER_LENGTH_SCALE
    for eq in runs:
        mat_name = run_mat.get(eq.name)
        if not mat_name:
            print(f"PlotEQCompareRuns: skip realtime {eq.name} (no mat)")
            continue
        pair = load_mat_mea_feedback(mat_name)
        if pair is None:
            print(f"PlotEQCompareRuns: skip realtime {eq.name} (no meaSigOS)")
            continue
        t_lab, u_lab = pair
        if t_lab.size == 0:
            continue
        _, trec, complete = run_duration(eq)
        t_proto = t_lab * TIME_SCALE_FROUDE
        u_proto_cm = (u_lab - u_lab[0]) * lam * 100.0
        td, ud = decimate(t_proto, u_proto_cm)
        color = COLORS[plotted % len(COLORS)]
        ls = "-" if complete else "--"
        ax.plot(
            td,
            ud,
            color=color,
            lw=1.0 if complete else 0.9,
            ls=ls,
            label=short_label(eq, incomplete=mark_incomplete and not complete),
        )
        plotted += 1
        meta = read_meta(eq)
        try:
            fv = float(meta.get("freeVibT", 0) or 0)
            if trec > fv > 0:
                t_eq_marks.append(trec - fv)
        except ValueError:
            pass
        if not complete:
            ax.axvline(
                float(t_lab[-1] * TIME_SCALE_FROUDE),
                color=color,
                lw=0.6,
                ls=":",
                alpha=0.5,
            )
    if plotted == 0:
        print(f"PlotEQCompareRuns: no series for {out.name}")
        plt.close(fig)
        return
    if t_eq_marks:
        te = float(np.median(t_eq_marks))
        ax.axvline(te, color="0.45", lw=0.8, ls="--", label="EQ end")
    dual_proto_model_axes(ax)
    n_leg = plotted + (1 if t_eq_marks else 0)
    place_legend_outside(fig, ax, n_leg)
    fig.savefig(out, dpi=DPI, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    print(f"PlotEQCompareRuns: wrote {out}  ({plotted} series)")


HELP = """\
usage: python plot/PlotEQCompareRuns.py [runDir ...]

  no args   dumps from TestMatrix_lab_runs.csv (grouped by model knobs)
  runDir    one or more dump folders (still grouped if in lab_runs)
  output    OSU_SSI_BRIDGE_DATA_LOCAL/plots/compare/<group>/hist_ux_*.png
  groups    soilMesh + soilProfile + soilEleType + constitutive +
            expElement + hold + non-default ξ  (not solver/np)
  complete  t_last >= Trec - 1 s
  realtime  meaSigOS (actuator feedback); dual proto/model axes
"""


def resolve_run_path(dump_name: str, root: Path | None) -> Path | None:
    if root is not None:
        p = root / dump_name
        if p.is_dir():
            return p
    local = REPO / "OSU_SSI_BRIDGE_DATA_LOCAL" / "opensees_data" / dump_name
    return local if local.is_dir() else None


def write_group_plots(
    slug: str,
    runs: list[Path],
    run_mat: dict[str, str],
    label: str,
) -> None:
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

    if len(sys.argv) > 1:
        # CLI dumps: bucket by lab_runs group when known, else "cli"
        buckets: dict[str, list[Path]] = {}
        labels: dict[str, str] = {}
        for a in sys.argv[1:]:
            p = Path(a).resolve()
            if not p.is_dir():
                print(f"PlotEQCompareRuns: skip (not a dir) {a}", file=sys.stderr)
                continue
            slug = dump_group.get(p.name, "cli")
            buckets.setdefault(slug, []).append(p)
            if slug not in labels:
                labels[slug] = slug
        if not buckets:
            print("PlotEQCompareRuns: no valid run dirs", file=sys.stderr)
            return 1
        for slug, runs in buckets.items():
            write_group_plots(slug, runs, run_mat, labels[slug])
        return 0

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

    n_written = 0
    for slug, rows in lab_groups.items():
        runs: list[Path] = []
        for row in rows:
            p = resolve_run_path(row["DumpFolder"], root)
            if p is None:
                print(f"PlotEQCompareRuns: missing dump {row['DumpFolder']}")
                continue
            if find_pier_top(p) is None:
                print(f"PlotEQCompareRuns: no pier_top {p.name}")
                continue
            runs.append(p)
        if not runs:
            continue
        write_group_plots(slug, runs, run_mat, group_label(rows[0]))
        n_written += 1
    print(f"\nPlotEQCompareRuns: {n_written} group folder(s) under compare/")
    return 0 if n_written else 1


if __name__ == "__main__":
    raise SystemExit(main())
