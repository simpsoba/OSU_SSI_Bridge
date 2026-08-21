#!/usr/bin/env python3
"""Overlay pier-top ux histories from lab dump folders.

  python plot/PlotEQCompareRuns.py
  python plot/PlotEQCompareRuns.py <runDir> <runDir> ...

Default: every subfolder under the Drive `opensees data` shortcut that has
`pier_top_disp.out` / `pier_top_disp.out.*`.

Writes (Δux from t0):
  OSU_SSI_PLOTS/compare/hist_ux_complete.png   # t_last >= Trec - 1 s
  OSU_SSI_PLOTS/compare/hist_ux_all.png        # complete + incomplete
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from paths import HERE

REPO = HERE.parent
PLOTS_ROOT = REPO / "OSU_SSI_PLOTS"
DRIVE_ROOT = REPO / "OSU_SSI_BRIDGE_DATA"
DPI = 140
TREC_TOL = 1.0  # s; match PlotEQ.truncated_end
COLORS = (
    "#1565c0",
    "#c45c12",
    "#2e7d32",
    "#6a1b9a",
    "#00838f",
    "#ad1457",
    "#5d4037",
    "#455a64",
)


def resolve_opensees_data() -> Path | None:
    """Folder that actually holds run dirs (often behind a Drive .lnk)."""
    direct = DRIVE_ROOT / "opensees data"
    if direct.is_dir():
        return direct
    shortcut_root = Path(r"G:\.shortcut-targets-by-id")
    if shortcut_root.is_dir():
        for p in sorted(shortcut_root.glob("*/opensees data")):
            if p.is_dir():
                return p
    return None


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
    meta = read_meta(eq)
    np_s = meta.get("np", "?")
    tag = f"{eq.name}  (np={np_s})"
    if incomplete:
        t_last, trec, _ = run_duration(eq)
        if t_last is not None and trec > 0:
            tag += f"  [inc {t_last:.0f}/{trec:.0f}s]"
        else:
            tag += "  [incomplete]"
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


def plot_compare(
    runs: list[Path],
    out: Path,
    title: str,
    mark_incomplete: bool,
) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(11.0, 4.6), constrained_layout=True)
    t_eq_marks: list[float] = []
    plotted = 0
    for i, eq in enumerate(runs):
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
        color = COLORS[plotted % len(COLORS)]
        ls = "-" if complete else "--"
        lw = 1.15 if complete else 1.0
        ax.plot(
            t,
            ux * 100.0,
            color=color,
            lw=lw,
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
            ax.axvline(float(t_last), color=color, lw=0.7, ls=":", alpha=0.55)
    if t_eq_marks:
        te = float(np.median(t_eq_marks))
        ax.axvline(te, color="#78909c", lw=1.0, ls=":", label=f"EQ end ~{te:.0f} s")
    ax.set_xlabel("t (s)")
    ax.set_ylabel(r"pier top $\Delta u_x$ (cm)")
    ax.set_title(title)
    ax.legend(fontsize=7.5, loc="best")
    ax.grid(True, ls=":", alpha=0.45)
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(f"PlotEQCompareRuns: wrote {out}  ({plotted} series)")


HELP = """\
usage: python plot/PlotEQCompareRuns.py [runDir ...]

  no args   all runs under OSU_SSI_BRIDGE_DATA/opensees data (.lnk ok)
  runDir    one or more dump folders
  output    OSU_SSI_PLOTS/compare/hist_ux_complete.png
            OSU_SSI_PLOTS/compare/hist_ux_all.png
  complete  t_last >= Trec - 1 s (same rule as PlotEQ truncated_end)
"""


def main() -> int:
    if any(a in ("-h", "--help") for a in sys.argv[1:]):
        print(HELP, end="")
        return 0
    if len(sys.argv) > 1:
        runs = [Path(a).resolve() for a in sys.argv[1:]]
    else:
        runs = discover_runs()
    if not runs:
        print("PlotEQCompareRuns: no runs with pier_top_disp found", file=sys.stderr)
        return 1

    complete_runs: list[Path] = []
    incomplete_runs: list[Path] = []
    print(f"PlotEQCompareRuns: {len(runs)} run(s)")
    for r in runs:
        t_last, trec, ok = run_duration(r)
        tag = "complete" if ok else "incomplete"
        print(f"  [{tag}] {r.name}  t_last={t_last}  Trec={trec}")
        (complete_runs if ok else incomplete_runs).append(r)

    out_dir = PLOTS_ROOT / "compare"
    if complete_runs:
        plot_compare(
            complete_runs,
            out_dir / "hist_ux_complete.png",
            "Storm Wave — pier top ux (complete runs only)",
            mark_incomplete=False,
        )
    else:
        print("PlotEQCompareRuns: no complete runs for hist_ux_complete.png")

    plot_compare(
        runs,
        out_dir / "hist_ux_all.png",
        "Storm Wave — pier top ux (all runs; dashed = incomplete)",
        mark_incomplete=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
