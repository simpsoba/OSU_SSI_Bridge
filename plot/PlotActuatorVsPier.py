#!/usr/bin/env python3
"""
Goals
-----
Overlay meaSigOS (actuator) and OpenSees pier UX on shared prototype axes
(full | D5–95), with amber vertical lines at each typeConv3==2 onset.

twoNodeLink runs use relative pier UX (inner top node minus inner base node,
nodes 4--2 for lumpedPlasticity) so the numerical line matches the actuator DOF.

  python plot/PlotActuatorVsPier.py
  python plot/PlotActuatorVsPier.py F07 F14

Writes ``plots/runs/<Test>/os/hist_ux_actuator_vs_pier.png``.
Type is 75% larger than the compare-plot paper style (``--font-scale 1`` to undo).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D

from PlotEQ import subplots_full_zoom
from PlotEQCompareRuns import pier_ux_legend_label
from PlotEQComparePairs import (
    COLOR_OTHER,
    COLOR_REF,
    LABEL_T_PROTO,
    SLOWDOWN_STATE,
    add_dual_time_xaxis,
    load_mea_ux_proto,
    load_pier_ux_mm,
    load_type_conv3,
)
from PlotEQCompareRuns import apply_paper_style
from gm_duration import arias_significant_duration
from lab_paths import (
    CYLINDER_LENGTH_SCALE,
    LOCAL_OPENSEES_DATA,
    TIME_SCALE_FROUDE,
    XLIM_FULL_PROTO_S,
    YLIM_DISP_PROTO_MM,
    load_lab_runs_rows,
    resolve_opensees_data,
    test_os_plots_dir,
)

COLOR_ACT = COLOR_OTHER  # #001F3F
COLOR_PIER = COLOR_REF  # #B0B0B0
COLOR_SLOW = "#FFC04D"
ALPHA_SLOW = 0.32
LW_SLOW = 1.0
LW_PIER = 1.25
PIER_HALO = [
    pe.Stroke(linewidth=LW_PIER + 2.5, foreground="white"),
    pe.Normal(),
]
OUT_NAME = "hist_ux_actuator_vs_pier.png"
DEFAULT_FONT_SCALE = 1.75
FONT_SCALE_KEYS = (
    "font.size",
    "axes.labelsize",
    "xtick.labelsize",
    "ytick.labelsize",
    "legend.fontsize",
)


def scale_paper_fonts(factor: float) -> None:
    """Multiply paper-style type sizes. factor=1.75 is 75% larger."""
    if factor == 1.0:
        return
    for key in FONT_SCALE_KEYS:
        plt.rcParams[key] = float(plt.rcParams[key]) * factor


def mat_dump_for_test(test_id: str) -> tuple[str, str] | None:
    """Return (MatFile, DumpFolder) for Test ID, or None if incomplete."""
    for row in load_lab_runs_rows():
        if (row.get("Test") or "").strip() != test_id:
            continue
        mat = (row.get("MatFile") or "").strip()
        dump = (row.get("DumpFolder") or "").strip()
        if mat and dump:
            return mat, dump
        return None
    return None


def slowdown_times_proto_s(mat_name: str) -> list[float]:
    """Prototype time at the start of each contiguous typeConv3==2 episode."""
    pair = load_type_conv3(mat_name)
    if pair is None:
        return []
    t_lab_s, state = pair
    times: list[float] = []
    in_span = False
    for i in range(len(state)):
        if int(state[i]) == SLOWDOWN_STATE and not in_span:
            in_span = True
            times.append(float(t_lab_s[i]) * TIME_SCALE_FROUDE)
        elif int(state[i]) != SLOWDOWN_STATE and in_span:
            in_span = False
    return times


def mark_slowdowns(ax, t_proto: list[float]) -> int:
    """Fixed-thickness amber vertical line at each slowdown onset."""
    for t in t_proto:
        ax.axvline(
            t,
            color=COLOR_SLOW,
            alpha=ALPHA_SLOW,
            lw=LW_SLOW,
            solid_capstyle="butt",
            zorder=1,
        )
    return len(t_proto)


def write_plot(
    test_id: str, *, data_root: Path | None = None, font_scale: float = DEFAULT_FONT_SCALE
) -> int:
    """
    Write one actuator vs pier PNG for a Test ID.

    Returns: 0 ok, 1 skip/error
    """
    apply_paper_style()
    scale_paper_fonts(font_scale)
    pair = mat_dump_for_test(test_id)
    if pair is None:
        print(f"PlotActuatorVsPier: skip {test_id} (no mat+dump)", file=sys.stderr)
        return 1
    mat, dump = pair
    root = data_root or resolve_opensees_data() or LOCAL_OPENSEES_DATA
    out = test_os_plots_dir(test_id) / OUT_NAME

    mea = load_mea_ux_proto(mat)
    dump_path = root / dump
    pier = load_pier_ux_mm(dump_path)
    if mea is None or pier is None:
        print(f"PlotActuatorVsPier: skip {test_id} (missing mea or pier)", file=sys.stderr)
        return 1
    t_act, u_act = mea
    t_pier, u_pier = pier
    pier_label = pier_ux_legend_label(dump_path)
    t_slow = slowdown_times_proto_s(mat)

    try:
        dur = arias_significant_duration()
        d595 = (float(dur.t5_s), float(dur.t95_s))
    except (OSError, ValueError):
        d595 = None

    fig_h = 4.2 * (0.65 + 0.35 * font_scale)
    fig, axes_f, axes_z = subplots_full_zoom(1, fig_h=fig_h, sharey=True, wspace=0.04)
    ax_f, ax_z = axes_f[0], axes_z[0]

    n_slow = mark_slowdowns(ax_f, t_slow)
    mark_slowdowns(ax_z, t_slow)
    for ax in (ax_f, ax_z):
        ax.plot(
            t_pier,
            u_pier,
            color=COLOR_PIER,
            lw=LW_PIER,
            label=pier_label,
            zorder=2,
            path_effects=PIER_HALO,
        )
        ax.plot(
            t_act,
            u_act,
            color=COLOR_ACT,
            lw=1.15,
            label="actuator (measured)",
            zorder=3,
        )
        ax.grid(True, ls=":", alpha=0.45)

    ax_f.set_xlim(*XLIM_FULL_PROTO_S)
    ax_f.set_ylim(*YLIM_DISP_PROTO_MM)
    if d595 is not None:
        ax_z.set_xlim(d595[0], d595[1])

    ax_f.set_xlabel(LABEL_T_PROTO)
    ax_z.set_xlabel(LABEL_T_PROTO)
    ax_z.tick_params(labelleft=False)
    add_dual_time_xaxis(ax_f, top=True)
    ax_z.secondary_xaxis(
        "top",
        functions=(
            lambda t_proto: t_proto / TIME_SCALE_FROUDE,
            lambda t_model: t_model * TIME_SCALE_FROUDE,
        ),
    )

    engine = fig.get_layout_engine()
    if engine is not None:
        engine.set(w_pad=0.015, h_pad=0.015, wspace=0.02, hspace=0.02)

    ax_f.set_ylabel(
        r"$\Delta u$ (mm) prototype scale"
        "\n"
        r"actuator / pier"
    )
    sec_y = ax_z.secondary_yaxis(
        "right",
        functions=(
            lambda u_proto: u_proto / CYLINDER_LENGTH_SCALE,
            lambda u_model: u_model * CYLINDER_LENGTH_SCALE,
        ),
    )
    sec_y.set_ylabel(r"$\Delta u/\lambda$ (mm) model scale")

    handles, labels = ax_f.get_legend_handles_labels()
    if n_slow:
        handles.append(
            Line2D(
                [0],
                [0],
                color=COLOR_SLOW,
                alpha=ALPHA_SLOW,
                lw=LW_SLOW,
                label="slowdown",
            )
        )
        labels.append("slowdown")
    ax_f.legend(
        handles,
        labels,
        loc="lower right",
        fontsize=plt.rcParams["legend.fontsize"],
        frameon=True,
        fancybox=False,
        edgecolor="#333333",
        facecolor="white",
        framealpha=1.0,
        handlelength=1.8,
    )

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"PlotActuatorVsPier: wrote {out}  (lines={n_slow})")
    return 0


def parse_argv(argv: list[str]) -> tuple[list[str], float]:
    """
    Split Test IDs from ``--font-scale``.

    Returns: (test_ids, font_scale)
    """
    tests: list[str] = []
    font_scale = DEFAULT_FONT_SCALE
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in ("-h", "--help"):
            print(__doc__)
            raise SystemExit(0)
        if a == "--font-scale":
            i += 1
            if i >= len(argv):
                raise SystemExit("PlotActuatorVsPier: --font-scale needs a number")
            font_scale = float(argv[i])
        elif a.startswith("--font-scale="):
            font_scale = float(a.split("=", 1)[1])
        elif a.startswith("-"):
            raise SystemExit(f"PlotActuatorVsPier: unknown option {a}")
        else:
            tests.append(a)
        i += 1
    return tests, font_scale


def main() -> int:
    tests, font_scale = parse_argv(sys.argv[1:])
    if not tests:
        for row in load_lab_runs_rows():
            tid = (row.get("Test") or "").strip()
            mat = (row.get("MatFile") or "").strip()
            dump = (row.get("DumpFolder") or "").strip()
            if tid and mat and dump:
                tests.append(tid)
    rc = 0
    for tid in tests:
        rc = max(rc, write_plot(tid, font_scale=font_scale))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
