#!/usr/bin/env python3
"""
Goals
-----
Overlay meaSigOS (actuator, native lab $t$) and OpenSees pier UX (mapped with
atTarget handshake) on shared lab axes (full | D5–95), with amber lines at each
mid-run typeConv3==2 onset.

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
from PlotEQCompareRuns import apply_paper_style, load_mat_mea_feedback, model_disp_to_proto_mm, pier_ux_legend_label
from PlotEQComparePairs import (
    COLOR_OTHER,
    COLOR_REF,
    load_pier_ux_mm,
)
from gm_duration import d595_proto_window, gm_start_time_s
from lab_paths import (
    CYLINDER_LENGTH_SCALE,
    LOCAL_OPENSEES_DATA,
    YLIM_DISP_PROTO_MM,
    full_xlim_model_s,
    load_lab_runs_rows,
    resolve_opensees_data,
    test_os_plots_dir,
)
from lab_time_map import (
    LABEL_T_LAB,
    add_dual_time_xaxis_lab,
    map_os_to_lab,
    os_window_to_lab,
    slowdown_lab_onsets,
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


def mark_slowdowns(ax, t_lab: list[float]) -> int:
    """Fixed-thickness amber vertical line at each slowdown onset (lab s)."""
    for t in t_lab:
        ax.axvline(
            t,
            color=COLOR_SLOW,
            alpha=ALPHA_SLOW,
            lw=LW_SLOW,
            solid_capstyle="butt",
            zorder=1,
        )
    return len(t_lab)


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

    mea = load_mat_mea_feedback(mat)
    dump_path = root / dump
    pier = load_pier_ux_mm(dump_path)
    if mea is None or pier is None:
        print(f"PlotActuatorVsPier: skip {test_id} (missing mea or pier)", file=sys.stderr)
        return 1
    t_act, u_m = mea
    u_act = model_disp_to_proto_mm(u_m)
    t_pier_os, u_pier = pier
    try:
        mapped = map_os_to_lab(t_pier_os, mat, dump_dir=dump_path)
    except RuntimeError as exc:
        print(f"PlotActuatorVsPier: skip {test_id} ({exc})", file=sys.stderr)
        return 1
    t_pier = mapped.t_lab
    pier_label = pier_ux_legend_label(dump_path)
    t_slow = slowdown_lab_onsets(mat)

    t0 = gm_start_time_s(dump_path)
    d595_os = d595_proto_window(gm_start_s=t0)
    d595 = (
        os_window_to_lab(d595_os[0], d595_os[1], mat, dump_dir=dump_path)
        if d595_os is not None
        else None
    )

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

    ax_f.set_xlim(*full_xlim_model_s(t_act if t_act.size >= t_pier.size else t_pier))
    ax_f.set_ylim(*YLIM_DISP_PROTO_MM)
    if d595 is not None:
        ax_z.set_xlim(d595[0], d595[1])

    ax_f.set_xlabel(LABEL_T_LAB)
    ax_z.set_xlabel(LABEL_T_LAB)
    ax_z.tick_params(labelleft=False)
    add_dual_time_xaxis_lab(ax_f, top=True)
    add_dual_time_xaxis_lab(ax_z, top=True)

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
