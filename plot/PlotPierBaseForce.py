#!/usr/bin/env python3
"""
Goals
-----
Plot numerical pier-base hinge section force (P, Mz) vs time (full | D5–95)
with amber lines at each typeConv3==2 onset. Dual axes: prototype and model
via Froude (P/λ³, M/λ⁴; t/√λ).

  python plot/PlotPierBaseForce.py --mesh-ladder
  python plot/PlotPierBaseForce.py F06 F14

Source: dump ``pier_hinge_force.out`` (OpenSees t; N, N·m; prototype).
Bottom zeroLength hinge for lumpedPlasticity (section 1 for forceBeamColumn).
Several Test IDs share one P axis and one M axis (max |·| over full records).
Writes ``plots/runs/<Test>/os/hist_pier_base_PM.png``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from PlotActuatorForce import (
    COLOR_WAVE,
    DEFAULT_FONT_SCALE,
    FONT_SCALE_KEYS,
    LW_WAVE,
    MESH_LADDER_TESTS,
    detect_wave_onset_proto_s,
    load_daq_force_kn,
    mark_slowdowns,
    mark_wave,
    mat_dump_for_test,
    run_title,
    scale_paper_fonts,
    slowdown_times_proto_s,
)
from PlotEQ import loadtxt_partial, subplots_full_zoom
from PlotEQComparePairs import COLOR_OTHER, LABEL_T_PROTO, add_dual_time_xaxis
from PlotEQCompareRuns import apply_paper_style
from gm_duration import d595_proto_window, gm_start_time_s
from lab_paths import (
    CYLINDER_LENGTH_SCALE,
    LOCAL_OPENSEES_DATA,
    TIME_SCALE_FROUDE,
    XLIM_FULL_PROTO_S,
    full_xlim_proto_s,
    resolve_opensees_data,
    test_os_plots_dir,
)

# Froude: F ~ λ³, M ~ F·L ~ λ⁴ (same density).
FORCE_SCALE_FROUDE = CYLINDER_LENGTH_SCALE**3
MOMENT_SCALE_FROUDE = CYLINDER_LENGTH_SCALE**4
N_TO_KN = 1.0e-3
NM_TO_KNM = 1.0e-3

COLOR_SIG = COLOR_OTHER
COLOR_SLOW = "#FFC04D"
ALPHA_SLOW = 0.32
LW_SLOW = 1.0
LW_SIG = 1.15
OUT_NAME = "hist_pier_base_PM.png"
HINGE_FORCE_NAME = "pier_hinge_force.out"


def d595_window(gm_start_s: float = 0.0) -> tuple[float, float] | None:
    """GM D5–95 in prototype seconds (incl. gmStartTime), or None if VT2 missing."""
    return d595_proto_window(gm_start_s=gm_start_s)


def resolve_hinge_force_path(dump_path: Path) -> Path | None:
    """Serial file or first parallel shard of pier_hinge_force.out."""
    plain = dump_path / HINGE_FORCE_NAME
    if plain.is_file():
        return plain
    shards = sorted(dump_path.glob(HINGE_FORCE_NAME + ".*"))
    return shards[0] if shards else None


def load_pier_base_pm(
    dump_path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """
    Bottom-hinge section force history.

    Args:    dump_path  LOCAL opensees_data/<DumpFolder>
    Returns: (t_proto_s, P_kN, M_kNm) or None
    """
    path = resolve_hinge_force_path(dump_path)
    if path is None:
        return None
    data = loadtxt_partial(path)
    if data.size == 0 or data.ndim != 2 or data.shape[1] < 3:
        return None
    t = np.asarray(data[:, 0], dtype=float)
    p_kn = np.asarray(data[:, 1], dtype=float) * N_TO_KN
    m_knm = np.asarray(data[:, 2], dtype=float) * NM_TO_KNM
    return t, p_kn, m_knm


def peak_abs(y: np.ndarray) -> float:
    """Peak |y| over finite samples."""
    if y.size == 0:
        return 0.0
    v = float(np.nanmax(np.abs(y)))
    return v if np.isfinite(v) else 0.0


def shared_ylims(
    test_ids: list[str],
    *,
    data_root: Path | None = None,
) -> tuple[tuple[float, float], tuple[float, float]] | None:
    """
    Symmetric (ylim_P, ylim_M) from max |P| and |M| among listed tests.

    Returns: ((-p, p), (-m, m)) in kN / kn·m prototype, or None
    """
    root = data_root or resolve_opensees_data() or LOCAL_OPENSEES_DATA
    p_max = 0.0
    m_max = 0.0
    for tid in test_ids:
        pair = mat_dump_for_test(tid)
        if pair is None:
            continue
        pm = load_pier_base_pm(root / pair[1])
        if pm is None:
            continue
        _, p_kn, m_knm = pm
        p_max = max(p_max, peak_abs(p_kn))
        m_max = max(m_max, peak_abs(m_knm))
    if p_max <= 0.0 and m_max <= 0.0:
        return None
    return (
        (-1.08 * max(p_max, 1.0e-9), 1.08 * max(p_max, 1.0e-9)),
        (-1.08 * max(m_max, 1.0e-9), 1.08 * max(m_max, 1.0e-9)),
    )


def add_model_time_top(ax) -> None:
    """Secondary x-axis: model time = t_proto / √λ."""
    ax.secondary_xaxis(
        "top",
        functions=(
            lambda t_proto: t_proto / TIME_SCALE_FROUDE,
            lambda t_model: t_model * TIME_SCALE_FROUDE,
        ),
    )


def write_plot(
    test_id: str,
    *,
    data_root: Path | None = None,
    font_scale: float = DEFAULT_FONT_SCALE,
    ylim_p: tuple[float, float] | None = None,
    ylim_m: tuple[float, float] | None = None,
) -> int:
    """
    Write one pier-base P/M PNG for a Test ID.

    Returns: 0 ok, 1 skip/error
    """
    apply_paper_style()
    scale_paper_fonts(font_scale)
    pair = mat_dump_for_test(test_id)
    if pair is None:
        print(f"PlotPierBaseForce: skip {test_id} (no mat+dump)", file=sys.stderr)
        return 1
    mat, dump = pair
    root = data_root or resolve_opensees_data() or LOCAL_OPENSEES_DATA
    out = test_os_plots_dir(test_id) / OUT_NAME

    pm = load_pier_base_pm(root / dump)
    if pm is None:
        print(f"PlotPierBaseForce: skip {test_id} (no {HINGE_FORCE_NAME})", file=sys.stderr)
        return 1
    t, p_kn, m_knm = pm
    t_slow = slowdown_times_proto_s(mat)
    # Wave onset from actuator force (same clock / dump as hinge).
    frc = load_daq_force_kn(root / dump)
    t_wave = (
        detect_wave_onset_proto_s(frc[0], frc[1]) if frc is not None else None
    )
    t0 = gm_start_time_s(root / dump)
    d595 = d595_window(t0)

    fig_h = 6.6 * (0.65 + 0.35 * font_scale)
    fig, axes_f, axes_z = subplots_full_zoom(2, fig_h=fig_h, sharey="row", wspace=0.04)
    series = (
        (axes_f[0], axes_z[0], p_kn, ylim_p, r"$P$ (kN) prototype", r"$P/\lambda^{3}$ (kN) model", FORCE_SCALE_FROUDE, "axial $P$"),
        (axes_f[1], axes_z[1], m_knm, ylim_m, r"$M$ (kN·m) prototype", r"$M/\lambda^{4}$ (kN·m) model", MOMENT_SCALE_FROUDE, "moment $M_z$"),
    )

    n_slow = 0
    has_wave = False
    for i, (ax_f, ax_z, y, ylim, ylab, ylab_m, scale, label) in enumerate(series):
        n_slow = mark_slowdowns(ax_f, t_slow)
        mark_slowdowns(ax_z, t_slow)
        has_wave = mark_wave(ax_f, t_wave) or has_wave
        mark_wave(ax_z, t_wave)
        for ax in (ax_f, ax_z):
            ax.plot(t, y, color=COLOR_SIG, lw=LW_SIG, label=label, zorder=5)
            ax.grid(True, ls=":", alpha=0.45)
            ax.axhline(0.0, color="#666666", lw=0.6, zorder=0)
        ax_f.set_xlim(*full_xlim_proto_s(t))
        if d595 is not None:
            ax_z.set_xlim(d595[0], d595[1])
        if ylim is not None:
            ax_f.set_ylim(*ylim)
        else:
            y_abs = peak_abs(y)
            if y_abs > 0.0:
                ax_f.set_ylim(-1.08 * y_abs, 1.08 * y_abs)
        ax_z.tick_params(labelleft=False)
        ax_f.set_ylabel(ylab)
        sec_y = ax_z.secondary_yaxis(
            "right",
            functions=(
                lambda v, s=scale: v / s,
                lambda v, s=scale: v * s,
            ),
        )
        sec_y.set_ylabel(ylab_m)
        if i == 0:
            add_dual_time_xaxis(ax_f, top=True)
            add_model_time_top(ax_z)
            ax_f.tick_params(labelbottom=False)
            ax_z.tick_params(labelbottom=False)
        else:
            ax_f.set_xlabel(LABEL_T_PROTO)
            ax_z.set_xlabel(LABEL_T_PROTO)
            add_model_time_top(ax_f)
            add_model_time_top(ax_z)

    engine = fig.get_layout_engine()
    if engine is not None:
        engine.set(w_pad=0.015, h_pad=0.015, wspace=0.02, hspace=0.06)

    handles, labels = axes_f[0].get_legend_handles_labels()
    handles2, labels2 = axes_f[1].get_legend_handles_labels()
    handles = handles + handles2
    labels = labels + labels2
    if n_slow:
        handles.append(
            Line2D(
                [0],
                [0],
                color=COLOR_SLOW,
                alpha=ALPHA_SLOW,
                lw=LW_SLOW,
                label=f"slowdown ({n_slow})",
            )
        )
        labels.append(f"slowdown ({n_slow})")
    if has_wave and t_wave is not None:
        handles.append(
            Line2D(
                [0],
                [0],
                color=COLOR_WAVE,
                lw=LW_WAVE,
                ls="--",
                label=rf"wave ($t/\sqrt{{\lambda}}={t_wave / TIME_SCALE_FROUDE:.0f}$ s)",
            )
        )
        labels.append(
            rf"wave ($t/\sqrt{{\lambda}}={t_wave / TIME_SCALE_FROUDE:.0f}$ s)"
        )
    axes_f[0].legend(
        handles,
        labels,
        loc="upper right",
        fontsize=plt.rcParams["legend.fontsize"],
        frameon=True,
        fancybox=False,
        edgecolor="#333333",
        facecolor="white",
        framealpha=1.0,
        handlelength=1.8,
    )

    fig.suptitle(
        run_title(test_id) + r"  —  pier base hinge",
        y=1.01,
        fontsize=plt.rcParams["axes.labelsize"],
    )

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    wave_txt = (
        f"wave={t_wave:.1f}s proto ({t_wave / TIME_SCALE_FROUDE:.1f}s model)"
        if t_wave is not None
        else "wave=none"
    )
    print(f"PlotPierBaseForce: wrote {out}  (lines={n_slow}, {wave_txt})")
    return 0


def parse_argv(argv: list[str]) -> tuple[list[str], float]:
    """Split Test IDs from ``--font-scale`` / ``--mesh-ladder``."""
    tests: list[str] = []
    font_scale = DEFAULT_FONT_SCALE
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in ("-h", "--help"):
            print(__doc__)
            raise SystemExit(0)
        if a == "--mesh-ladder":
            tests.extend(MESH_LADDER_TESTS)
        elif a == "--font-scale":
            i += 1
            if i >= len(argv):
                raise SystemExit("PlotPierBaseForce: --font-scale needs a number")
            font_scale = float(argv[i])
        elif a.startswith("--font-scale="):
            font_scale = float(a.split("=", 1)[1])
        elif a.startswith("-"):
            raise SystemExit(f"PlotPierBaseForce: unknown option {a}")
        else:
            tests.append(a)
        i += 1
    return tests, font_scale


def main() -> int:
    tests, font_scale = parse_argv(sys.argv[1:])
    if not tests:
        tests = list(MESH_LADDER_TESTS)
    ylims = shared_ylims(tests) if len(tests) > 1 else None
    ylim_p = ylim_m = None
    if ylims is not None:
        ylim_p, ylim_m = ylims
        print(
            f"PlotPierBaseForce: shared P [{ylim_p[0]:.3g}, {ylim_p[1]:.3g}] kN  "
            f"M [{ylim_m[0]:.3g}, {ylim_m[1]:.3g}] kN·m"
        )
    rc = 0
    for tid in tests:
        rc = max(
            rc,
            write_plot(tid, font_scale=font_scale, ylim_p=ylim_p, ylim_m=ylim_m),
        )
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
