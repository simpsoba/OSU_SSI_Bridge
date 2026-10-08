#!/usr/bin/env python3
"""
Goals
-----
Actuator hysteresis on lab-mapped OpenSees force (atTarget map).
F from daqForce (−); Δu / udot from meaSigOS on lab $t$ (udot = ∇ meaSig).

  python plot/PlotActuatorHyst.py
  python plot/PlotActuatorHyst.py F06 W05

Writes ``plots/runs/<Test>/os/hyst_actuator_Fu.png`` and
``hyst_actuator_Fu_by_udot.png``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from PlotActuatorForce import (
    COLOR_FRC,
    COLOR_WAVE,
    DEFAULT_FONT_SCALE,
    FORCE_SCALE_FROUDE,
    HYDRO_N_PERIODS,
    TIME_SCALE_FROUDE,
    detect_wave_hit_proto_s,
    hydro_force_zoom_proto,
    load_daq_force_kn,
    mat_dump_for_test,
    run_title,
    scale_paper_fonts,
    wave_period_proto_s,
)
from PlotEQCompareRuns import apply_paper_style, load_mat_mea_feedback, model_disp_to_proto_mm
from PlotHydroSpectra import default_tests
from gm_duration import d595_proto_window, gm_start_time_s
from lab_paths import (
    CYLINDER_LENGTH_SCALE,
    LOCAL_OPENSEES_DATA,
    M_TO_MM,
    XLIM_FULL_MODEL_S,
    resolve_opensees_data,
    test_os_plots_dir,
)
from lab_time_map import LABEL_T_LAB, map_os_to_lab

PAD = 0.08
WAVE_HALF_MODEL_S = 8.0
COLOR_ZONE = COLOR_WAVE
ALPHA_ZONE = 0.18
COLOR_RIGHT = "#1565c0"
COLOR_LEFT = "#c62828"
LW = 1.05
UDOT = r"\dot{(\Delta u)}"


def padded_lim(x: np.ndarray, pad: float = PAD) -> tuple[float, float] | None:
    if x.size == 0:
        return None
    lo, hi = float(np.nanmin(x)), float(np.nanmax(x))
    if not np.isfinite(lo) or not np.isfinite(hi):
        return None
    span = hi - lo
    if span <= 0.0:
        span = max(abs(hi), 1.0) * 0.05
    return lo - pad * span, hi + pad * span


def wave_window_lab(
    t_os: np.ndarray,
    f_os: np.ndarray,
    t_lab_os: np.ndarray,
    *,
    test_id: str,
    dump_path: Path,
) -> tuple[tuple[float, float], str] | None:
    t0 = gm_start_time_s(dump_path)
    d595 = d595_proto_window(t0)
    t_wave_period = wave_period_proto_s(test_id)
    wave_hit = detect_wave_hit_proto_s(t_os, f_os)
    t_wave_peak = None if wave_hit is None else wave_hit[1]
    if t_wave_period is not None:
        hydro = hydro_force_zoom_proto(
            t_os, f_os, t_wave_proto=t_wave_period, d595=d595, t_wave_peak=t_wave_peak
        )
        if hydro is None:
            return None
        m = (t_os >= hydro[0]) & (t_os <= hydro[1])
        if not np.any(m):
            return None
        return (
            (float(np.min(t_lab_os[m])), float(np.max(t_lab_os[m]))),
            rf"wave (${HYDRO_N_PERIODS:g}\,T_{{\mathrm{{w}}}}$)",
        )
    if wave_hit is None:
        return None
    t_on, t_pk = wave_hit
    half = WAVE_HALF_MODEL_S * TIME_SCALE_FROUDE
    m = (t_os >= t_on - 0.25 * half) & (t_os <= t_pk + half)
    if not np.any(m):
        return None
    return (
        (float(np.min(t_lab_os[m])), float(np.max(t_lab_os[m]))),
        r"wave $|F|$ peak",
    )


def colored_line(ax, x, y, c):
    pts = np.column_stack([x, y]).reshape(-1, 1, 2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    cols = np.where(c[:-1] >= 0.0, COLOR_RIGHT, COLOR_LEFT)
    ax.add_collection(LineCollection(segs, colors=cols, linewidths=LW, rasterized=True))
    ax.autoscale_view()


def write_fu(test_id: str) -> int:
    pair = mat_dump_for_test(test_id)
    if pair is None:
        return 1
    mat, dump = pair
    root = resolve_opensees_data() or LOCAL_OPENSEES_DATA
    dump_path = root / dump
    mea = load_mat_mea_feedback(mat)
    frc = load_daq_force_kn(dump_path)
    if mea is None or frc is None:
        return 1
    t_lab, u_m = mea
    u_lab = model_disp_to_proto_mm(u_m)
    t_os, f_kn = frc
    try:
        mapped = map_os_to_lab(t_os, mat, dump_dir=dump_path)
    except RuntimeError as exc:
        print(f"PlotActuatorHyst: skip {test_id} ({exc})", file=sys.stderr)
        return 1
    t_lab_os, f_end = mapped.t_lab, mapped.f_end
    # force onto dense lab grid for histories; hysteresis uses OS samples on mapped t
    f_on_lab = np.interp(t_lab, t_lab_os, f_kn, left=np.nan, right=np.nan)
    u_at_os = np.interp(t_lab_os, t_lab, u_lab)

    win = wave_window_lab(t_os, f_kn, t_lab_os, test_id=test_id, dump_path=dump_path)
    if win is None:
        print(f"skip {test_id}: no wave window", file=sys.stderr)
        return 1
    hydro_lab, wave_title = win
    m_w = (t_lab_os >= hydro_lab[0]) & (t_lab_os <= hydro_lab[1])

    apply_paper_style()
    scale_paper_fonts(DEFAULT_FONT_SCALE)
    lam = CYLINDER_LENGTH_SCALE
    fig = plt.figure(figsize=(10.8, 8.6), layout="constrained")
    gs = GridSpec(3, 1, figure=fig, height_ratios=[1.0, 1.0, 1.35])
    ax_f = fig.add_subplot(gs[0])
    ax_u = fig.add_subplot(gs[1], sharex=ax_f)
    ax_h = fig.add_subplot(gs[2])

    step = max(1, len(t_lab) // 12000)
    # top: histories (lab grid may be long — light decimate OK)
    ax_f.plot(t_lab[::step], f_on_lab[::step], color=COLOR_FRC, lw=1.05)
    ax_u.plot(t_lab[::step], u_lab[::step], color=COLOR_FRC, lw=1.05)
    for ax in (ax_f, ax_u):
        ax.axvspan(hydro_lab[0], hydro_lab[1], color=COLOR_ZONE, alpha=ALPHA_ZONE, zorder=0, lw=0)
        ax.grid(True, ls=":", alpha=0.45)
    ax_f.set_ylabel(r"$-F_{\mathrm{daq}}$ (kN) proto")
    ax_u.set_ylabel(r"$\Delta u$ (mm) proto")
    ax_u.set_xlabel(LABEL_T_LAB)
    ax_f.set_title(rf"$F(t)$, $\Delta u(t)$ on lab clock  ·  {wave_title} band", loc="left")
    # model dual
    ax_f.secondary_yaxis(
        "right",
        functions=(
            lambda f_p: f_p / FORCE_SCALE_FROUDE,
            lambda f_m: f_m * FORCE_SCALE_FROUDE,
        ),
    ).set_ylabel(r"$-F/\lambda^{3}$ (kN) model")
    ax_u.secondary_yaxis(
        "right",
        functions=(lambda u: u / lam, lambda u: u * lam),
    ).set_ylabel(r"$\Delta u/\lambda$ (mm) model")

    # hysteresis: grey full (every mapped OS sample), navy = wave band
    ax_h.plot(
        u_at_os,
        f_kn,
        color="#bdbdbd",
        lw=0.7,
        alpha=0.55,
        zorder=1,
        label="full (mapped)",
    )
    ax_h.plot(
        u_at_os[m_w],
        f_kn[m_w],
        color=COLOR_FRC,
        lw=1.15,
        zorder=3,
        label=wave_title,
    )
    ax_h.set_xlabel(r"$\Delta u$ (mm) prototype")
    ax_h.set_ylabel(r"$-F_{\mathrm{daq}}$ (kN) prototype")
    ax_h.grid(True, ls=":", alpha=0.45)
    ax_h.legend(loc="best", frameon=True, fancybox=False, edgecolor="#333")
    lim_u = padded_lim(u_at_os[m_w])
    lim_f = padded_lim(f_kn[m_w])
    if lim_u:
        ax_h.set_xlim(*lim_u)
    if lim_f:
        ax_h.set_ylim(*lim_f)
    ax_h.secondary_xaxis(
        "top", functions=(lambda x: x / lam, lambda x: x * lam)
    ).set_xlabel(r"$\Delta u/\lambda$ (mm) model")
    ax_h.secondary_yaxis(
        "right",
        functions=(
            lambda f_p: f_p / FORCE_SCALE_FROUDE,
            lambda f_m: f_m * FORCE_SCALE_FROUDE,
        ),
    ).set_ylabel(r"$-F/\lambda^{3}$ (kN) model")

    # xlim histories: full model span
    ax_f.set_xlim(*XLIM_FULL_MODEL_S)

    fig.suptitle(
        rf"{run_title(test_id)}  ·  atTarget map  ·  "
        rf"$f_{{\mathrm{{end}}}}={f_end*1e3:.1f}\,\mathrm{{ms}}$",
        y=1.01,
        fontsize=plt.rcParams["axes.labelsize"],
    )
    out = test_os_plots_dir(test_id) / "hyst_actuator_Fu.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"PlotActuatorHyst: wrote {out}")
    return 0


def write_udot(test_id: str) -> int:
    pair = mat_dump_for_test(test_id)
    if pair is None:
        return 1
    mat, dump = pair
    root = resolve_opensees_data() or LOCAL_OPENSEES_DATA
    dump_path = root / dump
    mea = load_mat_mea_feedback(mat)
    frc = load_daq_force_kn(dump_path)
    if mea is None or frc is None:
        return 1
    t_lab, u_m = mea
    u_lab = model_disp_to_proto_mm(u_m)
    t_os, f_kn = frc
    try:
        mapped = map_os_to_lab(t_os, mat, dump_dir=dump_path)
    except RuntimeError as exc:
        print(f"PlotActuatorHyst: skip {test_id} ({exc})", file=sys.stderr)
        return 1
    t_lab_os, f_end = mapped.t_lab, mapped.f_end
    f_on_lab = np.interp(t_lab, t_lab_os, f_kn, left=np.nan, right=np.nan)
    u_at_os = np.interp(t_lab_os, t_lab, u_lab)
    u_lab_mm0 = (u_m - u_m[0]) * M_TO_MM
    udot_proto = np.gradient(u_lab_mm0, t_lab) * TIME_SCALE_FROUDE
    udot_at_os = np.interp(t_lab_os, t_lab, udot_proto)

    win = wave_window_lab(t_os, f_kn, t_lab_os, test_id=test_id, dump_path=dump_path)
    if win is None:
        return 1
    hydro_lab, wave_title = win
    m = (t_lab_os >= hydro_lab[0]) & (t_lab_os <= hydro_lab[1])
    u_w, f_w, v_w = u_at_os[m], f_kn[m], udot_at_os[m]
    mean_r = float(np.nanmean(f_w[v_w > 0])) if np.any(v_w > 0) else float("nan")
    mean_l = float(np.nanmean(f_w[v_w < 0])) if np.any(v_w < 0) else float("nan")

    apply_paper_style()
    scale_paper_fonts(DEFAULT_FONT_SCALE)
    lam = CYLINDER_LENGTH_SCALE
    fig = plt.figure(figsize=(10.8, 11.0), layout="constrained")
    gs = GridSpec(4, 2, figure=fig, height_ratios=[1.0, 1.0, 1.0, 1.35])
    ax_f = fig.add_subplot(gs[0, :])
    ax_u = fig.add_subplot(gs[1, :], sharex=ax_f)
    ax_v = fig.add_subplot(gs[2, :], sharex=ax_f)
    ax_fu = fig.add_subplot(gs[3, 0])
    ax_fv = fig.add_subplot(gs[3, 1], sharey=ax_fu)

    step = max(1, len(t_lab) // 12000)
    for ax, y, ylab in (
        (ax_f, f_on_lab, r"$-F_{\mathrm{daq}}$ (kN) proto"),
        (ax_u, u_lab, r"$\Delta u$ (mm) proto"),
        (ax_v, udot_proto, rf"${UDOT}$ (mm/s) proto"),
    ):
        ax.plot(t_lab[::step], y[::step], color=COLOR_FRC, lw=1.05)
        ax.axvspan(hydro_lab[0], hydro_lab[1], color=COLOR_ZONE, alpha=ALPHA_ZONE, zorder=0, lw=0)
        ax.set_ylabel(ylab)
        ax.grid(True, ls=":", alpha=0.45)
        ax.set_xlim(*XLIM_FULL_MODEL_S)
    ax_v.axhline(0, color="#9e9e9e", lw=0.7)
    ax_v.set_xlabel(LABEL_T_LAB)
    for ax in (ax_f, ax_u):
        ax.tick_params(labelbottom=False)

    ax_fu.plot(u_at_os, f_kn, color="#bdbdbd", lw=0.65, alpha=0.5, zorder=1)
    colored_line(ax_fu, u_w, f_w, v_w)
    colored_line(ax_fv, v_w, f_w, v_w)
    ax_fu.axhline(0, color="#9e9e9e", lw=0.7)
    ax_fv.axhline(0, color="#9e9e9e", lw=0.7)
    ax_fv.axvline(0, color="#9e9e9e", lw=0.7)
    ax_fu.set_xlabel(r"$\Delta u$ (mm) prototype")
    ax_fv.set_xlabel(rf"${UDOT}$ (mm/s) prototype")
    ax_fu.set_ylabel(r"$-F_{\mathrm{daq}}$ (kN) prototype")
    for ax in (ax_fu, ax_fv):
        ax.grid(True, ls=":", alpha=0.45)
    lu, lf, lv = padded_lim(u_w), padded_lim(f_w), padded_lim(v_w)
    if lu:
        ax_fu.set_xlim(*lu)
    if lf:
        ax_fu.set_ylim(*lf)
        ax_fv.set_ylim(*lf)
    if lv:
        ax_fv.set_xlim(*lv)
    ax_fu.legend(
        handles=[
            Line2D([0], [0], color=COLOR_RIGHT, lw=2, label=rf"${UDOT}>0$"),
            Line2D([0], [0], color=COLOR_LEFT, lw=2, label=rf"${UDOT}<0$"),
            Patch(facecolor=COLOR_ZONE, alpha=ALPHA_ZONE, label=wave_title),
        ],
        loc="best",
        fontsize=8,
        frameon=True,
        fancybox=False,
        edgecolor="#333",
    )
    fig.suptitle(
        rf"{run_title(test_id)}  ·  atTarget map  ·  "
        rf"$f_{{\mathrm{{end}}}}={f_end*1e3:.1f}\,\mathrm{{ms}}$  ·  "
        rf"$\langle -F\rangle_{{\mathrm{{R}}}}={mean_r:.2f}$, "
        rf"$\langle -F\rangle_{{\mathrm{{L}}}}={mean_l:.2f}$ kN",
        y=1.01,
        fontsize=plt.rcParams["axes.labelsize"],
    )
    out = test_os_plots_dir(test_id) / "hyst_actuator_Fu_by_udot.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"PlotActuatorHyst: wrote {out}")
    return 0


def main() -> int:
    if any(a in ("-h", "--help") for a in sys.argv[1:]):
        print(__doc__)
        return 0
    tests = [a for a in sys.argv[1:] if not a.startswith("-")] or default_tests()
    if not tests:
        print("PlotActuatorHyst: no tests", file=sys.stderr)
        return 1
    rc = 0
    for tid in tests:
        rc = max(rc, write_fu(tid), write_udot(tid))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
