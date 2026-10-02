#!/usr/bin/env python3
"""
Goals
-----
Paper-style ground-motion figure for Tohoku FKSH19.NS1:
  (a) stacked ugddot, ugdot, ug time histories
  (b) elastic Sa(T) spectrum (5% damping)

Layout matches typical eq. engineering figures (time stack | spectrum).
Uses the cached spectrum CSV from PlotResponseSpectrum.py when present.

  python plot/PlotGroundMotionFigure.py
  python plot/PlotGroundMotionFigure.py path/to/file.VT2 [out_stem]
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

# Reuse spectrum helpers / font / paths.
import PlotResponseSpectrum as prs
OUT_DIR = HERE / "out" / "spectrum"
DEFAULT_VT2 = prs.DEFAULT_VT2
G = prs.G

FIG_W = 6.0  # in (= textwidth; save without tight bbox)
FIG_H = 3.0  # in — room for stacked TH + panel tags inside the canvas
FONT_SIZE = 9
LINE_COLOR = "black"
LINE_LW_TH = 0.45  # dense time history
LINE_LW_SA = 1.2
T_SPEC_MAX = 4.0  # s — match typical Sa panel range


def displacement_from_vel(vel: np.ndarray, dt: float) -> np.ndarray:
    """
    Integrate velocity to displacement (cumulative trapezoid).

    Args:    vel (m/s), dt (s)
    Returns: disp (m), not baseline-corrected
    """
    disp = np.empty_like(vel)
    disp[0] = 0.0
    disp[1:] = np.cumsum(0.5 * (vel[1:] + vel[:-1]) * dt)
    return disp


def baseline_correct_disp(disp: np.ndarray, t: np.ndarray, degree: int = 2) -> np.ndarray:
    """
    Remove a low-order polynomial trend from integrated displacement.

    PEER VT2 is processed velocity; raw ∫v still drifts from small residual
    bias. A quadratic (default) fit is enough for a paper time-history panel.

    Args:    disp (m), t (s), degree
    Returns: disp_bc (m)
    """
    coef = np.polyfit(t, disp, degree)
    return disp - np.polyval(coef, t)


def nice_ceil(x: float, pad: float = 1.08) -> float:
    """
    Round upward to a clean axis limit (1–1.5–2–2.5–3–4–5–6–8–10)·10^n.

    Args:    x  data peak (positive), pad  headroom before rounding
    Returns: limit >= pad*x
    """
    if x <= 0.0 or not math.isfinite(x):
        return 1.0
    target = float(x) * pad
    exp = math.floor(math.log10(target))
    mant = target / (10.0**exp)
    for step in (1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0):
        if mant <= step + 1e-12:
            return step * (10.0**exp)
    return 10.0 * (10.0**exp)


def _png_from_pdf(pdf_path: Path, png_path: Path, dpi: int = 300) -> None:
    prs._png_from_pdf(pdf_path, png_path, dpi=dpi)


def plot_gm_figure(
    t: np.ndarray,
    acc_g: np.ndarray,
    vel_cm_s: np.ndarray,
    disp_cm: np.ndarray,
    periods: np.ndarray,
    Sa_g: np.ndarray,
    pga_g: float,
    out_stem: Path,
) -> None:
    """
    (a) time histories | (b) Sa spectrum.

    Args:    t (s), acc_g (g), vel_cm_s (cm/s), disp_cm (cm),
             periods (s), Sa_g (g), pga_g (g), out_stem
    Returns: none (writes pdf + png)
    """
    prs.configure_font()

    fig = plt.figure(figsize=(FIG_W, FIG_H), layout="constrained")
    gs = GridSpec(3, 2, figure=fig, width_ratios=[1.15, 1.0], height_ratios=[1, 1, 1])

    ax_acc = fig.add_subplot(gs[0, 0])
    ax_vel = fig.add_subplot(gs[1, 0], sharex=ax_acc)
    ax_dsp = fig.add_subplot(gs[2, 0], sharex=ax_acc)
    ax_sa = fig.add_subplot(gs[:, 1])

    # --- (a) time histories ---
    ax_acc.plot(t, acc_g, color=LINE_COLOR, lw=LINE_LW_TH)
    # Unit g = gravity (math italic), not gram (\mathrm{g}).
    ax_acc.set_ylabel(r"$\ddot{u}_g(t)$ ($g$)")
    ax_acc.set_xlim(0.0, 300.0)

    ax_vel.plot(t, vel_cm_s, color=LINE_COLOR, lw=LINE_LW_TH)
    ax_vel.set_ylabel(r"$\dot{u}_g(t)$ (cm/s)")

    ax_dsp.plot(t, disp_cm, color=LINE_COLOR, lw=LINE_LW_TH)
    ax_dsp.set_ylabel(r"$u_g(t)$ (cm)")
    # Panel tag in the xlabel so constrained layout budgets the bottom margin.
    ax_dsp.set_xlabel(r"Time, $t$ (s)" + "\n" + r"\textbf{(a)}")
    ax_dsp.set_xlim(0.0, 300.0)
    ax_dsp.set_xticks([0, 50, 100, 150, 200, 250, 300])

    for ax in (ax_acc, ax_vel):
        plt.setp(ax.get_xticklabels(), visible=False)

    # Symmetric y-limits: nice 1–2–5 ceiling so ticks aren't glued to the data.
    for ax, y in (
        (ax_acc, acc_g),
        (ax_vel, vel_cm_s),
        (ax_dsp, disp_cm),
    ):
        lim = nice_ceil(float(np.nanmax(np.abs(y))))
        ax.set_ylim(-lim, lim)
        ax.grid(True, which="major")

    fig.align_ylabels([ax_acc, ax_vel, ax_dsp])

    # --- (b) Sa spectrum (linear, prepend PGA at T=0) ---
    mask = periods <= T_SPEC_MAX
    T_plot = np.concatenate([[0.0], periods[mask]])
    Sa_plot = np.concatenate([[pga_g], Sa_g[mask]])
    ax_sa.plot(T_plot, Sa_plot, color=LINE_COLOR, lw=LINE_LW_SA)
    ax_sa.set_xlabel(r"Period, $T_n$ (s)" + "\n" + r"\textbf{(b)}")
    ax_sa.set_ylabel(r"Spectral acceleration, $S_a$ ($g$)")
    ax_sa.set_xlim(0.0, T_SPEC_MAX)
    sa_ymax = nice_ceil(float(np.nanmax(Sa_plot)))
    ax_sa.set_ylim(0.0, sa_ymax)
    ax_sa.grid(True, which="major")

    # T1 / Sa(T1) — labels inside the axes, clear of the spectrum peak.
    T1 = prs.T1
    sa_t1 = float(np.interp(T1, periods, Sa_g))
    ref_c, ref_lw, ref_ls = "#888888", 0.7, (0, (3, 2))
    ax_sa.axvline(T1, color=ref_c, lw=ref_lw, ls=ref_ls, zorder=2)
    ax_sa.axhline(sa_t1, color=ref_c, lw=ref_lw, ls=ref_ls, zorder=2)

    # Sa(T1): left side, near the bottom.
    ax_sa.text(
        0.08,
        0.04 * sa_ymax,
        rf"$S_a(T_1) = {sa_t1:.3f}\,g$",
        ha="left",
        va="bottom",
        color="#555555",
        fontsize=FONT_SIZE,
        zorder=5,
    )
    # T1: just right of the vertical line, above the horizontal mark.
    ax_sa.text(
        T1 + 0.06,
        sa_t1 + 0.04 * sa_ymax,
        rf"$T_1 = {T1:.2f}$\,s",
        ha="left",
        va="bottom",
        color="#555555",
        fontsize=FONT_SIZE,
        zorder=5,
    )

    # Keep a hair of margin on every side so usetex labels are not clipped
    # at the page edge (exact figsize; no bbox_inches="tight").
    fig.get_layout_engine().set(h_pad=0.04, w_pad=0.04, rect=(0.01, 0.02, 0.99, 0.98))

    out_stem.parent.mkdir(parents=True, exist_ok=True)
    pdf_path = out_stem.with_suffix(".pdf")
    png_path = out_stem.with_suffix(".png")
    svg_path = out_stem.with_suffix(".svg")
    fig.savefig(pdf_path)
    fig.savefig(svg_path)
    plt.close(fig)
    _png_from_pdf(pdf_path, png_path, dpi=300)


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    vt2 = Path(args[0]) if args else DEFAULT_VT2
    if not vt2.is_file():
        print(f"VT2 not found: {vt2}", file=sys.stderr)
        return 1

    if len(args) >= 2:
        out_stem = Path(args[1])
    else:
        out_stem = OUT_DIR / "tohoku_FKSH19_NS1_gm_Sa"

    # Time histories from PEER VT2 (velocity file).
    vel_mps, dt = prs.read_peer_vt2(vt2)
    acc_mps2 = prs.accel_from_vel(vel_mps, dt)
    t = np.arange(len(vel_mps), dtype=float) * dt
    disp_raw = displacement_from_vel(vel_mps, dt)
    disp_m = baseline_correct_disp(disp_raw, t, degree=2)
    acc_g = acc_mps2 / G
    vel_cm_s = vel_mps * 100.0
    disp_cm = disp_m * 100.0
    pga_g = float(np.max(np.abs(acc_g)))

    # Spectrum: prefer CSV cache from PlotResponseSpectrum.py.
    spec_stem = OUT_DIR / "tohoku_FKSH19_NS1_SdSvSa"
    cached = prs.load_spectrum_csv(prs.spectrum_csv_path(spec_stem))
    if cached is not None:
        periods, _Sd, _Sv, Sa_g, pga_csv, method = cached
        if np.isfinite(pga_csv):
            pga_g = pga_csv
        print(f"spectrum cache  {prs.spectrum_csv_path(spec_stem)}")
    else:
        periods = np.logspace(np.log10(prs.T_MIN), np.log10(prs.T_MAX), prs.N_PERIOD)
        _Sd, _Sv, Sa_g, method = prs.compute_spectrum(acc_mps2, dt, periods, prs.ZETA)
        print(f"spectrum method {method}")

    plot_gm_figure(
        t, acc_g, vel_cm_s, disp_cm, periods, Sa_g, pga_g, out_stem
    )

    print(f"file   {vt2}")
    print(f"dt     {dt:g} s   npts={len(t)}   duration={t[-1]:.1f} s")
    print(f"PGA    {pga_g:.4f} g")
    print(f"PGV    {float(np.max(np.abs(vel_cm_s))):.2f} cm/s")
    print(f"PGD    {float(np.max(np.abs(disp_cm))):.2f} cm  (quadratic baseline correction)")
    print(f"       raw int(v) PGD={float(np.max(np.abs(disp_raw))*100):.2f} cm")
    print(f"Sa max {float(np.max(Sa_g)):.3f} g")
    sa_t1 = float(np.interp(prs.T1, periods, Sa_g))
    print(f"T_1    {prs.T1:.2f} s   Sa(T1)={sa_t1:.3f} g")
    print(f"wrote  {out_stem.with_suffix('.png')}")
    print(f"wrote  {out_stem.with_suffix('.pdf')}")
    print(f"wrote  {out_stem.with_suffix('.svg')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
