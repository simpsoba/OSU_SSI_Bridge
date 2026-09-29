#!/usr/bin/env python3
"""
Goals
-----
Prototype-scale power spectra (Welch PSD) for Wednesday / listed lab runs:

  - pier-top UX (mm)
  - pier-base hinge moment Mz (kN·m)
  - OpenFresco actuator daqForce (kN) — mainly hydro when outside EQ

  python plot/PlotHydroSpectra.py W01 W05
  python plot/PlotHydroSpectra.py          # all W## with waveT logged

Marks the flume wave frequency f_w = 1/(T_wave · √λ) when known.
Writes ``plots/runs/<Test>/os/psd_pier_ux_M_F.png``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import welch

from PlotActuatorForce import (
    DEFAULT_FONT_SCALE,
    detect_wave_start_proto_s,
    load_daq_force_kn,
    mat_dump_for_test,
    row_for_test,
    run_title,
    scale_paper_fonts,
    wave_period_proto_s,
)
from PlotEQCompareRuns import apply_paper_style
from PlotPierBaseForce import load_pier_base_pm
from gm_duration import d595_proto_window, gm_start_time_s
from lab_paths import (
    LOCAL_OPENSEES_DATA,
    TIME_SCALE_FROUDE,
    load_lab_runs_rows,
    resolve_opensees_data,
    test_os_plots_dir,
)

OUT_NAME = "psd_pier_ux_M_F.png"
COLOR_UX = "#1565C0"
COLOR_M = "#6A1B9A"
COLOR_F = "#001F3F"
COLOR_FW = "#C62828"


def load_pier_ux_mm(dump_path: Path) -> tuple[np.ndarray, np.ndarray] | None:
    """Pier-top UX (prototype mm) vs OpenSees t (prototype s)."""
    pier = dump_path / "pier_top_disp.out.0"
    if not pier.is_file():
        pier = dump_path / "pier_top_disp.out"
    if not pier.is_file():
        return None
    data = np.loadtxt(pier, ndmin=2)
    if data.size == 0 or data.shape[1] < 2:
        return None
    t = np.asarray(data[:, 0], dtype=float)
    ux_mm = np.asarray(data[:, 1], dtype=float) * 1000.0
    return t, ux_mm


def interpolate_uniform(
    t: np.ndarray,
    y: np.ndarray,
    *,
    t0: float,
    t1: float,
) -> tuple[np.ndarray, np.ndarray, float] | None:
    """
    Uniform samples on [t0, t1] for Welch.

    Returns: (t_u, y_u, dt) or None
    """
    m = (t >= t0) & (t <= t1) & np.isfinite(y)
    if int(np.count_nonzero(m)) < 64:
        return None
    tt = t[m]
    yy = y[m]
    dt = float(np.median(np.diff(tt)))
    if not np.isfinite(dt) or dt <= 0.0:
        return None
    t_u = np.arange(float(tt[0]), float(tt[-1]) + 0.5 * dt, dt)
    if t_u.size < 64:
        return None
    y_u = np.interp(t_u, tt, yy)
    y_u = y_u - float(np.mean(y_u))
    return t_u, y_u, dt


def welch_psd(
    y: np.ndarray,
    dt: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    One-sided Welch PSD vs frequency (Hz, prototype clock).

    Args:    y  zero-mean series; dt  sampling interval (s)
    Returns: (f_Hz, PSD)
    """
    nperseg = min(max(256, y.size // 4), y.size)
    if nperseg < 64:
        nperseg = y.size
    f_hz, pxx = welch(y, fs=1.0 / dt, nperseg=nperseg, detrend=False)
    return f_hz, pxx


def analysis_windows(
    t_end: float,
    gm_start: float,
    d595: tuple[float, float] | None,
    t_wave: float | None,
) -> list[tuple[str, float, float]]:
    """Named PSD windows on the prototype clock."""
    out: list[tuple[str, float, float]] = []
    if t_wave is not None and d595 is not None and t_wave < d595[0] - 5.0:
        out.append(("pre-EQ hydro", float(t_wave), float(d595[0])))
    if d595 is not None:
        out.append(("D5–95", float(d595[0]), min(float(d595[1]), t_end)))
        if t_end > d595[1] + 20.0:
            out.append(("post-EQ + wave", float(d595[1]), t_end))
    if not out:
        out.append(("full record", 0.0, t_end))
    return out


def write_plot(test_id: str, *, font_scale: float = DEFAULT_FONT_SCALE) -> int:
    """Write one PSD PNG for a Test ID."""
    apply_paper_style()
    scale_paper_fonts(font_scale)
    pair = mat_dump_for_test(test_id)
    if pair is None:
        print(f"PlotHydroSpectra: skip {test_id} (no mat+dump)", file=sys.stderr)
        return 1
    _, dump = pair
    root = resolve_opensees_data() or LOCAL_OPENSEES_DATA
    dump_path = root / dump

    frc = load_daq_force_kn(dump_path)
    pier = load_pier_ux_mm(dump_path)
    pm = load_pier_base_pm(dump_path)
    if frc is None or pier is None or pm is None:
        print(f"PlotHydroSpectra: skip {test_id} (missing F / ux / M)", file=sys.stderr)
        return 1

    t_f, f_kn = frc
    t_u, ux_mm = pier
    t_m, _p, m_knm = pm
    t_wave_period = wave_period_proto_s(test_id)
    t_wave = (
        detect_wave_start_proto_s(t_f, f_kn, t_wave_period)
        if t_wave_period is not None
        else None
    )
    gm0 = gm_start_time_s(dump_path)
    d595 = d595_proto_window(gm0)
    t_end = float(
        min(float(t_f[-1]), float(t_u[-1]), float(t_m[-1]))
    )
    windows = analysis_windows(t_end, gm0, d595, t_wave)
    f_wave = (1.0 / t_wave_period) if t_wave_period else None

    series = (
        ("pier top $u_x$", t_u, ux_mm, r"PSD (mm$^2$/Hz)", COLOR_UX),
        ("pier base $M_z$", t_m, m_knm, r"PSD ((kN·m)$^2$/Hz)", COLOR_M),
        ("actuator $F$", t_f, f_kn, r"PSD (kN$^2$/Hz)", COLOR_F),
    )

    n_win = len(windows)
    fig_h = 3.2 * n_win * (0.65 + 0.35 * font_scale)
    fig, axes = plt.subplots(
        n_win,
        3,
        figsize=(12.6, fig_h),
        sharex=True,
        layout="constrained",
        squeeze=False,
    )

    for i_win, (wname, t0, t1) in enumerate(windows):
        for j, (sname, t, y, ylab, color) in enumerate(series):
            ax = axes[i_win, j]
            uni = interpolate_uniform(t, y, t0=t0, t1=t1)
            if uni is None:
                ax.text(0.5, 0.5, "short", transform=ax.transAxes, ha="center")
                ax.set_title(f"{sname}\n{wname}", fontsize=9)
                continue
            _tu, yu, dt = uni
            f_hz, pxx = welch_psd(yu, dt)
            ax.semilogy(f_hz, pxx, color=color, lw=1.15)
            if f_wave is not None and f_wave < float(f_hz[-1]):
                ax.axvline(
                    f_wave,
                    color=COLOR_FW,
                    ls="--",
                    lw=1.2,
                    alpha=0.85,
                    label=rf"$f_w={f_wave:.3f}$ Hz",
                )
            ax.grid(True, ls=":", alpha=0.45)
            ax.set_xlim(0.0, min(1.0, float(f_hz[-1])))
            if i_win == 0:
                ax.set_title(sname, fontsize=10)
            if j == 0:
                ax.set_ylabel(f"{wname}\n{ylab}", fontsize=9)
            if i_win == n_win - 1:
                ax.set_xlabel(r"$f$ (Hz) prototype")
            if f_wave is not None and i_win == 0 and j == 2:
                ax.legend(loc="upper right", fontsize=8, frameon=True)

    out = test_os_plots_dir(test_id) / OUT_NAME
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.suptitle(
        run_title(test_id)
        + (
            rf"  ·  wave start $t={t_wave:.1f}$ s"
            if t_wave is not None
            else ""
        ),
        y=1.02,
        fontsize=plt.rcParams["axes.labelsize"],
    )
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(
        f"PlotHydroSpectra: wrote {out}  "
        f"windows={len(windows)}  wave_start="
        f"{t_wave if t_wave is not None else 'none'}"
    )
    return 0


def default_tests() -> list[str]:
    """All Test IDs with a logged model-scale wave period."""
    out: list[str] = []
    for row in load_lab_runs_rows():
        tid = (row.get("Test") or "").strip()
        if not tid:
            continue
        if (row.get("waveT_model_s") or "").strip():
            out.append(tid)
    return out


def main() -> int:
    args = [a for a in sys.argv[1:] if a not in ("-h", "--help")]
    if any(a in ("-h", "--help") for a in sys.argv[1:]):
        print(__doc__)
        return 0
    tests = args if args else default_tests()
    if not tests:
        print("PlotHydroSpectra: no tests", file=sys.stderr)
        return 1
    rc = 0
    for tid in tests:
        rc = max(rc, write_plot(tid))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
