#!/usr/bin/env python3
"""
Goals
-----
Plot OpenFresco actuator force (full | D5–95 | hydro) with amber lines at
each mid-run typeConv3==2 onset. Two clocks:

  hist_frc_actuator.png           — lab $t$ (daqForce $t_{\mathrm{int}}$
                                    → lab via atTarget map)
  hist_frc_actuator_opensees.png  — $t_{\mathrm{int}}$ (OpenSees domain)

Dual axes: prototype force, model via /λ³.

  python plot/PlotActuatorForce.py --mesh-ladder
  python plot/PlotActuatorForce.py F06 F08
  python plot/PlotActuatorForce.py

Source: dump ``ServerSetup_daqFrc.out`` (OpenSees $t_{\mathrm{int}}$, N).
Recorded daqForce is the experimental-element resisting force; plots use
F = −daqForce (external force on the numerical substructure).
Several Test IDs in one call share one F axis (max |F| over each full record).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from PlotEQComparePairs import (
    COLOR_OTHER,
    SLOWDOWN_STATE,
    load_type_conv3,
)
from PlotEQCompareRuns import apply_paper_style
from gm_duration import d595_proto_window, gm_start_time_s
from lab_paths import (
    CYLINDER_LENGTH_SCALE,
    LOCAL_OPENSEES_DATA,
    TIME_SCALE_FROUDE,
    full_xlim_model_s,
    full_xlim_proto_s,
    load_lab_runs_rows,
    resolve_opensees_data,
    test_os_plots_dir,
)
from lab_time_map import (
    LABEL_T_INT,
    LABEL_T_LAB,
    add_dual_time_xaxis_int,
    add_dual_time_xaxis_lab,
    handshake_map_for_dump,
    lab_times_to_os_handshake,
    map_os_to_lab,
    os_times_to_lab,
    os_window_to_lab,
    slowdown_lab_onsets,
)

# Froude force scale (same density): F_proto / F_model = λ³.
FORCE_SCALE_FROUDE = CYLINDER_LENGTH_SCALE**3
N_TO_KN = 1.0e-3

COLOR_FRC = COLOR_OTHER  # #001F3F
COLOR_SLOW = "#FFC04D"
COLOR_WAVE = "#C62828"  # red — distinct from navy history and amber slowdowns
COLOR_FY10 = "#546E7A"  # blue-grey reference
ALPHA_SLOW = 0.32
LW_SLOW = 1.0
LW_WAVE = 1.8
LW_FY10 = 1.2
LW_FRC = 1.15
OUT_NAME = "hist_frc_actuator.png"
OUT_NAME_OS = "hist_frc_actuator_opensees.png"
DAQ_FRC_NAME = "ServerSetup_daqFrc.out"
DEFAULT_FONT_SCALE = 1.75
# Lab wave arrives after free vib; search only for t_model >= this (Friday).
WAVE_SEARCH_MODEL_S = 200.0
# Legacy Friday spike zoom half-width (model s) when waveT is unknown.
WAVE_WINDOW_MODEL_S = 2.0
# Hydro panel: N model-scale wave periods outside strong EQ shaking.
HYDRO_N_PERIODS = 4.0
# Shared F axis for Fri+Wed: ±8 kN model → ±8·λ³ kN prototype.
F_LIM_MODEL_KN = 8.0
F_LIM_PROTO_KN = F_LIM_MODEL_KN * FORCE_SCALE_FROUDE
# Rough cantilever yield: Mn ≈ (2/π) Ast fy r ; F_y,eq = Mn / H_pier ; plot ±0.1 F_y,eq.
FY_EQ_FRAC = 0.10
FONT_SCALE_KEYS = (
    "font.size",
    "axes.labelsize",
    "xtick.labelsize",
    "ytick.labelsize",
    "legend.fontsize",
)

# Mesh ladder for the slowdown-vs-wave force story (few | many per mesh).
MESH_LADDER_TESTS: tuple[str, ...] = (
    "F05",  # Baseline, GPU, few
    "F22",  # Baseline, serial, many (stopped early)
    "F06",  # Moderate, GPU, few
    "F07",  # Moderate, CPU, many (stopped early)
    "F08",  # Moderate, CPU, many (stopped early)
    "F12",  # Large, GPU, few
    "F14",  # Large, CPU, many (complete)
    "F13",  # X-Large, GPU, many (complete)
)


def pier_fy_eq_kn() -> float:
    """
    Equivalent top force at cantilever yield (prototype kN).

    Mn ≈ (2/π) Ast fy r with Ast, fy, r = R_core from Parameters.tcl;
    F_y,eq = Mn / H_pier (H_pier = H_cyl·cylinderSF ≈ 7.25 m).

    Returns: F_y,eq in kN
    """
    inch = 0.0254
    foot = 12.0 * inch
    pi = np.pi
    d_pier = 4.0 * foot
    h_pier = 3.02 * ((4.0 * 12.0) / 20.0)  # H_cyl * cylinderSF
    fy = 470.0e6  # Pa
    cover = 2.0 * inch
    db_long = (10.0 / 8.0) * inch  # #10
    db_tran = (7.0 / 8.0) * inch  # #7 ties
    n_long = 28
    as_tot = n_long * (pi / 4.0) * db_long**2
    r_bar = 0.5 * d_pier - cover - 0.5 * db_long - db_tran  # R_core_pier
    m_n = (2.0 / pi) * as_tot * fy * r_bar  # N·m
    return float(m_n / h_pier) * N_TO_KN


def mark_fy10_ref(ax, f_kn: float) -> None:
    """± fraction of F_y,eq as light horizontal guides (under the history)."""
    for s in (+1.0, -1.0):
        ax.axhline(
            s * f_kn,
            color=COLOR_FY10,
            alpha=0.85,
            lw=LW_FY10,
            ls=":",
            zorder=2,
        )


def mark_force_exceedances(
    ax,
    t_proto: np.ndarray,
    f_kn: np.ndarray,
    lim_kn: float = F_LIM_PROTO_KN,
) -> int:
    """
    Triangles on ±lim where |F| exceeds the shared axis limit.

    Args:    ax; t_proto, f_kn; lim_kn  prototype kN half-range
    Returns: number of exceeded samples marked
    """
    if t_proto.size == 0 or f_kn.size != t_proto.size or lim_kn <= 0.0:
        return 0
    hi = np.asarray(f_kn, dtype=float) > lim_kn
    lo = np.asarray(f_kn, dtype=float) < -lim_kn
    n = 0
    if np.any(hi):
        ax.plot(
            t_proto[hi],
            np.full(int(np.count_nonzero(hi)), lim_kn),
            linestyle="none",
            marker="v",
            markersize=5.5,
            markerfacecolor="#C62828",
            markeredgecolor="#7F0000",
            markeredgewidth=0.4,
            zorder=6,
            label=rf"$|F|>{F_LIM_MODEL_KN:g}\,\mathrm{{kN}}$ model",
        )
        n += int(np.count_nonzero(hi))
    if np.any(lo):
        ax.plot(
            t_proto[lo],
            np.full(int(np.count_nonzero(lo)), -lim_kn),
            linestyle="none",
            marker="^",
            markersize=5.5,
            markerfacecolor="#C62828",
            markeredgecolor="#7F0000",
            markeredgewidth=0.4,
            zorder=6,
            label="_nolegend_" if np.any(hi) else rf"$|F|>{F_LIM_MODEL_KN:g}\,\mathrm{{kN}}$ model",
        )
        n += int(np.count_nonzero(lo))
    return n


def scale_paper_fonts(factor: float) -> None:
    """Multiply paper-style type sizes. factor=1.75 is 75% larger."""
    if factor == 1.0:
        return
    for key in FONT_SCALE_KEYS:
        plt.rcParams[key] = float(plt.rcParams[key]) * factor


def row_for_test(test_id: str) -> dict[str, str] | None:
    """TestMatrix row for a Test ID, or None."""
    for row in load_lab_runs_rows():
        if (row.get("Test") or "").strip() == test_id:
            return row
    return None


def mesh_label(soil_mesh: str) -> str:
    """soilMesh cell → Baseline | Moderate | Large | X-Large."""
    m = re.match(r"^(\d+)\s*\(([^)]+)\)", (soil_mesh or "").strip())
    if not m:
        return (soil_mesh or "?").strip() or "?"
    name = m.group(2).replace("XLARGE", "X-large").title()
    if name.lower() in ("production", "baseline"):
        return "Baseline"
    return name


def run_title(test_id: str) -> str:
    """Short figure title: Test · mesh · DOFs."""
    row = row_for_test(test_id)
    if row is None:
        return test_id
    mesh = mesh_label(row.get("soilMesh", ""))
    dofs = (row.get("DOFs") or "").strip()
    if dofs:
        return f"{test_id}  ·  {mesh}  ·  {dofs} DOF"
    return f"{test_id}  ·  {mesh}"


def mat_dump_for_test(test_id: str) -> tuple[str, str] | None:
    """Return (MatFile, DumpFolder) for Test ID, or None if incomplete."""
    row = row_for_test(test_id)
    if row is None:
        return None
    mat = (row.get("MatFile") or "").strip()
    dump = (row.get("DumpFolder") or "").strip()
    if mat and dump:
        return mat, dump
    return None


def slowdown_times_proto_s(mat_name: str) -> list[float]:
    """
    Deprecated name: lab onset times (model s), not prototype.

    Prefer ``slowdown_lab_onsets``. Kept for callers that still import this name.
    """
    return slowdown_lab_onsets(mat_name)


def mark_slowdowns(ax, t_proto: list[float]) -> int:
    """Fixed-thickness amber vertical line at each slowdown onset (plot x units)."""
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


def mark_wave(ax, t_proto: float | None) -> bool:
    """Red dashed line at detected wave onset (prototype s)."""
    if t_proto is None or not np.isfinite(t_proto):
        return False
    ax.axvline(
        float(t_proto),
        color=COLOR_WAVE,
        alpha=0.9,
        lw=LW_WAVE,
        ls="--",
        solid_capstyle="butt",
        zorder=2,
    )
    return True


def detect_wave_onset_proto_s(
    t_proto: np.ndarray,
    f_kn: np.ndarray,
    *,
    t_min_model_s: float = WAVE_SEARCH_MODEL_S,
) -> float | None:
    """Onset only; see ``detect_wave_hit_proto_s``."""
    hit = detect_wave_hit_proto_s(t_proto, f_kn, t_min_model_s=t_min_model_s)
    return None if hit is None else hit[0]


def detect_wave_hit_proto_s(
    t_proto: np.ndarray,
    f_kn: np.ndarray,
    *,
    t_min_model_s: float = WAVE_SEARCH_MODEL_S,
) -> tuple[float, float] | None:
    """
    Wave onset and |F| peak after free vib (prototype s).

    Search only t_model >= t_min_model_s. Confirm a real spike (peak ≫ quiet
    baseline), then walk back from that peak to the last sample still below a
    quiet threshold max(3·pre-wave p95, 2 kN) — the start of the rise.

    Args:    t_proto, f_kn  daqForce history; t_min_model_s  search floor (lab ~200 s)
    Returns: (t_onset, t_peak) prototype s, or None
    """
    if t_proto.size < 10 or f_kn.size != t_proto.size:
        return None
    t_min = float(t_min_model_s) * TIME_SCALE_FROUDE
    if float(t_proto[-1]) < t_min + 5.0:
        return None
    pre0 = t_min - 30.0 * TIME_SCALE_FROUDE
    m_pre = (t_proto >= pre0) & (t_proto < t_min)
    if not np.any(m_pre):
        m_pre = t_proto < t_min
    if not np.any(m_pre):
        return None
    f_base = float(np.nanpercentile(np.abs(f_kn[m_pre]), 95))
    m_post = t_proto >= t_min
    if not np.any(m_post):
        return None
    t_post = t_proto[m_post]
    f_post = np.abs(f_kn[m_post])
    i_pk = int(np.nanargmax(f_post))
    f_pk = float(f_post[i_pk])
    # Real wave spike, not free-vib chatter.
    if not np.isfinite(f_pk) or f_pk < max(15.0, 10.0 * max(f_base, 0.1)):
        return None
    # Onset = first rise out of quiet (near-zero), not mid-spike.
    thr_quiet = max(3.0 * f_base, 2.0)
    j = i_pk
    while j > 0 and float(f_post[j]) >= thr_quiet:
        j -= 1
    i_on = j + 1 if j < i_pk else i_pk
    return float(t_post[i_on]), float(t_post[i_pk])


def wave_period_proto_s(test_id: str) -> float | None:
    """
    Flume wave period on the prototype clock from the CSV, if logged.

    Args:    test_id  W## / F##
    Returns: T_wave · √λ (s), or None
    """
    row = row_for_test(test_id)
    if row is None:
        return None
    raw = (row.get("waveT_model_s") or "").strip()
    if not raw:
        return None
    try:
        t_model = float(raw)
    except ValueError:
        return None
    if t_model <= 0.0:
        return None
    return t_model * TIME_SCALE_FROUDE


def _wave_band_score(
    t_proto: np.ndarray,
    f_kn: np.ndarray,
    t0: float,
    t1: float,
    t_wave_proto: float,
) -> float:
    """FFT power near 1/T_wave in [t0, t1); 0 if the window is too short."""
    m = (t_proto >= t0) & (t_proto < t1)
    n = int(np.count_nonzero(m))
    if n < 200 or t_wave_proto <= 0.0:
        return 0.0
    y = np.asarray(f_kn[m], dtype=float)
    y = y - float(np.nanmean(y))
    dt = float(np.median(np.diff(t_proto[m])))
    if not np.isfinite(dt) or dt <= 0.0:
        return 0.0
    from numpy.fft import rfft, rfftfreq

    spectrum = np.abs(rfft(y))
    freqs = rfftfreq(n, d=dt)
    if spectrum.size < 2:
        return 0.0
    spectrum[0] = 0.0
    f0 = 1.0 / t_wave_proto
    band = (freqs >= 0.75 * f0) & (freqs <= 1.25 * f0)
    if not np.any(band):
        return 0.0
    return float(np.max(spectrum[band]))


def hydro_force_zoom_proto(
    t_proto: np.ndarray,
    f_kn: np.ndarray,
    *,
    t_wave_proto: float | None,
    d595: tuple[float, float] | None,
    t_wave_peak: float | None = None,
    n_periods: float = HYDRO_N_PERIODS,
) -> tuple[float, float] | None:
    """
    Prototype xlim for the hydro / wave panel.

    Prefer ``n_periods`` · T_wave outside D5–95 (post-EQ free vib first, then
    pre-EQ). If waveT is unknown, fall back to ±WAVE_WINDOW_MODEL_S about the
    Friday-style |F| peak.

    Args:    t_proto, f_kn; t_wave_proto  T·√λ; d595; t_wave_peak; n_periods
    Returns: (t_lo, t_hi) or None
    """
    t_end = float(t_proto[-1]) if t_proto.size else 0.0
    if t_wave_proto is not None and t_wave_proto > 0.0:
        width = float(n_periods) * float(t_wave_proto)
        if width <= 0.0 or t_end < width * 0.5:
            return None
        candidates: list[tuple[float, float]] = []
        if d595 is not None:
            t5, t95 = float(d595[0]), float(d595[1])
            # Post strong shaking (free vib + wave).
            if t_end >= t95 + 0.6 * width:
                t_lo = min(t95, t_end - width)
                t_hi = t_lo + width
                if t_hi > t_end:
                    t_hi = t_end
                    t_lo = max(t95, t_hi - width)
                candidates.append((t_lo, t_hi))
            # Late tail of a long record.
            if t_end > t95 + width:
                candidates.append((t_end - width, t_end))
            # Pre-EQ wave-only stretch.
            if t5 > width + 5.0:
                candidates.append((t5 - width, t5))
        else:
            candidates.append((max(0.0, t_end - width), t_end))
            if t_end > width + 5.0:
                candidates.append((0.0, width))
        if not candidates:
            return None
        best = max(
            candidates,
            key=lambda w: _wave_band_score(t_proto, f_kn, w[0], w[1], t_wave_proto),
        )
        # Require some hydro content; otherwise still show best post/pre window.
        return best
    # Friday-style: short window about |F| peak after free vib.
    if t_wave_peak is not None and np.isfinite(t_wave_peak):
        half = 0.5 * WAVE_WINDOW_MODEL_S * TIME_SCALE_FROUDE
        return float(t_wave_peak) - half, float(t_wave_peak) + half
    return None


def detect_wave_start_proto_s(
    t_proto: np.ndarray,
    f_kn: np.ndarray,
    t_wave_proto: float,
    *,
    n_periods: float = 3.0,
    f_p95_min_kn: float = 0.35,
    wave_ratio_min: float = 0.45,
) -> float | None:
    """
    First time the actuator force shows flume-wave content (prototype s).

    Sliding windows of ``n_periods`` · T_wave; require p95(|F|) and FFT power
    near 1/T_wave. Used for Wednesday (wave before EQ); Friday still uses the
    post–free-vib spike detector.

    Args:    t_proto, f_kn; t_wave_proto  T·√λ (s)
    Returns: onset prototype s, or None
    """
    if (
        t_proto.size < 50
        or f_kn.size != t_proto.size
        or t_wave_proto <= 0.0
    ):
        return None
    width = float(n_periods) * float(t_wave_proto)
    t_end = float(t_proto[-1])
    if t_end < width + 5.0:
        return None
    # Step ~0.25 period.
    step = max(float(t_wave_proto) * 0.25, 0.5)
    t0 = 0.0
    while t0 + width <= t_end + 1e-9:
        t1 = t0 + width
        m = (t_proto >= t0) & (t_proto < t1)
        if int(np.count_nonzero(m)) < 200:
            t0 += step
            continue
        f_win = np.asarray(f_kn[m], dtype=float)
        p95 = float(np.nanpercentile(np.abs(f_win), 95))
        if p95 < f_p95_min_kn:
            t0 += step
            continue
        score = _wave_band_score(t_proto, f_kn, t0, t1, t_wave_proto)
        # Peak PSD anywhere in the window for a ratio.
        y = f_win - float(np.nanmean(f_win))
        dt = float(np.median(np.diff(t_proto[m])))
        if not np.isfinite(dt) or dt <= 0.0:
            t0 += step
            continue
        from numpy.fft import rfft, rfftfreq

        spectrum = np.abs(rfft(y))
        spectrum[0] = 0.0
        peak = float(np.max(spectrum)) if spectrum.size else 0.0
        ratio = (score / peak) if peak > 0.0 else 0.0
        if ratio >= wave_ratio_min and score > 0.0:
            # Walk back from the detecting window like Friday's |F| spike:
            # quiet floor from samples before t0, then last rise out of quiet.
            m_pre = t_proto < t0
            if np.any(m_pre):
                f_base = float(np.nanpercentile(np.abs(f_kn[m_pre]), 95))
            else:
                f_base = 0.05
            thr_quiet = max(3.0 * f_base, 0.15)
            # Start from first sample in [t0, t1) above thr, walk earlier.
            m_win = (t_proto >= t0) & (t_proto < t1)
            idx = np.flatnonzero(m_win & (np.abs(f_kn) >= thr_quiet))
            if idx.size == 0:
                return float(t0)
            j = int(idx[0])
            while j > 0 and float(np.abs(f_kn[j])) >= thr_quiet:
                j -= 1
            i_on = j + 1 if j < idx[0] else int(idx[0])
            return float(t_proto[i_on])
        t0 += step
    return None


def load_daq_force_kn(dump_path: Path) -> tuple[np.ndarray, np.ndarray] | None:
    """
    Actuator force history as an external load on the numerical model.

    Reads OpenFresco ``daqForce`` (experimental resisting force) and returns
    F = −daqForce so plotted F is the force applied to the numerical side.

    Args:    dump_path  LOCAL opensees_data/<DumpFolder>
    Returns: (t_proto_s, F_proto_kN) or None
    """
    path = dump_path / DAQ_FRC_NAME
    if not path.is_file():
        shards = sorted(dump_path.glob(DAQ_FRC_NAME + ".*"))
        if not shards:
            return None
        path = shards[0]
    data = np.loadtxt(path, ndmin=2)
    if data.size == 0 or data.shape[1] < 2:
        return None
    t = np.asarray(data[:, 0], dtype=float)
    # daqForce = resisting; F_ext on numerical = −daqForce
    f_kn = -np.asarray(data[:, 1], dtype=float) * N_TO_KN
    return t, f_kn


def d595_window(gm_start_s: float = 0.0) -> tuple[float, float] | None:
    """GM D5–95 in prototype seconds (incl. gmStartTime), or None if VT2 missing."""
    return d595_proto_window(gm_start_s=gm_start_s)


def peak_abs_kn(f_kn: np.ndarray) -> float:
    """Peak |F| (kN proto) over the supplied samples."""
    if f_kn.size == 0:
        return 0.0
    y = float(np.nanmax(np.abs(f_kn)))
    return y if np.isfinite(y) else 0.0


def shared_ylim_kn(
    test_ids: list[str],
    *,
    data_root: Path | None = None,
) -> tuple[float, float] | None:
    """
    Symmetric ylim from the largest |F| among the listed tests (full records).

    Full-record peak so EQ/slowdown spikes and the later wave spike share one
    axis (slowdown artifact vs wave demand).

    Returns: (-pad, +pad) in kN prototype, or None if no force data
    """
    root = data_root or resolve_opensees_data() or LOCAL_OPENSEES_DATA
    y_max = 0.0
    for tid in test_ids:
        pair = mat_dump_for_test(tid)
        if pair is None:
            continue
        _, dump = pair
        frc = load_daq_force_kn(root / dump)
        if frc is None:
            continue
        y_max = max(y_max, peak_abs_kn(frc[1]))
    if y_max <= 0.0:
        return None
    pad = 1.08 * y_max
    return (-pad, pad)


def _write_force_figure(
    *,
    test_id: str,
    out: Path,
    t_x: np.ndarray,
    f_kn: np.ndarray,
    t_slow: list[float],
    t_wave: float | None,
    full_xlim: tuple[float, float] | None,
    d595: tuple[float, float] | None,
    hydro_xlim: tuple[float, float] | None,
    t_wave_period: float | None,
    ylim_kn: tuple[float, float],
    font_scale: float,
    clock: str,
    title_extra: str,
    wave_legend: str | None,
) -> tuple[int, int]:
    """
    One three-panel force figure on lab $t$ or $t_{\mathrm{int}}$.

    Args:    clock  "lab" | "int"; title_extra  after Test · …
    Returns: (n_slow, n_exceed)
    """
    apply_paper_style()
    scale_paper_fonts(font_scale)
    fig_h = 4.2 * (0.65 + 0.35 * font_scale)
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(15.6, fig_h),
        sharey=True,
        layout="constrained",
        gridspec_kw={"width_ratios": [1.35, 1.0, 1.0]},
    )
    ax_f, ax_z, ax_w = axes[0], axes[1], axes[2]

    n_slow = mark_slowdowns(ax_f, t_slow)
    mark_slowdowns(ax_z, t_slow)
    mark_slowdowns(ax_w, t_slow)
    has_wave = mark_wave(ax_f, t_wave)
    mark_wave(ax_z, t_wave)
    mark_wave(ax_w, t_wave)
    f_y10 = FY_EQ_FRAC * pier_fy_eq_kn()
    n_ex = 0
    for ax in (ax_f, ax_z, ax_w):
        mark_fy10_ref(ax, f_y10)
        ax.plot(
            t_x,
            f_kn,
            color=COLOR_FRC,
            lw=LW_FRC,
            label=r"$-F_{\mathrm{daq}}$ (on numerical)",
            zorder=5,
        )
        n_ex = mark_force_exceedances(ax, t_x, f_kn, F_LIM_PROTO_KN)
        ax.grid(True, ls=":", alpha=0.45)
        ax.axhline(0.0, color="#666666", lw=0.6, zorder=0)

    if full_xlim is not None:
        ax_f.set_xlim(*full_xlim)
    if d595 is not None:
        ax_z.set_xlim(d595[0], d595[1])
    if hydro_xlim is not None:
        ax_w.set_xlim(*hydro_xlim)
    else:
        if full_xlim is not None:
            ax_w.set_xlim(*full_xlim)
        ax_w.text(
            0.5,
            0.5,
            "no hydro window",
            transform=ax_w.transAxes,
            ha="center",
            va="center",
            color="#666666",
            fontsize=plt.rcParams["legend.fontsize"],
        )
    ax_f.set_ylim(*ylim_kn)

    xlab = LABEL_T_LAB if clock == "lab" else LABEL_T_INT
    for ax in (ax_f, ax_z, ax_w):
        ax.set_xlabel(xlab)
    ax_z.tick_params(labelleft=False)
    ax_w.tick_params(labelleft=False)
    if clock == "lab":
        for ax in (ax_f, ax_z, ax_w):
            add_dual_time_xaxis_lab(ax, top=True)
    else:
        for ax in (ax_f, ax_z, ax_w):
            add_dual_time_xaxis_int(ax, top=True)
    if hydro_xlim is not None and t_wave_period is not None:
        ax_w.set_title(
            rf"hydro (${HYDRO_N_PERIODS:g}\,T_{{\mathrm{{w}}}}$)",
            fontsize=10,
            pad=4,
        )
    elif hydro_xlim is not None:
        ax_w.set_title(r"wave $|F|$ peak", fontsize=10, pad=4)

    engine = fig.get_layout_engine()
    if engine is not None:
        engine.set(w_pad=0.015, h_pad=0.015, wspace=0.02, hspace=0.02)

    ax_f.set_ylabel(r"$-F_{\mathrm{daq}}$ (kN) prototype scale")
    sec_y = ax_w.secondary_yaxis(
        "right",
        functions=(
            lambda f_proto: f_proto / FORCE_SCALE_FROUDE,
            lambda f_model: f_model * FORCE_SCALE_FROUDE,
        ),
    )
    sec_y.set_ylabel(r"$-F_{\mathrm{daq}}/\lambda^{3}$ (kN) model scale")

    handles, labels = ax_f.get_legend_handles_labels()
    handles.append(
        Line2D(
            [0],
            [0],
            color=COLOR_FY10,
            lw=LW_FY10,
            ls=":",
            label=rf"$0.1\,F_{{y,\mathrm{{eq}}}}$ ($\pm${f_y10:.0f} kN)",
        )
    )
    labels.append(rf"$0.1\,F_{{y,\mathrm{{eq}}}}$ ($\pm${f_y10:.0f} kN)")
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
    if has_wave and t_wave is not None and wave_legend is not None:
        handles.append(
            Line2D(
                [0],
                [0],
                color=COLOR_WAVE,
                lw=LW_WAVE,
                ls="--",
                label=wave_legend,
            )
        )
        labels.append(wave_legend)
    ax_f.legend(
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
        rf"{run_title(test_id)}  ·  {title_extra}",
        y=1.02,
        fontsize=plt.rcParams["axes.labelsize"],
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return n_slow, n_ex


def write_plot(
    test_id: str,
    *,
    data_root: Path | None = None,
    font_scale: float = DEFAULT_FONT_SCALE,
    ylim_kn: tuple[float, float] | None = None,
) -> int:
    """
    Write lab-$t$ and $t_{\mathrm{int}}$ actuator-force PNGs for a Test ID.

    Panels: full history | D5–95 | hydro (N·T_wave outside strong EQ, or
    Friday |F|-peak zoom when waveT is unknown).

    Args:    ylim_kn  shared (ymin, ymax) in kN prototype; None = auto from this run
    Returns: 0 ok, 1 skip/error
    """
    pair = mat_dump_for_test(test_id)
    if pair is None:
        print(f"PlotActuatorForce: skip {test_id} (no mat+dump)", file=sys.stderr)
        return 1
    mat, dump = pair
    root = data_root or resolve_opensees_data() or LOCAL_OPENSEES_DATA
    dump_path = root / dump
    out_lab = test_os_plots_dir(test_id) / OUT_NAME
    out_os = test_os_plots_dir(test_id) / OUT_NAME_OS

    frc = load_daq_force_kn(dump_path)
    if frc is None:
        print(f"PlotActuatorForce: skip {test_id} (no {DAQ_FRC_NAME})", file=sys.stderr)
        return 1
    t_os, f_kn = frc
    try:
        mapped = map_os_to_lab(t_os, mat, dump_dir=dump_path)
    except RuntimeError as exc:
        print(f"PlotActuatorForce: skip {test_id} ({exc})", file=sys.stderr)
        return 1
    t_lab = mapped.t_lab
    t_slow_lab = slowdown_lab_onsets(mat)
    t0 = gm_start_time_s(dump_path)
    d595_os = d595_window(t0)
    t_wave_period = wave_period_proto_s(test_id)
    # Detect on t_int, then map markers / windows to lab (atTarget).
    t_wave_os = None
    t_wave_peak_os = None
    if t_wave_period is not None:
        t_wave_os = detect_wave_start_proto_s(t_os, f_kn, t_wave_period)
    wave_hit = detect_wave_hit_proto_s(t_os, f_kn)
    if wave_hit is not None:
        t_wave_peak_os = wave_hit[1]
        if t_wave_os is None:
            t_wave_os = wave_hit[0]
    hydro_os = hydro_force_zoom_proto(
        t_os,
        f_kn,
        t_wave_proto=t_wave_period,
        d595=d595_os,
        t_wave_peak=t_wave_peak_os,
    )
    d595_lab = (
        os_window_to_lab(d595_os[0], d595_os[1], mat, dump_dir=dump_path)
        if d595_os is not None
        else None
    )
    hydro_lab = (
        os_window_to_lab(hydro_os[0], hydro_os[1], mat, dump_dir=dump_path)
        if hydro_os is not None
        else None
    )
    t_wave_lab = (
        float(os_times_to_lab(t_wave_os, mat, dump_dir=dump_path)[0])
        if t_wave_os is not None
        else None
    )
    # Slowdown onsets on t_int via inverse atTarget map.
    t_slow_os: list[float] = []
    if t_slow_lab:
        try:
            hs = handshake_map_for_dump(dump_path, mat)
            t_slow_os = [
                float(x)
                for x in lab_times_to_os_handshake(
                    np.asarray(t_slow_lab, dtype=float), hs
                ).ravel()
            ]
        except RuntimeError:
            t_slow_os = []

    if ylim_kn is None:
        ylim_kn = (-F_LIM_PROTO_KN, F_LIM_PROTO_KN)

    n_slow, n_ex = _write_force_figure(
        test_id=test_id,
        out=out_lab,
        t_x=t_lab,
        f_kn=f_kn,
        t_slow=t_slow_lab,
        t_wave=t_wave_lab,
        full_xlim=full_xlim_model_s(t_lab),
        d595=d595_lab,
        hydro_xlim=hydro_lab,
        t_wave_period=t_wave_period,
        ylim_kn=ylim_kn,
        font_scale=font_scale,
        clock="lab",
        title_extra=(
            rf"atTarget map  ·  "
            rf"$f_{{\mathrm{{end}}}}={mapped.f_end*1e3:.1f}\,\mathrm{{ms}}$"
        ),
        wave_legend=(
            rf"wave start ($t={t_wave_lab:.0f}$ s lab)"
            if t_wave_lab is not None
            else None
        ),
    )
    wave_txt = f"wave={t_wave_lab:.1f}s lab" if t_wave_lab is not None else "wave=none"
    hydro_txt = (
        f"hydro=[{hydro_lab[0]:.1f},{hydro_lab[1]:.1f}] lab"
        if hydro_lab is not None
        else "hydro=none"
    )
    print(
        f"PlotActuatorForce: wrote {out_lab}  "
        f"(lines={n_slow}, exceed={n_ex}, {wave_txt}, {hydro_txt}, "
        f"t_end={float(t_lab[-1]):.1f}s lab)"
    )

    n_slow_os, n_ex_os = _write_force_figure(
        test_id=test_id,
        out=out_os,
        t_x=t_os,
        f_kn=f_kn,
        t_slow=t_slow_os,
        t_wave=float(t_wave_os) if t_wave_os is not None else None,
        full_xlim=full_xlim_proto_s(t_os),
        d595=d595_os,
        hydro_xlim=hydro_os,
        t_wave_period=t_wave_period,
        ylim_kn=ylim_kn,
        font_scale=font_scale,
        clock="int",
        title_extra=r"$t_{\mathrm{int}}$ (OpenSees domain)",
        wave_legend=(
            rf"wave start ($t_{{\mathrm{{int}}}}={float(t_wave_os):.0f}$ s)"
            if t_wave_os is not None
            else None
        ),
    )
    print(
        f"PlotActuatorForce: wrote {out_os}  "
        f"(lines={n_slow_os}, exceed={n_ex_os}, "
        f"t_end={float(t_os[-1]):.1f}s int)"
    )
    return 0


def parse_argv(argv: list[str]) -> tuple[list[str], float, bool]:
    """Split Test IDs from ``--font-scale`` / ``--mesh-ladder`` / ``--auto-ylim``."""
    tests: list[str] = []
    font_scale = DEFAULT_FONT_SCALE
    auto_ylim = False
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in ("-h", "--help"):
            print(__doc__)
            raise SystemExit(0)
        if a == "--mesh-ladder":
            tests.extend(MESH_LADDER_TESTS)
        elif a == "--auto-ylim":
            auto_ylim = True
        elif a == "--font-scale":
            i += 1
            if i >= len(argv):
                raise SystemExit("PlotActuatorForce: --font-scale needs a number")
            font_scale = float(argv[i])
        elif a.startswith("--font-scale="):
            font_scale = float(a.split("=", 1)[1])
        elif a.startswith("-"):
            raise SystemExit(f"PlotActuatorForce: unknown option {a}")
        else:
            tests.append(a)
        i += 1
    return tests, font_scale, auto_ylim


def main() -> int:
    tests, font_scale, auto_ylim = parse_argv(sys.argv[1:])
    if not tests:
        tests = list(MESH_LADDER_TESTS)
    if auto_ylim:
        ylim = shared_ylim_kn(tests) if len(tests) > 1 else None
        if ylim is not None:
            print(
                f"PlotActuatorForce: auto shared ylim "
                f"[{ylim[0]:.3g}, {ylim[1]:.3g}] kN proto"
            )
        else:
            print("PlotActuatorForce: per-run auto ylim")
    else:
        ylim = (-F_LIM_PROTO_KN, F_LIM_PROTO_KN)
        print(
            f"PlotActuatorForce: ylim [{ylim[0]:.3g}, {ylim[1]:.3g}] kN proto "
            f"(±{F_LIM_MODEL_KN:g} kN model)"
        )
    rc = 0
    for tid in tests:
        rc = max(rc, write_plot(tid, font_scale=font_scale, ylim_kn=ylim))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
