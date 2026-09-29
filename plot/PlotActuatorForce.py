#!/usr/bin/env python3
"""
Goals
-----
Plot OpenFresco daqForce vs time (full | D5–95) with amber lines at each
typeConv3==2 onset. Dual axes: prototype force / time, and model via /λ³
and /√λ.

  python plot/PlotActuatorForce.py --mesh-ladder
  python plot/PlotActuatorForce.py F06 F08
  python plot/PlotActuatorForce.py

Source: dump ``ServerSetup_daqFrc.out`` (OpenSees t, N, prototype).
Several Test IDs in one call share one F axis (max |F| over each full record).
Writes ``plots/runs/<Test>/os/hist_frc_actuator.png``.
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
    LABEL_T_PROTO,
    SLOWDOWN_STATE,
    add_dual_time_xaxis,
    load_type_conv3,
)
from PlotEQCompareRuns import apply_paper_style
from gm_duration import d595_proto_window, gm_start_time_s
from lab_paths import (
    CYLINDER_LENGTH_SCALE,
    LOCAL_OPENSEES_DATA,
    TIME_SCALE_FROUDE,
    XLIM_FULL_PROTO_S,
    full_xlim_proto_s,
    load_lab_runs_rows,
    resolve_opensees_data,
    test_os_plots_dir,
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
DAQ_FRC_NAME = "ServerSetup_daqFrc.out"
DEFAULT_FONT_SCALE = 1.75
# Lab wave arrives after free vib; search only for t_model >= this (Friday).
WAVE_SEARCH_MODEL_S = 200.0
# Legacy Friday spike zoom half-width (model s) when waveT is unknown.
WAVE_WINDOW_MODEL_S = 2.0
# Hydro panel: N model-scale wave periods outside strong EQ shaking.
HYDRO_N_PERIODS = 4.0
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


def load_daq_force_kn(dump_path: Path) -> tuple[np.ndarray, np.ndarray] | None:
    """
    OpenFresco daqForce history from a dump folder.

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
    f_kn = np.asarray(data[:, 1], dtype=float) * N_TO_KN
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


def write_plot(
    test_id: str,
    *,
    data_root: Path | None = None,
    font_scale: float = DEFAULT_FONT_SCALE,
    ylim_kn: tuple[float, float] | None = None,
) -> int:
    """
    Write one actuator-force PNG for a Test ID.

    Panels: full history | D5–95 | hydro (N·T_wave outside strong EQ, or
    Friday |F|-peak zoom when waveT is unknown).

    Args:    ylim_kn  shared (ymin, ymax) in kN prototype; None = auto from this run
    Returns: 0 ok, 1 skip/error
    """
    apply_paper_style()
    scale_paper_fonts(font_scale)
    pair = mat_dump_for_test(test_id)
    if pair is None:
        print(f"PlotActuatorForce: skip {test_id} (no mat+dump)", file=sys.stderr)
        return 1
    mat, dump = pair
    root = data_root or resolve_opensees_data() or LOCAL_OPENSEES_DATA
    out = test_os_plots_dir(test_id) / OUT_NAME

    frc = load_daq_force_kn(root / dump)
    if frc is None:
        print(f"PlotActuatorForce: skip {test_id} (no {DAQ_FRC_NAME})", file=sys.stderr)
        return 1
    t_frc, f_kn = frc
    t_slow = slowdown_times_proto_s(mat)
    wave_hit = detect_wave_hit_proto_s(t_frc, f_kn)
    t_wave = None if wave_hit is None else wave_hit[0]
    t_wave_peak = None if wave_hit is None else wave_hit[1]
    t0 = gm_start_time_s(root / dump)
    d595 = d595_window(t0)
    t_wave_period = wave_period_proto_s(test_id)
    hydro_xlim = hydro_force_zoom_proto(
        t_frc,
        f_kn,
        t_wave_proto=t_wave_period,
        d595=d595,
        t_wave_peak=t_wave_peak,
    )
    full_xlim = full_xlim_proto_s(t_frc)

    fig_h = 4.2 * (0.65 + 0.35 * font_scale)
    fig_w = 15.6
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(fig_w, fig_h),
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
    for ax in (ax_f, ax_z, ax_w):
        mark_fy10_ref(ax, f_y10)
        ax.plot(
            t_frc,
            f_kn,
            color=COLOR_FRC,
            lw=LW_FRC,
            label="actuator (daqForce)",
            zorder=5,
        )
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
    if ylim_kn is not None:
        ax_f.set_ylim(*ylim_kn)
    else:
        y_abs = peak_abs_kn(f_kn)
        if y_abs > 0.0:
            ax_f.set_ylim(-1.08 * y_abs, 1.08 * y_abs)

    ax_f.set_xlabel(LABEL_T_PROTO)
    ax_z.set_xlabel(LABEL_T_PROTO)
    ax_w.set_xlabel(LABEL_T_PROTO)
    ax_z.tick_params(labelleft=False)
    ax_w.tick_params(labelleft=False)
    add_dual_time_xaxis(ax_f, top=True)
    for ax in (ax_z, ax_w):
        ax.secondary_xaxis(
            "top",
            functions=(
                lambda t_proto: t_proto / TIME_SCALE_FROUDE,
                lambda t_model: t_model * TIME_SCALE_FROUDE,
            ),
        )
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

    ax_f.set_ylabel(r"$F$ (kN) prototype scale")
    sec_y = ax_w.secondary_yaxis(
        "right",
        functions=(
            lambda f_proto: f_proto / FORCE_SCALE_FROUDE,
            lambda f_model: f_model * FORCE_SCALE_FROUDE,
        ),
    )
    sec_y.set_ylabel(r"$F/\lambda^{3}$ (kN) model scale")

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

    fig.suptitle(run_title(test_id), y=1.02, fontsize=plt.rcParams["axes.labelsize"])

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    wave_txt = (
        f"wave={t_wave:.1f}s proto ({t_wave / TIME_SCALE_FROUDE:.1f}s model)"
        if t_wave is not None
        else "wave=none"
    )
    hydro_txt = (
        f"hydro=[{hydro_xlim[0]:.1f},{hydro_xlim[1]:.1f}]"
        if hydro_xlim is not None
        else "hydro=none"
    )
    print(
        f"PlotActuatorForce: wrote {out}  "
        f"(lines={n_slow}, {wave_txt}, {hydro_txt}, t_end={float(t_frc[-1]):.1f}s)"
    )
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
                raise SystemExit("PlotActuatorForce: --font-scale needs a number")
            font_scale = float(argv[i])
        elif a.startswith("--font-scale="):
            font_scale = float(a.split("=", 1)[1])
        elif a.startswith("-"):
            raise SystemExit(f"PlotActuatorForce: unknown option {a}")
        else:
            tests.append(a)
        i += 1
    return tests, font_scale


def main() -> int:
    tests, font_scale = parse_argv(sys.argv[1:])
    if not tests:
        tests = list(MESH_LADDER_TESTS)
    # Same F axis when several Test IDs are requested together.
    ylim = shared_ylim_kn(tests) if len(tests) > 1 else None
    if ylim is not None:
        print(
            f"PlotActuatorForce: shared ylim [{ylim[0]:.3g}, {ylim[1]:.3g}] kN proto"
        )
    rc = 0
    for tid in tests:
        rc = max(rc, write_plot(tid, font_scale=font_scale, ylim_kn=ylim))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
