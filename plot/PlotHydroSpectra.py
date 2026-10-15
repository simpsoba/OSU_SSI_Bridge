#!/usr/bin/env python3
"""
Goals
-----
Welch PSDs as the dissertation figures, Hann segments, segment means removed, density scaling;
no resample / interpolate. Friday (EQ then solitary wave): 40 s prototype segments, 75 % overlap
(plot/dissertation/fig_fri_psd_phases.py). Wednesday (periodic waves): three segments with 50 %
overlap tiling each window, which resolves the wave harmonics n f_w (plot/dissertation/fig_wed_psd.py).

Figure 1–2 (OpenSees domain \(t\)): structure + SSI recorder PSDs
(``psd_pier_ux_M_theta_F.png``, ``psd_pileM_py_soil.png``).

Figure 3 (``psd_tar_com_mea_lab_vs_pier_os.png``): tar/com/mea on lab \(t\) vs
pier \(u_x\) on \(t_{\mathrm{int}}\) (each native grid; lab \(f\) → proto /√λ).

Figure 4 (``psd_frc_actuator.png``): \(-F_{\mathrm{daq}}\) on \(t_{\mathrm{int}}\)
only (uniform recorder dt; no map / resample).

Columns = analysis windows (Wed: pre-EQ hydro | D5–95 | post-EQ + wave;
Fri: D5–95 | post-EQ free vib | wave) in the same time unit as each series.

Mode guides: eigenfrequencies f_1 and f_4 (fundamental and pier modes) of the numerical
model, from ``plot/out/eigen/mesh<soilMesh>/``. Higher soil modes (e.g. the 7.99 Hz
pair) have no interface participation and no Welch peak, so they are not marked.

  python plot/PlotHydroSpectra.py W05 F05
  python plot/PlotHydroSpectra.py --lab-vs-pier-os F06
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
    detect_wave_hit_proto_s,
    detect_wave_start_proto_s,
    load_daq_force_kn,
    mat_dump_for_test,
    row_for_test,
    run_title,
    scale_paper_fonts,
    wave_period_proto_s,
)
from PlotEQ import loadtxt_partial
from PlotEQCompareRuns import apply_paper_style, model_disp_to_proto_mm
from PlotMatOS import load_block, non_time_cols
from PlotPierBaseForce import load_pier_base_pm
from gm_duration import d595_proto_window, gm_start_time_s
from lab_paths import (
    CYLINDER_LENGTH_SCALE,
    LOCAL_OPENSEES_DATA,
    MAT_EXTRACT_DIR,
    TIME_SCALE_FROUDE,
    load_lab_runs_rows,
    resolve_opensees_data,
    test_os_plots_dir,
)
from lab_time_map import os_times_to_lab, os_window_to_lab

OUT_NAME = "psd_pier_ux_M_theta_F.png"
OUT_NAME_SSI = "psd_pileM_py_soil.png"
OUT_NAME_LAB_VS_PIER = "psd_tar_com_mea_lab_vs_pier_os.png"
OUT_NAME_FRC = "psd_frc_actuator.png"
COLOR_UX = "#1565C0"
COLOR_TAR = "#6A1B9A"
COLOR_COM = "#2A9D8F"
COLOR_MEA = "#001F3F"
COLOR_PIER_OS = "#B0B0B0"
COLOR_M = "#6A1B9A"
COLOR_F = "#001F3F"
COLOR_TH = "#00695C"
COLOR_PY = "#E65100"
COLOR_SOIL = "#5D4037"
COLOR_FW = "#C62828"
COLOR_F1 = "#2E7D32"
COLOR_FACT = "#6D4C41"
DEFAULT_PIER_T1_S = 2.09577  # gravity+holdPier eigen (mesh 0 and 19)
DEFAULT_PIER_T4_S = 0.3361
TSEG_PROTO_S = 40.0  # Welch segment, prototype s (as the dissertation PSD figures)
# Per-mesh archive from DumpEigenModes.tcl (avoids re-running eigen for PSD guides).
EIGEN_DIR = Path(__file__).resolve().parent / "out" / "eigen"
EIGEN_JSON_FALLBACK = Path(__file__).resolve().parent / "out" / "eigen" / "mesh0" / "eigen_modes.json"  # baseline
# Dense 3×N grid: keep near paper base (actuator force uses 1.75 for 1-row figs).
DEFAULT_FONT_SCALE_PSD = 1.15
F_MAX_PROTO_HZ = 100.0
F_MIN_PROTO_HZ = 3.0e-2
# Seki et al. (2026): transfer-system fn ≈ 24–28 Hz model → proto via √λ.
F_ACT_MODEL_HZ = (24.0, 28.0)
# Prefer this spring station when present (lean dumps include iy=10).
PY_IY_PREF = 10
HINGE_DEFO_NAME = "pier_hinge_defo.out"
LAMBDA = CYLINDER_LENGTH_SCALE
# S_model = S_proto / (amp_scale^2 · √λ); amp_scale = y_proto/y_model.
PSD_UX_TO_MODEL = 1.0 / (LAMBDA**2 * TIME_SCALE_FROUDE)  # mm²
PSD_M_TO_MODEL = 1.0 / (LAMBDA**8 * TIME_SCALE_FROUDE)  # (kN·m)²
PSD_F_TO_MODEL = 1.0 / (LAMBDA**6 * TIME_SCALE_FROUDE)  # kN²
# θ and γ are dimensionless; only the Hz axis Froude-scales.
PSD_THETA_TO_MODEL = 1.0 / TIME_SCALE_FROUDE
PSD_GAMMA_TO_MODEL = 1.0 / TIME_SCALE_FROUDE


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


def _recorder_paths(dump_path: Path, stem: str) -> list[Path]:
    """Serial ``stem`` or parallel ``stem.*`` shards (sorted)."""
    plain = dump_path / stem
    if plain.is_file():
        return [plain]
    return sorted(dump_path.glob(stem + ".*"))


def _eles_path_for_recorder(rec_path: Path, eles_stem: str) -> Path | None:
    """Matching ``eles_stem`` / ``eles_stem.<rank>`` next to a recorder shard."""
    name = rec_path.name
    if ".out." in name:
        rank = name.rsplit(".out.", 1)[1]
        cand = rec_path.with_name(f"{eles_stem}.{rank}")
        if cand.is_file():
            return cand
    cand = rec_path.with_name(eles_stem)
    return cand if cand.is_file() else None


def _read_ele_rows(path: Path) -> list[list[str]]:
    """Whitespace-split non-comment rows from an element list file."""
    rows: list[list[str]] = []
    with path.open(encoding="utf-8", errors="replace") as fd:
        for ln in fd:
            s = ln.strip()
            if not s or s.startswith("#"):
                continue
            rows.append(s.split())
    return rows


def load_pier_base_theta_rad(
    dump_path: Path,
) -> tuple[np.ndarray, np.ndarray] | None:
    """
    Bottom ZLS section rotation θ (rad) vs OpenSees t.

    ``pier_hinge_defo`` = (ε, θ) for lumpedPlasticity zeroLength section.
    """
    shards = _recorder_paths(dump_path, HINGE_DEFO_NAME)
    if not shards:
        return None
    data = loadtxt_partial(shards[0])
    if data.size == 0 or data.shape[1] < 3:
        return None
    t = np.asarray(data[:, 0], dtype=float)
    theta = np.asarray(data[:, 2], dtype=float)
    return t, theta


def load_center_pile_M_knm(
    dump_path: Path,
) -> tuple[np.ndarray, np.ndarray, str] | None:
    """
    Mid-depth center-pile i-end Mz (kN·m).

    Lean dumps keep the center shaft only; pick the mid iSeg.
    Returns: (t, M_kNm, label) or None
    """
    recs = _recorder_paths(dump_path, "pile_beam_globalForce.out")
    if not recs:
        return None
    rec = max(recs, key=lambda p: p.stat().st_size)
    eles_p = _eles_path_for_recorder(rec, "pile_beam_eles.txt")
    if eles_p is None:
        return None
    rows = _read_ele_rows(eles_p)
    if not rows:
        return None
    # rows: eleTag ip iSeg
    parsed = [(int(r[0]), int(r[1]), int(r[2])) for r in rows if len(r) >= 3]
    if not parsed:
        return None
    segs = sorted({s for _e, _ip, s in parsed})
    mid_seg = segs[len(segs) // 2]
    k = next(i for i, (_e, _ip, s) in enumerate(parsed) if s == mid_seg)
    data = loadtxt_partial(rec)
    if data.size == 0:
        return None
    n_ele = len(parsed)
    ncomp = (data.shape[1] - 1) // n_ele
    if ncomp < 6:
        return None
    # globalForce 2D: Fx Fy Mz | Fx Fy Mz  → i-end Mz = index 2
    m_nm = np.asarray(data[:, 1 + k * ncomp + 2], dtype=float)
    t = np.asarray(data[:, 0], dtype=float)
    label = f"pile $M_z$ (iSeg={mid_seg})"
    return t, m_nm * 1.0e-3, label


def load_py_force_kn(
    dump_path: Path,
) -> tuple[np.ndarray, np.ndarray, str] | None:
    """
    One center-pile p-y spring force (dir 1, kN), prefer iy=PY_IY_PREF.

    Returns: (t, F_kN, label) or None
    """
    best: tuple[int, Path, int, int] | None = None  # |iy-pref|, force_path, col, iy
    for eles_p in sorted(dump_path.glob("pile_springs_eles.txt*")):
        suffix = eles_p.name.replace("pile_springs_eles.txt", "")
        force_p = dump_path / f"pile_springs_force.out{suffix}"
        if not force_p.is_file():
            # serial: pile_springs_force.out next to pile_springs_eles.txt
            if suffix == "":
                force_p = dump_path / "pile_springs_force.out"
            if not force_p.is_file():
                continue
        rows = _read_ele_rows(eles_p)
        # eleTag kind ip iy
        local: list[tuple[int, int]] = []
        for r in rows:
            if len(r) < 4:
                continue
            if r[1] != "pile":
                continue
            local.append((int(r[2]), int(r[3])))  # ip, iy
        if not local:
            continue
        # column layout: per spring (F1, F2) = p-y, t-z
        for j, (_ip, iy) in enumerate(local):
            score = abs(iy - PY_IY_PREF)
            if best is None or score < best[0]:
                best = (score, force_p, j, iy)
    if best is None:
        return None
    _score, force_p, j, iy = best
    data = loadtxt_partial(force_p)
    if data.size == 0 or data.shape[1] < 1 + 2 * (j + 1):
        return None
    t = np.asarray(data[:, 0], dtype=float)
    f_n = np.asarray(data[:, 1 + 2 * j], dtype=float)
    return t, f_n * 1.0e-3, f"p-y $F$ (iy={iy})"


def load_soil_gamma(
    dump_path: Path,
) -> tuple[np.ndarray, np.ndarray, str] | None:
    """
    Center-column soil shear strain γ (mid list) vs OpenSees t.

    ``window_quad_strain`` stores (εx, εy, γ) per quad (1 GP for SSPquad).
    Returns: (t, gamma, label) or None
    """
    best: tuple[Path, Path] | None = None
    for quads_p in sorted(dump_path.glob("window_quads.txt*")):
        suffix = quads_p.name.replace("window_quads.txt", "")
        strain_p = dump_path / f"window_quad_strain.out{suffix}"
        if not strain_p.is_file():
            if suffix == "":
                strain_p = dump_path / "window_quad_strain.out"
            if not strain_p.is_file():
                continue
        rows = _read_ele_rows(quads_p)
        # Prefer center column if annotated.
        center = [r for r in rows if len(r) >= 3 and r[2] == "center"]
        use = center if center else rows
        if not use:
            continue
        # Keep the largest center-column shard.
        if best is None or strain_p.stat().st_size > best[0].stat().st_size:
            best = (strain_p, quads_p)
    if best is None:
        return None
    strain_p, quads_p = best
    rows = _read_ele_rows(quads_p)
    center = [r for r in rows if len(r) >= 3 and r[2] == "center"]
    use = center if center else rows
    k = len(use) // 2
    layer = use[k][1] if len(use[k]) > 1 else "?"
    # Map k into full row index for column offset when only center is selected
    # but file has all quads from this rank in file order (= rows order).
    if center and len(center) < len(rows):
        # strain file on this rank matches `rows` order for this shard
        full_idx = [i for i, r in enumerate(rows) if len(r) >= 3 and r[2] == "center"]
        if not full_idx:
            return None
        k_file = full_idx[len(full_idx) // 2]
        layer = rows[k_file][1] if len(rows[k_file]) > 1 else layer
    else:
        k_file = k
    data = loadtxt_partial(strain_p)
    if data.size == 0:
        return None
    n_ele = len(rows)
    ncomp = (data.shape[1] - 1) // max(n_ele, 1)
    if ncomp < 3 or k_file >= n_ele:
        return None
    # γ is component 2 of (εx, εy, γ)
    gamma = np.asarray(data[:, 1 + k_file * ncomp + 2], dtype=float)
    t = np.asarray(data[:, 0], dtype=float)
    return t, gamma, rf"soil $\gamma$ ({layer})"


def soil_mesh_id(test_id: str) -> int | None:
    """
    soilMesh integer from TestMatrix (e.g. ``0 (BASELINE)`` → 0).

    Args:    test_id
    Returns: mesh id, or None
    """
    row = row_for_test(test_id)
    if row is None:
        return None
    raw = (row.get("soilMesh") or "").strip()
    if not raw:
        return None
    try:
        return int(raw.split()[0])
    except ValueError:
        return None


def eigen_json_path(mesh: int | None) -> Path:
    """Preferred eigen JSON for a soilMesh id."""
    if mesh is not None:
        p = EIGEN_DIR / f"mesh{mesh}" / "eigen_modes.json"
        if p.is_file():
            return p
    return EIGEN_JSON_FALLBACK


def _eigen_period_s(mode: int, default: float, *, mesh: int | None) -> float:
    """Read T (s) for a mode from the mesh-keyed eigen JSON, else default."""
    path = eigen_json_path(mesh)
    if not path.is_file():
        return default
    try:
        import json

        data = json.loads(path.read_text(encoding="utf-8"))
        for m in data.get("modes_meta") or []:
            if int(m.get("mode", -1)) == mode:
                t = float(m["T"])
                if t > 0.0:
                    return t
    except (OSError, ValueError, KeyError, TypeError):
        pass
    return default


def pier_f1_hz(test_id: str) -> float:
    """
    Fundamental frequency from gravity+holdPier eigen for this run's mesh.

    ``rayleighT1`` in the CSV is the Rayleigh damping target (2.3 s), not the
    eigen period.

    Args:    test_id
    Returns: f1 (Hz, prototype / analysis clock)
    """
    mesh = soil_mesh_id(test_id)
    return 1.0 / _eigen_period_s(1, DEFAULT_PIER_T1_S, mesh=mesh)


def pier_f4_hz(test_id: str) -> float:
    """Mode-4 frequency for this run's mesh (Hz proto)."""
    mesh = soil_mesh_id(test_id)
    return 1.0 / _eigen_period_s(4, DEFAULT_PIER_T4_S, mesh=mesh)


def is_wed(test_id: str | None) -> bool:
    """Wednesday (periodic-wave) runs use the tiled Welch segmentation."""
    return bool(test_id) and str(test_id).upper().startswith("W")


def welch_psd(
    y: np.ndarray,
    dt: float,
    tseg: float | None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    One-sided Welch PSD, as the dissertation figures.

    Hann segments of ``tseg`` seconds (whole window if shorter) with 75 % overlap, or, with
    ``tseg=None``, three segments with 50 % overlap tiling the window (Wednesday periodic waves);
    segment means removed, density scaling: units of y²/Hz.

    Args:    y  series; dt  sampling interval (s); tseg  segment length (s, same clock as dt) or None
    Returns: (f_Hz, PSD)
    """
    y = np.asarray(y, dtype=float)
    if y.size < 64 or not np.isfinite(dt) or dt <= 0.0:
        return np.asarray([]), np.asarray([])
    if tseg is None:
        n = y.size // 2
        n -= n % 2
        nov = n // 2
    else:
        n = min(int(round(tseg / dt)), y.size)
        nov = 3 * n // 4
    return welch(
        y,
        fs=1.0 / dt,
        window="hann",
        nperseg=n,
        noverlap=nov,
        detrend="constant",
        return_onesided=True,
        scaling="density",
    )


def add_dual_freq_xaxis(ax, *, label: bool = True) -> None:
    """
    Top x-axis: model-scale frequency f_model = f_proto · √λ.

    Call only on the top row. Always set the xlabel (hide on side
    columns) so constrained layout keeps equal top margins.
    """
    sec = ax.secondary_xaxis(
        "top",
        functions=(
            lambda f_proto: f_proto * TIME_SCALE_FROUDE,
            lambda f_model: f_model / TIME_SCALE_FROUDE,
        ),
    )
    label_fs = float(plt.rcParams["axes.labelsize"])
    tick_fs = float(plt.rcParams["xtick.labelsize"])
    sec.tick_params(axis="x", labelsize=tick_fs, pad=2, length=3.5)
    sec.set_xlabel(r"$f\sqrt{\lambda}$ (Hz) model scale", labelpad=3)
    sec.xaxis.label.set_size(label_fs)
    sec.xaxis.label.set_visible(bool(label))


def mark_freq_guides(
    ax,
    *,
    f_wave: float | None,
    f_pier: float | None,
    f_mode4: float | None,
    f_max: float,
    with_labels: bool,
) -> None:
    """
    Vertical guides (prototype Hz): wave + superharmonics, the eigenfrequencies f1 and f4
    (fundamental and pier modes) of the numerical model, and the Seki transfer-system band.

    The measured pier-mode peak sits slightly below f4 during strong shaking (about 2.85 Hz,
    softening); the guides mark the numerical model, not the peaks.
    """
    # Actuation / transfer system (Seki et al. 2026): 24–28 Hz model.
    f_lo = F_ACT_MODEL_HZ[0] / TIME_SCALE_FROUDE
    f_hi = F_ACT_MODEL_HZ[1] / TIME_SCALE_FROUDE
    if f_lo < f_max:
        ax.axvspan(
            f_lo,
            min(f_hi, f_max),
            color=COLOR_FACT,
            alpha=0.18,
            lw=0,
            zorder=0,
            label=(rf"$f_\mathrm{{act}}$" if with_labels else None),
        )
    if f_wave is not None and f_wave > 0.0:
        for n, ls, alpha in (
            (1, "--", 0.95),
            (2, "-.", 0.85),
            (3, ":", 0.8),
        ):
            fn = n * f_wave
            if not (0.0 < fn < f_max):
                continue
            if n == 1:
                lab = rf"$f_w$" if with_labels else None
            else:
                lab = rf"${n}f_w$" if with_labels else None
            ax.axvline(
                fn,
                color=COLOR_FW,
                ls=ls,
                lw=1.15 if n == 1 else 1.0,
                alpha=alpha,
                label=lab,
            )
    if f_pier is not None and 0.0 < f_pier < f_max:
        ax.axvline(
            f_pier,
            color=COLOR_F1,
            ls="-",
            lw=1.35,
            alpha=0.9,
            label=(rf"$f_1$" if with_labels else None),
        )
    if f_mode4 is not None and 0.0 < f_mode4 < f_max:
        ax.axvline(
            f_mode4,
            color=COLOR_F1,
            ls="--",
            lw=1.05,
            alpha=0.75,
            label=(rf"$f_4$" if with_labels else None),
        )


def add_dual_psd_yaxis(ax, scale_to_model: float, ylabel: str) -> None:
    """
    Right y-axis: model-scale PSD density.

    Args:    scale_to_model  S_m = S_p · scale; ylabel  model-unit string
    """
    sec = ax.secondary_yaxis(
        "right",
        functions=(
            lambda s_p: s_p * scale_to_model,
            lambda s_m: s_m / scale_to_model,
        ),
    )
    label_fs = float(plt.rcParams["axes.labelsize"])
    tick_fs = float(plt.rcParams["ytick.labelsize"])
    sec.tick_params(axis="y", labelsize=tick_fs, pad=2)
    if ylabel:
        sec.set_ylabel(ylabel, labelpad=4)
        sec.yaxis.label.set_size(label_fs)


def analysis_windows(
    t_end: float,
    d595: tuple[float, float] | None,
    *,
    t_wave_wed: float | None,
    t_wave_fri: float | None,
) -> list[tuple[str, float, float]]:
    """
    Named PSD windows (same units as the series time axis).

    Wed (wave then EQ): pre-EQ hydro | D5–95 | post-EQ + wave.
    Fri (EQ then wave): D5–95 | post-EQ free vib | wave.
    """
    out: list[tuple[str, float, float]] = []
    if d595 is None:
        out.append(("full record", 0.0, t_end))
        return out
    t5, t95 = float(d595[0]), float(d595[1])
    t95_c = min(t95, t_end)

    # Wednesday: spectral wave onset before GM D5–95.
    if t_wave_wed is not None and t_wave_wed < t5 - 5.0:
        out.append(("pre-EQ hydro", float(t_wave_wed), t5))
        out.append(("D5–95", t5, t95_c))
        if t_end > t95 + 20.0:
            out.append(("post-EQ + wave", t95, t_end))
        return out

    # Friday (or Wed without early wave): EQ first.
    out.append(("D5–95", t5, t95_c))
    if t_wave_fri is not None and t_wave_fri > t95 + 5.0 and t_wave_fri < t_end - 5.0:
        out.append(("post-EQ free vib", t95, float(t_wave_fri)))
        out.append(("wave", float(t_wave_fri), t_end))
    elif t_end > t95 + 20.0:
        out.append(("post-EQ", t95, t_end))
    return out


def save_psd_figure(
    *,
    test_id: str,
    out_path: Path,
    windows: list[tuple[str, float, float]],
    series: list[tuple],
    f_wave: float | None,
    f_pier: float | None,
    f_mode4: float | None,
    t_wave_label: float | None,
    font_scale: float,
    title_extra: str = "",
) -> None:
    """
    Write one OpenSees-domain PSD grid PNG (Welch, no resample).

    series items: (sname, t_os, y, ylab_p, ylab_m, s2m, color)
    """
    n_win = len(windows)
    n_sig = len(series)
    if n_sig == 0 or n_win == 0:
        return
    fig_w = max(4.4 * n_win, 11.0)
    fig_h = 3.35 * n_sig * (0.7 + 0.3 * font_scale)
    fig, axes = plt.subplots(
        n_sig,
        n_win,
        figsize=(fig_w, fig_h),
        sharex=True,
        sharey="row",
        layout="constrained",
        squeeze=False,
    )
    fig.set_constrained_layout_pads(
        w_pad=0.02, h_pad=0.02, hspace=0.06, wspace=0.04
    )

    for i_sig, (sname, t, y, ylab_p, ylab_m, s2m, color) in enumerate(series):
        for j_win, (wname, t0, t1) in enumerate(windows):
            ax = axes[i_sig, j_win]
            got = _psd_proto_hz(t, y, t0, t1, lab_clock=False, tiled=is_wed(test_id))
            if got is None:
                ax.text(0.5, 0.5, "short", transform=ax.transAxes, ha="center")
            else:
                f_hz, pxx, f_hi = got
                ax.loglog(f_hz, pxx, color=color, lw=1.15, label=sname)
                mark_freq_guides(
                    ax,
                    f_wave=f_wave,
                    f_pier=f_pier,
                    f_mode4=f_mode4,
                    f_max=f_hi * 1.05,
                    with_labels=(i_sig == 0 and j_win == n_win - 1),
                )
                ax.set_xlim(F_MIN_PROTO_HZ, f_hi)
            ax.grid(True, which="both", ls=":", alpha=0.45)
            tick_fs = float(plt.rcParams["xtick.labelsize"])
            ax.tick_params(labelsize=tick_fs)
            if j_win == 0:
                ax.set_ylabel(f"{sname}\n{ylab_p}")
            if j_win == n_win - 1:
                add_dual_psd_yaxis(ax, s2m, ylab_m)
            if i_sig == 0:
                add_dual_freq_xaxis(ax, label=(j_win == n_win // 2))
            ax.text(
                0.02,
                0.96,
                wname,
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=tick_fs,
            )
            if i_sig == n_sig - 1 and j_win == n_win // 2:
                ax.set_xlabel(r"$f$ (Hz) prototype scale")
            if i_sig == 0 and j_win == n_win - 1:
                ax.legend(
                    loc="lower left",
                    frameon=True,
                    fontsize=float(plt.rcParams["legend.fontsize"]),
                )

    mesh = soil_mesh_id(test_id)
    mesh_bit = f"  ·  mesh{mesh}" if mesh is not None else ""
    fig.suptitle(
        run_title(test_id)
        + mesh_bit
        + title_extra
        + r"  ·  Welch PSD ($t_{\mathrm{int}}$)"
        + (
            rf"  ·  wave start $t={t_wave_label:.1f}$ s"
            if t_wave_label is not None
            else ""
        )
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"PlotHydroSpectra: wrote {out_path}  windows={len(windows)}")


def write_plot(test_id: str, *, font_scale: float = DEFAULT_FONT_SCALE_PSD) -> int:
    """Write structure + SSI PSD PNGs on $t_{\mathrm{int}}$ (no time map)."""
    apply_paper_style()
    scale_paper_fonts(font_scale)
    label_fs = float(plt.rcParams["axes.labelsize"])
    plt.rcParams["axes.titlesize"] = label_fs
    plt.rcParams["figure.titlesize"] = label_fs * 1.15
    pair = mat_dump_for_test(test_id)
    if pair is None:
        print(f"PlotHydroSpectra: skip {test_id} (no mat+dump)", file=sys.stderr)
        return 1
    _mat, dump = pair
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
    th = load_pier_base_theta_rad(dump_path)
    pile_m = load_center_pile_M_knm(dump_path)
    py = load_py_force_kn(dump_path)
    soil_g = load_soil_gamma(dump_path)

    t_wave_period = wave_period_proto_s(test_id)
    t_wave_wed = (
        detect_wave_start_proto_s(t_f, f_kn, t_wave_period)
        if t_wave_period is not None
        else None
    )
    wave_hit = detect_wave_hit_proto_s(t_f, f_kn)
    t_wave_fri = None if wave_hit is None else float(wave_hit[0])
    t_wave_label = t_wave_wed if t_wave_wed is not None else t_wave_fri
    gm0 = gm_start_time_s(dump_path)
    d595 = d595_proto_window(gm0)
    t_end = float(min(float(t_f[-1]), float(t_u[-1]), float(t_m[-1])))
    windows = analysis_windows(
        t_end,
        d595,
        t_wave_wed=t_wave_wed,
        t_wave_fri=t_wave_fri,
    )
    f_wave = (1.0 / t_wave_period) if t_wave_period else None
    f_pier = pier_f1_hz(test_id)
    f_mode4 = pier_f4_hz(test_id)
    mesh = soil_mesh_id(test_id)
    print(
        f"PlotHydroSpectra: {test_id}  mesh={mesh}  "
        f"f1={f_pier:.4f}  f4={f_mode4:.4f} Hz  "
        f"eigen={eigen_json_path(mesh).name}"
    )

    series_struct: list[tuple] = [
        (
            "pier top $u_x$",
            t_u,
            ux_mm,
            r"(mm$^2$/Hz) proto",
            r"(mm$^2$/Hz) model",
            PSD_UX_TO_MODEL,
            COLOR_UX,
        ),
        (
            "pier base $M_z$",
            t_m,
            m_knm,
            r"((kN$\cdot$m)$^2$/Hz) proto",
            r"((kN$\cdot$m)$^2$/Hz) model",
            PSD_M_TO_MODEL,
            COLOR_M,
        ),
    ]
    if th is not None:
        t_th, theta = th
        series_struct.append(
            (
                r"pier base $\theta$",
                t_th,
                theta,
                r"(rad$^2$/Hz) proto",
                r"(rad$^2$/Hz) model",
                PSD_THETA_TO_MODEL,
                COLOR_TH,
            )
        )
    else:
        print(f"PlotHydroSpectra: {test_id} no pier_hinge_defo — skip θ row")
    series_struct.append(
        (
            r"$-F_{\mathrm{daq}}$ (on numerical)",
            t_f,
            f_kn,
            r"(kN$^2$/Hz) proto",
            r"(kN$^2$/Hz) model",
            PSD_F_TO_MODEL,
            COLOR_F,
        )
    )

    out_dir = test_os_plots_dir(test_id)
    save_psd_figure(
        test_id=test_id,
        out_path=out_dir / OUT_NAME,
        windows=windows,
        series=series_struct,
        f_wave=f_wave,
        f_pier=f_pier,
        f_mode4=f_mode4,
        t_wave_label=t_wave_label,
        font_scale=font_scale,
    )

    series_ssi: list[tuple] = []
    if pile_m is not None:
        t_pm, m_pile, lab_pm = pile_m
        series_ssi.append(
            (
                lab_pm,
                t_pm,
                m_pile,
                r"((kN$\cdot$m)$^2$/Hz) proto",
                r"((kN$\cdot$m)$^2$/Hz) model",
                PSD_M_TO_MODEL,
                COLOR_M,
            )
        )
    if py is not None:
        t_py, f_py, lab_py = py
        series_ssi.append(
            (
                lab_py,
                t_py,
                f_py,
                r"(kN$^2$/Hz) proto",
                r"(kN$^2$/Hz) model",
                PSD_F_TO_MODEL,
                COLOR_PY,
            )
        )
    if soil_g is not None:
        t_g, gamma, lab_g = soil_g
        series_ssi.append(
            (
                lab_g,
                t_g,
                gamma,
                r"([-]$^2$/Hz) proto",
                r"([-]$^2$/Hz) model",
                PSD_GAMMA_TO_MODEL,
                COLOR_SOIL,
            )
        )
    if series_ssi:
        save_psd_figure(
            test_id=test_id,
            out_path=out_dir / OUT_NAME_SSI,
            windows=windows,
            series=series_ssi,
            f_wave=f_wave,
            f_pier=f_pier,
            f_mode4=f_mode4,
            t_wave_label=t_wave_label,
            font_scale=font_scale,
            title_extra="  ·  SSI",
        )
    else:
        print(f"PlotHydroSpectra: {test_id} no pile/p-y/soil recorders — skip SSI fig")
    return 0


def _window_uniform_dt(
    t: np.ndarray,
    y: np.ndarray,
    t0: float,
    t1: float,
) -> tuple[np.ndarray, float] | None:
    """
    Time-window cut on an already-uniform series (no interp, keep mean).

    Args:    t,y; t0,t1  window
    Returns: (y_win, dt) or None
    """
    m = (t >= t0) & (t <= t1) & np.isfinite(y)
    if int(np.count_nonzero(m)) < 64:
        return None
    yy = np.asarray(y[m], dtype=float)
    tt = np.asarray(t[m], dtype=float)
    dt = float(np.median(np.diff(tt)))
    if not np.isfinite(dt) or dt <= 0.0:
        return None
    return yy, dt


def _psd_proto_hz(
    t: np.ndarray,
    y: np.ndarray,
    t0: float,
    t1: float,
    *,
    lab_clock: bool,
    tiled: bool = False,
) -> tuple[np.ndarray, np.ndarray, float] | None:
    """
    Welch PSD on the native sample grid (no resample / interp); 40 s prototype segments.

    Args:    t,y  series; t0,t1  window; tiled  Wednesday segmentation; lab_clock  True → t is lab Time
             (convert f → proto via /√λ); False → OpenSees / proto t
    Returns: (f_proto_Hz, Pxx_proto, f_hi) or None
    """
    got = _window_uniform_dt(t, y, t0, t1)
    if got is None:
        return None
    yu, dt_use = got
    if lab_clock:
        f_lab, pxx_lab = welch_psd(yu, dt_use, None if tiled else TSEG_PROTO_S / TIME_SCALE_FROUDE)
        f_hz = f_lab / TIME_SCALE_FROUDE
        pxx = pxx_lab * TIME_SCALE_FROUDE
        nyquist = 0.5 / (dt_use * TIME_SCALE_FROUDE)
    else:
        f_hz, pxx = welch_psd(yu, dt_use, None if tiled else TSEG_PROTO_S)
        nyquist = 0.5 / dt_use
    m = (f_hz > 0.0) & (pxx > 0.0) & np.isfinite(pxx)
    if not np.any(m):
        return None
    f_hi = min(F_MAX_PROTO_HZ, 0.95 * nyquist)
    return f_hz[m], pxx[m], f_hi



def _load_os_disp_proto_mm(mat_name: str, key: str) -> tuple[np.ndarray, np.ndarray] | None:
    """
    Primary *OS displacement: lab Time (model s), Δu prototype mm.

    Args:    mat_name; key  tarSigOS | comSigOS | meaSigOS
    Returns: (t_lab, u_proto_mm) or None
    """
    npz = MAT_EXTRACT_DIR / f"{Path(mat_name).stem}.npz"
    if not npz.is_file():
        return None
    z = np.load(npz, allow_pickle=True)
    if f"{key}_time" in z.files and f"{key}_primary" in z.files:
        t = np.asarray(z[f"{key}_time"], dtype=float)
        u_m = np.asarray(z[f"{key}_primary"], dtype=float)
        return t, model_disp_to_proto_mm(u_m)
    block = load_block(z, key)
    if block is None:
        return None
    t, names, data = block
    cols = non_time_cols(names, data.shape[1])
    if not cols:
        return None
    return t, model_disp_to_proto_mm(np.asarray(data[:, cols[0]], dtype=float))


def write_force_psd(
    test_id: str, *, font_scale: float = DEFAULT_FONT_SCALE_PSD
) -> int:
    """
    Actuator force Welch PSD on $t_{\mathrm{int}}$ only (no map).

    Writes ``psd_frc_actuator.png``. Windows match structure/SSI PSDs.
    """
    apply_paper_style()
    scale_paper_fonts(font_scale)
    label_fs = float(plt.rcParams["axes.labelsize"])
    plt.rcParams["axes.titlesize"] = label_fs
    plt.rcParams["figure.titlesize"] = label_fs * 1.15
    pair = mat_dump_for_test(test_id)
    if pair is None:
        print(f"PlotHydroSpectra: skip {test_id} force PSD (no mat+dump)", file=sys.stderr)
        return 1
    _mat, dump = pair
    root = resolve_opensees_data() or LOCAL_OPENSEES_DATA
    dump_path = root / dump
    frc = load_daq_force_kn(dump_path)
    if frc is None:
        print(f"PlotHydroSpectra: skip {test_id} force PSD (no daqFrc)", file=sys.stderr)
        return 1
    t_f, f_kn = frc
    t_wave_period = wave_period_proto_s(test_id)
    t_wave_wed = (
        detect_wave_start_proto_s(t_f, f_kn, t_wave_period)
        if t_wave_period is not None
        else None
    )
    wave_hit = detect_wave_hit_proto_s(t_f, f_kn)
    t_wave_fri = None if wave_hit is None else float(wave_hit[0])
    t_wave_label = t_wave_wed if t_wave_wed is not None else t_wave_fri
    gm0 = gm_start_time_s(dump_path)
    d595 = d595_proto_window(gm0)
    windows = analysis_windows(
        float(t_f[-1]),
        d595,
        t_wave_wed=t_wave_wed,
        t_wave_fri=t_wave_fri,
    )
    if not windows:
        return 1
    f_wave = (1.0 / t_wave_period) if t_wave_period else None
    series = [
        (
            r"$-F_{\mathrm{daq}}$ (on numerical)",
            t_f,
            f_kn,
            r"(kN$^2$/Hz) proto",
            r"(kN$^2$/Hz) model",
            PSD_F_TO_MODEL,
            COLOR_F,
        )
    ]
    out = test_os_plots_dir(test_id) / OUT_NAME_FRC
    save_psd_figure(
        test_id=test_id,
        out_path=out,
        windows=windows,
        series=series,
        f_wave=f_wave,
        f_pier=pier_f1_hz(test_id),
        f_mode4=pier_f4_hz(test_id),
        t_wave_label=t_wave_label,
        font_scale=font_scale,
        title_extra=r"  ·  actuator force",
    )
    return 0


def write_lab_vs_pier_domain_psd(
    test_id: str, *, font_scale: float = DEFAULT_FONT_SCALE_PSD
) -> int:
    """
    PSD overlay: tar / com / mea on lab $t$ vs pier UX on $t_{\mathrm{int}}$.

    Writes ``psd_tar_com_mea_lab_vs_pier_os.png``. Freq axis = prototype Hz
    (lab PSD via /√λ; pier PSD on OpenSees t). Welch on each native grid
    only — no resample / interpolate.
    """
    apply_paper_style()
    scale_paper_fonts(font_scale)
    label_fs = float(plt.rcParams["axes.labelsize"])
    plt.rcParams["axes.titlesize"] = label_fs
    plt.rcParams["figure.titlesize"] = label_fs * 1.15
    pair = mat_dump_for_test(test_id)
    if pair is None:
        print(f"PlotHydroSpectra: skip {test_id} (no mat+dump)", file=sys.stderr)
        return 1
    mat, dump = pair
    root = resolve_opensees_data() or LOCAL_OPENSEES_DATA
    dump_path = root / dump

    tar = _load_os_disp_proto_mm(mat, "tarSigOS")
    com = _load_os_disp_proto_mm(mat, "comSigOS")
    mea = _load_os_disp_proto_mm(mat, "meaSigOS")
    pier = load_pier_ux_mm(dump_path)
    frc = load_daq_force_kn(dump_path)
    if mea is None or pier is None or frc is None:
        print(
            f"PlotHydroSpectra: skip {test_id} (need mea + pier + F)",
            file=sys.stderr,
        )
        return 1

    t_mea, u_mea = mea
    t_pier_os, u_pier = pier
    t_f_os, f_kn = frc

    t_wave_period = wave_period_proto_s(test_id)
    t_wave_wed_os = (
        detect_wave_start_proto_s(t_f_os, f_kn, t_wave_period)
        if t_wave_period is not None
        else None
    )
    wave_hit = detect_wave_hit_proto_s(t_f_os, f_kn)
    t_wave_fri_os = None if wave_hit is None else float(wave_hit[0])
    t_wave_wed_lab = (
        float(os_times_to_lab(t_wave_wed_os, mat, dump_dir=dump_path)[0])
        if t_wave_wed_os is not None
        else None
    )
    t_wave_fri_lab = (
        float(os_times_to_lab(t_wave_fri_os, mat, dump_dir=dump_path)[0])
        if t_wave_fri_os is not None
        else None
    )
    gm0 = gm_start_time_s(dump_path)
    d595_os = d595_proto_window(gm0)
    d595_lab = (
        os_window_to_lab(d595_os[0], d595_os[1], mat, dump_dir=dump_path)
        if d595_os is not None
        else None
    )
    windows_os = analysis_windows(
        float(t_pier_os[-1]),
        d595_os,
        t_wave_wed=t_wave_wed_os,
        t_wave_fri=t_wave_fri_os,
    )
    windows_lab = analysis_windows(
        float(t_mea[-1]),
        d595_lab,
        t_wave_wed=t_wave_wed_lab,
        t_wave_fri=t_wave_fri_lab,
    )
    n_win = min(len(windows_os), len(windows_lab))
    if n_win == 0:
        return 1

    lab_series: list[tuple[str, np.ndarray, np.ndarray, str]] = []
    if tar is not None:
        lab_series.append((r"tarSig (lab $t$)", tar[0], tar[1], COLOR_TAR))
    if com is not None:
        lab_series.append((r"comSig (lab $t$)", com[0], com[1], COLOR_COM))
    lab_series.append((r"meaSig (lab $t$)", t_mea, u_mea, COLOR_MEA))

    f_wave = (1.0 / t_wave_period) if t_wave_period else None
    f_pier = pier_f1_hz(test_id)
    f_mode4 = pier_f4_hz(test_id)

    fig_w = max(4.4 * n_win, 11.0)
    fig_h = 4.2 * (0.7 + 0.3 * font_scale)
    fig, axes = plt.subplots(
        1,
        n_win,
        figsize=(fig_w, fig_h),
        sharey=True,
        layout="constrained",
        squeeze=False,
    )
    fig.set_constrained_layout_pads(w_pad=0.02, h_pad=0.02, wspace=0.04)
    axes_1d = axes[0]

    for j in range(n_win):
        ax = axes_1d[j]
        wname, t0_os, t1_os = windows_os[j]
        _wn, t0_lab, t1_lab = windows_lab[j]
        f_hi = F_MAX_PROTO_HZ
        for label, t, u, color in lab_series:
            got = _psd_proto_hz(t, u, t0_lab, t1_lab, lab_clock=True, tiled=is_wed(test_id))
            if got is None:
                continue
            f_hz, pxx, f_hi_i = got
            ax.loglog(f_hz, pxx, color=color, lw=1.15, label=label, zorder=3)
            f_hi = min(f_hi, f_hi_i)
        got_p = _psd_proto_hz(t_pier_os, u_pier, t0_os, t1_os, lab_clock=False, tiled=is_wed(test_id))
        if got_p is not None:
            f_hz, pxx, f_hi_i = got_p
            ax.loglog(
                f_hz,
                pxx,
                color=COLOR_PIER_OS,
                lw=1.4,
                label=r"pier $u_x$ ($t_{\mathrm{int}}$)",
                zorder=2,
            )
            f_hi = min(f_hi, f_hi_i)
        mark_freq_guides(
            ax,
            f_wave=f_wave,
            f_pier=f_pier,
            f_mode4=f_mode4,
            f_max=f_hi * 1.05,
            with_labels=(j == n_win - 1),
        )
        ax.set_xlim(F_MIN_PROTO_HZ, f_hi)
        ax.grid(True, which="both", ls=":", alpha=0.45)
        ax.text(
            0.02,
            0.96,
            wname,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=float(plt.rcParams["xtick.labelsize"]),
        )
        add_dual_freq_xaxis(ax, label=(j == n_win // 2))
        if j == 0:
            ax.set_ylabel(r"$\Delta u$ PSD" "\n" r"(mm$^2$/Hz) proto")
        if j == n_win - 1:
            add_dual_psd_yaxis(ax, PSD_UX_TO_MODEL, r"(mm$^2$/Hz) model")
            ax.legend(
                loc="lower left",
                frameon=True,
                fontsize=float(plt.rcParams["legend.fontsize"]),
            )
        if j == n_win // 2:
            ax.set_xlabel(r"$f$ (Hz) prototype scale")

    mesh = soil_mesh_id(test_id)
    mesh_bit = f"  ·  mesh{mesh}" if mesh is not None else ""
    fig.suptitle(
        run_title(test_id)
        + mesh_bit
        + r"  ·  Welch PSD: tar/com/mea (lab $t$) vs pier ($t_{\mathrm{int}}$)"
    )
    out = test_os_plots_dir(test_id) / OUT_NAME_LAB_VS_PIER
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"PlotHydroSpectra: wrote {out}  windows={n_win}")
    return 0


def default_tests() -> list[str]:
    """Wed rows with waveT, plus Friday RTHS rows with a mat+dump pair."""
    out: list[str] = []
    for row in load_lab_runs_rows():
        tid = (row.get("Test") or "").strip()
        if not tid:
            continue
        if (row.get("waveT_model_s") or "").strip():
            out.append(tid)
            continue
        if tid.startswith("F") and mat_dump_for_test(tid) is not None:
            out.append(tid)
    return out


def main() -> int:
    raw = sys.argv[1:]
    if any(a in ("-h", "--help") for a in raw):
        print(__doc__)
        print(
            "  --lab-vs-pier-os     tar/com/mea (lab $t$) vs pier UX ($t_int$)"
        )
        return 0
    lab_pier = "--lab-vs-pier-os" in raw
    args = [a for a in raw if a not in ("-h", "--help", "--lab-vs-pier-os")]
    tests = args if args else default_tests()
    if not tests:
        print("PlotHydroSpectra: no tests", file=sys.stderr)
        return 1
    rc = 0
    for tid in tests:
        if lab_pier:
            rc = max(rc, write_lab_vs_pier_domain_psd(tid))
        else:
            rc = max(rc, write_plot(tid))
        rc = max(rc, write_force_psd(tid))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
