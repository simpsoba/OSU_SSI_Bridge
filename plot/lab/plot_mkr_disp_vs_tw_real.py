#!/usr/bin/env python3
"""
Real-data companion to plot_mkr_disp_vs_tw.py (mockup kept separate).

comSig + meaSig on lab Time. Host targets (tarSig = U_{n+α}) on the
MKR send grid mapped to lab:
  t_k = offset + (k + α_f) Δt_sim
  offset = t_ex - Δt_con - Δt_sim
  (first post-target 0→1 is one Δt_con past window end; that end is the
  lab's first completed Δt_sim — α_f only enters the send grid, not offset).

Top: a few response cycles around the detail window. Bottom: six-window
zoom (typeConv3 + comSig phase markers) of the shaded band.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "plot"))
from lab_paths import MAT_EXTRACT_DIR, M_TO_MM  # noqa: E402

OUT = Path(__file__).with_name("mkr_disp_vs_tw_real")
# F05: early Fri CudaMKR np=8 baseline mesh; only startup typeConv3==2
MAT_STEM = "0821_GusBridge_rowNeg4"
TEST_ID = "F05"
N_WIN = 6
N_COUNT = 10
N_SAMP = N_WIN * N_COUNT
# MKR ρ∞=0.5 → α_f = 2/3 (target at n+α)
ALPHA_F = 2.0 / 3.0
# top panel: ~4 cycles (T≈1.7 s model near this zoom)
N_CYCLES_TOP = 4
T_CYCLE_MODEL_S = 1.7
# ignore OpenFresco startup spikes before this when checking "no prior slowdown"
STARTUP_S = 1.0

COLOR_ACT = "#2a9d8f"
COLOR_MEA = "#6a1b9a"
COLOR_INIT = "#9e9e9e"
COLOR_OS = "#1f4e79"
COLOR_SLOW = "#c45c26"
COLOR_ZOOM = "#f4a261"

PHASE = {-1: "init", 0: "interp", 1: "extrap", 2: "slowdown"}


def slowdown_onset_times(t: np.ndarray, tc3: np.ndarray) -> np.ndarray:
    """Lab times (model s) of typeConv3 rising edges into 2."""
    prev = np.concatenate([[int(tc3[0])], tc3[:-1]])
    return t[(tc3 == 2) & (prev != 2)]


def target_lands(
    t: np.ndarray, tc3: np.ndarray, tar: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Lab times and tar values at typeConv3 1→0 with non-tiny tar."""
    land = np.concatenate([[False], (tc3[:-1] == 1) & (tc3[1:] == 0)])
    idx = [i for i in np.where(land)[0] if abs(float(tar[i])) > 1e-9]
    if not idx:
        raise RuntimeError("no real target lands")
    ii = np.asarray(idx, dtype=int)
    return t[ii], tar[ii]


def first_extrap_after_target(
    t: np.ndarray, tc3: np.ndarray, tar: np.ndarray
) -> float:
    """Lab time of first typeConv3 0→1 after the first real target land."""
    t_land_lab, _ = target_lands(t, tc3, tar)
    t_land0 = float(t_land_lab[0])
    prev = np.concatenate([[int(tc3[0])], tc3[:-1]])
    enter = np.where((tc3 == 1) & (prev == 0) & (t > t_land0))[0]
    if enter.size == 0:
        raise RuntimeError("no 0→1 extrapolate after first real target")
    return float(t[int(enter[0])])


def lab_clock_offset(t: np.ndarray, tc3: np.ndarray, tar: np.ndarray, dt_con: float, dt_sim: float) -> tuple[float, float, float]:
    """
    Offset from first completed lab Δt_sim (no α_f).

    t_ex = first 0→1 after first real target; window end = t_ex - Δt_con.
    offset = t_ex - Δt_con - Δt_sim.

    Returns: offset, t_goal, t_ex (model s).
    """
    t_ex = first_extrap_after_target(t, tc3, tar)
    t_goal = t_ex - dt_con
    return t_goal - dt_sim, t_goal, t_ex


def tar_on_os_grid(
    u_tar: np.ndarray, offset: float, dt_sim: float
) -> np.ndarray:
    """Place U_{n+α} on t_k = offset + (k+α_f) Δt_sim."""
    k = np.arange(u_tar.size, dtype=float)
    return offset + (k + ALPHA_F) * dt_sim


def pick_zoom_start(
    t: np.ndarray,
    tc1: np.ndarray,
    tc3: np.ndarray,
    com: np.ndarray,
) -> tuple[int, float, float]:
    """
    Index i0 for a healthy 6-window zoom in strong motion.

    Requires zero slowdown onsets after STARTUP_S and through the zoom end
    (startup typeConv3==2 spikes under 1 s are ignored).

    Returns: i0, frac_ex, amp_m
    """
    ons = slowdown_onset_times(t, tc3)
    prev = np.concatenate([[np.nan], tc1[:-1]])
    starts = np.where((tc1 == 1) & (prev == 10) & (t >= 55.0) & (t <= 80.0))[0]
    best = None
    for i0 in starts:
        i1 = i0 + N_SAMP
        if i1 > len(t):
            continue
        t_end = float(t[i1 - 1])
        if np.any((ons > STARTUP_S) & (ons <= t_end)):
            continue
        st = tc3[i0:i1]
        frac_ex = float((st == 1).mean())
        frac_sl = float((st == 2).mean())
        amp = float(com[i0:i1].max() - com[i0:i1].min())
        if frac_sl > 0.0 or not (0.30 <= frac_ex <= 0.50):
            continue
        score = (amp, -abs(frac_ex - 0.40))
        if best is None or score > best[0]:
            best = (score, i0, frac_ex, amp)
    if best is None:
        raise RuntimeError(
            f"no 6-window zoom on {TEST_ID} with zero post-startup slowdowns"
        )
    _, i0, frac_ex, amp = best
    return i0, frac_ex, amp


def plot_com_markers(ax, tw, u, tc3, *, zorder=4):
    """Unicolor comSig + phase markers (open extrap / filled interp)."""
    ax.plot(tw, u, color=COLOR_ACT, lw=1.85, zorder=zorder, label="comSig")
    kinds = np.array([PHASE.get(int(v), "interp") for v in tc3])
    style = {
        "init": dict(mfc="white", mec=COLOR_INIT, label="initialize"),
        "extrap": dict(mfc="white", mec=COLOR_ACT, label="extrapolate"),
        "interp": dict(mfc=COLOR_ACT, mec=COLOR_ACT, label="interpolate"),
        "slowdown": dict(mfc=COLOR_SLOW, mec=COLOR_SLOW, label="slowdown"),
    }
    for kind in ("init", "extrap", "interp", "slowdown"):
        m = kinds == kind
        if not np.any(m):
            continue
        ax.plot(tw[m], u[m], "o", ms=4.5, mew=1.1, zorder=zorder + 1, **style[kind])


def main() -> None:
    npz = MAT_EXTRACT_DIR / f"{MAT_STEM}.npz"
    z = np.load(npz, allow_pickle=True)
    t = np.asarray(z["stateOS_time"], dtype=float)
    tc3 = np.asarray(z["stateOS_data"][:, 0], dtype=int)
    tc1 = np.asarray(z["stateOS_data"][:, 1], dtype=float)
    com = np.asarray(z["comSigOS_data"][:, 0], dtype=float)
    mea = np.asarray(z["meaSigOS_data"][:, 0], dtype=float)
    tar = np.asarray(z["tarSigOS_data"][:, 0], dtype=float)
    dt = float(np.median(np.diff(t)))
    dt_sim = N_COUNT * dt

    t_land_lab, u_land = target_lands(t, tc3, tar)
    offset, t_goal, t_ex = lab_clock_offset(t, tc3, tar, dt, dt_sim)
    t_tar = tar_on_os_grid(u_land, offset, dt_sim)
    u_tar = u_land

    i0, frac_ex, amp = pick_zoom_start(t, tc1, tc3, com)
    i1 = i0 + N_SAMP
    t_zoom0 = float(t[i0])
    t_zoom1 = float(t[i1 - 1])
    ons = slowdown_onset_times(t, tc3)
    ons_prior = ons[ons <= t_zoom1]
    ons_post_startup = ons_prior[ons_prior > STARTUP_S]

    # global zero = first sample of each series
    com0 = float(com[0])
    mea0 = float(mea[0])
    tar0 = float(tar[0])
    com_mm = (com - com0) * M_TO_MM
    mea_mm = (mea - mea0) * M_TO_MM
    tar_mm = (u_tar - tar0) * M_TO_MM

    # pad ±Δt_sim so tar / com / mea polylines all span the zoom frame
    m_lab_z = (t >= t_zoom0 - dt_sim) & (t <= t_zoom1 + dt_sim)
    tw = (t[m_lab_z] - t_zoom0) / dt
    tc3_z = tc3[m_lab_z]
    com_z = com_mm[m_lab_z]
    mea_z = mea_mm[m_lab_z]
    m_tar_z = (t_tar >= t_zoom0 - dt_sim) & (t_tar <= t_zoom1 + dt_sim)
    tw_tar = (t_tar[m_tar_z] - t_zoom0) / dt
    tar_z = tar_mm[m_tar_z]
    m_tar_mark = (t_tar >= t_zoom0) & (t_tar <= t_zoom1)
    tw_tar_mark = (t_tar[m_tar_mark] - t_zoom0) / dt
    tar_mark = tar_mm[m_tar_mark]

    print(
        f"{TEST_ID} {MAT_STEM}: full t=[{t[0]:.3f},{t[-1]:.3f}] s model; "
        f"tar on OS grid, offset={offset*1e3:.1f} ms "
        f"(t_ex={t_ex*1e3:.1f} - dt_con - dt_sim; t_goal={t_goal*1e3:.1f}; "
        f"t_land0={t_land_lab[0]*1e3:.1f}); "
        f"n_targets={t_tar.size}; "
        f"zoom t=[{t_zoom0:.6f},{t_zoom1:.6f}] ({N_WIN} windows, "
        f"extrap={frac_ex:.0%}, amp={amp*M_TO_MM:.3f} mm), "
        f"tar zoom samples={tw_tar.size}; "
        f"slowdown onsets <= zoom: {list(map(float, ons_prior))} "
        f"(post-startup before zoom: {len(ons_post_startup)})"
    )

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 11,
            "figure.dpi": 150,
            "savefig.dpi": 200,
            "axes.grid": True,
            "grid.alpha": 0.35,
        }
    )
    fig, (ax_full, ax_tc, ax_zoom) = plt.subplots(
        3,
        1,
        figsize=(8.2, 8.0),
        gridspec_kw={"height_ratios": [1.35, 0.65, 1.45]},
    )

    half = 0.5 * N_CYCLES_TOP * T_CYCLE_MODEL_S
    t_mid = 0.5 * (t_zoom0 + t_zoom1)
    t_lo = max(float(t[0]), t_mid - half)
    t_hi = min(float(t[-1]), t_mid + half)
    m_tar = (t_tar >= t_lo - dt_sim) & (t_tar <= t_hi + dt_sim)
    m_com = (t >= t_lo - dt_sim) & (t <= t_hi + dt_sim)

    ax_full.plot(
        t_tar[m_tar],
        tar_mm[m_tar],
        color=COLOR_OS,
        lw=1.35,
        label=rf"$U_{{n+\alpha}}$ (tar @ $t_{{\mathrm{{OS}}}}+\mathrm{{offset}}$, "
        rf"offset$={offset*1e3:.1f}\,\mathrm{{ms}}$)",
        zorder=3,
    )
    ax_full.plot(
        t[m_com],
        com_mm[m_com],
        color=COLOR_ACT,
        lw=1.15,
        alpha=0.9,
        label="comSig",
        zorder=2,
    )
    ax_full.plot(
        t[m_com],
        mea_mm[m_com],
        color=COLOR_MEA,
        lw=1.1,
        alpha=0.9,
        label="meaSig",
        zorder=2,
    )
    y0, y1 = ax_full.get_ylim()
    ax_full.add_patch(
        Rectangle(
            (t_zoom0, y0),
            t_zoom1 - t_zoom0,
            y1 - y0,
            facecolor=COLOR_ZOOM,
            edgecolor="none",
            alpha=0.35,
            zorder=1,
            label="zoom window",
        )
    )
    ax_full.set_ylim(y0, y1)
    ax_full.set_xlim(t_lo, t_hi)
    ax_full.set_ylabel(r"$\Delta$disp (mm, model)")
    ax_full.legend(loc="upper right", frameon=False, fontsize=7.5)
    ax_full.set_title(
        rf"{TEST_ID} — tar $U_{{n+\alpha}}$ on "
        rf"$t=\mathrm{{offset}}+(k+\alpha_f)\Delta t_{{\mathrm{{sim}}}}$; "
        rf"offset$=t_{{\mathrm{{ex}}}}-\Delta t_{{\mathrm{{con}}}}-\Delta t_{{\mathrm{{sim}}}}$ "
        rf"$={offset*1e3:.1f}\,\mathrm{{ms}}$ / com / mea; "
        rf"orange = {N_WIN}-window detail",
        loc="left",
        fontsize=8.5,
    )
    ax_full.set_xlabel(r"$t$ lab (s, model)")

    ax_tc.plot(
        tw,
        tc3_z,
        color="#3d5a73",
        lw=1.8,
        drawstyle="steps-post",
        label=r"typeConv3",
    )
    ax_tc.set_yticks([-1, 0, 1, 2])
    ax_tc.set_yticklabels(["init", "interp", "extrap", "slow"])
    ax_tc.set_ylabel(r"typeConv3")
    ax_tc.set_ylim(-1.6, 2.6)
    ax_tc.legend(loc="upper right", frameon=False, fontsize=8.0)
    ax_tc.set_title(
        rf"Zoom: $t_{{\mathrm{{lab}}}}={t_zoom0:.3f}$–${t_zoom1:.3f}\,\mathrm{{s}}$ "
        rf"(extrap ${frac_ex*100:.0f}\%$)",
        loc="left",
        fontsize=10.0,
    )

    for k in range(N_WIN + 1):
        for ax in (ax_tc, ax_zoom):
            ax.axvline(k * N_COUNT, color="#555", ls="--", lw=0.75, alpha=0.45)

    ax_zoom.plot(
        tw_tar,
        tar_z,
        color=COLOR_OS,
        lw=2.0,
        label=rf"$U_{{n+\alpha}}$ (tar, offset$={offset*1e3:.1f}\,\mathrm{{ms}}$)",
        zorder=5,
    )
    ax_zoom.plot(
        tw_tar_mark,
        tar_mark,
        color=COLOR_OS,
        linestyle="None",
        marker="x",
        ms=8,
        markeredgewidth=1.5,
        zorder=6,
    )
    ax_zoom.plot(
        tw,
        mea_z,
        color=COLOR_MEA,
        lw=1.5,
        label="meaSig",
        zorder=3,
    )
    plot_com_markers(ax_zoom, tw, com_z, tc3_z, zorder=4)
    ax_zoom.set_xlabel(r"$T_w$ / $\Delta t_{\mathrm{con}}$  (from zoom start)")
    ax_zoom.set_ylabel(r"$\Delta$disp (mm, model)")
    ax_zoom.legend(loc="best", frameon=False, fontsize=8.0)
    ax_zoom.set_xlim(0, N_SAMP)
    ax_tc.set_xlim(0, N_SAMP)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT.with_suffix(f".{ext}"))
        print("wrote", OUT.with_suffix(f".{ext}"))


if __name__ == "__main__":
    main()
