#!/usr/bin/env python3
"""
Show the lab clock zero mismatch at session start (F05).

Lab Time ticks from 0 while the host finishes startup. First post-target
extrapolate entry is one Δt_con past end-of-window; that window end is the
lab's first completed Δt_sim (no α_f in the offset):
  offset = t_ex - Δt_con - Δt_sim.

Bottom: tar U_{n+α} on t_k = offset + (k+α_f) Δt_sim vs comSig on lab Time.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "plot"))
from lab_paths import MAT_EXTRACT_DIR, M_TO_MM  # noqa: E402

OUT = Path(__file__).with_name("clock_offset_start")
MAT_STEM = "0821_GusBridge_rowNeg4"
TEST_ID = "F05"
N_COUNT = 10
ALPHA_F = 2.0 / 3.0

COLOR_OS = "#1f4e79"
COLOR_COM = "#2a9d8f"
COLOR_MARK = "#c45c26"
COLOR_GOAL = "#2a9d8f"


def first_exceed(t: np.ndarray, u: np.ndarray, thr: float, u0: float = 0.0) -> tuple[float, float]:
    """First (t, u) with |u-u0| > thr."""
    i = int(np.argmax(np.abs(u - u0) > thr))
    if not (np.abs(u[i] - u0) > thr):
        raise RuntimeError(f"never exceeds {thr}")
    return float(t[i]), float(u[i])


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


def first_extrap_after_target(t: np.ndarray, tc3: np.ndarray, tar: np.ndarray) -> float:
    """Lab time of first typeConv3 0→1 after the first real target land."""
    t_land_lab, _ = target_lands(t, tc3, tar)
    t_land0 = float(t_land_lab[0])
    prev = np.concatenate([[int(tc3[0])], tc3[:-1]])
    enter = np.where((tc3 == 1) & (prev == 0) & (t > t_land0))[0]
    if enter.size == 0:
        raise RuntimeError("no 0→1 extrapolate after first real target")
    return float(t[int(enter[0])])


def main() -> None:
    z = np.load(MAT_EXTRACT_DIR / f"{MAT_STEM}.npz", allow_pickle=True)
    t = np.asarray(z["comSigOS_time"], dtype=float)
    tc3 = np.asarray(z["stateOS_data"][:, 0], dtype=int)
    com = np.asarray(z["comSigOS_data"][:, 0], dtype=float)
    tar = np.asarray(z["tarSigOS_data"][:, 0], dtype=float)

    dt = float(np.median(np.diff(t)))
    dt_sim = N_COUNT * dt
    t_land_lab, u_land = target_lands(t, tc3, tar)
    t_land = float(t_land_lab[0])
    t_ex = first_extrap_after_target(t, tc3, tar)
    t_goal = t_ex - dt
    offset = t_goal - dt_sim

    k = np.arange(u_land.size, dtype=float)
    t_tar = offset + (k + ALPHA_F) * dt_sim

    tar0 = float(tar[0])
    com0 = float(com[0])
    com_mm = (com - com0) * M_TO_MM
    tar_mm = (u_land - tar0) * M_TO_MM

    thr = 1.0
    t_com1, _ = first_exceed(t, com_mm, thr)
    t_tar1, _ = first_exceed(t_tar, tar_mm, thr)
    lag_1mm_ms = (t_com1 - t_tar1) * 1e3
    slow_ms = float(np.sum(tc3[t <= 2.0] == 2) * dt * 1e3)

    print(
        f"{TEST_ID}: t_land={t_land*1e3:.1f} ms; t_ex={t_ex*1e3:.1f} ms; "
        f"t_goal={t_goal*1e3:.1f} ms; offset={offset*1e3:.1f} ms "
        f"(t_ex - dt_con - dt_sim); "
        f"1 mm lag com-tar={lag_1mm_ms:.1f} ms; "
        f"slow samples in first 2 s={slow_ms:.1f} ms"
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
    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=(8.2, 6.4), gridspec_kw={"height_ratios": [1.15, 1.2]}
    )

    # --- top: lab state + first window complete ---
    t_hi = 0.35
    m = t <= t_hi
    ax0.plot(t[m], tc3[m], color="#3d5a73", lw=1.6, drawstyle="steps-post", label="typeConv3")
    ax0.axvline(
        t_land,
        color="#888",
        lw=1.2,
        ls="--",
        label=rf"first target land (${t_land*1e3:.0f}\,\mathrm{{ms}}$)",
    )
    ax0.axvline(
        t_goal,
        color=COLOR_GOAL,
        lw=1.6,
        label=rf"window end $t_{{\mathrm{{ex}}}}-\Delta t_{{\mathrm{{con}}}}$ (${t_goal*1e3:.1f}\,\mathrm{{ms}}$)",
    )
    ax0.axvline(
        t_ex,
        color=COLOR_MARK,
        lw=1.4,
        label=rf"first $0\rightarrow 1$ ($t_{{\mathrm{{ex}}}}={t_ex*1e3:.1f}\,\mathrm{{ms}}$)",
    )
    # first tar on OS grid
    t_tar0 = float(t_tar[0])
    ax0.axvline(
        t_tar0,
        color=COLOR_OS,
        lw=1.3,
        ls=":",
        label=rf"first $U_{{n+\alpha}}$ on OS grid (${t_tar0*1e3:.1f}\,\mathrm{{ms}}$)",
    )
    ax0.axvspan(0.0, t_goal, color=COLOR_MARK, alpha=0.10, label="before first lab $\\Delta t_{\\mathrm{sim}}$")
    ax0.set_yticks([-1, 0, 1, 2])
    ax0.set_yticklabels(["init", "interp", "extrap", "slow"])
    ax0.set_xlim(0, t_hi)
    ax0.set_ylim(-1.6, 2.6)
    ax0.set_ylabel(r"typeConv3")
    ax0.legend(loc="upper right", frameon=False, fontsize=6.5)
    ax0.set_title(
        rf"{TEST_ID} — offset $=t_{{\mathrm{{ex}}}}-\Delta t_{{\mathrm{{con}}}}-\Delta t_{{\mathrm{{sim}}}}$ "
        rf"$={offset*1e3:.1f}\,\mathrm{{ms}}$ (no $\alpha_f$ in offset; lab Time unchanged)",
        loc="left",
        fontsize=9.0,
    )

    # --- bottom: tar on OS grid vs com on lab time ---
    pad = 0.35
    t0z = min(t_tar1, t_com1) - pad
    t1z = max(t_tar1, t_com1) + pad
    m1 = (t >= t0z) & (t <= t1z)
    mt1 = (t_tar >= t0z) & (t_tar <= t1z)
    ax1.plot(
        t_tar[mt1],
        tar_mm[mt1],
        color=COLOR_OS,
        lw=2.0,
        marker="x",
        ms=6.0,
        markeredgewidth=1.3,
        label=rf"$U_{{n+\alpha}}$ (tar @ offset$+(k+\alpha_f)\Delta t_{{\mathrm{{sim}}}}$)",
    )
    ax1.plot(t[m1], com_mm[m1], color=COLOR_COM, lw=1.6, label="comSig (lab $t$)")
    ax1.axhline(thr, color="#555", ls=":", lw=0.9)
    ax1.axhline(-thr, color="#555", ls=":", lw=0.9)
    ax1.axvline(t_tar1, color=COLOR_OS, ls="--", lw=1.2)
    ax1.axvline(t_com1, color=COLOR_COM, ls="--", lw=1.2)
    ymin = min(float(np.min(tar_mm[mt1])), float(np.min(com_mm[m1])), -thr) - 0.4
    ymax = max(float(np.max(tar_mm[mt1])), float(np.max(com_mm[m1])), thr) + 0.9
    ax1.set_ylim(ymin, ymax)
    y_ann = thr + 0.25 * (ymax - thr)
    ax1.annotate(
        "",
        xy=(t_com1, y_ann),
        xytext=(t_tar1, y_ann),
        arrowprops=dict(arrowstyle="<->", color=COLOR_MARK, lw=1.5),
    )
    ax1.text(
        0.5 * (t_tar1 + t_com1),
        y_ann + 0.04 * (ymax - ymin),
        rf"com$-$tar at first 1 mm $\approx{lag_1mm_ms:.0f}\,\mathrm{{ms}}$"
        rf"  [$(1-\alpha_f)\Delta t_{{\mathrm{{sim}}}}\approx{(1-ALPHA_F)*dt_sim*1e3:.1f}$]",
        ha="center",
        va="bottom",
        color=COLOR_MARK,
        fontsize=8.5,
    )
    ax1.set_xlim(t0z, t1z)
    ax1.set_xlabel(r"$t$ lab (s, model)")
    ax1.set_ylabel(r"$\Delta$disp (mm, model)")
    ax1.legend(loc="upper left", frameon=False, fontsize=7.5)
    ax1.set_title(
        rf"offset ${offset*1e3:.1f}\,\mathrm{{ms}}$; tar on OS send grid "
        rf"(slow in first 2 s: ${slow_ms:.0f}\,\mathrm{{ms}}$)",
        loc="left",
        fontsize=9.5,
    )

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT.with_suffix(f".{ext}"))
        print("wrote", OUT.with_suffix(f".{ext}"))


if __name__ == "__main__":
    main()
