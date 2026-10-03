#!/usr/bin/env python3
"""
Three-window MKR–OpenFresco mockup: domain t and interface disp vs T_w.

Clocks from healthy lab proportions (F06 / STATEOS_SIGNALS):
  Δt_con = 1, Δt_sim = Δt_int = 10, α_f = 2/3, no slowdown.
  Host burst after force ≈ 4 Δt_con (~40% extrap): finish ≈ 3, begin ≈ 1.
  First begin ≈ 2 (no prior finish). Plot window: T_w / Δt_con ∈ [0, 30].

What matters for the hybrid handshake: force is acquired while domain t is
at the α-station. Disp panels show tar U_{n+α} (lab instruction / that
station) and comSig (extrap/interp). Committed U_{n+1} is not overlaid.

Targets U_{n+α} from a slow sine (peak at T_w=60).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "plot"))
from lab_paths import CYLINDER_LENGTH_SCALE  # noqa: E402

OUT = Path(__file__).with_name("mkr_disp_vs_tw")
DT = 10.0
AF = 2.0 / 3.0
FINISH = 3.0  # after atTarget: rd | sol | cmt (heavier)
BEGIN = 1.0  # pred | form | send
BURST = FINISH + BEGIN  # 4 → ~40% of window extrapolate
FIRST_BEGIN = 2.0
N_WIN = 3
TW_MAX = 30.0
LAMBDA = CYLINDER_LENGTH_SCALE

SINE_AMP = 2.0
SINE_T_PEAK = 60.0
SINE_OMEGA = 0.5 * np.pi / SINE_T_PEAK

COLOR_ACT = "#2a9d8f"
COLOR_INIT = "#9e9e9e"
COLOR_TAR = "#1f4e79"


def sine_u(t: float) -> float:
    return float(SINE_AMP * np.sin(SINE_OMEGA * t))


SEND = [FIRST_BEGIN] + [k * DT + BURST for k in range(1, N_WIN)]  # 2, 14, 24
FORCE_READY = [(k + 1) * DT for k in range(N_WIN)]  # 10, 20, 30
COMMIT = [fr + FINISH for fr in FORCE_READY]  # 13, 23, 33

U_ALPHA = np.array([sine_u(fr) for fr in FORCE_READY])


def reconstruct_U(u_alpha: np.ndarray, u0: float) -> np.ndarray:
    """U_0..U_N from U_{n+α}=(1-α)U_n+α U_{n+1}."""
    U = [float(u0)]
    for ua in u_alpha:
        U.append((float(ua) - (1.0 - AF) * U[-1]) / AF)
    return np.asarray(U)


U = reconstruct_U(U_ALPHA, sine_u(0.0))


def ualpha(n: int) -> float:
    return float(U_ALPHA[n])


def domain_events():
    """(T_w, domain t): α at send/hold through force+finish, n+1 at commit."""
    ev = [(0.0, 0.0)]
    ev.append((SEND[0], AF * DT))
    for k in range(N_WIN):
        ev.append((COMMIT[k], (k + 1) * DT))
        if k + 1 < N_WIN:
            ev.append((SEND[k + 1], (k + 1 + AF) * DT))
    return ev


def tar_events():
    """(T_w, tar U_{n+α}): steps at send; hold until next send."""
    ev = [(0.0, float(U[0]))]
    for k, ts in enumerate(SEND):
        ev.append((ts, ualpha(k)))
    return ev


def steps_from_events(tw: np.ndarray, events, col: int) -> np.ndarray:
    y = np.empty_like(tw)
    for i, x in enumerate(tw):
        val = events[0][col]
        for e in events:
            if x >= e[0]:
                val = e[col]
        y[i] = val
    return y


def lagrange(t: float, tp: np.ndarray, up: np.ndarray) -> float:
    tp = np.asarray(tp, dtype=float)
    up = np.asarray(up, dtype=float)
    n = tp.size
    if n == 1:
        return float(up[0])
    if n == 2:
        s = (t - tp[0]) / (tp[1] - tp[0])
        return float(up[0] + s * (up[1] - up[0]))
    tot = 0.0
    for i in range(n):
        num = 1.0
        den = 1.0
        for j in range(n):
            if i == j:
                continue
            num *= t - tp[j]
            den *= tp[i] - tp[j]
        tot += up[i] * num / den
    return float(tot)


def first_ramp(t: float) -> float:
    """Linear U_0 → U_{0+α}, same slope continued as extrap until SEND[1]."""
    v = (ualpha(0) - float(U[0])) / (FORCE_READY[0] - SEND[0])
    return float(U[0]) + v * (t - SEND[0])


def poly_after_send(k: int) -> tuple[np.ndarray, np.ndarray]:
    """
    2nd-order Lagrange after send of U_{(k+1)+α}: prior α-station, last
    extrapolated sample at send, new target at upcoming force ready.
    """
    t_send = SEND[k + 1]
    t_fr_next = FORCE_READY[k + 1]
    if k == 0:
        u_ex = first_ramp(t_send)
        t_prior, u_prior = FORCE_READY[0], ualpha(0)
    else:
        tp_prev, up_prev = poly_after_send(k - 1)
        u_ex = lagrange(t_send, tp_prev, up_prev)
        t_prior, u_prior = FORCE_READY[k], ualpha(k)
    tp = np.array([t_prior, t_send, t_fr_next])
    up = np.array([u_prior, u_ex, ualpha(k + 1)])
    return tp, up


def actuator(tw: float) -> float:
    if tw < SEND[0]:
        return float(U[0])
    if tw < SEND[1]:
        return first_ramp(tw)
    if tw < SEND[2]:
        tp, up = poly_after_send(0)
        return lagrange(tw, tp, up)
    tp, up = poly_after_send(1)
    return lagrange(tw, tp, up)


def extrap_ghost_segments():
    segs = []
    tt = np.linspace(SEND[1], FORCE_READY[1], 50)
    segs.append((tt, np.array([first_ramp(x) for x in tt])))
    tp, up = poly_after_send(0)
    tt = np.linspace(SEND[2], FORCE_READY[2], 50)
    segs.append((tt, np.array([lagrange(x, tp, up) for x in tt])))
    return segs


def comsig_phases():
    phases = [(0.0, SEND[0], "init")]
    phases.append((SEND[0], FORCE_READY[0], "interp"))
    for k in range(N_WIN - 1):
        t_fr = FORCE_READY[k]
        t_send = SEND[k + 1]
        phases.append((t_fr, t_send, "extrap"))
        phases.append((t_send, FORCE_READY[k + 1], "interp"))
    return phases


def phase_at(t: float) -> str:
    for t0, t1, kind in comsig_phases():
        if t0 - 1e-12 <= t <= t1 + 1e-12:
            return kind
    return "init"


def plot_comsig(ax, tw, u, *, scale=1.0, zorder=4, legend=True):
    (ln,) = ax.plot(
        tw,
        u * scale,
        color=COLOR_ACT,
        lw=1.85,
        zorder=zorder,
        label=r"comSig" if legend else None,
    )
    t_m = np.arange(0.0, TW_MAX + 0.5, 1.0)
    u_m = np.array([actuator(float(t)) for t in t_m]) * scale
    kinds = np.array([phase_at(float(t)) for t in t_m])
    style = {
        "init": dict(
            mfc="white",
            mec=COLOR_INIT,
            label="initialize" if legend else None,
        ),
        "extrap": dict(
            mfc="white",
            mec=COLOR_ACT,
            label="extrapolate" if legend else None,
        ),
        "interp": dict(
            mfc=COLOR_ACT,
            mec=COLOR_ACT,
            label="interpolate" if legend else None,
        ),
    }
    handles = [ln]
    for kind in ("init", "extrap", "interp"):
        m = kinds == kind
        if not np.any(m):
            continue
        (mk,) = ax.plot(
            t_m[m],
            u_m[m],
            "o",
            ms=5.0,
            mew=1.15,
            zorder=zorder + 1,
            **style[kind],
        )
        if legend:
            handles.append(mk)
    return handles


def plot_tar_step(ax, tw, u_tar, *, scale=1.0, zorder=5, legend=True):
    (ln,) = ax.plot(
        tw,
        u_tar * scale,
        color=COLOR_TAR,
        lw=2.0,
        drawstyle="steps-post",
        zorder=zorder,
        label=r"tar $U_{n+\alpha}$ (force at this station)" if legend else None,
    )
    t_send = np.asarray(SEND, dtype=float)
    u_send = np.array([ualpha(k) for k in range(len(SEND))]) * scale
    (mk,) = ax.plot(
        t_send,
        u_send,
        "o",
        ms=6.0,
        mfc="white",
        mec=COLOR_TAR,
        mew=1.4,
        zorder=zorder + 1,
        label=r"send" if legend else None,
    )
    return [ln, mk]


def main() -> None:
    tw = np.linspace(0.0, TW_MAX, 3000)
    t_dom = steps_from_events(tw, domain_events(), 1)
    u_tar = steps_from_events(tw, tar_events(), 1)
    u_act = np.array([actuator(x) for x in tw])

    t_mid = 0.5 * (SEND[2] + FORCE_READY[2])
    tp0, up0 = poly_after_send(0)
    u_ex_ghost = lagrange(t_mid, tp0, up0)
    u_in = actuator(t_mid)
    print(
        f"SEND={SEND} FORCE_READY={FORCE_READY} COMMIT={COMMIT}\n"
        f"U_alpha={U_ALPHA}  U_commit={U}\n"
        f"at t={t_mid:.1f}: extrap_ghost={u_ex_ghost:.3f}  interp={u_in:.3f}  "
        f"sine={sine_u(t_mid):.3f}  gap={u_ex_ghost - u_in:.3f}"
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
    fig, (ax0, ax1, ax2) = plt.subplots(
        3,
        1,
        figsize=(8.2, 7.2),
        sharex=False,
        gridspec_kw={"height_ratios": [1.0, 1.35, 1.2]},
    )

    ax0.plot(
        tw, t_dom, color="#3d5a73", lw=2.0, drawstyle="steps-post", label=r"domain $t$"
    )
    for fr in FORCE_READY:
        ax0.axvline(fr, color="#555", ls="--", lw=0.8, alpha=0.65)
    ax0.set_ylabel(r"domain $t$ / $\Delta t_{\mathrm{con}}$")
    ax0.set_ylim(-1.5, TW_MAX + AF * DT * 0.15)
    ax0.legend(loc="upper left", frameon=False)
    ax0.set_title(
        r"Three-window mockup (domain $t$ at $\alpha$ when force is read; "
        r"tar / comSig; finish${=}3$ + begin${=}1$)",
        loc="left",
        fontsize=10.0,
    )

    for i, (tt, uu) in enumerate(extrap_ghost_segments()):
        ax1.plot(
            tt,
            uu,
            color=COLOR_ACT,
            lw=1.0,
            ls="--",
            alpha=0.45,
            label=r"extrap (no new target)" if i == 0 else None,
            zorder=2,
        )

    plot_comsig(ax1, tw, u_act, zorder=3)
    plot_tar_step(ax1, tw, u_tar, zorder=5)
    for fr in FORCE_READY:
        ax1.axvline(fr, color="#555", ls="--", lw=0.8, alpha=0.45)

    ax1.set_xlabel(r"$T_w$ / $\Delta t_{\mathrm{con}}$  (model)")
    ax1.set_ylabel("interface displacement (arb. model)")
    ax1.legend(loc="upper left", frameon=False, fontsize=8.0)
    ax1.set_xlim(0, TW_MAX)
    ax0.set_xlim(0, TW_MAX)
    y_hi = max(float(np.max(u_act)), float(np.max(u_tar)), float(SINE_AMP)) + 0.35
    for _, uu in extrap_ghost_segments():
        y_hi = max(y_hi, float(np.max(uu)) + 0.15)
    ax1.set_ylim(-0.25, y_hi)

    # bottom: same tar / comSig on T_w, prototype scale
    plot_comsig(ax2, tw, u_act, scale=LAMBDA, zorder=3)
    plot_tar_step(ax2, tw, u_tar, scale=LAMBDA, zorder=5)
    ax2.set_xlabel(r"$T_w$ / $\Delta t_{\mathrm{con}}$")
    ax2.set_ylabel(rf"$\Delta u$ prototype (arb.$\times\lambda$, $\lambda={LAMBDA:g}$)")
    ax2.legend(loc="upper left", frameon=False, fontsize=8.0)
    ax2.set_xlim(0, TW_MAX)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT.with_suffix(f".{ext}"))
        print("wrote", OUT.with_suffix(f".{ext}"))


if __name__ == "__main__":
    main()
