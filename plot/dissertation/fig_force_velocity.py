"""Velocity-dependent part of the physical-subassembly force after removing the inertial term (dissertation Ch. 6
draft; read-only on data).

Same data and kinematics as fig_force_kinematics.py (W05 and F05, 5-95 % Arias window, 20 Hz zero-phase low-pass on
the displacement at model scale, central differences, prototype scale, resisting-force convention).
Step 1: inertial term only, F = m a + F0 (one-parameter regression on the acceleration; a and v are nearly
uncorrelated for oscillatory motion, so m is not traded against a damping term). The moving mass of Seki et al. (2026),
325 kg at model scale (cylinder, water inside it, and actuator piston; 2.74 m water depth), is 4.49 t at prototype
scale and is marked for reference.
Step 2: residual R = F - m a - F0 against v, with the mean of R in velocity bins (no functional form assumed) and its
interquartile range; the lines are least-squares fits of R = c v, R = c_D v|v|, and R = F_c sign(v) (Coulomb friction,
e.g., pin, guide rails, actuator) to the binned means; rms misfit of each printed.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SCRIPT = Path(__file__).parent / "fig_force_kinematics.py"
OUT = Path(sys.argv[1])
src = SCRIPT.read_text(encoding="utf-8")
ns = {"__file__": str(SCRIPT)}
exec(compile(src[:src.index("prs.configure_font()")], str(SCRIPT), "exec"), ns)
R = ns["R"]
M_SEKI = 325.0 * R.LAM ** 3 / 1e3  # t, prototype
EDGES = np.linspace(-0.30, 0.30, 25)  # m/s, prototype

ns["prs"].configure_font()
fig, A = plt.subplots(2, 2, figsize=(6.0, 4.6), layout="constrained", gridspec_kw=dict(width_ratios=[1.0, 1.35]))
for row, (tid, stem, t_eq, name) in enumerate(ns["CASES"]):
    d = ns["load"](stem, t_eq)
    a, v, F = d["a"], d["v"] / 1e3, d["F"]
    m, f0 = np.linalg.lstsq(np.c_[a, np.ones(a.size)], F, rcond=None)[0]
    res = F - m * a - f0
    print(f"{tid}: m = {m:.2f} t (Seki moving mass {M_SEKI:.2f} t); corr(a, v) = {np.corrcoef(a, v)[0, 1]:+.3f}; "
          f"residual std {res.std():.2f} kN (force std {F.std():.2f} kN)")
    ax = A[row, 0]
    ax.plot(a, F, color="0.6", lw=0.15, alpha=0.6, rasterized=True)
    xx = np.array([a.min(), a.max()])
    ax.plot(xx, m * xx + f0, color="#D55E00", lw=1.6, label=f"fit, $m$ = {m:.2f} t")
    ax.plot(xx, M_SEKI * xx + f0, color="#D55E00", lw=0.9, ls=(0, (3, 2)), label=f"Seki et al., {M_SEKI:.2f} t")
    ax.legend(loc="upper left", frameon=False, fontsize=6.3, handlelength=1.6)
    ax.set_xlabel(r"$\ddot u$ (m/s$^2$)")
    ax.set_ylabel(r"$\lambda_L^3 f_p$ (kN)")
    ax.set_title(rf"\textbf{{({'ac'[row]})}} {name}", loc="left")

    bx = A[row, 1]
    bx.plot(v * 1e3, res, color="0.7", lw=0.12, alpha=0.5, rasterized=True)
    idx = np.digitize(v, EDGES)
    vc, rm, q1, q3 = [], [], [], []
    for i in range(1, EDGES.size):
        s = idx == i
        if s.sum() < 2000:
            continue
        vc.append(v[s].mean())
        rm.append(res[s].mean())
        q1.append(np.percentile(res[s], 25))
        q3.append(np.percentile(res[s], 75))
    vc, rm = np.array(vc), np.array(rm)
    bx.fill_between(vc * 1e3, q1, q3, color="#009E73", alpha=0.18, lw=0)
    bx.plot(vc * 1e3, rm, "o-", color="#009E73", ms=3, lw=0.9, label="binned mean (IQR shaded)")
    c = np.sum(vc * rm) / np.sum(vc * vc)
    q = vc * np.abs(vc)
    cd = np.sum(q * rm) / np.sum(q * q)
    vv = np.linspace(vc.min(), vc.max(), 50)
    bx.plot(vv * 1e3, c * vv, color="k", lw=0.9, label=f"$c\\,\\dot u$, $c$ = {c:.2f} kN s/m")
    bx.plot(vv * 1e3, cd * vv * np.abs(vv), color="k", lw=0.9, ls=(0, (3, 2)),
            label=f"$c_D\\,\\dot u|\\dot u|$, $c_D$ = {cd:.1f} kN s$^2$/m$^2$")
    sg = np.sign(vc)
    fc = np.sum(sg * rm) / np.sum(sg * sg)
    vs = np.r_[np.linspace(vc.min(), -1e-4, 25), np.linspace(1e-4, vc.max(), 25)]
    bx.plot(vs * 1e3, fc * np.sign(vs), color="#CC79A7", lw=1.1, label=f"$F_c\,\mathrm{{sgn}}(\dot u)$, $F_c$ = {fc:.2f} kN")
    for lab_, pred in (("c v", c * vc), ("c_D v|v|", cd * q), ("F_c sgn(v)", fc * sg)):
        print(f"   {lab_:11}: rms misfit to binned means {np.sqrt(np.mean((rm - pred) ** 2)):.3f} kN")
    print(f"   F_c = {fc:.3f} kN (model scale {1e3 * fc / R.LAM ** 3:.1f} N)")
    print(f"   binned residual: c = {c:.2f} kN s/m, c_D = {cd:.1f} kN s2/m2; binned means range "
          f"{rm.min():+.2f} to {rm.max():+.2f} kN over v = {vc.min():+.2f} to {vc.max():+.2f} m/s")
    bx.axhline(0, color="0.85", lw=0.4, zorder=0)
    bx.axvline(0, color="0.85", lw=0.4, zorder=0)
    bx.set_ylim(-1.5, 1.5)
    bx.set_xlabel(r"$\dot u$ (mm/s)")
    bx.set_ylabel(r"$\lambda_L^3 f_p - m\ddot u - F_0$ (kN)")
    bx.legend(loc="upper left", frameon=False, fontsize=6.3, handlelength=1.6)
    bx.set_title(rf"\textbf{{({'bd'[row]})}}", loc="left")
fig.savefig(OUT.with_suffix(".png"), dpi=250)
fig.savefig(OUT.with_suffix(".pdf"), dpi=300)
print("wrote", OUT.with_suffix(".png"))
