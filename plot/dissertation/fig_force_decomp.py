"""Decomposition of the physical-subassembly force into stiffness, damping, and inertial terms (dissertation Ch. 6
draft; read-only on data).

Least-squares fit F = k u + c v + m a + F0 over the 5-95 % Arias window, with u, v, a from the measured displacement
(20 Hz zero-phase low-pass at model scale, central differences), as fig_force_kinematics.py; W05 and F05; prototype
scale, resisting-force convention. In each panel: grey = measured force against u, v, or a; color = the fitted term
of that column (k u, c v, or m a), plus F0; in the u panels, faint = m a against u. Text: coefficient, standard deviation of the term, and its share of the
force variance, cov(term, F) / var(F); the three shares add up to the R^2 of the fit.
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
COLS = [("u", r"$u$ (mm)", 1e-3, "#0072B2", "k", "kN/m"), ("v", r"$\dot u$ (mm/s)", 1e-3, "#009E73", "c", "kN s/m"),
        ("a", r"$\ddot u$ (m/s$^2$)", 1.0, "#D55E00", "m", "t")]

ns["prs"].configure_font()
fig, A = plt.subplots(2, 3, figsize=(6.0, 4.4), layout="constrained", sharey="row")
for row, (tid, stem, t_eq, name) in enumerate(ns["CASES"]):
    d = ns["load"](stem, t_eq)
    X = np.c_[d["u"] / 1e3, d["v"] / 1e3, d["a"], np.ones(d["F"].size)]
    coef = np.linalg.lstsq(X, d["F"], rcond=None)[0]
    fit = X @ coef
    r2 = 1 - np.sum((d["F"] - fit) ** 2) / np.sum((d["F"] - d["F"].mean()) ** 2)
    rua = np.corrcoef(d["u"], d["a"])[0, 1]
    print(f"{tid}: k = {coef[0]:.2f} kN/m, c = {coef[1]:.2f} kN s/m, m = {coef[2]:.2f} t, F0 = {coef[3]:.2f} kN; "
          f"R2 = {r2:.2f}; corr(u, a) = {rua:+.2f}")
    for j, (q, xl, sc, col, sym, unit) in enumerate(COLS):
        ax = A[row, j]
        term = coef[j] * d[q] * sc
        ax.plot(d[q], d["F"], color="0.6", lw=0.15, alpha=0.6, rasterized=True)
        if q == "u":  # the inertial term seen against u: the apparent u dependence (a ~ -w^2 u near the first mode)
            ax.plot(d["u"], coef[2] * d["a"] + coef[3], color="#D55E00", lw=0.15, alpha=0.35, rasterized=True)
            ax.text(0.96, 0.04, r"faint: $m\ddot u$ vs $u$", transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=6.3, color="#D55E00")
        o = np.argsort(d[q])
        ax.plot(d[q][o], term[o] + coef[3], color=col, lw=1.6)
        share = np.cov(term, d["F"])[0, 1] / np.var(d["F"], ddof=1)
        ax.text(0.04, 0.96, f"${sym}$ = {coef[j]:.2f} {unit}\nstd({sym}{q}) = {term.std():.2f} kN\nshare: "
                f"{100 * share:.0f}\\% of var.", transform=ax.transAxes, ha="left", va="top", fontsize=6.3,
                color=col, linespacing=1.1)
        ax.axhline(0, color="0.85", lw=0.4, zorder=0)
        ax.axvline(0, color="0.85", lw=0.4, zorder=0)
        ax.set_xlabel(xl)
        ax.set_title(rf"\textbf{{({'abcdef'[3 * row + j]})}}" + (f" {name}" if j == 0 else ""), loc="left")
        print(f"   {sym}{q}: std {term.std():.2f} kN, share {100 * share:.1f}% of variance")
    A[row, 0].set_ylabel(r"$\lambda_L^3 f_p$ (kN)")
    A[row, 2].text(0.96, 0.04, f"fit $R^2$ = {r2:.2f}\n$\\rho(u, \\ddot u)$ = {rua:+.2f}", transform=A[row, 2].transAxes,
                   ha="right", va="bottom", fontsize=6.3, color="0.3")
fig.savefig(OUT.with_suffix(".png"), dpi=250)
fig.savefig(OUT.with_suffix(".pdf"), dpi=300)
print("wrote", OUT.with_suffix(".png"))
