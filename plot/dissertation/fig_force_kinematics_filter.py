"""Force versus measured displacement, velocity, and acceleration, without and with the low-pass filter on the
displacement (dissertation Ch. 6 draft; read-only on data). Same data, window, and scaling as fig_force_kinematics.py
(W05 and F05, 5-95 % Arias window, prototype scale); rows alternate unfiltered / 20 Hz (model scale, zero phase).
Velocity and acceleration are central differences at 2048 Hz in both cases; the force is never filtered.
Option --negate: plot -lambda_L^3 f_p (the force the physical subassembly exerts on the numerical one, i.e., the
external force on the numerical subassembly) instead of the resisting force lambda_L^3 f_p.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
SCRIPT = Path(__file__).parent / "fig_force_kinematics.py"
OUT = Path(sys.argv[1])
SIGN = -1.0 if "--negate" in sys.argv[2:] else 1.0
YLAB = r"$-\lambda_L^3 f_p$ (kN)" if SIGN < 0 else r"$\lambda_L^3 f_p$ (kN)"
src = SCRIPT.read_text(encoding="utf-8")
ns = {"__file__": str(SCRIPT)}
exec(compile(src[:src.index("prs.configure_font()")], str(SCRIPT), "exec"), ns)

ns["prs"].configure_font()
fig, A = plt.subplots(4, 3, figsize=(6.0, 7.6), layout="constrained")
row = 0
for tid, stem, t_eq, name in ns["CASES"]:
    for fc, flab in ((1000.0, "unfiltered"), (20.0, "20 Hz low-pass")):
        ns["F_CUT"] = fc  # 1000 Hz is just below Nyquist (1024 Hz): no filtering in effect
        d = ns["load"](stem, t_eq)
        d["F"] = SIGN * d["F"]
        X = np.c_[d["u"] / 1e3, d["v"] / 1e3, d["a"], np.ones(d["F"].size)]
        coef = np.linalg.lstsq(X, d["F"], rcond=None)[0]
        for j, (q, xl) in enumerate((("u", r"$u$ (mm)"), ("v", r"$\dot u$ (mm/s)"), ("a", r"$\ddot u$ (m/s$^2$)"))):
            ax = A[row, j]
            ax.plot(d[q], d["F"], color="k", lw=0.15, alpha=0.6, rasterized=True)
            ax.axhline(0, color="0.8", lw=0.4, zorder=0)
            ax.axvline(0, color="0.8", lw=0.4, zorder=0)
            rho = np.corrcoef(d[q], d["F"])[0, 1]
            ax.text(0.04, 0.95, rf"$\rho = {rho:+.2f}$", transform=ax.transAxes, ha="left", va="top", fontsize=7)
            ax.set_xlabel(xl)
            if j > 0:
                ax.tick_params(labelleft=False)
            ax.set_title(rf"\textbf{{({'abcdefghijkl'[3 * row + j]})}}" + (f" {tid}, {flab}" if j == 0 else ""),
                         loc="left")
        A[row, 0].set_ylabel(YLAB)
        mtxt = f"fitted $m$ = {coef[2]:.2f} t" if SIGN > 0 else f"slope $-m$, $m$ = {-coef[2]:.2f} t"
        A[row, 2].text(0.96, 0.05, mtxt, transform=A[row, 2].transAxes, ha="right",
                       va="bottom", fontsize=6.5, color="0.3")
        print(f"{tid} {flab:15}: v range {np.ptp(d['v']):7.0f} mm/s, a range {np.ptp(d['a']):7.1f} m/s2, "
              f"rho(F,v) {np.corrcoef(d['v'], d['F'])[0, 1]:+.2f}, rho(F,a) {np.corrcoef(d['a'], d['F'])[0, 1]:+.2f}, "
              f"m {coef[2]:.2f} t")
        row += 1
for r in range(4):
    ylim = A[r, 0].get_ylim()
    for j in (1, 2):
        A[r, j].set_ylim(ylim)
fig.savefig(OUT.with_suffix(".png"), dpi=250)
fig.savefig(OUT.with_suffix(".pdf"), dpi=300)
print("wrote", OUT.with_suffix(".png"))
