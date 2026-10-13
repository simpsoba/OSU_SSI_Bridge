"""Solitary wave on the physical subassembly, F05 (dissertation Ch. 6 draft; read-only on data).

(a) interface force against the measured displacement from 5 s before to 30 s after the peak load-cell force; circle =
    peak force; arrows at half the peak force, beside the curve, show the sense of time
(b) displacement from 40 s before to 70 s after the peak force (circle = the peak force), with the amplitude at f_1
    (band 0.3-0.7 Hz) over the 10 s on each side of it
No derivatives of the displacement: during the wave impact the measured velocity and acceleration depend on the
low-pass cutoff (the transfer-system response), so they are not shown.
Force: load cell, lambda_L^3 f_p = -lambda^3 * topLRaw (cache_force.py); displacement: actuator meaSigOS; prototype scale
on the prototype-scaled laboratory clock. Usage: fig_fri_sol_kin.py <out.png> [test]
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, sosfiltfilt

sys.path.insert(0, str(Path(__file__).parent))
import fri_runs as R  # noqa: E402

sys.path.insert(0, str(R.ROOT / "plot"))
import PlotResponseSpectrum as prs  # noqa: E402

OUT = Path(sys.argv[1])
TID = sys.argv[2] if len(sys.argv) > 2 else "F05"
BLUE = "#0072B2"
stem = R.BY_ID[TID][1]
me = np.load(R.MAT_EXTRACT_DIR / f"{stem}.npz", allow_pickle=True)["meaSigOS_data"]
tm = me[:, -1]
dt = np.median(np.diff(tm))
tf, fr = R.load_force(stem)
F = -R.LAM ** 3 * np.interp(tm, tf, fr) / 1e3
tl = tm * R.SQL
u = me[:, 0] * R.LAM * 1e3
w = (tl > 300) & (tl < 380)
t_imp = tl[w][np.argmin(F[w])]
k = (tl >= t_imp - 5.0) & (tl <= t_imp + 30.0)
Fk, uk = F[k], u[k]
ip = np.argmin(Fk)
print(f"{TID}: peak {Fk[ip]:.1f} kN at {t_imp:.2f} s; u range {uk.min():.1f}..{uk.max():.1f} mm")

prs.configure_font()
fig, (ax, D) = plt.subplots(1, 2, figsize=(6.0, 2.5), layout="constrained", gridspec_kw=dict(width_ratios=[1, 2.3]))
ax.plot(uk, Fk, color=BLUE, lw=0.5)
ax.plot(uk[ip], Fk[ip], "o", ms=3.5, color=BLUE)
ax.axhline(0, color="0.85", lw=0.4, zorder=0)
ax.set_xlabel(r"$u$ (mm)")
ax.set_ylabel(r"$\lambda_L^3 f_p$ (kN)")
ax.set_title(r"\textbf{(a)}", loc="left")
# arrows beside the curve at half the peak force, before and after the peak, pointing the way time advances
half = 0.5 * Fk[ip]
ib = np.flatnonzero(Fk[:ip] > half)[-1]
ia = ip + np.flatnonzero(Fk[ip:] > half)[0]
fig.canvas.draw()
cen = ax.transData.transform(np.c_[uk[ib:ia + 1], Fk[ib:ia + 1]]).mean(axis=0)
bb = ax.get_window_extent()
for i in (ib, ia):
    p0, p1 = ax.transData.transform((uk[i - 3], Fk[i - 3])), ax.transData.transform((uk[i + 3], Fk[i + 3]))
    d = (p1 - p0) / np.hypot(*(p1 - p0))
    c = ax.transData.transform((uk[i], Fk[i]))
    n = np.array([-d[1], d[0]])
    n = n if np.dot(n, c - cen) > 0 else -n
    if not (bb.x0 + 12 < (c + 9 * n)[0] < bb.x1 - 12 and bb.y0 + 12 < (c + 9 * n)[1] < bb.y1 - 12):
        n = -n
    c = c + 9 * n
    tail, head = ax.transData.inverted().transform(c - 8 * d), ax.transData.inverted().transform(c + 8 * d)
    ax.annotate("", xy=head, xytext=tail, arrowprops=dict(arrowstyle="-|>", color="k", lw=1.0, mutation_scale=9,
                                                          shrinkA=0, shrinkB=0))
# (b) free vibration before and after the wave
kd = (tl >= t_imp - 40) & (tl <= t_imp + 70)
D.plot(tl[kd] - t_imp, u[kd], color=BLUE, lw=0.5, zorder=3)
D.plot(0.0, uk[ip], "o", ms=3.5, color=BLUE, zorder=4)
ub = sosfiltfilt(butter(4, [0.3, 0.7], btype="band", fs=1 / (dt * R.SQL), output="sos"), u)
amp = {}
for nm, (w0, w1) in (("before", (-12, -2)), ("after", (2, 12))):
    kw = (tl - t_imp >= w0) & (tl - t_imp <= w1)
    am, mu = np.abs(ub[kw]).max(), u[kw].mean()
    amp[nm] = am
    for yv in (mu - am, mu + am):
        D.plot([w0, w1], [yv, yv], color="k", lw=0.9, ls=(0, (4, 1.5)))
    D.text(0.5 * (w0 + w1), mu + am + 0.3, rf"$\pm${am:.1f} mm", ha="center", va="bottom", color="k")
print(f"amplitude at f1: before {amp['before']:.2f} mm, after {amp['after']:.2f} mm, ratio {amp['after'] / amp['before']:.1f}")
D.set_xlim(-40, 70)
D.set_ylim(D.get_ylim()[0], D.get_ylim()[1] + 1.0)
D.set_xlabel(r"Time from the peak force, $\sqrt{\lambda_L}\,\tau$ (s)")
D.set_ylabel(r"$u$ (mm)")
D.set_title(r"\textbf{(b)}", loc="left")
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
