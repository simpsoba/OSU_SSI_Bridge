"""Pier-base zero-length section of Test F30: moment against rotation, before and after the solitary wave
(dissertation Ch. 6 draft; read-only on data). Grey: from the start of the record to T_PRE before the peak load-cell
force of the solitary wave; red: from then to the end of the record. (b) zooms on the red part.
pier_hinge_defo/force.out: columns (time, axial, rotation) and (time, axial force, moment), as fig_wed_energy2.py.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fri_runs as R  # noqa: E402

sys.path.insert(0, str(R.ROOT / "plot"))
import PlotResponseSpectrum as prs  # noqa: E402

OUT = Path(sys.argv[1])
TEST, DUMP = "F30", R.DATA / "r-17_20260821_1649_Storm_Wave"
T_PRE = 2.0
RED = "#C0390B"


def rd(p):
    a = np.genfromtxt(p, invalid_raise=False)
    return a[np.isfinite(a).all(1)]


d = R.load(TEST)
m = sio.loadmat(R.DATA / "Simulink" / f"{R.BY_ID[TEST][1]}.mat", squeeze_me=True, struct_as_record=False)["data"]
raw = m.topLRaw.data
tl, F = raw[:, -1], -R.LAM ** 3 * raw[:, 0] / 1e3
w = (tl > 200) & (tl < 240)
t_imp = float(d["to_int"](tl[w][np.argmax(np.abs(F[w]))]))
hd, hf = rd(DUMP / "pier_hinge_defo.out.0"), rd(DUMP / "pier_hinge_force.out.0")
n = min(len(hd), len(hf))
t, th, M = hd[:n, 0], 1e3 * hd[:n, 2], 1e-6 * hf[:n, 2]
pre, post = t < t_imp - T_PRE, t >= t_imp - T_PRE
print(f"{TEST}: impact t_int {t_imp:.2f} s; rotation range before {th[pre].min():.2f} to {th[pre].max():.2f} mrad, "
      f"moment {M[pre].min():.2f} to {M[pre].max():.2f} MN m; after: rotation {th[post].min():.3f} to "
      f"{th[post].max():.3f} mrad, moment {M[post].min():.3f} to {M[post].max():.3f} MN m")

prs.configure_font()
fig, (ax, az) = plt.subplots(1, 2, figsize=(6.0, 2.7), layout="constrained")
ax.plot(th[pre], M[pre], color="0.7", lw=0.5, label="before the wave")
ax.plot(th[post], M[post], color=RED, lw=0.9, label="from the wave on")
x0, x1 = th[post].min(), th[post].max()
y0, y1 = M[post].min(), M[post].max()
px, py = 0.15 * (x1 - x0) + 0.05, 0.15 * (y1 - y0) + 0.02
ax.add_patch(plt.Rectangle((x0 - px, y0 - py), x1 - x0 + 2 * px, y1 - y0 + 2 * py, fc="none", ec="0.3", lw=0.6, ls="--"))
ax.axhline(0, color="0.85", lw=0.4, zorder=0)
ax.axvline(0, color="0.85", lw=0.4, zorder=0)
ax.set_xlabel(r"$\theta$ (mrad)")
ax.set_ylabel(r"$M$ (MN$\cdot$m)")
ax.legend(loc="upper left", frameon=False, fontsize=7, handlelength=1.4)
ax.set_title(r"\textbf{(a)}", loc="left")
az.plot(th[pre], M[pre], color="0.7", lw=0.5)
az.plot(th[post], M[post], color=RED, lw=0.9)
az.set_xlim(x0 - px, x1 + px)
az.set_ylim(y0 - py, y1 + py)
az.set_xlabel(r"$\theta$ (mrad)")
az.set_title(r"\textbf{(b)}", loc="left")
fig.savefig(OUT.with_suffix(".png"), dpi=250)
fig.savefig(OUT.with_suffix(".pdf"), dpi=300)
print("wrote", OUT.with_suffix(".png"))
