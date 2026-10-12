"""Interface force against measured displacement, velocity, and acceleration, Wednesday test W05 (16 ranks, no
slowdowns after start-up), dissertation Ch. 6; read-only on data.

Window: 5-95 % Arias intensity of the ground motion (t_int = 206.0-289.1 s), with the periodic waves still acting.
Force: load cell, lambda_L^3 f_p = -lambda^3 * topLRaw (resisting-force convention, as fig_wed_force.py).
Displacement: actuator meaSigOS; velocity and acceleration by central differences at 2048 Hz after a zero-phase
4th-order Butterworth low-pass at 20 Hz (model scale) on the displacement only (without it, differentiated sensor
noise dominates the acceleration; the fitted mass is stable for cutoffs of 5-20 Hz, and the unfiltered command
displacement gives nearly the same mass). Prototype scale (Froude).
Line in (c): F = m a + F0, least squares on the acceleration alone.
One row; bold panel tags; uniform font; constrained layout.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy.signal import butter, filtfilt

sys.path.insert(0, str(Path(__file__).parent))
import ctrl_error_wed as C  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "plot"))
import PlotResponseSpectrum as prs  # noqa: E402
from lab_paths import MAT_EXTRACT_DIR, CYLINDER_LENGTH_SCALE as LAM, TIME_SCALE_FROUDE as SQL  # noqa: E402

OUT = Path(sys.argv[1])
STEM = "0819_GusBridge_Parallel16_H0p5_T7p746_reg_Trial02"
T_LO, T_HI = 150.0 + 56.0, 150.0 + 139.05  # t_int, 5-95 % Arias
F_CUT = 20.0
ORANGE = "#C0390B"

z = np.load(MAT_EXTRACT_DIR / f"{STEM}.npz", allow_pickle=True)
st, me = z["stateOS_data"], z["meaSigOS_data"]
t, s1 = st[:, -1], np.rint(st[:, 2])
tk = t[1:][(s1[1:] > 0) & (s1[:-1] <= 0)]
tk = tk[1 if tk[1] - tk[0] < 3e-3 else 0:]
tm = me[:, -1]
ti = np.interp(tm, tk, np.arange(tk.size) * C.DT_SIM * SQL)
dt = np.median(np.diff(tm))
b, a = butter(4, F_CUT * 2 * dt, "low")
uf = filtfilt(b, a, me[:, 0])
vel = np.gradient(uf, dt)
acc = np.gradient(vel, dt)
raw = sio.loadmat(ROOT / "OSU_SSI_BRIDGE_DATA" / "Simulink" / f"{STEM}.mat", squeeze_me=True,
                  struct_as_record=False)["data"].topLRaw.data
F = -LAM ** 3 * np.interp(tm, raw[:, -1], raw[:, 0]) / 1e3
k = (ti >= T_LO) & (ti <= T_HI)
u, v, A, F = me[k, 0] * LAM * 1e3, vel[k] * SQL * 1e3, acc[k], F[k]
m, f0 = np.linalg.lstsq(np.c_[A, np.ones(A.size)], F, rcond=None)[0]
r2 = 1 - np.sum((F - m * A - f0) ** 2) / np.sum((F - F.mean()) ** 2)
print(f"W05: m = {m:.2f} t (prototype; {1e3 * m / LAM ** 3:.0f} kg model), F0 = {f0:.2f} kN, R^2 = {r2:.2f}; "
      f"corr(F, u/v/a) = {[round(np.corrcoef(x, F)[0, 1], 2) for x in (u, v, A)]}")

prs.configure_font()
fig, axs = plt.subplots(1, 3, figsize=(6.0, 2.3), layout="constrained", sharey=True)
for ax, x, xl, tg in zip(axs, (u, v, A), (r"$u$ (mm)", r"$\dot u$ (mm/s)", r"$\ddot u$ (m/s$^2$)"),
                         ("(a)", "(b)", "(c)")):
    ax.plot(x, F, color="0.45", lw=0.15, alpha=0.6, rasterized=True)
    ax.axhline(0, color="0.85", lw=0.4, zorder=0)
    ax.axvline(0, color="0.85", lw=0.4, zorder=0)
    ax.set_xlabel(xl)
    ax.set_title(rf"\textbf{{{tg}}}", loc="left")
xx = np.array([A.min(), A.max()])
axs[2].plot(xx, m * xx + f0, color=ORANGE, lw=1.4)
axs[2].annotate(rf"$m\ddot u$, $m$ = {m:.1f} t", (0.85 * xx[1], m * 0.85 * xx[1] + f0), xytext=(0.15, 7.5),
                color=ORANGE, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.6, shrinkA=1, shrinkB=1))
axs[0].set_ylabel(r"$\lambda_L^3 f_p$ (kN)")
fig.savefig(OUT.with_suffix(".png"), dpi=250)
fig.savefig(OUT.with_suffix(".pdf"), dpi=300)
print("wrote", OUT.with_suffix(".png"))
