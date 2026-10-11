"""Physical-subassembly force versus measured displacement, velocity, and acceleration (dissertation Ch. 6 draft;
read-only on data).

W02 (Wednesday, CPU MKR-alpha + ProfileSPD, 8 ranks, periodic waves then earthquake) and F05 (Friday, GPU MKR-alpha +
cuDSS, 8 ranks, no slowdowns, earthquake then solitary wave), over the 5-95 % Arias window of the ground motion.
Force: load cell, lambda_L^3 f_p = -lambda^3 * topLRaw (resisting-force convention, as fig_wed_force.py).
Displacement: actuator meaSigOS. Velocity and acceleration: central differences of the measured displacement at
2048 Hz after a zero-phase 4th-order Butterworth low-pass at F_CUT (model scale), applied to the displacement only.
All quantities at prototype scale (Froude): u_p = lambda u_m, v_p = sqrt(lambda) v_m, a_p = a_m, F_p = lambda^3 F_m.
Row 3 (printed and annotated): least-squares fit F = k u + c v + m a + F0 over the same window.
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
import fri_runs as R  # noqa: E402

sys.path.insert(0, str(R.ROOT / "plot"))
import PlotResponseSpectrum as prs  # noqa: E402
from lab_paths import MAT_EXTRACT_DIR  # noqa: E402

OUT = Path(sys.argv[1])
F_CUT = 40.0  # Hz, model scale
# test, mat stem, integration time of the start of the ground motion (prototype s)
CASES = [("W02", "0819_GusBridge_Parallel8_H0p5_T7p746_reg", 150.0, "W02 (Wed.), CPU, 8 ranks"),
         ("F05", "0821_GusBridge_rowNeg4", 0.0, "F05 (Fri.), GPU, 8 ranks")]


def load(stem, t_eq):
    z = np.load(MAT_EXTRACT_DIR / f"{stem}.npz", allow_pickle=True)
    st, me = z["stateOS_data"], z["meaSigOS_data"]
    t, s1 = st[:, -1], np.rint(st[:, 2])
    tk = t[1:][(s1[1:] > 0) & (s1[:-1] <= 0)]
    tk = tk[1 if tk[1] - tk[0] < 3e-3 else 0:]
    ti = np.interp(me[:, -1], tk, np.arange(tk.size) * R.DT_SIM * R.SQL)  # integration time of each sample
    m = sio.loadmat(R.DATA / "Simulink" / f"{stem}.mat", squeeze_me=True, struct_as_record=False)["data"]
    raw = m.topLRaw.data
    tm = me[:, -1]
    dt = np.median(np.diff(tm))
    b, a = butter(4, F_CUT * 2 * dt, "low")
    uf = filtfilt(b, a, me[:, 0])
    v = np.gradient(uf, dt)
    acc = np.gradient(v, dt)
    F = -R.LAM ** 3 * np.interp(tm, raw[:, -1], raw[:, 0]) / 1e3
    k = (ti >= t_eq + R.T_D5) & (ti <= t_eq + R.T_D95)
    return dict(u=me[k, 0] * R.LAM * 1e3, v=v[k] * R.SQL * 1e3, a=acc[k], F=F[k])


prs.configure_font()
fig, A = plt.subplots(2, 3, figsize=(6.0, 4.2), layout="constrained", sharey="row")
for row, (tid, stem, t_eq, name) in enumerate(CASES):
    d = load(stem, t_eq)
    X = np.c_[d["u"] / 1e3, d["v"] / 1e3, d["a"], np.ones(d["F"].size)]  # m, m/s, m/s^2
    coef, *_ = np.linalg.lstsq(X, d["F"], rcond=None)
    fit = X @ coef
    r2 = 1 - np.sum((d["F"] - fit) ** 2) / np.sum((d["F"] - d["F"].mean()) ** 2)
    rho = [np.corrcoef(d[q], d["F"])[0, 1] for q in ("u", "v", "a")]
    print(f"{tid}: k = {coef[0]:.1f} kN/m, c = {coef[1]:.1f} kN s/m, m = {coef[2]:.2f} kN s^2/m (= {coef[2]:.2f} t), "
          f"F0 = {coef[3]:.2f} kN; R^2 = {r2:.3f}; corr(F, u/v/a) = {np.round(rho, 3)}; "
          f"ranges u {np.ptp(d['u']):.0f} mm, v {np.ptp(d['v']):.0f} mm/s, a {np.ptp(d['a']):.1f} m/s^2, F {np.ptp(d['F']):.1f} kN")
    for j, (q, xl) in enumerate((("u", r"$u$ (mm)"), ("v", r"$\dot u$ (mm/s)"), ("a", r"$\ddot u$ (m/s$^2$)"))):
        ax = A[row, j]
        ax.plot(d[q], d["F"], color="k", lw=0.15, alpha=0.6, rasterized=True)
        ax.axhline(0, color="0.8", lw=0.4, zorder=0)
        ax.axvline(0, color="0.8", lw=0.4, zorder=0)
        ax.text(0.04, 0.95, rf"$\rho = {rho[j]:+.2f}$", transform=ax.transAxes, ha="left", va="top", fontsize=7)
        ax.set_xlabel(xl)
        ax.set_title(rf"\textbf{{({'abcdef'[3 * row + j]})}}" + (f" {name}" if j == 0 else ""), loc="left")
    A[row, 0].set_ylabel(r"$\lambda_L^3 f_p$ (kN)")
    A[row, 2].text(0.96, 0.05, f"$F \\approx m\\ddot u$: $m$ = {coef[2]:.1f} t\n$R^2$ = {r2:.2f} (with $k$, $c$)",
                   transform=A[row, 2].transAxes, ha="right", va="bottom", fontsize=6.5, color="0.3")
fig.savefig(OUT.with_suffix(".png"), dpi=250)
fig.savefig(OUT.with_suffix(".pdf"), dpi=300)
print("wrote", OUT.with_suffix(".png"))
