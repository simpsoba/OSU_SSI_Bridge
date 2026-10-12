"""Interface force against measured displacement, velocity, and acceleration over the whole record, colored by phase
(dissertation Ch. 6; read-only on data). One figure per test: this script makes the W05 figure; fig_fri_force_kin.py
calls make("F05", ...) for the Friday one.

W05 (Wednesday, CPU MKR-alpha + ProfileSPD, 16 ranks, no slowdowns after start-up): before the earthquake
(periodic waves only, t_int < 150 s), during it (150 s to 95 % Arias intensity, 289.1 s; waves still acting), after it.
F05 (Friday, GPU MKR-alpha + cuDSS, 8 ranks, no slowdowns): during the earthquake (0 to 139.05 s, still water),
after it (free vibration), and the solitary wave (from 5 s before to 30 s after the peak load-cell force, prototype
laboratory clock). The solitary-wave peak (about -77 kN) is off scale and marked.
Records start after start-up (t_int > 2 s) and end at the last new target.
Force: load cell, lambda_L^3 f_p = -lambda^3 * topLRaw (resisting-force convention, as fig_wed_force.py).
Displacement: actuator meaSigOS; velocity and acceleration by central differences at 2048 Hz after a zero-phase
4th-order Butterworth low-pass at 20 Hz (model scale) on the displacement only (without it, differentiated sensor
noise dominates the acceleration; the fitted mass is stable for cutoffs of 5-20 Hz, and the unfiltered command
displacement gives nearly the same mass). Prototype scale (Froude).
Line in (c): F = m a + F0, least squares on the acceleration alone during the earthquake (per-phase values
printed).
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

F_CUT = 20.0
RED = "#C0390B"
C_DUR, C_AFT, C_WAV, C_SOL = "0.35", "#E69F00", "#56B4E9", "#0072B2"
YLIM = (-11.0, 11.0)
CASES = [("W05", "0819_GusBridge_Parallel16_H0p5_T7p746_reg_Trial02", 150.0, "W05 (Wed.)"),
         ("F05", "0821_GusBridge_rowNeg4", 0.0, "F05 (Fri.)")]
T95 = 139.05  # s after the start of the ground motion


def load(stem):
    z = np.load(MAT_EXTRACT_DIR / f"{stem}.npz", allow_pickle=True)
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
    raw = sio.loadmat(ROOT / "OSU_SSI_BRIDGE_DATA" / "Simulink" / f"{stem}.mat", squeeze_me=True,
                      struct_as_record=False)["data"].topLRaw.data
    F = -LAM ** 3 * np.interp(tm, raw[:, -1], raw[:, 0]) / 1e3
    k = (ti > 2.0) & (tm <= tk[-1])
    return dict(tl=tm[k] * SQL, ti=ti[k], u=me[k, 0] * LAM * 1e3, v=vel[k] * SQL * 1e3, a=acc[k], F=F[k])


def phases(tid, d, t_eq):
    ti, tl = d["ti"], d["tl"]
    dur = (ti >= t_eq) & (ti < t_eq + T95)
    if tid.startswith("W"):
        return [("during the earthquake", dur, C_DUR), ("after", ti >= t_eq + T95, C_AFT),
                ("before (waves only)", ti < t_eq, C_WAV)]
    w = (tl > 300) & (tl < 380)
    t_imp = tl[w][np.argmin(d["F"][w])]
    sol = (tl >= t_imp - 5.0) & (tl <= t_imp + 30.0)
    return [("during the earthquake", dur, C_DUR), ("after", (ti >= t_eq + T95) & ~sol, C_AFT),
            ("solitary wave", sol, C_SOL)]


def make(test, out):
    tid, stem, t_eq, name = next(c for c in CASES if c[0] == test)
    prs.configure_font()
    fig, A = plt.subplots(1, 3, figsize=(6.0, 2.5), layout="constrained", sharey=True)
    A = A[None, :]
    handles = {}
    row = 0
    d = load(stem)
    P = phases(tid, d, t_eq)
    dur = P[0][1]
    m, f0 = np.linalg.lstsq(np.c_[d["a"][dur], np.ones(dur.sum())], d["F"][dur], rcond=None)[0]
    print(f"{tid}: m (earthquake) = {m:.2f} t ({1e3 * m / LAM ** 3:.0f} kg model)")
    for pname, s_, _ in P:
        mp = np.linalg.lstsq(np.c_[d["a"][s_], np.ones(s_.sum())], d["F"][s_], rcond=None)[0][0]
        print(f"  {pname:22}: t_int {d['ti'][s_].min():6.1f}-{d['ti'][s_].max():6.1f} s, m = {mp:5.2f} t, "
              f"corr(F, a) = {np.corrcoef(d['a'][s_], d['F'][s_])[0, 1]:+.2f}, |F| max {np.abs(d['F'][s_]).max():5.1f} kN")
    for j, (q, xl) in enumerate((("u", r"$u$ (mm)"), ("v", r"$\dot u$ (mm/s)"), ("a", r"$\ddot u$ (m/s$^2$)"))):
        ax = A[row, j]
        for pname, s_, col in P:
            thin = col in (C_WAV, C_SOL)
            ax.plot(d[q][s_], d["F"][s_], color=col, lw=0.4 if thin else 0.15, alpha=0.9 if thin else 0.6,
                    rasterized=True)
            handles.setdefault(pname, plt.Line2D([], [], color=col, lw=1.5))
        ax.axhline(0, color="0.85", lw=0.4, zorder=0)
        ax.axvline(0, color="0.85", lw=0.4, zorder=0)
        ax.set_xlabel(xl)
        ax.set_ylim(*YLIM)
        ax.set_title(rf"\textbf{{({'abc'[j]})}}" + (f" {name}" if j == 0 else ""), loc="left")
        if tid.startswith("F"):
            s_ = P[2][1]
            ax.text(0.97, 0.04, f"solitary-wave peak\n{d['F'][s_].min():.0f} kN, off scale", transform=ax.transAxes,
                    ha="right", va="bottom", fontsize=6, color=C_SOL)
    ax = A[row, 2]
    xx = np.array([d["a"][dur].min(), d["a"][dur].max()])
    ax.plot(xx, m * xx + f0, color=RED, lw=1.4)
    ax.annotate(rf"$m\ddot u$, $m$ = {m:.1f} t", (0.85 * xx[1], m * 0.85 * xx[1] + f0), xytext=(0.05, 0.93),
                textcoords="axes fraction", color=RED, ha="left", va="top",
                arrowprops=dict(arrowstyle="-", color=RED, lw=0.6, shrinkA=1, shrinkB=1))
    A[row, 0].set_ylabel(r"$\lambda_L^3 f_p$ (kN)")
    order = [n for n in ("before (waves only)", "during the earthquake", "after", "solitary wave") if n in handles]
    fig.legend([handles[n] for n in order], order, loc="outside upper center", ncol=4, frameon=False, handlelength=1.5,
               columnspacing=1.2)
    fig.savefig(out.with_suffix(".png"), dpi=250)
    fig.savefig(out.with_suffix(".pdf"), dpi=300)
    print("wrote", out.with_suffix(".png"))


if __name__ == "__main__":
    make("W05", Path(sys.argv[1]))
