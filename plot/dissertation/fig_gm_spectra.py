"""Spectra of the input ground motion, Tohoku 2011 FKSH19 NS (dissertation Ch. 6 draft; read-only on data).

Velocity record (PEER VT2, dt = 0.01 s, 300 s; prototype scale, scale factor 1) as applied at the base of the soil
column; acceleration by central differences, displacement by integration. Two signals: the whole record (0-300 s)
and its coda alone (150-300 s, oscillators starting from rest). Shaded: 0.1-0.3 Hz (T = 3.3-10 s), the long-period
band of fig_fri_psd_phases.py. Vertical lines: T_1 = 1/f_1 = 2.10 s and the period during strong shaking,
1/(0.70 f_1) = 2.98 s.
(a)-(c) 5 %-damped elastic response spectra: pseudo-acceleration S_a, pseudo-velocity S_v, displacement S_d
(d)     PSD of the ground displacement (Welch, Hann 40 s, 75 % overlap), against frequency
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import welch

sys.path.insert(0, str(Path(__file__).parent))
import fri_runs as R  # noqa: E402

sys.path.insert(0, str(R.ROOT / "plot"))
import PlotResponseSpectrum as prs  # noqa: E402

OUT = Path(sys.argv[1])
GM = R.ROOT / "ground-motion" / "Tohoku2011-FKSH" / "FKSH19.NS1.VT2"
F1 = 0.477
T1, T1S = 1 / F1, 1 / (0.336)
BAND_F = (0.1, 0.3)
SIG = [("whole record", 0.0, 300.0, "0.0", "-"), ("coda, 150-300 s", 150.0, 300.0, "#08519C", "-")]

prs.configure_font()
v, dt = prs.read_peer_vt2(GM)
a = prs.accel_from_vel(v, dt)
d = np.cumsum(v) * dt
t = np.arange(v.size) * dt
T = np.logspace(np.log10(0.05), np.log10(20.0), 90)
fig, A = plt.subplots(2, 2, figsize=(6.0, 4.8), layout="constrained")
for name, t0, t1, col, ls in SIG:
    k = (t >= t0) & (t < t1)
    Sd, Sv, Sa, how = prs.compute_spectrum(a[k], dt, T, 0.05)
    print(f"{name}: {how}; Sd at T1 {Sd[np.argmin(abs(T - T1))] * 1e3:.1f} mm, in band max Sd "
          f"{Sd[(T > 1 / BAND_F[1]) & (T < 1 / BAND_F[0])].max() * 1e3:.1f} mm")
    for ax, y, sc in ((A[0, 0], Sa, 1.0), (A[0, 1], Sv, 100.0), (A[1, 0], Sd, 1000.0)):
        ax.loglog(T, y * sc, color=col, ls=ls, lw=0.9)
    x = d[k] - d[k].mean()
    n = int(40 / dt)
    f, P = welch(x * 1e3, fs=1 / dt, window="hann", nperseg=n, noverlap=3 * n // 4)
    A[1, 1].loglog(f[1:], P[1:], color=col, ls=ls, lw=0.9)
    kb = (f >= BAND_F[0]) & (f <= BAND_F[1])
    print(f"   ground-displacement PSD: strongest in-band peak at {f[kb][np.argmax(P[kb])]:.3f} Hz; "
          f"overall peak (f > 0.05 Hz) at {f[f > 0.05][np.argmax(P[f > 0.05])]:.3f} Hz")
for ax, yl in ((A[0, 0], r"$S_a$ (g)"), (A[0, 1], r"$S_v$ (cm/s)"), (A[1, 0], r"$S_d$ (mm)")):
    ax.axvspan(1 / BAND_F[1], 1 / BAND_F[0], color="#DCE6F0", lw=0, zorder=-1)
    for x0, ls_ in ((T1, (0, (3, 2))), (T1S, (0, (1, 1.5)))):
        ax.axvline(x0, color="0.35", lw=0.6, ls=ls_, zorder=0)
    ax.set_xlim(T[0], T[-1])
    ax.set_xlabel("Period, $T$ (s)")
    ax.set_ylabel(yl)
A[1, 1].axvspan(*BAND_F, color="#DCE6F0", lw=0, zorder=-1)
for x0, ls_ in ((F1, (0, (3, 2))), (0.336, (0, (1, 1.5)))):
    A[1, 1].axvline(x0, color="0.35", lw=0.6, ls=ls_, zorder=0)
A[1, 1].set_xlim(0.02, 10)
A[1, 1].set_ylim(1e-6, 1e4)
A[1, 1].set_xlabel("Frequency, $f$ (Hz)")
A[1, 1].set_ylabel(r"PSD of $u_g$ (mm$^2$/Hz)")
for ax, s in zip(A.flat, ("(a)", "(b)", "(c)", "(d)")):
    ax.set_title(rf"\textbf{{{s}}}", loc="left")
A[0, 0].text(T1, 1.02, r"$T_1$", transform=A[0, 0].get_xaxis_transform(), ha="center", va="bottom", color="0.3")
A[0, 0].text(T1S * 1.08, 1.02, r"$T_1/0.70$", transform=A[0, 0].get_xaxis_transform(), ha="left", va="bottom", color="0.3")
A[1, 0].text(np.sqrt(10 / 0.3), 0.04, "long-period\nband", transform=A[1, 0].get_xaxis_transform(), ha="center",
             va="bottom", color="#3B5B7A", linespacing=0.95)
A[0, 1].text(0.04, 0.95, "whole record", transform=A[0, 1].transAxes, color="0.0")
A[0, 1].text(0.04, 0.85, "coda, 150-300 s", transform=A[0, 1].transAxes, color="#08519C")
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
