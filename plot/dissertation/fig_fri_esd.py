"""PSD and ESD of the interface force and displacement by phase, F03 (dissertation Ch. 6 draft; read-only on data).

F03: CPU MKR-alpha + ProfileSPD, baseline model, 8 ranks, one slowdown (so laboratory time ~ integration time).
Signals as fig_wed_psd.py (Simulink at 2048 Hz model, prototype-scaled clock; lambda_L^3 f_p in kN, u in mm).
Phases: earthquake (0-140 s), free vibration (140 s to 25 s before the solitary-wave peak), solitary wave and after
(to the end of the record). One spectrogram of the whole record (Hann segments of 40 s, 75 % overlap, so the sum of
w^2 over the segments is constant; each segment's mean removed); each segment is assigned to the phase that contains
its centre. ESD = (hop/fs) x sum of the segment PSDs, so the phase ESDs add exactly to the whole-record ESD; the PSD of
a phase is its ESD divided by its duration. The phase boundaries are blurred by half a segment (20 s).
(a), (b) PSD (mean square per Hz): whole record and each phase
(c), (d) ESD (energy per Hz, PSD x duration): whole record, each phase, and the sum of the phases
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy.signal import spectrogram

sys.path.insert(0, str(Path(__file__).parent))
import fri_runs as R  # noqa: E402

sys.path.insert(0, str(R.ROOT / "plot"))
import PlotResponseSpectrum as prs  # noqa: E402

OUT = Path(sys.argv[1])
TID = "F03"
TSEG = 40.0
PHASE_COL = ["0.55", "#08519C", "#6A3D9A"]


def load(tid):
    m = sio.loadmat(R.DATA / "Simulink" / f"{R.BY_ID[tid][1]}.mat", squeeze_me=True, struct_as_record=False)["data"]
    me, raw = m.meaSigOS.data, m.topLRaw.data
    return dict(tu=me[:, -1] * R.SQL, u=me[:, 0] * R.LAM * 1e3, tf=raw[:, -1] * R.SQL, F=-R.LAM ** 3 * raw[:, 0] / 1e3)


def spec(t, x):
    """Segment PSDs of the whole record: centres (s), f, PSD per segment (rows = f), and hop/fs."""
    fs = 1 / np.median(np.diff(t))
    n = int(TSEG * fs)
    n -= n % 4
    f, tc, S = spectrogram(x, fs=fs, window="hann", nperseg=n, noverlap=3 * n // 4, detrend="constant",
                           scaling="density", mode="psd")
    return t[0] + tc, f, S, (n // 4) / fs


prs.configure_font()
d = load(TID)
t_end = min(d["tf"][-1], d["tu"][-1])
w = (d["tf"] > 330) & (d["tf"] < 370)
t_sol = d["tf"][w][np.argmax(np.abs(d["F"][w]))]
PH = [("earthquake", 0.0, 140.0), ("free vibration", 140.0, t_sol - 25), ("solitary wave and after", t_sol - 25, t_end)]
print(f"{TID}: solitary peak at {t_sol:.1f} s; phases " + ", ".join(f"{n} {a:.0f}-{b:.0f} s" for n, a, b in PH))

fig, A = plt.subplots(2, 2, figsize=(6.0, 5.4), layout="constrained", sharex=True)
for c, (kt, kx, unit) in enumerate((("tf", "F", "kN"), ("tu", "u", "mm"))):
    k = d[kt] < t_end
    tc, f, S, dth = spec(d[kt][k], d[kx][k])
    E = S.sum(axis=1) * dth
    T = tc.size * dth
    A[0, c].loglog(f[1:], E[1:] / T, color="k", lw=1.2, label="whole record", zorder=3)
    A[1, c].loglog(f[1:], E[1:], color="k", lw=1.2, label="whole record", zorder=3)
    tot = np.zeros_like(E)
    for (name, a, b), col in zip(PH, PHASE_COL):
        sel = (tc >= a) & (tc < b)
        Ei, Ti = S[:, sel].sum(axis=1) * dth, sel.sum() * dth
        A[0, c].loglog(f[1:], Ei[1:] / Ti, color=col, lw=0.7, label=name)
        A[1, c].loglog(f[1:], Ei[1:], color=col, lw=0.7, label=name)
        tot += Ei
        print(f"  {kx} {name:24s}: T = {Ti:5.1f} s, energy = {np.trapz(Ei, f):10.4g} {unit}^2 s")
    A[1, c].loglog(f[1:], tot[1:], color="#E08214", lw=0.8, ls=(0, (2, 1.5)), label="sum of the phases", zorder=4)
    print(f"  {kx} whole record: energy = {np.trapz(E, f):10.4g}, sum of phases {np.trapz(tot, f):10.4g} {unit}^2 s")
    sym = r"$\lambda_L^3 f_p$" if kx == "F" else r"$u$"
    A[0, c].set_ylabel(f"PSD of {sym} ({unit}$^2$/Hz)")
    A[1, c].set_ylabel(f"ESD of {sym} ({unit}$^2$" + r"$\cdot$s/Hz)")
for ax in A.flat:
    ax.set_xlim(0.025, 100)
    for x0 in (1 / 12.0, 0.477, 2.975, 16.75, 60.0 / R.SQL):
        ax.axvline(x0, color="0.35", lw=0.5, ls=(0, (3, 2)), zorder=0)
for ax in A[1]:
    ax.set_xlabel("Real-time frequency, prototype (Hz)")
for ax, s in zip(A.flat, ("(a) force, PSD", "(b) displacement, PSD", "(c) force, ESD", "(d) displacement, ESD")):
    ax.set_title(rf"\textbf{{{s.split()[0]}}} {' '.join(s.split()[1:])}", loc="left")
A[0, 0].legend(loc="lower left", frameon=False, fontsize=7)
A[1, 0].legend(loc="lower left", frameon=False, fontsize=7)
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
