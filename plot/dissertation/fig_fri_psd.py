"""PSDs of the interface force and displacement, CPU-integrator runs only (dissertation Ch. 6 draft; read-only on data).

Same signals and Welch settings as fig_wed_psd.py (Simulink at 2048 Hz model, laboratory clock scaled to prototype;
lambda_L^3 f_p = -lambda^3 topLRaw in kN, u = meaSigOS lambda in mm; Hann, 50 % overlap, 3 segments per window).
Only runs integrated on the CPU (CPU MKR-alpha), so no GPU-integrator displacement enters the comparison.
(a), (b) wave loading, equal 70-s windows: W02 (Wednesday, regular waves, 8 ranks, 70-140 s, before the earthquake)
         and F03 (Friday, baseline, 8 ranks): solitary wave (330-400 s) and still water after the earthquake (250-320 s)
(c), (d) earthquake, Friday CPU-integrator runs over the interval common to all of them (from 0 to the end of the
         shortest record): F03 and F26 (baseline), F08 and F07 (moderate, NC), F15 (large)
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy.signal import welch

sys.path.insert(0, str(Path(__file__).parent))
import fri_runs as R  # noqa: E402

sys.path.insert(0, str(R.ROOT / "plot"))
import PlotResponseSpectrum as prs  # noqa: E402

OUT = Path(sys.argv[1])
MATD = R.DATA / "Simulink"
LAM, SQL = R.LAM, R.SQL
W02 = "0819_GusBridge_Parallel8_H0p5_T7p746_reg"
EQ = [("F03", "#C0390B", "-", "F03, ProfileSPD, baseline, 1 slowdown"),
      ("F26", "#8C510A", (0, (4, 1.5)), "F26, MUMPS, baseline, 22"),
      ("F15", "#2166AC", "-", "F15, cuDSS solves, large, 566"),
      ("F08", "#8C510A", "-", "F08, MUMPS, moderate, 475, NC"),
      ("F07", "#E8A090", "-", "F07, ProfileSPD, moderate, 2,216, NC")]
XLIM = (0.02, 100.0)


def load(mat):
    m = sio.loadmat(MATD / f"{mat}.mat", squeeze_me=True, struct_as_record=False)["data"]
    me, raw = m.meaSigOS.data, m.topLRaw.data
    return dict(tu=me[:, -1] * SQL, u=me[:, 0] * LAM * 1e3, tf=raw[:, -1] * SQL, F=-LAM ** 3 * raw[:, 0] / 1e3)


def psd(t, x, w0, w1, nseg=3):
    m = (t >= w0) & (t < w1)
    y = x[m] - x[m].mean()
    nper = int(2 * y.size // (nseg + 1))
    nper -= nper % 2
    return welch(y, fs=1 / np.median(np.diff(t[m])), window="hann", nperseg=nper, noverlap=nper // 2)


prs.configure_font()
fig, A = plt.subplots(2, 2, figsize=(6.0, 5.2), layout="constrained", sharex=True)
wd = load(W02)
f3 = load(R.BY_ID["F03"][1])
WAVES = [(wd, (70, 140), "#08519C", "-", "W02, regular waves (Wed.)"),
         (f3, (330, 400), "#000000", "-", "F03, solitary wave"),
         (f3, (250, 320), "0.6", "-", "F03, still water after the earthquake")]
for col, (kt, kx) in enumerate((("tf", "F"), ("tu", "u"))):
    ax = A[0, col]
    for d, (w0, w1), c, ls, lab in WAVES:
        f, P = psd(d[kt], d[kx], w0, w1)
        ax.loglog(f[1:], P[1:], color=c, ls=ls, lw=0.7, label=lab)
EQD = {tid: load(R.BY_ID[tid][1]) for tid, *_ in EQ}
t_end = min(R.load(tid)["t"][-1] * SQL for tid, *_ in EQ)
print(f"earthquake window 0-{t_end:.1f} s (prototype laboratory time)")
for col, (kt, kx) in enumerate((("tf", "F"), ("tu", "u"))):
    ax = A[1, col]
    for tid, c, ls, lab in EQ:
        f, P = psd(EQD[tid][kt], EQD[tid][kx], 0, t_end)
        ax.loglog(f[1:], P[1:], color=c, ls=ls, lw=0.7, label=lab, zorder=3 if tid == "F03" else 2)
for ax in A.flat:
    ax.set_xlim(*XLIM)
    for x0 in (1 / 12.0, 0.477, 2.870, 12.0, 60.0 / SQL):
        ax.axvline(x0, color="0.35", lw=0.5, ls=(0, (3, 2)), zorder=0)
for r in range(2):
    A[r, 0].set_ylabel(r"PSD of $\lambda_L^3 f_p$ (kN$^2$/Hz)")
    A[r, 1].set_ylabel(r"PSD of $u$ (mm$^2$/Hz)")
for ax in A[1]:
    ax.set_xlabel("Real-time frequency, prototype (Hz)")
for ax, s in zip(A.flat, ("(a) waves", "(b) waves", f"(c) earthquake, 0-{t_end:.0f} s", "(d) earthquake")):
    ax.set_title(rf"\textbf{{{s.split()[0]}}} {' '.join(s.split()[1:])}", loc="left")
A[0, 0].legend(loc="lower left", frameon=False, fontsize=7)
A[1, 0].legend(loc="lower left", frameon=False, fontsize=7)
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
