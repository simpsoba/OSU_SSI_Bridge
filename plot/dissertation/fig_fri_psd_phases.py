"""PSDs of the interface force and displacement over equal phase windows, F03 (dissertation Ch. 6 draft; read-only).

F03: CPU MKR-alpha + ProfileSPD, baseline model, 8 ranks, one slowdown (laboratory time ~ integration time), so no
GPU-integrator displacement enters. Signals as fig_wed_psd.py (prototype-scaled clock; lambda_L^3 f_p in kN, u in mm).
Three 120-s windows: earthquake (20-140 s; 5-95 % Arias 56-139 s), free vibration in still water (200-320 s), solitary
wave and after (330-450 s; peak near 349 s). Welch with the same Hann segments (40 s, 75 % overlap, segment means
removed) in every window, so the curves share resolution and averaging; whole record (faint) for reference.
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
TID = "F03"
TSEG = 40.0
PH = [("earthquake", 20.0, 140.0, "0.0"), ("free vibration", 200.0, 320.0, "#08519C"),
      ("solitary wave and after", 330.0, 450.0, "#C0390B")]


def load(tid):
    m = sio.loadmat(R.DATA / "Simulink" / f"{R.BY_ID[tid][1]}.mat", squeeze_me=True, struct_as_record=False)["data"]
    me, raw = m.meaSigOS.data, m.topLRaw.data
    return dict(tu=me[:, -1] * R.SQL, u=me[:, 0] * R.LAM * 1e3, tf=raw[:, -1] * R.SQL, F=-R.LAM ** 3 * raw[:, 0] / 1e3)


def psd(t, x, w0, w1):
    m = (t >= w0) & (t < w1)
    fs = 1 / np.median(np.diff(t[m]))
    n = int(TSEG * fs)
    n -= n % 4
    return welch(x[m], fs=fs, window="hann", nperseg=n, noverlap=3 * n // 4, detrend="constant")


prs.configure_font()
d = load(TID)
t_end = min(d["tf"][-1], d["tu"][-1])
fig, A = plt.subplots(1, 2, figsize=(6.0, 2.9), layout="constrained")
for c, (kt, kx, unit) in enumerate((("tf", "F", "kN"), ("tu", "u", "mm"))):
    ax = A[c]
    f, P = psd(d[kt], d[kx], 0.0, t_end)
    ax.loglog(f[1:], P[1:], color="0.75", lw=1.6, label="whole record", zorder=1)
    for name, a, b, col in PH:
        f, P = psd(d[kt], d[kx], a, b)
        ax.loglog(f[1:], P[1:], color=col, lw=0.8, label=f"{name}, {a:.0f}-{b:.0f} s", zorder=3)
        k = (f > 0.4) & (f < 0.6)
        print(f"{kx} {name:24s}: mean square {np.trapz(P, f):10.4g} {unit}^2, peak near f1 at {f[k][np.argmax(P[k])]:.3f} Hz")
    ax.set_xlim(0.025, 100)
    for x0 in (1 / 12.0, 0.477, 2.870, 12.0, 60.0 / R.SQL):
        ax.axvline(x0, color="0.35", lw=0.5, ls=(0, (3, 2)), zorder=0)
    ax.set_xlabel("Real-time frequency, prototype (Hz)")
    sym = r"$\lambda_L^3 f_p$" if kx == "F" else r"$u$"
    ax.set_ylabel(f"PSD of {sym} ({unit}$^2$/Hz)")
    ax.set_title(rf"\textbf{{({'ab'[c]})}} {'force' if kx == 'F' else 'displacement'}", loc="left")
A[0].legend(loc="lower left", frameon=False, fontsize=7)
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
