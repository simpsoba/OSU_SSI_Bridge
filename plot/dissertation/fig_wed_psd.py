"""Dissertation Ch. 6 (read-only on data): PSDs of the interface force and displacement, real time.

Simulink signals at 2048 Hz model, laboratory clock scaled to prototype (sqrt(lambda_L) t), so no aliasing:
  force        lambda_L^3 f_p = -lambda^3 * topLRaw (kN; verified against daqFrc: gain -13.83, corr 0.9999 < 2 Hz)
  displacement u = meaSigOS * lambda (mm)
Records end where the actuator stops tracking (1 rank). One window for every run: from 0 to where the 1-rank actuator stops tracking (272.4 s, real time). PSD: Welch, Hann segments, 50% overlap, even segment length chosen so the
segments tile the window (3 segments); one-sided; area = mean square.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy.signal import welch

ROOT = Path(__file__).resolve().parents[2]  # OSU_SSI_Bridge
sys.path.insert(0, str(ROOT / "plot"))
import PlotResponseSpectrum as prs  # noqa: E402
from lab_paths import CYLINDER_LENGTH_SCALE as LAM, TIME_SCALE_FROUDE as SQL  # noqa: E402

OUT = Path(sys.argv[1])
MATD = ROOT / "OSU_SSI_BRIDGE_DATA" / "Simulink"
ORANGE, BLUE, GRAY = "#C0390B", "#000000", "#8C8C8C"  # 1 rank deep vermillion, 8 ranks black, others gray
RUNS = [("W04", "0819_GusBridge_Serial_H0p5_T7p746_reg_Trial02", 1, ORANGE),
        ("W06", "0819_GusBridge_Parallel4_H0p5_T7p746_reg_Trial02", 4, GRAY),
        ("W02", "0819_GusBridge_Parallel8_H0p5_T7p746_reg", 8, BLUE),
        ("W05", "0819_GusBridge_Parallel16_H0p5_T7p746_reg_Trial02", 16, GRAY),
        ("W07", "0819_GusBridge_Parallel24_H0p5_T7p746_reg", 24, GRAY)]
F_WAVE = 1 / 12.0
REFS = [(F_WAVE, r"$f_w$"), (2 * F_WAVE, r"$2f_w$"), (3 * F_WAVE, r"$3f_w$"), (4 * F_WAVE, r"$4f_w$"),
        (5 * F_WAVE, r"$5f_w$"), (0.477, r"$f_1$"), (2.870, r"$f_2$"), (12.0, r"$f_{\mathrm{osc}}$"),
        (60.0 / SQL, r"$f_{\mathrm{mains}}$")]
XLIM = (0.02, 100.0)


def load(mat):
    m = sio.loadmat(MATD / f"{mat}.mat", squeeze_me=True, struct_as_record=False)["data"]
    me, tar, raw = m.meaSigOS.data, m.tarSigOS.data, m.topLRaw.data
    err = tar[:, 0] - me[:, 0]
    rms = np.sqrt(np.convolve(err ** 2, np.ones(512) / 512, mode="same"))
    bad = np.flatnonzero((rms > 2e-3) & (me[:, -1] > 50.0))
    t_end = me[max(bad[0] - 256, 0), -1] if bad.size else me[-1, -1]
    km, kr = me[:, -1] <= t_end, raw[:, -1] <= t_end
    return dict(tu=me[km, -1] * SQL, u=me[km, 0] * LAM * 1e3, tf=raw[kr, -1] * SQL, F=-LAM ** 3 * raw[kr, 0] / 1e3,
                t_end=t_end * SQL)


def psd(t, x, w1, nseg):
    m = t < w1
    y = x[m] - x[m].mean()
    nper = int(2 * y.size // (nseg + 1))
    nper -= nper % 2
    assert (nseg - 1) * (nper // 2) + nper <= y.size
    return welch(y, fs=1 / np.median(np.diff(t[m])), window="hann", nperseg=nper, noverlap=nper // 2)


def ref_labels(ax):
    for x0, _ in REFS:
        ax.axvline(x0, color="0.35", lw=0.5, ls=(0, (3, 2)), zorder=0)
    items = sorted(REFS, key=lambda it: it[0])
    lx0, lx1 = np.log10(ax.get_xlim())
    xf = [(np.log10(x) - lx0) / (lx1 - lx0) for x, _ in items]
    gap, xl = 0.062, list(xf)
    for i in range(1, len(xl)):
        xl[i] = max(xl[i], xl[i - 1] + gap)
    over = xl[-1] - 0.99
    if over > 0:
        xl = [v - over for v in xl]
        for i in range(len(xl) - 2, -1, -1):
            xl[i] = min(xl[i], xl[i + 1] - gap)
    for (x, lab), a_, b_ in zip(items, xf, xl):
        ax.plot([a_, a_, b_], [1.0, 1.025, 1.07], transform=ax.transAxes, color="0.45", lw=0.5, clip_on=False)
        ax.text(b_, 1.085, lab, transform=ax.transAxes, rotation=90, ha="center", va="bottom", color="0.25")


def label(ax, text, xy, xytext, color):
    ax.annotate(text, xy, xytext=xytext, color=color, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=color, lw=0.6, shrinkA=1, shrinkB=1))


prs.configure_font()
D = [load(r[1]) for r in RUNS]
w_all = D[0]["t_end"]
print(f"window: 0-{w_all:.1f} s (prototype real time)")


def tag(ax, s):
    ax.set_title(rf"\textbf{{{s}}}", loc="left", pad=34)


fig, A = plt.subplots(1, 2, figsize=(6.0, 2.7), layout="constrained")
for row, (kt, kx, ylab, ylim) in enumerate((("tf", "F", r"PSD of $\lambda_L^3 f_p$ (kN$^2$/Hz)", (1e-8, 1e2)),
                                             ("tu", "u", r"PSD of $u$ (mm$^2$/Hz)", (1e-9, 1e5)))):
    ax = A[row]
    for k in [0, 1, 3, 4, 2]:
        n, c = RUNS[k][2], RUNS[k][3]
        f, P = psd(D[k][kt], D[k][kx], w_all, 3)
        ax.loglog(f[1:], P[1:], color=c, lw=0.6, alpha=0.8 if c == GRAY else 1.0, zorder={1: 1, 8: 3}.get(n, 2))
    ax.set_xlim(*XLIM)
    ax.set_ylim(*ylim)
    ax.set_ylabel(ylab)
    tag(ax, f"({'ab'[row]})")
    ref_labels(ax)
for ax in A:
    ax.set_xlabel("Real-time frequency, prototype (Hz)")
label(A[0], "1 rank", (25, 3e-1), (1.2, 3e1), ORANGE)
label(A[0], "8 ranks", (40, 1e-5), (1.0, 1e-7), BLUE)
label(A[0], "4, 16, 24 ranks", (0.04, 6e-3), (0.023, 1e-6), GRAY)
label(A[1], "1 rank", (0.085, 0.08), (0.025, 3e-4), ORANGE)
label(A[1], "8 ranks", (0.75, 1.5), (1.5, 3e2), BLUE)
label(A[1], "4, 16, 24 ranks", (40.0, 2e-7), (0.8, 1e-8), GRAY)
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
