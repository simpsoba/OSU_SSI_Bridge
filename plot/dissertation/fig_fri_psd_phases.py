"""PSDs of the interface force and displacement over equal phase windows, F03 (dissertation Ch. 6 draft; read-only).

F03: CPU MKR-alpha + ProfileSPD, baseline model, 8 ranks, one slowdown (laboratory time ~ integration time), so no
GPU-integrator displacement enters. Signals as fig_wed_psd.py (prototype-scaled clock; lambda_L^3 f_p in kN, u in mm).
Three contiguous 150-s windows: earthquake (0-150 s; 5-95 % Arias 56-139 s), long-period tail of the ground motion
(150-300 s; the FKSH19 NS record ends at 300 s), free vibration and solitary wave (300-450 s; peak near 349 s); labels: earthquake, coda, solitary wave. Welch with the same Hann segments (40 s, 75 % overlap, segment means
removed) in every window, so the curves share resolution and averaging; whole record (light grey) for reference.
Direct labels instead of a legend; uniform 9 pt font.
Reference lines as fig_wed_psd.py: f_1, f_4 (fundamental and pier modes of the numerical model, eigen mesh0/mesh19), f_osc (nominal load-transfer band 15.5-18.1 Hz, Seki et al. 2026), f_mains; dotted:
first mode during strong shaking (0.336 Hz = 0.70 f_1, minimum of the 40-s sliding dominant frequency); shaded:
0.1-0.3 Hz long-period content of the FKSH19 NS ground motion.
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
TID = "F03"
TSEG = 40.0
PH = [("earthquake", 0.0, 150.0, "0.0"), ("coda", 150.0, 300.0, "#08519C"), ("solitary wave", 300.0, 450.0, "#C0390B")]
# direct labels (panel, x in Hz, y, ha) per curve; the colours are the same in both panels
LAB = {"earthquake": (1, 0.03, 4e4, "left"), "coda": (1, 0.03, 4.0, "left"), "solitary wave": (1, 0.03, 2.5e-3, "left"),
       "whole record": (1, 3.5, 3.0, "left")}  # with a leader to the curve at 1.3 Hz
GREY = "0.72"


def load(tid):
    stem = R.BY_ID[tid][1]
    me = np.load(R.MAT_EXTRACT_DIR / f"{stem}.npz", allow_pickle=True)["meaSigOS_data"]
    tf, fr = R.load_force(stem)  # cached load-cell force (cache_force.py)
    return dict(tu=me[:, -1] * R.SQL, u=me[:, 0] * R.LAM * 1e3, tf=tf * R.SQL, F=-R.LAM ** 3 * fr / 1e3)


def psd(t, x, w0, w1):
    m = (t >= w0) & (t < w1)
    fs = 1 / np.median(np.diff(t[m]))
    n = int(TSEG * fs)
    n -= n % 4
    return welch(x[m], fs=fs, window="hann", nperseg=n, noverlap=3 * n // 4, detrend="constant")

F1_SOFT = 0.336 / 0.477  # dominant frequency during strong shaking (40-s windows, 80-120 s) over f_1
REFS = [(0.336, rf"${F1_SOFT:.2f}f_1$"), (0.477, r"$f_1$"), (2.975, r"$f_4$"), (16.75, r"$f_{\mathrm{osc}}$"), (60.0 / R.SQL, r"$f_{\mathrm{mains}}$")]


def ref_labels(ax):
    ax.axvspan(15.5, 18.1, color="0.9", lw=0, zorder=-1)  # nominal load-transfer system band (Seki et al. 2026)
    for x0, _ in REFS:  # dotted: the softened first mode during strong shaking; dashed: nominal frequencies
        ax.axvline(x0, color="0.35", lw=0.6 if x0 < 0.4 else 0.5, ls=(0, (1, 1.5)) if x0 < 0.4 else (0, (3, 2)), zorder=0)
    ax.axvspan(0.1, 0.3, color="#DCE6F0", lw=0, zorder=-1)  # long-period ground motion (FKSH19 NS record coda)
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



prs.configure_font()
d = load(TID)
t_end = min(d["tf"][-1], d["tu"][-1])
fig, A = plt.subplots(1, 2, figsize=(6.0, 2.9), layout="constrained")
for c, (kt, kx, unit) in enumerate((("tf", "F", "kN"), ("tu", "u", "mm"))):
    ax = A[c]
    f, P = psd(d[kt], d[kx], 0.0, t_end)
    ax.loglog(f[1:], P[1:], color=GREY, lw=1.6, zorder=1)
    f_all, P_all = f, P
    for name, a, b, col in PH:
        f, P = psd(d[kt], d[kx], a, b)
        ax.loglog(f[1:], P[1:], color=col, lw=0.8, zorder=3)
        k = (f > 0.4) & (f < 0.6)
        print(f"{kx} {name:24s}: mean square {np.trapz(P, f):10.4g} {unit}^2, peak near f1 at {f[k][np.argmax(P[k])]:.3f} Hz")
    ax.set_xlim(0.025, 100)
    ref_labels(ax)
    for name, (pc, x, y, ha) in LAB.items():
        if pc == c:
            col = GREY if name == "whole record" else next(cc for nn, _, _, cc in PH if nn == name)
            if name == "whole record":
                ax.annotate(name, (1.3, np.interp(1.3, f_all, P_all)), xytext=(x, y), ha=ha, va="center", color="0.45",
                            arrowprops=dict(arrowstyle="-", color="0.55", lw=0.6, shrinkA=1, shrinkB=1))
            else:
                ax.text(x, y, name, ha=ha, va="center", color=col)
    if kx == "u":
        ax.text(np.sqrt(0.1 * 0.3), 0.03, "long-period\nband", transform=ax.get_xaxis_transform(), ha="center",
                va="bottom", color="#3B5B7A", linespacing=0.95)
    ax.set_xlabel("Real-time frequency, prototype (Hz)")
    sym = r"$\lambda_{\mathrm{L}}^3 f_{\mathrm{p}}$" if kx == "F" else r"$u$"
    ax.set_ylabel(f"PSD of {sym} ({unit}$^2$/Hz)")
    ax.set_title(rf"\textbf{{({'ab'[c]})}} {'force' if kx == 'F' else 'displacement'}", loc="left", pad=42)
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
