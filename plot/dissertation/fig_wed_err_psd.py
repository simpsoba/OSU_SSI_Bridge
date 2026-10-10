"""PSD of the actuator control error, Wednesday runs (dissertation Ch. 6; read-only on data).

e = command - measured (comSigOS - meaSigOS; AFC-compensated command vs feedback), prototype mm (x lambda_L),
laboratory clock scaled to prototype real time (sqrt(lambda_L) t). Simulink 2048 Hz model.
(a) all runs, one window from 0 to where the 1-rank actuator stops tracking (272.4 s real time)
(b) 8 ranks by phase on integrator time: before EQ (< 150 s), during (150-289.1 s, to the Arias 95% time),
    after (> 289.1 s); each PSD of its own window (mean square per Hz, so durations are comparable)
(c) measured vs command, all runs, same window as (a), with the 1:1 line
Welch, Hann, 50% overlap, even segment length tiling the window (3 segments in (a), 5 per phase in (b)).
Bold panel tags, no titles or legends (direct labels), uniform font, constrained layout.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import welch

sys.path.insert(0, str(Path(__file__).parent))
import ctrl_error_wed as C  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]  # OSU_SSI_Bridge
sys.path.insert(0, str(ROOT / "plot"))
import PlotResponseSpectrum as prs  # noqa: E402

LAM, SQL = C.LAM, C.SQL
OUT = Path(sys.argv[1])
ORANGE, BLACK, GRAY = "#C0390B", "#000000", "#8C8C8C"
COL = {1: ORANGE, 4: GRAY, 8: BLACK, 16: GRAY, 24: GRAY}
F_WAVE = 1 / 12.0
REFS = [(F_WAVE, r"$f_w$"), (0.477, r"$f_1$"), (2.870, r"$f_2$"), (12.0, r"$f_{\mathrm{osc}}$"),
        (60.0 / SQL, r"$f_{\mathrm{mains}}$")]
XLIM = (0.02, 100.0)


def psd(t, x, nseg):
    y = x - x.mean()
    nper = int(2 * y.size // (nseg + 1))
    nper -= nper % 2
    assert (nseg - 1) * (nper // 2) + nper <= y.size
    return welch(y, fs=1 / np.median(np.diff(t)), window="hann", nperseg=nper, noverlap=nper // 2)


def ref_labels(ax):
    for x0, _ in REFS:
        ax.axvline(x0, color="0.35", lw=0.5, ls=(0, (3, 2)), zorder=0)
    lx0, lx1 = np.log10(ax.get_xlim())
    for x, lab in REFS:
        ax.text((np.log10(x) - lx0) / (lx1 - lx0), 1.02, lab, transform=ax.transAxes, rotation=90,
                ha="center", va="bottom", color="0.25")


def label(ax, text, xy, xytext, color):
    ax.annotate(text, xy, xytext=xytext, color=color, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=color, lw=0.6, shrinkA=1, shrinkB=1))


def tag(ax, s):
    ax.set_title(rf"\textbf{{{s}}}", loc="left", pad=24)


prs.configure_font()
D = {r[3]: C.load(r[1], r[2]) for r in C.RUNS}
w_all = D[1]["t"][-1] * SQL
print(f"window (a): 0-{w_all:.1f} s prototype real time")

fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(6.0, 2.5), layout="constrained",
                                       gridspec_kw=dict(width_ratios=[1, 1, 0.9]))
ax_b.sharey(ax_a)
ax_b.tick_params(labelleft=False)
for n in (4, 16, 24, 8, 1):
    d = D[n]
    m = d["t"] * SQL <= w_all
    f, P = psd(d["t"][m] * SQL, (d["co"][m] - d["me"][m]) * LAM * 1e3, 3)
    ax_a.loglog(f[1:], P[1:], color=COL[n], lw=0.6, zorder={1: 3, 8: 2}.get(n, 1))
    ax_c.plot(d["co"][m] * LAM * 1e3, d["me"][m] * LAM * 1e3, color=COL[n], lw=0.4,
              zorder={1: 3, 8: 2}.get(n, 1))
d = D[8]
for (ph, a, b), c in zip(C.PHASES, ("0.65", BLACK, "#5B8DB8")):
    m = (d["ti"] >= a) & (d["ti"] < b)
    f, P = psd(d["t"][m] * SQL, (d["co"][m] - d["me"][m]) * LAM * 1e3, 5)
    ax_b.loglog(f[1:], P[1:], color=c, lw=0.6)
    print(ph, "rms e (proto mm)", np.sqrt(np.trapz(P, f)))
for ax, tg in ((ax_a, "(a)"), (ax_b, "(b)")):
    ax.set_xlim(*XLIM)
    ax.set_ylim(1e-9, 1e4)
    ax.set_xlabel("Frequency, prototype (Hz)")
    tag(ax, tg)
    ref_labels(ax)
ax_a.set_ylabel(r"PSD of $u_{\mathrm{com}} - u_{\mathrm{mea}}$ (mm$^2$/Hz)")
label(ax_a, "1 rank", (20.0, 2e-3), (6.0, 2e-1), ORANGE)
label(ax_a, "8 ranks", (25.0, 1e-6), (1.0, 3e-8), BLACK)
label(ax_a, "4, 16, 24 ranks", (0.05, 1.5e-3), (0.023, 2e-6), GRAY)
label(ax_b, "during", (0.6, 3e-1), (1.5, 3e1), BLACK)
label(ax_b, "after", (3.0, 1.5e-4), (0.7, 3e-7), "#5B8DB8")
label(ax_b, "before", (0.15, 5e-4), (0.023, 1e-6), "0.55")
ax_b.text(0.97, 0.95, "8 ranks", transform=ax_b.transAxes, ha="right", va="top")
lim = (-200, 220)
ax_c.plot(lim, lim, color="0.35", lw=0.6, ls=(0, (4, 2)), zorder=0)
ax_c.text(150, 195, "1:1", ha="right", va="bottom", color="0.35")
ax_c.set_xlim(*lim)
ax_c.set_ylim(*lim)
ax_c.set_aspect("equal")
ax_c.set_xticks([-200, -100, 0, 100, 200])
ax_c.set_yticks([-200, -100, 0, 100, 200])
ax_c.set_xlabel(r"$u_{\mathrm{com}}$ (mm)")
ax_c.set_ylabel(r"$u_{\mathrm{mea}}$ (mm)")
tag(ax_c, "(c)")
for n in (1, 8):
    d = D[n]
    m = d["t"] * SQL <= w_all
    print(n, "u_com range (proto mm)", (d["co"][m] * LAM * 1e3).min(), (d["co"][m] * LAM * 1e3).max())
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
