"""Interface force, Wednesday runs, all in real time (dissertation Ch. 6; read-only on data).

lambda_L^3 f_p = -lambda^3 * topLRaw (Simulink load cell, 2048 Hz model; verified against OpenFresco daqFrc,
gain -13.83, corr 0.9999 below 2 Hz), resisting-force convention, on the laboratory clock sqrt(lambda_L) t.
(a)        slowdown lane (one row per run, test ID in the free space) and the force history;
           time axis broken at 170 s, force axis broken at +-12 / +-15 kN
(b.1)-(b.3) zooms with a slowdown lane each: waves only, 24-rank slowdowns, 8-rank slowdown;
           multi-cycle period measurements between successive crests of the 8-rank trace
Records end where the 1-rank actuator stops tracking. Bold panel tags, no titles or legends, uniform font.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from matplotlib.patches import Rectangle
from scipy.signal import butter, filtfilt, find_peaks

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
ORDER = [0, 1, 3, 4, 2]
SEGX = [(0, 170), (170, 420)]
WRX = [1, 3.4]
BANDS = [(15, 110), (-12, 12), (-110, -15)]
HRY = [0.6, 1.7, 0.6]
ZOOMS = [((136, 168), (-5.5, 6.0), [1, 3, 4, 2]),
         ((216.2, 217.4), (-12, 15), [4, 2]),
         ((262.9, 263.9), (-15, 40), [4, 2])]


def fy_eq_kn():
    """Equivalent pier-top force at cantilever yield, F_y,eq = Mn / H_pier, Mn = (2/pi) Ast fy r
    (same estimate as plot/PlotActuatorForce.py: 28 #10 bars, fy = 470 MPa, 4 ft pier)."""
    inch = 0.0254
    d_pier, h_pier = 48 * inch, 3.02 * (48.0 / 20.0)
    as_tot = 28 * (np.pi / 4) * (1.25 * inch) ** 2
    r_bar = 0.5 * d_pier - 2 * inch - 0.5 * 1.25 * inch - 0.875 * inch
    return (2 / np.pi) * as_tot * 470e6 * r_bar / h_pier / 1e3


F_REF = 0.1 * fy_eq_kn()


def tag(ax, s):
    ax.set_title(rf"\textbf{{{s}}}", loc="left")


def load(mat):
    m = sio.loadmat(MATD / f"{mat}.mat", squeeze_me=True, struct_as_record=False)["data"]
    st, me, tar, raw = m.stateOS.data, m.meaSigOS.data, m.tarSigOS.data, m.topLRaw.data
    tl, ty = st[:, -1], np.rint(st[:, 0]).astype(int)
    on = tl[1:][(ty[1:] == 2) & (ty[:-1] != 2)]
    err = tar[:, 0] - me[:, 0]
    rms = np.sqrt(np.convolve(err ** 2, np.ones(512) / 512, mode="same"))
    bad = np.flatnonzero((rms > 2e-3) & (me[:, -1] > 50.0))
    t_end = me[max(bad[0] - 256, 0), -1] if bad.size else me[-1, -1]
    k = raw[:, -1] <= t_end
    return dict(t=raw[k, -1] * SQL, F=-LAM ** 3 * raw[k, 0] / 1e3, on=on[on > 0.7] * SQL, t_end=t_end * SQL)


def mark(ax, x, y, dx, dy):
    ax.plot([x - dx, x + dx], [y - dy, y + dy], transform=ax.transAxes, color="k", clip_on=False, lw=0.8)


def label(ax, text, xy, xytext, color, **kw):
    ax.annotate(text, xy, xytext=xytext, color=color, ha=kw.pop("ha", "left"), va=kw.pop("va", "center"),
                arrowprops=dict(arrowstyle="-", color=color, lw=0.6, shrinkA=1, shrinkB=1), **kw)


def crest_span(t, F, lo, hi, ncyc, lowpass=None):
    """First ncyc cycles between successive crests of F inside [lo, hi]."""
    dt = np.median(np.diff(t))
    x = F
    if lowpass:
        b, a = butter(2, lowpass * 2 * dt, "low")
        x = filtfilt(b, a, F)
    m = (t >= lo) & (t <= hi)
    dist = int((8.0 if lowpass else 0.05) / dt)
    ip, _ = find_peaks(x[m], prominence=0.5 if lowpass else 0.8, distance=dist)
    tp = t[m][ip]
    return tp[0], tp[ncyc]


prs.configure_font()
data = [load(r[1]) for r in RUNS]
dW = data[2]
MEAS = [crest_span(dW["t"], dW["F"], 140, 168, 2, lowpass=0.5),
        crest_span(dW["t"], dW["F"], 216.29, 216.7, 4),
        crest_span(dW["t"], dW["F"], 263.42, 263.9, 4)]
for (t0, t1), nc in zip(MEAS, (2, 4, 4)):
    print(f"measured {nc} cycles: {t0:.3f}-{t1:.3f} s -> {nc / (t1 - t0):.3f} Hz")

fig = plt.figure(figsize=(6.0, 6.3), layout="constrained")
sa, sb = fig.subfigures(2, 1, height_ratios=[2.3, 1.2])

# ---------------- (a) lane + history ----------------
g = sa.add_gridspec(4, 2, height_ratios=[0.95] + HRY, width_ratios=WRX, hspace=0.06, wspace=0.03)
LANE = [sa.add_subplot(g[0, j]) for j in range(2)]
AX = [[sa.add_subplot(g[i + 1, j]) for j in range(2)] for i in range(3)]
for i in range(3):
    for j in range(2):
        ax = AX[i][j]
        ax.set_xlim(*SEGX[j])
        ax.set_ylim(*BANDS[i])
        for k in ORDER:
            ax.plot(data[k]["t"], data[k]["F"], color=RUNS[k][3], lw=0.6 if RUNS[k][2] == 1 else 0.7)
        if i < 2:
            ax.spines["bottom"].set_visible(False)
            ax.tick_params(bottom=False, labelbottom=False)
        if j == 1:
            ax.spines["left"].set_visible(False)
            ax.tick_params(left=False, labelleft=False)
for i, sgn in ((0, 1), (2, -1)):  # 0.1 F_y,eq reference in the outer force bands
    for j in range(2):
        AX[i][j].axhline(sgn * F_REF, color="0.25", lw=0.7, ls=(0, (5, 3)), zorder=0)
AX[0][1].text(418, F_REF + 3, r"$0.1F_{y,\mathrm{eq}}$", ha="right", va="bottom", color="0.25")
AX[0][0].set_yticks([50, 100])
AX[1][0].set_yticks([-10, 0, 10])
AX[2][0].set_yticks([-100, -50])
AX[2][0].set_xticks([0, 100])
for i, ys in enumerate(((0,), (0, 1), (1,))):
    for y in ys:
        mark(AX[i][0], 0.0, y, 0.025, 0.05 * HRY[1] / HRY[i])
for j, x in ((0, 1.0), (1, 0.0)):
    mark(AX[2][j], x, 0.0, 0.036 * (1 if j == 0 else WRX[0] / WRX[1]), 0.12)
AX[1][0].set_ylabel(r"$\lambda_L^3 f_p$ (kN)")
AX[2][1].set_xlabel(r"Laboratory time, $\sqrt{\lambda_L}\,t$ (s)")
label(AX[0][1], "1 rank", (250, 70), (290, 80), ORANGE)
label(AX[1][1], "8 ranks", (330, 3.0), (340, 9.5), BLUE)
label(AX[1][1], "4, 16, 24 ranks", (395, -3.5), (355, -9.5), GRAY)
label(AX[1][1], "1 rank: test stopped", (data[0]["t_end"], -1.0), (290, -7.5), ORANGE)
for j in range(2):
    ax = LANE[j]
    lo, hi = SEGX[j]
    for r in range(5):
        n = RUNS[r][2]
        y = 4 - r
        tend = data[r]["t_end"]
        ax.add_patch(Rectangle((0, y - 0.32), tend, 0.64, fc="0.9", ec="none", zorder=1))
        on = data[r]["on"]
        on = on[(on >= lo) & (on <= min(hi, tend))]
        ax.vlines(on, y - 0.36, y + 0.36, color={1: ORANGE, 8: BLUE}.get(n, "0.15"), lw=1.0 if n != 1 else 0.4, zorder=2)
        if j == 1:
            ax.text(hi - 3, y, RUNS[r][0], ha="right", va="center", color="0.3", zorder=3)
    ax.set_xlim(lo, hi)
    ax.set_ylim(-0.6, 4.6)
    ax.grid(False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(bottom=False, labelbottom=False)
    if j == 1:
        ax.spines["left"].set_visible(False)
        ax.tick_params(left=False, labelleft=False)
LANE[0].set_yticks(range(5))
LANE[0].set_yticklabels([str(RUNS[r][2]) for r in range(5)][::-1])
LANE[0].set_ylabel("Slowdowns\n(MPI ranks)")
tag(LANE[0], "(a)")

# ---------------- (b) zooms, real time ----------------
gz = sb.add_gridspec(2, 3, height_ratios=[0.16, 1.0], hspace=0.04)
Z = [sb.add_subplot(gz[1, j]) for j in range(3)]
ZL = [sb.add_subplot(gz[0, j], sharex=Z[j]) for j in range(3)]
for j, ((z0, z1), yl, ks) in enumerate(ZOOMS):
    ax, lz = Z[j], ZL[j]
    for k in ks:
        n, c = RUNS[k][2], RUNS[k][3]
        d = data[k]
        m = (d["t"] >= z0 - 1) & (d["t"] <= z1 + 1)
        ax.plot(d["t"][m], d["F"][m], color=c, lw=1.0)
        on = d["on"][(d["on"] >= z0) & (d["on"] <= z1)]
        lz.vlines(on, 0, 1, color=BLUE if n == 8 else "0.15", lw=1.0)
    lz.add_patch(Rectangle((z0, 0.1), z1 - z0, 0.8, fc="0.9", ec="none", zorder=0))
    lz.set_ylim(0, 1)
    lz.grid(False)
    lz.set_yticks([])
    for s in ("left", "bottom"):
        lz.spines[s].set_visible(False)
    lz.tick_params(bottom=False, labelbottom=False)
    tag(lz, f"(b.{j + 1})")
    ax.set_xlim(z0, z1)
    ax.set_ylim(*yl)
    ax.set_xlabel(r"$\sqrt{\lambda_L}\,t$ (s)")
Z[0].set_ylabel(r"$\lambda_L^3 f_p$ (kN)")
Z[1].set_xticks([216.5, 217.0])
Z[1].set_yticks([-10, 0, 10])
Z[2].set_xticks([263.0, 263.5])
for j, ((t0, t1), nc) in enumerate(zip(MEAS, (2, 4, 4))):
    ax = Z[j]
    yl = ZOOMS[j][1]
    (z0_, z1_) = ZOOMS[j][0]
    yv = lambda tt: np.interp(tt, dW["t"], dW["F"])
    ya = max(yv(t0), yv(t1)) + 0.12 * (yl[1] - yl[0])
    for tt in (t0, t1):
        ax.plot([tt, tt], [yv(tt), ya], color="0.3", lw=0.5)
    ax.annotate("", (t0, ya), (t1, ya), arrowprops=dict(arrowstyle="<->", color="0.3", lw=0.7, shrinkA=0, shrinkB=0))
    lab = rf"${nc}/f_w$ = {t1 - t0:.1f} s" if j == 0 else rf"${nc}/f_{{\mathrm{{osc}}}}$ = {t1 - t0:.3f} s"
    xc = min(max(0.5 * (t0 + t1), z0_ + 0.24 * (z1_ - z0_)), z1_ - 0.24 * (z1_ - z0_))
    ax.text(xc, ya + 0.02 * (yl[1] - yl[0]), lab, ha="center", va="bottom", color="0.2", zorder=6,
            bbox=dict(fc="white", ec="none", alpha=0.9, pad=0.5))
label(Z[0], "8 ranks", (154.3, -2.7), (156.5, -4.6), BLUE)
label(Z[1], "24 ranks", (217.17, -9.0), (216.62, -9.5), GRAY)
label(Z[1], "8 ranks", (217.32, 1.6), (217.0, 12.0), BLUE)
label(Z[2], "8 ranks", (263.27, 30), (263.38, 26), BLUE)

fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
