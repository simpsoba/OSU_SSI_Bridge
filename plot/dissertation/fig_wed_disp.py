"""Interface displacement, Wednesday runs (dissertation Ch. 6; read-only on data).

u = meaSigOS x lambda (prototype mm), Simulink at 2048 Hz model. Records end where the serial actuator stops
tracking (0.25-s RMS of target - measured above 2 mm model). Integration time from the received-target count:
target k at t_int = (k + alpha_f) dt_int. Delay = laboratory time of each target minus k dt_sim (laboratory s).
(a)         slowdown lane (one row per run, test ID in the free space) and u on real time, sqrt(lambda_L) t
(b)         u on integration time; (b.1) zoom on the last cycles before the serial stop
(c)         delay behind real time after start-up, prototype-scaled laboratory s, symlog axis (as fig_fri_lanes.py)
Bold panel tags, no titles or legends (direct labels), uniform font, constrained layout.
"""
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[2]  # OSU_SSI_Bridge
sys.path.insert(0, str(ROOT / "plot"))
import PlotResponseSpectrum as prs  # noqa: E402
from lab_paths import CYLINDER_LENGTH_SCALE as LAM, MAT_EXTRACT_DIR, TIME_SCALE_FROUDE as SQL  # noqa: E402

OUT = Path(sys.argv[1])
DT_SIM = 10 / 2048
ALPHA_F = 1 / (1 + np.sqrt(0.5))
ORANGE, BLUE, GRAY = "#C0390B", "#000000", "#8C8C8C"  # serial deep vermillion, 8 ranks black, others gray
RUNS = [("W04", "0819_GusBridge_Serial_H0p5_T7p746_reg_Trial02", 1, ORANGE),
        ("W06", "0819_GusBridge_Parallel4_H0p5_T7p746_reg_Trial02", 4, GRAY),
        ("W02", "0819_GusBridge_Parallel8_H0p5_T7p746_reg", 8, BLUE),
        ("W05", "0819_GusBridge_Parallel16_H0p5_T7p746_reg_Trial02", 16, GRAY),
        ("W07", "0819_GusBridge_Parallel24_H0p5_T7p746_reg", 24, GRAY)]
WITH_STARTUP = os.environ.get("DELAY_STARTUP", "0") == "1"  # include the start-up offset in (c)
def fmt_delay(v):
    """Delay label at 0.1 ms resolution: 0 below 0.05 ms, 4 decimals below 1 s, else 0.1 s."""
    if v < 5e-5:
        return "0"
    return f"{v:.4f}" if v < 1.0 else f"{v:.1f}"

ORDER = [1, 3, 4, 2, 0]  # grays, 8 ranks, then serial on top
SEG = [(0, 170), (170, 420)]
WR = [1, 3.4]
ZOOMS = [(234, 245)]
YLIM = (-200, 200)


def tag(ax, s):
    ax.set_title(rf"\textbf{{{s}}}", loc="left")


def label(ax, text, xy, xytext, color, **kw):
    ax.annotate(text, xy, xytext=xytext, color=color, ha=kw.pop("ha", "left"), va=kw.pop("va", "center"),
                arrowprops=dict(arrowstyle="-", color=color, lw=0.6, shrinkA=1, shrinkB=1), **kw)


def mark(ax, x, y, dx, dy):
    ax.plot([x - dx, x + dx], [y - dy, y + dy], transform=ax.transAxes, color="k", clip_on=False, lw=0.8)


def load(stem):
    z = np.load(MAT_EXTRACT_DIR / f"{stem}.npz", allow_pickle=True)
    st, me = z["stateOS_data"], z["meaSigOS_data"]
    t, s1, typ = st[:, -1], np.rint(st[:, 2]), np.rint(st[:, 0]).astype(int)
    tk = t[1:][(s1[1:] > 0) & (s1[:-1] <= 0)]
    # delay behind real time from the rate-transition window starts (count = 1), which advance by exactly dt_sim
    # unless a slowdown stretches a window: free of the host-latency jitter carried by the target arrivals
    cnt = np.rint(st[:, 1]).astype(int)
    win = t[1:][(cnt[1:] == 1) & (cnt[:-1] != 1)]
    win = win[win >= tk[0]]
    lag = (win - win[0]) - DT_SIM * np.arange(win.size)
    k0 = 1 if tk[1] - tk[0] < 3e-3 else 0  # second start-up handshake flag
    if WITH_STARTUP:  # from the first target, start-up offset included
        su = win >= tk[k0]
    else:  # after the three start-up intervals, as fig_fri_lanes.py
        su = win >= tk[k0 + 3]
    win, lag = win[su], (lag[su] - lag[su][0]) * SQL  # delay, prototype-scaled laboratory s
    on = t[1:][(typ[1:] == 2) & (typ[:-1] != 2)]
    err = z["tarSigOS_data"][:, 0] - me[:, 0]
    rms = np.sqrt(np.convolve(err ** 2, np.ones(512) / 512, mode="same"))
    bad = np.flatnonzero((rms > 2e-3) & (me[:, -1] > 50.0))
    me = me[: (max(bad[0] - 256, 0) if bad.size else me.shape[0])]
    t_end = me[-1, -1]
    to_int = lambda x: (np.interp(x, tk, np.arange(tk.size)) + ALPHA_F) * DT_SIM * SQL
    on = on[(on > tk[k0 + 3]) & (on <= t_end)]  # after start-up
    kk = win <= t_end
    return dict(t=me[:, -1] * SQL, ti=to_int(me[:, -1]), u=me[:, 0] * LAM * 1e3, tk=win[kk] * SQL, lag=lag[kk],
                on=on * SQL, on_i=to_int(on), t_end=t_end * SQL, ti_end=to_int(t_end))


def broken_pair(sf, spec):
    sub = spec.subgridspec(1, 2, width_ratios=WR, wspace=0.03)
    al = sf.add_subplot(sub[0])
    ar = sf.add_subplot(sub[1], sharey=al)
    al.set_xlim(*SEG[0])
    ar.set_xlim(*SEG[1])
    ar.spines["left"].set_visible(False)
    ar.tick_params(left=False, labelleft=False)
    al.set_xticks([0, 100])
    for ax, x, w in ((al, 1.0, 0.036), (ar, 0.0, 0.036 * WR[0] / WR[1])):
        mark(ax, x, 0.0, w, 0.06)
    return al, ar


prs.configure_font()
data = [load(r[1]) for r in RUNS]

fig = plt.figure(figsize=(6.0, 6.6), layout="constrained")
s_a, s_b, s_z, s_c = fig.subfigures(4, 1, height_ratios=[1.85, 1.15, 1.25, 1.45])

# ---------------- (a) lane + u on real time ----------------
ga = s_a.add_gridspec(2, 1, height_ratios=[0.95, 1.6], hspace=0.05)
gl = ga[0].subgridspec(1, 2, width_ratios=WR, wspace=0.03)
LANE = [s_a.add_subplot(gl[j]) for j in range(2)]
A = broken_pair(s_a, ga[1])
for j in range(2):
    ax = LANE[j]
    lo, hi = SEG[j]
    for r in range(5):
        n = RUNS[r][2]
        y = 4 - r
        ax.add_patch(Rectangle((0, y - 0.32), data[r]["t_end"], 0.64, fc="0.9", ec="none", zorder=1))
        on = data[r]["on"]
        on = on[(on >= lo) & (on <= hi)]
        ax.vlines(on, y - 0.36, y + 0.36, color={1: ORANGE, 8: BLUE}.get(n, "0.15"), lw=1.0 if n != 1 else 0.4, zorder=2)
        if j == 1 and data[r]["on"].size == 0:  # only runs with no slowdowns at all, once, after the break
            ax.text(0.5 * (lo + min(hi, data[r]["t_end"])), y, "no slowdowns", ha="center", va="center", color="0.45", zorder=3)
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
LANE[0].set_yticklabels([("serial" if RUNS[r][2] == 1 else str(RUNS[r][2])) for r in range(5)][::-1])
LANE[0].set_ylabel("Ranks")
tag(LANE[0], "(a)")


def draw(ax_pair, key, lw):
    for ax in ax_pair:
        for k in ORDER:
            n = RUNS[k][2]
            d = data[k]
            ax.plot(d[key], d["u"], color=RUNS[k][3], lw=lw)
            if n == 1:
                ax.plot(d[key][-1], d["u"][-1], marker="x", color=ORANGE, ms=5, mew=1.0, zorder=5)
        ax.set_ylim(*YLIM)


draw(A, "t", 0.8)
A[0].set_ylabel(r"$u$ (mm)")
A[1].set_xlabel(r"Laboratory time, $\sqrt{\lambda_L}\,t$ (s)")
label(A[1], "serial", (271.5, -150), (290, -150), ORANGE)
label(A[1], "8 ranks", (330, 85), (340, 165), BLUE)
label(A[1], "4, 16, 24 ranks", (395, 62), (330, -120), GRAY)

# ---------------- (b) u on integrator time ----------------
B = broken_pair(s_b, s_b.add_gridspec(1, 1)[0])
draw(B, "ti", 0.8)
B[0].set_ylabel(r"$u$ (mm)")
B[1].set_xlabel(r"Integration time, $t_{\mathrm{int}}$ (s)")
tag(B[0], "(b)")
for (z0, z1), nm in zip(ZOOMS, ("(b.1)",)):
    B[1].add_patch(Rectangle((z0, YLIM[0] * 0.95), z1 - z0, (YLIM[1] - YLIM[0]) * 0.95, fill=False, ec="0.2", lw=0.6))
    B[1].text((z0 + z1) / 2, YLIM[1] * 0.98, rf"\textbf{{{nm}}}", ha="center", va="bottom")
label(B[1], "serial: test stopped", (data[0]["ti_end"], data[0]["u"][-1]), (300, -120), ORANGE)

# ---------------- (b.1) zoom on integration time ----------------
z0, z1 = ZOOMS[0]
ax = s_z.add_subplot(s_z.add_gridspec(1, 1)[0])
for k in ORDER:
    n, c = RUNS[k][2], RUNS[k][3]
    d = data[k]
    m = (d["ti"] >= z0 - 1) & (d["ti"] <= z1 + 1)
    ax.plot(d["ti"][m], d["u"][m], color=c, lw=1.1 if n in (1, 8) else 0.9)
    if n == 1 and d["ti"][-1] <= z1:
        ax.plot(d["ti"][-1], d["u"][-1], marker="x", color=ORANGE, ms=5, mew=1.0, zorder=5)
ax.set_xlim(z0, z1)
ax.set_ylim(*YLIM)
ax.set_xlabel(r"$t_{\mathrm{int}}$ (s)")
ax.set_ylabel(r"$u$ (mm)")
tag(ax, "(b.1)")
KR = {RUNS[k][2]: k for k in range(5)}
P8 = (236.6, float(np.interp(236.6, data[KR[8]]["ti"], data[KR[8]]["u"])))
PG = (238.3, float(np.interp(238.3, data[KR[24]]["ti"], data[KR[24]]["u"])))
LABELS = True
if LABELS:
    label(ax, "serial", (242.1, 106), (240.0, 165), ORANGE)
    label(ax, "8 ranks", P8, (P8[0] - 1.2, 150), BLUE)
    label(ax, "4, 16, 24 ranks", PG, (PG[0] - 1.0, -165), GRAY)

# ---------------- (c) delay behind real time (symlog, as fig_fri_lanes.py) ----------------
C = broken_pair(s_c, s_c.add_gridspec(1, 1)[0])
for ax in C:
    for k in ORDER:
        d = data[k]
        ax.plot(d["tk"], d["lag"], color=RUNS[k][3], lw=1.1)
    ax.set_yscale("symlog", linthresh=0.01, linscale=0.4)
    ax.set_ylim(-0.014, 60)  # room below zero for the label of the runs with no delay
C[0].set_yticks([0, 0.01, 0.1, 1, 10])
C[0].set_yticklabels(["0", "0.01", "0.1", "1", "10"])
C[0].set_ylabel("Delay behind\nreal time (s)")
C[1].set_xlabel(r"Laboratory time, $\sqrt{\lambda_L}\,t$ (s)")
tag(C[0], "(c)")
x_end = SEG[1][1]
d0 = data[0]
C[1].plot(d0["tk"][-1], d0["lag"][-1], marker="x", color=ORANGE, ms=5, mew=1.0, zorder=5)
C[1].text(d0["tk"][-1] + 4, d0["lag"][-1], f"serial: {fmt_delay(d0['lag'][-1])} s", ha="left", va="center", color=ORANGE)
end = {RUNS[k][2]: float(np.interp(x_end, data[k]["tk"], data[k]["lag"])) for k in (1, 2, 3, 4)}
print("delay at", x_end, "s:", {n: round(v, 5) for n, v in end.items()}, "; serial at its stop", round(d0["lag"][-1], 2))
# labels at the right edge, one per distinct level (runs within a factor of 1.6 share a label)
grp = []
for n in sorted(end, key=lambda n: end[n]):
    v = end[n]
    if grp and (v < 5e-5 and grp[-1][1] < 5e-5 or 0 < grp[-1][1] and v / grp[-1][1] < 1.6):
        grp[-1][0].append(n)
    else:
        grp.append([[n], v])
for ns, v in grp:
    vals = ", ".join(fmt_delay(end[n]) for n in ns)
    vals = "0" if all(end[n] < 5e-5 for n in ns) else vals
    c = "0.25" if 8 in ns else GRAY
    zero = all(end[n] < 5e-5 for n in ns)
    C[1].text(x_end - 3, -0.0015 if zero else v, f"{', '.join(str(n) for n in ns)} ranks: {vals} s", ha="right",
              va="top" if zero else "bottom", color=c)

# "Slowdowns" centred over both lane segments, on the baseline of the (a) tag
fig.canvas.draw()
fig.set_layout_engine("none")
p0, p1 = LANE[0].get_position(), LANE[1].get_position()
ttl = LANE[0].title
_, yb = fig.transFigure.inverted().transform(ttl.get_transform().transform(ttl.get_position()))
fig.text(0.5 * (p0.x0 + p1.x1), yb, "Slowdowns", ha="center", va="baseline")
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
