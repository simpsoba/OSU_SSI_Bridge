"""Interface displacement, Wednesday runs (dissertation Ch. 6; read-only on data).

u = meaSigOS x lambda (prototype mm), Simulink at 2048 Hz model. Records end where the 1-rank actuator stops
tracking (0.25-s RMS of target - measured above 2 mm model). Integrator time from the received-target count:
target k at t_int = (k + alpha_f) dt_int. Delay = laboratory time of each target minus k dt_sim (laboratory s).
(a)         slowdown lane (one row per run, test ID in the free space) and u on real time, sqrt(lambda_L) t
(b)         u on integrator time; (b.1), (b.2) zooms with a slowdown lane each
(c)         delay of the received targets behind hard real time; inset for the parallel runs
Bold panel tags, no titles or legends (direct labels), uniform font, constrained layout.
"""
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
ORANGE, BLUE, GRAY = "#C0390B", "#000000", "#8C8C8C"  # 1 rank deep vermillion, 8 ranks black, others gray
RUNS = [("W04", "0819_GusBridge_Serial_H0p5_T7p746_reg_Trial02", 1, ORANGE),
        ("W06", "0819_GusBridge_Parallel4_H0p5_T7p746_reg_Trial02", 4, GRAY),
        ("W02", "0819_GusBridge_Parallel8_H0p5_T7p746_reg", 8, BLUE),
        ("W05", "0819_GusBridge_Parallel16_H0p5_T7p746_reg_Trial02", 16, GRAY),
        ("W07", "0819_GusBridge_Parallel24_H0p5_T7p746_reg", 24, GRAY)]
ORDER = [1, 3, 4, 2, 0]  # grays, 8 ranks, then 1 rank on top
SEG = [(0, 170), (170, 420)]
WR = [1, 3.4]
ZOOMS = [(214, 222), (241, 249)]
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
    lag = (tk - tk[0]) - DT_SIM * np.arange(tk.size)
    on = t[1:][(typ[1:] == 2) & (typ[:-1] != 2)]
    err = z["tarSigOS_data"][:, 0] - me[:, 0]
    rms = np.sqrt(np.convolve(err ** 2, np.ones(512) / 512, mode="same"))
    bad = np.flatnonzero((rms > 2e-3) & (me[:, -1] > 50.0))
    me = me[: (max(bad[0] - 256, 0) if bad.size else me.shape[0])]
    t_end = me[-1, -1]
    to_int = lambda x: (np.interp(x, tk, np.arange(tk.size)) + ALPHA_F) * DT_SIM * SQL
    on = on[(on > 0.7) & (on <= t_end)]
    kk = tk <= t_end
    return dict(t=me[:, -1] * SQL, ti=to_int(me[:, -1]), u=me[:, 0] * LAM * 1e3, tk=tk[kk] * SQL, lag=lag[kk],
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
s_a, s_b, s_z, s_c = fig.subfigures(4, 1, height_ratios=[1.85, 1.15, 1.25, 1.15])

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
label(A[1], "1 rank", (271.5, -150), (290, -150), ORANGE)
label(A[1], "8 ranks", (330, 85), (340, 165), BLUE)
label(A[1], "4, 16, 24 ranks", (395, 62), (330, -120), GRAY)

# ---------------- (b) u on integrator time ----------------
B = broken_pair(s_b, s_b.add_gridspec(1, 1)[0])
draw(B, "ti", 0.8)
B[0].set_ylabel(r"$u$ (mm)")
B[1].set_xlabel(r"Integrator time, $t_{\mathrm{int}}$ (s)")
tag(B[0], "(b)")
for (z0, z1), nm in zip(ZOOMS, ("(b.1)", "(b.2)")):
    B[1].add_patch(Rectangle((z0, YLIM[0] * 0.95), z1 - z0, (YLIM[1] - YLIM[0]) * 0.95, fill=False, ec="0.2", lw=0.6))
    B[1].text((z0 + z1) / 2, YLIM[1] * 0.98, rf"\textbf{{{nm}}}", ha="center", va="bottom")
label(B[1], "1 rank: test stopped", (data[0]["ti_end"], data[0]["u"][-1]), (300, -120), ORANGE)

# ---------------- (b.1), (b.2) zooms on integrator time ----------------
gz = s_z.add_gridspec(2, 2, height_ratios=[0.16, 1.0], hspace=0.04)
Z = [s_z.add_subplot(gz[1, j]) for j in range(2)]
ZL = [s_z.add_subplot(gz[0, j], sharex=Z[j]) for j in range(2)]
for j, (z0, z1) in enumerate(ZOOMS):
    ax, lz = Z[j], ZL[j]
    for k in ORDER:
        n, c = RUNS[k][2], RUNS[k][3]
        d = data[k]
        m = (d["ti"] >= z0 - 1) & (d["ti"] <= z1 + 1)
        ax.plot(d["ti"][m], d["u"][m], color=c, lw=1.1 if n in (1, 8) else 0.9)
        if n == 1 and d["ti"][-1] <= z1:
            ax.plot(d["ti"][-1], d["u"][-1], marker="x", color=ORANGE, ms=5, mew=1.0, zorder=5)
        on = d["on_i"][(d["on_i"] >= z0) & (d["on_i"] <= z1)]
        if n != 1:
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
    ax.set_ylim(*YLIM)
    ax.set_xlabel(r"$t_{\mathrm{int}}$ (s)")
Z[0].set_ylabel(r"$u$ (mm)")
Z[1].tick_params(labelleft=False)
label(Z[0], "1 rank", (219.3, 45), (219.8, 140), ORANGE)
label(Z[0], "8 ranks", (217.0, 38), (215.0, 140), BLUE)
label(Z[0], "24 ranks", (216.8, 38), (215.0, -150), GRAY)
label(Z[1], "1 rank", (242.1, 106), (243.3, 160), ORANGE)
label(Z[1], "4, 8, 16, 24 ranks", (245.4, 80), (245.9, 150), "0.25")

# ---------------- (c) delay behind real time ----------------
C = broken_pair(s_c, s_c.add_gridspec(1, 1)[0])
for ax in C:
    for k in ORDER:
        d = data[k]
        ax.plot(d["tk"], d["lag"], color=RUNS[k][3], lw=1.1)
C[0].set_ylim(-1, 22)
C[0].set_ylabel("Delay,\nlaboratory (s)")
C[1].set_xlabel(r"Laboratory time, $\sqrt{\lambda_L}\,t$ (s)")
tag(C[0], "(c)")
label(C[1], "1 rank", (250, np.interp(250, data[0]["tk"], data[0]["lag"])), (215, 19), ORANGE)
ins = C[1].inset_axes([0.52, 0.40, 0.31, 0.52])
for k in (1, 3, 4, 2):
    d = data[k]
    ins.plot(d["tk"], d["lag"], color=RUNS[k][3], lw=1.0)
lag_end = {RUNS[k][2]: np.interp(SEG[1][1], data[k]["tk"], data[k]["lag"]) for k in (1, 2, 3, 4)}
blend = mtransforms.blended_transform_factory(ins.transAxes, ins.transData)
for lab, y, c in (("24 ranks", lag_end[24], GRAY), ("16 ranks", lag_end[16] + 0.025, GRAY),
                  ("4, 8 ranks", 0.5 * (lag_end[4] + lag_end[8]) - 0.025, "0.25")):
    ins.text(1.03, y, lab, transform=blend, ha="left", va="center", color=c)
ins.set_xlim(0, SEG[1][1])
ins.set_ylim(0, 0.35)
ins.set_xticks([0, 200, 400])
ins.grid(True)
for s in ("top", "right"):
    ins.spines[s].set_visible(True)

fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
