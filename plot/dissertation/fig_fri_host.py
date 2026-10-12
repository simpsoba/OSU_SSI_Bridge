"""Host timing of the Friday main-campaign tests (dissertation Ch. 6; read-only on data).

(a) host latency per dt_sim window (window start to new-target arrival), start-up windows excluded; boxes = quartiles,
    whiskers = 1st and 99th percentiles; dashed lines at 0.8 dt_sim (slowdown threshold) and dt_sim
(b) share of samples in each typeConv3 state (initialize excluded); numbers above the bars = slowdown episodes after
    start-up (k = thousands); runs stopped early are marked with a dagger
Tests grouped by model (baseline, moderate, large, extra-large) and by soil formulation; colors by configuration:
CPU MKR-alpha + ProfileSPD, CPU MKR-alpha + MUMPS, CPU MKR-alpha + cuDSS solves, GPU MKR-alpha + cuDSS.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import fri_runs as R  # noqa: E402

sys.path.insert(0, str(R.ROOT / "plot"))
import PlotResponseSpectrum as prs  # noqa: E402

OUT = Path(sys.argv[1])
CFG = {("CPU", "ProfileSPD"): ("#C0390B", "CPU MKR-$\\alpha$ + ProfileSPD"),
       ("CPU", "MUMPS"): ("#E69F00", "CPU MKR-$\\alpha$ + MUMPS"),
       ("CPU", "cuDSS"): ("#0072B2", "CPU MKR-$\\alpha$ + cuDSS solves"),
       ("GPU", "cuDSS"): ("#000000", "GPU MKR-$\\alpha$ + cuDSS")}
GROUPS = [("Baseline\n2028 DOFs", ["F03", "F04", "F26", "F22", "F05", "F27"]),
          ("Moderate\n3216", ["F07", "F08", "F06"]),
          ("Large\n4620", ["F14", "F15", "F09", "F10", "F12", "F11"]),
          ("XL\n5484", ["F13"]),
          ("Soil formulation\n2028", ["F24", "F23", "F20"])]
FORM = {"F24": "elastic", "F23": "stiff", "F20": "Quad\nPIMY"}
STATES = [(1, "extrapolate", "#1565c0"), (0, "interpolate", "#2e7d32"), (2, "slowdown", "#FF8C00")]
SLOW_TXT = "#B85C00"
DT = R.DT_SIM * 1e3


def cfg(r):
    return CFG[(r[5], r[6].split()[0])]


prs.configure_font()
ids, xs, gx = [], [], []
x = 0.0
for name, tl in GROUPS:
    gx.append((name, x, x + len(tl) - 1))
    for tid in tl:
        ids.append(tid)
        xs.append(x)
        x += 1
    x += 0.8
xs = np.array(xs)
D = {tid: R.load(tid) for tid in ids}

fig, (ax_a, ax_b) = plt.subplots(2, 1, figsize=(6.0, 4.6), layout="constrained", sharex=True,
                                 gridspec_kw=dict(height_ratios=[1.15, 1.0]))
lat = []
for tid in ids:
    v = D[tid]["lat"][3:] * 1e3
    lat.append(v[np.isfinite(v)])
bp = ax_a.boxplot(lat, positions=xs, widths=0.6, whis=(1, 99), showfliers=False, patch_artist=True,
                  medianprops=dict(lw=1.6))
for k, tid in enumerate(ids):
    c = cfg(R.BY_ID[tid])[0]
    bp["boxes"][k].set(facecolor=c, edgecolor=c)
    bp["medians"][k].set(color="w" if c == "#000000" else c, solid_capstyle="butt")
    for art in (bp["whiskers"][2 * k], bp["whiskers"][2 * k + 1], bp["caps"][2 * k], bp["caps"][2 * k + 1]):
        art.set(color=c, lw=0.9)
    ax_a.plot(xs[k], lat[k].mean(), marker="o", ms=2.6, mfc="white", mec=c, mew=0.7, zorder=5)
for yv, txt in ((0.8 * DT, r"$0.8\,\Delta t_{\mathrm{sim}}$"), (DT, r"$\Delta t_{\mathrm{sim}}$")):
    ax_a.axhline(yv, color="0.35", lw=0.6, ls=(0, (4, 2)), zorder=0)
    ax_a.text(xs[-1] + 0.6, yv, txt, ha="left", va="center", color="0.35")
ax_a.set_ylim(0, 6.2)
ax_a.set_ylabel("Host latency (ms)")
for name, x0, x1 in gx:
    ax_a.text(0.5 * (x0 + x1), 6.35, name, ha="center", va="bottom", linespacing=0.95)
    for ax in (ax_a, ax_b):
        if x0 > 0:
            ax.axvline(x0 - 0.9, color="0.75", lw=0.5, zorder=0)
handles = [plt.Rectangle((0, 0), 1, 1, fc=c, ec=c) for c, _ in CFG.values()]
fig.legend(handles, [l for _, l in CFG.values()], loc="outside upper center", ncol=2, frameon=False,
           handlelength=1.0, columnspacing=1.2)

bottom = np.zeros(len(ids))
for c, name, col in STATES:
    h = np.array([D[tid]["frac"][c] for tid in ids])
    ax_b.bar(xs, h, bottom=bottom, width=0.74, color=col, edgecolor="w", lw=0.3)
    for xi, b0, hj in zip(xs, bottom, h):
        if c != 2 and hj >= 12:
            ax_b.text(xi, b0 + 0.5 * hj, f"{hj:.0f}", ha="center", va="center", color="w", fontsize=6.5)
    bottom += h
for xi, tid in zip(xs, ids):
    d = D[tid]
    n = d["n_slow"]
    s = f"{n / 1000:.1f}k" if n >= 1000 else f"{n}"
    if not d["complete"]:
        s += r"$^\dagger$"
    ax_b.text(xi, 101.5, s, ha="center", va="bottom", color=SLOW_TXT, fontsize=6.5, rotation=90 if n >= 100 else 0)
ax_b.set_ylim(0, 100)
ax_b.set_ylabel(r"Share of record (\%)")
ax_b.set_xticks(xs)
ranks = {tid: ("serial" if R.BY_ID[tid][7] == 1 else R.BY_ID[tid][7]) for tid in ids}
ax_b.set_xticklabels([f"{tid}\n{FORM.get(tid, ranks[tid])}" for tid in ids], fontsize=6.5, linespacing=0.95)
ax_b.set_xlim(xs[0] - 0.6, xs[-1] + 2.6)
xr = xs[-1] + 0.45
bot = 0.0
for c, name, col in STATES[:2]:
    hl = D[ids[-1]]["frac"][c]
    ax_b.annotate(name, (xr, bot + 0.5 * hl), xytext=(xr + 0.35, bot + 0.5 * hl), ha="left", va="center",
                  color=col, fontsize=7, arrowprops=dict(arrowstyle="-", color=col, lw=0.6, shrinkA=1, shrinkB=0))
    bot += hl
ax_b.text(xr + 0.2, 101.5, "slowdown\nepisodes", ha="left", va="bottom", color=SLOW_TXT, fontsize=7,
          linespacing=0.95)
fig.text(0.01, 0.005, r"Tick labels: test ID and MPI ranks (serial: sequential interpreter). $^\dagger$Run stopped early.", fontsize=6.5,
         color="0.35", ha="left", va="bottom")
for ax, tg in ((ax_a, "(a)"), (ax_b, "(b)")):
    ax.set_title(rf"\textbf{{{tg}}}", loc="left", pad=22 if tg == "(a)" else 14)
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
