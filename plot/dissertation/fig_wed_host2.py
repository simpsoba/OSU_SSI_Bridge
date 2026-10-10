"""Host timing in the rate transition, Wednesday runs (dissertation Ch. 6; read-only on data).

stateOS at 2048 Hz, laboratory clock, model scale (dt_sim = 10/2048 s = 4.88 ms; one count tick = 0.49 ms).
Records end where the 1-rank actuator stops tracking. Start-up = first three target intervals (excluded in b-d).
(a) time to the first target: duration of the first integration step (MKR-alpha operator assembly on every rank,
    first ProfileSPD factorization of alpha on rank 0, first state determination); W02/W04 carry one extra 2-ms
    interval before it, skipped
(b) share of samples in each typeConv3 state (initialize excluded, as plot/PlotStateOSBars.py); numbers above the
    bars = slowdown episodes (entries into state 2) after start-up, k = thousands
(c) host latency per dt_sim window: window start (count resets to 1) to arrival of the new target (s1 rising edge)
    = extrapolation plus slowdown time in that window. Host critical path (force read, assembly, solves, predictor,
    update up to the experimental element) plus any overrun of the work after the target is sent; work that ends
    before the next force is held is not observed, so this is not a speedup of the whole step.
(d) interval between consecutive targets. interval_k = dt_sim + (latency_k - latency_k-1) + slowdown stretch
    (exact to one tick in every parallel step), so the mean is locked to dt_sim while the host keeps up.
(c), (d): box = quartiles, whiskers = 1st and 99th percentiles, x = maximum (triangle = off scale, value printed).
Bold panel tags, no titles or legends (direct labels with leaders), uniform font, constrained layout.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import ctrl_error_wed as C  # noqa: E402
from wed_timing_stats import timing  # noqa: E402
from host_time import host_times  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]  # OSU_SSI_Bridge
sys.path.insert(0, str(ROOT / "plot"))
import PlotResponseSpectrum as prs  # noqa: E402

OUT = Path(sys.argv[1])
ORANGE, BLACK, GRAY = "#C0390B", "#000000", "#8C8C8C"
COL = {1: ORANGE, 4: GRAY, 8: BLACK, 16: GRAY, 24: GRAY}
STATES = [(1, "extrapolate", "#1565c0"), (0, "interpolate", "#2e7d32"), (2, "slowdown", "#FF8C00")]
SLOW_TXT = "#B85C00"
DT = C.DT_SIM * 1e3


def tag(ax, s):
    ax.set_title(rf"\textbf{{{s}}}", loc="left", pad=14)


def boxes(ax, data, x, ymax):
    bp = ax.boxplot(data, positions=x, widths=0.5, whis=(1, 99), showfliers=False, patch_artist=True,
                    medianprops=dict(lw=2.4))
    for k, r in enumerate(R):
        c = COL[r["nr"]]
        bp["boxes"][k].set(facecolor=c, edgecolor=c)
        bp["medians"][k].set(color=c, solid_capstyle="butt")
        for art in (bp["whiskers"][2 * k], bp["whiskers"][2 * k + 1], bp["caps"][2 * k], bp["caps"][2 * k + 1]):
            art.set(color=c, lw=0.9)
        mx = data[k].max()
        ax.plot(x[k], min(mx, ymax), marker="x" if mx <= ymax else "^", color=c, ms=5, mew=1.0)
        if mx > ymax:
            ax.text(x[k] + 0.18, ymax, f"{mx:.0f}", ha="left", va="center", color=c)


def ref_lines(ax, items, x_text, ha="right"):
    for yv, lab in items:
        ax.axhline(yv, color="0.35", lw=0.6, ls=(0, (4, 2)), zorder=0)
        ax.text(x_text, yv + 0.06, lab, ha=ha, va="bottom", color="0.35")


prs.configure_font()
R = []
for wid, stem, dump, nr in C.RUNS:
    d = timing(stem)
    t0 = d["tk"][(1 if d["tk"][1] - d["tk"][0] < 3e-3 else 0) + 3]  # end of start-up
    typ, t = d["typ"], d["t"]
    use = typ != -1
    frac = [100 * np.mean(typ[use] == c) for c, *_ in STATES]
    n_slow = int(np.sum((typ[1:] == 2) & (typ[:-1] != 2) & (t[1:] > t0)))
    ta, h, _ = host_times(stem)
    m = ta > t0
    R.append(dict(wid=wid, nr=nr, first=d["first"] * 1e3, frac=frac, n_slow=n_slow, lat=h[m] * 1e3,
                  dt=np.diff(ta[m]) * 1e3))
for r in R:
    print(f"{r['wid']} {r['nr']:2d}: first {r['first']:.1f} ms  slowdowns {r['n_slow']}  latency mean {r['lat'].mean():.2f}"
          f"  interval mean {r['dt'].mean():.3f} (x{r['dt'].mean() / DT:.3f}) max {r['dt'].max():.1f}")
x = np.arange(len(R))
m1 = R[0]["lat"].mean()

fig, ((ax_a, ax_b), (ax_c, ax_d)) = plt.subplots(2, 2, figsize=(6.0, 4.6), layout="constrained")

# (a) time to first target
first = np.array([r["first"] for r in R])
ax_a.bar(x, first, width=0.6, color=[COL[r["nr"]] for r in R], edgecolor="none")
for xi, v in zip(x, first):
    ax_a.text(xi, v + 5, f"{v:.0f}", ha="center", va="bottom")
ax_a.set_ylim(0, 285)
ax_a.set_ylabel("Time to first target (ms)")

# (b) state shares; names on the right with leaders to the 24-rank bar; slowdown episodes on top
bottom = np.zeros(len(R))
seg = {}
for j, (c, name, col) in enumerate(STATES):
    h = np.array([r["frac"][j] for r in R])
    ax_b.bar(x, h, bottom=bottom, width=0.6, color=col, edgecolor="w", lw=0.4)
    for i in range(len(R)):
        if j < 2 and h[i] >= 12:
            ax_b.text(x[i], bottom[i] + 0.5 * h[i], f"{h[i]:.0f}", ha="center", va="center", color="w")
    seg[c] = (bottom.copy(), h.copy())
    bottom += h
for xi, r in zip(x, R):
    n = r["n_slow"]
    ax_b.text(xi, 102, f"{n / 1000:.0f}k" if n >= 1000 else f"{n}", ha="center", va="bottom", color=SLOW_TXT)
xr = x[-1] + 0.3  # right edge of the last bar
for c, name, col in STATES[:2]:
    b0, h = seg[c]
    yb = b0[-1] + 0.5 * h[-1]
    ax_b.annotate(name, (xr, yb), xytext=(xr + 0.55, yb), ha="left", va="center", color=col,
                  arrowprops=dict(arrowstyle="-", color=col, lw=0.6, shrinkA=1, shrinkB=0))
ax_b.text(xr + 0.25, 102, "slowdown\nepisodes", ha="left", va="bottom", color=SLOW_TXT, linespacing=1.0)
ax_b.set_ylim(0, 100)
ax_b.set_yticks([0, 25, 50, 75, 100])
ax_b.set_ylabel(r"Share of record (\%)")
ax_b.set_xlim(-0.6, len(R) + 1.05)
ax_b.spines["right"].set_visible(False)

# (c) host latency per window
boxes(ax_c, [r["lat"] for r in R], x, 8.3)
ref_lines(ax_c, ((0.8 * DT, r"$0.8\,\Delta t_{\mathrm{sim}}$"), (DT, r"$\Delta t_{\mathrm{sim}}$")), len(R) - 0.45)
ax_c.set_ylim(0, 8.7)
ax_c.set_ylabel("Host latency per window (ms)")

# (d) interval between targets
boxes(ax_d, [r["dt"] for r in R], x, 13.0)
ax_d.axhline(DT, color="0.35", lw=0.6, ls=(0, (4, 2)), zorder=0)
ax_d.annotate(r"$\Delta t_{\mathrm{sim}}$", (3.5, DT), xytext=(3.5, 10.2), ha="center", va="bottom", color="0.35",
              arrowprops=dict(arrowstyle="-", color="0.35", lw=0.6, shrinkA=1, shrinkB=0))
ax_d.set_ylim(2.5, 13.7)
ax_d.set_ylabel("Interval between targets (ms)")

for ax, tg in ((ax_a, "(a)"), (ax_b, "(b)"), (ax_c, "(c)"), (ax_d, "(d)")):
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r['nr']}\n{r['wid']}" for r in R])
    if ax is not ax_b:
        ax.set_xlim(-0.6, len(R) - 0.4)
    ax.set_xlabel("Ranks")
    tag(ax, tg)
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
