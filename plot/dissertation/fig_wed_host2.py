"""Host timing in the rate transition, Wednesday runs (dissertation Ch. 6; read-only on data).

stateOS at 2048 Hz, laboratory clock, model scale (dt_sim = 10/2048 s = 4.88 ms; one count tick = 0.49 ms).
Records end where the 1-rank actuator stops tracking. The controller logs one new-target flag per integration step
plus one (W05-W07) or two (W02, W04) start-up handshake flags within 4 ms of the start of the record (initialize
state, zero target); intervals are counted from the last handshake flag. Panels (b)-(d) include the whole record.
(a) the three start-up intervals of each test, log scale: (1) handshake to the first integration target: transient-
    analysis setup, assembly of M, A-tilde, and A on every rank, merge on rank 0, factorization of A-tilde at the
    predictor's first solve, first predictor and update pass; (2) first to second target: first force read,
    factorization of M and A at the corrector's first solves, second predictor; (3) first-pass warm-up
(d) share of samples in each typeConv3 state (initialize excluded, as plot/PlotStateOSBars.py); numbers above the
    bars = slowdown episodes (entries into state 2) after start-up + during start-up, k = thousands
(c) host latency per dt_sim window: window start (count resets to 1) to arrival of the new target (s1 rising edge)
    = extrapolation plus slowdown time in that window. Host critical path (force read, assembly, solves, predictor,
    update up to the experimental element) plus any overrun of the work after the target is sent; work that ends
    before the next force is held is not observed, so this is not a speedup of the whole step.
(b) interval between consecutive targets. interval_k = dt_sim + (latency_k - latency_k-1) + slowdown stretch
    (exact to one tick in every parallel step), so the mean is locked to dt_sim while the host keeps up.
(b), (c): box = quartiles, whiskers = 1st and 99th percentiles, x = maximum (triangle = off scale, value printed).
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
    tk = d["tk"]
    k0 = 1 if tk[1] - tk[0] < 3e-3 else 0  # last start-up handshake flag
    t_first, t0 = tk[k0 + 1], tk[k0 + 3]  # first integration target; end of the three start-up intervals
    typ, t = d["typ"], d["t"]
    use = typ != -1
    frac = [100 * np.mean(typ[use] == c) for c, *_ in STATES]
    on = t[1:][(typ[1:] == 2) & (typ[:-1] != 2)]
    n_slow, n_su = int(np.sum(on > t0)), int(np.sum((on >= tk[0]) & (on <= t0)))
    ta, h, _ = host_times(stem)
    iv = np.diff(tk[k0:]) * 1e3  # whole record, from the last handshake flag
    R.append(dict(wid=wid, nr=nr, first=d["first"] * 1e3, su=iv[:3], frac=frac, n_slow=n_slow, n_su=n_su,
                  lat=h[ta >= t_first] * 1e3, dt=iv))
for r in R:
    print(f"{r['wid']} {r['nr']:2d}: start-up {np.round(r['su'], 1)} ms  slowdowns {r['n_slow']}+{r['n_su']}  latency mean {r['lat'].mean():.2f}"
          f"  interval mean {r['dt'].mean():.3f} (x{r['dt'].mean() / DT:.3f}) max {r['dt'].max():.1f}")
x = np.arange(len(R))
m1 = R[0]["lat"].mean()

# Story order: (a) start-up intervals, (b) intervals between targets, (c) host latency, (d) rate-transition states.
# ax_d draws the intervals (top right, tag b); ax_b draws the states (bottom right, tag d).
fig, ((ax_a, ax_d), (ax_c, ax_b)) = plt.subplots(2, 2, figsize=(6.0, 4.6), layout="constrained")

# (a) the three start-up intervals of each test (log scale)
BW, SHADE = 0.26, (1.0, 0.6, 0.32)
for i, r in enumerate(R):
    for j in range(3):
        xj = x[i] + (j - 1) * BW
        ax_a.bar(xj, r["su"][j], width=BW * 0.9, color=COL[r["nr"]], alpha=SHADE[j], edgecolor="none")
        if i == 0:
            ax_a.text(xj, 1.15, f"{j + 1}", ha="center", va="bottom", color="w" if j == 0 else "0.2")
    ax_a.text(x[i] - BW, r["su"][0] * 1.08, f"{r['su'][0]:.0f}", ha="center", va="bottom")
ax_a.axhline(DT, color="0.35", lw=0.6, ls=(0, (4, 2)), zorder=0)
ax_a.text(len(R) - 0.45, DT * 0.92, r"$\Delta t_{\mathrm{sim}}$", ha="right", va="top", color="0.35")
ax_a.set_yscale("log")
ax_a.set_ylim(1.0, 600)
ax_a.set_ylabel("Start-up target interval (ms)")

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
    n, n_su = r["n_slow"], r["n_su"]
    lab = f"{n / 1000:.0f}k" if n >= 1000 else f"{n}+{n_su}"
    ax_b.text(xi, 102, lab, ha="center", va="bottom", color=SLOW_TXT)
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

for ax, tg in ((ax_a, "(a)"), (ax_d, "(b)"), (ax_c, "(c)"), (ax_b, "(d)")):
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r['nr']}\n{r['wid']}" for r in R])
    if ax is not ax_b:
        ax.set_xlim(-0.6, len(R) - 0.4)
    ax.set_xlabel("Ranks")
    tag(ax, tg)
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
