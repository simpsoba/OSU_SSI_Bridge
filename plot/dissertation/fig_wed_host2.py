"""Host timing in the rate transition, Wednesday runs (dissertation Ch. 6; read-only on data).

stateOS at 2048 Hz, laboratory clock, model scale (dt_sim = 10/2048 s = 4.88 ms; one count tick = 0.49 ms).
Records end where the 1-rank actuator stops tracking. The controller logs one new-target flag per integration step
plus two start-up handshake flags about 2 ms apart (in W05-W07 the first precedes the record); target intervals are
counted from the second handshake flag. All panels include the whole record.
(a) interval between consecutive targets; markers = the three start-up intervals: 1st (triangle) handshake to the
    first integration target (analysis setup, assembly of M, A-tilde, and A on every rank, merge on rank 0,
    factorization of A-tilde at the predictor's first solve, first predictor and update pass); 2nd (square) first
    force read and factorizations of M and A at the corrector's first solves; 3rd (diamond) first-pass overhead
(b) host latency per dt_sim window: window start (count resets to 1) to arrival of the new target (s1 rising edge)
    = extrapolation plus slowdown time in that window; markers = the first three windows. Host critical path plus
    any overrun of the work after the target is sent; not the total computational time of a step.
(a), (b): one log axis; box = quartiles, whiskers = 1st and 99th percentiles; dashed line dt_sim in (a),
    0.8 dt_sim (slowdown threshold) in (b)
(c) share of samples in each typeConv3 state (initialize excluded, as plot/PlotStateOSBars.py); numbers above the
    bars = slowdown episodes after start-up, with those during start-up below (+n; k = thousands)
One row; bold panel tags, no titles or legends (direct labels), uniform font, constrained layout.
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
MARK = (("^", "1st"), ("s", "2nd"), ("D", "3rd"))  # start-up intervals / windows

prs.configure_font()
R = []
for wid, stem, dump, nr in C.RUNS:
    d = timing(stem)
    tk, typ, t = d["tk"], d["typ"], d["t"]
    k0 = 1 if tk[1] - tk[0] < 3e-3 else 0  # second handshake flag
    use = typ != -1
    on = t[1:][(typ[1:] == 2) & (typ[:-1] != 2)]
    t0 = tk[k0 + 3]  # end of the three start-up intervals
    ta, h, _ = host_times(stem)
    R.append(dict(wid=wid, nr=nr, iv=np.diff(tk[k0:]) * 1e3, lat=h[ta >= tk[k0 + 1]] * 1e3,
                  frac=[100 * np.mean(typ[use] == c) for c, *_ in STATES],
                  n=int(np.sum(on > t0)), n_su=int(np.sum((on >= tk[0]) & (on <= t0)))))
for r in R:
    print(f"{r['wid']} {r['nr']:2d}: start-up intervals {np.round(r['iv'][:3], 1)} ms, windows {np.round(r['lat'][:3], 1)} ms;"
          f" mean interval {r['iv'].mean():.3f} ms, mean latency {r['lat'].mean():.2f} ms; slowdowns {r['n']} (+{r['n_su']})")
x = np.arange(len(R))


def panel(ax, key):
    bp = ax.boxplot([r[key] for r in R], positions=x, widths=0.5, whis=(1, 99), showfliers=False,
                    patch_artist=True, medianprops=dict(lw=2.0))
    for k, r in enumerate(R):
        c = COL[r["nr"]]
        bp["boxes"][k].set(facecolor=c, edgecolor=c)
        bp["medians"][k].set(color=c, solid_capstyle="butt")
        for art in (bp["whiskers"][2 * k], bp["whiskers"][2 * k + 1], bp["caps"][2 * k], bp["caps"][2 * k + 1]):
            art.set(color=c, lw=0.9)
        for j, v in enumerate(r[key][:3]):
            ax.plot(x[k], v, marker=MARK[j][0], ms=3.8 if j < 2 else 3.3, mfc="white", mec=c, mew=0.8, zorder=5)
            if key == "iv" and k == len(R) - 1:  # key once, beside the 24-rank column
                ax.text(x[k] + 0.2, v, MARK[j][1], ha="left", va="center", color="0.35")
    ax.set_yscale("log")
    ax.set_ylim(1.0, 400)
    ax.set_xticks(x)
    ax.set_xlim(-0.6, len(R) - 0.4)


fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(6.0, 2.6), layout="constrained",
                                       gridspec_kw=dict(width_ratios=[1.0, 1.0, 1.45]))
# (a) intervals, (b) latency: one log axis
panel(ax_a, "iv")
ax_a.axhline(DT, color="0.35", lw=0.6, ls=(0, (4, 2)), zorder=0)
ax_a.text(-0.5, DT * 0.86, r"$\Delta t_{\mathrm{sim}}$", ha="left", va="top", color="0.35")
ax_a.set_ylabel("Time (ms)")
panel(ax_b, "lat")
ax_b.sharey(ax_a)
ax_b.tick_params(labelleft=False)
ax_b.axhline(0.8 * DT, color="0.35", lw=0.6, ls=(0, (4, 2)), zorder=0)
ax_b.text(2.5, 0.8 * DT * 1.07, r"$0.8\,\Delta t_{\mathrm{sim}}$", ha="center", va="bottom", color="0.35")

# (c) states; slowdown episodes after start-up (start-up in parentheses)
bottom = np.zeros(len(R))
for j, (c, name, col) in enumerate(STATES):
    hh = np.array([r["frac"][j] for r in R])
    ax_c.bar(x, hh, bottom=bottom, width=0.74, color=col, edgecolor="w", lw=0.4)
    for xi, b0, hj in zip(x, bottom, hh):  # percentages inside the extrapolate / interpolate segments
        if c != 2 and hj >= 10:
            ax_c.text(xi, b0 + 0.5 * hj, f"{hj:.0f}", ha="center", va="center", color="w")
    bottom += hh
for xi, r in zip(x, R):
    n = f"{r['n'] / 1000:.0f}k" if r["n"] >= 1000 else f"{r['n']}"
    ax_c.text(xi, 101.5, f"{n}\n+{r['n_su']}", ha="center", va="bottom", color=SLOW_TXT, linespacing=0.95)
ax_c.set_ylim(0, 100)
ax_c.set_xticks(x)
ax_c.set_xlim(-0.6, len(R) + 1.05)
xr = x[-1] + 0.37
bot = 0.0
for j, (c, name, col) in enumerate(STATES[:2]):
    h_last = R[-1]["frac"][j]
    yb = bot + 0.5 * h_last
    ax_c.annotate(name, (xr, yb), xytext=(xr + 0.4, yb), ha="left", va="center", color=col,
                  arrowprops=dict(arrowstyle="-", color=col, lw=0.6, shrinkA=1, shrinkB=0))
    bot += h_last
ax_c.text(xr + 0.3, 101.5, "slowdown\nepisodes", ha="left", va="bottom", color=SLOW_TXT, linespacing=0.95)
ax_c.spines["right"].set_visible(False)
ax_c.set_ylabel(r"Share of record (\%)")

for ax, tg in ((ax_a, "(a)"), (ax_b, "(b)"), (ax_c, "(c)")):
    ax.set_xticklabels([str(r["nr"]) for r in R])
    ax.set_xlabel("Ranks")
    ax.set_title(rf"\textbf{{{tg}}}", loc="left", pad=24)
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
