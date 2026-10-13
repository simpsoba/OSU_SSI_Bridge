"""Rate-transition state shares of the main-campaign tests (dissertation Ch. 6; read-only on data).

One stacked bar per test: share of the stateOS samples (initialize excluded) in extrapolate, interpolate, and slowdown,
as in fig_wed_host2.py (c). Numbers above the bars: slowdown episodes after start-up, with those during start-up below
(+n; k = thousands). Tests grouped by model and ordered by configuration; tick labels = test ID and MPI ranks
(serial = sequential interpreter); configuration marker below each bar (open = not completed, NC).
Uniform 9 pt font; bold panel tag; constrained layout.
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
CFG = {"pspd": ("#C0390B", "s", r"CPU MKR-$\alpha$, ProfileSPD"), "mumps": ("#8C510A", "D", r"CPU MKR-$\alpha$, MUMPS"),
       "cpucudss": ("#2166AC", "^", r"CPU MKR-$\alpha$, cuDSS solves"), "gpu": ("#000000", "o", r"GPU MKR-$\alpha$, cuDSS")}
GROUPS = [("Baseline, 2,028 DOFs", [("F03", "pspd"), ("F04", "pspd"), ("F26", "mumps"), ("F22", "gpu"), ("F05", "gpu"),
                                    ("F27", "gpu")]),
          ("Moderate, 3,216", [("F07", "pspd"), ("F08", "mumps"), ("F06", "gpu")]),
          ("Large, 4,620", [("F14", "cpucudss"), ("F15", "cpucudss"), ("F09", "gpu"), ("F10", "gpu"), ("F12", "gpu"),
                            ("F11", "gpu")]),
          ("XL, 5,484", [("F13", "gpu")])]
STATES = [(1, "extrapolate", "#1565c0"), (0, "interpolate", "#2e7d32"), (2, "slowdown", "#FF8C00")]
SLOW_TXT = "#B85C00"

prs.configure_font()
fig, ax = plt.subplots(figsize=(6.0, 3.1), layout="constrained")
x, xs, info, gx = 0.0, [], [], []
for name, tl in GROUPS:
    x0 = x
    for tid, cfg in tl:
        xs.append(x)
        info.append((tid, cfg))
        x += 1
    gx.append((name, x0, x - 1))
    x += 0.7
xs = np.array(xs)
bottom = np.zeros(xs.size)
D = {tid: R.load(tid) for tid, _ in info}
for c_, name, col in STATES:
    h = np.array([D[tid]["frac"][c_] for tid, _ in info])
    ax.bar(xs, h, bottom=bottom, width=0.74, color=col, edgecolor="w", lw=0.4)
    for xi, b0, hj in zip(xs, bottom, h):
        if c_ != 2 and hj >= 14:
            ax.text(xi, b0 + 0.5 * hj, f"{hj:.0f}", ha="center", va="center", color="w", rotation=90)
    bottom += h
for xi, (tid, cfg) in zip(xs, info):
    d = D[tid]
    n = d["n_slow"]
    s = f"{n / 1000:.1f}k" if n >= 1000 else f"{n}"
    ax.text(xi, 101.5, f"{s}\n+{d['n_slow_su']}", ha="center", va="bottom", color=SLOW_TXT, linespacing=0.95)
    c, mk, _ = CFG[cfg]
    ax.plot(xi, -7.5, marker=mk, ms=4.6, color=c, mfc=c if d["complete"] else "white", mew=0.9, ls="none",
            clip_on=False)
    print(f"{tid}: extrap {d['frac'][1]:.1f}%, interp {d['frac'][0]:.1f}%, slow {d['frac'][2]:.2f}%, "
          f"episodes {n} (+{d['n_slow_su']})")
ax.set_ylim(0, 100)
ax.set_xticks(xs)
ax.set_xticklabels([f"{tid}\n" + ("serial" if R.BY_ID[tid][7] == 1 else f"{R.BY_ID[tid][7]}") +
                    ("\nNC" if not D[tid]["complete"] else "\n") for tid, _ in info], linespacing=0.95)
ax.tick_params(axis="x", length=0, pad=16)
ROW_HEADERS = True  # "Test" / "Ranks" left of the first tick label, aligned to it after the layout is fixed
for name, a, b in gx:
    ax.annotate(name, xy=(0.5 * (a + b), 0), xycoords=("data", "axes fraction"), xytext=(0, -50),
                textcoords="offset points", ha="center", va="top", fontweight="bold")
ax.set_xlim(xs[0] - 0.6, xs[-1] + 2.3)
xr = xs[-1] + 0.45
bot = 0.0
for c_, name, col in STATES[:2]:
    hl = D[info[-1][0]]["frac"][c_]
    ax.annotate(name, (xr, bot + 0.5 * hl), xytext=(xr + 0.35, bot + 0.5 * hl), ha="left", va="center", color=col,
                arrowprops=dict(arrowstyle="-", color=col, lw=0.6, shrinkA=1, shrinkB=0))
    bot += hl
ax.text(xr + 0.2, 101.5, "slowdown\nepisodes", ha="left", va="bottom", color=SLOW_TXT, linespacing=0.95)
ax.spines["right"].set_visible(False)
ax.set_ylabel(r"Share of record (\%)")
fig.legend([plt.Line2D([], [], marker=CFG[k][1], color=CFG[k][0], ls="none", ms=4.6) for k in CFG],
           [CFG[k][2] for k in CFG], loc="outside upper center", ncol=2, frameon=False, handletextpad=0.3,
           columnspacing=1.5)
if ROW_HEADERS:
    fig.canvas.draw()
    fig.set_layout_engine("none")
    t0 = ax.get_xticklabels()[0]
    bb = t0.get_window_extent()
    xf, yf = fig.transFigure.inverted().transform((bb.x0 - 6, bb.y1))
    fig.text(xf, yf, "Test\nRanks", ha="right", va="top", color="0.35", linespacing=0.95,
             multialignment="right")
fig.savefig(OUT.with_suffix(".png"), dpi=250)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
