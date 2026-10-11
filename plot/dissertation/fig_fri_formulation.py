"""Effect of the soil formulation on host latency, Friday baseline model, GPU MKR-alpha + cuDSS, 8 ranks
(dissertation Ch. 6; read-only on data).

Tests: F24 SSPQuad + elastic, F23 SSPQuad + PIMY stiff soil, F05 and F27 SSPQuad + PIMY (soft soil, reference),
F20 Quad + PIMY (soft soil).
(a) host latency per dt_sim window averaged over 2 s of integration time (start-up excluded), with the 5-95 % Arias
    band and the solitary-wave impact; dashed line at 0.8 dt_sim
(b) slowdown episodes after start-up (ticks), one lane per test
(c) share of windows at each latency step (stateOS resolution 0.49 ms), whole record and 5-95 % Arias
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
T_IMP = 347.4  # solitary-wave impact, t_int (median over completed tests; fig_fri_lanes.py)
TESTS = [("F24", "SSPQuad, elastic", "#56B4E9"), ("F23", "SSPQuad, PIMY stiff", "#009E73"),
         ("F05", "SSPQuad, PIMY soft", "#000000"), ("F27", "SSPQuad, PIMY soft (repeat)", "#8C8C8C"),
         ("F20", "Quad, PIMY soft", "#CC79A7")]
DT = R.DT_SIM * 1e3
BIN = 2.0

prs.configure_font()
D = {tid: R.load(tid) for tid, *_ in TESTS}
fig = plt.figure(figsize=(6.0, 5.2), layout="constrained")
gs = fig.add_gridspec(2, 2, height_ratios=[1.6, 0.8], width_ratios=[2.1, 1.0])
ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[1, 0], sharex=ax_a)
ax_c = fig.add_subplot(gs[:, 1])
for j, (tid, name, c) in enumerate(TESTS):
    d = D[tid]
    lat, ti = d["lat"][3:] * 1e3, d["ti_k"][3:]
    ok = np.isfinite(lat)
    edges = np.arange(0, ti[-1] + BIN, BIN)
    idx = np.digitize(ti[ok], edges)
    m = np.array([lat[ok][idx == i].mean() if np.any(idx == i) else np.nan for i in range(1, edges.size)])
    ax_a.plot(0.5 * (edges[1:] + edges[:-1]), m, color=c, lw=0.9 if tid != "F27" else 0.7, label=name, zorder=3)
    on = d["on"][d["on"] > d["t_su"]]
    ax_b.vlines(d["to_int"](on), j - 0.36, j + 0.36, color=c, lw=1.0)
    ax_b.add_patch(plt.Rectangle((0, j - 0.38), d["ti_k"][-1], 0.76, fc="0.92", ec="none", zorder=0))
    ax_b.text(d["ti_k"][-1] + 3, j, f"{d['n_slow']}", ha="left", va="center", fontsize=6.5, color="#B85C00")
    print(f"{tid} {name:28} mean latency {np.nanmean(lat):.2f} ms; slowdowns {d['n_slow']}")
for ax in (ax_a, ax_b):
    ax.axvspan(R.T_D5, R.T_D95, color="#d9e6f2", zorder=-1, lw=0)
    ax.axvline(T_IMP, color="0.3", lw=0.7, ls=(0, (4, 2)), zorder=1)
ax_a.axhline(0.8 * DT, color="0.35", lw=0.6, ls=(0, (4, 2)), zorder=0)
ax_a.text(5, 0.8 * DT - 0.05, r"$0.8\,\Delta t_{\mathrm{sim}}$", ha="left", va="top", color="0.35", fontsize=7)
ax_a.text(0.5 * (R.T_D5 + R.T_D95), 4.25, r"5--95\% Arias", ha="center", va="bottom", fontsize=7, color="0.3")
ax_a.text(T_IMP, 4.25, "solitary wave", ha="center", va="bottom", fontsize=7, color="0.3")
ax_a.set_ylim(1.0, 4.6)
ax_a.set_ylabel(f"Host latency, {BIN:.0f}-s mean (ms)")
ax_a.legend(loc="upper right", bbox_to_anchor=(1.0, 0.77), frameon=False, fontsize=7, handlelength=1.2)
ax_a.tick_params(labelbottom=False)
ax_b.set_yticks(range(len(TESTS)))
ax_b.set_yticklabels([t for t, *_ in TESTS], fontsize=6.5)
ax_b.set_ylim(len(TESTS) - 0.4, -0.6)
ax_b.set_xlim(0, 480)
ax_b.tick_params(axis="y", length=0)
ax_b.set_xlabel(r"Integration time, $t_{\mathrm{int}}$ (s)")
ax_b.set_ylabel("Slowdowns", fontsize=7)

TICK = 1e3 / R.FS  # stateOS resolution, 0.49 ms
EDGES = np.array([0, 2.5, 3.5, 4.5, 5.5, 6.5, 99]) * TICK  # bins centred on 2, 3, 4, 5, 6 ticks and above
BLAB = ["$\leq$0.98", "1.46", "1.95", "2.44", "2.93", "$\geq$3.42"]
BCOL = ["#deebf7", "#9ecae1", "#6baed6", "#3182bd", "#08519c", "#FF8C00"]
for j, (tid, name, c) in enumerate(TESTS):
    for k, (pn, f) in enumerate([("all", lambda t: t > 0), ("5--95\%", lambda t: (t >= R.T_D5) & (t < R.T_D95))]):
        d = D[tid]
        lat, ti = d["lat"][3:] * 1e3, d["ti_k"][3:]
        v = lat[f(ti) & np.isfinite(lat)]
        h = 100 * np.histogram(v, EDGES)[0] / v.size
        yy = j * 2.6 + k
        left = 0.0
        for b, (hb, cb) in enumerate(zip(h, BCOL)):
            ax_c.barh(yy, hb, left=left, height=0.8, color=cb, edgecolor="w", lw=0.3)
            left += hb
        ax_c.text(-2, yy, pn, ha="right", va="center", fontsize=6)
    ax_c.text(0, j * 2.6 - 0.62, f"{tid}: {name}", ha="left", va="bottom", fontsize=6.5, color=c)
ax_c.set_ylim(len(TESTS) * 2.6 - 1.1, -1.3)
ax_c.set_xlim(0, 100)
ax_c.set_yticks([])
ax_c.set_xlabel(r"Share of windows (\%)")
ax_c.legend([plt.Rectangle((0, 0), 1, 1, fc=cb, ec="none") for cb in BCOL], BLAB, title="Latency (ms)",
            loc="upper center", bbox_to_anchor=(0.45, -0.13), ncol=3, frameon=False, fontsize=6, title_fontsize=6.5,
            handlelength=0.9, columnspacing=0.8)
for sp in ("left", "right", "top"):
    ax_c.spines[sp].set_visible(False)
for ax, tg in ((ax_a, "(a)"), (ax_b, "(b)"), (ax_c, "(c)")):
    ax.set_title(rf"\textbf{{{tg}}}", loc="left")
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
