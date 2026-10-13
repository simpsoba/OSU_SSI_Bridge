"""Computational capacity of the main campaign (dissertation Ch. 6; read-only on data).

x = model (baseline, moderate, large, extra-large: 2,028, 3,216, 4,620, 5,484 DOFs); tests offset within each model by
configuration and rank count. Configurations: CPU MKR-alpha with ProfileSPD, with MUMPS, or with cuDSS solves
(integration on the CPU, linear solves on the GPU), and GPU MKR-alpha with cuDSS.
(a) host latency per dt_sim window after start-up, relative to dt_sim: marker = mean, thin line = 1st to 99th
    percentiles; dashed line = slowdown threshold 0.8 dt_sim
(b) slowdown episodes after start-up per minute of laboratory time (model scale), symlog axis with zero at the bottom
Filled markers: completed tests; open markers: tests stopped early. Numbers beside the markers: MPI ranks.
Key above the panels (four configurations, color and marker); bold panel tags; uniform font; constrained layout (style of the preliminary-campaign figures).
Soil-formulation variants and dry runs are not shown here.
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
CFG = {"pspd": ("#C0390B", "s", r"CPU MKR-$\alpha$, ProfileSPD"),
       "mumps": ("#8C510A", "D", r"CPU MKR-$\alpha$, MUMPS"),
       "cpucudss": ("#2166AC", "^", r"CPU MKR-$\alpha$, cuDSS solves"),
       "gpu": ("#000000", "o", r"GPU MKR-$\alpha$, cuDSS")}
# test: (configuration, x offset within its model)
TESTS = {"F03": ("pspd", -0.32), "F04": ("pspd", -0.20), "F26": ("mumps", -0.05),
         "F22": ("gpu", 0.10), "F05": ("gpu", 0.22), "F27": ("gpu", 0.34),
         "F07": ("pspd", -0.24), "F08": ("mumps", -0.06), "F06": ("gpu", 0.22),
         "F14": ("cpucudss", -0.30), "F15": ("cpucudss", -0.17), "F09": ("gpu", 0.00), "F10": ("gpu", 0.12),
         "F12": ("gpu", 0.24), "F11": ("gpu", 0.36),
         "F13": ("gpu", 0.10)}
SIZES = [2028, 3216, 4620, 5484]
NAMES = ["Baseline", "Moderate", "Large", "Extra-large"]
prs.configure_font()
fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(6.0, 2.9), layout="constrained")
pts = {}
for tid, (cfg, dx) in TESTS.items():
    r, d = R.BY_ID[tid], R.load(tid)
    c, mk, _ = CFG[cfg]
    x = SIZES.index(r[4]) + dx
    lat = d["lat"][3:]
    lat = lat[np.isfinite(lat)] / R.DT_SIM
    rate = d["n_slow"] / ((d["t"][-1] - d["t_su"]) / 60.0)
    pts[tid] = (x, lat.mean(), rate)
    mfc = c if d["complete"] else "white"
    ax_a.vlines(x, *np.percentile(lat, [1, 99]), color=c, lw=0.6, alpha=0.55, zorder=2)
    for ax, y in ((ax_a, lat.mean()), (ax_b, rate)):
        ax.plot(x, y, marker=mk, ms=4.6, color=c, mfc=mfc, mew=0.9, ls="none", zorder=4)
    rk = "serial" if r[7] == 1 else f"{r[7]}"
    ax_b.annotate(rk, (x, rate), xytext=(0, 5), textcoords="offset points", ha="center", va="bottom", fontsize=6,
                  color=c)
    print(f"{tid} {cfg:9} ranks {rk:>6}: latency/dt mean {lat.mean():.2f} (p1 {np.percentile(lat, 1):.2f}, "
          f"p99 {np.percentile(lat, 99):.2f}); slowdowns/min {rate:8.2f}; {'complete' if d['complete'] else 'stopped'}")
handles = [plt.Line2D([], [], marker=CFG[k][1], color=CFG[k][0], ls="none", ms=4.6, label=CFG[k][2]) for k in CFG]
fig.legend(handles=handles, loc="outside upper center", ncol=4, frameon=False, handletextpad=0.3, columnspacing=1.0,
           fontsize=7.5)
ax_a.axhline(0.8, color="0.35", lw=0.6, ls=(0, (4, 2)), zorder=0)
ax_a.text(3.5, 0.82, r"$0.8\,\Delta t_{\mathrm{sim}}$", ha="right", va="bottom", color="0.35")
ax_a.set_ylim(0, 1.05)
ax_a.set_ylabel(r"Host latency / $\Delta t_{\mathrm{sim}}$")
ax_b.set_yscale("symlog", linthresh=0.1, linscale=0.4)
ax_b.set_ylim(-0.02, 9000)
ax_b.set_yticks([0, 0.1, 1, 10, 100, 1000])
ax_b.set_yticklabels(["0", "0.1", "1", "10", "100", "1000"])
ax_b.set_ylabel("Slowdown episodes per minute")
ax_b.text(0.02, 0.98, "open markers: stopped early", transform=ax_b.transAxes, ha="left", va="top", fontsize=6.5,
          color="0.35")
for ax, tg in ((ax_a, "(a)"), (ax_b, "(b)")):
    ax.set_xticks(range(4))
    ax.set_xticklabels([f"{n}\n{s:,}" for n, s in zip(NAMES, SIZES)], linespacing=0.95)
    ax.set_xlim(-0.55, 3.6)
    ax.set_xlabel("Model (DOFs)")
    ax.set_title(rf"\textbf{{{tg}}}", loc="left")
fig.savefig(OUT.with_suffix(".png"), dpi=250)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
