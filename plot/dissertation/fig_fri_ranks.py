"""Rank sweeps with GPU MKR-alpha + cuDSS, Friday (dissertation Ch. 6; read-only on data).

Wet RTHS tests (filled) and lunch dry runs without actuator motion (open) of the moderate, large, and extra-large
models. (a) mean host latency after start-up relative to dt_sim, with the 1st-99th percentile range as a bar;
(b) slowdown episodes after start-up per minute of laboratory time (symlog). Vertical line: 24 physical cores.
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
SER = [("Moderate (3216 DOFs)", "#E69F00", "s", ["F06", "Fd01"]),
       ("Large (4620 DOFs)", "#0072B2", "^", ["F09", "F10", "Fd02", "Fd03", "F12", "Fd04", "F11"]),
       ("Extra-large (5484 DOFs)", "#000000", "o", ["Fd05", "Fd06", "F13", "Fd07", "Fd08", "Fd09"])]

prs.configure_font()
fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(6.0, 2.8), layout="constrained")
for name, c, mk, ids in SER:
    pts = []
    for tid in ids:
        r, d = R.BY_ID[tid], R.load(tid)
        lat = d["lat"][3:]
        lat = lat[np.isfinite(lat)] / R.DT_SIM
        rate = d["n_slow"] / ((d["t"][-1] - d["t_su"]) / 60.0)
        dry = r[8] == "dry"
        off = 0.25 if dry and tid == "Fd07" else 0.0
        pts.append((r[7], lat.mean(), rate))
        ax_a.vlines(r[7] + off, *np.percentile(lat, [1, 99]), color=c, lw=0.6, alpha=0.6)
        for ax, y in ((ax_a, lat.mean()), (ax_b, rate)):
            ax.plot(r[7] + off, y, marker=mk, ms=4.5, color=c, mfc="white" if dry else c, mew=0.9, ls="none")
        print(f"{tid} {name:24} np{r[7]:2d} {'dry' if dry else 'wet'} latency/dt {lat.mean():.2f} "
              f"(p99 {np.percentile(lat, 99):.2f}) slowdowns/min {rate:.2f}")
    pts.sort()
    for ax, k in ((ax_a, 1), (ax_b, 2)):
        ax.plot([p[0] for p in pts], [p[k] for p in pts], color=c, lw=0.6, alpha=0.5, zorder=0)
    ax_b.plot([], [], marker=mk, color=c, ls="-", lw=0.6, label=name)
ax_b.plot([], [], marker="o", color="0.4", mfc="0.4", ls="none", label="RTHS test")
ax_b.plot([], [], marker="o", color="0.4", mfc="white", ls="none", label="dry run")
ax_a.axhline(0.8, color="0.35", lw=0.6, ls=(0, (4, 2)))
ax_a.text(6.5, 0.82, r"$0.8\,\Delta t_{\mathrm{sim}}$", color="0.35", fontsize=7, va="bottom")
ax_a.set_ylim(0, 1.0)
ax_a.set_ylabel(r"Host latency / $\Delta t_{\mathrm{sim}}$")
ax_b.set_yscale("symlog", linthresh=0.1, linscale=0.4)
ax_b.set_ylim(-0.02, 300)
ax_b.set_yticks([0, 0.1, 1, 10, 100])
ax_b.set_yticklabels(["0", "0.1", "1", "10", "100"])
ax_b.set_ylabel("Slowdown episodes per minute")
ax_b.legend(loc="upper left", frameon=False, fontsize=6.5, handlelength=1.4)
for ax, tg in ((ax_a, "(a)"), (ax_b, "(b)")):
    ax.axvline(24, color="0.6", lw=0.6, ls=":")
    ax.set_xlim(5, 30)
    ax.set_xticks([8, 12, 16, 20, 24, 28])
    ax.set_xlabel("MPI ranks")
    ax.set_title(rf"\textbf{{{tg}}}", loc="left")
ax_a.text(24.3, 0.05, "24 cores", fontsize=6.5, color="0.4", rotation=90, va="bottom")
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
