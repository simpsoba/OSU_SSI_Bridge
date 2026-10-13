"""MPI-rank sweeps with GPU MKR-alpha + cuDSS in the main campaign (dissertation Ch. 6; read-only on data).

Moderate, large, and extra-large models; RTHS tests (filled) and lunch dry runs without actuator motion (open:
OpenSees, OpenFresco, and the rate transition running, actuator disconnected).
(a) host latency per dt_sim window after start-up, relative to dt_sim: marker = mean, capped bar = 1st to 99th
    percentiles; dashed line = slowdown threshold 0.8 dt_sim
(b) slowdown episodes after start-up per minute of laboratory time (model scale), symlog axis
Dotted vertical line: 24 physical cores (8 performance + 16 efficiency) of the host processor.
Direct model labels in (b); uniform 9 pt font; bold panel tags; constrained layout (style of fig_fri_capacity.py).
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
SER = [("Moderate, 3,216 DOFs", "#8C510A", "D", ["F06", "Fd01"]),
       ("Large, 4,620 DOFs", "#2166AC", "^", ["F09", "F10", "Fd02", "Fd03", "F12", "Fd04", "F11"]),
       ("Extra-large, 5,484 DOFs", "#000000", "o", ["Fd05", "Fd06", "F13", "Fd07", "Fd08", "Fd09"])]
OFF = {"Fd07": 0.35, "F13": -0.35, "Fd01": 0.35, "F06": 0.0}  # small x shifts where wet and dry share a rank count
LAB_B = {"Moderate, 3,216 DOFs": (10.0, 0.13), "Large, 4,620 DOFs": (10.0, 6.5), "Extra-large, 5,484 DOFs": (8.6, 45.0)}

prs.configure_font()
fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(6.0, 2.9), layout="constrained")
for name, c, mk, ids in SER:
    pts = []
    for tid in ids:
        r, d = R.BY_ID[tid], R.load(tid)
        lat = d["lat"][3:]
        lat = lat[np.isfinite(lat)] / R.DT_SIM
        rate = d["n_slow"] / ((d["t"][-1] - d["t_su"]) / 60.0)
        dry = r[8] == "dry"
        x = r[7] + OFF.get(tid, 0.0)
        p1, p99 = np.percentile(lat, [1, 99])
        pts.append((x, lat.mean(), rate))
        ax_a.errorbar(x, 0.5 * (p1 + p99), yerr=0.5 * (p99 - p1), fmt="none", ecolor=c, elinewidth=0.6, capsize=2.2,
                      capthick=0.6, alpha=0.6, zorder=2)
        for ax, y in ((ax_a, lat.mean()), (ax_b, rate)):
            ax.plot(x, y, marker=mk, ms=4.6, color=c, mfc="white" if dry else c, mew=0.9, ls="none", zorder=4)
        print(f"{tid} {name:24} {r[7]:2d} ranks {'dry' if dry else 'RTHS'}: latency/dt {lat.mean():.2f} "
              f"(p1 {p1:.2f}, p99 {p99:.2f}); slowdowns/min {rate:.2f}")
    pts.sort()
    for ax, k in ((ax_a, 1), (ax_b, 2)):
        ax.plot([p[0] for p in pts], [p[k] for p in pts], color=c, lw=0.6, alpha=0.5, zorder=1)
    ax_b.text(*LAB_B[name], name, color=c, ha="left", va="center")
ax_a.axhline(0.8, color="0.35", lw=0.6, ls=(0, (4, 2)), zorder=0)
ax_a.text(5.4, 0.82, r"$0.8\,\Delta t_{\mathrm{sim}}$", ha="left", va="bottom", color="0.35")
ax_a.set_ylim(0, 1.0)
ax_a.set_ylabel(r"Host latency / $\Delta t_{\mathrm{sim}}$")
ax_b.set_yscale("symlog", linthresh=0.1, linscale=0.4)
ax_b.set_ylim(-0.02, 400)
ax_b.set_yticks([0, 0.1, 1, 10, 100])
ax_b.set_yticklabels(["0", "0.1", "1", "10", "100"])
ax_b.set_ylabel("Slowdown episodes per minute")
ax_b.plot([], [], marker="o", color="0.3", mfc="0.3", ls="none", ms=4.6, label="RTHS test")
ax_b.plot([], [], marker="o", color="0.3", mfc="white", ls="none", ms=4.6, label="dry run")
ax_b.legend(loc="upper left", frameon=False, handletextpad=0.3, borderaxespad=0.2)
for ax, tg in ((ax_a, "(a)"), (ax_b, "(b)")):
    ax.axvline(24, color="0.55", lw=0.6, ls=":", zorder=0)
    ax.set_xlim(5, 30)
    ax.set_xticks([8, 12, 16, 20, 24, 28])
    ax.set_xlabel("MPI ranks")
    ax.set_title(rf"\textbf{{{tg}}}", loc="left")
ax_b.text(24.4, 0.15, "24 physical\ncores", ha="left", va="center", color="0.45", linespacing=0.95)
fig.savefig(OUT.with_suffix(".png"), dpi=250)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
