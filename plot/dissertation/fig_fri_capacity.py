"""Computational capacity map of the Friday main campaign (dissertation Ch. 6; read-only on data).

x = model size (DOFs; points offset by configuration for legibility).
(a) mean host latency after start-up, relative to dt_sim; dashed line at 0.8 (slowdown threshold)
(b) slowdown episodes after start-up per minute of laboratory time (model scale); symlog axis, zero at the bottom
Markers by configuration (filled: completed; open: stopped early, labeled); numbers next to markers = MPI ranks.
Soil-formulation variants (F20, F23, F24) are omitted here; see fig_fri_formulation.py.
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
CFG = {("CPU", "ProfileSPD"): ("#C0390B", "s", -0.30, "CPU MKR-$\\alpha$ + ProfileSPD"),
       ("CPU", "MUMPS"): ("#E69F00", "D", -0.15, "CPU MKR-$\\alpha$ + MUMPS"),
       ("CPU", "cuDSS"): ("#0072B2", "^", 0.0, "CPU MKR-$\\alpha$ + cuDSS solves"),
       ("GPU", "cuDSS"): ("#000000", "o", 0.20, "GPU MKR-$\\alpha$ + cuDSS")}
IDS = ["F03", "F04", "F26", "F22", "F05", "F27", "F07", "F08", "F06", "F09", "F10", "F12", "F11", "F14", "F15", "F13"]
SIZES = [2028, 3216, 4620, 5484]
NAMES = ["Baseline", "Moderate", "Large", "Extra-large"]

prs.configure_font()
fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(6.0, 2.9), layout="constrained")
seen = {}
rows = []
for tid in IDS:
    r = R.BY_ID[tid]
    d = R.load(tid)
    c, mk, dx, lab = CFG[(r[5], r[6].split()[0])]
    xi = SIZES.index(r[4])
    lat = d["lat"][3:]
    lat = lat[np.isfinite(lat)].mean() / R.DT_SIM
    mins = (d["t"][-1] - d["t_su"]) / 60.0
    rate = d["n_slow"] / mins
    rows.append((tid, r[4], lab, r[7], lat, rate))
    jit = {"F04": 0.07, "F27": 0.07, "F10": 0.08, "F12": 0.16, "F11": 0.24, "F15": 0.07}.get(tid, 0.0)
    for ax, y in ((ax_a, lat), (ax_b, rate)):
        h = ax.plot(xi + dx + jit, y, marker=mk, ms=4.6, color=c, mfc=c if d["complete"] else "white", mew=0.9,
                    ls="none", zorder=4)[0]
        ax.annotate(f"{r[7]}", (xi + dx + jit, y), xytext=(3, 0), textcoords="offset points", fontsize=6,
                    color=c, ha="left", va="center")
    if d["complete"]:
        seen.setdefault(lab, h)
ax_b.text(1.0, 1900, "stopped early (open)", ha="center", va="bottom", fontsize=6.5, color="0.35")
ax_a.axhline(0.8, color="0.35", lw=0.6, ls=(0, (4, 2)), zorder=0)
ax_a.text(-0.45, 0.82, r"slowdown threshold, $0.8\,\Delta t_{\mathrm{sim}}$", ha="left", va="bottom", color="0.35",
          fontsize=7)
ax_a.set_ylim(0, 1.0)
ax_a.set_ylabel(r"Mean host latency / $\Delta t_{\mathrm{sim}}$")
ax_b.set_yscale("symlog", linthresh=0.1, linscale=0.4)
ax_b.set_ylim(-0.02, 4000)
ax_b.set_yticks([0, 0.1, 1, 10, 100, 1000])
ax_b.set_yticklabels(["0", "0.1", "1", "10", "100", "1000"])
ax_b.set_ylabel("Slowdown episodes per minute")
for ax, tg in ((ax_a, "(a)"), (ax_b, "(b)")):
    ax.set_xticks(range(4))
    ax.set_xticklabels([f"{n}\n{s}" for n, s in zip(NAMES, SIZES)], fontsize=7, linespacing=0.95)
    ax.set_xlim(-0.55, 3.55)
    ax.set_xlabel("Model (DOFs)")
    ax.set_title(rf"\textbf{{{tg}}}", loc="left")
fig.legend(list(seen.values()), list(seen.keys()), loc="outside upper center", ncol=2, frameon=False,
           handletextpad=0.3, columnspacing=1.2)
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
for row in rows:
    print(f"{row[0]} {row[1]:5d} {row[2]:32} np{row[3]:2d} latency/dt_sim {row[4]:.2f} slowdowns/min {row[5]:7.2f}")
print("wrote", OUT.with_suffix(".png"))
