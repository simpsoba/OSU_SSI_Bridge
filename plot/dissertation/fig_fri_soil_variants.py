"""Soil-model variants on the baseline model, main campaign (dissertation Ch. 6 draft; read-only on data).

All GPU MKR-alpha + cuDSS, 8 ranks, baseline mesh: F05 (SSPquad + PIMY, the baseline), F20 (quad + PIMY),
F24 (SSPquad + elastic soil); PIMY = PressureIndependMultiYield (clay profile). Measured interface displacement u = meaSigOS * lambda (mm) on
integration time (slowdowns removed through the received-target count, as fri_runs.load). F23 (stiff soil) and F30 (2 % Rayleigh) not shown.
(a) u(t_int) during the earthquake; (b) peak |u| during the earthquake and the residual u (mean over 250-300 s);
(c) host latency per dt_sim window after start-up, relative to dt_sim, as fig_fri_capacity.py (marker = mean,
    faint box = interquartile range, capped bar = 1st to 99th percentiles; dashed = 0.8 dt_sim), and the slowdown
    episodes after start-up.
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
V = [("F05", "SSPquad + PIMY", "#000000", "-", 1.0), ("F20", "quad + PIMY", "0.5", (0, (3, 1.5)), 1.2),
     ("F24", "SSPquad + elastic", "#C0390B", (0, (1, 1.2)), 1.3)]
tg = np.arange(2, 420, 0.01)
D = {}
for tid, *_ in V:
    r, d = R.BY_ID[tid], R.load(tid)
    me = np.load(R.MAT_EXTRACT_DIR / f"{r[1]}.npz", allow_pickle=True)["meaSigOS_data"]
    u = np.interp(tg, d["to_int"](me[:, -1]), me[:, 0] * R.LAM * 1e3)
    D[tid] = dict(u=u, peak=np.abs(u[tg < 139]).max(), res=u[(tg > 250) & (tg < 300)].mean(),
                  lat=d["lat"][3:][np.isfinite(d["lat"][3:])] / R.DT_SIM, nsl=d["n_slow"])
    print(f"{tid}: peak {D[tid]['peak']:.1f} mm, residual {D[tid]['res']:.1f} mm, latency mean {D[tid]['lat'].mean():.3f}, "
          f"p99 {np.percentile(D[tid]['lat'], 99):.3f} dt_sim, "
          f"slowdowns {D[tid]['nsl']}")

prs.configure_font()
fig = plt.figure(figsize=(6.0, 3.9), layout="constrained")
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.8])
ax = fig.add_subplot(gs[0, :])
for z, (tid, lab, c, ls, lw) in enumerate(V):  # F05 first, so the dashed F20 shows on top of it
    ax.plot(tg, D[tid]["u"], color=c, ls=ls, lw=lw, zorder=2 + z, label=f"{tid}, {lab}")
ax.set_xlim(40, 160)
ax.set_ylim(-300, 215)
ax.set_xlabel(r"Integration time, $t_{\mathrm{int}}$ (s)")
ax.set_ylabel(r"$u$ (mm)")
ax.set_title(r"\textbf{(a)} interface displacement during the earthquake", loc="left")
lg = ax.legend(loc="lower right", ncol=2, frameon=False, handlelength=2.6, columnspacing=1.2, borderaxespad=0.2)
for t_, (tid, lab, c, *_) in zip(lg.get_texts(), V):
    t_.set_color(c)
bx = fig.add_subplot(gs[1, 0])
cx = fig.add_subplot(gs[1, 1])
y = np.arange(len(V))[::-1]
for yi, (tid, lab, c, *_) in zip(y, V):
    bx.plot(D[tid]["peak"], yi, "o", color=c, ms=5)
    bx.plot(D[tid]["res"], yi, "s", color=c, mfc="white", mew=1.0, ms=5)
    lat = D[tid]["lat"]  # as fig_fri_capacity.py: mean, interquartile box, 1st-99th percentile bar
    p1, p99 = np.percentile(lat, [1, 99])
    q1, q3 = np.percentile(lat, [25, 75])
    cx.errorbar(0.5 * (p1 + p99), yi, xerr=0.5 * (p99 - p1), fmt="none", ecolor=c, elinewidth=0.6, capsize=2.2,
                capthick=0.6, alpha=0.6, zorder=2)
    cx.add_patch(plt.Rectangle((q1, yi - 0.12), max(q3 - q1, 0.006), 0.24, fc=c, ec="none", alpha=0.22, zorder=3))
    cx.plot(lat.mean(), yi, "o", ms=4.6, color=c, zorder=4)
    cx.text(1.17, yi, f"{D[tid]['nsl']}", ha="right", va="center", color="#B85C00")
for a_ in (bx, cx):
    a_.set_yticks(y)
    a_.set_ylim(-0.6, len(V) - 0.4)
bx.set_yticklabels([tid for tid, *_ in V])
cx.set_yticklabels([])
bx.set_xlabel(r"$u$ (mm)")
bx.set_xlim(0, 210)
bx.text(0.98, 1.02, "filled: peak; open: residual", transform=bx.transAxes, ha="right", va="bottom", color="0.3")
bx.set_title(r"\textbf{(b)}", loc="left")
cx.axvline(0.8, color="0.35", lw=0.6, ls=(0, (4, 2)), zorder=0)
cx.text(0.79, len(V) - 0.5, r"$0.8\,\Delta t_{\mathrm{sim}}$", ha="right", va="bottom", color="0.35")
cx.set_xlabel(r"Host latency / $\Delta t_{\mathrm{sim}}$")
cx.set_xlim(0, 1.2)
cx.set_title(r"\textbf{(c)}", loc="left")
cx.text(0.97, 1.02, "slowdowns", transform=cx.transAxes, ha="right", va="bottom", color="#B85C00")
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
