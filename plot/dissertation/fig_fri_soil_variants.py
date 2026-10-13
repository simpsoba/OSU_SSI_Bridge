"""Soil-model variants on the baseline model, main campaign (dissertation Ch. 6 draft; read-only on data).

All GPU MKR-alpha + cuDSS, 8 ranks, baseline mesh: F05 (SSPquad, nonlinear soil), F20 (quad + PIMY), F23 (stiff
soil), F24 (elastic soil), F30 (2 % Rayleigh damping). Measured interface displacement u = meaSigOS * lambda (mm) on
integration time (slowdowns removed through the received-target count, as fri_runs.load).
(a) u(t_int) during the earthquake; (b) peak |u| during the earthquake and the residual u (mean over 250-300 s);
(c) mean host latency as a share of dt_sim, with the slowdown episodes after start-up.
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
V = [("F05", "SSPquad, nonlinear", "#000000", "-", 1.0), ("F20", "quad + PIMY", "0.55", (0, (4, 1.5)), 0.9),
     ("F30", r"2\% Rayleigh", "0.75", "-", 0.8), ("F23", "stiff soil", "#2166AC", "-", 0.9),
     ("F24", "elastic soil", "#C0390B", "-", 0.9)]
tg = np.arange(2, 420, 0.01)
D = {}
for tid, *_ in V:
    r, d = R.BY_ID[tid], R.load(tid)
    me = np.load(R.MAT_EXTRACT_DIR / f"{r[1]}.npz", allow_pickle=True)["meaSigOS_data"]
    u = np.interp(tg, d["to_int"](me[:, -1]), me[:, 0] * R.LAM * 1e3)
    D[tid] = dict(u=u, peak=np.abs(u[tg < 139]).max(), res=u[(tg > 250) & (tg < 300)].mean(),
                  lat=100 * np.nanmean(d["lat"][3:]) / R.DT_SIM, nsl=d["n_slow"])
    print(f"{tid}: peak {D[tid]['peak']:.1f} mm, residual {D[tid]['res']:.1f} mm, latency {D[tid]['lat']:.1f} %, "
          f"slowdowns {D[tid]['nsl']}")

prs.configure_font()
fig = plt.figure(figsize=(6.0, 3.9), layout="constrained")
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.8])
ax = fig.add_subplot(gs[0, :])
for tid, lab, c, ls, lw in V[::-1]:
    if tid == "F30":
        continue
    ax.plot(tg, D[tid]["u"], color=c, ls=ls, lw=lw)
ax.set_xlim(40, 160)
ax.set_ylim(-180, 290)
ax.set_xlabel(r"Integration time, $t_{\mathrm{int}}$ (s)")
ax.set_ylabel(r"$u$ (mm)")
ax.set_title(r"\textbf{(a)} interface displacement during the earthquake", loc="left")
for i, (tid, lab, c, ls, lw) in enumerate([v for v in V if v[0] != "F30"]):
    ax.text(0.02, 0.97 - 0.11 * i, f"{tid}, {lab}", transform=ax.transAxes, ha="left", va="top", color=c)
bx = fig.add_subplot(gs[1, 0])
cx = fig.add_subplot(gs[1, 1])
y = np.arange(len(V))[::-1]
for yi, (tid, lab, c, *_) in zip(y, V):
    bx.plot(D[tid]["peak"], yi, "o", color=c, ms=5)
    bx.plot(D[tid]["res"], yi, "s", color=c, mfc="white", mew=1.0, ms=5)
    cx.barh(yi, D[tid]["lat"], color=c, height=0.6)
    cx.text(D[tid]["lat"] + 0.8, yi, f"{D[tid]['nsl']}", va="center", color="#B85C00")
for a_ in (bx, cx):
    a_.set_yticks(y)
    a_.set_ylim(-0.6, len(V) - 0.4)
bx.set_yticklabels([f"{tid}, {lab}" for tid, lab, *_ in V])
cx.set_yticklabels([])
bx.set_xlabel(r"$u$ (mm)")
bx.set_xlim(0, 210)
bx.text(0.98, 1.02, "filled: peak; open: residual", transform=bx.transAxes, ha="right", va="bottom", color="0.3")
bx.set_title(r"\textbf{(b)}", loc="left")
cx.set_xlabel(r"Mean host latency / $\Delta t_{\mathrm{sim}}$ (\%)")
cx.set_xlim(0, 50)
cx.set_title(r"\textbf{(c)}", loc="left")
cx.text(0.97, 1.02, "slowdowns", transform=cx.transAxes, ha="right", va="bottom", color="#B85C00")
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
