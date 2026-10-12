"""Interface force of the baseline GPU tests F05, F27 (8 ranks) and F22 (serial, stopped), Friday (dissertation Ch. 6;
read-only on data). Check of whether the slowdowns of F22 produce force spikes as in the CPU tests.

(a) slowdown lanes and force histories on the prototype-scaled laboratory clock (F22 drawn on top)
(b) zoom on the largest |force| of F22 with its slowdown onsets
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

sys.path.insert(0, str(Path(__file__).parent))
import fri_runs as R  # noqa: E402

sys.path.insert(0, str(R.ROOT / "plot"))
import PlotResponseSpectrum as prs  # noqa: E402

OUT = Path(sys.argv[1])
TESTS = [("F05", "#000000", "F05 GPU, 8 ranks"), ("F27", "#8C8C8C", "F27 GPU, 8 ranks (repeat)"),
         ("F22", "#C0390B", "F22 GPU, serial (stopped)")]


def force(tid, t_end):
    m = sio.loadmat(R.DATA / "Simulink" / f"{R.BY_ID[tid][1]}.mat", squeeze_me=True, struct_as_record=False)["data"]
    raw = m.topLRaw.data
    k = raw[:, -1] <= t_end
    return raw[k, -1] * R.SQL, -R.LAM ** 3 * raw[k, 0] / 1e3


prs.configure_font()
D = {tid: R.load(tid) for tid, *_ in TESTS}
F = {tid: force(tid, D[tid]["t"][-1]) for tid, *_ in TESTS}
fig = plt.figure(figsize=(6.0, 4.8), layout="constrained")
sa, sb = fig.subfigures(2, 1, height_ratios=[1.5, 1.0])
la, fa = sa.subplots(2, 1, sharex=True, gridspec_kw=dict(height_ratios=[0.45, 2.0]))
for j, (tid, c, name) in enumerate(TESTS):
    d = D[tid]
    on = d["on"][d["on"] > d["t_su"]] * R.SQL
    la.add_patch(plt.Rectangle((0, j - 0.38), d["t"][-1] * R.SQL, 0.76, fc="0.92", ec="none", zorder=0))
    la.vlines(on, j - 0.36, j + 0.36, color=c, lw=0.5)
    la.text(452, j, f"{d['n_slow']}", ha="right", va="center", fontsize=6.5, color="#B85C00")
    t, f = F[tid]
    fa.plot(t, f, color=c, lw=0.5, label=name, zorder={"F22": 3, "F05": 2, "F27": 1}[tid])
    eq = t < 300
    print(f"{tid}: max |force| during EQ (lab < 300 s) {np.abs(f[eq]).max():.1f} kN; overall {np.abs(f).max():.1f} kN; "
          f"slowdowns {d['n_slow']}")
la.set_yticks(range(3))
la.set_yticklabels([t for t, *_ in TESTS], fontsize=6.5)
la.set_ylim(2.5, -0.5)
la.tick_params(axis="y", length=0)
la.set_title(r"\textbf{(a)}", loc="left")
fa.set_xlim(0, 455)
fa.set_ylim(-100, 130)
fa.set_ylabel(r"$\lambda_L^3 f_p$ (kN)")
fa.set_xlabel(r"Laboratory time, $\sqrt{\lambda_L}\,t$ (s)")
fa.legend(loc="upper right", frameon=False, fontsize=6.5, handlelength=1.2)

bx = sb.subplots()
t, f = F["F22"]
tp = t[np.argmax(np.abs(f))]
k = (t > tp - 3) & (t < tp + 3)
for tid, c, _ in TESTS[:2]:
    t2, f2 = F[tid]
    k2 = (t2 > tp - 3) & (t2 < tp + 3)
    bx.plot(t2[k2], f2[k2], color=c, lw=0.7)
bx.plot(t[k], f[k], color="#C0390B", lw=0.8)
on = D["F22"]["on"] * R.SQL
bx.vlines(on[np.abs(on - tp) < 3], 0.93, 1.0, transform=bx.get_xaxis_transform(), color="#FF8C00", lw=0.8)
bx.set_xlim(tp - 3, tp + 3)
bx.set_xlabel(r"$\sqrt{\lambda_L}\,t$ (s)")
bx.set_ylabel(r"$\lambda_L^3 f_p$ (kN)")
bx.set_title(r"\textbf{(b)} Largest F22 transient (ticks: F22 slowdown onsets)", loc="left")
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
