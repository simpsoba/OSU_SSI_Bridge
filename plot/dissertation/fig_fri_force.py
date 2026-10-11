"""Interface force in the Friday main campaign: hold transients versus the solitary wave (dissertation Ch. 6;
read-only on data).

lambda_L^3 f_p = -lambda^3 * topLRaw (load cell, as fig_wed_force.py), resisting-force convention, prototype scale,
on the prototype-scaled laboratory clock sqrt(lambda_L) t.
(a) slowdown lanes and force histories of the moderate model (F07 CPU ProfileSPD and F08 CPU MUMPS, both stopped
    early; F06 GPU) and of the large model with CPU integration and cuDSS solves (F15)
(b.1) solitary-wave impact, aligned at each test's peak, for GPU tests of increasing model size (F05, F06, F09, F13) and the CPU baseline (F03)
(b.2) an isolated slowdown of F11 (large model, GPU, 24 ranks) and the force transient it produced
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
CFG = {("CPU", "ProfileSPD"): "#C0390B", ("CPU", "MUMPS"): "#E69F00", ("CPU", "cuDSS"): "#0072B2",
       ("GPU", "cuDSS"): "#000000"}
MAIN = ["F07", "F08", "F15", "F06"]
SOL = [("F03", "#C0390B", "baseline, CPU"), ("F05", "#000000", "baseline"), ("F06", "#555555", "moderate"),
       ("F09", "#8C8C8C", "large"), ("F13", "#BBBBBB", "extra-large")]


def col(tid):
    r = R.BY_ID[tid]
    return CFG[(r[5], r[6].split()[0])]


def force(tid, t_end):
    m = sio.loadmat(R.DATA / "Simulink" / f"{R.BY_ID[tid][1]}.mat", squeeze_me=True, struct_as_record=False)["data"]
    raw = m.topLRaw.data
    k = raw[:, -1] <= t_end
    return raw[k, -1] * R.SQL, -R.LAM ** 3 * raw[k, 0] / 1e3


prs.configure_font()
ids = sorted(set(MAIN + [s[0] for s in SOL] + ["F11"]))
D = {tid: R.load(tid) for tid in ids}
F = {tid: force(tid, D[tid]["t"][-1]) for tid in ids}

fig = plt.figure(figsize=(6.0, 5.6), layout="constrained")
sa, sb = fig.subfigures(2, 1, height_ratios=[1.75, 1.0])
la, fa = sa.subplots(2, 1, sharex=True, gridspec_kw=dict(height_ratios=[0.5, 2.0]))
for j, tid in enumerate(MAIN):
    d = D[tid]
    on = d["on"][d["on"] > d["t_su"]] * R.SQL
    la.add_patch(plt.Rectangle((0, j - 0.38), d["t"][-1] * R.SQL, 0.76, fc="0.92", ec="none", zorder=0))
    la.vlines(on, j - 0.36, j + 0.36, color=col(tid), lw=0.4 if on.size > 200 else 0.9)
    la.text(482, j, f"{d['n_slow']}", ha="right", va="center", fontsize=6.5, color="#B85C00")
for tid in ["F06", "F15", "F08", "F07"]:
    t, f = F[tid]
    fa.plot(t, f, color=col(tid), lw=0.5 if tid != "F06" else 0.7, zorder=3 if tid == "F06" else 2)
la.set_yticks(range(len(MAIN)))
la.set_yticklabels(MAIN, fontsize=6.5)
for tl, tid in zip(la.get_yticklabels(), MAIN):
    tl.set_color(col(tid))
la.set_ylim(len(MAIN) - 0.5, -0.5)
la.tick_params(axis="y", length=0)
la.set_ylabel("Slowdowns", fontsize=7)
la.set_title(r"\textbf{(a)}", loc="left")
fa.set_xlim(0, 485)
fa.set_ylim(-145, 145)
fa.set_ylabel(r"$\lambda_L^3 f_p$ (kN)")
fa.set_xlabel(r"Laboratory time, $\sqrt{\lambda_L}\,t$ (s)")
for tid, xy, xyt in (("F07", (118, 105), (150, 125)), ("F08", (110, -124), (150, -132)),
                     ("F15", (44, -52), (150, -75)), ("F06", (347, -78), (380, -110))):
    r = R.BY_ID[tid]
    name = {"F07": "F07 CPU ProfileSPD (moderate, stopped)", "F08": "F08 CPU MUMPS (moderate, stopped)",
            "F15": "F15 CPU MKR + cuDSS (large)", "F06": "F06 GPU (moderate)"}[tid]
    fa.annotate(name, xy, xytext=xyt, color=col(tid), fontsize=6.5, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=col(tid), lw=0.5, shrinkA=1, shrinkB=1))
fa.text(347, 100, "solitary\nwave", ha="center", va="bottom", fontsize=6.5, color="0.3")

b1, b2 = sb.subplots(1, 2, gridspec_kw=dict(width_ratios=[1.4, 1.0]))
for tid, c, name in SOL:
    t, f = F[tid]
    w = (t > 335) & (t < 365)
    tp_ = t[w][np.argmax(np.abs(f[w]))]
    k = (t > tp_ - 1.5) & (t < tp_ + 3.0)
    b1.plot(t[k] - tp_, f[k], color=c, lw=0.8, label=f"{tid} {name}")
b1.legend(loc="lower right", frameon=False, fontsize=6, handlelength=1.0)
b1.set_xlim(-1.5, 3.0)
b1.set_xlabel(r"Time from the peak, $\sqrt{\lambda_L}\,\tau$ (s)")
b1.set_ylabel(r"$\lambda_L^3 f_p$ (kN)")
b1.set_title(r"\textbf{(b.1)}", loc="left")
d = D["F11"]
t, f = F["F11"]
on = d["on"][d["on"] > d["t_su"]] * R.SQL
tp = t[np.argmax(np.abs(f))]
k = (t > tp - 0.8) & (t < tp + 1.2)
b2.plot(t[k], f[k], color="#000000", lw=0.8)
for x in on[np.abs(on - tp) < 1.0]:
    b2.axvline(x, color="#FF8C00", lw=0.8, zorder=0)
b2.text(0.97, 0.05, "F11 GPU (large, 24 ranks)\norange: slowdown onset", transform=b2.transAxes, ha="right",
        va="bottom", fontsize=6.5)
b2.set_xlabel(r"$\sqrt{\lambda_L}\,t$ (s)")
b2.set_title(r"\textbf{(b.2)}", loc="left")
print(f"F11 peak {f[np.argmax(np.abs(f))]:.1f} kN at {tp:.2f} s; slowdown onsets nearby {np.round(on[np.abs(on - tp) < 1], 3)}")
for tid, c, name in SOL:
    t, f = F[tid]
    k = (t > 335) & (t < 365)
    print(f"{tid} solitary peak {f[k][np.argmax(np.abs(f[k]))]:.1f} kN")
for tid in MAIN:
    t, f = F[tid]
    print(f"{tid} max |force| {np.abs(f).max():.1f} kN")
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
