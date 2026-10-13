"""Interface force in the main campaign: force spikes from computational slowdowns, compared with the solitary-wave
load (dissertation Ch. 6; read-only on data).

lambda_L^3 f_p = -lambda^3 * topLRaw (load cell, resisting-force convention, prototype scale) on the prototype-scaled
laboratory clock, as fig_wed_force.py; dashed lines at +-0.1 F_y,eq (47.7 kN).
(a) baseline model, GPU MKR-alpha + cuDSS: F05 and F27 (8 ranks) and F22 (serial, NC), with slowdown lanes above
(b) moderate model: F07 (CPU MKR-alpha + ProfileSPD, NC) and F08 (CPU MKR-alpha + MUMPS, NC) and F06 (GPU, 8 ranks)
(c) isolated slowdown of F11 (large model, GPU, 24 ranks) and the force transient it produced; tick = slowdown onset
(d) solitary-wave pulse of the completed GPU tests on the four models, aligned at each peak (tau = 0)
Configuration colors and markers as fig_fri_capacity.py; uniform 9 pt font; bold panel tags; constrained layout.
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
CFG = {"pspd": ("#C0390B", "s"), "mumps": ("#8C510A", "D"), "cpucudss": ("#2166AC", "^"), "gpu": ("#000000", "o")}
ROWS = [("Baseline, GPU MKR-$\\alpha$", [("F05", "gpu", "0.45", "F05, 8"), ("F27", "gpu", "0.72", "F27, 8"),
                                         ("F22", "gpu", "#000000", "F22, serial, NC")]),
        ("Moderate", [("F06", "gpu", "0.0", r"F06, GPU MKR-$\alpha$, 8"),
                      ("F08", "mumps", "#8C510A", r"F08, CPU MKR-$\alpha$ + MUMPS, 8, NC"),
                      ("F07", "pspd", "#C0390B", r"F07, CPU MKR-$\alpha$ + ProfileSPD, 8, NC")])]
SOL = [("F05", "baseline", "0.0"), ("F06", "moderate", "#8C510A"), ("F09", "large", "#2166AC"),
       ("F13", "extra-large", "0.55")]
YL = 140.0


def fy_eq_kn():
    inch = 0.0254
    d_pier, h_pier = 48 * inch, 3.02 * (48.0 / 20.0)
    as_tot = 28 * (np.pi / 4) * (1.25 * inch) ** 2
    r_bar = 0.5 * d_pier - 2 * inch - 0.5 * 1.25 * inch - 0.875 * inch
    return (2 / np.pi) * as_tot * 470e6 * r_bar / h_pier / 1e3


F_REF = 0.1 * fy_eq_kn()


def force(tid, t_end):
    m = sio.loadmat(R.DATA / "Simulink" / f"{R.BY_ID[tid][1]}.mat", squeeze_me=True, struct_as_record=False)["data"]
    raw = m.topLRaw.data
    k = raw[:, -1] <= t_end
    return raw[k, -1] * R.SQL, -R.LAM ** 3 * raw[k, 0] / 1e3


prs.configure_font()
ids = sorted({t for _, rr in ROWS for t, *_ in rr} | {t for t, *_ in SOL} | {"F11"})
D = {tid: R.load(tid) for tid in ids}
F = {tid: force(tid, D[tid]["t"][-1]) for tid in ids}

fig = plt.figure(figsize=(6.0, 7.2), layout="constrained")
s_top, s_bot = fig.subfigures(2, 1, height_ratios=[2.2, 1.0])
gs = s_top.add_gridspec(4, 1, height_ratios=[0.3, 1.0, 0.3, 1.0])
for r_, (title, rows) in enumerate(ROWS):
    la = s_top.add_subplot(gs[2 * r_])
    fa = s_top.add_subplot(gs[2 * r_ + 1], sharex=la)
    for j, (tid, cfg, c, lab) in enumerate(rows):
        d = D[tid]
        on = d["on"][d["on"] > d["t_su"]] * R.SQL
        la.add_patch(plt.Rectangle((0, j - 0.36), d["t"][-1] * R.SQL, 0.72, fc="0.93", ec="none", zorder=0))
        la.vlines(on, j - 0.34, j + 0.34, color=c, lw=0.4 if on.size > 200 else 0.9)
        la.plot(d["t"][-1] * R.SQL, j, marker=CFG[cfg][1], ms=3.8, color=c, mfc=c if d["complete"] else "white",
                mew=0.8, ls="none", clip_on=False)
        t, f = F[tid]
        z = {"F05": 2, "F27": 1, "F22": 3, "F06": 3, "F07": 1, "F08": 2}[tid]
        fa.plot(t, f, color=c, lw=0.5, zorder=z, label=lab)
        print(f"{tid}: |F| max {np.abs(f).max():.1f} kN, earthquake phase (lab < 300 s) {np.abs(f[t < 300]).max():.1f}")
    la.set_yticks(range(len(rows)))
    la.set_yticklabels([lab.split(",")[0] for _, _, _, lab in rows])
    for tl, (_, _, c, _) in zip(la.get_yticklabels(), rows):
        tl.set_color(c)
    la.set_ylim(len(rows) - 0.5, -0.5)
    la.tick_params(axis="y", length=0, labelbottom=False)
    la.tick_params(axis="x", labelbottom=False)
    for s in ("left", "right", "top"):
        la.spines[s].set_visible(False)
    la.set_title(rf"\textbf{{({'ab'[r_]})}} {title}", loc="left")
    for yv in (F_REF, -F_REF):
        fa.axhline(yv, color="0.35", lw=0.6, ls=(0, (4, 2)), zorder=0)
    fa.text(200, F_REF + 3, r"$0.1F_{y,\mathrm{eq}}$", ha="center", va="bottom", color="0.35")
    fa.set_ylim(-YL, YL)
    fa.set_xlim(0, 475)
    fa.set_ylabel(r"$\lambda_L^3 f_p$ (kN)")
    lg = fa.legend(loc="upper right", frameon=True, handlelength=1.2, ncol=1, borderaxespad=0.2, framealpha=1.0)
    lg.get_frame().set(facecolor="white", edgecolor="none")
    if r_ == 1:
        fa.set_xlabel(r"Laboratory time, $\sqrt{\lambda_L}\,t$ (s)")
    else:
        fa.tick_params(labelbottom=False)
    fa.annotate("solitary wave", (347, -78), xytext=(300, -112), ha="right", va="center", color="0.35",
                arrowprops=dict(arrowstyle="-", color="0.5", lw=0.6, shrinkA=1, shrinkB=2))

bc, bd = s_bot.subplots(1, 2)
d = D["F11"]
t, f = F["F11"]
on = d["on"][d["on"] > d["t_su"]] * R.SQL
tp = t[np.argmax(np.abs(f))]
k = (t > tp - 0.6) & (t < tp + 1.4)
bc.plot(t[k] - tp, f[k], color="0.0", lw=0.8)
for x in on[np.abs(on - tp) < 1.5]:
    bc.axvline(x - tp, color="#FF8C00", lw=1.0, zorder=0)
    bc.text(x - tp + 0.03, 0.97, "slowdown", transform=bc.get_xaxis_transform(), ha="left", va="top", color="#B85C00")
bc.set_xlim(-0.6, 1.4)
bc.set_xlabel(r"$\sqrt{\lambda_L}\,(t - t_{\mathrm{peak}})$ (s)")
bc.set_ylabel(r"$\lambda_L^3 f_p$ (kN)")
bc.set_title(r"\textbf{(c)} isolated slowdown (F11)", loc="left")
print(f"F11 transient peak {f[np.argmax(np.abs(f))]:.1f} kN at {tp:.2f} s")
for tid, name, c in SOL:
    t, f = F[tid]
    w = (t > 335) & (t < 365)
    t0 = t[w][np.argmax(np.abs(f[w]))]
    kk = (t > t0 - 1.5) & (t < t0 + 3.0)
    bd.plot(t[kk] - t0, f[kk], color=c, lw=0.8, label=f"{tid}, {name}")
    print(f"{tid} solitary peak {f[w][np.argmax(np.abs(f[w]))]:.1f} kN")
bd.set_xlim(-1.5, 3.0)
bd.sharey(bc)
bc.set_ylim(-100, 120)
bd.tick_params(labelleft=False)
bd.set_xlabel(r"$\sqrt{\lambda_L}\,(t - t_{\mathrm{peak}})$ (s)")
bd.legend(loc="upper right", frameon=False, handlelength=1.2)
bd.set_title(r"\textbf{(d)} solitary wave", loc="left")
for ax in (bc, bd):
    ax.axhline(0, color="0.85", lw=0.4, zorder=0)
fig.savefig(OUT.with_suffix(".png"), dpi=250)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
