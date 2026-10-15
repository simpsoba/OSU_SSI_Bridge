"""Interface force in the main campaign: force spikes from computational slowdowns, compared with the solitary-wave
load (dissertation Ch. 6; read-only on data).

lambda_L^3 f_p = -lambda^3 * topLRaw (load cell, resisting-force convention, prototype scale) on the prototype-scaled
laboratory clock, as fig_wed_force.py; dashed lines at +-0.1 F_y,eq (47.7 kN).
(a)-(d) one row per model, with slowdown lanes above the force: (a) baseline, F05 (GPU, 8) and F22 (GPU, serial, NC);
    (b) moderate, F06 (GPU, 8), F08 (CPU + MUMPS, NC), F07 (CPU + ProfileSPD, NC); (c) large, F11 (GPU, 24) and F15
    (CPU MKR-alpha + cuDSS solves, 20); (d) extra-large, F13 (GPU, 24)
(e) the isolated slowdown of F11 at the strongest pulse against F12 (same model, no slowdown), each aligned at the
    arrival of the same integration step (tau = 0 at the slowdown onset in F11)
(f) solitary-wave pulse of the completed GPU tests on the four models, aligned at each peak (tau = 0)
No per-row legends: lane labels = test and ranks, configuration from the marker key (open = NC), as
fig_fri_capacity.py; uniform 9 pt font; bold panel tags; constrained layout.
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
ROWS = [("Baseline, 2,028 DOFs", [("F05", "gpu", "0.55", r"F05, GPU, 8"),
                                  ("F22", "gpu", "#000000", r"F22, GPU, serial, NC")]),
        ("Moderate, 3,216 DOFs", [("F06", "gpu", "#000000", r"F06, GPU, 8"),
                                  ("F08", "mumps", "#8C510A", r"F08, CPU + MUMPS, 8, NC"),
                                  ("F07", "pspd", "#C0390B", r"F07, CPU + ProfileSPD, 8, NC")]),
        ("Large, 4,620 DOFs", [("F11", "gpu", "#000000", r"F11, GPU, 24"),
                               ("F15", "cpucudss", "#2166AC", r"F15, CPU + cuDSS solves, 20")]),
        ("Extra-large, 5,484 DOFs", [("F13", "gpu", "#000000", r"F13, GPU, 24")])]
ZORD = {"F05": 1, "F22": 3, "F06": 3, "F07": 1, "F08": 2, "F11": 3, "F15": 2, "F13": 3}
SOL = [("F05", "baseline", "0.0"), ("F06", "moderate", "0.25"), ("F09", "large", "0.45"),
       ("F13", "extra-large", "0.65")]
YL = 100.0  # the largest spikes (110-133 kN, quoted in the text) are clipped


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
ids = sorted({t for _, rr in ROWS for t, *_ in rr} | {t for t, *_ in SOL} | {"F11", "F12"})
D = {tid: R.load(tid) for tid in ids}
F = {tid: force(tid, D[tid]["t"][-1]) for tid in ids}

fig = plt.figure(figsize=(6.0, 7.0), layout="constrained")
s_top, s_bot = fig.subfigures(2, 1, height_ratios=[3.6, 1.0])
hr = []
for _, rows in ROWS:
    hr += [0.15 * len(rows) + 0.06, 0.62]
gs = s_top.add_gridspec(2 * len(ROWS), 1, height_ratios=hr)
for r_, (title, rows) in enumerate(ROWS):
    la = s_top.add_subplot(gs[2 * r_])
    fa = s_top.add_subplot(gs[2 * r_ + 1], sharex=la)
    for j, (tid, cfg, c, lab) in enumerate(rows):
        d = D[tid]
        on = d["on"][d["on"] > d["t_su"]] * R.SQL
        la.add_patch(plt.Rectangle((0, j - 0.36), d["t"][-1] * R.SQL, 0.72, fc="0.93", ec="none", zorder=0))
        la.vlines(on, j - 0.34, j + 0.34, color=c, lw=0.4 if on.size > 200 else 0.9)
        if on.size == 0:
            la.text(0.5 * d["t"][-1] * R.SQL, j, "no slowdowns", ha="center", va="center", color="0.45", zorder=3)
        la.plot(d["t"][-1] * R.SQL, j, marker=CFG[cfg][1], ms=3.8, color=c, mfc=c if d["complete"] else "white",
                mew=0.8, ls="none", clip_on=False)
        t, f = F[tid]
        fa.plot(t, f, color=c, lw=0.5, zorder=ZORD[tid], label=lab)
        print(f"{tid}: |F| max {np.abs(f).max():.1f} kN, earthquake phase (lab < 300 s) {np.abs(f[t < 300]).max():.1f}")
    la.set_yticks(range(len(rows)))
    la.set_yticklabels([f"{tid}, {lab.split(', ')[2]}" for tid, _, _, lab in rows])
    for tl, (_, _, c, _) in zip(la.get_yticklabels(), rows):
        tl.set_color(c)
    la.set_ylim(len(rows) - 0.5, -0.5)
    la.tick_params(axis="y", length=0, labelbottom=False)
    la.tick_params(axis="x", labelbottom=False)
    for s in ("left", "right", "top"):
        la.spines[s].set_visible(False)
    la.set_title(rf"\textbf{{({'abcd'[r_]})}} {title}", loc="left")
    for yv in (F_REF, -F_REF):
        fa.axhline(yv, color="0.35", lw=0.6, ls=(0, (4, 2)), zorder=0)
    fa.text(472, -F_REF - 3, r"$-0.1F_{\mathrm{y,eq}}$", ha="right", va="top", color="0.35")
    fa.set_ylim(-YL, YL)
    fa.set_xlim(0, 475)
    fa.set_ylabel(r"$\lambda_{\mathrm{L}}^3 f_{\mathrm{p}}$ (kN)")
    if r_ == len(ROWS) - 1:
        fa.set_xlabel(r"Laboratory time, $\sqrt{\lambda_{\mathrm{L}}}\,t$ (s)")
    else:
        fa.tick_params(labelbottom=False)
    if r_ == 0:
        fa.annotate("solitary wave", (347, -60), xytext=(300, -72), ha="right", va="center", color="0.35",
                    arrowprops=dict(arrowstyle="-", color="0.5", lw=0.6, shrinkA=1, shrinkB=2))

KEY = [("gpu", r"GPU MKR-$\alpha$, cuDSS"), ("pspd", r"CPU MKR-$\alpha$, ProfileSPD"), ("mumps", r"CPU MKR-$\alpha$, MUMPS"),
       ("cpucudss", r"CPU MKR-$\alpha$, cuDSS solves")]
h = [plt.Line2D([], [], marker=CFG[k][1], color=CFG[k][0], ls="none", ms=4.6) for k, _ in KEY]
h.append(plt.Line2D([], [], marker="o", color="0.3", mfc="white", mew=0.8, ls="none", ms=4.6))
s_top.legend(h, [lab for _, lab in KEY] + ["not completed (NC)"], loc="outside upper center", ncol=3, frameon=False,
             handletextpad=0.3, columnspacing=1.2)

bc, bd = s_bot.subplots(1, 2)
d = D["F11"]
t, f = F["F11"]
tp = t[np.argmax(np.abs(f))] / R.SQL  # strongest transient (model-scale laboratory clock)
x0 = d["on"][d["on"] <= tp][-1]  # onset of the slowdown that produced it
ti0 = d["to_int"](x0)  # integration time carried by that step
for tid, c, lw, ls, z in (("F12", "0.5", 1.0, (0, (3, 1.5)), 3), ("F11", "0.0", 0.8, "-", 2)):
    t, f = F[tid]
    xs = x0 if tid == "F11" else np.interp(ti0, D[tid]["ti_k"], D[tid]["tk"])  # same step in the other test
    x = t - xs * R.SQL
    k = (x > -0.5) & (x < 1.0)
    bc.plot(x[k], f[k], color=c, lw=lw, ls=ls, zorder=z)
bc.set_xlim(-0.5, 0.5)
bc.set_xlabel(r"Time from the slowdown onset, $\sqrt{\lambda_{\mathrm{L}}}\,\tau$ (s)")
bc.set_ylabel(r"$\lambda_{\mathrm{L}}^3 f_{\mathrm{p}}$ (kN)")
bc.set_title(r"\textbf{(e)} isolated slowdown, F11 in (c)", loc="left")
bc.text(0.05, 104, "F11, 24 (slowdown)", ha="left", va="center", color="0.0")
bc.text(-0.46, -18, "F12, 20 (no slowdown)", ha="left", va="top", color="0.45")
print(f"F11 slowdown onset t_int {ti0:.2f} s")
for tid, name, c in SOL:
    t, f = F[tid]
    w = (t > 335) & (t < 365)
    t0 = t[w][np.argmax(np.abs(f[w]))]
    kk = (t > t0 - 1.0) & (t < t0 + 1.0)
    bd.plot(t[kk] - t0, f[kk], color=c, lw=0.8, label=f"{tid}, {name}")
    print(f"{tid} solitary peak {f[w][np.argmax(np.abs(f[w]))]:.1f} kN")
bd.set_xlim(-0.5, 0.5)
bd.sharey(bc)
bc.set_ylim(-100, 120)
bd.tick_params(labelleft=False)
bd.set_xlabel(r"Time from the peak force, $\sqrt{\lambda_{\mathrm{L}}}\,\tau$ (s)")
bd.text(0.06, -62, "F05, F06, F09, F13\n(four models)", ha="left", va="center", color="0.35", linespacing=0.95)
bd.set_title(r"\textbf{(f)} solitary wave", loc="left")
for ax in (bc, bd):
    ax.axhline(0, color="0.85", lw=0.4, zorder=0)
fig.savefig(OUT.with_suffix(".png"), dpi=250)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
