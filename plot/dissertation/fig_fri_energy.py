"""Energy injected by the slowdowns, Wednesday and Friday campaigns (dissertation Ch. 6; read-only on data).

E = int lambda_L^3 f_p du from the OpenSees recorders (Elmt101_ctrlDsp: commanded displacement; ServerSetup_daqFrc:
force returned to the integrator), as fig_wed_energy2.py, from the start of the earthquake (t_int = 150 s on
Wednesday, 0 s on Friday); tau = t_int - t_EQ.
(a) E(tau) for runs with frequent slowdowns (W04 CPU serial; F07 CPU ProfileSPD and F08 CPU MUMPS, moderate; F22 GPU
    serial, baseline) and reference runs without them (W02, F03, F05)
(b) energy change within 0.3 s after each slowdown onset (summed over the record) against the number of held
    slowdowns (slowdown episodes longer than one stateOS sample, i.e., past 0.9 dt_sim); line through the origin
    fitted by least squares
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
from lab_paths import MAT_EXTRACT_DIR  # noqa: E402

OUT = Path(sys.argv[1])
W04 = ("W04", "0819_GusBridge_Serial_H0p5_T7p746_reg_Trial02", "0819_GusBridge_Serial_H0p5_T7p746_reg_Trial02", 150.0)
W02 = ("W02", "0819_GusBridge_Parallel8_H0p5_T7p746_reg", "0819_GusBridge_Parallel8_H0p5_T7p746_reg", 150.0)
CASES = [(W04, "#7f2704", "W04 CPU, serial (Wed.)"),
         ((*R.BY_ID["F07"][:3], 0.0), "#C0390B", "F07 CPU ProfileSPD, moderate"),
         ((*R.BY_ID["F08"][:3], 0.0), "#E69F00", "F08 CPU MUMPS, moderate"),
         ((*R.BY_ID["F22"][:3], 0.0), "#CC79A7", "F22 GPU, serial, baseline"),
         (W02, "0.55", "W02, F03, F05 (no slowdowns)"),
         ((*R.BY_ID["F03"][:3], 0.0), "0.55", None),
         ((*R.BY_ID["F05"][:3], 0.0), "0.55", None)]
TAU_MAX = 125.0
WIN = 0.3  # s of integration time after each slowdown onset


def rd(p):
    a = np.genfromtxt(p, invalid_raise=False)
    return a[np.isfinite(a).all(1)]


def energy(dump, t0):
    d = R.DATA / dump
    c, f = rd(d / "Elmt101_ctrlDsp.out"), rd(d / "ServerSetup_daqFrc.out")
    t, u = c[:, 0], c[:, 1]
    F = np.interp(t, f[:, 0], f[:, 1])
    m = (t >= t0) & (t <= t0 + TAU_MAX)
    t, u, F = t[m], u[m], F[m]
    dE = 0.5 * (F[1:] + F[:-1]) * np.diff(u) / 1e3
    return t - t0, np.r_[0.0, np.cumsum(dE)], dE


def slowdowns(stem):
    st = np.load(MAT_EXTRACT_DIR / f"{stem}.npz", allow_pickle=True)["stateOS_data"]
    t, typ, s1 = st[:, -1], np.rint(st[:, 0]).astype(int), np.rint(st[:, 2])
    tk = t[1:][(s1[1:] > 0) & (s1[:-1] <= 0)]
    tk = tk[1 if tk[1] - tk[0] < 3e-3 else 0:]
    on = t[1:][(typ[1:] == 2) & (typ[:-1] != 2)]
    off = t[1:][(typ[1:] != 2) & (typ[:-1] == 2)]
    n = min(on.size, off.size)
    ti = np.interp(on[:n], tk, np.arange(tk.size) * R.DT_SIM * R.SQL)
    return ti[3:], (off[:n] - on[:n])[3:]  # drop start-up episodes


prs.configure_font()
fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(6.0, 2.9), layout="constrained",
                                     gridspec_kw=dict(width_ratios=[1.45, 1, 0.75]))
pts = []
for (tid, stem, dump, t0), c, name in CASES:
    tau, E, dE = energy(dump, t0)
    ax_a.plot(tau, E, color=c, lw=1.0 if name and "no slowdowns" not in name else 0.7, label=name)
    on, dur = slowdowns(stem)
    m = (on >= t0) & (on <= t0 + tau[-1])
    tm = 0.5 * (tau[1:] + tau[:-1]) + t0
    near = np.zeros(tm.size, bool)
    for i in np.searchsorted(tm, on[m]):
        near[max(i - 2, 0):i + int(WIN / (R.DT_SIM * R.SQL))] = True
    held = int(np.sum(dur[m] > 0.7e-3))
    pts.append((tid, c, held, dE[near].sum(), E[-1], tau[-1]))
    print(f"{tid}: {m.sum()} slowdowns, {held} held over tau 0-{tau[-1]:.0f} s; E end {E[-1]:.2f} kJ; "
          f"near slowdowns {dE[near].sum():.2f} kJ; elsewhere {dE[~near].sum():.2f} kJ")
ax_a.axhline(0, color="0.5", lw=0.5)
ax_a.set_xlim(0, TAU_MAX)
ax_a.set_xlabel(r"Time from the earthquake start, $t_{\mathrm{int}} - t_{\mathrm{EQ}}$ (s)")
ax_a.set_ylabel(r"$E$ (kJ)")
ax_a.legend(loc="lower left", frameon=False, fontsize=6.5, handlelength=1.2)
ax_a.set_title(r"\textbf{(a)}", loc="left")
# (b): every run, Friday and Wednesday, over the same window from the start of the earthquake; configuration markers
CFG = {"pspd": ("#C0390B", "s"), "mumps": ("#8C510A", "D"), "cpucudss": ("#2166AC", "^"), "gpu": ("#000000", "o")}
WED = [("W04", "0819_GusBridge_Serial_H0p5_T7p746_reg_Trial02"), ("W06", "0819_GusBridge_Parallel4_H0p5_T7p746_reg_Trial02"),
       ("W02", "0819_GusBridge_Parallel8_H0p5_T7p746_reg"), ("W05", "0819_GusBridge_Parallel16_H0p5_T7p746_reg_Trial02"),
       ("W07", "0819_GusBridge_Parallel24_H0p5_T7p746_reg")]
ALL = [(tid, stem, stem, 150.0, "pspd", True, "Wed.") for tid, stem in WED]
for r in R.RUNS:
    cfg = "gpu" if r[5] == "GPU" else {"ProfileSPD": "pspd", "MUMPS": "mumps", "cuDSS": "cpucudss"}[r[6]]
    ALL.append((r[0], r[1], r[2], 0.0, cfg, r[8] != "incomplete", "Fri."))
pts_b = []
for tid, stem, dump, t0, cfg, done, day in ALL:
    if not (R.DATA / dump / "Elmt101_ctrlDsp.out").exists():
        print(f"  (b) {tid}: no OpenSees dump {dump}, skipped")
        continue
    tau, E, dE = energy(dump, t0)
    on, dur = slowdowns(stem)
    m = (on >= t0) & (on <= t0 + tau[-1])
    tm = 0.5 * (tau[1:] + tau[:-1]) + t0
    near = np.zeros(tm.size, bool)
    for i in np.searchsorted(tm, on[m]):
        near[max(i - 2, 0):i + int(WIN / (R.DT_SIM * R.SQL))] = True
    held = int(np.sum(dur[m] > 0.7e-3))
    pts_b.append((tid, cfg, done, day, held, dE[near].sum()))
    print(f"  (b) {tid} ({day}, {cfg}{'' if done else ', NC'}): {held} held, dE near slowdowns {dE[near].sum():.3f} kJ")
for tid, cfg, done, day, held, en in pts_b:
    c, mk = CFG[cfg]
    ax_b.plot(held, en, marker=mk, ls="none", ms=4.6, color=c, mfc=c if done else "white", mew=0.9,
              alpha=0.55 if day == "Wed." else 1.0)
    if held >= 100:
        ax_b.annotate(tid, (held, en), xytext=(4, 2), textcoords="offset points", color=c)
ins = ax_c  # (c): zoom on the runs with few held slowdowns
for tid, cfg, done, day, held, en in pts_b:
    c, mk = CFG[cfg]
    ins.plot(held, 1e3 * en, marker=mk, ls="none", ms=4.0, color=c, mfc=c if done else "white", mew=0.8,
             alpha=0.55 if day == "Wed." else 1.0)
    if 0 < held < 100 and abs(en) > 0.015:
        ins.annotate(tid, (held, 1e3 * en), xytext=(3, 1), textcoords="offset points", color=c)
ins.axhline(0, color="0.5", lw=0.5)
ins.set_xlim(-1.5, 25)
ins.set_ylim(-170, 110)
ins.set_xlabel("Held slowdowns")
ins.set_ylabel(r"$\Delta E$ (J)")
ins.set_title(r"\textbf{(c)} few held", loc="left")
ins.grid(True)
ax_b.indicate_inset_zoom(ins, edgecolor="0.5", lw=0.5) if False else None
x = np.array([p[4] for p in pts_b], float)
y = np.array([p[5] for p in pts_b])
k = np.sum(x * y) / np.sum(x * x)
kf = np.sum(x[x < 1000] * y[x < 1000]) / np.sum(x[x < 1000] ** 2)
xx = np.linspace(0, x.max() * 1.05, 10)
ax_b.plot(xx, k * xx, color="0.5", lw=0.6, ls=(0, (4, 2)), zorder=0)
ax_b.text(0.97, 0.62, f"{1e3 * k:.1f} J per\nheld slowdown", ha="right", va="top", color="0.35", transform=ax_b.transAxes,
          linespacing=0.95)
print(f"fit (all {x.size} runs): {1e3 * k:.2f} J per held slowdown; without W04: {1e3 * kf:.2f} J")
ax_b.axhline(0, color="0.5", lw=0.5)
ax_b.set_xlabel("Held slowdowns")
ax_b.set_ylabel(r"$\Delta E$ within 0.3 s after slowdowns (kJ)")
ax_b.set_title(r"\textbf{(b)}", loc="left")
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
