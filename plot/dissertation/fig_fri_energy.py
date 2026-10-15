"""Energy exchanged at the interface over the whole record, main campaign (dissertation Ch. 6; read-only on data).

E = int lambda_L^3 f_p du from the OpenSees recorders (Elmt101_ctrlDsp: commanded displacement; ServerSetup_daqFrc:
force returned to the integrator), as fig_wed_energy2.py, from the start of the earthquake; tau = t_int - t_EQ.
Negative E: energy delivered to the numerical subassembly. Runs with frequent slowdowns (F07 CPU MKR-alpha +
ProfileSPD and F08 CPU MKR-alpha + MUMPS, moderate model, NC; F22 GPU MKR-alpha serial, baseline, NC; F15 CPU MKR-alpha
+ cuDSS solves, large) and references with 0-1 slowdowns (F03, F05, baseline). Configuration colors and end markers as
fig_fri_capacity.py (open = not completed); direct labels; uniform 9 pt font. The Wednesday serial run (W04) is in
fig_wed_energy2.py.
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
CFG = {"pspd": ("#C0390B", "s"), "mumps": ("#8C510A", "D"), "cpucudss": ("#2166AC", "^"), "gpu": ("#000000", "o")}
CASES = [("F07", "pspd", "F07, 8, NC"), ("F08", "mumps", "F08, 8, NC"), ("F22", "gpu", "F22, serial, NC"),
         ("F15", "cpucudss", "F15, 20"), ("F03", None, None), ("F05", None, None)]
TAU_MAX = 480.0  # whole record (the incomplete runs end near 125 s)
T_SOL = 348.0  # solitary-wave peak, t_int (F03, F05: no slowdowns there)
WIN = 0.3  # s of integration time after each slowdown onset
E_Y = 19.8  # kJ, M_y theta_y of the pier base hinge (same section as Wednesday; fig_wed_energy2.my_theta_y())


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
fig, ax = plt.subplots(figsize=(6.0, 2.6), layout="constrained")
ends = []
for tid, cfg, lab in CASES:
    r = R.BY_ID[tid]
    tau, E, _ = energy(r[2], 0.0)
    c, mk = CFG[cfg] if cfg else ("0.6", None)
    ax.plot(tau, E, color=c, lw=1.0 if cfg else 0.7, zorder=3 if cfg else 2)
    if cfg:
        ax.plot(tau[-1], E[-1], marker=mk, ms=4.6, color=c, mfc="white" if r[8] == "incomplete" else c, mew=0.9,
                ls="none", zorder=4)
        ends.append((tau[-1], E[-1], c, lab))
    else:
        ref_end = (tau[-1], E[-1])
    k125 = tau <= 125
    if tau[-1] > T_SOL + 20:
        e0, e1 = np.interp([T_SOL - 10, T_SOL + 20], tau, E)
        print(f"   change across the solitary wave ({T_SOL - 10:.0f}-{T_SOL + 20:.0f} s): {1e3 * (e1 - e0):+.1f} J")
    print(f"{tid}: E at 125 s {E[k125][-1]:.2f} kJ; at the end ({tau[-1]:.0f} s) {E[-1]:.2f} kJ")
for x, y, c, lab in ends:
    ax.text(x + 3, y + (0.28 if lab.startswith("F15") else 0.0), lab, ha="left", va="center", color=c)
ax.text(ref_end[0] + 6, ref_end[1] - 0.45, "F03, F05\n(0-1 slowdowns)", ha="left", va="center", color="0.45", linespacing=0.95)
ax.axhline(0, color="0.5", lw=0.5, zorder=1)
ax.axvline(T_SOL, color="0.35", lw=0.6, ls=(0, (4, 2)), zorder=1)
ax.text(T_SOL + 4, 0.04, "solitary wave", transform=ax.get_xaxis_transform(), ha="left", va="bottom", color="0.35")
ax.set_xlim(0, 540)
ax.set_xticks(range(0, 451, 50))
ax.spines["bottom"].set_bounds(0, 460)
ax.set_xlabel(r"Time from the start of the earthquake, $t_{\mathrm{int}} - t_{\mathrm{EQ}}$ (s)")
ax.set_ylabel(r"$E$ (kJ)")
sec = ax.secondary_yaxis("right", functions=(lambda e: e / E_Y, lambda r: r * E_Y))
sec.set_ylabel(r"$E/(M_{\mathrm{y}}\theta_{\mathrm{y}})$")
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
