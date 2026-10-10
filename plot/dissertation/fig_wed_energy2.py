"""Energy injected by the slowdowns and where it goes, Wednesday runs (dissertation Ch. 6; read-only on data).

(a) hybrid commanded pier-top displacement (OpenSees EE control displacement), all runs, integrator time
(b) offline pier-top displacement, 1 and 8 ranks (no physical subassembly; t_int = t + 150 s)
(c) net energy transferred across the interface, E = int lambda_L^3 f_p du (f_p resisting-force convention,
    u commanded); right axis in units of M_y theta_y of the pier base hinge
(d) pier base zero-length section: moment vs rotation
(e) center-pile section below the cap (first segment, first integration point): moment vs curvature
(f) p-y spring on the center pile at the cap soffit: force vs deformation
(d)-(f) over t_int = 150-243.6 s for every run (the 1-rank record ends where its actuator stops tracking).
Stations and the pile-moment interpolation follow plot/PlotLumpedHoldHysteresisPaper.py.
M_y = F_y,eq H_pier (cantilever estimate, as in plot/PlotActuatorForce.py; the M-theta secant stiffness drops
to 0.8 k_el at 2.9-3.4 MN m), theta_y = M_y / k_el with the model's initial base-hinge stiffness
k_el = E_c I_uncr / L_s,I (structure/PierSection.tcl, L_s,I = 0.74 H_pier; the 150-190 s M-theta fit is 2% higher).
Bold panel tags, no titles or legends (direct labels), uniform font, constrained layout.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]  # OSU_SSI_Bridge
sys.path.insert(0, str(ROOT / "plot"))
import PlotResponseSpectrum as prs  # noqa: E402

OUT = Path(sys.argv[1])
DATA = ROOT / "OSU_SSI_BRIDGE_DATA"
RC = Path(r"C:\Users\garaujor\OpenSees_Runs\OSU_SSI_Bridge_rankcheck\rankcheck")
T0, T1 = 150.0, 243.6  # common window (integrator time)
WIN = (226, 246)
ORANGE, BLACK, GRAY = "#C0390B", "#000000", "#8C8C8C"
RUNS = [("W04", "0819_GusBridge_Serial_H0p5_T7p746_reg_Trial02", 1, ORANGE),
        ("W06", "0819_GusBridge_Parallel4_H0p5_T7p746_reg_Trial02", 4, GRAY),
        ("W02", "0819_GusBridge_Parallel8_H0p5_T7p746_reg", 8, BLACK),
        ("W05", "0819_GusBridge_Parallel16_H0p5_T7p746_reg", 16, GRAY),
        ("W07", "0819_GusBridge_Parallel24_H0p5_T7p746_reg", 24, GRAY)]
ORDER = [1, 3, 4, 2, 0]  # grays, 8 ranks, 1 rank on top
COL = {r[2]: r[3] for r in RUNS}
SPRING_TAG, PILE_BEAM_TAG = 22021, 2120
S_IP1 = 0.5 * (1.0 - np.sqrt(0.6))


def rd(p):
    a = np.genfromtxt(p, invalid_raise=False)
    return a[np.isfinite(a).all(1)]


def tag(ax, s):
    ax.set_title(rf"\textbf{{{s}}}", loc="left")


def label(ax, text, xy, xytext, color):
    ax.annotate(text, xy, xytext=xytext, color=color, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=color, lw=0.6, shrinkA=1, shrinkB=1))


def first(d, name):
    return [p for p in (d / name, d / f"{name}.0") if p.exists()][0]


def ranked(d, eles, out, tagno):
    """(column index, recorder file) of element ``tagno`` in a per-rank (or serial) recorder set."""
    for e in sorted(d.glob(f"{eles}*")):
        tags = [int(float(l.split()[0])) for l in e.read_text().splitlines() if l.strip() and not l.startswith("#")]
        if tagno in tags:
            suffix = e.name[len(eles):]
            return tags.index(tagno), d / f"{out}{suffix}"
    raise FileNotFoundError(tagno)


def hybrid(dump):
    d = DATA / dump
    c, f = rd(d / "Elmt101_ctrlDsp.out"), rd(d / "ServerSetup_daqFrc.out")
    t, u = c[:, 0], c[:, 1]
    F = np.interp(t, f[:, 0], f[:, 1])
    m = (t >= T0) & (t <= T1)
    t, u, F = t[m], u[m], F[m]
    E = np.r_[0.0, np.cumsum(0.5 * (F[1:] + F[:-1]) * np.diff(u))] / 1e3
    out = dict(t=t, u=u * 1e3, E=E)
    hf, hd = rd(first(d, "pier_hinge_force.out")), rd(first(d, "pier_hinge_defo.out"))
    n = min(len(hf), len(hd))
    out["hinge"] = (hd[:n, 0], 1e3 * hd[:n, 2], 1e-6 * hf[:n, 2])
    k, ff = ranked(d, "pile_springs_eles.txt", "pile_springs_force.out", SPRING_TAG)
    _, fd = ranked(d, "pile_springs_eles.txt", "pile_springs_defo.out", SPRING_TAG)
    sf, sd = rd(ff), rd(fd)
    n = min(len(sf), len(sd))
    out["spring"] = (sd[:n, 0], 1e3 * sd[:n, 1 + 2 * k], 1e-3 * sf[:n, 1 + 2 * k])
    k, fg = ranked(d, "pile_beam_eles.txt", "pile_beam_globalForce.out", PILE_BEAM_TAG)
    _, fs = ranked(d, "pile_beam_eles.txt", "pile_beam_sec1_defo.out", PILE_BEAM_TAG)
    g, s = rd(fg), rd(fs)
    n = min(len(g), len(s))
    m_ip = -g[:n, 1 + 6 * k + 2] * (1 - S_IP1) + g[:n, 1 + 6 * k + 5] * S_IP1
    out["pile"] = (s[:n, 0], 1e3 * s[:n, 2 + 2 * k], 1e-6 * m_ip)
    return out


def offline(n):
    best = None
    for f in sorted((RC / f"out_np{n}").glob("pier_top_disp.out*")):
        a = rd(f)
        if a.ndim == 2 and (best is None or a.shape[0] > best.shape[0]):
            best = a
    return best[:, 0] + 150.0, best[:, 1] * 1e3


def sorted_xy(series):
    """(x, y) of a loop sorted by x, for placing a label on the curve."""
    t, x, y = series
    m = (t >= T0) & (t <= T1)
    o = np.argsort(x[m])
    return x[m][o], y[m][o]


def my_theta_y():
    """M_y theta_y (kJ) of the pier base hinge: M_y = F_y,eq H_pier, theta_y = M_y / k_el, k_el = E_c I_uncr / L_s,I."""
    inch = 0.0254
    h_pier = 3.02 * (48.0 / 20.0)
    as_tot = 28 * (np.pi / 4) * (1.25 * inch) ** 2
    r_bar = 24 * inch - 2 * inch - 0.625 * inch - 0.875 * inch
    m_y = (2 / np.pi) * as_tot * 470e6 * r_bar  # N m
    e_c, e_s, d_pier = 4700.0 * np.sqrt(28.0) * 1e6, 200e9, 48 * inch
    i_uncr = np.pi / 64 * d_pier ** 4 + (e_s / e_c - 1) * 0.5 * as_tot * r_bar ** 2  # m^4, transformed
    k_el = e_c * i_uncr / (0.74 * h_pier)  # N m / rad
    return m_y, m_y / k_el, m_y * (m_y / k_el) / 1e3


prs.configure_font()
H = {r[2]: hybrid(r[1]) for r in RUNS}
O = {n: offline(n) for n in (1, 8)}
M_Y, TH_Y, E_Y = my_theta_y()
print(f"M_y = {M_Y / 1e6:.2f} MN m, theta_y = {TH_Y * 1e3:.2f} mrad, M_y theta_y = {E_Y:.1f} kJ")

fig, A = plt.subplots(2, 3, figsize=(6.0, 4.5), layout="constrained")
(ax_a, ax_b, ax_c), (ax_d, ax_e, ax_f) = A

for k in ORDER:  # (a) hybrid
    n = RUNS[k][2]
    d = H[n]
    m = (d["t"] >= WIN[0]) & (d["t"] <= WIN[1])
    ax_a.plot(d["t"][m], d["u"][m], color=COL[n], lw=1.1 if n in (1, 8) else 0.8)
ax_a.plot(H[1]["t"][-1], H[1]["u"][-1], marker="x", color=ORANGE, ms=5, mew=1.0, zorder=5)
for n, ls in ((8, "-"), (1, (0, (4, 2)))):  # (b) offline
    t, u = O[n]
    m = (t >= WIN[0]) & (t <= WIN[1])
    ax_b.plot(t[m], u[m], color=COL[n], ls=ls, lw=1.1)
for ax in (ax_a, ax_b):
    ax.set_xlim(*WIN)
    ax.set_ylim(-200, 200)
    ax.set_xlabel(r"$t_{\mathrm{int}}$ (s)")
ax_a.set_ylabel(r"$u$ (mm)")
ax_b.tick_params(labelleft=False)
tag(ax_a, "(a)")
tag(ax_b, "(b)")
label(ax_a, "1 rank", (239.3, 99), (227.0, 170), ORANGE)
label(ax_a, "4, 8, 16, 24", (230.2, -52), (226.8, -170), "0.25")
label(ax_b, "1 and 8 ranks", (230.2, -52), (226.8, -170), "0.25")

for k in ORDER:  # (c) interface energy
    n = RUNS[k][2]
    ax_c.plot(H[n]["t"], H[n]["E"], color=COL[n], lw=1.1 if n in (1, 8) else 0.8)
ax_c.axhline(0, color="0.5", lw=0.6)
ax_c.set_xlim(T0, 250)
ax_c.set_ylim(-9, 3)
ax_c.set_xlabel(r"$t_{\mathrm{int}}$ (s)")
ax_c.set_ylabel(r"$E$ (kJ)")
sec = ax_c.secondary_yaxis("right", functions=(lambda e: e / E_Y, lambda r: r * E_Y))
sec.set_ylabel(r"$E/(M_y\theta_y)$")
tag(ax_c, "(c)")
label(ax_c, "1 rank", (238.0, np.interp(238.0, H[1]["t"], H[1]["E"])), (195.0, -6.5), ORANGE)
label(ax_c, "4, 8, 16, 24", (232.0, np.interp(232.0, H[8]["t"], H[8]["E"])), (160.0, 2.0), "0.25")

for ax, key, xl, yl, tg in ((ax_d, "hinge", r"$\theta$ (mrad)", r"$M$ (MN$\cdot$m)", "(d)"),
                            (ax_e, "pile", r"$\kappa$ ($10^{-3}$/m)", r"$M$ (MN$\cdot$m)", "(e)"),
                            (ax_f, "spring", r"$y$ (mm)", r"$p$ (kN)", "(f)")):
    for k in [0, 1, 3, 4, 2]:  # loops: 1 rank underneath, so it shows only where it goes beyond the others
        n = RUNS[k][2]
        t, x, y = H[n][key]
        m = (t >= T0) & (t <= T1)
        ax.plot(x[m], y[m], color=COL[n], lw=0.9 if n in (1, 8) else 0.6, alpha=1.0 if n in (1, 8) else 0.8)
    ax.axhline(0, color="0.8", lw=0.4, zorder=0)
    ax.axvline(0, color="0.8", lw=0.4, zorder=0)
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    tag(ax, tg)
ax_d.axhline(M_Y / 1e6, color="0.3", lw=0.6, ls=(0, (5, 3)), zorder=0)
ax_d.axhline(-M_Y / 1e6, color="0.3", lw=0.6, ls=(0, (5, 3)), zorder=0)
ax_d.text(ax_d.get_xlim()[0] + 0.5, M_Y / 1e6, r"$M_y$", ha="left", va="bottom", color="0.3")
th1, M1 = sorted_xy(H[1]["hinge"])
label(ax_d, "1 rank", (th1[-1], M1[-1]), (2.0, 5.0), ORANGE)
k1, Mp1 = sorted_xy(H[1]["pile"])
label(ax_e, "1 rank", (k1[-1], Mp1[-1]), (4.0, -2.0), ORANGE)
label(ax_e, "4, 8, 16, 24", (-5.0, np.interp(-5.0, *sorted_xy(H[8]["pile"]))), (-7.5, 2.0), "0.25")

fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
for n, d in H.items():
    t, th, M = d["hinge"]
    print(f"{n:>2} ranks: E(243.6) = {d['E'][-1]:6.2f} kJ = {d['E'][-1] / E_Y:+.2f} M_y theta_y; "
          f"max |theta| {np.abs(th[(t >= T0) & (t <= T1)]).max():.1f} mrad")
print("wrote", OUT.with_suffix(".png"))
