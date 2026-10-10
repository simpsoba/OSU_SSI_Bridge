"""Wave--structure interaction check, Wednesday runs (dissertation Ch. 6; read-only on data).

(a) pier-top displacement on integration time: hybrid Test W02 (8 ranks; commanded displacement of the
    experimental element, Elmt101_ctrlDsp) and the offline analysis with 8 ranks (no physical subassembly, so no
    waves; rank-check run, ground motion at t = 0, shifted to t_int = t + 150 s); NRMSE (offline as reference)
    over the earthquake, 150 s to the Arias 95% time (289.1 s)
(b) measured interface force of W02 under waves (wave period from successive crests of the force,
    low-passed at 0.5 Hz), before and at the start of the earthquake (115-175 s), mapped from the laboratory clock to
    integration time through the new-target flags (W02 runs about 0.06 s behind; invisible at this scale)
(c) commanded pier-top displacement of W02 over the same window, with the offline analysis (at rest before 150 s)
A wave-harmonic least-squares fit of the force was tried and dropped: on the offline displacement (no waves) it
reports harmonic amplitudes comparable to the remainder during strong shaking, i.e. the earthquake response
leaks into it with two-period windows.
One row; bold panel tags, no titles or legends (direct labels), uniform font, constrained layout.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy.signal import butter, filtfilt, find_peaks

ROOT = Path(__file__).resolve().parents[2]  # OSU_SSI_Bridge
sys.path.insert(0, str(ROOT / "plot"))
import PlotResponseSpectrum as prs  # noqa: E402
from lab_paths import CYLINDER_LENGTH_SCALE as LAM, TIME_SCALE_FROUDE as SQL  # noqa: E402

OUT = Path(sys.argv[1])
DATA = ROOT / "OSU_SSI_BRIDGE_DATA"
MATD = DATA / "Simulink"
RC = Path(r"C:\Users\garaujor\OpenSees_Runs\OSU_SSI_Bridge_rankcheck\rankcheck")
RUN = "0819_GusBridge_Parallel8_H0p5_T7p746_reg"
DT_SIM = 10 / 2048
T_EQ0, T_EQ95 = 150.0, 150.0 + 139.05
ZOOM = (115.0, 175.0)
BLACK = "#000000"


def rd(p):
    a = np.genfromtxt(p, invalid_raise=False)
    return a[np.isfinite(a).all(1)]


def tag(ax, s):
    ax.set_title(rf"\textbf{{{s}}}", loc="left")


def label(ax, text, xy, xytext, color, **kw):
    ax.annotate(text, xy, xytext=xytext, color=color, ha=kw.pop("ha", "left"), va=kw.pop("va", "center"),
                arrowprops=dict(arrowstyle="-", color=color, lw=0.6, shrinkA=1, shrinkB=1), **kw)


# hybrid and offline pier top (integration time, prototype)
c = rd(DATA / RUN / "Elmt101_ctrlDsp.out")
th, uh = c[:, 0], c[:, 1] * 1e3
best = None
for f in sorted((RC / "out_np8").glob("pier_top_disp.out*")):
    a = rd(f)
    if a.ndim == 2 and (best is None or a.shape[0] > best.shape[0]):
        best = a
to, uo = best[:, 0] + 150.0, best[:, 1] * 1e3
m = (th >= T_EQ0) & (th <= T_EQ95)
uoi = np.interp(th[m], to, uo)
nrmse = 100 * np.sqrt(np.mean((uh[m] - uoi) ** 2)) / np.ptp(uoi)
print(f"NRMSE over the earthquake {nrmse:.2f}%, max |diff| {np.abs(uh[m] - uoi).max():.2f} mm")

# measured force, laboratory clock -> integration time through the new-target flags
d = sio.loadmat(MATD / f"{RUN}.mat", squeeze_me=True, struct_as_record=False)["data"]
raw, st = d.topLRaw.data, d.stateOS.data
ts, s1 = st[:, -1], np.rint(st[:, 2])
tk = ts[1:][(s1[1:] > 0) & (s1[:-1] <= 0)]  # arrival of new target k (row k - 1 of ctrlDsp)
tf = (np.interp(raw[:, -1], tk, np.arange(tk.size)) - 1) * DT_SIM * SQL
F = -LAM ** 3 * raw[:, 0] / 1e3

prs.configure_font()
fig = plt.figure(figsize=(6.0, 2.5), layout="constrained")
gs = fig.add_gridspec(2, 2, width_ratios=[1.55, 1])
ax_a = fig.add_subplot(gs[:, 0])
ax_b = fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[1, 1], sharex=ax_b)

# (a) hybrid vs offline
ko = (to >= 140) & (to <= 300)
kh = (th >= 140) & (th <= 300)
ax_a.plot(to[ko], uo[ko], color="0.65", lw=1.4, solid_capstyle="butt", zorder=1)
ax_a.plot(th[kh], uh[kh], color=BLACK, lw=0.7, ls=(0, (3, 1.5)), zorder=2)
ax_a.set_xlim(140, 300)
ax_a.set_ylim(-200, 260)
ax_a.set_xlabel(r"Integrator time, $t_{\mathrm{int}}$ (s)")
ax_a.set_ylabel(r"$u$, pier top (mm)")
label(ax_a, "hybrid, W02", (247.0, np.interp(247.0, th, uh)), (175.0, 215.0), BLACK)
label(ax_a, "offline, no waves", (262.0, np.interp(262.0, to, uo)), (230.0, 238.0), "0.55")
ax_a.text(0.02, 0.04, rf"NRMSE = {nrmse:.2f}\%", transform=ax_a.transAxes, ha="left", va="bottom", color="0.25")
tag(ax_a, "(a)")

# (b) force, (c) displacement before the earthquake
kf = (tf >= ZOOM[0]) & (tf <= ZOOM[1])
ax_b.plot(tf[kf], F[kf], color=BLACK, lw=0.5)
ax_b.axvline(T_EQ0, color="0.35", lw=0.6, ls=(0, (4, 2)), zorder=0)
# wave period: successive crests of the low-passed force (0.5 Hz)
bb, aa = butter(2, 0.5 / (0.5 / np.median(np.diff(tf[kf]))))
Fl = filtfilt(bb, aa, F[kf])
pk, _ = find_peaks(Fl, distance=int(8.0 / np.median(np.diff(tf[kf]))), prominence=0.8)
tc = tf[kf][pk]
tc = tc[tc > 135.0]
print("force crests (t_int, s):", np.round(tc, 2), "spacing:", np.round(np.diff(tc), 2))
n_c = len(tc) - 1
y_ar = 3.4
ax_b.annotate("", (tc[0], y_ar), (tc[-1], y_ar), arrowprops=dict(arrowstyle="<->", color="0.25", lw=0.6, shrinkA=0, shrinkB=0))
for x0 in (tc[0], tc[-1]):
    ax_b.plot([x0, x0], [y_ar - 0.5, y_ar + 0.5], color="0.25", lw=0.6)
lab_T = (rf"{n_c}T_w" if n_c > 1 else "T_w")
ax_b.text(0.5 * (tc[0] + tc[-1]), y_ar + 0.25, rf"${lab_T}$ = {tc[-1] - tc[0]:.1f} s", ha="center", va="bottom",
          color="0.25", bbox=dict(fc="white", ec="none", pad=0.6))
ax_b.set_ylim(-4.5, 5.6)
ax_b.set_ylabel(r"$\lambda_L^3 f_p$ (kN)")
ax_b.tick_params(labelbottom=False)
tag(ax_b, "(b)")
kz = (th >= ZOOM[0]) & (th <= ZOOM[1])
ko2 = (to >= T_EQ0) & (to <= ZOOM[1])
ax_c.plot(np.r_[ZOOM[0], to[ko2]], np.r_[0.0, uo[ko2]], color="0.65", lw=1.4, solid_capstyle="butt", zorder=0)
ax_c.plot(th[kz], uh[kz], color=BLACK, lw=0.7, ls=(0, (3, 1.5)))
ax_c.axvline(T_EQ0, color="0.35", lw=0.6, ls=(0, (4, 2)), zorder=0)
ax_c.set_xlim(*ZOOM)
lim = 1.15 * max(np.abs(uh[kz]).max(), np.abs(uo[ko2]).max())
ax_c.set_ylim(-lim, lim)
ax_c.text(T_EQ0 + 1.0, 0.95 * lim, "earthquake", ha="left", va="top", color="0.35")
print(f"(c) max |u| hybrid {np.abs(uh[kz]).max():.2f} mm, offline {np.abs(uo[ko2]).max():.2f} mm")
ax_c.set_xlabel(r"$t_{\mathrm{int}}$ (s)")
ax_c.set_ylabel(r"$u$ (mm)")
tag(ax_c, "(c)")

fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
