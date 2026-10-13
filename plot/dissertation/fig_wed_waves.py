"""Wave--structure interaction check and interface force vs specimen motion, Wednesday Test W05 (dissertation Ch. 6,
figure prelim_wave_check; read-only on data). Row 2 reuses load() and phases() of fig_wed_force_kin.py.
W05 (16 ranks) replaces W02 (8 ranks) because it had no slowdowns after start-up (W02 had one, with a 35 kN transient).

Row 1, as fig_wed_waves.py but for W05 (16 ranks, no slowdowns after start-up):
(a) pier-top displacement on integration time: hybrid W05 (commanded displacement, Elmt101_ctrlDsp) and the offline
    analysis with 8 ranks (no physical subassembly, so no waves; t_int = t + 150 s); NRMSE (offline as reference)
    over the earthquake, 150 s to the Arias 95% time (289.1 s)
(b) measured interface force under waves before and at the start of the earthquake (115-175 s), mapped to
    integration time through the new-target flags; wave period from successive crests of the 0.5 Hz low-passed force
(c) commanded pier-top displacement over the same window, with the offline analysis
Row 2, as fig_wed_force_kin.py: interface force against measured u, u-dot, u-ddot (20 Hz zero-phase low-pass on the
displacement, central differences), whole record by phase (before the earthquake: waves only; during: 150 s to 95 %
Arias; after); line in (f): F = m a + F0 fitted during the earthquake.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy.signal import butter, filtfilt, find_peaks

sys.path.insert(0, str(Path(__file__).parent))
import fig_wed_force_kin as K  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "plot"))
import PlotResponseSpectrum as prs  # noqa: E402
from lab_paths import CYLINDER_LENGTH_SCALE as LAM, TIME_SCALE_FROUDE as SQL  # noqa: E402

OUT = Path(sys.argv[1])
DATA = ROOT / "OSU_SSI_BRIDGE_DATA"
RC = Path(r"C:\Users\garaujor\OpenSees_Runs\OSU_SSI_Bridge_rankcheck\rankcheck")
TEST, RUN, DUMP = "W05", "0819_GusBridge_Parallel16_H0p5_T7p746_reg_Trial02", "0819_GusBridge_Parallel16_H0p5_T7p746_reg"
DT_SIM = 10 / 2048
T_EQ0, T_EQ95 = 150.0, 150.0 + 139.05
ZOOM = (115.0, 175.0)
BLACK, RED = "#000000", "#C0390B"
PHASE_COL = {"before (waves only)": "#08519C", "during the earthquake": "0.55", "after": "#6A3D9A"}


def rd(p):
    a = np.genfromtxt(p, invalid_raise=False)
    return a[np.isfinite(a).all(1)]


def tag(ax, s):
    ax.set_title(rf"\textbf{{{s}}}", loc="left")


def label(ax, text, xy, xytext, color, **kw):
    ax.annotate(text, xy, xytext=xytext, color=color, ha=kw.pop("ha", "left"), va=kw.pop("va", "center"),
                arrowprops=dict(arrowstyle="-", color=color, lw=0.6, shrinkA=1, shrinkB=1), **kw)


# ---------------- row 1 data
c = rd(DATA / DUMP / "Elmt101_ctrlDsp.out")
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
print(f"{TEST}: NRMSE over the earthquake {nrmse:.2f}%, max |diff| {np.abs(uh[m] - uoi).max():.2f} mm, "
      f"peak |u| hybrid {np.abs(uh[m]).max():.0f} mm, offline {np.abs(uoi).max():.0f} mm")
d = sio.loadmat(DATA / "Simulink" / f"{RUN}.mat", squeeze_me=True, struct_as_record=False)["data"]
raw, st = d.topLRaw.data, d.stateOS.data
ts, s1 = st[:, -1], np.rint(st[:, 2])
tk = ts[1:][(s1[1:] > 0) & (s1[:-1] <= 0)]
tk = tk[1 if tk[1] - tk[0] < 3e-3 else 0:]
tf = (np.interp(raw[:, -1], tk, np.arange(tk.size)) - 1) * DT_SIM * SQL
F = -LAM ** 3 * raw[:, 0] / 1e3

prs.configure_font()
fig = plt.figure(figsize=(6.0, 5.2), layout="constrained")
s1f, s2f = fig.subfigures(2, 1, height_ratios=[1.0, 1.05])
gs = s1f.add_gridspec(2, 2, width_ratios=[1.55, 1])
ax_a = s1f.add_subplot(gs[:, 0])
ax_b = s1f.add_subplot(gs[0, 1])
ax_c = s1f.add_subplot(gs[1, 1], sharex=ax_b)

ko = (to >= 140) & (to <= 300)
kh = (th >= 140) & (th <= 300)
ax_a.plot(to[ko], uo[ko], color="0.65", lw=1.4, solid_capstyle="butt", zorder=1)
ax_a.plot(th[kh], uh[kh], color=BLACK, lw=0.7, ls=(0, (3, 1.5)), zorder=2)
ax_a.set_xlim(140, 300)
ax_a.set_ylim(-200, 260)
ax_a.set_xlabel(r"Integration time, $t_{\mathrm{int}}$ (s)")
ax_a.set_ylabel(r"$u$, pier top (mm)")
label(ax_a, f"hybrid, {TEST}", (205.0, np.interp(205.0, th, uh)), (160.0, 150.0), BLACK)
label(ax_a, "offline, no waves", (283.0, np.interp(283.0, to, uo)), (218.0, 235.0), "0.55")
ax_a.text(0.02, 0.04, rf"NRMSE = {nrmse:.2f}\%", transform=ax_a.transAxes, ha="left", va="bottom", color="0.25")
tag(ax_a, "(a)")

kf = (tf >= ZOOM[0]) & (tf <= ZOOM[1])
ax_b.plot(tf[kf], F[kf], color=BLACK, lw=0.5)
ax_b.axvline(T_EQ0, color="0.35", lw=0.6, ls=(0, (4, 2)), zorder=0)
bb, aa = butter(2, 0.5 / (0.5 / np.median(np.diff(tf[kf]))))
Fl = filtfilt(bb, aa, F[kf])
pk, _ = find_peaks(Fl, distance=int(8.0 / np.median(np.diff(tf[kf]))), prominence=0.8)
tc = tf[kf][pk]
tc = tc[tc > 140.0]  # the first W05 crest (137.4 s) precedes a developing wave train (10.2 s spacing)
print("force crests (t_int, s):", np.round(tc, 2), "spacing:", np.round(np.diff(tc), 2))
n_c = len(tc) - 1
y_ar = 3.4
ax_b.annotate("", (tc[0], y_ar), (tc[-1], y_ar),
              arrowprops=dict(arrowstyle="<->", color="0.25", lw=0.6, shrinkA=0, shrinkB=0))
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
print(f"(c) max |u| hybrid before EQ {np.abs(uh[(th >= ZOOM[0]) & (th < T_EQ0)]).max():.2f} mm")
ax_c.set_xlabel(r"$t_{\mathrm{int}}$ (s)")
ax_c.set_ylabel(r"$u$ (mm)")
tag(ax_c, "(c)")

# ---------------- row 2: force vs motion
dk = K.load(RUN)
P = [(n, s_, PHASE_COL[n]) for n, s_, _ in K.phases(TEST, dk, 150.0)]
dur = P[0][1]
mm, f0 = np.linalg.lstsq(np.c_[dk["a"][dur], np.ones(dur.sum())], dk["F"][dur], rcond=None)[0]
for n, s_, _ in P:
    print(f"  {n:22}: corr(F, u/v/a) = "
          f"{[round(np.corrcoef(dk[q][s_], dk['F'][s_])[0, 1], 2) for q in ('u', 'v', 'a')]}")
print(f"  m (earthquake) = {mm:.2f} t ({1e3 * mm / LAM ** 3:.0f} kg model)")
A = s2f.subplots(1, 3, sharey=True)
for j, (q, xl, tg) in enumerate((("u", r"$u$ (mm)", "(d)"), ("v", r"$\dot u$ (mm/s)", "(e)"),
                                 ("a", r"$\ddot u$ (m/s$^2$)", "(f)"))):
    ax = A[j]
    for n, s_, col in P:
        wav = n.startswith("before")
        ax.plot(dk[q][s_], dk["F"][s_], color=col, lw=0.4 if wav else 0.15, alpha=0.9 if wav else 0.55,
                rasterized=True)
    ax.axhline(0, color="0.85", lw=0.4, zorder=0)
    ax.axvline(0, color="0.85", lw=0.4, zorder=0)
    ax.set_xlabel(xl)
    ax.set_ylim(-9.5, 11)
    tag(ax, tg)
xx = np.array([dk["a"][dur].min(), dk["a"][dur].max()])
A[2].plot(xx, mm * xx + f0, color=RED, lw=1.4)
A[2].annotate(rf"$m\ddot u$, $m$ = {mm:.1f} t", (0.85 * xx[1], mm * 0.85 * xx[1] + f0), xytext=(0.05, 0.93),
              textcoords="axes fraction", color=RED, ha="left", va="top",
              arrowprops=dict(arrowstyle="-", color=RED, lw=0.6, shrinkA=1, shrinkB=1))
A[0].set_ylabel(r"$\lambda_L^3 f_p$ (kN)")
order = ["before (waves only)", "during the earthquake", "after"]
s2f.legend([plt.Line2D([], [], color=PHASE_COL[n], lw=1.6) for n in order],
           ["before the earthquake (waves only)", "during the earthquake", "after the earthquake"],
           loc="outside upper center", ncol=3, frameon=False, handlelength=1.5, columnspacing=1.2)
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"), dpi=300)
print("wrote", OUT.with_suffix(".png"))
