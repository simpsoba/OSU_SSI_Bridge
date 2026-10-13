"""F05 interface force against u, du/dt, d2u/dt2 during the solitary wave (a-c; arrows at half the peak force show the
sense of time) and the free vibration before and after it (d; amplitude at f_1, band 0.3-0.7 Hz, over the 10 s on
each side). Dissertation Ch. 6 draft; read-only on data. Usage: fig_fri_sol_kin.py <out.png> [test]."""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, filtfilt

sys.path.insert(0, r"C:\Users\garaujor\OpenSees_Runs\OSU_SSI_Bridge\plot\dissertation")
import fri_runs as R  # noqa: E402
import ctrl_error_wed as C  # noqa: E402

sys.path.insert(0, str(R.ROOT / "plot"))
import PlotResponseSpectrum as prs  # noqa: E402

TID = sys.argv[2] if len(sys.argv) > 2 else "F05"
stem = R.BY_ID[TID][1]
z = np.load(R.MAT_EXTRACT_DIR / f"{stem}.npz", allow_pickle=True)
st, me = z["stateOS_data"], z["meaSigOS_data"]
tm = me[:, -1]
dt = np.median(np.diff(tm))
b, a = butter(4, 20.0 * 2 * dt, "low")
uf = filtfilt(b, a, me[:, 0])
vel = np.gradient(uf, dt)
acc = np.gradient(vel, dt)
tf, fr = R.load_force(stem)
F = -R.LAM ** 3 * np.interp(tm, tf, fr) / 1e3
tl = tm * R.SQL
u, v = me[:, 0] * R.LAM * 1e3, vel * R.SQL * 1e3
w = (tl > 300) & (tl < 380)
t_imp = tl[w][np.argmin(F[w])]
k = (tl >= t_imp - 5.0) & (tl <= t_imp + 30.0)
tau = tl[k] - t_imp
M_EQ = 3.2  # t, fitted during the earthquake (fig_fri_force_kin)
print(f"{TID}: peak {F[k].min():.1f} kN at {t_imp:.2f} s; u range {u[k].min():.1f}..{u[k].max():.1f} mm; "
      f"|v| max {np.abs(v[k]).max():.1f} mm/s; |a| max {np.abs(acc[k]).max():.2f} m/s2; "
      f"m a at the largest |a|: {M_EQ * np.abs(acc[k]).max():.1f} kN")

prs.configure_font()
fig = plt.figure(figsize=(6.0, 4.4), layout="constrained")
gs = fig.add_gridspec(3, 3, height_ratios=[1.0, 0.58, 0.42])
A = [fig.add_subplot(gs[0, 0])]
A += [fig.add_subplot(gs[0, j], sharey=A[0]) for j in (1, 2)]
for ax in A[1:]:
    ax.tick_params(labelleft=False)
D = fig.add_subplot(gs[1, :])
E = fig.add_subplot(gs[2, :], sharex=D)
for ax, x, xl in zip(A, (u[k], v[k], acc[k]), (r"$u$ (mm)", r"$\dot u$ (mm/s)", r"$\ddot u$ (m/s$^2$)")):
    ax.plot(x, F[k], color="#0072B2", lw=0.5)
    ip = np.argmin(F[k])
    ax.plot(x[ip], F[k][ip], "o", ms=3.5, color="#0072B2")
    ax.axhline(0, color="0.85", lw=0.4, zorder=0)
    if xl != r"$u$ (mm)":
        ax.axvline(0, color="0.85", lw=0.4, zorder=0)
    ax.set_xlabel(xl)
# arrows on the branch before and after the peak, at half the peak force, pointing the way time advances
Fk = F[k]
ip = np.argmin(Fk)
half = 0.5 * Fk[ip]
ib = np.flatnonzero(Fk[:ip] > half)[-1]  # last sample above half the peak before it
ia = ip + np.flatnonzero(Fk[ip:] > half)[0]  # first sample above half the peak after it
fig.canvas.draw()
for ax, x in zip(A, (u[k], v[k], acc[k])):
    for i in (ib, ia):
        p0 = ax.transData.transform((x[i - 3], Fk[i - 3]))
        p1 = ax.transData.transform((x[i + 3], Fk[i + 3]))
        d = (p1 - p0) / np.hypot(*(p1 - p0))
        c = ax.transData.transform((x[i], Fk[i]))
        # beside the curve, on the side away from the loop, parallel to it
        cen = ax.transData.transform(np.c_[x[ib:ia + 1], Fk[ib:ia + 1]]).mean(axis=0)
        n = np.array([-d[1], d[0]])
        n = n if np.dot(n, c - cen) > 0 else -n
        bb = ax.get_window_extent()
        if not (bb.x0 + 12 < (c + 9 * n)[0] < bb.x1 - 12 and bb.y0 + 12 < (c + 9 * n)[1] < bb.y1 - 12):
            n = -n  # no room outside: use the inner side
        c = c + 9 * n
        tail, head = ax.transData.inverted().transform(c - 8 * d), ax.transData.inverted().transform(c + 8 * d)
        ax.annotate("", xy=head, xytext=tail, arrowprops=dict(arrowstyle="-|>", color="#000000", lw=1.0,
                                                              mutation_scale=9, shrinkA=0, shrinkB=0))
A[0].set_ylabel(r"$\lambda_L^3 f_p$ (kN)")
for ax, s in zip(A, ("(a)", "(b)", "(c)")):
    ax.set_title(rf"\textbf{{{s}}}" + (f" {TID}, solitary wave" if s == "(a)" else ""), loc="left")
# (d), (e) the free vibration before and after the wave: displacement and force histories and the amplitude at f_1 (band 0.3-0.7 Hz)
from scipy.signal import sosfiltfilt
kd = (tl >= t_imp - 40) & (tl <= t_imp + 70)
D.plot(tl[kd] - t_imp, u[kd], color="#0072B2", lw=0.5, zorder=3)
sos = butter(4, [0.3, 0.7], btype="band", fs=1 / (dt * R.SQL), output="sos")
ub = sosfiltfilt(sos, u)
amp = {}
for nm, (w0, w1) in (("before", (-12, -2)), ("after", (2, 12))):
    kw = (tl - t_imp >= w0) & (tl - t_imp <= w1)
    am, mu = np.abs(ub[kw]).max(), u[kw].mean()
    amp[nm] = am
    for yv in (mu - am, mu + am):
        D.plot([w0, w1], [yv, yv], color="#000000", lw=0.9, ls=(0, (4, 1.5)))
    D.text(0.5 * (w0 + w1), mu + am + 0.3, rf"$\pm${am:.1f} mm", ha="center", va="bottom", color="#000000")
print(f"amplitude at f1: before {amp['before']:.2f} mm, after {amp['after']:.2f} mm, ratio {amp['after'] / amp['before']:.1f}")
D.set_xlim(-40, 70)
D.set_ylim(D.get_ylim()[0], D.get_ylim()[1] + 1.0)
D.tick_params(labelbottom=False)
E.plot(tl[kd] - t_imp, F[kd], color="0.3", lw=0.6)
E.axhline(0, color="0.85", lw=0.4, zorder=0)
E.set_xlabel(r"Time from the peak force, $\sqrt{\lambda_L}\,\tau$ (s)")
E.set_ylabel(r"$\lambda_L^3 f_p$ (kN)")
E.set_title(r"\textbf{(e)} interface force", loc="left")
D.set_ylabel(r"$u$ (mm)")
D.set_title(r"\textbf{(d)} free vibration before and after the wave", loc="left")
fig.savefig(sys.argv[1], dpi=200)
print("wrote", sys.argv[1])
