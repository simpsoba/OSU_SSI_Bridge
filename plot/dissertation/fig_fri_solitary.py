"""Interface response across mesh refinement and under the solitary wave, Friday (dissertation Ch. 6; read-only on
data).

GPU MKR-alpha + cuDSS tests of the four models (same integrator implementation): F05 baseline (8 ranks), F06
moderate (8), F09 large (8), F13 extra-large (24); and F03 (baseline, CPU MKR-alpha + ProfileSPD, 8 ranks).
(a) measured interface displacement on integration time, with (a.1) a zoom on the strong shaking
(b) solitary wave, each test aligned at its peak load-cell force (tau = 0), prototype-scaled laboratory clock:
    (b.1) interface force (zoom on the pulse), (b.2) interface displacement (free vibration before and after the impact)
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy.signal import find_peaks

sys.path.insert(0, str(Path(__file__).parent))
import fri_runs as R  # noqa: E402

sys.path.insert(0, str(R.ROOT / "plot"))
import PlotResponseSpectrum as prs  # noqa: E402
from lab_paths import MAT_EXTRACT_DIR  # noqa: E402

OUT = Path(sys.argv[1])
MESH = [("F05", "baseline, 2028 DOFs", "#000000"), ("F06", "moderate, 3216", "#2b6ca3"),
        ("F09", "large, 4620", "#6aa8d8"), ("F13", "extra-large, 5484", "#E69F00")]
CPU = ("F03", "baseline, CPU MKR-$\\alpha$", "#C0390B")


def load(tid):
    d = R.load(tid)
    r = R.BY_ID[tid]
    z = np.load(MAT_EXTRACT_DIR / f"{r[1]}.npz", allow_pickle=True)
    me = z["meaSigOS_data"]
    t = me[:, -1]
    k = (t > d["tk"][0]) & (t <= d["t"][-1])
    m = sio.loadmat(R.DATA / "Simulink" / f"{r[1]}.mat", squeeze_me=True, struct_as_record=False)["data"]
    raw = m.topLRaw.data
    kf = raw[:, -1] <= d["t"][-1]
    tf, F = raw[kf, -1], -R.LAM ** 3 * raw[kf, 0] / 1e3
    w = (tf > 200) & (tf < 240)
    t_imp = tf[w][np.argmax(np.abs(F[w]))]
    return dict(t=t[k], ti=d["to_int"](t[k]), u=me[k, 0] * R.LAM * 1e3, tf=tf, F=F, t_imp=t_imp)


prs.configure_font()
L = {tid: load(tid) for tid, *_ in MESH + [CPU]}
fig = plt.figure(figsize=(6.0, 5.8), layout="constrained")
sa, sb = fig.subfigures(2, 1, height_ratios=[1.0, 1.15])
aa, az = sa.subplots(1, 2, gridspec_kw=dict(width_ratios=[2.2, 1.0]))
for tid, name, c in MESH:
    x = L[tid]
    aa.plot(x["ti"], x["u"], color=c, lw=0.6, label=f"{tid} {name}")
    k = (x["ti"] > 88) & (x["ti"] < 100)
    az.plot(x["ti"][k], x["u"][k], color=c, lw=0.9)
aa.set_xlim(0, 440)
aa.set_xlabel(r"Integration time, $t_{\mathrm{int}}$ (s)")
aa.set_ylabel(r"$\lambda_{\mathrm{L}} u_{\mathrm{m}}$ (mm)")
aa.legend(loc="upper right", frameon=False, fontsize=6.5, handlelength=1.2)
aa.add_patch(plt.Rectangle((88, -170), 12, 360, fc="none", ec="0.4", lw=0.5, ls="--"))
aa.set_title(r"\textbf{(a)} GPU MKR-$\alpha$, four models", loc="left")
az.set_xlim(88, 100)
az.set_xlabel(r"$t_{\mathrm{int}}$ (s)")
az.set_title(r"\textbf{(a.1)}", loc="left")

b1, b2 = sb.subplots(1, 2, gridspec_kw=dict(width_ratios=[1.0, 2.2]))
for tid, name, c in [CPU] + MESH:
    x = L[tid]
    tau_f = (x["tf"] - x["t_imp"]) * R.SQL
    kf = (tau_f > -1.5) & (tau_f < 3.0)
    b1.plot(tau_f[kf], x["F"][kf], color=c, lw=0.7, label=f"{tid} {name}", zorder=3 if tid == "F03" else 2)
    tau = (x["t"] - x["t_imp"]) * R.SQL
    k = (tau > -40) & (tau < 70)
    b2.plot(tau[k], x["u"][k], color=c, lw=0.7, zorder=3 if tid == "F03" else 2, label=f"{tid} {name}")
    ub = x["u"][(tau > -30) & (tau < -1)]
    ua = x["u"][(tau > 1) & (tau < 30)]
    pk = np.min(x["F"][kf])
    print(f"{tid}: impact at lab {x['t_imp'] * R.SQL:.1f} s (t_int {np.interp(x['t_imp'], x['t'], x['ti']):.1f}); "
          f"peak force {pk:.1f} kN; displacement p-p 30 s before {np.ptp(ub):.1f} mm, 30 s after {np.ptp(ua):.1f} mm, "
          f"mean before {ub.mean():.1f} after {ua.mean():.1f} mm")
b1.set_ylabel(r"$\lambda_{\mathrm{L}}^3 f_{\mathrm{p}}$ (kN)")
b1.set_xlim(-1.5, 3.0)
b1.set_xlabel(r"$\sqrt{\lambda_{\mathrm{L}}}\,\tau$ (s)")
b1.set_title(r"\textbf{(b.1)} Solitary wave", loc="left")
b2.legend(loc="upper right", frameon=False, fontsize=6.5, handlelength=1.2, ncol=2)
b2.set_title(r"\textbf{(b.2)} Aligned at the peak force ($\tau = 0$)", loc="left")
b2.set_ylabel(r"$\lambda_{\mathrm{L}} u_{\mathrm{m}}$ (mm)")
b2.set_ylim(60, 75)
b2.set_xlabel(r"Time from the peak force, $\sqrt{\lambda_{\mathrm{L}}}\,\tau$ (s)")
b2.set_xlim(-40, 70)
for ax in (b1, b2):
    ax.axvline(0, color="0.5", lw=0.5, ls=":")
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
