"""Effect of slowdowns on the interface displacement, Friday (dissertation Ch. 6; read-only on data).

Measured interface displacement lambda_L * u_m (prototype mm). Comparisons stay within one integrator
implementation, so that differences are due to the slowdowns (and run-to-run variability) only.
(a) moderate model on the prototype-scaled laboratory clock: F07 (CPU MKR-alpha + ProfileSPD) and F08 (CPU MKR-alpha
    + MUMPS), both stopped early, against F06 (GPU MKR-alpha); slowdown lanes above; (a.1) zoom at the end of F07
(b) baseline model on integration time: F22 (GPU MKR-alpha, serial, 836 slowdowns, stopped) against F05 and F27
    (GPU MKR-alpha, 8 ranks, 0 and 1 slowdowns); (b.1) difference from F05
(c) delay growth of the targets behind real time after start-up (as fig_fri_lanes.py)
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
COL = {"F07": "#C0390B", "F08": "#E69F00", "F06": "#000000", "F22": "#C0390B", "F05": "#000000", "F27": "#8C8C8C"}
NAME = {"F07": "F07 CPU + ProfileSPD, 8 ranks", "F08": "F08 CPU + MUMPS, 8 ranks", "F06": "F06 GPU, 8 ranks",
        "F22": "F22 GPU, serial", "F05": "F05 GPU, 8 ranks", "F27": "F27 GPU, 8 ranks (repeat)"}


def disp(tid):
    d = R.load(tid)
    z = np.load(MAT_EXTRACT_DIR / f"{R.BY_ID[tid][1]}.npz", allow_pickle=True)
    me = z["meaSigOS_data"]
    t = me[:, -1]
    k = (t > d["tk"][0]) & (t <= d["t"][-1])
    return dict(d=d, t=t[k] * R.SQL, ti=d["to_int"](t[k]), u=me[k, 0] * R.LAM * 1e3)


prs.configure_font()
U = {tid: disp(tid) for tid in COL}
fig = plt.figure(figsize=(6.0, 7.0), layout="constrained")
sa, sb, sc = fig.subfigures(3, 1, height_ratios=[1.45, 1.25, 0.8])

# (a) moderate model, laboratory clock
ga = sa.add_gridspec(2, 2, height_ratios=[0.35, 1.6], width_ratios=[2.2, 1.0])
la = sa.add_subplot(ga[0, 0])
aa = sa.add_subplot(ga[1, 0], sharex=la)
az = sa.add_subplot(ga[:, 1])
for j, tid in enumerate(["F07", "F08", "F06"]):
    d = U[tid]["d"]
    on = d["on"][d["on"] > d["t_su"]] * R.SQL
    la.add_patch(plt.Rectangle((0, j - 0.38), d["t"][-1] * R.SQL, 0.76, fc="0.92", ec="none", zorder=0))
    la.vlines(on, j - 0.36, j + 0.36, color=COL[tid], lw=0.4)
for tid in ["F06", "F08", "F07"]:
    u = U[tid]
    k = u["t"] < 140
    aa.plot(u["t"][k], u["u"][k], color=COL[tid], lw=0.8 if tid == "F06" else 0.7, label=NAME[tid])
    kz = (u["t"] > 112) & (u["t"] < 126)
    az.plot(u["t"][kz], u["u"][kz], color=COL[tid], lw=0.9)
la.set_yticks(range(3))
la.set_yticklabels(["F07", "F08", "F06"], fontsize=6.5)
for tl, tid in zip(la.get_yticklabels(), ["F07", "F08", "F06"]):
    tl.set_color(COL[tid])
la.set_ylim(2.5, -0.5)
la.tick_params(axis="y", length=0, labelbottom=False)
la.set_title(r"\textbf{(a)} Moderate model", loc="left")
aa.set_xlim(0, 140)
aa.set_xlabel(r"Laboratory time, $\sqrt{\lambda_L}\,t$ (s)")
aa.set_ylabel(r"$\lambda_L u_m$ (mm)")
aa.legend(loc="lower left", frameon=False, fontsize=6.5, handlelength=1.2)
aa.add_patch(plt.Rectangle((112, -150), 14, 300, fc="none", ec="0.4", lw=0.5, ls="--"))
az.set_xlim(112, 126)
az.set_xlabel(r"$\sqrt{\lambda_L}\,t$ (s)")
az.set_title(r"\textbf{(a.1)}", loc="left")

# (b) baseline model, integration time, GPU only
gb = sb.add_gridspec(2, 1, height_ratios=[1.5, 0.8])
bb = sb.add_subplot(gb[0])
bd = sb.add_subplot(gb[1], sharex=bb)
g = np.arange(5, 130, 0.005)
ref = np.interp(g, U["F05"]["ti"], U["F05"]["u"])
for tid in ["F05", "F27", "F22"]:
    u = U[tid]
    k = u["ti"] < 130
    bb.plot(u["ti"][k], u["u"][k], color=COL[tid], lw=0.8 if tid != "F27" else 0.6, label=NAME[tid],
            ls="-" if tid != "F22" else (0, (3, 1.5)), zorder=3 if tid == "F22" else 2)
    if tid != "F05":
        diff = np.interp(g, u["ti"], u["u"], right=np.nan) - ref
        bd.plot(g, diff, color=COL[tid], lw=0.6)
        print(f"{tid} - F05 over t_int 5-128 s: max |diff| {np.nanmax(np.abs(diff)):.1f} mm, "
              f"RMS {np.sqrt(np.nanmean(diff ** 2)):.2f} mm")
d22 = U["F22"]["d"]
on = d22["on"][d22["on"] > d22["t_su"]]
bd.vlines(d22["to_int"](on), 6.2, 7.6, color="#FF8C00", lw=0.3)
bd.text(129, -6.6, "ticks: F22 slowdowns", ha="right", va="center", fontsize=6.5, color="#B85C00")
bb.set_xlim(5, 130)
bb.set_ylabel(r"$\lambda_L u_m$ (mm)")
bb.legend(loc="lower left", frameon=False, fontsize=6.5, handlelength=1.6)
bb.tick_params(labelbottom=False)
bb.set_title(r"\textbf{(b)} Baseline model, GPU MKR-$\alpha$", loc="left")
bd.set_ylim(-8, 8)
bd.set_ylabel("Difference\nfrom F05 (mm)", fontsize=7)
bd.set_xlabel(r"Integration time, $t_{\mathrm{int}}$ (s)")
bd.set_title(r"\textbf{(b.1)}", loc="left")

# (c) delay growth
cx = sc.subplots()
for tid in ["F07", "F08", "F22", "F06", "F05"]:
    d = U[tid]["d"]
    gr = (d["lag"] - d["lag"][3]) * R.SQL
    cx.plot(d["tk"][3:] * R.SQL, gr[3:], color=COL[tid] if tid != "F22" else "#7f2704", lw=0.9,
            ls=(0, (3, 1.5)) if tid == "F22" else "-")
    cx.annotate(tid, (d["tk"][-1] * R.SQL, gr[-1]), xytext=(3, 0), textcoords="offset points", fontsize=6.5,
                va="center", color=COL[tid] if tid != "F22" else "#7f2704")
    print(f"{tid}: delay growth at end {gr[-1]:.2f} s (prototype), record {d['tk'][-1] * R.SQL:.0f} s")
cx.set_xlim(0, 480)
cx.set_ylim(-0.3, 6.5)
cx.set_xlabel(r"Laboratory time, $\sqrt{\lambda_L}\,t$ (s)")
cx.set_ylabel("Delay growth (s)")
cx.set_title(r"\textbf{(c)}", loc="left")
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
