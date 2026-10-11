"""When the slowdowns occurred in the Friday main campaign (dissertation Ch. 6; read-only on data).

(a) slowdown episodes after start-up (ticks), one lane per test, on integration time t_int (prototype); lanes grouped
    by model; shaded band = 5-95 % Arias intensity of the ground motion; dashed line = solitary-wave impact (peak
    load-cell force of each completed test, mapped to t_int); lanes end where the run ended
(b) growth of the delay of the targets behind real time after start-up (prototype-scaled laboratory clock),
    for selected tests; the fixed start-up offset (lag at the fourth target) is removed
Colors by configuration (as fig_fri_host.py).
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
CFG = {("CPU", "ProfileSPD"): "#C0390B", ("CPU", "MUMPS"): "#E69F00", ("CPU", "cuDSS"): "#0072B2",
       ("GPU", "cuDSS"): "#000000"}
GROUPS = [("Baseline", ["F03", "F04", "F26", "F22", "F05", "F27"]),
          ("Moderate", ["F07", "F08", "F06"]),
          ("Large", ["F14", "F15", "F09", "F10", "F12", "F11"]),
          ("XL", ["F13"]),
          ("Soil", ["F24", "F23", "F20"])]
FORM = {"F24": "elastic", "F23": "stiff", "F20": "Quad PIMY"}
DELAY = ["F07", "F08", "F22", "F15", "F14", "F13", "F09", "F06", "F05"]


def col(r):
    return CFG[(r[5], r[6].split()[0])]


def impact_lab(r, t_end):
    m = sio.loadmat(R.DATA / "Simulink" / f"{r[1]}.mat", squeeze_me=True, struct_as_record=False)["data"]
    raw = m.topLRaw.data
    t, F = raw[:, -1], -R.LAM ** 3 * raw[:, 0] / 1e3
    k = (t > 200) & (t < 240) & (t <= t_end)
    return t[k][np.argmax(np.abs(F[k]))] if k.any() else np.nan


prs.configure_font()
ids = [tid for _, tl in GROUPS for tid in tl]
D = {tid: R.load(tid) for tid in ids}
imp = [D[t]["to_int"](impact_lab(R.BY_ID[t], D[t]["t"][-1])) for t in ids if D[t]["complete"]]
t_imp = float(np.nanmedian(imp))
print(f"solitary-wave impact t_int: median {t_imp:.1f} s, range {np.nanmin(imp):.1f}-{np.nanmax(imp):.1f} s")

fig = plt.figure(figsize=(6.0, 6.4), layout="constrained")
sa, sb = fig.subfigures(2, 1, height_ratios=[2.0, 1.0])
ax = sa.subplots()
y, yt, yl = 0, [], []
for g, (name, tl) in enumerate(GROUPS):
    y0 = y
    for tid in tl:
        d, r = D[tid], R.BY_ID[tid]
        c = col(r)
        ti_end = d["ti_k"][-1]
        ax.add_patch(plt.Rectangle((0, y - 0.38), ti_end, 0.76, fc="0.92", ec="none", zorder=0))
        on = d["on"][d["on"] > d["t_su"]]
        ax.vlines(d["to_int"](on), y - 0.36, y + 0.36, color=c, lw=0.5 if on.size > 200 else 0.9, zorder=2)
        if not d["complete"]:
            ax.text(ti_end + 4, y, "stopped", ha="left", va="center", fontsize=6.5, color=c)
        ax.text(470, y, f"{d['n_slow']}", ha="right", va="center", fontsize=6.5, color="#B85C00")
        yt.append(y)
        yl.append(f"{tid}  {FORM.get(tid, r[6].split()[0])}, {r[7]}")
        y += 1
    ax.text(-82, 0.5 * (y0 + y - 1), name, ha="right", va="center", rotation=90, fontsize=7.5)
    if g < len(GROUPS) - 1:
        ax.axhline(y - 0.5 + 0.25, color="0.6", lw=0.5)
    y += 0.5
ax.axvspan(R.T_D5, R.T_D95, color="#d9e6f2", zorder=-1, lw=0)
ax.text(0.5 * (R.T_D5 + R.T_D95), -1.0, r"5--95\% Arias", ha="center", va="center", fontsize=7, color="0.3")
ax.axvline(t_imp, color="0.3", lw=0.7, ls=(0, (4, 2)), zorder=1)
ax.text(t_imp, -1.0, "solitary wave", ha="center", va="center", fontsize=7, color="0.3")
ax.text(470, -1.0, "episodes", ha="right", va="center", fontsize=6.5, color="#B85C00")
ax.set_yticks(yt)
ax.set_yticklabels(yl, fontsize=6.5)
for tl_, tid in zip(ax.get_yticklabels(), ids):
    tl_.set_color(col(R.BY_ID[tid]))
ax.set_ylim(y - 0.4, -1.6)
ax.set_xlim(0, 472)
ax.set_xlabel(r"Integration time, $t_{\mathrm{int}}$ (s)")
ax.tick_params(axis="y", length=0)
for s in ("left", "right", "top"):
    ax.spines[s].set_visible(False)
ax.set_title(r"\textbf{(a)}", loc="left")
NAMES = {("CPU", "ProfileSPD"): r"CPU MKR-$\alpha$ + ProfileSPD", ("CPU", "MUMPS"): r"CPU MKR-$\alpha$ + MUMPS",
         ("CPU", "cuDSS"): r"CPU MKR-$\alpha$ + cuDSS solves", ("GPU", "cuDSS"): r"GPU MKR-$\alpha$ + cuDSS"}
sa.legend([plt.Line2D([], [], color=c, lw=1.5) for c in CFG.values()], [NAMES[k] for k in CFG],
          loc="outside upper center", ncol=2, frameon=False, handlelength=1.2, columnspacing=1.2)

bx = sb.subplots()
for tid in DELAY:
    d, r = D[tid], R.BY_ID[tid]
    g = (d["lag"] - d["lag"][3]) * R.SQL
    bx.plot(d["tk"][3:] * R.SQL, g[3:], color=col(r), lw=1.0)
    xe, ye = d["tk"][-1] * R.SQL, g[-1]
    print(f"{tid}: start-up offset {d['lag'][3] * R.SQL:.3f} s, delay growth {g[-1]:.3f} s (prototype)")
    bx.annotate(f"{tid}", (xe, ye), xytext=(3, 0), textcoords="offset points", fontsize=6.5, color=col(r),
                va="center")
bx.set_yscale("symlog", linthresh=0.01, linscale=0.4)
bx.set_yticks([0, 0.01, 0.1, 1, 10])
bx.set_yticklabels(["0", "0.01", "0.1", "1", "10"])
bx.set_ylim(0, 12)
bx.set_xlim(0, 480)
bx.set_xlabel(r"Laboratory time, $\sqrt{\lambda_L}\,t$ (s)")
bx.set_ylabel("Delay growth after\nstart-up (s)")
bx.set_title(r"\textbf{(b)}", loc="left")
fig.savefig(OUT.with_suffix(".png"), dpi=200)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
