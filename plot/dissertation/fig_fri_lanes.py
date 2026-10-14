"""When the slowdowns occurred in the main campaign, and how far the tests fell behind real time (dissertation Ch. 6;
read-only on data).

(a) representative measured interface displacement (F05) on integration time
(b) slowdown episodes after start-up (ticks), one lane per test, on integration time t_int (prototype); lanes grouped
    by model and ordered by configuration; shaded band = D5-95 (5-95 % Arias intensity); dashed line =
    solitary-wave impact (median t_int of the peak load-cell force over the completed tests); lanes end where the
    test ended; numbers at the right = slowdown episodes after start-up; NC = not completed (stopped early)
(c) delay behind real time after start-up, from the rate-transition window starts (count = 1): windows start every
    dt_sim unless a slowdown stretches one, so the delay accumulates without the host-latency jitter of the target
    arrivals; start-up offset removed; prototype-scaled laboratory clock, symlog axis; marker at the end of each test
Configurations, colors, and markers as fig_fri_capacity.py. Uniform 9 pt font; bold panel tags; constrained layout.
Soil-formulation variants are not shown.
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
CFG = {"pspd": ("#C0390B", r"CPU MKR-$\alpha$, ProfileSPD", "s"),
       "mumps": ("#8C510A", r"CPU MKR-$\alpha$, MUMPS", "D"),
       "cpucudss": ("#2166AC", r"CPU MKR-$\alpha$, cuDSS solves", "^"),
       "gpu": ("#000000", r"GPU MKR-$\alpha$, cuDSS", "o")}
from lab_paths import MAT_EXTRACT_DIR  # noqa: E402


def window_delay(tid, t_su):
    """Delay behind real time from the window starts (model s), start-up offset removed."""
    st = np.load(MAT_EXTRACT_DIR / f"{R.BY_ID[tid][1]}.npz", allow_pickle=True)["stateOS_data"]
    t, cnt = st[:, -1], np.rint(st[:, 1]).astype(int)
    win = t[1:][(cnt[1:] == 1) & (cnt[:-1] != 1)]
    lag = (win - win[0]) - R.DT_SIM * np.arange(win.size)
    k = win >= t_su
    return win[k], lag[k] - lag[k][0]
GROUPS = [("Baseline", [("F03", "pspd"), ("F04", "pspd"), ("F26", "mumps"), ("F22", "gpu"), ("F05", "gpu"),
                        ("F27", "gpu")]),
          ("Moderate", [("F07", "pspd"), ("F08", "mumps"), ("F06", "gpu")]),
          ("Large", [("F14", "cpucudss"), ("F15", "cpucudss"), ("F09", "gpu"), ("F10", "gpu"), ("F12", "gpu"),
                     ("F11", "gpu")]),
          ("XL", [("F13", "gpu")])]
DELAY = [("F07", "pspd"), ("F08", "mumps"), ("F22", "gpu"), ("F15", "cpucudss"), ("F14", "cpucudss"),
         ("F13", "gpu"), ("F06", "gpu"), ("F05", "gpu")]
DLAB = {"F07": "F07", "F08": "F08", "F22": "F22, serial", "F15": "F15", "F14": "F14", "F13": "F13", "F06": "F06",
        "F05": "F05"}
T_END = 472.0
DY = {"F08": 4, "F22": -7, "F14": 5, "F13": -5}  # vertical label offsets (points) in (b) where end points crowd


def impact_lab(r, t_end):
    m = sio.loadmat(R.DATA / "Simulink" / f"{r[1]}.mat", squeeze_me=True, struct_as_record=False)["data"]
    raw = m.topLRaw.data
    t, F = raw[:, -1], -R.LAM ** 3 * raw[:, 0] / 1e3
    k = (t > 200) & (t < 240) & (t <= t_end)
    return t[k][np.argmax(np.abs(F[k]))] if k.any() else np.nan


prs.configure_font()
ids = [tid for _, tl in GROUPS for tid, _ in tl]
D = {tid: R.load(tid) for tid in ids}
imp = [D[t]["to_int"](impact_lab(R.BY_ID[t], D[t]["t"][-1])) for t in ids if D[t]["complete"]]
t_imp = float(np.nanmedian(imp))
print(f"solitary-wave impact t_int: median {t_imp:.1f} s, range {np.nanmin(imp):.1f}-{np.nanmax(imp):.1f} s")

fig = plt.figure(figsize=(6.0, 7.1), layout="constrained")
sa, sb = fig.subfigures(2, 1, height_ratios=[2.2, 1.0])
hx, ax = sa.subplots(2, 1, sharex=True, gridspec_kw=dict(height_ratios=[0.32, 1.75]))
zr = np.load(MAT_EXTRACT_DIR / f"{R.BY_ID['F05'][1]}.npz", allow_pickle=True)["meaSigOS_data"]
dz = D["F05"]
kz = (zr[:, -1] > dz["tk"][0]) & (zr[:, -1] <= dz["t"][-1])
hx.plot(dz["to_int"](zr[kz, -1]), zr[kz, 0] * R.LAM * 1e3, color="0.25", lw=0.5)
hx.axvspan(R.T_D5, R.T_D95, color="#d9e6f2", zorder=-1, lw=0)
hx.set_ylabel(r"$u$ (mm)")
hx.set_ylim(-230, 230)
hx.set_yticks([-200, 0, 200])
hx.text(0.995, 0.95, "interface displacement, F05", transform=hx.transAxes, ha="right", va="top", color="0.3")
hx.tick_params(labelbottom=False)
for s_ in ("right", "top"):
    hx.spines[s_].set_visible(False)
hx.set_title(r"\textbf{(a)}", loc="left")
y, yt, yl, ycol = 0, [], [], []
for g, (name, tl) in enumerate(GROUPS):
    y0 = y
    for tid, cfg in tl:
        d, r = D[tid], R.BY_ID[tid]
        c = CFG[cfg][0]
        ti_end = d["ti_k"][-1]
        ax.add_patch(plt.Rectangle((0, y - 0.36), ti_end, 0.72, fc="0.93", ec="none", zorder=0))
        on = d["on"][d["on"] > d["t_su"]]
        ax.vlines(d["to_int"](on), y - 0.34, y + 0.34, color=c, lw=0.45 if on.size > 200 else 0.9, zorder=2)
        if on.size == 0:
            ax.text(0.5 * ti_end, y, "no slowdowns", ha="center", va="center", color="0.45", zorder=3)
        ax.plot(ti_end, y, marker=CFG[cfg][2], ms=4.2, color=c, mfc=c if d["complete"] else "white", mew=0.9,
                ls="none", zorder=3)
        if not d["complete"]:
            ax.text(ti_end + 7, y, "NC", ha="left", va="center", color=c)
        ax.text(T_END + 30, y, f"{d['n_slow']:,}", ha="right", va="center", color=c)
        yt.append(y)
        yl.append(f"{tid}, " + ("serial" if r[7] == 1 else f"{r[7]}"))
        ycol.append(c)
        y += 1
    ax.text(-0.215, 0.5 * (y0 + y - 1), name, transform=ax.get_yaxis_transform(), ha="center", va="center",
            rotation=90)
    if g < len(GROUPS) - 1:
        ax.axhline(y - 0.25, color="0.7", lw=0.5, xmin=0, xmax=1)
    y += 0.5
ax.axvspan(R.T_D5, R.T_D95, color="#d9e6f2", zorder=-1, lw=0)
ax.text(0.5 * (R.T_D5 + R.T_D95), -1.0, r"$D_{5\text{--}95}$", ha="center", va="center", color="0.3")
ax.axvline(t_imp, color="0.3", lw=0.7, ls=(0, (4, 2)), zorder=1)
ax.text(t_imp, -1.0, "solitary wave", ha="center", va="center", color="0.3")
ax.text(T_END + 30, -1.0, "episodes", ha="right", va="center", color="0.3")
ax.set_yticks(yt)
ax.set_yticklabels(yl)
for tl_, c in zip(ax.get_yticklabels(), ycol):
    tl_.set_color(c)
ax.set_ylim(y - 0.4, -1.6)
ax.set_xlim(0, T_END + 32)
ax.set_xticks(range(0, 451, 50))
ax.set_xlabel(r"Integration time, $t_{\mathrm{int}}$ (s)")
ax.set_title(r"\textbf{(b)}", loc="left")
ax.tick_params(axis="y", length=0)
for s in ("left", "right", "top"):
    ax.spines[s].set_visible(False)
sa.legend([plt.Line2D([], [], color=CFG[k][0], marker=CFG[k][2], ms=4.6, lw=1.2) for k in CFG], [CFG[k][1] for k in CFG],
          loc="outside upper center", ncol=2, frameon=False, handlelength=1.2, columnspacing=1.5)

bx = sb.subplots()
for tid, cfg in DELAY:
    d = D[tid]
    c = CFG[cfg][0]
    tw, lw_ = window_delay(tid, d["t_su"])
    gr = lw_ * R.SQL
    bx.plot(tw * R.SQL, gr, color=c, lw=1.0)
    xe, ye = tw[-1] * R.SQL, gr[-1]
    bx.plot(xe, ye, marker=CFG[cfg][2], ms=4.2, color=c, mfc=c if d["complete"] else "white", mew=0.9,
            ls="none", zorder=4)
    print(f"{tid}: delay growth {gr[-1]:.3f} s (prototype), smallest step {1e3 * np.diff(gr).min():.2f} ms")
    lab = DLAB[tid] + (", NC" if not d["complete"] else "") + f": {max(ye, 0.0):.2g} s"  # final delay (log axis)
    bx.annotate(lab, (xe, ye), xytext=(6, DY.get(tid, 0)), textcoords="offset points", color=c, va="center")
bx.set_yscale("symlog", linthresh=0.01, linscale=0.4)
bx.set_yticks([0, 0.01, 0.1, 1, 10])
bx.set_yticklabels(["0", "0.01", "0.1", "1", "10"])
bx.set_ylim(-0.003, 12)
bx.set_xlim(0, 520)
bx.set_xlabel(r"Laboratory time, $\sqrt{\lambda_L}\,t$ (s)")
bx.set_ylabel("Delay behind\nreal time (s)")
bx.set_title(r"\textbf{(c)}", loc="left")
fig.savefig(OUT.with_suffix(".png"), dpi=250)
fig.savefig(OUT.with_suffix(".pdf"))
print("wrote", OUT.with_suffix(".png"))
