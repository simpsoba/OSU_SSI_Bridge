"""Solitary-wave impact: laboratory frames and the deformed numerical model (dissertation Ch. 6 draft; read-only on
data and videos). For visualization: the frames and the deformed shapes come from different tests with the same
solitary wave.

(a) wide view, video 2026-08-21-11h24m at 8.611 s (Test F10)
(b) close-up, video 2026-08-21-09h56m at 2:57.0, mirrored so that the wave travels in +x
(e) optional (SHOW_E): pier-base moment-rotation inset between (c) and (d); off by default
Layout 2 x 2: photos on top, model panels below without axes (5 m scale bar), cropped just below the pile tips.
(c) deformed near-field model of Test F30 (r-17, all-node recorders) at the peak load-cell force of the solitary wave,
    displacements amplified AMP; undeformed mesh in light grey
(d) post-earthquake state just before the wave (T_PRE before the peak force, amplified AMP) in light grey, and the
    change induced by the wave, u(t_w) - u(t_pre), added on top with amplification AMP_W, where t_w is the time of the
    largest wave-induced pier-top displacement within WIN_W after the peak force
Frames are extracted with the ffmpeg bundled with imageio-ffmpeg into out/friday/impact_frames (cached).
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from matplotlib.collections import LineCollection, PolyCollection

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import fri_runs as R  # noqa: E402
import PlotEQ as peq  # noqa: E402
import PlotEQParallel as pep  # noqa: E402
import PlotResponseSpectrum as prs  # noqa: E402

OUT = Path(sys.argv[1])
VID = Path(r"G:\Shared drives\Simpson team\Test Data\2026-OSU-SSI-Bridge\Videos")
FRAMES = [("2026-08-21-11h24m.MOV", 8.611, False, "wide"), ("2026-08-21-09h56m.mov", 177.0, True, "close")]
DUMP = R.DATA / "r-17_20260821_1649_Storm_Wave"
TEST = "F30"
AMP, AMP_W = 20.0, 400.0
T_PRE, WIN_W = 2.0, 4.0  # s, prototype integration time
RED, PURPLE = "#C0390B", "#6A3D9A"


def frame(name, t, flip, tag):
    import imageio_ffmpeg
    cache = HERE / "out" / "friday" / "impact_frames"
    cache.mkdir(parents=True, exist_ok=True)
    png = cache / f"{tag}_{t:.3f}{'_flip' if flip else ''}.png"
    if not png.exists():
        cmd = [imageio_ffmpeg.get_ffmpeg_exe(), "-hide_banner", "-loglevel", "error", "-ss", f"{t}", "-i",
               str(VID / name), "-frames:v", "1"] + (["-vf", "hflip"] if flip else []) + ["-y", str(png)]
        subprocess.run(cmd, check=True)
    return plt.imread(png)


def impact_t_int():
    d = R.load(TEST)
    m = sio.loadmat(R.DATA / "Simulink" / f"{R.BY_ID[TEST][1]}.mat", squeeze_me=True, struct_as_record=False)["data"]
    raw = m.topLRaw.data
    tl, F = raw[:, -1], -R.LAM ** 3 * raw[:, 0] / 1e3
    w = (tl > 200) & (tl < 240)
    i = np.argmax(np.abs(F[w]))
    return float(d["to_int"](tl[w][i])), float(F[w][i])


def model():
    np_run, metas = pep.load_np(DUMP)
    tmp = Path(tempfile.mkdtemp(prefix="eq_impact_"))
    try:
        pep.stitch(DUMP, tmp, np_run, metas)
        meta = peq.read_meta(tmp)
        tags, xy = peq.read_nodes(tmp)
        disp_tags = peq.read_disp_nodes(tmp)
        lines, quads = peq.read_eles(tmp)
        t, ux, uy = peq.load_window_disp(tmp, disp_tags, meta.get("dispFiles", "window_disp.out").split())
        js = peq.load_spring_json(meta)
        mesh = peq._frame_mesh_parts(disp_tags, xy, lines, quads, js, meta)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return mesh, t, ux, uy


def draw(ax, mesh, ux, uy, amp, face=True, grey=False, under=False):
    """Mesh with displacements ux, uy (m) amplified by amp; grey=True draws a light-grey reference."""
    X, Y = mesh["X0"] + amp * ux, mesh["Y0"] + amp * uy
    q = mesh["q_idx"]
    z = 0 if under else 2
    if q.size:
        polys = np.stack([np.c_[X[q[:, i]], Y[q[:, i]]] for i in range(4)], axis=1)
        if grey:
            pc = PolyCollection(polys, facecolors="none", edgecolors="0.8", linewidths=0.3, zorder=z)
        else:
            if face == "light":
                fc = [(r, g, b, 0.22) for r, g, b, _ in mesh["face"]]
                pc = PolyCollection(polys, facecolors=fc, edgecolors="0.62", linewidths=0.2, zorder=z)
            else:
                pc = PolyCollection(polys, facecolors=mesh["face"] if face else "none",
                                    edgecolors="0.35" if face else "0.55", linewidths=0.2, zorder=z)
        ax.add_collection(pc)
    for ij, col, lw in ((mesh["ij_s"], "0.72" if grey else RED, 1.0 if grey else 1.25), (mesh["ij_z"], PURPLE, 0.5)):
        if len(ij) and not (grey and col == PURPLE):
            seg = np.stack([np.c_[X[ij[:, 0]], Y[ij[:, 0]]], np.c_[X[ij[:, 1]], Y[ij[:, 1]]]], axis=1)
            ax.add_collection(LineCollection(seg, colors=col, linewidths=lw, zorder=z + 1))


t_imp, f_imp = impact_t_int()
mesh, t, ux, uy = model()
k_imp = int(np.argmin(np.abs(t - t_imp)))
k_pre = int(np.argmin(np.abs(t - (t_imp - T_PRE))))
# pier-top node: highest structural node at x ~ 0 (top of the pier, below the deck frame)
ij = mesh["ij_s"]
snodes = np.unique(ij) if len(ij) else np.arange(mesh["X0"].size)
top = snodes[np.argmax(np.where(np.abs(mesh["X0"][snodes]) < 0.05, mesh["Y0"][snodes], -np.inf))]
win = np.flatnonzero((t >= t_imp) & (t <= t_imp + WIN_W))
k_w = win[np.argmax(np.abs(ux[win, top] - ux[k_pre, top]))]
du_x, du_y = ux[k_w] - ux[k_pre], uy[k_w] - uy[k_pre]
print(f"{TEST}: peak force {f_imp:.1f} kN at t_int = {t_imp:.2f} s; pre-wave state at {t[k_pre]:.2f} s; "
      f"largest wave-induced pier-top change {1e3 * du_x[top]:+.2f} mm at {t[k_w]:.2f} s "
      f"(node at y = {mesh['Y0'][top]:.2f} m); max |du| over the mesh {1e3 * np.abs(du_x).max():.2f} mm; "
      f"pre-wave pier-top offset {1e3 * ux[k_pre, top]:.1f} mm")

prs.configure_font()
CROP = {"wide": (0.12, 0.86), "close": (0.0, 0.78)}  # kept fraction of the frame height (top, bottom)
BLUE = "#0072B2"
H_WATER = 2.4  # m, prototype water depth above the mudline
fig = plt.figure(figsize=(6.0, 8.0), layout="constrained")
gs = fig.add_gridspec(2, 2, height_ratios=[1.3, 1.35])


def tag_in(ax, text, dark=True):
    """Panel tag inside the axes: white on a translucent dark box (photos) or bold black (model panels)."""
    if dark:
        ax.text(0.025, 0.975, text, transform=ax.transAxes, ha="left", va="top", color="white", zorder=10,
                bbox=dict(boxstyle="square,pad=0.25", fc="0.1", ec="none", alpha=0.55))
    else:
        ax.text(0.0, 1.0, text, transform=ax.transAxes, ha="left", va="top", color="black", zorder=10)


for j, (name, tf, flip, tag) in enumerate(FRAMES):
    ax = fig.add_subplot(gs[0, j])
    img = frame(name, tf, flip, tag)
    h0, h1 = (int(f * img.shape[0]) for f in CROP[tag])
    ax.imshow(img[h0:h1])
    ax.set_axis_off()
    tag_in(ax, rf"\textbf{{({'ab'[j]})}}")
ax_c = fig.add_subplot(gs[1, 0])
ax_d = fig.add_subplot(gs[1, 1], sharex=ax_c, sharey=ax_c)
zeros = np.zeros_like(ux[0])
draw(ax_c, mesh, zeros, zeros, 0.0, grey=True, under=True)
draw(ax_c, mesh, ux[k_imp], uy[k_imp], AMP)
draw(ax_d, mesh, ux[k_pre], uy[k_pre], AMP, grey=True, under=True)
draw(ax_d, mesh, ux[k_pre] + (AMP_W / AMP) * du_x, uy[k_pre] + (AMP_W / AMP) * du_y, AMP, face="light")
for ax in (ax_c, ax_d):
    ax.autoscale_view()
    ax.set_ylim(-21.0, 15.5)  # cropped just below the pile tips; room for the labels above the deck
    ax.set_aspect("equal")
    ax.set_axis_off()
x0, x1 = ax_c.get_xlim()
qs = mesh["q_idx"].ravel()
xs0, xs1 = float(mesh["X0"][qs].min()), float(mesh["X0"][qs].max())  # soil block extent
for ax in (ax_c, ax_d):  # water above the mudline, and the direction of the wave
    ax.add_patch(plt.Rectangle((xs0, 0.0), xs1 - xs0, H_WATER, fc=BLUE, ec="none", alpha=0.10, zorder=-1))
    ax.plot([xs0, xs1], [H_WATER, H_WATER], color=BLUE, lw=0.6, alpha=0.6, zorder=-1)
    ax.annotate("", (xs0 + 3.4, H_WATER + 1.0), (xs0, H_WATER + 1.0),
                arrowprops=dict(arrowstyle="-|>,head_length=0.5,head_width=0.25", color=BLUE, lw=1.2))
    ax.text(xs0, H_WATER + 1.6, "wave", color=BLUE, fontsize=7, ha="left", va="bottom")
y0 = ax_c.get_ylim()[0]
xb, yb = x0 + 1.2, y0 + 0.8
ax_c.plot([xb, xb + 5.0], [yb, yb], color="0.15", lw=1.4, solid_capstyle="butt", zorder=6)
ax_c.text(xb + 2.5, yb + 0.5, "5 m", ha="center", va="bottom", fontsize=7, color="0.15", zorder=6)
# (d): dimension at the deck between the grey (before the wave) and red (with the wave) positions, true size in mm
xg = mesh["X0"][top] + AMP * ux[k_pre][top]
yg = mesh["Y0"][top] + AMP * uy[k_pre][top]
xr = xg + AMP_W * du_x[top]
yr = yg + AMP_W * du_y[top]
ydim = max(yg, yr) + 1.6
for xx, yy in ((xg, yg), (xr, yr)):
    ax_d.plot([xx, xx], [yy + 0.25, ydim + 0.4], color="0.35", lw=0.5, zorder=7)
ax_d.annotate("", (xr, ydim), (xg, ydim), arrowprops=dict(arrowstyle="<|-|>,head_length=0.3,head_width=0.15",
                                                          color="0.2", lw=0.7, shrinkA=0, shrinkB=0), zorder=7)
ax_d.text(0.5 * (xg + xr), ydim + 0.35, f"+{1e3 * du_x[top]:.1f} mm (wave)", ha="center", va="bottom", fontsize=7,
          color="0.15", zorder=7, bbox=dict(fc="white", ec="none", alpha=0.8, pad=0.4))
tag_in(ax_c, r"\textbf{(c)}", dark=False)
tag_in(ax_d, r"\textbf{(d)}", dark=False)
ax_c.text(0.0, 0.955, f"at the peak force, amplified {AMP:.0f}" + r"$\times$", transform=ax_c.transAxes,
          ha="left", va="top", fontsize=7, color="0.3")
ax_d.text(0.0, 0.955, f"grey: before the wave ({AMP:.0f}" + r"$\times$" + ")\n"
          f"red: + change by the wave ({AMP_W:.0f}" + r"$\times$" + ")", transform=ax_d.transAxes,
          ha="left", va="top", fontsize=7, color="0.3", linespacing=1.15)
SHOW_E = False  # (e) pier-base moment-rotation inset; off: it competed with (c) and (d) (see fig_fri_hinge.py)
if SHOW_E:
    # (e) pier-base zero-length section, moment against rotation, placed after the layout is fixed
    hd = np.genfromtxt(DUMP / "pier_hinge_defo.out.0", invalid_raise=False)
    hf = np.genfromtxt(DUMP / "pier_hinge_force.out.0", invalid_raise=False)
    hd, hf = hd[np.isfinite(hd).all(1)], hf[np.isfinite(hf).all(1)]
    nh = min(len(hd), len(hf))
    th, Mh, tt = 1e3 * hd[:nh, 2], 1e-6 * hf[:nh, 2], hd[:nh, 0]
    pre_h = tt < t[k_pre]
    fig.canvas.draw()
    fig.set_layout_engine("none")
    pc, pd = ax_c.get_position(), ax_d.get_position()
    w_e, h_e = 0.27, 0.19
    x_e = 0.5 * (pc.x1 + pd.x0) - 0.5 * w_e
    y_mid = fig.transFigure.inverted().transform(ax_d.transData.transform((0.0, -9.0)))[1]
    y_e = y_mid - 0.5 * h_e
    ax_e = fig.add_axes([x_e, y_e, w_e, h_e], zorder=20)
    ax_e.set_facecolor("white")
    ax_e.patch.set_alpha(1.0)
    fig.add_artist(plt.Rectangle((x_e - 0.062, y_e - 0.045), w_e + 0.072, h_e + 0.072, transform=fig.transFigure,
                                 fc="white", ec="0.75", lw=0.6, zorder=19))  # white card behind (e), with its labels
    ax_e.plot(th[pre_h], Mh[pre_h], color="0.68", lw=0.5)
    ax_e.plot(th[~pre_h], Mh[~pre_h], color=RED, lw=1.6, solid_capstyle="round")
    ax_e.axhline(0, color="0.88", lw=0.4, zorder=0)
    ax_e.axvline(0, color="0.88", lw=0.4, zorder=0)
    ax_e.tick_params(labelsize=7, length=2.5, pad=1.5)
    ax_e.set_xticks([-10, 0, 10])
    ax_e.set_yticks([-5, 0, 5])
    ax_e.set_xlabel(r"$\theta$ (mrad)", fontsize=7, labelpad=1)
    ax_e.set_ylabel(r"$M$ (MN$\cdot$m)", fontsize=7, labelpad=1)
    ax_e.text(-0.2, 1.08, r"\textbf{(e)} pier-base section", transform=ax_e.transAxes, ha="left", va="bottom", fontsize=7)
    ax_e.text(0.97, 0.05, "grey: earthquake\nred: wave", transform=ax_e.transAxes, ha="right", va="bottom",
              fontsize=6.5, color="0.3", linespacing=1.1)
    # leader from (e) to the pier base in (d) (deformed position as drawn there)
    base = int(np.argmin(np.hypot(mesh["X0"], mesh["Y0"])))
    xd = mesh["X0"][base] + AMP * ux[k_pre][base] + AMP_W * du_x[base]
    yd = mesh["Y0"][base] + AMP * uy[k_pre][base] + AMP_W * du_y[base]
    from matplotlib.patches import ConnectionPatch  # noqa: E402
    fig.add_artist(ConnectionPatch(xyA=(1.0, 1.0), coordsA=ax_e.transAxes, xyB=(xd, yd), coordsB=ax_d.transData,
                                   arrowstyle="-|>,head_length=0.4,head_width=0.2", color="0.3", lw=0.7,
                                   ls=(0, (1, 1.5)), shrinkB=1.5, zorder=21))
    print(f"(e): pier-base node {base} at ({mesh['X0'][base]:.2f}, {mesh['Y0'][base]:.2f}) m")

fig.savefig(OUT.with_suffix(".png"), dpi=250)
fig.savefig(OUT.with_suffix(".pdf"), dpi=300)
print("wrote", OUT.with_suffix(".png"))
