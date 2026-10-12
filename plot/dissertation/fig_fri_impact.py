"""Solitary-wave impact: laboratory frames and the deformed numerical model (dissertation Ch. 6 draft; read-only on
data and videos). For visualization: the frames and the deformed shapes come from different tests with the same
solitary wave.

(a) wide view, video 2026-08-21-11h24m at 8.611 s (Test F10)
(b) close-up, video 2026-08-21-09h56m at 2:57.0, mirrored so that the wave travels in +x
(e) pier-base zero-length section, moment against rotation (grey: before the wave; red: from the wave on),
    small panel between (c) and (d) with a leader to the pier base in (d)
Layout 2 x 2: photos on top (cropped toward square), model panels below without axes (5 m scale bar).
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
            pc = PolyCollection(polys, facecolors=mesh["face"] if face else "none",
                                edgecolors="0.35" if face else "0.55", linewidths=0.2, zorder=z)
        ax.add_collection(pc)
    for ij, col, lw in ((mesh["ij_s"], "0.72" if grey else RED, 1.0), (mesh["ij_z"], PURPLE, 0.5)):
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
CROP = {"wide": (0.20, 0.76), "close": (0.02, 0.62)}  # kept fraction of the frame height (top, bottom)
fig = plt.figure(figsize=(6.0, 8.3), layout="constrained")
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.75])
for j, (name, tf, flip, tag) in enumerate(FRAMES):
    ax = fig.add_subplot(gs[0, j])
    img = frame(name, tf, flip, tag)
    h0, h1 = (int(f * img.shape[0]) for f in CROP[tag])
    ax.imshow(img[h0:h1])
    ax.set_axis_off()
    ax.set_title(rf"\textbf{{({'ab'[j]})}}", loc="left")
ax_c = fig.add_subplot(gs[1, 0])
ax_d = fig.add_subplot(gs[1, 1], sharex=ax_c, sharey=ax_c)
zeros = np.zeros_like(ux[0])
draw(ax_c, mesh, zeros, zeros, 0.0, grey=True, under=True)
draw(ax_c, mesh, ux[k_imp], uy[k_imp], AMP)
draw(ax_d, mesh, ux[k_pre], uy[k_pre], AMP, grey=True, under=True)
draw(ax_d, mesh, ux[k_pre] + (AMP_W / AMP) * du_x, uy[k_pre] + (AMP_W / AMP) * du_y, AMP, face=False)
for ax in (ax_c, ax_d):
    ax.autoscale_view()
    ax.set_ylim(top=13.0)  # room for the labels above the deck
    ax.set_aspect("equal")
    ax.set_axis_off()
x0, x1 = ax_c.get_xlim()
y0 = ax_c.get_ylim()[0]
xb, yb = x0 + 0.3, y0 + 0.4
ax_c.plot([xb, xb + 5.0], [yb, yb], color="0.2", lw=1.2, solid_capstyle="butt")
ax_c.text(xb + 2.5, yb + 0.5, "5 m", ha="center", va="bottom", fontsize=7, color="0.2")
ax_c.set_title(r"\textbf{(c)}", loc="left")
ax_d.set_title(r"\textbf{(d)}", loc="left")
ax_c.text(0.02, 0.99, f"at the peak force\namplified {AMP:.0f}" + r"$\times$", transform=ax_c.transAxes,
          ha="left", va="top", fontsize=7, color="0.25")
ax_d.text(0.02, 0.99, f"grey: before the wave ({AMP:.0f}" + r"$\times$" + ")\n"
          f"red: + change by the wave ({AMP_W:.0f}" + r"$\times$" + ")", transform=ax_d.transAxes,
          ha="left", va="top", fontsize=7, color="0.25")
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
w_e, h_e = 0.16, 0.085
x_e = 0.5 * (pc.x1 + pd.x0) - 0.5 * w_e + 0.035
yb = ax_d.transData.transform((0.0, 0.0))
y_e = fig.transFigure.inverted().transform(yb)[1] + 0.03
ax_e = fig.add_axes([x_e, y_e, w_e, h_e])
ax_e.plot(th[pre_h], Mh[pre_h], color="0.72", lw=0.4)
ax_e.plot(th[~pre_h], Mh[~pre_h], color=RED, lw=1.1)
ax_e.axhline(0, color="0.85", lw=0.3, zorder=0)
ax_e.axvline(0, color="0.85", lw=0.3, zorder=0)
ax_e.tick_params(labelsize=6, length=2, pad=1)
ax_e.set_xticks([-10, 0, 10])
ax_e.set_yticks([-5, 0, 5])
ax_e.set_xlabel(r"$\theta$ (mrad)", fontsize=6.5, labelpad=1)
ax_e.set_ylabel(r"$M$ (MN$\cdot$m)", fontsize=6.5, labelpad=2)
ax_e.set_title(r"\textbf{(e)} pier base", loc="left", fontsize=7, pad=2)
ax_e.patch.set_alpha(0.9)
# leader from (e) to the pier base in (d) (deformed position as drawn there)
base = int(np.argmin(np.hypot(mesh["X0"], mesh["Y0"])))
xd = mesh["X0"][base] + AMP * ux[k_pre][base] + AMP_W * du_x[base]
yd = mesh["Y0"][base] + AMP * uy[k_pre][base] + AMP_W * du_y[base]
from matplotlib.patches import ConnectionPatch  # noqa: E402
fig.add_artist(ConnectionPatch(xyA=(1.0, 0.0), coordsA=ax_e.transAxes, xyB=(xd, yd), coordsB=ax_d.transData,
                               arrowstyle="-|>,head_length=0.35,head_width=0.18", color="0.35", lw=0.6,
                               ls=(0, (1, 1.5)), shrinkB=1.5))
print(f"(e): pier-base node {base} at ({mesh['X0'][base]:.2f}, {mesh['Y0'][base]:.2f}) m")
fig.savefig(OUT.with_suffix(".png"), dpi=250)
fig.savefig(OUT.with_suffix(".pdf"), dpi=300)
print("wrote", OUT.with_suffix(".png"))
