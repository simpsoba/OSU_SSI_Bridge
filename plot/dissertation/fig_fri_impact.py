"""Solitary-wave impact: laboratory frames and the deformed numerical model (dissertation Ch. 6 draft; read-only on
data and videos). For visualization: the frames and the deformed shape come from different tests with the same
solitary wave.

(a) wide view, video 2026-08-21-11h24m at 8.611 s (Test F10)
(b) close-up, video 2026-08-21-09h56m at 2:57.0, mirrored so that the wave travels in +x
(c) deformed near-field model of Test F30 (r-17, full-domain recorders) at the peak load-cell force of the solitary
    wave, displacements amplified by AMP; undeformed mesh in light grey
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
AMP = 20.0


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
    import scipy.io as sio
    m = sio.loadmat(R.DATA / "Simulink" / f"{R.BY_ID[TEST][1]}.mat", squeeze_me=True, struct_as_record=False)["data"]
    raw = m.topLRaw.data
    tl, F = raw[:, -1], -R.LAM ** 3 * raw[:, 0] / 1e3
    w = (tl > 200) & (tl < 240)
    tp = tl[w][np.argmax(np.abs(F[w]))]
    return float(d["to_int"](tp)), float(F[w][np.argmax(np.abs(F[w]))])


def mesh_at(t_target):
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
    k = int(np.argmin(np.abs(t - t_target)))
    return mesh, t[k], ux[k], uy[k]


t_imp, f_imp = impact_t_int()
mesh, t_k, ux, uy = mesh_at(t_imp)
print(f"{TEST}: impact at t_int = {t_imp:.2f} s (peak force {f_imp:.1f} kN); recorder frame t = {t_k:.2f} s; "
      f"max |ux| = {1e3 * np.abs(ux).max():.0f} mm")

prs.configure_font()
fig = plt.figure(figsize=(6.0, 3.9), layout="constrained")
gs = fig.add_gridspec(1, 3, width_ratios=[0.95, 0.95, 1.5])
for j, (name, t, flip, tag) in enumerate(FRAMES):
    ax = fig.add_subplot(gs[j])
    ax.imshow(frame(name, t, flip, tag))
    ax.set_axis_off()
    ax.set_title(rf"\textbf{{({'ab'[j]})}}", loc="left")
ax = fig.add_subplot(gs[2])
X0, Y0 = mesh["X0"], mesh["Y0"]
X, Y = X0 + AMP * ux, Y0 + AMP * uy
for xs, ys, alpha, z in ((X0, Y0, 0.0, 0), (X, Y, 1.0, 1)):
    q = mesh["q_idx"]
    if q.size:
        polys = np.stack([np.c_[xs[q[:, i]], ys[q[:, i]]] for i in range(4)], axis=1)
        if alpha:
            ax.add_collection(PolyCollection(polys, facecolors=mesh["face"], edgecolors="0.35", linewidths=0.2,
                                             zorder=z))
        else:
            ax.add_collection(PolyCollection(polys, facecolors="none", edgecolors="0.8", linewidths=0.3, zorder=z))
    for ij, col, lw in ((mesh["ij_s"], "#C0390B" if alpha else "0.75", 1.0), (mesh["ij_z"], "#6A3D9A", 0.5)):
        if len(ij) and (alpha or col != "#6A3D9A"):
            seg = np.stack([np.c_[xs[ij[:, 0]], ys[ij[:, 0]]], np.c_[xs[ij[:, 1]], ys[ij[:, 1]]]], axis=1)
            ax.add_collection(LineCollection(seg, colors=col, linewidths=lw, zorder=z + 1))
ax.autoscale_view()
ax.set_aspect("equal")
ax.set_xlabel(r"$x$ (m)")
ax.set_ylabel(r"$y$ (m)")
ax.set_title(r"\textbf{(c)}" + f" deformed shape, amplified {AMP:.0f}" + r"$\times$", loc="left")
fig.savefig(OUT.with_suffix(".png"), dpi=250)
fig.savefig(OUT.with_suffix(".pdf"), dpi=300)
print("wrote", OUT.with_suffix(".png"))
