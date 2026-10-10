#!/usr/bin/env python3
"""
Goals
-----
Paper figure for the offline lumped-plasticity, pier-base-hold case.
Each station is its own axes: Newmark (grey, solid) and MKR-α (black, dashed,
ρ_∞² = 0.5). The mesh is a left column; the histories are stacked to its
right on a shared 30–180 s time axis. Newmark stops early. A zoom column on
the far right magnifies each history around its peak |u| (box on the history
marks the zoomed window).

  (a) near-field mesh, stations marked with the panel letters
  (b) pier top
  (c–d) center pile at 1.0 m and 10.1 m depth
  (e) soil column at x = 5.8 m, mudline

u_x is the recorder value (m on file, mm on the figure). The relative version
plots Δu_x = u_x − u_x,base, where u_x,base is the Shin base primary node
(soil_base_primary.out; every base node is equalDOF to it in UX, so it carries
the input motion). One shared y-axis label. Pier top and soil have their own
limits; the pile panels share one in the relative version, while the pile head
(c) gets its own in the absolute version (it barely moves relative to the
ground there). Each limit is a round value that is also the outermost tick. Axes styling follows
PlotGroundMotionFigure.py (usetex CM Sans via PlotResponseSpectrum.configure_font,
left/bottom spines, light major grid).

Writes
------
  plot/out/eq_offline/compare/lumped_hold_integrator_every.{png,pdf}
  plot/out/eq_offline/compare/lumped_hold_integrator_every_rel.{png,pdf}

Usage
-----
  python plot/PlotLumpedHoldIntegratorPaper.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator

import math

import PlotResponseSpectrum as prs

from PlotModelSketchPaper import (
    PIER_COLOR,
    SOIL_ALPHA,
    SOIL_EDGE,
    SOIL_FILL,
    STIFF_COLOR,
)

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "plot" / "out" / "eq_offline"
# Serial dumps compared: UmfPack, dt = 0.00756 s, 180 s, recorded at every
# analysis step (eqRecDt 0). Newmark 0.5 0.25 and MKRAlphaExplicitMultiSOE 0.5.
# Earlier pairs (CudaMKRAlpha, 0.01 s recorders) wrote lumped_hold_integrator_cudamkr.*.
MKR_DUMP = "lumped_hold_mkr_every"
NM_DUMP = "lumped_hold_newmark_dt10_every"
OUT = ROOT / "compare" / "lumped_hold_integrator_every"

# Center pile: cap BC (head) and iy = 10 (mid-depth). See BuildPilesNodes.tcl.
PILE_TAGS = (1028, 2110)
PIER_TAG = 5
# Soil station off the shafts (piles sit at x = 0 and ±1.83 m).
SOIL_X = 5.7912
# Mudline only (the 10.1 m soil history tracks the 10.1 m pile one closely).
SOIL_Y = (0.0,)

# Newmark: thin mid-grey solid line under MKR-α: black dashes with open gaps, so the grey
# reads through the gaps where the two agree and as its own line where they differ (B&W
# safe). A wide light-grey line was hidden under the dashes in dense loops and zooms.
# Dash length is in points.
NM_COLOR = "#8c8c8c"
NM_LW = 1.2
MKR_LW = 0.7
MKR_LS = (0, (2.0, 2.2))
# Station letters on the mesh sketch (gray, regular weight) vs bold black panel tags.
STATION_GRAY = "#6e6e6e"

FIG_W = 6.0
# Mesh (true scale, no axes) is a left column, top-aligned with (b); the legend
# sits under it (LEG_DROP below the mesh, clear of the 5 m bar). Histories are stacked to the right and share one
# time axis; two zoom columns at the far right, without u axes (annotate_zoom
# marks the peak value and the largest gap instead).
ML, MR, MB, MT = 0.03, 0.06, 0.40, 0.08
COL_GAP = 0.66
LEG_DROP = 0.40
HIST_H = 0.75
HIST_GAP = 0.16
N_HIST = 4
STACK_H = N_HIST * HIST_H + (N_HIST - 1) * HIST_GAP
MESH_H = 1.80
ZOOM_GAP = 0.10
ZOOM_GAP2 = 0.06
ZOOM_W = 0.95
# Text size in the zooms (same 9 pt as the rest of the figure).
ZOOM_TICK_FS = 9
# Zoom window width (s). Windows start on an odd whole second, so every zoom
# has its two time ticks (t0 + 1, t0 + 3) at the same place on its axis.
ZOOM_SPAN = 4.0
# window_nodes extent used for the mesh (m), matching set_xlim / set_ylim below.
MESH_DX = 17.10
MESH_DY = 35.40
MESH_W = MESH_H * MESH_DX / MESH_DY
HIST_X = ML + MESH_W + COL_GAP
HIST_W = FIG_W - HIST_X - ZOOM_GAP - ZOOM_W - ZOOM_GAP2 - ZOOM_W - MR
ZOOM_X = HIST_X + HIST_W + ZOOM_GAP
ZOOM2_X = ZOOM_X + ZOOM_W + ZOOM_GAP2
FIG_H = MT + STACK_H + MB
# Time window shown on the histories (s).
T_MIN, T_MAX = 30.0, 180.0


def read_xy(path: Path) -> dict[int, tuple[float, float]]:
    """Node tag -> (x, y) in m from window_nodes.txt."""
    out: dict[int, tuple[float, float]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        tag_s, xs, ys = line.split()[:3]
        out[int(float(tag_s))] = (float(xs), float(ys))
    return out


def read_dt(eq: Path) -> float:
    """Analysis time step (s): dtAnalysis in window_meta.txt."""
    for line in (eq / "window_meta.txt").read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0] == "dtAnalysis":
            return float(parts[1])
    raise SystemExit(f"no dtAnalysis in {eq / 'window_meta.txt'}")


def read_order(path: Path) -> list[int]:
    """disp_nodes.txt tag order (column order of window_disp_*.out)."""
    tags: list[int] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        tags.append(int(float(line.split()[0])))
    return tags


def read_eles(path: Path) -> list[tuple[int, list[int]]]:
    """window_eles.txt as (eleTag, node tags)."""
    rows: list[tuple[int, list[int]]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        parts = [int(float(p)) for p in line.split()]
        rows.append((parts[0], parts[1:]))
    return rows


def read_layers(path: Path) -> dict[int, str]:
    """eleTag -> layer name from window_quads.txt."""
    layers: dict[int, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        ele_s, layer, _col = line.split()[:3]
        layers[int(ele_s)] = layer
    return layers


def chunk_of(index: int, chunk: int = 250) -> tuple[int, int]:
    """Return (file index, column index inside that file) for one node."""
    return index // chunk, index % chunk


def load_ux(path: Path, local_idx: list[int]) -> tuple[np.ndarray, np.ndarray]:
    """
    Horizontal displacement from one window_disp file.

    Args: path, local node indices inside this chunk (0-based)
    Returns: t (s), ux (n, k) in m, as recorded
    """
    cols = [0] + [1 + 2 * i for i in local_idx]
    arr = np.loadtxt(path, usecols=cols)
    t = arr[:, 0]
    return t, arr[:, 1:]


def load_pier_ux(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Pier-node recorder: time, UX, UY, RZ. UX in m, as recorded."""
    arr = np.loadtxt(path)
    return arr[:, 0], arr[:, 1]


def series_for(eq: Path, order: list[int], tags: list[int]) -> tuple[np.ndarray, np.ndarray]:
    """
    UX histories for ``tags``, concatenated across window_disp chunks.

    Returns: t (s), ux (n, n_tags) in m
    """
    loc = {tag: chunk_of(order.index(tag)) for tag in tags}
    files = sorted({fi for fi, _ in loc.values()})
    t_ref: np.ndarray | None = None
    cols: dict[int, np.ndarray] = {}
    for fi in files:
        want = [tag for tag, (f, _j) in loc.items() if f == fi]
        local = [loc[tag][1] for tag in want]
        t, ux = load_ux(eq / f"window_disp_{fi:02d}.out", local)
        if t_ref is None:
            t_ref = t
        for k, tag in enumerate(want):
            cols[tag] = ux[:, k]
    assert t_ref is not None
    return t_ref, np.column_stack([cols[tag] for tag in tags])


def soil_tags(xy: dict[int, tuple[float, float]]) -> list[int]:
    """Mudline and mid-depth nodes on the x ≈ 5.8 m soil line."""
    soil = {
        tag: (x, y)
        for tag, (x, y) in xy.items()
        if 10000 <= tag < 11000 and abs(x - SOIL_X) < 0.05
    }
    if len(soil) < len(SOIL_Y):
        raise SystemExit(f"soil column at x={SOIL_X} has {len(soil)} nodes")
    picked: list[int] = []
    for y_want in SOIL_Y:
        tag = min(soil, key=lambda t: abs(soil[t][1] - y_want))
        picked.append(tag)
    return picked


def member_kind(nodes: list[int], xy: dict[int, tuple[float, float]]) -> str | None:
    """
    Beam family for a 2-node element, or None to skip.

    Zero-length springs and the soil-spring duplicates (tags ≥ 20000) are omitted.
    """
    if len(nodes) != 2:
        return None
    n1, n2 = nodes
    if n1 >= 20000 or n2 >= 20000:
        return None
    x1, y1 = xy[n1]
    x2, y2 = xy[n2]
    if (x1 - x2) ** 2 + (y1 - y2) ** 2 < 1e-8:
        return None
    pair = {n1, n2}
    if pair <= {1, 2, 4, 5}:
        return "pier"
    if any(2000 <= n < 2300 for n in pair) and all(
        2000 <= n < 2300 or n in (1027, 1028, 1029) for n in pair
    ):
        return "pile"
    if any(n >= 3000 for n in pair):
        return "deck"
    if all(1000 <= n < 2000 or n in (1, 2) for n in pair):
        return "cap"
    return None


def draw_mesh(
    ax: plt.Axes,
    xy: dict[int, tuple[float, float]],
    eles: list[tuple[int, list[int]]],
    layers: dict[int, str],
    marks: list[tuple[int, str]],
) -> None:
    """Near-field mesh, structure, and the stations used in (b)–(f)."""
    polys = []
    faces = []
    for ele, nodes in eles:
        if len(nodes) != 4:
            continue
        polys.append([xy[n] for n in nodes])
        faces.append(SOIL_FILL.get(layers.get(ele, ""), "#d9d0c1"))
    quads = PolyCollection(
        polys,
        facecolors=faces,
        edgecolors=SOIL_EDGE,
        linewidths=0.15,
        alpha=SOIL_ALPHA,
        zorder=1,
    )
    ax.add_collection(quads)

    buckets: dict[str, list[list[tuple[float, float]]]] = {
        "pile": [],
        "cap": [],
        "deck": [],
        "pier": [],
    }
    for _ele, nodes in eles:
        kind = member_kind(nodes, xy)
        if kind is None:
            continue
        buckets[kind].append([xy[nodes[0]], xy[nodes[1]]])
    style = {
        "pile": (STIFF_COLOR, 1.15),
        "cap": (STIFF_COLOR, 1.05),
        "deck": (STIFF_COLOR, 1.05),
        "pier": (PIER_COLOR, 1.7),
    }
    for kind, segs in buckets.items():
        color, lw = style[kind]
        ax.add_collection(
            LineCollection(segs, colors=color, linewidths=lw, zorder=3, capstyle="round")
        )

    # Panel letters sit off the marker: piles to the left of the shafts, soil to the right.
    for tag, letter in marks:
        x, y = xy[tag]
        ax.scatter(
            [x],
            [y],
            s=10,
            c="#111111",
            marker="o",
            zorder=5,
            linewidths=0.3,
            edgecolors="white",
        )
        # Letters sit off the marker with a thin leader (as in the hysteresis figure).
        if tag == PIER_TAG:
            # Below-right of the soffit so the letter clears the deck chord.
            xytext, ha = (5.6, y - 2.2), "left"
        elif 2000 <= tag < 2300 or tag in (1027, 1028, 1029):
            xytext, ha = (-4.6, y), "right"
        else:
            # Soil: above-left of the marker, clear of the mesh edge.
            xytext, ha = (x - 1.8, y + 2.4), "right"
        textcoords = "data"
        ax.annotate(
            rf"({letter})",  # station reference: gray, regular weight (panel tags are bold black)
            xy=(x, y),
            xytext=xytext,
            textcoords=textcoords,
            ha=ha,
            va="center",
            color=STATION_GRAY,
            zorder=6,
            clip_on=False,
            arrowprops=dict(arrowstyle="-", lw=0.4, color=STATION_GRAY, shrinkA=1.5, shrinkB=1.5),
            # White backing box: a withStroke halo shrinks usetex glyphs.
            bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none", alpha=0.85),
        )

    ax.set_xlim(-8.55, 8.55)
    ax.set_ylim(-24.55, 10.85)
    ax.set_aspect("equal", adjustable="box")
    # No axes on the sketch: station coordinates are given in the history titles.
    ax.set_axis_off()
    # 5 m scale bar under the mesh (outside the data limits, so not clipped).
    x0, x1, yb = -7.5, -2.5, -26.3
    tick = 0.6
    ax.plot([x0, x1], [yb, yb], color="black", lw=0.9, clip_on=False, solid_capstyle="butt")
    for xe in (x0, x1):
        ax.plot([xe, xe], [yb - tick, yb + tick], color="black", lw=0.9, clip_on=False)
    ax.text(x1 + 0.8, yb, "5 m", ha="left", va="center", clip_on=False)
    # Panel letter inside the (axis-free) sketch, top left, clear of the deck.
    ax.text(0.0, 1.0, r"\textbf{(a)}", transform=ax.transAxes, ha="left", va="top")


# Space the in-axes panel title takes on a history (in), kept free of box tags.
TITLE_W_IN, TITLE_H_IN = 1.10, 0.18


def panel_title(ax: plt.Axes, letter: str, rest: str) -> None:
    """Bold panel letter and station name, inside the axes at the top left."""
    text = rf"\textbf{{({letter})}}\quad {rest}"
    ax.text(
        0.012,
        0.97,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        zorder=7,
        bbox=dict(boxstyle="square,pad=0.06", fc="white", ec="none", alpha=0.85),
    )


def plot_pair(ax: plt.Axes, t_m: np.ndarray, u_m: np.ndarray, t_n: np.ndarray, u_n: np.ndarray) -> None:
    """One station: Newmark solid grey, MKR-α dashed black. Arguments in m, drawn in mm."""
    ax.plot(t_n, 1000.0 * u_n, color=NM_COLOR, ls="-", lw=NM_LW, zorder=2, solid_capstyle="butt")
    ax.plot(t_m, 1000.0 * u_m, color="black", ls=MKR_LS, lw=MKR_LW, zorder=3)


def style_history(ax: plt.Axes) -> None:
    """Time-history axes as in PlotGroundMotionFigure.py: open box, light major grid."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, which="major", color="0.85", lw=0.5, alpha=0.25)
    ax.set_axisbelow(True)


def round_lim(peak: float) -> float:
    """
    Axis limit that is also the outermost tick: peak rounded up to half a decade.

    E.g. 18 -> 20, 199 -> 200, 230 -> 250, 305 -> 350 (mm).
    """
    if peak <= 0.0 or not math.isfinite(peak):
        return 1.0
    step = 0.5 * 10.0 ** math.floor(math.log10(peak))
    return step * math.ceil(peak / step - 1e-9)


def window_peak_mm(
    t_m: np.ndarray, u_m: np.ndarray, t_n: np.ndarray, u_n: np.ndarray
) -> float:
    """Peak |u| in mm for one station, from both integrators inside [T_MIN, T_MAX]."""
    in_m = (t_m >= T_MIN) & (t_m <= T_MAX)
    in_n = (t_n >= T_MIN) & (t_n <= T_MAX)
    peak_m = float(np.max(np.abs(u_m[in_m]))) if in_m.any() else 0.0
    peak_n = float(np.max(np.abs(u_n[in_n]))) if in_n.any() else 0.0
    return 1000.0 * max(peak_m, peak_n)


def pick_windows(
    rows: list[tuple[str, str, str, np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
    group_lim: dict[str, float],
) -> tuple[float, float]:
    """
    One window start per zoom column, shared by every station (so each column
    shares its time axis). Starts are odd whole seconds in [T_MIN, T_MAX - ZOOM_SPAN].

      column 1   largest sum over stations of peak |u_MKR| / station y-limit
      column 2   largest sum over stations of max |u_Newmark - u_MKR| / that
                 station's max over the whole record, among windows clear of column 1

    Newmark is interpolated onto the MKR record times only for this choice.
    Returns: (t0 column 1, t0 column 2), s.
    """
    starts = np.arange(T_MIN + 1.0, T_MAX - ZOOM_SPAN + 1e-9, 2.0)
    s1 = np.zeros_like(starts)
    s2 = np.zeros_like(starts)
    for group, _l, _r, t_m, u_m, t_n, u_n in rows:
        u = 1000.0 * u_m
        d = 1000.0 * np.abs(np.interp(t_m, t_n, u_n) - u_m)
        full = (t_m >= T_MIN) & (t_m <= T_MAX)
        d_full = max(float(d[full].max()), 1e-12)
        for i, t0 in enumerate(starts):
            m = (t_m >= t0) & (t_m <= t0 + ZOOM_SPAN)
            s1[i] += float(np.abs(u[m]).max()) / group_lim[group]
            s2[i] += float(d[m].max()) / d_full
    w1 = float(starts[int(np.argmax(s1))])
    clear = np.abs(starts - w1) >= ZOOM_SPAN
    w2 = float(starts[clear][int(np.argmax(s2[clear]))])
    return w1, w2


def zoom_window(
    t0: float, t_m: np.ndarray, u_m: np.ndarray, t_n: np.ndarray, u_n: np.ndarray
) -> tuple[float, float, float, float]:
    """
    ZOOM_SPAN window starting at t0 (an odd whole second, so the two ticks at
    t0 + 1 and t0 + 3 fall on even seconds at 1/4 and 3/4 of the axis).

    Returns: t0, t1 (s), y0, y1 (mm), y padded 10 % of the range.
    """
    t1 = t0 + ZOOM_SPAN
    in_m = (t_m >= t0) & (t_m <= t1)
    in_n = (t_n >= t0) & (t_n <= t1)
    vals = np.concatenate([1000.0 * u_m[in_m], 1000.0 * u_n[in_n]])
    lo, hi = float(vals.min()), float(vals.max())
    pad = 0.10 * (hi - lo)
    return t0, t1, lo - pad, hi + pad


def mm_text(v: float) -> str:
    """mm value for an annotation: one decimal below 20 mm, whole mm above."""
    return f"{v:.1f}" if abs(v) < 20.0 else f"{v:.0f}"


def annotate_zoom(
    ax_z: plt.Axes,
    t0: float,
    t1: float,
    t_m: np.ndarray,
    u_m: np.ndarray,
    t_n: np.ndarray,
    u_n: np.ndarray,
) -> tuple[float, float]:
    """
    Peak and valley values per zoom (no u axis). For the window maximum and the
    window minimum, each curve's own extreme is marked with a dot and the pair of
    values (Newmark grey, MKR-α black, upper value on top) is written in white
    space nearer its own point than the other one. When a pair finds no free spot,
    the zoom's u range is extended on that side (peak: up, valley: down) and the
    layout is redone, up to a limit; pairs still without room are left out.

    Returns: the final (y0, y1) of ax_z.
    """
    m = (t_m >= t0) & (t_m <= t1)
    tw, um = t_m[m], 1000.0 * u_m[m]
    mn = (t_n >= t0) & (t_n <= t1)
    twn, unn = t_n[mn], 1000.0 * u_n[mn]
    tf = np.linspace(t0, t1, 1600)
    uf = np.concatenate([np.interp(tf, tw, um), np.interp(tf, twn, unn)])
    tff = np.concatenate([tf, tf])
    fig = ax_z.figure
    rend = fig.canvas.get_renderer()
    bh = 0.26 / HIST_H
    tag_rect = (0.0, ZOOM_TAG_W_IN / ZOOM_W, 1.0 - ZOOM_TAG_H_IN / HIST_H, 1.0)

    def overlaps(a, b) -> bool:
        return not (a[1] <= b[0] or b[1] <= a[0] or a[3] <= b[2] or b[3] <= a[2])

    # The two extremes (index into MKR and Newmark samples) and their text.
    ext = []
    for pick in (np.argmax, np.argmin):
        km, kn = int(pick(um)), int(pick(unn))
        pts = [(unn[kn], "0.35", NM_COLOR, twn[kn]), (um[km], "black", "black", tw[km])]
        order = sorted(pts, key=lambda q: -q[0])
        w_px = 0.0
        for val, *_ in order:
            probe = ax_z.text(0, 0, rf"${mm_text(val)}$", fontsize=ZOOM_TICK_FS, transform=ax_z.transAxes)
            w_px = max(w_px, probe.get_window_extent(rend).width)
            probe.remove()
        bw = (w_px / fig.dpi + 0.04) / ZOOM_W
        tp = 0.5 * (twn[kn] + tw[km])
        up = 0.5 * (um[km] + unn[kn])
        ext.append((pts, order, bw, tp, up))

    def layout(y0: float, y1: float):
        fx = lambda t: (t - t0) / (t1 - t0)  # noqa: E731
        fy = lambda u: (u - y0) / (y1 - y0)  # noqa: E731
        xs, ys = fx(tff), fy(uf)
        taken = [tag_rect]
        out = []
        for ie, (_pts, _order, bw, tp, up) in enumerate(ext):
            px, py = fx(tp), fy(up)
            ox, oy = fx(ext[1 - ie][3]), fy(ext[1 - ie][4])
            best, best_n = None, None
            for lx in np.linspace(0.01, 0.99 - bw, 25):
                for bt in np.linspace(0.01, 0.99 - bh, 25):
                    r = (lx, lx + bw, bt, bt + bh)
                    if any(overlaps(r, q) for q in taken):
                        continue
                    cx, cy = 0.5 * (r[0] + r[1]), 0.5 * (r[2] + r[3])
                    dist = float(np.hypot((cx - px) * ZOOM_W, (cy - py) * HIST_H))
                    if dist > float(np.hypot((cx - ox) * ZOOM_W, (cy - oy) * HIST_H)):
                        continue
                    pad = 0.004
                    inside = (
                        (xs >= r[0] - pad) & (xs <= r[1] + pad)
                        & (ys >= r[2] - pad) & (ys <= r[3] + pad)
                    )
                    n = (int(inside.sum()), dist)
                    if best_n is None or n < best_n:
                        best, best_n = r, n
            ok = best is not None and best_n[0] == 0
            if ok:
                taken.append(best)
            out.append(best if ok else None)
        return out

    y0, y1 = ax_z.get_ylim()
    span0 = y1 - y0
    for _ in range(4):
        placed = layout(y0, y1)
        if all(r is not None for r in placed):
            break
        grow = 0.18 * (y1 - y0)
        if placed[0] is None and (y1 - y0) < 1.8 * span0:
            y1 += grow
        if placed[1] is None and (y1 - y0) < 1.8 * span0:
            y0 -= grow
        if (y1 - y0) >= 1.8 * span0:
            placed = layout(y0, y1)
            break
    ax_z.set_ylim(y0, y1)
    for (pts, order, _bw, tp, up), r in zip(ext, placed):
        if r is None:
            continue
        # Thin leader from the value pair to its point when they are not adjacent.
        px = (tp - t0) / (t1 - t0)
        py = (up - y0) / (y1 - y0)
        ex = min(max(px, r[0]), r[1])
        ey = min(max(py, r[2]), r[3])
        if np.hypot((ex - px) * ZOOM_W, (ey - py) * HIST_H) > 0.20:
            # Stop just short of the dot so it stays visible.
            k = 1.0 - 0.035 / max(np.hypot((ex - px) * ZOOM_W, (ey - py) * HIST_H), 1e-9)
            ax_z.plot(
                [ex, px + (ex - px) * (1.0 - k)],
                [ey, py + (ey - py) * (1.0 - k)],
                transform=ax_z.transAxes,
                color="0.3",
                lw=0.45,
                zorder=7,
            )
        for val, _tc, mc, tpp in pts:
            ax_z.plot([tpp], [val], marker="o", ms=3.0, mfc=mc, mec="black", mew=0.5, zorder=8)
        for line, (val, tcol, _mc, _tp2) in enumerate(order):
            ax_z.text(
                r[0],
                r[2] + bh * (0.75 if line == 0 else 0.25),
                rf"${mm_text(val)}$",
                transform=ax_z.transAxes,
                ha="left",
                va="center",
                fontsize=ZOOM_TICK_FS,
                color=tcol,
                zorder=9,
            )
    return y0, y1


# Zoom tag "(b.1)" inside the zoom at the top left (in), kept free of values.
ZOOM_TAG_W_IN, ZOOM_TAG_H_IN = 0.34, 0.17


def draw_zoom(
    ax_h: plt.Axes,
    ax_z: plt.Axes,
    t0: float,
    t_m: np.ndarray,
    u_m: np.ndarray,
    t_n: np.ndarray,
    u_n: np.ndarray,
    tag: str,
    time_labels: bool,
) -> tuple[float, float, float, float]:
    """
    Magnified window [t0, t0 + ZOOM_SPAN] in ax_z, tagged inside with ``tag``
    (e.g. "(b.1)"); the same window outlined on ax_h. Time labels only when
    ``time_labels`` (bottom row: each column shares one window). Returns the box
    (t0, t1, b0, b1); the caller tags the boxes on ax_h with place_box_tags.
    """
    t0, t1, y0, y1 = zoom_window(t0, t_m, u_m, t_n, u_n)
    # Extra headroom so the curves pass under the tag in the top-left corner.
    tag_w = ZOOM_TAG_W_IN / ZOOM_W + 0.03
    tag_top = 1.0 - ZOOM_TAG_H_IN / HIST_H - 0.04
    corner = []
    for t, u in ((t_m, u_m), (t_n, u_n)):
        m = (t >= t0) & (t <= t0 + tag_w * (t1 - t0))
        corner.append(1000.0 * u[m])
    v_corner = float(np.concatenate(corner).max())
    if (v_corner - y0) / (y1 - y0) > tag_top:
        y1 = y0 + (v_corner - y0) / tag_top
    plot_pair(ax_z, t_m, u_m, t_n, u_n)
    ax_z.set_xlim(t0, t1)
    ax_z.set_ylim(y0, y1)
    ax_z.grid(True, which="major", color="0.85", lw=0.5, alpha=0.25)
    ax_z.set_axisbelow(True)
    # Full frame so the inset reads as a magnifier of the boxed window.
    for side in ("top", "right"):
        ax_z.spines[side].set_visible(True)
    ax_z.set_xticks([t0 + 1.0, t0 + 3.0])
    # No u axis on the zooms (each has its own scale): both values are written at
    # the point of largest Newmark - MKR difference instead. Header: the tag.
    ax_z.yaxis.set_major_locator(MaxNLocator(nbins=3))
    ax_z.tick_params(axis="y", length=0, labelleft=False, labelright=False)
    ax_z.tick_params(axis="x", labelsize=ZOOM_TICK_FS, labelbottom=time_labels)
    y0, y1 = annotate_zoom(ax_z, t0, t1, t_m, u_m, t_n, u_n)
    ax_z.text(
        0.03,
        0.96,
        rf"\textbf{{{tag}}}",
        transform=ax_z.transAxes,
        ha="left",
        va="top",
        zorder=12,
        bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none", alpha=0.85),
    )
    # Outline clamped to the history limits, so the box stays closed.
    lo, hi = ax_h.get_ylim()
    b0, b1 = max(y0, lo), min(y1, hi)
    ax_h.add_patch(
        Rectangle(
            (t0, b0),
            t1 - t0,
            b1 - b0,
            fill=False,
            ec="0.35",
            lw=0.6,
            zorder=4,
        )
    )
    return t0, t1, b0, b1


# Box-tag size on a history, in axes fractions (8 pt bold "(b.1)").
TAG_W_IN, TAG_H_IN = 0.38, 0.15


def place_box_tags(
    ax_h: plt.Axes,
    boxes: list[tuple[str, float, float, float, float]],
    t_m: np.ndarray,
    u_m: np.ndarray,
    t_n: np.ndarray,
    u_n: np.ndarray,
) -> None:
    """
    Tag each zoom box on a history ("(b.1)", "(b.2)") where it covers the least
    of the curves: candidates sit beside or above/below the box at several
    heights; the score counts curve samples under the tag, and overlap with any
    box or with a tag already placed is ruled out.

    boxes: (tag, t0, t1, b0, b1) in data units.
    """
    (xa, xb), (ya, yb) = ax_h.get_xlim(), ax_h.get_ylim()
    fx = lambda t: (t - xa) / (xb - xa)  # noqa: E731
    fy = lambda u: (u - ya) / (yb - ya)  # noqa: E731
    w = TAG_W_IN / HIST_W
    h = TAG_H_IN / HIST_H
    # Curves resampled densely so steep, fast traces count as occupied everywhere.
    tf = np.linspace(xa, xb, 30000)
    pts = [
        np.column_stack([fx(tf), fy(1000.0 * np.interp(tf, t, u))])
        for t, u in ((t_m, u_m), (t_n, u_n))
    ]
    xy = np.vstack(pts)
    rects_box = [(fx(t0), fx(t1), fy(b0), fy(b1)) for _tag, t0, t1, b0, b1 in boxes]
    # The panel title sits inside the axes at the top left.
    rects_box.append((0.0, TITLE_W_IN / HIST_W, 1.0 - TITLE_H_IN / HIST_H, 1.0))
    placed: list[tuple[float, float, float, float]] = []
    gap = 0.01

    def overlaps(r, s) -> bool:
        return not (r[1] <= s[0] or s[1] <= r[0] or r[3] <= s[2] or s[3] <= r[2])

    for tag, t0, t1, b0, b1 in boxes:
        x0, x1 = fx(t0), fx(t1)
        y0b, y1b = fy(b0), fy(b1)
        xc = 0.5 * (x0 + x1)
        # A grid of spots around the box (within ~0.3 of the axes width); the
        # emptiest wins and gets a thin leader when it is not next to the box.
        # Tag always starts at the box's horizontal centre (extends to the right),
        # just above or below the box; only its height is chosen.
        cands = [(xc, xc + w, bt, bt + h) for bt in np.arange(0.01, 0.99 - h, 0.01)]
        best, best_score = None, None
        right = [(x1 + gap, x1 + gap + w, bt, bt + h) for bt in np.arange(0.01, 0.99 - h, 0.01)]
        # First above/below from the centre; if nothing fits, beside on the right.
        for pool in (cands, right):
            if best is not None:
                break
            for r in pool:
                if r[0] < 0.0 or r[1] > 1.0 or r[2] < 0.0 or r[3] > 1.0:
                    continue
                if any(overlaps(r, s) for s in rects_box + placed):
                    continue
                p = 0.005
                inside = (
                    (xy[:, 0] >= r[0] - p) & (xy[:, 0] <= r[1] + p)
                    & (xy[:, 1] >= r[2] - p) & (xy[:, 1] <= r[3] + p)
                )
                # Gap from the tag to the box, in axes fractions.
                dx = max(x0 - r[1], 0.0, r[0] - x1)
                dy = max(y0b - r[3], 0.0, r[2] - y1b)
                gap_own = float(np.hypot(dx, dy))
                # Only spots touching the own box, so each tag reads as its box's.
                if gap_own > 0.02:
                    continue
                # A tag must read as belonging to its own box: never closer to another one.
                closer_other = False
                for _t, s0, s1, c0, c1 in boxes:
                    if (s0, s1) == (t0, t1):
                        continue
                    ox = max(fx(s0) - r[1], 0.0, r[0] - fx(s1))
                    oy = max(fy(c0) - r[3], 0.0, r[2] - fy(c1))
                    if np.hypot(ox, oy) <= gap_own + 0.02:
                        closer_other = True
                if closer_other:
                    continue
                # Covered trace costs most; then distance from the own box.
                score = 50.0 * int(inside.sum()) + 400.0 * gap_own
                if best_score is None or score < best_score:
                    best, best_score = r, score
        if best is None:
            continue
        placed.append(best)
        # Leader to the nearest box edge when the tag is not touching the box.
        bx = min(max(0.5 * (best[0] + best[1]), x0), x1)
        by = min(max(0.5 * (best[2] + best[3]), y0b), y1b)
        tx = min(max(bx, best[0]), best[1])
        ty = min(max(by, best[2]), best[3])
        if np.hypot((tx - bx) * HIST_W, (ty - by) * HIST_H) > 0.06:
            ax_h.plot([tx, bx], [ty, by], transform=ax_h.transAxes, color=STATION_GRAY, lw=0.4, zorder=5)
        # Reference to a zoom panel: gray, regular weight (the zoom's own tag is bold black).
        ax_h.text(
            best[0],
            best[2],
            tag,
            transform=ax_h.transAxes,
            ha="left",
            va="bottom",
            fontsize=ZOOM_TICK_FS,
            color=STATION_GRAY,
            zorder=6,
        )


def nrmse(
    t_m: np.ndarray, x_m: np.ndarray, t_n: np.ndarray, x_n: np.ndarray, t_lo: float, t_hi: float
) -> float:
    """
    NRMSE of MKR-α against Newmark (reference), in %, over [t_lo, t_hi]:
    Newmark interpolated onto the MKR record times, RMS of the difference
    divided by the range (max - min) of the reference over the same times.
    """
    m = (t_m >= t_lo) & (t_m <= min(t_hi, float(t_n[-1])))
    ref = np.interp(t_m[m], t_n, x_n)
    rms = float(np.sqrt(np.mean((x_m[m] - ref) ** 2)))
    return 100.0 * rms / float(ref.max() - ref.min())


def place_note(ax: plt.Axes, text: str, curves: list[tuple[np.ndarray, np.ndarray]]) -> None:
    """
    Write ``text`` in the axes corner (top/bottom, right/left) that overlaps the
    fewest curve samples and no other text already on the axes.
    """
    fig = ax.figure
    rend = fig.canvas.get_renderer()
    others = [t.get_window_extent(rend) for t in ax.texts if t.get_visible()]
    if ax.title.get_text():
        others.append(ax.title.get_window_extent(rend))
    pts = np.vstack([ax.transData.transform(np.column_stack([x, y])) for x, y in curves])
    best, best_score = None, None
    # Corners first (preferred on ties), then a grid over the whole axes.
    cands = [
        (0.99, 0.97, "right", "top"),
        (0.99, 0.03, "right", "bottom"),
        (0.01, 0.03, "left", "bottom"),
        (0.01, 0.97, "left", "top"),
    ]
    cands += [
        (gx, gy, "left", "bottom")
        for gy in np.linspace(0.02, 0.85, 12)
        for gx in np.linspace(0.01, 0.70, 16)
    ]
    for (cx, cy, ha, va) in cands:
        tx = ax.text(cx, cy, text, transform=ax.transAxes, ha=ha, va=va, fontsize=ZOOM_TICK_FS)
        bb = tx.get_window_extent(rend).expanded(1.06, 1.15)
        tx.remove()
        ab = ax.get_window_extent(rend)
        if bb.x0 < ab.x0 or bb.x1 > ab.x1 or bb.y0 < ab.y0 or bb.y1 > ab.y1:
            continue
        if any(bb.overlaps(o) for o in others):
            continue
        inside = (
            (pts[:, 0] >= bb.x0) & (pts[:, 0] <= bb.x1) & (pts[:, 1] >= bb.y0) & (pts[:, 1] <= bb.y1)
        )
        score = int(inside.sum())
        if best_score is None or score < best_score:
            best, best_score = (cx, cy, ha, va), score
    if best is None:
        return
    cx, cy, ha, va = best
    ax.text(
        cx, cy, text, transform=ax.transAxes, ha=ha, va=va, fontsize=ZOOM_TICK_FS, zorder=8,
        bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none", alpha=0.85),
    )


def load_base_ux(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """soil_base_primary.out: time, UX, UY. UX in m, as recorded."""
    arr = np.loadtxt(path)
    return arr[:, 0], arr[:, 1]


def relative_to(t: np.ndarray, u: np.ndarray, t_b: np.ndarray, u_b: np.ndarray) -> np.ndarray:
    """u minus the base UX, interpolated onto t. u is (n,) or (n, k)."""
    ub = np.interp(t, t_b, u_b)
    return u - (ub[:, None] if u.ndim == 2 else ub)


def main() -> int:
    for relative in (False, True):
        render(relative)
    return 0


def render(relative: bool) -> None:
    """Draw one version: absolute u_x, or Δu_x relative to the base node."""
    mkr = ROOT / MKR_DUMP
    nm = ROOT / NM_DUMP
    for eq in (mkr, nm):
        if not (eq / "window_meta.txt").is_file():
            raise SystemExit(f"missing dump: {eq}")

    prs.configure_font()
    xy = read_xy(mkr / "window_nodes.txt")
    eles = read_eles(mkr / "window_eles.txt")
    layers = read_layers(mkr / "window_quads.txt")
    order_m = read_order(mkr / "disp_nodes.txt")
    order_n = read_order(nm / "disp_nodes.txt")
    soil = soil_tags(xy)

    t_pm, u_pm = load_pier_ux(mkr / "pier_node_5.out")
    t_pn, u_pn = load_pier_ux(nm / "pier_node_5.out")
    t_pile_m, u_pile_m = series_for(mkr, order_m, list(PILE_TAGS))
    t_pile_n, u_pile_n = series_for(nm, order_n, list(PILE_TAGS))
    t_soil_m, u_soil_m = series_for(mkr, order_m, soil)
    t_soil_n, u_soil_n = series_for(nm, order_n, soil)
    if relative:
        t_bm, u_bm = load_base_ux(mkr / "soil_base_primary.out")
        t_bn, u_bn = load_base_ux(nm / "soil_base_primary.out")
        u_pm = relative_to(t_pm, u_pm, t_bm, u_bm)
        u_pn = relative_to(t_pn, u_pn, t_bn, u_bn)
        u_pile_m = relative_to(t_pile_m, u_pile_m, t_bm, u_bm)
        u_pile_n = relative_to(t_pile_n, u_pile_n, t_bn, u_bn)
        u_soil_m = relative_to(t_soil_m, u_soil_m, t_bm, u_bm)
        u_soil_n = relative_to(t_soil_n, u_soil_n, t_bn, u_bn)

    # One pair per station. Letters match the marks on (a).
    # Station names only; where each one is shows on (a).
    rows: list[tuple[str, str, str, np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = [
        ("pier", "b", "Pier top", t_pm, u_pm, t_pn, u_pn),
    ]
    for k in range(len(PILE_TAGS)):
        rows.append(
            (
                "pile_head" if (k == 0 and not relative) else "pile",
                "cd"[k],
                "Center pile",
                t_pile_m,
                u_pile_m[:, k],
                t_pile_n,
                u_pile_n[:, k],
            )
        )
    for k in range(len(soil)):
        rows.append(
            (
                "soil",
                "e"[k],
                "Soil",
                t_soil_m,
                u_soil_m[:, k],
                t_soil_n,
                u_soil_n[:, k],
            )
        )
    marks = [(PIER_TAG, "b")]
    marks += [(tag, letter) for tag, letter in zip(PILE_TAGS, "cd")]
    marks += [(tag, letter) for tag, letter in zip(soil, "e")]

    fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=150)
    ax_a = fig.add_axes(
        (ML / FIG_W, (FIG_H - MT - MESH_H) / FIG_H, MESH_W / FIG_W, MESH_H / FIG_H)
    )
    y_top = FIG_H - MT - HIST_H
    hist_axes: list[plt.Axes] = []
    zoom_axes: list[plt.Axes] = []
    zoom2_axes: list[plt.Axes] = []
    for i in range(N_HIST):
        bottom = y_top - i * (HIST_H + HIST_GAP)
        hist_axes.append(
            fig.add_axes((HIST_X / FIG_W, bottom / FIG_H, HIST_W / FIG_W, HIST_H / FIG_H))
        )
        zoom_axes.append(
            fig.add_axes((ZOOM_X / FIG_W, bottom / FIG_H, ZOOM_W / FIG_W, HIST_H / FIG_H))
        )
        zoom2_axes.append(
            fig.add_axes((ZOOM2_X / FIG_W, bottom / FIG_H, ZOOM_W / FIG_W, HIST_H / FIG_H))
        )

    draw_mesh(ax_a, xy, eles, layers, marks)

    # One limit per group (pier, pile, soil): round ceiling of the group peak,
    # used as the outermost tick so the axes end exactly on it.
    group_lim: dict[str, float] = {}
    for group, _l, _r, t_m, u_m, t_n, u_n in rows:
        peak = window_peak_mm(t_m, u_m, t_n, u_n)
        group_lim[group] = max(group_lim.get(group, 0.0), peak)
    group_lim = {g: round_lim(pk) for g, pk in group_lim.items()}

    # One window per zoom column, shared by all stations.
    w1, w2 = pick_windows(rows, group_lim)
    nrmse_log: list[str] = []
    zoom_log = [f"column 1 {w1:.0f}-{w1 + ZOOM_SPAN:.0f} s, column 2 {w2:.0f}-{w2 + ZOOM_SPAN:.0f} s"]
    for i, (ax, ax_z, ax_z2, (group, letter, rest, t_m, u_m, t_n, u_n)) in enumerate(
        zip(hist_axes, zoom_axes, zoom2_axes, rows)
    ):
        last = i == len(rows) - 1
        plot_pair(ax, t_m, u_m, t_n, u_n)
        style_history(ax)
        ax.set_xlim(T_MIN, T_MAX)
        lim = group_lim[group]
        ax.set_ylim(-lim, lim)
        ax.set_yticks([-lim, 0.0, lim])
        # (x.1): strongest-response window; (x.2): largest-difference window.
        tag1, tag2 = f"({letter}.1)", f"({letter}.2)"
        box1 = draw_zoom(ax, ax_z, w1, t_m, u_m, t_n, u_n, tag1, last)
        box2 = draw_zoom(ax, ax_z2, w2, t_m, u_m, t_n, u_n, tag2, last)
        place_box_tags(ax, [(tag1, *box1), (tag2, *box2)], t_m, u_m, t_n, u_n)
        panel_title(ax, letter, rest)
        # NRMSE of MKR-α vs Newmark over the plotted window.
        e = nrmse(t_m, u_m, t_n, u_n, T_MIN, T_MAX)
        nrmse_log.append(f"({letter}) NRMSE {e:.3f} %")
        # Lower-left corner of every history (quiet early record).
        ax.text(
            0.012, 0.04, rf"NRMSE $= {e:.1f}\%$", transform=ax.transAxes, ha="left", va="bottom",
            fontsize=ZOOM_TICK_FS, zorder=8,
            bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none", alpha=0.85),
        )
        ax.set_xticks(np.arange(T_MIN, T_MAX + 1.0, 30.0))
    for ax in hist_axes[:-1]:
        ax.tick_params(labelbottom=False)
    hist_axes[-1].set_xlabel(r"Time, $t$ (s)")
    # One label for the whole stack, centered on it, left of the tick labels.
    ylab = (
        r"Relative horizontal displacement, $\Delta u_x$ (mm)"
        if relative
        else r"Horizontal displacement, $u_x$ (mm)"
    )
    fig.text(
        (HIST_X - 0.47) / FIG_W,
        (MB + 0.5 * STACK_H) / FIG_H,
        ylab,
        rotation=90,
        ha="center",
        va="center",
    )

    integ = [
        Line2D(
            [0],
            [0],
            color=NM_COLOR,
            lw=NM_LW,
            ls="-",
            label="Newmark\n" + r"$\gamma = 0.5$" + "\n" + r"$\beta = 0.25$",
        ),
        Line2D(
            [0],
            [0],
            color="black",
            lw=MKR_LW,
            ls=MKR_LS,
            label=r"MKR-$\alpha$" + "\n" + r"$\rho_\infty^{2} = 0.5$",
        ),
    ]
    # Legend under the mesh. Both runs share one dt, so it is stated once (title).
    dt_nm, dt_mkr = read_dt(nm), read_dt(mkr)
    if abs(dt_nm - dt_mkr) > 1e-12:
        raise SystemExit(f"legend assumes one dt: Newmark {dt_nm}, MKR {dt_mkr}")
    leg = fig.legend(
        handles=integ,
        title=rf"$\Delta t = {1e3 * dt_mkr:.2f}$ ms",
        loc="upper left",
        bbox_to_anchor=(0.0, (FIG_H - MT - MESH_H - LEG_DROP) / FIG_H),
        bbox_transform=fig.transFigure,
        handlelength=1.0,
        handletextpad=0.3,
        labelspacing=1.1,
        frameon=False,
        borderaxespad=0.0,
        alignment="left",
    )
    leg.get_title().set_fontsize(plt.rcParams["legend.fontsize"])

    out = OUT.with_name(OUT.name + "_rel") if relative else OUT
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".png"), dpi=300)
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".svg"))
    plt.close(fig)

    print(f"wrote {out.with_suffix('.png')}")
    print(f"wrote {out.with_suffix('.pdf')}")
    print(
        "stations pier",
        PIER_TAG,
        "pile",
        [(tag, xy[tag]) for tag in PILE_TAGS],
        "soil",
        [(tag, xy[tag]) for tag in soil],
    )
    print(f"figure {FIG_W:.2f} x {FIG_H:.2f} in")
    print(f"MKR t_end {t_pm[-1]:.2f} s   Newmark t_end {t_pn[-1]:.2f} s")
    print(f"pier u_x(0) {1000 * u_pm[0]:.3f} mm")
    print("y limits (mm)", group_lim)
    for line in nrmse_log:
        print(" ", line)
    for line in zoom_log:
        print("  zoom", line)


if __name__ == "__main__":
    sys.exit(main())
