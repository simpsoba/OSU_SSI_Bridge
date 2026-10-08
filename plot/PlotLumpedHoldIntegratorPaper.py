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
  plot/out/eq_offline/compare/lumped_hold_integrator_cudamkr.{png,pdf}
  plot/out/eq_offline/compare/lumped_hold_integrator_cudamkr_rel.{png,pdf}

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
# Serial dumps compared, both at dt = 0.00756 s and run to 180 s. Newmark: UmfPack
# (lumped_hold_newmark_dt10; the dt/2 runs lumped_hold_newmark[_t180] are kept).
# MKR-α: CudaMKRAlpha on CuDSS, run after the CudaMKRAlpha fix. The earlier figure from
# lumped_hold_mkr (MKRAlphaExplicitMultiSOE + UmfPack) stays as lumped_hold_integrator.*.
MKR_DUMP = "lumped_hold_cudamkr"
NM_DUMP = "lumped_hold_newmark_dt10"
OUT = ROOT / "compare" / "lumped_hold_integrator_cudamkr"

# Center pile: cap BC (head) and iy = 10 (mid-depth). See BuildPilesNodes.tcl.
PILE_TAGS = (1028, 2110)
PIER_TAG = 5
# Soil station off the shafts (piles sit at x = 0 and ±1.83 m).
SOIL_X = 5.7912
# Mudline only (the 10.1 m soil history tracks the 10.1 m pile one closely).
SOIL_Y = (0.0,)

# Newmark: wide light-grey solid line under MKR-α: thin black dashes, so the grey
# shows on both sides of the dash where the two agree (and still reads in B&W).
# Dash length is in points.
NM_COLOR = "#a8a8a8"
NM_LW = 2.0
MKR_LW = 0.8
MKR_LS = (0, (2.6, 1.8))

FIG_W = 6.0
# Mesh (true scale, no axes) is a left column, top-aligned with (b). The legend
# is one row across the top. Histories are stacked to the right and share one
# time axis; the zooms are a narrow column at the far right, with time labels
# only (their u range is the grey box on each history).
ML, MR, MB, MT = 0.08, 0.16, 0.40, 0.46
COL_GAP = 0.62
HIST_H = 0.62
HIST_GAP = 0.25
N_HIST = 4
STACK_H = N_HIST * HIST_H + (N_HIST - 1) * HIST_GAP
MESH_H = 2.20
ZOOM_GAP = 0.16
ZOOM_W = 1.00
# Zoom window width (s). Windows start on an odd whole second, so every zoom
# has its two time ticks (t0 + 1, t0 + 3) at the same place on its axis.
ZOOM_SPAN = 4.0
# window_nodes extent used for the mesh (m), matching set_xlim / set_ylim below.
MESH_DX = 17.10
MESH_DY = 35.40
MESH_W = MESH_H * MESH_DX / MESH_DY
HIST_X = ML + MESH_W + COL_GAP
HIST_W = FIG_W - HIST_X - ZOOM_GAP - ZOOM_W - MR
ZOOM_X = HIST_X + HIST_W + ZOOM_GAP
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
        if tag == PIER_TAG:
            # Below-right of the soffit so the letter clears the deck chord.
            xytext, ha = (8, -8), "left"
            textcoords = "offset points"
        elif 2000 <= tag < 2300 or tag in (1027, 1028, 1029):
            xytext, ha = (-4.6, y), "right"
            textcoords = "data"
        else:
            xytext, ha = (7.15, y), "left"
            textcoords = "data"
        ax.annotate(
            rf"\textbf{{({letter})}}",
            xy=(x, y),
            xytext=xytext,
            textcoords=textcoords,
            ha=ha,
            va="center",
            fontsize=7,
            color="#111111",
            zorder=6,
            clip_on=False,
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
    ax.set_title(r"\textbf{(a)}", loc="left", pad=5)


def coord_label(x_m: float, y_m: float) -> str:
    """'$(5.8, -1.0)$ m' = (x, y), with 0 for |v| < 0.05 (avoids a signed-zero '-0.0')."""

    def fmt(v: float) -> str:
        return "0" if abs(v) < 0.05 else f"{v:.1f}"

    return rf"$({fmt(x_m)}, {fmt(y_m)})$ m"


def panel_title(ax: plt.Axes, letter: str, rest: str) -> None:
    """Bold panel letter, normal-weight station text."""
    text = rf"\textbf{{({letter})}}" if not rest else rf"\textbf{{({letter})}}\quad {rest}"
    ax.set_title(text, loc="left", pad=5)


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


def zoom_window(
    t_m: np.ndarray, u_m: np.ndarray, t_n: np.ndarray, u_n: np.ndarray
) -> tuple[float, float, float, float]:
    """
    Zoom box around the peak |u| where both integrators exist.

    t0 is the odd whole second nearest t_pk - ZOOM_SPAN/2, so the peak falls in
    the middle half of the window.

    Returns: t0, t1 (s), y0, y1 (mm), y padded 10 % of the range.
    """
    t_end = min(T_MAX, float(t_n[-1]), float(t_m[-1]))
    inside = (t_m >= T_MIN) & (t_m <= t_end)
    k = int(np.argmax(np.abs(u_m[inside])))
    t_pk = float(t_m[inside][k])
    t0 = 2.0 * round((t_pk - 0.5 * ZOOM_SPAN - 1.0) / 2.0) + 1.0
    t1 = t0 + ZOOM_SPAN
    in_m = (t_m >= t0) & (t_m <= t1)
    in_n = (t_n >= t0) & (t_n <= t1)
    vals = np.concatenate([1000.0 * u_m[in_m], 1000.0 * u_n[in_n]])
    lo, hi = float(vals.min()), float(vals.max())
    pad = 0.10 * (hi - lo)
    return t0, t1, lo - pad, hi + pad


def draw_zoom(
    ax_h: plt.Axes,
    ax_z: plt.Axes,
    t_m: np.ndarray,
    u_m: np.ndarray,
    t_n: np.ndarray,
    u_n: np.ndarray,
) -> None:
    """Magnified peak window in ax_z; the same window outlined on the history ax_h."""
    t0, t1, y0, y1 = zoom_window(t_m, u_m, t_n, u_n)
    plot_pair(ax_z, t_m, u_m, t_n, u_n)
    ax_z.set_xlim(t0, t1)
    ax_z.set_ylim(y0, y1)
    ax_z.grid(True, which="major", color="0.85", lw=0.5, alpha=0.25)
    ax_z.set_axisbelow(True)
    # Full frame so the inset reads as a magnifier of the boxed window.
    for side in ("top", "right"):
        ax_z.spines[side].set_visible(True)
    ax_z.set_xticks([t0 + 1.0, t0 + 3.0])
    # No u labels on the zoom: its range is the grey box on the history.
    ax_z.yaxis.set_major_locator(MaxNLocator(nbins=3))
    ax_z.tick_params(axis="y", length=0, labelleft=False, labelright=False)
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


def label_zoom(
    ax_z: plt.Axes,
    text: str,
    t_m: np.ndarray,
    u_m: np.ndarray,
    t_n: np.ndarray,
    u_n: np.ndarray,
) -> None:
    """Bold label in the zoom corner with the fewest curve points (corner = 30 % x 30 %)."""
    (t0, t1), (y0, y1) = ax_z.get_xlim(), ax_z.get_ylim()
    xs, ys = [], []
    for t, u in ((t_m, u_m), (t_n, u_n)):
        m = (t >= t0) & (t <= t1)
        xs.append((t[m] - t0) / (t1 - t0))
        ys.append((1000.0 * u[m] - y0) / (y1 - y0))
    x, y = np.concatenate(xs), np.concatenate(ys)
    corners = {
        (0.04, 0.95, "left", "top"): (x < 0.30) & (y > 0.70),
        (0.96, 0.95, "right", "top"): (x > 0.70) & (y > 0.70),
        (0.04, 0.05, "left", "bottom"): (x < 0.30) & (y < 0.30),
        (0.96, 0.05, "right", "bottom"): (x > 0.70) & (y < 0.30),
    }
    (cx, cy, ha, va) = min(corners, key=lambda c: int(corners[c].sum()))
    ax_z.text(
        cx,
        cy,
        rf"\textbf{{{text}}}",
        transform=ax_z.transAxes,
        ha=ha,
        va=va,
        zorder=12,
        bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none", alpha=0.85),
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
    pile_where = [coord_label(*xy[tag]) for tag in PILE_TAGS]
    soil_where = [coord_label(*xy[tag]) for tag in soil]
    rows: list[tuple[str, str, str, np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = [
        ("pier", "b", f"Pier top, {coord_label(*xy[PIER_TAG])}", t_pm, u_pm, t_pn, u_pn),
    ]
    for k, where in enumerate(pile_where):
        rows.append(
            (
                "pile_head" if (k == 0 and not relative) else "pile",
                "cd"[k],
                f"Center pile, {where}",
                t_pile_m,
                u_pile_m[:, k],
                t_pile_n,
                u_pile_n[:, k],
            )
        )
    for k, where in enumerate(soil_where):
        rows.append(
            (
                "soil",
                "e"[k],
                f"Soil, {where}",
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
    for i in range(N_HIST):
        bottom = y_top - i * (HIST_H + HIST_GAP)
        hist_axes.append(
            fig.add_axes((HIST_X / FIG_W, bottom / FIG_H, HIST_W / FIG_W, HIST_H / FIG_H))
        )
        zoom_axes.append(
            fig.add_axes((ZOOM_X / FIG_W, bottom / FIG_H, ZOOM_W / FIG_W, HIST_H / FIG_H))
        )

    draw_mesh(ax_a, xy, eles, layers, marks)

    # One limit per group (pier, pile, soil): round ceiling of the group peak,
    # used as the outermost tick so the axes end exactly on it.
    group_lim: dict[str, float] = {}
    for group, _l, _r, t_m, u_m, t_n, u_n in rows:
        peak = window_peak_mm(t_m, u_m, t_n, u_n)
        group_lim[group] = max(group_lim.get(group, 0.0), peak)
    group_lim = {g: round_lim(pk) for g, pk in group_lim.items()}

    for ax, ax_z, (group, letter, rest, t_m, u_m, t_n, u_n) in zip(hist_axes, zoom_axes, rows):
        plot_pair(ax, t_m, u_m, t_n, u_n)
        style_history(ax)
        ax.set_xlim(T_MIN, T_MAX)
        lim = group_lim[group]
        ax.set_ylim(-lim, lim)
        ax.set_yticks([-lim, 0.0, lim])
        draw_zoom(ax, ax_z, t_m, u_m, t_n, u_n)
        panel_title(ax, letter, rest)
        # Zoom of the boxed window on the same row: (b.1), (c.1), ... Inside the
        # zoom (above it sit the time labels of the zoom above), in its emptiest corner.
        label_zoom(ax_z, f"({letter}.1)", t_m, u_m, t_n, u_n)
        ax.set_xticks(np.arange(T_MIN, T_MAX + 1.0, 30.0))
    for ax in hist_axes[:-1]:
        ax.tick_params(labelbottom=False)
    hist_axes[-1].set_xlabel(r"Time, $t$ (s)")
    # One label for the whole stack, centered on it, left of the tick labels.
    ylab = r"$\Delta u_x$ (mm)" if relative else r"$u_x$ (mm)"
    fig.text(
        (HIST_X - 0.50) / FIG_W,
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
            label=rf"Newmark, $\gamma = 0.5$, $\beta = 0.25$, $\Delta t = {1e3 * read_dt(nm):.2f}$ ms",
        ),
        Line2D(
            [0],
            [0],
            color="black",
            lw=MKR_LW,
            ls=MKR_LS,
            label=rf"MKR-$\alpha$, $\rho_\infty^{{2}} = 0.5$, $\Delta t = {1e3 * read_dt(mkr):.2f}$ ms",
        ),
    ]
    # One row across the top; a short line sample keeps each entry on one line.
    fig.legend(
        handles=integ,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.0 - 0.03 / FIG_H),
        bbox_transform=fig.transFigure,
        ncol=2,
        handlelength=1.6,
        handletextpad=0.5,
        columnspacing=1.6,
        frameon=False,
        borderaxespad=0.0,
    )

    out = OUT.with_name(OUT.name + "_rel") if relative else OUT
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".png"), dpi=300)
    fig.savefig(out.with_suffix(".pdf"))
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


if __name__ == "__main__":
    sys.exit(main())
