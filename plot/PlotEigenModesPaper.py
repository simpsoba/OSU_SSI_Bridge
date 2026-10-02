#!/usr/bin/env python3
"""
Goals
-----
Publication composite of the case-study eigenmodes (gravity + pier-base hold,
lumped plasticity, Baseline mesh). Same page language as
PlotModelSketchPaper.py: 6 in wide, 9 pt New Computer Modern Sans, stiff navy /
pier orange. Cubic Hermite frames from PlotEigenModes.

Default layout is full_grid: near-field zooms (a)--(e) for modes 1, 2, 4, 5, 25
in one row, 6 in wide. Shared window, equal max |u| scaling, BC pin triangles.

  python plot/PlotEigenModesPaper.py
  python plot/PlotEigenModesPaper.py --layout full_grid|column|stack|full

Default in:
  plot/out/eigen/compare_lumped_hold/eigen_modes.json
Default out:
  modal_analysis/figures/case_study_eigen_modes.{png,pdf,svg}
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.patches import Polygon

from PlotEigenModes import (
    NEP_HERMITE,
    ZOOM_GROUPS,
    deformed,
    domain_ylim,
    frame_segs,
    hermite_beam_xy,
    line_segs,
    load,
    node_xy,
    phi_maps,
    scale_for_mode,
    soil_polys,
    struct_node_xy,
    structure_xlim,
)
from PlotModelSketch import SPRING_GAP, draw_rot_spiral

# ------------------------------------------------------------
# 1. PAGE, FONT, COLORS (match PlotModelSketchPaper)
# ------------------------------------------------------------

FIG_W = 6.0  # in
# Zoom-only full_grid width (journal single-column budget)
FIG_W_FULL_GRID = 6.0
# FIG_H is chosen in build_figure from data aspect ratios (equal-aspect panels).
FONT_SIZE = 9
FONT_NAME = "NewComputerModernSans10"
FONT_DIR = Path.home() / "AppData" / "Local" / "Microsoft" / "Windows" / "Fonts"
FONT_FILES = (
    "NewCMSans10-Regular.otf",
    "NewCMSans10-Oblique.otf",
    "NewCMSans10-Bold.otf",
)

STIFF_COLOR = "#0b3c5d"
PIER_COLOR = "#cd6416"
PILE_COLOR = "#8B5A2B"
UNDEFORMED = "#b0b0b0"
UNDEFORMED_SOIL_FACE = "#e0e0e0"
UNDEFORMED_SOIL_EDGE = "#bdbdbd"
SOIL_FACE_DEF = "#f0e6c4"
SOIL_EDGE_DEF = "#8a7a64"
NODE_COLOR = "#263238"
# Shin thick free-field columns. Mode plots use a darker fill (not //// hatch):
# hatch on deformed quads follows surface waves and looks noisy at this size.
FF_FACE = "#c4b08a"
FF_EDGE = "#5a4a36"

STRUCT_COLOR = {
    "pier": PIER_COLOR,
    "deck": STIFF_COLOR,
    "cap": STIFF_COLOR,
    "pile": PILE_COLOR,
    "spring": "#6a1b9a",
    "ssi_spring": "#1565c0",  # same as schematic p-y (PlotModelSketch SSI_STYLE)
    "other": "#616161",
}
DRAW_GROUPS = ("pier", "deck", "cap", "pile", "spring", "ssi_spring", "other")

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
DEFAULT_JSON = HERE / "out" / "eigen" / "compare_lumped_hold" / "eigen_modes.json"
DEFAULT_OUT_DIR = REPO / "modal_analysis" / "figures"

# (mode_id, panel letter, view) — letters follow mode number (1,2,4,5,25)
# view: "near" = pier/deck/piles window; "full" = soil domain
# Column layout: (a) | (b)/(e) | (c) | (d)
PANELS_COLUMN = (
    (1, "a", "near"),
    (2, "b", "full"),
    (25, "e", "full"),
    (4, "c", "near"),
    (5, "d", "near"),
)
# Stack layout uses the same panel list as column (box order differs)
PANELS_STACK = PANELS_COLUMN
# All five modes, full soil domain — stacked review figure
PANELS_FULL = (
    (1, "a", "full"),
    (2, "b", "full"),
    (4, "c", "full"),
    (5, "d", "full"),
    (25, "e", "full"),
)
# Zoom-only row: (a,b,c) | (d,e). Letters follow mode number.
# Order matches _panel_boxes_full_grid.
PANELS_FULL_GRID = (
    (1, "a", "near"),
    (2, "b", "near"),
    (4, "c", "near"),
    (5, "d", "near"),
    (25, "e", "near"),
)
LAYOUTS = ("stack", "column", "full", "full_grid")
DEFAULT_LAYOUT = "full_grid"
PANELS_BY_LAYOUT = {
    "stack": PANELS_STACK,
    "column": PANELS_COLUMN,
    "full": PANELS_FULL,
    "full_grid": PANELS_FULL_GRID,
}
# Full-domain panels are small — exaggerate φ so soil waves still read
FULL_SF_BOOST = 2.5
# Near panels (a,c,d): tighter than PlotEigenModes.ZOOM_X_PAD (3.2 m)
NEAR_X_PAD = 1.6
# Zoom panels (full_grid): max |u| in the shared window = this × min(dx, dy)
NEAR_EQUAL_DISP_FRAC = 0.10
# Eigen BC glyphs: simple pin triangles (holdPier / soil-base fix UX+UY)
BC_PIN_S = 0.30  # m, pier-base triangle half-width
BC_SOIL_PIN_S = 0.20  # m, soil-base triangles
BC_FACE = "#ffffff"
BC_EDGE = "#9e9e9e"


def configure_font() -> None:
    """Register New CM Sans and set 9 pt text (same as elevation_paper)."""
    for name in FONT_FILES:
        path = FONT_DIR / name
        if path.is_file():
            font_manager.fontManager.addfont(str(path))
    plt.rcParams.update({
        "font.family": FONT_NAME,
        "font.size": FONT_SIZE,
        "axes.titlesize": FONT_SIZE,
        "axes.labelsize": FONT_SIZE,
        "xtick.labelsize": FONT_SIZE,
        "ytick.labelsize": FONT_SIZE,
        "mathtext.fontset": "custom",
        "mathtext.rm": FONT_NAME,
        "mathtext.it": f"{FONT_NAME}:oblique",
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "pdf.fonttype": 42,
    })


# ------------------------------------------------------------
# 2. FREE-FIELD HATCH + SPRING GLYPH + ONE PANEL
# ------------------------------------------------------------


def infer_L_half(
    xy: dict[int, tuple[float, float]],
    quads: list,
) -> float:
    """
    Near-field outer face |x| = L_half (Shin FF starts beyond this).

    Prefers sizes.L_half when present; else the largest soil |x| station
    strictly inside the mesh outer face.

    Args:    xy, quads
    Returns: L_half (m), or 0 if unknown
    """
    xs: list[float] = []
    for q in quads:
        for n in q:
            n = int(n)
            if n in xy:
                xs.append(abs(xy[n][0]))
    if not xs:
        return 0.0
    x_mesh = max(xs)
    uniq = sorted({round(x, 3) for x in xs})
    inner = [x for x in uniq if x < x_mesh - 0.5]
    return float(inner[-1]) if inner else 0.0


def soil_polys_nf_ff(
    quads: list,
    xy0: dict[int, tuple[float, float]],
    xy_draw: dict[int, tuple[float, float]],
    L_half: float,
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """
    Deformed soil polygons split by undeformed centre |x| ≷ L_half.

    Args:    quads, xy0  for FF test, xy_draw  plot coords, L_half (m)
    Returns: (nf_polys, ff_polys)
    """
    nf: list[np.ndarray] = []
    ff: list[np.ndarray] = []
    for q in quads:
        tags = [int(n) for n in q]
        if any(t not in xy0 or t not in xy_draw for t in tags):
            continue
        pts0 = np.array([xy0[t] for t in tags])
        xc = 0.5 * (float(pts0[:, 0].min()) + float(pts0[:, 0].max()))
        pts = np.array([xy_draw[t] for t in tags])
        if L_half > 0.0 and abs(xc) > L_half - 1.0e-6:
            ff.append(pts)
        else:
            nf.append(pts)
    return nf, ff


def draw_pier_rot_springs(
    ax,
    eles: list,
    xy1: dict[int, tuple[float, float]],
    color: str,
) -> None:
    """
    Rotational ZLS hinges as spirals at the deformed pier ends.

    Nodes are coincident (equalDOF UX/UY); spread them by SPRING_GAP like
    PlotModelSketch so the glyph reads (base 1--2, top 4--5).

    Args:    ax, eles, xy1  deformed coords, color
    Returns: none
    """
    for _e, ni, nj, grp in eles:
        if grp != "spring":
            continue
        ni, nj = int(ni), int(nj)
        if ni not in xy1 or nj not in xy1:
            continue
        tags = {ni, nj}
        if tags == {1, 2}:
            x, y = xy1[1]
            draw_rot_spiral(ax, x, y, x, y + SPRING_GAP, color)
        elif tags == {4, 5}:
            x, y = xy1[5]
            draw_rot_spiral(ax, x, y - SPRING_GAP, x, y, color)
        else:
            x, y = xy1[ni]
            draw_rot_spiral(ax, x, y, x, y + SPRING_GAP, color)


def draw_deformed_spring_coil(
    ax,
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    color: str,
) -> None:
    """
    Zigzag spring along a deformed chord (readable even when short).

    Args:    ax, endpoints (m), color
    Returns: none
    """
    dx = x1 - x0
    dy = y1 - y0
    L = float(np.hypot(dx, dy))
    if L < 1.0e-3:
        return
    nx, ny = -dy / L, dx / L
    # Larger amp than PlotModelSketch so short SSI gaps still read as coils
    n_zig = 5 if L < 0.35 else 7
    amp = min(0.22, max(0.06, 0.35 * L))
    xs = [x0]
    ys = [y0]
    for i in range(1, n_zig + 1):
        t = i / (n_zig + 1)
        side = 1.0 if (i % 2) else -1.0
        xs.append(x0 + t * dx + side * amp * nx)
        ys.append(y0 + t * dy + side * amp * ny)
    xs.append(x1)
    ys.append(y1)
    ax.plot(xs, ys, color=color, lw=0.45, zorder=5, solid_capstyle="round")


def draw_bc_triangle(
    ax,
    x: float,
    y: float,
    s: float,
    *,
    pointing: str = "down",
) -> None:
    """
    Simple pin triangle with tip at the node (no ground hatch).

    Args:    ax, x, y  tip (m), s  half-base (m),
             pointing  "down" | "left" | "right"
    Returns: none
    """
    h = 1.15 * s
    if pointing == "down":
        pts = [(x, y), (x - s, y - h), (x + s, y - h)]
    elif pointing == "left":
        pts = [(x, y), (x - h, y - s), (x - h, y + s)]
    elif pointing == "right":
        pts = [(x, y), (x + h, y - s), (x + h, y + s)]
    else:
        raise ValueError(f"pointing must be down|left|right, got {pointing!r}")
    ax.add_patch(
        Polygon(
            pts,
            closed=True,
            facecolor=BC_FACE,
            edgecolor=BC_EDGE,
            linewidth=0.2,  # same as deformed soil quad edges
            zorder=8,
        )
    )


def draw_eigen_bcs(
    ax,
    xy0: dict[int, tuple[float, float]],
    eles: list,
    quads: list,
    xlim: tuple[float, float],
) -> None:
    """
    Gravity + holdPier BC glyphs at undeformed nodes (eigen stage).

    Pier-base ZLS: both nodes pinned (UX+UY). Lower tip points down; upper
    (top of the rotational spring) points left so it clears the spiral.
    Soil base: downward triangles (fix 1 1 before Lysmer).

    Args:    ax, xy0, eles, quads, xlim  (clip soil pins to the window)
    Returns: none
    """
    # Pier-base hold: both ends of the base ZLS near y = 0.
    # Nodes are coincident in the mesh — spread like the spiral glyph.
    for _e, ni, nj, grp in eles:
        if grp != "spring":
            continue
        ni, nj = int(ni), int(nj)
        if ni not in xy0 or nj not in xy0:
            continue
        yi, yj = xy0[ni][1], xy0[nj][1]
        if abs(yi) > 1.0 and abs(yj) > 1.0:
            continue
        # Prefer cap TC (1) as the lower tip; inner (2) above
        if {ni, nj} == {1, 2}:
            x, y = xy0[1]
            draw_bc_triangle(ax, x, y, BC_PIN_S, pointing="down")
            draw_bc_triangle(ax, x, y + SPRING_GAP, BC_PIN_S, pointing="left")
        else:
            x, y = xy0[ni]
            draw_bc_triangle(ax, x, y, BC_PIN_S, pointing="down")
            draw_bc_triangle(ax, x, y + SPRING_GAP, BC_PIN_S, pointing="left")

    # Soil base nodes (lowest soil elevation)
    ys = [
        xy0[int(n)][1]
        for q in quads
        for n in q
        if int(n) in xy0
    ]
    if not ys:
        return
    yb = min(ys)
    xs = sorted({
        round(xy0[int(n)][0], 6)
        for q in quads
        for n in q
        if int(n) in xy0 and abs(xy0[int(n)][1] - yb) < 1.0e-6
    })
    xlo, xhi = xlim
    xs_win = [x for x in xs if xlo - 0.5 <= x <= xhi + 0.5]
    if not xs_win:
        return
    step = max(1, len(xs_win) // 10)
    for x in xs_win[::step]:
        draw_bc_triangle(ax, x, yb, BC_SOIL_PIN_S, pointing="down")


def draw_mode_panel(
    ax,
    *,
    xy0: dict[int, tuple[float, float]],
    xy1: dict[int, tuple[float, float]],
    eles: list,
    quads: list,
    phi: dict[int, tuple[float, float, float]],
    sf: float,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    letter: str,
    show_rot_springs: bool = True,
    show_ssi_springs: bool = True,
    show_bcs: bool = False,
    letter_loc: str = "upper left",
    L_half: float = 0.0,
    hatch_ff: bool = False,
) -> None:
    """
    Draw one mode: undeformed ghost, deformed Hermite frames + spring glyphs.

    Args:    ax, geometry, mode fields, window, letter, show_rot_springs,
             show_ssi_springs, show_bcs, letter_loc  "upper left" | "lower left",
             L_half, hatch_ff  Shin thick FF hatch (full-domain panels)
    Returns: none
    """
    # Undeformed soil ghost — near panels only (full-domain + large sf shears
    # leave empty wedges if a grey undeformed fill stays behind).
    if not hatch_ff:
        sp0 = soil_polys(quads, xy0)
        if sp0:
            ax.add_collection(
                PolyCollection(
                    sp0,
                    facecolors=UNDEFORMED_SOIL_FACE,
                    edgecolors=UNDEFORMED_SOIL_EDGE,
                    linewidths=0.15,
                    alpha=0.45,
                    zorder=0,
                )
            )
        sp1 = soil_polys(quads, xy1)
        if sp1:
            ax.add_collection(
                PolyCollection(
                    sp1,
                    facecolors=SOIL_FACE_DEF,
                    edgecolors=SOIL_EDGE_DEF,
                    linewidths=0.2,
                    alpha=0.40,
                    zorder=1,
                )
            )
    else:
        # Full domain: deformed NF + darker deformed FF. No interface line.
        # Keep soil translucent so the pier reads at this small scale.
        nf1, ff1 = soil_polys_nf_ff(quads, xy0, xy1, L_half)
        if nf1:
            ax.add_collection(
                PolyCollection(
                    nf1,
                    facecolors=SOIL_FACE_DEF,
                    edgecolors=SOIL_EDGE_DEF,
                    linewidths=0.2,
                    alpha=0.40,
                    zorder=1,
                )
            )
        if ff1:
            ax.add_collection(
                PolyCollection(
                    ff1,
                    facecolors=FF_FACE,
                    edgecolors=SOIL_EDGE_DEF,
                    linewidths=0.2,
                    alpha=0.50,
                    zorder=1.1,
                )
            )

    # Undeformed structure ghost
    if hatch_ff:
        # Full panels: pier + deck only (reference at rest; piles clutter the soil)
        for grp in ("pier", "deck"):
            segs = line_segs(eles, xy0, {grp})
            if segs:
                ax.add_collection(
                    LineCollection(
                        segs, colors=UNDEFORMED, linewidths=0.55, alpha=0.85, zorder=2
                    )
                )
    else:
        for grp in DRAW_GROUPS:
            if grp in ("spring", "ssi_spring"):
                continue
            segs = line_segs(eles, xy0, {grp})
            if segs:
                ax.add_collection(
                    LineCollection(
                        segs, colors=UNDEFORMED, linewidths=0.55, alpha=0.85, zorder=2
                    )
                )

    # Same pt widths on all panels (linewidth is paper units, not data units)
    pier_lw = 1.35
    frame_lw = 1.15
    draw_grps = ("pier", "deck", "cap", "pile", "other") if hatch_ff else DRAW_GROUPS
    for grp in draw_grps:
        if grp in ("spring", "ssi_spring"):
            continue
        segs = frame_segs(eles, xy0, phi, sf, {grp})
        if segs:
            ax.add_collection(
                LineCollection(
                    segs,
                    colors=STRUCT_COLOR.get(grp, "#333"),
                    linewidths=pier_lw if grp == "pier" else frame_lw,
                    zorder=3,
                )
            )

    # SSI translational springs (p-y / tip / soffit)
    if show_ssi_springs:
        ssi_color = STRUCT_COLOR["ssi_spring"]
        for seg in frame_segs(eles, xy0, phi, sf, {"ssi_spring"}):
            x0, y0 = float(seg[0, 0]), float(seg[0, 1])
            x1, y1 = float(seg[1, 0]), float(seg[1, 1])
            draw_deformed_spring_coil(ax, x0, y0, x1, y1, ssi_color)

    # Pier rotational ZLS (fiber hinges)
    if show_rot_springs:
        draw_pier_rot_springs(ax, eles, xy1, STRUCT_COLOR["spring"])

    # Node dots — near panels only (pile nodes look like a dotted centreline at full scale)
    if not hatch_ff:
        nx, ny = struct_node_xy(eles, xy1, ZOOM_GROUPS)
        if len(nx):
            ax.scatter(nx, ny, s=6, c=NODE_COLOR, alpha=0.45, linewidths=0, zorder=4)

    if show_bcs:
        draw_eigen_bcs(ax, xy0, eles, quads, xlim)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    # Limits padded to box aspect in build_figure so this does not letterbox.
    ax.set_aspect("equal", adjustable="box")
    # Mode shapes are arbitrarily scaled — no coordinate labels
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
        spine.set_color("0.55")

    # Panel tag only — mode / T / f go in the LaTeX caption
    if letter_loc == "lower left":
        tx, ty, va = 0.04, 0.04, "bottom"
    else:
        tx, ty, va = 0.03, 0.97, "top"
    ax.text(
        tx,
        ty,
        rf"({letter})",
        transform=ax.transAxes,
        ha="left",
        va=va,
        fontsize=FONT_SIZE,
        fontweight="bold",
        color="#222222",
        zorder=6,
    )


# ------------------------------------------------------------
# 3. FIGURE
# ------------------------------------------------------------


def pad_limits_to_aspect(
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    box_w: float,
    box_h: float,
) -> tuple[tuple[float, float], tuple[float, float]]:
    """
    Expand xlim or ylim so dy/dx equals the axes box aspect (h/w).

    Keeps equal-aspect plots from letterboxing inside their inch slot.

    Args:    xlim, ylim (m), box_w, box_h (in)
    Returns: (xlim, ylim)
    """
    x0, x1 = float(xlim[0]), float(xlim[1])
    y0, y1 = float(ylim[0]), float(ylim[1])
    dx = max(x1 - x0, 1.0e-6)
    dy = max(y1 - y0, 1.0e-6)
    target = box_h / max(box_w, 1.0e-6)
    if dy / dx < target:
        dy_new = target * dx
        mid = 0.5 * (y0 + y1)
        y0, y1 = mid - 0.5 * dy_new, mid + 0.5 * dy_new
    else:
        dx_new = dy / target
        mid = 0.5 * (x0 + x1)
        x0, x1 = mid - 0.5 * dx_new, mid + 0.5 * dx_new
    return (x0, x1), (y0, y1)


def window_unit_amp(
    xy0: dict[int, tuple[float, float]],
    phi: dict[int, tuple[float, float, float]],
    eles: list,
    quads: list,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
) -> float:
    """
    Max |u| = hypot(ux, uy) at unit scale for content inside the window.

    Includes Hermite pier/deck/cap/pile samples and soil nodes whose
    undeformed coordinates lie in xlim × ylim.

    Args:    xy0, phi, eles, quads, xlim, ylim (m)
    Returns: amplitude (m per unit sf)
    """
    xlo, xhi = xlim
    ylo, yhi = ylim
    max_a = 0.0

    def consider(x: float, y: float, ux: float, uy: float) -> None:
        nonlocal max_a
        if xlo <= x <= xhi and ylo <= y <= yhi:
            max_a = max(max_a, float(np.hypot(ux, uy)))

    for _e, ni, nj, grp in eles:
        if grp not in ("pier", "deck", "cap", "pile"):
            continue
        ni, nj = int(ni), int(nj)
        if ni not in xy0 or nj not in xy0 or ni not in phi or nj not in phi:
            continue
        xi, yi = xy0[ni]
        xj, yj = xy0[nj]
        uxi, uyi, rzi = phi[ni]
        uxj, uyj, rzj = phi[nj]
        pts = hermite_beam_xy(
            xi, yi, xj, yj, uxi, uyi, rzi, uxj, uyj, rzj, 1.0, NEP_HERMITE
        )
        x0 = np.linspace(xi, xj, NEP_HERMITE)
        y0 = np.linspace(yi, yj, NEP_HERMITE)
        for k in range(len(pts)):
            consider(
                float(x0[k]),
                float(y0[k]),
                float(pts[k, 0] - x0[k]),
                float(pts[k, 1] - y0[k]),
            )

    for q in quads:
        for n in q:
            n = int(n)
            if n not in xy0 or n not in phi:
                continue
            x, y = xy0[n]
            ux, uy, _rz = phi[n]
            consider(float(x), float(y), float(ux), float(uy))
    return max_a


def shared_near_window(
    xy0: dict[int, tuple[float, float]],
    phis: list,
    by_mode: dict[int, tuple],
    sf_by_mode: dict[int, float],
    eles: list,
    xz: tuple[float, float],
    near_mode_ids: list[int],
    box_w: float,
    box_h: float,
) -> tuple[tuple[float, float], tuple[float, float]]:
    """
    One near-field window for every zoom panel, centred on the undeformed
    structure so the undeformed pier sits at the same axes position.

    Expands just enough to fit deformed pier/deck/cap/pile of all near modes.

    Args:    geometry, mode fields, xz, near mode ids, panel size (in)
    Returns: (xlim, ylim)
    """
    ylim0 = domain_ylim(xy0, xy0)
    if ylim0 is None:
        ylim0 = (-25.0, 15.0)
    cx = 0.5 * (float(xz[0]) + float(xz[1]))
    cy = 0.5 * (float(ylim0[0]) + float(ylim0[1]))
    half_x = 0.5 * (float(xz[1]) - float(xz[0]))
    half_y = 0.5 * (float(ylim0[1]) - float(ylim0[0]))

    def grow(xy: dict[int, tuple[float, float]]) -> None:
        nonlocal half_x, half_y
        for grp in ("pier", "deck", "cap", "pile"):
            for seg in line_segs(eles, xy, {grp}):
                for px, py in seg:
                    half_x = max(half_x, abs(float(px) - cx))
                    half_y = max(half_y, abs(float(py) - cy))

    grow(xy0)
    for mid in near_mode_ids:
        if mid not in by_mode or mid not in sf_by_mode:
            continue
        i, _m = by_mode[mid]
        grow(deformed(xy0, phis[i], sf_by_mode[mid]))
        for seg in frame_segs(eles, xy0, phis[i], sf_by_mode[mid], {"pier", "deck", "cap", "pile"}):
            for px, py in seg:
                half_x = max(half_x, abs(float(px) - cx))
                half_y = max(half_y, abs(float(py) - cy))

    pad = 1.03
    xlim = (cx - pad * half_x, cx + pad * half_x)
    ylim = (cy - pad * half_y, cy + pad * half_y)
    return pad_limits_to_aspect(xlim, ylim, box_w, box_h)


def _panel_boxes_stack(
    fig_w: float, dy: float, dx_near: float, dx_full: float
) -> tuple[float, list[tuple[float, float, float, float]]]:
    """
    Layout: 1|4|5 on top; full-width 2; full-width 25.

    Args:    fig_w, dy, dx_near, dx_full  (in / m)
    Returns: fig_h (in), list of (x, y, w, h) inch boxes in PANELS order
    """
    ml = mr = 0.06
    mb = mt = 0.04
    gap_x = 0.10
    gap_y = 0.08

    w_near = (fig_w - ml - mr - 2.0 * gap_x) / 3.0
    h_near = w_near * (dy / dx_near)
    w_full = fig_w - ml - mr
    h_full = w_full * (dy / dx_full)
    fig_h = mb + 2.0 * h_full + gap_y + h_near + mt + gap_y

    y25 = mb
    y2 = y25 + h_full + gap_y
    y_near = y2 + h_full + gap_y
    # PANELS order: 1, 2, 25, 4, 5 — top row is 1|4|5
    boxes = [
        (ml + 0 * (w_near + gap_x), y_near, w_near, h_near),  # 1
        (ml, y2, w_full, h_full),  # 2
        (ml, y25, w_full, h_full),  # 25
        (ml + 1 * (w_near + gap_x), y_near, w_near, h_near),  # 4
        (ml + 2 * (w_near + gap_x), y_near, w_near, h_near),  # 5
    ]
    return fig_h, boxes


def _panel_boxes_column(
    fig_w: float, dy: float, dx_near: float, dx_full: float
) -> tuple[float, list[tuple[float, float, float, float]]]:
    """
    Layout: (a)1 | (b)2 above (e)25 | (c)4 | (d)5.

    Top of (b) aligns with top of (a); bottom of (e) with bottom of (a).
    w_full from that height + equal-aspect soil window.

    Args:    fig_w, dy, dx_near, dx_full  (in / m)
    Returns: fig_h (in), list of (x, y, w, h) inch boxes in PANELS order
    """
    ml = mr = 0.04
    mb = mt = 0.02
    # ~1/3 of the previous 0.08 / 0.06 gaps
    gap_x = 0.027
    gap_y = 0.02

    usable = fig_w - ml - mr - 3.0 * gap_x
    # Near width from sharing usable with one full column of equal-aspect soil
    # whose two stacked panels + gap_y fill h_near exactly.
    # h_near = w_near * dy/dx_near
    # 2*h_full + gap_y = h_near  →  h_full = (h_near - gap_y)/2
    # w_full = h_full * dx_full/dy
    # 3*w_near + w_full = usable
    aspect_n = dy / dx_near
    aspect_f = dy / dx_full
    # 3*w_near + ((w_near*aspect_n - gap_y)/2)* (dx_full/dy) = usable
    # 3*w_near + (w_near*aspect_n - gap_y)/(2*aspect_f) = usable
    coef = 3.0 + aspect_n / (2.0 * aspect_f)
    w_near = (usable + gap_y / (2.0 * aspect_f)) / coef
    h_near = w_near * aspect_n
    h_full = (h_near - gap_y) / 2.0
    w_full = h_full / aspect_f
    if w_near < 0.4 or h_full < 0.25:
        # fallback: clip and accept a short stack (should not hit at FIG_W=6)
        w_near = max(w_near, 0.4)
        h_near = w_near * aspect_n
        h_full = max((h_near - gap_y) / 2.0, 0.25)
        w_full = h_full / aspect_f
    fig_h = mb + h_near + mt

    x1 = ml
    x_col = x1 + w_near + gap_x
    x4 = x_col + w_full + gap_x
    x5 = x4 + w_near + gap_x
    y_near = mb
    # Flush: (e) on the bottom edge of (a), (b) on the top edge of (a)
    y25 = y_near
    y2 = y_near + h_near - h_full

    boxes = [
        (x1, y_near, w_near, h_near),  # 1 (a)
        (x_col, y2, w_full, h_full),  # 2 (b)
        (x_col, y25, w_full, h_full),  # 25 (e)
        (x4, y_near, w_near, h_near),  # 4 (c)
        (x5, y_near, w_near, h_near),  # 5 (d)
    ]
    return fig_h, boxes


def _panel_boxes_full(
    fig_w: float, dy: float, dx_full: float, n_panels: int = 5
) -> tuple[float, list[tuple[float, float, float, float]]]:
    """
    Layout: n full-domain panels stacked, mode order top → bottom.

    Args:    fig_w, dy, dx_full (in / m), n_panels
    Returns: fig_h (in), list of (x, y, w, h) inch boxes
    """
    ml = mr = 0.04
    mb = mt = 0.02
    gap_y = 0.03
    w = fig_w - ml - mr
    h = w * (dy / dx_full)
    fig_h = mb + mt + n_panels * h + (n_panels - 1) * gap_y
    boxes = []
    for i in range(n_panels):
        # top to bottom: panel 0 at top
        y = mb + (n_panels - 1 - i) * (h + gap_y)
        boxes.append((ml, y, w, h))
    return fig_h, boxes


def _panel_boxes_full_grid(
    fig_w: float, dy: float, dx_near: float, dx_full: float
) -> tuple[float, list[tuple[float, float, float, float]]]:
    """
    Zoom-only layout (PANELS_FULL_GRID order):

        (a,b,c) zooms | (d,e) zooms

    Five equal-width near panels (same height from equal aspect).

    Args:    fig_w, dy, dx_near, dx_full (in / m)
    Returns: fig_h (in), list of (x, y, w, h) inch boxes
    """
    _ = dx_full
    ml = mr = 0.04
    mb = mt = 0.02
    gap_x = 0.025
    gap_g = 0.04  # gap between the (a–c) and (d–e) groups
    n = 5
    aspect_n = dy / dx_near
    # One group gap + (n-1) panel gaps, but group gap replaces the mid panel gap
    usable = fig_w - ml - mr - (n - 2) * gap_x - gap_g
    w_near = usable / n
    h_near = w_near * aspect_n
    fig_h = mb + mt + h_near

    y = mb
    boxes: list[tuple[float, float, float, float]] = []
    x = ml
    for k in range(n):
        boxes.append((x, y, w_near, h_near))
        if k == 2:
            x += w_near + gap_g
        else:
            x += w_near + gap_x
    return fig_h, boxes



def build_figure(
    data: dict,
    out_stem: Path,
    layout: str = DEFAULT_LAYOUT,
    sf_mult: float = 1.0,
    show_bcs: bool = True,
) -> None:
    """
    Write one case-study eigenmode figure.

    Args:    data, out_stem without suffix,
             layout  "stack" | "column" | "full" | "full_grid",
             sf_mult  multiply every mode scale factor (1 = default),
             show_bcs  pier-base + soil-base pin glyphs (eigen / holdPier)
    Returns: none
    """
    if layout not in LAYOUTS:
        raise ValueError(f"layout must be one of {LAYOUTS}, got {layout!r}")

    panels = PANELS_BY_LAYOUT[layout]
    configure_font()
    xy0 = node_xy(data)
    phis = phi_maps(data)
    meta = data.get("modes_meta", [])
    eles = data.get("elements", [])
    quads = data.get("soil_quads", [])

    by_mode = {int(m["mode"]): (i, m) for i, m in enumerate(meta)}
    xz = structure_xlim(xy0, eles, pad=NEAR_X_PAD)
    if xz is None:
        xz = (-8.0, 8.0)

    xs_all = [p[0] for p in xy0.values()]
    x_full = (min(xs_all) - 2.0, max(xs_all) + 2.0)
    ys_all = [p[1] for p in xy0.values()]
    dy = 1.08 * (max(ys_all) - min(ys_all) + 1.0)
    dx_near = max(xz[1] - xz[0], 1.0)
    dx_full = max(x_full[1] - x_full[0], 1.0)
    sizes = data.get("sizes") or {}
    L_half = float(sizes["L_half"]) if sizes.get("L_half") else infer_L_half(xy0, quads)

    fig_w = FIG_W_FULL_GRID if layout == "full_grid" else FIG_W
    if layout == "stack":
        fig_h, panel_boxes = _panel_boxes_stack(fig_w, dy, dx_near, dx_full)
    elif layout == "full":
        fig_h, panel_boxes = _panel_boxes_full(fig_w, dy, dx_full, n_panels=len(panels))
    elif layout == "full_grid":
        fig_h, panel_boxes = _panel_boxes_full_grid(fig_w, dy, dx_near, dx_full)
    else:
        fig_h, panel_boxes = _panel_boxes_column(fig_w, dy, dx_near, dx_full)

    near_ids = [mid for mid, _let, view in panels if view == "near"]
    modes_with_full = {mid for mid, _let, view in panels if view == "full"}
    sf_by_mode: dict[int, float] = {}

    # Near panels: one shared window, sf so max |u| in that window matches
    shared_near: tuple[tuple[float, float], tuple[float, float]] | None = None
    near_box = next(
        ((w, h) for (_m, _l, v), (_x, _y, w, h) in zip(panels, panel_boxes) if v == "near"),
        None,
    )
    if near_ids and near_box is not None:
        w_n, h_n = near_box
        ylim0 = domain_ylim(xy0, xy0) or (-25.0, 15.0)
        xlim_p, ylim_p = pad_limits_to_aspect(xz, ylim0, w_n, h_n)
        amps = {}
        for mid in near_ids:
            if mid not in by_mode:
                continue
            i, _m = by_mode[mid]
            amps[mid] = window_unit_amp(
                xy0, phis[i], eles, quads, xlim_p, ylim_p
            )
        dx_w = xlim_p[1] - xlim_p[0]
        dy_w = ylim_p[1] - ylim_p[0]
        target_u = NEAR_EQUAL_DISP_FRAC * min(dx_w, dy_w) * float(sf_mult)
        for mid, amp in amps.items():
            sf_by_mode[mid] = (target_u / amp) if amp > 1.0e-30 else 0.0
        shared_near = shared_near_window(
            xy0,
            phis,
            by_mode,
            sf_by_mode,
            eles,
            xz,
            near_ids,
            w_n,
            h_n,
        )
        if show_bcs:
            # Room under the soil base for pin triangles + ground hatch
            x0, x1 = shared_near[0]
            y0, y1 = shared_near[1]
            drop = 2.2 * max(BC_PIN_S, BC_SOIL_PIN_S)
            shared_near = pad_limits_to_aspect((x0, x1), (y0 - drop, y1), w_n, h_n)
        # Re-fit sf on the final shared window so max |u| matches there
        xlim_f, ylim_f = shared_near
        dx_f = xlim_f[1] - xlim_f[0]
        dy_f = ylim_f[1] - ylim_f[0]
        target_u = NEAR_EQUAL_DISP_FRAC * min(dx_f, dy_f) * float(sf_mult)
        for mid in list(sf_by_mode):
            if mid not in by_mode:
                continue
            i, _m = by_mode[mid]
            amp = window_unit_amp(xy0, phis[i], eles, quads, xlim_f, ylim_f)
            sf_by_mode[mid] = (target_u / amp) if amp > 1.0e-30 else 0.0

    # Full-domain panels keep the usual scale_for_mode (+ boost)
    for mid, _let, _view in panels:
        if mid in sf_by_mode or mid not in by_mode:
            continue
        i, _m = by_mode[mid]
        sf_m, _amp, _H, _kind, _tgt = scale_for_mode(xy0, phis[i], eles)
        if mid in modes_with_full:
            sf_m *= FULL_SF_BOOST
        sf_by_mode[mid] = sf_m * float(sf_mult)

    fig = plt.figure(figsize=(fig_w, fig_h), dpi=300)

    def add_panel(x_in: float, y_in: float, w_in: float, h_in: float):
        return fig.add_axes(
            [x_in / fig_w, y_in / fig_h, w_in / fig_w, h_in / fig_h]
        )

    for (mode_id, letter, view), (x_in, y_in, w_in, h_in) in zip(
        panels, panel_boxes
    ):
        ax = add_panel(x_in, y_in, w_in, h_in)
        if mode_id not in by_mode:
            ax.text(0.5, 0.5, f"({letter}) missing", ha="center", va="center")
            ax.axis("off")
            continue
        i, _m = by_mode[mode_id]
        phi = phis[i]
        sf = sf_by_mode[mode_id]
        xy1 = deformed(xy0, phi, sf)
        if view == "near" and shared_near is not None:
            xlim, ylim = shared_near
        elif view == "near":
            ylim = domain_ylim(xy0, xy1)
            xlim = xz
            if ylim is None:
                ylim = (-25.0, 15.0)
            xlim, ylim = pad_limits_to_aspect(xlim, ylim, w_in, h_in)
        else:
            # X: deformed soil. Y: full height so pier/deck stay inside the axes
            xs = [xy1[int(n)][0] for q in quads for n in q if int(n) in xy1]
            if xs:
                pad_x = 0.02 * (max(xs) - min(xs) + 1.0)
                xlim = (min(xs) - pad_x, max(xs) + pad_x)
            else:
                xlim = x_full
            ylim = domain_ylim(xy0, xy1)
            if ylim is None:
                ylim = (-25.0, 15.0)
            xlim, ylim = pad_limits_to_aspect(xlim, ylim, w_in, h_in)
        # Springs on near panels and on the all-full review figures
        show_springs = (view == "near") or (layout in ("full", "full_grid"))
        draw_mode_panel(
            ax,
            xy0=xy0,
            xy1=xy1,
            eles=eles,
            quads=quads,
            phi=phi,
            sf=sf,
            xlim=xlim,
            ylim=ylim,
            letter=letter,
            show_rot_springs=show_springs,
            show_ssi_springs=show_springs,
            show_bcs=show_bcs,
            # Near panels: letter at bottom so it clears the deck
            letter_loc="lower left" if view == "near" else "upper left",
            L_half=L_half,
            hatch_ff=(view == "full"),
        )

    out_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".png", ".pdf", ".svg"):
        fig.savefig(out_stem.with_suffix(ext), dpi=300)
    plt.close(fig)
    print(
        f"PlotEigenModesPaper: wrote {out_stem}.png / .pdf / .svg  "
        f"({fig_w:.2f} x {fig_h:.2f} in, layout={layout}, sf_mult={sf_mult:g})"
    )


def main(argv: list[str] | None = None) -> int:
    """
    CLI: [--layout …] [--sf-mult f] [--no-bcs] [json] [out_stem].

    Args:    argv
    Returns: process status
    """
    args = list(sys.argv[1:] if argv is None else argv)
    layout = DEFAULT_LAYOUT
    sf_mult = 1.0
    show_bcs = True
    if "--layout" in args:
        i = args.index("--layout")
        if i + 1 >= len(args):
            print(
                "PlotEigenModesPaper: --layout needs stack|column|full|full_grid",
                file=sys.stderr,
            )
            return 1
        layout = args[i + 1]
        del args[i : i + 2]
    if "--sf-mult" in args:
        i = args.index("--sf-mult")
        if i + 1 >= len(args):
            print("PlotEigenModesPaper: --sf-mult needs a float", file=sys.stderr)
            return 1
        sf_mult = float(args[i + 1])
        del args[i : i + 2]
    if "--bcs" in args:
        args.remove("--bcs")
        show_bcs = True
    if "--no-bcs" in args:
        args.remove("--no-bcs")
        show_bcs = False
    json_path = Path(args[0]) if args else DEFAULT_JSON
    if not json_path.is_file():
        print(f"PlotEigenModesPaper: missing {json_path}", file=sys.stderr)
        return 1
    out_stem = (
        Path(args[1])
        if len(args) > 1
        else DEFAULT_OUT_DIR / "case_study_eigen_modes"
    )
    build_figure(
        load(json_path),
        out_stem,
        layout=layout,
        sf_mult=sf_mult,
        show_bcs=show_bcs,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
