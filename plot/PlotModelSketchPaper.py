#!/usr/bin/env python3
"""
Goals
-----
Publication version of the model elevation (6 x 6 in, 9 pt,
New Computer Modern Sans 10). Reads the same DumpModelSketch.tcl JSON as
PlotModelSketch.py and writes PNG + PDF next to elevation.png.

  (a) full soil domain, full-width strip, true scale
  (b) near field: deck, pier, cap, piles, SSI springs
  (c) zoom on the pile cap and pier base hinge
  (d) inset: pier top + λ_L^3 f_p (same force arrow on a, b)
      legend under (c)

Water: light fill + dashed free surface (physical flume, not meshed);
arrows = rho_w g h_w surcharge on the y=0 soil nodes.

  OpenSees PlotModel.tcl
  python plot/PlotModelSketchPaper.py [--fill none|gray|color] [in.json] [out_stem]

Default output (--fill none, lines only):
  plot/out/profile{N}/elevation/{Shin|ASDEA}/elevation_paper.{png,pdf}
  (--fill gray / color -> elevation_paper_gray / elevation_paper_color)
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from matplotlib import font_manager
from matplotlib.colors import to_rgba
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.legend_handler import HandlerBase
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Rectangle

from paths import elevation_dir
from pt39_outline import pt39_outline
from PlotModelSketch import (
    DEFAULT_JSON,
    SSI_STYLE as SSI_STYLE_BASE,
    SSI_TYPE_ORDER,
    STYLE as STYLE_BASE,
    display_nodes,
    draw_rot_spiral,
    draw_trans_coil,
    group_of_node,
    layer_style,
    load,
    node_map,
    parse_ssi_row,
)

# ------------------------------------------------------------
# 1. PAGE, FONT AND PANEL WINDOWS
# ------------------------------------------------------------

FIG_W = 6.0            # in
FIG_H = 6.0            # in
FONT_SIZE = 9          # pt
FONT_NAME = "NewComputerModernSans10"
FONT_DIR = Path.home() / "AppData" / "Local" / "Microsoft" / "Windows" / "Fonts"
FONT_FILES = ("NewCMSans10-Regular.otf", "NewCMSans10-Oblique.otf")

# Water is physical (flume); OpenSees only sees rho_w g h_w on the y=0 soil
# nodes (analysis/WaterSurfaceLoad.tcl), drawn as arrows onto the mudline.
# Stiff elastic frames share one color (solid); pier keeps a distinct orange.
# Hex matches the static-equilibrium / deck navy and a warmer pier orange.
STIFF_COLOR = "#0b3c5d"
PIER_COLOR = "#cd6416"
STYLE = {
    **STYLE_BASE,
    "deck": {**STYLE_BASE["deck"], "line": STIFF_COLOR, "node": STIFF_COLOR,
             "fill": "#90a4ae"},
    "cap": {**STYLE_BASE["cap"], "line": STIFF_COLOR, "node": STIFF_COLOR,
            "fill": "#90a4ae"},
    # stiff intermediate beam between fiber hinges
    "pier": {**STYLE_BASE["pier"], "line": PIER_COLOR, "node": PIER_COLOR,
             "fill": "#e0a070"},
}
SSI_STYLE = {**SSI_STYLE_BASE, "qz": {**SSI_STYLE_BASE["qz"], "line": "#d32f2f"}}

WATER = {"fill": "#b3e0f2", "alpha": 0.30, "line": "#1f78b4", "lw": 0.8,
         "ls": (0, (4, 2))}
SURCHARGE = {"color": "#0b3c5d", "lw": 0.5, "head": 3.5}

# Member extent shading (pier, cap, piles, deck outline):
#   color -- one tint per member (PlotModelSketch STYLE)
#   gray  -- one gray for all members
#   none  -- no shading; frame lines drawn thicker
FILL_MODES = ("color", "gray", "none")
FILL_DEFAULT = "none"
GRAY = {"fill": "#9e9e9e", "alpha": 0.35}

# Soil: earth tones, darker with depth; no blue so L2 never reads as water.
# Sand layers (profiles 1-3) get yellow-ochre.
SOIL_FILL = {
    "L2": "#f0e6c4",
    "L3": "#c2b577",
    "L3a": "#e3c46a",
    "L3b": "#d6ae4a",
    "L3c": "#c49a3a",
    "L5": "#b3a591",
}
SOIL_ALPHA = 0.48  # light face so pier / springs / BCs stay primary
SOIL_EDGE = "#8a7a64"  # softer mesh edges
# Shin free-field columns: same materials, t_FF = 10000 t_soil (Profiles.md)
FF_HATCH = {"pattern": "////", "edge": "#5a4a36"}

# SSI springs are zero length; draw each one hanging off its structure node.
#   p-y: horizontal, from the node outward
#   t-z, q-z: vertical, straight down when no member continues below;
#             otherwise a lead inward (away from p-y), then down, so the
#             coil does not sit on the pile line (q-z one lead past t-z)
# Shin base (soil/BuildShinLysmer.tcl): NF base nodes equalDOF UX to one
# primary node (NF dashpot); each FF column's inner base node equalDOF UX to
# its outer node (FF dashpot). Dashpots end on a fully fixed ghost node.
# Drawn in panel (a) below the base; sizes in m.
DASHPOT_LEN = 7.0      # m, rod + cylinder + rod
DASHPOT_H = 1.6        # m, cylinder height
DASHPOT_DROP = 3.4     # m, dashpot row below the base
# Path load on each Lysmer soil node: F = 2 c v(t) in UX (BuildShinLysmer.tcl)
# Arrow sits to the right of the dashpot (same elevation as the cylinder).
EQ_FORCE = {"len": 4.5, "gap": 1.2, "lw": 0.9, "head": 7.0}
EQDOF = {"color": "#444444", "lw": 0.6, "ls": (0, (1.5, 1.2)), "drop": 1.7}  # dense short dash
# HoldPierBase.tcl: pier-base UX+UY held at gravity; RZ free (pin). Soil base: UY stays fixed in EQ.
PIN = {"color": "#222222", "lw": 0.9}
UY_FIX = {"color": "#222222", "lw": 0.5, "h": 0.55}  # m, soil-base UY roller height
# OpenFresco expElement (generic UX at pier-top inner node): force feedback
# from the physical subassembly. Glyph = fixed -- actuator -- eyelet at pier.
# Panel (d) annotates the prototype force λ_L^3 f_p applied at that node.
EXP = {"color": "#ec2a10", "len": 1.55, "h": 0.34, "lw": 1.2}
# Force callout (Froude): prototype force = λ_L^3 f_p, with f_p the resisting
# force assembled on the physical DOFs (inertia, damping, hydro).
# (a) one-line; (b) arrow only; (d) formula by arrow, name under it inside axes.
FP_ARROW = {"len": 0.95, "gap": 0.06, "lw": 1.35, "head": 10.0,
            "label_dx": 0.08, "label_dy": 0.22, "name_dx": 1.85}
FP_TXT = r"Physical subassembly force, $\lambda_L^{3}\mathbf{f}_{p}$"
FP_TXT_NAME = "Physical subassembly force"
FP_TXT_MATH = r"$\lambda_L^{3}\mathbf{f}_{p}$"
# force_fp_label modes: "full" | "two" | "none"


SSI_SCALE = 0.7        # (-) overall glyph size
SSI_LEN = 0.55 * SSI_SCALE     # m, drawn p-y coil length
SSI_LEN_V = 0.42 * SSI_SCALE   # m, drawn t-z / q-z length (cap face stations 0.5 m apart)
PAIR_RING = {"grow": 3.0, "mew": 0.9}   # pt, soil ring = structure ms + grow
SSI_LEAD = 0.22 * SSI_SCALE    # m, loop offset for t-z beside a member (q-z: LOOP["qz_off"])
# Spring drawn as a closed loop: node -> coil -> turn -> return -> same node
LOOP = {
    "g": 0.075,        # m, half gap between coil track and return track
    "amp": 0.05,       # m, zigzag amplitude on the coil track
    "lead": 0.09,      # m, node to start of the tracks (along the spring)
    "n_zig": 7,
    "lw": 0.9,         # pt
    "qz_off": 0.40,    # m, q-z loop offset beside a member (clears the t-z loop)
}
LW_SCALE_NONE = 2.9    # (-) frame line width factor when fill = none
BOX = {"color": "#222222", "lw": 0.7, "ls": (0, (5, 2, 1, 2))}  # dash-dot zoom window

# Page margins and gaps (in)
LEFT = 0.48
RIGHT = 0.08
TOP = 0.22
GAP_TITLE = 0.20       # room for panel letter titles
GAP_XLABEL = 0.38      # room for tick labels + x label
GAP_COL = 0.50         # between (b) and (c), room for (c) y label

# Panel (c) window around the cap and pier-base hinge (m)
DETAIL_X = (-3.3, 3.3)
DETAIL_Y = (-2.5, 1.4)
# Panel (d) crop around pier-top inner node + soffit (m)
# Node 4 ~ (0, 6.55); room left for name under the force arrow.
TOP_X = (-5.00, 0.9)
TOP_Y = (5.70, 7.55)


def use_paper_font() -> None:
    """
    Register the per-user New CM Sans files and set 9 pt text + math.

    Args:    none
    Returns: none (updates rcParams)
    """
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
        "legend.fontsize": FONT_SIZE,
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


def soil_fill(name: str, profile: int | None) -> str:
    """
    Paper fill color for one soil layer.

    Args:    name  layer name (L2, L3, L5, ...), profile  soil profile number
    Returns: hex color
    """
    if profile == 1 and name == "L3":
        return SOIL_FILL["L3a"]
    if name.startswith("L5"):
        return SOIL_FILL["L5"]
    return SOIL_FILL.get(name, "#bdb3a4")


def nodes_with_member_below(data: dict) -> set[tuple[float, float]]:
    """
    Structure points that have a vertical pile or cap member going down.

    Args:    data
    Returns: set of rounded (x, y), m
    """
    xy = node_map(data)
    below = set()
    for _e, ni, nj, grp in data["elements"]:
        if grp not in ("pile", "cap") or int(ni) not in xy or int(nj) not in xy:
            continue
        (xa, ya), (xb, yb) = xy[int(ni)], xy[int(nj)]
        if abs(xa - xb) < 1.0e-6 and abs(ya - yb) > 1.0e-6:
            x_up, y_up = (xa, ya) if ya > yb else (xb, yb)
            below.add((round(x_up, 3), round(y_up, 3)))
    return below


def spring_station_xy(data: dict) -> set[tuple[float, float]]:
    """
    Points where a soil node, a zero-length spring, and a structure node coincide.

    Args:    data
    Returns: set of rounded (x, y), m
    """
    pts = set()
    for row in data.get("ssi_springs", []):
        xp, yp, _xi, _yi, stype, _kind = parse_ssi_row(row)
        if stype != "none":
            pts.add((round(xp, 6), round(yp, 6)))
    return pts


def loop_path(L: float, off: float, g: float, amp: float, lead: float,
              n_zig: int) -> np.ndarray:
    """
    Closed spring loop in local (s, n) coordinates, starting and ending at (0, 0).

    The coil runs on the track n = off + g, a half-circle turns at the far end,
    and a plain return runs on n = off - g back to the node.

    Args:    L  coil length along s, off  loop centre offset along n,
             g  half gap, amp  zigzag amplitude, lead  node to track start (all
             in the same length unit), n_zig  zigzag count
    Returns: (k, 2) array of (s, n)
    """
    s0, s1 = lead, lead + L
    n_top, n_bot = off + g, off - g
    pts = [(0.0, 0.0), (s0, n_top)]
    for i in range(1, n_zig + 1):
        t = i / (n_zig + 1)
        pts.append((s0 + t * L, n_top + (amp if i % 2 else -amp)))
    pts.append((s1, n_top))
    for phi in np.linspace(0.0, np.pi, 13)[1:-1]:
        pts.append((s1 + g * np.sin(phi), off + g * np.cos(phi)))
    pts += [(s1, n_bot), (s0, n_bot), (0.0, 0.0)]
    return np.asarray(pts)


def draw_loop_spring(a: plt.Axes, x: float, y: float, u: tuple, n: tuple,
                     L: float, off: float, color: str,
                     *, lw: float | None = None) -> None:
    """
    Draw one zero-length spring as a closed loop on its node.

    Args:    a, x, y  node (m), u  unit vector along the spring (coil direction),
             n  unit normal (side of the coil track and offset), L  coil length (m),
             off  loop offset along n (m), color, lw  line width (pt; default LOOP)
    Returns: none (updates a)
    """
    sn = loop_path(L, off, LOOP["g"], LOOP["amp"], LOOP["lead"], LOOP["n_zig"])
    xs = x + sn[:, 0] * u[0] + sn[:, 1] * n[0]
    ys = y + sn[:, 0] * u[1] + sn[:, 1] * n[1]
    a.plot(xs, ys, color=color, lw=LOOP["lw"] if lw is None else lw,
           zorder=5.5, solid_joinstyle="round")


def draw_ssi_attached(a: plt.Axes, data: dict) -> None:
    """
    Draw p-y / t-z / q-z coils attached to their structure node.

    Coincident duplicates (pile and cap springs at a pile head) draw once.

    Args:    a, data
    Returns: none (updates a)
    """
    below = nodes_with_member_below(data)
    done = set()
    for row in data.get("ssi_springs", []):
        xp, yp, _xi, _yi, stype, _kind = parse_ssi_row(row)
        if stype == "none":
            continue
        key = (round(xp, 3), round(yp, 3), stype)
        if key in done:
            continue
        done.add(key)
        color = SSI_STYLE[stype]["line"]
        out = 1.0 if xp >= -1.0e-6 else -1.0
        if stype in ("py", "pyliq", "py_elastic"):
            # outward, coil on the upper track
            draw_loop_spring(a, xp, yp, (out, 0.0), (0.0, 1.0), SSI_LEN, 0.0, color)
            continue
        # vertical: hang down; step inward when a pile or cap member runs below
        if (round(xp, 3), round(yp, 3)) in below:
            off = SSI_LEAD if stype in ("tz", "tzliq", "tz_elastic") else LOOP["qz_off"]
        else:
            off = 0.0
        draw_loop_spring(a, xp, yp, (0.0, -1.0), (-out, 0.0), SSI_LEN_V, off, color)


def draw_exp_element(a: plt.Axes, nodes: dict[int, tuple[float, float]],
                     *, label: bool = False, scale: float = 1.0,
                     lw: float | None = None, force_fp: bool = False,
                     force_fp_label: str = "full") -> None:
    """
    OpenFresco experimental element at the pier-top UX node.

    Default: fixed -- actuator -- eyelet (unused on the paper figure).
    force_fp: skip the glyph; draw the prototype force into the shared UX.

    Matches Run.tcl generic attachment: lumpedPlasticity node 4 (inner top),
    else node 5 (deck BC).

    Args:    a, nodes  display node map (m)
             label  unused
             scale  length multiplier; lw  line width override (pt)
             force_fp  arrow only (no actuator)
             force_fp_label  "full" | "two" | "none"
    Returns: none (updates a)
    """
    tag = 4 if 4 in nodes else (5 if 5 in nodes else None)
    if tag is None:
        return
    x0, y = nodes[tag]
    L = EXP["len"] * scale
    h = EXP["h"] * scale
    c = EXP["color"]
    line_w = EXP["lw"] if lw is None else lw
    r_eye = 0.38 * h
    z = 5.5

    if force_fp:
        # prototype force into the shared UX
        L_a = FP_ARROW["len"] * scale
        xa1 = x0 - FP_ARROW["gap"] * scale
        xa0 = xa1 - L_a
        a.add_patch(FancyArrowPatch(
            (xa0, y), (xa1, y), arrowstyle="-|>",
            mutation_scale=FP_ARROW["head"], lw=FP_ARROW["lw"],
            color=c, shrinkA=0.0, shrinkB=0.0, zorder=z + 3,
        ))
        xl = xa0 - FP_ARROW["label_dx"] * scale
        if force_fp_label == "full":
            a.text(
                xl, y, FP_TXT,
                color=c, fontsize=FONT_SIZE - 1, ha="right", va="center",
                zorder=z + 4, clip_on=False,
            )
        elif force_fp_label == "two":
            # formula left of shaft; name under arrow, shifted left ~"Physical"
            a.text(
                xl, y, FP_TXT_MATH,
                color=c, fontsize=FONT_SIZE - 1, ha="right", va="center",
                zorder=z + 4, clip_on=False,
            )
            a.text(
                0.5 * (xa0 + xa1) - FP_ARROW["name_dx"] * scale,
                y - FP_ARROW["label_dy"] * scale, FP_TXT_NAME,
                color=c, fontsize=FONT_SIZE - 1, ha="center", va="top",
                zorder=z + 4, clip_on=False,
            )
        return

    # Layout (left -> right): fixed | cylinder+caps | rod | eyelet@pier
    x_fix = x0 - L
    x_rod1 = x0 - r_eye                 # rod meets the eyelet
    x_cyl1 = x0 - 0.32 * L              # cylinder end toward pier
    x_cyl0 = x0 - 0.82 * L              # cylinder end toward fixed
    cap = 0.06 * L                      # end-cap thickness along x
    kw = {"color": c, "lw": line_w, "zorder": z, "solid_capstyle": "butt"}
    # cylinder body + slightly taller end caps (outline only)
    yb, yt = y - 0.45 * h, y + 0.45 * h
    yb_c, yt_c = y - 0.58 * h, y + 0.58 * h
    a.plot([x_cyl0, x_cyl1, x_cyl1, x_cyl0, x_cyl0],
           [yb, yb, yt, yt, yb], **kw)
    for x_cap in (x_cyl0, x_cyl1 - cap):
        a.plot([x_cap, x_cap + cap, x_cap + cap, x_cap, x_cap],
               [yb_c, yb_c, yt_c, yt_c, yb_c], **kw)
    # piston rod
    a.plot([x_cyl1, x_rod1], [y, y], color=c, lw=1.4 * line_w, zorder=z,
           solid_capstyle="butt")
    # clevis / eyelet at the pier (free end)
    a.add_patch(Circle(
        (x0, y), r_eye, facecolor="none", edgecolor=c, lw=line_w, zorder=z + 1,
    ))
    a.add_patch(Circle(
        (x0, y), 0.35 * r_eye, facecolor="none", edgecolor=c,
        lw=0.7 * line_w, zorder=z + 1.1,
    ))
    # fixed support on the left (replaces the icon's other eyelet)
    a.plot([x_fix, x_cyl0], [y, y], **kw)
    a.plot([x_fix, x_fix], [y - 0.9 * h, y + 0.9 * h], lw=1.0,
           color="#222222", zorder=z + 2)
    for k in range(4):
        yk = y - 0.75 * h + k * 0.5 * h
        a.plot([x_fix, x_fix - 0.55 * h], [yk, yk + 0.32 * h], lw=0.5,
               color="#222222", zorder=z + 2)



def draw_pin(a: plt.Axes, x: float, y: float, s: float,
             *, angle_deg: float = 0.0) -> None:
    """
    2D pin support (UX+UY fixed, RZ free) with tip at the node.

    Local glyph opens downward (angle 0). angle_deg rotates CCW; -90 puts
    the support to the left of the node (keeps the base ZLS spiral clear).

    Args:    a, x, y  tip (m), s  half-base width (m)
             angle_deg  CCW rotation from the default (down) orientation
    Returns: none (updates a)
    """
    h = 1.15 * s
    rad = np.radians(angle_deg)
    ca, sa = np.cos(rad), np.sin(rad)

    def xy(px: float, py: float) -> tuple[float, float]:
        return (x + ca * px - sa * py, y + sa * px + ca * py)

    a.add_patch(Polygon(
        [xy(0.0, 0.0), xy(-s, -h), xy(s, -h)],
        closed=True, facecolor="white", edgecolor=PIN["color"],
        lw=PIN["lw"], zorder=8,
    ))
    g0, g1 = xy(-1.35 * s, -h), xy(1.35 * s, -h)
    a.plot([g0[0], g1[0]], [g0[1], g1[1]], color=PIN["color"],
           lw=PIN["lw"], zorder=8)
    for k in range(5):
        p0 = xy(-1.1 * s + k * 0.55 * s, -h)
        p1 = xy(-1.1 * s + k * 0.55 * s - 0.25 * s, -h - 0.35 * s)
        a.plot([p0[0], p1[0]], [p0[1], p1[1]], color=PIN["color"],
               lw=0.5, zorder=8)


def draw_uy_roller(a: plt.Axes, x: float, y: float, s: float) -> None:
    """
    Vertical roller (UY fixed, UX free) under a soil-base node.

    Args:    a, x, y  contact on the base (m), s  glyph size (m)
    Returns: none (updates a)
    """
    c = UY_FIX["color"]
    a.add_patch(Circle((x, y - s), 0.45 * s, fill=False, edgecolor=c,
                       lw=UY_FIX["lw"], zorder=5))
    y0 = y - 1.55 * s
    a.plot([x - 0.9 * s, x + 0.9 * s], [y0, y0], color=c, lw=UY_FIX["lw"],
           zorder=5)
    for k in range(4):
        xk = x - 0.7 * s + k * 0.45 * s
        a.plot([xk, xk - 0.2 * s], [y0, y0 - 0.3 * s], color=c, lw=0.45,
               zorder=5)


def draw_dashpot(a: plt.Axes, x0: float, y: float, color: str,
                 *, eq_force: bool = True) -> None:
    """
    Horizontal Lysmer glyph matching BuildShinLysmer.tcl:

        fixed ghost -- Viscous -- soil base node  (+ Path load 2 c v)

    Soil attach is at x0 (end of the vertical lead). The dashpot extends
    to the left to a fixed support; the cylinder opens toward the free
    (soil) side; the 2cv arrow sits to the right of x0.

    Args:    a, x0, y  soil-side attach point (m), color
             eq_force  draw the earthquake force arrow
    Returns: none (updates a)
    """
    L, h = DASHPOT_LEN, DASHPOT_H
    # soil at x0; fixed support at x0 - L
    x_fix = x0 - L
    # cylinder: closed end toward fixed, open gap toward soil
    bx_c, bx_o = x0 - 0.7 * L, x0 - 0.3 * L
    xp = 0.5 * (bx_c + bx_o)  # piston at mid-cylinder
    kw = {"color": color, "zorder": 7, "solid_capstyle": "butt"}
    a.plot([x_fix, bx_c], [y, y], lw=1.2, **kw)
    a.plot([xp, x0], [y, y], lw=1.2, **kw)
    a.plot([bx_c, bx_c], [y - 0.5 * h, y + 0.5 * h], lw=1.0, **kw)
    a.plot([bx_c, bx_o], [y + 0.5 * h, y + 0.5 * h], lw=1.0, **kw)
    a.plot([bx_c, bx_o], [y - 0.5 * h, y - 0.5 * h], lw=1.0, **kw)
    a.plot([xp, xp], [y - 0.38 * h, y + 0.38 * h], lw=1.2, **kw)
    # fixed support on the left (ghost)
    a.plot([x_fix, x_fix], [y - 0.7 * h, y + 0.7 * h], lw=1.0, color="#222222",
           zorder=7)
    for k in range(4):
        yk = y - 0.7 * h + k * 0.47 * h
        a.plot([x_fix, x_fix - 0.9], [yk, yk + 0.35 * h], lw=0.5,
               color="#222222", zorder=7)
    if eq_force:
        # Path load on the soil node: arrow just to the right of the attach
        xa0 = x0 + EQ_FORCE["gap"]
        xa1 = xa0 + EQ_FORCE["len"]
        a.add_patch(FancyArrowPatch(
            (xa0, y), (xa1, y), arrowstyle="-|>",
            mutation_scale=EQ_FORCE["head"], lw=EQ_FORCE["lw"],
            color=color, shrinkA=0.0, shrinkB=0.0, zorder=8,
        ))


def draw_base_boundary(
    a: plt.Axes, data: dict, *, force_label: bool = True, force_under: bool = False,
    force_at: str = "NF",
) -> None:
    """
    Shin base: dashed equalDOF bus under the tied base nodes plus the
    three Lysmer dashpots. Each glyph is fixed (left) -- dashpot -- soil
    (right), with Path load F = 2 c v(t) arrowed at the soil end.

    Args:    a, data
             force_label  write the earthquake Path-load callout
             force_under  place that text under the glyph (else to the right)
             force_at  which dashpot gets the callout: "NF" | "leftmost"
    Returns: none (updates a)
    """
    dash = data.get("lysmer_dashpots", [])
    if not dash:
        return
    sz = data["sizes"]
    L_half = float(sz.get("L_half", 0.0))
    yb = min(float(L["y0"]) for L in data["soil_layers"])
    xs = sorted({
        round(q["xy"][i], 6)
        for q in data.get("soil_quads", [])
        for i in range(0, 8, 2)
        if abs(q["xy"][i + 1] - yb) < 1.0e-6
    })
    y_bus = yb - EQDOF["drop"]
    eq = {"color": EQDOF["color"], "lw": EQDOF["lw"], "ls": EQDOF["ls"], "zorder": 6}
    groups = [
        [x for x in xs if abs(x) < L_half - 1.0e-6],   # near field -> primary
        [x for x in xs if x <= -L_half + 1.0e-6],       # left FF column
        [x for x in xs if x >= L_half - 1.0e-6],        # right FF column
    ]
    for g in groups:
        if len(g) < 2:
            continue
        a.plot([g[0], g[-1]], [y_bus, y_bus], **eq)
        for x in g:
            a.plot([x, x], [yb, y_bus], **eq)

    # Soil base UY remains fixed in EQ (only UX is released for Lysmer)
    s_uy = UY_FIX["h"]
    step = max(1, len(xs) // 12)
    for x in xs[::step]:
        draw_uy_roller(a, x, yb, s_uy)

    color = STYLE["lysmer"]["line"]
    y_d = yb - DASHPOT_DROP
    eq_txt = r"Earthquake loading, $2c\,\dot{u}_g(t)$"
    if force_at == "leftmost":
        x_lab = min(float(d["x"]) for d in dash)
    else:
        x_lab = next(
            (float(d["x"]) for d in dash if d.get("role") == "NF"),
            float(dash[0]["x"]),
        )
    labeled = False
    for d in dash:
        x = float(d["x"])
        a.plot([x, x], [yb, y_d], color=color, lw=1.0, zorder=7)
        draw_dashpot(a, x, y_d, color)
        if force_label and (not labeled) and abs(x - x_lab) < 1.0e-9:
            if force_under:
                # center under dashpot + Path-load arrow, then nudge ~1 capital E
                # right so the leading "E" sits inside the axes
                xm = (x + 0.5 * (EQ_FORCE["gap"] + EQ_FORCE["len"] - DASHPOT_LEN)
                      + 0.65)
                a.text(xm, y_d - 0.55 * DASHPOT_H - 0.2, eq_txt, color=color,
                       fontsize=FONT_SIZE - 1, ha="center", va="top", zorder=8)
            else:
                xa = x + EQ_FORCE["gap"] + EQ_FORCE["len"]
                a.text(xa + 0.6, y_d, eq_txt, color=color,
                       fontsize=FONT_SIZE - 1, ha="left", va="center", zorder=8)
            labeled = True


# ------------------------------------------------------------
# 2. DRAW ONE PANEL
# ------------------------------------------------------------


def draw_model(
    a: plt.Axes,
    data: dict,
    nodes: dict[int, tuple[float, float]],
    *,
    xlim: tuple[float, float],
    lw: float,
    ms: float,
    quad_lw: float,
    coils: bool,
    dashpots: bool,
    fill: str = FILL_DEFAULT,
    soil_ms: float = 0.0,
    surcharge: bool = True,
    pin_s: float = 0.0,
    force_label: bool = True,
    force_under: bool = False,
    force_at: str = "NF",
    exp_element: bool = False,
    exp_label: bool = False,
    exp_scale: float = 1.0,
    exp_lw: float | None = None,
    exp_force_fp: bool = False,
    exp_force_fp_label: str = "full",
    equal_aspect: bool = True,
) -> None:
    """
    Draw soil, water, structure, hinges, and optional SSI coils / dashpots.

    Args:    a, data, nodes  display node map (m)
             xlim  panel x range (m); water band spans it, capped at the mesh
             lw, ms  frame line width (pt), node marker size (pt)
             quad_lw  soil mesh line width (pt)
             coils, dashpots  SSI springs / Lysmer glyphs
             fill  member shading
             soil_ms  soil-node marker size (pt; 0 = off)
             surcharge  water-pressure arrows on y=0
             pin_s  pier-base pin half-width (m; 0 = auto)
             force_label, force_under, force_at  earthquake Path-load callout
             exp_element  OpenFresco UX spring at the pier top
             exp_label, exp_scale, exp_lw, exp_force_fp, exp_force_fp_label
                 OpenFresco force-arrow options ("full"|"two"|"none")
             equal_aspect  True -> equal x/y scale (default); False for inset (d)
    Returns: none (updates a)
    """
    if fill == "none":
        lw *= LW_SCALE_NONE
    sz = data["sizes"]
    profile = data.get("soilProfile")
    H_pier = float(sz["H_pier"])
    springs = data.get("springs", [])
    tag_info = data.get("tags") or {}
    soil_base = int(tag_info.get("soil", 10000))
    spr_base = int(tag_info.get("spr", 20000))

    by_layer: dict[str, list] = {}
    for q in data.get("soil_quads", []):
        by_layer.setdefault(q["layer"], []).append(q)
    for nm, qs in by_layer.items():
        polys = [
            Polygon([(q["xy"][i], q["xy"][i + 1]) for i in range(0, 8, 2)])
            for q in qs
        ]
        # alpha on the face only so mesh edges stay opaque
        a.add_collection(PatchCollection(
            polys, facecolor=to_rgba(soil_fill(nm, profile), SOIL_ALPHA),
            edgecolor=SOIL_EDGE, linewidths=quad_lw, zorder=0,
        ))

    # Shin FF columns (|x| > L_half): hatch over the same fill (thick OOP)
    L_half = float(sz.get("L_half", 0.0) or 0.0)
    x_mesh = float(sz.get("xMeshHalf", 0.0) or 0.0)
    if L_half > 0.0 and x_mesh > L_half + 1.0e-6:
        ff_polys = []
        for q in data.get("soil_quads", []):
            xs = [q["xy"][i] for i in range(0, 8, 2)]
            if abs(0.5 * (min(xs) + max(xs))) > L_half - 1.0e-6:
                ff_polys.append(Polygon(
                    [(q["xy"][i], q["xy"][i + 1]) for i in range(0, 8, 2)]
                ))
        if ff_polys:
            a.add_collection(PatchCollection(
                ff_polys, facecolor="none", edgecolor=FF_HATCH["edge"],
                linewidths=0.0, hatch=FF_HATCH["pattern"], zorder=0.3,
            ))
            # NF | FF interface (helps when hatch is fine)
            yb = min(float(L["y0"]) for L in data["soil_layers"])
            yt = max(float(L["y1"]) for L in data["soil_layers"])
            for x_ff in (-L_half, L_half):
                a.plot([x_ff, x_ff], [yb, yt], color=FF_HATCH["edge"],
                       lw=0.5, ls=(0, (2, 1.5)), zorder=0.35)

    if soil_ms > 0.0:
        soil_xy = {
            (round(q["xy"][i], 6), round(q["xy"][i + 1], 6))
            for q in data.get("soil_quads", [])
            for i in range(0, 8, 2)
        }
        # soil nodes at spring stations sit under a structure node: ring them
        paired = spring_station_xy(data) & soil_xy if coils else set()
        xs, ys = zip(*(soil_xy - paired))
        a.plot(xs, ys, "o", ls="None", color=SOIL_EDGE, ms=soil_ms, zorder=0.5)
        if paired:
            xs, ys = zip(*paired)
            a.plot(xs, ys, "o", ls="None", ms=ms + PAIR_RING["grow"],
                   markerfacecolor="white", markeredgecolor=SOIL_EDGE,
                   markeredgewidth=PAIR_RING["mew"], zorder=4.9)

    h_water = float(sz.get("h_water", 0.0) or 0.0)
    if h_water > 0.0:
        x_mesh = float(sz.get("xMeshHalf", max(abs(xlim[0]), abs(xlim[1]))))
        x0_w, x1_w = max(xlim[0], -x_mesh), min(xlim[1], x_mesh)
        a.add_patch(Rectangle(
            (x0_w, 0.0), x1_w - x0_w, h_water,
            facecolor=WATER["fill"], edgecolor="none",
            alpha=WATER["alpha"], zorder=1,
        ))
        a.plot([x0_w, x1_w], [h_water, h_water], color=WATER["line"],
               lw=WATER["lw"], ls=WATER["ls"], zorder=1.5)
        # surcharge arrows: one per y=0 soil node, from the free surface down
        x_surf = sorted({
            round(q["xy"][i], 6)
            for q in data.get("soil_quads", [])
            for i in range(0, 8, 2)
            if abs(q["xy"][i + 1]) < 1.0e-6
        })
        for xs in x_surf if surcharge else []:
            if not (xlim[0] <= xs <= xlim[1]):
                continue
            a.add_patch(FancyArrowPatch(
                (xs, h_water), (xs, 0.0), arrowstyle="-|>",
                mutation_scale=SURCHARGE["head"], lw=SURCHARGE["lw"],
                color=SURCHARGE["color"], shrinkA=0.0, shrinkB=0.0,
                zorder=4,
            ))

    for fill_row in data["fills"] if fill != "none" else []:
        g = fill_row["group"]
        if g == "deck":
            continue
        st = STYLE[g] if fill == "color" else GRAY
        xc, w = float(fill_row["xc"]), float(fill_row["width"])
        y0, y1 = float(fill_row["y0"]), float(fill_row["y1"])
        if g == "pier" and springs and 2 in nodes and 4 in nodes:
            y0, y1 = nodes[2][1], nodes[4][1]
        a.add_patch(Rectangle(
            (xc - 0.5 * w, min(y0, y1)), w, abs(y1 - y0),
            facecolor=st["fill"], edgecolor="none",
            alpha=st["alpha"], zorder=2,
        ))

    if "dw_deck" in sz and fill != "none":
        st = STYLE["deck"] if fill == "color" else GRAY
        for poly in pt39_outline(
            float(sz["dw_deck"]), float(sz["dd_deck"]),
            float(sz["sw_deck"]), float(sz["cw_deck"]),
            float(sz["td_deck"]), float(sz["ts_deck"]), float(sz["tw_deck"]),
            y0=H_pier,
        ):
            poly.set_facecolor(st["fill"])
            poly.set_edgecolor(st.get("line", "none"))
            poly.set_alpha(st["alpha"])
            poly.set_linewidth(0.4)
            poly.set_zorder(2)
            a.add_patch(poly)

    segs, colors = [], []
    for _e, ni, nj, grp in data["elements"]:
        if grp in ("spring", "ssi_spring", "soil", "soil_bnd", "other"):
            continue
        if int(ni) not in nodes or int(nj) not in nodes:
            continue
        segs.append([nodes[int(ni)], nodes[int(nj)]])
        colors.append(STYLE.get(grp, STYLE["pile"])["line"])
    if segs:
        a.add_collection(LineCollection(segs, colors=colors, linewidths=lw,
                                        zorder=3))

    for _e, ni, nj, _sx, _sy in springs:
        if int(ni) in nodes and int(nj) in nodes:
            x0, y0 = nodes[int(ni)]
            x1, y1 = nodes[int(nj)]
            draw_rot_spiral(a, x0, y0, x1, y1, STYLE["spring"]["line"])

    if exp_element:
        draw_exp_element(a, nodes, label=exp_label, scale=exp_scale, lw=exp_lw,
                         force_fp=exp_force_fp, force_fp_label=exp_force_fp_label)

    if coils:
        draw_ssi_attached(a, data)

    if dashpots:
        draw_base_boundary(a, data, force_label=force_label,
                           force_under=force_under, force_at=force_at)

    # HoldPierBase: both ends of the base ZLS held in UX+UY (RZ free -> pin).
    # Pier-top hinge (nodes 4--5) is not held. Upper pin is rotated -90 deg
    # (support to the left) so the glyph does not cover the spiral.
    if pin_s > 0.0:
        for _e, ni, nj, _sx, _sy in springs:
            yi = nodes.get(int(ni), (0.0, 1.0e9))[1]
            yj = nodes.get(int(nj), (0.0, 1.0e9))[1]
            if abs(yi) > 1.0 and abs(yj) > 1.0:
                continue
            tags = [t for t in (int(ni), int(nj)) if t in nodes]
            if not tags:
                continue
            y_top = max(nodes[t][1] for t in tags)
            for tag in tags:
                # upper pin: support to the left (clear of the spiral)
                ang = -90.0 if nodes[tag][1] >= y_top - 1.0e-9 else 0.0
                draw_pin(a, nodes[tag][0], nodes[tag][1], pin_s,
                         angle_deg=ang)

    if ms > 0.0:
        for tag, (x, y) in nodes.items():
            if tag >= spr_base:
                continue
            g = group_of_node(tag, soil_base, spr_base)
            if g in STYLE and "node" in STYLE[g]:
                a.plot(
                    x, y, "o", color=STYLE[g]["node"], ms=ms, zorder=5,
                    markeredgecolor="white", markeredgewidth=0.25,
                )

    a.set_xlim(*xlim)
    if equal_aspect:
        a.set_aspect("equal", adjustable="box")
    else:
        a.set_aspect("auto")


def mark_window(a: plt.Axes, x: tuple, y: tuple, label: str,
                *, corner: str = "tr") -> None:
    """
    Dashed box showing where a zoom panel sits, with its panel letter.

    Args:    a, x (x0, x1), y (y0, y1) in m, label  e.g. "(b)"
             corner  text anchor on the box: "tr" | "br" | "tl" | "bl"
    Returns: none (updates a)
    """
    a.add_patch(Rectangle(
        (x[0], y[0]), x[1] - x[0], y[1] - y[0], fill=False,
        edgecolor=BOX["color"], lw=BOX["lw"], ls=BOX["ls"], zorder=9,
    ))
    # place the letter just outside the chosen corner
    loc = {
        "tr": (x[1], y[1], "left", "top"),
        # outside the box, just past the bottom-right corner
        "br": (x[1], y[0], "left", "top"),
        "tl": (x[0], y[1], "right", "top"),
        "bl": (x[0], y[0], "right", "bottom"),
    }
    tx, ty, ha, va = loc[corner]
    if corner == "br":
        a.text(tx + 0.04, ty - 0.02, f" {label}", ha=ha, va=va,
               fontsize=FONT_SIZE, zorder=9)
        return
    dx = -0.08 if ha == "right" else 0.08
    dy = 0.08 if va == "bottom" else -0.08
    pad = f" {label}" if ha == "left" else f"{label} "
    a.text(tx + dx, ty + dy, pad, ha=ha, va=va, fontsize=FONT_SIZE, zorder=9)


# ------------------------------------------------------------
# 3. LEGEND ENTRIES
# ------------------------------------------------------------


def legend_handles(data: dict, fill: str = FILL_DEFAULT) -> tuple[list, list]:
    """
    Build the legend list in reading order (structure, water, soil, springs).

    Short labels in many rows; the page layout places them full-width under
    panels (b)+(c) so long text does not run past the right edge.

    Args:    data, fill  member shading, one of FILL_MODES
    Returns: (handles, labels); gray mode pairs a gray patch with each line
    """
    sz = data["sizes"]
    profile = data.get("soilProfile")
    hs, labels = [], []
    lw_leg = 2.2 if fill == "none" else 1.1

    # Deck + pile cap share STIFF_COLOR → one legend row; pier keeps its color.
    members = [
        ("deck", "Deck and pile cap (stiff elastic frame)"),
        ("pier", "Pier (stiff elastic frame)"),
        ("pile", "Piles (disp.-based frame)"),
    ]
    for g, lab in members:
        if fill == "color":
            h = mpatches.Patch(facecolor=STYLE[g]["fill"], alpha=0.6,
                               edgecolor=STYLE[g]["line"], lw=0.5)
        else:
            # line + node dots at both ends (matches panel markers)
            h = MemberProxy(STYLE[g]["line"], STYLE[g]["node"], lw_leg,
                            ls=STYLE[g].get("ls", "-"))
            if fill == "gray":
                h = (mpatches.Patch(facecolor=GRAY["fill"], alpha=GRAY["alpha"],
                                    edgecolor="none"), h)
        hs.append(h)
        labels.append(lab)

    h_water = float(sz.get("h_water", 0.0) or 0.0)
    if h_water > 0.0:
        hs.append(WaterProxy())
        labels.append("Water surcharge")

    # Paper L1..Ln from the top (model tags L2 / L3 / L5 for profile 4)
    tops: dict[str, float] = {}
    for L in data.get("soil_layers", []):
        tops[L["name"]] = max(tops.get(L["name"], -1.0e9), float(L["y1"]))
    layer_names = sorted(tops, key=lambda n: -tops[n])
    if layer_names:
        fills = [soil_fill(nm, profile) for nm in layer_names]
        hs.append(SoilLayersProxy(fills))
        labels.append("Soil (SSPquad)")
    L_half = float(sz.get("L_half", 0.0) or 0.0)
    x_mesh = float(sz.get("xMeshHalf", 0.0) or 0.0)
    if L_half > 0.0 and x_mesh > L_half + 1.0e-6:
        hs.append(FFHatchProxy())
        labels.append("Free-field (stiffer SSPquad)")

    if data.get("springs"):
        hs.append(SpiralProxy(STYLE["spring"]["line"]))
        labels.append("Fiber section (zero-length)")

    names = {"py": "p-y", "tz": "t-z", "qz": "q-z"}
    present = {parse_ssi_row(r)[4] for r in data.get("ssi_springs", [])}
    for key in ("py", "tz", "qz"):
        if key in present or any(
            p.startswith(key) for p in present if p != "none"
        ):
            hs.append(SpringProxy(SSI_STYLE[key]["line"]))
            labels.append(f"{names[key]} springs")

    if data.get("lysmer_dashpots"):
        hs.append(DashpotProxy(STYLE["lysmer"]["line"]))
        labels.append("Lysmer dashpot")
        hs.append(Line2D([0], [0], color=EQDOF["color"], lw=EQDOF["lw"] + 0.2,
                         ls=EQDOF["ls"]))
        labels.append("equalDOF")
        hs.append(PinProxy())
        labels.append("Pier-base pin")
        hs.append(RollerProxy())
        labels.append("Soil-base roller")
    return hs, labels


def soil_desc(name: str, profile: int | None) -> str:
    """
    Short soil description for a legend entry.

    Args:    name  model layer tag, profile  soil profile number
    Returns: e.g. "soft clay"
    """
    if name == "L2" and profile == 4:
        return "soft clay"
    if name == "L3":
        return "sand" if profile == 1 else "medium clay"
    if name.startswith("L5"):
        return "stiff clay" if profile in (3, 4) else "dense sand"
    return layer_style(name, profile=profile)["label"].split(" ", 1)[-1]


class MemberProxy:
    """Legend stand-in for a frame member: short line with nodes at both ends."""

    def __init__(self, line: str, node: str, lw: float, ls: str | tuple = "-"):
        self.line = line
        self.node = node
        self.lw = lw
        self.ls = ls


class SpringProxy:
    """Legend stand-in for a translational spring (zigzag)."""

    def __init__(self, color: str):
        self.color = color
    """Legend stand-in for a translational spring (zigzag)."""

    def __init__(self, color: str):
        self.color = color


class SpiralProxy(SpringProxy):
    """Legend stand-in for a rotational fiber section (spiral)."""


class ExpProxy(SpringProxy):
    """Legend stand-in for OpenFresco (fixed -- actuator -- free pier)."""


class DashpotProxy(SpringProxy):
    """Legend stand-in for a Lysmer dashpot."""


class WaterProxy:
    """Legend stand-in for the water surcharge (blue box, arrow down)."""


class FrameTrioProxy:
    """Legend stand-in for deck / pier / cap as three parallel colored lines."""

    def __init__(self, colors: list[str], lw: float):
        self.colors = list(colors)
        self.lw = lw


class SpringsTrioProxy:
    """Legend stand-in for p-y / t-z / q-z as three stacked spring loops."""

    def __init__(self, colors: list[str]):
        self.colors = list(colors)


class SoilLayersProxy:
    """Legend stand-in for stacked soil fills (top to bottom, left to right)."""

    def __init__(self, fills: list[str]):
        self.fills = list(fills)


class FFHatchProxy:
    """Legend stand-in for Shin free-field column hatch."""


class PinProxy:
    """Legend stand-in for the pier-base pin (UX+UY held)."""


class RollerProxy:
    """Legend stand-in for soil-base UY rollers."""


class HandlerFrameTrio(HandlerBase):
    """Draw three short horizontal lines (deck, pier, cap colors)."""

    def create_artists(self, legend, orig, xd, yd, w, h, fontsize, trans):
        n = len(orig.colors)
        arts = []
        for i, c in enumerate(orig.colors):
            y = -yd + h * (n - i) / (n + 1)
            arts.append(Line2D([-xd, -xd + w], [y, y], color=c, lw=orig.lw,
                               transform=trans))
        return arts


class HandlerSpringsTrio(HandlerBase):
    """Draw up to three spring loops stacked in the handle box."""

    def create_artists(self, legend, orig, xd, yd, w, h, fontsize, trans):
        n = len(orig.colors)
        arts = []
        for i, c in enumerate(orig.colors):
            yc = -yd + h * (n - i) / (n + 1)
            g = 0.18 * h
            sn = loop_path(0.7 * w, 0.0, g, 0.12 * h, 0.06 * w, 5)
            arts.append(Line2D(-xd + 0.05 * w + sn[:, 0], yc + sn[:, 1],
                               color=c, lw=LOOP["lw"], solid_joinstyle="round",
                               transform=trans))
        return arts


class HandlerSoilLayers(HandlerBase):
    """Draw three (or n) soil color squares in one handle."""

    def create_artists(self, legend, orig, xd, yd, w, h, fontsize, trans):
        n = max(len(orig.fills), 1)
        gap = 0.08 * w
        side = min(h, (w - (n - 1) * gap) / n)
        y0 = -yd + 0.5 * (h - side)
        arts = []
        for i, fill in enumerate(orig.fills):
            x0 = -xd + i * (side + gap)
            arts.append(Rectangle(
                (x0, y0), side, side,
                facecolor=to_rgba(fill, SOIL_ALPHA), edgecolor=SOIL_EDGE,
                linewidth=0.4, transform=trans,
            ))
        return arts


class HandlerFFHatch(HandlerBase):
    """Draw a hatched square for the free-field columns."""

    def create_artists(self, legend, orig, xd, yd, w, h, fontsize, trans):
        return [Rectangle(
            (-xd, -yd), w, h,
            facecolor=to_rgba(SOIL_FILL["L3"], SOIL_ALPHA),
            edgecolor=FF_HATCH["edge"], linewidth=0.4,
            hatch=FF_HATCH["pattern"], transform=trans,
        )]


class HandlerPin(HandlerBase):
    """Draw a small pin triangle with ground hatch."""

    def create_artists(self, legend, orig, xd, yd, w, h, fontsize, trans):
        xc, tip = -xd + 0.5 * w, -yd + 0.95 * h
        s, hh = 0.35 * w, 0.7 * h
        y0 = tip - hh
        c = PIN["color"]
        return [
            Polygon([(xc, tip), (xc - s, y0), (xc + s, y0)], closed=True,
                    facecolor="white", edgecolor=c, lw=0.9, transform=trans),
            Line2D([xc - 1.2 * s, xc + 1.2 * s], [y0, y0], color=c, lw=0.9,
                   transform=trans),
        ]


class HandlerRoller(HandlerBase):
    """Draw a small UY roller (circle on a hatched ground line)."""

    def create_artists(self, legend, orig, xd, yd, w, h, fontsize, trans):
        xc, yc = -xd + 0.5 * w, -yd + 0.55 * h
        r = 0.22 * h
        y0 = -yd + 0.12 * h
        c = UY_FIX["color"]
        return [
            Circle((xc, yc), r, fill=False, edgecolor=c, lw=0.8, transform=trans),
            Line2D([xc - 1.6 * r, xc + 1.6 * r], [y0, y0], color=c, lw=0.8,
                   transform=trans),
        ]


class HandlerMember(HandlerBase):
    """Draw a short member line with filled node circles at both ends."""

    def create_artists(self, legend, orig, xd, yd, w, h, fontsize, trans):
        yc = -yd + 0.5 * h
        x0, x1 = -xd + 0.10 * w, -xd + 0.90 * w
        r = 0.30 * h
        return [
            Line2D([x0, x1], [yc, yc], color=orig.line, lw=orig.lw,
                   ls=getattr(orig, "ls", "-"),
                   transform=trans, solid_capstyle="butt"),
            Circle((x0, yc), r, facecolor=orig.node, edgecolor="white",
                   linewidth=0.45, transform=trans),
            Circle((x1, yc), r, facecolor=orig.node, edgecolor="white",
                   linewidth=0.45, transform=trans),
        ]


class HandlerSpring(HandlerBase):
    """Draw the closed spring loop (same shape as the panels) in the handle box."""

    def create_artists(self, legend, orig, xd, yd, w, h, fontsize, trans):
        x0, yc = -xd + 0.06 * w, -yd + 0.5 * h
        g = 0.32 * h
        sn = loop_path(0.62 * w, 0.0, g, 0.2 * h, 0.08 * w, LOOP["n_zig"])
        return [Line2D(x0 + sn[:, 0], yc + sn[:, 1], color=orig.color,
                       lw=LOOP["lw"], solid_joinstyle="round", transform=trans)]


class HandlerExp(HandlerBase):
    """Draw fixed -- actuator outline -- eyelet (OpenFresco) in the handle box."""

    def create_artists(self, legend, orig, xd, yd, w, h, fontsize, trans):
        xf, x_eye = -xd + 0.06 * w, -xd + 0.92 * w
        yc = -yd + 0.5 * h
        r = 0.22 * h
        x_cyl0, x_cyl1 = xf + 0.12 * w, x_eye - 0.28 * w
        cap = 0.06 * w
        c = orig.color
        yb, yt = yc - 0.28 * h, yc + 0.28 * h
        yb_c, yt_c = yc - 0.38 * h, yc + 0.38 * h
        return [
            Line2D([xf, x_cyl0], [yc, yc], color=c, lw=1.0, transform=trans),
            Line2D([x_cyl0, x_cyl1, x_cyl1, x_cyl0, x_cyl0],
                   [yb, yb, yt, yt, yb], color=c, lw=1.0, transform=trans),
            Line2D([x_cyl0, x_cyl0 + cap, x_cyl0 + cap, x_cyl0, x_cyl0],
                   [yb_c, yb_c, yt_c, yt_c, yb_c], color=c, lw=1.0,
                   transform=trans),
            Line2D([x_cyl1 - cap, x_cyl1, x_cyl1, x_cyl1 - cap, x_cyl1 - cap],
                   [yb_c, yb_c, yt_c, yt_c, yb_c], color=c, lw=1.0,
                   transform=trans),
            Line2D([x_cyl1, x_eye - r], [yc, yc], color=c, lw=1.4,
                   transform=trans),
            Circle((x_eye, yc), r, facecolor="none", edgecolor=c, lw=1.0,
                   transform=trans),
            Circle((x_eye, yc), 0.35 * r, facecolor="none", edgecolor=c,
                   lw=0.7, transform=trans),
            Line2D([xf, xf], [yc - 0.4 * h, yc + 0.4 * h], color="#222222",
                   lw=0.9, transform=trans),
        ]


class HandlerSpiral(HandlerBase):
    """Draw a two-turn spiral with short leads."""

    def create_artists(self, legend, orig, xd, yd, w, h, fontsize, trans):
        xc, yc = -xd + 0.5 * w, -yd + 0.5 * h
        th = np.linspace(0.0, 4.0 * np.pi, 80)
        r = 0.08 * h + 0.5 * h * th / th[-1]
        xs = xc + r * np.cos(th)
        ys = yc + r * np.sin(th)
        return [
            Line2D(xs, ys, color=orig.color, lw=1.1, transform=trans),
            Line2D([-xd, xc - 0.55 * h], [yc, yc], color=orig.color, lw=1.0,
                   transform=trans),
            Line2D([xs[-1], -xd + w], [ys[-1], yc], color=orig.color, lw=1.0,
                   transform=trans),
        ]


class HandlerDashpot(HandlerBase):
    """Draw a horizontal dashpot (gap toward free/right; no force arrow)."""

    def create_artists(self, legend, orig, xd, yd, w, h, fontsize, trans):
        x0, yc = -xd, -yd + 0.5 * h
        bx_c, bx_o = x0 + 0.25 * w, x0 + 0.65 * w
        xp = 0.5 * (bx_c + bx_o)
        c = orig.color
        return [
            Line2D([x0, bx_c], [yc, yc], color=c, lw=1.2, transform=trans),
            Line2D([xp, x0 + w], [yc, yc], color=c, lw=1.2, transform=trans),
            Line2D([bx_c, bx_c], [yc - 0.4 * h, yc + 0.4 * h],
                   color=c, lw=1.0, transform=trans),
            Line2D([bx_c, bx_o], [yc + 0.4 * h, yc + 0.4 * h],
                   color=c, lw=1.0, transform=trans),
            Line2D([bx_c, bx_o], [yc - 0.4 * h, yc - 0.4 * h],
                   color=c, lw=1.0, transform=trans),
            Line2D([xp, xp], [yc - 0.3 * h, yc + 0.3 * h],
                   color=c, lw=1.2, transform=trans),
        ]


class HandlerWater(HandlerBase):
    """Draw a light blue box, dashed free surface, and a down arrow."""

    def create_artists(self, legend, orig, xd, yd, w, h, fontsize, trans):
        x0, y0 = -xd, -yd
        box = Rectangle((x0, y0), w, h, facecolor=WATER["fill"], alpha=0.6,
                        edgecolor="none", transform=trans)
        top = Line2D([x0, x0 + w], [y0 + h, y0 + h], color=WATER["line"],
                     lw=WATER["lw"], ls=WATER["ls"], transform=trans)
        arrow = FancyArrowPatch(
            (x0 + 0.5 * w, y0 + h), (x0 + 0.5 * w, y0), arrowstyle="-|>",
            mutation_scale=5.0, lw=0.7, color=SURCHARGE["color"],
            shrinkA=0.0, shrinkB=0.0, transform=trans,
        )
        return [box, top, arrow]


LEGEND_HANDLERS = {
    MemberProxy: HandlerMember(),
    SpringProxy: HandlerSpring(),
    ExpProxy: HandlerExp(),
    SpiralProxy: HandlerSpiral(),
    DashpotProxy: HandlerDashpot(),
    WaterProxy: HandlerWater(),
    FrameTrioProxy: HandlerFrameTrio(),
    SpringsTrioProxy: HandlerSpringsTrio(),
    SoilLayersProxy: HandlerSoilLayers(),
    FFHatchProxy: HandlerFFHatch(),
    PinProxy: HandlerPin(),
    RollerProxy: HandlerRoller(),
}


# ------------------------------------------------------------
# 4. PAGE LAYOUT
# ------------------------------------------------------------


def plot(data: dict, out_stem: Path, fill: str = FILL_DEFAULT) -> None:
    """
    Lay out panels (a)-(d) and the legend on a 6 x 6 in page.

    Args:    data  decoded sketch dictionary; out_stem  path without suffix
             fill  member shading, one of FILL_MODES
    Returns: none (writes .png, .pdf, .svg)
    """
    use_paper_font()
    sz = data["sizes"]
    nodes = display_nodes(node_map(data), data.get("springs", []))

    y_bot = min(float(L["y0"]) for L in data["soil_layers"])
    y_top = max(y for _, y in nodes.values())

    # (a) full domain, true scale, full width
    # room on the left for fixed supports, on the right for 2cv arrows
    x_half_a = (float(sz.get("xMeshHalf", sz.get("L_half", 60.0)))
                + max(DASHPOT_LEN, EQ_FORCE["gap"] + EQ_FORCE["len"]) + 2.0)
    ylim_a = (y_bot - DASHPOT_DROP - DASHPOT_H - 1.0, y_top + 1.5)
    w_a = FIG_W - LEFT - RIGHT
    h_a = w_a * (ylim_a[1] - ylim_a[0]) / (2.0 * x_half_a)

    # (b) near field, full remaining height (same as pre-(d) layout)
    # pad past the Lysmer glyph so the EQ callout stays inside the axes
    x_half_b = max(
        0.5 * float(sz["dw_deck"]) + 0.8,
        DASHPOT_LEN + 1.4,
        EQ_FORCE["gap"] + EQ_FORCE["len"] + 0.8,
    )
    ylim_b = (y_bot - DASHPOT_DROP - DASHPOT_H - 1.6, y_top + 0.8)
    top_b = TOP + h_a + GAP_XLABEL + GAP_TITLE
    h_b = FIG_H - top_b - GAP_XLABEL
    w_b = h_b * (2.0 * x_half_b) / (ylim_b[1] - ylim_b[0])

    # (c) full column minus a right strip so (d) can overhang into white
    left_c = LEFT + w_b + GAP_COL
    right_room = 0.28  # in, white margin for (d) to sit partly outside (c)
    w_c = FIG_W - left_c - right_room
    h_c = w_c * (DETAIL_Y[1] - DETAIL_Y[0]) / (DETAIL_X[1] - DETAIL_X[0])

    def rect(left: float, top: float, w: float, h: float) -> list[float]:
        """Inches from the top-left corner -> figure fraction [l, b, w, h]."""
        return [left / FIG_W, 1.0 - (top + h) / FIG_H, w / FIG_W, h / FIG_H]

    fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=300)
    ax_a = fig.add_axes(rect(LEFT, TOP, w_a, h_a))
    ax_b = fig.add_axes(rect(LEFT, top_b, w_b, h_b))
    ax_c = fig.add_axes(rect(left_c, top_b, w_c, h_c))

    # (d) inset straddling (c)'s top-right; wide enough for name under the arrow
    pos_c = ax_c.get_position()
    dy_d = TOP_Y[1] - TOP_Y[0]
    dx_d = TOP_X[1] - TOP_X[0]
    w_d = 1.70 / FIG_W   # ~1.70 in wide (name sits left of the arrow)
    h_d = w_d * (dy_d / dx_d) * (FIG_W / FIG_H)
    left_d = 1.0 - 0.015 - w_d
    # one inset-height below a fully-raised seat, so (d) overlaps (c) again
    bot_d = pos_c.y1 - 0.55 * h_d
    bot_d = min(bot_d, 1.0 - 0.015 - h_d)
    ax_d = fig.add_axes([left_d, bot_d, w_d, h_d])
    ax_d.set_zorder(10)
    ax_d.patch.set_zorder(10)

    draw_model(ax_a, data, nodes, xlim=(-x_half_a, x_half_a),
               lw=0.5, ms=0.0, quad_lw=0.15, coils=False, dashpots=True,
               fill=fill, force_at="leftmost",
               exp_element=True, exp_force_fp=True, exp_scale=10.0,
               exp_force_fp_label="full")
    ax_a.set_ylim(*ylim_a)
    mark_window(ax_a, (-x_half_b, x_half_b), ylim_b, "(b)")
    ax_a.set_title("(a)", loc="left", pad=3)

    draw_model(ax_b, data, nodes, xlim=(-x_half_b, x_half_b),
               lw=0.8, ms=1.6, quad_lw=0.25, coils=True, dashpots=True,
               fill=fill, pin_s=0.28, force_label=True, force_under=True,
               exp_element=True, exp_force_fp=True, exp_scale=2.5,
               exp_force_fp_label="none")
    ax_b.set_ylim(*ylim_b)
    mark_window(ax_b, DETAIL_X, DETAIL_Y, "(c)")
    mark_window(ax_b, TOP_X, TOP_Y, "(d)", corner="br")
    ax_b.set_title("(b)", loc="left", pad=3)

    draw_model(ax_c, data, nodes, xlim=DETAIL_X,
               lw=0.8, ms=3.8, quad_lw=0.5, coils=True, dashpots=False,
               fill=fill, soil_ms=2.4, surcharge=False, pin_s=0.18,
               exp_element=False)
    ax_c.set_ylim(*DETAIL_Y)
    ax_c.set_title("(c)", loc="left", pad=3)

    draw_model(ax_d, data, nodes, xlim=TOP_X,
               lw=1.2, ms=4.0, quad_lw=0.3, coils=False, dashpots=False,
               fill=fill, soil_ms=0.0, surcharge=False, pin_s=0.0,
               exp_element=True, exp_label=False, exp_scale=1.05, exp_lw=1.5,
               exp_force_fp=True, exp_force_fp_label="two", equal_aspect=True)
    ax_d.set_ylim(*TOP_Y)
    ax_d.set_facecolor("white")
    ax_d.set_xticks([])
    ax_d.set_yticks([])
    ax_d.set_xlabel("")
    ax_d.set_ylabel("")
    ax_d.set_title("(d)", loc="left", pad=0.5, fontsize=FONT_SIZE - 1)
    for spine in ax_d.spines.values():
        spine.set_color("#222222")
        spine.set_linewidth(0.9)

    for a in (ax_a, ax_b, ax_c):
        a.set_xlabel(r"$x$ (m)", labelpad=1)
        a.set_ylabel(r"$y$ (m)", labelpad=1)
        a.tick_params(pad=1.5)
    ax_a.set_xticks(range(-75, 76, 25))
    ax_b.set_xticks([-5, 0, 5])
    ax_b.set_yticks(range(-20, 11, 5))
    ax_c.set_xticks(range(-3, 4, 1))
    ax_c.set_yticks(range(-2, 2, 1))

    top_leg = top_b + h_c + GAP_XLABEL
    handles, labels = legend_handles(data, fill)
    # 2 columns: ~8 rows in ~1.5 in under (c); pack to the page bottom
    fig.legend(
        handles, labels,
        loc="upper center",
        bbox_to_anchor=((left_c + 0.5 * w_c) / FIG_W, 1.0 - top_leg / FIG_H),
        bbox_transform=fig.transFigure,
        ncol=2, frameon=False, borderaxespad=0.0,
        handlelength=1.3, handleheight=0.85, handletextpad=0.35,
        columnspacing=1.0, labelspacing=0.22,
        handler_map=LEGEND_HANDLERS,
    )

    out_stem.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".png", ".pdf", ".svg"):
        fig.savefig(out_stem.with_suffix(suffix))
    plt.close(fig)
    print(f"PlotModelSketchPaper: wrote {out_stem}.png / .pdf / .svg  ({FIG_W} x {FIG_H} in)")


# ------------------------------------------------------------
# 5. COMMAND-LINE ENTRY POINT
# ------------------------------------------------------------


def main() -> int:
    """
    Resolve CLI paths and write the paper elevation.

    Args:    sys.argv  [--fill color|gray|none] [in.json] [out_stem]
    Returns: process status code
    """
    args = sys.argv[1:]
    fill = FILL_DEFAULT
    if "--fill" in args:
        i = args.index("--fill")
        fill = args[i + 1]
        del args[i:i + 2]
    if fill not in FILL_MODES:
        print(f"--fill must be one of {FILL_MODES} (got '{fill}')", file=sys.stderr)
        return 1
    json_path = Path(args[0]) if args else DEFAULT_JSON
    if not json_path.is_file():
        print(f"missing {json_path}; run OpenSees PlotModel.tcl first", file=sys.stderr)
        return 1
    data = load(json_path)
    if len(args) > 1:
        out_stem = Path(args[1]).with_suffix("")
    else:
        stem = "elevation_paper" if fill == FILL_DEFAULT else f"elevation_paper_{fill}"
        out_stem = elevation_dir(
            data.get("soilProfile") or 1, data.get("soilBoundary") or "Shin",
        ) / stem
    plot(data, out_stem, fill)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
