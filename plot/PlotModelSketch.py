#!/usr/bin/env python3
"""Elevation from plot/model_sketch.json (DumpModelSketch.tcl).

  OpenSees PlotModel.tcl
  python3 plot/PlotModelSketch.py [in.json] [out.png]
  # default → plot/out/profile{N}/elevation/{Shin|ASDEA}/elevation.png
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, Rectangle

from paths import HERE, elevation_png
from pt39_outline import pt39_outline

DEFAULT_JSON = HERE / "model_sketch.json"

SPRING_GAP = 0.70
# Visual length for coincident / short SSI springs (like rot. spiral gap)
SSI_VIS_LEN = 0.55

STYLE = {
    "pier": {"fill": "#f0a060", "line": "#c45c12", "node": "#c45c12", "alpha": 0.55},
    "deck": {"fill": "#90caf9", "line": "#1565c0", "node": "#0d47a1", "alpha": 0.30},
    "cap": {"fill": "#9e9e9e", "line": "#333333", "node": "#222222", "alpha": 0.45},
    "pile": {"fill": "#c4a484", "line": "#8B5A2B", "node": "#8B5A2B", "alpha": 0.50},
    "spring": {"line": "#6a1b9a"},
    "asdea": {"fill": "#ef9a9a", "line": "#c62828", "alpha": 0.35},
    "lysmer": {"line": "#ad1457"},
}

# One color per uniaxial SSI material (dir 1 lateral, dir 2 axial)
SSI_STYLE = {
    "py": {"line": "#1565c0", "label": "PySimple1"},
    "pyliq": {"line": "#00acc1", "label": "PyLiq1"},
    "py_elastic": {"line": "#90caf9", "label": "p-y elastic"},
    "tz": {"line": "#2e7d32", "label": "TzSimple1"},
    "tzliq": {"line": "#7b1fa2", "label": "TzLiq1"},
    "tz_elastic": {"line": "#a5d6a7", "label": "t-z elastic"},
    "qz": {"line": "#e65100", "label": "QzSimple1"},
    "qz_elastic": {"line": "#ffcc80", "label": "q-z elastic"},
    "none": {"line": "#9e9e9e", "label": "direct equalDOF"},
}
SSI_TYPE_ORDER = (
    "py", "pyliq", "py_elastic",
    "tz", "tzliq", "tz_elastic",
    "qz", "qz_elastic", "none",
)
# Axial (vertical) coils slightly offset from lateral so both read at one node
SSI_AXIAL_DX = 0.12


# Distinct colors per soil layer name
LAYER_STYLE = {
    "L2": {"fill": "#7e9bb5", "label": "L2 clay"},
    "L3": {"fill": "#6b8f71", "label": "L3 clay"},
    "L3a": {"fill": "#e8c547", "label": "L3 sand"},
    "L3b": {"fill": "#d4a017", "label": "L3 sand"},
    "L3c": {"fill": "#b8860b", "label": "L3 sand"},
    "L5": {"fill": "#8B4513", "label": "L5"},
}


def load(path: Path) -> dict:
    with path.open() as f:
        return json.load(f)


def node_map(data: dict) -> dict[int, tuple[float, float]]:
    return {int(t): (float(x), float(y)) for t, x, y in data["nodes"]}


def display_nodes(
    nodes: dict[int, tuple[float, float]],
    springs: list,
) -> dict[int, tuple[float, float]]:
    """Spread coincident lumpedPlasticity ZLS so both spirals read in elevation.

    Base: 1 (cap TC) stays, 2 (inner) up. Top: 5 (deck BC) stays, 4 (inner) down.
    """
    disp = dict(nodes)
    if not springs:
        return disp
    if 1 in disp and 2 in disp:
        x1, y1 = nodes[1]
        disp[1] = (x1, y1)
        disp[2] = (x1, y1 + SPRING_GAP)
    if 4 in disp and 5 in disp:
        x5, y5 = nodes[5]
        disp[5] = (x5, y5)
        disp[4] = (x5, y5 - SPRING_GAP)
    return disp


def group_of_node(tag: int) -> str:
    if tag < 1000:
        return "pier"
    if tag < 2000:
        return "cap"
    if tag < 3000:
        return "pile"
    if tag < 4000:
        return "deck"
    if tag < 10000:
        return "pile"
    if tag < 21000:
        return "ssi_spring"
    return "ssi_spring"


def layer_style(name: str, profile: int | None = None) -> dict:
    if profile == 1 and name == "L3":
        return {"fill": "#e8c547", "label": "L3 sand"}
    if name in LAYER_STYLE:
        return LAYER_STYLE[name]
    if name.startswith("L5"):
        return LAYER_STYLE["L5"]
    return {"fill": "#a0a0a0", "label": name}


def draw_rot_spiral(ax, x0: float, y0: float, x1: float, y1: float, color: str) -> None:
    xm = 0.5 * (x0 + x1)
    ym = 0.5 * (y0 + y1)
    gap = abs(y1 - y0)
    r_max = min(0.20, 0.32 * gap)
    n_turns = 2.0
    n = 90
    theta = np.linspace(0.0, 2.0 * np.pi * n_turns, n)
    r = 0.05 + (r_max - 0.05) * (theta / theta[-1])
    xs = xm + r * np.cos(theta)
    ys = ym + r * np.sin(theta)
    ax.plot([x0, xs[0]], [y0, ys[0]], color=color, lw=1.1, zorder=6)
    ax.plot([x1, xs[-1]], [y1, ys[-1]], color=color, lw=1.1, zorder=6)
    ax.plot(xs, ys, color=color, lw=1.35, zorder=6, solid_capstyle="round")


def draw_trans_coil(
    ax,
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    color: str,
) -> None:
    """Zigzag coil between two points (assumed visual length already applied)."""
    dx = x1 - x0
    dy = y1 - y0
    L = float(np.hypot(dx, dy))
    if L < 1.0e-9:
        return
    nx, ny = -dy / L, dx / L
    n_zig = 7
    amp = min(0.12, 0.25 * L)
    xs = [x0]
    ys = [y0]
    for i in range(1, n_zig + 1):
        t = i / (n_zig + 1)
        side = 1.0 if (i % 2) else -1.0
        xs.append(x0 + t * dx + side * amp * nx)
        ys.append(y0 + t * dy + side * amp * ny)
    xs.append(x1)
    ys.append(y1)
    ax.plot(xs, ys, color=color, lw=1.05, zorder=6, solid_capstyle="round")


def ssi_coil_ends(
    stype: str,
    xp: float,
    yp: float,
    xi: float,
    yi: float,
    kind: str = "pile",
) -> tuple[float, float, float, float]:
    """Place lateral coils horizontal at structure elev; axial vertical, mirrored for qz."""
    lat = stype in ("py", "pyliq", "py_elastic")
    dx = xi - xp
    # Outward from pile/cap center: iface−pile, else sign of x
    if abs(dx) > 1.0e-12:
        sign = 1.0 if dx > 0 else -1.0
    else:
        sign = 1.0 if (xi if abs(xi) > 1e-12 else xp) >= 0 else -1.0
    # Draw at structure elevation (yp) so cap mid/top/bot stay distinct
    y_draw = yp
    if lat:
        # Cap: always use visual length at face (soil may be coincident)
        if kind == "cap" or abs(dx) < 0.12:
            xm = xp if kind == "cap" else 0.5 * (xp + xi)
            return (
                xm - 0.5 * SSI_VIS_LEN * sign,
                y_draw,
                xm + 0.5 * SSI_VIS_LEN * sign,
                y_draw,
            )
        return xp, yp, xi, yi
    # Axial: tz/tzliq outward; qz mirrored inward so tips don't overlap
    if stype in ("qz", "qz_elastic"):
        xo = xi - sign * SSI_AXIAL_DX
    else:
        xo = xi + sign * SSI_AXIAL_DX
    return xo, y_draw - 0.5 * SSI_VIS_LEN, xo, y_draw + 0.5 * SSI_VIS_LEN


def parse_ssi_row(row: list) -> tuple[float, float, float, float, str, str]:
    """Return (xp, yp, xi, yi, stype, kind) from dump row (new or legacy)."""
    if len(row) >= 9:
        _e, _xs, _ys, xi, yi, xp, yp, kind, stype = row[:9]
        return float(xp), float(yp), float(xi), float(yi), str(stype), str(kind)
    if len(row) >= 8:
        _e, _xs, _ys, xi, yi, xp, yp, kind = row[:8]
        return float(xp), float(yp), float(xi), float(yi), "py", str(kind)
    x, y = float(row[1]), float(row[2])
    return x, y, x, y, "py", "pile"


def draw_lysmer_mark(ax, x: float, y: float, color: str) -> None:
    """Horizontal dashpot glyph at (x, y); sized for full-domain view."""
    w, h = 4.0, 1.2
    ax.plot([x - w, x - 0.7], [y, y], color=color, lw=1.6, zorder=7)
    ax.plot([x + 0.7, x + w], [y, y], color=color, lw=1.6, zorder=7)
    ax.add_patch(Rectangle(
        (x - 0.7, y - 0.5 * h), 1.4, h,
        fill=False, edgecolor=color, lw=1.5, zorder=7,
    ))
    ax.plot([x - 0.3, x - 0.3], [y - 0.35 * h, y + 0.35 * h],
            color=color, lw=1.3, zorder=8)
    ax.plot([x + 0.3, x + 0.3], [y - 0.35 * h, y + 0.35 * h],
            color=color, lw=1.3, zorder=8)


def plot(data: dict, out: Path) -> None:
    plt.rcParams["mathtext.fontset"] = "cm"
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.size"] = 9

    springs = data.get("springs", [])
    ssi = data.get("ssi_springs", [])
    quads = data.get("soil_quads", [])
    bnd_quads = data.get("bnd_quads", [])
    lysmer = data.get("lysmer_dashpots", [])
    layers = data.get("soil_layers", [])
    nodes = display_nodes(node_map(data), springs)
    sz = data["sizes"]
    H_pier = float(sz["H_pier"])
    W_cap = float(sz["W_cap"])
    H_cap = float(sz["H_cap"])
    L_pile = float(sz["L_pile"])
    nrow = int(sz.get("n_pile_row", 2))
    pier_type = data.get("pierEleType", "")
    pile_type = data.get("pileEleType", "")
    profile = data.get("soilProfile")
    boundary = data.get("soilBoundary") or ""

    out.parent.mkdir(parents=True, exist_ok=True)

    # Two panels if soil present: near-field + full domain
    has_soil = bool(quads)
    if has_soil:
        fig, (ax, axf) = plt.subplots(
            1, 2, figsize=(11.0, 9.5), dpi=160,
            gridspec_kw={"width_ratios": [1.0, 1.35]},
        )
        axes = [ax, axf]
    else:
        fig, ax = plt.subplots(figsize=(5.8, 10.2), dpi=160)
        axes = [ax]

    def draw_on(a: plt.Axes, *, full: bool) -> None:
        # Soil quads under structure
        if quads:
            by_layer: dict[str, list] = {}
            for q in quads:
                by_layer.setdefault(q["layer"], []).append(q)
            for nm, qs in by_layer.items():
                st = layer_style(nm, profile=profile)
                polys = []
                for q in qs:
                    xy = q["xy"]
                    pts = [(xy[i], xy[i + 1]) for i in range(0, 8, 2)]
                    polys.append(Polygon(pts, closed=True))
                coll = PatchCollection(
                    polys, facecolor=st["fill"], edgecolor="#333333",
                    linewidths=0.15, alpha=0.55, zorder=0,
                )
                a.add_collection(coll)

        # ASDEA ring (full-domain panel)
        if bnd_quads and full:
            st = STYLE["asdea"]
            polys = []
            for q in bnd_quads:
                xy = q["xy"]
                pts = [(xy[i], xy[i + 1]) for i in range(0, 8, 2)]
                polys.append(Polygon(pts, closed=True))
            a.add_collection(PatchCollection(
                polys, facecolor=st["fill"], edgecolor=st["line"],
                linewidths=0.45, alpha=st["alpha"], zorder=1,
                hatch="///",
            ))

        for fill in data["fills"]:
            g = fill["group"]
            if g == "deck":
                continue  # PT39 outline patches below
            st = STYLE[g]
            xc, w = float(fill["xc"]), float(fill["width"])
            y0, y1 = float(fill["y0"]), float(fill["y1"])
            if g == "pier" and springs and 2 in nodes:
                y0 = nodes[2][1]
                y1 = nodes[4][1] if 4 in nodes else y1
            a.add_patch(Rectangle(
                (xc - 0.5 * w, min(y0, y1)),
                w, abs(y1 - y0),
                facecolor=st["fill"], edgecolor="none",
                alpha=st["alpha"], zorder=2,
            ))

        # PT39 solid concrete outline (same patches as frame_concept)
        if "dw_deck" in sz and "dd_deck" in sz:
            dw = float(sz["dw_deck"])
            dd = float(sz["dd_deck"])
            sw = float(sz.get("sw_deck", 23.0 * 0.3048))
            cw = float(sz.get("cw_deck", 5.5 * 0.3048))
            td = float(sz.get("td_deck", 9.5 * 0.0254))
            ts = float(sz.get("ts_deck", 8.0 * 0.0254))
            tw = float(sz.get("tw_deck", 12.0 * 0.0254))
            y_soffit = H_pier
            st = STYLE["deck"]
            for poly in pt39_outline(dw, dd, sw, cw, td, ts, tw, y0=y_soffit):
                poly.set_facecolor(st["fill"])
                poly.set_edgecolor(st["line"])
                poly.set_alpha(st["alpha"])
                poly.set_linewidth(0.6)
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
            a.add_collection(LineCollection(
                segs, colors=colors, linewidths=1.35, zorder=3,
            ))

        for spr in springs:
            _e, ni, nj, _sx, _sy = spr
            if int(ni) not in nodes or int(nj) not in nodes:
                continue
            x0, y0 = nodes[int(ni)]
            x1, y1 = nodes[int(nj)]
            draw_rot_spiral(a, x0, y0, x1, y1, STYLE["spring"]["line"])

        if ssi and not full:
            for row in ssi:
                xp, yp, xi, yi, stype, kind = parse_ssi_row(row)
                if stype == "none":
                    continue
                st = SSI_STYLE.get(stype, SSI_STYLE["py"])
                x0, y0, x1, y1 = ssi_coil_ends(stype, xp, yp, xi, yi, kind)
                draw_trans_coil(a, x0, y0, x1, y1, st["line"])

        # Shin Lysmer at base (both panels; glyphs scale with full view)
        if lysmer and full:
            for d in lysmer:
                draw_lysmer_mark(
                    a, float(d["x"]), float(d["y"]), STYLE["lysmer"]["line"],
                )

        ms = 4.0
        for tag, (x, y) in nodes.items():
            if tag >= 20000:
                continue
            g = group_of_node(tag)
            if g not in STYLE:
                continue
            a.plot(
                x, y, "o", color=STYLE[g]["node"], ms=ms, zorder=5,
                markeredgecolor="white", markeredgewidth=0.35,
            )

        a.set_aspect("equal")
        a.set_xlabel(r"$x$ (m)")
        a.set_ylabel(r"$y$ (m)")
        a.grid(True, alpha=0.25, lw=0.5)

        y_bot = -H_cap - L_pile - 0.8
        if layers:
            y_bot = min(y_bot, min(float(L["y0"]) for L in layers) - 0.5)
        if bnd_quads:
            for q in bnd_quads:
                xy = q["xy"]
                ys = [xy[i + 1] for i in range(0, 8, 2)]
                y_bot = min(y_bot, min(ys) - 0.5)
        y_top = max(H_pier, max((y for _, y in nodes.values()), default=H_pier)) + 0.8

        if full and ("xMeshHalf" in sz or "L_half" in sz):
            Lh = float(sz.get("xMeshHalf", sz["L_half"]))
            pad = float(sz.get("w_FF", 0.0) or 0.0) + 2.0
            a.set_xlim(-Lh - pad, Lh + pad)
            ttl = "full soil domain"
            if boundary:
                ttl += f"  ({boundary})"
            a.set_title(ttl)
        else:
            half = 0.5 * W_cap + 3.2
            if "dw_deck" in sz:
                half = max(half, 0.5 * float(sz["dw_deck"]) + 0.8)
            a.set_xlim(-half, half)
            ttl = rf"pier={pier_type}  pile={pile_type}"
            if profile is not None:
                ttl += rf"  soil={profile}"
            if boundary:
                ttl += rf"  {boundary}"
            a.set_title(ttl)
        a.set_ylim(y_bot, y_top)

    draw_on(axes[0], full=False)
    if has_soil:
        draw_on(axes[1], full=True)

    handles = [
        mpatches.Patch(
            facecolor=STYLE["pier"]["fill"], alpha=0.55,
            edgecolor=STYLE["pier"]["line"], label="pier",
        ),
    ]
    if any(f.get("group") == "deck" for f in data.get("fills", [])) or any(
        e[3] == "deck" for e in data.get("elements", [])
    ) or "dw_deck" in sz:
        handles.append(
            mpatches.Patch(
                facecolor=STYLE["deck"]["fill"], alpha=0.35,
                edgecolor=STYLE["deck"]["line"], label="deck (PT39 outline)",
            )
        )
    handles.extend([
        mpatches.Patch(
            facecolor=STYLE["cap"]["fill"], alpha=0.45,
            edgecolor=STYLE["cap"]["line"], label="pile cap",
        ),
        mpatches.Patch(
            facecolor=STYLE["pile"]["fill"], alpha=0.5,
            edgecolor=STYLE["pile"]["line"],
            label=rf"piles ($\times {nrow}$)",
        ),
    ])
    seen = set()
    for L in layers:
        nm = L["name"]
        if nm in seen:
            continue
        seen.add(nm)
        st = layer_style(nm)
        lab = st["label"]
        if profile == 4 and nm == "L2":
            lab = "L2 soft clay"
        elif nm == "L3" and profile == 1:
            lab = "L3 sand"
        elif nm == "L3":
            lab = "L3 clay"
        if nm == "L5" and profile in (3, 4):
            lab = "L5 stiff clay"
        elif nm == "L5":
            lab = "L5 dense sand"
        handles.append(
            mpatches.Patch(facecolor=st["fill"], alpha=0.55, label=lab)
        )
    if springs:
        handles.append(
            Line2D(
                [0], [0], color=STYLE["spring"]["line"], lw=1.5,
                label="rot. spring (ZLS)",
            )
        )
    if ssi:
        present = {parse_ssi_row(r)[4] for r in ssi}
        for key in SSI_TYPE_ORDER:
            if key not in present:
                continue
            st = SSI_STYLE[key]
            handles.append(
                Line2D([0], [0], color=st["line"], lw=1.2, label=st["label"])
            )
    if bnd_quads:
        handles.append(
            mpatches.Patch(
                facecolor=STYLE["asdea"]["fill"], alpha=0.45,
                edgecolor=STYLE["asdea"]["line"], hatch="///",
                label="ASDEA boundary",
            )
        )
    if lysmer:
        handles.append(
            Line2D(
                [0], [0], color=STYLE["lysmer"]["line"], lw=1.5,
                label="Lysmer dashpot",
            )
        )
    fig.tight_layout()
    if has_soil:
        # One block above the full-domain panel; legend top = left axes top.
        pos_L = axes[0].get_position()
        pos_R = axes[1].get_position()
        xc = pos_R.x0 + 0.5 * pos_R.width
        fig.legend(
            handles=handles,
            loc="upper center",
            bbox_to_anchor=(xc, pos_L.y1),
            bbox_transform=fig.transFigure,
            ncol=3,
            fontsize=7.0,
            frameon=False,
            borderaxespad=0.0,
        )
    else:
        axes[0].legend(
            handles=handles, loc="upper left", fontsize=7.5, frameon=False,
        )
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(
        f"PlotModelSketch: wrote {out}  "
        f"({len(springs)} ZLS, {len(ssi)} ssi mats, {len(quads)} quads, "
        f"{len(bnd_quads)} ASDEA, {len(lysmer)} Lysmer)"
    )


def main() -> int:
    json_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_JSON
    if not json_path.is_file():
        print(
            f"missing {json_path}; run OpenSees PlotModel.tcl first",
            file=sys.stderr,
        )
        return 1
    data = load(json_path)
    if len(sys.argv) > 2:
        png_path = Path(sys.argv[2])
    else:
        prof = data.get("soilProfile") or 1
        bnd = data.get("soilBoundary") or "Shin"
        png_path = elevation_png(prof, bnd)
    plot(data, png_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
