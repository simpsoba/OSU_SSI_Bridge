#!/usr/bin/env python3
"""
Goals
-----
Show the PT39 solid-concrete outline and its stiff elastic frame idealization.
Annotate the translational mass and rotational inertia assigned at each node.

  python3 plot/PlotDeckFrameConcept.py

Writes:
  plot/out/deck/frame_concept.png

Prefers mass / geometry from plot/out/deck/deck_frame.json (BuildDeckNodes.tcl).
Falls back to Parameters-equivalent dims + length-weighted lumping.

Units
-----
Geometry is in m, nodal mass in kg, and rotational inertia in kg.m².
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon

from pt39_outline import pt39_outline

HERE = Path(__file__).resolve().parent
OUT = HERE / "out" / "deck" / "frame_concept.png"
JSON_PATH = HERE / "out" / "deck" / "deck_frame.json"

# ------------------------------------------------------------
# 1. PT39 GEOMETRY AND FRAME STYLE
# ------------------------------------------------------------

# Nodal inertia: mx, my (kg), Irot (kg·m²)
MassTriple = tuple[float, float, float]

# --- PT39 geometry fallback (Parameters.tcl / Deck_PT_39), m ---
foot = 0.3048
inch = 0.0254
DW = 39.0 * foot
DD = 6.0 * foot
SW = 23.0 * foot
CW = 5.5 * foot
TD = 9.5 * inch
TS = 8.0 * inch
TW = 12.0 * inch
BH = 32.0 * inch
YB = 43.6 * inch
A_DECK = 8869.0 * 6.4516e-4
DENS_DECK = 22.78e3 / 9.81
L_TRIB = 150.0 * foot
M_DECK = DENS_DECK * A_DECK * L_TRIB

MEMBERS_KIND = [
    ("BL", "BC", "soffit"),
    ("BC", "BR", "soffit"),
    ("TL", "TLi", "top"),
    ("TLi", "TC", "top"),
    ("TC", "TRi", "top"),
    ("TRi", "TR", "top"),
    ("BL", "TLi", "web"),
    ("BR", "TRi", "web"),
    ("BC", "TC", "center"),
    ("TL", "BarL", "barrier"),
    ("TR", "BarR", "barrier"),
]

STYLE = {
    "top": {"color": "#1565c0", "lw": 2.4},
    "soffit": {"color": "#2e7d32", "lw": 2.4},
    "web": {"color": "#6a1b9a", "lw": 2.2},
    "center": {"color": "#00838f", "lw": 2.0, "ls": "--"},
    "barrier": {"color": "#ef6c00", "lw": 1.8},
    "link": {"color": "#c62828", "lw": 2.8},
}


# ------------------------------------------------------------
# 2. FRAME GEOMETRY AND MASS
# ------------------------------------------------------------


def frame_nodes(dw: float, dd: float, sw: float, cw: float, bh: float, y0: float = 0.0):
    """
    Build named frame-node coordinates from PT39 dimensions.

    Args:    dw, dd, sw, cw, bh, y0  (m)
    Returns: node name → (x, y), m
    """
    x_top_outer = 0.5 * dw
    x_overhang_in = 0.5 * dw - cw
    x_soffit = 0.5 * sw
    y_top = y0 + dd
    y_bar = y_top + bh
    return {
        "BL": (-x_soffit, y0),
        "BC": (0.0, y0),
        "BR": (x_soffit, y0),
        "TL": (-x_top_outer, y_top),
        "TLi": (-x_overhang_in, y_top),
        "TC": (0.0, y_top),
        "TRi": (x_overhang_in, y_top),
        "TR": (x_top_outer, y_top),
        "BarL": (-x_top_outer, y_bar),
        "BarR": (x_top_outer, y_bar),
    }


def length_weighted_mass(
    nodes: dict[str, tuple[float, float]],
    members: list[tuple[str, str, str]],
    m_total: float,
    dw: float,
    dd: float,
    yb: float,
) -> dict[str, MassTriple]:
    """
    Reproduce BuildDeck length-weighted mass and rotational inertia.

    Args:    nodes, members, m_total  (kg), dw, dd, yb  (m)
    Returns: node name → (mx, my, Irot)
    """
    names = list(nodes.keys())
    w = {nm: 0.0 for nm in names}
    mem_L: list[tuple[str, str, float]] = []
    for a, b, _k in members:
        xa, ya = nodes[a]
        xb, yb_ = nodes[b]
        L = math.hypot(xb - xa, yb_ - ya)
        w[a] += 0.5 * L
        w[b] += 0.5 * L
        mem_L.append((a, b, L))
    wsum = sum(w.values())
    if wsum <= 0.0:
        raise RuntimeError("zero member-length weight sum")
    m = {nm: m_total * w[nm] / wsum for nm in names}
    irot = {nm: 0.0 for nm in names}
    for a, b, L in mem_L:
        m_mem = m_total * L / wsum
        i_end = m_mem * L * L / 105.0
        irot[a] += i_end
        irot[b] += i_end
    return {nm: (m[nm], m[nm], irot[nm]) for nm in names}


# ------------------------------------------------------------
# 3. EXPORTED DATA OR PARAMETER FALLBACK
# ------------------------------------------------------------


def load_or_fallback() -> tuple[
    dict[str, tuple[float, float]],
    list[tuple[str, str, str]],
    dict[str, MassTriple],
    dict,
    str,
]:
    """
    Read BuildDeck output, or reconstruct the same conceptual frame.

    Args:    none
    Returns: (nodes, members, nodal inertia, geometry, source label)
    """
    if JSON_PATH.is_file():
        with JSON_PATH.open() as f:
            data = json.load(f)
        y0 = float(data["y0"])
        nodes = {}
        mass: dict[str, MassTriple] = {}
        for n in data["nodes"]:
            nm = n["name"]
            nodes[nm] = (float(n["x"]), float(n["y"]) - y0)
            if "mx" in n and "my" in n:
                mx, my = float(n["mx"]), float(n["my"])
            else:
                mx = my = float(n["m"])
            irot = float(n.get("Irot", 0.0))
            mass[nm] = (mx, my, irot)
        members = []
        for m in data["members"]:
            kind = m["kind"]
            if kind == "center":
                members.append((m["i"], m["j"], "center"))
            elif m["i"].startswith("Bar") or m["j"].startswith("Bar"):
                members.append((m["i"], m["j"], "barrier"))
            elif float(nodes[m["i"]][1]) < 1e-9 and float(nodes[m["j"]][1]) < 1e-9:
                members.append((m["i"], m["j"], "soffit"))
            elif abs(float(nodes[m["i"]][1]) - float(nodes[m["j"]][1])) < 1e-9:
                members.append((m["i"], m["j"], "top"))
            else:
                members.append((m["i"], m["j"], "web"))
        geom = {
            "dw": float(data["dw"]),
            "dd": float(data["dd"]),
            "sw": float(data["sw"]),
            "cw": float(data["cw"]),
            "td": float(data["td"]),
            "ts": float(data["ts"]),
            "tw": float(data["tw"]),
            "bh": float(data["bh"]),
            "yb": float(data["yb"]),
            "m_deck": float(data["m_deck"]),
            "Iz_deck": float(data.get("Iz_deck", 0.0)),
            "I_steiner": float(data.get("I_steiner", 0.0)),
            "Irot_fill": float(data.get("Irot_fill", 0.0)),
        }
        return nodes, members, mass, geom, str(JSON_PATH)

    geom = {
        "dw": DW,
        "dd": DD,
        "sw": SW,
        "cw": CW,
        "td": TD,
        "ts": TS,
        "tw": TW,
        "bh": BH,
        "yb": YB,
        "m_deck": M_DECK,
        "Iz_deck": 0.0,
        "I_steiner": 0.0,
        "Irot_fill": 0.0,
    }
    nodes = frame_nodes(DW, DD, SW, CW, BH, y0=0.0)
    members = list(MEMBERS_KIND)
    mass = length_weighted_mass(nodes, members, M_DECK, DW, DD, YB)
    return nodes, members, mass, geom, "Parameters fallback"


def _fmt_mass_label(name: str, mx: float, my: float, irot: float) -> str:
    """
    Format one compact nodal-inertia annotation.

    Args:    name, mx, my  (kg), irot  (kg.m²)
    Returns: multiline label
    """
    return (
        f"{name}\n"
        f"$m_x$={mx:,.0f} kg\n"
        f"$m_y$={my:,.0f} kg\n"
        f"$I$={irot:,.0f} kg·m$^2$"
    )


# ------------------------------------------------------------
# 4. FRAME-CONCEPT FIGURE
# ------------------------------------------------------------


def main() -> int:
    """
    Draw the frame concept and write frame_concept.png.

    Args:    none
    Returns: process status code
    """
    nodes, members, mass, geom, src = load_or_fallback()
    dw, dd, sw, cw = geom["dw"], geom["dd"], geom["sw"], geom["cw"]
    td, ts, tw, bh, yb = geom["td"], geom["ts"], geom["tw"], geom["bh"], geom["yb"]
    m_deck = geom["m_deck"]

    # Pier stub for drawing only (not in mass dump)
    pier_y = -0.6 * dd
    draw_nodes = dict(nodes)
    draw_nodes["Pier"] = (0.0, pier_y)
    draw_members = list(members) + [("Pier", "BC", "link")]

    fig, ax = plt.subplots(figsize=(12.5, 6.2), constrained_layout=True)

    for poly in pt39_outline(dw, dd, sw, cw, td, ts, tw):
        poly.set_facecolor("#90a4ae")
        poly.set_edgecolor("#546e7a")
        poly.set_alpha(0.28)
        poly.set_linewidth(0.8)
        ax.add_patch(poly)

    for a, b, kind in draw_members:
        xa, ya = draw_nodes[a]
        xb, yb_ = draw_nodes[b]
        st = STYLE[kind]
        ax.plot(
            [xa, xb],
            [ya, yb_],
            color=st["color"],
            lw=st["lw"],
            ls=st.get("ls", "-"),
            solid_capstyle="round",
            zorder=3,
        )

    for name, (x, y) in draw_nodes.items():
        if name.startswith("Bar"):
            ms, mec, mfc = 6, "#ef6c00", "white"
        elif name == "Pier":
            ms, mec, mfc = 9, "#c62828", "#c62828"
        elif name in ("BC", "TC"):
            ms, mec, mfc = 8, "#00838f", "white"
        else:
            ms, mec, mfc = 7, "#37474f", "white"
        ax.plot(x, y, "o", ms=ms, mec=mec, mfc=mfc, mew=1.3, zorder=4)

    # Nodal inertia: mx, my (kg), I (kg·m²)
    for name, (x, y) in nodes.items():
        mx, my, irot = mass[name]
        dx = -0.18 if x < -1e-9 else (0.18 if x > 1e-9 else 0.0)
        ha = "right" if x < -1e-9 else ("left" if x > 1e-9 else "center")
        if name.startswith("Bar"):
            ay, va = y + 0.08, "bottom"
        elif y < 1e-9:
            ay, va = y - 0.08, "top"
        else:
            ay, va = y + 0.08, "bottom"
        ax.annotate(
            _fmt_mass_label(name, mx, my, irot),
            xy=(x, y),
            xytext=(x + dx, ay),
            fontsize=6.0,
            ha=ha,
            va=va,
            color="#37474f",
            zorder=6,
            linespacing=1.15,
        )

    ax.annotate(
        "pier top\n(node 5)",
        xy=(0.0, pier_y),
        xytext=(0.0, pier_y - 0.12),
        fontsize=7,
        ha="center",
        va="top",
        color="#37474f",
    )

    ax.plot(0.0, yb, "+", ms=12, mew=1.5, color="#6d4c41", zorder=5)
    ax.text(
        0.15,
        yb,
        f"CG  Σm={m_deck:,.0f} kg",
        fontsize=7.5,
        color="#6d4c41",
        va="center",
    )

    ax.set_aspect("equal")
    ax.set_xlabel("Transverse x (m)")
    ax.set_ylabel("Elevation y (m), soffit bottom = 0")
    ax.set_xlim(-0.55 * dw - 1.6, 0.55 * dw + 1.8)
    ax.set_ylim(-1.05 * dd, dd + bh + 0.95)
    ax.axhline(0.0, color="#bdbdbd", lw=0.6, zorder=0)
    ax.axvline(0.0, color="#bdbdbd", lw=0.6, ls=":", zorder=0)
    ax.grid(True, ls=":", alpha=0.35)

    handles = [
        mpatches.Patch(
            facecolor="#90a4ae", alpha=0.35, edgecolor="#546e7a",
            label="PT39 solid outline",
        ),
        Line2D([0], [0], color="#1565c0", lw=2.4, label="top slab chord"),
        Line2D([0], [0], color="#2e7d32", lw=2.4, label="soffit chord"),
        Line2D([0], [0], color="#6a1b9a", lw=2.2, label="web / strut"),
        Line2D([0], [0], color="#00838f", lw=2.0, ls="--", label="center web"),
        Line2D([0], [0], color="#ef6c00", lw=1.8, label="barrier"),
        Line2D([0], [0], color="#c62828", lw=2.8, label="pier–deck equalDOF"),
    ]
    ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=4,
        fontsize=8,
        framealpha=0.92,
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"PlotDeckFrameConcept: wrote {OUT}  (mass source: {src})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
