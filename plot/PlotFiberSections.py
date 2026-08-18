#!/usr/bin/env python3
"""Fiber section cuts from plot/fiber_sections.json (DumpFiberSections.tcl).

Draws what OpenSees integrates: graded horizontal strips, and (for the pier)
rebar fibers already merged onto z = 0 with combined area — same as
circularRebarYFibers / circularTubeFiberStripsGraded.

  python3 plot/PlotFiberSections.py [in.json] [out_dir]
  # default → plot/out/profile{N}/fibers/  (profile from soil_profile.json if present)
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.patches import Circle, Polygon

from paths import HERE, fibers_dir

DEFAULT_JSON = HERE / "fiber_sections.json"


def load(path: Path) -> dict:
    with path.open() as f:
        return json.load(f)


def y_breaks(Ro: float, n_fibers: int, n_edge: int, edge_frac: float = 1.0 / 6.0) -> list[float]:
    D = 2.0 * Ro
    Ledge = edge_frac * D
    Lmid = (1.0 - 2.0 * edge_frac) * D
    n_mid = n_fibers - 2 * n_edge
    y_bot, y_e1 = -Ro, -Ro + Ledge
    y_e2, y_top = Ro - Ledge, Ro
    br = [y_bot]
    for i in range(1, n_edge + 1):
        br.append(y_bot + Ledge * i / n_edge)
    br[-1] = y_e1
    for i in range(1, n_mid + 1):
        br.append(y_e1 + Lmid * i / n_mid)
    br[-1] = y_e2
    for i in range(1, n_edge + 1):
        br.append(y_e2 + Ledge * i / n_edge)
    br[-1] = y_top
    return br


def chord_half(R: float, y: float) -> float:
    v = R * R - y * y
    return math.sqrt(v) if v > 0.0 else 0.0


def annulus_strip_polys(Ro: float, Ri: float, y1: float, y2: float) -> list[np.ndarray]:
    """Horizontal strip of an annulus → polygons in (z, y)."""
    ys = np.linspace(y1, y2, 8)
    polys: list[np.ndarray] = []

    if abs(y1) >= Ri and abs(y2) >= Ri and y1 * y2 > 0:
        zo1 = chord_half(Ro, y1)
        zo2 = chord_half(Ro, y2)
        polys.append(np.array([
            [-zo1, y1], [zo1, y1], [zo2, y2], [-zo2, y2],
        ]))
        return polys

    for sign in (-1.0, 1.0):
        pts = []
        for y in ys:
            pts.append([sign * chord_half(Ro, y), y])
        for y in ys[::-1]:
            zi = chord_half(Ri, y) if abs(y) <= Ri else 0.0
            pts.append([sign * zi, y])
        polys.append(np.array(pts))
    return polys


def solid_strip_poly(R: float, y1: float, y2: float) -> np.ndarray:
    z1, z2 = chord_half(R, y1), chord_half(R, y2)
    return np.array([[-z1, y1], [z1, y1], [z2, y2], [-z2, y2]])


def plot_pier(sec: dict, out: Path) -> None:
    """RC circle as modeled: graded core/cover strips + rebar fibers on z=0."""
    R = float(sec["R"])
    Rc = float(sec["R_core"])
    nY = int(sec["nFiberY"])
    nE = int(sec["nFiberEdge"])
    ef = float(sec.get("edgeFrac", 1.0 / 6.0))
    br = y_breaks(R, nY, nE, ef)

    plt.rcParams["mathtext.fontset"] = "cm"
    fig, ax = plt.subplots(figsize=(5.2, 5.2), dpi=150)

    for i in range(len(br) - 1):
        y1, y2 = br[i], br[i + 1]
        for poly in annulus_strip_polys(R, Rc, y1, y2):
            ax.add_patch(Polygon(
                poly, closed=True, facecolor="#c4b59a", edgecolor="#8a7a60",
                linewidth=0.4, alpha=0.85, zorder=2,
            ))
        ax.add_patch(Polygon(
            solid_strip_poly(Rc, y1, y2), closed=True,
            facecolor="#9e9e9e", edgecolor="#666", linewidth=0.35, alpha=0.9, zorder=3,
        ))

    # Rebar as modeled: one fiber per unique y on z=0, circle area = fiber A
    for y, z, A in sec["rebar"]:
        r = math.sqrt(A / math.pi)
        ax.add_patch(Circle(
            (z, y), r, facecolor="#c62828", edgecolor="#7f0000",
            linewidth=0.5, zorder=5,
        ))

    ax.add_patch(Circle((0, 0), R, fill=False, ec="#333", lw=1.1, zorder=6))
    ax.add_patch(Circle((0, 0), Rc, fill=False, ec="#333", lw=0.8, ls="--", zorder=6))

    ax.set_aspect("equal")
    lim = 1.15 * R
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel(r"$z$ (m)")
    ax.set_ylabel(r"$y$ (m)")
    ax.set_title(r"Pier Fiber (ZLS / FBC) — as modeled")
    ax.legend(
        handles=[
            mpatches.Patch(facecolor="#9e9e9e", edgecolor="#666", label="core strips"),
            mpatches.Patch(facecolor="#c4b59a", edgecolor="#8a7a60", label="cover strips"),
            mpatches.Patch(
                facecolor="#c62828", edgecolor="#7f0000",
                label=r"rebar fiber ($A$ on $z=0$)",
            ),
        ],
        loc="upper right", fontsize=8, frameon=False,
    )
    ax.grid(True, alpha=0.25, lw=0.5)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"PlotFiberSections: wrote {out}")


def plot_pile(sec: dict, out: Path) -> None:
    """Steel tube as modeled: graded annular strips; fiber at strip centroid."""
    Ro = float(sec["Ro"])
    Ri = float(sec["Ri"])
    nY = int(sec["nFiberY"])
    nE = int(sec["nFiberEdge"])
    ef = float(sec.get("edgeFrac", 1.0 / 6.0))
    nrow = int(sec.get("n_pile_row", 1))
    br = y_breaks(Ro, nY, nE, ef)

    plt.rcParams["mathtext.fontset"] = "cm"
    fig, ax = plt.subplots(figsize=(5.2, 5.2), dpi=150)

    for i in range(len(br) - 1):
        y1, y2 = br[i], br[i + 1]
        for poly in annulus_strip_polys(Ro, Ri, y1, y2):
            ax.add_patch(Polygon(
                poly, closed=True, facecolor="#90caf9", edgecolor="#1565c0",
                linewidth=0.45, alpha=0.9, zorder=2,
            ))

    for y, z, A in sec["steel"]:
        ax.plot(
            z, y, "o", color="#0d47a1", ms=3.5, zorder=4,
            markeredgecolor="white", markeredgewidth=0.3,
        )

    ax.add_patch(Circle((0, 0), Ro, fill=False, ec="#333", lw=1.1, zorder=6))
    ax.add_patch(Circle((0, 0), Ri, fill=False, ec="#333", lw=0.9, ls="--", zorder=6))

    ax.set_aspect("equal")
    lim = 1.25 * Ro
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel(r"$z$ (m)")
    ax.set_ylabel(r"$y$ (m)")
    ax.set_title(rf"Pile Fiber tube — as modeled ($\times {nrow}$ area)")
    ax.legend(
        handles=[
            mpatches.Patch(facecolor="#90caf9", edgecolor="#1565c0", label="steel strips"),
            plt.Line2D(
                [0], [0], marker="o", color="w", markerfacecolor="#0d47a1",
                markersize=6, label="fiber centroid",
            ),
        ],
        loc="upper right", fontsize=8, frameon=False,
    )
    ax.grid(True, alpha=0.25, lw=0.5)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"PlotFiberSections: wrote {out}")


def main() -> int:
    json_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_JSON
    if not json_path.is_file():
        print(f"missing {json_path}", file=sys.stderr)
        return 1
    data = load(json_path)

    if len(sys.argv) > 2:
        out_dir = Path(sys.argv[2])
    else:
        # Fibers are not BC-specific; nest under active soil profile when known
        prof = data.get("soilProfile")
        if prof is None:
            sp = HERE / "soil_profile.json"
            if sp.is_file():
                prof = load(sp).get("soilProfile", 1)
            else:
                prof = 1
        out_dir = fibers_dir(prof)
    out_dir.mkdir(parents=True, exist_ok=True)

    n = 0
    if data.get("pier"):
        plot_pier(data["pier"], out_dir / "fiber_pier.png")
        n += 1
    if data.get("pile"):
        plot_pile(data["pile"], out_dir / "fiber_pile.png")
        n += 1
    if n == 0:
        print("PlotFiberSections: no Fiber sections in JSON (elastic-only?)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
