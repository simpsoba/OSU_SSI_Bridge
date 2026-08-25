#!/usr/bin/env python3
"""
Goals
-----
Merge per-rank JSON files from ExportPartitionMap.tcl.
Plot METIS element ownership for the bridge zoom and full soil domain.

  python3 plot/PlotPartition.py plot/out/profile4/partition/Shin/quad/lumpedPlasticity/np2
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Polygon

from PlotModelSketch import (
    SPRING_GAP,
    draw_rot_spiral,
    draw_trans_coil,
    ssi_coil_ends,
)

# ------------------------------------------------------------
# 1. RANK COLORS AND ELEMENT GROUPS
# ------------------------------------------------------------

# Okabe–Ito, rank 0..7
RANK_COLORS = (
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#CC79A7",
    "#F0E442",
    "#56B4E9",
    "#E69F00",
    "#332288",
)

QUAD_GRPS = {"soil", "soil_bnd"}
SKIP_LINE_GRPS = {"ssi_spring", "spring", "other"}


# ------------------------------------------------------------
# 2. READ AND MERGE RANK FILES
# ------------------------------------------------------------


def rank_color(pid: int) -> str:
    """
    Select the repeating color for one MPI rank.

    Args:    pid  rank number
    Returns: color string
    """
    return RANK_COLORS[int(pid) % len(RANK_COLORS)]


def load_rank(path: Path) -> dict:
    """
    Read one rank JSON file.

    Args:    path
    Returns: decoded rank dictionary
    """
    with path.open() as f:
        return json.load(f)


def merge_ranks(out_dir: Path) -> dict:
    """
    Merge rank.*.json files and verify unique element ownership.

    Args:    out_dir  partition output directory
    Returns: merged partition dictionary (also writes partition.json)
    """
    files = sorted(out_dir.glob("rank.*.json"), key=lambda p: int(p.suffixes[0].lstrip(".")))
    if not files:
        raise FileNotFoundError(f"no rank.*.json in {out_dir}")
    ranks = [load_rank(p) for p in files]
    np = int(ranks[0]["np"])
    if len(ranks) != np:
        raise ValueError(f"{out_dir}: expected {np} rank files, got {len(ranks)}")
    owner: dict[int, int] = {}
    elements: list[dict] = []
    for data in ranks:
        pid = int(data["pid"])
        for el in data.get("elements", []):
            e = int(el["e"])
            if e in owner:
                raise ValueError(f"element {e} owned by rank {owner[e]} and {pid}")
            owner[e] = pid
            row = dict(el)
            row["pid"] = pid
            elements.append(row)
    counts = [0] * np
    for pid in owner.values():
        counts[pid] += 1
    head = ranks[0]
    merged = {
        "np": np,
        "units": head.get("units", "m"),
        "pierEleType": head.get("pierEleType"),
        "pileEleType": head.get("pileEleType"),
        "soilEleType": head.get("soilEleType"),
        "soilProfile": head.get("soilProfile"),
        "soilBoundary": head.get("soilBoundary"),
        "sizes": head.get("sizes", {}),
        "counts": counts,
        "elements": sorted(elements, key=lambda r: int(r["e"])),
    }
    out_json = out_dir / "partition.json"
    with out_json.open("w") as f:
        json.dump(merged, f, indent=2)
        f.write("\n")
    return merged


# ------------------------------------------------------------
# 3. GEOMETRY AND SSI LOOKUPS
# ------------------------------------------------------------


def xy_pairs(xy: list) -> list[tuple[float, float]]:
    """
    Convert a flat coordinate list to (x, y) pairs.

    Args:    xy
    Returns: coordinate pairs
    """
    pts = []
    for i in range(0, len(xy) - 1, 2):
        pts.append((float(xy[i]), float(xy[i + 1])))
    return pts


def first_xy(el: dict) -> tuple[float, float] | None:
    """
    Read the first coordinate pair stored for an element.

    Args:    el
    Returns: (x, y), or None
    """
    pts = xy_pairs(el.get("xy") or [])
    if not pts:
        return None
    return pts[0]


def ssi_stypes(x: float, y: float, pile_xs: set[float], x_face: float,
               ymin_at: dict[float, float], H_cap: float) -> tuple[str, tuple[str, ...]]:
    """
    Infer station kind and uniaxial materials for one SSI element.

    Args:    x, y, pile_xs, x_face, ymin_at, H_cap  geometry (m)
    Returns: (station_kind, spring_types)
    """
    xr = round(x, 4)
    if xr in pile_xs:
        if abs(y - ymin_at.get(xr, y)) < 1.0e-6:
            return "pile", ("py", "qz")
        return "pile", ("py", "tz")
    if abs(abs(x) - x_face) < 0.08:
        return "cap", ("py", "tz")
    if abs(y + H_cap) < 0.12:
        return "cap", ("qz",)
    return "pile", ("py", "tz")


# ------------------------------------------------------------
# 4. PARTITION MAP
# ------------------------------------------------------------


def plot_map(data: dict, out: Path) -> None:
    """
    Draw the near-field and full-domain rank ownership maps.

    Args:    data  merged partition dictionary; out  destination PNG
    Returns: none (writes PNG)
    """
    plt.rcParams["mathtext.fontset"] = "cm"
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.size"] = 9

    np = int(data["np"])
    counts = data.get("counts") or [0] * np
    sz = data.get("sizes") or {}
    quads = [e for e in data["elements"] if e.get("grp") in QUAD_GRPS]
    lines = [
        e
        for e in data["elements"]
        if e.get("grp") not in QUAD_GRPS and e.get("grp") not in SKIP_LINE_GRPS
    ]
    nrot = sum(1 for e in data["elements"] if e.get("grp") == "spring")
    nssi = sum(1 for e in data["elements"] if e.get("grp") == "ssi_spring")

    fig, (ax, axf) = plt.subplots(
        1, 2, figsize=(11.0, 9.5), dpi=160,
        gridspec_kw={"width_ratios": [1.0, 1.35]},
    )

    def draw_on(a: plt.Axes, *, full: bool) -> None:
        """
        Draw rank-colored elements on one panel.

        Args:    a, full
        Returns: none (updates a)
        """
        by_pid: dict[int, list] = {i: [] for i in range(np)}
        for q in quads:
            pts = xy_pairs(q.get("xy") or [])
            if len(pts) < 3:
                continue
            by_pid[int(q["pid"])].append(Polygon(pts, closed=True))
        for pid, polys in by_pid.items():
            if not polys:
                continue
            a.add_collection(PatchCollection(
                polys,
                facecolor=rank_color(pid),
                edgecolor="#333333",
                linewidths=0.15,
                alpha=0.72 if full else 0.80,
                zorder=0,
            ))

        segs, colors = [], []
        for el in lines:
            pts = xy_pairs(el.get("xy") or [])
            if len(pts) < 2:
                continue
            x0, y0 = pts[0]
            x1, y1 = pts[1]
            if (x0 - x1) ** 2 + (y0 - y1) ** 2 < 1.0e-18:
                continue
            segs.append([(x0, y0), (x1, y1)])
            colors.append(rank_color(int(el["pid"])))
        if segs:
            a.add_collection(LineCollection(
                segs, colors=colors, linewidths=1.6, zorder=3,
            ))

        rot = [e for e in data["elements"] if e.get("grp") == "spring"]
        ssi = [e for e in data["elements"] if e.get("grp") == "ssi_spring"]
        for el in rot:
            xy0 = first_xy(el)
            if xy0 is None:
                continue
            x, y = xy0
            col = rank_color(int(el["pid"]))
            if y > 1.0:
                draw_rot_spiral(a, x, y - SPRING_GAP, x, y, col)
            else:
                draw_rot_spiral(a, x, y, x, y + SPRING_GAP, col)
        if ssi and not full:
            H_c = float(sz.get("H_cap") or 1.0)
            xs_all = []
            ymin_at: dict[float, float] = {}
            nx: Counter[float] = Counter()
            for el in ssi:
                xy0 = first_xy(el)
                if xy0 is None:
                    continue
                x, y = xy0
                xs_all.append(abs(x))
                xr = round(x, 4)
                nx[xr] += 1
                ymin_at[xr] = y if xr not in ymin_at else min(ymin_at[xr], y)
            pile_xs = {xr for xr, n in nx.items() if n >= 8}
            x_face = max(xs_all) if xs_all else 0.0
            for el in ssi:
                xy0 = first_xy(el)
                if xy0 is None:
                    continue
                x, y = xy0
                kind, stypes = ssi_stypes(x, y, pile_xs, x_face, ymin_at, H_c)
                col = rank_color(int(el["pid"]))
                for st in stypes:
                    x0, y0, x1, y1 = ssi_coil_ends(st, x, y, x, y, kind)
                    draw_trans_coil(a, x0, y0, x1, y1, col)

        a.set_xlabel(r"$x$ (m)")
        a.set_ylabel(r"$y$ (m)")
        a.grid(True, alpha=0.25, lw=0.5)

        ys = []
        xs = []
        for el in data["elements"]:
            for x, y in xy_pairs(el.get("xy") or []):
                xs.append(x)
                ys.append(y)
        H_pier = float(sz.get("H_pier") or (max(ys) if ys else 8.0))
        H_cap = float(sz.get("H_cap") or 1.0)
        L_pile = float(sz.get("L_pile") or 18.0)
        y_bot = -H_cap - L_pile - 0.8
        if ys:
            y_bot = min(y_bot, min(ys) - 0.5)
        y_top = H_pier + 2.0
        if sz.get("dd_deck"):
            y_top = max(y_top, H_pier + float(sz["dd_deck"]) + 2.0)
        if ys:
            y_top = max(y_top, max(ys) + 2.0)

        pier = data.get("pierEleType") or ""
        profile = data.get("soilProfile")
        boundary = data.get("soilBoundary") or ""
        if full:
            if "xMeshHalf" in sz:
                Lh = float(sz["xMeshHalf"])
            elif "L_half" in sz:
                Lh = float(sz["L_half"])
            elif xs:
                Lh = max(abs(min(xs)), abs(max(xs)))
            else:
                Lh = 60.0
            pad = float(sz.get("w_FF") or 0.0) + 2.0
            a.set_xlim(-Lh - pad, Lh + pad)
            ttl = rf"full domain  $n_p$={np}"
            if boundary:
                ttl += f"  ({boundary})"
            a.set_title(ttl)
            a.set_ylim(y_bot, y_top)
            a.set_aspect("equal", adjustable="box")
        else:
            W_cap = float(sz.get("W_cap") or 4.6)
            half = 0.5 * W_cap + 3.2
            if sz.get("dw_deck"):
                half = max(half, 0.5 * float(sz["dw_deck"]) + 0.8)
            a.set_xlim(-half, half)
            ttl = rf"pier={pier}"
            if profile is not None:
                ttl += rf"  soil={profile}"
            a.set_title(ttl)
            a.set_ylim(y_bot, y_top)
            a.set_aspect("equal", adjustable="box")

    draw_on(ax, full=False)
    draw_on(axf, full=True)

    handles = [
        mpatches.Patch(
            facecolor=rank_color(i), edgecolor="#333333", alpha=0.75,
            label=rf"rank {i}  ({counts[i]} eles)",
        )
        for i in range(np)
    ]
    fig.tight_layout()
    pos_L = ax.get_position()
    pos_R = axf.get_position()
    xc = pos_R.x0 + 0.5 * pos_R.width
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(xc, pos_L.y1),
        bbox_transform=fig.transFigure,
        ncol=min(4, np),
        fontsize=7.0,
        frameon=False,
        borderaxespad=0.0,
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(
        f"PlotPartition: wrote {out}  np={np}  "
        f"counts={counts}  quads={len(quads)}  beams={len(lines)}  "
        f"ZLS={nrot}  ssi={nssi}"
    )


# ------------------------------------------------------------
# 5. COMMAND-LINE ENTRY POINT
# ------------------------------------------------------------


def main() -> int:
    """
    Merge a partition directory and write its map.

    Args:    command-line arguments in sys.argv
    Returns: process status code
    """
    if len(sys.argv) < 2:
        print(
            "usage: python3 plot/PlotPartition.py <partitionOutDir>",
            file=sys.stderr,
        )
        return 1
    out_dir = Path(sys.argv[1])
    if not out_dir.is_dir():
        print(f"missing dir {out_dir}", file=sys.stderr)
        return 1
    data = merge_ranks(out_dir)
    if len(sys.argv) > 2:
        png_path = Path(sys.argv[2])
    else:
        png_path = out_dir / "partition.png"
    plot_map(data, png_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
