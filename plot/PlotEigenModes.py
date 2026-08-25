#!/usr/bin/env python3
"""
Goals
-----
Plot eigenmode shapes exported by DumpEigenModes.tcl.
Use one displacement scale on both components and both figure panels.

After:
  OpenSees run_gravity.tcl
  python3 plot/PlotEigenModes.py [in.json] [out.png]

Default output:
  plot/out/profile{N}/elevation/{BC}/modes/{soilEle}/{pier}/mode_XX.png

Each figure shows the pier/deck/pile zoom at left and full soil domain at right.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection

from paths import HERE, profile_root

DEFAULT_JSON = HERE / "eigen_modes.json"

# ------------------------------------------------------------
# 1. DISPLACEMENT SCALE AND FIGURE STYLE
# ------------------------------------------------------------

# Visual budget for plotted displacement (single sf on both ux and uy).
# Lateral-ish modes: max |ux| → SCALE_LATERAL · H
# Vertical-ish modes: max |uy| → SCALE_VERTICAL · H  (gentler; heave looks louder)
SCALE_LATERAL = 0.08
SCALE_VERTICAL = 0.03
STRUCT_GROUPS = ("pier", "deck", "cap", "pile", "spring", "other")
# Near-field zoom (left panel): pier / deck / piles — same idea as PlotModelSketch
ZOOM_GROUPS = frozenset({"pier", "deck", "cap", "pile", "spring"})
ZOOM_X_PAD = 3.2  # m beyond max |x| of structure (matches elevation sketch)
STRUCT_COLOR = {
    "pier": "#c45c12",
    "deck": "#1565c0",
    "cap": "#333333",
    "pile": "#8B5A2B",
    "spring": "#6a1b9a",
    "other": "#616161",
    "ssi_spring": "#90a4ae",
}


# ------------------------------------------------------------
# 2. INPUT AND MODE-SHAPE MAPS
# ------------------------------------------------------------


def load(path: Path) -> dict:
    """
    Read one DumpEigenModes JSON file.

    Args:    path
    Returns: decoded mode dictionary
    """
    with path.open() as f:
        return json.load(f)


def default_out_dir(data: dict) -> Path:
    """
    Build the default mode-figure directory from model metadata.

    Args:    data
    Returns: output directory (created if needed)
    """
    sp = data.get("soilProfile")
    sb = data.get("soilBoundary")
    pier = data.get("pierEleType") or "pier"
    sele = data.get("soilEleType") or "quad"
    if sp is not None and sb:
        d = (
            profile_root(sp)
            / "elevation"
            / str(sb)
            / "modes"
            / str(sele)
            / str(pier)
        )
    else:
        d = HERE / "out" / "modes" / str(sele) / str(pier)
    d.mkdir(parents=True, exist_ok=True)
    return d


def node_xy(data: dict) -> dict[int, tuple[float, float]]:
    """
    Convert dumped node rows to tag → (x, y).

    Args:    data
    Returns: node-coordinate mapping (m)
    """
    return {int(t): (float(x), float(y)) for t, x, y in data["nodes"]}


def phi_maps(data: dict) -> list[dict[int, tuple[float, float]]]:
    """
    Convert each mode table to tag → (ux, uy).

    Args:    data
    Returns: one displacement mapping per mode
    """
    out = []
    for mode_phi in data["phi"]:
        m = {}
        for row in mode_phi:
            tag, ux, uy = int(row[0]), float(row[1]), float(row[2])
            m[tag] = (ux, uy)
        out.append(m)
    return out


def domain_height(xy: dict[int, tuple[float, float]]) -> float:
    """
    Vertical extent used to set the visual displacement target.

    Args:    xy  node-coordinate mapping (m)
    Returns: domain height (m), at least 1
    """
    ymax = max((y for _, y in xy.values()), default=1.0)
    ymin = min((y for _, y in xy.values()), default=0.0)
    return max(ymax - ymin, 1.0)


# ------------------------------------------------------------
# 3. MODE SCALE AND DEFORMED GEOMETRY
# ------------------------------------------------------------


def phi_component_amps(
    phi: dict[int, tuple[float, float]],
    tags: set[int] | None = None,
) -> tuple[float, float, float]:
    """
    Find maximum mode-shape component and resultant amplitudes.

    Args:    phi, tags  optional subset of node tags
    Returns: (max|ux|, max|uy|, max resultant)
    """
    max_ux = max_uy = max_r = 0.0
    if tags is None:
        vals = phi.values()
    else:
        vals = (phi[t] for t in tags if t in phi)
    for ux, uy in vals:
        au, av = abs(ux), abs(uy)
        max_ux = max(max_ux, au)
        max_uy = max(max_uy, av)
        max_r = max(max_r, (ux * ux + uy * uy) ** 0.5)
    return max_ux, max_uy, max_r


def scale_for_mode(
    xy: dict[int, tuple[float, float]],
    phi: dict[int, tuple[float, float]],
) -> tuple[float, float, float, str, float]:
    """
    Choose one scale factor for ux and uy on both panels.

    Classify by max|ux| vs max|uy| over the full mesh.
    Lateral → sf so max|ux| = SCALE_LATERAL·H; vertical → max|uy| = SCALE_VERTICAL·H.

    Args:    xy, phi
    Returns: (scale_factor, controlling_amplitude, H, kind, target_fraction)
    """
    H = domain_height(xy)
    max_ux, max_uy, _ = phi_component_amps(phi, None)
    kind = "lateral" if max_ux >= max_uy else "vertical"
    target = SCALE_LATERAL if kind == "lateral" else SCALE_VERTICAL
    amp = max_ux if kind == "lateral" else max_uy
    if amp < 1.0e-30:
        return 0.0, 0.0, H, kind, target
    return target * H / amp, amp, H, kind, target


def deformed(
    xy: dict[int, tuple[float, float]],
    phi: dict[int, tuple[float, float]],
    sf: float,
) -> dict[int, tuple[float, float]]:
    """
    Apply a mode-shape scale factor to all node coordinates.

    Args:    xy, phi, sf
    Returns: deformed tag → (x, y) mapping
    """
    out = {}
    for t, (x, y) in xy.items():
        ux, uy = phi.get(t, (0.0, 0.0))
        out[t] = (x + sf * ux, y + sf * uy)
    return out


def line_segs(
    eles: list,
    xy: dict[int, tuple[float, float]],
    groups: set[str],
) -> list[np.ndarray]:
    """
    Build two-node line segments for selected element groups.

    Args:    eles, xy, groups
    Returns: coordinate arrays for plotting
    """
    segs = []
    for _e, ni, nj, grp in eles:
        if grp not in groups:
            continue
        if ni not in xy or nj not in xy:
            continue
        segs.append(np.array([xy[ni], xy[nj]]))
    return segs


def soil_polys(
    quads: list,
    xy: dict[int, tuple[float, float]],
) -> list[np.ndarray]:
    """
    Build quad polygons from node tags and a coordinate map.

    Args:    quads, xy
    Returns: valid polygon coordinate arrays
    """
    polys = []
    for q in quads:
        pts = []
        ok = True
        for n in q:
            n = int(n)
            if n not in xy:
                ok = False
                break
            pts.append(xy[n])
        if ok and len(pts) >= 3:
            polys.append(np.array(pts))
    return polys


def structure_xlim(
    xy: dict[int, tuple[float, float]],
    eles: list,
    pad: float = ZOOM_X_PAD,
) -> tuple[float, float] | None:
    """
    Find a symmetric x window around the bridge structure.

    Args:    xy, eles, pad  extra width (m)
    Returns: (xmin, xmax), or None when no structure is present
    """
    xs: list[float] = []
    for _e, ni, nj, grp in eles:
        if grp not in ZOOM_GROUPS:
            continue
        for n in (int(ni), int(nj)):
            if n in xy:
                xs.append(abs(xy[n][0]))
    if not xs:
        return None
    half = max(xs) + pad
    return (-half, half)


def domain_ylim(
    xy0: dict[int, tuple[float, float]],
    xy1: dict[int, tuple[float, float]],
) -> tuple[float, float] | None:
    """
    Find one y window containing undeformed and deformed nodes.

    Args:    xy0, xy1
    Returns: padded (ymin, ymax), or None
    """
    ys = [p[1] for p in xy0.values()] + [p[1] for p in xy1.values()]
    if not ys:
        return None
    pad = 0.05 * (max(ys) - min(ys) + 1.0)
    return (min(ys) - pad, max(ys) + pad)


# ------------------------------------------------------------
# 4. MODE PANELS
# ------------------------------------------------------------


def plot_panel(
    ax,
    xy0: dict[int, tuple[float, float]],
    xy1: dict[int, tuple[float, float]],
    eles: list,
    quads: list,
    title: str,
    bnd_quads: list | None = None,
    *,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    show_bnd: bool = True,
) -> None:
    """
    Draw undeformed and deformed soil and structure on one panel.

    Args:    ax, xy0, xy1, eles, quads, title, bnd_quads, xlim, ylim,
             show_bnd
    Returns: none (updates ax)
    """
    if bnd_quads is None:
        bnd_quads = []

    # Undeformed soil (light)
    sp0 = soil_polys(quads, xy0)
    if sp0:
        ax.add_collection(
            PolyCollection(
                sp0,
                facecolors="#cfd8dc",
                edgecolors="#90a4ae",
                linewidths=0.2,
                alpha=0.35,
            )
        )
    # Deformed soil
    sp1 = soil_polys(quads, xy1)
    if sp1:
        ax.add_collection(
            PolyCollection(
                sp1,
                facecolors="#ffe0b2",
                edgecolors="#ef6c00",
                linewidths=0.25,
                alpha=0.45,
            )
        )

    # ASDEA boundary quads (full-domain panel only)
    if show_bnd:
        b0 = soil_polys(bnd_quads, xy0)
        if b0:
            ax.add_collection(
                PolyCollection(
                    b0,
                    facecolors="#ef9a9a",
                    edgecolors="#c62828",
                    linewidths=0.35,
                    alpha=0.25,
                    zorder=1,
                )
            )
        b1 = soil_polys(bnd_quads, xy1)
        if b1:
            ax.add_collection(
                PolyCollection(
                    b1,
                    facecolors="#e53935",
                    edgecolors="#b71c1c",
                    linewidths=0.55,
                    alpha=0.40,
                    zorder=2,
                )
            )

    # Undeformed structure
    for grp in STRUCT_GROUPS:
        segs = line_segs(eles, xy0, {grp})
        if segs:
            ax.add_collection(
                LineCollection(
                    segs,
                    colors="#9e9e9e",
                    linewidths=0.8,
                    alpha=0.7,
                    zorder=3,
                )
            )
    # Deformed structure
    for grp in STRUCT_GROUPS:
        segs = line_segs(eles, xy1, {grp})
        if segs:
            ax.add_collection(
                LineCollection(
                    segs,
                    colors=STRUCT_COLOR.get(grp, "#333"),
                    linewidths=1.6,
                    zorder=4,
                )
            )

    ax.set_aspect("equal", adjustable="box")
    ax.set_title(title, fontsize=11)
    ax.tick_params(labelsize=9)
    ax.grid(False)
    ax.set_xlabel("x (m)", fontsize=10)
    ax.set_ylabel("y (m)", fontsize=10)

    if ylim is not None:
        ax.set_ylim(*ylim)
    else:
        yl = domain_ylim(xy0, xy1)
        if yl is not None:
            ax.set_ylim(*yl)

    if xlim is not None:
        ax.set_xlim(*xlim)
    else:
        xs = [p[0] for p in xy0.values()] + [p[0] for p in xy1.values()]
        if xs:
            pad_x = 0.05 * (max(xs) - min(xs) + 1.0)
            ax.set_xlim(min(xs) - pad_x, max(xs) + pad_x)


# ------------------------------------------------------------
# 5. COMMAND-LINE ENTRY POINT
# ------------------------------------------------------------


def main() -> int:
    """
    Read mode data and write one PNG per available mode.

    Args:    command-line arguments in sys.argv
    Returns: process status code
    """
    json_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_JSON
    if not json_path.is_file():
        print(f"PlotEigenModes: missing {json_path}; run OpenSees run_gravity.tcl first", file=sys.stderr)
        return 1

    data = load(json_path)
    xy0 = node_xy(data)
    phis = phi_maps(data)
    meta = data.get("modes_meta", [])
    eles = data.get("elements", [])
    quads = data.get("soil_quads", [])
    bnd_quads = data.get("bnd_quads", [])
    nModes = min(len(phis), len(meta), int(data.get("nModes", len(phis))))

    if len(sys.argv) > 2:
        out_arg = Path(sys.argv[2])
        out_dir = out_arg if out_arg.suffix == "" or out_arg.is_dir() else out_arg.parent
    else:
        out_dir = default_out_dir(data)
    out_dir.mkdir(parents=True, exist_ok=True)

    hdr = (
        f"pier={data.get('pierEleType')}  "
        f"soilEle={data.get('soilEleType', 'quad')}  "
        f"profile={data.get('soilProfile')}  BC={data.get('soilBoundary')}"
    )
    print(
        "PlotEigenModes: displacement scale  "
        f"u_plot = sf · φ  (same sf on ux, uy and on both panels);  "
        f"lateral → max|ux| = {SCALE_LATERAL}·H;  "
        f"vertical → max|uy| = {SCALE_VERTICAL}·H"
    )
    print(f"  mesh: {len(quads)} soil quads, {len(bnd_quads)} ASDEA bnd quads")
    written: list[Path] = []
    for i in range(nModes):
        m = meta[i]
        mode_id = int(m.get("mode", i + 1))
        T = m.get("T")
        tstr = f"T = {T:.4f} s" if T is not None else "T = —"
        sf, amp, H, kind, target = scale_for_mode(xy0, phis[i])
        xy1 = deformed(xy0, phis[i], sf)
        ylim = domain_ylim(xy0, xy1)
        xz = structure_xlim(xy0, eles)

        # Left: pier/deck/pile zoom; right: full soil domain (like elevation sketch)
        fig, (ax_z, ax_f) = plt.subplots(
            1,
            2,
            figsize=(13.5, 7.0),
            gridspec_kw={"width_ratios": [1.0, 1.35]},
        )
        plot_panel(
            ax_z,
            xy0,
            xy1,
            eles,
            quads,
            "pier / deck / piles",
            bnd_quads=bnd_quads,
            xlim=xz,
            ylim=ylim,
            show_bnd=False,
        )
        plot_panel(
            ax_f,
            xy0,
            xy1,
            eles,
            quads,
            "full soil domain",
            bnd_quads=bnd_quads,
            ylim=ylim,
            show_bnd=True,
        )
        fig.suptitle(
            f"Eigenmode {mode_id} — {tstr}   {kind}  "
            f"sf={sf:.4g}  (target {target}·H)   ({hdr})",
            fontsize=11,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        out_path = out_dir / f"mode_{mode_id:02d}.png"
        fig.savefig(out_path, dpi=160)
        plt.close(fig)
        written.append(out_path)
        print(
            f"  mode {mode_id:02d}: {kind:8s}  sf={sf:.6g}  amp={amp:.6g}  "
            f"H={H:.4g} m  → target |u|={target * H:.4g} m"
        )

    print(f"PlotEigenModes: wrote {len(written)} figures → {out_dir}")
    for p in written:
        print(f"  {p.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
