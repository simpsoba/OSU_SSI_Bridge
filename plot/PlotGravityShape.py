#!/usr/bin/env python3
"""
Goals
-----
Plot the post-gravity deformed shape exported by DumpGravityShape.tcl.
Show the bridge zoom and full soil domain with one fixed displacement scale.

The undeformed mesh is overlaid with displacement amplified by SCALE_FACTOR.
The same factor applies to ux and uy on both panels.

  python3 plot/PlotGravityShape.py [in.json] [out.png]
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt

from paths import HERE, gravity_dir
from PlotEigenModes import (
    deformed,
    domain_ylim,
    plot_panel,
    structure_xlim,
)

DEFAULT_JSON = HERE / "gravity_shape.json"

# ------------------------------------------------------------
# 1. DISPLAY SCALE
# ------------------------------------------------------------

# Fixed visual amplification of gravity displacements (same sf on ux, uy).
SCALE_FACTOR = 10.0


# ------------------------------------------------------------
# 2. INPUT AND OUTPUT PATHS
# ------------------------------------------------------------


def load(path: Path) -> dict:
    """
    Read one gravity-shape JSON file.

    Args:    path
    Returns: decoded gravity dictionary
    """
    with path.open() as f:
        return json.load(f)


def default_out_path(data: dict) -> Path:
    """
    Build the default gravity-figure path from model metadata.

    Args:    data
    Returns: destination PNG path
    """
    sp = data.get("soilProfile")
    sb = data.get("soilBoundary")
    sele = data.get("soilEleType") or "quad"
    pier = data.get("pierEleType") or "lumpedPlasticity"
    if sp is not None and sb:
        d = gravity_dir(sp, sb, sele, pier)
    else:
        d = HERE / "out" / "gravity" / str(sele) / str(pier)
    d.mkdir(parents=True, exist_ok=True)
    return d / "gravity_deformed.png"


# ------------------------------------------------------------
# 3. DISPLACEMENT MAP AND SCALE
# ------------------------------------------------------------


def node_xy_u(data: dict) -> tuple[
    dict[int, tuple[float, float]],
    dict[int, tuple[float, float]],
]:
    """
    Split dumped node rows into coordinate and displacement mappings.

    Args:    data
    Returns: (tag → (x, y), tag → (ux, uy)), m
    """
    xy: dict[int, tuple[float, float]] = {}
    u: dict[int, tuple[float, float]] = {}
    for row in data["nodes"]:
        tag, x, y, ux, uy = (
            int(row[0]),
            float(row[1]),
            float(row[2]),
            float(row[3]),
            float(row[4]),
        )
        xy[tag] = (x, y)
        u[tag] = (ux, uy)
    return xy, u


def scale_disp(
    xy: dict[int, tuple[float, float]],
    u: dict[int, tuple[float, float]],
    sf: float = SCALE_FACTOR,
) -> tuple[float, float, float]:
    """
    Report the fixed scale, true maximum displacement, and domain height.

    Args:    xy, u, sf
    Returns: (scale_factor, max|u|, domain_height)
    """
    ymax = max((y for _, y in xy.values()), default=1.0)
    ymin = min((y for _, y in xy.values()), default=0.0)
    H = max(ymax - ymin, 1.0)
    amp = 0.0
    for ux, uy in u.values():
        amp = max(amp, (ux * ux + uy * uy) ** 0.5)
    return sf, amp, H


# ------------------------------------------------------------
# 4. COMMAND-LINE ENTRY POINT
# ------------------------------------------------------------


def main() -> int:
    """
    Read gravity data and write the two-panel deformed-shape figure.

    Args:    command-line arguments in sys.argv
    Returns: process status code
    """
    json_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_JSON
    if not json_path.is_file():
        print(
            f"PlotGravityShape: missing {json_path}; "
            "run gravity + DumpGravityShape.tcl first",
            file=sys.stderr,
        )
        return 1

    if len(sys.argv) > 2:
        out_path = Path(sys.argv[2])
        out_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        out_path = default_out_path(load(json_path))

    data = load(json_path)
    xy0, u = node_xy_u(data)
    eles = data.get("elements", [])
    quads = data.get("soil_quads", [])
    bnd_quads = data.get("bnd_quads", [])
    sf, amp, H = scale_disp(xy0, u)
    xy1 = deformed(xy0, u, sf)

    # True (unscaled) settlement stats
    uys = [uy for _, uy in u.values()]
    uxs = [ux for ux, _ in u.values()]
    uy_min = min(uys) if uys else 0.0
    uy_max = max(uys) if uys else 0.0
    ux_abs = max((abs(v) for v in uxs), default=0.0)

    hdr = (
        f"pier={data.get('pierEleType')}  "
        f"soilEle={data.get('soilEleType', 'quad')}  "
        f"profile={data.get('soilProfile')}  BC={data.get('soilBoundary')}"
    )
    ylim = domain_ylim(xy0, xy1)
    xz = structure_xlim(xy0, eles)

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
        f"Post-gravity deformed (×{SCALE_FACTOR:g})  sf={sf:.4g}   ({hdr})\n"
        f"true: uy∈[{uy_min:.3e}, {uy_max:.3e}] m  max|ux|={ux_abs:.3e} m  "
        f"max|u|={amp:.3e} m  H={H:.4g} m",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(out_path, dpi=160)
    plt.close(fig)
    print(f"PlotGravityShape: wrote {out_path}")
    print(
        f"  sf={sf:.6g}  true max|u|={amp:.6g} m  "
        f"uy∈[{uy_min:.4e}, {uy_max:.4e}] m"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
