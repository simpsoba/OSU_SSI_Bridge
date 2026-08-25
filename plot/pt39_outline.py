"""
Goals
-----
Build reusable matplotlib polygons for the PT39 solid-concrete outline.
Use local y = 0 at the soffit bottom, with an optional global elevation shift.
"""

from __future__ import annotations

from matplotlib.patches import Polygon


# ------------------------------------------------------------
# 1. PT39 OUTLINE
# ------------------------------------------------------------


def pt39_outline(
    dw: float,
    dd: float,
    sw: float,
    cw: float,
    td: float,
    ts: float,
    tw: float,
    *,
    y0: float = 0.0,
) -> list[Polygon]:
    """
    Build top, soffit, web, and cantilever concrete patches.

    Args:    dw, dd, sw, cw, td, ts, tw  section dimensions (m);
             y0  soffit-bottom elevation (m)
    Returns: six matplotlib Polygon patches
    """
    x_top_outer = 0.5 * dw
    x_overhang_in = 0.5 * dw - cw
    x_soffit = 0.5 * sw

    def _poly(pts: list[tuple[float, float]]) -> Polygon:
        """
        Shift local section points to the requested elevation.

        Args:    pts  local (x, y) coordinates (m)
        Returns: closed Polygon
        """
        return Polygon([(x, y0 + y) for x, y in pts], closed=True)

    top = _poly([
        (-x_top_outer, dd - td),
        (x_top_outer, dd - td),
        (x_top_outer, dd),
        (-x_top_outer, dd),
    ])
    sof = _poly([
        (-x_soffit, 0.0),
        (x_soffit, 0.0),
        (x_soffit, ts),
        (-x_soffit, ts),
    ])
    x_top_web = x_overhang_in
    left = _poly([
        (-x_soffit, ts),
        (-x_soffit + tw, ts),
        (-x_top_web + tw, dd - td),
        (-x_top_web, dd - td),
    ])
    right = _poly([
        (x_soffit - tw, ts),
        (x_soffit, ts),
        (x_top_web, dd - td),
        (x_top_web - tw, dd - td),
    ])
    cant_L = _poly([
        (-x_top_outer, dd - 0.7 * td),
        (-x_overhang_in, dd - td),
        (-x_overhang_in, dd),
        (-x_top_outer, dd),
    ])
    cant_R = _poly([
        (x_overhang_in, dd - td),
        (x_top_outer, dd - 0.7 * td),
        (x_top_outer, dd),
        (x_overhang_in, dd),
    ])
    return [top, sof, left, right, cant_L, cant_R]
