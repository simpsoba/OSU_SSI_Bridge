"""PT39 box solid-concrete outline polygons (local y = 0 at soffit bottom)."""

from __future__ import annotations

from matplotlib.patches import Polygon


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
    """Solid concrete outline patches; y0 shifts soffit bottom to global elevation."""
    x_top_outer = 0.5 * dw
    x_overhang_in = 0.5 * dw - cw
    x_soffit = 0.5 * sw

    def _poly(pts: list[tuple[float, float]]) -> Polygon:
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
