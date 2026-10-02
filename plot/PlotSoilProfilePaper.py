#!/usr/bin/env python3
"""
Goals
-----
Paper-style soil + SSI spring properties versus depth (1x4).

  (a) Gr & Br   (b) su   (c) p_ult & t_ult   (d) y50, z50

  python plot/PlotSoilProfilePaper.py
  python plot/PlotSoilProfilePaper.py soil_profile.json pile_springs.json [out_stem]

Default: plot/out/profile{N}/soil_profile/soil_props_paper.{pdf,png,svg}
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from matplotlib.transforms import blended_transform_factory

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import PlotResponseSpectrum as prs
from paths import soil_profile_dir

DEFAULT_SOIL = HERE / "soil_profile.json"
DEFAULT_SPRINGS = HERE / "pile_springs.json"

# Single row — match suthesis3 textwidth (letter, 1.5 in + 1 in margins).
# Height above 2.4 in so 9 pt labels / Soft–Stiff are not cramped at 6 in wide.
FIG_W = 6.0  # in
FIG_H = 2.9  # in  (room for 2-line x-labels + panel tag)
FONT_SIZE = 9  # match spectrum / elevation paper figures
LINE_COLOR = "black"
LINE_LW = 1.0
MARKER_MS = 2.0
CONTACT_C = "#9e9e9e"

UNIT_LABEL = {
    "L2": r"Soft",
    "L3": r"Medium",
    "L5": r"Stiff",
}


def load(path: Path) -> dict:
    """
    Read a JSON dump.

    Args:    path
    Returns: decoded dict
    """
    with path.open(encoding="utf-8") as f:
        return json.load(f)


def step_xy(layers: list[dict], key: str, scale: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """
    Layer-wise step profile (value constant on each soil row).

    Args:    layers, key, scale
    Returns: (values, depths_m)
    """
    zs: list[float] = []
    vs: list[float] = []
    for L in layers:
        z0 = float(L["depthTop"])
        z1 = float(L["depthBot"])
        if z1 < z0:
            z0, z1 = z1, z0
        v = float(L.get(key, 0.0)) * scale
        zs.extend([z0, z1])
        vs.extend([v, v])
    return np.asarray(vs), np.asarray(zs)


def unit_spans_from_soil(layers: list[dict]) -> list[tuple[str, float, float]]:
    """
    Contiguous named units from soil_profile.json rows.

    Args:    layers
    Returns: list of (name, z_top, z_bot)
    """
    spans: list[tuple[str, float, float]] = []
    i = 0
    n = len(layers)
    while i < n:
        nm = str(layers[i]["name"])
        lo = float(layers[i]["depthTop"])
        j = i
        while j + 1 < n and layers[j + 1]["name"] == nm:
            j += 1
        hi = float(layers[j]["depthBot"])
        if hi < lo:
            lo, hi = hi, lo
        spans.append((nm, lo, hi))
        i = j + 1
    return spans


def mark_units(
    ax,
    spans: list[tuple[str, float, float]],
    *,
    label: bool,
    soft_x: float | None = None,
) -> None:
    """
    Layer contacts + Soft / Medium / Stiff labels (no fill).

    Args:    ax, spans, label
             soft_x  if set, Soft sits at this data-x (just right of Gr)
    Returns: none
    """
    for nm, lo, hi in spans:
        ax.axhline(hi, color=CONTACT_C, lw=0.45, ls=(0, (2, 2)), zorder=1)
        if not label:
            continue
        z_mid = 0.5 * (lo + hi)
        text = UNIT_LABEL.get(nm, nm)
        kw = dict(
            ha="left",
            va="center",
            fontsize=FONT_SIZE,
            color="#666666",
            zorder=2,
        )
        # Soft: data-x just past the shallow Gr series (white space, not piles).
        # Medium / Stiff stay on the left of the curve.
        if nm == "L2" and soft_x is not None:
            ax.text(soft_x, z_mid, text, **kw)
        else:
            ax.text(
                0.03, z_mid, text,
                transform=ax.get_yaxis_transform(), **kw,
            )


def draw_pile_group_faded(ax, springs: dict) -> None:
    """
    Faded cap + 3 in-plane piles on panel (a), depth-true.

    X is axes fraction (right side); Y is depth (m). Geometry from
    pile_springs.json (H_cap, L_pile, D_pile, W_cap). Outer piles at
    ±0.4 W_cap (Mackie 12 ft c/c on a 15 ft cap).

    Args:    ax  panel (a), springs  dump dict
    Returns: none
    """
    H_cap = float(springs.get("H_cap", 0.9906))
    L_pile = float(springs.get("L_pile", 18.288))
    D_pile = float(springs.get("D_pile", 0.6096))
    W_cap = float(springs.get("W_cap", 4.572))
    # Outer pile axes / cap width = 12/15 (Parameters s_pile_cap, W_cap).
    s_pile = 0.4 * W_cap

    # Cap width → ~0.34 of panel; group sits right of Soft/Medium/Stiff.
    x_c = 0.78
    span = 0.34
    scale = span / W_cap
    trans = blended_transform_factory(ax.transAxes, ax.transData)
    face = "#7a7a7a"
    edge = "#5a5a5a"
    kw = dict(
        transform=trans,
        facecolor=face,
        edgecolor=edge,
        lw=0.45,
        zorder=2,
        clip_on=True,
    )

    x_cap = x_c - 0.5 * W_cap * scale
    ax.add_patch(Rectangle((x_cap, 0.0), W_cap * scale, H_cap, alpha=0.28, **kw))

    tip_len = L_pile
    pile_w = max(D_pile * scale, 0.012)
    for dx in (-s_pile, 0.0, s_pile):
        px = x_c + dx * scale - 0.5 * pile_w
        ax.add_patch(Rectangle((px, H_cap), pile_w, tip_len, alpha=0.18, **kw))


# Shared depth window (m) and x-axis headroom for every panel.
DEPTH_MAX = 25.0
X_PAD = 1.08
N_X_INTERVALS = 4  # ticks at 0, xmax/4, …, xmax


def nice_ceil(x: float, pad: float = X_PAD) -> float:
    """
    Round padded peak up to (1–1.5–2–2.5–3–4–5–6–8–10)·10^n.

    Args:    x  data peak, pad
    Returns: xmax
    """
    if x <= 0.0 or not math.isfinite(x):
        return 1.0
    target = float(x) * pad
    exp = math.floor(math.log10(target))
    mant = target / (10.0**exp)
    for step in (1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0):
        if mant <= step + 1e-12:
            return step * (10.0**exp)
    return 10.0 * (10.0**exp)


def style_depth(ax, *, ylabel: bool) -> None:
    """
    Shared downward-positive depth axis (0 at top, DEPTH_MAX at bottom).

    Args:    ax, ylabel
    Returns: none
    """
    ax.set_ylim(DEPTH_MAX, 0.0)
    if ylabel:
        ax.set_ylabel(r"Depth (m)")
    else:
        ax.tick_params(labelleft=False)
    ax.set_yticks(np.arange(0.0, DEPTH_MAX + 0.1, 5.0))
    ax.grid(True, which="major")


def set_xlim0(ax, peak: float, *, n_intervals: int = N_X_INTERVALS) -> float:
    """
    X from 0; same pad→nice-ceil recipe and tick count on every panel.

    Args:    ax, peak  data max, n_intervals
    Returns: xmax
    """
    xmax = nice_ceil(float(peak))
    step = xmax / n_intervals
    ax.set_xlim(0.0, xmax)
    ax.set_xticks(np.arange(0.0, xmax + 0.5 * step, step))
    return xmax


def xlabel_with_tag(*lines: str) -> str:
    """
    Multi-line x-label via usetex shortstack (last line usually the panel tag).

    Args:    lines  e.g. r"Shear modulus,", r"$G_r$ (MPa)", r"\textbf{(a)}"
    Returns: LaTeX string for ax.set_xlabel
    """
    return r"\shortstack{" + r"\\".join(lines) + "}"


def stations_for_pile(data: dict, ip: int = 0) -> list[dict]:
    """
    Depth-sorted stations for one in-plane pile.

    Args:    data, ip
    Returns: station list
    """
    rows = [s for s in data.get("stations", []) if int(s["ip"]) == ip]
    rows.sort(key=lambda s: float(s["depth"]))
    return rows


def su_at_depth(layers: list[dict], depth: float) -> float:
    """
    Undrained strength c (Pa) on the soil row containing depth.

    Args:    layers  soil_profile rows, depth  (m below grade)
    Returns: c (Pa); last row if depth is past the mesh
    """
    for L in layers:
        z0 = float(L["depthTop"])
        z1 = float(L["depthBot"])
        if z1 < z0:
            z0, z1 = z1, z0
        if z0 - 1.0e-9 <= depth <= z1 + 1.0e-9:
            return float(L["c"])
    return float(layers[-1]["c"])


def plot_soil_props_paper(
    soil: dict,
    springs: dict,
    out_stem: Path,
) -> None:
    """
    1x4 paper figure: moduli, strength, spring capacities, deformations.

    Args:    soil, springs, out_stem
    Returns: none (writes pdf, svg, png)
    """
    layers = soil.get("layers", [])
    if not layers:
        raise ValueError("soil_profile.json has no layers")
    rows = stations_for_pile(springs, ip=0)
    if not rows:
        raise ValueError("pile_springs.json has no stations for ip=0")

    prs.configure_font()
    # Keep local FONT_SIZE if PlotResponseSpectrum ever diverges.
    plt.rcParams.update({
        "font.size": FONT_SIZE,
        "axes.labelsize": FONT_SIZE,
        "xtick.labelsize": FONT_SIZE,
        "ytick.labelsize": FONT_SIZE,
    })

    spans = unit_spans_from_soil(layers)

    Gr_MPa, zG = step_xy(layers, "Gr", scale=1e-6)
    Br_MPa, zB = step_xy(layers, "Br", scale=1e-6)
    su_kPa, zS = step_xy(layers, "c", scale=1e-3)

    # Per-length intensities (kN/m). Tip JSON tult is qult; reconstruct skin.
    D_pile = float(springs.get("D_pile", 0.6096))
    n_row = float(springs.get("n_pile_row", 2))
    pi = float(np.pi)

    depth = np.array([float(s["depth"]) for s in rows])
    trib = np.maximum(np.array([float(s["trib"]) for s in rows]), 1.0e-12)
    p_prime = np.array([float(s["pult"]) for s in rows]) / trib / 1.0e3

    t_prime = np.empty(len(rows), dtype=float)
    for i, s in enumerate(rows):
        if int(s.get("isTip", 0)) == 1 or str(s.get("axType", "")) == "qz":
            cu = su_at_depth(layers, float(s["depth"]))
            t_prime[i] = cu * pi * D_pile * n_row / 1.0e3
        else:
            t_prime[i] = float(s["tult"]) / trib[i] / 1.0e3

    y50_mm = np.array([float(s["y50"]) for s in rows]) * 1.0e3
    z50_mm = np.array([float(s["z50"]) for s in rows]) * 1.0e3
    is_tip = np.array([int(s.get("isTip", 0)) == 1 for s in rows], dtype=bool)
    is_shaft = ~is_tip

    # constrained layout; wspace is relative gutter (lower → wider axes).
    # Wrapped x-labels allow a tighter gutter than single-line titles did.
    fig = plt.figure(figsize=(FIG_W, FIG_H), layout="constrained")
    fig.get_layout_engine().set(h_pad=0.04, w_pad=0.02, wspace=0.06)
    gs = GridSpec(1, 4, figure=fig)

    # ---- (a) Gr / Br ----
    # Br = 50 Gr → Br[GPa] = Gr[MPa]/20; twin xlim locked so curves coincide.
    ax_a = fig.add_subplot(gs[0, 0])
    # Soft label: a little past max Gr in L2 (gap between series and piles).
    soft_Gr = max(
        (float(L["Gr"]) * 1.0e-6 for L in layers if str(L["name"]) == "L2"),
        default=0.0,
    )
    mark_units(ax_a, spans, label=True, soft_x=soft_Gr + 8.0)
    style_depth(ax_a, ylabel=True)
    ax_a.set_xlabel(
        xlabel_with_tag(r"Shear modulus,", r"$G_r$ (MPa)", r"\textbf{(a)}")
    )
    gr_xmax = set_xlim0(ax_a, float(np.nanmax(Gr_MPa)))
    draw_pile_group_faded(ax_a, springs)
    ax_a.plot(Gr_MPa, zG, color=LINE_COLOR, lw=LINE_LW, solid_capstyle="butt", zorder=3)
    ax_a2 = ax_a.twiny()
    Br_GPa = Br_MPa / 1000.0
    ax_a2.plot(
        Br_GPa, zB, color=LINE_COLOR, lw=LINE_LW, ls=(0, (4, 2)),
        solid_capstyle="butt", zorder=3,
    )
    ax_a2.set_xlabel(xlabel_with_tag(r"Bulk modulus,", r"$B_r$ (GPa)"))
    br_xmax = gr_xmax / 20.0
    br_step = br_xmax / N_X_INTERVALS
    ax_a2.set_xlim(0.0, br_xmax)
    ax_a2.set_xticks(np.arange(0.0, br_xmax + 0.5 * br_step, br_step))
    ax_a2.set_ylim(ax_a.get_ylim())
    ax_a2.spines["top"].set_visible(True)
    ax_a2.grid(False)

    # ---- (b) su ----
    ax_b = fig.add_subplot(gs[0, 1])
    mark_units(ax_b, spans, label=False)
    ax_b.plot(su_kPa, zS, color=LINE_COLOR, lw=LINE_LW, solid_capstyle="butt", zorder=3)
    style_depth(ax_b, ylabel=False)
    ax_b.set_xlabel(
        xlabel_with_tag(r"Undrained strength,", r"$s_u$ (kPa)", r"\textbf{(b)}")
    )
    set_xlim0(ax_b, float(np.nanmax(su_kPa)))

    # ---- (c) p_ult and t_ult (kN/m) ----
    ax_c = fig.add_subplot(gs[0, 2])
    mark_units(ax_c, spans, label=False)
    ax_c.plot(
        p_prime, depth, color=LINE_COLOR, lw=LINE_LW,
        marker="o", ms=MARKER_MS, zorder=3,
    )
    ax_c.plot(
        t_prime, depth, color=LINE_COLOR, lw=LINE_LW, ls=(0, (4, 2)),
        marker="^", ms=MARKER_MS + 0.4, zorder=3,
    )
    style_depth(ax_c, ylabel=False)
    ax_c.set_xlabel(xlabel_with_tag(r"Unit capacity", r"(kN/m)", r"\textbf{(c)}"))
    set_xlim0(ax_c, float(np.nanmax(p_prime)))
    z_mid = float(np.median(depth))
    ax_c.annotate(
        r"$p_\mathrm{ult}$",
        xy=(float(np.interp(z_mid, depth, p_prime)), z_mid),
        xytext=(4, 0),
        textcoords="offset points",
        ha="left",
        va="center",
        fontsize=FONT_SIZE,
        color=LINE_COLOR,
        zorder=5,
    )
    ax_c.annotate(
        r"$t_\mathrm{ult}$",
        xy=(float(np.interp(z_mid, depth, t_prime)), z_mid),
        xytext=(4, 0),
        textcoords="offset points",
        ha="left",
        va="center",
        fontsize=FONT_SIZE,
        color=LINE_COLOR,
        zorder=5,
    )

    # ---- (d) deformations ----
    ax_d = fig.add_subplot(gs[0, 3])
    mark_units(ax_d, spans, label=False)
    ax_d.plot(
        y50_mm, depth, color=LINE_COLOR, lw=LINE_LW,
        marker="o", ms=MARKER_MS, zorder=3,
    )
    ax_d.plot(
        z50_mm[is_shaft], depth[is_shaft], color=LINE_COLOR, lw=LINE_LW,
        ls=(0, (4, 2)), marker="^", ms=MARKER_MS + 0.4, zorder=3,
    )
    if np.any(is_tip):
        ax_d.plot(
            z50_mm[is_tip], depth[is_tip], color=LINE_COLOR, lw=0,
            marker="s", ms=MARKER_MS + 1.2, zorder=4,
        )
        ax_d.annotate(
            r"$z_{50}$ ($q$-$z$)",
            xy=(float(z50_mm[is_tip][0]), float(depth[is_tip][0])),
            xytext=(6, 0),
            textcoords="offset points",
            ha="left",
            va="center",
            fontsize=FONT_SIZE,
            color=LINE_COLOR,
            zorder=5,
        )
    style_depth(ax_d, ylabel=False)
    ax_d.set_xlabel(xlabel_with_tag(r"Deformation", r"(mm)", r"\textbf{(d)}"))
    set_xlim0(ax_d, max(float(np.nanmax(y50_mm)), float(np.nanmax(z50_mm))))
    if np.any(is_shaft):
        y50_x = float(np.median(y50_mm[is_shaft]))
        z50_x = float(np.median(z50_mm[is_shaft]))
        z_mid_s = float(np.median(depth[is_shaft]))
        ax_d.annotate(
            r"$y_{50}$",
            xy=(y50_x, z_mid_s),
            xytext=(-5, -12),
            textcoords="offset points",
            ha="right",
            va="center",
            fontsize=FONT_SIZE,
            color=LINE_COLOR,
            zorder=5,
        )
        ax_d.annotate(
            r"$z_{50}$ ($t$-$z$)",
            xy=(z50_x, z_mid_s),
            xytext=(5, 12),
            textcoords="offset points",
            ha="left",
            va="center",
            fontsize=FONT_SIZE,
            color=LINE_COLOR,
            zorder=5,
        )

    out_stem.parent.mkdir(parents=True, exist_ok=True)
    pdf_path = out_stem.with_suffix(".pdf")
    png_path = out_stem.with_suffix(".png")
    svg_path = out_stem.with_suffix(".svg")
    # Exact figsize (= textwidth); do not use bbox_inches="tight".
    fig.savefig(pdf_path)
    fig.savefig(svg_path)
    plt.close(fig)
    prs._png_from_pdf(pdf_path, png_path, dpi=300)
    print(f"PlotSoilProfilePaper: wrote {pdf_path}")
    print(f"PlotSoilProfilePaper: wrote {svg_path}")
    print(f"PlotSoilProfilePaper: wrote {png_path}")


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    soil_path = Path(args[0]) if len(args) >= 1 else DEFAULT_SOIL
    springs_path = Path(args[1]) if len(args) >= 2 else DEFAULT_SPRINGS

    if not soil_path.is_file():
        print(f"missing {soil_path}", file=sys.stderr)
        return 1
    if not springs_path.is_file():
        print(f"missing {springs_path}", file=sys.stderr)
        return 1

    soil = load(soil_path)
    springs = load(springs_path)

    if len(args) >= 3:
        out_stem = Path(args[2])
    else:
        out_stem = soil_profile_dir(soil.get("soilProfile", 4)) / "soil_props_paper"

    plot_soil_props_paper(soil, springs, out_stem)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
