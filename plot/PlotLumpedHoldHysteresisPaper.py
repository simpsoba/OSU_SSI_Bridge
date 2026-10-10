#!/usr/bin/env python3
"""
Goals
-----
Paper figure for the offline lumped-plasticity, pier-base-hold case: hysteresis
at four stations, Newmark (grey, solid) vs MKR-α (black, dashed), same runs
and styling as PlotLumpedHoldIntegratorPaper.py.

  (a) near-field mesh, stations marked with the panel letters
  (b) pier base zero-length section (zeroLengthSection, ZLS-I): moment vs rotation
  (c) p-y spring on the center pile at the cap soffit (y = -0.99 m): p vs y
  (d) center-pile section just below the cap (first segment, first IP):
      moment vs curvature
  (e) soil quad under the center pile, one row above the base boundary (the
      bottom row is tied to the base): tau_xy vs gamma_xy

(d) moment: dispBeamColumn sections were recorded as deformation only. The
moment at the first Gauss-Legendre point (s = (1 - sqrt(3/5))/2 of the length)
is interpolated from the element end moments (globalForce), M = -M_i (1 - s) +
M_j s. On the elastic segments below this reproduces the section curvature
exactly (correlation 1.000); after yield it is an approximation.

Pile and spring forces are for the model row (n_pile_row = 2 piles into the page).

Writes
------
  plot/out/eq_offline/compare/lumped_hold_hysteresis_every.{png,pdf}

Usage
-----
  python plot/PlotLumpedHoldHysteresisPaper.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

import PlotLumpedHoldIntegratorPaper as disp
import PlotResponseSpectrum as prs

ROOT = disp.ROOT
OUT = ROOT / "compare" / "lumped_hold_hysteresis_every"

# Stations (element tags in window_eles.txt).
SPRING_TAG = 22021  # p-y/t-z zeroLength, center pile, y = -0.99 m (dir 1 = p-y)
PILE_BEAM_TAG = 2120  # center pile, first segment below the cap (-0.99 to -1.91 m)
QUAD_TAG = 15485  # second row from the base, under the center pile (bottom row 15486 is boundary-dominated)
HINGE_NODE = 2  # pier base hinge sits between nodes 1 and 2 (top of cap)
QUADS_PER_FILE = 120  # window_quad_*_NN.out chunking
# First Gauss-Legendre point of the 3-point dispBeamColumn, as a fraction of L.
S_IP1 = 0.5 * (1.0 - np.sqrt(0.6))

ML, MR, MB, MT = 0.03, 0.14, 0.42, 0.24
MESH_H = disp.MESH_H  # same mesh size as the displacement figure
MESH_W = MESH_H * disp.MESH_DX / disp.MESH_DY
PW, PH = 1.45, 1.25  # loop panel size (in)
COL_GAP = 0.56  # mesh -> first panel column (y labels)
PCOL_GAP = 0.58  # between panel columns (y tick labels + label of the right column)
PROW_GAP = 0.66  # between panel rows (x tick labels + label of the top row + title below)
P_X0 = ML + MESH_W + COL_GAP
# Width follows the content; only the mesh size is tied to the displacement figure.
FIG_W = P_X0 + 2 * PW + PCOL_GAP + MR
FIG_H = MT + 2 * PH + PROW_GAP + MB
LEG_DROP = 0.30


def rows_of(path: Path) -> list[list[str]]:
    """Non-comment rows of a whitespace table."""
    out = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip() and not line.startswith("#"):
            out.append(line.split())
    return out


def column_index(path: Path, tag: int) -> int:
    """Position of element ``tag`` in a *_eles.txt / window_quads.txt list."""
    tags = [int(float(r[0])) for r in rows_of(path)]
    return tags.index(tag)


def load_station_series(eq: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """
    (x, y) hysteresis series per station, in plotting units.

    Returns: {"hinge": (theta mrad, M MN m), "spring": (y mm, p kN),
              "pile": (kappa 1e-3/m, M MN m), "quad": (gamma %, tau kPa)}
    """
    out: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    f = np.loadtxt(eq / "pier_hinge_force.out")
    d = np.loadtxt(eq / "pier_hinge_defo.out")
    out["time"] = (d[:, 0], d[:, 0])  # every recorder shares this time column
    out["hinge"] = (1e3 * d[:, 2], 1e-6 * f[:, 2])

    k = column_index(eq / "pile_springs_eles.txt", SPRING_TAG)
    f = np.loadtxt(eq / "pile_springs_force.out")
    d = np.loadtxt(eq / "pile_springs_defo.out")
    out["spring"] = (1e3 * d[:, 1 + 2 * k], 1e-3 * f[:, 1 + 2 * k])

    k = column_index(eq / "pile_beam_eles.txt", PILE_BEAM_TAG)
    s = np.loadtxt(eq / "pile_beam_sec1_defo.out")
    g = np.loadtxt(eq / "pile_beam_globalForce.out")
    m_i, m_j = g[:, 1 + 6 * k + 2], g[:, 1 + 6 * k + 5]
    m_ip = -m_i * (1.0 - S_IP1) + m_j * S_IP1
    out["pile"] = (1e3 * s[:, 2 + 2 * k], 1e-6 * m_ip)

    k = column_index(eq / "window_quads.txt", QUAD_TAG)
    fi, j = divmod(k, QUADS_PER_FILE)
    st = np.loadtxt(eq / f"window_quad_stress_{fi:02d}.out", usecols=[3 + 3 * j])
    sn = np.loadtxt(eq / f"window_quad_strain_{fi:02d}.out", usecols=[3 + 3 * j])
    assert len(sn) == len(out["time"][0]), "quad recorder time base differs"
    out["quad"] = (100.0 * sn, 1e-3 * st)
    return out


def station_points(xy: dict[int, tuple[float, float]], eles: dict[int, list[int]]) -> dict[str, tuple[float, float]]:
    """Mesh coordinates (m) of the four stations."""
    nd_spr = eles[SPRING_TAG][0]
    i_b, j_b = eles[PILE_BEAM_TAG][:2]
    q = eles[QUAD_TAG]
    return {
        "hinge": xy[HINGE_NODE],
        "spring": xy[nd_spr],
        "pile": (0.5 * (xy[i_b][0] + xy[j_b][0]), 0.5 * (xy[i_b][1] + xy[j_b][1])),
        "quad": (float(np.mean([xy[n][0] for n in q])), float(np.mean([xy[n][1] for n in q]))),
    }


def mark_stations(ax: plt.Axes, pts: dict[str, tuple[float, float]]) -> None:
    """Station dots on the mesh, letters off to the side with thin leaders (three sit near the cap)."""
    # (letter, station, label position in data units, ha)
    spec = [
        ("b", "hinge", (5.6, 3.0), "left"),
        ("c", "spring", (5.6, -2.6), "left"),
        ("d", "pile", (-5.0, -3.6), "right"),
        ("e", "quad", (5.6, -21.0), "left"),
    ]
    for letter, key, (lx, ly), ha in spec:
        x, y = pts[key]
        ax.scatter([x], [y], s=10, c="#111111", marker="o", zorder=5, linewidths=0.3, edgecolors="white")
        ax.annotate(
            rf"({letter})",  # station reference: gray, regular weight (panel tags are bold black)
            xy=(x, y),
            xytext=(lx, ly),
            textcoords="data",
            ha=ha,
            va="center",
            color=disp.STATION_GRAY,
            zorder=6,
            clip_on=False,
            arrowprops=dict(arrowstyle="-", lw=0.4, color=disp.STATION_GRAY, shrinkA=1.5, shrinkB=1.5),
            bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none", alpha=0.85),
        )


def loop_title(ax: plt.Axes, letter: str, name: str) -> None:
    """Panel letter and name above the axes (left aligned), clear of the loops."""
    ax.set_title(rf"\textbf{{({letter})}}~{name}", loc="left", pad=3)


def plot_loop(ax: plt.Axes, xm: np.ndarray, ym: np.ndarray, xn: np.ndarray, yn: np.ndarray) -> None:
    """Newmark (grey, under) and MKR-α (black dashed) loops."""
    ax.plot(xn, yn, color=disp.NM_COLOR, ls="-", lw=disp.NM_LW, zorder=2, solid_capstyle="round")
    ax.plot(xm, ym, color="black", ls=disp.MKR_LS, lw=disp.MKR_LW, zorder=3)


def nice_lim(peak: float) -> float:
    """Smallest {1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10} x 10^n at or above peak (axis limit = end tick)."""
    if peak <= 0.0:
        return 1.0
    exp = np.floor(np.log10(peak))
    for m in (1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0):
        if m * 10.0**exp >= peak * (1.0 - 1e-9):
            return float(m * 10.0**exp)
    return float(10.0 ** (exp + 1))


def sym_ticks(lim: float) -> list[float]:
    return [-lim, 0.0, lim]


def render() -> None:
    mkr = ROOT / disp.MKR_DUMP
    nm = ROOT / disp.NM_DUMP
    for eq in (mkr, nm):
        if not (eq / "window_meta.txt").is_file():
            raise SystemExit(f"missing dump: {eq}")

    prs.configure_font()
    xy = disp.read_xy(mkr / "window_nodes.txt")
    ele_rows = disp.read_eles(mkr / "window_eles.txt")
    eles = dict(ele_rows)
    layers = disp.read_layers(mkr / "window_quads.txt")
    s_m = load_station_series(mkr)
    s_n = load_station_series(nm)

    fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=150)
    ax_a = fig.add_axes((ML / FIG_W, (FIG_H - MT - MESH_H) / FIG_H, MESH_W / FIG_W, MESH_H / FIG_H))
    disp.draw_mesh(ax_a, xy, ele_rows, layers, [])
    pts = station_points(xy, eles)
    mark_stations(ax_a, pts)

    panels = [
        ("b", "hinge", "Pier base zero-length section", r"Rotation, $\theta$ (mrad)", r"Moment, $M$ (MN$\cdot$m)"),
        ("c", "spring", r"$p$--$y$ spring", r"Deformation, $y$ (mm)", r"Force, $p$ (kN)"),
        ("d", "pile", "Pile section", r"Curvature, $\kappa$ ($10^{-3}$/m)", r"Moment, $M$ (MN$\cdot$m)"),
        ("e", "quad", "Soil element", r"Shear strain, $\gamma_{xy}$ (\%)", r"Shear stress, $\tau_{xy}$ (kPa)"),
    ]
    log = []
    for idx, (letter, key, name, xlab, ylab) in enumerate(panels):
        r, c = divmod(idx, 2)
        x0 = P_X0 + c * (PW + PCOL_GAP)
        y0 = FIG_H - MT - (r + 1) * PH - r * PROW_GAP
        ax = fig.add_axes((x0 / FIG_W, y0 / FIG_H, PW / FIG_W, PH / FIG_H))
        xm, ym = s_m[key]
        xn, yn = s_n[key]
        plot_loop(ax, xm, ym, xn, yn)
        disp.style_history(ax)
        xl = nice_lim(max(np.abs(xm).max(), np.abs(xn).max()))
        yl = nice_lim(max(np.abs(ym).max(), np.abs(yn).max()))
        ax.set_xlim(-xl, xl)
        ax.set_ylim(-yl, yl)
        ax.set_xticks(sym_ticks(xl))
        ax.set_yticks(sym_ticks(yl))
        ax.axhline(0.0, color="0.8", lw=0.4, zorder=1)
        ax.axvline(0.0, color="0.8", lw=0.4, zorder=1)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab, labelpad=1.5)  # label tight to the tick labels
        loop_title(ax, letter, name)
        # NRMSE of the x-axis (deformation) quantity, MKR-α vs Newmark, whole record.
        e = disp.nrmse(s_m["time"][0], xm, s_n["time"][0], xn, 0.0, disp.T_MAX)
        # Top-left corner when it is clear of the loops, else the emptiest spot.
        note = "NRMSE\n" + rf"${e:.1f}\%$"
        tx = ax.text(
            0.02, 0.97, note, transform=ax.transAxes, ha="left", va="top",
            fontsize=disp.ZOOM_TICK_FS, zorder=8,
            bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none", alpha=0.85),
        )
        bb = tx.get_window_extent(fig.canvas.get_renderer())
        xy_px = np.vstack([ax.transData.transform(np.column_stack([xa, ya])) for xa, ya in ((xm, ym), (xn, yn))])
        if np.any((xy_px[:, 0] >= bb.x0) & (xy_px[:, 0] <= bb.x1) & (xy_px[:, 1] >= bb.y0) & (xy_px[:, 1] <= bb.y1)):
            tx.remove()
            disp.place_note(ax, note, [(xm, ym), (xn, yn)])
        log.append(f"({letter}) {key}: NRMSE of x {e:.3f} %")
        log.append(
            f"({letter}) {key}: |x| max MKR {np.abs(xm).max():.4g} NM {np.abs(xn).max():.4g}, "
            f"|y| max MKR {np.abs(ym).max():.4g} NM {np.abs(yn).max():.4g}"
        )

    integ = [
        Line2D([0], [0], color=disp.NM_COLOR, lw=disp.NM_LW, ls="-",
               label="Newmark\n" + r"$\gamma = 0.5$" + "\n" + r"$\beta = 0.25$"),
        Line2D([0], [0], color="black", lw=disp.MKR_LW, ls=disp.MKR_LS,
               label=r"MKR-$\alpha$" + "\n" + r"$\rho_\infty^{2} = 0.5$"),
    ]
    dt_nm, dt_mkr = disp.read_dt(nm), disp.read_dt(mkr)
    if abs(dt_nm - dt_mkr) > 1e-12:
        raise SystemExit(f"legend assumes one dt: Newmark {dt_nm}, MKR {dt_mkr}")
    leg = fig.legend(
        handles=integ,
        title=rf"$\Delta t = {1e3 * dt_mkr:.2f}$ ms",
        loc="upper left",
        bbox_to_anchor=(0.0, (FIG_H - MT - MESH_H - LEG_DROP) / FIG_H),
        bbox_transform=fig.transFigure,
        handlelength=1.0,
        handletextpad=0.3,
        labelspacing=1.1,
        frameon=False,
        borderaxespad=0.0,
        alignment="left",
    )
    leg.get_title().set_fontsize(plt.rcParams["legend.fontsize"])

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT.with_suffix(".png"), dpi=300)
    fig.savefig(OUT.with_suffix(".pdf"))
    fig.savefig(OUT.with_suffix(".svg"))
    plt.close(fig)
    print(f"wrote {OUT.with_suffix('.png')}")
    print(f"wrote {OUT.with_suffix('.pdf')}")
    print(f"figure {FIG_W:.2f} x {FIG_H:.2f} in")
    print("stations", {k: tuple(round(v, 3) for v in p) for k, p in pts.items()})
    for line in log:
        print(" ", line)


def main() -> int:
    render()
    return 0


if __name__ == "__main__":
    sys.exit(main())
