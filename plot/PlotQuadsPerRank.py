#!/usr/bin/env python3
"""
Goals
-----
Estimate near-field continuum quads per MPI rank for each soilMesh option.
Mark the even rank count at the lower edge of the 100--150 quads/rank band.

  python3 plot/PlotQuadsPerRank.py
  python3 plot/PlotQuadsPerRank.py [out.png]

Writes plot/out/quads_per_rank.png (or the path given).

Model assumptions
-----------------
These match soil/BuildSoilMesh.tcl, SoilDxBands.tcl, and the Shin default:

  * nQuad = (nX - 1) * (nY - 1). Continuum only; springs / beams / Lysmer
    are ~constant and small (~160 eles) so they do not change the shape.
  * nY = 28 (cap stations + 3 ft down to 15 ft below the pile tips).
    Independent of soilMesh.
  * nX from the same stepped bands + pile axes at 0, +/-6 ft + Shin FF
    column (L_half + 40 ft) on each side.
  * Perfect METIS balance: quads/rank = nQuad / np. Halo is ~nY nodes per
    vertical cut, so communication is ignored here.
  * RunParallel.tcl needs np >= 2. Serial (Run.tcl) is np = 1 on the plot
    only as a reference.
  * Sweet band 100-150 quads/rank. The ring on each curve is the even np
    whose quads/rank hits the lower edge (largest even np with
    nQuad/np >= 100, capped at 32). Elastic soil wants fewer ranks;
    Quad vs SSPQuad can take one step up.

Keep BANDS in sync with soil/SoilDxBands.tcl.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT_DEFAULT = HERE / "out" / "quads_per_rank.png"

# ------------------------------------------------------------
# 1. MESH AND PLOT KNOBS
# ------------------------------------------------------------

FOOT = 0.3048
S_PILE = 6.0 * FOOT
W_FF = 40.0 * FOOT
N_Y = 28

# soilMesh -> (label, [(dx_ft, x_end_ft), ...])
BANDS = {
    -2: ("-2 coarser", [(3, 6), (12, 30), (30, 90), (55, 200)]),
    -1: ("-1 coarse", [(3, 12), (14, 40), (20, 100), (50, 200)]),
    0: ("0 production", [(3, 12), (7, 40), (15, 100), (20, 140), (30, 200)]),
    1: ("1 moderate", [(3, 39), (7, 95), (15, 140), (20, 200)]),
    2: ("2 large", [(3, 84), (7, 140), (15, 200)]),
    3: ("3 x-large", [(3, 114), (7, 170), (15, 200)]),
    4: ("4 xx-large", [(3, 123), (7, 200)]),
}

NP_MIN = 2
NP_MAX = 32
SWEET_LO = 100.0
SWEET_HI = 150.0

# Okabe-Ito
COLORS = {
    -2: "#332288",
    -1: "#56B4E9",
    0: "#0072B2",
    1: "#009E73",
    2: "#E69F00",
    3: "#D55E00",
    4: "#CC79A7",
}


# ------------------------------------------------------------
# 2. SHIN MESH COUNTS
# ------------------------------------------------------------


def _push(xs: list[float], x: float, tol: float = 1e-6) -> None:
    """
    Append one coordinate unless an equal value is already present.

    Args:    xs, x, tol
    Returns: none (updates xs)
    """
    for v in xs:
        if abs(v - x) < tol:
            return
    xs.append(x)


def _fill_band(xs: list[float], x0: float, x1: float, dx: float) -> None:
    """
    Add one stepped mesh band, including both endpoints.

    Args:    xs, x0, x1, dx  coordinates and spacing (m)
    Returns: none (updates xs)
    """
    _push(xs, x0)
    if dx <= 0.0 or x1 <= x0 + 1e-9:
        _push(xs, x1)
        return
    x = x0
    while x + dx < x1 - 1e-7:
        x = x + dx
        _push(xs, x)
    _push(xs, x1)


def n_x_shin(bands_ft: list[tuple[float, float]]) -> int:
    """
    Count Shin x stations from near-field bands, pile axes, and FF columns.

    Args:    bands_ft  list of (dx, x_end) on one half-domain (ft)
    Returns: total mirrored x-station count
    """
    l_half = bands_ft[-1][1] * FOOT
    xs: list[float] = []
    for xp in (-S_PILE, 0.0, S_PILE):
        _push(xs, xp)
    x_prev = 0.0
    for dx_ft, x_end_ft in bands_ft:
        dx = dx_ft * FOOT
        x_end = min(x_end_ft * FOOT, l_half)
        if x_end <= x_prev + 1e-9:
            continue
        _fill_band(xs, x_prev, x_end, dx)
        x_prev = x_end
    _push(xs, l_half)
    _push(xs, l_half + W_FF)
    for x in list(xs):
        if x > 1e-9:
            _push(xs, -x)
    return len(xs)


def n_quad(n_x: int) -> int:
    """
    Continuum quad count for the fixed vertical mesh.

    Args:    n_x  number of x stations
    Returns: number of quadrilateral elements
    """
    return (n_x - 1) * (N_Y - 1)


def np_at_band_floor(n_quads: int) -> int:
    """
    Largest even rank count that keeps at least SWEET_LO quads/rank.

    Args:    n_quads  continuum element count
    Returns: even np in [NP_MIN, NP_MAX]
    """
    chosen = NP_MIN
    for n in range(NP_MIN, NP_MAX + 1, 2):
        if n_quads / float(n) >= SWEET_LO:
            chosen = n
    return chosen


def mesh_counts() -> list[tuple[int, str, int, int]]:
    """
    Calculate x-station and quad counts for every soilMesh option.

    Args:    none
    Returns: rows of (mesh_id, label, n_x, n_quads)
    """
    rows = []
    for mid in sorted(BANDS):
        label, bands = BANDS[mid]
        nx = n_x_shin(bands)
        nq = n_quad(nx)
        rows.append((mid, label, nx, nq))
    return rows


# ------------------------------------------------------------
# 3. PLOT AND REPORT
# ------------------------------------------------------------


def plot_quads_per_rank(out_path: Path) -> None:
    """
    Plot quads/rank curves and print the count summary.

    Args:    out_path  destination PNG
    Returns: none (writes PNG and prints counts)
    """
    np_vals = list(range(NP_MIN, NP_MAX + 1))
    rows = mesh_counts()

    fig, ax = plt.subplots(figsize=(9.4, 5.2))
    ax.axhspan(SWEET_LO, SWEET_HI, color="#eeeeee", zorder=0)
    ax.axhline(SWEET_LO, color="#bbbbbb", lw=0.8, zorder=1)
    ax.axhline(SWEET_HI, color="#bbbbbb", lw=0.8, zorder=1)
    ax.axvline(8, color="#888888", ls=":", lw=0.9, zorder=1)

    for mid, label, nx, nq in rows:
        y = [nq / float(n) for n in np_vals]
        ax.plot(
            np_vals,
            y,
            color=COLORS[mid],
            lw=1.6,
            marker="o",
            ms=3.2,
            label=f"{label}  ({nq} quads)",
            zorder=2,
        )
        rec = np_at_band_floor(nq)
        ax.plot(
            rec,
            nq / float(rec),
            marker="o",
            ms=9,
            color=COLORS[mid],
            markeredgecolor="k",
            markeredgewidth=0.7,
            zorder=3,
        )

    ax.set_xlim(2, 32)
    ax.set_xticks(list(range(2, 33, 2)))
    ax.set_ylim(0, 400)
    ax.set_xlabel("np  (MPI ranks)")
    ax.set_ylabel("quads / rank  (nQuad / np)")
    ax.set_title("Continuum quads per rank vs np")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.92)
    ax.text(
        0.02,
        0.04,
        "grey band: 100-150 quads/rank   dotted: np=8 (current matrix)   "
        "ring: even np at lower band edge   nY=28 Shin",
        transform=ax.transAxes,
        fontsize=7.5,
        color="#444",
    )
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

    print("nY=%d  nQuad=(nX-1)*(nY-1)" % N_Y)
    print("%6s %12s %4s %6s %s" % ("mesh", "name", "nX", "nQuad", "sweet np (quads/rank)"))
    for mid, label, nx, nq in rows:
        rec = np_at_band_floor(nq)
        print(
            "%6d %12s %4d %6d   np=%2d -> %.0f"
            % (mid, BANDS[mid][0], nx, nq, rec, nq / float(rec))
        )
    print("wrote", out_path)


# ------------------------------------------------------------
# 4. COMMAND-LINE ENTRY POINT
# ------------------------------------------------------------


if __name__ == "__main__":
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else OUT_DEFAULT
    plot_quads_per_rank(out)
