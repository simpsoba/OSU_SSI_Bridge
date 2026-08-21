#!/usr/bin/env python3
"""Near-field quads per MPI rank vs np, one curve per soilMesh.

  python3 plot/PlotQuadsPerRank.py
  python3 plot/PlotQuadsPerRank.py [out.png]

Writes plot/out/quads_per_rank.png (or the path given).

Assumptions (match soil/BuildSoilMesh.tcl + SoilDxBands.tcl, Shin default):
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
  * Sweet band 110-190 quads/rank is the PDMY / SSPQuad rule of thumb from
    the np discussion (expensive constitutive, cheap 2D halo). Elastic soil
    wants fewer ranks; Quad vs SSPQuad can take one step up.

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

# np that lands in ~110-190 quads/rank (see module docstring)
SWEET_NP = {-2: 4, -1: 4, 0: 8, 1: 12, 2: 16, 3: 16, 4: 16}
SWEET_LO = 110.0
SWEET_HI = 190.0

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


def _push(xs: list[float], x: float, tol: float = 1e-6) -> None:
    for v in xs:
        if abs(v - x) < tol:
            return
    xs.append(x)


def _fill_band(xs: list[float], x0: float, x1: float, dx: float) -> None:
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
    """Shin nX: NF bands 0 -> L_half, FF at L_half+w_FF, mirror, pile axes."""
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
    return (n_x - 1) * (N_Y - 1)


def mesh_counts() -> list[tuple[int, str, int, int]]:
    rows = []
    for mid in sorted(BANDS):
        label, bands = BANDS[mid]
        nx = n_x_shin(bands)
        nq = n_quad(nx)
        rows.append((mid, label, nx, nq))
    return rows


def plot_quads_per_rank(out_path: Path) -> None:
    np_vals = list(range(2, 25))
    rows = mesh_counts()

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
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
        rec = SWEET_NP[mid]
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

    ax.set_xlim(2, 24)
    ax.set_xticks(list(range(2, 25, 2)))
    ax.set_ylim(0, 400)
    ax.set_xlabel("np  (MPI ranks)")
    ax.set_ylabel("quads / rank  (nQuad / np)")
    ax.set_title("Continuum quads per rank vs np")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.92)
    ax.text(
        0.02,
        0.04,
        "grey band: 110-190 quads/rank   dotted: np=8 (current matrix)   "
        "ring: recommended np   nY=28 Shin",
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
        rec = SWEET_NP[mid]
        print(
            "%6d %12s %4d %6d   np=%2d -> %.0f"
            % (mid, BANDS[mid][0], nx, nq, rec, nq / float(rec))
        )
    print("wrote", out_path)


if __name__ == "__main__":
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else OUT_DEFAULT
    plot_quads_per_rank(out)
