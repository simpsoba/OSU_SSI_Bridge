#!/usr/bin/env python3
"""
Goals
-----
Plot pile p-y, t-z, and q-z spring properties versus depth.
Mark liquefiable stations and distinguish shaft response from the pile tip.

  python3 plot/PlotPileSprings.py [in.json] [out_dir]

Default output:
  plot/out/profile{N}/pile_springs/{pult,tult,y50z50}.png
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from paths import HERE, pile_springs_dir

DEFAULT_JSON = HERE / "pile_springs.json"

# ------------------------------------------------------------
# 1. LAYER STYLE
# ------------------------------------------------------------

# Soft layer bands (sand / clay)
LAYER_SAND = "#fff3e0"
LAYER_CLAY = "#e3f2fd"


# ------------------------------------------------------------
# 2. INPUT AND DEPTH-PLOT HELPERS
# ------------------------------------------------------------


def load(path: Path) -> dict:
    """
    Read one pile-spring JSON file.

    Args:    path
    Returns: decoded spring dictionary
    """
    with path.open() as f:
        return json.load(f)


def stations_for_pile(data: dict, ip: int = 0) -> list[dict]:
    """
    Select and depth-sort stations for one pile index.

    Args:    data, ip
    Returns: station dictionaries ordered by depth
    """
    rows = [s for s in data.get("stations", []) if int(s["ip"]) == ip]
    rows.sort(key=lambda s: float(s["depth"]))
    return rows


def shade_layers(ax, layers: list[dict], y_mode: str = "depth") -> None:
    """
    Shade and label soil layers on a depth or elevation axis.

    Args:    ax, layers, y_mode  "depth" or "elev"
    Returns: none (updates ax)
    """
    for L in layers:
        yt, yb = float(L["yTop"]), float(L["yBot"])
        if y_mode == "depth":
            z0, z1 = -yt, -yb
            lo, hi = min(z0, z1), max(z0, z1)
        else:
            lo, hi = min(yt, yb), max(yt, yb)
        color = LAYER_SAND if L.get("sand") else LAYER_CLAY
        ax.axhspan(lo, hi, facecolor=color, edgecolor="none", zorder=0)
        ym = 0.5 * (lo + hi)
        ax.text(
            0.99,
            ym,
            L["name"],
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="center",
            fontsize=7,
            color="#555",
            zorder=1,
        )


def style_depth_axis(ax, depth_max: float) -> None:
    """
    Apply the common downward-positive depth axis style.

    Args:    ax, depth_max  (m)
    Returns: none (updates ax)
    """
    ax.set_ylim(depth_max * 1.02, 0.0)  # depth down
    ax.set_ylabel("Depth below grade (m)")
    ax.grid(True, which="both", ls=":", alpha=0.45)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# ------------------------------------------------------------
# 3. P-Y CAPACITY
# ------------------------------------------------------------


def plot_pult(data: dict, rows: list[dict], out: Path) -> None:
    """
    Plot p-y ultimate and residual capacities.

    Args:    data, rows, out  destination PNG
    Returns: none (writes PNG)
    """
    depth = np.array([float(s["depth"]) for s in rows])
    pult = np.array([float(s["pult"]) for s in rows]) / 1.0e3  # kN
    pRes = np.array([float(s["pRes"]) for s in rows]) / 1.0e3
    frac = float(data.get("pRes_frac", 0.15))
    liq = np.array([int(s["useLiq"]) for s in rows], dtype=bool)

    fig, ax = plt.subplots(figsize=(5.2, 7.2), constrained_layout=True)
    shade_layers(ax, data.get("layers", []))
    ax.plot(pult, depth, "o-", color="#1565c0", ms=4.5, lw=1.6, label=r"$p_\mathrm{ult}$")
    ax.plot(
        pRes,
        depth,
        "s--",
        color="#c62828",
        ms=4.0,
        lw=1.3,
        label=rf"$p_\mathrm{{res}}$ (PyLiq1: ${frac:g}\,p_\mathrm{{ult}}$; else $=p_\mathrm{{ult}}$)",
    )
    if liq.any():
        ax.plot(
            pult[liq],
            depth[liq],
            "o",
            ms=8,
            mfc="none",
            mec="#ef6c00",
            mew=1.4,
            label="PyLiq1 station",
        )
    style_depth_axis(ax, float(depth.max()))
    ax.set_xlabel(r"$p$ capacity (kN)")
    ax.set_title("Pile p-y: $p_\\mathrm{ult}$ and residual")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.92)
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"PlotPileSprings: wrote {out}")


# ------------------------------------------------------------
# 4. T-Z AND Q-Z CAPACITY
# ------------------------------------------------------------


def plot_tult(data: dict, rows: list[dict], out: Path) -> None:
    """
    Plot shaft t-z capacity and annotate the q-z pile tip.

    Args:    data, rows, out  destination PNG
    Returns: none (writes PNG)
    """
    depth = np.array([float(s["depth"]) for s in rows])
    tult = np.array([float(s["tult"]) for s in rows]) / 1.0e3  # kN
    tRes = np.array([float(s["tRes"]) for s in rows]) / 1.0e3
    frac = float(data.get("pRes_frac", 0.15))
    tip = np.array([int(s["isTip"]) for s in rows], dtype=bool)
    liq = np.array([int(s["useLiq"]) for s in rows], dtype=bool)
    shaft = ~tip

    fig, ax = plt.subplots(figsize=(5.2, 7.2), constrained_layout=True)
    shade_layers(ax, data.get("layers", []))
    ax.plot(
        tult[shaft],
        depth[shaft],
        "o-",
        color="#2e7d32",
        ms=4.5,
        lw=1.6,
        label=r"$t_\mathrm{ult}$ (shaft)",
    )
    ax.plot(
        tRes[shaft],
        depth[shaft],
        "s--",
        color="#ad1457",
        ms=4.0,
        lw=1.3,
        label=rf"$t_\mathrm{{res}}$ (TzLiq1: ${frac:g}\,t_\mathrm{{ult}}$; else $=t_\mathrm{{ult}}$)",
    )
    if tip.any():
        # Tip q_ult ≫ shaft t_ult — keep shaft-scale xlim; annotate tip.
        q_tip = float(tult[tip][0])
        z_tip = float(depth[tip][0])
        ax.axhline(z_tip, color="#6a1b9a", ls=":", lw=1.0, alpha=0.55, zorder=2)
        ax.plot(
            [0.0],
            [z_tip],
            "D",
            color="#6a1b9a",
            ms=7,
            clip_on=False,
            label=rf"$q_\mathrm{{ult}}$ tip $= {q_tip:.0f}$ kN (off scale)",
        )
    if liq.any():
        ax.plot(
            tult[liq],
            depth[liq],
            "o",
            ms=8,
            mfc="none",
            mec="#ef6c00",
            mew=1.4,
            label="TzLiq1 station",
        )
    style_depth_axis(ax, float(depth.max()))
    # Shaft-only x-limits (tip annotated; do not stretch axis to q_ult)
    x_shaft = np.concatenate([tult[shaft], tRes[shaft]]) if shaft.any() else tult
    xmax = float(np.nanmax(x_shaft)) if x_shaft.size else 1.0
    ax.set_xlim(0.0, xmax * 1.15 if xmax > 0 else 1.0)
    ax.set_xlabel(r"$t$ capacity (kN)")
    ax.set_title("Pile t-z ultimate / residual (tip $q$ annotated)")
    ax.legend(loc="lower right", fontsize=7.5, framealpha=0.92)
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"PlotPileSprings: wrote {out}")


# ------------------------------------------------------------
# 5. REFERENCE DISPLACEMENTS
# ------------------------------------------------------------


def plot_y50_z50(data: dict, rows: list[dict], out: Path) -> None:
    """
    Plot p-y, t-z, and q-z reference displacements.

    Args:    data, rows, out  destination PNG
    Returns: none (writes PNG)
    """
    depth = np.array([float(s["depth"]) for s in rows])
    y50 = np.array([float(s["y50"]) for s in rows]) * 1.0e3  # mm
    z50 = np.array([float(s["z50"]) for s in rows]) * 1.0e3
    tip = np.array([int(s["isTip"]) for s in rows], dtype=bool)

    fig, ax = plt.subplots(figsize=(5.2, 7.2), constrained_layout=True)
    shade_layers(ax, data.get("layers", []))
    ax.plot(y50, depth, "o-", color="#1565c0", ms=4.5, lw=1.6, label=r"$y_{50}$ (p-y)")
    ax.plot(
        z50[~tip],
        depth[~tip],
        "s-",
        color="#2e7d32",
        ms=4.5,
        lw=1.6,
        label=r"$z_{50}$ (t-z shaft)",
    )
    if tip.any():
        ax.plot(
            z50[tip],
            depth[tip],
            "D",
            color="#6a1b9a",
            ms=7,
            label=r"$z_{50}$ (q-z tip)",
        )
    style_depth_axis(ax, float(depth.max()))
    ax.set_xlabel(r"$y_{50}$, $z_{50}$ (mm)")
    ax.set_title(r"Pile spring $y_{50}$ and $z_{50}$")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.92)
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"PlotPileSprings: wrote {out}")


# ------------------------------------------------------------
# 6. COMMAND-LINE ENTRY POINT
# ------------------------------------------------------------


def main() -> int:
    """
    Read pile-spring data and write all three figures.

    Args:    command-line arguments in sys.argv
    Returns: process status code
    """
    json_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_JSON
    if not json_path.is_file():
        print(
            f"missing {json_path}; run OpenSees PlotModel.tcl first",
            file=sys.stderr,
        )
        return 1

    data = load(json_path)
    rows = stations_for_pile(data, ip=0)
    if not rows:
        print("no stations for pile ip=0", file=sys.stderr)
        return 1

    if len(sys.argv) > 2:
        out_dir = Path(sys.argv[2])
    else:
        out_dir = pile_springs_dir(data.get("soilProfile", 1))
    out_dir.mkdir(parents=True, exist_ok=True)

    plot_pult(data, rows, out_dir / "pult.png")
    plot_tult(data, rows, out_dir / "tult.png")
    plot_y50_z50(data, rows, out_dir / "y50z50.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
