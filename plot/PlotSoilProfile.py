#!/usr/bin/env python3
"""
Goals
-----
Plot soil continuum parameters versus depth from soil_profile.json.
Separate the general material overview from the full PDMY02 argument card.

  python3 plot/PlotSoilProfile.py [in.json] [out_dir]

Default output:
  plot/out/profile{N}/soil_profile/{overview,pdmy02}.png
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from paths import HERE, soil_profile_dir

DEFAULT_JSON = HERE / "soil_profile.json"

# ------------------------------------------------------------
# 1. LAYER STYLE
# ------------------------------------------------------------

LAYER_SAND = "#fff3e0"
LAYER_CLAY = "#e3f2fd"


# ------------------------------------------------------------
# 2. INPUT AND DEPTH-PLOT HELPERS
# ------------------------------------------------------------


def load(path: Path) -> dict:
    """
    Read one soil-profile JSON file.

    Args:    path
    Returns: decoded profile dictionary
    """
    with path.open() as f:
        return json.load(f)


def shade_layers(ax, layers: list[dict]) -> None:
    """
    Shade contiguous soil layers and label them by name.

    Args:    ax, layers
    Returns: none (updates ax)
    """
    i = 0
    n = len(layers)
    while i < n:
        nm = layers[i]["name"]
        lo = float(layers[i]["depthTop"])
        j = i
        while j + 1 < n and layers[j + 1]["name"] == nm:
            j += 1
        hi = float(layers[j]["depthBot"])
        if hi < lo:
            lo, hi = hi, lo
        color = LAYER_SAND if layers[i].get("sand") else LAYER_CLAY
        ax.axhspan(lo, hi, facecolor=color, edgecolor="none", zorder=0)
        ax.text(
            0.98,
            0.5 * (lo + hi),
            nm,
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="center",
            fontsize=6.5,
            color="#555",
            zorder=1,
        )
        i = j + 1


def style_depth(ax, zmax: float) -> None:
    """
    Apply the common downward-positive depth axis style.

    Args:    ax, zmax  maximum depth (m)
    Returns: none (updates ax)
    """
    ax.set_ylim(zmax * 1.02, 0.0)
    ax.set_ylabel("Depth below grade (m)")
    ax.grid(True, which="both", ls=":", alpha=0.45)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def step_xy(layers: list[dict], key: str, scale: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """
    Expand one layer value into a vertical step profile.

    Args:    layers, key, scale
    Returns: (values, depths) arrays
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


def plot_step(ax, layers, key, *, scale=1.0, color="#1565c0", label=None, ls="-"):
    """
    Plot one layer-wise property as a depth step.

    Args:    ax, layers, key, scale, color, label, ls
    Returns: none (updates ax)
    """
    x, z = step_xy(layers, key, scale)
    ax.plot(x, z, ls, color=color, lw=1.8, label=label)


def twin_legend(ax, ax2, loc="lower right"):
    """
    Combine legends from a primary and twinned x axis.

    Args:    ax, ax2, loc
    Returns: none (updates ax)
    """
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc=loc, framealpha=0.9)


# ------------------------------------------------------------
# 3. SOIL-PARAMETER OVERVIEW
# ------------------------------------------------------------


def plot_overview(data: dict, layers: list[dict], out: Path) -> None:
    """
    Write the eight-panel soil-parameter overview.

    Args:    data, layers, out  destination PNG
    Returns: none (writes PNG)
    """
    zmax = max(float(L["depthBot"]) for L in layers)
    prof = int(data.get("soilProfile", 0))

    fig, axes = plt.subplots(2, 4, figsize=(14.0, 9.0), constrained_layout=True)
    fig.suptitle(f"Soil continuum parameters — profile {prof}", fontsize=13)

    ax = axes[0, 0]
    shade_layers(ax, layers)
    plot_step(ax, layers, "rho", color="#37474f", label=r"$\rho$")
    style_depth(ax, zmax)
    ax.set_xlabel(r"$\rho$ (kg/m³)")
    ax.set_title("Density")

    ax = axes[0, 1]
    shade_layers(ax, layers)
    plot_step(ax, layers, "Gr", scale=1e-6, color="#1565c0", label=r"$G_r$")
    plot_step(ax, layers, "Br", scale=1e-6, color="#c62828", label=r"$B_r$", ls="--")
    style_depth(ax, zmax)
    ax.set_xlabel("Modulus (MPa)")
    ax.set_title(r"Reference $G_r$, $B_r$")
    ax.legend(fontsize=8, loc="lower right", framealpha=0.9)

    ax = axes[0, 2]
    shade_layers(ax, layers)
    plot_step(ax, layers, "phi", color="#6a1b9a", label=r"$\phi$ (°)")
    ax2 = ax.twiny()
    x, z = step_xy(layers, "c", scale=1e-3)
    ax2.plot(x, z, ":", color="#2e7d32", lw=1.8, label=r"$c$ (kPa)")
    style_depth(ax, zmax)
    ax.set_xlabel(r"$\phi$ (°)")
    ax2.set_xlabel(r"$c$ (kPa)")
    ax.set_title(r"Strength $\phi$, $c$")
    twin_legend(ax, ax2)

    ax = axes[0, 3]
    shade_layers(ax, layers)
    plot_step(ax, layers, "Dr", color="#ef6c00", label=r"$D_r$ (%)")
    ax2 = ax.twiny()
    x, z = step_xy(layers, "k_pci")
    ax2.plot(x, z, "--", color="#00838f", lw=1.8, label=r"$k$ (pci)")
    style_depth(ax, zmax)
    ax.set_xlabel(r"$D_r$ (%)")
    ax2.set_xlabel(r"API $k$ (pci)")
    ax.set_title(r"$D_r$ and API $k$")
    twin_legend(ax, ax2)

    ax = axes[1, 0]
    shade_layers(ax, layers)
    plot_step(ax, layers, "contr1", color="#1565c0", label="contrac1")
    plot_step(ax, layers, "contr3", color="#c62828", label="contrac3", ls="--")
    style_depth(ax, zmax)
    ax.set_xlabel("(-)")
    ax.set_title("Contraction (PDMY02)")
    ax.legend(fontsize=8, loc="lower right", framealpha=0.9)

    ax = axes[1, 1]
    shade_layers(ax, layers)
    plot_step(ax, layers, "dilat1", color="#2e7d32", label="dilat1")
    plot_step(ax, layers, "dilat3", color="#ad1457", label="dilat3", ls="--")
    style_depth(ax, zmax)
    ax.set_xlabel("(-)")
    ax.set_title("Dilation (PDMY02)")
    ax.legend(fontsize=8, loc="lower right", framealpha=0.9)

    ax = axes[1, 2]
    shade_layers(ax, layers)
    plot_step(ax, layers, "PTA", color="#4527a0", label=r"PTAng (°)")
    ax2 = ax.twiny()
    x, z = step_xy(layers, "d")
    ax2.plot(x, z, "--", color="#546e7a", lw=1.8, label=r"$d$")
    style_depth(ax, zmax)
    ax.set_xlabel(r"PTAng (°)")
    ax2.set_xlabel(r"$d$ (-)")
    ax.set_title(r"Phase-transform angle, $d$")
    twin_legend(ax, ax2)

    ax = axes[1, 3]
    shade_layers(ax, layers)
    plot_step(ax, layers, "B_fsp", scale=1e-6, color="#00838f", label=r"$B_\mathrm{FSP}$")
    zs, ratios = [], []
    for L in layers:
        z0, z1 = float(L["depthTop"]), float(L["depthBot"])
        if z1 < z0:
            z0, z1 = z1, z0
        Gr = float(L["Gr"])
        Br = float(L["Br"])
        r = Br / Gr if Gr > 0 else 0.0
        zs.extend([z0, z1])
        ratios.extend([r, r])
    ax2 = ax.twiny()
    ax2.plot(ratios, zs, "--", color="#6d4c41", lw=1.8, label=r"$B_r/G_r$")
    style_depth(ax, zmax)
    # Clay is 50; wiki/PEER sand is ~1.8–2.3. Autoscale around 50+1e-6 looks like a staircase.
    ax2.set_xlim(0.0, 60.0)
    ax2.ticklabel_format(useOffset=False)
    bfx, _ = step_xy(layers, "B_fsp", 1e-6)
    bfmax = float(np.max(np.abs(bfx))) if len(bfx) else 0.0
    ax.set_xlim(0.0, max(bfmax * 1.15, 1.0))
    ax.set_xlabel(r"$B_\mathrm{FSP}$ (MPa)")
    ax2.set_xlabel(r"$B_r/G_r$")
    ax.set_title(r"FSP bulk and $B_r/G_r$")
    twin_legend(ax, ax2)

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"PlotSoilProfile: wrote {out}")


# ------------------------------------------------------------
# 4. PDMY02 MATERIAL CARD
# ------------------------------------------------------------


def plot_pdmy02(data: dict, layers: list[dict], out: Path) -> None:
    """
    Write the full PDMY02 argument card versus depth.

    Args:    data, layers, out  destination PNG
    Returns: none (writes PNG; skips when no sand is present)
    """
    zmax = max(float(L["depthBot"]) for L in layers)
    prof = int(data.get("soilProfile", 0))
    sands = [L for L in layers if L.get("sand")]
    if not sands:
        print("PlotSoilProfile: no sand layers — skip PDMY02 figure")
        return

    fig, axes = plt.subplots(3, 4, figsize=(14.5, 11.0), constrained_layout=True)
    fig.suptitle(
        f"PDMY02 material card — profile {prof}\n"
        r"(layer-varying args; defaults: contrac2=5, dilat2=3, liquefac1=1, liquefac2=0; "
        r"$e$ from wiki vs $D_r$)",
        fontsize=12,
    )

    panels = [
        (0, 0, "rho", 1.0, r"$\rho$ (kg/m³)", r"$\rho$", "#37474f"),
        (0, 1, "Gr", 1e-6, r"$G_r$ (MPa)", r"$G_r$", "#1565c0"),
        (0, 2, "Br", 1e-6, r"$B_r$ (MPa)", r"$B_r$", "#c62828"),
        (0, 3, "phi", 1.0, r"$\phi$ (°)", r"$\phi$", "#6a1b9a"),
        (1, 0, "gam_max", 1.0, r"$\gamma_\mathrm{max}$", r"$\gamma_\mathrm{max}$", "#00838f"),
        (1, 1, "pRef", 1e-3, r"$p'_r$ (kPa)", r"$p'_r$", "#5d4037"),
        (1, 2, "d", 1.0, r"$d$ (-)", r"$d$", "#546e7a"),
        (1, 3, "PTA", 1.0, r"PTAng (°)", r"PTAng", "#4527a0"),
    ]
    for r, c, key, scale, xlab, lab, color in panels:
        ax = axes[r, c]
        shade_layers(ax, layers)
        plot_step(ax, layers, key, scale=scale, color=color, label=lab)
        style_depth(ax, zmax)
        ax.set_xlabel(xlab)
        ax.set_title(lab)
        ax.legend(fontsize=8, loc="lower right", framealpha=0.9)

    ax = axes[2, 0]
    shade_layers(ax, layers)
    plot_step(ax, layers, "contr1", color="#1565c0", label="contrac1")
    plot_step(ax, layers, "contr3", color="#c62828", label="contrac3", ls=":")
    ax2 = ax.twiny()
    x, z = step_xy(layers, "contr2")
    ax2.plot(x, z, "--", color="#ef6c00", lw=1.6, label="contrac2 (def)")
    style_depth(ax, zmax)
    ax.set_xlabel("contrac1, contrac3 (-)")
    ax2.set_xlabel("contrac2 (-)")
    ax.set_title("Contraction")
    twin_legend(ax, ax2)

    ax = axes[2, 1]
    shade_layers(ax, layers)
    plot_step(ax, layers, "dilat1", color="#2e7d32", label="dilat1")
    plot_step(ax, layers, "dilat3", color="#ad1457", label="dilat3", ls=":")
    ax2 = ax.twiny()
    x, z = step_xy(layers, "dilat2")
    ax2.plot(x, z, "--", color="#ef6c00", lw=1.6, label="dilat2 (def)")
    style_depth(ax, zmax)
    ax.set_xlabel("dilat1, dilat3 (-)")
    ax2.set_xlabel("dilat2 (-)")
    ax.set_title("Dilation")
    twin_legend(ax, ax2)

    ax = axes[2, 2]
    shade_layers(ax, layers)
    plot_step(ax, layers, "liq1", color="#1565c0", label="liquefac1 (def)")
    plot_step(ax, layers, "liq2", color="#c62828", label="liquefac2 (def)", ls="--")
    plot_step(ax, layers, "e", color="#2e7d32", label=r"$e$ (wiki vs $D_r$)", ls=":")
    style_depth(ax, zmax)
    ax.set_xlabel("(-)")
    ax.set_title(r"liquefac + void ratio $e$")
    ax.legend(fontsize=7, loc="lower right", framealpha=0.9)

    ax = axes[2, 3]
    shade_layers(ax, layers)
    plot_step(ax, layers, "nYS", color="#455a64", label="nYS")
    ax2 = ax.twiny()
    x, z = step_xy(layers, "Dr")
    ax2.plot(x, z, "--", color="#ef6c00", lw=1.8, label=r"$D_r$ (%)")
    style_depth(ax, zmax)
    ax.set_xlabel("nYS (-)")
    ax2.set_xlabel(r"$D_r$ (%)")
    ax.set_title(r"Yield surfaces + $D_r$")
    twin_legend(ax, ax2)

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"PlotSoilProfile: wrote {out}")


# ------------------------------------------------------------
# 5. COMMAND-LINE ENTRY POINT
# ------------------------------------------------------------


def main() -> int:
    """
    Read profile data and write both soil-profile figures.

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
    layers = data.get("layers", [])
    if not layers:
        print("no layers in JSON", file=sys.stderr)
        return 1

    if len(sys.argv) > 2:
        out_dir = Path(sys.argv[2])
    else:
        out_dir = soil_profile_dir(data.get("soilProfile", 1))
    out_dir.mkdir(parents=True, exist_ok=True)

    plot_overview(data, layers, out_dir / "overview.png")
    plot_pdmy02(data, layers, out_dir / "pdmy02.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
