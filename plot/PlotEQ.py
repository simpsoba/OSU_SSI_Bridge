#!/usr/bin/env python3
"""Goals
-----
Post-process one serial OpenSees earthquake-window recorder dump.
Use the ``DO_*`` switches to select history, depth-stacked pile/soil ux,
envelope, hysteresis, quad, spring, pile-section, and frame plots. Write
PNGs to ``LOCAL/plots/runs/<Test>/eq/`` for lab dumps, else
``<eqOutDir>/plots/``. Use ``PlotEQParallel.py`` for OpenSeesMP dumps
whose files end in ``.$pid``.

Units: N, m, s.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.gridspec import GridSpec
from matplotlib.tri import Triangulation
from matplotlib.colors import to_rgba
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, MaxNLocator

from lab_paths import (
    CYLINDER_LENGTH_SCALE,
    M_TO_MM,
    TIME_SCALE_FROUDE,
    XLIM_FULL_PROTO_S,
    YLIM_DISP_PROTO_MM,
    is_lab_dump,
    run_eq_plots_dir,
)
from paths import HERE, elevation_dir, eq_compare_dir, eq_dir, pile_springs_dir
from PlotModelSketch import layer_style
from gm_duration import arias_significant_duration

# ------------------------------------------------------------
# EDIT
# ------------------------------------------------------------
EQ_OUT = eq_dir(3, "Shin", "quad", "forceBeamColumn")

DO_HIST = 1
DO_DEPTH_HIST = 1     # pile + soil ux vs t, one axes per depth
DO_ENVELOPE = 1       # pile-node ux min/max (symmetric)
DO_SPRING_ENV = 1     # spring defo vs y50/z50, force vs pult/tult/qult
DO_HYST = 1           # all p-y, t-z, and q-z loops
DO_QUAD_PEAK = 1      # window peak |tau_xy| and |gamma_xy|
DO_QUAD_HYST = 1      # tau_xy vs gamma_xy vs depth (center + near-FF columns)
DO_HINGE = 1          # pier hinge hist + M-rot, P-axial, P-M, axial-rot
DO_PILE_SEC = 1       # pile M-kappa hyst; peak M and kappa vs depth
DO_FRAMES = 0         # 0=auto (on when every window node has disp); 1=force; -1=off
DO_FRAME_HIST = 1     # ux history side panel on animation frames

N_FRAMES = 0          # >0 = that many equally spaced; 0 = use FRAME_FPS
FRAME_FPS = 30        # max PNGs per second of analysis (0 = every recorded sample)
SCALE = 20.0
SUBTRACT_T0 = 1       # nodal plots/frames relative to first sample
HYST_QUAD_X = None    # fallback if no center/ff tags; None = outermost |x|
PIER_TOP = 5
PIER_BOT = 1
DPI = 140
FRAME_DPI = 80
MOVIE_FPS = 30        # encode at this fps; duration = model-scale window length
FRAME_T0_MODEL_S = 16.0  # skip frames before this lab/model time (s)
SPRING_MINLEN = 0.30  # m, glyph floor so coincident zeroLength springs stay visible
FRAME_HIST_DEPTH_M = (0.0, 2.5, 5.0, 10.0)  # pile/soil stations below grade
FRAME_FS = 8          # frame tick / annotation / legend
FRAME_FS_AXIS = 9     # frame axis labels and mesh title
FRAME_FS_SUP = 10     # frame history supylabel
PEAK_GAMMA_VMAX = 10.0  # window_peak_gamma_xy colorbar top (×10⁻³ units)
# ------------------------------------------------------------

GRAY = "#90a4ae"
ORANGE = "#c45c12"
BLUE = "#1565c0"
BROWN = "#8B5A2B"
PURPLE = "#6a1b9a"
GREEN = "#2e7d32"
NAMES = ("L", "C", "R")
COLORS = {"L": BLUE, "C": BROWN, "R": PURPLE}


# ------------------------------------------------------------
# 1. INPUT AND RECORDER I/O
# ------------------------------------------------------------

def _skip_hash(path: Path) -> list[str]:
    """Read nonblank, noncomment lines from a text file.

    Args:
        path: Input text-file path.
    Returns:
        Lines with blanks and ``#`` comments removed.
    """
    lines = []
    for ln in path.read_text().splitlines():
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        lines.append(s)
    return lines


def read_meta(eq: Path) -> dict[str, str]:
    """Read key-value pairs from the serial window metadata file.

    Args:
        eq: Serial recorder directory.
    Returns:
        Metadata values keyed by their Tcl names.
    """
    meta = {}
    for ln in _skip_hash(eq / "window_meta.txt"):
        k, _, rest = ln.partition(" ")
        meta[k] = rest.strip()
    return meta


def node_tag_bases(meta: dict | None) -> tuple[int, int, int, int, int]:
    """Get node-tag ranges used to classify model parts.

    Args:
        meta: Window metadata, or ``None`` for legacy defaults.
    Returns:
        Soil base, spring base, soffit offset, boundary base, and last soil tag.
    """
    meta = meta or {}

    def gi(key: str, default: int) -> int:
        """Convert one metadata value to an integer.

        Args:
            key: Metadata key.
            default: Value used when the key is absent or invalid.
        Returns:
            Parsed integer value.
        """
        try:
            return int(float(meta.get(key, default)))
        except (TypeError, ValueError):
            return default

    soil = gi("tagShift_soil", 10000)
    spr = gi("nodeTag_sprSoil_base", 20000)
    return (
        soil,
        spr,
        gi("sprSoffitOff", 920),
        gi("nodeTag_bnd_base", 30000),
        gi("soilNodeLast", spr - 1),
    )


def read_node_file(path: Path) -> tuple[list[int], dict[int, tuple[float, float]]]:
    """Read node tags and planar coordinates.

    Args:
        path: Node-list file path; coordinates are in m.
    Returns:
        Ordered tags and a tag-to-``(x, y)`` coordinate map in m.
    """
    tags: list[int] = []
    xy: dict[int, tuple[float, float]] = {}
    for ln in _skip_hash(path):
        a, b, c = ln.split()[:3]
        t = int(a)
        tags.append(t)
        xy[t] = (float(b), float(c))
    return tags, xy


def read_nodes(eq: Path) -> tuple[list[int], dict[int, tuple[float, float]]]:
    """Read every node listed in ``window_nodes.txt``.

    Args:
        eq: Serial recorder directory.
    Returns:
        Ordered tags and a tag-to-``(x, y)`` coordinate map in m.
    """
    return read_node_file(eq / "window_nodes.txt")


def read_disp_nodes(eq: Path) -> list[int]:
    """Read the node order used by displacement recorder columns.

    Args:
        eq: Serial recorder directory.
    Returns:
        Node tags in ``window_disp*.out`` column order.

    Lean dumps record only pier and center-pile displacements. Legacy dumps
    without ``disp_nodes.txt`` use every window node.
    """
    p = eq / "disp_nodes.txt"
    if not p.is_file():
        return read_nodes(eq)[0]
    return read_node_file(p)[0]


def read_eles(eq: Path) -> tuple[list[list[int]], list[list[int]]]:
    """Read line and quadrilateral connectivity from the window element list.

    Args:
        eq: Serial recorder directory.
    Returns:
        Two-node line connectivity and four-node quad connectivity.
    """
    lines = []
    quads = []
    for ln in _skip_hash(eq / "window_eles.txt"):
        toks = [int(x) for x in ln.split()]
        nodes = toks[1:]
        if len(nodes) >= 4:
            quads.append(nodes[:4])
        elif len(nodes) == 2:
            lines.append(nodes)
    return lines, quads


def load_window_disp(
    eq: Path, tags: list[int], disp_files: list[str]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load and stitch serial nodal displacement recorder files.

    Args:
        eq: Serial recorder directory.
        tags: Nodes expected in recorder-column order.
        disp_files: Recorder filenames to stitch.
    Returns:
        Time in s, horizontal displacement in m, and vertical displacement in m.
    """
    ux_parts = []
    uy_parts = []
    t = None
    for fn in disp_files:
        a = loadtxt_partial(eq / fn)
        if a.size == 0:
            continue
        if a.ndim == 1:
            a = a.reshape(1, -1)
        ti = a[:, 0]
        if t is None:
            t = ti
        data = a[:, 1:]
        ux_parts.append(data[:, 0::2])
        uy_parts.append(data[:, 1::2])
    if not ux_parts:
        raise SystemExit(f"PlotEQ: no data in disp files {disp_files}")
    n = min(p.shape[0] for p in ux_parts)
    if t is not None and len(t) != n:
        print(f"PlotEQ: WARNING window disp row mismatch; trim to {n}")
        t = t[:n]
    ux = np.hstack([p[:n] for p in ux_parts])
    uy = np.hstack([p[:n] for p in uy_parts])
    if ux.shape[1] != len(tags):
        raise SystemExit(
            f"PlotEQ: window disp cols {ux.shape[1]} != nDispNodes {len(tags)}"
        )
    return t, ux, uy


def loadtxt_partial(path: Path) -> np.ndarray:
    """Load numeric rows while tolerating one truncated final recorder line.

    Args:
        path: Numeric recorder-file path.
    Returns:
        Two-dimensional floating-point data array.
    """
    try:
        a = np.loadtxt(path)
    except ValueError:
        rows: list[list[float]] = []
        ncol: int | None = None
        with path.open() as f:
            for ln in f:
                s = ln.strip()
                if not s or s.startswith("#"):
                    continue
                try:
                    row = [float(x) for x in s.split()]
                except ValueError:
                    break
                if ncol is None:
                    ncol = len(row)
                if len(row) != ncol:
                    break
                rows.append(row)
        a = np.asarray(rows, dtype=float)
    if a.size == 0:
        return np.empty((0, 0))
    if a.ndim == 1:
        a = a.reshape(1, -1)
    return a


def load_spring_pt(
    eq: Path, force_name: str, defo_name: str, n_ele: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load paired spring force and deformation recorders.

    Args:
        eq: Serial recorder directory.
        force_name: Spring-force recorder filename; forces are in N.
        defo_name: Spring-deformation recorder filename; deformations are in m.
        n_ele: Number of recorded spring elements.
    Returns:
        Time in s, force array in N, and deformation array in m.
    """
    f = loadtxt_partial(eq / force_name)
    d = loadtxt_partial(eq / defo_name)
    if f.ndim == 1:
        f = f.reshape(1, -1)
        d = d.reshape(1, -1)
    ncomp = (f.shape[1] - 1) // n_ele
    F = f[:, 1:].reshape(f.shape[0], n_ele, ncomp)
    U = d[:, 1:].reshape(d.shape[0], n_ele, ncomp)
    return f[:, 0], F, U


# cap BL/BC/BR (Parameters.tcl): pile heads, shared with the cap soffit
PILE_HEAD_TAGS = {1027, 1028, 1029}


# ------------------------------------------------------------
# 2. GEOMETRY AND COMMON PLOT HELPERS
# ------------------------------------------------------------

def pile_groups(xy: dict[int, tuple[float, float]]) -> dict[str, list[int]]:
    """Group pile nodes by shaft and order each shaft from head to tip.

    Args:
        xy: Node coordinates in m keyed by tag.
    Returns:
        Pile names mapped to node tags ordered by decreasing y.
    """
    piles = [t for t in xy if 2000 <= t < 3000 or t in PILE_HEAD_TAGS]
    if not piles:
        return {}
    xs = sorted({round(xy[t][0], 2) for t in piles})
    if len(xs) == 1:
        names = ["C"]
    elif len(xs) <= 3:
        names = list(NAMES[: len(xs)])
    else:
        names = [f"p{i}" for i in range(len(xs))]
    out: dict[str, list[int]] = {n: [] for n in names}
    for t in piles:
        x = round(xy[t][0], 2)
        i = min(range(len(xs)), key=lambda k: abs(xs[k] - x))
        out[names[i]].append(t)
    for n in out:
        out[n].sort(key=lambda t: xy[t][1], reverse=True)
    return out


def maybe_t0(arr: np.ndarray) -> np.ndarray:
    """Subtract the first sample when ``SUBTRACT_T0`` is enabled.

    Args:
        arr: Time-history array in its native units.
    Returns:
        Original or first-sample-relative array in the same units.
    """
    if not SUBTRACT_T0 or arr.shape[0] < 1:
        return arr
    return arr - arr[0]


def sym_xlim(ax, *vals: np.ndarray, pad: float = 1.05) -> None:
    """Set symmetric x limits around all finite input values.

    Args:
        ax: Matplotlib axes to update.
        *vals: Arrays in the plotted x-axis units.
        pad: Multiplicative limit padding.
    Returns:
        None.
    """
    m = 0.0
    for v in vals:
        if v is None or len(np.atleast_1d(v)) == 0:
            continue
        a = np.asarray(v, dtype=float)
        if not np.isfinite(a).any():
            continue
        m = max(m, float(np.nanmax(np.abs(a))))
    if not np.isfinite(m) or m <= 0:
        m = 1.0
    ax.set_xlim(-pad * m, pad * m)


def layer_at_y(y: float, layers: list[dict]) -> str:
    """Find the soil-layer name at one elevation.

    Args:
        y: Elevation in m.
        layers: Layer dictionaries with vertical bounds in m.
    Returns:
        Soil-layer name.
    """
    if not layers:
        return "soil"

    def span(L: dict) -> tuple[float, float]:
        """Normalize one layer's vertical bounds.

        Args:
            L: Soil-layer dictionary with bounds in m.
        Returns:
            Bottom and top elevations in m.
        """
        if "yBot" in L and "yTop" in L:
            a, b = float(L["yBot"]), float(L["yTop"])
        else:
            a, b = float(L["y0"]), float(L["y1"])
        return min(a, b), max(a, b)

    for L in layers:
        y0, y1 = span(L)
        if y0 - 1e-6 <= y <= y1 + 1e-6:
            return str(L["name"])
    if y > max(span(L)[1] for L in layers):
        return str(layers[0]["name"])
    return str(layers[-1]["name"])


def _pair_idx(pairs: list[tuple[int, int]], idx: dict[int, int]) -> np.ndarray:
    """Map node-tag pairs to displacement-column index pairs.

    Args:
        pairs: Node-tag pairs.
        idx: Node tag to recorder-column index.
    Returns:
        Integer ``(n, 2)`` index array.
    """
    rows = []
    for a, b in pairs:
        if a in idx and b in idx:
            rows.append((idx[a], idx[b]))
    if not rows:
        return np.zeros((0, 2), dtype=int)
    return np.asarray(rows, dtype=int)


def _spring_segments(
    X: np.ndarray, Y: np.ndarray, ij: np.ndarray, fallback: np.ndarray, minlen: float
) -> np.ndarray:
    """Build visible spring segments, including ticks for zero-length springs.

    Args:
        X: Deformed x coordinates in m.
        Y: Deformed y coordinates in m.
        ij: Spring endpoint index pairs.
        fallback: Unit directions for coincident endpoints.
        minlen: Minimum displayed spring length in m.
    Returns:
        Segment endpoints with shape ``(n, 2, 2)`` in m.
    """
    if len(ij) == 0:
        return np.zeros((0, 2, 2))
    p0 = np.column_stack((X[ij[:, 0]], Y[ij[:, 0]]))
    p1 = np.column_stack((X[ij[:, 1]], Y[ij[:, 1]]))
    d = p1 - p0
    L = np.linalg.norm(d, axis=1)
    mid = 0.5 * (p0 + p1)
    u = np.empty_like(d)
    ok = L >= 1.0e-12
    u[ok] = d[ok] / L[ok, None]
    u[~ok] = fallback[~ok]
    half = np.where(L[:, None] < minlen, 0.5 * minlen * u, 0.5 * d)
    return np.stack([mid - half, mid + half], axis=1)


def _line_segments(X: np.ndarray, Y: np.ndarray, ij: np.ndarray) -> np.ndarray:
    """Build line-segment endpoints from coordinate arrays.

    Args:
        X: Deformed x coordinates in m.
        Y: Deformed y coordinates in m.
        ij: Endpoint index pairs.
    Returns:
        Segment endpoints with shape ``(n, 2, 2)`` in m.
    """
    if len(ij) == 0:
        return np.zeros((0, 2, 2))
    a = np.column_stack((X[ij[:, 0]], Y[ij[:, 0]]))
    b = np.column_stack((X[ij[:, 1]], Y[ij[:, 1]]))
    return np.stack([a, b], axis=1)


# ------------------------------------------------------------
# 3. HISTORIES AND PILE ENVELOPES
# ------------------------------------------------------------

def eq_end_time(meta: dict, t: np.ndarray) -> float | None:
    """Infer the earthquake end from record duration and free vibration.

    Args:
        meta: Window metadata with ``Trec`` and ``freeVibT`` in s.
        t: Recorder times in s.
    Returns:
        Earthquake end time in s, or ``None`` when unavailable.
    """
    try:
        fv = float(meta.get("freeVibT", 0) or 0)
    except ValueError:
        return None
    if fv <= 0 or len(t) < 2:
        return None
    try:
        trec = float(meta.get("Trec", t[-1]))
    except ValueError:
        trec = float(t[-1])
    te = trec - fv
    if te <= 0 or te >= float(t[-1]):
        return None
    return te


def truncated_end(meta: dict, t: np.ndarray | None) -> float | None:
    """Find the last sample when a recorder dump stops before ``Trec``.

    Args:
        meta: Window metadata with ``Trec`` in s.
        t: Recorder times in s, or ``None``.
    Returns:
        Truncated end time in s, or ``None`` for a complete/unknown run.
    """
    if t is None or len(t) < 2:
        return None
    try:
        trec = float(meta.get("Trec", 0) or 0)
    except ValueError:
        return None
    t_last = float(t[-1])
    if trec > 0.0 and t_last < trec - 1.0e-3:
        return t_last
    return None


def mark_eq_end(ax, t_eq: float | None) -> None:
    """Mark the earthquake end on a time-history axes.

    Args:
        ax: Matplotlib axes to update.
        t_eq: Earthquake end time in s, or ``None``.
    Returns:
        None.
    """
    if t_eq is None:
        return
    ax.axvline(t_eq, color="#78909c", lw=1.0, ls=":", label="EQ end")


def mark_last_sample(
    ax, t, t_eq: float | None, t_cut: float | None = None
) -> None:
    """Mark earthquake end and an incomplete recorder's last sample.

    Args:
        ax: Matplotlib axes to update.
        t: Recorder times in s.
        t_eq: Earthquake end time in s, or ``None``.
        t_cut: Known truncated end time in s, or ``None``.
    Returns:
        None.
    """
    mark_eq_end(ax, t_eq)
    if t_cut is not None:
        ax.axvline(float(t_cut), color="#c62828", lw=1.0, ls="--",
                   label=f"last sample t={float(t_cut):.2f} s")
        return
    if t_eq is None and t is not None and len(t) > 2 and float(t[-1]) < 80.0:
        ax.axvline(float(t[-1]), color="#c62828", lw=1.0, ls="--",
                   label=f"last sample t={float(t[-1]):.2f} s")


def load_d595_proto() -> tuple[float, float] | None:
    """GM D5–95 bounds on the prototype / OpenSees ``t_num`` clock (gmStart≈0).

    Returns:
        ``(t5_s, t95_s)`` or ``None`` if the VT2 is missing.
    """
    try:
        d = arias_significant_duration()
    except (OSError, ValueError) as exc:
        print(f"PlotEQ: D5-95 unavailable ({exc})")
        return None
    return float(d.t5_s), float(d.t95_s)


def finish_full_zoom_pair(
    ax_f,
    ax_z,
    d595: tuple[float, float] | None,
    t,
    t_eq: float | None,
    t_cut: float | None,
    *,
    xlabel: bool = True,
    zoom_title: bool = True,
    full_xlim: tuple[float, float] | None = XLIM_FULL_PROTO_S,
    full_ylim: tuple[float, float] | None = None,
) -> None:
    """Shared grid/markers; full-panel window + D5–95 xlim on the zoom axes.

    Args:
        ax_f: Full-history axes.
        ax_z: Zoom axes.
        d595: ``(t5, t95)`` in s, or ``None``.
        t, t_eq, t_cut: Recorder time and markers.
        xlabel: Label both bottom x axes.
        zoom_title: Title the zoom panel.
        full_xlim: Prototype-time window for the full panel (default 0–300 s model).
        full_ylim: Optional shared y limits (e.g. ±200 mm for ux).
    Returns:
        None.
    """
    for ax in (ax_f, ax_z):
        mark_last_sample(ax, t, t_eq, t_cut)
        ax.grid(True, ls=":", alpha=0.45)
    if full_xlim is not None:
        ax_f.set_xlim(*full_xlim)
    if full_ylim is not None:
        ax_f.set_ylim(*full_ylim)
    if d595 is not None:
        ax_z.set_xlim(d595[0], d595[1])
        if zoom_title:
            ax_z.set_title(r"D5–95 zoom", fontsize=9, pad=4)
    if xlabel:
        ax_f.set_xlabel(r"$t_\mathrm{num}$ (s)")
        ax_z.set_xlabel(r"$t_\mathrm{num}$ (s)")


def subplots_full_zoom(
    n_rows: int,
    *,
    fig_h: float,
    fig_w: float = 12.8,
    sharey: str | bool = "row",
    wspace: float = 0.04,
    hspace: float = 0.04,
):
    """Create ``n_rows × 2`` axes: full history | D5–95 zoom.

    Uses ``layout='constrained'`` (not a separate ``constrained_layout=``
    flag). Do not put ``wspace`` in ``gridspec_kw`` — that fights the
    layout engine and leaves a wide gutter between columns.

    Args:
        n_rows: Number of signal rows.
        fig_h: Figure height in inches.
        fig_w: Figure width in inches.
        sharey: Matplotlib ``sharey`` for the grid.
        wspace, hspace: Constrained-layout gaps as a fraction of subplot size.
    Returns:
        ``(fig, axes_full, axes_zoom)`` — lists of length ``n_rows``.
    """
    fig, axes = plt.subplots(
        n_rows,
        2,
        figsize=(fig_w, fig_h),
        sharex="col",
        sharey=sharey,
        layout="constrained",
        gridspec_kw={"width_ratios": [1.35, 1.0]},
    )
    engine = fig.get_layout_engine()
    if engine is not None and hasattr(engine, "set"):
        engine.set(w_pad=0.02, h_pad=0.02, wspace=wspace, hspace=hspace)
    if n_rows == 1:
        return fig, [axes[0]], [axes[1]]
    return fig, list(axes[:, 0]), list(axes[:, 1])


def hyst_loop(ax, x, y, xlabel: str, ylabel: str, title: str) -> None:
    """Draw one labeled hysteresis loop.

    Args:
        ax: Matplotlib axes to update.
        x: Deformation-like history in label-defined units.
        y: Force-like history in label-defined units.
        xlabel: Horizontal-axis label.
        ylabel: Vertical-axis label.
        title: Panel title.
    Returns:
        None.
    """
    ax.plot(x, y, color=BROWN, lw=0.85, rasterized=True)
    ax.axhline(0.0, color="#9e9e9e", lw=0.6)
    ax.axvline(0.0, color="#9e9e9e", lw=0.6)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, ls=":", alpha=0.45)


def plot_hist(
    out: Path,
    t,
    ux,
    uy,
    idx,
    groups,
    t_eq: float | None = None,
    t_cut: float | None = None,
    d595: tuple[float, float] | None = None,
) -> None:
    """Plot pier and pile-head displacement histories (full | D5–95 zoom).

    Args:
        out: Plot output directory.
        t: Recorder times in s.
        ux: Horizontal nodal displacements in m.
        uy: Vertical nodal displacements in m.
        idx: Node tag to displacement-column index.
        groups: Pile names mapped to ordered node tags.
        t_eq: Earthquake end time in s, or ``None``.
        t_cut: Truncated recorder time in s, or ``None``.
        d595: ``(t5, t95)`` prototype s for the zoom panel, or ``None``.
    Returns:
        None; writes ``hist_ux.png`` and ``hist_uy.png``.
    """

    def _draw_ux(ax) -> None:
        if PIER_TOP in idx:
            ax.plot(
                t,
                to_mm(ux[:, idx[PIER_TOP]]),
                color=ORANGE,
                lw=1.4,
                label=f"pier top ({PIER_TOP}) ux",
            )
        if PIER_BOT in idx:
            ax.plot(
                t,
                to_mm(ux[:, idx[PIER_BOT]]),
                color=ORANGE,
                lw=1.0,
                ls="--",
                label=f"pier bot ({PIER_BOT}) ux",
            )
        for name, tags in groups.items():
            if not tags:
                continue
            ax.plot(
                t,
                to_mm(ux[:, idx[tags[0]]]),
                color=COLORS.get(name, "#333"),
                lw=1.0,
                label=f"pile {name} head ux",
            )

    def _draw_uy(ax) -> None:
        if PIER_TOP in idx:
            ax.plot(t, to_mm(uy[:, idx[PIER_TOP]]), color=ORANGE, lw=1.4, label="pier top uy")
        if PIER_BOT in idx:
            ax.plot(
                t,
                to_mm(uy[:, idx[PIER_BOT]]),
                color=ORANGE,
                lw=1.0,
                ls="--",
                label="pier bot uy",
            )

    if d595 is None:
        fig, ax = plt.subplots(figsize=(10.4, 4.2), constrained_layout=True)
        _draw_ux(ax)
        mark_last_sample(ax, t, t_eq, t_cut)
        ax.set_xlabel(r"$t_\mathrm{num}$ (s)")
        ax.set_ylabel("ux (mm)")
        ax.legend(fontsize=8, ncol=2)
        ax.grid(True, ls=":", alpha=0.45)
        fig.savefig(out / "hist_ux.png", dpi=DPI)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(10.4, 3.6), constrained_layout=True)
        _draw_uy(ax)
        mark_last_sample(ax, t, t_eq, t_cut)
        ax.set_xlabel(r"$t_\mathrm{num}$ (s)")
        ax.set_ylabel("uy (mm)")
        ax.legend(fontsize=8)
        ax.grid(True, ls=":", alpha=0.45)
        fig.savefig(out / "hist_uy.png", dpi=DPI)
        plt.close(fig)
        return

    fig, axes_f, axes_z = subplots_full_zoom(1, fig_h=4.2, sharey=True)
    ax_f, ax_z = axes_f[0], axes_z[0]
    _draw_ux(ax_f)
    _draw_ux(ax_z)
    finish_full_zoom_pair(
        ax_f, ax_z, d595, t, t_eq, t_cut, full_ylim=YLIM_DISP_PROTO_MM
    )
    ax_f.set_ylabel("ux (mm)")
    ax_f.legend(fontsize=8, ncol=2)
    fig.savefig(out / "hist_ux.png", dpi=DPI)
    plt.close(fig)

    fig, axes_f, axes_z = subplots_full_zoom(1, fig_h=3.6, sharey=True)
    ax_f, ax_z = axes_f[0], axes_z[0]
    _draw_uy(ax_f)
    _draw_uy(ax_z)
    finish_full_zoom_pair(
        ax_f, ax_z, d595, t, t_eq, t_cut, full_ylim=YLIM_DISP_PROTO_MM
    )
    ax_f.set_ylabel("uy (mm)")
    ax_f.legend(fontsize=8)
    fig.savefig(out / "hist_uy.png", dpi=DPI)
    plt.close(fig)


def preferred_pile_tags(groups: dict[str, list[int]]) -> tuple[str, list[int]]:
    """Pick the center pile shaft when present, else the first nonempty group.

    Args:
        groups: Pile names mapped to head-to-tip node tags.
    Returns:
        ``(name, tags)`` or ``("", [])`` when empty.
    """
    if groups.get("C"):
        return "C", groups["C"]
    for name in NAMES:
        if groups.get(name):
            return name, groups[name]
    for name, tags in groups.items():
        if tags:
            return name, tags
    return "", []


def soil_column_nodes(
    tags: list[int],
    xy: dict[int, tuple[float, float]],
    meta: dict,
    x_tgt: float = 0.0,
) -> list[tuple[int, float]]:
    """Soil continuum nodes along one vertical column (nearest ``x_tgt``).

    Args:
        tags: Candidate node tags (usually ``disp_nodes``).
        xy: Coordinates in m keyed by tag.
        meta: Window metadata for soil tag bounds.
        x_tgt: Target column x in m (center pile ≈ 0).
    Returns:
        ``(tag, y)`` pairs ordered head-to-tip (decreasing y).
    """
    soil_base, _spr, _soff, _bnd, soil_last = node_tag_bases(meta)
    cands = [
        t
        for t in tags
        if t in xy and soil_base <= t <= soil_last
    ]
    if not cands:
        return []
    # Prefer a narrow band around x_tgt; widen if the lean/full window is coarse.
    for x_tol in (0.35, 1.0, 3.0):
        near = [t for t in cands if abs(xy[t][0] - x_tgt) <= x_tol]
        if near:
            cands = near
            break
    by_y: dict[float, int] = {}
    for t in cands:
        x, y = xy[t]
        yk = round(y, 3)
        prev = by_y.get(yk)
        if prev is None or abs(x - x_tgt) < abs(xy[prev][0] - x_tgt):
            by_y[yk] = t
    return [
        (by_y[yk], float(yk))
        for yk in sorted(by_y.keys(), reverse=True)
    ]


def _interp_ux(t_ref: np.ndarray, t_src: np.ndarray, u_src: np.ndarray) -> np.ndarray:
    """Interpolate one displacement history onto a reference time grid.

    Args:
        t_ref: Target times in s.
        t_src: Source times in s.
        u_src: Source displacements in m.
    Returns:
        Displacements in m on ``t_ref``.
    """
    if len(t_src) == len(t_ref) and np.allclose(t_src, t_ref, atol=1e-6, rtol=0):
        return np.asarray(u_src, dtype=float)
    return np.interp(t_ref, t_src, u_src)


def estimate_soil_ux_column(
    eq: Path,
    meta: dict,
    js: dict | None,
    t_pile: np.ndarray,
    pile_tags: list[int],
    xy: dict[int, tuple[float, float]],
    idx: dict[int, int],
    ux: np.ndarray,
    ip_prefer: int = 1,
) -> list[tuple[float, np.ndarray, str]] | None:
    """Estimate soil ux(t) at SSI horizons: u_soil ≈ u_pile − u_py.

    Spring is ``zeroLength soil→dup`` with dup equalDOF'd to the pile, so
    recorded py deformation is u_dup − u_soil ≈ u_pile − u_soil. UX has no
    fold datum; SUBTRACT_T0 on both series removes a constant uy-style offset
    if the caller already zeroed the pile history.

    Args:
        eq: Serial (or stitched) dump directory.
        meta: Window metadata.
        js: Spring JSON, or ``None``.
        t_pile: Pile recorder times in s.
        pile_tags: Center-pile tags head-to-tip.
        xy: Coordinates in m.
        idx: Disp-column index by node tag.
        ux: Nodal ux in m (already t0-relative if SUBTRACT_T0).
        ip_prefer: Zero-based pile index (1 = center).
    Returns:
        ``(y, u_soil, label)`` rows tipward, or ``None`` when unavailable.
    """
    rows = read_pile_spring_eles(eq)
    if not rows:
        return None
    stations = stations_for_plot(rows, js, meta)
    if not stations or len(stations) != len(rows):
        return None
    try:
        t_spr, _F, U = load_spring_pt(
            eq, "pile_springs_force.out", "pile_springs_defo.out", len(rows)
        )
    except Exception:
        return None
    if U.ndim != 3 or U.shape[2] < 1 or len(t_spr) < 2:
        return None

    pile_pts = [
        (t, xy[t][1], idx[t])
        for t in pile_tags
        if t in xy and t in idx
    ]
    if not pile_pts:
        return None

    out_rows: list[tuple[float, np.ndarray, str]] = []
    for i, st in enumerate(stations):
        if int(st.get("ip", -1)) != ip_prefer:
            continue
        if is_qz_station(st):
            # Tip q-z is axial (dir 2); skip for lateral soil estimate.
            continue
        y_s = float(st["y"])
        # Nearest pile node in elevation.
        tag, y_p, icol = min(pile_pts, key=lambda p: abs(p[1] - y_s))
        if abs(y_p - y_s) > 0.25:
            continue
        u_py = maybe_t0(U[:, i, 0])
        u_py_i = _interp_ux(t_pile, t_spr, u_py)
        u_pile = ux[:, icol]
        u_soil = u_pile - u_py_i
        out_rows.append((y_s, u_soil, f"y={y_s:.2f}"))
    if not out_rows:
        return None
    out_rows.sort(key=lambda r: r[0], reverse=True)
    return out_rows


def format_depth_label(y: float) -> str:
    """Depth below grade ``|y|`` with three significant figures.

    Args:
        y: Elevation in m (down is negative).
    Returns:
        Label such as ``0.00``, ``0.99``, ``1.91``, ``10.1``.
    """
    from math import floor, log10

    d = abs(float(y))
    if d < 0.05:
        return "0.00"
    if d < 1.0:
        return f"{d:.2f}"
    order = int(floor(log10(d)))
    decimals = max(0, 2 - order)  # 3 significant figures in fixed point
    return f"{d:.{decimals}f}"


def foundation_cap_ux(
    idx: dict[int, int],
    xy: dict[int, tuple[float, float]],
    ux: np.ndarray,
    pile_tags: list[int],
    meta: dict,
) -> np.ndarray | None:
    """Lateral motion at the pile head / cap BC (tag 1028 when present).

    Args:
        idx, xy, ux: Disp maps / histories (m).
        pile_tags: Center-pile tags head-to-tip.
        meta: Window metadata for soil tag bound.
    Returns:
        ux history in m, or ``None``.
    """
    soil_base, *_ = node_tag_bases(meta)
    for tag in (1028, 1025):
        if tag in idx and tag < soil_base:
            return ux[:, idx[tag]]
    if pile_tags:
        head = pile_tags[0]
        if head in idx:
            return ux[:, idx[head]]
    return None


def plot_stacked_depth_hist(
    out: Path,
    fname: str,
    t: np.ndarray,
    rows: list[tuple[str, np.ndarray]],
    color: str,
    t_eq: float | None,
    t_cut: float | None,
    *,
    qty_label: str,
    d595: tuple[float, float] | None = None,
) -> None:
    """Depth-stacked time histories: full | D5–95, no figure titles.

    Every axes shows ytick values and a depth ylabel; ``fig.supylabel`` is the
    plotted quantity (e.g. ``$u_x$``).

    Args:
        out: Plot output directory.
        fname: Output PNG name.
        t: Times in s.
        rows: ``(depth_label, series)`` head-to-tip.
        color: Trace color.
        t_eq, t_cut: Markers.
        qty_label: Shared quantity (e.g. ``$u_x$ (mm)``).
        d595: Zoom window, or ``None`` for a single column.
    Returns:
        None; writes ``fname``.
    """
    if not rows:
        return
    n = len(rows)
    fig_h = max(3.0, 0.48 * n + 0.8)
    all_u = [u for _, u in rows]
    m = 0.0
    for u in all_u:
        a = np.asarray(u, dtype=float)
        if np.isfinite(a).any():
            m = max(m, float(np.nanmax(np.abs(a))))
    ylim_u = m * M_TO_MM * 1.05 if m > 0 else 1e-4 * M_TO_MM

    if d595 is None:
        fig, axes = plt.subplots(n, 1, figsize=(9.2, fig_h), sharex=True)
        axes_f = [axes] if n == 1 else list(axes)
        axes_z = None
        fig.subplots_adjust(
            left=0.16, right=0.99, top=0.995, bottom=0.045, hspace=0.08
        )
    else:
        fig, axes = plt.subplots(
            n,
            2,
            figsize=(11.0, fig_h),
            sharex="col",
            sharey="row",
            gridspec_kw={"width_ratios": [1.45, 1.0]},
        )
        # Room for supylabel + depth ylabel + tick numbers on both columns.
        fig.subplots_adjust(
            left=0.14,
            right=0.995,
            top=0.995,
            bottom=0.045,
            wspace=0.22,
            hspace=0.08,
        )
        axes_f = [axes[0]] if n == 1 else list(axes[:, 0])
        axes_z = [axes[1]] if n == 1 else list(axes[:, 1])

    for i, (lab, u) in enumerate(rows):
        ax_list = (axes_f[i],) if axes_z is None else (axes_f[i], axes_z[i])
        for ax in ax_list:
            ax.plot(t, to_mm(u), color=color, lw=0.9)
            ax.axhline(0.0, color="#bbb", lw=0.6)
            ax.set_ylim(-ylim_u, ylim_u)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
            ax.tick_params(labelsize=7, labelleft=True, pad=1)
            ax.set_ylabel(
                lab, fontsize=8, rotation=0, va="center", ha="right", labelpad=6
            )
            mark_last_sample(ax, t, t_eq, t_cut)
            ax.grid(True, ls=":", alpha=0.4)
        if axes_z is not None and d595 is not None:
            axes_z[i].set_xlim(d595[0], d595[1])
        axes_f[i].set_xlim(*XLIM_FULL_PROTO_S)
        if i == n - 1:
            axes_f[i].set_xlabel(r"$t_\mathrm{num}$ (s)")
            if axes_z is not None:
                axes_z[i].set_xlabel(r"$t_\mathrm{num}$ (s)")

    fig.supylabel(qty_label, fontsize=10)
    fig.savefig(out / fname, dpi=DPI)
    plt.close(fig)


# Back-compat name used by older call sites.
def plot_stacked_ux_depth(
    out: Path,
    fname: str,
    t: np.ndarray,
    rows: list[tuple[str, np.ndarray]],
    color: str,
    t_eq: float | None,
    t_cut: float | None,
    note: str = "",
    d595: tuple[float, float] | None = None,
) -> None:
    """Deprecated wrapper: depth-stacked ``$u_x$`` (``note`` ignored)."""
    del note
    plot_stacked_depth_hist(
        out,
        fname,
        t,
        rows,
        color,
        t_eq,
        t_cut,
        qty_label=r"$u_x$ (mm)",
        d595=d595,
    )


def load_pier_node_rz(
    eq: Path, node_tag: int
) -> tuple[np.ndarray, np.ndarray] | None:
    """Load pier-node RZ from ``pier_node_<tag>.out`` (time, ux, uy, rz).

    Args:
        eq: Dump directory (serial or stitched).
        node_tag: Pier node tag.
    Returns:
        ``(t_s, rz_rad)`` or ``None``.
    """
    cands = [eq / f"pier_node_{node_tag}.out"]
    cands.extend(sorted(eq.glob(f"pier_node_{node_tag}.out.*")))
    path = next((p for p in cands if p.is_file()), None)
    if path is None:
        return None
    a = loadtxt_partial(path)
    if a.size == 0:
        return None
    if a.ndim == 1:
        a = a.reshape(1, -1)
    if a.shape[1] < 4:
        return None
    return a[:, 0], a[:, 3]


def pile_chord_rotations(
    pile_tags: list[int],
    xy: dict[int, tuple[float, float]],
    idx: dict[int, int],
    ux: np.ndarray,
) -> list[tuple[int, np.ndarray]]:
    """Approximate nodal RZ from adjacent pile chord slopes (``Δu_x / Δy``).

    ``window_disp`` records only UX,UY. True RZ is available for pier nodes
    via ``pier_node_*.out``; shaft nodes use this chord estimate.

    Args:
        pile_tags: Head-to-tip tags.
        xy, idx, ux: Geometry and lateral disp (m).
    Returns:
        ``(tag, theta_rad)`` in the same order as ``pile_tags`` (subset).
    """
    tags = [t for t in pile_tags if t in idx and t in xy]
    if len(tags) < 2:
        return []
    ys = np.array([xy[t][1] for t in tags], dtype=float)
    cols = [idx[t] for t in tags]
    out: list[tuple[int, np.ndarray]] = []
    for i, tag in enumerate(tags):
        if i == 0:
            i0, i1 = 0, 1
        elif i == len(tags) - 1:
            i0, i1 = len(tags) - 2, len(tags) - 1
        else:
            i0, i1 = i - 1, i + 1
        dy = float(ys[i0] - ys[i1])
        if abs(dy) < 1e-9:
            th = np.zeros(ux.shape[0], dtype=float)
        else:
            th = (ux[:, cols[i0]] - ux[:, cols[i1]]) / dy
        out.append((tag, th))
    return out


def plot_depth_histories(
    out: Path,
    eq: Path,
    meta: dict,
    js: dict | None,
    t: np.ndarray,
    ux: np.ndarray,
    idx: dict[int, int],
    groups: dict[str, list[int]],
    xy: dict[int, tuple[float, float]],
    disp_tags: list[int],
    t_eq: float | None,
    t_cut: float | None,
    d595: tuple[float, float] | None = None,
) -> None:
    """Write pile/soil ux and pile rotation depth-stacked histories.

    Soil: recorded continuum column when present; else ``u_pile − u_py``.
    Rotation: pier-base RZ from ``pier_node_1`` at grade when present; shaft
    from chord slopes of ``window_disp`` UX (no RZ in that recorder).

    Args:
        out, eq, meta, js, t, ux, idx, groups, xy, disp_tags, t_eq, t_cut, d595
    Returns:
        None; writes depth-stack PNGs when data allow.
    """
    pname, pile_tags = preferred_pile_tags(groups)
    pile_rows: list[tuple[str, np.ndarray]] = []
    # Cap top / pier base at grade (depth 0) when recorded — keep pile head too.
    if 1 in idx and 1 in xy:
        pile_rows.append((format_depth_label(xy[1][1]), ux[:, idx[1]]))
    for tag in pile_tags:
        if tag not in idx or tag not in xy:
            continue
        pile_rows.append((format_depth_label(xy[tag][1]), ux[:, idx[tag]]))
    if pile_rows:
        plot_stacked_depth_hist(
            out,
            "hist_pile_ux_depth.png",
            t,
            pile_rows,
            COLORS.get(pname, BROWN),
            t_eq,
            t_cut,
            qty_label=r"$u_x$ (mm)",
            d595=d595,
        )
        print(f"PlotEQ: wrote {out / 'hist_pile_ux_depth.png'}  ({len(pile_rows)} depths)")
    else:
        print("PlotEQ: no pile shaft in disp_nodes -- skip hist_pile_ux_depth")

    soil_nodes = soil_column_nodes(disp_tags, xy, meta, x_tgt=0.0)
    soil_rows: list[tuple[str, np.ndarray]] = []
    if soil_nodes:
        for tag, y in soil_nodes:
            if tag not in idx:
                continue
            soil_rows.append((format_depth_label(y), ux[:, idx[tag]]))
    if not soil_rows and pile_tags:
        est = estimate_soil_ux_column(
            eq, meta, js, t, pile_tags, xy, idx, ux, ip_prefer=1
        )
        if est:
            soil_rows = [(format_depth_label(y), u) for y, u, _lab in est]
    has_grade = any(lab == "0.00" for lab, _ in soil_rows)
    if soil_rows and not has_grade:
        if 1 in idx:
            soil_rows = [(format_depth_label(xy.get(1, (0.0, 0.0))[1]), ux[:, idx[1]]), *soil_rows]
        else:
            u_cap = foundation_cap_ux(idx, xy, ux, pile_tags, meta)
            if u_cap is not None:
                soil_rows = [("0.00", u_cap), *soil_rows]
    if soil_rows:
        plot_stacked_depth_hist(
            out,
            "hist_soil_ux_depth.png",
            t,
            soil_rows,
            GREEN,
            t_eq,
            t_cut,
            qty_label=r"$u_x$ (mm)",
            d595=d595,
        )
        print(f"PlotEQ: wrote {out / 'hist_soil_ux_depth.png'}  ({len(soil_rows)} depths)")
    else:
        print("PlotEQ: no soil disp / spring estimate -- skip hist_soil_ux_depth")

    # --- rotation vs depth (cap + pile) ---
    rz_rows: list[tuple[str, np.ndarray]] = []
    pier_rz = load_pier_node_rz(eq, 1)
    if pier_rz is not None:
        t_p, rz_p = pier_rz
        rz_grade = maybe_t0(_interp_ux(t, t_p, rz_p))
        rz_rows.append(("0.00", rz_grade))
    chords = pile_chord_rotations(pile_tags, xy, idx, ux)
    for tag, th in chords:
        th = maybe_t0(th)
        # Skip duplicating head if we already placed pier RZ at 0.00 and head
        # elevation is the soffit (~1 m) — still include head chord at its depth.
        rz_rows.append((format_depth_label(xy[tag][1]), th))
    # Deduplicate identical depth labels (keep first).
    seen: set[str] = set()
    rz_uniq: list[tuple[str, np.ndarray]] = []
    for lab, series in rz_rows:
        if lab in seen:
            continue
        seen.add(lab)
        rz_uniq.append((lab, series))
    if rz_uniq:
        plot_stacked_depth_hist(
            out,
            "hist_pile_rz_depth.png",
            t,
            rz_uniq,
            PURPLE,
            t_eq,
            t_cut,
            qty_label=r"$\theta$ (rad)",
            d595=d595,
        )
        print(
            f"PlotEQ: wrote {out / 'hist_pile_rz_depth.png'}  ({len(rz_uniq)} depths)"
            "  [grade=pier_node RZ; shaft=chord dux/dy]"
        )
    else:
        print("PlotEQ: no pile rotation estimate -- skip hist_pile_rz_depth")


def plot_pier_hinge(
    out: Path,
    eq: Path,
    meta: dict,
    t_eq: float | None,
    t_cut: float | None = None,
    d595: tuple[float, float] | None = None,
) -> None:
    """Plot pier-base hinge histories and hysteresis panels.

    Args:
        out: Plot output directory.
        eq: Serial recorder directory.
        meta: Window metadata describing the hinge recorder.
        t_eq: Earthquake end time in s, or ``None``.
        t_cut: Truncated recorder time in s, or ``None``.
        d595: ``(t5, t95)`` for history zoom, or ``None``.
    Returns:
        None; writes hinge PNGs when recorder files exist.
    """
    kind = meta.get("pierHinge", "")
    if kind not in ("lumpedPlasticity", "forceBeamColumn"):
        return
    fn_f = eq / meta.get("hingeForceFile", "pier_hinge_force.out")
    fn_d = eq / meta.get("hingeDefoFile", "pier_hinge_defo.out")
    if not fn_f.is_file() or not fn_d.is_file():
        print("PlotEQ: no pier_hinge_*.out -- skip hinge plots")
        return
    F = loadtxt_partial(fn_f)
    D = loadtxt_partial(fn_d)
    if F.size == 0 or D.size == 0:
        return
    n = min(len(F), len(D))
    t = F[:n, 0]
    P = F[:n, 1]
    M = F[:n, 2] if F.shape[1] > 2 else np.zeros(n)
    ax_d = D[:n, 1]
    rot_d = D[:n, 2] if D.shape[1] > 2 else np.zeros(n)
    if t_eq is None:
        t_eq = eq_end_time(meta, t)
    # ZLS sectionDeformation = axial displacement, rotation (not strain/kappa).
    is_zls = kind == "lumpedPlasticity"
    if is_zls:
        ax_lab = r"$\Delta u_\mathrm{ax}$ (mm)" if SUBTRACT_T0 else r"$u_\mathrm{ax}$ (mm)"
        rot_lab = r"$\Delta\theta$ (rad)" if SUBTRACT_T0 else r"$\theta$ (rad)"
        m_title = r"Mz vs $\theta$"
        p_title = r"P vs $u_\mathrm{ax}$"
        ar_title = r"$u_\mathrm{ax}$ vs $\theta$"
    else:
        ax_lab = r"$\Delta\varepsilon$" if SUBTRACT_T0 else r"$\varepsilon$"
        rot_lab = r"$\Delta\kappa$ (1/m)" if SUBTRACT_T0 else r"$\kappa$ (1/m)"
        m_title = r"Mz vs $\kappa$"
        p_title = r"P vs $\varepsilon$"
        ar_title = r"$\varepsilon$ vs $\kappa$"
    ax_plot = ax_d - ax_d[0] if SUBTRACT_T0 else ax_d
    rot_plot = rot_d - rot_d[0] if SUBTRACT_T0 else rot_d
    if is_zls:
        ax_plot = to_mm(ax_plot)
    P_kN = P / 1.0e3
    M_kNm = M / 1.0e3

    channels = [
        (ax_plot, ax_lab, ORANGE),
        (rot_plot, rot_lab, ORANGE),
        (P_kN, "P (kN)", BLUE),
        (M_kNm, "Mz (kN·m)", BLUE),
    ]
    if d595 is None:
        fig, axes = plt.subplots(
            2, 2, figsize=(10.4, 6.2), sharex="col", constrained_layout=True
        )
        flat = axes.ravel()
        for ax, (y, ylab, col) in zip(flat, channels):
            ax.plot(t, y, color=col, lw=1.2)
            ax.set_ylabel(ylab)
            mark_last_sample(ax, t, t_eq, t_cut)
            ax.grid(True, ls=":", alpha=0.45)
        axes[1, 0].set_xlabel(r"$t_\mathrm{num}$ (s)")
        axes[1, 1].set_xlabel(r"$t_\mathrm{num}$ (s)")
        axes[0, 0].set_title(f"Pier base hinge  ({kind})")
        fig.savefig(out / "hist_hinge.png", dpi=DPI)
        plt.close(fig)
    else:
        fig, axes_f, axes_z = subplots_full_zoom(4, fig_h=8.4, sharey="row")
        for i, (y, ylab, col) in enumerate(channels):
            for ax in (axes_f[i], axes_z[i]):
                ax.plot(t, y, color=col, lw=1.2)
            axes_f[i].set_ylabel(ylab)
            finish_full_zoom_pair(
                axes_f[i],
                axes_z[i],
                d595,
                t,
                t_eq,
                t_cut,
                xlabel=(i == 3),
                zoom_title=(i == 0),
            )
        axes_f[0].set_title(f"Pier base hinge  ({kind})")
        fig.savefig(out / "hist_hinge.png", dpi=DPI)
        plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 8.4), constrained_layout=True)
    hyst_loop(axes[0, 0], rot_plot, M_kNm, rot_lab, "Mz (kN·m)", m_title)
    hyst_loop(axes[0, 1], ax_plot, P_kN, ax_lab, "P (kN)", p_title)
    hyst_loop(axes[1, 0], M_kNm, P_kN, "Mz (kN·m)", "P (kN)", "P vs Mz")
    hyst_loop(axes[1, 1], rot_plot, ax_plot, rot_lab, ax_lab, ar_title)
    fig.suptitle(f"Pier base hinge  ({kind})", fontsize=11)
    fig.savefig(out / "hyst_hinge.png", dpi=DPI)
    plt.close(fig)
    print(f"PlotEQ: wrote {out / 'hist_hinge.png'}  {out / 'hyst_hinge.png'}")


def plot_envelope(out: Path, ux, idx, groups, xy) -> None:
    """Plot pile displacement extrema against depth.

    Args:
        out: Plot output directory.
        ux: Horizontal nodal displacements in m.
        idx: Node tag to displacement-column index.
        groups: Pile names mapped to ordered node tags.
        xy: Node coordinates in m keyed by tag.
    Returns:
        None; writes ``pile_envelope_ux.png``.
    """
    fig, ax = plt.subplots(figsize=(5.6, 7.4), constrained_layout=True)
    xs = []
    for name, tags in groups.items():
        y = np.array([xy[t][1] for t in tags])
        col = np.array([ux[:, idx[t]] for t in tags])
        umin = to_mm(col.min(axis=1))
        umax = to_mm(col.max(axis=1))
        c = COLORS.get(name, "#333")
        ax.plot(umin, y, color=c, lw=1.5, label=f"{name} min")
        ax.plot(umax, y, color=c, lw=1.5, ls="--", label=f"{name} max")
        xs.extend([umin, umax])
    ax.axvline(0.0, color="#bbb", lw=0.8)
    sym_xlim(ax, *xs)
    ax.set_xlabel("ux min / max (mm)")
    ax.set_ylabel("y (m), down is negative")
    ax.legend(fontsize=8)
    ax.grid(True, ls=":", alpha=0.45)
    fig.savefig(out / "pile_envelope_ux.png", dpi=DPI)
    plt.close(fig)


def read_pile_beam_eles(eq: Path) -> list[tuple[int, int, int]]:
    """Read pile beam element, pile, and vertical-station indices.

    Args:
        eq: Serial recorder directory.
    Returns:
        ``(element tag, pile index, station index)`` rows.
    """
    p = eq / "pile_beam_eles.txt"
    if not p.is_file():
        return []
    rows: list[tuple[int, int, int]] = []
    for ln in _skip_hash(p):
        a = ln.split()
        rows.append((int(a[0]), int(a[1]), int(a[2])))
    return rows


def plot_pile_depth_env(
    out: Path,
    fname: str,
    groups: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
    xlabel: str,
) -> None:
    """Plot pile-section extrema against elevation.

    Args:
        out: Plot output directory.
        fname: Output PNG filename.
        groups: Pile names mapped to elevation and min/max arrays.
        xlabel: Horizontal-axis label with units.
    Returns:
        None; writes the named PNG.
    """
    fig, ax = plt.subplots(figsize=(5.6, 7.4), constrained_layout=True)
    xs: list[np.ndarray] = []
    for name, (y, vmin, vmax) in groups.items():
        c = COLORS.get(name, "#333")
        ax.plot(vmin, y, color=c, lw=1.5, label=f"{name} min")
        ax.plot(vmax, y, color=c, lw=1.5, ls="--", label=f"{name} max")
        xs.extend([vmin, vmax])
    ax.axvline(0.0, color="#bbb", lw=0.8)
    sym_xlim(ax, *xs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("y (m), down is negative")
    ax.legend(fontsize=8)
    ax.grid(True, ls=":", alpha=0.45)
    fig.savefig(out / fname, dpi=DPI)
    plt.close(fig)


def plot_pile_section(out: Path, eq: Path, meta: dict, xy: dict) -> None:
    """Plot pile moment/curvature envelopes and hysteresis.

    Args:
        out: Plot output directory.
        eq: Serial recorder directory.
        meta: Window metadata describing pile recorders.
        xy: Node coordinates in m keyed by tag.
    Returns:
        None; writes available pile-section PNGs.
    """
    rows = read_pile_beam_eles(eq)
    fn_g = eq / meta.get("pileBeamGlobalForceFile", "pile_beam_globalForce.out")
    if not rows or not fn_g.is_file():
        print("PlotEQ: no pile_beam_globalForce.out -- skip pile section plots")
        return
    n_ele = len(rows)
    G = loadtxt_partial(fn_g)
    if G.size == 0:
        return
    ncomp = (G.shape[1] - 1) // n_ele
    if ncomp < 6:
        print(f"PlotEQ: pile globalForce ncomp={ncomp}, need 6 -- skip")
        return
    G3 = G[:, 1 : 1 + n_ele * ncomp].reshape(len(G), n_ele, ncomp)
    M_i = G3[:, :, 2] / 1.0e3
    M_j = G3[:, :, 5] / 1.0e3
    ev = read_ele_nodes(eq)
    y_i = np.full(n_ele, np.nan)
    y_j = np.full(n_ele, np.nan)
    ip_of = np.array([r[1] for r in rows], dtype=int)
    iy_of = np.array([r[2] for r in rows], dtype=int)
    for k, (e, ip, iy) in enumerate(rows):
        nd = ev.get(e, [])
        if len(nd) >= 2 and nd[0] in xy and nd[1] in xy:
            y_i[k] = xy[nd[0]][1]
            y_j[k] = xy[nd[1]][1]
    fn_s = eq / meta.get("pileBeamSec1DefoFile", "pile_beam_sec1_defo.out")
    kap = None
    if meta.get("pileEleType") == "dispBeamColumn" and fn_s.is_file():
        S = loadtxt_partial(fn_s)
        if S.size:
            ns = min(len(G), len(S))
            M_i = M_i[:ns]
            M_j = M_j[:ns]
            ncs = (S.shape[1] - 1) // n_ele
            if ncs >= 2:
                S3 = S[:ns, 1 : 1 + n_ele * ncs].reshape(ns, n_ele, ncs)
                kap = S3[:, :, 1]
                if float(np.nanmean(M_i * kap)) < 0.0:
                    M_i = -M_i
                    M_j = -M_j
    groups_m: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for ip in sorted(set(ip_of.tolist())):
        ks = [k for k in range(n_ele) if ip_of[k] == ip]
        ks.sort(key=lambda k: iy_of[k])
        if not ks or not np.isfinite(y_i[ks[0]]):
            continue
        ys = [y_i[ks[0]]]
        mt = [M_i[:, ks[0]]]
        for k in ks:
            ys.append(y_j[k])
            mt.append(M_j[:, k])
        y = np.array(ys)
        vmin = np.array([m.min() for m in mt])
        vmax = np.array([m.max() for m in mt])
        groups_m[pile_name(int(ip))] = (y, vmin, vmax)
    if groups_m:
        plot_pile_depth_env(out, "pile_envelope_M.png", groups_m, "Mz min / max (kN·m)")
        print(f"PlotEQ: wrote {out / 'pile_envelope_M.png'}")

    if kap is None:
        return
    ns = kap.shape[0]
    groups_k: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for ip in sorted(set(ip_of.tolist())):
        ks = [k for k in range(n_ele) if ip_of[k] == ip and np.isfinite(y_i[k])]
        ks.sort(key=lambda k: iy_of[k])
        if not ks:
            continue
        y = y_i[ks]
        vmin = kap[:, ks].min(axis=0)
        vmax = kap[:, ks].max(axis=0)
        groups_k[pile_name(int(ip))] = (y, vmin, vmax)
    if groups_k:
        plot_pile_depth_env(
            out, "pile_envelope_kappa.png", groups_k, r"$\kappa$ min / max (1/m), i-end"
        )
        print(f"PlotEQ: wrote {out / 'pile_envelope_kappa.png'}")

    titles = []
    stations = []
    for k, (_e, ip, iy) in enumerate(rows):
        yi = float(y_i[k]) if np.isfinite(y_i[k]) else 0.0
        titles.append(f"{pile_name(ip)}  y={yi:.2f}")
        stations.append({"ip": ip, "iy": iy})
    U = np.zeros((ns, n_ele, 1))
    F = np.zeros_like(U)
    U[:, :, 0] = maybe_t0(kap)
    F[:, :, 0] = maybe_t0(M_i)
    hyst_grid(
        out, "hyst_pile_mk.png", U, F, 0, titles,
        r"$\Delta\kappa$ (1/m)" if SUBTRACT_T0 else r"$\kappa$ (1/m)",
        "Mz (kN·m)", 3,
        which=pile_hyst_order(stations, n_ele),
        share_x=True,
        share_y=True,
    )
    print(f"PlotEQ: wrote {out / 'hyst_pile_mk.png'}")


# ------------------------------------------------------------
# 4. QUAD STRESS AND STRAIN
# ------------------------------------------------------------

def read_window_quad_list(eq: Path) -> list[int]:
    """Read window quad element tags in recorder-column order.

    Args:
        eq: Serial recorder directory.
    Returns:
        Ordered quad element tags.
    """
    return [t for t, _, _ in read_window_quads(eq)]


def read_window_quads(eq: Path) -> list[tuple[int, str, str]]:
    """Read quad tags, layer names, and recorder-column groups.

    Args:
        eq: Serial recorder directory.
    Returns:
        ``(element tag, layer, column)`` rows; column may be center, ff, or window.
    """
    rows: list[tuple[int, str, str]] = []
    p = eq / "window_quads.txt"
    if not p.is_file():
        return rows
    for ln in _skip_hash(p):
        a = ln.split()
        tag = int(a[0])
        layer = a[1] if len(a) > 1 else ""
        col = a[2] if len(a) > 2 else ""
        rows.append((tag, layer, col))
    return rows


def read_ele_nodes(eq: Path) -> dict[int, list[int]]:
    """Read element connectivity keyed by element tag.

    Args:
        eq: Serial recorder directory.
    Returns:
        Element tags mapped to node-tag lists.
    """
    m: dict[int, list[int]] = {}
    p = eq / "window_eles.txt"
    if not p.is_file():
        return m
    for ln in _skip_hash(p):
        toks = [int(x) for x in ln.split()]
        m[toks[0]] = toks[1:]
    return m


def _file_ncols(path: Path) -> int:
    """Count fields in the first numeric recorder row.

    Args:
        path: Recorder-file path.
    Returns:
        Number of columns, or zero for an empty file.
    """
    with path.open() as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith("#"):
                continue
            return len(s.split())
    return 0


def peak_abs_gp_comp(
    eq: Path, files: list[str], n_ele: int, n_gp: int, n_comp: int = 3, icomp: int = 2
) -> np.ndarray:
    """Find each element's peak absolute component over time and Gauss points.

    Args:
        eq: Serial recorder directory.
        files: Recorder filenames covering all elements.
        n_ele: Total number of recorded elements.
        n_gp: Gauss points per element.
        n_comp: Components recorded at each Gauss point.
        icomp: Zero-based component index.
    Returns:
        Peak absolute value for each element in recorder units.
    """
    peak = np.zeros(n_ele)
    i0 = 0
    blk = n_gp * n_comp
    for fn in files:
        path = eq / fn
        if not path.is_file():
            raise FileNotFoundError(path)
        ncol = _file_ncols(path)
        nstr = ncol - 1
        if nstr < blk or nstr % blk != 0:
            raise SystemExit(
                f"PlotEQ: {fn} has {nstr} data cols, expected multiples of {blk}"
            )
        ne = nstr // blk
        cols = [
            1 + i * blk + g * n_comp + icomp
            for i in range(ne)
            for g in range(n_gp)
        ]
        a = np.loadtxt(path, usecols=cols)
        if a.ndim == 1:
            a = a.reshape(1, -1)
        v = np.abs(a.reshape(a.shape[0], ne, n_gp))
        peak[i0 : i0 + ne] = np.max(v, axis=(0, 2))
        i0 += ne
        del a, v
    if i0 != n_ele:
        raise SystemExit(f"PlotEQ: quad files covered {i0} eles, window_quads has {n_ele}")
    return peak


def load_quad_gp_mean_keep(
    eq: Path,
    files: list[str],
    n_ele: int,
    n_gp: int,
    keep: list[int],
    n_comp: int = 3,
    icomp: int = 2,
) -> tuple[np.ndarray, np.ndarray]:
    """Load Gauss-point means for selected global quad indices.

    Args:
        eq: Serial recorder directory.
        files: Recorder filenames covering all elements.
        n_ele: Total number of recorded elements.
        n_gp: Gauss points per element.
        keep: Global quad indices to retain.
        n_comp: Components recorded at each Gauss point.
        icomp: Zero-based component index.
    Returns:
        Time in s and selected component histories in recorder units.
    """
    pos = {int(g): j for j, g in enumerate(keep)}
    n_keep = len(keep)
    blk = n_gp * n_comp
    t: np.ndarray | None = None
    out: np.ndarray | None = None
    i0 = 0
    for fn in files:
        path = eq / fn
        if not path.is_file():
            raise FileNotFoundError(path)
        ncol = _file_ncols(path)
        nstr = ncol - 1
        if nstr < blk or nstr % blk != 0:
            raise SystemExit(
                f"PlotEQ: {fn} has {nstr} data cols, expected multiples of {blk}"
            )
        ne = nstr // blk
        local = [gi for gi in range(i0, i0 + ne) if gi in pos]
        if t is None:
            t = np.loadtxt(path, usecols=(0,))
            if t.ndim == 0:
                t = np.array([float(t)])
            out = np.zeros((len(t), n_keep))
        if local:
            cols = []
            dest = []
            for gi in local:
                i = gi - i0
                dest.append(pos[gi])
                for g in range(n_gp):
                    cols.append(1 + i * blk + g * n_comp + icomp)
            a = np.loadtxt(path, usecols=cols)
            if a.ndim == 1:
                a = a.reshape(1, -1)
            a = a.reshape(a.shape[0], len(local), n_gp).mean(axis=2)
            for k, j in enumerate(dest):
                out[:, j] = a[:, k]
            del a
        i0 += ne
    if i0 != n_ele:
        raise SystemExit(f"PlotEQ: quad files covered {i0} eles, window_quads has {n_ele}")
    if t is None or out is None:
        raise SystemExit("PlotEQ: no quad recorder files")
    return t, out


def _quad_triangulation(
    xy: dict[int, tuple[float, float]], quads: list[list[int]]
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[int, int]]:
    """Split quadrilaterals into triangles for field plotting.

    Args:
        xy: Node coordinates in m keyed by tag.
        quads: Four-node quad connectivity.
    Returns:
        Points in m, triangle indices, quad-face indices, and node point indices.
    """
    pts: list[tuple[float, float]] = []
    imap: dict[int, int] = {}

    def ix(n: int) -> int:
        """Get or create a triangulation point index.

        Args:
            n: Node tag.
        Returns:
            Point-array index.
        """
        if n not in imap:
            imap[n] = len(pts)
            pts.append(xy[n])
        return imap[n]

    tri: list[tuple[int, int, int]] = []
    face: list[int] = []
    for iq, q in enumerate(quads):
        a, b, c, d = (ix(n) for n in q)
        tri.append((a, b, c))
        tri.append((a, c, d))
        face.extend((iq, iq))
    return np.asarray(pts, float), np.asarray(tri, int), np.asarray(face, int), imap


def plot_window_peak_field(
    out: Path,
    name: str,
    xy: dict[int, tuple[float, float]],
    quads: list[list[int]],
    values: np.ndarray,
    lines: list[list[int]],
    cbar_label: str,
    title: str,
    soil_base: int = 10000,
    vmax: float | None = None,
) -> None:
    """Plot one peak quad field over the soil window.

    Args:
        out: Plot output directory.
        name: Output PNG filename.
        xy: Node coordinates in m keyed by tag.
        quads: Four-node quad connectivity.
        values: One peak value per quad in colorbar units.
        lines: Structural line connectivity.
        cbar_label: Colorbar label with units.
        title: Plot title.
        soil_base: First soil node tag.
        vmax: Colorbar maximum; ``None`` = data max.
    Returns:
        None; writes the named PNG when values are finite.
    """
    pts, triangles, face, imap = _quad_triangulation(xy, quads)
    tri = Triangulation(pts[:, 0], pts[:, 1], triangles)
    znode = np.zeros(len(pts))
    cnt = np.zeros(len(pts))
    for q, v in zip(quads, values):
        for n in q:
            i = imap[n]
            znode[i] += float(v)
            cnt[i] += 1.0
    znode /= np.maximum(cnt, 1.0)

    fig, ax = plt.subplots(figsize=(6.4, 8.0), constrained_layout=True)
    vals = np.asarray(values, dtype=float)
    face_v = vals[face]
    if not np.isfinite(face_v).any():
        print(f"PlotEQ: skip {name} (all-NaN field)")
        plt.close(fig)
        return
    face_v = np.nan_to_num(face_v, nan=0.0)
    zplot = np.nan_to_num(znode, nan=0.0)
    tpc = ax.tripcolor(
        tri,
        facecolors=face_v,
        cmap="inferno",
        edgecolors="none",
        shading="flat",
    )
    vmax_use = float(vmax) if vmax is not None else float(np.nanmax(face_v))
    if not np.isfinite(vmax_use) or vmax_use <= 0:
        vmax_use = 1.0
    tpc.set_clim(0.0, vmax_use)
    if np.isfinite(zplot).any() and float(np.nanmax(zplot) - np.nanmin(zplot)) > 0:
        ax.tricontour(tri, zplot, levels=8, colors="k", linewidths=0.35, alpha=0.45)
    segs = []
    for a, b in lines:
        if a in xy and b in xy and max(a, b) < soil_base:
            segs.append([xy[a], xy[b]])
    if segs:
        ax.add_collection(
            LineCollection(segs, colors=ORANGE, linewidths=0.9, zorder=4)
        )
    pad = 0.4
    ax.set_xlim(float(pts[:, 0].min()) - pad, float(pts[:, 0].max()) + pad)
    ax.set_ylim(float(pts[:, 1].min()) - pad, float(pts[:, 1].max()) + pad)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m), down is negative")
    ax.set_title(title, fontsize=10)
    ax.grid(True, ls=":", alpha=0.35)
    fig.colorbar(tpc, ax=ax, label=cbar_label, fraction=0.046, pad=0.04)
    fig.savefig(out / name, dpi=DPI)
    plt.close(fig)


def plot_quad_shear_peaks(
    out: Path,
    eq: Path,
    meta: dict,
    xy: dict[int, tuple[float, float]],
    lines: list[list[int]],
) -> None:
    """Plot peak quad shear stress and strain fields.

    Args:
        out: Plot output directory.
        eq: Serial recorder directory.
        meta: Window metadata describing quad recorders.
        xy: Node coordinates in m keyed by tag.
        lines: Structural line connectivity.
    Returns:
        None; writes available peak shear PNGs.
    """
    qtags = read_window_quad_list(eq)
    if not qtags:
        print("PlotEQ: no window_quads.txt -- skip shear peaks")
        return
    ev = read_ele_nodes(eq)
    quads: list[list[int]] = []
    keep: list[int] = []
    for i, t in enumerate(qtags):
        nn = ev.get(t, [])
        if len(nn) < 4 or any(n not in xy for n in nn[:4]):
            continue
        quads.append(nn[:4])
        keep.append(i)
    if len(keep) < len(qtags):
        print(f"PlotEQ: skipped {len(qtags) - len(keep)} window quads missing nodes")
    if not quads:
        print("PlotEQ: no window quad connectivity -- skip shear peaks")
        return
    n_gp = int(float(meta.get("quadNgp", 4)))
    n_ele = len(qtags)
    idx = np.asarray(keep, dtype=int)
    soil_base = node_tag_bases(meta)[0]
    sig_files = meta.get("quadStressFiles", "").split()
    eps_files = meta.get("quadStrainFiles", "").split()
    if sig_files:
        tau = peak_abs_gp_comp(eq, sig_files, n_ele, n_gp)[idx]
        plot_window_peak_field(
            out, "window_peak_tau_xy.png", xy, quads, tau / 1.0e3, lines,
            r"peak $|\tau_{xy}|$ (kPa)",
            r"peak $|\tau_{xy}|$  (max over $t$ and Gauss pts)",
            soil_base,
        )
        print(
            f"PlotEQ: wrote {out / 'window_peak_tau_xy.png'}  "
            f"max={float(tau.max()) / 1e3:.3g} kPa"
        )
    else:
        print("PlotEQ: no quadStressFiles -- skip tau contour")
    if eps_files:
        gam = peak_abs_gp_comp(eq, eps_files, n_ele, n_gp)[idx]
        plot_window_peak_field(
            out, "window_peak_gamma_xy.png", xy, quads, gam * 1.0e3, lines,
            r"peak $|\gamma_{xy}|$ ($\times 10^{-3}$)",
            r"peak $|\gamma_{xy}|$  (max over $t$ and Gauss pts)",
            soil_base,
            vmax=PEAK_GAMMA_VMAX,
        )
        print(
            f"PlotEQ: wrote {out / 'window_peak_gamma_xy.png'}  "
            f"max={float(gam.max()) * 1e3:.3g}e-3"
        )
    else:
        print("PlotEQ: no quadStrainFiles -- skip gamma contour")


def quad_centroids(
    qrows: list[tuple[int, str, str]],
    ev: dict[int, list[int]],
    xy: dict[int, tuple[float, float]],
) -> list[tuple[int, str, float, float, int, str]]:
    """Compute centroids for window quads with complete connectivity.

    Args:
        qrows: Quad tag, layer, and column rows.
        ev: Element connectivity keyed by tag.
        xy: Node coordinates in m keyed by tag.
    Returns:
        Quad tag, layer, centroid in m, global index, and column rows.
    """
    out = []
    for i, (tag, nm, col) in enumerate(qrows):
        nn = ev.get(tag, [])
        if len(nn) < 4 or any(n not in xy for n in nn[:4]):
            continue
        xc = 0.25 * sum(xy[n][0] for n in nn[:4])
        yc = 0.25 * sum(xy[n][1] for n in nn[:4])
        out.append((tag, nm, xc, yc, i, col))
    return out


def quad_depth_column(
    cents: list[tuple[int, str, float, float, int, str]],
    x_target: float | None,
    column: str | None = None,
) -> list[tuple[int, str, float, float, int, str]]:
    """Select one quad per soil row nearest a target x coordinate.

    Args:
        cents: Quad centroid rows with coordinates in m.
        x_target: Target x coordinate in m; ``None`` selects outermost absolute x.
        column: Optional recorder-column name to filter first.
    Returns:
        Selected centroid rows ordered from top to bottom.
    """
    if column:
        cents = [r for r in cents if r[5] == column]
    if not cents:
        return []
    if x_target is None:
        x_target = max(cents, key=lambda r: abs(r[2]))[2]
    by_y: dict[float, list[tuple[int, str, float, float, int, str]]] = {}
    for row in cents:
        yk = round(row[3], 6)
        by_y.setdefault(yk, []).append(row)
    col = []
    for yk in sorted(by_y.keys(), reverse=True):
        col.append(min(by_y[yk], key=lambda r: abs(r[2] - x_target)))
    return col


def plot_quad_shear_hyst(
    out: Path,
    eq: Path,
    meta: dict,
    xy: dict[int, tuple[float, float]],
) -> None:
    """Plot shear stress-strain loops through selected soil columns.

    Args:
        out: Plot output directory.
        eq: Serial recorder directory.
        meta: Window metadata describing quad recorders.
        xy: Node coordinates in m keyed by tag.
    Returns:
        None; writes available tau-gamma hysteresis PNGs.
    """
    qrows = read_window_quads(eq)
    if not qrows:
        print("PlotEQ: no window_quads.txt -- skip tau-gamma hyst")
        return
    sig_files = meta.get("quadStressFiles", "").split()
    eps_files = meta.get("quadStrainFiles", "").split()
    if not sig_files or not eps_files:
        print("PlotEQ: need quad stress and strain files for tau-gamma hyst")
        return
    ev = read_ele_nodes(eq)
    cents = quad_centroids(qrows, ev, xy)
    cols_present = sorted({r[5] for r in cents if r[5] in ("center", "ff")})
    if not cols_present:
        cols_present = [""]
    n_gp = int(float(meta.get("quadNgp", 4)))
    n_ele = len(qrows)
    for col_name in cols_present:
        if col_name == "center":
            x_tgt = 0.0
            fname = "hyst_tau_gamma_center.png"
            label = "center"
        elif col_name == "ff":
            x_tgt = float(meta["eqFFColumnX"]) if meta.get("eqFFColumnX") else None
            fname = "hyst_tau_gamma_ff.png"
            label = "near-FF"
        else:
            x_tgt = HYST_QUAD_X
            fname = "hyst_tau_gamma.png"
            label = "window"
        col = quad_depth_column(cents, x_tgt, column=col_name or None)
        if not col:
            print(f"PlotEQ: no {label} column for tau-gamma hyst")
            continue
        keep = [r[4] for r in col]
        t, tau = load_quad_gp_mean_keep(eq, sig_files, n_ele, n_gp, keep)
        t2, gam = load_quad_gp_mean_keep(eq, eps_files, n_ele, n_gp, keep)
        if len(t2) != len(t):
            n = min(len(t), len(t2))
            t, tau, gam = t[:n], tau[:n], gam[:n]
        if SUBTRACT_T0:
            tau = maybe_t0(tau)
            gam = maybe_t0(gam)
        U = np.zeros((len(t), len(col), 1))
        F = np.zeros_like(U)
        U[:, :, 0] = gam * 1.0e3
        F[:, :, 0] = tau / 1.0e3
        xc = float(np.mean([r[2] for r in col]))
        titles = [f"{nm or 'quad'}  y={yc:.2f} m" for _, nm, _, yc, _, _ in col]
        hyst_grid(
            out, fname, U, F, 0, titles,
            r"$\gamma_{xy}$ ($\times 10^{-3}$)", r"$\tau_{xy}$ (kPa)", 3,
            share_x=True,
        )
        print(
            f"PlotEQ: wrote {out / fname}  "
            f"n={len(col)}  x~{xc:.2f} m  ({label}, GP mean, d from t0={SUBTRACT_T0})"
        )
        # Keep legacy name as a copy of the center (or only) column.
        if col_name in ("", "center") and fname != "hyst_tau_gamma.png":
            shutil.copy2(out / fname, out / "hyst_tau_gamma.png")


# ------------------------------------------------------------
# 5. SSI SPRINGS AND HYSTERESIS
# ------------------------------------------------------------

def load_spring_json(meta: dict) -> dict | None:
    """Load spring capacities and station data for the active soil profile.

    Args:
        meta: Window metadata identifying the soil profile.
    Returns:
        Parsed spring JSON, or ``None`` when no file exists.
    """
    sp = meta.get("soilProfile", "")
    cands = []
    if sp:
        cands.append(pile_springs_dir(sp) / "pile_springs.json")
    cands.append(HERE / "pile_springs.json")
    for p in cands:
        if p.is_file():
            with p.open() as f:
                return json.load(f)
    return None


def cap_totals(js: dict, profile: int) -> tuple[float, float, float, float, float]:
    """Get or calculate cap spring capacities and reference deformations.

    Args:
        js: Spring input data in N and m.
        profile: Soil-profile number.
    Returns:
        Lateral, shaft, and soffit capacities in N; y50 and z50 in m.
    """
    H = float(js.get("H_cap", 0.9906))
    y50 = float(js.get("y50_cap", 0.01))
    z50 = float(js.get("z50_cap", 0.01))
    if "PultCap" in js and "TultCap" in js and "QultSoffit" in js:
        return float(js["PultCap"]), float(js["TultCap"]), float(js["QultSoffit"]), y50, z50
    if profile == 4:
        cu, rho = 1.80e4, 1300.0
    else:
        cu, rho = 3.59e4, 1488.0
    foot = 0.3048
    W = float(js.get("W_cap", 15.0 * foot))
    L = float(js.get("L_cap", 10.0 * foot))
    g = 9.81
    rho_w = 1000.0
    alpha = 0.75
    gam = (rho - rho_w) * g
    P = cu * L * H / 2.0 * (4.0 + gam * H / cu + 0.25 * H / L + 2.0 * alpha)
    T = cu * 2.0 * H * (W + L)
    Q = 9.0 * cu * W * L
    return P, T, Q, y50, z50


def load_sketch_sizes(meta: dict) -> dict:
    """Load model-sketch dimensions for spring plot fallbacks.

    Args:
        meta: Window metadata identifying soil profile and boundary.
    Returns:
        Model size values in m, or an empty dictionary.
    """
    sp = meta.get("soilProfile", "")
    bnd = meta.get("soilBoundary", "Shin")
    cands = []
    if sp:
        cands.append(elevation_dir(sp, bnd) / "model_sketch.json")
    cands.append(HERE / "model_sketch.json")
    for pth in cands:
        if pth.is_file():
            with pth.open() as f:
                return json.load(f).get("sizes") or {}
    return {}


def trib_from_x(xs: np.ndarray) -> np.ndarray:
    """Compute one-dimensional tributary widths.

    Args:
        xs: Ordered x coordinates in m.
    Returns:
        Tributary widths in m.
    """
    n = len(xs)
    trib = np.ones(n)
    if n == 1:
        return trib
    trib[0] = 0.5 * (xs[1] - xs[0])
    trib[-1] = 0.5 * (xs[-1] - xs[-2])
    for i in range(1, n - 1):
        trib[i] = 0.5 * (xs[i + 1] - xs[i - 1])
    return trib


def soffit_x_default(nsof: int, sizes: dict) -> np.ndarray:
    """Build fallback x coordinates for cap-soffit springs.

    Args:
        nsof: Number of soffit springs.
        sizes: Model dimensions in m.
    Returns:
        Soffit spring x coordinates in m.
    """
    s = float(sizes.get("s_pile_cap", 1.8288))
    wsoil = float(sizes.get("W_cap_soil", 2.0 * s))
    xf = 0.5 * wsoil
    xs = np.array([-s, -0.5 * s, 0.0, 0.5 * s, s])
    if len(xs) != nsof:
        return np.linspace(-xf, xf, nsof)
    return xs


def split_cap_eles(eq: Path) -> tuple[list[int], list[int]]:
    """Separate cap-face and cap-soffit spring element tags.

    Args:
        eq: Serial recorder directory.
    Returns:
        Cap-face tags and cap-soffit tags.
    """
    face_e, sof_e = [], []
    path = eq / "cap_springs_eles.txt"
    if not path.is_file():
        return face_e, sof_e
    for ln in _skip_hash(path):
        parts = ln.split()
        e = int(parts[0])
        k = parts[1] if len(parts) > 1 else "cap"
        if k == "cap_soffit":
            sof_e.append(e)
        else:
            face_e.append(e)
    return face_e, sof_e


def read_ele_list(path: Path) -> list[int]:
    """Read the first element tag from each data line.

    Args:
        path: Element-list file path.
    Returns:
        Ordered element tags.
    """
    out = []
    for ln in _skip_hash(path):
        out.append(int(ln.split()[0]))
    return out


def read_pile_spring_eles(eq: Path) -> list[dict]:
    """Read pile spring identity and station indices.

    Args:
        eq: Serial recorder directory.
    Returns:
        Rows with element tag, kind, and optional pile/station indices.
    """
    path = eq / "pile_springs_eles.txt"
    out: list[dict] = []
    for ln in _skip_hash(path):
        parts = ln.split()
        row: dict = {"e": int(parts[0]), "kind": parts[1] if len(parts) > 1 else "pile"}
        if len(parts) >= 4:
            row["ip"] = int(parts[2])
            row["iy"] = int(parts[3])
        out.append(row)
    return out


def parse_lean_stations(meta: dict) -> list[dict]:
    """Parse lean-dump SSI stations from window metadata.

    Args:
        meta: Window metadata containing ``leanStations``.
    Returns:
        Station rows with index, elevation in m, layer, tip flag, and segment.
    """
    raw = meta.get("leanStations", "").strip()
    if not raw:
        return []
    out: list[dict] = []
    for inner in re.findall(r"\{([^}]*)\}", raw):
        a = inner.split()
        if len(a) < 5:
            continue
        out.append({
            "iy": int(float(a[0])),
            "y": float(a[1]),
            "layer": a[2],
            "isTip": int(float(a[3])),
            "iSeg": int(float(a[4])),
        })
    return out


def stations_for_plot(
    rows: list[dict], js: dict | None, meta: dict
) -> list[dict] | None:
    """Align spring stations with recorder columns.

    Args:
        rows: Recorder element rows.
        js: Spring JSON, or ``None``.
        meta: Window metadata with lean-station fallback data.
    Returns:
        Stations in recorder-column order, or ``None`` when alignment fails.
    """
    aligned = align_pile_stations(rows, (js or {}).get("stations") or [])
    if aligned is not None:
        return aligned
    lean = parse_lean_stations(meta)
    if not lean or not rows:
        return None
    by_iy = {int(s["iy"]): s for s in lean}
    out: list[dict] = []
    for r in rows:
        if "iy" not in r or "ip" not in r:
            return None
        s = by_iy.get(int(r["iy"]))
        if s is None:
            return None
        out.append({
            "ip": int(r["ip"]),
            "iy": int(r["iy"]),
            "y": float(s["y"]),
            "isTip": int(s["isTip"]),
            "layer": s.get("layer", ""),
        })
    return out


def align_pile_stations(
    rows: list[dict], stations: list[dict]
) -> list[dict] | None:
    """Match station data to recorder rows by pile and vertical index.

    Args:
        rows: Recorder element rows.
        stations: Candidate station dictionaries.
    Returns:
        Stations in recorder-column order, or ``None`` when unmatched.
    """
    if not rows:
        return []
    if not stations:
        return None
    by = {(int(s["ip"]), int(s["iy"])): s for s in stations}
    if rows and all("ip" in r and int(r["ip"]) >= 0 for r in rows):
        out = []
        for r in rows:
            s = by.get((int(r["ip"]), int(r["iy"])))
            if s is None:
                return None
            out.append(s)
        return out
    if len(stations) == len(rows):
        return list(stations)
    return None


def env_minmax(arr: np.ndarray, icomp: int) -> tuple[np.ndarray, np.ndarray]:
    """Find time-history extrema for one recorder component.

    Args:
        arr: Array shaped ``(time, element, component)`` in native units.
        icomp: Zero-based component index.
    Returns:
        Minimum and maximum arrays in the input units.
    """
    v = arr[:, :, icomp]
    return v.min(axis=0), v.max(axis=0)


def pile_name(ip: int) -> str:
    """Convert a zero-based pile index to a plot label.

    Args:
        ip: Zero-based pile index.
    Returns:
        ``L``, ``C``, ``R``, or a generated pile label.
    """
    return NAMES[ip] if ip < len(NAMES) else f"p{ip}"


def is_qz_station(s: dict) -> bool:
    """Test whether a station represents a pile-tip q-z spring.

    Args:
        s: Spring station dictionary.
    Returns:
        ``True`` for a tip/q-z station.
    """
    return bool(s.get("isTip")) or s.get("axType") == "qz"


def pile_hyst_order(stations: list, n_ele: int, qz: bool | None = None) -> list[int]:
    """Order pile hysteresis panels by depth and pile.

    Args:
        stations: Spring or beam station dictionaries.
        n_ele: Number of recorded elements.
        qz: ``True`` for tips, ``False`` for shafts, or ``None`` for all.
    Returns:
        Recorder-column indices ordered by station then pile.
    """
    if len(stations) != n_ele:
        return list(range(n_ele))
    by_iy: dict[int, list[tuple[int, int]]] = {}
    for i, s in enumerate(stations):
        is_q = is_qz_station(s)
        if qz is True and not is_q:
            continue
        if qz is False and is_q:
            continue
        iy = int(s["iy"])
        by_iy.setdefault(iy, []).append((int(s["ip"]), i))
    out = []
    for iy in sorted(by_iy):
        recs = sorted(by_iy[iy], key=lambda t: t[0])
        out.extend(i for _, i in recs)
    return out or list(range(n_ele))


def plot_depth_env(
    out: Path,
    name: str,
    y_cap: np.ndarray,
    cap_pos: np.ndarray,
    xlabel: str,
    cap_label: str,
    groups: list[tuple[str, np.ndarray, np.ndarray, np.ndarray]],
    cap_neg: np.ndarray | None = None,
) -> None:
    """Plot spring extrema and capacities against depth.

    Args:
        out: Plot output directory.
        name: Output PNG filename.
        y_cap: Capacity elevations in m.
        cap_pos: Positive capacities/reference deformations in axis units.
        xlabel: Horizontal-axis label with units.
        cap_label: Capacity curve label.
        groups: Name, elevation in m, minimum, and maximum arrays.
        cap_neg: Optional negative capacity in axis units.
    Returns:
        None; writes the named PNG.
    """
    cap_pos = np.asarray(cap_pos, dtype=float)
    if cap_neg is None:
        cap_neg = -cap_pos
    else:
        cap_neg = np.asarray(cap_neg, dtype=float)
    fig, ax = plt.subplots(figsize=(5.8, 7.4), constrained_layout=True)
    xs: list[np.ndarray] = [cap_pos, cap_neg]
    ax.plot(cap_neg, y_cap, color=GRAY, lw=1.2, label=cap_label)
    ax.plot(cap_pos, y_cap, color=GRAY, lw=1.2)
    for lab, yi, vn, vp in groups:
        c = COLORS.get(lab, "#333")
        ax.plot(vn, yi, color=c, lw=1.5, label=f"{lab} min")
        ax.plot(vp, yi, color=c, lw=1.5, ls="--", label=f"{lab} max")
        xs.extend([vn, vp])
    ax.axvline(0.0, color="#bbb", lw=0.8)
    sym_xlim(ax, *xs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("y (m), down is negative")
    ax.legend(fontsize=7, loc="best")
    ax.grid(True, ls=":", alpha=0.45)
    fig.savefig(out / name, dpi=DPI)
    plt.close(fig)


def plot_x_env(
    out: Path,
    name: str,
    x: np.ndarray,
    vmin: np.ndarray,
    vmax: np.ndarray,
    cap_pos: np.ndarray,
    ylabel: str,
    cap_label: str,
    compression_only: bool = False,
) -> None:
    """Plot spring extrema and capacities along x.

    Args:
        out: Plot output directory.
        name: Output PNG filename.
        x: Spring x coordinates in m.
        vmin: Minimum response in axis units.
        vmax: Maximum response in axis units.
        cap_pos: Capacity/reference values in axis units.
        ylabel: Vertical-axis label with units.
        cap_label: Capacity curve label.
        compression_only: Show only the negative capacity branch when true.
    Returns:
        None; writes the named PNG.
    """
    cap_pos = np.asarray(cap_pos, dtype=float)
    fig, ax = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
    ax.plot(x, vmin, color=BROWN, lw=1.5, marker="o", ms=4, label="min")
    ax.plot(x, vmax, color=BROWN, lw=1.5, ls="--", marker="o", ms=4, label="max")
    if compression_only:
        cap_c = -np.abs(cap_pos)
        ax.plot(x, cap_c, color=GRAY, lw=1.2, label=cap_label)
        ylim_vals = (vmin, vmax, cap_c)
    else:
        ax.plot(x, cap_pos, color=GRAY, lw=1.2, label=cap_label)
        ax.plot(x, -cap_pos, color=GRAY, lw=1.2)
        ylim_vals = (vmin, vmax, cap_pos)
    ax.axhline(0.0, color="#bbb", lw=0.8)
    m = 0.0
    for v in ylim_vals:
        a = np.asarray(v, dtype=float)
        if a.size and np.isfinite(a).any():
            m = max(m, float(np.nanmax(np.abs(a))))
    if not np.isfinite(m) or m <= 0.0:
        m = 1.0
    ax.set_ylim(-1.05 * m, 1.05 * m)
    ax.set_xlabel("x (m)")
    ax.set_ylabel(ylabel)
    ax.legend(fontsize=8)
    ax.grid(True, ls=":", alpha=0.45)
    fig.savefig(out / name, dpi=DPI)
    plt.close(fig)


def plot_spring_envelopes(
    out: Path,
    eq: Path,
    js: dict | None,
    meta: dict,
    groups_xy: dict[str, list[int]],
    xy: dict[int, tuple[float, float]],
) -> None:
    """Plot pile and cap spring deformation/force envelopes.

    Args:
        out: Plot output directory.
        eq: Serial recorder directory.
        js: Spring capacities and stations, or ``None``.
        meta: Window metadata.
        groups_xy: Pile names mapped to ordered node tags.
        xy: Node coordinates in m keyed by tag.
    Returns:
        None; writes available spring envelope PNGs.
    """
    rows = read_pile_spring_eles(eq)
    pile_eles = [r["e"] for r in rows]
    _, Fp, Up = load_spring_pt(
        eq, "pile_springs_force.out", "pile_springs_defo.out", len(pile_eles)
    )
    n = len(pile_eles)
    stations = stations_for_plot(rows, js, meta)
    if stations is None or len(stations) != n:
        print(
            f"PlotEQ: pile stations {0 if stations is None else len(stations)}"
            f" != nEle {n} -- skip spring env"
        )
        return
    if n and ("pult" not in stations[0] or "y50" not in stations[0]):
        print("PlotEQ: no pult/y50 (need pile_springs.json) -- skip spring env")
        return

    ips_idx: dict[int, list[int]] = {}
    for i, s in enumerate(stations):
        ips_idx.setdefault(int(s["ip"]), []).append(i)
    g = []
    for ip in sorted(ips_idx):
        idx = np.array(ips_idx[ip], dtype=int)
        # MP stitch order follows rank, not depth (lean: iy 13,14,17,0,3,...).
        order = np.argsort([int(stations[i]["iy"]) for i in idx])
        idx = idx[order]
        y = np.array([float(stations[i]["y"]) for i in idx])
        g.append((pile_name(ip), idx, y))
    y0 = g[0][2]
    i0 = g[0][1]
    y50 = np.array([float(stations[i]["y50"]) for i in i0])
    z50 = np.array([float(stations[i]["z50"]) for i in i0])
    pult = np.array([float(stations[i]["pult"]) for i in i0])
    tult = np.array([float(stations[i]["tult"]) for i in i0])

    def as_groups_u(minmax):
        """Spring deflection extrema (m → mm) by pile group."""
        umin, umax = minmax
        umin, umax = to_mm(umin), to_mm(umax)
        return [(name, y, umin[idx], umax[idx]) for name, idx, y in g]

    def as_groups(minmax):
        """Force extrema by pile group (native N)."""
        umin, umax = minmax
        return [(name, y, umin[idx], umax[idx]) for name, idx, y in g]

    u_py = env_minmax(Up, 0)
    f_py = env_minmax(Fp, 0)
    plot_depth_env(
        out, "spring_env_py_defo.png", y0, y50,
        "u_py min / max (mm)", r"$\pm y_{50}$", as_groups_u(u_py),
    )
    plot_depth_env(
        out, "spring_env_py_force.png", y0, pult,
        "p min / max (N)", r"$\pm p_\mathrm{ult}$", as_groups(f_py),
    )
    if Up.shape[2] > 1:
        u_tz = env_minmax(Up, 1)
        f_tz = env_minmax(Fp, 1)
        shaft_g = []
        for name, idx, y in g:
            keep = [j for j, i in enumerate(idx) if not is_qz_station(stations[i])]
            if not keep:
                continue
            idx_s = idx[keep]
            y_s = y[keep]
            shaft_g.append((name, idx_s, y_s))
        if shaft_g:
            y_s0 = shaft_g[0][2]
            i_s0 = shaft_g[0][1]
            z50_s = np.array([float(stations[i]["z50"]) for i in i_s0])
            tult_s = np.array([float(stations[i]["tult"]) for i in i_s0])

            def as_shaft_u(minmax):
                umin, umax = minmax
                umin, umax = to_mm(umin), to_mm(umax)
                return [(nm, yi, umin[idx], umax[idx]) for nm, idx, yi in shaft_g]

            def as_shaft(minmax):
                umin, umax = minmax
                return [(nm, yi, umin[idx], umax[idx]) for nm, idx, yi in shaft_g]

            plot_depth_env(
                out, "spring_env_tz_defo.png", y_s0, z50_s,
                "u_tz min / max (mm)", r"$\pm z_{50}$ (shaft t-z)",
                as_shaft_u(u_tz),
            )
            plot_depth_env(
                out, "spring_env_tz_force.png", y_s0, tult_s,
                "t min / max (N)", r"$\pm t_\mathrm{ult}$ (shaft t-z)",
                as_shaft(f_tz),
            )
        tip = [i for i, s in enumerate(stations) if is_qz_station(s)]
        if tip:
            xp = []
            for i in tip:
                ip = int(stations[i]["ip"])
                nm = pile_name(ip)
                tags = groups_xy.get(nm, [])
                xp.append(xy[tags[0]][0] if tags else float(ip))
            xp = np.array(xp)
            z50q = np.array([float(stations[i]["z50"]) for i in tip])
            qult = np.array([float(stations[i]["tult"]) for i in tip])
            plot_x_env(
                out, "spring_env_pile_qz_defo.png", xp,
                to_mm(u_tz[0][tip]), to_mm(u_tz[1][tip]), z50q,
                "u_qz min / max (mm)", r"$z_{50}$ compression",
                compression_only=True,
            )
            plot_x_env(
                out, "spring_env_pile_qz_force.png", xp,
                f_tz[0][tip], f_tz[1][tip], qult,
                "q min / max (N)", r"$q_\mathrm{ult}$ compression",
                compression_only=True,
            )

    face_e, sof_e = split_cap_eles(eq)
    P, T, Q, y50c, z50c = cap_totals(js or {}, int(meta.get("soilProfile", 3)))
    H = float((js or {}).get("H_cap", 0.9906))
    sizes = load_sketch_sizes(meta)

    if face_e and (eq / "cap_springs_force.out").is_file():
        _, Fc, Uc = load_spring_pt(
            eq, "cap_springs_force.out", "cap_springs_defo.out", len(face_e)
        )
        face = (js or {}).get("cap_face") or []
        hfrac = np.array([0.25, 0.5, 0.25, 0.25, 0.5, 0.25])[: len(face_e)]
        if len(face) == len(face_e):
            yf = np.array([float(r["y"]) for r in face])
            pcap = np.array([float(r["pult"]) for r in face])
            tcap = np.array([float(r["tult"]) for r in face])
            y50c_v = np.array([float(r["y50"]) for r in face])
            z50c_v = np.array([float(r["z50"]) for r in face])
        else:
            yf = np.array([0.0, -0.5 * H, -H, 0.0, -0.5 * H, -H])[: len(face_e)]
            pcap = P * hfrac / 2.0
            tcap = T * hfrac / 2.0
            y50c_v = np.full(len(face_e), y50c)
            z50c_v = np.full(len(face_e), z50c)
        n2 = len(face_e) // 2
        yL = yf[:n2] if n2 else yf

        def cap_groups_u(minmax, yf=yf, n2=n2):
            umin, umax = minmax
            umin, umax = to_mm(umin), to_mm(umax)
            if not n2:
                return [("cap", yf, umin, umax)]
            return [
                ("L", yf[:n2], umin[:n2], umax[:n2]),
                ("R", yf[n2:], umin[n2:], umax[n2:]),
            ]

        def cap_groups(minmax, yf=yf, n2=n2):
            """Arrange cap-face extrema into left/right depth groups.

            Args:
                minmax: Minimum and maximum arrays in native units.
                yf: Cap-face elevations in m.
                n2: Number of stations on one cap face.
            Returns:
                Cap-face labels, elevations in m, and grouped extrema.
            """
            umin, umax = minmax
            if not n2:
                return [("cap", yf, umin, umax)]
            return [
                ("L", yf[:n2], umin[:n2], umax[:n2]),
                ("R", yf[n2:], umin[n2:], umax[n2:]),
            ]

        plot_depth_env(
            out, "spring_env_cap_py_defo.png", yL, y50c_v[:n2] if n2 else y50c_v,
            "u_py min / max (mm)", r"$\pm y_{50}$ cap", cap_groups_u(env_minmax(Uc, 0)),
        )
        plot_depth_env(
            out, "spring_env_cap_py_force.png", yL, pcap[:n2] if n2 else pcap,
            "p min / max (N)", r"$\pm p_\mathrm{ult}$ cap", cap_groups(env_minmax(Fc, 0)),
        )
        if Uc.shape[2] > 1:
            plot_depth_env(
                out, "spring_env_cap_tz_defo.png", yL, z50c_v[:n2] if n2 else z50c_v,
                "u_tz min / max (mm)", r"$\pm z_{50}$ cap", cap_groups_u(env_minmax(Uc, 1)),
            )
            plot_depth_env(
                out, "spring_env_cap_tz_force.png", yL, tcap[:n2] if n2 else tcap,
                "t min / max (N)", r"$\pm t_\mathrm{ult}$ cap", cap_groups(env_minmax(Fc, 1)),
            )

    if sof_e and (eq / "cap_springs_soffit_force.out").is_file():
        nsof = len(sof_e)
        _, Fq, Uq = load_spring_pt(
            eq, "cap_springs_soffit_force.out", "cap_springs_soffit_defo.out", nsof
        )
        ic = 0 if Fq.shape[2] == 1 else 1
        umin, umax = env_minmax(Uq, ic)
        fmin, fmax = env_minmax(Fq, ic)
        soff = (js or {}).get("cap_soffit") or []
        if len(soff) == nsof:
            xq = np.array([float(r["x"]) for r in soff])
            qult = np.array([float(r["qult"]) for r in soff])
            z50q = np.array([float(r.get("z50", z50c)) for r in soff])
        else:
            xq = soffit_x_default(nsof, sizes)
            trib = trib_from_x(xq)
            qult = Q * trib / trib.sum()
            z50q = np.full(nsof, z50c)
        plot_x_env(
            out, "spring_env_qz_defo.png", xq, to_mm(umin), to_mm(umax), z50q,
            "u_qz min / max (mm)", r"$z_{50}$ compression",
            compression_only=True,
        )
        plot_x_env(
            out, "spring_env_qz_force.png", xq, fmin, fmax, qult,
            "q min / max (N)", r"$q_\mathrm{ult}$ compression",
            compression_only=True,
        )


def hyst_grid(
    out: Path,
    fname: str,
    U: np.ndarray,
    F: np.ndarray,
    icomp: int,
    titles: list[str],
    xlabel: str,
    ylabel: str,
    ncol: int,
    which: list[int] | None = None,
    share_x: bool = False,
    share_y: bool = False,
) -> None:
    """Plot many element hysteresis loops on a panel grid.

    Args:
        out: Plot output directory.
        fname: Output PNG filename.
        U: Deformation histories in label-defined units.
        F: Force histories in label-defined units.
        icomp: Zero-based component index.
        titles: Element panel titles.
        xlabel: Horizontal-axis label with units.
        ylabel: Vertical-axis label with units.
        ncol: Number of panel columns.
        which: Optional recorder-column order/subset.
        share_x: Use one symmetric x range.
        share_y: Use one symmetric y range.
    Returns:
        None; writes the named PNG.
    """
    if which is None:
        which = list(range(U.shape[1]))
    n = len(which)
    if n == 0:
        return
    xlim = ylim = None
    if share_x:
        xm = max(float(np.nanmax(np.abs(U[:, i, icomp]))) for i in which)
        xlim = (-1.05 * (xm if xm > 0 else 1.0), 1.05 * (xm if xm > 0 else 1.0))
    if share_y:
        ym = max(float(np.nanmax(np.abs(F[:, i, icomp]))) for i in which)
        ylim = (-1.05 * (ym if ym > 0 else 1.0), 1.05 * (ym if ym > 0 else 1.0))
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(
        nrow, ncol, figsize=(3.0 * ncol, 1.55 * nrow), squeeze=False
    )
    for k, i in enumerate(which):
        r, c = divmod(k, ncol)
        ax = axes[r][c]
        ax.plot(U[:, i, icomp], F[:, i, icomp], color=BROWN, lw=0.55, rasterized=True)
        ax.set_title(titles[i] if i < len(titles) else f"ele {i}", fontsize=8)
        ax.grid(True, ls=":", alpha=0.35)
        if xlim is not None:
            ax.set_xlim(*xlim)
        if ylim is not None:
            ax.set_ylim(*ylim)
        if r == nrow - 1:
            ax.set_xlabel(xlabel, fontsize=8)
        if c == 0:
            ax.set_ylabel(ylabel, fontsize=8)
        ax.tick_params(labelsize=7)
    for k in range(n, nrow * ncol):
        r, c = divmod(k, ncol)
        axes[r][c].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / fname, dpi=110)
    plt.close(fig)


def plot_hyst(out: Path, eq: Path, js: dict | None, meta: dict) -> None:
    """Plot pile and cap spring hysteresis grids.

    Args:
        out: Plot output directory.
        eq: Serial recorder directory.
        js: Spring capacities and stations, or ``None``.
        meta: Window metadata.
    Returns:
        None; writes available spring hysteresis PNGs.
    """
    rows = read_pile_spring_eles(eq)
    pile_eles = [r["e"] for r in rows]
    _, Fp, Up = load_spring_pt(
        eq, "pile_springs_force.out", "pile_springs_defo.out", len(pile_eles)
    )
    stations = stations_for_plot(rows, js, meta) or []
    sizes = load_sketch_sizes(meta)
    titles = []
    for i, e in enumerate(pile_eles):
        if i < len(stations):
            s = stations[i]
            titles.append(f"{pile_name(int(s['ip']))}  y={float(s['y']):.2f}")
        else:
            titles.append(f"ele {e}")
    Up_mm = Up * M_TO_MM
    hyst_grid(
        out, "hyst_pile_py.png", Up_mm, Fp, 0, titles, "u_py (mm)", "p (N)", 3,
        which=pile_hyst_order(stations, Up.shape[1]),
        share_x=True,
    )
    if Up.shape[2] > 1:
        hyst_grid(
            out, "hyst_pile_tz.png", Up_mm, Fp, 1, titles,
            "u_tz (mm)", "t (N)", 3, which=pile_hyst_order(stations, Up.shape[1], qz=False),
            share_x=True,
        )
        qz_i = pile_hyst_order(stations, Up.shape[1], qz=True)
        if qz_i:
            hyst_grid(
                out, "hyst_pile_qz.png", Up_mm, Fp, 1, titles,
                "u_qz (mm)", "q (N)", 3, which=qz_i,
                share_x=True,
            )

    face_e, sof_e = split_cap_eles(eq)
    face_labs = ["L top", "L mid", "L bot", "R top", "R mid", "R bot"]
    if face_e and (eq / "cap_springs_force.out").is_file():
        _, Fc, Uc = load_spring_pt(
            eq, "cap_springs_force.out", "cap_springs_defo.out", len(face_e)
        )
        face = (js or {}).get("cap_face") or []
        tcap = []
        for i, e in enumerate(face_e):
            if i < len(face):
                tcap.append(f"cap y={float(face[i]['y']):.2f}")
            elif i < len(face_labs):
                tcap.append(face_labs[i])
            else:
                tcap.append(f"cap ele {e}")
        cap_ord = [0, 3, 1, 4, 2, 5] if len(face_e) == 6 else None
        Uc_mm = Uc * M_TO_MM
        hyst_grid(
            out, "hyst_cap_py.png", Uc_mm, Fc, 0, tcap, "u_py (mm)", "p (N)", 2,
            which=cap_ord,
            share_x=True,
        )
        if Uc.shape[2] > 1:
            hyst_grid(
                out, "hyst_cap_tz.png", Uc_mm, Fc, 1, tcap, "u_tz (mm)", "t (N)", 2,
                which=cap_ord,
                share_x=True,
            )
    if sof_e and (eq / "cap_springs_soffit_force.out").is_file():
        _, Fq, Uq = load_spring_pt(
            eq, "cap_springs_soffit_force.out",
            "cap_springs_soffit_defo.out", len(sof_e),
        )
        ic = 0 if Fq.shape[2] == 1 else 1
        soff = (js or {}).get("cap_soffit") or []
        xq = (
            np.array([float(r["x"]) for r in soff])
            if len(soff) == len(sof_e)
            else soffit_x_default(len(sof_e), sizes)
        )
        tsof = [f"qz x={x:.2f}" for x in xq]
        Uq_mm = Uq * M_TO_MM
        hyst_grid(
            out, "hyst_soffit_qz.png", Uq_mm, Fq, ic, tsof, "u_qz (mm)", "q (N)", 4,
            share_x=True,
        )


# ------------------------------------------------------------
# 6. DEFORMED FRAMES AND MOVIES
# ------------------------------------------------------------


def _nearest_frame_depth(
    column: list[tuple[int, float]],
    targets: tuple[float, ...],
    *,
    y_grade: float,
) -> list[tuple[int, float, float]]:
    """Pick nodes nearest each target depth; always include tip.

    Args:
        column: ``(tag, y)`` head-to-tip.
        targets: Depths below grade in m.
        y_grade: Grade elevation in m.
    Returns:
        ``(tag, y, depth)`` unique, shallow to deep.
    """
    if not column:
        return []
    picked: list[tuple[int, float, float]] = []
    used: set[int] = set()

    def add_nearest(d_tgt: float) -> None:
        best = None
        best_err = 1e99
        for tg, y in column:
            if tg in used:
                continue
            d = abs(float(y) - y_grade)
            err = abs(d - d_tgt)
            if err < best_err:
                best_err = err
                best = (tg, float(y), d)
        if best is not None:
            used.add(best[0])
            picked.append(best)

    for d_tgt in targets:
        add_nearest(d_tgt)
    tip_tg, tip_y = column[-1]
    if tip_tg not in used:
        picked.append((tip_tg, float(tip_y), abs(float(tip_y) - y_grade)))
    picked.sort(key=lambda r: r[2])
    return picked


def pick_frame_traces(
    tags: list[int],
    xy: dict[int, tuple[float, float]],
    meta: dict,
    idx: dict[int, int],
) -> list[tuple[str, int, str]]:
    """Pier top/base + pile and soil ux stations for frame side panels.

    Args:
        tags, xy, meta, idx: Window maps.
    Returns:
        ``(label, tag, color)`` in plot order.
    """
    cmap = plt.get_cmap("tab10")
    out: list[tuple[str, int, str]] = []
    ci = 0
    has_pier_base = False

    def push(lab: str, tg: int) -> None:
        nonlocal ci
        if tg not in idx or tg not in xy:
            return
        out.append((lab, tg, cmap(ci % 10)))
        ci += 1

    if PIER_TOP in idx:
        push("pier top", PIER_TOP)
    if PIER_BOT in idx:
        push("pier base", PIER_BOT)
        has_pier_base = True

    groups = pile_groups({tg: xy[tg] for tg in tags if tg in xy})
    _name, pile_tags = preferred_pile_tags(groups)
    pile_col = [(tg, xy[tg][1]) for tg in pile_tags if tg in idx and tg in xy]
    y_grade = float(pile_col[0][1]) if pile_col else 0.0
    if 1 in xy:
        y_grade = max(y_grade, float(xy[1][1]))

    pile_targets = FRAME_HIST_DEPTH_M[1:] if has_pier_base else FRAME_HIST_DEPTH_M
    for tg, y, d in _nearest_frame_depth(pile_col, pile_targets, y_grade=y_grade):
        is_tip = pile_col and tg == pile_col[-1][0]
        if is_tip:
            push(f"pile base ({format_depth_label(y)} m)", tg)
        elif abs(d) < 0.15 and not has_pier_base:
            push("pile 0.0 m", tg)
        else:
            tgt = min(FRAME_HIST_DEPTH_M, key=lambda t: abs(t - d))
            if abs(tgt - d) < 1.0:
                push(f"pile {tgt:.1f} m", tg)
            else:
                push(f"pile {format_depth_label(y)} m", tg)

    soil_col = soil_column_nodes(tags, xy, meta, x_tgt=0.0)
    for tg, y, d in _nearest_frame_depth(soil_col, FRAME_HIST_DEPTH_M, y_grade=y_grade):
        is_tip = soil_col and tg == soil_col[-1][0]
        if is_tip:
            push(f"soil base ({format_depth_label(y)} m)", tg)
        elif abs(d) < 0.15:
            push("soil 0.0 m", tg)
        else:
            tgt = min(FRAME_HIST_DEPTH_M, key=lambda t: abs(t - d))
            if abs(tgt - d) < 1.0:
                push(f"soil {tgt:.1f} m", tg)
            else:
                push(f"soil {format_depth_label(y)} m", tg)
    return out


def _frame_mesh_parts(
    tags: list[int],
    xy: dict[int, tuple[float, float]],
    lines: list[list[int]],
    quads: list[list[int]],
    js: dict | None,
    meta: dict | None,
) -> dict:
    """Soil patches, structure, and spring segments for deformed frames.

    Returns:
        Dict with idx, X0, Y0, q_idx, ij_s, ij_z, fb, face, names_used, prof,
        xlim, ylim (limits filled by caller from amp).
    """
    idx = {tg: i for i, tg in enumerate(tags)}
    X0 = np.array([xy[tg][0] for tg in tags])
    Y0 = np.array([xy[tg][1] for tg in tags])
    soil_base, spr_base, soffit_off, _bnd, soil_last = node_tag_bases(meta)

    soil_quads = []
    for q in quads:
        if len(q) < 4:
            continue
        nn = q[:4]
        if (
            min(nn) >= soil_base
            and max(nn) <= soil_last
            and all(n in idx for n in nn)
        ):
            soil_quads.append(nn)
    q_idx = (
        np.array([[idx[n] for n in q] for q in soil_quads], dtype=int)
        if soil_quads
        else np.zeros((0, 4), dtype=int)
    )

    layers = list((js or {}).get("soil_layers") or (js or {}).get("layers") or [])
    prof = (js or {}).get("soilProfile")
    try:
        prof = int(prof) if prof is not None else None
    except (TypeError, ValueError):
        prof = None
    names_used: list[str] = []
    face = []
    for q in soil_quads:
        yc = 0.25 * sum(xy[n][1] for n in q)
        nm = layer_at_y(yc, layers)
        st = layer_style(nm, profile=prof)
        if nm not in names_used:
            names_used.append(nm)
        face.append(to_rgba(st["fill"], alpha=0.55))

    struct, spr = [], []
    for ln in lines:
        if len(ln) < 2:
            continue
        a, b = ln[0], ln[1]
        hi, lo = max(a, b), min(a, b)
        if hi >= spr_base:
            spr.append((lo, hi) if lo < spr_base else (a, b))
        elif hi < soil_base:
            struct.append((a, b))
    ij_s = _pair_idx(struct, idx)
    ij_z = _pair_idx(spr, idx)
    fb = np.zeros((len(ij_z), 2))
    fb[:, 0] = 1.0
    spr_found = min((n for n in tags if n >= spr_base), default=spr_base)
    for kk, (a, b) in enumerate(spr):
        dup = b if b >= spr_base else a
        if dup - spr_found >= soffit_off:
            fb[kk] = (0.0, 1.0)

    return {
        "idx": idx,
        "X0": X0,
        "Y0": Y0,
        "q_idx": q_idx,
        "ij_s": ij_s,
        "ij_z": ij_z,
        "fb": fb,
        "face": face,
        "names_used": names_used,
        "prof": prof,
    }


def to_mm(u) -> np.ndarray:
    """Recorder displacement (m) → plot units (mm)."""
    return np.asarray(u, dtype=float) * M_TO_MM


def t_proto_to_model(t_proto: float | np.ndarray) -> np.ndarray:
    """Prototype recorder time → Froude model (lab) time (s)."""
    return np.asarray(t_proto, dtype=float) / TIME_SCALE_FROUDE


def t_model_to_proto(t_model: float | np.ndarray) -> np.ndarray:
    """Froude model (lab) time → prototype recorder time (s)."""
    return np.asarray(t_model, dtype=float) * TIME_SCALE_FROUDE


def frame_t0_proto() -> float:
    """First prototype time included in animation frames (s)."""
    return FRAME_T0_MODEL_S * TIME_SCALE_FROUDE


def frame_time_title(t_proto: float, ifr: int, nfr: int) -> str:
    """Mesh title with prototype and model clock plus deform scale."""
    t_model = float(t_proto_to_model(t_proto))
    return (
        f"t = {t_model:.3f} s model  ({float(t_proto):.3f} s proto)   "
        f"deform x{SCALE:g}   {ifr + 1}/{nfr}"
    )


def _frame_hist_mm_formatter() -> FuncFormatter:
    """Integer mm tick labels (1 mm precision)."""
    return FuncFormatter(lambda v, _: f"{int(round(v))}")


def _finish_frame_hist_dual_axes(
    fig: plt.Figure,
    axes_h: list[plt.Axes],
    ylim_u: float,
) -> None:
    """Model-scale time (top of history column) and ux/λ (per-panel secondary y).

    Each history row keeps prototype ``u_x`` (mm) on the left. ``secondary_yaxis``
    on the right mirrors primary tick positions with labels ÷ λ; one row carries
    the shared y-label so ticks are not clipped by the figure edge.
    """
    mm_fmt = _frame_hist_mm_formatter()
    visible = [ax for ax in axes_h if ax.get_visible()]
    if not visible:
        return

    ax_top = visible[0]
    label_row = len(visible) // 2  # shared y-label only once

    sec_x = ax_top.secondary_xaxis(
        "top",
        functions=(
            lambda t_p: t_proto_to_model(t_p),
            lambda t_m: t_model_to_proto(t_m),
        ),
    )
    sec_x.set_xlabel(
        r"$t/\sqrt{\lambda}$ (s) model scale", fontsize=FRAME_FS_AXIS, labelpad=2
    )
    sec_x.tick_params(labelsize=FRAME_FS)

    for i, ax in enumerate(visible):
        ax.yaxis.set_major_locator(MaxNLocator(nbins=3, integer=True))
        ax.yaxis.set_major_formatter(mm_fmt)
        sec_y = ax.secondary_yaxis(
            "right",
            functions=(
                lambda u_p: np.asarray(u_p, dtype=float) / CYLINDER_LENGTH_SCALE,
                lambda u_m: np.asarray(u_m, dtype=float) * CYLINDER_LENGTH_SCALE,
            ),
        )
        ylo, yhi = ax.get_ylim()
        yticks = [
            y
            for y in ax.get_yticks()
            if ylo - 1e-9 <= float(y) <= yhi + 1e-9
        ]
        if yticks:
            sec_y.set_yticks(
                np.asarray(yticks, dtype=float) / CYLINDER_LENGTH_SCALE
            )
        sec_y.yaxis.set_major_formatter(mm_fmt)
        sec_y.tick_params(
            labelsize=FRAME_FS,
            pad=1,
            length=3,
            labelright=True,
        )
        if i == label_row:
            sec_y.set_ylabel(
                r"$u_x/\lambda$ (mm) model scale",
                fontsize=FRAME_FS_AXIS,
                labelpad=4,
            )

def _create_frame_hist_figure(
    t: np.ndarray,
    ux: np.ndarray,
    mesh: dict,
    traces: list[tuple[str, int, str]],
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    nfr: int,
) -> dict:
    """Build mesh + ux-history figure for animation or preview.

    Returns:
        Artist handles and arrays for ``_update_frame_hist_figure``.
    """
    idx = mesh["idx"]
    ntr = max(len(traces), 1)
    fig = plt.figure(figsize=(12.5, 8.2))
    gs = GridSpec(
        ntr,
        2,
        figure=fig,
        width_ratios=[1.05, 1.35],
        wspace=0.22,
        hspace=0.06,
        left=0.08,
        right=0.90,
        top=0.94,
        bottom=0.07,
    )
    ax_m = fig.add_subplot(gs[:, 0])
    ax_m.set_aspect("equal")
    ax_m.set_xlabel(r"$x$ (m) prototype scale", fontsize=FRAME_FS_AXIS)
    ax_m.set_ylabel(r"$y$ (m) prototype scale", fontsize=FRAME_FS_AXIS)
    ax_m.tick_params(labelsize=FRAME_FS)
    ax_m.grid(True, ls=":", alpha=0.35)
    ax_m.set_xlim(*xlim)
    ax_m.set_ylim(*ylim)
    ttl = ax_m.set_title("", fontsize=FRAME_FS_AXIS)

    q_idx = mesh["q_idx"]
    face = mesh["face"]
    pc = PolyCollection(
        np.zeros((max(len(q_idx), 1), 4, 2)),
        facecolors=face if face else "#cfd8dc",
        edgecolors="#333333",
        linewidths=0.18,
        zorder=0,
    )
    if len(q_idx) == 0:
        pc.set_visible(False)
    ax_m.add_collection(pc)
    lc_s = LineCollection(
        np.zeros((max(len(mesh["ij_s"]), 1), 2, 2)),
        colors=ORANGE, linewidths=1.25, zorder=3,
    )
    if len(mesh["ij_s"]) == 0:
        lc_s.set_visible(False)
    ax_m.add_collection(lc_s)
    lc_z = LineCollection(
        np.zeros((max(len(mesh["ij_z"]), 1), 2, 2)),
        colors=PURPLE, linewidths=1.55, zorder=4,
    )
    if len(mesh["ij_z"]) == 0:
        lc_z.set_visible(False)
    ax_m.add_collection(lc_z)

    prof = mesh["prof"]
    handles = [
        Patch(
            facecolor=to_rgba(layer_style(nm, profile=prof)["fill"], alpha=0.55),
            edgecolor="#333333",
            label=layer_style(nm, profile=prof).get("label", nm),
        )
        for nm in mesh["names_used"]
    ]
    handles.append(plt.Line2D([0], [0], color=ORANGE, lw=1.4, label="structure"))
    handles.append(plt.Line2D([0], [0], color=PURPLE, lw=1.6, label="springs"))
    ax_m.legend(handles=handles, loc="lower left", fontsize=FRAME_FS, framealpha=0.88)

    umax = 1e-4
    for _lab, tg, _c in traces:
        if tg not in idx:
            continue
        a = np.asarray(ux[:, idx[tg]], dtype=float)
        if np.isfinite(a).any():
            umax = max(umax, float(np.nanmax(np.abs(a))))
    ylim_mm = umax * M_TO_MM * 1.08

    markers = []
    for _lab, tg, col in traces:
        if tg not in idx:
            markers.append(None)
            continue
        ln, = ax_m.plot([], [], "o", ms=5, color=col, zorder=6,
                        markeredgecolor="k", markeredgewidth=0.35)
        markers.append(ln)

    hist = []
    axes_h = [fig.add_subplot(gs[i, 1]) for i in range(ntr)]
    for i, (lab, tg, col) in enumerate(traces):
        ax = axes_h[i]
        if tg not in idx:
            ax.set_visible(False)
            hist.append(None)
            continue
        u = ux[:, idx[tg]] * M_TO_MM
        ax.plot(t, u, color="#b0b0b0", lw=0.9, zorder=1)
        run_ln, = ax.plot([], [], color=col, lw=1.15, zorder=2)
        vln = ax.axvline(float(t[0]), color="#555", ls="--", lw=0.7, zorder=3)
        ax.axhline(0.0, color="#bbb", lw=0.5)
        ax.set_ylim(-ylim_mm, ylim_mm)
        ax.tick_params(labelsize=FRAME_FS, pad=2)
        ax.grid(True, ls=":", alpha=0.35)
        ax.text(
            0.97, 0.08, lab, transform=ax.transAxes,
            ha="right", va="bottom", fontsize=FRAME_FS, color=col, zorder=5,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      alpha=0.82, edgecolor="none"),
        )
        if i < ntr - 1:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel(r"$t$ (s) prototype scale", fontsize=FRAME_FS_AXIS)
        ax.set_xlim(float(t[0]), float(t[-1]))
        hist.append({"run": run_ln, "vline": vln})

    _finish_frame_hist_dual_axes(fig, axes_h, ylim_mm)

    pos = axes_h[ntr // 2].get_position()
    fig.supylabel(r"$u_x$ (mm) prototype scale", fontsize=FRAME_FS_SUP, x=pos.x0 - 0.045)

    return {
        "fig": fig,
        "ax_m": ax_m,
        "ttl": ttl,
        "pc": pc,
        "lc_s": lc_s,
        "lc_z": lc_z,
        "q_idx": q_idx,
        "ij_s": mesh["ij_s"],
        "ij_z": mesh["ij_z"],
        "fb": mesh["fb"],
        "markers": markers,
        "hist": hist,
        "traces": traces,
        "idx": idx,
        "X0": mesh["X0"],
        "Y0": mesh["Y0"],
        "t": t,
        "ux": ux,
        "nfr": nfr,
    }


def _update_frame_hist_figure(ctx: dict, k: int, ifr: int, ux_k: np.ndarray, uy_k: np.ndarray) -> None:
    """Refresh one animation frame (mesh deformation + running histories)."""
    X = ctx["X0"] + SCALE * ux_k
    Y = ctx["Y0"] + SCALE * uy_k
    q_idx = ctx["q_idx"]
    if len(q_idx):
        ctx["pc"].set_verts(np.stack((X[q_idx], Y[q_idx]), axis=-1))
    if len(ctx["ij_s"]):
        ctx["lc_s"].set_segments(_line_segments(X, Y, ctx["ij_s"]))
    if len(ctx["ij_z"]):
        ctx["lc_z"].set_segments(_spring_segments(X, Y, ctx["ij_z"], ctx["fb"], SPRING_MINLEN))
    t_now = float(ctx["t"][k])
    ctx["ttl"].set_text(frame_time_title(t_now, ifr, ctx["nfr"]))
    for (lab, tg, _col), mk, ha in zip(ctx["traces"], ctx["markers"], ctx["hist"]):
        if tg not in ctx["idx"] or mk is None or ha is None:
            continue
        i = ctx["idx"][tg]
        mk.set_data([X[i]], [Y[i]])
        u = ctx["ux"][:, i] * M_TO_MM
        ha["run"].set_data(ctx["t"][: k + 1], u[: k + 1])
        ha["vline"].set_xdata([t_now, t_now])


def frame_steps(t: np.ndarray) -> np.ndarray:
    """Select recorder samples for deformed-shape frames.

    Window starts at ``FRAME_T0_MODEL_S`` on the model clock. Frame count
    follows model-scale duration at ``FRAME_FPS`` (or ``N_FRAMES``).

    Args:
        t: Recorder times in s (prototype / ``t_num``).
    Returns:
        Sample indices; ``N_FRAMES`` takes precedence over ``FRAME_FPS``.
    """
    nt = len(t)
    if nt < 1:
        return np.zeros(0, dtype=int)
    t0 = frame_t0_proto()
    i0 = int(np.searchsorted(t, t0, side="left"))
    i0 = min(max(i0, 0), nt - 1)
    if i0 >= nt - 1:
        return np.array([nt - 1], dtype=int)
    win = np.arange(i0, nt, dtype=int)
    nw = len(win)
    if N_FRAMES and N_FRAMES > 0:
        nfr = min(int(N_FRAMES), nw)
        local = np.unique(np.linspace(0, nw - 1, nfr).round().astype(int))
        return win[local]
    fps = float(FRAME_FPS)
    if fps <= 0 or nw < 2:
        return win
    T_model = (float(t[-1]) - float(t[i0])) / TIME_SCALE_FROUDE
    if T_model <= 0:
        return win
    nfr = min(nw, max(2, int(round(T_model * fps)) + 1))
    local = np.unique(np.linspace(0, nw - 1, nfr).round().astype(int))
    return win[local]


def plot_frames(
    out: Path,
    t: np.ndarray,
    ux: np.ndarray,
    uy: np.ndarray,
    tags: list[int],
    xy: dict[int, tuple[float, float]],
    lines: list[list[int]],
    quads: list[list[int]],
    js: dict | None,
    meta: dict | None = None,
) -> None:
    """Render deformed soil, structure, spring frames and an MP4.

    When ``DO_FRAME_HIST`` and trace nodes exist, each frame includes running
    ``u_x(t)`` histories (gray full trace + colored segment to current time).

    Args:
        out: Plot output directory.
        t: Recorder times in s.
        ux: Horizontal nodal displacements in m.
        uy: Vertical nodal displacements in m.
        tags: Nodes in displacement-column order.
        xy: Undeformed node coordinates in m keyed by tag.
        lines: Structural and spring line connectivity.
        quads: Soil quad connectivity.
        js: Spring/model JSON, or ``None``.
        meta: Window metadata, or ``None``.
    Returns:
        None; writes frame PNGs and, when available, ``eq_window.mp4``.
    """
    d = out / "frames"
    d.mkdir(parents=True, exist_ok=True)
    for old in d.glob("*.png"):
        old.unlink()

    steps = frame_steps(t)
    nfr = len(steps)
    mesh = _frame_mesh_parts(tags, xy, lines, quads, js, meta)
    idx = mesh["idx"]
    amp = float(np.max(np.abs(ux[steps]))) if len(steps) else 0.0
    amp = max(amp, float(np.max(np.abs(uy[steps]))) if len(steps) else 0.0)
    pad = SCALE * amp + 0.5
    xlim = (float(mesh["X0"].min()) - pad, float(mesh["X0"].max()) + pad)
    ylim = (float(mesh["Y0"].min()) - pad, float(mesh["Y0"].max()) + pad)

    traces: list[tuple[str, int, str]] = []
    if DO_FRAME_HIST and meta is not None:
        traces = pick_frame_traces(tags, xy, meta, idx)

    ndig = max(4, len(str(max(nfr - 1, 0))))

    if traces:
        ctx = _create_frame_hist_figure(t, ux, mesh, traces, xlim, ylim, nfr)
        fig = ctx["fig"]
        for ifr, k in enumerate(steps):
            _update_frame_hist_figure(ctx, k, ifr, ux[k], uy[k])
            fig.savefig(
                d / f"frame_{ifr:0{ndig}d}.png",
                dpi=FRAME_DPI,
                facecolor="white",
                pil_kwargs={"compress_level": 1},
            )
            if ifr == 0 or (ifr + 1) % 400 == 0 or ifr + 1 == nfr:
                print(f"PlotEQ: frames {ifr + 1}/{nfr}", flush=True)
        plt.close(fig)
        print(f"PlotEQ: {nfr} frames (+ histories) -> {d}")
    else:
        fig, ax = plt.subplots(figsize=(6.0, 8.0))
        ax.set_aspect("equal")
        ax.set_xlabel(r"$x$ (m) prototype scale")
        ax.set_ylabel(r"$y$ (m) prototype scale")
        ax.grid(True, ls=":", alpha=0.35)
        ttl = ax.set_title("", fontsize=FRAME_FS_AXIS)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)

        q_idx = mesh["q_idx"]
        pc = PolyCollection(
            np.zeros((max(len(q_idx), 1), 4, 2)),
            facecolors=mesh["face"] if mesh["face"] else "#cfd8dc",
            edgecolors="#333333",
            linewidths=0.18,
            zorder=0,
        )
        if len(q_idx) == 0:
            pc.set_visible(False)
        ax.add_collection(pc)
        lc_s = LineCollection(
            np.zeros((max(len(mesh["ij_s"]), 1), 2, 2)),
            colors=ORANGE, linewidths=1.25, zorder=3,
        )
        if len(mesh["ij_s"]) == 0:
            lc_s.set_visible(False)
        ax.add_collection(lc_s)
        lc_z = LineCollection(
            np.zeros((max(len(mesh["ij_z"]), 1), 2, 2)),
            colors=PURPLE, linewidths=1.55, zorder=4,
        )
        if len(mesh["ij_z"]) == 0:
            lc_z.set_visible(False)
        ax.add_collection(lc_z)

        prof = mesh["prof"]
        handles = [
            Patch(facecolor=to_rgba(layer_style(nm, profile=prof)["fill"], alpha=0.55),
                  edgecolor="#333333", label=layer_style(nm, profile=prof).get("label", nm))
            for nm in mesh["names_used"]
        ]
        handles.append(plt.Line2D([0], [0], color=ORANGE, lw=1.4, label="structure"))
        handles.append(plt.Line2D([0], [0], color=PURPLE, lw=1.6, label="springs"))
        ax.legend(handles=handles, loc="lower left", fontsize=FRAME_FS, framealpha=0.88)
        fig.subplots_adjust(left=0.12, right=0.98, bottom=0.07, top=0.95)

        X0, Y0 = mesh["X0"], mesh["Y0"]
        for ifr, k in enumerate(steps):
            X = X0 + SCALE * ux[k]
            Y = Y0 + SCALE * uy[k]
            if len(q_idx):
                pc.set_verts(np.stack((X[q_idx], Y[q_idx]), axis=-1))
            if len(mesh["ij_s"]):
                lc_s.set_segments(_line_segments(X, Y, mesh["ij_s"]))
            if len(mesh["ij_z"]):
                lc_z.set_segments(_spring_segments(X, Y, mesh["ij_z"], mesh["fb"], SPRING_MINLEN))
            ttl.set_text(frame_time_title(float(t[k]), ifr, nfr))
            fig.savefig(
                d / f"frame_{ifr:0{ndig}d}.png",
                dpi=FRAME_DPI,
                facecolor="white",
                pil_kwargs={"compress_level": 1},
            )
            if ifr == 0 or (ifr + 1) % 400 == 0 or ifr + 1 == nfr:
                print(f"PlotEQ: frames {ifr + 1}/{nfr}", flush=True)
        plt.close(fig)
        print(f"PlotEQ: {nfr} frames -> {d}")

    if mux_frame_mp4(out, d, t, steps, ndig):
        prune_movie_frames(d, steps, t, ux, uy, tags, ndig)


def ffmpeg_bin() -> str | None:
    """Locate an ffmpeg executable.

    Args:
        None.
    Returns:
        Executable path, or ``None`` when ffmpeg is unavailable.
    """
    p = shutil.which("ffmpeg")
    if p:
        return p
    try:
        import imageio_ffmpeg
        return imageio_ffmpeg.get_ffmpeg_exe()
    except Exception:
        return None


def mux_frame_mp4(
    out: Path, d: Path, t: np.ndarray, steps: np.ndarray, ndig: int
) -> bool:
    """Encode frames so MP4 wall-clock duration matches model-scale window.

    Args:
        out: Plot output directory.
        d: Frame PNG directory.
        t: Recorder times in s (prototype).
        steps: Recorder sample indices used for frames.
        ndig: Zero-padding width in frame filenames.
    Returns:
        ``True`` when ``eq_window.mp4`` is encoded successfully.
    """
    nfr = len(steps)
    if nfr < 2:
        return False
    ff = ffmpeg_bin()
    if not ff:
        print("PlotEQ: no ffmpeg -- PNGs only")
        return False
    T_model = (float(t[steps[-1]]) - float(t[steps[0]])) / TIME_SCALE_FROUDE
    if T_model <= 0:
        T_model = float(nfr - 1) * 0.02
    in_fps = nfr / T_model
    out_fps = min(float(MOVIE_FPS), in_fps)
    mp4 = out / "eq_window.mp4"
    cmd = [
        ff, "-y", "-loglevel", "error",
        "-framerate", f"{in_fps:.8f}",
        "-i", str(d / f"frame_%0{ndig}d.png"),
        "-vf", f"fps={out_fps:.4f}",
        "-c:v", "libx264", "-pix_fmt", "yuv420p",
        "-crf", "18", "-movflags", "+faststart",
        str(mp4),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode == 0:
        t0m = t_proto_to_model(float(t[steps[0]]))
        t1m = t_proto_to_model(float(t[steps[-1]]))
        print(
            f"PlotEQ: movie {mp4}  ({out_fps:.1f} fps, {T_model:.2f} s model; "
            f"{t0m:.2f}–{t1m:.2f} s model)"
        )
        return True
    print(f"PlotEQ: ffmpeg skipped ({r.stderr.strip() or r.returncode})")
    return False


def prune_movie_frames(
    d: Path,
    steps: np.ndarray,
    t: np.ndarray,
    ux: np.ndarray,
    uy: np.ndarray,
    tags: list[int],
    ndig: int,
) -> None:
    """Keep key stills and remove intermediate PNGs after MP4 encoding.

    Args:
        d: Frame PNG directory.
        steps: Recorder sample indices used for frames.
        t: Recorder times in s.
        ux: Horizontal nodal displacements in m.
        uy: Vertical nodal displacements in m.
        tags: Nodes in displacement-column order.
        ndig: Zero-padding width in frame filenames.
    Returns:
        None; keeps start/end and pier-top extrema stills.
    """
    nfr = len(steps)
    if nfr < 1:
        return
    keep: dict[str, int] = {"t0": 0, "tend": nfr - 1}
    idx = {tg: i for i, tg in enumerate(tags)}
    if PIER_TOP in idx:
        col = idx[PIER_TOP]

        def near(k: int) -> int:
            """Find the rendered frame nearest one recorder sample.

            Args:
                k: Recorder sample index.
            Returns:
                Nearest frame index.
            """
            return int(np.argmin(np.abs(steps.astype(np.int64) - int(k))))

        keep["ux_max"] = near(int(np.argmax(ux[:, col])))
        keep["ux_min"] = near(int(np.argmin(ux[:, col])))
        keep["uy_max"] = near(int(np.argmax(uy[:, col])))
        keep["uy_min"] = near(int(np.argmin(uy[:, col])))
    kept = []
    for label, ifr in keep.items():
        src = d / f"frame_{ifr:0{ndig}d}.png"
        if not src.is_file():
            continue
        k = int(steps[ifr])
        dst = d / f"{label}_t{t[k]:.3f}s.png"
        shutil.copy2(src, dst)
        kept.append(dst.name)
    for ifr in range(nfr):
        p = d / f"frame_{ifr:0{ndig}d}.png"
        if p.is_file():
            p.unlink()
    print(f"PlotEQ: kept {len(kept)} frames ({', '.join(kept)})")


# ------------------------------------------------------------
# 7. COMPARISON, CLI, AND MAIN
# ------------------------------------------------------------

def load_pier_top(eq: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """Load the dedicated pier-top displacement recorder.

    Args:
        eq: Serial recorder directory.
    Returns:
        Time in s and pier-top ux/uy in m, or ``None`` when absent.
    """
    p = eq / "pier_top_disp.out"
    if not p.is_file():
        return None
    a = np.loadtxt(p)
    if a.ndim == 1:
        a = a.reshape(1, -1)
    t, ux, uy = a[:, 0], a[:, 1], a[:, 2]
    if SUBTRACT_T0:
        ux = ux - ux[0]
        uy = uy - uy[0]
    return t, ux, uy


def overlay_bcs(
    profile: int,
    soil_ele: str = "quad",
    pier_ele: str = "lumpedPlasticity",
) -> int:
    """Plot Shin and ASDEA pier-top history overlays.

    Args:
        profile: Soil-profile number.
        soil_ele: Soil element type.
        pier_ele: Pier element/hinge type.
    Returns:
        Process status: zero on success, one when an input recorder is missing.
    """
    shin = eq_dir(profile, "Shin", soil_ele, pier_ele, "serial")
    asd = eq_dir(profile, "ASDEA", soil_ele, pier_ele, "serial")
    a = load_pier_top(shin)
    b = load_pier_top(asd)
    if a is None or b is None:
        missing = []
        if a is None:
            missing.append(str(shin / "pier_top_disp.out"))
        if b is None:
            missing.append(str(asd / "pier_top_disp.out"))
        print("PlotEQ: overlay needs both BCs:\n  " + "\n  ".join(missing),
              file=sys.stderr)
        return 1
    t_s, ux_s, uy_s = a
    t_a, ux_a, uy_a = b
    out = eq_compare_dir(profile, pier_ele)
    out.mkdir(parents=True, exist_ok=True)
    meta = read_meta(shin) if (shin / "window_meta.txt").is_file() else {}
    t_eq = eq_end_time(meta, t_s)

    fig, ax = plt.subplots(figsize=(10.4, 4.2), constrained_layout=True)
    ax.plot(t_s, to_mm(ux_s), color="#1565c0", lw=1.35, label="Shin  pier top ux")
    ax.plot(t_a, to_mm(ux_a), color="#c45c12", lw=1.15, ls="--", label="ASDEA  pier top ux")
    mark_eq_end(ax, t_eq)
    ax.set_xlabel(r"$t_\mathrm{num}$ (s)")
    ax.set_ylabel("ux (mm)")
    ax.set_title(f"Profile {profile}  {pier_ele}  Shin vs ASDEA  (pier top, Δ from t0)")
    ax.legend(fontsize=8)
    ax.grid(True, ls=":", alpha=0.45)
    fig.savefig(out / "hist_ux_overlay.png", dpi=DPI)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10.4, 3.6), constrained_layout=True)
    ax.plot(t_s, to_mm(uy_s), color="#1565c0", lw=1.35, label="Shin  pier top uy")
    ax.plot(t_a, to_mm(uy_a), color="#c45c12", lw=1.15, ls="--", label="ASDEA  pier top uy")
    mark_eq_end(ax, t_eq)
    ax.set_xlabel(r"$t_\mathrm{num}$ (s)")
    ax.set_ylabel("uy (mm)")
    ax.set_title(f"Profile {profile}  {pier_ele}  Shin vs ASDEA  (pier top, Δ from t0)")
    ax.legend(fontsize=8)
    ax.grid(True, ls=":", alpha=0.45)
    fig.savefig(out / "hist_uy_overlay.png", dpi=DPI)
    plt.close(fig)
    print(f"PlotEQ: overlay -> {out}")
    print(f"  Shin  n={len(t_s)}  t={float(t_s[-1]):.2f} s")
    print(f"  ASDEA n={len(t_a)}  t={float(t_a[-1]):.2f} s")
    return 0


HELP = """\
usage: python3 plot/PlotEQ.py [eqOutDir]
       python3 plot/PlotEQ.py --overlay PROFILE [pierEleType]

  eqOutDir   serial recorder folder (default: EQ_OUT in this file)
  --overlay  Shin vs ASDEA overlay for that soil profile
  plots      lab dumps → LOCAL/plots/runs/<Test>/eq/; else <eqOutDir>/plots/
  MP dumps   python3 plot/PlotEQParallel.py [eqOutDir]
"""


def main() -> int:
    """Run the serial plot workflow selected by the command line and switches.

    Args:
        None; reads ``sys.argv`` without changing the existing CLI.
    Returns:
        Process status code.
    """
    if any(a in ("-h", "--help") for a in sys.argv[1:]):
        print(HELP, end="")
        return 0
    if len(sys.argv) >= 3 and sys.argv[1] == "--overlay":
        pier = sys.argv[3] if len(sys.argv) > 3 else "lumpedPlasticity"
        return overlay_bcs(int(sys.argv[2]), pier_ele=pier)
    eq = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(EQ_OUT)
    if not (eq / "window_nodes.txt").is_file():
        if (eq / "window_nodes.txt.0").is_file():
            print(
                f"PlotEQ: OpenSeesMP dump in {eq} (window_nodes.txt.0).\n"
                f"  python3 plot/PlotEQParallel.py {eq}",
                file=sys.stderr,
            )
        else:
            print(f"PlotEQ: missing window_nodes.txt in {eq}", file=sys.stderr)
        return 1

    meta = read_meta(eq)
    tags, xy = read_nodes(eq)
    disp_tags = read_disp_nodes(eq)
    idx = {t: i for i, t in enumerate(disp_tags)}
    lines, quads = read_eles(eq)
    disp_files = meta.get("dispFiles", "window_disp.out").split()
    # Frames: 0=auto when every window node has a disp column; 1=force; -1=off.
    full_window_disp = len(tags) > 0 and len(disp_tags) >= len(tags)
    if DO_FRAMES < 0:
        want_frames = False
    elif DO_FRAMES > 0:
        want_frames = True
    else:
        want_frames = full_window_disp
    need_disp = DO_HIST or DO_DEPTH_HIST or DO_ENVELOPE or want_frames
    if need_disp:
        t, ux, uy = load_window_disp(eq, disp_tags, disp_files)
        if SUBTRACT_T0:
            ux = maybe_t0(ux)
            uy = maybe_t0(uy)
    else:
        t = ux = uy = None

    if is_lab_dump(eq):
        out = run_eq_plots_dir(eq.name).resolve()
        print(f"PlotEQ: lab dump -> plots-out {out}")
    else:
        out = eq / "plots"
    out.mkdir(parents=True, exist_ok=True)
    # Pile groups drive the ux history and envelope, so key them off the nodes
    # that own a displacement column.
    groups = pile_groups({tg: xy[tg] for tg in disp_tags if tg in xy})
    js = load_spring_json(meta)
    if js is None:
        print("PlotEQ: pile_springs.json not found -- spring capacity overlays skipped")

    t_eq = eq_end_time(meta, t) if t is not None else None
    t_cut = truncated_end(meta, t)
    d595 = load_d595_proto()
    if d595 is not None:
        print(f"PlotEQ: D5-95 zoom  [{d595[0]:.1f}, {d595[1]:.1f}] s  (prototype / t_num)")
    if t is not None and len(t) > 1:
        print(
            f"PlotEQ: n={len(t)}  t={float(t[0]):.3g}..{float(t[-1]):.3g} s"
            + (f"  Trec={meta.get('Trec', '?')} s (incomplete)" if t_cut else "")
        )
    if DO_HIST:
        plot_hist(out, t, ux, uy, idx, groups, t_eq, t_cut, d595=d595)
        print(f"PlotEQ: wrote {out / 'hist_ux.png'}")
    if DO_DEPTH_HIST and t is not None:
        plot_depth_histories(
            out,
            eq,
            meta,
            js,
            t,
            ux,
            idx,
            groups,
            xy,
            disp_tags,
            t_eq,
            t_cut,
            d595=d595,
        )
    if DO_HINGE:
        plot_pier_hinge(out, eq, meta, t_eq, t_cut, d595=d595)
    if DO_PILE_SEC:
        plot_pile_section(out, eq, meta, xy)
    if DO_ENVELOPE:
        plot_envelope(out, ux, idx, groups, xy)
        print(f"PlotEQ: wrote {out / 'pile_envelope_ux.png'}")
    if DO_SPRING_ENV:
        plot_spring_envelopes(out, eq, js, meta, groups, xy)
        print(f"PlotEQ: spring envelopes -> {out}")
    if DO_HYST:
        plot_hyst(out, eq, js, meta)
        print(f"PlotEQ: hysteresis -> {out}")
    if DO_QUAD_PEAK:
        plot_quad_shear_peaks(out, eq, meta, xy, lines)
    if DO_QUAD_HYST:
        plot_quad_shear_hyst(out, eq, meta, xy)
    if want_frames:
        if not full_window_disp:
            print(
                f"PlotEQ: {len(tags) - len(disp_tags)} window nodes have no disp"
                " column (lean dump) -- skip frames"
            )
        elif t is None:
            print("PlotEQ: no window_disp -- skip frames")
        else:
            plot_frames(out, t, ux, uy, disp_tags, xy, lines, quads, js, meta)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
