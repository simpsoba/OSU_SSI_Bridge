#!/usr/bin/env python3
"""EQ window postprocess (serial dumps: window_nodes.txt, disp_nodes.txt, … ).

  python3 plot/PlotEQ.py
  python3 plot/PlotEQ.py /path/to/eqOutDir
  python3 plot/PlotEQ.py --overlay 3
  python3 plot/PlotEQ.py --overlay 3 lumpedPlasticity

OpenSeesMP dumps (name.$pid): python3 plot/PlotEQParallel.py [eqOutDir]

Writes PNGs to <eqOutDir>/plots/. Edit the block below; leave a switch at 0 to skip.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.tri import Triangulation
from matplotlib.colors import to_rgba
from matplotlib.patches import Patch

from paths import HERE, elevation_dir, eq_compare_dir, eq_dir, pile_springs_dir
from PlotModelSketch import layer_style

# ------------------------------------------------------------
# EDIT
# ------------------------------------------------------------
EQ_OUT = eq_dir(3, "Shin", "quad", "forceBeamColumn")

DO_HIST = 1
DO_ENVELOPE = 1       # pile-node ux min/max (symmetric)
DO_SPRING_ENV = 1     # spring defo vs y50/z50, force vs pult/tult/qult
DO_HYST = 1           # all p-y, t-z, and q-z loops
DO_QUAD_PEAK = 1      # window peak |tau_xy| and |gamma_xy|
DO_QUAD_HYST = 1      # tau_xy vs gamma_xy vs depth (one window column)
DO_HINGE = 1          # pier hinge hist + M-rot, P-axial, P-M, axial-rot
DO_PILE_SEC = 1       # pile M-kappa hyst; peak M and kappa vs depth
DO_FRAMES = 1         # window deform snapshots + MP4

N_FRAMES = 0          # >0 = that many equally spaced; 0 = use FRAME_FPS
FRAME_FPS = 30        # max PNGs per second of analysis (0 = every recorded sample)
SCALE = 20.0
SUBTRACT_T0 = 1       # nodal plots/frames relative to first sample
HYST_QUAD_X = None    # m, column centroid; None = outermost |x| in the window
PIER_TOP = 5
PIER_BOT = 1
DPI = 140
FRAME_DPI = 80
MOVIE_FPS = 30        # encode so MP4 duration = Trec (1:1)
SPRING_MINLEN = 0.30  # m, glyph floor so coincident zeroLength springs stay visible
# ------------------------------------------------------------

GRAY = "#90a4ae"
ORANGE = "#c45c12"
BLUE = "#1565c0"
BROWN = "#8B5A2B"
PURPLE = "#6a1b9a"
GREEN = "#2e7d32"
NAMES = ("L", "C", "R")
COLORS = {"L": BLUE, "C": BROWN, "R": PURPLE}


def _skip_hash(path: Path) -> list[str]:
    lines = []
    for ln in path.read_text().splitlines():
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        lines.append(s)
    return lines


def read_meta(eq: Path) -> dict[str, str]:
    meta = {}
    for ln in _skip_hash(eq / "window_meta.txt"):
        k, _, rest = ln.partition(" ")
        meta[k] = rest.strip()
    return meta


def read_node_file(path: Path) -> tuple[list[int], dict[int, tuple[float, float]]]:
    tags: list[int] = []
    xy: dict[int, tuple[float, float]] = {}
    for ln in _skip_hash(path):
        a, b, c = ln.split()[:3]
        t = int(a)
        tags.append(t)
        xy[t] = (float(b), float(c))
    return tags, xy


def read_nodes(eq: Path) -> tuple[list[int], dict[int, tuple[float, float]]]:
    """Geometry: every node with coordinates (window_nodes.txt)."""
    return read_node_file(eq / "window_nodes.txt")


def read_disp_nodes(eq: Path) -> list[int]:
    """Column order of window_disp*.out.

    recordersON=2 records displacement for the pier and the center pile only, so
    disp_nodes.txt is a subset of window_nodes.txt. Older dumps have no such
    file and every window node owns a column.
    """
    p = eq / "disp_nodes.txt"
    if not p.is_file():
        return read_nodes(eq)[0]
    return read_node_file(p)[0]


def read_eles(eq: Path) -> tuple[list[list[int]], list[list[int]]]:
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
    ux_parts = []
    uy_parts = []
    t = None
    for fn in disp_files:
        a = np.loadtxt(eq / fn)
        if a.ndim == 1:
            a = a.reshape(1, -1)
        ti = a[:, 0]
        if t is None:
            t = ti
        data = a[:, 1:]
        ux_parts.append(data[:, 0::2])
        uy_parts.append(data[:, 1::2])
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
    """np.loadtxt, but drop a truncated last line (live recorder)."""
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
    f = loadtxt_partial(eq / force_name)
    d = loadtxt_partial(eq / defo_name)
    if f.ndim == 1:
        f = f.reshape(1, -1)
        d = d.reshape(1, -1)
    ncomp = (f.shape[1] - 1) // n_ele
    F = f[:, 1:].reshape(f.shape[0], n_ele, ncomp)
    U = d[:, 1:].reshape(d.shape[0], n_ele, ncomp)
    return f[:, 0], F, U


def pile_groups(xy: dict[int, tuple[float, float]]) -> dict[str, list[int]]:
    piles = [t for t in xy if 2000 <= t < 3000]
    if not piles:
        return {}
    xs = sorted({round(xy[t][0], 2) for t in piles})
    names = list(NAMES[: len(xs)]) if len(xs) <= 3 else [f"p{i}" for i in range(len(xs))]
    out: dict[str, list[int]] = {n: [] for n in names}
    for t in piles:
        x = round(xy[t][0], 2)
        i = min(range(len(xs)), key=lambda k: abs(xs[k] - x))
        out[names[i]].append(t)
    for n in out:
        out[n].sort(key=lambda t: xy[t][1], reverse=True)
    return out


def maybe_t0(arr: np.ndarray) -> np.ndarray:
    if not SUBTRACT_T0 or arr.shape[0] < 1:
        return arr
    return arr - arr[0]


def sym_xlim(ax, *vals: np.ndarray, pad: float = 1.05) -> None:
    m = 0.0
    for v in vals:
        if v is None or len(np.atleast_1d(v)) == 0:
            continue
        m = max(m, float(np.nanmax(np.abs(v))))
    if m <= 0:
        m = 1.0
    ax.set_xlim(-pad * m, pad * m)


def layer_at_y(y: float, layers: list[dict]) -> str:
    if not layers:
        return "soil"

    def span(L: dict) -> tuple[float, float]:
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
    """(n, 2, 2) segments from soil->dup; short ZLs get a min-length tick."""
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
    if len(ij) == 0:
        return np.zeros((0, 2, 2))
    a = np.column_stack((X[ij[:, 0]], Y[ij[:, 0]]))
    b = np.column_stack((X[ij[:, 1]], Y[ij[:, 1]]))
    return np.stack([a, b], axis=1)


def eq_end_time(meta: dict, t: np.ndarray) -> float | None:
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


def mark_eq_end(ax, t_eq: float | None) -> None:
    if t_eq is None:
        return
    ax.axvline(t_eq, color="#78909c", lw=1.0, ls=":", label="EQ end")


def mark_last_sample(ax, t, t_eq: float | None) -> None:
    mark_eq_end(ax, t_eq)
    if t_eq is None and t is not None and len(t) > 2 and float(t[-1]) < 80.0:
        ax.axvline(float(t[-1]), color="#c62828", lw=1.0, ls="--",
                   label=f"last sample t={float(t[-1]):.2f} s")


def hyst_loop(ax, x, y, xlabel: str, ylabel: str, title: str) -> None:
    ax.plot(x, y, color=BROWN, lw=0.85, rasterized=True)
    ax.axhline(0.0, color="#9e9e9e", lw=0.6)
    ax.axvline(0.0, color="#9e9e9e", lw=0.6)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, ls=":", alpha=0.45)


def plot_hist(out: Path, t, ux, uy, idx, groups, t_eq: float | None = None) -> None:
    fig, ax = plt.subplots(figsize=(10.4, 4.2), constrained_layout=True)
    if PIER_TOP in idx:
        ax.plot(t, ux[:, idx[PIER_TOP]], color=ORANGE, lw=1.4,
                label=f"pier top ({PIER_TOP}) ux")
    if PIER_BOT in idx:
        ax.plot(t, ux[:, idx[PIER_BOT]], color=ORANGE, lw=1.0, ls="--",
                label=f"pier bot ({PIER_BOT}) ux")
    for name, tags in groups.items():
        if not tags:
            continue
        ax.plot(t, ux[:, idx[tags[0]]], color=COLORS.get(name, "#333"),
                lw=1.0, label=f"pile {name} head ux")
    mark_eq_end(ax, t_eq)
    if t_eq is None and t is not None and len(t) > 2 and float(t[-1]) < 80.0:
        ax.axvline(float(t[-1]), color="#c62828", lw=1.0, ls="--",
                   label=f"last sample t={float(t[-1]):.2f} s")
    ax.set_xlabel("t (s)")
    ax.set_ylabel("ux (m)")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, ls=":", alpha=0.45)
    fig.savefig(out / "hist_ux.png", dpi=DPI)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10.4, 3.6), constrained_layout=True)
    if PIER_TOP in idx:
        ax.plot(t, uy[:, idx[PIER_TOP]], color=ORANGE, lw=1.4, label="pier top uy")
    if PIER_BOT in idx:
        ax.plot(t, uy[:, idx[PIER_BOT]], color=ORANGE, lw=1.0, ls="--",
                label="pier bot uy")
    mark_eq_end(ax, t_eq)
    if t_eq is None and t is not None and len(t) > 2 and float(t[-1]) < 80.0:
        ax.axvline(float(t[-1]), color="#c62828", lw=1.0, ls="--",
                   label=f"last sample t={float(t[-1]):.2f} s")
    ax.set_xlabel("t (s)")
    ax.set_ylabel("uy (m)")
    ax.legend(fontsize=8)
    ax.grid(True, ls=":", alpha=0.45)
    fig.savefig(out / "hist_uy.png", dpi=DPI)
    plt.close(fig)


def plot_pier_hinge(
    out: Path, eq: Path, meta: dict, t_eq: float | None
) -> None:
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
        ax_lab = r"$\Delta u_\mathrm{ax}$ (m)" if SUBTRACT_T0 else r"$u_\mathrm{ax}$ (m)"
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
    P_kN = P / 1.0e3
    M_kNm = M / 1.0e3

    fig, axes = plt.subplots(2, 2, figsize=(10.4, 6.2), sharex="col",
                             constrained_layout=True)
    axes[0, 0].plot(t, ax_plot, color=ORANGE, lw=1.2)
    axes[1, 0].plot(t, rot_plot, color=ORANGE, lw=1.2)
    axes[0, 1].plot(t, P_kN, color=BLUE, lw=1.2)
    axes[1, 1].plot(t, M_kNm, color=BLUE, lw=1.2)
    axes[0, 0].set_ylabel(ax_lab)
    axes[1, 0].set_ylabel(rot_lab)
    axes[0, 1].set_ylabel("P (kN)")
    axes[1, 1].set_ylabel("Mz (kN·m)")
    axes[1, 0].set_xlabel("t (s)")
    axes[1, 1].set_xlabel("t (s)")
    axes[0, 0].set_title(f"Pier base hinge  ({kind})")
    for ax in axes.ravel():
        mark_last_sample(ax, t, t_eq)
        ax.grid(True, ls=":", alpha=0.45)
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
    fig, ax = plt.subplots(figsize=(5.6, 7.4), constrained_layout=True)
    xs = []
    for name, tags in groups.items():
        y = np.array([xy[t][1] for t in tags])
        col = np.array([ux[:, idx[t]] for t in tags])
        umin = col.min(axis=1)
        umax = col.max(axis=1)
        c = COLORS.get(name, "#333")
        ax.plot(umin, y, color=c, lw=1.5, label=f"{name} min")
        ax.plot(umax, y, color=c, lw=1.5, ls="--", label=f"{name} max")
        xs.extend([umin, umax])
    ax.axvline(0.0, color="#bbb", lw=0.8)
    sym_xlim(ax, *xs)
    ax.set_xlabel("ux min / max (m)")
    ax.set_ylabel("y (m), down is negative")
    ax.legend(fontsize=8)
    ax.grid(True, ls=":", alpha=0.45)
    fig.savefig(out / "pile_envelope_ux.png", dpi=DPI)
    plt.close(fig)


def read_pile_beam_eles(eq: Path) -> list[tuple[int, int, int]]:
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


def read_window_quad_list(eq: Path) -> list[int]:
    return [t for t, _ in read_window_quads(eq)]


def read_window_quads(eq: Path) -> list[tuple[int, str]]:
    rows: list[tuple[int, str]] = []
    p = eq / "window_quads.txt"
    if not p.is_file():
        return rows
    for ln in _skip_hash(p):
        a = ln.split()
        rows.append((int(a[0]), a[1] if len(a) > 1 else ""))
    return rows


def read_ele_nodes(eq: Path) -> dict[int, list[int]]:
    m: dict[int, list[int]] = {}
    p = eq / "window_eles.txt"
    if not p.is_file():
        return m
    for ln in _skip_hash(p):
        toks = [int(x) for x in ln.split()]
        m[toks[0]] = toks[1:]
    return m


def _file_ncols(path: Path) -> int:
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
    """Max over time and Gauss pts of |component icomp| (0-based in each GP)."""
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
    """GP-mean of component icomp vs time for global window_quads indices `keep`."""
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
    pts: list[tuple[float, float]] = []
    imap: dict[int, int] = {}

    def ix(n: int) -> int:
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
) -> None:
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
    tpc = ax.tripcolor(
        tri,
        facecolors=values[face],
        cmap="inferno",
        edgecolors="none",
        shading="flat",
    )
    vmax = float(np.nanmax(values))
    if vmax <= 0:
        vmax = 1.0
    tpc.set_clim(0.0, vmax)
    ax.tricontour(tri, znode, levels=8, colors="k", linewidths=0.35, alpha=0.45)
    segs = []
    for a, b in lines:
        if a in xy and b in xy and max(a, b) < 10000:
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
    sig_files = meta.get("quadStressFiles", "").split()
    eps_files = meta.get("quadStrainFiles", "").split()
    if sig_files:
        tau = peak_abs_gp_comp(eq, sig_files, n_ele, n_gp)[idx]
        plot_window_peak_field(
            out, "window_peak_tau_xy.png", xy, quads, tau / 1.0e3, lines,
            r"peak $|\tau_{xy}|$ (kPa)",
            r"peak $|\tau_{xy}|$  (max over $t$ and Gauss pts)",
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
        )
        print(
            f"PlotEQ: wrote {out / 'window_peak_gamma_xy.png'}  "
            f"max={float(gam.max()) * 1e3:.3g}e-3"
        )
    else:
        print("PlotEQ: no quadStrainFiles -- skip gamma contour")


def quad_centroids(
    qrows: list[tuple[int, str]],
    ev: dict[int, list[int]],
    xy: dict[int, tuple[float, float]],
) -> list[tuple[int, str, float, float, int]]:
    """(eleTag, layer, xc, yc, global index) for window quads with 4 nodes in xy."""
    out = []
    for i, (tag, nm) in enumerate(qrows):
        nn = ev.get(tag, [])
        if len(nn) < 4 or any(n not in xy for n in nn[:4]):
            continue
        xc = 0.25 * sum(xy[n][0] for n in nn[:4])
        yc = 0.25 * sum(xy[n][1] for n in nn[:4])
        out.append((tag, nm, xc, yc, i))
    return out


def quad_depth_column(
    cents: list[tuple[int, str, float, float, int]],
    x_target: float | None,
) -> list[tuple[int, str, float, float, int]]:
    """One ele per soil row, nearest x_target (default: outermost |x|)."""
    if not cents:
        return []
    if x_target is None:
        x_target = max(cents, key=lambda r: abs(r[2]))[2]
    by_y: dict[float, list[tuple[int, str, float, float, int]]] = {}
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
    col = quad_depth_column(cents, HYST_QUAD_X)
    if not col:
        print("PlotEQ: no window column for tau-gamma hyst")
        return
    n_gp = int(float(meta.get("quadNgp", 4)))
    n_ele = len(qrows)
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
    titles = [f"{nm or 'quad'}  y={yc:.2f} m" for _, nm, _, yc, _ in col]
    hyst_grid(
        out, "hyst_tau_gamma.png", U, F, 0, titles,
        r"$\gamma_{xy}$ ($\times 10^{-3}$)", r"$\tau_{xy}$ (kPa)", 3,
        share_x=True,
    )
    print(
        f"PlotEQ: wrote {out / 'hyst_tau_gamma.png'}  "
        f"n={len(col)}  x≈{xc:.2f} m  (GP mean, Δ from t0={SUBTRACT_T0})"
    )


def load_spring_json(meta: dict) -> dict | None:
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
    """PultCap, TultCap, QultSoffit, y50_cap, z50_cap (Mokwa, same as BuildSoilSprings)."""
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
    s = float(sizes.get("s_pile_cap", 1.8288))
    wsoil = float(sizes.get("W_cap_soil", 2.0 * s))
    xf = 0.5 * wsoil
    xs = np.array([-s, -0.5 * s, 0.0, 0.5 * s, s])
    if len(xs) != nsof:
        return np.linspace(-xf, xf, nsof)
    return xs


def split_cap_eles(eq: Path) -> tuple[list[int], list[int]]:
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
    out = []
    for ln in _skip_hash(path):
        out.append(int(ln.split()[0]))
    return out


def read_pile_spring_eles(eq: Path) -> list[dict]:
    """Rows from pile_springs_eles.txt: e, kind, optional ip iy."""
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


def align_pile_stations(
    rows: list[dict], stations: list[dict]
) -> list[dict] | None:
    """Stations in recorder-column order. None if they cannot be matched."""
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
    """arr (nt, n, ncomp) -> (umin, umax) length n."""
    v = arr[:, :, icomp]
    return v.min(axis=0), v.max(axis=0)


def pile_name(ip: int) -> str:
    return NAMES[ip] if ip < len(NAMES) else f"p{ip}"


def is_qz_station(s: dict) -> bool:
    return bool(s.get("isTip")) or s.get("axType") == "qz"


def pile_hyst_order(stations: list, n_ele: int, qz: bool | None = None) -> list[int]:
    """Row = depth (iy), cols L/C/R. qz True = tips only, False = shaft only."""
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
    m = max(float(np.nanmax(np.abs(v))) for v in ylim_vals)
    if m <= 0:
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
    rows = read_pile_spring_eles(eq)
    pile_eles = [r["e"] for r in rows]
    _, Fp, Up = load_spring_pt(
        eq, "pile_springs_force.out", "pile_springs_defo.out", len(pile_eles)
    )
    n = len(pile_eles)
    stations = align_pile_stations(rows, (js or {}).get("stations") or [])
    if stations is None or len(stations) != n:
        print(
            f"PlotEQ: pile stations {0 if stations is None else len(stations)}"
            f" != nEle {n} -- skip spring env"
        )
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

    def as_groups(minmax):
        umin, umax = minmax
        return [(name, y, umin[idx], umax[idx]) for name, idx, y in g]

    u_py = env_minmax(Up, 0)
    f_py = env_minmax(Fp, 0)
    plot_depth_env(
        out, "spring_env_py_defo.png", y0, y50,
        "u_py min / max (m)", r"$\pm y_{50}$", as_groups(u_py),
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

            def as_shaft(minmax):
                umin, umax = minmax
                return [(nm, yi, umin[idx], umax[idx]) for nm, idx, yi in shaft_g]

            plot_depth_env(
                out, "spring_env_tz_defo.png", y_s0, z50_s,
                "u_tz min / max (m)", r"$\pm z_{50}$ (shaft t-z)",
                as_shaft(u_tz),
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
                u_tz[0][tip], u_tz[1][tip], z50q,
                "u_qz min / max (m)", r"$z_{50}$ compression",
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

        def cap_groups(minmax, yf=yf, n2=n2):
            umin, umax = minmax
            if not n2:
                return [("cap", yf, umin, umax)]
            return [
                ("L", yf[:n2], umin[:n2], umax[:n2]),
                ("R", yf[n2:], umin[n2:], umax[n2:]),
            ]

        plot_depth_env(
            out, "spring_env_cap_py_defo.png", yL, y50c_v[:n2] if n2 else y50c_v,
            "u_py min / max (m)", r"$\pm y_{50}$ cap", cap_groups(env_minmax(Uc, 0)),
        )
        plot_depth_env(
            out, "spring_env_cap_py_force.png", yL, pcap[:n2] if n2 else pcap,
            "p min / max (N)", r"$\pm p_\mathrm{ult}$ cap", cap_groups(env_minmax(Fc, 0)),
        )
        if Uc.shape[2] > 1:
            plot_depth_env(
                out, "spring_env_cap_tz_defo.png", yL, z50c_v[:n2] if n2 else z50c_v,
                "u_tz min / max (m)", r"$\pm z_{50}$ cap", cap_groups(env_minmax(Uc, 1)),
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
            out, "spring_env_qz_defo.png", xq, umin, umax, z50q,
            "u_qz min / max (m)", r"$z_{50}$ compression",
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
    rows = read_pile_spring_eles(eq)
    pile_eles = [r["e"] for r in rows]
    _, Fp, Up = load_spring_pt(
        eq, "pile_springs_force.out", "pile_springs_defo.out", len(pile_eles)
    )
    stations = align_pile_stations(rows, (js or {}).get("stations") or []) or []
    sizes = load_sketch_sizes(meta)
    titles = []
    for i, e in enumerate(pile_eles):
        if i < len(stations):
            s = stations[i]
            titles.append(f"{pile_name(int(s['ip']))}  y={float(s['y']):.2f}")
        else:
            titles.append(f"ele {e}")
    hyst_grid(
        out, "hyst_pile_py.png", Up, Fp, 0, titles, "u_py (m)", "p (N)", 3,
        which=pile_hyst_order(stations, Up.shape[1]),
        share_x=True,
    )
    if Up.shape[2] > 1:
        hyst_grid(
            out, "hyst_pile_tz.png", Up, Fp, 1, titles,
            "u_tz (m)", "t (N)", 3, which=pile_hyst_order(stations, Up.shape[1], qz=False),
            share_x=True,
        )
        qz_i = pile_hyst_order(stations, Up.shape[1], qz=True)
        if qz_i:
            hyst_grid(
                out, "hyst_pile_qz.png", Up, Fp, 1, titles,
                "u_qz (m)", "q (N)", 3, which=qz_i,
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
        hyst_grid(
            out, "hyst_cap_py.png", Uc, Fc, 0, tcap, "u_py (m)", "p (N)", 2,
            which=cap_ord,
            share_x=True,
        )
        if Uc.shape[2] > 1:
            hyst_grid(
                out, "hyst_cap_tz.png", Uc, Fc, 1, tcap, "u_tz (m)", "t (N)", 2,
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
        hyst_grid(
            out, "hyst_soffit_qz.png", Uq, Fq, ic, tsof, "u_qz (m)", "q (N)", 4,
            share_x=True,
        )


def frame_steps(t: np.ndarray) -> np.ndarray:
    """Indices into the recorder. N_FRAMES wins; else FRAME_FPS of analysis time."""
    nt = len(t)
    if nt < 1:
        return np.zeros(0, dtype=int)
    if N_FRAMES and N_FRAMES > 0:
        nfr = min(int(N_FRAMES), nt)
        return np.unique(np.linspace(0, nt - 1, nfr).round().astype(int))
    fps = float(FRAME_FPS)
    if fps <= 0 or nt < 2:
        return np.arange(nt, dtype=int)
    T = float(t[-1] - t[0])
    if T <= 0:
        return np.arange(nt, dtype=int)
    nfr = min(nt, max(2, int(round(T * fps)) + 1))
    return np.unique(np.linspace(0, nt - 1, nfr).round().astype(int))


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
) -> None:
    d = out / "frames"
    d.mkdir(parents=True, exist_ok=True)
    for old in d.glob("*.png"):
        old.unlink()

    steps = frame_steps(t)
    nfr = len(steps)

    idx = {tg: i for i, tg in enumerate(tags)}
    X0 = np.array([xy[tg][0] for tg in tags])
    Y0 = np.array([xy[tg][1] for tg in tags])

    soil_quads = []
    for q in quads:
        if len(q) < 4:
            continue
        nn = q[:4]
        if min(nn) >= 10000 and max(nn) < 20000 and all(n in idx for n in nn):
            soil_quads.append(nn)
    q_idx = (
        np.array([[idx[n] for n in q] for q in soil_quads], dtype=int)
        if soil_quads
        else np.zeros((0, 4), dtype=int)
    )

    layers = list(
        (js or {}).get("soil_layers") or (js or {}).get("layers") or []
    )
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

    struct = []
    spr = []
    for ln in lines:
        if len(ln) < 2:
            continue
        a, b = ln[0], ln[1]
        hi, lo = max(a, b), min(a, b)
        if 20000 <= hi < 30000:
            # soil node first so fallback/segment order is soil -> dup
            spr.append((lo, hi) if lo < 20000 else (a, b))
        elif hi < 10000:
            struct.append((a, b))
    ij_s = _pair_idx(struct, idx)
    ij_z = _pair_idx(spr, idx)
    fb = np.zeros((len(ij_z), 2))
    fb[:, 0] = 1.0
    spr_base = min((n for n in tags if n >= 20000), default=20000)
    for k, (a, b) in enumerate(spr):
        dup = b if b >= 20000 else a
        if dup - spr_base >= 920:
            fb[k] = (0.0, 1.0)

    amp = float(np.max(np.abs(ux[steps])))
    amp = max(amp, float(np.max(np.abs(uy[steps]))))
    pad = SCALE * amp + 0.5
    xlim = (float(X0.min()) - pad, float(X0.max()) + pad)
    ylim = (float(Y0.min()) - pad, float(Y0.max()) + pad)

    fig, ax = plt.subplots(figsize=(6.0, 8.0))
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.grid(True, ls=":", alpha=0.35)
    ttl = ax.set_title("", fontsize=10)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

    pc = PolyCollection(
        np.zeros((max(len(q_idx), 1), 4, 2)),
        facecolors=face if face else "#cfd8dc",
        edgecolors="#333333",
        linewidths=0.18,
        zorder=0,
    )
    if len(q_idx) == 0:
        pc.set_visible(False)
    ax.add_collection(pc)
    lc_s = LineCollection(
        np.zeros((max(len(ij_s), 1), 2, 2)),
        colors=ORANGE, linewidths=1.25, zorder=3,
    )
    if len(ij_s) == 0:
        lc_s.set_visible(False)
    ax.add_collection(lc_s)
    lc_z = LineCollection(
        np.zeros((max(len(ij_z), 1), 2, 2)),
        colors=PURPLE, linewidths=1.55, zorder=4,
    )
    if len(ij_z) == 0:
        lc_z.set_visible(False)
    ax.add_collection(lc_z)

    handles = [
        Patch(facecolor=to_rgba(layer_style(nm, profile=prof)["fill"], alpha=0.55),
              edgecolor="#333333", label=layer_style(nm, profile=prof).get("label", nm))
        for nm in names_used
    ]
    handles.append(plt.Line2D([0], [0], color=ORANGE, lw=1.4, label="structure"))
    handles.append(plt.Line2D([0], [0], color=PURPLE, lw=1.6, label="springs"))
    ax.legend(handles=handles, loc="lower left", fontsize=7, framealpha=0.88)
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.07, top=0.95)

    ndig = max(4, len(str(max(nfr - 1, 0))))
    for ifr, k in enumerate(steps):
        X = X0 + SCALE * ux[k]
        Y = Y0 + SCALE * uy[k]
        if len(q_idx):
            verts = np.stack((X[q_idx], Y[q_idx]), axis=-1)
            pc.set_verts(verts)
        if len(ij_s):
            lc_s.set_segments(_line_segments(X, Y, ij_s))
        if len(ij_z):
            lc_z.set_segments(_spring_segments(X, Y, ij_z, fb, SPRING_MINLEN))
        ttl.set_text(f"t = {t[k]:.3f} s   scale x{SCALE:g}   {ifr + 1}/{nfr}")
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
    nfr = len(steps)
    if nfr < 2:
        return False
    ff = ffmpeg_bin()
    if not ff:
        print("PlotEQ: no ffmpeg -- PNGs only")
        return False
    T = float(t[steps[-1]] - t[steps[0]])
    if T <= 0:
        T = float(nfr - 1) * 0.02
    in_fps = nfr / T
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
        print(f"PlotEQ: movie {mp4}  ({out_fps:.1f} fps, {T:.2f} s)")
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
    """After a successful MP4, keep t0 / tend / pier-top ux,uy extrema; drop the rest."""
    nfr = len(steps)
    if nfr < 1:
        return
    keep: dict[str, int] = {"t0": 0, "tend": nfr - 1}
    idx = {tg: i for i, tg in enumerate(tags)}
    if PIER_TOP in idx:
        col = idx[PIER_TOP]

        def near(k: int) -> int:
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


def load_pier_top(eq: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
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
    """Shin vs ASDEA pier-top histories for one soil profile and pier type."""
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
    ax.plot(t_s, ux_s, color="#1565c0", lw=1.35, label="Shin  pier top ux")
    ax.plot(t_a, ux_a, color="#c45c12", lw=1.15, ls="--", label="ASDEA  pier top ux")
    mark_eq_end(ax, t_eq)
    ax.set_xlabel("t (s)")
    ax.set_ylabel("ux (m)")
    ax.set_title(f"Profile {profile}  {pier_ele}  Shin vs ASDEA  (pier top, Δ from t0)")
    ax.legend(fontsize=8)
    ax.grid(True, ls=":", alpha=0.45)
    fig.savefig(out / "hist_ux_overlay.png", dpi=DPI)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10.4, 3.6), constrained_layout=True)
    ax.plot(t_s, uy_s, color="#1565c0", lw=1.35, label="Shin  pier top uy")
    ax.plot(t_a, uy_a, color="#c45c12", lw=1.15, ls="--", label="ASDEA  pier top uy")
    mark_eq_end(ax, t_eq)
    ax.set_xlabel("t (s)")
    ax.set_ylabel("uy (m)")
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
  plots      <eqOutDir>/plots/   (panels: DO_* switches above)
  MP dumps   python3 plot/PlotEQParallel.py [eqOutDir]
"""


def main() -> int:
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
    need_disp = DO_HIST or DO_ENVELOPE or DO_FRAMES
    if need_disp:
        t, ux, uy = load_window_disp(eq, disp_tags, disp_files)
        if SUBTRACT_T0:
            ux = maybe_t0(ux)
            uy = maybe_t0(uy)
    else:
        t = ux = uy = None

    out = eq / "plots"
    out.mkdir(parents=True, exist_ok=True)
    # Pile groups drive the ux history and envelope, so key them off the nodes
    # that own a displacement column.
    groups = pile_groups({tg: xy[tg] for tg in disp_tags if tg in xy})
    js = load_spring_json(meta)
    if js is None:
        print("PlotEQ: pile_springs.json not found -- spring capacity overlays skipped")

    t_eq = eq_end_time(meta, t) if t is not None else None
    if DO_HIST:
        plot_hist(out, t, ux, uy, idx, groups, t_eq)
        print(f"PlotEQ: wrote {out / 'hist_ux.png'}")
    if DO_HINGE:
        plot_pier_hinge(out, eq, meta, t_eq)
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
    if DO_FRAMES:
        # Soil patches need a displacement per corner; without them the frames
        # would mix a deformed structure with an undeformed mesh.
        if len(disp_tags) < len(tags):
            print(
                f"PlotEQ: {len(tags) - len(disp_tags)} window nodes have no disp"
                " column (recordersON=2) -- skip frames"
            )
        else:
            plot_frames(out, t, ux, uy, disp_tags, xy, lines, quads, js)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
