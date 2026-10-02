#!/usr/bin/env python3
"""
Goals
-----
Compare offline baseline EQ dumps under ``plot/out/eq_offline/``:

  * pier top/base UX histories (full + D5–95 zoom)
  * hinge M–θ histories and hysteresis (θ ≈ κ·Lp for forceBeamColumn)
  * deformed-window snapshots at free-FBC peak time (all cases) and at each
    case's own pier-top peak
  * pointer plots for pile springs / pile shafts (per-case PlotEQ outputs)

Writes only to ``plot/out/eq_offline/compare/``. Never touches lab data trees.

Usage
-----
  python plot/PlotEQOfflineCompare.py
  python plot/PlotEQOfflineCompare.py --integrators newmark
  python plot/PlotEQOfflineCompare.py --root plot/out/eq_offline
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import PlotEQ as peq
from gm_duration import d595_proto_window

REPO = Path(__file__).resolve().parents[1]
DEFAULT_ROOT = REPO / "plot" / "out" / "eq_offline"

STRUCT_SLUGS = ("lumped_hold", "lumped_free", "fbc_hold", "fbc_free")
COLORS = {
    "lumped_hold": "#1f77b4",
    "lumped_free": "#ff7f0e",
    "fbc_hold": "#2ca02c",
    "fbc_free": "#d62728",
}
LABELS = {
    "lumped_hold": "lumped + hold",
    "lumped_free": "lumped + free",
    "fbc_hold": "FBC + hold",
    "fbc_free": "FBC + free",
}

# Priestley Lp when meta lacks Lp_pier (same formula as Parameters.tcl).
_INCH = 0.0254
_H_PIER = 3.02 * (48.0 / 20.0)
_FY_MPA = 470.0
_DB = (10.0 / 8.0) * _INCH
DEFAULT_LP = max(0.08 * _H_PIER + 0.022 * _FY_MPA * _DB, 0.044 * _FY_MPA * _DB)


def case_dir(root: Path, struct: str, integ: str) -> Path:
    return root / f"{struct}_{integ}"


def load_pier_ux(eq: Path, tag: int) -> tuple[np.ndarray, np.ndarray] | None:
    """Return (t, ux) from pier_node_<tag>.out, or None if missing."""
    path = eq / f"pier_node_{tag}.out"
    if not path.is_file():
        return None
    a = peq.loadtxt_partial(path)
    if a.size == 0 or a.shape[1] < 2:
        return None
    t = a[:, 0]
    ux = a[:, 1]
    if peq.SUBTRACT_T0:
        ux = ux - ux[0]
    return t, ux


def load_pier_series(
    eq: Path, tag: int, kind: str
) -> tuple[np.ndarray, np.ndarray] | None:
    """Return (t, channel_x) for disp/vel/accel on pier node ``tag``.

    kind: ``disp`` -> pier_node_<tag>.out
          ``vel``  -> pier_node_<tag>_vel.out
          ``acc``  -> pier_node_<tag>_acc.out
    """
    if kind == "disp":
        path = eq / f"pier_node_{tag}.out"
    elif kind == "vel":
        path = eq / f"pier_node_{tag}_vel.out"
    elif kind == "acc":
        path = eq / f"pier_node_{tag}_acc.out"
    else:
        raise ValueError(kind)
    if not path.is_file():
        return None
    a = peq.loadtxt_partial(path)
    if a.size == 0 or a.shape[1] < 2:
        return None
    t = a[:, 0]
    x = a[:, 1]
    if kind == "disp" and peq.SUBTRACT_T0:
        x = x - x[0]
    return t, x


def peak_abs_ux(t: np.ndarray, ux: np.ndarray) -> tuple[float, float, int]:
    """Return (t_peak, ux_at_peak, index) for max |ux|."""
    k = int(np.argmax(np.abs(ux)))
    return float(t[k]), float(ux[k]), k


def lp_from_meta(meta: dict) -> float:
    raw = meta.get("Lp_pier")
    if raw is None or raw == "":
        return float(DEFAULT_LP)
    return float(raw)


def load_hinge(
    eq: Path, meta: dict, *, top: bool = False
) -> tuple[np.ndarray, np.ndarray, np.ndarray, str] | None:
    """
    Load hinge Mz and rotation-like deformation.

    Returns: (t, M_Nm, theta_rad, kind) where kind is 'theta' or 'kappa*Lp'.
    For forceBeamColumn, θ ≈ κ · Lp (approximate plastic rotation).
    """
    kind = meta.get("pierHinge", "")
    if top:
        if kind != "lumpedPlasticity":
            return None
        fn_f = eq / meta.get("hingeTopForceFile", "pier_hinge_top_force.out")
        fn_d = eq / meta.get("hingeTopDefoFile", "pier_hinge_top_defo.out")
    else:
        fn_f = eq / meta.get("hingeForceFile", "pier_hinge_force.out")
        fn_d = eq / meta.get("hingeDefoFile", "pier_hinge_defo.out")
    if not fn_f.is_file() or not fn_d.is_file():
        return None
    F = peq.loadtxt_partial(fn_f)
    D = peq.loadtxt_partial(fn_d)
    if F.size == 0 or D.size == 0:
        return None
    n = min(len(F), len(D))
    t = F[:n, 0]
    M = F[:n, 2] if F.shape[1] > 2 else np.zeros(n)
    rot = D[:n, 2] if D.shape[1] > 2 else np.zeros(n)
    if peq.SUBTRACT_T0:
        rot = rot - rot[0]
    if kind == "forceBeamColumn":
        lp = lp_from_meta(meta)
        rot = rot * lp
        lab = "kappa*Lp"
    else:
        lab = "theta"
    return t, M, rot, lab


def _interp_align(
    t_a: np.ndarray, y_a: np.ndarray, t_b: np.ndarray, y_b: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Align two series onto the overlapping time grid of ``t_a``."""
    t0 = max(float(t_a[0]), float(t_b[0]))
    t1 = min(float(t_a[-1]), float(t_b[-1]))
    mask = (t_a >= t0) & (t_a <= t1)
    t = t_a[mask]
    ya = y_a[mask]
    yb = np.interp(t, t_b, y_b)
    return t, ya, yb


def plot_pier_histories(
    out: Path,
    root: Path,
    integ: str,
    d595: tuple[float, float] | None,
) -> None:
    """Overlay pier top/base UX (and VX/AX), plus UX drift = top − base.

    Also writes per-case ``plots/hist_ux.png`` (pier only) and
    ``plots/hist_drift_ux.png`` so the busy pile-head PlotEQ hist is replaced
    for these offline dumps.
    """
    # --- per-case pier-only UX + drift (no re-analysis) ---
    for struct in STRUCT_SLUGS:
        eq = case_dir(root, struct, integ)
        if not (eq / "window_meta.txt").is_file():
            continue
        top = load_pier_series(eq, 5, "disp")
        bot = load_pier_series(eq, 1, "disp")
        if top is None or bot is None:
            continue
        t, u_top, u_bot = _interp_align(top[0], top[1], bot[0], bot[1])
        drift = u_top - u_bot
        plots = eq / "plots"
        plots.mkdir(parents=True, exist_ok=True)
        meta = peq.read_meta(eq)
        t_eq = peq.eq_end_time(meta, t)
        t_cut = peq.truncated_end(meta, t)

        fig, axes_f, axes_z = peq.subplots_full_zoom(1, fig_h=4.2, sharey=True)
        for ax in (axes_f[0], axes_z[0]):
            ax.plot(t, peq.to_mm(u_top), color=peq.ORANGE, lw=1.4, label="pier top")
            ax.plot(
                t, peq.to_mm(u_bot), color=peq.ORANGE, lw=1.0, ls="--", label="pier bot"
            )
        peq.finish_full_zoom_pair(
            axes_f[0], axes_z[0], d595, t, t_eq, t_cut,
            full_ylim=peq.YLIM_DISP_PROTO_MM,
        )
        axes_f[0].set_ylabel(r"$\Delta u_x$ (mm)")
        axes_f[0].legend(fontsize=8)
        axes_f[0].set_title(f"{LABELS[struct]} — pier UX")
        p_ux = plots / "hist_ux.png"
        fig.savefig(p_ux, dpi=peq.DPI)
        plt.close(fig)
        print(f"PlotEQOfflineCompare: wrote {p_ux}  (pier top/bot only)")

        fig, axes_f, axes_z = peq.subplots_full_zoom(1, fig_h=4.2, sharey=True)
        for ax in (axes_f[0], axes_z[0]):
            ax.plot(t, peq.to_mm(drift), color=peq.BLUE, lw=1.3, label=r"$u_\mathrm{top}-u_\mathrm{bot}$")
            ax.axhline(0.0, color="#bbb", lw=0.8)
        peq.finish_full_zoom_pair(
            axes_f[0], axes_z[0], d595, t, t_eq, t_cut,
            full_ylim=peq.YLIM_DISP_PROTO_MM,
        )
        axes_f[0].set_ylabel(r"drift $\Delta u_x$ (mm)")
        axes_f[0].legend(fontsize=8)
        axes_f[0].set_title(f"{LABELS[struct]} — pier UX drift")
        p_dr = plots / "hist_drift_ux.png"
        fig.savefig(p_dr, dpi=peq.DPI)
        plt.close(fig)
        print(f"PlotEQOfflineCompare: wrote {p_dr}")

    # --- compare overlays: UX / VX / AX ---
    series = (
        ("disp", r"$\Delta u_x$ (mm)", True, "ux"),
        ("vel", r"$v_x$ (m/s)", False, "vx"),
        ("acc", r"$a_x$ (m/s$^2$)", False, "ax"),
    )
    for kind, ylab, to_mm, stem in series:
        fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.2), constrained_layout=True)
        ax_top_f, ax_top_z = axes[0, 0], axes[0, 1]
        ax_bot_f, ax_bot_z = axes[1, 0], axes[1, 1]
        n_ok = 0
        for struct in STRUCT_SLUGS:
            eq = case_dir(root, struct, integ)
            if not (eq / "window_meta.txt").is_file():
                continue
            top = load_pier_series(eq, 5, kind)
            bot = load_pier_series(eq, 1, kind)
            col = COLORS[struct]
            lab = LABELS[struct]
            for ax_f, ax_z, ser in (
                (ax_top_f, ax_top_z, top),
                (ax_bot_f, ax_bot_z, bot),
            ):
                if ser is None:
                    continue
                t, x = ser
                y = peq.to_mm(x) if to_mm else x
                ax_f.plot(t, y, color=col, lw=1.1, label=lab)
                ax_z.plot(t, y, color=col, lw=1.1, label=lab)
                n_ok += 1
        if n_ok == 0:
            plt.close(fig)
            continue
        for ax, title in (
            (ax_top_f, rf"Pier top {stem} (node 5)"),
            (ax_bot_f, rf"Pier base {stem} (node 1)"),
        ):
            ax.set_title(title)
            ax.set_ylabel(ylab)
            ax.grid(True, ls=":", alpha=0.45)
            ax.legend(fontsize=8, loc="best")
        if d595 is not None:
            for ax in (ax_top_z, ax_bot_z):
                ax.set_xlim(d595[0], d595[1])
                ax.set_title("D5–95 zoom")
                ax.grid(True, ls=":", alpha=0.45)
        ax_bot_f.set_xlabel(r"$t$ (s)")
        ax_bot_z.set_xlabel(r"$t$ (s)")
        fig.suptitle(f"Baseline offline EQ — pier {stem} ({integ})", fontsize=11)
        path = out / f"hist_pier_{stem}_{integ}.png"
        fig.savefig(path, dpi=peq.DPI)
        plt.close(fig)
        print(f"PlotEQOfflineCompare: wrote {path}")

    # --- compare drift overlay ---
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.0), constrained_layout=True)
    ax_f, ax_z = axes
    n_ok = 0
    for struct in STRUCT_SLUGS:
        eq = case_dir(root, struct, integ)
        if not (eq / "window_meta.txt").is_file():
            continue
        top = load_pier_series(eq, 5, "disp")
        bot = load_pier_series(eq, 1, "disp")
        if top is None or bot is None:
            continue
        t, u_top, u_bot = _interp_align(top[0], top[1], bot[0], bot[1])
        y = peq.to_mm(u_top - u_bot)
        ax_f.plot(t, y, color=COLORS[struct], lw=1.1, label=LABELS[struct])
        ax_z.plot(t, y, color=COLORS[struct], lw=1.1, label=LABELS[struct])
        n_ok += 1
    if n_ok:
        for ax in (ax_f, ax_z):
            ax.axhline(0.0, color="#bbb", lw=0.8)
            ax.grid(True, ls=":", alpha=0.45)
            ax.legend(fontsize=8, loc="best")
        ax_f.set_ylabel(r"drift $\Delta u_x$ (mm)")
        ax_f.set_title(r"Pier UX drift $u_\mathrm{top}-u_\mathrm{bot}$")
        ax_f.set_xlabel(r"$t$ (s)")
        ax_z.set_xlabel(r"$t$ (s)")
        if d595 is not None:
            ax_z.set_xlim(d595[0], d595[1])
            ax_z.set_title("D5–95 zoom")
        fig.suptitle(f"Baseline offline EQ — pier UX drift ({integ})", fontsize=11)
        path = out / f"hist_drift_ux_{integ}.png"
        fig.savefig(path, dpi=peq.DPI)
        plt.close(fig)
        print(f"PlotEQOfflineCompare: wrote {path}")
    else:
        plt.close(fig)


def plot_hinge_compare(out: Path, root: Path, integ: str) -> None:
    """Hinge M–θ history and hysteresis for base (and top when lumped)."""
    # Base: history + hyst
    fig_h, ax_mh = plt.subplots(figsize=(9.5, 4.2), constrained_layout=True)
    fig_y, ax_hy = plt.subplots(figsize=(6.4, 5.6), constrained_layout=True)
    for struct in STRUCT_SLUGS:
        eq = case_dir(root, struct, integ)
        meta_p = eq / "window_meta.txt"
        if not meta_p.is_file():
            continue
        meta = peq.read_meta(eq)
        h = load_hinge(eq, meta, top=False)
        if h is None:
            continue
        t, M, th, kind = h
        col = COLORS[struct]
        lab = f"{LABELS[struct]} ({kind})"
        ax_mh.plot(t, M / 1e3, color=col, lw=1.0, label=lab)
        ax_hy.plot(th, M / 1e3, color=col, lw=0.9, label=lab)
    ax_mh.set_xlabel(r"$t$ (s)")
    ax_mh.set_ylabel(r"$M_z$ (kN·m)")
    ax_mh.set_title(f"Pier base hinge moment ({integ})")
    ax_mh.grid(True, ls=":", alpha=0.45)
    ax_mh.legend(fontsize=8)
    p1 = out / f"hist_hinge_base_Mz_{integ}.png"
    fig_h.savefig(p1, dpi=peq.DPI)
    plt.close(fig_h)

    ax_hy.set_xlabel(r"$\theta$ (rad)  [FBC: $\kappa\cdot L_p$]")
    ax_hy.set_ylabel(r"$M_z$ (kN·m)")
    ax_hy.set_title(f"Pier base hinge hysteresis ({integ})")
    ax_hy.grid(True, ls=":", alpha=0.45)
    ax_hy.legend(fontsize=8)
    p2 = out / f"hyst_hinge_base_{integ}.png"
    fig_y.savefig(p2, dpi=peq.DPI)
    plt.close(fig_y)
    print(f"PlotEQOfflineCompare: wrote {p1}  {p2}")

    # Top hinges (lumped only)
    fig_h, ax_mh = plt.subplots(figsize=(9.5, 4.2), constrained_layout=True)
    fig_y, ax_hy = plt.subplots(figsize=(6.4, 5.6), constrained_layout=True)
    n_top = 0
    for struct in ("lumped_hold", "lumped_free"):
        eq = case_dir(root, struct, integ)
        if not (eq / "window_meta.txt").is_file():
            continue
        meta = peq.read_meta(eq)
        h = load_hinge(eq, meta, top=True)
        if h is None:
            continue
        t, M, th, kind = h
        col = COLORS[struct]
        lab = LABELS[struct]
        ax_mh.plot(t, M / 1e3, color=col, lw=1.0, label=lab)
        ax_hy.plot(th, M / 1e3, color=col, lw=0.9, label=lab)
        n_top += 1
    if n_top:
        ax_mh.set_xlabel(r"$t$ (s)")
        ax_mh.set_ylabel(r"$M_z$ (kN·m)")
        ax_mh.set_title(f"Pier top hinge moment ({integ})")
        ax_mh.grid(True, ls=":", alpha=0.45)
        ax_mh.legend(fontsize=8)
        p3 = out / f"hist_hinge_top_Mz_{integ}.png"
        fig_h.savefig(p3, dpi=peq.DPI)
        ax_hy.set_xlabel(r"$\theta$ (rad)")
        ax_hy.set_ylabel(r"$M_z$ (kN·m)")
        ax_hy.set_title(f"Pier top hinge hysteresis ({integ})")
        ax_hy.grid(True, ls=":", alpha=0.45)
        ax_hy.legend(fontsize=8)
        p4 = out / f"hyst_hinge_top_{integ}.png"
        fig_y.savefig(p4, dpi=peq.DPI)
        print(f"PlotEQOfflineCompare: wrote {p3}  {p4}")
    plt.close(fig_h)
    plt.close(fig_y)

    # Rotation histories (base)
    fig, ax = plt.subplots(figsize=(9.5, 4.2), constrained_layout=True)
    for struct in STRUCT_SLUGS:
        eq = case_dir(root, struct, integ)
        if not (eq / "window_meta.txt").is_file():
            continue
        meta = peq.read_meta(eq)
        h = load_hinge(eq, meta, top=False)
        if h is None:
            continue
        t, _M, th, kind = h
        ax.plot(t, th, color=COLORS[struct], lw=1.0, label=f"{LABELS[struct]} ({kind})")
    ax.set_xlabel(r"$t$ (s)")
    ax.set_ylabel(r"$\theta$ (rad)")
    ax.set_title(f"Pier base hinge rotation ({integ}; FBC uses κ·Lp)")
    ax.grid(True, ls=":", alpha=0.45)
    ax.legend(fontsize=8)
    p5 = out / f"hist_hinge_base_theta_{integ}.png"
    fig.savefig(p5, dpi=peq.DPI)
    plt.close(fig)
    print(f"PlotEQOfflineCompare: wrote {p5}")


def write_snapshot(
    eq: Path,
    t_target: float,
    out_png: Path,
    title: str,
) -> bool:
    """One deformed-window + history frame at t_target (prototype s)."""
    meta = peq.read_meta(eq)
    tags, xy = peq.read_nodes(eq)
    disp_tags = peq.read_disp_nodes(eq)
    lines, quads = peq.read_eles(eq)
    disp_files = meta.get("dispFiles", "window_disp.out").split()
    t, ux, uy = peq.load_window_disp(eq, disp_tags, disp_files)
    if t.size == 0:
        return False
    if peq.SUBTRACT_T0:
        ux = peq.maybe_t0(ux)
        uy = peq.maybe_t0(uy)
    idx = {tg: i for i, tg in enumerate(disp_tags)}
    js = peq.load_spring_json(meta)
    traces = peq.pick_frame_traces(disp_tags, xy, meta, idx)
    if not traces:
        return False
    k = int(np.argmin(np.abs(t - t_target)))
    mesh = peq._frame_mesh_parts(disp_tags, xy, lines, quads, js, meta)
    amp = max(float(np.max(np.abs(ux))), float(np.max(np.abs(uy))), 1e-6)
    pad = peq.SCALE * amp + 0.5
    xlim = (float(mesh["X0"].min()) - pad, float(mesh["X0"].max()) + pad)
    ylim = (float(mesh["Y0"].min()) - pad, float(mesh["Y0"].max()) + pad)
    ctx = peq._create_frame_hist_figure(t, ux, mesh, traces, xlim, ylim, 1)
    peq._update_frame_hist_figure(ctx, k, 0, ux[k], uy[k])
    ctx["ttl"].set_text(f"{title}\n{peq.frame_time_title(float(t[k]), 0, 1)}")
    out_png.parent.mkdir(parents=True, exist_ok=True)
    ctx["fig"].savefig(out_png, dpi=120, facecolor="white")
    plt.close(ctx["fig"])
    print(f"PlotEQOfflineCompare: wrote {out_png}")
    return True


def plot_snapshots(out: Path, root: Path, integ: str) -> None:
    """Sync free-FBC peak + each case's own pier-top peak."""
    free = case_dir(root, "fbc_free", integ)
    if not (free / "window_meta.txt").is_file():
        print(f"PlotEQOfflineCompare: no fbc_free_{integ} for sync snapshots")
        return
    top = load_pier_ux(free, 5)
    if top is None:
        print("PlotEQOfflineCompare: no pier_node_5 on fbc_free")
        return
    t_sync, ux_sync, _ = peak_abs_ux(*top)
    sync_dir = out / f"snap_sync_fbc_free_peak_{integ}"
    own_dir = out / f"snap_own_peak_{integ}"
    note = sync_dir / "peak_time.txt"
    note.parent.mkdir(parents=True, exist_ok=True)
    note.write_text(
        f"integrator {integ}\n"
        f"ref_case fbc_free_{integ}\n"
        f"t_peak_s {t_sync:.6g}\n"
        f"ux_peak_m {ux_sync:.6g}\n"
        f"note sync snapshots use this prototype time for all cases\n",
        encoding="utf-8",
    )

    peaks = []
    for struct in STRUCT_SLUGS:
        eq = case_dir(root, struct, integ)
        if not (eq / "window_meta.txt").is_file():
            continue
        top_c = load_pier_ux(eq, 5)
        if top_c is None:
            continue
        t_pk, ux_pk, _ = peak_abs_ux(*top_c)
        peaks.append((struct, t_pk, ux_pk))
        write_snapshot(
            eq,
            t_sync,
            sync_dir / f"{struct}.png",
            f"{LABELS[struct]} @ free-FBC peak t={t_sync:.2f}s",
        )
        write_snapshot(
            eq,
            t_pk,
            own_dir / f"{struct}.png",
            f"{LABELS[struct]} @ own peak t={t_pk:.2f}s",
        )

    own_note = own_dir / "peak_times.txt"
    own_dir.mkdir(parents=True, exist_ok=True)
    lines = [f"integrator {integ}", "struct t_peak_s ux_peak_m"]
    lines.extend(f"{s} {tp:.6g} {ux:.6g}" for s, tp, ux in peaks)
    own_note.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot_integrator_pair(
    out: Path,
    root: Path,
    struct: str,
    integ_a: str,
    integ_b: str,
    d595: tuple[float, float] | None,
) -> None:
    """Overlay two integrators for one structural case (incomplete dumps OK)."""
    eq_a = case_dir(root, struct, integ_a)
    eq_b = case_dir(root, struct, integ_b)
    if not (eq_a / "pier_node_5.out").is_file() or not (eq_b / "pier_node_5.out").is_file():
        print(f"PlotEQOfflineCompare: missing pier dumps for {struct} {integ_a}/{integ_b}")
        return
    dest = out / f"{struct}_{integ_a}_vs_{integ_b}"
    dest.mkdir(parents=True, exist_ok=True)
    styles = {
        integ_a: dict(color="#1f77b4", lw=1.3),
        integ_b: dict(color="#d62728", lw=1.1, ls="--"),
    }
    labels = {integ_a: integ_a, integ_b: integ_b}

    def _series(eq: Path, tag: int, kind: str):
        s = load_pier_series(eq, tag, kind)
        return s

    # UX top / bot / drift
    fig, axes = plt.subplots(3, 2, figsize=(11.2, 9.0), constrained_layout=True)
    rows = (
        (5, "disp", r"Pier top $\Delta u_x$ (mm)", True),
        (1, "disp", r"Pier base $\Delta u_x$ (mm)", True),
        (None, "drift", r"Drift $u_\mathrm{top}-u_\mathrm{bot}$ (mm)", True),
    )
    t_last = {}
    for i_row, (tag, kind, ylab, to_mm) in enumerate(rows):
        ax_f, ax_z = axes[i_row]
        for integ, eq in ((integ_a, eq_a), (integ_b, eq_b)):
            if kind == "drift":
                top = _series(eq, 5, "disp")
                bot = _series(eq, 1, "disp")
                if top is None or bot is None:
                    continue
                t, u_top, u_bot = _interp_align(top[0], top[1], bot[0], bot[1])
                y = peq.to_mm(u_top - u_bot) if to_mm else (u_top - u_bot)
            else:
                ser = _series(eq, tag, kind)
                if ser is None:
                    continue
                t, x = ser
                y = peq.to_mm(x) if to_mm else x
            t_last[integ] = float(t[-1])
            ax_f.plot(t, y, label=labels[integ], **styles[integ])
            ax_z.plot(t, y, label=labels[integ], **styles[integ])
        ax_f.set_ylabel(ylab)
        ax_f.grid(True, ls=":", alpha=0.45)
        ax_z.grid(True, ls=":", alpha=0.45)
        if d595 is not None:
            ax_z.set_xlim(d595[0], min(d595[1], max(t_last.values()) if t_last else d595[1]))
        if i_row == 0:
            ax_f.set_title(f"{LABELS.get(struct, struct)}")
            ax_z.set_title("D5–95 zoom (to last overlapping sample)")
            ax_f.legend(fontsize=8)
    axes[-1, 0].set_xlabel(r"$t$ (s)")
    axes[-1, 1].set_xlabel(r"$t$ (s)")
    note = ", ".join(f"{k} to {v:.1f} s" for k, v in t_last.items())
    fig.suptitle(f"{struct} UX  ({note})", fontsize=11)
    path = dest / "hist_ux.png"
    fig.savefig(path, dpi=peq.DPI)
    plt.close(fig)
    print(f"PlotEQOfflineCompare: wrote {path}")

    # VX / AX top
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 6.4), constrained_layout=True)
    for i_row, (kind, ylab) in enumerate(
        (("vel", r"Pier top $v_x$ (m/s)"), ("acc", r"Pier top $a_x$ (m/s$^2$)"))
    ):
        ax_f, ax_z = axes[i_row]
        for integ, eq in ((integ_a, eq_a), (integ_b, eq_b)):
            ser = _series(eq, 5, kind)
            if ser is None:
                continue
            t, x = ser
            ax_f.plot(t, x, label=labels[integ], **styles[integ])
            ax_z.plot(t, x, label=labels[integ], **styles[integ])
        ax_f.set_ylabel(ylab)
        ax_f.grid(True, ls=":", alpha=0.45)
        ax_z.grid(True, ls=":", alpha=0.45)
        ax_f.legend(fontsize=8)
        if d595 is not None and t_last:
            ax_z.set_xlim(d595[0], min(d595[1], max(t_last.values())))
    axes[0, 0].set_title("full")
    axes[0, 1].set_title("D5–95 zoom")
    axes[-1, 0].set_xlabel(r"$t$ (s)")
    axes[-1, 1].set_xlabel(r"$t$ (s)")
    fig.suptitle(f"{struct} top VX / AX", fontsize=11)
    path = dest / "hist_vx_ax.png"
    fig.savefig(path, dpi=peq.DPI)
    plt.close(fig)
    print(f"PlotEQOfflineCompare: wrote {path}")

    # Hinges if both exist
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.0), constrained_layout=True)
    n_h = 0
    for integ, eq in ((integ_a, eq_a), (integ_b, eq_b)):
        if not (eq / "window_meta.txt").is_file():
            continue
        meta = peq.read_meta(eq)
        h = load_hinge(eq, meta, top=False)
        if h is None:
            continue
        t, M, th, kind = h
        n_h += 1
        axes[0, 0].plot(t, M / 1e3, label=f"{integ} ({kind})", **styles[integ])
        axes[0, 1].plot(th, M / 1e3, label=f"{integ} ({kind})", **styles[integ])
        ht = load_hinge(eq, meta, top=True)
        if ht is not None:
            tt, Mt, tht, _ = ht
            axes[1, 0].plot(tt, Mt / 1e3, label=integ, **styles[integ])
            axes[1, 1].plot(tht, Mt / 1e3, label=integ, **styles[integ])
    if n_h:
        axes[0, 0].set_title("base Mz vs t")
        axes[0, 1].set_title(r"base hyst $M$–$\theta$")
        axes[1, 0].set_title("top Mz vs t")
        axes[1, 1].set_title(r"top hyst $M$–$\theta$")
        for ax in axes.ravel():
            ax.grid(True, ls=":", alpha=0.45)
            ax.legend(fontsize=8)
        axes[0, 0].set_ylabel(r"$M_z$ (kN·m)")
        axes[1, 0].set_ylabel(r"$M_z$ (kN·m)")
        axes[1, 0].set_xlabel(r"$t$ (s)")
        axes[1, 1].set_xlabel(r"$\theta$ (rad)")
        path = dest / "hist_hinge.png"
        fig.savefig(path, dpi=peq.DPI)
        plt.close(fig)
        print(f"PlotEQOfflineCompare: wrote {path}")
    else:
        plt.close(fig)
        print("PlotEQOfflineCompare: no hinge files for pair (killed dump?)")


def write_index(out: Path, root: Path, integs: list[str]) -> None:
    """Short index of per-case PlotEQ spring/pile PNGs."""
    lines = [
        "# Offline baseline EQ compare",
        f"root: {root}",
        "",
        "Per-case PlotEQ plots (springs, piles, window hist) live in:",
    ]
    for integ in integs:
        for struct in STRUCT_SLUGS:
            eq = case_dir(root, struct, integ)
            plots = eq / "plots"
            lines.append(f"  {plots.relative_to(REPO) if plots.is_dir() else eq.name + ' (missing)'}")
            if plots.is_dir():
                for name in sorted(
                    p.name
                    for p in plots.glob("*.png")
                    if any(
                        k in p.name
                        for k in (
                            "spring",
                            "pile",
                            "hinge",
                            "hist_ux",
                            "envelope",
                        )
                    )
                ):
                    lines.append(f"    - {name}")
    path = out / "INDEX.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"PlotEQOfflineCompare: wrote {path}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    ap.add_argument(
        "--integrators",
        nargs="+",
        default=("mkr",),
        choices=("newmark", "mkr"),
    )
    ap.add_argument(
        "--pair",
        nargs=3,
        metavar=("STRUCT", "INTEG_A", "INTEG_B"),
        default=None,
        help="Overlay two integrators for one case, e.g. lumped_hold mkr newmark",
    )
    args = ap.parse_args()
    root = args.root.resolve()
    out = root / "compare"
    out.mkdir(parents=True, exist_ok=True)

    # Safety: stay under eq_offline
    try:
        root.relative_to(DEFAULT_ROOT.resolve())
    except ValueError:
        if root != DEFAULT_ROOT.resolve():
            # allow explicit root only if still named eq_offline leaf
            if root.name != "eq_offline":
                raise SystemExit(f"refusing compare root outside eq_offline: {root}")

    d595 = d595_proto_window(0.0)
    if args.pair is not None:
        struct, a, b = args.pair
        plot_integrator_pair(out, root, struct, a, b, d595)
        return 0
    for integ in args.integrators:
        present = [
            s
            for s in STRUCT_SLUGS
            if (case_dir(root, s, integ) / "window_meta.txt").is_file()
        ]
        if not present:
            print(f"PlotEQOfflineCompare: no dumps for {integ}")
            continue
        print(f"PlotEQOfflineCompare: {integ} cases={present}")
        plot_pier_histories(out, root, integ, d595)
        plot_hinge_compare(out, root, integ)
        plot_snapshots(out, root, integ)
    write_index(out, root, list(args.integrators))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
