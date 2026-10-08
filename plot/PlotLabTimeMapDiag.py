#!/usr/bin/env python3
"""
Goals
-----
Lab-time map diagnostics for RTHS. Map from the start: ctrlDisp
recorder \(t_{i+1}\) ↔ lab atTarget; \(f=t_{\mathrm{lab}}-t_{\mathrm{OS}}/\sqrt{\lambda}\):

  diag_lab_time_map.png   before/after, map curve, f
  diag_offset_align.png   ctrlDisp vs tar/com/mea (start + peak)

  python plot/PlotLabTimeMapDiag.py
  python plot/PlotLabTimeMapDiag.py W02 W04
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

from PlotActuatorForce import (
    COLOR_FRC,
    TIME_SCALE_FROUDE,
    mat_dump_for_test,
    run_title,
    scale_paper_fonts,
)
from PlotEQComparePairs import load_type_conv3
from PlotEQCompareRuns import (
    apply_paper_style,
    load_mat_mea_feedback,
    model_disp_to_proto_mm,
)
from gm_duration import d595_proto_window, gm_start_time_s
from lab_paths import (
    LOCAL_OPENSEES_DATA,
    MAT_EXTRACT_DIR,
    M_TO_MM,
    resolve_opensees_data,
    test_os_plots_dir,
)
from lab_time_map import (
    handshake_map_for_dump,
    lab_times_to_os_handshake,
    map_os_times_handshake,
    slowdown_lab_onsets,
)

N_COUNT = 10
SLOWDOWN_STATE = 2
T_EARLY_S = 1.0
COLOR_LAB = COLOR_FRC
COLOR_MEA = COLOR_FRC
COLOR_COM = "#2a9d8f"
COLOR_TAR = "#6a1b9a"
COLOR_OS = "#c62828"
COLOR_MAP = "#2e7d32"
COLOR_REF = "#666666"
COLOR_K = "#6a1b9a"
COLOR_SLOW = "#FFC04D"

# Half-width (lab s) for motion panels centered on the handshake peak.
PEAK_ZOOM_HALF_S = 7.5


def lab_clock_parts(mat_name: str) -> tuple[float, float, float, float, float]:
    z = np.load(MAT_EXTRACT_DIR / f"{Path(mat_name).stem}.npz", allow_pickle=True)
    t = np.asarray(z["stateOS_time"], dtype=float)
    tc3 = np.asarray(z["stateOS_data"][:, 0], dtype=int)
    dt_con = float(np.median(np.diff(t)))
    dt_sim = float(N_COUNT) * dt_con
    land = np.concatenate([[False], (tc3[:-1] == 1) & (tc3[1:] == 0)])
    t_land = float(t[int(np.where(land)[0][0])])
    prev = np.concatenate([[int(tc3[0])], tc3[:-1]])
    t_ex = float(t[np.where((tc3 == 1) & (prev == 0) & (t > t_land))[0][0]])
    t_goal = t_ex - dt_con
    k = t_goal - dt_sim
    return k, t_land, t_ex, t_goal, dt_sim


def midrun_cum_slow(t_lab, tc3, t_ex):
    is_slow = (tc3 == SLOWDOWN_STATE).astype(float)
    d_slow = np.zeros_like(t_lab)
    d_slow[1:] = is_slow[:-1] * np.diff(t_lab)
    d_slow[t_lab <= t_ex] = 0.0
    cum = np.cumsum(d_slow)
    prev = np.concatenate([[int(tc3[0])], tc3[:-1]])
    ons = t_lab[(tc3 == SLOWDOWN_STATE) & (prev != SLOWDOWN_STATE) & (t_lab > t_ex)]
    return cum, np.asarray(ons, dtype=float)


def map_os_with_slow_f(t_os, k, t_lab, cum_mid):
    target = k + t_os / TIME_SCALE_FROUDE
    g = np.asarray(t_lab - cum_mid, dtype=float)
    for i in range(1, g.size):
        if g[i] < g[i - 1]:
            g[i] = g[i - 1]
    idx = np.clip(np.searchsorted(g, target, side="left"), 0, t_lab.size - 1)
    f = np.maximum(np.asarray(cum_mid[idx], dtype=float), 0.0)
    return target + f, f


def lab_times_to_os(t_lab_events, t_os, t_mapped):
    if t_lab_events.size == 0:
        return np.asarray([], dtype=float)
    tm = np.array(t_mapped, dtype=float, copy=True)
    for i in range(1, tm.size):
        if tm[i] < tm[i - 1]:
            tm[i] = tm[i - 1]
    idx = np.clip(np.searchsorted(tm, t_lab_events, side="left"), 0, t_os.size - 1)
    return np.asarray(t_os[idx], dtype=float)


def load_os_block(mat_name: str, key: str):
    npz = MAT_EXTRACT_DIR / f"{Path(mat_name).stem}.npz"
    z = np.load(npz, allow_pickle=True)
    if f"{key}_time" in z.files and f"{key}_primary" in z.files:
        return np.asarray(z[f"{key}_time"], dtype=float), np.asarray(
            z[f"{key}_primary"], dtype=float
        )
    data = np.asarray(z[f"{key}_data"], dtype=float)
    names = [str(x) for x in z[f"{key}_signalNames"].tolist()]
    time_col = next(i for i, n in enumerate(names) if n.strip().lower() == "time")
    signal_col = next(i for i in range(data.shape[1]) if i != time_col)
    return data[:, time_col], data[:, signal_col]


def write_map(test_id: str) -> int:
    pair = mat_dump_for_test(test_id)
    if pair is None:
        return 1
    mat, dump = pair
    root = resolve_opensees_data() or LOCAL_OPENSEES_DATA
    dump_path = root / dump
    mea = load_mat_mea_feedback(mat)
    tar = load_os_block(mat, "tarSigOS")
    if mea is None or tar is None:
        return 1
    t_lab_mea, u_model = mea
    u_lab = model_disp_to_proto_mm(u_model)
    t_tar, u_tar = tar[0], model_disp_to_proto_mm(tar[1])
    d = np.loadtxt(dump_path / "Elmt101_ctrlDsp.out", ndmin=2)
    t_os = np.asarray(d[:, 0], dtype=float)
    u_os = (np.asarray(d[:, 1], dtype=float) - d[0, 1]) * M_TO_MM
    t_naive = t_os / TIME_SCALE_FROUDE

    try:
        hs = handshake_map_for_dump(dump_path, mat, t_tar, u_tar)
    except RuntimeError as exc:
        print(f"PlotLabTimeMapDiag: {test_id} handshake failed: {exc}", file=sys.stderr)
        return 1
    t_mapped = map_os_times_handshake(t_os, hs)
    f_nl = t_mapped - t_naive

    t0 = gm_start_time_s(dump_path)
    d595 = d595_proto_window(t0)
    if d595 is not None:
        z0, z1 = d595[0] / TIME_SCALE_FROUDE, d595[1] / TIME_SCALE_FROUDE
    else:
        z0, z1 = 36.0, 90.0
    t_peak_lab = float(hs.t_lab[hs.i_peak])
    z1 = max(z1, t_peak_lab + 5.0)

    apply_paper_style()
    scale_paper_fonts(1.35)
    fig = plt.figure(figsize=(11.4, 9.2), layout="constrained")
    gs = GridSpec(3, 2, figure=fig, height_ratios=[1.0, 1.0, 1.1])
    ax_b = fig.add_subplot(gs[0, :])
    ax_a = fig.add_subplot(gs[1, :])
    ax_m = fig.add_subplot(gs[2, 0])
    ax_f = fig.add_subplot(gs[2, 1])

    step_lab = max(1, len(t_lab_mea) // 8000)
    step_os = max(1, len(t_os) // 8000)
    step_m = max(1, len(t_os) // 12000)
    step_tar = max(1, len(t_tar) // 8000)

    ax_b.plot(
        t_tar[::step_tar],
        u_tar[::step_tar],
        color=COLOR_TAR,
        lw=0.9,
        ls="--",
        label="tarSigOS (lab $t$)",
        zorder=3,
    )
    ax_b.plot(
        t_lab_mea[::step_lab],
        u_lab[::step_lab],
        color=COLOR_LAB,
        lw=0.85,
        alpha=0.75,
        label="meaSigOS (lab $t$)",
        zorder=2,
    )
    ax_b.plot(
        t_naive[::step_os],
        u_os[::step_os],
        color=COLOR_OS,
        lw=0.9,
        label=r"ctrlDisp ($t_{\mathrm{int}}/\sqrt{\lambda}$)",
        zorder=4,
    )
    ax_b.set_xlim(z0, z1)
    ax_b.set_ylabel(r"$\Delta u$ (mm, prototype)")
    ax_b.set_title(r"before: naive Froude ($t_{\mathrm{int}}/\sqrt{\lambda}$)")
    ax_b.legend(loc="upper right", frameon=True, fancybox=False, edgecolor="#333")
    ax_b.grid(True, ls=":", alpha=0.45)

    ax_a.plot(
        t_tar[::step_tar],
        u_tar[::step_tar],
        color=COLOR_TAR,
        lw=0.9,
        ls="--",
        label="tarSigOS (lab $t$)",
        zorder=3,
    )
    ax_a.plot(
        t_lab_mea[::step_lab],
        u_lab[::step_lab],
        color=COLOR_LAB,
        lw=0.85,
        alpha=0.75,
        label="meaSigOS (lab $t$)",
        zorder=2,
    )
    ax_a.plot(
        t_mapped[::step_os],
        u_os[::step_os],
        color=COLOR_OS,
        lw=0.9,
        label=r"ctrlDisp on $t=$atTarget$(t_{\mathrm{int}})$",
        zorder=4,
    )
    ax_a.set_xlim(z0, z1)
    ax_a.set_xlabel(r"$t$ (s) lab / model scale")
    ax_a.set_ylabel(r"$\Delta u$ (mm, prototype)")
    ax_a.set_title(
        rf"after: $t_{{\mathrm{{int}}}}\leftrightarrow$atTarget  ·  "
        rf"$f_{{\mathrm{{end}}}}={float(hs.f_at[-1])*1e3:.1f}\,\mathrm{{ms}}$"
    )
    ax_a.legend(loc="upper right", frameon=True, fancybox=False, edgecolor="#333")
    ax_a.grid(True, ls=":", alpha=0.45)

    ax_m.plot(
        t_os[::step_m],
        t_naive[::step_m],
        color="#bdbdbd",
        lw=1.0,
        ls=":",
        label=r"$t_{\mathrm{int}}/\sqrt{\lambda}$",
    )
    ax_m.plot(
        t_os[::step_m],
        t_mapped[::step_m],
        color=COLOR_MAP,
        lw=1.05,
        label=r"$t=$atTarget$(t_{\mathrm{int}})$",
    )
    ax_m.set_xlabel(r"$t_{\mathrm{int}}$ (s, prototype)")
    ax_m.set_ylabel(r"$t$ (s, lab / model)")
    ax_m.set_title(r"map: recorder $t_{i+1}$ $\leftrightarrow$ atTarget")
    ax_m.legend(loc="lower right", frameon=True, fancybox=False, edgecolor="#333")
    ax_m.grid(True, ls=":", alpha=0.45)
    ax_m.set_xlim(0.0, float(t_os[-1]))
    ax_m.set_ylim(0.0, max(float(t_mapped[-1]), float(t_naive[-1])) * 1.02)

    ax_f.plot(t_os[::step_m], f_nl[::step_m] * 1e3, color=COLOR_MAP, lw=1.2)
    ax_f.axhline(0.0, color=COLOR_REF, lw=0.7, ls="--")
    # typeConv3==2 onsets (lab) → OS t via inverse atTarget map
    ons_lab = slowdown_lab_onsets(mat)
    n_slow = len(ons_lab)
    if n_slow:
        t_slow_os = lab_times_to_os_handshake(np.asarray(ons_lab, dtype=float), hs)
        step_mark = max(1, n_slow // 40)
        for tos in np.asarray(t_slow_os, dtype=float)[::step_mark]:
            if 0.0 <= float(tos) <= float(t_os[-1]):
                ax_f.axvline(float(tos), color=COLOR_SLOW, alpha=0.35, lw=0.7, zorder=0)
    ax_f.set_xlabel(r"$t_{\mathrm{int}}$ (s, prototype)")
    ax_f.set_ylabel(r"$f(t_{\mathrm{int}})$ (ms model)")
    ax_f.set_title(
        rf"$f=t-t_{{\mathrm{{int}}}}/\sqrt{{\lambda}}$  "
        rf"($f_{{\mathrm{{end}}}}={float(f_nl[-1])*1e3:.1f}\,\mathrm{{ms}}$; "
        rf"amber$=$slowdown onset, $n={n_slow}$)"
    )
    y_lo = min(0.0, float(np.min(f_nl) * 1e3))
    y_hi = max(10.0, float(np.max(f_nl) * 1e3) * 1.15)
    ax_f.set_ylim(y_lo - 0.05 * (y_hi - y_lo + 1e-9), y_hi)
    ax_f.grid(True, ls=":", alpha=0.45)

    fig.suptitle(
        rf"{run_title(test_id)}  ·  "
        rf"atTarget map (from start)  ·  "
        rf"$n_{{\mathrm{{pair}}}}={hs.t_os.size}$",
        fontsize=plt.rcParams["axes.labelsize"],
        y=1.01,
    )
    out = test_os_plots_dir(test_id) / "diag_lab_time_map.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(
        f"PlotLabTimeMapDiag: wrote {out}  "
        f"f_end={float(f_nl[-1])*1e3:.2f} ms  n_pair={hs.t_os.size}"
    )
    return 0


def annotate_k_slip(ax, t_ref, u_ref, t_naive, u_os, k, m_ref, m_n):
    if not (np.any(m_ref) and np.any(m_n)):
        return
    i_pk = int(np.argmax(np.abs(u_ref[m_ref])))
    t_pk = float(t_ref[m_ref][i_pk])
    u_pk = float(u_ref[m_ref][i_pk])
    sign = 1.0 if u_pk >= 0 else -1.0
    u_n = u_os[m_n]
    t_n = t_naive[m_n]
    cand = np.where(sign * u_n > 0.5 * max(abs(u_pk), 1e-9))[0]
    if cand.size:
        j = int(cand[np.argmin(np.abs(t_n[cand] - t_pk))])
    else:
        j = int(np.argmin(np.abs(t_n - (t_pk - k))))
    t_os_pk = float(t_n[j])
    u_os_pk = float(u_n[j])
    y0 = min(float(np.nanmin(u_ref[m_ref])), float(np.nanmin(u_n)))
    y1 = max(float(np.nanmax(u_ref[m_ref])), float(np.nanmax(u_n)))
    y_arr = max(u_pk, u_os_pk) + 0.12 * (y1 - y0 + 1e-9)
    ax.annotate(
        "",
        xy=(t_pk, y_arr),
        xytext=(t_os_pk, y_arr),
        arrowprops=dict(arrowstyle="<->", color=COLOR_K, lw=1.6),
    )
    ax.text(
        0.5 * (t_pk + t_os_pk),
        y_arr,
        rf"$k={k*1e3:.1f}\,\mathrm{{ms}}$",
        color=COLOR_K,
        ha="center",
        va="bottom",
        fontsize=10,
        fontweight="bold",
    )
    pad = 0.18 * (y1 - y0 + 1e-9)
    ax.set_ylim(y0 - 0.08 * (y1 - y0 + 1e-9), y_arr + pad)


def write_align(test_id: str) -> int:
    pair = mat_dump_for_test(test_id)
    if pair is None:
        return 1
    mat, dump = pair
    root = resolve_opensees_data() or LOCAL_OPENSEES_DATA
    dump_path = root / dump
    tar = load_os_block(mat, "tarSigOS")
    com = load_os_block(mat, "comSigOS")
    mea = load_os_block(mat, "meaSigOS")
    t_tar, u_tar = tar[0], model_disp_to_proto_mm(tar[1])
    t_com, u_com = com[0], model_disp_to_proto_mm(com[1])
    t_mea, u_mea = mea[0], model_disp_to_proto_mm(mea[1])
    d = np.loadtxt(dump_path / "Elmt101_ctrlDsp.out", ndmin=2)
    t_os = np.asarray(d[:, 0], dtype=float)
    u_os = (np.asarray(d[:, 1], dtype=float) - d[0, 1]) * M_TO_MM
    try:
        hs = handshake_map_for_dump(dump_path, mat, t_tar, u_tar)
    except RuntimeError as exc:
        print(f"PlotLabTimeMapDiag: {test_id} handshake failed: {exc}", file=sys.stderr)
        return 1
    t_mapped = map_os_times_handshake(t_os, hs)
    t_naive = t_os / TIME_SCALE_FROUDE

    # Start zoom: span naive vs mapped first motion; after = lab [0, T_EARLY].
    t_start_lab = float(hs.t_lab[0])
    t_start_naive = float(hs.t_os[0]) / TIME_SCALE_FROUDE
    z0_start_b = 0.0
    z1_start_b = max(T_EARLY_S, t_start_lab + 0.2, t_start_naive + 0.2)
    z0_start_a = 0.0
    z1_start_a = T_EARLY_S

    t_peak_lab = float(hs.t_lab[hs.i_peak])
    t_peak_naive = float(hs.t_os[hs.i_peak]) / TIME_SCALE_FROUDE
    # Same xlim for peak before/after (span both clocks so lag is visible).
    z0_peak = max(0.0, min(t_peak_lab, t_peak_naive) - PEAK_ZOOM_HALF_S)
    z1_peak = max(t_peak_lab, t_peak_naive) + PEAK_ZOOM_HALF_S

    apply_paper_style()
    scale_paper_fonts(1.25)
    fig, axes = plt.subplots(4, 1, figsize=(11.2, 10.4), layout="constrained")
    ax_sb, ax_sa, ax_pb, ax_pa = axes
    step_lab = max(1, len(t_com) // 10000)
    step_os = max(1, len(t_os) // 12000)

    def draw(ax, z0, z1, t_os_plot, label_os):
        m_tar = (t_tar >= z0) & (t_tar <= z1)
        m_com = (t_com >= z0) & (t_com <= z1)
        m_mea = (t_mea >= z0) & (t_mea <= z1)
        m_os = (t_os_plot >= z0) & (t_os_plot <= z1)
        ax.plot(
            t_mea[m_mea][::step_lab],
            u_mea[m_mea][::step_lab],
            color=COLOR_MEA,
            lw=1.0,
            alpha=0.85,
            label="meaSigOS",
            zorder=2,
        )
        ax.plot(
            t_com[m_com][::step_lab],
            u_com[m_com][::step_lab],
            color=COLOR_COM,
            lw=1.15,
            label="comSigOS",
            zorder=3,
        )
        st = max(1, int(np.count_nonzero(m_tar) // 4000))
        ax.plot(
            t_tar[m_tar][::st],
            u_tar[m_tar][::st],
            color=COLOR_TAR,
            lw=1.0,
            ls="--",
            label="tarSigOS",
            zorder=4,
        )
        ax.plot(
            t_os_plot[m_os][::step_os],
            u_os[m_os][::step_os],
            color=COLOR_OS,
            lw=1.2,
            label=label_os,
            zorder=5,
        )
        ax.set_xlim(z0, z1)
        ax.set_ylabel(r"$\Delta u$ (mm, proto)")
        ax.legend(loc="best", fontsize=7.5, frameon=True, fancybox=False, edgecolor="#333")
        ax.grid(True, ls=":", alpha=0.45)
        u_all = np.concatenate(
            [u_mea[m_mea], u_com[m_com], u_tar[m_tar], u_os[m_os] if np.any(m_os) else [0.0]]
        )
        y0, y1 = float(np.nanmin(u_all)), float(np.nanmax(u_all))
        pad = 0.12 * (y1 - y0 + 1e-9)
        ax.set_ylim(y0 - pad, y1 + pad)

    draw(ax_sb, z0_start_b, z1_start_b, t_naive, r"ctrlDisp ($t_{\mathrm{int}}/\sqrt{\lambda}$)")
    ax_sb.set_title(
        rf"start before: naive Froude  "
        rf"(zoom $[{z0_start_b:.2f},{z1_start_b:.2f}]$ s)",
        loc="left",
    )
    draw(ax_sa, z0_start_a, z1_start_a, t_mapped, r"ctrlDisp on $t=$atTarget")
    ax_sa.set_title(
        rf"start after: atTarget map  "
        rf"(zoom $[{z0_start_a:.2f},{z1_start_a:.2f}]$ s · "
        rf"$f_0={float(hs.f_at[0])*1e3:.1f}\,\mathrm{{ms}}$)",
        loc="left",
    )
    ax_sa.set_xlabel(r"$t$ (s) lab / model")

    draw(ax_pb, z0_peak, z1_peak, t_naive, r"ctrlDisp ($t_{\mathrm{int}}/\sqrt{\lambda}$)")
    ax_pb.set_title(
        rf"peak before: naive Froude  "
        rf"(zoom $[{z0_peak:.1f},{z1_peak:.1f}]$ s · "
        rf"atTarget peak ${t_peak_lab:.1f}$ s, naive ${t_peak_naive:.1f}$ s)",
        loc="left",
    )
    draw(ax_pa, z0_peak, z1_peak, t_mapped, r"ctrlDisp on $t=$atTarget")
    ax_pa.set_title(
        rf"peak after: atTarget map  "
        rf"(zoom $[{z0_peak:.1f},{z1_peak:.1f}]$ s · peak ${t_peak_lab:.1f}$ s · "
        rf"$f={float(hs.f_at[hs.i_peak])*1e3:.1f}\,\mathrm{{ms}}$)",
        loc="left",
    )
    ax_pa.set_xlabel(r"$t$ (s) lab / model scale")

    fig.suptitle(
        rf"{run_title(test_id)}  ·  atTarget map: ctrlDisp vs tar / com / mea",
        fontsize=plt.rcParams["axes.labelsize"],
        y=1.01,
    )
    engine = fig.get_layout_engine()
    if engine is not None:
        engine.set(h_pad=0.03, hspace=0.07)

    out = test_os_plots_dir(test_id) / "diag_offset_align.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"PlotLabTimeMapDiag: wrote {out}")
    return 0


def main() -> int:
    if any(a in ("-h", "--help") for a in sys.argv[1:]):
        print(__doc__)
        return 0
    from PlotHydroSpectra import default_tests

    tests = [a for a in sys.argv[1:] if not a.startswith("-")] or default_tests()
    if not tests:
        print("PlotLabTimeMapDiag: no tests", file=sys.stderr)
        return 1
    rc = 0
    for tid in tests:
        rc = max(rc, write_map(tid), write_align(tid))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
