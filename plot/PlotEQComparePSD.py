#!/usr/bin/env python3
"""
Goals
-----
Pairwise FFT PSD overlays within each physical-model compare group.

Same interim reference as PlotEQComparePairs (fewest typeConv3→2 in GM
D5–95). One figure per non-ref dump:

  psd_ux_pair_<Test>.png

Five rows on a shared prototype-Hz axis (Hann FFT; lab \(f\) → proto via /√λ):
  tarSig / comSig / meaSig  — lab \(t\)
  pier \(u_x\) / ctrlDisp   — \(t_{\mathrm{int}}\)

Layout: FFT window history | PSD (all groups).

Window:
  Friday / campaign — GM D5–95 (OS proto; lab via atTarget map)
  Wed0819 — fair end from shortest OS progress
  \(T^\star=\min(t_{\mathrm{end}}^{\mathrm{OS}}-\mathrm{gmStart})\):
  pier / ctrlDisp FFT on OS \([0,\mathrm{gmStart}+T^\star]\);
  tar/com/mea on lab \([0,t_{\mathrm{lab}}^\star]\) with \(t_{\mathrm{lab}}^\star\)
  = atTarget map of that OS cut (ctrlDisp↔stateOS handshakes).

  python plot/PlotEQComparePSD.py
  python plot/PlotEQComparePSD.py <runDir> ...

Writes under OSU_SSI_BRIDGE_DATA_LOCAL/plots/compare/<Mesh>/<variant>/pairs/
(excluded runs → compare/_excluded/…). Ref = grey; other = navy.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from compare_groups import (
    analysis_skip_reason,
    dump_to_group,
    dump_to_row,
    dump_to_test_id,
    group_label,
    groups_by_dump,
    legend_labels_for_dumps,
)
from gm_duration import arias_significant_duration, d595_proto_window, gm_start_time_s
from lab_paths import compare_plots_dir, resolve_opensees_data
from lab_time_map import (
    LABEL_T_INT,
    LABEL_T_LAB,
    add_dual_time_xaxis_lab,
    handshake_map_for_dump,
    load_ctrl_disp_proto_mm,
    map_os_times_handshake,
)
from PlotEQComparePairs import (
    COLOR_OTHER,
    COLOR_REF,
    PAIR_DPI,
    pair_file_slug,
    pick_reference,
)
from PlotEQCompareRuns import apply_paper_style, resolve_run_path, run_duration, run_to_mat_name
from PlotHydroSpectra import (
    F_MIN_PROTO_HZ,
    _fft_proto_hz,
    _load_os_disp_proto_mm,
    add_dual_freq_xaxis,
    load_pier_ux_mm,
    mark_freq_guides,
    pier_f1_hz,
    pier_f4_hz,
    pier_f5_hz,
    ssi_guide_freqs,
)
from PlotActuatorForce import scale_paper_fonts

HELP = """\
usage: python plot/PlotEQComparePSD.py [runDir ...]

  no args   dumps from TestMatrix_lab_runs.csv (same groups as hist pairs)
  runDir    one or more dump folders
  output    LOCAL/plots/compare/<Mesh>/<variant>/pairs/psd_ux_pair_<Test>.png
  clocks    tar/com/mea = lab $t$; pier/ctrlDisp = $t_{\mathrm{int}}$
  layout    history (FFT window) | PSD
  fft       Hann periodogram, mean kept, no resample
  Fri       D5-95 on $t_{\mathrm{int}}$; lab cuts via atTarget map
  Wed       fair end: T*=min(t_end_OS-gmStart); pier/ctrl on $t_{\mathrm{int}}$;
            lab [0, atTarget(cut)] via handshake map
"""

FONT_SCALE = 1.15
# (key, row label, lab_clock) — OS-domain rows use fair/D5–95 $t_{\mathrm{int}}$ windows
ROW_SPECS = (
    ("tarSigOS", r"tarSig (lab $t$)", True),
    ("comSigOS", r"comSig (lab $t$)", True),
    ("meaSigOS", r"meaSig (lab $t$)", True),
    ("pier", r"pier $u_x$ ($t_{\mathrm{int}}$)", False),
    ("ctrlDisp", r"ctrlDisp ($t_{\mathrm{int}}$)", False),
)


def _series_for_dump(
    dump: Path,
    mat_name: str,
) -> dict[str, tuple[np.ndarray, np.ndarray] | None]:
    """
    Load tar/com/mea (lab t), pier ux and ctrlDisp (OS t) for one dump.

    Args:    dump  eqOutDir; mat_name  Simulink mat
    Returns: {key: (t, y_proto_mm) or None}
    """
    out: dict[str, tuple[np.ndarray, np.ndarray] | None] = {}
    for key in ("tarSigOS", "comSigOS", "meaSigOS"):
        out[key] = _load_os_disp_proto_mm(mat_name, key)
    out["pier"] = load_pier_ux_mm(dump)
    out["ctrlDisp"] = load_ctrl_disp_proto_mm(dump)
    return out


def _os_to_lab_cut(dump: Path, mat_name: str, t_os: float) -> float | None:
    """
    Lab atTarget time for an OpenSees recorder instant.

    Args:    dump; mat_name; t_os  proto s
    Returns: t_lab (s) or None
    """
    try:
        hs = handshake_map_for_dump(dump, mat_name)
        return float(np.asarray(map_os_times_handshake(t_os, hs)).reshape(-1)[0])
    except RuntimeError:
        return None


def _d595_windows(
    dump: Path,
    mat_name: str,
) -> tuple[tuple[float, float] | None, tuple[float, float] | None]:
    """
    D5–95 on OpenSees proto t and on lab Time (atTarget map).

    Args:    dump; mat_name
    Returns: (d595_os, d595_lab) each (t0, t1) or None
    """
    gm0 = gm_start_time_s(dump)
    d595_os = d595_proto_window(gm0)
    if d595_os is None:
        return None, None
    t0 = _os_to_lab_cut(dump, mat_name, float(d595_os[0]))
    t1 = _os_to_lab_cut(dump, mat_name, float(d595_os[1]))
    if t0 is None or t1 is None or t1 <= t0:
        return d595_os, None
    return d595_os, (t0, t1)


def _is_wed0819_group(slug: str) -> bool:
    """True for Wednesday 2026-08-19 compare folder."""
    return "Wed0819" in slug.replace("\\", "/")


def _os_end_progress(
    dump: Path,
    pier: tuple[np.ndarray, np.ndarray] | None,
) -> tuple[float, float] | None:
    """
    OpenSees stop time and post-gmStart progress.

    Args:    dump; pier  (t_os, ux) or None
    Returns: (t_end_os, T = t_end - gmStart) or None
    """
    if pier is None or pier[0].size < 64:
        return None
    t_end = float(pier[0][-1])
    gm0 = float(gm_start_time_s(dump))
    prog = t_end - gm0
    if not np.isfinite(prog) or prog <= 0.0:
        return None
    return t_end, prog


def _fair_windows(
    ref: Path,
    other: Path,
    mat_ref: str,
    mat_other: str,
    ser_r: dict[str, tuple[np.ndarray, np.ndarray] | None],
    ser_o: dict[str, tuple[np.ndarray, np.ndarray] | None],
) -> (
    tuple[
        tuple[float, float],
        tuple[float, float],
        tuple[float, float],
        tuple[float, float],
        float,
    ]
    | None
):
    """
    Fair FFT windows from shortest OS post-gmStart progress.

    Pier: [0, gmStart+T★] on OS t.
    Lab: [0, atTarget(gmStart+T★)] via handshake map.

    Args:    ref, other, mats, series dicts
    Returns: (win_lab_r, win_lab_o, win_os_r, win_os_o, T_star) or None
    """
    pr = _os_end_progress(ref, ser_r.get("pier"))
    po = _os_end_progress(other, ser_o.get("pier"))
    if pr is None or po is None:
        return None
    _end_r, prog_r = pr
    _end_o, prog_o = po
    t_star = float(min(prog_r, prog_o))
    pier_r = ser_r.get("pier")
    pier_o = ser_o.get("pier")
    if pier_r is None or pier_o is None:
        return None
    gm_r = float(gm_start_time_s(ref))
    gm_o = float(gm_start_time_s(other))
    t_cut_os_r = min(gm_r + t_star, float(pier_r[0][-1]))
    t_cut_os_o = min(gm_o + t_star, float(pier_o[0][-1]))
    t_cut_lab_r = _os_to_lab_cut(ref, mat_ref, t_cut_os_r)
    t_cut_lab_o = _os_to_lab_cut(other, mat_other, t_cut_os_o)
    if t_cut_lab_r is None or t_cut_lab_o is None:
        return None
    if t_cut_lab_r <= 0.0 or t_cut_lab_o <= 0.0:
        return None
    return (
        (0.0, t_cut_lab_r),
        (0.0, t_cut_lab_o),
        (0.0, t_cut_os_r),
        (0.0, t_cut_os_o),
        t_star,
    )


def plot_psd_pair(
    ref: Path,
    other: Path,
    out: Path,
    *,
    mat_ref: str,
    mat_other: str,
    labels: dict[str, str],
    test_ref: str,
    test_other: str,
    fair_lab_end: bool = False,
) -> bool:
    """
    Write one 5-row PSD pair (ref grey / other navy).

    Args:    ref, other dumps; mats; legend labels; Test IDs;
             fair_lab_end  Wed: OS progress -> lab [0, t_cut]; else D5-95
    Returns: True if a figure was written
    """
    ser_r = _series_for_dump(ref, mat_ref)
    ser_o = _series_for_dump(other, mat_other)
    if ser_r.get("meaSigOS") is None or ser_o.get("meaSigOS") is None:
        print(f"PlotEQComparePSD: skip {other.name} (need meaSigOS)", file=sys.stderr)
        return False

    if fair_lab_end:
        fair = _fair_windows(ref, other, mat_ref, mat_other, ser_r, ser_o)
        if fair is None:
            print(
                f"PlotEQComparePSD: skip {other.name} (fair lab end)",
                file=sys.stderr,
            )
            return False
        win_lab_r, win_lab_o, win_os_r, win_os_o, t_star = fair
        win_label = (
            rf"fair $[0,t_{{\mathrm{{cut}}}}]$ "
            rf"($T^\star$={t_star:.1f} s int; pier/ctrl on $t_{{\mathrm{{int}}}}$)"
        )
    else:
        win_os_r, win_lab_r = _d595_windows(ref, mat_ref)
        win_os_o, win_lab_o = _d595_windows(other, mat_other)
        win_label = "D5–95"
        if win_os_r is None or win_os_o is None:
            print(f"PlotEQComparePSD: skip {other.name} (no D5–95)", file=sys.stderr)
            return False

    f_pier = pier_f1_hz(test_other or test_ref)
    f_mode4 = pier_f4_hz(test_other or test_ref)
    f_mode5 = pier_f5_hz(test_other or test_ref)
    f_ssi = ssi_guide_freqs(test_other or test_ref)

    apply_paper_style()
    scale_paper_fonts(FONT_SCALE)
    label_fs = float(plt.rcParams["axes.labelsize"])
    plt.rcParams["axes.titlesize"] = label_fs
    plt.rcParams["figure.titlesize"] = label_fs * 1.15

    n_row = len(ROW_SPECS)
    # History (FFT window) | PSD for every group.
    n_col = 2
    fig_w = 12.8
    fig, axes = plt.subplots(
        n_row,
        n_col,
        figsize=(fig_w, 2.55 * n_row * (0.7 + 0.3 * FONT_SCALE)),
        sharex=False,
        layout="constrained",
        squeeze=False,
    )
    fig.set_constrained_layout_pads(h_pad=0.02, hspace=0.05, wspace=0.06)

    lab_ref = labels.get(ref.name, ref.name)
    lab_oth = labels.get(other.name, other.name)
    f_hi_all = 100.0
    any_curve = False

    for i, (key, row_name, lab_clock) in enumerate(ROW_SPECS):
        ax_psd = axes[i, 1]
        ax_hist = axes[i, 0]
        # tar/com/mea = lab $t$; pier/ctrlDisp = $t_{\mathrm{int}}$.
        use_lab = lab_clock
        series_r, series_o = ser_r.get(key), ser_o.get(key)
        if use_lab:
            win_r, win_o = win_lab_r, win_lab_o
        else:
            win_r, win_o = win_os_r, win_os_o
        if win_r is None or win_o is None:
            ax_psd.text(
                0.5, 0.5, "no window", transform=ax_psd.transAxes, ha="center"
            )
            ax_hist.text(
                0.5, 0.5, "no window", transform=ax_hist.transAxes, ha="center"
            )
            ax_psd.set_ylabel(row_name)
            ax_hist.set_ylabel(row_name)
            continue

        got_r = got_o = None
        if series_r is not None:
            t_r, y_r = series_r
            got_r = _fft_proto_hz(t_r, y_r, win_r[0], win_r[1], lab_clock=use_lab)
            m = (t_r >= win_r[0]) & (t_r <= win_r[1]) & np.isfinite(y_r)
            if np.any(m):
                ax_hist.plot(
                    t_r[m],
                    y_r[m],
                    color=COLOR_REF,
                    lw=1.0,
                    label=lab_ref,
                    zorder=2,
                )
        if series_o is not None:
            t_o, y_o = series_o
            got_o = _fft_proto_hz(t_o, y_o, win_o[0], win_o[1], lab_clock=use_lab)
            m = (t_o >= win_o[0]) & (t_o <= win_o[1]) & np.isfinite(y_o)
            if np.any(m):
                ax_hist.plot(
                    t_o[m],
                    y_o[m],
                    color=COLOR_OTHER,
                    lw=0.95,
                    label=lab_oth,
                    zorder=3,
                )

        if got_r is not None:
            f_r, p_r, f_hi_r = got_r
            ax_psd.loglog(f_r, p_r, color=COLOR_REF, lw=1.35, label=lab_ref, zorder=2)
            f_hi_all = min(f_hi_all, f_hi_r)
            any_curve = True
        if got_o is not None:
            f_o, p_o, f_hi_o = got_o
            ax_psd.loglog(f_o, p_o, color=COLOR_OTHER, lw=1.15, label=lab_oth, zorder=3)
            f_hi_all = min(f_hi_all, f_hi_o)
            any_curve = True
        if got_r is None and got_o is None:
            ax_psd.text(
                0.5, 0.5, "short / missing", transform=ax_psd.transAxes, ha="center"
            )

        ax_psd.grid(True, which="both", ls=":", alpha=0.45)
        ax_psd.set_ylabel(f"{row_name}\n" r"(mm$^2$/Hz) proto")
        if i == 0:
            add_dual_freq_xaxis(ax_psd, label=True)
            ax_psd.legend(
                loc="lower left",
                frameon=True,
                fontsize=float(plt.rcParams["legend.fontsize"]),
            )
            ax_hist.set_title(r"FFT window", fontsize=label_fs)
            ax_psd.set_title(r"FFT PSD", fontsize=label_fs)
            ax_hist.legend(
                loc="upper right",
                frameon=True,
                fontsize=float(plt.rcParams["legend.fontsize"]),
            )
        if i == n_row - 1:
            ax_psd.set_xlabel(r"$f$ (Hz) prototype scale")

        ax_hist.set_ylabel(rf"{row_name}" "\n" r"$\Delta u$ (mm) proto")
        ax_hist.grid(True, which="both", ls=":", alpha=0.45)
        if use_lab and i == 0:
            add_dual_time_xaxis_lab(ax_hist, top=True)

    # Shared proto-Hz xlim on every PSD row (min Nyquist across series).
    if any_curve:
        f_max_guide = f_hi_all * 1.05
        for i in range(n_row):
            ax_psd = axes[i, 1]
            mark_freq_guides(
                ax_psd,
                f_wave=None,
                f_pier=f_pier,
                f_mode4=f_mode4,
                f_mode5=f_mode5,
                f_ssi=f_ssi,
                f_max=f_max_guide,
                with_labels=(i == 0),
            )
            ax_psd.set_xlim(F_MIN_PROTO_HZ, f_hi_all)

    # History xlim = FFT window (lab rows share lab window; OS rows share OS).
    i_last_lab = i_last_os = -1
    for i, (_k, _n, lab_clock) in enumerate(ROW_SPECS):
        ax_h = axes[i, 0]
        if lab_clock and win_lab_r is not None and win_lab_o is not None:
            t0 = min(float(win_lab_r[0]), float(win_lab_o[0]))
            t1 = max(float(win_lab_r[1]), float(win_lab_o[1]))
            if t1 > t0:
                ax_h.set_xlim(t0, t1)
                i_last_lab = i
        elif (not lab_clock) and win_os_r is not None and win_os_o is not None:
            t0 = min(float(win_os_r[0]), float(win_os_o[0]))
            t1 = max(float(win_os_r[1]), float(win_os_o[1]))
            if t1 > t0:
                ax_h.set_xlim(t0, t1)
                i_last_os = i
    if i_last_lab >= 0:
        axes[i_last_lab, 0].set_xlabel(LABEL_T_LAB)
    if i_last_os >= 0:
        axes[i_last_os, 0].set_xlabel(LABEL_T_INT)

    if not any_curve:
        plt.close(fig)
        return False

    tag_r = test_ref or pair_file_slug(ref.name)
    tag_o = test_other or pair_file_slug(other.name)
    fig.suptitle(
        rf"FFT PSD {win_label}  ·  ref {tag_r} (grey) vs {tag_o} (navy)"
        r"  ·  lab $t$ / $t_{\mathrm{int}}$"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=PAIR_DPI, bbox_inches="tight", pad_inches=0.22)
    plt.close(fig)
    print(f"PlotEQComparePSD: wrote {out}")
    return True


def write_group_psd_pairs(
    slug: str,
    runs: list[Path],
    run_mat: dict[str, str],
    label: str,
) -> None:
    """
    Pick interim ref and write psd_ux_pair_*.png for each other dump.

    Args:    slug, runs, run_mat, label
    Returns: none
    """
    duration = arias_significant_duration()
    out_dir = compare_plots_dir(slug) / "pairs"
    rows_by_dump = dump_to_row()
    test_ids = dump_to_test_id()
    in_excluded = slug.replace("\\", "/").startswith("_excluded/")
    if not in_excluded:
        runs = [
            r
            for r in runs
            if analysis_skip_reason(
                rows_by_dump.get(r.name), run_mat.get(r.name, "") or ""
            )
            is None
        ]
    if not runs:
        print(f"\nPlotEQComparePSD: group {slug} — all excluded, skip")
        return

    ref, stats = pick_reference(
        runs, run_mat, duration, allow_excluded=in_excluded
    )
    print(f"\nPlotEQComparePSD: group {slug}")
    print(f"  {label}")
    if ref is None:
        print("  no stateOS + meaSigOS — skip")
        return

    incomplete = {eq_dir.name: not run_duration(eq_dir)[2] for eq_dir in runs}
    labels = legend_labels_for_dumps(
        [eq_dir.name for eq_dir in runs],
        incomplete=incomplete,
    )
    mat_ref = run_mat[ref.name]
    test_ref = test_ids.get(ref.name, "")

    others = [r for r in runs if r.name != ref.name]
    if not others:
        print("  only one dump — no pairs")
        return

    n_ok = 0
    for other in others:
        if other.name not in stats:
            print(f"  skip {other.name} (no mat / stateOS)")
            continue
        mat_other = run_mat.get(other.name)
        if not mat_other:
            continue
        tag = pair_file_slug(other.name, test_ids)
        ok = plot_psd_pair(
            ref,
            other,
            out_dir / f"psd_ux_pair_{tag}.png",
            mat_ref=mat_ref,
            mat_other=mat_other,
            labels=labels,
            test_ref=test_ref,
            test_other=test_ids.get(other.name, ""),
            fair_lab_end=_is_wed0819_group(slug),
        )
        if ok:
            n_ok += 1
    print(f"  wrote {n_ok} PSD pair(s) -> {out_dir}")


def main() -> int:
    if any(a in ("-h", "--help") for a in sys.argv[1:]):
        print(HELP, end="")
        return 0

    root = resolve_opensees_data()
    run_mat = run_to_mat_name()
    dump_group = dump_to_group()
    lab_groups = groups_by_dump()

    if len(sys.argv) > 1:
        buckets: dict[str, list[Path]] = {}
        labels: dict[str, str] = {}
        for arg in sys.argv[1:]:
            path = Path(arg).resolve()
            if not path.is_dir():
                print(f"PlotEQComparePSD: skip (not a dir) {arg}", file=sys.stderr)
                continue
            slug = dump_group.get(path.name, "cli")
            buckets.setdefault(slug, []).append(path)
            labels.setdefault(slug, slug)
        if not buckets:
            return 1
        for slug, runs in buckets.items():
            write_group_psd_pairs(slug, runs, run_mat, labels[slug])
        return 0

    if not lab_groups:
        print("PlotEQComparePSD: empty TestMatrix_lab_runs.csv", file=sys.stderr)
        return 1

    n_groups = 0
    for slug, rows in lab_groups.items():
        runs: list[Path] = []
        for row in rows:
            path = resolve_run_path(row["DumpFolder"], root)
            if path is None:
                continue
            runs.append(path)
        if len(runs) < 2:
            continue
        write_group_psd_pairs(slug, runs, run_mat, group_label(rows[0]))
        n_groups += 1

    print(f"\nPlotEQComparePSD: {n_groups} group folder(s)")
    return 0 if n_groups else 1


if __name__ == "__main__":
    raise SystemExit(main())
