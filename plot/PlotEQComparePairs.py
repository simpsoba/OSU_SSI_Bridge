#!/usr/bin/env python3
"""
Goals
-----
Pairwise Simulink compare plots within each physical-model group.

  Reference (interim) = dump with the fewest typeConv3==2 rising edges inside
  the ground-motion D5–95 window (Arias on FKSH19.NS1.VT2). Prefer complete
  analyses when counts tie; require meaSigOS; exclude testBaseline mats,
  single-precision CuDSS (unstable OpenSees recorders), and twoNodeLink
  (those land under plots/compare/_excluded/…).

  True reference (later) = offline OpenSees for the same folder with
  realTimeON 0 (no OpenFresco / Simulink). Until those exist, fewest
  D5–95 slowdowns is the interim reference.

  python plot/PlotEQComparePairs.py
  python plot/PlotEQComparePairs.py <runDir> ...

Writes under OSU_SSI_BRIDGE_DATA_LOCAL/plots/compare/<Mesh>/<variant>/pairs/ :

  hist_ux_pair_F##.png            Simulink only (mandatory; meaSigOS + stateOS)
  hist_ux_pair_F##_opensees.png   OpenSees only (optional; pier_top vs t_num)
  reference.txt                       chosen ref + D5–95 slowdown counts

Mesh folders: Baseline / Moderate / Large / X-Large (excluded → _excluded/).

Do not mix OpenSees and Simulink on one figure — hybrid clocks desync.

  Ref = grey (#B0B0B0); other = navy (#001F3F) — BRB-Calibration pair.
  Axis notation (matched): prototype $t$, $\Delta u$; model $t/\sqrt{\lambda}$,
  $\Delta u/\lambda$. Simulink $\Delta u$ = actuator (meaSigOS); OpenSees =
  pier-top nodal.

Units
-----
  Simulink Time: t_lab (s, model). Primary x = t_lab√λ (prototype).
  meaSigOS disp: model m → primary Δu (mm, prototype) via ·λ·1000.
  stateOS: same t_lab clock as meaSigOS.
  OpenSees companion: t_num (s, prototype), pier Δux (mm).
  D5–95 on GM / prototype clock: [t5, t95] from gm_duration (gmStart≈0).
  Lab D5–95 window for counts: t_lab ∈ [t5/√λ, t95/√λ].

Loading note (2026-08-21): folder nicknames may say Storm_Wave; the campaign
is earthquake followed by tsunami (not storm-wave-only).
"""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MaxNLocator
import numpy as np

from compare_groups import (
    analysis_skip_reason,
    dump_to_group,
    dump_to_row,
    dump_to_test_id,
    group_label,
    groups_by_dump,
    legend_labels_for_dumps,
    test_file_slug,
)
from gm_duration import SignificantDuration, arias_significant_duration
from lab_paths import (
    CYLINDER_LENGTH_SCALE,
    MAT_EXTRACT_DIR,
    M_TO_MM,
    TIME_SCALE_FROUDE,
    compare_plots_dir,
    resolve_opensees_data,
)
from paths import HERE
from PlotEQCompareRuns import (
    DPI,
    UX_ABS_MAX_MM,
    apply_paper_style,
    find_pier_top,
    load_mat_mea_feedback,
    loadtxt_partial,
    model_disp_to_proto_mm,
    resolve_run_path,
    run_duration,
    run_to_mat_name,
    sym_ylim,
)
from PlotMatOS import load_block, non_time_cols

# ------------------------------------------------------------
# knobs
# ------------------------------------------------------------
REPO = HERE.parent
SLOWDOWN_STATE = 2  # typeConv3 wait flag (OpenFresco)
PAIR_DPI = DPI

# Pair colors (same pair as BRB-Calibration plot_dimensions):
# grey = reference / baseline; navy = comparison series.
COLOR_REF = "#B0B0B0"
COLOR_OTHER = "#001F3F"
# Fixed prototype-time window for pair figures (matches campaign Trec ≈ 450 s).
XLIM_PROTO_S = (0.0, 450.0)

HELP = """\
usage: python plot/PlotEQComparePairs.py [runDir ...]

  no args   dumps from TestMatrix_lab_runs.csv (grouped by mesh + model knobs)
  runDir    one or more dump folders (still grouped if in lab_runs)
  output    LOCAL/plots/compare/<Baseline|Moderate|Large|X-Large>/<variant>/pairs/
            (excluded runs → compare/_excluded/<reason>/<variant>/pairs/)
            hist_ux_pair_F##.png           Simulink (primary)
            hist_ux_pair_F##_opensees.png  OpenSees companion
  ref       fewest typeConv3→2 events in GM D5-95 (complete preferred)
  note      offline OpenSees (realTimeON 0) is the eventual true reference
"""


# ------------------------------------------------------------
# 1. STATEOS SLOWDOWNS IN D5–95
# ------------------------------------------------------------


@dataclass
class SlowdownStats:
    """Slowdown rising edges for one mat/dump on the lab clock."""

    mat_name: str
    n_d595: int
    n_all: int
    t_onset_lab_s: np.ndarray  # all rising edges (lab s)
    t_onset_d595_lab_s: np.ndarray  # those inside D5–95


def load_type_conv3(mat_name: str) -> tuple[np.ndarray, np.ndarray] | None:
    """
    Lab Time and integer typeConv3 from a mat extract.

    Args:    mat_name  Simulink .mat file name
    Returns: (t_lab_s, typeConv3) or None
    """
    npz_path = MAT_EXTRACT_DIR / f"{Path(mat_name).stem}.npz"
    if not npz_path.is_file():
        return None
    z = np.load(npz_path, allow_pickle=True)
    block = load_block(z, "stateOS")
    if block is None:
        return None
    t_lab_s, names, data = block
    cols = non_time_cols(names, data.shape[1])
    if not cols:
        return None
    return t_lab_s, np.rint(data[:, cols[0]]).astype(int)


def rising_edges_to_state(
    t_lab_s: np.ndarray,
    state: np.ndarray,
    target: int = SLOWDOWN_STATE,
) -> np.ndarray:
    """
    Lab times of rising edges into `target` (e.g. typeConv3 → 2).

    Args:    t_lab_s, state, target
    Returns: 1-d array of onset times (lab s)
    """
    if state.size < 2:
        return np.asarray([], dtype=float)
    enter = (state[1:] == target) & (state[:-1] != target)
    return np.asarray(t_lab_s[1:][enter], dtype=float)


def slowdown_stats_for_mat(
    mat_name: str,
    duration: SignificantDuration,
) -> SlowdownStats | None:
    """
    Count slowdown onsets inside the lab-mapped D5–95 window.

    Lab window: [t5/√λ, t95/√λ] so it matches GM [t5, t95] under Froude.

    Args:    mat_name, duration  Arias D5–95 on GM / prototype clock
    Returns: SlowdownStats or None if no stateOS
    """
    pair = load_type_conv3(mat_name)
    if pair is None:
        return None
    t_lab_s, state = pair
    onsets = rising_edges_to_state(t_lab_s, state)
    t5_lab = duration.t5_s / TIME_SCALE_FROUDE
    t95_lab = duration.t95_s / TIME_SCALE_FROUDE
    in_window = (onsets >= t5_lab) & (onsets <= t95_lab)
    return SlowdownStats(
        mat_name=mat_name,
        n_d595=int(np.count_nonzero(in_window)),
        n_all=int(onsets.size),
        t_onset_lab_s=onsets,
        t_onset_d595_lab_s=onsets[in_window],
    )


def is_single_precision_dump(dump_name: str, rows_by_dump: dict[str, dict[str, str]]) -> bool:
    """
    True if the as-run matrix marks CuDSS single precision (dFFI).

    Those runs often blew up / NaN'd in OpenSees; do not use as interim ref.

    Args:    dump_name; rows_by_dump  from dump_to_row()
    Returns: bool
    """
    row = rows_by_dump.get(dump_name) or {}
    blob = " ".join(
        [
            row.get("Goal", ""),
            row.get("Name", ""),
            row.get("postPartitionSystem", ""),
            row.get("prePartitionSystem", ""),
        ]
    ).lower()
    return ("-precision" in blob) or ("single precision" in blob)


def pick_reference(
    runs: list[Path],
    run_mat: dict[str, str],
    duration: SignificantDuration,
    *,
    allow_excluded: bool = False,
) -> tuple[Path | None, dict[str, SlowdownStats]]:
    """
    Interim reference: fewest D5–95 slowdowns; complete preferred on ties.

    Requires plottable Simulink meaSigOS. Unless ``allow_excluded``, skip
    dry/PID baseline mats, single-precision CuDSS, and twoNodeLink.

    Args:    runs, run_mat, duration; allow_excluded for ``_excluded/`` groups
    Returns: (ref_path or None, {dump_name: SlowdownStats})
    """
    stats: dict[str, SlowdownStats] = {}
    rows_by_dump = dump_to_row()
    # sort key: n_d595, incomplete, is_baseline, dump name
    candidates: list[tuple[int, int, int, str, Path]] = []

    for eq_dir in runs:
        mat_name = run_mat.get(eq_dir.name)
        if not mat_name:
            continue
        row = rows_by_dump.get(eq_dir.name)
        if not allow_excluded and analysis_skip_reason(row, mat_name):
            continue
        s = slowdown_stats_for_mat(mat_name, duration)
        if s is None:
            continue
        stats[eq_dir.name] = s
        if load_mea_ux_proto(mat_name) is None:
            continue
        _, _, complete = run_duration(eq_dir)
        candidates.append(
            (s.n_d595, 0 if complete else 1, 0, eq_dir.name, eq_dir)
        )

    if not candidates:
        return None, stats
    candidates.sort(key=lambda row: (row[0], row[1], row[2], row[3]))
    return candidates[0][4], stats


# ------------------------------------------------------------
# 2. SERIES HELPERS
# ------------------------------------------------------------


def load_mea_ux_proto(mat_name: str) -> tuple[np.ndarray, np.ndarray] | None:
    """
    Simulink meaSigOS Δu (mm, prototype) vs t_lab√λ (s, prototype).

    Args:    mat_name
    Returns: (t_proto_s, ux_proto_mm) or None
    """
    pair = load_mat_mea_feedback(mat_name)
    if pair is None:
        return None
    t_lab_s, u_model_m = pair
    if t_lab_s.size == 0:
        return None
    t_proto_s = t_lab_s * TIME_SCALE_FROUDE
    ux_mm = model_disp_to_proto_mm(u_model_m)
    if not np.all(np.isfinite(ux_mm)) or float(np.nanmax(np.abs(ux_mm))) > UX_ABS_MAX_MM:
        return None
    return t_proto_s, ux_mm


def load_pier_ux_mm(eq_dir: Path) -> tuple[np.ndarray, np.ndarray] | None:
    """
    OpenSees pier-top Δux (mm) vs t_num (s).

    Args:    eq_dir
    Returns: (t_num_s, ux_mm) or None
    """
    pier = find_pier_top(eq_dir)
    if pier is None:
        return None
    data = loadtxt_partial(pier)
    if data.size == 0 or data.shape[1] < 2:
        return None
    t_s = data[:, 0]
    ux_mm = (data[:, 1] - data[0, 1]) * M_TO_MM
    if not np.all(np.isfinite(ux_mm)) or float(np.nanmax(np.abs(ux_mm))) > UX_ABS_MAX_MM:
        return None
    return t_s, ux_mm


def apply_pair_style() -> None:
    """
    Paper style with ~1.5× base font sizes for pair figures.

    Args:    none
    Returns: none
    """
    apply_paper_style()
    plt.rcParams.update(
        {
            "font.size": 14,
            "axes.labelsize": 15,
            "axes.titlesize": 14,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "legend.fontsize": 10,
            "lines.linewidth": 1.15,
        }
    )


def add_dual_disp_yaxis(
    ax: plt.Axes,
    *,
    source: str,
    primary_label: bool = True,
) -> None:
    """
    Left = Δu prototype; right = Δu/λ model (same slash notation as time).

    Args:    ax; source  "actuator" (Simulink meaSigOS) or "pier" (OpenSees);
             primary_label  if False, omit left ylabel (shared row)
    Returns: none
    """
    if primary_label:
        if source == "actuator":
            ax.set_ylabel(
                r"$\Delta u$ (mm) prototype scale"
                "\n"
                r"actuator (meaSigOS)"
            )
        else:
            ax.set_ylabel(
                r"$\Delta u$ (mm) prototype scale"
                "\n"
                r"pier top (OpenSees)"
            )
    sec_y = ax.secondary_yaxis(
        "right",
        functions=(
            lambda u_proto: u_proto / CYLINDER_LENGTH_SCALE,
            lambda u_model: u_model * CYLINDER_LENGTH_SCALE,
        ),
    )
    sec_y.set_ylabel(r"$\Delta u/\lambda$ (mm) model scale")


def add_dual_time_xaxis(ax: plt.Axes, *, top: bool) -> None:
    """
    Top secondary: $t/\sqrt{\lambda}$ (s) model scale.

    Args:    ax; top  if True, add secondary model-time axis above
    Returns: none
    """
    if top:
        sec_x = ax.secondary_xaxis(
            "top",
            functions=(
                lambda t_proto: t_proto / TIME_SCALE_FROUDE,
                lambda t_model: t_model * TIME_SCALE_FROUDE,
            ),
        )
        sec_x.set_xlabel(r"$t/\sqrt{\lambda}$ (s) model scale")


LABEL_T_PROTO = r"$t$ (s) prototype scale"


def format_state_axis(ax: plt.Axes, y_stack: list[np.ndarray]) -> None:
    """
    Integer typeConv3 limits and ticks.

    Args:    ax; y_stack  state series used for ylim
    Returns: none
    """
    if y_stack:
        y_all = np.concatenate(y_stack)
        y_min = int(np.floor(float(np.nanmin(y_all))))
        y_max = int(np.ceil(float(np.nanmax(y_all))))
        if y_max <= y_min:
            y_max = y_min + 1
        ax.set_ylim(y_min - 0.25, y_max + 0.35)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, min_n_ticks=2))
    ax.yaxis.set_major_formatter(
        FuncFormatter(lambda v, _: f"{int(round(v))}")
    )
    ax.set_ylabel("typeConv3")


def pair_file_slug(dump_name: str, test_ids: dict[str, str] | None = None) -> str:
    """
    File tag for one dump: Test ID when known, else row+HHMM fallback.

    Args:    dump_name  DumpFolder; test_ids  optional dump→Test map
    Returns: e.g. F26 or r-02_1601
    """
    ids = test_ids if test_ids is not None else dump_to_test_id()
    tid = ids.get(dump_name)
    if tid is not None:
        return test_file_slug(tid)
    match = re.match(r"^(r[+-]?\d+)_(\d{8})_(\d{4})_", dump_name, flags=re.IGNORECASE)
    if match:
        return f"{match.group(1)}_{match.group(3)}"
    return re.sub(r"[^A-Za-z0-9.+-]+", "_", dump_name)[:48]


def dump_slug(name: str) -> str:
    """Deprecated alias for pair_file_slug (kept for callers)."""
    return pair_file_slug(name)


def _pair_line_styles(ref: Path, other: Path) -> tuple[str, str]:
    """
    Solid if OpenSees analysis reached Trec; dashed if incomplete.

    Args:    ref, other
    Returns: (ls_ref, ls_other)
    """
    _, _, complete_ref = run_duration(ref)
    _, _, complete_other = run_duration(other)
    return ("-" if complete_ref else "--", "-" if complete_other else "--")


# ------------------------------------------------------------
# 3. PAIRWISE FIGURES (Simulink primary; OpenSees companion)
# ------------------------------------------------------------


def plot_pair_simulink(
    ref: Path,
    other: Path,
    out: Path,
    duration: SignificantDuration,
    stats: dict[str, SlowdownStats],
    labels: dict[str, str],
) -> None:
    """
    All-Simulink 2×2 figure:

      Δu full | Δu D5–95 zoom
      typeConv3 full | typeConv3 D5–95 zoom

    Args:    ref, other, out, duration, stats, labels
    Returns: none (writes PNG)
    """
    mat_ref = stats[ref.name].mat_name
    mat_other = stats[other.name].mat_name
    mea_ref = load_mea_ux_proto(mat_ref)
    mea_other = load_mea_ux_proto(mat_other)
    if mea_ref is None or mea_other is None:
        print(f"PlotEQComparePairs: skip Simulink pair {other.name} (no meaSigOS)")
        return

    state_ref = load_type_conv3(mat_ref)
    state_other = load_type_conv3(mat_other)

    apply_pair_style()
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(11.2, 6.2),
        sharex="col",
        sharey="row",
        gridspec_kw={
            "height_ratios": [2.2, 1.3],
            "width_ratios": [1.35, 1.0],
            "hspace": 0.18,
            "wspace": 0.28,
        },
    )
    ax_ux_f, ax_ux_z = axes[0]
    ax_st_f, ax_st_z = axes[1]

    ls_ref, ls_other = _pair_line_styles(ref, other)
    label_ref = f"ref  {labels[ref.name]}"
    label_other = labels[other.name]
    t_ref, u_ref = mea_ref
    t_oth, u_oth = mea_other
    t5, t95 = duration.t5_s, duration.t95_s

    for ax in (ax_ux_f, ax_ux_z):
        ax.plot(t_ref, u_ref, color=COLOR_REF, lw=1.2, ls=ls_ref, zorder=2)
        ax.plot(t_oth, u_oth, color=COLOR_OTHER, lw=1.1, ls=ls_other, zorder=3)
    sym_ylim(ax_ux_f)

    y_stack: list[np.ndarray] = []
    for ax in (ax_st_f, ax_st_z):
        if state_ref is not None:
            t_lab, st = state_ref
            ax.plot(
                t_lab * TIME_SCALE_FROUDE,
                st,
                color=COLOR_REF,
                lw=1.1,
                drawstyle="steps-post",
                zorder=2,
            )
            if ax is ax_st_f:
                y_stack.append(st)
        if state_other is not None:
            t_lab, st = state_other
            ax.plot(
                t_lab * TIME_SCALE_FROUDE,
                st,
                color=COLOR_OTHER,
                lw=1.05,
                drawstyle="steps-post",
                zorder=3,
            )
            if ax is ax_st_f:
                y_stack.append(st)
        ax.axhline(SLOWDOWN_STATE, color="0.65", lw=0.7, ls=":", zorder=1)
    format_state_axis(ax_st_f, y_stack)
    ax_st_z.set_ylabel("")

    add_dual_disp_yaxis(ax_ux_f, source="actuator", primary_label=True)
    add_dual_disp_yaxis(ax_ux_z, source="actuator", primary_label=False)
    add_dual_time_xaxis(ax_ux_f, top=True)
    add_dual_time_xaxis(ax_ux_z, top=True)
    ax_ux_z.set_title(r"D5–95 zoom", fontsize=14, pad=6)
    ax_st_f.set_xlabel(LABEL_T_PROTO)
    ax_st_z.set_xlabel(LABEL_T_PROTO)

    ax_ux_f.set_xlim(*XLIM_PROTO_S)
    ax_st_f.set_xlim(*XLIM_PROTO_S)
    ax_ux_z.set_xlim(t5, t95)
    ax_st_z.set_xlim(t5, t95)

    fig.legend(
        handles=[
            Line2D([0], [0], color=COLOR_REF, lw=1.2, ls=ls_ref, label=label_ref),
            Line2D(
                [0], [0], color=COLOR_OTHER, lw=1.1, ls=ls_other, label=label_other
            ),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.02),
        ncol=2,
        fontsize=10,
        frameon=False,
        columnspacing=1.2,
        handlelength=1.6,
        handletextpad=0.4,
    )

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=PAIR_DPI, bbox_inches="tight", pad_inches=0.22)
    plt.close(fig)
    print(f"PlotEQComparePairs: wrote {out}  (Simulink)")


def plot_pair_opensees(
    ref: Path,
    other: Path,
    out: Path,
    duration: SignificantDuration,
    labels: dict[str, str],
) -> None:
    """
    Optional OpenSees-only companion: pier Δu full | D5–95 zoom.

    Args:    ref, other, out, duration, labels
    Returns: none (writes PNG when both piers exist)
    """
    pier_ref = load_pier_ux_mm(ref)
    pier_other = load_pier_ux_mm(other)
    if pier_ref is None or pier_other is None:
        print(f"PlotEQComparePairs: skip OpenSees pair {other.name} (no pier ux)")
        return

    apply_pair_style()
    fig, (ax_f, ax_z) = plt.subplots(
        1,
        2,
        figsize=(11.2, 3.8),
        sharey=True,
        gridspec_kw={"width_ratios": [1.35, 1.0], "wspace": 0.28},
    )
    ls_ref, ls_other = _pair_line_styles(ref, other)
    label_ref = f"ref  {labels[ref.name]}"
    label_other = labels[other.name]
    t_ref, u_ref = pier_ref
    t_oth, u_oth = pier_other
    t5, t95 = duration.t5_s, duration.t95_s

    for ax in (ax_f, ax_z):
        ax.plot(t_ref, u_ref, color=COLOR_REF, lw=1.2, ls=ls_ref, zorder=2)
        ax.plot(t_oth, u_oth, color=COLOR_OTHER, lw=1.1, ls=ls_other, zorder=3)
    sym_ylim(ax_f)

    add_dual_disp_yaxis(ax_f, source="pier", primary_label=True)
    add_dual_disp_yaxis(ax_z, source="pier", primary_label=False)
    add_dual_time_xaxis(ax_f, top=True)
    add_dual_time_xaxis(ax_z, top=True)
    ax_f.set_xlabel(LABEL_T_PROTO)
    ax_z.set_xlabel(LABEL_T_PROTO)
    ax_z.set_title(r"D5–95 zoom", fontsize=14, pad=6)
    ax_f.set_xlim(*XLIM_PROTO_S)
    ax_z.set_xlim(t5, t95)

    fig.legend(
        handles=[
            Line2D([0], [0], color=COLOR_REF, lw=1.2, ls=ls_ref, label=label_ref),
            Line2D(
                [0], [0], color=COLOR_OTHER, lw=1.1, ls=ls_other, label=label_other
            ),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.1),
        ncol=2,
        fontsize=10,
        frameon=False,
        columnspacing=1.2,
        handlelength=1.6,
        handletextpad=0.4,
    )

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=PAIR_DPI, bbox_inches="tight", pad_inches=0.22)
    plt.close(fig)
    print(f"PlotEQComparePairs: wrote {out}  (OpenSees)")


def write_group_pairs(
    slug: str,
    runs: list[Path],
    run_mat: dict[str, str],
    label: str,
    duration: SignificantDuration,
) -> None:
    """
    Pick interim ref and write Simulink (+ optional OpenSees) PNGs per other.

    Campaign groups drop baseline / precision / twoNodeLink. Groups already
    under ``_excluded/`` keep those dumps so pairs still land there.

    Args:    slug, runs, run_mat, label, duration
    Returns: none
    """
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
        print(f"\nPlotEQComparePairs: group {slug}")
        print(f"  {label}")
        print("  all dumps excluded (baseline / precision / twoNodeLink) — skip")
        return

    ref, stats = pick_reference(
        runs, run_mat, duration, allow_excluded=in_excluded
    )
    print(f"\nPlotEQComparePairs: group {slug}")
    print(f"  {label}")
    if ref is None:
        print("  no stateOS + meaSigOS — skip pairs")
        return

    incomplete = {
        eq_dir.name: not run_duration(eq_dir)[2] for eq_dir in runs
    }
    labels = legend_labels_for_dumps(
        [eq_dir.name for eq_dir in runs],
        incomplete=incomplete,
    )

    lines = [
        f"group: {slug}",
        f"label: {label}",
        (
            f"D5-95 (GM / prototype): t5={duration.t5_s:.3f}s  "
            f"t95={duration.t95_s:.3f}s  D={duration.d5_95_s:.3f}s  "
            f"({duration.path.name})"
        ),
        (
            "interim_reference: fewest typeConv3→2 rising edges in D5-95 "
            "(complete preferred; meaSigOS required; "
            "exclude single-precision CuDSS and twoNodeLink "
            "from campaign folders — those go under _excluded/; "
            "dry OpenFresco baseline F01 is included)"
        ),
        (
            "primary_figures: Simulink meaSigOS + stateOS "
            "(same lab clock; do not mix with OpenSees)"
        ),
        "companion_figures: *_opensees.png = pier_top vs t_num only",
        "pair_file_tag: Test ID (F##/W##) from TestMatrix_lab_runs.csv",
        (
            "true_reference_later: offline OpenSees per folder "
            "(realTimeON 0; no OpenFresco / Simulink)"
        ),
        (
            "loading: 2026-08-21 campaign is earthquake followed by tsunami "
            "(folder Storm_Wave is a matrix nickname, not storm-wave-only)"
        ),
        f"reference_dump: {ref.name}",
        f"reference_test: {test_ids.get(ref.name, '')}",
        f"reference_label: {labels.get(ref.name, ref.name)}",
        f"reference_n_d595: {stats[ref.name].n_d595}",
        f"reference_n_all: {stats[ref.name].n_all}",
        "",
        "dump\ttest\tn_d595\tn_all\tcomplete",
    ]
    for eq_dir in sorted(runs, key=lambda p: p.name):
        s = stats.get(eq_dir.name)
        _, _, complete = run_duration(eq_dir)
        n_d = s.n_d595 if s else -1
        n_a = s.n_all if s else -1
        tid = test_ids.get(eq_dir.name, "")
        lines.append(
            f"{eq_dir.name}\t{tid}\t{n_d}\t{n_a}\t{int(complete)}"
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    ref_txt = out_dir / "reference.txt"
    ref_txt.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"  ref = {ref.name}  (n_d595={stats[ref.name].n_d595}) -> {ref_txt}")

    others = [r for r in runs if r.name != ref.name]
    if not others:
        print("  only one dump with stateOS — no pairs")
        return

    for other in others:
        if other.name not in stats:
            print(f"  skip {other.name} (no mat / stateOS)")
            continue
        other_tag = pair_file_slug(other.name, test_ids)
        plot_pair_simulink(
            ref,
            other,
            out_dir / f"hist_ux_pair_{other_tag}.png",
            duration,
            stats,
            labels,
        )
        plot_pair_opensees(
            ref,
            other,
            out_dir / f"hist_ux_pair_{other_tag}_opensees.png",
            duration,
            labels,
        )


# ------------------------------------------------------------
# 4. DISCOVER / MAIN
# ------------------------------------------------------------


def main() -> int:
    if any(a in ("-h", "--help") for a in sys.argv[1:]):
        print(HELP, end="")
        return 0

    duration = arias_significant_duration()
    print(
        f"PlotEQComparePairs: D5-95 = {duration.d5_95_s:.2f} s "
        f"[{duration.t5_s:.1f}, {duration.t95_s:.1f}] from {duration.path.name}"
    )

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
                print(f"PlotEQComparePairs: skip (not a dir) {arg}", file=sys.stderr)
                continue
            slug = dump_group.get(path.name, "cli")
            buckets.setdefault(slug, []).append(path)
            labels.setdefault(slug, slug)
        if not buckets:
            return 1
        for slug, runs in buckets.items():
            write_group_pairs(slug, runs, run_mat, labels[slug], duration)
        return 0

    if not lab_groups:
        print("PlotEQComparePairs: empty TestMatrix_lab_runs.csv", file=sys.stderr)
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
        write_group_pairs(slug, runs, run_mat, group_label(rows[0]), duration)
        n_groups += 1

    print(f"\nPlotEQComparePairs: {n_groups} group folder(s)")
    return 0 if n_groups else 1


if __name__ == "__main__":
    raise SystemExit(main())
