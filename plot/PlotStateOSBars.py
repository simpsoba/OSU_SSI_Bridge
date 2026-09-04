#!/usr/bin/env python3
"""
Goals
-----
Stacked bar charts: fraction of stateOS samples in each OpenFresco
typeConv3 stage (Seki §2.1) per campaign **Test** ID (F## / W## / Fd##).

  python plot/PlotStateOSBars.py

Reads:
  plot/lab/TestMatrix_lab_runs.csv  (Test column = F/W/Fd/Fx IDs)
  LOCAL/mat_extract/*.npz           (stateOS via TestMatrix_lab_runs.csv)

Writes:
  OSU_SSI_BRIDGE_DATA_LOCAL/plots/compare/stateos/bar_typeconv3_eq35_90_wed.png
  OSU_SSI_BRIDGE_DATA_LOCAL/plots/compare/stateos/bar_typeconv3_eq35_90_fri.png
  OSU_SSI_BRIDGE_DATA_LOCAL/plots/compare/stateos/bar_typeconv3_full_wed.png
  OSU_SSI_BRIDGE_DATA_LOCAL/plots/compare/stateos/bar_typeconv3_full_fri.png
  …_with_dry.png companions (Fri includes Fd## dry mat-only)

Windows:
  eq35_90 — fixed lab-clock band [35, 90] s (semi–D5–95). Not Arias /
            Froude-mapped GM time: slowdowns stretch wall clock, so a fixed
            lab window is the comparable “most of the EQ” band. Early stops
            show an unfinished (hatched) top segment so every bar is 100% of
            that window.
  full    — every stateOS sample in the mat extract (lab clock); no unfinished.
Basis: sample count (DAQ rate) for stages; typeConv3==2 split into brief (1-tick)
vs sustained (≥2 ticks) episodes; unfinished by time in the eq window.
See STATEOS_SIGNALS.md for stage codes.
"""

from __future__ import annotations

import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patheffects
import numpy as np

from compare_groups import (
    abbreviate_integrator,
    abbreviate_solver,
    analysis_skip_reason,
    dump_to_row,
    has_hybrid_execute,
    is_dry_mat_only,
    mesh_neqn,
)
from lab_paths import lab_runs_csv_path, stateos_compare_plots_dir
from PlotEQComparePairs import load_type_conv3
from PlotEQCompareRuns import DPI, run_to_mat_name

# Stack bottom→top = extrapolate → interpolate → brief slowdown → sustained slowdown.
# typeConv3==2 at count 9 (0.9 Δt_sim): 1-sample episode = brief; ≥2 = sustained hold.
# Initialize (−1) omitted from the denominator.
STAGE_KEY = int | str
STAGES: tuple[tuple[STAGE_KEY, str, str, str], ...] = (
    (1, "extrapolate", "#1565c0", "#ffffff"),
    (0, "interpolate", "#2e7d32", "#ffffff"),
    ("slowdown_brief", "slowdown (brief)", "#FFC04D", "#996600"),
    ("slowdown_sustained", "slowdown (sustained)", "#FF8C00", "#5C2200"),
)
INIT_CODE = -1
SLOWDOWN_KEYS = frozenset({"slowdown_brief", "slowdown_sustained"})

# Fixed lab-clock window (s): expected “most of EQ” band when the schedule
# is healthy. Rounder stand-in for mapped Arias ≈ 36–90 s.
EQ_WINDOW_LAB_S = (35.0, 90.0)
UNFINISHED_FACE = "#616161"
UNFINISHED_EDGE = "#212121"
UNFINISHED_HATCH = "///"

# Show in-bar % when segment is at least this tall (always label slowdown > 0).
PCT_LABEL_MIN = 3.0
SLOWDOWN_LABEL_COLOR = "#ffffff"
FS = 8  # one size for ticks, legend, in-bar %, and key table

# Footnote key table columns (blank Notes = campaign default).
KEY_HEADERS = (
    "Test",
    "Mesh",
    "DOFs",
    "Soil",
    "Ele",
    "Integrator",
    "Solver",
    "np",
    "Notes",
)

# Bar order: mesh small→large; within mesh standard Soft/SSPQuad by np, variants last.
# twoNodeLink / single-precision omitted (analysis_skip_reason).
MESH_RANK: dict[str, int] = {
    "Wed0819": -1,
    "Baseline": 0,
    "Moderate": 1,
    "Large": 2,
    "X-Large": 3,
}
BAR_ORDER_BY_MESH: dict[str, tuple[str, ...]] = {
    "Wed0819": ("W01", "W02", "W03", "W04", "W05", "W06", "W07"),
    # F01 = dry OpenFresco baseline (0836); kept in the main campaign set.
    "Baseline": (
        "F01",
        "F03",
        "F04",
        "F05",
        "F27",
        "F22",
        "F26",
        "F20",
        "F23",
        "F24",
        "F30",
    ),
    "Moderate": ("F06", "F17", "F07", "F08"),
    "Large": ("F09", "F10", "F11", "F12", "F14", "F15"),
    "X-Large": ("F13",),
}
MESH_GROUP_GAP = 0.55  # x-axis spacer between mesh blocks

# W## vs F## / Fd## / Fx## — separate PNGs per lab day.
CAMPAIGN_DAYS: tuple[tuple[str, str], ...] = (
    ("wed", "Wed"),
    ("fri", "Fri"),
)


def campaign_day(test_id: str) -> str | None:
    """
    Lab-day bucket for one Test ID.

    Args:    test_id  W## | F## | Fd## | Fx##
    Returns: ``wed`` | ``fri`` | None
    """
    tid = (test_id or "").strip()
    if tid.startswith("W"):
        return "wed"
    if tid.startswith("F"):
        return "fri"
    return None


def filter_trials_by_day(
    trials: list[TrialFractions],
    day: str,
) -> list[TrialFractions]:
    """Keep trials whose Test ID belongs to one lab-day bucket."""
    return [t for t in trials if campaign_day(t.test_id) == day]


@dataclass
class TrialFractions:
    """typeConv3 sample fractions for one Test ID."""

    test_id: str
    dump_folder: str
    mat_name: str
    fracs: dict[STAGE_KEY, float]  # stage key → fraction in [0, 1]
    n_samples: int
    unfinished: float = 0.0  # fraction of eq window not recorded (0 on full)


def _empty_fracs() -> dict[STAGE_KEY, float]:
    return {code: 0.0 for code, _, _, _ in STAGES}


def _slowdown_episode_masks(state: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Brief vs sustained masks for typeConv3==2 samples.

    Sustained = episode length ≥2 DAQ samples (stayed in state 2).

    Args:    state  integer typeConv3 series (full length)
    Returns: (brief_mask, sustained_mask) booleans, True only where state==2
    """
    n = state.size
    brief = np.zeros(n, dtype=bool)
    sustained = np.zeros(n, dtype=bool)
    i = 0
    while i < n:
        if state[i] != 2:
            i += 1
            continue
        j = i
        while j < n and state[j] == 2:
            j += 1
        if j - i >= 2:
            sustained[i:j] = True
        else:
            brief[i:j] = True
        i = j
    return brief, sustained


def fractions_by_code(state: np.ndarray) -> tuple[dict[STAGE_KEY, float], int]:
    """
    Sample fractions of each plotted typeConv3 value on a state slice.

    Initialize (−1) samples are dropped from the denominator.

    Args:    state  integer typeConv3 series
    Returns: ({stage: fraction}, n_samples used)
    """
    state = np.rint(state).astype(int)
    use = state != INIT_CODE
    n = int(use.sum())
    if n == 0:
        return _empty_fracs(), 0
    brief_m, sust_m = _slowdown_episode_masks(state)
    return {
        1: float(np.count_nonzero(use & (state == 1))) / n,
        0: float(np.count_nonzero(use & (state == 0))) / n,
        "slowdown_brief": float(np.count_nonzero(use & brief_m)) / n,
        "slowdown_sustained": float(np.count_nonzero(use & sust_m)) / n,
    }, n


def load_lab_rows() -> list[dict[str, str]]:
    """
    As-run CSV rows sorted by DateTime (campaign Test = F/W/Fd/Fx).

    Args:    none
    Returns: list of row dicts
    """
    path = lab_runs_csv_path()
    if not path.is_file():
        raise FileNotFoundError(f"missing {path}")
    with path.open(newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return []
    if "Test" not in rows[0]:
        raise KeyError(f"{path}: missing Test column (F## / W## / Fd## / Fx##)")
    out: list[dict[str, str]] = []
    for row in rows:
        test_s = (row.get("Test") or "").strip()
        if not test_s:
            continue
        row["Test"] = test_s
        out.append(row)

    def _sort_key(r: dict[str, str]) -> tuple:
        dt = (r.get("DateTime") or "").strip()
        return (dt or "9999", r["Test"])

    out.sort(key=_sort_key)
    return out


def is_excluded_row(row: dict[str, str], mat_name: str) -> str | None:
    """
    Skip reason for twoNodeLink / single-precision, or None if OK.

    Args:    row, mat_name
    Returns: reason string or None
    """
    return analysis_skip_reason(row, mat_name)


def mesh_short(text: str) -> str:
    """
    soilMesh cell → Baseline | Moderate | Large | X-Large.

    Args:    text  e.g. ``0 (BASELINE)``
    Returns: short mesh name (legacy PRODUCTION → Baseline)
    """
    m = re.match(r"^(\d+)\s*\(([^)]+)\)", (text or "").strip())
    if not m:
        return (text or "?").strip() or "?"
    name = m.group(2).replace("XLARGE", "X-large").title()
    if name.lower() in ("production", "baseline"):
        return "Baseline"
    if name.upper() in ("WED0819", "WED"):
        return "Wed0819"
    return name


def soil_short(text: str) -> str:
    """
    soilProfile cell → SOFT | STIFF | …

    Args:    text  e.g. ``4 (SOFT)``
    Returns: parenthetical tag or raw cell
    """
    m = re.match(r"^(\d+)\s*\(([^)]+)\)", (text or "").strip())
    return m.group(2) if m else ((text or "?").strip() or "?")


def run_time_hhmm(row: dict[str, str]) -> str:
    """
    Lab wall-clock start time for one as-run row.

    Args:    row  TestMatrix_lab_runs.csv dict
    Returns: ``HH:MM`` or empty if unknown
    """
    dt = (row.get("DateTime") or "").strip()
    if dt:
        m = re.match(r"^\d{4}-\d{2}-\d{2}\s+(\d{1,2}):(\d{2})", dt)
        if m:
            return f"{int(m.group(1)):02d}:{m.group(2)}"
    dump = (row.get("DumpFolder") or "").strip()
    m = re.match(r"^r[+-]?\d+_\d{8}_(\d{4})_", dump)
    if m:
        hhmm = m.group(1)
        return f"{hhmm[:2]}:{hhmm[2:]}"
    return ""


def test_key_cells(row: dict[str, str], test_id: str) -> list[str]:
    """
    One footnote-table row for a Test (fixed columns).

    Args:    row  TestMatrix_lab_runs.csv dict; test_id
    Returns: cells matching KEY_HEADERS
    """
    notes: list[str] = []
    t_run = run_time_hhmm(row)
    if t_run:
        notes.append(t_run)
    ch = (row.get("constraintsHandler") or "").strip()
    if ch and ch.lower() != "transformation":
        notes.append(f"constraints={ch}")
    mat_l = (row.get("MatFile") or "").lower()
    note_l = (row.get("Note") or "").lower()
    if is_dry_mat_only(row) or "testbaseline" in mat_l or (
        "dry" in note_l and "baseline" in note_l
    ):
        notes.append("dry")
    if has_hybrid_execute(row):
        notes.append("Hybrid Execute")
    const = (row.get("soilConstitutive") or "").strip()
    if const and const.lower() != "inelastic":
        notes.append(const)
    exp = (row.get("expElementType") or "").strip()
    if exp and exp != "generic":
        notes.append(exp)
    if (row.get("holdPierON") or "").strip() == "0":
        notes.append("noHold")
    xi = (row.get("rayleighXi1") or "").strip()
    if xi and xi not in ("0.03", "0.030"):
        notes.append(f"ξ={xi}")
    return [
        str(test_id),
        mesh_short(row.get("soilMesh", "")),
        mesh_neqn(row),
        soil_short(row.get("soilProfile", "")),
        (row.get("soilEleType") or "").strip() or "?",
        abbreviate_integrator(row.get("eqIntegrator", "")),
        abbreviate_solver(row.get("postPartitionSystem", "")),
        (row.get("Number of Procs") or "?").strip(),
        " ".join(notes),
    ]


def fmt_pct(pct: float) -> str:
    """Format a percentage for an in-bar label (extrapolate / interpolate)."""
    if pct >= 10.0:
        return f"{pct:.0f}"
    if pct >= 1.0:
        return f"{pct:.0f}"
    return f"{pct:.1f}"


def fmt_pct_slowdown(pct: float) -> str:
    """Slowdown % — three decimals so tiny shares stay visible."""
    return f"{pct:.3f}"


def fractions_in_eq_window(
    t_lab_s: np.ndarray,
    state: np.ndarray,
    t0_s: float,
    t1_s: float,
) -> tuple[dict[STAGE_KEY, float], int, float]:
    """
    Stage fractions inside a fixed lab window, plus unfinished time share.

    Recorded stages use sample counts on [t0, min(t_end, t1)]. If the mat
    ends before t1, unfinished = (t1 − t_end) / (t1 − t0) and stage fractions
    are scaled by (1 − unfinished) so the bar sums to 1 over the window.

    Args:    t_lab_s, state, t0_s, t1_s  (lab s)
    Returns: ({typeConv3: fraction}, n_samples in covered window, unfinished)
    """
    empty = _empty_fracs()
    if t_lab_s.size == 0:
        return empty, 0, 1.0
    t_end = float(t_lab_s[-1])
    width = t1_s - t0_s
    if width <= 0.0:
        raise ValueError(f"bad eq window [{t0_s}, {t1_s}]")
    if t_end <= t0_s:
        return empty, 0, 1.0
    unfinished = max(0.0, (t1_s - t_end) / width) if t_end < t1_s else 0.0
    covered = 1.0 - unfinished
    t_hi = min(t_end, t1_s)
    mask = (t_lab_s >= t0_s) & (t_lab_s <= t_hi)
    raw, n = fractions_by_code(state[mask])
    if n == 0:
        return empty, 0, unfinished if unfinished > 0.0 else 1.0
    fracs = {code: raw[code] * covered for code in raw}
    return fracs, n, unfinished


def fractions_full_history(
    state: np.ndarray,
) -> tuple[dict[STAGE_KEY, float], int, float]:
    """
    Sample fractions of each typeConv3 value over the full stateOS record.

    Args:    state
    Returns: ({typeConv3: fraction}, n_samples, unfinished=0)
    """
    fracs, n = fractions_by_code(state)
    return fracs, n, 0.0


def collect_trials(
    fraction_fn,
    *,
    include_dry: bool = False,
) -> tuple[list[TrialFractions], list[str]]:
    """
    Build typeConv3 fractions for every plottable Test row.

    Args:    fraction_fn  (t_lab_s, state) -> (fracs, n, unfinished)
             include_dry  include dry mat-only rows (no DumpFolder)
    Returns: (trials, skip messages)
    """
    run_mat = run_to_mat_name()
    trials: list[TrialFractions] = []
    skips: list[str] = []
    for row in load_lab_rows():
        test_id = (row.get("Test") or "").strip()
        dump = (row.get("DumpFolder") or "").strip()
        mat_name = (row.get("MatFile") or "").strip() or run_mat.get(dump, "")
        if is_dry_mat_only(row) and not include_dry:
            skips.append(f"{test_id}: dry mat-only (use _with_dry)")
            continue
        if not mat_name:
            skips.append(f"{test_id}: no mapped mat ({dump or 'no dump'})")
            continue
        reason = is_excluded_row(row, mat_name)
        if reason:
            skips.append(f"{test_id}: skip ({reason})")
            continue
        pair = load_type_conv3(mat_name)
        if pair is None:
            skips.append(f"{test_id}: no stateOS extract ({mat_name})")
            continue
        t_lab_s, state = pair
        fracs, n, unfinished = fraction_fn(t_lab_s, state)
        if n == 0 and unfinished <= 0.0:
            skips.append(f"{test_id}: empty state slice")
            continue
        trials.append(
            TrialFractions(
                test_id=test_id,
                dump_folder=dump or f"mat:{Path(mat_name).stem}",
                mat_name=mat_name,
                fracs=fracs,
                n_samples=n,
                unfinished=unfinished,
            )
        )
    return trials, skips


def _row_for_trial(
    t: TrialFractions,
    rows_by_dump: dict[str, dict[str, str]],
    rows_by_test: dict[str, dict[str, str]],
) -> dict[str, str]:
    """Resolve CSV row for a trial (dump key or Test ID)."""
    row = rows_by_dump.get(t.dump_folder) or {}
    if row:
        return row
    return rows_by_test.get(t.test_id) or {}


def trial_mesh_name(
    t: TrialFractions,
    rows_by_dump: dict[str, dict[str, str]],
    rows_by_test: dict[str, dict[str, str]] | None = None,
) -> str:
    """
    Mesh label for one trial (Baseline | Moderate | Large | X-Large).

    Args:    t, rows_by_dump, rows_by_test
    Returns: short mesh name
    """
    row = _row_for_trial(t, rows_by_dump, rows_by_test or {})
    return mesh_short(row.get("soilMesh", ""))


def sort_trials_by_mesh(
    trials: list[TrialFractions],
    rows_by_dump: dict[str, dict[str, str]],
    rows_by_test: dict[str, dict[str, str]] | None = None,
) -> list[TrialFractions]:
    """
    Order bars mesh-first (Baseline → X-Large), variants last within each block.

    Args:    trials, rows_by_dump, rows_by_test
    Returns: reordered trials (same objects)
    """
    by_test = rows_by_test or {}
    by_id = {t.test_id: t for t in trials}
    ordered: list[TrialFractions] = []
    seen: set[str] = set()
    for mesh in sorted(MESH_RANK, key=lambda m: MESH_RANK[m]):
        for test_id in BAR_ORDER_BY_MESH.get(mesh, ()):
            t = by_id.get(test_id)
            if t is None:
                continue
            ordered.append(t)
            seen.add(test_id)
    leftovers = [t for t in trials if t.test_id not in seen]
    leftovers.sort(
        key=lambda t: (
            MESH_RANK.get(trial_mesh_name(t, rows_by_dump, by_test), 99),
            t.test_id,
        )
    )
    ordered.extend(leftovers)
    return ordered


def bar_x_positions(
    trials: list[TrialFractions],
    rows_by_dump: dict[str, dict[str, str]],
    rows_by_test: dict[str, dict[str, str]] | None = None,
) -> tuple[np.ndarray, list[float]]:
    """
    Bar centers with gaps between mesh groups; divider x between blocks.

    Args:    trials, rows_by_dump, rows_by_test
    Returns: (x positions, vertical guide positions between mesh blocks)
    """
    by_test = rows_by_test or {}
    x_vals: list[float] = []
    dividers: list[float] = []
    pos = 0.0
    prev_mesh: str | None = None
    for t in trials:
        mesh = trial_mesh_name(t, rows_by_dump, by_test)
        if prev_mesh is not None and mesh != prev_mesh:
            dividers.append(pos - 0.5 * MESH_GROUP_GAP)
            pos += MESH_GROUP_GAP
        x_vals.append(pos)
        pos += 1.0
        prev_mesh = mesh
    return np.asarray(x_vals, dtype=float), dividers


def _label_segment(
    ax,
    x_i: float,
    h: float,
    y0: float,
    *,
    code: STAGE_KEY | None,
    color: str,
    edge: str = "#333333",
) -> None:
    """In-bar % text for one stack segment."""
    if h <= 0:
        return
    is_slow = code in SLOWDOWN_KEYS
    show = (is_slow and h > 0) or h >= PCT_LABEL_MIN
    if not show:
        return
    if is_slow:
        if code == "slowdown_brief":
            y_txt = max(y0 - 0.5, 0.35)
            va = "top"
        else:
            y_txt = y0 + h + 0.5
            va = "bottom"
        txt = ax.text(
            x_i,
            y_txt,
            fmt_pct_slowdown(h),
            ha="center",
            va=va,
            fontsize=FS,
            color=SLOWDOWN_LABEL_COLOR,
            clip_on=True,
            zorder=5,
        )
        txt.set_path_effects(
            [
                patheffects.Stroke(linewidth=2.0, foreground=edge),
                patheffects.Normal(),
            ]
        )
        return
    ax.text(
        x_i,
        y0 + 0.5 * h,
        fmt_pct(h),
        ha="center",
        va="center",
        fontsize=FS,
        color=color,
    )


def plot_stacked_bars(
    trials: list[TrialFractions],
    out: Path,
    *,
    ylabel: str,
) -> None:
    """
    One stacked bar per Test ID; footnote key maps IDs to model/solver.

    Args:    trials, out, ylabel
    Returns: none (writes PNG)
    """
    if not trials:
        print(f"PlotStateOSBars: no trials to plot ({out.name})")
        return

    rows_by_dump = dump_to_row()
    rows_by_test = {r["Test"]: r for r in load_lab_rows()}
    trials = sort_trials_by_mesh(trials, rows_by_dump, rows_by_test)
    labels = [str(t.test_id) for t in trials]
    x, mesh_dividers = bar_x_positions(trials, rows_by_dump, rows_by_test)
    show_unfinished = any(t.unfinished > 0.0 for t in trials)

    n_gaps = len(mesh_dividers)
    fig_w = max(12.0, 0.32 * (len(trials) + n_gaps * MESH_GROUP_GAP) + 2.5)
    # Key table: header + one row per Test; ~0.22 in per row.
    key_h = 0.22 * (len(trials) + 1) + 0.20
    fig_h = 5.8 + key_h
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(
        2,
        1,
        height_ratios=[5.2, key_h],
        hspace=0.08,
        left=0.07,
        right=0.99,
        top=0.90,
        bottom=0.02,
    )
    ax = fig.add_subplot(gs[0])
    ax_key = fig.add_subplot(gs[1])
    ax_key.set_axis_off()

    bottoms = np.zeros(len(trials), dtype=float)
    for code, name, face, edge in STAGES:
        heights = np.array([t.fracs.get(code, 0.0) for t in trials], dtype=float) * 100.0
        lw = 0.75 if code in SLOWDOWN_KEYS else 0.4
        ax.bar(
            x,
            heights,
            bottom=bottoms,
            width=0.78,
            label=name,
            color=face,
            edgecolor=edge,
            linewidth=lw,
        )
        for i, (h, y0) in enumerate(zip(heights, bottoms)):
            _label_segment(ax, x[i], h, y0, code=code, color="white", edge=edge)
        bottoms = bottoms + heights

    if show_unfinished:
        heights_u = np.array([t.unfinished for t in trials], dtype=float) * 100.0
        ax.bar(
            x,
            heights_u,
            bottom=bottoms,
            width=0.78,
            label="unfinished",
            facecolor=UNFINISHED_FACE,
            edgecolor=UNFINISHED_EDGE,
            linewidth=0.5,
            hatch=UNFINISHED_HATCH,
        )

    n_legend = len(STAGES) + (1 if show_unfinished else 0)
    for x_div in mesh_dividers:
        ax.axvline(x_div, color="#cccccc", linewidth=0.8, linestyle="--", zorder=0)
    ax.set_xlim(float(x[0]) - 0.6, float(x[-1]) + 0.4)
    ax.set_ylim(0, 110)  # headroom for sustained-slowdown labels above 100%
    ax.set_ylabel(ylabel, fontsize=FS)
    ax.set_xlabel("Test ID", fontsize=FS)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=FS)
    ax.tick_params(axis="both", labelsize=FS)
    ax.grid(axis="y", alpha=0.25, linewidth=0.6)
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.01),
        ncol=min(n_legend, 5),
        fontsize=FS,
        frameon=False,
        borderaxespad=0.0,
    )

    # Footnote key table: Test → mesh / soil / ele / solver / …
    cell_text: list[list[str]] = []
    for t in trials:
        row = _row_for_trial(t, rows_by_dump, rows_by_test)
        if row:
            cell_text.append(test_key_cells(row, t.test_id))
        else:
            pad = [""] * (len(KEY_HEADERS) - 2)
            cell_text.append([str(t.test_id), t.dump_folder, *pad])

    table = ax_key.table(
        cellText=cell_text,
        colLabels=list(KEY_HEADERS),
        loc="upper center",
        cellLoc="left",
        colLoc="left",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(FS)
    table.scale(1.0, 1.15)
    for (r, c), cell in table.get_celld().items():
        cell.set_linewidth(0.4)
        cell.set_edgecolor("#bdbdbd")
        cell.PAD = 0.02
        # Match chart FS (table cells otherwise pick a larger default).
        if r == 0:
            cell.set_text_props(weight="bold", fontsize=FS)
            cell.set_facecolor("#eeeeee")
        else:
            cell.set_text_props(fontsize=FS)
            cell.set_facecolor("#f7f7f7" if r % 2 == 0 else "#ffffff")
        # Test, Mesh, DOFs, Soil, Ele, Integrator, Solver, np, Notes
        if c in (0, 7):  # Test, np
            cell.set_width(0.045)
        elif c == 2:  # DOFs
            cell.set_width(0.07)
        elif c == 3:  # Soil
            cell.set_width(0.07)
        elif c == 1:  # Mesh
            cell.set_width(0.09)
        elif c == 4:  # Ele
            cell.set_width(0.09)
        elif c == 5:  # Integrator (mathtext)
            cell.set_width(0.20)
        elif c == 6:  # Solver
            cell.set_width(0.12)
        elif c == 8:  # Notes
            cell.set_width(0.18)

    # Nudge key table slightly below the gridspec slot (extra gap under the bars).
    fig.canvas.draw()
    box = ax_key.get_position()
    ax_key.set_position([box.x0, box.y0 - 0.06 * box.height, box.width, box.height])

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=DPI, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print(f"PlotStateOSBars: wrote {out}  ({len(trials)} tests)")


def main() -> int:
    """Build stacked typeConv3 bar charts (eq 35–90 s lab window and full)."""
    out_dir = stateos_compare_plots_dir()
    t0, t1 = EQ_WINDOW_LAB_S

    def eq_fn(t_lab_s: np.ndarray, state: np.ndarray) -> tuple[dict[STAGE_KEY, float], int, float]:
        return fractions_in_eq_window(t_lab_s, state, t0, t1)

    def full_fn(_t_lab_s: np.ndarray, state: np.ndarray) -> tuple[dict[STAGE_KEY, float], int, float]:
        return fractions_full_history(state)

    def run_suite(*, include_dry: bool, day: str, day_label: str) -> bool:
        trials_eq, skips = collect_trials(eq_fn, include_dry=include_dry)
        trials_eq = filter_trials_by_day(trials_eq, day)
        print(f"\n--- stateOS bars {day_label} ---")
        for msg in skips:
            if include_dry or "dry mat-only" not in msg:
                print(f"  skip: {msg}")
        dry_suffix = "_with_dry" if include_dry else ""
        plot_stacked_bars(
            trials_eq,
            out_dir / f"bar_typeconv3_eq35_90_{day}{dry_suffix}.png",
            ylabel=f"fraction of lab window [{t0:.0f}, {t1:.0f}] s (%)",
        )
        trials_full, skips_full = collect_trials(full_fn, include_dry=include_dry)
        trials_full = filter_trials_by_day(trials_full, day)
        for msg in skips_full:
            if msg not in skips and (include_dry or "dry mat-only" not in msg):
                print(f"  skip: {msg}")
        plot_stacked_bars(
            trials_full,
            out_dir / f"bar_typeconv3_full_{day}{dry_suffix}.png",
            ylabel="fraction of full-record samples (%)",
        )
        return bool(trials_eq or trials_full)

    ok = False
    for day, day_label in CAMPAIGN_DAYS:
        ok_main = run_suite(
            include_dry=False,
            day=day,
            day_label=f"{day_label} (wet + dry+OS baseline)",
        )
        ok_dry = run_suite(
            include_dry=True,
            day=day,
            day_label=f"{day_label} with dry mat-only",
        )
        ok = ok or ok_main or ok_dry

    for obsolete in (
        "bar_typeconv3_eq35_90.png",
        "bar_typeconv3_full.png",
        "bar_typeconv3_eq35_90_with_dry.png",
        "bar_typeconv3_full_with_dry.png",
        "bar_typeconv3_gm_d595.png",
    ):
        old = out_dir / obsolete
        if old.is_file():
            old.unlink()
            print(f"PlotStateOSBars: removed obsolete {old.name}")

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
