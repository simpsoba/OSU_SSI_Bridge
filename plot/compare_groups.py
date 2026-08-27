#!/usr/bin/env python3
"""
Goals
-----
Bucket as-run lab dumps for fair overlays under LOCAL/plots/compare/.

  Campaign folders (mesh size):
    Baseline / Moderate / Large / X-Large
      /<profile>_<ele>_<const>_<exp>[_noHold][_xi…]/pairs/

  Excluded (testBaseline, twoNodeLink, single-precision):
    _excluded/<reason>/<same variant>/pairs/

  Same compare set = same mesh folder + same variant knobs.
  Not in key = solver, np, integrator, hybridExecuteMode, …

Reads: plot/lab/TestMatrix_lab_runs.csv
       (fallback: OSU_SSI_BRIDGE_DATA_LOCAL/TestMatrix_lab_runs.csv)
"""

from __future__ import annotations

import csv
import re
from collections import OrderedDict
from pathlib import Path

from lab_paths import lab_runs_csv_path

# Columns that define the structural / soil model (not numerics).
MODEL_KEYS = (
    "soilMesh",
    "soilProfile",
    "soilEleType",
    "soilConstitutive",
    "expElementType",
    "holdPierON",
    "rayleighXi1",
)

# Default campaign ξ1 — only append to the slug when different.
DEFAULT_XI1 = ("0.03", "0.030")

# Mesh index → compare top folder (matches PlotStateOSBars mesh_short).
MESH_FOLDER: dict[str, str] = {
    "0": "Baseline",
    "1": "Moderate",
    "2": "Large",
    "3": "X-Large",
}

SKIP_FOLDER: dict[str, str] = {
    "twoNodeLink": "twoNodeLink",
    "single-precision CuDSS": "singlePrecision",
}

# r±NN_YYYYMMDD_HHMM_…
_FOLDER_RE = re.compile(r"^r([+-]?\d+)_(\d{8})_(\d{4})_", re.IGNORECASE)


# ------------------------------------------------------------
# slug pieces
# ------------------------------------------------------------


def _slug_part(text: str) -> str:
    """
    One path-safe token from a CSV cell (drop parenthetical tags).

    Args:    text  e.g. "4 (SOFT)" or "SSPQuad"
    Returns: "4" or "SSPQuad"
    """
    s = (text or "").strip()
    s = re.sub(r"\s*\([^)]*\)\s*", "", s)
    s = s.replace("%", "pct")
    s = re.sub(r"[^A-Za-z0-9.+-]+", "_", s)
    return s.strip("_") or "x"


def mesh_folder_name(row: dict[str, str]) -> str:
    """
    Top compare folder from soilMesh (Baseline / Moderate / …).

    Args:    row  lab_runs CSV dict
    Returns: folder name
    """
    mesh = _slug_part(row.get("soilMesh", ""))
    mesh_match = re.match(r"^(\d+)", mesh)
    if not mesh_match:
        return f"mesh_{mesh}" if mesh else "mesh_x"
    return MESH_FOLDER.get(mesh_match.group(1), f"mesh{mesh_match.group(1)}")


def variant_slug(row: dict[str, str]) -> str:
    """
    Model-knob subfolder under a mesh (or under _excluded/<reason>/).

    Args:    row
    Returns: e.g. 4_SSPQuad_Inelastic_generic
    """
    parts = [
        _slug_part(row.get("soilProfile", "")),
        _slug_part(row.get("soilEleType", "")),
        _slug_part(row.get("soilConstitutive", "")),
        _slug_part(row.get("expElementType", "")),
    ]
    if (row.get("holdPierON") or "").strip() == "0":
        parts.append("noHold")
    xi = (row.get("rayleighXi1") or "").strip()
    if xi and xi not in DEFAULT_XI1:
        parts.append(f"xi{_slug_part(xi)}")
    return "_".join(parts)


def analysis_skip_reason(
    row: dict[str, str] | None,
    mat_name: str = "",
) -> str | None:
    """
    Skip reason for campaign bars / pairwise compare, or None if OK.

    Excludes single-precision CuDSS and twoNodeLink. Dry OpenFresco
    baseline (testBaseline / F01) is kept — useful reference without
    wet flume motion. Dry mat-only rows (no DumpFolder) are gated by
    PlotStateOSBars ``include_dry``, not here.

    Args:    row  TestMatrix_lab_runs dict (may be None); mat_name
    Returns: reason string or None
    """
    row = row or {}
    exp = (row.get("expElementType") or "").strip().lower()
    if exp == "twonodelink":
        return "twoNodeLink"
    blob = " ".join(
        [
            row.get("Goal", ""),
            row.get("postPartitionSystem", ""),
            row.get("eqIntegrator", ""),
            mat_name or "",
        ]
    ).lower()
    if "-precision" in blob or "single precision" in blob:
        return "single-precision CuDSS"
    return None


def is_dry_mat_only(row: dict[str, str]) -> bool:
    """True when the row has a mat but no OpenSees dump folder."""
    return bool((row.get("MatFile") or "").strip()) and not (
        row.get("DumpFolder") or ""
    ).strip()


def group_slug(row: dict[str, str], mat_name: str | None = None) -> str:
    """
    Relative path under LOCAL/plots/compare/.

    Campaign: Baseline/4_SSPQuad_Inelastic_generic
    Excluded: _excluded/twoNodeLink/4_SSPQuad_Inelastic_twoNodeLink_noHold

    Args:    row; mat_name  optional (defaults to row MatFile)
    Returns: slash-separated relative path
    """
    mat = mat_name if mat_name is not None else (row.get("MatFile") or "")
    variant = variant_slug(row)
    reason = analysis_skip_reason(row, mat)
    if reason:
        folder = SKIP_FOLDER.get(reason, _slug_part(reason))
        return f"_excluded/{folder}/{variant}"
    return f"{mesh_folder_name(row)}/{variant}"


def group_label(row: dict[str, str]) -> str:
    """
    Short human string for logs (not the folder name).

    Args:    row
    Returns: comma-separated knob summary
    """
    bits = [
        mesh_folder_name(row),
        row.get("soilProfile", ""),
        row.get("soilEleType", ""),
        row.get("soilConstitutive", ""),
        row.get("expElementType", ""),
    ]
    if (row.get("holdPierON") or "").strip() == "0":
        bits.append("hold=0")
    xi = (row.get("rayleighXi1") or "").strip()
    if xi and xi not in DEFAULT_XI1:
        bits.append(f"xi={xi}")
    reason = analysis_skip_reason(row, row.get("MatFile", ""))
    if reason:
        bits.append(f"excluded={SKIP_FOLDER.get(reason, reason)}")
    return ", ".join(b for b in bits if b)


# ------------------------------------------------------------
# CSV → groups
# ------------------------------------------------------------


def load_lab_rows(path: Path | None = None) -> list[dict[str, str]]:
    """
    Read the curated as-run matrix.

    Args:    path  optional override (default: lab_runs_csv_path())
    Returns: list of row dicts (empty if missing)
    """
    csv_path = path or lab_runs_csv_path()
    if not csv_path.is_file():
        return []
    with csv_path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def groups_by_dump(
    path: Path | None = None,
) -> OrderedDict[str, list[dict[str, str]]]:
    """
    Compare-group slug → lab_runs rows (CSV order within each group).

    Args:    path  optional CSV override
    Returns: OrderedDict
    """
    out: OrderedDict[str, list[dict[str, str]]] = OrderedDict()
    for row in load_lab_rows(path):
        dump = (row.get("DumpFolder") or "").strip()
        if not dump:
            continue
        slug = group_slug(row)
        out.setdefault(slug, []).append(row)
    return out


def dump_to_group(path: Path | None = None) -> dict[str, str]:
    """
    DumpFolder name → compare group slug.

    Args:    path  optional CSV override
    Returns: {dump_folder: slug}
    """
    mapping: dict[str, str] = {}
    for slug, rows in groups_by_dump(path).items():
        for row in rows:
            mapping[row["DumpFolder"]] = slug
    return mapping


def dump_to_row(path: Path | None = None) -> dict[str, dict[str, str]]:
    """
    DumpFolder name → lab_runs CSV row.

    Args:    path  optional CSV override
    Returns: {dump_folder: row}
    """
    return {
        row["DumpFolder"]: row
        for row in load_lab_rows(path)
        if (row.get("DumpFolder") or "").strip()
    }


def dump_to_test_id(path: Path | None = None) -> dict[str, str]:
    """
    DumpFolder → campaign Test ID (W## / F## / Fd## / Fx##).

    Args:    path  optional CSV override
    Returns: {dump_folder: test_id}
    """
    out: dict[str, str] = {}
    for row in load_lab_rows(path):
        dump = (row.get("DumpFolder") or "").strip()
        test_s = (row.get("Test") or "").strip()
        if not dump or not test_s:
            continue
        out[dump] = test_s
    return out


def test_file_slug(test_id: str | int) -> str:
    """
    File-safe Test ID tag for pair PNGs and plots/runs/.

    Args:    test_id  CSV Test cell (F27) or legacy int
    Returns: e.g. F27 (or T25 for legacy ints)
    """
    s = str(test_id).strip()
    if re.match(r"^(W|F|Fd|Fx)\d+$", s, flags=re.IGNORECASE):
        # Normalize letter prefix case: Fd / Fx / F / W + zero-padded digits as given
        m = re.match(r"^(W|F|Fd|Fx)(\d+)$", s, flags=re.IGNORECASE)
        assert m is not None
        prefix = m.group(1)
        prefix = {"w": "W", "f": "F", "fd": "Fd", "fx": "Fx"}[prefix.lower()]
        return f"{prefix}{m.group(2)}"
    try:
        return f"T{int(s)}"
    except ValueError:
        return s or "Tx"

def row_tag(run_cell: str) -> str:
    """
    Zero-padded row tag: r-08, r+01, r-120.

    Args:    run_cell  TestMatrix Run column
    Returns: legend row token
    """
    n = int(str(run_cell).strip())
    if abs(n) >= 100:
        return f"r+{n}" if n > 0 else f"r{n}"
    if n > 0:
        return f"r+{n:02d}"
    if n == 0:
        return "r+00"
    return f"r-{abs(n):02d}"


def abbreviate_integrator(text: str) -> str:
    """
    Shorten eqIntegrator for legends / tables (mathtext).

    MKRAlphaExplicitMultiSOE 0.5 → MKR-$\\alpha$ $\\rho_{\\infty}=0.5$
    CudaMKRAlpha 0.5 → CUDA-MKR-$\\alpha$ $\\rho_{\\infty}=0.5$
    """
    s = (text or "").strip()
    s_clean = re.sub(r"\s*-incrementalAccel\b", "", s, flags=re.IGNORECASE)
    m = re.match(
        r"^(Cuda)?MKRAlpha(?:ExplicitMultiSOE)?\s+([0-9.]+)",
        s_clean,
        flags=re.IGNORECASE,
    )
    if m:
        rho = m.group(2)
        if m.group(1):
            return rf"CUDA-MKR-$\alpha$ $\rho_{{\infty}}={rho}$"
        return rf"MKR-$\alpha$ $\rho_{{\infty}}={rho}$"
    return re.sub(r"\s+", " ", s_clean).strip()


def abbreviate_solver(text: str) -> str:
    """
    Shorten postPartitionSystem for legends.

    DistributedCuDSS → CuDSS; ParallelProfileSPD → ProfileSPD.
    Strips ``-hybridExecuteMode N`` (shown as Hybrid Execute in Notes).
    """
    s = (text or "").strip()
    s = re.sub(r"\s*-hybridExecuteMode\s*\d+\b", "", s, flags=re.IGNORECASE)
    s = s.replace("DistributedCuDSS", "CuDSS")
    s = s.replace("ParallelProfileSPD", "ProfileSPD")
    return re.sub(r"\s+", " ", s).strip()


def has_hybrid_execute(row: dict[str, str]) -> bool:
    """True when CuDSS hybridExecuteMode 1 (or Goal says hybrid execute)."""
    sol = (row.get("postPartitionSystem") or "") + " " + (row.get("Goal") or "")
    return bool(
        re.search(r"hybridExecuteMode\s*1", sol, flags=re.IGNORECASE)
        or re.search(r"\bhybrid\s+execute\b", sol, flags=re.IGNORECASE)
    )


# Fallback when a row has no DOFs cell: OpenSees systemSize after gravity
# (Run.tcl / RunParallel.tcl "Model size … DOFs"), serial probe 2026-08-27.
# Profile 4, Shin, SSPquad/quad, inelastic, lumpedPlasticity + disp piles.
# (quad vs SSPquad and expElementType twoNodeLink do not change this count.)
MESH_NEQN: dict[int, int] = {
    -2: 1200,  # still estimate (not in campaign; not re-probed)
    -1: 1500,  # still estimate
    0: 2028,
    1: 3216,
    2: 4620,
    3: 5484,
    4: 5900,  # still estimate
}


def mesh_neqn(row: dict[str, str]) -> str:
    """
    DOF count for footnote tables.

    Prefer the as-run ``DOFs`` column (OpenSees ``systemSize``). Fall back to
    ``MESH_NEQN`` from ``soilMesh`` when the cell is blank.
    """
    raw = (row.get("DOFs") or "").strip().replace(",", "")
    if raw.isdigit():
        return f"{int(raw):,}"
    mesh = (row.get("soilMesh") or "").strip()
    m = re.match(r"^(-?\d+)", mesh)
    if not m:
        return "?"
    n = MESH_NEQN.get(int(m.group(1)))
    return f"{n:,}" if n is not None else "?"


def format_run_legend(
    row: dict[str, str],
    *,
    include_hhmm: bool = False,
    incomplete: bool = False,
) -> str:
    """
    Legend line: F## [| HH:MM] | integrator | solver | [np=N][*].

    Args:    row  lab_runs CSV dict
             include_hhmm  append wall-clock when labels collide
             incomplete  append *
    Returns: legend string
    """
    test_s = (row.get("Test") or "").strip()
    tag = test_file_slug(test_s) if test_s else row_tag(row.get("Run", "0"))
    dump = (row.get("DumpFolder") or "").strip()
    match = _FOLDER_RE.match(dump)
    if include_hhmm and match:
        hhmm = match.group(3)
        tag = f"{tag} {hhmm[:2]}:{hhmm[2:]}"
    integ = abbreviate_integrator(row.get("eqIntegrator", ""))
    sol = abbreviate_solver(row.get("postPartitionSystem", ""))
    np_procs = (row.get("Number of Procs") or "?").strip()
    label = f"{tag} | {integ} | {sol} | [np={np_procs}]"
    if incomplete:
        label += "*"
    return label


def legend_labels_for_dumps(
    dump_names: list[str],
    *,
    incomplete: dict[str, bool] | None = None,
    path: Path | None = None,
) -> dict[str, str]:
    """
    Build legend labels for a set of dumps; add HH:MM when bases collide.

    Args:    dump_names  DumpFolder names in plot order
             incomplete  dump → append *
             path  optional CSV override
    Returns: {dump_name: legend}
    """
    rows = dump_to_row(path)
    base: dict[str, str] = {}
    for name in dump_names:
        row = rows.get(name)
        if row is None:
            base[name] = name
        else:
            base[name] = format_run_legend(row, include_hhmm=False)

    counts: dict[str, int] = {}
    for label in base.values():
        counts[label] = counts.get(label, 0) + 1
    need_hhmm = {name for name, label in base.items() if counts[label] > 1}

    flags = incomplete or {}
    out: dict[str, str] = {}
    for name in dump_names:
        row = rows.get(name)
        if row is None:
            out[name] = name + ("*" if flags.get(name) else "")
            continue
        out[name] = format_run_legend(
            row,
            include_hhmm=name in need_hhmm,
            incomplete=bool(flags.get(name)),
        )
    return out
