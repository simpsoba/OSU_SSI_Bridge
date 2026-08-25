#!/usr/bin/env python3
"""
Goals
-----
Bucket as-run lab dumps by *model* knobs so PlotEQCompareRuns / PlotMatOS
put fair overlays in one folder.

  Same group  = same mesh, soil profile/element/constitutive, exp element,
                holdPier, and (if non-default) Rayleigh ξ1.
  Not in key  = solver, np, integrator, precision, hybridExecuteMode, …

Reads: plot/opensees_data/TestMatrix_lab_runs.csv
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


def group_slug(row: dict[str, str]) -> str:
    """
    Folder name under LOCAL/plots/compare/.

    Args:    row  one TestMatrix_lab_runs.csv dict
    Returns: e.g. mesh0_4_SSPQuad_Inelastic_generic
    """
    mesh = _slug_part(row.get("soilMesh", ""))
    mesh_match = re.match(r"^(\d+)", mesh)
    mesh_tag = f"mesh{mesh_match.group(1)}" if mesh_match else f"mesh_{mesh}"

    parts = [
        mesh_tag,
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


def group_label(row: dict[str, str]) -> str:
    """
    Short human string for logs (not the folder name).

    Args:    row
    Returns: comma-separated knob summary
    """
    bits = [
        row.get("soilMesh", ""),
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
