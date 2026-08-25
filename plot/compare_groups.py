#!/usr/bin/env python3
"""Group as-run lab dumps by *model* knobs for fair compare overlays.

Solver / np / integrator changes stay in the same group (same mesh, soil,
element, constitutive, hybrid setup). Physics changes get their own folder.

Reads: plot/opensees_data/TestMatrix_lab_runs.csv (fallback: LOCAL copy).
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


def _slug_part(s: str) -> str:
    s = (s or "").strip()
    s = re.sub(r"\s*\([^)]*\)\s*", "", s)  # drop " (SOFT)" style tags
    s = s.replace("%", "pct")
    s = re.sub(r"[^A-Za-z0-9.+-]+", "_", s)
    return s.strip("_") or "x"


def group_slug(row: dict[str, str]) -> str:
    """Folder name under LOCAL/plots/compare/."""
    mesh = _slug_part(row.get("soilMesh", ""))
    # Prefer short mesh tags: 0, 1, 2, 3
    m0 = re.match(r"^(\d+)", mesh)
    mesh_tag = f"mesh{m0.group(1)}" if m0 else f"mesh_{mesh}"
    parts = [
        mesh_tag,
        _slug_part(row.get("soilProfile", "")),
        _slug_part(row.get("soilEleType", "")),
        _slug_part(row.get("soilConstitutive", "")),
        _slug_part(row.get("expElementType", "")),
    ]
    hold = (row.get("holdPierON") or "").strip()
    if hold == "0":
        parts.append("noHold")
    xi = (row.get("rayleighXi1") or "").strip()
    if xi and xi not in ("0.03", "0.030"):
        parts.append(f"xi{_slug_part(xi)}")
    return "_".join(parts)


def group_label(row: dict[str, str]) -> str:
    """Short human label for logs / optional titles."""
    bits = [
        row.get("soilMesh", ""),
        row.get("soilProfile", ""),
        row.get("soilEleType", ""),
        row.get("soilConstitutive", ""),
        row.get("expElementType", ""),
    ]
    hold = (row.get("holdPierON") or "").strip()
    if hold == "0":
        bits.append("hold=0")
    xi = (row.get("rayleighXi1") or "").strip()
    if xi and xi not in ("0.03", "0.030"):
        bits.append(f"xi={xi}")
    return ", ".join(b for b in bits if b)


def load_lab_rows(path: Path | None = None) -> list[dict[str, str]]:
    p = path or lab_runs_csv_path()
    if not p.is_file():
        return []
    with p.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def groups_by_dump(
    path: Path | None = None,
) -> OrderedDict[str, list[dict[str, str]]]:
    """slug → list of lab_runs rows (same DumpFolder order as CSV)."""
    out: OrderedDict[str, list[dict[str, str]]] = OrderedDict()
    for row in load_lab_rows(path):
        dump = (row.get("DumpFolder") or "").strip()
        if not dump:
            continue
        slug = group_slug(row)
        out.setdefault(slug, []).append(row)
    return out


def dump_to_group(path: Path | None = None) -> dict[str, str]:
    """DumpFolder name → compare group slug."""
    m: dict[str, str] = {}
    for slug, rows in groups_by_dump(path).items():
        for r in rows:
            m[r["DumpFolder"]] = slug
    return m
