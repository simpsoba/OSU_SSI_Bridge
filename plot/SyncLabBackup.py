#!/usr/bin/env python3
"""
Goals
-----
Ingest lab dumps into the local working mirror, then archive to Shared Drive.

  python plot/SyncLabBackup.py
  python plot/SyncLabBackup.py --extract-only
  python plot/SyncLabBackup.py --no-extract
  python plot/SyncLabBackup.py --force-extract
  python plot/SyncLabBackup.py --no-archive

Flow
----
1. Read new/changed files from the live upload folder
   ``G:/.shortcut-targets-by-id/…/opensees data/`` (ingest only).
2. Copy into ``OSU_SSI_BRIDGE_DATA_LOCAL/opensees_data/`` (working mirror).
3. Copy LOCAL dumps into the Shared Drive safety archive under date folders
   ``…/2026-OSU-SSI-Bridge/<YYYY-MM-DD>/opensees_data/`` (skips ``mat_extract/``,
   ``plots/``).
4. Extract Simulink ``data`` + *OS channels into LOCAL ``mat_extract/``.

Post-process (PlotEQ*, PlotMatOS, compare) must use the LOCAL mirror only —
never the ingest shortcut or Shared Drive. Skip every ``plots/`` directory so
generated figures never enter dump trees.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import scipy.io as sio

from lab_paths import (
    ARCHIVE_SKIP_TOP,
    LOCAL_MAT_RUN_MAP_PATH,
    LOCAL_OPENSEES_DATA,
    LOCAL_ROOT,
    MANIFEST_PATH,
    MAT_EXTRACT_DIR,
    MAT_KEYS,
    LAB_RUNS_CSV,
    MAT_RUN_MAP_PATH,
    SHARED_DRIVE_ARCHIVE,
    archive_day_from_name,
    archive_opensees_data_for_day,
    cell_str,
    migrate_legacy_local_layout,
    migrate_legacy_plots,
    migrate_plots_layout,
    resolve_lab_ingest_source,
    resolve_simulink_dir,
    signal_names,
    time_column,
)


# ------------------------------------------------------------
# 1. MANIFEST AND FILE DISCOVERY
# ------------------------------------------------------------


def utc_now() -> str:
    """
    Current UTC time in the manifest timestamp format.

    Args:    none
    Returns: UTC timestamp as YYYY-MM-DDTHH:MM:SSZ
    """
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def file_fingerprint(path: Path) -> dict[str, int]:
    """
    Size and modification time used to detect a changed source file.

    Args:    path
    Returns: {"size": bytes, "mtime_ns": nanoseconds}
    """
    stat = path.stat()
    return {"size": int(stat.st_size), "mtime_ns": int(stat.st_mtime_ns)}


def load_manifest() -> dict:
    """
    Read the LOCAL backup manifest.

    Args:    none
    Returns: manifest dict, or an empty manifest when the file is missing
    """
    if not MANIFEST_PATH.is_file():
        return {"files": {}, "updated": None}
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


def save_manifest(manifest: dict) -> None:
    """
    Write the LOCAL backup manifest with a fresh UTC timestamp.

    Args:    manifest
    Returns: none
    """
    MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    manifest["updated"] = utc_now()
    MANIFEST_PATH.write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )


def iter_source_files(
    source_root: Path,
    *,
    skip_top: frozenset[str] | None = None,
) -> list[Path]:
    """
    Files below a source root, excluding ``plots/`` and optional top dirs.

    Args:    source_root  ingest, LOCAL, or archive folder
             skip_top     extra top-level names to skip (e.g. mat_extract)
    Returns: file paths in pathlib traversal order
    """
    skip_top = skip_top or frozenset()
    source_files: list[Path] = []
    for path in source_root.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(source_root)
        if rel.parts and rel.parts[0] in skip_top:
            continue
        relative_parts = {part.lower() for part in rel.parts}
        if "plots" in relative_parts:
            continue
        source_files.append(path)
    return source_files


def needs_copy(
    source_path: Path,
    destination_path: Path,
    old_fingerprint: dict | None,
) -> bool:
    """
    Whether a source file is new or changed since the last sync.

    Args:    source_path, destination_path, old_fingerprint
    Returns: True when the file must be copied
    """
    if not destination_path.is_file():
        return True
    source_fingerprint = file_fingerprint(source_path)
    if old_fingerprint is None:
        return source_fingerprint != file_fingerprint(destination_path)
    return (
        source_fingerprint["size"] != old_fingerprint.get("size")
        or source_fingerprint["mtime_ns"] != old_fingerprint.get("mtime_ns")
    )


# ------------------------------------------------------------
# 2. TREE SYNC (ingest → LOCAL → archive)
# ------------------------------------------------------------


def sync_tree(
    source_root: Path,
    destination_root: Path,
    manifest: dict,
    *,
    manifest_key: str = "files",
    skip_top: frozenset[str] | None = None,
    label: str = "copy",
) -> tuple[int, int, list[Path]]:
    """
    Copy new or changed files from source into destination.

    Args:    source_root, destination_root, manifest
             manifest_key  dict key under manifest for fingerprints
             skip_top      top-level dirs to omit
             label         log prefix
    Returns: (copied count, skipped count, copied Simulink .mat paths)
    """
    destination_root.mkdir(parents=True, exist_ok=True)
    manifest_files = manifest.setdefault(manifest_key, {})
    copied_count = 0
    skipped_count = 0
    copied_mat_paths: list[Path] = []

    for source_path in iter_source_files(source_root, skip_top=skip_top):
        relative_path = source_path.relative_to(source_root).as_posix()
        destination_path = destination_root / relative_path
        old_fingerprint = manifest_files.get(relative_path)

        if not needs_copy(source_path, destination_path, old_fingerprint):
            skipped_count += 1
            continue

        destination_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_path, destination_path)
        manifest_files[relative_path] = file_fingerprint(source_path)
        copied_count += 1

        if source_path.suffix.lower() == ".mat" and "Simulink" in source_path.parts:
            copied_mat_paths.append(destination_path)
        print(f"  {label} {relative_path}")

    return copied_count, skipped_count, copied_mat_paths


def sync_archive_by_day(
    local_root: Path,
    manifest: dict,
) -> tuple[int, int]:
    """
    Copy LOCAL opensees_data into ``…/2026-OSU-SSI-Bridge/<day>/opensees_data/``.

    Day is inferred from the dump folder name or Simulink mat stem
    (``0819_…`` → 2026-08-19, ``r-*_20260821_…`` / ``0821_…`` → 2026-08-21).

    Args:    local_root  LOCAL opensees_data; manifest
    Returns: (copied count, skipped count)
    """
    manifest_files = manifest.setdefault("archive_files", {})
    copied_count = 0
    skipped_count = 0

    for source_path in iter_source_files(local_root, skip_top=ARCHIVE_SKIP_TOP):
        relative = source_path.relative_to(local_root)
        parts = relative.parts
        if parts[0] == "Simulink" and len(parts) > 1:
            day = archive_day_from_name(parts[1])
        else:
            day = archive_day_from_name(parts[0])
        destination_root = archive_opensees_data_for_day(day)
        destination_path = destination_root / relative
        # Manifest key includes day so the same relative path can exist twice.
        manifest_key = f"{day}/{relative.as_posix()}"
        old_fingerprint = manifest_files.get(manifest_key)

        if not needs_copy(source_path, destination_path, old_fingerprint):
            skipped_count += 1
            continue

        destination_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_path, destination_path)
        manifest_files[manifest_key] = file_fingerprint(source_path)
        copied_count += 1
        print(f"  archive {day}/{relative.as_posix()}")

    return copied_count, skipped_count


# ------------------------------------------------------------
# 3. MATLAB STRUCT HELPERS
# ------------------------------------------------------------


def unwrap_struct(value: Any) -> Any:
    """
    Unwrap a one-item MATLAB object array.

    Args:    value  scipy.io MATLAB value
    Returns: contained value, or the original value
    """
    if isinstance(value, np.ndarray) and value.dtype == object and value.size == 1:
        return value.flat[0]
    return value


def extract_os_block(block: Any) -> dict:
    """
    Read signal names and numeric data from one OpenFresco *OS struct.

    Args:    block  MATLAB struct with signalNames and data
    Returns: {"signalNames": list[str], "data": float array}
    """
    names = signal_names(block)
    data = np.asarray(block.data, dtype=float)
    return {"signalNames": names, "data": data}


def extract_data_container(block: Any) -> dict:
    """
    Keep `data` metadata and nested *OS arrays, excluding topL and AFC.

    Args:    block  MATLAB `data` container
    Returns: extracted metadata, nested *OS blocks, and value shapes
    """
    extracted: dict = {
        "fileName": [cell_str(value) for value in block.fileName.flatten()],
        "sigNames": [cell_str(value) for value in block.sigNames.flatten()],
    }

    for key in ("tarSigOS", "comSigOS", "meaSigOS", "stateOS"):
        if not hasattr(block, key):
            continue
        nested_block = unwrap_struct(getattr(block, key))
        if hasattr(nested_block, "signalNames") and hasattr(nested_block, "data"):
            extracted[key] = extract_os_block(nested_block)

    if hasattr(block, "values"):
        value_shapes: list[list[int] | None] = []
        for value in block.values.flatten():
            unwrapped_value = unwrap_struct(value)
            if hasattr(unwrapped_value, "data"):
                value_shapes.append(list(np.asarray(unwrapped_value.data).shape))
            elif isinstance(unwrapped_value, np.ndarray):
                value_shapes.append(list(unwrapped_value.shape))
            else:
                value_shapes.append(None)
        extracted["values_shapes"] = value_shapes

    return extracted


# ------------------------------------------------------------
# 4. MAT-TO-NPZ EXTRACT
# ------------------------------------------------------------


def extract_mat(mat_path: Path, output_npz: Path) -> dict:
    """
    Extract MAT_KEYS channels from one Simulink file into compressed NPZ.

    Args:    mat_path, output_npz
    Returns: summary dict with source, output, keys, target rows, and last time (s)
    """
    print(f"  extract {mat_path.name}")
    raw_mat = sio.loadmat(
        str(mat_path),
        variable_names=list(MAT_KEYS),
        squeeze_me=False,
        struct_as_record=False,
    )
    payload: dict = {"mat": mat_path.name, "extracted": utc_now()}
    metadata: dict = {}

    for key in MAT_KEYS:
        if key not in raw_mat:
            continue
        block = unwrap_struct(raw_mat[key])

        if key == "data":
            data_container = extract_data_container(block)
            metadata["data_fileName"] = np.array(
                data_container["fileName"],
                dtype=object,
            )
            metadata["data_sigNames"] = np.array(
                data_container["sigNames"],
                dtype=object,
            )
            if "values_shapes" in data_container:
                metadata["data_values_shapes"] = np.array(
                    data_container["values_shapes"],
                    dtype=object,
                )
            for os_key in ("tarSigOS", "comSigOS", "meaSigOS", "stateOS"):
                if os_key not in data_container:
                    continue
                payload[f"data_{os_key}_data"] = data_container[os_key]["data"]
                metadata[f"data_{os_key}_signalNames"] = np.array(
                    data_container[os_key]["signalNames"],
                    dtype=object,
                )
            continue

        if hasattr(block, "signalNames") and hasattr(block, "data"):
            os_block = extract_os_block(block)
            payload[f"{key}_data"] = os_block["data"]
            metadata[f"{key}_signalNames"] = np.array(
                os_block["signalNames"],
                dtype=object,
            )

            names = os_block["signalNames"]
            try:
                time_index = time_column(names)
                payload[f"{key}_time"] = os_block["data"][:, time_index]
                for signal_index, name in enumerate(names):
                    if signal_index != time_index:
                        payload[f"{key}_primary"] = os_block["data"][:, signal_index]
                        metadata[f"{key}_primaryName"] = np.array(
                            [name],
                            dtype=object,
                        )
                        break
            except KeyError:
                pass

    output_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_npz, **payload, **metadata)

    return {
        "mat": mat_path.name,
        "npz": output_npz.name,
        "keys": sorted(
            key for key in payload if key.endswith("_data") or key == "mat"
        ),
        "n_tar": (
            int(payload["tarSigOS_data"].shape[0])
            if "tarSigOS_data" in payload
            else 0
        ),
        "t_real_last": (
            float(payload["tarSigOS_time"][-1])
            if "tarSigOS_time" in payload
            else None
        ),
    }


def extract_mats(mat_paths: list[Path], force: bool = False) -> list[dict]:
    """
    Extract MAT files whose NPZ output is missing or older.

    Args:    mat_paths; force  re-extract even when output is current
    Returns: summary dicts for files extracted during this call
    """
    MAT_EXTRACT_DIR.mkdir(parents=True, exist_ok=True)
    summaries: list[dict] = []

    for mat_path in mat_paths:
        output_npz = MAT_EXTRACT_DIR / f"{mat_path.stem}.npz"
        if output_npz.is_file() and not force:
            source_mtime_ns = mat_path.stat().st_mtime_ns
            output_mtime_ns = output_npz.stat().st_mtime_ns
            if output_mtime_ns >= source_mtime_ns:
                print(f"  skip extract (up to date) {mat_path.name}")
                continue
        summaries.append(extract_mat(mat_path, output_npz))

    return summaries


def list_local_or_drive_mats(source_root: Path) -> list[Path]:
    """
    List Simulink MAT files under one opensees_data root.

    Args:    source_root
    Returns: sorted .mat paths, or an empty list if Simulink/ is missing
    """
    simulink_dir = resolve_simulink_dir(source_root)
    if simulink_dir is None:
        return []
    return sorted(simulink_dir.glob("*.mat"))


# ------------------------------------------------------------
# 5. MAT-TO-RUN MAP REPORT
# ------------------------------------------------------------


def report_map_coverage(source_root: Path) -> None:
    """
    Print unmapped MAT files, mapped dumps, and runs without MAT files.

    Args:    source_root  opensees_data folder used to check dump presence
    Returns: none
    """
    from lab_paths import build_mat_run_catalog, load_mat_orphans

    catalog = build_mat_run_catalog()
    orphans = load_mat_orphans()
    mat_names = {path.name for path in list_local_or_drive_mats(source_root)}
    known = set(catalog.get("mats", {}))

    print("as-run index:", LAB_RUNS_CSV)
    print("orphan registry:", MAT_RUN_MAP_PATH)
    for mat_name in sorted(mat_names - known):
        print(f"  UNMAPPED mat on disk: {mat_name}")

    for mat_name, info in sorted(catalog.get("mats", {}).items()):
        run_name = info.get("run")
        if run_name is None:
            print(f"  {mat_name} -> (no OpenSees run)  [{info.get('note', '')}]")
            continue
        run_path = source_root / run_name
        status = "ok" if run_path.is_dir() else "MISSING DUMP"
        print(f"  {mat_name} -> {run_name}  [{status}]")

    for run_name in catalog.get("runs_without_mat", []):
        print(f"  run without mat: {run_name}")

    for mat_name in orphans.get("pending_mat_upload", []):
        print(f"  pending mat upload: {mat_name}")
    for mat_name in orphans.get("pending_dump_upload", []):
        print(f"  pending dump for mat: {mat_name}")


# ------------------------------------------------------------
# 6. CLI
# ------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    """
    Parse sync and extraction controls.

    Args:    none (reads sys.argv)
    Returns: argparse namespace
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--extract-only",
        action="store_true",
        help="only extract mats already on the local mirror",
    )
    parser.add_argument(
        "--no-extract",
        action="store_true",
        help="sync files only (ingest + archive)",
    )
    parser.add_argument(
        "--no-archive",
        action="store_true",
        help="skip LOCAL → Shared Drive archive copy",
    )
    parser.add_argument(
        "--force-extract",
        action="store_true",
        help="re-extract all mats",
    )
    return parser.parse_args()


def main() -> int:
    """
    Migrate layouts, ingest → LOCAL → archive, extract mats, report coverage.

    Args:    none (reads CLI flags)
    Returns: process status, 0 on success and 1 when ingest/local missing
    """
    args = parse_args()

    if migrate_legacy_local_layout():
        print("SyncLabBackup: migrated legacy opensees data/ → opensees_data/")
    if migrate_legacy_plots():
        print("SyncLabBackup: migrated legacy OSU_SSI_PLOTS/ → LOCAL/plots/")
    layout_stats = migrate_plots_layout()
    if any(layout_stats.values()):
        print(
            "SyncLabBackup: plots layout → runs/{eq,os}, mats/*/os, compare/ "
            f"({layout_stats})"
        )

    ingest = resolve_lab_ingest_source()
    archive = SHARED_DRIVE_ARCHIVE
    LOCAL_ROOT.mkdir(parents=True, exist_ok=True)
    LOCAL_OPENSEES_DATA.mkdir(parents=True, exist_ok=True)

    print(f"SyncLabBackup: ingest  {ingest or '(missing)'}")
    print(f"SyncLabBackup: local   {LOCAL_OPENSEES_DATA}")
    print(f"SyncLabBackup: archive {archive}  (<day>/opensees_data)")

    manifest = load_manifest()
    copied_mat_paths: list[Path] = []

    if not args.extract_only:
        if ingest is None:
            print(
                "SyncLabBackup: no ingest folder "
                "(shortcut-targets …/opensees data); skip ingest copy",
                file=sys.stderr,
            )
        else:
            copied_count, skipped_count, copied_mat_paths = sync_tree(
                ingest,
                LOCAL_OPENSEES_DATA,
                manifest,
                manifest_key="files",
                label="ingest",
            )
            print(
                f"SyncLabBackup: ingest copied={copied_count} "
                f"skipped={skipped_count}"
            )

        if not args.no_archive:
            if not _has_local_runs():
                print(
                    "SyncLabBackup: local mirror has no run dirs; skip archive",
                    file=sys.stderr,
                )
            else:
                ac, as_ = sync_archive_by_day(LOCAL_OPENSEES_DATA, manifest)
                print(
                    f"SyncLabBackup: archive copied={ac} skipped={as_}"
                )

        save_manifest(manifest)

        if MAT_RUN_MAP_PATH.is_file():
            LOCAL_MAT_RUN_MAP_PATH.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(MAT_RUN_MAP_PATH, LOCAL_MAT_RUN_MAP_PATH)

    if not args.no_extract:
        local_mat_paths = list_local_or_drive_mats(LOCAL_OPENSEES_DATA)
        if args.force_extract:
            extract_mats(local_mat_paths, force=True)
        elif copied_mat_paths:
            extract_mats(copied_mat_paths, force=False)
        else:
            extract_mats(local_mat_paths, force=False)

    report_map_coverage(LOCAL_OPENSEES_DATA)
    return 0


def _has_local_runs() -> bool:
    """True when LOCAL_OPENSEES_DATA has at least one r± run folder."""
    if not LOCAL_OPENSEES_DATA.is_dir():
        return False
    for child in LOCAL_OPENSEES_DATA.iterdir():
        if child.is_dir() and child.name.startswith(("r+", "r-")):
            return True
    return False


if __name__ == "__main__":
    raise SystemExit(main())
