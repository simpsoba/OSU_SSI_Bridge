#!/usr/bin/env python3
"""
Goals
-----
Mirror the read-only Shared Drive lab dumps into
`OSU_SSI_BRIDGE_DATA_LOCAL/opensees_data/`.

  python plot/SyncLabBackup.py
  python plot/SyncLabBackup.py --extract-only
  python plot/SyncLabBackup.py --no-extract
  python plot/SyncLabBackup.py --force-extract

Copy new or changed files into the LOCAL mirror. Skip every `plots/`
directory so generated figures never enter the backup tree. Extract the
Simulink `data` and *OS channels named by `lab_paths.MAT_KEYS` into
`LOCAL/opensees_data/mat_extract/<stem>.npz`.

The Shared Drive is read-only to this script. All writes go to LOCAL.
Legacy local mirror and plot layouts are migrated before sync or extract.
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
    LOCAL_MAT_RUN_MAP_PATH,
    LOCAL_OPENSEES_DATA,
    LOCAL_ROOT,
    MANIFEST_PATH,
    MAT_EXTRACT_DIR,
    MAT_KEYS,
    MAT_RUN_MAP_PATH,
    cell_str,
    load_mat_run_map,
    migrate_legacy_local_layout,
    migrate_legacy_plots,
    resolve_opensees_data,
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


def iter_source_files(source_root: Path) -> list[Path]:
    """
    Files below a source root, excluding anything inside `plots/`.

    Args:    source_root  Drive, junction, or LOCAL opensees_data folder
    Returns: file paths in pathlib traversal order
    """
    source_files: list[Path] = []
    for path in source_root.rglob("*"):
        if not path.is_file():
            continue
        relative_parts = {
            part.lower() for part in path.relative_to(source_root).parts
        }
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
    Returns: True when the file must be copied into LOCAL
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
# 2. DRIVE-TO-LOCAL SYNC
# ------------------------------------------------------------


def sync_tree(
    source_root: Path,
    destination_root: Path,
    manifest: dict,
) -> tuple[int, int, list[Path]]:
    """
    Copy new or changed files from the source tree into LOCAL.

    Args:    source_root, destination_root, manifest
    Returns: (copied count, skipped count, copied Simulink .mat paths)
    """
    destination_root.mkdir(parents=True, exist_ok=True)
    manifest_files = manifest.setdefault("files", {})
    copied_count = 0
    skipped_count = 0
    copied_mat_paths: list[Path] = []

    for source_path in iter_source_files(source_root):
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
        print(f"  copy {relative_path}")

    return copied_count, skipped_count, copied_mat_paths


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
    mat_run_map = load_mat_run_map()
    mat_names = {
        path.name for path in list_local_or_drive_mats(source_root)
    }
    known_mat_names = set(mat_run_map.get("mats", {}))

    print("mat_run_map:", MAT_RUN_MAP_PATH)
    for mat_name in sorted(mat_names - known_mat_names):
        print(f"  UNMAPPED mat: {mat_name}")

    for mat_name, info in sorted(mat_run_map.get("mats", {}).items()):
        run_name = info.get("run")
        if run_name is None:
            print(f"  {mat_name} -> (no OpenSees run)")
            continue
        run_path = source_root / run_name
        status = "ok" if run_path.is_dir() else "MISSING DUMP"
        print(f"  {mat_name} -> {run_name}  [{status}]")

    for run_name in mat_run_map.get("runs_without_mat", []):
        print(f"  run without mat: {run_name}")


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
        help="only extract mats already on the local mirror (or Drive if no mirror)",
    )
    parser.add_argument(
        "--no-extract",
        action="store_true",
        help="sync files only",
    )
    parser.add_argument(
        "--force-extract",
        action="store_true",
        help="re-extract all mats",
    )
    return parser.parse_args()


def main() -> int:
    """
    Migrate local layouts, sync Drive files, extract MAT files, and report.

    Args:    none (reads CLI flags)
    Returns: process status, 0 on success and 1 when no data root is found
    """
    args = parse_args()

    if migrate_legacy_local_layout():
        print("SyncLabBackup: migrated legacy opensees data/ → opensees_data/")
    if migrate_legacy_plots():
        print("SyncLabBackup: migrated legacy OSU_SSI_PLOTS/ → LOCAL/plots/")

    source_root = resolve_opensees_data()
    if source_root is None:
        print("SyncLabBackup: no opensees data folder found", file=sys.stderr)
        return 1

    print(f"SyncLabBackup: source {source_root}")
    print(f"SyncLabBackup: local  {LOCAL_OPENSEES_DATA}")
    LOCAL_ROOT.mkdir(parents=True, exist_ok=True)

    manifest = load_manifest()
    copied_mat_paths: list[Path] = []

    if not args.extract_only:
        if source_root.resolve() == LOCAL_OPENSEES_DATA.resolve():
            print("SyncLabBackup: source is already the local mirror; skip copy")
        else:
            copied_count, skipped_count, copied_mat_paths = sync_tree(
                source_root,
                LOCAL_OPENSEES_DATA,
                manifest,
            )
            save_manifest(manifest)
            print(
                f"SyncLabBackup: copied={copied_count} skipped={skipped_count}"
            )

        if MAT_RUN_MAP_PATH.is_file():
            LOCAL_MAT_RUN_MAP_PATH.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(MAT_RUN_MAP_PATH, LOCAL_MAT_RUN_MAP_PATH)

    if not args.no_extract:
        local_mat_paths = list_local_or_drive_mats(LOCAL_OPENSEES_DATA)
        source_mat_paths = list_local_or_drive_mats(source_root)
        all_mat_paths = local_mat_paths if local_mat_paths else source_mat_paths

        if args.force_extract:
            extract_mats(all_mat_paths, force=True)
        elif copied_mat_paths:
            extract_mats(copied_mat_paths, force=False)
        else:
            extract_mats(all_mat_paths, force=False)

    coverage_root = (
        source_root if source_root.is_dir() else LOCAL_OPENSEES_DATA
    )
    report_map_coverage(coverage_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
