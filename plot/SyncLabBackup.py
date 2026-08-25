#!/usr/bin/env python3
"""Mirror Drive lab dumps locally and extract Simulink *OS / data channels.

  python plot/SyncLabBackup.py
  python plot/SyncLabBackup.py --extract-only
  python plot/SyncLabBackup.py --no-extract

Copies new/changed files under the Drive `opensees data` tree into
`OSU_SSI_BRIDGE_DATA_LOCAL/opensees_data/` (skips any `plots/` dirs).
When new or larger `.mat` files appear, extracts `data` + *OS blocks to
`OSU_SSI_BRIDGE_DATA_LOCAL/opensees_data/mat_extract/<stem>.npz`.

Drive stays read-only; this never writes back to Drive.
On first run after an upgrade, moves legacy `opensees data/` into `opensees_data/`.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

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


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def file_fingerprint(path: Path) -> dict:
    st = path.stat()
    return {"size": int(st.st_size), "mtime_ns": int(st.st_mtime_ns)}


def load_manifest() -> dict:
    if not MANIFEST_PATH.is_file():
        return {"files": {}, "updated": None}
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


def save_manifest(man: dict) -> None:
    MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    man["updated"] = utc_now()
    MANIFEST_PATH.write_text(json.dumps(man, indent=2) + "\n", encoding="utf-8")


def iter_source_files(src_root: Path) -> list[Path]:
    """All files under src_root except anything inside a plots/ directory."""
    out: list[Path] = []
    for p in src_root.rglob("*"):
        if not p.is_file():
            continue
        parts = {x.lower() for x in p.relative_to(src_root).parts}
        if "plots" in parts:
            continue
        out.append(p)
    return out


def needs_copy(src: Path, dst: Path, old: dict | None) -> bool:
    if not dst.is_file():
        return True
    fp = file_fingerprint(src)
    if old is None:
        return fp != file_fingerprint(dst)
    return fp["size"] != old.get("size") or fp["mtime_ns"] != old.get("mtime_ns")


def sync_tree(src_root: Path, dst_root: Path, man: dict) -> tuple[int, int, list[Path]]:
    """Copy new/changed files. Returns (copied, skipped, new_or_updated_mats)."""
    dst_root.mkdir(parents=True, exist_ok=True)
    files = man.setdefault("files", {})
    copied = 0
    skipped = 0
    mats: list[Path] = []
    for src in iter_source_files(src_root):
        rel = src.relative_to(src_root).as_posix()
        dst = dst_root / rel
        old = files.get(rel)
        if not needs_copy(src, dst, old):
            skipped += 1
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        files[rel] = file_fingerprint(src)
        copied += 1
        if src.suffix.lower() == ".mat" and "Simulink" in src.parts:
            mats.append(dst)
        print(f"  copy {rel}")
    return copied, skipped, mats


def unwrap_struct(arr):
    if isinstance(arr, np.ndarray) and arr.dtype == object and arr.size == 1:
        return arr.flat[0]
    return arr


def extract_os_block(block) -> dict:
    names = signal_names(block)
    data = np.asarray(block.data, dtype=float)
    return {"signalNames": names, "data": data}


def extract_data_container(block) -> dict:
    """Keep `data` metadata + nested *OS arrays (not topL / AFC)."""
    out: dict = {
        "fileName": [cell_str(x) for x in block.fileName.flatten()],
        "sigNames": [cell_str(x) for x in block.sigNames.flatten()],
    }
    for key in ("tarSigOS", "comSigOS", "meaSigOS", "stateOS"):
        if not hasattr(block, key):
            continue
        nested = unwrap_struct(getattr(block, key))
        if hasattr(nested, "signalNames") and hasattr(nested, "data"):
            out[key] = extract_os_block(nested)
    # values cell: store shapes only if not already covered by nested *OS
    if hasattr(block, "values"):
        shapes = []
        for v in block.values.flatten():
            x = unwrap_struct(v)
            if hasattr(x, "data"):
                shapes.append(list(np.asarray(x.data).shape))
            elif isinstance(x, np.ndarray):
                shapes.append(list(x.shape))
            else:
                shapes.append(None)
        out["values_shapes"] = shapes
    return out


def extract_mat(mat_path: Path, out_npz: Path) -> dict:
    print(f"  extract {mat_path.name}")
    raw = sio.loadmat(
        str(mat_path),
        variable_names=list(MAT_KEYS),
        squeeze_me=False,
        struct_as_record=False,
    )
    payload: dict = {"mat": mat_path.name, "extracted": utc_now()}
    meta: dict = {}
    for key in MAT_KEYS:
        if key not in raw:
            continue
        block = unwrap_struct(raw[key])
        if key == "data":
            dc = extract_data_container(block)
            # flatten nested OS into npz keys under data_
            meta["data_fileName"] = np.array(dc["fileName"], dtype=object)
            meta["data_sigNames"] = np.array(dc["sigNames"], dtype=object)
            if "values_shapes" in dc:
                meta["data_values_shapes"] = np.array(dc["values_shapes"], dtype=object)
            for os_key in ("tarSigOS", "comSigOS", "meaSigOS", "stateOS"):
                if os_key not in dc:
                    continue
                payload[f"data_{os_key}_data"] = dc[os_key]["data"]
                meta[f"data_{os_key}_signalNames"] = np.array(
                    dc[os_key]["signalNames"], dtype=object
                )
            continue
        if hasattr(block, "signalNames") and hasattr(block, "data"):
            osb = extract_os_block(block)
            payload[f"{key}_data"] = osb["data"]
            meta[f"{key}_signalNames"] = np.array(osb["signalNames"], dtype=object)
            # convenience: Time + primary signal for tar/mea
            names = osb["signalNames"]
            try:
                ti = time_column(names)
                payload[f"{key}_time"] = osb["data"][:, ti]
                # first non-Time column as primary
                for j, n in enumerate(names):
                    if j != ti:
                        payload[f"{key}_primary"] = osb["data"][:, j]
                        meta[f"{key}_primaryName"] = np.array([n], dtype=object)
                        break
            except KeyError:
                pass
    out_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out_npz, **payload, **meta)
    summary = {
        "mat": mat_path.name,
        "npz": out_npz.name,
        "keys": sorted(k for k in payload if k.endswith("_data") or k == "mat"),
        "n_tar": int(payload["tarSigOS_data"].shape[0])
        if "tarSigOS_data" in payload
        else 0,
        "t_real_last": float(payload["tarSigOS_time"][-1])
        if "tarSigOS_time" in payload
        else None,
    }
    return summary


def extract_mats(mat_paths: list[Path], force: bool = False) -> list[dict]:
    MAT_EXTRACT_DIR.mkdir(parents=True, exist_ok=True)
    summaries: list[dict] = []
    for mat_path in mat_paths:
        out_npz = MAT_EXTRACT_DIR / f"{mat_path.stem}.npz"
        if out_npz.is_file() and not force:
            src_m = mat_path.stat().st_mtime_ns
            dst_m = out_npz.stat().st_mtime_ns
            if dst_m >= src_m:
                print(f"  skip extract (up to date) {mat_path.name}")
                continue
        summaries.append(extract_mat(mat_path, out_npz))
    return summaries


def list_local_or_drive_mats(src_root: Path) -> list[Path]:
    sim = resolve_simulink_dir(src_root)
    if sim is None:
        return []
    return sorted(sim.glob("*.mat"))


def report_map_coverage(src_root: Path) -> None:
    mmap = load_mat_run_map()
    mats = {p.name for p in list_local_or_drive_mats(src_root)}
    known = set(mmap.get("mats", {}))
    print("mat_run_map:", MAT_RUN_MAP_PATH)
    for name in sorted(mats - known):
        print(f"  UNMAPPED mat: {name}")
    for name, info in sorted(mmap.get("mats", {}).items()):
        run = info.get("run")
        if run is None:
            print(f"  {name} -> (no OpenSees run)")
            continue
        run_path = src_root / run
        ok = run_path.is_dir()
        print(f"  {name} -> {run}  [{'ok' if ok else 'MISSING DUMP'}]")
    for run in mmap.get("runs_without_mat", []):
        print(f"  run without mat: {run}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--extract-only",
        action="store_true",
        help="only extract mats already on the local mirror (or Drive if no mirror)",
    )
    ap.add_argument("--no-extract", action="store_true", help="sync files only")
    ap.add_argument("--force-extract", action="store_true", help="re-extract all mats")
    args = ap.parse_args()

    if migrate_legacy_local_layout():
        print("SyncLabBackup: migrated legacy opensees data/ → opensees_data/")
    if migrate_legacy_plots():
        print("SyncLabBackup: migrated legacy OSU_SSI_PLOTS/ → LOCAL/plots/")

    src = resolve_opensees_data()
    if src is None:
        print("SyncLabBackup: no opensees data folder found", file=sys.stderr)
        return 1
    # Prefer Drive as source when both exist
    if LOCAL_OPENSEES_DATA in src.parents or src == LOCAL_OPENSEES_DATA:
        # only local available
        pass

    print(f"SyncLabBackup: source {src}")
    print(f"SyncLabBackup: local  {LOCAL_OPENSEES_DATA}")
    LOCAL_ROOT.mkdir(parents=True, exist_ok=True)

    man = load_manifest()
    new_mats: list[Path] = []

    if not args.extract_only:
        if src.resolve() == LOCAL_OPENSEES_DATA.resolve():
            print("SyncLabBackup: source is already the local mirror; skip copy")
        else:
            copied, skipped, new_mats = sync_tree(src, LOCAL_OPENSEES_DATA, man)
            save_manifest(man)
            print(f"SyncLabBackup: copied={copied} skipped={skipped}")
        if MAT_RUN_MAP_PATH.is_file():
            LOCAL_MAT_RUN_MAP_PATH.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(MAT_RUN_MAP_PATH, LOCAL_MAT_RUN_MAP_PATH)

    if not args.no_extract:
        # Prefer local Simulink mirror; fall back to Drive if empty / missing
        local_mats = list_local_or_drive_mats(LOCAL_OPENSEES_DATA)
        drive_mats = list_local_or_drive_mats(src)
        mats = local_mats if local_mats else drive_mats
        if args.force_extract:
            extract_mats(mats, force=True)
        elif new_mats:
            extract_mats(new_mats, force=False)
        else:
            # first run / extract-only: extract any missing npz
            extract_mats(mats, force=False)

    report_map_coverage(src if src.is_dir() else LOCAL_OPENSEES_DATA)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
