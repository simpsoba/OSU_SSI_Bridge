#!/usr/bin/env python3
"""EQ postprocess for OpenSeesMP dumps (name.$pid).

  python3 plot/PlotEQParallel.py
  python3 plot/PlotEQParallel.py /path/to/eqOutDir

Reads np from window_meta.txt.0 (must match ranks 0..np-1). Stitches a serial
folder, then calls PlotEQ.py. Lean window quads get corners from
model_sketch.json here (PlotEQ.py stays serial-only).
Panels: the DO_* switches in PlotEQ.py. Plots go to <eqOutDir>/plots/.
"""

from __future__ import annotations

import json
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

import PlotEQ as peq
from paths import HERE, elevation_dir, eq_dir

EQ_OUT = eq_dir(4, "Shin", "SSPquad", "lumpedPlasticity", "parallel")

# Recorders that live on one rank (copy the first non-empty shard).
# pier_node_$tag.out is one file per pier node, so each one is also single-rank.
ONE_RANK = (
    "pier_hinge_force.out",
    "pier_hinge_defo.out",
    "pier_hinge_top_force.out",
    "pier_hinge_top_defo.out",
    "pier_top_disp.out",
    "pier_top_acc.out",
    "soil_base_primary.out",
)

# Time-series whose columns follow concat of the matching *_eles.txt ranks.
HSTACK = (
    "pile_beam_globalForce.out",
    "pile_beam_sec1_defo.out",
    "pile_springs_force.out",
    "pile_springs_defo.out",
    "cap_springs_force.out",
    "cap_springs_defo.out",
    "cap_springs_soffit_force.out",
    "cap_springs_soffit_defo.out",
)

TEXT_UNIQUE_ELE = (
    "window_eles.txt",
    "window_quads.txt",
    "pile_beam_eles.txt",
    "pile_springs_eles.txt",
    "cap_springs_eles.txt",
)

# window_nodes.txt and disp_nodes.txt are stitched with the disp columns instead.


def rank_files(eq: Path, name: str) -> dict[int, Path]:
    out: dict[int, Path] = {}
    prefix = name + "."
    for p in eq.iterdir():
        if not p.is_file() or not p.name.startswith(prefix):
            continue
        suf = p.name[len(prefix) :]
        if suf.isdigit():
            out[int(suf)] = p
    return dict(sorted(out.items()))


def parse_meta(path: Path) -> dict[str, str]:
    meta: dict[str, str] = {}
    for ln in peq._skip_hash(path):
        k, _, rest = ln.partition(" ")
        meta[k] = rest.strip()
    return meta


def load_np(eq: Path) -> tuple[int, dict[int, dict[str, str]]]:
    files = rank_files(eq, "window_meta.txt")
    if 0 not in files:
        raise SystemExit(
            f"PlotEQParallel: missing window_meta.txt.0 in {eq}"
        )
    metas = {pid: parse_meta(p) for pid, p in files.items()}
    np_run = int(float(metas[0].get("np", 0) or 0))
    if np_run < 2:
        raise SystemExit(
            f"PlotEQParallel: meta np={np_run}; serial dump -> PlotEQ.py"
        )
    expect = set(range(np_run))
    got = set(metas)
    extra = sorted(got - expect)
    missing = sorted(expect - got)
    if extra or missing:
        raise SystemExit(
            f"PlotEQParallel: meta np={np_run} but ranks {sorted(got)}"
            + (f"  extra={extra}" if extra else "")
            + (f"  missing={missing}" if missing else "")
            + "\n  leftover name.$pid from another -np? EQRecorders rank 0 "
            "wipes files in eqOutDir before a new run."
        )
    return np_run, metas


def one_rank_names(eq: Path) -> list[str]:
    """ONE_RANK plus the pier_node_$tag.out files this run happened to write."""
    names = list(ONE_RANK)
    for p in sorted(eq.glob("pier_node_*.out.*")):
        stem = p.name.rsplit(".", 1)[0]
        if stem not in names:
            names.append(stem)
    return names


def shard_or_none(eq: Path, name: str, pid: int) -> Path | None:
    p = eq / f"{name}.{pid}"
    return p if p.is_file() else None


def concat_text(
    eq: Path,
    dest: Path,
    name: str,
    np_run: int,
    unique_col0: bool,
) -> int:
    header = None
    rows: list[str] = []
    seen: set[str] = set()
    for pid in range(np_run):
        p = shard_or_none(eq, name, pid)
        if p is None:
            continue
        for ln in p.read_text().splitlines():
            s = ln.strip()
            if not s:
                continue
            if s.startswith("#"):
                if header is None:
                    header = ln.rstrip()
                continue
            key = s.split()[0] if unique_col0 else s
            if unique_col0:
                if key in seen:
                    continue
                seen.add(key)
            rows.append(s)
    out = dest / name
    with out.open("w") as f:
        if header:
            f.write(header + "\n")
        for s in rows:
            f.write(s + "\n")
    return len(rows)


def load_rank_disp(
    eq: Path, tags: list[int], disp_files: list[str]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return peq.load_window_disp(eq, tags, disp_files)


def parse_node_shard(path: Path) -> tuple[str | None, list[int], list[tuple[str, str]]]:
    """One window_nodes.txt / disp_nodes.txt shard: header, tags, (x, y) strings."""
    header = None
    tags: list[int] = []
    xy: list[tuple[str, str]] = []
    for ln in path.read_text().splitlines():
        s = ln.strip()
        if not s:
            continue
        if s.startswith("#"):
            if header is None:
                header = ln.rstrip()
            continue
        a, b, c = s.split()[:3]
        tags.append(int(a))
        xy.append((b, c))
    return header, tags, xy


def write_node_file(
    path: Path, header: str | None, tags: list[int], xy: dict[int, tuple[str, str]]
) -> None:
    with path.open("w") as f:
        if header:
            f.write(header + "\n")
        for tg in tags:
            b, c = xy[tg]
            f.write(f"{tg} {b} {c}\n")


def stitch_nodes_disp(
    eq: Path, dest: Path, np_run: int, metas: dict[int, dict[str, str]]
) -> tuple[int, int, int]:
    """Unique nodes (drop ghosts): geometry -> window_nodes.txt, columns ->
    disp_nodes.txt + window_disp.out. recordersON=2 records displacement for a
    subset of the geometry, so the two lists differ."""
    geom_header = None
    disp_header = None
    geom_tags: list[int] = []
    disp_tags: list[int] = []
    geom_xy: dict[int, tuple[str, str]] = {}
    disp_xy: dict[int, tuple[str, str]] = {}
    n_raw = 0
    ux_keep: list[np.ndarray] = []
    uy_keep: list[np.ndarray] = []
    t = None
    for pid in range(np_run):
        gpth = shard_or_none(eq, "window_nodes.txt", pid)
        if gpth is None:
            raise SystemExit(f"PlotEQParallel: missing window_nodes.txt.{pid}")
        hdr, gtags, gxy = parse_node_shard(gpth)
        if geom_header is None:
            geom_header = hdr
        n_raw += len(gtags)
        for i, tg in enumerate(gtags):
            if tg in geom_xy:
                continue
            geom_xy[tg] = gxy[i]
            geom_tags.append(tg)
        dpth = shard_or_none(eq, "disp_nodes.txt", pid) or gpth
        dhdr, dtags, dxy = parse_node_shard(dpth)
        if disp_header is None:
            disp_header = dhdr
        dfs = metas[pid].get("dispFiles", "").split()
        dfs = [f for f in dfs if (eq / f).is_file()]
        keep_idx: list[int] = []
        for i, tg in enumerate(dtags):
            if tg in disp_xy:
                continue
            disp_xy[tg] = dxy[i]
            disp_tags.append(tg)
            keep_idx.append(i)
        if not dtags or not dfs or not keep_idx:
            continue
        ti, uxr, uyr = load_rank_disp(eq, dtags, dfs)
        if t is None:
            t = ti
        n = min(len(t), len(ti), uxr.shape[0])
        t = t[:n]
        sel = np.asarray(keep_idx, dtype=int)
        ux_keep.append(uxr[:n][:, sel])
        uy_keep.append(uyr[:n][:, sel])
    write_node_file(dest / "window_nodes.txt", geom_header, geom_tags, geom_xy)
    write_node_file(dest / "disp_nodes.txt", disp_header, disp_tags, disp_xy)
    if not ux_keep:
        raise SystemExit("PlotEQParallel: no window disp on any rank")
    n = min(p.shape[0] for p in ux_keep)
    ux = np.hstack([p[:n] for p in ux_keep])
    uy = np.hstack([p[:n] for p in uy_keep])
    t = t[:n]
    n_disp = len(disp_tags)
    if ux.shape[1] != n_disp:
        raise SystemExit(
            f"PlotEQParallel: disp cols {ux.shape[1]} != unique disp nodes {n_disp}"
        )
    mixed = np.empty((n, 1 + 2 * n_disp))
    mixed[:, 0] = t
    mixed[:, 1::2] = ux
    mixed[:, 2::2] = uy
    np.savetxt(dest / "window_disp.out", mixed, fmt="%.10g")
    return n_raw, len(geom_tags), n_disp


def hstack_recorders(eq: Path, dest: Path, name: str, np_run: int) -> None:
    parts: list[np.ndarray] = []
    t = None
    for pid in range(np_run):
        p = shard_or_none(eq, name, pid)
        if p is None or p.stat().st_size < 100:
            continue
        a = peq.loadtxt_partial(p)
        if a.size == 0:
            continue
        if t is None:
            t = a[:, 0]
        n = min(len(t), a.shape[0])
        t = t[:n]
        parts.append(a[:n, 1:])
    if not parts:
        return
    n = min(p.shape[0] for p in parts)
    data = np.hstack([p[:n] for p in parts])
    out = np.column_stack([t[:n], data])
    np.savetxt(dest / name, out, fmt="%.10g")


def copy_one_rank(eq: Path, dest: Path, name: str, np_run: int) -> None:
    for pid in range(np_run):
        p = shard_or_none(eq, name, pid)
        if p is None or p.stat().st_size < 100:
            continue
        shutil.copy2(p, dest / name)
        return


def stitch_quads(
    eq: Path, dest: Path, np_run: int, metas: dict[int, dict[str, str]]
) -> tuple[list[str], list[str], int]:
    nq = concat_text(eq, dest, "window_quads.txt", np_run, unique_col0=True)
    sig: list[str] = []
    eps: list[str] = []
    for pid in range(np_run):
        sig.extend(
            f for f in metas[pid].get("quadStressFiles", "").split() if (eq / f).is_file()
        )
        eps.extend(
            f for f in metas[pid].get("quadStrainFiles", "").split() if (eq / f).is_file()
        )
    if sig:
        hstack_named_files(eq, dest, "window_quad_stress.out", sig)
    if eps:
        hstack_named_files(eq, dest, "window_quad_strain.out", eps)
    return (
        ["window_quad_stress.out"] if sig else [],
        ["window_quad_strain.out"] if eps else [],
        nq,
    )


def hstack_named_files(
    eq: Path, dest: Path, out_name: str, files: list[str]
) -> None:
    parts: list[np.ndarray] = []
    t = None
    for fn in files:
        a = peq.loadtxt_partial(eq / fn)
        if a.size == 0:
            continue
        if t is None:
            t = a[:, 0]
        n = min(len(t), a.shape[0])
        t = t[:n]
        parts.append(a[:n, 1:])
    if not parts:
        return
    n = min(p.shape[0] for p in parts)
    data = np.hstack([p[:n] for p in parts])
    np.savetxt(dest / out_name, np.column_stack([t[:n], data]), fmt="%.10g")


def write_meta(
    dest: Path,
    base: dict[str, str],
    n_nodes: int,
    n_eles: int,
    n_quads: int,
    sig: list[str],
    eps: list[str],
    np_run: int,
    pier_files: list[str] | None = None,
    n_disp: int = 0,
) -> None:
    skip = {
        "dispFiles",
        "pierNodeFiles",
        "nDispNodes",
        "pid",
        "np",
        "fileSuffix",
        "nWindowNodes",
        "nWindowEles",
        "nWindowQuads",
        "quadStressFiles",
        "quadStrainFiles",
        "quadStrainRsp",
        "quadStressRsp",
        "quadEleFile",
    }
    lines = []
    for k, v in base.items():
        if k in skip:
            continue
        lines.append(f"{k} {v}")
    lines.append(f"nWindowNodes {n_nodes}")
    lines.append(f"nDispNodes {n_disp}")
    lines.append(f"nWindowEles {n_eles}")
    lines.append(f"nWindowQuads {n_quads}")
    lines.append("dispNodesFile disp_nodes.txt")
    lines.append("dispFiles window_disp.out")
    if pier_files:
        lines.append(f"pierNodeFiles {' '.join(pier_files)}")
    lines.append(f"np {np_run}")
    if sig:
        lines.append("quadEleFile window_quads.txt")
        if "quadNgp" not in base:
            lines.append("quadNgp 1")
        lines.append(f"quadStressFiles {' '.join(sig)}")
        lines.append(f"quadStressRsp {base.get('quadStressRsp', 'stress2D3')}")
        lines.append("quadStrainRsp " + base.get("quadStrainRsp", "strain2D3"))
        lines.append(f"quadStrainFiles {' '.join(eps)}")
    (dest / "window_meta.txt").write_text("\n".join(lines) + "\n")


def sketch_quad_xy(meta: dict) -> dict[int, list[tuple[float, float]]]:
    """Undeformed 4-node corners from DumpModelSketch. Keys = ele tags."""
    sp = meta.get("soilProfile", "")
    bnd = meta.get("soilBoundary", "Shin")
    cands = []
    if sp:
        cands.append(elevation_dir(sp, bnd) / "model_sketch.json")
    cands.append(HERE / "model_sketch.json")
    js = {}
    for pth in cands:
        if pth.is_file():
            with pth.open() as f:
                js = json.load(f)
            break
    out: dict[int, list[tuple[float, float]]] = {}
    for row in js.get("soil_quads") or []:
        raw = row.get("xy") or []
        if len(raw) < 8:
            continue
        pts = [(float(raw[2 * k]), float(raw[2 * k + 1])) for k in range(4)]
        out[int(row["e"])] = pts
    return out


def fill_lean_quad_geom(dest: Path, meta: dict) -> int:
    """Add sketch corners so PlotEQ can find 4 nodes per window quad.

    Only needed for dumps whose window_nodes.txt misses a quad corner. Geometry
    only: these tags get no displacement column.
    """
    sketch = sketch_quad_xy(meta)
    qtags = peq.read_window_quad_list(dest)
    ev = peq.read_ele_nodes(dest)
    tags, xy = peq.read_nodes(dest)
    extra_nodes: list[tuple[int, float, float]] = []
    extra_eles: list[str] = []
    nid = 910000
    n_fill = 0
    have = set(xy)
    for t in qtags:
        nn = ev.get(t, [])
        if len(nn) >= 4 and all(n in have for n in nn[:4]):
            continue
        pts = sketch.get(t)
        if not pts:
            continue
        ids = []
        for px, py in pts:
            nid += 1
            extra_nodes.append((nid, px, py))
            ids.append(nid)
            have.add(nid)
        extra_eles.append(f"{t} {ids[0]} {ids[1]} {ids[2]} {ids[3]}")
        n_fill += 1
    if not extra_nodes:
        return 0
    with (dest / "window_nodes.txt").open("a") as f:
        for tg, px, py in extra_nodes:
            f.write(f"{tg} {px} {py}\n")
    with (dest / "window_eles.txt").open("a") as f:
        for ln in extra_eles:
            f.write(ln + "\n")
    return n_fill


def stitch(eq: Path, dest: Path, np_run: int, metas: dict[int, dict[str, str]]) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    n_raw, n_u, n_disp = stitch_nodes_disp(eq, dest, np_run, metas)
    print(
        f"PlotEQParallel: unique nodes {n_u} (dropped {n_raw - n_u} ghosts), "
        f"disp columns for {n_disp}"
    )
    n_eles = concat_text(eq, dest, "window_eles.txt", np_run, unique_col0=True)
    for name in TEXT_UNIQUE_ELE:
        if name == "window_eles.txt" or name == "window_quads.txt":
            continue
        concat_text(eq, dest, name, np_run, unique_col0=True)
    sig, eps, n_quads = stitch_quads(eq, dest, np_run, metas)
    for name in HSTACK:
        hstack_recorders(eq, dest, name, np_run)
    pier_files = []
    for name in one_rank_names(eq):
        copy_one_rank(eq, dest, name, np_run)
        if name.startswith("pier_node_") and (dest / name).is_file():
            pier_files.append(name)
    write_meta(
        dest, metas[0], n_u, n_eles, n_quads, sig, eps, np_run, pier_files, n_disp
    )
    n_fill = fill_lean_quad_geom(dest, metas[0])
    if n_fill:
        print(f"PlotEQParallel: sketch corners for {n_fill} window quads")


HELP = """\
usage: python3 plot/PlotEQParallel.py [eqOutDir]

  eqOutDir   MP recorder folder (default: EQ_OUT in this file)
  np         from window_meta.txt.0; ranks must be 0..np-1
  plots      <eqOutDir>/plots/   (PlotEQ.py DO_* switches)
"""


def main() -> int:
    if any(a in ("-h", "--help") for a in sys.argv[1:]):
        print(HELP, end="")
        return 0
    eq = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(EQ_OUT)
    if (eq / "window_nodes.txt").is_file() and not (eq / "window_nodes.txt.0").is_file():
        print(
            f"PlotEQParallel: serial dump in {eq}\n  python3 plot/PlotEQ.py {eq}",
            file=sys.stderr,
        )
        return 1
    np_run, metas = load_np(eq)
    print(f"PlotEQParallel: np={np_run}  {eq}")
    tmp = Path(tempfile.mkdtemp(prefix="eqmp_"))
    try:
        stitch(eq, tmp, np_run, metas)
        argv = sys.argv
        sys.argv = [argv[0], str(tmp)]
        try:
            rc = peq.main()
        finally:
            sys.argv = argv
        src = tmp / "plots"
        dst = eq / "plots"
        if src.is_dir():
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(src, dst)
            print(f"PlotEQParallel: plots -> {dst}")
        return rc
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
