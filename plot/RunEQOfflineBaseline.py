#!/usr/bin/env python3
"""
Goals
-----
Run the four baseline structural cases (fixed/free × lumped/FBC) with
``MKRAlphaExplicitMultiSOE 0.5``, writing recorder dumps only under
``plot/out/eq_offline/``.

Does not touch ``OSU_SSI_BRIDGE_DATA*``, Shared Drive, ``plot/out/eigen/``,
``modal_analysis/``, or root ``Overrides.tcl``.

Usage
-----
  python plot/RunEQOfflineBaseline.py --smoke
  python plot/RunEQOfflineBaseline.py
  python plot/RunEQOfflineBaseline.py --only fbc_free_mkr lumped_hold_mkr
  python plot/RunEQOfflineBaseline.py --skip-plot
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OV_DIR = REPO / "analysis" / "eq_offline" / "overrides"
OUT_ROOT = REPO / "plot" / "out" / "eq_offline"
LOG_DIR = OUT_ROOT / "logs"
PY = sys.executable

# Builds that ship MKRAlphaExplicitMultiSOE (stock OpenSees PATH often does not).
MKR_OPENSEES_CANDIDATES = (
    Path(r"C:\Users\garaujor\source\repos\simpsoba\OpenSees\build-cuda\Release\OpenSees.exe"),
    Path(r"C:\Users\garaujor\source\repos\OpenSees-CUDA\build\Release\OpenSees.exe"),
    Path(r"C:\Users\garaujor\source\repos\simpsoba\RTHS-CUDA\OpenSees\build\Release\OpenSees.exe"),
    Path(r"C:\Users\garaujor\OpenSees_Runs\RTHS-CUDA-ELAINE\RTHS-CUDA\OpenSees\build\Release\OpenSees.exe"),
)

# Do not write here (lab ingest / archive / other agent).
FORBIDDEN_PREFIXES = (
    "OSU_SSI_BRIDGE_DATA",
    "OSU_SSI_BRIDGE_DATA_LOCAL",
    "plot/out/eigen",
    "modal_analysis",
)

CASES: list[dict] = []
for pier, hold, pier_slug, hold_slug in (
    ("lumpedPlasticity", 1, "lumped", "hold"),
    ("lumpedPlasticity", 0, "lumped", "free"),
    ("forceBeamColumn", 1, "fbc", "hold"),
    ("forceBeamColumn", 0, "fbc", "free"),
):
    slug = f"{pier_slug}_{hold_slug}_mkr"
    CASES.append(
        {
            "slug": slug,
            "pierEleType": pier,
            "holdPierON": hold,
            "eqIntegrator": "MKRAlphaExplicitMultiSOE 0.5",
        }
    )
# Newmark companion for hold+lumped (implicit; may need smaller DT_FACTOR).
CASES.append(
    {
        "slug": "lumped_hold_newmark",
        "pierEleType": "lumpedPlasticity",
        "holdPierON": 1,
        "eqIntegrator": "Newmark 0.5 0.25",
    }
)


def _safe_outdir(rel: str) -> Path:
    """Resolve a repo-relative outDIR and refuse lab / eigen trees."""
    rel_norm = rel.replace("\\", "/").lstrip("./")
    for bad in FORBIDDEN_PREFIXES:
        if rel_norm == bad or rel_norm.startswith(bad + "/"):
            raise SystemExit(f"refusing outDIR under forbidden tree: {rel}")
    out = (REPO / rel_norm).resolve()
    try:
        out.relative_to(OUT_ROOT.resolve())
    except ValueError as exc:
        raise SystemExit(f"outDIR must stay under {OUT_ROOT}: {out}") from exc
    return out


def write_override(case: dict) -> Path:
    """Write one argv Overrides.tcl for an offline baseline EQ case."""
    OV_DIR.mkdir(parents=True, exist_ok=True)
    slug = case["slug"]
    out_rel = f"plot/out/eq_offline/{slug}"
    _safe_outdir(out_rel)
    path = OV_DIR / f"{slug}.tcl"
    # Newmark average-accel: smaller dt + slightly looser NormDispIncr (hold+lumped
    # previously died near t~54 s at DT_FACTOR=10/5 with 1e-8).
    dt_line = ""
    if "Newmark" in case["eqIntegrator"]:
        dt_line = (
            "set DT_FACTOR 5\n"
            "set eqTestTol 1.0e-6\n"
            "set eqTestIter 50\n"
            "set eqRecoverTol 1.0e-5\n"
            "set eqRecoverIter 50\n"
        )
    text = f"""# Auto-written by plot/RunEQOfflineBaseline.py — offline baseline EQ only.
# Dump: {out_rel}/  (never lab data / eigen / modal_analysis)
set runEQ 1
set plotFigures 0
set realTimeON 0
set recordersON 1
set eqPrintON 1
set eqPrintDt 10.0
set holdPierON {case["holdPierON"]}
set holdPierRZON 0
set pierEleType "{case["pierEleType"]}"
set soilMesh 0
set soilProfile 4
set soilBoundary "Shin"
set soilEleType "SSPquad"
set soilConstitutive "inelastic"
set eqIntegrator "{case["eqIntegrator"]}"
set prePartitionSystem "UmfPack"
set postPartitionSystem "UmfPack"
set constraintsHandler "Transformation"
{dt_line}set gmDir [file join $gmRoot Tohoku2011-FKSH]
set gmVelFile [file join $gmDir FKSH19.NS1.VT2]
set gmStartTime 0.0
set outDIR "{out_rel}"
"""
    path.write_text(text, encoding="utf-8", newline="\n")
    return path


def find_opensees() -> str:
    """Prefer a CUDA/MKR build; PATH stock OpenSees often lacks MKRAlpha."""
    env = os.environ.get("OPENSEES_EXE", "").strip()
    if env:
        p = Path(env)
        if not p.is_file():
            raise SystemExit(f"OPENSEES_EXE not found: {env}")
        return str(p)
    for cand in MKR_OPENSEES_CANDIDATES:
        if cand.is_file():
            return str(cand)
    wh = shutil.which("OpenSees") or shutil.which("OpenSees.exe")
    if wh:
        return wh
    raise SystemExit(
        "OpenSees with MKRAlphaExplicitMultiSOE not found. "
        "Set OPENSEES_EXE to a CUDA build."
    )


def run_case(case: dict, *, smoke: bool, skip_plot: bool) -> int:
    """Run one OpenSees case and optionally PlotEQ."""
    ov = write_override(case)
    slug = case["slug"]
    out = OUT_ROOT / slug
    out.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = LOG_DIR / f"{slug}.log"
    exe = find_opensees()
    env = os.environ.copy()
    if smoke:
        env["REGEN_EQ_TMAX"] = env.get("REGEN_EQ_TMAX", "20")
        env["REGEN_FREE_VIB"] = env.get("REGEN_FREE_VIB", "0")
    else:
        env.pop("REGEN_EQ_TMAX", None)
        env.pop("REGEN_FREE_VIB", None)

    cmd = [exe, str(REPO / "Run.tcl"), str(ov)]
    print(f"=== {slug} ===", flush=True)
    print(f"  OpenSees: {' '.join(cmd)}", flush=True)
    print(f"  log: {log}", flush=True)
    t0 = time.time()
    with log.open("w", encoding="utf-8", errors="replace") as fd:
        fd.write(f"# cmd: {cmd}\n# cwd: {REPO}\n")
        if smoke:
            fd.write(
                f"# smoke REGEN_EQ_TMAX={env.get('REGEN_EQ_TMAX')} "
                f"REGEN_FREE_VIB={env.get('REGEN_FREE_VIB')}\n"
            )
        fd.flush()
        proc = subprocess.run(
            cmd,
            cwd=str(REPO),
            env=env,
            stdout=fd,
            stderr=subprocess.STDOUT,
            text=True,
        )
    elapsed = time.time() - t0
    print(f"  exit={proc.returncode}  wall={elapsed / 60.0:.1f} min", flush=True)
    log_txt = log.read_text(encoding="utf-8", errors="replace")
    if (
        proc.returncode != 0
        or "analyze failed" in log_txt
        or "No Integrator type exists" in log_txt
    ):
        print(f"  FAILED — see {log}", flush=True)
        return proc.returncode if proc.returncode != 0 else 3

    meta = out / "window_meta.txt"
    if not meta.is_file():
        print(f"  missing {meta}", flush=True)
        return 2

    if not skip_plot:
        # Hist / hinge / spring / pile plots only — no movie frames (compare
        # script writes the peak snapshots). Full-window movies are huge.
        plot_py = (
            "import sys\n"
            f"sys.path.insert(0, r'{REPO / 'plot'}')\n"
            "import PlotEQ as peq\n"
            "peq.DO_FRAMES = -1\n"
            f"sys.argv = ['PlotEQ.py', r'{out}']\n"
            "raise SystemExit(peq.main())\n"
        )
        print(f"  PlotEQ (no frames): {out}", flush=True)
        pr = subprocess.run([PY, "-c", plot_py], cwd=str(REPO))
        if pr.returncode != 0:
            return pr.returncode
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--smoke",
        action="store_true",
        help="Short EQ (REGEN_EQ_TMAX=20, no free vib) to verify recorders",
    )
    ap.add_argument(
        "--only",
        nargs="+",
        default=None,
        help="Case slug(s) to run (default: all four MKR cases)",
    )
    ap.add_argument(
        "--skip-plot",
        action="store_true",
        help="Run OpenSees only; skip per-case PlotEQ",
    )
    ap.add_argument(
        "--list",
        action="store_true",
        help="Print case slugs and exit",
    )
    args = ap.parse_args()

    if args.list:
        for c in CASES:
            print(c["slug"])
        return 0

    wanted = {c["slug"] for c in CASES}
    if args.only:
        bad = [s for s in args.only if s not in wanted]
        if bad:
            raise SystemExit(f"unknown slug(s): {bad}\nknown: {sorted(wanted)}")
        cases = [c for c in CASES if c["slug"] in set(args.only)]
    else:
        cases = list(CASES)

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    for key in ("REGEN_OUTDIR",):
        if key in os.environ:
            _safe_outdir(os.environ[key])

    print(f"OpenSees exe: {find_opensees()}", flush=True)
    rc_all = 0
    for case in cases:
        rc = run_case(case, smoke=args.smoke, skip_plot=args.skip_plot)
        if rc != 0:
            rc_all = rc
            print(f"stopping after failure on {case['slug']}", flush=True)
            break
    return rc_all


if __name__ == "__main__":
    raise SystemExit(main())
