#!/usr/bin/env python3
"""
Goals
-----
OpenSeesMP lumped-hold companions of the serial offline EQ dumps.
Does not write into the existing serial folders (lumped_hold_mkr,
lumped_hold_newmark, or any other slug already under eq_offline).

Five runs, 4 ranks, DT_FACTOR 10, first 180 s of the earthquake, no free
vibration. Each process is killed at 30 min.

  MKRAlphaExplicitMultiSOE 0.5 + Mumps
  MKRAlphaExplicitMultiSOE 0.5 + DistributedCuDSS
  CudaMKRAlpha 0.5          + DistributedCuDSS
  Newmark 0.5 0.25          + Mumps            (NormDispIncr 1e-6/50, recover 1e-5/50)
  Newmark 0.5 0.25          + DistributedCuDSS (same test)

OpenSeesMP is the 2026-10-05 RTHS-CUDA build.

Usage
-----
  python plot/RunEQOfflineMPLumpedHold.py
"""

from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OV_DIR = REPO / "analysis" / "eq_offline" / "overrides"
OUT_ROOT = REPO / "plot" / "out" / "eq_offline"
LOG_DIR = OUT_ROOT / "logs"

# Serial results. The recorder stage deletes every file in outDIR.
PROTECTED = {
    "lumped_hold_mkr",
    "lumped_hold_newmark",
    "lumped_free_mkr",
    "fbc_hold_mkr",
    "fbc_free_mkr",
    "fbc_free_newmark",
    "fbc_hold_newmark",
    "lumped_free_newmark",
    "compare",
    "logs",
}

OPENSEES_MP = Path(
    r"C:\Users\garaujor\source\repos\simpsoba\RTHS-CUDA\OpenSees\build\Release\OpenSeesMP.exe"
)
TCL_LIBRARY = Path(
    r"C:\Users\garaujor\source\repos\simpsoba\RTHS-CUDA\OpenSees\build\Release\lib\tcl8.6"
)
MPIEXEC = Path(r"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\bin\mpiexec.exe")
SETVARS = Path(r"C:\Program Files (x86)\Intel\oneAPI\setvars.bat")

NP = 4
WALL_LIMIT_S = 30 * 60
EQ_TMAX_S = 180.0

CASES = [
    {
        "slug": "mp4_lumped_hold_mkr_mumps",
        "eqIntegrator": "MKRAlphaExplicitMultiSOE 0.5",
        "postPartitionSystem": "Mumps",
        "newmark": False,
    },
    {
        "slug": "mp4_lumped_hold_mkr_cudss",
        "eqIntegrator": "MKRAlphaExplicitMultiSOE 0.5",
        "postPartitionSystem": "DistributedCuDSS",
        "newmark": False,
    },
    {
        "slug": "mp4_lumped_hold_cudamkr_cudss",
        "eqIntegrator": "CudaMKRAlpha 0.5",
        "postPartitionSystem": "DistributedCuDSS",
        "newmark": False,
    },
    {
        "slug": "mp4_lumped_hold_newmark_mumps",
        "eqIntegrator": "Newmark 0.5 0.25",
        "postPartitionSystem": "Mumps",
        "newmark": True,
    },
    {
        "slug": "mp4_lumped_hold_newmark_cudss",
        "eqIntegrator": "Newmark 0.5 0.25",
        "postPartitionSystem": "DistributedCuDSS",
        "newmark": True,
    },
]


def write_override(case: dict) -> Path:
    """Write one argv file. Refuses any slug that already holds a serial dump."""
    slug = case["slug"]
    if slug in PROTECTED or not slug.startswith("mp4_"):
        raise SystemExit(f"refusing to write protected or non-mp slug: {slug}")
    out_rel = f"plot/out/eq_offline/{slug}"
    out = (REPO / out_rel).resolve()
    if out.name in PROTECTED:
        raise SystemExit(f"refusing outDIR {out}")
    OV_DIR.mkdir(parents=True, exist_ok=True)
    dt_factor = int(case.get("dt_factor", 10))
    if case["newmark"]:
        test = (
            f"set DT_FACTOR {dt_factor}\n"
            "set eqTestTol 1.0e-6\n"
            "set eqTestIter 50\n"
            "set eqRecoverTol 1.0e-5\n"
            "set eqRecoverIter 50\n"
        )
    else:
        test = f"set DT_FACTOR {dt_factor}\n"
    text = f"""# OpenSeesMP lumped-hold, first {EQ_TMAX_S:.0f} s. Does not touch serial dumps.
# Dump: {out_rel}/
set runEQ 1
set plotFigures 0
set realTimeON 0
set recordersON 1
set eqPrintON 1
set eqPrintDt 10.0
set holdPierON 1
set holdPierRZON 0
set pierEleType "lumpedPlasticity"
set soilMesh 0
set soilProfile 4
set soilBoundary "Shin"
set soilEleType "SSPquad"
set soilConstitutive "inelastic"
set eqIntegrator "{case["eqIntegrator"]}"
set prePartitionSystem "UmfPack"
set postPartitionSystem "{case["postPartitionSystem"]}"
set constraintsHandler "Transformation"
{test}set eqTmax {EQ_TMAX_S}
set eqFreeVibT 0
set gmDir [file join $gmRoot Tohoku2011-FKSH]
set gmVelFile [file join $gmDir FKSH19.NS1.VT2]
set gmStartTime 0.0
set outDIR "{out_rel}"
"""
    path = OV_DIR / f"{slug}.tcl"
    path.write_text(text, encoding="utf-8", newline="\n")
    return path


def mpi_env() -> dict[str, str]:
    """Intel MPI pin plus the RTHS-CUDA Tcl library. One thread per rank.

    impi.dll sits next to mpiexec and is not on PATH, so the ranks exit
    0xC0000135 (DLL not found) unless that directory is prepended.
    """
    env = os.environ.copy()
    env["TCL_LIBRARY"] = str(TCL_LIBRARY)
    env["I_MPI_PIN"] = "on"
    env["I_MPI_PIN_CELL"] = "core"
    env["I_MPI_FABRICS"] = "shm"
    env["OMP_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    mpi_bin = str(MPIEXEC.parent)
    env["PATH"] = mpi_bin + os.pathsep + env.get("PATH", "")
    return env


def kill_tree(pid: int, slug: str = "") -> None:
    """Stop mpiexec and the ranks it spawned.

    Intel MPI's hydra_pmi_proxy can outlive mpiexec and keep recorder
    files open, so also stop any OpenSeesMP still running this slug.
    """
    subprocess.run(
        ["taskkill", "/F", "/T", "/PID", str(pid)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    if not slug:
        return
    ps = (
        "Get-CimInstance Win32_Process -Filter \"Name='OpenSeesMP.exe'\" | "
        f"Where-Object {{ $_.CommandLine -like '*{slug}*' }} | "
        "ForEach-Object { taskkill /F /T /PID $_.ParentProcessId }"
    )
    subprocess.run(
        ["powershell", "-NoProfile", "-Command", ps],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )


def clear_out(out: Path) -> None:
    """Drop leftover recorder files so EQRecorders can recreate the folder."""
    out.mkdir(parents=True, exist_ok=True)
    for path in out.iterdir():
        if path.is_dir():
            continue
        path.unlink()


def run_case(case: dict) -> dict:
    """Run one case. Returns a short record (exit, wall, timed_out)."""
    ov = write_override(case)
    slug = case["slug"]
    out = OUT_ROOT / slug
    clear_out(out)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"{slug}.log"
    # setvars.bat only affects the cmd process that calls it. A .bat file
    # avoids re-quoting that path when Python builds the CreateProcess line.
    bat = LOG_DIR / f"{slug}.bat"
    bat.write_text(
        "\r\n".join(
            [
                "@echo off",
                f'call "{SETVARS}" intel64 >nul 2>&1',
                "if errorlevel 1 exit /b 1",
                "set OMP_NUM_THREADS=1",
                "set MKL_NUM_THREADS=1",
                "set I_MPI_PIN=on",
                "set I_MPI_PIN_CELL=core",
                "set I_MPI_FABRICS=shm",
                (
                    f'"{MPIEXEC}" -n {NP} "{OPENSEES_MP}" '
                    f'"{REPO / "RunParallel.tcl"}" "{ov}"'
                ),
                "exit /b %ERRORLEVEL%",
                "",
            ]
        ),
        encoding="utf-8",
    )
    cmd = ["cmd", "/d", "/c", str(bat)]
    print(f"=== {slug} ===", flush=True)
    print(f"  {' '.join(cmd)}", flush=True)
    print(f"  log: {log_path}", flush=True)
    t0 = time.time()
    timed_out = False
    with log_path.open("w", encoding="utf-8", errors="replace") as fd:
        fd.write(f"# cmd: {cmd}\n# cwd: {REPO}\n")
        fd.flush()
        proc = subprocess.Popen(
            cmd,
            cwd=str(REPO),
            env=mpi_env(),
            stdout=fd,
            stderr=subprocess.STDOUT,
        )
        try:
            proc.wait(timeout=WALL_LIMIT_S)
        except subprocess.TimeoutExpired:
            timed_out = True
            kill_tree(proc.pid, slug)
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                kill_tree(proc.pid, slug)
    wall = time.time() - t0
    code = proc.returncode
    print(
        f"  exit={code}  wall={wall / 60.0:.1f} min  timeout={timed_out}",
        flush=True,
    )
    return {
        "slug": slug,
        "newmark": case["newmark"],
        "code": code,
        "wall_s": wall,
        "timed_out": timed_out,
        "log": log_path,
        "out": out,
    }


def log_text(path: Path) -> str:
    if not path.is_file():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def reached_eq(text: str) -> bool:
    """True once rank 0 has printed the earthquake header."""
    return "----- EQ (OpenSeesMP" in text


def finished_180(text: str) -> bool:
    """True when the driver printed EQ done after a 180 s record."""
    return "EQ done" in text and "EQ 180" in text


def dll_missing(text: str) -> bool:
    return "c0000135" in text.lower() or "STATUS_DLL_NOT_FOUND" in text


def main() -> int:
    for path in (OPENSEES_MP, TCL_LIBRARY, MPIEXEC, SETVARS):
        if not path.exists():
            raise SystemExit(f"missing {path}")
    rc = 0
    newmark_retry: list[dict] = []
    for case in CASES:
        rec = run_case(case)
        text = log_text(rec["log"])
        if dll_missing(text) or not reached_eq(text):
            raise SystemExit(
                f"{rec['slug']}: OpenSeesMP did not enter the EQ. See {rec['log']}"
            )
        if case["newmark"] and not finished_180(text):
            print(f"  {rec['slug']} did not finish 180 s — will retry DT_FACTOR 5", flush=True)
            newmark_retry.append(case)
        elif rec["code"] not in (0, None) and not rec["timed_out"]:
            print(f"  {rec['slug']} returned {rec['code']}", flush=True)
            rc = rec["code"] or rc
    for case in newmark_retry:
        retry = dict(case)
        retry["slug"] = case["slug"] + "_dt5"
        retry["dt_factor"] = 5
        print(f"  retry {retry['slug']}", flush=True)
        rec = run_case(retry)
        text = log_text(rec["log"])
        if dll_missing(text):
            raise SystemExit(f"{rec['slug']}: OpenSeesMP did not start (missing DLL). See {rec['log']}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
