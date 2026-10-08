#!/usr/bin/env python3
"""
Goals
-----
One-shot regen of campaign post-process figures under
``OSU_SSI_BRIDGE_DATA_LOCAL/plots/runs/*/os/`` and ``plots/compare/``.

Does **not** replot ``eq/`` (PlotEQ*) or paper sketches.

Clocks: lab ``$t$``; OpenSees ``$t_{int}$``; OS→lab via atTarget map
(``lab_time_map.py``). See ``plot/lab/LAB_RUN_MAP.md`` § Clocks and map.

  python plot/RegenCampaignOsPlots.py
  python plot/RegenCampaignOsPlots.py --sync
  python plot/RegenCampaignOsPlots.py --compare-only
  python plot/RegenCampaignOsPlots.py W04 F08

Args:
  --sync          run SyncLabBackup.py first (ingest + archive + mat extract)
  --runs-only     per-Test OS suite only
  --compare-only  hist/PSD pairs + stateOS bars only
  Test IDs        optional subset; default = all CSV rows with MatFile+DumpFolder
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from lab_paths import load_lab_runs_rows
from paths import HERE

REPO = HERE.parent
PY = sys.executable
PLOT = HERE


def mat_dump_tests() -> list[str]:
    """Test IDs with both MatFile and DumpFolder."""
    out: list[str] = []
    for row in load_lab_runs_rows():
        tid = (row.get("Test") or "").strip()
        mat = (row.get("MatFile") or "").strip()
        dump = (row.get("DumpFolder") or "").strip()
        if tid and mat and dump:
            out.append(tid)
    return out


def run_script(script: str, args: list[str] | None = None) -> int:
    """
    Run ``plot/<script>`` with the current interpreter.

    Args:    script  basename; args  CLI after script name
    Returns: process return code
    """
    cmd = [PY, str(PLOT / script)] + (args or [])
    print(f"\n=== {' '.join(cmd)} ===", flush=True)
    return int(subprocess.call(cmd, cwd=str(REPO)))


def parse_argv(argv: list[str]) -> tuple[list[str], bool, bool, bool]:
    """
    Split flags and Test IDs.

    Returns: (tests, do_sync, runs_only, compare_only)
    """
    tests: list[str] = []
    do_sync = False
    runs_only = False
    compare_only = False
    for a in argv:
        if a in ("-h", "--help"):
            print(__doc__)
            raise SystemExit(0)
        if a == "--sync":
            do_sync = True
        elif a == "--runs-only":
            runs_only = True
        elif a == "--compare-only":
            compare_only = True
        elif a.startswith("-"):
            raise SystemExit(f"RegenCampaignOsPlots: unknown option {a}")
        else:
            tests.append(a)
    if runs_only and compare_only:
        raise SystemExit("RegenCampaignOsPlots: use only one of --runs-only / --compare-only")
    return tests, do_sync, runs_only, compare_only


def main() -> int:
    tests, do_sync, runs_only, compare_only = parse_argv(sys.argv[1:])
    if not tests:
        tests = mat_dump_tests()
    if not tests and not compare_only:
        print("RegenCampaignOsPlots: no mat+dump Test IDs", file=sys.stderr)
        return 1

    rc = 0
    if do_sync:
        rc = max(rc, run_script("SyncLabBackup.py"))

    if not compare_only:
        print(f"RegenCampaignOsPlots: {len(tests)} Test ID(s)", flush=True)
        for script in (
            "PlotActuatorForce.py",
            "PlotPierBaseForce.py",
            "PlotActuatorVsPier.py",
            "PlotActuatorHyst.py",
            "PlotHydroSpectra.py",  # structure + SSI + psd_frc_actuator ($t_int$)
            "PlotLabTimeMapDiag.py",
        ):
            rc = max(rc, run_script(script, tests))
        rc = max(rc, run_script("PlotHydroSpectra.py", ["--lab-vs-pier-os"] + tests))
        rc = max(rc, run_script("PlotMatOS.py"))

    if not runs_only:
        rc = max(rc, run_script("PlotEQComparePairs.py"))
        rc = max(rc, run_script("PlotEQComparePSD.py"))
        rc = max(rc, run_script("PlotStateOSBars.py"))

    print(f"\nRegenCampaignOsPlots: done (rc={rc})", flush=True)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
