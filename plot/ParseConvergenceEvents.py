#!/usr/bin/env python3
"""
Goals
-----
Scan the offline EQ logs once and save every convergence event, so plots can
mark them without re-reading the logs.

Events (one row each, convergence_events.csv):
  fail         "the Algorithm failed at time T", with the norm / tolerance of
               the CTestNormDispIncr report just before it (one failed attempt)
  recover_tol  "recover NormDispIncr TOL at step N  t~T" (relaxed tolerance)
  recover_dt   "recover dt=DT s at step N  t~T"          (dt/4 sub-steps)

Bad steps (one row each, convergence_steps.csv), classified by the rung of the
RunParallel.tcl / Run.tcl recovery ladder that finally got the step through:
  tol    retry at dt with the relaxed tolerance converged
  dt2    relaxed tolerance at dt failed; two dt/2 sub-steps converged
         (silent in the log: inferred from a relaxed-tol failure with no
         "recover dt" line after it)
  dt4    dt/2 failed too; four dt/4 sub-steps converged
  abort  last bad step of a run whose ladder was exhausted (the driver's
         "analyze failed at step" error is in the log)

A block starts at "recover at step N  t~T". Failures reported against the
relaxed tolerance belong to that block; a failure against the normal tolerance
is the first attempt of the next bad step and closes the block. t_s is the
committed time at the start of the bad step.

Writes
------
  plot/out/eq_offline/logs/convergence_events.csv
      run, event, t_s, step, norm, tol, dt_s
  plot/out/eq_offline/logs/convergence_steps.csv
      run, step, t_s, outcome, n_fail
  plot/out/eq_offline/logs/convergence_summary.csv
      run, completed, n_fail, n_recover_tol, n_recover_dt, n_bad_steps,
      n_tol, n_dt2, n_dt4, n_abort, first_fail_s, last_event_s

``completed`` is True when the log has "EQ done". Runs with no events still get
a summary row.

Usage
-----
  python plot/ParseConvergenceEvents.py [logs_dir]
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

LOGS = Path(__file__).resolve().parent / "out" / "eq_offline" / "logs"

RE_NORM = re.compile(r"after:\s*(\d+) iterations\s+current Norm:\s*([0-9.eE+-]+)\s*\(max:\s*([0-9.eE+-]+)")
RE_FAIL = re.compile(r"Algorithm failed at time\s+([0-9.eE+-]+)")
RE_START = re.compile(r"recover at step\s+(\d+)\s+t~([0-9.eE+-]+)")
RE_TOL = re.compile(r"recover NormDispIncr\s+([0-9.eE+-]+)\s+at step\s+(\d+)\s+t~([0-9.eE+-]+)")
RE_DT = re.compile(r"recover dt=([0-9.eE+-]+)\s*s\s+at step\s+(\d+)\s+t~([0-9.eE+-]+)")


def parse_log(path: Path) -> tuple[list[dict], list[dict], bool]:
    """Events in file order, bad steps in file order, and whether the EQ stage finished."""
    events: list[dict] = []
    steps: list[dict] = []
    norm = tol = ""
    done = False
    aborted = False
    block: dict | None = None
    relaxed = ""  # relaxed tolerance of the current block, as printed

    def close() -> None:
        nonlocal block
        if block is not None:
            if block["cut_dt4"]:
                block["outcome"] = "dt4"
            elif block["n_retry_fail"] > 0:
                block["outcome"] = "dt2"
            else:
                block["outcome"] = "tol"
            steps.append(block)
            block = None

    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "EQ done" in line:
            done = True
        if "analyze failed at step" in line:
            aborted = True
        m = RE_NORM.search(line)
        if m:
            norm, tol = m.group(2), m.group(3)
            continue
        m = RE_FAIL.search(line)
        if m:
            events.append({"event": "fail", "t_s": m.group(1), "step": "", "norm": norm, "tol": tol, "dt_s": ""})
            if block is not None:
                if relaxed and tol and float(tol) == float(relaxed):
                    block["n_fail"] += 1
                    block["n_retry_fail"] += 1
                else:
                    close()  # normal-tolerance failure: next bad step's first attempt
            norm = tol = ""
            continue
        m = RE_START.search(line)
        if m:
            close()
            # The normal-tolerance failure that triggered this block was just printed.
            block = {"step": m.group(1), "t_s": m.group(2), "n_fail": 1, "n_retry_fail": 0, "cut_dt4": False}
            continue
        m = RE_TOL.search(line)
        if m:
            relaxed = m.group(1)
            events.append({"event": "recover_tol", "t_s": m.group(3), "step": m.group(2), "norm": "", "tol": m.group(1), "dt_s": ""})
            continue
        m = RE_DT.search(line)
        if m:
            events.append({"event": "recover_dt", "t_s": m.group(3), "step": m.group(2), "norm": "", "tol": "", "dt_s": m.group(1)})
            if block is not None:
                block["cut_dt4"] = True
    close()
    # Ladder exhausted: the driver raised its "analyze failed at step" error. A run
    # that merely stops (killed, wall-clock cap) keeps its last outcome.
    if steps and not done and aborted:
        steps[-1]["outcome"] = "abort"
    return events, steps, done


def main(argv: list[str]) -> int:
    logs = Path(argv[0]) if argv else LOGS
    rows: list[dict] = []
    step_rows: list[dict] = []
    summary: list[dict] = []
    for path in sorted(logs.glob("*.log")):
        if path.name.startswith("_"):
            continue
        run = path.stem
        events, steps, done = parse_log(path)
        for e in events:
            rows.append({"run": run, **e})
        for s in steps:
            step_rows.append({"run": run, "step": s["step"], "t_s": s["t_s"], "outcome": s["outcome"], "n_fail": s["n_fail"]})
        fails = [float(e["t_s"]) for e in events if e["event"] == "fail"]
        outcomes = [s["outcome"] for s in steps]
        summary.append({
            "run": run,
            "completed": done,
            "n_fail": len(fails),
            "n_recover_tol": sum(e["event"] == "recover_tol" for e in events),
            "n_recover_dt": sum(e["event"] == "recover_dt" for e in events),
            "n_bad_steps": len(steps),
            "n_tol": outcomes.count("tol"),
            "n_dt2": outcomes.count("dt2"),
            "n_dt4": outcomes.count("dt4"),
            "n_abort": outcomes.count("abort"),
            "first_fail_s": f"{min(fails):.4f}" if fails else "",
            "last_event_s": f"{max(float(e['t_s']) for e in events):.4f}" if events else "",
        })

    ev_path = logs / "convergence_events.csv"
    with ev_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["run", "event", "t_s", "step", "norm", "tol", "dt_s"])
        w.writeheader()
        w.writerows(rows)
    st_path = logs / "convergence_steps.csv"
    with st_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["run", "step", "t_s", "outcome", "n_fail"])
        w.writeheader()
        w.writerows(step_rows)
    sm_path = logs / "convergence_summary.csv"
    with sm_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        w.writeheader()
        w.writerows(summary)

    for s in summary:
        print(
            f"{s['run']:36s} done={str(s['completed']):5s} attempts_failed={s['n_fail']:4d} "
            f"bad_steps={s['n_bad_steps']:4d} (tol {s['n_tol']:4d}, dt/2 {s['n_dt2']:3d}, "
            f"dt/4 {s['n_dt4']:3d}, abort {s['n_abort']})"
        )
    print(f"wrote {ev_path}")
    print(f"wrote {st_path}")
    print(f"wrote {sm_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
