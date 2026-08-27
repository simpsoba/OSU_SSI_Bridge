#!/usr/bin/env python3
"""Regenerate eq_window.mp4 for recordersON=1 lab dumps (frames only)."""

from __future__ import annotations

import sys
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "plot"))

import PlotEQ as peq  # noqa: E402
import PlotEQParallel as pep  # noqa: E402

STATIC_SWITCHES = (
    "DO_HIST",
    "DO_DEPTH_HIST",
    "DO_ENVELOPE",
    "DO_SPRING_ENV",
    "DO_HYST",
    "DO_QUAD_PEAK",
    "DO_QUAD_HYST",
    "DO_HINGE",
    "DO_PILE_SEC",
)


def recorders_on_1(dump: Path) -> bool:
    meta = dump / "window_meta.txt.0"
    if not meta.is_file():
        meta = dump / "window_meta.txt"
    if not meta.is_file():
        return False
    for ln in meta.read_text(encoding="utf-8", errors="replace").splitlines():
        s = ln.strip()
        if s.startswith("recordersON"):
            parts = s.split()
            return len(parts) >= 2 and parts[1] == "1"
    return False


def runs_ok_in_log(log_path: Path) -> set[str]:
    """Runs that already finished with exit 0 in movie_regen.log."""
    if not log_path.is_file():
        return set()
    ok: set[str] = set()
    current: str | None = None
    for ln in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if ln.startswith("--- ") and ln.endswith(" ---"):
            current = ln[4:-4].strip()
        elif ln.startswith("exit ") and current is not None:
            if ln.strip() == "exit 0":
                ok.add(current)
            current = None
    return ok


def main() -> int:
    for name in STATIC_SWITCHES:
        setattr(peq, name, 0)
    peq.DO_FRAMES = 1
    peq.DO_FRAME_HIST = 1
    peq.N_FRAMES = 0
    peq.FRAME_FPS = 30

    root = REPO / "OSU_SSI_BRIDGE_DATA_LOCAL" / "opensees_data"
    plots_root = REPO / "OSU_SSI_BRIDGE_DATA_LOCAL" / "plots"
    log_path = plots_root / "movie_regen.log"
    all_runs = sorted(d.name for d in root.iterdir() if d.is_dir() and recorders_on_1(d))
    done = runs_ok_in_log(log_path)
    runs = [r for r in all_runs if r not in done]
    if not all_runs:
        print("RegenEqMovies: no recordersON=1 dumps")
        return 1
    if not runs:
        print(f"RegenEqMovies: all {len(all_runs)} recordersON=1 runs already ok")
        return 0
    if done:
        print(f"RegenEqMovies: skip ok {sorted(done)}")
    print(f"RegenEqMovies: remaining {runs}")

    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    with log_path.open("a", encoding="utf-8") as log:
        log.write(f"\n=== RegenEqMovies resume {stamp}  runs={runs}\n")
        rc_all = 0
        for run in runs:
            eq = root / run
            # Drop stale partial frame dumps from an interrupted pass.
            frames_dir = plots_root / "runs" / run / "eq" / "frames"
            if frames_dir.is_dir():
                n_old = sum(1 for _ in frames_dir.glob("frame_*.png"))
                if n_old > 12:
                    print(f"RegenEqMovies: clear stale frames ({n_old}) {frames_dir}", flush=True)
                    for p in frames_dir.glob("frame_*.png"):
                        p.unlink(missing_ok=True)
            print(f"\n===== {run} =====", flush=True)
            log.write(f"\n--- {run} ---\n")
            log.flush()
            sys.argv = ["RegenEqMovies.py", str(eq)]
            rc = pep.main()
            log.write(f"exit {rc}\n")
            log.flush()
            print(f"RegenEqMovies: {run} exit {rc}", flush=True)
            rc_all |= rc
        log.write(f"=== done rc={rc_all}\n")
    return rc_all


if __name__ == "__main__":
    raise SystemExit(main())
