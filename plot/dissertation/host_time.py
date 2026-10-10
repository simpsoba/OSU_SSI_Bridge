"""Host turnaround per step from stateOS: time from the start of each dt_sim window (count resets to 1) to the
arrival of the new target (s1 rising edge). Read-only."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np, ctrl_error_wed as C
from wed_timing_stats import timing


def host_times(stem):
    z = np.load(C.MAT_EXTRACT_DIR / f"{stem}.npz", allow_pickle=True)
    st = z["stateOS_data"]
    d = timing(stem)
    n = d["t"].size
    t, cnt, s1 = st[:n, -1], np.rint(st[:n, 1]).astype(int), np.rint(st[:n, 2])
    win = t[1:][(cnt[1:] == 1) & (cnt[:-1] != 1)]          # window starts
    arr = t[1:][(s1[1:] > 0) & (s1[:-1] <= 0)]              # target arrivals
    i = np.searchsorted(win, arr) - 1                       # latest window start before each arrival
    ok = i >= 0
    h = arr[ok] - win[i[ok]]
    return arr[ok], h, cnt


if __name__ == "__main__":
    med1 = None
    for wid, stem, dump, nr in C.RUNS:
        ta, h, cnt = host_times(stem)
        h = h[ta > 2.0] * 1e3
        q = np.percentile(h, [5, 25, 50, 75, 95, 99])
        med1 = med1 or np.median(h)
        print(f"{wid} {nr:2d}: n={h.size}  host time pct(5,25,50,75,95,99) ms = {np.round(q,2)}  mean {h.mean():.2f}  "
              f"speedup(median) {med1/np.median(h):.2f}  speedup(mean) vs 1: -  | count values seen {np.unique(cnt)[:12]}")
