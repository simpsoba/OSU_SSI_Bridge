"""Shared loader: new-target times, state series, tracking cut (Wednesday runs; dissertation Ch. 6; read-only)."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
import ctrl_error_wed as C  # noqa: E402


def timing(stem):
    z = np.load(C.MAT_EXTRACT_DIR / f"{stem}.npz", allow_pickle=True)
    st, ta, me = z["stateOS_data"], z["tarSigOS_data"], z["meaSigOS_data"]
    n = min(len(st), len(ta), len(me))
    t, s1, typ = st[:n, -1], np.rint(st[:n, 2]), np.rint(st[:n, 0]).astype(int)
    rms = np.sqrt(np.convolve((ta[:n, 0] - me[:n, 0]) ** 2, np.ones(512) / 512, mode="same"))
    bad = np.flatnonzero((rms > 2e-3) & (t > 50.0))
    end = max(bad[0] - 256, 0) if bad.size else n
    t, s1, typ = t[:end], s1[:end], typ[:end]
    tk = t[1:][(s1[1:] > 0) & (s1[:-1] <= 0)]
    k_first = 1 if tk[1] - tk[0] < 3e-3 else 0  # W02/W04: one 2-ms interval precedes the long first step
    return dict(t=t, typ=typ, tk=tk, first=tk[k_first + 1] - tk[k_first], dt=np.diff(tk[k_first + 4:]),
                startup=np.diff(tk[k_first:k_first + 4]))


if __name__ == "__main__":
    for wid, stem, dump, nr in C.RUNS:
        d = timing(stem)
        q = np.percentile(d["dt"] * 1e3, [1, 25, 50, 75, 99, 99.9])
        use = d["typ"] != -1
        fr = [100 * np.mean(d["typ"][use] == c) for c in (0, 1, 2)]
        print(f"{wid} {nr:2d}: first {d['first']*1e3:6.1f} ms  startup {np.round(d['startup']*1e3,1)}  intervals n={d['dt'].size}"
              f"  pct(1,25,50,75,99,99.9)={np.round(q,2)}  max {d['dt'].max()*1e3:.1f}  mean {d['dt'].mean()*1e3:.3f}"
              f"  >dt_sim+1 {100*np.mean(d['dt']>C.DT_SIM+1e-3):.2f}%  | interp {fr[0]:.2f}% extrap {fr[1]:.2f}% slow {fr[2]:.3f}%")
