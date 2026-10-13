"""One-time cache of the load-cell force (topLRaw) of the Simulink records (dissertation Ch. 6; read-only on data).

The Simulink .mat files are large once decompressed; the dissertation figures only need topLRaw (column 0 = load-cell
force, model scale; last column = laboratory time, model s). This writes cache/force_<mat stem>.npz with t (float64)
and f (float64, as recorded) for each run, one at a time, skipping files already cached.
  python plot/dissertation/cache_force.py [stem ...]     # default: every Friday run and the five Wednesday runs
fri_runs.load_force(stem) reads the cache (and falls back to the .mat file).
"""
import gc
import sys
from pathlib import Path

import numpy as np
import scipy.io as sio

sys.path.insert(0, str(Path(__file__).parent))
import fri_runs as R  # noqa: E402

WED = ["0819_GusBridge_Serial_H0p5_T7p746_reg_Trial02", "0819_GusBridge_Parallel4_H0p5_T7p746_reg_Trial02",
       "0819_GusBridge_Parallel8_H0p5_T7p746_reg", "0819_GusBridge_Parallel16_H0p5_T7p746_reg_Trial02",
       "0819_GusBridge_Parallel24_H0p5_T7p746_reg"]

if __name__ == "__main__":
    stems = sys.argv[1:] or [r[1] for r in R.RUNS] + WED
    R.FORCE_CACHE.mkdir(exist_ok=True)
    for stem in stems:
        out = R.FORCE_CACHE / f"force_{stem}.npz"
        if out.exists():
            print("cached  ", stem)
            continue
        try:
            raw = sio.loadmat(R.DATA / "Simulink" / f"{stem}.mat", squeeze_me=True, struct_as_record=False)["data"].topLRaw.data
        except MemoryError:
            print("MEMORY  ", stem, "(retry later)")
            gc.collect()
            continue
        np.savez(out, t=raw[:, -1].astype(np.float64), f=raw[:, 0].astype(np.float64))
        print("wrote   ", stem, raw.shape)
        del raw
        gc.collect()
