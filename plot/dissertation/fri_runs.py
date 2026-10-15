"""Friday (2026-08-21) main-campaign runs: metadata and a shared timing loader (dissertation Ch. 6; read-only on data).

Earthquake from t_int = 0 (gmStartTime 0), followed by the solitary wave near the end of the record.
stateOS at 2048 Hz on the laboratory clock (model scale); dt_sim = 10/2048 s = 4.88 ms; integration time t_int
(prototype) advances by dt_sim*sqrt(lambda) = 7.57 ms per target. Columns: typeConv3 (-1 init, 0 interpolate,
1 extrapolate, 2 slowdown), count (1..10 per dt_sim), s1 (new-target flag).
As on Wednesday, the start-up handshake raises two flags ~2 ms apart; intervals are counted from the second one.
Records are used up to the last new target (incomplete runs F07, F08, F22 stop early).
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]  # OSU_SSI_Bridge
sys.path.insert(0, str(ROOT / "plot"))
from lab_paths import MAT_EXTRACT_DIR, CYLINDER_LENGTH_SCALE as LAM, TIME_SCALE_FROUDE as SQL  # noqa: E402

DATA = ROOT / "OSU_SSI_BRIDGE_DATA"
FS, DT_SIM = 2048.0, 10 / 2048
T_D5, T_D95 = 56.0, 139.05  # Arias 5 % and 95 % of the ground motion, t_int (prototype s) from its start

# test, mat stem, dump, model, DOFs, integration processor, solver, ranks, note
RUNS = [
    ("F03", "0821_GusBridge_rowPos1", "r+01_20260821_0930_Storm_Wave", "Baseline", 2028, "CPU", "ProfileSPD", 8, ""),
    ("F04", "0821_GusBridge_rowPos1_16Core", "r+01_20260821_0949_Storm_Wave", "Baseline", 2028, "CPU", "ProfileSPD", 16, ""),
    ("F26", "0821_GusBridge_rowNeg2_8Core", "r-02_20260821_1601_Storm_Wave", "Baseline", 2028, "CPU", "MUMPS", 8, ""),
    ("F22", "0821_GusBridge_rowNeg6_Serial", "r-06_20260821_1528_Storm_Wave", "Baseline", 2028, "GPU", "cuDSS", 1, "incomplete"),
    ("F05", "0821_GusBridge_rowNeg4", "r-04_20260821_1010_Storm_Wave", "Baseline", 2028, "GPU", "cuDSS", 8, ""),
    ("F27", "0821_GusBridge_rowNeg4_8Core_P1", "r-04_20260821_1625_Storm_Wave", "Baseline", 2028, "GPU", "cuDSS", 8, ""),
    ("F07", "0821_GusBridge_rowNeg80", "r-80_20260821_1054_Storm_Wave", "Moderate", 3216, "CPU", "ProfileSPD", 8, "incomplete"),
    ("F08", "0821_GusBridge_rowNeg81", "r-81_20260821_1058_Storm_Wave", "Moderate", 3216, "CPU", "MUMPS", 8, "incomplete"),
    ("F06", "0821_GusBridge_rowNeg8", "r-08_20260821_1034_Storm_Wave", "Moderate", 3216, "GPU", "cuDSS", 8, ""),
    ("F17", "0821_GusBridge_rowNeg11_8Core", "r-11_20260821_1451_Storm_Wave", "Moderate", 3216, "GPU", "cuDSS hybrid", 8, ""),
    ("F09", "0821_GusBridge_rowNeg9", "r-09_20260821_1102_Storm_Wave", "Large", 4620, "GPU", "cuDSS", 8, ""),
    ("F10", "0821_GusBridge_rowNeg9_16Core", "r-09_20260821_1108_Storm_Wave", "Large", 4620, "GPU", "cuDSS", 16, ""),
    ("F12", "0821_GusBridge_rowNeg9_20Core", "r-09_20260821_1303_Storm_Wave", "Large", 4620, "GPU", "cuDSS", 20, ""),
    ("F11", "0821_GusBridge_rowNeg9_24Core", "r-09_20260821_1132_Storm_Wave", "Large", 4620, "GPU", "cuDSS", 24, ""),
    ("F14", "0821_GusBridge_rowNeg91_16Core", "r-91_20260821_1342_Storm_Wave", "Large", 4620, "CPU", "cuDSS", 16, ""),  # host WiFi on (G. Araujo)
    ("F15", "0821_GusBridge_rowNeg91_20Core", "r-91_20260821_1358_Storm_Wave", "Large", 4620, "CPU", "cuDSS", 20, ""),
    ("F13", "0821_GusBridge_rowNeg10_24Core_duplicate", "r-10_20260821_1312_Storm_Wave", "Extra-large", 5484, "GPU", "cuDSS", 24, ""),
    ("F20", "0821_GusBridge_rowNeg15_8Core", "r-15_20260821_1521_Storm_Wave", "Quad PIMY", 2028, "GPU", "cuDSS", 8, ""),
    ("F24", "0821_GusBridge_rowNeg14_8Core", "r-14_20260821_1544_Storm_Wave", "Elastic soil", 2028, "GPU", "cuDSS", 8, ""),
    ("F23", "0821_GusBridge_rowNeg13_8Core", "r-13_20260821_1535_Storm_Wave", "Stiff soil", 2028, "GPU", "cuDSS", 8, ""),
    ("F30", "0821_GusBridge_rowNeg17_8Core_P1", "r-17_20260821_1649_Storm_Wave", "2% Rayleigh", 2028, "GPU", "cuDSS", 8, ""),
]
# lunch dry runs: OpenSees + OpenFresco + rate transition running and commands generated, but the actuator did not
# move (measured displacement within 0.1 mm, load-cell force ~0) and no OpenSees dump; records end near t_int = 310-330 s,
# before the solitary wave. Timing only.
DRY = [
    ("Fd01", "0821_GusBridge_rowNeg8_16Core", "", "Moderate", 3216, "GPU", "cuDSS", 16, "dry"),
    ("Fd02", "0821_GusBridge_rowNeg9_18Core", "", "Large", 4620, "GPU", "cuDSS", 18, "dry"),
    ("Fd03", "0821_GusBridge_rowNeg9_19Core", "", "Large", 4620, "GPU", "cuDSS", 19, "dry"),
    ("Fd04", "0821_GusBridge_rowNeg9_22Core", "", "Large", 4620, "GPU", "cuDSS", 22, "dry"),
    ("Fd05", "0821_GusBridge_rowNeg10_20Core", "", "Extra-large", 5484, "GPU", "cuDSS", 20, "dry"),
    ("Fd06", "0821_GusBridge_rowNeg10_22Core", "", "Extra-large", 5484, "GPU", "cuDSS", 22, "dry"),
    ("Fd07", "0821_GusBridge_rowNeg10_24Core", "", "Extra-large", 5484, "GPU", "cuDSS", 24, "dry"),
    ("Fd08", "0821_GusBridge_rowNeg10_26Core", "", "Extra-large", 5484, "GPU", "cuDSS", 26, "dry"),
    ("Fd09", "0821_GusBridge_rowNeg10_28Core", "", "Extra-large", 5484, "GPU", "cuDSS", 28, "dry"),
]
BY_ID = {r[0]: r for r in RUNS + DRY}
FORCE_CACHE = Path(__file__).resolve().parent / "cache"


def load_force(stem):
    """Load-cell force record of one Simulink file: laboratory time (model s) and topLRaw column 0 (as recorded).
    Reads cache/force_<stem>.npz written by cache_force.py; falls back to the (large) .mat file."""
    c = FORCE_CACHE / f"force_{stem}.npz"
    if c.exists():
        z = np.load(c)
        return z["t"], z["f"]
    import scipy.io as sio
    raw = sio.loadmat(DATA / "Simulink" / f"{stem}.mat", squeeze_me=True, struct_as_record=False)["data"].topLRaw.data
    return raw[:, -1], raw[:, 0]


def label(r):
    tid, _, _, model, dofs, proc, solver, nr, _ = r
    integ = "GPU MKR" if proc == "GPU" else "CPU MKR"
    return f"{integ} + {solver}, " + ("serial" if nr == 1 else f"{nr} ranks")


def load(tid):
    """Timing record of one test: target arrivals, latency per window, states, slowdown episodes, delay."""
    r = BY_ID[tid]
    z = np.load(MAT_EXTRACT_DIR / f"{r[1]}.npz", allow_pickle=True)
    st = z["stateOS_data"]
    t, typ = st[:, -1], np.rint(st[:, 0]).astype(int)
    cnt, s1 = np.rint(st[:, 1]).astype(int), np.rint(st[:, 2])
    tk_all = t[1:][(s1[1:] > 0) & (s1[:-1] <= 0)]
    end = np.searchsorted(t, tk_all[-1] + DT_SIM)  # record ends one window after the last target
    t, typ, cnt = t[:end], typ[:end], cnt[:end]
    k0 = 1 if tk_all[1] - tk_all[0] < 3e-3 else 0  # second handshake flag
    tk = tk_all[k0:]
    win = t[1:][(cnt[1:] == 1) & (cnt[:-1] != 1)]  # dt_sim window starts
    i = np.searchsorted(win, tk) - 1
    ok = i >= 0
    lat = np.full(tk.size, np.nan)
    lat[ok] = tk[ok] - win[i[ok]]
    on = t[1:][(typ[1:] == 2) & (typ[:-1] != 2)]  # slowdown onsets
    off = t[1:][(typ[1:] != 2) & (typ[:-1] == 2)]
    if off.size < on.size:
        off = np.append(off, t[-1])
    dur = off[:on.size] - on
    t_su = tk[3] if tk.size > 3 else t[-1]  # end of the three start-up intervals
    ti_k = np.arange(tk.size) * DT_SIM * SQL  # integration time carried by each target (prototype s)
    use = typ != -1
    return dict(r=r, t=t, typ=typ, tk=tk, iv=np.diff(tk), lat=lat, on=on, dur=dur, t_su=t_su, ti_k=ti_k,
                frac={c: 100 * np.mean(typ[use] == c) for c in (0, 1, 2)},
                n_slow=int(np.sum(on > t_su)), n_slow_su=int(np.sum(on <= t_su)),
                to_int=lambda x: np.interp(x, tk, ti_k), complete=(r[8] not in ("incomplete", "dry")))


def window_delay(tid, t_su, t0, with_startup=False):
    """Delay behind real time from the dt_sim window starts (count == 1), model s, rebased at the end of start-up
    (t_su, after the three start-up intervals) or, with_startup, at the first target t0. Multiply by SQL for the
    prototype-scaled clock. Window starts carry no host-latency jitter, unlike the target arrivals."""
    st = np.load(MAT_EXTRACT_DIR / f"{BY_ID[tid][1]}.npz", allow_pickle=True)["stateOS_data"]
    t, cnt = st[:, -1], np.rint(st[:, 1]).astype(int)
    win = t[1:][(cnt[1:] == 1) & (cnt[:-1] != 1)]
    lag = (win - win[0]) - DT_SIM * np.arange(win.size)
    k = win >= (t0 if with_startup else t_su)
    return win[k], lag[k] - lag[k][0]


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    print(f"{'test':4} {'model':12} {'DOF':>5} {'config':30} {'lab_s':>6} {'t_int':>6} {'mean_iv':>7} {'lat_mean':>8}"
          f" {'lat_p50':>7} {'lat_p99':>7} {'extrap':>6} {'slow%':>6} {'n_slow':>6} {'+su':>3} {'slow_s':>7} {'delay':>7}")
    for r in RUNS + DRY:
        d = load(r[0])
        lat = d["lat"][3:] * 1e3
        lat = lat[np.isfinite(lat)]
        print(f"{r[0]:4} {r[3]:12} {r[4]:5d} {label(r):30} {d['t'][-1]:6.1f} {d['ti_k'][-1]:6.1f} "
              f"{1e3 * d['iv'][3:].mean():7.3f} {lat.mean():8.2f} {np.median(lat):7.2f} {np.percentile(lat, 99):7.2f} "
              f"{d['frac'][1]:6.1f} {d['frac'][2]:6.2f} {d['n_slow']:6d} {d['n_slow_su']:3d} {d['dur'].sum():7.2f} "
              f"{window_delay(r[0], d['t_su'], d['tk'][0])[1][-1] * SQL:7.4f}")  # prototype s, after start-up
