"""Control error of the Wednesday runs by phase (dissertation Ch. 6; read-only on data).

Phases on integrator time t_int (prototype): before EQ < 150 s (waves only); during 150 s to the Arias 95%
time (289.1 s; Arias 5% at 206.0 s); after > 289.1 s. 1-rank record cut where its actuator stops tracking.
NRMSE = RMS(y - x) / (max x - min x) within the phase (Chae et al. 2013 / Dong et al. 2015 convention), x the
reference (target or compensated command), y measured. Delays: least-squares shift (see xdelay), model scale.
  command -> measured   comSigOS vs meaSigOS (actuator + AFC residual)
  target -> measured    tarSigOS vs meaSigOS
  integrator -> target  ctrlDsp/lambda placed at hard real time vs tarSigOS on the laboratory clock: rate
                        transition plus lag behind real time. Hard real time is referenced to the arrival of the
                        row at t_int = 5 s (row r due at t_r0 + (r - r0) dt_sim), which removes the start-up offset
                        (first targets 40-290 ms late at t < 0.3 s, model at rest); the offset is reported apart.
Row r of ctrlDsp is carried by new target r + 1; tarSigOS reaches row r when target r + 2 arrives (checked).
"""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]  # OSU_SSI_Bridge
sys.path.insert(0, str(ROOT / "plot"))
from lab_paths import MAT_EXTRACT_DIR, CYLINDER_LENGTH_SCALE as LAM, TIME_SCALE_FROUDE as SQL  # noqa: E402
DATA = ROOT / "OSU_SSI_BRIDGE_DATA"
RUNS = [("W04", "0819_GusBridge_Serial_H0p5_T7p746_reg_Trial02", "0819_GusBridge_Serial_H0p5_T7p746_reg_Trial02", 1),
        ("W06", "0819_GusBridge_Parallel4_H0p5_T7p746_reg_Trial02", "0819_GusBridge_Parallel4_H0p5_T7p746_reg_Trial02", 4),
        ("W02", "0819_GusBridge_Parallel8_H0p5_T7p746_reg", "0819_GusBridge_Parallel8_H0p5_T7p746_reg", 8),
        ("W05", "0819_GusBridge_Parallel16_H0p5_T7p746_reg_Trial02", "0819_GusBridge_Parallel16_H0p5_T7p746_reg", 16),
        ("W07", "0819_GusBridge_Parallel24_H0p5_T7p746_reg", "0819_GusBridge_Parallel24_H0p5_T7p746_reg", 24)]
FS, DT_SIM = 2048.0, 10 / 2048
T_EQ0, T_EQ95 = 150.0, 150.0 + 139.05
T_REF = 5.0  # s, integrator time at which hard real time is referenced
PHASES = [("before", -1.0, T_EQ0), ("during", T_EQ0, T_EQ95), ("after", T_EQ95, 1e9)]


def load(stem, dump):
    z = np.load(MAT_EXTRACT_DIR / f"{stem}.npz", allow_pickle=True)
    st, ta, co, me = z["stateOS_data"], z["tarSigOS_data"], z["comSigOS_data"], z["meaSigOS_data"]
    n = min(len(st), len(ta), len(co), len(me))
    t, s1 = st[:n, -1], np.rint(st[:n, 2])
    tk = t[1:][(s1[1:] > 0) & (s1[:-1] <= 0)]
    ta, co, me = ta[:n, 0], co[:n, 0], me[:n, 0]
    rms = np.sqrt(np.convolve((ta - me) ** 2, np.ones(512) / 512, mode="same"))
    bad = np.flatnonzero((rms > 2e-3) & (t > 50.0))
    end = max(bad[0] - 256, 0) if bad.size else n
    c = np.loadtxt(DATA / dump / "Elmt101_ctrlDsp.out")
    c = c[:min(len(c), tk.size - 1)]
    ti = (np.interp(t, tk, np.arange(tk.size)) - 1) * DT_SIM * SQL  # integrator time of the row reached
    r0 = int(round(T_REF / (DT_SIM * SQL)))
    th = tk[r0 + 1] + (np.arange(len(c)) - r0) * DT_SIM  # hard real time of row r (model, laboratory clock)
    lag = tk[1:len(c) + 1] - th  # arrival lag of each row behind hard real time (s)
    ui = np.interp(t, th, c[:, 1] / LAM, right=np.nan)  # integrator displacement on hard real time, resampled at 2048 Hz
    keep = (t > 1.0) & (np.arange(n) < end)
    tr = np.arange(len(c)) * DT_SIM * SQL
    kr = tk[1:len(c) + 1] <= t[keep][-1]
    return dict(t=t[keep], ti=ti[keep], ta=ta[keep], co=co[keep], me=me[keep], ui=ui[keep],
                tr=tr[kr], lag=lag[kr], offset=(tk[r0 + 1] - tk[1]) - r0 * DT_SIM)


def xdelay(x, y, maxlag_s=0.2):
    """Lag (ms) by which y trails x: the shift minimizing mean (y[k + L] - x[k])^2 over a fixed overlap
    (equivalent to the peak of the cross-correlation normalized for overlap and energy), parabolic sub-sample
    refinement. The raw cross-correlation peak is flat and biased on these nonstationary, low-frequency-dominated
    records (8 ranks, during: 3.5 ms on the whole phase, 0.1 ms on 10-s windows, against 4.9 ms here)."""
    M = int(maxlag_s * FS)
    n = x.size - 2 * M
    xs = x[M:M + n]
    L = np.arange(-M, M + 1)
    e = np.array([np.mean((y[M + l:M + l + n] - xs) ** 2) for l in L])
    i = int(np.argmin(e))
    d = 0.0 if i in (0, L.size - 1) else 0.5 * (e[i - 1] - e[i + 1]) / (e[i - 1] - 2 * e[i] + e[i + 1])
    return 1e3 * (L[i] + d) / FS


def nrmse(x, y):
    return 100 * np.sqrt(np.mean((y - x) ** 2)) / (x.max() - x.min())


if __name__ == "__main__":
    print(f"{'run':4} {'ranks':>5} {'phase':7} {'dur_s':>6} {'rng_mm':>7} {'RMS_tm':>7} {'NRMSE_tm':>8} {'NRMSE_cm':>8}"
          f" {'d_cm':>6} {'d_tm':>6} {'d_it':>8} {'lag_med':>8} {'lag_max':>8}")
    for wid, stem, dump, nr in RUNS:
        d = load(stem, dump)
        print(f"{wid} start-up offset {1e3 * d['offset']:.1f} ms")
        for ph, a, b in PHASES:
            m = (d["ti"] >= a) & (d["ti"] < b)
            if m.sum() < 4 * FS:
                continue
            x, y, cm, ui = d["ta"][m], d["me"][m], d["co"][m], d["ui"][m]
            mr = (d["tr"] >= a) & (d["tr"] < b)
            ok = np.isfinite(ui)
            print(f"{wid:4} {nr:5d} {ph:7} {m.sum() / FS:6.1f} {1e3 * (x.max() - x.min()):7.1f} "
                  f"{1e3 * np.sqrt(np.mean((y - x) ** 2)):7.3f} {nrmse(x, y):7.2f}% {nrmse(cm, y):7.2f}% "
                  f"{xdelay(cm, y):6.2f} {xdelay(x, y):6.2f} {(xdelay(ui[ok], x[ok]) if nr > 1 else np.nan):8.2f} "
                  f"{1e3 * np.median(d['lag'][mr]):8.2f} {1e3 * d['lag'][mr].max():8.2f}")
    print("dur = phase duration on the laboratory clock (model s); rng = target range (model mm); RMS_tm (model mm);"
          " delays d_* in model ms (cm command->measured, tm target->measured, it integrator->target);"
          " lag_med/lag_max = arrival lag of new targets behind hard real time (model ms).")
