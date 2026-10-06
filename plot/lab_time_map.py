#!/usr/bin/env python3
"""
Goals
-----
Map OpenSees recorder time onto the Simulink lab clock for RTHS plots:

  t_lab = k + t_OS / √λ + f

- k = first typeConv3 1→0 → t_ex − Δt_con − Δt_sim (startup waits inside k).
- f = cumulative mid-run typeConv3==2 only (after t_ex).
- Lab Time is ground truth; never place lab events with naive t_lab·√λ.

See plot/lab/MKR_OPENFRESCO_TIMING.md §4.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from lab_paths import MAT_EXTRACT_DIR, TIME_SCALE_FROUDE

N_COUNT = 10  # Δt_sim = N_COUNT · Δt_con
SLOWDOWN_STATE = 2  # typeConv3 wait flag (OpenFresco)

LABEL_T_LAB = r"$t$ (s) lab / model scale"
LABEL_T_PROTO_TOP = r"$t\sqrt{\lambda}$ (s) prototype scale"


def load_type_conv3(mat_name: str) -> tuple[np.ndarray, np.ndarray] | None:
    """
    Lab Time and integer typeConv3 from a mat extract.

    Args:    mat_name  Simulink .mat file name
    Returns: (t_lab_s, typeConv3) or None
    """
    npz_path = MAT_EXTRACT_DIR / f"{Path(mat_name).stem}.npz"
    if not npz_path.is_file():
        return None
    z = np.load(npz_path, allow_pickle=True)
    if "stateOS_time" in z.files and "stateOS_data" in z.files:
        t = np.asarray(z["stateOS_time"], dtype=float)
        tc3 = np.rint(np.asarray(z["stateOS_data"][:, 0], dtype=float)).astype(int)
        return t, tc3
    return None


@dataclass(frozen=True)
class LabTimeMap:
    """Mapped OS samples plus clock knobs."""

    t_lab: np.ndarray
    k: float
    t_ex: float
    f_end: float
    f_at: np.ndarray


def lab_clock_offset_s(mat_name: str) -> tuple[float, float]:
    """
    Lab-clock offset k and first mid-run execute time t_ex (model s).

    Args:    mat_name  MatFile basename (with or without .mat)
    Returns: (k, t_ex) lab s
    """
    z = np.load(MAT_EXTRACT_DIR / f"{Path(mat_name).stem}.npz", allow_pickle=True)
    t = np.asarray(z["stateOS_time"], dtype=float)
    tc3 = np.asarray(z["stateOS_data"][:, 0], dtype=int)
    dt_con = float(np.median(np.diff(t)))
    dt_sim = float(N_COUNT) * dt_con
    land = np.concatenate([[False], (tc3[:-1] == 1) & (tc3[1:] == 0)])
    lands = np.flatnonzero(land)
    if lands.size == 0:
        raise RuntimeError(f"no typeConv3 1→0 landing in {mat_name}")
    t_land = float(t[int(lands[0])])
    prev = np.concatenate([[int(tc3[0])], tc3[:-1]])
    ex = np.flatnonzero((tc3 == 1) & (prev == 0) & (t > t_land))
    if ex.size == 0:
        raise RuntimeError(f"no typeConv3 0→1 after landing in {mat_name}")
    t_ex = float(t[int(ex[0])])
    return t_ex - dt_con - dt_sim, t_ex


def midrun_cum_slow(
    t_lab: np.ndarray,
    tc3: np.ndarray,
    t_ex: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Cumulative mid-run slowdown duration and onset times (lab s).

    Args:    t_lab, tc3  stateOS; t_ex  first execute after landing
    Returns: (cum_f, onset_lab_times)
    """
    is_slow = (np.asarray(tc3, dtype=int) == SLOWDOWN_STATE).astype(float)
    d_slow = np.zeros_like(t_lab, dtype=float)
    d_slow[1:] = is_slow[:-1] * np.diff(t_lab)
    d_slow[t_lab <= t_ex] = 0.0
    cum = np.cumsum(d_slow)
    prev = np.concatenate([[int(tc3[0])], np.asarray(tc3[:-1], dtype=int)])
    ons = t_lab[
        (np.asarray(tc3, dtype=int) == SLOWDOWN_STATE)
        & (prev != SLOWDOWN_STATE)
        & (t_lab > t_ex)
    ]
    return cum, np.asarray(ons, dtype=float)


def slowdown_lab_onsets(mat_name: str) -> list[float]:
    """
    Lab times at the start of each mid-run typeConv3==2 episode.

    Args:    mat_name
    Returns: onset times (lab s), empty if no stateOS
    """
    pair = load_type_conv3(mat_name)
    if pair is None:
        return []
    try:
        _k, t_ex = lab_clock_offset_s(mat_name)
    except RuntimeError:
        return []
    t_lab, tc3 = pair
    _cum, ons = midrun_cum_slow(t_lab, tc3, t_ex)
    return [float(x) for x in ons]


def map_os_to_lab(t_os_proto: np.ndarray, mat_name: str) -> LabTimeMap:
    """
    Map OpenSees prototype times onto lab Time.

    Args:    t_os_proto  recorder t (s, prototype); mat_name
    Returns: LabTimeMap
    """
    t_os = np.asarray(t_os_proto, dtype=float)
    k, t_ex = lab_clock_offset_s(mat_name)
    pair = load_type_conv3(mat_name)
    if pair is None:
        raise RuntimeError(f"no stateOS for {mat_name}")
    t_st, tc3 = pair
    cum, _ons = midrun_cum_slow(t_st, tc3, t_ex)
    target = k + t_os / TIME_SCALE_FROUDE
    g = np.asarray(t_st - cum, dtype=float)
    for i in range(1, g.size):
        if g[i] < g[i - 1]:
            g[i] = g[i - 1]
    idx = np.clip(np.searchsorted(g, target, side="left"), 0, t_st.size - 1)
    f_at = np.maximum(cum[idx], 0.0)
    t_lab = target + f_at
    return LabTimeMap(
        t_lab=t_lab,
        k=float(k),
        t_ex=float(t_ex),
        f_end=float(f_at[-1]) if f_at.size else 0.0,
        f_at=f_at,
    )


def os_times_to_lab(t_os_proto: np.ndarray | float, mat_name: str) -> np.ndarray:
    """
    Map one or more OS times to lab (same rule as ``map_os_to_lab``).

    Args:    t_os_proto  scalar or array (prototype s); mat_name
    Returns: lab s (array, same shape)
    """
    arr = np.atleast_1d(np.asarray(t_os_proto, dtype=float))
    mapped = map_os_to_lab(arr.ravel(), mat_name).t_lab.reshape(arr.shape)
    return mapped


def os_window_to_lab(
    t0_os: float,
    t1_os: float,
    mat_name: str,
) -> tuple[float, float]:
    """
    Map an OpenSees time window [t0, t1] to lab Time.

    Args:    t0_os, t1_os  prototype s; mat_name
    Returns: (t0_lab, t1_lab)
    """
    t_lab = os_times_to_lab(np.asarray([t0_os, t1_os], dtype=float), mat_name)
    return float(t_lab[0]), float(t_lab[1])


def resample_on_lab(
    t_lab: np.ndarray,
    y: np.ndarray,
    *,
    t0: float | None = None,
    t1: float | None = None,
    dt: float | None = None,
) -> tuple[np.ndarray, np.ndarray, float] | None:
    """
    Uniform lab-time samples for Welch / FFT (after k+f map).

    Args:    t_lab, y  mapped series; t0,t1  window (lab s); dt  step (default median)
    Returns: (t_u, y_u, dt) or None if too short
    """
    t = np.asarray(t_lab, dtype=float)
    yy = np.asarray(y, dtype=float)
    if t.size < 8 or yy.size != t.size:
        return None
    lo = float(t[0]) if t0 is None else float(t0)
    hi = float(t[-1]) if t1 is None else float(t1)
    m = (t >= lo) & (t <= hi) & np.isfinite(yy)
    if int(np.count_nonzero(m)) < 64:
        return None
    tt = t[m]
    vv = yy[m]
    # Enforce increasing t for interp (slowdowns can flat-line mapped t).
    order = np.argsort(tt, kind="mergesort")
    tt = tt[order]
    vv = vv[order]
    uniq = np.concatenate([[True], np.diff(tt) > 0.0])
    tt = tt[uniq]
    vv = vv[uniq]
    if tt.size < 64:
        return None
    dt_u = float(dt) if dt is not None else float(np.median(np.diff(tt)))
    if not np.isfinite(dt_u) or dt_u <= 0.0:
        return None
    t_u = np.arange(float(tt[0]), float(tt[-1]) + 0.5 * dt_u, dt_u)
    if t_u.size < 64:
        return None
    y_u = np.interp(t_u, tt, vv)
    y_u = y_u - float(np.mean(y_u))
    return t_u, y_u, dt_u


def add_dual_time_xaxis_lab(ax, *, top: bool = True) -> None:
    """
    Primary x = lab / model s; optional top = t√λ prototype s.

    Args:    ax; top  add secondary proto-time axis
    Returns: none
    """
    if not top:
        return
    sec = ax.secondary_xaxis(
        "top",
        functions=(
            lambda t_lab: t_lab * TIME_SCALE_FROUDE,
            lambda t_proto: t_proto / TIME_SCALE_FROUDE,
        ),
    )
    sec.set_xlabel(LABEL_T_PROTO_TOP)
