#!/usr/bin/env python3
"""
Goals
-----
Map OpenSees recorder time onto the Simulink lab clock for RTHS plots.

Handshake map (preferred)
  ctrlDisp recorder stamps completed steps as \(t_{i+1}\) with value
  \(u_{i+\alpha_f}\). Pair from the start 1:1 with lab **atTarget**
  (end of that \(\Delta t_{\mathrm{sim}}\) window — when force is valid).
  Index windows by enter-interpolate (target arrival); map uses atTarget.

  t_lab(t_OS) from those pairs already includes startup and waits.
  Define  f = t_lab − t_OS/√λ  (no separate k).

Legacy: typeConv3==2-only wait proxy (under-counts long extrapolate).

Lab Time is ground truth; never place lab events with naive t_lab·√λ.
See plot/lab/MKR_OPENFRESCO_TIMING.md §4.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from lab_paths import MAT_EXTRACT_DIR, M_TO_MM, TIME_SCALE_FROUDE

N_COUNT = 10  # Δt_sim = N_COUNT · Δt_con
SLOWDOWN_STATE = 2  # typeConv3 wait flag (OpenFresco)
INTERP_STATE = 0
EXTRAP_STATE = 1

# Lab / wall clock (Simulink Time). OpenSees integrator = $t_{\mathrm{int}}$.
LABEL_T_LAB = r"$t$ (s) lab / model scale"
LABEL_T_PROTO_TOP = r"$t\sqrt{\lambda}$ (s) prototype scale"
LABEL_T_INT = r"$t_{\mathrm{int}}$ (s) OpenSees / prototype scale"
LABEL_T_INT_MODEL_TOP = r"$t_{\mathrm{int}}/\sqrt{\lambda}$ (s) model scale"
LABEL_T_INT_SINCE_GM = r"$t_{\mathrm{int}}-\mathrm{gmStart}$ (s) prototype scale"


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


@dataclass(frozen=True)
class HandshakeMap:
    """
    ctrlDisp recorder \(t_{i+1}\) ↔ lab atTarget pairs (from start).

    ``t_lab`` is atTarget (window end). ``t_arrive`` is enter-interpolate
    (value appears). ``u_mm`` is ctrlDisp (proto mm).
    ``f_at = t_lab − t_OS/√λ`` (no separate k).
    """

    t_os: np.ndarray
    t_lab: np.ndarray
    t_arrive: np.ndarray
    u_mm: np.ndarray
    f_at: np.ndarray
    i_peak: int


def load_state_os(
    mat_name: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """
    stateOS time, typeConv3, count (typeConv1).

    Args:    mat_name
    Returns: (t_lab, typeConv3, count) or None
    """
    npz_path = MAT_EXTRACT_DIR / f"{Path(mat_name).stem}.npz"
    if not npz_path.is_file():
        return None
    z = np.load(npz_path, allow_pickle=True)
    if "stateOS_time" not in z.files or "stateOS_data" not in z.files:
        return None
    t = np.asarray(z["stateOS_time"], dtype=float)
    d = np.asarray(z["stateOS_data"], dtype=float)
    if d.ndim != 2 or d.shape[1] < 2 or t.size < 8:
        return None
    tc3 = np.rint(d[:, 0]).astype(int)
    count = np.rint(d[:, 1]).astype(int)
    return t, tc3, count


def enter_interpolate_indices(tc3: np.ndarray) -> np.ndarray:
    """
    Sample indices where typeConv3 enters interpolate (0) from nonzero.

    Args:    tc3  state enum
    Returns: indices (int array)
    """
    prev = np.concatenate([[int(tc3[0])], np.asarray(tc3[:-1], dtype=int)])
    return np.flatnonzero((prev != INTERP_STATE) & (np.asarray(tc3) == INTERP_STATE))


def at_target_indices(tc3: np.ndarray, arrive: np.ndarray) -> np.ndarray:
    """
    For each arrive index, lab atTarget = last sample before next extrapolate.

    If the target lands on the last substep, atTarget == arrive.

    Args:    tc3; arrive  enter-interpolate indices
    Returns: atTarget indices (same length as arrive)
    """
    tc3 = np.asarray(tc3, dtype=int)
    prev = np.concatenate([[int(tc3[0])], tc3[:-1]])
    # starts of extrapolate windows
    extrap_start = np.flatnonzero(
        (prev != EXTRAP_STATE) & (tc3 == EXTRAP_STATE)
    )
    at = np.empty(arrive.size, dtype=int)
    for i, j in enumerate(np.asarray(arrive, dtype=int)):
        # first extrap start strictly after arrive
        p = int(np.searchsorted(extrap_start, j + 1, side="left"))
        if p < extrap_start.size:
            at[i] = int(extrap_start[p]) - 1
        else:
            at[i] = int(tc3.size) - 1
        if at[i] < j:
            at[i] = j
    return at


def build_handshake_map(
    t_os: np.ndarray,
    u_os_mm: np.ndarray,
    t_state: np.ndarray,
    tc3: np.ndarray,
    *,
    u_corr_min: float = 0.99,
) -> HandshakeMap:
    """
    Pair ctrlDisp from the start with lab atTarget (1:1 handshakes).

    Args:    t_os, u_os_mm  ctrlDisp; t_state, tc3  stateOS
    Returns: HandshakeMap
    """
    t_os = np.asarray(t_os, dtype=float)
    u_os = np.asarray(u_os_mm, dtype=float)
    t_state = np.asarray(t_state, dtype=float)
    tc3 = np.asarray(tc3, dtype=int)
    if t_os.size < 8 or t_state.size < 8:
        raise RuntimeError("handshake map: ctrlDisp/stateOS too short")

    arrive = enter_interpolate_indices(tc3)
    if arrive.size < 8:
        raise RuntimeError("handshake map: too few enter-interpolate events")
    at_idx = at_target_indices(tc3, arrive)
    t_arrive = t_state[arrive]
    t_at = t_state[at_idx]

    n = min(int(t_os.size), int(arrive.size))
    if abs(int(t_os.size) - int(arrive.size)) > 2:
        raise RuntimeError(
            f"handshake map: ctrl={t_os.size} vs arrive={arrive.size}"
        )
    t_os_p = t_os[:n].copy()
    u_p = u_os[:n].copy()
    t_lab_p = t_at[:n].copy()
    t_arr_p = t_arrive[:n].copy()
    for i in range(1, t_lab_p.size):
        if t_lab_p[i] < t_lab_p[i - 1]:
            t_lab_p[i] = t_lab_p[i - 1]

    i_peak = int(np.nanargmin(u_p))
    f_at = t_lab_p - t_os_p / TIME_SCALE_FROUDE
    _ = u_corr_min  # reserved for tar QA in handshake_map_for_dump
    return HandshakeMap(
        t_os=t_os_p,
        t_lab=t_lab_p,
        t_arrive=t_arr_p,
        u_mm=u_p,
        f_at=f_at,
        i_peak=i_peak,
    )


def _handshake_monotone_pairs(
    hs: HandshakeMap,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Monotone (t_os, t_lab) for interpolation either way.

    Args:    hs
    Returns: (t_os, t_lab)
    """
    t_os = np.array(hs.t_os, dtype=float, copy=True)
    t_lab = np.array(hs.t_lab, dtype=float, copy=True)
    for i in range(1, t_os.size):
        if t_os[i] < t_os[i - 1]:
            t_os[i] = t_os[i - 1]
        if t_lab[i] < t_lab[i - 1]:
            t_lab[i] = t_lab[i - 1]
    return t_os, t_lab


def map_os_times_handshake(
    t_os_proto: np.ndarray | float,
    hs: HandshakeMap,
) -> np.ndarray:
    """
    Interpolate handshake pairs onto query OS times.

    Args:    t_os_proto  proto s; hs  HandshakeMap
    Returns: lab s (same shape)
    """
    arr = np.atleast_1d(np.asarray(t_os_proto, dtype=float))
    t_os, t_lab = _handshake_monotone_pairs(hs)
    return np.interp(arr.ravel(), t_os, t_lab).reshape(arr.shape)


def lab_times_to_os_handshake(
    t_lab_s: np.ndarray | float,
    hs: HandshakeMap,
) -> np.ndarray:
    """
    Inverse map: lab Time → last completed OpenSees recorder t (proto s).

    Piecewise-constant on atTarget knots (not linear in a wait gap).

    Args:    t_lab_s  lab s; hs  HandshakeMap
    Returns: t_os proto s (same shape)
    """
    arr = np.atleast_1d(np.asarray(t_lab_s, dtype=float))
    t_os, t_lab = _handshake_monotone_pairs(hs)
    uniq = np.concatenate([[True], np.diff(t_lab) > 0.0])
    t_lab_u = t_lab[uniq]
    t_os_u = t_os[uniq]
    idx = np.searchsorted(t_lab_u, arr.ravel(), side="right") - 1
    idx = np.clip(idx, 0, t_os_u.size - 1)
    return t_os_u[idx].reshape(arr.shape)


def load_ctrl_disp_proto_mm(
    dump_dir: Path,
) -> tuple[np.ndarray, np.ndarray] | None:
    """
    OpenFresco ctrlDisp recorder (every analysis step).

    Args:    dump_dir  eqOutDir
    Returns: (t_os_proto_s, u_proto_mm) or None
    """
    dump_dir = Path(dump_dir)
    cands = sorted(dump_dir.glob("Elmt*_ctrlDsp.out*"))
    if not cands:
        return None
    # prefer unsuffixed / .out.0
    cands.sort(key=lambda p: (p.suffix != ".out", p.name))
    try:
        d = np.loadtxt(cands[0], ndmin=2)
    except Exception:
        return None
    if d.ndim != 2 or d.shape[1] < 2 or d.shape[0] < 8:
        return None
    t = np.asarray(d[:, 0], dtype=float)
    u = (np.asarray(d[:, 1], dtype=float) - float(d[0, 1])) * M_TO_MM
    return t, u


def handshake_map_for_dump(
    dump_dir: Path,
    mat_name: str,
    t_tar: np.ndarray | None = None,
    u_tar_mm: np.ndarray | None = None,
    *,
    u_corr_min: float = 0.99,
) -> HandshakeMap:
    """
    Build atTarget handshake map for a dump + mat.

    Args:    dump_dir; mat_name; optional tar for u-correlation QA
    Returns: HandshakeMap
    """
    ctrl = load_ctrl_disp_proto_mm(dump_dir)
    if ctrl is None:
        raise RuntimeError(f"no Elmt*_ctrlDsp.out under {dump_dir}")
    st = load_state_os(mat_name)
    if st is None:
        raise RuntimeError(f"no stateOS for {mat_name}")
    t_os, u_os = ctrl
    t_state, tc3, _count = st
    hs = build_handshake_map(t_os, u_os, t_state, tc3)
    if t_tar is not None and u_tar_mm is not None:
        u_ref = np.interp(hs.t_arrive, np.asarray(t_tar, dtype=float), np.asarray(u_tar_mm, dtype=float))
        a = hs.u_mm - float(np.mean(hs.u_mm))
        b = u_ref - float(np.mean(u_ref))
        den = float(np.linalg.norm(a) * np.linalg.norm(b))
        corr = float(np.dot(a, b) / den) if den > 0.0 else 0.0
        if corr < float(u_corr_min):
            raise RuntimeError(
                f"handshake map: ctrl vs tar@arrive corr={corr:.4f} < {u_corr_min}"
            )
    return hs


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


def _dump_dir_for_mat(mat_name: str) -> Path | None:
    """
    LOCAL dump folder for a MatFile (TestMatrix pairing).

    Args:    mat_name
    Returns: Path or None
    """
    from lab_paths import (
        LOCAL_OPENSEES_DATA,
        mat_to_run_from_csv,
        resolve_opensees_data,
    )

    stem = Path(mat_name).name
    mats = mat_to_run_from_csv()
    dump = mats.get(stem)
    if dump is None:
        stem_key = Path(stem).stem
        for k, v in mats.items():
            if Path(k).stem == stem_key:
                dump = v
                break
    if not dump:
        return None
    root = resolve_opensees_data() or LOCAL_OPENSEES_DATA
    path = Path(root) / dump
    return path if path.is_dir() else None


def _map_os_to_lab_typeconv3(t_os: np.ndarray, mat_name: str) -> LabTimeMap:
    """
    Legacy map: k + t_OS/√λ + cum(typeConv3==2). Under-counts waits.

    Args:    t_os  proto s; mat_name
    Returns: LabTimeMap
    """
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


def map_os_to_lab(
    t_os_proto: np.ndarray,
    mat_name: str,
    dump_dir: Path | str | None = None,
) -> LabTimeMap:
    """
    Map OpenSees prototype times onto lab Time (atTarget handshake).

    Prefers ctrlDisp↔stateOS atTarget pairs. Falls back to typeConv3==2
    wait proxy if the dump/ctrlDisp is missing.

    Args:    t_os_proto  recorder t (s, prototype); mat_name; dump_dir
    Returns: LabTimeMap  (k=0 for handshake; f = t_lab − t_OS/√λ)
    """
    t_os = np.asarray(t_os_proto, dtype=float)
    dump = Path(dump_dir) if dump_dir is not None else _dump_dir_for_mat(mat_name)
    if dump is not None and dump.is_dir():
        try:
            hs = handshake_map_for_dump(dump, mat_name)
            t_lab = np.asarray(map_os_times_handshake(t_os, hs), dtype=float)
            f_at = t_lab - t_os / TIME_SCALE_FROUDE
            try:
                _k, t_ex = lab_clock_offset_s(mat_name)
            except RuntimeError:
                t_ex = 0.0
            return LabTimeMap(
                t_lab=t_lab,
                k=0.0,
                t_ex=float(t_ex),
                f_end=float(f_at[-1]) if f_at.size else 0.0,
                f_at=f_at,
            )
        except RuntimeError:
            pass
    return _map_os_to_lab_typeconv3(t_os, mat_name)


def os_times_to_lab(
    t_os_proto: np.ndarray | float,
    mat_name: str,
    dump_dir: Path | str | None = None,
) -> np.ndarray:
    """
    Map one or more OS times to lab (same rule as ``map_os_to_lab``).

    Args:    t_os_proto  scalar or array (prototype s); mat_name; dump_dir
    Returns: lab s (array, same shape)
    """
    arr = np.atleast_1d(np.asarray(t_os_proto, dtype=float))
    mapped = map_os_to_lab(arr.ravel(), mat_name, dump_dir=dump_dir).t_lab
    return mapped.reshape(arr.shape)


def os_window_to_lab(
    t0_os: float,
    t1_os: float,
    mat_name: str,
    dump_dir: Path | str | None = None,
) -> tuple[float, float]:
    """
    Map an OpenSees time window [t0, t1] to lab Time.

    Args:    t0_os, t1_os  prototype s; mat_name; dump_dir
    Returns: (t0_lab, t1_lab)
    """
    t_lab = os_times_to_lab(
        np.asarray([t0_os, t1_os], dtype=float), mat_name, dump_dir=dump_dir
    )
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
    Uniform lab-time samples for FFT (after k+f map). Mean is kept.

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
    return t_u, y_u, dt_u


def add_dual_time_xaxis_lab(ax, *, top: bool = True) -> None:
    """
    Primary x = lab $t$ (model s); optional top = $t\sqrt{\lambda}$ proto s.

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


def add_dual_time_xaxis_int(ax, *, top: bool = True) -> None:
    """
    Primary x = $t_{\mathrm{int}}$ (proto s); optional top = model via /√λ.

    Args:    ax; top  add secondary model-time axis
    Returns: none
    """
    if not top:
        return
    sec = ax.secondary_xaxis(
        "top",
        functions=(
            lambda t_int: t_int / TIME_SCALE_FROUDE,
            lambda t_model: t_model * TIME_SCALE_FROUDE,
        ),
    )
    sec.set_xlabel(LABEL_T_INT_MODEL_TOP)
