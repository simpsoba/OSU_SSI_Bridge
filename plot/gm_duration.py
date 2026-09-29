#!/usr/bin/env python3
"""
Goals
-----
Significant duration (D5–95) from a PEER velocity VT2 via Arias intensity.

  Classical Arias uses acceleration. We differentiate the VT2 (cm/s → m/s)
  with a central difference, then:

    I_A(t) = (π / 2g) ∫ a² dt
    t_5, t_95 = times when I_A / I_A∞ reaches 5% and 95%
    D5–95 = t_95 − t_5

Units: returned times are in the record's own clock (s). On the OpenSees /
prototype analysis clock that is ``t5 + gmStartTime`` … ``t95 + gmStartTime``
(Path ``-startTime``). Use ``d595_proto_window(gm_start_s=…)`` for zoom panels.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from paths import HERE

REPO = HERE.parent
DEFAULT_VT2 = (
    REPO / "ground-motion" / "Tohoku2011-FKSH" / "FKSH19.NS1.VT2"
)
G_ACCEL = 9.81  # m/s²


@dataclass(frozen=True)
class SignificantDuration:
    """Arias-based D5–95 on the ground-motion clock."""

    t5_s: float
    t95_s: float
    d5_95_s: float
    ia_total: float
    dt_s: float
    npts: int
    path: Path

    @property
    def duration_s(self) -> float:
        return (self.npts - 1) * self.dt_s


def load_peer_vt2_velocity_mps(path: Path) -> tuple[np.ndarray, float]:
    """
    Read a PEER VT2 velocity series.

    Args:    path  PEER NGA VT2 (header + values in cm/s)
    Returns: (v_mps, dt_s)
    """
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 5:
        raise ValueError(f"short PEER file: {path}")
    header = lines[3]
    match = re.search(
        r"NPTS\s*=\s*(\d+).*DT\s*=\s*([0-9.eE+-]+)",
        header,
        flags=re.IGNORECASE,
    )
    if not match:
        raise ValueError(f"no NPTS/DT in {path}: {header!r}")
    npts = int(match.group(1))
    dt_s = float(match.group(2))
    values_cm_s: list[float] = []
    for line in lines[4:]:
        values_cm_s.extend(float(token) for token in line.split())
    if len(values_cm_s) < npts:
        raise ValueError(
            f"{path}: expected {npts} samples, got {len(values_cm_s)}"
        )
    velocity_mps = np.asarray(values_cm_s[:npts], dtype=float) / 100.0
    return velocity_mps, dt_s


def arias_significant_duration(
    path: Path | None = None,
) -> SignificantDuration:
    """
    D5–95 from VT2 velocity via differentiated acceleration.

    Args:    path  default FKSH19.NS1.VT2 in the repo
    Returns: SignificantDuration on the GM clock (s)
    """
    vt2 = path or DEFAULT_VT2
    velocity_mps, dt_s = load_peer_vt2_velocity_mps(vt2)
    accel_mps2 = np.gradient(velocity_mps, dt_s)
    # I_A = (π/2g) ∫ a² dt  (m/s)
    arias = np.cumsum(accel_mps2**2) * dt_s * (np.pi / (2.0 * G_ACCEL))
    arias_total = float(arias[-1])
    if arias_total <= 0.0:
        raise ValueError(f"zero Arias intensity for {vt2}")
    fraction = arias / arias_total
    time_s = np.arange(len(fraction), dtype=float) * dt_s
    t5_s = float(np.interp(0.05, fraction, time_s))
    t95_s = float(np.interp(0.95, fraction, time_s))
    return SignificantDuration(
        t5_s=t5_s,
        t95_s=t95_s,
        d5_95_s=t95_s - t5_s,
        ia_total=arias_total,
        dt_s=dt_s,
        npts=len(velocity_mps),
        path=vt2,
    )


def gm_start_time_s(
    source: Path | Mapping[str, str] | None = None,
) -> float:
    """
    OpenSees Path ``-startTime`` (prototype s) from dump meta.

    Args:    source  dump folder, ``window_meta.txt[.0]``, or meta dict
    Returns: gmStartTime ≥ 0, or 0.0 if unset / missing
    """
    if source is None:
        return 0.0
    if isinstance(source, Mapping):
        raw = source.get("gmStartTime") or source.get("gmStart") or ""
        try:
            return max(0.0, float(str(raw).split()[0]))
        except (TypeError, ValueError, IndexError):
            return 0.0
    path = Path(source)
    if path.is_dir():
        for name in ("window_meta.txt", "window_meta.txt.0"):
            cand = path / name
            if cand.is_file():
                return gm_start_time_s(cand)
        return 0.0
    if not path.is_file():
        return 0.0
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        key, _, rest = s.partition(" ")
        if key in ("gmStartTime", "gmStart"):
            try:
                return max(0.0, float(rest.split()[0]))
            except (ValueError, IndexError):
                return 0.0
    return 0.0


def d595_proto_window(
    gm_start_s: float = 0.0,
    path: Path | None = None,
) -> tuple[float, float] | None:
    """
    D5–95 on the OpenSees / prototype analysis clock.

    Args:    gm_start_s  Path ``-startTime`` (s); path  optional VT2 override
    Returns: (t5 + gmStart, t95 + gmStart), or None if the VT2 is missing
    """
    try:
        dur = arias_significant_duration(path)
    except (OSError, ValueError):
        return None
    t0 = max(0.0, float(gm_start_s))
    return float(dur.t5_s) + t0, float(dur.t95_s) + t0


def d595_lab_window(
    gm_start_s: float = 0.0,
    *,
    time_scale_froude: float,
    path: Path | None = None,
) -> tuple[float, float] | None:
    """
    D5–95 on the Simulink lab / model clock (prototype ÷ √λ).

    Args:    gm_start_s  Path ``-startTime`` on the prototype clock (s)
             time_scale_froude  √λ (proto / model)
             path  optional VT2 override
    Returns: (t5_lab, t95_lab), or None
    """
    proto = d595_proto_window(gm_start_s, path=path)
    if proto is None:
        return None
    scale = float(time_scale_froude)
    if scale <= 0.0:
        raise ValueError(f"time_scale_froude must be > 0 (got {scale})")
    return proto[0] / scale, proto[1] / scale


if __name__ == "__main__":
    duration = arias_significant_duration()
    print(
        f"{duration.path.name}: t5={duration.t5_s:.3f}s  "
        f"t95={duration.t95_s:.3f}s  D5-95={duration.d5_95_s:.3f}s  "
        f"IA={duration.ia_total:.4f} m/s"
    )
