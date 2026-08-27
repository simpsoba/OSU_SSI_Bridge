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

Units: returned times are in the record's own clock (s), matching OpenSees
t_num when gmStartTime = 0 on a prototype-scale schedule.
"""

from __future__ import annotations

import re
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


if __name__ == "__main__":
    duration = arias_significant_duration()
    print(
        f"{duration.path.name}: t5={duration.t5_s:.3f}s  "
        f"t95={duration.t95_s:.3f}s  D5-95={duration.d5_95_s:.3f}s  "
        f"IA={duration.ia_total:.4f} m/s"
    )
