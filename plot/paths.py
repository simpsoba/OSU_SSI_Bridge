#!/usr/bin/env python3
"""
Goals
-----
Shared output layout for local (non-lab) plot/*.py dumps under plot/out/.

  plot/out/profile{N}/
    eq/{serial|parallel}/{Shin|ASDEA}/{quad|SSPquad}/{pierEleType}/
    eq/compare/{pierEleType}/
    elevation/{Shin|ASDEA}/elevation.png
    pile_springs/{pult,tult,y50z50}.png
    soil_profile/{overview,pdmy02}.png
    fibers/fiber_*.png

Lab campaign PNGs do *not* use this tree — see lab_paths.plots_root().
"""

from __future__ import annotations

from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT_ROOT = HERE / "out"


def profile_root(soil_profile: int | str) -> Path:
    """
    Args:    soil_profile  matrix soilProfile index
    Returns: plot/out/profile{N}/
    """
    return OUT_ROOT / f"profile{int(soil_profile)}"


def elevation_dir(soil_profile: int | str, soil_boundary: str) -> Path:
    """
    Args:    soil_profile, soil_boundary  (Shin | ASDEA)
    Returns: …/elevation/{boundary}/
    """
    return profile_root(soil_profile) / "elevation" / str(soil_boundary)


def elevation_png(soil_profile: int | str, soil_boundary: str) -> Path:
    """Args/Returns: path to elevation.png under elevation_dir()."""
    return elevation_dir(soil_profile, soil_boundary) / "elevation.png"


def pile_springs_dir(soil_profile: int | str) -> Path:
    """Args/Returns: …/pile_springs/."""
    return profile_root(soil_profile) / "pile_springs"


def soil_profile_dir(soil_profile: int | str) -> Path:
    """Args/Returns: …/soil_profile/."""
    return profile_root(soil_profile) / "soil_profile"


def fibers_dir(soil_profile: int | str) -> Path:
    """Args/Returns: …/fibers/."""
    return profile_root(soil_profile) / "fibers"


def eq_dir(
    soil_profile: int | str,
    soil_boundary: str,
    soil_ele_type: str = "quad",
    pier_ele_type: str = "lumpedPlasticity",
    run_kind: str = "serial",
) -> Path:
    """
    EQ recorder dump folder for a local plot/out run.

    Args:    soil_profile, soil_boundary, soil_ele_type, pier_ele_type,
             run_kind  ("serial" | "parallel")
    Returns: …/eq/{run_kind}/{boundary}/{ele}/{pier}/
    """
    return (
        profile_root(soil_profile)
        / "eq"
        / str(run_kind)
        / str(soil_boundary)
        / str(soil_ele_type)
        / pier_ele_type
    )


def eq_compare_dir(
    soil_profile: int | str,
    pier_ele_type: str = "lumpedPlasticity",
) -> Path:
    """Args/Returns: …/eq/compare/{pierEleType}/."""
    return profile_root(soil_profile) / "eq" / "compare" / pier_ele_type


def gravity_dir(
    soil_profile: int | str,
    soil_boundary: str,
    soil_ele_type: str = "quad",
    pier_ele_type: str = "lumpedPlasticity",
) -> Path:
    """Args/Returns: …/elevation/…/gravity/{ele}/{pier}/."""
    return (
        elevation_dir(soil_profile, soil_boundary)
        / "gravity"
        / str(soil_ele_type)
        / pier_ele_type
    )


def modes_dir(
    soil_profile: int | str,
    soil_boundary: str,
    pier_ele_type: str = "lumpedPlasticity",
    soil_ele_type: str = "quad",
) -> Path:
    """Args/Returns: …/elevation/…/modes/{ele}/{pier}/."""
    return (
        elevation_dir(soil_profile, soil_boundary)
        / "modes"
        / str(soil_ele_type)
        / pier_ele_type
    )


def modes_png(
    soil_profile: int | str,
    soil_boundary: str,
    pier_ele_type: str = "lumpedPlasticity",
    soil_ele_type: str = "quad",
) -> Path:
    """Args/Returns: mode_01.png under modes_dir()."""
    return modes_dir(
        soil_profile, soil_boundary, pier_ele_type, soil_ele_type
    ) / "mode_01.png"


def partition_dir(
    soil_profile: int | str,
    soil_boundary: str,
    soil_ele_type: str = "quad",
    pier_ele_type: str = "lumpedPlasticity",
    np: int = 2,
) -> Path:
    """
    Args:    np  MPI rank count (folder tag np{N})
    Returns: …/partition/{boundary}/{ele}/{pier}/np{N}/
    """
    return (
        profile_root(soil_profile)
        / "partition"
        / str(soil_boundary)
        / str(soil_ele_type)
        / pier_ele_type
        / f"np{int(np)}"
    )
