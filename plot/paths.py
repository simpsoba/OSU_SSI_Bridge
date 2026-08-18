# Shared output layout for plot/*.py
#
#   plot/out/profile{N}/
#     eq/{serial|parallel}/{Shin|ASDEA}/{quad|SSPquad}/{pierEleType}/
#     eq/compare/{pierEleType}/
#     elevation/{Shin|ASDEA}/elevation.png
#     pile_springs/{pult,tult,y50z50}.png
#     soil_profile/{overview,pdmy02}.png
#     fibers/fiber_*.png

from __future__ import annotations

from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT_ROOT = HERE / "out"


def profile_root(soil_profile: int | str) -> Path:
    return OUT_ROOT / f"profile{int(soil_profile)}"


def elevation_dir(soil_profile: int | str, soil_boundary: str) -> Path:
    return profile_root(soil_profile) / "elevation" / str(soil_boundary)


def elevation_png(soil_profile: int | str, soil_boundary: str) -> Path:
    return elevation_dir(soil_profile, soil_boundary) / "elevation.png"


def pile_springs_dir(soil_profile: int | str) -> Path:
    return profile_root(soil_profile) / "pile_springs"


def soil_profile_dir(soil_profile: int | str) -> Path:
    return profile_root(soil_profile) / "soil_profile"


def fibers_dir(soil_profile: int | str) -> Path:
    return profile_root(soil_profile) / "fibers"


def eq_dir(
    soil_profile: int | str,
    soil_boundary: str,
    soil_ele_type: str = "quad",
    pier_ele_type: str = "lumpedPlasticity",
    run_kind: str = "serial",
) -> Path:
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
    return profile_root(soil_profile) / "eq" / "compare" / pier_ele_type


def gravity_dir(
    soil_profile: int | str,
    soil_boundary: str,
    soil_ele_type: str = "quad",
    pier_ele_type: str = "lumpedPlasticity",
) -> Path:
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
    return (
        profile_root(soil_profile)
        / "partition"
        / str(soil_boundary)
        / str(soil_ele_type)
        / pier_ele_type
        / f"np{int(np)}"
    )
