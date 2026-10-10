"""Regenerate the Ch. 6 figures of the dissertation that come from OSU_SSI_Bridge.

  python plot/dissertation/make_figures.py                     # preliminary campaign -> plot/dissertation/out
  python plot/dissertation/make_figures.py --copy-to <figdir>  # ... and copy the PDFs under their Ch. 6 names
  python plot/dissertation/make_figures.py --with-upstream --copy-to <figdir>
                                                               # also rerun the upstream paper scripts in plot/
  python plot/dissertation/make_figures.py --list              # print the script -> figure map

Preliminary-campaign scripts live in this folder. The other Ch. 6 figures are made by the paper scripts in plot/,
which stay where they are (they read JSON and recorder outputs relative to plot/ and write to their usual
outputs); --with-upstream runs them unchanged and copies their PDFs. Read-only on the lab data.
Hand-drawn Ch. 6 figures (bridge_schematic, physical_setup_schematic, static_equilibrium, three-loop-rths, the
CPU/GPU flowcharts, mkr_openfresco_timing_strip) have their sources in chapters/06-case-study/figures/src.
"""
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PLOT = HERE.parent
REPO = PLOT.parent
OUT = HERE / "out"

# script in this folder -> Ch. 6 figure name
PRELIM = {
    "fig_offline_combined.py": "offline_integration_combined",
    "fig_wed_host2.py": "prelim_host_timing",
    "fig_wed_disp.py": "prelim_interface_displacement",
    "fig_wed_force.py": "prelim_interface_force",
    "fig_wed_psd.py": "prelim_interface_psd",
    "fig_wed_waves.py": "prelim_wave_check",
    "fig_wed_err_psd.py": "prelim_control_error",
    "fig_wed_energy2.py": "prelim_interface_energy",
}
# Ch. 6 figure name -> (script in plot/ with arguments, PDF it writes); checked byte-identical to the Ch. 6 PDFs
UPSTREAM = {
    "numerical_subassembly_schematic": (["PlotModelSketchPaper.py"],
                                        PLOT / "out/profile4/elevation/Shin/elevation_paper.pdf"),
    "soil_props_paper": (["PlotSoilProfilePaper.py"], PLOT / "out/profile4/soil_profile/soil_props_paper.pdf"),
    "tohoku_FKSH19_NS1_gm_Sa": (["PlotGroundMotionFigure.py"], PLOT / "out/spectrum/tohoku_FKSH19_NS1_gm_Sa.pdf"),
    "case_study_eigen_modes": (["PlotEigenModesPaper.py"], REPO / "modal_analysis/figures/case_study_eigen_modes.pdf"),
    "lumped_hold_integrator_every": (["PlotLumpedHoldIntegratorPaper.py"],
                                     PLOT / "out/eq_offline/compare/lumped_hold_integrator_every.pdf"),
    "lumped_hold_integrator_every_rel": (None, PLOT / "out/eq_offline/compare/lumped_hold_integrator_every_rel.pdf"),
    "lumped_hold_hysteresis_every": (["PlotLumpedHoldHysteresisPaper.py"],
                                     PLOT / "out/eq_offline/compare/lumped_hold_hysteresis_every.pdf"),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--copy-to", type=Path, help="figure folder of the dissertation chapter")
    ap.add_argument("--only", nargs="*", help="subset of Ch. 6 figure names")
    ap.add_argument("--with-upstream", action="store_true", help="also rerun the paper scripts in plot/")
    ap.add_argument("--list", action="store_true", help="print the script -> figure map and exit")
    args = ap.parse_args()
    if args.list:
        for s, n in PRELIM.items():
            print(f"{n:34s} plot/dissertation/{s}")
        for n, (cmd, pdf) in UPSTREAM.items():
            print(f"{n:34s} plot/{cmd[0] if cmd else '(written with the entry above)'} -> {pdf.relative_to(REPO)}")
        return
    want = lambda n: not args.only or n in args.only  # noqa: E731
    OUT.mkdir(exist_ok=True)
    for script, name in PRELIM.items():
        if not want(name):
            continue
        print(f"--- {script} -> {name}", flush=True)
        subprocess.run([sys.executable, "-I", str(HERE / script), str(OUT / name)], check=True, cwd=HERE)
        if args.copy_to:
            shutil.copy2(OUT / f"{name}.pdf", args.copy_to / f"{name}.pdf")
            print("copied", args.copy_to / f"{name}.pdf")
    if not args.with_upstream:
        return
    for name, (cmd, pdf) in UPSTREAM.items():
        if not want(name):
            continue
        if cmd:
            print(f"--- plot/{cmd[0]} -> {name}", flush=True)
            subprocess.run([sys.executable, str(PLOT / cmd[0]), *cmd[1:]], check=True, cwd=REPO)
        if args.copy_to:
            shutil.copy2(pdf, args.copy_to / f"{name}.pdf")
            print("copied", args.copy_to / f"{name}.pdf")


if __name__ == "__main__":
    main()
