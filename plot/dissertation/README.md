# Dissertation figures (Ch. 6)

Entry point for the Ch. 6 figures of the dissertation that come from this repository.

```
python plot/dissertation/make_figures.py --list
python plot/dissertation/make_figures.py --copy-to <dissertation-repo>/chapters/06-case-study/figures
python plot/dissertation/make_figures.py --with-upstream --copy-to <dissertation-repo>/chapters/06-case-study/figures
```

## Preliminary campaign (scripts in this folder)

Wednesday 2026-08-19 runs W04 (1 rank), W06 (4), W02 (8), W05 (16), and W07 (24): CPU MKR-alpha
(`MKRAlphaExplicitMultiSOE -incrementalAccel`) with `ProfileSPD`, 2136-DOF model. The scripts only read data.

| Script | Ch. 6 figure (PDF name) | Content |
|---|---|---|
| `fig_wed_host2.py` | `prelim_host_timing` | time to first target; rate-transition state shares and slowdown episodes; host latency per window; interval between targets |
| `fig_wed_disp.py` | `prelim_interface_displacement` | interface displacement on real and integrator time; delay of the targets behind real time |
| `fig_wed_force.py` | `prelim_interface_force` | interface force on real time, with zooms |
| `fig_wed_psd.py` | `prelim_interface_psd` | PSDs of the interface force and displacement |
| `fig_wed_waves.py` | `prelim_wave_check` | hybrid W02 vs offline (no waves) pier top; force and pier-top motion under waves before the earthquake |
| `fig_wed_err_psd.py` | `prelim_control_error` | PSD of command - measured; measured vs command |
| `fig_wed_energy2.py` | `prelim_interface_energy` | energy across the interface; offline 1 vs 8 ranks; hysteresis of pier base, pile, p-y spring |

Helpers: `ctrl_error_wed.py` (run list, phases, loaders; prints NRMSE and delays by phase when run),
`wed_timing_stats.py` (new-target times, states, tracking cut), `host_time.py` (host latency per window).

Data: `LOCAL/mat_extract/*.npz` and `OSU_SSI_BRIDGE_DATA/Simulink/*.mat` (via `plot/lab_paths.py`),
`OSU_SSI_BRIDGE_DATA/<run>/` OpenSees recorders, and the offline rank-check runs in
`C:\Users\garaujor\OpenSees_Runs\OSU_SSI_Bridge_rankcheck\rankcheck\out_np{1,8}` (energy figure, panel b).
Definitions are in each script's docstring (tracking cut of the 1-rank run, phases, latency, k_el, PSD windows).

## Offline Newmark vs MKR-alpha (merged figure)

`fig_offline_combined.py` -> `offline_integration_combined` (Ch. 6 figure `fig:case-study-offline`): the two paper
figures below (displacement histories and hysteresis loops) in one figure with one mesh, legend, and caption. It
imports `plot/PlotLumpedHoldIntegratorPaper.py` and `plot/PlotLumpedHoldHysteresisPaper.py` for data, styling,
zoom windows, and NRMSE, so style changes there carry over.

## Other Ch. 6 figures (paper scripts in `plot/`, left in place)

These read JSON and recorder outputs relative to `plot/`, so they are not copied here; `--with-upstream` runs
them unchanged. Each output below was checked byte-identical to the PDF in the chapter (2026-10-09).

| Ch. 6 figure | Script | Output |
|---|---|---|
| `numerical_subassembly_schematic` | `plot/PlotModelSketchPaper.py` (default `--fill none`) | `plot/out/profile4/elevation/Shin/elevation_paper.pdf` |
| `soil_props_paper` | `plot/PlotSoilProfilePaper.py` | `plot/out/profile4/soil_profile/soil_props_paper.pdf` |
| `tohoku_FKSH19_NS1_gm_Sa` | `plot/PlotGroundMotionFigure.py` (spectrum cache from `PlotResponseSpectrum.py`) | `plot/out/spectrum/tohoku_FKSH19_NS1_gm_Sa.pdf` |
| `case_study_eigen_modes` | `plot/PlotEigenModesPaper.py` | `modal_analysis/figures/case_study_eigen_modes.pdf` |
| `lumped_hold_integrator_every`, `_rel` | `plot/PlotLumpedHoldIntegratorPaper.py` | `plot/out/eq_offline/compare/lumped_hold_integrator_every{,_rel}.pdf` |
| `lumped_hold_hysteresis_every` | `plot/PlotLumpedHoldHysteresisPaper.py` | `plot/out/eq_offline/compare/lumped_hold_hysteresis_every.pdf` |

Not from this repository: `bridge_schematic`, `physical_setup_schematic`, `static_equilibrium`, `three-loop-rths`,
and `mkr_openfresco_timing_strip` (Inkscape), and the CPU/GPU flowcharts (`gen_*.py`); all in
`chapters/06-case-study/figures/src` of the dissertation repository.
