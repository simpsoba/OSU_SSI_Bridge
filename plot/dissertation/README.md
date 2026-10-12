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
| `fig_wed_force_kin.py` | `prelim_force_kinematics` (not yet in the chapter) | W05 interface force vs measured u, u-dot, u-ddot (20 Hz low-pass), with the inertia fit F = m u-ddot |
| `fig_wed_energy2.py` | `prelim_interface_energy` | energy across the interface; offline 1 vs 8 ranks; hysteresis of pier base, pile, p-y spring |

Helpers: `ctrl_error_wed.py` (run list, phases, loaders; prints NRMSE and delays by phase when run),
`wed_timing_stats.py` (new-target times, states, tracking cut), `host_time.py` (host latency per window).

Data: `LOCAL/mat_extract/*.npz` and `OSU_SSI_BRIDGE_DATA/Simulink/*.mat` (via `plot/lab_paths.py`),
`OSU_SSI_BRIDGE_DATA/<run>/` OpenSees recorders, and the offline rank-check runs in
`C:\Users\garaujor\OpenSees_Runs\OSU_SSI_Bridge_rankcheck\rankcheck\out_np{1,8}` (energy figure, panel b).
Definitions are in each script's docstring (tracking cut of the 1-rank run, phases, latency, k_el, PSD windows).

## Main campaign (Friday 2026-08-21; drafts, not yet in the chapter)

Earthquake followed by the solitary wave; baseline/moderate/large/extra-large models (2028/3216/4620/5484 DOFs);
CPU MKR-alpha (ProfileSPD, MUMPS, or cuDSS solves) vs GPU MKR-alpha (`CudaMKRAlpha`) + cuDSS. Timing comes from
stateOS (`LOCAL/mat_extract`), force from the Simulink load cell. `python plot/dissertation/make_figures.py --main`
writes them to `out/friday`.

| Script | Draft figure | Content |
|---|---|---|
| `fig_fri_host.py` | `fri_host_timing` | host latency per window and rate-transition state shares for every wet test, grouped by model |
| `fig_fri_capacity.py` | `fri_capacity` | mean latency and slowdown rate vs model size, by configuration |
| `fig_fri_lanes.py` | `fri_slowdown_lanes` | slowdown episodes on integration time (Arias band, solitary wave); delay growth behind real time |
| `fig_fri_ranks.py` | `fri_rank_sweep` | GPU rank sweeps of the moderate, large, and extra-large models, wet tests and lunch dry runs |
| `fig_fri_formulation.py` | `fri_formulation` | elastic / stiff / soft SSPQuad and Quad-PIMY soil on the baseline model: latency over time, slowdowns |
| `fig_fri_force.py` | `fri_interface_force` | hold-induced force transients vs the solitary-wave impact |
| `fig_fri_disp.py` | `fri_interface_displacement` | slowdowns in the displacement history: lab-clock lag (moderate, CPU vs GPU) and GPU 1 vs 8 ranks on integration time |
| `fig_fri_solitary.py` | `fri_solitary_mesh` | displacement across the four meshes (GPU); solitary-wave pulse and the free vibration it re-excites |

Helper: `fri_runs.py` (run list with configuration, dry runs, timing loader; prints a stats table when run).

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
