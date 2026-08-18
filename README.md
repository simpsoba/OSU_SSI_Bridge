# OSU SSI Bridge

2D pier-bridge SSI model, built incrementally.

- Pier geometry: Shin et al. (2007); Mackie et al. (2008)
- Lab cylinder / hydro-RTHS: Neumann (2021); Neumann et al. (2023)
- Plane: transverse
- Above grade: base rotational spring + stiff beam
- Soil: layered continuum (`soil/Profiles.md`; L3 sand or clay)

See **[NOTES.md](NOTES.md)**. Cite keys: **[reference/references.bib](reference/references.bib)** (PDFs stay local, not in git).

Clone (private; org members have write):

```bash
git clone git@github.com:simpsoba/OSU_SSI_Bridge.git
```

Run from the repo root (`file join` / `pathlib`, so Windows and Linux checkouts are the same). Checkout is LF (`.gitattributes`); do not convert Tcl to CRLF. `OpenSees` / `OpenSees.exe` on PATH. Plotters use `python3` or `python`. Do not commit `plot/out/`. Lab workflow: [simpsoba/onboarding](https://github.com/simpsoba/onboarding).

Edit knobs in `Parameters.tcl` (grep `EDIT` / `USER INPUT`; IDs under TAGS CONVENTION). Driver: `Run.tcl` (`runEQ` 0 = gravity+eigen, 1 = gravity+EQ). Analysis settings (`numberer`, `system`, `analyze`) are in `Run.tcl`; staging and recorders live in `analysis/`.
Figures (build + dump + Python plotters): `OpenSees PlotModel.tcl`.
Mode shapes after gravity+eigen: `plot/out/profile{N}/elevation/{BC}/modes/{soilEle}/{pierEleType}/mode_01.png` …
Gravity deformed: `plot/out/profile{N}/elevation/{BC}/gravity/{soilEle}/{pierEleType}/gravity_deformed.png`
EQ dumps: `plot/out/profile{N}/eq/{serial|parallel}/{BC}/{soilEle}/{pierEleType}/`
- `pierEleType`: `elasticBeamColumn` | `forceBeamColumn` | `lumpedPlasticity`
- `pierCrackedFactor`: elastic pier \(I\) × factor; also stretches lumped \(L_s\) (default 0.5; Fiber \(E_c,E_s\) unchanged)
- `pileEleType`: `elasticBeamColumn` | `dispBeamColumn`
- `soilProfile`: `1` | `2` | `3` | `4`
- `soilBoundary`: `Shin` | `ASDEA`
- `soilEleType`: `quad` | `SSPquad` (SSPquad mass from material ρ; no ele ρ)
- `soilConstitutive`: `inelastic` (PIMY/PDMY02) | `elastic` (`ElasticIsotropic3D` from skeleton \(G_r,B_r\); FSP on sands)
- `pileSpring`: `inelastic` | `elastic` | `none` (p-y / t-z / q-z together)

Plots (`plot/out/profile{N}/`) after `PlotModel.tcl`:
```bash
OpenSees PlotModel.tcl
# or JSON only:  set ::plotSkipPython 1; source PlotModel.tcl
# then optionally re-run one plotter:
python3 plot/PlotModelSketch.py      # → elevation/{Shin|ASDEA}/elevation.png
python3 plot/PlotFiberSections.py    # → fibers/fiber_*.png
python3 plot/PlotPileSprings.py      # → pile_springs/{pult,tult,y50z50}.png
python3 plot/PlotSoilProfile.py      # → soil_profile/{overview,pdmy02}.png
```
