# OSU SSI Bridge — Design Notes

Living notes for a **2D pier + soil SSI** model, built piece by piece.

---

## Intent

Single-pier (or pier-focused) bridge SSI model:

- Pier geometry from Shin et al. (2007), Mackie et al. (2008), and related
  PEER reports (`reference/`).
- Hybrid simulation compatible with Neumann (2021) and Neumann et al. (2023);
  similitude notes in Seki et al. (2026).
- Soil: layered continuum; profiles in `soil/Profiles.md`.
- Motion in the **transverse** plane.
- Above-grade pier: one base rotational spring with a stiff beam to the top.
- Transverse deck: stiff PT39 frame (`structure/BuildDeckNodes.tcl`).

Prefer a small driver (`Run.tcl`) + sourced modules over one large script.

---

## Hard constraints (agreed)

| Constraint | Detail |
|---|---|
| Dimension | **2D** (`BasicBuilder -ndm 2 -ndf 3`) |
| Units | **N, m, s** |
| EQ direction | **Transverse** |
| Above-ground pier | **one base rotational spring + stiff beam to top** |
| Physical pier | Lab cylinder: **D_cyl = 20.0 in**, **H_cyl = 3.02 m** (Neumann `Hc_m`) |
| Prototype pier D | **D_pier = 4.0 ft** (Shin et al. 2007); **cylinderSF = D_pier/D_cyl**; **H_pier = H_cyl·cylinderSF** |
| Build style | Incremental; each piece understandable alone |
| Legibility | Comment node/element maps; named tags; gravity from `structNodeTags`; switches in `Run.tcl`; knobs in `Parameters.tcl` |

### Above-ground pier (lumped plasticity)

```
  mudline / base                              top of pier
  node 1 ──spring── node 2 ════ stiff beam ════ node 5
           (rot)              (η · EI)
```

- Base spring: `zeroLengthSection`, rotation (dir 3)
- Stiff beam: `elasticBeamColumn` with `η · I` (nodes 2→5)
- `equalDOF` on UX/UY at the spring interface; rotation through the spring
- Spring elastic scale: `Ls = H/(3 α_mon α_I · pierCrackedFactor)` with `α_mon = 1/(1−1/η)` (Fiber full \(E_c,E_s\); \(L_s\) matches cracked \(I\))
- Pier \(EI\) (elastic members): uncracked transformed \(I_\mathrm{uncr}=I_g+(E_s/E_c-1)I_s\) (ring \(I_s=\tfrac12 A_s R_\mathrm{bar}^2\)), then \(I=\) `pierCrackedFactor` \(\times I_\mathrm{uncr}\) (default 0.5). Fiber hinges keep full \(E_c,E_s\). Set factor to 1 for uncracked elastic \(I\).

---

## Reference documents

Local PDFs in `reference/` (gitignored). BibTeX: `reference/references.bib`.

### Bridge geometry

| File | Cite key | Role |
|---|---|---|
| `Shin et al. - 2007 - … Liquefiable Soils.pdf` | `shinPerformanceBasedEvaluationBridges2007` | Bridge geometry |
| `Kramer - Using OpenSees … Liquefiable Soils.pdf` | `kramerUsingOpenSeesPerformanceBased` | Same bridge class; OpenSees / PBEE (PEER 2008/07) |
| `Mackie et al. - 2008 - … Benchmark Reinforced Concrete Bridges.pdf` | `mackieIntegratedProbabilisticPerformanceBased2008` | Benchmark RC bridges (PEER 2007/09) |
| `Mackie and Stojadinović - 2003 - Seismic Demands ….pdf` | `mackieSeismicDemandsPerformanceBased2003` | Seismic demands / PBEE (PEER 2003/16) |

### Hybrid simulation setup

| File | Cite key | Role |
|---|---|---|
| `Neumann - 2021 - … Real-Time Hybrid Simulati.pdf` | `neumannFluidStructureInteraction2021` | OSU MS thesis; RTHS / FSI |
| `Neumann et al. - 2023 - … Cascading Seismic and Tsunami Events.pdf` | `neumannHydrodynamicRealTimeHybrid2023` | Hydrodynamic RTHS demonstration |
| `Seki et al. - 2026 - … monopile offshore.pdf` | `sekiHydrodynamicRealtimeHybrid2026` | Hydro-RTHS; similitude |

In the Neumann RTHS setup the **physically tested pier** is represented numerically as a
very stiff `elasticBeamColumn` (`Iy × 1000` on that shaft) so flexure sits in
the base **rotational SSI spring** (`zeroLength` + `Steel01`) — same idea as
our `lumpedPlasticity` base hinge (η·EI beam + ZLS). Paper: deformations
mainly in the SSI spring; stiff link between base and top for equilibrium with
the pinned physical specimen.

### SSI boundary conditions

Switch: `soilBoundary` in `Parameters.tcl` — **`Shin`** | **`ASDEA`**.  
Full write-up (theory refs, Lysmer \(2c\,v\), ASDEA stages, Abell `setTime`): **`soil/Boundary.md`**.

| Mode | Scripts | Short form |
|------|---------|------------|
| Shin | `BuildSoilMesh` + `ActivateEQBoundary` → `BuildShinLysmer` | Thick FF + 3 rock Lysmers, \(F=2c\,v\) |
| ASDEA | `BuildASDEABoundary` (Stage 0) + `ActivateEQBoundary` (Stage 1) | `ASDAbsorbingBoundary2D` L/BL/B/BR/R; rock bottom; layer G0 on sides |

**Gravity BCs:** Shin = base `fix 1 1` (+ FF `equalDOF`); ASDEA = Stage 0 ring. Soil first (quad γ'/γ + Linear ponding). Pile/cap **nodes** exist (beams later). Ride-along MP is pile/cap retained, soil constrained. pile/cap→dup `equalDOF` stays during soil gravity, then is dropped and re-added after folding beam nodes only (offset frozen). Then Linear structure weight (`analysis/StructureGravityLoads.tcl`). Driver: `Run.tcl` (`runEQ` 0). After last `loadConst`, call `ActivateEQBoundary` (Shin: unfix UX + Lysmer; ASDEA: Stage 1; PyLiq/TzLiq → stage 1).

Shared velocity Path: `analysis/BuildVelSeries.tcl` (`gmVelFile`, `gmVelDT`).

#### Key references

**Shin / Lysmer**

- Shin et al. (2007); Kramer et al. (PEER 2008/07) — FF columns / liquefiable-bridge SSI context (`reference/`)
- Lysmer & Kuhlemeyer (1969) — dashpot \(c=\rho V_s A\)
- Joyner & Chen (1975) — outcrop input \(\Rightarrow F=2c\,v\)
- Elgamal et al. (2008); Zhang et al. (2008) — OpenSees Viscous + Path Lysmer pattern (some scripts scale by \(c\); we use \(2c\))
- Rock half-space under L5: \(V_s\approx 760\) m/s, \(\rho=2400\) kg/m³, force scale \(2c\)

**ASDEA**

- OpenSees *ASDAbsorbingBoundary* (Petracca & Camata, ASDEA); theory note Nielsen (2006)
- OpenSees manual: Stage 0 → Stage 1 via `setParameter … stage`; Path without `-startTime` unless scheduled otherwise

`-fx` is attached at element creation but **not applied in Stage 0** (only Stage 1 calls `addBaseActions`).

### Pier geometry (`Parameters.tcl`)

| Quantity | Value | Notes |
|---|---|---|
| `D_cyl` | 0.508 m (20.0 in) | Lab specimen (Neumann 2021; Neumann et al. 2023) |
| `H_cyl` | 3.02 m | Lab specimen (Neumann `Hc_m`; model pier 118.9 in) |
| `D_pier` | 1.2192 m (4.0 ft) | Shin et al. (2007); Mackie et al. (2008) |
| `cylinderSF` | **2.4** | `(4 ft = 48 in) / 20 in` |
| `H_pier` | **7.248 m** | `H_cyl · cylinderSF` |

Units: **N, m, s**. Elastic section: `structure/PierSection.tcl`.

---

## Proposed file layout

```
OSU_SSI_Bridge/
  NOTES.md
  README.md
  Parameters.tcl           # knobs (# <-- EDIT); TAGS CONVENTION
  Run.tcl                  # driver (runEQ 0|1); optional Overrides.tcl argv
  RunParallel.tcl          # OpenSeesMP driver; optional Overrides.tcl argv
  TestMatrix.csv           # campaign rows
  RunTestMatrix.py         # --row N -> Overrides.tcl (gitignored)
  BuildModel.tcl           # structure nodes + soil quads + SSI springs as needed;
                           # pier/deck/cap/pile beam-columns after soil gravity
  reference/
  analysis/                # README.md: gravity stages; settings in Run.tcl
  structure/               # README.md: Nodes vs Elements
  soil/                    # Profiles.md, Boundary.md, Build*, analysis
  plot/out/                # figures + EQ dumps (serial|parallel)
  tmp/                     # local scratch; safe to delete
```

Knobs in `Parameters.tcl` (`# <-- EDIT`); IDs in the TAGS CONVENTION section.
`tmp/` is gitignored; `Run.tcl` and `PlotModel.tcl` do not source it.
Matrix runs: `python RunTestMatrix.py --row N` then pass `Overrides.tcl` as argv to `Run.tcl` / `RunParallel.tcl` (`overridesON` defaults to 1; no file → forced off).
Wave Name (Storm Wave / Big Tsunami) prototype vs lab depths/heights/periods: `WaveCatalog.csv` (documentation only for now).

---

## Build stages

1. **Scaffold** — notes, units  
2. **Pier geometry + section** — done  
3. **Above-ground pier base hinge** — done (`lumpedPlasticity`)  
4. **Pile cap** — done (stiff frame)  
5. **Piles** — done (3 ×2 pipe shafts)  
6. **Soil layers** — done (`soil/Profiles.md` + mesh/materials)  
7. **SSI coupling** — done (coincident p-y/t-z/q-z, cap faces + soffit q-z; `Ge=0.69`)  
8. **Recorders + verification** — in progress (`Run.tcl`: gravity or gravity+EQ)

Driver: `Run.tcl` (`BuildModel.tcl` assembly). Knobs in `Parameters.tcl` (`# <-- EDIT`).
`runEQ` 0 = gravity+eigen; 1 = gravity then EQ (`Run.tcl`).
Analysis: soil gravity in `analysis/SoilGravity.tcl`; fold in `FoldStructNodes.tcl`;
`numberer` / `analyze` in `Run.tcl`. Recorders: `EQRecorders.tcl`.
Status: pier + deck + pile cap + piles + soil + springs.
Switch in `Parameters.tcl`:
- `pierEleType`: `elasticBeamColumn` | `forceBeamColumn` | `lumpedPlasticity`
- `pileEleType`: `elasticBeamColumn` | `dispBeamColumn`
- `soilProfile`: `1` | `2` | `3` | `4`
- `soilBoundary`: `Shin` | `ASDEA`
- `soilEleType`: `quad` | `SSPquad`
- `soilConstitutive`: `inelastic` | `elastic` (`ElasticIsotropic3D` from skeleton \(G_r,B_r\); sands still FSP)
- `pileSpring`: `inelastic` | `elastic` (\(k\approx p_\mathrm{ult}/y_{50}\), …) | `none` (equalDOF pile↔soil)

Edit knobs in `Parameters.tcl` (`# <-- EDIT`), then re-run.
Mesh / figures: `OpenSees PlotModel.tcl` (build + JSON dumps + Python plotters).
JSON only: `set ::plotSkipPython 1` then source `PlotModel.tcl`.
Plot: `python3 plot/PlotModelSketch.py` → `plot/out/profile{N}/elevation/{Shin|ASDEA}/elevation.png`.

Gravity: `runEQ 0` (elastic soil gravity → stage 1 → plastic → fold structure nodes → structure weight; Shin base `fix` or ASDEA Stage 0; then `loadConst` + eigen in `Run.tcl`). With `runEQ 1`, then `ActivateEQBoundary` (Lysmer/ASDEA + PyLiq/TzLiq stage 1) and the transient loop in `Run.tcl`. Details: `soil/Boundary.md`.

### Pier mesh (current)

```
elastic / forceBeamColumn:
  node 1 (0,0) ── ele 2 ── node 5 (0, H_pier)

lumpedPlasticity (base ZLS + η·EI beam):
  1 ──ZLS-I (ele 1)── 2 ════ η·EI (ele 2) ════ 5
  equalDOF UX/UY: 1→2; SteelMPF
  ZLS scale: Ls = H/(3 α_mon α_I · pierCrackedFactor) on E and through peak;
             (ε_u−ε_peak)·Lp post-peak (Priestley)
```

Fiber path: graded strips + rebar on section y (`CircleStripFibers.tcl`).
Pier elastic \(EI\) (elasticBeamColumn / FBC mid / \(\eta I\)): uncracked transformed × `pierCrackedFactor` (default 0.5). Fiber materials use full \(E_c,E_s\); lumped \(L_s\) includes the cracked factor.
Priestley / Caltrans \(L_p = 0.08 H + 0.022 f_y d_b\) (\(f_y\) in MPa, lengths in m; \(0.022\) and \(0.044\) carry units \(\mathrm{MPa}^{-1}\); same as Mackie 2003 \(0.15/0.3\) with ksi·in).

### Nodal rotary mass

No element `-mass`. Translation is lumped; RZ uses the following.

**Pier** (ends 1 and 5; lumpedPlasticity splits base mass across the spring pair):

\[
m = \tfrac12 \rho_L H,\qquad
I_\mathrm{rot} = \frac{\rho_L H^3}{105}
\]

\(I_\mathrm{rot}\) is the diagonal rotational entry of an elastic consistent-mass matrix (\(\rho L\cdot 4L^2/420\)). Order \(\sim m\) (kg·m² with \(k\sim 1\,\mathrm{m}\)); small vs \(m H^2\).

**Pile cap** (11 frame nodes: 3×3 pile grid + soffit mid-bay). Tributary rectangles scaled to \(m_\mathrm{cap}\); \(I_{\mathrm{rot},i}=m_i(dx_i^2+dy_i^2)/12\). Outer \(dx\) reaches \(\pm W/2\) (1.5 ft overhang beyond the frame at \(\pm s\)).

Shared nodes stack mass with `IncrMass` (`structure/IncrMass.tcl`): deck soffit CL on pier top; cap TC on pier base; pile heads on cap bottom.

### Pile cap

Stiff `elasticBeamColumn` frame (steel \(E\), \(A=HW\), \(I=HW^3/12\)).
Mackie et al. (2008) Type 1A / Ketchum 3×2 under 4 ft column:
\(H=3.25\,\mathrm{ft}\), \(W=15\,\mathrm{ft}\) (transverse), \(L=10\,\mathrm{ft}\)
(out-of-plane), \(s=6\,\mathrm{ft}\) (3\(D\) c/c). Top center **is** pier base
(node 1); frame and face springs at \(\pm s\); soffit mid-bay at \(\pm s/2\); eles **1101+**. No `equalDOF`. Face Py/Tz on TL/ML/BL and TR/MR/BR; soffit q-z on 5 bottom stations.

Names: T/M/B = top / mid / bot; L/C/R = left / center / right; BML/BMR = bot mid-bay.

```
  y=0     TL ---------- 1(=TC) ---------- TR
           |           / | \              |
  y=-H/2  ML ---------- MC ---------- MR
           |         /    |    \           |
  y=-H    BL -- BML -- BC -- BMR -- BR   pile heads at ±s, 0
          -s   -s/2     0    s/2     s
```

### Piles

Open steel pipe: \(D=2\,\mathrm{ft}\), \(t=0.5\,\mathrm{in}\), three shafts at
cap bottom nodes. Each shaft × `n_pile_row=2` (out-of-plane row condensed into
\(A\) and \(I\)). Length \(L_\mathrm{pile}=60\,\mathrm{ft}\) (Mackie); 20 segments
→ \(\Delta y=3\,\mathrm{ft}\) (even on 60 ft; SI run still N·m·s).
Plots: `python3 plot/PlotModelSketch.py`, `python3 plot/PlotFiberSections.py`,
`python3 plot/PlotPileSprings.py`, `python3 plot/PlotSoilProfile.py`
→ `plot/out/profile{N}/` (elevation split by Shin|ASDEA; springs/soil/fibers by profile only).

| `pileEleType` | Section |
|---|---|
| `elasticBeamColumn` | \(A\), \(E_s\), \(I\) × `n_pile_row` |
| `dispBeamColumn` | Fiber: graded tube strips (`circularTubeFiberStripsGraded`); strip areas × `n_pile_row` |

Tags: tip/interior nodes **2001+** (left), **2101+** (center), **2201+** (right);
eles **2100+**.

Nodal mass on piles: half segment on tip and head (`IncrMass` onto cap \(m_i\));
full mid. \(I_\mathrm{rot}=\rho_L\,\mathrm{d}y^3/105\) per segment end.

### Model sketch / fiber plots (`plot/`)

OpenSees build dumps JSON; Python only reads those files.

| Dump (Tcl) | Plot (Python) | Output |
|---|---|---|
| `plot/DumpModelSketch.tcl` | `plot/PlotModelSketch.py` | `plot/out/profileN/elevation/{BC}/elevation.png` |
| `plot/DumpFiberSections.tcl` | `plot/PlotFiberSections.py` | `plot/out/profileN/fibers/fiber_*.png` |
| `plot/DumpPileSprings.tcl` | `plot/PlotPileSprings.py` | `plot/out/profileN/pile_springs/*.png` |
| `plot/DumpSoilProfile.tcl` | `plot/PlotSoilProfile.py` | `plot/out/profileN/soil_profile/*.png` |

Fiber figures show the section **as modeled**: graded strips; pier rebar fibers
on \(z=0\) with combined \(A\) (merged ±z), not individual bar positions on the ring.

When soil arrives: extend the dumps. Plotters stay dumb.

---

## Open decisions

- [x] Pier diameter / scale — **D_pier = 4 ft**, **cylinderSF = 2.4**, **H_cyl = 3.02 m**, **H_pier = 7.248 m**
- [ ] Sandwich `η` / spring `α`’s
- [ ] Single pier only vs pier + abbreviated deck
- [x] Soil: 2D continuum; profiles 1–4 (L3 sand / clay+sand / all clay / soft L2+L3a) — `soil/Profiles.md`
- [x] Units: **N, m, s**
- [ ] Tag ranges (structure vs soil)

---

## Changelog

| Date | Note |
|---|---|
| 2026-08-10 | Initial notes |
| 2026-08-11 | Reference PDFs + `references.bib` |
| 2026-08-11 | Writing-style rule |
| 2026-08-11 | Pier geometry in `Parameters.tcl` (D, cylinderSF, H); elastic section |
| 2026-08-11 | Dropped cross-links to other local OpenSees models in code/notes |
| 2026-08-11 | 2-node pier; pierEleType elastic | force; graded fiber strips |
| 2026-08-11 | forceBeamColumn → ConcentratedCurvature + Priestley Lp |
| 2026-08-11 | lumpedPlasticity base ZLS + η·EI; SteelMPF |
| 2026-08-11 | ZLS hybrid scale: Ls elastic through peak; Lp post-peak |
| 2026-08-11 | Pile cap 3×3 stiff frame (tagShift 1000); base spring mass split |
| 2026-08-11 | Piles: 3×2 pipes, elastic / dispBeam Fiber strips; s=6 ft |
| 2026-08-11 | Sketch from OpenSees dump JSON (Dump/PlotModelSketch) |
| 2026-08-11 | Plots moved to `plot/`; fiber figs as-modeled (rebar on z=0) |
| 2026-08-11 | Cap 15×10×3.25 ft, L_pile=60 ft, dy=3 ft (Mackie Type 1A) |
| 2026-08-11 | Soil profiles: L2/L5 + L3 sand(3 PDMY) or clay; `soil/Profiles.md` |
| 2026-08-13 | lumpedPlasticity: single base spring (dropped top ZLS) |
| 2026-08-13 | Model switches: soilConstitutive, pileSpring; PyLiq/TzLiq stage 1 at ActivateEQBoundary |
| 2026-08-13 | Pier EI: uncracked transformed × `pierCrackedFactor` (default 0.5; elastic only; Ls ÷ factor) |
| 2026-08-15 | `Run.tcl` named stages; gravity split (SoilGravity / FoldStructNodes); settings in driver |
| 2026-08-14 | Dual-face pile springs (±R, ½ p-y/t-z); ZL soil↔dup(ndf2) twin of iface/cap; soffit q-z on BL/BC/BR |
| 2026-08-14 | Cap frame to ±(s+R); mid chord; soffit meshed (7 q-z stations, tributary \(Q_\mathrm{ult}\)) |
| 2026-08-17 | Cap frame to ±s (drop outer-face ring); 5 soffit q-z; soil must-hit pile axes only (3 ft inner band) |
| 2026-08-14 | Cap face py/tz by trib height: Mokwa \(b=L\) (½/side); tz \(c\cdot 2H(W+L)\) |
| 2026-08-15 | Gravity: soil Linear 0→1, hold+stage 1, then Linear structure |
| 2026-08-15 | Gravity ride-along: pile retained, soil constrained (star with pile→dup) |
| 2026-08-15 | 2D piles: coincident ZL p-y/t-z only (no radial stubs or shaft voids) |
| 2026-08-17 | `recordersON 2`: soil nodes are geometry only (`disp_nodes.txt` holds the `window_disp` columns) |
| 2026-08-17 | `uy_dup - uy_pile` = 0.2--0.8 mm is the fold datum shift, not a broken tie: the pile keeps its gravity settlement in `nodeCoord`, the unfolded dup keeps it in `nodeDisp`. coord + disp agree to 6e-8 m |
| 2026-08-17 | EQ `constraints Auto` (both drivers): ZLS `deformation` = relative nodal disp to 4e-8 m, same as `Transformation`, at `Plain`'s runtime. `Plain` was 1e-3 m off in UY |
| 2026-08-17 | Why `Plain` fails: no constraint FE_Element, `PlainNumberer` just shares the retained DOF's equation number, so a tied pair cannot hold two different total displacements. The 0.79 mm fold offset then kicks a vertical transient at EQ start (`uy_pile` 1.7e-3 m for the first 0.2--0.5 s, then ~5e-5 m). UX, which has no offset, is clean. The `sp -const -subtractInit` holds fine under all three |
