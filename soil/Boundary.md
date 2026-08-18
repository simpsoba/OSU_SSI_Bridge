# Soil domain boundaries — Shin vs ASDEA

Switch: `soilBoundary` in `Parameters.tcl` (`Shin` | `ASDEA`).  
Same **near-field** continuum (`|x| ≤ L_half`) and same base velocity Path (`analysis/BuildVelSeries.tcl`).  
Shin adds thick FF columns of width `w_FF` out to `L_half+w_FF`; ASDEA stops at `L_half` and puts the ASD ring outside that face.

Gravity driver: `Run.tcl` (`runEQ` 0) — soil gravity (`SoilGravity.tcl`), fold (`FoldStructNodes.tcl`), structure weight + `loadConst` in the driver (Shin base `fix` or ASDEA Stage 0).  
After gravity: `source soil/ActivateEQBoundary.tcl` (ASDEA/Lysmer + PyLiq/TzLiq → stage 1).

### Gravity BCs (by FF approach)

| | **Shin** | **ASDEA** |
|---|----------|-----------|
| Base | `fix 1 1` on all base nodes (`BuildSoilMesh`) | ASD bottom + BL/BR in **Stage 0** (penalty “fix”) |
| Sides | Thick FF columns + face `equalDOF`; no lateral fix except base | ASD **L** / **R** Stage 0 |
| Surface | Free | Free (no top ASD) |
| EQ switch | `ActivateEQBoundary` → unfix UX, 3 Lysmers + \(2c\,v\) | `ActivateEQBoundary` → `setParameter … stage` **1** |

Body forces stay on for EQ (`loadConst`); only the artificial boundary changes.

### Transformation constraint handler

Gravity and EQ both use `constraints Transformation`. OpenSees builds the retained/constrained DOF map when the **analysis is created**, so:

1. `wipeAnalysis` before any `remove sp` / new `equalDOF` (Shin) or before relying on a new analysis after ASDEA stage change.
2. `equalDOF $retained $constrained dof…` — primary = retained; put Lysmer / \(2c\,v\) loads on retained nodes only. Shin FF faces: retained = outermost column; base elevation UX only (UY on `fix`); above base UX+UY.
3. Lysmer: one fully fixed ghost — Viscous — soil base (both ndf=2; no mate/`equalDOF`).
4. Create the EQ `analysis` **after** `ActivateEQBoundary` so Transformation sees the final SP/MP set.

---

## Shin (`soilBoundary "Shin"`)

### What we do

1. Thick free-field (FF) columns on L/R (`t_FF = t_FF_factor · t_soil`, width `w_FF` beyond `L_half`), face `equalDOF` to the near field.
2. Gravity: base `fix 1 1` (see table above). Scripts: `analysis/SoilGravity.tcl`, then fold + structure weight in `Run.tcl`.
3. EQ (`BuildShinLysmer.tcl`): release base UX; three Lysmer–Kuhlemeyer dashpots + force \(F = 2c\,v(t)\):
   - near-field base primary (all other NF base nodes `equalDOF` in UX), \(c = \rho_r V_{s,r}\,L_{\mathrm{nf}}\,t_{\mathrm{soil}}\)
   - left FF outer base, \(c = \rho_r V_{s,r}\,w_{FF}\,t_{FF}\)
   - right FF outer base, same \(c\)
4. \(v(t)\): integrated, baseline-corrected rock outcrop velocity (`gmVelFile`).

Rock props: `rockVs`, `rockRho` (half-space under L5; \(V_s \approx 760\) m/s, \(\rho = 2400\) kg/m³).

### References (theory + practice)

| Role | Citation | Notes |
|------|----------|--------|
| FF column idea / liquefiable-bridge SSI context | Shin et al. (2007); Kramer et al. (PEER 2008/07) | Thick outer columns as free-field; geometry also from Shin/Mackie |
| Lysmer dashpot | Lysmer & Kuhlemeyer (1969) | \(c = \rho V_s A\) radiation damper |
| Outcrop force \(F = 2c\,v\) | Joyner & Chen (1975); common OpenSees soil-column practice | Incident + reflected; ASD element uses the same factor internally |
| OpenSees Lysmer + Path pattern | Elgamal et al. (2008); Zhang et al. (2008); OpenSees Viscous / Path | Base: free UX, `equalDOF` to one retained primary, Viscous + Path; FF columns with their own Lysmers. Some scripts scale Path by \(c\); **we use \(2c\)** (outcrop) |
| Rock \(V_s,\rho\) | Firm rock half-space under L5 | `rockVs ≈ 760` m/s, `rockRho = 2400` kg/m³; force scale \(2\cdot\rho\cdot V_s\cdot A\) |

Local PDFs (geometry / SSI context): `reference/` — Shin 2007, Kramer PEER 2008/07, Elgamal et al. 2008, Zhang et al. 2008.  
BibTeX: `reference/references.bib`.

---

## ASDEA (`soilBoundary "ASDEA"`)

### What we do

1. No thick Shin FF; continuum ends at \(\pm L_\mathrm{half}\) (same NF as Shin). ASD ring outside.
2. Extruded `ASDAbsorbingBoundary2D` ring: **L, BL, B, BR, R** (`BuildASDEABoundary.tcl`).
3. Props: **bottom / corners = rock** (`rockG`, `rockNu`, `rockRho`); **L/R = layer** `soilG0`, `asdeaNu`, `soilRho` at face mid-height.
4. Bottom / BL / BR created with `-fx $tsTag_velBase` (same velocity series as Shin).
5. Created in **Stage 0** (penalty “fix” for gravity). After gravity: `setParameter -val 1 -ele … stage` → Stage 1 (absorbing + base traction).

### References

| Role | Citation | Notes |
|------|----------|--------|
| Element (OpenSees) | Petracca & Camata (ASDEA); OpenSees manual *ASDAbsorbingBoundary* | 2D/3D; Stage 0 → Stage 1 via `setParameter … stage` |
| Theory basis cited by manual | Nielsen (2006) | Absorbing BC assembly: free-field + Lysmer + traction transfer |
| Usage note (Path scheduling) | OpenSees manual; optional Path `-startTime` + `setTime` | Bottom ASD with `-fx`; rock-like \(G,\nu,\rho\) on ASD; `setParameter stage`. Offset Path only if the series is scheduled that way |

---

## `-fx` in Stage 0 (ASDEA)

**Yes — ignored during Stage 0.**

In `ASDAbsorbingBoundary2D::getResistingForce` / `…IncInertia`, Stage 0 only calls `addRPenaltyStage0`.  
`addBaseActions` (reads `-fx`/`-fy` via `TimeSeries::getFactor(domainTime)`) runs only in **Stage 1** (absorbing), together with Lysmer and free-field terms.

The series must still exist at element construction (OpenSees copies the TimeSeries into the element). It is simply not sampled until Stage 1.

---

## Path `-startTime` and `setTime`

**Not an ASD requirement.** The element only does `TimeSeries::getFactor(domain->getCurrentTime())` in Stage 1. Official OpenSees ASD example: `loadConst -time 0.0` → `setParameter … stage` → transient from \(t=0\) (Path has no `-startTime`). No special `setTime`.

Some workflows offset the Path for scheduling only:

```tcl
timeSeries Path 11 -dt 0.0025 -values {…} -factor 0.5 -startTime 7.0
…
loadConst -time 1.0          ;# after gravity
# activate ASD stage → 1
setTime 7.0                  ;# align domain clock with Path -startTime
```

OpenSees `Path` with `-startTime 7.0` maps `values[0]` to domain time **7.0**. For \(t < 7\), `getFactor(t)` is 0. `setTime 7.0` so the first EQ step hits the start of that series.

**Our default:** Path without `-startTime` (starts at 0) and `loadConst -time 0.0` after gravity — Stage 1 + transient at \(t = 0\). Use `-startTime` + `setTime` only if the velocity file is scheduled that way.

---

## ASDEA forum notes

Local copy: `soil/ASD-absorbing-boundary-forum-2026.pdf`  
Online: [asdeasoft.net … p=11058](https://asdeasoft.net/forum/viewtopic.php?p=11058) (thread *ASD absorb Boundary condition*, Aug 2025–Jul 2026).

- **Linear elastic by design.** Standard ASD uses fixed \(G,\nu,\rho\). That matches our rock bottom + layer `soilG0` sides (elastic FF/Lysmer props, not live PDMY \(G\)).
- **Optional `-mat` (nonlinear G update):** newer feature; STKO says assign an nDMaterial instead of \(G,\nu,\rho\). Implemented/tested first on **3D**; needs a recent solver (they pointed to **OpenSees 3.7.2** — older installers reject `-mat`). **2D support was still unclear** as of the last post (Jul 2026); unanswered whether free-field becomes fully nonlinear or only \(G\) is refreshed for dashpots.
- **Depth-varying \(G\):** forum asked average vs sublayers; no detailed formula reply. Our practice (per-elevation side props from the adjacent layer) is the practical answer without waiting on `-mat`.
- **Gravity / initial stress (manual ASD):** asked whether it follows a conventional model; STKO did not answer beyond pointing to 3.7.2 for `-mat`. We still follow the OpenSees manual: Stage 0 for gravity, then `setParameter … stage` → 1 before EQ.
- Nothing in that thread requires special `setTime` or changes Stage 0/`-fx` behavior.
