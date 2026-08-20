# analysis/

Gravity and EQ. Analysis settings sit next to each `analyze`: soil gravity in
`SoilGravity.tcl`; structure weight, eigen, and EQ in `Run.tcl`. Recorders
stay in this folder.

| File | Role |
|---|---|
| `SoilGravity.tcl` | Elastic soil gravity, stage 1, plastic gravity (LoadControl here) |
| `FoldStructNodes.tcl` | Fold beam nodes to soil settlement, then `BuildStructElements.tcl` |
| `StructureGravityLoads.tcl` | Nodal -mg on the structure |
| `HoldPierBase.tcl` | After gravity: freeze pier-base UX/UY (`-subtractInit`); gated by `holdPierON` |
| `RayleighDamping.tcl` | Regions + `αM`, `βKinit` (pier / hinge last) |
| `GravityHelpers.tcl` | Spring kine print; nearest soil node |
| `WaterSurfaceLoad.tcl` | Ponding if `h_water > 0` |
| `BuildVelSeries.tcl` | Outcrop velocity Path (PEER VT2 or dummy zeros) |
| `EQRecorders.tcl` | `recordersON` 0 off / 1 full window / 2 center+near-FF columns / 3 nine SSI (+ same columns). Near-FF at `eqFFColumnFrac*L_half` (not Shin thick FF). `-dT` = GM DT. `window_nodes.txt` = geometry, `disp_nodes.txt` = `window_disp` columns |
| `EigenAfterGravity.tcl` | Mode JSON/PNG after `eigen` in `Run.tcl` |
