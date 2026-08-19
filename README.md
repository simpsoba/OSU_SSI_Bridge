# OSU SSI Bridge

2D transverse pier-bridge SSI. Geometry: Shin et al. (2007), Mackie et al. (2008). Lab cylinder: Neumann (2021); Neumann et al. (2023). Units **N, m, s**.

Parameters of the model are mostly changed in `Parameters.tcl` (`# <-- EDIT`). 
The main file to run is: `Run.tcl`

To obtain the latest version of the model from GitHub
```bash
git clone git@github.com:simpsoba/OSU_SSI_Bridge.git
OpenSees.exe Run.tcl # for serial runs
mpiexec -n 4 OpenSeesMP.exe RunParallel.tcl # for parallel runs
```

Before `mpiexec` (Windows `cmd`), set these so Intel MPI pins one rank per core and does not oversubscribe OpenMP/MKL. Point `TCL_LIBRARY` and `OpenSeesMP.exe` at the OpenSees-CUDA build:

```bat
set TCL_LIBRARY=C:\projects\RTHS-CUDA\OpenSees\build-mp\lib\tcl8.6
set I_MPI_PIN=on
set I_MPI_PIN_CELL=core
set I_MPI_FABRICS=shm
set OMP_NUM_THREADS=1
set MKL_NUM_THREADS=1

cd /d C:\path\to\OSU_SSI_Bridge
mpiexec -n 4 C:\projects\RTHS-CUDA\OpenSees\build-mp\Release\OpenSeesMP.exe RunParallel.tcl
```

Ranks other than 0 send OpenSees `opserr` (C++ warnings/errors) to `opensees.rankN.log` via `logFile … -noEcho`, and Tcl `puts` to stdout is dropped. Rank 0 still prints to the console. Tcl `error` traces from a failing rank can still show up through MPI.

Named tags live in `Parameters.tcl` (TAGS CONVENTION). Shared nodes (same tag, stacked mass, no `equalDOF`): cap TC = pier base (1); deck soffit BC = pier top (5); pile heads = cap BL / BC / BR.

### Structure — pier (default `lumpedPlasticity`)

Node 3 is unused.

```
  5  pier top / deck soffit BC     y = H_pier
  |  ele 3  top rotational spring (ZLS-J)
  4  ZLS-J inner                   (same coords as 5)
  ║  ele 2  eta*EI beam
  2  ZLS-I inner                   (same coords as 1)
  |  ele 1  base rotational spring (ZLS-I)
  1  pier base / cap TC            y = 0
```

```
# Nodes
#   1 : pier base, cap top center
#   2 : base-spring inner (lumpedPlasticity only)
#   4 : top-spring inner  (lumpedPlasticity only)
#   5 : pier top, deck soffit center
#
# Elements
#   1 : base rotational spring     (nodes 1 -> 2)
#   2 : eta*EI pier beam           (nodes 2 -> 4)
#   3 : top rotational spring      (nodes 4 -> 5)
```

`elasticBeamColumn` / `forceBeamColumn`: nodes 1 and 5 only; ele 2 is `1 -> 5`. Eles 1 and 3 are not created.

### Structure — deck (PT39 frame)

```
  BarL 3009                                    BarR 3010
     |                                            |
  TL 3004 -- TLi 3005 -- TC 3006 -- TRi 3007 -- TR 3008     y = H_pier+dd
     |          |          |          |          |
             BL 3001 ---- BC 5 ---- BR 3003               y = H_pier
```

```
# Nodes
#   3001 : soffit BL     5 : soffit BC (pier top)     3003 : soffit BR
#   3004 : top TL     3005 : top TLi     3006 : top TC
#   3007 : top TRi    3008 : top TR
#   3009 : barrier L  3010 : barrier R
#
# Elements (elasticBeamColumn, 3100--3110)
#   3100 : soffit BL-BC     (3001 -> 5)
#   3101 : soffit BC-BR     (5 -> 3003)
#   3102 : top TL-TLi       (3004 -> 3005)
#   3103 : top TLi-TC       (3005 -> 3006)
#   3104 : top TC-TRi       (3006 -> 3007)
#   3105 : top TRi-TR       (3007 -> 3008)
#   3106 : outer web L      (3001 -> 3005)
#   3107 : outer web R      (3003 -> 3007)
#   3108 : center web       (5 -> 3006)
#   3109 : barrier L        (3004 -> 3009)
#   3110 : barrier R        (3008 -> 3010)
```

### Structure — pile cap

```
  y=0     TL 1021 ---------- TC 1 ---------- TR 1023
           |                / | \              |
  y=-H/2  ML 1024 ---------- MC 1025 ---------- MR 1026
           |              /    |    \           |
  y=-H    BL 1027 -- BML 1036 -- BC 1028 -- BMR 1037 -- BR 1029
          -s        -s/2         0          s/2         s
```

```
# Nodes
#   1021 : TL     1 : TC (pier base)     1023 : TR
#   1024 : ML     1025 : MC              1026 : MR
#   1027 : BL (left pile head)
#   1028 : BC (center pile head)
#   1029 : BR (right pile head)
#   1036 : BML    1037 : BMR
#
# Elements (elasticBeamColumn, 1101--1116)
#   1101--1102 : left  vertical   TL-ML, ML-BL
#   1103--1104 : center vertical  TC-MC, MC-BC
#   1105--1106 : right vertical   TR-MR, MR-BR
#   1107--1108 : top horizontal   TL-TC, TC-TR
#   1109--1110 : mid horizontal   ML-MC, MC-MR
#   1111--1114 : soffit           BL-BML-BC-BMR-BR
#   1115--1116 : diagonals        TL-BR, TR-BL
```

### Piles (3 shafts × `nSeg_pile` = 20)

Heads are the cap bottom nodes. Below the head: `node = 2000 + ip*100 + iy` (`ip` = 0 left, 1 center, 2 right; `iy` = 1..20, tip at 20). Elements: `ele = 2100 + ip*20 + (iy-1)`, head → first station, then down the shaft.

```
# Nodes (head is a cap node; tip = iy 20)
#   left:   head 1027,  2001 .. 2020 (tip)
#   center: head 1028,  2101 .. 2120 (tip)
#   right:  head 1029,  2201 .. 2220 (tip)
#
# Elements
#   left:   2100 .. 2119     (1027 -> 2001 -> ... -> 2020)
#   center: 2120 .. 2139     (1028 -> 2101 -> ... -> 2120)
#   right:  2140 .. 2159     (1029 -> 2201 -> ... -> 2220)
```

---

## Where to edit

### Rayleigh damping

`analysis/RayleighDamping.tcl`

Damping is applied by `region`. Soil, far field boundary, SSI springs, piles, cap, and deck get near-zero damping. The pier beam, then the pier hinges (lumped plasticity), get the full `αM`, `βKcomm`. The last region that owns a node wins `αM`, so keep the pier groups last.

To turn a group on or off, swap `$aOff`/`$bOff` with `$alphaM`/`$betaKcomm` on that `-rayleigh` line.

### Soil mesh

`Parameters.tcl`, there is a list called `soilDxBands` that defines the soil mesh. Two coarser lists sit commented below it. Uncomment one list and comment the others.

Vertical size is `dy_soil` (keep equal to pile `dy`). Layer materials: `soil/Profiles.md`. Builder: `soil/BuildSoilMesh.tcl`.

### Ground motion

`Parameters.tcl`, Ground motion / EQ. Two records in `ground-motion/`:

| Record | Folder / file |
|---|---|
| Tohoku 2011, KiK-net FKSH19 borehole NS | `Tohoku2011-FKSH/FKSH19.NS1.VT2` (default) |
| El Centro 1940 Array #9, 180 (NS) | `ImperialValley1940-ElCentro/RSN6_IMPVALL.I_I-ELC180.VT2` |

Uncomment one `gmDir` / `gmVelFile` pair. Scale: `gmScaleFactor`. The Path is built in `analysis/BuildVelSeries.tcl` (PEER VT2 cm/s → m/s).

### Gravity and mass on the structure and piles

Mass is nodal (no element `-mass`). Densities and sizes are in `Parameters.tcl`; the lump masses are written in:

| Piece | File |
|---|---|
| pier | `structure/BuildPierNodes.tcl` |
| deck | `structure/BuildDeckNodes.tcl` |
| pile cap | `structure/BuildPileCapNodes.tcl` |
| piles | `structure/BuildPilesNodes.tcl` |


Structure gravity load: `analysis/StructureGravityLoads.tcl` walks `structNodeTags` and applies −mg from that mass. Soil self-weight is separate (`analysis/SoilGravity.tcl`).

After gravity, `Run.tcl` adds 1 kg / 0.1 kg·m² on every node so **M** is not singular.

### SP at the pier base (before EQ)

`analysis/HoldPierBase.tcl`, sourced from `Run.tcl` after `loadConst`. `sp -const -subtractInit` holds pier-base UX and UY at the gravity displacement (and the lumped-plasticity inner ZLS node). RZ stays free.

To drop the hold, stop sourcing that file (or comment the `sp` lines).

### Elastic vs inelastic

Switches in `Parameters.tcl` (Model switches, plus pier / pile type):

| Piece | Knob | Builder |
|---|---|---|
| soil quads | `soilConstitutive` `inelastic` \| `elastic` | `soil/BuildSoilMaterials.tcl` |
| SSI springs (p-y / t-z / q-z together) | `pileSpring` `inelastic` \| `elastic` \| `none` | `soil/BuildSoilSprings.tcl` |
| piles | `pileEleType` `elasticBeamColumn` \| `dispBeamColumn` | `structure/BuildPileSection.tcl` |
| pier | `pierEleType` `elasticBeamColumn` \| `forceBeamColumn` \| `lumpedPlasticity` | `structure/PierSection.tcl` |

### Fiber sections (pier springs and piles)

Strip counts and materials in `Parameters.tcl` (`nFiberY_pier`, `nFiberEdge_pier`, `nFiberY_pile`, `nFiberEdge_pile`; \(f'_c\), \(f_y\), bar layout).

| Piece | File |
|---|---|
| pier ZLS / forceBeamColumn hinges | `structure/PierSection.tcl` (`pierBuildFiberHinge`) |
| pile tube | `structure/BuildPileSection.tcl` |
| strip geometry | `structure/CircleStripFibers.tcl` |

Figures: `OpenSees PlotModel.tcl` → `plot/out/profile{N}/fibers/`.

### Free-vibration duration

`eqFreeVibT` (seconds) in `Parameters.tcl` defines the number of seconds that are run in free vibration after the earthquake ends.

### EQ progress / timings

`eqPrintON` in `Run.tcl` and `RunParallel.tcl` (not `Parameters.tcl`). Default `1` prints analysis time, wall-clock elapsed, and pier-top disp every `eqPrintDt` s (debug). Set `eqPrintON 0` to silence the loop. Ignored when `realTimeON 1`. The one-line `EQ done` summary still prints elapsed at the end.

### Recorder folder / GM start / realTimeON

Same files (`Run.tcl`, `RunParallel.tcl`):

- `outDIR` — recorder folder (`trial1`, `runA`, …). Relative to the process cwd. `""` uses the auto path `plot/out/profile{N}/eq/{serial|parallel}/...`.
- `gmStartTime` — Path `-startTime` (s). `0` omits it. Domain clock stays at 0 after gravity; the GM is silent until that time. When `realTimeON` is 0, `eqNstepsAll` covers `gmStartTime + Trec + eqFreeVibT`.
- `realTimeON 1` — OpenFresco, no recovery, `realTimeNsteps`, no `eqPrintON`. Default `0` is the usual EQ loop. Needs `pierEleType lumpedPlasticity`. `eleTag_exp` (default 101) is the generic experimental element.

### OpenFresco (`expElement`)

Gated by `realTimeON`. Create the experimental element **before** `numberer` / `system` / `analysis Transient`. Recorders go after `EQRecorders.tcl`.

**Parallel:** pin the pier on rank 0, then create `expElement` on rank 0 after `partition`, then `barrier`:

```tcl
partition -keepOnRank 0 3 \
	$eleTag_pier_botSpr $eleTag_pier $eleTag_pier_topSpr
# rank 0: expElement on $nodeTag_pierTopZeroLengthInner
# barrier
# then numberer / system / constraints / analysis Transient
```

Without OpenFresco (`realTimeON 0`), `RunParallel.tcl` just calls `partition`. Do not create the analysis first: `partition` only sees elements already in the domain, and the numberer is built from the mesh that remains after the split.

---

## Run / figures

`OpenSees PlotModel.tcl` (build + dump + Python). JSON only: `set ::plotSkipPython 1` then source `PlotModel.tcl`.

EQ dumps: `outDIR` if set, else `plot/out/profile{N}/eq/{serial|parallel}/{BC}/{soilEle}/{pierEleType}/`
Modes: `plot/out/profile{N}/elevation/{BC}/modes/{soilEle}/{pierEleType}/`
Gravity deformed: `…/gravity/{soilEle}/{pierEleType}/gravity_deformed.png`
