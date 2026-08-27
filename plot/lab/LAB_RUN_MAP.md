# OSU SSI Bridge — lab run map (Wed 2026-08-19 + Fri 2026-08-21)

Narrative companion to the machine-readable as-run index. **Do not duplicate the
pair table here** — use the CSV.

## Canonical sources

| File | Role |
|------|------|
| **`plot/lab/TestMatrix_lab_runs.csv`** | One row per as-run. Campaign/model knobs first (incl. **`DOFs`** = OpenSees `systemSize` after gravity); trailing **`DateTime`**, **`DumpFolder`**, **`MatFile`**, **`LabTrial`**, **`Note`**. Leading **`Test`** = `W##` (Wed) / `F##` (Fri RTHS Trial) / `Fd##` (dry lunch) / `Fx##` (abort / no Trial). Git-tracked; working copy: `OSU_SSI_BRIDGE_DATA_LOCAL/TestMatrix_lab_runs.csv`. |
| **`plot/lab/1_Monopile_matrix.xlsx`** | Lab schedule workbook (**Run Log** = wall-clock and Trial IDs). |
| **`plot/lab/mat_run_map.json`** | **Orphans only:** mats without a dump, pending uploads, duplicate mat aliases. Not paired runs. |
| **`plot/lab/STATEOS_SIGNALS.md`** | Seki §2.1 / `typeConv3` field guide for `hist_os_state.png`. |

**Archive:** `G:\Shared drives\Simpson team\Test Data\2026-OSU-SSI-Bridge\`  
(`2026-08-19/opensees_data/`, `2026-08-21/opensees_data/`)  
**Working mirror:** `OSU_SSI_BRIDGE_DATA_LOCAL/opensees_data/`  
**Plots:** `OSU_SSI_BRIDGE_DATA_LOCAL/plots/` — `runs/<Test>/{eq,os}/`, `compare/…`

## Test ID scheme

| Prefix | Meaning |
|--------|---------|
| **W##** | Wednesday 2026-08-19 bridge RTHS Trial (Run Log). Knobs ≈ `2026-08-19` tag + **constraints Auto** + **ParallelProfileSPD** (same family as Fri `r+01`). |
| **F##** | Friday 2026-08-21 bridge RTHS Trial (gaps OK: no F02, F16, …). |
| **Fd##** | Friday lunch dry (no OpenSees dump). |
| **Fx##** | Friday dump without a clean Trial number (abort / early stop). |

Wed (4 pm+ bridge Trials 1–7): **W01–W07**. Mats + OpenSees dumps live under
`GustavoModel/data/08192026` (copied to LOCAL and Simpson archive). Filenames are
all `H0p5_T7p746` (Run Log `H0p156` labels were outdated). Orphan probe:
`0819_testGus.mat`.

## Loading

**Friday (2026-08-21):** earthquake followed by wave (folder `Storm_Wave` is a
matrix nickname only). **Wednesday (2026-08-19):** EQ together with deep-water
waves. Infer loading from the dump and lab notes, not from the folder name alone.

## How mat ↔ dump pairs were chosen

1. Lab trial log (spreadsheet) — wall-clock, trial number, dry/NAN flags.
2. Row tag — `rowNeg9` → `r-09_…`.
3. Duration — mat `Time_last × √2.4` ≈ OpenSees pier record length (ratio ≈ 1.0).
4. Early shape — pier Δux vs `meaSigOS` when clocks desync.

**Clocks:** mat `Time` = lab DAQ (model scale). OpenSees `t` = prototype. λ = 2.4.

**Slowdowns:** `stateOS` / `typeConv3` = 2 marks OpenFresco waits. Many brief events
can still give duration ratio ≈ 1.0 if sample fraction in state 2 stays small.
Campaign bars: `python plot/PlotStateOSBars.py` → `plots/compare/stateos/`.

## Compare reference

Eventual true reference = offline OpenSees per folder (`realTimeON 0`, no OpenFresco).
Until then, pairwise overlays (`PlotEQComparePairs.py`) use the fewest `typeConv3→2`
rising edges in GM **D5–95** as interim reference. Skip single-precision CuDSS and
twoNodeLink. Keep dry OpenFresco baseline **F01**.

## Pending / incomplete (see CSV `Note` + orphan JSON)

- **`pending_mat_upload`:** cleared for `rowNeg9_24Core` (paired as **F11** → `r-09_1132`).
- **Wed mats + dumps:** ingested from Gustavo `…/08192026` → LOCAL + Simpson archive.
- **Dumps without mat** (empty `MatFile`): `Fx##` aborts + W01/W02 pending mats — notes in CSV.
- **Mats without dump** (`mats_without_dump` / **Fd##**): lunch np sweeps.

**Compare layout:** `plots/compare/{Baseline,Moderate,Large,X-Large}/<variant>/pairs/`
with `hist_ux_pair_F##.png`. Excluded (twoNodeLink / single-precision)
→ `plots/compare/_excluded/<reason>/…`. Per-run plots: `plots/runs/<Test>/{eq,os}/`.

**Dry OpenFresco baseline (F01, `r+01_…0836`):** included in campaign bars/pairs.
**Dry (Fd01–Fd09):** no dump; stateOS bars in `bar_typeconv3_*_with_dry.png` only.

## Scripts

- Ingest + mat extract: `python plot/SyncLabBackup.py`
- OpenSees EQ plots: `python plot/PlotEQParallel.py <LOCAL dump>`
- Simulink *OS: `python plot/PlotMatOS.py`
- Pairwise pier compare: `python plot/PlotEQComparePairs.py`
