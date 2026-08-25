# OSU SSI Bridge — lab run map (2026-08-21)

Maps each Simulink `.mat` to an OpenSees dump folder (or says why it has no partner).
Companion to `mat_run_map.json` (used by the plot scripts).

**Archive:** `G:\Shared drives\Simpson team\Test Data\2026-08-21-OSU-SSI-Bridge\opensees_data\`  
**Machine mirror:** `OSU_SSI_BRIDGE_DATA_LOCAL/opensees_data/` (gitignored)  
**Plots:** `OSU_SSI_BRIDGE_DATA_LOCAL/plots/`

---

## How pairs were chosen

1. **Lab trial log** (spreadsheet, Aug 21) — wall-clock HHMM in notes, trial number, dry/NAN flags.
2. **Row tag** — `rowNeg9` → `r-09_…`, etc.
3. **Duration** — mat `Time_last × √2.4` ≈ OpenSees `pier_top_disp` last `t` (ratio ≈ 1.0).
4. **Early shape** — pier Δux vs `meaSigOS` on first ~50 s prototype time when slowdowns desync clocks.

**Clocks:** mat `Time` = lab DAQ (model scale). OpenSees `t` = numerical prototype time.  
Froude scale λ = 2.4 → multiply lab time by √λ for prototype overlays.

**Slowdowns:** `stateOS` column 0 (`typeConv3`) = 2 marks OpenFresco waits. Many brief events
(100–500+) can still give duration ratio ≈ 1.0 if wall time in state 2 stays below ~1%.

---

## Confirmed pairs (mat ↔ dump)

| Trial | Lab time | Mat | OpenSees dump | np | Slowdowns / holds (lab log) | Status |
|------:|----------|-----|---------------|---:|------------------------------|--------|
| 01 | 08:36 | `0821_testBaseline_Run02.mat` | `r+01_20260821_0836_Storm_Wave` | 8 | PID dry baseline | paired |
| 03 | 09:30 | `0821_GusBridge_rowPos1.mat` | `r+01_20260821_0930_Storm_Wave` | 8 | 1 slowdown | paired |
| 04 | 09:49 | `0821_GusBridge_rowPos1_16Core.mat` | `r+01_20260821_0949_Storm_Wave` | 16 | 5 slowdowns | paired |
| 05 | 10:10 | `0821_GusBridge_rowNeg4.mat` | `r-04_20260821_1010_Storm_Wave` | 8 | 0 slowdowns | paired |
| 06 | 10:34 | `0821_GusBridge_rowNeg8.mat` | `r-08_20260821_1034_Storm_Wave` | 8 | — | paired |
| 07 | 10:54 | `0821_GusBridge_rowNeg80.mat` | `r-80_20260821_1054_Storm_Wave` | 8 | all slowdown | paired, incomplete OS |
| 08 | 10:58 | `0821_GusBridge_rowNeg81.mat` | `r-81_20260821_1058_Storm_Wave` | 8 | heavy slowdown | paired, incomplete OS |
| 09 | 11:02 | `0821_GusBridge_rowNeg9.mat` | `r-09_20260821_1102_Storm_Wave` | 8 | ≥4 holds | paired |
| 10 | 11:08 | `0821_GusBridge_rowNeg9_16Core.mat` | `r-09_20260821_1108_Storm_Wave` | 16 | ignore first 20 s | paired |
| 11 | 11:32 | `0821_GusBridge_rowNeg9_24Core.mat` | `r-09_20260821_1132_Storm_Wave` | 24 | 1 hold | **mat pending upload** |
| — | 13:03 | `0821_GusBridge_rowNeg9_20Core.mat` | `r-09_20260821_1303_Storm_Wave` | 20 | ~5 slowdowns | paired |
| — | 13:12 | `0821_GusBridge_rowNeg10_24Core_duplicate.mat` | `r-10_20260821_1312_Storm_Wave` | 24 | heavy slowdowns | paired |
| — | 13:42 | `0821_GusBridge_rowNeg91_16Core.mat` | `r-91_20260821_1342_Storm_Wave` | 16 | heavy slowdowns | paired |
| — | 13:58 | `0821_GusBridge_rowNeg91_20Core.mat` | `r-91_20260821_1358_Storm_Wave` | 20 | heavy slowdowns | paired |
| 17 | 14:51 | `0821_GusBridge_rowNeg11_8Core.mat` | `r-11_20260821_1451_Storm_Wave` | 8 | 2 slowdowns | paired |
| 20 | 15:21 | `0821_GusBridge_rowNeg15_8Core.mat` | `r-15_20260821_1521_Storm_Wave` | 8 | 3 holds, 9 slowdowns | paired |
| 22 | 15:28 | `0821_GusBridge_rowNeg6_Serial.mat` | `r-06_20260821_1528_Storm_Wave` | 1 | ≥35 holds; ~2/3 slowdown; ~80 s lab | paired, incomplete OS |
| 23 | 15:35 | `0821_GusBridge_rowNeg13_8Core.mat` | `r-13_20260821_1535_Storm_Wave` | 8 | no slowdowns | paired |
| 24 | 15:44 | `0821_GusBridge_rowNeg14_8Core.mat` | `r-14_20260821_1544_Storm_Wave` | 8 | 1 slowdown @ 10 s | paired |
| — | 15:12 | `0821_GusBridge_rowNeg120_20Core.mat` | `r-120_20260821_1512_Storm_Wave` | 20 | retry after 1503 NAN | paired, very short OS |
| — | 15:16 | `0821_GusBridge_rowNeg121_20Core.mat` | `r-121_20260821_1516_Storm_Wave` | 20 | retry after 1503 NAN | paired, very short OS |
| 26 | 16:01 | `0821_GusBridge_rowNeg2_8Core.mat` | `r-02_20260821_1601_Storm_Wave` | 8 | ~10 slowdowns, 1–2 holds | paired |
| 27 | 16:25 | `0821_GusBridge_rowNeg4_8Core_P1.mat` | `r-04_20260821_1625_Storm_Wave` | 8 | 1 slowdown | paired |
| 28 | 16:39 | `0821_GusBridge_rowNeg16_8Core_P1.mat` | `r-16_20260821_1639_Storm_Wave` | 8 | crazy signal; 1 slowdown @ 25 s | paired, incomplete OS |
| 29 | 16:44 | `0821_GusBridge_rowNeg16_8Core_P1_alpha0.mat` | `r-16_20260821_1644_Storm_Wave` | 8 | CudaMKRAlpha 0.0 | paired |
| 29 | 16:44 | `0821_GusBridge_rowNeg16_8Core_P1_alpha0_dry.mat` | `r-16_20260821_1644_Storm_Wave` | 8 | same session as row above | duplicate mat |
| 30 | 16:49 | `0821_GusBridge_rowNeg17_8Core_P1.mat` | `r-17_20260821_1649_Storm_Wave` | 8 | no slowdown | paired |

Trial numbers without a row in the morning/afternoon sheets (np sweeps, extra sessions) are marked “—”.

---

## Pending upload

| Item | Expected partner | Notes |
|------|------------------|-------|
| `0821_GusBridge_rowNeg9_24Core.mat` | `r-09_20260821_1132_Storm_Wave` | Trial 11; dump on Drive, mat not in Simulink yet |
| `0821_GusBridge_rowNeg8_16Core.mat` | new dump ~12:45 wall clock | 16-core `-08`; folder not on Drive yet |

---

## OpenSees dumps without mat

| Dump | t_last (OS) | Explanation |
|------|------------:|-------------|
| `r-04_20260821_1515_Storm_Wave` | ~2.6 s | TestMatrix row; abort / dry start between afternoon `-04` attempts. No Simulink save. |
| `r-04_20260821_1520_Storm_Wave` | ~2.6 s | Same — second abort before `-15` / `-06` block. |
| `r-04_20260821_1615_Storm_Wave` | ~194 s | Incomplete `-04` run between Trial 26 (`1601`) and Trial 27 (`1625`). No mat on Drive; likely discarded dry / failed hybrid. **Not** `rowNeg2_16Core`. |
| `r-12_20260821_1503_Storm_Wave` | ~9.6 s | Trial 18: `rowNeg12_20Core` — “didn’t work, NAN”. Mat not archived. |
| `r-16_20260821_1622_Storm_Wave` | ~156 s | `-16` twoNodeLink, no hold. Pier ux matches `1639` early (same physics); no full-length mat. Likely restart between P1 trials; mat not saved or discarded. |

---

## Simulink mats without dump

| Mat | t_lab (s) | Explanation |
|-----|----------:|-------------|
| `0821_testBaseline_Run01.mat` | ~107 | Short baseline; no OpenSees dump. |
| `0821_GusBridge_rowNeg2_16Core.mat` | ~121 | 16-core `-02` probe; only `r-02_1601` (8-core) on Drive. Dry / np check or incomplete. |
| `0821_GusBridge_rowNeg8_16Core.mat` | ~201 | 16-core `-08`; morning dump is 8-core `1034`. Awaiting ~1245 folder. |
| `0821_GusBridge_rowNeg9_18Core.mat`, `_19Core`, `_22Core` | ~202–207 | `-09` np sweeps; no dump (duration ratio ~0.73 vs full runs). Dry trials. |
| `0821_GusBridge_rowNeg10_20Core.mat` … `_28Core` | ~200–214 | `-10` np sweeps; no dump. Use `Neg10_24Core_duplicate` for `r-10_1312`. |

---

## Failed / NAN trials (lab log)

| Trial | Mat attempted | Dump | Outcome |
|------:|---------------|------|---------|
| 18 | `rowNeg12_20Core` | `r-12_1503` | NAN; dump only (~9.6 s OS) |
| 19 | `rowNeg120_20Core` | (1503 slot) | NAN at 1503; **retried** → `r-120_1512` + mat |
| — | `rowNeg121_20Core` | (1503 slot) | NAN note; **retried** → `r-121_1516` + mat |

---

## Slowdown reference (from lab log + extracts)

Heavy event count but duration ratio still ≈ 1.0:

- `rowNeg91_16Core` (~149 events), `rowNeg91_20Core` (~568)
- `rowNeg10_24Core_duplicate` (~97)
- `rowNeg9_20Core` (~5–7)

Clock desync matters for **realtime** overlays (`hist_ux_*_realtime.png`), not for rejecting
these pairs. Use OpenSees-time plots or early-window correlation when unsure.

---

## Overwrites

Run folders use unique `r±NN_YYYYMMDD_HHMM_` names. Coexisting variants (`alpha0` + `alpha0_dry`,
`Neg10_24Core` + `_duplicate`, multiple `-04` / `-16` timestamps) argue against silent overwrites.
Gaps look like **dry trials or NAN aborts** where one side was kept.

---

## Maintenance

- Machine map: `plot/opensees_data/mat_run_map.json` (used by `PlotMatOS.py`, `PlotEQCompareRuns.py`).
- After Drive adds files: `python plot/SyncLabBackup.py`, then refresh plots.
- Update this file when new pairs are confirmed or uploads land.
