# Soil profiles

Units in the model: **N, m, s**. Tables also list ft / kPa where that matches pile `dy` or the OpenSees wiki.

Coordinate: **z = 0 at cap top** (same as the structure). +z up. Cap occupies `0 → −H_cap` (`H_cap = 3.25 ft`). Pile head at `z = −H_cap`, tip at `z = −H_cap − 60 ft = −63.25 ft`. Pile / soil vertical step: **`dy = 3 ft = 0.9144 m`** (`nSeg_pile = 20`).

Layer interfaces sit on **pile nodes**: `z = −3.25 − k·3` ft.

Three named units on every option: **L2 / L3 / L5**. What changes is the L2 thickness (opt. 4 only) and the L3/L5 material.

| Option | L2 | L3 | L5 |
|---|---|---|---|
| **1** | crust clay, 9.25 ft | liquefiable sand, one PDMY02 \(D_r\) ramp | dense sand (PDMY02 + FSP), constant |
| **2** | crust clay, 9.25 ft | medium stiff clay, PIMY \(c\)-\(G_r\) ramp | dense sand (PDMY02 + FSP), constant |
| **3** | crust clay, 9.25 ft | medium stiff clay (same as opt. 2) | **stiff clay**, wiki PIMY ramp |
| **4** | **soft clay, 27.25 ft** (same top; deeper L2/L3 contact) | medium stiff clay, shorter remainder | stiff clay (same as opt. 3) |

PIMY on all clays (**no** `FluidSolidPorous`). PDMY02 + FSP on sands only.

Layer thicknesses follow the Shin et al. (2007) / Kramer PEER 2008/07 liquefiable-bridge class (prototype P3-style column), snapped to the pile grid. Constitutive **ends** mix Kramer Table 5.1 clay (\(c\), \(G_r\), \(\rho\)) with OpenSees wiki PDMY02 columns vs \(D_r\). Within each named unit, props vary linearly in \(z\) (clay-on-clay contacts match).

`BuildSoilMaterials.tcl` builds one nDMaterial per 3 ft row (FSP wrap on sands). Mesh and springs read `soilMatRow($iy)` / row `c`, \(\phi\), \(\rho\).

---

## Mesh height

Soil quads follow the pile stations. Uniform **3 ft** quads (same as pile `dy`): one soil row per pile segment, plus two thinner rows through the cap. Near-surface PEER meshes often use ~2–2.5 ft; dense sand along the pile ~3.3 ft; below tip sometimes coarser. Three feet sits in that range.

---

## Geometry (all options)

Prototype thicknesses from grade, rounded onto the pile grid. Grade ≈ cap top.

| Unit | Prototype (m / ft) | Rounded | z top (ft) | z bot (ft) | z top (m) | z bot (m) | Pile segs |
|---|---|---|---|---|---|---|---|
| L2 opts 1–3 | 2.43 / 8.0 | 9.25 ft from z=0 | 0 | −9.25 | 0 | −2.819 | 2 below cap |
| L2 opt. 4 | thicker crust (same grade) | 27.25 ft | 0 | −27.25 | 0 | −8.306 | 8 below cap |
| L3 opts 1–3 | 9.73 / 31.9 | 33 ft | −9.25 | −42.25 | −2.819 | −12.878 | 11 |
| L3 opt. 4 | remainder under thicker L2 | 15 ft | −27.25 | −42.25 | −8.306 | −12.878 | 5 |
| L5 (on pile) | ~6.1 / 20 to tip | 21 ft | −42.25 | −63.25 | −12.878 | −19.279 | 7 |
| L5 below tip | ~4.6 / 15 | 15 ft | −63.25 | −78.25 | −19.279 | −23.851 | 5 extra |

L5 top is always −42.25 ft. Opt. 4 only moves the **L2/L3** contact down 18 ft. Same quad mesh.

Cap occupies 3.25 ft of L2 in the structure model. The continuum still runs through that footprint (same overlap as the pile shafts). L2 along the pile below the soffit is **6 ft** (opts 1–3) or **24 ft** (opt. 4). Face springs sit on the soil line at `|x| = s` (outer pile axes); soffit q-z at `z = −H_cap`.

27 soil rows (2 cap + 25 × 3 ft). One material per row, evaluated at the **row centroid**.

---

## In-layer interpolation

Linear in \(z\) between unit **control points**. \(\xi = (z_\mathrm{top} - z)/(z_\mathrm{top} - z_\mathrm{bot})\). \(\xi=0\) at the top of the unit, \(1\) at the bottom. Clay: \(\rho\), \(G_r\), \(c\) share that \(\xi\); \(B_r = 50\,G_r\) (undrained \(\nu \approx 0.49\)). Sand: \(D_r\) linear in \(z\), then the wiki PDMY02 columns linear in \(D_r\) between 50–60% and 60–75%.

Shin et al. (2007) Table 1 ranges, Kramer (2008) Table 5.1 clay ends used here:

| Knot | Source | \(\rho\) (kg/m³) | \(G_r\) (MPa) | \(c\) (kPa) | \(\phi\) (°) | \(D_r\) |
|---|---|---|---|---|---|---|
| crust clay, low | Kramer T5.1 / Shin crust | 1488 | 57.28 | 35.9 | 0 | — |
| medium clay, low | Kramer T5.1 deep-clay top | 1521 | 66.39 | 39.7 | 0 | — |
| medium clay, high | Kramer T5.1 deep-clay base | 1669 | 110.68 | 58.4 | 0 | — |
| wiki soft clay | OpenSees PIMY | 1300 | 13.0 | 18.0 | 0 | — |
| wiki stiff clay | OpenSees PIMY | 1800 | 150.0 | 75.0 | 0 | — |
| wiki sand 50% | OpenSees PDMY02 | 1900 | 100 | 0 | 33.5 | 50% |
| wiki sand 60% | OpenSees PDMY02 | 2000 | 110 | 0 | 35.0 | 60% |
| wiki sand 75% | OpenSees PDMY02 | 2100 | 130 | 0 | 36.5 | 75% |
| PEER dense sand | L5 sand | 2260 | 154.8 | 0 | 39.3 | ~95% |

Clay-on-clay contacts **meet**: L2 bot = L3 top (\(c=39.7\) kPa on opts 2–4); L3 bot = L5 top on opts 3–4 (\(c=58.4\)). Opt. 1 L2 uses the full Shin crust 36–58 kPa (L3 below is sand). Kramer’s crust in Table 5.1 dips in the middle of a thicker surface clay; we keep a monotonic ramp.

L5 dense sand is flat (Shin \(\phi=40^\circ\); PEER \(G_r\) already ~155 MPa). Stiff and soft clay ends are wiki, not Shin Table 1.

Centroid samples never sit on the interface, so adjacent rows still differ by one half-step of the local slope. Sand/clay contacts also change the constitutive model.

---

## Option 1 — L3 = liquefiable sand

L2: crust \(c=35.9\to58.4\) kPa. L3: one sand ramp \(D_r=50\to75\%\) over 33 ft (Shin \(\phi=33\sim36^\circ\)). L5: PEER dense sand, constant.

Control points (interfaces):

| z (ft) | Face | \(\rho\) | \(G_r\) (MPa) | \(c\) (kPa) | \(\phi\) (°) | \(D_r\) (%) |
|---|---|---:|---:|---:|---:|---:|
| 0 | L2 top | 1488 | 57.3 | 35.9 | 0 | — |
| −9.25 | L2 bot / L3 top | 1669 / 1900 | 110.7 / 100 | 58.4 / 0 | 0 / 33.5 | — / 50 |
| −42.25 | L3 bot / L5 top | 2100 / 2260 | 130 / 154.8 | 0 | 36.5 / 39.3 | 75 / 95 |
| −78.25 | L5 bot | 2260 | 154.8 | 0 | 39.3 | 95 |

Per-row (centroid). \(G_r\) ratio at L3/L5 ≈ **1.21**. L2/L3 \(G_r\) almost 1.0 (clay slightly stiffer than the first sand row).

| # | z top | z bot | z mid | unit | type | ρ (kg/m³) | G_r (MPa) | B_r (MPa) | c (kPa) | φ (°) | D_r (%) |
|---:|---:|---:|---:|---|---|---:|---:|---:|---:|---:|---:|
| 1 | 0 | −1.625 | −0.8125 | L2 | clay | 1504 | 62.0 | 3099 | 37.9 | 0 | — |
| 2 | −1.625 | −3.25 | −2.4375 | L2 | clay | 1536 | 71.4 | 3568 | 41.8 | 0 | — |
| 3 | −3.25 | −6.25 | −4.75 | L2 | clay | 1581 | 84.7 | 4235 | 47.5 | 0 | — |
| 4 | −6.25 | −9.25 | −7.75 | L2 | clay | 1640 | 102.0 | 5101 | 54.8 | 0 | — |
| 5 | −9.25 | −12.25 | −10.75 | L3 | sand | 1911 | 101.1 | 234 | 0 | 33.7 | 51.1 |
| 6 | −12.25 | −15.25 | −13.75 | L3 | sand | 1934 | 103.4 | 235 | 0 | 34.0 | 53.4 |
| 7 | −15.25 | −18.25 | −16.75 | L3 | sand | 1957 | 105.7 | 237 | 0 | 34.4 | 55.7 |
| 8 | −18.25 | −21.25 | −19.75 | L3 | sand | 1980 | 108.0 | 239 | 0 | 34.7 | 58.0 |
| 9 | −21.25 | −24.25 | −22.75 | L3 | sand | 2002 | 110.3 | 240 | 0 | 35.0 | 60.2 |
| 10 | −24.25 | −27.25 | −25.75 | L3 | sand | 2017 | 113.3 | 243 | 0 | 35.2 | 62.5 |
| 11 | −27.25 | −30.25 | −28.75 | L3 | sand | 2032 | 116.4 | 246 | 0 | 35.5 | 64.8 |
| 12 | −30.25 | −33.25 | −31.75 | L3 | sand | 2047 | 119.4 | 249 | 0 | 35.7 | 67.0 |
| 13 | −33.25 | −36.25 | −34.75 | L3 | sand | 2062 | 122.4 | 252 | 0 | 35.9 | 69.3 |
| 14 | −36.25 | −39.25 | −37.75 | L3 | sand | 2077 | 125.5 | 255 | 0 | 36.2 | 71.6 |
| 15 | −39.25 | −42.25 | −40.75 | L3 | sand | 2092 | 128.5 | 258 | 0 | 36.4 | 73.9 |
| 16 | −42.25 | −45.25 | −43.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |
| 17 | −45.25 | −48.25 | −46.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |
| 18 | −48.25 | −51.25 | −49.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |
| 19 | −51.25 | −54.25 | −52.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |
| 20 | −54.25 | −57.25 | −55.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |
| 21 | −57.25 | −60.25 | −58.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |
| 22 | −60.25 | −63.25 | −61.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |
| 23 | −63.25 | −66.25 | −64.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |
| 24 | −66.25 | −69.25 | −67.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |
| 25 | −69.25 | −72.25 | −70.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |
| 26 | −72.25 | −75.25 | −73.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |
| 27 | −75.25 | −78.25 | −76.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |

Sand \(B_r\) is the wiki drained-ish bulk (not \(50 G_r\)). API \(k\) for p-y follows \(D_r\) between 60, 81, and 121 pci.

---

## Option 2 — L3 = medium stiff clay

Same L5 as Option 1. L2 meets L3 at \(c=39.7\) kPa (Shin/Kramer medium-clay low), so the crust only uses 36–40 kPa of the 36–58 crust range. L3: \(c=39.7\to58.4\) kPa over 33 ft.

Control points:

| z (ft) | Face | \(\rho\) | \(G_r\) (MPa) | \(c\) (kPa) | \(\phi\) (°) | \(D_r\) (%) |
|---|---|---:|---:|---:|---:|---:|
| 0 | L2 top | 1488 | 57.3 | 35.9 | 0 | — |
| −9.25 | L2 bot / L3 top | 1521 | 66.4 | 39.7 | 0 | — |
| −42.25 | L3 bot / L5 top | 1669 / 2260 | 110.7 / 154.8 | 58.4 / 0 | 0 / 39.3 | — / 95 |
| −78.25 | L5 bot | 2260 | 154.8 | 0 | 39.3 | 95 |

L3/L5 \(G_r\) ratio ≈ **1.42** (plus PIMY→PDMY). L2/L3 ≈ **1.05**.

| # | z top | z bot | z mid | unit | type | ρ (kg/m³) | G_r (MPa) | B_r (MPa) | c (kPa) | φ (°) | D_r (%) |
|---:|---:|---:|---:|---|---|---:|---:|---:|---:|---:|---:|
| 1 | 0 | −1.625 | −0.8125 | L2 | clay | 1491 | 58.1 | 2904 | 36.2 | 0 | — |
| 2 | −1.625 | −3.25 | −2.4375 | L2 | clay | 1497 | 59.7 | 2984 | 36.9 | 0 | — |
| 3 | −3.25 | −6.25 | −4.75 | L2 | clay | 1505 | 62.0 | 3098 | 37.9 | 0 | — |
| 4 | −6.25 | −9.25 | −7.75 | L2 | clay | 1516 | 64.9 | 3246 | 39.1 | 0 | — |
| 5 | −9.25 | −12.25 | −10.75 | L3 | clay | 1528 | 68.4 | 3420 | 40.6 | 0 | — |
| 6 | −12.25 | −15.25 | −13.75 | L3 | clay | 1541 | 72.4 | 3622 | 42.2 | 0 | — |
| 7 | −15.25 | −18.25 | −16.75 | L3 | clay | 1555 | 76.5 | 3823 | 44.0 | 0 | — |
| 8 | −18.25 | −21.25 | −19.75 | L3 | clay | 1568 | 80.5 | 4024 | 45.6 | 0 | — |
| 9 | −21.25 | −24.25 | −22.75 | L3 | clay | 1582 | 84.5 | 4226 | 47.4 | 0 | — |
| 10 | −24.25 | −27.25 | −25.75 | L3 | clay | 1595 | 88.5 | 4427 | 49.0 | 0 | — |
| 11 | −27.25 | −30.25 | −28.75 | L3 | clay | 1608 | 92.6 | 4628 | 50.8 | 0 | — |
| 12 | −30.25 | −33.25 | −31.75 | L3 | clay | 1622 | 96.6 | 4830 | 52.5 | 0 | — |
| 13 | −33.25 | −36.25 | −34.75 | L3 | clay | 1635 | 100.6 | 5031 | 54.1 | 0 | — |
| 14 | −36.25 | −39.25 | −37.75 | L3 | clay | 1649 | 104.6 | 5232 | 55.8 | 0 | — |
| 15 | −39.25 | −42.25 | −40.75 | L3 | clay | 1662 | 108.7 | 5434 | 57.5 | 0 | — |
| 16 | −42.25 | −45.25 | −43.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |
| 17 | −45.25 | −48.25 | −46.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |
| 18 | −48.25 | −51.25 | −49.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |
| 19 | −51.25 | −54.25 | −52.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |
| 20 | −54.25 | −57.25 | −55.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |
| 21 | −57.25 | −60.25 | −58.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |
| 22 | −60.25 | −63.25 | −61.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |
| 23 | −63.25 | −66.25 | −64.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |
| 24 | −66.25 | −69.25 | −67.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |
| 25 | −69.25 | −72.25 | −70.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |
| 26 | −72.25 | −75.25 | −73.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |
| 27 | −75.25 | −78.25 | −76.75 | L5 | sand | 2260 | 154.8 | 279 | 0 | 39.3 | 95.0 |

This is the non-liquefiable bound on the same mesh (not a full Pier-1 column with sand between L2 and bottom clay).

---

## Option 3 — L3 medium clay + L5 stiff clay

Same L2 / L3 as Option 2. L5 is wiki [PIMY stiff clay](https://opensees.berkeley.edu/wiki/index.php/PressureIndependMultiYield_Material) instead of dense sand. Ramp meets L3 at \(c=58.4\) kPa and reaches wiki \(c=75\) kPa at the base. No sand in the column.

Control points:

| z (ft) | Face | \(\rho\) | \(G_r\) (MPa) | \(c\) (kPa) |
|---|---|---:|---:|---:|
| 0 | L2 top | 1488 | 57.3 | 35.9 |
| −9.25 | L2 bot / L3 top | 1521 | 66.4 | 39.7 |
| −42.25 | L3 bot / L5 top | 1669 | 110.7 | 58.4 |
| −78.25 | L5 bot | 1800 | 150.0 | 75.0 |

L3/L5 \(G_r\) ratio ≈ **1.03** (was 2.3 with two constants). Same row 1–15 as Option 2.

| # | z top | z bot | z mid | unit | type | ρ (kg/m³) | G_r (MPa) | B_r (MPa) | c (kPa) | φ (°) | D_r (%) |
|---:|---:|---:|---:|---|---|---:|---:|---:|---:|---:|---:|
| 1 | 0 | −1.625 | −0.8125 | L2 | clay | 1491 | 58.1 | 2904 | 36.2 | 0 | — |
| 2 | −1.625 | −3.25 | −2.4375 | L2 | clay | 1497 | 59.7 | 2984 | 36.9 | 0 | — |
| 3 | −3.25 | −6.25 | −4.75 | L2 | clay | 1505 | 62.0 | 3098 | 37.9 | 0 | — |
| 4 | −6.25 | −9.25 | −7.75 | L2 | clay | 1516 | 64.9 | 3246 | 39.1 | 0 | — |
| 5 | −9.25 | −12.25 | −10.75 | L3 | clay | 1528 | 68.4 | 3420 | 40.6 | 0 | — |
| 6 | −12.25 | −15.25 | −13.75 | L3 | clay | 1541 | 72.4 | 3622 | 42.2 | 0 | — |
| 7 | −15.25 | −18.25 | −16.75 | L3 | clay | 1555 | 76.5 | 3823 | 44.0 | 0 | — |
| 8 | −18.25 | −21.25 | −19.75 | L3 | clay | 1568 | 80.5 | 4024 | 45.6 | 0 | — |
| 9 | −21.25 | −24.25 | −22.75 | L3 | clay | 1582 | 84.5 | 4226 | 47.4 | 0 | — |
| 10 | −24.25 | −27.25 | −25.75 | L3 | clay | 1595 | 88.5 | 4427 | 49.0 | 0 | — |
| 11 | −27.25 | −30.25 | −28.75 | L3 | clay | 1608 | 92.6 | 4628 | 50.8 | 0 | — |
| 12 | −30.25 | −33.25 | −31.75 | L3 | clay | 1622 | 96.6 | 4830 | 52.5 | 0 | — |
| 13 | −33.25 | −36.25 | −34.75 | L3 | clay | 1635 | 100.6 | 5031 | 54.1 | 0 | — |
| 14 | −36.25 | −39.25 | −37.75 | L3 | clay | 1649 | 104.6 | 5232 | 55.8 | 0 | — |
| 15 | −39.25 | −42.25 | −40.75 | L3 | clay | 1662 | 108.7 | 5434 | 57.5 | 0 | — |
| 16 | −42.25 | −45.25 | −43.75 | L5 | clay | 1674 | 112.3 | 5616 | 59.1 | 0 | — |
| 17 | −45.25 | −48.25 | −46.75 | L5 | clay | 1685 | 115.6 | 5780 | 60.5 | 0 | — |
| 18 | −48.25 | −51.25 | −49.75 | L5 | clay | 1696 | 118.9 | 5944 | 61.9 | 0 | — |
| 19 | −51.25 | −54.25 | −52.75 | L5 | clay | 1707 | 122.2 | 6108 | 63.2 | 0 | — |
| 20 | −54.25 | −57.25 | −55.75 | L5 | clay | 1718 | 125.4 | 6271 | 64.6 | 0 | — |
| 21 | −57.25 | −60.25 | −58.75 | L5 | clay | 1729 | 128.7 | 6435 | 66.0 | 0 | — |
| 22 | −60.25 | −63.25 | −61.75 | L5 | clay | 1740 | 132.0 | 6599 | 67.4 | 0 | — |
| 23 | −63.25 | −66.25 | −64.75 | L5 | clay | 1751 | 135.3 | 6763 | 68.8 | 0 | — |
| 24 | −66.25 | −69.25 | −67.75 | L5 | clay | 1762 | 138.5 | 6927 | 70.2 | 0 | — |
| 25 | −69.25 | −72.25 | −70.75 | L5 | clay | 1773 | 141.8 | 7090 | 71.5 | 0 | — |
| 26 | −72.25 | −75.25 | −73.75 | L5 | clay | 1784 | 145.1 | 7254 | 72.9 | 0 | — |
| 27 | −75.25 | −78.25 | −76.75 | L5 | clay | 1795 | 148.4 | 7418 | 74.3 | 0 | — |

---

## Option 4 — thicker soft L2, medium L3, stiff L5

Same L5 as Option 3. L2 is the same crust contact moved down 18 ft: soft clay from grade to **−27.25 ft** (~14 D for D = 2 ft), wiki [PIMY soft](https://opensees.berkeley.edu/wiki/index.php/PressureIndependMultiYield_Material) at the top meeting medium-clay low at the L2/L3 contact. L3 is then the 15 ft remainder (\(c=39.7\to58.4\) kPa). No extra unit name.

Control points:

| z (ft) | Face | \(\rho\) | \(G_r\) (MPa) | \(c\) (kPa) |
|---|---|---:|---:|---:|
| 0 | L2 top (soft) | 1300 | 13.0 | 18.0 |
| −27.25 | L2 bot / L3 top | 1521 | 66.4 | 39.7 |
| −42.25 | L3 bot / L5 top | 1669 | 110.7 | 58.4 |
| −78.25 | L5 bot | 1800 | 150.0 | 75.0 |

L2/L3 \(G_r\) ratio at the −27.25 ft centroids ≈ **1.12** (was 5.1). L3/L5 ≈ **1.06**.

| # | z top | z bot | z mid | unit | type | ρ (kg/m³) | G_r (MPa) | B_r (MPa) | c (kPa) | φ (°) | D_r (%) |
|---:|---:|---:|---:|---|---|---:|---:|---:|---:|---:|---:|
| 1 | 0 | −1.625 | −0.8125 | L2 | clay | 1307 | 14.6 | 730 | 18.6 | 0 | — |
| 2 | −1.625 | −3.25 | −2.4375 | L2 | clay | 1320 | 17.8 | 889 | 19.9 | 0 | — |
| 3 | −3.25 | −6.25 | −4.75 | L2 | clay | 1339 | 22.3 | 1115 | 21.8 | 0 | — |
| 4 | −6.25 | −9.25 | −7.75 | L2 | clay | 1363 | 28.2 | 1409 | 24.2 | 0 | — |
| 5 | −9.25 | −12.25 | −10.75 | L2 | clay | 1387 | 34.1 | 1703 | 26.6 | 0 | — |
| 6 | −12.25 | −15.25 | −13.75 | L2 | clay | 1412 | 39.9 | 1997 | 28.9 | 0 | — |
| 7 | −15.25 | −18.25 | −16.75 | L2 | clay | 1436 | 45.8 | 2291 | 31.3 | 0 | — |
| 8 | −18.25 | −21.25 | −19.75 | L2 | clay | 1460 | 51.7 | 2585 | 33.7 | 0 | — |
| 9 | −21.25 | −24.25 | −22.75 | L2 | clay | 1485 | 57.6 | 2879 | 36.1 | 0 | — |
| 10 | −24.25 | −27.25 | −25.75 | L2 | clay | 1509 | 63.5 | 3173 | 38.5 | 0 | — |
| 11 | −27.25 | −30.25 | −28.75 | L3 | clay | 1536 | 70.8 | 3541 | 41.6 | 0 | — |
| 12 | −30.25 | −33.25 | −31.75 | L3 | clay | 1565 | 79.7 | 3984 | 45.3 | 0 | — |
| 13 | −33.25 | −36.25 | −34.75 | L3 | clay | 1595 | 88.5 | 4427 | 49.0 | 0 | — |
| 14 | −36.25 | −39.25 | −37.75 | L3 | clay | 1625 | 97.4 | 4870 | 52.8 | 0 | — |
| 15 | −39.25 | −42.25 | −40.75 | L3 | clay | 1654 | 106.3 | 5313 | 56.5 | 0 | — |
| 16 | −42.25 | −45.25 | −43.75 | L5 | clay | 1674 | 112.3 | 5616 | 59.1 | 0 | — |
| 17 | −45.25 | −48.25 | −46.75 | L5 | clay | 1685 | 115.6 | 5780 | 60.5 | 0 | — |
| 18 | −48.25 | −51.25 | −49.75 | L5 | clay | 1696 | 118.9 | 5944 | 61.9 | 0 | — |
| 19 | −51.25 | −54.25 | −52.75 | L5 | clay | 1707 | 122.2 | 6108 | 63.2 | 0 | — |
| 20 | −54.25 | −57.25 | −55.75 | L5 | clay | 1718 | 125.4 | 6271 | 64.6 | 0 | — |
| 21 | −57.25 | −60.25 | −58.75 | L5 | clay | 1729 | 128.7 | 6435 | 66.0 | 0 | — |
| 22 | −60.25 | −63.25 | −61.75 | L5 | clay | 1740 | 132.0 | 6599 | 67.4 | 0 | — |
| 23 | −63.25 | −66.25 | −64.75 | L5 | clay | 1751 | 135.3 | 6763 | 68.8 | 0 | — |
| 24 | −66.25 | −69.25 | −67.75 | L5 | clay | 1762 | 138.5 | 6927 | 70.2 | 0 | — |
| 25 | −69.25 | −72.25 | −70.75 | L5 | clay | 1773 | 141.8 | 7090 | 71.5 | 0 | — |
| 26 | −72.25 | −75.25 | −73.75 | L5 | clay | 1784 | 145.1 | 7254 | 72.9 | 0 | — |
| 27 | −75.25 | −78.25 | −76.75 | L5 | clay | 1795 | 148.4 | 7418 | 74.3 | 0 | — |

Pile springs would follow the row `soilC` / `soilRho` / `soilPhi` once per-row mats exist.

---

## Constitutive parameters (N, m, s)

OpenSees wiki lists use **kN, m, s** (`ρ` in t/m³, stress in kPa). Convert: `ρ_SI = 1000 ρ`, `G_SI = 1000 G_kPa`, same for B and c.

Shared: `nd = 2`, `γ_max = 0.1`, `noYieldSurf = 20` (default). Clay `φ = 0`, `d = 0`, `p'_r = 100 kPa`. Sand `d = 0.5`, `p'_r = 101 kPa`.

### PIMY (clays)

`nDMaterial PressureIndependMultiYield tag 2 rho Gr Br c gam_max phi Pr d`

Knots for the ramps. Interpolation in \(z\) as above; \(B_r = 50 G_r\) on every clay row (Kramer Table 5.1 lists \(B/G=5\); we do not use that).

| Knot | Source | ρ (kg/m³) | G_r (Pa) | B_r (Pa) | c (Pa) | γ_max | φ | p'_r (Pa) | d |
|---|---|---|---|---|---|---|---|---|---|
| crust low | Kramer T5.1 | 1488 | 5.728e7 | 2.864e9 | 3.59e4 | 0.1 | 0 | 1.00e5 | 0 |
| medium low | Kramer T5.1 | 1521 | 6.639e7 | 3.320e9 | 3.97e4 | 0.1 | 0 | 1.00e5 | 0 |
| medium high | Kramer T5.1 | 1669 | 1.107e8 | 5.534e9 | 5.84e4 | 0.1 | 0 | 1.00e5 | 0 |
| wiki soft | OpenSees PIMY | 1300 | 1.30e7 | 6.50e8 | 1.80e4 | 0.1 | 0 | 1.00e5 | 0 |
| wiki stiff | OpenSees PIMY | 1800 | 1.50e8 | 7.50e9 | 7.50e4 | 0.1 | 0 | 1.00e5 | 0 |

Clay **B_r = 50 G_r** (not 60), i.e. undrained Poisson ν ≈ 0.49:

\[
\frac{B}{G}=\frac{2(1+\nu)}{3(1-2\nu)}\Big|_{\nu=0.49}=50
\qquad\Rightarrow\qquad
B_\mathrm{r}=50\,G_\mathrm{r}
\]

Wiki soft/medium list smaller B_r (ν ≈ 0.41); we keep undrained B. Wiki stiff lists B_r = 7.5e8; we replace it with **50 G_r**.

**No `FluidSolidPorous` on any PIMY** (L2, soft, L3 clay, L5 stiff). Undrained response is in B_r, not an FSP wrap.

### PDMY02 (sands)

`nDMaterial PressureDependMultiYield02 tag 2 rho Gr Br phi gam_max Pr d PTA contr1 contr3 dilat1 dilat3`

Wiki knots. Option 1 L3 interpolates in \(D_r\) between these; L5 stays on the PEER dense row.

| Knot | Source | ρ | G_r | B_r | φ | γ_max | p'_r | d | PTAng | contr1 | contr3 | dilat1 | dilat3 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 50% | wiki | 1900 | 1.00e8 | 2.33e8 | 33.5 | 0.1 | 1.01e5 | 0.5 | 25.5 | 0.045 | 0.15 | 0.06 | 0.15 |
| 60% | wiki | 2000 | 1.10e8 | 2.40e8 | 35.0 | 0.1 | 1.01e5 | 0.5 | 26 | 0.028 | 0.05 | 0.10 | 0.05 |
| 75% | wiki | 2100 | 1.30e8 | 2.60e8 | 36.5 | 0.1 | 1.01e5 | 0.5 | 26 | 0.013 | 0.0 | 0.30 | 0.0 |
| L5 | PEER dense sand | 2260 | 1.548e8 | 2.792e8 | 39.3 | 0.1 | 1.01e5 | 0.5 | 26 | 0.0042 | 0.0 | 0.564 | 0.0 |

ρ in kg/m³; G_r, B_r, p'_r in Pa. Leave `contrac2`, `dilat2`, `liquefac*` at OpenSees defaults.
Void ratio **`e` vs Dr** (wiki): 50% → 0.70, 60% → 0.65, 75% → 0.55, L5 0.43 — linear in \(D_r\) on L3, passed explicitly in `BuildSoilMaterials.tcl` (not the OpenSees default `e=0.6`).

### FluidSolidPorous (sands only)

`nDMaterial FluidSolidPorousMaterial tagFSP 2 soilMatTag combinedBulkModul`

| | This model | Wiki note |
|---|---|---|
| `combinedBulkModul` | **2.2e9 Pa** (2.2e6 kPa) | wiki: B_c ≈ B_f / n with B_f = 2.2e6 kPa; we use B_f (common PEER OpenSees sand practice) |
| `pa` | 1.01e5 Pa (default) | |
| Wrap | L3 and L5 when those layers are sand (opt. 1; opt. 2 L5 only) | **never** on PIMY clays |
| Element body force | **buoyant** γ' = (ρ − 1000) g on FSP sands | wiki FSP note; sand does ρ−1 t/m³ in kN·m·s |
| Clay body force | total γ = ρ g (no FSP) | undrained clay |
| Ponding `h_water` | Consistent \(F_y\) on soil **y=0** top edges: \(p=\rho_w g h_w\), half to each node × **`t_soil`** (never Shin \(t_{FF}\)) | `WaterSurfaceLoad.tcl`; no structure hydro; body forces unchanged |

Quads reference the **FSP tag**, not the bare PDMY tag. `updateMaterialStage` on both the solid and the FSP wrapper (stage 0 gravity → stage 1 EQ).

---

## Switch

One knob, e.g. `soilProfile 1` | `2` | `3` | `4`:

- **1** → L3 sand (PDMY + FSP, \(D_r\) ramp), L5 dense sand (PDMY + FSP, constant)
- **2** → L3 medium clay (PIMY ramp), L5 dense sand (PDMY + FSP)
- **3** → L3 medium clay (PIMY ramp), L5 stiff clay (PIMY ramp, no FSP)
- **4** → thicker soft L2 (PIMY ramp to −27.25 ft), medium L3 (PIMY ramp), stiff L5 (PIMY ramp)

L2 is crust clay (opts 1–3) or a deeper soft crust (opt. 4). Same nodes and quad connectivity.

---

## Lateral domain and free-field columns (Shin et al. 2007)

Shin et al. (2007): embankment-side soil extended **73.2 m (240 ft)** outward from the slope crest; outermost soil columns thickened out-of-plane and constrained so nodes at the same elevation share the same horizontal motion. Here the reference is the **pier center** (no abutment slope).

| Quantity | Meaning | Value |
|---|---|---|
| Near-field OOP thickness `t` | Plane-strain thickness of near-field quads | **45.72 m (150 ft)** = pier spacing |
| Free-field OOP thickness `t_FF` | Outermost columns (Shin) | **`10000 t` = 4.572×10⁵ m** |
| Near-field half-extent `L_half` | Pier center → NF outer face | **200 ft (61.0 m)** |
| FF column width `w_FF` | Shin BC strip beyond NF | **40 ft (12.19 m)** |
| Shin outer face | `L_half + w_FF` | **240 ft (73.2 m)** |

So each side: near field out to `L_half`, then (Shin only) a `w_FF` column. Full Shin soil width = 2 × (`L_half+w_FF`) = **480 ft**. ASDEA continuum stops at `L_half`; ASD sits outside.

```
  x = −(L_half+w_FF)              x = 0              x = +(L_half+w_FF)   [Shin]
  x = −L_half                     x = 0                     x = +L_half   [ASDEA NF]
  |← w_FF →|← near field L_half →|pier|← near field L_half →|← w_FF →|
  t_FF = 10000 t       t = 45.72 m            t = 45.72 m            t_FF = 10000 t
  equalDOF on column faces (retained = outer): UX+UY above base; base UX only
```

**Why not Shin’s 100 m for `t`?** Shin’s cut is **longitudinal**: out-of-plane is the **transverse embankment** (`embank_width × factor` → 100 m). We are **transverse** 2D: out-of-plane is **along the bridge**. Then `t` is the soil/deck tributary of **one pier**, so it should not exceed the longitudinal pier spacing.

Piers are at **150 ft (45.72 m)** centers (Mackie et al. 2008 / Shin geometry). One pier owns the strip from midspan to midspan → thickness **`t = S = 150 ft`**, not `S/2`. Half-spacing (75 ft) is the distance to the next midspan; the full strip is both sides.

Upper bound: `t ≤ 150 ft`. Cap length `L_cap = 10 ft` is only the foundation footprint — too small for soil mass. Use **`t = 150 ft`**.

End treatment: thick outer column + `equalDOF` on its two vertical faces → 1D free-field shear column. Retained = outermost (Lysmer). Shin’s wording is horizontal only; we tie UX and UY above the base. At the base, UX only so UY `fix` is not doubled with an MP.

Base: `fix` UX/UY on bottom nodes (near field and FF).

Horizontal mesh: finer near piles, coarsen toward the FF. Vertical mesh on the 3 ft pile stations (profiles above).

### Knobs (`Parameters.tcl`)

| Name | Meaning | Value |
|---|---|---|
| `t_soil` | Near-field OOP thickness | **45.72 m (150 ft)** |
| `t_FF_factor` | `t_FF / t_soil` | 10000 |
| `L_half` | Pier center → near-field outer face | **200 ft (61.0 m)** |
| `w_FF` | Shin FF column width beyond NF | **40 ft** |

ASDEA continuum = \(\pm L_\mathrm{half}\). Shin continuum = \(\pm(L_\mathrm{half}+w_{FF})\).
