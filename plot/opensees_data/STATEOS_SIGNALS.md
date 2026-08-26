# stateOS — OpenFresco signal generation (Seki §2.1)

How to read `hist_os_state.png` from `plot/PlotMatOS.py`. Companion to
`LAB_RUN_MAP.md` (slowdown counts) and `mat_run_map.json`.

**Paper:** Seki et al. (2026), *Computational capacity in hydrodynamic RTHS*,
§2.1 Signal generation task — Stateflow rate-transition on the **target**
(Simulink + OpenFresco), between the OpenSees **host** and the flume
**controller / DAQ**.

**Three-loop roles (Seki Fig. 2):**

| Loop | Machine | This campaign |
|------|---------|---------------|
| Integrator | Host | OpenSees MP (`np` cores); irregular **target** displacements |
| Rate transition | Target | Signal generation + sync; **reference** commands at Δt_con |
| Sensor–control | Controller / DAQ | Actuator + `*OS` recorders |

OpenSees does not implement extrapolation or slowdown — that logic runs on
the DSP so the actuator still gets deterministic references when the host
is late.

---

## `stateOS` columns (Simulink → mat extract)

Path prefix: `OpenSees//OpenFresco1/Subsystem4/…`

| Column | Campaign name | Seki §2.1 role |
|--------|---------------|----------------|
| **`typeConv3`** | state enum | Stateflow state (below) |
| **`typeConv1`** | substep / count | **count** — 0 … N_up×N each Δt_sim (here 0–10) |
| **`typeConv2/s1`, `/s2`** | flag s1, s2 | Duplicate **flag** outputs in our charts (identical on r-17) |
| **`typeConv2/s3`** | flag s3 | **flag** — 1 when a new host target has arrived |

Verified on `0821_GusBridge_rowNeg17_8Core_P1` (r-17): s3 = 1 iff
`typeConv3` ∈ {1, 2}; s1 and s2 track each other.

### `typeConv3` values

| Value | Seki state | Meaning |
|------:|------------|---------|
| **−1** | initialize | Before first target; zero references (startup only) |
| **0** | extrapolate | Host step in flight; Lagrange extrapolation of last targets |
| **1** | interpolate | New target arrived; correct reference over remaining substeps |
| **2** | **slowdown** | Host late (~ past 0.8 Δt_sim); actuator slows / holds |

**Slowdown metric:** rising edge into `typeConv3 == 2` (see
`plot/PlotEQComparePairs.py`, `SLOWDOWN_STATE = 2`). Pairwise compare uses
the fewest such events inside GM **D5–95** as the interim reference.

During healthy runs the top panel flickers **0 ↔ 1** (extrapolate most of
each OpenSees step, brief interpolate when the answer lands). Sustained **2**
means the host could not keep real-time pace.

### Example: r-17 (Trial 30, np=8, no slowdown)

- ~63% extrapolate (0), ~37% interpolate (1), 0.06% slowdown (2).
- Two startup spikes to 2 at lab **t ≈ 0.01 s** and **0.19 s**; **zero**
  slowdown onsets inside D5–95 (lab 36.1–89.8 s).
- **`typeConv1`** cycles 0→10 each Δt_sim window; count ≈ 7 during
  extrapolate, ≈ 2 during interpolate, ≈ 9 during slowdown.

Heavy-slowdown mats (e.g. r-80, r-81) spend most of the record in state 2.

---

## Related *OS signals (`hist_os_tar_com_mea.png`)

| Signal | Role |
|--------|------|
| **`tarSigOS`** | Irregular **targets** from OpenSees when they arrive |
| **`comSigOS`** | Smoothed **references** from signal generation (§2.1 output) |
| **`meaSigOS`** | Measured physical response |

§2.1 is the block that turns **tar → com**. Slowdown freezes **com** on lab
time while the host catches up.

---

## Clocks on `hist_os_state.png`

- **Bottom x:** lab **model** time (`Time` in the mat).
- **Top x:** prototype-scaled lab time, t_lab √λ (λ = 2.4, Froude).

OpenSees recorder time in `.out` files is **prototype** time. Do not overlay
OpenSees and Simulink traces without √λ (or use separate axes as in
`PlotEQComparePairs.py`).

---

## References

- Seki et al. (2026), WEER — `reference/references.bib` →
  `sekiComputationalCapacityHydrodynamic2026`
- Schellenberg et al. [47], Stojadinovic et al. [53] — Stateflow implementation
  cited in Seki §2.1
- Neumann et al. (2023) — same three-loop hydro-RTHS architecture for the
  OSU cascading EQ–tsunami bridge demo
