# stateOS — OpenFresco signal generation (Seki §2.1)

How to read `hist_os_state.png` from `plot/PlotMatOS.py`. Companion to
Companion to `TestMatrix_lab_runs.csv` (as-run index) and `mat_run_map.json`
(orphan mats only).

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

| Column | Seki §2.1 role | What it does on these mats |
|--------|----------------|----------------------------|
| **`typeConv3`** | State enum | −1 / 0 / 1 / 2 (table below) |
| **`typeConv1`** | **count** | Substep index within one Δt_sim: cycles **1…10** (then resets). Slowdown sits at **9**. |
| **`typeConv2/s1`** | **flag** (arrival pulse) | Rises **1→0** with extrapolate→interpolate (host target landed). Brief pulse (~2 DAQ samples). |
| **`typeConv2/s2`** | Delayed flag copy | Same pulse as s1, delayed by ~1 count tick during interpolate entry. |
| **`typeConv2/s3`** | Wait / active-extrap marker | **1** for the whole extrapolate (and slowdown) stretch; drops after interpolate begins. Not the paper “flag” itself. |

One healthy Δt_sim cycle on r-17 (lab ~5.08 s; count 1…10):

```text
count  typeConv3   s1 s2 s3
  10   interpolate  0  0  0   ← end of previous window
   1   extrapolate  0  0  1   ← new window; waiting on host
   2   extrapolate  0  0  1
   3   extrapolate  0  0  1
   4   extrapolate  0  0  1
   5   interpolate  1  0  1   ← target in; s1 pulse starts
   6   interpolate  1  1  0
   7   interpolate  0  1  0
   8…10 interpolate 0  0  0   ← pure interpolate to end of window
```

Verified on `0821_GusBridge_rowNeg17_8Core_P1`: during extrapolate/slowdown,
s3≡1 and s1≡s2≡0; s1 rises only on **1→0** (into interpolate).

### `typeConv3` values

| Value | Seki state | Meaning |
|------:|------------|---------|
| **−1** | initialize | Before first target; zero references (startup only) |
| **0** | **interpolate** | New host target arrived; correct reference over remaining substeps |
| **1** | **extrapolate** | Host step in flight; Lagrange extrapolation of last targets |
| **2** | **slowdown** | Host late (~ past 0.8 Δt_sim); actuator slows / holds |

**Stateflow transitions (Seki §2.1):** extrapolate → interpolate *or* slowdown,
depending on whether the flag arrives before the count exceeds a fraction of
Δt_sim (~0.8). Slowdown recovers into interpolate when the flag eventually
turns 1. On the mats: every entry into **2** is from **1**; every exit from
**2** goes to **0** (checked on r-17, r-08, r+01, r-80, r-81, r-91, r-06).
The paper’s “flag” lines up with **s1** (and the delayed **s2**), not s3.

**Slowdown metric:** rising edge into `typeConv3 == 2` (see
`plot/PlotEQComparePairs.py`, `SLOWDOWN_STATE = 2`). Pairwise compare uses
the fewest such events inside GM **D5–95** as the interim reference.

Campaign stacked bars (`plot/PlotStateOSBars.py`): fixed **lab** window
**[35, 90] s** (semi–D5–95; not Froude-mapped Arias — slowdowns stretch wall
clock). Early stops show a hatched **unfinished** top. Full-record bars are
separate (`bar_typeconv3_full.png`). **Slowdown (2)** is split into **brief**
(1 DAQ sample per episode) and **sustained** (≥2 samples); labels on the
orange segments show combined slowdown % to three decimals.

During healthy runs the top panel flickers **1 ↔ 0** (extrapolate while the
host is integrating, then interpolate once the target lands). Sustained **2**
means the host could not keep real-time pace.

### Example: r-17 (Trial 30, np=8, no slowdown)

- ~63% interpolate (0), ~37% extrapolate (1), 0.06% slowdown (2) over the full
  record (sample fractions; interpolate often wins on DAQ samples when the
  target arrives early in Δt_sim).
- Two startup spikes to 2 at lab **t ≈ 0.01 s** and **0.19 s**; **zero**
  slowdown onsets inside D5–95 (lab 36.1–89.8 s).
- **`typeConv1`** cycles 1→10 each Δt_sim; mean count ≈ 2–4 during
  extrapolate (1), ≈ 7 during interpolate (0), = 9 during slowdown (2).

Heavy-slowdown mats (e.g. r-80, r-81) enter **2** often from **1**, each
spike brief on the sample clock.

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
