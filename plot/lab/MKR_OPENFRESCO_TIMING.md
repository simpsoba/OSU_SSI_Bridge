# MKR-α + OpenFresco: what the actuator does each step

Walkthrough of how `MKRAlphaExplicitMultiSOE` / `CudaMKRAlpha` couples to
`expElement generic` + SCRAMNetGT + Simulink predictor–corrector in this repo.
Companion to `STATEOS_SIGNALS.md` (how to read `stateOS` on the mats).

Wiring in `Run.tcl` / `RunParallel.tcl` (`realTimeON 1`):

- trial → lab: UX displacement only (`expControlPoint` disp)
- lab → host: force only (`daqForce`)
- no time / \(\alpha_f\) signal on SCRAMNet

---

## 1. What MKR asks for

MKR forms equilibrium at the weighted station \(t_{\mathrm{int},n} + \alpha_f\Delta t\), not at
\(t_{\mathrm{int},n+1}\). For `MKRAlphaExplicitMultiSOE` / `CudaMKRAlpha` with
\(\rho_\infty^{\mathrm{eq}}=0.5\),

\[
\alpha_f=\frac{1}{1+\sqrt{\rho_\infty^{\mathrm{eq}}}}=2-\sqrt{2}\approx 0.5858
\]

(not the KR-α map \(\alpha_f=1/(1+\rho)=2/3\)).

Inside `newStep` (ExplicitAlpha / MKR family):

1. **Compute** \(\mathbf{u}_{n+1}\) and \(\mathbf{v}_{n+1}\) from the previous acceleration
   (explicit — no lab call yet).
2. Form \(\mathbf{u}_{n+\alpha_f} = (1-\alpha_f)\mathbf{u}_n + \alpha_f \mathbf{u}_{n+1}\).
3. Put that on the nodes; advance domain time to \(t_{\mathrm{int},n} + \alpha_f\Delta t\).
4. `updateDomain` → `EEGeneric::update` → `setTrialResponse` (ctrlDisp =
   \(u_{n+\alpha_f}\)).
5. Residual / `getResistingForce` → `acquire()` waits for `atTarget`, then
   reads `daqForce`.
6. That force is used to **solve \(\mathbf{a}_{n+1}\)**.
7. `commit`: domain time \(+= (1-\alpha_f)\Delta t\) → \(t_{\mathrm{int},n+1}\). No second
   lab move (`ECSCRAMNetGT::commitState` is a no-op).

So the physical force updates **acceleration**. Displacement \(\mathbf{u}_{n+1}\) was
already predicted before the handshake.

Numerical elements evaluate \(P_r(\mathbf{u}_{n+\alpha_f})\) constitutively. The
generic element does not: it commands \(u_{n+\alpha_f}\) and returns the
measured force at `atTarget`.

---

## 2. What the lab knows

Simulink never sees \(n\), \(n+\alpha_f\), or \(n+1\). Each handshake is:

1. Host writes target disp, raises `newTarget` (`switchPC` ack).
2. Over the current \(\Delta t_\mathrm{sim}\) window: **extrapolate** until the
   target lands, then **interpolate** toward it (`typeConv3` 1 → 0).
3. End of window: `atTarget`; host reads force.

Arrival deadline = end of that \(\Delta t_\mathrm{sim}\), not “OpenSees
\(t_{\mathrm{int},n+1}\)”. The value tracked is whatever OpenSees sent — for MKR,
\(u_{n+\alpha_f}\).

Default `-updateElemDisp` off: after the force read, commit does not send
\(u_{n+1}\). The actuator walks successive **\(\alpha\)-stations**, not the
pure nodal \(\mathbf{u}_n, \mathbf{u}_{n+1}\) (unless \(\alpha_f = 1\)).

“Behind” relative to host end-of-step: when OpenSees has committed \(\mathbf{u}_{n+1}\),
the last imposed target is still \(u_{n+\alpha_f}\). The DOF does not sit
frozen — extrapolation / interpolation at \(\Delta t_\mathrm{con}\) keeps
`comSig` moving (unless slowdown).

---

## 3. Mockup (no similitude, with solve time)

Assumptions (illustrative sub-ms stamps):

| Quantity | Value |
|----------|--------|
| \(\Delta t = \Delta t_\mathrm{int} = \Delta t_\mathrm{sim}\) | \(0.010\,\mathrm{s}\) |
| \(\Delta t_\mathrm{con}\) | \(0.001\,\mathrm{s}\) → \(N = 10\) counts / window |
| \(\alpha_f\) | \(2-\sqrt{2}\approx 0.5858\) → \(\alpha_f\Delta t \approx 0.0059\,\mathrm{s}\) (mock \(\Delta t=0.010\)) |
| Host work after each force | ~\(4\,\mathrm{ms}\) / count ~4 (read / solve / commit / predict / form / send) — ~40% of window extrapolate, ~60% interpolate (healthy-run order; cf. `STATEOS_SIGNALS.md` ~37% extrap) |
| First begin-step (before any force) | ~\(2\,\mathrm{ms}\) (compute \(\mathbf{u}_1, \mathbf{v}_1\) / form \(\mathbf{u}_{0+\alpha_f}\) / send) — not instantaneous; no prior solve |
| Healthy run | later targets land ~count 4; no slowdown |

Clocks: \(t\) = lab / wall time; \(t_\mathrm{int}\) = OpenSees domain (integration) time.
“Actuator” below means the command path (`comSig`); measured disp lags a bit.

### Start

| \(t\) | Integrator \(t_\mathrm{int}\) | Lab | Host |
|---:|---:|---|---|
| 0.0000 | 0.0000 | at \(u_0\) | committed \(\mathbf{u}_0,\mathbf{v}_0,\mathbf{a}_0\) |

### Enter step \(0\to1\) (first begin-step — a bit expensive)

No prior force/solve; still more host work than the short predict–form–send
tail after later commits (cold start / first `newStep`).

| \(t\) | Integrator \(t_\mathrm{int}\) | Lab | Host |
|---:|---:|---|---|
| 0.0000 | 0.0000 | window 0, count 1 | **begin step:** compute \(\mathbf{u}_1, \mathbf{v}_1\) from \(\mathbf{a}_0\) |
| 0.0008 | 0.0000 | count 1–2 | form \(\mathbf{u}_{0+\alpha_f}\) (in progress) |
| 0.0015 | **0.0059** | count 2 | set domain time to \(0+\alpha_f\) |
| **0.0020** | 0.0059 | receives target | **send** \(u_{0+\alpha_f}\) |
| 0.0021 | 0.0059 | extrap → interp | **block** on `atTarget` |

### Window 0 (host waiting)

| \(t\) | count | Lab | Integrator \(t_\mathrm{int}\) | Host |
|---:|---:|---|---:|---|
| 0.0021–0.0060 | 3–7 | interp toward \(u_{0+\alpha_f}\) | 0.0059 | blocked |
| 0.0070–0.0090 | 8–10 | interp → \(u_{0+\alpha_f}\) | 0.0059 | blocked |
| **0.0100** | end | **`atTarget`** | 0.0059 | unblock |

Force for the residual is **not** mid-flight: `acquire()` waits for
`atTarget`, then reads daq.

### After force 0 (~4 ms into window 1 — ~40% extrap)

| \(t\) | count | Lab | Integrator \(t_\mathrm{int}\) | Host |
|---:|---:|---|---:|---|
| 0.0100 | 1 | window 1 starts; extrapolate | 0.0059 | **read** \(F\) at \(u_{0+\alpha_f}\) |
| 0.0108 | 1–2 | extrapolate | 0.0059 | **solve** \(\mathbf{a}_1\) (in progress) |
| 0.0114 | 2 | extrapolate | 0.0059 | **solve** \(\mathbf{a}_1\) (done) |
| 0.0117 | 2 | extrapolate | **0.0100** | **commit** \(n=1\) |
| 0.0123 | 3 | extrapolate | 0.0100 | **compute** \(\mathbf{u}_2, \mathbf{v}_2\) from \(\mathbf{a}_1\) |
| 0.0129 | 3 | extrapolate | 0.0100 | **form** \(\mathbf{u}_{1+\alpha_f}\) |
| 0.0135 | 4 | extrapolate | **0.0159** | set domain time to \(1+\alpha_f\) |
| **0.0140** | 4 | flag → **interpolate** | 0.0159 | **send** \(u_{1+\alpha_f}\) |
| 0.0141 | 4 | interp toward \(u_{1+\alpha_f}\) | 0.0159 | **block** on `atTarget` |

### Rest of window 1

| \(t\) | count | Lab | Integrator \(t_\mathrm{int}\) | Host |
|---:|---:|---|---:|---|
| 0.0141–0.0190 | 4–10 | interp → \(u_{1+\alpha_f}\) | 0.0159 | blocked |
| **0.0200** | end | **`atTarget`** | 0.0159 | unblock |

### After force 1 (same ~4 ms burst into window 2)

| \(t\) | count | Lab | Integrator \(t_\mathrm{int}\) | Host |
|---:|---:|---|---:|---|
| 0.0200 | 1 | window 2; extrapolate | 0.0159 | **read** \(F\) at \(u_{1+\alpha_f}\) |
| 0.0208 | 1–2 | extrapolate | 0.0159 | **solve** \(\mathbf{a}_2\) (in progress) |
| 0.0214 | 2 | extrapolate | 0.0159 | **solve** \(\mathbf{a}_2\) (done) |
| 0.0217 | 2 | extrapolate | **0.0200** | **commit** \(n=2\) |
| 0.0223 | 3 | extrapolate | 0.0200 | **compute** \(\mathbf{u}_3, \mathbf{v}_3\) |
| 0.0229 | 3 | extrapolate | 0.0200 | **form** \(\mathbf{u}_{2+\alpha_f}\) |
| 0.0235 | 4 | extrapolate | **0.0259** | set domain time to \(2+\alpha_f\) |
| **0.0240** | 4 | → **interpolate** | 0.0259 | **send** \(u_{2+\alpha_f}\) |
| 0.0241 | 4 | interp toward \(u_{2+\alpha_f}\) | 0.0259 | **block** on `atTarget` |

### Strip

![MKR-α + OpenFresco timing strip](mkr_openfresco_timing_strip.svg)

*MKR-α + OpenFresco handshake timing (6 in, 9 pt). Wall time \(t\) in units of \(\Delta t_\mathrm{con}\): 0 at the start of init, \(k\) at the end of the cold start (factorize, not to scale), then (k+N), (k+2N), (k+3N), with (Delta t_mathrm{sim}=N,Delta t_mathrm{con}) (drawn with (N=10)). Integrator time \(t_\mathrm{int}\) (OpenSees domain time) is in units of \(\Delta t_\mathrm{int}=\sqrt{\lambda}\,\Delta t_\mathrm{sim}\); circles mark \(t_\mathrm{int}=n+\alpha_f\) at each at-target instant (MKR \(\alpha_f=1/(1+\sqrt{\rho_\infty})\approx 0.586\) for \(\rho_\infty=0.5\)). The target machine sends the reference to the controller every \(\Delta t_\mathrm{con}\), interpolating toward the latest target and extrapolating while the host works. The zoom shows the burst after the first at target: finish step 0 (\(\approx 3\): read force, solve \(\mathbf{a}_1\), short commit of \(\mathbf{u}_1, \mathbf{v}_1, \mathbf{a}_1\)), then begin step 1 (\(\approx 1\): compute \(\mathbf{u}_2\) and \(\mathbf{v}_2\) from \(\mathbf{a}_1\), form \(\mathbf{u}_{1+\alpha_f}\) and \(\mathbf{v}_{1+\alpha_f}\), advance \(t_\mathrm{int}\), send). Interface \(u\) (scalar, the actuated DOF): the target \(u_{n+\alpha_f}\) steps at each send; reference samples every \(\Delta t_\mathrm{con}\) (circles = interpolate, squares aligned with the curve = extrapolate; open = on target, solid = off target) reach the target only at the at-target instants (open circles). Dashed blue: where extrapolation would have continued had no new target arrived; dashed green: the unused start of each interpolating 2nd-order Lagrange fit, from the previous target to the send. Targets are illustrative (\(u_{0+\alpha_f}=1.02\), \(u_{1+\alpha_f}=1.29\), \(u_{2+\alpha_f}=0.27\)), not the slow sine of `plot_mkr_disp_vs_tw.py`.*

Labels ↔ code (three-loop terms): target = `tar` (target signal); reference = `comSig` (rate-transition output to the controller every \(\Delta t_\mathrm{con}\)); measured force = `daqForce` (measured signal); at target = SCRAMNet `atTarget`; read force = `getForce`; begin step = `newStep`; advance \(t_\mathrm{int}\) = `updateDomain`; factorize = `formOperators` (assemble and factorize M, the α-operator, and A).

### Integrator \(t_\mathrm{int}\) and displacements vs \(t\)

Three-window **mockup** with the same clocks as the strip
(\(\Delta t_\mathrm{sim}=10\,\Delta t_\mathrm{con}\), first begin \(\approx 2\),
later burst = finish \(\approx 3\) + begin \(\approx 1\) / ~40% extrap,
MKR \(\alpha_f=1/(1+\sqrt{0.5})\approx 0.5858\), no slowdown). \(u_{n+\alpha_f}\)
from a slow sine (peak at \(t=60\)). Script: `plot_mkr_disp_vs_tw.py`.

![Domain time and tar / comSig vs wall time](mkr_disp_vs_tw.png)

Top: integrator time \(t_\mathrm{int}\) holds at \(n+\alpha_f\) through at target (force is read at
that station), then jumps at commit and at the next send. Middle / bottom:
**tar** \(u_{n+\alpha_f}\) (stepwise at send — lab instruction) and **comSig**
(**2nd-order Lagrange**, Seki et al. 2026 §2.1: one fit through interpolate,
continues as extrapolate past at target; refit when the next target lands).
Axis \(t/\Delta t_\mathrm{con}\in[0,30]\); bottom scaled to prototype
(\(\times\lambda\), \(\lambda=2.4\)).

Glossary: begin step = OpenSees `newStep`; read force = `getForce`; advance
\(t_\mathrm{int}\) = `updateDomain`; at target = SCRAMNet `atTarget`; send = SCRAMNet
target write.

Pattern after the first send: **wait → read \(F\) → solve \(\mathbf{a}\) →
commit → compute → form \(\mathbf{u}_{\cdot+\alpha_f}\) → send → wait**. Lab
**extrapolates** during the host burst, then **interpolates** once the new
target lands.

### Actuator vs host stations (ideal tracking)

| Host station | Integrator \(t_\mathrm{int}\) | Wall when committed / forced | Actuator ≈ |
|---|---:|---|---|
| \(n=0\) | 0.0000 | \(t=0\) | \(u_0\) |
| send \(u_{0+\alpha_f}\) | 0.0059 | \(t\approx 0.002\) (first begin) | leaving \(u_0\) toward \(u_{0+\alpha_f}\) |
| \(0+\alpha_f\) (force) | 0.0059 | `atTarget` at \(t=0.010\) | \(u_{0+\alpha_f}\) |
| \(n=1\) (commit) | 0.0100 | \(t\approx 0.011\) | still \(u_{0+\alpha_f}\) (extrap starting) |
| send \(u_{1+\alpha_f}\) | 0.0159 | \(t=0.014\) | leaving \(u_{0+\alpha_f}\) toward \(u_{1+\alpha_f}\) |
| \(1+\alpha_f\) (force) | 0.0159 | \(t=0.020\) | \(u_{1+\alpha_f}\) |
| \(n=2\) (commit) | 0.0200 | \(t\approx 0.021\) | still \(u_{1+\alpha_f}\) |

### Real-time budget (this mockup)

- Window = 10 counts; host burst ≈ 4 → ~40% extrapolate, ~60% interpolate.
- Slowdown threshold ~0.8 \(\Delta t_\mathrm{sim}\) (`typeConv3 = 2`) — count 4
  is comfortable; ~count 9 would land near the deadline.
- During the wait, OpenSees is mostly blocked on `atTarget`. Solve overlaps
  the **start** of the next window (lab already extrapolating), not a serial
  “move, then solve, then move” with a frozen actuator.

With Froude similitude, \(\Delta t_\mathrm{int}\) (prototype) and
\(\Delta t_\mathrm{sim}\) (model) differ by \(\sqrt{\lambda}\); the handshake
logic is the same, only the rate mapping changes (`HybridCtrlParameters.rate`
in the Speedgoat `initializeSimulation.m`).

---

## 4. Comparing `tarSig` to `comSig` / `meaSig` (lab clock offset)

`tarSig`, `comSig`, and `meaSig` share **lab Time**. OpenSees domain time does
not. At session start the lab clock already ticks while the host finishes
startup (first factorization / first send); that pad is a one-shot sync debt,
not mid-run slowdown. Leave lab Time alone; shift the host target stations
when overlaying.

### What `tarSig` is

For MKR, each host target is \(u_{n+\alpha_f}\) (not committed \(\mathbf{u}_{n+1}\)).
On the mats, new values appear with `typeConv3` **1→0** (and the `typeConv2/s1`
pulse). In OpenSees time, sends sit at

\[
t_{\mathrm{int},k}=(k+\alpha_f)\,\Delta t_{\mathrm{sim}},\qquad k=0,1,2,\ldots
\]

(\(\Delta t_{\mathrm{sim}}=N\,\Delta t_{\mathrm{con}}\); on F05, \(N=10\),
\(\Delta t_{\mathrm{con}}\approx 0.488\,\mathrm{ms}\)).

### Offset (no \(\alpha_f\) in the definition)

Use the first **completed lab window** after the first real target:

1. \(t_{\mathrm{land}}\) — first `typeConv3` 1→0 with a non-tiny `tar`.
2. \(t_{\mathrm{ex}}\) — first 0→1 (into extrapolate) after that land.
3. Window end \(t_{\mathrm{goal}}=t_{\mathrm{ex}}-\Delta t_{\mathrm{con}}\)
   (the 0→1 sample is already one \(\Delta t_{\mathrm{con}}\) into the next
   window; at target (`atTarget`) is the count-10 end of the prior window).
4. Pin OpenSees \(t_{\mathrm{int}}=\Delta t_{\mathrm{sim}}\) to that window end:

\[
\mathrm{offset}
= t_{\mathrm{ex}}-\Delta t_{\mathrm{con}}-\Delta t_{\mathrm{sim}}
= t_{\mathrm{goal}}-\Delta t_{\mathrm{sim}}.
\]

\(\alpha_f\) does **not** enter this formula. On F05 (rowNeg4 / `r-04_…1010`):
\(\mathrm{offset}\approx 187.5\,\mathrm{ms}\).

### Where to put each \(u_{n+\alpha_f}\) on the lab axis

\[
t_k=\mathrm{offset}+(k+\alpha_f)\,\Delta t_{\mathrm{sim}}.
\]

Plot those stations against `comSig` / `meaSig` on raw lab Time. Scripts:
`plot/lab/plot_clock_offset_start.py`, `plot/lab/plot_mkr_disp_vs_tw_real.py`
(figures `clock_offset_start.*`, `mkr_disp_vs_tw_real.*`).

### Expected gap when there is no slowdown

Lab finishes the trial at window end (\(t_{\mathrm{int},n}+\Delta t_{\mathrm{sim}}\) in the
synced OS clock). The send station is at \(t_{\mathrm{int},n}+\alpha_f\Delta t_{\mathrm{sim}}\).
So `comSig` should reach that tar level about

\[
(1-\alpha_f)\,\Delta t_{\mathrm{sim}}
\]

later (\(\approx 2.02\,\mathrm{ms}\) on F05 with MKR \(\alpha_f\approx 0.5858\)).
Checked on F05 no-slowdown windows with that \(\alpha_f\) on the OS send grid.
First-1 mm com−tar is the same order.

Do **not** use pier-top UX for this check: the pier recorder is committed
\(u_{n+1}\) samples, not the \(\alpha\)-stations on SCRAMNet. Use `tarSig`.

### Related (not the same)

| Quantity | Role |
|----------|------|
| Mat − pier record-length gap | ≈ startup pad; rough peer of `offset`, not the definition above |
| `t_land - α_f Δt_sim` | Alternate offset that *builds in* \(\alpha_f\); prefer window-end form so \(\alpha_f\) is tested, not assumed |
| Mid-run `typeConv3→2` | Real-time slowdowns; keep out of the startup offset |

### Post-process: map OpenSees recorders onto lab Time

For RTHS history / PSD / hyst figures under `plots/runs/<Test>/os/`, use

\[
t_{\mathrm{lab}}=k+\frac{t_{\mathrm{OS}}}{\sqrt{\lambda}}+f
\]

with \(k\) from the window-end offset above (no `|tar|` cut) and \(f\) = cumulative
**mid-run** `typeConv3==2` only. Helper: `plot/lab_time_map.py`. Used by
`PlotActuatorForce`, `PlotPierBaseForce`, `PlotActuatorVsPier`, `PlotHydroSpectra`,
`PlotActuatorHyst`, `PlotLabTimeMapDiag`, and compare `*_opensees.png` companions.
Leave `eq/` on domain time; PlotMatOS / Simulink pairs stay on native lab Time.

---

## 5. Code anchors

| Piece | Where |
|-------|--------|
| OpenFresco UX disp / force, `-checkTime` | `Run.tcl`, `RunParallel.tcl` (`realTimeON`) |
| MKR `newStep` / \(\mathbf{u}_{n+\alpha_f}\) / commit | simpsoba `ExplicitAlphaMultiSOE.cpp` (MKR family) |
| EEGeneric send trial / get force | OpenFresco `EEGeneric.cpp` |
| SCRAMNet wait on `atTarget` | `ECSCRAMNetGT::acquire()` |
| `stateOS` extrap / interp / slowdown | `plot/lab/STATEOS_SIGNALS.md` |
| Lab↔OS offset; tar vs com/mea | §4; `plot_clock_offset_start.py`, `plot_mkr_disp_vs_tw_real.py` |
| OS→lab map for `os/` plotters | §4 post-process; `plot/lab_time_map.py` |

Stock OpenFresco HybridController / PredictorCorrector (no \(\alpha_f\) input):
`simpsoba/RTHS-CUDA/OpenFresco/SRC/experimentalControl/Simulink/`.
Seki campaign init scripts (same handshake map): Shared Drive
`2024 - FakeWind - rths - seki/model/speedgoat/`,
`2022 - Monopile OWT  - rths - Seki/RTActualTestModels/`.
