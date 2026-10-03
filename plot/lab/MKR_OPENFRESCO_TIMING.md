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

MKR forms equilibrium at the weighted station \(t_n + \alpha_f\Delta t\), not at
\(t_{n+1}\). With \(\rho_\infty = 0.5\), \(\alpha_f = 2/3\).

Inside `newStep` (ExplicitAlpha / MKR family):

1. **Predict** \(U_{n+1}\) and \(\dot U_{n+1}\) from the previous acceleration
   (explicit — no lab call yet).
2. Form \(U_{n+\alpha_f} = (1-\alpha_f)U_n + \alpha_f U_{n+1}\).
3. Put that on the nodes; advance domain time to \(t_n + \alpha_f\Delta t\).
4. `updateDomain` → `EEGeneric::update` → `setTrialResponse` (ctrlDisp =
   \(U_{n+\alpha_f}\)).
5. Residual / `getResistingForce` → `acquire()` waits for `atTarget`, then
   reads `daqForce`.
6. That force is used to **solve \(\ddot U_{n+1}\)**.
7. `commit`: domain time \(+= (1-\alpha_f)\Delta t\) → \(t_{n+1}\). No second
   lab move (`ECSCRAMNetGT::commitState` is a no-op).

So the physical force updates **acceleration**. Displacement \(U_{n+1}\) was
already predicted before the handshake.

Numerical elements evaluate \(P_r(U_{n+\alpha_f})\) constitutively. The
generic element does not: it commands \(U_{n+\alpha_f}\) and returns the
measured force at `atTarget`.

---

## 2. What the lab knows

Simulink never sees \(n\), \(n+\alpha_f\), or \(n+1\). Each handshake is:

1. Host writes target disp, raises `newTarget` (`switchPC` ack).
2. Over the current \(\Delta t_\mathrm{sim}\) window: **extrapolate** until the
   target lands, then **interpolate** toward it (`typeConv3` 1 → 0).
3. End of window: `atTarget`; host reads force.

Arrival deadline = end of that \(\Delta t_\mathrm{sim}\), not “OpenSees
\(t_{n+1}\)”. The value tracked is whatever OpenSees sent — for MKR,
\(U_{n+\alpha_f}\).

Default `-updateElemDisp` off: after the force read, commit does not send
\(U_{n+1}\). The actuator walks successive **\(\alpha\)-stations**, not the
pure nodal \(U_n, U_{n+1}\) (unless \(\alpha_f = 1\)).

“Behind” relative to host end-of-step: when OpenSees has committed \(U_{n+1}\),
the last imposed target is still \(U_{n+\alpha_f}\). The DOF does not sit
frozen — extrapolation / interpolation at \(\Delta t_\mathrm{con}\) keeps
`comSig` moving (unless slowdown).

---

## 3. Mockup (no similitude, with solve time)

Assumptions (illustrative sub-ms stamps):

| Quantity | Value |
|----------|--------|
| \(\Delta t = \Delta t_\mathrm{int} = \Delta t_\mathrm{sim}\) | \(0.010\,\mathrm{s}\) |
| \(\Delta t_\mathrm{con}\) | \(0.001\,\mathrm{s}\) → \(N = 10\) counts / window |
| \(\alpha_f\) | \(2/3\) → \(\alpha_f\Delta t \approx 0.0067\,\mathrm{s}\) |
| Host work after each force | ~\(3\,\mathrm{ms}\) (read / solve / commit / predict / form / send) |
| Healthy run | target lands ~count 4; no slowdown |

Clocks: \(T_w\) = lab / wall time; \(t\) = OpenSees domain time.
“Actuator” below means the command path (`comSig`); measured disp lags a bit.

### Start

| \(T_w\) | Domain \(t\) | Lab | Host |
|---:|---:|---|---|
| 0.0000 | 0.0000 | at \(U_0\) | committed \(U_0,\dot U_0,\ddot U_0\) |

### Enter step \(0\to1\)

| \(T_w\) | Domain \(t\) | Lab | Host |
|---:|---:|---|---|
| 0.0000 | 0.0000 | window 0, count 1 | predict \(U_1\) from \(\ddot U_0\) |
| 0.0001 | 0.0000 | count 1 | form \(U_{0+\alpha}\) |
| 0.0002 | **0.0067** | count 1 | set domain time to \(0+\alpha\) |
| 0.0003 | 0.0067 | receives target | **send** \(U_{0+\alpha}\) |
| 0.0004 | 0.0067 | | **block** on `atTarget` |

### Window 0 (host waiting)

| \(T_w\) | count | Lab | Domain \(t\) | Host |
|---:|---:|---|---:|---|
| 0.0004–0.0040 | 1–5 | extrap → interp | 0.0067 | blocked |
| 0.0050–0.0090 | 6–10 | interp → \(U_{0+\alpha}\) | 0.0067 | blocked |
| **0.0100** | end | **`atTarget`** | 0.0067 | unblock |

Force for the residual is **not** mid-flight: `acquire()` waits for
`atTarget`, then reads daq.

### After force 0 (~3 ms into window 1)

| \(T_w\) | count | Lab | Domain \(t\) | Host |
|---:|---:|---|---:|---|
| 0.0100 | 1 | window 1 starts; extrapolate | 0.0067 | **read** \(F\) at \(U_{0+\alpha}\) |
| 0.0104 | 1 | extrapolate | 0.0067 | **solve** \(\ddot U_1\) (in progress) |
| 0.0108 | 2 | extrapolate | 0.0067 | **solve** \(\ddot U_1\) (done) |
| 0.0111 | 2 | extrapolate | **0.0100** | **commit** \(n=1\) |
| 0.0114 | 2 | extrapolate | 0.0100 | **predict** \(U_2\) from \(\ddot U_1\) |
| 0.0117 | 2 | extrapolate | 0.0100 | **form** \(U_{1+\alpha}\) |
| 0.0120 | 3 | extrapolate | **0.0167** | set domain time to \(1+\alpha\) |
| 0.0125 | 3 | extrapolate | 0.0167 | prep SCRAMNet write |
| **0.0130** | 4 | flag → **interpolate** | 0.0167 | **send** \(U_{1+\alpha}\) |
| 0.0131 | 4 | interp toward \(U_{1+\alpha}\) | 0.0167 | **block** on `atTarget` |

### Rest of window 1

| \(T_w\) | count | Lab | Domain \(t\) | Host |
|---:|---:|---|---:|---|
| 0.0131–0.0190 | 4–10 | interp → \(U_{1+\alpha}\) | 0.0167 | blocked |
| **0.0200** | end | **`atTarget`** | 0.0167 | unblock |

### After force 1 (same burst into window 2)

| \(T_w\) | count | Lab | Domain \(t\) | Host |
|---:|---:|---|---:|---|
| 0.0200 | 1 | window 2; extrapolate | 0.0167 | **read** \(F\) at \(U_{1+\alpha}\) |
| 0.0204 | 1 | extrapolate | 0.0167 | **solve** \(\ddot U_2\) (in progress) |
| 0.0208 | 2 | extrapolate | 0.0167 | **solve** \(\ddot U_2\) (done) |
| 0.0211 | 2 | extrapolate | **0.0200** | **commit** \(n=2\) |
| 0.0214 | 2 | extrapolate | 0.0200 | **predict** \(U_3\) |
| 0.0217 | 2 | extrapolate | 0.0200 | **form** \(U_{2+\alpha}\) |
| 0.0220 | 3 | extrapolate | **0.0267** | set domain time to \(2+\alpha\) |
| 0.0225 | 3 | extrapolate | 0.0267 | prep SCRAMNet write |
| **0.0230** | 4 | → **interpolate** | 0.0267 | **send** \(U_{2+\alpha}\) |
| 0.0231 | 4 | interp toward \(U_{2+\alpha}\) | 0.0267 | **block** on `atTarget` |

### Strip

```text
T_w  0.000          0.010  0.011  0.013           0.020  0.021  0.023
     |-- window 0 --|---- window 1 ----|---- window 2 ----
lab  [  → U_{0+α}   ][extrap | → U_{1+α} ][extrap | → U_{2+α}
host [==== wait ====][rd sol cmt pred frm send][==== wait ====]...
                      10.0 10.8 11.1 11.4 11.7 13.0  ms
t    0 → 0.0067 ----→ 0.010 ------→ 0.0167 ----→ 0.020 → 0.0267
```

Pattern after the first send: **wait → read \(F\) → solve \(\ddot U\) →
commit → predict → form \(U_{\cdot+\alpha}\) → send → wait**. Lab
**extrapolates** during the host burst, then **interpolates** once the new
target lands.

### Actuator vs host stations (ideal tracking)

| Host station | Domain \(t\) | Wall when committed / forced | Actuator ≈ |
|---|---:|---|---|
| \(n=0\) | 0.0000 | \(T_w=0\) | \(U_0\) |
| \(0+\alpha\) (force) | 0.0067 | `atTarget` at \(T_w=0.010\) | \(U_{0+\alpha}\) |
| \(n=1\) (commit) | 0.0100 | \(T_w\approx 0.011\) | still \(U_{0+\alpha}\) (extrap starting) |
| send \(U_{1+\alpha}\) | 0.0167 | \(T_w=0.013\) | leaving \(U_{0+\alpha}\) toward \(U_{1+\alpha}\) |
| \(1+\alpha\) (force) | 0.0167 | \(T_w=0.020\) | \(U_{1+\alpha}\) |
| \(n=2\) (commit) | 0.0200 | \(T_w\approx 0.021\) | still \(U_{1+\alpha}\) |

### Real-time budget (this mockup)

- Window = 10 ms; host burst = 3 ms → target ~30% into the next window.
- Slowdown threshold ~0.8 \(\Delta t_\mathrm{sim}\) (`typeConv3 = 2`) — 3 ms is
  comfortable; ~9 ms would land near the deadline.
- During the wait, OpenSees is mostly blocked on `atTarget`. Solve overlaps
  the **start** of the next window (lab already extrapolating), not a serial
  “move, then solve, then move” with a frozen actuator.

With Froude similitude, \(\Delta t_\mathrm{int}\) (prototype) and
\(\Delta t_\mathrm{sim}\) (model) differ by \(\sqrt{\lambda}\); the handshake
logic is the same, only the rate mapping changes (`HybridCtrlParameters.rate`
in the Speedgoat `initializeSimulation.m`).

---

## 4. Code anchors

| Piece | Where |
|-------|--------|
| OpenFresco UX disp / force, `-checkTime` | `Run.tcl`, `RunParallel.tcl` (`realTimeON`) |
| MKR `newStep` / \(U_{n+\alpha_f}\) / commit | simpsoba `ExplicitAlphaMultiSOE.cpp` (MKR family) |
| EEGeneric send trial / get force | OpenFresco `EEGeneric.cpp` |
| SCRAMNet wait on `atTarget` | `ECSCRAMNetGT::acquire()` |
| `stateOS` extrap / interp / slowdown | `plot/lab/STATEOS_SIGNALS.md` |

Stock OpenFresco HybridController / PredictorCorrector (no \(\alpha_f\) input):
`simpsoba/RTHS-CUDA/OpenFresco/SRC/experimentalControl/Simulink/`.
Seki campaign init scripts (same handshake map): Shared Drive
`2024 - FakeWind - rths - seki/model/speedgoat/`,
`2022 - Monopile OWT  - rths - Seki/RTActualTestModels/`.
