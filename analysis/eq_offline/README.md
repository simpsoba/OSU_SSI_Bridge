# Offline baseline EQ compare (temp)

Four structural cases on mesh 0 (profile 4, Shin, SSPquad), integrator
`MKRAlphaExplicitMultiSOE 0.5`, constraints `Transformation`:

| slug | holdPierON | pierEleType |
|---|---|---|
| lumped_hold_mkr | 1 | lumpedPlasticity |
| lumped_free_mkr | 0 | lumpedPlasticity |
| fbc_hold_mkr | 1 | forceBeamColumn |
| fbc_free_mkr | 0 | forceBeamColumn |

**Outputs only under** `plot/out/eq_offline/<slug>/`. Does not write to
`OSU_SSI_BRIDGE_DATA*`, Shared Drive, `plot/out/eigen/`, or `modal_analysis/`.

Needs an OpenSees build that ships MKR (stock PATH often does not). The runner
prefers `simpsoba/OpenSees/build-cuda/Release/OpenSees.exe`, or set
`OPENSEES_EXE`.

```bat
cd <repo>
python plot/RunEQOfflineBaseline.py --smoke
python plot/RunEQOfflineBaseline.py
python plot/PlotEQOfflineCompare.py --integrators mkr
```
