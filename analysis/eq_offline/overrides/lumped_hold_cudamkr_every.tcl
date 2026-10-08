# Serial lumped-hold, DT_FACTOR 10, first 180 s, every analysis step.
# Dump: plot/out/eq_offline/lumped_hold_cudamkr_every/  Does not set dTRecorder.
set runEQ 1
set plotFigures 0
set realTimeON 0
set recordersON 1
set eqPrintON 1
set eqPrintDt 10.0
set holdPierON 1
set holdPierRZON 0
set pierEleType "lumpedPlasticity"
set soilMesh 0
set soilProfile 4
set soilBoundary "Shin"
set soilEleType "SSPquad"
set soilConstitutive "inelastic"
set eqIntegrator "CudaMKRAlpha 0.5"
set prePartitionSystem "UmfPack"
set postPartitionSystem "CuDSS"
set constraintsHandler "Transformation"
set DT_FACTOR 10
set eqTmax 180.0
set eqFreeVibT 0
set gmDir [file join $gmRoot Tohoku2011-FKSH]
set gmVelFile [file join $gmDir FKSH19.NS1.VT2]
set gmStartTime 0.0
set outDIR "plot/out/eq_offline/lumped_hold_cudamkr_every"
