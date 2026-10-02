# Auto-written by plot/RunEQOfflineBaseline.py — offline baseline EQ only.
# Dump: plot/out/eq_offline/fbc_free_mkr/  (never lab data / eigen / modal_analysis)
set runEQ 1
set plotFigures 0
set realTimeON 0
set recordersON 1
set eqPrintON 1
set eqPrintDt 10.0
set holdPierON 0
set holdPierRZON 0
set pierEleType "forceBeamColumn"
set soilMesh 0
set soilProfile 4
set soilBoundary "Shin"
set soilEleType "SSPquad"
set soilConstitutive "inelastic"
set eqIntegrator "MKRAlphaExplicitMultiSOE 0.5"
set prePartitionSystem "UmfPack"
set postPartitionSystem "UmfPack"
set constraintsHandler "Transformation"
set gmDir [file join $gmRoot Tohoku2011-FKSH]
set gmVelFile [file join $gmDir FKSH19.NS1.VT2]
set gmStartTime 0.0
set outDIR "plot/out/eq_offline/fbc_free_mkr"
