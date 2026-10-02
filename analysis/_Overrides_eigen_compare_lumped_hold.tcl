# Baseline mesh: lumpedPlasticity + holdPier (UX+UY). Temp compare only.
set runEQ 0
set plotFigures 1
set realTimeON 0
set recordersON 0
set holdPierON 1
set holdPierRZON 0
set pierEleType "lumpedPlasticity"
set nModesEigen 40
set soilMesh 0
set soilProfile 4
set soilEleType "SSPquad"
set soilConstitutive "inelastic"
set soilBoundary "Shin"
set outDIR "_eigen_compare_lumped_hold"
set eigenOutPath [file join $plotDir out eigen compare_lumped_hold eigen_modes.json]
file mkdir [file dirname $eigenOutPath]
