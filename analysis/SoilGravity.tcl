# analysis/SoilGravity.tcl
# Goals:
#   Apply soil self-weight (elastic, then plastic). No structure -mg yet.
#   Tie soil UX,UY to the coincident pile/cap so the continuum settles with
#   the shafts. Leave absorbing BCs in the gravity state.
# Call after WaterSurfaceLoad. Do not ActivateEQBoundary yet.

if {![info exists soilBoundary]} {
	error "SoilGravity.tcl: source Parameters.tcl first"
}
if {![info exists eleTag_soil_last]} {
	error "SoilGravity.tcl: BuildSoilMesh.tcl first"
}
if {![info exists soilDir]} {
	error "SoilGravity.tcl: set soilDir first"
}

# Pairs from BuildSoilSprings: {soilNd retainNd}. retainNd is the pile or
# cap at the same (x,y). equalDOF retain soil 1 2: soil UX,UY follow the
# pile/cap under soil self-weight (first node retained, second constrained).
# gravMpRemove / gravFoldNodes are used in FoldStructNodes.tcl after this.
if {![info exists gravSoilLockPairs]} {
	error "SoilGravity.tcl: gravSoilLockPairs missing (BuildSoilSprings.tcl first)"
}
set nGravLock 0
set gravMpRemove {}
set gravFoldNodes {}
array unset soilLocked
foreach pair $gravSoilLockPairs {
	lassign $pair soilNd retainNd
	if {[info exists soilLocked($soilNd)]} {
		error [format "SoilGravity.tcl: soil %d already constrained to %d (also %d)" \
			$soilNd $soilLocked($soilNd) $retainNd]
	}
	equalDOF $retainNd $soilNd 1 2
	set soilLocked($soilNd) $retainNd
	lappend gravMpRemove $soilNd
	lappend gravFoldNodes $retainNd
	incr nGravLock
}

# Pier/deck/cap/pile beam-columns are not in the model yet. Pin RZ so
# structure nodes cannot spin.
# Ride pier/deck UX,UY with a pile/cap already locked to soil (gravLockSeen
# from BuildSoilSprings; cap BC if it is in that set, else first retainNd).
# The retain node itself ties to nearest soil. Do not tie structure to soil
# nodes that already follow a pile: those continuum nodes are constrained.
if {![info exists structNodeTags] || [llength $structNodeTags] < 1} {
	error "SoilGravity.tcl: structNodeTags empty (BuildModel.tcl first)"
}
set gravStructFixRZ {}
if {![array exists gravLockSeen]} {
	array set gravLockSeen {}
}
set gravRetainNd ""
if {$pileSpring ne "none"} {
	if {[info exists nodeTag_cap_BC] && [info exists gravLockSeen($nodeTag_cap_BC)]} {
		set gravRetainNd $nodeTag_cap_BC
	} elseif {[llength $gravSoilLockPairs] > 0} {
		set gravRetainNd [lindex [lindex $gravSoilLockPairs 0] 1]
	}
}
foreach n $structNodeTags {
	fix $n 0 0 1
	lappend gravStructFixRZ $n
	if {[info exists gravLockSeen($n)]} { continue }
	if {$gravRetainNd ne "" && $n != $gravRetainNd} {
		equalDOF $gravRetainNd $n 1 2
	} else {
		set soilNd [soilNdNearestExisting 0.0 0.0]
		if {$soilNd < 0} {
			error "SoilGravity.tcl: no soil node to lock structure $n"
		}
		equalDOF $soilNd $n 1 2
	}
	lappend gravMpRemove $n
	lappend gravFoldNodes $n
	set gravLockSeen($n) 1
	incr nGravLock
}

# Continuum elastic (Liq stay stage 0)
set updateLiqSpringStage 0
set soilMatStageWanted 0
source [file join $soilDir UpdateSoilStage.tcl]

setTime 0.0
set dLambda 0.1
set nStep [expr {int(round(1.0/$dLambda))}]

wipeAnalysis
constraints Transformation
# constraints Plain
numberer Plain
# numberer RCM
system UmfPack
# system BandGeneral
test NormDispIncr 1.0e-8 50 0
algorithm Newton
# algorithm KrylovNewton
integrator LoadControl $dLambda
analysis Static

puts [format "----- gravity  stage %d  %s  %dx%.3g  lock=%d -----" \
	$soilMatStage $soilBoundary $nStep $dLambda $nGravLock]
set ok [analyze $nStep]
if {$ok != 0} {
	error "SoilGravity.tcl: gravity application (elastic) failed (ok=$ok)"
}

wipeAnalysis
loadConst -time 0.0

# Continuum -> stage 1; Liq stay 0 until ActivateEQBoundary
set updateLiqSpringStage 0
set soilMatStageWanted 1
source [file join $soilDir UpdateSoilStage.tcl]

wipeAnalysis
constraints Transformation
# constraints Plain
numberer Plain
# numberer RCM
system UmfPack
# system BandGeneral
test NormDispIncr 5.0e-8 100 0
algorithm KrylovNewton
# algorithm Newton
integrator LoadControl $dLambda
analysis Static
puts [format "----- gravity  stage %d  %s  %dx%.3g -----" \
	$soilMatStage $soilConstitutive $nStep $dLambda]
set ok [analyze $nStep]
if {$ok != 0} {
	error "SoilGravity.tcl: gravity application (plastic) failed (ok=$ok)"
}

wipeAnalysis
loadConst -time 0.0
