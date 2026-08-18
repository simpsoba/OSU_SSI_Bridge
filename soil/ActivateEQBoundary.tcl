# soil/ActivateEQBoundary.tcl
# Goals: Lysmer / ASDEA stage 1 + PyLiq/TzLiq stage 1.
# Call after gravity loadConst. Do not call before SoilGravity.
# Knobs: Parameters.tcl. BuildVelSeries must exist (ASDEA also at mesh time).
#
# =====================================================================
# 5. BOUNDARY CONDITIONS / CONSTRAINTS
# =====================================================================

if {![info exists soilBoundary]} {
	error "ActivateEQBoundary.tcl: source Parameters.tcl first"
}

set scriptDir [file dirname [file normalize [info script]]]

# Ensure no live Transformation analysis while BCs change
wipeAnalysis

if {$soilBoundary eq "Shin"} {
	if {![info exists gmVelNPTS]} {
		source [file join [file dirname $scriptDir] analysis BuildVelSeries.tcl]
	}
	source [file join $scriptDir BuildShinLysmer.tcl]
} elseif {$soilBoundary eq "ASDEA"} {
	if {![info exists soilEleBndTags] || [llength $soilEleBndTags] < 1} {
		error "ActivateEQBoundary.tcl: no ASDEA elements (BuildASDEABoundary missing?)"
	}
	# setParameter -val $value -ele $eleTag1 ... $paramName
	setParameter -val 1 -ele {*}$soilEleBndTags stage
	set soilBndStage 1
	puts [format "----- ASDEA stage -> 1 (%d elements) -----" \
		[llength $soilEleBndTags]]
} else {
	error "ActivateEQBoundary.tcl: soilBoundary must be Shin or ASDEA (got $soilBoundary)"
}

# Liq springs: stage 1 only for EQ (same gate as ASDEA / Lysmer)
set updateLiqSpringStage 1
set soilMatStageWanted 1
# Continuum already stage 1 from mid-gravity; UpdateSoilStage re-asserts it
source [file join $scriptDir UpdateSoilStage.tcl]
