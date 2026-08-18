# analysis/StructureGravityLoads.tcl
# Goals: nodal -mg on pier, cap, deck, piles (structNodeTags).
# Call after FoldStructNodes. Analyze / loadConst are in Run.tcl.
#
# =====================================================================
# 6. LOADS
# =====================================================================

if {![info exists gravity_accel]} {
	error "StructureGravityLoads.tcl: source Parameters.tcl first"
}
if {![info exists structNodeTags] || [llength $structNodeTags] < 1} {
	error "StructureGravityLoads.tcl: structNodeTags empty (BuildModel.tcl first)"
}

# timeSeries Linear $tag
timeSeries Linear $tsTag_gravStruct

set nLoad 0
set wSum 0.0
# pattern Plain $patternTag $tsTag { load $nodeTag $Fx $Fy $Mz }
# Y-weight from nodal mass already on pier / cap / deck / piles
pattern Plain $patternTag_gravStruct $tsTag_gravStruct {
	foreach n $structNodeTags {
		set m [nodeMass $n 2]
		if {$m <= 0.0} { continue }
		load $n 0.0 [expr {-$m*$gravity_accel}] 0.0
		incr nLoad
		set wSum [expr {$wSum + $m*$gravity_accel}]
	}
}

puts [format "----- Structure gravity  %d nodes  SumW=%.3e N -----" \
	$nLoad $wSum]
