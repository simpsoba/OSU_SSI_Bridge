# analysis/GravityHelpers.tcl
# Goals: print SSI spring kinematics; find the nearest soil node.
# Sourced from Run.tcl before SoilGravity.

# Max relative UY between pile/cap, dup, and soil. label is a stage name.
# Args:    label (string)
# Returns: none (puts)
proc gravPrintSpringKine {label} {
	global springEqualDOFPairs
	if {[llength $springEqualDOFPairs] < 1} {
		return
	}
	set maxRel 0.0
	set maxPile 0.0
	foreach pr $springEqualDOFPairs {
		lassign $pr retainNd dupNd soilNd
		set uyR [nodeDisp $retainNd 2]
		set uyD [nodeDisp $dupNd 2]
		set uyS [nodeDisp $soilNd 2]
		set rel [expr {abs($uyD - $uyS)}]
		if {$rel > $maxRel} { set maxRel $rel }
		set ap [expr {abs($uyR)}]
		if {$ap > $maxPile} { set maxPile $ap }
	}
	puts [format "  SSI kine (%s): max |uy_dup-uy_soil|=%.4e  max |uy_pile/cap|=%.4e  n=%d" \
		$label $maxRel $maxPile [llength $springEqualDOFPairs]]
}

# Nearest existing continuum node (tags in [nodeTag_soil_base, soilNodeLast]).
# Args:    xT yT (m)
# Returns: node tag, or -1 if none
proc soilNdNearestExisting {xT yT} {
	global nodeTag_soil_base soilNodeLast
	set best -1
	set bd 1.0e99
	foreach n [getNodeTags] {
		if {![isSoilContinuumNode $n]} { continue }
		set xy [nodeCoord $n]
		set d [expr {hypot([lindex $xy 0] - $xT, [lindex $xy 1] - $yT)}]
		if {$d < $bd} {
			set bd $d
			set best $n
		}
	}
	return $best
}
