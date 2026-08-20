# soil/TagHelpers.tcl
# Soil / pile tag formulas + small helpers so dense meshes do not collide.
# Sourced from Parameters.tcl. BuildSoilMesh sets soilNodeStride = nY.
# Spring / boundary bases move to clear existing tags when those builders run.

# Next thousand (or align) strictly above `after`.
proc nextTagBlock {after {align 1000}} {
	expr {int((($after / $align) + 1) * $align)}
}

# If $varName <= $after, set it to nextTagBlock($after). Returns 1 if moved.
proc ensureAbove {varName after} {
	upvar 1 $varName v
	if {$v <= $after} {
		set v [nextTagBlock $after]
		return 1
	}
	return 0
}

# Soil continuum node. ix, iy 0-based. Stride = nY (set in BuildSoilMesh).
proc soilNodeTag {ix iy} {
	global nodeTag_soil_base soilNodeStride
	expr {$nodeTag_soil_base + $ix*$soilNodeStride + $iy}
}

# Pile shaft node below the head (iy = 1..nSeg). Heads are cap BL/BC/BR.
proc pileNodeTag {ip iy} {
	global nodeTag_pile_base pileNodeStride
	expr {$nodeTag_pile_base + $ip*$pileNodeStride + $iy}
}

# Continuum soil node (uses soilNodeLast from BuildSoilMesh).
proc isSoilContinuumNode {n} {
	global nodeTag_soil_base soilNodeLast
	if {![info exists soilNodeLast]} { return 0 }
	expr {$n >= $nodeTag_soil_base && $n <= $soilNodeLast}
}

# Element group from live tag ranges. Empty if not soil / spring / boundary.
proc modelEleGroup {e} {
	global eleTag_spr_base eleTag_spr_last eleTag_soil_base eleTag_soil_last \
		eleTag_bnd_base eleTag_bnd_last
	if {[info exists eleTag_spr_base]} {
		set hi $eleTag_spr_base
		if {[info exists eleTag_spr_last]} {
			set hi $eleTag_spr_last
		} elseif {[info exists eleTag_bnd_base]} {
			set hi [expr {$eleTag_bnd_base - 1}]
		} else {
			set hi [expr {$eleTag_spr_base + 99999}]
		}
		if {$e >= $eleTag_spr_base && $e <= $hi} {
			return "ssi_spring"
		}
	}
	if {[info exists eleTag_bnd_base]} {
		set hi $eleTag_bnd_base
		if {[info exists eleTag_bnd_last]} {
			set hi $eleTag_bnd_last
		} else {
			set hi [expr {$eleTag_bnd_base + 99999}]
		}
		if {$e >= $eleTag_bnd_base && $e <= $hi} {
			return "soil_bnd"
		}
	}
	if {[info exists eleTag_soil_base]} {
		set hi $eleTag_soil_base
		if {[info exists eleTag_soil_last]} {
			set hi $eleTag_soil_last
		} elseif {[info exists eleTag_spr_base]} {
			set hi [expr {$eleTag_spr_base - 1}]
		} else {
			set hi [expr {$eleTag_soil_base + 99999}]
		}
		if {$e >= $eleTag_soil_base && $e <= $hi} {
			return "soil"
		}
	}
	return ""
}
