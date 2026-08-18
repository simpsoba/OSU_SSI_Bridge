# BuildPileSection.tcl
# Units: N, m, s
#
# Goals: pipe section for pileEleType.
#   elasticBeamColumn -- A_pile, Es_pile, I_pile already in Parameters (x n_pile_row)
#   dispBeamColumn    -- Fiber: graded tube strips; area x n_pile_row
# Expects Parameters.tcl and a model domain.

if {![info exists pileEleType] || ![info exists A_pile] || ![info exists I_pile]} {
	error "BuildPileSection.tcl: source Parameters.tcl first"
}

if {$pileEleType eq "elasticBeamColumn"} {

	# A, E, I used directly by BuildPilesNodes.tcl

} elseif {$pileEleType eq "dispBeamColumn"} {

	source [file join [file dirname [info script]] CircleStripFibers.tcl]

	# uniaxialMaterial Steel01 $matTag $Fy $E0 $b
	uniaxialMaterial Steel01 $matTag_steel_pile $fy_pile $Es_pile $b_pile

	set pileFibers [circularTubeFiberStripsGraded \
		$Ro_pile $Ri_pile $nFiberY_pile $nFiberEdge_pile]

	# section Fiber $secTag { fiber $y $z $A $matTag ... }
	section Fiber $secTag_pile {
		foreach fiberData $pileFibers {
			lassign $fiberData yLoc zLoc Af
			fiber $yLoc $zLoc [expr {$n_pile_row*$Af}] $matTag_steel_pile
		}
	}

} else {
	error "BuildPileSection.tcl: unknown pileEleType '$pileEleType' (elasticBeamColumn|dispBeamColumn)"
}
