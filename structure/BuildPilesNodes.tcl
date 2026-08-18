# BuildPilesNodes.tcl
# Units: N, m, s
#
# Goals: three steel-pipe shafts under the pile cap. Heads = cap bottom
# (1027, 1028, 1029). L_pile in nSeg_pile segments, downward.
#
# pileEleType (Parameters.tcl):
#   elasticBeamColumn -- A, E, I from Parameters (already x n_pile_row)
#   dispBeamColumn    -- Fiber sec from BuildPileSection.tcl
#
# Mass: nodal only (no element -mass).
#   Translation: half segment on tip and head; full on mids.
#   Rotation:    Irot = rhoL*dy^3/105 per segment end (mids get 2x).
#   Heads = cap bottom nodes; pile half-seg via IncrMass (keeps cap lump).
#
# Expects: Parameters, BuildPileCapNodes, BuildPileSection already sourced.
# Does not apply tip fixity (soil / PlotModel / a caller may).
# Shafts: BuildPilesElements.tcl after nodes.

if {![info exists L_pile] || ![info exists nSeg_pile] || ![info exists pileEleType]} {
	error "BuildPilesNodes.tcl: source Parameters.tcl first"
}
if {![info exists structNodeTags]} { set structNodeTags {} }
if {![info exists nodeTag_cap_BL] || ![info exists H_cap] || ![info exists s_pile_cap]} {
	error "BuildPilesNodes.tcl: source structure/BuildPileCapNodes.tcl first"
}

set scriptDir [file dirname [file normalize [info script]]]
source [file join $scriptDir IncrMass.tcl]

if {$pileEleType eq "dispBeamColumn"} {
	if {![info exists secTag_pile]} {
		error "BuildPilesNodes.tcl: source structure/BuildPileSection.tcl first"
	}
}
if {$nSeg_pile < 1} {
	error "BuildPilesNodes.tcl: nSeg_pile must be >= 1"
}

set dy [expr {$L_pile/double($nSeg_pile)}];           # m
set m_pile_seg  [expr {$rhoL_pile*$dy}];              # kg
set m_pile_half [expr {0.5*$m_pile_seg}];             # kg
set Irot_seg    [expr {$rhoL_pile*pow($dy,3)/105.0}]; # kg*m^2

set yHead [expr {-$H_cap}];  # m, cap bottom
set s $s_pile_cap

# Head tags and x: left, center, right
# tag (below head) = tagShift_pile + ip*100 + iy  (iy=1..nSeg; heads = cap BL/BC/BR)
set pileHeads [list \
	[list $nodeTag_cap_BL [expr {-$s}]] \
	[list $nodeTag_cap_BC 0.0] \
	[list $nodeTag_cap_BR $s] \
]

# =====================================================================
# 2. MODEL BUILDER / NODES
# =====================================================================
set tipTags {}

for {set ip 0} {$ip < $n_pile} {incr ip} {
	lassign [lindex $pileHeads $ip] headTag xP

	# Head: stack half-segment + Irot onto existing cap mass
	IncrMass $headTag $m_pile_half $m_pile_half $Irot_seg

	# Nodes below head: iy = 1 .. nSeg (tip = nSeg)
	# tag = tagShift_pile + ip*100 + iy  (heads = cap BL/BC/BR, IncrMass only)
	# node $tag $x $y -mass $mx $my $mRz
	for {set iy 1} {$iy <= $nSeg_pile} {incr iy} {
		set nTag [expr {$nodeTag_pile_base + $ip*100 + $iy}]
		set y [expr {$yHead - $iy*$dy}]

		if {$iy == $nSeg_pile} {
			# tip
			node $nTag $xP $y \
				-mass $m_pile_half $m_pile_half $Irot_seg
			lappend tipTags $nTag
		} else {
			# mid: two adjoining segments
			node $nTag $xP $y \
				-mass $m_pile_seg $m_pile_seg [expr {2.0*$Irot_seg}]
		}
		lappend structNodeTags $nTag
	}
}

set nodeTag_pile_tips $tipTags
