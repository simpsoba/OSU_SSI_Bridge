# soil/SoilDxBands.tcl
# Units: N, m, s
#
# Goals: map soilMesh (integer) to soilDxBands + L_half.
# Sourced from Parameters.tcl. RefreshDerivedKnobs calls ApplySoilDxBands
# after set soilMesh and again after Overrides.tcl.
#
# Each band row: {mesh size    x end}
# mesh size -- horizontal quad width in this ring
# x end     -- outer |x| where this width stops (next row starts there)
# x end = previous x end + n * mesh size (integer n), or the last cell is skinny.
# Last x end is L_half (near-field outer face).

# Fill soilDxBands (m) and L_half from the current soilMesh knob.
# Args: none (reads soilMesh, foot)
# Returns: none (sets soilDxBands, L_half)
proc ApplySoilDxBands {} {
	global soilMesh foot soilDxBands L_half
	if {![info exists soilMesh]} {
		error "SoilDxBands.tcl: soilMesh missing"
	}
	if {![info exists foot]} {
		error "SoilDxBands.tcl: foot missing (source Parameters.tcl first)"
	}
	if {$soilMesh == 0} {
		# production (2026-08-19 tag): 3 ft SSI to 12 ft; outer 30 ft to 200 ft
		set soilDxBands [list \
			[list [expr { 3.0*$foot}] [expr { 12.0*$foot}]] \
			[list [expr { 7.0*$foot}] [expr { 40.0*$foot}]] \
			[list [expr {15.0*$foot}] [expr {100.0*$foot}]] \
			[list [expr {20.0*$foot}] [expr {140.0*$foot}]] \
			[list [expr {30.0*$foot}] [expr {200.0*$foot}]] \
			]
	} elseif {$soilMesh == 1} {
		# moderate SSI: 3 ft to 39 ft (13 cells), 7 ft to 95 ft, then graded to 200 ft
		set soilDxBands [list \
			[list [expr { 3.0*$foot}] [expr { 39.0*$foot}]] \
			[list [expr { 7.0*$foot}] [expr { 95.0*$foot}]] \
			[list [expr {15.0*$foot}] [expr {140.0*$foot}]] \
			[list [expr {20.0*$foot}] [expr {200.0*$foot}]] \
			]
	} elseif {$soilMesh == 2} {
		# large SSI: each size one ring out; drop 20 ft (~83 x-stations)
		set soilDxBands [list \
			[list [expr { 3.0*$foot}] [expr { 84.0*$foot}]] \
			[list [expr { 7.0*$foot}] [expr {140.0*$foot}]] \
			[list [expr {15.0*$foot}] [expr {200.0*$foot}]] \
			]
	} elseif {$soilMesh == 3} {
		# x-large SSI: 3 ft to 114 ft, 7 ft to 170 ft, 15 ft to 200 ft (~99 x-stations)
		set soilDxBands [list \
			[list [expr { 3.0*$foot}] [expr {114.0*$foot}]] \
			[list [expr { 7.0*$foot}] [expr {170.0*$foot}]] \
			[list [expr {15.0*$foot}] [expr {200.0*$foot}]] \
			]
	} elseif {$soilMesh == 4} {
		# xx-large SSI: drop 15 ft; 3 ft to 123 ft, 7 ft to 200 ft (~107 x-stations)
		set soilDxBands [list \
			[list [expr { 3.0*$foot}] [expr {123.0*$foot}]] \
			[list [expr { 7.0*$foot}] [expr {200.0*$foot}]] \
			]
	} elseif {$soilMesh == -1} {
		# coarse: same inner SSI; ~25 x-stations
		set soilDxBands [list \
			[list [expr { 3.0*$foot}] [expr { 12.0*$foot}]] \
			[list [expr {14.0*$foot}] [expr { 40.0*$foot}]] \
			[list [expr {20.0*$foot}] [expr {100.0*$foot}]] \
			[list [expr {50.0*$foot}] [expr {200.0*$foot}]] \
			]
	} elseif {$soilMesh == -2} {
		# coarser: 3 ft only to outer pile (±s); ~19 x-stations
		set soilDxBands [list \
			[list [expr { 3.0*$foot}] [expr {  6.0*$foot}]] \
			[list [expr {12.0*$foot}] [expr { 30.0*$foot}]] \
			[list [expr {30.0*$foot}] [expr { 90.0*$foot}]] \
			[list [expr {55.0*$foot}] [expr {200.0*$foot}]] \
			]
	} else {
		error "SoilDxBands.tcl: soilMesh must be -2, -1, 0, 1, 2, 3, or 4 (got '$soilMesh')"
	}
	set L_half [lindex [lindex $soilDxBands end] 1]
}
