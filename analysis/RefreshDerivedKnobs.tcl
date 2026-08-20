# RefreshDerivedKnobs.tcl
# Units: N, m, s
#
# Goals: after Overrides.tcl (or any late change to soilMesh / DT_FACTOR),
# rebuild the derived mesh bands and analysis step so they match the knobs.
# Requires: foot, cylinderSF already set (Parameters.tcl).

# ---
# 1. soilDxBands + L_half from soilMesh
# ---
if {![info exists soilMesh]} {
	error "RefreshDerivedKnobs.tcl: soilMesh missing"
}
if {$soilMesh == 0} {
	# production: 3 ft SSI to 12 ft
	set soilDxBands [list \
		[list [expr { 3.0*$foot}] [expr { 12.0*$foot}]] \
		[list [expr { 7.0*$foot}] [expr { 40.0*$foot}]] \
		[list [expr {15.0*$foot}] [expr {100.0*$foot}]] \
		[list [expr {20.0*$foot}] [expr {140.0*$foot}]] \
		[list [expr {20.0*$foot}] [expr {200.0*$foot}]] \
		]
} elseif {$soilMesh == 1} {
	# fine: 3 ft bands to 201 ft NF
	set soilDxBands [list \
		[list [expr { 3.0*$foot}] [expr { 12.0*$foot}]] \
		[list [expr { 3.0*$foot}] [expr { 39.0*$foot}]] \
		[list [expr { 3.0*$foot}] [expr { 99.0*$foot}]] \
		[list [expr { 3.0*$foot}] [expr {141.0*$foot}]] \
		[list [expr { 3.0*$foot}] [expr {201.0*$foot}]] \
		]
} elseif {$soilMesh == 2} {
	# finer: 3 ft bands to 270 ft NF
	set soilDxBands [list \
		[list [expr { 3.0*$foot}] [expr { 12.0*$foot}]] \
		[list [expr { 3.0*$foot}] [expr { 39.0*$foot}]] \
		[list [expr { 3.0*$foot}] [expr { 99.0*$foot}]] \
		[list [expr { 3.0*$foot}] [expr {141.0*$foot}]] \
		[list [expr { 3.0*$foot}] [expr {270.0*$foot}]] \
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
	error "RefreshDerivedKnobs.tcl: soilMesh must be -2, -1, 0, 1, or 2 (got '$soilMesh')"
}
set L_half [lindex [lindex $soilDxBands end] 1]

# ---
# 2. dtAnalysis from DT_FACTOR
# ---
if {![info exists DT_FACTOR]} {
	error "RefreshDerivedKnobs.tcl: DT_FACTOR missing"
}
if {![info exists cylinderSF]} {
	error "RefreshDerivedKnobs.tcl: cylinderSF missing"
}
set dtAnalysis [expr {$DT_FACTOR/2048.0*sqrt($cylinderSF)}]
