# IncrMass.tcl
# Units: N, m, s
#
# Goals: add to existing nodal mass without wiping what is already there.
# OpenSees `mass` replaces the diagonal; use this when stacking contributions
# on a shared node (deck soffit CL <-> pier top; cap TC <-> pier base;
# pile head <-> cap bottom).
#
# Args:    nodeTag  dmx ?dmy? ?dIrot?   (kg, kg.m2)
# Returns: none (updates the node)

proc IncrMass {nodeTag args} {
	set n [llength $args]
	if {$n < 1} {
		error "IncrMass: want nodeTag mx ?my? ?Irot? ..."
	}
	set existing {}
	for {set i 1} {$i <= $n} {incr i} {
		lappend existing [nodeMass $nodeTag $i]
	}
	set total {}
	for {set i 0} {$i < $n} {incr i} {
		lappend total [expr {[lindex $existing $i] + [lindex $args $i]}]
	}
	mass $nodeTag {*}$total
}
