# analysis/FoldStructNodes.tcl
# Goals:
#   Drop soil-gravity equalDOF, move pile/cap/pier/deck coords onto the
#   current settlement, create pier/deck/cap/pile beam-columns, then
#   re-tie pile/cap -> spring dups.
# Call after SoilGravity.tcl.

if {![info exists structDir]} {
	error "FoldStructNodes.tcl: set structDir first"
}

# remove mp on the constrained node of each gravity equalDOF
foreach nd $gravMpRemove {
	remove mp $nd
}

# Drop pile/cap -> dup ties so fold does not move the spring soil-side node
foreach pr $springEqualDOFPairs {
	lassign $pr retainNd dupNd soilNd
	remove mp $dupNd
}

foreach n $gravStructFixRZ {
	remove sp $n 3
}

# Move each unique structure node to undeformed + displacement, then zero disp
set nFold 0
array unset foldDone
foreach n $gravFoldNodes {
	if {[info exists foldDone($n)]} { continue }
	set foldDone($n) 1
	set xy [nodeCoord $n]
	set x0 [lindex $xy 0]
	set y0 [lindex $xy 1]
	set ux [nodeDisp $n 1]
	set uy [nodeDisp $n 2]
	setNodeCoord $n 1 [expr {$x0 + $ux}]
	setNodeCoord $n 2 [expr {$y0 + $uy}]
	setNodeDisp $n 1 0.0 -commit
	setNodeDisp $n 2 0.0 -commit
	setNodeDisp $n 3 0.0 -commit
	incr nFold
}
source [file join $structDir BuildStructElements.tcl]

set nSprTie 0
foreach pr $springEqualDOFPairs {
	lassign $pr retainNd dupNd soilNd
	equalDOF $retainNd $dupNd 1 2
	incr nSprTie
}
puts [format "----- fold  %d nodes  restore spr MP %d -----" $nFold $nSprTie]
if {$nSprTie > 0} {
	gravPrintSpringKine "after fold"
}
