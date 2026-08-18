# DumpFiberSections.tcl
# Units: N, m, s
#
# Goals: write plot/fiber_sections.json from fiber lists left by PierSection /
# BuildPileSection (same strip data registered in OpenSees).
# Plot: python3 plot/PlotFiberSections.py
# Expects Parameters + PierSection (+ BuildPileSection if piles are built).

set scriptDir [file dirname [file normalize [info script]]]
if {![info exists fiberOutPath]} {
	set fiberOutPath [file join $scriptDir fiber_sections.json]
}

proc fibNum {x} {
	return [format "%.8g" [expr {double($x)}]]
}

set outFd [open $fiberOutPath w]
puts $outFd "\{"
puts $outFd [format "  \"pierEleType\": \"%s\"," $pierEleType]
puts $outFd [format "  \"pileEleType\": \"%s\"," $pileEleType]
if {[info exists soilProfile]} {
	puts $outFd [format "  \"soilProfile\": %d," $soilProfile]
} else {
	puts $outFd "  \"soilProfile\": null,"
}

# ---- Pier Fiber (forceBeamColumn hinges / lumpedPlasticity ZLS) ----
set hasPier 0
if {[info exists coreFibers] && [info exists coverFibers] && [info exists rebarFibers]} {
	set hasPier 1
	puts $outFd "  \"pier\": \{"
	puts $outFd "    \"kind\": \"RC_circle\","
	puts $outFd [format "    \"R\": %s," [fibNum $R_pier]]
	puts $outFd [format "    \"R_core\": %s," [fibNum $R_core_pier]]
	puts $outFd [format "    \"R_bar\": %s," [fibNum $R_core_pier]]
	puts $outFd [format "    \"nFiberY\": %d," $nFiberY_pier]
	puts $outFd [format "    \"nFiberEdge\": %d," $nFiberEdge_pier]
	puts $outFd "    \"edgeFrac\": 0.1666666667,"
	puts $outFd "    \"core\": \["
	set nc [llength $coreFibers]
	set i 0
	foreach f $coreFibers {
		incr i
		lassign $f y z A
		set c [expr {($i < $nc) ? "," : ""}]
		puts $outFd [format "      \[%s, %s, %s\]%s" [fibNum $y] [fibNum $z] [fibNum $A] $c]
	}
	puts $outFd "    \],"
	puts $outFd "    \"cover\": \["
	set nc [llength $coverFibers]
	set i 0
	foreach f $coverFibers {
		incr i
		lassign $f y z A
		set c [expr {($i < $nc) ? "," : ""}]
		puts $outFd [format "      \[%s, %s, %s\]%s" [fibNum $y] [fibNum $z] [fibNum $A] $c]
	}
	puts $outFd "    \],"
	puts $outFd "    \"rebar\": \["
	set nc [llength $rebarFibers]
	set i 0
	foreach f $rebarFibers {
		incr i
		lassign $f y z A
		set c [expr {($i < $nc) ? "," : ""}]
		puts $outFd [format "      \[%s, %s, %s\]%s" [fibNum $y] [fibNum $z] [fibNum $A] $c]
	}
	puts $outFd "    \]"
	puts $outFd "  \},"
} else {
	puts $outFd "  \"pier\": null,"
}

# ---- Pile Fiber (dispBeamColumn tube strips) ----
if {[info exists pileFibers]} {
	puts $outFd "  \"pile\": \{"
	puts $outFd "    \"kind\": \"steel_tube\","
	puts $outFd [format "    \"Ro\": %s," [fibNum $Ro_pile]]
	puts $outFd [format "    \"Ri\": %s," [fibNum $Ri_pile]]
	puts $outFd [format "    \"nFiberY\": %d," $nFiberY_pile]
	puts $outFd [format "    \"nFiberEdge\": %d," $nFiberEdge_pile]
	puts $outFd "    \"edgeFrac\": 0.1666666667,"
	puts $outFd [format "    \"n_pile_row\": %d," $n_pile_row]
	puts $outFd "    \"steel\": \["
	set nc [llength $pileFibers]
	set i 0
	foreach f $pileFibers {
		incr i
		lassign $f y z A
		set c [expr {($i < $nc) ? "," : ""}]
		puts $outFd [format "      \[%s, %s, %s\]%s" \
			[fibNum $y] [fibNum $z] [fibNum [expr {$n_pile_row*$A}]] $c]
	}
	puts $outFd "    \]"
	puts $outFd "  \}"
} else {
	puts $outFd "  \"pile\": null"
}

puts $outFd "\}"
close $outFd

