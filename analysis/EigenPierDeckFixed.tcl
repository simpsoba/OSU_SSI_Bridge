# analysis/EigenPierDeckFixed.tcl
# Goals: eigen of elastic pier + deck frame only (no soil/cap/piles).
# Fixed pier base UX+UY+RZ. Compare to hand T = 2 pi sqrt(M/(3EI/H^3)).
#
#   OpenSees analysis/EigenPierDeckFixed.tcl
#
# Units: N, m, s. Temp check only — does not write recorders.

wipe

set root [file dirname [file dirname [file normalize [info script]]]]
set structDir [file join $root structure]
set analysisDir [file join $root analysis]
source [file join $root Parameters.tcl]

# Force elastic pier for this check (ignore Parameters lumpedPlasticity)
set pierEleType "elasticBeamColumn"

model BasicBuilder -ndm 2 -ndf 3
set structNodeTags {}
source [file join $structDir PierSection.tcl]
source [file join $structDir BuildPierNodes.tcl]
source [file join $structDir BuildDeckNodes.tcl]
source [file join $structDir BuildPierElements.tcl]
source [file join $structDir BuildDeckElements.tcl]

# Fixed base (true cantilever)
fix $nodeTag_pierBase_capTC 1 1 1

# Optional: gravity then PDelta eigen (matches campaign geomTransf)
timeSeries Linear 1
pattern Plain 1 1 {
	foreach n $structNodeTags {
		set m [nodeMass $n 1]
		if {$m > 0.0} {
			load $n 0.0 [expr {-$m*$gravity_accel}] 0.0
		}
	}
}
constraints Transformation
numberer Plain
system UmfPack
algorithm Linear
integrator LoadControl 0.1
analysis Static
analyze 10
loadConst -time 0.0

wipeAnalysis
constraints Transformation
system UmfPack
set nModes 8
set lambdas [eigen -fullGenLapack $nModes]

set pi 3.141592653589793
set Mhand [expr {$m_deck + 0.5*$rhoL_pier*$H_pier}]
set khand [expr {3.0*$Ec_pier*$I_pier/pow($H_pier,3)}]
set Thand [expr {2.0*$pi*sqrt($Mhand/$khand)}]

puts "----- Pier + deck only (elasticBeamColumn, base fixed) -----"
puts [format "  hand SDOF  M=%.3e kg  k=%.3e N/m  T=%.5f s" $Mhand $khand $Thand]
puts [format "  I_pier=%.5f m4  H=%.4f m  Ec=%.4e Pa  m_deck=%.3e kg" \
	$I_pier $H_pier $Ec_pier $m_deck]
puts "| mode | T (s) | f (Hz) |"
set i 1
foreach lam $lambdas {
	if {$lam <= 0.0} {
		puts [format "| %4d | -- | -- |" $i]
	} else {
		set w [expr {sqrt($lam)}]
		set T [expr {2.0*$pi/$w}]
		set f [expr {1.0/$T}]
		puts [format "| %4d | %8.5f | %8.5f |" $i $T $f]
	}
	incr i
}
