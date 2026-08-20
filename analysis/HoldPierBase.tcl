# analysis/HoldPierBase.tcl
# Goals: after gravity loadConst, optionally freeze pier-base UX/UY at the
# current displacement. RZ stays free. Eigen/EQ then see a pin at the cap
# center (TC). Gated by holdPierON in Parameters.tcl (default on).
#
# sp $nodeTag $dof $disp -const -subtractInit  ->  hold gravity disp, not snap to 0

if {![info exists holdPierON]} {
	set holdPierON 1
}
if {!$holdPierON} {
	puts "----- Hold pier base skipped (holdPierON 0) -----"
	return
}

timeSeries Constant $tsTag_holdPier

if {$pierEleType eq "lumpedPlasticity"} {
	# pierBase_capTC / ZLS outer; pierBaseZeroLengthInner / stiff-column foot
	remove mp $nodeTag_pierBaseZeroLengthInner
	pattern Plain $patternTag_holdPier $tsTag_holdPier "
		sp $nodeTag_pierBase_capTC 1 0.0 -const -subtractInit
		sp $nodeTag_pierBase_capTC 2 0.0 -const -subtractInit
		sp $nodeTag_pierBaseZeroLengthInner 1 0.0 -const -subtractInit
		sp $nodeTag_pierBaseZeroLengthInner 2 0.0 -const -subtractInit
	"
} else {
	pattern Plain $patternTag_holdPier $tsTag_holdPier "
		sp $nodeTag_pierBase_capTC 1 0.0 -const -subtractInit
		sp $nodeTag_pierBase_capTC 2 0.0 -const -subtractInit
	"
}

puts [format "----- Hold pier base UX+UY (%s) -----" $pierEleType]
