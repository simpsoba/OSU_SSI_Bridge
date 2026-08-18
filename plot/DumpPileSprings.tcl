# DumpPileSprings.tcl
# After BuildSoilSprings.tcl: write plot/pile_springs.json for PlotPileSprings.py.
#
# Optional: set pileSpringOutPath before sourcing.

if {![info exists pileSpringPropsDump]} {
	error "DumpPileSprings.tcl: BuildSoilSprings.tcl first (need pileSpringPropsDump)"
}

set scriptDir [file dirname [file normalize [info script]]]
if {![info exists pileSpringOutPath]} {
	set pileSpringOutPath [file join $scriptDir pile_springs.json]
}

proc psNum {x} {
	return [format "%.8g" $x]
}

set outFd [open $pileSpringOutPath w]
puts $outFd "\{"
puts $outFd "  \"units\": \"N, m, s\","
puts $outFd [format "  \"soilProfile\": %d," $soilProfile]
if {[info exists pileSpring]} {
	puts $outFd [format "  \"pileSpring\": \"%s\"," $pileSpring]
}
puts $outFd [format "  \"pRes_frac\": %s," [psNum $pRes_frac]]
puts $outFd [format "  \"Ge_pile\": %s," [psNum $Ge_pile]]
puts $outFd [format "  \"n_pile_row\": %d," $n_pile_row]
puts $outFd [format "  \"D_pile\": %s," [psNum $D_pile]]
puts $outFd [format "  \"H_cap\": %s," [psNum $H_cap]]
puts $outFd [format "  \"W_cap\": %s," [psNum $W_cap]]
puts $outFd [format "  \"L_cap\": %s," [psNum $L_cap]]
if {[info exists PultCap]} {
	puts $outFd [format "  \"PultCap\": %s," [psNum $PultCap]]
	puts $outFd [format "  \"TultCap\": %s," [psNum $TultCap]]
	puts $outFd [format "  \"QultSoffit\": %s," [psNum $QultSoffit]]
}
puts $outFd [format "  \"L_pile\": %s," [psNum $L_pile]]
puts $outFd [format "  \"z50_tz\": %s," [psNum $z50_tz]]

# Layer bands (elevation y)
puts $outFd "  \"layers\": \["
set nLayer [llength $soilLayerNames]
set idx 0
foreach nm $soilLayerNames {
	incr idx
	set jsonComma [expr {($idx < $nLayer) ? "," : ""}]
	set sand [expr {$soilIsSand($nm) ? "true" : "false"}]
	puts $outFd [format "    \{\"name\": \"%s\", \"yTop\": %s, \"yBot\": %s, \"sand\": %s\}%s" \
		$nm [psNum $soilYTop($nm)] [psNum $soilYBot($nm)] $sand $jsonComma]
}
puts $outFd "  \],"

# Stations: one row per pile spring side
# ip iy y depth layer isTip useLiq pult pRes tAx tRes y50 z50 latType axType trib [side]
puts $outFd "  \"stations\": \["
set nSpring [llength $pileSpringPropsDump]
set idx 0
foreach rec $pileSpringPropsDump {
	incr idx
	set side 0.0
	if {[llength $rec] >= 17} {
		lassign $rec ip iyP y depth nm isTip useLiq pult pRes tAx tRes y50 z50Ax latType axType trib side
	} else {
		lassign $rec ip iyP y depth nm isTip useLiq pult pRes tAx tRes y50 z50Ax latType axType trib
	}
	set jsonComma [expr {($idx < $nSpring) ? "," : ""}]
	puts $outFd [format "    \{\"ip\": %d, \"iy\": %d, \"side\": %s, \"y\": %s, \"depth\": %s, \"layer\": \"%s\", \"isTip\": %d, \"useLiq\": %d, \"pult\": %s, \"pRes\": %s, \"tult\": %s, \"tRes\": %s, \"y50\": %s, \"z50\": %s, \"latType\": \"%s\", \"axType\": \"%s\", \"trib\": %s\}%s" \
		$ip $iyP [psNum $side] [psNum $y] [psNum $depth] $nm $isTip $useLiq \
		[psNum $pult] [psNum $pRes] [psNum $tAx] [psNum $tRes] \
		[psNum $y50] [psNum $z50Ax] $latType $axType [psNum $trib] $jsonComma]
}
puts $outFd "  ],"

puts $outFd "  \"y50_cap\": [format %.8g $y50_cap],"
puts $outFd "  \"z50_cap\": [format %.8g $z50_cap],"
if {![info exists capFacePropsDump]} { set capFacePropsDump {} }
if {![info exists capSoffitPropsDump]} { set capSoffitPropsDump {} }
puts $outFd "  \"cap_face\": \["
set nCap [llength $capFacePropsDump]
set idx 0
foreach rec $capFacePropsDump {
	incr idx
	lassign $rec e x y pult tult y50 z50
	set jsonComma [expr {($idx < $nCap) ? "," : ""}]
	puts $outFd [format "    \{\"ele\": %d, \"x\": %s, \"y\": %s, \"pult\": %s, \"tult\": %s, \"y50\": %s, \"z50\": %s\}%s" \
		$e [psNum $x] [psNum $y] [psNum $pult] [psNum $tult] \
		[psNum $y50] [psNum $z50] $jsonComma]
}
puts $outFd "  ],"
puts $outFd "  \"cap_soffit\": \["
set nSoffitOut [llength $capSoffitPropsDump]
set idx 0
foreach rec $capSoffitPropsDump {
	incr idx
	lassign $rec e x y qult z50
	set jsonComma [expr {($idx < $nSoffitOut) ? "," : ""}]
	puts $outFd [format "    \{\"ele\": %d, \"x\": %s, \"y\": %s, \"qult\": %s, \"z50\": %s\}%s" \
		$e [psNum $x] [psNum $y] [psNum $qult] [psNum $z50] $jsonComma]
}
puts $outFd "  ]"
puts $outFd "\}"
close $outFd

