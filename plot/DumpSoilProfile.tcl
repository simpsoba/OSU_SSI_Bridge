# DumpSoilProfile.tcl
# After BuildSoilMaterials.tcl: write plot/soil_profile.json for PlotSoilProfile.py.
#
# Optional: set soilProfileOutPath before sourcing.

if {![info exists nSoilRows] || ![info exists soilGr]} {
	error "DumpSoilProfile.tcl: BuildSoilMaterials.tcl first"
}

set scriptDir [file dirname [file normalize [info script]]]
if {![info exists soilProfileOutPath]} {
	set soilProfileOutPath [file join $scriptDir soil_profile.json]
}

proc spNum {x} {
	return [format "%.8g" $x]
}

set outFd [open $soilProfileOutPath w]
puts $outFd "\{"
puts $outFd "  \"units\": \"N, m, s\","
puts $outFd [format "  \"soilProfile\": %d," $soilProfile]
puts $outFd [format "  \"gam_max\": %s," [spNum $gam_max]]
puts $outFd [format "  \"pRef_clay\": %s," [spNum $pRef_c]]
puts $outFd [format "  \"pRef_sand\": %s," [spNum $pRef_s]]
puts $outFd [format "  \"B_fsp\": %s," [spNum $B_fsp]]
puts $outFd "  \"layers\": \["
set nLayer $nSoilRows
for {set iy 0} {$iy < $nSoilRows} {incr iy} {
	set nm $soilRowLayer($iy)
	set yT [lindex $soilYs $iy]
	set yB [lindex $soilYs [expr {$iy + 1}]]
	set jsonComma [expr {($iy < $nSoilRows - 1) ? "," : ""}]
	set sand [expr {$soilIsSand($iy) ? "true" : "false"}]
	set fsp  [expr {$soilFSP($iy) ? "true" : "false"}]
	puts $outFd [format "    \{"]
	puts $outFd [format "      \"name\": \"%s\"," $nm]
	puts $outFd [format "      \"type\": \"%s\"," $soilType($iy)]
	puts $outFd [format "      \"sand\": %s," $sand]
	puts $outFd [format "      \"fsp\": %s," $fsp]
	puts $outFd [format "      \"yTop\": %s," [spNum $yT]]
	puts $outFd [format "      \"yBot\": %s," [spNum $yB]]
	puts $outFd [format "      \"depthTop\": %s," [spNum [expr {-$yT}]]]
	puts $outFd [format "      \"depthBot\": %s," [spNum [expr {-$yB}]]]
	puts $outFd [format "      \"rho\": %s," [spNum $soilRho($iy)]]
	puts $outFd [format "      \"Gr\": %s," [spNum $soilGr($iy)]]
	puts $outFd [format "      \"Br\": %s," [spNum $soilBr($iy)]]
	puts $outFd [format "      \"phi\": %s," [spNum $soilPhi($iy)]]
	puts $outFd [format "      \"c\": %s," [spNum $soilC($iy)]]
	puts $outFd [format "      \"Dr\": %s," [spNum $soilDr($iy)]]
	puts $outFd [format "      \"k_pci\": %s," [spNum $soilKpci($iy)]]
	puts $outFd [format "      \"PTA\": %s," [spNum $soilPTA($iy)]]
	puts $outFd [format "      \"contr1\": %s," [spNum $soilContr1($iy)]]
	puts $outFd [format "      \"contr3\": %s," [spNum $soilContr3($iy)]]
	puts $outFd [format "      \"dilat1\": %s," [spNum $soilDilat1($iy)]]
	puts $outFd [format "      \"dilat3\": %s," [spNum $soilDilat3($iy)]]
	puts $outFd [format "      \"d\": %s," [spNum $soilD($iy)]]
	puts $outFd [format "      \"gam_max\": %s," [spNum $soilGamMax($iy)]]
	puts $outFd [format "      \"pRef\": %s," [spNum $soilPRef($iy)]]
	puts $outFd [format "      \"nYS\": %s," [spNum $soilNYS($iy)]]
	puts $outFd [format "      \"contr2\": %s," [spNum $soilContr2($iy)]]
	puts $outFd [format "      \"dilat2\": %s," [spNum $soilDilat2($iy)]]
	puts $outFd [format "      \"liq1\": %s," [spNum $soilLiq1($iy)]]
	puts $outFd [format "      \"liq2\": %s," [spNum $soilLiq2($iy)]]
	puts $outFd [format "      \"e\": %s," [spNum $soilE($iy)]]
	puts $outFd [format "      \"cs1\": %s," [spNum $soilCs1($iy)]]
	puts $outFd [format "      \"cs2\": %s," [spNum $soilCs2($iy)]]
	puts $outFd [format "      \"cs3\": %s," [spNum $soilCs3($iy)]]
	puts $outFd [format "      \"pa\": %s," [spNum $soilPa($iy)]]
	puts $outFd [format "      \"B_fsp\": %s" [spNum $soilB_fsp($iy)]]
	puts $outFd [format "    \}%s" $jsonComma]
}
puts $outFd "  \]"
puts $outFd "\}"
close $outFd

