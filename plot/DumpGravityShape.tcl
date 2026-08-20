# DumpGravityShape.tcl
# Units: N, m, s
#
# Goals: after gravity (+ loadConst), write undeformed coords + nodeDisp
# for PlotGravityShape.py. Call before eigen / ActivateEQBoundary / EQ.
# Optional: set gravityShapeOutPath first.

if {![info exists H_pier]} {
	error "DumpGravityShape.tcl: model + Parameters required"
}

set scriptDir [file dirname [file normalize [info script]]]
if {![info exists gravityShapeOutPath]} {
	set gsDir $scriptDir
	if {[info exists plotDir] && [info exists soilProfile] && [info exists soilBoundary]} {
		set eleType "quad"
		if {[info exists soilEleType]} { set eleType $soilEleType }
		set gsDir [file join $plotDir out profile$soilProfile elevation \
			$soilBoundary gravity $eleType $pierEleType]
	}
	file mkdir $gsDir
	set gravityShapeOutPath [file join $gsDir gravity_shape.json]
}

proc gsNum {x} {
	return [format "%.8g" [expr {double($x)}]]
}

proc gsClassifyEle {e ni nj} {
	set g [modelEleGroup $e]
	if {$g ne ""} { return $g }
	set ci [nodeCoord $ni]
	set cj [nodeCoord $nj]
	set dx [expr {[lindex $ci 0] - [lindex $cj 0]}]
	set dy [expr {[lindex $ci 1] - [lindex $cj 1]}]
	if {[expr {sqrt($dx*$dx + $dy*$dy)}] < 1.0e-9} { return "spring" }
	if {[info exists ::eleTag_cap_base] && [info exists ::eleTag_pile_base]} {
		if {$e >= $::eleTag_cap_base && $e < $::eleTag_pile_base} {
			return "cap"
		}
	}
	set hi [expr {($ni > $nj) ? $ni : $nj}]
	set lo [expr {($ni < $nj) ? $ni : $nj}]
	if {$hi >= 3000 && $hi < 4000} { return "deck" }
	if {$hi >= 2000 && $hi < 3000} { return "pile" }
	if {$lo >= 1000 && $hi < 2000} { return "cap" }
	if {$ni < 1000 && $nj < 1000} { return "pier" }
	return "other"
}

set outFd [open $gravityShapeOutPath w]
puts $outFd "\{"
puts $outFd "  \"units\": \"m\","
puts $outFd "  \"title\": \"OSU SSI Bridge -- post-gravity\","
puts $outFd [format "  \"pierEleType\": \"%s\"," $pierEleType]
if {[info exists soilEleType]} {
	puts $outFd [format "  \"soilEleType\": \"%s\"," $soilEleType]
} else {
	puts $outFd "  \"soilEleType\": \"quad\","
}
if {[info exists soilProfile]} {
	puts $outFd [format "  \"soilProfile\": %d," $soilProfile]
} else {
	puts $outFd "  \"soilProfile\": null,"
}
if {[info exists soilBoundary]} {
	puts $outFd [format "  \"soilBoundary\": \"%s\"," $soilBoundary]
} else {
	puts $outFd "  \"soilBoundary\": null,"
}
puts $outFd [format "  \"time\": %s," [gsNum [getTime]]]

# Nodes: [tag, x0, y0, ux, uy]
set nTags [lsort -integer [getNodeTags]]
set nNode [llength $nTags]
puts $outFd "  \"nodes\": \["
set idx 0
foreach n $nTags {
	incr idx
	set xy [nodeCoord $n]
	set ux 0.0
	set uy 0.0
	catch { set ux [nodeDisp $n 1] }
	catch { set uy [nodeDisp $n 2] }
	set jsonComma [expr {($idx < $nNode) ? "," : ""}]
	puts $outFd [format "    \[%d, %s, %s, %s, %s\]%s" \
		$n [gsNum [lindex $xy 0]] [gsNum [lindex $xy 1]] \
		[gsNum $ux] [gsNum $uy] $jsonComma]
}
puts $outFd "  \],"

puts $outFd "  \"elements\": \["
set beamEles {}
foreach e [lsort -integer [getEleTags]] {
	if {[info exists soilEleTags] && [lsearch -exact $soilEleTags $e] >= 0} {
		continue
	}
	set en [eleNodes $e]
	if {[llength $en] < 2} { continue }
	set ni [lindex $en 0]
	set nj [lindex $en 1]
	set grp [gsClassifyEle $e $ni $nj]
	if {$grp eq "soil" || $grp eq "soil_bnd"} { continue }
	lappend beamEles [list $e $ni $nj $grp]
}
set nEle [llength $beamEles]
set idx 0
foreach row $beamEles {
	incr idx
	lassign $row e ni nj grp
	set jsonComma [expr {($idx < $nEle) ? "," : ""}]
	puts $outFd [format "    \[%d, %d, %d, \"%s\"]%s" $e $ni $nj $grp $jsonComma]
}
puts $outFd "  \],"

puts $outFd "  \"soil_quads\": \["
set soilQuads {}
if {[info exists soilEleTags]} {
	foreach e $soilEleTags {
		set en [eleNodes $e]
		if {[llength $en] < 4} { continue }
		lappend soilQuads [lrange $en 0 3]
	}
}
set nSoilQuad [llength $soilQuads]
set idx 0
foreach q $soilQuads {
	incr idx
	set jsonComma [expr {($idx < $nSoilQuad) ? "," : ""}]
	puts $outFd [format "    \[%d, %d, %d, %d\]%s" \
		[lindex $q 0] [lindex $q 1] [lindex $q 2] [lindex $q 3] $jsonComma]
}
puts $outFd "  \],"

puts $outFd "  \"bnd_quads\": \["
set bndQuads {}
if {[info exists asdeaEleTags]} {
	foreach e $asdeaEleTags {
		set en [eleNodes $e]
		if {[llength $en] < 4} { continue }
		lappend bndQuads [lrange $en 0 3]
	}
}
set nBndQuad [llength $bndQuads]
set idx 0
foreach q $bndQuads {
	incr idx
	set jsonComma [expr {($idx < $nBndQuad) ? "," : ""}]
	puts $outFd [format "    \[%d, %d, %d, %d\]%s" \
		[lindex $q 0] [lindex $q 1] [lindex $q 2] [lindex $q 3] $jsonComma]
}
puts $outFd "  \]"
puts $outFd "\}"
close $outFd

