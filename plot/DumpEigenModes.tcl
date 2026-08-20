# DumpEigenModes.tcl
# Units: N, m, s
#
# Goals: after eigen, write mode shapes (nodeEigenvector) to JSON for PlotEigenModes.py.
# Expects eigenLambdas / nModesEigen. Optional: set eigenOutPath first.
# nodeEigenvector $nodeTag $mode <$dof>  (mode and dof 1-based)

if {![info exists eigenLambdas]} {
	error "DumpEigenModes.tcl: run eigen first (need eigenLambdas)"
}
if {![info exists nModesEigen]} {
	set nModesEigen [llength $eigenLambdas]
}

set scriptDir [file dirname [file normalize [info script]]]
if {![info exists root]} {
	set root [file dirname $scriptDir]
}
if {![info exists plotDir]} {
	set plotDir [file join $root plot]
}
if {![info exists eigenOutPath]} {
	if {[info exists soilProfile] && [info exists soilBoundary]} {
		if {![info exists soilEleType]} { set soilEleType "quad" }
		set eigDir [file join $plotDir out profile$soilProfile elevation $soilBoundary \
			modes $soilEleType $pierEleType]
	} else {
		set eigDir [file join $plotDir out modes $pierEleType]
	}
	file mkdir $eigDir
	set eigenOutPath [file join $eigDir eigen_modes.json]
}

proc eigNum {x} {
	return [format "%.8g" [expr {double($x)}]]
}

proc eigClassifyEle {e ni nj} {
	set g [modelEleGroup $e]
	if {$g ne ""} { return $g }
	set ci [nodeCoord $ni]
	set cj [nodeCoord $nj]
	set dx [expr {[lindex $ci 0] - [lindex $cj 0]}]
	set dy [expr {[lindex $ci 1] - [lindex $cj 1]}]
	if {[expr {sqrt($dx*$dx + $dy*$dy)}] < 1.0e-9} {
		return "spring"
	}
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

set outFd [open $eigenOutPath w]
puts $outFd "\{"
puts $outFd "  \"units\": \"m\","
puts $outFd "  \"title\": \"OSU SSI Bridge -- eigenmodes\","
puts $outFd [format "  \"pierEleType\": \"%s\"," $pierEleType]
if {[info exists soilEleType]} {
	puts $outFd [format "  \"soilEleType\": \"%s\"," $soilEleType]
} else {
	puts $outFd "  \"soilEleType\": \"quad\","
}
if {[info exists soilMatStage]} {
	puts $outFd [format "  \"soilMatStage\": %d," $soilMatStage]
} else {
	puts $outFd "  \"soilMatStage\": null,"
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
puts $outFd [format "  \"pierCrackedFactor\": %s," [eigNum $pierCrackedFactor]]
puts $outFd [format "  \"nModes\": %d," $nModesEigen]

# Periods
set piVal 3.141592653589793
puts $outFd "  \"modes_meta\": \["
set iMode 0
foreach lam $eigenLambdas {
	incr iMode
	set jsonComma [expr {($iMode < $nModesEigen) ? "," : ""}]
	if {$lam <= 0.0} {
		puts $outFd [format "    \{\"mode\": %d, \"lambda\": %s, \"T\": null, \"f\": null\}%s" \
			$iMode [eigNum $lam] $jsonComma]
	} else {
		set w [expr {sqrt($lam)}]
		set T [expr {2.0*$piVal/$w}]
		set f [expr {1.0/$T}]
		puts $outFd [format "    \{\"mode\": %d, \"lambda\": %s, \"T\": %s, \"f\": %s\}%s" \
			$iMode [eigNum $lam] [eigNum $T] [eigNum $f] $jsonComma]
	}
}
puts $outFd "  \],"

# Nodes: [tag, x, y]
set nTags [lsort -integer [getNodeTags]]
set nNode [llength $nTags]
puts $outFd "  \"nodes\": \["
set idx 0
foreach n $nTags {
	incr idx
	set xy [nodeCoord $n]
	set jsonComma [expr {($idx < $nNode) ? "," : ""}]
	puts $outFd [format "    \[%d, %s, %s\]%s" \
		$n [eigNum [lindex $xy 0]] [eigNum [lindex $xy 1]] $jsonComma]
}
puts $outFd "  \],"

# 2-node elements (structure + springs); skip soil continuum quads here
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
	set grp [eigClassifyEle $e $ni $nj]
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

# Soil quads as node tags [n1,n2,n3,n4] for deformed mesh
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

# ASDEA absorbing boundary quads (outer ring; node tags -> deformed via phi)
puts $outFd "  \"bnd_quads\": \["
set bndQuads {}
if {[info exists soilEleBndTags]} {
	foreach e $soilEleBndTags {
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
puts $outFd "  \],"

# Mode shapes: phi[mode-1] = [[tag, ux, uy], ...]
puts $outFd "  \"phi\": \["
set iMode 0
foreach lam $eigenLambdas {
	incr iMode
	puts $outFd "    \["
	set idx 0
	foreach n $nTags {
		incr idx
		# nodeEigenvector $tag $mode $dof -- 1-based mode and dof
		set ux 0.0
		set uy 0.0
		if {[catch {set ux [nodeEigenvector $n $iMode 1]}]} {
			set ux 0.0
		}
		if {[catch {set uy [nodeEigenvector $n $iMode 2]}]} {
			set uy 0.0
		}
		set jsonComma [expr {($idx < $nNode) ? "," : ""}]
		puts $outFd [format "      \[%d, %s, %s\]%s" \
			$n [eigNum $ux] [eigNum $uy] $jsonComma]
	}
	set modeComma [expr {($iMode < $nModesEigen) ? "," : ""}]
	puts $outFd "    \]$modeComma"
}
puts $outFd "  \]"

puts $outFd "\}"
close $outFd

# Convenience copy next to other plot JSON
file copy -force $eigenOutPath [file join $plotDir eigen_modes.json]

set dumpEigenModesDone 1
