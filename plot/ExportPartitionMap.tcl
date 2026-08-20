# plot/ExportPartitionMap.tcl
# Units: N, m, s
#
# Goals: after METIS `partition`, write this rank's elements/nodes, then
# (rank 0) merge + PNG. Call from RunParallel.tcl when exportPartitionMap=1.
# Plot:  python3 plot/PlotPartition.py $partitionOutDir
#
# Optional: set partitionOutDir before sourcing.

if {![info exists np] || ![info exists pid]} {
	error "ExportPartitionMap.tcl: need np / pid from OpenSeesMP"
}
if {![info exists H_pier]} {
	error "ExportPartitionMap.tcl: model + Parameters required"
}

set scriptDir [file dirname [file normalize [info script]]]
if {![info exists root]} {
	set root [file dirname $scriptDir]
}
if {![info exists partitionOutDir]} {
	set eleType "quad"
	if {[info exists soilEleType]} { set eleType $soilEleType }
	if {![info exists plotDir]} { set plotDir $scriptDir }
	set partitionOutDir [file join $plotDir out profile$soilProfile partition \
		$soilBoundary $eleType $pierEleType np$np]
}

# Format a float for JSON.
# Args: x (any numeric)
# Returns: string
proc partMapNum {x} {
	return [format "%.8g" [expr {double($x)}]]
}

# Element group from tag / coincidence (same bins as DumpModelSketch).
# Args: e ni nj (int)
# Returns: grp string
proc partMapClassify {e ni nj} {
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

# ---
# 1. OUTPUT DIR (rank 0 mkdir; then every rank writes)
# ---
if {$pid == 0} {
	file mkdir $partitionOutDir
}
barrier

set rankPath [file join $partitionOutDir [format "rank.%d.json" $pid]]
set outFd [open $rankPath w]
puts $outFd "\{"
puts $outFd [format "  \"np\": %d," $np]
puts $outFd [format "  \"pid\": %d," $pid]
puts $outFd "  \"units\": \"m\","
puts $outFd [format "  \"pierEleType\": \"%s\"," $pierEleType]
if {[info exists pileEleType]} {
	puts $outFd [format "  \"pileEleType\": \"%s\"," $pileEleType]
} else {
	puts $outFd "  \"pileEleType\": \"\","
}
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
	puts $outFd "  \"soilBoundary\": \"\","
}

puts $outFd "  \"sizes\": \{"
puts $outFd [format "    \"D_pier\": %s, \"H_pier\": %s," [partMapNum $D_pier] [partMapNum $H_pier]]
puts $outFd [format "    \"W_cap\": %s, \"H_cap\": %s, \"L_pile\": %s," \
	[partMapNum $W_cap] [partMapNum $H_cap] [partMapNum $L_pile]]
if {[info exists L_half]} {
	puts $outFd [format "    \"L_half\": %s," [partMapNum $L_half]]
}
if {[info exists w_FF]} {
	puts $outFd [format "    \"w_FF\": %s," [partMapNum $w_FF]]
}
if {[info exists xMeshHalf]} {
	puts $outFd [format "    \"xMeshHalf\": %s," [partMapNum $xMeshHalf]]
}
if {[info exists dw_deck]} {
	puts $outFd [format "    \"dw_deck\": %s," [partMapNum $dw_deck]]
} else {
	puts $outFd "    \"dw_deck\": null,"
}
if {[info exists dd_deck]} {
	puts $outFd [format "    \"dd_deck\": %s" [partMapNum $dd_deck]]
} else {
	puts $outFd "    \"dd_deck\": null"
}
puts $outFd "  \},"

# ---
# 2. LOCAL NODES + ELEMENTS (owned eles; nodes include ghosts)
# ---
set nTags [lsort -integer [getNodeTags]]
set nNode [llength $nTags]
puts $outFd "  \"nodes\": \["
set idx 0
foreach n $nTags {
	incr idx
	set xy [nodeCoord $n]
	set jsonComma [expr {($idx < $nNode) ? "," : ""}]
	puts $outFd [format "    \[%d, %s, %s\]%s" \
		$n [partMapNum [lindex $xy 0]] [partMapNum [lindex $xy 1]] $jsonComma]
}
puts $outFd "  \],"

set dumpedEles {}
foreach e [lsort -integer [getEleTags]] {
	set en [eleNodes $e]
	if {[llength $en] < 2} { continue }
	set ni [lindex $en 0]
	set nj [lindex $en 1]
	set grp [partMapClassify $e $ni $nj]
	set xyParts {}
	foreach nd $en {
		set xy [nodeCoord $nd]
		lappend xyParts [partMapNum [lindex $xy 0]] [partMapNum [lindex $xy 1]]
	}
	lappend dumpedEles [list $e $grp [join $en ", "] [join $xyParts ", "]]
}
set nEle [llength $dumpedEles]
puts $outFd "  \"elements\": \["
set idx 0
foreach row $dumpedEles {
	incr idx
	lassign $row e grp nCsv xyCsv
	set jsonComma [expr {($idx < $nEle) ? "," : ""}]
	puts $outFd [format "    \{\"e\": %d, \"grp\": \"%s\", \"n\": \[%s\], \"xy\": \[%s\]\}%s" \
		$e $grp $nCsv $xyCsv $jsonComma]
}
puts $outFd "  \]"
puts $outFd "\}"
close $outFd

barrier

# ---
# 3. MERGE + PNG (rank 0)
# ---
if {$pid == 0} {
	set python3bin [FindPython3]
	if {$python3bin eq ""} {
		puts "ExportPartitionMap: WARNING Python not found; rank JSON only"
	} else {
		if {![info exists root]} {
			set root [file dirname $scriptDir]
		}
		if {[catch {set pyOut [exec {*}$python3bin [file join $root plot PlotPartition.py] \
			$partitionOutDir]} err]} {
			puts "ExportPartitionMap: WARNING PlotPartition.py failed:\n$err"
		} elseif {$pyOut ne ""} {
			puts $pyOut
		}
	}
	puts [format "ExportPartitionMap: dir=%s" $partitionOutDir]
}
barrier
