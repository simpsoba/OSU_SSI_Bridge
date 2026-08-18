# DumpModelSketch.tcl
# Units: N, m, s
#
# Goals: write plot/model_sketch.json after the structure (+ optional soil) is built.
# Plot: python3 plot/PlotModelSketch.py
# Expects a built model and Parameters.tcl. Optional: set sketchOutPath first.

if {![info exists H_pier] || ![info exists D_pier]} {
	error "DumpModelSketch.tcl: source Parameters.tcl (and build) first"
}

set scriptDir [file dirname [file normalize [info script]]]
if {![info exists sketchOutPath]} {
	set sketchOutPath [file join $scriptDir model_sketch.json]
}

# Format a float for JSON.
# Args: x (numeric)
# Returns: string
proc sketchNum {x} {
	return [format "%.8g" [expr {double($x)}]]
}

# Element group from tag / coincidence.
# Args: e ni nj (int)
# Returns: grp string
proc sketchClassifyEle {e ni nj} {
	# Tag-based when soil/springs present
	if {[info exists eleTag_spr_base] && $e >= $eleTag_spr_base && $e < 30000} {
		return "ssi_spring"
	}
	if {$e >= 22000 && $e < 30000} {
		return "ssi_spring"
	}
	if {$e >= 35000 && $e < 40000} {
		return "soil_bnd"
	}
	if {$e >= 15000 && $e < 20000} {
		return "soil"
	}
	set ci [nodeCoord $ni]
	set cj [nodeCoord $nj]
	set dx [expr {[lindex $ci 0] - [lindex $cj 0]}]
	set dy [expr {[lindex $ci 1] - [lindex $cj 1]}]
	if {[expr {sqrt($dx*$dx + $dy*$dy)}] < 1.0e-9} {
		return "spring"
	}
	set hi [expr {($ni > $nj) ? $ni : $nj}]
	set lo [expr {($ni < $nj) ? $ni : $nj}]
	# Cap ele tags (TC may be pier base tag < 1000)
	if {[info exists ::eleTag_cap_base] && [info exists ::eleTag_pile_base]} {
		if {$e >= $::eleTag_cap_base && $e < $::eleTag_pile_base} {
			return "cap"
		}
	}
	if {$hi >= 3000 && $hi < 4000} {
		return "deck"
	}
	if {$hi >= 2000 && $hi < 3000} {
		return "pile"
	}
	if {$lo >= 1000 && $hi < 2000} {
		return "cap"
	}
	if {$ni < 1000 && $nj < 1000} {
		return "pier"
	}
	return "other"
}

set outFd [open $sketchOutPath w]
puts $outFd "\{"
puts $outFd "  \"units\": \"m\","
puts $outFd "  \"title\": \"OSU SSI Bridge\","
puts $outFd [format "  \"pierEleType\": \"%s\"," $pierEleType]
puts $outFd [format "  \"pileEleType\": \"%s\"," $pileEleType]
if {[info exists soilConstitutive]} {
	puts $outFd [format "  \"soilConstitutive\": \"%s\"," $soilConstitutive]
}
if {[info exists pileSpring]} {
	puts $outFd [format "  \"pileSpring\": \"%s\"," $pileSpring]
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
puts $outFd [format "  \"dy_pile\": %s," [sketchNum [expr {$L_pile/double($nSeg_pile)}]]]

# Fills: structure widths
puts $outFd "  \"fills\": \["
puts $outFd [format "    \{\"group\": \"pier\", \"xc\": 0, \"y0\": 0, \"y1\": %s, \"width\": %s\}," \
	[sketchNum $H_pier] [sketchNum $D_pier]]
# Deck solid outline drawn in PlotModelSketch from sizes (not a rectangle fill)
# Physical W_cap (overhang beyond the frame at +/-s)
set capFillW $W_cap
puts $outFd [format "    \{\"group\": \"cap\", \"xc\": 0, \"y0\": %s, \"y1\": 0, \"width\": %s\}," \
	[sketchNum [expr {-$H_cap}]] [sketchNum $capFillW]]
set fillI 0
foreach xP [list [expr {-$s_pile_cap}] 0.0 $s_pile_cap] {
	incr fillI
	set jsonComma [expr {($fillI < 3) ? "," : ""}]
	puts $outFd [format "    \{\"group\": \"pile\", \"xc\": %s, \"y0\": %s, \"y1\": %s, \"width\": %s\}%s" \
		[sketchNum $xP] \
		[sketchNum [expr {-$H_cap - $L_pile}]] \
		[sketchNum [expr {-$H_cap}]] \
		[sketchNum $D_pile] $jsonComma]
}
puts $outFd "  \],"

# Soil layer bands (for legend / optional strip)
puts $outFd "  \"soil_layers\": \["
if {[info exists soilLayerNames]} {
	set nLayer [llength $soilLayerNames]
	set idx 0
	foreach nm $soilLayerNames {
		incr idx
		set jsonComma [expr {($idx < $nLayer) ? "," : ""}]
		set sand [expr {[info exists soilIsSand($nm)] ? $soilIsSand($nm) : 0}]
		puts $outFd [format "    \{\"name\": \"%s\", \"y0\": %s, \"y1\": %s, \"sand\": %d\}%s" \
			$nm [sketchNum $soilYBot($nm)] [sketchNum $soilYTop($nm)] $sand $jsonComma]
	}
}
puts $outFd "  \],"

# Soil quads as 4-node polygons + layer name
puts $outFd "  \"soil_quads\": \["
set soilQuads {}
if {[info exists soilEleTags]} {
	foreach e $soilEleTags {
		set en [eleNodes $e]
		if {[llength $en] < 4} { continue }
		set nm "L2"
		if {[info exists soilEleLayer($e)]} { set nm $soilEleLayer($e) }
		set pts {}
		foreach n $en {
			set xy [nodeCoord $n]
			lappend pts [lindex $xy 0] [lindex $xy 1]
		}
		lappend soilQuads [list $e $nm $pts]
	}
}
set nSoilQuadsOut [llength $soilQuads]
set idx 0
foreach q $soilQuads {
	incr idx
	lassign $q e nm pts
	set jsonComma [expr {($idx < $nSoilQuadsOut) ? "," : ""}]
	puts $outFd [format "    \{\"e\": %d, \"layer\": \"%s\", \"xy\": \[%s, %s, %s, %s, %s, %s, %s, %s\]\}%s" \
		$e $nm \
		[sketchNum [lindex $pts 0]] [sketchNum [lindex $pts 1]] \
		[sketchNum [lindex $pts 2]] [sketchNum [lindex $pts 3]] \
		[sketchNum [lindex $pts 4]] [sketchNum [lindex $pts 5]] \
		[sketchNum [lindex $pts 6]] [sketchNum [lindex $pts 7]] \
		$jsonComma]
}
puts $outFd "  \],"

puts $outFd "  \"sizes\": \{"
puts $outFd [format "    \"D_pier\": %s, \"H_pier\": %s," [sketchNum $D_pier] [sketchNum $H_pier]]
puts $outFd [format "    \"W_cap\": %s, \"H_cap\": %s, \"s_pile_cap\": %s," \
	[sketchNum $W_cap] [sketchNum $H_cap] [sketchNum $s_pile_cap]]
if {[info exists dw_deck]} {
	puts $outFd [format "    \"dw_deck\": %s, \"dd_deck\": %s," \
		[sketchNum $dw_deck] [sketchNum $dd_deck]]
	puts $outFd [format "    \"sw_deck\": %s, \"cw_deck\": %s," \
		[sketchNum $sw_deck] [sketchNum $cw_deck]]
	puts $outFd [format "    \"td_deck\": %s, \"ts_deck\": %s, \"tw_deck\": %s," \
		[sketchNum $td_deck] [sketchNum $ts_deck] [sketchNum $tw_deck]]
}
if {[info exists W_cap_soil]} {
	puts $outFd [format "    \"W_cap_soil\": %s," [sketchNum $W_cap_soil]]
}
puts $outFd [format "    \"D_pile\": %s, \"L_pile\": %s, \"n_pile_row\": %d," \
	[sketchNum $D_pile] [sketchNum $L_pile] $n_pile_row]
if {[info exists L_half]} {
	puts $outFd [format "    \"L_half\": %s, \"t_soil\": %s," [sketchNum $L_half] [sketchNum $t_soil]]
}
if {[info exists w_FF]} {
	puts $outFd [format "    \"w_FF\": %s," [sketchNum $w_FF]]
}
if {[info exists xMeshHalf]} {
	puts $outFd [format "    \"xMeshHalf\": %s," [sketchNum $xMeshHalf]]
}
puts $outFd [format "    \"foot\": %s" [sketchNum $foot]]
puts $outFd "  \},"

# Structure + spring nodes only (skip continuum soil nodes for clarity)
puts $outFd "  \"nodes\": \["
set nTags {}
foreach n [lsort -integer [getNodeTags]] {
	if {$n >= 10000 && $n < 20000} { continue }
	lappend nTags $n
}
set nNode [llength $nTags]
set idx 0
foreach n $nTags {
	incr idx
	set xy [nodeCoord $n]
	set jsonComma [expr {($idx < $nNode) ? "," : ""}]
	puts $outFd [format "    \[%d, %s, %s\]%s" \
		$n [sketchNum [lindex $xy 0]] [sketchNum [lindex $xy 1]] $jsonComma]
}
puts $outFd "  \],"

# Beam / spring elements (not soil quads)
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
	set grp [sketchClassifyEle $e $ni $nj]
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

# Rotational ZLS (lumpedPlasticity: base ele 1 and top ele 3)
puts $outFd "  \"springs\": \["
set sprDump {}
foreach row $beamEles {
	lassign $row e ni nj grp
	if {$grp eq "spring"} {
		set xy [nodeCoord $ni]
		lappend sprDump [list $e $ni $nj [lindex $xy 0] [lindex $xy 1]]
	}
}
set nSpring [llength $sprDump]
set idx 0
foreach s $sprDump {
	incr idx
	lassign $s e ni nj x y
	set jsonComma [expr {($idx < $nSpring) ? "," : ""}]
	puts $outFd [format "    \[%d, %d, %d, %s, %s\]%s" \
		$e $ni $nj [sketchNum $x] [sketchNum $y] $jsonComma]
}
puts $outFd "  \],"

# SSI springs: [e, xs,ys, xi,yi, xp,yp, kind, type]
#   type = py | tz | qz | pyliq | tzliq  (one row per uniaxial mat)
puts $outFd "  \"ssi_springs\": \["
set ssiSprings {}
if {[info exists ssiSpringDump] && [llength $ssiSpringDump] > 0} {
	set ssiSprings $ssiSpringDump
} else {
	foreach row $beamEles {
		lassign $row e ni nj grp
		if {$grp eq "ssi_spring"} {
			set xy [nodeCoord $ni]
			set x [lindex $xy 0]
			set y [lindex $xy 1]
			lappend ssiSprings [list $e $x $y $x $y $x $y "pile" "py"]
		}
	}
}
set nSsi [llength $ssiSprings]
set idx 0
foreach s $ssiSprings {
	incr idx
	lassign $s e xs ys xi yi xp yp kind stype
	set jsonComma [expr {($idx < $nSsi) ? "," : ""}]
	puts $outFd [format "    \[%d, %s, %s, %s, %s, %s, %s, \"%s\", \"%s\"]%s" \
		$e [sketchNum $xs] [sketchNum $ys] \
		[sketchNum $xi] [sketchNum $yi] \
		[sketchNum $xp] [sketchNum $yp] $kind $stype $jsonComma]
}
puts $outFd "  \],"

# ASDEA absorbing quads (4-node) and/or Shin Lysmer dashpots
puts $outFd "  \"bnd_quads\": \["
set bndQuads {}
if {[info exists soilEleBndTags]} {
	foreach e $soilEleBndTags {
		set en [eleNodes $e]
		if {[llength $en] < 4} { continue }
		set pts {}
		foreach n $en {
			set xy [nodeCoord $n]
			lappend pts [lindex $xy 0] [lindex $xy 1]
		}
		lappend bndQuads [list $e $pts]
	}
}
set nBndQuad [llength $bndQuads]
set idx 0
foreach q $bndQuads {
	incr idx
	lassign $q e pts
	set jsonComma [expr {($idx < $nBndQuad) ? "," : ""}]
	puts $outFd [format "    \{\"e\": %d, \"xy\": \[%s, %s, %s, %s, %s, %s, %s, %s\]\}%s" \
		$e \
		[sketchNum [lindex $pts 0]] [sketchNum [lindex $pts 1]] \
		[sketchNum [lindex $pts 2]] [sketchNum [lindex $pts 3]] \
		[sketchNum [lindex $pts 4]] [sketchNum [lindex $pts 5]] \
		[sketchNum [lindex $pts 6]] [sketchNum [lindex $pts 7]] \
		$jsonComma]
}
puts $outFd "  \],"

puts $outFd "  \"lysmer_dashpots\": \["
set loadDump {}
if {[info exists lysmerDashpotDump]} {
	set loadDump $lysmerDashpotDump
}
set nLoad [llength $loadDump]
set idx 0
foreach row $loadDump {
	incr idx
	lassign $row e x y role c
	set jsonComma [expr {($idx < $nLoad) ? "," : ""}]
	puts $outFd [format "    \{\"e\": %d, \"x\": %s, \"y\": %s, \"role\": \"%s\", \"c\": %s\}%s" \
		$e [sketchNum $x] [sketchNum $y] $role [sketchNum $c] $jsonComma]
}
puts $outFd "  \],"

puts $outFd "  \"pile_tips\": \["
if {[info exists nodeTag_pile_tips]} {
	set nTip [llength $nodeTag_pile_tips]
	set idx 0
	foreach t $nodeTag_pile_tips {
		incr idx
		set jsonComma [expr {($idx < $nTip) ? "," : ""}]
		puts $outFd [format "    %d%s" $t $jsonComma]
	}
}
puts $outFd "  \]"

puts $outFd "\}"
close $outFd

