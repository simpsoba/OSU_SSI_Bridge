# BuildPilesElements.tcl
# Units: N, m, s
#
# Goals: shaft beams under the cap. Call after BuildPilesNodes.
# tag = eleTag_pile_base + ip*nSeg + (iy-1); nodes: pileNodeTag ip iy

if {[info exists pileElementsDone] && $pileElementsDone} {
	return
}
if {![info exists pileHeads] || ![info exists nSeg_pile]} {
	error "BuildPilesElements.tcl: source BuildPilesNodes.tcl first"
}

# =====================================================================
# 4. ELEMENTS
# =====================================================================
# tag = eleTag_pile_base + ip*nSeg + (iy-1); nodes: pileNodeTag ip iy
geomTransf $pileGeoTransf $transfTag_pile

set e [expr {$eleTag_pile_base - 1}]

for {set ip 0} {$ip < $n_pile} {incr ip} {
	lassign [lindex $pileHeads $ip] headTag xP
	set prev $headTag
	for {set iy 1} {$iy <= $nSeg_pile} {incr iy} {
		set nTag [pileNodeTag $ip $iy]
		incr e
		if {$pileEleType eq "elasticBeamColumn"} {
			# element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
			element elasticBeamColumn $e \
				$prev $nTag $A_pile $Es_pile $I_pile $transfTag_pile
		} elseif {$pileEleType eq "dispBeamColumn"} {
			# element dispBeamColumn $eleTag $iNode $jNode $numIntgrPts $secTag $transfTag
			element dispBeamColumn $e \
				$prev $nTag $nIP_pile $secTag_pile $transfTag_pile
		} else {
			error "BuildPilesElements.tcl: unknown pileEleType '$pileEleType'"
		}
		set prev $nTag
	}
}

set eleTag_pile_last $e
set pileElementsDone 1
