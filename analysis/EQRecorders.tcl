# analysis/EQRecorders.tcl
# Goals: node and element recorders for the EQ run, sourced after
# `analysis Transient`. Units N, m, s. Every recorder samples at the ground
# motion step (-dT $gmVelDT), not at every dtAnalysis step.
#
# recordersON
#   0  off
#   1  full window: nodes with |x| <= eqWindowX, the quads inside it, every
#      pile beam, every SSI spring (pile, cap face, cap soffit)
#   2  center column: pier nodes (UX UY RZ), both rotational springs, the
#      soil-base primary node, the whole center pile, every center-pile
#      spring, and every x=0 soil quad (grade to base). No cap springs.
#      Quad corners are geometry only (no soil UX).
#   3  nine SSI horizons (old 2): first / mid / last station of L2, L3, L5
#      on the center pile, each with its spring, pile segment, and x=0 quad.
#
# Serial (getNP = 1) writes $name; OpenSeesMP writes $name.$pid after rank 0
# clears eqOutDir. plot/PlotEQ.py reads the serial names, plot/PlotEQParallel.py
# stitches the shards and then calls PlotEQ.py.
# Drivers call `record` once after this file (and after OpenFresco recorders).
#
# 1. RANK, OUTPUT FOLDER, HELPERS
# 2. MONITOR SET      -- which nodes, quads, pile beams and springs
# 3. INDEX FILES      -- one line per recorder column block, in column order
# 4. RECORDERS
# 5. window_meta.txt

# ---
# 1. RANK, OUTPUT FOLDER, HELPERS
# ---
set eqNP 1
set eqPID 0
if {[llength [info commands getNP]]} {
	set eqNP [getNP]
	set eqPID [getPID]
}

if {![info exists recordersON]} {
	set recordersON 1
}
if {![string is integer -strict $recordersON] || $recordersON < 0 || $recordersON > 3} {
	error "EQRecorders.tcl: recordersON must be an integer 0..3 (got '$recordersON')"
}
set eqRecLean [expr {$recordersON == 2 || $recordersON == 3}]
set eqRecNine [expr {$recordersON == 3}]
if {$eqNP > 1} {
	set eqRecSuf [format ".%d" $eqPID]
} else {
	set eqRecSuf ""
}

# Output path for one recorder file.
# Args: name (e.g. window_disp.out)
# Returns: absolute path ($name serial, $name.$pid under OpenSeesMP)
proc eqRecPath {name} {
	global eqOutDir eqRecSuf
	return [file join $eqOutDir ${name}${eqRecSuf}]
}

# Does this rank hold the node / element? Serial: every tag in the model.
# Args: tag
# Returns: 1 | 0
proc eqOwnsNode {n} {
	global eqLocalNode
	return [info exists eqLocalNode($n)]
}
proc eqOwnsEle {e} {
	global eqLocalEle
	return [info exists eqLocalEle($e)]
}

# Element tags this rank holds, same order as the input.
# Args: tags (list; a single tag is a one-element list)
# Returns: subset
proc eqOwnedEles {tags} {
	set out {}
	foreach e $tags {
		if {[eqOwnsEle $e]} { lappend out $e }
	}
	return $out
}

# Element tags this rank holds in [$e0, $e1] inclusive.
# Args: e0 e1
# Returns: ascending subset
proc eqOwnedEleRange {e0 e1} {
	set out {}
	for {set e $e0} {$e <= $e1} {incr e} {
		if {[eqOwnsEle $e]} { lappend out $e }
	}
	return $out
}

# Sample at the ground-motion step, so file length follows the record and not
# dtAnalysis. 0 = every analysis step (no -dT).
set eqRecDt 0.0
if {[info exists gmVelDT] && $gmVelDT > 0.0} {
	set eqRecDt $gmVelDT
}

# Time flags shared by every recorder below.
# Args: none
# Returns: {-time -dT $eqRecDt}, or {-time} when eqRecDt is 0
proc eqRecTime {} {
	global eqRecDt
	if {$eqRecDt > 0.0} {
		return [list -time -dT $eqRecDt]
	}
	return [list -time]
}

# recorder Node -file $fileName -time -dT $dT -node $nodeTags -dof $dofs $rsp
# Args: name (file), tags (node list), dofs (list), rsp (disp | accel)
# Returns: none (writes nothing when tags is empty, e.g. no local nodes)
proc eqRecNode {name tags dofs rsp} {
	if {[llength $tags] < 1} { return }
	recorder Node -file [eqRecPath $name] {*}[eqRecTime] \
		-node {*}$tags -dof {*}$dofs $rsp
}

# recorder Element -file $fileName -time -dT $dT -ele $eleTags $args
# Args: name (file), tags (element list), args (e.g. section 1 deformation)
# Returns: none (writes nothing when tags is empty)
proc eqRecEle {name tags args} {
	if {[llength $tags] < 1} { return }
	recorder Element -file [eqRecPath $name] {*}[eqRecTime] \
		-ele {*}$tags {*}$args
}

# Recorders live under plot/out next to the figures. Serial and parallel keep
# separate folders so the two dumps cannot overwrite each other.
# outDIR in Run.tcl / RunParallel.tcl overrides; "" keeps this auto path.
if {[info exists outDIR] && $outDIR ne ""} {
	set eqOutDir $outDIR
}
if {![info exists eqOutDir]} {
	if {$eqNP > 1} {
		set eqRunKind "parallel"
	} else {
		set eqRunKind "serial"
	}
	if {![info exists plotDir]} {
		set plotDir [file join $root plot]
	}
	set eqOutDir [file join $plotDir out profile$soilProfile eq \
		$eqRunKind $soilBoundary $soilEleType $pierEleType]
}

if {$recordersON == 0} {
	if {$eqNP <= 1 || $eqPID == 0} {
		puts "EQRecorders: recordersON=0 -- no recorders"
	}
	return
}

# Clear the folder first: leftover name.$pid from another -np would be stitched.
if {$eqNP <= 1 || $eqPID == 0} {
	if {[file isdirectory $eqOutDir]} {
		foreach oldOut [glob -nocomplain [file join $eqOutDir *]] {
			if {[file isdirectory $oldOut]} { continue }
			file delete -force $oldOut
		}
	}
	file mkdir $eqOutDir
}
if {$eqNP > 1} {
	barrier
}

array unset eqLocalEle
array unset eqLocalNode
foreach e [getEleTags] { set eqLocalEle($e) 1 }
foreach n [getNodeTags] { set eqLocalNode($n) 1 }

# ---
# 2. MONITOR SET
# ---
# Center-pile SSI stations for lean dumps. iy comes from the spring dump
# (dy_soil and the layer thicknesses), not a hardcoded tag list.
# Args: nine (1 = first/mid/last of L2, L3, L5; 0 = every center-pile station)
# Returns: none (sets eqLeanStations {iy y layer isTip iSeg}, and the lookups
#          eqLeanStationKey(ip,iy) and eqLeanSegKeep(iSeg))
proc eqLeanBuildStations {nine} {
	global pileSpringPropsDump n_pile nSeg_pile
	global eqLeanStations eqLeanStationKey eqLeanSegKeep
	array unset eqLeanStationKey
	array unset eqLeanSegKeep
	set eqLeanStations {}
	if {![info exists pileSpringPropsDump]} { return }
	set ipC [expr {($n_pile - 1)/2}]
	set rows {}
	if {$nine} {
		array unset byLayer
		foreach rec $pileSpringPropsDump {
			lassign $rec ip iy y depth nm isTip
			if {$ip != $ipC} { continue }
			lappend byLayer($nm) [list $iy $y $nm $isTip]
		}
		foreach nm {L2 L3 L5} {
			if {![info exists byLayer($nm)]} { continue }
			set layerRows [lsort -integer -index 0 $byLayer($nm)]
			set nRow [llength $layerRows]
			set pick [list 0 [expr {($nRow - 1)/2}] [expr {$nRow - 1}]]
			foreach i [lsort -integer -unique $pick] {
				lappend rows [lindex $layerRows $i]
			}
		}
	} else {
		foreach rec $pileSpringPropsDump {
			lassign $rec ip iy y depth nm isTip
			if {$ip != $ipC} { continue }
			lappend rows [list $iy $y $nm $isTip]
		}
	}
	foreach rec [lsort -integer -index 0 $rows] {
		lassign $rec iy y nm isTip
		# station iy is the i-end of segment iy+1; the tip is the j-end of
		# the last segment, so it reuses iSeg = nSeg_pile
		set iSeg [expr {$iy + 1}]
		if {$iSeg > $nSeg_pile} { set iSeg $nSeg_pile }
		if {[info exists eqLeanStationKey($ipC,$iy)]} { continue }
		set eqLeanStationKey($ipC,$iy) 1
		set eqLeanSegKeep($iSeg) 1
		lappend eqLeanStations [list $iy $y $nm $isTip $iSeg]
	}
}

# Soil row whose top or bottom face sits at the pile station y and whose
# material is still $layer. The two candidate rows differ at an L3/L5 contact;
# at the tip the row above is the one the pile still loads.
# Args: yPile (m), layer (L2 | L3 | L5), isTip (1 | 0)
# Returns: row index iy, or -1
proc eqQuadRowForPile {yPile layer isTip} {
	global soilYs nSoilRows soilRowLayer
	set yTol 1.0e-4
	set iyTop -1
	set iyBot -1
	for {set iy 0} {$iy < $nSoilRows} {incr iy} {
		if {abs([lindex $soilYs $iy] - $yPile) < $yTol} { set iyTop $iy }
		if {abs([lindex $soilYs [expr {$iy + 1}]] - $yPile) < $yTol} { set iyBot $iy }
	}
	if {$isTip && $iyBot >= 0 && $soilRowLayer($iyBot) eq $layer} { return $iyBot }
	if {$iyTop >= 0 && $soilRowLayer($iyTop) eq $layer} { return $iyTop }
	if {$iyBot >= 0 && $soilRowLayer($iyBot) eq $layer} { return $iyBot }
	return -1
}

if {$eqRecLean} {
	eqLeanBuildStations $eqRecNine
}

# Pier: nodes 2 and 4 are the zeroLength inner nodes, so they exist only for
# pierEleType lumpedPlasticity. Cap TC aliases node 1, deck BC aliases node 5.
set eqPierNodes {}
foreach nameVar {nodeTag_pierBase_capTC nodeTag_pierBaseZeroLengthInner \
		nodeTag_pierTopZeroLengthInner nodeTag_pierTop_deckBC} {
	if {![info exists $nameVar]} { continue }
	set n [set $nameVar]
	if {[eqOwnsNode $n]} { lappend eqPierNodes $n }
}

# Pile beams: eleTag_pile_base + ip*nSeg_pile + (iSeg - 1), shaft by shaft.
# recordersON=2: every center-pile segment. 3: the nine station segments.
set eqPileBeamRows {}
if {[info exists eleTag_pile_base] && [info exists eleTag_pile_last] \
		&& $eleTag_pile_last >= $eleTag_pile_base} {
	set ipCenter [expr {($n_pile - 1)/2}]
	for {set ip 0} {$ip < $n_pile} {incr ip} {
		for {set iSeg 1} {$iSeg <= $nSeg_pile} {incr iSeg} {
			if {$eqRecLean && $ip != $ipCenter} { continue }
			if {$eqRecNine && ![info exists eqLeanSegKeep($iSeg)]} {
				continue
			}
			set e [expr {$eleTag_pile_base + $ip*$nSeg_pile + $iSeg - 1}]
			if {![eqOwnsEle $e]} { continue }
			lappend eqPileBeamRows [list $e $ip $iSeg]
		}
	}
	set ePileEnd [expr {$eleTag_pile_base + $n_pile*$nSeg_pile - 1}]
	if {$ePileEnd != $eleTag_pile_last} {
		puts [format "EQRecorders: WARNING pile beams end at %d, eleTag_pile_last=%d" \
			$ePileEnd $eleTag_pile_last]
	}
}

# Pile SSI springs: one zeroLength per station (dir 1 = p-y, dir 2 = t-z, or q-z
# at the tip), built head to tip, so pileSprEleRec order is the column order.
set eqPileSprRows {}
if {[info exists pileSprEleRec]} {
	foreach rec $pileSprEleRec {
		lassign $rec e ip iy isTip
		if {$eqRecLean && $ip != [expr {($n_pile - 1)/2}]} { continue }
		if {$eqRecNine && ![info exists eqLeanStationKey($ip,$iy)]} { continue }
		if {![eqOwnsEle $e]} { continue }
		lappend eqPileSprRows [list $e $ip $iy]
	}
}
# Cap tags follow the pile springs: the six face springs, then the soffit q-z row.
set eqCapFaceEles {}
set eqCapSoffitEles {}
if {!$eqRecLean && [info exists nPileSprings] && [info exists eleTag_spr_last]} {
	set nSof 0
	if {[info exists nSoffit]} { set nSof $nSoffit }
	set eCapFirst [expr {$eleTag_spr_base + $nPileSprings}]
	set eCapLast [expr {$eleTag_spr_last - $nSof}]
	set eqCapFaceEles [eqOwnedEleRange $eCapFirst $eCapLast]
	set eqCapSoffitEles [eqOwnedEleRange [expr {$eCapLast + 1}] $eleTag_spr_last]
}

# Quads, geometry nodes (window_nodes.txt) and displacement nodes (window_disp).
#   1: nodes with |x| <= eqWindowX, then the quads with all four nodes there.
#      Geometry and displacement are the same set.
#   2: every x=0 soil quad (grade to base). Geometry carries the quad corners;
#      only the pier and the center-pile nodes get a displacement channel.
#   3: one x=0 quad per nine-horizon station (same geometry/disp split as 2).
#      Quad peak and hysteresis plots read stress and strain, so the corners
#      need coordinates only. Soil UX at a station follows from
#      u_pile - (spring deformation dir 1); the vertical needs the dup/pile
#      gravity offset too (see NOTES.md).
# soilEleTags index = ix*nSoilRows + iy (BuildSoilMesh.tcl loops ix, then iy).
set eqQuadEles {}
set eqGeomNodeTags {}
set eqDispNodeTags {}
array unset eqInWindow
if {$eqRecLean} {
	set ixCenter -1
	for {set ix 0} {$ix < [llength $soilXs] - 1} {incr ix} {
		if {abs([lindex $soilXs $ix]) < 1.0e-9} {
			set ixCenter $ix
			break
		}
	}
	if {$eqRecNine} {
		array unset rowSeen
		foreach st $eqLeanStations {
			lassign $st iy y nm isTip iSeg
			set iyQ [eqQuadRowForPile $y $nm $isTip]
			if {$ixCenter < 0 || $iyQ < 0 || [info exists rowSeen($iyQ)]} {
				continue
			}
			set rowSeen($iyQ) 1
			set e [lindex $soilEleTags [expr {$ixCenter*$nSoilRows + $iyQ}]]
			if {![eqOwnsEle $e]} { continue }
			lappend eqQuadEles $e
		}
	} else {
		for {set iyQ 0} {$iyQ < $nSoilRows} {incr iyQ} {
			if {$ixCenter < 0} { break }
			set e [lindex $soilEleTags [expr {$ixCenter*$nSoilRows + $iyQ}]]
			if {![eqOwnsEle $e]} { continue }
			lappend eqQuadEles $e
		}
	}
	set dispNodes $eqPierNodes
	foreach row $eqPileBeamRows {
		lappend dispNodes {*}[eleNodes [lindex $row 0]]
	}
	set eqDispNodeTags [lsort -integer -unique $dispNodes]
	set geomNodes $eqDispNodeTags
	foreach e $eqQuadEles {
		lappend geomNodes {*}[eleNodes $e]
	}
	set eqGeomNodeTags [lsort -integer -unique $geomNodes]
} else {
	foreach n [lsort -integer [getNodeTags]] {
		if {abs([lindex [nodeCoord $n] 0]) > $eqWindowX + 1.0e-9} { continue }
		lappend eqGeomNodeTags $n
	}
	set eqDispNodeTags $eqGeomNodeTags
	foreach n $eqGeomNodeTags { set eqInWindow($n) 1 }
	foreach e $soilEleTags {
		if {![eqOwnsEle $e]} { continue }
		if {[catch {set enodes [eleNodes $e]}]} { continue }
		if {[llength $enodes] < 4} { continue }
		set inside 1
		foreach en $enodes {
			if {![info exists eqInWindow($en)]} { set inside 0; break }
		}
		if {$inside} { lappend eqQuadEles $e }
	}
}
array unset eqInWindow
foreach n $eqGeomNodeTags { set eqInWindow($n) 1 }

# ---
# 3. INDEX FILES
# ---
# Geometry: every node a plot script needs coordinates for.
set nodesFd [open [eqRecPath window_nodes.txt] w]
if {$eqRecNine} {
	puts $nodesFd "# tag x y  (pier, nine SSI stations, station quad corners)"
} elseif {$eqRecLean} {
	puts $nodesFd "# tag x y  (pier, center pile, x=0 soil column)"
} else {
	puts $nodesFd "# tag x y  (|x| <= $eqWindowX m)"
}
foreach n $eqGeomNodeTags {
	set xy [nodeCoord $n]
	puts $nodesFd [format "%d %.8g %.8g" $n [lindex $xy 0] [lindex $xy 1]]
}
close $nodesFd

# One line per (ux uy) column pair in window_disp.out -- a subset of the above.
set dispFd [open [eqRecPath disp_nodes.txt] w]
puts $dispFd "# tag x y  (column order in window_disp*.out)"
foreach n $eqDispNodeTags {
	set xy [nodeCoord $n]
	puts $dispFd [format "%d %.8g %.8g" $n [lindex $xy 0] [lindex $xy 1]]
}
close $dispFd

# Any local element with all nodes in the list above: pier and pile segments
# (2 nodes) draw as lines, quads (4 nodes) as patches.
set elesFd [open [eqRecPath window_eles.txt] w]
puts $elesFd {# eleTag n1 n2 [n3 n4]  (all nodes in window_nodes.txt)}
set nWinEle 0
foreach e [getEleTags] {
	if {[catch {set enodes [eleNodes $e]}]} { continue }
	set inside 1
	foreach en $enodes {
		if {![info exists eqInWindow($en)]} { set inside 0; break }
	}
	if {!$inside} { continue }
	puts $elesFd [concat $e $enodes]
	incr nWinEle
}
close $elesFd

set quadsFd [open [eqRecPath window_quads.txt] w]
if {$eqRecNine} {
	puts $quadsFd {# eleTag layer  (x=0 quad at each SSI station)}
} elseif {$eqRecLean} {
	puts $quadsFd {# eleTag layer  (x=0 column, grade to base)}
} else {
	puts $quadsFd {# eleTag layer  (all four nodes inside the window)}
}
foreach e $eqQuadEles {
	set nm "?"
	if {[info exists soilEleLayer($e)]} { set nm $soilEleLayer($e) }
	puts $quadsFd [format "%d %s" $e $nm]
}
close $quadsFd

set pileBeamFd [open [eqRecPath pile_beam_eles.txt] w]
puts $pileBeamFd "# eleTag ip iSeg   (section 1 = first IP = i-end)"
foreach row $eqPileBeamRows {
	puts $pileBeamFd [format "%d %d %d" {*}$row]
}
close $pileBeamFd

set pileSprFd [open [eqRecPath pile_springs_eles.txt] w]
puts $pileSprFd "# eleTag kind ip iy   (kind: pile)"
foreach row $eqPileSprRows {
	puts $pileSprFd [format "%d pile %d %d" {*}$row]
}
close $pileSprFd

if {[llength $eqCapFaceEles] > 0 || [llength $eqCapSoffitEles] > 0} {
	set capSprFd [open [eqRecPath cap_springs_eles.txt] w]
	puts $capSprFd "# eleTag kind   (cap = face p-y/t-z; cap_soffit = q-z dir 2)"
	foreach e $eqCapFaceEles {
		puts $capSprFd [format "%d cap" $e]
	}
	foreach e $eqCapSoffitEles {
		puts $capSprFd [format "%d cap_soffit" $e]
	}
	close $capSprFd
}

# ---
# 4. RECORDERS
# ---
# UX, UY: two columns per node in disp_nodes.txt order. Chunks keep one file
# from growing past a few hundred columns.
set nodeChunk 250
set nGeomNode [llength $eqGeomNodeTags]
set nDispNode [llength $eqDispNodeTags]
set windowDispFiles {}
for {set i0 0} {$i0 < $nDispNode} {incr i0 $nodeChunk} {
	if {$nDispNode <= $nodeChunk} {
		set name window_disp.out
	} else {
		set name [format "window_disp_%02d.out" [expr {$i0/$nodeChunk}]]
	}
	eqRecNode $name [lrange $eqDispNodeTags $i0 [expr {$i0 + $nodeChunk - 1}]] \
		{1 2} disp
	lappend windowDispFiles [file tail [eqRecPath $name]]
}

# Quad (4 Gauss pts) or SSPquad (1 IP at the centroid).
# PIMY 2D getStress is (sxx, syy, sxy) so sxy = tau_xy; getStrain is
# (exx, eyy, gxy) so gxy = gamma_xy. PDMY02 packs 2D the same way.
# material $ip stress would add sigma_zz and eta_r; keep the 3-comp element query.
if {$soilEleType eq "SSPquad"} {
	set sigRsp "stress2D3"
	set epsRsp "strain2D3"
	set nGPq 1
	set quadChunk 120
} else {
	set sigRsp "stresses"
	set epsRsp "strains"
	set nGPq 4
	set quadChunk 40
}
set nQuadRec [llength $eqQuadEles]
set quadStressFiles {}
set quadStrainFiles {}
for {set i0 0} {$i0 < $nQuadRec} {incr i0 $quadChunk} {
	if {$nQuadRec <= $quadChunk} {
		set sigName window_quad_stress.out
		set epsName window_quad_strain.out
	} else {
		set iC [expr {$i0/$quadChunk}]
		set sigName [format "window_quad_stress_%02d.out" $iC]
		set epsName [format "window_quad_strain_%02d.out" $iC]
	}
	set batch [lrange $eqQuadEles $i0 [expr {$i0 + $quadChunk - 1}]]
	eqRecEle $sigName $batch $sigRsp
	eqRecEle $epsName $batch $epsRsp
	lappend quadStressFiles [file tail [eqRecPath $sigName]]
	lappend quadStrainFiles [file tail [eqRecPath $epsName]]
}

# Pier: one file per node (UX UY RZ) so a rank that owns part of the pier still
# writes labelled columns. pier_top_disp.out stays for the Shin/ASDEA overlay.
set pierNodeFiles {}
foreach n $eqPierNodes {
	set name [format "pier_node_%d.out" $n]
	eqRecNode $name $n {1 2 3} disp
	lappend pierNodeFiles [file tail [eqRecPath $name]]
}
if {[eqOwnsNode $nodeTag_pierTop_deckBC]} {
	eqRecNode pier_top_disp.out $nodeTag_pierTop_deckBC {1 2 3} disp
	if {!$eqRecLean} {
		eqRecNode pier_top_acc.out $nodeTag_pierTop_deckBC {1 2 3} accel
	}
}

# Shin near-field base: all other base nodes equalDOF to nPrimary in UX, and the
# Lysmer dashpot plus the 2cv load sit on it, so this is the input motion check.
if {[info exists nPrimary] && [eqOwnsNode $nPrimary]} {
	eqRecNode soil_base_primary.out $nPrimary {1 2} disp
}

# Pier hinges: the zeroLength sections (lumpedPlasticity) or the i-end IP of the
# forceBeamColumn. 2D fiber section: force = (P, Mz), deformation = (eps, kappa).
# For a ZLS those two are axial displacement (m) and rotation theta (rad).
if {$pierEleType eq "lumpedPlasticity"} {
	set botSpr [eqOwnedEles $eleTag_pier_botSpr]
	set topSpr [eqOwnedEles $eleTag_pier_topSpr]
	eqRecEle pier_hinge_force.out $botSpr section force
	eqRecEle pier_hinge_defo.out $botSpr section deformation
	eqRecEle pier_hinge_top_force.out $topSpr section force
	eqRecEle pier_hinge_top_defo.out $topSpr section deformation
} elseif {$pierEleType eq "forceBeamColumn"} {
	set pierEle [eqOwnedEles $eleTag_pier]
	eqRecEle pier_hinge_force.out $pierEle section 1 force
	eqRecEle pier_hinge_defo.out $pierEle section 1 deformation
}

# Pile shafts: global end forces (2D: Px Py Mz at i and j).
# dispBeamColumn section 1 = first IP (Lobatto i-end: eps, kappa).
set pileBeamEles {}
foreach row $eqPileBeamRows {
	lappend pileBeamEles [lindex $row 0]
}
eqRecEle pile_beam_globalForce.out $pileBeamEles globalForce
if {$pileEleType eq "dispBeamColumn"} {
	eqRecEle pile_beam_sec1_defo.out $pileBeamEles section 1 deformation
}

# Interface zeroLength: localForce and deformation, dir 1 = p-y, dir 2 = t-z.
# The soffit q-z row is dir 2 only, so it gets its own files and the column
# count of each file stays even.
set pileSprEles {}
foreach row $eqPileSprRows {
	lappend pileSprEles [lindex $row 0]
}
eqRecEle pile_springs_force.out $pileSprEles localForce
eqRecEle pile_springs_defo.out $pileSprEles deformation
eqRecEle cap_springs_force.out $eqCapFaceEles localForce
eqRecEle cap_springs_defo.out $eqCapFaceEles deformation
eqRecEle cap_springs_soffit_force.out $eqCapSoffitEles localForce
eqRecEle cap_springs_soffit_defo.out $eqCapSoffitEles deformation

# ---
# 5. window_meta.txt
# ---
set metaFd [open [eqRecPath window_meta.txt] w]
puts $metaFd "soilProfile $soilProfile"
puts $metaFd "soilBoundary $soilBoundary"
puts $metaFd "soilEleType $soilEleType"
puts $metaFd "pierHinge $pierEleType"
puts $metaFd "pileEleType $pileEleType"
puts $metaFd "recordersON $recordersON"
puts $metaFd "eqWindowX $eqWindowX"
puts $metaFd "nWindowNodes $nGeomNode"
puts $metaFd "nDispNodes $nDispNode"
puts $metaFd "nWindowEles $nWinEle"
puts $metaFd "nWindowQuads $nQuadRec"
puts $metaFd "dtAnalysis $dtAnalysis"
puts $metaFd "recDt $eqRecDt"
puts $metaFd "eqNsteps $eqNsteps"
set tMeta $Trec
if {[info exists gmStartTime] && $gmStartTime ne "" && $gmStartTime > 0} {
	puts $metaFd "gmStartTime $gmStartTime"
	set tMeta [expr {$tMeta + $gmStartTime}]
}
if {[info exists fvNsteps]} {
	puts $metaFd "freeVibT $eqFreeVibT"
	puts $metaFd "freeVibNsteps $fvNsteps"
	puts $metaFd "eqNstepsAll $eqNstepsAll"
	puts $metaFd "Trec [expr {$tMeta + $eqFreeVibT}]"
} else {
	puts $metaFd "Trec $tMeta"
}
puts $metaFd "gmVelFile $gmVelFile"
puts $metaFd "gmScaleFactor $gmScaleFactor"
puts $metaFd "tagShift_soil $tagShift_soil"
puts $metaFd "nodeTag_sprSoil_base $nodeTag_sprSoil_base"
puts $metaFd "nodeTag_bnd_base $nodeTag_bnd_base"
puts $metaFd "sprSoffitOff $sprSoffitOff"
if {[info exists soilNodeLast]} {
	puts $metaFd "soilNodeLast $soilNodeLast"
}
if {[info exists soilNodeStride]} {
	puts $metaFd "soilNodeStride $soilNodeStride"
}
if {$eqNP > 1} {
	puts $metaFd "pid $eqPID"
	puts $metaFd "np $eqNP"
	puts $metaFd "fileSuffix $eqRecSuf"
}
puts $metaFd "dispNodesFile disp_nodes.txt"
puts $metaFd "dispFiles $windowDispFiles"
puts $metaFd "dispFormat time then (ux uy) per node, disp_nodes.txt order"
puts $metaFd "pierNodeFiles $pierNodeFiles"
puts $metaFd "pierNodeFormat time then (ux uy rz) for the node in the file name"
if {$eqRecLean} {
	puts $metaFd "leanStations $eqLeanStations"
	puts $metaFd "leanStationFormat {iy y layer isTip iSeg} per SSI horizon"
}
if {$nQuadRec > 0} {
	puts $metaFd "quadEleFile window_quads.txt"
	puts $metaFd "quadNgp $nGPq"
	puts $metaFd "quadStressRsp $sigRsp"
	puts $metaFd "quadStrainRsp $epsRsp"
	puts $metaFd "quadStressFiles $quadStressFiles"
	puts $metaFd "quadStrainFiles $quadStrainFiles"
	puts $metaFd "quadStressFormat time then nGP*(sxx syy sxy) per ele; sxy=tau_xy"
	puts $metaFd "quadStrainFormat time then nGP*(exx eyy gxy) per ele; gxy=gamma_xy"
}
if {$pierEleType eq "lumpedPlasticity" || $pierEleType eq "forceBeamColumn"} {
	puts $metaFd "hingeForceFile pier_hinge_force.out"
	puts $metaFd "hingeDefoFile pier_hinge_defo.out"
}
if {$pierEleType eq "lumpedPlasticity"} {
	puts $metaFd "hingeTopForceFile pier_hinge_top_force.out"
	puts $metaFd "hingeTopDefoFile pier_hinge_top_defo.out"
}
if {[llength $eqPileBeamRows] > 0} {
	puts $metaFd "pileBeamGlobalForceFile pile_beam_globalForce.out"
	if {$pileEleType eq "dispBeamColumn"} {
		puts $metaFd "pileBeamSec1DefoFile pile_beam_sec1_defo.out"
		if {[info exists nIP_pile]} {
			puts $metaFd "nIP_pile $nIP_pile"
		}
	}
}
if {[llength $eqPileSprRows] > 0} {
	puts $metaFd "pileSpringsForceFile pile_springs_force.out"
	puts $metaFd "pileSpringsDefoFile pile_springs_defo.out"
}
if {[llength $eqCapFaceEles] > 0} {
	puts $metaFd "capSpringsForceFile cap_springs_force.out"
	puts $metaFd "capSpringsDefoFile cap_springs_defo.out"
}
if {[llength $eqCapSoffitEles] > 0} {
	puts $metaFd "capSoffitForceFile cap_springs_soffit_force.out"
	puts $metaFd "capSoffitDefoFile cap_springs_soffit_defo.out"
}
if {[info exists nPrimary]} {
	puts $metaFd "soilBasePrimaryFile soil_base_primary.out"
	puts $metaFd "soilBasePrimaryNode $nPrimary"
}
close $metaFd

if {$eqRecNine} {
	set recKind "nine SSI"
} elseif {$eqRecLean} {
	set recKind "center column"
} else {
	set recKind [format "|x|<=%.3g m" $eqWindowX]
}
set recCounts [format "nodes=%d (disp %d) eles=%d quads=%d pileBeams=%d pileSprings=%d" \
	$nGeomNode $nDispNode $nWinEle $nQuadRec \
	[llength $eqPileBeamRows] [llength $eqPileSprRows]]
if {$eqNP <= 1} {
	puts [format "----- EQ recorders  %s  dT=%g s  %s -> %s -----" \
		$recKind $eqRecDt $recCounts $eqOutDir]
	if {$eqRecNine} {
		puts "  SSI stations (iy y layer isTip iSeg):"
		foreach st $eqLeanStations {
			puts [format "    %s" $st]
		}
	} elseif {$eqRecLean} {
		puts [format "  center-pile stations: %d" [llength $eqLeanStations]]
	}
} else {
	puts [format "EQRecorders rank %d: %s  %s" $eqPID $recKind $recCounts]
}
