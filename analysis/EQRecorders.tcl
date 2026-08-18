# analysis/EQRecorders.tcl
# Goals: node/element recorders after `analysis Transient`.
# recordersON: 0 off; 1 full |x|<=eqWindowX window; 2 lean (center pile, x=0 soil).
# Serial (getNP=1): same files as before. OpenSeesMP (getNP>1): local tags;
# files name.$pid after rank-0 mkdir.

set eqNP 1
set eqPID 0
if {[llength [info commands getNP]]} {
	set eqNP [getNP]
	set eqPID [getPID]
}

if {![info exists recordersON]} {
	set recordersON 1
}
if {$recordersON != 0 && $recordersON != 1 && $recordersON != 2} {
	error "EQRecorders.tcl: recordersON must be 0, 1, or 2 (got '$recordersON')"
}

if {$recordersON == 0} {
	if {$eqNP <= 1 || $eqPID == 0} {
		puts "EQRecorders: recordersON=0 -- no node recorders"
	}
	return
}
set eqRecLean [expr {$recordersON == 2}]

if {$eqNP > 1} {
	set eqRecSuf [format ".%d" $eqPID]
} else {
	set eqRecSuf ""
}

# Path under eqOutDir. Serial: $name. MP: $name.$pid
# Args: name (e.g. window_disp.out)
# Returns: abs path (string)
proc eqRecPath {name} {
	global eqOutDir eqRecSuf
	return [file join $eqOutDir ${name}${eqRecSuf}]
}

# Local element tags in [$e0, $e1] inclusive.
# Args: e0 e1 (int)
# Returns: list of tags this rank owns
proc eqRecRange {e0 e1} {
	global eqLocalEle
	set out {}
	for {set e $e0} {$e <= $e1} {incr e} {
		if {[info exists eqLocalEle($e)]} {
			lappend out $e
		}
	}
	return $out
}

# Tags from a list that this rank owns.
# Args: tags (int list)
# Returns: subset
proc eqRecOwned {tags} {
	global eqLocalEle eqNP
	if {$eqNP <= 1} { return $tags }
	set out {}
	foreach e $tags {
		if {[info exists eqLocalEle($e)]} { lappend out $e }
	}
	return $out
}

# True if x is the x=0 soil column (lean window).
# Args: x (m)
# Returns: 1 | 0
proc eqOnCenterline {x} {
	return [expr {abs($x) < 1.0e-6}]
}

# Keep this pile-spring station when recordersON=2: center shaft only.
# Args: ip nPile
# Returns: 1 | 0
proc eqKeepPileSpr {ip nPile} {
	set ipC [expr {($nPile - 1) / 2}]
	return [expr {$ip == $ipC}]
}

# Recorders live under plot/out (same tree as figures). serial | parallel so
# OpenSees and OpenSeesMP dumps do not overwrite each other.
if {![info exists eqOutDir]} {
	if {$eqNP > 1} {
		set eqRunKind "parallel"
	} else {
		set eqRunKind "serial"
	}
	set eleType $soilEleType
	if {![info exists plotDir]} {
		set plotDir [file join $root plot]
	}
	set eqOutDir [file join $plotDir out profile$soilProfile eq \
		$eqRunKind $soilBoundary $eleType $pierEleType]
}
if {$eqNP <= 1} {
	if {[file isdirectory $eqOutDir]} {
		foreach oldOut [glob -nocomplain [file join $eqOutDir *]] {
			if {[file isdirectory $oldOut]} { continue }
			file delete -force $oldOut
		}
	}
	file mkdir $eqOutDir
} else {
	if {$eqPID == 0} {
		if {[file isdirectory $eqOutDir]} {
			foreach oldOut [glob -nocomplain [file join $eqOutDir *]] {
				if {[file isdirectory $oldOut]} { continue }
				file delete -force $oldOut
			}
		}
		file mkdir $eqOutDir
	}
	barrier
}

array unset eqLocalEle
array unset eqLocalNode
foreach e [getEleTags] { set eqLocalEle($e) 1 }
foreach n [getNodeTags] { set eqLocalNode($n) 1 }

array set eqIsStruct {}
if {[info exists structNodeTags]} {
	foreach n $structNodeTags { set eqIsStruct($n) 1 }
}
set nSoilLo -1
set nSoilHi -1
if {[info exists tagShift_soil] && [info exists tagShift_spr]} {
	set nSoilLo $tagShift_soil
	set nSoilHi [expr {$tagShift_spr - 1}]
}

# ---
# 1. WINDOW NODES
# ---
set eqWindowNodeTags {}
set nodesFd [open [eqRecPath window_nodes.txt] w]
if {$eqRecLean} {
	puts $nodesFd "# tag x y  (structure + x=0 soil column; recordersON=2)"
} else {
	puts $nodesFd "# tag x y   (|x|<=$eqWindowX m)"
}
foreach n [lsort -integer [getNodeTags]] {
	set xy [nodeCoord $n]
	set x [lindex $xy 0]
	set y [lindex $xy 1]
	set keep 0
	if {$eqRecLean} {
		if {[info exists eqIsStruct($n)]} {
			set keep 1
		} elseif {$nSoilLo >= 0 && $n >= $nSoilLo && $n <= $nSoilHi} {
			if {[eqOnCenterline $x]} { set keep 1 }
		}
	} elseif {abs($x) <= $eqWindowX + 1.0e-9} {
		set keep 1
	}
	if {!$keep} { continue }
	lappend eqWindowNodeTags $n
	puts $nodesFd [format "%d %.8g %.8g" $n $x $y]
}
close $nodesFd

array set inWindow {}
foreach n $eqWindowNodeTags { set inWindow($n) 1 }
set elesFd [open [eqRecPath window_eles.txt] w]
puts $elesFd {# eleTag n1 n2 [n3 n4...]  (all nodes in window)}
set nWinEle 0
foreach e [getEleTags] {
	if {[catch {set enodes [eleNodes $e]}]} { continue }
	set ok 1
	foreach en $enodes {
		if {![info exists inWindow($en)]} { set ok 0; break }
	}
	if {!$ok} { continue }
	puts $elesFd [concat $e $enodes]
	incr nWinEle
}
close $elesFd

# ---
# 2. WINDOW QUADS
# ---
# 1: all four nodes in |x|<=eqWindowX.  2: column with left face at x=0.
set winQuads {}
set quadsFd [open [eqRecPath window_quads.txt] w]
if {$eqRecLean} {
	puts $quadsFd {# eleTag layer  (centerline column, left face x=0; recordersON=2)}
} else {
	puts $quadsFd {# eleTag layer  (near-field quads; all four nodes in |x|<=eqWindowX)}
}
if {[info exists soilEleTags]} {
	if {$eqRecLean && [info exists soilXs] && [info exists soil_nY]} {
		set ixC -1
		set nXq [llength $soilXs]
		for {set ix 0} {$ix < $nXq} {incr ix} {
			if {abs([lindex $soilXs $ix]) < 1.0e-9} {
				set ixC $ix
				break
			}
		}
		set nRow [expr {$soil_nY - 1}]
		if {$ixC >= 0 && $ixC < $nXq - 1 && $nRow >= 1} {
			for {set iy 0} {$iy < $nRow} {incr iy} {
				set idx [expr {$ixC * $nRow + $iy}]
				if {$idx >= [llength $soilEleTags]} { break }
				set e [lindex $soilEleTags $idx]
				if {![info exists eqLocalEle($e)]} { continue }
				lappend winQuads $e
				set nm "?"
				if {[info exists soilEleLayer($e)]} { set nm $soilEleLayer($e) }
				puts $quadsFd [format "%d %s" $e $nm]
			}
		}
	} else {
		foreach e $soilEleTags {
			if {![info exists eqLocalEle($e)]} { continue }
			if {[catch {set enodes [eleNodes $e]}]} { continue }
			if {[llength $enodes] < 4} { continue }
			set ok 1
			foreach en $enodes {
				if {![info exists inWindow($en)]} { set ok 0; break }
			}
			if {!$ok} { continue }
			lappend winQuads $e
			set nm "?"
			if {[info exists soilEleLayer($e)]} { set nm $soilEleLayer($e) }
			puts $quadsFd [format "%d %s" $e $nm]
		}
	}
}
close $quadsFd

set nodeChunk 250
set windowDispFiles {}
if {[llength $eqWindowNodeTags] < 1} {
	if {$eqNP <= 1} {
		puts "EQRecorders: WARNING no window nodes -- skip window recorder"
	} else {
		puts [format "EQRecorders rank %d: no window nodes -- skip" $eqPID]
	}
} else {
	set iChunk 0
	set tagBatch {}
	set nWin [llength $eqWindowNodeTags]
	foreach n $eqWindowNodeTags {
		lappend tagBatch $n
		if {[llength $tagBatch] >= $nodeChunk} {
			if {$nWin <= $nodeChunk} {
				set outPath [eqRecPath window_disp.out]
			} else {
				set outPath [eqRecPath [format "window_disp_%02d.out" $iChunk]]
			}
			# recorder Node -file $fileName -time -node $nodeTags -dof 1 2 disp
			eval recorder Node -file $outPath -time -node $tagBatch -dof 1 2 disp
			lappend windowDispFiles [file tail $outPath]
			incr iChunk
			set tagBatch {}
		}
	}
	if {[llength $tagBatch] > 0} {
		if {$nWin <= $nodeChunk && $iChunk == 0} {
			set outPath [eqRecPath window_disp.out]
		} else {
			set outPath [eqRecPath [format "window_disp_%02d.out" $iChunk]]
		}
		eval recorder Node -file $outPath -time -node $tagBatch -dof 1 2 disp
		lappend windowDispFiles [file tail $outPath]
		incr iChunk
	}
}

# Quad (4 Gauss pts) or SSPquad (1 IP at centroid).
# Quad wiki: stresses/strains. PIMY 2D getStress is (sxx, syy, sxy) so sxy = tau_xy;
# getStrain is (exx, eyy, gxy) so gxy = gamma_xy. Same packing PDMY02 uses in 2D.
# material $ip stress would add sigma_zz and eta_r; we keep the 3-comp element query.
set quadStressFiles {}
set quadStrainFiles {}
set nGPq 0
set sigRsp ""
set epsRsp ""
if {[llength $winQuads] >= 1} {
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
	set nWinQuad [llength $winQuads]
	set iC 0
	set tagBatch {}
	foreach e $winQuads {
		lappend tagBatch $e
		if {[llength $tagBatch] >= $quadChunk} {
			if {$nWinQuad <= $quadChunk && $iC == 0} {
				set stressPath [eqRecPath window_quad_stress.out]
				set strainPath [eqRecPath window_quad_strain.out]
			} else {
				set stressPath [eqRecPath [format "window_quad_stress_%02d.out" $iC]]
				set strainPath [eqRecPath [format "window_quad_strain_%02d.out" $iC]]
			}
			eval recorder Element -file $stressPath -time -ele $tagBatch $sigRsp
			eval recorder Element -file $strainPath -time -ele $tagBatch $epsRsp
			lappend quadStressFiles [file tail $stressPath]
			lappend quadStrainFiles [file tail $strainPath]
			incr iC
			set tagBatch {}
		}
	}
	if {[llength $tagBatch] > 0} {
		if {$nWinQuad <= $quadChunk && $iC == 0} {
			set stressPath [eqRecPath window_quad_stress.out]
			set strainPath [eqRecPath window_quad_strain.out]
		} else {
			set stressPath [eqRecPath [format "window_quad_stress_%02d.out" $iC]]
			set strainPath [eqRecPath [format "window_quad_strain_%02d.out" $iC]]
		}
		eval recorder Element -file $stressPath -time -ele $tagBatch $sigRsp
		eval recorder Element -file $strainPath -time -ele $tagBatch $epsRsp
		lappend quadStressFiles [file tail $stressPath]
		lappend quadStrainFiles [file tail $strainPath]
		incr iC
	}
}

set nWinQ 0
if {[info exists nWinQuad]} { set nWinQ $nWinQuad }
if {[llength $eqWindowNodeTags] >= 1} {
	if {$eqRecLean} {
		set recKind "lean"
	} else {
		set recKind [format "|x|<=%.3g m" $eqWindowX]
	}
	if {$eqNP <= 1} {
		puts [format "----- EQ recorders  %s  nodes=%d  eles=%d  quads=%d -> %s -----" \
			$recKind [llength $eqWindowNodeTags] $nWinEle $nWinQ $eqOutDir]
	} else {
		puts [format "EQRecorders rank %d: %s  nodes=%d  eles=%d  quads=%d -> %s" \
			$eqPID $recKind [llength $eqWindowNodeTags] $nWinEle $nWinQ $eqOutDir]
	}
}

set metaFd [open [eqRecPath window_meta.txt] w]
puts $metaFd "soilProfile $soilProfile"
puts $metaFd "soilBoundary $soilBoundary"
puts $metaFd "soilEleType $soilEleType"
puts $metaFd "recordersON $recordersON"
puts $metaFd "eqWindowX $eqWindowX"
puts $metaFd "nWindowNodes [llength $eqWindowNodeTags]"
puts $metaFd "nWindowEles $nWinEle"
puts $metaFd "nWindowQuads [llength $winQuads]"
puts $metaFd "dtAnalysis $dtAnalysis"
puts $metaFd "eqNsteps $eqNsteps"
if {[info exists fvNsteps]} {
	puts $metaFd "freeVibT $eqFreeVibT"
	puts $metaFd "freeVibNsteps $fvNsteps"
	puts $metaFd "eqNstepsAll $eqNstepsAll"
	puts $metaFd "Trec [expr {$Trec + $eqFreeVibT}]"
} else {
	puts $metaFd "Trec $Trec"
}
puts $metaFd "gmVelFile $gmVelFile"
puts $metaFd "gmScaleFactor $gmScaleFactor"
puts $metaFd "dispFiles $windowDispFiles"
puts $metaFd "dispFormat each file: time then (ux uy) for that chunk; chunks follow window_nodes.txt order"
if {$eqNP > 1} {
	puts $metaFd "pid $eqPID"
	puts $metaFd "np $eqNP"
	puts $metaFd "fileSuffix $eqRecSuf"
}
if {[llength $quadStressFiles] > 0} {
	puts $metaFd "quadEleFile window_quads.txt"
	puts $metaFd "quadNgp $nGPq"
	puts $metaFd "quadStressRsp $sigRsp"
	puts $metaFd "quadStrainRsp $epsRsp"
	puts $metaFd "quadStressFiles $quadStressFiles"
	puts $metaFd "quadStrainFiles $quadStrainFiles"
	puts $metaFd "quadStressFormat time then nGP*(sxx syy sxy) per ele; sxy=tau_xy; ele order window_quads.txt"
	puts $metaFd "quadStrainFormat time then nGP*(exx eyy gxy) per ele; gxy=gamma_xy"
}
puts $metaFd "pierHinge $pierEleType"
puts $metaFd "pileEleType $pileEleType"
if {$pierEleType eq "lumpedPlasticity" || $pierEleType eq "forceBeamColumn"} {
	puts $metaFd "hingeForceFile pier_hinge_force.out"
	puts $metaFd "hingeDefoFile pier_hinge_defo.out"
}
if {$pierEleType eq "lumpedPlasticity"} {
	puts $metaFd "hingeTopForceFile pier_hinge_top_force.out"
	puts $metaFd "hingeTopDefoFile pier_hinge_top_defo.out"
}
if {[info exists eleTag_pile_base] && [info exists eleTag_pile_last] \
		&& $eleTag_pile_last >= $eleTag_pile_base} {
	puts $metaFd "pileBeamGlobalForceFile pile_beam_globalForce.out"
	if {$pileEleType eq "dispBeamColumn"} {
		puts $metaFd "pileBeamSec1DefoFile pile_beam_sec1_defo.out"
		if {[info exists nIP_pile]} {
			puts $metaFd "nIP_pile $nIP_pile"
		}
	}
}
if {[info exists nSprings] && $nSprings >= 1} {
	puts $metaFd "pileSpringsForceFile pile_springs_force.out"
	puts $metaFd "pileSpringsDefoFile pile_springs_defo.out"
	puts $metaFd "capSpringsForceFile cap_springs_force.out"
	puts $metaFd "capSpringsDefoFile cap_springs_defo.out"
	if {[info exists nSoffit] && $nSoffit > 0} {
		puts $metaFd "capSoffitForceFile cap_springs_soffit_force.out"
		puts $metaFd "capSoffitDefoFile cap_springs_soffit_defo.out"
	}
}
close $metaFd

if {$eqNP <= 1 || [info exists eqLocalNode($nodeTag_pierTop_deckBC)]} {
	recorder Node -file [eqRecPath pier_top_disp.out] -time \
		-node $nodeTag_pierTop_deckBC -dof 1 2 3 disp
	recorder Node -file [eqRecPath pier_top_acc.out] -time \
		-node $nodeTag_pierTop_deckBC -dof 1 2 3 accel
}
if {[info exists nPrimary]} {
	if {$eqNP <= 1 || [info exists eqLocalNode($nPrimary)]} {
		recorder Node -file [eqRecPath soil_base_primary.out] -time \
			-node $nPrimary -dof 1 2 disp
	}
}

# Base hinge: ZLS (lumpedPlasticity) or forceBeamColumn i-end IP (section 1).
# 2D Fiber: force = P Mz; deformation = eps kappa.
# ZLS: those are axial displacement (m) and rotation theta (rad).
if {$pierEleType eq "lumpedPlasticity"} {
	if {$eqNP <= 1 || [info exists eqLocalEle($eleTag_pier_botSpr)]} {
		recorder Element -file [eqRecPath pier_hinge_force.out] -time \
			-ele $eleTag_pier_botSpr section force
		recorder Element -file [eqRecPath pier_hinge_defo.out] -time \
			-ele $eleTag_pier_botSpr section deformation
	}
	if {$eqNP <= 1 || [info exists eqLocalEle($eleTag_pier_topSpr)]} {
		recorder Element -file [eqRecPath pier_hinge_top_force.out] -time \
			-ele $eleTag_pier_topSpr section force
		recorder Element -file [eqRecPath pier_hinge_top_defo.out] -time \
			-ele $eleTag_pier_topSpr section deformation
	}
} elseif {$pierEleType eq "forceBeamColumn"} {
	if {$eqNP <= 1 || [info exists eqLocalEle($eleTag_pier)]} {
		recorder Element -file [eqRecPath pier_hinge_force.out] -time \
			-ele $eleTag_pier section 1 force
		recorder Element -file [eqRecPath pier_hinge_defo.out] -time \
			-ele $eleTag_pier section 1 deformation
	}
}

# Pile shafts: global end forces (2D: Px Py Mz at i and j).
# dispBeamColumn: section 1 = first IP (Lobatto i-end: eps, kappa).
# recordersON=2: center pile all segments (no outer piles).
if {[info exists eleTag_pile_base] && [info exists eleTag_pile_last] \
		&& $eleTag_pile_last >= $eleTag_pile_base} {
	set pileBeamFd [open [eqRecPath pile_beam_eles.txt] w]
	puts $pileBeamFd "# eleTag ip iy   (section 1 = first IP = i-end)"
	set eleCur [expr {$eleTag_pile_base - 1}]
	set pileBeamKeep {}
	set ipC [expr {($n_pile - 1) / 2}]
	for {set iPile 0} {$iPile < $n_pile} {incr iPile} {
		for {set iSeg 1} {$iSeg <= $nSeg_pile} {incr iSeg} {
			incr eleCur
			if {$eqRecLean && $iPile != $ipC} { continue }
			lappend pileBeamKeep $eleCur
			if {$eqNP <= 1 || [info exists eqLocalEle($eleCur)]} {
				puts $pileBeamFd [format "%d %d %d" $eleCur $iPile $iSeg]
			}
		}
	}
	close $pileBeamFd
	if {!$eqRecLean && $eleCur != $eleTag_pile_last} {
		puts [format "EQRecorders: WARNING pile_beam_eles last=%d  eleTag_pile_last=%d" \
			$eleCur $eleTag_pile_last]
	}
	set pileLoc [eqRecOwned $pileBeamKeep]
	if {[llength $pileLoc] >= 1} {
		eval recorder Element -file [eqRecPath pile_beam_globalForce.out] -time \
			-ele $pileLoc globalForce
		if {$pileEleType eq "dispBeamColumn"} {
			eval recorder Element -file [eqRecPath pile_beam_sec1_defo.out] -time \
				-ele $pileLoc section 1 deformation
		}
	}
}

# Interface zeroLength: dir 1 = p-y, dir 2 = t-z (localForce / deformation).
# Soffit q-z is dir 2 only -- own files so column counts stay even.
if {[info exists nSprings] && $nSprings >= 1} {
	set nPileSpr 0
	if {[info exists nPileSprings]} { set nPileSpr $nPileSprings }
	set nSof 0
	if {[info exists nSoffit]} { set nSof $nSoffit }
	set ePile0 $eleTag_spr_base
	set ePile1 [expr {$ePile0 + $nPileSpr - 1}]
	set eCap0 [expr {$ePile1 + 1}]
	set eCapFace [expr {$eleTag_spr_last - $nSof}]
	set eSof0 [expr {$eCapFace + 1}]
	set eSof1 $eleTag_spr_last

	array unset eqKeepPileSprEle
	array unset eqPileSprIp
	array unset eqPileSprIy
	set pileSprKeep {}
	if {[info exists pileSprEleRec]} {
		foreach rec $pileSprEleRec {
			lassign $rec e ip iyP isTip
			if {$eqRecLean && ![eqKeepPileSpr $ip $n_pile]} {
				continue
			}
			set eqKeepPileSprEle($e) 1
			set eqPileSprIp($e) $ip
			set eqPileSprIy($e) $iyP
			lappend pileSprKeep $e
		}
	}

	if {[info exists ssiSpringDump]} {
		set pileSprFd [open [eqRecPath pile_springs_eles.txt] w]
		set capSprFd [open [eqRecPath cap_springs_eles.txt] w]
		puts $pileSprFd "# eleTag  kind  ip  iy"
		puts $capSprFd "# eleTag  kind  (cap face 2-dir; cap_soffit q-z dir 2)"
		array unset sprSeen
		foreach rec $ssiSpringDump {
			set e [lindex $rec 0]
			if {$e < 0 || [info exists sprSeen($e)]} { continue }
			if {$eqNP > 1 && ![info exists eqLocalEle($e)]} { continue }
			set sprSeen($e) 1
			set kind [lindex $rec 7]
			if {$kind eq "pile"} {
				if {[info exists eqKeepPileSprEle] \
						&& [array size eqKeepPileSprEle] > 0 \
						&& ![info exists eqKeepPileSprEle($e)]} {
					continue
				}
				set ipW -1
				set iyW -1
				if {[info exists eqPileSprIp($e)]} {
					set ipW $eqPileSprIp($e)
					set iyW $eqPileSprIy($e)
				}
				puts $pileSprFd [format "%d %s %d %d" $e $kind $ipW $iyW]
			} else {
				puts $capSprFd [format "%d %s" $e $kind]
			}
		}
		close $pileSprFd
		close $capSprFd
	}

	if {$nPileSpr >= 1} {
		if {$eqRecLean && [llength $pileSprKeep] >= 1} {
			set pileSprLoc [eqRecOwned $pileSprKeep]
		} else {
			set pileSprLoc [eqRecRange $ePile0 $ePile1]
		}
		if {[llength $pileSprLoc] >= 1} {
			eval recorder Element -file [eqRecPath pile_springs_force.out] -time \
				-ele $pileSprLoc localForce
			eval recorder Element -file [eqRecPath pile_springs_defo.out] -time \
				-ele $pileSprLoc deformation
		}
	}
	if {$eCapFace >= $eCap0} {
		if {$eqNP <= 1} {
			recorder Element -file [eqRecPath cap_springs_force.out] -time \
				-eleRange $eCap0 $eCapFace localForce
			recorder Element -file [eqRecPath cap_springs_defo.out] -time \
				-eleRange $eCap0 $eCapFace deformation
		} else {
			set capSprLoc [eqRecRange $eCap0 $eCapFace]
			if {[llength $capSprLoc] >= 1} {
				eval recorder Element -file [eqRecPath cap_springs_force.out] -time \
					-ele $capSprLoc localForce
				eval recorder Element -file [eqRecPath cap_springs_defo.out] -time \
					-ele $capSprLoc deformation
			}
		}
	}
	if {$nSof > 0} {
		if {$eqNP <= 1} {
			recorder Element -file [eqRecPath cap_springs_soffit_force.out] -time \
				-eleRange $eSof0 $eSof1 localForce
			recorder Element -file [eqRecPath cap_springs_soffit_defo.out] -time \
				-eleRange $eSof0 $eSof1 deformation
		} else {
			set sofSprLoc [eqRecRange $eSof0 $eSof1]
			if {[llength $sofSprLoc] >= 1} {
				eval recorder Element -file [eqRecPath cap_springs_soffit_force.out] -time \
					-ele $sofSprLoc localForce
				eval recorder Element -file [eqRecPath cap_springs_soffit_defo.out] -time \
					-ele $sofSprLoc deformation
			}
		}
	}
}

record
