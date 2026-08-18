# soil/BuildSoilSprings.tcl
# Goals: pile and cap p-y / t-z / q-z. Call AFTER soil mesh (ndf=3; dups use ndf=2).
# Knobs: Parameters.tcl. Fold + re-tie: analysis/FoldStructNodes.tcl.
#
# Layout:
#   helpers (lookup, API pu/tu/qu)
#   Pass 1: full-station capacities
#   Pass 2: ndf=2 dups
#   pile shaft + tip springs
#   cap face + soffit
#
# =====================================================================
# 3. MATERIALS AND SECTIONS
# =====================================================================
# =====================================================================
# 4. ELEMENTS
# =====================================================================
# =====================================================================
# 5. BOUNDARY CONDITIONS / CONSTRAINTS
# =====================================================================

if {![info exists soil_nX] || ![info exists Ge_pile]} {
	error "BuildSoilSprings.tcl: source BuildSoilMesh.tcl / Parameters first"
}
if {![info exists nSoilRows] || ![info exists soilRowLayer]} {
	error "BuildSoilSprings.tcl: soil rows from BuildSoilMaterials.tcl required"
}
if {![info exists nodeTag_pile_tips]} {
	error "BuildSoilSprings.tcl: source BuildPilesNodes.tcl first"
}
if {$pileSpring ne "inelastic" && $pileSpring ne "elastic" && $pileSpring ne "none"} {
	error "BuildSoilSprings.tcl: pileSpring must be inelastic|elastic|none (got '$pileSpring')"
}

set gravSoilLockPairs {}
set springEqualDOFPairs {}
array unset gravLockSeen
# Record one soil <-> pile/cap pair for SoilGravity (soil follows retain).
# Args:    soilNd retainNd
# Returns: none (skips invalid or duplicate retainNd)
proc gravLockPush {soilNd rNd} {
	global gravSoilLockPairs gravLockSeen
	if {$soilNd < 0 || $rNd < 0 || $soilNd == $rNd} { return }
	if {[info exists gravLockSeen($rNd)]} { return }
	set gravLockSeen($rNd) 1
	lappend gravSoilLockPairs [list $soilNd $rNd]
}
# pile/cap retained; dup (spring soil-side node) follows UX,UY
# equalDOF $rNodeTag $cNodeTag $dof1 $dof2
# Args:    rNd (retain)  dup  soilNd
proc springTie {rNd dup soilNd} {
	global springEqualDOFPairs
	equalDOF $rNd $dup 1 2
	lappend springEqualDOFPairs [list $rNd $dup $soilNd]
}

set pi_s 3.141592653589793
set pci_to_Pa_m [expr {271.447 / 0.0254}]

# ---- helpers (ix,iy 0-based; node tag = nodeTag_soil_base + ix*100 + iy) ----
proc soilIxAtX {xT} {
	global soilXs
	set best 0
	set bd 1.0e99
	set i 0
	foreach x $soilXs {
		set d [expr {abs($x - $xT)}]
		if {$d < $bd} { set bd $d; set best $i }
		incr i
	}
	return $best
}

proc soilIyAtY {yT} {
	global soilYs
	set best 0
	set bd 1.0e99
	set i 0
	foreach y $soilYs {
		set d [expr {abs($y - $yT)}]
		if {$d < $bd} { set bd $d; set best $i }
		incr i
	}
	return $best
}

proc soilNodeAtIxIy {ix iy} {
	global nodeTag_soil_base
	return [expr {$nodeTag_soil_base + $ix*100 + $iy}]
}

# Soil node at (x, y) m, or -1 if that tag is missing.
proc soilNdAt {x y} {
	set ix [soilIxAtX $x]
	set iy [soilIyAtY $y]
	set nd [soilNodeAtIxIy $ix $iy]
	if {[lsearch -exact [getNodeTags] $nd] < 0} { return -1 }
	return $nd
}

proc sigVprime {yMid} {
	global soilYs nSoilRows soilRho soilIsSand rho_w gravity_accel
	set sig 0.0
	for {set iy 0} {$iy < $nSoilRows} {incr iy} {
		set yt [lindex $soilYs $iy]
		set yb [lindex $soilYs [expr {$iy + 1}]]
		if {$yMid >= $yt - 1.0e-9} { break }
		set yLo $yb
		if {$yMid > $yb} { set yLo $yMid }
		set rho $soilRho($iy)
		if {$soilIsSand($iy)} {
			set gam [expr {($rho - $rho_w)*$gravity_accel}]
		} else {
			set gam [expr {$rho*$gravity_accel}]
		}
		set yA [expr {($yt < 0) ? $yt : 0.0}]
		if {$yLo >= $yA} { continue }
		set th [expr {$yA - $yLo}]
		if {$th > 0} {
			set sig [expr {$sig + $gam*$th}]
		}
	}
	if {$sig < 1.0} { set sig 1.0 }
	return $sig
}

proc puClay {cu D J zDepth gam} {
	set p1 [expr {(3.0*$cu + $gam*$zDepth + $J*$cu*$zDepth/$D)*$D}]
	set p2 [expr {9.0*$cu*$D}]
	return [expr {($p1 < $p2) ? $p1 : $p2}]
}

proc puSand {phi D zDepth sigV k_Pa_m} {
	if {$phi < 25} {
		set C1 1.0; set C2 1.5; set C3 10.0
	} elseif {$phi < 30} {
		set C1 1.5; set C2 2.0; set C3 20.0
	} elseif {$phi < 35} {
		set C1 2.0; set C2 2.5; set C3 30.0
	} elseif {$phi < 40} {
		set C1 2.5; set C2 3.5; set C3 50.0
	} else {
		set C1 3.0; set C2 4.0; set C3 70.0
	}
	set gamP [expr {$sigV/max($zDepth,0.1)}]
	set ps [expr {($C1*$zDepth + $C2*$D)*$gamP*$zDepth}]
	set pd [expr {$C3*$D*$gamP*$zDepth}]
	return [expr {($ps < $pd) ? $ps : $pd}]
}

proc y50Clay {D} {
	return [expr {2.5*0.01*$D}]
}

proc y50Sand {pu k_Pa_m zDepth} {
	set kz [expr {$k_Pa_m*max($zDepth,0.1)}]
	return [expr {0.549*$pu/$kz}]
}

proc tuClay {cu D trib} {
	global pi_s
	return [expr {$cu*$pi_s*$D*$trib}]
}

proc tuSand {sigV D trib phi} {
	global pi_s
	set delta [expr {30.0*$pi_s/180.0}]
	set fmax 95.8e3
	set K 0.8
	if {$phi >= 39.0} { set K 1.25 }
	set f [expr {$K*$sigV*tan($delta)}]
	if {$f > $fmax} { set f $fmax }
	return [expr {$f*$pi_s*$D*$trib}]
}

proc quClay {cu D} {
	global pi_s
	return [expr {9.0*$cu*$pi_s*0.25*$D*$D}]
}

proc quSand {sigV D phi} {
	global pi_s
	set Area [expr {$pi_s*0.25*$D*$D}]
	set q [expr {40.0*$sigV}]
	if {$q > 9580.0e3} { set q 9580.0e3 }
	return [expr {$q*$Area}]
}

# Pass 1: full-station capacities
set pileHeads [list \
	[list $nodeTag_cap_BL [expr {-$s_pile_cap}]] \
	[list $nodeTag_cap_BC 0.0] \
	[list $nodeTag_cap_BR $s_pile_cap] \
]

set maxKpy 1.0e6
set springRec {}
for {set ip 0} {$ip < $n_pile} {incr ip} {
	lassign [lindex $pileHeads $ip] headTag xP
	for {set iyP 0} {$iyP <= $nSeg_pile} {incr iyP} {
		if {$iyP == 0} {
			set pileNd $headTag
			set y [expr {-$H_cap}]
			set trib [expr {0.5*$dy_soil}]
			set isTip 0
		} elseif {$iyP == $nSeg_pile} {
			set pileNd [expr {$nodeTag_pile_base + $ip*100 + $iyP}]
			set y [expr {-$H_cap - $iyP*$dy_soil}]
			set trib [expr {0.5*$dy_soil}]
			set isTip 1
		} else {
			set pileNd [expr {$nodeTag_pile_base + $ip*100 + $iyP}]
			set y [expr {-$H_cap - $iyP*$dy_soil}]
			set trib $dy_soil
			set isTip 0
		}
		set iyS [soilRowAtY $y]
		set nm $soilRowLayer($iyS)
		set zDepth [expr {max(-$y, 0.05)}]
		set sig [sigVprime $y]
		set D $D_pile
		if {$soilIsSand($iyS)} {
			set k [expr {$soilKpci($iyS)*$pci_to_Pa_m}]
			set puL [puSand $soilPhi($iyS) $D $zDepth $sig $k]
			set y50 [y50Sand $puL $k $zDepth]
			set tult [expr {[tuSand $sig $D $trib $soilPhi($iyS)]*$n_pile_row}]
			set tzType 2
			set useLiq [expr {$soilProfile == 1 && $nm eq "L3"}]
			set pyType 2
		} else {
			set cu $soilC($iyS)
			set gam [expr {$soilRho($iyS)*$gravity_accel}]
			set puL [puClay $cu $D 0.5 $zDepth $gam]
			set y50 [y50Clay $D]
			set tult [expr {[tuClay $cu $D $trib]*$n_pile_row}]
			set tzType 1
			set useLiq 0
			set pyType 2
			set k 0.0
		}
		if {$y50 < 1.0e-5} { set y50 1.0e-5 }
		set pult [expr {$puL*$trib*$Ge_pile*$n_pile_row}]
		set Kpy [expr {$pult/$y50}]
		if {$Kpy > $maxKpy} { set maxKpy $Kpy }

		set qult 0.0
		set z50q $z50_tz
		set qzType 1
		if {$isTip} {
			if {$soilIsSand($iyS)} {
				set qult [expr {[quSand $sig $D $soilPhi($iyS)]*$n_pile_row}]
				set z50q [expr {0.01*$D}]
				set qzType 2
			} else {
				set qult [expr {[quClay $soilC($iyS) $D]*$n_pile_row}]
				set z50q [expr {0.00625*$D}]
				set qzType 1
			}
		}
		lappend springRec [list $ip $iyP $headTag $xP $pileNd $y $trib $isTip \
			$nm $pult $y50 $pyType $tult $tzType $useLiq $qult $z50q $qzType]
	}
}

# Pass 2: nodes + materials + elements
set eSpr [expr {$eleTag_spr_base - 1}]
set nSpr 0
set ssiSpringDump {}
set pileSprEleRec {}
set pileSpringPropsDump {}
set capFacePropsDump {}
set capSoffitPropsDump {}
set liqSpringMatTags {}

# ---- pileSpring none: equalDOF only ----
if {$pileSpring eq "none"} {
	foreach rec $springRec {
		lassign $rec ip iyP headTag xP pileNd y trib isTip \
			nm pult y50 pyType tult tzType useLiq qult z50q qzType
		set soilNd [soilNdAt $xP $y]
		if {$soilNd < 0} {
			puts "WARNING: no soil node at pile ip=$ip y=$y -- skip"
			continue
		}
		equalDOF $soilNd $pileNd 1 2
		set gravLockSeen($pileNd) 1
		set xyS [nodeCoord $soilNd]
		set xs [lindex $xyS 0]
		set ys [lindex $xyS 1]
		if {$isTip} {
			set latType "none"; set axType "none"
			set tAx $qult; set z50Ax $z50q; set pRes $pult; set tRes $qult
		} else {
			set latType "none"; set axType "none"
			set tAx $tult; set z50Ax $z50_tz; set pRes $pult; set tRes $tult
		}
		lappend ssiSpringDump [list -1 $xs $ys $xP $y $xP $y "pile" $latType]
		lappend ssiSpringDump [list -1 $xs $ys $xP $y $xP $y "pile" $axType]
		lappend pileSpringPropsDump [list $ip $iyP $y [expr {-$y}] $nm $isTip 0 \
			$pult $pRes $tAx $tRes $y50 $z50Ax $latType $axType $trib 0.0]
	}
}

# Cap capacities (needed before work lists)
# Face p-y: Mokwa one-face wall (b = L_cap OOP); bilateral -> 1/2 per side, then x h/H.
# Face t-z: full vertical skin c*2H(W+L); x h/(2H) onto the six face springs (OOP lumped).
# Soffit q-z: 9 c W L, split by tributary width on meshed bottom chord.
set cu2 $soilC(L2)
set Hcap $H_cap
set Lcap $L_cap
set Wcap $W_cap
set gam2 [expr {($soilRho(L2) - $rho_w)*$gravity_accel}]
set alpha $alpha_cap_py
# Mokwa phi=0 wall; width = out-of-plane length
set PultCap [expr {$cu2*$Lcap*$Hcap/2.0*(4.0 + $gam2*$Hcap/$cu2 + 0.25*$Hcap/$Lcap + 2.0*$alpha)}]
set Aside [expr {2.0*$Hcap*($Wcap + $Lcap)}]
set TultCap [expr {$cu2*$Aside}]
set QultSoffit [expr {9.0*$cu2*$Wcap*$Lcap}]
# Tributary heights at grade / mid / soffit
set hCapTop [expr {0.25*$Hcap}]
set hCapMid [expr {0.5*$Hcap}]
set hCapBot [expr {0.25*$Hcap}]

# ---- Build work lists + all ndf=2 dups in one batch ----
set sprWork {}
set dupNodes {}
if {$pileSpring ne "none"} {
	foreach rec $springRec {
		lassign $rec ip iyP headTag xP pileNd y trib isTip \
			nm pult y50 pyType tult tzType useLiq qult z50q qzType
		set soilNd [soilNdAt $xP $y]
		if {$soilNd < 0} {
			puts "WARNING: no soil at pile ip=$ip y=$y -- skip"
			continue
		}
		set dup [expr {$nodeTag_sprSoil_base + $ip*100 + $iyP}]
		lappend dupNodes [list $dup $xP $y]
		lappend sprWork [list $ip $iyP $xP $pileNd $y $trib $isTip $nm \
			$pult $y50 $pyType $tult $tzType $useLiq $qult $z50q $qzType \
			$soilNd $dup]
	}
}

# Face stations: [capNd x y hTrib] -- pult = PultCap*h/(2H), tult = TultCap*h/(2H)
# On the outer pile axes (frame edge = shaft centerline)
set capFaceNodes [list \
	[list $nodeTag_cap_TL [expr {-$s_pile_cap}] 0.0 $hCapTop] \
	[list $nodeTag_cap_ML [expr {-$s_pile_cap}] [expr {-0.5*$H_cap}] $hCapMid] \
	[list $nodeTag_cap_BL [expr {-$s_pile_cap}] [expr {-$H_cap}] $hCapBot] \
	[list $nodeTag_cap_TR $s_pile_cap 0.0 $hCapTop] \
	[list $nodeTag_cap_MR $s_pile_cap [expr {-0.5*$H_cap}] $hCapMid] \
	[list $nodeTag_cap_BR $s_pile_cap [expr {-$H_cap}] $hCapBot] \
]
set capFaceWork {}
set iCap 0
foreach row $capFaceNodes {
	lassign $row capNd xC yC hTrib
	set pOne [expr {$PultCap * $hTrib / (2.0*$Hcap)}]
	set tOne [expr {$TultCap * $hTrib / (2.0*$Hcap)}]
	set ix [soilIxAtX $xC]
	set iy [soilIyAtY $yC]
	set soilNd [soilNodeAtIxIy $ix $iy]
	if {[lsearch -exact [getNodeTags] $soilNd] < 0} {
		puts "WARNING: no soil for cap spring at $capNd (x=$xC y=$yC) -- skip"
		incr iCap
		continue
	}
	if {$pileSpring eq "none"} {
		equalDOF $soilNd $capNd 1 2
		set gravLockSeen($capNd) 1
		set xy [nodeCoord $soilNd]
		lappend ssiSpringDump [list -1 \
			[lindex $xy 0] [lindex $xy 1] $xC $yC $xC $yC "cap" "none"]
		lappend ssiSpringDump [list -1 \
			[lindex $xy 0] [lindex $xy 1] $xC $yC $xC $yC "cap" "none"]
		incr iCap
		continue
	}
	set dup [expr {$nodeTag_sprSoil_base + 900 + $iCap}]
	lappend dupNodes [list $dup $xC $yC]
	lappend capFaceWork [list $capNd $xC $yC $soilNd $dup $pOne $tOne]
	incr iCap
}

# Soffit q-z: meshed bottom chord; Qult = 9 c W_cap L_cap, split by tributary length
if {![info exists capSoffitStations]} {
	error "BuildSoilSprings.tcl: BuildPileCapNodes.tcl must export capSoffitStations"
}
set sofUse {}
foreach row $capSoffitStations {
	lassign $row capNd xS
	lappend sofUse [list $capNd $xS]
}
if {[llength $sofUse] == 0} {
	error "BuildSoilSprings.tcl: no soffit soil stations (check capSoffitStations)"
}
set sofXs {}
foreach row $sofUse {
	lappend sofXs [lindex $row 1]
}
set sofTrib {}
set nSof [llength $sofXs]
for {set i 0} {$i < $nSof} {incr i} {
	set xi [lindex $sofXs $i]
	if {$nSof == 1} {
		set trib 1.0
	} elseif {$i == 0} {
		set trib [expr {0.5*([lindex $sofXs 1] - $xi)}]
	} elseif {$i == $nSof - 1} {
		set trib [expr {0.5*($xi - [lindex $sofXs [expr {$i - 1}]])}]
	} else {
		set trib [expr {0.5*([lindex $sofXs [expr {$i + 1}]] - [lindex $sofXs [expr {$i - 1}]])}]
	}
	lappend sofTrib $trib
}
set sofTribSum 0.0
foreach t $sofTrib { set sofTribSum [expr {$sofTribSum + $t}] }

set capSoffitWork {}
set iSof 0
foreach row $sofUse {
	lassign $row capNd xS
	set yS [expr {-$H_cap}]
	set trib [lindex $sofTrib $iSof]
	set qOne [expr {$QultSoffit * $trib / $sofTribSum}]
	set ix [soilIxAtX $xS]
	set iy [soilIyAtY $yS]
	set soilNd [soilNodeAtIxIy $ix $iy]
	if {[lsearch -exact [getNodeTags] $soilNd] < 0} {
		puts "WARNING: no soil for soffit q-z at $capNd (x=$xS) -- skip"
		incr iSof
		continue
	}
	if {$pileSpring eq "none"} {
		incr iSof
		continue
	}
	set dup [expr {$nodeTag_sprSoil_base + 920 + $iSof}]
	lappend dupNodes [list $dup $xS $yS]
	lappend capSoffitWork [list $capNd $xS $yS $soilNd $dup $qOne]
	incr iSof
}

# One ndf=2 switch for all dups, then ndf=3 for ZL
if {[llength $dupNodes] > 0} {
	model BasicBuilder -ndm 2 -ndf 2
	foreach dn $dupNodes {
		lassign $dn dup xD yD
		node $dup $xD $yD
	}
	model BasicBuilder -ndm 2 -ndf 3
}

# ---- pile shaft springs ----
foreach w $sprWork {
	lassign $w ip iyP xP pileNd y trib isTip nm \
		pultS y50 pyType tultS tzType useLiq qultS z50q qzType \
		soilNd dup
	set hasQz $isTip

	springTie $pileNd $dup $soilNd
	gravLockPush $soilNd $pileNd

	set mPy [expr {$matTag_py_base + $nSpr}]
	set mTz [expr {$matTag_tz_base + $nSpr}]
	incr eSpr

	if {$isTip} {
		set tAx $qultS
		set z50Ax $z50q
		set pRes $pultS
		set tRes $qultS
		if {$pileSpring eq "elastic"} {
			set kP [expr {$pultS/$y50}]
			uniaxialMaterial Elastic $mPy $kP
			set latType "py_elastic"
			if {$hasQz} {
				set kQ [expr {$qultS/$z50q}]
				uniaxialMaterial Elastic $mTz $kQ
				set axType "qz_elastic"
			} else {
				set axType "none"
			}
		} else {
			uniaxialMaterial PySimple1 $mPy $pyType $pultS $y50 $Cd_py $c_dash_py
			set latType "py"
			if {$hasQz} {
				uniaxialMaterial QzSimple1 $mTz $qzType $qultS $z50q 0.0 $c_dash_py
				set axType "qz"
			} else {
				set axType "none"
			}
		}
		set useLiqS 0
		if {$hasQz} {
			element zeroLength $eSpr $soilNd $dup -mat $mPy $mTz -dir 1 2 -doRayleigh 1
		} else {
			element zeroLength $eSpr $soilNd $dup -mat $mPy -dir 1 -doRayleigh 1
		}
	} elseif {$pileSpring eq "elastic"} {
		set tAx $tultS
		set z50Ax $z50_tz
		set pRes $pultS
		set tRes $tultS
		set kP [expr {$pultS/$y50}]
		set kT [expr {$tultS/$z50_tz}]
		uniaxialMaterial Elastic $mPy $kP
		uniaxialMaterial Elastic $mTz $kT
		set latType "py_elastic"
		set axType "tz_elastic"
		set useLiqS 0
		element zeroLength $eSpr $soilNd $dup -mat $mPy $mTz -dir 1 2 -doRayleigh 1
	} elseif {$useLiq} {
		set tAx $tultS
		set z50Ax $z50_tz
		set pRes [expr {$pRes_frac*$pultS}]
		set tRes [expr {$pRes_frac*$tultS}]
		set sEle1 [lindex $soilEleTags 0]
		set sEle2 $sEle1
		foreach ee $soilEleTags {
			if {$soilEleLayer($ee) eq $nm} {
				set sEle1 $ee
				set sEle2 $ee
				break
			}
		}
		uniaxialMaterial PyLiq1 $mPy $pyType $pultS $y50 $Cd_py $c_dash_py $pRes $sEle1 $sEle2
		uniaxialMaterial TzLiq1 $mTz $tzType $tultS $z50_tz $tRes $sEle1 $sEle2
		lappend liqSpringMatTags $mPy $mTz
		set latType "pyliq"
		set axType "tzliq"
		set useLiqS 1
		element zeroLength $eSpr $soilNd $dup -mat $mPy $mTz -dir 1 2 -doRayleigh 1
	} else {
		set tAx $tultS
		set z50Ax $z50_tz
		set pRes $pultS
		set tRes $tultS
		uniaxialMaterial PySimple1 $mPy $pyType $pultS $y50 $Cd_py $c_dash_py
		uniaxialMaterial TzSimple1 $mTz $tzType $tultS $z50_tz $c_dash_py
		set latType "py"
		set axType "tz"
		set useLiqS 0
		element zeroLength $eSpr $soilNd $dup -mat $mPy $mTz -dir 1 2 -doRayleigh 1
	}

	set xyS [nodeCoord $soilNd]
	set xs [lindex $xyS 0]
	set ys [lindex $xyS 1]
	lappend ssiSpringDump [list $eSpr $xs $ys $xP $y $xP $y "pile" $latType]
	lappend ssiSpringDump [list $eSpr $xs $ys $xP $y $xP $y "pile" $axType]
	lappend pileSprEleRec [list $eSpr $ip $iyP $isTip]
	lappend pileSpringPropsDump [list $ip $iyP $y [expr {-$y}] $nm $isTip $useLiqS \
		$pultS $pRes $tAx $tRes $y50 $z50Ax $latType $axType $trib 0.0]
	incr nSpr
}

set nPileSprings $nSpr

# ---- cap face springs (tributary height; py 1/2-per-side already in pOne) ----
foreach w $capFaceWork {
	lassign $w capNd xC yC soilNd dup pOne tOne
	springTie $capNd $dup $soilNd
	gravLockPush $soilNd $capNd
	set mPy [expr {$matTag_py_base + $nSpr}]
	set mTz [expr {$matTag_tz_base + $nSpr}]
	if {$pileSpring eq "elastic"} {
		set kP [expr {$pOne/$y50_cap}]
		set kT [expr {$tOne/$z50_cap}]
		uniaxialMaterial Elastic $mPy $kP
		uniaxialMaterial Elastic $mTz $kT
		set latType "py_elastic"
		set axType "tz_elastic"
	} else {
		uniaxialMaterial PySimple1 $mPy 2 $pOne $y50_cap $Cd_py $c_dash_py
		uniaxialMaterial TzSimple1 $mTz 1 $tOne $z50_cap $c_dash_py
		set latType "py"
		set axType "tz"
	}
	incr eSpr
	element zeroLength $eSpr $soilNd $dup -mat $mPy $mTz -dir 1 2 -doRayleigh 1
	lappend capFacePropsDump [list $eSpr $xC $yC $pOne $tOne $y50_cap $z50_cap]
	set xy [nodeCoord $soilNd]
	lappend ssiSpringDump [list $eSpr \
		[lindex $xy 0] [lindex $xy 1] $xC $yC $xC $yC "cap" $latType]
	lappend ssiSpringDump [list $eSpr \
		[lindex $xy 0] [lindex $xy 1] $xC $yC $xC $yC "cap" $axType]
	incr nSpr
}

# ---- cap soffit q-z ----
set nSoffit 0
foreach w $capSoffitWork {
	lassign $w capNd xS yS soilNd dup qOne
	springTie $capNd $dup $soilNd
	gravLockPush $soilNd $capNd
	set mQz [expr {$matTag_tz_base + $nSpr}]
	if {$pileSpring eq "elastic"} {
		set kQ [expr {$qOne/$z50_cap}]
		uniaxialMaterial Elastic $mQz $kQ
		set axType "qz_elastic"
	} else {
		uniaxialMaterial QzSimple1 $mQz 1 $qOne $z50_cap 0.0 $c_dash_py
		set axType "qz"
	}
	incr eSpr
	element zeroLength $eSpr $soilNd $dup -mat $mQz -dir 2 -doRayleigh 1
	lappend capSoffitPropsDump [list $eSpr $xS $yS $qOne $z50_cap]
	set xy [nodeCoord $soilNd]
	lappend ssiSpringDump [list $eSpr \
		[lindex $xy 0] [lindex $xy 1] $xS $yS $xS $yS "cap_soffit" $axType]
	incr nSpr
	incr nSoffit
}

set eleTag_spr_last $eSpr
set nSprings $nSpr
