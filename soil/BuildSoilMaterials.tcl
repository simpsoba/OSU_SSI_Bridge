# soil/BuildSoilMaterials.tcl
# Goals: PIMY clay / PDMY02+FSP sand, one mat per soil row.
# Knobs: Parameters.tcl. Layers / ramps: soil/Profiles.md.
#
# Sets lists used by BuildSoilMesh.tcl / DumpSoilProfile.tcl:
#   soilLayerNames         -- L2 L3 L5
#   soilYTop / soilYBot    -- unit bounds (m, z=0 at cap top)
#   soilYs / nSoilRows     -- vertical stations / row count (mesh reuses)
#   soilMatRow($iy)        -- quad material (FSP if sand)
#   soilSolidTags / soilFspTags -- updateMaterialStage lists
#   soilRho / soilC / ... keyed by row iy and by unit name
#
# =====================================================================
# 3. MATERIALS AND SECTIONS
# =====================================================================

if {![info exists soilProfile] || ![info exists matTag_soil_base] || ![info exists soilConstitutive]} {
	error "BuildSoilMaterials.tcl: source Parameters.tcl first"
}
if {![info exists H_cap] || ![info exists dy_soil] || ![info exists nSeg_below_tip]} {
	error "BuildSoilMaterials.tcl: need H_cap, dy_soil, nSeg_below_tip"
}
if {$soilConstitutive ne "elastic" && $soilConstitutive ne "inelastic"} {
	error "BuildSoilMaterials.tcl: soilConstitutive must be elastic|inelastic (got '$soilConstitutive')"
}

set nd_soil 2
set gam_max 0.1
set nYS 20
set pRef_c 1.0e5
set pRef_s 1.01e5
set d_clay 0.0
set d_sand 0.5

# PDMY02 trailing args (explicit so e can vary with Dr; OpenSees default e=0.6)
set pdmy_contr2 5.0
set pdmy_dilat2 3.0
set pdmy_liq1 1.0
set pdmy_liq2 0.0
set pdmy_cs1 0.9
set pdmy_cs2 0.02
set pdmy_cs3 0.7
set pdmy_pa 101.0

# ---- Kramer T5.1 / wiki knots (Profiles.md). Clay Br = 50 Gr. ----
set cLo_rho 1488.0; set cLo_Gr 5.728e7;  set cLo_c 3.59e4
set mLo_rho 1521.0; set mLo_Gr 6.639e7;  set mLo_c 3.97e4
set mHi_rho 1669.0; set mHi_Gr 1.10684e8; set mHi_c 5.84e4
set soft_rho 1300.0; set soft_Gr 1.30e7;  set soft_c 1.80e4
set stiff_rho 1800.0; set stiff_Gr 1.50e8; set stiff_c 7.50e4

# Wiki PDMY02 knots vs Dr: rho Gr Br phi PTA c1 c3 d1 d3 e k_pci
set pdmy50 [list 50.0 1900.0 1.00e8 2.33e8 33.5 25.5 0.045 0.15 0.06 0.15 0.70 $k_pci_L3a]
set pdmy60 [list 60.0 2000.0 1.10e8 2.40e8 35.0 26.0 0.028 0.05 0.10 0.05 0.65 $k_pci_L3b]
set pdmy75 [list 75.0 2100.0 1.30e8 2.60e8 36.5 26.0 0.013 0.0  0.30 0.0  0.55 $k_pci_L3c]

set L5s_rho 2260.0; set L5s_Gr 1.548e8; set L5s_Br 2.792e8; set L5s_phi 39.3
set L5s_PTA 26.0; set L5s_c1 0.0042; set L5s_c3 0.0; set L5s_d1 0.564; set L5s_d3 0.0
set L5s_Dr 95.0; set L5s_e 0.43

proc soilLerp {a b xi} {
	return [expr {$a + $xi*($b - $a)}]
}

proc soilXi {y yt yb} {
	set den [expr {$yt - $yb}]
	if {abs($den) < 1.0e-12} { return 0.0 }
	return [expr {($yt - $y)/$den}]
}

proc soilClampXi {xi} {
	if {$xi < 0.0} { return 0.0 }
	if {$xi > 1.0} { return 1.0 }
	return $xi
}

proc soilBlendList {A B t} {
	set out {}
	foreach a $A b $B {
		lappend out [expr {$a + $t*($b - $a)}]
	}
	return $out
}

# return: Dr rho Gr Br phi PTA c1 c3 d1 d3 e kpci
proc soilSandAtDr {Dr} {
	global pdmy50 pdmy60 pdmy75
	if {$Dr <= 50.0} { return $pdmy50 }
	if {$Dr >= 75.0} {
		set p $pdmy75
		lset p 0 $Dr
		return $p
	}
	if {$Dr <= 60.0} {
		set t [expr {($Dr - 50.0)/10.0}]
		set p [soilBlendList $pdmy50 $pdmy60 $t]
	} else {
		set t [expr {($Dr - 60.0)/15.0}]
		set p [soilBlendList $pdmy60 $pdmy75 $t]
	}
	lset p 0 $Dr
	return $p
}

proc soilClayPair {y yt yb rho0 Gr0 c0 rho1 Gr1 c1} {
	set xi [soilClampXi [soilXi $y $yt $yb]]
	set rho [soilLerp $rho0 $rho1 $xi]
	set Gr  [soilLerp $Gr0 $Gr1 $xi]
	set c   [soilLerp $c0 $c1 $xi]
	set Br  [expr {50.0*$Gr}]
	return [list $rho $Gr $Br $c]
}

proc soilGKToENu {G K} {
	set nu [expr {(3.0*$K - 2.0*$G)/(2.0*(3.0*$K + $G))}]
	if {$nu < 0.0} { set nu 0.0 }
	if {$nu > 0.49} { set nu 0.49 }
	set E [expr {2.0*$G*(1.0 + $nu)}]
	return [list $E $nu]
}

proc soilStoreENu {nm E nu} {
	global soilEIso soilNuIso
	set soilEIso($nm) $E
	set soilNuIso($nm) $nu
}

proc soilClearProps {} {
	global soilMatEle soilIsSand soilRho soilC soilPhi soilKpci soilDr soilG0
	global soilType soilGr soilBr soilPTA soilContr1 soilContr3 soilDilat1 soilDilat3
	global soilD soilFSP soilB_fsp soilGamMax soilPRef soilNYS
	global soilContr2 soilDilat2 soilLiq1 soilLiq2 soilE soilCs1 soilCs2 soilCs3 soilPa
	global soilEIso soilNuIso soilRowLayer soilMatRow soilSolidTag
	foreach a {soilMatEle soilIsSand soilRho soilC soilPhi soilKpci soilDr soilG0 \
		soilType soilGr soilBr soilPTA soilContr1 soilContr3 soilDilat1 soilDilat3 \
		soilD soilFSP soilB_fsp soilGamMax soilPRef soilNYS \
		soilContr2 soilDilat2 soilLiq1 soilLiq2 soilE soilCs1 soilCs2 soilCs3 soilPa \
		soilEIso soilNuIso soilRowLayer soilMatRow soilSolidTag} {
		array unset $a
	}
}

proc soilStoreClay {nm matTag rho Gr Br c} {
	global soilMatEle soilIsSand soilRho soilC soilPhi soilKpci soilDr soilG0
	global soilType soilGr soilBr soilPTA soilContr1 soilContr3 soilDilat1 soilDilat3
	global soilD soilFSP soilB_fsp soilGamMax soilPRef soilNYS
	global soilContr2 soilDilat2 soilLiq1 soilLiq2 soilE soilCs1 soilCs2 soilCs3 soilPa
	global d_clay gam_max pRef_c nYS
	set soilMatEle($nm) $matTag
	set soilIsSand($nm) 0
	set soilType($nm) "clay"
	set soilRho($nm) $rho
	set soilGr($nm) $Gr
	set soilBr($nm) $Br
	set soilC($nm) $c
	set soilPhi($nm) 0.0
	set soilKpci($nm) 0.0
	set soilDr($nm) 0.0
	set soilG0($nm) $Gr
	set soilPTA($nm) 0.0
	set soilContr1($nm) 0.0
	set soilContr3($nm) 0.0
	set soilDilat1($nm) 0.0
	set soilDilat3($nm) 0.0
	set soilD($nm) $d_clay
	set soilFSP($nm) 0
	set soilB_fsp($nm) 0.0
	set soilGamMax($nm) $gam_max
	set soilPRef($nm) $pRef_c
	set soilNYS($nm) $nYS
	set soilContr2($nm) 0.0
	set soilDilat2($nm) 0.0
	set soilLiq1($nm) 0.0
	set soilLiq2($nm) 0.0
	set soilE($nm) 0.0
	set soilCs1($nm) 0.0
	set soilCs2($nm) 0.0
	set soilCs3($nm) 0.0
	set soilPa($nm) 0.0
}

proc soilStoreSand {nm matTag rho Gr Br phi PTA c1 c3 d1 d3 Dr kpci e} {
	global soilMatEle soilIsSand soilRho soilC soilPhi soilKpci soilDr soilG0
	global soilType soilGr soilBr soilPTA soilContr1 soilContr3 soilDilat1 soilDilat3
	global soilD soilFSP soilB_fsp soilGamMax soilPRef soilNYS
	global soilContr2 soilDilat2 soilLiq1 soilLiq2 soilE soilCs1 soilCs2 soilCs3 soilPa
	global d_sand B_fsp gam_max pRef_s nYS
	global pdmy_contr2 pdmy_dilat2 pdmy_liq1 pdmy_liq2 pdmy_cs1 pdmy_cs2 pdmy_cs3 pdmy_pa
	set soilMatEle($nm) $matTag
	set soilIsSand($nm) 1
	set soilType($nm) "sand"
	set soilRho($nm) $rho
	set soilGr($nm) $Gr
	set soilBr($nm) $Br
	set soilC($nm) 0.0
	set soilPhi($nm) $phi
	set soilKpci($nm) $kpci
	set soilDr($nm) $Dr
	set soilG0($nm) $Gr
	set soilPTA($nm) $PTA
	set soilContr1($nm) $c1
	set soilContr3($nm) $c3
	set soilDilat1($nm) $d1
	set soilDilat3($nm) $d3
	set soilD($nm) $d_sand
	set soilFSP($nm) 1
	set soilB_fsp($nm) $B_fsp
	set soilGamMax($nm) $gam_max
	set soilPRef($nm) $pRef_s
	set soilNYS($nm) $nYS
	set soilContr2($nm) $pdmy_contr2
	set soilDilat2($nm) $pdmy_dilat2
	set soilLiq1($nm) $pdmy_liq1
	set soilLiq2($nm) $pdmy_liq2
	set soilE($nm) $e
	set soilCs1($nm) $pdmy_cs1
	set soilCs2($nm) $pdmy_cs2
	set soilCs3($nm) $pdmy_cs3
	set soilPa($nm) $pdmy_pa
}

# ---- unit geometry (m); z=0 at cap top ----
set soilLayerNames {L2 L3 L5}
set soilYTop(L2) 0.0
if {$soilProfile == 4} {
	set soilYBot(L2) [expr {-27.25*$foot}]
} else {
	set soilYBot(L2) [expr {-9.25*$foot}]
}
set soilYTop(L3) $soilYBot(L2)
set soilYBot(L3) [expr {-42.25*$foot}]
set soilYTop(L5) $soilYBot(L3)
set soilYBot(L5) [expr {-$H_cap - $L_pile - $nSeg_below_tip*$dy_soil}]

# ---- vertical stations (same list BuildSoilMesh.tcl uses) ----
set soilYs {}
lappend soilYs 0.0
lappend soilYs [expr {-0.5*$H_cap}]
lappend soilYs [expr {-$H_cap}]
set yBot $soilYBot(L5)
set yCoord [expr {-$H_cap}]
while {$yCoord > $yBot + 1.0e-9} {
	set yCoord [expr {$yCoord - $dy_soil}]
	if {$yCoord < $yBot} { set yCoord $yBot }
	lappend soilYs $yCoord
}
set soilYs [lsort -real -decreasing -unique $soilYs]
set nY [llength $soilYs]
set nSoilRows [expr {$nY - 1}]

proc soilLayerAtY {yc} {
	global soilLayerNames soilYTop soilYBot
	foreach nm $soilLayerNames {
		if {$yc <= $soilYTop($nm) + 1.0e-6 && $yc >= $soilYBot($nm) - 1.0e-6} {
			return $nm
		}
	}
	return [lindex $soilLayerNames end]
}

proc soilRowAtY {y} {
	global soilYs nSoilRows
	set best 0
	set bd 1.0e99
	for {set iy 0} {$iy < $nSoilRows} {incr iy} {
		set yt [lindex $soilYs $iy]
		set yb [lindex $soilYs [expr {$iy + 1}]]
		set yc [expr {0.5*($yt + $yb)}]
		set d [expr {abs($yc - $y)}]
		if {$d < $bd} {
			set bd $d
			set best $iy
		}
	}
	return $best
}

soilClearProps

set soilIsSand(L2) 0
set soilType(L2) "clay"
if {$soilProfile == 1} {
	set soilIsSand(L3) 1
	set soilType(L3) "sand"
	set soilIsSand(L5) 1
	set soilType(L5) "sand"
} elseif {$soilProfile == 2} {
	set soilIsSand(L3) 0
	set soilType(L3) "clay"
	set soilIsSand(L5) 1
	set soilType(L5) "sand"
} elseif {$soilProfile == 3 || $soilProfile == 4} {
	set soilIsSand(L3) 0
	set soilType(L3) "clay"
	set soilIsSand(L5) 0
	set soilType(L5) "clay"
} else {
	error "BuildSoilMaterials.tcl: soilProfile must be 1, 2, 3, or 4 (got $soilProfile)"
}

# First-row mat of each unit (DumpPileSprings / sketches)
array unset unitFirst
set soilSolidTags {}
set soilFspTags {}
array unset soilRhoByMat

for {set iy 0} {$iy < $nSoilRows} {incr iy} {
	set yT [lindex $soilYs $iy]
	set yB [lindex $soilYs [expr {$iy + 1}]]
	set yc [expr {0.5*($yT + $yB)}]
	set nm [soilLayerAtY $yc]
	set soilRowLayer($iy) $nm
	set solidTag [expr {$matTag_soil_base + 1 + $iy}]
	set fspTag   [expr {$matTag_soil_base + 100 + $iy}]
	set isSand $soilIsSand($nm)

	if {$nm eq "L2"} {
		if {$soilProfile == 1} {
			lassign [soilClayPair $yc $soilYTop(L2) $soilYBot(L2) \
				$cLo_rho $cLo_Gr $cLo_c $mHi_rho $mHi_Gr $mHi_c] \
				rho Gr Br c
		} elseif {$soilProfile == 4} {
			lassign [soilClayPair $yc $soilYTop(L2) $soilYBot(L2) \
				$soft_rho $soft_Gr $soft_c $mLo_rho $mLo_Gr $mLo_c] \
				rho Gr Br c
		} else {
			lassign [soilClayPair $yc $soilYTop(L2) $soilYBot(L2) \
				$cLo_rho $cLo_Gr $cLo_c $mLo_rho $mLo_Gr $mLo_c] \
				rho Gr Br c
		}
	} elseif {$nm eq "L3"} {
		if {$soilProfile == 1} {
			set isSand 1
			set xi [soilClampXi [soilXi $yc $soilYTop(L3) $soilYBot(L3)]]
			set Dr [soilLerp 50.0 75.0 $xi]
			lassign [soilSandAtDr $Dr] DrNow rho Gr Br phi PTA c1 c3 d1 d3 e kpci
		} else {
			lassign [soilClayPair $yc $soilYTop(L3) $soilYBot(L3) \
				$mLo_rho $mLo_Gr $mLo_c $mHi_rho $mHi_Gr $mHi_c] \
				rho Gr Br c
		}
	} else {
		# L5
		if {$soilProfile == 1 || $soilProfile == 2} {
			set isSand 1
			set rho $L5s_rho; set Gr $L5s_Gr; set Br $L5s_Br
			set phi $L5s_phi; set PTA $L5s_PTA
			set c1 $L5s_c1; set c3 $L5s_c3; set d1 $L5s_d1; set d3 $L5s_d3
			set Dr $L5s_Dr; set e $L5s_e; set kpci $k_pci_L5
		} else {
			lassign [soilClayPair $yc $soilYTop(L5) $soilYBot(L5) \
				$mHi_rho $mHi_Gr $mHi_c $stiff_rho $stiff_Gr $stiff_c] \
				rho Gr Br c
		}
	}

	if {$isSand} {
		if {$soilConstitutive eq "elastic"} {
			lassign [soilGKToENu $Gr $Br] Emod nuVal
			nDMaterial ElasticIsotropic3D $solidTag $Emod $nuVal $rho
			soilStoreENu $iy $Emod $nuVal
		} else {
			nDMaterial PressureDependMultiYield02 $solidTag $nd_soil \
				$rho $Gr $Br $phi $gam_max $pRef_s $d_sand \
				$PTA $c1 $c3 $d1 $d3 $nYS \
				$pdmy_contr2 $pdmy_dilat2 $pdmy_liq1 $pdmy_liq2 $e
		}
		nDMaterial FluidSolidPorous $fspTag $nd_soil $solidTag $B_fsp
		soilStoreSand $iy $fspTag $rho $Gr $Br $phi $PTA $c1 $c3 $d1 $d3 $Dr $kpci $e
		set soilMatRow($iy) $fspTag
		set soilRhoByMat($fspTag) $rho
		lappend soilFspTags $fspTag
	} else {
		if {$soilConstitutive eq "elastic"} {
			lassign [soilGKToENu $Gr $Br] Emod nuVal
			nDMaterial ElasticIsotropic3D $solidTag $Emod $nuVal $rho
			soilStoreENu $iy $Emod $nuVal
		} else {
			nDMaterial PressureIndependMultiYield $solidTag $nd_soil \
				$rho $Gr $Br $c $gam_max 0.0 $pRef_c $d_clay $nYS
		}
		soilStoreClay $iy $solidTag $rho $Gr $Br $c
		set soilMatRow($iy) $solidTag
		set soilRhoByMat($solidTag) $rho
	}
	set soilSolidTag($iy) $solidTag
	lappend soilSolidTags $solidTag
	if {![info exists unitFirst($nm)]} {
		set unitFirst($nm) $iy
		set soilMatEle($nm) $soilMatRow($iy)
	}
}

# Unit-level props = first row of that unit (DumpPileSprings y-bands).
# Cap Mokwa uses L2 at mid-cap, set below.
foreach nm $soilLayerNames {
	if {![info exists unitFirst($nm)]} { continue }
	set iy0 $unitFirst($nm)
	set soilRho($nm) $soilRho($iy0)
	set soilC($nm) $soilC($iy0)
	set soilGr($nm) $soilGr($iy0)
	set soilBr($nm) $soilBr($iy0)
	set soilPhi($nm) $soilPhi($iy0)
	set soilDr($nm) $soilDr($iy0)
	set soilKpci($nm) $soilKpci($iy0)
	set soilG0($nm) $soilG0($iy0)
	set soilFSP($nm) $soilFSP($iy0)
	set soilB_fsp($nm) $soilB_fsp($iy0)
	set soilGamMax($nm) $soilGamMax($iy0)
}

# Cap springs: L2 at mid-cap (not the first thin row if it is tiny)
set iyCap [soilRowAtY [expr {-0.5*$H_cap}]]
set soilC(L2) $soilC($iyCap)
set soilRho(L2) $soilRho($iyCap)
