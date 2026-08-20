# soil/BuildSoilMesh.tcl
# Goals: near-field quads + Shin FF or ASDEA ring (stage 0).
# Knobs: Parameters.tcl. Tags: TAGS CONVENTION.
# Call after structure nodes (ndf=3), with: model BasicBuilder -ndm 2 -ndf 2
# Layout: soil/Profiles.md. Pile shafts and cap: continuum through, no excavation.
#
# =====================================================================
# 2. MODEL BUILDER / NODES
# =====================================================================

if {![info exists soilLayerNames] || ![info exists t_soil]} {
	error "BuildSoilMesh.tcl: source BuildSoilMaterials.tcl first"
}
if {![info exists H_cap] || ![info exists s_pile_cap]} {
	error "BuildSoilMesh.tcl: source structure builders first"
}
if {$soilBoundary ne "Shin" && $soilBoundary ne "ASDEA"} {
	error "BuildSoilMesh.tcl: soilBoundary must be Shin or ASDEA (got $soilBoundary)"
}
if {$soilEleType ne "quad" && $soilEleType ne "SSPquad"} {
	error "BuildSoilMesh.tcl: soilEleType must be quad|SSPquad (got '$soilEleType')"
}
# Continuum half-width: NF only (ASDEA) or NF+FF (Shin)
if {$soilBoundary eq "Shin"} {
	set xMeshHalf [expr {$L_half + $w_FF}]
} else {
	set xMeshHalf $L_half
}
set meshDir [file dirname [file normalize [info script]]]

# ---- vertical stations (m) ----
# Built in BuildSoilMaterials.tcl (one mat per row). Do not rebuild.
if {![info exists soilYs] || [llength $soilYs] < 2} {
	error "BuildSoilMesh.tcl: soilYs from BuildSoilMaterials.tcl required"
}
if {![info exists nSoilRows] || ![info exists soilMatRow]} {
	error "BuildSoilMesh.tcl: soilMatRow / nSoilRows from BuildSoilMaterials.tcl required"
}
set nY [llength $soilYs]
if {$nY != $nSoilRows + 1} {
	error [format "BuildSoilMesh.tcl: nY=%d but nSoilRows=%d" $nY $nSoilRows]
}

# ---- horizontal stations (m), symmetric stepped bands ----
# soilDxBands: near field {mesh size, x end}. Last x end = L_half.
# Shin: one extra station at L_half+w_FF. ASDEA: continuum stops at L_half.
if {![info exists soilDxBands] || [llength $soilDxBands] < 1} {
	error "BuildSoilMesh.tcl: set soilDxBands in Parameters.tcl"
}
set xTol 1.0e-4
set ib 0
foreach row $soilDxBands {
	incr ib
	if {[llength $row] != 2} {
		error [format "BuildSoilMesh.tcl: soilDxBands row %d needs {mesh size  x end}" $ib]
	}
	set dxRow [lindex $row 0]
	set xEndRow [lindex $row 1]
	if {$dxRow <= 0.0} {
		error [format "BuildSoilMesh.tcl: soilDxBands row %d mesh size must be > 0" $ib]
	}
	if {$ib > 1 && $xEndRow <= $xEndPrev + $xTol} {
		error [format "BuildSoilMesh.tcl: soilDxBands row %d x end must increase" $ib]
	}
	set xEndPrev $xEndRow
}
set nBand [llength $soilDxBands]
set lastEnd [lindex [lindex $soilDxBands end] 1]
if {![info exists L_half]} {
	set L_half $lastEnd
}
if {abs($lastEnd - $L_half) > $xTol} {
	error [format "BuildSoilMesh.tcl: last soilDxBands x end (%.4g m) must equal L_half (%.4g m)" \
		$lastEnd $L_half]
}
if {![info exists w_FF] || $w_FF <= 0.0} {
	error "BuildSoilMesh.tcl: w_FF must be > 0 (Shin column / ASDEA ring)"
}

set soilXs {}
proc pushX {lst x} {
	upvar 1 $lst L
	set tol 1.0e-6
	foreach v $L {
		if {abs($v - $x) < $tol} { return }
	}
	lappend L $x
}
proc fillBand {lst x0 x1 dx} {
	# inclusive endpoints; step ~ dx (shrink last step to land on x1)
	upvar 1 $lst L
	pushX L $x0
	if {$dx <= 0 || $x1 <= $x0 + 1.0e-9} {
		pushX L $x1
		return
	}
	set x $x0
	while {$x + $dx < $x1 - 1.0e-7} {
		set x [expr {$x + $dx}]
		pushX L $x
	}
	pushX L $x1
}

# Must-hit: pile axes (inner band also lands on 0, +/-s_pile_cap)
foreach xP [list [expr {-$s_pile_cap}] 0.0 $s_pile_cap] {
	pushX soilXs $xP
}

# Positive half: NF bands (0 -> L_half), then Shin FF station at L_half+w_FF
set soilMeshLeftover {}
set xPrev 0.0
for {set ib 0} {$ib < $nBand} {incr ib} {
	set dx [lindex [lindex $soilDxBands $ib] 0]
	set xEnd [lindex [lindex $soilDxBands $ib] 1]
	if {$xEnd > $L_half} { set xEnd $L_half }
	if {$xEnd <= $xPrev + 1.0e-9} { continue }
	set span [expr {$xEnd - $xPrev}]
	set nFull [expr {int(floor(($span + 1.0e-9)/$dx))}]
	set rem [expr {$span - $nFull*$dx}]
	if {$rem > $xTol && $rem < $dx - $xTol} {
		lappend soilMeshLeftover [format "%.3g ft at %.3g-%.3g ft" \
			[expr {$rem/$foot}] \
			[expr {($xEnd - $rem)/$foot}] \
			[expr {$xEnd/$foot}]]
	}
	fillBand soilXs $xPrev $xEnd $dx
	set xPrev $xEnd
}
if {[llength $soilMeshLeftover] > 0} {
	puts [format "----- Soil mesh leftover  %s -----" [join $soilMeshLeftover "; "]]
}
pushX soilXs $L_half
if {$soilBoundary eq "Shin"} {
	pushX soilXs [expr {$L_half + $w_FF}]
}

# Mirror to negative half
set posList $soilXs
foreach x $posList {
	if {$x > 1.0e-9} {
		pushX soilXs [expr {-$x}]
	}
}
set soilXs [lsort -real $soilXs]
set nX [llength $soilXs]

# ---- nodes: tag = nodeTag_soil_base + ix*nY + iy (stride = nY, no wrap) ----
# Spring / ASDEA bases move later via ensureAbove max(getNodeTags|getEleTags).
set soilNodeStride $nY
set soilNodeLast [soilNodeTag [expr {$nX - 1}] [expr {$nY - 1}]]
puts [format "----- Soil tags  nX=%d nY=%d stride=%d  nodes %d..%d -----" \
	$nX $nY $soilNodeStride $nodeTag_soil_base $soilNodeLast]

# node $nodeTag $xCoord $yCoord
for {set ix 0} {$ix < $nX} {incr ix} {
	set x [lindex $soilXs $ix]
	for {set iy 0} {$iy < $nY} {incr iy} {
		set y [lindex $soilYs $iy]
		node [soilNodeTag $ix $iy] $x $y
	}
}

# =====================================================================
# 4. ELEMENTS
# =====================================================================
# ---- quads ----
# Syntax (node order CCW from top-left: n1 TL, n4 BL, n3 BR, n2 TR):
#   element quad    $eleTag $n1 $n4 $n3 $n2 $thick $type $matTag <$pressure $rho $b1 $b2>
#   element SSPquad $eleTag $n1 $n4 $n3 $n2 $matTag $type $thick <$b1 $b2>
# quad: eleDensity = soilRhoByMat($mat) (explicit; avoid 0 -> getRho fallback).
# SSPquad: no rho argument -- mass from nDMaterial::getRho() (PIMY/PDMY/ElasticIsotropic3D
#   already store rho; FSP sands use the wrapped solid's rho via getRho when available).
# b1, b2 = body forces (horizontal 0, vertical self-weight)
set e $eleTag_soil_base
set soilEleTags {}
array unset soilEleLayer
array set soilEleLayer {}
array unset soilEleThick
array set soilEleThick {}
set nQuad 0
set pressure 0.0
set bfx 0.0

for {set ix 0} {$ix < $nX - 1} {incr ix} {
	set xL [lindex $soilXs $ix]
	set xR [lindex $soilXs [expr {$ix + 1}]]
	set xc [expr {0.5*($xL + $xR)}]
	# thickness: Shin FF columns (beyond L_half) use t_FF; else t_soil
	set inFF [expr {$soilBoundary eq "Shin" && \
		abs($xc) > $L_half - 1.0e-6}]
	if {$inFF} {
		set thick $t_FF
	} else {
		set thick $t_soil
	}
	for {set iy 0} {$iy < $nY - 1} {incr iy} {
		set yT [lindex $soilYs $iy]
		set yB [lindex $soilYs [expr {$iy + 1}]]
		set yc [expr {0.5*($yT + $yB)}]
		set n1 [soilNodeTag $ix $iy]
		set n2 [soilNodeTag [expr {$ix + 1}] $iy]
		set n3 [soilNodeTag [expr {$ix + 1}] [expr {$iy + 1}]]
		set n4 [soilNodeTag $ix [expr {$iy + 1}]]
		# all four nodes must exist
		set ok 1
		foreach nn [list $n1 $n2 $n3 $n4] {
			if {[lsearch -exact [getNodeTags] $nn] < 0} { set ok 0; break }
		}
		if {!$ok} { continue }

		set nm $soilRowLayer($iy)
		set mat $soilMatRow($iy)
		if {![info exists soilRhoByMat($mat)]} {
			error "BuildSoilMesh.tcl: no soilRhoByMat for mat=$mat (row $iy $nm)"
		}
		set rho $soilRhoByMat($mat)
		if {$soilIsSand($iy)} {
			set bfy [expr {-($rho - $rho_w)*$gravity_accel}]
		} else {
			set bfy [expr {-$rho*$gravity_accel}]
		}
		incr e
		if {$soilEleType eq "SSPquad"} {
			# element SSPquad $tag $i $j $k $l $matTag $type $thick <$b1 $b2>
			element SSPquad $e $n1 $n4 $n3 $n2 $mat "PlaneStrain" $thick \
				$bfx $bfy
		} else {
			# element quad $tag $i $j $k $l $thick $type $matTag <$p $rho $b1 $b2>
			element quad $e $n1 $n4 $n3 $n2 $thick PlaneStrain $mat \
				$pressure $rho $bfx $bfy
		}
		lappend soilEleTags $e
		set soilEleLayer($e) $nm
		set soilEleThick($e) $thick
		incr nQuad
	}
}
set eleTag_soil_last $e

# Node -> continuum elements (for SSI springs / partition -samePart)
array unset soilNodeEles
array set soilNodeEles {}
foreach ee $soilEleTags {
	foreach nn [eleNodes $ee] {
		lappend soilNodeEles($nn) $ee
	}
}

# =====================================================================
# 5. BOUNDARY CONDITIONS / CONSTRAINTS
# =====================================================================
# ---- base fixity (Shin only; ASDEA Stage 0 provides restraints) ----
# fix $nodeTag $dx $dy  (1 = fixed)
set iyBot [expr {$nY - 1}]
if {$soilBoundary eq "Shin"} {
	for {set ix 0} {$ix < $nX} {incr ix} {
		set nTag [soilNodeTag $ix $iyBot]
		if {[lsearch -exact [getNodeTags] $nTag] >= 0} {
			fix $nTag 1 1
		}
	}
}

# ---- Shin free-field equalDOF on outer column faces ----
# equalDOF $rNodeTag $cNodeTag $dof1 ... -- retained = outermost (Lysmer node).
# Base iy: UX only (UY remains on fix SP -> no SP+MP clash under Plain).
# Above base: UX and UY (1D FF column).
if {$soilBoundary eq "Shin"} {
	set ixOuterL 0
	set ixInnerL 1
	set ixInnerR [expr {$nX - 2}]
	set ixOuterR [expr {$nX - 1}]
	set iyBot [expr {$nY - 1}]
	foreach side [list \
		[list $ixOuterL $ixInnerL] \
		[list $ixOuterR $ixInnerR]] {
		lassign $side ixOuter ixInner
		for {set iy 0} {$iy < $nY} {incr iy} {
			set nOuter [soilNodeTag $ixOuter $iy]
			set nInner [soilNodeTag $ixInner $iy]
			if {[lsearch -exact [getNodeTags] $nOuter] < 0} { continue }
			if {[lsearch -exact [getNodeTags] $nInner] < 0} { continue }
			if {$iy == $iyBot} {
				equalDOF $nOuter $nInner 1
			} else {
				equalDOF $nOuter $nInner 1 2
			}
		}
	}
}

# Export mesh meta for springs / sketch / EQ boundary
set soil_nX $nX
set soil_nY $nY
set soil_nQuad $nQuad

# ASDEA ring (Stage 0). Velocity series required for bottom -fx.
if {$soilBoundary eq "ASDEA"} {
	source [file join [file dirname $meshDir] analysis BuildVelSeries.tcl]
	source [file join $meshDir BuildASDEABoundary.tcl]
}
