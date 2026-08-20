# analysis/RayleighDamping.tcl
# Goals: ξ = αM/(2ω) + β ω/2 at T1, T2. After ActivateEQBoundary, before EQ.
# OpenSees skips -rayleigh if all four coefs are 0, so "off" = offFac * (αM, β).
#
# Knobs may already be set by Overrides.tcl (TestMatrix). Defaults apply only
# when a variable is missing. Region *ON: 1 = full damping, 0 = near-zero.
# rayleighStiff: committed | initial  (which K gets β).
# Last region that owns a node overwrites αM. Keep pier / hinge last.

# ------------------------------------------------------------
# EDIT (defaults if not already set)
# ------------------------------------------------------------
if {![info exists rayleighT1]} {
	set rayleighT1 2.3;                   # <-- EDIT  s
}
if {![info exists rayleighT2]} {
	set rayleighT2 [expr {sqrt($cylinderSF)/20.0}];  # s  20 Hz lab -> proto
}
if {![info exists rayleighXi1]} {
	set rayleighXi1 0.03;                 # <-- EDIT  (-) at T1
}
if {![info exists rayleighXi2]} {
	set rayleighXi2 0.20;                 # <-- EDIT  (-) at T2
}
if {![info exists rayleighOffFac]} {
	set rayleighOffFac 1.0e-8;            # <-- EDIT  (-) near-zero scale
}
if {![info exists rayleighStiff]} {
	set rayleighStiff "committed";        # <-- EDIT  committed | initial
}
# Region full (1) vs off (0). Defaults match the old hard-coded layout.
if {![info exists rayleighSoilON]}       { set rayleighSoilON 0 }
if {![info exists rayleighBoundON]}      { set rayleighBoundON 0 }
if {![info exists rayleighSprPileON]}    { set rayleighSprPileON 0 }
if {![info exists rayleighSprCapFaceON]} { set rayleighSprCapFaceON 0 }
if {![info exists rayleighSprSoffitON]}  { set rayleighSprSoffitON 0 }
if {![info exists rayleighPilesON]}      { set rayleighPilesON 0 }
if {![info exists rayleighCapON]}        { set rayleighCapON 0 }
if {![info exists rayleighDeckON]}       { set rayleighDeckON 0 }
if {![info exists rayleighPierON]}       { set rayleighPierON 1 }
if {![info exists rayleighPierHingeON]}  { set rayleighPierHingeON 1 }
# ------------------------------------------------------------

if {$rayleighT1 <= 0.0 || $rayleighT2 <= 0.0 || $rayleighT1 == $rayleighT2} {
	error "RayleighDamping.tcl: rayleighT1 and rayleighT2 must be > 0 and distinct"
}
if {$rayleighStiff ne "committed" && $rayleighStiff ne "initial"} {
	error "RayleighDamping.tcl: rayleighStiff must be committed|initial (got '$rayleighStiff')"
}
foreach _onName {
	rayleighSoilON rayleighBoundON rayleighSprPileON rayleighSprCapFaceON
	rayleighSprSoffitON rayleighPilesON rayleighCapON rayleighDeckON
	rayleighPierON rayleighPierHingeON
} {
	set _on [set $_onName]
	if {$_on != 0 && $_on != 1} {
		error "RayleighDamping.tcl: $_onName must be 0 or 1 (got '$_on')"
	}
}
unset _onName _on

set w1 [expr {2.0*$pi/$rayleighT1}]
set w2 [expr {2.0*$pi/$rayleighT2}]
set w2w1den [expr {$w2*$w2 - $w1*$w1}]
set betaKcomm [expr {2.0*($rayleighXi2*$w2 - $rayleighXi1*$w1)/$w2w1den}]
set alphaM [expr {2.0*$rayleighXi1*$w1 - $betaKcomm*$w1*$w1}]
unset w2w1den
set aOff [expr {$rayleighOffFac*$alphaM}]
set bOff [expr {$rayleighOffFac*$betaKcomm}]

# Pick α and β for a region: on=1 full, on=0 near-zero.
# Returns list {aM bK bKinit bKcomm} for -rayleigh.
# Args: on (0|1) betaFull (the β to use when on; usually $betaKcomm)
proc rayleighFour {on betaFull} {
	global alphaM aOff bOff rayleighStiff
	if {$on} {
		set a $alphaM
		set b $betaFull
	} else {
		set a $aOff
		set b $bOff
	}
	if {$rayleighStiff eq "initial"} {
		return [list $a 0.0 $b 0.0]
	}
	return [list $a 0.0 0.0 $b]
}

# ξ(T) = αM/(2ω) + β ω/2, ω = 2π/T
# Args: T (s). Returns: xi (-)
proc rayleighXiAtT {T} {
	global alphaM betaKcomm pi
	set w [expr {2.0*$pi/$T}]
	return [expr {$alphaM/(2.0*$w) + 0.5*$betaKcomm*$w}]
}

puts [format "----- Rayleigh  T1=%.3g s xi=%.3g  T2=%.4g s xi=%.3g  stiff=%s  xi(0.3s)=%.3g -----" \
	$rayleighT1 $rayleighXi1 $rayleighT2 $rayleighXi2 $rayleighStiff [rayleighXiAtT 0.3]]

# region IDs (not knobs)
set reg_soil       801
set reg_bound      802
set reg_sprPile    803
set reg_sprCapFace 804
set reg_sprSoffit  805
set reg_piles      806
set reg_cap        807
set reg_deck       808
set reg_pier       809
set reg_pierHinge  810

# region $tag -eleRange $a $b -rayleigh $alphaM $betaK $betaKinit $betaKcomm

# soil: all continuum quads / SSPquads (near-field + FF, tags interleaved)
region $reg_soil -eleRange $eleTag_soil_base $eleTag_soil_last \
	-rayleigh {*}[rayleighFour $rayleighSoilON $betaKcomm]

# bound: Shin Lysmer dashpots or ASDEA ring
if {[info exists eleTag_bnd_last] && $eleTag_bnd_last >= $eleTag_bnd_base} {
	region $reg_bound -eleRange $eleTag_bnd_base $eleTag_bnd_last \
		-rayleigh {*}[rayleighFour $rayleighBoundON $betaKcomm]
}

# SSI springs: one contiguous ele range, pile | cap-face | soffit
if {$pileSpring ne "none" && [info exists nPileSprings] && $nPileSprings >= 1} {
	set eleSprPile0 $eleTag_spr_base
	set eleSprPile1 [expr {$eleSprPile0 + $nPileSprings - 1}]
	# sprPile: pile zeroLengths (p-y + t-z; tip q-z)
	region $reg_sprPile -eleRange $eleSprPile0 $eleSprPile1 \
		-rayleigh {*}[rayleighFour $rayleighSprPileON $betaKcomm]

	set nSoffitHere 0
	if {[info exists nSoffit]} { set nSoffitHere $nSoffit }
	set eleSprFace0 [expr {$eleSprPile1 + 1}]
	set eleSprFace1 [expr {$eleTag_spr_last - $nSoffitHere}]
	if {$eleSprFace1 >= $eleSprFace0} {
		# sprCapFace: six cap-face springs (p-y + t-z)
		region $reg_sprCapFace -eleRange $eleSprFace0 $eleSprFace1 \
			-rayleigh {*}[rayleighFour $rayleighSprCapFaceON $betaKcomm]
	}
	if {$nSoffitHere > 0} {
		# sprSoffit: cap soffit q-z
		region $reg_sprSoffit -eleRange [expr {$eleSprFace1 + 1}] $eleTag_spr_last \
			-rayleigh {*}[rayleighFour $rayleighSprSoffitON $betaKcomm]
	}
	unset eleSprPile0 eleSprPile1 nSoffitHere eleSprFace0 eleSprFace1
}

# piles: three shaft beams
if {[info exists eleTag_pile_last] && $eleTag_pile_last >= $eleTag_pile_base} {
	region $reg_piles -eleRange $eleTag_pile_base $eleTag_pile_last \
		-rayleigh {*}[rayleighFour $rayleighPilesON $betaKcomm]
}

# cap: pile-cap frame
region $reg_cap -eleRange $eleTag_cap_base $eleTag_cap_last \
	-rayleigh {*}[rayleighFour $rayleighCapON $betaKcomm]

# deck: box girder + barriers
region $reg_deck -eleRange $eleTag_deck_base $eleTag_deck_last \
	-rayleigh {*}[rayleighFour $rayleighDeckON $betaKcomm]

# pier: β/η (accounting for eta stiffness factor). FBC / elastic: full β.
set betaPier $betaKcomm
if {$pierEleType eq "lumpedPlasticity"} {
	if {![info exists eta_pier] || $eta_pier <= 0.0} {
		error "RayleighDamping.tcl: lumpedPlasticity needs eta_pier > 0"
	}
	set betaPier [expr {$betaKcomm/$eta_pier}]
}
region $reg_pier -ele $eleTag_pier \
	-rayleigh {*}[rayleighFour $rayleighPierON $betaPier]

# pierHinge: both ZLS. Last so αM wins at the hinge nodes.
if {$pierEleType eq "lumpedPlasticity"} {
	region $reg_pierHinge -ele $eleTag_pier_botSpr $eleTag_pier_topSpr \
		-rayleigh {*}[rayleighFour $rayleighPierHingeON $betaKcomm]
}

unset aOff bOff
