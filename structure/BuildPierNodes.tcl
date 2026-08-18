# BuildPierNodes.tcl
# Units: N, m, s
#
# Goals: pier nodes on +Y. Tags: Parameters.tcl
#   pierBase_capTC, pierBaseZeroLengthInner, pierTopZeroLengthInner, pierTop_deckBC
# Height: H_pier (Parameters.tcl).
#
# pierEleType (Parameters.tcl):
#   elasticBeamColumn -- A, E, I from PierSection
#   forceBeamColumn   -- ConcentratedCurvature (Fiber hinges + Elastic mid)
#   lumpedPlasticity  -- ZLS-I + eta*EI beam + ZLS-J (Ls = ratio * H; Lp plastic)
#
# Mass: nodal lump only (no element -mass).
#   Translation: m = rhoL*H/2 per end.
#   Rotation:    Irot = rhoL*H^3/105  (= consistent-mass diagonal rhoL*4L^2/420).
#   lumpedPlasticity: split m and Irot evenly on each ZLS pair (base 1-2, top 4-5).
#                equalDOF: base UX; top UX+UY (BuildPierElements.tcl).
#
# Expects: Parameters.tcl, model domain, PierSection.tcl already sourced.
# Does not apply boundary conditions (caller does).
# Beams/ZLS: BuildPierElements.tcl (after nodes; gravity waits until soil settle).

if {![info exists structNodeTags]} { set structNodeTags {} }

if {![info exists H_pier] || ![info exists pierEleType]} {
	error "BuildPierNodes.tcl: source Parameters.tcl first"
}
if {![info exists A_pier] || ![info exists Ec_pier] || ![info exists I_pier]} {
	error "BuildPierNodes.tcl: source structure/PierSection.tcl first"
}
if {![info exists rhoL_pier]} {
	error "BuildPierNodes.tcl: rhoL_pier missing (source Parameters.tcl)"
}
if {$pierEleType eq "forceBeamColumn"} {
	if {![info exists Lp_pier] || ![info exists secTag_elastic_pier] || ![info exists intTag_pier]} {
		error "BuildPierNodes.tcl: need Lp_pier, secTag_elastic_pier, intTag_pier from Parameters / PierSection"
	}
	if {[expr {2.0*$Lp_pier}] >= $H_pier} {
		error [format "BuildPier: 2*Lp=%.4f m >= H_pier=%.4f m" [expr {2.0*$Lp_pier}] $H_pier]
	}
}
if {![info exists nodeTag_pierBase_capTC] || ![info exists nodeTag_pierTop_deckBC]} {
	error "BuildPierNodes.tcl: source Parameters.tcl first (pier node tags)"
}
if {$pierEleType eq "lumpedPlasticity"} {
	if {![info exists Lp_pier] || ![info exists eta_pier] \
			|| ![info exists Ls_I_pier] || ![info exists Ls_J_pier]} {
		error "BuildPierNodes.tcl: need Lp, eta, Ls_I, Ls_J for lumpedPlasticity"
	}
	if {![info exists secTag_pier_I] || ![info exists secTag_pier_J]} {
		error "BuildPierNodes.tcl: need secTag_pier_I and secTag_pier_J"
	}
	if {![info exists eleTag_pier_botSpr] || ![info exists eleTag_pier_topSpr]} {
		error "BuildPierNodes.tcl: need eleTag_pier_botSpr and eleTag_pier_topSpr"
	}
	if {![info exists nodeTag_pierBaseZeroLengthInner] \
			|| ![info exists nodeTag_pierTopZeroLengthInner]} {
		error "BuildPierNodes.tcl: need both ZLS inner node tags"
	}
}

# Nodal mass: m = rhoL*H/2 per end; Irot = rhoL*H^3/105 (consistent-mass RZ diagonal)
set m_end_pier  [expr {0.5*$rhoL_pier*$H_pier}];              # kg
set Irot_pier   [expr {$rhoL_pier*pow($H_pier,3)/105.0}];      # kg*m^2

# =====================================================================
# 2. MODEL BUILDER / NODES
# =====================================================================
if {$pierEleType eq "elasticBeamColumn" || $pierEleType eq "forceBeamColumn"} {

	# node $tag $x $y -mass $mx $my $mRz
	# pierBase_capTC at y=0; pierTop_deckBC at y=H_pier
	node $nodeTag_pierBase_capTC 0.0 0.0 \
		-mass $m_end_pier $m_end_pier $Irot_pier
	node $nodeTag_pierTop_deckBC 0.0 $H_pier \
		-mass $m_end_pier $m_end_pier $Irot_pier
	lappend structNodeTags $nodeTag_pierBase_capTC $nodeTag_pierTop_deckBC

} elseif {$pierEleType eq "lumpedPlasticity"} {

	# Each ZLS pair splits that end's m and Irot evenly.
	set m_h   [expr {0.5*$m_end_pier}];   # kg
	set Irot_h [expr {0.5*$Irot_pier}];   # kg*m^2

	# node $tag $x $y -mass $mx $my $mRz
	node $nodeTag_pierBase_capTC 0.0 0.0 \
		-mass $m_h $m_h $Irot_h
	node $nodeTag_pierBaseZeroLengthInner 0.0 0.0 \
		-mass $m_h $m_h $Irot_h
	node $nodeTag_pierTopZeroLengthInner 0.0 $H_pier \
		-mass $m_h $m_h $Irot_h
	node $nodeTag_pierTop_deckBC 0.0 $H_pier \
		-mass $m_h $m_h $Irot_h
	lappend structNodeTags $nodeTag_pierBase_capTC \
		$nodeTag_pierBaseZeroLengthInner \
		$nodeTag_pierTopZeroLengthInner $nodeTag_pierTop_deckBC

} else {
	error "BuildPier: pierEleType must be elasticBeamColumn, forceBeamColumn, or lumpedPlasticity (got '$pierEleType')"
}
