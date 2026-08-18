# BuildPierElements.tcl
# Units: N, m, s
#
# Goals:
#   Build the pier element(s) on nodes that already exist.
#   Call after BuildPierNodes (SSI: after FoldStructNodes).
#
# pierEleType (Parameters.tcl):
#   elasticBeamColumn -- one beam, pierBase_capTC -> pierTop_deckBC
#   forceBeamColumn   -- same span; Fiber hinges + Elastic mid (Priestley Lp)
#   lumpedPlasticity  -- ZLS-I 1--2, eta*EI 2--4, ZLS-J 4--5

if {[info exists pierElementsDone] && $pierElementsDone} {
	return
}
if {![info exists pierEleType] || ![info exists nodeTag_pierBase_capTC]} {
	error "BuildPierElements.tcl: source Parameters.tcl first"
}

# =====================================================================
# 4. ELEMENTS
# =====================================================================
geomTransf $pierGeoTransf $transfTag_pier

if {$pierEleType eq "elasticBeamColumn"} {

	# element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
	element elasticBeamColumn $eleTag_pier \
		$nodeTag_pierBase_capTC $nodeTag_pierTop_deckBC \
		$A_pier $Ec_pier $I_pier $transfTag_pier

} elseif {$pierEleType eq "forceBeamColumn"} {

	# beamIntegration ConcentratedCurvature $tag $secI $LpI $secJ $LpJ $secE
	beamIntegration ConcentratedCurvature $intTag_pier \
		$secTag_pier $Lp_pier \
		$secTag_pier $Lp_pier \
		$secTag_elastic_pier

	# element forceBeamColumn $eleTag $iNode $jNode $transfTag $integrationTag
	element forceBeamColumn $eleTag_pier \
		$nodeTag_pierBase_capTC $nodeTag_pierTop_deckBC \
		$transfTag_pier $intTag_pier

} elseif {$pierEleType eq "lumpedPlasticity"} {

	# Fiber ZLS is P+M, no V. equalDOF $rNodeTag $cNodeTag $dof1 ...
	# base: UX only; UY free so gravity axial goes through ZLS-I (EA/Ls)
	equalDOF $nodeTag_pierBase_capTC $nodeTag_pierBaseZeroLengthInner 1
	# top: UX and UY; ZLS-J is rotation only, deck follows node 4 in translation
	equalDOF $nodeTag_pierTop_deckBC $nodeTag_pierTopZeroLengthInner 1 2
	# equalDOF $nodeTag_pierTop_deckBC $nodeTag_pierTopZeroLengthInner 1

	# element zeroLengthSection $eleTag $iNode $jNode $secTag <-orient $x1 $x2 $x3 $yp1 $yp2 $yp3>
	# Local x along pier (+Y); local y along -X -> fiber y in transverse plane
	element zeroLengthSection $eleTag_pier_botSpr \
		$nodeTag_pierBase_capTC $nodeTag_pierBaseZeroLengthInner $secTag_pier_I \
		-orient 0.0 1.0 0.0  -1.0 0.0 0.0 -doRayleigh 1

	# element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
	element elasticBeamColumn $eleTag_pier \
		$nodeTag_pierBaseZeroLengthInner $nodeTag_pierTopZeroLengthInner \
		$A_pier $Ec_pier [expr {$eta_pier*$I_pier}] $transfTag_pier

	element zeroLengthSection $eleTag_pier_topSpr \
		$nodeTag_pierTopZeroLengthInner $nodeTag_pierTop_deckBC $secTag_pier_J \
		-orient 0.0 1.0 0.0  -1.0 0.0 0.0 -doRayleigh 1

} else {
	error "BuildPierElements: unknown pierEleType '$pierEleType'"
}

set pierElementsDone 1
