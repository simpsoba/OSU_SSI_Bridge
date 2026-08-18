# BuildDeckElements.tcl
# Units: N, m, s
#
# Fixed PT39 frame. Call after BuildDeckNodes.tcl.
#
#   BarL 3009                                    BarR 3010
#      |                                            |
#   TL 3004 -- TLi 3005 -- TC 3006 -- TRi 3007 -- TR 3008     y = H_pier+dd
#      |          |          |          |          |
#              BL 3001 ---- BC 5 ---- BR 3003               y = H_pier
#
# Frame members: A_deck_frame, I_deck_frame.  geomTransf 3001.
# Center web BC-TC (pierTop_deckBC to 3006): A_link, I_link (pier A,I x eta_deckLink).
# Eles 3100-3110.

if {[info exists deckElementsDone] && $deckElementsDone} {
	return
}
if {![info exists A_link]} {
	error "BuildDeckElements.tcl: source BuildDeckNodes.tcl first"
}

# =====================================================================
# 4. ELEMENTS
# =====================================================================
geomTransf Linear 3001

# element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
# Soffit
element elasticBeamColumn 3100 \
	3001 $nodeTag_deck_BC \
	$A_deck_frame $Ec_deck $I_deck_frame 3001
element elasticBeamColumn 3101 \
	$nodeTag_deck_BC 3003 \
	$A_deck_frame $Ec_deck $I_deck_frame 3001
# Top flange
element elasticBeamColumn 3102 \
	3004 3005 \
	$A_deck_frame $Ec_deck $I_deck_frame 3001
element elasticBeamColumn 3103 \
	3005 3006 \
	$A_deck_frame $Ec_deck $I_deck_frame 3001
element elasticBeamColumn 3104 \
	3006 3007 \
	$A_deck_frame $Ec_deck $I_deck_frame 3001
element elasticBeamColumn 3105 \
	3007 3008 \
	$A_deck_frame $Ec_deck $I_deck_frame 3001
# Outer webs
element elasticBeamColumn 3106 \
	3001 3005 \
	$A_deck_frame $Ec_deck $I_deck_frame 3001
element elasticBeamColumn 3107 \
	3003 3007 \
	$A_deck_frame $Ec_deck $I_deck_frame 3001
# Center web (stiff)
element elasticBeamColumn 3108 \
	$nodeTag_deck_BC 3006 \
	$A_link $Ec_deck $I_link 3001
# Barriers
element elasticBeamColumn 3109 \
	3004 3009 \
	$A_deck_frame $Ec_deck $I_deck_frame 3001
element elasticBeamColumn 3110 \
	3008 3010 \
	$A_deck_frame $Ec_deck $I_deck_frame 3001

set eleTag_deck_last 3110
set deckElementsDone 1
