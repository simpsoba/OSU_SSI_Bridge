# BuildDeckNodes.tcl
# Units: N, m, s
#
# Goals: PT39 box as a stiff elasticBeamColumn frame (2D transverse).
# Center web uses pier A,I x eta_deckLink; other members use A_deck_frame,
# I_deck_frame (eta on geometry only; E = Ec_deck).
#
#   BarL 3009                                    BarR 3010
#      |                                            |
#   TL 3004 -- TLi 3005 -- TC 3006 -- TRi 3007 -- TR 3008     y = H_pier+dd
#      |          |          |          |          |
#              BL 3001 ---- BC 5 ---- BR 3003               y = H_pier
#
# BC = pier top (same node, deck mass stacked with IncrMass, no equalDOF).
# TLi/TRi = inner web (overhang).  BarL/BarR = jersey top.
#
# Mass: m = dens*A_deck*L_trib; length-weighted half-member lump -> scale to m;
# Irot: pier-style consistent-mass diagonal per member, I_end = m_mem L^2/105.
#
# Expects: Parameters.tcl, PierSection.tcl (A_pier, I_pier), BuildPierNodes.tcl.
# Frame beams: BuildDeckElements.tcl after nodes.

if {![info exists A_deck] || ![info exists m_deck] || ![info exists H_pier]} {
	error "BuildDeckNodes.tcl: source Parameters.tcl first"
}
if {![info exists nodeTag_deck_BC]} {
	error "BuildDeckNodes.tcl: need nodeTag_deck_BC from Parameters.tcl"
}
if {![info exists A_pier] || ![info exists I_pier]} {
	error "BuildDeckNodes.tcl: source PierSection.tcl first (need A_pier, I_pier)"
}
if {![info exists structNodeTags]} { set structNodeTags {} }

set scriptDir [file dirname [file normalize [info script]]]
source [file join $scriptDir IncrMass.tcl]

set y0 $H_pier
set xTopOuter [expr {0.5*$dw_deck}]
set xOverIn   [expr {0.5*$dw_deck - $cw_deck}]
set xSoffit   [expr {0.5*$sw_deck}]
set yTop      [expr {$y0 + $dd_deck}]
set yBar      [expr {$yTop + $bh_deck}]
set yCG       [expr {$y0 + $yb_deck}]

# Node names -> tags (BC shares pierTop_deckBC)
set deckNodeTag(BL)   3001
set deckNodeTag(BC)   $nodeTag_deck_BC
set deckNodeTag(BR)   3003
set deckNodeTag(TL)   3004
set deckNodeTag(TLi)  3005
set deckNodeTag(TC)   3006
set deckNodeTag(TRi)  3007
set deckNodeTag(TR)   3008
set deckNodeTag(BarL) 3009
set deckNodeTag(BarR) 3010

set deckNodeXY(BL)   [list [expr {-$xSoffit}] $y0]
set deckNodeXY(BC)   [list 0.0 $y0]
set deckNodeXY(BR)   [list $xSoffit $y0]
set deckNodeXY(TL)   [list [expr {-$xTopOuter}] $yTop]
set deckNodeXY(TLi)  [list [expr {-$xOverIn}] $yTop]
set deckNodeXY(TC)   [list 0.0 $yTop]
set deckNodeXY(TRi)  [list $xOverIn $yTop]
set deckNodeXY(TR)   [list $xTopOuter $yTop]
set deckNodeXY(BarL) [list [expr {-$xTopOuter}] $yBar]
set deckNodeXY(BarR) [list $xTopOuter $yBar]

# Members: {nameI nameJ kind}  kind=frame|center
set deckMembers [list \
	[list BL BC frame] [list BC BR frame] \
	[list TL TLi frame] [list TLi TC frame] [list TC TRi frame] [list TRi TR frame] \
	[list BL TLi frame] [list BR TRi frame] \
	[list BC TC center] \
	[list TL BarL frame] [list TR BarR frame] \
]

# Length-weighted half-member lump
set deckNames [list BL BC BR TL TLi TC TRi TR BarL BarR]
foreach nm $deckNames { set deckW($nm) 0.0 }

# Distance between two deck {x y} lists (m).
proc deckDist {xy1 xy2} {
	lassign $xy1 x1 y1
	lassign $xy2 x2 y2
	return [expr {hypot($x2 - $x1, $y2 - $y1)}]
}

foreach mem $deckMembers {
	lassign $mem a b kind
	set L [deckDist $deckNodeXY($a) $deckNodeXY($b)]
	set deckW($a) [expr {$deckW($a) + 0.5*$L}]
	set deckW($b) [expr {$deckW($b) + 0.5*$L}]
}

set Wsum 0.0
foreach nm $deckNames { set Wsum [expr {$Wsum + $deckW($nm)}] }
if {$Wsum <= 0.0} {
	error "BuildDeckNodes.tcl: zero member-length weight sum"
}
foreach nm $deckNames {
	set deckM($nm) [expr {$m_deck*$deckW($nm)/$Wsum}]
	set deckIrot($nm) 0.0
}

# Member lengths: translational weight uses half-L; rotary uses pier-style
# consistent-mass diagonal I_end = (m_mem) L^2/105 with m_mem = m_deck*L/SumL.
set Lsum $Wsum
foreach mem $deckMembers {
	lassign $mem a b kind
	set L [deckDist $deckNodeXY($a) $deckNodeXY($b)]
	set m_mem [expr {$m_deck*$L/$Lsum}]
	set Irot_end [expr {$m_mem*$L*$L/105.0}]
	set deckIrot($a) [expr {$deckIrot($a) + $Irot_end}]
	set deckIrot($b) [expr {$deckIrot($b) + $Irot_end}]
}

# Diagnostics: solid-block Iz vs Steiner (not used for assignment)
set Iz_deck [expr {$m_deck*($dw_deck*$dw_deck + $dd_deck*$dd_deck)/12.0}]
set I_steiner_deck 0.0
foreach nm $deckNames {
	lassign $deckNodeXY($nm) x y
	set dx $x
	set dy [expr {$y - $yCG}]
	set I_steiner_deck [expr {$I_steiner_deck + $deckM($nm)*($dx*$dx + $dy*$dy)}]
}
set nDeck [llength $deckNames]
set Irot_deck_sum 0.0
foreach nm $deckNames {
	set Irot_deck_sum [expr {$Irot_deck_sum + $deckIrot($nm)}]
}

# =====================================================================
# 2. MODEL BUILDER / NODES
# =====================================================================
# node $tag $x $y -mass $mx $my $mRz
# BC = nodeTag_deck_BC (pierTop_deckBC): IncrMass below
# --- soffit ---
node 3001 [expr {-$xSoffit}] $y0 \
	-mass $deckM(BL) $deckM(BL) $deckIrot(BL)
node 3003 $xSoffit $y0 \
	-mass $deckM(BR) $deckM(BR) $deckIrot(BR)
# --- top flange ---
node 3004 [expr {-$xTopOuter}] $yTop \
	-mass $deckM(TL) $deckM(TL) $deckIrot(TL)
node 3005 [expr {-$xOverIn}] $yTop \
	-mass $deckM(TLi) $deckM(TLi) $deckIrot(TLi)
node 3006 0.0 $yTop \
	-mass $deckM(TC) $deckM(TC) $deckIrot(TC)
node 3007 $xOverIn $yTop \
	-mass $deckM(TRi) $deckM(TRi) $deckIrot(TRi)
node 3008 $xTopOuter $yTop \
	-mass $deckM(TR) $deckM(TR) $deckIrot(TR)
# --- barriers ---
node 3009 [expr {-$xTopOuter}] $yBar \
	-mass $deckM(BarL) $deckM(BarL) $deckIrot(BarL)
node 3010 $xTopOuter $yBar \
	-mass $deckM(BarR) $deckM(BarR) $deckIrot(BarR)
lappend structNodeTags 3001 3003 3004 3005 3006 3007 3008 3009 3010

# Stack deck BC mass onto pier top
set mPierTopBefore [nodeMass $nodeTag_deck_BC 1]
IncrMass $nodeTag_deck_BC $deckM(BC) $deckM(BC) $deckIrot(BC)
set mPierTopAfter [nodeMass $nodeTag_deck_BC 1]

set A_link [expr {$eta_deckLink_A*$A_pier}]
set I_link [expr {$eta_deckLink_I*$I_pier}]

# Dump for PlotDeckFrameConcept.py
set plotRoot [file join [file dirname [file dirname [file normalize [info script]]]] plot]
file mkdir [file join $plotRoot out deck]
set deckFrameOutPath [file join $plotRoot out deck deck_frame.json]
set jsonFd [open $deckFrameOutPath w]
puts $jsonFd "\{"
puts $jsonFd "  \"units\": \"N, m, s\","
puts $jsonFd [format "  \"m_deck\": %.8g," $m_deck]
puts $jsonFd [format "  \"Iz_deck\": %.8g," $Iz_deck]
puts $jsonFd [format "  \"I_steiner\": %.8g," $I_steiner_deck]
puts $jsonFd [format "  \"Irot_fill\": %.8g," 0.0]
puts $jsonFd [format "  \"Irot_sum\": %.8g," $Irot_deck_sum]
puts $jsonFd "  \"Irot_method\": \"member consistent-mass diagonal m_mem L^2/105\","
puts $jsonFd [format "  \"y0\": %.8g," $y0]
puts $jsonFd [format "  \"yCG\": %.8g," $yCG]
puts $jsonFd [format "  \"dw\": %.8g, \"dd\": %.8g, \"sw\": %.8g, \"cw\": %.8g," \
	$dw_deck $dd_deck $sw_deck $cw_deck]
puts $jsonFd [format "  \"td\": %.8g, \"ts\": %.8g, \"tw\": %.8g, \"bh\": %.8g, \"yb\": %.8g," \
	$td_deck $ts_deck $tw_deck $bh_deck $yb_deck]
puts $jsonFd [format "  \"A_frame\": %.8g, \"I_frame\": %.8g, \"Ec\": %.8g," \
	$A_deck_frame $I_deck_frame $Ec_deck]
puts $jsonFd [format "  \"A_link\": %.8g, \"I_link\": %.8g," $A_link $I_link]
puts $jsonFd [format "  \"soffitCL_node\": %d," $nodeTag_deck_BC]
puts $jsonFd "  \"nodes\": \["
set idx 0
foreach nm $deckNames {
	incr idx
	lassign $deckNodeXY($nm) x y
	set jsonComma [expr {($idx < $nDeck) ? "," : ""}]
	puts $jsonFd [format "    \{\"name\": \"%s\", \"tag\": %d, \"x\": %.8g, \"y\": %.8g, \"m\": %.8g, \"mx\": %.8g, \"my\": %.8g, \"Irot\": %.8g\}%s" \
		$nm $deckNodeTag($nm) $x $y $deckM($nm) $deckM($nm) $deckM($nm) $deckIrot($nm) $jsonComma]
}
puts $jsonFd "  \],"
puts $jsonFd "  \"members\": \["
set nMem [llength $deckMembers]
set idx 0
foreach mem $deckMembers {
	incr idx
	lassign $mem a b kind
	set jsonComma [expr {($idx < $nMem) ? "," : ""}]
	puts $jsonFd [format "    \{\"i\": \"%s\", \"j\": \"%s\", \"kind\": \"%s\"\}%s" \
		$a $b $kind $jsonComma]
}
puts $jsonFd "  \]"
puts $jsonFd "\}"
close $jsonFd
