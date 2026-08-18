# BuildPileCapNodes.tcl
# Units: N, m, s
#
# Goals: 2D transverse pile cap -- stiff elasticBeamColumn frame.
# Half-spacing s = s_pile_cap. Frame and face springs stop at the outer
# pile axes (+/-s). Physical W_cap overhang is tributary mass only.
#
# T/M/B = top / mid / bot (y = 0, -H/2, -H)
# L/C/R = left / center / right (pile axes at x = +/-s, 0)
# BML/BMR = bot mid-bay, x = +/-s/2
#
#   y=0     TL ---------- TC ---------- TR
#            |           / | \           |
#   y=-H/2  ML ---------- MC ---------- MR
#            |         /    |    \        |
#   y=-H    BL -- BML -- BC -- BMR -- BR
#           -s   -s/2     0    s/2     s
#
# TC = pier base (IncrMass).  BL, BC, BR = pile heads.
# Face Py/Tz on TL/ML/BL and TR/MR/BR.
#
# Mass (nodal only): tributary m_i = dens*dx*dy*L_cap, then scale to m_cap.
# Irot_i = m_i (dx^2 + dy^2)/12  (that rectangle about the node).
# dx from midpoints, outer to +/- W/2. dy = H/4, H/2, H/4.
#
# Expects: Parameters.tcl, model domain, BuildPierNodes.tcl already sourced.
# Frame beams: BuildPileCapElements.tcl after nodes.

if {![info exists H_cap] || ![info exists nodeTag_cap_TC] || ![info exists A_cap]} {
	error "BuildPileCapNodes.tcl: source Parameters.tcl first"
}
if {![info exists structNodeTags]} { set structNodeTags {} }
if {$nodeTag_cap_TC != $nodeTag_pierBase_capTC} {
	error "BuildPileCapNodes.tcl: nodeTag_cap_TC must equal nodeTag_pierBase_capTC"
}
if {![info exists n_cap_nodes]} {
	error "BuildPileCapNodes.tcl: need n_cap_nodes from Parameters"
}

set scriptDir [file dirname [file normalize [info script]]]
source [file join $scriptDir IncrMass.tcl]

set s $s_pile_cap
set yMid [expr {-0.5*$H_cap}]
set yBot [expr {-$H_cap}]
set xEdgeL [expr {-0.5*$W_cap}]
set xEdgeR [expr { 0.5*$W_cap}]
if {$s > $xEdgeR + 1.0e-9} {
	error "BuildPileCapNodes.tcl: s_pile_cap > W_cap/2 (outer trib dx < 0)"
}

# Tributary dx for each station on a sorted x-row. Outer to xEdgeL / xEdgeR.
# Args: xRow (m, increasing) xEdgeL xEdgeR (m)
# Returns: list of dx (m), same order as xRow
proc capTribWidths {xRow xEdgeL xEdgeR} {
	set n [llength $xRow]
	set dxs {}
	for {set i 0} {$i < $n} {incr i} {
		set x [lindex $xRow $i]
		if {$i == 0} {
			set xL $xEdgeL
		} else {
			set xL [expr {0.5*([lindex $xRow [expr {$i - 1}]] + $x)}]
		}
		if {$i == $n - 1} {
			set xR $xEdgeR
		} else {
			set xR [expr {0.5*($x + [lindex $xRow [expr {$i + 1}]])}]
		}
		set dx [expr {$xR - $xL}]
		if {$dx <= 0.0} {
			error [format "capTribWidths: dx=%.4e at x=%.4f" $dx $x]
		}
		lappend dxs $dx
	}
	return $dxs
}

set xRowTop [list [expr {-$s}] 0.0 $s]
set xRowBot [list [expr {-$s}] [expr {-0.5*$s}] 0.0 \
	[expr {0.5*$s}] $s]
set dxTop [capTribWidths $xRowTop $xEdgeL $xEdgeR]
set dxBot [capTribWidths $xRowBot $xEdgeL $xEdgeR]
set dyT [expr {0.25*$H_cap}];             # m, top row
set dyM [expr {0.50*$H_cap}];             # m, mid row
set dyB [expr {0.25*$H_cap}];             # m, bot row

array unset capM
array unset capIrot
set capMassTags {}

# m = dens*dx*dy*L ; Irot = m (dx^2+dy^2)/12
proc capSetTrib {tag dx dy} {
	global capM capIrot capMassTags dens_cap L_cap
	set capM($tag) [expr {$dens_cap*$dx*$dy*$L_cap}]
	set capIrot($tag) [expr {$capM($tag)*($dx*$dx + $dy*$dy)/12.0}]
	lappend capMassTags $tag
}

foreach tag [list $nodeTag_cap_TL $nodeTag_cap_TC $nodeTag_cap_TR] \
		dx $dxTop {
	capSetTrib $tag $dx $dyT
}
foreach tag [list $nodeTag_cap_ML $nodeTag_cap_MC $nodeTag_cap_MR] \
		dx $dxTop {
	capSetTrib $tag $dx $dyM
}
foreach tag [list $nodeTag_cap_BL $nodeTag_cap_BML $nodeTag_cap_BC \
		$nodeTag_cap_BMR $nodeTag_cap_BR] \
		dx $dxBot {
	capSetTrib $tag $dx $dyB
}
if {[llength $capMassTags] != $n_cap_nodes} {
	error "BuildPileCapNodes.tcl: n_cap_nodes=$n_cap_nodes but trib list has [llength $capMassTags]"
}

set mSum 0.0
foreach tag $capMassTags { set mSum [expr {$mSum + $capM($tag)}] }
if {$mSum <= 0.0} {
	error "BuildPileCapNodes.tcl: zero tributary mass sum"
}
set capMscale [expr {$m_cap/$mSum}]
set I_steiner 0.0
set Irot_sum 0.0
set mMin 1.0e99
set mMax 0.0
set IrotMin 1.0e99
set IrotMax 0.0
foreach tag $capMassTags {
	set capM($tag) [expr {$capMscale*$capM($tag)}]
	set capIrot($tag) [expr {$capMscale*$capIrot($tag)}]
	set Irot_sum [expr {$Irot_sum + $capIrot($tag)}]
	if {$capM($tag) < $mMin} { set mMin $capM($tag) }
	if {$capM($tag) > $mMax} { set mMax $capM($tag) }
	if {$capIrot($tag) < $IrotMin} { set IrotMin $capIrot($tag) }
	if {$capIrot($tag) > $IrotMax} { set IrotMax $capIrot($tag) }
}

set capXY [list \
	[list $nodeTag_cap_TL [expr {-$s}] 0.0] \
	[list $nodeTag_cap_TC 0.0 0.0] \
	[list $nodeTag_cap_TR $s 0.0] \
	[list $nodeTag_cap_ML [expr {-$s}] $yMid] \
	[list $nodeTag_cap_MC 0.0 $yMid] \
	[list $nodeTag_cap_MR $s $yMid] \
	[list $nodeTag_cap_BL [expr {-$s}] $yBot] \
	[list $nodeTag_cap_BML [expr {-0.5*$s}] $yBot] \
	[list $nodeTag_cap_BC 0.0 $yBot] \
	[list $nodeTag_cap_BMR [expr {0.5*$s}] $yBot] \
	[list $nodeTag_cap_BR $s $yBot] \
]
foreach row $capXY {
	lassign $row tag x y
	set dy [expr {$y - $yMid}]
	set I_steiner [expr {$I_steiner + $capM($tag)*($x*$x + $dy*$dy)}]
}
set Iz_cap [expr {$m_cap*($W_cap*$W_cap + $H_cap*$H_cap)/12.0}];  # kg*m^2

# =====================================================================
# 2. MODEL BUILDER / NODES
# =====================================================================
# node $tag $x $y -mass $mx $my $mRz
# --- top y=0: TL pile line, TC pier base, TR ---
node $nodeTag_cap_TL [expr {-$s}] 0.0 \
	-mass $capM($nodeTag_cap_TL) $capM($nodeTag_cap_TL) $capIrot($nodeTag_cap_TL)
IncrMass $nodeTag_cap_TC $capM($nodeTag_cap_TC) $capM($nodeTag_cap_TC) \
	$capIrot($nodeTag_cap_TC)
node $nodeTag_cap_TR $s 0.0 \
	-mass $capM($nodeTag_cap_TR) $capM($nodeTag_cap_TR) $capIrot($nodeTag_cap_TR)

# --- mid y=-H/2: ML, MC, MR ---
node $nodeTag_cap_ML [expr {-$s}] $yMid \
	-mass $capM($nodeTag_cap_ML) $capM($nodeTag_cap_ML) $capIrot($nodeTag_cap_ML)
node $nodeTag_cap_MC 0.0 $yMid \
	-mass $capM($nodeTag_cap_MC) $capM($nodeTag_cap_MC) $capIrot($nodeTag_cap_MC)
node $nodeTag_cap_MR $s $yMid \
	-mass $capM($nodeTag_cap_MR) $capM($nodeTag_cap_MR) $capIrot($nodeTag_cap_MR)

# --- bot y=-H: BL pile, BML mid-bay, BC, BMR, BR pile ---
node $nodeTag_cap_BL [expr {-$s}] $yBot \
	-mass $capM($nodeTag_cap_BL) $capM($nodeTag_cap_BL) $capIrot($nodeTag_cap_BL)
node $nodeTag_cap_BML [expr {-0.5*$s}] $yBot \
	-mass $capM($nodeTag_cap_BML) $capM($nodeTag_cap_BML) $capIrot($nodeTag_cap_BML)
node $nodeTag_cap_BC 0.0 $yBot \
	-mass $capM($nodeTag_cap_BC) $capM($nodeTag_cap_BC) $capIrot($nodeTag_cap_BC)
node $nodeTag_cap_BMR [expr {0.5*$s}] $yBot \
	-mass $capM($nodeTag_cap_BMR) $capM($nodeTag_cap_BMR) $capIrot($nodeTag_cap_BMR)
node $nodeTag_cap_BR $s $yBot \
	-mass $capM($nodeTag_cap_BR) $capM($nodeTag_cap_BR) $capIrot($nodeTag_cap_BR)

# Soffit stations for springs (tag, x) left->right -- exported for BuildSoilSprings
set capSoffitStations [list \
	[list $nodeTag_cap_BL [expr {-$s}]] \
	[list $nodeTag_cap_BML [expr {-0.5*$s}]] \
	[list $nodeTag_cap_BC 0.0] \
	[list $nodeTag_cap_BMR [expr {0.5*$s}]] \
	[list $nodeTag_cap_BR $s] \
]

# TC = pier bot: already on structNodeTags
lappend structNodeTags \
	$nodeTag_cap_TL $nodeTag_cap_TR \
	$nodeTag_cap_ML $nodeTag_cap_MC $nodeTag_cap_MR \
	$nodeTag_cap_BL $nodeTag_cap_BML $nodeTag_cap_BC \
	$nodeTag_cap_BMR $nodeTag_cap_BR
