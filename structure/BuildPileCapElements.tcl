# BuildPileCapElements.tcl -- stiff cap frame (explicit beams)
# Call after BuildPileCapNodes.tcl. Tags: eleTag_cap_base.
#
# T/M/B top/mid/bot; L/C/R left/center/right; BML/BMR mid-bay
#   TL------TC------TR
#   ML------MC------MR
#   BL--BML--BC--BMR--BR

if {[info exists capElementsDone] && $capElementsDone} {
	return
}
if {![info exists nodeTag_cap_TC] || ![info exists A_cap]} {
	error "BuildPileCapElements.tcl: source BuildPileCapNodes.tcl first"
}

# =====================================================================
# 4. ELEMENTS
# =====================================================================
geomTransf Linear $transfTag_cap

# element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
set e0 $eleTag_cap_base

# Verticals -- pile axes (two segments each)
element elasticBeamColumn $e0 \
	$nodeTag_cap_TL $nodeTag_cap_ML $A_cap $E_cap $I_cap $transfTag_cap
element elasticBeamColumn [expr {$e0 + 1}] \
	$nodeTag_cap_ML $nodeTag_cap_BL $A_cap $E_cap $I_cap $transfTag_cap
element elasticBeamColumn [expr {$e0 + 2}] \
	$nodeTag_cap_TC $nodeTag_cap_MC $A_cap $E_cap $I_cap $transfTag_cap
element elasticBeamColumn [expr {$e0 + 3}] \
	$nodeTag_cap_MC $nodeTag_cap_BC $A_cap $E_cap $I_cap $transfTag_cap
element elasticBeamColumn [expr {$e0 + 4}] \
	$nodeTag_cap_TR $nodeTag_cap_MR $A_cap $E_cap $I_cap $transfTag_cap
element elasticBeamColumn [expr {$e0 + 5}] \
	$nodeTag_cap_MR $nodeTag_cap_BR $A_cap $E_cap $I_cap $transfTag_cap
# Horizontals -- top, mid, bottom (pile -> ... -> pile)
element elasticBeamColumn [expr {$e0 + 6}] \
	$nodeTag_cap_TL $nodeTag_cap_TC $A_cap $E_cap $I_cap $transfTag_cap
element elasticBeamColumn [expr {$e0 + 7}] \
	$nodeTag_cap_TC $nodeTag_cap_TR $A_cap $E_cap $I_cap $transfTag_cap
element elasticBeamColumn [expr {$e0 + 8}] \
	$nodeTag_cap_ML $nodeTag_cap_MC $A_cap $E_cap $I_cap $transfTag_cap
element elasticBeamColumn [expr {$e0 + 9}] \
	$nodeTag_cap_MC $nodeTag_cap_MR $A_cap $E_cap $I_cap $transfTag_cap
element elasticBeamColumn [expr {$e0 + 10}] \
	$nodeTag_cap_BL $nodeTag_cap_BML $A_cap $E_cap $I_cap $transfTag_cap
element elasticBeamColumn [expr {$e0 + 11}] \
	$nodeTag_cap_BML $nodeTag_cap_BC $A_cap $E_cap $I_cap $transfTag_cap
element elasticBeamColumn [expr {$e0 + 12}] \
	$nodeTag_cap_BC $nodeTag_cap_BMR $A_cap $E_cap $I_cap $transfTag_cap
element elasticBeamColumn [expr {$e0 + 13}] \
	$nodeTag_cap_BMR $nodeTag_cap_BR $A_cap $E_cap $I_cap $transfTag_cap
# Diagonals (X-brace on pile-axis rectangle)
element elasticBeamColumn [expr {$e0 + 14}] \
	$nodeTag_cap_TL $nodeTag_cap_BR $A_cap $E_cap $I_cap $transfTag_cap
element elasticBeamColumn [expr {$e0 + 15}] \
	$nodeTag_cap_TR $nodeTag_cap_BL $A_cap $E_cap $I_cap $transfTag_cap

set eleTag_cap_last [expr {$e0 + 15}]
set capElementsDone 1
