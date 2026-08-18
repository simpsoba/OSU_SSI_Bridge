# soil/BuildShinLysmer.tcl
# Goals: three Lysmer dashpots + 2c v (Joyner & Chen 1975).
# Call from ActivateEQBoundary after gravity. Do not call during SoilGravity.
# Create the EQ analysis afterward (constraints Transformation).
#
# =====================================================================
# 5. BOUNDARY CONDITIONS / CONSTRAINTS
# =====================================================================
# =====================================================================
# 6. LOADS
# =====================================================================

if {![info exists soilBoundary] || $soilBoundary ne "Shin"} {
	error "BuildShinLysmer.tcl: soilBoundary must be Shin"
}
if {![info exists soil_nX] || ![info exists tsTag_velBase]} {
	error "BuildShinLysmer.tcl: mesh + BuildVelSeries first"
}

# Drop any live analysis so SP/MP edits are not under a Transformation map
wipeAnalysis

set nX $soil_nX
set nY $soil_nY
set iyBot [expr {$nY - 1}]
set xFF_inner $L_half

# Classify base nodes: FF vs near-field
set baseNF {}
set nL_FF ""
set nR_FF ""
for {set ix 0} {$ix < $nX} {incr ix} {
	set nTag [expr {$nodeTag_soil_base + $ix*100 + $iyBot}]
	if {[lsearch -exact [getNodeTags] $nTag] < 0} { continue }
	set x [lindex $soilXs $ix]
	if {$x <= -$xFF_inner + 1.0e-6} {
		# left FF column (prefer outermost)
		if {$nL_FF eq "" || $x < [lindex [nodeCoord $nL_FF] 0]} {
			set nL_FF $nTag
		}
	} elseif {$x >= $xFF_inner - 1.0e-6} {
		if {$nR_FF eq "" || $x > [lindex [nodeCoord $nR_FF] 0]} {
			set nR_FF $nTag
		}
	} else {
		lappend baseNF $nTag
	}
}
if {$nL_FF eq "" || $nR_FF eq "" || [llength $baseNF] < 1} {
	error "BuildShinLysmer.tcl: could not find FF / near-field base nodes"
}

# Primary NF base = closest to x=0 (retained for equalDOF; Lysmer attaches here)
set nPrimary [lindex $baseNF 0]
set xBest 1.0e99
foreach n $baseNF {
	set ax [expr {abs([lindex [nodeCoord $n] 0])}]
	if {$ax < $xBest} {
		set xBest $ax
		set nPrimary $n
	}
}

# Release UX on all base nodes (UY SP from gravity remains)
set nRemoved 0
for {set ix 0} {$ix < $nX} {incr ix} {
	set nTag [expr {$nodeTag_soil_base + $ix*100 + $iyBot}]
	if {[lsearch -exact [getNodeTags] $nTag] < 0} { continue }
	if {[catch {remove sp $nTag 1} err]} {
		error "BuildShinLysmer.tcl: remove sp $nTag 1 failed: $err"
	}
	incr nRemoved
}

# Near-field base UX: retained = nPrimary, constrained = other NF base nodes
# equalDOF $rNodeTag $cNodeTag $dof1 ...
foreach n $baseNF {
	if {$n != $nPrimary} {
		equalDOF $nPrimary $n 1
	}
}

set L_nf [expr {2.0*$xFF_inner}]
set c_nf [expr {$rockRho*$rockVs*$L_nf*$t_soil}]
set c_FF [expr {$rockRho*$rockVs*$w_FF*$t_FF}]
set F_nf [expr {2.0*$c_nf}]
set F_FF [expr {2.0*$c_FF}]

set matNF [expr {$matTag_lysmer_base + 1}]
set matFF [expr {$matTag_lysmer_base + 2}]
# uniaxialMaterial Viscous $matTag $C $alpha
uniaxialMaterial Viscous $matNF $c_nf 1.0
uniaxialMaterial Viscous $matFF $c_FF 1.0

# Dashpot: fully fixed ghost -- Viscous dir 1 -- soil base node.
# Both ends ndf=2; no intermediate free node / equalDOF.
# element zeroLength $eleTag $iNode $jNode -mat $matTag -dir $dir
proc addLysmerDashpot {eleTag matTag soilNode x y} {
	upvar 1 nodeTag_bnd_base nodeTag_bnd_base
	set nFix [expr {$nodeTag_bnd_base + $eleTag}]
	node $nFix $x $y
	fix $nFix 1 1
	element zeroLength $eleTag $nFix $soilNode -mat $matTag -dir 1
}

set e $eleTag_bnd_base
set xy [nodeCoord $nPrimary]
addLysmerDashpot [incr e] $matNF $nPrimary [lindex $xy 0] [lindex $xy 1]
set eleLysmer_nf $e

set xy [nodeCoord $nL_FF]
addLysmerDashpot [incr e] $matFF $nL_FF [lindex $xy 0] [lindex $xy 1]
set eleLysmer_L $e

set xy [nodeCoord $nR_FF]
addLysmerDashpot [incr e] $matFF $nR_FF [lindex $xy 0] [lindex $xy 1]
set eleLysmer_R $e
set eleTag_bnd_last $e

if {![info exists patternTag_lysmer]} {
	error "BuildShinLysmer.tcl: patternTag_lysmer missing (source Parameters.tcl)"
}
catch {remove loadPattern $patternTag_lysmer}
# pattern Plain $patternTag $tsTag { load $nodeTag $Fx $Fy ... }
eval "pattern Plain $patternTag_lysmer $tsTag_velBase {
	load $nPrimary $F_nf 0.0
	load $nL_FF    $F_FF 0.0
	load $nR_FF    $F_FF 0.0
}"

set soilLysmerON 1
# Sketch dump: [ele, x, y, role, c]  role = NF | LFF | RFF
set lysmerDashpotDump {}
lappend lysmerDashpotDump [list $eleLysmer_nf \
	[lindex [nodeCoord $nPrimary] 0] [lindex [nodeCoord $nPrimary] 1] NF $c_nf]
lappend lysmerDashpotDump [list $eleLysmer_L \
	[lindex [nodeCoord $nL_FF] 0] [lindex [nodeCoord $nL_FF] 1] LFF $c_FF]
lappend lysmerDashpotDump [list $eleLysmer_R \
	[lindex [nodeCoord $nR_FF] 0] [lindex [nodeCoord $nR_FF] 1] RFF $c_FF]

puts [format "----- Shin Lysmer  NF c=%.3e  FF c=%.3e  nBase=%d -----" \
	$c_nf $c_FF $nRemoved]
