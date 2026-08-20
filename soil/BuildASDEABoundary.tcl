# soil/BuildASDEABoundary.tcl
# Goals: ASDAbsorbingBoundary2D L/BL/B/BR/R at stage 0.
# Call after near-field mesh. Stage 1: ActivateEQBoundary after gravity.
#
# =====================================================================
# 2. MODEL BUILDER / NODES
# =====================================================================
# =====================================================================
# 4. ELEMENTS
# =====================================================================
# =====================================================================
# 5. BOUNDARY CONDITIONS / CONSTRAINTS
# =====================================================================

if {![info exists soilBoundary] || $soilBoundary ne "ASDEA"} {
	error "BuildASDEABoundary.tcl: soilBoundary must be ASDEA"
}
if {![info exists soil_nX] || ![info exists tsTag_velBase]} {
	error "BuildASDEABoundary.tcl: mesh + BuildVelSeries first"
}
if {![info exists rockG] || ![info exists asdeaNu]} {
	error "BuildASDEABoundary.tcl: rock / asdeaNu from Parameters"
}

# Keep ASDEA outer tags above the continuum mesh.
set maxN [tcl::mathfunc::max {*}[getNodeTags]]
if {[ensureAbove nodeTag_bnd_base $maxN]} {
	puts [format "----- ASDEA tags  nodes -> %d (above max node %d) -----" \
		$nodeTag_bnd_base $maxN]
}
set maxE [tcl::mathfunc::max {*}[getEleTags]]
if {[ensureAbove eleTag_bnd_base $maxE]} {
	puts [format "----- ASDEA tags  eles -> %d (above max ele %d) -----" \
		$eleTag_bnd_base $maxE]
}

set nX $soil_nX
set nY $soil_nY
set ixL 0
set ixR [expr {$nX - 1}]
set iyBot [expr {$nY - 1}]
set xL [lindex $soilXs $ixL]
set xR [lindex $soilXs $ixR]
set yBot [lindex $soilYs $iyBot]

if {$w_FF <= 0} {
	error "BuildASDEABoundary.tcl: w_FF must be > 0 (ASD ring width = Shin FF strip)"
}
set hext $w_FF

# Outer node tags packed from nY / nX so a dense mesh cannot collide:
#   left  iy:  nodeTag_bnd_base + iy
#   right iy:  nodeTag_bnd_base + nY + iy
#   bottom ix: nodeTag_bnd_base + 2*nY + ix
#   BL: + 2*nY + nX; BR: +1 after BL
proc bndNodeExists {tag} {
	expr {[lsearch -exact [getNodeTags] $tag] >= 0}
}

set bndLeft0 $nodeTag_bnd_base
set bndRight0 [expr {$nodeTag_bnd_base + $nY}]
set bndBot0 [expr {$nodeTag_bnd_base + 2*$nY}]
set nBL [expr {$bndBot0 + $nX}]
set nBR [expr {$nBL + 1}]

for {set iy 0} {$iy < $nY} {incr iy} {
	set y [lindex $soilYs $iy]
	node [expr {$bndLeft0 + $iy}] [expr {$xL - $hext}] $y
	node [expr {$bndRight0 + $iy}] [expr {$xR + $hext}] $y
}
for {set ix 0} {$ix < $nX} {incr ix} {
	set x [lindex $soilXs $ix]
	node [expr {$bndBot0 + $ix}] $x [expr {$yBot - $hext}]
}
node $nBL [expr {$xL - $hext}] [expr {$yBot - $hext}]
node $nBR [expr {$xR + $hext}] [expr {$yBot - $hext}]

set soilEleBndTags {}
set e $eleTag_bnd_base
set thick $t_soil

# ---- Left / Right vertical faces ----
for {set iy 0} {$iy < $nY - 1} {incr iy} {
	set yT [lindex $soilYs $iy]
	set yB [lindex $soilYs [expr {$iy + 1}]]
	set yc [expr {0.5*($yT + $yB)}]
	set G $soilG0($iy)
	set rho $soilRho($iy)
	set nu $asdeaNu

	# Left: N1=outer BL, N2=soil BL, N3=soil TL, N4=outer TL
	# element ASDAbsorbingBoundary2D $eleTag $n1 $n2 $n3 $n4 $G $v $rho $thickness $bType <-fx $tsTag>
	set n1 [expr {$bndLeft0 + ($iy + 1)}]
	set n2 [soilNodeTag $ixL [expr {$iy + 1}]]
	set n3 [soilNodeTag $ixL $iy]
	set n4 [expr {$bndLeft0 + $iy}]
	incr e
	element ASDAbsorbingBoundary2D $e $n1 $n2 $n3 $n4 $G $nu $rho $thick L
	lappend soilEleBndTags $e

	# Right: N1=soil BL, N2=outer BR, N3=outer TR, N4=soil TL
	set n1 [soilNodeTag $ixR [expr {$iy + 1}]]
	set n2 [expr {$bndRight0 + ($iy + 1)}]
	set n3 [expr {$bndRight0 + $iy}]
	set n4 [soilNodeTag $ixR $iy]
	incr e
	element ASDAbsorbingBoundary2D $e $n1 $n2 $n3 $n4 $G $nu $rho $thick R
	lappend soilEleBndTags $e
}

# ---- Bottom face (rock) + corners ----
set G $rockG
set nu $rockNu
set rho $rockRho

for {set ix 0} {$ix < $nX - 1} {incr ix} {
	set n1 [expr {$bndBot0 + $ix}]
	set n2 [expr {$bndBot0 + ($ix + 1)}]
	set n3 [soilNodeTag [expr {$ix + 1}] $iyBot]
	set n4 [soilNodeTag $ix $iyBot]
	# skip if soil corner nodes missing (should not happen on base)
	set ok 1
	foreach nn [list $n3 $n4] {
		if {![bndNodeExists $nn]} { set ok 0; break }
	}
	if {!$ok} { continue }
	incr e
	element ASDAbsorbingBoundary2D $e $n1 $n2 $n3 $n4 $G $nu $rho $thick B \
		-fx $tsTag_velBase
	lappend soilEleBndTags $e
}

# BL corner
set n1 $nBL
set n2 [expr {$bndBot0 + $ixL}]
set n3 [soilNodeTag $ixL $iyBot]
set n4 [expr {$bndLeft0 + $iyBot}]
incr e
element ASDAbsorbingBoundary2D $e $n1 $n2 $n3 $n4 $G $nu $rho $thick BL \
	-fx $tsTag_velBase
lappend soilEleBndTags $e

# BR corner
set n1 [expr {$bndBot0 + $ixR}]
set n2 $nBR
set n3 [expr {$bndRight0 + $iyBot}]
set n4 [soilNodeTag $ixR $iyBot]
incr e
element ASDAbsorbingBoundary2D $e $n1 $n2 $n3 $n4 $G $nu $rho $thick BR \
	-fx $tsTag_velBase
lappend soilEleBndTags $e

set eleTag_bnd_last $e
set soilBndStage 0
