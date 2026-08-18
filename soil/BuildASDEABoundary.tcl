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

# Outer node tags: left/right columns, bottom row, two corners
#   left  iy:  nodeTag_bnd_base + 1000 + iy
#   right iy:  nodeTag_bnd_base + 2000 + iy
#   bottom ix: nodeTag_bnd_base + 3000 + ix
#   BL: +4000; BR: +4001
proc bndNodeExists {tag} {
	expr {[lsearch -exact [getNodeTags] $tag] >= 0}
}

for {set iy 0} {$iy < $nY} {incr iy} {
	set y [lindex $soilYs $iy]
	set nL [expr {$nodeTag_bnd_base + 1000 + $iy}]
	set nR [expr {$nodeTag_bnd_base + 2000 + $iy}]
	node $nL [expr {$xL - $hext}] $y
	node $nR [expr {$xR + $hext}] $y
}
for {set ix 0} {$ix < $nX} {incr ix} {
	set x [lindex $soilXs $ix]
	set nB [expr {$nodeTag_bnd_base + 3000 + $ix}]
	node $nB $x [expr {$yBot - $hext}]
}
set nBL [expr {$nodeTag_bnd_base + 4000}]
set nBR [expr {$nodeTag_bnd_base + 4001}]
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
	set n1 [expr {$nodeTag_bnd_base + 1000 + ($iy + 1)}]
	set n2 [expr {$nodeTag_soil_base + $ixL*100 + ($iy + 1)}]
	set n3 [expr {$nodeTag_soil_base + $ixL*100 + $iy}]
	set n4 [expr {$nodeTag_bnd_base + 1000 + $iy}]
	incr e
	element ASDAbsorbingBoundary2D $e $n1 $n2 $n3 $n4 $G $nu $rho $thick L
	lappend soilEleBndTags $e

	# Right: N1=soil BL, N2=outer BR, N3=outer TR, N4=soil TL
	set n1 [expr {$nodeTag_soil_base + $ixR*100 + ($iy + 1)}]
	set n2 [expr {$nodeTag_bnd_base + 2000 + ($iy + 1)}]
	set n3 [expr {$nodeTag_bnd_base + 2000 + $iy}]
	set n4 [expr {$nodeTag_soil_base + $ixR*100 + $iy}]
	incr e
	element ASDAbsorbingBoundary2D $e $n1 $n2 $n3 $n4 $G $nu $rho $thick R
	lappend soilEleBndTags $e
}

# ---- Bottom face (rock) + corners ----
set G $rockG
set nu $rockNu
set rho $rockRho

for {set ix 0} {$ix < $nX - 1} {incr ix} {
	set n1 [expr {$nodeTag_bnd_base + 3000 + $ix}]
	set n2 [expr {$nodeTag_bnd_base + 3000 + ($ix + 1)}]
	set n3 [expr {$nodeTag_soil_base + ($ix + 1)*100 + $iyBot}]
	set n4 [expr {$nodeTag_soil_base + $ix*100 + $iyBot}]
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
set n2 [expr {$nodeTag_bnd_base + 3000 + $ixL}]
set n3 [expr {$nodeTag_soil_base + $ixL*100 + $iyBot}]
set n4 [expr {$nodeTag_bnd_base + 1000 + $iyBot}]
incr e
element ASDAbsorbingBoundary2D $e $n1 $n2 $n3 $n4 $G $nu $rho $thick BL \
	-fx $tsTag_velBase
lappend soilEleBndTags $e

# BR corner
set n1 [expr {$nodeTag_bnd_base + 3000 + $ixR}]
set n2 $nBR
set n3 [expr {$nodeTag_bnd_base + 2000 + $iyBot}]
set n4 [expr {$nodeTag_soil_base + $ixR*100 + $iyBot}]
incr e
element ASDAbsorbingBoundary2D $e $n1 $n2 $n3 $n4 $G $nu $rho $thick BR \
	-fx $tsTag_velBase
lappend soilEleBndTags $e

set eleTag_bnd_last $e
set soilBndStage 0
