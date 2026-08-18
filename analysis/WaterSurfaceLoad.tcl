# analysis/WaterSurfaceLoad.tcl
# Goals: ponding on free soil top (y~0). No-op if h_water<=0.
# Call after BuildSoilMesh, before SoilGravity.
#
# =====================================================================
# 6. LOADS
# =====================================================================

if {![info exists h_water] || ![info exists rho_w] || ![info exists gravity_accel]} {
	error "WaterSurfaceLoad.tcl: source Parameters.tcl first"
}
if {![info exists t_soil]} {
	error "WaterSurfaceLoad.tcl: source Parameters.tcl first (t_soil)"
}
if {![info exists soilEleTags]} {
	error "WaterSurfaceLoad.tcl: BuildSoilMesh.tcl first"
}

if {$h_water <= 0.0} {
	set soilWaterLoadON 0
	return
}

set p [expr {$rho_w*$gravity_accel*$h_water}]
array unset fyWater
array set fyWater {}
set nEdge 0
set Fsum 0.0
set yTol 1.0e-6

foreach e $soilEleTags {
	set nodes [eleNodes $e]
	# Construction order: n1(TL) n4(BL) n3(BR) n2(TR) -> top edge n1-n2
	set nL [lindex $nodes 0]
	set nR [lindex $nodes 3]
	set yL [lindex [nodeCoord $nL] 1]
	set yR [lindex [nodeCoord $nR] 1]
	if {abs($yL) > $yTol || abs($yR) > $yTol} { continue }

	set xL [lindex [nodeCoord $nL] 0]
	set xR [lindex [nodeCoord $nR] 0]
	set L [expr {abs($xR - $xL)}]
	if {$L < 1.0e-12} { continue }

	set Fy [expr {-0.5*$p*$L*$t_soil}]
	if {![info exists fyWater($nL)]} { set fyWater($nL) 0.0 }
	if {![info exists fyWater($nR)]} { set fyWater($nR) 0.0 }
	set fyWater($nL) [expr {$fyWater($nL) + $Fy}]
	set fyWater($nR) [expr {$fyWater($nR) + $Fy}]
	set Fsum [expr {$Fsum + 2.0*$Fy}]
	incr nEdge
}

if {$nEdge < 1} {
	error "WaterSurfaceLoad.tcl: no y=0 top edges found"
}

timeSeries Linear $tsTag_hWater
# pattern Plain $patternTag $tsTag { load $nodeTag $Fx $Fy $Mz }
# ponding: vertical force on y=0 top-edge nodes
pattern Plain $patternTag_hWater $tsTag_hWater {
	foreach n [array names fyWater] {
		load $n 0.0 $fyWater($n)
	}
}

set soilWaterLoadON 1
puts [format "----- Water surface load: h=%.3f m  p=%.3e Pa  edges=%d  SumFy=%.3e N (t=t_soil) -----" \
	$h_water $p $nEdge $Fsum]
