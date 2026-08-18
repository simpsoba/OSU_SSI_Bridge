# CircleStripFibers.tcl
# Goals: circular / annular section -> horizontal fiber strips for a 2D
# Fiber section (in-plane bending; strain varies with section y).
# Graded default: edgeFrac=1/6 each side (e.g. nFibers=21, nEdge=5 -> 5+11+5).

# Strip area and first moment about y=0 for a circular slice [y1, y2].
# Args: radius y1 y2 (m)
# Returns: {area firstMoment} (m2, m3)
proc circleStripProperties {radius y1 y2} {
	set R [expr {double($radius)}]

	if {$R <= 0.0} {
		return [list 0.0 0.0]
	}

	# Clip the strip to the circle so we never integrate outside +/-R
	set a [expr {max(double($y1), -$R)}]
	set b [expr {min(double($y2),  $R)}]

	if {$b <= $a} {
		return [list 0.0 0.0]
	}

	set R2 [expr {$R*$R}]

	# Half-chord lengths at the strip edges: sqrt(R^2 - y^2)
	set sa2 [expr {max(0.0, $R2 - $a*$a)}]
	set sb2 [expr {max(0.0, $R2 - $b*$b)}]
	set sa [expr {sqrt($sa2)}]
	set sb [expr {sqrt($sb2)}]

	# asin args must stay in [-1, 1] (floating-point guard)
	set argA [expr {max(-1.0, min(1.0, $a/$R))}]
	set argB [expr {max(-1.0, min(1.0, $b/$R))}]

	# Antiderivative of the circle chord width -> strip area = F(b) - F(a)
	set FA [expr {$a*$sa + $R2*asin($argA)}]
	set FB [expr {$b*$sb + $R2*asin($argB)}]
	set area [expr {$FB - $FA}]

	# First moment int y dA about y = 0  ->  later centroid = firstMoment / area
	set firstMoment [expr {
		(2.0/3.0) *
		(pow($sa2, 1.5) - pow($sb2, 1.5))
	}]

	return [list $area $firstMoment]
}

# Uniform strips from -Ro to +Ro. Tube: outer strip minus inner (Ri=0 -> solid).
# Args: outerRadius innerRadius nFibers ?yCenter? ?zCenter? (m)
# Returns: list of {yCentroid zCentroid area}
proc circularTubeFiberStrips {
	outerRadius innerRadius nFibers
	{yCenter 0.0} {zCenter 0.0}
} {
	set Ro [expr {double($outerRadius)}]
	set Ri [expr {double($innerRadius)}]

	if {$Ro <= 0.0} {
		error "outerRadius must be greater than zero"
	}
	if {$Ri < 0.0} {
		error "innerRadius cannot be negative"
	}
	if {$Ri >= $Ro} {
		error "innerRadius must be smaller than outerRadius"
	}
	if {$nFibers < 1 || $nFibers != int($nFibers)} {
		error "nFibers must be a positive integer"
	}

	set dy [expr {2.0*$Ro/double($nFibers)}];	# strip thickness
	set fibers {}

	for {set i 0} {$i < $nFibers} {incr i} {
		set y1 [expr {-$Ro + $i*$dy}]
		set y2 [expr {-$Ro + ($i + 1)*$dy}]

		# Pin the first/last edges exactly to +/-Ro (avoid tiny gaps from float)
		if {$i == 0} {
			set y1 [expr {-$Ro}]
		}
		if {$i == $nFibers - 1} {
			set y2 $Ro
		}

		# Solid outer disk strip
		lassign [circleStripProperties $Ro $y1 $y2] outerArea outerMoment

		# Hollow out the inner disk (Ri = 0 -> solid section)
		if {$Ri > 0.0} {
			lassign [circleStripProperties $Ri $y1 $y2] innerArea innerMoment
		} else {
			set innerArea   0.0
			set innerMoment 0.0
		}

		set area [expr {$outerArea - $innerArea}]
		set firstMoment [expr {$outerMoment - $innerMoment}]

		if {$area <= 0.0} {
			continue;	# strip missed the wall (near the equator of a thin tube, rare)
		}

		# Strip centroid in section coords; z stays 0 for 2D bending
		set yLocal [expr {$firstMoment/$area}]
		set yCentroid [expr {$yCenter + $yLocal}]
		set zCentroid $zCenter

		lappend fibers [list $yCentroid $zCentroid $area]
	}

	return $fibers
}

# Graded strip spacing; denser at extreme fibers (edgeFrac each side).
# Args: outerRadius innerRadius nFibers ?nEdge? ?edgeFrac? ?yCenter? ?zCenter?
# Returns: list of {yCentroid zCentroid area}
proc circularTubeFiberStripsGraded {
	outerRadius innerRadius nFibers
	{nEdge 5}
	{edgeFrac 0.1666666667}
	{yCenter 0.0} {zCenter 0.0}
} {
	set Ro [expr {double($outerRadius)}]
	set Ri [expr {double($innerRadius)}]
	set N  [expr {int($nFibers)}]
	set nEdgeEach [expr {int($nEdge)}]
	set ef [expr {double($edgeFrac)}]

	if {$Ro <= 0.0} {
		error "outerRadius must be greater than zero"
	}
	if {$Ri < 0.0 || $Ri >= $Ro} {
		error "innerRadius must satisfy 0 <= Ri < Ro"
	}
	if {$nEdgeEach < 1} {
		error "nEdge must be >= 1"
	}
	set nMid [expr {$N - 2*$nEdgeEach}]
	if {$nMid < 1} {
		error "nFibers=$N too small for nEdge=$nEdgeEach on each side (need mid >= 1)"
	}
	if {$ef <= 0.0 || $ef >= 0.5} {
		error "edgeFrac must be in (0, 0.5)"
	}

	# Depth spans of the three zones along the diameter
	set D [expr {2.0*$Ro}]
	set Ledge [expr {$ef*$D}];			# each outer band (default 1/6 of D)
	set Lmid  [expr {(1.0 - 2.0*$ef)*$D}];	# middle band (default 2/3 of D)

	set yBot  [expr {-$Ro}]
	set yE1   [expr {-$Ro + $Ledge}];		# end of bottom edge zone
	set yE2   [expr { $Ro - $Ledge}];		# start of top edge zone
	set yTop  $Ro

	# yBreaks = sorted strip boundaries from bottom to top
	set yBreaks [list $yBot]

	# Bottom edge zone: nEdge equal steps over Ledge
	for {set i 1} {$i <= $nEdgeEach} {incr i} {
		lappend yBreaks [expr {$yBot + $Ledge*double($i)/double($nEdgeEach)}]
	}
	set yBreaks [lreplace $yBreaks end end $yE1];	# land exactly on yE1

	# Middle zone: nMid equal steps over Lmid
	for {set i 1} {$i <= $nMid} {incr i} {
		lappend yBreaks [expr {$yE1 + $Lmid*double($i)/double($nMid)}]
	}
	set yBreaks [lreplace $yBreaks end end $yE2]

	# Top edge zone: nEdge equal steps over Ledge
	for {set i 1} {$i <= $nEdgeEach} {incr i} {
		lappend yBreaks [expr {$yE2 + $Ledge*double($i)/double($nEdgeEach)}]
	}
	set yBreaks [lreplace $yBreaks end end $yTop]

	# One fiber per [yBreaks[i], yBreaks[i+1]] interval
	set fibers {}
	set nSeg [expr {[llength $yBreaks] - 1}]
	for {set i 0} {$i < $nSeg} {incr i} {
		set y1 [lindex $yBreaks $i]
		set y2 [lindex $yBreaks [expr {$i + 1}]]

		lassign [circleStripProperties $Ro $y1 $y2] outerArea outerMoment
		if {$Ri > 0.0} {
			lassign [circleStripProperties $Ri $y1 $y2] innerArea innerMoment
		} else {
			set innerArea 0.0
			set innerMoment 0.0
		}
		set area [expr {$outerArea - $innerArea}]
		set firstMoment [expr {$outerMoment - $innerMoment}]
		if {$area <= 0.0} {
			continue
		}
		set yLocal [expr {$firstMoment/$area}]
		lappend fibers [list [expr {$yCenter + $yLocal}] $zCenter $area]
	}
	return $fibers
}

# Same y-breaks as graded tube; split core (r<=Rc) vs cover.
# Args: outerRadius coreRadius nFibers ?nEdge? ?edgeFrac? ?yCenter? ?zCenter?
# Returns: {coreFibers coverFibers} each {yCentroid zCentroid area}
proc circularCoreCoverFiberStripsGraded {
	outerRadius coreRadius nFibers
	{nEdge 5}
	{edgeFrac 0.1666666667}
	{yCenter 0.0} {zCenter 0.0}
} {
	set Ro [expr {double($outerRadius)}]
	set Rc [expr {double($coreRadius)}]
	set N  [expr {int($nFibers)}]
	set nEdgeEach [expr {int($nEdge)}]
	set ef [expr {double($edgeFrac)}]

	if {$Ro <= 0.0 || $Rc <= 0.0 || $Rc >= $Ro} {
		error "need 0 < coreRadius < outerRadius"
	}
	if {$nEdgeEach < 1} {
		error "nEdge must be >= 1"
	}
	set nMid [expr {$N - 2*$nEdgeEach}]
	if {$nMid < 1} {
		error "nFibers=$N too small for nEdge=$nEdgeEach on each side (need mid >= 1)"
	}
	if {$ef <= 0.0 || $ef >= 0.5} {
		error "edgeFrac must be in (0, 0.5)"
	}

	set D [expr {2.0*$Ro}]
	set Ledge [expr {$ef*$D}]
	set Lmid  [expr {(1.0 - 2.0*$ef)*$D}]

	set yBot [expr {-$Ro}]
	set yE1  [expr {-$Ro + $Ledge}]
	set yE2  [expr { $Ro - $Ledge}]
	set yTop $Ro

	set yBreaks [list $yBot]
	for {set i 1} {$i <= $nEdgeEach} {incr i} {
		lappend yBreaks [expr {$yBot + $Ledge*double($i)/double($nEdgeEach)}]
	}
	set yBreaks [lreplace $yBreaks end end $yE1]
	for {set i 1} {$i <= $nMid} {incr i} {
		lappend yBreaks [expr {$yE1 + $Lmid*double($i)/double($nMid)}]
	}
	set yBreaks [lreplace $yBreaks end end $yE2]
	for {set i 1} {$i <= $nEdgeEach} {incr i} {
		lappend yBreaks [expr {$yE2 + $Ledge*double($i)/double($nEdgeEach)}]
	}
	set yBreaks [lreplace $yBreaks end end $yTop]

	set coreFibers {}
	set coverFibers {}
	set nSeg [expr {[llength $yBreaks] - 1}]
	for {set i 0} {$i < $nSeg} {incr i} {
		set y1 [lindex $yBreaks $i]
		set y2 [lindex $yBreaks [expr {$i + 1}]]

		lassign [circleStripProperties $Ro $y1 $y2] Ao Mo
		lassign [circleStripProperties $Rc $y1 $y2] Ac Mc

		if {$Ac > 0.0} {
			set yc [expr {$yCenter + $Mc/$Ac}]
			lappend coreFibers [list $yc $zCenter $Ac]
		}
		set Acover [expr {$Ao - $Ac}]
		set Mcover [expr {$Mo - $Mc}]
		if {$Acover > 0.0} {
			set yv [expr {$yCenter + $Mcover/$Acover}]
			lappend coverFibers [list $yv $zCenter $Acover]
		}
	}
	return [list $coreFibers $coverFibers]
}

# Rebar ring as fibers on z=0; bars with the same section-y are merged.
# Angles match OpenSees layer circ: y = R*cos(theta), startAng at the first bar.
# Args: radius numBars barArea ?yCenter? ?zCenter? ?startAng? (m, -, m2, deg)
# Returns: list of {yCentroid zCentroid area}
proc circularRebarYFibers {
	radius numBars barArea
	{yCenter 0.0} {zCenter 0.0} {startAng 0.0}
} {
	set R [expr {double($radius)}]
	set As [expr {double($barArea)}]
	set n [expr {int($numBars)}]

	if {$R <= 0.0 || $As <= 0.0 || $n < 1} {
		error "circularRebarYFibers: need positive radius, barArea, and numBars"
	}

	set pi [expr {acos(-1.0)}]
	set dAng [expr {360.0/double($n)}]

	# Bucket bars by y; round the key so +/-z pairs land in the same bucket
	array unset areaAtY
	array set areaAtY {}
	array unset yExact
	array set yExact {}

	for {set i 0} {$i < $n} {incr i} {
		set angDeg [expr {double($startAng) + $i*$dAng}]
		set ang [expr {$angDeg*$pi/180.0}]
		set y [expr {$yCenter + $R*cos($ang)}]

		set yKey [format %.8f $y]
		if {[info exists areaAtY($yKey)]} {
			set areaAtY($yKey) [expr {$areaAtY($yKey) + $As}]
		} else {
			set areaAtY($yKey) $As
			set yExact($yKey) $y
		}
	}

	# Emit one fiber per unique y, bottom to top
	set fibers {}
	foreach yKey [lsort -real [array names yExact]] {
		lappend fibers [list $yExact($yKey) 0.0 $areaAtY($yKey)]
	}
	return $fibers
}
