# analysis/BuildVelSeries.tcl
# Goals: outcrop velocity Path (PEER VT2 or dummy zeros). PEER CM/S -> m/s
# before Path -factor. Knobs: gmVelFile, gmVelDT, gmScaleFactor.
# Call from Run.tcl (EQ) or the ASDEA mesh.
#
# =====================================================================
# 6. LOADS
# =====================================================================
# (timeSeries Path used by Shin Lysmer loads and ASDEA -fx)

if {![info exists tsTag_velBase]} {
	error "BuildVelSeries.tcl: source Parameters.tcl first"
}

catch {remove timeSeries $tsTag_velBase}

set velOk 0
if {[info exists gmVelFile] && $gmVelFile ne ""} {
	if {![file exists $gmVelFile]} {
		error "BuildVelSeries.tcl: gmVelFile not found: $gmVelFile"
	}
	set inFd [open $gmVelFile r]
	set line1 [gets $inFd]
	close $inFd

	set isPeer 0
	if {[string match "*PEER*" $line1]} {
		set isPeer 1
	}

	if {$isPeer} {
		set inFd [open $gmVelFile r]
		set hdr [gets $inFd]
		if {![string match "*PEER*" $hdr]} {
			close $inFd
			error "BuildVelSeries.tcl: expected PEER header in $gmVelFile"
		}
		set hdrMeta [gets $inFd]
		set lineUnits [gets $inFd]
		set nPtsLine [gets $inFd]
		if {![regexp -nocase {DT\s*=\s*([0-9]*\.?[0-9]+)} $nPtsLine -> dtStr]} {
			close $inFd
			error "BuildVelSeries.tcl: failed to parse DT from: $nPtsLine"
		}
		set gmVelDT [expr {double($dtStr)}]
		if {$gmVelDT <= 0.0} {
			close $inFd
			error "BuildVelSeries.tcl: DT must be > 0 (got $gmVelDT)"
		}
		if {![regexp -nocase {NPTS\s*=\s*([0-9]+)} $nPtsLine -> nStr]} {
			close $inFd
			error "BuildVelSeries.tcl: failed to parse NPTS from: $nPtsLine"
		}
		set nPts [expr {int($nStr)}]
		if {$nPts < 2} {
			close $inFd
			error "BuildVelSeries.tcl: NPTS must be >= 2 (got $nPts)"
		}
		set rawLine [read $inFd]
		close $inFd

		set unitTok 1.0
		if {[string match "*CM/S*" [string toupper $lineUnits]]} {
			set unitTok 0.01
		} elseif {[string match "*M/S*" [string toupper $lineUnits]] \
			&& ![string match "*CM/S*" [string toupper $lineUnits]]} {
			set unitTok 1.0
		} else {
			error "BuildVelSeries.tcl: unrecognized units line: $lineUnits"
		}

		set velVals {}
		foreach tok [regexp -all -inline -- {[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?} $rawLine] {
			lappend velVals [expr {$unitTok*double($tok)}]
		}
		if {[llength $velVals] != $nPts} {
			error [format "BuildVelSeries.tcl: NPTS=%d but read %d values from %s" \
				$nPts [llength $velVals] $gmVelFile]
		}
		# timeSeries Path $tag -dt $dt -values {...} -factor $gmScaleFactor
		timeSeries Path $tsTag_velBase -dt $gmVelDT -values $velVals -factor $gmScaleFactor
		set gmVelNPTS [llength $velVals]
		set gmVelDuration [expr {($gmVelNPTS - 1)*$gmVelDT}]
		set velOk 1
		set vMax 0.0
		foreach velSamp $velVals {
			set tokA [expr {abs($velSamp)}]
			if {$tokA > $vMax} { set vMax $tokA }
		}
		set vMax [expr {$vMax*$gmScaleFactor}]
		puts [format "----- Vel series  %s  T=%.4g s  |v|nMax=%.4g m/s -----" \
			[file tail $gmVelFile] $gmVelDuration $vMax]
	} else {
		# Filename looks like PEER but header was not
		if {[string match -nocase "*.VT2" $gmVelFile] \
			|| [string match -nocase "*.AT2" $gmVelFile] \
			|| [string match -nocase "*.DT2" $gmVelFile]} {
			error "BuildVelSeries.tcl: $gmVelFile looks like PEER but line 1 is not a PEER header"
		}
		if {![info exists gmVelDT] || $gmVelDT <= 0} {
			error "BuildVelSeries.tcl: gmVelDT must be > 0 when gmVelFile is set"
		}
		timeSeries Path $tsTag_velBase -dt $gmVelDT -filePath $gmVelFile -factor $gmScaleFactor
		set velOk 1
		if {![info exists gmVelNPTS] || $gmVelNPTS < 2} {
			error "BuildVelSeries.tcl: set gmVelNPTS (>=2) for plain -filePath series"
		}
		set gmVelDuration [expr {($gmVelNPTS - 1)*$gmVelDT}]
		puts [format "----- Vel series  %s  T=%.4g s -----" \
			[file tail $gmVelFile] $gmVelDuration]
	}
}

if {!$velOk} {
	if {![info exists gmVelDT] || $gmVelDT <= 0} {
		set gmVelDT 0.005
	}
	timeSeries Path $tsTag_velBase -dt $gmVelDT -values {0.0 0.0} -factor $gmScaleFactor
	set gmVelNPTS 2
	set gmVelDuration $gmVelDT
	puts [format "----- Vel series  dummy zeros  dt=%.4g s -----" $gmVelDT]
}
