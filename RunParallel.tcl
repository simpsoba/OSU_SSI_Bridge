# RunParallel.tcl
# Units: N, m, s
#
# OpenSeesMP driver: same gravity as Run.tcl (replicated on every rank),
# then METIS partition and Mumps + MKRAlphaExplicitMultiSOE for EQ.
# Eigen stays on serial Run.tcl. runEQ 0 stops after partition.
#
#   mpirun -np N OpenSeesMP RunParallel.tcl
#
# Knobs: Parameters.tcl. Switches below.

set runDir [file dirname [file normalize [info script]]]
set root $runDir
set structDir [file join $root structure]
set soilDir [file join $root soil]
set analysisDir [file join $root analysis]
set plotDir [file join $root plot]

set np [getNP]
set pid [getPID]
if {$np < 2} {
	error "RunParallel.tcl: need OpenSeesMP with np>=2 (got np=$np). Serial: OpenSees Run.tcl"
}

# Rank 0 only writes to the terminal. File channels stay on every rank.
if {$pid != 0} {
	rename puts _putsAll
	proc puts {args} {
		set start 0
		if {[lindex $args 0] eq "-nonewline"} {
			set start 1
		}
		set n [expr {[llength $args] - $start}]
		if {$n >= 2} {
			set ch [lindex $args $start]
			if {$ch eq "stdout" || $ch eq "stderr"} {
				return
			}
			_putsAll {*}$args
			return
		}
	}
}

wipe

# ------------------------------------------------------------
# Switches
# ------------------------------------------------------------
# runEQ = 0: gravity, freeze EQ mesh, partition, stop (no transient)
# runEQ = 1: then EQ
# exportPartitionMap = 1: after partition, plot/ExportPartitionMap.tcl
# plotFigures = 1: gravity shape PNG (rank 0)
# eqPrintON = 1: print analysis t, elapsed, pier top every eqPrintDt s (debug)
# eqPrintON = 0: silent EQ loop
# outDIR: recorder folder (EQ + OpenFresco). "" = auto plot/out/...
# gmStartTime: Path -startTime (s). 0 = series starts at t=0
# realTimeON = 0: EQ to eqNstepsAll, recovery on, no OpenFresco
# realTimeON = 1: OpenFresco, no recovery, realTimeNsteps (needs lumpedPlasticity; ignores eqPrintON)
# eleTag_exp: OpenFresco generic element (UX on ZLS-J inner)
set runEQ 1;                              # <-- EDIT  0 | 1
set exportPartitionMap 0;                 # <-- EDIT  0 | 1
set plotFigures 0;                        # <-- EDIT  0 | 1  (rank 0)
set eqPrintON 1;                          # <-- EDIT  0 | 1  (rank that owns pier top)
set eqPrintDt 5.0;                        # <-- EDIT  s, analysis time between progress lines
set outDIR trial1;                        # <-- EDIT  folder ("" = auto)
set gmStartTime 0.0;                      # <-- EDIT  s (0 = omit -startTime)
set realTimeON 0;                         # <-- EDIT  0 | 1
set realTimeNsteps 1000000000000;         # <-- EDIT  steps when realTimeON 1
set eleTag_exp 101;                       # <-- EDIT  OpenFresco generic ele tag

source [file join $root Parameters.tcl]
if {[info exists env(REGEN_PROFILE)] && $env(REGEN_PROFILE) ne ""} {
	set soilProfile $env(REGEN_PROFILE)
}
if {[info exists env(REGEN_BOUNDARY)] && $env(REGEN_BOUNDARY) ne ""} {
	set soilBoundary $env(REGEN_BOUNDARY)
}
if {[info exists env(REGEN_EQ_TMAX)] && $env(REGEN_EQ_TMAX) ne ""} {
	set eqTmax $env(REGEN_EQ_TMAX)
}
if {[info exists env(REGEN_FREE_VIB)] && $env(REGEN_FREE_VIB) ne ""} {
	set eqFreeVibT $env(REGEN_FREE_VIB)
}
if {[info exists env(REGEN_PLOTFIGURES)] && $env(REGEN_PLOTFIGURES) ne ""} {
	set plotFigures $env(REGEN_PLOTFIGURES)
}
if {[info exists env(REGEN_RUNEQ)] && $env(REGEN_RUNEQ) ne ""} {
	set runEQ $env(REGEN_RUNEQ)
}
if {[info exists env(REGEN_EXPORT_PARTITION)] && $env(REGEN_EXPORT_PARTITION) ne ""} {
	set exportPartitionMap $env(REGEN_EXPORT_PARTITION)
}
if {$runEQ != 0 && $runEQ != 1} {
	error "RunParallel.tcl: runEQ must be 0 (stop after partition) or 1 (EQ) (got '$runEQ')"
}
if {$realTimeON != 0 && $realTimeON != 1} {
	error "RunParallel.tcl: realTimeON must be 0 or 1 (got '$realTimeON')"
}
if {$gmStartTime < 0} {
	error "RunParallel.tcl: gmStartTime must be >= 0 (got '$gmStartTime')"
}
if {$realTimeON && $pierEleType ne "lumpedPlasticity"} {
	error "RunParallel.tcl: realTimeON=1 needs pierEleType lumpedPlasticity (got '$pierEleType')"
}
if {$realTimeON && $realTimeNsteps < 1} {
	error "RunParallel.tcl: realTimeNsteps must be >= 1 when realTimeON=1 (got '$realTimeNsteps')"
}
if {$exportPartitionMap != 0 && $exportPartitionMap != 1} {
	error "RunParallel.tcl: exportPartitionMap must be 0 or 1 (got '$exportPartitionMap')"
}
if {$recordersON != 0 && $recordersON != 1 && $recordersON != 2} {
	error "RunParallel.tcl: recordersON must be 0, 1, or 2 (got '$recordersON')"
}

if {$pid == 0} {
	puts [format "RunParallel: np=%d  runEQ=%d  realTimeON=%d  recordersON=%d  exportPartitionMap=%d  pier=%s  pile=%s  profile=%s  boundary=%s  constitutive=%s  springs=%s  soilEle=%s" \
		$np $runEQ $realTimeON $recordersON $exportPartitionMap $pierEleType $pileEleType $soilProfile $soilBoundary $soilConstitutive $pileSpring $soilEleType]
	if {$runEQ} {
		puts [format "  dt=%.6g s  DT_FACTOR=%d  cylinderSF=%.4g  gmStartTime=%.4g s  outDIR=%s" \
			$dtAnalysis $DT_FACTOR $cylinderSF $gmStartTime $outDIR]
	}
}

# ------------------------------------------------------------
# Build (structure nodes + soil quads + SSI springs as needed;
# pier/deck/cap/pile beam-columns after soil gravity)
# ------------------------------------------------------------
source [file join $root BuildModel.tcl]

# Ponding on y=0 if h_water > 0 (no-op if 0).
source [file join $analysisDir WaterSurfaceLoad.tcl]

# ------------------------------------------------------------
# Soil gravity (elastic -> stage 1 -> plastic)
# ------------------------------------------------------------
source [file join $analysisDir GravityHelpers.tcl]
source [file join $analysisDir SoilGravity.tcl]

# ------------------------------------------------------------
# Move structure-node coordinates onto the soil settlement, then create
# pier, deck, cap, and pile beam-columns
# ------------------------------------------------------------
source [file join $analysisDir FoldStructNodes.tcl]

# ------------------------------------------------------------
# Structure weight
# ------------------------------------------------------------
source [file join $analysisDir StructureGravityLoads.tcl]

wipeAnalysis
constraints Transformation
numberer Plain
system UmfPack
test NormDispIncr 5.0e-8 100 0
algorithm KrylovNewton
set dLambda 0.1
set nStep [expr {int(round(1.0/$dLambda))}]
integrator LoadControl $dLambda
analysis Static

set ok [analyze $nStep]
if {$ok != 0} {
	error [format "RunParallel.tcl: structure weight failed rank %d (ok=%d)" $pid $ok]
}
if {$pid == 0} {
	gravPrintSpringKine "after structure weight"
}

set nNodesModel [llength [getNodeTags]]
set nElesModel  [llength [getEleTags]]
set nEqnModel   [systemSize]
if {$pid == 0} {
	puts [format "----- Model size: %d nodes  %d elements  %d DOFs (systemSize) -----" \
		$nNodesModel $nElesModel $nEqnModel]
}

wipeAnalysis
loadConst -time 0.0
set soilGravityDone 1

if {$pid == 0} {
	puts "----- Gravity done  loadConst t=0 -----"
}

if {$plotFigures && $pid == 0} {
	source [file join $plotDir DumpGravityShape.tcl]
	set python3bin [FindPython3]
	if {$python3bin eq ""} {
		puts "RunParallel: WARNING Python not found; gravity JSON only"
	} elseif {[info exists gravityShapeOutPath]} {
		if {[catch {exec {*}$python3bin [file join $root plot PlotGravityShape.py] \
			$gravityShapeOutPath} err]} {
			puts "RunParallel: WARNING PlotGravityShape.py failed:\n$err"
		}
	}
}
barrier

source [file join $analysisDir HoldPierBase.tcl]

# IncrMass $nodeTag $dmx $dmy ?$dIrot?
if {[info procs IncrMass] eq ""} {
	source [file join $structDir IncrMass.tcl]
}
set dmTrans 1.0;                          # kg
set dmRot   0.1;                          # kg.m2
set nMassIncr 0
foreach n [getNodeTags] {
	set ndf [getNDF $n]
	if {$ndf >= 3} {
		IncrMass $n $dmTrans $dmTrans $dmRot
	} else {
		IncrMass $n $dmTrans $dmTrans
	}
	incr nMassIncr
}

# ------------------------------------------------------------
# EQ boundary + partition + transient
# ------------------------------------------------------------
model BasicBuilder -ndm 2 -ndf 2
if {![info exists gmVelNPTS]} {
	source [file join $analysisDir BuildVelSeries.tcl]
}
source [file join $soilDir ActivateEQBoundary.tcl]

# EDIT analysis/RayleighDamping.tcl to change damping model!
source [file join $analysisDir RayleighDamping.tcl]

if {![info exists gmVelDuration]} {
	set gmVelDuration [expr {($gmVelNPTS - 1)*$gmVelDT}]
}
set Trec $gmVelDuration
if {$eqTmax ne "" && $eqTmax > 0 && $eqTmax < $Trec} {
	set Trec $eqTmax
}
set eqNsteps [expr {int(ceil($Trec/$dtAnalysis))}]
if {$eqNsteps < 1} { set eqNsteps 1 }
if {![info exists eqFreeVibT] || $eqFreeVibT eq "" || $eqFreeVibT <= 0} {
	set eqFreeVibT 0.0
	set fvNsteps 0
} else {
	set fvNsteps [expr {int(ceil($eqFreeVibT/$dtAnalysis))}]
	if {$fvNsteps < 1} { set fvNsteps 1 }
}
set eqNstepsAll [expr {$eqNsteps + $fvNsteps}]

wipeAnalysis


# partition splits the frozen EQ mesh (METIS). Gravity stayed replicated.
# realTimeON: pin ZLS-I, eta beam, ZLS-J on rank 0 so node 4 stays local.
if {$realTimeON} {
	# partition -keepOnRank $rank $nEle $ele1 ...
	partition -keepOnRank 0 3 \
		$eleTag_pier_botSpr $eleTag_pier $eleTag_pier_topSpr
} else {
	partition
}
if {$pid == 0} {
	puts [format "----- partition  rank 0/%d  local nodes=%d  eles=%d -----" \
		$np [llength [getNodeTags]] [llength [getEleTags]]]
}
if {$exportPartitionMap} {
	source [file join $plotDir ExportPartitionMap.tcl]
}
if {!$runEQ} {
	if {$pid == 0} {
		puts "RunParallel: gravity + partition done (runEQ=0, no EQ)"
	}
} else {
	# ---
	# OpenFresco (realTimeON): rank 0 only, after keepOnRank, before numberer
	# ---
	if {$realTimeON} {
		if {$pid == 0} {
			if {[lsearch -exact [getNodeTags] $nodeTag_pierTopZeroLengthInner] < 0} {
				error "RunParallel.tcl: realTimeON node $nodeTag_pierTopZeroLengthInner missing on rank 0 after keepOnRank"
			}
			# model BasicBuilder -ndm $ndm -ndf $ndf
			model BasicBuilder -ndm 2 -ndf 3
			loadPackage OpenFrescoTcl
			puts "\n-------------------"
			puts "experimental element on"

			set Kexp 1e-2;                  # N/m  -initStif (1 dof)

			# expControlPoint $tag $dof rsp
			expControlPoint 1 1 disp
			expControlPoint 2 1 force

			# expControl SCRAMNetGT $tag -nodeID $id $memSize -trialCP $cp -outCP $cp
			expControl SCRAMNetGT 1 -nodeID 3 4096 -trialCP 1 -outCP 2

			# expSetup NoTransformation $tag -control $ctrlTag -dir $dof -sizeTrialOut $nTrial $nOut
			expSetup NoTransformation 1 -control 1 \
				-dir 1 \
				-sizeTrialOut 1 1

			# expSite LocalSite $tag $setupTag
			expSite LocalSite 1 1

			# expElement generic $eleTag -node $nodeTag -dof $dof -site $site -initStif $K
			# node 4 = $nodeTag_pierTopZeroLengthInner: top node of the stiff element, UX
			expElement generic $eleTag_exp -node $nodeTag_pierTopZeroLengthInner \
				-dof 1 -site 1 -initStif $Kexp -noRayleigh -checkTime
		}
		barrier
	}

	# numberer $type
	numberer ParallelPlain
	# numberer ParallelRCM
	# system $type
	# system DistributedCuDSS
	# system ParallelProfileSPD
	system Mumps
	# constraints $type
	# Note: because of the sp-roller condition, Plain constraint handler cannot be used for the EQ Analysis
	# off the relative displacement of its own two nodes.
	# constraints Plain
	constraints Auto
	# constraints Penalty 1.0e18 1.0e18
	# test $type $tol $maxIter $flag
	test EnergyIncr 1e-8 25 0
	# algorithm $type
	algorithm Linear
	# integrator MKRAlphaExplicitMultiSOE $rhoInf
	# integrator MKRAlphaExplicitMultiSOE 0.5
	integrator MKRAlphaExplicitMultiSOE 0.5 -incrementalAccel
	# integrator AlphaOSGeneralized 0.0
	# integrator CudaMKRAlpha 0.5
	# analysis $type
	analysis Transient

	source [file join $analysisDir EQRecorders.tcl]

	if {$realTimeON} {
		if {$pid == 0} {
			if {![file isdirectory $eqOutDir]} {
				file mkdir $eqOutDir
			}
			# recorder Element -file $fileName -time -ele $eleTag $rsp
			recorder Element -file [file join $eqOutDir Elmt${eleTag_exp}_ctrlDsp.out] \
				-time -ele $eleTag_exp ctrlDisp
			recorder Element -file [file join $eqOutDir Elmt${eleTag_exp}_GlbFrc.out] \
				-time -ele $eleTag_exp forces
			recorder Element -file [file join $eqOutDir Elmt${eleTag_exp}_BscFrc.out] \
				-time -ele $eleTag_exp basicForces
			# expRecorder Setup -file $fileName -time -setup $tag $rsp
			expRecorder Setup -file [file join $eqOutDir ServerSetup_daqFrc.out] \
				-time -setup 1 daqForce
		}
		barrier
	}
	record
	
	if {$pid == 0} {
		if {$realTimeON} {
			puts [format "----- EQ (OpenSeesMP + Mumps)  realTimeON  dt=%.6g s  nSteps=%s  rec=%d -----" \
				$dtAnalysis $realTimeNsteps $recordersON]
		} else {
			puts [format "----- EQ (OpenSeesMP + Mumps)  dt=%.6g s  T=%.4g+%.4g s  nSteps=%d  rec=%d -----" \
				$dtAnalysis $Trec $eqFreeVibT $eqNsteps $recordersON]
		}
	}

	set ownsPierTop 0
	if {[lsearch -exact [getNodeTags] $nodeTag_pierTop_deckBC] >= 0} {
		set ownsPierTop 1
	}
	if {$pid == 0 && !$ownsPierTop} {
		puts "  pier top not on rank 0 -- progress is t/step only"
	}

	set t0 [clock microseconds]
	set ok 0
	set nFail 0
	if {$realTimeON} {
		for {set i 1} {$i <= $realTimeNsteps} {incr i} {
			analyze 1 $dtAnalysis
		}
	} else {
		set dtHalf [expr {0.5*$dtAnalysis}]
		set dtQ [expr {0.25*$dtAnalysis}]
		set dtPrint $eqPrintDt
		set tPrint $dtPrint
		for {set i 1} {$i <= $eqNstepsAll} {incr i} {
			if {$fvNsteps > 0 && $i == [expr {$eqNsteps + 1}] && $pid == 0} {
				puts [format "----- free vibration %.4g s (%d steps) -----" \
					$eqFreeVibT $fvNsteps]
			}
			set ok [analyze 1 $dtAnalysis]
			if {$ok != 0} {
				incr nFail
				if {$pid == 0} {
					puts [format "  recover at step %d  t~%.4g s" $i [getTime]]
				}
				set ok [analyze 1 $dtHalf]
				if {$ok == 0} {
					set ok [analyze 1 $dtHalf]
				}
				if {$ok != 0} {
					if {$pid == 0} {
						puts [format "  recover dt=%.6g s at step %d  t~%.4g s" $dtQ $i [getTime]]
					}
					set ok 0
					for {set k 1} {$k <= 4} {incr k} {
						set ok [analyze 1 $dtQ]
						if {$ok != 0} { break }
					}
				}
				if {$ok != 0} {
					error [format "RunParallel.tcl: analyze failed at step %d / %d (t~%.4g s) after dt/2, dt/4x4" \
						$i $eqNstepsAll [getTime]]
				}
			}
			if {$eqPrintON && $pid == 0} {
				set tNow [getTime]
				set elapsed [expr {([clock microseconds] - $t0)/1.0e6}]
				while {$tNow >= $tPrint} {
					if {$ownsPierTop} {
						puts [format "  analysis t=%.3f s  elapsed=%.2f s  pier top ux=%.4e m  uy=%.4e m  step %d / %d" \
							$tNow $elapsed \
							[nodeDisp $nodeTag_pierTop_deckBC 1] [nodeDisp $nodeTag_pierTop_deckBC 2] \
							$i $eqNstepsAll]
					} else {
						puts [format "  analysis t=%.3f s  elapsed=%.2f s  step %d / %d" \
							$tNow $elapsed $i $eqNstepsAll]
					}
					set tPrint [expr {$tPrint + $dtPrint}]
				}
			}
		}
	}
	set elapsed [expr {([clock microseconds] - $t0)/1.0e6}]
	set soilEQDone 1

	if {$pid == 0} {
		if {$realTimeON} {
			puts [format "----- EQ done  realTimeON  nSteps=%s -----" $realTimeNsteps]
		} else {
			puts [format "----- EQ done  recoveries=%d -----" $nFail]
		}
		puts [format "  analysis t=%.4g s  elapsed=%.2f s  (EQ %.4g s + freeVib %.4g s)" \
			[getTime] $elapsed $Trec $eqFreeVibT]
	}
	if {$pid == 0 && $ownsPierTop} {
		puts [format "  pier top ux=%.4e m  uy=%.4e m" \
			[nodeDisp $nodeTag_pierTop_deckBC 1] [nodeDisp $nodeTag_pierTop_deckBC 2]]
	}
	if {$pid == 0} {
		if {$recordersON == 0 && !$realTimeON} {
			puts "  recordersON=0"
		} else {
			puts [format "  outDIR=%s  (per-rank files *.\$pid)" $eqOutDir]
		}
		puts "RunParallel: gravity + EQ done"
	}
}
