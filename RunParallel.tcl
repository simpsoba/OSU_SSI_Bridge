# RunParallel.tcl
# Units: N, m, s
#
# OpenSeesMP driver: same gravity as Run.tcl (replicated on every rank),
# then METIS partition and parallel SOE + EQ integrator.
# Eigen stays on serial Run.tcl. runEQ 0 stops after partition.
#
#   mpirun -np N OpenSeesMP RunParallel.tcl
#
# Knobs: Parameters.tcl. Switches + analysis knobs below.

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
	# logFile $fileName ?-noEcho? ?-append?
	# opserr (C++ warnings/errors) to this file only; Tcl puts hijacked below.
	logFile [file join $runDir [format "opensees.rank%d.log" $pid]] -noEcho
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
# gmStartTime: Path -startTime (s). 0 = series starts at t=0.
#   When realTimeON is 0, eqNstepsAll covers startTime + Trec + eqFreeVibT.
# realTimeON = 0: EQ to eqNstepsAll, recovery on, no OpenFresco
# realTimeON = 1: OpenFresco, no recovery, realTimeNsteps (ignores eqPrintON)
# expElementType: OpenFresco element on the pier beam ends (UX)
#   "generic"     -- single-node force at the beam top (current)
#   "twoNodeLink" -- equal-opposite UX forces between beam i/j nodes
#     lumpedPlasticity: nodes 2--4; elastic/forceBeam: nodes 1--5
# eleTag_exp: experimental element tag
set runEQ 1;                              # <-- EDIT  0 | 1
set exportPartitionMap 0;                 # <-- EDIT  0 | 1
set plotFigures 0;                        # <-- EDIT  0 | 1  (rank 0)
set eqPrintON 1;                          # <-- EDIT  0 | 1  (rank that owns pier top)
set eqPrintDt 5.0;                        # <-- EDIT  s, analysis time between progress lines
set outDIR trial1;                        # <-- EDIT  folder ("" = auto)
set gmStartTime 0.0;                      # <-- EDIT  s (0 = omit -startTime)
set realTimeON 0;                         # <-- EDIT  0 | 1
set realTimeNsteps 1000000000000;         # <-- EDIT  steps when realTimeON 1
set expElementType "generic";             # <-- EDIT  generic | twoNodeLink
set eleTag_exp 101;                       # <-- EDIT  OpenFresco ele tag

# ------------------------------------------------------------
# Analysis knobs (EDIT these strings — include any args)
# ------------------------------------------------------------
# Gravity + structure weight BEFORE partition (full mesh on every rank).
#   "UmfPack"             recommended (each rank solves its own copy)
#   "CuDSS"               GPU; each rank may open the device
#   "BandGeneral" | "ProfileSPD"   (numberer RCM)
#   "Mumps"               allowed; see NOTE
#   "DistributedCuDSS"    allowed; see NOTE (only rank 0 uses the GPU)
#   "ParallelProfileSPD"  allowed; see NOTE (numberer ParallelRCM)
#
# NOTE — Mumps / DistributedCuDSS / ParallelProfileSPD before partition:
#   Every rank still has the full mesh. The parallel solver adds every rank's
#   K and F, so you solve (np*K) x = np*F. Same x, ~np times the work.
#   DistributedCuDSS is mainly so only rank 0 talks to the GPU under MPI.
#
# EQ AFTER partition (mesh split):
#   "Mumps" | "DistributedCuDSS" | "ParallelProfileSPD"
#
set prePartitionSystem  "UmfPack";                              # <-- EDIT
set postPartitionSystem "DistributedCuDSS";                                # <-- EDIT

# EQ constraints (gravity always uses Transformation).
# ASDEA forces Transformation for EQ no matter what you put here.
#   "Auto" | "Transformation" | "Plain" | "Penalty 1.0e18 1.0e18"
set constraintsHandler  "Transformation";                                 # <-- EDIT

# EQ integrator string (name + args). Algorithm/test follow from the name.
#   "MKRAlphaExplicitMultiSOE 0.5 -incrementalAccel"  -> Linear, no test
#   "MKRAlphaExplicitMultiSOE 0.5"
#   "CudaMKRAlpha 0.5 -incrementalAccel"              -> forces DistributedCuDSS
#   "CudaMKRAlpha 0.5"
#   "AlphaOSGeneralized 0.5"
#   "TRBDF2"                                          -> KrylovNewton + test
#   "Newmark 0.5 0.25"
set eqIntegrator        "CudaMKRAlpha 0.5";  # <-- EDIT

# --- apply system + matching numberer (used for pre and post) ---
proc applySystem {sysStr} {
	set name [lindex $sysStr 0]
	if {$name eq "Mumps" || $name eq "DistributedCuDSS"} {
		# Before partition (replicated mesh): works via (np*K)x = np*F; see NOTE above.
		numberer ParallelPlain
		system {*}$sysStr
	} elseif {$name eq "ParallelProfileSPD"} {
		numberer ParallelRCM
		system {*}$sysStr
	} elseif {$name eq "BandGeneral" || $name eq "ProfileSPD"} {
		numberer RCM
		system {*}$sysStr
	} elseif {$name eq "UmfPack" || $name eq "CuDSS"} {
		numberer Plain
		system {*}$sysStr
	} else {
		error "applySystem: unknown system '$sysStr'"
	}
}

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
if {$expElementType ne "generic" && $expElementType ne "twoNodeLink"} {
	error "RunParallel.tcl: expElementType must be generic or twoNodeLink (got '$expElementType')"
}
if {$gmStartTime < 0} {
	error "RunParallel.tcl: gmStartTime must be >= 0 (got '$gmStartTime')"
}
if {$realTimeON && $realTimeNsteps < 1} {
	error "RunParallel.tcl: realTimeNsteps must be >= 1 when realTimeON=1 (got '$realTimeNsteps')"
}
if {$exportPartitionMap != 0 && $exportPartitionMap != 1} {
	error "RunParallel.tcl: exportPartitionMap must be 0 or 1 (got '$exportPartitionMap')"
}
if {![string is integer -strict $recordersON] || $recordersON < 0 || $recordersON > 3} {
	error "RunParallel.tcl: recordersON must be an integer 0..3 (got '$recordersON')"
}

if {$pid == 0} {
	puts [format "RunParallel: np=%d  runEQ=%d  realTimeON=%d  expElementType=%s  recordersON=%d  exportPartitionMap=%d  pier=%s  pile=%s  profile=%s  boundary=%s  constitutive=%s  springs=%s  soilEle=%s  soilMesh=%d" \
		$np $runEQ $realTimeON $expElementType $recordersON $exportPartitionMap $pierEleType $pileEleType $soilProfile $soilBoundary $soilConstitutive $pileSpring $soilEleType $soilMesh]
	puts [format "  prePartitionSystem=%s  postPartitionSystem=%s  constraints=%s  integrator=%s" \
		$prePartitionSystem $postPartitionSystem $constraintsHandler $eqIntegrator]
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
applySystem $prePartitionSystem
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

# Hold pier-base UX/UY at gravity disp (holdPierON; see HoldPierBase.tcl).
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
# Path is silent until t >= gmStartTime; fold that pad into eqNsteps
set tWait 0.0
if {[info exists gmStartTime] && $gmStartTime ne "" && $gmStartTime > 0} {
	set tWait $gmStartTime
}
set eqNsteps [expr {int(ceil(($tWait + $Trec)/$dtAnalysis))}]
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
# Profiles 1–2: -samePart keeps each SSI spring with continuum at its soil node
# (needed for PyLiq1/TzLiq1 and sand-column SSI under MPI).
# realTimeON: keep the pier beam (and ZLS hinges) on rank 0 for OpenFresco.
set partitionArgs {}
if {($soilProfile == 1 || $soilProfile == 2) \
		&& [info exists ssiPartitionSamePart] \
		&& [llength $ssiPartitionSamePart] > 0} {
	foreach grp $ssiPartitionSamePart {
		set nSame [llength $grp]
		if {$nSame >= 2} {
			# partition -samePart $n $e1 $e2 ...
			lappend partitionArgs -samePart $nSame {*}$grp
		}
	}
	if {$pid == 0} {
		puts [format "----- partition -samePart  %d spring/continuum groups (profile %d) -----" \
			[llength $ssiPartitionSamePart] $soilProfile]
	}
}
if {$realTimeON} {
	# partition ... -keepOnRank $rank $nEle $ele1 ...
	if {$pierEleType eq "lumpedPlasticity"} {
		partition {*}$partitionArgs -keepOnRank 0 3 \
			$eleTag_pier_botSpr $eleTag_pier $eleTag_pier_topSpr
	} else {
		partition {*}$partitionArgs -keepOnRank 0 1 $eleTag_pier
	}
} else {
	partition {*}$partitionArgs
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
	#   generic / twoNodeLink on ends of eleTag_pier (see Run.tcl comments)
	# ---
	if {$realTimeON} {
		if {$pierEleType eq "lumpedPlasticity"} {
			set expNodeI $nodeTag_pierBaseZeroLengthInner
			set expNodeJ $nodeTag_pierTopZeroLengthInner
		} else {
			set expNodeI $nodeTag_pierBase_capTC
			set expNodeJ $nodeTag_pierTop_deckBC
		}
		if {$pid == 0} {
			foreach nCheck [list $expNodeI $expNodeJ] {
				if {[lsearch -exact [getNodeTags] $nCheck] < 0} {
					error "RunParallel.tcl: realTimeON node $nCheck missing on rank 0 after keepOnRank"
				}
			}
			# model BasicBuilder -ndm $ndm -ndf $ndf
			model BasicBuilder -ndm 2 -ndf 3
			loadPackage OpenFrescoTcl
			puts "\n-------------------"
			puts "experimental element on ($expElementType)"

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

			if {$expElementType eq "twoNodeLink"} {
				expElement twoNodeLink $eleTag_exp \
					$expNodeI $expNodeJ \
					-dir 1 -site 1 -initStif $Kexp -noRayleigh \
					-orient 1.0 0.0 0.0  0.0 1.0 0.0
				puts [format "  twoNodeLink ele %d  nodes %d--%d (UX)" \
					$eleTag_exp $expNodeI $expNodeJ]
			} else {
				expElement generic $eleTag_exp -node $expNodeJ \
					-dof 1 -site 1 -initStif $Kexp -noRayleigh -checkTime
				puts [format "  generic ele %d  node %d (UX)" $eleTag_exp $expNodeJ]
			}
		}
		barrier
	}
	
	# --- EQ analysis objects (knobs: postPartitionSystem, constraintsHandler, eqIntegrator) ---
	set intName [lindex $eqIntegrator 0]

	# system of equations + solver (+ matching numberer inside applySystem).
	# CudaMKRAlpha needs DistributedCuDSS so only rank 0 uses the GPU under MPI.
	if {$intName eq "CudaMKRAlpha" && [lindex $postPartitionSystem 0] ne "DistributedCuDSS"} {
		puts "CudaMKRAlpha -> forcing postPartitionSystem DistributedCuDSS (was $postPartitionSystem)"
		set postPartitionSystem "DistributedCuDSS"
	}
	applySystem $postPartitionSystem

	# constraint handler (ASDEA sp-rollers need Transformation)
	if {[info exists soilBoundary] && $soilBoundary eq "ASDEA"} {
		constraints Transformation
	} else {
		constraints {*}$constraintsHandler
	}

	# time integrator
	integrator {*}$eqIntegrator

	# solution algorithm (+ convergence test when the integrator is implicit)
	if {$intName eq "MKRAlphaExplicitMultiSOE" || $intName eq "CudaMKRAlpha" \
			|| $intName eq "AlphaOSGeneralized"} {
		algorithm Linear
	} elseif {$intName eq "TRBDF2" || $intName eq "Newmark"} {
		test NormDispIncr 1.0e-8 25 0
		algorithm KrylovNewton
	} else {
		error "RunParallel.tcl: unknown eqIntegrator '$eqIntegrator'"
	}

	# analysis type
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
			puts [format "----- EQ (OpenSeesMP + %s + %s)  realTimeON  dt=%.6g s  nSteps=%s  rec=%d -----" \
				$postPartitionSystem $eqIntegrator $dtAnalysis $realTimeNsteps $recordersON]
		} else {
			puts [format "----- EQ (OpenSeesMP + %s + %s)  dt=%.6g s  T=%.4g+%.4g+%.4g s (start+EQ+freeVib)  nSteps=%d  rec=%d -----" \
				$postPartitionSystem $eqIntegrator $dtAnalysis $tWait $Trec $eqFreeVibT $eqNstepsAll $recordersON]
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
				# test $type $tol $maxIter $flag
				test NormDispIncr 1.0e-6 25 0
				if {$pid == 0} {
					puts [format "  recover NormDispIncr 1e-6 at step %d  t~%.4g s" $i [getTime]]
				}
				set ok [analyze 1 $dtAnalysis]
				if {$ok != 0} {
					set ok [analyze 1 $dtHalf]
					if {$ok == 0} {
						set ok [analyze 1 $dtHalf]
					}
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
				# test $type $tol $maxIter $flag
				test NormDispIncr 1.0e-8 25 0
				if {$ok != 0} {
					error [format "RunParallel.tcl: analyze failed at step %d / %d (t~%.4g s) after NormDispIncr 1e-8/1e-6, dt/2, dt/4x4" \
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
		puts [format "  analysis t=%.4g s  elapsed=%.2f s  (start %.4g s + EQ %.4g s + freeVib %.4g s)" \
			[getTime] $elapsed $tWait $Trec $eqFreeVibT]
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
