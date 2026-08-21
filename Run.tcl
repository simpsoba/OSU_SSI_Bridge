# Run.tcl
# Units: N, m, s
#
# Driver: parameters -> build (structure nodes, soil quads, SSI springs as needed)
# -> soil gravity -> move structure nodes onto settlement ->
# pier/deck/cap/pile beam-columns ->
# structure weight -> hold pier-base UX/UY -> (eigen | EQ).
# Papers: Shin et al. (2007); Mackie et al. (2008); Kramer PEER 2008/07;
# Neumann (2021); Neumann et al. (2023). See reference/ and NOTES.md.
#
# Knobs: Parameters.tcl (TAGS CONVENTION for IDs). Switches + analysis knobs below.
# Optional: OpenSees Run.tcl Overrides.tcl  (see overridesON; RunTestMatrix.py --row N)

# Folders next to this file (structure/, soil/, analysis/, plot/).
set runDir [file dirname [file normalize [info script]]]
set root $runDir
set structDir [file join $root structure]
set soilDir [file join $root soil]
set analysisDir [file join $root analysis]
set plotDir [file join $root plot]

wipe

# ------------------------------------------------------------
# Switches
# ------------------------------------------------------------
# runEQ = 0: gravity + eigen (absorbing BCs still in gravity state)
# runEQ = 1: gravity then EQ (Lysmer / ASDEA stage 1)
# plotFigures = 0: analysis only
# plotFigures = 1: gravity shape PNG; modes too if runEQ 0
# eqPrintON = 0: silent EQ loop (no per-interval wall-clock timings)
# eqPrintON = 1: print analysis t, elapsed, pier top every eqPrintDt s
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
# overridesON: apply argv Overrides.tcl (forced to 0 if no file is passed)
set runEQ 1;                              # <-- EDIT  0 | 1
set plotFigures 0;                        # <-- EDIT  0 | 1
set eqPrintON 1;                          # <-- EDIT  0 | 1
set eqPrintDt 5.0;                        # <-- EDIT  s, analysis time between progress lines
set outDIR trial1;                        # <-- EDIT  folder ("" = auto)
set gmStartTime 0.0;                      # <-- EDIT  s (0 = omit -startTime)
set realTimeON 0;                         # <-- EDIT  0 | 1
set realTimeNsteps 1000000000000;         # <-- EDIT  steps when realTimeON 1
set expElementType "generic";             # <-- EDIT  generic | twoNodeLink
set eleTag_exp 101;                       # <-- EDIT  OpenFresco ele tag
set overridesON 1;                        # <-- EDIT  0 | 1  (forced to 0 if no argv file)

# ------------------------------------------------------------
# Analysis knobs (EDIT these strings — include any args)
# ------------------------------------------------------------
# Gravity + structure weight (no MPI partition here).
#   "UmfPack" | "CuDSS" | "FullGeneral"
#   "BandGeneral" | "ProfileSPD"   (numberer RCM)
# (DistributedCuDSS / Mumps need OpenSeesMP — use RunParallel.tcl.)
#
set prePartitionSystem  "UmfPack";                              # <-- EDIT
set postPartitionSystem "UmfPack";                              # <-- EDIT

# EQ constraints (gravity always uses Transformation).
# ASDEA forces Transformation for EQ no matter what you put here.
#   "Auto" | "Transformation" | "Plain" | "Penalty 1.0e18 1.0e18"
set constraintsHandler  "Auto";                                 # <-- EDIT

# EQ integrator string (name + args). Algorithm/test follow from the name.
#   "MKRAlphaExplicitMultiSOE 0.5 -incrementalAccel"  -> Linear, no test
#   "MKRAlphaExplicitMultiSOE 0.5"
#   "CudaMKRAlpha 0.5 -incrementalAccel"              -> forces CuDSS
#   "CudaMKRAlpha 0.5"
#   "AlphaOSGeneralized 0.5"
#   "TRBDF2"                                          -> KrylovNewton + test
#   "Newmark 0.5 0.25"
set eqIntegrator        "MKRAlphaExplicitMultiSOE 0.5 -incrementalAccel";  # <-- EDIT

# --- apply system + matching numberer (used for pre and post) ---
proc applySystem {sysStr} {
	set name [lindex $sysStr 0]
	if {$name eq "BandGeneral" || $name eq "ProfileSPD"} {
		numberer RCM
		system {*}$sysStr
	} elseif {$name eq "UmfPack" || $name eq "CuDSS" || $name eq "FullGeneral"} {
		numberer Plain
		system {*}$sysStr
	} else {
		error "applySystem: unknown system '$sysStr' (UmfPack, CuDSS, BandGeneral, FullGeneral, ProfileSPD)"
	}
}

source [file join $root Parameters.tcl]
if {[info exists env(REGEN_PROFILE)] && $env(REGEN_PROFILE) ne ""} {
	set soilProfile $env(REGEN_PROFILE)
}
if {[info exists env(REGEN_BOUNDARY)] && $env(REGEN_BOUNDARY) ne ""} {
	set soilBoundary $env(REGEN_BOUNDARY)
}
# Same short-run overrides RunParallel.tcl takes (smoke tests without editing
# Parameters.tcl).
if {[info exists env(REGEN_EQ_TMAX)] && $env(REGEN_EQ_TMAX) ne ""} {
	set eqTmax $env(REGEN_EQ_TMAX)
}
if {[info exists env(REGEN_FREE_VIB)] && $env(REGEN_FREE_VIB) ne ""} {
	set eqFreeVibT $env(REGEN_FREE_VIB)
}

# Optional overrides file: OpenSees Run.tcl Overrides.tcl
set overridesFile ""
if {[info exists argv] && [llength $argv] >= 1} {
	set overridesFile [lindex $argv 0]
}
if {$overridesFile eq ""} {
	set overridesON 0
}
if {$overridesON && $overridesFile ne ""} {
	if {![file exists $overridesFile]} {
		error "Run.tcl: overrides file not found: $overridesFile"
	}
	source $overridesFile
	RefreshDerivedKnobs
	puts "----- Overrides ON  $overridesFile -----"
} elseif {$overridesFile ne "" && !$overridesON} {
	puts "----- Overrides file ignored (overridesON=0)  $overridesFile -----"
}

if {$runEQ != 0 && $runEQ != 1} {
	error "Run.tcl: runEQ must be 0 (gravity) or 1 (gravity+EQ) (got '$runEQ')"
}
if {$realTimeON != 0 && $realTimeON != 1} {
	error "Run.tcl: realTimeON must be 0 or 1 (got '$realTimeON')"
}
if {$expElementType ne "generic" && $expElementType ne "twoNodeLink"} {
	error "Run.tcl: expElementType must be generic or twoNodeLink (got '$expElementType')"
}
if {$gmStartTime < 0} {
	error "Run.tcl: gmStartTime must be >= 0 (got '$gmStartTime')"
}
if {$realTimeON && $realTimeNsteps < 1} {
	error "Run.tcl: realTimeNsteps must be >= 1 when realTimeON=1 (got '$realTimeNsteps')"
}

if {![string is integer -strict $recordersON] || $recordersON < 0 || $recordersON > 4} {
	error "Run.tcl: recordersON must be an integer 0..4 (got '$recordersON')"
}

puts [format "Run: runEQ=%d  realTimeON=%d  expElementType=%s  recordersON=%d  pier=%s  pile=%s  profile=%s  boundary=%s  constitutive=%s  springs=%s  soilEle=%s  soilMesh=%d" \
	$runEQ $realTimeON $expElementType $recordersON $pierEleType $pileEleType $soilProfile $soilBoundary $soilConstitutive $pileSpring $soilEleType $soilMesh]
puts [format "  prePartitionSystem=%s  postPartitionSystem=%s  constraints=%s  integrator=%s" \
	$prePartitionSystem $postPartitionSystem $constraintsHandler $eqIntegrator]
if {$runEQ} {
	puts [format "  dt=%.6g s  DT_FACTOR=%d  cylinderSF=%.4g  gmStartTime=%.4g s  outDIR=%s" \
		$dtAnalysis $DT_FACTOR $cylinderSF $gmStartTime $outDIR]
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
# algorithm Newton
set dLambda 0.1
set nStep [expr {int(round(1.0/$dLambda))}]
integrator LoadControl $dLambda
analysis Static

set ok [analyze $nStep]
if {$ok != 0} {
	error "Run.tcl: structure weight failed (ok=$ok)"
}
gravPrintSpringKine "after structure weight"

set nNodesModel [llength [getNodeTags]]
set nElesModel  [llength [getEleTags]]
set nEqnModel   [systemSize]
puts [format "----- Model size: %d nodes  %d elements  %d DOFs (systemSize) -----" \
	$nNodesModel $nElesModel $nEqnModel]

wipeAnalysis
loadConst -time 0.0
set soilGravityDone 1

puts "----- Gravity done  loadConst t=0 -----"

# ------------------------------------------------------------
# Figures (optional)
# ------------------------------------------------------------
if {$plotFigures} {
	source [file join $plotDir DumpGravityShape.tcl]
	set python3bin [FindPython3]
	if {$python3bin eq ""} {
		puts "Run: WARNING Python not found; gravity JSON only"
	} elseif {[info exists gravityShapeOutPath]} {
		if {[catch {exec {*}$python3bin [file join $root plot PlotGravityShape.py] \
			$gravityShapeOutPath} err]} {
			puts "Run: WARNING PlotGravityShape.py failed:\n$err"
		}
	}
}

# Hold pier-base UX/UY at gravity disp (holdPierON; see HoldPierBase.tcl).
source [file join $analysisDir HoldPierBase.tcl]

# IncrMass $nodeTag $dmx $dmy ?$dIrot?
# Small lump on every node so M is not singular. Element rho stays.
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

if {!$runEQ} {
	# ------------------------------------------------------------
	# Eigen (gravity BC; before dashpots / ASDEA stage 1)
	# ------------------------------------------------------------
	wipeAnalysis
	constraints Transformation
	applySystem $prePartitionSystem

	puts [format "----- Eigen (%d modes) -----" $nModesEigen]
	set lambdas [eigen $nModesEigen]

	set i 1
	puts "| mode | lambda | omega (rad/s) | T (s) | f (Hz) |"
	foreach lam $lambdas {
		if {$lam <= 0.0} {
			puts [format "| %4d | %10.4e | -- | -- | -- |  (non-positive)" $i $lam]
		} else {
			set w [expr {sqrt($lam)}]
			set T [expr {2.0*$pi/$w}]
			set f [expr {1.0/$T}]
			puts [format "| %4d | %10.4e | %12.5f | %8.5f | %8.5f |" $i $lam $w $T $f]
		}
		incr i
	}
	set eigenLambdas $lambdas
	set soilEigenDone 1
	puts "----- Eigen done -----"

	source [file join $analysisDir EigenAfterGravity.tcl]
	puts "Run: gravity + eigen done"
} else {
	# ------------------------------------------------------------
	# EQ boundary + transient
	# ------------------------------------------------------------
	model BasicBuilder -ndm 2 -ndf 2
	if {![info exists gmVelNPTS]} {
		source [file join $analysisDir BuildVelSeries.tcl]
	}
	source [file join $soilDir ActivateEQBoundary.tcl]
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

	# ---
	# OpenFresco (realTimeON): expElement on pier beam ends (UX)
	#   generic: force at beam top; twoNodeLink: equal-opposite on i/j
	#   attachment = ends of eleTag_pier (2--4 lumped; 1--5 otherwise)
	# ---
	if {$realTimeON} {
		# model BasicBuilder -ndm $ndm -ndf $ndf
		model BasicBuilder -ndm 2 -ndf 3
		loadPackage OpenFrescoTcl
		puts "\n-------------------"
		puts "experimental element on ($expElementType)"

		set Kexp 1e-2;                      # N/m  -initStif (1 dof)

		# Ends of the numerical pier beam (eleTag_pier).
		if {$pierEleType eq "lumpedPlasticity"} {
			set expNodeI $nodeTag_pierBaseZeroLengthInner
			set expNodeJ $nodeTag_pierTopZeroLengthInner
		} else {
			set expNodeI $nodeTag_pierBase_capTC
			set expNodeJ $nodeTag_pierTop_deckBC
		}

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
			# Parallel UX actuator: equal-opposite forces on beam i/j.
			expElement twoNodeLink $eleTag_exp \
				$expNodeI $expNodeJ \
				-dir 1 -site 1 -initStif $Kexp -noRayleigh \
				-orient 1.0 0.0 0.0  0.0 1.0 0.0
			puts [format "  twoNodeLink ele %d  nodes %d--%d (UX)" \
				$eleTag_exp $expNodeI $expNodeJ]
		} else {
			# Single-node experimental force at beam top, UX.
			expElement generic $eleTag_exp -node $expNodeJ \
				-dof 1 -site 1 -initStif $Kexp -noRayleigh -checkTime
			puts [format "  generic ele %d  node %d (UX)" $eleTag_exp $expNodeJ]
		}
	}
	
	wipeAnalysis
	# --- EQ analysis objects (knobs: postPartitionSystem, constraintsHandler, eqIntegrator) ---
	set intName [lindex $eqIntegrator 0]

	# system of equations + solver (+ matching numberer inside applySystem).
	# CudaMKRAlpha needs CuDSS.
	if {$intName eq "CudaMKRAlpha" && [lindex $postPartitionSystem 0] ne "CuDSS"} {
		puts "CudaMKRAlpha -> forcing postPartitionSystem CuDSS (was $postPartitionSystem)"
		set postPartitionSystem "CuDSS"
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
		error "Run.tcl: unknown eqIntegrator '$eqIntegrator'"
	}

	# analysis type
	analysis Transient

	source [file join $analysisDir EQRecorders.tcl]

	if {$realTimeON} {
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
	record
	
	
	if {$realTimeON} {
		puts [format "----- EQ (%s + %s)  realTimeON  dt=%.6g s  nSteps=%s  rec=%d -----" \
			$postPartitionSystem $eqIntegrator $dtAnalysis $realTimeNsteps $recordersON]
	} else {
		puts [format "----- EQ (%s + %s)  dt=%.6g s  T=%.4g+%.4g+%.4g s (start+EQ+freeVib)  nSteps=%d  rec=%d -----" \
			$postPartitionSystem $eqIntegrator $dtAnalysis $tWait $Trec $eqFreeVibT $eqNstepsAll $recordersON]
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
			if {$fvNsteps > 0 && $i == [expr {$eqNsteps + 1}]} {
				puts [format "----- free vibration %.4g s (%d steps) -----" \
					$eqFreeVibT $fvNsteps]
			}
			set ok [analyze 1 $dtAnalysis]
			if {$ok != 0} {
				incr nFail
				puts [format "  recover at step %d  t~%.4g s" $i [getTime]]
				# test $type $tol $maxIter $flag
				test NormDispIncr 1.0e-6 25 0
				puts [format "  recover NormDispIncr 1e-6 at step %d  t~%.4g s" $i [getTime]]
				set ok [analyze 1 $dtAnalysis]
				if {$ok != 0} {
					set ok [analyze 1 $dtHalf]
					if {$ok == 0} {
						set ok [analyze 1 $dtHalf]
					}
				}
				if {$ok != 0} {
					puts [format "  recover dt=%.6g s at step %d  t~%.4g s" $dtQ $i [getTime]]
					set ok 0
					for {set k 1} {$k <= 4} {incr k} {
						set ok [analyze 1 $dtQ]
						if {$ok != 0} { break }
					}
				}
				# test $type $tol $maxIter $flag
				test NormDispIncr 1.0e-8 25 0
				if {$ok != 0} {
					error [format "Run.tcl: analyze failed at step %d / %d (t~%.4g s) after NormDispIncr 1e-8/1e-6, dt/2, dt/4x4" \
						$i $eqNstepsAll [getTime]]
				}
			}
			if {$eqPrintON} {
				set tNow [getTime]
				set elapsed [expr {([clock microseconds] - $t0)/1.0e6}]
				while {$tNow >= $tPrint} {
					puts [format "  analysis t=%.3f s  elapsed=%.2f s  pier top ux=%.4e m  uy=%.4e m  step %d / %d" \
						$tNow $elapsed \
						[nodeDisp $nodeTag_pierTop_deckBC 1] [nodeDisp $nodeTag_pierTop_deckBC 2] \
						$i $eqNstepsAll]
					set tPrint [expr {$tPrint + $dtPrint}]
				}
			}
		}
	}
	set elapsed [expr {([clock microseconds] - $t0)/1.0e6}]
	set soilEQDone 1

	if {$realTimeON} {
		puts [format "----- EQ done  realTimeON  nSteps=%s -----" $realTimeNsteps]
	} else {
		puts [format "----- EQ done  recoveries=%d -----" $nFail]
	}
	puts [format "  analysis t=%.4g s  elapsed=%.2f s  (start %.4g s + EQ %.4g s + freeVib %.4g s)" \
		[getTime] $elapsed $tWait $Trec $eqFreeVibT]
	puts [format "  pier top ux=%.4e m  uy=%.4e m" \
		[nodeDisp $nodeTag_pierTop_deckBC 1] [nodeDisp $nodeTag_pierTop_deckBC 2]]
	if {$recordersON == 0 && !$realTimeON} {
		puts "  recordersON=0"
	} else {
		puts [format "  outDIR=%s" $eqOutDir]
	}
	puts "Run: gravity + EQ done"
}
