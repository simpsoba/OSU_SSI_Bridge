# PlotModel.tcl
# Units: N, m, s
#
# Build full SSI mesh, dump JSON, run Python plotters. No gravity / EQ.
#   OpenSees PlotModel.tcl
#
# Optional: set ::plotSkipPython 1  (JSON only). Knobs: Parameters.tcl.

set scriptDir [file dirname [file normalize [info script]]]
set root $scriptDir
set structDir [file join $root structure]
set soilDir [file join $root soil]
set analysisDir [file join $root analysis]
set plotDir [file join $root plot]

wipe

source [file join $root Parameters.tcl]

source [file join $root BuildModel.tcl]
source [file join $structDir BuildStructElements.tcl]

# Shin Lysmer for sketch; ASDEA ring already Stage 0 on mesh
model BasicBuilder -ndm 2 -ndf 2
if {$soilBoundary eq "Shin"} {
	source [file join $analysisDir BuildVelSeries.tcl]
	source [file join $soilDir ActivateEQBoundary.tcl]
}

file mkdir $plotDir
file mkdir [file join $plotDir out]
set caseDir [file join $plotDir out profile$soilProfile]
set elevDir [file join $caseDir elevation $soilBoundary]
set sprDir  [file join $caseDir pile_springs]
set soilProfDir [file join $caseDir soil_profile]
set fibDir  [file join $caseDir fibers]
foreach dirPath [list $caseDir $elevDir $sprDir $soilProfDir $fibDir] {
	file mkdir $dirPath
}
set sketchOutPath      [file join $elevDir model_sketch.json]
set fiberOutPath       [file join $fibDir fiber_sections.json]
set pileSpringOutPath  [file join $sprDir pile_springs.json]
set soilProfileOutPath [file join $soilProfDir soil_profile.json]
source [file join $plotDir DumpModelSketch.tcl]
source [file join $plotDir DumpFiberSections.tcl]
source [file join $plotDir DumpPileSprings.tcl]
source [file join $plotDir DumpSoilProfile.tcl]
foreach {src dst} [list \
	$sketchOutPath     [file join $plotDir model_sketch.json] \
	$fiberOutPath      [file join $plotDir fiber_sections.json] \
	$pileSpringOutPath [file join $plotDir pile_springs.json] \
	$soilProfileOutPath [file join $plotDir soil_profile.json] \
] {
	file copy -force $src $dst
}
puts [format "PlotModel: dumps -> %s  springs=%s" $caseDir $nSprings]

if {![info exists ::plotSkipPython] || !$::plotSkipPython} {
	set python3bin [FindPython3]
	if {$python3bin eq ""} {
		puts "PlotModel: WARNING Python not found; JSON dumps only"
	} else {
		cd $root
		foreach plotScript {
			plot/PlotModelSketch.py
			plot/PlotFiberSections.py
			plot/PlotPileSprings.py
			plot/PlotSoilProfile.py
		} {
			puts "PlotModel: python3 $plotScript"
			if {[catch {exec {*}$python3bin $plotScript} err]} {
				puts "PlotModel: WARNING $plotScript failed:\n$err"
			}
		}
	}
}
puts "PlotModel: done"
