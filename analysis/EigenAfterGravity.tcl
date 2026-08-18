# analysis/EigenAfterGravity.tcl
# Goals: mode-shape JSON + PNG after `eigen`. No-op if plotFigures=0.

if {!$plotFigures} {
	return
}
if {![info exists plotDir]} {
	set plotDir [file join $root plot]
}
source [file join $plotDir DumpEigenModes.tcl]
set python3bin [FindPython3]
if {$python3bin eq ""} {
	puts "EigenAfterGravity: WARNING Python not found; JSON only"
} elseif {[info exists eigenOutPath]} {
	if {[catch {exec {*}$python3bin [file join $root plot PlotEigenModes.py] \
		$eigenOutPath} err]} {
		puts "EigenAfterGravity: WARNING PlotEigenModes.py failed:\n$err"
	}
}
