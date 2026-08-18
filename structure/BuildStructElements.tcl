# structure/BuildStructElements.tcl
# Goals: create pier, deck, cap, and pile beam-columns (nodes already exist).
# Call after FoldStructNodes (Run.tcl) or immediately from PlotModel.tcl.

set structDirHere [file dirname [file normalize [info script]]]
source [file join $structDirHere BuildPierElements.tcl]
source [file join $structDirHere BuildDeckElements.tcl]
source [file join $structDirHere BuildPileCapElements.tcl]
source [file join $structDirHere BuildPilesElements.tcl]
