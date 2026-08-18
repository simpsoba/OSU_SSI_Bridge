# BuildModel.tcl
# Goals: structure nodes + soil quads + SSI springs as needed.
# Pier/deck/cap/pile beam-columns wait until after soil gravity (FoldStructNodes).
# Knobs: Parameters.tcl. Tags: TAGS CONVENTION.
if {![info exists root] || ![info exists structDir] || ![info exists soilDir]} {
	error "BuildModel.tcl: set root, structDir, soilDir first"
}
if {![info exists H_pier] || ![info exists pierEleType]} {
	error "BuildModel.tcl: source Parameters.tcl first"
}

# =====================================================================
# 2. MODEL BUILDER / NODES
# 3. MATERIALS AND SECTIONS
# 4. ELEMENTS
# =====================================================================

# model BasicBuilder -ndm $ndm -ndf $ndf
model BasicBuilder -ndm 2 -ndf 3
set structNodeTags {}
source [file join $structDir PierSection.tcl]
source [file join $structDir BuildPierNodes.tcl]
source [file join $structDir BuildDeckNodes.tcl]
source [file join $structDir BuildPileCapNodes.tcl]
source [file join $structDir BuildPileSection.tcl]
source [file join $structDir BuildPilesNodes.tcl]

model BasicBuilder -ndm 2 -ndf 2
source [file join $soilDir BuildSoilMaterials.tcl]
source [file join $soilDir BuildSoilMesh.tcl]

model BasicBuilder -ndm 2 -ndf 3
source [file join $soilDir BuildSoilSprings.tcl]

set buildModelDone 1
set nSprR 0
if {[info exists nSprings]} { set nSprR $nSprings }
puts [format "----- Build  structNodes=%d  nX=%d nY=%d nQuad=%d  spr=%d  %s %s %s -----" \
	[llength $structNodeTags] $soil_nX $soil_nY $soil_nQuad $nSprR \
	$pierEleType $soilBoundary $soilEleType]
