# soil/UpdateSoilStage.tcl
# Goals: updateMaterialStage 0 (elastic) or 1 (plastic) on continuum + FSP.
# Call from SoilGravity / ActivateEQBoundary. Liq springs stay 0 until EQ
# (updateLiqSpringStage).
#
# =====================================================================
# 3. MATERIALS AND SECTIONS
# =====================================================================
# (stage update only; materials created in BuildSoilMaterials / BuildSoilSprings)

if {![info exists soilMatStageWanted]} {
	error "UpdateSoilStage.tcl: set soilMatStageWanted to 0 or 1 first"
}
if {$soilMatStageWanted != 0 && $soilMatStageWanted != 1} {
	error "UpdateSoilStage.tcl: soilMatStageWanted must be 0 or 1"
}
if {![info exists soilProfile] || ![info exists soilConstitutive]} {
	error "UpdateSoilStage.tcl: source Parameters + BuildSoilMaterials first"
}

set nSolid 0
set nFsp 0
set nLiq 0

# Continuum stage only for multi-yield (+ FSP). ElasticIsotropic3D has no stage.
if {$soilConstitutive eq "inelastic"} {
	if {![info exists soilSolidTags]} {
		error "UpdateSoilStage.tcl: soilSolidTags from BuildSoilMaterials.tcl required"
	}
	foreach m $soilSolidTags {
		# updateMaterialStage -material $matTag -stage $stage
		updateMaterialStage -material $m -stage $soilMatStageWanted
		incr nSolid
	}
	if {[info exists soilFspTags]} {
		foreach m $soilFspTags {
			updateMaterialStage -material $m -stage $soilMatStageWanted
			incr nFsp
		}
	}
}

# PyLiq1 / TzLiq1: only when updateLiqSpringStage is set (ActivateEQBoundary).
# Default off so accidental source during gravity does not flip them.
if {![info exists updateLiqSpringStage]} {
	set updateLiqSpringStage 0
}
if {$updateLiqSpringStage && [info exists liqSpringMatTags] && [llength $liqSpringMatTags] > 0} {
	foreach m $liqSpringMatTags {
		updateMaterialStage -material $m -stage $soilMatStageWanted
		incr nLiq
	}
}

set soilMatStage $soilMatStageWanted
if {$updateLiqSpringStage} {
	puts [format "----- soil stage %d  solids=%d  FSP=%d  liqSpr=%d -----" \
		$soilMatStageWanted $nSolid $nFsp $nLiq]
}
