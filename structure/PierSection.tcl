# PierSection.tcl
# Units: N, m, s
#
# Goals: pier section(s) from Parameters.tcl.
#   elasticBeamColumn -- Elastic A, I (uncracked transformed x pierCrackedFactor)
#   forceBeamColumn   -- Fiber hinge + Elastic mid-span
#   lumpedPlasticity  -- Fiber ZLS at each end
#       elastic: moduli / Ls; strains through peak x Ls
#                 Ls_I = zlsLeqRatio_I * H ; Ls_J = zlsLeqRatio_J * H
#       plastic: (eps_u - eps_peak) x Lp (Priestley) at both ends
#
# Fiber hinges keep full Ec, Es. Concrete02IS + Mander (CircularColumn), Pa.
# SteelMPF. Ec = 4700 sqrt(f'c[MPa]) MPa.
# Expects Parameters.tcl. Needs a model domain.

if {![info exists D_pier] || ![info exists pierEleType] || ![info exists secTag_pier]} {
	error "PierSection.tcl: source Parameters.tcl first"
}
if {![info exists pierCrackedFactor] || ![info exists As_tot_pier] || ![info exists R_core_pier]} {
	error "PierSection.tcl: need pierCrackedFactor, As_tot_pier, R_core_pier from Parameters.tcl"
}
if {![info exists pi]} {
	set pi 3.141592653589793
}
if {$pierCrackedFactor <= 0.0} {
	error "PierSection: pierCrackedFactor must be > 0"
}

set Ec_pier [expr {4700.0*sqrt($fc_pier/1.0e6)*1.0e6}];  # Pa
set A_pier  [expr {$pi/4.0*$D_pier*$D_pier}];              # m^2
set I_g_pier [expr {$pi/64.0*pow($D_pier,4)}];             # m^4, gross concrete
set R_pier  [expr {0.5*$D_pier}];                          # m

# Uncracked transformed I about Ec: I_g + (Es/Ec - 1) I_s
# Ring bars at R_bar = R_core: I_s = 1/2 As R_bar^2
set I_s_pier [expr {0.5*$As_tot_pier*$R_core_pier*$R_core_pier}];  # m^4
set I_uncr_pier [expr {$I_g_pier + ($Es_pier/$Ec_pier - 1.0)*$I_s_pier}];  # m^4
set I_pier [expr {$pierCrackedFactor*$I_uncr_pier}];       # m^4, elastic EI only

# One Fiber hinge: materials + section.
# Ls -- elastic scale (1 for forceBeamColumn)
# Lp -- plastic scale (1 for forceBeamColumn; Priestley for ZLS)
# Split at peak: eps_peak' = eps_peak*Ls ; eps_u' = eps_peak*Ls + (eps_u-eps_peak)*Lp
# Args: secTag matCover matCore matSteel Ls Lp label Ec Es fy b R0 cR1 cR2
#       fc fcc fcu eps0 epsu epsc epsuC fct fctC Ets coreFibers coverFibers rebarFibers
# Returns: none (creates materials + Fiber section)
proc pierBuildFiberHinge {secTag matCover matCore matSteel Ls Lp \
	label Ec Es fy b R0 cR1 cR2 \
	fc fcc fcu eps0 epsu epsc epsuC fct fctC Ets \
	coreFibers coverFibers rebarFibers} {

	if {$Ls <= 0.0 || $Lp <= 0.0} {
		error "pierBuildFiberHinge: Ls and Lp must be > 0"
	}

	# Ec, Es = material moduli (full); Ls scales for ZLS only
	set Ec_mat  [expr {$Ec/$Ls}]
	set Ets_mat [expr {$Ets/$Ls}]
	set Es_mat  [expr {$Es/$Ls}]
	set eps0_mat  [expr {$eps0*$Ls}]
	set epsc_mat  [expr {$epsc*$Ls}]
	set epsu_mat  [expr {$eps0*$Ls + ($epsu - $eps0)*$Lp}]
	set epsuC_mat [expr {$epsc*$Ls + ($epsuC - $epsc)*$Lp}]

	uniaxialMaterial Concrete02IS $matCover \
		$Ec_mat \
		[expr {-$fc}] [expr {-$eps0_mat}] \
		0.0 [expr {-$epsu_mat}] \
		0.1 $fct $Ets_mat
	uniaxialMaterial Concrete02IS $matCore \
		$Ec_mat \
		[expr {-$fcc}] [expr {-$epsc_mat}] \
		[expr {-$fcu}] [expr {-$epsuC_mat}] \
		0.1 $fctC $Ets_mat
	uniaxialMaterial SteelMPF $matSteel \
		$fy $fy $Es_mat \
		$b $b $R0 $cR1 $cR2

	# Build fiber body with tags already substituted (Fiber eval is not proc-local)
	set fibBody ""
	foreach fiberData $coreFibers {
		lassign $fiberData yLoc zLoc Af
		append fibBody [format "fiber %.12e %.12e %.12e %d\n" $yLoc $zLoc $Af $matCore]
	}
	foreach fiberData $coverFibers {
		lassign $fiberData yLoc zLoc Af
		append fibBody [format "fiber %.12e %.12e %.12e %d\n" $yLoc $zLoc $Af $matCover]
	}
	foreach fiberData $rebarFibers {
		lassign $fiberData yLoc zLoc Af
		append fibBody [format "fiber %.12e %.12e %.12e %d\n" $yLoc $zLoc $Af $matSteel]
	}
	section Fiber $secTag $fibBody
}

if {$pierEleType eq "elasticBeamColumn"} {

	section Elastic $secTag_pier $Ec_pier $A_pier $I_pier

} elseif {$pierEleType eq "forceBeamColumn" || $pierEleType eq "lumpedPlasticity"} {

	set here [file dirname [file normalize [info script]]]
	source [file join $here CircleStripFibers.tcl]

	if {$R_core_pier <= 0.0 || $R_core_pier >= $R_pier} {
		error "PierSection: bad R_core_pier=$R_core_pier (check cover / bars)"
	}
	set R_bar $R_core_pier;                               # m, long. bar centerline

	# --- Mander confinement (CircularColumn formulas, Pa) ---
	set As_tot  [expr {$n_long_pier*$As_long_pier}];       # m^2
	set sprime  [expr {$s_tran_pier - $db_tran_pier}];     # m
	set rho_cc  [expr {$As_tot/($pi/4.0*$dcs_pier*$dcs_pier)}];  # (-)
	set ke      [expr {pow(1.0 - $sprime/2.0/$dcs_pier, 2)/(1.0 - $rho_cc)}]
	set fl      [expr {0.5*$ke*$rho_t_pier*$fy_pier}];     # Pa
	set fcc_pier [expr {$fc_pier*(-1.254 + 2.254*sqrt(1.0 + 7.94*$fl/$fc_pier) - 2.0*$fl/$fc_pier)}]
	if {$fcc_pier > 2.23*$fc_pier} {
		set fcc_pier [expr {2.23*$fc_pier}]
	}
	set epsc_core [expr {$eps0_pier*(1.0 + 5.0*($fcc_pier/$fc_pier - 1.0))}]
	set ecr       [expr {$Ec_pier/($Ec_pier - $fcc_pier/$epsc_core)}]
	set epsu_core [expr {0.004 + 1.4*0.12*$fy_pier/$fc_pier*$rho_t_pier}]
	if {$epsu_core > 0.025} {
		set epsu_core 0.025
	}
	set fcu_core [expr {$fcc_pier*$epsu_core/$epsc_core*$ecr/($ecr - 1.0 + pow($epsu_core/$epsc_core, $ecr))}]
	set fct_pier  [expr {0.14*$fc_pier}];                  # Pa
	set fct_core  [expr {0.14*$fcc_pier}];                 # Pa
	set Ets_pier  [expr {$fct_pier/$eps0_pier}];           # Pa

	lassign [circularCoreCoverFiberStripsGraded \
		$R_pier $R_core_pier $nFiberY_pier $nFiberEdge_pier] \
		coreFibers coverFibers
	set rebarFibers [circularRebarYFibers $R_bar $n_long_pos_pier $As_bundle_pier]

	if {$pierEleType eq "forceBeamColumn"} {

		pierBuildFiberHinge $secTag_pier \
			$matTag_cover_pier $matTag_core_pier $matTag_steel_pier \
			1.0 1.0 "FBC" \
			$Ec_pier $Es_pier $fy_pier $b_steel_pier \
			$R0_steel_pier $cR1_steel_pier $cR2_steel_pier \
			$fc_pier $fcc_pier $fcu_core \
			$eps0_pier $epsu_pier $epsc_core $epsu_core \
			$fct_pier $fct_core $Ets_pier \
			$coreFibers $coverFibers $rebarFibers

		section Elastic $secTag_elastic_pier $Ec_pier $A_pier $I_pier

	} else {

		if {![info exists Ls_I_pier] || ![info exists Ls_J_pier] || ![info exists Lp_pier]} {
			error "PierSection: lumpedPlasticity needs Ls_I_pier, Ls_J_pier, Lp_pier"
		}

		pierBuildFiberHinge $secTag_pier_I \
			$matTag_cover_pier_I $matTag_core_pier_I $matTag_steel_pier_I \
			$Ls_I_pier $Lp_pier "ZLS-I" \
			$Ec_pier $Es_pier $fy_pier $b_steel_pier \
			$R0_steel_pier $cR1_steel_pier $cR2_steel_pier \
			$fc_pier $fcc_pier $fcu_core \
			$eps0_pier $epsu_pier $epsc_core $epsu_core \
			$fct_pier $fct_core $Ets_pier \
			$coreFibers $coverFibers $rebarFibers

		pierBuildFiberHinge $secTag_pier_J \
			$matTag_cover_pier_J $matTag_core_pier_J $matTag_steel_pier_J \
			$Ls_J_pier $Lp_pier "ZLS-J" \
			$Ec_pier $Es_pier $fy_pier $b_steel_pier \
			$R0_steel_pier $cR1_steel_pier $cR2_steel_pier \
			$fc_pier $fcc_pier $fcu_core \
			$eps0_pier $epsu_pier $epsc_core $epsu_core \
			$fct_pier $fct_core $Ets_pier \
			$coreFibers $coverFibers $rebarFibers
	}

} else {
	error "PierSection: pierEleType must be elasticBeamColumn, forceBeamColumn, or lumpedPlasticity (got '$pierEleType')"
}
