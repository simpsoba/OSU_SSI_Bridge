# Parameters.tcl
# Edit knobs here. IDs: TAGS CONVENTION. Driver: Run.tcl.
# Units: N, m, s
#
# Where to edit (Rayleigh groups, soil mesh, GM, mass, pier SP, elastic,
# fibers, free vib): README.md.
# User-changeable knobs: `# <-- EDIT` (physics, switches, block bases).
# Grep:  grep -n "EDIT\|USER INPUT" Parameters.tcl analysis/RayleighDamping.tcl
# Run.tcl also marks `runEQ` / `plotFigures`. Derived values are not EDIT.

# ------------------------------------------------------------
# Units
# ------------------------------------------------------------
set inch 0.0254;                          # m
set foot [expr {12.0*$inch}];             # m
set pi   3.141592653589793;               # (-)
set gravity_accel 9.81;                   # <-- EDIT  m/s^2

# ---- Lab cylinder -- Neumann (2021); Neumann et al. (2023) ----
set D_cyl [expr {20.0*$inch}];            # m, specimen diameter (20 in)  # <-- EDIT
set H_cyl 3.02;                           # m, specimen height (Neumann Hc_m)  # <-- EDIT

# ---- Prototype pier -- Shin et al. (2007); Mackie et al. (2008) ----
set D_pier [expr {4.0*$foot}];            # m, pier diameter (4 ft)  # <-- EDIT

# Length scale (prototype / lab) and pier height: (4 ft = 48 in) / 20 in
set cylinderSF [expr {(4.0*12.0)/20.0}];  # (-) D_pier / D_cyl
set H_pier [expr {$H_cyl*$cylinderSF}];   # m

# ------------------------------------------------------------
# Model switches
# ------------------------------------------------------------
#   pierEleType / pileEleType / soilProfile / soilBoundary -- below
#   soilConstitutive: inelastic = PIMY/PDMY02 (+ FSP sands)
#                     elastic   = ElasticIsotropic (skeleton Gr,Br; FSP on sands)
#   pileSpring:       inelastic | elastic | none  (py + tz + qz together)
#   soilEleType:      quad | SSPquad  (near-field continuum; ASDEA ring unchanged)
set soilConstitutive "inelastic";         # <-- EDIT  inelastic | elastic
# elastic -> ElasticIsotropic3D from skeleton Gr,Br (FSP on sands; needs getCopy)
set pileSpring "inelastic";               # <-- EDIT  inelastic | elastic | none
set soilEleType "SSPquad";                   # <-- EDIT  quad | SSPquad
# SSPquad: no ele rho -- mass from nDMaterial getRho() (PIMY/PDMY/ElasticIsotropic3D rho).
# Body forces still set on the element (b1,b2).

# ------------------------------------------------------------
# Pier
# ------------------------------------------------------------
#   pierEleType:
#     elasticBeamColumn  -- single elastic pier
#     forceBeamColumn    -- ConcentratedCurvature (Fiber hinges + Elastic mid)
#     lumpedPlasticity   -- two Fiber ZLS + eta*EI beam (hybrid)
set pierEleType "lumpedPlasticity";       # <-- EDIT  elasticBeamColumn | forceBeamColumn | lumpedPlasticity
set pierGeoTransf "PDelta";               # <-- EDIT  geomTransf type

# lumpedPlasticity: two Fiber ZLS + eta*I beam (hybrid: beam ~ 100x EI)
#   1 --ZLS-I-- 2 ==== eta EI ==== 4 --ZLS-J-- 5
#   K_theta ~ Ec I_uncr / Ls ; Ls = zlsLeqRatio * H_pier
#   Priestley Lp still scales post-peak (eps_u - eps_peak)
#   Profile 3 Shin quad eigen vs FBC: T1=2.06 s, extra pier ~0.34 s
#   at zlsLeqRatio_I=0.74, zlsLeqRatio_J=0.25 (eta=100). Lp still Priestley.
set eta_pier        100.0;                # <-- EDIT  (-) stiff-beam factor
set zlsLeqRatio_I   0.74;                 # <-- EDIT  (-) Ls_I / H_pier  (base)
set zlsLeqRatio_J   0.25;                 # <-- EDIT  (-) Ls_J / H_pier  (top)

# Pier flexural EI (elastic members): uncracked transformed, then x cracked factor.
#   I_uncr = I_g + (Es/Ec - 1) I_s   (I_s = 1/2 As R_bar^2 for ring bars)
#   I_pier = pierCrackedFactor * I_uncr   -> elasticBeamColumn / FBC mid / eta*I beam
# Fiber hinge Ec, Es are full; lumped ZLS: Ls = zlsLeqRatio*H (elastic), Lp post-peak.
set pierCrackedFactor 0.5;                # <-- EDIT  (-) 1 = uncracked transformed

# Fiber strips along section depth (fiber pier types); denser at extremes
set nFiberY_pier    21;                   # <-- EDIT  (-) total concrete strips
set nFiberEdge_pier 5;                    # <-- EDIT  (-) strips in each extreme band

# ---- Pier materials ----
# Concrete
set fc_pier   28.0e6;                     # <-- EDIT  Pa, compressive strength
set eps0_pier 0.002;                      # <-- EDIT  (-) strain at peak stress
set epsu_pier 0.005;                      # <-- EDIT  (-) cover spalling / crushing strain
set dens_c_pier 2400.0;                   # <-- EDIT  kg/m^3, concrete mass density

# Steel (SteelMPF)
set fy_pier      470.0e6;                 # <-- EDIT  Pa, yield stress
set Es_pier      200.0e9;                 # <-- EDIT  Pa, elastic modulus
set b_steel_pier 0.01;                    # <-- EDIT  (-) post-yield stiffness ratio
set R0_steel_pier  20.0;                  # <-- EDIT  (-) SteelMPF transition
set cR1_steel_pier 0.925;                 # <-- EDIT  (-) SteelMPF
set cR2_steel_pier 0.15;                  # <-- EDIT  (-) SteelMPF
set dens_s_pier 7850.0;                   # <-- EDIT  kg/m^3, steel mass density

# ---- Pier reinforcement ----
set cover_pier   [expr {2.0*$inch}];      # <-- EDIT  m, clear cover
set longbar_pier 10;                      # <-- EDIT  (-) US longitudinal bar size (#)
set tranbar_pier 7;                       # <-- EDIT  (-) US tie / hoop bar size (#)
set s_tran_pier  [expr {3.5*$inch}];      # <-- EDIT  m, tie spacing (#7 @ 3.5 in)

# 28 x #10 in pairs -> 14 bundles around the ring
set n_long_pier     28;                   # <-- EDIT  (-) total longitudinal bars
set n_bundle_pier   2;                    # <-- EDIT  (-) bars per bundle
set n_long_pos_pier [expr {$n_long_pier/$n_bundle_pier}];  # (-) ring positions

set db_long_pier   [expr {($longbar_pier/8.0)*$inch}];           # m, long. bar diameter
set db_tran_pier   [expr {($tranbar_pier/8.0)*$inch}];           # m, tie diameter
set As_long_pier   [expr {$pi/4.0*$db_long_pier*$db_long_pier}];  # m^2, one #10
set As_bundle_pier [expr {$n_bundle_pier*$As_long_pier}];         # m^2, fiber area at each ring point

# Derived ratios / core size (not inputs)
set A_g_pier   [expr {$pi/4.0*$D_pier*$D_pier}];               # m^2, gross area
set rho_l_pier [expr {$n_long_pier*$As_long_pier/$A_g_pier}];  # (-) long. steel ratio
# d_cs as in CircularColumn (for rho_t / Mander); R_core to inside of cover
set dcs_pier    [expr {$D_pier - 2.0*$cover_pier - $db_tran_pier}];  # m
set R_core_pier [expr {0.5*$D_pier - $cover_pier - 0.5*$db_long_pier - $db_tran_pier}];  # m
set rho_t_pier  [expr {$pi*$db_tran_pier*$db_tran_pier/$s_tran_pier/$dcs_pier}];  # (-)
# Linear mass density: concrete on (A_g - As) + steel on As
set As_tot_pier [expr {$n_long_pier*$As_long_pier}];           # m^2
set rhoL_pier   [expr {$dens_c_pier*($A_g_pier - $As_tot_pier) + $dens_s_pier*$As_tot_pier}];  # kg/m

# Priestley / Caltrans plastic hinge length (cantilever: L = H_pier)
#   SI:  Lp = 0.08 L + 0.022 fy db  >=  0.044 fy db
#   with fy in MPa, db and L in m; 0.022 and 0.044 have units MPa^-1
#   (same as Caltrans/Mackie 2003: 0.08 L' + 0.15 fy dbl >= 0.3 fy dbl, fy ksi, dbl in)
set Lp_pier [expr {0.08*$H_pier + 0.022*($fy_pier*1.0e-6)*$db_long_pier}];  # m
set Lp_min_pier [expr {0.044*($fy_pier*1.0e-6)*$db_long_pier}];               # m
if {$Lp_pier < $Lp_min_pier} {
	set Lp_pier $Lp_min_pier
}

# ZLS equivalent lengths (lumpedPlasticity). Lp is Priestley, unchanged.
if {$eta_pier <= 1.0} {
	error "Parameters: eta_pier must be > 1"
}
if {$pierCrackedFactor <= 0.0} {
	error "Parameters: pierCrackedFactor must be > 0"
}
if {$zlsLeqRatio_I <= 0.0 || $zlsLeqRatio_J <= 0.0} {
	error "Parameters: zlsLeqRatio_I and zlsLeqRatio_J must be > 0"
}
set Ls_I_pier [expr {$zlsLeqRatio_I*$H_pier}];  # m, base hinge elastic scale
set Ls_J_pier [expr {$zlsLeqRatio_J*$H_pier}];  # m, top hinge elastic scale

# ------------------------------------------------------------
# Pile cap
# ------------------------------------------------------------
# Mackie et al. (2008) Type 1A / Ketchum 3x2 under 4 ft column
set H_cap [expr {3.25*$foot}];            # <-- EDIT  m, thickness (3.25 ft)
set W_cap [expr {15.0*$foot}];            # <-- EDIT  m, plan width transverse (15 ft, long)
set L_cap [expr {10.0*$foot}];            # <-- EDIT  m, plan length out-of-plane (10 ft, short)
set E_cap 2.0e11;                         # <-- EDIT  Pa, steel E (near-rigid frame)
set s_pile_cap [expr {6.0*$foot}];        # <-- EDIT  m, center->outer (6 ft c/c = 3 D_pile)
set dens_cap $dens_c_pier;                # kg/m^3

# Cap mass / A / I use W_cap (prototype). Frame and face-spring x stop at
# the outer pile axes (s); W_cap_soil = 2s. Overhang to +/- W/2 is tributary
# mass only (no extra soil columns).
set A_cap [expr {$H_cap*$W_cap}];         # m^2
set I_cap [expr {$H_cap*pow($W_cap,3)/12.0}];  # m^4
set m_cap [expr {$dens_cap*$H_cap*$W_cap*$L_cap}];  # kg
# Cap frame: 3x3 pile grid + 2 soffit mid-bay = 11
# m_i and Irot_i: tributary rectangle in BuildPileCapNodes.tcl
set n_cap_nodes 11;                       # (-) node count

# ------------------------------------------------------------
# Piles
# ------------------------------------------------------------
#   pileEleType:
#     elasticBeamColumn -- A, E, I (x n_pile_row)
#     dispBeamColumn    -- Fiber tube strips (graded), areas x n_pile_row
set pileEleType "dispBeamColumn";         # <-- EDIT  elasticBeamColumn | dispBeamColumn
set pileGeoTransf "PDelta";               # <-- EDIT  geomTransf type
set D_pile [expr {2.0*$foot}];            # <-- EDIT  m, outer diameter (2 ft)
set t_pile [expr {0.5*$inch}];            # <-- EDIT  m, wall thickness (0.5 in)
set n_pile 3;                             # <-- EDIT  (-) shafts in plane
set n_pile_row 2;                         # <-- EDIT  (-) piles into the page (area x2)
set L_pile [expr {60.0*$foot}];           # <-- EDIT  m, length below cap (60 ft, Mackie)
set nSeg_pile 20;                         # <-- EDIT  (-) segments -> dy = 3 ft (exact on 60 ft)
set nIP_pile 3;                           # <-- EDIT  (-) dispBeamColumn integration points

# Fiber strips (dispBeamColumn); denser at extreme fibers
set nFiberY_pile    21;                   # <-- EDIT  (-) total tube strips
set nFiberEdge_pile 5;                    # <-- EDIT  (-) strips in each extreme band

# Steel (Steel01 for pile Fiber)
set fy_pile 340.0e6;                      # <-- EDIT  Pa, yield (ASTM A572-class)
set Es_pile 200.0e9;                      # <-- EDIT  Pa
set b_pile  0.01;                         # <-- EDIT  (-) post-yield ratio
set dens_s_pile $dens_s_pier;             # kg/m^3

# Derived pipe section (one pipe, then x n_pile_row)
set Ro_pile [expr {0.5*$D_pile}];         # m
set Ri_pile [expr {$Ro_pile - $t_pile}];  # m
# Mesh / face-spring x: outer pile axes (no +/-R face ring)
set xCapSoilHalf $s_pile_cap;             # m
set W_cap_soil [expr {2.0*$xCapSoilHalf}];         # m
set A_pipe_one [expr {$pi*($Ro_pile*$Ro_pile - $Ri_pile*$Ri_pile)}];  # m^2
set I_pipe_one [expr {$pi/4.0*(pow($Ro_pile,4) - pow($Ri_pile,4))}];  # m^4
set A_pile [expr {$n_pile_row*$A_pipe_one}];  # m^2
set I_pile [expr {$n_pile_row*$I_pipe_one}];  # m^4
set rhoL_pile [expr {$dens_s_pile*$A_pile}];  # kg/m (steel wall only)

# ------------------------------------------------------------
# Deck
# ------------------------------------------------------------
#   No fiber. Factors on A,I only (E = Ec). Mass: rho A L_trib, length-weighted
#   nodal lump; Irot from member consistent-mass diagonal (see BuildDeckNodes).
set dw_deck [expr {39.0*$foot}];          # <-- EDIT  m, top width out-to-out
set dd_deck [expr {6.0*$foot}];           # <-- EDIT  m, depth out-to-out
set sw_deck [expr {23.0*$foot}];          # <-- EDIT  m, soffit width
set cw_deck [expr {5.5*$foot}];           # <-- EDIT  m, cantilever
set td_deck [expr {9.5*$inch}];           # <-- EDIT  m, top slab
set ts_deck [expr {8.0*$inch}];           # <-- EDIT  m, soffit slab
set tw_deck [expr {12.0*$inch}];          # <-- EDIT  m, web (nominal)
set bh_deck [expr {32.0*$inch}];          # <-- EDIT  m, jersey barrier height
set yb_deck [expr {43.6*$inch}];          # <-- EDIT  m, CG above soffit bottom
# Ketchum A / Iy (Deck_PT_39); cracked flexure target = 0.5 Iy (not assigned per member)
set A_deck [expr {8869.0*6.4516e-4}];     # <-- EDIT  m^2
set Iy_deck [expr {325.5*0.008631}];      # <-- EDIT  m^4
set fc_deck 34473.8e3;                    # <-- EDIT  Pa, deck f'c
set Ec_deck [expr {4700.0*sqrt($fc_deck/1.0e6)*1.0e6}];  # Pa
set dens_deck [expr {22.78e3/9.81}];      # <-- EDIT  kg/m^3 (gamma=22.78 kN/m^3)
set L_trib_deck [expr {150.0*$foot}];     # <-- EDIT  m, one-pier OOP tributary (= t_soil)
set m_deck [expr {$dens_deck*$A_deck*$L_trib_deck}];  # kg
# Mass assignment (BuildDeckNodes): length-weighted mx=my; Irot = Sum m_mem L^2/105
# (pier-style). Cap-style Iz-Steiner fill is not used for the hollow box.

# Member base props (strip-scale); eta multiplies A and I only
set A_deck_mem 0.75;                      # <-- EDIT  m^2, base frame member area
set I_deck_mem 0.02;                      # <-- EDIT  m^4, base frame member I
set eta_deck_A 100.0;                     # <-- EDIT  (-) -> A_mem = eta A
set eta_deck_I 1000.0;                    # <-- EDIT  (-) -> I_mem = eta I
# Center web / pier link factors on pier A,I (stiff link; E stays Ec_deck)
set eta_deckLink_A 10.0;                  # <-- EDIT  (-)
set eta_deckLink_I 100.0;                 # <-- EDIT  (-)

set A_deck_frame [expr {$eta_deck_A*$A_deck_mem}];  # m^2
set I_deck_frame [expr {$eta_deck_I*$I_deck_mem}];  # m^4

# ------------------------------------------------------------
# Soil domain -- Profiles.md
# ------------------------------------------------------------
#   soilProfile: 1 = L2 crust clay + L3 liq sand (Dr ramp) + L5 dense sand
#                2 = L2 crust clay + L3 med clay + L5 dense sand
#                3 = L2 crust clay + L3 med clay + L5 stiff clay
#                4 = thicker soft L2 (to -27.25 ft) + med L3 + stiff L5
#   soilBoundary:
#     Shin  -- thick FF columns + 3 Lysmer dashpots (rock, F=2c v)
#     ASDEA -- ASDAbsorbingBoundary2D ring (L/BL/B/BR/R); rock bottom,
#             layer G0 on sides; setParameter stage 1 after gravity
#   Refs / Stage0 -fx / Abell setTime: soil/Boundary.md
set soilProfile 4;                        # <-- EDIT  1 | 2 | 3 | 4
set soilBoundary "Shin";                  # <-- EDIT  Shin | ASDEA
# soilConstitutive / pileSpring: see Model switches above
set rho_w 1000.0;                         # <-- EDIT  kg/m^3
set h_water 2.4;                          # <-- EDIT  m, free water above y=0 (0 = none)
# Ponding: consistent Fy on soil top edges only (WaterSurfaceLoad.tcl).
# No structure hydrostatic. Body forces unchanged (gamma'/gamma).
set t_soil [expr {150.0*$foot}];          # <-- EDIT  m, near-field OOP thickness

# Vertical mesh (keep dy_soil = pile dy = L_pile / nSeg_pile)
set dy_soil [expr {3.0*$foot}];           # <-- EDIT  m, vertical quad size (= pile dy)
set nSeg_below_tip 5;                     # <-- EDIT  (-) L5 below tip (15 ft)

# Horizontal mesh, |x| from the pier (near field). Each row:
#   {mesh size    x end for that size}
# mesh size -- horizontal quad width in this ring
# x end     -- outer |x| where this width stops (next row starts there)
# x end = previous x end + n * mesh size (integer n), or the last cell is skinny.
# Last x end is L_half (near-field outer face).
#
# soilMesh: pick one band list
#    0  production (~35 x-stations with Shin)
#    1  fine      (3 ft bands to 201 ft NF; ~137 x-stations)
#    2  finer     (3 ft bands to 270 ft NF)
#   -1  coarse    (~25 x-stations)
#   -2  coarser   (~19 x-stations)
set soilMesh 1;                           # <-- EDIT  -2 | -1 | 0 | 1 | 2

if {$soilMesh == 0} {
	# production: 3 ft SSI to 12 ft
	set soilDxBands [list \
		[list [expr { 3.0*$foot}] [expr { 12.0*$foot}]] \
		[list [expr { 7.0*$foot}] [expr { 40.0*$foot}]] \
		[list [expr {15.0*$foot}] [expr {100.0*$foot}]] \
		[list [expr {20.0*$foot}] [expr {140.0*$foot}]] \
		[list [expr {30.0*$foot}] [expr {200.0*$foot}]] \
		]
} elseif {$soilMesh == 1} {
	# fine: 3 ft bands to 201 ft NF
	set soilDxBands [list \
		[list [expr { 3.0*$foot}] [expr { 12.0*$foot}]] \
		[list [expr { 3.0*$foot}] [expr { 39.0*$foot}]] \
		[list [expr { 3.0*$foot}] [expr { 99.0*$foot}]] \
		[list [expr { 3.0*$foot}] [expr {141.0*$foot}]] \
		[list [expr { 3.0*$foot}] [expr {201.0*$foot}]] \
		]
} elseif {$soilMesh == 2} {
	# finer: 3 ft bands to 270 ft NF
	set soilDxBands [list \
		[list [expr { 3.0*$foot}] [expr { 12.0*$foot}]] \
		[list [expr { 3.0*$foot}] [expr { 39.0*$foot}]] \
		[list [expr { 3.0*$foot}] [expr { 99.0*$foot}]] \
		[list [expr { 3.0*$foot}] [expr {141.0*$foot}]] \
		[list [expr { 3.0*$foot}] [expr {270.0*$foot}]] \
		]
} elseif {$soilMesh == -1} {
	# coarse: same inner SSI; ~25 x-stations
	set soilDxBands [list \
		[list [expr { 3.0*$foot}] [expr { 12.0*$foot}]] \
		[list [expr {14.0*$foot}] [expr { 40.0*$foot}]] \
		[list [expr {20.0*$foot}] [expr {100.0*$foot}]] \
		[list [expr {50.0*$foot}] [expr {200.0*$foot}]] \
		]
} elseif {$soilMesh == -2} {
	# coarser: 3 ft only to outer pile (±s); ~19 x-stations
	set soilDxBands [list \
		[list [expr { 3.0*$foot}] [expr {  6.0*$foot}]] \
		[list [expr {12.0*$foot}] [expr { 30.0*$foot}]] \
		[list [expr {30.0*$foot}] [expr { 90.0*$foot}]] \
		[list [expr {55.0*$foot}] [expr {200.0*$foot}]] \
		]
} else {
	error "Parameters.tcl: soilMesh must be -2, -1, 0, 1, or 2 (got '$soilMesh')"
}
set L_half [lindex [lindex $soilDxBands end] 1];  # m, NF outer face

# Far-field column: Shin uses one thick column of width w_FF; ASDEA uses
# an absorbing ring of the same width. OOP t_FF is Shin-only.
set w_FF [expr {40.0*$foot}];             # <-- EDIT  m, FF column / ring width
set t_FF_factor 10000.0;                  # <-- EDIT  (-) t_FF / t_soil (Shin only)
set t_FF [expr {$t_FF_factor*$t_soil}];   # m

set B_fsp 2.2e9;                          # <-- EDIT  Pa, FSP combined bulk (OpenSees wiki B_f)

# Rock half-space under L5 (Lysmer / ASDEA bottom; Joyner & Chen outcrop scale)
set rockVs 760.0;                         # <-- EDIT  m/s
set rockRho 2400.0;                       # <-- EDIT  kg/m^3
set rockNu 0.25;                          # <-- EDIT  (-) ASDEA elastic
set rockG [expr {$rockRho*$rockVs*$rockVs}];  # Pa
set asdeaNu 0.25;                         # <-- EDIT  (-) sides (not continuum bulk nu)

# ------------------------------------------------------------
# Ground motion / EQ
# ------------------------------------------------------------
# Base velocity Path (PEER VT2 or plain m/s file). Empty -> dummy zeros.
# VT2 is cm/s -- analysis/BuildVelSeries converts to m/s. PEER header sets DT.
# Uncomment one pair (Tohoku | El Centro).
set gmRoot [file join [file dirname [file normalize [info script]]] ground-motion]
# Tohoku 2011, KiK-net FKSH19 borehole NS
set gmDir [file join $gmRoot Tohoku2011-FKSH]
set gmVelFile [file join $gmDir FKSH19.NS1.VT2];  # <-- EDIT  Path file (or "")
# El Centro 1940 Array #9, 180 (NS)
# set gmDir [file join $gmRoot ImperialValley1940-ElCentro]
# set gmVelFile [file join $gmDir RSN6_IMPVALL.I_I-ELC180.VT2]
set gmVelDT 0.01;                         # <-- EDIT  s (PEER header overrides if present)
set gmScaleFactor 1.0;                    # <-- EDIT  (-) Path scale (units converted first)

# EQ transient step (TRBDF2). Time scale is sqrt(cylinderSF) (length scale).
# dtAnalysis = DT_FACTOR / 2048 * sqrt(cylinderSF)
set DT_FACTOR 10;                         # <-- EDIT  (-) skip vs 2048 Hz
set dtAnalysis [expr {$DT_FACTOR/2048.0*sqrt($cylinderSF)}];  # s
# Optional truncate (s); empty -> full velocity record length
set eqTmax "";                            # <-- EDIT  s (empty = full record)
# Free vibration after the last GM sample (Path is 0 past the record).
# Lysmer/ASDEA then only radiate. 0 -> stop when the earthquake record ends.
set eqFreeVibT 60.0;                      # <-- EDIT  s
set nModesEigen 10;                       # <-- EDIT  (-) modes after gravity (runEQ 0)
# recordersON (EQ; both Run.tcl and RunParallel.tcl). Every recorder samples at
# -dT gmVelDT (the PEER step), not at every dtAnalysis step:
#   0  off
#   1  full window: |x|<=eqWindowX nodes + quads; all pile beams, all SSI springs
#   2  center column: pier nodes 1/2/4/5 (UX UY RZ), both rotational springs,
#      soil-base primary, the whole center pile, every center-pile spring,
#      and every x=0 soil quad (grade to base). No cap springs, no pier accel.
#   3  nine SSI horizons (old 2): first / mid / last station of L2, L3, L5
set recordersON 2;                        # <-- EDIT  0 | 1 | 2 | 3
set eqWindowX 10.0;                       # <-- EDIT  m, |x| <= this for deformed-shape dump (recordersON 1)

# ------------------------------------------------------------
# Springs
# ------------------------------------------------------------
# Group effect on p-y only (Mokwa chart @ 3D; EQ average of 3 in-plane rows)
#   f_m(3D) = 0.82, 0.67, 0.58 -> Ge = 0.69; OOP side-by-side @ 3D -> 1.0
set Ge_pile 0.69;                         # <-- EDIT  (-)
set num_spring_pCap 6;                    # <-- EDIT  (-) 3 elev x 2 faces

# API sand k (pci) from soil/api-k-vs-Dr.csv (below-WT), PDMY knots
set k_pci_L3a 60.0;                       # <-- EDIT  Dr 50% knot
set k_pci_L3b 81.0;                       # <-- EDIT  Dr 60% knot
set k_pci_L3c 121.0;                      # <-- EDIT  Dr 75% knot
set k_pci_L5  170.0;                      # <-- EDIT  Dr ~95%

# Spring / dashpot
set Cd_py 1.0;                            # <-- EDIT  (-) PySimple1 drag
set c_dash_py 0.01;                       # <-- EDIT  (-) dashpot
set pRes_frac 0.15;                       # <-- EDIT  (-) PyLiq1 / TzLiq residual fraction
set z50_tz 0.001;                         # <-- EDIT  m, shaft t-z
set y50_cap 0.01;                         # <-- EDIT  m, cap p-y / t-z
set z50_cap 0.01;                         # <-- EDIT  m
set alpha_cap_py 0.75;                    # <-- EDIT  (-) Mokwa phi=0 wall adhesion factor

# =====================================================================
# TAGS CONVENTION
# =====================================================================
# Mesh encodings stay in the builders; only the bases are here.
# Soil: node = tagShift_soil + ix*nY + iy  (BuildSoilMesh sets stride = nY).
# Springs / ASDEA / Lysmer: if a default base is already taken, that builder
# moves it to the next thousand above max(getNodeTags) / max(getEleTags).
#
# Shared nodes (same tag, stacked mass, no equalDOF):
#   cap TC = nodeTag_pierBase_capTC; deck soffit BC = nodeTag_pierTop_deckBC;
#   pile heads = cap BL / BC / BR.

# ---- Block bases ----
#   pier     1-999
#   cap      1000    eles 1101+
#   piles    2000    nodes 2001+/2101+/2201+; eles 2100+
#   deck     3000    nodes 3001+; eles 3100+
#   soil     10000   nodes 10000+; quads 15000+
#   springs  20000   dups + zeroLength; Py/Tz mats  (may move)
#   ASDEA/Lysmer mates  30000 / eles 35000          (may move)
set tagShift_cap  1000;                   # <-- EDIT
set tagShift_pile 2000;                   # <-- EDIT
set tagShift_deck 3000;                   # <-- EDIT
set tagShift_soil 10000;                  # <-- EDIT
set tagShift_spr  20000;                  # <-- EDIT
set nodeTag_bnd_base 30000;               # <-- EDIT  ASDEA outer / Lysmer mates
set eleTag_bnd_base  35000;               # <-- EDIT  ASDEA + Lysmer zeroLength

# ---- Pier (edit these; everything else follows) ----
# lumpedPlasticity:
#   pierBase_capTC --ZLS-I-- pierBaseZeroLengthInner
#     ==== eta*EI ====
#   pierTopZeroLengthInner --ZLS-J-- pierTop_deckBC
# elastic / forceBeamColumn: pierBase_capTC and pierTop_deckBC only.
set nodeTag_pierBase_capTC             1; # y=0; pile-cap top center
set nodeTag_pierBaseZeroLengthInner    2; # same as base; lumpedPlasticity
set nodeTag_pierTopZeroLengthInner     4; # same as top; lumpedPlasticity
set nodeTag_pierTop_deckBC             5; # y=H_pier; deck soffit center
set eleTag_pier_botSpr   1;               # zeroLengthSection base 1->2
set eleTag_pier          2;               # eta beam 2->4 (or 1->5)
set eleTag_pier_topSpr   3;               # zeroLengthSection top 4->5
set transfTag_pier       1;               # geomTransf
set intTag_pier          1;               # ConcentratedCurvature
set secTag_pier        101;               # Fiber (forceBeamColumn hinges)
set secTag_pier_I      101;               # Fiber ZLS base
set secTag_pier_J      103;               # Fiber ZLS top
set secTag_elastic_pier 102;              # Elastic mid-span
set matTag_cover_pier  201;
set matTag_core_pier   202;
set matTag_steel_pier  203;
set matTag_cover_pier_I 211;
set matTag_core_pier_I  212;
set matTag_steel_pier_I 213;
set matTag_cover_pier_J 221;
set matTag_core_pier_J  222;
set matTag_steel_pier_J 223;

# ---- Pile cap ----
# T/M/B = top / mid / bot (y = 0, -H/2, -H)
# L/C/R = left / center / right (pile axes at x = +/-s, 0)
# BML/BMR = bot mid-bay, x = +/-s/2
#
#   y=0     TL ---------- TC ---------- TR
#            |           / | \           |
#   y=-H/2  ML ---------- MC ---------- MR
#            |         /    |    \        |
#   y=-H    BL -- BML -- BC -- BMR -- BR
#           -s   -s/2     0    s/2     s
#
# TC = nodeTag_pierBase_capTC.  BL, BC, BR = pile heads.
# Face Py/Tz on TL/ML/BL and TR/MR/BR.
set nodeTag_cap_TL  [expr {$tagShift_cap + 21}];  # 1021  top left @ -s
set nodeTag_cap_TC  $nodeTag_pierBase_capTC;      # pier base (shared)
set nodeTag_cap_TR  [expr {$tagShift_cap + 23}];  # 1023  top right @ +s
set nodeTag_cap_ML  [expr {$tagShift_cap + 24}];  # 1024  mid left
set nodeTag_cap_MC  [expr {$tagShift_cap + 25}];  # 1025  mid center
set nodeTag_cap_MR  [expr {$tagShift_cap + 26}];  # 1026  mid right
set nodeTag_cap_BL  [expr {$tagShift_cap + 27}];  # 1027  pile head left
set nodeTag_cap_BC  [expr {$tagShift_cap + 28}];  # 1028  pile head center
set nodeTag_cap_BR  [expr {$tagShift_cap + 29}];  # 1029  pile head right
set nodeTag_cap_BML [expr {$tagShift_cap + 36}];  # 1036  bot mid-bay -s/2
set nodeTag_cap_BMR [expr {$tagShift_cap + 37}];  # 1037  bot mid-bay +s/2
set transfTag_cap   [expr {$tagShift_cap + 1}];   # 1001
set eleTag_cap_base [expr {$tagShift_cap + 101}]; # 1101  first cap element

# ---- Piles ----
# BuildPilesNodes: node = tagShift_pile + ipile*pileNodeStride + iy  (iy=1..nSeg)
#   left head 1027, nodes 2001...; center 1028 / 2101...; right 1029 / 2201...
#   pileNodeStride starts at 100 and grows if nSeg_pile >= 100
set nodeTag_pile_base $tagShift_pile
set pileNodeStride 100
set eleTag_pile_base  [expr {$tagShift_pile + 100}];  # 2100
set transfTag_pile    [expr {$tagShift_pile + 1}];    # 2001
set secTag_pile       301;                # Fiber or Elastic
set matTag_steel_pile 401;                # Steel01

# ---- Deck (fixed PT39 frame) ----
#   BarL 3009                                    BarR 3010
#      |                                            |
#   TL 3004 -- TLi 3005 -- TC 3006 -- TRi 3007 -- TR 3008     y = H_pier+dd
#      |          |          |          |          |
#              BL 3001 ---- BC 5 ---- BR 3003               y = H_pier
#
# BC = nodeTag_pierTop_deckBC.  TLi/TRi = inner web.  BarL/BarR = jersey top.
set nodeTag_deck_BL   3001
set nodeTag_deck_BC   $nodeTag_pierTop_deckBC
set nodeTag_deck_BR   3003
set nodeTag_deck_TL   3004
set nodeTag_deck_TLi  3005
set nodeTag_deck_TC   3006
set nodeTag_deck_TRi  3007
set nodeTag_deck_TR   3008
set nodeTag_deck_BarL 3009
set nodeTag_deck_BarR 3010
set transfTag_deck    3001
set eleTag_deck_base  3100

# ---- Soil continuum ----
# BuildSoilMesh: node = tagShift_soil + ix*soilNodeStride + iy  (ix,iy 0-based)
# Layer mats: one solid per soil row at matTag_soil_base+1+iy, FSP +100+iy
set nodeTag_soil_base $tagShift_soil
set eleTag_soil_base  [expr {$tagShift_soil + 5000}];  # 15000
set matTag_soil_base  501;                # solid + FSP wrappers
set matTag_lysmer_base [expr {$tagShift_soil + 50}];   # 10050  Viscous
set soilNodeStride 100;                   # overwritten to nY in BuildSoilMesh

# ---- Interface springs ----
# BuildSoilSprings: pile dup = tagShift_spr + ip*pileNodeStride + iyP;
#   cap face +sprCapFaceOff+i; soffit +sprSoffitOff+i
set nodeTag_sprSoil_base $tagShift_spr
set eleTag_spr_base      [expr {$tagShift_spr + 2000}];  # 22000
set matTag_py_base       [expr {$tagShift_spr + 3000}];  # 23000
set matTag_tz_base       [expr {$tagShift_spr + 4000}];  # 24000
set sprCapFaceOff 900
set sprSoffitOff  920

source [file join [file dirname [file normalize [info script]]] soil TagHelpers.tcl]

# ---- Time series / load patterns ----
set tsTag_velBase         9001;           # -fx / Shin Path
set tsTag_gravStruct      9101;           # structure weight Linear
set tsTag_hWater          9102;           # ponding Linear
set tsTag_holdPier        9103;           # Constant; post-gravity UX/UY hold
set patternTag_gravStruct 110;            # structure weight
set patternTag_hWater     111;            # ponding
set patternTag_holdPier   112;            # pier base UX/UY (after gravity)
set patternTag_lysmer     116;            # Shin 2c v

# Resolve python3 (Linux) or python (Windows) for plotters.
# Args: none
# Returns: executable path, or "" if neither is on PATH
proc FindPython3 {} {
	foreach name {python3 python} {
		set bin [auto_execok $name]
		if {$bin ne ""} {
			return $bin
		}
	}
	return ""
}
