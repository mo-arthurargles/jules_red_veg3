! (c) Crown copyright, Met Office. All rights reserved.
PROGRAM test_red_canopy

USE ancil_info, ONLY: land_pts, nmasst, dim_cslayer, dim_cs1, ainfo_type
USE jules_surface_types_mod, ONLY: nnpft, npft, soil
USE jules_vegetation_mod, ONLY: l_red, l_triffid, triffid_period,              &
  phenol_period, frac_min, frac_seed, l_limit_canhc
USE jules_soil_mod, ONLY: dzsoil
USE timestep_mod, ONLY: timestep
USE pftparm, ONLY: g_leaf_0, kext
USE trif, ONLY: g_root, g_wood
USE red_io, ONLY: alpha_recrt, crwn_area0, dom_order, height0, lai_bal0,       &
  mass0, massi, mclass, mort_base, phi_a, phi_g, phi_h, phi_l, k_cai
USE veg3_parm_mod, ONLY: veg3_parm_allocate, veg3_parm_init, veg3_ctrl,        &
  red_parms, litter_parms
USE veg3_field_mod, ONLY: veg3_field_allocate, veg_state, red_state,         &
  red_veg3_couple
USE veg3_canopy_mod, ONLY: veg3_canopy_height_order,                         &
  veg3_canopy_maintain_frac_min
USE veg3_litter_mod, ONLY: veg3_litter
USE veg3_red_assimilate_parition_mod, ONLY: veg3_red_assimilate_partition
USE veg3_red_dynamic_mod, ONLY: veg3_red_dynamic
USE cancap_mod, ONLY: cancap
USE jules_surface_mod, ONLY: hleaf, hwood
USE conversions_mod, ONLY: rsec_per_day

IMPLICIT NONE

TYPE(ainfo_type) :: ainfo
REAL :: qtop, qmid, qsmall, wtop, wsmall, available, expected, year
REAL :: old_frac(1,2), old_carbon, npp(1,2), turnover(1,2), leaf_rate(1,2)
REAL :: old_leaf_litter(1,2), expected_npp, explicit_litter
REAL :: heat_capacity(1), cover(1), crop_dvi(1,0)
REAL :: min_test(2), k_test(2), excl_test(2), tile_excl_test(1)
REAL :: pnd_test(1,2,2), frac_test(1,2)
INTEGER :: n

land_pts = 1
nnpft = 2
npft = 2
nmasst = 2
soil = 3
dim_cslayer = 1
dim_cs1 = 4
l_red = .TRUE.
l_triffid = .TRUE.
l_limit_canhc = .FALSE.
triffid_period = 1
phenol_period = 1
timestep = 86400.0
frac_min = 0.001
frac_seed = 0.001
ALLOCATE(g_leaf_0(2), g_root(2), g_wood(2), dzsoil(1), kext(2))
g_leaf_0 = 0.5
g_root = 0.2
g_wood = 0.01
dzsoil = 1.0
kext = 0.5
alpha_recrt = 0.2
crwn_area0 = 1.0
dom_order = 1.0
height0 = 10.0
lai_bal0 = 1.0
mass0 = 10.0
massi = 100.0
mclass = 2
mort_base = 0.0
phi_a = 0.5
phi_g = 0.75
phi_h = 0.25
phi_l = 0.25
k_cai = 1.0
CALL veg3_parm_allocate(1,3,2,2)
CALL veg3_parm_init(1,3,2,2,2)
CALL veg3_field_allocate(1,3,2,2,4)

red_state%mass_cohort(:,1) = 10.0
red_state%mass_cohort(:,2) = 100.0
red_state%ht_mass(1,:) = (/10.0,30.0/)
red_state%ht_mass(2,:) = (/5.0,20.0/)
red_state%crwn_area_mass(:,1) = 1.0
red_state%crwn_area_mass(:,2) = 2.0
red_state%lai_bal_mass(:,1) = 1.0
red_state%lai_bal_mass(:,2) = 4.0
red_state%leaf_mass(:,1) = 0.2
red_state%leaf_mass(:,2) = 1.6
red_state%froot_mass = red_state%leaf_mass
red_state%wood_mass = red_state%mass_cohort - 2.0 * red_state%leaf_mass
red_state%g_mass_scale(:,1) = 1.0
red_state%g_mass_scale(:,2) = 3.0
CALL veg3_canopy_height_order(2,2,4,red_parms%mclass,red_state%ht_mass,        &
  red_state%order_pft_desc,red_state%order_mclass_desc,                       &
  red_state%order_pft_asc,red_state%order_mclass_asc)
veg_state%phen = 0.5
red_state%plantNumDensity = 0.0
red_state%plantNumDensity(1,1,2) = 0.6
CALL red_veg3_couple(ainfo)
CALL near('single-class height',veg_state%canht(1,1),30.0)
CALL near('single-class LAI',veg_state%lai_bal(1,1),4.0)
CALL near('phenological LAI',veg_state%lai(1,1),2.0)
CALL near('single-class canopy wood',veg_state%woodC_canopy(1,1),48.4)
CALL near('empty PFT carbon',veg_state%vegCpft(1,2),0.0)
CALL near('empty PFT physical fallback',veg_state%canht(1,2),5.0)

red_state%plantNumDensity(1,1,1) = 0.8
CALL red_veg3_couple(ainfo)
qtop = 1.0 - EXP(-1.2)
qsmall = 1.0 - EXP(-0.8)
expected = (qtop * 30.0 + qsmall * 10.0) / (qtop + qsmall)
CALL near('intra overlap retained in height',veg_state%canht(1,1),expected)

red_state%plantNumDensity(1,2,2) = 0.4
CALL red_veg3_couple(ainfo)
available = 1.0 - red_parms%frac_excl(1)
qmid = 1.0 - EXP(-0.8)
wtop = available * qtop
wsmall = (available - (available - wtop) * qmid) * qsmall
expected = (wtop * 30.0 + wsmall * 10.0) / (wtop + wsmall)
CALL near('inter overlap excluded from height',veg_state%canht(1,1),expected)
CALL near('exposed plus intra class weight',red_state%frac_intra_mass(1,1,1),wsmall)
expected = (wtop * 4.0 + wsmall) / (wtop + wsmall)
CALL near('intra mean LAI',veg_state%lai_bal(1,1),expected)
CALL near('whole population gridbox carbon',veg_state%vegC(1),108.0)
DO n = 1,2
  expected = SUM(red_state%plantNumDensity(1,n,:) * red_state%wood_mass(n,:))
  CALL near('all shaded wood retained',veg_state%woodC(1,n) * veg_state%frac(1,n),expected)
END DO

CALL cancap(1,1,(/1/),3,1,veg_state,veg_state%canht(:,1),                   &
  veg_state%lai(:,1),heat_capacity,cover,crop_dvi)
expected = hleaf * veg_state%leafC_canopy(1,1) + hwood * veg_state%woodC_canopy(1,1)
CALL near('physical heat capacity',heat_capacity(1),expected)

year = rsec_per_day * 360.0
leaf_rate = 0.5 / year
CALL veg3_litter(1,(/1/),veg3_ctrl,1,2,litter_parms,leaf_rate,veg_state,turnover)
DO n = 1,2
  expected = SUM(red_state%plantNumDensity(1,n,:) *                          &
    (leaf_rate(1,n) * 0.5 * red_state%leaf_mass(n,:) +                       &
     litter_parms%g_root(n) * red_state%froot_mass(n,:) +                    &
     litter_parms%g_wood(n) * red_state%wood_mass(n,:)))
  CALL near('whole population turnover',turnover(1,n) * veg_state%frac(1,n),expected)
END DO

npp = 0.0
CALL veg3_red_assimilate_partition(1,2,red_parms,npp,turnover,               &
  veg_state%frac(:,1:2),red_state)
DO n = 1,2
  CALL near('negative assimilate has no recruits',red_state%P_s(1,n),0.0)
  CALL near('all deficit becomes shrinkage',red_state%G_struct(1,n),red_state%P_a(1,n))
END DO

npp = turnover + 1.0e-6
CALL veg3_red_assimilate_partition(1,2,red_parms,npp,turnover,               &
  veg_state%frac(:,1:2),red_state)
DO n = 1,2
  CALL near('only exposed area assimilates',red_state%P_a(1,n),1.0e-6 * veg_state%frac(1,n))
  expected = SUM(red_state%plantNumDensity(1,n,:) * red_state%g_mass(1,n,:))
  CALL near('MST growth conserves assimilate',expected,red_state%G_struct(1,n))
  CALL near('MST class scaling retained',red_state%g_mass(1,n,2) / red_state%g_mass(1,n,1),3.0)
END DO

old_frac = veg_state%frac(:,1:2)
old_carbon = veg_state%vegC(1)
old_leaf_litter = veg_state%leaf_litC * old_frac
expected_npp = SUM(npp * old_frac)
veg_state%npp_dr_out = npp * year
red_state%mort = 0.0
CALL veg3_red_dynamic(veg3_ctrl%dt_red,1,(/1/),veg3_ctrl,1,2,2,             &
  red_parms,veg_state,red_state)
explicit_litter = SUM(turnover * old_frac) +                               &
  SUM(veg_state%mort_litC * veg_state%frac(:,1:2))
CALL red_veg3_couple(ainfo,old_frac)
IF (MAXVAL(ABS(veg_state%frac(:,1:2) - old_frac)) < 1.0e-6) THEN
  PRINT *, 'Test must exercise a changing top-down PFT fraction'
  STOP 1
END IF
CALL near('NPP gridbox flux survives area change',veg_state%npp_gb(1),expected_npp)
CALL near('litter components survive area change',                         &
  SUM(veg_state%leaf_litC * veg_state%frac(:,1:2)),SUM(old_leaf_litter))
CALL near('whole population budget',veg_state%vegC(1) - old_carbon,         &
  (expected_npp - veg_state%litC(1) / year) * veg3_ctrl%dt_red)
CALL near('implicit and explicit litter agree',veg_state%litC(1) / year,explicit_litter,2.0e-9)
CALL near('soil receives reported litter',                                &
  SUM(veg_state%litCpft * veg_state%frac(:,1:2)),veg_state%litC(1))

k_test = 1.0
tile_excl_test = 0.0
excl_test = 0.0
min_test = (/0.2, 0.0/)
pnd_test = 0.0
pnd_test(1,1,2) = -LOG(0.9) / red_state%crwn_area_mass(1,2)
CALL veg3_canopy_maintain_frac_min(1,2,2,4,red_state%order_pft_desc,          &
  red_state%order_mclass_desc,min_test,k_test,excl_test,tile_excl_test,       &
  red_state%crwn_area_mass,red_state%mass_cohort(:,1),pnd_test,frac_test)
CALL near('minimum includes existing taller cohort',frac_test(1,1),0.2)
CALL near('minimum preserves existing taller cohort density',                &
  pnd_test(1,1,2),-LOG(0.9) / red_state%crwn_area_mass(1,2))

min_test = 0.2
excl_test = 0.2
pnd_test = 0.0
CALL veg3_canopy_maintain_frac_min(1,2,2,4,red_state%order_pft_desc,          &
  red_state%order_mclass_desc,min_test,k_test,excl_test,tile_excl_test,       &
  red_state%crwn_area_mass,red_state%mass_cohort(:,1),pnd_test,frac_test)
CALL near('minimum retained after taller PFT recruitment',frac_test(1,1),0.2)
CALL near('minimum retained after shorter PFT recruitment',frac_test(1,2),0.2)

k_test(1) = 0.0
min_test = 0.0
excl_test = 0.0
pnd_test = 0.0
CALL veg3_canopy_maintain_frac_min(1,2,2,4,red_state%order_pft_desc,          &
  red_state%order_mclass_desc,min_test,k_test,excl_test,tile_excl_test,       &
  red_state%crwn_area_mass,red_state%mass_cohort(:,1),pnd_test,frac_test)
CALL near('zero k_cai with zero minimum',frac_test(1,1),0.0)

PRINT *, 'RED canopy and carbon coupling tests passed'

CONTAINS

SUBROUTINE near(label,actual,wanted,absolute_tolerance)
CHARACTER(LEN=*), INTENT(IN) :: label
REAL, INTENT(IN) :: actual,wanted
REAL, INTENT(IN), OPTIONAL :: absolute_tolerance
REAL :: tolerance

tolerance = 2.0e-5 * MAX(ABS(wanted),1.0e-8)
IF (PRESENT(absolute_tolerance)) tolerance = absolute_tolerance
IF (.NOT. (ABS(actual - wanted) <= tolerance)) THEN
  PRINT *, label, ': actual=',actual, ', expected=',wanted
  STOP 1
END IF
END SUBROUTINE near

END PROGRAM test_red_canopy
