
! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in Veg3 Ecosystem Demography
! *****************************COPYRIGHT****************************************
!
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise.

MODULE veg3_red_dynamic_mod

IMPLICIT NONE

PRIVATE
PUBLIC :: veg3_red_dynamic

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='VEG3_RED_DYNAMIC_MOD'

CONTAINS

!-----------------------------------------------------------------------------
SUBROUTINE veg3_red_dynamic(                                                   &
                !IN Control vars
                dt_red,veg_index_pts,veg_index,veg3_ctrl,land_pts,             &
                nnpft,nmasst,                                                  &
                !IN red_parms
                red_parms,                                                     &
                !IN state
                veg_state,red_state                                            &
                !OUT Diagnostics
                )

!Only get the data structures - the data comes through the calling tree
USE veg3_parm_mod,ONLY:  red_parm_type, veg3_ctrl_type
USE veg3_field_mod,ONLY:  veg_state_type, red_state_type
USE veg3_canopy_mod,ONLY: veg3_canopy_frac, veg3_canopy_maintain_frac_min

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Objects with INTENT IN
!-----------------------------------------------------------------------------
TYPE(red_state_type)  :: red_state
TYPE(red_parm_type)   :: red_parms
TYPE(veg_state_type)  :: veg_state
TYPE(veg3_ctrl_type)  :: veg3_ctrl

!----------------------------------------------------------------------------
! Integers with INTENT IN
!----------------------------------------------------------------------------
INTEGER, INTENT(IN) :: land_pts,nnpft,veg_index(land_pts),veg_index_pts,nmasst

!----------------------------------------------------------------------------
! Reals with INTENT IN
!----------------------------------------------------------------------------
REAL, INTENT(IN)   ::                                                          &
dt_red
              !  Dynamic vegetation time-step (s)

!-----------------------------------------------------------------------------
!Local Vars
!-----------------------------------------------------------------------------
INTEGER                ::l,n,k,j

REAL                   ::                                                      &
frac_shade(land_pts,nnpft),                                                    &
              !  Competitive shading of seedlings in each PFT. (-)
frac_above_mclass1(land_pts,nnpft),                                            &
              !  frac_above at the moment each PFT's lowest mass class was
              !  processed (see veg3_canopy_frac).
growth_up(nmasst),                                                             &
              !  Positive structural growth of each mass class, summed over
              !  the intra- and inter-PFT groups. [kg C m-2 s-1]
growth_down(nmasst)
              !  Magnitude of negative structural growth of each mass class,
              !  summed over the intra- and inter-PFT groups. [kg C m-2 s-1]


!End of headers

! Initialise vars
veg_state%mort_litC(:,:)    = 0.0
frac_shade(:,:)             = 0.0

! Pass 1: update each PFT's size structure (growth/mortality) independently.
! mort_litC remains gridbox-normalised at this point.
DO l = 1,land_pts
  DO n = 1,nnpft
    !Estimate the inter-PFT competition.
    DO j = 1, nnpft
      frac_shade(l,n) = MIN(1.0,frac_shade(l,n) + red_parms%comp_coef(n,j)     &
          * veg_state%frac(l,j))
    END DO

    ! Growth and shrinkage of the intra- and inter-PFT groups are kept apart
    ! so that opposite signs are not averaged before forming class fluxes.
    growth_up(:)   = 0.0
    growth_down(:) = 0.0
    DO k = 1, red_parms%mclass(n)
      growth_up(k) = red_state%plantNumDensity_intra(l,n,k)                    &
                     * MAX(0.0, red_state%g_mass_intra(l,n,k))                 &
                     + red_state%plantNumDensity_inter(l,n,k)                  &
                     * MAX(0.0, red_state%g_mass_inter(l,n,k))
      growth_down(k) = red_state%plantNumDensity_intra(l,n,k)                  &
                       * MAX(0.0, -red_state%g_mass_intra(l,n,k))              &
                       + red_state%plantNumDensity_inter(l,n,k)                &
                       * MAX(0.0, -red_state%g_mass_inter(l,n,k))
    END DO

    ! Update the PFT number density, using growth/P_s/mort precomputed by
    ! veg3_red_assimilate_partition/veg3_mortality_partition.
    CALL update_pft_size_structure(                                            &
      !IN sizing and control
      red_parms%mclass(n), dt_red,                                             &
      !IN fields
      red_state%P_s(l,n), frac_shade(l,n),                                     &
      red_state%mort(l,n,1:red_parms%mclass(n)),                               &
      growth_up(1:red_parms%mclass(n)), growth_down(1:red_parms%mclass(n)),    &
      red_state%mass_cohort(n,1:red_parms%mclass(n)),                          &
      !INOUT state
      red_state%plantNumDensity(l,n,1:red_parms%mclass(n)),                    &
      !OUT diagnostics
      veg_state%mort_litC(l,n)                                                 &
      )

  END DO
END DO

! Pass 2: enforce each PFT's minimum vegetation fraction using the
! overlap-corrected frac across all PFTs, debiting the implied recruitment
! from mort_litC (still gridbox-normalised - see item 8).
CALL veg3_canopy_frac(land_pts, nnpft, nmasst, veg3_ctrl%npft_totmclass,       &
  red_state%order_pft_desc, red_state%order_mclass_desc, red_parms%k_cai,      &
  red_parms%frac_excl, veg_state%frac_tile_excl, red_state%crwn_area_mass,     &
  red_state%plantNumDensity, red_state%frac_mass, red_state%CAI_mass,          &
  veg_state%frac(:,1:nnpft), veg_state%cai, frac_above_mclass1,                &
  red_state%CAI_intra_mass)

CALL veg3_canopy_maintain_frac_min(land_pts, nnpft, nmasst,                    &
  veg3_ctrl%npft_totmclass, red_state%order_pft_desc,                          &
  red_state%order_mclass_desc, red_parms%frac_min, red_parms%k_cai,            &
  red_parms%frac_excl, veg_state%frac_tile_excl, red_state%crwn_area_mass,     &
  red_state%mass_cohort(:,1), red_state%plantNumDensity,                       &
  veg_state%frac(:,1:nnpft), dt_red, veg_state%mort_litC)

! Pass 3: normalise mort_litC per unit PFT canopy area.
DO n = 1,nnpft
  DO l = 1,land_pts
    IF (veg_state%frac(l,n) > 0.0) THEN
      veg_state%mort_litC(l,n) = veg_state%mort_litC(l,n) /                    &
        veg_state%frac(l,n)
    ELSE
      veg_state%mort_litC(l,n) = 0.0
    END IF
  END DO
END DO

END SUBROUTINE veg3_red_dynamic
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
SUBROUTINE update_pft_size_structure(                                          &
                !IN sizing and control
                mclass,dt_red,                                                 &
                !IN fields
                P_s,frac_shade,mort,growth_up,growth_down,mass_cohort,         &
                !INOUT state
                plantNumDensity,                                               &
                !OUT diagnostics
                mort_litC                                                      &
                )

! Description:
!   Explicit update of one PFT's number density across its mass classes from
!   growth, shrinkage, recruitment and mortality.
!
! Method:
!   Upward and downward fluxes are formed from the summed positive and
!   negative growth of the overlap groups, so a class can grow and shrink at
!   once. All outflows from a class are scaled together so that the class
!   cannot be exhausted, conserving plants. The demographic litter is the
!   part of the assimilate not retained in the population.

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Integers with INTENT IN
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) :: mclass
              !  Number of mass classes for this PFT.

!-----------------------------------------------------------------------------
! Reals with INTENT IN
!-----------------------------------------------------------------------------
REAL, INTENT(IN)     ::                                                        &
dt_red,                                                                        &
              !  Dynamic vegetation time-step. [s]
P_s,                                                                           &
              !  Gridbox carbon assimilate devoted to recruitment.
              !  [kg C m-2 s-1]
frac_shade,                                                                    &
              !  Competitive shading of seedlings in this PFT. [-]
mort(mclass),                                                                  &
              !  Total mortality rate across the mass cohorts, precomputed by
              !  veg3_mortality_partition. [s-1]
growth_up(mclass),                                                             &
              !  Positive structural growth of each mass class. [kg C m-2 s-1]
growth_down(mclass),                                                           &
              !  Magnitude of negative structural growth of each mass class.
              !  [kg C m-2 s-1]
mass_cohort(mclass)
              !  Mass of an individual member of each mass cohort. [kg C]

!-----------------------------------------------------------------------------
! Reals with INTENT INOUT
!-----------------------------------------------------------------------------
REAL, INTENT(IN OUT)  ::                                                       &
plantNumDensity(mclass)
              !  Population density within each mass cohort. [m-2]

!-----------------------------------------------------------------------------
! Reals with INTENT OUT
!-----------------------------------------------------------------------------
REAL, INTENT(OUT)    :: mort_litC
              !  Mortality/demographic litter for this PFT, normalised per
              !  unit gridbox area, veg3_red_dynamics renormalises to
              !  PFT canopy area. [kg C m-2 s-1]

!-----------------------------------------------------------------------------
!Local Vars
!-----------------------------------------------------------------------------
INTEGER              :: k

REAL                 ::                                                        &
flux_up(mclass),                                                               &
              !  Plants growing from each class into the class above.
              !  [m-2 s-1]
flux_down(mclass),                                                             &
              !  Plants shrinking from each class into the class below.
              !  [m-2 s-1]
flux_mort(mclass),                                                             &
              !  Plants lost from each class by mortality. [m-2 s-1]
flux_recruit,                                                                  &
              !  Successful recruits into the lowest class. [m-2 s-1]
outflow,                                                                       &
              !  Total loss rate of plants from a class. [m-2 s-1]
dplantNumDensity_dt,                                                           &
              !  Net rate of change of a class's number density. [m-2 s-1]
carbon_stored
              !  Rate of change of the population's carbon. [kg C m-2 s-1]

!End of headers

flux_recruit = P_s * (1.0 - frac_shade) / mass_cohort(1)

DO k = 1, mclass
  ! Growth beyond the top class and shrinkage below the lowest class cannot
  ! move plants; the former goes to litter and the latter is met by
  ! starvation mortality (veg3_red_assimilate_partition).
  flux_up(k)   = 0.0
  flux_down(k) = 0.0
  IF (k < mclass) flux_up(k) = growth_up(k)                                    &
                               / (mass_cohort(k+1) - mass_cohort(k))
  IF (k > 1) flux_down(k) = growth_down(k)                                     &
                            / (mass_cohort(k) - mass_cohort(k-1))
  flux_mort(k) = mort(k) * plantNumDensity(k)

  ! Prevent the class from being exhausted over a timestep.
  outflow = flux_up(k) + flux_down(k) + flux_mort(k)
  IF (outflow * dt_red > plantNumDensity(k)) THEN
    flux_up(k)   = flux_up(k)   * plantNumDensity(k) / (outflow * dt_red)
    flux_down(k) = flux_down(k) * plantNumDensity(k) / (outflow * dt_red)
    flux_mort(k) = flux_mort(k) * plantNumDensity(k) / (outflow * dt_red)
  END IF
END DO

carbon_stored = 0.0
DO k = 1, mclass
  dplantNumDensity_dt = - flux_up(k) - flux_down(k) - flux_mort(k)
  IF (k == 1) dplantNumDensity_dt = dplantNumDensity_dt + flux_recruit
  IF (k > 1) dplantNumDensity_dt = dplantNumDensity_dt + flux_up(k-1)
  IF (k < mclass) dplantNumDensity_dt = dplantNumDensity_dt + flux_down(k+1)

  carbon_stored = carbon_stored + dplantNumDensity_dt * mass_cohort(k)
  plantNumDensity(k) = MAX(0.0,                                                &
                           plantNumDensity(k) + dplantNumDensity_dt * dt_red)
END DO

! Assimilate (P_s plus net structural growth) not retained in the population.
mort_litC = P_s + SUM(growth_up) - SUM(growth_down) - carbon_stored

END SUBROUTINE update_pft_size_structure
!-----------------------------------------------------------------------------

END MODULE veg3_red_dynamic_mod
