
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
MODULE veg3_litter_mod

IMPLICIT NONE

PRIVATE
PUBLIC :: veg3_litter

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='VEG3_LITTER_MOD'

CONTAINS
SUBROUTINE veg3_Litter(                                                        &
                !IN Control vars
                veg_index_pts,veg_index,land_pts,nnpft,                        &
                !IN parms
                litter_parms,                                                  &
                !IN fields
                g_leaf_dr,                                                     &
                !INOUT state
                veg_state,red_state                                            &
                )

! Description:
!   Local litterfall (leaf, root and wood turnover, excluding mortality) for
!   an individual in each mass class, aggregated to the PFT litter pools.
!   Turnover includes all plants, both intra- and inter-PFT overlapped.

!Only get the data structures - the data comes through the calling tree
USE conversions_mod, ONLY: rsec_per_day

USE veg3_parm_mod, ONLY:  litter_parm_type
USE veg3_field_mod, ONLY:  veg_state_type, red_state_type

IMPLICIT NONE

!----------------------------------------------------------------------------
! Integers with INTENT IN
!----------------------------------------------------------------------------
INTEGER, INTENT(IN) :: land_pts,nnpft,veg_index(land_pts),veg_index_pts

!-----------------------------------------------------------------------------
! Objects with INTENT IN
!-----------------------------------------------------------------------------
TYPE(litter_parm_type),INTENT(IN) :: litter_parms

!-----------------------------------------------------------------------------
! Reals with INTENT IN
!-----------------------------------------------------------------------------
REAL, INTENT(IN)       ::                                                      &
g_leaf_dr(land_pts,nnpft)
        ! Mean phenology-driven leaf turnover rate for driving vegetation
        ! dynamics, accumulated and averaged over the vegetation dynamics
        ! timestep in veg3_run_ctrl. [s-1]

!-----------------------------------------------------------------------------
! Objects with INTENT INOUT
!-----------------------------------------------------------------------------
TYPE(veg_state_type),INTENT(IN OUT)   :: veg_state
TYPE(red_state_type),INTENT(IN OUT)   :: red_state

!-----------------------------------------------------------------------------
!Local Vars
!-----------------------------------------------------------------------------
INTEGER                :: l,n,k,i,nmass

REAL                   ::                                                      &
leaf_ind,                                                                      &
        ! Leaf turnover of an individual. [kg C s-1]
root_ind,                                                                      &
        ! Fine root turnover of an individual. [kg C s-1]
wood_ind,                                                                      &
        ! Wood turnover of an individual. [kg C s-1]
leaf_litter,                                                                   &
        ! Leaf turnover per gridbox area. [kg C m-2 s-1]
root_litter,                                                                   &
        ! Fine root turnover per gridbox area. [kg C m-2 s-1]
wood_litter
        ! Wood turnover per gridbox area. [kg C m-2 s-1]

!End of headers

nmass = SIZE(red_state%leaf_mass,2)
red_state%litter_mass(:,:,:) = 0.0

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) COLLAPSE(2)                   &
!$OMP PRIVATE(l,n,k,i,leaf_ind,root_ind,wood_ind,leaf_litter,root_litter,      &
!$OMP         wood_litter)                                                     &
!$OMP SHARED(litter_parms,veg_state,red_state,veg_index,veg_index_pts,nnpft,   &
!$OMP        nmass,g_leaf_dr)
DO n = 1, nnpft
  DO i = 1, veg_index_pts
    l = veg_index(i)

    leaf_litter = 0.0
    root_litter = 0.0
    wood_litter = 0.0

    ! Unused mass classes have zero tissue mass and number density.
    DO k = 1, nmass
      leaf_ind = g_leaf_dr(l,n) * veg_state%phen(l,n) *                        &
                 red_state%leaf_mass(n,k)
      root_ind = litter_parms%g_root(n) * red_state%froot_mass(n,k)
      wood_ind = litter_parms%g_wood(n) * red_state%wood_mass(n,k)
      red_state%litter_mass(l,n,k) = leaf_ind + root_ind + wood_ind

      leaf_litter = leaf_litter + red_state%plantNumDensity(l,n,k) * leaf_ind
      root_litter = root_litter + red_state%plantNumDensity(l,n,k) * root_ind
      wood_litter = wood_litter + red_state%plantNumDensity(l,n,k) * wood_ind
    END DO

    ! Normalise by top-down PFT fraction and convert to kg C m-2 (360d)-1.
    IF (veg_state%frac(l,n) > 0.0) THEN
      veg_state%leaf_litC(l,n) = leaf_litter / veg_state%frac(l,n)             &
                                 * rsec_per_day * 360.0
      veg_state%root_litC(l,n) = root_litter / veg_state%frac(l,n)             &
                                 * rsec_per_day * 360.0
      veg_state%wood_litC(l,n) = wood_litter / veg_state%frac(l,n)             &
                                 * rsec_per_day * 360.0
    ELSE
      veg_state%leaf_litC(l,n) = 0.0
      veg_state%root_litC(l,n) = 0.0
      veg_state%wood_litC(l,n) = 0.0
    END IF

  END DO
END DO
!$OMP END PARALLEL DO

RETURN
END SUBROUTINE veg3_Litter

!-----------------------------------------------------------------------------

END MODULE veg3_litter_mod
