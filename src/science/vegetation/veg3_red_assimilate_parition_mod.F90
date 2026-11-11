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

MODULE veg3_red_assimilate_parition_mod

IMPLICIT NONE

PRIVATE
PUBLIC :: veg3_red_assimilate_partition

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName=                            &
  'VEG3_RED_ASSIMILATE_PARITION_MOD'

CONTAINS

!-----------------------------------------------------------------------------
SUBROUTINE veg3_red_assimilate_partition(                                      &
                !IN Control vars
                land_pts,nnpft,                                                &
                !IN red_parms
                red_parms,                                                     &
                !IN fields
                npp_dr,local_litter,CAI,                                       &
                !INOUT state
                red_state                                                      &
                )

! Estimates each PFT's gridbox growth (P_a) from npp_dr/local_litter/CAI,
! splits it into recruitment (P_s) and structural growth (G_struct), where
! P_a = P_s + G_struct, then scales G_struct across the mass class
! structure (g_mass), ahead of veg3_red_dynamic.

!Only get the data structures - the data comes through the calling tree
USE veg3_parm_mod,  ONLY: red_parm_type
USE veg3_field_mod, ONLY: red_state_type

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Objects with INTENT IN OUT
!-----------------------------------------------------------------------------
TYPE(red_state_type)  :: red_state
TYPE(red_parm_type)   :: red_parms

!----------------------------------------------------------------------------
! Integers with INTENT IN
!----------------------------------------------------------------------------
INTEGER, INTENT(IN) :: land_pts,nnpft

!----------------------------------------------------------------------------
! Reals with INTENT IN
!----------------------------------------------------------------------------
REAL, INTENT(IN)   ::                                                          &
npp_dr(land_pts,nnpft),                                                        &
              !  Mean NPP driving vegetation, per unit CAI. (kgC m-2 s-1)
local_litter(land_pts,nnpft),                                                  &
              !  Litter production, per unit CAI. (kgC m-2 s-1)
CAI(land_pts,nnpft)
              !  PFT crown area index. (m2 m-2)

!-----------------------------------------------------------------------------
!Local Vars
!-----------------------------------------------------------------------------
INTEGER                ::l,n,k,mclass

REAL                   ::                                                      &
g0,                                                                            &
              !  Boundary growth for an individual member of the smallest mass
              !  cohort. (kgC s-1)
plantNumDensity_g_sum
              !  Summation of the relative cohort contribution towards the
              !  total PFT assimilate (m-2)

!End of headers

red_state%P_a(:,:)      = 0.0
red_state%G_struct(:,:) = 0.0
red_state%P_s(:,:)      = 0.0
red_state%g_mass(:,:,:) = 0.0

DO l = 1,land_pts
  DO n = 1,nnpft
    mclass = red_parms%mclass(n)

    ! Sum product of the number density and the allometric scaling
    plantNumDensity_g_sum = 0.0
    DO k = 1, mclass
      plantNumDensity_g_sum = plantNumDensity_g_sum                            &
        + red_state%plantNumDensity(l,n,k) * red_state%g_mass_scale(n,k)
    END DO

    ! Total gridbox growth (P_a), split into recruitment (P_s) and
    ! structural growth (G_struct): P_a = P_s + G_struct.
    red_state%P_a(l,n) = (npp_dr(l,n) - local_litter(l,n)) * CAI(l,n)
    red_state%P_s(l,n) = red_parms%alpha_recrt(n) * red_state%P_a(l,n)
    red_state%G_struct(l,n) = (1.0 - red_parms%alpha_recrt(n))                 &
      * red_state%P_a(l,n)

    g0 = 0.0
    IF (plantNumDensity_g_sum > 0) g0 = red_state%G_struct(l,n)                &
      / plantNumDensity_g_sum

    ! Recruitment does not apply when growth is negative - the associated
    ! shrinkage is instead applied via g0/g_mass_scale metabolic (phi_g)
    ! scaling in update_pft_size_structure.
    IF (red_state%P_a(l,n) < 0.0) red_state%P_s(l,n) = 0.0

    ! Metabolic (phi_g) scaling of growth across the mass classes.
    red_state%g_mass(l,n,1:mclass) = g0 * red_state%g_mass_scale(n,1:mclass)
  END DO
END DO

END SUBROUTINE veg3_red_assimilate_partition
!-----------------------------------------------------------------------------

END MODULE veg3_red_assimilate_parition_mod
