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
                npp_dr,frac,                                                   &
                !INOUT state
                red_state                                                      &
                )

! Description:
!   Partitions each PFT's gridbox assimilate across its mass classes and its
!   intra- and inter-PFT overlap groups, ahead of veg3_red_dynamic.
!
! Method:
!   NPP is shared across intra-PFT plants by metabolic (phi_g) scaling;
!   inter-PFT plants receive none. Every plant loses its own local
!   litterfall (litter_mass, from veg3_litter). Recruitment
!   P_s = alpha*MAX(0,P_a) is debited from the intra group only. A deficit
!   in the lowest class cannot shrink, so becomes starvation mortality.

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
              !  Mean NPP driving vegetation, per top-down tile area.
              !  [kg C m-2 s-1]
frac(land_pts,nnpft)
              !  Exposed top-down PFT fraction. [-]

!-----------------------------------------------------------------------------
!Local Vars
!-----------------------------------------------------------------------------
INTEGER                :: l,n,k,mclass

REAL                   ::                                                      &
npp_gb,                                                                        &
              !  NPP per gridbox area available to the intra group.
              !  [kg C m-2 s-1]
litter_gb,                                                                     &
              !  Local litterfall of the whole population per gridbox area.
              !  [kg C m-2 s-1]
w_intra_sum,                                                                   &
              !  Intra-group number density weighted by the metabolic
              !  growth scaling. [m-2]
deficit
              !  Carbon deficit of the lowest mass class per gridbox area.
              !  [kg C m-2 s-1]

!End of headers

red_state%P_a(:,:)               = 0.0
red_state%G_struct(:,:)          = 0.0
red_state%P_s(:,:)               = 0.0
red_state%npp_mass(:,:,:)        = 0.0
red_state%g_mass_intra(:,:,:)    = 0.0
red_state%g_mass_inter(:,:,:)    = 0.0
red_state%g_mass(:,:,:)          = 0.0
red_state%mort_starv(:,:,:)      = 0.0

DO l = 1,land_pts
  DO n = 1,nnpft
    mclass = red_parms%mclass(n)

    w_intra_sum = 0.0
    litter_gb   = 0.0
    DO k = 1, mclass
      w_intra_sum = w_intra_sum + red_state%plantNumDensity_intra(l,n,k)       &
                    * red_state%g_mass_scale(n,k)
      litter_gb = litter_gb + red_state%plantNumDensity(l,n,k)                 &
                  * red_state%litter_mass(l,n,k)
    END DO

    ! A positive exposed fraction implies a non-empty intra group, so no NPP
    ! is lost by this guard.
    npp_gb = 0.0
    IF (w_intra_sum > 0.0) npp_gb = npp_dr(l,n) * frac(l,n)

    ! Total gridbox growth (P_a), split into recruitment (P_s) and
    ! structural growth (G_struct): P_a = P_s + G_struct.
    red_state%P_a(l,n) = npp_gb - litter_gb
    red_state%P_s(l,n) = red_parms%alpha_recrt(n)                              &
                         * MAX(0.0, red_state%P_a(l,n))
    red_state%G_struct(l,n) = red_state%P_a(l,n) - red_state%P_s(l,n)

    DO k = 1, mclass
      IF (w_intra_sum > 0.0) THEN
        red_state%npp_mass(l,n,k) = npp_gb * red_state%g_mass_scale(n,k)       &
                                    / w_intra_sum
        red_state%g_mass_intra(l,n,k) = red_state%npp_mass(l,n,k)              &
          - red_state%litter_mass(l,n,k)                                       &
          - red_state%P_s(l,n) * red_state%g_mass_scale(n,k) / w_intra_sum
      ELSE
        red_state%g_mass_intra(l,n,k) = -red_state%litter_mass(l,n,k)
      END IF
      red_state%g_mass_inter(l,n,k) = -red_state%litter_mass(l,n,k)

      ! Population-mean individual growth, diagnostic only.
      IF (red_state%plantNumDensity(l,n,k) > 0.0) THEN
        red_state%g_mass(l,n,k) = (red_state%plantNumDensity_intra(l,n,k)      &
          * red_state%g_mass_intra(l,n,k)                                      &
          + red_state%plantNumDensity_inter(l,n,k)                             &
          * red_state%g_mass_inter(l,n,k))                                     &
          / red_state%plantNumDensity(l,n,k)
      ELSE
        red_state%g_mass(l,n,k) = red_state%g_mass_intra(l,n,k)
      END IF
    END DO

    ! Plants in the lowest class cannot shrink, so their carbon deficit is
    ! met by the death of plants (carbon starvation).
    deficit = red_state%plantNumDensity_intra(l,n,1)                           &
              * MAX(0.0, -red_state%g_mass_intra(l,n,1))                       &
              + red_state%plantNumDensity_inter(l,n,1)                         &
              * MAX(0.0, -red_state%g_mass_inter(l,n,1))
    IF (red_state%plantNumDensity(l,n,1) > 0.0) THEN
      red_state%mort_starv(l,n,1) = deficit                                    &
        / (red_state%plantNumDensity(l,n,1) * red_state%mass_cohort(n,1))
    END IF
  END DO
END DO

END SUBROUTINE veg3_red_assimilate_partition
!-----------------------------------------------------------------------------

END MODULE veg3_red_assimilate_parition_mod
