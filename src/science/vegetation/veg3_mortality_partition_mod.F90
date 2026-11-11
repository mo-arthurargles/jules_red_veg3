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

MODULE veg3_mortality_partition_mod

IMPLICIT NONE

PRIVATE
PUBLIC :: veg3_mortality_partition

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName=                            &
  'VEG3_MORTALITY_PARTITION_MOD'

CONTAINS

!-----------------------------------------------------------------------------
SUBROUTINE veg3_mortality_partition(                                           &
                !IN Control vars
                land_pts,nnpft,nmasst,                                         &
                !IN red_parms
                red_parms,                                                     &
                !IN fields
                mort_add,                                                      &
                !INOUT state
                red_state                                                      &
                )

! Aggregates mortality sources (disturbance mort_dstrb, baseline) into the
! total mortality rate (mort) across each PFT's mass classes, ahead of
! veg3_red_dynamic: mort = mort_base + mort_dstrb.

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
INTEGER, INTENT(IN) :: land_pts,nnpft,nmasst

!----------------------------------------------------------------------------
! Reals with INTENT IN
!----------------------------------------------------------------------------
REAL, INTENT(IN)   :: mort_add(land_pts,nnpft,nmasst)
              !  Additional (disturbance) plant mortality across plant mass
              !  (s-1)

!-----------------------------------------------------------------------------
!Local Vars
!-----------------------------------------------------------------------------
INTEGER                ::l,n,mclass

!End of headers

red_state%mort_dstrb(:,:,:) = 0.0
red_state%mort(:,:,:)       = 0.0

DO l = 1,land_pts
  DO n = 1,nnpft
    mclass = red_parms%mclass(n)

    red_state%mort_dstrb(l,n,1:mclass) = mort_add(l,n,1:mclass)
    red_state%mort(l,n,1:mclass) = red_parms%mort_base(n)                      &
      + red_state%mort_dstrb(l,n,1:mclass)
  END DO
END DO

END SUBROUTINE veg3_mortality_partition
!-----------------------------------------------------------------------------

END MODULE veg3_mortality_partition_mod
