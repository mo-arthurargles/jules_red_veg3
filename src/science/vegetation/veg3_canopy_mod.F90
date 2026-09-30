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
! Description:
!   Closed-canopy crown overlap scheme for veg3/RED. Replaces the naive
!   frac=SUM(plantNumDensity*crwn_area_mass) calculation with a
!   height-ordered gap-fraction correction based on the crown area index
!   (CAI), so that taller cohorts/PFTs shade and reduce the ground area
!   available to shorter ones.

MODULE veg3_canopy_mod

IMPLICIT NONE

PRIVATE
PUBLIC :: veg3_canopy_height_order, veg3_canopy_frac,                          &
          veg3_canopy_maintain_frac_min

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='VEG3_CANOPY_MOD'

CONTAINS

!-------------------------------------------------------------------------------
SUBROUTINE veg3_canopy_height_order(nnpft, nmasst, npft_totmclass, mclass,     &
                                     ht_mass, order_pft_desc,                  &
                                     order_mclass_desc, order_pft_asc,         &
                                     order_mclass_asc)
! Flattens the (PFT, mass-class) pairs and sorts them by static allometric
! height, tallest-to-shortest (desc) and shortest-to-tallest (asc). Computed
! once at initialisation, since ht_mass does not vary with plantNumDensity.
! Unused trailing slots (beyond the number of valid PFT/mass-class pairs)
! are zero-padded; consumers should stop at the first zero PFT index.

IMPLICIT NONE

INTEGER, INTENT(IN) :: nnpft, nmasst, npft_totmclass
INTEGER, INTENT(IN) :: mclass(nnpft)
REAL, INTENT(IN)    :: ht_mass(nnpft,nmasst)

INTEGER, INTENT(OUT) :: order_pft_desc(npft_totmclass)
INTEGER, INTENT(OUT) :: order_mclass_desc(npft_totmclass)
INTEGER, INTENT(OUT) :: order_pft_asc(npft_totmclass)
INTEGER, INTENT(OUT) :: order_mclass_asc(npft_totmclass)

INTEGER :: n, k, i, j, nvalid
INTEGER :: pft_tmp, mclass_tmp
REAL    :: ht_tmp

!End of header

order_pft_desc(:)    = 0
order_mclass_desc(:) = 0
order_pft_asc(:)     = 0
order_mclass_asc(:)  = 0

! Flatten the valid (n,k) pairs into the ascending order arrays.
nvalid = 0
DO n = 1, nnpft
  DO k = 1, mclass(n)
    nvalid = nvalid + 1
    order_pft_asc(nvalid)    = n
    order_mclass_asc(nvalid) = k
  END DO
END DO

! Simple insertion sort by height, ascending (nvalid is small: nnpft*nmasst).
DO i = 2, nvalid
  pft_tmp    = order_pft_asc(i)
  mclass_tmp = order_mclass_asc(i)
  ht_tmp     = ht_mass(pft_tmp,mclass_tmp)
  j = i - 1
  DO WHILE (j >= 1)
    IF (ht_mass(order_pft_asc(j),order_mclass_asc(j)) <= ht_tmp) EXIT
    order_pft_asc(j+1)    = order_pft_asc(j)
    order_mclass_asc(j+1) = order_mclass_asc(j)
    j = j - 1
  END DO
  order_pft_asc(j+1)    = pft_tmp
  order_mclass_asc(j+1) = mclass_tmp
END DO

! Descending order is simply the reverse of the valid ascending entries.
DO i = 1, nvalid
  order_pft_desc(i)    = order_pft_asc(nvalid - i + 1)
  order_mclass_desc(i) = order_mclass_asc(nvalid - i + 1)
END DO

RETURN
END SUBROUTINE veg3_canopy_height_order

!-------------------------------------------------------------------------------
SUBROUTINE veg3_canopy_frac(land_pts, nnpft, nmasst, npft_totmclass,           &
                             order_pft_desc, order_mclass_desc, k_cai,         &
                             frac_excl, frac_tile_excl, crwn_area_mass,        &
                             plantNumDensity, frac_mass, CAI_mass,             &
                             CAI_overlapped, frac, CAI, frac_above_mclass1)
! Height-ordered gap-fraction overlap calculation.
! Processes all (PFT, mass-class) cohorts tallest-to-shortest at each land
! point, so each cohort's fraction is only allowed to claim ground not
! already claimed by taller cohorts (frac_above) or reserved for other
! land cover/competing PFTs (frac_tile_excl/frac_excl).

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, nnpft, nmasst, npft_totmclass
INTEGER, INTENT(IN) :: order_pft_desc(npft_totmclass)
INTEGER, INTENT(IN) :: order_mclass_desc(npft_totmclass)
REAL, INTENT(IN)    :: k_cai(nnpft)
              ! PFT crown overlap coefficient. (-)
REAL, INTENT(IN)    :: frac_excl(nnpft)
              ! Fraction reserved for other PFTs' own frac_min. (-)
REAL, INTENT(IN)    :: frac_tile_excl(land_pts)
              ! Fraction reserved for non-vegetated tiles (urban/lake/ice). (-)
REAL, INTENT(IN)    :: crwn_area_mass(nnpft,nmasst)
REAL, INTENT(IN)    :: plantNumDensity(land_pts,nnpft,nmasst)

REAL, INTENT(OUT) :: frac_mass(land_pts,nnpft,nmasst)
REAL, INTENT(OUT) :: CAI_mass(land_pts,nnpft,nmasst)
REAL, INTENT(OUT) :: CAI_overlapped(land_pts,nnpft,nmasst)
REAL, INTENT(OUT) :: frac(land_pts,nnpft)
REAL, INTENT(OUT) :: CAI(land_pts,nnpft)
REAL, INTENT(OUT) :: frac_above_mclass1(land_pts,nnpft)
              ! frac_above at the moment the PFT's lowest mass class (k=1)
              ! was processed. Used by veg3_canopy_maintain_frac_min.

INTEGER :: l, i, n, k
REAL    :: frac_above, frac_open

!End of header

frac(:,:)               = 0.0
CAI(:,:)                = 0.0
frac_mass(:,:,:)        = 0.0
CAI_mass(:,:,:)         = 0.0
CAI_overlapped(:,:,:)   = 0.0
frac_above_mclass1(:,:) = 0.0

DO l = 1, land_pts
  frac_above = 0.0

  DO i = 1, npft_totmclass
    n = order_pft_desc(i)
    k = order_mclass_desc(i)
    IF (n == 0) EXIT ! Zero indicates unused memory element in ordered arrays.

    IF (k == 1) frac_above_mclass1(l,n) = frac_above

    ! EQ1: crown area index of this cohort.
    CAI_mass(l,n,k) = plantNumDensity(l,n,k) * crwn_area_mass(n,k)

    ! EQ2: gap-fraction overlap, restricted to the ground not already
    ! claimed by taller cohorts or excluded tiles/competitor PFTs.
    frac_open = MAX(0.0, 1.0 - frac_above - frac_tile_excl(l) - frac_excl(n))
    frac_mass(l,n,k) = frac_open * (1.0 - EXP(-k_cai(n) * CAI_mass(l,n,k)))

    ! EQ3: overlap-corrected CAI per unit fraction claimed.
    IF (CAI_mass(l,n,k) > 0.0) THEN
      CAI_overlapped(l,n,k) = frac_mass(l,n,k) / CAI_mass(l,n,k)
    END IF

    frac(l,n) = frac(l,n) + frac_mass(l,n,k)
    CAI(l,n)  = CAI(l,n)  + CAI_mass(l,n,k)
    frac_above = frac_above + frac_mass(l,n,k)
  END DO
END DO

RETURN
END SUBROUTINE veg3_canopy_frac

!-------------------------------------------------------------------------------
SUBROUTINE veg3_canopy_maintain_frac_min(land_pts, nnpft, frac_min, k_cai,     &
                                          frac_excl, frac_tile_excl,           &
                                          crwn_area_mass_1, mass_mass_1,       &
                                          frac_above_mclass1, frac,            &
                                          plantNumDensity_1, dt, mort_litC)
! Maintains each PFT's minimum vegetation fraction by adding recruits to its
! lowest mass class, inverting the to find the plantNumDensity needed. If
! 
! When dt/mort_litC are supplied (runtime use), the implied recruitment is
! debited from mort_litC to conserve carbon; omitted at initialisation.

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, nnpft
REAL, INTENT(IN) :: frac_min(nnpft), k_cai(nnpft), frac_excl(nnpft)
REAL, INTENT(IN) :: frac_tile_excl(land_pts)
REAL, INTENT(IN) :: crwn_area_mass_1(nnpft), mass_mass_1(nnpft)
REAL, INTENT(IN) :: frac_above_mclass1(land_pts,nnpft)

REAL, INTENT(IN OUT) :: frac(land_pts,nnpft)
REAL, INTENT(IN OUT) :: plantNumDensity_1(land_pts,nnpft)

REAL, INTENT(IN), OPTIONAL :: dt
REAL, INTENT(IN OUT), OPTIONAL :: mort_litC(land_pts,nnpft)

INTEGER :: l, n
REAL    :: frac_open, cai_new, plantnumdensity_new, delta_pnd

!End of header

DO n = 1, nnpft
  DO l = 1, land_pts
    IF (frac(l,n) < frac_min(n)) THEN
      frac_open = 1.0 - frac_above_mclass1(l,n) - frac_tile_excl(l)           &
                  - frac_excl(n)

      IF (frac_open > frac_min(n)) THEN
        ! Algebraic inverse of EQ2 for the lowest mass class.
        cai_new = -LOG( (frac_open - frac_min(n)) / frac_open ) / k_cai(n)
        plantnumdensity_new = cai_new / crwn_area_mass_1(n)

        IF (PRESENT(mort_litC) .AND. PRESENT(dt)) THEN
          delta_pnd = MAX(0.0, plantnumdensity_new - plantNumDensity_1(l,n))
          mort_litC(l,n) = mort_litC(l,n) - delta_pnd * mass_mass_1(n) / dt
        END IF

        plantNumDensity_1(l,n) = MAX(plantNumDensity_1(l,n),                  &
                                      plantnumdensity_new)
        frac(l,n) = frac_min(n)
      END IF
    END IF
  END DO
END DO

RETURN
END SUBROUTINE veg3_canopy_maintain_frac_min

END MODULE veg3_canopy_mod
