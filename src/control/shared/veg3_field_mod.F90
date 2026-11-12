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

MODULE veg3_field_mod

! Description:
!   Holds the veg3/RED vegetation state variables/fields (veg_state,
!   red_state) and provides routines to allocate, deallocate, associate and
!   initialise them.
!
!   Note: the soil biogeochemistry state (soil_state) previously held here
!   now lives in soil_bgc_4pool_field_mod.F90, and the veg3/RED-soil
!   carbon coupling routine (veg3_soil_couple) now lives in
!   soil_bgc_4pool_control_mod.F90. This keeps soil biogeochemistry
!   modular and independent of the choice of vegetation dynamics model.
!   Its allocation/association is therefore called separately from
!   init_mod.F90, not from within this module.

USE veg3_parm_mod, ONLY: red_parms
USE um_types, ONLY: real_jlslsm
USE jules_vegetation_mod, ONLY: l_red

!Use at module level
USE ancil_info,    ONLY: ainfo_type
USE prognostics,   ONLY: progs_type

IMPLICIT NONE

! Structure to keep the vegetation state variables and fields
TYPE :: veg_state_type
  REAL, POINTER ::                                                             &
      leafC(:,:),                                                              &
              ! PFT leaf carbon per PFT area fraction. (kg C m-2)
      rootC(:,:),                                                              &
              ! PFT root carbon per PFT area fraction. (kg C m-2)
      woodC(:,:),                                                              &
              ! PFT woody carbon per PFT area fraction. (kg C m-2)
      lai_bal(:,:),                                                            &
              ! Balanced LAI, per unit open (overlap-corrected,
              ! frac_mass-weighted) PFT ground fraction. (m2 m-2)
      vegCpft(:,:),                                                            &
              ! Total PFT carbon density per PFT area fraction. (kg C m-2)
      vegC(:),                                                                 &
              ! Gridbox mean vegetation carbon. (kg C m-2)
      npp_gb(:),                                                               &
              ! Gridbox mean NPP driving RED. (kg C m-2 s-1)
      npp_n_gb(:),                                                             &
              ! Gridbox mean NPP after nitrogen limitation. Nitrogen is not
              ! yet coupled to veg3/RED, so this is npp_gb converted to
              ! (360d)-1 units, matching trif_vars_data%npp_n_gb.
              ! (kg C m-2 (360d)-1)
      npp_acc(:,:),                                                            &
              ! Accumulated NPP. (kg C m-2 s-1)
      npp_dr_out(:,:),                                                         &
              ! A diagnostic NPP variable driving RED. (kg C m-2 (360d)-1)
      frac(:,:),                                                               &
              ! Gridbox area fraction of each tile type including PFTs.
      phen(:,:),                                                               &
              ! PFT phenology state variable, diagnosed from lai and lai_bal.
      lai(:,:),                                                                &
              ! Leaf area index, per unit open (overlap-corrected,
              ! frac_mass-weighted) PFT ground fraction. Associated with
              ! trif_vars%lai_pft (m2 m-2).
      canht(:,:),                                                              &
              ! Canopy height. Associated with trif_vars%canht_pft. (m)
      g_leaf_phen(:,:),                                                        &
              ! PFT leaf carbon mass density turnover. ((360d)-1)
      g_leaf(:,:),                                                             &
              ! Pointer necessary to get the diagnosed g_leaf rate from the
              ! rest of JULES. ((360d)-1)
      g_leaf_day(:,:),                                                         &
              ! PFT mean leaf turnover rate for input to phenology. ((360d)-1)
      g_leaf_dr_out(:,:),                                                      &
              ! PFT mean leaf turnover rate for driving vegetation dynamics.
              ! ((360d)-1)
      lai_phen(:,:),                                                           &
              ! Diagnostic LAI immediately following the phenology update.
              ! Associated with trifctl%lai_phen_pft. (m2 m-2)
      g_leaf_acc(:,:),                                                         &
              ! Accumulated leaf turnover rate since the last phenology call.
              ! ((360d)-1)
      g_leaf_phen_acc(:,:),                                                    &
              ! Accumulated mean phenological leaf turnover rate since the
              ! last vegetation dynamics call. ((360d)-1)
      leaf_litC(:,:),                                                          &
              ! Leaf litter carbon flux per PFT fraction. (kg C m-2 (360d)-1)
      root_litC(:,:),                                                          &
              ! Root litter carbon flux per PFT fraction. (kg C m-2 (360d)-1)
      wood_litC(:,:),                                                          &
              ! Wood litter carbon flux per PFT fraction. (kg C m-2 (360d)-1)
      litCpft(:,:),                                                            &
              ! Total litter carbon flux per PFT fraction, also includes any
              ! litter from dynamics (e.g., mortality). (kg C m-2 (360d)-1)
      litC(:)
              ! Total litter carbon flux per gridbox. (kg C m-2 (360d)-1)

  REAL, ALLOCATABLE ::                                                         &
    mort_litC(:,:),                                                            &
            ! Mortality/demographic litter carbon flux from vegetation
            ! dynamics, normalised per unit PFT canopy area (kg C m-2 s-1).
    frac_tile_excl(:),                                                         &
            ! Fraction of the gridbox excluded from vegetation (non-soil,
            ! non-PFT tiles: urban/lake/ice). Static. (-)
    cai(:,:),                                                                  &
            ! PFT crown area index (CAI) (m2 m-2)
    vegCpft_CAI(:,:),                                                          &
            ! Total PFT carbon density per unit crown area index (CAI).
            ! (kg C m-2)
    leafC_CAI(:,:),                                                            &
            ! Open-canopy PFT leaf carbon per unit crown area index
            ! (CAI). (kg C m-2)
    rootC_CAI(:,:),                                                            &
            ! Open-canopy PFT root carbon per unit crown area index
            ! (CAI). (kg C m-2)
    woodC_CAI(:,:)
            ! Open-canopy PFT woody carbon per unit crown area index
            ! (CAI). (kg C m-2)

END TYPE veg_state_type

! Structure to keep the RED state variables and fields
TYPE :: red_state_type
  REAL, ALLOCATABLE ::                                                         &
    mass_cohort(:,:),                                                          &
              !  PFT plant mass tiles/classes. (kg C)
    ht_mass(:,:),                                                              &
              !  PFT height across plant mass. (m)
    lai_bal_mass(:,:),                                                         &
              !  PFT balanced leaf area index across plant mass. (m2/m2)
    crwn_area_mass(:,:),                                                       &
              !  PFT crown area across plant mass. (m2)
    g_mass_scale(:,:),                                                         &
              !  PFT plant growth scaling wrt metabolic scaling
              !  theory across mass. (kg C /kg C)
    leaf_mass(:,:),                                                            &
              !  Individual leaf carbon mass, by mass class. (kg C)
    froot_mass(:,:),                                                           &
              !  Individual fine root carbon mass, by mass class. (kg C)
    wood_mass(:,:),                                                            &
              !  Individual wood carbon mass, by mass class. (kg C)
    mclass_geom_mult(:),                                                       &
              !  PFT geometric scaling coefficent for binning mass classes (-)
    frac_mass(:,:,:),                                                          &
              !  PFT closed-canopy overlap-corrected area fraction, by mass
              !  class. (-)
    CAI_mass(:,:,:),                                                           &
              !  PFT crown area index, by mass class. (m2 m-2)
    mort(:,:,:),                                                               &
              !  Total PFT mortality rate across plant mass (mort_base +
              !  mort_dstrb), precomputed by veg3_mortality_partition_mod.
              !  (/s)
    mort_dstrb(:,:,:),                                                         &
              !  Disturbance PFT mortality rate across plant mass. (/s)
    g_mass(:,:,:),                                                             &
              !  Individual growth across plant mass, precomputed by
              !  veg3_red_assimilate_parition_mod. (kgC s-1)
    P_s(:,:),                                                                  &
              !  Gridbox carbon assimilate devoted to recruitment, precomputed
              !  by veg3_red_assimilate_parition_mod. (kgC m-2 s-1)
    P_a(:,:),                                                                  &
              !  Total gridbox carbon assimilate (P_a = P_s + G_struct).
              !  (kgC m-2 s-1)
    G_struct(:,:)
              !  Gridbox carbon assimilate devoted to structural growth
              !  (P_a = P_s + G_struct). (kgC m-2 s-1)
  REAL, POINTER ::                                                             &
    plantNumDensity(:,:,:)
              !  PFT number density across plant mass. (/m2)
  INTEGER, ALLOCATABLE ::                                                      &
    order_pft_desc(:), order_mclass_desc(:),                                   &
              !  Flattened (PFT, mass-class) pairs ordered tallest-to-shortest
              !  by static allometric height.
    order_pft_asc(:), order_mclass_asc(:)
              !  As above, ordered shortest-to-tallest.
END TYPE red_state_type

TYPE(veg_state_type)   :: veg_state
TYPE(red_state_type)   :: red_state

!Private by default
PRIVATE

!Expose routines
PUBLIC :: veg3_field_init, veg3_field_allocate, veg3_field_deallocate,         &
          veg3_field_assoc, red_veg3_couple

!Expose data
PUBLIC :: veg_state, red_state

!Expose data structures
PUBLIC :: veg_state_type, red_state_type

!Allow external code to read but not write
!PROTECTED ::

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='VEG3_FIELD_MOD'

CONTAINS
!-------------------------------------------------------------------------------

SUBROUTINE veg3_field_allocate(land_pts,nsurft,nnpft,nmasst,npft_totmclass)

IMPLICIT NONE
INTEGER, INTENT(IN) :: land_pts, nsurft, nnpft, nmasst, npft_totmclass

!End of Header

ALLOCATE(veg_state%leafC        ( land_pts, nnpft) )
ALLOCATE(veg_state%rootC        ( land_pts, nnpft) )
ALLOCATE(veg_state%woodC        ( land_pts, nnpft) )
ALLOCATE(veg_state%vegCpft      ( land_pts, nnpft) )
ALLOCATE(veg_state%lai_bal      ( land_pts, nnpft) )
ALLOCATE(veg_state%canht        ( land_pts, nnpft) )
ALLOCATE(veg_state%lai          ( land_pts, nnpft) )
ALLOCATE(veg_state%phen         ( land_pts, nnpft) )
ALLOCATE(veg_state%npp_acc      ( land_pts, nnpft) )
ALLOCATE(veg_state%npp_dr_out   ( land_pts, nnpft) )
ALLOCATE(veg_state%g_leaf_phen  ( land_pts, nnpft) )
ALLOCATE(veg_state%g_leaf       ( land_pts, nnpft) )
ALLOCATE(veg_state%g_leaf_day   ( land_pts, nnpft) )
ALLOCATE(veg_state%g_leaf_dr_out( land_pts, nnpft) )
ALLOCATE(veg_state%lai_phen     ( land_pts, nnpft) )
ALLOCATE(veg_state%g_leaf_acc   ( land_pts, nnpft) )
ALLOCATE(veg_state%g_leaf_phen_acc ( land_pts, nnpft) )
ALLOCATE(veg_state%frac         ( land_pts, nsurft) )
ALLOCATE(veg_state%vegC         ( land_pts) )
ALLOCATE(veg_state%npp_gb       ( land_pts) )
ALLOCATE(veg_state%npp_n_gb     ( land_pts) )
ALLOCATE(veg_state%leaf_litC    ( land_pts, nnpft) )
ALLOCATE(veg_state%root_litC    ( land_pts, nnpft) )
ALLOCATE(veg_state%wood_litC    ( land_pts, nnpft) )
ALLOCATE(veg_state%litCpft     ( land_pts, nnpft) )
ALLOCATE(veg_state%litC         ( land_pts) )
ALLOCATE(veg_state%mort_litC    ( land_pts, nnpft) )
ALLOCATE(veg_state%cai          ( land_pts, nnpft) )
ALLOCATE(veg_state%vegCpft_CAI  ( land_pts, nnpft) )
ALLOCATE(veg_state%leafC_CAI    ( land_pts, nnpft) )
ALLOCATE(veg_state%rootC_CAI    ( land_pts, nnpft) )
ALLOCATE(veg_state%woodC_CAI    ( land_pts, nnpft) )
ALLOCATE(veg_state%frac_tile_excl( land_pts) )

!Initialise
veg_state%leafC(:,:)           = 0.0
veg_state%rootC(:,:)           = 0.0
veg_state%woodC(:,:)           = 0.0
veg_state%vegCpft(:,:)         = 0.0
veg_state%lai_bal(:,:)         = 0.0
veg_state%canht(:,:)           = 0.0
veg_state%lai(:,:)             = 0.0
veg_state%phen(:,:)            = 0.0
veg_state%npp_acc(:,:)         = 0.0
veg_state%npp_dr_out(:,:)      = 0.0
veg_state%g_leaf_phen(:,:)     = 0.0
veg_state%g_leaf(:,:)          = 0.0
veg_state%g_leaf_day(:,:)      = 0.0
veg_state%g_leaf_dr_out(:,:)   = 0.0
veg_state%lai_phen(:,:)        = 0.0
veg_state%g_leaf_acc(:,:)      = 0.0
veg_state%g_leaf_phen_acc(:,:) = 0.0
veg_state%frac(:,:)            = 0.0
veg_state%vegC(:)              = 0.0
veg_state%npp_gb(:)            = 0.0
veg_state%npp_n_gb(:)          = 0.0
veg_state%leaf_litC(:,:)       = 0.0
veg_state%root_litC(:,:)       = 0.0
veg_state%wood_litC(:,:)       = 0.0
veg_state%litCpft(:,:)         = 0.0
veg_state%litC(:)              = 0.0
veg_state%mort_litC(:,:)       = 0.0
veg_state%cai(:,:)             = 0.0
veg_state%vegCpft_CAI(:,:)     = 0.0
veg_state%leafC_CAI(:,:)       = 0.0
veg_state%rootC_CAI(:,:)       = 0.0
veg_state%woodC_CAI(:,:)       = 0.0
veg_state%frac_tile_excl(:)    = 0.0

! RED

! Allocate red_data_type
ALLOCATE(red_state%mass_cohort          (nnpft, nmasst ))
ALLOCATE(red_state%ht_mass            (nnpft, nmasst ))
ALLOCATE(red_state%lai_bal_mass       (nnpft, nmasst ))
ALLOCATE(red_state%crwn_area_mass     (nnpft, nmasst ))
ALLOCATE(red_state%g_mass_scale       (nnpft, nmasst ))
ALLOCATE(red_state%leaf_mass           (nnpft, nmasst ))
ALLOCATE(red_state%froot_mass          (nnpft, nmasst ))
ALLOCATE(red_state%wood_mass           (nnpft, nmasst ))
ALLOCATE(red_state%plantNumDensity    (land_pts, nnpft, nmasst ))
ALLOCATE(red_state%mort               (land_pts, nnpft, nmasst ))
ALLOCATE(red_state%mort_dstrb         (land_pts, nnpft, nmasst ))
ALLOCATE(red_state%frac_mass          (land_pts, nnpft, nmasst ))
ALLOCATE(red_state%CAI_mass           (land_pts, nnpft, nmasst ))
ALLOCATE(red_state%g_mass             (land_pts, nnpft, nmasst ))
ALLOCATE(red_state%P_s                (land_pts, nnpft ))
ALLOCATE(red_state%P_a                 (land_pts, nnpft ))
ALLOCATE(red_state%G_struct            (land_pts, nnpft ))
ALLOCATE(red_state%order_pft_desc     (npft_totmclass))
ALLOCATE(red_state%order_mclass_desc  (npft_totmclass))
ALLOCATE(red_state%order_pft_asc      (npft_totmclass))
ALLOCATE(red_state%order_mclass_asc   (npft_totmclass))

! Initialise red_data_type
red_state%mass_cohort(:,:)          = 0.0
red_state%ht_mass(:,:)            = 0.0
red_state%lai_bal_mass(:,:)       = 0.0
red_state%crwn_area_mass(:,:)     = 0.0
red_state%g_mass_scale(:,:)       = 0.0
red_state%leaf_mass(:,:)           = 0.0
red_state%froot_mass(:,:)          = 0.0
red_state%wood_mass(:,:)           = 0.0
red_state%plantNumDensity(:,:,:)  = 0.0
red_state%mort(:,:,:)             = 0.0
red_state%mort_dstrb(:,:,:)       = 0.0
red_state%frac_mass(:,:,:)        = 0.0
red_state%CAI_mass(:,:,:)         = 0.0
red_state%g_mass(:,:,:)           = 0.0
red_state%P_s(:,:)                = 0.0
red_state%P_a(:,:)                = 0.0
red_state%G_struct(:,:)           = 0.0
red_state%order_pft_desc(:)       = 0
red_state%order_mclass_desc(:)    = 0
red_state%order_pft_asc(:)        = 0
red_state%order_mclass_asc(:)     = 0

RETURN
END SUBROUTINE veg3_field_allocate

!-------------------------------------------------------------------------------

SUBROUTINE veg3_field_deallocate()

! Deallocates the veg_state and red_state POINTER components that were
! allocated locally by veg3_field_allocate but which veg3_field_assoc later
! re-associates onto external targets (e.g. progs, trif_vars_data and
! trifctl_data fields). This must be called before that re-association is
! done, otherwise the original local allocations become orphaned (leaked) -
! this can be a significant leak for red_state%plantNumDensity in particular.
!
! Note: veg_state%phen, veg_state%npp_acc, veg_state%mort_litC and the
! red_state ALLOCATABLE components (mass_cohort, ht_mass, lai_bal_mass,
! crwn_area_mass, g_mass_scale, leaf_mass, froot_mass, wood_mass, mort) are
! deliberately NOT included here. None of these are re-associated in
! veg3_field_assoc (the red_state ones are ALLOCATABLE rather than
! POINTER, so cannot be re-associated at all) and they remain the arrays
! allocated in veg3_field_allocate for the lifetime of the run. They are
! read/written every timestep (e.g. in red_veg3_couple and
! veg3_red_dynamic_mod), so deallocating them here - before that ongoing
! use - would break the run.

IMPLICIT NONE

!End of header

DEALLOCATE(veg_state%leafC)
DEALLOCATE(veg_state%rootC)
DEALLOCATE(veg_state%woodC)
DEALLOCATE(veg_state%vegCpft)
DEALLOCATE(veg_state%lai_bal)
DEALLOCATE(veg_state%canht)
DEALLOCATE(veg_state%lai)
DEALLOCATE(veg_state%npp_dr_out)
DEALLOCATE(veg_state%frac)
DEALLOCATE(veg_state%g_leaf_phen)
DEALLOCATE(veg_state%g_leaf)
DEALLOCATE(veg_state%g_leaf_day)
DEALLOCATE(veg_state%g_leaf_dr_out)
DEALLOCATE(veg_state%lai_phen)
DEALLOCATE(veg_state%g_leaf_acc)
DEALLOCATE(veg_state%g_leaf_phen_acc)
DEALLOCATE(veg_state%leaf_litC)
DEALLOCATE(veg_state%root_litC)
DEALLOCATE(veg_state%wood_litC)
DEALLOCATE(veg_state%litCpft)
DEALLOCATE(veg_state%litC)
DEALLOCATE(veg_state%vegC)
DEALLOCATE(veg_state%npp_gb)
DEALLOCATE(veg_state%npp_n_gb)

DEALLOCATE(red_state%plantNumDensity)

RETURN
END SUBROUTINE veg3_field_deallocate

!-------------------------------------------------------------------------------
SUBROUTINE veg3_field_assoc(progs, ainfo, trifctl_data, trif_vars_data)

! Initial code to associate the veg3 and red fields to the rest of JULES
! This new routine moves out the pointers from veg3_set_fields to here
! to keep clean and more in line with in init.F90.

! We still need to run with l_triffid, but this will need to be addressed in a
! future revision to allow for veg3 to run fully independently of the switch.

! Note: trifctl_data and trif_vars_data are passed in as arguments (rather than
! USE-associated from jules_fields_mod) because this module lives in
! src/control/shared and must also build for the UM, where jules_fields_mod
! (a standalone-only module) is not available.

USE jules_vegetation_mod,     ONLY: l_triffid
USE trifctl,                  ONLY: trifctl_data_type
USE trif_vars_mod,            ONLY: trif_vars_data_type

IMPLICIT NONE

TYPE(progs_type), INTENT(IN) :: progs
TYPE(ainfo_type), INTENT(IN) :: ainfo
TYPE(trifctl_data_type), INTENT(IN), TARGET :: trifctl_data
TYPE(trif_vars_data_type), INTENT(IN), TARGET :: trif_vars_data
! End of header
!-------------------------------------------------------------------------------
IF (l_red .AND. l_triffid) THEN

  ! Deallocate the local pointer targets set up in veg3_field_allocate before
  ! they are re-associated onto external targets below, to avoid leaking the
  ! original allocations (see veg3_field_deallocate for details).
  CALL veg3_field_deallocate()

  ! Set pointers to the prognostic fields
  ! Note: veg_state%phen is not associated to an external target here - it is
  ! not an independent prognostic. It is diagnosed each phenology call as
  ! lai / lai_bal (see veg3_phenol_couple in next_gen_biogeochem_mod), so it
  ! remains the array allocated locally in veg3_field_allocate.
  red_state%plantNumDensity => progs%plantNumDensity
  veg_state%frac        => ainfo%frac_surft
  veg_state%canht       => progs%canht_pft
  veg_state%lai         => progs%lai_pft

  ! Additional TRIFFID fields mimicking the pointer associations set up in
  ! trif_vars_assoc and trifctl_assoc, so that veg3/RED diagnostics feed
  ! straight back into the standard JULES carbon diagnostics and outputs.

  ! Firstly trif_vars_data
  veg_state%lai_bal => trif_vars_data%lai_bal_pft
  veg_state%leafC   => trif_vars_data%leafc_pft
  veg_state%rootC   => trif_vars_data%rootc_pft
  veg_state%woodC   => trif_vars_data%woodc_pft
  veg_state%leaf_litC => trif_vars_data%leaf_litc_pft
  veg_state%root_litC => trif_vars_data%root_litc_pft
  veg_state%wood_litC => trif_vars_data%wood_litc_pft

  ! Next trifctl_data
  veg_state%npp_dr_out => trifctl_data%npp_dr_out_pft
  veg_state%g_leaf => trifctl_data%g_leaf_pft
  veg_state%g_leaf_phen => trifctl_data%g_leaf_phen_pft
  veg_state%g_leaf_day => trifctl_data%g_leaf_day_pft
  veg_state%g_leaf_dr_out => trifctl_data%g_leaf_dr_out_pft
  veg_state%lai_phen => trifctl_data%lai_phen_pft
  veg_state%g_leaf_acc => trifctl_data%g_leaf_acc_pft
  veg_state%g_leaf_phen_acc => trifctl_data%g_leaf_phen_acc_pft
  veg_state%vegCpft => trifctl_data%c_veg_pft
  veg_state%litCpft => trifctl_data%lit_c_pft
  veg_state%vegC => trifctl_data%cv_gb
  veg_state%litC => trifctl_data%lit_c_mn_gb
  veg_state%npp_gb => trifctl_data%npp_gb
  veg_state%npp_n_gb => trif_vars_data%npp_n_gb

END IF

RETURN

END SUBROUTINE veg3_field_assoc


!-------------------------------------------------------------------------------
SUBROUTINE veg3_set_fields(land_pts,nsurft,nnpft,nmasst,ainfo,progs)

!Source parms, etc from io modules - these should be available on and offline

! Above only allocated if triffid on - needs to be addressed
USE jules_vegetation_mod,     ONLY: l_triffid, triffid_period
USE conversions_mod,          ONLY: rsec_per_day

USE jules_surface_types_mod,  ONLY: soil, urban, urban_canyon, urban_roof,     &
                                    lake, ice
USE veg3_canopy_mod,          ONLY: veg3_canopy_frac, veg3_canopy_maintain_frac_min
USE veg3_parm_mod,            ONLY: veg3_ctrl
USE pftparm,                   ONLY: lma

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, nsurft, nnpft, nmasst

TYPE(ainfo_type), INTENT(IN OUT) :: ainfo
TYPE(progs_type), INTENT(IN) :: progs

INTEGER :: l,n

REAL :: frac_above_mclass1(land_pts,nnpft)
        ! frac_above at the moment each PFT's lowest mass class was
        ! processed (see veg3_canopy_frac).
REAL :: frac_min_init(nnpft)
        ! Minimum fraction to apply at initialisation: frac_seed on a
        ! genuine cold start (no plants), else frac_min on restart.

!End of header
!-----------------------------------------------------------------------------

IF (l_red .AND. l_triffid) THEN
  veg_state%phen(:,:) = 1.0

  ! Static fraction excluded from vegetation (non-soil, non-PFT tiles).
  veg_state%frac_tile_excl(:) = 0.0
  IF (urban        > 0) veg_state%frac_tile_excl(:) =                          &
    veg_state%frac_tile_excl(:) + ainfo%frac_surft(:,urban)
  IF (urban_canyon  > 0) veg_state%frac_tile_excl(:) =                         &
    veg_state%frac_tile_excl(:) + ainfo%frac_surft(:,urban_canyon)
  IF (urban_roof    > 0) veg_state%frac_tile_excl(:) =                         &
    veg_state%frac_tile_excl(:) + ainfo%frac_surft(:,urban_roof)
  IF (lake          > 0) veg_state%frac_tile_excl(:) =                         &
    veg_state%frac_tile_excl(:) + ainfo%frac_surft(:,lake)
  IF (ice           > 0) veg_state%frac_tile_excl(:) =                         &
    veg_state%frac_tile_excl(:) + ainfo%frac_surft(:,ice)

  ! Preliminary overlap-corrected fraction from the dump-read
  ! plantNumDensity (may still be below frac_min at this point).
  CALL veg3_canopy_frac(land_pts, nnpft, nmasst, veg3_ctrl%npft_totmclass,     &
    red_state%order_pft_desc, red_state%order_mclass_desc, red_parms%k_cai,    &
    red_parms%frac_excl, veg_state%frac_tile_excl, red_state%crwn_area_mass,   &
    red_state%plantNumDensity, red_state%frac_mass, red_state%CAI_mass,        &
    veg_state%frac(:,1:nnpft), veg_state%cai, frac_above_mclass1)

  ! Top up any PFT below its minimum fraction (no litter adjustment at
  ! initialisation - see item 9). Use frac_seed on a genuine cold start
  ! (no plants present); otherwise maintain the restarted frac_min.
  IF (SUM(progs%plantNumDensity) <= 0.0) THEN
    frac_min_init(:) = red_parms%frac_seed(:)
  ELSE
    frac_min_init(:) = red_parms%frac_min(:)
  END IF

  CALL veg3_canopy_maintain_frac_min(land_pts, nnpft, frac_min_init,           &
    red_parms%k_cai, red_parms%frac_excl, veg_state%frac_tile_excl,            &
    red_state%crwn_area_mass(:,1), red_state%mass_cohort(:,1),                 &
    frac_above_mclass1, veg_state%frac(:,1:nnpft),                             &
    red_state%plantNumDensity(:,:,1))

  DO l = 1, land_pts
    DO n = 1, nnpft

      ! red_veg3_couple needs to know the phen state variable, which is
      ! diagnosed from lai_bal. Both are not prognostic variables, but lai
      ! is, so we can infer the phen state from this and plantNumDensity
      ! for lai_bal.

      CALL pft_mean_from_mass_class(                                           &
        !IN sizing
        red_parms%mclass(n),                                                   &
        !IN mass-cohort properties
        red_state%plantNumDensity(l,n,1:red_parms%mclass(n)),                  &
        red_state%mass_cohort(n,1:red_parms%mclass(n)),                        &
        red_state%lai_bal_mass(n,1:red_parms%mclass(n)),                       &
        red_state%ht_mass(n,1:red_parms%mclass(n)),                            &
        red_state%leaf_mass(n,1:red_parms%mclass(n)),                          &
        red_state%froot_mass(n,1:red_parms%mclass(n)),                         &
        red_state%wood_mass(n,1:red_parms%mclass(n)),                          &
        red_state%crwn_area_mass(n,1:red_parms%mclass(n)),                     &
        !IN PFT trait parameter used to invert leafC for lai_bal
        lma(n),                                                                &
        !IN overlap-corrected PFT ground area fraction, by mass class
        red_state%frac_mass(l,n,1:red_parms%mclass(n)),                        &
        !IN total PFT fraction and crown area index
        veg_state%frac(l,n),veg_state%cai(l,n),                                &
        !OUT fields (all recomputed properly below by the red_veg3_couple
        ! call, not needed here)
        veg_state%vegCpft(l,n),                                                &
        veg_state%leafC(l,n),veg_state%rootC(l,n),veg_state%woodC(l,n),        &
        veg_state%vegCpft_CAI(l,n),veg_state%leafC_CAI(l,n),                   &
        veg_state%rootC_CAI(l,n),veg_state%woodC_CAI(l,n),                     &
        veg_state%canht(l,n),veg_state%lai_bal(l,n)                            &
        )


      ! Here we estimate the phenology diagnosed from the lai and lai_bal
      ! from the dump/intialisation.
      IF (veg_state%lai_bal(l,n) > 0.0 .AND.                                   &
          veg_state%lai_bal(l,n) >= veg_state%lai(l,n)) THEN
        veg_state%phen(l,n) = veg_state%lai(l,n) / veg_state%lai_bal(l,n)

      ELSE
        veg_state%phen(l,n) = TINY(0.0)

      END IF

    END DO
  END DO

  ! Initialise veg state from red prognostic
  CALL red_veg3_couple(ainfo)

  !-----------------------------------------------------------------------------
  ! AJW <<END
  !-----------------------------------------------------------------------------

END IF

RETURN
END SUBROUTINE veg3_set_fields

!-------------------------------------------------------------------------------
SUBROUTINE veg3_red_set_fields(nnpft,nmasst)

!Source parms, etc from io modules - these should be available on and offline

USE veg3_canopy_mod, ONLY: veg3_canopy_height_order
USE veg3_parm_mod, ONLY: veg3_ctrl
USE jules_surface_mod, ONLY: cmass
USE pftparm, ONLY: lma
USE ereport_mod, ONLY: ereport

IMPLICIT NONE

INTEGER, INTENT(IN) :: nnpft, nmasst

INTEGER :: n,k
INTEGER :: errorstatus
CHARACTER(LEN=100) :: message_text
CHARACTER(LEN=*), PARAMETER :: RoutineName='VEG3_RED_SET_FIELDS'

!End of header

! Setup Allometry
DO k = 1,nmasst
  DO n = 1,nnpft
    IF (k == 1) THEN
      red_state%mass_cohort(n,k) = red_parms%mass0(n)
      red_state%ht_mass(n,k) = red_parms%height0(n)
      red_state%crwn_area_mass(n,k) = red_parms%crwn_area0(n)
      red_state%g_mass_scale(n,k) = 1.0
      red_state%lai_bal_mass(n,k) = red_parms%lai_bal0(n)
    ELSE IF (k > 1 .AND. k <= red_parms%mclass(n)) THEN
      red_state%mass_cohort(n,k) = red_state%mass_cohort(n,k-1)                &
        * red_parms%mclass_geom_mult(n)
      red_state%ht_mass(n,k) = red_parms%height0(n)                            &
        * (red_state%mass_cohort(n,k) / red_parms%mass0(n))** red_parms%phi_h(n)
      red_state%crwn_area_mass(n,k) = red_parms%crwn_area0(n)                  &
        * (red_state%mass_cohort(n,k) / red_parms%mass0(n))** red_parms%phi_a(n)
      red_state%g_mass_scale(n,k) = (red_state%mass_cohort(n,k)                &
        / red_parms%mass0(n))** red_parms%phi_g(n)
      red_state%lai_bal_mass(n,k) = red_parms%lai_bal0(n)                      &
        * (red_state%mass_cohort(n,k) / red_parms%mass0(n))** red_parms%phi_l(n)
    END IF

    IF (k <= red_parms%mclass(n)) THEN
      ! Individual leaf/fine-root/wood carbon mass within this mass class,
      ! partitioned out of the individual's total mass (mass_cohort). Fine
      ! root mass is assumed equal to leaf mass (as for the PFT-area-mean
      ! leafC/rootC previously calculated directly from lai_bal).
      red_state%leaf_mass(n,k) = cmass * lma(n) * red_state%lai_bal_mass(n,k)  &
        * red_state%crwn_area_mass(n,k)
      red_state%froot_mass(n,k) = red_state%leaf_mass(n,k)
      red_state%wood_mass(n,k) = red_state%mass_cohort(n,k)                    &
        - red_state%leaf_mass(n,k) - red_state%froot_mass(n,k)

      IF (red_state%wood_mass(n,k) < 0.0) THEN
        errorstatus = 101
        WRITE(message_text, '(A,I0,A,I0)') 'Negative wood_mass diagnosed ' //  &
          'for pft=', n, ', mass class=', k
        CALL ereport(RoutineName, errorstatus, message_text)
      END IF
    END IF
  END DO
END DO

! Height-order the (PFT, mass-class) pairs once (static allometry).
CALL veg3_canopy_height_order(nnpft, nmasst, veg3_ctrl%npft_totmclass,         &
  red_parms%mclass, red_state%ht_mass, red_state%order_pft_desc,               &
  red_state%order_mclass_desc, red_state%order_pft_asc,                        &
  red_state%order_mclass_asc)

RETURN
END SUBROUTINE veg3_red_set_fields

!-------------------------------------------------------------------------------

SUBROUTINE veg3_field_init(land_pts,nsurft,nnpft,npft,nmasst,ainfo,progs)

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, nnpft, npft, nsurft,nmasst

TYPE(ainfo_type), INTENT(IN OUT) :: ainfo
TYPE(progs_type), INTENT(IN) :: progs

!End of header

IF (l_red) THEN
  ! Set up all the necessary fields in RED
  CALL veg3_red_set_fields(nnpft,nmasst)
  ! Set up the rest of the fields for VEG3
  CALL veg3_set_fields(land_pts,nsurft,nnpft,nmasst,ainfo,progs)
END IF

RETURN
END SUBROUTINE veg3_field_init
!-----------------------------------------------------------------------------

SUBROUTINE red_veg3_couple(ainfo)
!-----------------------------------------------------------------------------
! Main coupling routine to aggregate between mass tiles to pft properties
! provides diagnostics and the main coupling fields back to JULES
! Same routine called during initialisation and runtime
!-----------------------------------------------------------------------------

USE ancil_info,                   ONLY: land_pts, nmasst, ainfo_type
USE jules_surface_types_mod,      ONLY: nnpft, soil
USE gridbox_mean_mod,             ONLY: pfttiles_to_gbm,                       &
                                        masstiles_to_pfttiles
USE conversions_mod, ONLY: rsec_per_day
USE veg3_parm_mod, ONLY: veg3_ctrl
USE veg3_canopy_mod, ONLY: veg3_canopy_frac
USE pftparm, ONLY: lma

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(INOUT).
!-----------------------------------------------------------------------------
TYPE(ainfo_type),INTENT(IN OUT) :: ainfo

!Local
INTEGER :: l,n  ! Index variables.

REAL :: frac_old(land_pts,nnpft)
        ! PFT fraction before this call updates veg_state%frac.
REAL :: vegCpft_old(land_pts,nnpft)
        ! PFT carbon density before this call updates veg_state%vegCpft.
        ! Used to derive the implicit litter flux below.
REAL :: frac_flux(land_pts,nnpft)
        ! Representative PFT fraction used to convert per-PFT-area fluxes
        ! to/from the gridbox mean over this coupling step, taken as the
        ! midpoint of the old and new PFT fraction.
REAL :: frac_above_mclass1(land_pts,nnpft)
        ! frac_above at the moment each PFT's lowest mass class was
        ! processed (see veg3_canopy_frac). Unused here.
REAL :: cai_frac_scale(land_pts,nnpft)
        ! Ratio of CAI to frac, used to convert the per-CAI NPP/litter
        ! fluxes output by photosynthesis/litter calculations to a per-PFT
        ! ground area (frac) basis for output.

!-----------------------------------------------------------------------------
!end of header

frac_old(:,:) = veg_state%frac(:,1:nnpft)
vegCpft_old(:,:) = veg_state%vegCpft(:,:)

veg_state%vegCpft(:,:)  = 0.0
veg_state%lai_bal(:,:)  = 0.0
veg_state%canht(:,:)    = 0.0
veg_state%lai_bal(:,:) = 0.0
veg_state%canht(:,:)   = 0.0
veg_state%vegC(:) = 0.0
veg_state%litC(:) = 0.0
veg_state%npp_n_gb(:) = 0.0
veg_state%npp_gb(:) = 0.0

! Estimate the vegetation fraction and CAI.
CALL veg3_canopy_frac(land_pts, nnpft, nmasst, veg3_ctrl%npft_totmclass,       &
  red_state%order_pft_desc, red_state%order_mclass_desc, red_parms%k_cai,      &
  red_parms%frac_excl, veg_state%frac_tile_excl, red_state%crwn_area_mass,     &
  red_state%plantNumDensity, red_state%frac_mass, red_state%CAI_mass,          &
  veg_state%frac(:,1:nnpft), veg_state%cai, frac_above_mclass1)

DO n = 1,nnpft
  DO l = 1,land_pts

    ! Estimate the PFT mean physical properties (per-frac and per-CAI, see
    ! pft_mean_from_mass_class) by aggregating across the mass class
    ! structure.
    CALL pft_mean_from_mass_class(                                             &
      !IN sizing
      red_parms%mclass(n),                                                     &
      !IN mass-cohort properties
      red_state%plantNumDensity(l,n,1:red_parms%mclass(n)),                    &
      red_state%mass_cohort(n,1:red_parms%mclass(n)),                          &
      red_state%lai_bal_mass(n,1:red_parms%mclass(n)),                         &
      red_state%ht_mass(n,1:red_parms%mclass(n)),                              &
      red_state%leaf_mass(n,1:red_parms%mclass(n)),                            &
      red_state%froot_mass(n,1:red_parms%mclass(n)),                           &
      red_state%wood_mass(n,1:red_parms%mclass(n)),                            &
      red_state%crwn_area_mass(n,1:red_parms%mclass(n)),                       &
      !IN PFT trait parameter used to invert leafC for lai_bal
      lma(n),                                                                  &
      !IN overlap-corrected PFT ground area fraction, by mass class
      red_state%frac_mass(l,n,1:red_parms%mclass(n)),                          &
      !IN total PFT fraction and crown area index
      veg_state%frac(l,n),veg_state%cai(l,n),                                  &
      !OUT fields
      veg_state%vegCpft(l,n),                                                  &
      veg_state%leafC(l,n),veg_state%rootC(l,n),veg_state%woodC(l,n),          &
      veg_state%vegCpft_CAI(l,n),veg_state%leafC_CAI(l,n),                     &
      veg_state%rootC_CAI(l,n),veg_state%woodC_CAI(l,n),                       &
      veg_state%canht(l,n),veg_state%lai_bal(l,n)                              &
      )

    veg_state%lai(l,n)   = veg_state%phen(l,n) * veg_state%lai_bal(l,n)

    ! NPP and its derived litter fluxes are normalised per unit PFT canopy
    ! area. frac_flux is the representative PFT fraction used to convert
    ! between per-PFT-area and gridbox-mean quantities over this coupling
    ! step, taken as the midpoint of the old and new PFT fraction (mirroring
    ! TRIFFID's own frac_flux normalisation of lit_c).
    frac_flux(l,n) = 0.5 * (frac_old(l,n) + veg_state%frac(l,n))

    ! Rescale NPP and litterfall fluxes by frac_old/frac_flux so the total
    ! gridbox flux is retained under the new PFT area. Left unchanged if
    ! frac_flux is zero.
    IF (frac_flux(l,n) > 0.0) THEN
      veg_state%npp_dr_out(l,n) = veg_state%npp_dr_out(l,n) * frac_old(l,n)    &
                                   / frac_flux(l,n)
      veg_state%leaf_litC(l,n)  = veg_state%leaf_litC(l,n)  * frac_old(l,n)    &
                                   / frac_flux(l,n)
      veg_state%root_litC(l,n)  = veg_state%root_litC(l,n)  * frac_old(l,n)    &
                                   / frac_flux(l,n)
      veg_state%wood_litC(l,n)  = veg_state%wood_litC(l,n)  * frac_old(l,n)    &
                                   / frac_flux(l,n)
      veg_state%mort_litC(l,n)  = veg_state%mort_litC(l,n)  * frac_old(l,n)    &
                                   / frac_flux(l,n)
    END IF

    ! npp_dr_out and the litter fluxes above are on a per-CAI basis (see
    ! sf_stom/cancap and veg3_litter_mod); convert to per-frac for output
    ! and for the litCpft_implicit call below.
    IF (veg_state%frac(l,n) > 0.0) THEN
      cai_frac_scale(l,n) = veg_state%cai(l,n) / veg_state%frac(l,n)
    ELSE
      cai_frac_scale(l,n) = 0.0
    END IF

    ! Cap npp_dr_out at CAI=1 (as for growth - see next_gen_biogeochem_mod)
    ! before converting from per-CAI to per-PFT ground area (frac) density.
    veg_state%npp_dr_out(l,n) =  veg_state%cai(l,n) *                          &
                                 veg_state%npp_dr_out(l,n)
    veg_state%npp_dr_out(l,n) = veg_state%npp_dr_out(l,n) * cai_frac_scale(l,n)
    veg_state%leaf_litC(l,n)  = veg_state%leaf_litC(l,n)  * cai_frac_scale(l,n)
    veg_state%root_litC(l,n)  = veg_state%root_litC(l,n)  * cai_frac_scale(l,n)
    veg_state%wood_litC(l,n)  = veg_state%wood_litC(l,n)  * cai_frac_scale(l,n)

    ! Derive litCpft implicitly from carbon conservation (NPP minus the
    ! change in standing PFT carbon), which guarantees the litter flux is
    ! consistent with the actual change in vegetation carbon.
    CALL litCpft_implicit(veg_state%npp_dr_out(l,n), frac_old(l,n),            &
                           veg_state%frac(l,n), frac_flux(l,n),                &
                           vegCpft_old(l,n), veg_state%vegCpft(l,n),           &
                           veg3_ctrl%dt_red, veg_state%litCpft(l,n))

    ! Update bare soil, also subtracts from urban, lake and snow fractions
    veg_state%frac(l,soil) = MAX(0.0, 1.0 - SUM(veg_state%frac(l,1:nnpft)) -   &
                             veg_state%frac_tile_excl(l))

    ! Aggregate the per-PFT contributions to update the gridbox-level fields
    veg_state%vegC(l) = veg_state%vegC(l) + veg_state%vegCpft(l,n) *           &
                        veg_state%frac(l,n)
    veg_state%npp_n_gb(l) = veg_state%npp_n_gb(l) + veg_state%npp_dr_out(l,n) *&
                            frac_flux(l,n)
    veg_state%npp_gb(l) = veg_state%npp_gb(l) + veg_state%npp_dr_out(l,n) *    &
                          frac_flux(l,n) / (rsec_per_day * 360.0)
    veg_state%litC(l) = veg_state%litC(l) + veg_state%litCpft(l,n) *           &
                        frac_flux(l,n)

  END DO

END DO

! Litter must be positive for the 4-pool soil carbon model (soilcarb_layers
! calculates the DPM:RPM ratio as a fraction of the total litter, which
! would otherwise divide by zero when a point has no litter). If there is
! no litter, make it a very small positive value.
DO l = 1,land_pts
  IF (SUM(veg_state%litCpft(l,:)) == 0.0) THEN
    veg_state%litCpft(l,1) = TINY(0.0)
  END IF
END DO


RETURN
END SUBROUTINE red_veg3_couple

!-----------------------------------------------------------------------------
SUBROUTINE pft_mean_from_mass_class(                                           &
                !IN sizing
                mclass,                                                        &
                !IN mass-cohort properties
                plantNumDensity,mass_cohort,lai_bal_mass,ht_mass,              &
                leaf_mass,froot_mass,wood_mass,crwn_area_mass,                 &
                !IN PFT trait parameter used to invert leafC for lai_bal
                lma_pft,                                                       &
                !IN overlap-corrected PFT ground area fraction, by mass class
                frac_mass,                                                     &
                !IN total PFT fraction and crown area index (CAI)
                frac,CAI_pft,                                                  &
                !OUT fields, per unit frac (whole-population)
                vegCpft,leafC,rootC,woodC,                                     &
                !OUT fields, per unit CAI (open-canopy)
                vegCpft_CAI,leafC_CAI,rootC_CAI,woodC_CAI,canht,lai_bal        &
                )
!-----------------------------------------------------------------------------
! Aggregates the plant number density across the mass class structure of a
! single PFT at a single point into the PFT mean physical properties.
!
! vegCpft/leafC/rootC/woodC weight by the whole population
! (plantNumDensity), normalised per unit frac. canht, lai_bal and the
! _CAI fields weight by frac_mass (open-canopy, overlap-corrected) and
! normalise per unit CAI_pft instead, for use by the land-surface/
! photosynthesis scheme.
!-----------------------------------------------------------------------------

USE jules_surface_mod, ONLY: cmass

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
plantNumDensity(mclass),                                                       &
              !  Population density within each mass cohort. (m-2)
mass_cohort(mclass),                                                           &
              !  Mass of an individual member of each mass cohort. (kgC)
lai_bal_mass(mclass),                                                          &
              !  Balanced leaf area index of an individual member of each mass
              !  cohort. (m2 m-2)
ht_mass(mclass),                                                               &
              !  Height of an individual member of each mass cohort. (m)
leaf_mass(mclass),                                                             &
              !  Individual leaf carbon mass of each mass cohort. (kg C)
froot_mass(mclass),                                                            &
              !  Individual fine root carbon mass of each mass cohort. (kg C)
wood_mass(mclass),                                                             &
              !  Individual wood carbon mass of each mass cohort. (kg C)
crwn_area_mass(mclass),                                                        &
              !  Crown area of an individual member of each mass cohort.
              !  (m2)
frac_mass(mclass),                                                             &
              !  PFT closed-canopy overlap-corrected ground area fraction
              !  contributed by each mass cohort. (-)
frac,                                                                          &
              !  Total PFT ground area fraction, SUM(frac_mass). (-)
CAI_pft,                                                                       &
              !  Total PFT crown area index, SUM(CAI_mass). (m2 m-2)
lma_pft
              !  PFT leaf mass per unit area (1/SLA), used to invert leafC for
              !  the effective balanced LAI. (kg leaf m-2)

!-----------------------------------------------------------------------------
! Reals with INTENT OUT
!-----------------------------------------------------------------------------
REAL, INTENT(OUT)    ::                                                        &
vegCpft,                                                                       &
              !  Total PFT carbon density per unit PFT ground area. (kg C m-2)
leafC,                                                                         &
              !  PFT leaf carbon density per unit PFT ground area. (kg C m-2)
rootC,                                                                         &
              !  PFT fine root carbon density per unit PFT ground area.
              !  (kg C m-2)
woodC,                                                                         &
              !  PFT wood carbon density per unit PFT ground area. (kg C m-2)
vegCpft_CAI,                                                                   &
              !  Total PFT carbon density per unit crown area index (CAI).
              !  (kg C m-2)
leafC_CAI,                                                                     &
              !  Open-canopy PFT leaf carbon density per unit crown area
              !  index (CAI). (kg C m-2)
rootC_CAI,                                                                     &
              !  Open-canopy PFT fine root carbon density per unit crown
              !  area index (CAI). (kg C m-2)
woodC_CAI,                                                                     &
              !  Open-canopy PFT wood carbon density per unit crown area
              !  index (CAI). (kg C m-2)
canht,                                                                         &
              !  Canopy height per unit crown area index (CAI). (m)
lai_bal
              !  Open-canopy (bounded) balanced LAI per unit crown area
              !  index (CAI). (m2 m-2)

!-----------------------------------------------------------------------------
!Local Vars
!-----------------------------------------------------------------------------
INTEGER              :: k
REAL                 :: vegCpft_raw, canht_raw, leafC_raw, rootC_raw,          &
                        woodC_raw,                                             &
                        leafC_CAI_raw, rootC_CAI_raw, woodC_CAI_raw
              !  Unnormalised sums.

!End of headers

! Initialise vars
vegCpft_raw    = 0.0
canht_raw      = 0.0
leafC_raw      = 0.0
rootC_raw      = 0.0
woodC_raw      = 0.0
leafC_CAI_raw  = 0.0
rootC_CAI_raw  = 0.0
woodC_CAI_raw  = 0.0

! leafC/rootC/woodC aggregate the same way as vegCpft (whole-population,
! per unit frac). canht and the open-canopy carbon densities weight by
! frac_mass instead, converting leaf/root/wood mass to a density via
! crown area.
DO k = 1, mclass
  vegCpft_raw   = vegCpft_raw   + plantNumDensity(k) * mass_cohort(k)
  canht_raw     = canht_raw     + frac_mass(k) * ht_mass(k)
  leafC_raw     = leafC_raw     + plantNumDensity(k) * leaf_mass(k)
  rootC_raw     = rootC_raw     + plantNumDensity(k) * froot_mass(k)
  woodC_raw     = woodC_raw     + plantNumDensity(k) * wood_mass(k)
  leafC_CAI_raw = leafC_CAI_raw + frac_mass(k) * leaf_mass(k)                  &
                                   / crwn_area_mass(k)
  rootC_CAI_raw = rootC_CAI_raw + frac_mass(k) * froot_mass(k)                 &
                                   / crwn_area_mass(k)
  woodC_CAI_raw = woodC_CAI_raw + frac_mass(k) * wood_mass(k)                  &
                                   / crwn_area_mass(k)
END DO

! Normalise the whole-population fields per unit frac.
IF (frac > 0.0) THEN
  vegCpft = vegCpft_raw / frac
  leafC   = leafC_raw   / frac
  rootC   = rootC_raw   / frac
  woodC   = woodC_raw   / frac
ELSE

  ! If frac is zero, set mean to lowest mass class value
  vegCpft = mass_cohort(1)
  leafC   = leaf_mass(1)
  rootC   = froot_mass(1)
  woodC   = wood_mass(1)
END IF

! Normalise the open-canopy fields per unit CAI_pft, inverting
! leafC_CAI = cmass * lma * lai_bal to get the effective balanced LAI.
IF (CAI_pft > 0.0) THEN
  vegCpft_CAI = vegCpft_raw / CAI_pft
  canht       = canht_raw   / CAI_pft
  leafC_CAI   = leafC_CAI_raw / CAI_pft
  rootC_CAI   = rootC_CAI_raw / CAI_pft
  woodC_CAI   = woodC_CAI_raw / CAI_pft
  ! Think about lai_bal here.
  lai_bal     = leafC_CAI / (cmass * lma_pft)
ELSE

  ! If CAI_pft is zero, set mean to lowest mass class value
  vegCpft_CAI = mass_cohort(1)
  canht       = ht_mass(1)
  leafC_CAI   = leaf_mass(1)  / crwn_area_mass(1)
  rootC_CAI   = froot_mass(1) / crwn_area_mass(1)
  woodC_CAI   = wood_mass(1)  / crwn_area_mass(1)
  lai_bal     = lai_bal_mass(1)
END IF

END SUBROUTINE pft_mean_from_mass_class
!-----------------------------------------------------------------------------

SUBROUTINE litCpft_implicit(npp_dr, frac_old, frac_new, frac_flux,             &
                             vegCpft_old, vegCpft_new, dt, litCpft)
!-----------------------------------------------------------------------------
! Derives the PFT litter flux implicitly from carbon conservation (NPP minus
! the change in standing PFT carbon), rather than summing the explicit
! leaf/root/wood/mortality litter fluxes:
!   litCpft = npp_dr - d(vegCpft*frac) / (frac_flux*dt) * rsec_per_day*360
! frac_flux (the midpoint of frac_old and frac_new) normalises the fluxes,
! mirroring TRIFFID's own frac_flux normalisation of lit_c.
!-----------------------------------------------------------------------------
USE conversions_mod, ONLY: rsec_per_day

IMPLICIT NONE

REAL, INTENT(IN)  :: npp_dr
              !  PFT NPP driving RED, normalised per unit PFT area.
              !  (kg C m-2 (360d)-1)
REAL, INTENT(IN)  :: frac_old, frac_new
              !  PFT fraction before/after this call. (-)
REAL, INTENT(IN)  :: frac_flux
              !  Representative PFT fraction used to normalise the litter
              !  flux, taken as the midpoint of frac_old and frac_new. (-)
REAL, INTENT(IN)  :: vegCpft_old, vegCpft_new
              !  PFT carbon density per PFT area before/after this call.
              !  (kg C m-2)
REAL, INTENT(IN)  :: dt
              !  Vegetation dynamics timestep over which growth was applied
              !  and vegCpft changed. (s)
REAL, INTENT(OUT) :: litCpft
              !  Implicit litter flux consistent with carbon conservation.
              !  (kg C m-2 (360d)-1)

!End of header

IF (frac_flux > 0.0) THEN
  litCpft = npp_dr - (vegCpft_new * frac_new - vegCpft_old * frac_old)         &
            / (frac_flux * dt) * rsec_per_day * 360.0
ELSE
  litCpft = 0.0
END IF

END SUBROUTINE litCpft_implicit
!-----------------------------------------------------------------------------

END MODULE veg3_field_mod
