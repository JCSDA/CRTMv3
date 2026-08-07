!
! ODPS_Predictor
!
! Module containing routines to compute the optical depth predictors for
! the Optical Depth in Pressure Space (ODPS) algorithm.
!
!
! CREATION HISTORY:
!       Written by:     Yong Han, 29-Aug-2006
!                       yong.han@noaa.gov
!
!       Modified by:    Tong Zhu, 18-Nov-2008
!                       tong.zhu@noaa.gov
!
!       Modified by:    James Rosinski, 08-Feb-2019
!                       Rosinski@ucar.edu
!
! (C) Copyright 2019 UCAR
!

MODULE ODPS_Predictor

  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module use
  USE Type_Kinds              , ONLY: fp
  USE CRTM_Parameters         , ONLY: MINIMUM_ABSORBER_AMOUNT, &
                                      RECIPROCAL_GRAVITY
  USE CRTM_Atmosphere_Define  , ONLY: CRTM_Atmosphere_type
  USE CRTM_GeometryInfo_Define, ONLY: CRTM_GeometryInfo_type, &
                                      CRTM_GeometryInfo_GetValue
  USE ODPS_Predictor_Define   , ONLY: ODPS_Predictor_type, &
                                      PAFV_type          , &
                                      PAFV_Associated    , &
                                      MAX_OPTRAN_ORDER
  USE ODPS_CoordinateMapping  , ONLY: Map_Input           , &
                                      Map_Input_TL        , &
                                      Map_Input_AD        , &
                                      Compute_Interp_Index
  USE ODPS_Define             , ONLY: ODPS_type
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Datatypes
  PUBLIC :: iVar_type
  ! Procedures
  PUBLIC :: ODPS_Assemble_Predictors
  PUBLIC :: ODPS_Assemble_Predictors_TL
  PUBLIC :: ODPS_Assemble_Predictors_AD
  PUBLIC :: ODPS_Compute_Predictor
  PUBLIC :: ODPS_Compute_Predictor_TL
  PUBLIC :: ODPS_Compute_Predictor_AD
  PUBLIC :: ODPS_Compute_Predictor_ODAS
  PUBLIC :: ODPS_Compute_Predictor_ODAS_TL
  PUBLIC :: ODPS_Compute_Predictor_ODAS_AD
  PUBLIC :: ODPS_Get_max_n_Predictors
  PUBLIC :: ODPS_Get_n_Components
  PUBLIC :: ODPS_Get_n_Absorbers
  PUBLIC :: ODPS_Get_Component_ID
  PUBLIC :: ODPS_Get_Absorber_ID
  PUBLIC :: ODPS_Get_Ozone_Component_ID
  PUBLIC :: ODPS_Get_SaveFWVFlag
  PUBLIC :: ODPS_Validate_Group
  PUBLIC :: ODPS_Kernel_n_Predictors
  PUBLIC :: ODPS_Max_n_Predictors_For
  ! Parameters
  PUBLIC :: TOT_ComID
  PUBLIC :: WLO_ComID
  PUBLIC :: WET_ComID
  PUBLIC :: CO2_ComID
  PUBLIC :: GROUP_1
  PUBLIC :: GROUP_2
  PUBLIC :: GROUP_3
  PUBLIC :: GROUP_MW_O3
  PUBLIC :: GROUP_UV_NO2
  PUBLIC :: GROUP_IR_QDRY2
  PUBLIC :: GROUP_IR_QDRY1
  PUBLIC :: RESERVED_ZSSMIS_GROUP
  PUBLIC :: RESERVED_ZAMSUA_GROUP
  PUBLIC :: ODPS_INVALID_ID
  PUBLIC :: ALLOW_OPTRAN


  ! -----------------
  ! Module parameters
  ! -----------------
  ! The ODPS groups. The complete definition of every group (its predictor
  ! basis, component roster, absorber roster, and per-component predictor
  ! counts) lives in the single GROUP_REGISTRY table declared below, after
  ! the component and absorber ID constants it references.
  ! Group 7 is the MW+ozone variant of group 3 (indexes 4 - 6 are Zeeman);
  ! group 8 is the UV/VIS variant of group 2 with an added scene-NO2
  ! component; groups 9 and 10 are the humidity-aware-dry (QDRY) variants
  ! of groups 2 and 1, whose dry component carries the 13 water-vapor
  ! predictor formulas in addition to the 7 heritage dry predictors.
  INTEGER, PARAMETER  :: N_G = 10
  INTEGER, PARAMETER  :: MAX_COMPONENTS_ANY_GROUP = 8
  INTEGER, PARAMETER  :: MAX_ABSORBERS_ANY_GROUP  = 6
  ! Predictor basis classes: which shared per-layer formulation a group uses
  INTEGER, PARAMETER :: BASIS_RESERVED = 0  ! Zeeman-reserved; never dispatched here
  INTEGER, PARAMETER :: BASIS_IR       = 1  ! IR/VIS/UV formulation (groups 1, 2, 8)
  INTEGER, PARAMETER :: BASIS_MW       = 2  ! MW formulation (groups 3, 7)
  ! Group index. Group indexes 4 - 6 are RESERVED for the Zeeman sub-algorithms
  ! and are never valid in a standard ODPS TauCoeff file: 4 is Zeeman SSMIS,
  ! 5 is Zeeman AMSU-A, 6 is an unassigned Zeeman reserve. This module is the
  ! single owner of the group index space; the ODZeeman code derives its
  ! ODPS_gINDEX_* constants from the RESERVED_* parameters below. The zero
  ! entries at positions 4 - 6 in the dimension tables above are placeholders
  ! for these reserved indexes (Zeeman has its own predictor module and never
  ! reaches these tables).
  INTEGER, PARAMETER :: GROUP_1 = 1
  INTEGER, PARAMETER :: GROUP_2 = 2
  INTEGER, PARAMETER :: GROUP_3 = 3
  INTEGER, PARAMETER :: RESERVED_ZSSMIS_GROUP = 4  ! Zeeman SSMIS (ODZeeman only)
  INTEGER, PARAMETER :: RESERVED_ZAMSUA_GROUP = 5  ! Zeeman AMSU-A (ODZeeman only)
  INTEGER, PARAMETER :: GROUP_MW_O3 = 7   ! MW with a scene-ozone component
  INTEGER, PARAMETER :: GROUP_UV_NO2 = 8  ! UV/VIS with a scene-NO2 component
  INTEGER, PARAMETER :: GROUP_IR_QDRY2 = 9   ! group 2 with the QDRY dry component
  INTEGER, PARAMETER :: GROUP_IR_QDRY1 = 10  ! group 1 with the QDRY dry component
  ! The groups a standard ODPS TauCoeff file may legitimately carry
  INTEGER, PARAMETER :: VALID_GROUPS(7) = &
    (/ GROUP_1, GROUP_2, GROUP_3, GROUP_MW_O3, GROUP_UV_NO2, &
       GROUP_IR_QDRY2, GROUP_IR_QDRY1 /)

  ! Component IDs.
  !
  ! REGISTRY NOTE: component IDs are inherited from the heritage transmittance
  ! production "molecule set" numbering (see CRTM_coef,
  ! src/apps/TauProd/Infrared/Check_ProcessControl_File/Tau_Production_Parameters.f90):
  !   1 - 7   individual molecules (HITRAN order; 7 is O2)
  !   8/9/10  all-no-continua / continua-only / all-with-continua
  !   12      wet (water vapor, line and continua)
  !   13      dry
  !   14      ozone
  !   15      wco (water vapor continua only)
  !   20      dry, group-2 formulation
  !   101     effective molecule 1 (water vapor line only, "wlo")
  !   112-121 effective single-gas components (112 wet, 113 dry, 114 ozone,
  !           118 CH4, 119 CO, 120 N2O, 121 CO2)
  !   122     NO2 (CRTM extension, added with GROUP_UV_NO2)
  !   123/124 QDRY dry components (CRTM extension, added with the QDRY
  !           groups 9/10): the group-2/group-1 effective-dry TARGETS
  !           unchanged, but fitted and applied with 20 predictors (the 7
  !           heritage dry predictors plus water-vapor predictors 1 - 13
  !           of the WLO formulation). The band-averaged effective-dry
  !           transmittance carries the spectral-correlation term of the
  !           in-band water vs dry-gas line overlap, which responds to
  !           humidity; a temperature-only dry basis cannot represent
  !           that response (P3 investigation, 2026-08).
  ! New component IDs must be coordinated with the coefficient generation
  ! package (crtm-coeffgen) and recorded here.
  INTEGER,  PARAMETER :: TOT_ComID = 10    ! total tau
  INTEGER,  PARAMETER :: DRY_ComID_G1 =  7   ! dry gas for Group-1 sensors
  INTEGER,  PARAMETER :: DRY_ComID_G2 = 20   ! dry gas, for Gorup-2 sensors
  INTEGER,  PARAMETER :: DRY_ComID_Q2 = 123  ! QDRY dry, group-9 (group-2 target)
  INTEGER,  PARAMETER :: DRY_ComID_Q1 = 124  ! QDRY dry, group-10 (group-1 target)
  INTEGER,  PARAMETER :: WLO_ComID = 101  ! water vapor line only, no continua
  INTEGER,  PARAMETER :: WCO_ComID = 15   ! water vapor continua only, no line absorption
  INTEGER,  PARAMETER :: OZO_ComID = 114  ! ozone
  INTEGER,  PARAMETER :: CO2_ComID = 121  ! CO2
  INTEGER,  PARAMETER :: N2O_ComID = 120  ! N2O
  INTEGER,  PARAMETER :: CO_ComID  = 119  ! CO
  INTEGER,  PARAMETER :: CH4_ComID = 118  ! CH4
  INTEGER,  PARAMETER :: NO2_ComID = 122  ! NO2 (scene component, UV/VIS group 8)
  ! Sentinel returned by the ID accessors for an out-of-range query
  INTEGER,  PARAMETER :: ODPS_INVALID_ID = -1

  ! Microwave sensors
  INTEGER,  PARAMETER :: EDRY_ComID = 113  ! Effective dry
  INTEGER,  PARAMETER :: WET_ComID = 12  ! water vapor line & no continua

  ! IR sensor Component indexes (sequence index in an array)
  INTEGER, PARAMETER :: COMP_DRY_IR = 1
  INTEGER, PARAMETER :: COMP_WLO_IR = 2
  INTEGER, PARAMETER :: COMP_WCO_IR = 3
  INTEGER, PARAMETER :: COMP_OZO_IR = 4
  INTEGER, PARAMETER :: COMP_CO2_IR = 5
  INTEGER, PARAMETER :: COMP_N2O_IR = 6
  INTEGER, PARAMETER :: COMP_CO_IR  = 7
  INTEGER, PARAMETER :: COMP_CH4_IR = 8

  ! MW sensor Component indexes
  INTEGER, PARAMETER :: COMP_DRY_MW = 1
  INTEGER, PARAMETER :: COMP_WET_MW = 2
  INTEGER, PARAMETER :: COMP_OZO_MW = 3   ! GROUP_MW_O3 only

  ! UV/VIS group-8 NO2 component index (6th component of the group-2-based set)
  INTEGER, PARAMETER :: COMP_NO2_G8 = 6   ! GROUP_UV_NO2 only

  ! Absorber IDs (HITRAN)
  INTEGER, PARAMETER ::   H2O_ID =  1
  INTEGER, PARAMETER ::   CO2_ID =  2
  INTEGER, PARAMETER ::    O3_ID =  3
  INTEGER, PARAMETER ::   N2O_ID =  4
  INTEGER, PARAMETER ::    CO_ID =  5
  INTEGER, PARAMETER ::   CH4_ID =  6
  INTEGER, PARAMETER ::   NO2_ID = 10
  ! All gases CRTM's ODPS kernels know how to consume
  INTEGER, PARAMETER :: KNOWN_GAS_IDS(7) = &
    (/ H2O_ID, CO2_ID, O3_ID, N2O_ID, CO_ID, CH4_ID, NO2_ID /)

  ! Absorber (Molecule) indexes for accessing absorber profile array
  INTEGER,  PARAMETER :: ABS_H2O_IR = 1
  INTEGER,  PARAMETER :: ABS_O3_IR  = 2
  INTEGER,  PARAMETER :: ABS_CO2_IR = 3
  INTEGER,  PARAMETER :: ABS_N2O_IR = 4
  INTEGER,  PARAMETER :: ABS_CO_IR  = 5
  INTEGER,  PARAMETER :: ABS_CH4_IR = 6

  INTEGER,  PARAMETER :: ABS_H2O_MW = 1
  INTEGER,  PARAMETER :: ABS_O3_MW  = 2   ! GROUP_MW_O3 only

  ! UV/VIS group-8 absorber array is [H2O,O3,CO2,NO2]; the first three reuse the
  ! IR indexes (ABS_H2O_IR/ABS_O3_IR/ABS_CO2_IR = 1/2/3), NO2 is the 4th.
  INTEGER,  PARAMETER :: ABS_NO2_G8 = 4   ! GROUP_UV_NO2 only

  ! ---------------------------------------------------------------------
  ! THE GROUP REGISTRY: one row per ODPS group, the complete definition in
  ! one place (basis class, component roster, absorber roster, and the
  ! per-component predictor counts). Rows 4 to 6 are the Zeeman-reserved
  ! placeholders (BASIS_RESERVED, all-zero rosters). The rosters are padded
  ! with zeros to the fixed component/absorber maximums; only the first
  ! n_Components / n_Absorbers entries are meaningful.
  !
  ! To add a group: add one row here, extend N_G, add the group's named
  ! index constant above, add it to VALID_GROUPS, and provide predictor
  ! kernels for any component ID not already handled (see the kernel
  ! dispatch in ODPS_Compute_Predictor and its TL/AD companions).
  ! ---------------------------------------------------------------------
  TYPE :: ODPS_Group_Spec_type
    CHARACTER(12) :: Name
    INTEGER       :: Basis
    INTEGER       :: n_Components
    INTEGER       :: n_Absorbers
    INTEGER       :: Max_n_Predictors
    INTEGER       :: Component_ID(MAX_COMPONENTS_ANY_GROUP)
    INTEGER       :: Absorber_ID(MAX_ABSORBERS_ANY_GROUP)
    INTEGER       :: n_Predictors(MAX_COMPONENTS_ANY_GROUP)
  END TYPE ODPS_Group_Spec_type

  TYPE(ODPS_Group_Spec_type), PARAMETER :: GROUP_REGISTRY(N_G) = (/ &
    ODPS_Group_Spec_type( 'IR_HIRES    ', BASIS_IR, 8, 6, 18, &
      (/ DRY_ComID_G1, WLO_ComID, WCO_ComID, OZO_ComID, CO2_ComID, N2O_ComID, CO_ComID, CH4_ComID /), &
      (/ H2O_ID, O3_ID, CO2_ID, N2O_ID, CO_ID, CH4_ID /), &
      (/ 7, 18, 7, 11, 11, 14, 10, 11 /) ), &
    ODPS_Group_Spec_type( 'IR_BROAD    ', BASIS_IR, 5, 3, 15, &
      (/ DRY_ComID_G2, WLO_ComID, WCO_ComID, OZO_ComID, CO2_ComID, 0, 0, 0 /), &
      (/ H2O_ID, O3_ID, CO2_ID, 0, 0, 0 /), &
      (/ 7, 15, 7, 11, 10, 0, 0, 0 /) ), &
    ODPS_Group_Spec_type( 'MW          ', BASIS_MW, 2, 1, 14, &
      (/ EDRY_ComID, WET_ComID, 0, 0, 0, 0, 0, 0 /), &
      (/ H2O_ID, 0, 0, 0, 0, 0 /), &
      (/ 7, 14, 0, 0, 0, 0, 0, 0 /) ), &
    ODPS_Group_Spec_type( 'RSVD_ZSSMIS ', BASIS_RESERVED, 0, 0, 0, &
      (/ 0, 0, 0, 0, 0, 0, 0, 0 /), (/ 0, 0, 0, 0, 0, 0 /), &
      (/ 0, 0, 0, 0, 0, 0, 0, 0 /) ), &
    ODPS_Group_Spec_type( 'RSVD_ZAMSUA ', BASIS_RESERVED, 0, 0, 0, &
      (/ 0, 0, 0, 0, 0, 0, 0, 0 /), (/ 0, 0, 0, 0, 0, 0 /), &
      (/ 0, 0, 0, 0, 0, 0, 0, 0 /) ), &
    ODPS_Group_Spec_type( 'RSVD_ZEEMAN3', BASIS_RESERVED, 0, 0, 0, &
      (/ 0, 0, 0, 0, 0, 0, 0, 0 /), (/ 0, 0, 0, 0, 0, 0 /), &
      (/ 0, 0, 0, 0, 0, 0, 0, 0 /) ), &
    ODPS_Group_Spec_type( 'MW_O3       ', BASIS_MW, 3, 2, 14, &
      (/ EDRY_ComID, WET_ComID, OZO_ComID, 0, 0, 0, 0, 0 /), &
      (/ H2O_ID, O3_ID, 0, 0, 0, 0 /), &
      (/ 7, 14, 11, 0, 0, 0, 0, 0 /) ), &
    ODPS_Group_Spec_type( 'UV_NO2      ', BASIS_IR, 6, 4, 15, &
      (/ DRY_ComID_G2, WLO_ComID, WCO_ComID, OZO_ComID, CO2_ComID, NO2_ComID, 0, 0 /), &
      (/ H2O_ID, O3_ID, CO2_ID, NO2_ID, 0, 0 /), &
      (/ 7, 15, 7, 11, 10, 3, 0, 0 /) ), &
    ODPS_Group_Spec_type( 'IR_BROAD_QD ', BASIS_IR, 5, 3, 20, &
      (/ DRY_ComID_Q2, WLO_ComID, WCO_ComID, OZO_ComID, CO2_ComID, 0, 0, 0 /), &
      (/ H2O_ID, O3_ID, CO2_ID, 0, 0, 0 /), &
      (/ 20, 15, 7, 11, 10, 0, 0, 0 /) ), &
    ODPS_Group_Spec_type( 'IR_HIRES_QD ', BASIS_IR, 8, 6, 20, &
      (/ DRY_ComID_Q1, WLO_ComID, WCO_ComID, OZO_ComID, CO2_ComID, N2O_ComID, CO_ComID, CH4_ComID /), &
      (/ H2O_ID, O3_ID, CO2_ID, N2O_ID, CO_ID, CH4_ID /), &
      (/ 20, 18, 7, 11, 11, 14, 10, 11 /) ) /)

  ! Literal constants
  REAL(fp), PARAMETER :: ZERO      = 0.0_fp
  REAL(fp), PARAMETER :: ONE       = 1.0_fp
  REAL(fp), PARAMETER :: TWO       = 2.0_fp
  REAL(fp), PARAMETER :: THREE     = 3.0_fp
  REAL(fp), PARAMETER :: FOUR      = 4.0_fp
  REAL(fp), PARAMETER :: TEN       = 10.0_fp
  REAL(fp), PARAMETER :: POINT_25  = 0.25_fp
  REAL(fp), PARAMETER :: POINT_5   = 0.5_fp
  REAL(fp), PARAMETER :: POINT_75  = 0.75_fp
  REAL(fp), PARAMETER :: ONE_POINT_5   = 1.5_fp
  REAL(fp), PARAMETER :: ONE_POINT_25  = 1.25_fp
  REAL(fp), PARAMETER :: ONE_POINT_75  = 1.75_fp


  LOGICAL, PARAMETER  :: ALLOW_OPTRAN = .TRUE.


  ! ------------------------------------------
  ! Structure definition to hold forward model
  ! variables across FWD, TL, and AD calls
  ! ------------------------------------------
  TYPE :: iVar_type
    PRIVATE
    INTEGER :: dummy
  END TYPE iVar_type

CONTAINS


!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

!--------------------------------------------------------------------------------
!
! NAME:
!       ODPS_Assemble_Predictors
!
! PURPOSE:
!       Subroutine to assemble all the gas absorption model predictors
!       for the ODPS algorithm.
!
! CALLING SEQUENCE:
!       CALL ODPS_Assemble_Predictors( &
!              TC       , &  ! Input
!              Atm      , &  ! Input
!              GeoInfo  , &  ! Input
!              Predictor  )  ! Output
!
! INPUT ARGUMENTS:
!       TC:           ODPS structure holding tau coefficients
!                        UNITS:      N/A
!                        TYPE:       ODPS_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Atm       :     CRTM Atmosphere structure containing the atmospheric
!                       state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       GeoInfo     :   CRTM_GeometryInfo structure containing the
!                       view geometry information.
!                       UNITS:      N/A
!                       TYPE:       CRTM_GeometryInfo_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       Predictor:      Predictor structure containing the integrated absorber
!                       and predictor profiles.
!                       UNITS:      N/A
!                       TYPE:       ODPS_Predictor_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
!--------------------------------------------------------------------------------

  SUBROUTINE ODPS_Assemble_Predictors( &
    TC       , &
    Atm      , &
    GeoInfo  , &
    Predictor  )
    ! Arguments
    TYPE(ODPS_type)             , INTENT(IN)     :: TC
    TYPE(CRTM_Atmosphere_type)  , INTENT(IN)     :: Atm
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)     :: GeoInfo
    TYPE(ODPS_Predictor_type)   , INTENT(IN OUT) :: Predictor
    ! Local variables
    REAL(fp) :: Temperature(Predictor%n_Layers)
    REAL(fp) :: Absorber(Predictor%n_Layers, TC%n_Absorbers)
    INTEGER  :: H2O_idx
    REAL(fp) :: Secant_Sensor_Zenith


    ! Map data from user to internal fixed pressure layers/levels
    CALL Map_Input( &
      Atm                            , &
      TC                             , &
      GeoInfo                        , &
      Temperature                    , &
      Absorber                       , &
      Predictor%User_Level_LnPressure, &
      Predictor%Ref_Level_LnPressure , &
      Predictor%Secant_Zenith        , &
      H2O_idx                        , &
      Predictor%PAFV)


    ! ...Store the surface secant zenith angle
    CALL CRTM_GeometryInfo_GetValue( GeoInfo, Secant_Trans_Zenith = Secant_Sensor_Zenith )
    Predictor%Secant_Zenith_Surface = Secant_Sensor_Zenith


    ! Compute predictor
    CALL ODPS_Compute_Predictor( &
      TC%Group_index         , &
      TC%Component_ID        , &
      TC%Absorber_ID         , &
      Temperature            , &
      Absorber               , &
      TC%Ref_Level_Pressure  , &
      TC%Ref_Temperature     , &
      TC%Ref_Absorber        , &
      Predictor%Secant_Zenith, &
      Predictor                )
    ! ...Optional ODAS for water vapour lines
    IF( ALLOW_OPTRAN .AND. TC%n_OCoeffs > 0 )THEN
      CALL ODPS_Compute_Predictor_ODAS( &
        Temperature            , &
        Absorber(:,H2O_idx)    , &
        TC%Ref_Level_Pressure  , &
        TC%Ref_Pressure        , &
        Predictor%Secant_Zenith, &
        TC%Alpha               , &
        TC%Alpha_C1            , &
        TC%Alpha_C2            , &
        Predictor                )
    END IF
    ! ...Save the interpolation indices
    IF ( PAFV_Associated(Predictor%PAFV) ) THEN
      CALL Compute_Interp_Index( &
        Predictor%Ref_Level_LnPressure ,  &
        Predictor%User_Level_LnPressure,  &
        Predictor%PAFV%ODPS2User_Idx)
    END IF

  END SUBROUTINE ODPS_Assemble_Predictors

!--------------------------------------------------------------------------------
!
! NAME:
!       ODPS_Assemble_Predictors_TL
!
! PURPOSE:
!       Subroutine to assemble all the tangent-linear gas absorption model
!       predictors.
!       It first interpolates the user temperature and absorber profiles on the
!       internal pressure grids and then calls the predictor computation routine
!       to compute the predictors
!
! CALLING SEQUENCE:
!       CALL ODPS_Assemble_Predictors_TL( &
!         TC          , &
!         Predictor   , &
!         Atm_TL      , &
!         Predictor_TL  )
!
! INPUT ARGUMENTS:
!          TC:           ODPS structure holding tau coefficients
!                        UNITS:      N/A
!                        TYPE:       ODPS_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Atm_TL    :     CRTM Atmosphere structure containing the atmospheric
!                       state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Predictor:      Predictor structure containing the integrated absorber
!                       and predictor profiles.
!                       UNITS:      N/A
!                       TYPE:       ODPS_Predictor_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       Predictor_TL:   Predictor structure containing the integrated absorber
!                       and predictor profiles.
!                       UNITS:      N/A
!                       TYPE:       ODPS_Predictor_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
!--------------------------------------------------------------------------------

  SUBROUTINE ODPS_Assemble_Predictors_TL( &
    TC          , &  ! Input
    Predictor   , &  ! Input
    Atm_TL      , &  ! Input
    Predictor_TL  )  ! Output
    ! Arguments
    TYPE(ODPS_type)           , INTENT(IN)     :: TC
    TYPE(ODPS_Predictor_type) , INTENT(IN)     :: Predictor
    TYPE(CRTM_Atmosphere_type), INTENT(IN)     :: Atm_TL
    TYPE(ODPS_Predictor_type) , INTENT(IN OUT) :: Predictor_TL
    ! Local variables
    REAL(fp) :: Absorber_TL(Predictor%n_Layers, TC%n_Absorbers)
    REAL(fp) :: Temperature_TL(Predictor%n_Layers)


    ! Map data from user to internal fixed pressure layers/levels
    CALL Map_Input_TL( &
      TC            , &
      Atm_TL        , &
      Temperature_TL, &
      Absorber_TL   , &
      Predictor%PAFV  )


    ! Compute predictor
    CALL ODPS_Compute_Predictor_TL( &
      TC%Group_index            , &
      TC%Component_ID           , &
      TC%Absorber_ID            , &
      Predictor%PAFV%Temperature, &
      Predictor%PAFV%Absorber   , &
      TC%Ref_Temperature        , &
      TC%Ref_Absorber           , &
      Predictor%Secant_Zenith   , &
      Predictor                 , &
      Temperature_TL            , &
      Absorber_TL               , &
      Predictor_TL                )
    ! ...Optional ODAS for water vapour lines
    IF ( ALLOW_OPTRAN .AND. TC%n_OCoeffs > 0 ) THEN
      CALL ODPS_Compute_Predictor_ODAS_TL( &
        Predictor%PAFV%Temperature                       , &
        Predictor%PAFV%Absorber(:,Predictor%PAFV%H2O_idx), &
        TC%Ref_Pressure                                  , &
        Predictor%Secant_Zenith                          , &
        TC%Alpha                                         , &
        TC%Alpha_C2                                      , &
        Predictor                                        , &
        Temperature_TL                                   , &
        Absorber_TL(:,Predictor%PAFV%H2O_idx)            , &
        Predictor_TL                                       )
    END IF

  END SUBROUTINE ODPS_Assemble_Predictors_TL


!--------------------------------------------------------------------------------
!
! NAME:
!       ODPS_Assemble_Predictors_AD
!
! PURPOSE:
!       Subroutine to assemble the adjoint of the gas absorption model
!       predictors.
!       It first calss the adjoint of the predictor computation routine and
!       then performs the adjoint interpolation of the user temperature and
!       absorber profiles on the internal pressure grid.
!
! CALLING SEQUENCE:
!       CALL ODPS_Assemble_Predictors_AD( &
!         TC          , &
!         Predictor   , &
!         Predictor_AD, &
!         Atm_AD        )
!
! INPUT ARGUMENTS:
!          TC:           ODPS structure holding tau coefficients
!                        UNITS:      N/A
!                        TYPE:       ODPS_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Predictor:      Predictor structure containing the integrated absorber
!                       and predictor profiles.
!                       UNITS:      N/A
!                       TYPE:       ODPS_Predictor_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Predictor_AD:   Predictor structure containing the integrated absorber
!                       and predictor profiles.
!                       UNITS:      N/A
!                       TYPE:       ODPS_Predictor_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!
!       Atm_AD    :     CRTM Atmosphere structure containing the atmospheric
!                       state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
!--------------------------------------------------------------------------------

  SUBROUTINE ODPS_Assemble_Predictors_AD( &
    TC          , &
    Predictor   , &
    Predictor_AD, &
    Atm_AD        )
    ! Arguments
    TYPE(ODPS_type)           , INTENT(IN)     :: TC
    TYPE(ODPS_Predictor_type) , INTENT(IN)     :: Predictor
    TYPE(ODPS_Predictor_type) , INTENT(IN OUT) :: predictor_AD
    TYPE(CRTM_Atmosphere_type), INTENT(IN OUT) :: Atm_AD
    ! Local variables
    REAL(fp) :: Absorber_AD(Predictor%n_Layers, TC%n_Absorbers)
    REAL(fp) :: Temperature_AD(Predictor%n_Layers)

    ! Initialise local adjoint variables
    Temperature_AD = ZERO
    Absorber_AD    = ZERO


    ! Compute predictor
    ! ...Optional ODAS for water vapour lines
    IF ( ALLOW_OPTRAN .AND. TC%n_OCoeffs > 0 ) THEN
      CALL ODPS_Compute_Predictor_ODAS_AD( &
        Predictor%PAFV%Temperature                       , &
        Predictor%PAFV%Absorber(:,Predictor%PAFV%H2O_idx), &
        TC%Ref_Pressure                                  , &
        Predictor%Secant_Zenith                          , &
        TC%Alpha                                         , &
        TC%Alpha_C2                                      , &
        Predictor                                        , &
        Predictor_AD                                     , &
        Temperature_AD                                   , &
        Absorber_AD(:,Predictor%PAFV%H2O_idx)              )
    END IF
    ! ...The main ODPS predictor
    CALL ODPS_Compute_Predictor_AD( &
      TC%Group_index            , &
      TC%Component_ID           , &
      TC%Absorber_ID            , &
      Predictor%PAFV%Temperature, &
      Predictor%PAFV%Absorber   , &
      TC%Ref_Temperature        , &
      TC%Ref_Absorber           , &
      Predictor%Secant_Zenith   , &
      Predictor                 , &
      Predictor_AD              , &
      Temperature_AD            , &
      Absorber_AD                 )


    ! Map data from user to internal fixed pressure layers/levels
    CALL Map_Input_AD( &
      TC            , &
      Temperature_AD, &
      Absorber_AD   , &
      Atm_AD        , &
      Predictor%PAFV  )

  END SUBROUTINE ODPS_Assemble_Predictors_AD


!------------------------------------------------------------------------------
!
! NAME:
!       ODPS_Compute_Predictor
!
! PURPOSE:
!       Subroutine to predictors
!
! CALLING SEQUENCE:
!       CALL ODPS_Compute_Predictor( &
!         Group_ID,           &  ! Input
!         Temperature,        &  ! Input
!         Absorber,           &  ! Input
!         Ref_Level_Pressure, &  ! Input
!         Ref_Temperature,    &  ! Input
!         Ref_Absorber,       &  ! Input
!         secang,             &  ! Input
!         Predictor )            ! Output
!
! INPUT ARGUMENTS:
!       Group_ID   :     The ID of predictor group
!                        UNITS:      N?A
!                        TYPE:       INTEGER
!                        DIMENSION:  scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Temperature:     Temperature profile
!                        UNITS:      K
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!       Absorber   :     Absorber profiles
!                        UNITS:      vary
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-2(n_Layers x n_Absorbers) array
!                        ATTRIBUTES: INTENT(IN)
!
!       Ref_Level_Pressure : Reference level pressure profile
!                        UNITS:      hPa
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(0:n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!       Ref_Temperature : Reference layer temperature profile
!                        UNITS:      K
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!       Ref_Absorber :   Reference absorber profiles
!                        UNITS:      vary
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-2(n_Layers x n_Absorbers) array
!                        ATTRIBUTES: INTENT(IN)
!
!       secang       :   Secont sensor zenith angle profile
!                        UNITS:      N/A
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       Predictor:      Predictor structure containing the integrated absorber
!                       and predictor profiles.
!                       UNITS:      N/A
!                       TYPE:       ODPS_Predictor_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
!------------------------------------------------------------------------------

  SUBROUTINE ODPS_Compute_Predictor( &
    Group_ID,           &
    Component_ID,       &
    Absorber_ID,        &
    Temperature,        &
    Absorber,           &
    Ref_Level_Pressure, &
    Ref_Temperature,    &
    Ref_Absorber,       &
    secang,             &
    Predictor )

    INTEGER,                   INTENT(IN)     :: Group_ID
    INTEGER,                   INTENT(IN)     :: Component_ID(:)
    INTEGER,                   INTENT(IN)     :: Absorber_ID(:)
    REAL(fp),                  INTENT(IN)     :: Temperature(:)
    REAL(fp),                  INTENT(IN)     :: Absorber(:, :)
    REAL(fp),                  INTENT(IN)     :: Ref_Level_Pressure(0:)
    REAL(fp),                  INTENT(IN)     :: Ref_Temperature(:)
    REAL(fp),                  INTENT(IN)     :: Ref_Absorber(:, :)
    REAL(fp),                  INTENT(IN)     :: Secang(:)
    TYPE(ODPS_Predictor_type), INTENT(IN OUT) :: Predictor

    ! ---------------
    ! Local variables
    ! ---------------
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'ODPS_Compute_Predictor'
    INTEGER    ::    n_layers
    INTEGER    ::    k ! n_Layers, n_Levels
    INTEGER    ::    j ! n_absorbers
    REAL(fp) ::    PDP
    REAL(fp) ::    Tzp_ref
    REAL(fp) ::    Tzp_sum
    REAL(fp) ::    Tzp(SIZE(Absorber, DIM=1))
    REAL(fp) ::    Tz_ref
    REAL(fp) ::    Tz_sum
    REAL(fp) ::    Tz(SIZE(Absorber, DIM=1))
    REAL(fp) ::    GAz_ref(SIZE(Absorber, DIM=2))
    REAL(fp) ::    GAz_sum(SIZE(Absorber, DIM=2))
    REAL(fp) ::    GAz(SIZE(Absorber, DIM=1), SIZE(Absorber, DIM=2))
    REAL(fp) ::    GAzp_ref(SIZE(Absorber, DIM=2))
    REAL(fp) ::    GAzp_sum(SIZE(Absorber, DIM=2))
    REAL(fp) ::    GAzp(SIZE(Absorber, DIM=1), SIZE(Absorber, DIM=2))
    REAL(fp) ::    GATzp_ref(SIZE(Absorber, DIM=2))
    REAL(fp) ::    GATzp_sum (SIZE(Absorber, DIM=2))
    REAL(fp) ::    GATzp(SIZE(Absorber, DIM=1), SIZE(Absorber, DIM=2))
    ! Kernel dispatch bookkeeping and shared per-layer variables
    INTEGER  :: ic, np
    INTEGER  :: ja_h2o, ja_o3, ja_co2, ja_n2o, ja_co, ja_ch4, ja_no2
    LOGICAL  :: has_trace
    REAL(fp) :: DT, T, T2, DT2
    REAL(fp) :: H2O, H2O_A, H2O_R, H2O_S, H2O_R4, H2OdH2OTzp
    REAL(fp) :: CO2, O3, O3_A, O3_R
    REAL(fp) :: NO2, NO2_A
    REAL(fp) :: CO, CO_A, CO_R, CO_S, CO_ACOdCOzp
    REAL(fp) :: N2O, N2O_A, N2O_R, N2O_S
    REAL(fp) :: CH4, CH4_A, CH4_R, CH4_ACH4zp

    n_Layers = Predictor%n_Layers

    Predictor%Secant_Zenith = Secang

    !------------------------------------
    ! Compute integrated variables
    !------------------------------------

    Tzp_ref   = ZERO
    Tzp_sum   = ZERO
    Tz_ref    = ZERO
    Tz_sum    = ZERO
    GAz_ref   = ZERO
    GAz_sum   = ZERO
    GAzp_ref  = ZERO
    GAzp_sum  = ZERO
    GATzp_ref = ZERO
    GATzp_sum = ZERO

    Layer_Loop : DO k = 1, n_Layers

      ! weight for integrated variables
      if(k == 1)then
            PDP = Ref_Level_Pressure(0) * &
                 ( Ref_Level_Pressure(1) - Ref_Level_Pressure(0) )

      else
            PDP = Ref_Level_Pressure(k) * &
                 ( Ref_Level_Pressure(k) - Ref_Level_Pressure(k-1) )

      endif

      ! Temperature
      Tz_ref  = Tz_ref + Ref_Temperature(k)
      Tz_sum  = Tz_sum + Temperature(k)
      Tz(k)   = Tz_sum / Tz_ref
      Tzp_ref = Tzp_ref + PDP * Ref_Temperature(k)
      Tzp_sum = Tzp_sum + PDP*Temperature(k)
      Tzp(k)  = Tzp_sum/Tzp_ref

      ! absorbers
      DO j = 1, SIZE(Absorber, DIM=2)
        GAz_ref(j)   = GAz_ref(j) + Ref_absorber(k, j)
        GAz_sum(j)   = GAz_sum(j) + Absorber(k, j)
        GAz(k, j)    = GAz_sum(j) / GAz_ref(j)
        GAzp_ref(j)  = GAzp_ref(j) + PDP*Ref_absorber(k, j)
        GAzp_sum(j)  = GAzp_sum(j) + PDP*Absorber(k, j)
        GAzp(k, j)   = GAzp_sum(j) / GAzp_ref(j)
        GATzp_ref(j) = GATzp_ref(j) + PDP*Ref_absorber(k, j)*Ref_Temperature(k)
        GATzp_sum(j) = GATzp_sum(j) + PDP*Absorber(k, j)*Temperature(k)
        GATzp(k, j)  = GATzp_sum(j) / GATzp_ref(j)
      END DO

      ! save FW variables for TL and AD routines
      IF ( PAFV_Associated(Predictor%PAFV) ) THEN
        Predictor%PAFV%PDP(k)          = PDP
        Predictor%PAFV%Tz_ref(k)       = Tz_ref
        Predictor%PAFV%Tz(k)           = Tz(k)
        Predictor%PAFV%Tzp_ref(k)      = Tzp_ref
        Predictor%PAFV%Tzp(k)          = Tzp(k)
        Predictor%PAFV%GAz_ref(k, :)   = GAz_ref
        Predictor%PAFV%GAz_sum(k, :)   = GAz_sum
        Predictor%PAFV%GAz(k, :)       = GAz(k, :)
        Predictor%PAFV%GAzp_ref(k, :)  = GAzp_ref
        Predictor%PAFV%GAzp_sum(k, :)  = GAzp_sum
        Predictor%PAFV%GAzp(k, :)      = GAzp(k, :)
        Predictor%PAFV%GATzp_ref(k, :) = GATzp_ref
        Predictor%PAFV%GATzp_sum(k, :) = GATzp_sum
        Predictor%PAFV%GATzp(k, :)     = GATzp(k, :)
      END IF

    END DO Layer_Loop

    !----------------------------------------------------------------
    ! Per-component predictor computation. The registry supplies the
    ! roster; kernels are dispatched per component ID within the basis
    ! layer loop. Forward X assignments are independent of one another,
    ! so kernel order does not affect results.
    !----------------------------------------------------------------

    ! Number of predictors per component (kernel capability). has_trace
    ! selects the group-1 style WLO/CO2 variants and must be set first.
    has_trace = ANY( Component_ID == CO_ComID )   ! validation guarantees the trio
    DO ic = 1, SIZE(Component_ID)
      Predictor%n_CP(ic) = ODPS_Kernel_n_Predictors( &
        GROUP_REGISTRY(Group_ID)%Basis, Component_ID(ic), has_trace )
    END DO

    ! Resolve each gas's position in this group's absorber roster
    ! (0 when the gas is not carried; its kernel is then never dispatched)
    ja_h2o = Absorber_Position(H2O_ID)
    ja_o3  = Absorber_Position(O3_ID)
    ja_co2 = Absorber_Position(CO2_ID)
    ja_n2o = Absorber_Position(N2O_ID)
    ja_co  = Absorber_Position(CO_ID)
    ja_ch4 = Absorber_Position(CH4_ID)
    ja_no2 = Absorber_Position(NO2_ID)

    ! Silence gfortran complaints about maybe-used-uninit by init to HUGE()
    N2O_S       = HUGE(N2O_S)
    N2O_R       = HUGE(N2O_R)
    N2O_A       = HUGE(N2O_A)
    N2O         = HUGE(N2O)
    CO_S        = HUGE(CO_S)
    CO_R        = HUGE(CO_R)
    CO_ACOdCOzp = HUGE(CO_ACOdCOzp)
    CO_A        = HUGE(CO_A)
    CH4_R       = HUGE(CH4_R)
    CH4_ACH4zp  = HUGE(CH4_ACH4zp)
    CH4_A       = HUGE(CH4_A)
    CH4         = HUGE(CH4)

    Basis_Select: IF ( GROUP_REGISTRY(Group_ID)%Basis == BASIS_IR ) THEN

      IR_Layer_Loop : DO k = 1, n_Layers

        !------------------------------------------
        !  Relative Temperature
        !------------------------------------------
        dT = Temperature(k) - Ref_Temperature(k)
        T = Temperature(k) / Ref_Temperature(k)

        !-------------------------------------------
        !  Abosrber amount scalled by the reference
        !-------------------------------------------
        IF ( ja_h2o > 0 ) H2O = Absorber(k,ja_h2o)/Ref_Absorber(k, ja_h2o)
        IF ( ja_o3  > 0 ) O3  = Absorber(k,ja_o3)/Ref_absorber(k,ja_o3)
        IF ( ja_co2 > 0 ) CO2 = Absorber(k,ja_co2)/Ref_absorber(k,ja_co2)

        ! Combinations of variables common to all predictor groups
        T2   = T*T
        DT2  = DT*ABS( DT )

        IF ( ja_h2o > 0 ) THEN
          H2O_A = SECANG(k)*H2O
          H2O_R  = SQRT( H2O_A )
          H2O_S  = H2O_A*H2O_A
          H2O_R4 = SQRT( H2O_R )
          H2OdH2OTzp = H2O/GATzp(k, ja_h2o)
        END IF

        IF ( ja_o3 > 0 ) THEN
          O3_A = SECANG(k)*O3
          O3_R = SQRT( O3_A )
        END IF

        IF( has_trace )THEN
          CO  = Absorber(k,ja_co)/Ref_absorber(k, ja_co)
          N2O = Absorber(k,ja_n2o)/Ref_absorber(k,ja_n2o)
          CH4 = Absorber(k,ja_ch4)/Ref_absorber(k,ja_ch4)

          N2O_A = SECANG(k)*N2O
          N2O_R = SQRT( N2O_A )
          N2O_S = N2O_A*N2O_A

          CO_A = SECANG(k)*CO
          CO_R = SQRT( CO_A )
          CO_S = CO_A*CO_A
          CO_ACOdCOzp = CO_A*CO/GAzp(k, ja_co)

          CH4_A = SECANG(k)*CH4
          CH4_R = SQRT(CH4_A)
          CH4_ACH4zp = SECANG(k)*GAzp(k, ja_ch4)
        END IF

        IR_Component_Loop : DO ic = 1, SIZE(Component_ID)
          np = Predictor%n_CP(ic)
          SELECT CASE ( Component_ID(ic) )
            CASE ( DRY_ComID_G1, DRY_ComID_G2, DRY_ComID_Q1, DRY_ComID_Q2 )
              CALL FWD_Kernel_DRY(k, ic, np)
            CASE ( WLO_ComID )
              CALL FWD_Kernel_WLO(k, ic, np)
            CASE ( WCO_ComID )
              CALL FWD_Kernel_WCO(k, ic)
            CASE ( OZO_ComID )
              CALL FWD_Kernel_OZO(k, ic)
            CASE ( CO2_ComID )
              CALL FWD_Kernel_CO2(k, ic, np)
            CASE ( N2O_ComID )
              CALL FWD_Kernel_N2O(k, ic)
            CASE ( CO_ComID )
              CALL FWD_Kernel_CO(k, ic)
            CASE ( CH4_ComID )
              CALL FWD_Kernel_CH4(k, ic)
            CASE ( NO2_ComID )
              CALL FWD_Kernel_NO2(k, ic)
          END SELECT
        END DO IR_Component_Loop

      END DO IR_Layer_Loop

    ELSE Basis_Select   ! BASIS_MW

      MW_Layer_Loop : DO k = 1, n_Layers

        !------------------------------------------
        !  Relative Temperature
        !------------------------------------------
        dT = Temperature(k) - Ref_Temperature(k)
        T = Temperature(k) / Ref_Temperature(k)

        !-------------------------------------------
        !  Abosrber amount scalled by the reference
        !-------------------------------------------
        ! Combinations of variables common to all predictor groups
        T2  = T*T
        DT2 = DT*ABS( DT )

        IF ( ja_h2o > 0 ) THEN
          H2O = Absorber(k,ja_h2o)/Ref_Absorber(k, ja_h2o)
          H2O_A = SECANG(k)*H2O
          H2O_R  = SQRT( H2O_A )
          H2O_S  = H2O_A*H2O_A
          H2O_R4 = SQRT( H2O_R )
          H2OdH2OTzp = H2O/GATzp(k, ja_h2o)
        END IF

        IF ( ja_o3 > 0 ) THEN
          O3   = Absorber(k,ja_o3)/Ref_Absorber(k, ja_o3)
          O3_A = SECANG(k)*O3
          O3_R = SQRT( O3_A )
        END IF

        MW_Component_Loop : DO ic = 1, SIZE(Component_ID)
          np = Predictor%n_CP(ic)
          SELECT CASE ( Component_ID(ic) )
            CASE ( EDRY_ComID )
              CALL FWD_Kernel_DRY(k, ic, np)
            CASE ( WET_ComID )
              CALL FWD_Kernel_WET_MW(k, ic)
            CASE ( OZO_ComID )
              CALL FWD_Kernel_OZO(k, ic)
          END SELECT
        END DO MW_Component_Loop

      END DO MW_Layer_Loop

    END IF Basis_Select

CONTAINS

    ! Position of a HITRAN absorber ID in the file's absorber roster
    PURE FUNCTION Absorber_Position( Gas_ID ) RESULT( Position )
      INTEGER, INTENT(IN) :: Gas_ID
      INTEGER :: Position
      INTEGER :: ja
      Position = 0
      DO ja = 1, SIZE(Absorber_ID)
        IF ( Absorber_ID(ja) == Gas_ID ) THEN
          Position = ja
          RETURN
        END IF
      END DO
    END FUNCTION Absorber_Position

    !  ----------------------
    !   Fixed (Dry) predictors (IR and MW use the same formulation)
    !  ----------------------
    SUBROUTINE FWD_Kernel_DRY( k, ic, np )
      INTEGER, INTENT(IN) :: k, ic, np
      Predictor%X(k, 1, ic)  = SECANG(k)
      Predictor%X(k, 2, ic)  = SECANG(k) * T
      Predictor%X(k, 3, ic)  = SECANG(k) * T2
      Predictor%X(k, 4, ic)  = T
      Predictor%X(k, 5, ic)  = SECANG(k) * SECANG(k)
      Predictor%X(k, 6, ic)  = T2
      Predictor%X(k, 7, ic)  = Tz(k)
      ! QDRY dry components (123/124) carry water-vapor predictors 1 - 13
      ! of the WLO formulation in slots 8 - 20: the band-averaged
      ! effective-dry target responds to humidity through the in-band
      ! water/dry-gas line-overlap correlation term
      IF ( np >= 20 ) THEN
        Predictor%X(k, 8, ic)  = H2O_A
        Predictor%X(k, 9, ic)  = H2O_A*DT
        Predictor%X(k,10, ic)  = H2O_S
        Predictor%X(k,11, ic)  = H2O_A*DT2
        Predictor%X(k,12, ic)  = H2O_R4
        Predictor%X(k,13, ic)  = H2O_S*H2O_A
        Predictor%X(k,14, ic)  = H2O_R
        Predictor%X(k,15, ic)  = H2O_R*DT
        Predictor%X(k,16, ic)  = H2O_S*H2O_S
        Predictor%X(k,17, ic)  = H2OdH2OTzp
        Predictor%X(k,18, ic)  = H2O_R*H2OdH2OTzp
        Predictor%X(k,19, ic)  = (SECANG(k)*GAzp(k,ja_h2o))**2
        Predictor%X(k,20, ic)  = SECANG(k)*GAzp(k,ja_h2o)
      END IF
    END SUBROUTINE FWD_Kernel_DRY

    !  --------------------------
    !  Water vapor continuum predictors
    !  --------------------------
    SUBROUTINE FWD_Kernel_WCO( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      Predictor%X(k, 1, ic) = H2O_A/T
      Predictor%X(k, 2, ic) = H2O_A/T * H2O
      Predictor%X(k, 3, ic) = H2O_A/T2 * H2O/T2
      Predictor%X(k, 4, ic) = H2O_A/T2
      Predictor%X(k, 5, ic) = H2O_A/T2 * H2O
      Predictor%X(k, 6, ic) = H2O_A/T2**2
      Predictor%X(k, 7, ic) = H2O_A
    END SUBROUTINE FWD_Kernel_WCO

    !  -----------------------
    !  Ozone predictors (same formulation for the IR and MW_O3 groups)
    !  -----------------------
    SUBROUTINE FWD_Kernel_OZO( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      Predictor%X(k, 1, ic)  = O3_A
      Predictor%X(k, 2, ic)  = O3_A*DT
      Predictor%X(k, 3, ic)  = O3_A*O3*GAzp(k,ja_o3)
      Predictor%X(k, 4, ic)  = O3_A*O3_A
      Predictor%X(k, 5, ic)  = O3_A*GAzp(k,ja_o3)
      Predictor%X(k, 6, ic)  = O3_A*SQRT(SECANG(k)*GAzp(k,ja_o3))
      Predictor%X(k, 7, ic)  = O3_R*DT   !T*T*T
      Predictor%X(k, 8, ic)  = O3_R
      Predictor%X(k, 9, ic)  = O3_R*O3/GAzp(k,ja_o3)
      Predictor%X(k,10, ic)  = SECANG(k)*GAzp(k,ja_o3)
      Predictor%X(k,11, ic)  = (SECANG(k)*GAzp(k,ja_o3))**2
    END SUBROUTINE FWD_Kernel_OZO

    !  -----------------------
    !  Carbon dioxide predictors; predictor 11 (CO amount) is carried only
    !  by rosters that request 11 predictors for CO2 (group 1)
    !  -----------------------
    SUBROUTINE FWD_Kernel_CO2( k, ic, np )
      INTEGER, INTENT(IN) :: k, ic, np
      Predictor%X(k, 1, ic)  = SECANG(k) * T
      Predictor%X(k, 2, ic)  = SECANG(k) * T2
      Predictor%X(k, 3, ic)  = T
      Predictor%X(k, 4, ic)  = T2
      Predictor%X(k, 5, ic)  = SECANG(k)
      Predictor%X(k, 6, ic)  = SECANG(k)*CO2
      Predictor%X(k, 7, ic)  = SECANG(k) * Tzp(k)
      Predictor%X(k, 8, ic)  = (SECANG(k) * GAzp(k, ja_co2))**2
      Predictor%X(k, 9, ic)  = Tzp(k)**3
      Predictor%X(k, 10, ic) = SECANG(k) * Tzp(k) * SQRT(T)
      IF ( np >= 11 ) THEN
        Predictor%X(k, 11, ic) = CO_A
      END IF
    END SUBROUTINE FWD_Kernel_CO2

    !  --------------------------
    !  Water-line predictors; predictors 16 - 18 (CH4/CO cross terms) are
    !  carried only by rosters that request 18 predictors for WLO (group 1)
    !  --------------------------
    SUBROUTINE FWD_Kernel_WLO( k, ic, np )
      INTEGER, INTENT(IN) :: k, ic, np
      Predictor%X(k, 1, ic) = H2O_A
      Predictor%X(k, 2, ic) = H2O_A*DT
      Predictor%X(k, 3, ic) = H2O_S
      Predictor%X(k, 4, ic) = H2O_A*DT2
      Predictor%X(k, 5, ic) = H2O_R4
      Predictor%X(k, 6, ic) = H2O_S*H2O_A
      Predictor%X(k, 7, ic) = H2O_R
      Predictor%X(k, 8, ic) = H2O_R*DT
      Predictor%X(k, 9, ic) = H2O_S*H2O_S
      Predictor%X(k,10, ic) = H2OdH2OTzp
      Predictor%X(k,11, ic) = H2O_R*H2OdH2OTzp
      Predictor%X(k,12, ic) = (SECANG(k)*GAzp(k,ja_h2o))**2
      Predictor%X(k,13, ic) = SECANG(k)*GAzp(k,ja_h2o)
      Predictor%X(k,14, ic) = SECANG(k)
      Predictor%X(k,15, ic) = SECANG(k) * CO2
      IF ( np >= 18 ) THEN
        Predictor%X(k,16, ic) = CH4_A
        Predictor%X(k,17, ic) = CH4_A*CH4_A*DT
        Predictor%X(k,18, ic) = CO_A
      END IF
    END SUBROUTINE FWD_Kernel_WLO

    !  -----------------------
    !  Carbon monoxide
    !  -----------------------
    SUBROUTINE FWD_Kernel_CO( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      Predictor%X(k, 1,  ic)   = CO_A
      Predictor%X(k, 2,  ic)   = CO_A*DT
      Predictor%X(k, 3,  ic)   = SQRT( CO_R )
      Predictor%X(k, 4,  ic)   = CO_R*DT
      Predictor%X(k, 5,  ic)   = CO_S
      Predictor%X(k, 6,  ic)   = CO_R
      Predictor%X(k, 7,  ic)   = CO_A*DT2
      Predictor%X(k, 8,  ic)   = CO_ACOdCOzp
      Predictor%X(k, 9,  ic)   = CO_ACOdCOzp/CO_R
      Predictor%X(k, 10, ic)   = CO_ACOdCOzp * SQRT( GAzp(k, ja_co) )
    END SUBROUTINE FWD_Kernel_CO

    !  -----------------------
    !  Methane predictors
    !  -----------------------
    SUBROUTINE FWD_Kernel_CH4( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      Predictor%X(k, 1,  ic)  = CH4_A*DT
      Predictor%X(k, 2,  ic)  = CH4_R
      Predictor%X(k, 3,  ic)  = CH4_A*CH4_A
      Predictor%X(k, 4,  ic)  = CH4_A
      Predictor%X(k, 5,  ic)  = CH4*DT
      Predictor%X(k, 6,  ic)  = CH4_ACH4zp
      Predictor%X(k, 7,  ic)  = CH4_ACH4zp**2
      Predictor%X(k, 8,  ic)  = SQRT(CH4_R)
      Predictor%X(k, 9,  ic)  = GATzp(k, ja_ch4)
      Predictor%X(k, 10, ic)  = SECANG(k)*GATzp(k, ja_ch4)
      Predictor%X(k, 11, ic)  = CH4_R * CH4/GAzp(k, ja_ch4)
    END SUBROUTINE FWD_Kernel_CH4

    !  -----------------------
    !    N2O predictors
    !  -----------------------
    SUBROUTINE FWD_Kernel_N2O( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      Predictor%X(k, 1, ic)   = N2O_A*DT
      Predictor%X(k, 2, ic)   = N2O_R
      Predictor%X(k, 3, ic)   = N2O*DT
      Predictor%X(k, 4, ic)   = N2O_A**POINT_25
      Predictor%X(k, 5, ic)   = N2O_A
      Predictor%X(k, 6, ic)   = SECANG(k) * GAzp(k, ja_n2o)
      Predictor%X(k, 7, ic)   = SECANG(k) * GATzp(k, ja_n2o)
      Predictor%X(k, 8, ic)   = N2O_S
      Predictor%X(k, 9, ic)   = GATzp(k, ja_n2o)
      Predictor%X(k,10, ic)   = N2O_R*N2O / GAzp(k, ja_n2o)
      Predictor%X(k,11, ic)   = CH4_A
      Predictor%X(k,12, ic)   = CH4_A*GAzp(k, ja_ch4)
      Predictor%X(k,13, ic)   = CO_A
      Predictor%X(k,14, ic)   = CO_A*SECANG(k)*GAzp(k, ja_co)
    END SUBROUTINE FWD_Kernel_N2O

    !  -----------------------
    !  NO2 predictors (GROUP_UV_NO2). Scene NO2 in the UV/VIS is pure
    !  Beer-Lambert electronic-cross-section extinction: layer OD = sigma(T)*N,
    !  exactly linear in amount. A compact set suffices: amount*secant, plus
    !  DT and DT2 terms for the smooth (~3-12%) sigma(T) temperature dependence.
    !  -----------------------
    SUBROUTINE FWD_Kernel_NO2( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      NO2   = Absorber(k,ja_no2)/Ref_Absorber(k, ja_no2)
      NO2_A = SECANG(k)*NO2
      Predictor%X(k, 1, ic) = NO2_A
      Predictor%X(k, 2, ic) = NO2_A*DT
      Predictor%X(k, 3, ic) = NO2_A*DT2
    END SUBROUTINE FWD_Kernel_NO2

    !  --------------------------------
    !  Water vapor, MW (line and continuum together)
    !  --------------------------------
    SUBROUTINE FWD_Kernel_WET_MW( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      Predictor%X(k, 1, ic) = H2O_A/T
      Predictor%X(k, 2, ic) = H2O_A/T * H2O
      Predictor%X(k, 3, ic) = H2O_A/T2 * H2O/T2
      Predictor%X(k, 4, ic) = H2O_A/T2
      Predictor%X(k, 5, ic) = H2O_A/T2 * H2O
      Predictor%X(k, 6, ic) = H2O_A/T2**2
      Predictor%X(k, 7, ic) = H2O_A
      Predictor%X(k, 8, ic) = H2O_A*DT
      Predictor%X(k, 9, ic) = (SECANG(k)*GAzp(k,ja_h2o))**2
      Predictor%X(k, 10,ic) = SECANG(k)*GAzp(k,ja_h2o)
      Predictor%X(k, 11,ic) = SECANG(k)
      Predictor%X(k, 12,ic) = H2O_S*H2O_A
      Predictor%X(k, 13,ic) = H2O_S*H2O_S
      Predictor%X(k, 14,ic) = H2OdH2OTzp
    END SUBROUTINE FWD_Kernel_WET_MW

  END SUBROUTINE  ODPS_Compute_Predictor

!============================== TL
!------------------------------------------------------------------------------
!
! NAME:
!       ODPS_Compute_Predictor_TL
!
! PURPOSE:
!       Subroutine to predictors
!
! CALLING SEQUENCE:
!       CALL ODPS_Compute_Predictor_TL( &
!         Group_ID,           &  ! Input
!         Temperature,        &  ! Input
!         Absorber,           &  ! Input
!         Ref_Temperature,    &  ! Input
!         Ref_Absorber,       &  ! Input
!         secang,             &  ! Input
!         Predictor           &  ! Input
!         Temperature_TL,     &  ! Input
!         Absorber_TL,        &  ! Input
!         Predictor_TL )         ! Output
!
! INPUT ARGUMENTS:
!       Group_ID   :     The ID of predictor group
!                        UNITS:      N?A
!                        TYPE:       INTEGER
!                        DIMENSION:  scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Temperature:     Temperature profile
!                        UNITS:      K
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!       Absorber   :     Absorber profiles
!                        UNITS:      vary
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-2(n_Layers x n_Absorbers) array
!                        ATTRIBUTES: INTENT(IN)
!
!       Ref_Temperature : Reference layer temperature profile
!                        UNITS:      K
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!       Ref_Absorber :   Reference absorber profiles
!                        UNITS:      vary
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-2(n_Layers x n_Absorbers) array
!                        ATTRIBUTES: INTENT(IN)
!
!       secang       :   Secont sensor zenith angle profile
!                        UNITS:      N/A
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!       Predictor:      Predictor structure containing the integrated absorber
!                       and predictor profiles.
!                       UNITS:      N/A
!                       TYPE:       ODPS_Predictor_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Temperature_TL:  Temperature_TL profile
!                        UNITS:      K
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!       Absorber_TL:     Absorber_TL profiles
!                        UNITS:      vary
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-2(n_Layers x n_Absorbers) array
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!
!       Predictor_TL:   Predictor_TL structure containing the integrated absorber_TL
!                       and predictor_TL profiles.
!                       UNITS:      N/A
!                       TYPE:       ODPS_Predictor_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
!------------------------------------------------------------------------------

  SUBROUTINE ODPS_Compute_Predictor_TL( &
    Group_ID,           &
    Component_ID,       &
    Absorber_ID,        &
    Temperature,        &
    Absorber,           &
    Ref_Temperature,    &
    Ref_Absorber,       &
    secang,             &
    Predictor,          &
    Temperature_TL,     &
    Absorber_TL,        &
    Predictor_TL )

    INTEGER,                           INTENT(IN)     :: Group_ID
    INTEGER,                           INTENT(IN)     :: Component_ID(:)
    INTEGER,                           INTENT(IN)     :: Absorber_ID(:)
    REAL(fp),                          INTENT(IN)     :: Temperature(:)
    REAL(fp),                          INTENT(IN)     :: Absorber(:, :)
    REAL(fp),                          INTENT(IN)     :: Ref_Temperature(:)
    REAL(fp),                          INTENT(IN)     :: Ref_Absorber(:, :)
    REAL(fp),                          INTENT(IN)     :: Secang(:)
    REAL(fp),                          INTENT(IN)     :: Temperature_TL(:)
    REAL(fp),                          INTENT(IN)     :: Absorber_TL(:, :)
    TYPE(ODPS_Predictor_type), TARGET, INTENT(IN)     :: Predictor
    TYPE(ODPS_Predictor_type),         INTENT(IN OUT) :: Predictor_TL

    ! ---------------
    ! Local variables
    ! ---------------
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'ODPS_Compute_Predictor_TL'
    INTEGER    ::    n_layers
    INTEGER    ::    k ! n_Layers, n_Levels
    INTEGER    ::    j ! n_absorbers
    REAL(fp) ::    Tzp_sum_TL
    REAL(fp) ::    Tzp_TL(SIZE(Absorber, DIM=1))
    REAL(fp) ::    Tz_sum_TL
    REAL(fp) ::    Tz_TL(SIZE(Absorber, DIM=1))
    REAL(fp) ::    GAz_sum_TL(SIZE(Absorber, DIM=2))
    REAL(fp) ::    GAz_TL(SIZE(Absorber, DIM=1), SIZE(Absorber, DIM=2))
    REAL(fp) ::    GAzp_sum_TL(SIZE(Absorber, DIM=2))
    REAL(fp) ::    GAzp_TL(SIZE(Absorber, DIM=1), SIZE(Absorber, DIM=2))
    REAL(fp) ::    GATzp_sum_TL(SIZE(Absorber, DIM=2))
    REAL(fp) ::    GATzp_TL(SIZE(Absorber, DIM=1), SIZE(Absorber, DIM=2))
    ! Kernel dispatch bookkeeping and shared per-layer variables
    INTEGER  :: ic, np
    INTEGER  :: ja_h2o, ja_o3, ja_co2, ja_n2o, ja_co, ja_ch4, ja_no2
    LOGICAL  :: has_trace
    REAL(fp) :: DT, DT_TL, T, T_TL, T2, T2_TL, DT2, DT2_TL
    REAL(fp) :: H2O, H2O_TL, H2O_A, H2O_A_TL, H2O_R, H2O_R_TL
    REAL(fp) :: H2O_S, H2O_S_TL, H2O_R4, H2O_R4_TL, H2OdH2OTzp, H2OdH2OTzp_TL
    REAL(fp) :: CO2, CO2_TL, O3, O3_TL, O3_A, O3_A_TL, O3_R, O3_R_TL
    REAL(fp) :: NO2, NO2_TL, NO2_A, NO2_A_TL
    REAL(fp) :: CO, CO_TL, CO_A, CO_A_TL, CO_R, CO_R_TL, CO_S, CO_S_TL
    REAL(fp) :: CO_ACOdCOzp, CO_ACOdCOzp_TL
    REAL(fp) :: N2O, N2O_TL, N2O_A, N2O_A_TL, N2O_R, N2O_R_TL, N2O_S, N2O_S_TL
    REAL(fp) :: CH4, CH4_TL, CH4_A, CH4_A_TL, CH4_R, CH4_R_TL
    REAL(fp) :: CH4_ACH4zp, CH4_ACH4zp_TL
!JR Static initialization means only 1 copy of the variable. OpenMP over profiles
!JR means $OPENMP_NUM_THREADS copies are needed. So change to run-time initialization
!JR    TYPE(PAFV_type), POINTER  :: PAFV => NULL()
    TYPE(PAFV_type), POINTER  :: PAFV

    NULLIFY(PAFV)

    ! use short name
    PAFV => Predictor%PAFV

    n_Layers = Predictor%n_Layers

    Predictor_TL%Secant_Zenith = Secang

    !------------------------------------
    ! Compute integrated variables
    !------------------------------------

    Tzp_sum_TL   = ZERO
    Tz_sum_TL    = ZERO
    GAz_sum_TL   = ZERO
    GAzp_sum_TL  = ZERO
    GATzp_sum_TL = ZERO

    Layer_Loop : DO k = 1, n_Layers

      ! Temperature
      Tz_sum_TL  = Tz_sum_TL + Temperature_TL(k)
      Tz_TL(k)   = Tz_sum_TL / PAFV%Tz_ref(k)
      Tzp_sum_TL = Tzp_sum_TL + PAFV%PDP(k)*Temperature_TL(k)
      Tzp_TL(k)  = Tzp_sum_TL/PAFV%Tzp_ref(k)

      ! absorbers
      DO j = 1, SIZE(Absorber, DIM=2)
        GAz_sum_TL(j)   = GAz_sum_TL(j) + Absorber_TL(k, j)
        GAz_TL(k, j)    = GAz_sum_TL(j) / PAFV%GAz_ref(k,j)
        GAzp_sum_TL(j)  = GAzp_sum_TL(j) + PAFV%PDP(k)*Absorber_TL(k, j)
        GAzp_TL(k, j)   = GAzp_sum_TL(j) / PAFV%GAzp_ref(k,j)
        GATzp_sum_TL(j) = GATzp_sum_TL(j) + PAFV%PDP(k)*Absorber_TL(k, j)*Temperature(k) + &
                          PAFV%PDP(k)*Absorber(k, j)*Temperature_TL(k)
        GATzp_TL(k, j)  = GATzp_sum_TL(j) / PAFV%GATzp_ref(k,j)
      END DO

    END DO Layer_Loop

    !----------------------------------------------------------------
    ! Per-component tangent-linear predictor computation; mirrors the
    ! forward kernel dispatch (TL X assignments are independent of one
    ! another, so kernel order does not affect results).
    !----------------------------------------------------------------

    ! Number of predictors per component (kernel capability). has_trace
    ! selects the group-1 style WLO/CO2 variants and must be set first.
    has_trace = ANY( Component_ID == CO_ComID )   ! validation guarantees the trio
    DO ic = 1, SIZE(Component_ID)
      Predictor_TL%n_CP(ic) = ODPS_Kernel_n_Predictors( &
        GROUP_REGISTRY(Group_ID)%Basis, Component_ID(ic), has_trace )
    END DO

    ! Resolve each gas's position in this group's absorber roster
    ja_h2o = Absorber_Position(H2O_ID)
    ja_o3  = Absorber_Position(O3_ID)
    ja_co2 = Absorber_Position(CO2_ID)
    ja_n2o = Absorber_Position(N2O_ID)
    ja_co  = Absorber_Position(CO_ID)
    ja_ch4 = Absorber_Position(CH4_ID)
    ja_no2 = Absorber_Position(NO2_ID)

    ! Silence gfortran complaints about maybe-used-uninit by init to HUGE()
    N2O_TL        = HUGE(N2O_TL)
    N2O_S_TL      = HUGE(N2O_S_TL)
    N2O_S         = HUGE(N2O_S)
    N2O_R         = HUGE(N2O_R)
    N2O_R_TL      = HUGE(N2O_R_TL)
    N2O           = HUGE(N2O)
    N2O_A         = HUGE(N2O_A)
    N2O_A_TL      = HUGE(N2O_A_TL)
    CO_S_TL       = HUGE(CO_S_TL)
    CO_S          = HUGE(CO_S)
    CO_R          = HUGE(CO_R)
    CO_R_TL       = HUGE(CO_R_TL)
    CO_A_TL       = HUGE(CO_A_TL)
    CO_A          = HUGE(CO_A)
    CO_ACODCOZP_TL= HUGE(CO_ACODCOZP_TL)
    CO_ACODCOZP   = HUGE(CO_ACODCOZP)
    CH4_TL        = HUGE(CH4_TL)
    CH4_R_TL      = HUGE(CH4_R_TL)
    CH4_A_TL      = HUGE(CH4_A_TL)
    CH4_A         = HUGE(CH4_A)
    CH4_R         = HUGE(CH4_R)
    CH4           = HUGE(CH4)
    CH4_ACH4ZP_TL = HUGE(CH4_ACH4ZP_TL)
    CH4_ACH4ZP    = HUGE(CH4_ACH4ZP)

    Basis_Select: IF ( GROUP_REGISTRY(Group_ID)%Basis == BASIS_IR ) THEN

      IR_Layer_Loop : DO k = 1, n_Layers

        !------------------------------------------
        !  Relative Temperature
        !------------------------------------------
        dT = Temperature(k) - Ref_Temperature(k)
        T = Temperature(k) / Ref_Temperature(k)
        dT_TL = Temperature_TL(k)
        T_TL = Temperature_TL(k) / Ref_Temperature(k)

        !-------------------------------------------
        !  Abosrber amount scalled by the reference
        !-------------------------------------------
        IF ( ja_h2o > 0 ) THEN
          H2O    = Absorber(k,ja_h2o)/Ref_Absorber(k, ja_h2o)
          H2O_TL = Absorber_TL(k,ja_h2o)/Ref_Absorber(k, ja_h2o)
        END IF
        IF ( ja_o3 > 0 ) THEN
          O3    = Absorber(k,ja_o3)/Ref_absorber(k,ja_o3)
          O3_TL = Absorber_TL(k,ja_o3)/Ref_absorber(k,ja_o3)
        END IF
        IF ( ja_co2 > 0 ) THEN
          CO2    = Absorber(k,ja_co2)/Ref_absorber(k,ja_co2)
          CO2_TL = Absorber_TL(k,ja_co2)/Ref_absorber(k,ja_co2)
        END IF

        ! Combinations of variables common to all predictor groups
        T2   = T*T
        DT2  = DT*ABS( DT )

        T2_TL  = TWO*T*T_TL
        IF( DT > ZERO) THEN
          DT2_TL = TWO*DT*DT_TL
        ELSE
          DT2_TL = - TWO*DT*DT_TL
        ENDIF

        IF ( ja_h2o > 0 ) THEN
          H2O_A  = SECANG(k)*H2O
          H2O_R  = SQRT( H2O_A )
          H2O_S  = H2O_A*H2O_A
          H2O_R4 = SQRT( H2O_R )
          H2OdH2OTzp = H2O/PAFV%GATzp(k, ja_h2o)

          H2O_A_TL  = SECANG(k)*H2O_TL
          H2O_R_TL  = (POINT_5 / SQRT(H2O_A)) * H2O_A_TL
          H2O_S_TL  = TWO * H2O_A * H2O_A_TL
          H2O_R4_TL = (POINT_5 / SQRT(H2O_R)) * H2O_R_TL
          H2OdH2OTzp_TL = H2O_TL/PAFV%GATzp(k, ja_h2o) - &
                          H2O * GATzp_TL(k, ja_h2o)/PAFV%GATzp(k, ja_h2o)**2
        END IF

        IF ( ja_o3 > 0 ) THEN
          O3_A = SECANG(k)*O3
          O3_R = SQRT( O3_A )

          O3_A_TL = SECANG(k)*O3_TL
          O3_R_TL = (POINT_5 / SQRT(O3_A)) * O3_A_TL
        END IF

        IF( has_trace )THEN
          CO  = Absorber(k,ja_co)/Ref_absorber(k, ja_co)
          N2O = Absorber(k,ja_n2o)/Ref_absorber(k,ja_n2o)
          CH4 = Absorber(k,ja_ch4)/Ref_absorber(k,ja_ch4)

          CO_TL  = Absorber_TL(k,ja_co)/Ref_absorber(k, ja_co)
          N2O_TL = Absorber_TL(k,ja_n2o)/Ref_absorber(k,ja_n2o)
          CH4_TL = Absorber_TL(k,ja_ch4)/Ref_absorber(k,ja_ch4)

          N2O_A = SECANG(k)*N2O
          N2O_R = SQRT( N2O_A )
          N2O_S = N2O_A*N2O_A

          N2O_A_TL = SECANG(k) * N2O_TL
          N2O_R_TL = (POINT_5 / SQRT(N2O_A)) * N2O_A_TL
          N2O_S_TL = TWO * N2O_A * N2O_A_TL

          CO_A = SECANG(k)*CO
          CO_R = SQRT( CO_A )
          CO_S = CO_A*CO_A
          CO_ACOdCOzp = CO_A*CO/PAFV%GAzp(k, ja_co)

          CO_A_TL = SECANG(k)*CO_TL
          CO_R_TL = (POINT_5 / SQRT(CO_A)) * CO_A_TL
          CO_S_TL = TWO * CO_A * CO_A_TL
          CO_ACOdCOzp_TL = CO_A_TL*CO/PAFV%GAzp(k, ja_co) + CO_A*CO_TL/PAFV%GAzp(k, ja_co) &
                           - CO_A*CO*GAzp_TL(k, ja_co)/PAFV%GAzp(k, ja_co)**2

          CH4_A = SECANG(k)*CH4
          CH4_R = SQRT(CH4_A)
          CH4_ACH4zp = SECANG(k)*PAFV%GAzp(k, ja_ch4)

          CH4_A_TL = SECANG(k)*CH4_TL
          CH4_R_TL = (POINT_5 / SQRT(CH4_A)) * CH4_A_TL
          CH4_ACH4zp_TL = SECANG(k)*GAzp_TL(k, ja_ch4)
        END IF

        IR_Component_Loop : DO ic = 1, SIZE(Component_ID)
          np = Predictor_TL%n_CP(ic)
          SELECT CASE ( Component_ID(ic) )
            CASE ( DRY_ComID_G1, DRY_ComID_G2, DRY_ComID_Q1, DRY_ComID_Q2 )
              CALL TL_Kernel_DRY(k, ic, np)
            CASE ( WLO_ComID )
              CALL TL_Kernel_WLO(k, ic, np)
            CASE ( WCO_ComID )
              CALL TL_Kernel_WCO(k, ic)
            CASE ( OZO_ComID )
              CALL TL_Kernel_OZO(k, ic)
            CASE ( CO2_ComID )
              CALL TL_Kernel_CO2(k, ic, np)
            CASE ( N2O_ComID )
              CALL TL_Kernel_N2O(k, ic)
            CASE ( CO_ComID )
              CALL TL_Kernel_CO(k, ic)
            CASE ( CH4_ComID )
              CALL TL_Kernel_CH4(k, ic)
            CASE ( NO2_ComID )
              CALL TL_Kernel_NO2(k, ic)
          END SELECT
        END DO IR_Component_Loop

      END DO IR_Layer_Loop

    ELSE Basis_Select   ! BASIS_MW

      MW_Layer_Loop : DO k = 1, n_Layers

        !------------------------------------------
        !  Relative Temperature
        !------------------------------------------
        dT = Temperature(k) - Ref_Temperature(k)
        T = Temperature(k) / Ref_Temperature(k)

        dT_TL = Temperature_TL(k)
        T_TL  = Temperature_TL(k) / Ref_Temperature(k)

        !-------------------------------------------
        !  Abosrber amount scalled by the reference
        !-------------------------------------------
        IF ( ja_h2o > 0 ) THEN
          H2O = Absorber(k,ja_h2o)/Ref_Absorber(k, ja_h2o)
          H2O_TL = Absorber_TL(k,ja_h2o)/Ref_Absorber(k, ja_h2o)
        END IF

        ! Combinations of variables common to all predictor groups
        T2  = T*T
        DT2 = DT*ABS( DT )

        T2_TL  = TWO*T*T_TL
        IF( DT > ZERO) THEN
          DT2_TL = TWO*DT*DT_TL
        ELSE
          DT2_TL = - TWO*DT*DT_TL
        ENDIF

        IF ( ja_h2o > 0 ) THEN
          H2O_A = SECANG(k)*H2O
          H2O_R  = SQRT( H2O_A )
          H2O_S  = H2O_A*H2O_A
          H2O_R4 = SQRT( H2O_R )
          H2OdH2OTzp = H2O/PAFV%GATzp(k, ja_h2o)

          H2O_A_TL  = SECANG(k)*H2O_TL
          H2O_R_TL  = (POINT_5 / SQRT(H2O_A)) * H2O_A_TL
          H2O_S_TL  = TWO * H2O_A * H2O_A_TL
          H2O_R4_TL = (POINT_5 / SQRT(H2O_R)) * H2O_R_TL
          H2OdH2OTzp_TL = H2O_TL/PAFV%GATzp(k, ja_h2o) - &
                          H2O * GATzp_TL(k, ja_h2o)/PAFV%GATzp(k, ja_h2o)**2
        END IF

        IF ( ja_o3 > 0 ) THEN
          O3      = Absorber(k,ja_o3)/Ref_Absorber(k, ja_o3)
          O3_TL   = Absorber_TL(k,ja_o3)/Ref_Absorber(k, ja_o3)
          O3_A    = SECANG(k)*O3
          O3_A_TL = SECANG(k)*O3_TL
          O3_R    = SQRT( O3_A )
          O3_R_TL = (POINT_5 / SQRT(O3_A)) * O3_A_TL
        END IF

        MW_Component_Loop : DO ic = 1, SIZE(Component_ID)
          np = Predictor_TL%n_CP(ic)
          SELECT CASE ( Component_ID(ic) )
            CASE ( EDRY_ComID )
              CALL TL_Kernel_DRY(k, ic, np)
              Predictor_TL%X(:, 8:, ic) = ZERO
            CASE ( WET_ComID )
              CALL TL_Kernel_WET_MW(k, ic)
            CASE ( OZO_ComID )
              CALL TL_Kernel_OZO(k, ic)
          END SELECT
        END DO MW_Component_Loop

      END DO MW_Layer_Loop

    END IF Basis_Select

    NULLIFY(PAFV)

CONTAINS

    ! Position of a HITRAN absorber ID in the file's absorber roster
    PURE FUNCTION Absorber_Position( Gas_ID ) RESULT( Position )
      INTEGER, INTENT(IN) :: Gas_ID
      INTEGER :: Position
      INTEGER :: ja
      Position = 0
      DO ja = 1, SIZE(Absorber_ID)
        IF ( Absorber_ID(ja) == Gas_ID ) THEN
          Position = ja
          RETURN
        END IF
      END DO
    END FUNCTION Absorber_Position

    !  ----------------------
    !   Fixed (Dry) predictors (IR and MW use the same formulation)
    !  ----------------------
    SUBROUTINE TL_Kernel_DRY( k, ic, np )
      INTEGER, INTENT(IN) :: k, ic, np
      Predictor_TL%X(k, 1, ic)  = ZERO
      Predictor_TL%X(k, 2, ic)  = SECANG(k) * T_TL
      Predictor_TL%X(k, 3, ic)  = SECANG(k) * T2_TL
      Predictor_TL%X(k, 4, ic)  = T_TL
      Predictor_TL%X(k, 5, ic)  = ZERO
      Predictor_TL%X(k, 6, ic)  = T2_TL
      Predictor_TL%X(k, 7, ic)  = Tz_TL(k)
      ! QDRY water predictors 8 - 20: TL of the WLO formulas 1 - 13
      IF ( np >= 20 ) THEN
        Predictor_TL%X(k, 8, ic)  = H2O_A_TL
        Predictor_TL%X(k, 9, ic)  = H2O_A_TL*DT + H2O_A*DT_TL
        Predictor_TL%X(k,10, ic)  = H2O_S_TL
        Predictor_TL%X(k,11, ic)  = H2O_A_TL*DT2 + H2O_A*DT2_TL
        Predictor_TL%X(k,12, ic)  = H2O_R4_TL
        Predictor_TL%X(k,13, ic)  = H2O_S_TL*H2O_A + H2O_S*H2O_A_TL
        Predictor_TL%X(k,14, ic)  = H2O_R_TL
        Predictor_TL%X(k,15, ic)  = H2O_R_TL*DT + H2O_R*DT_TL
        Predictor_TL%X(k,16, ic)  = TWO*H2O_S*H2O_S_TL
        Predictor_TL%X(k,17, ic)  = H2OdH2OTzp_TL
        Predictor_TL%X(k,18, ic)  = H2O_R_TL*H2OdH2OTzp + H2O_R*H2OdH2OTzp_TL
        Predictor_TL%X(k,19, ic)  = TWO*SECANG(k)**2 * PAFV%GAzp(k,ja_h2o)*GAzp_TL(k,ja_h2o)
        Predictor_TL%X(k,20, ic)  = SECANG(k)*GAzp_TL(k,ja_h2o)
      END IF
    END SUBROUTINE TL_Kernel_DRY

    !  --------------------------
    !  Water vapor continuum predictors
    !  --------------------------
    SUBROUTINE TL_Kernel_WCO( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      Predictor_TL%X(k, 1, ic) = H2O_A_TL/T - H2O_A * T_TL/T**2
      Predictor_TL%X(k, 2, ic) = H2O_A_TL*H2O/T + H2O_A*H2O_TL/T - H2O_A*H2O*T_TL/T**2
      Predictor_TL%X(k, 3, ic) = H2O_A_TL*H2O/T2**2 + H2O_A*H2O_TL/T2**2 - &
                                 Two*H2O_A*H2O*T2_TL/T2**3
      Predictor_TL%X(k, 4, ic) = H2O_A_TL/T2 - H2O_A * T2_TL/T2**2
      Predictor_TL%X(k, 5, ic) = H2O_A_TL*H2O/T2 + H2O_A*H2O_TL/T2 - H2O_A*H2O*T2_TL/T2**2
      Predictor_TL%X(k, 6, ic) = H2O_A_TL/T2**2 - TWO*H2O_A*T2_TL/T2**3
      Predictor_TL%X(k, 7, ic) = H2O_A_TL
    END SUBROUTINE TL_Kernel_WCO

    !  -----------------------
    !  Ozone predictors (same formulation for the IR and MW_O3 groups)
    !  -----------------------
    SUBROUTINE TL_Kernel_OZO( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      Predictor_TL%X(k, 1, ic)  = O3_A_TL
      Predictor_TL%X(k, 2, ic)  = O3_A_TL*DT + O3_A*DT_TL
      Predictor_TL%X(k, 3, ic)  = O3_A_TL*O3*PAFV%GAzp(k,ja_o3) + O3_A*O3_TL*PAFV%GAzp(k,ja_o3) &
                                   + O3_A*O3*GAzp_TL(k,ja_o3)
      Predictor_TL%X(k, 4, ic)  = TWO*O3_A*O3_A_TL
      Predictor_TL%X(k, 5, ic)  = O3_A_TL*PAFV%GAzp(k,ja_o3) + O3_A*GAzp_TL(k,ja_o3)
      Predictor_TL%X(k, 6, ic)  = O3_A_TL*SQRT(SECANG(k)*PAFV%GAzp(k,ja_o3)) + &
                                   POINT_5*O3_A*SQRT(SECANG(k)/PAFV%GAzp(k,ja_o3))* &
                                   GAzp_TL(k,ja_o3)
      Predictor_TL%X(k, 7, ic)  = O3_R_TL*DT + O3_R*DT_TL  !T*T*T
      Predictor_TL%X(k, 8, ic)  = O3_R_TL
      Predictor_TL%X(k, 9, ic)  = O3_R_TL*O3/PAFV%GAzp(k,ja_o3) + O3_R*O3_TL/PAFV%GAzp(k,ja_o3) &
                                   - O3_R*O3*GAzp_TL(k,ja_o3)/PAFV%GAzp(k,ja_o3)**2
      Predictor_TL%X(k,10, ic)  = SECANG(k)*GAzp_TL(k,ja_o3)
      Predictor_TL%X(k,11, ic)  = TWO*SECANG(k)**2 * PAFV%GAzp(k,ja_o3)*GAzp_TL(k,ja_o3)
    END SUBROUTINE TL_Kernel_OZO

    !  -----------------------
    !  Carbon dioxide predictors; predictor 11 (CO amount) is carried only
    !  by rosters that request 11 predictors for CO2 (group 1)
    !  -----------------------
    SUBROUTINE TL_Kernel_CO2( k, ic, np )
      INTEGER, INTENT(IN) :: k, ic, np
      Predictor_TL%X(k, 1, ic)  = SECANG(k) * T_TL
      Predictor_TL%X(k, 2, ic)  = SECANG(k) * T2_TL
      Predictor_TL%X(k, 3, ic)  = T_TL
      Predictor_TL%X(k, 4, ic)  = T2_TL
      Predictor_TL%X(k, 5, ic)  = ZERO
      Predictor_TL%X(k, 6, ic)  = SECANG(k)*CO2_TL
      Predictor_TL%X(k, 7, ic)  = SECANG(k)*Tzp_TL(k)
      Predictor_TL%X(k, 8, ic)  = TWO*SECANG(k)**2 * PAFV%GAzp(k, ja_co2)* GAzp_TL(k, ja_co2)
      Predictor_TL%X(k, 9, ic)  = THREE*PAFV%Tzp(k)**2*Tzp_TL(k)
      Predictor_TL%X(k, 10, ic) = SECANG(k)*( SQRT(T)*Tzp_TL(k) + (POINT_5*PAFV%Tzp(k)/SQRT(T))*T_TL )
      IF ( np >= 11 ) THEN
        Predictor_TL%X(k, 11, ic) = CO_A_TL
      END IF
    END SUBROUTINE TL_Kernel_CO2

    !  --------------------------
    !  Water-line predictors; predictors 16 - 18 (CH4/CO cross terms) are
    !  carried only by rosters that request 18 predictors for WLO (group 1)
    !  --------------------------
    SUBROUTINE TL_Kernel_WLO( k, ic, np )
      INTEGER, INTENT(IN) :: k, ic, np
      Predictor_TL%X(k, 1, ic) = H2O_A_TL
      Predictor_TL%X(k, 2, ic) = H2O_A_TL*DT + H2O_A*DT_TL
      Predictor_TL%X(k, 3, ic) = H2O_S_TL
      Predictor_TL%X(k, 4, ic) = H2O_A_TL*DT2 + H2O_A*DT2_TL
      Predictor_TL%X(k, 5, ic) = H2O_R4_TL
      Predictor_TL%X(k, 6, ic) = H2O_S_TL*H2O_A + H2O_S*H2O_A_TL
      Predictor_TL%X(k, 7, ic) = H2O_R_TL
      Predictor_TL%X(k, 8, ic) = H2O_R_TL*DT + H2O_R*DT_TL
      Predictor_TL%X(k, 9, ic) = TWO*H2O_S*H2O_S_TL
      Predictor_TL%X(k,10, ic) = H2OdH2OTzp_TL
      Predictor_TL%X(k,11, ic) = H2O_R_TL*H2OdH2OTzp + H2O_R*H2OdH2OTzp_TL
      Predictor_TL%X(k,12, ic) = TWO*SECANG(k)**2 * PAFV%GAzp(k,ja_h2o)*GAzp_TL(k,ja_h2o)
      Predictor_TL%X(k,13, ic) = SECANG(k)*GAzp_TL(k,ja_h2o)
      Predictor_TL%X(k,14, ic) = ZERO
      Predictor_TL%X(k,15, ic) = SECANG(k)*CO2_TL
      IF ( np >= 18 ) THEN
        Predictor_TL%X(k,16, ic) = CH4_A_TL
        Predictor_TL%X(k,17, ic) = TWO*CH4_A*CH4_A_TL*DT + CH4_A*CH4_A*DT_TL
        Predictor_TL%X(k,18, ic) = CO_A_TL
      END IF
    END SUBROUTINE TL_Kernel_WLO

    !  -----------------------
    !  Carbon monoxide
    !  -----------------------
    SUBROUTINE TL_Kernel_CO( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      Predictor_TL%X(k, 1,  ic)   = CO_A_TL
      Predictor_TL%X(k, 2,  ic)   = CO_A_TL*DT + CO_A*DT_TL
      Predictor_TL%X(k, 3,  ic)   = (POINT_5/SQRT(CO_R))*CO_R_TL
      Predictor_TL%X(k, 4,  ic)   = CO_R_TL*DT + CO_R*DT_TL
      Predictor_TL%X(k, 5,  ic)   = CO_S_TL
      Predictor_TL%X(k, 6,  ic)   = CO_R_TL
      Predictor_TL%X(k, 7,  ic)   = CO_A_TL*DT2 + CO_A*DT2_TL
      Predictor_TL%X(k, 8,  ic)   = CO_ACOdCOzp_TL
      Predictor_TL%X(k, 9,  ic)   = CO_ACOdCOzp_TL/CO_R - CO_ACOdCOzp*CO_R_TL/CO_R**2
      Predictor_TL%X(k, 10, ic)   = CO_ACOdCOzp_TL * SQRT( PAFV%GAzp(k, ja_co) ) + &
                                    (POINT_5*CO_ACOdCOzp/SQRT(PAFV%GAzp(k, ja_co))) * &
                                    GAzp_TL(k, ja_co)
    END SUBROUTINE TL_Kernel_CO

    !  -----------------------
    !  Methane predictors
    !  -----------------------
    SUBROUTINE TL_Kernel_CH4( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      Predictor_TL%X(k, 1,  ic)  = CH4_A_TL*DT + CH4_A*DT_TL
      Predictor_TL%X(k, 2,  ic)  = CH4_R_TL
      Predictor_TL%X(k, 3,  ic)  = TWO*CH4_A*CH4_A_TL
      Predictor_TL%X(k, 4,  ic)  = CH4_A_TL
      Predictor_TL%X(k, 5,  ic)  = CH4_TL*DT + CH4*DT_TL
      Predictor_TL%X(k, 6,  ic)  = CH4_ACH4zp_TL
      Predictor_TL%X(k, 7,  ic)  = TWO*CH4_ACH4zp*CH4_ACH4zp_TL
      Predictor_TL%X(k, 8,  ic)  = (POINT_5/SQRT(CH4_R))*CH4_R_TL
      Predictor_TL%X(k, 9,  ic)  = GATzp_TL(k, ja_ch4)
      Predictor_TL%X(k, 10, ic)  = SECANG(k)*GATzp_TL(k, ja_ch4)
      Predictor_TL%X(k, 11, ic)  = CH4_R_TL*CH4/PAFV%GAzp(k, ja_ch4) + &
                                   CH4_R*CH4_TL/PAFV%GAzp(k, ja_ch4) - &
                                   CH4_R*CH4*GAzp_TL(k, ja_ch4)/PAFV%GAzp(k, ja_ch4)**2
    END SUBROUTINE TL_Kernel_CH4

    !  -----------------------
    !    N2O predictors
    !  -----------------------
    SUBROUTINE TL_Kernel_N2O( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      Predictor_TL%X(k, 1, ic)   = N2O_A_TL*DT + N2O_A*DT_TL
      Predictor_TL%X(k, 2, ic)   = N2O_R_TL
      Predictor_TL%X(k, 3, ic)   = N2O_TL*DT + N2O*DT_TL
      Predictor_TL%X(k, 4, ic)   = POINT_25*N2O_A**(-POINT_75) * N2O_A_TL
      Predictor_TL%X(k, 5, ic)   = N2O_A_TL
      Predictor_TL%X(k, 6, ic)   = SECANG(k) * GAzp_TL(k, ja_n2o)
      Predictor_TL%X(k, 7, ic)   = SECANG(k) * GATzp_TL(k, ja_n2o)
      Predictor_TL%X(k, 8, ic)   = N2O_S_TL
      Predictor_TL%X(k, 9, ic)   = GATzp_TL(k, ja_n2o)
      Predictor_TL%X(k,10, ic)   = N2O_R_TL*N2O / PAFV%GAzp(k, ja_n2o) + &
                                   N2O_R*N2O_TL / PAFV%GAzp(k, ja_n2o) - &
                                   N2O_R*N2O*GAzp_TL(k, ja_n2o)/PAFV%GAzp(k, ja_n2o)**2
      Predictor_TL%X(k,11, ic)   = CH4_A_TL
      Predictor_TL%X(k,12, ic)   = CH4_A_TL*PAFV%GAzp(k, ja_ch4) + CH4_A*GAzp_TL(k, ja_ch4)
      Predictor_TL%X(k,13, ic)   = CO_A_TL
      Predictor_TL%X(k,14, ic)   = CO_A_TL*SECANG(k)*PAFV%GAzp(k, ja_co) + &
                                   CO_A*SECANG(k)*GAzp_TL(k, ja_co)
    END SUBROUTINE TL_Kernel_N2O

    !  -----------------------
    !  NO2 predictors TL (GROUP_UV_NO2)
    !  -----------------------
    SUBROUTINE TL_Kernel_NO2( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      NO2      = Absorber(k,ja_no2)/Ref_Absorber(k, ja_no2)
      NO2_TL   = Absorber_TL(k,ja_no2)/Ref_Absorber(k, ja_no2)
      NO2_A    = SECANG(k)*NO2
      NO2_A_TL = SECANG(k)*NO2_TL
      Predictor_TL%X(k, 1, ic) = NO2_A_TL
      Predictor_TL%X(k, 2, ic) = NO2_A_TL*DT + NO2_A*DT_TL
      Predictor_TL%X(k, 3, ic) = NO2_A_TL*DT2 + NO2_A*DT2_TL
    END SUBROUTINE TL_Kernel_NO2

    !  --------------------------------
    !  Water vapor, MW (line and continuum together)
    !  --------------------------------
    SUBROUTINE TL_Kernel_WET_MW( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      Predictor_TL%X(k, 1, ic) = H2O_A_TL/T - H2O_A*T_TL/T**2
      Predictor_TL%X(k, 2, ic) = H2O_A_TL*H2O/T + H2O_A*H2O_TL/T - H2O_A*H2O*T_TL/T**2
      Predictor_TL%X(k, 3, ic) = H2O_A_TL*H2O/T2**2 + H2O_A*H2O_TL/T2**2 - &
                                 Two*H2O_A*H2O*T2_TL/T2**3
      Predictor_TL%X(k, 4, ic) = H2O_A_TL/T2 - H2O_A * T2_TL/T2**2
      Predictor_TL%X(k, 5, ic) = H2O_A_TL*H2O/T2 + H2O_A*H2O_TL/T2 - H2O_A*H2O*T2_TL/T2**2
      Predictor_TL%X(k, 6, ic) = H2O_A_TL/T2**2 - TWO*H2O_A*T2_TL/T2**3
      Predictor_TL%X(k, 7, ic) = H2O_A_TL
      Predictor_TL%X(k, 8, ic) = H2O_A_TL*DT + H2O_A*DT_TL
      Predictor_TL%X(k, 9, ic) = TWO*SECANG(k)**2 * PAFV%GAzp(k,ja_h2o)*GAzp_TL(k,ja_h2o)
      Predictor_TL%X(k, 10,ic) = SECANG(k)*GAzp_TL(k,ja_h2o)
      Predictor_TL%X(k, 11,ic) = ZERO
      Predictor_TL%X(k, 12,ic) = H2O_S_TL*H2O_A + H2O_S*H2O_A_TL
      Predictor_TL%X(k, 13,ic) = TWO*H2O_S*H2O_S_TL
      Predictor_TL%X(k, 14,ic) = H2OdH2OTzp_TL
    END SUBROUTINE TL_Kernel_WET_MW

  END SUBROUTINE  ODPS_Compute_Predictor_TL

!============================== END OF TL

!============================== AD

!------------------------------------------------------------------------------
!
! NAME:
!       ODPS_Compute_Predictor_AD
!
! PURPOSE:
!       Subroutine to predictors
!
! CALLING SEQUENCE:
!       CALL ODPS_Compute_Predictor_AD( &
!         Group_ID,           &  ! Input
!         Temperature,        &  ! Input
!         Absorber,           &  ! Input
!         Ref_Temperature,    &  ! Input
!         Ref_Absorber,       &  ! Input
!         secang,             &  ! Input
!         Predictor,          &  ! Input
!         Predictor_AD,       &  ! Input
!         Temperature_AD      &  ! Output
!         Absorber_AD         )  ! Output
!
! INPUT ARGUMENTS:
!       Group_ID   :     The ID of predictor group
!                        UNITS:      N?A
!                        TYPE:       INTEGER
!                        DIMENSION:  scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Temperature:     Temperature profile
!                        UNITS:      K
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!       Absorber   :     Absorber profiles
!                        UNITS:      vary
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-2(n_Layers x n_Absorbers) array
!                        ATTRIBUTES: INTENT(IN)
!
!       Ref_Temperature : Reference layer temperature profile
!                        UNITS:      K
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!       Ref_Absorber :   Reference absorber profiles
!                        UNITS:      vary
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-2(n_Layers x n_Absorbers) array
!                        ATTRIBUTES: INTENT(IN)
!
!       secang       :   Secont sensor zenith angle profile
!                        UNITS:      N/A
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!       Predictor:      Predictor structure containing the integrated absorber
!                       and predictor profiles.
!                       UNITS:      N/A
!                       TYPE:       ODPS_Predictor_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Predictor_AD:   Predictor_AD structure containing the integrated absorber_AD
!                       and predictor_AD profiles.
!                       UNITS:      N/A
!                       TYPE:       ODPS_Predictor_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!
!       Temperature_AD:  Temperature_AD profile
!                        UNITS:      K
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN OUT)
!
!       Absorber_AD:     Absorber_AD profiles
!                        UNITS:      vary
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-2(n_Layers x n_Absorbers) array
!                        ATTRIBUTES: INTENT(IN OUT)
!
!------------------------------------------------------------------------------

  SUBROUTINE ODPS_Compute_Predictor_AD( &
    Group_ID,           &
    Component_ID,       &
    Absorber_ID,        &
    Temperature,        &
    Absorber,           &
    Ref_Temperature,    &
    Ref_Absorber,       &
    secang,             &
    Predictor,          &
    Predictor_AD,       &
    Temperature_AD,     &
    Absorber_AD )

    INTEGER,                           INTENT(IN)     :: Group_ID
    INTEGER,                           INTENT(IN)     :: Component_ID(:)
    INTEGER,                           INTENT(IN)     :: Absorber_ID(:)
    REAL(fp),                          INTENT(IN)     :: Temperature(:)
    REAL(fp),                          INTENT(IN)     :: Absorber(:, :)
    REAL(fp),                          INTENT(IN)     :: Ref_Temperature(:)
    REAL(fp),                          INTENT(IN)     :: Ref_Absorber(:, :)
    REAL(fp),                          INTENT(IN)     :: Secang(:)
    TYPE(ODPS_Predictor_type), TARGET, INTENT(IN)     :: Predictor
    TYPE(ODPS_Predictor_type),         INTENT(IN OUT) :: Predictor_AD
    REAL(fp),                          INTENT(IN OUT) :: Temperature_AD(:)
    REAL(fp),                          INTENT(IN OUT) :: Absorber_AD(:, :)

    ! ---------------
    ! Local variables
    ! ---------------
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'ODPS_Compute_Predictor_AD'
    INTEGER    ::    n_layers
    INTEGER    ::    k ! n_Layers, n_Levels
    INTEGER    ::    j ! n_absorbers
    REAL(fp) ::    Tzp_sum_AD
    REAL(fp) ::    Tzp_AD(SIZE(Absorber, DIM=1))
    REAL(fp) ::    Tz_sum_AD
    REAL(fp) ::    Tz_AD(SIZE(Absorber, DIM=1))
    REAL(fp) ::    GAz_sum_AD(SIZE(Absorber, DIM=2))
    REAL(fp) ::    GAz_AD(SIZE(Absorber, DIM=1), SIZE(Absorber, DIM=2))
    REAL(fp) ::    GAzp_sum_AD(SIZE(Absorber, DIM=2))
    REAL(fp) ::    GAzp_AD(SIZE(Absorber, DIM=1), SIZE(Absorber, DIM=2))
    REAL(fp) ::    GATzp_sum_AD(SIZE(Absorber, DIM=2))
    REAL(fp) ::    GATzp_AD(SIZE(Absorber, DIM=1), SIZE(Absorber, DIM=2))
    ! Kernel dispatch bookkeeping and shared per-layer variables
    INTEGER  :: ic
    INTEGER  :: ic_dry, ic_wlo, ic_wco, ic_ozo, ic_co2
    INTEGER  :: ic_n2o, ic_co, ic_ch4, ic_no2, ic_wet
    INTEGER  :: ja_h2o, ja_o3, ja_co2, ja_n2o, ja_co, ja_ch4, ja_no2
    LOGICAL  :: has_trace
    REAL(fp) :: DT, DT_AD, T, T_AD, T2, T2_AD, DT2, DT2_AD
    REAL(fp) :: H2O, H2O_AD, H2O_A, H2O_A_AD, H2O_R, H2O_R_AD
    REAL(fp) :: H2O_S, H2O_S_AD, H2O_R4, H2O_R4_AD, H2OdH2OTzp, H2OdH2OTzp_AD
    REAL(fp) :: CO2, CO2_AD, O3, O3_AD, O3_A, O3_A_AD, O3_R, O3_R_AD
    REAL(fp) :: NO2, NO2_AD, NO2_A, NO2_A_AD
    REAL(fp) :: CO, CO_AD, CO_A, CO_A_AD, CO_R, CO_R_AD, CO_S, CO_S_AD
    REAL(fp) :: CO_ACOdCOzp, CO_ACOdCOzp_AD
    REAL(fp) :: N2O, N2O_AD, N2O_A, N2O_A_AD, N2O_R, N2O_R_AD, N2O_S, N2O_S_AD
    REAL(fp) :: CH4, CH4_AD, CH4_A, CH4_A_AD, CH4_R, CH4_R_AD
    REAL(fp) :: CH4_ACH4zp, CH4_ACH4zp_AD
!JR Static initialization means only 1 copy of the variable. OpenMP over profiles
!JR means $OPENMP_NUM_THREADS copies are needed. So change to run-time initialization
!    TYPE(PAFV_type), POINTER  :: PAFV => NULL()
    TYPE(PAFV_type), POINTER  :: PAFV

!JR Looks silly since next stmt sets PAFV, but nullify here to guard agains future code changes
    NULLIFY(PAFV)
    ! use short name
    PAFV => Predictor%PAFV

    n_Layers = Predictor%n_Layers


    !------------------------------------
    ! Compute integrated variables
    !------------------------------------

    Tzp_sum_AD   = ZERO
    Tz_sum_AD    = ZERO
    GAz_sum_AD   = ZERO
    GAzp_sum_AD  = ZERO
    GATzp_sum_AD = ZERO

    GATzp_AD = ZERO
    GAzp_AD  = ZERO
    GAz_AD   = ZERO
    Tzp_AD   = ZERO
    Tz_AD    = ZERO

    !----------------------------------------------------------------
    ! Per-component adjoint predictor computation. Unlike the forward
    ! and tangent-linear cases, the ORDER of the adjoint kernels is
    ! bit-significant: adjoint contributions accumulate into shared
    ! variables, and floating-point sums depend on evaluation order.
    ! The kernel call sequence below reproduces the heritage monolith
    ! order exactly: DRY, WCO, NO2, OZO, CO2, WLO (with the group-1
    ! extension), CO, CH4, N2O, then the trace-gas variable chains,
    ! then the common variable chains. Do not reorder.
    !----------------------------------------------------------------

    ! Resolve each gas's position in this group's absorber roster
    ja_h2o = Absorber_Position(H2O_ID)
    ja_o3  = Absorber_Position(O3_ID)
    ja_co2 = Absorber_Position(CO2_ID)
    ja_n2o = Absorber_Position(N2O_ID)
    ja_co  = Absorber_Position(CO_ID)
    ja_ch4 = Absorber_Position(CH4_ID)
    ja_no2 = Absorber_Position(NO2_ID)
    has_trace = ANY( Component_ID == CO_ComID )   ! validation guarantees the trio

    ! Resolve each component's position in this group's roster
    ic_dry = Component_Position(DRY_ComID_G1)
    IF ( ic_dry == 0 ) ic_dry = Component_Position(DRY_ComID_G2)
    IF ( ic_dry == 0 ) ic_dry = Component_Position(DRY_ComID_Q1)
    IF ( ic_dry == 0 ) ic_dry = Component_Position(DRY_ComID_Q2)
    IF ( ic_dry == 0 ) ic_dry = Component_Position(EDRY_ComID)
    ic_wlo = Component_Position(WLO_ComID)
    ic_wco = Component_Position(WCO_ComID)
    ic_ozo = Component_Position(OZO_ComID)
    ic_co2 = Component_Position(CO2_ComID)
    ic_n2o = Component_Position(N2O_ComID)
    ic_co  = Component_Position(CO_ComID)
    ic_ch4 = Component_Position(CH4_ComID)
    ic_no2 = Component_Position(NO2_ComID)
    ic_wet = Component_Position(WET_ComID)

    ! Zero the adjoint accumulators (once per call, as in the heritage code)
    DT_AD          = ZERO
    T_AD           = ZERO
    T2_AD          = ZERO
    DT2_AD         = ZERO
    H2O_AD         = ZERO
    H2O_A_AD       = ZERO
    H2O_R_AD       = ZERO
    H2O_S_AD       = ZERO
    H2O_R4_AD      = ZERO
    H2OdH2OTzp_AD  = ZERO
    CO2_AD         = ZERO
    O3_AD          = ZERO
    O3_A_AD        = ZERO
    O3_R_AD        = ZERO
    NO2_AD         = ZERO
    NO2_A_AD       = ZERO
    CO_AD          = ZERO
    CO_A_AD        = ZERO
    CO_R_AD        = ZERO
    CO_S_AD        = ZERO
    CO_ACOdCOzp_AD = ZERO
    N2O_AD         = ZERO
    N2O_A_AD       = ZERO
    N2O_R_AD       = ZERO
    N2O_S_AD       = ZERO
    CH4_AD         = ZERO
    CH4_A_AD       = ZERO
    CH4_R_AD       = ZERO
    CH4_ACH4zp_AD  = ZERO

    Basis_Select: IF ( GROUP_REGISTRY(Group_ID)%Basis == BASIS_IR ) THEN

      IR_Layer_Loop : DO k = n_Layers, 1, -1

        !------------------------------------------
        !  Relative Temperature
        !------------------------------------------
        dT = Temperature(k) - Ref_Temperature(k)
        T = Temperature(k) / Ref_Temperature(k)

        !-------------------------------------------
        !  Abosrber amount scalled by the reference
        !-------------------------------------------
        IF ( ja_h2o > 0 ) H2O = Absorber(k,ja_h2o)/Ref_Absorber(k, ja_h2o)
        IF ( ja_o3  > 0 ) O3  = Absorber(k,ja_o3)/Ref_absorber(k,ja_o3)
        IF ( ja_co2 > 0 ) CO2 = Absorber(k,ja_co2)/Ref_absorber(k,ja_co2)

        ! Combinations of variables common to all predictor groups
        T2   = T*T
        DT2  = DT*ABS( DT )

        IF ( ja_h2o > 0 ) THEN
          H2O_A = SECANG(k)*H2O
          H2O_R  = SQRT( H2O_A )
          H2O_S  = H2O_A*H2O_A
          H2O_R4 = SQRT( H2O_R )
          H2OdH2OTzp = H2O/PAFV%GATzp(k, ja_h2o)
        END IF

        IF ( ja_o3 > 0 ) THEN
          O3_A = SECANG(k)*O3
          O3_R = SQRT( O3_A )
        END IF

        IF( has_trace )THEN
          CO  = Absorber(k,ja_co)/Ref_absorber(k, ja_co)
          N2O = Absorber(k,ja_n2o)/Ref_absorber(k,ja_n2o)
          CH4 = Absorber(k,ja_ch4)/Ref_absorber(k,ja_ch4)

          N2O_A = SECANG(k)*N2O
          N2O_R = SQRT( N2O_A )
          N2O_S = N2O_A*N2O_A

          CO_A = SECANG(k)*CO
          CO_R = SQRT( CO_A )
          CO_S = CO_A*CO_A
          CO_ACOdCOzp = CO_A*CO/PAFV%GAzp(k, ja_co)

          CH4_A = SECANG(k)*CH4
          CH4_R = SQRT(CH4_A)
          CH4_ACH4zp = SECANG(k)*PAFV%GAzp(k, ja_ch4)

        END IF

        ! Adjoint kernels in the heritage monolith order (see note above)
        IF ( ic_dry > 0 ) CALL AD_Kernel_DRY(k, ic_dry, ODPS_Kernel_n_Predictors( &
          GROUP_REGISTRY(Group_ID)%Basis, Component_ID(ic_dry), has_trace ))
        IF ( ic_wco > 0 ) CALL AD_Kernel_WCO(k, ic_wco)
        IF ( ic_no2 > 0 ) CALL AD_Kernel_NO2(k, ic_no2)
        IF ( ic_ozo > 0 ) CALL AD_Kernel_OZO_IR(k, ic_ozo)
        IF ( ic_co2 > 0 ) CALL AD_Kernel_CO2(k, ic_co2)
        IF ( ic_wlo > 0 ) CALL AD_Kernel_WLO(k, ic_wlo, ODPS_Kernel_n_Predictors( &
          GROUP_REGISTRY(Group_ID)%Basis, WLO_ComID, has_trace ))
        IF ( has_trace ) THEN
          CALL AD_Kernel_CO(k, ic_co)
          CALL AD_Kernel_CH4(k, ic_ch4)
          CALL AD_Kernel_N2O(k, ic_n2o)

          ! Trace-gas variable chains (group 1)
          GAzp_AD(k, ja_ch4) = GAzp_AD(k, ja_ch4) + SECANG(k)*CH4_ACH4zp_AD
          CH4_A_AD = CH4_A_AD + (POINT_5/SQRT(CH4_A)) * CH4_R_AD
          CH4_AD = CH4_AD + SECANG(k)*CH4_A_AD
          CH4_ACH4zp_AD = ZERO
          CH4_R_AD =      ZERO
          CH4_A_AD =      ZERO

          GAzp_AD(k, ja_co) = GAzp_AD(k, ja_co)                                      &
                              - CO_ACOdCOzp_AD*CO_A*CO/PAFV%GAzp(k, ja_co)**2
          CO_A_AD = CO_A_AD + CO_ACOdCOzp_AD*CO/PAFV%GAzp(k, ja_co)                  &
                            + CO_S_AD*TWO * CO_A                                     &
                            + CO_R_AD*POINT_5/SQRT(CO_A)
          CO_AD = CO_AD + CO_ACOdCOzp_AD*CO_A/PAFV%GAzp(k, ja_co)                    &
                        + CO_A_AD*SECANG(k)
          CO_ACOdCOzp_AD = ZERO
          CO_S_AD        = ZERO
          CO_R_AD        = ZERO
          CO_A_AD        = ZERO

          N2O_A_AD = N2O_A_AD + N2O_R_AD * POINT_5 / SQRT(N2O_A)                    &
                              + N2O_S_AD * TWO * N2O_A
          N2O_AD = N2O_AD + N2O_A_AD * SECANG(k)
          N2O_A_AD = ZERO
          N2O_R_AD = ZERO
          N2O_S_AD = ZERO

          Absorber_AD(k,ja_ch4) = Absorber_AD(k,ja_ch4)                              &
                                    + CH4_AD/Ref_absorber(k,ja_ch4)
          Absorber_AD(k,ja_n2o) = Absorber_AD(k,ja_n2o)                              &
                                    + N2O_AD/Ref_absorber(k,ja_n2o)
          Absorber_AD(k,ja_co)  = Absorber_AD(k,ja_co)                               &
                                    + CO_AD/Ref_absorber(k, ja_co)
          CO_AD  = ZERO
          N2O_AD = ZERO
          CH4_AD = ZERO

        END IF

        ! Combinations of variables common to all predictor groups

        IF ( ja_o3 > 0 ) THEN
          O3_A_AD = O3_A_AD + O3_R_AD * POINT_5 / SQRT(O3_A)
          O3_AD = O3_AD + O3_A_AD * SECANG(k)
          O3_A_AD = ZERO
          O3_R_AD = ZERO
        END IF

        IF ( ja_h2o > 0 ) THEN
          GATzp_AD(k, ja_h2o) = GATzp_AD(k, ja_h2o)                                  &
                                  - H2OdH2OTzp_AD*H2O/PAFV%GATzp(k, ja_h2o)**2
          H2O_R_AD = H2O_R_AD + H2O_R4_AD * POINT_5 / SQRT(H2O_R)
          H2O_A_AD = H2O_A_AD + H2O_S_AD * TWO * H2O_A                               &
                   + H2O_R_AD * POINT_5 / SQRT(H2O_A)
          H2O_AD = H2O_AD + H2O_A_AD * SECANG(k)                                     &
                 + H2OdH2OTzp_AD / PAFV%GATzp(k, ja_h2o)

          H2O_A_AD      = ZERO
          H2O_R_AD      = ZERO
          H2O_S_AD      = ZERO
          H2O_R4_AD     = ZERO
          H2OdH2OTzp_AD = ZERO
        END IF

        IF( DT > ZERO) THEN
          DT_AD = DT_AD + DT2_AD*TWO*DT
        ELSE
          DT_AD = DT_AD - DT2_AD*TWO*DT
        ENDIF
        T_AD = T_AD + T2_AD*TWO*T
        T2_AD   = ZERO
        DT2_AD  = ZERO

        !-------------------------------------------
        !  Abosrber amount scalled by the reference
        !-------------------------------------------
        IF ( ja_co2 > 0 ) THEN
          Absorber_AD(k,ja_co2) = Absorber_AD(k,ja_co2)            &
                                    + CO2_AD / Ref_absorber(k,ja_co2)
        END IF
        IF ( ja_o3 > 0 ) THEN
          Absorber_AD(k,ja_o3)  = Absorber_AD(k,ja_o3)             &
                                    + O3_AD / Ref_absorber(k,ja_o3)
        END IF
        IF ( ja_h2o > 0 ) THEN
          Absorber_AD(k,ja_h2o) = Absorber_AD(k,ja_h2o)            &
                                    + H2O_AD / Ref_Absorber(k, ja_h2o)
        END IF
        H2O_AD = ZERO
        O3_AD  = ZERO
        CO2_AD = ZERO

        !------------------------------------------
        !  Relative Temperature
        !------------------------------------------
        Temperature_AD(k) = Temperature_AD(k) + T_AD/Ref_Temperature(k) &
                          + dT_AD
        dT_AD = ZERO
        T_AD  = ZERO

      END DO IR_Layer_Loop

    ELSE Basis_Select   ! BASIS_MW

      ! set number of predictors
      DO ic = 1, SIZE(Component_ID)
        Predictor_AD%n_CP(ic) = ODPS_Kernel_n_Predictors( &
          GROUP_REGISTRY(Group_ID)%Basis, Component_ID(ic), has_trace )
      END DO

      MW_Layer_Loop : DO k = n_Layers, 1, -1

        !------------------------------------------
        !  Relative Temperature
        !------------------------------------------
        dT = Temperature(k) - Ref_Temperature(k)
        T = Temperature(k) / Ref_Temperature(k)

        !-------------------------------------------
        !  Abosrber amount scalled by the reference
        !-------------------------------------------
        ! Combinations of variables common to all predictor groups
        T2  = T*T
        DT2 = DT*ABS( DT )

        IF ( ja_h2o > 0 ) THEN
          H2O = Absorber(k,ja_h2o)/Ref_Absorber(k, ja_h2o)
          H2O_A = SECANG(k)*H2O
          H2O_R  = SQRT( H2O_A )
          H2O_S  = H2O_A*H2O_A
          H2O_R4 = SQRT( H2O_R )
          H2OdH2OTzp = H2O/PAFV%GATzp(k, ja_h2o)
        END IF

        IF ( ja_o3 > 0 ) THEN
          O3   = Absorber(k,ja_o3)/Ref_Absorber(k, ja_o3)
          O3_A = SECANG(k)*O3
          O3_R = SQRT( O3_A )
        END IF

        ! Adjoint kernels in the heritage monolith order (see note above)
        IF ( ic_dry > 0 ) CALL AD_Kernel_DRY(k, ic_dry, ODPS_Kernel_n_Predictors( &
          GROUP_REGISTRY(Group_ID)%Basis, Component_ID(ic_dry), has_trace ))
        IF ( ic_wet > 0 ) CALL AD_Kernel_WET_MW(k, ic_wet)
        IF ( ic_ozo > 0 ) CALL AD_Kernel_OZO_MW(k, ic_ozo)

        !-------------------------------------------
        !  Abosrber amount scalled by the reference
        !-------------------------------------------

        ! Combinations of variables common to all predictor groups

        IF ( ja_h2o > 0 ) THEN
          GATzp_AD(k, ja_h2o) = GATzp_AD(k, ja_h2o)                                  &
                                  - H2OdH2OTzp_AD*H2O/PAFV%GATzp(k, ja_h2o)**2
          H2O_R_AD = H2O_R_AD + H2O_R4_AD * POINT_5 / SQRT(H2O_R)
          H2O_A_AD = H2O_A_AD + H2O_S_AD * TWO * H2O_A                               &
                   + H2O_R_AD * POINT_5 / SQRT(H2O_A)
          H2O_AD = H2O_AD + H2O_A_AD * SECANG(k)                                     &
                 + H2OdH2OTzp_AD / PAFV%GATzp(k, ja_h2o)
          H2O_A_AD      = ZERO
          H2O_R_AD      = ZERO
          H2O_S_AD      = ZERO
          H2O_R4_AD     = ZERO
          H2OdH2OTzp_AD = ZERO
        END IF

        IF( DT > ZERO) THEN
          DT_AD = DT_AD + DT2_AD*TWO*DT
        ELSE
          DT_AD = DT_AD - DT2_AD*TWO*DT
        ENDIF
        T_AD = T_AD + T2_AD*TWO*T
        T2_AD   = ZERO
        DT2_AD  = ZERO

        IF ( ja_h2o > 0 ) THEN
          Absorber_AD(k,ja_h2o) = Absorber_AD(k,ja_h2o)            &
                                    + H2O_AD / Ref_Absorber(k, ja_h2o)
        END IF
        H2O_AD = ZERO

        !------------------------------------------
        !  Relative Temperature
        !------------------------------------------
        Temperature_AD(k) = Temperature_AD(k) + T_AD/Ref_Temperature(k) &
                          + dT_AD
        dT_AD = ZERO
        T_AD  = ZERO

      END DO MW_Layer_Loop

    END IF Basis_Select

    Adjoint_Layer_Loop : DO k = n_Layers, 1, -1

      ! absorbers
      DO j = SIZE(Absorber, DIM=2), 1, -1

        GATzp_sum_AD(j) = GATzp_sum_AD(j) + GATzp_AD(k, j)/PAFV%GATzp_ref(k,j)
        GAzp_sum_AD(j) = GAzp_sum_AD(j) + GAzp_AD(k, j)/PAFV%GAzp_ref(k,j)
        GAz_sum_AD(j) = GAz_sum_AD(j) + GAz_AD(k, j)/PAFV%GAz_ref(k,j)
        Temperature_AD(k) = Temperature_AD(k) + GATzp_sum_AD(j)*PAFV%PDP(k)*Absorber(k, j)
        Absorber_AD(k, j) = Absorber_AD(k, j) + GAz_sum_AD(j) + GAzp_sum_AD(j)*PAFV%PDP(k) &
                                              + GATzp_sum_AD(j)*PAFV%PDP(k)*Temperature(k)
        GATzp_AD(k, j) = ZERO
        GAzp_AD(k, j)  = ZERO
        GAz_AD(k, j)   = ZERO

      END DO

      ! Temperature
      Tzp_sum_AD = Tzp_sum_AD + Tzp_AD(k)/PAFV%Tzp_ref(k)
      Tz_sum_AD = Tz_sum_AD + Tz_AD(k)/PAFV%Tz_ref(k)
      Temperature_AD(k) = Temperature_AD(k) + Tz_sum_AD + PAFV%PDP(k)*Tzp_sum_AD
      Tzp_AD(k) = ZERO
      Tz_AD(k)  = ZERO

    END DO Adjoint_Layer_Loop

    NULLIFY(PAFV)

CONTAINS

    ! Position of a HITRAN absorber ID in the file's absorber roster
    PURE FUNCTION Absorber_Position( Gas_ID ) RESULT( Position )
      INTEGER, INTENT(IN) :: Gas_ID
      INTEGER :: Position
      INTEGER :: ja
      Position = 0
      DO ja = 1, SIZE(Absorber_ID)
        IF ( Absorber_ID(ja) == Gas_ID ) THEN
          Position = ja
          RETURN
        END IF
      END DO
    END FUNCTION Absorber_Position

    ! Position of a component ID in the file's component roster
    PURE FUNCTION Component_Position( Com_ID ) RESULT( Position )
      INTEGER, INTENT(IN) :: Com_ID
      INTEGER :: Position
      INTEGER :: jc
      Position = 0
      DO jc = 1, SIZE(Component_ID)
        IF ( Component_ID(jc) == Com_ID ) THEN
          Position = jc
          RETURN
        END IF
      END DO
    END FUNCTION Component_Position

    !  ----------------------
    !   Fixed (Dry) predictors (IR and MW use the same formulation)
    !  ----------------------
    SUBROUTINE AD_Kernel_DRY( k, ic, np )
      INTEGER, INTENT(IN) :: k, ic, np
      T_AD     = T_AD                                    &
                 + Predictor_AD%X(k, 2, ic) * SECANG(k)  &
                 + Predictor_AD%X(k, 4, ic)
      T2_AD    = T2_AD                                   &
                 + Predictor_AD%X(k, 3, ic) * SECANG(k)  &
                 + Predictor_AD%X(k, 6, ic)
      Tz_AD(k) = Tz_AD(k) + Predictor_AD%X(k, 7, ic)

      Predictor_AD%X(k, 1, ic)  = ZERO
      Predictor_AD%X(k, 2, ic)  = ZERO
      Predictor_AD%X(k, 3, ic)  = ZERO
      Predictor_AD%X(k, 4, ic)  = ZERO
      Predictor_AD%X(k, 5, ic)  = ZERO
      Predictor_AD%X(k, 6, ic)  = ZERO
      Predictor_AD%X(k, 7, ic)  = ZERO

      ! QDRY water predictors 8 - 20: adjoint of the WLO formulas 1 - 13.
      ! Accumulates into the same water adjoint variables the WLO/WCO
      ! kernels use; the common variable chains at the end of the layer
      ! iteration convert them to Absorber_AD / Temperature_AD.
      IF ( np >= 20 ) THEN
        H2O_A_AD = H2O_A_AD                                &
                   + Predictor_AD%X(k, 8, ic)              &
                   + Predictor_AD%X(k, 9, ic)*DT           &
                   + Predictor_AD%X(k,11, ic)*DT2          &
                   + Predictor_AD%X(k,13, ic)*H2O_S

        DT_AD    = DT_AD                                   &
                   + Predictor_AD%X(k, 9, ic)*H2O_A        &
                   + Predictor_AD%X(k,15, ic)*H2O_R

        H2O_S_AD = H2O_S_AD                                &
                   + Predictor_AD%X(k,10, ic)              &
                   + Predictor_AD%X(k,13, ic)*H2O_A        &
                   + Predictor_AD%X(k,16, ic)*TWO*H2O_S

        DT2_AD   = DT2_AD + Predictor_AD%X(k,11, ic)*H2O_A

        H2O_R4_AD= H2O_R4_AD + Predictor_AD%X(k,12, ic)

        H2O_R_AD = H2O_R_AD                                &
                   + Predictor_AD%X(k,14, ic)              &
                   + Predictor_AD%X(k,15, ic)*DT           &
                   + Predictor_AD%X(k,18, ic)*H2OdH2OTzp

        H2OdH2OTzp_AD = H2OdH2OTzp_AD                      &
                   + Predictor_AD%X(k,17, ic)              &
                   + Predictor_AD%X(k,18, ic)*H2O_R

        GAzp_AD(k,ja_h2o) = GAzp_AD(k,ja_h2o)                                       &
                   + Predictor_AD%X(k,19, ic)*TWO*SECANG(k)**2*PAFV%GAzp(k,ja_h2o)  &
                   + Predictor_AD%X(k,20, ic)*SECANG(k)

        Predictor_AD%X(k, 8:20, ic) = ZERO
      END IF
    END SUBROUTINE AD_Kernel_DRY

    !  --------------------------
    !  Water vapor continuum predictors
    !  --------------------------
    SUBROUTINE AD_Kernel_WCO( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      H2O_A_AD = H2O_A_AD                                &
                 + Predictor_AD%X(k, 1, ic)/T            &
                 + Predictor_AD%X(k, 2, ic)*H2O/T        &
                 + Predictor_AD%X(k, 3, ic)*H2O/T2**2    &
                 + Predictor_AD%X(k, 4, ic)/T2           &
                 + Predictor_AD%X(k, 5, ic)*H2O/T2       &
                 + Predictor_AD%X(k, 6, ic)/T2**2        &
                 + Predictor_AD%X(k, 7, ic)
      T_AD     = T_AD                                       &
                 - Predictor_AD%X(k, 1, ic)*H2O_A/T**2      &
                 - Predictor_AD%X(k, 2, ic)*H2O_A*H2O/T**2
      H2O_AD   = H2O_AD                                 &
                 + Predictor_AD%X(k, 2, ic)*H2O_A/T     &
                 + Predictor_AD%X(k, 3, ic)*H2O_A/T2**2 &
                 + Predictor_AD%X(k, 5, ic)*H2O_A/T2
      T2_AD    = T2_AD                                          &
                 - Predictor_AD%X(k, 3, ic)*TWO*H2O_A*H2O/T2**3 &
                 - Predictor_AD%X(k, 4, ic)*H2O_A/T2**2         &
                 - Predictor_AD%X(k, 5, ic)*H2O_A*H2O/T2**2     &
                 - Predictor_AD%X(k, 6, ic)*TWO*H2O_A/T2**3

      Predictor_AD%X(k, 1, ic) = ZERO
      Predictor_AD%X(k, 2, ic) = ZERO
      Predictor_AD%X(k, 3, ic) = ZERO
      Predictor_AD%X(k, 4, ic) = ZERO
      Predictor_AD%X(k, 5, ic) = ZERO
      Predictor_AD%X(k, 6, ic) = ZERO
      Predictor_AD%X(k, 7, ic) = ZERO
    END SUBROUTINE AD_Kernel_WCO

    !  -----------------------
    !  NO2 predictors AD (GROUP_UV_NO2); adjoint of the compact set
    !  X1=NO2_A, X2=NO2_A*DT, X3=NO2_A*DT2, NO2_A=secang*NO2, NO2=Abs/Ref.
    !  Self-contained through to Absorber_AD, as in the heritage code.
    !  -----------------------
    SUBROUTINE AD_Kernel_NO2( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      NO2   = Absorber(k,ja_no2)/Ref_Absorber(k, ja_no2)
      NO2_A = SECANG(k)*NO2

      NO2_A_AD = NO2_A_AD                            &
                 + Predictor_AD%X(k, 1, ic)          &
                 + Predictor_AD%X(k, 2, ic)*DT       &
                 + Predictor_AD%X(k, 3, ic)*DT2
      DT_AD    = DT_AD  + Predictor_AD%X(k, 2, ic)*NO2_A
      DT2_AD   = DT2_AD + Predictor_AD%X(k, 3, ic)*NO2_A
      Predictor_AD%X(k, 1, ic) = ZERO
      Predictor_AD%X(k, 2, ic) = ZERO
      Predictor_AD%X(k, 3, ic) = ZERO

      NO2_AD   = NO2_AD + NO2_A_AD*SECANG(k)
      NO2_A_AD = ZERO
      Absorber_AD(k,ja_no2) = Absorber_AD(k,ja_no2)     &
                                + NO2_AD/Ref_Absorber(k, ja_no2)
      NO2_AD   = ZERO
    END SUBROUTINE AD_Kernel_NO2

    !  -----------------------
    !  Ozone predictors (IR groups; chain propagation happens in the
    !  common epilogue)
    !  -----------------------
    SUBROUTINE AD_Kernel_OZO_IR( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      O3_A_AD = O3_A_AD                                                       &
                + Predictor_AD%X(k, 1, ic)                                    &
                + Predictor_AD%X(k, 2, ic)*DT                                 &
                + Predictor_AD%X(k, 3, ic)*O3*PAFV%GAzp(k,ja_o3)              &
                + Predictor_AD%X(k, 4, ic)*TWO*O3_A                           &
                + Predictor_AD%X(k, 5, ic)*PAFV%GAzp(k,ja_o3)                 &
                + Predictor_AD%X(k, 6, ic)*SQRT(SECANG(k)*PAFV%GAzp(k,ja_o3))

      DT_AD   = DT_AD                                                         &
                + Predictor_AD%X(k, 2, ic)*O3_A                               &
                + Predictor_AD%X(k, 7, ic)*O3_R

      O3_AD   = O3_AD                                                         &
                + Predictor_AD%X(k, 3, ic)*O3_A*PAFV%GAzp(k,ja_o3)            &
                + Predictor_AD%X(k, 9, ic)*O3_R/PAFV%GAzp(k,ja_o3)

      GAzp_AD(k,ja_o3) = GAzp_AD(k,ja_o3)                                     &
                + Predictor_AD%X(k, 3, ic)*O3_A*O3                            &
                + Predictor_AD%X(k, 5, ic)*O3_A                               &
                + Predictor_AD%X(k, 6, ic)*POINT_5*O3_A*SQRT(SECANG(k)/PAFV%GAzp(k,ja_o3)) &
                - Predictor_AD%X(k, 9, ic)*O3_R*O3/PAFV%GAzp(k,ja_o3)**2      &
                + Predictor_AD%X(k,10, ic)*SECANG(k)                          &
                + Predictor_AD%X(k,11, ic)*TWO*SECANG(k)**2*PAFV%GAzp(k,ja_o3)

      O3_R_AD = O3_R_AD                                                       &
                + Predictor_AD%X(k, 7, ic)*DT                                 &
                + Predictor_AD%X(k, 8, ic)                                    &
                + Predictor_AD%X(k, 9, ic)*O3/PAFV%GAzp(k,ja_o3)

      Predictor_AD%X(k, 1, ic)  = ZERO
      Predictor_AD%X(k, 2, ic)  = ZERO
      Predictor_AD%X(k, 3, ic)  = ZERO
      Predictor_AD%X(k, 4, ic)  = ZERO
      Predictor_AD%X(k, 5, ic)  = ZERO
      Predictor_AD%X(k, 6, ic)  = ZERO
      Predictor_AD%X(k, 7, ic)  = ZERO
      Predictor_AD%X(k, 8, ic)  = ZERO
      Predictor_AD%X(k, 9, ic)  = ZERO
      Predictor_AD%X(k,10, ic)  = ZERO
      Predictor_AD%X(k,11, ic)  = ZERO

      Predictor_AD%X(k,12, ic)  = ZERO
      Predictor_AD%X(k,13, ic)  = ZERO
    END SUBROUTINE AD_Kernel_OZO_IR

    !  -----------------------
    !  Ozone predictors (GROUP_MW_O3; self-contained through to
    !  Absorber_AD, as in the heritage code)
    !  -----------------------
    SUBROUTINE AD_Kernel_OZO_MW( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      O3_A_AD = O3_A_AD                                                  &
                + Predictor_AD%X(k, 1, ic)                               &
                + Predictor_AD%X(k, 2, ic)*DT                            &
                + Predictor_AD%X(k, 3, ic)*O3*PAFV%GAzp(k,ja_o3)         &
                + Predictor_AD%X(k, 4, ic)*TWO*O3_A                      &
                + Predictor_AD%X(k, 5, ic)*PAFV%GAzp(k,ja_o3)            &
                + Predictor_AD%X(k, 6, ic)*SQRT(SECANG(k)*PAFV%GAzp(k,ja_o3))

      DT_AD   = DT_AD                                                    &
                + Predictor_AD%X(k, 2, ic)*O3_A                          &
                + Predictor_AD%X(k, 7, ic)*O3_R

      O3_AD   = O3_AD                                                    &
                + Predictor_AD%X(k, 3, ic)*O3_A*PAFV%GAzp(k,ja_o3)       &
                + Predictor_AD%X(k, 9, ic)*O3_R/PAFV%GAzp(k,ja_o3)

      GAzp_AD(k,ja_o3) = GAzp_AD(k,ja_o3)                                &
                + Predictor_AD%X(k, 3, ic)*O3_A*O3                       &
                + Predictor_AD%X(k, 5, ic)*O3_A                          &
                + Predictor_AD%X(k, 6, ic)*POINT_5*O3_A*SQRT(SECANG(k)/PAFV%GAzp(k,ja_o3)) &
                - Predictor_AD%X(k, 9, ic)*O3_R*O3/PAFV%GAzp(k,ja_o3)**2 &
                + Predictor_AD%X(k,10, ic)*SECANG(k)                     &
                + Predictor_AD%X(k,11, ic)*TWO*SECANG(k)**2*PAFV%GAzp(k,ja_o3)

      O3_R_AD = O3_R_AD                                                  &
                + Predictor_AD%X(k, 7, ic)*DT                            &
                + Predictor_AD%X(k, 8, ic)                               &
                + Predictor_AD%X(k, 9, ic)*O3/PAFV%GAzp(k,ja_o3)

      Predictor_AD%X(k, 1:11, ic) = ZERO

      O3_A_AD = O3_A_AD + O3_R_AD * POINT_5 / SQRT(O3_A)
      O3_AD   = O3_AD + O3_A_AD * SECANG(k)
      O3_R_AD = ZERO
      O3_A_AD = ZERO
      Absorber_AD(k,ja_o3) = Absorber_AD(k,ja_o3)            &
                               + O3_AD / Ref_Absorber(k, ja_o3)
      O3_AD   = ZERO
    END SUBROUTINE AD_Kernel_OZO_MW

    !  -----------------------
    !  Carbon dioxide predictors (1 - 10; the group-1 predictor 11 adjoint
    !  is handled with the WLO extension to preserve the heritage
    !  accumulation order)
    !  -----------------------
    SUBROUTINE AD_Kernel_CO2( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      T_AD      = T_AD                                          &
                  + Predictor_AD%X(k, 1, ic)*SECANG(k)          &
                  + Predictor_AD%X(k, 3, ic)                    &
                  + Predictor_AD%X(k, 10, ic)*SECANG(k)*(POINT_5*PAFV%Tzp(k)/SQRT(T))

      T2_AD     = T2_AD                                         &
                  + Predictor_AD%X(k, 2, ic)*SECANG(k)          &
                  + Predictor_AD%X(k, 4, ic)

      CO2_AD    = CO2_AD + Predictor_AD%X(k, 6, ic)*SECANG(k)

      Tzp_AD(k) = Tzp_AD(k)                                     &
                  + Predictor_AD%X(k, 7, ic)*SECANG(k)          &
                  + Predictor_AD%X(k, 9, ic)*THREE*PAFV%Tzp(k)**2 &
                  + Predictor_AD%X(k, 10, ic)*SECANG(k)*SQRT(T)

      GAzp_AD(k, ja_co2) = GAzp_AD(k, ja_co2)                   &
                  + Predictor_AD%X(k, 8, ic)*TWO*SECANG(k)**2*PAFV%GAzp(k, ja_co2)


      Predictor_AD%X(k, 1, ic)  = ZERO
      Predictor_AD%X(k, 2, ic)  = ZERO
      Predictor_AD%X(k, 3, ic)  = ZERO
      Predictor_AD%X(k, 4, ic)  = ZERO
      Predictor_AD%X(k, 5, ic)  = ZERO
      Predictor_AD%X(k, 6, ic)  = ZERO
      Predictor_AD%X(k, 7, ic)  = ZERO
      Predictor_AD%X(k, 8, ic)  = ZERO
      Predictor_AD%X(k, 9, ic)  = ZERO
      Predictor_AD%X(k,10, ic)  = ZERO
    END SUBROUTINE AD_Kernel_CO2

    !  --------------------------
    !  Water-line predictors (1 - 15), plus the group-1 extension
    !  (WLO predictors 16 - 18 AND the CO2 component's predictor 11)
    !  when the roster requests 18 WLO predictors. The CO2 predictor-11
    !  adjoint lives here, not in AD_Kernel_CO2, because the heritage
    !  monolith accumulates it at exactly this point in the sequence.
    !  --------------------------
    SUBROUTINE AD_Kernel_WLO( k, ic, np )
      INTEGER, INTENT(IN) :: k, ic, np
      H2O_A_AD = H2O_A_AD                                 &
                 + Predictor_AD%X(k, 1, ic)               &
                 + Predictor_AD%X(k, 2, ic)*DT            &
                 + Predictor_AD%X(k, 4, ic)*DT2           &
                 + Predictor_AD%X(k, 6, ic)*H2O_S

      DT_AD    = DT_AD                                    &
                 + Predictor_AD%X(k, 2, ic)*H2O_A         &
                 + Predictor_AD%X(k, 8, ic)*H2O_R

      H2O_S_AD = H2O_S_AD                                 &
                 + Predictor_AD%X(k, 3, ic)               &
                 + Predictor_AD%X(k, 6, ic)*H2O_A         &
                 + Predictor_AD%X(k, 9, ic)*TWO*H2O_S

      DT2_AD   = DT2_AD + Predictor_AD%X(k, 4, ic)*H2O_A

      H2O_R4_AD= H2O_R4_AD + Predictor_AD%X(k, 5, ic)

      H2O_R_AD = H2O_R_AD                                 &
                 + Predictor_AD%X(k, 7, ic)               &
                 + Predictor_AD%X(k, 8, ic)*DT            &
                 + Predictor_AD%X(k,11, ic)*H2OdH2OTzp

      H2OdH2OTzp_AD = H2OdH2OTzp_AD                       &
                 + Predictor_AD%X(k,10, ic)               &
                 + Predictor_AD%X(k,11, ic)*H2O_R

      GAzp_AD(k,ja_h2o) = GAzp_AD(k,ja_h2o)                                            &
                 + Predictor_AD%X(k,12, ic)*TWO*SECANG(k)**2*PAFV%GAzp(k,ja_h2o)  &
                 + Predictor_AD%X(k,13, ic)*SECANG(k)

      CO2_AD   = CO2_AD + Predictor_AD%X(k,15, ic)*SECANG(k)

      Predictor_AD%X(k, 1, ic) = ZERO
      Predictor_AD%X(k, 2, ic) = ZERO
      Predictor_AD%X(k, 3, ic) = ZERO
      Predictor_AD%X(k, 4, ic) = ZERO
      Predictor_AD%X(k, 5, ic) = ZERO
      Predictor_AD%X(k, 6, ic) = ZERO
      Predictor_AD%X(k, 7, ic) = ZERO
      Predictor_AD%X(k, 8, ic) = ZERO
      Predictor_AD%X(k, 9, ic) = ZERO
      Predictor_AD%X(k,10, ic) = ZERO
      Predictor_AD%X(k,11, ic) = ZERO
      Predictor_AD%X(k,12, ic) = ZERO
      Predictor_AD%X(k,13, ic) = ZERO
      Predictor_AD%X(k,14, ic) = ZERO
      Predictor_AD%X(k,15, ic) = ZERO

      ! Addtional predictors for group 1
      IF ( np >= 18 ) THEN

        CO_A_AD  = CO_A_AD   + Predictor_AD%X(k, 18, ic)          &
                             + Predictor_AD%X(k, 11, ic_co2)
        CH4_A_AD = CH4_A_AD                                       &
                   + Predictor_AD%X(k, 16, ic)                    &
                   + Predictor_AD%X(k, 17, ic)*TWO*CH4_A*DT
        DT_AD    = DT_AD + Predictor_AD%X(k, 17, ic)*CH4_A*CH4_A

        Predictor_AD%X(k, 16, ic) = ZERO
        Predictor_AD%X(k, 17, ic) = ZERO
        Predictor_AD%X(k, 18, ic) = ZERO
        Predictor_AD%X(k, 11, ic_co2) = ZERO

      END IF
    END SUBROUTINE AD_Kernel_WLO

    !  -----------------------
    !  Carbon monoxide
    !  -----------------------
    SUBROUTINE AD_Kernel_CO( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      CO_A_AD = CO_A_AD                                               &
                + Predictor_AD%X(k, 1,  ic)                           &
                + Predictor_AD%X(k, 2,  ic)*DT                        &
                + Predictor_AD%X(k, 7,  ic)*DT2

      DT_AD   = DT_AD                                                 &
                + Predictor_AD%X(k, 2,  ic)*CO_A                      &
                + Predictor_AD%X(k, 4,  ic)*CO_R

      CO_R_AD = CO_R_AD                                               &
                + Predictor_AD%X(k, 3,  ic)*POINT_5/SQRT(CO_R)        &
                + Predictor_AD%X(k, 4,  ic)*DT                        &
                + Predictor_AD%X(k, 6,  ic)                           &
                - Predictor_AD%X(k, 9,  ic)*CO_ACOdCOzp/CO_R**2

      CO_S_AD = CO_S_AD + Predictor_AD%X(k, 5,  ic)

      DT2_AD  = DT2_AD + Predictor_AD%X(k, 7,  ic)*CO_A

      CO_ACOdCOzp_AD = CO_ACOdCOzp_AD                                 &
                 + Predictor_AD%X(k, 8,  ic)                          &
                 + Predictor_AD%X(k, 9,  ic)/CO_R                     &
                 + Predictor_AD%X(k,10,  ic)*SQRT(PAFV%GAzp(k, ja_co))

      GAzp_AD(k, ja_co) = GAzp_AD(k, ja_co)                           &
                 + Predictor_AD%X(k, 10, ic)*                         &
                   POINT_5*CO_ACOdCOzp/SQRT(PAFV%GAzp(k, ja_co))

      Predictor_AD%X(k, 1,  ic)   = ZERO
      Predictor_AD%X(k, 2,  ic)   = ZERO
      Predictor_AD%X(k, 3,  ic)   = ZERO
      Predictor_AD%X(k, 4,  ic)   = ZERO
      Predictor_AD%X(k, 5,  ic)   = ZERO
      Predictor_AD%X(k, 6,  ic)   = ZERO
      Predictor_AD%X(k, 7,  ic)   = ZERO
      Predictor_AD%X(k, 8,  ic)   = ZERO
      Predictor_AD%X(k, 9,  ic)   = ZERO
      Predictor_AD%X(k, 10, ic)   = ZERO
    END SUBROUTINE AD_Kernel_CO

    !  -----------------------
    !  Methane predictors
    !  -----------------------
    SUBROUTINE AD_Kernel_CH4( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      CH4_A_AD = CH4_A_AD                                                &
                 + Predictor_AD%X(k, 1,  ic)*DT                          &
                 + Predictor_AD%X(k, 3,  ic)*TWO*CH4_A                   &
                 + Predictor_AD%X(k, 4,  ic)

      DT_AD    = DT_AD                                                   &
                 + Predictor_AD%X(k, 1,  ic)*CH4_A                       &
                 + Predictor_AD%X(k, 5,  ic)*CH4

      CH4_R_AD = CH4_R_AD                                                &
                 + Predictor_AD%X(k, 2,  ic)                             &
                 + Predictor_AD%X(k, 8,  ic)*POINT_5/SQRT(CH4_R)         &
                 + Predictor_AD%X(k, 11, ic)*CH4/PAFV%GAzp(k, ja_ch4)

      CH4_AD   = CH4_AD                                                  &
                 + Predictor_AD%X(k, 5,  ic)*DT                          &
                 + Predictor_AD%X(k, 11, ic)*CH4_R/PAFV%GAzp(k, ja_ch4)

      CH4_ACH4zp_AD = CH4_ACH4zp_AD                                      &
                 + Predictor_AD%X(k, 6,  ic)                             &
                 + Predictor_AD%X(k, 7,  ic)*TWO*CH4_ACH4zp

      GATzp_AD(k, ja_ch4) = GATzp_AD(k, ja_ch4)                          &
                 + Predictor_AD%X(k, 9,  ic)                             &
                 + Predictor_AD%X(k,10,  ic)*SECANG(k)

      GAzp_AD(k, ja_ch4) = GAzp_AD(k, ja_ch4)                            &
                 - Predictor_AD%X(k, 11, ic)*                            &
                   CH4_R*CH4/PAFV%GAzp(k, ja_ch4)**2

      Predictor_AD%X(k, 1,  ic)  = ZERO
      Predictor_AD%X(k, 2,  ic)  = ZERO
      Predictor_AD%X(k, 3,  ic)  = ZERO
      Predictor_AD%X(k, 4,  ic)  = ZERO
      Predictor_AD%X(k, 5,  ic)  = ZERO
      Predictor_AD%X(k, 6,  ic)  = ZERO
      Predictor_AD%X(k, 7,  ic)  = ZERO
      Predictor_AD%X(k, 8,  ic)  = ZERO
      Predictor_AD%X(k, 9,  ic)  = ZERO
      Predictor_AD%X(k, 10, ic)  = ZERO
      Predictor_AD%X(k, 11, ic)  = ZERO
    END SUBROUTINE AD_Kernel_CH4

    !  -----------------------
    !    N2O predictors
    !  -----------------------
    SUBROUTINE AD_Kernel_N2O( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      N2O_A_AD = N2O_A_AD                                                 &
                 + Predictor_AD%X(k, 1, ic)*DT                            &
                 + Predictor_AD%X(k, 4, ic)*POINT_25*N2O_A**(-POINT_75)   &
                 + Predictor_AD%X(k, 5, ic)

      DT_AD    = DT_AD                                                    &
                 + Predictor_AD%X(k, 1, ic)*N2O_A                         &
                 + Predictor_AD%X(k, 3, ic)*N2O

      N2O_R_AD = N2O_R_AD                                                 &
                 + Predictor_AD%X(k, 2, ic)                               &
                 + Predictor_AD%X(k,10, ic)*N2O/PAFV%GAzp(k, ja_n2o)

      N2O_AD   = N2O_AD                                                   &
                 + Predictor_AD%X(k, 3, ic)*DT                            &
                 + Predictor_AD%X(k,10, ic)*N2O_R/PAFV%GAzp(k, ja_n2o)

      GAzp_AD(k, ja_n2o) = GAzp_AD(k, ja_n2o)                             &
                 + Predictor_AD%X(k, 6, ic)*SECANG(k)                     &
                 - Predictor_AD%X(k,10, ic)*N2O_R*N2O/PAFV%GAzp(k, ja_n2o)**2

      GATzp_AD(k, ja_n2o) = GATzp_AD(k, ja_n2o)                           &
                 + Predictor_AD%X(k, 7, ic)*SECANG(k)                     &
                 + Predictor_AD%X(k, 9, ic)

      N2O_S_AD = N2O_S_AD + Predictor_AD%X(k, 8, ic)

      CH4_A_AD = CH4_A_AD                                                 &
                 + Predictor_AD%X(k,11, ic)                               &
                 + Predictor_AD%X(k,12, ic)*PAFV%GAzp(k, ja_ch4)

      GAzp_AD(k, ja_ch4) = GAzp_AD(k, ja_ch4)                             &
                 + Predictor_AD%X(k,12, ic)*CH4_A

      CO_A_AD = CO_A_AD                                                   &
                + Predictor_AD%X(k,13, ic)                                &
                + Predictor_AD%X(k,14, ic)*SECANG(k)*PAFV%GAzp(k, ja_co)

      GAzp_AD(k, ja_co) = GAzp_AD(k, ja_co)                               &
                + Predictor_AD%X(k,14, ic)*CO_A*SECANG(k)

      Predictor_AD%X(k, 1, ic)   = ZERO
      Predictor_AD%X(k, 2, ic)   = ZERO
      Predictor_AD%X(k, 3, ic)   = ZERO
      Predictor_AD%X(k, 4, ic)   = ZERO
      Predictor_AD%X(k, 5, ic)   = ZERO
      Predictor_AD%X(k, 6, ic)   = ZERO
      Predictor_AD%X(k, 7, ic)   = ZERO
      Predictor_AD%X(k, 8, ic)   = ZERO
      Predictor_AD%X(k, 9, ic)   = ZERO
      Predictor_AD%X(k,10, ic)   = ZERO

      Predictor_AD%X(k,11, ic)   = ZERO
      Predictor_AD%X(k,12, ic)   = ZERO
      Predictor_AD%X(k,13, ic)   = ZERO
      Predictor_AD%X(k,14, ic)   = ZERO
    END SUBROUTINE AD_Kernel_N2O

    !  --------------------------------
    !  Water vapor, MW (line and continuum together)
    !  --------------------------------
    SUBROUTINE AD_Kernel_WET_MW( k, ic )
      INTEGER, INTENT(IN) :: k, ic
      H2O_A_AD = H2O_A_AD                                            &
                 + Predictor_AD%X(k, 1, ic)/T                        &
                 + Predictor_AD%X(k, 2, ic)*H2O/T                    &
                 + Predictor_AD%X(k, 3, ic)*H2O/T2**2                &
                 + Predictor_AD%X(k, 4, ic)/T2                       &
                 + Predictor_AD%X(k, 5, ic)*H2O/T2                   &
                 + Predictor_AD%X(k, 6, ic)/T2**2                    &
                 + Predictor_AD%X(k, 7, ic)                          &
                 + Predictor_AD%X(k, 8, ic)*DT                       &
                 + Predictor_AD%X(k,12, ic)*H2O_S

      T_AD     = T_AD                                                &
                 - Predictor_AD%X(k, 1, ic)*H2O_A/T**2               &
                 - Predictor_AD%X(k, 2, ic)*H2O_A*H2O/T**2

      H2O_AD   = H2O_AD                                              &
                 + Predictor_AD%X(k, 2, ic)*H2O_A/T                  &
                 + Predictor_AD%X(k, 3, ic)*H2O_A/T2**2              &
                 + Predictor_AD%X(k, 5, ic)*H2O_A/T2

      T2_AD    = T2_AD                                               &
                 - Predictor_AD%X(k, 3, ic)*Two*H2O_A*H2O/T2**3      &
                 - Predictor_AD%X(k, 4, ic)*H2O_A/T2**2              &
                 - Predictor_AD%X(k, 5, ic)*H2O_A*H2O/T2**2          &
                 - Predictor_AD%X(k, 6, ic)*TWO*H2O_A/T2**3

      DT_AD    = DT_AD + Predictor_AD%X(k, 8, ic)*H2O_A

      GAzp_AD(k,ja_h2o) = GAzp_AD(k,ja_h2o)                                       &
                 + Predictor_AD%X(k, 9, ic)*TWO*SECANG(k)**2*PAFV%GAzp(k,ja_h2o)  &
                 + Predictor_AD%X(k, 10,ic)*SECANG(k)

      H2O_S_AD = H2O_S_AD                                            &
                 + Predictor_AD%X(k, 12,ic)*H2O_A                    &
                 + Predictor_AD%X(k, 13,ic)*TWO*H2O_S

      H2OdH2OTzp_AD = H2OdH2OTzp_AD + Predictor_AD%X(k, 14,ic)

      Predictor_AD%X(k, 1, ic) = ZERO
      Predictor_AD%X(k, 2, ic) = ZERO
      Predictor_AD%X(k, 3, ic) = ZERO
      Predictor_AD%X(k, 4, ic) = ZERO
      Predictor_AD%X(k, 5, ic) = ZERO
      Predictor_AD%X(k, 6, ic) = ZERO
      Predictor_AD%X(k, 7, ic) = ZERO
      Predictor_AD%X(k, 8, ic) = ZERO
      Predictor_AD%X(k, 9, ic) = ZERO
      Predictor_AD%X(k, 10,ic) = ZERO
      Predictor_AD%X(k, 11,ic) = ZERO
      Predictor_AD%X(k, 12,ic) = ZERO
      Predictor_AD%X(k, 13,ic) = ZERO
      Predictor_AD%X(k, 14,ic) = ZERO
    END SUBROUTINE AD_Kernel_WET_MW

  END SUBROUTINE  ODPS_Compute_Predictor_AD

!============================== END OF AD

!------------------------------------------------------------------------------
!
! NAME:
!       ODPS_Compute_Predictor_ODAS
!
! PURPOSE:
!       Subroutine to compute ODAS predictors for water vapor line
!       absorption.
!
! CALLING SEQUENCE:
!       CALL ODPS_Compute_Predictor_ODAS( &
!         Temperature,    &  ! Input
!         vapor,          &  ! Inpt
!         Level_Pressure, &  ! Input
!         Pressure,       &  ! Input
!         secant_angle,   &  ! Input
!         Alpha,          &  ! Input
!         Alpha_C1,       &  ! Input
!         Alpha_C2,       &  ! Input
!         Predictor)         ! In/Output
!
! INPUT ARGUMENTS:
!
!       Temperature:     Temperature profile
!                        UNITS:      K
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!         Vapor    :     Water vapor mixing ratio profile
!                        UNITS:      g/kg
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!       Level_Pressure : level pressure profile
!                        UNITS:      hPa
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(0:n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!            Pressure :  Layer pressure profile
!                        UNITS:      hPa
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!       secang_angle  :  Secant sensor zenith angle profile
!                        UNITS:      N/A
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!      Alpha, Alpha_C1, Alpha_C2  :  Coefficients for converting water vapor integrated amount
!                        to ODAS water vapor regression space
!                        UNITS:      N/A
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
! IN/OUTPUT ARGUMENTS:
!       Predictor:      Predictor structure containing the integrated absorber
!                       and predictor profiles.
!                       UNITS:      N/A
!                       TYPE:       ODPS_Predictor_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
!------------------------------------------------------------------------------

  SUBROUTINE ODPS_Compute_Predictor_ODAS(&
    Temperature,    &
    Vapor,          &
    Level_Pressure, &
    Pressure,       &
    secant_angle,   &
    Alpha,          &
    Alpha_C1,       &
    Alpha_C2,       &
    Predictor)
     REAL(fp),                  INTENT(IN)     :: Temperature(:)
     REAL(fp),                  INTENT(IN)     :: Vapor(:)
     REAL(fp),                  INTENT(IN)     :: Level_Pressure(0:)
     REAL(fp),                  INTENT(IN)     :: Pressure(:)
     REAL(fp),                  INTENT(IN)     :: secant_angle(:)
     REAL(fp),                  INTENT(IN)     :: Alpha
     REAL(fp),                  INTENT(IN)     :: Alpha_C1
     REAL(fp),                  INTENT(IN)     :: Alpha_C2
     TYPE(ODPS_Predictor_type), INTENT(IN OUT) :: Predictor

     !Local
     REAL(fp), DIMENSION(0:Predictor%n_Layers) :: xL_t, xL_p
     REAL(fp) :: t2, p2, s_t, s_p, Inverse, dPonG, d_Absorber
     REAL(fp) :: Int_vapor_prev, Int_vapor, AveA, ap1
     INTEGER  :: i, k

     ! Regular predictors
     DO k = 1, Predictor%n_Layers
       t2 = Temperature(k)*Temperature(k)
       p2 = Pressure(k)*Pressure(k)
       Predictor%OX(k, 1)  = Temperature(k)
       Predictor%OX(k, 2)  = Pressure(k)
       Predictor%OX(k, 3)  = t2
       Predictor%OX(k, 4)  = p2
       Predictor%OX(k, 5)  = Temperature(k) * Pressure(k)
       Predictor%OX(k, 6)  = t2 * Pressure(k)
       Predictor%OX(k, 7)  = Temperature(k) * p2
       Predictor%OX(k, 8)  = t2 * p2
       Predictor%OX(k, 9)  = Pressure(k)**POINT_25
       Predictor%OX(k, 10) = Vapor(k)
       Predictor%OX(k, 11) = Vapor(k)/t2
       Predictor%OX(k, 14) = secant_angle(k)
     END DO

     ! Integrated predictors
     Int_vapor_prev = ZERO
     s_t = ZERO
     s_p = ZERO
     xL_t(0) = ZERO
     xL_p(0) = ZERO

     DO k = 1, Predictor%n_Layers

       dPonG = RECIPROCAL_GRAVITY * (Level_Pressure(k) - Level_Pressure(k-1))
       d_Absorber = dPonG*Vapor(k)*secant_angle(k)
       Int_vapor = Int_vapor_prev + d_Absorber
       AveA = POINT_5 * (Int_vapor_prev + Int_vapor)

       Predictor%dA(k) = d_Absorber

       s_t = s_t + ( Temperature( k ) * d_Absorber )  ! T*
       s_p = s_p + ( Pressure( k )    * d_Absorber )  ! P*

       IF ( Int_vapor > MINIMUM_ABSORBER_AMOUNT ) THEN
         Inverse = ONE / Int_vapor
       ELSE
         Inverse = ZERO
       END IF

       xL_t(k) = POINT_5 * s_t * Inverse
       xL_p(k) = POINT_5 * s_p * Inverse

       Predictor%OX(k, 12) = xL_t(k) + xL_t(k-1)
       Predictor%OX(k, 13) = xL_p(k) + xL_p(k-1)

       Ap1 =  LOG((aveA - Alpha_C2) / Alpha_C1) / &
          !  ----------------------------------------------
                             Alpha

       Predictor%Ap(k, 1) = Ap1
       DO i = 2, MAX_OPTRAN_ORDER
         Predictor%Ap(k, i) = Predictor%Ap(k, i-1) * Ap1
       END DO

       Int_vapor_prev = Int_vapor

       ! Save variables for TL and AD routines
       IF(Predictor%PAFV%OPTRAN)THEN
         Predictor%PAFV%dPonG(k)      = dPonG
         Predictor%PAFV%d_Absorber(k) = d_Absorber
         Predictor%PAFV%Int_vapor(k)  = Int_vapor
         Predictor%PAFV%AveA(k)       = AveA
         Predictor%PAFV%Inverse(k)    = Inverse
         Predictor%PAFV%s_t(k)        = s_t
         Predictor%PAFV%s_p(k)        = s_p
         Predictor%PAFV%Ap1(k)        = Ap1
       END IF

     END DO

  END SUBROUTINE ODPS_Compute_Predictor_ODAS

!------------------------------------------------------------------------------
!
! NAME:
!       ODPS_Compute_Predictor_ODAS_TL
!
! PURPOSE:
!       Subroutine to compute TL ODAS predictors for water vapor line
!       absorption.
!
! CALLING SEQUENCE:
!    CALL ODPS_Compute_Predictor_ODAS_TL( &
!      Temperature,    &  ! Input
!      vapor,          &  ! Inpt
!      Pressure,       &  ! Input
!      secant_angle,   &  ! Input
!      Alpha,          &  ! Input
!      Alpha_C2,       &  ! Input
!      Predictor,      &  ! Input
!      Temperature_TL, &  ! Input
!      Vapor_TL,       &  ! Input
!      Predictor_TL)      ! Output
!
! INPUT ARGUMENTS:
!
!       Temperature:     Temperature profile
!                        UNITS:      K
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!         Vapor    :     Water vapor mixing ratio profile
!                        UNITS:      g/kg
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!            Pressure :  Layer pressure profile
!                        UNITS:      hPa
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!       secang_angle  :  Secant sensor zenith angle profile
!                        UNITS:      N/A
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!      Alpha, Alpha_C2  :  Coefficients for converting water vapor integrated amount
!                        to ODAS water vapor regression space
!                        UNITS:      N/A
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Predictor:      Predictor structure containing the integrated absorber
!                       and predictor profiles.
!                       UNITS:      N/A
!                       TYPE:       ODPS_Predictor_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Temperature_TL:  TL Temperature profile
!                        UNITS:      K
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!         Vapor_TL    :  TL water vapor mixing ratio profile
!                        UNITS:      g/kg
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
! IN/OUTPUT ARGUMENTS:
!       Predictor_TL:   TL Predictor structure containing the integrated absorber
!                       and predictor profiles.
!                       UNITS:      N/A
!                       TYPE:       ODPS_Predictor_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
!------------------------------------------------------------------------------

  SUBROUTINE ODPS_Compute_Predictor_ODAS_TL( &
    Temperature,    &
    Vapor,          &
    Pressure,       &
    secant_angle,   &
    Alpha,          &
    Alpha_C2,       &
    Predictor,      &
    Temperature_TL, &
    Vapor_TL,       &
    Predictor_TL)

     REAL(fp),                          INTENT(IN)     :: Temperature(:)
     REAL(fp),                          INTENT(IN)     :: Vapor(:)
     REAL(fp),                          INTENT(IN)     :: Pressure(:)
     REAL(fp),                          INTENT(IN)     :: secant_angle(:)
     REAL(fp),                          INTENT(IN)     :: Alpha
     REAL(fp),                          INTENT(IN)     :: Alpha_C2
     TYPE(ODPS_Predictor_type), TARGET, INTENT(IN)     :: Predictor
     REAL(fp),                          INTENT(IN)     :: Temperature_TL(:)
     REAL(fp),                          INTENT(IN)     :: Vapor_TL(:)
     TYPE(ODPS_Predictor_type),         INTENT(IN OUT) :: Predictor_TL

     !Local
     REAL(fp), DIMENSION(0:Predictor%n_Layers) :: xL_t_TL, xL_p_TL
     REAL(fp) :: t2,   p2, &
                 t2_TL,    s_t_TL, s_p_TL, Inverse_TL, d_Absorber_TL
     REAL(fp) :: Int_vapor_prev_TL, Int_vapor_TL, AveA_TL, ap1_TL
     INTEGER  :: i, k
!JR Static initialization means only 1 copy of the variable. OpenMP over profiles
!JR means $OPENMP_NUM_THREADS copies are needed. So change to run-time initialization
!JR     TYPE(PAFV_type), POINTER  :: PAFV => NULL()
     TYPE(PAFV_type), POINTER  :: PAFV

!JR Looks silly since next stmt sets PAFV, but nullify here to guard agains future code changes
     NULLIFY(PAFV)
     ! short name
     PAFV => Predictor%PAFV

     ! Regular predictors
     DO k = 1, Predictor%n_Layers
       t2 = Temperature(k)*Temperature(k)
       p2 = Pressure(k)*Pressure(k)
       t2_TL = TWO*Temperature(k)*Temperature_TL(k)
       Predictor_TL%OX(k, 1)  = Temperature_TL(k)
       Predictor_TL%OX(k, 2)  = ZERO
       Predictor_TL%OX(k, 3)  = t2_TL
       Predictor_TL%OX(k, 4)  = ZERO
       Predictor_TL%OX(k, 5)  = Temperature_TL(k) * Pressure(k)
       Predictor_TL%OX(k, 6)  = t2_TL * Pressure(k)
       Predictor_TL%OX(k, 7)  = Temperature_TL(k) * p2
       Predictor_TL%OX(k, 8)  = t2_TL * p2
       Predictor_TL%OX(k, 9)  = ZERO
       Predictor_TL%OX(k, 10) = Vapor_TL(k)
       Predictor_TL%OX(k, 11) = Vapor_TL(k)/t2 - (Vapor(k)/t2**2)*t2_TL
       Predictor_TL%OX(k, 14) = ZERO
     END DO

     ! Integrated predictors

     Int_vapor_prev_TL = ZERO
     s_t_TL = ZERO
     s_p_TL = ZERO
     xL_t_TL(0) = ZERO
     xL_p_TL(0) = ZERO

     DO k = 1, Predictor%n_Layers

       d_Absorber_TL = PAFV%dPonG(k)*Vapor_TL(k)*secant_angle(k)

       Int_vapor_TL  = Int_vapor_prev_TL + d_Absorber_TL
       AveA_TL = POINT_5 * (Int_vapor_prev_TL + Int_vapor_TL)

       Predictor_TL%dA(k) = d_Absorber_TL

       s_t_TL = s_t_TL + ( Temperature_TL( k )*PAFV%d_Absorber(k) +  Temperature( k )*d_Absorber_TL)
       s_p_TL = s_p_TL + ( Pressure( k )*d_Absorber_TL )

       IF ( PAFV%Int_vapor(k) > MINIMUM_ABSORBER_AMOUNT ) THEN
         Inverse_TL = -(ONE/PAFV%Int_vapor(k)**2)*Int_vapor_TL
       ELSE
         Inverse_TL = ZERO
       END IF

       xL_t_TL(k) = POINT_5 * (s_t_TL*PAFV%Inverse(k) + PAFV%s_t(k)*Inverse_TL)
       xL_p_TL(k) = POINT_5 * (s_p_TL*PAFV%Inverse(k) + PAFV%s_p(k)*Inverse_TL)

       Predictor_TL%OX(k, 12) = xL_t_TL(k) + xL_t_TL(k-1)
       Predictor_TL%OX(k, 13) = xL_p_TL(k) + xL_p_TL(k-1)

       Ap1_TL =           aveA_TL / &
        !        -----------------------------------
                 ( Alpha * (PAFV%aveA(k) - Alpha_C2 ) )

       Predictor_TL%Ap(k, 1) = Ap1_TL

       DO i = 2, MAX_OPTRAN_ORDER
         Predictor_TL%Ap(k, i) = Predictor_TL%Ap(k, i-1)*PAFV%Ap1(k) + Predictor%Ap(k, i-1)*Ap1_TL
       END DO

       Int_vapor_prev_TL = Int_vapor_TL

     END DO

  END SUBROUTINE ODPS_Compute_Predictor_ODAS_TL

!------------------------------------------------------------------------------
!
! NAME:
!       ODPS_Compute_Predictor_ODAS_AD
!
! PURPOSE:
!       Subroutine to compute AD ODAS predictors for water vapor line
!       absorption.
!
! CALLING SEQUENCE:
!    CALL ODPS_Compute_Predictor_ODAS_AD( &
!      Temperature,    &  ! Input
!      vapor,          &  ! Inpt
!      Pressure,       &  ! Input
!      secant_angle,   &  ! Input
!      Alpha,          &  ! Input
!      Alpha_C2,       &  ! Input
!      Predictor,      &  ! Input
!      Predictor_AD,   &  ! Input
!      Temperature_AD, &  ! Input
!      Vapor_AD)          ! Input
!
! INPUT ARGUMENTS:
!
!       Temperature:     Temperature profile
!                        UNITS:      K
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!         Vapor    :     Water vapor mixing ratio profile
!                        UNITS:      g/kg
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!            Pressure :  Layer pressure profile
!                        UNITS:      hPa
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!       secang_angle  :  Secant sensor zenith angle profile
!                        UNITS:      N/A
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(IN)
!
!      Alpha, Alpha_C2  :  Coefficients for converting water vapor integrated amount
!                        to ODAS water vapor regression space
!                        UNITS:      N/A
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Predictor:      Predictor structure containing the integrated absorber
!                       and predictor profiles.
!                       UNITS:      N/A
!                       TYPE:       ODPS_Predictor_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Predictor_aD:   AD Predictor structure containing the integrated absorber
!                       and predictor profiles.
!                       UNITS:      N/A
!                       TYPE:       ODPS_Predictor_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
! IN/OUTPUT ARGUMENTS:
!       Temperature_AD:  AD Temperature profile
!                        UNITS:      K
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(INoUT)
!
!         Vapor_AD    :  AD water vapor mixing ratio profile
!                        UNITS:      g/kg
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Rank-1(n_Layers) array
!                        ATTRIBUTES: INTENT(INOUT)
!
!------------------------------------------------------------------------------

  SUBROUTINE ODPS_Compute_Predictor_ODAS_AD( &
    Temperature,    &
    Vapor,          &
    Pressure,       &
    secant_angle,   &
    Alpha,          &
    Alpha_C2,       &
    Predictor,      &
    Predictor_AD,   &
    Temperature_AD, &
    Vapor_AD)

     REAL(fp),                          INTENT(IN)     :: Temperature(:)
     REAL(fp),                          INTENT(IN)     :: Vapor(:)
     REAL(fp),                          INTENT(IN)     :: Pressure(:)
     REAL(fp),                          INTENT(IN)     :: secant_angle(:)
     REAL(fp),                          INTENT(IN)     :: Alpha
     REAL(fp),                          INTENT(IN)     :: Alpha_C2
     TYPE(ODPS_Predictor_type), TARGET, INTENT(IN)     :: Predictor
     TYPE(ODPS_Predictor_type),         INTENT(IN OUT) :: Predictor_AD
     REAL(fp),                          INTENT(IN OUT) :: Temperature_AD(:)
     REAL(fp),                          INTENT(IN OUT) :: Vapor_AD(:)

     !Local
     REAL(fp), DIMENSION(0:Predictor%n_Layers) :: xL_t_AD, xL_p_AD
     REAL(fp) :: t2, p2
     REAL(fp) :: t2_AD, s_t_AD, s_p_AD, Inverse_AD, d_Absorber_AD
     REAL(fp) :: Int_vapor_prev_AD, Int_vapor_AD, AveA_AD, ap1_AD
     INTEGER  :: i, k
!JR Static initialization means only 1 copy of the variable. OpenMP over profiles
!JR means $OPENMP_NUM_THREADS copies are needed. So change to run-time initialization
!JR     TYPE(PAFV_type), POINTER  :: PAFV => NULL()
     TYPE(PAFV_type), POINTER  :: PAFV

!JR Looks silly since next stmt sets PAFV, but nullify here to guard agains future code changes
     NULLIFY(PAFV)
     ! short name
     PAFV => Predictor%PAFV

     ! Integrated predictors

     ! --- AD part
     Int_vapor_prev_AD = ZERO
     Int_vapor_AD      = ZERO
     Ap1_AD            = ZERO
     AveA_AD           = ZERO
     xL_t_AD(Predictor%n_Layers) = ZERO
     xL_p_AD(Predictor%n_Layers) = ZERO
     s_t_AD            = ZERO
     s_p_AD            = ZERO
     Inverse_AD        = ZERO
     d_Absorber_AD     = ZERO
     DO k = Predictor%n_Layers, 1, -1

       Int_vapor_AD = Int_vapor_AD + Int_vapor_prev_AD
       Int_vapor_prev_AD = ZERO

       DO i = MAX_OPTRAN_ORDER, 2, -1
         Ap1_AD = Ap1_AD + Predictor%Ap(k, i-1)*Predictor_AD%Ap(k, i)
         Predictor_AD%Ap(k, i-1) = Predictor_AD%Ap(k, i-1) + PAFV%Ap1(k)*Predictor_AD%Ap(k, i)
         Predictor_AD%Ap(k, i) = ZERO
       END DO

       Ap1_AD = Ap1_AD + Predictor_AD%Ap(k, 1)
       Predictor_AD%Ap(k, 1) = ZERO

       aveA_AD = aveA_AD + &
                            Ap1_AD / &
        !        -----------------------------------
                 ( Alpha * (PAFV%aveA(k) - Alpha_C2 ) )
       Ap1_AD = ZERO

       xL_t_AD(k)   = xL_t_AD(k)   + Predictor_AD%OX(k, 12)
       xL_t_AD(k-1) =                Predictor_AD%OX(k, 12)  ! combine with initialization for xL_t_AD(k-1)
       Predictor_AD%OX(k, 12) = ZERO
       xL_p_AD(k)   = xL_p_AD(k)   + Predictor_AD%OX(k, 13)
       xL_p_AD(k-1) =                Predictor_AD%OX(k, 13) ! combine with initialization for xL_p_AD(k-1)
       Predictor_AD%OX(k, 13) = ZERO

       s_p_AD     = s_p_AD     + POINT_5*PAFV%Inverse(k)*xL_p_AD(k)
       Inverse_AD = Inverse_AD + POINT_5*PAFV%s_p(k)*xL_p_AD(k)
       s_t_AD     = s_t_AD     + POINT_5*PAFV%Inverse(k)*xL_t_AD(k)
       Inverse_AD = Inverse_AD + POINT_5*PAFV%s_t(k)*xL_t_AD(k)
       xL_t_AD(k) = ZERO
       xL_p_AD(k) = ZERO

       IF ( PAFV%Int_vapor(k) > MINIMUM_ABSORBER_AMOUNT ) THEN
         Int_vapor_AD = Int_vapor_AD -(ONE/PAFV%Int_vapor(k)**2)*Inverse_AD
         Inverse_AD = ZERO
       ELSE
         Inverse_AD = ZERO
       END IF

       d_Absorber_AD = d_Absorber_AD + Pressure( k )*s_p_AD &
                                     + Temperature( k )*s_t_AD
       Temperature_AD( k ) = Temperature_AD( k ) + PAFV%d_Absorber(k)*s_t_AD

       d_Absorber_AD = d_Absorber_AD + Predictor_AD%dA(k)
       Predictor_AD%dA(k) = ZERO

       Int_vapor_prev_AD = Int_vapor_prev_AD + POINT_5 * AveA_AD
       Int_vapor_AD = Int_vapor_AD + POINT_5 *  AveA_AD
       AveA_AD = ZERO

       Int_vapor_prev_AD = Int_vapor_prev_AD + Int_vapor_AD
       d_Absorber_AD = d_Absorber_AD + Int_vapor_AD
       Int_vapor_AD = ZERO

       Vapor_AD(k) = Vapor_AD(k) + PAFV%dPonG(k)*d_Absorber_AD*secant_angle(k)
       d_Absorber_AD = ZERO

     END DO


     ! AD Regular predictors
     DO k = Predictor%n_Layers, 1, -1
       t2 = Temperature(k)*Temperature(k)
       p2 = Pressure(k)*Pressure(k)

       Temperature_AD(k) = Temperature_AD(k) &
                           +              Predictor_AD%OX(k, 1) &
                           + Pressure(k)* Predictor_AD%OX(k, 5) &
                           +          p2* Predictor_AD%OX(k, 7)
       t2_AD =                            Predictor_AD%OX(k, 3) &
               +             Pressure(k)* Predictor_AD%OX(k, 6) &
               +                      p2* Predictor_AD%OX(k, 8) &
               -        (Vapor(k)/t2**2)* Predictor_AD%OX(k, 11)

       Vapor_AD(k) = Vapor_AD(k) &
                     + Predictor_AD%OX(k, 10) &
                     + Predictor_AD%OX(k, 11)/t2

       Predictor_AD%OX(k, 1:11) = ZERO
       Predictor_AD%OX(k,   14) = ZERO

       Temperature_AD(k) = Temperature_AD(k) + TWO*Temperature(k)*t2_AD

     END DO

  END SUBROUTINE ODPS_Compute_Predictor_ODAS_AD


  ! Registry accessors. Out-of-range or Zeeman-reserved group queries return
  ! 0 (dimensions) or ODPS_INVALID_ID (IDs); ODPS_Validate_Group is the
  ! load-time gate that makes such queries unreachable in normal operation.

  PURE FUNCTION ODPS_Get_max_n_Predictors( Group_Index ) RESULT( max_n_Predictors )
    INTEGER, INTENT( IN ) :: Group_Index
    INTEGER :: max_n_Predictors
    IF ( Group_Index >= 1 .AND. Group_Index <= N_G ) THEN
      max_n_Predictors = GROUP_REGISTRY(Group_Index)%Max_n_Predictors
    ELSE
      max_n_Predictors = 0
    END IF
  END FUNCTION ODPS_Get_max_n_Predictors


  PURE FUNCTION ODPS_Get_n_Components( Group_Index ) RESULT( n_Components )
    INTEGER, INTENT( IN ) :: Group_Index
    INTEGER :: n_Components
    IF ( Group_Index >= 1 .AND. Group_Index <= N_G ) THEN
      n_Components = GROUP_REGISTRY(Group_Index)%n_Components
    ELSE
      n_Components = 0
    END IF
  END FUNCTION ODPS_Get_n_Components


  PURE FUNCTION ODPS_Get_n_Absorbers( Group_Index ) RESULT( n_Absorbers )
    INTEGER, INTENT( IN ) :: Group_Index
    INTEGER :: n_Absorbers
    IF ( Group_Index >= 1 .AND. Group_Index <= N_G ) THEN
      n_Absorbers = GROUP_REGISTRY(Group_Index)%n_Absorbers
    ELSE
      n_Absorbers = 0
    END IF
  END FUNCTION ODPS_Get_n_Absorbers


  PURE FUNCTION ODPS_Get_Component_ID(Component_Index, Group_Index) RESULT( Component_ID )
    INTEGER, INTENT( IN ) :: Component_Index
    INTEGER, INTENT( IN ) :: Group_Index
    INTEGER :: Component_ID
    IF ( Group_Index >= 1 .AND. Group_Index <= N_G ) THEN
      IF ( Component_Index >= 1 .AND. &
           Component_Index <= GROUP_REGISTRY(Group_Index)%n_Components ) THEN
        Component_ID = GROUP_REGISTRY(Group_Index)%Component_ID(Component_Index)
        RETURN
      END IF
    END IF
    Component_ID = ODPS_INVALID_ID  ! Out-of-range query; see ODPS_Validate_Group
  END FUNCTION ODPS_Get_Component_ID


  PURE FUNCTION ODPS_Get_Absorber_ID(Absorber_Index, Group_Index) RESULT(  Absorber_ID )
    INTEGER, INTENT( IN ) :: Absorber_Index
    INTEGER, INTENT( IN ) :: Group_Index
    INTEGER :: Absorber_ID
    IF ( Group_Index >= 1 .AND. Group_Index <= N_G ) THEN
      IF ( Absorber_Index >= 1 .AND. &
           Absorber_Index <= GROUP_REGISTRY(Group_Index)%n_Absorbers ) THEN
        Absorber_ID = GROUP_REGISTRY(Group_Index)%Absorber_ID(Absorber_Index)
        RETURN
      END IF
    END IF
    Absorber_ID = ODPS_INVALID_ID  ! Out-of-range query; see ODPS_Validate_Group
  END FUNCTION ODPS_Get_Absorber_ID


  PURE FUNCTION ODPS_Get_Ozone_Component_ID(Group_Index) RESULT( Ozone_Component_ID )
    INTEGER, INTENT(IN) :: Group_Index
    INTEGER :: Ozone_Component_ID
    IF( Group_Index == GROUP_1 .OR. Group_Index == GROUP_2 .OR. &
        Group_Index == GROUP_MW_O3 .OR. Group_Index == GROUP_UV_NO2 )THEN
      Ozone_Component_ID = OZO_ComID
    ELSE
      Ozone_Component_ID = ODPS_INVALID_ID
    END IF
  END FUNCTION ODPS_Get_Ozone_Component_ID


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       ODPS_Validate_Group
!
! PURPOSE:
!       Validate the group metadata of a loaded ODPS TauCoeff structure
!       against the supported group definitions: the Group_Index must be a
!       supported (non-Zeeman-reserved) group, and the file's Component_ID
!       and Absorber_ID rosters must match that group's compiled maps
!       exactly, in content and order (the predictor code addresses
!       components positionally).
!
! CALLING SEQUENCE:
!       Is_Valid = ODPS_Validate_Group( Group_Index,  &
!                                       Component_ID, &
!                                       Absorber_ID,  &
!                                       Message       )
!
! INPUTS:
!       Group_Index:   The file's ODPS group index.
!       Component_ID:  The file's component ID roster.
!       Absorber_ID:   The file's absorber ID roster.
!
! OUTPUTS:
!       Message:       Explanation of the failure (blank when valid).
!
! FUNCTION RESULT:
!       Is_Valid:      .TRUE. when the metadata is consistent with a
!                      supported group; .FALSE. otherwise.
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION ODPS_Validate_Group( Group_Index, Component_ID, Absorber_ID, Message ) &
    RESULT( Is_Valid )
    INTEGER,      INTENT(IN)  :: Group_Index
    INTEGER,      INTENT(IN)  :: Component_ID(:)
    INTEGER,      INTENT(IN)  :: Absorber_ID(:)
    CHARACTER(*), INTENT(OUT) :: Message
    LOGICAL :: Is_Valid
    ! Local variables
    INTEGER :: nc, na, i, n_trace
    LOGICAL :: Has_Trace

    Is_Valid = .FALSE.
    Message  = ''

    ! Reserved (Zeeman) indexes get a specific explanation
    IF ( Group_Index >= RESERVED_ZSSMIS_GROUP .AND. &
         Group_Index <= RESERVED_ZAMSUA_GROUP + 1 ) THEN
      WRITE( Message, '("Group_Index ",i0," is reserved for the Zeeman ", &
        &"sub-algorithms. Zeeman companion coefficients (z*.TauCoeff) are ", &
        &"loaded via the ODZeeman path for SSMIS/AMSU-A sensors only; a ", &
        &"standard ODPS TauCoeff must use group 1, 2, 3, 7, 8, 9, or 10.")' ) &
        Group_Index
      RETURN
    END IF

    ! Unknown group index
    IF ( .NOT. ANY( VALID_GROUPS == Group_Index ) ) THEN
      WRITE( Message, '("Group_Index ",i0," is not a supported ODPS group ", &
        &"(supported: 1, 2, 3, 7, 8, 9, 10).")' ) Group_Index
      RETURN
    END IF

    ! Kernel-capability validation. A file's roster is valid when every
    ! component ID maps to a compiled predictor kernel for this group's
    ! basis and every gas those kernels consume is present in the file's
    ! absorber roster. Order and subsetting are free: the compute path
    ! dispatches by the file's own roster. Unknown component IDs (for
    ! example raw molecule-set 13/14 from externally trained files) are
    ! rejected: a kernel is a trained CRTM predictor formulation, not
    ! just a gas label.
    nc = SIZE(Component_ID)
    na = SIZE(Absorber_ID)
    IF ( nc < 1 .OR. na < 1 ) THEN
      WRITE( Message, '("Empty roster: ",i0," components, ",i0, &
        &" absorbers.")' ) nc, na
      RETURN
    END IF

    ! Duplicates
    DO i = 1, nc
      IF ( ANY( Component_ID(1:i-1) == Component_ID(i) ) ) THEN
        WRITE( Message, '("Duplicate Component_ID ",i0,".")' ) Component_ID(i)
        RETURN
      END IF
    END DO
    DO i = 1, na
      IF ( ANY( Absorber_ID(1:i-1) == Absorber_ID(i) ) ) THEN
        WRITE( Message, '("Duplicate Absorber_ID ",i0,".")' ) Absorber_ID(i)
        RETURN
      END IF
    END DO

    ! Known absorbers only
    DO i = 1, na
      IF ( .NOT. ANY( KNOWN_GAS_IDS == Absorber_ID(i) ) ) THEN
        WRITE( Message, '("Absorber_ID ",i0," is not a gas CRTM knows ", &
          &"(known HITRAN IDs: 1, 2, 3, 4, 5, 6, 10).")' ) Absorber_ID(i)
        RETURN
      END IF
    END DO

    ! The group-1 trace components travel together, and their extension
    ! predictors live in the WLO and CO2 kernels
    n_trace = COUNT( (/ ANY(Component_ID == CO_ComID),  &
                        ANY(Component_ID == CH4_ComID), &
                        ANY(Component_ID == N2O_ComID) /) )
    Has_Trace = ( n_trace == 3 )
    IF ( n_trace > 0 .AND. .NOT. Has_Trace ) THEN
      Message = 'The trace components (CO 119, CH4 118, N2O 120) must '// &
                'appear together or not at all.'
      RETURN
    END IF
    IF ( Has_Trace .AND. .NOT. ( ANY(Component_ID == WLO_ComID) .AND. &
                                 ANY(Component_ID == CO2_ComID) ) ) THEN
      Message = 'A roster with the trace components also requires the '// &
                'WLO (101) and CO2 (121) components (their kernels carry '// &
                'the trace cross-term predictors).'
      RETURN
    END IF

    ! Every component must have a kernel for this basis, and every gas its
    ! kernel consumes must be in the absorber roster
    DO i = 1, nc
      IF ( ODPS_Kernel_n_Predictors( GROUP_REGISTRY(Group_Index)%Basis, &
                                     Component_ID(i), Has_Trace ) == 0 ) THEN
        WRITE( Message, '("Component_ID ",i0," has no predictor kernel for ", &
          &"the ",a," basis (supported IR: 7, 20, 101, 15, 114, 121, 120, ", &
          &"119, 118, 122, 123, 124; MW: 113, 12, 114).")' ) Component_ID(i), &
          TRIM(Basis_Name(GROUP_REGISTRY(Group_Index)%Basis))
        RETURN
      END IF
      IF ( .NOT. Required_Gases_Present( Component_ID(i), Absorber_ID ) ) THEN
        WRITE( Message, '("Component_ID ",i0," requires a gas that is not ", &
          &"in the file''s absorber roster.")' ) Component_ID(i)
        RETURN
      END IF
    END DO

    Is_Valid = .TRUE.

  CONTAINS

    PURE FUNCTION Basis_Name( Basis ) RESULT( Name )
      INTEGER, INTENT(IN) :: Basis
      CHARACTER(8) :: Name
      SELECT CASE ( Basis )
        CASE ( BASIS_IR ); Name = 'IR/UV'
        CASE ( BASIS_MW ); Name = 'MW'
        CASE DEFAULT;      Name = 'RESERVED'
      END SELECT
    END FUNCTION Basis_Name

    PURE FUNCTION Required_Gases_Present( Com_ID, Gases ) RESULT( Present )
      INTEGER, INTENT(IN) :: Com_ID
      INTEGER, INTENT(IN) :: Gases(:)
      LOGICAL :: Present
      SELECT CASE ( Com_ID )
        CASE ( WLO_ComID )
          Present = ANY(Gases == H2O_ID) .AND. ANY(Gases == CO2_ID)
        CASE ( WCO_ComID, WET_ComID )
          Present = ANY(Gases == H2O_ID)
        CASE ( OZO_ComID )
          Present = ANY(Gases == O3_ID)
        CASE ( CO2_ComID )
          Present = ANY(Gases == CO2_ID)
        CASE ( N2O_ComID )
          Present = ANY(Gases == N2O_ID) .AND. ANY(Gases == CH4_ID) &
                    .AND. ANY(Gases == CO_ID)
        CASE ( CO_ComID )
          Present = ANY(Gases == CO_ID)
        CASE ( CH4_ComID )
          Present = ANY(Gases == CH4_ID)
        CASE ( NO2_ComID )
          Present = ANY(Gases == NO2_ID)
        CASE ( DRY_ComID_Q1, DRY_ComID_Q2 )   ! QDRY dry carries water predictors
          Present = ANY(Gases == H2O_ID)
        CASE DEFAULT   ! heritage dry components consume no variable gas
          Present = .TRUE.
      END SELECT
    END FUNCTION Required_Gases_Present

  END FUNCTION ODPS_Validate_Group


!------------------------------------------------------------------------------
! Kernel capability: the number of predictors CRTM's compiled kernel
! computes for a component ID under a given basis (0 = no kernel; the
! component is unsupported). Has_Trace selects the group-1 style variants
! (WLO 18 vs 15, CO2 11 vs 10) that add trace-gas cross terms.
!------------------------------------------------------------------------------
  PURE FUNCTION ODPS_Kernel_n_Predictors( Basis, Component_ID, Has_Trace ) &
    RESULT( n_Predictors )
    INTEGER, INTENT(IN) :: Basis
    INTEGER, INTENT(IN) :: Component_ID
    LOGICAL, INTENT(IN) :: Has_Trace
    INTEGER :: n_Predictors
    n_Predictors = 0
    SELECT CASE ( Basis )
      CASE ( BASIS_IR )
        SELECT CASE ( Component_ID )
          CASE ( DRY_ComID_G1, DRY_ComID_G2 ); n_Predictors = 7
          CASE ( DRY_ComID_Q1, DRY_ComID_Q2 ); n_Predictors = 20
          CASE ( WLO_ComID );  n_Predictors = MERGE( 18, 15, Has_Trace )
          CASE ( WCO_ComID );  n_Predictors = 7
          CASE ( OZO_ComID );  n_Predictors = 11
          CASE ( CO2_ComID );  n_Predictors = MERGE( 11, 10, Has_Trace )
          CASE ( N2O_ComID );  n_Predictors = 14
          CASE ( CO_ComID );   n_Predictors = 10
          CASE ( CH4_ComID );  n_Predictors = 11
          CASE ( NO2_ComID );  n_Predictors = 3
        END SELECT
      CASE ( BASIS_MW )
        SELECT CASE ( Component_ID )
          CASE ( EDRY_ComID ); n_Predictors = 7
          CASE ( WET_ComID );  n_Predictors = 14
          CASE ( OZO_ComID );  n_Predictors = 11
        END SELECT
    END SELECT
  END FUNCTION ODPS_Kernel_n_Predictors


!------------------------------------------------------------------------------
! The predictor-array capacity a file's roster needs: the maximum kernel
! predictor count over its components. Returns 0 for an unsupported group
! or any unsupported component (the load-time validation reports the
! specific reason).
!------------------------------------------------------------------------------
  PURE FUNCTION ODPS_Max_n_Predictors_For( Group_Index, Component_ID ) &
    RESULT( Max_n_Predictors )
    INTEGER, INTENT(IN) :: Group_Index
    INTEGER, INTENT(IN) :: Component_ID(:)
    INTEGER :: Max_n_Predictors
    INTEGER :: i, np
    LOGICAL :: Has_Trace
    Max_n_Predictors = 0
    IF ( Group_Index < 1 .OR. Group_Index > N_G ) RETURN
    Has_Trace = ANY( Component_ID == CO_ComID )
    DO i = 1, SIZE(Component_ID)
      np = ODPS_Kernel_n_Predictors( GROUP_REGISTRY(Group_Index)%Basis, &
                                     Component_ID(i), Has_Trace )
      IF ( np == 0 ) THEN
        Max_n_Predictors = 0
        RETURN
      END IF
      Max_n_Predictors = MAX( Max_n_Predictors, np )
    END DO
  END FUNCTION ODPS_Max_n_Predictors_For


  ! This function gets a flag (true or false) indicating the
  ! need for saveing the FWD variables
  PURE FUNCTION ODPS_Get_SaveFWVFlag() RESULT(Flag)
    LOGICAL :: Flag
    Flag = .TRUE.
  END FUNCTION ODPS_Get_SaveFWVFlag

END MODULE ODPS_Predictor
