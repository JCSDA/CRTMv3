!
! test_CONST_MIXED_Polarization.f90
!
! Unit test for the CONST_MIXED_POLARIZATION (=13) surface-optics polarization
! mixing in CRTM_Compute_SfcOptics.
!
! Background
! ----------
! For a "constant mixed" polarization channel the surface emissivity and
! reflectivity are mixed between the vertical (V) and horizontal (H) components
! using the *fixed* per-channel polarization angle PolAngle:
!
!     e = eV * sin^2(PolAngle) + eH * (1 - sin^2(PolAngle))
!
! A prior bug scaled PolAngle by GeometryInfo%Distance_Ratio
! (= sin(scan)/sin(zenith)), i.e. SIN2_Angle = (Distance_Ratio*sin(PolAngle))^2.
! Distance_Ratio is only meaningful for the V/H-mixed cases, where it converts
! the local zenith angle to the scan angle. PolAngle is a constant, so the mix
! must NOT depend on Distance_Ratio (hence "CONST" mixed). This test pins that
! behaviour for TMS (TROPICS / tomorrow.io) sensors, the only sensors that use
! polarization type 13.
!
! Strategy (per channel, sea-water surface so only the FASTEM MW-water model
! supplies eV/eH):
!   (A) Distance_Ratio invariance  -- the regression gate for the bug. Holding
!       the zenith angle fixed and varying ONLY Distance_Ratio must leave the
!       mixed emissivity/reflectivity unchanged. (The buggy code changed with
!       Distance_Ratio; the corrected code does not.)
!   (B) PolAngle sensitivity       -- non-vacuity guard. Overriding PolAngle
!       0 deg vs 90 deg must change the result, proving PolAngle is actually
!       used and that eV /= eH (which is what makes (A) a meaningful gate).
!   (C) Formula reconstruction     -- structural correctness. The mixed value
!       must equal eV*sin^2(PolAngle) + eH*(1-sin^2(PolAngle)), with eV obtained
!       at PolAngle=90 deg and eH at PolAngle=0 deg.
!
PROGRAM test_CONST_MIXED_Polarization

  ! ============================================================================
  ! **** ENVIRONMENT SETUP ****
  USE UnitTest_Define,          ONLY: UnitTest_type
  USE Type_Kinds,               ONLY: fp
  USE Message_Handler,          ONLY: SUCCESS
  USE CRTM_Parameters,          ONLY: MAX_N_STOKES
  USE CRTM_Surface_Define,      ONLY: CRTM_Surface_type
  USE CRTM_GeometryInfo_Define, ONLY: CRTM_GeometryInfo_type,     &
                                      CRTM_GeometryInfo_SetValue, &
                                      CRTM_GeometryInfo_Destroy
  USE CRTM_SpcCoeff,            ONLY: SC,                          &
                                      CRTM_SpcCoeff_Load,          &
                                      CRTM_SpcCoeff_Destroy,       &
                                      SpcCoeff_IsMicrowaveSensor,  &
                                      CONST_MIXED_POLARIZATION
  USE CRTM_SfcOptics_Define,    ONLY: CRTM_SfcOptics_type,   &
                                      CRTM_SfcOptics_Create, &
                                      CRTM_SfcOptics_Destroy
  USE CRTM_MWwaterCoeff,        ONLY: CRTM_MWwaterCoeff_Load_FASTEM
  USE CRTM_SfcOptics,           ONLY: iVar_type, CRTM_Compute_SfcOptics

  IMPLICIT NONE

  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME      = 'test_CONST_MIXED_Polarization'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  CHARACTER(*), PARAMETER :: SENSOR_ID         = 'tms_tropics-01'  ! polarization type 13
  LOGICAL,      PARAMETER :: QUIET             = .TRUE.

  REAL(fp), PARAMETER :: DEG2RAD = ACOS(-1.0_fp) / 180.0_fp

  ! Geometry. Distance_Ratio is varied between two well-separated values while
  ! the zenith angle (which alone drives eV/eH) is held fixed.
  REAL(fp), PARAMETER :: ZENITH_ANGLE   = 30.0_fp
  REAL(fp), PARAMETER :: DISTANCE_RATIO_1 = 0.50_fp
  REAL(fp), PARAMETER :: DISTANCE_RATIO_2 = 1.00_fp

  ! Tolerances
  REAL(fp), PARAMETER :: INV_TOL     = 1.0e-12_fp  ! (A),(C): expected bit-exact
  REAL(fp), PARAMETER :: SENS_TOL    = 1.0e-4_fp   ! (B): eV and eH must differ by > this

  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256), DIMENSION(1) :: Sensor_Id_Arr
  INTEGER  :: Error_Status
  INTEGER  :: n_Channels
  INTEGER  :: SensorIndex, ChannelIndex
  REAL(fp) :: e_d1, e_d2            ! mixed emissivity at Distance_Ratio 1 / 2
  REAL(fp) :: r_d1, r_d2            ! mixed reflectivity at Distance_Ratio 1 / 2
  REAL(fp) :: e_real, e_v, e_h      ! emissivity at real / 90deg / 0deg PolAngle
  REAL(fp) :: pa_deg, sin2          ! channel PolAngle [deg] and sin^2(PolAngle)
  REAL(fp) :: e_formula             ! reconstructed mix
  REAL(fp) :: pa_save               ! original PolAngle for restore

  TYPE(UnitTest_type)          :: test
  TYPE(iVar_type)              :: iVar
  TYPE(CRTM_Surface_type)      :: Sfc(1)
  TYPE(CRTM_GeometryInfo_type) :: gInfo(1)
  TYPE(CRTM_SfcOptics_type)    :: SfcOptics

  ! ============================================================================
  ! 1. **** INITIALISE UNIT TEST ****
  CALL test%Init(.TRUE.)
  CALL test%Setup(PROGRAM_NAME, PROGRAM_NAME, .TRUE.)

  ! ============================================================================
  ! 2. **** LOAD COEFFICIENTS ****
  SensorIndex      = 1
  Sensor_Id_Arr(1) = SENSOR_ID
  Error_Status = CRTM_SpcCoeff_Load( Sensor_Id_Arr,                 &
                                     File_Path = COEFFICIENTS_PATH, &
                                     netCDF    = .TRUE.,            &
                                     Quiet     = QUIET )
  CALL test%Assert(Error_Status == SUCCESS)
  IF ( Error_Status /= SUCCESS ) THEN
    WRITE(*,'(/5x,"Could not load SpcCoeff for ",a,"; aborting.")') TRIM(SENSOR_ID)
    STOP 1
  END IF

  ! This test only makes sense for a microwave sensor that uses CONST_MIXED.
  CALL test%Assert( SpcCoeff_IsMicrowaveSensor( SC(SensorIndex) ) )

  Error_Status = CRTM_MWwaterCoeff_Load_FASTEM( 'FASTEM6', Quiet = QUIET )
  CALL test%Assert(Error_Status == SUCCESS)

  n_Channels = SC(SensorIndex)%n_Channels
  WRITE(*,'(/5x,"Sensor ",a," : ",i0," channels")') TRIM(SENSOR_ID), n_Channels

  ! ============================================================================
  ! 3. **** SET UP SURFACE + GEOMETRY ****
  ! 100% sea-water so only Compute_MW_Water_SfcOptics (FASTEM) is invoked, and
  ! eV /= eH (water is strongly polarizing).
  Sfc(1)%Water_Coverage    = 1.0_fp
  Sfc(1)%Water_Type        = 1          ! sea water
  Sfc(1)%Water_Temperature = 290.0_fp
  Sfc(1)%Wind_Speed        = 5.0_fp
  Sfc(1)%Wind_Direction    = 0.0_fp

  CALL CRTM_GeometryInfo_SetValue( gInfo,                          &
                                   Source_Azimuth_Angle = 0.0_fp,  &
                                   Sensor_Azimuth_Angle = 0.0_fp,  &
                                   Sensor_Scan_Angle    = 0.0_fp,  &
                                   Sensor_Zenith_Angle  = ZENITH_ANGLE )

  ! Allocate with the full Stokes dimension (room for the V/H emissivity
  ! components the MW-water model fills), then set the n_Stokes *flag* to 1.
  ! n_Stokes==1 is what selects the "decoupled polarization" branch in
  ! CRTM_Compute_SfcOptics where the CONST_MIXED_POLARIZATION mixing lives --
  ! this mirrors the real forward-model setup (CRTM_Forward_Module: Create with
  ! MAX_N_STOKES, then SfcOptics%n_Stokes = RTV%n_Stokes).
  CALL CRTM_SfcOptics_Create( SfcOptics, 1, MAX_N_STOKES )  ! n_Angles, n_Stokes
  SfcOptics%n_Stokes  = 1
  SfcOptics%Angle(1)  = ZENITH_ANGLE
  SfcOptics%Weight(1) = 1.0_fp

  ! ============================================================================
  ! 4. **** PER-CHANNEL TESTS ****
  ChannelLoop: DO ChannelIndex = 1, n_Channels

    ! Guard: every channel of this sensor must be CONST_MIXED_POLARIZATION.
    CALL test%Assert( SC(SensorIndex)%Polarization(ChannelIndex) &
                      == CONST_MIXED_POLARIZATION )

    ! ------------------------------------------------------------------
    ! (A) Distance_Ratio invariance  -- the regression gate for the bug.
    ! ------------------------------------------------------------------
    gInfo(1)%Distance_Ratio = DISTANCE_RATIO_1
    Error_Status = CRTM_Compute_SfcOptics( Sfc(1), gInfo(1), SensorIndex, &
                                           ChannelIndex, SfcOptics, iVar )
    CALL test%Assert(Error_Status == SUCCESS)
    e_d1 = SfcOptics%Emissivity(1,1)
    r_d1 = SfcOptics%Reflectivity(1,1,1,1)

    gInfo(1)%Distance_Ratio = DISTANCE_RATIO_2
    Error_Status = CRTM_Compute_SfcOptics( Sfc(1), gInfo(1), SensorIndex, &
                                           ChannelIndex, SfcOptics, iVar )
    CALL test%Assert(Error_Status == SUCCESS)
    e_d2 = SfcOptics%Emissivity(1,1)
    r_d2 = SfcOptics%Reflectivity(1,1,1,1)

    ! Mixed emissivity / reflectivity must be independent of Distance_Ratio.
    CALL test%Assert_EqualWithin( e_d1, e_d2, INV_TOL )
    CALL test%Assert_EqualWithin( r_d1, r_d2, INV_TOL )

    ! Sanity: a physical emissivity.
    CALL test%Assert( e_d2 > 0.0_fp .AND. e_d2 <= 1.0_fp )

    ! ------------------------------------------------------------------
    ! (B)/(C) PolAngle sensitivity + formula reconstruction.
    !         Distance_Ratio held fixed (proven irrelevant by (A)).
    ! ------------------------------------------------------------------
    pa_save = SC(SensorIndex)%PolAngle(ChannelIndex)
    pa_deg  = pa_save
    sin2    = SIN(pa_deg*DEG2RAD)**2

    ! Mixed value at the channel's real PolAngle.
    e_real = e_d2

    ! eH : PolAngle = 0 deg  -> sin^2 = 0 -> e = eH
    SC(SensorIndex)%PolAngle(ChannelIndex) = 0.0_fp
    Error_Status = CRTM_Compute_SfcOptics( Sfc(1), gInfo(1), SensorIndex, &
                                           ChannelIndex, SfcOptics, iVar )
    CALL test%Assert(Error_Status == SUCCESS)
    e_h = SfcOptics%Emissivity(1,1)

    ! eV : PolAngle = 90 deg -> sin^2 = 1 -> e = eV
    SC(SensorIndex)%PolAngle(ChannelIndex) = 90.0_fp
    Error_Status = CRTM_Compute_SfcOptics( Sfc(1), gInfo(1), SensorIndex, &
                                           ChannelIndex, SfcOptics, iVar )
    CALL test%Assert(Error_Status == SUCCESS)
    e_v = SfcOptics%Emissivity(1,1)

    ! Restore the real PolAngle.
    SC(SensorIndex)%PolAngle(ChannelIndex) = pa_save

    ! (B) The V and H limits must differ (else the mix is vacuous).
    CALL test%Refute_EqualWithin( e_h, e_v, SENS_TOL )

    ! (C) The real mix must reconstruct from eV, eH and sin^2(PolAngle).
    e_formula = e_v*sin2 + e_h*(1.0_fp - sin2)
    CALL test%Assert_EqualWithin( e_real, e_formula, INV_TOL )

  END DO ChannelLoop

  ! ============================================================================
  ! 5. **** REPORT & CLEAN UP ****
  CALL test%Report()

  CALL CRTM_SfcOptics_Destroy(SfcOptics)
  CALL CRTM_GeometryInfo_Destroy(gInfo)
  Error_Status = CRTM_SpcCoeff_Destroy()

  IF ( test%n_Failed() == 0 ) THEN
    STOP 0
  ELSE
    STOP 1
  END IF

END PROGRAM test_CONST_MIXED_Polarization
