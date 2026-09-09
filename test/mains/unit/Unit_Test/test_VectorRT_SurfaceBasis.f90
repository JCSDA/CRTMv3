!
! test_VectorRT_SurfaceBasis
!
! Pins the surface (V,H) -> Stokes (I,Q) basis conversion on the vector
! radiative-transfer path (Options%n_Stokes > 1).
!
! Background
! ----------
! The microwave surface models return emissivity in the (V,H) basis:
! component 1 is the vertical emissivity eV and component 2 the horizontal
! emissivity eH. The vector solver, however, consumes the Stokes vector: the
! flattened source built by Reshape_Surf_Opt feeds CRTM_ADA / CRTM_Emission
! directly as (I,Q,U,V) per angle. CRTM's own scalar branch states the
! relationship between the two bases explicitly:
!
!     UNPOLARIZED / FIRST_STOKES  ->  ( eV + eH ) / 2      (= Stokes I)
!     SECOND_STOKES_COMPONENT     ->  ( eV - eH ) / 2      (= Stokes Q)
!     VL_POLARIZATION             ->    eV
!     HL_POLARIZATION             ->    eH
!
! so on the n_Stokes > 1 path the first two components handed to the solver
! must be (eV+eH)/2 and (eV-eH)/2, not eV and eH.
!
! What this test does
! -------------------
! It calls CRTM_Compute_SfcOptics directly for a single microwave channel over
! an ocean surface, three times on the same scene:
!
!   1. scalar (n_Stokes=1) with the channel temporarily set to VL_POLARIZATION,
!      which by the table above returns exactly eV;
!   2. scalar with the channel set to HL_POLARIZATION, returning exactly eH;
!   3. vector (n_Stokes=2) with the channel's real polarization restored.
!
! and then asserts, to machine precision,
!
!     Emissivity(:,1) == ( eV + eH ) / 2
!     Emissivity(:,2) == ( eV - eH ) / 2
!
! together with the matching reflectivity identities
!
!     R(1,1) == R(2,2) == ( rV + rH ) / 2
!     R(1,2) == R(2,1) == ( rV - rH ) / 2 .
!
! On the unconverted code the first two emissivity components come back as eV
! and eH, so the test fails by the full polarization difference (order 0.2 in
! emissivity over ocean), not by a tolerance margin.
!
! The test needs no cloud lookup table and no reference radiances: it is a
! direct algebraic statement about the surface handoff, which is what the
! self-consistency tests (TL vs finite difference, adjoint dot product,
! K vs AD) structurally cannot check.
!

PROGRAM test_VectorRT_SurfaceBasis

  ! -----------------
  ! Environment setup
  ! -----------------
  USE CRTM_Module
  USE CRTM_SpcCoeff           , ONLY: SC
  USE CRTM_SfcOptics_Define   , ONLY: CRTM_SfcOptics_type      , &
                                      CRTM_SfcOptics_Create    , &
                                      CRTM_SfcOptics_Destroy   , &
                                      CRTM_SfcOptics_Associated
  USE CRTM_SfcOptics          , ONLY: CRTM_Compute_SfcOptics, iVar_type
  USE CRTM_GeometryInfo_Define, ONLY: CRTM_GeometryInfo_type, &
                                      CRTM_GeometryInfo_SetValue
  USE CRTM_GeometryInfo       , ONLY: CRTM_GeometryInfo_Compute
  USE SensorInfo_Parameters   , ONLY: VL_POLARIZATION, HL_POLARIZATION
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_VectorRT_SurfaceBasis'
  CHARACTER(*), PARAMETER :: SENSOR       = 'amsua_n19'
  CHARACTER(*), PARAMETER :: PATH         = 'testinput/'

  ! Scene: open ocean, moderate wind, mid-latitude SST.
  REAL(fp), PARAMETER :: WIND_SPEED  = 7.0_fp      ! m/s
  REAL(fp), PARAMETER :: WATER_TEMP  = 285.0_fp    ! K
  REAL(fp), PARAMETER :: SALINITY    = 33.0_fp     ! ppmv
  REAL(fp), PARAMETER :: ZENITH      = 45.0_fp     ! deg, well away from nadir so
                                                   ! eV and eH are clearly distinct
  INTEGER , PARAMETER :: N_ANGLES    = 1
  INTEGER , PARAMETER :: CHANNEL     = 1           ! 23.8 GHz, strong V/H contrast

  ! Machine-precision assertion: the identities are exact algebra, not physics
  ! approximations, so the only slack needed is floating-point round-off.
  REAL(fp), PARAMETER :: TOL = 1.0e-12_fp

  CHARACTER(256) :: Version
  INTEGER :: Error_Status, n_Channels, i, saved_pol
  LOGICAL :: ok_eI, ok_eQ, ok_rII, ok_rIQ, all_ok
  REAL(fp) :: eV(N_ANGLES), eH(N_ANGLES)
  REAL(fp) :: rV(N_ANGLES), rH(N_ANGLES)
  REAL(fp) :: d_eI, d_eQ, d_rII, d_rIQ

  TYPE(CRTM_ChannelInfo_type)  :: ChannelInfo(1)
  TYPE(CRTM_Surface_type)      :: Sfc
  TYPE(CRTM_GeometryInfo_type) :: gInfo
  TYPE(CRTM_SfcOptics_type)    :: SfcOptics_s, SfcOptics_v
  TYPE(iVar_type)              :: iVar

  CALL CRTM_Version(Version)
  WRITE(*,'(/5x,a)') 'Vector-RT surface (V,H) -> Stokes (I,Q) basis verification'
  WRITE(*,'(5x,a/)') 'CRTM Version: '//TRIM(Version)

  ! --------------
  ! Initialize
  ! --------------
  Error_Status = CRTM_Init( (/ SENSOR /), ChannelInfo, &
                            File_Path = PATH, Quiet = .TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Init failed', FAILURE ); STOP 1
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
  IF ( n_Channels < CHANNEL ) THEN
    CALL Display_Message( PROGRAM_NAME, 'sensor has too few channels', FAILURE ); STOP 1
  END IF

  ! ------------------------
  ! Ocean surface + geometry
  ! ------------------------
  Sfc%Water_Coverage    = ONE
  Sfc%Land_Coverage     = ZERO
  Sfc%Snow_Coverage     = ZERO
  Sfc%Ice_Coverage      = ZERO
  Sfc%Water_Type        = 1                ! sea water
  Sfc%Water_Temperature = WATER_TEMP
  Sfc%Wind_Speed        = WIND_SPEED
  Sfc%Wind_Direction    = ZERO
  Sfc%Salinity          = SALINITY

  CALL CRTM_GeometryInfo_SetValue( gInfo, Sensor_Zenith_Angle = ZENITH )
  CALL CRTM_GeometryInfo_Compute( gInfo )

  ! Two SfcOptics containers on the same scene: scalar and 2-Stokes vector.
  CALL CRTM_SfcOptics_Create( SfcOptics_s, N_ANGLES, MAX_N_STOKES )
  CALL CRTM_SfcOptics_Create( SfcOptics_v, N_ANGLES, MAX_N_STOKES )
  IF ( .NOT. CRTM_SfcOptics_Associated(SfcOptics_s) .OR. &
       .NOT. CRTM_SfcOptics_Associated(SfcOptics_v) ) THEN
    CALL Display_Message( PROGRAM_NAME, 'SfcOptics_Create failed', FAILURE ); STOP 1
  END IF

  SfcOptics_s%Angle(1)      = ZENITH
  SfcOptics_s%Weight(1)     = ONE
  SfcOptics_s%Index_Sat_Ang = 1
  SfcOptics_s%n_Angles      = N_ANGLES
  SfcOptics_v = SfcOptics_s

  saved_pol = SC(1)%Polarization(CHANNEL)

  ! ------------------------------------------------------------------
  ! 1. eV: force the channel to pure vertical polarization, scalar path
  ! ------------------------------------------------------------------
  SfcOptics_s%n_Stokes             = 1
  SC(1)%Polarization(CHANNEL)      = VL_POLARIZATION
  Error_Status = CRTM_Compute_SfcOptics( Sfc, gInfo, 1, CHANNEL, SfcOptics_s, iVar )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Compute_SfcOptics (V) failed', FAILURE ); STOP 1
  END IF
  eV(1:N_ANGLES) = SfcOptics_s%Emissivity(1:N_ANGLES,1)
  DO i = 1, N_ANGLES
    rV(i) = SfcOptics_s%Reflectivity(i,1,i,1)
  END DO

  ! --------------------------------------------------------------------
  ! 2. eH: force the channel to pure horizontal polarization, scalar path
  ! --------------------------------------------------------------------
  SC(1)%Polarization(CHANNEL)      = HL_POLARIZATION
  Error_Status = CRTM_Compute_SfcOptics( Sfc, gInfo, 1, CHANNEL, SfcOptics_s, iVar )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Compute_SfcOptics (H) failed', FAILURE ); STOP 1
  END IF
  eH(1:N_ANGLES) = SfcOptics_s%Emissivity(1:N_ANGLES,1)
  DO i = 1, N_ANGLES
    rH(i) = SfcOptics_s%Reflectivity(i,1,i,1)
  END DO

  SC(1)%Polarization(CHANNEL) = saved_pol

  ! Sanity: the scene must actually be polarized, otherwise the test proves
  ! nothing (eV == eH would satisfy both the right and the wrong conversion).
  IF ( ABS(eV(1)-eH(1)) < 0.05_fp ) THEN
    WRITE(*,'(5x,a,f8.4)') 'FAIL: scene is not polarized enough, eV-eH = ', eV(1)-eH(1)
    STOP 1
  END IF

  ! --------------------------------------------
  ! 3. Vector path on the identical scene
  ! --------------------------------------------
  SfcOptics_v%n_Stokes = 2
  Error_Status = CRTM_Compute_SfcOptics( Sfc, gInfo, 1, CHANNEL, SfcOptics_v, iVar )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Compute_SfcOptics (vector) failed', FAILURE ); STOP 1
  END IF

  ! ------------------------------------------------
  ! Assertions: the (V,H) -> (I,Q) basis conversion
  ! ------------------------------------------------
  d_eI  = ZERO ; d_eQ  = ZERO
  d_rII = ZERO ; d_rIQ = ZERO
  DO i = 1, N_ANGLES
    d_eI  = MAX(d_eI,  ABS( SfcOptics_v%Emissivity(i,1) - POINT_5*(eV(i)+eH(i)) ))
    d_eQ  = MAX(d_eQ,  ABS( SfcOptics_v%Emissivity(i,2) - POINT_5*(eV(i)-eH(i)) ))
    d_rII = MAX(d_rII, ABS( SfcOptics_v%Reflectivity(i,1,i,1) - POINT_5*(rV(i)+rH(i)) ))
    d_rII = MAX(d_rII, ABS( SfcOptics_v%Reflectivity(i,2,i,2) - POINT_5*(rV(i)+rH(i)) ))
    d_rIQ = MAX(d_rIQ, ABS( SfcOptics_v%Reflectivity(i,1,i,2) - POINT_5*(rV(i)-rH(i)) ))
    d_rIQ = MAX(d_rIQ, ABS( SfcOptics_v%Reflectivity(i,2,i,1) - POINT_5*(rV(i)-rH(i)) ))
  END DO

  ok_eI  = ( d_eI  < TOL )
  ok_eQ  = ( d_eQ  < TOL )
  ok_rII = ( d_rII < TOL )
  ok_rIQ = ( d_rIQ < TOL )

  WRITE(*,'(5x,a,i0,a,f6.2,a)') 'Channel ', CHANNEL, ' at ', ZENITH, ' deg over ocean'
  WRITE(*,'(5x,a,f10.6,a,f10.6)') 'eV = ', eV(1), '   eH = ', eH(1)
  WRITE(*,'(5x,a,f10.6,a,f10.6)') 'expected  I = ', POINT_5*(eV(1)+eH(1)), &
                                  '   Q = ', POINT_5*(eV(1)-eH(1))
  WRITE(*,'(5x,a,f10.6,a,f10.6)') 'computed  I = ', SfcOptics_v%Emissivity(1,1), &
                                  '   Q = ', SfcOptics_v%Emissivity(1,2)
  WRITE(*,'(/5x,a,es12.4,a,l1)') 'emissivity I  max|diff| = ', d_eI,  '   pass = ', ok_eI
  WRITE(*,'(5x,a,es12.4,a,l1)')  'emissivity Q  max|diff| = ', d_eQ,  '   pass = ', ok_eQ
  WRITE(*,'(5x,a,es12.4,a,l1)')  'reflect. I,Q diag       = ', d_rII, '   pass = ', ok_rII
  WRITE(*,'(5x,a,es12.4,a,l1)')  'reflect. I,Q off-diag   = ', d_rIQ, '   pass = ', ok_rIQ

  all_ok = ok_eI .AND. ok_eQ .AND. ok_rII .AND. ok_rIQ

  CALL CRTM_SfcOptics_Destroy( SfcOptics_s )
  CALL CRTM_SfcOptics_Destroy( SfcOptics_v )
  Error_Status = CRTM_Destroy( ChannelInfo )

  IF ( all_ok ) THEN
    WRITE(*,'(/5x,a/)') 'PASS: surface Stokes basis conversion verified'
    STOP 0
  ELSE
    WRITE(*,'(/5x,a/)') 'FAIL: surface emissivity/reflectivity are not in the Stokes basis'
    STOP 1
  END IF

END PROGRAM test_VectorRT_SurfaceBasis
