!
! test_VectorRT_Unsupported
!
! Asserts that combinations the vector (n_Stokes > 1) path cannot honour are
! refused rather than silently substituted.
!
! Why this exists
! ---------------
! The polarimetric path accumulated its defects by being quietly wrong rather
! than loudly broken, so every combination known to be unsupported should fail
! where a user can see it.
!
! The case covered here is the radiative transfer algorithm selector. In
! CRTM_RTSolution the n_Stokes > 1 branch is taken before RT_Algorithm_Id is
! ever consulted, so a caller who asks for SOI receives ADA instead. SOI has no
! vector solver, and the two algorithms do not agree, so the caller is handed
! another algorithm's answer under the name of the one they asked for. That is
! a silent substitution, not a graceful fallback.
!
! What this test does
! -------------------
! It runs one microwave scene twice:
!
!   1. RT_ADA at n_Stokes = 2, which must SUCCEED. This is the control: without
!      it the test would also pass against a build that rejected every vector
!      run for some unrelated reason.
!   2. RT_SOI at n_Stokes = 2, which must FAIL.
!
! Against the unguarded code the second call returns SUCCESS, having quietly
! run ADA, so the test fails.
!
! No cloud lookup table is needed: the guard is reached before any solver runs,
! and the ADA control is a clear-sky vector run, which is now a supported
! configuration in its own right.
!

PROGRAM test_VectorRT_Unsupported

  USE CRTM_Module
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_VectorRT_Unsupported'
  CHARACTER(*), PARAMETER :: PATH   = './testinput/'
  CHARACTER(*), PARAMETER :: SENSOR = 'amsua_n19'

  ! Load_ECMWF84_Atm_Data fills atm(1) AND atm(2), so two profiles are
  ! mandatory; asking for one writes out of bounds and segfaults.
  INTEGER,  PARAMETER :: N_PROFILES  = 2
  INTEGER,  PARAMETER :: N_LAYERS    = 100
  INTEGER,  PARAMETER :: N_ABSORBERS = 6
  INTEGER,  PARAMETER :: N_CLOUDS    = 0
  INTEGER,  PARAMETER :: N_AEROSOLS  = 0
  REAL(fp), PARAMETER :: ZENITH      = 53.0_fp

  CHARACTER(256) :: Version
  INTEGER :: Error_Status, Allocate_Status, n_Channels, m
  INTEGER :: stat_ada, stat_soi
  LOGICAL :: ok_ada, ok_soi, all_ok

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(1)
  TYPE(CRTM_Geometry_type)    :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc(N_PROFILES)
  TYPE(CRTM_Options_type)     :: Options(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)

  CALL CRTM_Version(Version)
  WRITE(*,'(/5x,a)') 'Vector-RT unsupported-combination refusal'
  WRITE(*,'(5x,a/)') 'CRTM Version: '//TRIM(Version)

  Error_Status = CRTM_Init( (/ SENSOR /), ChannelInfo, &
                            File_Path = PATH, Quiet = .TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Init failed', FAILURE ); STOP 1
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))

  ALLOCATE( RTSolution(n_Channels,N_PROFILES), STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN; WRITE(*,*) 'Alloc error'; STOP 1; END IF
  CALL CRTM_RTSolution_Create( RTSolution, N_LAYERS )
  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Atmosphere_Create failed', FAILURE ); STOP 1
  END IF

  CALL Load_ECMWF84_Atm_Data()          ! fills Atm(1) and Atm(2)

  DO m = 1, N_PROFILES
    Sfc(m)%Water_Coverage    = ONE
    Sfc(m)%Water_Type        = 1
    Sfc(m)%Water_Temperature = 290.0_fp
    Sfc(m)%Wind_Speed        = 6.0_fp
    Sfc(m)%Salinity          = 33.0_fp
    CALL CRTM_Geometry_SetValue( Geometry(m), Sensor_Zenith_Angle = ZENITH )
  END DO

  ! ------------------------------------------------------------------
  ! 1. Control: ADA at n_Stokes = 2 must succeed
  ! ------------------------------------------------------------------
  DO m = 1, N_PROFILES
    Options(m)%n_Stokes        = 2
    Options(m)%RT_Algorithm_Id = RT_ADA
  END DO
  stat_ada = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution, Options=Options )
  ok_ada = ( stat_ada == SUCCESS )

  ! ------------------------------------------------------------------
  ! 2. SOI at n_Stokes = 2 must be refused
  ! ------------------------------------------------------------------
  DO m = 1, N_PROFILES
    Options(m)%n_Stokes        = 2
    Options(m)%RT_Algorithm_Id = RT_SOI
  END DO
  WRITE(*,'(5x,a)') 'The following error message is EXPECTED:'
  stat_soi = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution, Options=Options )
  ok_soi = ( stat_soi /= SUCCESS )

  WRITE(*,'(/5x,a,i0,a,l1)') 'RT_ADA + n_Stokes=2 status = ', stat_ada, &
        '   (expect SUCCESS)   pass = ', ok_ada
  WRITE(*,'(5x,a,i0,a,l1)')  'RT_SOI + n_Stokes=2 status = ', stat_soi, &
        '   (expect FAILURE)   pass = ', ok_soi

  all_ok = ok_ada .AND. ok_soi

  Error_Status = CRTM_Destroy( ChannelInfo )

  IF ( all_ok ) THEN
    WRITE(*,'(/5x,a/)') 'PASS: unsupported vector combinations are refused, supported ones run'
    STOP 0
  ELSE
    IF ( .NOT. ok_soi ) THEN
      WRITE(*,'(/5x,a)') 'FAIL: SOI with n_Stokes>1 returned SUCCESS. It silently ran ADA'
      WRITE(*,'(5x,a)')  '      and reported another algorithm''s answer as SOI''s.'
    ELSE
      WRITE(*,'(/5x,a)') 'FAIL: the supported ADA vector control did not run.'
    END IF
    WRITE(*,'(a)') ''
    STOP 1
  END IF

CONTAINS

  INCLUDE 'Load_ECMWF84_Atm_Data.inc'

END PROGRAM test_VectorRT_Unsupported
