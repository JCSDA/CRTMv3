!
! test_TELSEM2_MWland
!
! Integration test for the TELSEM2 microwave land surface emissivity atlas.
!
! For a pure-land MW scene (amsua_n19) at a known land cell it checks two things:
!
!   1. When the TELSEM2 atlas is loaded it is actually used: the surface
!      emissivity differs from the NESDIS_LandEM fallback (run with the atlas
!      absent) and stays physically reasonable.
!
!   2. The atlas path is independent of the surface state. The atlas depends only
!      on lat/lon/month/frequency/angle, so the K-matrix and tangent-linear
!      sensitivities of Tb to LAI, Vegetation_Fraction and Soil_Moisture_Content
!      must be exactly zero (in contrast to the NESDIS path exercised by
!      test_Land_Jacobian, where they are non-zero).
!
! The atlas is loaded explicitly via the MWlandCoeff_File argument pointing at a
! test-only staged copy, so the default drop-in behaviour of the other land
! tests (which use NESDIS_LandEM) is left unchanged.
!
! Exit status: STOP 0 = success, STOP 1 = failure.
!

PROGRAM test_TELSEM2_MWland

  USE CRTM_Module
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME      = 'test_TELSEM2_MWland'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  CHARACTER(*), PARAMETER :: ATLAS_FILE        = 'TELSEM2.MWland.test.nc'  ! test-only staged name
  CHARACTER(*), PARAMETER :: SENSOR_ID         = 'amsua_n19'

  INTEGER,  PARAMETER :: N_PROFILES  = 2   ! matches Load_Atm_Data.inc
  INTEGER,  PARAMETER :: N_LAYERS    = 92
  INTEGER,  PARAMETER :: N_ABSORBERS = 2
  INTEGER,  PARAMETER :: N_CLOUDS    = 0
  INTEGER,  PARAMETER :: N_AEROSOLS  = 0
  INTEGER,  PARAMETER :: N_SENSORS   = 1

  REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp
  ! Land cell with TELSEM2 climatology (northern Argentina), September
  REAL(fp), PARAMETER :: LATITUDE  = -30.0_fp
  REAL(fp), PARAMETER :: LONGITUDE = 302.0_fp
  INTEGER,  PARAMETER :: MONTH     = 9
  ! Base land state
  REAL(fp), PARAMETER :: LAI0 = 2.0_fp, VEG0 = 0.5_fp, SMC0 = 0.2_fp
  ! Zero-sensitivity tolerance for the atlas path (should be exact zero)
  REAL(fp), PARAMETER :: TOL_ZERO = 1.0e-10_fp
  ! Minimum emissivity change vs NESDIS to prove the atlas is active
  REAL(fp), PARAMETER :: MIN_DIFF = 1.0e-3_fp

  CHARACTER(256) :: Message, Version
  INTEGER :: Error_Status, Alloc_Status, n_Channels, l, m
  LOGICAL :: failed
  REAL(fp) :: maxdiff, maxjac

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)    :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm(N_PROFILES), Atm_TL(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc(N_PROFILES), Sfc_TL(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:), RTSolution_TL(:,:), RTSolution_K(:,:)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atmosphere_K(:,:)
  TYPE(CRTM_Surface_type),    ALLOCATABLE :: Surface_K(:,:)
  REAL(fp), ALLOCATABLE :: emis_atlas(:), emis_nesdis(:)

  CALL CRTM_Version( Version )
  CALL Program_Message( PROGRAM_NAME, &
    'Validate the TELSEM2 MW land emissivity atlas integration.', &
    'CRTM Version: '//TRIM(Version) )

  failed = .FALSE.

  ! ------------------------------------------------------------------
  ! 1. Initialise WITH the TELSEM2 atlas
  ! ------------------------------------------------------------------
  Error_Status = CRTM_Init( (/SENSOR_ID/), ChannelInfo, &
                            File_Path=COEFFICIENTS_PATH, MWlandCoeff_File=ATLAS_FILE )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error initializing CRTM with atlas', FAILURE ); STOP 1
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))

  ALLOCATE( RTSolution(n_Channels,N_PROFILES), RTSolution_TL(n_Channels,N_PROFILES), &
            RTSolution_K(n_Channels,N_PROFILES), Atmosphere_K(n_Channels,N_PROFILES), &
            Surface_K(n_Channels,N_PROFILES), emis_atlas(n_Channels), emis_nesdis(n_Channels), &
            STAT = Alloc_Status )
  IF ( Alloc_Status /= 0 ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error allocating arrays', FAILURE ); STOP 1
  END IF
  CALL CRTM_Atmosphere_Create( Atm,    N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atm_TL, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atmosphere_K, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )

  CALL Load_Atm_Data()
  CALL Load_Land_Surface()
  CALL CRTM_Geometry_SetValue( Geometry, &
                               Sensor_Zenith_Angle = ZENITH_ANGLE, &
                               Sensor_Scan_Angle   = SCAN_ANGLE, &
                               Latitude  = LATITUDE, &
                               Longitude = LONGITUDE, &
                               Month     = MONTH )

  ! 1a. Forward -> atlas emissivity
  Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error in CRTM Forward (atlas)', FAILURE ); STOP 1
  END IF
  emis_atlas = RTSolution(:,1)%Surface_Emissivity
  DO l = 1, n_Channels
    IF ( emis_atlas(l) < 0.4_fp .OR. emis_atlas(l) > 1.0_fp ) THEN
      WRITE(Message,'("TELSEM2 emissivity out of range at channel ",i0,": ",es13.5)') l, emis_atlas(l)
      CALL Display_Message( PROGRAM_NAME, TRIM(Message), FAILURE ); failed = .TRUE.
    END IF
  END DO

  ! 1b. K-matrix -> surface Jacobians must be zero on the atlas path
  CALL CRTM_Atmosphere_Zero( Atmosphere_K )
  CALL CRTM_Surface_Zero( Surface_K )
  RTSolution_K%Radiance               = ZERO
  RTSolution_K%Brightness_Temperature = ONE
  Error_Status = CRTM_K_Matrix( Atm, Sfc, RTSolution_K, Geometry, ChannelInfo, &
                                Atmosphere_K, Surface_K, RTSolution )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error in CRTM K_Matrix (atlas)', FAILURE ); STOP 1
  END IF
  maxjac = ZERO
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      maxjac = MAX( maxjac, ABS(Surface_K(l,m)%Lai), &
                            ABS(Surface_K(l,m)%Vegetation_Fraction), &
                            ABS(Surface_K(l,m)%Soil_Moisture_Content) )
    END DO
  END DO
  WRITE(*,'(5x,"Max |K-matrix surface Jacobian| on atlas path = ",es12.4)') maxjac
  IF ( maxjac > TOL_ZERO ) THEN
    CALL Display_Message( PROGRAM_NAME, &
      'K-matrix surface Jacobian non-zero on atlas path', FAILURE ); failed = .TRUE.
  END IF

  ! 1c. Tangent-linear -> Tb sensitivity must be zero for each surface direction
  maxjac = ZERO
  CALL TL_Direction( 'LAI' )
  CALL TL_Direction( 'VEG' )
  CALL TL_Direction( 'SMC' )
  WRITE(*,'(5x,"Max |tangent-linear dTb| on atlas path        = ",es12.4)') maxjac
  IF ( maxjac > TOL_ZERO ) THEN
    CALL Display_Message( PROGRAM_NAME, &
      'Tangent-linear Tb non-zero on atlas path', FAILURE ); failed = .TRUE.
  END IF

  Error_Status = CRTM_Destroy( ChannelInfo )

  ! ------------------------------------------------------------------
  ! 2. Initialise WITHOUT the atlas -> NESDIS_LandEM fallback
  ! ------------------------------------------------------------------
  Error_Status = CRTM_Init( (/SENSOR_ID/), ChannelInfo, File_Path=COEFFICIENTS_PATH )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error initializing CRTM (NESDIS)', FAILURE ); STOP 1
  END IF
  Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error in CRTM Forward (NESDIS)', FAILURE ); STOP 1
  END IF
  emis_nesdis = RTSolution(:,1)%Surface_Emissivity
  Error_Status = CRTM_Destroy( ChannelInfo )

  ! 2a. The atlas must change the emissivity (i.e. it was actually used)
  maxdiff = ZERO
  WRITE(*,'(/5x,a)') 'Chan   NESDIS-emis   TELSEM2-emis   diff'
  DO l = 1, n_Channels
    WRITE(*,'(5x,i4,3f14.6)') RTSolution(l,1)%Sensor_Channel, emis_nesdis(l), emis_atlas(l), &
                              emis_atlas(l)-emis_nesdis(l)
    maxdiff = MAX( maxdiff, ABS(emis_atlas(l)-emis_nesdis(l)) )
  END DO
  WRITE(*,'(5x,"Max |TELSEM2 - NESDIS| = ",es12.4)') maxdiff
  IF ( maxdiff < MIN_DIFF ) THEN
    CALL Display_Message( PROGRAM_NAME, &
      'Atlas did not change the emissivity (not active?)', FAILURE ); failed = .TRUE.
  END IF

  ! 3. Clean up
  CALL CRTM_Atmosphere_Destroy( Atm )
  CALL CRTM_Atmosphere_Destroy( Atm_TL )
  CALL CRTM_Atmosphere_Destroy( Atmosphere_K )
  DEALLOCATE( RTSolution, RTSolution_TL, RTSolution_K, Atmosphere_K, Surface_K, &
              emis_atlas, emis_nesdis )

  IF ( failed ) THEN
    CALL Display_Message( PROGRAM_NAME, 'FAILED', FAILURE ); STOP 1
  ELSE
    CALL Display_Message( PROGRAM_NAME, &
      'PASSED: atlas active, surface Jacobians zero on atlas path', INFORMATION ); STOP 0
  END IF

CONTAINS

  ! Run the tangent-linear model with a unit perturbation in one surface variable
  ! and accumulate the largest |dTb| (which must be zero on the atlas path).
  SUBROUTINE TL_Direction( which )
    CHARACTER(*), INTENT(IN) :: which
    INTEGER :: stat
    CALL CRTM_Atmosphere_Zero( Atm_TL )
    CALL CRTM_Surface_Zero( Sfc_TL )
    SELECT CASE ( which )
      CASE ('LAI'); Sfc_TL%Lai                   = ONE
      CASE ('VEG'); Sfc_TL%Vegetation_Fraction   = ONE
      CASE ('SMC'); Sfc_TL%Soil_Moisture_Content = ONE
    END SELECT
    stat = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, ChannelInfo, &
                                RTSolution, RTSolution_TL )
    IF ( stat /= SUCCESS ) THEN
      CALL Display_Message( PROGRAM_NAME, 'Error in CRTM Tangent_Linear ('//which//')', FAILURE ); STOP 1
    END IF
    maxjac = MAX( maxjac, MAXVAL(ABS(RTSolution_TL%Brightness_Temperature)) )
  END SUBROUTINE TL_Direction

  ! Pure-land surface for all profiles
  SUBROUTINE Load_Land_Surface()
    INTEGER :: mm
    DO mm = 1, N_PROFILES
      Sfc(mm)%Land_Coverage         = 1.0_fp
      Sfc(mm)%Water_Coverage        = 0.0_fp
      Sfc(mm)%Snow_Coverage         = 0.0_fp
      Sfc(mm)%Ice_Coverage          = 0.0_fp
      Sfc(mm)%Land_Type             = 1
      Sfc(mm)%Soil_Type             = 1
      Sfc(mm)%Vegetation_Type       = 7
      Sfc(mm)%Land_Temperature      = 290.0_fp
      Sfc(mm)%Soil_Temperature      = 290.0_fp
      Sfc(mm)%Soil_Moisture_Content = SMC0
      Sfc(mm)%Lai                   = LAI0
      Sfc(mm)%Vegetation_Fraction   = VEG0
    END DO
  END SUBROUTINE Load_Land_Surface

  INCLUDE 'Load_Atm_Data.inc'

END PROGRAM test_TELSEM2_MWland
