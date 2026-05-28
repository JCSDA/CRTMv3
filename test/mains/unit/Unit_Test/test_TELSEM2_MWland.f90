!
! test_TELSEM2_MWland
!
! Integration test for the TELSEM2 microwave land surface emissivity atlas.
!
! For pure-land MW scenes (amsua_n19) it checks:
!
!   1. The atlas is actually used when loaded: surface emissivity differs from
!      the NESDIS_LandEM fallback (run with the atlas absent) and stays physical.
!
!   2. The atlas path is independent of the surface state: the K-matrix and
!      tangent-linear sensitivities of Tb to LAI, Vegetation_Fraction and
!      Soil_Moisture_Content are exactly zero (in contrast to test_Land_Jacobian,
!      where the NESDIS path gives non-zero values).
!
!   3. The geometry inputs drive the lookup:
!        - a second land cell gives a different emissivity (spatial dependence);
!        - a different month at the same cell gives a different emissivity
!          (seasonal dependence);
!        - an ocean point (no land climatology) falls back to NESDIS_LandEM
!          through the full CRTM_Forward path, bit-identical to a NESDIS-only run.
!
! The atlas is loaded explicitly via the MWlandCoeff_File argument pointing at a
! test-only staged copy, so the default drop-in behaviour of the other land
! tests is left unchanged.
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

  ! Land cell A with TELSEM2 climatology (northern Argentina)
  REAL(fp), PARAMETER :: LAT_A = -30.0_fp, LON_A = 302.0_fp
  ! Land cell B (central North America) -- expect a different emissivity
  REAL(fp), PARAMETER :: LAT_B =  35.0_fp, LON_B = 260.0_fp
  ! Ocean point (equatorial mid-Pacific) -- no land climatology -> NESDIS fallback
  REAL(fp), PARAMETER :: LAT_O =   0.0_fp, LON_O = 200.0_fp
  INTEGER,  PARAMETER :: MON_SEP = 9, MON_JAN = 1

  ! Base land state
  REAL(fp), PARAMETER :: LAI0 = 2.0_fp, VEG0 = 0.5_fp, SMC0 = 0.2_fp
  ! Tolerances
  REAL(fp), PARAMETER :: TOL_ZERO = 1.0e-10_fp  ! "exact zero" for atlas-path Jacobians / fallback equality
  REAL(fp), PARAMETER :: MIN_DIFF = 1.0e-3_fp   ! atlas-vs-NESDIS and spatial difference
  REAL(fp), PARAMETER :: MIN_SEAS = 1.0e-4_fp   ! seasonal (month) difference

  CHARACTER(256) :: Message, Version
  INTEGER :: Error_Status, Alloc_Status, n_Channels, l, m
  LOGICAL :: failed
  REAL(fp) :: maxjac

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)    :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm(N_PROFILES), Atm_TL(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc(N_PROFILES), Sfc_TL(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:), RTSolution_TL(:,:), RTSolution_K(:,:)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atmosphere_K(:,:)
  TYPE(CRTM_Surface_type),    ALLOCATABLE :: Surface_K(:,:)
  ! Surface emissivity (profile 1) for each scenario
  REAL(fp), ALLOCATABLE :: emis_A(:), emis_B(:), emis_A_jan(:), emis_O_atlas(:)
  REAL(fp), ALLOCATABLE :: emis_nesdis_A(:), emis_nesdis_O(:)

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
            Surface_K(n_Channels,N_PROFILES), &
            emis_A(n_Channels), emis_B(n_Channels), emis_A_jan(n_Channels), &
            emis_O_atlas(n_Channels), emis_nesdis_A(n_Channels), emis_nesdis_O(n_Channels), &
            STAT = Alloc_Status )
  IF ( Alloc_Status /= 0 ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error allocating arrays', FAILURE ); STOP 1
  END IF
  CALL CRTM_Atmosphere_Create( Atm,    N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atm_TL, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atmosphere_K, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )

  CALL Load_Atm_Data()
  CALL Load_Land_Surface()

  ! 1a. Forward at cell A (Sep) -> atlas emissivity
  CALL Run_Forward_At( LAT_A, LON_A, MON_SEP, emis_A, 'A/Sep' )
  DO l = 1, n_Channels
    IF ( emis_A(l) < 0.4_fp .OR. emis_A(l) > 1.0_fp ) THEN
      WRITE(Message,'("TELSEM2 emissivity out of range at channel ",i0,": ",es13.5)') l, emis_A(l)
      CALL Display_Message( PROGRAM_NAME, TRIM(Message), FAILURE ); failed = .TRUE.
    END IF
  END DO

  ! 1b. K-matrix at cell A -> surface Jacobians must be zero on the atlas path
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
    CALL Display_Message( PROGRAM_NAME, 'K-matrix surface Jacobian non-zero on atlas path', FAILURE )
    failed = .TRUE.
  END IF

  ! 1c. Tangent-linear at cell A -> Tb sensitivity must be zero for each direction
  maxjac = ZERO
  CALL TL_Direction( 'LAI' )
  CALL TL_Direction( 'VEG' )
  CALL TL_Direction( 'SMC' )
  WRITE(*,'(5x,"Max |tangent-linear dTb| on atlas path        = ",es12.4)') maxjac
  IF ( maxjac > TOL_ZERO ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Tangent-linear Tb non-zero on atlas path', FAILURE )
    failed = .TRUE.
  END IF

  ! 1d. Spatial dependence: a different land cell must give a different emissivity
  CALL Run_Forward_At( LAT_B, LON_B, MON_SEP, emis_B, 'B/Sep' )
  CALL Require_Different( emis_A, emis_B, MIN_DIFF, 'spatial (cell A vs cell B)' )

  ! 1e. Seasonal dependence: a different month at cell A must differ
  CALL Run_Forward_At( LAT_A, LON_A, MON_JAN, emis_A_jan, 'A/Jan' )
  CALL Require_Different( emis_A, emis_A_jan, MIN_SEAS, 'seasonal (cell A Sep vs Jan)' )

  ! 1f. Ocean point: atlas has no data -> falls back through CRTM_Forward
  CALL Run_Forward_At( LAT_O, LON_O, MON_SEP, emis_O_atlas, 'ocean (atlas loaded)' )

  Error_Status = CRTM_Destroy( ChannelInfo )

  ! ------------------------------------------------------------------
  ! 2. Initialise WITHOUT the atlas -> NESDIS_LandEM fallback
  ! ------------------------------------------------------------------
  Error_Status = CRTM_Init( (/SENSOR_ID/), ChannelInfo, File_Path=COEFFICIENTS_PATH )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error initializing CRTM (NESDIS)', FAILURE ); STOP 1
  END IF
  CALL Run_Forward_At( LAT_A, LON_A, MON_SEP, emis_nesdis_A, 'A/Sep (NESDIS)' )
  CALL Run_Forward_At( LAT_O, LON_O, MON_SEP, emis_nesdis_O, 'ocean (NESDIS)' )
  Error_Status = CRTM_Destroy( ChannelInfo )

  ! 2a. Atlas must change the emissivity at cell A (i.e. it was actually used)
  CALL Require_Different( emis_A, emis_nesdis_A, MIN_DIFF, 'atlas active (cell A: TELSEM2 vs NESDIS)' )

  ! 2b. Ocean fallback: atlas-loaded run must equal the NESDIS-only run (bit-identical)
  CALL Require_Same( emis_O_atlas, emis_nesdis_O, TOL_ZERO, 'ocean fallback (atlas run == NESDIS run)' )

  ! ------------------------------------------------------------------
  ! 3. Report and clean up
  ! ------------------------------------------------------------------
  WRITE(*,'(/5x,a)') 'Chan   NESDIS(A)   TELSEM2(A)  TELSEM2(B)  TELSEM2(A,Jan)   ocean(atlas/NESDIS)'
  DO l = 1, n_Channels
    WRITE(*,'(5x,i4,4f12.6,2x,2f12.6)') RTSolution(l,1)%Sensor_Channel, &
          emis_nesdis_A(l), emis_A(l), emis_B(l), emis_A_jan(l), emis_O_atlas(l), emis_nesdis_O(l)
  END DO

  CALL CRTM_Atmosphere_Destroy( Atm )
  CALL CRTM_Atmosphere_Destroy( Atm_TL )
  CALL CRTM_Atmosphere_Destroy( Atmosphere_K )
  DEALLOCATE( RTSolution, RTSolution_TL, RTSolution_K, Atmosphere_K, Surface_K, &
              emis_A, emis_B, emis_A_jan, emis_O_atlas, emis_nesdis_A, emis_nesdis_O )

  IF ( failed ) THEN
    CALL Display_Message( PROGRAM_NAME, 'FAILED', FAILURE ); STOP 1
  ELSE
    CALL Display_Message( PROGRAM_NAME, &
      'PASSED: atlas active; spatial/seasonal dependence; ocean fallback; zero surface Jacobians', &
      INFORMATION ); STOP 0
  END IF

CONTAINS

  ! Set the geometry to (lat,lon,month) for all profiles and run the forward
  ! model, returning profile 1's per-channel surface emissivity.
  SUBROUTINE Run_Forward_At( lat, lon, mon, emis, label )
    REAL(fp),     INTENT(IN)  :: lat, lon
    INTEGER,      INTENT(IN)  :: mon
    REAL(fp),     INTENT(OUT) :: emis(:)
    CHARACTER(*), INTENT(IN)  :: label
    INTEGER :: stat
    CALL CRTM_Geometry_SetValue( Geometry, &
                                 Sensor_Zenith_Angle = ZENITH_ANGLE, &
                                 Sensor_Scan_Angle   = SCAN_ANGLE, &
                                 Latitude  = lat, Longitude = lon, Month = mon )
    stat = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution )
    IF ( stat /= SUCCESS ) THEN
      CALL Display_Message( PROGRAM_NAME, 'Error in CRTM Forward ('//label//')', FAILURE ); STOP 1
    END IF
    emis = RTSolution(:,1)%Surface_Emissivity
  END SUBROUTINE Run_Forward_At

  ! Tangent-linear with a unit perturbation in one surface variable; accumulate
  ! the largest |dTb| (must be zero on the atlas path). Uses the current Geometry.
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

  ! Require that two emissivity spectra differ by at least 'tol' on some channel.
  SUBROUTINE Require_Different( a, b, tol, label )
    REAL(fp),     INTENT(IN) :: a(:), b(:), tol
    CHARACTER(*), INTENT(IN) :: label
    REAL(fp) :: d
    d = MAXVAL(ABS(a-b))
    WRITE(*,'(5x,"Max |diff| for ",a," = ",es12.4)') label, d
    IF ( d < tol ) THEN
      CALL Display_Message( PROGRAM_NAME, 'Expected a difference but found none: '//label, FAILURE )
      failed = .TRUE.
    END IF
  END SUBROUTINE Require_Different

  ! Require that two emissivity spectra are equal to within 'tol' on all channels.
  SUBROUTINE Require_Same( a, b, tol, label )
    REAL(fp),     INTENT(IN) :: a(:), b(:), tol
    CHARACTER(*), INTENT(IN) :: label
    REAL(fp) :: d
    d = MAXVAL(ABS(a-b))
    WRITE(*,'(5x,"Max |diff| for ",a," = ",es12.4)') label, d
    IF ( d > tol ) THEN
      CALL Display_Message( PROGRAM_NAME, 'Expected equality but found a difference: '//label, FAILURE )
      failed = .TRUE.
    END IF
  END SUBROUTINE Require_Same

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
