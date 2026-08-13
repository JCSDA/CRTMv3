!
! test_TELSEM2_SnowIce
!
! Integration test for the class-consistent TELSEM2 atlas path in the MW SNOW
! and ICE surface drivers (Phase G of the TELSEM2 DA program).
!
! Contract under test: when the atlas is loaded and its classification agrees
! with the declared surface -- snow declared over a cell the January/September
! climatology classes as snow/continental ice (class2 17..22), or ice declared
! over a climatological sea-ice cell (class2 11..16) -- the drivers use the
! atlas emissivity and its Release-2 uncertainty. When the declaration
! contradicts the climatology (snow/ice declared over a vegetated cell), the
! drivers fall back to the NESDIS models bit-for-bit.
!
! Scenes (amsua_n19, clear sky):
!   1. Snow declared over Greenland (72N, Jan; class2 20)   -> atlas, std > 0
!   2. Snow declared over Argentina (30S, Sep; class2 1)    -> NESDIS fallback
!   3. Ice declared over Antarctic sea ice (72.5S, Jan; 12) -> atlas, std > 0
!   4. Ice declared over Argentina (30S, Sep; class2 1)     -> NESDIS fallback
!   5. Land declared over Greenland == snow declared over Greenland (the same
!      atlas cell serves every fraction -> identical emissivity)
!   6. Fractional land/snow (0.5/0.5) over Greenland: emissivity and std equal
!      the pure-coverage runs (identical contributions add linearly in coverage)
!
! No SensorData is attached, so the NESDIS snow baseline is the deterministic
! LandEM-with-snow-depth branch and the ice baseline the model default -- fixed
! references for the bit-identical fallback comparisons.
!
! Exit status: STOP 0 = success, STOP 1 = failure.
!
PROGRAM test_TELSEM2_SnowIce

  USE CRTM_Module
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME      = 'test_TELSEM2_SnowIce'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  CHARACTER(*), PARAMETER :: ATLAS_FILE        = 'TELSEM2.MWland.test.nc'
  CHARACTER(*), PARAMETER :: SENSOR_ID         = 'amsua_n19'

  INTEGER,  PARAMETER :: N_PROFILES  = 2   ! matches Load_Atm_Data.inc
  INTEGER,  PARAMETER :: N_LAYERS    = 92
  INTEGER,  PARAMETER :: N_ABSORBERS = 2
  REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp

  ! Greenland continental ice (class2 20/21 all year)
  REAL(fp), PARAMETER :: LAT_G = 72.0_fp,   LON_G = 320.0_fp
  ! Antarctic sea ice in January (class2 12)
  REAL(fp), PARAMETER :: LAT_I = -72.5_fp,  LON_I = 281.5_fp
  ! Vegetated cell (class2 1) -- inconsistent with declared snow/ice
  REAL(fp), PARAMETER :: LAT_A = -30.0_fp,  LON_A = 302.0_fp
  INTEGER,  PARAMETER :: MON_JAN = 1, MON_SEP = 9

  REAL(fp), PARAMETER :: TOL_ZERO = 1.0e-10_fp
  REAL(fp), PARAMETER :: MIN_DIFF = 1.0e-3_fp

  CHARACTER(256) :: Version
  INTEGER :: Error_Status, n_Channels, l
  LOGICAL :: failed

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(1)
  TYPE(CRTM_Geometry_type)    :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)
  REAL(fp), ALLOCATABLE :: e_snow_grl(:), e_snow_arg(:), e_ice_ant(:), e_ice_arg(:)
  REAL(fp), ALLOCATABLE :: e_land_grl(:), e_frac_grl(:), s_snow_grl(:), s_frac_grl(:)
  REAL(fp), ALLOCATABLE :: n_snow_grl(:), n_snow_arg(:), n_ice_ant(:), n_ice_arg(:)

  CALL CRTM_Version( Version )
  CALL Program_Message( PROGRAM_NAME, &
    'Validate the class-consistent TELSEM2 path in the MW snow/ice drivers.', &
    'CRTM Version: '//TRIM(Version) )
  failed = .FALSE.

  ! ------------------------------------------------------------------
  ! 1. Runs WITH the atlas
  ! ------------------------------------------------------------------
  Error_Status = CRTM_Init( (/SENSOR_ID/), ChannelInfo, &
                            File_Path=COEFFICIENTS_PATH, MWlandCoeff_File=ATLAS_FILE )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error initializing CRTM with atlas', FAILURE ); STOP 1
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
  ALLOCATE( RTSolution(n_Channels,N_PROFILES), &
            e_snow_grl(n_Channels), e_snow_arg(n_Channels), e_ice_ant(n_Channels), &
            e_ice_arg(n_Channels), e_land_grl(n_Channels), e_frac_grl(n_Channels), &
            s_snow_grl(n_Channels), s_frac_grl(n_Channels), &
            n_snow_grl(n_Channels), n_snow_arg(n_Channels), n_ice_ant(n_Channels), &
            n_ice_arg(n_Channels) )
  CALL CRTM_RTSolution_Create( RTSolution, N_LAYERS )
  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, 0, 0 )
  CALL Load_Atm_Data()

  CALL Set_Snow_Surface( 1.0_fp )
  CALL Run_At( LAT_G, LON_G, MON_JAN, e_snow_grl, 'snow/Greenland (atlas)' )
  s_snow_grl = RTSolution(:,1)%Surface_Emissivity_Std
  CALL Run_At( LAT_A, LON_A, MON_SEP, e_snow_arg, 'snow/Argentina (atlas)' )
  ! Inconsistent declaration -> fallback -> no atlas uncertainty
  IF ( ANY( ABS(RTSolution(:,1)%Surface_Emissivity_Std) > ZERO ) ) THEN
    CALL Display_Message( PROGRAM_NAME, &
      'Std non-zero for snow declared over a vegetated cell', FAILURE ); failed = .TRUE.
  END IF

  CALL Set_Ice_Surface( 1.0_fp )
  CALL Run_At( LAT_I, LON_I, MON_JAN, e_ice_ant, 'ice/Antarctic (atlas)' )
  IF ( ANY( RTSolution(:,1)%Surface_Emissivity_Std <= ZERO ) ) THEN
    CALL Display_Message( PROGRAM_NAME, &
      'Std not positive for ice declared over climatological sea ice', FAILURE ); failed = .TRUE.
  END IF
  CALL Run_At( LAT_A, LON_A, MON_SEP, e_ice_arg, 'ice/Argentina (atlas)' )

  CALL Set_Land_Surface( 1.0_fp )
  CALL Run_At( LAT_G, LON_G, MON_JAN, e_land_grl, 'land/Greenland (atlas)' )

  ! Fractional land/snow over the same Greenland cell
  CALL Set_Land_Surface( 0.5_fp )
  Sfc(:)%Snow_Coverage    = 0.5_fp
  Sfc(:)%Snow_Temperature = 260.0_fp
  Sfc(:)%Snow_Depth       = 100.0_fp
  CALL Run_At( LAT_G, LON_G, MON_JAN, e_frac_grl, 'land+snow/Greenland (atlas)' )
  s_frac_grl = RTSolution(:,1)%Surface_Emissivity_Std

  Error_Status = CRTM_Destroy( ChannelInfo )

  ! ------------------------------------------------------------------
  ! 2. Baseline runs WITHOUT the atlas (NESDIS snow/ice paths)
  ! ------------------------------------------------------------------
  Error_Status = CRTM_Init( (/SENSOR_ID/), ChannelInfo, File_Path=COEFFICIENTS_PATH )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error initializing CRTM (NESDIS)', FAILURE ); STOP 1
  END IF
  CALL Set_Snow_Surface( 1.0_fp )
  CALL Run_At( LAT_G, LON_G, MON_JAN, n_snow_grl, 'snow/Greenland (NESDIS)' )
  CALL Run_At( LAT_A, LON_A, MON_SEP, n_snow_arg, 'snow/Argentina (NESDIS)' )
  CALL Set_Ice_Surface( 1.0_fp )
  CALL Run_At( LAT_I, LON_I, MON_JAN, n_ice_ant, 'ice/Antarctic (NESDIS)' )
  CALL Run_At( LAT_A, LON_A, MON_SEP, n_ice_arg, 'ice/Argentina (NESDIS)' )
  Error_Status = CRTM_Destroy( ChannelInfo )

  ! ------------------------------------------------------------------
  ! 3. Assertions
  ! ------------------------------------------------------------------
  ! Consistent declarations use the atlas (emissivity changes, std positive)
  CALL Require_Different( e_snow_grl, n_snow_grl, MIN_DIFF, &
       'atlas active for consistent snow (Greenland)' )
  IF ( ANY( s_snow_grl <= ZERO ) ) THEN
    CALL Display_Message( PROGRAM_NAME, &
      'Std not positive for snow declared over Greenland', FAILURE ); failed = .TRUE.
  END IF
  CALL Require_Different( e_ice_ant, n_ice_ant, MIN_DIFF, &
       'atlas active for consistent ice (Antarctic)' )

  ! Inconsistent declarations fall back bit-for-bit
  CALL Require_Same( e_snow_arg, n_snow_arg, TOL_ZERO, &
       'NESDIS fallback for snow declared over vegetation' )
  CALL Require_Same( e_ice_arg, n_ice_arg, TOL_ZERO, &
       'NESDIS fallback for ice declared over vegetation' )

  ! One atlas cell serves every fraction
  CALL Require_Same( e_land_grl, e_snow_grl, TOL_ZERO, &
       'land-declared == snow-declared emissivity over Greenland' )
  CALL Require_Same( e_frac_grl, e_snow_grl, TOL_ZERO, &
       'fractional land+snow emissivity == pure coverage' )
  CALL Require_Same( s_frac_grl, s_snow_grl, TOL_ZERO, &
       'fractional land+snow std == pure coverage (linear in coverage)' )

  CALL CRTM_Atmosphere_Destroy( Atm )
  IF ( failed ) THEN
    CALL Display_Message( PROGRAM_NAME, 'FAILED', FAILURE ); STOP 1
  ELSE
    CALL Display_Message( PROGRAM_NAME, &
      'PASSED: class-consistent snow/ice use the atlas with uncertainty; '// &
      'inconsistent declarations fall back to NESDIS bit-for-bit; '// &
      'fractional coverages combine linearly', INFORMATION ); STOP 0
  END IF

CONTAINS

  SUBROUTINE Run_At( lat, lon, mon, emis, label )
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
  END SUBROUTINE Run_At

  SUBROUTINE Reset_Surface()
    INTEGER :: mm
    DO mm = 1, N_PROFILES
      Sfc(mm)%Land_Coverage  = 0.0_fp
      Sfc(mm)%Water_Coverage = 0.0_fp
      Sfc(mm)%Snow_Coverage  = 0.0_fp
      Sfc(mm)%Ice_Coverage   = 0.0_fp
    END DO
  END SUBROUTINE Reset_Surface

  SUBROUTINE Set_Land_Surface( coverage )
    REAL(fp), INTENT(IN) :: coverage
    INTEGER :: mm
    CALL Reset_Surface()
    DO mm = 1, N_PROFILES
      Sfc(mm)%Land_Coverage         = coverage
      Sfc(mm)%Land_Type             = 1
      Sfc(mm)%Soil_Type             = 1
      Sfc(mm)%Vegetation_Type       = 7
      Sfc(mm)%Land_Temperature      = 270.0_fp
      Sfc(mm)%Soil_Temperature      = 270.0_fp
      Sfc(mm)%Soil_Moisture_Content = 0.2_fp
      Sfc(mm)%Lai                   = 2.0_fp
      Sfc(mm)%Vegetation_Fraction   = 0.5_fp
    END DO
  END SUBROUTINE Set_Land_Surface

  SUBROUTINE Set_Snow_Surface( coverage )
    REAL(fp), INTENT(IN) :: coverage
    INTEGER :: mm
    CALL Reset_Surface()
    DO mm = 1, N_PROFILES
      Sfc(mm)%Snow_Coverage    = coverage
      Sfc(mm)%Snow_Type        = 1
      Sfc(mm)%Snow_Temperature = 260.0_fp
      Sfc(mm)%Snow_Depth       = 100.0_fp
      ! The NESDIS fallback branch reads these through the LandEM call
      Sfc(mm)%Soil_Type        = 1
      Sfc(mm)%Vegetation_Type  = 7
      Sfc(mm)%Lai              = 2.0_fp
    END DO
  END SUBROUTINE Set_Snow_Surface

  SUBROUTINE Set_Ice_Surface( coverage )
    REAL(fp), INTENT(IN) :: coverage
    INTEGER :: mm
    CALL Reset_Surface()
    DO mm = 1, N_PROFILES
      Sfc(mm)%Ice_Coverage    = coverage
      Sfc(mm)%Ice_Type        = 1
      Sfc(mm)%Ice_Temperature = 263.0_fp
      Sfc(mm)%Ice_Thickness   = 50.0_fp
    END DO
  END SUBROUTINE Set_Ice_Surface

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

  INCLUDE 'Load_Atm_Data.inc'

END PROGRAM test_TELSEM2_SnowIce
