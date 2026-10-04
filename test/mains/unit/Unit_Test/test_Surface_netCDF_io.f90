!
! test_Surface_netCDF_io
!
! Round-trip unit test for the CRTM Surface netCDF file I/O added for the
! REL-3.2.0 baseline-format conversion. Builds a rank-2 Surface(L x M) array
! with distinct nonzero values in every field, writes it with NetCDF=.TRUE.,
! inquires the dimensions, reads it back, and verifies that every serialized
! field round-trips exactly (the packed-schema field map is index-sensitive,
! so each field is checked explicitly), plus an overall CRTM_Surface_Compare.
!
! STOP 0 = PASS, STOP 1 = FAIL.
!

PROGRAM test_Surface_netCDF_io

  ! -----------------
  ! Environment setup
  ! -----------------
  USE Type_Kinds        , ONLY: fp
  USE Message_Handler   , ONLY: SUCCESS, Display_Message
  USE CRTM_Surface_Define, ONLY: CRTM_Surface_type   , &
                                 CRTM_Surface_WriteFile , &
                                 CRTM_Surface_ReadFile  , &
                                 CRTM_Surface_InquireFile, &
                                 CRTM_Surface_Compare   , &
                                 CRTM_Surface_Destroy
  IMPLICIT NONE

  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_Surface_netCDF_io'
  CHARACTER(*), PARAMETER :: FILENAME     = 'test_Surface_netCDF_io.nc'
  INTEGER     , PARAMETER :: N_CHANNELS   = 3
  INTEGER     , PARAMETER :: N_PROFILES   = 2

  ! ---------
  ! Variables
  ! ---------
  TYPE(CRTM_Surface_type)              :: sfc_in(N_CHANNELS,N_PROFILES)
  TYPE(CRTM_Surface_type), ALLOCATABLE :: sfc_out(:,:)
  INTEGER :: err_stat
  INTEGER :: l, m
  INTEGER :: n_File_Channels, n_File_Profiles
  INTEGER :: n_fail

  n_fail = 0

  ! Build the input array with distinct, nonzero values
  DO m = 1, N_PROFILES
    DO l = 1, N_CHANNELS
      CALL Make_Surface( sfc_in(l,m), l, m )
    END DO
  END DO

  ! Write it out in netCDF format
  err_stat = CRTM_Surface_WriteFile( FILENAME, sfc_in, NetCDF=.TRUE., Quiet=.TRUE. )
  IF ( err_stat /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error writing netCDF Surface file', err_stat )
    STOP 1
  END IF

  ! Inquire the dimensions
  err_stat = CRTM_Surface_InquireFile( FILENAME, &
               n_Channels = n_File_Channels, &
               n_Profiles = n_File_Profiles, &
               NetCDF     = .TRUE. )
  IF ( err_stat /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error inquiring netCDF Surface file', err_stat )
    STOP 1
  END IF
  IF ( n_File_Channels /= N_CHANNELS .OR. n_File_Profiles /= N_PROFILES ) THEN
    WRITE(*,'("FAIL: inquired dims (",i0,",",i0,") /= expected (",i0,",",i0,")")') &
      n_File_Channels, n_File_Profiles, N_CHANNELS, N_PROFILES
    n_fail = n_fail + 1
  END IF

  ! Read it back
  err_stat = CRTM_Surface_ReadFile( FILENAME, sfc_out, NetCDF=.TRUE., Quiet=.TRUE. )
  IF ( err_stat /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error reading netCDF Surface file', err_stat )
    STOP 1
  END IF

  ! Check the returned shape
  IF ( .NOT. ALLOCATED(sfc_out) ) THEN
    WRITE(*,'("FAIL: sfc_out not allocated after read")')
    STOP 1
  END IF
  IF ( SIZE(sfc_out,DIM=1) /= N_CHANNELS .OR. SIZE(sfc_out,DIM=2) /= N_PROFILES ) THEN
    WRITE(*,'("FAIL: read shape (",i0,",",i0,") /= expected (",i0,",",i0,")")') &
      SIZE(sfc_out,DIM=1), SIZE(sfc_out,DIM=2), N_CHANNELS, N_PROFILES
    STOP 1
  END IF

  ! Field-by-field exact round-trip check + overall compare
  DO m = 1, N_PROFILES
    DO l = 1, N_CHANNELS
      CALL Compare_Element( sfc_in(l,m), sfc_out(l,m), l, m )
      IF ( .NOT. CRTM_Surface_Compare( sfc_in(l,m), sfc_out(l,m) ) ) THEN
        WRITE(*,'("FAIL: CRTM_Surface_Compare false at (",i0,",",i0,")")') l, m
        n_fail = n_fail + 1
      END IF
    END DO
  END DO

  ! Clean up
  CALL CRTM_Surface_Destroy( sfc_in )
  CALL CRTM_Surface_Destroy( sfc_out )

  IF ( n_fail == 0 ) THEN
    WRITE(*,'(/,"SURFACE NETCDF ROUNDTRIP PASS")')
    STOP 0
  ELSE
    WRITE(*,'(/,"SURFACE NETCDF ROUNDTRIP FAIL: ",i0," mismatch(es)")') n_fail
    STOP 1
  END IF

CONTAINS

  ! Populate every field of a Surface element with a distinct nonzero value
  ! derived from its (l,m) indices.
  SUBROUTINE Make_Surface( sfc, l, m )
    TYPE(CRTM_Surface_type), INTENT(IN OUT) :: sfc
    INTEGER,                 INTENT(IN)     :: l, m
    REAL(fp) :: r
    INTEGER  :: i
    r = REAL( 100*l + 10*m, fp )
    i = 10*l + m
    ! Coverage fractions (validity is not enforced by the netCDF path)
    sfc%Land_Coverage         = 0.001_fp * r + 0.11_fp
    sfc%Water_Coverage        = 0.001_fp * r + 0.12_fp
    sfc%Snow_Coverage         = 0.001_fp * r + 0.13_fp
    sfc%Ice_Coverage          = 0.001_fp * r + 0.14_fp
    ! Surface-type-independent
    sfc%Wind_Speed            = r + 1.0_fp
    ! Land
    sfc%Land_Temperature      = r + 2.0_fp
    sfc%Soil_Moisture_Content = r + 3.0_fp
    sfc%Canopy_Water_Content  = r + 4.0_fp
    sfc%Vegetation_Fraction   = r + 5.0_fp
    sfc%Soil_Temperature      = r + 6.0_fp
    sfc%LAI                   = r + 7.0_fp
    sfc%Land_Type             = i + 1
    sfc%Soil_Type             = i + 2
    sfc%Vegetation_Type       = i + 3
    ! Water
    sfc%Water_Temperature     = r + 8.0_fp
    sfc%Wind_Direction        = r + 9.0_fp
    sfc%Salinity              = r + 10.0_fp
    sfc%Water_Type            = i + 4
    ! Snow
    sfc%Snow_Temperature      = r + 11.0_fp
    sfc%Snow_Depth            = r + 12.0_fp
    sfc%Snow_Density          = r + 13.0_fp
    sfc%Snow_Grain_Size       = r + 14.0_fp
    sfc%Snow_Type             = i + 5
    ! Ice
    sfc%Ice_Temperature       = r + 15.0_fp
    sfc%Ice_Thickness         = r + 16.0_fp
    sfc%Ice_Density           = r + 17.0_fp
    sfc%Ice_Roughness         = r + 18.0_fp
    sfc%Ice_Type              = i + 6
  END SUBROUTINE Make_Surface

  ! Exact field-by-field comparison; increments host n_fail per mismatch.
  SUBROUTINE Compare_Element( a, b, l, m )
    TYPE(CRTM_Surface_type), INTENT(IN) :: a, b
    INTEGER,                 INTENT(IN) :: l, m
    CALL CheckR( 'Land_Coverage'        , a%Land_Coverage        , b%Land_Coverage        , l, m )
    CALL CheckR( 'Water_Coverage'       , a%Water_Coverage       , b%Water_Coverage       , l, m )
    CALL CheckR( 'Snow_Coverage'        , a%Snow_Coverage        , b%Snow_Coverage        , l, m )
    CALL CheckR( 'Ice_Coverage'         , a%Ice_Coverage         , b%Ice_Coverage         , l, m )
    CALL CheckR( 'Wind_Speed'           , a%Wind_Speed           , b%Wind_Speed           , l, m )
    CALL CheckR( 'Land_Temperature'     , a%Land_Temperature     , b%Land_Temperature     , l, m )
    CALL CheckR( 'Soil_Moisture_Content', a%Soil_Moisture_Content, b%Soil_Moisture_Content, l, m )
    CALL CheckR( 'Canopy_Water_Content' , a%Canopy_Water_Content , b%Canopy_Water_Content , l, m )
    CALL CheckR( 'Vegetation_Fraction'  , a%Vegetation_Fraction  , b%Vegetation_Fraction  , l, m )
    CALL CheckR( 'Soil_Temperature'     , a%Soil_Temperature     , b%Soil_Temperature     , l, m )
    CALL CheckR( 'LAI'                  , a%LAI                  , b%LAI                  , l, m )
    CALL CheckR( 'Water_Temperature'    , a%Water_Temperature    , b%Water_Temperature    , l, m )
    CALL CheckR( 'Wind_Direction'       , a%Wind_Direction       , b%Wind_Direction       , l, m )
    CALL CheckR( 'Salinity'             , a%Salinity             , b%Salinity             , l, m )
    CALL CheckR( 'Snow_Temperature'     , a%Snow_Temperature     , b%Snow_Temperature     , l, m )
    CALL CheckR( 'Snow_Depth'           , a%Snow_Depth           , b%Snow_Depth           , l, m )
    CALL CheckR( 'Snow_Density'         , a%Snow_Density         , b%Snow_Density         , l, m )
    CALL CheckR( 'Snow_Grain_Size'      , a%Snow_Grain_Size      , b%Snow_Grain_Size      , l, m )
    CALL CheckR( 'Ice_Temperature'      , a%Ice_Temperature      , b%Ice_Temperature      , l, m )
    CALL CheckR( 'Ice_Thickness'        , a%Ice_Thickness        , b%Ice_Thickness        , l, m )
    CALL CheckR( 'Ice_Density'          , a%Ice_Density          , b%Ice_Density          , l, m )
    CALL CheckR( 'Ice_Roughness'        , a%Ice_Roughness        , b%Ice_Roughness        , l, m )
    CALL CheckI( 'Land_Type'            , a%Land_Type            , b%Land_Type            , l, m )
    CALL CheckI( 'Soil_Type'            , a%Soil_Type            , b%Soil_Type            , l, m )
    CALL CheckI( 'Vegetation_Type'      , a%Vegetation_Type      , b%Vegetation_Type      , l, m )
    CALL CheckI( 'Water_Type'           , a%Water_Type           , b%Water_Type           , l, m )
    CALL CheckI( 'Snow_Type'            , a%Snow_Type            , b%Snow_Type            , l, m )
    CALL CheckI( 'Ice_Type'             , a%Ice_Type             , b%Ice_Type             , l, m )
  END SUBROUTINE Compare_Element

  SUBROUTINE CheckR( name, a, b, l, m )
    CHARACTER(*), INTENT(IN) :: name
    REAL(fp),     INTENT(IN) :: a, b
    INTEGER,      INTENT(IN) :: l, m
    IF ( a /= b ) THEN
      WRITE(*,'("FAIL: ",a," (",i0,",",i0,") in=",es24.16," out=",es24.16)') &
        TRIM(name), l, m, a, b
      n_fail = n_fail + 1
    END IF
  END SUBROUTINE CheckR

  SUBROUTINE CheckI( name, a, b, l, m )
    CHARACTER(*), INTENT(IN) :: name
    INTEGER,      INTENT(IN) :: a, b
    INTEGER,      INTENT(IN) :: l, m
    IF ( a /= b ) THEN
      WRITE(*,'("FAIL: ",a," (",i0,",",i0,") in=",i0," out=",i0)') &
        TRIM(name), l, m, a, b
      n_fail = n_fail + 1
    END IF
  END SUBROUTINE CheckI

END PROGRAM test_Surface_netCDF_io
