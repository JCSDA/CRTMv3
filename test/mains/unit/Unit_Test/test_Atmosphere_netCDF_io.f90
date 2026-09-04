!
! test_Atmosphere_netCDF_io
!
! Round-trip unit test for the CRTM Atmosphere netCDF file I/O added for the
! REL-3.2.0 baseline-format conversion. Builds a rank-2 Atmosphere(L x M) array
! (with clouds AND aerosols populated) with distinct nonzero values, writes it
! with NetCDF=.TRUE., inquires the dimensions, reads it back, and verifies that
! every serialized field (including nested Cloud/Aerosol arrays) round-trips
! exactly, plus an overall CRTM_Atmosphere_Compare.
!
! STOP 0 = PASS, STOP 1 = FAIL.
!

PROGRAM test_Atmosphere_netCDF_io

  ! -----------------
  ! Environment setup
  ! -----------------
  USE Type_Kinds          , ONLY: fp
  USE Message_Handler     , ONLY: SUCCESS, Display_Message
  USE CRTM_Atmosphere_Define, ONLY: CRTM_Atmosphere_type    , &
                                    CRTM_Atmosphere_Create  , &
                                    CRTM_Atmosphere_Destroy , &
                                    CRTM_Atmosphere_Associated, &
                                    CRTM_Atmosphere_WriteFile , &
                                    CRTM_Atmosphere_ReadFile  , &
                                    CRTM_Atmosphere_InquireFile, &
                                    CRTM_Atmosphere_Compare
  IMPLICIT NONE

  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_Atmosphere_netCDF_io'
  CHARACTER(*), PARAMETER :: FILENAME     = 'test_Atmosphere_netCDF_io.nc'
  INTEGER     , PARAMETER :: N_CHANNELS   = 2
  INTEGER     , PARAMETER :: N_PROFILES   = 2
  INTEGER     , PARAMETER :: N_LAYERS     = 4
  INTEGER     , PARAMETER :: N_ABSORBERS  = 2
  INTEGER     , PARAMETER :: N_CLOUDS     = 2
  INTEGER     , PARAMETER :: N_AEROSOLS   = 2

  ! ---------
  ! Variables
  ! ---------
  TYPE(CRTM_Atmosphere_type)              :: atm_in(N_CHANNELS,N_PROFILES)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm_out(:,:)
  INTEGER :: err_stat
  INTEGER :: l, m
  INTEGER :: n_File_Channels, n_File_Profiles
  INTEGER :: n_fail

  n_fail = 0

  ! Build the input array with distinct, nonzero values
  CALL CRTM_Atmosphere_Create( atm_in, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(atm_in)) ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error creating input Atmosphere array', 1 )
    STOP 1
  END IF
  DO m = 1, N_PROFILES
    DO l = 1, N_CHANNELS
      CALL Make_Atm( atm_in(l,m), l, m )
    END DO
  END DO

  ! Write it out in netCDF format
  err_stat = CRTM_Atmosphere_WriteFile( FILENAME, atm_in, NetCDF=.TRUE., Quiet=.TRUE. )
  IF ( err_stat /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error writing netCDF Atmosphere file', err_stat )
    STOP 1
  END IF

  ! Inquire the dimensions
  err_stat = CRTM_Atmosphere_InquireFile( FILENAME, &
               n_Channels = n_File_Channels, &
               n_Profiles = n_File_Profiles, &
               NetCDF     = .TRUE. )
  IF ( err_stat /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error inquiring netCDF Atmosphere file', err_stat )
    STOP 1
  END IF
  IF ( n_File_Channels /= N_CHANNELS .OR. n_File_Profiles /= N_PROFILES ) THEN
    WRITE(*,'("FAIL: inquired dims (",i0,",",i0,") /= expected (",i0,",",i0,")")') &
      n_File_Channels, n_File_Profiles, N_CHANNELS, N_PROFILES
    n_fail = n_fail + 1
  END IF

  ! Read it back
  err_stat = CRTM_Atmosphere_ReadFile( FILENAME, atm_out, NetCDF=.TRUE., Quiet=.TRUE. )
  IF ( err_stat /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error reading netCDF Atmosphere file', err_stat )
    STOP 1
  END IF

  ! Check the returned shape
  IF ( .NOT. ALLOCATED(atm_out) ) THEN
    WRITE(*,'("FAIL: atm_out not allocated after read")')
    STOP 1
  END IF
  IF ( SIZE(atm_out,DIM=1) /= N_CHANNELS .OR. SIZE(atm_out,DIM=2) /= N_PROFILES ) THEN
    WRITE(*,'("FAIL: read shape (",i0,",",i0,") /= expected (",i0,",",i0,")")') &
      SIZE(atm_out,DIM=1), SIZE(atm_out,DIM=2), N_CHANNELS, N_PROFILES
    STOP 1
  END IF

  ! Field-by-field exact round-trip check + overall compare
  DO m = 1, N_PROFILES
    DO l = 1, N_CHANNELS
      CALL Compare_Element( atm_in(l,m), atm_out(l,m), l, m )
      IF ( .NOT. CRTM_Atmosphere_Compare( atm_in(l,m), atm_out(l,m) ) ) THEN
        WRITE(*,'("FAIL: CRTM_Atmosphere_Compare false at (",i0,",",i0,")")') l, m
        n_fail = n_fail + 1
      END IF
    END DO
  END DO

  ! Clean up
  CALL CRTM_Atmosphere_Destroy( atm_in )
  CALL CRTM_Atmosphere_Destroy( atm_out )

  IF ( n_fail == 0 ) THEN
    WRITE(*,'(/,"ATMOSPHERE NETCDF ROUNDTRIP PASS")')
    STOP 0
  ELSE
    WRITE(*,'(/,"ATMOSPHERE NETCDF ROUNDTRIP FAIL: ",i0," mismatch(es)")') n_fail
    STOP 1
  END IF

CONTAINS

  ! Populate every serialized field of an Atmosphere element with distinct
  ! nonzero values derived from its (l,m) indices. Height / n_Added_Layers /
  ! Add_Extra_Layers are intentionally left at their Create defaults (they are
  ! not serialized, mirroring the binary format).
  SUBROUTINE Make_Atm( atm, l, m )
    TYPE(CRTM_Atmosphere_type), INTENT(IN OUT) :: atm
    INTEGER,                    INTENT(IN)     :: l, m
    REAL(fp) :: r
    INTEGER  :: i, j, k, c, a
    r = REAL( 100*l + 10*m, fp )
    i = 10*l + m
    atm%Climatology = MOD(i,6) + 1
    DO j = 1, N_ABSORBERS
      atm%Absorber_ID(j)    = i + j
      atm%Absorber_Units(j) = i + 10 + j
    END DO
    DO k = 0, N_LAYERS
      atm%Level_Pressure(k) = r + 0.5_fp*REAL(k,fp) + 1.0_fp
    END DO
    DO k = 1, N_LAYERS
      atm%Pressure(k)          = r + REAL(k,fp) + 2.0_fp
      atm%Temperature(k)       = r + REAL(k,fp) + 3.0_fp
      atm%Relative_Humidity(k) = r + 0.01_fp*REAL(k,fp) + 0.1_fp
      atm%Cloud_Fraction(k)    = 0.001_fp*r + 0.01_fp*REAL(k,fp)
      DO j = 1, N_ABSORBERS
        atm%Absorber(k,j) = r + REAL(k,fp) + 0.25_fp*REAL(j,fp) + 4.0_fp
      END DO
    END DO
    DO c = 1, N_CLOUDS
      atm%Cloud(c)%Type = i + c
      DO k = 1, N_LAYERS
        atm%Cloud(c)%Effective_Radius(k)   = r + 10.0_fp*REAL(c,fp) + REAL(k,fp) + 5.0_fp
        atm%Cloud(c)%Effective_Variance(k) = r + 10.0_fp*REAL(c,fp) + REAL(k,fp) + 6.0_fp
        atm%Cloud(c)%Water_Content(k)      = r + 10.0_fp*REAL(c,fp) + REAL(k,fp) + 7.0_fp
        atm%Cloud(c)%Water_Density(k)      = r + 10.0_fp*REAL(c,fp) + REAL(k,fp) + 8.0_fp
      END DO
    END DO
    DO a = 1, N_AEROSOLS
      atm%Aerosol(a)%Type = i + 100 + a
      DO k = 1, N_LAYERS
        atm%Aerosol(a)%Effective_Radius(k)   = r + 20.0_fp*REAL(a,fp) + REAL(k,fp) + 9.0_fp
        atm%Aerosol(a)%Effective_Variance(k) = r + 20.0_fp*REAL(a,fp) + REAL(k,fp) + 10.0_fp
        atm%Aerosol(a)%Concentration(k)      = r + 20.0_fp*REAL(a,fp) + REAL(k,fp) + 11.0_fp
      END DO
    END DO
  END SUBROUTINE Make_Atm

  ! Exact field-by-field comparison; increments host n_fail per mismatch.
  SUBROUTINE Compare_Element( a, b, l, m )
    TYPE(CRTM_Atmosphere_type), INTENT(IN) :: a, b
    INTEGER,                    INTENT(IN) :: l, m
    INTEGER :: j, k, c, ae
    CALL CheckI( 'Climatology', a%Climatology, b%Climatology, l, m )
    CALL CheckI( 'n_Layers'   , a%n_Layers   , b%n_Layers   , l, m )
    CALL CheckI( 'n_Absorbers', a%n_Absorbers, b%n_Absorbers, l, m )
    CALL CheckI( 'n_Clouds'   , a%n_Clouds   , b%n_Clouds   , l, m )
    CALL CheckI( 'n_Aerosols' , a%n_Aerosols , b%n_Aerosols , l, m )
    DO j = 1, N_ABSORBERS
      CALL CheckI( 'Absorber_ID'   , a%Absorber_ID(j)   , b%Absorber_ID(j)   , l, m )
      CALL CheckI( 'Absorber_Units', a%Absorber_Units(j), b%Absorber_Units(j), l, m )
    END DO
    DO k = 0, N_LAYERS
      CALL CheckR( 'Level_Pressure', a%Level_Pressure(k), b%Level_Pressure(k), l, m )
    END DO
    DO k = 1, N_LAYERS
      CALL CheckR( 'Pressure'         , a%Pressure(k)         , b%Pressure(k)         , l, m )
      CALL CheckR( 'Temperature'      , a%Temperature(k)      , b%Temperature(k)      , l, m )
      CALL CheckR( 'Relative_Humidity', a%Relative_Humidity(k), b%Relative_Humidity(k), l, m )
      CALL CheckR( 'Cloud_Fraction'   , a%Cloud_Fraction(k)   , b%Cloud_Fraction(k)   , l, m )
      DO j = 1, N_ABSORBERS
        CALL CheckR( 'Absorber', a%Absorber(k,j), b%Absorber(k,j), l, m )
      END DO
    END DO
    DO c = 1, N_CLOUDS
      CALL CheckI( 'Cloud%Type'    , a%Cloud(c)%Type    , b%Cloud(c)%Type    , l, m )
      CALL CheckI( 'Cloud%n_Layers', a%Cloud(c)%n_Layers, b%Cloud(c)%n_Layers, l, m )
      DO k = 1, N_LAYERS
        CALL CheckR( 'Cloud%Effective_Radius'  , a%Cloud(c)%Effective_Radius(k)  , b%Cloud(c)%Effective_Radius(k)  , l, m )
        CALL CheckR( 'Cloud%Effective_Variance', a%Cloud(c)%Effective_Variance(k), b%Cloud(c)%Effective_Variance(k), l, m )
        CALL CheckR( 'Cloud%Water_Content'     , a%Cloud(c)%Water_Content(k)     , b%Cloud(c)%Water_Content(k)     , l, m )
        CALL CheckR( 'Cloud%Water_Density'     , a%Cloud(c)%Water_Density(k)     , b%Cloud(c)%Water_Density(k)     , l, m )
      END DO
    END DO
    DO ae = 1, N_AEROSOLS
      CALL CheckI( 'Aerosol%Type'    , a%Aerosol(ae)%Type    , b%Aerosol(ae)%Type    , l, m )
      CALL CheckI( 'Aerosol%n_Layers', a%Aerosol(ae)%n_Layers, b%Aerosol(ae)%n_Layers, l, m )
      DO k = 1, N_LAYERS
        CALL CheckR( 'Aerosol%Effective_Radius'  , a%Aerosol(ae)%Effective_Radius(k)  , b%Aerosol(ae)%Effective_Radius(k)  , l, m )
        CALL CheckR( 'Aerosol%Effective_Variance', a%Aerosol(ae)%Effective_Variance(k), b%Aerosol(ae)%Effective_Variance(k), l, m )
        CALL CheckR( 'Aerosol%Concentration'     , a%Aerosol(ae)%Concentration(k)     , b%Aerosol(ae)%Concentration(k)     , l, m )
      END DO
    END DO
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

END PROGRAM test_Atmosphere_netCDF_io
