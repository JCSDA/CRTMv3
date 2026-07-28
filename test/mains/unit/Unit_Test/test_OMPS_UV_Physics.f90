!
! test_OMPS_UV_Physics
!
! Baseline-independent physics verification of the four per-platform OMPS UV
! products (u.omps-np_n20, u.omps-np_n21, u.omps-tc_n20, u.omps-tc_n21),
! all loaded in a single CRTM_Init call and run as one multi-sensor forward.
!
! On the ECMWF84 ocean column (clear sky, daytime solar geometry, scene NO2
! from the GEOS-CF climatology used by test_UV_NO2_TLAD) the test asserts,
! with tolerances 2-4x wider than the values measured at product acceptance:
!
!   1. Radiances positive and finite on every channel.
!   2. BUV spectral shape: the nadir profilers' normalized radiance
!      (radiance / band solar irradiance) collapses from 310 nm into the
!      Hartley band; the total-column mappers peak near 340 nm.
!   3. Cross-platform consistency: NOAA-20 vs NOAA-21 normalized radiance at
!      matched wavelengths (NP pair and TC pair) agrees to rms < 1%,
!      max < 3%. A wavelength-registration or channel-numbering error of
!      even a fraction of the 0.42 nm channel spacing breaks this through
!      the Fraunhofer structure.
!   4. Dichroic-range consistency: NP vs TC on the same platform over
!      302-310 nm agrees to < 2.5% (the operational NM/NP standard is ~2%
!      on the real instruments).
!   5. Ozone weighting functions of both profilers peak monotonically
!      deeper with increasing wavelength (260 -> 310 nm anchors).
!   6. Scene-NO2 doubling lowers radiance on every channel, with the peak
!      response on the long-wavelength side for the mappers.
!   7. An ozone increase lowers radiance strongly near 307 nm.
!   8. Adjoint dot-product closure over T + H2O + O3 + NO2, and K == AD on
!      the most NO2-sensitive channel of each sensor.
!
! The four coefficient pairs are pre-release: the test is registered only
! when they are present (symlinked) in the source testinput directory.
!
! Exit: STOP 0 if every check passes, STOP 1 otherwise.
!
! CREATION HISTORY:
!       Written by:     Benjamin Johnson, 28-Jul-2026
!                       Promoted from the OMPS product-acceptance driver
!                       (coeff_consistency_2026-07-26/OMPS_PRODUCT_VERIFICATION.md).
!
PROGRAM test_OMPS_UV_Physics

  USE CRTM_Module
  USE SpcCoeff_Define   , ONLY: SpcCoeff_type, SpcCoeff_Destroy
  USE SpcCoeff_netCDF_IO, ONLY: SpcCoeff_netCDF_ReadFile
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_OMPS_UV_Physics'
  CHARACTER(*), PARAMETER :: PATH = './testinput/'
  INTEGER,      PARAMETER :: N_SENSORS = 4
  CHARACTER(13), PARAMETER :: SENSORS(N_SENSORS) = &
    (/ 'u.omps-np_n20', 'u.omps-np_n21', 'u.omps-tc_n20', 'u.omps-tc_n21' /)
  INTEGER, PARAMETER :: S_NP20 = 1, S_NP21 = 2, S_TC20 = 3, S_TC21 = 4
  INTEGER, PARAMETER :: MAX_CH = 198

  INTEGER,  PARAMETER :: N_PROFILES  = 2
  INTEGER,  PARAMETER :: N_LAYERS    = 100
  INTEGER,  PARAMETER :: N_ABSORBERS = 7
  INTEGER,  PARAMETER :: N_CLOUDS    = 0
  INTEGER,  PARAMETER :: N_AEROSOLS  = 0
  REAL(fp), PARAMETER :: ZENITH      = 45.0_fp
  REAL(fp), PARAMETER :: SOLAR_ZEN   = 30.0_fp
  INTEGER,  PARAMETER :: IDX_H2O = 1, IDX_O3 = 3, IDX_NO2 = 7

  REAL(fp), PARAMETER :: TOL_ADJ = 1.0e-10_fp   ! solar-RT roundoff bound
  REAL(fp), PARAMETER :: TOL_K   = 1.0e-9_fp

  CHARACTER(256) :: Version
  INTEGER  :: Error_Status, Allocate_Status
  INTEGER  :: i, l, m, ii, ntot, kpeak, l0
  INTEGER  :: n_per(N_SENSORS), off(N_SENSORS)
  REAL(fp) :: lam(MAX_CH,N_SENSORS), esun(MAX_CH,N_SENSORS)
  REAL(fp) :: no2_clim(N_LAYERS)
  REAL(fp) :: LHS, RHS, dy, rel, mx, rms, pk_prev, pk
  REAL(fp), ALLOCATABLE :: R0(:,:), R_no2(:), R_o3(:)
  REAL(fp) :: wf(N_LAYERS)
  LOGICAL  :: all_ok
  TYPE(SpcCoeff_type) :: sc

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)    :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm_TL(N_PROFILES), Atm_AD(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc_TL(N_PROFILES), Sfc_AD(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:), RTSolution_pert(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_TL(:,:), RTSolution_AD(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_K(:,:)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atm_K(:,:)
  TYPE(CRTM_Surface_type),    ALLOCATABLE :: Sfc_K(:,:)

  all_ok = .TRUE.

  CALL CRTM_Version(Version)
  WRITE(*,'(/5x,a)') 'OMPS UV physics verification (4 per-platform products)'
  WRITE(*,'(5x,a/)') 'CRTM Version: '//TRIM(Version)

  Error_Status = CRTM_Init( SENSORS, ChannelInfo, File_Path=PATH, Quiet=.TRUE. )
  CALL judge( 'CRTM_Init (all four sensors, one call)', Error_Status == SUCCESS )
  IF ( Error_Status /= SUCCESS ) STOP 1

  ntot = 0
  DO i = 1, N_SENSORS
    n_per(i) = CRTM_ChannelInfo_n_Channels(ChannelInfo(i))
    off(i)   = ntot
    ntot     = ntot + n_per(i)
  END DO
  CALL judge( 'channel counts 151/158/196/198', &
              ALL( n_per == (/151, 158, 196, 198/) ) )

  ! Band centers and solar irradiance from the loaded SpcCoeff files.
  lam = ZERO ; esun = ZERO
  DO i = 1, N_SENSORS
    Error_Status = SpcCoeff_netCDF_ReadFile( PATH//TRIM(SENSORS(i))//'.SpcCoeff.nc', sc, Quiet=.TRUE. )
    IF ( Error_Status /= SUCCESS .OR. SIZE(sc%Wavenumber) /= n_per(i) ) THEN
      CALL judge( 'SpcCoeff re-read for '//TRIM(SENSORS(i)), .FALSE. )
    ELSE
      lam(1:n_per(i),i)  = 1.0e7_fp / sc%Wavenumber
      esun(1:n_per(i),i) = sc%Solar_Irradiance
    END IF
    CALL SpcCoeff_Destroy( sc )
  END DO

  ALLOCATE( RTSolution(ntot,N_PROFILES), RTSolution_pert(ntot,N_PROFILES), &
            RTSolution_TL(ntot,N_PROFILES), RTSolution_AD(ntot,N_PROFILES), &
            RTSolution_K(ntot,N_PROFILES), &
            Atm_K(ntot,N_PROFILES), Sfc_K(ntot,N_PROFILES), &
            R0(ntot,N_PROFILES), R_no2(ntot), R_o3(ntot), STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN; WRITE(*,*) 'Alloc error'; STOP 1; END IF

  CALL CRTM_RTSolution_Create( RTSolution,      N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_pert, N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_TL,   N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_AD,   N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_K,    N_LAYERS )

  CALL CRTM_Atmosphere_Create( Atm,    N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atm_TL, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atm_AD, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atm_K,  N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )

  CALL Load_ECMWF84_Atm_Data()
  CALL Set_NO2_Climatology()
  DO m = 1, N_PROFILES
    Atm(m)%Absorber_Id(IDX_NO2)    = NO2_ID
    Atm(m)%Absorber_Units(IDX_NO2) = VOLUME_MIXING_RATIO_UNITS
    Atm(m)%Absorber(:,IDX_NO2)     = no2_clim
  END DO
  Atm(2)%Absorber(:,IDX_NO2) = 1.3_fp * Atm(2)%Absorber(:,IDX_NO2)

  DO m = 1, N_PROFILES
    Atm_TL(m)%Climatology = Atm(m)%Climatology
    Atm_TL(m)%Absorber_ID = Atm(m)%Absorber_ID ; Atm_TL(m)%Absorber_Units = Atm(m)%Absorber_Units
    Atm_AD(m)%Climatology = Atm(m)%Climatology
    Atm_AD(m)%Absorber_ID = Atm(m)%Absorber_ID ; Atm_AD(m)%Absorber_Units = Atm(m)%Absorber_Units
    DO l = 1, ntot
      Atm_K(l,m)%Climatology = Atm(m)%Climatology
      Atm_K(l,m)%Absorber_ID = Atm(m)%Absorber_ID ; Atm_K(l,m)%Absorber_Units = Atm(m)%Absorber_Units
    END DO
  END DO

  DO m = 1, N_PROFILES
    Sfc(m)%Water_Coverage    = ONE
    Sfc(m)%Water_Type        = 1
    Sfc(m)%Water_Temperature = 290.0_fp
    Sfc(m)%Wind_Speed        = 6.0_fp
    Sfc(m)%Salinity          = 33.0_fp
    CALL CRTM_Geometry_SetValue( Geometry(m), Sensor_Zenith_Angle = ZENITH, &
                                 Source_Zenith_Angle = SOLAR_ZEN )
  END DO

  ! ---- base forward, all four sensors in one call ----
  Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution )
  CALL judge( 'multi-sensor forward', Error_Status == SUCCESS )
  IF ( Error_Status /= SUCCESS ) STOP 1
  DO m = 1, N_PROFILES
    DO l = 1, ntot
      R0(l,m) = RTSolution(l,m)%Radiance
    END DO
  END DO
  CALL judge( 'all radiances positive and finite', &
              ALL( R0 > ZERO ) .AND. ALL( ABS(R0) < HUGE(ONE) ) )

  ! ---- BUV spectral shape (normalized radiance = R / E_sun) ----
  ! Profilers: Hartley-band collapse relative to 310 nm.
  CALL judge( 'np_n20 Hartley cutoff (NR 250 nm < 0.1 x NR 310 nm)', &
              nr_at(S_NP20, 250.0_fp) < 0.1_fp * nr_at(S_NP20, 310.0_fp) )
  CALL judge( 'np_n21 Hartley cutoff (NR 250 nm < 0.1 x NR 310 nm)', &
              nr_at(S_NP21, 250.0_fp) < 0.1_fp * nr_at(S_NP21, 310.0_fp) )
  ! Mappers: Huggins rise to a broad maximum near 340 nm.
  CALL judge( 'tc_n20 band shape (NR 340 > 3 x NR 305; NR 340 > NR 380)', &
              nr_at(S_TC20, 340.0_fp) > 3.0_fp * nr_at(S_TC20, 305.0_fp) .AND. &
              nr_at(S_TC20, 340.0_fp) > nr_at(S_TC20, 380.0_fp) )
  CALL judge( 'tc_n21 band shape (NR 340 > 3 x NR 305; NR 340 > NR 380)', &
              nr_at(S_TC21, 340.0_fp) > 3.0_fp * nr_at(S_TC21, 305.0_fp) .AND. &
              nr_at(S_TC21, 340.0_fp) > nr_at(S_TC21, 380.0_fp) )

  ! ---- cross-platform consistency at matched wavelengths ----
  CALL cross_platform( S_NP20, S_NP21, mx, rms )
  CALL judge( 'NP cross-platform NR (rms < 1%, max < 3%)', &
              rms < 0.01_fp .AND. mx < 0.03_fp )
  WRITE(*,'(9x,"NP  n20 vs n21: rms ",f6.3,"%  max ",f6.3,"%")') 100.0_fp*rms, 100.0_fp*mx
  CALL cross_platform( S_TC20, S_TC21, mx, rms )
  CALL judge( 'TC cross-platform NR (rms < 1%, max < 3%)', &
              rms < 0.01_fp .AND. mx < 0.03_fp )
  WRITE(*,'(9x,"TC  n20 vs n21: rms ",f6.3,"%  max ",f6.3,"%")') 100.0_fp*rms, 100.0_fp*mx

  ! ---- dichroic overlap, same platform (302-310 nm) ----
  CALL dichroic( S_NP20, S_TC20, mx )
  CALL judge( 'n20 NP vs TC dichroic overlap (max < 2.5%)', mx < 0.025_fp )
  WRITE(*,'(9x,"n20 dichroic max ",f6.3,"%")') 100.0_fp*mx
  CALL dichroic( S_NP21, S_TC21, mx )
  CALL judge( 'n21 NP vs TC dichroic overlap (max < 2.5%)', mx < 0.025_fp )
  WRITE(*,'(9x,"n21 dichroic max ",f6.3,"%")') 100.0_fp*mx

  ! ---- NO2 doubled ----
  DO m = 1, N_PROFILES
    Atm(m)%Absorber(:,IDX_NO2) = 2.0_fp * Atm(m)%Absorber(:,IDX_NO2)
  END DO
  Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution_pert )
  IF ( Error_Status /= SUCCESS ) THEN; WRITE(*,*) 'forward fail'; STOP 1; END IF
  DO l = 1, ntot
    R_no2(l) = RTSolution_pert(l,1)%Radiance
  END DO
  DO m = 1, N_PROFILES
    Atm(m)%Absorber(:,IDX_NO2) = 0.5_fp * Atm(m)%Absorber(:,IDX_NO2)
  END DO
  CALL judge( 'NO2 x2 lowers radiance on every channel', &
              ALL( (R_no2 - R0(:,1)) / R0(:,1) <= 1.0e-9_fp ) )
  DO i = 3, 4   ! mappers: response window and long-wavelength peak
    mx = ZERO ; l0 = 1
    DO l = off(i)+1, off(i)+n_per(i)
      rel = ABS( (R_no2(l) - R0(l,1)) / R0(l,1) )
      IF ( rel > mx ) THEN; mx = rel; l0 = l - off(i); END IF
    END DO
    CALL judge( TRIM(SENSORS(i))//' NO2 response in [0.1%,2%], peak beyond 350 nm', &
                mx > 0.001_fp .AND. mx < 0.02_fp .AND. lam(l0,i) > 350.0_fp )
  END DO

  ! ---- O3 + 5% ----
  DO m = 1, N_PROFILES
    Atm(m)%Absorber(:,IDX_O3) = 1.05_fp * Atm(m)%Absorber(:,IDX_O3)
  END DO
  Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution_pert )
  IF ( Error_Status /= SUCCESS ) THEN; WRITE(*,*) 'forward fail'; STOP 1; END IF
  DO l = 1, ntot
    R_o3(l) = RTSolution_pert(l,1)%Radiance
  END DO
  DO m = 1, N_PROFILES
    Atm(m)%Absorber(:,IDX_O3) = Atm(m)%Absorber(:,IDX_O3) / 1.05_fp
  END DO
  DO i = 1, N_SENSORS
    mx = ZERO
    DO l = off(i)+1, off(i)+n_per(i)
      mx = MIN( mx, (R_o3(l) - R0(l,1)) / R0(l,1) )
    END DO
    CALL judge( TRIM(SENSORS(i))//' O3 +5% strongest response in [-20%,-5%]', &
                mx < -0.05_fp .AND. mx > -0.20_fp )
  END DO
  CALL judge( 'O3 +5% never raises radiance above +0.1%', &
              ALL( (R_o3 - R0(:,1)) / R0(:,1) < 1.0e-3_fp ) )

  ! ---- K-Matrix, all channels ----
  CALL CRTM_Atmosphere_Zero( Atm_K ) ; CALL CRTM_Surface_Zero( Sfc_K )
  CALL CRTM_RTSolution_Zero( RTSolution_K )
  DO m = 1, N_PROFILES
    DO l = 1, ntot
      RTSolution_K(l,m)%Radiance = ONE
    END DO
  END DO
  Error_Status = CRTM_K_Matrix( Atm, Sfc, RTSolution_K, Geometry, ChannelInfo, &
                                Atm_K, Sfc_K, RTSolution )
  CALL judge( 'multi-sensor K-Matrix', Error_Status == SUCCESS )
  IF ( Error_Status /= SUCCESS ) STOP 1

  ! Profiler ozone weighting functions descend with wavelength.
  DO i = 1, 2
    pk_prev = ZERO
    ok_block: BLOCK
      LOGICAL :: mono
      REAL(fp), PARAMETER :: ANCHORS(6) = &
        (/ 260.0_fp, 280.0_fp, 290.0_fp, 300.0_fp, 305.0_fp, 310.0_fp /)
      mono = .TRUE.
      DO ii = 1, SIZE(ANCHORS)
        l0 = nearest_ch( i, ANCHORS(ii) )
        wf = Atm_K(off(i)+l0,1)%Absorber(:,IDX_O3) * Atm(1)%Absorber(:,IDX_O3)
        kpeak = MAXLOC( ABS(wf), DIM=1 )
        pk = Atm(1)%Pressure(kpeak)
        IF ( ii > 1 .AND. pk < pk_prev ) mono = .FALSE.
        IF ( ii > 2 .AND. pk <= pk_prev ) mono = .FALSE.   ! strict from 280 on
        pk_prev = pk
      END DO
      CALL judge( TRIM(SENSORS(i))//' O3 weighting-function peaks descend 260->310 nm', mono )
    END BLOCK ok_block
  END DO

  ! ---- adjoint dot-product closure (T + H2O + O3 + NO2) ----
  CALL CRTM_Atmosphere_Zero( Atm_TL ) ; CALL CRTM_Surface_Zero( Sfc_TL )
  DO m = 1, N_PROFILES
    DO ii = 1, N_LAYERS
      Atm_TL(m)%Temperature(ii)      = 0.5_fp * SIN( 0.7_fp*REAL(ii,fp) + 1.3_fp*REAL(m,fp) )
      Atm_TL(m)%Absorber(ii,IDX_H2O) = 0.05_fp * Atm(m)%Absorber(ii,IDX_H2O) &
                                        * COS( 0.9_fp*REAL(ii,fp) + 0.4_fp*REAL(m,fp) )
      Atm_TL(m)%Absorber(ii,IDX_O3)  = 0.05_fp * Atm(m)%Absorber(ii,IDX_O3) &
                                        * SIN( 0.5_fp*REAL(ii,fp) + 0.9_fp*REAL(m,fp) )
      Atm_TL(m)%Absorber(ii,IDX_NO2) = 0.05_fp * Atm(m)%Absorber(ii,IDX_NO2) &
                                        * SIN( 1.1_fp*REAL(ii,fp) + 0.8_fp*REAL(m,fp) )
    END DO
  END DO
  Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, ChannelInfo, &
                                      RTSolution, RTSolution_TL )
  IF ( Error_Status /= SUCCESS ) THEN; WRITE(*,*) 'TL fail'; STOP 1; END IF
  LHS = ZERO
  CALL CRTM_RTSolution_Zero( RTSolution_AD )
  DO m = 1, N_PROFILES
    DO l = 1, ntot
      dy = RTSolution_TL(l,m)%Radiance
      LHS = LHS + dy*dy
      RTSolution_AD(l,m)%Radiance = dy
    END DO
  END DO
  CALL CRTM_Atmosphere_Zero( Atm_AD ) ; CALL CRTM_Surface_Zero( Sfc_AD )
  Error_Status = CRTM_Adjoint( Atm, Sfc, RTSolution_AD, Geometry, ChannelInfo, &
                               Atm_AD, Sfc_AD, RTSolution )
  IF ( Error_Status /= SUCCESS ) THEN; WRITE(*,*) 'AD fail'; STOP 1; END IF
  RHS = ZERO
  DO m = 1, N_PROFILES
    DO ii = 1, N_LAYERS
      RHS = RHS + Atm_TL(m)%Temperature(ii)      * Atm_AD(m)%Temperature(ii)
      RHS = RHS + Atm_TL(m)%Absorber(ii,IDX_H2O) * Atm_AD(m)%Absorber(ii,IDX_H2O)
      RHS = RHS + Atm_TL(m)%Absorber(ii,IDX_O3)  * Atm_AD(m)%Absorber(ii,IDX_O3)
      RHS = RHS + Atm_TL(m)%Absorber(ii,IDX_NO2) * Atm_AD(m)%Absorber(ii,IDX_NO2)
    END DO
  END DO
  rel = ABS(LHS-RHS) / MAX(ABS(LHS), TINY(ONE))
  CALL judge( 'adjoint dot-product closure (T+H2O+O3+NO2)', rel < TOL_ADJ )
  WRITE(*,'(9x,"<dy,dy>=",es16.9,"  <x,gx>=",es16.9,"  rel=",es10.3)') LHS, RHS, rel

  ! ---- K vs AD on each sensor's most NO2-sensitive channel ----
  DO i = 1, N_SENSORS
    mx = -ONE ; l0 = off(i) + 1
    DO l = off(i)+1, off(i)+n_per(i)
      rel = ABS( R_no2(l) - R0(l,1) )
      IF ( rel > mx ) THEN; mx = rel; l0 = l; END IF
    END DO
    CALL CRTM_Atmosphere_Zero( Atm_AD ) ; CALL CRTM_Surface_Zero( Sfc_AD )
    CALL CRTM_RTSolution_Zero( RTSolution_AD )
    RTSolution_AD(l0,1)%Radiance = ONE
    Error_Status = CRTM_Adjoint( Atm, Sfc, RTSolution_AD, Geometry, ChannelInfo, &
                                 Atm_AD, Sfc_AD, RTSolution )
    IF ( Error_Status /= SUCCESS ) THEN; WRITE(*,*) 'AD fail'; STOP 1; END IF
    mx = MAX( MAXVAL( ABS( Atm_K(l0,1)%Temperature - Atm_AD(1)%Temperature ) ), &
              MAXVAL( ABS( Atm_K(l0,1)%Absorber(:,IDX_O3)  - Atm_AD(1)%Absorber(:,IDX_O3) ) ), &
              MAXVAL( ABS( Atm_K(l0,1)%Absorber(:,IDX_NO2) - Atm_AD(1)%Absorber(:,IDX_NO2) ) ) )
    rel = mx / MAX( MAXVAL(ABS(Atm_K(l0,1)%Temperature)), &
                    MAXVAL(ABS(Atm_K(l0,1)%Absorber(:,IDX_O3))), &
                    MAXVAL(ABS(Atm_K(l0,1)%Absorber(:,IDX_NO2))), TINY(ONE) )
    CALL judge( TRIM(SENSORS(i))//' K == AD on most NO2-sensitive channel', rel < TOL_K )
  END DO

  Error_Status = CRTM_Destroy( ChannelInfo )

  WRITE(*,'(/5x,a)') '====================================================='
  IF ( all_ok ) THEN
    WRITE(*,'(5x,a)') 'ALL CHECKS PASSED'
    STOP 0
  ELSE
    WRITE(*,'(5x,a)') 'CHECKS FAILED'
    STOP 1
  END IF

CONTAINS

  SUBROUTINE judge( name, ok )
    CHARACTER(*), INTENT(IN) :: name
    LOGICAL,      INTENT(IN) :: ok
    WRITE(*,'(5x,"[",a,"] ",a)') MERGE('PASS','FAIL',ok), name
    IF ( .NOT. ok ) all_ok = .FALSE.
  END SUBROUTINE judge

  ! Channel of sensor i nearest to wavelength target (nm).
  INTEGER FUNCTION nearest_ch( i, target )
    INTEGER,  INTENT(IN) :: i
    REAL(fp), INTENT(IN) :: target
    nearest_ch = MINLOC( ABS( lam(1:n_per(i),i) - target ), DIM=1 )
  END FUNCTION nearest_ch

  ! Normalized radiance of sensor i at the channel nearest to target (nm).
  REAL(fp) FUNCTION nr_at( i, target )
    INTEGER,  INTENT(IN) :: i
    REAL(fp), INTENT(IN) :: target
    INTEGER :: k
    k = nearest_ch( i, target )
    nr_at = R0(off(i)+k,1) / esun(k,i)
  END FUNCTION nr_at

  ! Normalized radiance of sensor i linearly interpolated to wavelength wl.
  REAL(fp) FUNCTION nr_interp( i, wl )
    INTEGER,  INTENT(IN) :: i
    REAL(fp), INTENT(IN) :: wl
    INTEGER  :: k
    REAL(fp) :: w, nr_lo, nr_hi
    k = 1
    DO WHILE ( k < n_per(i)-1 .AND. lam(k+1,i) < wl )
      k = k + 1
    END DO
    w = ( wl - lam(k,i) ) / ( lam(k+1,i) - lam(k,i) )
    nr_lo = R0(off(i)+k,  1) / esun(k,  i)
    nr_hi = R0(off(i)+k+1,1) / esun(k+1,i)
    nr_interp = (ONE - w)*nr_lo + w*nr_hi
  END FUNCTION nr_interp

  ! Relative normalized-radiance difference of sensors ia vs ib over their
  ! common wavelength range, evaluated on ia's grid with ib interpolated.
  SUBROUTINE cross_platform( ia, ib, max_rel, rms_rel )
    INTEGER,  INTENT(IN)  :: ia, ib
    REAL(fp), INTENT(OUT) :: max_rel, rms_rel
    REAL(fp) :: lo, hi, nra, nrb, r, s
    INTEGER  :: k, n
    lo = MAX( lam(1,ia), lam(1,ib) )
    hi = MIN( lam(n_per(ia),ia), lam(n_per(ib),ib) )
    max_rel = ZERO ; s = ZERO ; n = 0
    DO k = 1, n_per(ia)
      IF ( lam(k,ia) < lo .OR. lam(k,ia) > hi ) CYCLE
      nra = R0(off(ia)+k,1) / esun(k,ia)
      nrb = nr_interp( ib, lam(k,ia) )
      r   = (nra - nrb) / nrb
      max_rel = MAX( max_rel, ABS(r) )
      s = s + r*r ; n = n + 1
    END DO
    rms_rel = SQRT( s / REAL(MAX(n,1),fp) )
  END SUBROUTINE cross_platform

  ! Max relative NP-vs-TC normalized-radiance difference over 302-310 nm.
  SUBROUTINE dichroic( i_np, i_tc, max_rel )
    INTEGER,  INTENT(IN)  :: i_np, i_tc
    REAL(fp), INTENT(OUT) :: max_rel
    REAL(fp) :: nra, nrb, r
    INTEGER  :: k
    max_rel = ZERO
    DO k = 1, n_per(i_np)
      IF ( lam(k,i_np) < 302.0_fp .OR. lam(k,i_np) > 310.0_fp ) CYCLE
      nra = R0(off(i_np)+k,1) / esun(k,i_np)
      nrb = nr_interp( i_tc, lam(k,i_np) )
      r   = (nra - nrb) / nrb
      max_rel = MAX( max_rel, ABS(r) )
    END DO
  END SUBROUTINE dichroic

  ! Daytime NO2 reference (GEOS-CF mean over the TEMPO O-B validation domain)
  ! on the ECMWF84 100-layer grid; ppmv. Same profile family as the group-8
  ! TauCoeff Ref_Absorber (see test_UV_NO2_TLAD).
  SUBROUTINE Set_NO2_Climatology()
    no2_clim = (/ &
      7.019022e-09_fp, 9.475711e-09_fp, 6.922833e-08_fp, 5.354402e-07_fp, 2.481293e-06_fp,  &
      8.153599e-06_fp, 2.185519e-05_fp, 5.164744e-05_fp, 1.178604e-04_fp, 2.551657e-04_fp,  &
      5.161727e-04_fp, 9.991564e-04_fp, 1.839313e-03_fp, 2.953800e-03_fp, 4.047471e-03_fp,  &
      4.924563e-03_fp, 5.398932e-03_fp, 5.551236e-03_fp, 5.462271e-03_fp, 5.151474e-03_fp,  &
      4.615022e-03_fp, 4.000497e-03_fp, 3.426450e-03_fp, 2.866598e-03_fp, 2.270770e-03_fp,  &
      1.883071e-03_fp, 1.616375e-03_fp, 1.534265e-03_fp, 1.441712e-03_fp, 1.342972e-03_fp,  &
      1.174026e-03_fp, 1.006428e-03_fp, 8.623019e-04_fp, 7.306015e-04_fp, 6.184407e-04_fp,  &
      4.961050e-04_fp, 3.598377e-04_fp, 2.714093e-04_fp, 2.404335e-04_fp, 2.145286e-04_fp,  &
      1.967208e-04_fp, 1.784723e-04_fp, 1.568922e-04_fp, 1.358346e-04_fp, 1.082782e-04_fp,  &
      8.063631e-05_fp, 5.949495e-05_fp, 4.488859e-05_fp, 3.060448e-05_fp, 2.444995e-05_fp,  &
      1.899578e-05_fp, 1.452061e-05_fp, 1.313692e-05_fp, 1.178142e-05_fp, 1.138935e-05_fp,  &
      1.277731e-05_fp, 1.413816e-05_fp, 1.577888e-05_fp, 1.818921e-05_fp, 2.055447e-05_fp,  &
      2.287432e-05_fp, 2.285868e-05_fp, 2.284186e-05_fp, 2.282535e-05_fp, 2.269303e-05_fp,  &
      2.247449e-05_fp, 2.225972e-05_fp, 2.208659e-05_fp, 2.194654e-05_fp, 2.184646e-05_fp,  &
      2.192552e-05_fp, 2.200356e-05_fp, 2.202075e-05_fp, 2.203063e-05_fp, 2.206344e-05_fp,  &
      2.210567e-05_fp, 2.228723e-05_fp, 2.251860e-05_fp, 2.303798e-05_fp, 2.363930e-05_fp,  &
      2.498871e-05_fp, 2.630772e-05_fp, 2.630885e-05_fp, 2.561118e-05_fp, 2.342883e-05_fp,  &
      2.000029e-05_fp, 1.652792e-05_fp, 1.373170e-05_fp, 1.182672e-05_fp, 1.070237e-05_fp,  &
      1.041749e-05_fp, 1.086241e-05_fp, 1.220802e-05_fp, 1.701331e-05_fp, 2.478315e-05_fp,  &
      3.423267e-05_fp, 4.403548e-05_fp, 4.735015e-05_fp, 4.735015e-05_fp, 4.735015e-05_fp /)
  END SUBROUTINE Set_NO2_Climatology

  INCLUDE 'Load_ECMWF84_Atm_Data.inc'

END PROGRAM test_OMPS_UV_Physics
