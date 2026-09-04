!
! test_VectorRT_TLADK
!
! Baseline-independent TL/AD/K correctness check for the vector-RT
! (n_Stokes > 1) cloud-scattering path.
!
! The n_Stokes > 1 ADA branch is reachable only with a >= 6-phase-element
! cloud LUT (the experimental 'CRTM-Exp' scheme; stock LUTs carry a single
! phase element and are hard-rejected by the forward guard), so none of the
! standard regression tests exercise it. This test initializes mwr_aws with
! Cloud_Model='CRTM-Exp' + CloudCoeff_Exp_Full6.nc, runs an overcast snow
! column over ocean at Options%n_Stokes=2 (overcast so the solver is
! isolated from the fractional-cloud combine, whose n_Stokes>1 adjoint is a
! known deferred item), and verifies:
!   1. TL vs central finite-difference of the forward model, for BOTH
!      Stokes components:
!        - d Stokes(1:2) / d Cloud%Water_Content  (phase-matrix /
!          Normalize_Phase chain, incl. the polarized-block D2 mirror)
!        - d Stokes(1:2) / d Temperature          (AMOM thermal-source
!          intensity-slot guard + Kirchhoff sum)
!   2. Adjoint dot-product over the full Stokes vector
!      <dy,dy> == <x, AD(dy)>   with x spanning Temperature AND
!      Water_Content on every layer/profile  (AD = TL^T ?)
!   3. K-Matrix vs Adjoint Jacobian equality (Temperature and
!      Water_Content columns)
! and runs the same scene at n_Stokes=1 (scalar control) to validate the
! harness against the long-verified scalar path.
!
! Exit: STOP 0 if every check passes, STOP 1 otherwise.
!
! CREATION HISTORY:
!       Written by:     Benjamin Johnson, 11-Jun-2026
!                       Setup adapted from test_CloudCoeff_Exp_Forward;
!                       verification machinery from test_Downwelling_TLADK.
!
PROGRAM test_VectorRT_TLADK

  USE CRTM_Module
  USE CRTM_MWwaterCoeff, ONLY: CRTM_MWwaterCoeff_Load_FASTEM
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_VectorRT_TLADK'
  CHARACTER(*), PARAMETER :: PATH   = './testinput/'
  CHARACTER(*), PARAMETER :: SENSOR = 'mwr_aws'
  CHARACTER(*), PARAMETER :: LUT    = 'CloudCoeff_Exp_Full6.nc'

  ! Profile / column setup (ECMWF84 ocean column, as in the Exp forward test)
  INTEGER,  PARAMETER :: N_PROFILES  = 2
  INTEGER,  PARAMETER :: N_LAYERS    = 100
  INTEGER,  PARAMETER :: N_ABSORBERS = 6
  INTEGER,  PARAMETER :: N_CLOUDS    = 1
  INTEGER,  PARAMETER :: N_AEROSOLS  = 0
  REAL(fp), PARAMETER :: ZENITH      = 53.0_fp
  INTEGER,  PARAMETER :: KC1 = 78, KC2 = 86       ! cloud vertical band (layers)
  INTEGER,  PARAMETER :: KP  = 82                 ! perturbed layer (mid-band)
  REAL(fp), PARAMETER :: REFF_S = 500.0_fp        ! snow effective radius (microns)
  REAL(fp), PARAMETER :: WC_S   = 1.0_fp          ! kg/m^2 per layer
  ! Small enough that every layer's single-scatter albedo falls below CRTM's
  ! scattering threshold, so the solve leaves ADA for the emission path.
  REAL(fp), PARAMETER :: WC_CLEAR = 1.0e-8_fp
  REAL(fp), PARAMETER :: WIND_DIR   = 100.0_fp   ! relative azimuth = 60 deg
  REAL(fp), PARAMETER :: SENSOR_AZI =  40.0_fp

  ! The adjoint dot-product tolerance is deliberately tight (the correct code
  ! achieves ~1e-15): a one-sided TL/AD inconsistency in the phase-normalization
  ! polarized blocks shows up at ~5e-11, which 1e-9 would let through.
  REAL(fp), PARAMETER :: TOL_FD  = 1.0e-3_fp      ! TL vs finite difference
  REAL(fp), PARAMETER :: TOL_ADJ = 1.0e-12_fp     ! adjoint dot-product
  REAL(fp), PARAMETER :: TOL_K   = 1.0e-9_fp      ! K vs AD

  ! Perturbation-variable selectors for the FD check
  INTEGER, PARAMETER :: VAR_WC = 1, VAR_T = 2, VAR_WSP = 3, VAR_WDIR = 4

  CHARACTER(256) :: Version
  INTEGER :: Error_Status, Allocate_Status, n_Channels
  INTEGER :: l, m
  LOGICAL :: ok_s1_fd, ok_s1_adj, ok_s1_k
  LOGICAL :: ok_v_fd_wc1, ok_v_fd_wc2, ok_v_fd_t1, ok_v_fd_t2, ok_v_adj, ok_v_k
  LOGICAL :: ok_c_fd_t1, ok_c_fd_t2, ok_c_adj, ok_c_k
  LOGICAL :: ok_f_fd, ok_f_adj, ok_f_k
  LOGICAL :: ok_4_fd1, ok_4_fd2, ok_4_fd3, ok_4_fd4, ok_4_adj, ok_4_k
  LOGICAL :: ok_4_fdR, ok_4_adjR
  LOGICAL :: ok_4f_fdR, ok_4f_adjR
  LOGICAL :: ok_4_wdirU, ok_4_wdirV, ok_4_wspQ

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(1)
  TYPE(CRTM_Geometry_type)    :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm_TL(N_PROFILES), Atm_AD(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc_TL(N_PROFILES), Sfc_AD(N_PROFILES)
  TYPE(CRTM_Options_type)     :: Options(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:), RTSolution_pert(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_TL(:,:), RTSolution_AD(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_K(:,:)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atm_K(:,:)
  TYPE(CRTM_Surface_type),    ALLOCATABLE :: Sfc_K(:,:)

  CALL CRTM_Version(Version)
  WRITE(*,'(/5x,a)') 'Vector-RT (n_Stokes>1) cloud-scattering TL/AD/K verification'
  WRITE(*,'(5x,a/)') 'CRTM Version: '//TRIM(Version)

  ! --------------------------------------------------------------------------
  ! Initialize CRTM with the experimental cloud-optics scheme (6 phase
  ! elements -> the n_Stokes>1 scattering guard admits the run)
  ! --------------------------------------------------------------------------
  Error_Status = CRTM_Init( (/ SENSOR /), ChannelInfo,      &
                            Cloud_Model       = 'CRTM-Exp', &
                            CloudCoeff_File   = LUT,        &
                            CloudCoeff_Format = 'netCDF',   &
                            File_Path         = PATH,       &
                            Quiet             = .TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Init (Cloud_Model=CRTM-Exp) failed', FAILURE )
    STOP 1
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
  IF ( n_Channels < 1 ) THEN
    CALL Display_Message( PROGRAM_NAME, 'no channels loaded for '//SENSOR, FAILURE )
    STOP 1
  END IF

  ALLOCATE( RTSolution(n_Channels,N_PROFILES), RTSolution_pert(n_Channels,N_PROFILES), &
            RTSolution_TL(n_Channels,N_PROFILES), RTSolution_AD(n_Channels,N_PROFILES), &
            RTSolution_K(n_Channels,N_PROFILES), &
            Atm_K(n_Channels,N_PROFILES), Sfc_K(n_Channels,N_PROFILES), &
            STAT=Allocate_Status )
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
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Atmosphere_Create failed', FAILURE )
    STOP 1
  END IF

  ! Base column for every profile, with an overcast snow band
  CALL Load_ECMWF84_Atm_Data()         ! fills Atm(1)
  DO m = 2, N_PROFILES
    Atm(m) = Atm(1)
  END DO
  DO m = 1, N_PROFILES
    Atm(m)%n_Clouds                           = 1
    Atm(m)%Cloud_Fraction                     = ZERO
    Atm(m)%Cloud_Fraction(KC1:KC2)            = ONE      ! overcast: isolate the solver
    Atm(m)%Cloud(1)%Type                      = SNOW_CLOUD
    Atm(m)%Cloud(1)%Effective_Radius          = ZERO
    Atm(m)%Cloud(1)%Water_Content             = ZERO
    Atm(m)%Cloud(1)%Effective_Radius(KC1:KC2) = REFF_S
    Atm(m)%Cloud(1)%Water_Content(KC1:KC2)    = WC_S * (ONE + 0.2_fp*REAL(m-1,fp))
  END DO

  ! Congruent TL/AD/K input atmospheres
  DO m = 1, N_PROFILES
    Atm_TL(m)%Climatology = Atm(m)%Climatology
    Atm_TL(m)%Absorber_ID = Atm(m)%Absorber_ID ; Atm_TL(m)%Absorber_Units = Atm(m)%Absorber_Units
    Atm_TL(m)%Cloud(1)%Type = Atm(m)%Cloud(1)%Type
    Atm_AD(m)%Climatology = Atm(m)%Climatology
    Atm_AD(m)%Absorber_ID = Atm(m)%Absorber_ID ; Atm_AD(m)%Absorber_Units = Atm(m)%Absorber_Units
    Atm_AD(m)%Cloud(1)%Type = Atm(m)%Cloud(1)%Type
    DO l = 1, n_Channels
      Atm_K(l,m)%Climatology = Atm(m)%Climatology
      Atm_K(l,m)%Absorber_ID = Atm(m)%Absorber_ID ; Atm_K(l,m)%Absorber_Units = Atm(m)%Absorber_Units
      Atm_K(l,m)%Cloud(1)%Type = Atm(m)%Cloud(1)%Type
    END DO
  END DO

  ! Ocean surface + geometry
  DO m = 1, N_PROFILES
    Sfc(m)%Water_Coverage    = ONE
    Sfc(m)%Water_Type        = 1
    Sfc(m)%Water_Temperature = 290.0_fp
    Sfc(m)%Wind_Speed        = 12.0_fp
    Sfc(m)%Wind_Direction    = WIND_DIR
    Sfc(m)%Salinity          = 33.0_fp
    ! Non-zero relative wind azimuth. At zero the surface third and fourth
    ! Stokes components vanish identically (they are odd harmonics), which
    ! would make every U and V check below vacuous.
    CALL CRTM_Geometry_SetValue( Geometry(m), Sensor_Zenith_Angle  = ZENITH, &
                                              Sensor_Azimuth_Angle = SENSOR_AZI )
  END DO

  ! --------------------------------------------------------------------------
  ! Scalar control (n_Stokes = 1): validates the harness on the proven path
  ! --------------------------------------------------------------------------
  CALL set_options( 1 )
  WRITE(*,'(/5x,"=========== n_Stokes = 1 scalar control (ADA, overcast snow) ===========")')
  CALL check_fd ( 0, VAR_WC, ok_s1_fd )
  CALL check_adj( 1, ok_s1_adj )
  CALL check_k  ( 1, ok_s1_k )

  ! --------------------------------------------------------------------------
  ! Vector RT (n_Stokes = 2)
  ! --------------------------------------------------------------------------
  CALL set_options( 2 )
  WRITE(*,'(/5x,"=========== n_Stokes = 2 vector RT (ADA, overcast snow) ===========")')
  CALL check_fd ( 1, VAR_WC, ok_v_fd_wc1 )   ! dI/dWC
  CALL check_fd ( 2, VAR_WC, ok_v_fd_wc2 )   ! dQ/dWC (polarized phase chain)
  CALL check_fd ( 1, VAR_T , ok_v_fd_t1  )   ! dI/dT  (thermal source)
  CALL check_fd ( 2, VAR_T , ok_v_fd_t2  )   ! dQ/dT
  CALL check_adj( 2, ok_v_adj )              ! full-Stokes dot product
  CALL check_k  ( 2, ok_v_k )

  ! --------------------------------------------------------------------------
  ! Vector RT with scattering switched OFF (n_Stokes = 2). Dropping the water
  ! content below CRTM's scattering trigger routes the solve to CRTM_Emission
  ! plus CRTM_Emission_Stokes instead of ADA, which is a different code path
  ! with its own tangent-linear and adjoint. Without this block that path has
  ! no Jacobian coverage at all, and it is the path a clear-sky polarimetric
  ! run takes, which is the main use for ocean wind-vector work.
  ! Water content is the wrong control variable here (there is no cloud left to
  ! perturb), so the checks drive temperature, which reaches Stokes Q through
  ! both the surface Planck term and the reflected downwelling.
  ! --------------------------------------------------------------------------
  CALL set_wc( WC_CLEAR )
  WRITE(*,'(/5x,"=========== n_Stokes = 2 vector RT (no scattering, Emission) ===========")')
  CALL check_fd ( 1, VAR_T, ok_c_fd_t1 )     ! dI/dT
  CALL check_fd ( 2, VAR_T, ok_c_fd_t2 )     ! dQ/dT (polarized surface chain)
  CALL check_adj( 2, ok_c_adj )              ! full-Stokes dot product
  CALL check_k  ( 2, ok_c_k )
  CALL set_wc( WC_S )                        ! restore the scattering column

  ! --------------------------------------------------------------------------
  ! Vector RT with FRACTIONAL cloud cover (n_Stokes = 2). The blocks above are
  ! deliberately overcast so the solver is isolated from the clear/cloudy
  ! combine. This one exercises that combine: the forward model blends every
  ! Stokes component of the clear and cloudy columns, so its adjoint has to
  ! seed every Stokes component of both. Seeding %Radiance alone leaves the
  ! clear-sky half of the vector Jacobian unseeded, which the dot-product
  ! identity detects and a K-vs-AD check cannot.
  ! --------------------------------------------------------------------------
  CALL set_cfrac( 0.5_fp )
  WRITE(*,'(/5x,"=========== n_Stokes = 2 vector RT (fractional cloud) ===========")')
  CALL check_fd ( 2, VAR_WC, ok_f_fd )       ! dQ/dWC through the combine
  CALL check_adj( 2, ok_f_adj )              ! full-Stokes dot product
  CALL check_k  ( 2, ok_f_k )
  CALL set_cfrac( ONE )                      ! restore

  ! --------------------------------------------------------------------------
  ! FULL Stokes vector (n_Stokes = 4). Everything above runs at n_Stokes = 2,
  ! so the polarized phase-matrix blocks that only exist beyond two components,
  ! (1,3) (3,1) (2,3) (3,2) (3,3) (2,4) (4,2) (3,4) (4,3) (4,4), have never been
  ! differentiated, and neither has the U/V chain from the surface through the
  ! azimuthal accumulation to the reported Stokes vector.
  !
  ! FASTEM4 is loaded here because the FASTEM6 default has no third or fourth
  ! Stokes azimuth model at all, so U and V would be identically zero and every
  ! check below would pass vacuously.
  ! --------------------------------------------------------------------------
  Error_Status = CRTM_MWwaterCoeff_Load_FASTEM( 'FASTEM4', Quiet=.TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'FASTEM4 load failed', FAILURE ); STOP 1
  END IF
  CALL set_options( 4 )
  WRITE(*,'(/5x,"=========== n_Stokes = 4 full Stokes (ADA, overcast snow) ===========")')
  CALL check_fd ( 1, VAR_WC, ok_4_fd1 )      ! dI/dWC
  CALL check_fd ( 2, VAR_WC, ok_4_fd2 )      ! dQ/dWC
  CALL check_fd ( 3, VAR_WC, ok_4_fd3 )      ! dU/dWC
  CALL check_fd ( 4, VAR_WC, ok_4_fd4 )      ! dV/dWC
  ! The observable polarimetric microwave exists for: the third and fourth
  ! Stokes components respond to WIND DIRECTION, through the odd harmonics of
  ! the surface azimuth model. Nothing had ever differentiated the surface on
  ! the vector path, so this is the Jacobian the whole capability rests on.
  CALL check_fd ( 3, VAR_WDIR, ok_4_wdirU )  ! dU/d(wind direction)
  CALL check_fd ( 4, VAR_WDIR, ok_4_wdirV )  ! dV/d(wind direction)
  CALL check_fd ( 2, VAR_WSP , ok_4_wspQ  )  ! dQ/d(wind speed)
  CALL check_adj( 4, ok_4_adj )              ! four-component dot product
  CALL check_k  ( 4, ok_4_k )
  ! The reported Radiance is now the Stokes vector projected onto the channel
  ! polarization, so it has its own tangent linear and adjoint. Passing ks_out=0
  ! selects %Radiance rather than a Stokes component, and check_adj_radiance
  ! seeds %Radiance rather than %Stokes, which is the only way to exercise the
  ! transpose of the projection.
  CALL check_fd ( 0, VAR_WC, ok_4_fdR )      ! d(projected Radiance)/dWC
  CALL check_adj_radiance( ok_4_adjR )

  ! Same two checks with FRACTIONAL cloud. The overcast block above has a total
  ! cloud cover of one, which makes the clear/cloudy split of the reported
  ! radiance degenerate: the cloudy column gets everything either way. Only a
  ! cover strictly between zero and one distinguishes a correct split from an
  ! absent one, and only the reported radiance exercises it, since the Stokes
  ! seeds are split separately.
  CALL set_cfrac( 0.5_fp )
  WRITE(*,'(/5x,"=========== n_Stokes = 4, fractional, reported radiance ===========")')
  CALL check_fd ( 0, VAR_WC, ok_4f_fdR )
  CALL check_adj_radiance( ok_4f_adjR )
  CALL set_cfrac( ONE )

  Error_Status = CRTM_Destroy( ChannelInfo )

  WRITE(*,'(/5x,a)') '====================================================='
  WRITE(*,'(5x,"scalar control  TL vs FD (dI/dWC)     : ",a)') MERGE('PASS','FAIL',ok_s1_fd)
  WRITE(*,'(5x,"scalar control  adjoint dot-product   : ",a)') MERGE('PASS','FAIL',ok_s1_adj)
  WRITE(*,'(5x,"scalar control  K vs AD               : ",a)') MERGE('PASS','FAIL',ok_s1_k)
  WRITE(*,'(5x,"n_Stokes=2      TL vs FD (dI/dWC)     : ",a)') MERGE('PASS','FAIL',ok_v_fd_wc1)
  WRITE(*,'(5x,"n_Stokes=2      TL vs FD (dQ/dWC)     : ",a)') MERGE('PASS','FAIL',ok_v_fd_wc2)
  WRITE(*,'(5x,"n_Stokes=2      TL vs FD (dI/dT)      : ",a)') MERGE('PASS','FAIL',ok_v_fd_t1)
  WRITE(*,'(5x,"n_Stokes=2      TL vs FD (dQ/dT)      : ",a)') MERGE('PASS','FAIL',ok_v_fd_t2)
  WRITE(*,'(5x,"n_Stokes=2      adjoint dot-product   : ",a)') MERGE('PASS','FAIL',ok_v_adj)
  WRITE(*,'(5x,"n_Stokes=2      K vs AD               : ",a)') MERGE('PASS','FAIL',ok_v_k)
  WRITE(*,'(5x,"n_Stokes=2 clear TL vs FD (dI/dT)     : ",a)') MERGE('PASS','FAIL',ok_c_fd_t1)
  WRITE(*,'(5x,"n_Stokes=2 clear TL vs FD (dQ/dT)     : ",a)') MERGE('PASS','FAIL',ok_c_fd_t2)
  WRITE(*,'(5x,"n_Stokes=2 clear adjoint dot-product  : ",a)') MERGE('PASS','FAIL',ok_c_adj)
  WRITE(*,'(5x,"n_Stokes=2 clear K vs AD              : ",a)') MERGE('PASS','FAIL',ok_c_k)
  WRITE(*,'(5x,"n_Stokes=2 frac  TL vs FD (dQ/dWC)    : ",a)') MERGE('PASS','FAIL',ok_f_fd)
  WRITE(*,'(5x,"n_Stokes=2 frac  adjoint dot-product  : ",a)') MERGE('PASS','FAIL',ok_f_adj)
  WRITE(*,'(5x,"n_Stokes=2 frac  K vs AD              : ",a)') MERGE('PASS','FAIL',ok_f_k)
  WRITE(*,'(5x,"n_Stokes=4      TL vs FD (dI/dWC)     : ",a)') MERGE('PASS','FAIL',ok_4_fd1)
  WRITE(*,'(5x,"n_Stokes=4      TL vs FD (dQ/dWC)     : ",a)') MERGE('PASS','FAIL',ok_4_fd2)
  WRITE(*,'(5x,"n_Stokes=4      TL vs FD (dU/dWC)     : ",a)') MERGE('PASS','FAIL',ok_4_fd3)
  WRITE(*,'(5x,"n_Stokes=4      TL vs FD (dV/dWC)     : ",a)') MERGE('PASS','FAIL',ok_4_fd4)
  WRITE(*,'(5x,"n_Stokes=4      TL vs FD (dU/dWindDir)  : ",a)') MERGE('PASS','FAIL',ok_4_wdirU)
  WRITE(*,'(5x,"n_Stokes=4      TL vs FD (dV/dWindDir)  : ",a)') MERGE('PASS','FAIL',ok_4_wdirV)
  WRITE(*,'(5x,"n_Stokes=4      TL vs FD (dQ/dWindSpd)  : ",a)') MERGE('PASS','FAIL',ok_4_wspQ)
  WRITE(*,'(5x,"n_Stokes=4      adjoint dot-product   : ",a)') MERGE('PASS','FAIL',ok_4_adj)
  WRITE(*,'(5x,"n_Stokes=4      K vs AD               : ",a)') MERGE('PASS','FAIL',ok_4_k)
  WRITE(*,'(5x,"n_Stokes=4      TL vs FD (dRadiance)   : ",a)') MERGE('PASS','FAIL',ok_4_fdR)
  WRITE(*,'(5x,"n_Stokes=4      adjoint via %Radiance  : ",a)') MERGE('PASS','FAIL',ok_4_adjR)
  WRITE(*,'(5x,"n_Stokes=4 frac TL vs FD (dRadiance)   : ",a)') MERGE('PASS','FAIL',ok_4f_fdR)
  WRITE(*,'(5x,"n_Stokes=4 frac adjoint via %Radiance  : ",a)') MERGE('PASS','FAIL',ok_4f_adjR)
  IF ( ok_s1_fd .AND. ok_s1_adj .AND. ok_s1_k .AND. &
       ok_v_fd_wc1 .AND. ok_v_fd_wc2 .AND. ok_v_fd_t1 .AND. ok_v_fd_t2 .AND. &
       ok_v_adj .AND. ok_v_k .AND. &
       ok_c_fd_t1 .AND. ok_c_fd_t2 .AND. ok_c_adj .AND. ok_c_k .AND. &
       ok_f_fd .AND. ok_f_adj .AND. ok_f_k .AND. &
       ok_4_fd1 .AND. ok_4_fd2 .AND. ok_4_fd3 .AND. ok_4_fd4 .AND. &
       ok_4_adj .AND. ok_4_k .AND. ok_4_fdR .AND. ok_4_adjR .AND. &
       ok_4f_fdR .AND. ok_4f_adjR .AND. &
       ok_4_wdirU .AND. ok_4_wdirV .AND. ok_4_wspQ ) THEN
    WRITE(*,'(5x,a)') 'ALL CHECKS PASSED'
    STOP 0
  ELSE
    WRITE(*,'(5x,a)') 'CHECKS FAILED'
    STOP 1
  END IF

CONTAINS

  SUBROUTINE set_options( ns )
    INTEGER, INTENT(IN) :: ns
    INTEGER :: mm
    DO mm = 1, N_PROFILES
      Options(mm)%n_Stokes        = ns
      Options(mm)%RT_Algorithm_Id = RT_ADA
    END DO
  END SUBROUTINE set_options

  ! Reset the cloud water content on every profile, preserving the per-profile
  ! spread the scattering blocks rely on.
  SUBROUTINE set_wc( wc )
    REAL(fp), INTENT(IN) :: wc
    INTEGER :: mm
    DO mm = 1, N_PROFILES
      Atm(mm)%Cloud(1)%Water_Content(KC1:KC2) = wc * (ONE + 0.2_fp*REAL(mm-1,fp))
    END DO
  END SUBROUTINE set_wc

  ! Cloud fraction in the cloud band. ONE is overcast; anything strictly between
  ! zero and one routes the run through the clear/cloudy combine.
  SUBROUTINE set_cfrac( cf )
    REAL(fp), INTENT(IN) :: cf
    INTEGER :: mm
    DO mm = 1, N_PROFILES
      Atm(mm)%Cloud_Fraction(KC1:KC2) = cf
    END DO
  END SUBROUTINE set_cfrac

  ! Selected output: Stokes component ks (n_Stokes>1) or the scalar Radiance (ks=0)
  REAL(fp) FUNCTION get_out( rts, ks )
    TYPE(CRTM_RTSolution_type), INTENT(IN) :: rts
    INTEGER, INTENT(IN) :: ks
    IF ( ks > 0 ) THEN
      get_out = rts%Stokes(ks)
    ELSE
      get_out = rts%Radiance
    END IF
  END FUNCTION get_out

  SUBROUTINE set_seed( rts, ks, val )
    TYPE(CRTM_RTSolution_type), INTENT(INOUT) :: rts
    INTEGER,  INTENT(IN) :: ks
    REAL(fp), INTENT(IN) :: val
    IF ( ks > 0 ) THEN
      rts%Stokes(ks) = val
    ELSE
      rts%Radiance = val
    END IF
  END SUBROUTINE set_seed

  ! Access the perturbed forward variable (profile 1, layer KP)
  REAL(fp) FUNCTION get_var( var )
    INTEGER, INTENT(IN) :: var
    SELECT CASE ( var )
      CASE ( VAR_WC )   ; get_var = Atm(1)%Cloud(1)%Water_Content(KP)
      CASE ( VAR_WSP )  ; get_var = Sfc(1)%Wind_Speed
      CASE ( VAR_WDIR ) ; get_var = Sfc(1)%Wind_Direction
      CASE DEFAULT      ; get_var = Atm(1)%Temperature(KP)
    END SELECT
  END FUNCTION get_var

  SUBROUTINE set_var( var, val )
    INTEGER,  INTENT(IN) :: var
    REAL(fp), INTENT(IN) :: val
    SELECT CASE ( var )
      CASE ( VAR_WC )   ; Atm(1)%Cloud(1)%Water_Content(KP) = val
      CASE ( VAR_WSP )  ; Sfc(1)%Wind_Speed     = val
      CASE ( VAR_WDIR ) ; Sfc(1)%Wind_Direction = val
      CASE DEFAULT      ; Atm(1)%Temperature(KP) = val
    END SELECT
  END SUBROUTINE set_var

  ! ----------------------------------------------------------------
  ! Check 1 : TL vs central finite difference, output Stokes(ks_out)
  !           (ks_out=0 -> Radiance), perturbing variable `var`.
  ! ----------------------------------------------------------------
  SUBROUTINE check_fd( ks_out, var, ok )
    INTEGER, INTENT(IN)  :: ks_out, var
    LOGICAL, INTENT(OUT) :: ok
    CHARACTER(16) :: vname, oname
    REAL(fp) :: tl, fd, Rp, Rm, ratio, best, delta, X0
    REAL(fp) :: fd_all(n_Channels)
    INTEGER  :: ii, kk, ch

    SELECT CASE ( var )
      CASE ( VAR_WC )   ; vname = 'Water_Content'
      CASE ( VAR_WSP )  ; vname = 'Wind_Speed'
      CASE ( VAR_WDIR ) ; vname = 'Wind_Direction'
      CASE DEFAULT      ; vname = 'Temperature'
    END SELECT
    IF ( ks_out > 0 ) THEN
      WRITE(oname,'("Stokes(",i0,")")') ks_out
    ELSE
      oname = 'Radiance'
    END IF

    ! TL with a unit perturbation of the variable
    CALL CRTM_Atmosphere_Zero( Atm_TL ) ; CALL CRTM_Surface_Zero( Sfc_TL )
    SELECT CASE ( var )
      CASE ( VAR_WC )   ; Atm_TL(1)%Cloud(1)%Water_Content(KP) = ONE
      CASE ( VAR_WSP )  ; Sfc_TL(1)%Wind_Speed     = ONE
      CASE ( VAR_WDIR ) ; Sfc_TL(1)%Wind_Direction = ONE
      CASE DEFAULT      ; Atm_TL(1)%Temperature(KP) = ONE
    END SELECT
    Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, ChannelInfo, &
                                        RTSolution, RTSolution_TL, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'TL fail' ; ok=.FALSE. ; RETURN ; END IF

    ! Channel selection by FD probe (not max|TL|: a broken zero TL must not
    ! hide itself), restricted to scattering channels.
    X0 = get_var( var )
    delta = ABS(X0) * 0.1_fp / 256.0_fp
    CALL set_var( var, X0 + delta )
    Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution_pert, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'FD forward fail' ; ok=.FALSE. ; RETURN ; END IF
    DO ii = 1, n_Channels ; fd_all(ii) = get_out(RTSolution_pert(ii,1),ks_out) ; END DO
    CALL set_var( var, X0 - delta )
    Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution_pert, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'FD forward fail' ; ok=.FALSE. ; RETURN ; END IF
    DO ii = 1, n_Channels
      fd_all(ii) = ( fd_all(ii) - get_out(RTSolution_pert(ii,1),ks_out) ) / ( 2.0_fp*delta )
    END DO
    CALL set_var( var, X0 )
    ch = 0 ; best = ZERO
    DO ii = 1, n_Channels
      ! Surface variables reach every channel; the scattering restriction only
      ! makes sense for the cloud/temperature probes.
      IF ( var /= VAR_WSP .AND. var /= VAR_WDIR ) THEN
        IF ( .NOT. RTSolution(ii,1)%Scattering_Flag ) CYCLE
      END IF
      IF ( ABS(fd_all(ii)) >= best ) THEN ; best = ABS(fd_all(ii)) ; ch = ii ; END IF
    END DO
    IF ( ch == 0 ) ch = 1
    tl = get_out(RTSolution_TL(ch,1),ks_out)

    best = HUGE(ONE)
    WRITE(*,'(/7x,"[FD] d ",a," / d ",a,"(",i0,")   channel ",i0,"   TL=",es13.6)') &
          TRIM(oname), TRIM(vname), KP, RTSolution(ch,1)%Sensor_Channel, tl
    DO kk = 4, 14
      delta = ABS(X0) * 0.1_fp / (2.0_fp**kk)
      CALL set_var( var, X0 + delta )
      Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution_pert, Options=Options )
      IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'FD forward fail' ; ok=.FALSE. ; RETURN ; END IF
      Rp = get_out(RTSolution_pert(ch,1),ks_out)
      CALL set_var( var, X0 - delta )
      Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution_pert, Options=Options )
      IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'FD forward fail' ; ok=.FALSE. ; RETURN ; END IF
      Rm = get_out(RTSolution_pert(ch,1),ks_out)
      CALL set_var( var, X0 )
      fd = ( Rp - Rm ) / ( 2.0_fp*delta )
      ratio = fd / tl
      IF ( ABS(ratio-ONE) < best ) best = ABS(ratio-ONE)
      WRITE(*,'(9x,"delta=",es10.3,"  FD=",es16.9,"  FD/TL=",f14.10)') delta, fd, ratio
    END DO
    ok = ( best < TOL_FD )
    WRITE(*,'(7x,"-> best |FD/TL - 1| = ",es11.4,"   ",a)') best, MERGE('PASS','FAIL',ok)
  END SUBROUTINE check_fd

  ! ----------------------------------------------------------------
  ! Check 2 : adjoint dot-product over the full Stokes vector,
  !           x spanning Temperature + Water_Content everywhere.
  ! ----------------------------------------------------------------
  SUBROUTINE check_adj( ns, ok )
    INTEGER, INTENT(IN)  :: ns
    LOGICAL, INTENT(OUT) :: ok
    REAL(fp) :: LHS, RHS, dy, rel_adj, RHS_sfc
    INTEGER  :: ii, ks, mm

    ! The surface is perturbed alongside the atmosphere. Without it the identity
    ! never touches the surface adjoint at all, and a completely broken
    ! SfcOptics_AD would satisfy it: Sfc_TL was zeroed and Sfc_AD never read.
    ! Wind direction is the one that matters, because it is the only route to
    ! the third and fourth Stokes components and so the observable that
    ! polarimetric microwave exists for.
    CALL CRTM_Atmosphere_Zero( Atm_TL ) ; CALL CRTM_Surface_Zero( Sfc_TL )
    DO mm = 1, N_PROFILES
      DO ii = 1, N_LAYERS
        Atm_TL(mm)%Temperature(ii)            = 0.5_fp * SIN( 0.7_fp*REAL(ii,fp) + 1.3_fp*REAL(mm,fp) )
        Atm_TL(mm)%Cloud(1)%Water_Content(ii) = 0.1_fp * COS( 0.9_fp*REAL(ii,fp) + 0.4_fp*REAL(mm,fp) )
      END DO
      Sfc_TL(mm)%Wind_Speed        = 0.30_fp + 0.10_fp*REAL(mm,fp)
      Sfc_TL(mm)%Wind_Direction    = 2.00_fp - 0.50_fp*REAL(mm,fp)
      Sfc_TL(mm)%Water_Temperature = 0.20_fp + 0.05_fp*REAL(mm,fp)
      Sfc_TL(mm)%Salinity          = 0.10_fp
    END DO
    Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, ChannelInfo, &
                                        RTSolution, RTSolution_TL, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'TL fail' ; ok=.FALSE. ; RETURN ; END IF

    LHS = ZERO
    CALL CRTM_RTSolution_Zero( RTSolution_AD )
    DO mm = 1, N_PROFILES
      DO l = 1, n_Channels
        IF ( ns > 1 ) THEN
          DO ks = 1, ns
            dy = RTSolution_TL(l,mm)%Stokes(ks)
            LHS = LHS + dy*dy
            RTSolution_AD(l,mm)%Stokes(ks) = dy
          END DO
        ELSE
          dy = RTSolution_TL(l,mm)%Radiance
          LHS = LHS + dy*dy
          RTSolution_AD(l,mm)%Radiance = dy
        END IF
      END DO
    END DO

    CALL CRTM_Atmosphere_Zero( Atm_AD ) ; CALL CRTM_Surface_Zero( Sfc_AD )
    Error_Status = CRTM_Adjoint( Atm, Sfc, RTSolution_AD, Geometry, ChannelInfo, &
                                 Atm_AD, Sfc_AD, RTSolution, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'AD fail' ; ok=.FALSE. ; RETURN ; END IF

    RHS = ZERO ; RHS_sfc = ZERO
    DO mm = 1, N_PROFILES
      DO ii = 1, N_LAYERS
        RHS = RHS + Atm_TL(mm)%Temperature(ii)            * Atm_AD(mm)%Temperature(ii)
        RHS = RHS + Atm_TL(mm)%Cloud(1)%Water_Content(ii) * Atm_AD(mm)%Cloud(1)%Water_Content(ii)
      END DO
      RHS_sfc = RHS_sfc &
              + Sfc_TL(mm)%Wind_Speed        * Sfc_AD(mm)%Wind_Speed        &
              + Sfc_TL(mm)%Wind_Direction    * Sfc_AD(mm)%Wind_Direction    &
              + Sfc_TL(mm)%Water_Temperature * Sfc_AD(mm)%Water_Temperature &
              + Sfc_TL(mm)%Salinity          * Sfc_AD(mm)%Salinity
    END DO
    RHS = RHS + RHS_sfc
    rel_adj = ABS(LHS-RHS) / MAX(ABS(LHS), TINY(ONE))
    ok = ( rel_adj < TOL_ADJ )
    WRITE(*,'(/7x,"[ADJ] <dy,dy>=",es16.9,"   <x,gx>=",es16.9)') LHS, RHS
    ! Surface share, so it is visible whether the surface adjoint is actually
    ! being tested or merely swamped by the atmospheric terms.
    WRITE(*,'(7x,"surface share of <x,gx> = ",f8.4," %")') 100.0_fp*RHS_sfc/MAX(ABS(RHS),TINY(ONE))
    WRITE(*,'(7x,"-> relative difference = ",es11.4,"   ",a)') rel_adj, MERGE('PASS','FAIL',ok)
  END SUBROUTINE check_adj

  ! Adjoint dot product seeded through %Radiance instead of %Stokes. This is
  ! the transpose of the channel-polarization projection, and nothing else
  ! exercises it: every other check seeds Stokes components directly, which
  ! bypasses the projection entirely.
  SUBROUTINE check_adj_radiance( ok )
    LOGICAL, INTENT(OUT) :: ok
    REAL(fp) :: LHS, RHS, dy, rel_adj
    INTEGER  :: ii, mm

    CALL CRTM_Atmosphere_Zero( Atm_TL ) ; CALL CRTM_Surface_Zero( Sfc_TL )
    DO mm = 1, N_PROFILES
      DO ii = 1, N_LAYERS
        Atm_TL(mm)%Temperature(ii)            = 0.5_fp * SIN( 0.7_fp*REAL(ii,fp) + 1.3_fp*REAL(mm,fp) )
        Atm_TL(mm)%Cloud(1)%Water_Content(ii) = 0.1_fp * COS( 0.9_fp*REAL(ii,fp) + 0.4_fp*REAL(mm,fp) )
      END DO
    END DO
    Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, ChannelInfo, &
                                        RTSolution, RTSolution_TL, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'TL fail' ; ok=.FALSE. ; RETURN ; END IF

    LHS = ZERO
    CALL CRTM_RTSolution_Zero( RTSolution_AD )
    DO mm = 1, N_PROFILES
      DO l = 1, n_Channels
        dy = RTSolution_TL(l,mm)%Radiance
        LHS = LHS + dy*dy
        RTSolution_AD(l,mm)%Radiance = dy
      END DO
    END DO

    CALL CRTM_Atmosphere_Zero( Atm_AD ) ; CALL CRTM_Surface_Zero( Sfc_AD )
    Error_Status = CRTM_Adjoint( Atm, Sfc, RTSolution_AD, Geometry, ChannelInfo, &
                                 Atm_AD, Sfc_AD, RTSolution, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'AD fail' ; ok=.FALSE. ; RETURN ; END IF

    RHS = ZERO
    DO mm = 1, N_PROFILES
      DO ii = 1, N_LAYERS
        RHS = RHS + Atm_TL(mm)%Temperature(ii)            * Atm_AD(mm)%Temperature(ii)
        RHS = RHS + Atm_TL(mm)%Cloud(1)%Water_Content(ii) * Atm_AD(mm)%Cloud(1)%Water_Content(ii)
      END DO
    END DO
    rel_adj = ABS(LHS-RHS) / MAX(ABS(LHS), TINY(ONE))
    ok = ( rel_adj < TOL_ADJ )
    WRITE(*,'(/7x,"[ADJ-R] <dy,dy>=",es16.9,"   <x,gx>=",es16.9)') LHS, RHS
    WRITE(*,'(7x,"-> relative difference = ",es11.4,"   ",a)') rel_adj, MERGE('PASS','FAIL',ok)
  END SUBROUTINE check_adj_radiance

  ! ----------------------------------------------------------------
  ! Check 3 : K-Matrix vs Adjoint Jacobian (Stokes(1) seed; Temperature
  !           and Water_Content columns), one channel/profile.
  ! ----------------------------------------------------------------
  SUBROUTINE check_k( ns, ok )
    INTEGER, INTENT(IN)  :: ns
    LOGICAL, INTENT(OUT) :: ok
    REAL(fp) :: maxdiff, scal, rel_k
    INTEGER  :: ks0, l0, m0, mm

    ks0 = MERGE( 1, 0, ns > 1 )   ! Stokes(1) seed for vector, Radiance for scalar
    l0 = 1 ; m0 = 1

    CALL CRTM_Atmosphere_Zero( Atm_K ) ; CALL CRTM_Surface_Zero( Sfc_K )
    CALL CRTM_RTSolution_Zero( RTSolution_K )
    DO mm = 1, N_PROFILES
      DO l = 1, n_Channels
        CALL set_seed( RTSolution_K(l,mm), ks0, ONE )
      END DO
    END DO
    Error_Status = CRTM_K_Matrix( Atm, Sfc, RTSolution_K, Geometry, ChannelInfo, &
                                  Atm_K, Sfc_K, RTSolution, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'K fail' ; ok=.FALSE. ; RETURN ; END IF

    CALL CRTM_Atmosphere_Zero( Atm_AD ) ; CALL CRTM_Surface_Zero( Sfc_AD )
    CALL CRTM_RTSolution_Zero( RTSolution_AD )
    CALL set_seed( RTSolution_AD(l0,m0), ks0, ONE )
    Error_Status = CRTM_Adjoint( Atm, Sfc, RTSolution_AD, Geometry, ChannelInfo, &
                                 Atm_AD, Sfc_AD, RTSolution, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'AD fail' ; ok=.FALSE. ; RETURN ; END IF

    maxdiff = MAX( MAXVAL( ABS( Atm_K(l0,m0)%Temperature - Atm_AD(m0)%Temperature ) ),  &
                   MAXVAL( ABS( Atm_K(l0,m0)%Cloud(1)%Water_Content                     &
                                - Atm_AD(m0)%Cloud(1)%Water_Content ) ),                &
                   ABS( Sfc_K(l0,m0)%Wind_Speed        - Sfc_AD(m0)%Wind_Speed ),       &
                   ABS( Sfc_K(l0,m0)%Wind_Direction    - Sfc_AD(m0)%Wind_Direction ),   &
                   ABS( Sfc_K(l0,m0)%Water_Temperature - Sfc_AD(m0)%Water_Temperature ),&
                   ABS( Sfc_K(l0,m0)%Salinity          - Sfc_AD(m0)%Salinity ) )
    scal    = MAX( MAXVAL(ABS(Atm_K(l0,m0)%Temperature)), &
                   MAXVAL(ABS(Atm_K(l0,m0)%Cloud(1)%Water_Content)), &
                   ABS(Sfc_K(l0,m0)%Wind_Speed), ABS(Sfc_K(l0,m0)%Wind_Direction), &
                   ABS(Sfc_K(l0,m0)%Water_Temperature), TINY(ONE) )
    rel_k   = maxdiff / scal
    ok = ( rel_k < TOL_K )
    WRITE(*,'(/7x,"[K] K vs AD (channel ",i0,"):  max|K-AD|/max|K| = ",es11.4,"   ",a)') &
          RTSolution(l0,m0)%Sensor_Channel, rel_k, MERGE('PASS','FAIL',ok)
  END SUBROUTINE check_k

  INCLUDE 'Load_ECMWF84_Atm_Data.inc'

END PROGRAM test_VectorRT_TLADK
