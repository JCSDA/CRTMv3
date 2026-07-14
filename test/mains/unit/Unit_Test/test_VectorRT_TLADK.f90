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

  ! The adjoint dot-product tolerance is deliberately tight (the correct code
  ! achieves ~1e-15): a one-sided TL/AD inconsistency in the phase-normalization
  ! polarized blocks shows up at ~5e-11, which 1e-9 would let through.
  REAL(fp), PARAMETER :: TOL_FD  = 1.0e-3_fp      ! TL vs finite difference
  REAL(fp), PARAMETER :: TOL_ADJ = 1.0e-12_fp     ! adjoint dot-product
  REAL(fp), PARAMETER :: TOL_K   = 1.0e-9_fp      ! K vs AD

  ! Perturbation-variable selectors for the FD check
  INTEGER, PARAMETER :: VAR_WC = 1, VAR_T = 2

  CHARACTER(256) :: Version
  INTEGER :: Error_Status, Allocate_Status, n_Channels
  INTEGER :: l, m
  LOGICAL :: ok_s1_fd, ok_s1_adj, ok_s1_k
  LOGICAL :: ok_v_fd_wc1, ok_v_fd_wc2, ok_v_fd_t1, ok_v_fd_t2, ok_v_adj, ok_v_k

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
    Sfc(m)%Wind_Speed        = 6.0_fp
    Sfc(m)%Salinity          = 33.0_fp
    CALL CRTM_Geometry_SetValue( Geometry(m), Sensor_Zenith_Angle = ZENITH )
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
  IF ( ok_s1_fd .AND. ok_s1_adj .AND. ok_s1_k .AND. &
       ok_v_fd_wc1 .AND. ok_v_fd_wc2 .AND. ok_v_fd_t1 .AND. ok_v_fd_t2 .AND. &
       ok_v_adj .AND. ok_v_k ) THEN
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
    IF ( var == VAR_WC ) THEN
      get_var = Atm(1)%Cloud(1)%Water_Content(KP)
    ELSE
      get_var = Atm(1)%Temperature(KP)
    END IF
  END FUNCTION get_var

  SUBROUTINE set_var( var, val )
    INTEGER,  INTENT(IN) :: var
    REAL(fp), INTENT(IN) :: val
    IF ( var == VAR_WC ) THEN
      Atm(1)%Cloud(1)%Water_Content(KP) = val
    ELSE
      Atm(1)%Temperature(KP) = val
    END IF
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

    IF ( var == VAR_WC ) THEN ; vname = 'Water_Content' ; ELSE ; vname = 'Temperature' ; END IF
    IF ( ks_out > 0 ) THEN
      WRITE(oname,'("Stokes(",i0,")")') ks_out
    ELSE
      oname = 'Radiance'
    END IF

    ! TL with a unit perturbation of the variable
    CALL CRTM_Atmosphere_Zero( Atm_TL ) ; CALL CRTM_Surface_Zero( Sfc_TL )
    IF ( var == VAR_WC ) THEN
      Atm_TL(1)%Cloud(1)%Water_Content(KP) = ONE
    ELSE
      Atm_TL(1)%Temperature(KP) = ONE
    END IF
    Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, ChannelInfo, &
                                        RTSolution, RTSolution_TL, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'TL fail' ; ok=.FALSE. ; RETURN ; END IF

    ! Channel selection by FD probe (not max|TL|: a broken zero TL must not
    ! hide itself), restricted to scattering channels.
    X0 = get_var( var )
    delta = ABS(X0) * 0.1_fp / 256.0_fp
    CALL set_var( var, X0 + delta )
    Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution_pert, Options=Options )
    DO ii = 1, n_Channels ; fd_all(ii) = get_out(RTSolution_pert(ii,1),ks_out) ; END DO
    CALL set_var( var, X0 - delta )
    Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution_pert, Options=Options )
    DO ii = 1, n_Channels
      fd_all(ii) = ( fd_all(ii) - get_out(RTSolution_pert(ii,1),ks_out) ) / ( 2.0_fp*delta )
    END DO
    CALL set_var( var, X0 )
    ch = 0 ; best = ZERO
    DO ii = 1, n_Channels
      IF ( .NOT. RTSolution(ii,1)%Scattering_Flag ) CYCLE
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
      Rp = get_out(RTSolution_pert(ch,1),ks_out)
      CALL set_var( var, X0 - delta )
      Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution_pert, Options=Options )
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
    REAL(fp) :: LHS, RHS, dy, rel_adj
    INTEGER  :: ii, ks, mm

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

    RHS = ZERO
    DO mm = 1, N_PROFILES
      DO ii = 1, N_LAYERS
        RHS = RHS + Atm_TL(mm)%Temperature(ii)            * Atm_AD(mm)%Temperature(ii)
        RHS = RHS + Atm_TL(mm)%Cloud(1)%Water_Content(ii) * Atm_AD(mm)%Cloud(1)%Water_Content(ii)
      END DO
    END DO
    rel_adj = ABS(LHS-RHS) / MAX(ABS(LHS), TINY(ONE))
    ok = ( rel_adj < TOL_ADJ )
    WRITE(*,'(/7x,"[ADJ] <dy,dy>=",es16.9,"   <x,gx>=",es16.9)') LHS, RHS
    WRITE(*,'(7x,"-> relative difference = ",es11.4,"   ",a)') rel_adj, MERGE('PASS','FAIL',ok)
  END SUBROUTINE check_adj

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
                                - Atm_AD(m0)%Cloud(1)%Water_Content ) ) )
    scal    = MAX( MAXVAL(ABS(Atm_K(l0,m0)%Temperature)), &
                   MAXVAL(ABS(Atm_K(l0,m0)%Cloud(1)%Water_Content)), TINY(ONE) )
    rel_k   = maxdiff / scal
    ok = ( rel_k < TOL_K )
    WRITE(*,'(/7x,"[K] K vs AD (channel ",i0,"):  max|K-AD|/max|K| = ",es11.4,"   ",a)') &
          RTSolution(l0,m0)%Sensor_Channel, rel_k, MERGE('PASS','FAIL',ok)
  END SUBROUTINE check_k

  INCLUDE 'Load_ECMWF84_Atm_Data.inc'

END PROGRAM test_VectorRT_TLADK
