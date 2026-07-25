!
! test_UV_NO2_TLAD
!
! TL/AD/K parity check for the UV/VIS scene-NO2 ODPS component
! (GROUP_UV_NO2, Group_Index = 8).
!
! The group-8 NO2 predictor blocks in ODPS_Predictor.f90 are a compact
! 3-predictor set {NO2_A, NO2_A*DT, NO2_A*DT2} for the analytic Beer-Lambert
! NO2 extinction (layer OD = sigma(T)*N, exactly linear in amount); the
! forward path builds clean but no test exercises the TL/AD/K consistency of
! the new blocks. This test initializes the TEMPO UV product on a clear-sky
! ECMWF84 ocean column (daytime solar geometry - the UV signal is reflected
! solar) with an added NO2 absorber (UARS climatology) and verifies:
!   1. TL vs central finite difference for column perturbations of
!      NO2 (relative), H2O (relative) and Temperature (additive) -- the NO2
!      check probes the new predictor block end-to-end; with a zero-base
!      parity file the T response also rides the NO2 DT/DT2 terms.
!   2. Adjoint dot-product  <dy,dy> == <x, AD(dy)>  with x spanning
!      Temperature, H2O and NO2 on every layer of every profile.
!   3. K-Matrix vs Adjoint Jacobian equality (Temperature, H2O and NO2
!      columns) on the most NO2-sensitive channel.
!
! The test adapts to the loaded TauCoeff per variable: if the forward FD
! probe shows no response (e.g. H2O against the parity-gate file whose base
! components are zero-predictor, or NO2 against a plain group-2 file) the
! TL must be IDENTICALLY ZERO -- that run is the consistency/backward-
! compatibility control. With a responding file the TL must converge to the
! FD derivative.
!
! Exit: STOP 0 if every check passes, STOP 1 otherwise.
!
! CREATION HISTORY:
!       Written by:     Benjamin Johnson, 25-Jul-2026
!                       Adapted from test_MW_O3_TLAD (b9c525a).
!
PROGRAM test_UV_NO2_TLAD

  USE CRTM_Module
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_UV_NO2_TLAD'
  CHARACTER(*), PARAMETER :: PATH   = './testinput/'
  CHARACTER(*), PARAMETER :: SENSOR = 'u.tempo_is40e'

  ! Profile / column setup (ECMWF84 ocean column; clear sky; daytime)
  INTEGER,  PARAMETER :: N_PROFILES  = 2
  INTEGER,  PARAMETER :: N_LAYERS    = 100
  INTEGER,  PARAMETER :: N_ABSORBERS = 7
  INTEGER,  PARAMETER :: N_CLOUDS    = 0
  INTEGER,  PARAMETER :: N_AEROSOLS  = 0
  REAL(fp), PARAMETER :: ZENITH      = 45.0_fp
  REAL(fp), PARAMETER :: SOLAR_ZEN   = 30.0_fp
  INTEGER,  PARAMETER :: IDX_H2O = 1, IDX_NO2 = 7  ! absorber slots (7 = added NO2)

  REAL(fp), PARAMETER :: TOL_FD   = 1.0e-3_fp     ! TL vs finite difference
  ! The adjoint dot-product tolerance is looser than the MW test's 1e-12: any
  ! dy here traverses the UV solar scattering RT (azimuth Fourier loop, 1028
  ! channels), whose TL and AD sums accumulate float64 roundoff in different
  ! orders (~5e-12 relative observed). The predictor mapping itself is pinned
  ! at machine precision by test_ODPS_NO2_Predictor_TLAD (no RT in the loop);
  ! a real transpose error would fail both tests by many orders.
  REAL(fp), PARAMETER :: TOL_ADJ      = 1.0e-10_fp ! combined T+H2O+NO2
  REAL(fp), PARAMETER :: TOL_ADJ_NO2  = 1.0e-10_fp ! NO2-only perturbation
  REAL(fp), PARAMETER :: TOL_K    = 1.0e-9_fp     ! K vs AD
  ! Forward UV radiances are O(1e0..1e2) mW/(m2.sr.cm-1); a column FD response
  ! below this is numerically indistinguishable from zero -> blind variable.
  REAL(fp), PARAMETER :: FD_ZERO  = 1.0e-9_fp

  ! Perturbation-variable selectors
  INTEGER, PARAMETER :: VAR_T = 1, VAR_H2O = 2, VAR_NO2 = 3

  CHARACTER(256) :: Version
  INTEGER :: Error_Status, Allocate_Status, n_Channels
  INTEGER :: l, m
  INTEGER :: ch_no2               ! most NO2-sensitive channel (0 if none)
  LOGICAL :: no2_active           ! loaded file responds to scene NO2
  LOGICAL :: ok_fd_no2, ok_fd_h2o, ok_fd_t, ok_adj, ok_adj_no2, ok_k

  REAL(fp) :: no2_clim(N_LAYERS)  ! UARS NO2 climatology, ppmv, 100-layer grid

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(1)
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

  CALL CRTM_Version(Version)
  WRITE(*,'(/5x,a)') 'UV scene-NO2 (GROUP_UV_NO2) TL/AD/K verification'
  WRITE(*,'(5x,a/)') 'CRTM Version: '//TRIM(Version)

  Error_Status = CRTM_Init( (/ SENSOR /), ChannelInfo, File_Path=PATH, Quiet=.TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Init failed', FAILURE )
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

  ! Base clear-sky column for every profile, plus the added NO2 absorber
  ! (slot 7, UARS climatology). Second profile gets a scaled NO2 column so
  ! both profiles contribute distinct NO2 signals to the adjoint dot-product.
  CALL Load_ECMWF84_Atm_Data()         ! fills Atm(1) and Atm(2), absorbers 1-6
  CALL Set_NO2_Climatology()
  DO m = 1, N_PROFILES
    Atm(m)%Absorber_Id(IDX_NO2)    = NO2_ID
    Atm(m)%Absorber_Units(IDX_NO2) = VOLUME_MIXING_RATIO_UNITS
    Atm(m)%Absorber(:,IDX_NO2)     = no2_clim
  END DO
  Atm(2)%Absorber(:,IDX_NO2) = 1.3_fp * Atm(2)%Absorber(:,IDX_NO2)

  ! Congruent TL/AD/K input atmospheres
  DO m = 1, N_PROFILES
    Atm_TL(m)%Climatology = Atm(m)%Climatology
    Atm_TL(m)%Absorber_ID = Atm(m)%Absorber_ID ; Atm_TL(m)%Absorber_Units = Atm(m)%Absorber_Units
    Atm_AD(m)%Climatology = Atm(m)%Climatology
    Atm_AD(m)%Absorber_ID = Atm(m)%Absorber_ID ; Atm_AD(m)%Absorber_Units = Atm(m)%Absorber_Units
    DO l = 1, n_Channels
      Atm_K(l,m)%Climatology = Atm(m)%Climatology
      Atm_K(l,m)%Absorber_ID = Atm(m)%Absorber_ID ; Atm_K(l,m)%Absorber_Units = Atm(m)%Absorber_Units
    END DO
  END DO

  ! Ocean surface + daytime geometry (solar zenith above horizon: the UV
  ! radiance is reflected solar, so a sun below the horizon would make every
  ! check vacuously zero)
  DO m = 1, N_PROFILES
    Sfc(m)%Water_Coverage    = ONE
    Sfc(m)%Water_Type        = 1
    Sfc(m)%Water_Temperature = 290.0_fp
    Sfc(m)%Wind_Speed        = 6.0_fp
    Sfc(m)%Salinity          = 33.0_fp
    CALL CRTM_Geometry_SetValue( Geometry(m), Sensor_Zenith_Angle = ZENITH, &
                                 Source_Zenith_Angle = SOLAR_ZEN )
  END DO

  ! --------------------------------------------------------------------------
  ! Checks. check_fd(VAR_NO2) also determines no2_active/ch_no2 from the
  ! forward FD probe, so it must run first.
  ! --------------------------------------------------------------------------
  CALL check_fd ( VAR_NO2, ok_fd_no2 )
  CALL check_fd ( VAR_H2O, ok_fd_h2o )
  CALL check_fd ( VAR_T  , ok_fd_t   )
  CALL check_adj( ok_adj )
  CALL check_adj_no2( ok_adj_no2 )
  CALL check_k  ( ok_k )

  Error_Status = CRTM_Destroy( ChannelInfo )

  WRITE(*,'(/5x,a)') '====================================================='
  IF ( no2_active ) THEN
    WRITE(*,'(5x,a)') 'TauCoeff mode : scene-NO2 ACTIVE (group-8 NO2 component)'
  ELSE
    WRITE(*,'(5x,a)') 'TauCoeff mode : NO2-blind control (no NO2 component)'
  END IF
  WRITE(*,'(5x,"TL vs FD (NO2 column)                 : ",a)') MERGE('PASS','FAIL',ok_fd_no2)
  WRITE(*,'(5x,"TL vs FD (H2O column)                 : ",a)') MERGE('PASS','FAIL',ok_fd_h2o)
  WRITE(*,'(5x,"TL vs FD (Temperature)                : ",a)') MERGE('PASS','FAIL',ok_fd_t)
  WRITE(*,'(5x,"adjoint dot-product (T+H2O+NO2)       : ",a)') MERGE('PASS','FAIL',ok_adj)
  WRITE(*,'(5x,"adjoint dot-product (NO2 only)        : ",a)') MERGE('PASS','FAIL',ok_adj_no2)
  WRITE(*,'(5x,"K vs AD (T/H2O/NO2 columns)           : ",a)') MERGE('PASS','FAIL',ok_k)
  IF ( ok_fd_no2 .AND. ok_fd_h2o .AND. ok_fd_t .AND. ok_adj .AND. ok_adj_no2 .AND. ok_k ) THEN
    WRITE(*,'(5x,a)') 'ALL CHECKS PASSED'
    STOP 0
  ELSE
    WRITE(*,'(5x,a)') 'CHECKS FAILED'
    STOP 1
  END IF

CONTAINS

  ! Apply a column perturbation of size eps to profile 1 of the forward state:
  ! multiplicative (1+eps) on the absorber columns, additive eps [K] on T.
  SUBROUTINE perturb( var, base, eps )
    INTEGER,  INTENT(IN) :: var
    REAL(fp), INTENT(IN) :: base(N_LAYERS)
    REAL(fp), INTENT(IN) :: eps
    SELECT CASE ( var )
      CASE ( VAR_T )   ; Atm(1)%Temperature        = base + eps
      CASE ( VAR_H2O ) ; Atm(1)%Absorber(:,IDX_H2O) = base * (ONE + eps)
      CASE ( VAR_NO2 ) ; Atm(1)%Absorber(:,IDX_NO2) = base * (ONE + eps)
    END SELECT
  END SUBROUTINE perturb

  SUBROUTINE get_base( var, base )
    INTEGER,  INTENT(IN)  :: var
    REAL(fp), INTENT(OUT) :: base(N_LAYERS)
    SELECT CASE ( var )
      CASE ( VAR_T )   ; base = Atm(1)%Temperature
      CASE ( VAR_H2O ) ; base = Atm(1)%Absorber(:,IDX_H2O)
      CASE ( VAR_NO2 ) ; base = Atm(1)%Absorber(:,IDX_NO2)
    END SELECT
  END SUBROUTINE get_base

  ! ----------------------------------------------------------------
  ! Check 1 : TL vs central finite difference for a whole-column
  !           perturbation of variable `var` on profile 1. The TL
  !           direction matches the FD perturbation shape, so the
  !           FD/TL ratio must -> 1. A variable with no forward
  !           response (blind against the loaded file) must show a
  !           consistently zero TL instead.
  ! ----------------------------------------------------------------
  SUBROUTINE check_fd( var, ok )
    INTEGER, INTENT(IN)  :: var
    LOGICAL, INTENT(OUT) :: ok
    CHARACTER(16) :: vname
    REAL(fp) :: base(N_LAYERS), fd_all(n_Channels)
    REAL(fp) :: tl, fd, Rp, Rm, ratio, best, eps, tl_max
    INTEGER  :: ii, kk, ch

    SELECT CASE ( var )
      CASE ( VAR_T )   ; vname = 'Temperature'
      CASE ( VAR_H2O ) ; vname = 'H2O column'
      CASE ( VAR_NO2 ) ; vname = 'NO2 column'
    END SELECT
    CALL get_base( var, base )

    ! TL along the same direction as the FD perturbation
    CALL CRTM_Atmosphere_Zero( Atm_TL ) ; CALL CRTM_Surface_Zero( Sfc_TL )
    SELECT CASE ( var )
      CASE ( VAR_T )   ; Atm_TL(1)%Temperature        = ONE
      CASE ( VAR_H2O ) ; Atm_TL(1)%Absorber(:,IDX_H2O) = base
      CASE ( VAR_NO2 ) ; Atm_TL(1)%Absorber(:,IDX_NO2) = base
    END SELECT
    Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, ChannelInfo, &
                                        RTSolution, RTSolution_TL )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'TL fail' ; ok=.FALSE. ; RETURN ; END IF

    ! Channel selection by FD probe (not max|TL|: a broken zero TL must not
    ! hide itself).
    eps = 1.0e-3_fp
    CALL perturb( var, base, +eps )
    Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution_pert )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'FD forward fail' ; ok=.FALSE. ; RETURN ; END IF
    DO ii = 1, n_Channels ; fd_all(ii) = RTSolution_pert(ii,1)%Radiance ; END DO
    CALL perturb( var, base, -eps )
    Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution_pert )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'FD forward fail' ; ok=.FALSE. ; RETURN ; END IF
    DO ii = 1, n_Channels
      fd_all(ii) = ( fd_all(ii) - RTSolution_pert(ii,1)%Radiance ) / ( 2.0_fp*eps )
    END DO
    CALL perturb( var, base, ZERO )

    ! Blind-variable detection: a variable the loaded file cannot see (NO2
    ! against a plain group-2 file; H2O against the zero-base parity file)
    ! must show ZERO forward response AND zero TL -- that is the
    ! consistency/backward-compatibility contract.
    IF ( MAXVAL(ABS(fd_all)) <= FD_ZERO ) THEN
      tl_max = ZERO
      DO ii = 1, n_Channels
        tl_max = MAX( tl_max, ABS(RTSolution_TL(ii,1)%Radiance) )
      END DO
      ok = ( tl_max <= FD_ZERO )
      WRITE(*,'(/7x,"[FD] ",a,": no forward response (blind for this file).  max|TL| = ",es11.4)') &
            TRIM(vname), tl_max
      WRITE(*,'(7x,"-> TL consistently zero : ",a)') MERGE('PASS','FAIL',ok)
      IF ( var == VAR_NO2 ) THEN
        no2_active = .FALSE.
        ch_no2 = 0
      END IF
      RETURN
    END IF
    IF ( var == VAR_NO2 ) THEN
      no2_active = .TRUE.
      ch_no2 = MAXLOC( ABS(fd_all), DIM=1 )
    END IF

    ch = MAXLOC( ABS(fd_all), DIM=1 )
    tl = RTSolution_TL(ch,1)%Radiance
    IF ( ABS(tl) <= TINY(ONE) ) THEN
      WRITE(*,'(/7x,"[FD] d Radiance / d ",a,"   channel ",i0,": TL is ZERO but FD=",es13.6)') &
            TRIM(vname), RTSolution(ch,1)%Sensor_Channel, fd_all(ch)
      ok = .FALSE.
      RETURN
    END IF

    best = HUGE(ONE)
    WRITE(*,'(/7x,"[FD] d Radiance / d ",a,"   channel ",i0,"   TL=",es13.6)') &
          TRIM(vname), RTSolution(ch,1)%Sensor_Channel, tl
    DO kk = 4, 14
      eps = 0.1_fp / (2.0_fp**kk)
      CALL perturb( var, base, +eps )
      Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution_pert )
      IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'FD forward fail' ; ok=.FALSE. ; RETURN ; END IF
      Rp = RTSolution_pert(ch,1)%Radiance
      CALL perturb( var, base, -eps )
      Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution_pert )
      IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'FD forward fail' ; ok=.FALSE. ; RETURN ; END IF
      Rm = RTSolution_pert(ch,1)%Radiance
      CALL perturb( var, base, ZERO )
      fd = ( Rp - Rm ) / ( 2.0_fp*eps )
      ratio = fd / tl
      IF ( ABS(ratio-ONE) < best ) best = ABS(ratio-ONE)
      WRITE(*,'(9x,"eps=",es10.3,"  FD=",es16.9,"  FD/TL=",f14.10)') eps, fd, ratio
    END DO
    ok = ( best < TOL_FD )
    WRITE(*,'(7x,"-> best |FD/TL - 1| = ",es11.4,"   ",a)') best, MERGE('PASS','FAIL',ok)
  END SUBROUTINE check_fd

  ! ----------------------------------------------------------------
  ! Check 2 : adjoint dot-product  <dy,dy> == <x, AD(dy)>,
  !           x spanning Temperature + H2O + NO2 everywhere.
  ! ----------------------------------------------------------------
  SUBROUTINE check_adj( ok )
    LOGICAL, INTENT(OUT) :: ok
    REAL(fp) :: LHS, RHS, dy, rel_adj
    INTEGER  :: ii, mm

    CALL CRTM_Atmosphere_Zero( Atm_TL ) ; CALL CRTM_Surface_Zero( Sfc_TL )
    DO mm = 1, N_PROFILES
      DO ii = 1, N_LAYERS
        Atm_TL(mm)%Temperature(ii)      = 0.5_fp * SIN( 0.7_fp*REAL(ii,fp) + 1.3_fp*REAL(mm,fp) )
        Atm_TL(mm)%Absorber(ii,IDX_H2O) = 0.05_fp * Atm(mm)%Absorber(ii,IDX_H2O) &
                                          * COS( 0.9_fp*REAL(ii,fp) + 0.4_fp*REAL(mm,fp) )
        Atm_TL(mm)%Absorber(ii,IDX_NO2) = 0.05_fp * Atm(mm)%Absorber(ii,IDX_NO2) &
                                          * SIN( 1.1_fp*REAL(ii,fp) + 0.8_fp*REAL(mm,fp) )
      END DO
    END DO
    Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, ChannelInfo, &
                                        RTSolution, RTSolution_TL )
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
                                 Atm_AD, Sfc_AD, RTSolution )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'AD fail' ; ok=.FALSE. ; RETURN ; END IF

    RHS = ZERO
    DO mm = 1, N_PROFILES
      DO ii = 1, N_LAYERS
        RHS = RHS + Atm_TL(mm)%Temperature(ii)      * Atm_AD(mm)%Temperature(ii)
        RHS = RHS + Atm_TL(mm)%Absorber(ii,IDX_H2O) * Atm_AD(mm)%Absorber(ii,IDX_H2O)
        RHS = RHS + Atm_TL(mm)%Absorber(ii,IDX_NO2) * Atm_AD(mm)%Absorber(ii,IDX_NO2)
      END DO
    END DO
    rel_adj = ABS(LHS-RHS) / MAX(ABS(LHS), TINY(ONE))
    ok = ( rel_adj < TOL_ADJ )
    WRITE(*,'(/7x,"[ADJ] <dy,dy>=",es16.9,"   <x,gx>=",es16.9)') LHS, RHS
    WRITE(*,'(7x,"-> relative difference = ",es11.4,"   ",a)') rel_adj, MERGE('PASS','FAIL',ok)
  END SUBROUTINE check_adj

  ! ----------------------------------------------------------------
  ! Check 2b: adjoint dot-product with x spanning NO2 ONLY, so the
  !           dy signal is generated entirely by the new group-8
  !           chain (the RT legs are still traversed, so the bound
  !           is the same solar-RT roundoff as the combined check).
  ! ----------------------------------------------------------------
  SUBROUTINE check_adj_no2( ok )
    LOGICAL, INTENT(OUT) :: ok
    REAL(fp) :: LHS, RHS, dy, rel_adj
    INTEGER  :: ii, mm

    CALL CRTM_Atmosphere_Zero( Atm_TL ) ; CALL CRTM_Surface_Zero( Sfc_TL )
    DO mm = 1, N_PROFILES
      DO ii = 1, N_LAYERS
        Atm_TL(mm)%Absorber(ii,IDX_NO2) = 0.05_fp * Atm(mm)%Absorber(ii,IDX_NO2) &
                                          * SIN( 1.1_fp*REAL(ii,fp) + 0.8_fp*REAL(mm,fp) )
      END DO
    END DO
    Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, ChannelInfo, &
                                        RTSolution, RTSolution_TL )
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
                                 Atm_AD, Sfc_AD, RTSolution )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'AD fail' ; ok=.FALSE. ; RETURN ; END IF

    RHS = ZERO
    DO mm = 1, N_PROFILES
      DO ii = 1, N_LAYERS
        RHS = RHS + Atm_TL(mm)%Absorber(ii,IDX_NO2) * Atm_AD(mm)%Absorber(ii,IDX_NO2)
      END DO
    END DO
    ! An NO2-blind file gives LHS = RHS = 0 exactly; that degenerate pass is
    ! the backward-compatibility control.
    rel_adj = ABS(LHS-RHS) / MAX(ABS(LHS), TINY(ONE))
    ok = ( rel_adj < TOL_ADJ_NO2 )
    WRITE(*,'(/7x,"[ADJ-NO2] <dy,dy>=",es16.9,"   <x,gx>=",es16.9)') LHS, RHS
    WRITE(*,'(7x,"-> relative difference = ",es11.4,"   ",a)') rel_adj, MERGE('PASS','FAIL',ok)
  END SUBROUTINE check_adj_no2

  ! ----------------------------------------------------------------
  ! Check 3 : K-Matrix vs Adjoint Jacobian (Temperature, H2O and NO2
  !           columns) on the most NO2-sensitive channel (channel 1
  !           for an NO2-blind file).
  ! ----------------------------------------------------------------
  SUBROUTINE check_k( ok )
    LOGICAL, INTENT(OUT) :: ok
    REAL(fp) :: maxdiff, scal, rel_k
    INTEGER  :: l0, m0, mm

    l0 = MAX( ch_no2, 1 ) ; m0 = 1

    CALL CRTM_Atmosphere_Zero( Atm_K ) ; CALL CRTM_Surface_Zero( Sfc_K )
    CALL CRTM_RTSolution_Zero( RTSolution_K )
    DO mm = 1, N_PROFILES
      DO l = 1, n_Channels
        RTSolution_K(l,mm)%Radiance = ONE
      END DO
    END DO
    Error_Status = CRTM_K_Matrix( Atm, Sfc, RTSolution_K, Geometry, ChannelInfo, &
                                  Atm_K, Sfc_K, RTSolution )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'K fail' ; ok=.FALSE. ; RETURN ; END IF

    CALL CRTM_Atmosphere_Zero( Atm_AD ) ; CALL CRTM_Surface_Zero( Sfc_AD )
    CALL CRTM_RTSolution_Zero( RTSolution_AD )
    RTSolution_AD(l0,m0)%Radiance = ONE
    Error_Status = CRTM_Adjoint( Atm, Sfc, RTSolution_AD, Geometry, ChannelInfo, &
                                 Atm_AD, Sfc_AD, RTSolution )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'AD fail' ; ok=.FALSE. ; RETURN ; END IF

    maxdiff = MAX( MAXVAL( ABS( Atm_K(l0,m0)%Temperature - Atm_AD(m0)%Temperature ) ),          &
                   MAXVAL( ABS( Atm_K(l0,m0)%Absorber(:,IDX_H2O) - Atm_AD(m0)%Absorber(:,IDX_H2O) ) ), &
                   MAXVAL( ABS( Atm_K(l0,m0)%Absorber(:,IDX_NO2) - Atm_AD(m0)%Absorber(:,IDX_NO2) ) ) )
    scal    = MAX( MAXVAL(ABS(Atm_K(l0,m0)%Temperature)),          &
                   MAXVAL(ABS(Atm_K(l0,m0)%Absorber(:,IDX_H2O))),  &
                   MAXVAL(ABS(Atm_K(l0,m0)%Absorber(:,IDX_NO2))), TINY(ONE) )
    rel_k   = maxdiff / scal
    ok = ( rel_k < TOL_K )
    WRITE(*,'(/7x,"[K] K vs AD (channel ",i0,"):  max|K-AD|/max|K| = ",es11.4,"   ",a)') &
          RTSolution(l0,m0)%Sensor_Channel, rel_k, MERGE('PASS','FAIL',ok)
    IF ( no2_active ) THEN
      WRITE(*,'(7x,"max|dR/dNO2| (K, channel ",i0,") = ",es11.4)') &
            RTSolution(l0,m0)%Sensor_Channel, MAXVAL(ABS(Atm_K(l0,m0)%Absorber(:,IDX_NO2)))
    END IF
  END SUBROUTINE check_k

  ! UARS NO2 climatology (Kerr/Louisnard, NO2_vmr_UARS.txt) log-P interpolated
  ! onto the ECMWF84 100-layer grid; ppmv. Same profile family as the
  ! Ref_Absorber written into the group-8 TauCoeff.
  SUBROUTINE Set_NO2_Climatology()
    no2_clim = (/ &
      3.793332e-06_fp, 2.282579e-05_fp, 7.604145e-05_fp, 2.653026e-04_fp, 1.458149e-03_fp,  &
      2.313445e-03_fp, 2.860727e-03_fp, 3.223867e-03_fp, 3.684018e-03_fp, 4.757346e-03_fp,  &
      5.907983e-03_fp, 7.026512e-03_fp, 6.319363e-03_fp, 3.891646e-03_fp, 2.760278e-03_fp,  &
      1.903978e-03_fp, 1.352134e-03_fp, 1.132655e-03_fp, 1.005111e-03_fp, 9.050189e-04_fp,  &
      7.515740e-04_fp, 6.018733e-04_fp, 5.057622e-04_fp, 4.186505e-04_fp, 3.793749e-04_fp,  &
      3.669254e-04_fp, 3.555738e-04_fp, 3.496207e-04_fp, 3.438747e-04_fp, 3.385942e-04_fp,  &
      3.340939e-04_fp, 3.297371e-04_fp, 3.252582e-04_fp, 3.197292e-04_fp, 3.143633e-04_fp,  &
      3.091517e-04_fp, 3.061413e-04_fp, 3.041805e-04_fp, 3.022723e-04_fp, 3.004143e-04_fp,  &
      3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp,  &
      3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp,  &
      3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp,  &
      3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp,  &
      3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp,  &
      3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp,  &
      3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp,  &
      3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp,  &
      3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp,  &
      3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp,  &
      3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp,  &
      3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp, 3.000000e-04_fp /)
  END SUBROUTINE Set_NO2_Climatology

  INCLUDE 'Load_ECMWF84_Atm_Data.inc'

END PROGRAM test_UV_NO2_TLAD
