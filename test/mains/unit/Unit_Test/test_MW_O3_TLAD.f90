!
! test_MW_O3_TLAD
!
! TL/AD/K parity check for the MW scene-ozone ODPS component
! (GROUP_MW_O3, Group_Index = 7).
!
! The group-7 ozone predictor blocks in ODPS_Predictor.f90 are transcriptions
! of the validated IR ozone formulation; the forward path has been verified
! (zero-coefficient plumbing identity + synthetic single-channel isolation)
! but no test exercises their TL/AD/K consistency. This test initializes
! mwr_aws on a clear-sky ECMWF84 ocean column and verifies:
!   1. TL vs central finite difference for column perturbations of
!      O3 (relative), H2O (relative) and Temperature (additive) -- the O3
!      check probes the new predictor block end-to-end.
!   2. Adjoint dot-product  <dy,dy> == <x, AD(dy)>  with x spanning
!      Temperature, H2O and O3 on every layer of every profile.
!   3. K-Matrix vs Adjoint Jacobian equality (Temperature, H2O and O3
!      columns) on the most O3-sensitive channel.
!
! The test adapts to the loaded TauCoeff: with a 2-component group-3 file
! (no ozone absorber) the O3 response must be IDENTICALLY ZERO in both the
! forward FD probe and the TL -- that run is the backward-compatibility
! control. With a group-7 file the O3 TL must converge to the FD derivative.
!
! Exit: STOP 0 if every check passes, STOP 1 otherwise.
!
! CREATION HISTORY:
!       Written by:     Benjamin Johnson, 15-Jul-2026
!                       Setup and verification machinery adapted from
!                       test_VectorRT_TLADK.
!
PROGRAM test_MW_O3_TLAD

  USE CRTM_Module
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_MW_O3_TLAD'
  CHARACTER(*), PARAMETER :: PATH   = './testinput/'
  CHARACTER(*), PARAMETER :: SENSOR = 'mwr_aws'

  ! Profile / column setup (ECMWF84 ocean column; clear sky)
  INTEGER,  PARAMETER :: N_PROFILES  = 2
  INTEGER,  PARAMETER :: N_LAYERS    = 100
  INTEGER,  PARAMETER :: N_ABSORBERS = 6
  INTEGER,  PARAMETER :: N_CLOUDS    = 0
  INTEGER,  PARAMETER :: N_AEROSOLS  = 0
  REAL(fp), PARAMETER :: ZENITH      = 53.0_fp
  INTEGER,  PARAMETER :: IDX_H2O = 1, IDX_O3 = 3   ! ECMWF84 absorber slots

  REAL(fp), PARAMETER :: TOL_FD   = 1.0e-3_fp     ! TL vs finite difference
  REAL(fp), PARAMETER :: TOL_ADJ  = 1.0e-12_fp    ! adjoint dot-product
  REAL(fp), PARAMETER :: TOL_K    = 1.0e-9_fp     ! K vs AD
  ! Forward radiances are O(1e0..1e2) mW/(m2.sr.cm-1); a column FD response
  ! below this is numerically indistinguishable from zero -> O3-blind file.
  REAL(fp), PARAMETER :: FD_ZERO  = 1.0e-9_fp

  ! Perturbation-variable selectors
  INTEGER, PARAMETER :: VAR_T = 1, VAR_H2O = 2, VAR_O3 = 3

  CHARACTER(256) :: Version
  INTEGER :: Error_Status, Allocate_Status, n_Channels
  INTEGER :: l, m
  INTEGER :: ch_o3                 ! most O3-sensitive channel (0 if none)
  LOGICAL :: o3_active             ! loaded file responds to scene O3
  LOGICAL :: ok_fd_o3, ok_fd_h2o, ok_fd_t, ok_adj, ok_k

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
  WRITE(*,'(/5x,a)') 'MW scene-ozone (GROUP_MW_O3) TL/AD/K verification'
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

  ! Base clear-sky column for every profile; second profile gets a scaled O3
  ! column so both profiles contribute distinct ozone signals to the adjoint
  ! dot-product.
  CALL Load_ECMWF84_Atm_Data()         ! fills Atm(1) and Atm(2)
  Atm(2)%Absorber(:,IDX_O3) = 1.3_fp * Atm(2)%Absorber(:,IDX_O3)

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
  ! Checks. check_fd(VAR_O3) also determines o3_active/ch_o3 from the forward
  ! FD probe, so it must run first.
  ! --------------------------------------------------------------------------
  CALL check_fd ( VAR_O3 , ok_fd_o3  )
  CALL check_fd ( VAR_H2O, ok_fd_h2o )
  CALL check_fd ( VAR_T  , ok_fd_t   )
  CALL check_adj( ok_adj )
  CALL check_k  ( ok_k )

  Error_Status = CRTM_Destroy( ChannelInfo )

  WRITE(*,'(/5x,a)') '====================================================='
  IF ( o3_active ) THEN
    WRITE(*,'(5x,a)') 'TauCoeff mode : scene-O3 ACTIVE (group-7 ozone component)'
  ELSE
    WRITE(*,'(5x,a)') 'TauCoeff mode : O3-blind control (2-component MW file)'
  END IF
  WRITE(*,'(5x,"TL vs FD (O3 column)                  : ",a)') MERGE('PASS','FAIL',ok_fd_o3)
  WRITE(*,'(5x,"TL vs FD (H2O column)                 : ",a)') MERGE('PASS','FAIL',ok_fd_h2o)
  WRITE(*,'(5x,"TL vs FD (Temperature)                : ",a)') MERGE('PASS','FAIL',ok_fd_t)
  WRITE(*,'(5x,"adjoint dot-product (T+H2O+O3)        : ",a)') MERGE('PASS','FAIL',ok_adj)
  WRITE(*,'(5x,"K vs AD (T/H2O/O3 columns)            : ",a)') MERGE('PASS','FAIL',ok_k)
  IF ( ok_fd_o3 .AND. ok_fd_h2o .AND. ok_fd_t .AND. ok_adj .AND. ok_k ) THEN
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
      CASE ( VAR_O3 )  ; Atm(1)%Absorber(:,IDX_O3)  = base * (ONE + eps)
    END SELECT
  END SUBROUTINE perturb

  SUBROUTINE get_base( var, base )
    INTEGER,  INTENT(IN)  :: var
    REAL(fp), INTENT(OUT) :: base(N_LAYERS)
    SELECT CASE ( var )
      CASE ( VAR_T )   ; base = Atm(1)%Temperature
      CASE ( VAR_H2O ) ; base = Atm(1)%Absorber(:,IDX_H2O)
      CASE ( VAR_O3 )  ; base = Atm(1)%Absorber(:,IDX_O3)
    END SELECT
  END SUBROUTINE get_base

  ! ----------------------------------------------------------------
  ! Check 1 : TL vs central finite difference for a whole-column
  !           perturbation of variable `var` on profile 1. The TL
  !           direction matches the FD perturbation shape, so the
  !           FD/TL ratio must -> 1.
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
      CASE ( VAR_O3 )  ; vname = 'O3 column'
    END SELECT
    CALL get_base( var, base )

    ! TL along the same direction as the FD perturbation
    CALL CRTM_Atmosphere_Zero( Atm_TL ) ; CALL CRTM_Surface_Zero( Sfc_TL )
    SELECT CASE ( var )
      CASE ( VAR_T )   ; Atm_TL(1)%Temperature        = ONE
      CASE ( VAR_H2O ) ; Atm_TL(1)%Absorber(:,IDX_H2O) = base
      CASE ( VAR_O3 )  ; Atm_TL(1)%Absorber(:,IDX_O3)  = base
    END SELECT
    Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, ChannelInfo, &
                                        RTSolution, RTSolution_TL )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'TL fail' ; ok=.FALSE. ; RETURN ; END IF

    ! Channel selection by FD probe (not max|TL|: a broken zero TL must not
    ! hide itself).
    eps = 1.0e-3_fp
    CALL perturb( var, base, +eps )
    Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution_pert )
    DO ii = 1, n_Channels ; fd_all(ii) = RTSolution_pert(ii,1)%Radiance ; END DO
    CALL perturb( var, base, -eps )
    Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution_pert )
    DO ii = 1, n_Channels
      fd_all(ii) = ( fd_all(ii) - RTSolution_pert(ii,1)%Radiance ) / ( 2.0_fp*eps )
    END DO
    CALL perturb( var, base, ZERO )

    ! O3 mode detection: an O3-blind (group-3) file must show ZERO forward
    ! response AND zero TL to an O3 perturbation -- that is the
    ! backward-compatibility contract.
    IF ( var == VAR_O3 ) THEN
      tl_max = ZERO
      DO ii = 1, n_Channels
        tl_max = MAX( tl_max, ABS(RTSolution_TL(ii,1)%Radiance) )
      END DO
      o3_active = ( MAXVAL(ABS(fd_all)) > FD_ZERO )
      ch_o3 = MAXLOC( ABS(fd_all), DIM=1 )
      IF ( .NOT. o3_active ) THEN
        ok = ( tl_max <= FD_ZERO )
        WRITE(*,'(/7x,"[FD] O3: no forward response (O3-blind file).  max|TL| = ",es11.4)') tl_max
        WRITE(*,'(7x,"-> TL consistently zero : ",a)') MERGE('PASS','FAIL',ok)
        ch_o3 = 0
        RETURN
      END IF
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
      Rp = RTSolution_pert(ch,1)%Radiance
      CALL perturb( var, base, -eps )
      Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution_pert )
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
  !           x spanning Temperature + H2O + O3 everywhere.
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
        Atm_TL(mm)%Absorber(ii,IDX_O3)  = 0.05_fp * Atm(mm)%Absorber(ii,IDX_O3) &
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
        RHS = RHS + Atm_TL(mm)%Absorber(ii,IDX_O3)  * Atm_AD(mm)%Absorber(ii,IDX_O3)
      END DO
    END DO
    rel_adj = ABS(LHS-RHS) / MAX(ABS(LHS), TINY(ONE))
    ok = ( rel_adj < TOL_ADJ )
    WRITE(*,'(/7x,"[ADJ] <dy,dy>=",es16.9,"   <x,gx>=",es16.9)') LHS, RHS
    WRITE(*,'(7x,"-> relative difference = ",es11.4,"   ",a)') rel_adj, MERGE('PASS','FAIL',ok)
  END SUBROUTINE check_adj

  ! ----------------------------------------------------------------
  ! Check 3 : K-Matrix vs Adjoint Jacobian (Temperature, H2O and O3
  !           columns) on the most O3-sensitive channel (channel 1
  !           for an O3-blind file).
  ! ----------------------------------------------------------------
  SUBROUTINE check_k( ok )
    LOGICAL, INTENT(OUT) :: ok
    REAL(fp) :: maxdiff, scal, rel_k
    INTEGER  :: l0, m0, mm

    l0 = MAX( ch_o3, 1 ) ; m0 = 1

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
                   MAXVAL( ABS( Atm_K(l0,m0)%Absorber(:,IDX_O3)  - Atm_AD(m0)%Absorber(:,IDX_O3)  ) ) )
    scal    = MAX( MAXVAL(ABS(Atm_K(l0,m0)%Temperature)),          &
                   MAXVAL(ABS(Atm_K(l0,m0)%Absorber(:,IDX_H2O))),  &
                   MAXVAL(ABS(Atm_K(l0,m0)%Absorber(:,IDX_O3))), TINY(ONE) )
    rel_k   = maxdiff / scal
    ok = ( rel_k < TOL_K )
    WRITE(*,'(/7x,"[K] K vs AD (channel ",i0,"):  max|K-AD|/max|K| = ",es11.4,"   ",a)') &
          RTSolution(l0,m0)%Sensor_Channel, rel_k, MERGE('PASS','FAIL',ok)
    IF ( o3_active ) THEN
      WRITE(*,'(7x,"max|dR/dO3| (K, channel ",i0,") = ",es11.4)') &
            RTSolution(l0,m0)%Sensor_Channel, MAXVAL(ABS(Atm_K(l0,m0)%Absorber(:,IDX_O3)))
    END IF
  END SUBROUTINE check_k

  INCLUDE 'Load_ECMWF84_Atm_Data.inc'

END PROGRAM test_MW_O3_TLAD
