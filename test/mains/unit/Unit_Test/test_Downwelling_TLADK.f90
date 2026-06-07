!
! test_Downwelling_TLADK
!
! Baseline-independent correctness check for the always-on surface downwelling
! radiance output (RTSolution%Down_Radiance) in the Tangent-Linear, Adjoint and
! K-Matrix models.
!
! Downwelling radiance at the surface is a first-class, always-computed output.
! This test verifies its Jacobians without relying on stored baselines, for both
! the standard TOA upwelling radiance (control) and the surface Down_Radiance:
!   1. TL vs central finite-difference of the forward model   (TL = dF/dx ?)
!   2. Adjoint dot-product test  <TL.x,TL.x> == <x,AD.(TL.x)> (AD = TL^T ?)
!   3. K-Matrix vs Adjoint Jacobian equality                  (K wiring ok ?)
!
! TOA control is seeded via RTSolution%Radiance; downwelling via %Down_Radiance.
! Exit: STOP 0 if every check passes, STOP 1 otherwise.
!
PROGRAM test_Downwelling_TLADK

  USE CRTM_Module
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME      = 'test_Downwelling_TLADK'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'

  ! Clear-sky profile/sensor setup (Emission solver path)
  INTEGER,  PARAMETER :: N_PROFILES  = 2
  INTEGER,  PARAMETER :: N_LAYERS    = 92
  INTEGER,  PARAMETER :: N_ABSORBERS = 2
  INTEGER,  PARAMETER :: N_CLOUDS    = 1    ! allocate a cloud slot; content toggled per scene
  INTEGER,  PARAMETER :: N_AEROSOLS  = 0
  INTEGER,  PARAMETER :: N_SENSORS   = 1
  REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp
  INTEGER,  PARAMETER :: PERT_LAYER = 85          ! perturbed temperature layer (within the low cloud)

  REAL(fp), PARAMETER :: TOL_FD  = 1.0e-3_fp      ! TL vs finite difference
  REAL(fp), PARAMETER :: TOL_ADJ = 1.0e-9_fp      ! adjoint dot-product
  REAL(fp), PARAMETER :: TOL_K   = 1.0e-9_fp      ! K vs AD

  CHARACTER(256) :: Version, Sensor_Id
  INTEGER :: Error_Status, Allocate_Status, n_Channels
  INTEGER :: l, m
  LOGICAL :: ok_toa, ok_dwn, ok_toa_s, ok_dwn_s, ok_toa_a, ok_dwn_a

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(N_SENSORS)
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
  WRITE(*,'(/5x,a)') 'Downwelling (surface Down_Radiance) TL/AD/K verification'
  WRITE(*,'(5x,a/)') 'CRTM Version: '//TRIM(Version)

  Sensor_Id = 'atms_npp'

  Error_Status = CRTM_Init( (/Sensor_Id/), ChannelInfo, File_Path=COEFFICIENTS_PATH )
  IF ( Error_Status /= SUCCESS ) THEN
    WRITE(*,*) 'Error initializing CRTM'; STOP 1
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))

  ALLOCATE( RTSolution(n_Channels,N_PROFILES), RTSolution_pert(n_Channels,N_PROFILES), &
            RTSolution_TL(n_Channels,N_PROFILES), RTSolution_AD(n_Channels,N_PROFILES), &
            RTSolution_K(n_Channels,N_PROFILES), &
            Atm_K(n_Channels,N_PROFILES), Sfc_K(n_Channels,N_PROFILES), &
            STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN; WRITE(*,*) 'Alloc error'; STOP 1; END IF

  CALL CRTM_Atmosphere_Create( Atm,    N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atm_TL, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atm_AD, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atm_K,  N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )

  CALL Load_Atm_Data()
  CALL Load_Sfc_Data()

  ! Make the TL/AD/K input atmospheres layer-independently congruent with the FWD
  ! atmosphere (required by the fractional-cloud ClearSkyCopy path).
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

  CALL CRTM_Geometry_SetValue( Geometry, &
                               Sensor_Zenith_Angle = ZENITH_ANGLE, &
                               Sensor_Scan_Angle   = SCAN_ANGLE )

  CALL verify( .FALSE., .FALSE., RT_ADA, ZERO,   ok_toa   )  ! clear-sky: TOA radiance (control)
  CALL verify( .TRUE. , .FALSE., RT_ADA, ZERO,   ok_dwn   )  ! clear-sky: surface Down_Radiance (emission)
  CALL verify( .FALSE., .TRUE. , RT_SOI, 0.5_fp, ok_toa_s )  ! SOI scattering (fractional: + combine)
  CALL verify( .TRUE. , .TRUE. , RT_SOI, 0.5_fp, ok_dwn_s )  ! SOI Down_Radiance (fractional: + combine)
  CALL verify( .FALSE., .TRUE. , RT_ADA, ONE,    ok_toa_a )  ! ADA scattering (overcast: isolate solver)
  CALL verify( .TRUE. , .TRUE. , RT_ADA, ONE,    ok_dwn_a )  ! ADA Down_Radiance (overcast: isolate solver)

  Error_Status = CRTM_Destroy( ChannelInfo )

  WRITE(*,'(/5x,a)') '====================================================='
  WRITE(*,'(5x,"clear-sky TOA control       : ",a)') MERGE('PASS','FAIL',ok_toa)
  WRITE(*,'(5x,"clear-sky Down_Radiance     : ",a)') MERGE('PASS','FAIL',ok_dwn)
  WRITE(*,'(5x,"SOI scattering TOA control  : ",a)') MERGE('PASS','FAIL',ok_toa_s)
  WRITE(*,'(5x,"SOI scattering Down_Radiance: ",a)') MERGE('PASS','FAIL',ok_dwn_s)
  WRITE(*,'(5x,"ADA scattering TOA control  : ",a)') MERGE('PASS','FAIL',ok_toa_a)
  WRITE(*,'(5x,"ADA scattering Down_Radiance: ",a)') MERGE('PASS','FAIL',ok_dwn_a)
  IF ( ok_toa .AND. ok_dwn .AND. ok_toa_s .AND. ok_dwn_s .AND. ok_toa_a .AND. ok_dwn_a ) THEN
    WRITE(*,'(5x,a)') 'ALL CHECKS PASSED'
    STOP 0
  ELSE
    WRITE(*,'(5x,a)') 'CHECKS FAILED'
    STOP 1
  END IF

CONTAINS

  ! Read the selected output (TOA radiance or surface downwelling) from an RTSolution
  REAL(fp) FUNCTION get_out( rts, downwelling )
    TYPE(CRTM_RTSolution_type), INTENT(IN) :: rts
    LOGICAL, INTENT(IN) :: downwelling
    IF ( downwelling ) THEN ; get_out = rts%Down_Radiance ; ELSE ; get_out = rts%Radiance ; END IF
  END FUNCTION get_out

  ! Seed the selected adjoint/K output to a value
  SUBROUTINE set_seed( rts, downwelling, val )
    TYPE(CRTM_RTSolution_type), INTENT(INOUT) :: rts
    LOGICAL,  INTENT(IN) :: downwelling
    REAL(fp), INTENT(IN) :: val
    IF ( downwelling ) THEN ; rts%Down_Radiance = val ; ELSE ; rts%Radiance = val ; END IF
  END SUBROUTINE set_seed

  SUBROUTINE verify( downwelling, scattering, rt_alg, cfrac, all_ok )
    LOGICAL,  INTENT(IN)  :: downwelling, scattering
    INTEGER,  INTENT(IN)  :: rt_alg
    REAL(fp), INTENT(IN)  :: cfrac
    LOGICAL,  INTENT(OUT) :: all_ok
    CHARACTER(48) :: tag
    REAL(fp) :: tl, fd, R0, Rp, Rm, ratio, best, delta, T0
    REAL(fp) :: LHS, RHS, rel_adj, dy
    REAL(fp) :: maxdiff, scal, rel_k
    REAL(fp) :: fd_all(n_Channels)
    INTEGER  :: ii, kk, ch, l0, m0, nscat
    LOGICAL  :: ok1, ok2, ok3

    ! Configure the scene: clear-sky (emission) or SOI scattering (opt-in downwelling)
    DO m = 1, N_PROFILES
      IF ( scattering ) THEN
        Options(m)%RT_Algorithm_Id       = rt_alg
        Options(m)%Compute_Down_Radiance = .TRUE.
        Atm(m)%n_Clouds                  = 1
        ! Thick low cloud so the cloudy (scattering) contribution dominates the
        ! SURFACE downwelling for the sensitive channels (not masked by clear emission).
        Atm(m)%Cloud_Fraction            = ZERO
        Atm(m)%Cloud_Fraction(70:90)     = cfrac  ! cfrac=1 overcast (isolate solver); 0<cfrac<1 exercises combine
        Atm(m)%Cloud(1)%Type             = SNOW_CLOUD   ! mm-size particles -> strong MW scattering
        Atm(m)%Cloud(1)%Effective_Radius = ZERO
        Atm(m)%Cloud(1)%Water_Content    = ZERO
        Atm(m)%Cloud(1)%Effective_Radius(70:90) = 500.0_fp
        Atm(m)%Cloud(1)%Water_Content(70:90)    = 5.0_fp
      ELSE
        Options(m)%RT_Algorithm_Id       = RT_ADA
        Options(m)%Compute_Down_Radiance = .FALSE.
        Atm(m)%n_Clouds                  = 0
        Atm(m)%Cloud(1)%Water_Content    = ZERO
      END IF
    END DO

    IF ( downwelling ) THEN ; tag = 'Down_Radiance' ; ELSE ; tag = 'TOA Radiance' ; END IF
    IF ( scattering ) THEN
      IF ( rt_alg == RT_SOI ) THEN ; tag = TRIM(tag)//' [SOI scattering]' ; ELSE ; tag = TRIM(tag)//' [ADA scattering]' ; END IF
    ELSE ; tag = TRIM(tag)//' [clear-sky]' ; END IF
    WRITE(*,'(/5x,"============== Output: ",a," ==============")') TRIM(tag)

    ! ----------------------------------------------------------------
    ! Check 1 : TL vs central finite-difference of the forward model
    ! ----------------------------------------------------------------
    CALL CRTM_Atmosphere_Zero( Atm_TL ) ; CALL CRTM_Surface_Zero( Sfc_TL )
    Atm_TL(1)%Temperature(PERT_LAYER) = ONE
    Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, ChannelInfo, &
                                        RTSolution, RTSolution_TL, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'TL fail' ; all_ok=.FALSE. ; RETURN ; END IF

    IF ( scattering ) THEN
      nscat = COUNT( RTSolution(:,1)%Scattering_Flag )
      WRITE(*,'(7x,"(scattering channels in profile 1: ",i0," of ",i0,")")') nscat, n_Channels
    END IF

    ! Pick the channel most sensitive to T(PERT_LAYER) via a finite-difference PROBE,
    ! NOT max|TL|: a broken (zero) TL would otherwise hide itself from selection.
    T0 = Atm(1)%Temperature(PERT_LAYER)
    delta = ABS(T0) * 0.1_fp / 256.0_fp
    Atm(1)%Temperature(PERT_LAYER) = T0 + delta
    Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution_pert, Options=Options )
    DO ii = 1, n_Channels ; fd_all(ii) = get_out(RTSolution_pert(ii,1),downwelling) ; END DO
    Atm(1)%Temperature(PERT_LAYER) = T0 - delta
    Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution_pert, Options=Options )
    DO ii = 1, n_Channels
      fd_all(ii) = ( fd_all(ii) - get_out(RTSolution_pert(ii,1),downwelling) ) / ( 2.0_fp*delta )
    END DO
    Atm(1)%Temperature(PERT_LAYER) = T0
    ch = 0 ; best = ZERO
    DO ii = 1, n_Channels
      IF ( scattering .AND. .NOT. RTSolution(ii,1)%Scattering_Flag ) CYCLE
      IF ( ABS(fd_all(ii)) >= best ) THEN ; best = ABS(fd_all(ii)) ; ch = ii ; END IF
    END DO
    IF ( ch == 0 ) ch = 1
    tl = get_out(RTSolution_TL(ch,1),downwelling)

    T0 = Atm(1)%Temperature(PERT_LAYER)
    best = HUGE(ONE)
    WRITE(*,'(7x,"[1] TL vs finite-difference  (channel ",i0,", d/dT(",i0,"), TL=",es13.6,")")') &
          RTSolution(ch,1)%Sensor_Channel, PERT_LAYER, tl
    DO kk = 4, 16
      delta = ABS(T0) * 0.1_fp / (2.0_fp**kk)
      Atm(1)%Temperature(PERT_LAYER) = T0 + delta
      Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution_pert, Options=Options )
      Rp = get_out(RTSolution_pert(ch,1),downwelling)
      Atm(1)%Temperature(PERT_LAYER) = T0 - delta
      Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution_pert, Options=Options )
      Rm = get_out(RTSolution_pert(ch,1),downwelling)
      Atm(1)%Temperature(PERT_LAYER) = T0
      fd = ( Rp - Rm ) / ( 2.0_fp*delta )
      ratio = fd / tl
      IF ( ABS(ratio-ONE) < best ) best = ABS(ratio-ONE)
      WRITE(*,'(9x,"delta=",es10.3,"  FD=",es16.9,"  FD/TL=",f14.10)') delta, fd, ratio
    END DO
    ok1 = ( best < TOL_FD )
    WRITE(*,'(7x,"-> best |FD/TL - 1| = ",es11.4,"   ",a)') best, MERGE('PASS','FAIL',ok1)

    ! ----------------------------------------------------------------
    ! Check 2 : Adjoint dot-product test (T perturbation, all layers/profiles)
    ! ----------------------------------------------------------------
    CALL CRTM_Atmosphere_Zero( Atm_TL ) ; CALL CRTM_Surface_Zero( Sfc_TL )
    DO m = 1, N_PROFILES
      DO ii = 1, N_LAYERS
        Atm_TL(m)%Temperature(ii) = 0.5_fp * SIN( 0.7_fp*REAL(ii,fp) + 1.3_fp*REAL(m,fp) )
      END DO
    END DO
    Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, ChannelInfo, &
                                        RTSolution, RTSolution_TL, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'TL fail' ; all_ok=.FALSE. ; RETURN ; END IF

    LHS = ZERO
    CALL CRTM_RTSolution_Zero( RTSolution_AD )
    DO m = 1, N_PROFILES
      DO l = 1, n_Channels
        dy = get_out(RTSolution_TL(l,m),downwelling)
        LHS = LHS + dy*dy
        CALL set_seed( RTSolution_AD(l,m), downwelling, dy )
      END DO
    END DO

    CALL CRTM_Atmosphere_Zero( Atm_AD ) ; CALL CRTM_Surface_Zero( Sfc_AD )
    Error_Status = CRTM_Adjoint( Atm, Sfc, RTSolution_AD, Geometry, ChannelInfo, &
                                 Atm_AD, Sfc_AD, RTSolution, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'AD fail' ; all_ok=.FALSE. ; RETURN ; END IF

    RHS = ZERO
    DO m = 1, N_PROFILES
      DO ii = 1, N_LAYERS
        RHS = RHS + Atm_TL(m)%Temperature(ii) * Atm_AD(m)%Temperature(ii)
      END DO
    END DO
    rel_adj = ABS(LHS-RHS) / MAX(ABS(LHS), TINY(ONE))
    ok2 = ( rel_adj < TOL_ADJ )
    WRITE(*,'(7x,"[2] Adjoint dot-product:  <dy,dy>=",es16.9,"  <x,gx>=",es16.9)') LHS, RHS
    WRITE(*,'(7x,"-> relative difference = ",es11.4,"   ",a)') rel_adj, MERGE('PASS','FAIL',ok2)

    ! ----------------------------------------------------------------
    ! Check 3 : K-Matrix vs Adjoint Jacobian for one channel
    ! ----------------------------------------------------------------
    l0 = ch ; m0 = 1
    CALL CRTM_Atmosphere_Zero( Atm_K ) ; CALL CRTM_Surface_Zero( Sfc_K )
    CALL CRTM_RTSolution_Zero( RTSolution_K )
    DO l = 1, n_Channels
      CALL set_seed( RTSolution_K(l,1), downwelling, ONE )
      CALL set_seed( RTSolution_K(l,2), downwelling, ONE )
    END DO
    Error_Status = CRTM_K_Matrix( Atm, Sfc, RTSolution_K, Geometry, ChannelInfo, &
                                  Atm_K, Sfc_K, RTSolution, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'K fail' ; all_ok=.FALSE. ; RETURN ; END IF

    CALL CRTM_Atmosphere_Zero( Atm_AD ) ; CALL CRTM_Surface_Zero( Sfc_AD )
    CALL CRTM_RTSolution_Zero( RTSolution_AD )
    CALL set_seed( RTSolution_AD(l0,m0), downwelling, ONE )
    Error_Status = CRTM_Adjoint( Atm, Sfc, RTSolution_AD, Geometry, ChannelInfo, &
                                 Atm_AD, Sfc_AD, RTSolution, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'AD fail' ; all_ok=.FALSE. ; RETURN ; END IF

    maxdiff = MAXVAL( ABS( Atm_K(l0,m0)%Temperature - Atm_AD(m0)%Temperature ) )
    scal    = MAX( MAXVAL(ABS(Atm_K(l0,m0)%Temperature)), TINY(ONE) )
    rel_k   = maxdiff / scal
    ok3 = ( rel_k < TOL_K )
    WRITE(*,'(7x,"[3] K vs AD Jacobian (channel ",i0,"):  max|K-AD|/max|K| = ",es11.4,"   ",a)') &
          RTSolution(l0,m0)%Sensor_Channel, rel_k, MERGE('PASS','FAIL',ok3)

    all_ok = ( ok1 .AND. ok2 .AND. ok3 )
  END SUBROUTINE verify

  INCLUDE 'Load_Atm_Data.inc'
  INCLUDE 'Load_Sfc_Data.inc'

END PROGRAM test_Downwelling_TLADK
