!
! test_VectorRT_PARMIO_TLAD
!
! Jacobian coverage for the PARMIO polarimetric surface, through the full
! radiative transfer chain.
!
! Why a separate test
! -------------------
! PARMIO is the microwave-water backend at and above PARMIO_FREQ_THRESHOLD
! (200 GHz) whenever its lookup table is loaded, which by default it is. It
! carries a four-Stokes azimuth model of its own
! (PARMIO_Azimuth_Module.f90), entirely independent of the FASTEM one, with its
! own tangent-linear and adjoint. Every other polarimetric test loads FASTEM4
! and so never touches it.
!
! It was briefly believed that PARMIO could not be reached through the
! atmosphere at all, because only two shipped sensors have channels above the
! gate and the mwr_aws one is at 325.15 GHz, on a water-vapour line, where the
! surface is invisible: the measured Stokes U there is around 1e-10 of the
! intensity. That reasoning was wrong for the other one. TROPICS channel 12 sits
! at 204.783 GHz, between the 183 and 325 GHz lines, and is transparent enough
! that the polarimetric surface reaches the top of the atmosphere with
! U/I = 1.6e-3, comparable to the best FASTEM window channels. So PARMIO is
! fully testable end to end, and this is that test.
!
! What it checks
! --------------
! Clear sky over ocean at n_Stokes = 4, so no cloud lookup table is involved and
! a failure cannot be blamed on coefficient quality:
!
!   1. the PARMIO channel really does carry a polarimetric signal to the top of
!      the atmosphere, otherwise everything below is vacuous;
!   2. dU/d(wind direction) against central finite differences. This is the
!      observable polarimetric microwave exists for, on the PARMIO backend;
!   3. the adjoint dot-product identity with the surface control variables
!      perturbed, which is the only check that tests the transpose;
!   4. K against AD on the surface Jacobians.
!
! No coefficient gating beyond the sensor itself: clear sky needs no cloud LUT.
!

PROGRAM test_VectorRT_PARMIO_TLAD

  USE CRTM_Module
  USE CRTM_SpcCoeff, ONLY: SC
  USE CRTM_MW_Water_SfcOptics, ONLY: PARMIO_Is_Active_At
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_VectorRT_PARMIO_TLAD'
  CHARACTER(*), PARAMETER :: PATH   = './testinput/'
  CHARACTER(*), PARAMETER :: SENSOR = 'tms_tropics-01'
  ! Drive PARMIO down into a band the LUT actually covers.
  !
  ! The default 200 GHz floor lands in a hole: the production table's
  ! sss_nominal_m group stops at 183.31 GHz and sss_nominal_h starts at 229,
  ! so this sensor's only above-floor channel, 204.783 GHz, has no data and
  ! the dispatcher correctly declines it. Before coverage was checked it was
  ! served anyway, silently evaluated at 229 GHz.
  !
  ! Lowering the floor puts the 91.319 GHz channel on PARMIO, which sits well
  ! inside sss_nominal_m and is a window channel, so the surface polarimetric
  ! signal actually reaches the top of the atmosphere. That makes this a
  ! stronger test than it was: the signal is real rather than clamped.
  LOGICAL,      PARAMETER :: USE_PARMIO_EVERYWHERE = .TRUE.

  INTEGER,  PARAMETER :: N_PROFILES  = 2      ! the ECMWF84 loader fills atm(1) and atm(2)
  INTEGER,  PARAMETER :: N_LAYERS    = 100
  INTEGER,  PARAMETER :: N_ABSORBERS = 6
  REAL(fp), PARAMETER :: ZENITH      = 53.0_fp
  REAL(fp), PARAMETER :: SENSOR_AZI  = 40.0_fp
  REAL(fp), PARAMETER :: WIND_DIR    = 100.0_fp   ! relative azimuth 60 deg
  REAL(fp), PARAMETER :: WIND_SPEED  = 12.0_fp

  REAL(fp), PARAMETER :: TOL_FD   = 1.0e-3_fp
  REAL(fp), PARAMETER :: TOL_ADJ  = 1.0e-12_fp
  REAL(fp), PARAMETER :: TOL_K    = 1.0e-9_fp
  REAL(fp), PARAMETER :: UI_FLOOR = 1.0e-6_fp    ! |U/I| a PARMIO channel must clear

  CHARACTER(256) :: Version
  INTEGER  :: Error_Status, Allocate_Status, n_Channels, l, m, ch_parmio, kk
  LOGICAL  :: ok_signal, ok_fd, ok_adj, ok_k, all_ok
  REAL(fp) :: ui, best_ui, tl, fd, Rp, Rm, delta, X0, best
  REAL(fp) :: LHS, RHS, RHS_sfc, dy, rel_adj, maxdiff, scal, rel_k

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(1)
  TYPE(CRTM_Geometry_type)    :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm(N_PROFILES), Atm_TL(N_PROFILES), Atm_AD(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc(N_PROFILES), Sfc_TL(N_PROFILES), Sfc_AD(N_PROFILES)
  TYPE(CRTM_Options_type)     :: Options(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTS(:,:), RTS_pert(:,:), RTS_TL(:,:), RTS_AD(:,:), RTS_K(:,:)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atm_K(:,:)
  TYPE(CRTM_Surface_type),    ALLOCATABLE :: Sfc_K(:,:)

  CALL CRTM_Version(Version)
  WRITE(*,'(/5x,a)') 'PARMIO polarimetric surface: Jacobians through the full RT chain'
  WRITE(*,'(5x,a/)') 'CRTM Version: '//TRIM(Version)

  Error_Status = CRTM_Init( (/ SENSOR /), ChannelInfo, File_Path=PATH, Quiet=.TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Init failed', FAILURE ); STOP 1
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))

  ALLOCATE( RTS(n_Channels,N_PROFILES), RTS_pert(n_Channels,N_PROFILES), &
            RTS_TL(n_Channels,N_PROFILES), RTS_AD(n_Channels,N_PROFILES), &
            RTS_K(n_Channels,N_PROFILES), Atm_K(n_Channels,N_PROFILES), &
            Sfc_K(n_Channels,N_PROFILES), STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN; WRITE(*,*) 'Alloc error'; STOP 1; END IF
  CALL CRTM_RTSolution_Create( RTS,      N_LAYERS )
  CALL CRTM_RTSolution_Create( RTS_pert, N_LAYERS )
  CALL CRTM_RTSolution_Create( RTS_TL,   N_LAYERS )
  CALL CRTM_RTSolution_Create( RTS_AD,   N_LAYERS )
  CALL CRTM_RTSolution_Create( RTS_K,    N_LAYERS )
  CALL CRTM_Atmosphere_Create( Atm,    N_LAYERS, N_ABSORBERS, 0, 0 )
  CALL CRTM_Atmosphere_Create( Atm_TL, N_LAYERS, N_ABSORBERS, 0, 0 )
  CALL CRTM_Atmosphere_Create( Atm_AD, N_LAYERS, N_ABSORBERS, 0, 0 )
  CALL CRTM_Atmosphere_Create( Atm_K,  N_LAYERS, N_ABSORBERS, 0, 0 )

  CALL Load_ECMWF84_Atm_Data()

  DO m = 1, N_PROFILES
    Atm_TL(m)%Climatology = Atm(m)%Climatology
    Atm_TL(m)%Absorber_Id = Atm(m)%Absorber_Id ; Atm_TL(m)%Absorber_Units = Atm(m)%Absorber_Units
    Atm_AD(m)%Climatology = Atm(m)%Climatology
    Atm_AD(m)%Absorber_Id = Atm(m)%Absorber_Id ; Atm_AD(m)%Absorber_Units = Atm(m)%Absorber_Units
    DO l = 1, n_Channels
      Atm_K(l,m)%Climatology = Atm(m)%Climatology
      Atm_K(l,m)%Absorber_Id = Atm(m)%Absorber_Id ; Atm_K(l,m)%Absorber_Units = Atm(m)%Absorber_Units
    END DO
    Sfc(m)%Water_Coverage    = ONE
    Sfc(m)%Water_Type        = 1
    Sfc(m)%Water_Temperature = 290.0_fp
    Sfc(m)%Wind_Speed        = WIND_SPEED
    Sfc(m)%Wind_Direction    = WIND_DIR
    Sfc(m)%Salinity          = 33.0_fp
    CALL CRTM_Geometry_SetValue( Geometry(m), Sensor_Zenith_Angle  = ZENITH, &
                                              Sensor_Azimuth_Angle = SENSOR_AZI )
    Options(m)%n_Stokes          = 4
    Options(m)%Use_PARMIO_MWSSEM = USE_PARMIO_EVERYWHERE
    Options(m)%RT_Algorithm_Id = RT_ADA
  END DO

  ! ------------------------------------------------------------------
  ! 1. Find a PARMIO channel that actually carries U to the top
  ! ------------------------------------------------------------------
  Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTS, Options=Options )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Forward failed', FAILURE ); STOP 1
  END IF
  ch_parmio = 0 ; best_ui = ZERO
  WRITE(*,'(5x,a)') 'ch   freq(GHz)  backend      Stokes I         Stokes U        |U/I|'
  DO l = 1, n_Channels
    ui = ABS(RTS(l,1)%Stokes(3)) / MAX(ABS(RTS(l,1)%Stokes(1)), TINY(ONE))
    WRITE(*,'(4x,i2,2x,f10.3,2x,a,2(2x,es15.6),2x,es11.3)') &
      l, SC(1)%Frequency(l), MERGE('PARMIO','FASTEM', PARMIO_Is_Active_At(SC(1)%Frequency(l), USE_PARMIO_EVERYWHERE)), &
      RTS(l,1)%Stokes(1), RTS(l,1)%Stokes(3), ui
    IF ( PARMIO_Is_Active_At(SC(1)%Frequency(l), USE_PARMIO_EVERYWHERE) .AND. ui > best_ui ) THEN
      best_ui = ui ; ch_parmio = l
    END IF
  END DO
  ok_signal = ( ch_parmio > 0 .AND. best_ui > UI_FLOOR )
  IF ( .NOT. ok_signal ) THEN
    WRITE(*,'(/5x,a)') 'FAIL: no PARMIO channel carries a polarimetric signal to the top'
    WRITE(*,'(5x,a)')  '      of the atmosphere, so the Jacobian checks would be vacuous.'
    STOP 1
  END IF
  WRITE(*,'(/5x,a,i0,a,f8.3,a,es11.3)') 'PARMIO channel under test: ', ch_parmio, &
        '  at ', SC(1)%Frequency(ch_parmio), ' GHz,  |U/I| = ', best_ui

  ! ------------------------------------------------------------------
  ! 2. dU / d(wind direction), tangent linear against finite differences
  ! ------------------------------------------------------------------
  CALL CRTM_Atmosphere_Zero( Atm_TL ) ; CALL CRTM_Surface_Zero( Sfc_TL )
  Sfc_TL(1)%Wind_Direction = ONE
  Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, ChannelInfo, &
                                      RTS, RTS_TL, Options=Options )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'TL failed', FAILURE ); STOP 1
  END IF
  tl = RTS_TL(ch_parmio,1)%Stokes(3)
  X0 = Sfc(1)%Wind_Direction
  best = HUGE(ONE)
  WRITE(*,'(/5x,a,es14.6)') 'TL  dU/d(wind direction) = ', tl
  DO kk = 4, 14
    delta = ABS(X0) * 0.1_fp / (2.0_fp**kk)
    Sfc(1)%Wind_Direction = X0 + delta
    Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTS_pert, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN; WRITE(*,*) 'FD fail'; STOP 1; END IF
    Rp = RTS_pert(ch_parmio,1)%Stokes(3)
    Sfc(1)%Wind_Direction = X0 - delta
    Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTS_pert, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN; WRITE(*,*) 'FD fail'; STOP 1; END IF
    Rm = RTS_pert(ch_parmio,1)%Stokes(3)
    Sfc(1)%Wind_Direction = X0
    fd = ( Rp - Rm ) / ( 2.0_fp*delta )
    IF ( ABS(fd/tl - ONE) < best ) best = ABS(fd/tl - ONE)
  END DO
  ok_fd = ( best < TOL_FD )
  WRITE(*,'(5x,a,es11.4,a,l1)') 'best |FD/TL - 1|         = ', best, '   pass = ', ok_fd

  ! ------------------------------------------------------------------
  ! 3. Adjoint dot product, surface control variables perturbed
  ! ------------------------------------------------------------------
  CALL CRTM_Atmosphere_Zero( Atm_TL ) ; CALL CRTM_Surface_Zero( Sfc_TL )
  DO m = 1, N_PROFILES
    DO l = 1, N_LAYERS
      Atm_TL(m)%Temperature(l) = 0.5_fp * SIN( 0.7_fp*REAL(l,fp) + 1.3_fp*REAL(m,fp) )
    END DO
    Sfc_TL(m)%Wind_Speed        = 0.30_fp + 0.10_fp*REAL(m,fp)
    Sfc_TL(m)%Wind_Direction    = 2.00_fp - 0.50_fp*REAL(m,fp)
    Sfc_TL(m)%Water_Temperature = 0.20_fp + 0.05_fp*REAL(m,fp)
    Sfc_TL(m)%Salinity          = 0.10_fp
  END DO
  Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, ChannelInfo, &
                                      RTS, RTS_TL, Options=Options )
  IF ( Error_Status /= SUCCESS ) THEN; WRITE(*,*) 'TL fail'; STOP 1; END IF

  LHS = ZERO
  CALL CRTM_RTSolution_Zero( RTS_AD )
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      DO kk = 1, 4
        dy = RTS_TL(l,m)%Stokes(kk)
        LHS = LHS + dy*dy
        RTS_AD(l,m)%Stokes(kk) = dy
      END DO
    END DO
  END DO
  CALL CRTM_Atmosphere_Zero( Atm_AD ) ; CALL CRTM_Surface_Zero( Sfc_AD )
  Error_Status = CRTM_Adjoint( Atm, Sfc, RTS_AD, Geometry, ChannelInfo, &
                               Atm_AD, Sfc_AD, RTS, Options=Options )
  IF ( Error_Status /= SUCCESS ) THEN; WRITE(*,*) 'AD fail'; STOP 1; END IF

  RHS = ZERO ; RHS_sfc = ZERO
  DO m = 1, N_PROFILES
    DO l = 1, N_LAYERS
      RHS = RHS + Atm_TL(m)%Temperature(l) * Atm_AD(m)%Temperature(l)
    END DO
    RHS_sfc = RHS_sfc &
            + Sfc_TL(m)%Wind_Speed        * Sfc_AD(m)%Wind_Speed        &
            + Sfc_TL(m)%Wind_Direction    * Sfc_AD(m)%Wind_Direction    &
            + Sfc_TL(m)%Water_Temperature * Sfc_AD(m)%Water_Temperature &
            + Sfc_TL(m)%Salinity          * Sfc_AD(m)%Salinity
  END DO
  RHS = RHS + RHS_sfc
  rel_adj = ABS(LHS-RHS) / MAX(ABS(LHS), TINY(ONE))
  ok_adj = ( rel_adj < TOL_ADJ )
  WRITE(*,'(/5x,a,es16.9,a,es16.9)') '<dy,dy> = ', LHS, '   <x,gx> = ', RHS
  WRITE(*,'(5x,a,f8.4,a)')      'surface share of <x,gx>  = ', 100.0_fp*RHS_sfc/MAX(ABS(RHS),TINY(ONE)), ' %'
  WRITE(*,'(5x,a,es11.4,a,l1)') 'adjoint dot product      = ', rel_adj, '   pass = ', ok_adj

  ! ------------------------------------------------------------------
  ! 4. K against AD on the surface Jacobians
  ! ------------------------------------------------------------------
  CALL CRTM_Atmosphere_Zero( Atm_K ) ; CALL CRTM_Surface_Zero( Sfc_K )
  CALL CRTM_RTSolution_Zero( RTS_K )
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      RTS_K(l,m)%Stokes(3) = ONE          ! seed U, the polarimetric component
    END DO
  END DO
  Error_Status = CRTM_K_Matrix( Atm, Sfc, RTS_K, Geometry, ChannelInfo, &
                                Atm_K, Sfc_K, RTS, Options=Options )
  IF ( Error_Status /= SUCCESS ) THEN; WRITE(*,*) 'K fail'; STOP 1; END IF

  CALL CRTM_Atmosphere_Zero( Atm_AD ) ; CALL CRTM_Surface_Zero( Sfc_AD )
  CALL CRTM_RTSolution_Zero( RTS_AD )
  RTS_AD(ch_parmio,1)%Stokes(3) = ONE
  Error_Status = CRTM_Adjoint( Atm, Sfc, RTS_AD, Geometry, ChannelInfo, &
                               Atm_AD, Sfc_AD, RTS, Options=Options )
  IF ( Error_Status /= SUCCESS ) THEN; WRITE(*,*) 'AD fail'; STOP 1; END IF

  maxdiff = MAX( ABS( Sfc_K(ch_parmio,1)%Wind_Speed     - Sfc_AD(1)%Wind_Speed ),     &
                 ABS( Sfc_K(ch_parmio,1)%Wind_Direction - Sfc_AD(1)%Wind_Direction ), &
                 ABS( Sfc_K(ch_parmio,1)%Water_Temperature - Sfc_AD(1)%Water_Temperature ) )
  scal    = MAX( ABS(Sfc_K(ch_parmio,1)%Wind_Speed), ABS(Sfc_K(ch_parmio,1)%Wind_Direction), &
                 ABS(Sfc_K(ch_parmio,1)%Water_Temperature), TINY(ONE) )
  rel_k   = maxdiff / scal
  ok_k    = ( rel_k < TOL_K )
  WRITE(*,'(5x,a,es11.4,a,l1)') 'K vs AD, surface         = ', rel_k, '   pass = ', ok_k
  WRITE(*,'(5x,a,es14.6)')      'dU/d(wind direction), K  = ', Sfc_K(ch_parmio,1)%Wind_Direction

  all_ok = ok_signal .AND. ok_fd .AND. ok_adj .AND. ok_k
  Error_Status = CRTM_Destroy( ChannelInfo )

  IF ( all_ok ) THEN
    WRITE(*,'(/5x,a/)') 'PASS: PARMIO polarimetric surface Jacobians verified end to end'
    STOP 0
  ELSE
    WRITE(*,'(/5x,a/)') 'FAIL: PARMIO polarimetric surface Jacobians'
    STOP 1
  END IF

CONTAINS

  INCLUDE 'Load_ECMWF84_Atm_Data.inc'

END PROGRAM test_VectorRT_PARMIO_TLAD
