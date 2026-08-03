!
! test_VectorRT_IR_Surface
!
! Verification of the polarized INFRARED water surface, which is the first
! polarimetric capability in CRTM outside the microwave.
!
! Background
! ----------
! Before this, an infrared run at n_Stokes > 1 returned Q = U = V = exactly
! zero on every channel at every geometry, because no infrared surface model
! produced anything but the first Stokes component. The polarization was not
! missing from the physics, only from the code: Fresnel reflection from
! seawater is strongly polarized in the infrared, and CRTM_IR_Water_SfcOptics
! already computed rv and rh separately from the Wieliczka and Friedman
! refractive index, then averaged them away on the following line.
!
! The model now takes the V/H SHAPE from Fresnel and the MAGNITUDE from the
! validated WuSmith or Nalli emissivity lookup table, so that
! e_I = (e_V + e_H)/2 is the lookup table value by construction and the
! scalar path cannot move.
!
! What is asserted, and why each is a check the unfixed code fails or that
! guards against a plausible way of getting it wrong
! ----------------------------------------------------------------------
!   1. Q is EXACTLY zero at nadir. At normal incidence rv == rh identically,
!      so any nonzero Q there is an indexing or geometry error, not physics.
!   2. Q grows strictly monotonically with view angle on a window channel.
!      This is the assertion the unfixed code fails: it returned zero at
!      every angle.
!   3. I is bit-identical between n_Stokes = 1 and n_Stokes = 4. This is the
!      construction check. Any arithmetic that perturbs the intensity, for
!      instance re-deriving it as (e_V+e_H)/2 instead of leaving the lookup
!      table value alone, shows up here immediately.
!   4. U and V are exactly zero. Thermal emission from an azimuthally
!      symmetric surface has no third or fourth Stokes component and no
!      infrared azimuth model exists, so a nonzero value would be leakage.
!   5. An OPAQUE channel carries essentially no Q. A channel that cannot see
!      the surface must not acquire a surface polarization signal, which
!      catches a Q injected somewhere other than the surface.
!   6. The polarization bound I^2 >= Q^2 + U^2 + V^2 holds, with the
!      polarized part non-vacuously nonzero.
!   7. Tangent linear against central finite differences for dQ/d(wind speed)
!      and dQ/d(water temperature), the two state variables the infrared
!      surface emissivity depends on.
!   8. The adjoint dot-product identity over the full Stokes vector, which is
!      the only instrument that has ever caught a non-transpose on this path.
!   9. K-matrix equals the adjoint Jacobian.
!
! Clear sky throughout, so no cloud lookup table is involved and a failure
! cannot be blamed on coefficient quality. Night-time geometry, so the solar
! BRDF does not contribute.
!
! Exit: STOP 0 if every check passes, STOP 1 otherwise.
!
! CREATION HISTORY:
!       Written by:     Benjamin Johnson, 03-Aug-2026
!
PROGRAM test_VectorRT_IR_Surface

  USE CRTM_Module
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_VectorRT_IR_Surface'
  CHARACTER(*), PARAMETER :: PATH   = './testinput/'
  CHARACTER(*), PARAMETER :: SENSOR = 'abi_g18'

  INTEGER,  PARAMETER :: N_PROFILES  = 2     ! the ECMWF84 loader fills atm(1) and atm(2)
  INTEGER,  PARAMETER :: N_LAYERS    = 100
  INTEGER,  PARAMETER :: N_ABSORBERS = 6
  INTEGER,  PARAMETER :: N_CLOUDS    = 0
  INTEGER,  PARAMETER :: N_AEROSOLS  = 0

  ! abi_g18 carries ABI bands 7 to 16. Index 8 is the 11.2 micron window,
  ! which sees the surface; index 2 is the 6.2 micron water-vapour band,
  ! which is opaque and must not acquire a surface signal.
  INTEGER,  PARAMETER :: CH_WINDOW = 8
  INTEGER,  PARAMETER :: CH_OPAQUE = 2

  INTEGER,  PARAMETER :: N_ANG = 5
  REAL(fp), PARAMETER :: ZEN(N_ANG) = (/ 0.0_fp, 15.0_fp, 30.0_fp, 45.0_fp, 60.0_fp /)
  REAL(fp), PARAMETER :: ZEN_J      = 55.0_fp   ! Jacobian tests, well away from nadir

  ! Wind speed sits deliberately INSIDE a lookup-table cell, not on a knot.
  ! The WuSmith/Nalli emissivity table is tabulated at 0, 0.5, ... 10, 12.5,
  ! 15 m/s and interpolated piecewise, so a central finite difference taken
  ! at a knot straddles two different polynomial pieces and disagrees with
  ! the analytic tangent linear by about 2e-5 relative. That is an artifact
  ! of differencing across the knot, not a defect: measured at 6.0 m/s the
  ! disagreement is 1.8e-5, and at 6.25 m/s, mid-cell, it is 4.3e-8 with
  ! everything else unchanged. Do not "fix" a future failure by loosening
  ! the tolerance without first checking whether the wind speed has landed
  ! on a knot.
  REAL(fp), PARAMETER :: WIND  = 6.25_fp
  REAL(fp), PARAMETER :: WTEMP = 290.0_fp
  REAL(fp), PARAMETER :: SALIN = 33.0_fp
  REAL(fp), PARAMETER :: NIGHT = 100.0_fp      ! source zenith below the horizon

  ! Finite-difference steps, large enough to clear roundoff on a lookup-table
  ! interpolation and small enough to stay linear.
  REAL(fp), PARAMETER :: DWIND  = 1.0e-3_fp
  REAL(fp), PARAMETER :: DWTEMP = 1.0e-3_fp

  REAL(fp), PARAMETER :: TOL_FD    = 1.0e-6_fp   ! relative, TL against finite difference
  REAL(fp), PARAMETER :: TOL_ADJ   = 1.0e-12_fp  ! relative, dot-product identity
  REAL(fp), PARAMETER :: TOL_K     = 1.0e-12_fp  ! relative, K against AD
  REAL(fp), PARAMETER :: FLOOR_SIG = 1.0e-6_fp   ! non-degeneracy floor for Q

  CHARACTER(256) :: Version
  INTEGER  :: Error_Status, Allocate_Status, n_Channels, l, m, ia
  LOGICAL  :: all_ok
  REAL(fp) :: LHS, RHS, rel, polfrac, worst_bound, maxU, maxV, maxdI
  REAL(fp) :: dQ_ws_tl, dQ_ws_fd, dQ_wt_tl, dQ_wt_fd
  REAL(fp), ALLOCATABLE :: S(:,:,:), Rs(:,:)
  REAL(fp), ALLOCATABLE :: Sp(:,:), Sm(:,:)

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(1)
  TYPE(CRTM_Geometry_type)    :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm_TL(N_PROFILES), Atm_AD(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc_TL(N_PROFILES), Sfc_AD(N_PROFILES)
  TYPE(CRTM_Options_type)     :: Options(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTS(:,:), RTS_TL(:,:), RTS_AD(:,:), RTS_K(:,:)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atm_K(:,:)
  TYPE(CRTM_Surface_type),    ALLOCATABLE :: Sfc_K(:,:)

  all_ok = .TRUE.

  CALL CRTM_Version(Version)
  WRITE(*,'(/5x,a)') 'Polarized infrared water surface (n_Stokes > 1), clear sky'
  WRITE(*,'(5x,a/)') 'CRTM Version: '//TRIM(Version)

  Error_Status = CRTM_Init( (/ SENSOR /), ChannelInfo, File_Path = PATH, Quiet = .TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Init failed', FAILURE ); STOP 1
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
  IF ( n_Channels < CH_WINDOW ) THEN
    CALL Display_Message( PROGRAM_NAME, 'unexpected channel count', FAILURE ); STOP 1
  END IF

  ALLOCATE( RTS(n_Channels,N_PROFILES), RTS_TL(n_Channels,N_PROFILES), &
            RTS_AD(n_Channels,N_PROFILES), RTS_K(n_Channels,N_PROFILES), &
            Atm_K(n_Channels,N_PROFILES), Sfc_K(n_Channels,N_PROFILES), &
            S(4,n_Channels,N_ANG), Rs(n_Channels,N_ANG), &
            Sp(4,n_Channels), Sm(4,n_Channels), STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN; WRITE(*,*) 'Alloc error'; STOP 1; END IF
  CALL CRTM_RTSolution_Create( RTS,    N_LAYERS )
  CALL CRTM_RTSolution_Create( RTS_TL, N_LAYERS )
  CALL CRTM_RTSolution_Create( RTS_AD, N_LAYERS )
  CALL CRTM_RTSolution_Create( RTS_K,  N_LAYERS )
  CALL CRTM_Atmosphere_Create( Atm,    N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atm_TL, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atm_AD, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atm_K,  N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )

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
  END DO

  ! ==================================================================
  ! Forward: view-angle sweep, vector against scalar
  ! ==================================================================
  DO ia = 1, N_ANG
    CALL Set_Scene( ZEN(ia), WIND, WTEMP )
    CALL Run_Forward( 1 )
    DO l = 1, n_Channels ; Rs(l,ia) = RTS(l,1)%Radiance ; END DO
    CALL Run_Forward( 4 )
    DO l = 1, n_Channels ; S(1:4,l,ia) = RTS(l,1)%Stokes(1:4) ; END DO
  END DO

  ! 1. Q exactly zero at nadir
  CALL judge( 'Q is exactly zero at nadir (rv == rh there)', &
              MAXVAL(ABS(S(2,1:n_Channels,1))) == ZERO )

  ! 2. Q strictly increasing with view angle on the window channel.
  !    This is the assertion the unfixed code fails: it returned zero at
  !    every angle, so the strict inequality could never hold.
  CALL judge( 'Q grows strictly with view angle (window channel)', &
              ALL( (/ (S(2,CH_WINDOW,ia+1) > S(2,CH_WINDOW,ia), ia=1,N_ANG-1) /) ) )
  WRITE(*,'(7x,a)') 'window-channel Q by zenith angle:'
  DO ia = 1, N_ANG
    WRITE(*,'(9x,f6.1,a,es24.17,a,es12.4)') ZEN(ia), ' deg  Q = ', S(2,CH_WINDOW,ia), &
          '   Q/I = ', S(2,CH_WINDOW,ia)/S(1,CH_WINDOW,ia)
  END DO

  ! 3. Intensity bit-identical between the scalar and vector runs
  maxdI = MAXVAL(ABS(S(1,1:n_Channels,:) - Rs(1:n_Channels,:)))
  CALL judge( 'I is bit-identical between n_Stokes 1 and 4', maxdI == ZERO )
  WRITE(*,'(7x,a,es24.17)') 'max |I(vector) - I(scalar)| = ', maxdI

  ! 4. U and V exactly zero
  maxU = MAXVAL(ABS(S(3,1:n_Channels,:)))
  maxV = MAXVAL(ABS(S(4,1:n_Channels,:)))
  CALL judge( 'U and V are exactly zero (no IR azimuth model)', &
              maxU == ZERO .AND. maxV == ZERO )

  ! 5. Opaque channel acquires no surface polarization. Compared against the
  !    window channel at the same angle rather than an absolute threshold, so
  !    the check scales with the scene.
  CALL judge( 'opaque channel carries negligible Q versus the window channel', &
              ABS(S(2,CH_OPAQUE,N_ANG)) < 1.0e-6_fp * ABS(S(2,CH_WINDOW,N_ANG)) )
  WRITE(*,'(7x,a,es12.4,a,es12.4)') 'opaque Q = ', S(2,CH_OPAQUE,N_ANG), &
        '   window Q = ', S(2,CH_WINDOW,N_ANG)

  ! 6. Polarization bound, and non-vacuous
  worst_bound = MINVAL( S(1,1:n_Channels,:)**2 - S(2,1:n_Channels,:)**2 &
                      - S(3,1:n_Channels,:)**2 - S(4,1:n_Channels,:)**2 )
  polfrac = MAXVAL(ABS(S(2,1:n_Channels,:))/MAX(ABS(S(1,1:n_Channels,:)),1.0e-300_fp))
  CALL judge( 'polarization bound I^2 >= Q^2+U^2+V^2 holds', worst_bound >= ZERO )
  CALL judge( 'the polarized part is non-vacuously nonzero', polfrac > FLOOR_SIG )
  WRITE(*,'(7x,a,es12.4,a,es12.4)') 'worst bound margin = ', worst_bound, &
        '   max Q/I = ', polfrac

  ! ==================================================================
  ! 7. Tangent linear against central finite differences
  ! ==================================================================
  CALL Set_Scene( ZEN_J, WIND, WTEMP )

  ! ...wind speed
  CALL CRTM_Atmosphere_Zero( Atm_TL ) ; CALL CRTM_Surface_Zero( Sfc_TL )
  DO m = 1, N_PROFILES ; Sfc_TL(m)%Wind_Speed = ONE ; END DO
  CALL Run_TL( 4 )
  dQ_ws_tl = RTS_TL(CH_WINDOW,1)%Stokes(2)
  CALL Set_Scene( ZEN_J, WIND+DWIND, WTEMP ) ; CALL Run_Forward( 4 )
  Sp(1:4,CH_WINDOW) = RTS(CH_WINDOW,1)%Stokes(1:4)
  CALL Set_Scene( ZEN_J, WIND-DWIND, WTEMP ) ; CALL Run_Forward( 4 )
  Sm(1:4,CH_WINDOW) = RTS(CH_WINDOW,1)%Stokes(1:4)
  dQ_ws_fd = (Sp(2,CH_WINDOW) - Sm(2,CH_WINDOW))/(TWO*DWIND)
  rel = ABS(dQ_ws_tl - dQ_ws_fd)/MAX(ABS(dQ_ws_fd),1.0e-300_fp)
  CALL judge( 'TL dQ/d(wind speed) against finite difference', &
              rel < TOL_FD .AND. ABS(dQ_ws_fd) > ZERO )
  WRITE(*,'(7x,a,es24.17,a,es12.4)') 'dQ/dU10  TL = ', dQ_ws_tl, '   rel = ', rel

  ! ...water temperature
  CALL Set_Scene( ZEN_J, WIND, WTEMP )
  CALL CRTM_Atmosphere_Zero( Atm_TL ) ; CALL CRTM_Surface_Zero( Sfc_TL )
  DO m = 1, N_PROFILES ; Sfc_TL(m)%Water_Temperature = ONE ; END DO
  CALL Run_TL( 4 )
  dQ_wt_tl = RTS_TL(CH_WINDOW,1)%Stokes(2)
  CALL Set_Scene( ZEN_J, WIND, WTEMP+DWTEMP ) ; CALL Run_Forward( 4 )
  Sp(1:4,CH_WINDOW) = RTS(CH_WINDOW,1)%Stokes(1:4)
  CALL Set_Scene( ZEN_J, WIND, WTEMP-DWTEMP ) ; CALL Run_Forward( 4 )
  Sm(1:4,CH_WINDOW) = RTS(CH_WINDOW,1)%Stokes(1:4)
  dQ_wt_fd = (Sp(2,CH_WINDOW) - Sm(2,CH_WINDOW))/(TWO*DWTEMP)
  rel = ABS(dQ_wt_tl - dQ_wt_fd)/MAX(ABS(dQ_wt_fd),1.0e-300_fp)
  CALL judge( 'TL dQ/d(water temperature) against finite difference', &
              rel < TOL_FD .AND. ABS(dQ_wt_fd) > ZERO )
  WRITE(*,'(7x,a,es24.17,a,es12.4)') 'dQ/dTs   TL = ', dQ_wt_tl, '   rel = ', rel

  ! ==================================================================
  ! 8. Adjoint dot-product identity over the full Stokes vector
  !
  !    <TL(x), y> == <x, AD(y)>, with x spanning both surface variables and
  !    y seeded on every Stokes component so the polarized rows of the
  !    transpose are actually exercised.
  ! ==================================================================
  CALL Set_Scene( ZEN_J, WIND, WTEMP )
  CALL CRTM_Atmosphere_Zero( Atm_TL ) ; CALL CRTM_Surface_Zero( Sfc_TL )
  DO m = 1, N_PROFILES
    Sfc_TL(m)%Wind_Speed        = ONE
    Sfc_TL(m)%Water_Temperature = 0.5_fp
  END DO
  CALL Run_TL( 4 )

  LHS = ZERO
  CALL CRTM_RTSolution_Zero( RTS_AD )
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      RTS_AD(l,m)%Stokes(1:4) = RTS_TL(l,m)%Stokes(1:4)
      LHS = LHS + SUM( RTS_TL(l,m)%Stokes(1:4)**2 )
    END DO
  END DO
  CALL CRTM_Atmosphere_Zero( Atm_AD ) ; CALL CRTM_Surface_Zero( Sfc_AD )
  DO m = 1, N_PROFILES
    Options(m)%n_Stokes        = 4
    Options(m)%RT_Algorithm_Id = RT_ADA
  END DO
  Error_Status = CRTM_Adjoint( Atm, Sfc, RTS_AD, Geometry, ChannelInfo, &
                               Atm_AD, Sfc_AD, RTS, Options=Options )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL judge( 'CRTM_Adjoint call', .FALSE. )
  ELSE
    RHS = ZERO
    DO m = 1, N_PROFILES
      RHS = RHS + ONE     * Sfc_AD(m)%Wind_Speed        &
                + 0.5_fp  * Sfc_AD(m)%Water_Temperature
    END DO
    rel = ABS(LHS - RHS)/MAX(ABS(LHS),1.0e-300_fp)
    CALL judge( 'adjoint dot-product identity over all four Stokes components', &
                rel < TOL_ADJ .AND. LHS > ZERO )
    WRITE(*,'(7x,a,es24.17)') '<TL,y>  = ', LHS
    WRITE(*,'(7x,a,es24.17,a,es12.4)') '<x,ADy> = ', RHS, '   rel = ', rel
  END IF

  ! ==================================================================
  ! 9. K-matrix against the adjoint
  ! ==================================================================
  CALL CRTM_RTSolution_Zero( RTS_K )
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      RTS_K(l,m)%Stokes(1:4) = RTS_TL(l,m)%Stokes(1:4)
    END DO
    CALL CRTM_Surface_Zero( Sfc_K(:,m) )
  END DO
  Error_Status = CRTM_K_Matrix( Atm, Sfc, RTS_K, Geometry, ChannelInfo, &
                                Atm_K, Sfc_K, RTS, Options=Options )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL judge( 'CRTM_K_Matrix call', .FALSE. )
  ELSE
    LHS = ZERO ; RHS = ZERO
    DO m = 1, N_PROFILES
      DO l = 1, n_Channels
        LHS = LHS + Sfc_K(l,m)%Wind_Speed + Sfc_K(l,m)%Water_Temperature
      END DO
      RHS = RHS + Sfc_AD(m)%Wind_Speed + Sfc_AD(m)%Water_Temperature
    END DO
    rel = ABS(LHS - RHS)/MAX(ABS(RHS),1.0e-300_fp)
    CALL judge( 'K-matrix summed over channels equals the adjoint', &
                rel < TOL_K .AND. ABS(RHS) > ZERO )
    WRITE(*,'(7x,a,es24.17,a,es12.4)') 'sum(K) = ', LHS, '   rel = ', rel
  END IF

  ! ------------------------------------------------------------------
  CALL CRTM_Atmosphere_Destroy(Atm)
  CALL CRTM_RTSolution_Destroy(RTS)
  Error_Status = CRTM_Destroy( ChannelInfo )

  IF ( all_ok ) THEN
    WRITE(*,'(/5x,a/)') 'ALL CHECKS PASSED'
    STOP 0
  ELSE
    WRITE(*,'(/5x,a/)') 'ONE OR MORE CHECKS FAILED'
    STOP 1
  END IF

CONTAINS

  SUBROUTINE Set_Scene( zenith, wind_speed, water_temp )
    REAL(fp), INTENT(IN) :: zenith, wind_speed, water_temp
    INTEGER :: mm
    DO mm = 1, N_PROFILES
      Sfc(mm) = CRTM_Surface_type()
      Sfc(mm)%Water_Coverage    = ONE
      Sfc(mm)%Water_Type        = 1
      Sfc(mm)%Water_Temperature = water_temp
      Sfc(mm)%Wind_Speed        = wind_speed
      Sfc(mm)%Salinity          = SALIN
      CALL CRTM_Geometry_SetValue( Geometry(mm), Sensor_Zenith_Angle = zenith, &
                                   Source_Zenith_Angle = NIGHT )
    END DO
  END SUBROUTINE Set_Scene

  SUBROUTINE Run_Forward( ns )
    INTEGER, INTENT(IN) :: ns
    INTEGER :: mm
    DO mm = 1, N_PROFILES
      Options(mm)%n_Stokes        = ns
      Options(mm)%RT_Algorithm_Id = RT_ADA
    END DO
    Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTS, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN
      CALL Display_Message( PROGRAM_NAME, 'CRTM_Forward failed', FAILURE ); STOP 1
    END IF
  END SUBROUTINE Run_Forward

  SUBROUTINE Run_TL( ns )
    INTEGER, INTENT(IN) :: ns
    INTEGER :: mm
    DO mm = 1, N_PROFILES
      Options(mm)%n_Stokes        = ns
      Options(mm)%RT_Algorithm_Id = RT_ADA
    END DO
    Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, &
                                        ChannelInfo, RTS, RTS_TL, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN
      CALL Display_Message( PROGRAM_NAME, 'CRTM_Tangent_Linear failed', FAILURE ); STOP 1
    END IF
  END SUBROUTINE Run_TL

  SUBROUTINE judge( label, ok )
    CHARACTER(*), INTENT(IN) :: label
    LOGICAL,      INTENT(IN) :: ok
    IF ( ok ) THEN
      WRITE(*,'(5x,"PASS  ",a)') label
    ELSE
      WRITE(*,'(5x,"FAIL  ",a)') label
      all_ok = .FALSE.
    END IF
  END SUBROUTINE judge

  INCLUDE 'Load_ECMWF84_Atm_Data.inc'

END PROGRAM test_VectorRT_IR_Surface
