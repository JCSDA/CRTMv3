!
! test_MW_Sounder_Physics
!
! Baseline-independent physics verification for the microwave sounder family,
! anchored to shipped AMSU-A heritage inside the same run. Five sensors in
! one CRTM_Init and one multi-sensor forward:
!
!   amsua_n19       - the heritage anchor: its 57.290344 GHz line-splitting
!                     channels (10-14) must peak at the operationally
!                     validated pressures (33.9 / 17.5 / 7.6 / 3.7 / 1.9 hPa).
!                     If the anchor itself moves, the suite must say so.
!   mwts3_fy3e      - carries the SAME line-splitting design (ch 13-17);
!                     the coefficients (crtm-coeffgen, NWP-SAF passbands,
!                     fixed-LBLRTM path per crtm-coeffgen#75) must reproduce
!                     the anchor ladder. This is the check that caught both
!                     the double-offset conversion defect (four channels
!                     collapsed to near-copies, all peaking ~49 hPa) and the
!                     multi-band convolution defect (#71).
!   mwhs2_fy3e      - the 118.75 GHz bank: the line-center channel must peak
!                     high (mesosphere-adjacent), the 166 GHz window low.
!   mwrirm_fy3g     - conical imager: brightness-temperature sanity.
!   gems2_amethyst  - 118.75 GHz line-splitting smallsat sounder: monotone
!                     weighting-function descent up the line.
!
! All five: BT within physical bounds, adjoint dot-product closure over
! T + H2O, K equal to AD on a probe channel.
!
! The non-heritage coefficient pairs are pre-release: the test is registered
! only when all four are present (symlinked) in the source testinput
! directory, following the OMPS/TEMPO gating pattern.
!
! Exit: STOP 0 if every check passes, STOP 1 otherwise.
!
! CREATION HISTORY:
!       Written by:     Benjamin Johnson, 28-Jul-2026
!                       Companion to test_OMPS_UV_Physics and
!                       test_TEMPO_UVVIS_Physics.
!
PROGRAM test_MW_Sounder_Physics

  USE CRTM_Module
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_MW_Sounder_Physics'
  CHARACTER(*), PARAMETER :: PATH = './testinput/'
  INTEGER, PARAMETER :: N_SENSORS = 5
  CHARACTER(14), PARAMETER :: SENSORS(N_SENSORS) = &
    (/ 'amsua_n19     ', 'mwts3_fy3e    ', 'mwhs2_fy3e    ', &
       'mwrirm_fy3g   ', 'gems2_amethyst' /)
  INTEGER, PARAMETER :: S_AMSUA = 1, S_MWTS3 = 2, S_MWHS2 = 3, &
                        S_MWRIRM = 4, S_GEMS2 = 5

  INTEGER,  PARAMETER :: N_PROFILES  = 2
  INTEGER,  PARAMETER :: N_LAYERS    = 100
  INTEGER,  PARAMETER :: N_ABSORBERS = 6
  REAL(fp), PARAMETER :: ZENITH = 45.0_fp
  INTEGER,  PARAMETER :: IDX_H2O = 1

  ! AMSU-A heritage line-splitting ladder (hPa), channels 10-14, measured on
  ! the shipped amsua_n19 with this driver convention. The same instrument
  ! design flies as MWTS-2 ch 9-13 and MWTS-3 ch 13-17.
  REAL(fp), PARAMETER :: LADDER(5) = &
    (/ 33.93_fp, 17.49_fp, 7.55_fp, 3.70_fp, 1.91_fp /)
  REAL(fp), PARAMETER :: LADDER_RTOL = 0.35_fp   ! covers backend evolution

  INTEGER  :: Error_Status, Allocate_Status
  INTEGER  :: i, l, m, ii, ntot, kpeak, l0
  INTEGER  :: n_per(N_SENSORS), off(N_SENSORS)
  REAL(fp) :: LHS, RHS, dy, rel, pk, pk_prev
  REAL(fp), ALLOCATABLE :: BT(:,:)
  REAL(fp) :: wf(N_LAYERS), peaks(5)
  LOGICAL  :: all_ok, ok

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)    :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm_TL(N_PROFILES), Atm_AD(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc_TL(N_PROFILES), Sfc_AD(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_TL(:,:), RTSolution_AD(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_K(:,:)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atm_K(:,:)
  TYPE(CRTM_Surface_type),    ALLOCATABLE :: Sfc_K(:,:)

  all_ok = .TRUE.
  WRITE(*,'(/5x,a)') 'MW sounder physics verification (AMSU-A anchored)'

  Error_Status = CRTM_Init( SENSORS, ChannelInfo, File_Path=PATH, Quiet=.TRUE. )
  CALL judge( 'CRTM_Init (five sensors, one call)', Error_Status == SUCCESS )
  IF ( Error_Status /= SUCCESS ) STOP 1
  ntot = 0
  DO i = 1, N_SENSORS
    n_per(i) = CRTM_ChannelInfo_n_Channels(ChannelInfo(i))
    off(i)   = ntot
    ntot     = ntot + n_per(i)
  END DO
  CALL judge( 'channel counts 15/17/15/26/24', &
              ALL( n_per == (/15, 17, 15, 26, 24/) ) )

  ALLOCATE( RTSolution(ntot,N_PROFILES), RTSolution_TL(ntot,N_PROFILES), &
            RTSolution_AD(ntot,N_PROFILES), RTSolution_K(ntot,N_PROFILES), &
            Atm_K(ntot,N_PROFILES), Sfc_K(ntot,N_PROFILES), BT(ntot,N_PROFILES), &
            STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN; WRITE(*,*) 'Alloc error'; STOP 1; END IF
  CALL CRTM_RTSolution_Create( RTSolution, N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_TL, N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_AD, N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_K, N_LAYERS )
  CALL CRTM_Atmosphere_Create( Atm,    N_LAYERS, N_ABSORBERS, 0, 0 )
  CALL CRTM_Atmosphere_Create( Atm_TL, N_LAYERS, N_ABSORBERS, 0, 0 )
  CALL CRTM_Atmosphere_Create( Atm_AD, N_LAYERS, N_ABSORBERS, 0, 0 )
  CALL CRTM_Atmosphere_Create( Atm_K,  N_LAYERS, N_ABSORBERS, 0, 0 )

  CALL Load_ECMWF84_Atm_Data()
  DO m = 1, N_PROFILES
    Atm_TL(m)%Climatology = Atm(m)%Climatology
    Atm_TL(m)%Absorber_ID = Atm(m)%Absorber_ID ; Atm_TL(m)%Absorber_Units = Atm(m)%Absorber_Units
    Atm_AD(m)%Climatology = Atm(m)%Climatology
    Atm_AD(m)%Absorber_ID = Atm(m)%Absorber_ID ; Atm_AD(m)%Absorber_Units = Atm(m)%Absorber_Units
    DO l = 1, ntot
      Atm_K(l,m)%Climatology = Atm(m)%Climatology
      Atm_K(l,m)%Absorber_ID = Atm(m)%Absorber_ID ; Atm_K(l,m)%Absorber_Units = Atm(m)%Absorber_Units
    END DO
    Sfc(m)%Water_Coverage = ONE ; Sfc(m)%Water_Type = 1
    Sfc(m)%Water_Temperature = 290.0_fp ; Sfc(m)%Wind_Speed = 6.0_fp ; Sfc(m)%Salinity = 33.0_fp
    CALL CRTM_Geometry_SetValue( Geometry(m), Sensor_Zenith_Angle = ZENITH )
  END DO

  Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution )
  CALL judge( 'multi-sensor forward', Error_Status == SUCCESS )
  IF ( Error_Status /= SUCCESS ) STOP 1
  DO m = 1, N_PROFILES
    DO l = 1, ntot
      BT(l,m) = RTSolution(l,m)%Brightness_Temperature
    END DO
  END DO
  CALL judge( 'BT within 100-350 K everywhere', &
              ALL( BT > 100.0_fp ) .AND. ALL( BT < 350.0_fp ) )

  CALL CRTM_Atmosphere_Zero( Atm_K ) ; CALL CRTM_Surface_Zero( Sfc_K )
  CALL CRTM_RTSolution_Zero( RTSolution_K )
  DO m = 1, N_PROFILES
    DO l = 1, ntot
      RTSolution_K(l,m)%Brightness_Temperature = ONE
    END DO
  END DO
  Error_Status = CRTM_K_Matrix( Atm, Sfc, RTSolution_K, Geometry, ChannelInfo, &
                                Atm_K, Sfc_K, RTSolution )
  CALL judge( 'multi-sensor K-Matrix', Error_Status == SUCCESS )
  IF ( Error_Status /= SUCCESS ) STOP 1

  ! ---- heritage anchor: AMSU-A ch 10-14 ladder ----
  DO ii = 1, 5
    peaks(ii) = t_peak_hpa( S_AMSUA, 9 + ii )
  END DO
  CALL check_ladder( 'amsua_n19 ch10-14 heritage ladder', peaks )

  ! ---- MWTS-3: same design, generated coefficients must match the anchor ----
  DO ii = 1, 5
    peaks(ii) = t_peak_hpa( S_MWTS3, 12 + ii )
  END DO
  CALL check_ladder( 'mwts3_fy3e ch13-17 matches the AMSU-A ladder', peaks )

  ! ---- MWHS-2: 118 GHz line center high, 166 GHz window low ----
  pk = t_peak_hpa( S_MWHS2, 2 )
  CALL judge( 'mwhs2_fy3e 118.75 GHz line-center channel peaks above 15 hPa', &
              pk < 15.0_fp )
  WRITE(*,'(9x,"line-center peak ",f8.2," hPa")') pk
  pk = t_peak_hpa( S_MWHS2, 10 )
  CALL judge( 'mwhs2_fy3e 166 GHz window channel peaks below 700 hPa', &
              pk > 700.0_fp )

  ! ---- GEMS2: monotone descent up the 118.75 GHz line ----
  pk_prev = t_peak_hpa( S_GEMS2, 5 )
  pk = t_peak_hpa( S_GEMS2, 10 )
  ok = pk < pk_prev
  pk_prev = pk
  pk = t_peak_hpa( S_GEMS2, 15 )
  ok = ok .AND. ( pk < pk_prev ) .AND. ( pk < 10.0_fp )
  CALL judge( 'gems2_amethyst weighting functions climb the 118 GHz line '// &
              '(ch5 > ch10 > ch15, line center above 10 hPa)', ok )

  ! ---- adjoint closure over everything ----
  CALL CRTM_Atmosphere_Zero( Atm_TL ) ; CALL CRTM_Surface_Zero( Sfc_TL )
  DO m = 1, N_PROFILES
    DO ii = 1, N_LAYERS
      Atm_TL(m)%Temperature(ii)      = 0.5_fp * SIN( 0.7_fp*REAL(ii,fp) + 1.3_fp*REAL(m,fp) )
      Atm_TL(m)%Absorber(ii,IDX_H2O) = 0.05_fp * Atm(m)%Absorber(ii,IDX_H2O) &
                                        * COS( 0.9_fp*REAL(ii,fp) + 0.4_fp*REAL(m,fp) )
    END DO
  END DO
  Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, ChannelInfo, &
                                      RTSolution, RTSolution_TL )
  IF ( Error_Status /= SUCCESS ) THEN; WRITE(*,*) 'TL fail'; STOP 1; END IF
  LHS = ZERO
  CALL CRTM_RTSolution_Zero( RTSolution_AD )
  DO m = 1, N_PROFILES
    DO l = 1, ntot
      dy = RTSolution_TL(l,m)%Brightness_Temperature
      LHS = LHS + dy*dy
      RTSolution_AD(l,m)%Brightness_Temperature = dy
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
    END DO
  END DO
  rel = ABS(LHS-RHS) / MAX(ABS(LHS), TINY(ONE))
  CALL judge( 'adjoint dot-product closure (T+H2O, all five sensors)', &
              rel < 1.0e-10_fp )
  WRITE(*,'(9x,"rel = ",es10.3)') rel

  ! ---- K vs AD on the top line-splitting channel of each sounder ----
  DO i = 1, N_SENSORS
    SELECT CASE ( i )
      CASE ( S_AMSUA )  ; l0 = off(i) + 14   ! top line-splitting channel
      CASE ( S_MWTS3 )  ; l0 = off(i) + 17   ! top line-splitting channel
      CASE ( S_MWHS2 )  ; l0 = off(i) + 2    ! 118.75 GHz line center
      CASE ( S_MWRIRM ) ; l0 = off(i) + 1    ! 10.65 GHz window
      CASE ( S_GEMS2 )  ; l0 = off(i) + 15   ! 118.75 GHz line center
    END SELECT
    CALL CRTM_Atmosphere_Zero( Atm_AD ) ; CALL CRTM_Surface_Zero( Sfc_AD )
    CALL CRTM_RTSolution_Zero( RTSolution_AD )
    RTSolution_AD(l0,1)%Brightness_Temperature = ONE
    Error_Status = CRTM_Adjoint( Atm, Sfc, RTSolution_AD, Geometry, ChannelInfo, &
                                 Atm_AD, Sfc_AD, RTSolution )
    IF ( Error_Status /= SUCCESS ) THEN; WRITE(*,*) 'AD fail'; STOP 1; END IF
    pk = MAX( MAXVAL( ABS( Atm_K(l0,1)%Temperature - Atm_AD(1)%Temperature ) ), &
              MAXVAL( ABS( Atm_K(l0,1)%Absorber(:,IDX_H2O) - Atm_AD(1)%Absorber(:,IDX_H2O) ) ) )
    rel = pk / MAX( MAXVAL(ABS(Atm_K(l0,1)%Temperature)), &
                    MAXVAL(ABS(Atm_K(l0,1)%Absorber(:,IDX_H2O))), TINY(ONE) )
    CALL judge( TRIM(SENSORS(i))//' K == AD on probe channel', rel < 1.0e-9_fp )
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

  SUBROUTINE judge( name, okv )
    CHARACTER(*), INTENT(IN) :: name
    LOGICAL,      INTENT(IN) :: okv
    WRITE(*,'(5x,"[",a,"] ",a)') MERGE('PASS','FAIL',okv), name
    IF ( .NOT. okv ) all_ok = .FALSE.
  END SUBROUTINE judge

  ! Peak pressure (hPa) of the temperature Jacobian of sensor i, channel ch.
  REAL(fp) FUNCTION t_peak_hpa( i, ch )
    INTEGER, INTENT(IN) :: i, ch
    wf = Atm_K(off(i)+ch,1)%Temperature
    kpeak = MAXLOC( ABS(wf), DIM=1 )
    t_peak_hpa = Atm(1)%Pressure(kpeak)
  END FUNCTION t_peak_hpa

  ! Five-rung line-splitting ladder: each rung within LADDER_RTOL of the
  ! heritage value AND strictly descending. The double-offset conversion
  ! defect made rungs 2-5 near-equal at ~49 hPa; the multi-band convolution
  ! defect (#71) put them at 41-49 hPa; both fail here loudly.
  SUBROUTINE check_ladder( name, p )
    CHARACTER(*), INTENT(IN) :: name
    REAL(fp),     INTENT(IN) :: p(5)
    LOGICAL :: okv
    INTEGER :: k
    okv = .TRUE.
    DO k = 1, 5
      IF ( ABS(p(k) - LADDER(k)) > LADDER_RTOL * LADDER(k) ) okv = .FALSE.
      IF ( k > 1 ) THEN
        IF ( p(k) >= p(k-1) ) okv = .FALSE.
      END IF
    END DO
    CALL judge( name, okv )
    WRITE(*,'(9x,"peaks hPa: ",5f9.2,"   (heritage ",5f7.2,")")') p, LADDER
  END SUBROUTINE check_ladder

  INCLUDE 'Load_ECMWF84_Atm_Data.inc'

END PROGRAM test_MW_Sounder_Physics
