!
! test_Emissivity_Jacobian
!
! Validate the surface-emissivity Jacobian output d(Tb)/d(emissivity) in
! RTSolution_K%Surface_Emissivity against central finite differences of the
! forward model, on both branches of Assign_Common_Output_AD:
!
!   1. Computed-emissivity branch (TELSEM2 atlas active, clear-sky emission RT):
!      the G2 fix folds the reflected-downwelling term (via the specular
!      closure r = 1-e) into the reported value. The finite-difference
!      reference perturbs Options%Emissivity around the atlas emissivity --
!      for the single-stream emission solver the user-emissivity forward
!      reproduces the atlas forward exactly, so FD and K describe the same
!      derivative and must agree tightly.
!
!   2. User-emissivity branch, clear sky (the branch historically flagged
!      "!! need to check !!" -- the G3 validation).
!
!   3. User-emissivity branch under scattering (overcast cloud, ADA): the G1
!      fix restores the emission term that was hardcoded to zero since V2.4.
!
!   4. Same under the SOI solver, all channels: the G6 fix zeroes the SOI
!      adjoint surface duals on entry; without it, residue from channel n's
!      seed leaks into channel n+1's Jacobian (the user-emissivity specular
!      branch reads but never zeroes the reflectivity adjoint diagonals).
!
!   5. Fractional cloud (profile 2, Cloud_Fraction = 0.5): the G4 fix sums
!      the coverage-weighted clear- and cloudy-column captures, so the
!      reported d(Tb)/d(emissivity) is the total-column derivative and must
!      match the finite difference of the (1-cc)/cc-combined forward.
!      Before the fix the value was cloudy-column-only and fails this test.
!
PROGRAM test_Emissivity_Jacobian

  USE CRTM_Module
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME      = 'test_Emissivity_Jacobian'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  CHARACTER(*), PARAMETER :: ATLAS_FILE        = 'TELSEM2.MWland.test.nc'  ! test-only staged name
  CHARACTER(*), PARAMETER :: SENSOR_ID         = 'amsua_n19'

  INTEGER,  PARAMETER :: N_PROFILES  = 2   ! matches Load_Atm_Data.inc
  INTEGER,  PARAMETER :: N_LAYERS    = 92
  INTEGER,  PARAMETER :: N_ABSORBERS = 2
  INTEGER,  PARAMETER :: N_SENSORS   = 1

  REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp
  ! Land cell with TELSEM2 climatology (northern Argentina), September
  REAL(fp), PARAMETER :: LAT_A = -30.0_fp, LON_A = 302.0_fp
  INTEGER,  PARAMETER :: MON_SEP = 9

  ! Central-difference step in emissivity
  REAL(fp), PARAMETER :: DELTA = 1.0E-3_fp
  ! Relative tolerances: the clear-sky emission radiance is exactly linear in
  ! emissivity, so FD truncation error is set only by the Planck inversion --
  ! tight. Scattering solvers add iteration/quadrature noise -- looser.
  REAL(fp), PARAMETER :: TOL_CLEAR = 1.0E-4_fp
  REAL(fp), PARAMETER :: TOL_SCAT  = 1.0E-3_fp
  ! Sanity: user-emissivity forward at the atlas emissivity must reproduce
  ! the atlas forward (same state, same solver inputs).
  REAL(fp), PARAMETER :: TOL_TB    = 1.0E-6_fp
  ! Base emissivity for the pure user-emissivity scattering sections
  REAL(fp), PARAMETER :: EMIS0 = 0.9_fp

  CHARACTER(256) :: Message, Version
  INTEGER  :: Error_Status, Alloc_Status, n_Channels, l
  LOGICAL  :: failed

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)    :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc(N_PROFILES)
  TYPE(CRTM_Options_type)     :: Opt(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:), RTSolution_K(:,:)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atmosphere_K(:,:)
  TYPE(CRTM_Surface_type),    ALLOCATABLE :: Surface_K(:,:)

  REAL(fp), ALLOCATABLE :: emis_b(:), tb0(:), tbchk(:), tbp(:), tbm(:)
  REAL(fp), ALLOCATABLE :: fd(:), k_atlas(:), k_user(:), eps(:)
  REAL(fp), ALLOCATABLE :: tbp2(:), tbm2(:), fd2(:), k_user2(:)

  CALL CRTM_Version( Version )
  CALL Program_Message( PROGRAM_NAME, &
    'Validate d(Tb)/d(emissivity) in RTSolution_K%Surface_Emissivity '//&
    'against central finite differences.', &
    'CRTM Version: '//TRIM(Version) )

  failed = .FALSE.

  Error_Status = CRTM_Init( (/SENSOR_ID/), ChannelInfo, &
                            File_Path=COEFFICIENTS_PATH, MWlandCoeff_File=ATLAS_FILE )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error initializing CRTM with atlas', FAILURE ); STOP 1
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))

  ALLOCATE( RTSolution(n_Channels,N_PROFILES), RTSolution_K(n_Channels,N_PROFILES), &
            Atmosphere_K(n_Channels,N_PROFILES), Surface_K(n_Channels,N_PROFILES), &
            emis_b(n_Channels), tb0(n_Channels), tbchk(n_Channels), &
            tbp(n_Channels), tbm(n_Channels), fd(n_Channels), &
            k_atlas(n_Channels), k_user(n_Channels), eps(n_Channels), &
            tbp2(n_Channels), tbm2(n_Channels), fd2(n_Channels), &
            k_user2(n_Channels), &
            STAT = Alloc_Status )
  IF ( Alloc_Status /= 0 ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error allocating arrays', FAILURE ); STOP 1
  END IF
  CALL CRTM_RTSolution_Create( RTSolution,   N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_K, N_LAYERS )
  CALL CRTM_Options_Create( Opt, n_Channels )

  ! Clear-sky atmosphere and pure-land surface
  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, 0, 0 )
  CALL CRTM_Atmosphere_Create( Atmosphere_K, N_LAYERS, N_ABSORBERS, 0, 0 )
  CALL Load_Atm_Data()
  CALL Load_Land_Surface()
  CALL CRTM_Geometry_SetValue( Geometry, &
                               Sensor_Zenith_Angle = ZENITH_ANGLE, &
                               Sensor_Scan_Angle   = SCAN_ANGLE, &
                               Latitude  = LAT_A, Longitude = LON_A, Month = MON_SEP )

  ! ------------------------------------------------------------------
  ! 1. Computed-emissivity branch (atlas), clear-sky emission RT
  ! ------------------------------------------------------------------
  CALL Run_Forward( tb0, use_emis=.FALSE. )
  emis_b = RTSolution(:,1)%Surface_Emissivity

  CALL Run_K( k_atlas, use_emis=.FALSE. )

  ! User-emissivity forward at the atlas emissivity must reproduce the run
  CALL Run_Forward( tbchk, use_emis=.TRUE., emis=emis_b )
  CALL Assert_Close( tbchk, tb0, TOL_TB, 'sanity: user-emissivity forward reproduces atlas Tb', &
                     absolute=.TRUE. )

  ! Central difference around the atlas emissivity
  CALL Run_Forward( tbp, use_emis=.TRUE., emis=emis_b+DELTA )
  CALL Run_Forward( tbm, use_emis=.TRUE., emis=emis_b-DELTA )
  fd = ( tbp - tbm ) / ( 2.0_fp*DELTA )

  CALL Assert_Close( k_atlas, fd, TOL_CLEAR, 'G2: computed-path (atlas) K vs FD, clear' )

  ! ------------------------------------------------------------------
  ! 2. User-emissivity branch, clear sky (G3, the "need to check" branch)
  ! ------------------------------------------------------------------
  CALL Run_K( k_user, use_emis=.TRUE., emis=emis_b )
  CALL Assert_Close( k_user, fd, TOL_CLEAR, 'G3: user-emissivity K vs FD, clear' )

  ! ------------------------------------------------------------------
  ! 3. User-emissivity branch under scattering (overcast cloud, ADA) -- G1
  ! ------------------------------------------------------------------
  CALL CRTM_Atmosphere_Destroy( Atm )
  CALL CRTM_Atmosphere_Destroy( Atmosphere_K )
  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, 1, 1 )
  CALL CRTM_Atmosphere_Create( Atmosphere_K, N_LAYERS, N_ABSORBERS, 1, 1 )
  CALL Load_Atm_Data()   ! profile 1: overcast (Cloud_Fraction=1) snow cloud

  eps = EMIS0
  CALL Run_Forward( tbp, use_emis=.TRUE., emis=eps+DELTA )
  CALL Run_Forward( tbm, use_emis=.TRUE., emis=eps-DELTA )
  fd = ( tbp - tbm ) / ( 2.0_fp*DELTA )
  CALL Run_K( k_user, use_emis=.TRUE., emis=eps )
  CALL Assert_Close( k_user, fd, TOL_SCAT, 'G1: user-emissivity K vs FD, overcast ADA' )

  ! ------------------------------------------------------------------
  ! 4. Same under the SOI solver, all channels -- G6 (channel carryover)
  ! ------------------------------------------------------------------
  Opt(:)%RT_Algorithm_ID = RT_SOI
  CALL Run_Forward( tbp, use_emis=.TRUE., emis=eps+DELTA )
  CALL Run_Forward( tbm, use_emis=.TRUE., emis=eps-DELTA )
  fd = ( tbp - tbm ) / ( 2.0_fp*DELTA )
  CALL Run_K( k_user, use_emis=.TRUE., emis=eps )
  CALL Assert_Close( k_user, fd, TOL_SCAT, 'G6: user-emissivity K vs FD, overcast SOI, all channels' )
  Opt(:)%RT_Algorithm_ID = RT_ADA

  ! ------------------------------------------------------------------
  ! 5. Fractional cloud (profile 2, Cloud_Fraction = 0.5) -- G4: the
  !    reported d(Tb)/d(emissivity) must be the total-column derivative
  ! ------------------------------------------------------------------
  CALL Run_Forward( tbp, use_emis=.TRUE., emis=eps+DELTA, tb2=tbp2 )
  CALL Run_Forward( tbm, use_emis=.TRUE., emis=eps-DELTA, tb2=tbm2 )
  fd2 = ( tbp2 - tbm2 ) / ( 2.0_fp*DELTA )
  CALL Run_K( k_user, use_emis=.TRUE., emis=eps, k2=k_user2 )
  CALL Assert_Close( k_user2, fd2, TOL_SCAT, 'G4: user-emissivity K vs FD, fractional cloud (cc=0.5)' )

  ! ------------------------------------------------------------------
  ! Report and clean up
  ! ------------------------------------------------------------------
  Error_Status = CRTM_Destroy( ChannelInfo )
  CALL CRTM_Atmosphere_Destroy( Atm )
  CALL CRTM_Atmosphere_Destroy( Atmosphere_K )
  CALL CRTM_Options_Destroy( Opt )
  DEALLOCATE( RTSolution, RTSolution_K, Atmosphere_K, Surface_K, &
              emis_b, tb0, tbchk, tbp, tbm, fd, k_atlas, k_user, eps, &
              tbp2, tbm2, fd2, k_user2 )

  IF ( failed ) THEN
    CALL Display_Message( PROGRAM_NAME, 'FAILED', FAILURE ); STOP 1
  ELSE
    CALL Display_Message( PROGRAM_NAME, &
      'PASSED: d(Tb)/d(emissivity) matches finite differences on the '//&
      'computed (atlas) and user-emissivity branches, clear and scattering, '//&
      'ADA and SOI, overcast and fractional cloud', &
      INFORMATION ); STOP 0
  END IF

CONTAINS

  ! Forward run; optionally with user-specified emissivity (same value vector
  ! for every profile). Returns profile 1's brightness temperatures, and
  ! optionally profile 2's (the fractional-cloud profile).
  SUBROUTINE Run_Forward( tb, use_emis, emis, tb2 )
    REAL(fp),           INTENT(OUT) :: tb(:)
    LOGICAL,            INTENT(IN)  :: use_emis
    REAL(fp), OPTIONAL, INTENT(IN)  :: emis(:)
    REAL(fp), OPTIONAL, INTENT(OUT) :: tb2(:)
    INTEGER :: stat, m
    DO m = 1, N_PROFILES
      Opt(m)%Use_Emissivity = use_emis
      IF ( use_emis ) CALL CRTM_Options_SetEmissivity( Opt(m), emis )
    END DO
    stat = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution, Options=Opt )
    IF ( stat /= SUCCESS ) THEN
      CALL Display_Message( PROGRAM_NAME, 'Error in CRTM Forward', FAILURE ); STOP 1
    END IF
    tb = RTSolution(:,1)%Brightness_Temperature
    IF ( PRESENT(tb2) ) tb2 = RTSolution(:,2)%Brightness_Temperature
  END SUBROUTINE Run_Forward

  ! K-matrix run with a unit brightness-temperature seed; returns profile 1's
  ! d(Tb)/d(emissivity) from RTSolution_K%Surface_Emissivity, and optionally
  ! profile 2's.
  SUBROUTINE Run_K( k, use_emis, emis, k2 )
    REAL(fp),           INTENT(OUT) :: k(:)
    LOGICAL,            INTENT(IN)  :: use_emis
    REAL(fp), OPTIONAL, INTENT(IN)  :: emis(:)
    REAL(fp), OPTIONAL, INTENT(OUT) :: k2(:)
    INTEGER :: stat, m
    DO m = 1, N_PROFILES
      Opt(m)%Use_Emissivity = use_emis
      IF ( use_emis ) CALL CRTM_Options_SetEmissivity( Opt(m), emis )
    END DO
    CALL CRTM_Atmosphere_Zero( Atmosphere_K )
    CALL CRTM_Surface_Zero( Surface_K )
    RTSolution_K%Radiance               = ZERO
    RTSolution_K%Brightness_Temperature = ONE
    stat = CRTM_K_Matrix( Atm, Sfc, RTSolution_K, Geometry, ChannelInfo, &
                          Atmosphere_K, Surface_K, RTSolution, Options=Opt )
    IF ( stat /= SUCCESS ) THEN
      CALL Display_Message( PROGRAM_NAME, 'Error in CRTM K_Matrix', FAILURE ); STOP 1
    END IF
    k = RTSolution_K(:,1)%Surface_Emissivity
    IF ( PRESENT(k2) ) k2 = RTSolution_K(:,2)%Surface_Emissivity
  END SUBROUTINE Run_K

  ! Assert per-channel closeness, relative by default (denominator floored at
  ! 1 K/unit-emissivity so near-zero-sensitivity channels are judged absolutely).
  SUBROUTINE Assert_Close( a, b, tol, label, absolute )
    REAL(fp),          INTENT(IN) :: a(:), b(:), tol
    CHARACTER(*),      INTENT(IN) :: label
    LOGICAL, OPTIONAL, INTENT(IN) :: absolute
    REAL(fp) :: err, maxerr
    LOGICAL  :: abs_mode
    INTEGER  :: lc
    abs_mode = .FALSE.
    IF ( PRESENT(absolute) ) abs_mode = absolute
    maxerr = ZERO
    DO lc = 1, SIZE(a)
      IF ( abs_mode ) THEN
        err = ABS( a(lc) - b(lc) )
      ELSE
        err = ABS( a(lc) - b(lc) ) / MAX( ABS(b(lc)), ONE )
      END IF
      maxerr = MAX( maxerr, err )
    END DO
    WRITE(*,'(5x,"Max err for ",a," = ",es12.4," (tol ",es9.2,")")') label, maxerr, tol
    IF ( maxerr > tol ) THEN
      WRITE(Message,'("Mismatch: ",a," max err ",es12.4," exceeds ",es9.2)') label, maxerr, tol
      CALL Display_Message( PROGRAM_NAME, TRIM(Message), FAILURE )
      failed = .TRUE.
      WRITE(*,'(7x,"chan",12x,"K",18x,"ref")')
      DO lc = 1, SIZE(a)
        WRITE(*,'(7x,i4,2es19.10)') lc, a(lc), b(lc)
      END DO
    END IF
  END SUBROUTINE Assert_Close

  ! Pure-land surface for all profiles (matches test_TELSEM2_MWland)
  SUBROUTINE Load_Land_Surface()
    INTEGER :: mm
    DO mm = 1, N_PROFILES
      Sfc(mm)%Land_Coverage         = 1.0_fp
      Sfc(mm)%Water_Coverage        = 0.0_fp
      Sfc(mm)%Snow_Coverage         = 0.0_fp
      Sfc(mm)%Ice_Coverage          = 0.0_fp
      Sfc(mm)%Land_Type             = 1
      Sfc(mm)%Soil_Type             = 1
      Sfc(mm)%Vegetation_Type       = 7
      Sfc(mm)%Land_Temperature      = 290.0_fp
      Sfc(mm)%Soil_Temperature      = 290.0_fp
      Sfc(mm)%Soil_Moisture_Content = 0.2_fp
      Sfc(mm)%Lai                   = 2.0_fp
      Sfc(mm)%Vegetation_Fraction   = 0.5_fp
    END DO
  END SUBROUTINE Load_Land_Surface

  INCLUDE 'Load_Atm_Data.inc'

END PROGRAM test_Emissivity_Jacobian
