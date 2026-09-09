!
! test_Land_Jacobian
!
! Validation oracle for the analytic microwave LAND surface Jacobians
! (issue #281, Phases 1/2/3). For a pure-land MW scene it checks the analytic
! d(Tb)/d{LAI, Vegetation_Fraction, Soil_Moisture_Content, Soil_Temperature,
! Land_Temperature} obtained from both
!   * the K-matrix (adjoint path), and
!   * the tangent-linear model
! against central finite differences of the forward model. It also asserts that
! Canopy_Water_Content has an exactly-zero Jacobian (the LandEM forward never
! consumes it, so a valid analytic Jacobian is zero by construction).
!
! Temperature note: Land_Temperature carries TWO contributions -- the dominant
! skin-T Planck emission term (frequency-independent, via CRTM_Compute_SurfaceT_AD)
! plus the emissivity part through the LandEM thermal ratio gsect0. Soil_Temperature
! carries only the emissivity part (it never enters the surface emission), so it is
! the clean oracle for the Phase-3 dielectric+gsect0 derivative. The FD of the
! forward captures the full total in each case.
!
! A microwave window sensor is used (amsua_n19: channels 1-2 at 23.8/31.4 GHz
! sit below the 80 GHz cutoff of the NESDIS land emissivity model, so they are
! surface sensitive). Channels above the cutoff use a constant default
! emissivity, so both their analytic and finite-difference sensitivities are
! ~0 and must still agree.
!
! Exit status: STOP 0 = success, STOP 1 = failure.
!

PROGRAM test_Land_Jacobian

  ! ============================================================================
  USE CRTM_Module
  IMPLICIT NONE
  ! ============================================================================

  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME      = 'test_Land_Jacobian'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  CHARACTER(*), PARAMETER :: SENSOR_ID         = 'amsua_n19'

  INTEGER, PARAMETER :: N_PROFILES  = 2   ! matches Load_Atm_Data.inc
  INTEGER, PARAMETER :: N_LAYERS    = 92
  INTEGER, PARAMETER :: N_ABSORBERS = 2
  INTEGER, PARAMETER :: N_CLOUDS    = 0
  INTEGER, PARAMETER :: N_AEROSOLS  = 0
  INTEGER, PARAMETER :: N_SENSORS   = 1

  REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp

  ! Base land state (fractions kept inside (0,1) so they are not clipped)
  REAL(fp), PARAMETER :: LAI0   = 2.0_fp
  REAL(fp), PARAMETER :: VEG0   = 0.5_fp
  REAL(fp), PARAMETER :: SMC0   = 0.2_fp
  REAL(fp), PARAMETER :: TSOIL0 = 290.0_fp   ! in [100,350] -> not aliased to skin
  REAL(fp), PARAMETER :: TLAND0 = 290.0_fp
  ! Central finite-difference perturbations
  REAL(fp), PARAMETER :: DLAI = 1.0e-3_fp
  REAL(fp), PARAMETER :: DVEG = 1.0e-3_fp
  REAL(fp), PARAMETER :: DSMC = 1.0e-4_fp
  REAL(fp), PARAMETER :: DTS  = 1.0e-2_fp    ! soil/land temperature perturbation (K)
  ! Agreement tolerances: |analytic - FD| <= TOL_ABS + TOL_REL*|FD|
  REAL(fp), PARAMETER :: TOL_REL = 5.0e-3_fp
  REAL(fp), PARAMETER :: TOL_ABS = 1.0e-4_fp
  ! A channel counts as surface-sensitive when |FD| exceeds this (K per unit)
  REAL(fp), PARAMETER :: ACTIVE_THRESHOLD = 1.0e-2_fp

  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: Message, Version
  INTEGER :: Error_Status, Alloc_Status
  INTEGER :: n_Channels, l, m
  INTEGER :: n_active
  LOGICAL :: failed

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)    :: Geometry(N_PROFILES)

  TYPE(CRTM_Atmosphere_type)  :: Atm(N_PROFILES), Atm_TL(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc(N_PROFILES), Sfc_TL(N_PROFILES)

  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:), RTSolution_TL(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_K(:,:)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atmosphere_K(:,:)
  TYPE(CRTM_Surface_type)   , ALLOCATABLE :: Surface_K(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTS_p(:,:), RTS_m(:,:)

  ! Analytic and finite-difference Jacobians, dims (n_Channels, N_PROFILES)
  REAL(fp), ALLOCATABLE :: ad_lai(:,:), ad_veg(:,:), ad_smc(:,:)   ! K-matrix (adjoint)
  REAL(fp), ALLOCATABLE :: tl_lai(:,:), tl_veg(:,:), tl_smc(:,:)   ! tangent-linear
  REAL(fp), ALLOCATABLE :: fd_lai(:,:), fd_veg(:,:), fd_smc(:,:)   ! finite difference
  REAL(fp), ALLOCATABLE :: ad_tsoil(:,:), ad_tland(:,:), ad_cwc(:,:)
  REAL(fp), ALLOCATABLE :: tl_tsoil(:,:), tl_tland(:,:)
  REAL(fp), ALLOCATABLE :: fd_tsoil(:,:), fd_tland(:,:)


  ! Header
  CALL CRTM_Version( Version )
  CALL Program_Message( PROGRAM_NAME, &
    'Validate analytic MW land LAI/Vegetation_Fraction Jacobians vs finite differences.', &
    'CRTM Version: '//TRIM(Version) )


  ! 1. Initialize the CRTM
  ! ----------------------
  WRITE( *,'(/5x,"Initializing the CRTM (",a,")...")' ) SENSOR_ID
  Error_Status = CRTM_Init( (/SENSOR_ID/), ChannelInfo, File_Path=COEFFICIENTS_PATH )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error initializing CRTM', FAILURE ); STOP 1
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))


  ! 2. Allocate arrays and structures
  ! ---------------------------------
  ALLOCATE( RTSolution(n_Channels,N_PROFILES), RTSolution_TL(n_Channels,N_PROFILES), &
            RTSolution_K(n_Channels,N_PROFILES), Atmosphere_K(n_Channels,N_PROFILES), &
            Surface_K(n_Channels,N_PROFILES), RTS_p(n_Channels,N_PROFILES), &
            RTS_m(n_Channels,N_PROFILES), &
            ad_lai(n_Channels,N_PROFILES), ad_veg(n_Channels,N_PROFILES), ad_smc(n_Channels,N_PROFILES), &
            tl_lai(n_Channels,N_PROFILES), tl_veg(n_Channels,N_PROFILES), tl_smc(n_Channels,N_PROFILES), &
            fd_lai(n_Channels,N_PROFILES), fd_veg(n_Channels,N_PROFILES), fd_smc(n_Channels,N_PROFILES), &
            ad_tsoil(n_Channels,N_PROFILES), ad_tland(n_Channels,N_PROFILES), ad_cwc(n_Channels,N_PROFILES), &
            tl_tsoil(n_Channels,N_PROFILES), tl_tland(n_Channels,N_PROFILES), &
            fd_tsoil(n_Channels,N_PROFILES), fd_tland(n_Channels,N_PROFILES), &
            STAT = Alloc_Status )
  IF ( Alloc_Status /= 0 ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error allocating arrays', FAILURE ); STOP 1
  END IF

  CALL CRTM_Atmosphere_Create( Atm   , N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atm_TL, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atmosphere_K, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error allocating Atmosphere', FAILURE ); STOP 1
  END IF


  ! 3. Assign input data
  ! --------------------
  CALL Load_Atm_Data()       ! US standard atmosphere (fills Atm(1:2))
  CALL Load_Land_Surface()   ! pure-land surface for all profiles

  CALL CRTM_Geometry_SetValue( Geometry, &
                               Sensor_Zenith_Angle = ZENITH_ANGLE, &
                               Sensor_Scan_Angle   = SCAN_ANGLE )


  ! 4. Analytic Jacobian from the K-matrix (adjoint path)
  ! -----------------------------------------------------
  CALL CRTM_Atmosphere_Zero( Atmosphere_K )
  CALL CRTM_Surface_Zero( Surface_K )
  RTSolution_K%Radiance               = ZERO
  RTSolution_K%Brightness_Temperature = ONE

  Error_Status = CRTM_K_Matrix( Atm, Sfc, RTSolution_K, Geometry, ChannelInfo, &
                                Atmosphere_K, Surface_K, RTSolution )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error in CRTM K_Matrix', FAILURE ); STOP 1
  END IF
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      ad_lai(l,m)   = Surface_K(l,m)%Lai
      ad_veg(l,m)   = Surface_K(l,m)%Vegetation_Fraction
      ad_smc(l,m)   = Surface_K(l,m)%Soil_Moisture_Content
      ad_tsoil(l,m) = Surface_K(l,m)%Soil_Temperature
      ad_tland(l,m) = Surface_K(l,m)%Land_Temperature
      ad_cwc(l,m)   = Surface_K(l,m)%Canopy_Water_Content
    END DO
  END DO


  ! 5. Analytic Jacobian from the tangent-linear model
  ! --------------------------------------------------
  ! ...LAI direction
  CALL CRTM_Atmosphere_Zero( Atm_TL )
  CALL CRTM_Surface_Zero( Sfc_TL )
  Sfc_TL%Lai = ONE
  Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, ChannelInfo, &
                                      RTSolution, RTSolution_TL )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error in CRTM Tangent_Linear (LAI)', FAILURE ); STOP 1
  END IF
  tl_lai = RTSolution_TL%Brightness_Temperature

  ! ...Vegetation_Fraction direction
  CALL CRTM_Atmosphere_Zero( Atm_TL )
  CALL CRTM_Surface_Zero( Sfc_TL )
  Sfc_TL%Vegetation_Fraction = ONE
  Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, ChannelInfo, &
                                      RTSolution, RTSolution_TL )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error in CRTM Tangent_Linear (Veg)', FAILURE ); STOP 1
  END IF
  tl_veg = RTSolution_TL%Brightness_Temperature

  ! ...Soil_Moisture_Content direction
  CALL CRTM_Atmosphere_Zero( Atm_TL )
  CALL CRTM_Surface_Zero( Sfc_TL )
  Sfc_TL%Soil_Moisture_Content = ONE
  Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, ChannelInfo, &
                                      RTSolution, RTSolution_TL )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error in CRTM Tangent_Linear (SMC)', FAILURE ); STOP 1
  END IF
  tl_smc = RTSolution_TL%Brightness_Temperature

  ! ...Soil_Temperature direction
  CALL CRTM_Atmosphere_Zero( Atm_TL )
  CALL CRTM_Surface_Zero( Sfc_TL )
  Sfc_TL%Soil_Temperature = ONE
  Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, ChannelInfo, &
                                      RTSolution, RTSolution_TL )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error in CRTM Tangent_Linear (Tsoil)', FAILURE ); STOP 1
  END IF
  tl_tsoil = RTSolution_TL%Brightness_Temperature

  ! ...Land_Temperature direction (emission + emissivity parts)
  CALL CRTM_Atmosphere_Zero( Atm_TL )
  CALL CRTM_Surface_Zero( Sfc_TL )
  Sfc_TL%Land_Temperature = ONE
  Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, ChannelInfo, &
                                      RTSolution, RTSolution_TL )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error in CRTM Tangent_Linear (Tland)', FAILURE ); STOP 1
  END IF
  tl_tland = RTSolution_TL%Brightness_Temperature


  ! 6. Central finite differences of the forward model
  ! --------------------------------------------------
  ! ...LAI
  Sfc%Lai = LAI0 + DLAI
  CALL Run_Forward( RTS_p, 'LAI+' )
  Sfc%Lai = LAI0 - DLAI
  CALL Run_Forward( RTS_m, 'LAI-' )
  fd_lai = (RTS_p%Brightness_Temperature - RTS_m%Brightness_Temperature)/(2.0_fp*DLAI)
  Sfc%Lai = LAI0
  ! ...Vegetation_Fraction
  Sfc%Vegetation_Fraction = VEG0 + DVEG
  CALL Run_Forward( RTS_p, 'VEG+' )
  Sfc%Vegetation_Fraction = VEG0 - DVEG
  CALL Run_Forward( RTS_m, 'VEG-' )
  fd_veg = (RTS_p%Brightness_Temperature - RTS_m%Brightness_Temperature)/(2.0_fp*DVEG)
  Sfc%Vegetation_Fraction = VEG0
  ! ...Soil_Moisture_Content
  Sfc%Soil_Moisture_Content = SMC0 + DSMC
  CALL Run_Forward( RTS_p, 'SMC+' )
  Sfc%Soil_Moisture_Content = SMC0 - DSMC
  CALL Run_Forward( RTS_m, 'SMC-' )
  fd_smc = (RTS_p%Brightness_Temperature - RTS_m%Brightness_Temperature)/(2.0_fp*DSMC)
  Sfc%Soil_Moisture_Content = SMC0
  ! ...Soil_Temperature
  Sfc%Soil_Temperature = TSOIL0 + DTS
  CALL Run_Forward( RTS_p, 'Tsoil+' )
  Sfc%Soil_Temperature = TSOIL0 - DTS
  CALL Run_Forward( RTS_m, 'Tsoil-' )
  fd_tsoil = (RTS_p%Brightness_Temperature - RTS_m%Brightness_Temperature)/(2.0_fp*DTS)
  Sfc%Soil_Temperature = TSOIL0
  ! ...Land_Temperature
  Sfc%Land_Temperature = TLAND0 + DTS
  CALL Run_Forward( RTS_p, 'Tland+' )
  Sfc%Land_Temperature = TLAND0 - DTS
  CALL Run_Forward( RTS_m, 'Tland-' )
  fd_tland = (RTS_p%Brightness_Temperature - RTS_m%Brightness_Temperature)/(2.0_fp*DTS)
  Sfc%Land_Temperature = TLAND0


  ! 7. Compare
  ! ----------
  failed   = .FALSE.
  n_active = 0
  WRITE(*,'(/5x,a)') 'Per-channel comparison (K per unit). Columns: AD / TL / FD'
  WRITE(*,'(5x,a)')  '  m  ch       dTb/dLAI (AD/TL/FD)         dTb/dVeg (AD/TL/FD)        dTb/dSMC (AD/TL/FD)'
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      WRITE(*,'(5x,i3,i4,2x,3es12.3,2x,3es12.3,2x,3es12.3)') &
            m, RTSolution(l,m)%Sensor_Channel, &
            ad_lai(l,m), tl_lai(l,m), fd_lai(l,m), &
            ad_veg(l,m), tl_veg(l,m), fd_veg(l,m), &
            ad_smc(l,m), tl_smc(l,m), fd_smc(l,m)
      IF ( ABS(fd_lai(l,m)) > ACTIVE_THRESHOLD .AND. ABS(fd_smc(l,m)) > ACTIVE_THRESHOLD ) &
        n_active = n_active + 1
      CALL Check( 'AD dTb/dLAI', m, RTSolution(l,m)%Sensor_Channel, ad_lai(l,m), fd_lai(l,m), failed )
      CALL Check( 'TL dTb/dLAI', m, RTSolution(l,m)%Sensor_Channel, tl_lai(l,m), fd_lai(l,m), failed )
      CALL Check( 'AD dTb/dVeg', m, RTSolution(l,m)%Sensor_Channel, ad_veg(l,m), fd_veg(l,m), failed )
      CALL Check( 'TL dTb/dVeg', m, RTSolution(l,m)%Sensor_Channel, tl_veg(l,m), fd_veg(l,m), failed )
      CALL Check( 'AD dTb/dSMC', m, RTSolution(l,m)%Sensor_Channel, ad_smc(l,m), fd_smc(l,m), failed )
      CALL Check( 'TL dTb/dSMC', m, RTSolution(l,m)%Sensor_Channel, tl_smc(l,m), fd_smc(l,m), failed )
    END DO
  END DO

  ! ...temperatures (separate table) + Canopy_Water_Content structural zero
  WRITE(*,'(/5x,a)') 'Temperature Jacobians (K per K). Columns: AD / TL / FD'
  WRITE(*,'(5x,a)')  '  m  ch     dTb/dTsoil (AD/TL/FD)       dTb/dTland (AD/TL/FD)     K(Canopy_Water)'
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      WRITE(*,'(5x,i3,i4,2x,3es12.3,2x,3es12.3,2x,es11.3)') &
            m, RTSolution(l,m)%Sensor_Channel, &
            ad_tsoil(l,m), tl_tsoil(l,m), fd_tsoil(l,m), &
            ad_tland(l,m), tl_tland(l,m), fd_tland(l,m), ad_cwc(l,m)
      CALL Check( 'AD dTb/dTsoil', m, RTSolution(l,m)%Sensor_Channel, ad_tsoil(l,m), fd_tsoil(l,m), failed )
      CALL Check( 'TL dTb/dTsoil', m, RTSolution(l,m)%Sensor_Channel, tl_tsoil(l,m), fd_tsoil(l,m), failed )
      CALL Check( 'AD dTb/dTland', m, RTSolution(l,m)%Sensor_Channel, ad_tland(l,m), fd_tland(l,m), failed )
      CALL Check( 'TL dTb/dTland', m, RTSolution(l,m)%Sensor_Channel, tl_tland(l,m), fd_tland(l,m), failed )
      ! Canopy_Water_Content is never consumed by the LandEM forward -> Jacobian
      ! must be exactly zero (guarded so a future forward change that wires it in
      ! is flagged here rather than silently producing an unvalidated Jacobian).
      IF ( ABS(ad_cwc(l,m)) > TOL_ABS ) THEN
        WRITE(Message,'("Canopy_Water_Content Jacobian expected 0 but got ",es13.5, &
              &" (profile ",i0," channel ",i0,")")') ad_cwc(l,m), m, RTSolution(l,m)%Sensor_Channel
        CALL Display_Message( PROGRAM_NAME, TRIM(Message), FAILURE ); failed = .TRUE.
      END IF
    END DO
  END DO

  WRITE(*,'(/5x,"Surface-sensitive channel evaluations (|FD|>",es8.1,"): ",i0)') ACTIVE_THRESHOLD, n_active
  IF ( n_active < 1 ) THEN
    CALL Display_Message( PROGRAM_NAME, &
      'No surface-sensitive channels exercised -- test is not meaningful', FAILURE )
    failed = .TRUE.
  END IF


  ! 8. Clean up
  ! -----------
  Error_Status = CRTM_Destroy( ChannelInfo )
  CALL CRTM_Atmosphere_Destroy( Atm )
  CALL CRTM_Atmosphere_Destroy( Atm_TL )
  CALL CRTM_Atmosphere_Destroy( Atmosphere_K )
  DEALLOCATE( RTSolution, RTSolution_TL, RTSolution_K, Atmosphere_K, Surface_K, &
              RTS_p, RTS_m, ad_lai, ad_veg, ad_smc, tl_lai, tl_veg, tl_smc, &
              fd_lai, fd_veg, fd_smc, ad_tsoil, ad_tland, ad_cwc, &
              tl_tsoil, tl_tland, fd_tsoil, fd_tland )

  IF ( failed ) THEN
    CALL Display_Message( PROGRAM_NAME, 'FAILED: analytic Jacobians disagree with finite differences', FAILURE )
    STOP 1
  ELSE
    CALL Display_Message( PROGRAM_NAME, 'PASSED: analytic Jacobians match finite differences', INFORMATION )
    STOP 0
  END IF


CONTAINS


  ! Run the forward model into the supplied RTSolution array
  SUBROUTINE Run_Forward( RTS, label )
    TYPE(CRTM_RTSolution_type), INTENT(IN OUT) :: RTS(:,:)
    CHARACTER(*),               INTENT(IN)     :: label
    INTEGER :: stat
    stat = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTS )
    IF ( stat /= SUCCESS ) THEN
      CALL Display_Message( PROGRAM_NAME, 'Error in CRTM Forward ('//label//')', FAILURE ); STOP 1
    END IF
  END SUBROUTINE Run_Forward


  ! Compare one analytic value against its finite-difference reference
  SUBROUTINE Check( name, m, ch, analytic, fd, failed )
    CHARACTER(*), INTENT(IN)    :: name
    INTEGER,      INTENT(IN)    :: m, ch
    REAL(fp),     INTENT(IN)    :: analytic, fd
    LOGICAL,      INTENT(IN OUT):: failed
    REAL(fp) :: tol
    tol = TOL_ABS + TOL_REL*ABS(fd)
    IF ( ABS(analytic - fd) > tol ) THEN
      WRITE(Message,'(a," mismatch: profile ",i0," channel ",i0,&
            &": analytic=",es13.5," FD=",es13.5," |diff|=",es11.3," tol=",es11.3)') &
            TRIM(name), m, ch, analytic, fd, ABS(analytic-fd), tol
      CALL Display_Message( PROGRAM_NAME, TRIM(Message), FAILURE )
      failed = .TRUE.
    END IF
  END SUBROUTINE Check


  ! Pure-land surface for all profiles
  SUBROUTINE Load_Land_Surface()
    INTEGER :: mm
    DO mm = 1, N_PROFILES
      Sfc(mm)%Land_Coverage        = 1.0_fp
      Sfc(mm)%Water_Coverage       = 0.0_fp
      Sfc(mm)%Snow_Coverage        = 0.0_fp
      Sfc(mm)%Ice_Coverage         = 0.0_fp
      Sfc(mm)%Land_Type            = 1            ! valid NPOESS land type (IR/VIS only)
      Sfc(mm)%Soil_Type            = 1            ! COARSE          (MW land model)
      Sfc(mm)%Vegetation_Type      = 7            ! GROUNDCOVER     (MW land model)
      Sfc(mm)%Land_Temperature     = TLAND0
      Sfc(mm)%Soil_Temperature     = TSOIL0
      Sfc(mm)%Soil_Moisture_Content= SMC0
      Sfc(mm)%Lai                  = LAI0
      Sfc(mm)%Vegetation_Fraction  = VEG0
    END DO
  END SUBROUTINE Load_Land_Surface


  INCLUDE 'Load_Atm_Data.inc'

END PROGRAM test_Land_Jacobian
