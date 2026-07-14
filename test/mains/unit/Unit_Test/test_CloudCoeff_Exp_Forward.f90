!
! test_CloudCoeff_Exp_Forward
!
! Forward-mode coverage for the experimental ('CRTM-Exp') cloud-optics scheme.
!
! Initializes CRTM with Cloud_Model='CRTM-Exp' and the complete 6-habit
! experimental LUT (CloudCoeff_Exp_Full6.nc), then runs the mwr_aws microwave
! sensor over an ocean US-Standard column for three profiles:
!     1 = clear, 2 = thin graupel, 3 = heavy graupel
! and asserts that the experimental scattering path is physical:
!     * the forward model runs and returns physical TBs for EVERY channel and
!       profile (incl. optically THIN frozen cloud, WC = 0.05 kg/m^2/layer),
!     * a graupel cloud produces a brightness-temperature DEPRESSION, and
!     * that depression GROWS with water content.
!
! The thin profile is included deliberately: it is the regime that exposed a
! surface-reflectivity bug (the FASTEM-fit reflection correction extrapolating
! to a non-physical value at the near-grazing Gaussian quadrature angles the
! scattering RT uses, >= 200 GHz / PARMIO) which produced -1e15 K TBs at the
! 325 GHz AWS sideband channels. Guarded in CRTM_PARMIO (reflectivity clamped
! to [0,1]); this test locks that in.
!
! STOP 0 on success, STOP 1 on failure.
!
! CREATION HISTORY:
!       Written by:     Benjamin Johnson, 05-Jun-2026
!                       Adapted from the exp_aws_scatter validation driver.
!

PROGRAM test_CloudCoeff_Exp_Forward

  USE CRTM_Module
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_CloudCoeff_Exp_Forward'
  CHARACTER(*), PARAMETER :: PATH   = './testinput/'
  CHARACTER(*), PARAMETER :: SENSOR = 'mwr_aws'
  CHARACTER(*), PARAMETER :: LUT    = 'CloudCoeff_Exp_Full6.nc'

  ! Profile / column setup
  INTEGER,  PARAMETER :: N_PROFILES  = 3       ! 1=clear, 2=thin graupel, 3=heavy graupel
  INTEGER,  PARAMETER :: N_LAYERS    = 100
  INTEGER,  PARAMETER :: N_ABSORBERS = 6
  INTEGER,  PARAMETER :: N_CLOUDS    = 1
  INTEGER,  PARAMETER :: N_AEROSOLS  = 0
  REAL(fp), PARAMETER :: ZENITH      = 53.0_fp  ! AWS conical scan ~53 deg
  INTEGER,  PARAMETER :: KC1 = 78, KC2 = 86     ! cloud vertical band (layers)
  REAL(fp), PARAMETER :: REFF_G      = 500.0_fp ! graupel effective radius (microns)
  REAL(fp), PARAMETER :: WC_THIN      = 0.05_fp  ! kg/m^2 per layer (optically thin)
  REAL(fp), PARAMETER :: WC_HEAVY    = 1.00_fp  ! kg/m^2 per layer (heavy)

  ! Pass/fail thresholds (conservative; the path saturates well above these)
  REAL(fp), PARAMETER :: MIN_DEPRESSION = 1.0_fp   ! K
  REAL(fp), PARAMETER :: TB_LO = 50.0_fp, TB_HI = 330.0_fp

  TYPE(CRTM_ChannelInfo_type)             :: chinfo(1)
  TYPE(CRTM_Geometry_type)                :: geo(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:)

  INTEGER  :: err, nch, m
  REAL(fp) :: dep_thin, dep_heavy, tb_min, tb_max
  LOGICAL  :: ok

  ok = .TRUE.

  ! --------------------------------------------------------------------------
  ! Initialize CRTM with the experimental cloud-optics scheme
  ! --------------------------------------------------------------------------
  err = CRTM_Init( (/ SENSOR /), chinfo,           &
                   Cloud_Model       = 'CRTM-Exp', &
                   CloudCoeff_File   = LUT,        &
                   CloudCoeff_Format = 'netCDF',   &
                   File_Path         = PATH,       &
                   Quiet             = .TRUE. )
  IF ( err /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Init (Cloud_Model=CRTM-Exp) failed', FAILURE )
    STOP 1
  END IF
  nch = SUM( CRTM_ChannelInfo_n_Channels(chinfo) )
  IF ( nch < 1 ) THEN
    CALL Display_Message( PROGRAM_NAME, 'no channels loaded for '//SENSOR, FAILURE )
    STOP 1
  END IF

  ! --------------------------------------------------------------------------
  ! Build the atmosphere/surface/geometry
  ! --------------------------------------------------------------------------
  ALLOCATE( rts(nch, N_PROFILES) )
  CALL CRTM_Atmosphere_Create( atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Atmosphere_Create failed', FAILURE )
    STOP 1
  END IF

  CALL Load_ECMWF84_Atm_Data()        ! fills atm(1) (US-Standard ocean column)
  DO m = 2, N_PROFILES                 ! identical base column for every profile
    atm(m) = atm(1)
  END DO

  ! Profile 1: clear sky
  atm(1)%n_Clouds       = 0
  atm(1)%Cloud_Fraction = ZERO
  ! Profiles 2 & 3: graupel cloud in the band (thin / heavy loading)
  CALL Set_Graupel( atm(2), WC_THIN )
  CALL Set_Graupel( atm(3), WC_HEAVY )

  ! Ocean surface + geometry, identical for all profiles
  DO m = 1, N_PROFILES
    sfc(m)%Water_Coverage    = 1.0_fp
    sfc(m)%Water_Type        = 1          ! SEA_WATER
    sfc(m)%Water_Temperature = 290.0_fp
    sfc(m)%Wind_Speed        = 6.0_fp
    sfc(m)%Salinity          = 33.0_fp
    CALL CRTM_Geometry_SetValue( geo(m), Sensor_Zenith_Angle = ZENITH )
  END DO

  ! --------------------------------------------------------------------------
  ! Forward model
  ! --------------------------------------------------------------------------
  err = CRTM_Forward( atm, sfc, geo, chinfo, rts )
  IF ( err /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Forward failed', FAILURE )
    STOP 1
  END IF

  ! --------------------------------------------------------------------------
  ! Physical checks
  ! --------------------------------------------------------------------------
  ! Every TB (all profiles, all channels) must be physical
  tb_min = MINVAL( rts%Brightness_Temperature )
  tb_max = MAXVAL( rts%Brightness_Temperature )
  ! Max TB depression (clear - cloudy) over all channels, per loading
  dep_thin   = MAXVAL( rts(:,1)%Brightness_Temperature - rts(:,2)%Brightness_Temperature )
  dep_heavy = MAXVAL( rts(:,1)%Brightness_Temperature - rts(:,3)%Brightness_Temperature )

  WRITE(*,'(/a)')       ' CRTM-Exp forward scattering check (mwr_aws, graupel over ocean):'
  WRITE(*,'(a,i0)')     '   channels                 : ', nch
  WRITE(*,'(a,2f8.2)')  '   all-profile TB range (K) : ', tb_min, tb_max
  WRITE(*,'(a,f8.3)')   '   max depression, thin     : ', dep_thin
  WRITE(*,'(a,f8.3)')   '   max depression, heavy    : ', dep_heavy

  ! 1) all radiances physical (forward model produced sensible TBs everywhere)
  IF ( tb_min < TB_LO .OR. tb_max > TB_HI ) THEN
    WRITE(*,'(a,2f10.2)') ' FAIL: a TB is outside the physical range ', tb_min, tb_max
    ok = .FALSE.
  END IF
  ! 2) a graupel cloud scatters -> measurable cold depression
  IF ( dep_heavy < MIN_DEPRESSION ) THEN
    WRITE(*,'(a,f8.3,a,f6.2,a)') ' FAIL: heavy-graupel depression ', dep_heavy, &
         ' K is below the ', MIN_DEPRESSION, ' K threshold'
    ok = .FALSE.
  END IF
  ! 3) depression grows with water content (monotonic scattering response)
  IF ( dep_heavy <= dep_thin ) THEN
    WRITE(*,'(a)') ' FAIL: depression did not increase with water content'
    ok = .FALSE.
  END IF

  ! --------------------------------------------------------------------------
  ! Clean up + verdict
  ! --------------------------------------------------------------------------
  DEALLOCATE( rts )
  CALL CRTM_Atmosphere_Destroy( atm )
  err = CRTM_Destroy( chinfo )

  IF ( ok ) THEN
    WRITE(*,'(/a)') ' PASS: CRTM-Exp cloud-optics forward path produces physical scattering.'
    STOP 0
  ELSE
    WRITE(*,'(/a)') ' FAIL: CRTM-Exp forward checks failed.'
    STOP 1
  END IF

CONTAINS

  ! Put a graupel cloud of the given per-layer water content into the band.
  SUBROUTINE Set_Graupel( a, wc )
    TYPE(CRTM_Atmosphere_type), INTENT(IN OUT) :: a
    REAL(fp),                   INTENT(IN)     :: wc
    a%n_Clouds                            = 1
    a%Cloud_Fraction                      = ZERO
    a%Cloud_Fraction(KC1:KC2)             = 1.0_fp
    a%Cloud(1)%Type                       = GRAUPEL_CLOUD
    a%Cloud(1)%Effective_Radius           = ZERO
    a%Cloud(1)%Water_Content              = ZERO
    a%Cloud(1)%Effective_Radius(KC1:KC2)  = REFF_G
    a%Cloud(1)%Water_Content(KC1:KC2)     = wc
  END SUBROUTINE Set_Graupel

  INCLUDE 'Load_ECMWF84_Atm_Data.inc'

END PROGRAM test_CloudCoeff_Exp_Forward
