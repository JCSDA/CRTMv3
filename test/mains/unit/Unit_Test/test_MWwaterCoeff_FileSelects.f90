!
! test_MWwaterCoeff_FileSelects
!
! Proves that MWwaterCoeff_File actually selects the microwave water emissivity
! model, and that a caller asking for one model does not silently receive
! another.
!
! The defect
! ----------
! The file-based MWwaterCoeff load in CRTM_LifeCycle is commented out; it was
! only ever needed for the FASTEM5 binary lookup tables. The model is chosen by
! the scheme string instead. That left MWwaterCoeff_File accepted, stored,
! echoed in the "Loading MW water emissivity coefficients" message, and
! otherwise ignored, so CRTM_Init(MWwaterCoeff_File='FASTEM4...') returned
! SUCCESS with FASTEM6 loaded and said nothing.
!
! That is not a harmless no-op. FASTEM6 has no third or fourth Stokes azimuth
! model and returns U and V as identically zero, while FASTEM4 carries all
! four. A polarimetric caller who selected a polarimetric backend therefore got
! an unpolarised surface, and U = 0 is indistinguishable from a scene with no
! polarimetric signal. JEDI/UFO reaches CRTM through exactly this argument: it
! builds TRIM(MWwaterCoeff)//".MWwater.EmisCoeff.nc" from its own yaml key, so
! a JEDI user writing "MWwaterCoeff: FASTEM4" was silently given FASTEM6.
!
! What this test does
! -------------------
! It uses the one unambiguous observable that separates the two models: over
! ocean at a nonzero relative azimuth, FASTEM4 produces a nonzero third and
! fourth Stokes emissivity and FASTEM6 produces exactly zero.
!
!   1. CRTM_Init with MWwaterCoeff_File naming FASTEM4 and no scheme argument,
!      then assert U and V are nonzero. This is the assertion that fails
!      against the unfixed code, where FASTEM6 is loaded and both are exactly
!      zero.
!   2. CRTM_Init with MWwaterCoeff_File naming FASTEM6, then assert U and V are
!      exactly zero. Without this the test could pass by always loading
!      FASTEM4, so it is what makes step 1 evidence of selection rather than of
!      a new hardcoded default.
!
! Note the filename never has to exist on disk for this to work, and that is
! the point rather than an oversight: nothing opens it. It is a model selector
! wearing a filename, which is why the argument was so easy to ignore.
!

PROGRAM test_MWwaterCoeff_FileSelects

  ! -----------------
  ! Environment setup
  ! -----------------
  USE CRTM_Module
  USE CRTM_SfcOptics_Define   , ONLY: CRTM_SfcOptics_type      , &
                                      CRTM_SfcOptics_Create    , &
                                      CRTM_SfcOptics_Destroy   , &
                                      CRTM_SfcOptics_Associated
  USE CRTM_SfcOptics          , ONLY: CRTM_Compute_SfcOptics, iVar_type
  USE CRTM_GeometryInfo_Define, ONLY: CRTM_GeometryInfo_type, &
                                      CRTM_GeometryInfo_SetValue
  USE CRTM_GeometryInfo       , ONLY: CRTM_GeometryInfo_Compute
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_MWwaterCoeff_FileSelects'
  ! amsua_n19 channel 1 is 23.8 GHz, below the PARMIO frequency gate, so the
  ! microwave water dispatcher uses the FASTEM backend under test.
  CHARACTER(*), PARAMETER :: SENSORS(1) = (/ 'amsua_n19' /)
  CHARACTER(*), PARAMETER :: PATH       = 'testinput/'
  INTEGER     , PARAMETER :: CHANNEL    = 1

  CHARACTER(*), PARAMETER :: FILE_F4 = 'FASTEM4.MWwater.EmisCoeff.nc'
  CHARACTER(*), PARAMETER :: FILE_F6 = 'FASTEM6.MWwater.EmisCoeff.nc'

  ! Ocean, brisk wind, well off nadir and away from the azimuths where the odd
  ! harmonics vanish, so a real FASTEM4 signal is comfortably above round-off.
  REAL(fp), PARAMETER :: WIND_SPEED = 12.0_fp
  REAL(fp), PARAMETER :: WATER_TEMP = 285.0_fp
  REAL(fp), PARAMETER :: SALINITY   = 33.0_fp
  REAL(fp), PARAMETER :: ZENITH     = 45.0_fp
  REAL(fp), PARAMETER :: WIND_DIR   = 100.0_fp
  REAL(fp), PARAMETER :: SENSOR_AZI =  40.0_fp   ! relative azimuth = +60

  INTEGER , PARAMETER :: N_ANGLES = 1
  ! FASTEM6 returns the polarimetric components as an untouched ZERO, so the
  ! null side is exact rather than approximate.
  REAL(fp), PARAMETER :: TOL          = 0.0_fp
  REAL(fp), PARAMETER :: SIGNAL_FLOOR = 1.0e-8_fp

  CHARACTER(256) :: Version
  REAL(fp) :: eU_f4, eV_f4, eU_f6, eV_f6
  LOGICAL  :: ok_selects_f4, ok_selects_f6, all_ok

  CALL CRTM_Version(Version)
  WRITE(*,'(/5x,a)') 'MWwaterCoeff_File model-selection verification'
  WRITE(*,'(5x,a/)') 'CRTM Version: '//TRIM(Version)

  CALL surface_UV_for_file( FILE_F4, eU_f4, eV_f4 )
  CALL surface_UV_for_file( FILE_F6, eU_f6, eV_f6 )

  ! Asking for FASTEM4 must give the model that carries U and V ...
  ok_selects_f4 = ( ABS(eU_f4) > SIGNAL_FLOOR ) .AND. ( ABS(eV_f4) > SIGNAL_FLOOR )
  ! ... and asking for FASTEM6 must give the model that does not, so that the
  ! first assertion demonstrates selection rather than a new fixed default.
  ok_selects_f6 = ( ABS(eU_f6) <= TOL ) .AND. ( ABS(eV_f6) <= TOL )

  WRITE(*,'(5x,a,a)')              'MWwaterCoeff_File = ', FILE_F4
  WRITE(*,'(5x,a,es14.6,a,es14.6)')'    surface  U = ', eU_f4, '   V = ', eV_f4
  WRITE(*,'(5x,a,a)')              'MWwaterCoeff_File = ', FILE_F6
  WRITE(*,'(5x,a,es14.6,a,es14.6)')'    surface  U = ', eU_f6, '   V = ', eV_f6

  WRITE(*,'(/5x,a,l1)') 'FASTEM4 requested and delivered (U,V nonzero) ...  pass = ', ok_selects_f4
  WRITE(*,'(5x,a,l1)')  'FASTEM6 requested and delivered (U,V zero) ......  pass = ', ok_selects_f6

  all_ok = ok_selects_f4 .AND. ok_selects_f6

  WRITE(*,'(/5x,a)') '=================================================='
  IF ( all_ok ) THEN
    WRITE(*,'(5x,a)') 'RESULT: PASS - MWwaterCoeff_File selects the model'
  ELSE
    WRITE(*,'(5x,a)') 'RESULT: FAIL - MWwaterCoeff_File did not select the model'
  END IF
  WRITE(*,'(5x,a/)') '=================================================='

  IF ( all_ok ) THEN
    STOP 0
  ELSE
    STOP 1
  END IF

CONTAINS

  ! Initialise CRTM selecting the water model by filename alone, evaluate the
  ! ocean surface optics once, and return the polarimetric components.
  SUBROUTINE surface_UV_for_file( MWfile, eU, eV )
    CHARACTER(*), INTENT(IN)  :: MWfile
    REAL(fp)    , INTENT(OUT) :: eU, eV

    INTEGER                      :: err
    TYPE(CRTM_ChannelInfo_type)  :: ChannelInfo(1)
    TYPE(CRTM_Surface_type)      :: Sfc
    TYPE(CRTM_GeometryInfo_type) :: gInfo
    TYPE(CRTM_SfcOptics_type)    :: SfcOptics
    TYPE(iVar_type)              :: iVar

    ! Deliberately no MWwaterCoeff_Scheme: the filename is the only selector.
    err = CRTM_Init( SENSORS, ChannelInfo,             &
                     File_Path         = PATH,         &
                     MWwaterCoeff_File = MWfile,       &
                     Quiet             = .TRUE.        )
    IF ( err /= SUCCESS ) THEN
      CALL Display_Message( PROGRAM_NAME, 'CRTM_Init failed for '//MWfile, FAILURE ); STOP 1
    END IF

    Sfc%Water_Coverage    = ONE
    Sfc%Land_Coverage     = ZERO
    Sfc%Snow_Coverage     = ZERO
    Sfc%Ice_Coverage      = ZERO
    Sfc%Water_Type        = 1
    Sfc%Water_Temperature = WATER_TEMP
    Sfc%Wind_Speed        = WIND_SPEED
    Sfc%Wind_Direction    = WIND_DIR
    Sfc%Salinity          = SALINITY

    CALL CRTM_GeometryInfo_SetValue( gInfo, Sensor_Zenith_Angle  = ZENITH, &
                                            Sensor_Azimuth_Angle = SENSOR_AZI )
    CALL CRTM_GeometryInfo_Compute( gInfo )

    CALL CRTM_SfcOptics_Create( SfcOptics, N_ANGLES, MAX_N_STOKES )
    IF ( .NOT. CRTM_SfcOptics_Associated(SfcOptics) ) THEN
      CALL Display_Message( PROGRAM_NAME, 'SfcOptics_Create failed', FAILURE ); STOP 1
    END IF
    SfcOptics%Angle(1)      = ZENITH
    SfcOptics%Weight(1)     = ONE
    SfcOptics%Index_Sat_Ang = 1
    SfcOptics%n_Angles      = N_ANGLES
    ! Scalar branch: it writes component 1 only, so components 3 and 4 still
    ! hold exactly what the surface model produced.
    SfcOptics%n_Stokes      = 1

    err = CRTM_Compute_SfcOptics( Sfc, gInfo, 1, CHANNEL, SfcOptics, iVar )
    IF ( err /= SUCCESS ) THEN
      CALL Display_Message( PROGRAM_NAME, 'Compute_SfcOptics failed', FAILURE ); STOP 1
    END IF

    eU = SfcOptics%Emissivity(1,3)
    eV = SfcOptics%Emissivity(1,4)

    CALL CRTM_SfcOptics_Destroy( SfcOptics )
    err = CRTM_Destroy( ChannelInfo )

  END SUBROUTINE surface_UV_for_file

END PROGRAM test_MWwaterCoeff_FileSelects
