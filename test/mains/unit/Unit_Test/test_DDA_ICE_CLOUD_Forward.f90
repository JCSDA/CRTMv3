!
! test_DDA_ICE_CLOUD_Forward
!
! Regression coverage for the v3.2.0 DDA-ARTS behavior change: ICE_CLOUD now goes
! through the full scattering branch (the legacy non-scattering shortcut applies
! only to Mie-TAMU tables), with the default DDA habit IconCloudIce. This altered
! radiances/Jacobians for DDA-ARTS users, and nothing pinned the new values.
!
! Loads a DDA-ARTS CloudCoeff (CloudCoeff_DDA_Moradi_2022.nc) and runs atms_n21
! (183 GHz channels) over an ocean US-Standard column for three profiles:
!     1 = clear, 2 = thin ice, 3 = heavy ice
! and asserts the ICE_CLOUD scattering path is physical and ACTIVE:
!     * every TB (all channels/profiles) is physical,
!     * an ice cloud produces a brightness-temperature DEPRESSION (scattering),
!     * the depression GROWS with ice water content.
! If ICE_CLOUD ever reverts to the non-scattering shortcut under DDA-ARTS, the
! depression collapses and this test fails.
!
! STOP 0 on success, STOP 1 on failure.
!
PROGRAM test_DDA_ICE_CLOUD_Forward

  USE CRTM_Module
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_DDA_ICE_CLOUD_Forward'
  CHARACTER(*), PARAMETER :: PATH   = './testinput/'
  CHARACTER(*), PARAMETER :: SENSOR = 'atms_n21'
  CHARACTER(*), PARAMETER :: LUT    = 'CloudCoeff_DDA_Moradi_2022.nc'

  INTEGER,  PARAMETER :: N_PROFILES  = 3       ! 1=clear, 2=thin ice, 3=heavy ice
  INTEGER,  PARAMETER :: N_LAYERS    = 100
  INTEGER,  PARAMETER :: N_ABSORBERS = 6
  INTEGER,  PARAMETER :: N_CLOUDS    = 1
  INTEGER,  PARAMETER :: N_AEROSOLS  = 0
  REAL(fp), PARAMETER :: ZENITH      = 30.0_fp
  INTEGER,  PARAMETER :: KC1 = 60, KC2 = 72     ! upper-tropospheric (cold) ice band
  REAL(fp), PARAMETER :: REFF_I      = 100.0_fp ! ice effective radius (microns)
  REAL(fp), PARAMETER :: WC_THIN     = 0.05_fp  ! kg/m^2 per layer (thin)
  REAL(fp), PARAMETER :: WC_HEAVY    = 1.00_fp  ! kg/m^2 per layer (heavy)

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

  ! Init CRTM with the DDA-ARTS cloud table (scheme is read from the file)
  err = CRTM_Init( (/ SENSOR /), chinfo,         &
                   CloudCoeff_File   = LUT,       &
                   CloudCoeff_Format = 'netCDF',  &
                   File_Path         = PATH,      &
                   Quiet             = .TRUE. )
  IF ( err /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Init (DDA-ARTS CloudCoeff) failed', FAILURE ); STOP 1
  END IF
  nch = SUM( CRTM_ChannelInfo_n_Channels(chinfo) )
  IF ( nch < 1 ) THEN
    CALL Display_Message( PROGRAM_NAME, 'no channels loaded for '//SENSOR, FAILURE ); STOP 1
  END IF

  ALLOCATE( rts(nch, N_PROFILES) )
  CALL CRTM_Atmosphere_Create( atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Atmosphere_Create failed', FAILURE ); STOP 1
  END IF

  CALL Load_ECMWF84_Atm_Data()
  DO m = 2, N_PROFILES
    atm(m) = atm(1)
  END DO
  atm(1)%n_Clouds       = 0
  atm(1)%Cloud_Fraction = ZERO
  CALL Set_Ice( atm(2), WC_THIN )
  CALL Set_Ice( atm(3), WC_HEAVY )

  DO m = 1, N_PROFILES
    sfc(m)%Water_Coverage    = 1.0_fp
    sfc(m)%Water_Type        = 1          ! SEA_WATER
    sfc(m)%Water_Temperature = 290.0_fp
    sfc(m)%Wind_Speed        = 6.0_fp
    sfc(m)%Salinity          = 33.0_fp
    CALL CRTM_Geometry_SetValue( geo(m), Sensor_Zenith_Angle = ZENITH )
  END DO

  err = CRTM_Forward( atm, sfc, geo, chinfo, rts )
  IF ( err /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Forward failed', FAILURE ); STOP 1
  END IF

  tb_min = MINVAL( rts%Brightness_Temperature )
  tb_max = MAXVAL( rts%Brightness_Temperature )
  dep_thin  = MAXVAL( rts(:,1)%Brightness_Temperature - rts(:,2)%Brightness_Temperature )
  dep_heavy = MAXVAL( rts(:,1)%Brightness_Temperature - rts(:,3)%Brightness_Temperature )

  WRITE(*,'(/a)')       ' DDA-ARTS ICE_CLOUD forward scattering check (atms_n21, ice over ocean):'
  WRITE(*,'(a,i0)')     '   channels                 : ', nch
  WRITE(*,'(a,2f8.2)')  '   all-profile TB range (K) : ', tb_min, tb_max
  WRITE(*,'(a,f8.3)')   '   max depression, thin     : ', dep_thin
  WRITE(*,'(a,f8.3)')   '   max depression, heavy    : ', dep_heavy

  IF ( tb_min < TB_LO .OR. tb_max > TB_HI ) THEN
    WRITE(*,'(a,2f10.2)') ' FAIL: a TB is outside the physical range ', tb_min, tb_max
    ok = .FALSE.
  END IF
  IF ( dep_heavy < MIN_DEPRESSION ) THEN
    WRITE(*,'(a,f8.3,a,f6.2,a)') ' FAIL: heavy-ice depression ', dep_heavy, &
         ' K is below the ', MIN_DEPRESSION, ' K threshold (ICE_CLOUD not scattering?)'
    ok = .FALSE.
  END IF
  IF ( dep_heavy <= dep_thin ) THEN
    WRITE(*,'(a)') ' FAIL: depression did not increase with ice water content'
    ok = .FALSE.
  END IF

  DEALLOCATE( rts )
  CALL CRTM_Atmosphere_Destroy( atm )
  err = CRTM_Destroy( chinfo )

  IF ( ok ) THEN
    WRITE(*,'(/a)') ' PASS: DDA-ARTS ICE_CLOUD forward path produces physical scattering.'
    STOP 0
  ELSE
    WRITE(*,'(/a)') ' FAIL: DDA-ARTS ICE_CLOUD forward checks failed.'
    STOP 1
  END IF

CONTAINS

  ! Put an ice cloud of the given per-layer water content into the band.
  SUBROUTINE Set_Ice( a, wc )
    TYPE(CRTM_Atmosphere_type), INTENT(IN OUT) :: a
    REAL(fp),                   INTENT(IN)     :: wc
    a%n_Clouds                            = 1
    a%Cloud_Fraction                      = ZERO
    a%Cloud_Fraction(KC1:KC2)             = 1.0_fp
    a%Cloud(1)%Type                       = ICE_CLOUD
    a%Cloud(1)%Effective_Radius           = ZERO
    a%Cloud(1)%Water_Content              = ZERO
    a%Cloud(1)%Effective_Radius(KC1:KC2)  = REFF_I
    a%Cloud(1)%Water_Content(KC1:KC2)     = wc
  END SUBROUTINE Set_Ice

  INCLUDE 'Load_ECMWF84_Atm_Data.inc'

END PROGRAM test_DDA_ICE_CLOUD_Forward
