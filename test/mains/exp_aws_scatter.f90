!
! exp_aws_scatter
!
! Validation driver for the experimental ('CRTM-Exp') Ren-snow CloudCoeff LUT.
! Runs the MWR_AWS sensor over an ocean US-Standard profile and reports the
! frozen-cloud scattering depression (TB_clear - TB_cloudy) for ICE / SNOW /
! GRAUPEL across the AWS scattering channels (89, 165.5, 183, 325 GHz).
!
! Build:  see tools/cloudcoeff_exp/build_and_run_exp_aws.sh
! Run from the repo test/ directory (needs ./testinput/ coeffs + the LUT).
!
PROGRAM exp_aws_scatter

  USE CRTM_Module
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'exp_aws_scatter'
  CHARACTER(*), PARAMETER :: PATH    = './testinput/'
  CHARACTER(*), PARAMETER :: SENSOR  = 'mwr_aws'
  CHARACTER(*), PARAMETER :: LUT     = 'CloudCoeff_Exp_Full6.nc'

  INTEGER, PARAMETER :: N_PROFILES = 2          ! 1 = clear, 2 = cloudy
  INTEGER, PARAMETER :: N_LAYERS   = 100
  INTEGER, PARAMETER :: N_ABSORBERS= 6
  INTEGER, PARAMETER :: N_CLOUDS   = 1
  INTEGER, PARAMETER :: N_AEROSOLS = 0
  REAL(fp), PARAMETER :: ZENITH = 53.0_fp       ! AWS conical scan ~53 deg

  ! cloud vertical band (warm-ish band so liquid habits are realistic; frozen still valid)
  INTEGER, PARAMETER :: KC1 = 78, KC2 = 86

  TYPE(CRTM_ChannelInfo_type)             :: chinfo(1)
  TYPE(CRTM_Geometry_type)                :: geo(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:)

  INTEGER :: err, nch, m, jc
  REAL(fp) :: dep(4), tbclr
  ! AWS scattering channels of interest (ch 9=89, 10=165.5, 11=176.3, 16=325.2 GHz)
  INTEGER, PARAMETER :: NSHOW = 4
  INTEGER, PARAMETER :: SHOW_CH(NSHOW) = (/ 9, 10, 11, 16 /)

  ! experiment sweep — all 6 habits of the complete LUT
  INTEGER  :: hk, wk
  INTEGER,  PARAMETER :: NHAB = 6
  INTEGER,  PARAMETER :: HID(NHAB)  = (/ WATER_CLOUD, ICE_CLOUD, RAIN_CLOUD, &
                                         SNOW_CLOUD, GRAUPEL_CLOUD, HAIL_CLOUD /)
  CHARACTER(6), PARAMETER :: HNM(NHAB) = (/ 'WATER ','ICE   ','RAIN  ','SNOW  ','GRAUP ','HAIL  '/)
  REAL(fp), PARAMETER :: REFF(NHAB) = (/ 15.0_fp, 60.0_fp, 1000.0_fp, &
                                         500.0_fp, 500.0_fp, 4000.0_fp /)   ! microns -> Dm
  REAL(fp), PARAMETER :: WC(3)   = (/ 0.05_fp, 0.2_fp, 0.5_fp /)       ! kg/m^2 per layer

  ! --------------------------------------------------------------------------
  WRITE(*,'(/a)') '=== Initializing CRTM (mwr_aws, Cloud_Model=CRTM-Exp) ==='
  err = CRTM_Init( (/ SENSOR /), chinfo,            &
                   Cloud_Model       = 'CRTM-Exp',  &
                   CloudCoeff_File   = LUT,         &
                   CloudCoeff_Format = 'netCDF',    &
                   File_Path         = PATH,        &
                   Quiet             = .TRUE. )
  IF ( err /= SUCCESS ) THEN
    CALL Display_Message(PROGRAM_NAME,'CRTM_Init failed',FAILURE); STOP 1
  END IF
  nch = SUM(CRTM_ChannelInfo_n_Channels(chinfo))
  WRITE(*,'(a,i0,a)') ' loaded ', nch, ' AWS channels'

  ALLOCATE( rts(nch, N_PROFILES) )
  CALL CRTM_Atmosphere_Create( atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
    CALL Display_Message(PROGRAM_NAME,'Atmosphere_Create failed',FAILURE); STOP 1
  END IF

  CALL Load_ECMWF84_Atm_Data()        ! fills atm(1) and atm(2) gas profiles (US-Std)
  atm(2) = atm(1)                      ! make the two profiles identical (clean clear vs cloudy)

  ! clear profile: no cloud
  atm(1)%n_Clouds = 0
  atm(1)%Cloud_Fraction = ZERO

  ! ocean surface, both profiles
  DO m = 1, N_PROFILES
    sfc(m)%Water_Coverage    = 1.0_fp
    sfc(m)%Water_Type        = 1          ! SEA_WATER
    sfc(m)%Water_Temperature = 290.0_fp
    sfc(m)%Wind_Speed        = 6.0_fp
    sfc(m)%Salinity          = 33.0_fp
    CALL CRTM_Geometry_SetValue( geo(m), Sensor_Zenith_Angle = ZENITH )
  END DO

  WRITE(*,'(a,i0,a,i0,a,2(f6.1),a)') ' cloud band layers ', KC1, '-', KC2, &
       '  (T = ', atm(1)%Temperature(KC1), atm(1)%Temperature(KC2), ' K)'

  ! --------------------------------------------------------------------------
  WRITE(*,'(/a)') '=== Frozen-cloud scattering depression  TB_clear - TB_cloudy  (K) ==='
  WRITE(*,'(a)')  '  habit   Reff   WC/lyr |   89.0   165.5   176.3   325.2 GHz   | TBclr@89'
  WRITE(*,'(a)')  '  -----   ----   ------ | ------  ------  ------  ------       | --------'

  DO hk = 1, NHAB
    DO wk = 1, 3
      ! reset cloudy profile = clear profile + a frozen cloud in the band
      atm(2)            = atm(1)
      atm(2)%n_Clouds   = 1
      atm(2)%Cloud_Fraction = ZERO
      atm(2)%Cloud_Fraction(KC1:KC2)          = 1.0_fp
      atm(2)%Cloud(1)%Type                    = HID(hk)
      atm(2)%Cloud(1)%Effective_Radius        = ZERO
      atm(2)%Cloud(1)%Water_Content           = ZERO
      atm(2)%Cloud(1)%Effective_Radius(KC1:KC2) = REFF(hk)
      atm(2)%Cloud(1)%Water_Content(KC1:KC2)    = WC(wk)

      err = CRTM_Forward( atm, sfc, geo, chinfo, rts )
      IF ( err /= SUCCESS ) THEN
        CALL Display_Message(PROGRAM_NAME,'CRTM_Forward failed',FAILURE); STOP 1
      END IF

      DO jc = 1, NSHOW
        dep(jc) = rts(SHOW_CH(jc),1)%Brightness_Temperature &
                - rts(SHOW_CH(jc),2)%Brightness_Temperature
      END DO
      tbclr = rts(SHOW_CH(1),1)%Brightness_Temperature
      WRITE(*,'(2x,a6,f6.0,f8.3," |",4f8.2,"       |",f8.2)') &
        HNM(hk), REFF(hk), WC(wk), dep, tbclr
    END DO
  END DO

  WRITE(*,'(/a)') '=== done ==='
  DEALLOCATE( rts )
  CALL CRTM_Atmosphere_Destroy( atm )
  err = CRTM_Destroy( chinfo )

CONTAINS

  INCLUDE 'Load_ECMWF84_Atm_Data.inc'

END PROGRAM exp_aws_scatter
