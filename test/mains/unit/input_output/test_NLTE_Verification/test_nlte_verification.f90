!-------------------------------------------------------
!
! Description:
!       Verification test for the netCDF SpcCoeff sibling-file
!       load introduced in CRTM v3.1.5 (backport of the
!       REL-3.2.0 fix).
!
!       The netCDF SpcCoeff format stores the non-LTE correction
!       and antenna correction coefficients in separate
!       <sensor>.NLTECoeff.nc / <sensor>.ACCoeff.nc sibling files
!       (the binary format embeds them inline).  Prior to v3.1.5
!       the netCDF reader never loaded the siblings, silently
!       disabling both corrections.
!
!       This test loads iasi616_metop-b (IR, 27 NLTE channels) and
!       amsua_n19 (MW, antenna correction) from two staged
!       directories -- ./testinput/ with the siblings co-located,
!       ./testinput_no_siblings/ without them -- and asserts:
!
!         1. With siblings present, NLTECoeff/ACCoeff are loaded.
!         2. Without them, the load still succeeds (graceful
!            degradation, matching pre-fix behaviour).
!         3. Forward model, sibling on vs off: daytime BTs differ
!            by > 0.5 K on NLTE channels ONLY; night-time BTs are
!            identical (the NLTE correction is solar-driven).
!
!       Assertion 1 and the daytime part of assertion 3 fail on an
!       unfixed 3.1.x library, guarding against regression of the
!       sibling load.
!
!       Date: 2026-08-23       Author: J. Benjamin
!
!-------------------------------------------------------

PROGRAM test_nlte_verification

  ! Module usage
  USE CRTM_Module
  USE CRTM_SpcCoeff,    ONLY: SC, CRTM_SpcCoeff_Load, CRTM_SpcCoeff_Destroy
  USE NLTECoeff_Define, ONLY: NLTECoeff_Associated
  USE ACCoeff_Define,   ONLY: ACCoeff_Associated
  USE UnitTest_Define,  ONLY: UnitTest_type,     &
                              UnitTest_Init,     &
                              UnitTest_Setup,    &
                              UnitTest_Assert,   &
                              UnitTest_Report,   &
                              UnitTest_n_Failed

  ! Disable all implicit typing
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: Program_Name = 'test_nlte_verification'
  CHARACTER(*), PARAMETER :: IR_SENSOR    = 'iasi616_metop-b'
  CHARACTER(*), PARAMETER :: MW_SENSOR    = 'amsua_n19'
  CHARACTER(*), PARAMETER :: PATH_SIB     = './testinput/'
  CHARACTER(*), PARAMETER :: PATH_NOSIB   = './testinput_no_siblings/'
  INTEGER,      PARAMETER :: N_IASI616_NLTE_CHANNELS = 27

  ! Profile and geometry setup (matches the forward regression tests)
  INTEGER,  PARAMETER :: N_PROFILES  = 2
  INTEGER,  PARAMETER :: N_LAYERS    = 92
  INTEGER,  PARAMETER :: N_ABSORBERS = 2
  INTEGER,  PARAMETER :: N_CLOUDS    = 1
  INTEGER,  PARAMETER :: N_AEROSOLS  = 1
  REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp
  REAL(fp), PARAMETER :: SZA_DAY   =  45.0_fp
  REAL(fp), PARAMETER :: SZA_NIGHT = 120.0_fp

  ! With these profiles the NLTE band-head channels shift by ~4 K at
  ! SZA=45; anything above this floor proves the correction is active.
  REAL(fp), PARAMETER :: NLTE_MIN_DAY_DBT = 0.5_fp
  ! Sibling on/off must be bit-comparable wherever NLTE is inactive.
  REAL(fp), PARAMETER :: TOL_IDENTICAL = 1.0e-9_fp

  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(1)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)

  TYPE(UnitTest_type) :: utest
  INTEGER  :: err_stat, n_Channels, l
  LOGICAL  :: nlte_loaded
  LOGICAL,  ALLOCATABLE :: nlte_mask(:)
  REAL(fp), ALLOCATABLE :: bt_on_day(:,:),  bt_on_night(:,:)
  REAL(fp), ALLOCATABLE :: bt_off_day(:,:), bt_off_night(:,:)
  REAL(fp) :: max_night, max_day_nlte, max_day_other

  CALL UnitTest_Init(utest, .TRUE.)
  CALL UnitTest_Setup(utest, 'NLTE_Verification_Test', Program_Name, .TRUE.)


  ! ============================================================
  ! Part 1: read-level checks via CRTM_SpcCoeff_Load
  ! ============================================================

  ! IR sensor, sibling present: NLTECoeff must be loaded
  err_stat = CRTM_SpcCoeff_Load( (/IR_SENSOR/), &
                                 File_Path = PATH_SIB, &
                                 netCDF    = .TRUE.,   &
                                 Quiet     = .TRUE. )
  WRITE(*,'(a)') 'ASSERT: '//IR_SENSOR//' netCDF load (sibling present) succeeds'
  CALL UnitTest_Assert(utest, err_stat == SUCCESS)
  WRITE(*,'(a)') 'ASSERT: NLTECoeff loaded from co-located sibling file'
  CALL UnitTest_Assert(utest, NLTECoeff_Associated(SC(1)%NC))
  IF ( NLTECoeff_Associated(SC(1)%NC) ) THEN
    WRITE(*,'(a,i0)') '        n_NLTE_Channels = ', SC(1)%NC%n_NLTE_Channels
    CALL UnitTest_Assert(utest, SC(1)%NC%n_NLTE_Channels == N_IASI616_NLTE_CHANNELS)
  ELSE
    CALL UnitTest_Assert(utest, .FALSE.)
  END IF
  err_stat = CRTM_SpcCoeff_Destroy()

  ! IR sensor, sibling absent: load succeeds, NLTECoeff stays empty
  err_stat = CRTM_SpcCoeff_Load( (/IR_SENSOR/), &
                                 File_Path = PATH_NOSIB, &
                                 netCDF    = .TRUE.,     &
                                 Quiet     = .TRUE. )
  WRITE(*,'(a)') 'ASSERT: '//IR_SENSOR//' netCDF load (sibling absent) still succeeds'
  CALL UnitTest_Assert(utest, err_stat == SUCCESS)
  WRITE(*,'(a)') 'ASSERT: NLTECoeff not associated when sibling absent'
  CALL UnitTest_Assert(utest, .NOT. NLTECoeff_Associated(SC(1)%NC))
  err_stat = CRTM_SpcCoeff_Destroy()

  ! MW sensor: same contract for the ACCoeff sibling
  err_stat = CRTM_SpcCoeff_Load( (/MW_SENSOR/), &
                                 File_Path = PATH_SIB, &
                                 netCDF    = .TRUE.,   &
                                 Quiet     = .TRUE. )
  WRITE(*,'(a)') 'ASSERT: '//MW_SENSOR//' ACCoeff loaded from co-located sibling file'
  CALL UnitTest_Assert(utest, err_stat == SUCCESS .AND. ACCoeff_Associated(SC(1)%AC))
  err_stat = CRTM_SpcCoeff_Destroy()

  err_stat = CRTM_SpcCoeff_Load( (/MW_SENSOR/), &
                                 File_Path = PATH_NOSIB, &
                                 netCDF    = .TRUE.,     &
                                 Quiet     = .TRUE. )
  WRITE(*,'(a)') 'ASSERT: '//MW_SENSOR//' ACCoeff not associated when sibling absent'
  CALL UnitTest_Assert(utest, err_stat == SUCCESS .AND. .NOT. ACCoeff_Associated(SC(1)%AC))
  err_stat = CRTM_SpcCoeff_Destroy()


  ! ============================================================
  ! Part 2: forward-model check that the NLTE correction is
  !         active by day and inert by night
  ! ============================================================

  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
    CALL Display_Message( Program_Name, 'Error creating Atmosphere', FAILURE )
    STOP 1
  END IF
  CALL Load_Atm_Data()
  CALL Load_Sfc_Data()

  CALL Run_Forward( PATH_SIB,   bt_on_day,  bt_on_night,  nlte_loaded )
  WRITE(*,'(a)') 'ASSERT: NLTECoeff associated after CRTM_Init (sibling present)'
  CALL UnitTest_Assert(utest, nlte_loaded)

  CALL Run_Forward( PATH_NOSIB, bt_off_day, bt_off_night, nlte_loaded )
  WRITE(*,'(a)') 'ASSERT: NLTECoeff not associated after CRTM_Init (sibling absent)'
  CALL UnitTest_Assert(utest, .NOT. nlte_loaded)

  IF ( .NOT. ALLOCATED(nlte_mask) ) THEN
    CALL Display_Message( Program_Name, &
      'NLTE channel mask never built -- sibling load is broken; '// &
      'skipping BT comparison', FAILURE )
    CALL UnitTest_Assert(utest, .FALSE.)
  ELSE
    max_night     = MAXVAL( ABS(bt_on_night - bt_off_night) )
    max_day_nlte  = MAXVAL( ABS(bt_on_day - bt_off_day), &
                            MASK=SPREAD(nlte_mask, DIM=2, NCOPIES=N_PROFILES) )
    max_day_other = MAXVAL( ABS(bt_on_day - bt_off_day), &
                            MASK=SPREAD(.NOT. nlte_mask, DIM=2, NCOPIES=N_PROFILES) )

    WRITE(*,'(a,es13.6,a)') 'Sibling on/off max |dBT|, night, all channels  : ', max_night, ' K'
    WRITE(*,'(a,es13.6,a)') 'Sibling on/off max |dBT|, day, NLTE channels   : ', max_day_nlte, ' K'
    WRITE(*,'(a,es13.6,a)') 'Sibling on/off max |dBT|, day, other channels  : ', max_day_other, ' K'

    WRITE(*,'(a)') 'ASSERT: night-time BTs identical with and without sibling'
    CALL UnitTest_Assert(utest, max_night < TOL_IDENTICAL)
    WRITE(*,'(a)') 'ASSERT: daytime NLTE-channel BTs shifted by the correction'
    CALL UnitTest_Assert(utest, max_day_nlte > NLTE_MIN_DAY_DBT)
    WRITE(*,'(a)') 'ASSERT: daytime non-NLTE-channel BTs unaffected'
    CALL UnitTest_Assert(utest, max_day_other < TOL_IDENTICAL)
  END IF

  CALL CRTM_Atmosphere_Destroy(Atm)

  CALL UnitTest_Report(utest)
  IF ( UnitTest_n_Failed(utest) > 0 ) STOP 1


CONTAINS


  ! Initialize the CRTM for the IR sensor from Path with netCDF
  ! SpcCoeff, report whether the NLTE coefficients were picked up,
  ! and return day/night brightness temperatures.
  SUBROUTINE Run_Forward( Path, BT_Day, BT_Night, NLTE_Loaded )
    CHARACTER(*),          INTENT(IN)  :: Path
    REAL(fp), ALLOCATABLE, INTENT(OUT) :: BT_Day(:,:), BT_Night(:,:)
    LOGICAL,               INTENT(OUT) :: NLTE_Loaded
    INTEGER :: err, alloc_stat

    err = CRTM_Init( (/IR_SENSOR/), &
                     ChannelInfo, &
                     File_Path       = Path, &
                     SpcCoeff_Format = 'netCDF', &
                     Quiet           = .TRUE. )
    IF ( err /= SUCCESS ) THEN
      CALL Display_Message( Program_Name, &
        'Error initializing CRTM from '//Path, FAILURE )
      STOP 1
    END IF
    n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))

    NLTE_Loaded = NLTECoeff_Associated( SC(1)%NC )
    IF ( NLTE_Loaded .AND. .NOT. ALLOCATED(nlte_mask) ) THEN
      ALLOCATE( nlte_mask(n_Channels) )
      DO l = 1, n_Channels
        nlte_mask(l) = ANY( SC(1)%NC%NLTE_Channel(:) == SC(1)%Sensor_Channel(l) )
      END DO
      WRITE(*,'(a,i0,a)') 'NLTE channel mask: ', COUNT(nlte_mask), ' channels flagged'
    END IF

    ALLOCATE( BT_Day(n_Channels, N_PROFILES), &
              BT_Night(n_Channels, N_PROFILES), &
              RTSolution(n_Channels, N_PROFILES), STAT=alloc_stat )
    IF ( alloc_stat /= 0 ) STOP 1
    CALL CRTM_RTSolution_Create( RTSolution, N_LAYERS )

    CALL Run_One_SZA( SZA_DAY,   BT_Day )
    CALL Run_One_SZA( SZA_NIGHT, BT_Night )

    DEALLOCATE( RTSolution )
    err = CRTM_Destroy( ChannelInfo )
    IF ( err /= SUCCESS ) THEN
      CALL Display_Message( Program_Name, 'Error destroying CRTM', FAILURE )
      STOP 1
    END IF
  END SUBROUTINE Run_Forward


  SUBROUTINE Run_One_SZA( SZA, BT )
    REAL(fp), INTENT(IN)  :: SZA
    REAL(fp), INTENT(OUT) :: BT(:,:)
    INTEGER :: err, m

    CALL CRTM_Geometry_SetValue( Geometry, &
                                 Sensor_Zenith_Angle = ZENITH_ANGLE, &
                                 Sensor_Scan_Angle   = SCAN_ANGLE, &
                                 Source_Zenith_Angle = SZA )
    err = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution )
    IF ( err /= SUCCESS ) THEN
      CALL Display_Message( Program_Name, 'Error in CRTM Forward Model', FAILURE )
      STOP 1
    END IF
    DO m = 1, N_PROFILES
      DO l = 1, n_Channels
        BT(l,m) = RTSolution(l,m)%Brightness_Temperature
      END DO
    END DO
  END SUBROUTINE Run_One_SZA


  INCLUDE 'Load_Atm_Data.inc'
  INCLUDE 'Load_Sfc_Data.inc'

END PROGRAM test_nlte_verification
