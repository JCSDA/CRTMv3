!
! Test_CRTM_V30
!
! This code may be used for testing the CRTM code as a bench mark and the consistency.
! Each run, the code performs the test for one sensor and one algorithm.  
!
! One hundred profiles over various surface types and different Sun/sensor zenith and azimuth
!  angles are used. The Sun zenith angles are less than 70 degrees here. Both clear and cloudy
!  with one dust aerosol are used. One runs this by using any number of profiles by using the
!  parameters "n1" and "n2". 
!
!  The CRTM V3.0 can be used for microwave, infrared sensors for n_Stokes =1. Same as previous versions,
!    only Fourier zeroth component for azimuth angle is considered. The surface emissivity/
!    reflectivity is NOT for the zeroth component and it may be "actual" value that depends
!    on incident/outgoing directions.
!
!  For visible/ultraviolet sensors, the CRTM 3.0 can be performed in either scalar model for n_Stokes = 1
!  (using RT_ADA) or fully polarized mode for n_Stokes = 4 (using RT_VMOM).
!
!
! CREATION HISTORY:
!       Written by:     David Groff, 9-Aug-2010
!                       david.groff@noaa.gov
!       Updated by:     Quanhua (Mark) Liu, 12-Jan-2022
!                       Quanhua.Liu@noaa.gov

PROGRAM Test_CRTM_V30

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module usage
  USE Type_Kinds
  USE Message_Handler
  USE Sort_Utility
  USE CRTM_Module
  USE UnitTest_Define
  USE Timing_Utility
  USE CRTM_CloudCover_Define, ONLY: DEFAULT_OVERLAP_ID, &
                                    CloudCover_Maximum_Overlap, &
                                    CloudCover_Random_Overlap , &
                                    CloudCover_MaxRan_Overlap , &
                                    CloudCover_Average_Overlap, &
                                    CloudCover_Overcast_Overlap,&
                                    CloudCover_Overlap_IsValid, &
                                    CloudCover_Overlap_Name
  USE CRTM_SpcCoeff,          ONLY: SC, &
                                    SpcCoeff_IsInfraredSensor, &
                                    SpcCoeff_IsMicrowaveSensor
  ! Disable all implicit typing
  IMPLICIT NONE
!   INTEGER :: 
  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*),  PARAMETER :: PROGRAM_NAME = 'Test_CRTM_V30'
  CHARACTER(*),  PARAMETER :: PROGRAM_VERSION_ID = &
    '$Id: Test_CRTM_V30.f90  $'

!    ----------------  Setting ------------------------------------------------------
  INTEGER , PARAMETER :: Test_Case = 1  ! 1: MW atms_n20, 2: IR modis_aqua, !
           ! 3: VIS v.viirs-m_j1 (scalar solver), 4: VIS v.viirs-m_j1 (vector solver)
  LOGICAL :: FWD_check = .true.  
  LOGICAL :: TL_check = .true.  
  LOGICAL :: AD_check = .true.
  
  LOGICAL :: FWD_TL = .false.
  LOGICAL :: TL_AD = .false.
  LOGICAL :: AD_KM = .false.  
  INTEGER, PARAMETER :: n1 = 1, n2 = 100, n_profile_step = 1 !start and end profile
  INTEGER, PARAMETER :: istoke=1  ! The first STokes component
  INTEGER, PARAMETER :: EMISSION_INDEX = 3
!    ----------------  Setting  END-----------------------------------------------
  LOGICAL :: IR_MW_Sensor = .false.
  INTEGER :: n_Stokes
  INTEGER     , PARAMETER :: N_RT_ALGORITHMS = 1
  CHARACTER(len=8) :: RT_ALGORITHM_NAME(N_RT_ALGORITHMS)
  INTEGER :: RT_ALGORITHM_ID(N_RT_ALGORITHMS)

  INTEGER     , PARAMETER :: TEST_N_SENSORS = 1
  CHARACTER(len=100) :: SENSOR_ID = 'v.modis_aqua'

!  INTEGER, PARAMETER :: used_n_channel = 9
!  INTEGER :: used_channel_index(used_n_channel)

  ! The test inputs
  CHARACTER(*), PARAMETER :: ATMOSPHERE_FILENAME = './CRTM_atm100.bin'
  CHARACTER(*), PARAMETER :: SURFACE_FILENAME    = './CRTM_sfc100.bin'
  CHARACTER(*), PARAMETER :: GEOMETRY_FILENAME   = './CRTM_geo100.bin'
  ! Tolerances
  REAL(fp), PARAMETER :: FWD_ULP   = 5.0e+02_fp  ! 2.5 ULPs
  REAL(fp), PARAMETER :: TL_ULP    = 5.0e+02_fp  ! 2.5 ULPs
  REAL(fp), PARAMETER :: TLAD_ULP  = 1.0e+06_fp
  INTEGER,  PARAMETER :: AD_SIGFIG = 6
  REAL(fp), PARAMETER :: DX = 0.001_fp
  ! Number of NL perturbations
  INTEGER , PARAMETER :: N_ALPHA = 2
  ! Scaling factors for NL perturbations
  REAL(fp), PARAMETER :: ALPHA(N_ALPHA) = [ 0.1000_fp, &
                                            0.0001_fp  ]

  REAL(fp), PARAMETER :: FWDTL_TOLERANCE(N_ALPHA) = [ 2.0e-04_fp, &  
                                                      2.0e-7_fp  ] 
  ! Generic
  LOGICAL , PARAMETER :: QUIET = .TRUE.

  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: err_msg, alloc_msg
  INTEGER :: err_stat, alloc_stat
  INTEGER :: i, ic, n, k, m 
  INTEGER :: n_aerosols, n_clouds, n_absorbers, n_layers
  INTEGER :: n_sensors, n_Channels
  INTEGER :: n_profiles
  INTEGER :: isensor1, isensor2, dsensor
  TYPE(UnitTest_type)  :: utest
  TYPE(Timing_type) :: Timing
  TYPE(CRTM_Atmosphere_type) , ALLOCATABLE :: atm(:)
  TYPE(CRTM_Surface_type)    , ALLOCATABLE :: sfc(:)
  TYPE(CRTM_Geometry_type)   , ALLOCATABLE :: geo(:)
  TYPE(CRTM_Options_type),    ALLOCATABLE :: opt(:)
  TYPE(CRTM_ChannelInfo_type) :: chinfo(TEST_N_SENSORS)

  ! Output header
  CALL Program_Message( PROGRAM_NAME, &
                        'A comprehensive CRTM test. The test is intended to account '//&
                        'for the full range of input conditions and instrument types.', &
                        '$Revision: 92324 $' )
 IF( Test_Case < 4 ) THEN
   n_Stokes = 1
   IF( Test_Case == 1 ) Sensor_Id = 'atms_n20'
   IF( Test_Case == 2 ) Sensor_Id = 'modis_aqua'   
   IF( Test_Case == 3 ) Sensor_Id = 'v.viirs-m_j1'
!  err_stat = CRTM_Init( (/Sensor_Id/), &  ! Input... must be an array, hence the (/../)
!         chinfo  , &  ! Output
!         File_Path='../crtm_coeff_2.4.1/')
  err_stat = CRTM_Init( (/Sensor_Id/), &  ! Input... must be an array, hence the (/../)
         chinfo  , &  ! Output
         CloudCoeff_File='Cloud_V3.bin', &
         AerosolCoeff_File='Aerosol_V3.bin', & 
         File_Path='../crtm_coeff_3.0_test/') 


     RT_ALGORITHM_NAME(N_RT_ALGORITHMS) = 'ADA     '
     RT_ALGORITHM_ID(N_RT_ALGORITHMS) = RT_ADA

 ELSE IF( Test_Case == 4 ) THEN
   n_Stokes = 4
   Sensor_Id = 'v.viirs-m_j1'
   err_stat = CRTM_Init( (/Sensor_Id/), &  ! Input... must be an array, hence the (/../)
         chinfo  , &  ! Output
         CloudCoeff_File='Cloud_V3.bin', &
         AerosolCoeff_File='Aerosol_V3.bin', & 
         File_Path='../crtm_coeff_3.0_test/') 
     RT_ALGORITHM_NAME(N_RT_ALGORITHMS) = 'VMOM    '
     RT_ALGORITHM_ID(N_RT_ALGORITHMS) = RT_VMOM
 END IF

! 
   IF( SpcCoeff_IsInfraredSensor(SC(1)).or.SpcCoeff_IsMicrowaveSensor(SC(1)) ) THEN
     IR_MW_Sensor = .true.
   END IF
 
  IF ( err_stat /= SUCCESS ) THEN
    err_msg = 'Error initialising the CRTM'
    CALL Display_Message( PROGRAM_NAME, err_msg, err_stat ); STOP
  END IF
 
  ! Allocate the structure arrays
  ! ...Determine the number of profiles
  err_stat = CRTM_Atmosphere_InquireFile( ATMOSPHERE_FILENAME, n_Profiles=n_profiles )
  IF ( err_stat /= SUCCESS ) THEN
    err_msg = 'Error inquiring Atmosphere datafile '//ATMOSPHERE_FILENAME
    CALL Display_Message( PROGRAM_NAME, err_msg, err_stat ); STOP
  END IF
  ! ...Determine the number of sensors
  n_sensors = TEST_N_SENSORS
  ! ...Allocate input data structures
  ALLOCATE( atm(n_profiles), sfc(n_profiles), geo(n_profiles), &
            STAT=alloc_stat, ERRMSG=alloc_msg )
  IF ( alloc_stat /= 0 ) THEN
    err_msg = 'Error allocating structure arrays - '//TRIM(alloc_msg)
    CALL Display_Message( PROGRAM_NAME, err_msg, err_stat ); STOP
  END IF

  ! Read the input datafiles
  ! ...Atmosphere
  print *,' n_Profiles ',n_Profiles
  err_stat = CRTM_Atmosphere_ReadFile( ATMOSPHERE_FILENAME, atm, Quiet=QUIET )

  IF ( err_stat /= SUCCESS ) THEN
    err_msg = 'Error reading Atmosphere datafile '//ATMOSPHERE_FILENAME
    CALL Display_Message( PROGRAM_NAME, err_msg, err_stat ); STOP
  END IF
  ! ...Surface
  err_stat = CRTM_Surface_ReadFile( SURFACE_FILENAME, sfc, Quiet=QUIET )
  IF ( err_stat /= SUCCESS ) THEN
    err_msg = 'Error reading Surface datafile '//SURFACE_FILENAME
    CALL Display_Message( PROGRAM_NAME, err_msg, err_stat ); STOP
  END IF
  ! ...Geometry
  err_stat = CRTM_Geometry_ReadFile( GEOMETRY_FILENAME, geo, Quiet=QUIET )
  
  IF ( err_stat /= SUCCESS ) THEN
    err_msg = 'Error reading Geometry datafile '//GEOMETRY_FILENAME
    CALL Display_Message( PROGRAM_NAME, err_msg, err_stat ); STOP
  END IF
  ! Unit test initialization
  CALL UnitTest_Init(utest)

  ALLOCATE( Opt(n_profiles) )
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(chinfo(:)))
  CALL CRTM_Options_Create( Opt, n_Channels )

!  Opt(1)%RT_Algorithm_Id = RT_VMOM !RT_VMOM  !RT_ADA
  DO k = 1, n_profiles
  Opt(:)%Include_Scattering = .TRUE.
  Opt(:)%Use_n_Streams = .TRUE.
  Opt(:)%n_Streams = 8
  Opt(:)%n_Stokes = n_Stokes
  END DO

 IF(FWD_check) CALL Test_CRTM_FWD(utest, atm(n1:n2), sfc(n1:n2), geo(n1:n2), chinfo)
 IF(TL_check) CALL Test_CRTM_TL(utest, atm(::n_profile_step), sfc(::n_profile_step), geo(::n_profile_step), chinfo)
 IF(AD_check) CALL Test_CRTM_AD(utest, atm(::n_profile_step), sfc(::n_profile_step), geo(::n_profile_step), chinfo)

  IF( FWD_TL ) THEN
    CALL Test_CRTM_FWDTL(utest, atm(n1:n2), sfc(n1:n2), geo(n1:n2), chinfo)
  END IF
  IF( TL_AD ) THEN  
    CALL Test_CRTM_TLAD(utest, atm(n1:n2), sfc(n1:n2), geo(n1:n2), chinfo) 
  END IF
  IF( AD_KM ) THEN
    CALL Test_CRTM_ADK(utest, atm(n1:n2), sfc(n1:n2), geo(n1:n2), chinfo)
  END IF

  ! Test Summary
  CALL UnitTest_Summary(utest)

  ! Cleanup
  err_stat = CRTM_Destroy(chinfo)
  DEALLOCATE( atm, sfc, geo )

CONTAINS

  ! ==================
  ! Forward model test
  ! ==================
  SUBROUTINE Test_CRTM_FWD(utest, atm, sfc, geo, chinfo)
    ! Arguments
    TYPE(UnitTest_type)        , INTENT(IN OUT) :: utest
    TYPE(CRTM_Atmosphere_type) , INTENT(IN)     :: atm(:)
    TYPE(CRTM_Surface_type)    , INTENT(IN)     :: sfc(:)
    TYPE(CRTM_Geometry_type)   , INTENT(IN)     :: geo(:)
    TYPE(CRTM_ChannelInfo_type), INTENT(IN)     :: chinfo(:)
    ! Parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Test_CRTM_FWD'
    ! Variables
    CHARACTER(256) :: err_msg, alloc_msg, io_msg, utest_msg
    CHARACTER(256) :: file_Ref
    CHARACTER(256) :: test_failure_file
    CHARACTER(7) :: file_status
    LOGICAL :: ref_data_exists
    INTEGER :: fid
    INTEGER :: err_stat, alloc_stat, io_stat
    INTEGER :: i, m
    INTEGER :: n_profiles
    INTEGER :: n_channels
    INTEGER :: dtb_maxloc
    INTEGER :: n_failed
    REAL(fp) :: fwd_tolerance
    REAL(fp), ALLOCATABLE :: tb_Ref(:)
    REAL(fp), ALLOCATABLE :: tb_FWD(:)
    REAL(fp), ALLOCATABLE :: dtb(:)
    TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_Ref(:,:)
    TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_FWD(:,:)
    TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm_FWD(:)
    TYPE(CRTM_Surface_type),    ALLOCATABLE :: sfc_FWD(:)


    ! Get test dimensions
    n_channels = SUM(CRTM_ChannelInfo_n_Channels(chinfo))
    n_profiles = SIZE(atm)
    n_sensors  = SIZE(chinfo)


    ! Perform all the allocations
    ALLOCATE( rts_Ref(n_channels, n_Profiles), &
              rts_FWD( n_channels, n_Profiles), &
              atm_FWD(n_Profiles), &
              sfc_FWD(n_Profiles), &
              tb_FWD(n_channels) , &
              dtb(n_channels)    , &
              STAT   = alloc_stat, &
              ERRMSG = alloc_msg   )
    IF ( alloc_stat /= 0 ) THEN
      err_msg = 'Error allocating data structure arrays - '//TRIM(alloc_msg)
      CALL Display_Message( ROUTINE_NAME, err_msg, FAILURE ); STOP
    END IF


    ! Loop over types of radiative transfer algorithms
    rt_algorithm_loop: DO i = 1, N_RT_ALGORITHMS


      ! Setup for test failure reporting
      file_status = 'REPLACE'
      test_failure_file = TRIM(RT_ALGORITHM_NAME(i))//'.test_failure_report'
      

      ! Output info
      ! ...Algorithm identifier
      WRITE(*,'(30("*"),1x,"FWD Comparisons for RT Algorithm ",a,1x,30("*"))') TRIM(RT_ALGORITHM_NAME(i))
      ! ...Sensors to process
      WRITE(*,'(4x,"- Sensors: ",99(a,:))') chinfo%sensor_id


      ! Copy the inputs so they can be modified if necessary
      atm_FWD = atm
      sfc_FWD = sfc


      ! Specify the RT algorithm to be used
      opt(:)%RT_Algorithm_Id = RT_ALGORITHM_ID(i)
      ! ...For emission algorithm calculations
      IF ( i == EMISSION_INDEX ) THEN
        atm_FWD%n_Clouds   = 0
        atm_FWD%n_Aerosols = 0
      END IF


      ! Perform forward calculations
      WRITE(*, '(4x,"- Running forward model...")')
      CALL Timing_Begin(timing)
      err_stat = CRTM_Forward( atm_FWD      , &
                               sfc_FWD      , &
                               geo          , &
                               chinfo       , &
                               rts_FWD      , &
                               Options = opt  )
      CALL Timing_End(timing)
      IF ( err_stat /= SUCCESS ) THEN
        err_msg = 'Error when calling CRTM_Forward'
        CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
      END IF
      CALL Timing_Display(timing)


      ! Read the algorithm reference data
      IF( n_Stokes == 1 ) THEN
      file_Ref = 'Results/'//TRIM(Sensor_ID)//'.'//TRIM(RT_ALGORITHM_NAME(i))//'.RTSolution.bin'
      ELSE
      file_Ref = 'Results/'//TRIM(Sensor_ID)//'.'//TRIM(RT_ALGORITHM_NAME(i))//'.VectorRTSolution.bin'      
      END IF
      ref_data_exists = File_Exists(file_Ref)

      ! ...If the reference data file doesn't exist, create it.
      !    This should only happen once, the first time this procedure is called.
      !    However, the write call is left in here in case the reference files
      !    need to be recreated.
      IF ( .NOT. ref_data_exists ) THEN
        err_msg = 'Reference datafile '//TRIM(file_Ref)//' does not exist. Creating it...'
        CALL Display_Message( ROUTINE_NAME, err_msg, WARNING )
        err_stat = CRTM_RTSolution_WriteFile( file_Ref, rts_FWD, Quiet=QUIET )
        IF ( err_stat /= SUCCESS ) THEN
          err_msg = 'Error writing '//TRIM(file_Ref)
          CALL Display_Message( ROUTINE_NAME, err_msg, err_stat ); STOP
        END IF
      END IF
      ! ...Now read the guaranteed to exist reference data file.
      err_stat = CRTM_RTSolution_ReadFile( file_Ref, rts_Ref, Quiet=QUIET )
      IF ( err_stat /= SUCCESS ) THEN
        err_msg = 'Error reading '//TRIM(file_Ref)
        CALL Display_Message( ROUTINE_NAME, err_msg, err_stat ); STOP
      END IF


      ! Initialise test
      WRITE(utest_msg,'("FWD comparison test | ",&
                       &"Algorithm: ",a)') &
            TRIM(RT_ALGORITHM_NAME(i))
      CALL UnitTest_Setup(utest,TRIM(utest_msg))


      ! Perform the tests profile by profile
      profile_loop: DO m = 1, n_profiles


        ! Compute the test quantities
        tb_Ref = rts_Ref(:,m)%Brightness_Temperature
        tb_FWD = rts_FWD(:,m)%Brightness_Temperature

!        write(6,'(6E14.6)') tb_FWD
!        write(6,'(6E14.6)') tb_Ref
        fwd_tolerance = SPACING(MAX(MAXVAL(tb_Ref),MAXVAL(tb_FWD))) * FWD_ULP
        dtb = ABS(tb_Ref - tb_FWD)


        ! Apply the test
        CALL UnitTest_Assert( utest, ALL(dtb < fwd_tolerance) )


        ! Output info for failed test
        IF ( UnitTest_Failed(utest) ) THEN
          ! Open failure report file
          fid = Get_Lun()
          OPEN( fid, FILE     = test_failure_file, &
                     FORM     = 'FORMATTED', &
                     STATUS   = file_status, &
                     POSITION = 'APPEND', &
                     IOSTAT   = io_stat, &
                     IOMSG    = io_msg )
          IF ( io_stat /= 0 ) THEN
            err_msg = 'Error opening '//TRIM(test_failure_file)//' - '//TRIM(io_msg)
            CALL Display_Message( ROUTINE_NAME, err_msg, WARNING )
          ELSE
            ! Update file status
            file_status = 'OLD'
            ! Report failure
            n_failed = COUNT(.NOT. (dtb < fwd_tolerance))
            dtb_maxloc = MAXLOC(dtb, DIM=1)
            WRITE(fid,'(40("="))')
            WRITE(fid,'("*** dtb, FWD_TOLERANCE test ", &
                       &"failed for profile #",i0)') m
            WRITE(fid,'("Number of failures: ",i0," of ",i0)') n_failed, n_channels
            WRITE(fid,'("Values for largest magnitude failure:")')
            WRITE(fid,'("tb_Ref                ",f19.15)') tb_Ref(dtb_maxloc)
            WRITE(fid,'("tb_FWD              - ",f19.15)') tb_FWD(dtb_maxloc)
            WRITE(fid,'("                      ",19("-"))')
            WRITE(fid,'("***--->>> delta     = ",f19.15,2x,"(",es22.15,")")') dtb(dtb_maxloc), dtb(dtb_maxloc)
            WRITE(fid,'("***--->>> threshold = ",f19.15,2x,"(",es22.15,")")') fwd_tolerance  , fwd_tolerance
            CALL CRTM_RTSolution_Inspect(rts_Ref(dtb_maxloc,m)-rts_FWD(dtb_maxloc,m), Unit=fid)
            WRITE(fid,'("*** dtb, FWD_TOLERANCE test ", &
                     &"failed for profile #",i0)') m
            WRITE(fid,'(40("="),/)')
            CLOSE(fid)
          END IF
        END IF

      END DO profile_loop

      CALL UnitTest_Report(utest)

    END DO rt_algorithm_loop


    ! Cleanup
    DEALLOCATE( rts_Ref, &
                rts_FWD, &
                atm_FWD, &
                sfc_FWD, &
                tb_Ref , &
                tb_FWD , &
                dtb    , &
                STAT   = alloc_stat )

  END SUBROUTINE Test_CRTM_FWD


  ! ====================
  ! Tangent-linear model
  ! ====================
  SUBROUTINE Test_CRTM_TL(utest, atm, sfc, geo, chinfo)
    ! Arguments
    TYPE(UnitTest_type)        , INTENT(IN OUT) :: utest
    TYPE(CRTM_Atmosphere_type) , INTENT(IN)     :: atm(:)
    TYPE(CRTM_Surface_type)    , INTENT(IN)     :: sfc(:)
    TYPE(CRTM_Geometry_type)   , INTENT(IN)     :: geo(:)
    TYPE(CRTM_ChannelInfo_type), INTENT(IN)     :: chinfo(:)
    ! Parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Test_CRTM_TL'
    ! Variables
    CHARACTER(256) :: err_msg, alloc_msg, io_msg, utest_msg
    CHARACTER(256) :: file_Ref
    CHARACTER(256) :: test_failure_file
    CHARACTER(7) :: file_status
    LOGICAL :: ref_data_exists
    INTEGER :: fid
    INTEGER :: err_stat, alloc_stat, io_stat
    INTEGER :: i, m
    INTEGER :: n_profiles
    INTEGER :: n_channels
    INTEGER :: dtb_tl_maxloc
    INTEGER :: n_failed
    REAL(fp) :: tl_tolerance
    REAL(fp), ALLOCATABLE :: tb_TL_Ref(:)
    REAL(fp), ALLOCATABLE :: tb_TL(:)
    REAL(fp), ALLOCATABLE :: dtb_TL(:)
    TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_TL_Ref(:,:)
    TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_FWD(:,:), rts_TL(:,:)
    TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm_FWD(:), atm_TL(:)
    TYPE(CRTM_Surface_type),    ALLOCATABLE :: sfc_FWD(:), sfc_TL(:)


    ! Get test dimensions
    n_channels = SUM(CRTM_ChannelInfo_n_Channels(chinfo))
    n_profiles = SIZE(atm)
    n_sensors  = SIZE(chinfo)


    ! Perform all the allocations
    ALLOCATE( rts_TL_Ref(n_channels, n_Profiles), &
              rts_FWD(   n_channels, n_Profiles), &
              rts_TL(    n_channels, n_Profiles), &
              atm_FWD(n_Profiles)  , &
              atm_TL( n_Profiles)  , &
              sfc_FWD(n_Profiles)  , &
              sfc_TL( n_Profiles)  , &
              tb_TL_Ref(n_channels), &
              tb_TL( n_channels)   , &
              dtb_TL(n_channels)   , &
              STAT   = alloc_stat  , &
              ERRMSG = alloc_msg     )
    IF ( alloc_stat /= 0 ) THEN
      err_msg = 'Error allocating data structure arrays - '//TRIM(alloc_msg)
      CALL Display_Message( ROUTINE_NAME, err_msg, FAILURE ); STOP
    END IF


    ! Loop over types of radiative transfer algorithms
    rt_algorithm_loop: DO i = 1, N_RT_ALGORITHMS


      ! Setup for test failure reporting
      file_status = 'REPLACE'
      test_failure_file = TRIM(RT_ALGORITHM_NAME(i))//'.TL.test_failure_report'
      
      
      ! Output info
      ! ...Algorithm identifier
      WRITE(*,'(30("*"),1x,"TL Comparisons for RT Algorithm ",a,1x,30("*"))') TRIM(RT_ALGORITHM_NAME(i))
      ! ...Sensors to process
      WRITE(*, '(4x,"- Sensors: ",99(a,:))') chinfo%sensor_id


      ! Copy the inputs so they can be modified if necessary
      atm_FWD = atm
      sfc_FWD = sfc


      ! Specify the RT algorithm to be used
      opt(:)%RT_Algorithm_Id = RT_ALGORITHM_ID(i)
      ! ...For emission algorithm calculations
      IF ( i == EMISSION_INDEX ) THEN
        atm_FWD%n_Clouds   = 0
        atm_FWD%n_Aerosols = 0
      END IF


      ! Assign the tangent-linear perturbations
      CALL Assign_TL_Atmosphere( DX, atm_FWD, atm_TL )
      CALL Assign_TL_Surface( DX, sfc_FWD, sfc_TL )


      ! Perform tangent-linear calculations
      WRITE(*, '(4x,"- Running tangent-linear model...")')
      CALL Timing_Begin(timing)
      err_stat = CRTM_Tangent_Linear( atm_FWD      , &
                                      sfc_FWD      , &
                                      atm_TL       , &
                                      sfc_TL       , &
                                      geo          , &
                                      chinfo       , &
                                      rts_FWD      , &
                                      rts_TL       , &
                                      Options = opt  )
      CALL Timing_End(timing)
      IF ( err_stat /= SUCCESS ) THEN
        err_msg = 'Error when calling CRTM_Tangent_Linear'
        CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
      END IF
      CALL Timing_Display(timing)


      ! Read the algorithm reference data
      IF( n_Stokes == 1 ) THEN
      file_Ref = 'Results/'//TRIM(Sensor_ID)//'.'//TRIM(RT_ALGORITHM_NAME(i))//'.TL.RTSolution.bin'
      ELSE
      file_Ref = 'Results/'//TRIM(Sensor_ID)//'.'//TRIM(RT_ALGORITHM_NAME(i))//'.TL.VectorRTSolution.bin'      
      END IF

      ref_data_exists = File_Exists(file_Ref)
      ! ...If the reference data file doesn't exist, create it.
      !    This should only happen once, the first time this procedure is called.
      !    However, the write call is left in here in case the reference files
      !    need to be recreated.
      IF ( .NOT. ref_data_exists ) THEN
        err_msg = 'Reference datafile '//TRIM(file_Ref)//' does not exist. Creating it...'
        CALL Display_Message( ROUTINE_NAME, err_msg, WARNING )
        err_stat = CRTM_RTSolution_WriteFile( file_Ref, rts_TL, Quiet=QUIET )
        IF ( err_stat /= SUCCESS ) THEN
          err_msg = 'Error writing '//TRIM(file_Ref)
          CALL Display_Message( ROUTINE_NAME, err_msg, err_stat ); STOP
        END IF
      END IF
      ! ...Now read the guaranteed to exist reference data file.
      err_stat = CRTM_RTSolution_ReadFile( file_Ref, rts_TL_Ref, Quiet=QUIET )
      IF ( err_stat /= SUCCESS ) THEN
        err_msg = 'Error reading '//TRIM(file_Ref)
        CALL Display_Message( ROUTINE_NAME, err_msg, err_stat ); STOP
      END IF


      ! Initialise test
      WRITE(utest_msg,'("TL comparison test | ",&
                       &"Algorithm: ",a)') &
            TRIM(RT_ALGORITHM_NAME(i))
      CALL UnitTest_Setup(utest,TRIM(utest_msg))


      ! Perform the tests profile by profile
      profile_loop: DO m = 1, n_profiles


        ! Compute the test quantities
        tb_TL_Ref = rts_TL_Ref(:,m)%Brightness_Temperature
        tb_TL     = rts_TL(:,m)%Brightness_Temperature
        tl_tolerance = SPACING(MAX(MAXVAL(tb_TL_Ref),MAXVAL(tb_TL))) * TL_ULP
        dtb_TL = ABS(tb_TL_Ref - tb_TL)


        ! Apply the test
        CALL UnitTest_Assert( utest, ALL(dtb_TL < tl_tolerance) )


        ! Output info for failed test
        IF ( UnitTest_Failed(utest) ) THEN
          ! Open failure report file
          fid = Get_Lun()
          OPEN( fid, FILE     = test_failure_file, &
                     FORM     = 'FORMATTED', &
                     STATUS   = file_status, &
                     POSITION = 'APPEND', &
                     IOSTAT   = io_stat, &
                     IOMSG    = io_msg )
          IF ( io_stat /= 0 ) THEN
            err_msg = 'Error opening '//TRIM(test_failure_file)//' - '//TRIM(io_msg)
            CALL Display_Message( ROUTINE_NAME, err_msg, WARNING )
          ELSE
            ! Update file status
            file_status = 'OLD'
            ! Report failure
            n_failed = COUNT(.NOT. (dtb_TL < tl_tolerance))
            dtb_TL_maxloc = MAXLOC(dtb_TL, DIM=1)
            WRITE(fid,'(40("="))')
            WRITE(fid,'("*** dtb_TL, TL_TOLERANCE test ", &
                     &"failed for profile #",i0)') m
            WRITE(fid,'("Number of failures: ",i0," of ",i0)') n_failed, n_channels
            WRITE(fid,'("Values for largest magnitude failure:")')
            WRITE(fid,'("tb_TL_Ref             ",f19.15)') tb_TL_Ref(dtb_TL_maxloc)
            WRITE(fid,'("tb_TL               - ",f19.15)') tb_TL(dtb_TL_maxloc)
            WRITE(fid,'("                      ",19("-"))')
            WRITE(fid,'("***--->>> delta     = ",f19.15,2x,"(",es22.15,")")') dtb_TL(dtb_TL_maxloc), dtb_TL(dtb_TL_maxloc)
            WRITE(fid,'("***--->>> threshold = ",f19.15,2x,"(",es22.15,")")') tl_tolerance  , tl_tolerance
            CALL CRTM_RTSolution_Inspect(rts_TL_Ref(dtb_TL_maxloc,m)-rts_TL(dtb_TL_maxloc,m), Unit=fid)
            WRITE(fid,'("*** dtb_TL, TL_TOLERANCE test ", &
                     &"failed for profile #",i0)') m
            WRITE(fid,'(40("="),/)')
            CLOSE(fid)
          END IF
        END IF

      END DO profile_loop

      CALL UnitTest_Report(utest)

    END DO rt_algorithm_loop


    ! Cleanup
    DEALLOCATE( rts_TL_Ref, &
                rts_FWD   , &
                rts_TL    , &
                atm_FWD   , &
                atm_TL    , &
                sfc_FWD   , &
                sfc_TL    , &
                tb_TL_Ref , &
                tb_TL     , &
                dtb_TL    , &
                STAT = alloc_stat )

  END SUBROUTINE Test_CRTM_TL


  ! =============
  ! Adjoint model
  ! =============
  SUBROUTINE Test_CRTM_AD(utest, atm, sfc, geo, chinfo)
    ! Arguments
    TYPE(UnitTest_type)        , INTENT(IN OUT) :: utest
    TYPE(CRTM_Atmosphere_type) , INTENT(IN)     :: atm(:)
    TYPE(CRTM_Surface_type)    , INTENT(IN)     :: sfc(:)
    TYPE(CRTM_Geometry_type)   , INTENT(IN)     :: geo(:)
    TYPE(CRTM_ChannelInfo_type), INTENT(IN)     :: chinfo(:)
    ! Parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Test_CRTM_AD'
    ! Variables
    CHARACTER(256) :: err_msg, alloc_msg, io_msg, utest_msg
    CHARACTER(256) :: atm_file_Ref, sfc_file_Ref
    CHARACTER(256) :: atm_failure_file, sfc_failure_file
    CHARACTER(7) :: atm_file_status, sfc_file_status
    LOGICAL :: atm_ref_data_exists, sfc_ref_data_exists
    INTEGER :: fid
    INTEGER :: err_stat, alloc_stat, io_stat
    INTEGER :: i, m
    INTEGER :: n_channels
    INTEGER :: n_profiles
    INTEGER :: n_sensors
    TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_FWD(:,:)
    TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_AD(:,:)
    TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm_FWD(:)
    TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm_AD(:)
    TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm_AD_Ref(:)
    TYPE(CRTM_Surface_type),    ALLOCATABLE :: sfc_FWD(:)
    TYPE(CRTM_Surface_type),    ALLOCATABLE :: sfc_AD(:)
    TYPE(CRTM_Surface_type),    ALLOCATABLE :: sfc_AD_Ref(:)

    ! Get test dimensions
    n_channels = SUM(CRTM_ChannelInfo_n_Channels(chinfo))
    n_profiles = SIZE(atm)
    n_sensors  = SIZE(chinfo)


    ! Perform all the allocations
    ALLOCATE( rts_FWD(n_channels, n_Profiles), &
              rts_AD( n_channels, n_Profiles), &
              atm_FWD(   n_Profiles), &
              atm_AD(    n_Profiles), &
              atm_AD_Ref(   n_Profiles), &
              sfc_FWD(   n_Profiles), &
              sfc_AD(    n_Profiles), &
              sfc_AD_Ref(n_Profiles), &
              STAT   = alloc_stat, &
              ERRMSG = alloc_msg   )
    IF ( alloc_stat /= 0 ) THEN
      err_msg = 'Error allocating data structure arrays - '//TRIM(alloc_msg)
      CALL Display_Message( ROUTINE_NAME, err_msg, FAILURE ); STOP
    END IF


    ! Loop over types of radiative transfer algorithms
    rt_algorithm_loop: DO i = 1, N_RT_ALGORITHMS


      ! Setup for test failure reporting
      atm_file_status = 'REPLACE'
      atm_failure_file = TRIM(RT_ALGORITHM_NAME(i))//'.AD.atmosphere.test_failure_report'
      sfc_file_status = 'REPLACE'
      sfc_failure_file = TRIM(RT_ALGORITHM_NAME(i))//'.AD.surface.test_failure_report'
      

      ! Output info
      ! ...Algorithm identifier
      WRITE(*,'(30("*"),1x,"AD Comparisons for RT Algorithm ",a,1x,30("*"))') TRIM(RT_ALGORITHM_NAME(i))
      ! ...Sensors to process
      WRITE(*, '(4x,"- Sensors: ",99(a,:))') chinfo%sensor_id


      ! Copy the inputs so they can be modified if necessary
      atm_FWD = atm
      sfc_FWD = sfc


      ! Specify the RT algorithm to be used
      opt(:)%RT_Algorithm_Id = RT_ALGORITHM_ID(i)
      ! ...For emission algorithm calculations
      IF ( i == EMISSION_INDEX ) THEN
        atm_FWD%n_Clouds   = 0
        atm_FWD%n_Aerosols = 0
      END IF


      ! Assign the adjoint data
      ! ...The input
      CALL CRTM_RTSolution_Zero( rts_AD ); rts_AD%Brightness_Temperature = ONE
      ! ...The output
      atm_AD = atm_FWD; CALL CRTM_Atmosphere_Zero( atm_AD )
      sfc_AD = sfc_FWD; CALL CRTM_Surface_Zero( sfc_AD )


      ! Perform the adjoint calculations
      WRITE(*, '(4x,"- Running adjoint model...")')
      CALL Timing_Begin(timing)
      err_stat = CRTM_Adjoint( atm_FWD      , &
                               sfc_FWD      , &
                               rts_AD       , &
                               geo          , &
                               chinfo       , &
                               atm_AD       , &
                               sfc_AD       , &
                               rts_FWD      , &
                               Options = opt  )
      CALL Timing_End(timing)
      IF ( err_stat /= SUCCESS ) THEN
        err_msg = 'Error when calling CRTM_Adjoint'
        CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
      END IF
      CALL Timing_Display(timing)


      ! Read the algorithm ATMOSPHERE reference data
      IF( n_Stokes == 1 ) THEN
      atm_file_Ref = 'Results/'//TRIM(Sensor_ID)//'.'//TRIM(RT_ALGORITHM_NAME(i))//'.AD.Atmosphere.bin'
      ELSE
      atm_file_Ref = 'Results/'//TRIM(Sensor_ID)//'.'//TRIM(RT_ALGORITHM_NAME(i))//'.AD.VectorAtmosphere.bin'      
      END IF


      atm_ref_data_exists = File_Exists(atm_file_Ref)
      ! ...If the reference data file doesn't exist, create it.
      !    This should only happen once, the first time this procedure is called.
      !    However, the write call is left in here in case the reference files
      !    need to be recreated.
      IF ( .NOT. atm_ref_data_exists ) THEN
        err_msg = 'Reference datafile '//TRIM(atm_file_Ref)//' does not exist. Creating it...'
        CALL Display_Message( ROUTINE_NAME, err_msg, WARNING )
        err_stat = CRTM_Atmosphere_WriteFile( atm_file_Ref, atm_AD, Quiet=QUIET )
        IF ( err_stat /= SUCCESS ) THEN
          err_msg = 'Error writing '//TRIM(atm_file_Ref)
          CALL Display_Message( ROUTINE_NAME, err_msg, err_stat ); STOP
        END IF
      END IF
      ! ...Now read the guaranteed to exist reference data file.
      err_stat = CRTM_Atmosphere_ReadFile( atm_file_Ref, atm_AD_Ref, Quiet=QUIET )
      IF ( err_stat /= SUCCESS ) THEN
        err_msg = 'Error reading '//TRIM(atm_file_Ref)
        CALL Display_Message( ROUTINE_NAME, err_msg, err_stat ); STOP
      END IF


      ! Read the algorithm surface reference data
      IF( n_Stokes == 1 ) THEN
      sfc_file_Ref = 'Results/'//TRIM(Sensor_ID)//'.'//TRIM(RT_ALGORITHM_NAME(i))//'.AD.Surface.bin'
      ELSE
      sfc_file_Ref = 'Results/'//TRIM(Sensor_ID)//'.'//TRIM(RT_ALGORITHM_NAME(i))//'.AD.VectorSurface.bin'      
      END IF
      sfc_ref_data_exists = File_Exists(sfc_file_Ref)
      ! ...If the reference data file doesn't exist, create it.
      !    This should only happen once, the first time this procedure is called.
      !    However, the write call is left in here in case the reference files
      !    need to be recreated.
      IF ( .NOT. sfc_ref_data_exists ) THEN
        err_msg = 'Reference datafile '//TRIM(sfc_file_Ref)//' does not exist. Creating it...'
        CALL Display_Message( ROUTINE_NAME, err_msg, WARNING )
        err_stat = CRTM_Surface_WriteFile( sfc_file_Ref, sfc_AD, Quiet=QUIET )
        IF ( err_stat /= SUCCESS ) THEN
          err_msg = 'Error writing '//TRIM(sfc_file_Ref)
          CALL Display_Message( ROUTINE_NAME, err_msg, err_stat ); STOP
        END IF
      END IF
      ! ...Now read the guaranteed to exist reference data file.
      err_stat = CRTM_Surface_ReadFile( sfc_file_Ref, sfc_AD_Ref, Quiet=QUIET )
      IF ( err_stat /= SUCCESS ) THEN
        err_msg = 'Error reading '//TRIM(sfc_file_Ref)
        CALL Display_Message( ROUTINE_NAME, err_msg, err_stat ); STOP
      END IF


      ! Initialise test
      WRITE(utest_msg,'("AD comparison test | ",&
                       &"Algorithm: ",a)') &
            TRIM(RT_ALGORITHM_NAME(i))
      CALL UnitTest_Setup(utest,TRIM(utest_msg))


      ! Perform the tests profile by profile
      profile_loop: DO m = 1, n_profiles


        ! Perform the atmosphere test
        CALL UnitTest_Assert( utest, CRTM_Atmosphere_Compare(atm_AD_Ref(m), atm_AD(m), n_SigFig=AD_SIGFIG) )


        ! Output info for failed test
        IF ( UnitTest_Failed(utest) ) THEN
          ! Open failure report file
          fid = Get_Lun()
          OPEN( fid, FILE     = atm_failure_file, &
                     FORM     = 'FORMATTED', &
                     STATUS   = atm_file_status, &
                     POSITION = 'APPEND', &
                     IOSTAT   = io_stat, &
                     IOMSG    = io_msg )
          IF ( io_stat /= 0 ) THEN
            err_msg = 'Error opening '//TRIM(atm_failure_file)//' - '//TRIM(io_msg)
            CALL Display_Message( ROUTINE_NAME, err_msg, WARNING )
          ELSE
            ! Update file status
            atm_file_status = 'OLD'
            ! Report failure
            WRITE(fid,'(40("="))')
            WRITE(fid,'("*** ATM AD test failed for profile #",i0)') m
            WRITE(fid,'("No. of significant figures used in comparison: ", i0)') AD_SIGFIG
            CALL CRTM_Atmosphere_Inspect(atm_AD_Ref(m)-atm_AD(m), Unit=fid)
            WRITE(fid,'("*** ATM AD test failed for profile #",i0)') m
            WRITE(fid,'(40("="),/)')
            CLOSE(fid)
          END IF
        END IF


        ! Perform the surface test
        CALL UnitTest_Assert( utest, CRTM_Surface_Compare(sfc_AD_Ref(m), sfc_AD(m), n_SigFig=AD_SIGFIG) )


        ! Output info for failed test
        IF ( UnitTest_Failed(utest) ) THEN
          ! Open failure report file
          fid = Get_Lun()
          OPEN( fid, FILE     = sfc_failure_file, &
                     FORM     = 'FORMATTED', &
                     STATUS   = sfc_file_status, &
                     POSITION = 'APPEND', &
                     IOSTAT   = io_stat, &
                     IOMSG    = io_msg )
          IF ( io_stat /= 0 ) THEN
            err_msg = 'Error opening '//TRIM(sfc_failure_file)//' - '//TRIM(io_msg)
            CALL Display_Message( ROUTINE_NAME, err_msg, WARNING )
          ELSE
            ! Update file status
            sfc_file_status = 'OLD'
            ! Report failure
            WRITE(fid,'(40("="))')
            WRITE(fid,'("*** SFC AD test failed for profile #",i0)') m
            WRITE(fid,'("No. of significant figures used in comparison: ", i0)') AD_SIGFIG
            CALL CRTM_Surface_Inspect(sfc_AD_Ref(m)-sfc_AD(m), Unit=fid)
            WRITE(fid,'("*** SFC AD test failed for profile #",i0)') m
            WRITE(fid,'(40("="),/)')
            CLOSE(fid)
          END IF
        END IF

      END DO profile_loop

      CALL UnitTest_Report(utest)

    END DO rt_algorithm_loop


    ! Cleanup
    DEALLOCATE( rts_FWD   , &
                rts_AD    , &
                atm_FWD   , &
                atm_AD    , &
                atm_AD_Ref, &
                sfc_FWD   , &
                sfc_AD    , &
                sfc_AD_Ref, &
                STAT = alloc_stat )

  END SUBROUTINE Test_CRTM_AD


  ! ==================================
  ! Forward/tangent-linear consistency
  ! ==================================
  SUBROUTINE Test_CRTM_FWDTL(utest, atm, sfc, geo, chinfo)
    ! Arguments
    TYPE(UnitTest_type)        , INTENT(IN OUT) :: utest
    TYPE(CRTM_Atmosphere_type) , INTENT(IN)     :: atm(:)
    TYPE(CRTM_Surface_type)    , INTENT(IN)     :: sfc(:)
    TYPE(CRTM_Geometry_type)   , INTENT(IN)     :: geo(:)
    TYPE(CRTM_ChannelInfo_type), INTENT(IN)     :: chinfo(:)
    ! Parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Test_CRTM_FWDTL'
    ! ...Components to test
    INTEGER     , PARAMETER :: N_ATM_COMPONENTS = 3
    CHARACTER(*), PARAMETER :: ATM_COMPONENT_NAME(N_ATM_COMPONENTS) = &
      [ 'Temperature', &
        'Water Vapor', &
        'Ozone      '  ]
    INTEGER     , PARAMETER :: N_CLOUD_COMPONENTS = 2
    CHARACTER(*), PARAMETER :: CLOUD_COMPONENT_NAME(N_CLOUD_COMPONENTS) = &
      [ 'Effective Radius', &
        'Water Content   '  ]
    INTEGER     , PARAMETER :: N_AEROSOL_COMPONENTS = 2
    CHARACTER(*), PARAMETER :: AEROSOL_COMPONENT_NAME(N_AEROSOL_COMPONENTS) = &
      [ 'Effective Radius', &
        'Concentration   '  ]
!  There are no TL/AD part in land surface model
    INTEGER     , PARAMETER :: N_LAND_SFC_COMPONENTS = 1
    CHARACTER(*), PARAMETER :: LAND_SFC_COMPONENT_NAME(N_LAND_SFC_COMPONENTS) = &
      [ 'Land Temperature     ' ]

    INTEGER     , PARAMETER :: N_WATER_SFC_COMPONENTS = 4
    CHARACTER(*), PARAMETER :: WATER_SFC_COMPONENT_NAME(N_WATER_SFC_COMPONENTS) = &
      [ 'Water Temperature', &
        'Wind Speed       ', &
        'Wind Direction   ', &
        'Salinity         '  ]
!  There are no TL/AD part in snow surface model
    INTEGER     , PARAMETER :: N_SNOW_SFC_COMPONENTS = 1
    CHARACTER(*), PARAMETER :: SNOW_SFC_COMPONENT_NAME(N_SNOW_SFC_COMPONENTS) = &
      [ 'Snow_Temperature' ]   
!  There are no TL/AD part in ice surface model
    INTEGER     , PARAMETER :: N_ICE_SFC_COMPONENTS = 1
    CHARACTER(*), PARAMETER :: ICE_SFC_COMPONENT_NAME(N_ICE_SFC_COMPONENTS) = &
      [ 'Ice_Temperature' ]

    ! ...Atmosphere layer skip value. No need to test all layers.
    INTEGER, PARAMETER :: LAYER_STEP = 20
    ! ...Cloud layer begin and step value
    INTEGER, PARAMETER :: CLOUD_LAYER_BEGIN = 40
    INTEGER, PARAMETER :: CLOUD_LAYER_STEP  = 5
    ! ...Aerosol layer begin and step value
    INTEGER, PARAMETER :: AEROSOL_LAYER_BEGIN = 40
    INTEGER, PARAMETER :: AEROSOL_LAYER_STEP  = 5
    ! Variables
    CHARACTER(256) :: err_msg, alloc_msg, io_msg, utest_msg
    CHARACTER(256) :: atm_failure_file, cloud_failure_file, aerosol_failure_file
    CHARACTER(256) :: landsfc_failure_file, watersfc_failure_file, snowsfc_failure_file, icesfc_failure_file
    CHARACTER(7) :: atm_file_status, cloud_file_status, aerosol_file_status
    CHARACTER(7) :: landsfc_file_status, watersfc_file_status, snowsfc_file_status, icesfc_file_status
    INTEGER :: fid
    INTEGER :: err_stat, alloc_stat, io_stat
    INTEGER :: i, k, m, ic, ia
    INTEGER :: ialpha
    INTEGER :: icomponent
    INTEGER :: n_channels
    INTEGER :: n_profiles
    INTEGER :: n_sensors
    INTEGER :: h2o_idx, o3_idx
    INTEGER :: n_failed
    INTEGER :: dtb_maxloc
    REAL(fp) :: delta
    REAL(fp), ALLOCATABLE :: tb_NLm(:), tb_NLp(:), dtb_NL(:)
    REAL(fp), ALLOCATABLE :: tb_TL(:)
    REAL(fp), ALLOCATABLE :: dtb_delta(:)
    TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_Base(:,:)
    TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_FWD(:,:)
    TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_TL(:,:)
    TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_NLm(:,:)
    TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_NLp(:,:)
    TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm_FWD(:)
    TYPE(CRTM_Surface_type)   , ALLOCATABLE :: sfc_FWD(:)
    TYPE(CRTM_Atmosphere_type) :: atm_TL(1)
    TYPE(CRTM_Atmosphere_type) :: atm_NLm(1)
    TYPE(CRTM_Atmosphere_type) :: atm_NLp(1)
    TYPE(CRTM_Surface_type)    :: sfc_TL(1)
    TYPE(CRTM_Surface_type)    :: sfc_NLp(1)
    TYPE(CRTM_Surface_type)    :: sfc_NLm(1)


    ! Get test dimensions
    n_channels = SUM(CRTM_ChannelInfo_n_Channels(chinfo))
    n_profiles = SIZE(atm)
    n_sensors  = SIZE(chinfo)

              dtb_maxloc = 1
    ! Perform all the allocations
    ALLOCATE( rts_Base(n_channels, n_Profiles), &
              rts_FWD( n_channels, n_Profiles), &
              rts_TL(  n_channels, n_Profiles), &
              rts_NLp( n_channels, n_Profiles), &
              rts_NLm( n_channels, n_Profiles), &
              tb_NLm(   n_channels), &
              tb_NLp(   n_channels), &
              dtb_NL(   n_channels), &
              tb_TL(    n_channels), &
              dtb_delta(n_channels), &
              atm_FWD(n_Profiles), &
              sfc_FWD(n_Profiles), &
              STAT   = alloc_stat, &
              ERRMSG = alloc_msg   )
    IF ( alloc_stat /= 0 ) THEN
      err_msg = 'Error allocating data structure arrays - '//TRIM(alloc_msg)
      CALL Display_Message( ROUTINE_NAME, err_msg, FAILURE ); STOP
    END IF


    ! Loop over types of radiative transfer algorithms
    rt_algorithm_loop: DO i = 1, N_RT_ALGORITHMS            


      ! Setup for test failure reporting
      atm_file_status      = 'REPLACE'
      cloud_file_status    = 'REPLACE'
      aerosol_file_status  = 'REPLACE'
      landsfc_file_status  = 'REPLACE'
      watersfc_file_status = 'REPLACE'
      snowsfc_file_status  = 'REPLACE'
      icesfc_file_status   = 'REPLACE'
      
      atm_failure_file      = TRIM(RT_ALGORITHM_NAME(i))//'.FWDTL.atmosphere.test_failure_report'
      cloud_failure_file    = TRIM(RT_ALGORITHM_NAME(i))//'.FWDTL.cloud.test_failure_report'
      aerosol_failure_file  = TRIM(RT_ALGORITHM_NAME(i))//'.FWDTL.aerosol.test_failure_report'
      landsfc_failure_file  = TRIM(RT_ALGORITHM_NAME(i))//'.FWDTL.landsfc.test_failure_report'
      watersfc_failure_file = TRIM(RT_ALGORITHM_NAME(i))//'.FWDTL.watersfc.test_failure_report'
      snowsfc_failure_file  = TRIM(RT_ALGORITHM_NAME(i))//'.FWDTL.snowsfc.test_failure_report'
      icesfc_failure_file   = TRIM(RT_ALGORITHM_NAME(i))//'.FWDTL.icesfc.test_failure_report'
      

      ! Output info
      ! ...Algorithm identifier
      WRITE(*,'(30("*"),1x,"FWD/TL Comparisons for RT Algorithm ",a,1x,30("*"))') TRIM(RT_ALGORITHM_NAME(i))
      ! ...Sensors to process
      WRITE(*, '(4x,"- Sensors: ",99(a,:))') chinfo%sensor_id


      ! Copy the inputs so they can be modified if necessary
      atm_FWD = atm
      sfc_FWD = sfc


      ! Specify the RT algorithm to be used
      opt(:)%RT_Algorithm_Id = RT_ALGORITHM_ID(i)
      ! ...For emission algorithm calculations
      IF ( i == EMISSION_INDEX ) THEN
        atm_FWD%n_Clouds   = 0
        atm_FWD%n_Aerosols = 0
      END IF


      ! Loop over perturbations to be tested
      alpha_loop: DO ialpha = 1, N_ALPHA


        ! The perturbation value
        delta = ALPHA(ialpha) * DX


        ! Loop over profiles
        profile_loop: DO m = 1, n_profiles


          ! ====================================
          ! =====     ATMOSPHERE TESTS     =====
          ! ====================================

          ! Get the absorber indices for the current profile
          h2o_idx = CRTM_Get_AbsorberIdx(atm_FWD(m), H2O_ID)
          o3_idx  = CRTM_Get_AbsorberIdx(atm_FWD(m), O3_ID)


          ! Initialise the surface TL data (only has to be done once here)
          sfc_TL = sfc_FWD(m)
          CALL CRTM_Surface_Zero(sfc_TL)


          ! Loop over the components to test
          atm_component_loop: DO icomponent = 1, N_ATM_COMPONENTS


            ! Initialise test
            WRITE(utest_msg,'("ATM FWD/TL test | ",&
                             &"Profile: ",i0," | ",&
                             &"Component: ",a," | ",&
                             &"Alpha: ",es9.2," | ",&
                             &"Algorithm: ",a)') &
                  m, TRIM(ATM_COMPONENT_NAME(icomponent)), &
                  ALPHA(ialpha), &
                  TRIM(RT_ALGORITHM_NAME(i))
            CALL UnitTest_Setup(utest,TRIM(utest_msg))


            ! Loop over the atmospheric layers
            layer_loop: DO k = 1, atm_FWD(m)%n_Layers, LAYER_STEP
              WRITE(*, '(4x,"- Running tangent- and non-linear model for layer ",i0," perturbation...")') k


              ! Initialise the atm TL and NL data
              atm_NLm(1) = atm_FWD(m)
              atm_NLp(1) = atm_FWD(m)
              atm_TL(1)  = atm_FWD(m); CALL CRTM_Atmosphere_Zero(atm_TL)


              ! Select atmosphere component to test
              SELECT CASE(TRIM(ATM_COMPONENT_NAME(icomponent)))

                CASE('Temperature')
                  delta = ALPHA(ialpha) * DX
                  atm_TL(1)%Temperature(k)  = delta
                  atm_NLp(1)%Temperature(k) = atm_FWD(m)%Temperature(k) + (delta/TWO)
                  atm_NLm(1)%Temperature(k) = atm_FWD(m)%Temperature(k) - (delta/TWO)

                CASE('Water Vapor')
                  delta = ALPHA(ialpha) * DX * atm_FWD(m)%Absorber(k,h2o_idx)
                  atm_NLm(1)%Absorber(k,h2o_idx) = atm_FWD(m)%Absorber(k,h2o_idx) - (delta/TWO)
                  IF( atm_NLm(1)%Absorber(k,h2o_idx) < ZERO ) THEN
                    delta = atm_FWD(m)%Absorber(k,h2o_idx) * 0.3_fp
                    atm_NLm(1)%Absorber(k,h2o_idx) = atm_FWD(m)%Absorber(k,h2o_idx) - (delta/TWO)
                  END IF
                  atm_TL(1)%Absorber(k,h2o_idx)  = delta
                  atm_NLp(1)%Absorber(k,h2o_idx) = atm_FWD(m)%Absorber(k,h2o_idx) + (delta/TWO)

                CASE('Ozone')
                  delta = ALPHA(ialpha) * DX * atm_FWD(m)%Absorber(k,o3_idx)
                  atm_NLm(1)%Absorber(k,o3_idx) = atm_FWD(m)%Absorber(k,o3_idx) - (delta/TWO)
                  IF( atm_NLm(1)%Absorber(k,o3_idx) < ZERO ) THEN
                    delta = atm_FWD(m)%Absorber(k,o3_idx) * 0.3_fp
                    atm_NLm(1)%Absorber(k,o3_idx) = atm_FWD(m)%Absorber(k,o3_idx) - (delta/TWO)
                  END IF
                  atm_TL(1)%Absorber(k,o3_idx)  = delta
                  atm_NLp(1)%Absorber(k,o3_idx) = atm_FWD(m)%Absorber(k,o3_idx) + (delta/TWO)


                CASE DEFAULT
                  err_msg = 'Unrecognised ATM component name: '//TRIM(ATM_COMPONENT_NAME(icomponent))
                  CALL Display_Message(ROUTINE_NAME, err_msg, FAILURE ); STOP

              END SELECT


              ! Perform tangent-linear calculations
              err_stat = CRTM_Tangent_Linear( atm_FWD(m:m)      , &
                                              sfc_FWD(m:m)      , &
!                                              [atm_TL]          , &
!                                              [sfc_TL]          , &
                                              atm_TL          , &
                                              sfc_TL          , &
                                              geo(m:m)          , &
                                              chinfo            , &
                                              rts_FWD(:,m:m)    , &
                                              rts_TL(:,m:m)     , &
                                              Options = opt(m:m)  )

              IF ( err_stat /= SUCCESS ) THEN
                err_msg = 'Error when calling CRTM_Tangent_Linear'
                CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
              END IF

              ! Perform non-linear calculations
              ! ...Negative perturbation
!!             CALL CRTM_RTSolution_ZERO(rts_NLm(:,m:m))
              err_stat = CRTM_Forward( atm_NLm    , &
                                       sfc_FWD(m:m) , &
                                       geo(m:m)     , &
                                       chinfo       , &
                                       rts_NLm(:,m:m), &
                                       Options = opt(m:m) )

              IF ( err_stat /= SUCCESS ) THEN
                err_msg = 'Error when calling CRTM_Forward for negative perturbation'
                CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
              END IF
              ! ...Positive perturbation
!!              CALL CRTM_RTSolution_ZERO(rts_NLp(:,m:m))
              err_stat = CRTM_Forward( atm_NLp    , &
                                       sfc_FWD(m:m) , &
                                       geo(m:m)     , &
                                       chinfo       , &
                                       rts_NLp(:,m:m), &
                                       Options = opt(m:m) )
              IF ( err_stat /= SUCCESS ) THEN
                err_msg = 'Error when calling CRTM_Forward for positive perturbation'
                CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
              END IF


              ! Compute the test quantities
           IF( IR_MW_Sensor ) THEN
              tb_NLm = rts_NLm(:,m)%Brightness_Temperature
              tb_NLp = rts_NLp(:,m)%Brightness_Temperature
              tb_TL  = rts_TL(:,m)%Brightness_Temperature
           ELSE
             IF( n_Stokes == 1 ) THEN
               tb_NLm = rts_NLm(:,m)%Radiance
               tb_NLp = rts_NLp(:,m)%Radiance
               tb_TL  = rts_TL(:,m)%Radiance
             ELSE
               tb_NLm = rts_NLm(:,m)%Stokes(iStoke)
               tb_NLp = rts_NLp(:,m)%Stokes(iStoke)
               tb_TL  = rts_TL(:,m)%Stokes(iStoke)
             END IF
           END IF
              dtb_NL = tb_NLp - tb_NLm
              dtb_delta = ABS(dtb_NL - tb_TL)
         
              ! Apply the test
              CALL UnitTest_Assert( utest, ALL(dtb_delta < FWDTL_TOLERANCE(ialpha)) )


              ! Output info for failed test
              IF ( UnitTest_Failed(utest) ) THEN
                 print *,' TL failed ',TRIM(ATM_COMPONENT_NAME(icomponent))
                 write(6,'(6E15.8)') dtb_delta,FWDTL_TOLERANCE(ialpha),delta
                ! Open failure report file
                fid = Get_Lun()
                OPEN( fid, FILE     = atm_failure_file, &
                           FORM     = 'FORMATTED', &
                           STATUS   = atm_file_status, &
                           POSITION = 'APPEND', &
                           IOSTAT   = io_stat, &
                           IOMSG    = io_msg )

                IF ( io_stat /= 0 ) THEN
                  err_msg = 'Error opening '//TRIM(atm_failure_file)//' - '//TRIM(io_msg)
                  CALL Display_Message( ROUTINE_NAME, err_msg, WARNING )
                ELSE
                  ! Update file status
                  atm_file_status = 'OLD'
                  ! Report failure
                  n_failed = COUNT(.NOT. (dtb_delta < FWDTL_TOLERANCE(ialpha)))
                  dtb_maxloc = MAXLOC(dtb_delta, DIM=1)
                  WRITE(fid,'(40("="))')
                  WRITE(fid,'("*** ATM: dtb_delta, FWDTL_TOLERANCE test ", &
                           &"failed for profile #",i0," layer ",i0)') m, k
                  WRITE(fid,'("Number of failures: ",i0," of ",i0)') n_failed, n_channels
                  WRITE(fid,'("alpha = ",es9.2)') ALPHA(ialpha)
                  WRITE(fid,'("Values for largest magnitude failure:")')
                  WRITE(fid,'("tb_NLp                ",f19.15)') tb_NLp(dtb_maxloc)
                  WRITE(fid,'("tb_NLm              - ",f19.15)') tb_NLm(dtb_maxloc)
                  WRITE(fid,'("                      ",19("-"))')
                  WRITE(fid,'("dtb_NL              = ",f19.15,2x,"(",es22.15,")")') dtb_NL(dtb_maxloc)     , dtb_NL(dtb_maxloc)
                  WRITE(fid,'("tb_TL               = ",f19.15,2x,"(",es22.15,")")') tb_TL(dtb_maxloc)      , tb_TL(dtb_maxloc)
                  WRITE(fid,'("***--->>> delta     = ",f19.15,2x,"(",es22.15,")")') dtb_delta(dtb_maxloc)  , dtb_delta(dtb_maxloc)
                  WRITE(fid,'("***--->>> threshold = ",f19.15,2x,"(",es22.15,")")') FWDTL_TOLERANCE(ialpha), FWDTL_TOLERANCE(ialpha)
                  WRITE(fid,'("*** ATM: dtb_delta, FWDTL_TOLERANCE test ", &
                           &"failed for profile #",i0," layer ",i0)') m, k
                  WRITE(fid,'(40("="),/)')
                  CLOSE(fid)
                END IF
              END IF

            END DO layer_loop

            CALL UnitTest_Report(utest)

          END DO atm_component_loop


          ! ===============================
          ! =====     CLOUD TESTS     =====
          ! ===============================

          ! Initialise the surface TL data (only has to be done once here)
          sfc_TL = sfc_FWD(m)
          CALL CRTM_Surface_Zero(sfc_TL)


          ! Loop over the clouds
          cloud_loop: DO ic = 1, atm_FWD(m)%n_Clouds


            ! Loop over the components to test
            cloud_component_loop: DO icomponent = 1, N_CLOUD_COMPONENTS


              ! Initialise test
              WRITE(utest_msg,'("CLOUD FWD/TL test | ",&
                               &"Profile: ",i0," | ",&
                               &"Cloud: ",i0," | ",&
                               &"Component: ",a," | ",&
                               &"Alpha: ",es9.2," | ",&
                               &"Algorithm: ",a)') &
                    m, ic, TRIM(CLOUD_COMPONENT_NAME(icomponent)), &
                    ALPHA(ialpha), &
                    TRIM(RT_ALGORITHM_NAME(i))
              CALL UnitTest_Setup(utest,TRIM(utest_msg))


              ! Loop over the atmospheric layers
              cloud_layer_loop: DO k = CLOUD_LAYER_BEGIN, atm_FWD(m)%n_Layers, CLOUD_LAYER_STEP


                ! No cloud in this layer
                IF ( atm_FWD(m)%Cloud(ic)%Effective_Radius(k) == ZERO .OR. &
                     atm_FWD(m)%Cloud(ic)%Water_Content(k)    == ZERO ) CYCLE cloud_layer_loop


                WRITE(*, '(4x,"- Running tangent- and non-linear model for cloud ",i0,&
                             &", layer ",i0," perturbation...")') ic, k


                ! Initialise the atm TL and NL data
                atm_NLm(1) = atm_FWD(m)
                atm_NLp(1) = atm_FWD(m)
                atm_TL(1)  = atm_FWD(m); CALL CRTM_Atmosphere_Zero(atm_TL)


                ! Select cloud component to test
                SELECT CASE(TRIM(CLOUD_COMPONENT_NAME(icomponent)))

                  CASE('Effective Radius')
                    delta = ALPHA(ialpha) * DX
                    atm_NLm(1)%Cloud(ic)%Effective_Radius(k) = atm_FWD(m)%Cloud(ic)%Effective_Radius(k) - (delta/TWO)
                    IF( atm_NLm(1)%Cloud(ic)%Effective_Radius(k) < ZERO ) THEN
                      delta = atm_FWD(m)%Cloud(ic)%Effective_Radius(k) * 0.3_fp
                      atm_NLm(1)%Cloud(ic)%Effective_Radius(k) = atm_FWD(m)%Cloud(ic)%Effective_Radius(k) - (delta/TWO)
                    END IF
                    atm_TL(1)%Cloud(ic)%Effective_Radius(k)  = delta
                    atm_NLp(1)%Cloud(ic)%Effective_Radius(k) = atm_FWD(m)%Cloud(ic)%Effective_Radius(k) + (delta/TWO)
                  CASE('Water Content')
                    delta = ALPHA(ialpha) * DX * atm_FWD(m)%Cloud(ic)%Water_Content(k)
                    atm_NLm(1)%Cloud(ic)%Water_Content(k) = atm_FWD(m)%Cloud(ic)%Water_Content(k) - (delta/TWO)
                    IF( atm_NLm(1)%Cloud(ic)%Water_Content(k) < ZERO ) THEN
                      delta = atm_FWD(m)%Cloud(ic)%Water_Content(k) * 0.3_fp
                      atm_NLm(1)%Cloud(ic)%Water_Content(k) = atm_FWD(m)%Cloud(ic)%Water_Content(k) - (delta/TWO)
                    END IF                    
                    atm_TL(1)%Cloud(ic)%Water_Content(k)  = delta
                    atm_NLp(1)%Cloud(ic)%Water_Content(k) = atm_FWD(m)%Cloud(ic)%Water_Content(k) + (delta/TWO)


                  CASE DEFAULT
                    err_msg = 'Unrecognised CLOUD component name: '//TRIM(CLOUD_COMPONENT_NAME(icomponent))
                    CALL Display_Message(ROUTINE_NAME, err_msg, FAILURE ); STOP

                END SELECT


                ! Perform tangent-linear calculations
                err_stat = CRTM_Tangent_Linear( atm_FWD(m:m)      , &
                                                sfc_FWD(m:m)      , &
                                                [atm_TL]          , &
                                                [sfc_TL]          , &
                                                geo(m:m)          , &
                                                chinfo            , &
                                                rts_FWD(:,m:m)    , &
                                                rts_TL(:,m:m)     , &
                                                Options = opt(m:m)  )
                IF ( err_stat /= SUCCESS ) THEN
                  err_msg = 'Error when calling CRTM_Tangent_Linear'
                  CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
                END IF


                ! Perform non-linear calculations
                ! ...Negative perturbation
                err_stat = CRTM_Forward( atm_NLm    , &
                                         sfc_FWD(m:m) , &
                                         geo(m:m)     , &
                                         chinfo       , &
                                         rts_NLm(:,m:m), &
                                         Options = opt(m:m) )
                IF ( err_stat /= SUCCESS ) THEN
                  err_msg = 'Error when calling CRTM_Forward for negative perturbation'
                  CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
                END IF
                ! ...Positive perturbation
                err_stat = CRTM_Forward( atm_NLp    , &
                                         sfc_FWD(m:m) , &
                                         geo(m:m)     , &
                                         chinfo       , &
                                         rts_NLp(:,m:m), &
                                         Options = opt(m:m) )
                IF ( err_stat /= SUCCESS ) THEN
                  err_msg = 'Error when calling CRTM_Forward for positive perturbation'
                  CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
                END IF
                ! Compute the test quantities
           IF( IR_MW_Sensor ) THEN
              tb_NLm = rts_NLm(:,m)%Brightness_Temperature
              tb_NLp = rts_NLp(:,m)%Brightness_Temperature
              tb_TL  = rts_TL(:,m)%Brightness_Temperature
           ELSE
             IF( n_Stokes == 1 ) THEN
               tb_NLm = rts_NLm(:,m)%Radiance
               tb_NLp = rts_NLp(:,m)%Radiance
               tb_TL  = rts_TL(:,m)%Radiance
             ELSE
               tb_NLm = rts_NLm(:,m)%Stokes(iStoke)
               tb_NLp = rts_NLp(:,m)%Stokes(iStoke)
               tb_TL  = rts_TL(:,m)%Stokes(iStoke)
             END IF
           END IF

              dtb_NL = tb_NLp - tb_NLm
              dtb_delta = ABS(dtb_NL - tb_TL)


                ! Apply the test
                CALL UnitTest_Assert( utest, ALL(dtb_delta < FWDTL_TOLERANCE(ialpha)) )


                ! Output info for failed test
                IF ( UnitTest_Failed(utest) ) THEN
                  ! Open failure report file
                  fid = Get_Lun()
                  OPEN( fid, FILE     = cloud_failure_file, &
                             FORM     = 'FORMATTED', &
                             STATUS   = cloud_file_status, &
                             POSITION = 'APPEND', &
                             IOSTAT   = io_stat, &
                             IOMSG    = io_msg )
                  IF ( io_stat /= 0 ) THEN
                    err_msg = 'Error opening '//TRIM(cloud_failure_file)//' - '//TRIM(io_msg)
                    CALL Display_Message( ROUTINE_NAME, err_msg, WARNING )
                  ELSE
                    ! Update file status
                    cloud_file_status = 'OLD'
                    ! Report failure
                    n_failed = COUNT(.NOT. (dtb_delta < FWDTL_TOLERANCE(ialpha)))
                    dtb_maxloc = MAXLOC(dtb_delta, DIM=1)
                    WRITE(fid,'(40("="))')
                    WRITE(fid,'("*** CLOUD: dtb_delta, FWDTL_TOLERANCE test ", &
                             &"failed for profile #",i0," cloud #",i0," and layer ",i0)') m, ic, k
                    WRITE(fid,'("Number of failures: ",i0," of ",i0)') n_failed, n_channels
                    WRITE(fid,'("alpha = ",es9.2)') ALPHA(ialpha)
                    WRITE(fid,'("Values for largest magnitude failure:")')
                    WRITE(fid,'("tb_NLp                ",f19.15)') tb_NLp(dtb_maxloc)
                    WRITE(fid,'("tb_NLm              - ",f19.15)') tb_NLm(dtb_maxloc)
                    WRITE(fid,'("                      ",19("-"))')
                    WRITE(fid,'("dtb_NL              = ",f19.15,2x,"(",es22.15,")")') dtb_NL(dtb_maxloc)   , dtb_NL(dtb_maxloc)
                    WRITE(fid,'("tb_TL               = ",f19.15,2x,"(",es22.15,")")') tb_TL(dtb_maxloc)    , tb_TL(dtb_maxloc)
                    WRITE(fid,'("***--->>> delta     = ",f19.15,2x,"(",es22.15,")")') dtb_delta(dtb_maxloc), dtb_delta(dtb_maxloc)
                    WRITE(fid,'("***--->>> threshold = ",f19.15,2x,"(",es22.15,")")') FWDTL_TOLERANCE(ialpha), &
                                                                                      FWDTL_TOLERANCE(ialpha)
                    WRITE(fid,'("*** CLOUD: dtb_delta, FWDTL_TOLERANCE test ", &
                             &"failed for profile #",i0," layer ",i0)') m, k
                    WRITE(fid,'(40("="),/)')
                    CLOSE(fid)

!                  IF( m > 0 ) THEN
!                    print *,' qliu -02 '
!                    STOP
!                  END IF

                  END IF
                END IF

              END DO cloud_layer_loop

              CALL UnitTest_Report(utest)

            END DO cloud_component_loop

          END DO cloud_loop


          ! =================================
          ! =====     AEROSOL TESTS     =====
          ! =================================

          ! Initialise the surface TL data (only has to be done once here)
          sfc_TL = sfc_FWD(m)
          CALL CRTM_Surface_Zero(sfc_TL)


          ! Loop over the aerosols
          aerosol_loop: DO ia = 1, atm_FWD(m)%n_Aerosols


            ! Loop over the components to test
            aerosol_component_loop: DO icomponent = 1, N_AEROSOL_COMPONENTS


              ! Initialise test
              WRITE(utest_msg,'("AEROSOL FWD/TL test | ",&
                               &"Profile: ",i0," | ",&
                               &"Aerosol: ",i0," | ",&
                               &"Component: ",a," | ",&
                               &"Alpha: ",es9.2," | ",&
                               &"Algorithm: ",a)') &
                    m, ic, TRIM(AEROSOL_COMPONENT_NAME(icomponent)), &
                    ALPHA(ialpha), &
                    TRIM(RT_ALGORITHM_NAME(i))
              CALL UnitTest_Setup(utest,TRIM(utest_msg))


              ! Loop over the atmospheric layers
              aerosol_layer_loop: DO k = AEROSOL_LAYER_BEGIN, atm_FWD(m)%n_Layers, AEROSOL_LAYER_STEP


                ! No aerosol in this layer
                IF ( atm_FWD(m)%Aerosol(ia)%Effective_Radius(k) == ZERO .OR. &
                     atm_FWD(m)%Aerosol(ia)%Concentration(k)    == ZERO ) CYCLE aerosol_layer_loop


                WRITE(*, '(4x,"- Running tangent- and non-linear model for aerosol ",i0,&
                             &", layer ",i0," perturbation...")') ic, k


                ! Initialise the atm TL and NL data
                atm_NLm(1) = atm_FWD(m)
                atm_NLp(1) = atm_FWD(m)
                atm_TL(1)  = atm_FWD(m); CALL CRTM_Atmosphere_Zero(atm_TL)


                ! Select aerosol component to test
                SELECT CASE(TRIM(AEROSOL_COMPONENT_NAME(icomponent)))

                  CASE('Effective Radius')
                    delta = ALPHA(ialpha) * DX
                    atm_TL(1)%Aerosol(ia)%Effective_Radius(k)  = delta
                    atm_NLp(1)%Aerosol(ia)%Effective_Radius(k) = atm_FWD(m)%Aerosol(ia)%Effective_Radius(k) + (delta/TWO)
                    atm_NLm(1)%Aerosol(ia)%Effective_Radius(k) = atm_FWD(m)%Aerosol(ia)%Effective_Radius(k) - (delta/TWO)

                  CASE('Concentration')
                    delta = ALPHA(ialpha) * DX * atm_FWD(m)%Aerosol(ia)%Concentration(k)
                    atm_TL(1)%Aerosol(ia)%Concentration(k)  = delta
                    atm_NLp(1)%Aerosol(ia)%Concentration(k) = atm_FWD(m)%Aerosol(ia)%Concentration(k) + (delta/TWO)
                    atm_NLm(1)%Aerosol(ia)%Concentration(k) = atm_FWD(m)%Aerosol(ia)%Concentration(k) - (delta/TWO)
!           print *,' aerosol con ',k,delta,atm_NLm(1)%Aerosol(ia)%Concentration(k), &
!               atm_NLp(1)%Aerosol(ia)%Concentration(k)
                  CASE DEFAULT
                    err_msg = 'Unrecognised AEROSOL component name: '//TRIM(AEROSOL_COMPONENT_NAME(icomponent))
                    CALL Display_Message(ROUTINE_NAME, err_msg, FAILURE ); STOP

                END SELECT


                ! Perform tangent-linear calculations
                err_stat = CRTM_Tangent_Linear( atm_FWD(m:m)      , &
                                                sfc_FWD(m:m)      , &
                                                [atm_TL]          , &
                                                [sfc_TL]          , &
                                                geo(m:m)          , &
                                                chinfo            , &
                                                rts_FWD(:,m:m)    , &
                                                rts_TL(:,m:m)     , &
                                                Options = opt(m:m)  )
                IF ( err_stat /= SUCCESS ) THEN
                  err_msg = 'Error when calling CRTM_Tangent_Linear'
                  CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
                END IF


                ! Perform non-linear calculations
                ! ...Negative perturbation
                err_stat = CRTM_Forward( atm_NLm    , &
                                         sfc_FWD(m:m) , &
                                         geo(m:m)     , &
                                         chinfo       , &
                                         rts_NLm(:,m:m), &
                                         Options = opt(m:m) )
                IF ( err_stat /= SUCCESS ) THEN
                  err_msg = 'Error when calling CRTM_Forward for negative perturbation'
                  CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
                END IF
                ! ...Positive perturbation
                err_stat = CRTM_Forward( atm_NLp    , &
                                         sfc_FWD(m:m) , &
                                         geo(m:m)     , &
                                         chinfo       , &
                                         rts_NLp(:,m:m), &
                                         Options = opt(m:m) )
                IF ( err_stat /= SUCCESS ) THEN
                  err_msg = 'Error when calling CRTM_Forward for positive perturbation'
                  CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
                END IF


                ! Compute the test quantities
           IF( IR_MW_Sensor ) THEN
              tb_NLm = rts_NLm(:,m)%Brightness_Temperature
              tb_NLp = rts_NLp(:,m)%Brightness_Temperature
              tb_TL  = rts_TL(:,m)%Brightness_Temperature
           ELSE
             IF( n_Stokes == 1 ) THEN
               tb_NLm = rts_NLm(:,m)%Radiance
               tb_NLp = rts_NLp(:,m)%Radiance
               tb_TL  = rts_TL(:,m)%Radiance
             ELSE
               tb_NLm = rts_NLm(:,m)%Stokes(iStoke)
               tb_NLp = rts_NLp(:,m)%Stokes(iStoke)
               tb_TL  = rts_TL(:,m)%Stokes(iStoke)
             END IF
           END IF
              dtb_NL = tb_NLp - tb_NLm
              dtb_delta = ABS(dtb_NL - tb_TL)


                ! Apply the test
                CALL UnitTest_Assert( utest, ALL(dtb_delta < FWDTL_TOLERANCE(ialpha)) )


                ! Output info for failed test
                IF ( UnitTest_Failed(utest) ) THEN
                  ! Open failure report file
                  fid = Get_Lun()
                  OPEN( fid, FILE     = aerosol_failure_file, &
                             FORM     = 'FORMATTED', &
                             STATUS   = aerosol_file_status, &
                             POSITION = 'APPEND', &
                             IOSTAT   = io_stat, &
                             IOMSG    = io_msg )
                  IF ( io_stat /= 0 ) THEN
                    err_msg = 'Error opening '//TRIM(aerosol_failure_file)//' - '//TRIM(io_msg)
                    CALL Display_Message( ROUTINE_NAME, err_msg, WARNING )
                  ELSE
                    ! Update file status
                    aerosol_file_status = 'OLD'
                    ! Report failure
                    n_failed = COUNT(.NOT. (dtb_delta < FWDTL_TOLERANCE(ialpha)))
                    dtb_maxloc = MAXLOC(dtb_delta, DIM=1)
                    WRITE(fid,'(40("="))')
                    WRITE(fid,'("*** AEROSOL: dtb_delta, FWDTL_TOLERANCE test ", &
                             &"failed for profile #",i0," aerosol #",i0," and layer ",i0)') m, ic, k
                    WRITE(fid,'("Number of failures: ",i0," of ",i0)') n_failed, n_channels
                    WRITE(fid,'("alpha = ",es9.2)') ALPHA(ialpha)
                    WRITE(fid,'("Values for largest magnitude failure:")')
                    WRITE(fid,'("tb_NLp                ",f19.15)') tb_NLp(dtb_maxloc)
                    WRITE(fid,'("tb_NLm              - ",f19.15)') tb_NLm(dtb_maxloc)
                    WRITE(fid,'("                      ",19("-"))')
                    WRITE(fid,'("dtb_NL              = ",f19.15,2x,"(",es22.15,")")') dtb_NL(dtb_maxloc)   , dtb_NL(dtb_maxloc)
                    WRITE(fid,'("tb_TL               = ",f19.15,2x,"(",es22.15,")")') tb_TL(dtb_maxloc)    , tb_TL(dtb_maxloc)
                    WRITE(fid,'("***--->>> delta     = ",f19.15,2x,"(",es22.15,")")') dtb_delta(dtb_maxloc), dtb_delta(dtb_maxloc)
                    WRITE(fid,'("***--->>> threshold = ",f19.15,2x,"(",es22.15,")")') FWDTL_TOLERANCE(ialpha), &
                                                                                      FWDTL_TOLERANCE(ialpha)
                    WRITE(fid,'("*** AEROSOL: dtb_delta, FWDTL_TOLERANCE test ", &
                             &"failed for profile #",i0," layer ",i0)') m, k
                    WRITE(fid,'(40("="),/)')
!                    STOP
                    
                    
                    
                    CLOSE(fid)
                  END IF
                END IF

              END DO aerosol_layer_loop

              CALL UnitTest_Report(utest)

            END DO aerosol_component_loop

          END DO aerosol_loop


          ! ======================================
          ! =====     LAND SURFACE TESTS     =====
          ! ======================================

          ! Initialise the atmosphere TL data (only has to be done once here)
          atm_TL = atm_FWD(m)
          CALL CRTM_Atmosphere_Zero(atm_TL)


          ! Loop over the components to test
          land_sfc_component_loop: DO icomponent = 1, N_LAND_SFC_COMPONENTS


            ! No land, so no need to loop
            IF ( .NOT. (sfc_FWD(m)%Land_Coverage > ZERO) ) EXIT land_sfc_component_loop


            ! Initialise test
            WRITE(utest_msg,'("LAND SFC FWD/TL test | ",&
                             &"Profile: ",i0," | ",&
                             &"Component: ",a," | ",&
                             &"Alpha: ",es9.2," | ",&
                             &"Algorithm: ",a)') &
                  m, TRIM(LAND_SFC_COMPONENT_NAME(icomponent)), &
                  ALPHA(ialpha), &
                  TRIM(RT_ALGORITHM_NAME(i))
            CALL UnitTest_Setup(utest,TRIM(utest_msg))


            ! Initialise the sfc TL and NL data
            sfc_NLm = sfc_FWD(m)
            sfc_NLp = sfc_FWD(m)
            sfc_TL  = sfc_FWD(m); CALL CRTM_Surface_Zero(sfc_TL)


            ! Select surface component to test
            WRITE(*, '(4x,"- Running tangent- and non-linear model for land surface perturbation...")')
            SELECT CASE(TRIM(LAND_SFC_COMPONENT_NAME(icomponent)))

              CASE('Land Temperature')
                delta = ALPHA(ialpha) * DX
                sfc_TL%Land_Temperature  = delta
                sfc_NLp%Land_Temperature = sfc_FWD(m)%Land_Temperature + (delta/TWO)
                sfc_NLm%Land_Temperature = sfc_FWD(m)%Land_Temperature - (delta/TWO)

              CASE('Soil Moisture Content')
                delta = ALPHA(ialpha) * DX * sfc_FWD(m)%Soil_Moisture_Content
                sfc_TL%Soil_Moisture_Content  = delta
                sfc_NLp%Soil_Moisture_Content = sfc_FWD(m)%Soil_Moisture_Content + (delta/TWO)
                sfc_NLm%Soil_Moisture_Content = sfc_FWD(m)%Soil_Moisture_Content - (delta/TWO)

              CASE('Canopy Water Content')
                delta = ALPHA(ialpha) * DX * sfc_FWD(m)%Canopy_Water_Content
                sfc_TL%Canopy_Water_Content  = delta
                sfc_NLp%Canopy_Water_Content = sfc_FWD(m)%Canopy_Water_Content + (delta/TWO)
                sfc_NLm%Canopy_Water_Content = sfc_FWD(m)%Canopy_Water_Content - (delta/TWO)

              CASE('Vegetation Fraction')
                delta = ALPHA(ialpha) * DX * sfc_FWD(m)%Vegetation_Fraction
                sfc_TL%Vegetation_Fraction  = delta
                sfc_NLp%Vegetation_Fraction = sfc_FWD(m)%Vegetation_Fraction + (delta/TWO)
                sfc_NLm%Vegetation_Fraction = sfc_FWD(m)%Vegetation_Fraction - (delta/TWO)

              CASE('Soil Temperature')
                delta = ALPHA(ialpha) * DX * sfc_FWD(m)%Soil_Temperature
                sfc_TL%Soil_Temperature  = delta
                sfc_NLp%Soil_Temperature = sfc_FWD(m)%Soil_Temperature + (delta/TWO)
                sfc_NLm%Soil_Temperature = sfc_FWD(m)%Soil_Temperature - (delta/TWO)

              CASE('LAI')
                delta = ALPHA(ialpha) * DX * sfc_FWD(m)%LAI
                sfc_TL%LAI  = delta
                sfc_NLp%LAI = sfc_FWD(m)%LAI + (delta/TWO)
                sfc_NLm%LAI = sfc_FWD(m)%LAI - (delta/TWO)

              CASE DEFAULT
                err_msg = 'Unrecognised LAND SFC component name: '//TRIM(LAND_SFC_COMPONENT_NAME(icomponent))
                CALL Display_Message(ROUTINE_NAME, err_msg, FAILURE ); STOP

            END SELECT


            ! Perform tangent-linear calculations
            err_stat = CRTM_Tangent_Linear( atm_FWD(m:m)      , &
                                            sfc_FWD(m:m)      , &
                                            [atm_TL]          , &
                                            [sfc_TL]          , &
                                            geo(m:m)          , &
                                            chinfo            , &
                                            rts_FWD(:,m:m)    , &
                                            rts_TL(:,m:m)     , &
                                            Options = opt(m:m)  )
            IF ( err_stat /= SUCCESS ) THEN
              err_msg = 'Error when calling CRTM_Tangent_Linear'
              CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
            END IF


            ! Perform non-linear calculations
            ! ...Negative perturbation
            err_stat = CRTM_Forward( atm_FWD(m:m)      , &
                                     [sfc_NLm]         , &
                                     geo(m:m)          , &
                                     chinfo            , &
                                     rts_NLm(:,m:m)    , &
                                     Options = opt(m:m)  )
            IF ( err_stat /= SUCCESS ) THEN
              err_msg = 'Error when calling CRTM_Forward for negative perturbation'
              CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
            END IF
            ! ...Positive perturbation
            err_stat = CRTM_Forward( atm_FWD(m:m)      , &
                                     [sfc_NLp]         , &
                                     geo(m:m)          , &
                                     chinfo            , &
                                     rts_NLp(:,m:m)    , &
                                     Options = opt(m:m)  )
            IF ( err_stat /= SUCCESS ) THEN
              err_msg = 'Error when calling CRTM_Forward for positive perturbation'
              CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
            END IF


            ! Compute the test quantities
           IF( IR_MW_Sensor ) THEN
              tb_NLm = rts_NLm(:,m)%Brightness_Temperature
              tb_NLp = rts_NLp(:,m)%Brightness_Temperature
              tb_TL  = rts_TL(:,m)%Brightness_Temperature
           ELSE
             IF( n_Stokes == 1 ) THEN
               tb_NLm = rts_NLm(:,m)%Radiance
               tb_NLp = rts_NLp(:,m)%Radiance
               tb_TL  = rts_TL(:,m)%Radiance
             ELSE
               tb_NLm = rts_NLm(:,m)%Stokes(iStoke)
               tb_NLp = rts_NLp(:,m)%Stokes(iStoke)
               tb_TL  = rts_TL(:,m)%Stokes(iStoke)
             END IF
           END IF
              dtb_NL = tb_NLp - tb_NLm
              dtb_delta = ABS(dtb_NL - tb_TL)

            ! Apply the test
            CALL UnitTest_Assert( utest, ALL(dtb_delta < FWDTL_TOLERANCE(ialpha)) )


            ! Output info for failed test
            IF ( UnitTest_Failed(utest) ) THEN
              ! Open failure report file
              fid = Get_Lun()
              OPEN( fid, FILE     = landsfc_failure_file, &
                         FORM     = 'FORMATTED', &
                         STATUS   = landsfc_file_status, &
                         POSITION = 'APPEND', &
                         IOSTAT   = io_stat, &
                         IOMSG    = io_msg )
              IF ( io_stat /= 0 ) THEN
                err_msg = 'Error opening '//TRIM(landsfc_failure_file)//' - '//TRIM(io_msg)
                CALL Display_Message( ROUTINE_NAME, err_msg, WARNING )
              ELSE
                ! Update file status
                landsfc_file_status = 'OLD'
                ! Report failure
                n_failed = COUNT(.NOT. (dtb_delta < FWDTL_TOLERANCE(ialpha)))
                dtb_maxloc = MAXLOC(dtb_delta, DIM=1)
                WRITE(fid,'(40("="))')
                WRITE(fid,'("*** LAND SFC: dtb_delta, FWDTL_TOLERANCE test ", &
                         &"failed for profile #",i0)') m
                WRITE(fid,'("Number of failures: ",i0," of ",i0)') n_failed, n_channels
                WRITE(fid,'("alpha = ",es9.2)') ALPHA(ialpha)
                WRITE(fid,'("Input perturbation = ",es13.6)') delta
                WRITE(fid,'("Values for largest magnitude failure:")')
                WRITE(fid,'("tb_NLp                ",f19.15)') tb_NLp(dtb_maxloc)
                WRITE(fid,'("tb_NLm              - ",f19.15)') tb_NLm(dtb_maxloc)
                WRITE(fid,'("                      ",19("-"))')
                WRITE(fid,'("dtb_NL              = ",f19.15,2x,"(",es22.15,")")') dtb_NL(dtb_maxloc)     , dtb_NL(dtb_maxloc)
                WRITE(fid,'("tb_TL               = ",f19.15,2x,"(",es22.15,")")') tb_TL(dtb_maxloc)      , tb_TL(dtb_maxloc)
                WRITE(fid,'("***--->>> delta     = ",f19.15,2x,"(",es22.15,")")') dtb_delta(dtb_maxloc)  , dtb_delta(dtb_maxloc)
                WRITE(fid,'("***--->>> threshold = ",f19.15,2x,"(",es22.15,")")') FWDTL_TOLERANCE(ialpha), FWDTL_TOLERANCE(ialpha)
                CALL CRTM_Surface_Inspect(sfc_FWD(m), Unit=fid)
                WRITE(fid,'("*** LAND SFC: dtb_delta, FWDTL_TOLERANCE test ", &
                         &"failed for profile #",i0)') m
                WRITE(fid,'(40("="),/)')
                CLOSE(fid)
              END IF
            END IF

            CALL UnitTest_Report(utest)

          END DO land_sfc_component_loop


          ! =======================================
          ! =====     WATER SURFACE TESTS     =====
          ! =======================================

          ! Initialise the atmosphere TL data (only has to be done once here)
          atm_TL = atm_FWD(m)
          CALL CRTM_Atmosphere_Zero(atm_TL)


          ! Loop over the components to test
          water_sfc_component_loop: DO icomponent = 1, N_WATER_SFC_COMPONENTS


            ! No water, so no need to loop
            IF ( .NOT. (sfc_FWD(m)%Water_Coverage > ZERO) ) EXIT water_sfc_component_loop


            ! Initialise test
            WRITE(utest_msg,'("WATER SFC FWD/TL test | ",&
                             &"Profile: ",i0," | ",&
                             &"Component: ",a," | ",&
                             &"Alpha: ",es9.2," | ",&
                             &"Algorithm: ",a)') &
                  m, TRIM(WATER_SFC_COMPONENT_NAME(icomponent)), &
                  ALPHA(ialpha), &
                  TRIM(RT_ALGORITHM_NAME(i))
            CALL UnitTest_Setup(utest,TRIM(utest_msg))


            ! Initialise the sfc TL and NL data
            sfc_NLm = sfc_FWD(m)
            sfc_NLp = sfc_FWD(m)
            sfc_TL  = sfc_FWD(m); CALL CRTM_Surface_Zero(sfc_TL)


            ! Select surface component to test
            WRITE(*, '(4x,"- Running tangent- and non-linear model for water surface perturbation...")')
            SELECT CASE(TRIM(WATER_SFC_COMPONENT_NAME(icomponent)))

              CASE('Water Temperature')
                delta = ALPHA(ialpha) * DX
                sfc_TL%Water_Temperature  = delta
                sfc_NLp%Water_Temperature = sfc_FWD(m)%Water_Temperature + (delta/TWO)
                sfc_NLm%Water_Temperature = sfc_FWD(m)%Water_Temperature - (delta/TWO)

              CASE('Wind Speed')
                delta = ALPHA(ialpha) * DX * sfc_FWD(m)%Wind_Speed
                sfc_TL%Wind_Speed  = delta
                sfc_NLp%Wind_Speed = sfc_FWD(m)%Wind_Speed + (delta/TWO)
                sfc_NLm%Wind_Speed = sfc_FWD(m)%Wind_Speed - (delta/TWO)

              CASE('Wind Direction')
                delta = ALPHA(ialpha) * DX * sfc_FWD(m)%Wind_Direction
                sfc_TL%Wind_Direction  = delta
                sfc_NLp%Wind_Direction = sfc_FWD(m)%Wind_Direction + (delta/TWO)
                sfc_NLm%Wind_Direction = sfc_FWD(m)%Wind_Direction - (delta/TWO)

              CASE('Salinity')
                delta = ALPHA(ialpha) * DX * sfc_FWD(m)%Salinity
                sfc_TL%Salinity  = delta
                sfc_NLp%Salinity = sfc_FWD(m)%Salinity + (delta/TWO)
                sfc_NLm%Salinity = sfc_FWD(m)%Salinity - (delta/TWO)

              CASE DEFAULT
                err_msg = 'Unrecognised WATER SFC component name: '//TRIM(WATER_SFC_COMPONENT_NAME(icomponent))
                CALL Display_Message(ROUTINE_NAME, err_msg, FAILURE ); STOP

            END SELECT


            ! Perform tangent-linear calculations
            err_stat = CRTM_Tangent_Linear( atm_FWD(m:m)      , &
                                            sfc_FWD(m:m)      , &
                                            [atm_TL]          , &
                                            [sfc_TL]          , &
                                            geo(m:m)          , &
                                            chinfo            , &
                                            rts_FWD(:,m:m)    , &
                                            rts_TL(:,m:m)     , &
                                            Options = opt(m:m)  )
            IF ( err_stat /= SUCCESS ) THEN
              err_msg = 'Error when calling CRTM_Tangent_Linear'
              CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
            END IF


            ! Perform non-linear calculations
            ! ...Negative perturbation
            err_stat = CRTM_Forward( atm_FWD(m:m)      , &
                                     [sfc_NLm]         , &
                                     geo(m:m)          , &
                                     chinfo            , &
                                     rts_NLm(:,m:m)    , &
                                     Options = opt(m:m)  )
            IF ( err_stat /= SUCCESS ) THEN
              err_msg = 'Error when calling CRTM_Forward for negative perturbation'
              CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
            END IF
            ! ...Positive perturbation
            err_stat = CRTM_Forward( atm_FWD(m:m)      , &
                                     [sfc_NLp]         , &
                                     geo(m:m)          , &
                                     chinfo            , &
                                     rts_NLp(:,m:m)    , &
                                     Options = opt(m:m)  )
            IF ( err_stat /= SUCCESS ) THEN
              err_msg = 'Error when calling CRTM_Forward for positive perturbation'
              CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
            END IF


            ! Compute the test quantities
           IF( IR_MW_Sensor ) THEN
              tb_NLm = rts_NLm(:,m)%Brightness_Temperature
              tb_NLp = rts_NLp(:,m)%Brightness_Temperature
              tb_TL  = rts_TL(:,m)%Brightness_Temperature
           ELSE
             IF( n_Stokes == 1 ) THEN
               tb_NLm = rts_NLm(:,m)%Radiance
               tb_NLp = rts_NLp(:,m)%Radiance
               tb_TL  = rts_TL(:,m)%Radiance
             ELSE
               tb_NLm = rts_NLm(:,m)%Stokes(iStoke)
               tb_NLp = rts_NLp(:,m)%Stokes(iStoke)
               tb_TL  = rts_TL(:,m)%Stokes(iStoke)
             END IF
           END IF
              dtb_NL = tb_NLp - tb_NLm
              dtb_delta = ABS(dtb_NL - tb_TL)


            ! Apply the test
            CALL UnitTest_Assert( utest, ALL(dtb_delta < FWDTL_TOLERANCE(ialpha)) )


            ! Output info for failed test
            IF ( UnitTest_Failed(utest) ) THEN
              ! Open failure report file
              fid = Get_Lun()
              OPEN( fid, FILE     = watersfc_failure_file, &
                         FORM     = 'FORMATTED', &
                         STATUS   = watersfc_file_status, &
                         POSITION = 'APPEND', &
                         IOSTAT   = io_stat, &
                         IOMSG    = io_msg )
              IF ( io_stat /= 0 ) THEN
                err_msg = 'Error opening '//TRIM(watersfc_failure_file)//' - '//TRIM(io_msg)
                CALL Display_Message( ROUTINE_NAME, err_msg, WARNING )
              ELSE
                ! Update file status
                watersfc_file_status = 'OLD'
                ! Report failure
                n_failed = COUNT(.NOT. (dtb_delta < FWDTL_TOLERANCE(ialpha)))
                dtb_maxloc = MAXLOC(dtb_delta, DIM=1)
                WRITE(fid,'(40("="))')
                WRITE(fid,'("*** WATER SFC: dtb_delta, FWDTL_TOLERANCE test ", &
                         &"failed for profile #",i0)') m
                WRITE(fid,'("Number of failures: ",i0," of ",i0)') n_failed, n_channels
                WRITE(fid,'("alpha = ",es9.2)') ALPHA(ialpha)
                WRITE(fid,'("Input perturbation = ",es13.6)') delta
                WRITE(fid,'("Values for largest magnitude failure:")')
                WRITE(fid,'("tb_NLp                ",f19.15)') tb_NLp(dtb_maxloc)
                WRITE(fid,'("tb_NLm              - ",f19.15)') tb_NLm(dtb_maxloc)
                WRITE(fid,'("                      ",19("-"))')
                WRITE(fid,'("dtb_NL              = ",f19.15,2x,"(",es22.15,")")') dtb_NL(dtb_maxloc)     , dtb_NL(dtb_maxloc)
                WRITE(fid,'("tb_TL               = ",f19.15,2x,"(",es22.15,")")') tb_TL(dtb_maxloc)      , tb_TL(dtb_maxloc)
                WRITE(fid,'("***--->>> delta     = ",f19.15,2x,"(",es22.15,")")') dtb_delta(dtb_maxloc)  , dtb_delta(dtb_maxloc)
                WRITE(fid,'("***--->>> threshold = ",f19.15,2x,"(",es22.15,")")') FWDTL_TOLERANCE(ialpha), FWDTL_TOLERANCE(ialpha)
                CALL CRTM_Surface_Inspect(sfc_FWD(m), Unit=fid)
                WRITE(fid,'("*** WATER SFC: dtb_delta, FWDTL_TOLERANCE test ", &
                         &"failed for profile #",i0)') m
                WRITE(fid,'(40("="),/)')
                CLOSE(fid)
              END IF
            END IF

            CALL UnitTest_Report(utest)

          END DO water_sfc_component_loop


          ! ======================================
          ! =====     SNOW SURFACE TESTS     =====
          ! ======================================

          ! Initialise the atmosphere TL data (only has to be done once here)
          atm_TL = atm_FWD(m)
          CALL CRTM_Atmosphere_Zero(atm_TL)


          ! Loop over the components to test
!          snow_sfc_component_loop: DO icomponent = 1, N_SNOW_SFC_COMPONENTS
          snow_sfc_component_loop: DO icomponent = 1, 1   ! No TL/AD for Snow now N_SNOW_SFC_COMPONENTS

            ! No snow, so no need to loop
            IF ( .NOT. (sfc_FWD(m)%Snow_Coverage > ZERO) ) EXIT snow_sfc_component_loop


            ! Initialise test
            WRITE(utest_msg,'("SNOW SFC FWD/TL test | ",&
                             &"Profile: ",i0," | ",&
                             &"Component: ",a," | ",&
                             &"Alpha: ",es9.2," | ",&
                             &"Algorithm: ",a)') &
                  m, TRIM(SNOW_SFC_COMPONENT_NAME(icomponent)), &
                  ALPHA(ialpha), &
                  TRIM(RT_ALGORITHM_NAME(i))
            CALL UnitTest_Setup(utest,TRIM(utest_msg))


            ! Initialise the sfc TL and NL data
            sfc_NLm = sfc_FWD(m)
            sfc_NLp = sfc_FWD(m)
            sfc_TL  = sfc_FWD(m); CALL CRTM_Surface_Zero(sfc_TL)


            ! Select surface component to test
            WRITE(*, '(4x,"- Running tangent- and non-linear model for snow surface perturbation...")')
               print *,' name = ',TRIM(SNOW_SFC_COMPONENT_NAME(icomponent))
            SELECT CASE(TRIM(SNOW_SFC_COMPONENT_NAME(icomponent)))
              CASE('Snow_Temperature')
                delta = ALPHA(ialpha) * DX
                sfc_TL%Snow_Temperature  = delta
                sfc_NLp%Snow_Temperature = sfc_FWD(m)%Snow_Temperature + (delta/TWO)
                sfc_NLm%Snow_Temperature = sfc_FWD(m)%Snow_Temperature - (delta/TWO)
 
              CASE('Snow_Depth')
                delta = ALPHA(ialpha) * DX * sfc_FWD(m)%Snow_Depth
                sfc_TL%Snow_Depth  = delta
                sfc_NLp%Snow_Depth = sfc_FWD(m)%Snow_Depth + (delta/TWO)
                sfc_NLm%Snow_Depth = sfc_FWD(m)%Snow_Depth - (delta/TWO)
            print *,' snow_depth ',delta,sfc_NLm(1)%Snow_Depth,sfc_NLp(1)%Snow_Depth
             

              CASE('Snow_Density')
                delta = ALPHA(ialpha) * DX * sfc_FWD(m)%Snow_Density
                sfc_TL%Snow_Density  = delta
                sfc_NLp%Snow_Density = sfc_FWD(m)%Snow_Density + (delta/TWO)
                sfc_NLm%Snow_Density = sfc_FWD(m)%Snow_Density - (delta/TWO)

              CASE('Snow_Grain_Size')
                delta = ALPHA(ialpha) * DX * sfc_FWD(m)%Snow_Grain_Size
                sfc_TL%Snow_Grain_Size  = delta
                sfc_NLp%Snow_Grain_Size = sfc_FWD(m)%Snow_Grain_Size + (delta/TWO)
                sfc_NLm%Snow_Grain_Size = sfc_FWD(m)%Snow_Grain_Size - (delta/TWO)

              CASE DEFAULT
                err_msg = 'Unrecognised SNOW SFC component name: '//TRIM(SNOW_SFC_COMPONENT_NAME(icomponent))
                CALL Display_Message(ROUTINE_NAME, err_msg, FAILURE ); STOP

            END SELECT


            ! Perform tangent-linear calculations
            err_stat = CRTM_Tangent_Linear( atm_FWD(m:m)      , &
                                            sfc_FWD(m:m)      , &
                                            [atm_TL]          , &
                                            [sfc_TL]          , &
                                            geo(m:m)          , &
                                            chinfo            , &
                                            rts_FWD(:,m:m)    , &
                                            rts_TL(:,m:m)     , &
                                            Options = opt(m:m)  )
            IF ( err_stat /= SUCCESS ) THEN
              err_msg = 'Error when calling CRTM_Tangent_Linear'
              CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
            END IF


            ! Perform non-linear calculations
            ! ...Negative perturbation
            err_stat = CRTM_Forward( atm_FWD(m:m)      , &
                                     [sfc_NLm]         , &
                                     geo(m:m)          , &
                                     chinfo            , &
                                     rts_NLm(:,m:m)    , &
                                     Options = opt(m:m)  )
            IF ( err_stat /= SUCCESS ) THEN
              err_msg = 'Error when calling CRTM_Forward for negative perturbation'
              CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
            END IF
            ! ...Positive perturbation
            err_stat = CRTM_Forward( atm_FWD(m:m)      , &
                                     [sfc_NLp]         , &
                                     geo(m:m)          , &
                                     chinfo            , &
                                     rts_NLp(:,m:m)    , &
                                     Options = opt(m:m)  )
            IF ( err_stat /= SUCCESS ) THEN
              err_msg = 'Error when calling CRTM_Forward for positive perturbation'
              CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
            END IF


            ! Compute the test quantities
           IF( IR_MW_Sensor ) THEN
              tb_NLm = rts_NLm(:,m)%Brightness_Temperature
              tb_NLp = rts_NLp(:,m)%Brightness_Temperature
              tb_TL  = rts_TL(:,m)%Brightness_Temperature
           ELSE
             IF( n_Stokes == 1 ) THEN
               tb_NLm = rts_NLm(:,m)%Radiance
               tb_NLp = rts_NLp(:,m)%Radiance
               tb_TL  = rts_TL(:,m)%Radiance
             ELSE
               tb_NLm = rts_NLm(:,m)%Stokes(iStoke)
               tb_NLp = rts_NLp(:,m)%Stokes(iStoke)
               tb_TL  = rts_TL(:,m)%Stokes(iStoke)
             END IF
           END IF
              dtb_NL = tb_NLp - tb_NLm
              dtb_delta = ABS(dtb_NL - tb_TL)


            ! Apply the test
            CALL UnitTest_Assert( utest, ALL(dtb_delta < FWDTL_TOLERANCE(ialpha)) )


            ! Output info for failed test
            IF ( UnitTest_Failed(utest) ) THEN
              ! Open failure report file
              fid = Get_Lun() 
              OPEN( fid, FILE     = snowsfc_failure_file, &
                         FORM     = 'FORMATTED', &
                         STATUS   = snowsfc_file_status, &
                         POSITION = 'APPEND', &
                         IOSTAT   = io_stat, &
                         IOMSG    = io_msg )
              IF ( io_stat /= 0 ) THEN
                err_msg = 'Error opening '//TRIM(snowsfc_failure_file)//' - '//TRIM(io_msg)
                CALL Display_Message( ROUTINE_NAME, err_msg, WARNING )
              ELSE
                ! Update file status
                snowsfc_file_status = 'OLD'
                ! Report failure
                n_failed = COUNT(.NOT. (dtb_delta < FWDTL_TOLERANCE(ialpha)))
                dtb_maxloc = MAXLOC(dtb_delta, DIM=1)
                WRITE(fid,'(40("="))')
                WRITE(fid,'("*** SNOW SFC: dtb_delta, FWDTL_TOLERANCE test ", &
                         &"failed for profile #",i0)') m
                WRITE(fid,'("Number of failures: ",i0," of ",i0)') n_failed, n_channels
                WRITE(fid,'("alpha = ",es9.2)') ALPHA(ialpha)
                WRITE(fid,'("Input perturbation = ",es13.6)') delta
                WRITE(fid,'("Values for largest magnitude failure:")')
                WRITE(fid,'("tb_NLp                ",f19.15)') tb_NLp(dtb_maxloc)
                WRITE(fid,'("tb_NLm              - ",f19.15)') tb_NLm(dtb_maxloc)
                WRITE(fid,'("                      ",19("-"))')
                WRITE(fid,'("dtb_NL              = ",f19.15,2x,"(",es22.15,")")') dtb_NL(dtb_maxloc)     , dtb_NL(dtb_maxloc)
                WRITE(fid,'("tb_TL               = ",f19.15,2x,"(",es22.15,")")') tb_TL(dtb_maxloc)      , tb_TL(dtb_maxloc)
                WRITE(fid,'("***--->>> delta     = ",f19.15,2x,"(",es22.15,")")') dtb_delta(dtb_maxloc)  , dtb_delta(dtb_maxloc)
                WRITE(fid,'("***--->>> threshold = ",f19.15,2x,"(",es22.15,")")') FWDTL_TOLERANCE(ialpha), FWDTL_TOLERANCE(ialpha)
                CALL CRTM_Surface_Inspect(sfc_FWD(m), Unit=fid)
                WRITE(fid,'("*** SNOW SFC: dtb_delta, FWDTL_TOLERANCE test ", &
                         &"failed for profile #",i0)') m
                WRITE(fid,'(40("="),/)')
                CLOSE(fid)
              END IF
            END IF

            CALL UnitTest_Report(utest)

          END DO snow_sfc_component_loop


          ! ======================================
          ! =====     ICE SURFACE TESTS     =====
          ! ======================================

          ! Initialise the atmosphere TL data (only has to be done once here)
          atm_TL = atm_FWD(m)
          CALL CRTM_Atmosphere_Zero(atm_TL)


          ! Loop over the components to test
          ice_sfc_component_loop: DO icomponent = 1, N_ICE_SFC_COMPONENTS

            ! No ice, so no need to loop
            IF ( .NOT. (sfc_FWD(m)%Ice_Coverage > ZERO) ) EXIT ice_sfc_component_loop


            ! Initialise test
            WRITE(utest_msg,'("ICE SFC FWD/TL test | ",&
                             &"Profile: ",i0," | ",&
                             &"Component: ",a," | ",&
                             &"Alpha: ",es9.2," | ",&
                             &"Algorithm: ",a)') &
                  m, TRIM(ICE_SFC_COMPONENT_NAME(icomponent)), &
                  ALPHA(ialpha), &
                  TRIM(RT_ALGORITHM_NAME(i))
            CALL UnitTest_Setup(utest,TRIM(utest_msg))


            ! Initialise the sfc TL and NL data
            sfc_NLm = sfc_FWD(m)
            sfc_NLp = sfc_FWD(m)
            sfc_TL  = sfc_FWD(m); CALL CRTM_Surface_Zero(sfc_TL)


            ! Select surface component to test
            WRITE(*, '(4x,"- Running tangent- and non-linear model for ice surface perturbation...")')
            SELECT CASE(TRIM(ICE_SFC_COMPONENT_NAME(icomponent)))

              CASE('Ice_Temperature')
                delta = ALPHA(ialpha) * DX
                sfc_TL%Ice_Temperature  = delta
                sfc_NLp%Ice_Temperature = sfc_FWD(m)%Ice_Temperature + (delta/TWO)
                sfc_NLm%Ice_Temperature = sfc_FWD(m)%Ice_Temperature - (delta/TWO)

              CASE('Ice_Thickness')
                delta = ALPHA(ialpha) * DX * sfc_FWD(m)%Ice_Thickness
                sfc_TL%Ice_Thickness  = delta
                sfc_NLp%Ice_Thickness = sfc_FWD(m)%Ice_Thickness + (delta/TWO)
                sfc_NLm%Ice_Thickness = sfc_FWD(m)%Ice_Thickness - (delta/TWO)

              CASE('Ice_Density')
                delta = ALPHA(ialpha) * DX * sfc_FWD(m)%Ice_Density
                sfc_TL%Ice_Density  = delta
                sfc_NLp%Ice_Density = sfc_FWD(m)%Ice_Density + (delta/TWO)
                sfc_NLm%Ice_Density = sfc_FWD(m)%Ice_Density - (delta/TWO)

              CASE('Ice_Roughness')
                delta = ALPHA(ialpha) * DX * sfc_FWD(m)%Ice_Roughness
                sfc_TL%Ice_Roughness  = delta
                sfc_NLp%Ice_Roughness = sfc_FWD(m)%Ice_Roughness + (delta/TWO)
                sfc_NLm%Ice_Roughness = sfc_FWD(m)%Ice_Roughness - (delta/TWO)

              CASE DEFAULT
                err_msg = 'Unrecognised ICE SFC component name: '//TRIM(ICE_SFC_COMPONENT_NAME(icomponent))
                CALL Display_Message(ROUTINE_NAME, err_msg, FAILURE ); STOP

            END SELECT


            ! Perform tangent-linear calculations
            err_stat = CRTM_Tangent_Linear( atm_FWD(m:m)      , &
                                            sfc_FWD(m:m)      , &
                                            [atm_TL]          , &
                                            [sfc_TL]          , &
                                            geo(m:m)          , &
                                            chinfo            , &
                                            rts_FWD(:,m:m)    , &
                                            rts_TL(:,m:m)     , &
                                            Options = opt(m:m)  )
            IF ( err_stat /= SUCCESS ) THEN
              err_msg = 'Error when calling CRTM_Tangent_Linear'
              CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
            END IF


            ! Perform non-linear calculations
            ! ...Negative perturbation
            err_stat = CRTM_Forward( atm_FWD(m:m)      , &
                                     [sfc_NLm]         , &
                                     geo(m:m)          , &
                                     chinfo            , &
                                     rts_NLm(:,m:m)    , &
                                     Options = opt(m:m)  )
            IF ( err_stat /= SUCCESS ) THEN
              err_msg = 'Error when calling CRTM_Forward for negative perturbation'
              CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
            END IF
            ! ...Positive perturbation
            err_stat = CRTM_Forward( atm_FWD(m:m)      , &
                                     [sfc_NLp]         , &
                                     geo(m:m)          , &
                                     chinfo            , &
                                     rts_NLp(:,m:m)    , &
                                     Options = opt(m:m)  )
            IF ( err_stat /= SUCCESS ) THEN
              err_msg = 'Error when calling CRTM_Forward for positive perturbation'
              CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
            END IF


            ! Compute the test quantities
           IF( IR_MW_Sensor ) THEN
              tb_NLm = rts_NLm(:,m)%Brightness_Temperature
              tb_NLp = rts_NLp(:,m)%Brightness_Temperature
              tb_TL  = rts_TL(:,m)%Brightness_Temperature
           ELSE
             IF( n_Stokes == 1 ) THEN
               tb_NLm = rts_NLm(:,m)%Radiance
               tb_NLp = rts_NLp(:,m)%Radiance
               tb_TL  = rts_TL(:,m)%Radiance
             ELSE
               tb_NLm = rts_NLm(:,m)%Stokes(iStoke)
               tb_NLp = rts_NLp(:,m)%Stokes(iStoke)
               tb_TL  = rts_TL(:,m)%Stokes(iStoke)
             END IF
           END IF
              dtb_NL = tb_NLp - tb_NLm
              dtb_delta = ABS(dtb_NL - tb_TL)


            ! Apply the test
            CALL UnitTest_Assert( utest, ALL(dtb_delta < FWDTL_TOLERANCE(ialpha)) )

            ! Output info for failed test
            IF ( UnitTest_Failed(utest) ) THEN
              ! Open failure report file
              fid = Get_Lun()
              OPEN( fid, FILE     = icesfc_failure_file, &
                         FORM     = 'FORMATTED', &
                         STATUS   = icesfc_file_status, &
                         POSITION = 'APPEND', &
                         IOSTAT   = io_stat, &
                         IOMSG    = io_msg )
              IF ( io_stat /= 0 ) THEN
                err_msg = 'Error opening '//TRIM(icesfc_failure_file)//' - '//TRIM(io_msg)
                CALL Display_Message( ROUTINE_NAME, err_msg, WARNING )
              ELSE
                ! Update file status
                icesfc_file_status = 'OLD'
                ! Report failure
                n_failed = COUNT(.NOT. (dtb_delta < FWDTL_TOLERANCE(ialpha)))
                dtb_maxloc = MAXLOC(dtb_delta, DIM=1)
                WRITE(fid,'(40("="))')
                WRITE(fid,'("*** ICE SFC: dtb_delta, FWDTL_TOLERANCE test ", &
                         &"failed for profile #",i0)') m
                WRITE(fid,'("Number of failures: ",i0," of ",i0)') n_failed, n_channels
                WRITE(fid,'("alpha = ",es9.2)') ALPHA(ialpha)
                WRITE(fid,'("Input perturbation = ",es13.6)') delta
                WRITE(fid,'("Values for largest magnitude failure:")')
                WRITE(fid,'("tb_NLp                ",f19.15)') tb_NLp(dtb_maxloc)
                WRITE(fid,'("tb_NLm              - ",f19.15)') tb_NLm(dtb_maxloc)
                WRITE(fid,'("                      ",19("-"))')
                WRITE(fid,'("dtb_NL              = ",f19.15,2x,"(",es22.15,")")') dtb_NL(dtb_maxloc)     , dtb_NL(dtb_maxloc)
                WRITE(fid,'("tb_TL               = ",f19.15,2x,"(",es22.15,")")') tb_TL(dtb_maxloc)      , tb_TL(dtb_maxloc)
                WRITE(fid,'("***--->>> delta     = ",f19.15,2x,"(",es22.15,")")') dtb_delta(dtb_maxloc)  , dtb_delta(dtb_maxloc)
                WRITE(fid,'("***--->>> threshold = ",f19.15,2x,"(",es22.15,")")') FWDTL_TOLERANCE(ialpha), FWDTL_TOLERANCE(ialpha)
                CALL CRTM_Surface_Inspect(sfc_FWD(m), Unit=fid)
                WRITE(fid,'("*** ICE SFC: dtb_delta, FWDTL_TOLERANCE test ", &
                         &"failed for profile #",i0)') m
                WRITE(fid,'(40("="),/)')
                CLOSE(fid)
              END IF
            END IF

            CALL UnitTest_Report(utest)

          END DO ice_sfc_component_loop

        END DO profile_loop

      END DO alpha_loop

    END DO rt_algorithm_loop


    ! Destroy the data structures
    CALL CRTM_Atmosphere_Destroy(atm_FWD)
    CALL CRTM_Atmosphere_Destroy(atm_NLm)
    CALL CRTM_Atmosphere_Destroy(atm_NLp)
    CALL CRTM_Atmosphere_Destroy(atm_TL)
    CALL CRTM_Surface_Destroy(sfc_FWD)
    CALL CRTM_Surface_Destroy(sfc_NLm)
    CALL CRTM_Surface_Destroy(sfc_NLp)
    CALL CRTM_Surface_Destroy(sfc_TL)
    CALL CRTM_RTSolution_Destroy(rts_Base)
    CALL CRTM_RTSolution_Destroy(rts_NLm)
    CALL CRTM_RTSolution_Destroy(rts_NLp)
    CALL CRTM_RTSolution_Destroy(rts_TL)
    CALL CRTM_Options_Destroy(opt)


    ! Deallocate the structure arrays
    DEALLOCATE( rts_Base , &
                rts_NLm  , &
                rts_NLp  , &
                rts_TL   , &
                tb_NLm   , &
                tb_NLp   , &
                dtb_NL   , &
                tb_TL    , &
                dtb_delta, &
                atm_FWD  , &
                sfc_FWD  , &
                STAT = alloc_stat )


  END SUBROUTINE Test_CRTM_FWDTL


  ! ==================================
  ! Tangent-linear/adjoint consistency
  ! ==================================

  SUBROUTINE Test_CRTM_TLAD(utest, atm, sfc, geo, chinfo)
    ! Arguments
    TYPE(UnitTest_type)        , INTENT(IN OUT) :: utest
    TYPE(CRTM_Atmosphere_type) , INTENT(IN)     :: atm(:)
    TYPE(CRTM_Surface_type)    , INTENT(IN)     :: sfc(:)
    TYPE(CRTM_Geometry_type)   , INTENT(IN)     :: geo(:)
    TYPE(CRTM_ChannelInfo_type), INTENT(IN)     :: chinfo(:)
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Test_CRTM_TLAD'
    ! Local variable
    CHARACTER(256) :: err_msg, alloc_msg, io_msg, utest_msg
    CHARACTER(256) :: test_failure_file
    CHARACTER(7) :: file_status
    INTEGER :: fid
    INTEGER :: err_stat, alloc_stat, io_stat
    INTEGER :: i, m, nc
    INTEGER :: ii, ic, k, n, kk 
    INTEGER :: n_channels, n_profiles, n_sensors
    INTEGER :: n_aerosols, n_clouds, n_absorbers, n_layers 
    REAL(fp) :: TLtTL, dxtAD
    REAL(fp) :: delta
    REAL(fp) :: TLAD_tolerance
    TYPE(CRTM_RTSolution_type),  ALLOCATABLE :: rts_Base(:,:,:)
    TYPE(CRTM_RTSolution_type),  ALLOCATABLE :: rts_TL(:,:,:)
    TYPE(CRTM_RTSolution_type),  ALLOCATABLE :: rts_AD(:,:,:)
    TYPE(CRTM_RTSolution_type),  ALLOCATABLE :: rts_AD_in(:,:,:)
    TYPE(CRTM_Atmosphere_type),  ALLOCATABLE :: atm_FWD(:)
    TYPE(CRTM_Atmosphere_type),  ALLOCATABLE :: atm_TL(:)
    TYPE(CRTM_Atmosphere_type),  ALLOCATABLE :: atm_AD(:)
    TYPE(CRTM_Surface_type),     ALLOCATABLE :: sfc_FWD(:)
    TYPE(CRTM_Surface_type),     ALLOCATABLE :: sfc_TL(:)
    TYPE(CRTM_Surface_type),     ALLOCATABLE :: sfc_AD(:)
    TYPE(Timing_type) :: timing


    ! Get test dimensions
    n_channels = SUM(CRTM_ChannelInfo_n_Channels(chinfo))
    n_profiles = SIZE(atm)
    n_sensors  = SIZE(chinfo)
    n_layers   = SIZE(atm(1)%Temperature(:)) 
    
    ! Perform all the allocations
    ALLOCATE( rts_Base (n_channels, n_profiles, N_RT_ALGORITHMS), &
              rts_TL   (n_channels, n_profiles, N_RT_ALGORITHMS), &
              rts_AD   (n_channels, n_profiles, N_RT_ALGORITHMS), &
              rts_AD_in(n_channels, n_profiles, N_RT_ALGORITHMS), & 
              atm_FWD  (n_Profiles), &
              atm_TL   (n_Profiles), &
              atm_AD   (n_Profiles), &
              sfc_FWD  (n_Profiles), &
              sfc_TL   (n_Profiles), &
              sfc_AD   (n_Profiles), &
              STAT   = alloc_stat, &
              ERRMSG = alloc_msg   )
    IF ( alloc_stat /= 0 ) THEN
      err_msg = 'Error allocating data structure arrays - '//TRIM(alloc_msg)
      CALL Display_Message(ROUTINE_NAME, err_msg, FAILURE ); STOP
    END IF

    call CRTM_RTSolution_Create(rts_Base, n_layers) 
    call CRTM_RTSolution_Create(rts_TL, n_layers) 
    call CRTM_RTSolution_Create(rts_AD, n_layers) 
    call CRTM_RTSolution_Create(rts_AD_in, n_layers) 

 
    ! Loop over types of radiative transfer algorithms
    rt_algorithm_loop: DO i = 1, N_RT_ALGORITHMS               
      opt(:)%RT_Algorithm_Id = RT_ALGORITHM_ID(i)
    CALL CRTM_RTSolution_Zero( rts_Base(:,:,i) )
    CALL CRTM_RTSolution_Zero( rts_TL(:,:,i) )
    CALL CRTM_RTSolution_Zero( rts_AD(:,:,i) )   
    
      ! Setup for test failure reporting
      file_status = 'REPLACE'
      test_failure_file = TRIM(RT_ALGORITHM_NAME(i))//'.TLAD.test_failure_report'
      
  !   opt(:)%Overlap_Id = CloudCover_Maximum_Overlap()      
  !   opt(:)%Overlap_Id = CloudCover_Random_Overlap()      
  !   opt(:)%Overlap_Id = CloudCover_MaxRan_Overlap()      
  !   opt(:)%Overlap_Id = CloudCover_Overcast_Overlap()      

      ! Output info
      ! ...Algorithm identifier
      WRITE(*,'(30("*"),1x,"TL/AD Comparisons for RT Algorithm ",a,1x,30("*"))') TRIM(RT_ALGORITHM_NAME(i))
      ! ...Sensors to process
      WRITE(*, '(4x,"- Sensors: ",99(a,:))') chinfo%sensor_id


      ! Initialise adjoint structures
      atm_AD = atm
      CALL CRTM_Atmosphere_Zero( atm_AD )
      sfc_AD = sfc
      CALL CRTM_Surface_Zero( sfc_AD )


      ! Copy the inputs so they can be modified if necessary
      atm_FWD = atm
      sfc_FWD = sfc


      ! For emission algorithm calculations
      IF ( i == EMISSION_INDEX ) THEN
        atm_FWD%n_Clouds   = 0
        atm_FWD%n_Aerosols = 0
      END IF

      ! Perform baseline calculations for
      ! SOI, ADA and Emission RT algorithms
      WRITE(*, '(4x,"- Running forward model...")')
      CALL Timing_Begin(timing)
      err_stat = CRTM_Forward( atm_FWD        , &
                               sfc_FWD        , &
                               geo            , &
                               chinfo         , &
                               rts_Base(:,:,i), &
                               Options = opt(n1:n2)    )

      CALL Timing_End(timing)
      IF ( err_stat /= SUCCESS ) THEN
        err_msg = 'Error when calling CRTM_Forward'
        CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
      END IF
      CALL Timing_Display(timing)

      ! Assign perturbation for tangent-linear calculation

    CALL CRTM_RTSolution_Zero( rts_Base(:,:,i) )
    CALL CRTM_RTSolution_Zero( rts_TL(:,:,i) )
      CALL CRTM_Atmosphere_ZERO(atm_TL)   
      CALL CRTM_Surface_ZERO(sfc_TL)
      
      CALL Assign_TL_Atmosphere( DX, atm_FWD, atm_TL )
      CALL Assign_TL_Surface(    DX, sfc_FWD, sfc_TL )

      ! Perform TL calculations around the FWD baseline calc
      WRITE(*, '(4x,"- Running tangent-linear model...")')

 
      CALL Timing_Begin(timing)
      err_stat = CRTM_Tangent_Linear( atm_FWD        , &
                                      sfc_FWD        , &
                                      atm_TL         , &
                                      sfc_TL         , &
                                      geo            , &
                                      chinfo         , &
                                      rts_Base(:,:,i), &
                                      rts_TL(:,:,i)  , &
                                      Options = opt(n1:n2)    )
            print *,' TL err_stat = ',err_stat
   
      CALL Timing_End(timing)
      IF ( err_stat /= SUCCESS ) THEN
        err_msg = 'Error when calling CRTM_Tangent_Linear'
        CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
      END IF
      CALL Timing_Display(timing)

      ! Initialise adjoint structure
    CALL CRTM_RTSolution_Zero( rts_AD(:,:,i) )
    CALL CRTM_RTSolution_Zero( rts_Base(:,:,i) )
      CALL CRTM_Atmosphere_ZERO(atm_AD)   
      CALL CRTM_Surface_ZERO(sfc_AD)
          
      IF( IR_MW_Sensor ) THEN
        rts_AD(:,:,i)%Brightness_Temperature    = rts_TL(:,:,i)%Brightness_Temperature   
        rts_AD_in(:,:,i)%Brightness_Temperature = rts_TL(:,:,i)%Brightness_Temperature
        rts_AD(:,:,i)%Radiance = ZERO
      ELSE
        IF( n_Stokes == 1 ) THEN
          rts_AD(:,:,i)%Radiance    = rts_TL(:,:,i)%Radiance
          rts_AD_in(:,:,i)%Radiance = rts_TL(:,:,i)%Radiance
          rts_AD(:,:,i)%Brightness_Temperature = ZERO
        ELSE
          rts_AD(:,:,i)%Stokes(iStoke)    = rts_TL(:,:,i)%Stokes(iStoke)
          rts_AD_in(:,:,i)%Stokes(iStoke) = rts_TL(:,:,i)%Stokes(iStoke)
          rts_AD(:,:,i)%Brightness_Temperature = ZERO
        END IF
      END IF

      
      ! Perform AD calculations
      WRITE(*, '(4x,"- Running adjoint model...")')
      CALL Timing_Begin(timing)
      err_stat = CRTM_Adjoint( atm_FWD        , &
                               sfc_FWD        , &
                               rts_AD(:,:,i)  , &
                               geo            , &
                               chinfo         , &
                               atm_AD         , &
                               sfc_AD         , &
                               rts_Base(:,:,i), &
                               Options = opt(n1:n2)    )
            print *,' AD err_stat = ',err_stat
      CALL Timing_End(timing)
      IF ( err_stat /= SUCCESS ) THEN
        err_msg = 'Error when calling CRTM_Adjoint'
        CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
      END IF
      CALL Timing_Display(timing)

      ! Initialise test
      WRITE(utest_msg,'("TL/AD comparison test | ",&
                       &"Algorithm: ",a)') &
            TRIM(RT_ALGORITHM_NAME(i))
      CALL UnitTest_Setup(utest,TRIM(utest_msg))


      ! Perform test on each profile separately
      profile_loop: DO m = 1, n_profiles

        ! Output some info
        IF( n_profiles > 10 ) THEN
        
        IF ( MOD(m,n_profiles/10) == 0 ) &
          WRITE(*, '(6x,"Testing profile #",i0," of ",i0,"...")') m, n_profiles
          
        END IF
        ! Compute the TL test quantity
                            
      IF( IR_MW_Sensor ) THEN
        TLtTL = SUM( rts_TL(:,m,i)%Brightness_Temperature * rts_AD_in(:,m,i)%Brightness_Temperature )
      ELSE
        IF( n_Stokes == 1 ) THEN
        TLtTL = SUM( rts_TL(:,m,i)%Radiance * rts_AD_in(:,m,i)%Radiance )
        ELSE
        TLtTL = SUM( rts_TL(:,m,i)%Stokes(iStoke) * rts_AD_in(:,m,i)%Stokes(iStoke) )
        END IF
      END IF

 
      print *,' TLDL ',TLtTL
        ! Compute the AD test quantity
        ! ...Atmospheric part
        dxtAD = SUM(atm_TL(m)%Level_Pressure * atm_AD(m)%Level_Pressure) + &
                SUM(atm_TL(m)%Pressure       * atm_AD(m)%Pressure      ) + &
                SUM(atm_TL(m)%Temperature    * atm_AD(m)%Temperature   ) + &
                SUM(atm_TL(m)%Cloud_Fraction * atm_AD(m)%Cloud_Fraction) + &          
                SUM(atm_TL(m)%Absorber       * atm_AD(m)%Absorber      )

        ! ...Cloud part
        IF ( i /= EMISSION_INDEX ) THEN
          DO nc = 1, atm_FWD(m)%n_Clouds
            dxtAD = dxtAD + &
                    SUM(atm_TL(m)%Cloud(nc)%Water_Content    * atm_AD(m)%Cloud(nc)%Water_Content   ) +  &
                    SUM(atm_TL(m)%Cloud(nc)%Effective_Radius * atm_AD(m)%Cloud(nc)%Effective_Radius)
          END DO
        END IF
        ! ...Aerosol part
        IF ( i /= EMISSION_INDEX ) THEN
          DO nc = 1, atm_FWD(m)%n_Aerosols
            dxtAD = dxtAD + &
                    SUM(atm_TL(m)%Aerosol(nc)%Concentration    * atm_AD(m)%Aerosol(nc)%Concentration   ) +  &
                    SUM(atm_TL(m)%Aerosol(nc)%Effective_Radius * atm_AD(m)%Aerosol(nc)%Effective_Radius)
          END DO
        END IF
        ! ...sfc part
        IF ( sfc_FWD(m)%Land_Coverage > ZERO ) THEN
          dxtAD = dxtAD + &
                  ( sfc_TL(m)%Land_Temperature      * sfc_AD(m)%Land_Temperature      ) + &
                  ( sfc_TL(m)%Soil_Moisture_Content * sfc_AD(m)%Soil_Moisture_Content ) + &
                  ( sfc_TL(m)%Canopy_Water_Content  * sfc_AD(m)%Canopy_Water_Content  ) + &
                  ( sfc_TL(m)%Vegetation_Fraction   * sfc_AD(m)%Vegetation_Fraction   ) + &
                  ( sfc_TL(m)%Soil_Temperature      * sfc_AD(m)%Soil_Temperature      )
        END IF
        IF ( sfc_FWD(m)%Water_Coverage > ZERO ) THEN
          dxtAD = dxtAD + &
                  ( sfc_TL(m)%Water_Temperature     * sfc_AD(m)%Water_Temperature     ) + &
                  ( sfc_TL(m)%Wind_Speed            * sfc_AD(m)%Wind_Speed            ) + &
                  ( sfc_TL(m)%Wind_Direction        * sfc_AD(m)%Wind_Direction        ) + &
                  ( sfc_TL(m)%Salinity              * sfc_AD(m)%Salinity              )
        END IF
        IF ( sfc_FWD(m)%Snow_Coverage > ZERO ) THEN
          dxtAD = dxtAD + &
                  ( sfc_TL(m)%Snow_Temperature      * sfc_AD(m)%Snow_Temperature      ) + &
                  ( sfc_TL(m)%Snow_Depth            * sfc_AD(m)%Snow_Depth            ) + &
                  ( sfc_TL(m)%Snow_Density          * sfc_AD(m)%Snow_Density          ) + &
                  ( sfc_TL(m)%Snow_Grain_Size       * sfc_AD(m)%Snow_Grain_Size       )
        END IF
        IF ( sfc_FWD(m)%Ice_Coverage > ZERO ) THEN
          dxtAD = dxtAD + &
                  ( sfc_TL(m)%Ice_Temperature       * sfc_AD(m)%Ice_Temperature       ) + &
                  ( sfc_TL(m)%Ice_Thickness         * sfc_AD(m)%Ice_Thickness         ) + &
                  ( sfc_TL(m)%Ice_Density           * sfc_AD(m)%Ice_Density           ) + &
                  ( sfc_TL(m)%Ice_Roughness         * sfc_AD(m)%Ice_Roughness         )
        END IF

       print *,' dxtAD ',dxtAD,TLtTL,dxtAD-TLtTL
        ! Compute the test quantities
        TLAD_tolerance = SPACING( MAX(ABS(TLtTL),ABS(dxtAD)) ) * TLAD_ULP
        
        !  qliu, 01/19/2022
        TLAD_tolerance = 2.0_fp * TLAD_tolerance
        
        delta = ABS(TLtTL - dxtAD)


        ! Apply the test
        CALL UnitTest_Assert( utest, delta < TLAD_tolerance )
        write(6,'(ES25.18,5x,ES25.18,5x,ES25.18,5x,ES25.18)') TLtTL, dxtAD, delta, TLAD_tolerance


        ! Output info for failed test
        IF ( UnitTest_Failed(utest) ) THEN
          ! Open failure report file
          fid = Get_Lun()
          OPEN( fid, FILE     = test_failure_file, &
                     FORM     = 'FORMATTED', &
                     STATUS   = file_status, &
                     POSITION = 'APPEND', &
                     IOSTAT   = io_stat, &
                     IOMSG    = io_msg )
          IF ( io_stat /= 0 ) THEN
            err_msg = 'Error opening '//TRIM(test_failure_file)//' - '//TRIM(io_msg)
            CALL Display_Message( ROUTINE_NAME, err_msg, WARNING )
          ELSE
            ! Update file status
            file_status = 'OLD'
            ! Report failure
            WRITE(fid,'(40("="))')
            WRITE(fid,'("*** TLtTL, dxtAD equality test failed for profile #",i0)') m
            WRITE(fid,'("TLtTL                = ",es22.15)') TLtTL
            WRITE(fid,'("dxtAD                = ",es22.15)') dxtAD
            WRITE(fid,'("***--->>> delta      = ",es22.15)') delta
            WRITE(fid,'("***--->>> threshold  = ",es22.15)') TLAD_tolerance
!            CALL CRTM_Atmosphere_Inspect(atm_FWD(m), Unit=fid)
!            CALL CRTM_Surface_Inspect(sfc_FWD(m), Unit=fid)
            WRITE(fid,'("*** TLtTL, dxtAD equality test failed for profile #",i0)') m
            WRITE(fid,'(40("="),/)')
            CLOSE(fid)
            
          END IF
        END IF

      END DO profile_loop

      CALL UnitTest_Report(utest)

    END DO rt_algorithm_loop


    ! Destroy the data structures
    CALL CRTM_Atmosphere_Destroy(atm_FWD)
    CALL CRTM_Atmosphere_Destroy(atm_TL)
    CALL CRTM_Atmosphere_Destroy(atm_AD)
    CALL CRTM_Surface_Destroy(sfc_FWD)
    CALL CRTM_Surface_Destroy(sfc_TL)
    CALL CRTM_Surface_Destroy(sfc_AD)
    CALL CRTM_RTSolution_Destroy(rts_Base)
    CALL CRTM_RTSolution_Destroy(rts_TL)
    CALL CRTM_RTSolution_Destroy(rts_AD)
    CALL CRTM_Options_Destroy(opt)


    ! Deallocate the structure arrays
    DEALLOCATE( rts_Base, &
                rts_AD  , &
                rts_TL  , &
                atm_FWD , &
                atm_TL  , &
                atm_AD  , &
                sfc_FWD , &
                sfc_TL  , &
                sfc_AD  , &
                STAT = alloc_stat )

  END SUBROUTINE Test_CRTM_TLAD


  ! ==================================
  ! K-matrix/adjoint consistency
  ! ==================================
  SUBROUTINE Test_CRTM_ADK(utest, atm, sfc, geo, chinfo)
    ! Arguments
    TYPE(UnitTest_type)        , INTENT(IN OUT) :: utest
    TYPE(CRTM_Atmosphere_type) , INTENT(IN)     :: atm(:)
    TYPE(CRTM_Surface_type)    , INTENT(IN)     :: sfc(:)
    TYPE(CRTM_Geometry_type)   , INTENT(IN)     :: geo(:)
    TYPE(CRTM_ChannelInfo_type), INTENT(IN)     :: chinfo(:)
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Test_CRTM_ADK'
    CHARACTER(256) :: atm_file_Ref, sfc_file_Ref
    CHARACTER(256) :: atm_failure_file, sfc_failure_file
    CHARACTER(7) :: atm_file_status, sfc_file_status
    ! Local variable
    CHARACTER(256) :: err_msg, alloc_msg, io_msg, utest_msg
    CHARACTER(256) :: test_failure_file
    CHARACTER(7) :: file_status
    INTEGER :: fid
    INTEGER :: err_stat, alloc_stat, io_stat
    INTEGER :: i, m, nc
    INTEGER :: ii, ic, k, n, kk 
    INTEGER :: n_channels, n_profiles, n_sensors
    INTEGER :: n_aerosols, n_clouds, n_absorbers, n_layers 
    REAL(fp) :: TLtTL, dxtAD, Tolerance_x
    REAL(fp) :: delta
    REAL(fp) :: TLAD_tolerance
    TYPE(CRTM_RTSolution_type),  ALLOCATABLE :: rts_Base(:,:,:)
    TYPE(CRTM_RTSolution_type),  ALLOCATABLE :: rts_K(:,:,:)
    TYPE(CRTM_RTSolution_type),  ALLOCATABLE :: rts_AD(:,:,:)
    TYPE(CRTM_RTSolution_type),  ALLOCATABLE :: rts_TL(:,:,:)
    TYPE(CRTM_Atmosphere_type),  ALLOCATABLE :: atm_FWD(:)
    TYPE(CRTM_Atmosphere_type),  ALLOCATABLE :: atm_TL(:)
    TYPE(CRTM_Atmosphere_type),  ALLOCATABLE :: atm_K(:,:)
    TYPE(CRTM_Atmosphere_type),  ALLOCATABLE :: atm_AD(:), atm_AD_Ref(:)
    TYPE(CRTM_Surface_type),     ALLOCATABLE :: sfc_FWD(:)
    TYPE(CRTM_Surface_type),     ALLOCATABLE :: sfc_TL(:)
    TYPE(CRTM_Surface_type),     ALLOCATABLE :: sfc_K(:,:)
    TYPE(CRTM_Surface_type),     ALLOCATABLE :: sfc_AD(:), sfc_AD_Ref(:)
    TYPE(Timing_type) :: timing


    ! Get test dimensions
    n_channels = SUM(CRTM_ChannelInfo_n_Channels(chinfo))
    n_profiles = SIZE(atm)
    n_sensors  = SIZE(chinfo)
    n_layers   = SIZE(atm(1)%Temperature(:)) 
    
    ! Perform all the allocations
    ALLOCATE( rts_Base (n_channels, n_profiles, N_RT_ALGORITHMS), &
              rts_K   (n_channels, n_profiles, N_RT_ALGORITHMS), &
              rts_AD   (n_channels, n_profiles, N_RT_ALGORITHMS), &
              rts_TL   (n_channels, n_profiles, N_RT_ALGORITHMS), &
              atm_FWD  (n_Profiles), &
              atm_K   (n_channels, n_Profiles), &
              atm_AD   (n_Profiles), &
              atm_TL   (n_Profiles), &
              atm_AD_Ref   (n_Profiles), &
              sfc_FWD  (n_Profiles), &
              sfc_K   (n_channels, n_Profiles), &
              sfc_AD   (n_Profiles), &
              sfc_TL   (n_Profiles), &
              sfc_AD_Ref   (n_Profiles), &
              STAT   = alloc_stat, &
              ERRMSG = alloc_msg   )
    IF ( alloc_stat /= 0 ) THEN
      err_msg = 'Error allocating data structure arrays - '//TRIM(alloc_msg)
      CALL Display_Message(ROUTINE_NAME, err_msg, FAILURE ); STOP
    END IF

    call CRTM_RTSolution_Create(rts_Base, n_layers) 
    call CRTM_RTSolution_Create(rts_AD, n_layers) 
    call CRTM_RTSolution_Create(rts_K, n_layers) 
    call CRTM_RTSolution_Create(rts_TL, n_layers) 
    
    ! Loop over types of radiative transfer algorithms
    rt_algorithm_loop: DO i = 1, N_RT_ALGORITHMS               
      opt(:)%RT_Algorithm_Id = RT_ALGORITHM_ID(i)
    CALL CRTM_RTSolution_Zero( rts_Base(:,:,i) )
    CALL CRTM_RTSolution_Zero( rts_AD(:,:,i) )   
    CALL CRTM_RTSolution_Zero( rts_K(:,:,i) )  
    CALL CRTM_RTSolution_Zero( rts_TL(:,:,i) )       

       ! Setup for test failure reporting
      atm_file_status = 'REPLACE'
      atm_failure_file = TRIM(RT_ALGORITHM_NAME(i))//'.ADK.atmosphere.test_failure_report'
      sfc_file_status = 'REPLACE'
      sfc_failure_file = TRIM(RT_ALGORITHM_NAME(i))//'.ADK.surface.test_failure_report'
      

      ! Output info
      ! ...Algorithm identifier
      WRITE(*,'(30("*"),1x,"AD Comparisons for RT Algorithm ",a,1x,30("*"))') TRIM(RT_ALGORITHM_NAME(i))
      ! ...Sensors to process
      WRITE(*, '(4x,"- Sensors: ",99(a,:))') chinfo%sensor_id
 
      ! Output info
      ! ...Algorithm identifier
      WRITE(*,'(30("*"),1x,"K/AD Comparisons for RT Algorithm ",a,1x,30("*"))') TRIM(RT_ALGORITHM_NAME(i))
      ! ...Sensors to process
      WRITE(*, '(4x,"- Sensors: ",99(a,:))') chinfo%sensor_id


      ! Initialise adjoint structures
      atm_AD = atm
      CALL CRTM_Atmosphere_Zero( atm_AD )
      sfc_AD = sfc
      CALL CRTM_Surface_Zero( sfc_AD )

      DO kk = 1, n_channels
      DO ii = 1, n_profiles
        atm_K(kk,ii) = atm_AD(ii)
        sfc_K(kk,ii) = sfc_AD(ii)
      END DO
      END DO

      ! Copy the inputs so they can be modified if necessary
      atm_FWD = atm
      sfc_FWD = sfc

      atm_TL = atm
      sfc_TL = sfc

      ! For emission algorithm calculations
      IF ( i == EMISSION_INDEX ) THEN
        atm_FWD%n_Clouds   = 0
        atm_FWD%n_Aerosols = 0
      END IF
      
 771  FORMAT(I5,3f12.4,8ES13.4)

      CALL CRTM_Atmosphere_ZERO(atm_AD)
      CALL CRTM_Surface_ZERO(sfc_AD)
      CALL CRTM_RTSolution_Zero( rts_AD(:,:,i) )


      IF( IR_MW_Sensor ) THEN
        rts_AD(:,:,i)%Brightness_Temperature    = ONE
        rts_AD(:,:,i)%Radiance = ZERO
      ELSE
        IF( n_Stokes == 1 ) THEN
          rts_AD(:,:,i)%Radiance    = ONE
          rts_AD(:,:,i)%Brightness_Temperature = ZERO
        ELSE
          rts_AD(:,:,i)%Stokes(iStoke)    = ONE
        END IF
      END IF

      print *,' *********  AD ****************** '

      ! Perform AD calculations
      WRITE(*, '(4x,"- Running adjoint model...")')
      CALL Timing_Begin(timing)
      err_stat = CRTM_Adjoint( atm_FWD        , &
                               sfc_FWD        , &
                               rts_AD(:,:,i)  , &
                               geo            , &
                               chinfo         , &
                               atm_AD         , &
                               sfc_AD         , &
                               rts_Base(:,:,i), &
                               Options = opt(n1:n2)    )
            print *,' AD err_stat = ',err_stat

      CALL Timing_End(timing)
      IF ( err_stat /= SUCCESS ) THEN
        err_msg = 'Error when calling CRTM_Adjoint'
        CALL Display_Message(ROUTINE_NAME, err_msg, err_stat ); STOP
      END IF
      CALL Timing_Display(timing)

      CALL CRTM_Atmosphere_ZERO(atm_K)
      CALL CRTM_Surface_ZERO(sfc_K)
      CALL CRTM_RTSolution_Zero( rts_K(:,:,i) )

      IF( IR_MW_Sensor ) THEN
        rts_K(:,:,i)%Brightness_Temperature    = ONE
        rts_K(:,:,i)%Radiance = ZERO
      ELSE
        IF( n_Stokes == 1 ) THEN
          rts_K(:,:,i)%Radiance    = ONE
          rts_K(:,:,i)%Brightness_Temperature = ZERO
        ELSE
          rts_K(:,:,i)%Stokes(iStoke)    = ONE
        END IF
      END IF

      print *,' *********  K-Matrix ****************** '
     err_stat = CRTM_K_Matrix( Atm_FWD  , &  ! FWD Input
                           Sfc_FWD , &  ! FWD Input
                           rts_K(:,:,i), &  ! K   Input
                           geo , &  ! Input
                           chinfo , &  ! Input
                           atm_K  , &  ! K   Output
                           sfc_K  , &  ! K   Output
                           rts_Base(:,:,i), &  ! FWD Output
                           Options = Opt(n1:n2) ) 

      DO ii = 1, n_profiles
        atm_AD_Ref(ii) = atm_K(1,ii)
        sfc_AD_Ref(ii) = sfc_K(1,ii)
      END DO
      
      IF( n_channels > 1 ) THEN
      DO ii = 1, n_profiles
        DO kk = 2, n_channels
          atm_AD_Ref(ii) = atm_AD_Ref(ii) + atm_K(kk,ii)
          sfc_AD_Ref(ii) = sfc_AD_Ref(ii) + sfc_K(kk,ii)
        END DO
      END DO    
      END IF


      ! Initialise test
      WRITE(utest_msg,'("K/AD comparison test | ",&
                       &"Algorithm: ",a)') &
            TRIM(RT_ALGORITHM_NAME(i))
      CALL UnitTest_Setup(utest,TRIM(utest_msg))


      ! Perform test on each profile separately
      profile_loop: DO m = 1, n_profiles      
        CALL UnitTest_Assert( utest, CRTM_Atmosphere_Compare(atm_AD_Ref(m), atm_AD(m), n_SigFig=AD_SIGFIG) )

        ! Output info for failed test
        IF ( UnitTest_Failed(utest) ) THEN
          ! Open failure report file
          fid = Get_Lun()
          OPEN( fid, FILE     = atm_failure_file, &
                     FORM     = 'FORMATTED', &
                     STATUS   = atm_file_status, &
                     POSITION = 'APPEND', &
                     IOSTAT   = io_stat, &
                     IOMSG    = io_msg )
          IF ( io_stat /= 0 ) THEN
            err_msg = 'Error opening '//TRIM(atm_failure_file)//' - '//TRIM(io_msg)
            CALL Display_Message( ROUTINE_NAME, err_msg, WARNING )
          ELSE
            ! Update file status
            atm_file_status = 'OLD'
            ! Report failure
            WRITE(fid,'(40("="))')
            WRITE(fid,'("*** ATM AD test failed for profile #",i0)') m
            WRITE(fid,'("No. of significant figures used in comparison: ", i0)') AD_SIGFIG
            CALL CRTM_Atmosphere_Inspect(atm_AD_Ref(m)-atm_AD(m), Unit=fid)
            WRITE(fid,'("*** ATM AD test failed for profile #",i0)') m
            WRITE(fid,'(40("="),/)')
            print *,' qliu failed m = ',m
            STOP
            
            CLOSE(fid)
          END IF
        END IF


        ! Perform the surface test
        CALL UnitTest_Assert( utest, CRTM_Surface_Compare(sfc_AD_Ref(m), sfc_AD(m), n_SigFig=AD_SIGFIG) )


        ! Output info for failed test
        IF ( UnitTest_Failed(utest) ) THEN
          ! Open failure report file
          fid = Get_Lun()
          OPEN( fid, FILE     = sfc_failure_file, &
                     FORM     = 'FORMATTED', &
                     STATUS   = sfc_file_status, &
                     POSITION = 'APPEND', &
                     IOSTAT   = io_stat, &
                     IOMSG    = io_msg )
          IF ( io_stat /= 0 ) THEN
            err_msg = 'Error opening '//TRIM(sfc_failure_file)//' - '//TRIM(io_msg)
            CALL Display_Message( ROUTINE_NAME, err_msg, WARNING )
          ELSE
            ! Update file status
            sfc_file_status = 'OLD'
            ! Report failure
            WRITE(fid,'(40("="))')
            WRITE(fid,'("*** SFC AD test failed for profile #",i0)') m
            WRITE(fid,'("No. of significant figures used in comparison: ", i0)') AD_SIGFIG
            CALL CRTM_Surface_Inspect(sfc_AD_Ref(m)-sfc_AD(m), Unit=fid)
            WRITE(fid,'("*** SFC AD test failed for profile #",i0)') m
            WRITE(fid,'(40("="),/)')
            CLOSE(fid)
          END IF
        END IF

      END DO profile_loop

      CALL UnitTest_Report(utest)

    END DO rt_algorithm_loop


    ! Cleanup

    ! Destroy the data structures
    CALL CRTM_Atmosphere_Destroy(atm_FWD)
    CALL CRTM_Atmosphere_Destroy(atm_AD)
    CALL CRTM_Atmosphere_Destroy(atm_K)
    CALL CRTM_Atmosphere_Destroy(atm_AD_Ref)
    CALL CRTM_Surface_Destroy(sfc_FWD)
    CALL CRTM_Surface_Destroy(sfc_AD)
    CALL CRTM_Surface_Destroy(sfc_K)
    CALL CRTM_Surface_Destroy(sfc_AD_Ref)
    CALL CRTM_RTSolution_Destroy(rts_Base)
    CALL CRTM_RTSolution_Destroy(rts_AD)
   CALL CRTM_RTSolution_Destroy(rts_K)
    CALL CRTM_Options_Destroy(opt)


    ! Deallocate the structure arrays
    DEALLOCATE( rts_Base, &
                rts_AD  , &
                rts_K  , &
                atm_FWD , &
                atm_K  , &
                atm_AD  , &
                atm_AD_Ref  , &
                sfc_FWD , &
                sfc_K  , &
                sfc_AD  , &
                sfc_AD_Ref  , &
                STAT = alloc_stat )

  END SUBROUTINE Test_CRTM_ADK


  ! ========================================
  ! Subroutine to Assign TL Atmosphere input
  ! ========================================
  ELEMENTAL SUBROUTINE Assign_TL_Atmosphere( &
    dx           , & ! Input
    Atmosphere   , & ! Input
    Atmosphere_TL )  ! Output
    ! Arguments
    REAL(fp)                  , INTENT(IN)  :: dx
    TYPE(CRTM_Atmosphere_Type), INTENT(IN)  :: Atmosphere
    TYPE(CRTM_Atmosphere_Type), INTENT(OUT) :: Atmosphere_TL
    ! Local variables
    INTEGER :: j, n

    ! Initialize TL structure fields
    Atmosphere_TL = Atmosphere
    CALL CRTM_Atmosphere_Zero(Atmosphere_TL)

    ! Assign TL pressure, temperature values
    Atmosphere_TL%Temperature = Atmosphere%Temperature * DX
    Atmosphere_TL%Pressure    = Atmosphere%Pressure * DX
    ! Assign TL pressure, temperature values
    Atmosphere_TL%Cloud_Fraction = Atmosphere%Cloud_Fraction * DX
    ! assign TL absorber values
    DO j = 1, Atmosphere%n_Absorbers
      Atmosphere_TL%Absorber(:,j) = Atmosphere%Absorber(:,j) * DX
    END DO
    ! assign TL cloud values
    DO n = 1, Atmosphere%n_Clouds
      Atmosphere_TL%Cloud(n)%Water_Content    = Atmosphere%Cloud(n)%Water_Content * DX
      Atmosphere_TL%Cloud(n)%Effective_Radius = Atmosphere%Cloud(n)%Effective_Radius * DX
    END DO
    ! assign TL aerosol values
    DO n = 1, Atmosphere%n_Aerosols
      Atmosphere_TL%Aerosol(n)%Concentration    = Atmosphere%Aerosol(n)%Concentration * DX
      Atmosphere_TL%Aerosol(n)%Effective_Radius = Atmosphere%Aerosol(n)%Effective_Radius * DX
    END DO
  END SUBROUTINE Assign_TL_Atmosphere


  ! =====================================
  ! Subroutine to assign TL Surface input
  ! =====================================
  ELEMENTAL SUBROUTINE Assign_TL_Surface( &
    dx         , & ! Input
    Surface    , & ! Input
    Surface_TL )   ! Output
    ! Arguments
    REAL(fp)                , INTENT(IN)  :: dx
    TYPE(CRTM_Surface_Type) , INTENT(IN)  :: Surface
    TYPE(CRTM_Surface_Type) , INTENT(OUT) :: Surface_TL

    ! Initialize
    Surface_TL = Surface
    CALL CRTM_Surface_Zero(Surface_TL)

    ! Perturb the components
    IF ( Surface%Land_Coverage > ZERO ) THEN
      Surface_TL%Land_Coverage         = Surface%Land_Coverage
      Surface_TL%Land_Temperature      = dx*Surface%Land_Temperature
      Surface_TL%Soil_Moisture_Content = dx*Surface%Soil_Moisture_Content
      Surface_TL%Canopy_Water_Content  = dx*Surface%Canopy_Water_Content
      Surface_TL%Soil_Temperature      = dx*Surface%Soil_Temperature
    END IF
    IF ( Surface%Water_Coverage > ZERO ) THEN
      Surface_TL%Water_Coverage    = Surface%Water_Coverage
      Surface_TL%Water_Temperature = dx*Surface%Water_Temperature
      Surface_TL%Wind_Speed        = dx*Surface%Wind_Speed
      Surface_TL%Wind_Direction    = dx*Surface%Wind_Direction
      Surface_TL%Salinity          = dx*Surface%Salinity
    END IF
    IF ( Surface%Snow_Coverage > ZERO ) THEN
      Surface_TL%Snow_Coverage    = Surface%Snow_Coverage
      Surface_TL%Snow_Temperature = dx*Surface%Snow_Temperature
      Surface_TL%Snow_Depth       = dx*Surface%Snow_Depth
      Surface_TL%Snow_Density     = dx*Surface%Snow_Density
      Surface_TL%Snow_Grain_Size  = dx*Surface%Snow_Grain_Size
    END IF
    IF ( Surface%Ice_Coverage > ZERO ) THEN
      Surface_TL%Ice_Coverage    = Surface%Ice_Coverage
      Surface_TL%Ice_Temperature = dx*Surface%Ice_Temperature
      Surface_TL%Ice_Thickness   = dx*Surface%Ice_Thickness
      Surface_TL%Ice_Density     = dx*Surface%Ice_Density
      Surface_TL%Ice_Roughness   = dx*Surface%Ice_Roughness
    END IF
  END SUBROUTINE Assign_TL_Surface
!
!
!
!  ! Function to assign perturbed FWD Surface inputs
  ELEMENTAL SUBROUTINE Assign_Perturbed_Surface( &
    Surface    , & ! Input
    Surface_TL , & ! Input
    alpha      , & ! Input
    Surface_NLp, & ! Output
    Surface_NLm)   ! Output
    ! Arguments
    TYPE(CRTM_Surface_Type) , INTENT(IN)  :: Surface
    TYPE(CRTM_Surface_Type) , INTENT(IN)  :: Surface_TL
    REAL(fp)                , INTENT(IN)  :: alpha
    TYPE(CRTM_Surface_Type) , INTENT(OUT) :: Surface_NLp
    TYPE(CRTM_Surface_Type) , INTENT(OUT) :: Surface_NLm
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Assign_Perturbed_Surface'

    ! Initialize
    Surface_NLp = Surface
    Surface_NLm = Surface

    ! Perturb the components
    IF ( Surface%Land_Coverage > ZERO ) THEN

      Surface_NLm%Land_Temperature = Surface%Land_Temperature - &
        Surface_TL%Land_Temperature*alpha

      Surface_NLp%Land_Temperature = Surface%Land_Temperature + &
        Surface_TL%Land_Temperature*alpha

      Surface_NLm%Soil_Moisture_Content = Surface%Soil_Moisture_Content - &
        Surface_TL%Soil_Moisture_Content*alpha

      Surface_NLp%Soil_Moisture_Content = Surface%Soil_Moisture_Content + &
        Surface_TL%Soil_Moisture_Content*alpha

      Surface_NLm%Canopy_Water_Content = Surface%Canopy_Water_Content - &
        Surface_TL%Canopy_Water_Content*alpha

      Surface_NLp%Canopy_Water_Content = Surface%Canopy_Water_Content + &
        Surface_TL%Canopy_Water_Content*alpha

      Surface_NLm%Soil_Temperature = Surface%Soil_Temperature - &
        Surface_TL%Soil_Temperature*alpha

      Surface_NLp%Soil_Temperature = Surface%Soil_Temperature + &
        Surface_TL%Soil_Temperature*alpha
    END IF
    IF ( Surface%Water_Coverage > ZERO ) THEN
      Surface_NLm%Water_Temperature = Surface%Water_Temperature - &
        Surface_TL%Water_Temperature*alpha

      Surface_NLp%Water_Temperature = Surface%Water_Temperature + &
        Surface_TL%Water_Temperature*alpha

      Surface_NLm%Wind_Speed = Surface%Wind_Speed - &
        Surface_TL%Wind_Speed*alpha

      Surface_NLp%Wind_Speed = Surface%Wind_Speed + &
        Surface_TL%Wind_Speed*alpha

      Surface_NLm%Wind_Direction = Surface%Wind_Direction - &
        Surface_TL%Wind_Direction*alpha

      Surface_NLp%Wind_Direction = Surface%Wind_Direction + &
        Surface_TL%Wind_Direction*alpha

      Surface_NLm%Salinity = Surface%Salinity - &
        Surface_TL%Salinity*alpha

      Surface_NLp%Salinity = Surface%Salinity + &
        Surface_TL%Salinity*alpha
    END IF
    IF ( Surface%Snow_Coverage > ZERO ) THEN
      Surface_NLm%Snow_Temperature = Surface%Snow_Temperature - &
        Surface_TL%Snow_Temperature*alpha

      Surface_NLp%Snow_Temperature = Surface%Snow_Temperature + &
        Surface_TL%Snow_Temperature*alpha

      Surface_NLm%Snow_Depth = Surface%Snow_Depth - &
        Surface_TL%Snow_Depth*alpha

      Surface_NLp%Snow_Depth = Surface%Snow_Depth + &
        Surface_TL%Snow_Depth*alpha

      Surface_NLm%Snow_Density = Surface%Snow_Density - &
        Surface_TL%Snow_Density*alpha

      Surface_NLp%Snow_Density = Surface%Snow_Density + &
        Surface_TL%Snow_Density*alpha

      Surface_NLm%Snow_Grain_Size = Surface%Snow_Grain_Size - &
        Surface_TL%Snow_Grain_Size*alpha

      Surface_NLp%Snow_Grain_Size = Surface%Snow_Grain_Size + &
        Surface_TL%Snow_Grain_Size*alpha
    END IF
    IF ( Surface%Ice_Coverage > ZERO ) THEN
      Surface_NLm%Ice_Temperature = Surface%Ice_Temperature - &
        Surface_TL%Ice_Temperature*alpha

      Surface_NLp%Ice_Temperature = Surface%Ice_Temperature + &
        Surface_TL%Ice_Temperature*alpha

      Surface_NLm%Ice_Thickness = Surface%Ice_Thickness - &
        Surface_TL%Ice_Thickness*alpha

      Surface_NLp%Ice_Thickness = Surface%Ice_Thickness + &
        Surface_TL%Ice_Thickness*alpha

      Surface_NLm%Ice_Density = Surface%Ice_Density - &
        Surface_TL%Ice_Density*alpha

      Surface_NLp%Ice_Density = Surface%Ice_Density + &
        Surface_TL%Ice_Density*alpha

      Surface_NLm%Ice_Roughness = Surface%Ice_Roughness - &
        Surface_TL%Ice_Roughness*alpha

      Surface_NLp%Ice_Roughness = Surface%Ice_Roughness + &
        Surface_TL%Ice_Roughness*alpha
    END IF
  END SUBROUTINE Assign_Perturbed_Surface

END PROGRAM Test_CRTM_V30
