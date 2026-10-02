!-------------------------------------------------------
!
! Description:
!       Verification test for the AOD Skip_Profile handling
!       backported to CRTM v3.1.6 (#372).
!
!       Before v3.1.6 only CRTM_AOD honoured Options%Skip_Profile.
!       CRTM_AOD_TL, CRTM_AOD_AD and CRTM_AOD_K computed a flagged
!       profile anyway, and returned FAILURE when its input was
!       invalid.
!
!       This test runs all four AOD functions on three aerosol
!       profiles (v.abi_g18), first with every profile valid and
!       none flagged (the reference), then with profile 2 given a
!       negative pressure and flagged with Skip_Profile, and asserts:
!
!         1. The corrupted profile really is invalid.
!         2. CRTM_AOD, _TL, _AD and _K each return SUCCESS.
!         3. Profiles 1 and 3 are bit-identical to the reference
!            for every function.
!         4. The outputs for the skipped profile are left exactly
!            as the caller set them.
!
!       Assertion 2 fails for _TL, _AD and _K on an unfixed 3.1.x
!       library, guarding against regression of the backport.
!
!       Date: 2026-10-01
!
!-------------------------------------------------------

PROGRAM test_aod_skip_profile

  ! Module usage
  USE CRTM_Module
  USE UnitTest_Define, ONLY: UnitTest_type,     &
                             UnitTest_Init,     &
                             UnitTest_Setup,    &
                             UnitTest_Assert,   &
                             UnitTest_Report,   &
                             UnitTest_n_Failed

  ! Disable all implicit typing
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: Program_Name = 'test_aod_skip_profile'
  CHARACTER(*), PARAMETER :: SENSOR_ID    = 'v.abi_g18'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'

  ! Profile setup (matches the AOD regression tests, plus a third profile)
  INTEGER,  PARAMETER :: N_PROFILES  = 3
  INTEGER,  PARAMETER :: N_LAYERS    = 92
  INTEGER,  PARAMETER :: N_ABSORBERS = 2
  INTEGER,  PARAMETER :: N_CLOUDS    = 0
  INTEGER,  PARAMETER :: N_AEROSOLS  = 3
  INTEGER,  PARAMETER :: SKIPPED     = 2
  ! Marker the caller leaves in the skipped profile's outputs
  REAL(fp), PARAMETER :: SENTINEL    = -999.0_fp

  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(1)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES), Atm_Ref(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm_TL(N_PROFILES), Atm_AD(N_PROFILES)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atm_K(:,:)
  TYPE(CRTM_Options_type)                 :: Options(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:), rts_TL(:,:), rts_AD(:,:), rts_K(:,:)

  ! Reference outputs: AOD and TL layer optical depths, AD and K aerosol sensitivities
  REAL(fp), ALLOCATABLE :: ref_fwd(:,:,:), ref_tl(:,:,:), ref_ad(:,:,:), ref_k(:,:,:,:)
  REAL(fp), ALLOCATABLE :: out_fwd(:,:,:), out_tl(:,:,:), out_ad(:,:,:), out_k(:,:,:,:)

  TYPE(UnitTest_type) :: utest
  INTEGER :: err, stat(4), n_Channels, l, m
  LOGICAL :: valid, untouched

  CALL UnitTest_Init(utest, .TRUE.)
  CALL UnitTest_Setup(utest, 'AOD_Skip_Profile_Test', Program_Name, .TRUE.)

  err = CRTM_Init( (/SENSOR_ID/), ChannelInfo, File_Path=COEFFICIENTS_PATH, Quiet=.TRUE. )
  IF ( err /= SUCCESS ) THEN
    CALL Display_Message( Program_Name, 'Error initializing CRTM', FAILURE )
    STOP 1
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))

  ! Profiles 1 and 2 from the AOD regression loader; profile 3 is
  ! profile 1 with half the aerosol loading, so all three differ.
  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
    CALL Display_Message( Program_Name, 'Error creating Atmosphere', FAILURE )
    STOP 1
  END IF
  CALL Load_Atm_Data()
  Atm(3) = Atm(1)
  DO l = 1, N_AEROSOLS
    Atm(3)%Aerosol(l)%Concentration = 0.5_fp * Atm(1)%Aerosol(l)%Concentration
  END DO
  Atm_Ref = Atm

  ! Perturbation, adjoint and K-matrix structures
  CALL CRTM_Atmosphere_Create( Atm_TL, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atm_AD, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  ALLOCATE( Atm_K(n_Channels, N_PROFILES), &
            rts(n_Channels, N_PROFILES), rts_TL(n_Channels, N_PROFILES), &
            rts_AD(n_Channels, N_PROFILES), rts_K(n_Channels, N_PROFILES) )
  CALL CRTM_Atmosphere_Create( Atm_K, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_RTSolution_Create( rts,    N_LAYERS )
  CALL CRTM_RTSolution_Create( rts_TL, N_LAYERS )
  CALL CRTM_RTSolution_Create( rts_AD, N_LAYERS )
  CALL CRTM_RTSolution_Create( rts_K,  N_LAYERS )
  ALLOCATE( ref_fwd(N_LAYERS, n_Channels, N_PROFILES), ref_tl(N_LAYERS, n_Channels, N_PROFILES), &
            ref_ad(N_LAYERS, N_AEROSOLS, N_PROFILES),  ref_k(N_LAYERS, N_AEROSOLS, n_Channels, N_PROFILES), &
            out_fwd(N_LAYERS, n_Channels, N_PROFILES), out_tl(N_LAYERS, n_Channels, N_PROFILES), &
            out_ad(N_LAYERS, N_AEROSOLS, N_PROFILES),  out_k(N_LAYERS, N_AEROSOLS, n_Channels, N_PROFILES) )


  ! 1. Reference: all profiles valid, none flagged
  ! ----------------------------------------------
  CALL Run_All_AOD( stat, ref_fwd, ref_tl, ref_ad, ref_k )
  IF ( ANY(stat /= SUCCESS) ) THEN
    CALL Display_Message( Program_Name, 'Reference AOD calls failed', FAILURE )
    STOP 1
  END IF


  ! 2. Profile SKIPPED invalid and flagged with Skip_Profile
  ! --------------------------------------------------------
  Atm = Atm_Ref
  Atm(SKIPPED)%Pressure(10)       = -100.0_fp
  Atm(SKIPPED)%Level_Pressure(10) = -100.0_fp
  Options(SKIPPED)%Skip_Profile   = .TRUE.

  WRITE(*,'(/a)') 'ASSERT: the flagged profile has invalid input'
  valid = CRTM_Atmosphere_IsValid( Atm(SKIPPED) )
  CALL UnitTest_Assert(utest, .NOT. valid)

  CALL Run_All_AOD( stat, out_fwd, out_tl, out_ad, out_k )

  WRITE(*,'(/a)') 'ASSERT: CRTM_AOD returns SUCCESS'
  CALL UnitTest_Assert(utest, stat(1) == SUCCESS)
  WRITE(*,'(/a)') 'ASSERT: CRTM_AOD_TL returns SUCCESS'
  CALL UnitTest_Assert(utest, stat(2) == SUCCESS)
  WRITE(*,'(/a)') 'ASSERT: CRTM_AOD_AD returns SUCCESS'
  CALL UnitTest_Assert(utest, stat(3) == SUCCESS)
  WRITE(*,'(/a)') 'ASSERT: CRTM_AOD_K returns SUCCESS'
  CALL UnitTest_Assert(utest, stat(4) == SUCCESS)

  WRITE(*,'(/a)') 'ASSERT: CRTM_AOD unflagged profiles identical to reference'
  CALL UnitTest_Assert(utest, ALL(out_fwd(:,:,1) == ref_fwd(:,:,1)) .AND. ALL(out_fwd(:,:,3) == ref_fwd(:,:,3)))
  WRITE(*,'(/a)') 'ASSERT: CRTM_AOD_TL unflagged profiles identical to reference'
  CALL UnitTest_Assert(utest, ALL(out_tl(:,:,1) == ref_tl(:,:,1)) .AND. ALL(out_tl(:,:,3) == ref_tl(:,:,3)))
  WRITE(*,'(/a)') 'ASSERT: CRTM_AOD_AD unflagged profiles identical to reference'
  CALL UnitTest_Assert(utest, ALL(out_ad(:,:,1) == ref_ad(:,:,1)) .AND. ALL(out_ad(:,:,3) == ref_ad(:,:,3)))
  WRITE(*,'(/a)') 'ASSERT: CRTM_AOD_K unflagged profiles identical to reference'
  CALL UnitTest_Assert(utest, ALL(out_k(:,:,:,1) == ref_k(:,:,:,1)) .AND. ALL(out_k(:,:,:,3) == ref_k(:,:,:,3)))

  WRITE(*,'(/a)') 'ASSERT: skipped profile outputs left as the caller set them'
  untouched = ALL(out_fwd(:,:,SKIPPED) == SENTINEL) .AND. ALL(out_tl(:,:,SKIPPED) == SENTINEL) .AND. &
              ALL(out_ad(:,:,SKIPPED)  == SENTINEL) .AND. ALL(out_k(:,:,:,SKIPPED) == SENTINEL)
  CALL UnitTest_Assert(utest, untouched)


  ! 3. Clean up and report
  ! ----------------------
  err = CRTM_Destroy( ChannelInfo )
  CALL CRTM_Atmosphere_Destroy( Atm )
  CALL CRTM_Atmosphere_Destroy( Atm_Ref )
  CALL CRTM_Atmosphere_Destroy( Atm_TL )
  CALL CRTM_Atmosphere_Destroy( Atm_AD )
  CALL CRTM_Atmosphere_Destroy( Atm_K )

  CALL UnitTest_Report(utest)
  IF ( UnitTest_n_Failed(utest) > 0 ) STOP 1

CONTAINS

  ! Call CRTM_AOD, _TL, _AD and _K with the current Atm and Options, after
  ! setting every output for the flagged profile (if any) to SENTINEL.
  SUBROUTINE Run_All_AOD( Stat, Fwd, TL, AD, K )
    INTEGER,  INTENT(OUT) :: Stat(4)
    REAL(fp), INTENT(OUT) :: Fwd(:,:,:), TL(:,:,:), AD(:,:,:), K(:,:,:,:)
    INTEGER :: ll, mm, na

    ! TL input: 1% perturbation of aerosol concentration and radius
    CALL CRTM_Atmosphere_Zero( Atm_TL )
    DO mm = 1, N_PROFILES
      DO na = 1, N_AEROSOLS
        Atm_TL(mm)%Aerosol(na)%Type             = Atm(mm)%Aerosol(na)%Type
        Atm_TL(mm)%Aerosol(na)%Concentration    = 0.01_fp * Atm_Ref(mm)%Aerosol(na)%Concentration
        Atm_TL(mm)%Aerosol(na)%Effective_Radius = 0.01_fp * Atm_Ref(mm)%Aerosol(na)%Effective_Radius
      END DO
    END DO

    ! Adjoint and K-matrix inputs; the flagged profile's outputs get the sentinel
    CALL CRTM_Atmosphere_Zero( Atm_AD )
    CALL CRTM_Atmosphere_Zero( Atm_K )
    CALL CRTM_RTSolution_Zero( rts )
    CALL CRTM_RTSolution_Zero( rts_TL )
    DO mm = 1, N_PROFILES
      DO ll = 1, n_Channels
        rts_AD(ll,mm)%Layer_Optical_Depth = ONE
        rts_K(ll,mm)%Layer_Optical_Depth  = ONE
      END DO
      IF ( Options(mm)%Skip_Profile ) CALL Set_Sentinel( mm )
    END DO

    Stat(1) = CRTM_AOD( Atm, ChannelInfo, rts, Options=Options )
    DO ll = 1, n_Channels
      Fwd(:,ll,:) = RESHAPE([(rts(ll,mm)%Layer_Optical_Depth, mm=1,N_PROFILES)], [N_LAYERS,N_PROFILES])
    END DO

    Stat(2) = CRTM_AOD_TL( Atm, Atm_TL, ChannelInfo, rts, rts_TL, Options=Options )
    DO ll = 1, n_Channels
      TL(:,ll,:) = RESHAPE([(rts_TL(ll,mm)%Layer_Optical_Depth, mm=1,N_PROFILES)], [N_LAYERS,N_PROFILES])
    END DO

    Stat(3) = CRTM_AOD_AD( Atm, rts_AD, ChannelInfo, rts, Atm_AD, Options=Options )
    DO mm = 1, N_PROFILES
      DO na = 1, N_AEROSOLS
        AD(:,na,mm) = Atm_AD(mm)%Aerosol(na)%Concentration
      END DO
    END DO

    Stat(4) = CRTM_AOD_K( Atm, rts_K, ChannelInfo, rts, Atm_K, Options=Options )
    DO mm = 1, N_PROFILES
      DO ll = 1, n_Channels
        DO na = 1, N_AEROSOLS
          K(:,na,ll,mm) = Atm_K(ll,mm)%Aerosol(na)%Concentration
        END DO
      END DO
    END DO
  END SUBROUTINE Run_All_AOD


  ! Mark every output the AOD functions could write for profile mm
  SUBROUTINE Set_Sentinel( mm )
    INTEGER, INTENT(IN) :: mm
    INTEGER :: ll, na
    DO ll = 1, n_Channels
      rts(ll,mm)%Layer_Optical_Depth    = SENTINEL
      rts_TL(ll,mm)%Layer_Optical_Depth = SENTINEL
      DO na = 1, N_AEROSOLS
        Atm_K(ll,mm)%Aerosol(na)%Concentration = SENTINEL
      END DO
    END DO
    DO na = 1, N_AEROSOLS
      Atm_AD(mm)%Aerosol(na)%Concentration = SENTINEL
    END DO
  END SUBROUTINE Set_Sentinel


  INCLUDE 'Load_Atm_Data.inc'

END PROGRAM test_aod_skip_profile
