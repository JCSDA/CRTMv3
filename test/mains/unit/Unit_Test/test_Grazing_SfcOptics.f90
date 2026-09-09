!
! test_Grazing_SfcOptics
!
! Regression test for the catastrophic-reflectivity guard in the microwave
! water-surface optics: CRTM_FastemX and CRTM_PARMIO.
!
! Both backends derive V/H reflectivity by applying a polynomial reflection
! correction to the bare (1 - emissivity). That polynomial is fit for typical
! view angles; at the near-grazing Gaussian quadrature angles the scattering RT
! uses, and at high frequency, it extrapolates to a wildly non-physical
! reflectivity (~1e35) that blows the adding-doubling radiance up to ~ -1e15 K.
! Both backends now clamp grossly out-of-range reflectivity to the bare
! (1 - emissivity), which is physical by construction.
!
! This test drives each backend DIRECTLY at a sweep of grazing zenith angles at
! 325 GHz (where the uncaught correction blows up) and asserts:
!   * the returned reflectivity stays within the guard band [R_LO,R_HI]
!     -- i.e. no blow-up; this fails if the clamp is removed (raw ~1e35),
!   * the guard actually fired (at >=1 grazing angle the output equals the bare
!     (1 - emissivity) fall-back, which only happens when the clamp engages), and
!   * emissivity is physical.
! It also checks a non-grazing view angle, where the reflectivity must be
! physical [0,1] (normal behaviour preserved).
!
! Note on n_Angles: FASTEM caps its correction angle at 60 deg when n_Angles>1,
! so in the real (multi-stream) RT it is protected at grazing; its catastrophic
! regime is only reachable at n_Angles==1, which this test uses to exercise the
! FASTEM guard. PARMIO has no such cap and blows up at grazing in the real RT
! (the bug that produced -1e15 K AWS 325 GHz TBs).
!
! STOP 0 on success, STOP 1 on failure.
!
! CREATION HISTORY:
!       Written by:     Benjamin Johnson, 06-Jun-2026
!

PROGRAM test_Grazing_SfcOptics

  USE CRTM_Module
  USE CRTM_MWwaterCoeff, ONLY: MWwaterC
  USE CRTM_PARMIOCoeff,  ONLY: PARMIOC, CRTM_PARMIOCoeff_IsLoaded
  USE CRTM_FastemX,      ONLY: Fastem_iVar => iVar_type, Compute_FastemX
  USE CRTM_PARMIO,       ONLY: Parmio_iVar => iVar_type, Compute_PARMIO
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_Grazing_SfcOptics'
  CHARACTER(*), PARAMETER :: PATH   = './testinput/'
  CHARACTER(*), PARAMETER :: SENSOR = 'mwr_aws'

  ! Surface state + the high frequency that drives the correction out of range
  REAL(fp), PARAMETER :: FREQ  = 325.0_fp        ! GHz
  REAL(fp), PARAMETER :: SST   = 290.0_fp        ! K
  REAL(fp), PARAMETER :: SSS   = 33.0_fp         ! psu
  REAL(fp), PARAMETER :: WIND  = 6.0_fp          ! m/s
  REAL(fp), PARAMETER :: TRANS = 0.5_fp          ! atmospheric transmittance (enables the correction)
  ! Guard band -- must match R_PHYS_LO/R_PHYS_HI in CRTM_FastemX / CRTM_PARMIO
  REAL(fp), PARAMETER :: R_LO  = -0.5_fp, R_HI = 1.5_fp
  REAL(fp), PARAMETER :: FB_TOL = 1.0e-8_fp      ! fall-back match tolerance
  ! Grazing-angle sweep (up to ~86 deg, the largest quadrature angle the
  ! scattering RT actually uses) + one well-behaved view angle. 82 deg is below
  ! the blow-up threshold (guard idle); 84-86 deg trigger it.
  INTEGER,  PARAMETER :: NANG = 3
  REAL(fp), PARAMETER :: ZA(NANG) = (/ 82.0_fp, 84.0_fp, 86.0_fp /)
  REAL(fp), PARAMETER :: ZA_VIEW = 50.0_fp

  TYPE(CRTM_ChannelInfo_type) :: chinfo(1)
  TYPE(Fastem_iVar) :: fvar
  TYPE(Parmio_iVar) :: pvar
  REAL(fp) :: emis(4), refl(4)
  INTEGER  :: err, i
  LOGICAL  :: ok, fastem_fired, parmio_fired

  ok = .TRUE.; fastem_fired = .FALSE.; parmio_fired = .FALSE.

  ! Load the MW-water surface coefficients (FASTEM always; PARMIO LUT if staged)
  err = CRTM_Init( (/ SENSOR /), chinfo, File_Path=PATH, Quiet=.TRUE. )
  IF ( err /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Init failed', FAILURE ); STOP 1
  END IF
  IF ( .NOT. CRTM_PARMIOCoeff_IsLoaded() ) THEN
    CALL Display_Message( PROGRAM_NAME, 'PARMIO LUT not loaded; cannot test PARMIO path', FAILURE )
    STOP 1
  END IF

  ! --------------------------------------------------------------------------
  ! Grazing-angle sweep -- both backends must stay bounded; the guard must fire
  ! --------------------------------------------------------------------------
  WRITE(*,'(/a)') ' Grazing-angle reflectivity guard (325 GHz, ocean):'
  WRITE(*,'(a)')  '   backend  za      Rv          Rh        fell_back'
  DO i = 1, NANG
    ! FASTEM (n_Angles=1 bypasses the 60-deg correction-angle cap)
    emis = ZERO; refl = ZERO
    CALL Compute_FastemX( MWwaterC, FREQ, 1, ZA(i), SST, SSS, WIND, fvar, &
                          emis, refl, Transmittance=TRANS )
    CALL Check( 'FASTEM', ZA(i), emis, refl, fastem_fired )

    ! PARMIO at the same angle (no grazing cap; blows up in the real RT)
    emis = ZERO; refl = ZERO
    CALL Compute_PARMIO( PARMIOC, FREQ, 1, ZA(i), SST, SSS, WIND, pvar, &
                         emis, refl, Transmittance=TRANS )
    CALL Check( 'PARMIO', ZA(i), emis, refl, parmio_fired )
  END DO

  ! --------------------------------------------------------------------------
  ! Well-behaved view angle -- reflectivity must be physical [0,1]
  ! --------------------------------------------------------------------------
  emis = ZERO; refl = ZERO
  CALL Compute_FastemX( MWwaterC, FREQ, 1, ZA_VIEW, SST, SSS, WIND, fvar, &
                        emis, refl, Transmittance=TRANS )
  IF ( ANY(refl(1:2) < ZERO) .OR. ANY(refl(1:2) > ONE) ) THEN
    WRITE(*,'(a,2es12.3)') ' FAIL: FASTEM view-angle reflectivity not physical: ', refl(1:2); ok=.FALSE.
  END IF
  emis = ZERO; refl = ZERO
  CALL Compute_PARMIO( PARMIOC, FREQ, 1, ZA_VIEW, SST, SSS, WIND, pvar, &
                       emis, refl, Transmittance=TRANS )
  IF ( ANY(refl(1:2) < ZERO) .OR. ANY(refl(1:2) > ONE) ) THEN
    WRITE(*,'(a,2es12.3)') ' FAIL: PARMIO view-angle reflectivity not physical: ', refl(1:2); ok=.FALSE.
  END IF

  ! --------------------------------------------------------------------------
  ! The guard must have actually fired for both backends (else the test is not
  ! exercising the catastrophic regime).
  ! --------------------------------------------------------------------------
  IF ( .NOT. fastem_fired ) THEN
    WRITE(*,'(a)') ' FAIL: FASTEM reflectivity guard never engaged at grazing'; ok=.FALSE.
  END IF
  IF ( .NOT. parmio_fired ) THEN
    WRITE(*,'(a)') ' FAIL: PARMIO reflectivity guard never engaged at grazing'; ok=.FALSE.
  END IF

  err = CRTM_Destroy( chinfo )

  IF ( ok ) THEN
    WRITE(*,'(/a)') ' PASS: grazing-angle reflectivity guard holds for FASTEM and PARMIO.'
    STOP 0
  ELSE
    WRITE(*,'(/a)') ' FAIL: grazing-angle reflectivity guard test failed.'
    STOP 1
  END IF

CONTAINS

  ! Assert reflectivity bounded (guard band) and emissivity physical; set `fired`
  ! when the guard engaged -- detected by the output reflectivity matching the
  ! bare (1 - emissivity) fall-back (only produced when the clamp triggers).
  SUBROUTINE Check( tag, za_deg, e, r, fired )
    CHARACTER(*), INTENT(IN)    :: tag
    REAL(fp),     INTENT(IN)    :: za_deg, e(:), r(:)
    LOGICAL,      INTENT(INOUT) :: fired
    LOGICAL :: fell_back
    fell_back = ( ABS(r(1) - (ONE-e(1))) < FB_TOL ) .AND. &
                ( ABS(r(2) - (ONE-e(2))) < FB_TOL )
    IF ( ANY(r(1:2) < R_LO) .OR. ANY(r(1:2) > R_HI) ) THEN
      WRITE(*,'(a,a,a,f5.1,a,2es12.3)') ' FAIL: ',tag,' reflectivity out of band at za=',za_deg,' : ',r(1:2)
      ok = .FALSE.
    END IF
    IF ( ANY(e(1:2) < ZERO) .OR. ANY(e(1:2) > ONE) ) THEN
      WRITE(*,'(a,a,a,f5.1,a,2es12.3)') ' FAIL: ',tag,' emissivity not physical at za=',za_deg,' : ',e(1:2)
      ok = .FALSE.
    END IF
    IF ( fell_back ) fired = .TRUE.
    WRITE(*,'(3x,a,2x,f5.1,2es12.3,4x,l1)') tag, za_deg, r(1:2), fell_back
  END SUBROUTINE Check

END PROGRAM test_Grazing_SfcOptics
