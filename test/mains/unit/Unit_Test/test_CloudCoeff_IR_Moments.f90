!
! test_CloudCoeff_IR_Moments
!
! Consistency test for the IR cloud phase coefficients of CloudCoeff_MC6.nc
! (JCSDA/CRTMv3#394). Each IR stream set is a delta-fit: slot 0 holds 0.5, slots
! 1..n-1 hold (l+1/2) chi'_l with chi'_l = (chi_l - f)/(1 - f), and slot n holds the
! truncation factor f. The fit keeps the first moments of the phase function, so the
! implied unscaled first moment chi_1 = c_1/1.5 (1 - f) + f must equal the table's own
! asymmetry parameter g_IR. In the default CloudCoeff.nc this fails for IR class 4
! (ICE_CLOUD, HAIL_CLOUD), whose sets have chi_1 = chi_2 = 0 at every point (#394);
! CloudCoeff_MC6.nc replaces that class with the MC6 ice optics.
!
! Checks, for every IR density class and every stored stream set of CloudCoeff_MC6.nc
! (the 2-stream entry shares the 4-stream set's storage at offset 0, so its slot 2 is
! c_2, not a truncation factor, and it is skipped):
!   1. c_0 = 0.5 at the set's Legendre_Offset
!   2. median |chi_1 - g_IR| < 1e-3
!   3. fewer than 5% of the frequency/radius points have |chi_1 - g_IR| > 0.02
!      (the delta-fit quadrature misses the narrowest UV forward peaks at a few points)
!
! STOP 0 on success, STOP 1 on failure.
!
PROGRAM test_CloudCoeff_IR_Moments

  USE CRTM_CloudCoeff, ONLY: CRTM_CloudCoeff_Load, CloudC
  USE Message_Handler, ONLY: SUCCESS, FAILURE, Display_Message
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_CloudCoeff_IR_Moments'
  CHARACTER(*), PARAMETER :: CLOUDCOEFF_FILE = 'CloudCoeff_MC6.nc'
  CHARACTER(*), PARAMETER :: FILE_PATH = './testinput/'
  INTEGER, PARAMETER :: DP = KIND(1.0d0)
  REAL(DP), PARAMETER :: TOLERANCE = 1.0d-5
  REAL(DP), PARAMETER :: MEDIAN_LIMIT = 1.0d-3
  REAL(DP), PARAMETER :: OUTLIER_LIMIT = 2.0d-2
  REAL(DP), PARAMETER :: OUTLIER_FRACTION = 0.05d0
  CHARACTER(256) :: msg
  REAL(DP), ALLOCATABLE :: f(:,:), chi1(:,:), e(:)
  REAL(DP) :: med, frac
  INTEGER :: err_stat, k, m, n, offset, n_fail

  err_stat = CRTM_CloudCoeff_Load( 'CRTM', CLOUDCOEFF_FILE, FILE_PATH, netCDF = .TRUE., Quiet = .TRUE. )
  IF ( err_stat /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error loading '//CLOUDCOEFF_FILE, FAILURE )
    STOP 1
  END IF

  n_fail = 0
  DO k = 1, SIZE(CloudC%g_IR, 3)
    DO m = 1, CloudC%n_Stream_Sets
      offset = CloudC%Legendre_Offset(m)
      n      = CloudC%n_Streams(m)
      ! Skip a set stored inside a larger one (same offset, more streams)
      IF ( ANY(CloudC%Legendre_Offset(1:CloudC%n_Stream_Sets) == offset .AND. &
               CloudC%n_Streams(1:CloudC%n_Stream_Sets) > n) ) CYCLE
      ! 1. Normalisation slot
      IF ( ANY(ABS(CloudC%pcoeff_IR(:,:,k,offset,1) - 0.5_DP) > TOLERANCE) ) THEN
        WRITE( msg, '("IR class ",i0,", ",i0,"-stream set: c_0 is not 0.5")' ) k, n
        CALL Display_Message( PROGRAM_NAME, msg, FAILURE )
        n_fail = n_fail + 1
      END IF
      ! 2.-3. Implied first moment against g_IR
      f    = CloudC%pcoeff_IR(:,:,k,offset+n,1)
      chi1 = CloudC%pcoeff_IR(:,:,k,offset+1,1) / 1.5_DP * (1.0_DP - f) + f
      e    = PACK(ABS(chi1 - CloudC%g_IR(:,:,k)), .TRUE.)
      med  = Median(e)
      frac = REAL(COUNT(e > OUTLIER_LIMIT), DP) / REAL(SIZE(e), DP)
      IF ( med >= MEDIAN_LIMIT .OR. frac >= OUTLIER_FRACTION ) THEN
        WRITE( msg, '("IR class ",i0,", ",i0,"-stream set: median |chi_1 - g| = ",es9.2, &
                     &", fraction above ",f4.2," = ",f5.3)' ) k, n, med, OUTLIER_LIMIT, frac
        CALL Display_Message( PROGRAM_NAME, msg, FAILURE )
        n_fail = n_fail + 1
      END IF
    END DO
  END DO

  IF ( n_fail > 0 ) STOP 1
  WRITE(*,'(a)') PROGRAM_NAME//': every IR stream set of '//CLOUDCOEFF_FILE//' is consistent with g_IR'
  STOP 0

CONTAINS

  ! Median by insertion sort (a few hundred values per set)
  FUNCTION Median(x) RESULT(xm)
    REAL(DP), INTENT(IN) :: x(:)
    REAL(DP) :: xm, a(SIZE(x)), t
    INTEGER :: i, j
    a = x
    DO i = 2, SIZE(a)
      t = a(i)
      j = i - 1
      DO WHILE ( j >= 1 )
        IF ( a(j) <= t ) EXIT
        a(j+1) = a(j)
        j = j - 1
      END DO
      a(j+1) = t
    END DO
    xm = a((SIZE(a)+1)/2)
  END FUNCTION Median

END PROGRAM test_CloudCoeff_IR_Moments
