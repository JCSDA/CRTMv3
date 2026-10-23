!
! test_CloudCoeff_IR_Legendre_Layout
!
! Regression test for JCSDA/CRTMv3#384: the IR cloud phase coefficients must be
! laid out like the MW ones. CRTM_CloudScatter reads pcoeff(l + Legendre_Offset)
! for l >= 1 from every table, which assumes a 0-based Legendre dimension in which
! each stream set starts with its normalised c_0 = 0.5 at Legendre_Offset. With a
! 1-based pcoeff_IR (v3.1.0 to v3.1.6) every IR set was read one Legendre order low
! (c_1 = 0.5 for all particles).
!
! Checks, for the loaded CloudCoeff:
!   1. pcoeff_IR and pcoeff_S_MW have the same (0) Legendre lower bound
!   2. every IR stream set has c_0 = 0.5 at its Legendre_Offset
!   3. c_1 of the IR sets is not the normalisation 0.5 everywhere
!
! STOP 0 on success, STOP 1 on failure.
!
PROGRAM test_CloudCoeff_IR_Legendre_Layout

  USE CRTM_CloudCoeff, ONLY: CRTM_CloudCoeff_Load, CloudC
  USE Message_Handler, ONLY: SUCCESS, FAILURE, Display_Message
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_CloudCoeff_IR_Legendre_Layout'
  CHARACTER(*), PARAMETER :: CLOUDCOEFF_FILE = 'CloudCoeff.nc'
  CHARACTER(*), PARAMETER :: FILE_PATH = './testinput/'
  REAL(KIND(1.0d0)), PARAMETER :: TOLERANCE = 1.0d-5
  CHARACTER(256) :: msg
  INTEGER :: err_stat, m, offset, n_fail

  err_stat = CRTM_CloudCoeff_Load( 'CRTM', CLOUDCOEFF_FILE, FILE_PATH, netCDF = .TRUE., Quiet = .TRUE. )
  IF ( err_stat /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error loading '//CLOUDCOEFF_FILE, FAILURE )
    STOP 1
  END IF

  n_fail = 0

  ! 1. Same Legendre indexing for IR and MW (CRTM_CloudScatter indexes both alike)
  IF ( LBOUND(CloudC%pcoeff_IR, 4) /= LBOUND(CloudC%pcoeff_S_MW, 4) ) THEN
    WRITE( msg, '("pcoeff_IR Legendre lower bound ",i0," differs from pcoeff_S_MW ",i0)' ) &
      LBOUND(CloudC%pcoeff_IR, 4), LBOUND(CloudC%pcoeff_S_MW, 4)
    CALL Display_Message( PROGRAM_NAME, msg, FAILURE )
    n_fail = n_fail + 1
  END IF

  DO m = 1, CloudC%n_Stream_Sets
    offset = CloudC%Legendre_Offset(m)
    IF ( offset < LBOUND(CloudC%pcoeff_IR, 4) .OR. offset + 1 > UBOUND(CloudC%pcoeff_IR, 4) ) THEN
      WRITE( msg, '(i0,"-stream set offset ",i0," outside the pcoeff_IR Legendre bounds")' ) &
        CloudC%n_Streams(m), offset
      CALL Display_Message( PROGRAM_NAME, msg, FAILURE )
      n_fail = n_fail + 1
      CYCLE
    END IF
    ! 2. Each IR stream set starts with c_0 = 0.5 at its offset
    IF ( ANY(ABS(CloudC%pcoeff_IR(:,:,:,offset,1) - 0.5d0) > TOLERANCE) ) THEN
      WRITE( msg, '(i0,"-stream IR set: c_0 at Legendre index ",i0," is not 0.5")' ) CloudC%n_Streams(m), offset
      CALL Display_Message( PROGRAM_NAME, msg, FAILURE )
      n_fail = n_fail + 1
    END IF
    ! 3. ...and its c_1 is a real first moment, not the normalisation term
    IF ( ALL(ABS(CloudC%pcoeff_IR(:,:,:,offset+1,1) - 0.5d0) <= TOLERANCE) ) THEN
      WRITE( msg, '(i0,"-stream IR set: c_1 is 0.5 everywhere (Legendre index shifted)")' ) CloudC%n_Streams(m)
      CALL Display_Message( PROGRAM_NAME, msg, FAILURE )
      n_fail = n_fail + 1
    END IF
  END DO

  IF ( n_fail > 0 ) STOP 1
  WRITE(*,'(a)') PROGRAM_NAME//': IR Legendre layout matches MW; all stream sets start with c_0 = 0.5'
  STOP 0

END PROGRAM test_CloudCoeff_IR_Legendre_Layout
