!
! test_Long_Path_Init
!
! Regression coverage for issue #238: coefficient file paths must not be
! silently truncated by fixed-length buffers. CRTM is initialized through a
! deliberately long File_Path (a two-level symlink to the test-data tree,
! about 300 characters) that exceeds the old fixed caps (80/128/256) in the
! coefficient loaders. With the deferred-length path handling the load must
! succeed; with the old fixed-length buffers the path was clipped and the
! files were not found.
!
! Builds ./<150 a>/<150 b> as a symlink to ../testinput, then initializes
! atms_n21 from that path (default loads: SpcCoeff, TauCoeff, CloudCoeff,
! AerosolCoeff, and the surface emissivity coefficients, exercising the
! converted loaders).
!
! STOP 0 on success, STOP 1 on failure.
!
PROGRAM test_Long_Path_Init

  USE CRTM_Module
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_Long_Path_Init'
  CHARACTER(*), PARAMETER :: SENSOR       = 'atms_n21'
  INTEGER,      PARAMETER :: PAD = 150     ! length of each long path component

  TYPE(CRTM_ChannelInfo_type) :: chinfo(1)
  CHARACTER(:), ALLOCATABLE   :: padA, padB, long_path
  INTEGER :: err, nch, estat, cstat

  padA = REPEAT('a', PAD)
  padB = REPEAT('b', PAD)

  ! Build a long directory path that resolves to the real test-data tree:
  !   ./<padA>/<padB>  is a symlink to  ../testinput  (i.e. ./testinput)
  CALL EXECUTE_COMMAND_LINE( 'rm -rf ./'//padA, WAIT=.TRUE. )
  CALL EXECUTE_COMMAND_LINE( 'mkdir -p ./'//padA, WAIT=.TRUE. )
  CALL EXECUTE_COMMAND_LINE( 'ln -sfn ../testinput ./'//padA//'/'//padB, &
                             WAIT=.TRUE., EXITSTAT=estat, CMDSTAT=cstat )
  IF ( cstat /= 0 .OR. estat /= 0 ) THEN
    CALL Display_Message( PROGRAM_NAME, 'could not build the long symlink path', FAILURE )
    STOP 1
  END IF
  long_path = './'//padA//'/'//padB//'/'

  WRITE(*,'(/a)')     ' Long-path CRTM_Init check (issue #238):'
  WRITE(*,'(a,i0)')   '   File_Path length (characters) : ', LEN(long_path)
  IF ( LEN(long_path) <= 256 ) THEN
    CALL Display_Message( PROGRAM_NAME, 'test path is not longer than the old 256 cap', FAILURE )
    CALL Cleanup(); STOP 1
  END IF

  err = CRTM_Init( (/ SENSOR /), chinfo, File_Path=long_path, Quiet=.TRUE. )
  IF ( err /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, &
      'CRTM_Init failed through a '//itoa(LEN(long_path))//'-character path '// &
      '(coefficient path was truncated?)', FAILURE )
    CALL Cleanup(); STOP 1
  END IF

  nch = SUM( CRTM_ChannelInfo_n_Channels(chinfo) )
  WRITE(*,'(a,i0)')   '   channels loaded              : ', nch
  IF ( nch < 1 ) THEN
    CALL Display_Message( PROGRAM_NAME, 'no channels loaded through the long path', FAILURE )
    err = CRTM_Destroy( chinfo ); CALL Cleanup(); STOP 1
  END IF

  err = CRTM_Destroy( chinfo )
  CALL Cleanup()

  WRITE(*,'(/a)') ' PASS: CRTM initialized from a long path without truncation.'
  STOP 0

CONTAINS

  SUBROUTINE Cleanup()
    CALL EXECUTE_COMMAND_LINE( 'rm -rf ./'//padA, WAIT=.TRUE. )
  END SUBROUTINE Cleanup

  PURE FUNCTION itoa(i) RESULT(s)
    INTEGER, INTENT(IN) :: i
    CHARACTER(:), ALLOCATABLE :: s
    CHARACTER(16) :: buf
    WRITE(buf,'(i0)') i
    s = TRIM(buf)
  END FUNCTION itoa

END PROGRAM test_Long_Path_Init
