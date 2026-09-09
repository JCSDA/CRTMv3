!
! CloudCoeff_Exp_netCDF_IO
!
! netCDF reader (and minimal inquire) for the EXPERIMENTAL CloudCoeff (v1).
! See CloudCoeff_Experimental_Schema_v1.md. Legacy CloudCoeff IO is unchanged.
!
! The CloudCoeff_Exp arrays are declared in canonical Fortran order
! (fastest index first); the file stores the reversed (C/CDL) order, so nf90_get_var
! fills the Fortran arrays directly with no manual transpose.
!

MODULE CloudCoeff_Exp_netCDF_IO

  USE Type_Kinds       , ONLY: Long, Double
  USE Message_Handler  , ONLY: SUCCESS, FAILURE, INFORMATION, Display_Message
  USE File_Utility     , ONLY: File_Exists
  USE String_Utility   , ONLY: StrClean
  USE CloudCoeff_Exp_Define, ONLY: CloudCoeff_Exp_type, &
                                   CloudCoeff_Exp_Create, &
                                   CloudCoeff_Exp_Destroy, &
                                   CloudCoeff_Exp_Associated
  USE netcdf
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: CloudCoeff_Exp_netCDF_ReadFile
  PUBLIC :: CloudCoeff_Exp_netCDF_InquireFile

  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: CloudCoeff_Exp_netCDF_IO.f90  v1  2026-06-01 $'
  INTEGER, PARAMETER :: ML = 512
  INTEGER, PARAMETER :: SL = 32      ! must match CloudCoeff_Exp_Define SL / file nchar

CONTAINS

  !----------------------------------------------------------------------------
  ! Read the experimental CloudCoeff netCDF into a CloudCoeff_Exp_type structure.
  !----------------------------------------------------------------------------
  FUNCTION CloudCoeff_Exp_netCDF_ReadFile( Filename, CloudCoeff_Exp, Quiet ) RESULT( err_stat )
    CHARACTER(*),              INTENT(IN)  :: Filename
    TYPE(CloudCoeff_Exp_type), INTENT(OUT) :: CloudCoeff_Exp
    LOGICAL,         OPTIONAL, INTENT(IN)  :: Quiet
    INTEGER :: err_stat
    ! Local
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CloudCoeff_Exp_netCDF_ReadFile'
    CHARACTER(ML) :: msg
    LOGICAL :: noisy, ok
    INTEGER :: fid, vid, did
    INTEGER :: nF, nT, nMu, nDm, nH, nL, nP
    INTEGER :: i

    err_stat = SUCCESS
    noisy = .TRUE.; IF ( PRESENT(Quiet) ) noisy = .NOT. Quiet

    IF ( .NOT. File_Exists(Filename) ) THEN
      CALL Display_Message( ROUTINE_NAME, 'File not found: '//TRIM(Filename), FAILURE )
      err_stat = FAILURE; RETURN
    END IF

    ok = chk( nf90_open( Filename, NF90_NOWRITE, fid ), 'open '//TRIM(Filename) )
    IF ( .NOT. ok ) THEN; err_stat = FAILURE; RETURN; END IF

    ! --- dimensions ---
    nF=0; nT=0; nMu=0; nDm=0; nH=0; nL=0; nP=0
    CALL get_dim( 'n_Frequency'     , nF  )
    CALL get_dim( 'n_Temperature'   , nT  )
    CALL get_dim( 'n_Mu'            , nMu )
    CALL get_dim( 'n_Dm'            , nDm )
    CALL get_dim( 'n_Habit'         , nH  )
    CALL get_dim( 'n_Legendre'      , nL  )
    CALL get_dim( 'n_Phase_Elements', nP  )
    IF ( .NOT. ok ) GO TO 900

    CALL CloudCoeff_Exp_Create( CloudCoeff_Exp, nF, nT, nMu, nDm, nH, nL, nP )
    IF ( .NOT. CloudCoeff_Exp_Associated( CloudCoeff_Exp ) ) THEN
      msg = 'Allocation failed'; ok = .FALSE.; GO TO 900
    END IF

    ! --- global attributes (best-effort) ---
    IF ( nf90_get_att( fid, NF90_GLOBAL, 'Scheme',  CloudCoeff_Exp%Scheme  ) == NF90_NOERR ) &
      CALL StrClean( CloudCoeff_Exp%Scheme )
    i = nf90_get_att( fid, NF90_GLOBAL, 'Release', CloudCoeff_Exp%Release )
    i = nf90_get_att( fid, NF90_GLOBAL, 'Version', CloudCoeff_Exp%Version )

    ! --- axes ---
    CALL get_var_d( 'Frequency'  , CloudCoeff_Exp%Frequency   )
    CALL get_var_d( 'Temperature', CloudCoeff_Exp%Temperature )
    CALL get_var_d( 'Mu'         , CloudCoeff_Exp%Mu          )
    CALL get_var_d( 'Dm'         , CloudCoeff_Exp%Dm          )
    ! --- habit metadata ---
    CALL get_var_i( 'Habit_Id'   , CloudCoeff_Exp%Habit_Id    )
    CALL get_var_i( 'Habit_Phase', CloudCoeff_Exp%Habit_Phase )
    CALL get_var_d( 'mD_a'       , CloudCoeff_Exp%mD_a        )
    CALL get_var_d( 'mD_b'       , CloudCoeff_Exp%mD_b        )
    ! Optional per-habit Effective_Radius -> Dm multiplier. Absent in older LUTs ->
    ! stays at the 1.0 identity set by CloudCoeff_Exp_Create (no behaviour change).
    IF ( ok ) THEN
      IF ( nf90_inq_varid( fid, 'Reff_to_Dm', vid ) == NF90_NOERR ) &
        ok = chk( nf90_get_var( fid, vid, CloudCoeff_Exp%Reff_to_Dm ), 'get Reff_to_Dm' )
    END IF
    ! Habit_Name : char(n_Habit, nchar) -> CHARACTER(SL) array
    IF ( ok ) ok = chk( nf90_inq_varid( fid, 'Habit_Name', vid ), 'varid Habit_Name' )
    IF ( ok ) ok = chk( nf90_get_var( fid, vid, CloudCoeff_Exp%Habit_Name ), 'get Habit_Name' )
    IF ( ok ) THEN
      DO i = 1, nH; CALL StrClean( CloudCoeff_Exp%Habit_Name(i) ); END DO
    END IF
    ! --- bulk optics + truncation ---
    CALL get_var_d5( 'ke', CloudCoeff_Exp%ke )
    CALL get_var_d5( 'ka', CloudCoeff_Exp%ka )
    CALL get_var_d5( 'g' , CloudCoeff_Exp%g  )
    CALL get_var_d5( 'kb', CloudCoeff_Exp%kb )
    IF ( ok ) ok = chk( nf90_inq_varid( fid, 'n_Legendre_Eff', vid ), 'varid n_Legendre_Eff' )
    IF ( ok ) ok = chk( nf90_get_var( fid, vid, CloudCoeff_Exp%n_Legendre_Eff ), 'get n_Legendre_Eff' )
    ! --- phase expansion ---
    IF ( ok ) ok = chk( nf90_inq_varid( fid, 'pcoeff', vid ), 'varid pcoeff' )
    IF ( ok ) ok = chk( nf90_get_var( fid, vid, CloudCoeff_Exp%pcoeff ), 'get pcoeff' )

900 CONTINUE
    i = nf90_close( fid )
    IF ( .NOT. ok ) THEN
      CALL Display_Message( ROUTINE_NAME, 'Read failed: '//TRIM(msg), FAILURE )
      CALL CloudCoeff_Exp_Destroy( CloudCoeff_Exp )
      err_stat = FAILURE; RETURN
    END IF
    IF ( noisy ) THEN
      WRITE(msg,'("CloudCoeff_Exp read OK: ",a," nFreq=",i0," nDm=",i0," nMu=",i0,&
                &" nT=",i0," nHabit=",i0," Lmax=",i0," nPhase=",i0)') &
        TRIM(CloudCoeff_Exp%Scheme), nF, nDm, nMu, nT, nH, nL, nP
      CALL Display_Message( ROUTINE_NAME, TRIM(msg), INFORMATION )
    END IF

  CONTAINS

    LOGICAL FUNCTION chk( status, what )
      INTEGER,      INTENT(IN) :: status
      CHARACTER(*), INTENT(IN) :: what
      chk = ( status == NF90_NOERR )
      IF ( .NOT. chk ) msg = TRIM(what)//': '//TRIM(nf90_strerror(status))
    END FUNCTION chk

    SUBROUTINE get_dim( name, val )
      CHARACTER(*), INTENT(IN)  :: name
      INTEGER,      INTENT(OUT) :: val
      val = 0
      IF ( .NOT. ok ) RETURN
      ok = chk( nf90_inq_dimid( fid, name, did ), 'dimid '//name )
      IF ( ok ) ok = chk( nf90_inquire_dimension( fid, did, len=val ), 'dimlen '//name )
    END SUBROUTINE get_dim

    SUBROUTINE get_var_d( name, arr )
      CHARACTER(*), INTENT(IN)  :: name
      REAL(Double), INTENT(OUT) :: arr(:)
      IF ( .NOT. ok ) RETURN
      ok = chk( nf90_inq_varid( fid, name, vid ), 'varid '//name )
      IF ( ok ) ok = chk( nf90_get_var( fid, vid, arr ), 'get '//name )
    END SUBROUTINE get_var_d

    SUBROUTINE get_var_i( name, arr )
      CHARACTER(*),  INTENT(IN)  :: name
      INTEGER(Long), INTENT(OUT) :: arr(:)
      IF ( .NOT. ok ) RETURN
      ok = chk( nf90_inq_varid( fid, name, vid ), 'varid '//name )
      IF ( ok ) ok = chk( nf90_get_var( fid, vid, arr ), 'get '//name )
    END SUBROUTINE get_var_i

    SUBROUTINE get_var_d5( name, arr )
      CHARACTER(*), INTENT(IN)  :: name
      REAL(Double), INTENT(OUT) :: arr(:,:,:,:,:)
      IF ( .NOT. ok ) RETURN
      ok = chk( nf90_inq_varid( fid, name, vid ), 'varid '//name )
      IF ( ok ) ok = chk( nf90_get_var( fid, vid, arr ), 'get '//name )
    END SUBROUTINE get_var_d5

  END FUNCTION CloudCoeff_Exp_netCDF_ReadFile


  !----------------------------------------------------------------------------
  ! Minimal inquire: dimensions + scheme string.
  !----------------------------------------------------------------------------
  FUNCTION CloudCoeff_Exp_netCDF_InquireFile( Filename, n_Frequency, n_Temperature, &
       n_Mu, n_Dm, n_Habit, n_Legendre, n_Phase_Elements, Scheme ) RESULT( err_stat )
    CHARACTER(*),           INTENT(IN)  :: Filename
    INTEGER,      OPTIONAL, INTENT(OUT) :: n_Frequency, n_Temperature, n_Mu, n_Dm, &
                                           n_Habit, n_Legendre, n_Phase_Elements
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: Scheme
    INTEGER :: err_stat
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CloudCoeff_Exp_netCDF_InquireFile'
    INTEGER :: fid, did, i, v
    LOGICAL :: ok

    err_stat = SUCCESS
    IF ( .NOT. File_Exists(Filename) ) THEN
      CALL Display_Message( ROUTINE_NAME, 'File not found: '//TRIM(Filename), FAILURE )
      err_stat = FAILURE; RETURN
    END IF
    ok = ( nf90_open( Filename, NF90_NOWRITE, fid ) == NF90_NOERR )
    IF ( .NOT. ok ) THEN; err_stat = FAILURE; RETURN; END IF
    CALL gd( 'n_Frequency'     , n_Frequency      )
    CALL gd( 'n_Temperature'   , n_Temperature    )
    CALL gd( 'n_Mu'            , n_Mu             )
    CALL gd( 'n_Dm'            , n_Dm             )
    CALL gd( 'n_Habit'         , n_Habit          )
    CALL gd( 'n_Legendre'      , n_Legendre       )
    CALL gd( 'n_Phase_Elements', n_Phase_Elements )
    IF ( PRESENT(Scheme) ) THEN
      Scheme = ' '
      IF ( nf90_get_att( fid, NF90_GLOBAL, 'Scheme', Scheme ) == NF90_NOERR ) CALL StrClean( Scheme )
    END IF
    i = nf90_close( fid )
    IF ( .NOT. ok ) err_stat = FAILURE
  CONTAINS
    SUBROUTINE gd( name, out )
      CHARACTER(*),      INTENT(IN)  :: name
      INTEGER, OPTIONAL, INTENT(OUT) :: out
      IF ( .NOT. PRESENT(out) ) RETURN
      out = 0
      IF ( ok .AND. nf90_inq_dimid( fid, name, did ) == NF90_NOERR ) THEN
        IF ( nf90_inquire_dimension( fid, did, len=v ) == NF90_NOERR ) out = v
      END IF
    END SUBROUTINE gd
  END FUNCTION CloudCoeff_Exp_netCDF_InquireFile

END MODULE CloudCoeff_Exp_netCDF_IO
