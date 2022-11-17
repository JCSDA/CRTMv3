!----------------------------------------------------------------------------------
!:sdoc+:
!
!
! FASTEM_Coeff_Reader
!
! Module containing the load/destruction routines to handel
! the shared CSEM microwave water surface emissivity
! model data in NetCDF format. 
!
! PUBLIC DATA:
!   MWwaterC:  Data structure containing the microwave water surface
!              emissivity data.
!
! SIDE EFFECTS:
!       Routines in this module modify the contents of the public
!       data structure MWwaterC.
!
! RESTRICTIONS:
!       Routines in this module should only be called during the
!       CSEM initialisation.
!
! CREATION HISTORY:
!       Written by:     Ming  Chen, 12-Nov-2015
!                       ming.chen@noaa.gov
!

MODULE FASTEM_Coeff_Reader
  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module use
  USE CSEM_Type_Kinds,           ONLY: fp => CSEM_fp
  USE CSEM_Exception_Handler
  USE CRTM_MWwaterCoeff_Define,  ONLY: CRTM_MWwaterCoeff_type, &
                                 CRTM_MWwaterCoeff_Associated, &
                                 CRTM_MWwaterCoeff_Create, &
                                 CRTM_MWwaterCoeff_Destroy

  USE netcdf
  ! Disable implicit typing
  IMPLICIT NONE
 
  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! The shared data
  PUBLIC :: CSEM_MWwaterC
  ! Procedures
  PUBLIC :: CSEM_MWwaterCoeff_Load
  PUBLIC :: CSEM_MWwaterCoeff_CleanUp
  PUBLIC :: CSEM_MWwaterCoeff_IsLoaded
   
  ! -----------------
  ! Module parameters
  ! -----------------
  ! Message string length
  INTEGER, PARAMETER :: ML = 512
  
  
  ! Group name of the target set 
  INTEGER,      PARAMETER :: N_SUBGRP = 6
  INTEGER,      PARAMETER :: MAX_VAR_DIMS = 4
  CHARACTER(*), PARAMETER :: GROUP_NAME         = 'MWwaterCoeff' 
  CHARACTER(10) :: SUB_GROUP_NAME(N_SUBGRP) =  &
                            (/'FCcoeff ', 'FRcoeff ','RCcoeff ', &
                              'AZcoeff ', 'SSCcoeff','LSCcoeff'/)


  ! Dimension names
  CHARACTER(*), PARAMETER :: ANGLE_DIMNAME      = 'n_Angles'
  CHARACTER(*), PARAMETER :: FREQUENCY_DIMNAME  = 'n_Frequencies'
  CHARACTER(*), PARAMETER :: WIND_SPEED_DIMNAME = 'n_Wind_Speeds'

  ! Variable names. Case sensitive.
  CHARACTER(*), PARAMETER :: RELEASE_VARNAME    = 'Release'
  CHARACTER(*), PARAMETER :: VERSION_VARNAME    = 'Version'
  CHARACTER(*), PARAMETER :: DATA_TYPE_VARNAME  = 'Data_Type'
  CHARACTER(*), PARAMETER :: ANGLE_VARNAME      = 'Angle'
  CHARACTER(*), PARAMETER :: FREQUENCY_VARNAME  = 'Frequency'
  CHARACTER(*), PARAMETER :: WIND_SPEED_VARNAME = 'Wind_Speed'
  CHARACTER(*), PARAMETER :: EMISSIVITY_VARNAME = 'Emissivity'


  CHARACTER(*), PARAMETER :: ANGLE_UNITS      = 'degrees from vertical'
  CHARACTER(*), PARAMETER :: FREQUENCY_UNITS  = 'Inverse centimetres (cm^-1)'
  CHARACTER(*), PARAMETER :: WIND_SPEED_UNITS = 'metres per second (m.s^-1)'
  CHARACTER(*), PARAMETER :: EMISSIVITY_UNITS = 'None.'

  ! Variable _FillValue attribute.
  CHARACTER(*), PARAMETER :: FILLVALUE_ATTNAME = '_FillValue'

  INTEGER,      PARAMETER :: IP_FILLVALUE = -1
  REAL(fp),     PARAMETER :: FP_FILLVALUE = -999.0_fp
  
 
 
  ! --------------------------------------------------
  ! Status Control variable 
  ! -------------------------------------------------
  LOGICAL,SAVE :: MWwaterC_Loaded  = .FALSE.
  ! --------------------------------------------------
  ! The shared microwave water surface emissivity data
  ! --------------------------------------------------
  TYPE(CRTM_MWwaterCoeff_type),  SAVE :: CSEM_MWwaterC


CONTAINS


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CSEM_MWwaterCoeff_Load
!
! PURPOSE:
!       Function to load the microwave water surface emissivity data into
!       the public data structure MWwaterC
!
! CALLING SEQUENCE:
!       Error_Status = CSEM_MWwaterCoeff_Load( &
!                        Filename,                              &
!                        File_Path         = File_Path        , &
!                        Quiet             = Quiet            , &
!                        Process_ID        = Process_ID       , &
!                        Output_Process_ID = Output_Process_ID  )
!
! INPUT ARGUMENTS:
!       Filename:           Name of the MWwaterCoeff file.
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
!
!
! OPTIONAL INPUT ARGUMENTS:
!       File_Path:          Character string specifying a file path for the
!                           input data file. If not specified, the current
!                           directory is the default.
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Quiet:              Set this logical argument to suppress INFORMATION
!                           messages being printed to stdout
!                           If == .FALSE., INFORMATION messages are OUTPUT [DEFAULT].
!                              == .TRUE.,  INFORMATION messages are SUPPRESSED.
!                           If not specified, default is .FALSE.
!                           UNITS:      N/A
!                           TYPE:       LOGICAL
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Process_ID:         Set this argument to the MPI process ID that this
!                           function call is running under. This value is used
!                           solely for controlling INFORMATIOn message output.
!                           If MPI is not being used, ignore this argument.
!                           This argument is ignored if the Quiet argument is set.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Output_Process_ID:  Set this argument to the MPI process ID in which
!                           all INFORMATION messages are to be output. If
!                           the passed Process_ID value agrees with this value
!                           the INFORMATION messages are output. 
!                           This argument is ignored if the Quiet argument
!                           is set.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status:       The return value is an integer defining the error
!                           status. The error codes are defined in the
!                           Message_Handler module.
!                           If == SUCCESS the data load was successful
!                              == FAILURE an unrecoverable error occurred.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!
! SIDE EFFECTS:
!       This function modifies the contents of the public data
!       structure MWwaterC.
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CSEM_MWwaterCoeff_Load( &
    Filename         , &  ! Input
    File_Path        , &  ! Optional input
    Quiet            , &  ! Optional input
    Version          , &  ! Optional input
    Process_ID       , &  ! Optional input
    Output_Process_ID) &  ! Optional input
  RESULT( Error_Status )
    ! Arguments
    CHARACTER(*),           INTENT(IN) :: Filename
    CHARACTER(*), OPTIONAL, INTENT(IN) :: File_Path
    LOGICAL     , OPTIONAL, INTENT(IN) :: Quiet             
    INTEGER     , OPTIONAL, INTENT(IN) :: Version             
    INTEGER     , OPTIONAL, INTENT(IN) :: Process_ID
    INTEGER     , OPTIONAL, INTENT(IN) :: Output_Process_ID
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CSEM_MWwaterCoeff_Load'
    ! Local variables
    CHARACTER(ML) :: Message, pid_msg
    CHARACTER(ML) :: MWwaterCoeff_File
    LOGICAL :: noisy
    
    INTEGER :: File_Group_ID(N_SUBGRP+1)
    INTEGER :: ndim_subgrp(N_SUBGRP)
    INTEGER :: dims_subgrp(MAX_VAR_DIMS,N_SUBGRP)
    INTEGER :: varid, i, dim1(1),dim3(3)


    ! Setup 
    Error_Status = SUCCESS
    ! ...Assign the filename to local variable
    MWwaterCoeff_File = ADJUSTL(Filename)
    ! ...Add the file path
    IF ( PRESENT(File_Path) ) MWwaterCoeff_File = TRIM(ADJUSTL(File_Path))//TRIM(MWwaterCoeff_File)
    ! ...Check Quiet argument
    noisy = .TRUE.
    IF ( PRESENT(Quiet) ) noisy = .NOT. Quiet
    ! ...Check the MPI Process Ids
    IF ( noisy .AND. PRESENT(Process_ID) .AND. PRESENT(Output_Process_ID) ) THEN
      IF ( Process_Id /= Output_Process_Id ) noisy = .FALSE.
    END IF
    ! ...Create a process ID message tag for error messages
    IF ( PRESENT(Process_Id) ) THEN
      WRITE( pid_msg,'("; Process ID: ",i0)' ) Process_ID
    ELSE
      pid_msg = ''
    END IF

    IF(MWwaterC_Loaded)  THEN
       PRINT*, 'MWwaterCoeff already loaded, reloading ...'
       Error_Status = CSEM_MWwaterCoeff_CleanUP()
    ENDIF
    Error_Status = INQ_EmisCoeff_File( &
          MWwaterCoeff_File,     &
          ndim_subgrp,           &
          dims_subgrp,           &
          File_Group_ID)  

    ! Allocate the output structure
    CALL  CRTM_MWwaterCoeff_Create(CSEM_MWwaterC,ndim_subgrp,dims_subgrp)
 
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error occurred allocating MWwaterC structure.'
      PRINT*,Message
      Error_Status = check( nf90_close(File_Group_ID(1)) )
      RETURN
    END IF

    ! Read the MWwaterCoeff data file
    IF(PRESENT(Version)) CSEM_MWwaterC%Version = Version
    i = 1 ; dim1 = dims_subgrp(1:1,i)
    Error_Status = check(nf90_inq_varid(File_Group_ID(i+1), "C", varid))
    Error_Status = check(nf90_get_var(File_Group_ID(i+1), varid,  &
              CSEM_MWwaterC%FCCoeff%C,  count=dim1))
    i = 2 ; dim1 =    dims_subgrp(1:1,i)
    Error_Status = check(nf90_inq_varid(File_Group_ID(i+1), "C", varid))      
    Error_Status = check(nf90_get_var(File_Group_ID(i+1),   varid,  &
              CSEM_MWwaterC%FRCoeff%C,  count=dim1))
    i = 3 ; dim3 =    dims_subgrp(1:3,i)      
    Error_Status = check(nf90_inq_varid(File_Group_ID(i+1), "C", varid))           
    Error_Status = check(nf90_get_var(File_Group_ID(i+1),   varid,  &
              CSEM_MWwaterC%RCCoeff%C, count=dim3))
    i = 4 ; dim3 =    dims_subgrp(1:3,i)      
    Error_Status = check(nf90_inq_varid(File_Group_ID(i+1), "C", varid))           
    Error_Status = check(nf90_get_var(File_Group_ID(i+1),   varid,  &
             CSEM_MWwaterC%AZCoeff%C, count=dim3))
    i = 5  ; dim1 = dims_subgrp(1:1,i)       
    Error_Status = check(nf90_inq_varid(File_Group_ID(i+1), "C", varid))           
    Error_Status = check(nf90_get_var(File_Group_ID(i+1),   varid,  &
             CSEM_MWwaterC%SSCCoeff%C, count=dim1))
    i = 6  ; dim3 =    dims_subgrp(1:3,i)          
    Error_Status = check(nf90_inq_varid(File_Group_ID(i+1), "C", varid))           
    Error_Status = check(nf90_get_var(File_Group_ID(i+1),   varid,  &
             CSEM_MWwaterC%LSCCoeff%C, count=dim3))


    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error occurred allocating MWwaterC structure.'
      PRINT*,Message
      Error_Status = CSEM_MWwaterCoeff_CleanUp()    
      Error_Status = check(nf90_close(File_Group_ID(1)) )
      RETURN
    END IF
   
    MWwaterC_Loaded = .TRUE. 
 
 
  END FUNCTION CSEM_MWwaterCoeff_Load


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_MWwaterCoeff_Destroy
!
! PURPOSE:
!       Function to deallocate the public data structure MWwaterC containing
!       the CSEM microwave water surface emissivity data.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_MWwaterCoeff_Destroy( Process_ID = Process_ID )
!
! OPTIONAL INPUTS:
!       Process_ID:       Set this argument to the MPI process ID that this
!                         function call is running under. This value is used
!                         solely for controlling message output. If MPI is not
!                         being used, ignore this argument.
!                         UNITS:      N/A
!                         TYPE:       INTEGER
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status:     The return value is an integer defining the error
!                         status. The error codes are defined in the
!                         Message_Handler module.
!                         If == SUCCESS the deallocation of the public data
!                                       structure was successful
!                            == FAILURE an unrecoverable error occurred.
!                         UNITS:      N/A
!                         TYPE:       INTEGER
!                         DIMENSION:  Scalar
!
! SIDE EFFECTS:
!       This function modifies the contents of the public data
!       structure MWwaterC.
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CSEM_MWwaterCoeff_CleanUp( Process_ID ) RESULT( err_stat )
    ! Arguments
    INTEGER, OPTIONAL, INTENT(IN) :: Process_ID
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_MWwaterCoeff_Destroy'
    ! Local variables
    CHARACTER(ML) :: msg, pid_msg

    ! Setup
    err_stat = SUCCESS
    ! ...Create a process ID message tag for error messages
    IF ( PRESENT(Process_Id) ) THEN
      WRITE( pid_msg,'("; Process ID: ",i0)' ) Process_ID
    ELSE
      pid_msg = ''
    END IF
    
 
    ! Destroy the structure
    IF(MWwaterC_Loaded) CALL CRTM_MWwaterCoeff_Destroy(CSEM_MWwaterC )
    IF ( CRTM_MWwaterCoeff_Associated( CSEM_MWwaterC ) ) THEN
      err_stat = FAILURE
      msg = 'Error deallocating MWwaterCoeff shared data structure'//TRIM(pid_msg)
      CALL Display_Message( ROUTINE_NAME, msg, err_stat ); RETURN
    END IF
    
    MWwaterC_Loaded = .FALSE.

  END FUNCTION CSEM_MWwaterCoeff_CleanUp

  
!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CSEM_MWwaterCoeff_IsLoaded
!
! PURPOSE:
!       Function to test if microwave water surface emissivity data has
!       been loaded into the public data structure MWwaterC.
!
! CALLING SEQUENCE:
!       status = CSEM_MWwaterCoeff_IsLoaded()
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CSEM_MWwaterCoeff_IsLoaded() RESULT( IsLoaded )
    LOGICAL :: IsLoaded
    IsLoaded = MWwaterC_Loaded
  END FUNCTION CSEM_MWwaterCoeff_IsLoaded


!##################################################################################
!##################################################################################
!##                                                                              ##
!##                          ## PRIVATE MODULE ROUTINES ##                       ##
!##                                                                              ##
!##################################################################################
!##################################################################################


  FUNCTION INQ_EmisCoeff_File( &
       FILE_NAME,     &
       ndim_subgrp,   &
       dims_subgrp,   &
       File_Group_ID) & 
    RESULT ( Error_status)
  
    CHARACTER(*), INTENT(IN)  :: FILE_NAME
    INTEGER,      INTENT(OUT) :: ndim_subgrp(N_SUBGRP)
    INTEGER,      INTENT(OUT) :: dims_subgrp(MAX_VAR_DIMS,N_SUBGRP)
    INTEGER,      INTENT(OUT),   OPTIONAL  :: File_Group_ID(N_SUBGRP+1)
    INTEGER   :: Error_status

    INTEGER :: ncid, grpid, varid
    LOGICAL :: Existence
    INTEGER :: i, j, len, ndims, numgrps
    INTEGER :: ncids(N_SUBGRP+1)
    INTEGER :: grpids(N_SUBGRP+1)
    INTEGER :: dimids(MAX_VAR_DIMS)
  
    Error_Status = SUCCESS
    INQUIRE( FILE = TRIM( FILE_NAME), EXIST = Existence )
    IF ( .NOT. Existence ) THEN
      PRINT*,'File '//TRIM( FILE_NAME )//' not found.'
      Error_Status = FAILURE
      RETURN
    END IF

    !Open the NetCDF file:
    Error_Status = check(nf90_open( TRIM( FILE_NAME ), nf90_nowrite, ncid))
    Error_Status = check(nf90_inq_grps(ncid, numgrps, ncids))
    grpids(1) = ncid

    IF(numgrps >= 1) THEN
       Error_Status = check(nf90_inq_grp_ncid(ncid, TRIM(GROUP_NAME), grpid))
    END IF
    IF(Error_Status /= SUCCESS) THEN
       ncid = grpids(1)
       Error_Status = check( nf90_close(ncid) ) 
       PRINT*, 'NC file does not have the target dataset group.....' 
       RETURN    
    END IF
    Error_Status = check(nf90_inq_grps(grpid, numgrps, ncids))
    IF(numgrps /= N_SUBGRP) THEN
       PRINT*, 'Number of sub-groups is not equal to the defined '
    END IF
    DO i = 1, N_SUBGRP
        Error_Status = check(nf90_inq_grp_ncid(grpid, TRIM(SUB_GROUP_NAME(i)), grpids(i+1)))
        IF(Error_Status /= SUCCESS) THEN 
          Error_Status = check( nf90_close(ncid) ) 
          print*, TRIM(SUB_GROUP_NAME(i)) // '  is not the sub-group name...'
          RETURN
        END IF 
        !Error_Status = check(nf90_inq_dimid(grpids(i+1), "dim1",  dimid))
        !Error_Status = check(nf90_inquire_dimension(grpids(i+1), dimid, len=len))
        Error_Status = check(nf90_inq_varid(grpids(i+1), "C", varid))
        Error_Status = check(nf90_inquire_variable(grpids(i+1), varid,  ndims=ndims, dimids=dimids))
        ndim_subgrp(i) = ndims
        DO j = 1, ndims
          Error_Status = check(nf90_inquire_dimension(grpids(i+1), dimids(j), len=len))
          dims_subgrp(j,i) = len
        END DO

    END DO
 
    IF(PRESENT(File_Group_ID)) THEN
      File_Group_ID =  grpids 
      RETURN
    ENDIF
    
    Error_Status = check( nf90_close(ncid) )

  END FUNCTION  INQ_EmisCoeff_File


  
  FUNCTION check(status) RESULT ( Err_Status )

    INTEGER, INTENT ( IN) :: status
    INTEGER :: err_status
    err_status = SUCCESS
    IF (status /= nf90_noerr) THEN
      PRINT *, trim(nf90_strerror(status))
      err_status = FAILURE
    END IF
  END FUNCTION check


END MODULE FASTEM_Coeff_Reader
