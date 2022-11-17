!
! IRSSEM_EmisCoeff_V2_Reader
!
! Module containing routines to read the netCDF format EmisCoeff V2 files 
! of the NESDIS physical Infrared ocean surface models
!       
!
! CREATION HISTORY:
!       Written by:     Ming Chen Apr-08-2020
!                       ming.chen@noaa.gov
!

MODULE IRSSEM_EmisCoeff_V2_Reader

  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module use
  USE CSEM_Type_Kinds,         ONLY: fp => CSEM_fp
  USE CSEM_Exception_Handler
  USE IRSSEM_EmisCoeff_V2_Define, ONLY:                    &
                               EmisCoeff_V2_Type,          &
                               SPECTRAL_EMISCOEFF_V2_TYPE, &
                               SENSOR_EMISCOEFF_V2_TYPE,   &
                               Associated_EmisCoeff_V2,    &
                               Allocate_EmisCoeff_V2,      &
                               Destroy_EmisCoeff_V2

  USE netcdf
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  PUBLIC :: Load_IRSSEM_V2_LUT
  PUBLIC :: Close_IRSSEM_V2_LUT
  PUBLIC :: IRwaterC_V2
  
  ! -----------------
  ! Module parameters
  ! -----------------
    
  ! Group name of the target set 
  CHARACTER(*), PARAMETER :: GROUP_NAME         = 'IRwaterCoeff' 

  ! Dimension names
  CHARACTER(*), PARAMETER :: ANGLE_DIMNAME      = 'n_Angles'
  CHARACTER(*), PARAMETER :: FREQUENCY_DIMNAME  = 'n_Frequencies'
  CHARACTER(*), PARAMETER :: WIND_SPEED_DIMNAME = 'n_Wind_Speeds'
  CHARACTER(*), PARAMETER :: TEMPERATURE_DIMNAME = 'n_Temperature'

  ! Variable names. Case sensitive.
  CHARACTER(*), PARAMETER :: RELEASE_VARNAME    = 'Release'
  CHARACTER(*), PARAMETER :: VERSION_VARNAME    = 'Version'
  CHARACTER(*), PARAMETER :: DATA_TYPE_VARNAME  = 'Data_Type'
  CHARACTER(*), PARAMETER :: ANGLE_VARNAME      = 'Angle'
  CHARACTER(*), PARAMETER :: FREQUENCY_VARNAME  = 'Frequency'
  CHARACTER(*), PARAMETER :: WIND_SPEED_VARNAME = 'Wind_Speed'
  CHARACTER(*), PARAMETER :: TEMPERATURE_VARNAME = 'Temperature'
  CHARACTER(*), PARAMETER :: EMISSIVITY_VARNAME = 'Emissivity'

! Variable units attribute.
  CHARACTER(*), PARAMETER :: UNITS_ATTNAME      = 'units'
  CHARACTER(*), PARAMETER :: ANGLE_UNITS        = 'degrees from vertical'
  CHARACTER(*), PARAMETER :: FREQUENCY_UNITS    = 'Inverse centimetres (cm^-1)'
  CHARACTER(*), PARAMETER :: WIND_SPEED_UNITS   = 'metres per second (m.s^-1)'
  CHARACTER(*), PARAMETER :: TEMPERATURE_UNITS   = 'degrees (k)'
  CHARACTER(*), PARAMETER :: EMISSIVITY_UNITS   = 'None.'

  ! Variable _FillValue attribute.
  CHARACTER(*), PARAMETER :: FILLVALUE_ATTNAME = '_FillValue'

  INTEGER,      PARAMETER :: IP_FILLVALUE = -1
  REAL(fp),     PARAMETER :: FP_FILLVALUE = -999.0_fp

  
  TYPE(EmisCoeff_V2_type), SAVE :: IRwaterC_V2
  LOGICAL, SAVE :: IRwaterC_Loaded =.FALSE.

CONTAINS

!--------------------------------------------------------------------------------
!
! NAME:
!       Load_IRSSEM_V2_LUT
!
! PURPOSE:
!       Function to read data from a netCDF format Spectral EmisCoeff file.
!
! CALLING SEQUENCE:
!       Error_Status = Load_IRSSEM_V2_LUT( NC_Filename)                     ! Input
!
! INPUT ARGUMENTS:
!       NC_Filename:     Character string specifying the name of the netCDF EmisCoeff
!                        format data file.
!                        UNITS:      N/A
!                        TYPE:       CHARACTER(*)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!
! FUNCTION RESULT:
!       Error_Status:    The return value is an integer defining the error status.
!                        The error codes are defined in the Message_Handler module.
!                        If == SUCCESS the netCDF file read was successful
!                           == FAILURE an unrecoverable read error occurred.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!
! SIDE EFFECTS:
!
! COMMENTS:
!
!--------------------------------------------------------------------------------

  FUNCTION Load_IRSSEM_V2_LUT( NC_Filename)     &  ! input
    RESULT ( Error_Status )
    ! Arguments
    CHARACTER(*),           INTENT(IN)     :: NC_Filename
   ! Function result
    INTEGER :: Error_Status
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Read_EmisCoeff_netCDF'
    ! Function variables
    CHARACTER(256) :: Message
    INTEGER :: NC_FileID(2)
    INTEGER :: NC_Frequency_ID
    INTEGER :: NC_Angle_ID
    INTEGER :: NC_WindSpeed_ID
    INTEGER :: NC_Temperature_ID
    INTEGER :: NC_Emissivity_ID
    INTEGER :: n_Frequencies
    INTEGER :: n_Wind_Speeds
    INTEGER :: n_Temperatures
    INTEGER :: n_Angles
    INTEGER :: i,j,k
    REAL(fp), ALLOCATABLE :: ee(:,:,:,:)
    ! Set up
    Error_Status = SUCCESS
    IF(IRwaterC_Loaded)  Error_Status = Close_IRSSEM_V2_LUT()

    CALL INQ_EmisCoeff_V2_File(NC_Filename,n_Frequencies,n_Angles,n_Wind_Speeds,n_Temperatures,NC_FileID)
     
    ! Check the structure data type
    IF ( IRwaterC_V2%Data_Type /= SPECTRAL_EMISCOEFF_V2_TYPE) THEN
      Message = 'EmisCoeff structure data type invalid'
      PRINT*,Message
      Error_Status = check( nf90_close(NC_FileID(1)) )
      RETURN
    END IF

    ! Allocate the output structure
    Error_Status = Allocate_EmisCoeff_V2( n_Angles, &
                                       n_Frequencies, &
                                       n_Wind_Speeds, &
                                       n_Temperatures, &
                                       IRwaterC_V2)

   
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error occurred allocating IRwaterC structure.'
      PRINT*,Message
      Error_Status = check( nf90_close(NC_FileID(1)) )
      RETURN
    END IF
    IRwaterC_Loaded = .TRUE.    
    Error_Status = check(nf90_inq_varid(NC_FileID(2), ANGLE_VARNAME,      NC_Angle_ID))
    Error_Status = check(nf90_inq_varid(NC_FileID(2), FREQUENCY_VARNAME,  NC_Frequency_ID))
    Error_Status = check(nf90_inq_varid(NC_FileID(2), WIND_SPEED_VARNAME, NC_WindSpeed_ID))
    Error_Status = check(nf90_inq_varid(NC_FileID(2), TEMPERATURE_VARNAME, NC_TEMPERATURE_ID))
    Error_Status = check(nf90_inq_varid(NC_FileID(2), EMISSIVITY_VARNAME, NC_Emissivity_ID))
   
    ALLOCATE(ee(n_Temperatures, n_Wind_Speeds, n_Frequencies, n_Angles))
    Error_Status = check(nf90_get_var( NC_FileID(2),   NC_Angle_ID,     &
        IRwaterC_V2%Angle,           count=(/n_Angles/)))
    Error_Status = check(nf90_get_var(NC_FileID(2),    NC_Frequency_ID, &
        IRwaterC_V2%Frequency,   count=(/n_Frequencies/)))
    Error_Status = check(nf90_get_var(NC_FileID(2),    NC_Windspeed_ID, &
        IRwaterC_V2%Wind_Speed,  count=(/n_Wind_Speeds/)))
    Error_Status = check(nf90_get_var(NC_FileID(2),    NC_Temperature_ID, &
        IRwaterC_V2%Temperature,  count=(/n_Temperatures/)))
    !Error_Status = check(nf90_get_var(NC_FileID(2),    NC_Emissivity_ID,&
    !    IRwaterC_V2%Emissivity,  count=(/n_Angles, n_Frequencies, n_Wind_Speeds, n_Temperatures/)))
    Error_Status = check(nf90_get_var(NC_FileID(2),    NC_Emissivity_ID,&
        ee,  count=(/n_Temperatures, n_Wind_Speeds, n_Frequencies, n_Angles/)))
    DO i = 1, n_Temperatures
      DO j = 1, n_Wind_Speeds
        DO k = 1, n_Angles
          IRwaterC_V2%Emissivity(k,:,j,i) = ee(i,j,:,k)
        END DO
      END DO
    END DO
  
    DEALLOCATE(ee)
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error occurred allocating IRwaterC structure.'
      PRINT*,Message
      Error_Status = close_irssem_V2_lut()      
      Error_Status = check(nf90_close(NC_FileID(1)) )
      RETURN
    END IF

  END FUNCTION Load_IRSSEM_V2_LUT
  

  FUNCTION Close_IRSSEM_V2_LUT() RESULT ( Err_Status )
    INTEGER :: err_status
    PRINT*, 'Destroy IRSSEM LUT ..'
    err_status =1
    IF(IRwaterC_Loaded) err_status = Destroy_EmisCoeff_V2( IRwaterC_V2)
    IRwaterC_Loaded = .FALSE.
  END FUNCTION Close_IRSSEM_V2_LUT

!##################################################################################
!##################################################################################
!##                                                                              ##
!##                          ## PRIVATE MODULE ROUTINES ##                       ##
!##                                                                              ##
!##################################################################################
!##################################################################################


  SUBROUTINE INQ_EmisCoeff_V2_File(FILE_NAME,n_Frequencies,n_Angles,n_Wind_Speeds,n_Temperatures,File_Group_ID)
  
    CHARACTER(*),  INTENT(IN) :: FILE_NAME
    INTEGER,  INTENT(OUT)  :: n_Angles
    INTEGER,  INTENT(OUT)  :: n_Frequencies
    INTEGER,  INTENT(OUT)  :: n_Wind_Speeds
    INTEGER,  INTENT(OUT)  :: n_Temperatures
    INTEGER,  INTENT(OUT), OPTIONAL:: File_Group_ID(2)
  
    INTEGER :: ncid
    LOGICAL :: Existence
    INTEGER :: Angle_dimID
    INTEGER :: Frequency_dimID
    INTEGER :: Wind_Speed_dimID
    INTEGER :: Temperature_dimID
    INTEGER :: Error_status
    INTEGER :: numgrps, grpid
    INTEGER :: ncids(3)
  
  
    INQUIRE( FILE = TRIM( FILE_NAME), EXIST = Existence )
    IF ( .NOT. Existence ) THEN
      PRINT*,'File '//TRIM( FILE_NAME )//' not found.'
      n_Angles = -1; n_Frequencies = -1; n_Wind_Speeds = -1
      RETURN
    END IF

    !Open the NetCDF file:
    Error_Status = check(nf90_open( TRIM( FILE_NAME ), nf90_nowrite, ncid))
    Error_Status = check(nf90_inq_grps(ncid, numgrps, ncids))
    grpid = ncid
    IF(numgrps >= 1) THEN
       Error_Status = check(nf90_inq_grp_ncid(ncid, TRIM(GROUP_NAME), grpid))
    END IF

    Error_Status = check(nf90_inq_dimid(grpid, ANGLE_DIMNAME,  Angle_dimID))
    Error_Status = check(nf90_inq_dimid(grpid, FREQUENCY_DIMNAME, Frequency_dimID))
    Error_Status = check(nf90_inq_dimid(grpid, WIND_SPEED_DIMNAME,Wind_Speed_dimID))
    Error_Status = check(nf90_inq_dimid(grpid, TEMPERATURE_DIMNAME,Temperature_dimID))
   
    Error_Status = check(nf90_inquire_DIMENSION(grpid, Angle_dimID,  len = n_Angles))
    Error_Status = check(nf90_inquire_DIMENSION(grpid, Frequency_dimID,  len = n_Frequencies))
    Error_Status = check(nf90_inquire_DIMENSION(grpid, Wind_Speed_dimID, len = n_Wind_Speeds))
    Error_Status = check(nf90_inquire_DIMENSION(grpid, Temperature_dimID, len = n_Temperatures))
 
    IF(PRESENT(File_Group_ID)) THEN
      File_Group_ID = (/ ncid, grpid /)
      RETURN
    ENDIF
    
    Error_Status = check( nf90_close(ncid) )

  END SUBROUTINE INQ_EmisCoeff_V2_File


  
  FUNCTION check(status) RESULT ( Err_Status )

    INTEGER, INTENT ( IN) :: status
    INTEGER :: err_status
    err_status = SUCCESS
    IF (status /= nf90_noerr) THEN
      PRINT *, trim(nf90_strerror(status))
      err_status = FAILURE
    END IF
  END FUNCTION check

END MODULE IRSSEM_EmisCoeff_V2_Reader
