!
! TauCoeff_Inspect
!
! Program to inspect the contents of a CRTM TauCoeff file (binary or netCDF)
!
!
! CREATION HISTORY:
!       Written by:     Gemini Assistant, 12-Jan-2026
!                       modeled after SpcCoeff_Inspect
!

PROGRAM TauCoeff_Inspect

  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module usage
  USE File_Utility      , ONLY: File_Exists
  USE Message_Handler   , ONLY: SUCCESS, FAILURE, WARNING, INFORMATION, &
                                Program_Message, Display_Message
  USE Type_Kinds        , ONLY: Long
  USE ODPS_Define       , ONLY: ODPS_type, Destroy_ODPS, &
                                Info_ODPS, Equal_ODPS, Associated_ODPS, &
                                SENSOR_TYPE_NAME, ODPS_ALGORITHM
  USE ODPS_Binary_IO    , ONLY: Read_ODPS_Binary
  USE ODPS_netCDF_IO    , ONLY: Read_ODPS_netCDF
  USE ODAS_Define       , ONLY: ODAS_type, Destroy_ODAS, Info_ODAS, &
                                ODAS_ALGORITHM, Associated_ODAS
  USE ODAS_Binary_IO    , ONLY: Read_ODAS_Binary
  USE ODAS_netCDF_IO    , ONLY: Read_ODAS_netCDF
  USE Binary_File_Utility, ONLY: Open_Binary_File
  USE NETCDF, ONLY: NF90_NOWRITE, NF90_GLOBAL, NF90_OPEN, NF90_GET_ATT, NF90_CLOSE, &
                    NF90_NOERR
  
  ! Disable implicit typing
  IMPLICIT NONE

  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'TauCoeff_Inspect'
  CHARACTER(*), PARAMETER :: PROGRAM_VERSION_ID = '1.1'
  INTEGER,      PARAMETER :: ML = 256

  ! ---------
  ! Variables
  ! ---------
  INTEGER :: err_stat
  CHARACTER(ML) :: filename, msg, info_msg
  INTEGER :: n_args
  TYPE(ODPS_type) :: odps
  TYPE(ODAS_type) :: odas
  LOGICAL :: is_nc, is_bin
  INTEGER :: Algorithm_ID

  ! Output program header
  CALL Program_Message( PROGRAM_NAME, &
                        'Program to display the contents of a CRTM '//&
                        'Binary/netCDF format TauCoeff file to stdout.', &
                        PROGRAM_VERSION_ID )

  ! Get the filename
  n_args = COMMAND_ARGUMENT_COUNT() 
  IF ( n_args > 0 ) THEN 
    CALL GET_COMMAND_ARGUMENT(1, filename) 
  ELSE 
    WRITE( *,FMT='(/5x,"Enter the TauCoeff filename: ")',ADVANCE='NO' ) 
    READ( *,'(a)' ) filename 
  END IF
  filename = ADJUSTL(filename)
  IF ( .NOT. File_Exists( TRIM(filename) ) ) THEN
    msg = 'File '//TRIM(filename)//' not found.'
    CALL Display_Message( PROGRAM_NAME, msg, FAILURE ); STOP
  END IF

  ! Check file extension
  is_nc = (INDEX(TRIM(filename), '.nc') == LEN_TRIM(filename) - 2)
  is_bin = (INDEX(TRIM(filename), '.bin') == LEN_TRIM(filename) - 3)
  
  ! Peek at the algorithm ID
  Algorithm_ID = -1
  IF (is_bin) THEN
     err_stat = Get_AlgorithmID_Bin( TRIM(filename), Algorithm_ID )
  ELSE IF (is_nc) THEN
     err_stat = Get_AlgorithmID_NC( TRIM(filename), Algorithm_ID )
  END IF

  IF ( Algorithm_ID == ODPS_ALGORITHM ) THEN
     IF (is_bin) err_stat = Read_ODPS_Binary( TRIM(filename), odps, Quiet=1 )
     IF (is_nc)  err_stat = Read_ODPS_netCDF( TRIM(filename), odps, Quiet=1 )
     IF ( err_stat == SUCCESS ) THEN
        CALL ODPS_Inspect( odps )
        err_stat = Destroy_ODPS( odps )
        STOP
     END IF
  ELSE IF ( Algorithm_ID == ODAS_ALGORITHM ) THEN
     IF (is_bin) err_stat = Read_ODAS_Binary( TRIM(filename), odas, Quiet=1 )
     IF (is_nc)  err_stat = Read_ODAS_netCDF( TRIM(filename), odas, Quiet=1 )
     IF ( err_stat == SUCCESS ) THEN
        CALL ODAS_Inspect( odas )
        err_stat = Destroy_ODAS( odas )
        STOP
     END IF
  END IF

  ! Fallback if algorithm not identified or read failed
  msg = 'Error reading TauCoeff file '//TRIM(filename)//'. Unknown algorithm or format.'
  CALL Display_Message( PROGRAM_NAME, msg, FAILURE )

CONTAINS

  FUNCTION Get_AlgorithmID_Bin( Filename, Algorithm_ID ) RESULT( Error_Status )
    CHARACTER(*), INTENT(IN)  :: Filename
    INTEGER,      INTENT(OUT) :: Algorithm_ID
    INTEGER :: Error_Status
    INTEGER :: FileID, IO_Status
    INTEGER(Long) :: Dummy, Alg

    Error_Status = FAILURE
    IF ( Open_Binary_File( Filename, FileID ) /= SUCCESS ) RETURN
    READ( FileID, IOSTAT=IO_Status ) Dummy, Dummy ! Release, Version
    IF ( IO_Status == 0 ) THEN
       READ( FileID, IOSTAT=IO_Status ) Alg
       IF ( IO_Status == 0 ) THEN
          Algorithm_ID = Alg
          Error_Status = SUCCESS
       END IF
    END IF
    CLOSE( FileID )
  END FUNCTION Get_AlgorithmID_Bin

  FUNCTION Get_AlgorithmID_NC( Filename, Algorithm_ID ) RESULT( Error_Status )
    CHARACTER(*), INTENT(IN)  :: Filename
    INTEGER,      INTENT(OUT) :: Algorithm_ID
    INTEGER :: Error_Status
    INTEGER :: ncid, IO_Status, Alg

    Error_Status = FAILURE
    IF ( NF90_OPEN(Filename, NF90_NOWRITE, ncid) /= NF90_NOERR ) RETURN
    IF ( NF90_GET_ATT(ncid, NF90_GLOBAL, "Algorithm", Alg) == NF90_NOERR ) THEN
       Algorithm_ID = Alg
       Error_Status = SUCCESS
    END IF
    IO_Status = NF90_CLOSE(ncid)
  END FUNCTION Get_AlgorithmID_NC

  SUBROUTINE ODPS_Inspect( ODPS )
    TYPE(ODPS_type), INTENT(IN) :: ODPS
    
    WRITE(*,'(1x,"ODPS OBJECT")')
    ! Release/version info
    WRITE(*,'(3x,"Release.Version  :",1x,i0,".",i0)') ODPS%Release, ODPS%Version
    ! Dimensions
    WRITE(*,'(3x,"n_Layers         :",1x,i0)') ODPS%n_Layers
    WRITE(*,'(3x,"n_Components     :",1x,i0)') ODPS%n_Components
    WRITE(*,'(3x,"n_Absorbers      :",1x,i0)') ODPS%n_Absorbers
    WRITE(*,'(3x,"n_Channels       :",1x,i0)') ODPS%n_Channels
    WRITE(*,'(3x,"n_Coeffs         :",1x,i0)') ODPS%n_Coeffs
    
    IF ( .NOT. Associated_ODPS(ODPS) ) RETURN
    
    ! Scalar info
    WRITE(*,'(3x,"Sensor_Id        :",1x,a )') TRIM(ODPS%Sensor_Id)
    WRITE(*,'(3x,"WMO_Satellite_ID :",1x,i0)') ODPS%WMO_Satellite_ID 
    WRITE(*,'(3x,"WMO_Sensor_ID    :",1x,i0)') ODPS%WMO_Sensor_ID
    WRITE(*,'(3x,"Sensor_Type      :",1x,a )') TRIM(SENSOR_TYPE_NAME(ODPS%Sensor_Type))
    WRITE(*,'(3x,"Group_Index      :",1x,i0)') ODPS%Group_Index
    
    ! OPTRAN info
    IF ( ODPS%n_OCoeffs > 0 ) THEN
       WRITE(*,'(3x,"OPTRAN Data      :")')
       WRITE(*,'(5x,"n_OCoeffs        :",1x,i0)') ODPS%n_OCoeffs
       WRITE(*,'(5x,"Alpha            :",1x,es13.6)') ODPS%Alpha
       WRITE(*,'(5x,"Alpha_C1         :",1x,es13.6)') ODPS%Alpha_C1
       WRITE(*,'(5x,"Alpha_C2         :",1x,es13.6)') ODPS%Alpha_C2
       WRITE(*,'(5x,"OComponent_Index :",1x,i0)') ODPS%OComponent_Index
    END IF

    ! Arrays summary
    WRITE(*,'(3x,"Sensor_Channel   :")')
    WRITE(*,'(10(1x,i5,:))') ODPS%Sensor_Channel
    
    WRITE(*,'(3x,"Component_ID     :")')
    WRITE(*,'(10(1x,i5,:))') ODPS%Component_ID

    WRITE(*,'(3x,"Absorber_ID      :")')
    WRITE(*,'(10(1x,i5,:))') ODPS%Absorber_ID

  END SUBROUTINE ODPS_Inspect

  SUBROUTINE ODAS_Inspect( ODAS )
    TYPE(ODAS_type), INTENT(IN) :: ODAS

    WRITE(*,'(1x,"ODAS OBJECT")')
    ! Release/version info
    WRITE(*,'(3x,"Release.Version  :",1x,i0,".",i0)') ODAS%Release, ODAS%Version
    ! Dimensions
    WRITE(*,'(3x,"n_Predictors     :",1x,i0)') ODAS%n_Predictors
    WRITE(*,'(3x,"n_Absorbers      :",1x,i0)') ODAS%n_Absorbers
    WRITE(*,'(3x,"n_Channels       :",1x,i0)') ODAS%n_Channels
    WRITE(*,'(3x,"n_Alphas         :",1x,i0)') ODAS%n_Alphas
    WRITE(*,'(3x,"n_Coeffs         :",1x,i0)') ODAS%n_Coeffs

    IF ( .NOT. Associated_ODAS(ODAS) ) RETURN

    ! Scalar info
    WRITE(*,'(3x,"Sensor_Id        :",1x,a )') TRIM(ODAS%Sensor_Id)
    WRITE(*,'(3x,"WMO_Satellite_ID :",1x,i0)') ODAS%WMO_Satellite_ID
    WRITE(*,'(3x,"WMO_Sensor_ID    :",1x,i0)') ODAS%WMO_Sensor_ID
    WRITE(*,'(3x,"Sensor_Type      :",1x,a )') TRIM(SENSOR_TYPE_NAME(ODAS%Sensor_Type))

    ! Arrays summary
    WRITE(*,'(3x,"Sensor_Channel   :")')
    WRITE(*,'(10(1x,i5,:))') ODAS%Sensor_Channel

    WRITE(*,'(3x,"Absorber_ID      :")')
    WRITE(*,'(10(1x,i5,:))') ODAS%Absorber_ID

  END SUBROUTINE ODAS_Inspect

END PROGRAM TauCoeff_Inspect
