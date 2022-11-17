!
! MWwaterLUT_Define
!
! Module defining the MWwaterLUT object containing the
! Look-Up Table (LUT) for the microWave (MW) sea surface emissivity
! model.
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 10-Nov-2011
!                       paul.vandelst@noaa.gov
!

MODULE CRTM_MWwaterLUT_Define

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds,        ONLY: Long   => CSEM_Long, fp => CSEM_fp, &
                                    Double => CSEM_Double
  USE CSEM_Exception_Handler, ONLY: SUCCESS, FAILURE, INFORMATION, Display_Message
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Datatypes
  PUBLIC :: MWwaterLUT_type

  ! Procedures
  PUBLIC :: MWwaterLUT_Associated
  PUBLIC :: MWwaterLUT_Destroy
  PUBLIC :: MWwaterLUT_Create
  PUBLIC :: MWwaterLUT_Inspect
  PUBLIC :: MWwaterLUT_ValidRelease
  PUBLIC :: MWwaterLUT_Info
  PUBLIC :: MWwaterLUT_DefineVersion


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
    '$Id: MWwaterLUT_Define.f90 21141 2012-09-14 17:40:43Z paul.vandelst@noaa.gov $'
  ! Release and version
  INTEGER, PARAMETER :: MWWATERLUT_RELEASE = 1  ! This determines structure and file formats.
  INTEGER, PARAMETER :: MWWATERLUT_VERSION = 1  ! This is just the default data version.
  ! Close status for write errors
  CHARACTER(*), PARAMETER :: WRITE_ERROR_STATUS = 'DELETE'
  ! Literal constants
  REAL(fp), PARAMETER :: ZERO = 0.0_fp
  REAL(fp), PARAMETER :: ONE  = 1.0_fp
  ! String lengths
  INTEGER,  PARAMETER :: ML = 256 ! Message length
  INTEGER,  PARAMETER :: SL =  80 ! String length


  ! ----------------------------------
  ! MWwaterLUT data type definitions
  ! ----------------------------------
  !:tdoc+:
  TYPE :: MWwaterLUT_type
    ! Allocation indicator
    LOGICAL :: Is_Allocated = .FALSE.
    ! Release and version information
    INTEGER(Long) :: Release = MWwaterLUT_RELEASE
    INTEGER(Long) :: Version = MWwaterLUT_VERSION
    ! Dimensions
    INTEGER(Long) :: n_Angles       = 0   ! I1 dimension
    INTEGER(Long) :: n_Frequencies  = 0   ! I2 dimension
    INTEGER(Long) :: n_Temperatures = 0   ! I3 dimension
    INTEGER(Long) :: n_Wind_Speeds  = 0   ! I4 dimension
    ! Dimensional vectors
    REAL(Double),  ALLOCATABLE :: Angle(:)        ! I1
    REAL(Double),  ALLOCATABLE :: Frequency(:)    ! I2
    REAL(Double),  ALLOCATABLE :: Temperature(:)  ! I3
    REAL(Double),  ALLOCATABLE :: Wind_Speed(:)   ! I4
    ! Large-scale correction emissivity data
    REAL(Double),  ALLOCATABLE :: ev(:,:,:,:)  ! I1 x I2 x I3 x I4
    REAL(Double),  ALLOCATABLE :: eh(:,:,:,:)  ! I1 x I2 x I3 x I4
  END TYPE MWwaterLUT_type
  !:tdoc-:


CONTAINS


!################################################################################
!################################################################################
!##                                                                            ##
!##                           ## PUBLIC PROCEDURES ##                          ##
!##                                                                            ##
!################################################################################
!################################################################################

!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       MWwaterLUT_Associated
!
! PURPOSE:
!       Pure function to test the status of the allocatable components
!       of the MWwaterLUT structure.
!
! CALLING SEQUENCE:
!       Status = MWwaterLUT_Associated( MWwaterLUT )
!
! OBJECTS:
!       MWwaterLUT:    Structure which is to have its member's
!                      status tested.
!                      UNITS:      N/A
!                      TYPE:       MWwaterLUT_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN)
!
! FUNCTION RESULT:
!       Status:        The return value is a logical value indicating the
!                      status of the components.
!                       .TRUE.  - if ANY of the MWwaterLUT allocatable members
!                                 are in use.
!                       .FALSE. - if ALL of the MWwaterLUT allocatable members
!                                 are not in use.
!                      UNITS:      N/A
!                      TYPE:       LOGICAL
!                      DIMENSION:  Same as input
!
!:sdoc-:
!--------------------------------------------------------------------------------

  PURE FUNCTION MWwaterLUT_Associated( self ) RESULT( Status )
    TYPE(MWwaterLUT_type), INTENT(IN) :: self
    LOGICAL :: Status
    Status = self%Is_Allocated
  END FUNCTION MWwaterLUT_Associated


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       MWwaterLUT_Destroy
!
! PURPOSE:
!       Pure subroutine to re-initialize MWwaterLUT objects.
!
! CALLING SEQUENCE:
!       CALL MWwaterLUT_Destroy( MWwaterLUT )
!
! OBJECTS:
!       MWwaterLUT:    Re-initialized MWwaterLUT structure.
!                      UNITS:      N/A
!                      TYPE:       MWwaterLUT_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  PURE SUBROUTINE MWwaterLUT_Destroy( self )
    TYPE(MWwaterLUT_type), INTENT(OUT) :: self
    self%Is_Allocated = .FALSE.
    self%n_Angles       = 0
    self%n_Frequencies  = 0
    self%n_Temperatures = 0
    self%n_Wind_Speeds  = 0
  END SUBROUTINE MWwaterLUT_Destroy


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       MWwaterLUT_Create
!
! PURPOSE:
!       Pure subroutine to create an instance of an MWwaterLUT object.
!
! CALLING SEQUENCE:
!       CALL MWwaterLUT_Create( MWwaterLUT    , &
!                               n_Angles      , &
!                               n_Frequencies , &
!                               n_Temperatures, &
!                               n_Wind_Speeds   )
!
! OBJECTS:
!       MWwaterLUT:         MWwaterLUT object structure.
!                           UNITS:      N/A
!                           TYPE:       MWwaterLUT_type
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(OUT)
!
! INPUTS:
!       n_Angles:           Number of zenith angles for which is are data.
!                           Must be > 0.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
!
!       n_Frequencies:      Number of spectral frequencies for which there are
!                           data.
!                           Must be > 0.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
!
!       n_Temperatures:     Number of surface temperatures for which there are
!                           data.
!                           Must be > 0.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
!
!       n_Wind_Speeds:      Number of surface wind speeds for which there are
!                           data.
!                           Must be > 0.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  PURE SUBROUTINE MWwaterLUT_Create( &
    self          , &  ! Output
    n_Angles      , &  ! Input
    n_Frequencies , &  ! Input
    n_Temperatures, &  ! Input
    n_Wind_Speeds   )  ! Input
    ! Arguments
    TYPE(MWwaterLUT_type), INTENT(OUT) :: self
    INTEGER              , INTENT(IN)  :: n_Angles
    INTEGER              , INTENT(IN)  :: n_Frequencies
    INTEGER              , INTENT(IN)  :: n_Temperatures
    INTEGER              , INTENT(IN)  :: n_Wind_Speeds
    ! Local variables
    INTEGER :: alloc_stat

    ! Check input
    IF ( n_Angles       < 1 .OR. &
         n_Frequencies  < 1 .OR. &
         n_Temperatures < 1 .OR. &
         n_Wind_Speeds  < 1 ) RETURN


    ! Perform the allocation
    ALLOCATE( self%Angle( n_Angles ), &
              self%Frequency( n_Frequencies ), &
              self%Temperature( n_Temperatures ), &
              self%Wind_Speed( n_Wind_Speeds ), &
              self%ev( n_Angles, n_Frequencies, n_Temperatures, n_Wind_Speeds ), &
              self%eh( n_Angles, n_Frequencies, n_Temperatures, n_Wind_Speeds ), &
              STAT = alloc_stat )
    IF ( alloc_stat /= 0 ) RETURN


    ! Initialise
    ! ...Dimensions
    self%n_Angles       = n_Angles
    self%n_Frequencies  = n_Frequencies
    self%n_Temperatures = n_Temperatures
    self%n_Wind_Speeds  = n_Wind_Speeds
    ! ...Arrays
    self%Angle       = ZERO
    self%Frequency   = ZERO
    self%Temperature = ZERO
    self%Wind_Speed  = ZERO
    self%ev          = ZERO
    self%eh          = ZERO

    ! Set allocation indicator
    self%Is_Allocated = .TRUE.

  END SUBROUTINE MWwaterLUT_Create


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       MWwaterLUT_Inspect
!
! PURPOSE:
!       Subroutine to print the contents of a MWwaterLUT object to stdout.
!
! CALLING SEQUENCE:
!       CALL MWwaterLUT_Inspect( MWwaterLUT )
!
! OBJECTS:
!       MWwaterLUT:     MWwaterLUT object to display.
!                       UNITS:      N/A
!                       TYPE:       MWwaterLUT_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE MWwaterLUT_Inspect( self, pause )
    TYPE(MWwaterLUT_type), INTENT(IN) :: self
    LOGICAL,     OPTIONAL, INTENT(IN) :: pause
    LOGICAL :: wait
    INTEGER :: i2, i3, i4

    wait = .FALSE.
    IF ( PRESENT(pause) ) wait = pause

    WRITE(*,'(1x,"MWwaterLUT OBJECT")')
    ! Release/version info
    WRITE(*,'(3x,"Release.Version : ",i0,".",i0)') self%Release, self%Version
    ! Dimensions
    WRITE(*,'(3x,"n_Angles        : ",i0)') self%n_Angles
    WRITE(*,'(3x,"n_Frequencies   : ",i0)') self%n_Frequencies
    WRITE(*,'(3x,"n_Temperatures  : ",i0)') self%n_Temperatures
    WRITE(*,'(3x,"n_Wind_Speeds   : ",i0)') self%n_Wind_Speeds
    IF ( .NOT. MWwaterLUT_Associated(self) ) RETURN
    ! Dimension arrays
    WRITE(*,'(3x,"Angle :")')
    WRITE(*,'(5(1x,es13.6,:))') self%Angle
    WRITE(*,'(3x,"Frequency :")')
    WRITE(*,'(5(1x,es13.6,:))') self%Frequency
    WRITE(*,'(3x,"Temperature :")')
    WRITE(*,'(5(1x,es13.6,:))') self%Temperature
    WRITE(*,'(3x,"Wind_Speed :")')
    WRITE(*,'(5(1x,es13.6,:))') self%Wind_Speed

    ! Emissivity arrays
    WRITE(*,'(/3x,"Emissivity(vertical polarisation) :")')
    IF ( wait ) THEN
      WRITE(*,FMT='(/1x,"Paused. Press <ENTER> to continue...")',ADVANCE='NO')
      READ(*,*)
    END IF

    DO i4 = 1, self%n_Wind_Speeds
      WRITE(*,'(5x,"WIND_SPEED  :",es13.6)') self%Wind_Speed(i4)
      DO i3 = 1, self%n_Temperatures
        WRITE(*,'(5x,"TEMPERATURE :",es13.6)') self%Temperature(i3)
        DO i2 = 1, self%n_Frequencies
          WRITE(*,'(5x,"FREQUENCY   :",es13.6)') self%Frequency(i2)
          WRITE(*,'(5(1x,es13.6,:))') self%ev(:,i2,i3,i4)
        END DO
      END DO
    END DO

    WRITE(*,'(/3x,"Emissivity(horizontal polarisation) :")')
    IF ( wait ) THEN
      WRITE(*,FMT='(/1x,"Paused. Press <ENTER> to continue...")',ADVANCE='NO')
      READ(*,*)
    END IF

    DO i4 = 1, self%n_Wind_Speeds
      WRITE(*,'(5x,"WIND_SPEED  :",es13.6)') self%Wind_Speed(i4)
      DO i3 = 1, self%n_Temperatures
        WRITE(*,'(5x,"TEMPERATURE :",es13.6)') self%Temperature(i3)
        DO i2 = 1, self%n_Frequencies
          WRITE(*,'(5x,"FREQUENCY   :",es13.6)') self%Frequency(i2)
          WRITE(*,'(5(1x,es13.6,:))') self%eh(:,i2,i3,i4)
        END DO
      END DO
    END DO
  END SUBROUTINE MWwaterLUT_Inspect



!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       MWwaterLUT_ValidRelease
!
! PURPOSE:
!       Function to check the MWwaterLUT Release value.
!
! CALLING SEQUENCE:
!       IsValid = MWwaterLUT_ValidRelease( MWwaterLUT )
!
! INPUTS:
!       MWwaterLUT:     MWwaterLUT object for which the Release component
!                       is to be checked.
!                       UNITS:      N/A
!                       TYPE:       MWwaterLUT_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
! FUNCTION RESULT:
!       IsValid:        Logical value defining the release validity.
!                       UNITS:      N/A
!                       TYPE:       LOGICAL
!                       DIMENSION:  Scalar
!
!:sdoc-:
!----------------------------------------------------------------------------------

  FUNCTION MWwaterLUT_ValidRelease( self ) RESULT( IsValid )
    ! Arguments
    TYPE(MWwaterLUT_type), INTENT(IN) :: self
    ! Function result
    LOGICAL :: IsValid
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'MWwaterLUT_ValidRelease'
    ! Local variables
    CHARACTER(ML) :: msg

    ! Set up
    IsValid = .TRUE.


    ! Check release is not too old
    IF ( self%Release < MWWATERLUT_RELEASE ) THEN
      IsValid = .FALSE.
      WRITE( msg,'("An MWwaterLUT data update is needed. ", &
                  &"MWwaterLUT release is ",i0,". Valid release is ",i0,"." )' ) &
                  self%Release, MWWATERLUT_RELEASE
      CALL Display_Message( ROUTINE_NAME, msg, INFORMATION ); RETURN
    END IF


    ! Check release is not too new
    IF ( self%Release > MWWATERLUT_RELEASE ) THEN
      IsValid = .FALSE.
      WRITE( msg,'("An MWwaterLUT software update is needed. ", &
                  &"MWwaterLUT release is ",i0,". Valid release is ",i0,"." )' ) &
                  self%Release, MWWATERLUT_RELEASE
      CALL Display_Message( ROUTINE_NAME, msg, INFORMATION ); RETURN
    END IF

  END FUNCTION MWwaterLUT_ValidRelease


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       MWwaterLUT_Info
!
! PURPOSE:
!       Subroutine to return a string containing version and dimension
!       information about a MWwaterLUT object.
!
! CALLING SEQUENCE:
!       CALL MWwaterLUT_Info( MWwaterLUT, Info )
!
! OBJECTS:
!       MWwaterLUT:      MWwaterLUT object about which info is required.
!                        UNITS:      N/A
!                        TYPE:       MWwaterLUT_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       Info:            String containing version and dimension information
!                        about the MWwaterLUT object.
!                        UNITS:      N/A
!                        TYPE:       CHARACTER(*)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE MWwaterLUT_Info( self, Info )
    ! Arguments
    TYPE(MWwaterLUT_type), INTENT(IN)  :: self
    CHARACTER(*)         , INTENT(OUT) :: Info
    ! Parameters
    INTEGER, PARAMETER :: CARRIAGE_RETURN = 13
    INTEGER, PARAMETER :: LINEFEED = 10
    ! Local variables
    CHARACTER(2000) :: Long_String

    ! Write the required data to the local string
    WRITE( Long_String, &
           '(a,1x,"MWwaterLUT RELEASE.VERSION: ",i2,".",i2.2,a,3x, &
           &"N_ANGLES=",i0,2x,&
           &"N_FREQUENCIES=",i0,2x,&
           &"N_TEMPERATURES=",i0,2x,&
           &"N_WIND_SPEEDS=",i0 )' ) &
           ACHAR(CARRIAGE_RETURN)//ACHAR(LINEFEED), &
           self%Release, self%Version, &
           ACHAR(CARRIAGE_RETURN)//ACHAR(LINEFEED), &
           self%n_Angles      , &
           self%n_Frequencies , &
           self%n_Temperatures, &
           self%n_Wind_Speeds

    ! Trim the output based on the
    ! dummy argument string length
    Info = Long_String(1:MIN(LEN(Info), LEN_TRIM(Long_String)))

  END SUBROUTINE MWwaterLUT_Info


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       MWwaterLUT_DefineVersion
!
! PURPOSE:
!       Subroutine to return the module version information.
!
! CALLING SEQUENCE:
!       CALL MWwaterLUT_DefineVersion( Id )
!
! OUTPUTS:
!       Id:    Character string containing the version Id information
!              for the module.
!              UNITS:      N/A
!              TYPE:       CHARACTER(*)
!              DIMENSION:  Scalar
!              ATTRIBUTES: INTENT(OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE MWwaterLUT_DefineVersion( Id )
    CHARACTER(*), INTENT(OUT) :: Id
    Id = MODULE_VERSION_ID
  END SUBROUTINE MWwaterLUT_DefineVersion



END MODULE CRTM_MWwaterLUT_Define
