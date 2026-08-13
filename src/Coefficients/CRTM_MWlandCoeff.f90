!
! CRTM_MWlandCoeff
!
! Module containing the shared CRTM microwave land surface emissivity atlas
! (TELSEM2) data and its load/destruction routines.
!
! PUBLIC DATA:
!   MWlandC:  TELSEM2Atlas structure containing the microwave land surface
!             emissivity climatology atlas.
!
! SIDE EFFECTS:
!       The load/destroy routines in this module modify the contents of the
!       public data structure MWlandC.
!
! RESTRICTIONS:
!       The load/destroy routines should only be called during CRTM
!       initialisation/destruction. CRTM_TELSEM2_Emissivity_Uncertainty is a
!       read-only query and may be called at any time after CRTM_Init
!       (including from within OpenMP regions).
!

MODULE CRTM_MWlandCoeff

  ! -----------------
  ! Environment setup
  ! -----------------
  USE Type_Kinds            , ONLY: fp
  USE Message_Handler       , ONLY: SUCCESS, FAILURE, INFORMATION, Display_Message
  USE File_Utility          , ONLY: Join_Path
  USE TELSEM2Atlas_Define   , ONLY: TELSEM2Atlas_type     , &
                                    TELSEM2Atlas_Associated, &
                                    TELSEM2Atlas_Destroy
  USE TELSEM2Atlas_netCDF_IO, ONLY: TELSEM2Atlas_netCDF_ReadFile
  USE TELSEM2_Atlas_Module  , ONLY: TELSEM2_Setup_Grid, &
                                    TELSEM2_Emissivity_Cov
  ! Disable all implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  ! The shared data
  PUBLIC :: MWlandC
  ! Procedures
  PUBLIC :: CRTM_MWlandCoeff_Load
  PUBLIC :: CRTM_MWlandCoeff_Destroy
  PUBLIC :: CRTM_MWlandCoeff_IsLoaded
  PUBLIC :: CRTM_TELSEM2_Emissivity_Uncertainty


  ! -----------------
  ! Module parameters
  ! -----------------
  INTEGER, PARAMETER :: ML = 512


  ! ----------------------------------------------------
  ! The shared microwave land surface emissivity atlas
  ! ----------------------------------------------------
  TYPE(TELSEM2Atlas_type), SAVE :: MWlandC


CONTAINS


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_MWlandCoeff_Load
!
! PURPOSE:
!       Function to load the TELSEM2 microwave land surface emissivity atlas
!       into the public data structure MWlandC.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_MWlandCoeff_Load( &
!                        Filename,                              &
!                        File_Path         = File_Path        , &
!                        netCDF            = netCDF           , &
!                        Quiet             = Quiet            , &
!                        Process_ID        = Process_ID       , &
!                        Output_Process_ID = Output_Process_ID  )
!
! INPUT ARGUMENTS:
!       Filename:           Name of the TELSEM2 atlas coefficient file.
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       File_Path:          File path for the input data file. Default is the
!                           current directory.
!       netCDF:             Present for interface consistency with the other
!                           coefficient loaders. The TELSEM2 atlas is only
!                           distributed in netCDF, so this is effectively always
!                           treated as netCDF.
!       Quiet:              Suppress INFORMATION messages if .TRUE.
!       Process_ID:         MPI process ID (message control only).
!       Output_Process_ID:  MPI process ID that emits messages.
!
! FUNCTION RESULT:
!       Error_Status:       SUCCESS or FAILURE.
!
!:sdoc-:
!------------------------------------------------------------------------------
  FUNCTION CRTM_MWlandCoeff_Load( &
    Filename         , &  ! Input
    File_Path        , &  ! Optional input
    netCDF           , &  ! Optional input
    Quiet            , &  ! Optional input
    Process_ID       , &  ! Optional input
    Output_Process_ID) &  ! Optional input
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),           INTENT(IN) :: Filename
    CHARACTER(*), OPTIONAL, INTENT(IN) :: File_Path
    LOGICAL,      OPTIONAL, INTENT(IN) :: netCDF
    LOGICAL,      OPTIONAL, INTENT(IN) :: Quiet
    INTEGER,      OPTIONAL, INTENT(IN) :: Process_ID
    INTEGER,      OPTIONAL, INTENT(IN) :: Output_Process_ID
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_MWlandCoeff_Load'
    ! Local variables
    CHARACTER(ML) :: msg, pid_msg
    CHARACTER(:), ALLOCATABLE :: MWlandCoeff_File
    LOGICAL :: noisy

    ! Setup
    err_stat = SUCCESS
    ! ...Assign the filename to local variable
    MWlandCoeff_File = TRIM(ADJUSTL(Filename))
    ! ...Add the file path
    IF ( PRESENT(File_Path) ) MWlandCoeff_File = Join_Path(File_Path, MWlandCoeff_File)
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

    ! Read the TELSEM2 atlas file (netCDF)
    err_stat = TELSEM2Atlas_netCDF_ReadFile( TRIM(MWlandCoeff_File), MWlandC )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error reading MWlandCoeff TELSEM2 atlas file '//TRIM(MWlandCoeff_File)//TRIM(pid_msg)
      CALL Load_Cleanup(); RETURN
    END IF

    ! Build the equal-area grid geometry and reverse lookup
    CALL TELSEM2_Setup_Grid( MWlandC )

    IF ( noisy ) THEN
      WRITE( msg,'("TELSEM2 MW land emissivity atlas loaded: ",i0," cells, ",i0," months")' ) &
             MWlandC%n_Data, MWlandC%n_Months
      CALL Display_Message( ROUTINE_NAME, TRIM(msg)//TRIM(pid_msg), INFORMATION )
    END IF

  CONTAINS

    SUBROUTINE Load_CleanUp()
      CALL TELSEM2Atlas_Destroy( MWlandC )
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE Load_CleanUp

  END FUNCTION CRTM_MWlandCoeff_Load


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_MWlandCoeff_Destroy
!
! PURPOSE:
!       Function to deallocate the public data structure MWlandC.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_MWlandCoeff_Destroy( Process_ID = Process_ID )
!
!:sdoc-:
!------------------------------------------------------------------------------
  FUNCTION CRTM_MWlandCoeff_Destroy( Process_ID ) RESULT( err_stat )
    ! Arguments
    INTEGER, OPTIONAL, INTENT(IN) :: Process_ID
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_MWlandCoeff_Destroy'
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
    CALL TELSEM2Atlas_Destroy( MWlandC )
    IF ( TELSEM2Atlas_Associated( MWlandC ) ) THEN
      err_stat = FAILURE
      msg = 'Error deallocating MWlandCoeff shared data structure'//TRIM(pid_msg)
      CALL Display_Message( ROUTINE_NAME, msg, err_stat ); RETURN
    END IF

  END FUNCTION CRTM_MWlandCoeff_Destroy


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_MWlandCoeff_IsLoaded
!
! PURPOSE:
!       Function to test if the TELSEM2 MW land emissivity atlas has been
!       loaded into the public data structure MWlandC.
!
!:sdoc-:
!------------------------------------------------------------------------------
  FUNCTION CRTM_MWlandCoeff_IsLoaded() RESULT( IsLoaded )
    LOGICAL :: IsLoaded
    IsLoaded = TELSEM2Atlas_Associated( MWlandC )
  END FUNCTION CRTM_MWlandCoeff_IsLoaded


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_TELSEM2_Emissivity_Uncertainty
!
! PURPOSE:
!       Query the loaded TELSEM2 atlas for the V/H emissivity background and
!       the full inter-channel emissivity error covariance at a location and
!       month, for an arbitrary set of channel frequencies. This is the
!       data-assimilation access point for the TELSEM2 uncertainty content
!       (emissivity background error for emissivity-sink assimilation,
!       observation-error inflation, correlated-R modelling, spectral
!       transfer of retrieved emissivities).
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_TELSEM2_Emissivity_Uncertainty( &
!                        Latitude, Longitude, Month, Frequency, &
!                        Emissivity_V, Emissivity_H, Covariance, Valid, &
!                        Zenith_Angle = Zenith_Angle, &
!                        Resolution   = Resolution )
!
! INPUTS:
!       Latitude:     degrees, -90..90
!       Longitude:    degrees (any range; reduced modulo 360)
!       Month:        1..12
!       Frequency:    channel frequencies, GHz. Rank-1, size n_chan.
!
! OUTPUTS:
!       Emissivity_V: V-pol emissivity background (n_chan)
!       Emissivity_H: H-pol emissivity background (n_chan)
!       Covariance:   emissivity error covariance (2*n_chan x 2*n_chan);
!                     rows/columns 1..n_chan are V-pol, n_chan+1..2*n_chan
!                     are H-pol. Frequency-interpolated from the atlas
!                     channels; no angular dependence.
!       Valid:        .FALSE. (with zeroed outputs) when the atlas has no
!                     climatology at the requested cell.
!
! OPTIONAL INPUTS:
!       Zenith_Angle: degrees; affects the emissivity outputs only
!                     (default 0 = nadir).
!       Resolution:   spatial-averaging resolution in degrees (>= the native
!                     0.25); emissivities and covariance are averaged over
!                     the contributing cells.
!
! FUNCTION RESULT:
!       Error_Status: SUCCESS, or FAILURE when the atlas is not loaded, the
!                     loaded coefficient file carries no uncertainty content
!                     (Release 1), or the output array dimensions do not
!                     match SIZE(Frequency).
!
! RESTRICTIONS:
!       Read-only; callable any time after CRTM_Init with the MW land atlas
!       opted in, including from within OpenMP regions.
!
!:sdoc-:
!------------------------------------------------------------------------------
  FUNCTION CRTM_TELSEM2_Emissivity_Uncertainty( &
    Latitude    , &  ! Input
    Longitude   , &  ! Input
    Month       , &  ! Input
    Frequency   , &  ! Input
    Emissivity_V, &  ! Output
    Emissivity_H, &  ! Output
    Covariance  , &  ! Output
    Valid       , &  ! Output
    Zenith_Angle, &  ! Optional input
    Resolution  ) &  ! Optional input
  RESULT( err_stat )
    ! Arguments
    REAL(fp),           INTENT(IN)  :: Latitude
    REAL(fp),           INTENT(IN)  :: Longitude
    INTEGER,            INTENT(IN)  :: Month
    REAL(fp),           INTENT(IN)  :: Frequency(:)
    REAL(fp),           INTENT(OUT) :: Emissivity_V(:)
    REAL(fp),           INTENT(OUT) :: Emissivity_H(:)
    REAL(fp),           INTENT(OUT) :: Covariance(:,:)
    LOGICAL,            INTENT(OUT) :: Valid
    REAL(fp), OPTIONAL, INTENT(IN)  :: Zenith_Angle
    REAL(fp), OPTIONAL, INTENT(IN)  :: Resolution
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_TELSEM2_Emissivity_Uncertainty'
    ! Local variables
    CHARACTER(ML) :: msg
    INTEGER  :: n_chan
    REAL(fp) :: theta

    err_stat = SUCCESS
    Emissivity_V(:)  = 0.0_fp
    Emissivity_H(:)  = 0.0_fp
    Covariance(:,:)  = 0.0_fp
    Valid            = .FALSE.

    ! Checks
    IF ( .NOT. CRTM_MWlandCoeff_IsLoaded() ) THEN
      err_stat = FAILURE
      msg = 'TELSEM2 MW land emissivity atlas is not loaded'
      CALL Display_Message( ROUTINE_NAME, msg, err_stat ); RETURN
    END IF
    IF ( .NOT. MWlandC%Has_Uncertainty ) THEN
      err_stat = FAILURE
      msg = 'Loaded TELSEM2 atlas carries no uncertainty content (Release 1 file)'
      CALL Display_Message( ROUTINE_NAME, msg, err_stat ); RETURN
    END IF
    n_chan = SIZE(Frequency)
    IF ( SIZE(Emissivity_V) /= n_chan .OR. SIZE(Emissivity_H) /= n_chan .OR. &
         SIZE(Covariance,1) /= 2*n_chan .OR. SIZE(Covariance,2) /= 2*n_chan ) THEN
      err_stat = FAILURE
      msg = 'Output array dimensions inconsistent with SIZE(Frequency)'
      CALL Display_Message( ROUTINE_NAME, msg, err_stat ); RETURN
    END IF
    IF ( Month < 1 .OR. Month > 12 ) THEN
      err_stat = FAILURE
      msg = 'Month must be in 1..12'
      CALL Display_Message( ROUTINE_NAME, msg, err_stat ); RETURN
    END IF

    theta = 0.0_fp
    IF ( PRESENT(Zenith_Angle) ) theta = Zenith_Angle

    IF ( PRESENT(Resolution) ) THEN
      CALL TELSEM2_Emissivity_Cov( MWlandC, Latitude, Longitude, Month, theta, &
                                   Frequency, Emissivity_V, Emissivity_H, Valid, &
                                   Covariance=Covariance, Resolution=Resolution )
    ELSE
      CALL TELSEM2_Emissivity_Cov( MWlandC, Latitude, Longitude, Month, theta, &
                                   Frequency, Emissivity_V, Emissivity_H, Valid, &
                                   Covariance=Covariance )
    END IF
  END FUNCTION CRTM_TELSEM2_Emissivity_Uncertainty

END MODULE CRTM_MWlandCoeff
