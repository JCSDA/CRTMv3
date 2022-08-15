!
! CRTM_MoleculeScatter
!
! Module to compute molecule optical properties.
!       
!
! CREATION HISTORY:
!       Written by:     Quanhua Liu, 03-Oct-2008
!                       Quanhua.Liu@noaa.gov
!

MODULE CRTM_MoleculeScatter

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE Type_Kinds            , ONLY: fp
  USE Message_Handler       , ONLY: SUCCESS, FAILURE, Display_Message
  USE CRTM_Parameters       , ONLY: ZERO, ONE, TWO
  USE CRTM_Atmosphere_Define, ONLY: CRTM_Atmosphere_type
  USE CRTM_AtmOptics_Define , ONLY: CRTM_AtmOptics_type
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  PUBLIC :: CRTM_Compute_MoleculeScatter
  PUBLIC :: CRTM_Compute_MoleculeScatter_TL
  PUBLIC :: CRTM_Compute_MoleculeScatter_AD


  ! -----------------
  ! Module parameters
  ! -----------------
  ! Rayleigh factor
  REAL(fp),     PARAMETER :: RFACTOR = 27.0363_fp  ! = 287.0/9.8*923.1907/1000.0


CONTAINS


!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Compute_MoleculeScatter
!
! PURPOSE:
!       Function to compute molecular scattering and extinction
!
! CALLING SEQUENCE:
!       Error_Status =  CRTM_Compute_MoleculeScatter( Wavenumber             , &  ! Input
!                                                     Atmosphere             , &  ! Input
!                                                     AtmOptics              , &  ! In/Output
!                                                     Message_Log=Message_Log  )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Wavenumber:      Spectral frequency
!                        UNITS:      Inverse centimetres (cm^-1)
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Atmosphere:      Structure containing the atmospheric
!                        profile data.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_Atmosphere_type)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       AtmOptics:       Structure containing the atmospheric optics data
!                        to which the molecular scattering and extinction
!                        component is added.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_AtmOptics_type)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
! OPTIONAL INPUT ARGUMENTS:
!       Message_Log:     Character string specifying a filename in which any
!                        messages will be logged. If not specified, or if an
!                        error occurs opening the log file, the default action
!                        is to output messages to standard output.
!                        UNITS:      N/A
!                        TYPE:       CHARACTER(*)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status:    The return value is an integer defining the error status.
!                        The error codes are defined in the ERROR_HANDLER module.
!                        If == SUCCESS the computation was sucessful
!                           == FAILURE an unrecoverable error occurred
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_Compute_MoleculeScatter( Wavenumber,     &  ! Input
                                         Atmosphere,     &  ! Input
                                         AtmOptics,      &  ! In/Output
                                         Message_Log )   &  ! Error messaging
                                      RESULT( Error_Status )
    ! Arguments
    REAL(fp),                   INTENT(IN)     :: Wavenumber
    TYPE(CRTM_Atmosphere_type), INTENT(IN)     :: Atmosphere
    TYPE(CRTM_AtmOptics_type), INTENT(IN OUT) :: AtmOptics 
    CHARACTER(*), OPTIONAL,     INTENT(IN)     :: Message_Log
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Compute_MoleculeScatter'

    ! Local variables
    INTEGER :: k
    REAL(fp) :: Wavelength , Opt_unit, Optical_Depth, dpF, dpF2

    ! Setup
    ! -----
    Error_Status = SUCCESS
    
    ! Check input
    IF( Wavenumber > ZERO ) THEN
      Wavelength = Compute_Wavelength( Wavenumber )
    ELSE
      Error_Status = FAILURE
      CALL Display_Message(ROUTINE_NAME,'Invalid wavenumber',Error_Status,Message_Log=Message_Log)
      RETURN
    END IF


    ! Calculate the scattering parameters
    ! -----------------------------------
    ! Compute optical scaling unit                                                                   
    CALL RAYLO(Wavelength, Opt_unit, AtmOptics%depolarization)
!  Opt_unit = Opt_unit*(6.0_fp + AtmOptics%depolarization*3.0_fp)/(6.0_fp-7.0_fp*AtmOptics%depolarization)
    dpF = (ONE - AtmOptics%depolarization)/(ONE+AtmOptics%depolarization/TWO)
    dpF2 = (ONE-TWO*AtmOptics%depolarization)/(ONE+AtmOptics%depolarization/TWO)
    
    ! Loop over atmospheric layers
    DO k = 1, Atmosphere%n_Layers
      Optical_Depth = RFACTOR*Opt_unit*(Atmosphere%Level_Pressure(k)-Atmosphere%Level_Pressure(k-1))
      AtmOptics%Optical_Depth(k)         = AtmOptics%Optical_Depth(k)         + Optical_Depth 
      AtmOptics%Single_Scatter_Albedo(k) = AtmOptics%Single_Scatter_Albedo(k) + Optical_Depth 

      ! The Rayleigh spherical expansion coefficients are constant.
      AtmOptics%Phase_Coefficient(2,1,k) = AtmOptics%Phase_Coefficient(2,1,k) + &
                                       dpF * 0.25_fp * Optical_Depth 
      IF( AtmOptics%n_Stokes > 1 ) THEN
         AtmOptics%Phase_Coefficient(2,2,k) = AtmOptics%Phase_Coefficient(2,2,k) + &
                              dpF * 1.5_fp * Optical_Depth
      END IF
      IF( AtmOptics%n_Stokes > 2 ) THEN                                             
        AtmOptics%Phase_Coefficient(1,4,k) = AtmOptics%Phase_Coefficient(1,4,k) + &
                                    dpF2 * 0.75_fp * Optical_Depth 
        AtmOptics%Phase_Coefficient(2,5,k) = AtmOptics%Phase_Coefficient(2,5,k) - &
                                      dpF * 0.612372_fp * Optical_Depth
      END IF 

    END DO
  END FUNCTION CRTM_Compute_MoleculeScatter


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Compute_MoleculeScatter_TL
!
! PURPOSE:
!       Function to compute the tangent-linear molecular scattering and
!       extinction
!
! CALLING SEQUENCE:
!       Error_Status =  CRTM_Compute_MoleculeScatter_TL( Wavenumber             , &  ! Input
!                                                        Atmosphere_TL          , &  ! Input
!                                                        AtmOptics_TL           , &  ! In/Output
!                                                        Message_Log=Message_Log  )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Wavenumber:      Spectral frequency
!                        UNITS:      Inverse centimetres (cm^-1)
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Atmosphere_TL:   Structure containing the tangent-linear atmospheric
!                        profile data.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_Atmosphere_type)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       AtmOptics_TL:    Structure containing the tangent-linear atmospheric
!                        optics data to which the molecular scattering and
!                        extinction component is added.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_AtmOptics_type)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
! OPTIONAL INPUT ARGUMENTS:
!       Message_Log:     Character string specifying a filename in which any
!                        messages will be logged. If not specified, or if an
!                        error occurs opening the log file, the default action
!                        is to output messages to standard output.
!                        UNITS:      N/A
!                        TYPE:       CHARACTER(*)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status:    The return value is an integer defining the error status.
!                        The error codes are defined in the ERROR_HANDLER module.
!                        If == SUCCESS the computation was sucessful
!                           == FAILURE an unrecoverable error occurred
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_Compute_MoleculeScatter_TL( Wavenumber,     &  ! Input
                                            Atmosphere_TL,  &  ! TL Input
                                            AtmOptics_TL,   &  ! TL Output
                                            Message_Log )   &  ! Error messaging
                                         RESULT( Error_Status )
    ! Arguments
    REAL(fp),                   INTENT(IN)     :: Wavenumber
    TYPE(CRTM_Atmosphere_type), INTENT(IN)     :: Atmosphere_TL
    TYPE(CRTM_AtmOptics_type), INTENT(IN OUT) :: AtmOptics_TL
    CHARACTER(*), OPTIONAL,     INTENT(IN)     :: Message_Log
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Compute_MoleculeScatter_TL'
    ! Local variables
    INTEGER :: k
    REAL(fp) :: Wavelength , Opt_unit, Optical_Depth_TL, dpF, dpF2

    ! Setup
    ! -----
    Error_Status = SUCCESS
    
    ! Check input
    IF( Wavenumber > ZERO ) THEN
      Wavelength = Compute_Wavelength( Wavenumber )
    ELSE
      Error_Status = FAILURE
      CALL Display_Message(ROUTINE_NAME,'Invalid wavenumber',Error_Status,Message_Log=Message_Log)
      RETURN
    END IF

    dpF = (ONE - AtmOptics_TL%depolarization)/(ONE+AtmOptics_TL%depolarization/TWO)
    dpF2 = (ONE-TWO*AtmOptics_TL%depolarization)/(ONE+AtmOptics_TL%depolarization/TWO)
    ! Calculate the TL scattering parameters
    ! --------------------------------------
    CALL RAYLO(Wavelength, Opt_unit, AtmOptics_TL%depolarization)
    DO k = 1, Atmosphere_TL%n_Layers
      Optical_Depth_TL = RFACTOR*Opt_unit*(Atmosphere_TL%Level_Pressure(k)-Atmosphere_TL%Level_Pressure(k-1))
      AtmOptics_TL%Optical_Depth(k)         = AtmOptics_TL%Optical_Depth(k)         + Optical_Depth_TL
      AtmOptics_TL%Single_Scatter_Albedo(k) = AtmOptics_TL%Single_Scatter_Albedo(k) + Optical_Depth_TL
      ! The Rayleigh spherical expansion coefficients are constant.
      AtmOptics_TL%Phase_Coefficient(2,1,k) = AtmOptics_TL%Phase_Coefficient(2,1,k) + &
                                             dpF * 0.25_fp * Optical_Depth_TL

    IF( AtmOptics_TL%n_Stokes > 1 ) THEN
       AtmOptics_TL%Phase_Coefficient(2,2,k) = AtmOptics_TL%Phase_Coefficient(2,2,k) + &
                                         dpF * 1.5_fp * Optical_Depth_TL
    END IF
    IF( AtmOptics_TL%n_Stokes > 2 ) THEN                                             
      AtmOptics_TL%Phase_Coefficient(1,4,k) = AtmOptics_TL%Phase_Coefficient(1,4,k) + &
                                        dpF2 * 0.75_fp * Optical_Depth_TL
      AtmOptics_TL%Phase_Coefficient(2,5,k) = AtmOptics_TL%Phase_Coefficient(2,5,k) - &
                                        dpF * 0.612372_fp * Optical_Depth_TL
    END IF 
   END DO

  END FUNCTION CRTM_Compute_MoleculeScatter_TL


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Compute_MoleculeScatter_AD
!
! PURPOSE:
!       Function to compute the molecular scattering and extinction adjoint.
!
! CALLING SEQUENCE:
!       Error_Status =  CRTM_Compute_MoleculeScatter_AD( Wavenumber             , &  ! Input
!                                                        AtmOptics_AD           , &  ! In/Output
!                                                        Atmosphere_AD          , &  ! Input
!                                                        Message_Log=Message_Log  )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Wavenumber:      Spectral frequency
!                        UNITS:      Inverse centimetres (cm^-1)
!                        TYPE:       REAL(fp)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       AtmOptics_AD:    Structure containing the adjoint atmospheric
!                        optics data from which the molecular scattering and
!                        extinction component is taken.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_AtmOptics_type)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
! OUTPUT ARGUMENTS:
!       Atmosphere_AD:   Structure containing the adjoint atmospheric
!                        profile data.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_Atmosphere_type)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
! OPTIONAL INPUT ARGUMENTS:
!       Message_Log:     Character string specifying a filename in which any
!                        messages will be logged. If not specified, or if an
!                        error occurs opening the log file, the default action
!                        is to output messages to standard output.
!                        UNITS:      N/A
!                        TYPE:       CHARACTER(*)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status:    The return value is an integer defining the error status.
!                        The error codes are defined in the ERROR_HANDLER module.
!                        If == SUCCESS the computation was sucessful
!                           == FAILURE an unrecoverable error occurred
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_Compute_MoleculeScatter_AD( Wavenumber,     &  ! Input
                                            AtmOptics_AD,   &  ! AD Input
                                            Atmosphere_AD,  &  ! AD Output
                                            Message_Log )   &  ! Error messaging
                                         RESULT( Error_Status )
    ! Arguments
    REAL(fp),                   INTENT(IN)     :: Wavenumber
    TYPE(CRTM_AtmOptics_type), INTENT(IN OUT) :: AtmOptics_AD
    TYPE(CRTM_Atmosphere_type), INTENT(IN OUT) :: Atmosphere_AD
    CHARACTER(*), OPTIONAL,     INTENT(IN)     :: Message_Log
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Compute_MoleculeScatter_AD'
    INTEGER :: k
    REAL(fp) :: Wavelength , Opt_unit, Optical_Depth_AD, dpF, dpF2

    ! Setup
    ! -----
    Error_Status = SUCCESS
    
    ! Check input
    IF( Wavenumber > ZERO ) THEN
      Wavelength = Compute_Wavelength( Wavenumber )
    ELSE
      Error_Status = FAILURE
      CALL Display_Message(ROUTINE_NAME,'Invalid wavenumber',Error_Status,Message_Log=Message_Log)
      RETURN
    END IF

    dpF = (ONE - AtmOptics_AD%depolarization)/(ONE+AtmOptics_AD%depolarization/TWO)
    dpF2 = (ONE-TWO*AtmOptics_AD%depolarization)/(ONE+AtmOptics_AD%depolarization/TWO)
    ! Calculate the AD scattering parameters
    ! --------------------------------------
    CALL RAYLO(Wavelength, Opt_unit, AtmOptics_AD%depolarization)
    DO k = 1, Atmosphere_AD%n_Layers
     ! Compute atmospheric polarisatin component
      IF( AtmOptics_AD%n_Stokes > 2 ) THEN
        Optical_Depth_AD=Optical_Depth_AD+dpF2*0.75_fp*AtmOptics_AD%Phase_Coefficient(1,4,k)                                      
        Optical_Depth_AD=Optical_Depth_AD-dpF*0.612372_fp*AtmOptics_AD%Phase_Coefficient(2,5,k)  
      END IF 
      IF( AtmOptics_AD%n_Stokes > 1 ) THEN
        Optical_Depth_AD=Optical_Depth_AD+dpF*1.5_fp*AtmOptics_AD%Phase_Coefficient(2,2,k)  
      END IF
      Optical_Depth_AD = AtmOptics_AD%Single_Scatter_Albedo(k)
      Optical_Depth_AD = Optical_Depth_AD +dpF*0.25_fp*AtmOptics_AD%Phase_Coefficient(2,1,k)
      Atmosphere_AD%Level_Pressure(k)   = Atmosphere_AD%Level_Pressure(k)   + &
                                            RFACTOR*Opt_unit*Optical_Depth_AD
      Atmosphere_AD%Level_Pressure(k-1) = Atmosphere_AD%Level_Pressure(k-1) - &
                                            RFACTOR*Opt_unit*Optical_Depth_AD
    END DO
    
  END FUNCTION CRTM_Compute_MoleculeScatter_AD


!##################################################################################
!##################################################################################
!##                                                                              ##
!##                          ## PRIVATE MODULE ROUTINES ##                       ##
!##                                                                              ##
!##################################################################################
!##################################################################################

  ! ---------------------------------------------------------------------
  ! Simple function to convert wavenumber (cm^-1) to wavelength (microns)
  ! ---------------------------------------------------------------------
  FUNCTION Compute_Wavelength( Wavenumber ) RESULT( Wavelength )
    REAL(fp), INTENT(IN) :: Wavenumber
    REAL(fp) :: Wavelength
    REAL(fp), PARAMETER :: WFACTOR = 10000.0_fp
    Wavelength = WFACTOR/Wavenumber
  END FUNCTION Compute_Wavelength

  
  !--------------------------------------------------------------------
  !   OPTICAL DEPTH FOR Molecule SCATTERING
  ! ARGUMENTS -
  !     WE      R*8    IN      WAVELENGTH ( MICRO METER)
  !     RAYLO   R*8    OUT     OPTICAL DEPTH PER KM
  !--------------------------------------------------------------------
  SUBROUTINE RAYLO(WE, OPT, DELT)
     REAL(fp), INTENT(IN)  :: WE
     REAL(fp), INTENT(OUT) :: OPT
     REAL(fp) :: DELT
     REAL(fp) :: X1, DY, X2, AS
     X1=1.0_fp/(WE*WE)
     AS=(6432.8_fp+2949810.0_fp/(146.0_fp-X1)+25540.0_fp/(41.0_fp-X1))*1.0e-08_fp + 1.0_fp
     X2=(AS*AS - 1.0_fp)**2

     DY = (6.0_fp+3.0_fp*DELT)/(6.0_fp-7.0_fp*DELT)

     OPT = X2*DY/(WE**4)
  END SUBROUTINE RAYLO

END MODULE CRTM_MoleculeScatter
