!
! NESDIS_WaterIR_Emiss_V2_Module
!
! Module containing function to invoke the CSEM Infrared
! Sea Surface Emissivity Model (IRSSEM).
!
! This module was modfified and renamed from CRTM_IRSSEM for use in CSEM.
! The model coefficient files are loaded from the added model constructor 
! functions.

!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 22-Jun-2005
!                       paul.vandelst@noaa.gov
!
!       Modfied by:     Ming Chen, May-06-2020
!                       ming.chen@noaa.gov
!
MODULE NESDIS_WaterIR_Emiss_V2_Module

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  USE CSEM_Exception_Handler,   ONLY: SUCCESS, FAILURE, Display_Message

  USE CSEM_Interpolation, ONLY: NPTS        ,      &
                                CSEM_LPoly_type  , &
                                CSEM_find_index  , &
                                CSEM_interp_1D   , &
                                CSEM_interp_3D   , &
                                CSEM_interp_4D   , &
                                CSEM_Clear_LPoly , &
                                CSEM_LPoly , &      
                                CSEM_LPoly_TL, &
                                CSEM_LPoly_AD, &
                                CSEM_Interp_3D_TL, &
                                CSEM_Interp_3D_AD, &
                                CSEM_Interp_4D_TL, &
                                CSEM_Interp_4D_AD
                               
  USE IRSSEM_EmisCoeff_V2_Reader

  ! Disable implicit typing
  IMPLICIT NONE

  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Procedures
  PUBLIC :: NESDIS_WaterIR_Emiss_V2
  PUBLIC :: NESDIS_WaterIR_Emiss_V2_TL
  PUBLIC :: NESDIS_WaterIR_Emiss_V2_AD
  PUBLIC :: Einterp_V2_type
  PUBLIC :: IRSSEM_V2_Setup
  PUBLIC :: IRSSEM_V2_CleanUp
  PUBLIC :: IRSSEM_V2_Initialized
  
  ! -----------------
  ! Module parameters
  ! -----------------
  ! Version Id for the module
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: CSEM_IRSSEM.f90 21141 2012-09-14 17:40:43Z paul.vandelst@noaa.gov $'
  ! Message string length
  INTEGER, PARAMETER :: ML = 256
  INTEGER ,PARAMETER :: N_Angles    = 18
  REAL(fp),PARAMETER :: ONE = 1.0_fp, ZERO = 0.0_fp
  
  LOGICAL,SAVE :: IRSSEM_EmisCoeff_Init
  
  ! -------------------------------
  ! Structure definition to hold
  ! forward interpolating variables
  ! across fwd, tl and adjoint
  ! -------------------------------
  ! The interpolation routine structure
  TYPE :: Einterp_V2_type
    ! The dimensions
    INTEGER :: n_Pts    = NPTS
    ! The interpolating polynomials
    INTEGER :: N_Angles = 0   
    TYPE(CSEM_LPoly_type)    :: wlp(N_Angles)! Angle
    TYPE(CSEM_LPoly_type)    :: xlp       ! Frequency
    TYPE(CSEM_LPoly_type)    :: ylp       ! Wind Speed
    TYPE(CSEM_LPoly_type)    :: tlp       ! Temperature
    ! The LUT interpolation indices
    INTEGER           :: i1(N_Angles),i2(N_Angles)! Angle
    INTEGER           :: j1    , j2       ! Frequency
    INTEGER           :: k1    , k2       ! Wind Speed
    INTEGER           :: t1    , t2       ! Temperature
    ! The LUT interpolation boundary check
    LOGICAL           :: a_outbound(N_Angles)! Angle
    LOGICAL           :: f_outbound       ! Frequency
    LOGICAL           :: v_outbound       ! Wind Speed
    LOGICAL           :: t_outbound       ! Temperature
    ! The interpolation input
    REAL(fp)          :: a_int(N_Angles)  ! Angle
    REAL(fp)          :: f_int            ! Frequency
    REAL(fp)          :: v_int            ! Wind Speed
    REAL(fp)          :: t_int            ! Temperature
    ! The data to be interpolated
    REAL(fp)          :: a(NPTS,N_Angles) ! Angle
    REAL(fp)          :: f(NPTS)          ! Frequency
    REAL(fp)          :: v(NPTS)          ! Wind Speed
    REAL(fp)          :: t(NPTS)          ! Temperature
    REAL(fp)          :: e(NPTS, NPTS, NPTS, NPTS, N_Angles)! Emissivity
  END TYPE Einterp_V2_type

CONTAINS


!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

  
!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       NESDIS_WaterIR_Emiss_V2
!
! PURPOSE:
!       Function to compute the CSEM infrared sea surface emissivity (IRSSE)
!       for input wind speed, frequency, and angles.
!
! CALLING SEQUENCE:
!       Error_Status = NESDIS_WaterIR_Emiss_V2( Wind_Speed, &  ! Input
!                                           Frequency,  &  ! Input
!                                           Angle,      &  ! Input
!                                           Emissivity) &  ! Output
!                                               
!
! INPUTS:
!       Wind_Speed:     Wind speed.
!                       UNITS:      metres per second (m.s^-1)
!                       TYPE:       REAL(fp)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Frequency:      Infrared frequency.
!                       UNITS:      inverse centimetres (cm^-1)
!                       TYPE:       REAL(fp)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Angle:          Surface zenith angle.
!                       UNITS:      Degrees
!                       TYPE:       REAL(fp)
!                       DIMENSION:  Rank-1 (n_Angles)
!                       ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       Emissivity:     Sea surface emissivities for the
!                       requested wind speed, frequency, and angles.
!                       UNITS:      N/A
!                       TYPE:       REAL(fp)
!                       DIMENSION:  Same as input ANGLE argument.
!                       ATTRIBUTES: INTENT(OUT)
!
!
! FUNCTION RESULT:
!       Error_Status:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the emissivity computation was successful.
!                          == FAILURE an unrecoverable error occurred.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!
!:sdoc-:
!--------------------------------------------------------------------------------

  FUNCTION NESDIS_WaterIR_Emiss_V2( &
    Wind_Speed,  &  ! Input
    Temperature, &  ! Input
    Frequency ,  &  ! Input
    Angle     ,  &  ! Input
    Emissivity,  &  ! Output
    iVar)        &  ! Output
  RESULT( Error_Status )
    ! Arguments
    REAL(fp),          INTENT(IN)  :: Wind_Speed     ! v
    REAL(fp),          INTENT(IN)  :: Temperature    ! t
    REAL(fp),          INTENT(IN)  :: Frequency      ! f
    REAL(fp),          INTENT(IN)  :: Angle(:)       ! a
    REAL(fp),          INTENT(OUT) :: Emissivity(:)
    TYPE(Einterp_V2_type),INTENT(OUT) :: iVar
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'NESDIS_WaterIR_Emiss'
    ! Local variables
    !CHARACTER(ML) :: msg
    INTEGER :: i

    ! Set up
    Error_Status = SUCCESS
    IF (.NOT. IRSSEM_EmisCoeff_Init) THEN
      ! use default IRSSEM
      PRINT*, 'IRSSEM EmisCoff not loaded ...'
      RETURN
    ENDIF
    
    Error_Status = Compute_IRSSEM_Einterp( &
                   Wind_Speed,  &  ! Input
                   Temperature, &  ! Input
                   Frequency ,  &  ! Input
                   Angle     ,  &  ! Input
                   iVar      ) 
    ! Compute the interpolated emissivity
    DO i = 1, size(Angle)
      CALL CSEM_Interp_4D(iVar%e(:,:,:,:,i), & ! Input
                      iVar%wlp(i)        , & ! Input
                      iVar%xlp           , & ! Input
                      iVar%ylp           , & ! Input
                      iVar%tlp           , & ! Input
                      Emissivity(i)   )      ! Output
    END DO

  END FUNCTION NESDIS_WaterIR_Emiss_V2
!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       NESDIS_WaterIR_Emiss_V2_TL
!
! PURPOSE:
!       Function to compute the tangent linear of the CSEM infrared sea surface
!       emissivity (IRSSE) for input wind speed, frequency, and angles.
!
! CALLING SEQUENCE:
!       Error_Status = NESDIS_WaterIR_Emiss_V2_TL( &
!                                              Emissivity_TL, &  ! Input
!                                              Wind_Speed_TL)    ! Output
!
! INPUTS:
!
!!       Wind_Speed_TL:  Tangent linear wind speed.
!                       *** MUST HAVE VALUE ON ENTRY ***
!                       UNITS:      per metres per second, (m.s^-1)^-1
!                       TYPE:       REAL(fp)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!       iVar:           Structure containing internal variables required for
!                       subsequent tangent-linear or adjoint model calls.
!                       The contents of this structure are NOT accessible
!                       outside of the module.
!                       UNITS:      N/A
!                       TYPE:       iVar_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(INOUT)
!
! OUTPUTS:
!       Emissivity_TL:  Tangent linear sea surface emissivity.
!                       UNITS:      N/A
!                       TYPE:       REAL(fp)
!                       DIMENSION:  Rank-1 (n_Angles)
!                       ATTRIBUTES: INTENT(IN)
!
! FUNCTION RESULT:
!       Error_Status:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the computation was successful.
!                          == FAILURE an unrecoverable error occurred.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!
!:sdoc-:
!--------------------------------------------------------------------------------

  FUNCTION NESDIS_WaterIR_Emiss_V2_TL( &
    Wind_Speed_TL,  &  ! Input
    Temperature_TL, &  ! Input
    Emissivity_TL,  &  ! Output
    iVar        )   &  ! Input
  RESULT ( Error_Status )
    ! Arguments
    REAL(fp),   INTENT(IN)  ::  Wind_Speed_TL
    REAL(fp),   INTENT(IN)  ::  Temperature_TL
    REAL(fp),   INTENT(OUT) ::  Emissivity_TL(:)
    TYPE(Einterp_V2_type),   INTENT(IN) :: iVar 
    
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'NESDIS_WaterIR_Emiss_TL'
    ! Local variables
    !CHARACTER(ML) :: msg
    REAL(fp) :: v_TL(NPTS)
    REAL(fp) :: t_TL(NPTS)
    REAL(fp) :: e_TL(NPTS,NPTS,NPTS,NPTS)
    TYPE(CSEM_LPoly_Type) :: ylp_TL, xlp_TL, wlp_TL, tlp_TL
    INTEGER :: i

    ! Set up
    Error_Status = SUCCESS
    IF (.NOT. IRSSEM_EmisCoeff_Init) THEN
      ! use default IRSSEM
      PRINT*, 'IRSSEM EmisCoff not loaded ...'
      RETURN
   ENDIF
    
    ! ...Initialise local TL variables
    v_TL = ZERO
    t_TL = ZERO
    e_TL = ZERO
    CALL CSEM_Clear_LPoly(wlp_TL)
    CALL CSEM_Clear_LPoly(xlp_TL)
    CALL CSEM_Clear_LPoly(tlp_TL)
    IF(size(Emissivity_TL) /= iVar%n_Angles) THEN
       PRINT*, 'Inconsistent stream angles in IRSSEM ...'
       PRINT*, 'NESDIS_WaterIR_Emiss_TL',size(Emissivity_TL),iVar%n_Angles
       RETURN
    END IF 
 
    ! ...No TL if wind speed is out of bounds
    IF ( iVar%v_outbound ) THEN
       Emissivity_TL = ZERO
       RETURN
    END IF


    ! polynomials for wind speed
    CALL CSEM_LPoly_TL( iVar%v, iVar%v_int,   & ! FWD Input
                   iVar%ylp,                  & ! FWD Input
                   v_TL, Wind_Speed_TL,       & ! TL  Input
                   ylp_TL                    )  ! TL  Output
    ! polynomials for Temperature
    CALL CSEM_LPoly_TL( iVar%t, iVar%t_int,   & ! FWD Input
                   iVar%tlp,                  & ! FWD Input
                   t_TL, Temperature_TL,      & ! TL  Input
                   tlp_TL                    )  ! TL  Output

    DO i=1,iVar%n_Angles

      ! Perform interpolation
      CALL CSEM_interp_4D_TL( iVar%e(:,:,:,:,i), & ! FWD Emissivity input
                        iVar%wlp(i), & ! FWD polynomial input
                        iVar%xlp   , & ! FWD polynomial input
                        iVar%ylp   , & ! FWD polynomial input
                        iVar%tlp   , & ! FWD polynomial input
                        e_TL, wlp_TL, xlp_TL, ylp_TL, tlp_TL, & ! TL input
                        Emissivity_TL(i)              ) ! Output
    END DO

  END FUNCTION NESDIS_WaterIR_Emiss_V2_TL

  
 
!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       NESDIS_WaterIR_Emiss_V2_AD
!
! PURPOSE:
!       Function to compute the adjoint of the CSEM infrared sea surface
!       emissivity (IRSSE) for input wind speed, frequency, and angles.
!
! CALLING SEQUENCE:
!       Error_Status = NESDIS_WaterIR_Emiss_V2_AD( &
!                                              Emissivity_AD, &  ! Input
!                                              Wind_Speed_AD)    ! Output
!
! INPUTS:
!
!       Emissivity_AD:  Adjoint sea surface emissivity.
!                       *** SET TO ZERO ON EXIT ***
!                       UNITS:      N/A
!                       TYPE:       REAL(fp)
!                       DIMENSION:  Rank-1 (n_Angles)
!                       ATTRIBUTES: INTENT(IN OUT)
!
!       iVar:           Structure containing internal variables required for
!                       subsequent tangent-linear or adjoint model calls.
!                       The contents of this structure are NOT accessible
!                       outside of the CRTM_IRSSEM module.
!                       UNITS:      N/A
!                       TYPE:       iVar_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(INOUT)
!
! OUTPUTS:
!       Wind_Speed_AD:  Adjoint wind speed.
!                       *** MUST HAVE VALUE ON ENTRY ***
!                       UNITS:      per metres per second, (m.s^-1)^-1
!                       TYPE:       REAL(fp)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
! FUNCTION RESULT:
!       Error_Status:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the computation was successful.
!                          == FAILURE an unrecoverable error occurred.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!
!:sdoc-:
!--------------------------------------------------------------------------------

  FUNCTION NESDIS_WaterIR_Emiss_V2_AD ( &
    Emissivity_AD,  &  ! Input
    Wind_Speed_AD,  &  ! Output
    Temperature_AD, &  ! Output
    iVar        )  &  ! Input
  RESULT ( Error_Status )
    ! Arguments
    REAL(fp)          , INTENT(INOUT)  :: Emissivity_AD(:)
    REAL(fp)          , INTENT(INOUT)  :: Wind_Speed_AD
    REAL(fp)          , INTENT(INOUT)  :: Temperature_AD
    TYPE(Einterp_V2_type), INTENT(INOUT)  :: iVar
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'NESDIS_WaterIR_Emiss_AD'
    ! Local variables
    !CHARACTER(ML) :: msg
    INTEGER  :: i
    REAL(fp) :: v_AD(NPTS)
    REAL(fp) :: t_AD(NPTS)
    REAL(fp) :: e_AD(NPTS,NPTS,NPTS,NPTS)
    TYPE(CSEM_LPoly_Type) :: ylp_AD, xlp_AD, wlp_AD, tlp_AD

    ! Set up
    Error_Status = SUCCESS
    IF (.NOT. IRSSEM_EmisCoeff_Init) THEN
      ! use default IRSSEM
      PRINT*, 'IRSSEM EmisCoff not loaded ...'
      RETURN
    ENDIF
    ! ...Initialise local AD variables
    v_AD = ZERO
    t_AD = ZERO
    e_AD = ZERO
    !Wind_Speed_AD = ZERO !!!
    
     ! ...Initialize local variables
    CALL CSEM_Clear_LPoly(wlp_AD)
    CALL CSEM_Clear_LPoly(xlp_AD)
    CALL CSEM_Clear_LPoly(ylp_AD)
    CALL CSEM_Clear_LPoly(tlp_AD)
    IF(size(Emissivity_AD) /= iVar%n_Angles) THEN
       PRINT*, 'Inconsistent stream angles in IRSSEM ...'
       PRINT*, 'NESDIS_WaterIR_Emiss_AD',size(Emissivity_AD),iVar%n_Angles
      RETURN
    END IF 
    DO i=1,iVar%n_Angles
       IF ( iVar%v_outbound ) THEN
         Wind_Speed_AD  = Wind_Speed_AD + ZERO
         Temperature_AD =Temperature_AD + ZERO
         RETURN
       END IF
    
      ! Perform interpolation
       CALL CSEM_Interp_4D_AD(iVar%e(:,:,:,:,i), & ! FWD Emissivity input
                   iVar%wlp(i)  , & ! FWD Input
                   iVar%xlp     , & ! FWD Input
                   iVar%ylp     , & ! FWD Input
                   iVar%tlp     , & ! FWD Input
                   Emissivity_AD(i), & ! AD Input
                   e_AD, wlp_AD, xlp_AD, ylp_AD, tlp_AD) ! AD Output
    END DO
    ! Compute the wind speed adjoint
    CALL CSEM_Lpoly_AD(iVar%v,     & ! FWD Input
                  iVar%v_int,      & ! FWD Input
                  iVar%ylp,        & ! FWD Input
                  ylp_AD,          & ! AD  Input
                  v_AD,            & ! AD  Output
                  Wind_Speed_AD)     ! AD  Output
    ! Compute the Temperature adjoint
    CALL CSEM_Lpoly_AD(iVar%t,     & ! FWD Input
                  iVar%t_int,      & ! FWD Input
                  iVar%tlp,        & ! FWD Input
                  tlp_AD,          & ! AD  Input
                  t_AD,            & ! AD  Output
                  Temperature_AD  )  ! AD  Output

  END FUNCTION NESDIS_WaterIR_Emiss_V2_AD
   
  

!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       Compute_IRSSEM_Einterp
!
! PURPOSE:
!       Function to interpolate the CRTM infrared sea surface emissivity (IRSSE)
!       for input wind speed, frequency, and angles.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_Compute_IRSSEM( Wind_Speed, &  ! Input
!                                           Frequency,  &  ! Input
!                                           Angle,      &  ! Input
!                                           Emissivity, &  ! Output
!                                           iVar        )  ! Internal Variable Output
!
! INPUTS:
!       Wind_Speed:     Wind speed.
!                       UNITS:      metres per second (m.s^-1)
!                       TYPE:       REAL(fp)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Frequency:      Infrared frequency.
!                       UNITS:      inverse centimetres (cm^-1)
!                       TYPE:       REAL(fp)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Angle:          Surface zenith angle.
!                       UNITS:      Degrees
!                       TYPE:       REAL(fp)
!                       DIMENSION:  Rank-1 (n_Angles)
!                       ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!
!       iVar:           Structure containing internal variables required for
!                       subsequent tangent-linear or adjoint model calls.
!                       The contents of this structure are NOT accessible
!                       outside of the CRTM_IRSSEM module.
!                       UNITS:      N/A
!                       TYPE:       iVar_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(OUT)
!
! FUNCTION RESULT:
!       Error_Status:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the emissivity computation was successful.
!                          == FAILURE an unrecoverable error occurred.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!
!:sdoc-:
!--------------------------------------------------------------------------------

  FUNCTION Compute_IRSSEM_Einterp( &
    Wind_Speed,  &  ! Input
    Temperature, &  ! Input
    Frequency ,  &  ! Input
    Angle     ,  &  ! Input
    iVar      )  &  ! Internal variable output
  RESULT( Error_Status )
    ! Arguments
    REAL(fp)       , INTENT(IN)  :: Wind_Speed     ! v
    REAL(fp)       , INTENT(IN)  :: Temperature    ! t
    REAL(fp)       , INTENT(IN)  :: Frequency      ! f
    REAL(fp)       , INTENT(IN)  :: Angle(:)       ! a
    TYPE(Einterp_V2_type), INTENT(OUT) :: iVar
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Interpolate_IRSSEM'
    ! Local variables
    !CHARACTER(256):: msg
    REAL(fp)     , SAVE :: Wind_Speed_Saved           = -999.0      ! v
    REAL(fp)     , SAVE :: Frequency_Saved            = -999.0      ! f
    REAL(fp)     , SAVE :: Temperature_Saved          = -999.0      ! t
    REAL(fp)     , SAVE :: Angle_Saved(N_Angles)      = -999.0      ! a
    TYPE(Einterp_V2_type), SAVE :: iVar_Saved
    
    INTEGER :: i
    ! Set up
    Error_Status = SUCCESS
    iVar  = iVar_Saved
 
    ! Compute the wind speed interpolating polynomial
    ! ...Find the LUT indices and check if input is out of bounds
    IF(Wind_Speed /= Wind_Speed_Saved) THEN  
       iVar%v_int = Wind_Speed
       CALL CSEM_find_index(IRwaterC_V2%Wind_Speed, iVar%v_int, &
            iVar%k1, iVar%k2, iVar%v_outbound)
       iVar%v = IRwaterC_V2%Wind_Speed(iVar%k1:iVar%k2)
       ! ...Compute the polynomial
       CALL CSEM_LPoly(iVar%v,  & ! Input
                   iVar%v_int,  & ! Input
                   iVar%ylp    )  ! Output
       Wind_Speed_Saved = Wind_Speed
    ENDIF
    IF(Temperature /= Temperature_Saved) THEN  
       iVar%t_int = Temperature
       CALL CSEM_find_index(IRwaterC_V2%Temperature, iVar%t_int, &
            iVar%t1, iVar%t2, iVar%t_outbound)
       iVar%t = IRwaterC_V2%Temperature(iVar%t1:iVar%t2)
       ! ...Compute the polynomial
       CALL CSEM_LPoly(iVar%t,  & ! Input
                   iVar%t_int,  & ! Input
                   iVar%tlp    )  ! Output
       Temperature_Saved = Temperature
    ENDIF

    ! Compute the frequency interpolating polynomial
    ! ...Find the LUT indices and check if input is out of bounds
    IF( Frequency /= Frequency_Saved) THEN
       iVar%f_int = Frequency
       CALL CSEM_find_index(IRwaterC_V2%Frequency, iVar%f_int, &
             iVar%j1, iVar%j2, iVar%f_outbound)
       iVar%f = IRwaterC_V2%Frequency(iVar%j1:iVar%j2)
       ! ...Compute the polynomial
       CALL CSEM_LPoly(iVar%f, & ! Input
                   iVar%f_int, & ! Input
                   iVar%xlp    ) ! Output
       Frequency_Saved = Frequency
    ENDIF

    ! Compute the angle interpolating polynomials
    ! ...Find the LUT indices and check if input is out of bounds
    iVar%n_Angles = size(Angle)
    DO i=1, iVar%n_Angles 
       iVar%a_int(i) = ABS(Angle(i))
       IF(iVar%a_int(i) /= Angle_Saved(i)) THEN
         CALL CSEM_find_index(IRwaterC_V2%Angle, iVar%a_int(i), &
              iVar%i1(i), iVar%i2(i), iVar%a_outbound(i))
         iVar%a(:,i)  = IRwaterC_V2%Angle(iVar%i1(i):iVar%i2(i))
         ! ...Compute the polynomial
         CALL CSEM_LPoly(iVar%a(:,i),   & ! Input
                  iVar%a_int(i),        & ! Input
                  iVar%wlp(i)    )        ! Output
         Angle_Saved(i) = iVar%a_int(i)
       END IF  
       ! Compute the interpolated emissivity 
       iVar%e(:,:,:,:,i) = IRwaterC_V2%Emissivity( iVar%i1(i) :iVar%i2(i) , &
                      iVar%j1   :iVar%j2   , &
                      iVar%k1   :iVar%k2, iVar%t1   :iVar%t2 )

    ENDDO
   
    iVar_Saved = iVar

  END FUNCTION Compute_IRSSEM_Einterp


  FUNCTION IRSSEM_V2_Setup( Coeff_File_Name) RESULT( Error_Status )

    CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: Coeff_File_Name         ! IRSSEM netcdf coeff file name
    
    ! Function result
    INTEGER :: Error_Status
     

    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'IRSSEM_Setup'
    ! Local variables
    CHARACTER(256) :: Message
    CHARACTER(LEN=256) :: IRWaterCoeff_File

    IF (PRESENT(Coeff_File_Name)) THEN
       IRWaterCoeff_File= TRIM(Coeff_File_Name)
    ELSE
       IRWaterCoeff_File='Nalli.IRwater.EmisCoeff.nc'
    END IF
    
    ! ------  
    ! Set up  
    ! ------  
    Error_Status = SUCCESS
    
    IF (IRSSEM_EmisCoeff_Init) THEN
        Message = 'IRSSEM_LUT already initialised.'
        PRINT*,Message
        RETURN
    ENDIF
    
    !PRINT*,'Initialising IRSSEM_EmisCoeff...'
    !print*,IRWaterCoeff_File
    
    Error_Status=Load_IRSSEM_V2_LUT(TRIM(IRWaterCoeff_File))
    IF (Error_Status /= SUCCESS) THEN
        Error_Status = FAILURE
        PRINT*, 'Error initialising IRSSEM_LUT..'
        RETURN
    END IF
    IRSSEM_EmisCoeff_Init= .TRUE.


  END FUNCTION IRSSEM_V2_Setup
  
  FUNCTION IRSSEM_V2_Initialized() RESULT( LUT_Status )
    LOGICAL :: LUT_Status
    LUT_Status = IRSSEM_EmisCoeff_Init
  END FUNCTION IRSSEM_V2_Initialized

  FUNCTION IRSSEM_V2_CleanUP()  RESULT(Error_Status)
    INTEGER :: Error_Status
    Error_Status = Close_IRSSEM_V2_LUT()
    IRSSEM_EmisCoeff_Init= .FALSE.
  END FUNCTION IRSSEM_V2_CleanUP
  
END MODULE NESDIS_WaterIR_Emiss_V2_Module
