!
! CRTM_MW_Land_SfcOptics
!
! Module to compute the surface optical properties for LAND surfaces at
! microwave frequencies required for determining the LAND surface
! contribution to the radiative transfer.
!
! This module is provided to allow developers to "wrap" their existing
! codes inside the provided functions to simplify integration into
! the main CRTM_SfcOptics module.
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 23-Jun-2005
!                       paul.vandelst@noaa.gov
!

MODULE CRTM_MW_Land_SfcOptics

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE Type_Kinds,               ONLY: fp
  USE Message_Handler,          ONLY: SUCCESS, FAILURE, Display_Message
  USE CRTM_Parameters,          ONLY: ZERO, ONE, MAX_N_ANGLES
  USE CRTM_SpcCoeff,            ONLY: SC
  USE CRTM_Surface_Define,      ONLY: CRTM_Surface_type
  USE CRTM_GeometryInfo_Define, ONLY: CRTM_GeometryInfo_type, &
                                      CRTM_GeometryInfo_GetValue
  USE CRTM_SfcOptics_Define,    ONLY: CRTM_SfcOptics_type
  USE NESDIS_LandEM_Module,     ONLY: NESDIS_LandEM
  USE CRTM_MWlandCoeff,         ONLY: MWlandC, CRTM_MWlandCoeff_IsLoaded
  USE TELSEM2_Atlas_Module,     ONLY: TELSEM2_Emissivity
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Data types
  PUBLIC :: iVar_type
  ! Science routines
  PUBLIC :: Compute_MW_Land_SfcOptics
  PUBLIC :: Compute_MW_Land_SfcOptics_TL
  PUBLIC :: Compute_MW_Land_SfcOptics_AD


  ! -----------------
  ! Module parameters
  ! -----------------
  ! Message length
  INTEGER, PARAMETER :: ML = 256
  ! Valid type indices for the microwave land emissivity model
  ! ...The soil types
  INTEGER, PARAMETER :: N_VALID_SOIL_TYPES = 8
  INTEGER, PARAMETER :: INVALID_SOIL    =  0
  INTEGER, PARAMETER :: COARSE          =  1
  INTEGER, PARAMETER :: MEDIUM          =  2
  INTEGER, PARAMETER :: FINE            =  3
  INTEGER, PARAMETER :: COARSE_MEDIUM   =  4
  INTEGER, PARAMETER :: COARSE_FINE     =  5
  INTEGER, PARAMETER :: MEDIUM_FINE     =  6
  INTEGER, PARAMETER :: COARSE_MED_FINE =  7
  INTEGER, PARAMETER :: ORGANIC         =  8
  ! ...The vegetation types
  INTEGER, PARAMETER :: N_VALID_VEGETATION_TYPES       = 12
  INTEGER, PARAMETER :: INVALID_VEGETATION             =  0
  INTEGER, PARAMETER :: BROADLEAF_EVERGREEN_TREES      =  1
  INTEGER, PARAMETER :: BROADLEAF_DECIDUOUS_TREES      =  2
  INTEGER, PARAMETER :: BROADLEAF_NEEDLELEAF_TREES     =  3
  INTEGER, PARAMETER :: NEEDLELEAF_EVERGREEN_TREES     =  4
  INTEGER, PARAMETER :: NEEDLELEAF_DECIDUOUS_TREES     =  5
  INTEGER, PARAMETER :: BROADLEAF_TREES_GROUNDCOVER    =  6
  INTEGER, PARAMETER :: GROUNDCOVER                    =  7
  INTEGER, PARAMETER :: GROADLEAF_SHRUBS_GROUNDCOVER   =  8
  INTEGER, PARAMETER :: BROADLEAF_SHRUBS_BARE_SOIL     =  9
  INTEGER, PARAMETER :: DWARF_TREES_SHRUBS_GROUNDCOVER = 10
  INTEGER, PARAMETER :: BARE_SOIL                      = 11
  INTEGER, PARAMETER :: CULTIVATIONS                   = 12


  ! --------------------------------------
  ! Structure definition to hold forward
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
  TYPE :: iVar_type
    PRIVATE
    ! Whether the low-frequency canopy model path was used. When .FALSE.
    ! (high-frequency default emissivity, or an invalid-type early return)
    ! the TL/AD results are zero.
    LOGICAL  :: Compute     = .FALSE.
    ! Whether the input vegetation fraction / soil moisture were clipped to [0,1]
    LOGICAL  :: Veg_Clipped = .FALSE.
    LOGICAL  :: Smc_Clipped = .FALSE.
    ! Forward values used to form vlai = Lai * Veg_Frac (Veg_Frac is post-clip)
    REAL(fp) :: Lai      = ZERO
    REAL(fp) :: Veg_Frac = ZERO
    ! Cached d(emissivity)/d(vlai) per angle: V is index 1, H is index 2
    REAL(fp), DIMENSION(MAX_N_ANGLES) :: dEV_dvlai = ZERO
    REAL(fp), DIMENSION(MAX_N_ANGLES) :: dEH_dvlai = ZERO
    ! Cached d(emissivity)/d(Soil_Moisture_Content) per angle
    REAL(fp), DIMENSION(MAX_N_ANGLES) :: dEV_dmv = ZERO
    REAL(fp), DIMENSION(MAX_N_ANGLES) :: dEH_dmv = ZERO
  END TYPE iVar_type


CONTAINS



!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       Compute_MW_Land_SfcOptics
!
! PURPOSE:
!       Function to compute the surface emissivity and reflectivity at microwave
!       frequencies over a land surface.
!
!       This function is a wrapper for third party code.
!
! CALLING SEQUENCE:
!       Error_Status = Compute_MW_Land_SfcOptics( &
!                        Surface     , &
!                        SensorIndex , &
!                        ChannelIndex, &
!                        SfcOptics     )
!
! INPUTS:
!       Surface:         CRTM_Surface structure containing the surface state
!                        data.
!                        UNITS:      N/A
!                        TYPE:       CRTM_Surface_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       GeometryInfo:    CRTM_GeometryInfo structure containing the
!                        view geometry information.
!                        UNITS:      N/A
!                        TYPE:       CRTM_GeometryInfo_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       SensorIndex:     Sensor index id. This is a unique index associated
!                        with a (supported) sensor used to access the
!                        shared coefficient data for a particular sensor.
!                        See the ChannelIndex argument.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:    Channel index id. This is a unique index associated
!                        with a (supported) sensor channel used to access the
!                        shared coefficient data for a particular sensor's
!                        channel.
!                        See the SensorIndex argument.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       SfcOptics:       CRTM_SfcOptics structure containing the surface
!                        optical properties required for the radiative
!                        transfer calculation. On input the Angle component
!                        is assumed to contain data.
!                        UNITS:      N/A
!                        TYPE:       CRTM_SfcOptics_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
! FUNCTION RESULT:
!       Error_Status:    The return value is an integer defining the error status.
!                        The error codes are defined in the Message_Handler module.
!                        If == SUCCESS the computation was sucessful
!                           == FAILURE an unrecoverable error occurred
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!
! COMMENTS:
!       Note the INTENT on the output SfcOptics argument is IN OUT rather
!       than just OUT as it is assumed to contain some data upon input.
!
!:sdoc-:
!----------------------------------------------------------------------------------

  FUNCTION Compute_MW_Land_SfcOptics( &
    Surface     , &  ! Input
    GeometryInfo, &  ! Input
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    SfcOptics   , &  ! Output
    iVar        ) &  ! Internal variable output
  RESULT ( err_stat )
    ! Arguments
    TYPE(CRTM_Surface_type),      INTENT(IN)     :: Surface
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)     :: GeometryInfo
    INTEGER,                      INTENT(IN)     :: SensorIndex
    INTEGER,                      INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_SfcOptics_type),    INTENT(IN OUT) :: SfcOptics
    TYPE(iVar_type),              INTENT(OUT)    :: iVar
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Compute_MW_Land_SfcOptics'
    REAL(fp),     PARAMETER :: FREQUENCY_CUTOFF   = 80.0_fp  ! GHz
    REAL(fp),     PARAMETER :: DEFAULT_EMISSIVITY = 0.95_fp
    ! Local variables
    CHARACTER(ML) :: msg
    INTEGER :: i
    INTEGER  :: month
    REAL(fp) :: lat, lon, ev, eh
    LOGICAL  :: atlas_valid, valid


    ! Set up
    err_stat = SUCCESS

    ! ----------------------------------------------------------------------
    ! TELSEM2 atlas path. When the microwave land emissivity atlas is loaded
    ! and has land climatology at this location/month, use it for all angles.
    ! The atlas depends only on lat/lon/month/frequency/angle (no CRTM control
    ! variable), so iVar%Compute is left .FALSE. and the TL/AD results are zero.
    ! A no-data cell (e.g. open water / permanent ice) falls through to the
    ! NESDIS_LandEM model below.
    ! ----------------------------------------------------------------------
    IF ( CRTM_MWlandCoeff_IsLoaded() ) THEN
      CALL CRTM_GeometryInfo_GetValue( GeometryInfo, &
                                       Latitude  = lat, &
                                       Longitude = lon, &
                                       Month     = month )
      atlas_valid = .TRUE.
      DO i = 1, SfcOptics%n_Angles
        CALL TELSEM2_Emissivity( MWlandC, lat, lon, month, &
                                 SC(SensorIndex)%Frequency(ChannelIndex), &
                                 SfcOptics%Angle(i), ev, eh, valid )
        IF ( .NOT. valid ) THEN
          atlas_valid = .FALSE.
          EXIT
        END IF
        SfcOptics%Emissivity(i,1) = ev
        SfcOptics%Emissivity(i,2) = eh
        ! Assume specular surface
        SfcOptics%Reflectivity(i,1,i,1) = ONE - ev
        SfcOptics%Reflectivity(i,2,i,2) = ONE - eh
      END DO
      IF ( atlas_valid ) RETURN  ! err_stat=SUCCESS; iVar%Compute stays .FALSE.
    END IF

    ! ...Check the soil type...
    IF ( Surface%Soil_Type < 1 .OR. &
         Surface%Soil_Type > N_VALID_SOIL_TYPES ) THEN
      SfcOptics%Emissivity   = ZERO
      SfcOptics%Reflectivity = ZERO
      err_stat = FAILURE
      msg = 'Invalid soil type index specified'
      CALL Display_Message( ROUTINE_NAME, msg, err_stat ); RETURN
    END IF
    ! ...and the vegetation type
    IF ( Surface%Vegetation_Type < 1 .OR. &
         Surface%Vegetation_Type > N_VALID_VEGETATION_TYPES ) THEN
      SfcOptics%Emissivity   = ZERO
      SfcOptics%Reflectivity = ZERO
      err_stat = FAILURE
      msg = 'Invalid vegetation type index specified'
      CALL Display_Message( ROUTINE_NAME, msg, err_stat ); RETURN
    END IF


    ! Compute the surface optical parameters
    IF ( SC(SensorIndex)%Frequency(ChannelIndex) < FREQUENCY_CUTOFF ) THEN
      ! Frequency is low enough for the model.
      ! ...Cache the forward state needed for the LAI/vegetation Jacobian.
      !    vlai = Lai*Vegetation_Fraction, with the vegetation fraction clipped
      !    to [0,1] inside NESDIS_LandEM (replicated here for the chain rule).
      iVar%Compute     = .TRUE.
      iVar%Lai         = Surface%Lai
      iVar%Veg_Frac    = MAX(MIN(Surface%Vegetation_Fraction,ONE),ZERO)
      iVar%Veg_Clipped = (Surface%Vegetation_Fraction < ZERO) .OR. &
                         (Surface%Vegetation_Fraction > ONE)
      iVar%Smc_Clipped = (Surface%Soil_Moisture_Content < ZERO) .OR. &
                         (Surface%Soil_Moisture_Content > ONE)
      DO i = 1, SfcOptics%n_Angles
        CALL NESDIS_LandEM(SfcOptics%Angle(i),            & ! Input, Degree
                           SC(SensorIndex)%Frequency(ChannelIndex),   & ! Input, GHz
                           Surface%Soil_Moisture_Content, & ! Input, g.cm^-3
                           Surface%Vegetation_Fraction,   & ! Input
                           Surface%Soil_Temperature,      & ! Input, K
                           Surface%Land_Temperature,      & ! Input, K
                           Surface%Lai,                   & ! Input, Leaf Area Index
                           Surface%Soil_Type,             & ! Input, Soil Type (1 -  9)
                           Surface%Vegetation_Type,       & ! Input, Vegetation Type (1 - 13)
                           ZERO,                          & ! Input, Snow depth, mm
                           SfcOptics%Emissivity(i,2),     & ! Output, H component
                           SfcOptics%Emissivity(i,1),     & ! Output, V component
                           dEV_dvlai = iVar%dEV_dvlai(i), & ! Optional output, V
                           dEH_dvlai = iVar%dEH_dvlai(i), & ! Optional output, H
                           dEV_dmv   = iVar%dEV_dmv(i),   & ! Optional output, V
                           dEH_dmv   = iVar%dEH_dmv(i)    ) ! Optional output, H
        ! Assume specular surface
        SfcOptics%Reflectivity(i,1,i,1) = ONE-SfcOptics%Emissivity(i,1)
        SfcOptics%Reflectivity(i,2,i,2) = ONE-SfcOptics%Emissivity(i,2)
      END DO
    ELSE
      ! Frequency is too high for model. Use default.
      DO i = 1, SfcOptics%n_Angles
        SfcOptics%Emissivity(i,1:2)         = DEFAULT_EMISSIVITY
        SfcOptics%Reflectivity(i,1:2,i,1:2) = ONE-DEFAULT_EMISSIVITY
      END DO
    END IF

  END FUNCTION Compute_MW_Land_SfcOptics


!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       Compute_MW_Land_SfcOptics_TL
!
! PURPOSE:
!       Function to compute the tangent-linear surface emissivity and
!       reflectivity at microwave frequencies over a land surface.
!
!       This function is a wrapper for third party code.
!
!       NB: CURRENTLY THIS IS A STUB FUNCTION AS THERE ARE NO TL
!           COMPONENTS IN THE MW LAND SFCOPTICS COMPUTATIONS.
!
! CALLING SEQUENCE:
!       Error_Status = Compute_MW_Land_SfcOptics_TL( SfcOptics_TL )
!
! OUTPUTS:
!       SfcOptics_TL:    Structure containing the tangent-linear surface
!                        optical properties required for the tangent-
!                        linear radiative transfer calculation.
!                        UNITS:      N/A
!                        TYPE:       CRTM_SfcOptics_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
! FUNCTION RESULT:
!       Error_Status:    The return value is an integer defining the error status.
!                        The error codes are defined in the Message_Handler module.
!                        If == SUCCESS the computation was sucessful
!                           == FAILURE an unrecoverable error occurred
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!
! COMMENTS:
!       Note the INTENT on the output SfcOptics_TL argument is IN OUT rather
!       than just OUT. This is necessary because the argument may be defined
!       upon input.
!
!:sdoc-:
!----------------------------------------------------------------------------------

  FUNCTION Compute_MW_Land_SfcOptics_TL( &
    SfcOptics   , &  ! FWD Input
    Surface_TL  , &  ! TL  Input
    SfcOptics_TL, &  ! TL  Output
    iVar        ) &  ! Internal variable input
  RESULT ( err_stat )
    ! Arguments
    TYPE(CRTM_SfcOptics_type), INTENT(IN)     :: SfcOptics
    TYPE(CRTM_Surface_type),   INTENT(IN)     :: Surface_TL
    TYPE(CRTM_SfcOptics_type), INTENT(IN OUT) :: SfcOptics_TL
    TYPE(iVar_type),           INTENT(IN)     :: iVar
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Compute_MW_Land_SfcOptics_TL'
    ! Local variables
    INTEGER  :: i
    REAL(fp) :: veg_frac_TL, vlai_TL, smc_TL


    ! Set up
    err_stat = SUCCESS
    SfcOptics_TL%Reflectivity = ZERO
    SfcOptics_TL%Emissivity   = ZERO

    ! No sensitivity unless the low-frequency canopy model path was taken
    IF ( .NOT. iVar%Compute ) RETURN

    ! Tangent-linear of vlai = Lai*Veg_Frac, with Veg_Frac clipped to [0,1]
    IF ( iVar%Veg_Clipped ) THEN
      veg_frac_TL = ZERO
    ELSE
      veg_frac_TL = Surface_TL%Vegetation_Fraction
    END IF
    vlai_TL = iVar%Veg_Frac*Surface_TL%Lai + iVar%Lai*veg_frac_TL

    ! Tangent-linear of soil moisture (clipped to [0,1] in the forward)
    IF ( iVar%Smc_Clipped ) THEN
      smc_TL = ZERO
    ELSE
      smc_TL = Surface_TL%Soil_Moisture_Content
    END IF

    ! Propagate to the surface emissivity/reflectivity (specular: r = 1 - e)
    DO i = 1, SfcOptics%n_Angles
      SfcOptics_TL%Emissivity(i,1) = iVar%dEV_dvlai(i)*vlai_TL + iVar%dEV_dmv(i)*smc_TL  ! V
      SfcOptics_TL%Emissivity(i,2) = iVar%dEH_dvlai(i)*vlai_TL + iVar%dEH_dmv(i)*smc_TL  ! H
      SfcOptics_TL%Reflectivity(i,1,i,1) = -SfcOptics_TL%Emissivity(i,1)
      SfcOptics_TL%Reflectivity(i,2,i,2) = -SfcOptics_TL%Emissivity(i,2)
    END DO

  END FUNCTION Compute_MW_Land_SfcOptics_TL



!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       Compute_MW_Land_SfcOptics_AD
!
! PURPOSE:
!       Function to compute the adjoint surface emissivity and
!       reflectivity at microwave frequencies over a land surface.
!
!       This function is a wrapper for third party code.
!
!       NB: CURRENTLY THIS IS A STUB FUNCTION AS THERE ARE NO AD
!           COMPONENTS IN THE MW LAND SFCOPTICS COMPUTATIONS.
!
! CALLING SEQUENCE:
!       Error_Status = Compute_MW_Land_SfcOptics_AD( SfcOptics_AD )
!
! INPUTS:
!       SfcOptics_AD:    Structure containing the adjoint surface optical
!                        properties required for the adjoint radiative
!                        transfer calculation.
!                        *** COMPONENTS MODIFIED UPON OUTPUT ***
!                        UNITS:      N/A
!                        TYPE:       CRTM_SfcOptics_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
! FUNCTION RESULT:
!       Error_Status:    The return value is an integer defining the error status.
!                        The error codes are defined in the Message_Handler module.
!                        If == SUCCESS the computation was sucessful
!                           == FAILURE an unrecoverable error occurred
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!
! COMMENTS:
!       Note the INTENT on the input adjoint arguments are IN OUT regardless
!       of their specification as "input" or "output". This is because these
!       arguments may contain information on input, or need to be zeroed on
!       output (or both).
!
!:sdoc-:
!----------------------------------------------------------------------------------

  FUNCTION Compute_MW_Land_SfcOptics_AD( &
    SfcOptics   , &  ! FWD Input
    SfcOptics_AD, &  ! AD  Input
    Surface_AD  , &  ! AD  Output
    iVar        ) &  ! Internal variable input
  RESULT( err_stat )
    ! Arguments
    TYPE(CRTM_SfcOptics_type),    INTENT(IN)     :: SfcOptics
    TYPE(CRTM_SfcOptics_type),    INTENT(IN OUT) :: SfcOptics_AD
    TYPE(CRTM_Surface_type),      INTENT(IN OUT) :: Surface_AD
    TYPE(iVar_type),              INTENT(IN)     :: iVar
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Compute_MW_Land_SfcOptics_AD'
    ! Local variables
    INTEGER  :: i
    REAL(fp) :: vlai_AD, smc_AD


    ! Set up
    err_stat = SUCCESS

    ! No sensitivity unless the low-frequency canopy model path was taken; the
    ! incoming adjoints are still consumed (zeroed) so they are not double counted.
    IF ( .NOT. iVar%Compute ) THEN
      SfcOptics_AD%Reflectivity = ZERO
      SfcOptics_AD%Emissivity   = ZERO
      RETURN
    END IF

    ! Adjoint of the emissivity/reflectivity -> vlai and soil moisture
    vlai_AD = ZERO
    smc_AD  = ZERO
    DO i = 1, SfcOptics%n_Angles
      ! Adjoint of specular reflectivity (r = 1 - e): e_AD += -r_AD, then zero r_AD
      SfcOptics_AD%Emissivity(i,1) = SfcOptics_AD%Emissivity(i,1) - SfcOptics_AD%Reflectivity(i,1,i,1)
      SfcOptics_AD%Emissivity(i,2) = SfcOptics_AD%Emissivity(i,2) - SfcOptics_AD%Reflectivity(i,2,i,2)
      SfcOptics_AD%Reflectivity(i,1,i,1) = ZERO
      SfcOptics_AD%Reflectivity(i,2,i,2) = ZERO
      ! Adjoint of emissivity = dE/dvlai * vlai + dE/dmv * soil_moisture
      vlai_AD = vlai_AD + iVar%dEV_dvlai(i)*SfcOptics_AD%Emissivity(i,1) &
                        + iVar%dEH_dvlai(i)*SfcOptics_AD%Emissivity(i,2)
      smc_AD  = smc_AD  + iVar%dEV_dmv(i)*SfcOptics_AD%Emissivity(i,1) &
                        + iVar%dEH_dmv(i)*SfcOptics_AD%Emissivity(i,2)
      SfcOptics_AD%Emissivity(i,1) = ZERO
      SfcOptics_AD%Emissivity(i,2) = ZERO
    END DO

    ! Adjoint of vlai = Lai*Veg_Frac (Veg_Frac clipped to [0,1])
    Surface_AD%Lai = Surface_AD%Lai + iVar%Veg_Frac*vlai_AD
    IF ( .NOT. iVar%Veg_Clipped ) THEN
      Surface_AD%Vegetation_Fraction = Surface_AD%Vegetation_Fraction + iVar%Lai*vlai_AD
    END IF

    ! Adjoint of soil moisture (clipped to [0,1] in the forward)
    IF ( .NOT. iVar%Smc_Clipped ) THEN
      Surface_AD%Soil_Moisture_Content = Surface_AD%Soil_Moisture_Content + smc_AD
    END IF

    ! Ensure no residual surface-optics adjoints leak downstream
    SfcOptics_AD%Reflectivity = ZERO
    SfcOptics_AD%Emissivity   = ZERO

  END FUNCTION Compute_MW_Land_SfcOptics_AD

END MODULE CRTM_MW_Land_SfcOptics
