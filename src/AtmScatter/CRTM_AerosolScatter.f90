
! CRTM_AerosolScatter
!
!
! Module to compute the aerosol absorption and scattering properties
! required for radiative transfer in an atmosphere with aerosols.
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 15-Feb-2005
!                       paul.vandelst@noaa.gov
!       Modified by     Quanhua Liu, 03-Oct-2006
!                       Quanhua.Liu@noaa.gov
!       Modified by:    James Rosinski, 08-Feb-2019
!                       Rosinski@ucar.edu
!       Modified by     Yingtao Ma, 2020/6/11; yingtao.ma@noaa.gov
!                       Implemented CMAQ aerosol (based on Quanhua's v2.0.2).
!       Modified by:    Yingtao Ma, 2020/8/13; yingtao.ma@noaa.gov
!                       Rewrote the TL/AD code to correctly handle size variance.
!       Modified by     Cheng Dang, 20-Nov-2020;
!                       dangch@ucar.edu
!                       Update FWD/TL/AD to include RH dimension and GOCART-GEOS5 table
!                       14-Mar-2021; Update to include NAAPS table
! (C) Copyright 2019 UCAR
!


MODULE CRTM_AerosolScatter

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE Type_Kinds,               ONLY: fp
  USE Message_Handler,          ONLY: SUCCESS, FAILURE, Display_Message
  USE String_Utility,           ONLY: StrUpCase
  USE CRTM_Parameters,          ONLY: ZERO, ONE, POINT_5, ONEpointFIVE, &
                                      AEROSOL_CONTENT_THRESHOLD, &
                                      HGPHASE  ! <<< NEED TO REMOVE THIS IN FUTURE
  USE CRTM_SpcCoeff,            ONLY: SC, &
                                      SpcCoeff_IsMicrowaveSensor , &
                                      SpcCoeff_IsInfraredSensor  , &
                                      SpcCoeff_IsVisibleSensor   , &
                                      SpcCoeff_IsUltravioletSensor
  USE CRTM_AerosolCoeff,        ONLY: AeroC
  USE AerosolCoeff_Define,      ONLY: AerosolCoeff_typeID_to_index, &
                                      AerosolCoeff_typeName_to_index
  USE CRTM_Atmosphere_Define,   ONLY: CRTM_Atmosphere_type, &
                                      AerosolCoeff_BYPASS_AEROSOL
  USE CRTM_GeometryInfo_Define, ONLY: CRTM_GeometryInfo_type
  USE CRTM_Interpolation,       ONLY: NPTS        , &
                                      LPoly_type  , &
                                      find_index  , &
                                      interp_1D   , &
                                      interp_2D   , &
                                      interp_3D   , &
                                      interp_4D   , &
                                      interp_2D_TL, &
                                      interp_3D_TL, &
                                      interp_4D_TL, &
                                      interp_2D_AD, &
                                      interp_3D_AD, &
                                      interp_4D_AD, &
                                      Clear_LPoly , &
                                      LPoly       , &
                                      LPoly_TL    , &
                                      LPoly_AD
  USE CRTM_AtmOptics_Define,    ONLY: CRTM_AtmOptics_type

  ! Internal variable definition module
  USE ASvar_Define, ONLY: ASvar_type, &
                          ASinterp_type


  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Procedures
  PUBLIC :: CRTM_Compute_AerosolScatter
  PUBLIC :: CRTM_Compute_AerosolScatter_TL
  PUBLIC :: CRTM_Compute_AerosolScatter_AD


  ! -----------------
  ! Module parameters
  ! -----------------
  ! Message string length
  INTEGER, PARAMETER :: ML = 256
  ! Number of stream angle definitions
  INTEGER, PARAMETER :: FOUR_STREAMS      =  4
  INTEGER, PARAMETER :: SIX_STREAMS       =  6
  INTEGER, PARAMETER :: EIGHT_STREAMS     =  8
  INTEGER, PARAMETER :: SIXTEEN_STREAMS   = 16


CONTAINS


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Compute_AerosolScatter
!
! PURPOSE:
!       Function to compute the aerosol absorption and scattering properties
!       and populate the output AerosolScatter structure for a single channel.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_Compute_AerosolScatter( Atmosphere    , &
!                                                   SensorIndex   , &
!                                                   ChannelIndex  , &
!                                                   AerosolScatter, &
!                                                   ASvar           )
!
! INPUT ARGUMENTS:
!       Atmosphere:      CRTM_Atmosphere structure containing the atmospheric
!                        profile data.
!                        UNITS:      N/A
!                        TYPE:       CRTM_Atmosphere_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:    Channel index id. This is a unique index associated
!                        with a (supported) sensor channel used to access the
!                        shared coefficient data.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       SensorIndex:     Sensor index id. This ia a unique index associated
!                        with a sensor.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       AerosolScatter:  CRTM_AtmOptics structure containing the aerosol
!                        absorption and scattering properties required by
!                        the radiative transfer.
!                        UNITS:      N/A
!                        TYPE:       CRTM_AtmOptics_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
!       ASvar:           Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        UNITS:      N/A
!                        TYPE:       ASvar_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(OUT)
!
! FUNCTION RESULT:
!       Error_Status:    The return value is an integer defining the error status.
!                        The error codes are defined in the ERROR_HANDLER module.
!                        If == SUCCESS the computation was sucessful
!                           == FAILURE an unrecoverable error occurred
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_Compute_AerosolScatter( &
    Atm         , &  ! Input
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    AScat       , &  ! Output
    ASV         ) &  ! Internal variable output
  RESULT( Error_Status )
    !Arguments
    TYPE(CRTM_Atmosphere_type), INTENT(IN)     :: Atm
    INTEGER                   , INTENT(IN)     :: ChannelIndex
    INTEGER                   , INTENT(IN)     :: SensorIndex
    TYPE(CRTM_AtmOptics_type) , INTENT(IN OUT) :: AScat
    TYPE(ASvar_type)          , INTENT(IN OUT) :: ASV
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Compute_AerosolScatter'
    ! Local Variables
    CHARACTER(ML) :: Message
    INTEGER :: k, ka, l, m, n
    REAL(fp) :: Frequency
    LOGICAL  :: Layer_Mask(Atm%n_Layers)
    INTEGER  :: Layer_Index(Atm%n_Layers)
    INTEGER  :: nAerosol_Layers
    REAL(fp) :: bs

   ! ------
   ! Set up
   ! ------
   Error_Status = SUCCESS
   IF ( SpcCoeff_IsMicrowaveSensor(SC(SensorIndex)) ) RETURN
   IF ( Atm%n_Aerosols == 0 ) RETURN
   ASV%Total_bs = ZERO
   ! Frequency in inverse centimetres
   Frequency = SC(SensorIndex)%Wavenumber(ChannelIndex)
   ! Determine offset for Legendre coefficients in the
   ! AeroC lookup table corresponding to the
   ! number of streams
   SELECT CASE(AScat%n_Legendre_Terms)
      CASE (FOUR_STREAMS)   ; AScat%lOffset = 0
      CASE (SIX_STREAMS)    ; AScat%lOffset = 5
      CASE (EIGHT_STREAMS)  ; AScat%lOffset = 12
      CASE (SIXTEEN_STREAMS); AScat%lOffset = 21
      CASE DEFAULT
        AScat%lOffset = 0
        ! Use two-stream model or HG and RAYLEIGH Phase function
        IF( HGPHASE ) THEN
          AScat%n_Legendre_Terms = 0
        ELSE
          Error_Status = FAILURE
          WRITE(Message,'("The n_Legendre_Terms in AerosolScatter, ",i0,", do not fit model")') &
                        AScat%n_Legendre_Terms
          CALL Display_Message(ROUTINE_NAME, Message, Error_Status )
          RETURN
        END IF
   END SELECT

    ! -----------------------------------------------
    ! Loop over the different Aerosols in the profile
    ! -----------------------------------------------
    Aerosol_loop: DO n = 1, Atm%n_Aerosols

      ! Only process aerosols with more
      ! than the threshold Aerosol amount
      Layer_Mask  = Atm%Aerosol(n)%Concentration > AEROSOL_CONTENT_THRESHOLD
      nAerosol_Layers = COUNT(Layer_Mask)
      IF ( nAerosol_Layers == 0 ) CYCLE Aerosol_loop

      ! --------------------------------------
      ! Loop over the current Aerosol's Layers
      ! --------------------------------------
      Layer_Index(1:nAerosol_Layers) = PACK((/(k,k=1,Atm%Aerosol(n)%n_Layers)/), Layer_Mask)
      Aerosol_Layer_loop: DO k = 1, nAerosol_Layers
        ka = Layer_Index(k)

        ! Get the aerosol optical properties
        CALL Get_Aerosol_Opt( AScat                                 , & ! Input
                              Frequency                             , & ! Input
                              Atm%Aerosol(n)%Type                   , & ! Input
                              Atm%Aerosol(n)%Effective_Radius(ka)   , & ! Input
                              Atm%Aerosol(n)%Effective_Variance(ka) , & ! Input
                              Atm%Relative_Humidity(ka)             , & ! Input
                              ASV%ke(ka,n)                          , & ! Output
                              ASV%w(ka,n)                           , & ! Output
                              ASV%g(ka,n)                           , & ! Output
                              ASV%kb(ka,n)                          , & ! Output
                              ASV%pcoeff(:,:,ka,n)                  , & ! Output
                              ASV%asi(ka,n)                        ) ! Interpolation

        ! interpolation quality control
        IF( ASV%ke(ka,n) <= ZERO ) THEN
          ASV%ke(ka,n) = ZERO
          ASV%w(ka,n)  = ZERO
        END IF
        IF( ASV%w(ka,n) <= ZERO ) THEN
          ASV%w(ka,n) = ZERO
          ASV%pcoeff(:,:,ka,n) = ZERO
        END IF
        IF( ASV%w(ka,n) >= ONE ) THEN
          ASV%w(ka,n) = ONE
        END IF

        ! Compute the optical depth (absorption + scattering)
        !   tau = rho.ke
        ! where
        !   rho = Integrated Aerosol Concentration for a layer(kg/m^2) [M.L^-2]
        !   ke  = mass extintion coefficient (m^2/kg) [L^2.M^-1]
        ! Note that since all these computations are done for a given
        ! layer, the optical depth is the same as the volume extinction
        ! coefficient, be. Usually,
        !   tau = be.d(z)
        ! but we are working with height/thickness independent quantities
        ! so that
        !   tau = be
        ! This is why the optical depth is used in the denominator to
        ! compute the single scatter albedo in the Layer_loop below.
        AScat%Optical_Depth(ka) = AScat%Optical_Depth(ka) + &
                                  (ASV%ke(ka,n)*Atm%Aerosol(n)%Concentration(ka))

        ! Compute the phase matrix coefficients
        ! p = p + p(LUT)*bs
        ! where
        !   p(LUT) = the phase coefficient from the LUT
        IF( AScat%n_Phase_Elements > 0 .and. AScat%Include_Scattering ) THEN
        ! Compute the volume scattering coefficient for the current
        ! aerosol layer and accumulate it for the layer total for the
        ! profile (i.e. all aerosols)
        !   bs = rho.w.ke
        ! where
        !   bs  = volume scattering coefficient for a layer [dimensionless]
        !   rho = integrated aerosol concentration for a layer (kg/m^2) [M.L^-2]
        !   w   = single scatter albedo [dimensionless]
        !   ke  = mass extintion coefficient (m^2/kg) [L^2.M^-1]
        !  qliu
        bs = Atm%Aerosol(n)%Concentration(ka) * ASV%ke(ka,n) * ASV%w(ka,n)
        ASV%Total_bs(ka) = ASV%Total_bs(ka) + bs
        AScat%Single_Scatter_Albedo(ka) = AScat%Single_Scatter_Albedo(ka) + bs
        AScat%Asymmetry_Factor(ka) = AScat%Asymmetry_Factor(ka) + (ASV%g(ka,n) * bs)
        AScat%Backscat_Coefficient(ka) = AScat%Backscat_Coefficient(ka) + &
                                         (ASV%kb(ka,n) * Atm%Aerosol(n)%Concentration(ka))

          DO m = 1, AScat%n_Phase_Elements
            DO l = 0, AScat%n_Legendre_Terms
              AScat%Phase_Coefficient(l,m,ka) = AScat%Phase_Coefficient(l,m,ka) + &
                                                (ASV%pcoeff(l,m,ka,n) * bs)
            END DO
          END DO
        END IF
      END DO Aerosol_Layer_loop
    END DO Aerosol_loop


  END FUNCTION CRTM_Compute_AerosolScatter


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Compute_AerosolScatter_TL
!
! PURPOSE:
!       Function to compute the tangent-linear aerosol absorption and
!       scattering properties and populate the output AerosolScatter_TL
!       structure for a single channel.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_Compute_AerosolScatter_TL( Atmosphere       , &
!                                                      AerosolScatter   , &
!                                                      Atmosphere_TL    , &
!                                                      SensorIndex      , &
!                                                      ChannelIndex     , &
!                                                      AerosolScatter_TL, &
!                                                      ASvar              )
!
! INPUT ARGUMENTS:
!       Atmosphere:         CRTM_Atmosphere structure containing the atmospheric
!                           profile data.
!                           UNITS:      N/A
!                           TYPE:       CRTM_Atmosphere_type
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
!
!       AerosolScatter:     CRTM_AtmOptics structure containing the forward model
!                           aerosol absorption and scattering properties required
!                           for radiative transfer.
!                           UNITS:      N/A
!                           TYPE:       CRTM_AtmOptics_type
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
!
!       Atmosphere_TL:      CRTM Atmosphere structure containing the tangent-linear
!                           atmospheric state data.
!                           UNITS:      N/A
!                           TYPE:       CRTM_Atmosphere_type
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
!
!       SensorIndex:        Sensor index id. This ia a unique index associated
!                           with a sensor.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
!
!
!       Channel_Index:      Channel index id. This is a unique index associated
!                           with a (supported) sensor channel used to access the
!                           shared coefficient data.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
!
!       ASvar:              Structure containing internal variables required for
!                           subsequent tangent-linear or adjoint model calls.
!                           UNITS:      N/A
!                           TYPE:       ASvar_type
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       AerosolScatter_TL:  CRTM_AtmOptics structure containing the tangent-linear
!                           aerosol absorption and scattering properties required
!                           for radiative transfer.
!                           UNITS:      N/A
!                           TYPE:       CRTM_AtmOptics_type
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN OUT)
!
! FUNCTION RESULT:
!       Error_Status:       The return value is an integer defining the error status.
!                           The error codes are defined in the ERROR_HANDLER module.
!                           If == SUCCESS the computation was sucessful
!                              == FAILURE an unrecoverable error occurred
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_Compute_AerosolScatter_TL( &
    Atm         , & ! FWD Input
    AScat       , & ! FWD Input
    Atm_TL      , & ! TL  Input
    SensorIndex , & ! Input
    ChannelIndex, & ! Input
    AScat_TL    , & ! TL  Input
    ASV         ) & ! Internal variable input
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_Atmosphere_type) , INTENT(IN)     :: Atm
    TYPE(CRTM_AtmOptics_type)  , INTENT(IN)     :: AScat
    TYPE(CRTM_Atmosphere_type) , INTENT(IN)     :: Atm_TL
    INTEGER                    , INTENT(IN)     :: SensorIndex
    INTEGER                    , INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_AtmOptics_type)  , INTENT(IN OUT) :: AScat_TL
    TYPE(ASvar_type)           , INTENT(IN)     :: ASV
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Compute_AerosolScatter_TL'
    ! Local variables
    INTEGER  :: k, ka, l, m, n
    INTEGER  :: n_Legendre_Terms, n_Phase_Elements
    REAL(fp) :: Frequency
    LOGICAL  :: Layer_Mask(Atm%n_Layers)
    INTEGER  :: Layer_Index(Atm%n_Layers)
    INTEGER  :: nAerosol_Layers
    REAL(fp) :: ke_TL, w_TL, g_TL, kb_TL
    REAL(fp) :: pcoeff_TL(0:AScat%n_Legendre_Terms, AScat%n_Phase_Elements)
    REAL(fp) :: bs, bs_TL

    ! ------
    ! Set up
    ! ------
    Error_Status = SUCCESS
    IF ( SpcCoeff_IsMicrowaveSensor(SC(SensorIndex)) ) RETURN
    IF ( Atm%n_Aerosols == 0 ) RETURN
    ! Frequency
    Frequency = SC(SensorIndex)%Wavenumber(ChannelIndex)
    ! Phase matrix dimensions
    n_Legendre_Terms = AScat_TL%n_Legendre_Terms
    n_Phase_Elements = AScat_TL%n_Phase_Elements
    AScat_TL%lOffset = AScat%lOffset


    ! -----------------------------------------------
    ! Loop over the different Aerosols in the profile
    ! -----------------------------------------------
    Aerosol_loop: DO n = 1, Atm%n_Aerosols

      ! Only process aerosols with more
      ! than the threshold aerosol amount
      Layer_Mask = Atm%Aerosol(n)%Concentration > AEROSOL_CONTENT_THRESHOLD
      nAerosol_Layers = COUNT(Layer_Mask)
      IF (nAerosol_Layers == 0) CYCLE Aerosol_loop


      ! --------------------------------------
      ! Loop over the current aerosol's layers
      ! --------------------------------------
      Layer_Index(1:nAerosol_Layers) = PACK((/(k,k=1,Atm%Aerosol(n)%n_Layers)/), Layer_Mask)
      Aerosol_Layer_loop: DO k = 1, nAerosol_Layers
        ka = Layer_Index(k)

        CALL Get_Aerosol_Opt_TL(AScat_TL                                , & ! Input
                                Atm%Aerosol(n)%Type                     , & ! Input
                                ASV%ke(ka,n)                            , & ! Input
                                ASV%w(ka,n)                             , & ! Input
                                Atm_TL%Aerosol(n)%Effective_Radius(ka)  , & ! TL  Input
                                Atm_TL%Aerosol(n)%Effective_Variance(ka), & ! TL  Input
                                Atm_TL%Relative_Humidity(ka)            , & ! TL  Input
                                ke_TL                                   , & ! TL  Output
                                w_TL                                    , & ! TL  Output
                                g_TL                                    , & ! TL  Output
                                kb_TL                                   , & ! TL  Output
                                pcoeff_TL                               , & ! TL  Output
                                ASV%asi(ka,n)                           ) ! Interpolation

        ! interpolation quality control
        IF( ASV%ke(ka,n) <= ZERO ) THEN
          ke_TL = ZERO
          w_TL = ZERO
        END IF
        IF( ASV%w(ka,n) <= ZERO ) THEN
          w_TL  = ZERO
          pcoeff_TL = ZERO
        END IF
        IF( ASV%w(ka,n) >= ONE ) THEN
          w_TL  = ZERO
        END IF

        ! Compute the optical depth (absorption + scattering)
        AScat_TL%Optical_Depth(ka) = AScat_TL%Optical_Depth(ka) + &
                                     (ke_TL        * Atm%Aerosol(n)%Concentration(ka)) + &
                                     (ASV%ke(ka,n) * Atm_TL%Aerosol(n)%Concentration(ka))

        ! Compute the phase matrix coefficients
        IF( n_Phase_Elements > 0 .and. AScat%Include_Scattering ) THEN
          ! Compute the volume scattering coefficient
          bs = Atm%Aerosol(n)%Concentration(ka) * ASV%ke(ka,n) * ASV%w(ka,n)
          bs_TL = (Atm_TL%Aerosol(n)%Concentration(ka) * ASV%ke(ka,n) * ASV%w(ka,n) ) + &
                   (Atm%Aerosol(n)%Concentration(ka) * ke_TL * ASV%w(ka,n) ) + &
                   (Atm%Aerosol(n)%Concentration(ka) * ASV%ke(ka,n) * w_TL )
          AScat_TL%Single_Scatter_Albedo(ka) = AScat_TL%Single_Scatter_Albedo(ka) + bs_TL
          AScat_TL%Asymmetry_Factor(ka) = AScat_TL%Asymmetry_Factor(ka) + &
                                          (g_TL * bs) + (ASV%g(ka,n) * bs_TL)
          AScat_TL%Backscat_Coefficient(ka) = AScat_TL%Backscat_Coefficient(ka) + &
                                              (kb_TL * Atm%Aerosol(n)%Concentration(ka)) + &
                                              (ASV%kb(ka,n) * Atm_TL%Aerosol(n)%Concentration(ka))

          DO m = 1, n_Phase_Elements
            DO l = 0, n_Legendre_Terms
              AScat_TL%Phase_Coefficient(l,m,ka) = AScat_TL%Phase_Coefficient(l,m,ka) + &
                                                   (pcoeff_TL(l,m)       * bs   ) + &
                                                   (ASV%pcoeff(l,m,ka,n) * bs_TL)
            END DO
          END DO
        END IF
      END DO Aerosol_Layer_loop
    END DO Aerosol_loop

  END FUNCTION CRTM_Compute_AerosolScatter_TL


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Compute_AerosolScatter_AD
!
! PURPOSE:
!       Function to compute the adjoint aerosol absorption and scattering
!       properties for a single channel.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_Compute_AerosolScatter_AD( Atmosphere       , &
!                                                      AerosolScatter   , &
!                                                      AerosolScatter_AD, &
!                                                      SensorIndex      , &
!                                                      ChannelIndex     , &
!                                                      Atmosphere_AD    , &
!                                                      ASVariables        )
!
! INPUT ARGUMENTS:
!       Atmosphere:         CRTM_Atmosphere structure containing the atmospheric
!                           profile data.
!                           UNITS:      N/A
!                           TYPE:       CRTM_Atmosphere_type
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
!
!       AerosolScatter:     CRTM_AtmOptics structure containing the forward model
!                           aerosol absorption and scattering properties required
!                           for radiative transfer.
!                           UNITS:      N/A
!                           TYPE:       CRTM_AtmOptics_type
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
!
!       AerosolScatter_AD:  CRTM_AtmOptics structure containing the adjoint
!                           aerosol absorption and scattering properties.
!                           **NOTE: On EXIT from this function, the contents of
!                                   this structure may be modified (e.g. set to
!                                   zero.)
!                           UNITS:      N/A
!                           TYPE:       CRTM_AtmOptics_type
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN OUT)
!
!       SensorIndex:        Sensor index id. This ia a unique index associated
!                           with a sensor.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
!
!       Channel_Index:      Channel index id. This is a unique index associated
!                           with a (supported) sensor channel used to access the
!                           shared coefficient data.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
!
!       ASvar:              Structure containing internal variables required for
!                           subsequent tangent-linear or adjoint model calls.
!                           UNITS:      N/A
!                           TYPE:       ASvar_type
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       Atmosphere_AD:      CRTM Atmosphere structure containing the adjoint
!                           atmospheric state data.
!                           UNITS:      N/A
!                           TYPE:       CRTM_Atmosphere_type
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN OUT)
!
!
! FUNCTION RESULT:
!       Error_Status:       The return value is an integer defining the error
!                           status. The error codes are defined in the
!                           ERROR_HANDLER module.
!                           If == SUCCESS the computation was sucessful
!                              == FAILURE an unrecoverable error occurred
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_Compute_AerosolScatter_AD( &
    Atm         , &  ! FWD Input
    AScat       , &  ! FWD Input
    AScat_AD    , &  ! AD  Input
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    Atm_AD      , &  ! AD  Output
    ASV         ) &  ! Internal Variable input
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_Atmosphere_type), INTENT(IN)     :: Atm
    TYPE(CRTM_AtmOptics_type) , INTENT(IN)     :: AScat
    TYPE(CRTM_AtmOptics_type) , INTENT(IN OUT) :: AScat_AD
    INTEGER                   , INTENT(IN)     :: SensorIndex
    INTEGER                   , INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_Atmosphere_type), INTENT(IN OUT) :: Atm_AD
    TYPE(ASvar_type)          , INTENT(IN)     :: ASV
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Compute_AerosolScatter_AD'
    ! Local variables
    INTEGER   :: k, ka, l, m, n
    INTEGER   :: n_Legendre_Terms, n_Phase_Elements
    REAL(fp)  :: Frequency
    LOGICAL   :: Layer_Mask(Atm%n_Layers)
    INTEGER   :: Layer_Index(Atm%n_Layers)
    INTEGER   :: nAerosol_Layers
    REAL(fp)  :: ke_AD, w_AD, g_AD, kb_AD
    REAL(fp)  :: pcoeff_AD(0:AScat%n_Legendre_Terms, AScat%n_Phase_Elements)
    REAL(fp)  :: bs, bs_AD

    ! ------
    ! Set up
    ! ------
    Error_Status = SUCCESS
    IF ( SpcCoeff_IsMicrowaveSensor(SC(SensorIndex)) ) RETURN
    IF ( Atm%n_Aerosols == 0 ) RETURN
    ! Frequency
    Frequency = SC(SensorIndex)%Wavenumber(ChannelIndex)
    ! Phase matrix dimensions
    n_Legendre_Terms = AScat_AD%n_Legendre_Terms
    n_Phase_Elements = AScat_AD%n_Phase_Elements
    AScat_AD%lOffset = AScat%lOffset


    ! ----------------------------------------------------------
    ! Adjoint of accumulated optical properties for all aerosols
    ! ----------------------------------------------------------
    ! Shorten variable name
    l = n_Legendre_Terms

    ! ------------------------------------------------
    ! Loop over different types of aerosols in profile
    ! ------------------------------------------------
    Aerosol_loop: DO n = 1, Atm%n_Aerosols

      ! Only process aerosols with more than
      ! the threshold aerosol concentration
      Layer_Mask = Atm%Aerosol(n)%Concentration > AEROSOL_CONTENT_THRESHOLD
      nAerosol_Layers = COUNT(Layer_Mask)
      IF ( nAerosol_Layers == 0 ) CYCLE Aerosol_loop

      ! --------------------------------------
      ! Loop over the current aerosol's layers
      ! --------------------------------------
      Layer_Index(1:nAerosol_Layers) = PACK((/(k,k=1,Atm%Aerosol(n)%n_Layers)/), Layer_Mask)
      Aerosol_Layer_loop: DO k = 1, nAerosol_Layers
        ka = Layer_Index(k)

        ! Initialize the individual
        ! Aerosol adjoint variables
        bs_AD = ZERO
        pcoeff_AD = ZERO
        ke_AD     = ZERO
        w_AD      = ZERO
        g_AD      = ZERO
        kb_AD     = ZERO

        ! Compute the adjoint of the
        ! phase matrix coefficients
        IF( n_Phase_Elements > 0 .and. AScat%Include_Scattering ) THEN
          ! Recompute the forward model volume scattering
          ! coefficient for the current aerosol type ONLY
          bs = Atm%Aerosol(n)%Concentration(ka) * ASV%ke(ka,n) * ASV%w(ka,n)
          DO m = 1, n_Phase_Elements
            DO l = 0, n_Legendre_Terms
              bs_AD = bs_AD + (ASV%pcoeff(l,m,ka,n) * AScat_AD%Phase_Coefficient(l,m,ka))
              pcoeff_AD(l,m) = pcoeff_AD(l,m) + (bs * AScat_AD%Phase_Coefficient(l,m,ka))
            END DO
          END DO
          ! Adjoint of asymmetry factor
          g_AD  = g_AD  + (bs * AScat_AD%Asymmetry_Factor(ka))
          bs_AD = bs_AD + (ASV%g(ka,n) * AScat_AD%Asymmetry_Factor(ka))

          ! NOTE: bs_AD is not reinitialized after this
          !       point since it is reinitialized at the
          !       start of the Aerosol_Layer_loop
          bs_AD = bs_AD + AScat_AD%Single_Scatter_Albedo(ka)
          w_AD  = w_AD  + (Atm%Aerosol(n)%Concentration(ka) * ASV%ke(ka,n)* bs_AD )
        END IF

        ! Compute the adjoint of the optical
        ! depth (absorption + scattering)
        Atm_AD%Aerosol(n)%Concentration(ka) = Atm_AD%Aerosol(n)%Concentration(ka) + &
                                              (ASV%ke(ka,n) * AScat_AD%Optical_Depth(ka))
        ke_AD = ke_AD + (Atm%Aerosol(n)%Concentration(ka) * AScat_AD%Optical_Depth(ka))

        ! Compute the adjoint of the backscatter coefficient
        Atm_AD%Aerosol(n)%Concentration(ka) = Atm_AD%Aerosol(n)%Concentration(ka) + &
                                              (ASV%kb(ka,n) * AScat_AD%Backscat_Coefficient(ka))
        kb_AD = kb_AD + (Atm%Aerosol(n)%Concentration(ka) * AScat_AD%Backscat_Coefficient(ka))

        ! Compute the adjoint of the volume
        ! scattering coefficient.

        ke_AD = ke_AD + (Atm%Aerosol(n)%Concentration(ka) * bs_AD * ASV%w(ka,n) )
        Atm_AD%Aerosol(n)%Concentration(ka) = Atm_AD%Aerosol(n)%Concentration(ka) + &
                                              ( bs_AD * ASV%ke(ka,n) * ASV%w(ka,n) )

        ! interpolation quality control
        IF( ASV%w(ka,n) >= ONE ) THEN
          w_AD = ZERO
        END IF
        IF( ASV%ke(ka,n) <= ZERO ) THEN
          ke_AD = ZERO
          w_AD = ZERO
        END IF
        IF( ASV%w(ka,n) <= ZERO ) THEN
          w_AD = ZERO
          pcoeff_AD = ZERO
        END IF

        CALL Get_Aerosol_Opt_AD(AScat_AD                                , & ! Input
                                Atm%Aerosol(n)%Type                     , & ! Input
                                ASV%ke(ka,n)                            , & ! input
                                ASV%w(ka,n)                             , & ! input
                                ke_AD                                   , & ! AD Input
                                w_AD                                    , & ! AD Input
                                g_AD                                    , & ! AD Input
                                kb_AD                                   , & ! AD Input
                                pcoeff_AD                               , & ! AD Input
                                Atm_AD%Aerosol(n)%Effective_Radius(ka)  , & ! AD  Output
                                Atm_AD%Aerosol(n)%Effective_Variance(ka), & ! AD  Output
                                Atm_AD%Relative_Humidity(ka)            , & ! AD  Output
                                ASV%asi(ka,n)                           )   ! Interpolation

      END DO Aerosol_Layer_loop
    END DO Aerosol_loop

  END FUNCTION CRTM_Compute_AerosolScatter_AD



!################################################################################
!################################################################################
!##                                                                            ##
!##                        ## PRIVATE MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################


  ! ----------------------------------------
  ! Subroutine to obtain the bulk
  ! optical properties of aerosol
  ! extinction coefficient (ke),
  ! scattering coefficient (w)
  ! asymmetry factor (g), and
  ! spherical Legendre coefficients (pcoeff)
  ! ----------------------------------------
  SUBROUTINE Get_Aerosol_Opt( AerosolScatter, &  ! Input AerosolScatter structure
                              Frequency     , &  ! Input in cm^-1
                              Aerosol_Type  , &  ! Input see CRTM_Aerosol_Define.f90
                              Reff          , &  ! Input effective radius (mm)
                              Rsig          , &  ! Input effective variance, unitless
                              RH            , &  ! Input relative humidity, unitless
                              ke            , &  ! Output extinction coefficient (=~ optical depth)
                              w             , &  ! Output single scattering albedo
                              g             , &  ! Output asymmetry factor
                              kb            , &  ! Output backscatter coefficient
                              pcoeff        , &  ! Output spherical Legendre coefficients
                              asi             )  ! Output interpolation data
    ! Arguments
    TYPE(CRTM_AtmOptics_type) , INTENT(IN)     :: AerosolScatter
    REAL(fp)                  , INTENT(IN)     :: Frequency
    INTEGER                   , INTENT(IN)     :: Aerosol_Type
    REAL(fp)                  , INTENT(IN)     :: Reff, Rsig
    REAL(fp)                  , INTENT(IN)     :: RH
    REAL(fp)                  , INTENT(OUT)    :: ke
    REAL(fp)                  , INTENT(OUT)    :: w
    REAL(fp)                  , INTENT(OUT)    :: g
    REAL(fp)                  , INTENT(OUT)    :: kb
    REAL(fp)                  , INTENT(IN OUT) :: pcoeff(0:,:)
    TYPE(ASinterp_type)       , INTENT(IN OUT) :: asi
    ! Local variables
    INTEGER  :: k, l, m

    ! Get the aerosol type LUT index
    ! ------------------------------
    ! Bypass aerosol
    IF ( Aerosol_Type == AerosolCoeff_BYPASS_AEROSOL ) THEN
      ke = ZERO
      w  = ZERO
      g  = ZERO
      kb = ZERO
      pcoeff = ZERO
      RETURN
    END IF

    k = AerosolCoeff_typeID_to_index( AeroC, Aerosol_Type )

    ! Find the Frequency indices for interpolation (Always Dim 1)
    asi%f_int = MAX( MIN( AeroC%Frequency( AeroC%n_Wavelengths ), Frequency ), AeroC%Frequency(1) )
    CALL find_index( AeroC%Frequency(:), asi%f_int, asi%i1, asi%i2, asi%f_outbound )
    asi%f = AeroC%Frequency( asi%i1:asi%i2 )
    CALL LPoly( asi%f, asi%f_int, asi%wlp )

    ! Relative Humidity (Dim 2)
    IF ( AeroC%n_RH > 1 ) THEN
      asi%h_int = MAX( MIN( AeroC%RH( AeroC%n_RH ), RH ), AeroC%RH(1) )
      CALL find_index( AeroC%RH(:), asi%h_int, asi%h1, asi%h2, asi%h_outbound )
      asi%h = AeroC%RH( asi%h1:asi%h2 )
      CALL LPoly( asi%h, asi%h_int, asi%hlp )
    ELSE
      asi%h1 = 1 ; asi%h2 = 1 ; asi%h_outbound = .FALSE.
    END IF

    ! Effective Radius (Dim 3)
    IF ( AeroC%n_Radii > 1 ) THEN
      asi%r_int = MAX( MIN( AeroC%Reff( AeroC%n_Radii, k ), Reff ), AeroC%Reff(1, k) )
      CALL find_index( AeroC%Reff(:, k), asi%r_int, asi%j1, asi%j2, asi%r_outbound )
      asi%r = AeroC%Reff( asi%j1:asi%j2, k )
      CALL LPoly( asi%r, asi%r_int, asi%xlp )
    ELSE
      asi%j1 = 1 ; asi%j2 = 1 ; asi%r_outbound = .FALSE.
    END IF

    ! Size Variance (Dim 4)
    IF ( AeroC%n_Sigma > 1 ) THEN
      asi%v_int = MAX( MIN( AeroC%Rsig( AeroC%n_Sigma, k ), Rsig ), AeroC%Rsig(1, k) )
      CALL find_index( AeroC%Rsig(:, k), asi%v_int, asi%k1, asi%k2, asi%v_outbound )
      asi%v = AeroC%Rsig( asi%k1:asi%k2, k )
      CALL LPoly( asi%v, asi%v_int, asi%vlp )
    ELSE
      asi%k1 = 1 ; asi%k2 = 1 ; asi%v_outbound = .FALSE.
    END IF

    ! Perform 4D Interpolation
    CALL Interp_4D( AeroC%ke(asi%i1:asi%i2, asi%h1:asi%h2, asi%j1:asi%j2, asi%k1:asi%k2, k), &
                    asi%wlp, asi%hlp, asi%xlp, asi%vlp, ke )
    CALL Interp_4D( AeroC%w(asi%i1:asi%i2, asi%h1:asi%h2, asi%j1:asi%j2, asi%k1:asi%k2, k), &
                    asi%wlp, asi%hlp, asi%xlp, asi%vlp, w )
    CALL Interp_4D( AeroC%g(asi%i1:asi%i2, asi%h1:asi%h2, asi%j1:asi%j2, asi%k1:asi%k2, k), &
                    asi%wlp, asi%hlp, asi%xlp, asi%vlp, g )
    CALL Interp_4D( AeroC%kb(asi%i1:asi%i2, asi%h1:asi%h2, asi%j1:asi%j2, asi%k1:asi%k2, k), &
                    asi%wlp, asi%hlp, asi%xlp, asi%vlp, kb )

    IF ( AerosolScatter%n_Phase_Elements > 0 .AND. AerosolScatter%Include_Scattering ) THEN
      pcoeff(0,1) = POINT_5
      DO m = 1, AerosolScatter%n_Phase_Elements
        DO l = 1, AerosolScatter%n_Legendre_Terms
          CALL Interp_4D( AeroC%pcoeff(asi%i1:asi%i2, asi%h1:asi%h2, asi%j1:asi%j2, asi%k1:asi%k2, k, l+AerosolScatter%lOffset, m), &
                          asi%wlp, asi%hlp, asi%xlp, asi%vlp, pcoeff(l,m) )
        END DO
      END DO
    ELSE
      ! Absorption coefficient
      ke = ke * (ONE - w)
    END IF

  END SUBROUTINE Get_Aerosol_Opt


  ! ---------------------------------------------
  ! Subroutine to obtain the tangent-linear
  !  bulk optical properties of a aerosol:
  !   extinction coefficient (ke_TL),
  !   scattereing coefficient (w_TL)
  !   asymmetry factor (g_TL), and
  !   spherical Legendre coefficients (pcoeff_TL)
  ! ---------------------------------------------

  SUBROUTINE Get_Aerosol_Opt_TL(AerosolScatter_TL, &  ! Input  AerosolScatterTL structure
                                Aerosol_Type     , &  ! Input  see CRTM_Aerosol_Define.f90
                                ke               , &  ! Input
                                w                , &  ! Input
                                Reff_TL          , &  ! Input  TL effective radius (mm)
                                Rsig_TL          , &  ! Input  TL effective radius (mm)
                                RH_TL            , &  ! Input  TL effective radius (mm)
                                ke_TL            , &  ! Output TL extinction coefficient (=~ optical depth)
                                w_TL             , &  ! Output TL single scattering albedo
                                g_TL             , &  ! Output TL asymmetry factor
                                kb_TL            , &  ! Output TL backscatter coefficient
                                pcoeff_TL        , &  ! TL  Output spherical Legendre coefficients
                                asi                )  ! Input interpolation data

    ! Arguments
    TYPE(CRTM_AtmOptics_type),  INTENT(IN)     :: AerosolScatter_TL
    INTEGER ,                   INTENT(IN)     :: Aerosol_Type
    REAL(fp),                   INTENT(IN)     :: ke, w
    REAL(fp),                   INTENT(IN)     :: Reff_TL, Rsig_TL, RH_TL
    REAL(fp),                   INTENT(OUT)    :: ke_TL, w_TL, g_TL, kb_TL
    REAL(fp),                   INTENT(IN OUT) :: pcoeff_TL(0:,:)
    TYPE(ASinterp_type),        INTENT(IN)     :: asi
    ! Local variables
    INTEGER  :: k, l, m
    REAL(fp) :: f_int_TL, r_int_TL, v_int_TL, h_int_TL
    REAL(fp) :: f_TL(NPTS), r_TL(NPTS), v_TL(NPTS), h_TL(NPTS)
    REAL(fp) :: zg_TL(NPTS,NPTS,NPTS,NPTS)
    TYPE(LPoly_type) :: wlp_TL, xlp_TL, vlp_TL, hlp_TL
    REAL(fp), POINTER :: z4(:,:,:,:)

    ! Bypass aerosol
    IF ( Aerosol_Type == AerosolCoeff_BYPASS_AEROSOL ) THEN
      ke_TL     = ZERO
      w_TL      = ZERO
      g_TL      = ZERO
      kb_TL     = ZERO
      pcoeff_TL = ZERO
      RETURN
    END IF

    ! Setup
    ! -----
    ! No TL output when all dimensions are outside LUT bounds
    IF ( asi%f_outbound .AND. asi%h_outbound .AND. asi%r_outbound .AND. asi%v_outbound ) THEN
      ke_TL     = ZERO
      w_TL      = ZERO
      g_TL      = ZERO
      kb_TL     = ZERO
      pcoeff_TL = ZERO
      RETURN
    END IF

    ! The TL inputs
    f_int_TL = ZERO ; f_TL = ZERO
    h_int_TL = RH_TL ; h_TL = ZERO
    r_int_TL = Reff_TL ; r_TL = ZERO
    v_int_TL = Rsig_TL ; v_TL = ZERO
    zg_TL = ZERO

    ! Calculate the TL interpolating polynomials
    CALL LPoly_TL( asi%f, asi%f_int, asi%wlp, f_TL, f_int_TL, wlp_TL )
    CALL LPoly_TL( asi%h, asi%h_int, asi%hlp, h_TL, h_int_TL, hlp_TL )
    CALL LPoly_TL( asi%r, asi%r_int, asi%xlp, r_TL, r_int_TL, xlp_TL )
    CALL LPoly_TL( asi%v, asi%v_int, asi%vlp, v_TL, v_int_TL, vlp_TL )

    ! Get the aerosol type LUT index
    k = AerosolCoeff_typeID_to_index( AeroC, Aerosol_Type )

    ! Extinction coefficient
    z4 => AeroC%ke(asi%i1:asi%i2, asi%h1:asi%h2, asi%j1:asi%j2, asi%k1:asi%k2, k)
    CALL Interp_4D_TL( z4, asi%wlp, asi%hlp, asi%xlp, asi%vlp, &
                       zg_TL, wlp_TL, hlp_TL, xlp_TL, vlp_TL, ke_TL )
    ! Single scatter albedo
    z4 => AeroC%w(asi%i1:asi%i2, asi%h1:asi%h2, asi%j1:asi%j2, asi%k1:asi%k2, k)
    CALL Interp_4D_TL( z4, asi%wlp, asi%hlp, asi%xlp, asi%vlp, &
                       zg_TL, wlp_TL, hlp_TL, xlp_TL, vlp_TL, w_TL )
    ! Asymmetry factor
    z4 => AeroC%g(asi%i1:asi%i2, asi%h1:asi%h2, asi%j1:asi%j2, asi%k1:asi%k2, k)
    CALL Interp_4D_TL( z4, asi%wlp, asi%hlp, asi%xlp, asi%vlp, &
                       zg_TL, wlp_TL, hlp_TL, xlp_TL, vlp_TL, g_TL )
    ! Backscatter coefficient
    z4 => AeroC%kb(asi%i1:asi%i2, asi%h1:asi%h2, asi%j1:asi%j2, asi%k1:asi%k2, k)
    CALL Interp_4D_TL( z4, asi%wlp, asi%hlp, asi%xlp, asi%vlp, &
                       zg_TL, wlp_TL, hlp_TL, xlp_TL, vlp_TL, kb_TL )

    ! Phase matrix coefficients
    IF (AerosolScatter_TL%n_Phase_Elements > 0 .AND. AerosolScatter_TL%Include_Scattering) THEN
      pcoeff_TL(0,1) = ZERO
      DO m = 1, AerosolScatter_TL%n_Phase_Elements
        DO l = 1, AerosolScatter_TL%n_Legendre_Terms
          z4 => AeroC%pcoeff(asi%i1:asi%i2, asi%h1:asi%h2, asi%j1:asi%j2, asi%k1:asi%k2, k, l+AerosolScatter_TL%lOffset, m)
          CALL Interp_4D_TL( z4, asi%wlp, asi%hlp, asi%xlp, asi%vlp, &
                             zg_TL, wlp_TL, hlp_TL, xlp_TL, vlp_TL, pcoeff_TL(l,m) )
        END DO
      END DO
    ELSE
      ! Absorption coefficient
      IF ( w < ONE ) THEN
        ke_TL = ke_TL * (ONE - w) - (ke / (ONE - w)) * w_TL
      ELSE
        ke_TL = ZERO
      END IF
    END IF
    NULLIFY(z4)

  END SUBROUTINE Get_Aerosol_Opt_TL


  ! ---------------------------------------------
  ! Subroutine to obtain the adjoint
  !  bulk optical properties of a aerosol:
  !   extinction coefficient (ke_AD),
  !   scattereing coefficient (w_AD)
  !   asymmetry factor (g_AD), and
  !   spherical Legendre coefficients (pcoeff_AD)
  ! ---------------------------------------------
  SUBROUTINE Get_Aerosol_Opt_AD( AerosolScatter_AD, & ! Input AerosolScatter AD structure
                                 Aerosol_Type     , & ! Input see CRTM_Aerosol_Define.f90
                                 ke               , & ! Input
                                 w                , & ! Input
                                 ke_AD            , & ! AD Input extinction cross section
                                 w_AD             , & ! AD Input single scatter albedo
                                 g_AD             , & ! AD Input asymmetry factor
                                 kb_AD            , & ! AD Input backscatter coefficient
                                 pcoeff_AD        , & ! AD Input spherical Legendre coefficients
                                 Reff_AD          , & ! AD Outputmode radius
                                 Rsig_AD          , & ! AD Output radius deviation
                                 RH_AD            , & ! AD Output relative humidity
                                 asi                ) ! Input interpolation data
    ! Arguments
    TYPE(CRTM_AtmOptics_type),  INTENT(IN)     :: AerosolScatter_AD
    INTEGER ,                   INTENT(IN)     :: Aerosol_Type
    REAL(fp),                   INTENT(IN)     :: ke, w
    REAL(fp),                   INTENT(IN OUT) :: ke_AD            ! AD Input
    REAL(fp),                   INTENT(IN OUT) :: w_AD             ! AD Input
    REAL(fp),                   INTENT(IN OUT) :: g_AD             ! AD Input
    REAL(fp),                   INTENT(IN OUT) :: kb_AD            ! AD Input
    REAL(fp),                   INTENT(IN OUT) :: pcoeff_AD(0:,:)  ! AD Input
    REAL(fp),                   INTENT(IN OUT) :: Reff_AD          ! AD Output
    REAL(fp),                   INTENT(IN OUT) :: Rsig_AD          ! AD Output
    REAL(fp),                   INTENT(IN OUT) :: RH_AD            ! AD Output
    TYPE(ASinterp_type),        INTENT(IN)     :: asi
    ! Local variables
    INTEGER  :: k, l, m
    REAL(fp) :: f_int_AD, r_int_AD, v_int_AD, h_int_AD
    REAL(fp) :: f_AD(NPTS), r_AD(NPTS), v_AD(NPTS), h_AD(NPTS)
    REAL(fp) :: z4_AD(NPTS,NPTS,NPTS,NPTS)
    TYPE(LPoly_type) :: wlp_AD, xlp_AD, vlp_AD, hlp_AD
    REAL(fp), POINTER :: z4(:,:,:,:)

    ! Bypass aerosol
    IF ( Aerosol_Type == AerosolCoeff_BYPASS_AEROSOL ) THEN
      Reff_AD   = ZERO ; Rsig_AD   = ZERO ; RH_AD     = ZERO
      ke_AD     = ZERO ; w_AD      = ZERO ; g_AD      = ZERO ; kb_AD     = ZERO
      pcoeff_AD = ZERO
      RETURN
    END IF

    ! Setup
    ! -----
    ! No AD output all dimensions are outside LUT bounds
    IF ( asi%f_outbound .AND. asi%h_outbound .AND. asi%r_outbound .AND. asi%v_outbound ) THEN
      Reff_AD   = ZERO ; Rsig_AD   = ZERO ; RH_AD     = ZERO
      ke_AD     = ZERO ; w_AD      = ZERO ; g_AD      = ZERO ; kb_AD     = ZERO
      pcoeff_AD = ZERO
      RETURN
    END IF

    ! Initialise local adjoint variables
    f_int_AD = ZERO ; r_int_AD = ZERO ; v_int_AD = ZERO ; h_int_AD = ZERO
    f_AD = ZERO ; r_AD = ZERO ; v_AD = ZERO ; h_AD = ZERO
    z4_AD = ZERO
    CALL Clear_LPoly(wlp_AD)
    CALL Clear_LPoly(hlp_AD)
    CALL Clear_LPoly(xlp_AD)
    CALL Clear_LPoly(vlp_AD)

    ! Get the aerosol type LUT index
    k = AerosolCoeff_typeID_to_index( AeroC, Aerosol_Type )

    ! Phase matrix coefficients
    IF (AerosolScatter_AD%n_Phase_Elements > 0 .AND. AerosolScatter_AD%Include_Scattering ) THEN
      DO m = 1, AerosolScatter_AD%n_Phase_Elements
        DO l = 1, AerosolScatter_AD%n_Legendre_Terms
          z4 => AeroC%pcoeff(asi%i1:asi%i2, asi%h1:asi%h2, asi%j1:asi%j2, asi%k1:asi%k2, k, l+AerosolScatter_AD%lOffset, m)
          CALL Interp_4D_AD( z4, asi%wlp, asi%hlp, asi%xlp, asi%vlp, pcoeff_AD(l,m), &
                             z4_AD, wlp_AD, hlp_AD, xlp_AD, vlp_AD )
        END DO
      END DO
      pcoeff_AD(0,1) = ZERO
    ELSE
      ! Absorption coefficient
      IF ( w < ONE ) THEN
        w_AD  = w_AD - (ke / (ONE - w)) * ke_AD
        ke_AD = ke_AD * (ONE - w)
      ELSE
        ke_AD = ZERO
      END IF
    END IF

    ! Backscatter coefficient
    z4 => AeroC%kb(asi%i1:asi%i2, asi%h1:asi%h2, asi%j1:asi%j2, asi%k1:asi%k2, k)
    CALL Interp_4D_AD( z4, asi%wlp, asi%hlp, asi%xlp, asi%vlp, kb_AD, &
                       z4_AD, wlp_AD, hlp_AD, xlp_AD, vlp_AD )
    kb_AD = ZERO

    ! Asymmetry factor
    z4 => AeroC%g(asi%i1:asi%i2, asi%h1:asi%h2, asi%j1:asi%j2, asi%k1:asi%k2, k)
    CALL Interp_4D_AD( z4, asi%wlp, asi%hlp, asi%xlp, asi%vlp, g_AD, &
                       z4_AD, wlp_AD, hlp_AD, xlp_AD, vlp_AD )
    g_AD = ZERO

    ! Single scatter albedo
    z4 => AeroC%w(asi%i1:asi%i2, asi%h1:asi%h2, asi%j1:asi%j2, asi%k1:asi%k2, k)
    CALL Interp_4D_AD( z4, asi%wlp, asi%hlp, asi%xlp, asi%vlp, w_AD, &
                       z4_AD, wlp_AD, hlp_AD, xlp_AD, vlp_AD )
    w_AD = ZERO

    ! Extinction coefficient
    z4 => AeroC%ke(asi%i1:asi%i2, asi%h1:asi%h2, asi%j1:asi%j2, asi%k1:asi%k2, k)
    CALL Interp_4D_AD( z4, asi%wlp, asi%hlp, asi%xlp, asi%vlp, ke_AD, &
                       z4_AD, wlp_AD, hlp_AD, xlp_AD, vlp_AD )
    ke_AD = ZERO
    NULLIFY(z4)

    ! Compute the AD of the interpolating polynomials
    CALL LPoly_AD( asi%v, asi%v_int, asi%vlp, vlp_AD, v_AD, v_int_AD )
    CALL LPoly_AD( asi%r, asi%r_int, asi%xlp, xlp_AD, r_AD, r_int_AD )
    CALL LPoly_AD( asi%h, asi%h_int, asi%hlp, hlp_AD, h_AD, h_int_AD )
    CALL LPoly_AD( asi%f, asi%f_int, asi%wlp, wlp_AD, f_AD, f_int_AD )

    ! The AD outputs
    RH_AD   = RH_AD + h_int_AD
    Reff_AD = Reff_AD + r_int_AD
    Rsig_AD = Rsig_AD + v_int_AD

  END SUBROUTINE Get_Aerosol_Opt_AD

END MODULE CRTM_AerosolScatter
