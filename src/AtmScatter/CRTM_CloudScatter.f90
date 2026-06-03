!
! CRTM_CloudScatter
!
! Module to compute the cloud particle absorption and scattering properties
! required for radiative transfer in a cloudy atmosphere.
!
!
! CREATION HISTORY
!        Written by:     Quanhua Liu,    quanhua.liu@noaa.gov
!                        Yong Han,       yong.han@noaa.gov
!                        Paul van Delst, paul.vandelst@noaa.gov
!                        02-July-2005
!
!       Modified by:    James Rosinski, 08-Feb-2019
!                       Rosinski@ucar.edu
!
!       Modified by:    Isaac Moradi, 18-Sept-2021
!                       Isaac.Moradi@nasa.gov
!
! (C) Copyright 2019 UCAR
!

MODULE CRTM_CloudScatter

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE Type_Kinds,               ONLY: fp
  USE Message_Handler,          ONLY: SUCCESS, FAILURE, Display_Message
  USE CRTM_Parameters,          ONLY: ZERO, ONE, ONE_THOUSAND, &
                                      POINT_5, ONEpointFIVE, &
                                      MAX_N_LAYERS, &
                                      MAX_N_CLOUDS, &
                                      WATER_CONTENT_THRESHOLD, &
                                      BS_THRESHOLD, &
                                      MAX_N_LEGENDRE_TERMS, &
                                      MAX_N_PHASE_ELEMENTS, &
                                      HGPHASE  ! <<< NEED TO REMOVE THIS IN FUTURE
  USE CRTM_SpcCoeff,            ONLY: SC, &
                                      SpcCoeff_IsMicrowaveSensor , &
                                      SpcCoeff_IsInfraredSensor  , &
                                      SpcCoeff_IsVisibleSensor   , &
                                      SpcCoeff_IsUltravioletSensor
  USE CRTM_CloudCoeff,          ONLY: CloudC, &
                                      INVALID_CLOUDCOEFF, &
                                      MIE_TAMU_CLOUDCOEFF, &
                                      DDA_ARTS_CLOUDCOEFF, &
                                      CloudC_Exp, &
                                      Active_Cloud_Scheme, &
                                      CRTM_EXP_CLOUDCOEFF

  USE CRTM_Atmosphere_Define,   ONLY: CRTM_Atmosphere_type, &
                                      WATER_CLOUD, &
                                      ICE_CLOUD, &
                                      RAIN_CLOUD, &
                                      SNOW_CLOUD, &
                                      GRAUPEL_CLOUD, &
                                      HAIL_CLOUD, &
                                      PlateType1 , &
                                      ColumnType1 , &
                                      SixBulletRosette, &
                                      Perpendicular4_BulletRosette, &
                                      Flat3_BulletRosette, &
                                      IconCloudIce, &
                                      SectorSnowflake, &
                                      EvansSnowAggregate, &
                                      EightColumnAggregate, &
                                      LargePlateAggregate, &
                                      LargeColumnAggregate, &
                                      LargeBlockAggregate, &
                                      IconSnow, &
                                      IconHail, &
                                      GemGraupel, &
                                      GemSnow, &
                                      GemHail, &
                                      IceSphere, &
                                      LiquidSphere, &
                                      ! Parameters used for simplifying interpolation
                                      ! of cloud optical properties
                                      CLOUD_TYPE_MIE_TAMU, &
                                      CLOUD_INDEX_MIE_TAMU, &
                                      CLOUD_STATE_MIE_TAMU, &
                                      CLOUD_TYPE_DDA_ARTS, &
                                      CLOUD_INDEX_DDA_ARTS, &
                                      CLOUD_STATE_DDA_ARTS, &
                                      N_VALID_CLOUDS_MIE_TAMU, &
                                      N_VALID_CLOUDS_DDA_ARTS, &
                                      LIQUID, &
                                      FROZEN

  USE CRTM_GeometryInfo_Define, ONLY: CRTM_GeometryInfo_type
  USE CRTM_Interpolation,       ONLY: NPTS        , &
                                      LPoly_type  , &
                                      find_index  , &
                                      interp_1D   , &
                                      interp_2D   , &
                                      interp_3D   , &
                                      interp_2D_TL, &
                                      interp_3D_TL, &
                                      interp_2D_AD, &
                                      interp_3D_AD, &
                                      Clear_LPoly , &
                                      LPoly       , &
                                      LPoly_TL    , &
                                      LPoly_AD
  USE CRTM_AtmOptics_Define,    ONLY: CRTM_AtmOptics_type

  ! Internal variable definition module
  USE CSvar_Define, ONLY: CSvar_type, &
                          CSinterp_type, &
                          CSvar_Associated, &
                          CSvar_Destroy   , &
                          CSvar_Create

  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Procedures
  PUBLIC :: CRTM_Compute_CloudScatter
  PUBLIC :: CRTM_Compute_CloudScatter_TL
  PUBLIC :: CRTM_Compute_CloudScatter_AD


  ! -----------------
  ! Module parameters
  ! -----------------
  ! Version Id for the module
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: CRTM_CloudScatter.f90 99117 2017-11-27 18:37:14Z tong.zhu@noaa.gov $'
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
!       CRTM_Compute_CloudScatter
!
! PURPOSE:
!       Function to compute the cloud particle absorption and scattering
!       properties and populate the output CloudScatter structure for a
!       single channel.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_Compute_CloudScatter( Atmosphere  , &
!                                                 SensorIndex , &
!                                                 ChannelIndex, &
!                                                 CloudScatter, &
!                                                 CSvar         )
!
! INPUT ARGUMENTS:
!       Atmosphere:      CRTM_Atmosphere structure containing the atmospheric
!                        profile data.
!                        UNITS:      N/A
!                        TYPE:       CRTM_Atmosphere_type
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
! OUTPUT ARGUMENTS:
!        CloudScatter:   CRTM_AtmOptics structure containing the cloud particle
!                        absorption and scattering properties required for
!                        radiative transfer.
!                        UNITS:      N/A
!                        TYPE:       CRTM_AtmOptics_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
!        CSvar:          Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        UNITS:      N/A
!                        TYPE:       CSvar_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(OUT)
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
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_Compute_CloudScatter( &
    Atm         , &  ! Input
    GeometryInfo, &  ! Input
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    CScat       , &  ! Output
    CSV         ) &  ! Internal variable output
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_Atmosphere_type), INTENT(IN)     :: Atm
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)   :: GeometryInfo
    INTEGER                   , INTENT(IN)     :: SensorIndex
    INTEGER                   , INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_AtmOptics_type) , INTENT(IN OUT) :: CScat
    TYPE(CSvar_type)          , INTENT(IN OUT) :: CSV
    ! Function result
    INTEGER :: Error_Status
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Compute_CloudScatter'
    ! Local variables
    CHARACTER(ML) :: Message
    INTEGER  :: k, kc, l, m, n, j
    INTEGER  :: Le                       ! experimental-scheme effective truncation order
    REAL(fp) :: Frequency_MW, Frequency_IR
    LOGICAL  :: Layer_Mask(Atm%n_Layers)
    INTEGER  :: Layer_Index(Atm%n_Layers)
    INTEGER  :: nCloud_Layers
    REAL(fp) :: bs

    ! ------
    ! Set up
    ! ------
    Error_Status = SUCCESS

    IF (Atm%n_Clouds == 0) RETURN
    CSV%Total_bs = ZERO
    ! Spectral variables
    Frequency_MW = SC(SensorIndex)%Frequency(ChannelIndex)
    Frequency_IR = SC(SensorIndex)%Wavenumber(ChannelIndex)
    ! Determine the phase-function truncation order.
    IF ( Active_Cloud_Scheme == CRTM_EXP_CLOUDCOEFF ) THEN
      ! Experimental scheme: truncation order is taken from the LUT (running
      ! maximum accumulated in the cloud loop below), DECOUPLED from the RT
      ! stream count. The legacy {4,6,8,16} lOffset block packing is not used.
      CScat%n_Legendre_Terms = 0
      CScat%lOffset          = 0
    ELSE
    ! Determine offset for Legendre coefficients in
    ! the CloudC lookup table corresponding to the
    ! number of streams
    SELECT CASE(CScat%n_Legendre_Terms)
      CASE (FOUR_STREAMS)   ; CScat%lOffset = 0
      CASE (SIX_STREAMS)    ; CScat%lOffset = 5
      CASE (EIGHT_STREAMS)  ; CScat%lOffset = 12
      CASE (SIXTEEN_STREAMS); CScat%lOffset = 21
      CASE DEFAULT
        CScat%lOffset = 0  ! Is this correct?
        ! Use two-stream model or HG and RAYLEIGH Phase function
        IF( HGPHASE ) THEN
          CScat%n_Legendre_Terms = 0
        ELSE
          Error_Status = FAILURE
          WRITE(Message,'("The n_Legendre_Terms in CloudScatter, ",i0,", do not fit model")') &
                        CScat%n_Legendre_Terms
          CALL Display_Message( ROUTINE_NAME,Message,Error_Status )
          RETURN
        END IF
    END SELECT
    END IF

    ! ---------------------------------------------
    ! Loop over the different clouds in the profile
    ! ---------------------------------------------
    Cloud_loop: DO n = 1, Atm%n_Clouds
      ! Only process clouds with more
      ! than the threshold water amount
      Layer_Mask    = Atm%Cloud(n)%Water_Content > WATER_CONTENT_THRESHOLD
      nCloud_Layers = COUNT(Layer_Mask)
      IF ( nCloud_Layers == 0 ) CYCLE Cloud_loop

      ! ------------------------------------
      ! Loop over the current cloud's layers
      ! ------------------------------------
      Layer_Index(1:nCloud_Layers) = PACK((/(k, k=1,Atm%Cloud(n)%n_Layers)/), Layer_Mask)

      Cloud_Layer_loop: DO k = 1, nCloud_Layers
        kc = Layer_Index(k)

        ! Call sensor specific routines
        IF ( Active_Cloud_Scheme == CRTM_EXP_CLOUDCOEFF ) THEN
          ! Experimental 'CRTM-Exp' scheme (microwave only in v1)
          IF ( SpcCoeff_IsMicrowaveSensor(SC(SensorIndex)) ) THEN
            CALL Get_Cloud_Opt_MW_Exp(CScat                            , & ! Input
                                      Frequency_MW                     , & ! Input
                                      Atm%Cloud(n)%Type                , & ! Input
                                      Atm%Cloud(n)%Effective_Radius(kc), & ! Input (mapped to Dm)
                                      Atm%Temperature(kc)              , & ! Input
                                      CSV%ke(kc,n)                     , & ! Output
                                      CSV%kb(kc,n)                     , & ! Output
                                      CSV%w(kc,n)                      , & ! Output
                                      CSV%pcoeff(:,:,kc,n)             , & ! Output
                                      Le                               , & ! Output (effective truncation)
                                      CSV%csi(kc,n)                      ) ! Interpolation
            CScat%n_Legendre_Terms = MAX( CScat%n_Legendre_Terms, Le )
          ELSE
            CSV%ke(kc,n)         = ZERO
            CSV%kb(kc,n)         = ZERO
            CSV%w(kc,n)          = ZERO
            CSV%pcoeff(:,:,kc,n) = ZERO
          END IF
        ELSE IF ( SpcCoeff_IsMicrowaveSensor(SC(SensorIndex)) ) THEN
          CALL Get_Cloud_Opt_MW(CScat                            , & ! Input
                                Frequency_MW                     , & ! Input
                                Atm%Cloud(n)%Type                , & ! Input
                                Atm%Cloud(n)%Water_Density(kc)   , & ! Input
                                Atm%Cloud(n)%Effective_Radius(kc), & ! Input
                                Atm%Temperature(kc)              , & ! Input
                                CSV%ke(kc,n)                     , & ! Output
                                CSV%kb(kc,n)                     , & ! Output
                                CSV%w(kc,n)                      , & ! Output
                                CSV%pcoeff(:,:,kc,n)             , & ! Output
                                CSV%csi(kc,n)                      ) ! Interpolation
        ELSE IF ( SpcCoeff_IsInfraredSensor(SC(SensorIndex)) .OR. &
                  SpcCoeff_IsUltravioletSensor(SC(SensorIndex)) .OR. &
                  SpcCoeff_IsVisibleSensor( SC(SensorIndex))      ) THEN
          ! IR and visible use the same cloud optical data file, but distingished with Frequency
          CALL Get_Cloud_Opt_IR(CScat                            , & ! Input
                                Frequency_IR                     , & ! Input
                                Atm%Cloud(n)%Type                , & ! Input
                                Atm%Cloud(n)%Effective_Radius(kc), & ! Input
                                CSV%ke(kc,n)                     , & ! Output
                                CSV%kb(kc,n)                     , & ! Output
                                CSV%w(kc,n)                      , & ! Output
                                CSV%pcoeff(:,:,kc,n)             , & ! Output
                                CSV%csi(kc,n)                      ) ! Interpolation
        ELSE
          CSV%ke(kc,n)         = ZERO
          CSV%kb(kc,n)         = ZERO
          CSV%w(kc,n)          = ZERO
          CSV%pcoeff(:,:,kc,n) = ZERO
        END IF

        ! interpolation quality control
        IF( CSV%ke(kc,n) <= ZERO ) THEN
          CSV%ke(kc,n) = ZERO
          CSV%kb(kc,n) = ZERO
          CSV%w(kc,n)  = ZERO
        END IF
        IF( CSV%w(kc,n) <= ZERO ) THEN
          CSV%w(kc,n) = ZERO
          CSV%kb(kc,n) = ZERO
          CSV%pcoeff(:,:,kc,n) = ZERO
        END IF

        IF( CSV%w(kc,n) >= ONE ) THEN
          CSV%w(kc,n) = ONE
        END IF


        ! Compute the optical depth (absorption + scattering)
        !   tau = rho.ke
        ! where
        !   rho = integrated cloud water density for a layer (kg/m^2) [M.L^-2]
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
        CScat%Optical_Depth(kc) = CScat%Optical_Depth(kc) + &
                                  (CSV%ke(kc,n)*Atm%Cloud(n)%Water_Content(kc))
        IF ( SC(SensorIndex)%Is_Active_Sensor .and. CScat%Include_Scattering ) THEN
           CScat%Backscat_Coefficient(kc) = CScat%Backscat_Coefficient(kc) + &
                                     (CSV%kb(kc,n)*Atm%Cloud(n)%Water_Content(kc))
        END IF

        ! Compute the phase matrix coefficients
        !   p = p + p(LUT)*bs
        ! where
        !   p(LUT) = the phase coefficient from the LUT
        IF( CScat%n_Phase_Elements > 0  .and. CScat%Include_Scattering ) THEN
        ! Compute the volume scattering coefficient for the current
        ! cloud layer and accumulate it for the layer total for the
        ! profile (i.e. all clouds)
        !   bs = rho.w.ke
        ! where
        !   bs  = volume scattering coefficient for a layer [dimensionless]
        !   rho = integrated cloud water density for a layer (kg/m^2) [M.L^-2]
        !   w   = single scatter albedo [dimensionless]
        !   ke  = mass extintion coefficient (m^2/kg) [L^2.M^-1]
        bs = Atm%Cloud(n)%Water_Content(kc) * CSV%ke(kc,n) * CSV%w(kc,n)
        CSV%Total_bs(kc) = CSV%Total_bs(kc) + bs
        CScat%Single_Scatter_Albedo(kc) = CScat%Single_Scatter_Albedo(kc) + bs

          DO m = 1, CScat%n_Phase_Elements
            DO l = 1, CScat%n_Legendre_Terms
              CScat%Phase_Coefficient(l,m,kc) = CScat%Phase_Coefficient(l,m,kc) + &
                                                (CSV%pcoeff(l,m,kc,n) * bs)
            END DO
          END DO
        END IF
      END DO Cloud_Layer_loop
    END DO Cloud_loop

  END FUNCTION CRTM_Compute_CloudScatter


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Compute_CloudScatter_TL
!
! PURPOSE:
!       Function to compute the tangent-linear cloud particle absorption and
!       scattering properties and populate the output CloudScatter_TL structure
!       for a single channel.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_Compute_CloudScatter_TL( Atmosphere     , &
!                                                    CloudScatter   , &
!                                                    Atmosphere_TL  , &
!                                                    SensorIndex    , &
!                                                    ChannelIndex   , &
!                                                    CloudScatter_TL, &
!                                                    CSvar            )
!
! INPUT ARGUMENTS:
!       Atmosphere:       CRTM_Atmosphere structure containing the atmospheric
!                         profile data.
!                         UNITS:      N/A
!                         TYPE:       CRTM_Atmosphere_type
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN)
!
!       CloudScatter:     CRTM_AtmOptics structure containing the forward model
!                         cloud particle absorption and scattering properties
!                         required for radiative transfer.
!                         UNITS:      N/A
!                         TYPE:       CRTM_AtmOptics_type
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN)
!
!       Atmosphere_TL:    CRTM Atmosphere structure containing the tangent-linear
!                         atmospheric state data.
!                         UNITS:      N/A
!                         TYPE:       CRTM_Atmosphere_type
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN)
!
!       SensorIndex:      Sensor index id. This is a unique index associated
!                         with a (supported) sensor used to access the
!                         shared coefficient data for a particular sensor.
!                         See the ChannelIndex argument.
!                         UNITS:      N/A
!                         TYPE:       INTEGER
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:     Channel index id. This is a unique index associated
!                         with a (supported) sensor channel used to access the
!                         shared coefficient data for a particular sensor's
!                         channel.
!                         See the SensorIndex argument.
!                         UNITS:      N/A
!                         TYPE:       INTEGER
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN)
!
!       CSvar:            Structure containing internal variables required for
!                         subsequent tangent-linear or adjoint model calls.
!                         UNITS:      N/A
!                         TYPE:       CSvar_type
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN)
! OUTPUT ARGUMENTS:
!       CloudScatter_TL:  CRTM_AtmOptics structure containing the tangent-linear
!                         cloud particle absorption and scattering properties
!                         required for radiative transfer.
!                         UNITS:      N/A
!                         TYPE:       CRTM_AtmOptics_type
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN OUT)
!
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
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_Compute_CloudScatter_TL( &
    Atm         , &  ! FWD Input
    CScat       , &  ! FWD Input
    Atm_TL      , &  ! TL  Input
    GeometryInfo, &  ! Input
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    CScat_TL    , &  ! TL  Output
    CSV         ) &  ! Internal variable input
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_Atmosphere_type), INTENT(IN)     :: Atm
    TYPE(CRTM_AtmOptics_type) , INTENT(IN)     :: CScat
    TYPE(CRTM_Atmosphere_type), INTENT(IN)     :: Atm_TL
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)   :: GeometryInfo
    INTEGER                   , INTENT(IN)     :: SensorIndex
    INTEGER                   , INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_AtmOptics_type) , INTENT(IN OUT) :: CScat_TL
    TYPE(CSvar_type)          , INTENT(IN)     :: CSV
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Compute_CloudScatter'
    ! Local variables
    INTEGER  :: k, kc, l, m, n
    INTEGER  :: n_Legendre_Terms, n_Phase_Elements
    REAL(fp) :: Frequency_MW, Frequency_IR
    LOGICAL  :: Layer_Mask(Atm%n_Layers)
    INTEGER  :: Layer_Index(Atm%n_Layers)
    INTEGER  :: nCloud_Layers
    REAL(fp) :: ke_TL, w_TL, kb_TL
    REAL(fp) :: pcoeff_TL(0:CScat%n_Legendre_Terms, CScat%n_Phase_Elements)
    REAL(fp) :: bs, bs_TL

    ! ------
    ! Set up
    ! ------
    Error_Status = SUCCESS
    IF (Atm%n_Clouds == 0) RETURN
    ! Experimental scheme sets n_Legendre_Terms dynamically in the forward; mirror it
    ! onto AtmOptics_TL so the local n_Legendre_Terms (below) and the downstream
    ! RT/Combine/clear-sky-copy TL stay congruent with the forward.
    IF ( Active_Cloud_Scheme == CRTM_EXP_CLOUDCOEFF ) CScat_TL%n_Legendre_Terms = CScat%n_Legendre_Terms
    ! Spectral variables
    Frequency_MW = SC(SensorIndex)%Frequency(ChannelIndex)
    Frequency_IR = SC(SensorIndex)%Wavenumber(ChannelIndex)
    ! Phase matrix dimensions
    n_Legendre_Terms = CScat_TL%n_Legendre_Terms
    n_Phase_Elements = CScat_TL%n_Phase_Elements
    CScat_TL%lOffset = CScat%lOffset

    ! ---------------------------------------------
    ! Loop over the different clouds in the profile
    ! ---------------------------------------------
    Cloud_loop: DO n = 1, Atm%n_Clouds
      ! Only process clouds with more
      ! than the threshold water amount
      Layer_Mask    = Atm%Cloud(n)%Water_Content > WATER_CONTENT_THRESHOLD
      nCloud_Layers = COUNT(Layer_Mask)
      IF ( nCloud_Layers == 0 ) CYCLE Cloud_loop

      ! ------------------------------------
      ! Loop over the current cloud's layers
      ! ------------------------------------
      Layer_Index(1:nCloud_Layers) = PACK((/(k, k=1,Atm%Cloud(n)%n_Layers)/), Layer_Mask)
      Cloud_Layer_loop: DO k = 1, nCloud_Layers
        kc = Layer_Index(k)

        ! Call sensor specific routines
        IF ( Active_Cloud_Scheme == CRTM_EXP_CLOUDCOEFF ) THEN
          IF ( SpcCoeff_IsMicrowaveSensor(SC(SensorIndex)) ) THEN
            CALL Get_Cloud_Opt_MW_Exp_TL( CScat_TL, Atm%Cloud(n)%Type, CSV%ke(kc,n), CSV%w(kc,n), &
                                          Atm_TL%Cloud(n)%Effective_Radius(kc), Atm_TL%Temperature(kc), &
                                          ke_TL, kb_TL, w_TL, pcoeff_TL, CSV%csi(kc,n) )
          ELSE
            ke_TL = ZERO ; kb_TL = ZERO ; w_TL = ZERO ; pcoeff_TL = ZERO
          END IF
        ELSE IF ( SpcCoeff_IsMicrowaveSensor(SC(SensorIndex)) ) THEN
          CALL Get_Cloud_Opt_MW_TL(CScat_TL                            , & ! Input
                                   Atm%Cloud(n)%Type                   , & ! Input
                                   CSV%ke(kc,n)                        , & ! Input
                                   CSV%kb(kc,n)                        , & ! Input
                                   CSV%w(kc,n)                         , & ! Input
                                   SC(SensorIndex)%Is_Active_Sensor    , & ! Input
                                   Atm_TL%Cloud(n)%Water_Density(kc)   , & ! TL  Input
                                   Atm_TL%Cloud(n)%Effective_Radius(kc), & ! TL  Input
                                   Atm_TL%Temperature(kc)              , & ! TL  Input
                                   ke_TL                               , & ! TL  Output
                                   kb_TL                               , & ! TL  Output
                                   w_TL                                , & ! TL  Output
                                   pcoeff_TL                           , & ! TL  Output
                                   CSV%csi(kc,n)                         ) ! Interpolation
        ELSE IF ( SpcCoeff_IsInfraredSensor(SC(SensorIndex)) .OR. &
                  SpcCoeff_IsUltravioletSensor(SC(SensorIndex)) .OR. &
                  SpcCoeff_IsVisibleSensor( SC(SensorIndex))      ) THEN
          CALL Get_Cloud_Opt_IR_TL(CScat_TL                            , & ! Input
                                   Atm%Cloud(n)%Type                   , & ! Input
                                   CSV%ke(kc,n)                        , & ! Input
                                   CSV%kb(kc,n)                        , & ! Input
                                   CSV%w(kc,n)                         , & ! Input
                                   SC(SensorIndex)%Is_Active_Sensor    , & ! Input
                                   Atm_TL%Cloud(n)%Effective_Radius(kc), & ! TL  Input
                                   ke_TL                               , & ! TL  Output
                                   kb_TL                               , & ! TL  Output
                                   w_TL                                , & ! TL  Output
                                   pcoeff_TL                           , & ! TL  Output
                                   CSV%csi(kc,n)                         ) ! Interpolation
        ELSE
          ke_TL     = ZERO
          kb_TL     = ZERO
          w_TL      = ZERO
          pcoeff_TL = ZERO
        END IF

        ! interpolation quality control
        IF( CSV%ke(kc,n) <= ZERO ) THEN
          ke_TL = ZERO
          kb_TL = ZERO
          w_TL  = ZERO
        END IF
        IF( CSV%w(kc,n) <= ZERO ) THEN
          w_TL = ZERO
          kb_TL = ZERO
          pcoeff_TL = ZERO
        END IF
        IF( CSV%w(kc,n) >= ONE ) THEN
          w_TL = ZERO
        END IF

        ! Compute the optical depth (absorption + scattering)
        CScat_TL%Optical_Depth(kc) = CScat_TL%Optical_Depth(kc) + &
                                     (ke_TL        * Atm%Cloud(n)%Water_Content(kc)) + &
                                     (CSV%ke(kc,n) * Atm_TL%Cloud(n)%Water_Content(kc))

        ! Compute the Backscatering coefficient
        CScat_TL%Backscat_Coefficient(kc) = CScat_TL%Backscat_Coefficient(kc) + &
                                     (kb_TL        * Atm%Cloud(n)%Water_Content(kc)) + &
                                     (CSV%kb(kc,n) * Atm_TL%Cloud(n)%Water_Content(kc))

        ! Compute the phase matrix coefficients
        IF( n_Phase_Elements > 0 .and. CScat%Include_Scattering ) THEN
        ! Compute the volume scattering coefficient
        bs = Atm%Cloud(n)%Water_Content(kc) * CSV%ke(kc,n) * CSV%w(kc,n)
        bs_TL = (Atm_TL%Cloud(n)%Water_Content(kc) * CSV%ke(kc,n) * CSV%w(kc,n) ) + &
                (Atm%Cloud(n)%Water_Content(kc) * ke_TL * CSV%w(kc,n) ) + &
                (Atm%Cloud(n)%Water_Content(kc) * CSV%ke(kc,n) * w_TL )

        CScat_TL%Single_Scatter_Albedo(kc) = CScat_TL%Single_Scatter_Albedo(kc) + bs_TL

          DO m = 1, n_Phase_Elements
            DO l = 0, n_Legendre_Terms
              CScat_TL%Phase_Coefficient(l,m,kc) = CScat_TL%Phase_Coefficient(l,m,kc) + &
                                                   (pcoeff_TL(l,m)       * bs   ) + &
                                                   (CSV%pcoeff(l,m,kc,n) * bs_TL)
            END DO
          END DO
        END IF
      END DO Cloud_Layer_loop
    END DO Cloud_loop


  END FUNCTION CRTM_Compute_CloudScatter_TL


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Compute_CloudScatter_AD
!
! PURPOSE:
!       Function to compute the adjoint of the cloud particle absorption and
!       scattering properties for a single channel.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_Compute_CloudScatter_AD(  Atmosphere     , &
!                                                     CloudScatter   , &
!                                                     CloudScatter_AD, &
!                                                     SensorIndex    , &
!                                                     ChannelIndex   , &
!                                                     Atmosphere_AD  , &
!                                                     CSvar            )
!
! INPUT ARGUMENTS:
!       Atmosphere:       CRTM_Atmosphere structure containing the atmospheric
!                         profile data.
!                         UNITS:      N/A
!                         TYPE:       CRTM_Atmosphere_type
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN)
!
!       CloudScatter:     CRTM_AtmOptics structure containing the forward model
!                         cloud particle absorption and scattering properties
!                         required for radiative transfer.
!                         UNITS:      N/A
!                         TYPE:       CRTM_AtmOptics_type
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN)
!
!       CloudScatter_AD:  CRTM_AtmOptics structure containing the adjoint
!                         of the cloud particle absorption and scattering
!                         properties required for radiative transfer.
!                         **NOTE: On EXIT from this function, the contents of
!                                 this structure may be modified (e.g. set to
!                                 zero.)
!                         UNITS:      N/A
!                         TYPE:       CRTM_AtmOptics_type
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN OUT)
!
!       SensorIndex:      Sensor index id. This is a unique index associated
!                         with a (supported) sensor used to access the
!                         shared coefficient data for a particular sensor.
!                         See the ChannelIndex argument.
!                         UNITS:      N/A
!                         TYPE:       INTEGER
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:     Channel index id. This is a unique index associated
!                         with a (supported) sensor channel used to access the
!                         shared coefficient data for a particular sensor's
!                         channel.
!                         See the SensorIndex argument.
!                         UNITS:      N/A
!                         TYPE:       INTEGER
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN)
!
!       CSvar:            Structure containing internal variables required for
!                         subsequent tangent-linear or adjoint model calls.
!                         UNITS:      N/A
!                         TYPE:       CSvar_type
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       Atmosphere_AD:    CRTM Atmosphere structure containing the adjoint
!                         atmospheric state data.
!                         UNITS:      N/A
!                         TYPE:       CRTM_Atmosphere_type
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN OUT)
!
!
! FUNCTION RESULT:
!       Error_Status:     The return value is an integer defining the error status.
!                         The error codes are defined in the Message_Handler module.
!                         If == SUCCESS the computation was sucessful
!                            == FAILURE an unrecoverable error occurred
!                         UNITS:      N/A
!                         TYPE:       INTEGER
!                         DIMENSION:  Scalar
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_Compute_CloudScatter_AD( &
    Atm         , &  ! FWD Input
    CScat       , &  ! FWD Input
    CScat_AD    , &  ! AD  Input
    GeometryInfo, &  ! Input
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    Atm_AD      , &  ! AD  Output
    CSV         ) &  ! Internal variable input
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_Atmosphere_type), INTENT(IN)     :: Atm
    TYPE(CRTM_AtmOptics_type) , INTENT(IN)     :: CScat
    TYPE(CRTM_AtmOptics_type) , INTENT(IN OUT) :: CScat_AD
    TYPE(CRTM_GeometryInfo_type), INTENT(IN) :: GeometryInfo
    INTEGER                   , INTENT(IN)     :: SensorIndex
    INTEGER                   , INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_Atmosphere_type), INTENT(IN OUT) :: Atm_AD
    TYPE(CSvar_type)          , INTENT(IN)     :: CSV
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Compute_CloudScatter_AD'
    ! Local variables
    INTEGER  :: k, kc, l, m, n
    INTEGER  :: n_Legendre_Terms, n_Phase_Elements
    REAL(fp) :: Frequency_MW, Frequency_IR
    LOGICAL  :: Layer_Mask(Atm%n_Layers)
    INTEGER  :: Layer_Index(Atm%n_Layers)
    INTEGER  :: nCloud_Layers
    REAL(fp) :: ke_AD, w_AD, kb_AD
    REAL(fp) :: pcoeff_AD(0:CScat%n_Legendre_Terms, CScat%n_Phase_Elements)
    REAL(fp) :: bs, bs_AD

    ! ------
    ! Set up
    ! ------
    Error_Status = SUCCESS
    IF ( Atm%n_Clouds == 0 ) RETURN
    ! Experimental scheme: mirror the dynamic forward Legendre count onto AtmOptics_AD
    ! so the local n_Legendre_Terms (below) matches the forward.
    IF ( Active_Cloud_Scheme == CRTM_EXP_CLOUDCOEFF ) CScat_AD%n_Legendre_Terms = CScat%n_Legendre_Terms
    ! Spectral variables
    Frequency_MW = SC(SensorIndex)%Frequency(ChannelIndex)
    Frequency_IR = SC(SensorIndex)%Wavenumber(ChannelIndex)
    ! Phase matrix dimensions
    n_Legendre_Terms = CScat_AD%n_Legendre_Terms
    n_Phase_Elements = CScat_AD%n_Phase_Elements
    CScat_AD%lOffset = CScat%lOffset

    ! ---------------------------------------------
    ! Loop over the different clouds in the profile
    ! ---------------------------------------------
    Cloud_loop: DO n = 1, Atm%n_Clouds
      ! Only process clouds with more
      ! than the threshold water amount
      Layer_Mask    = Atm%Cloud(n)%Water_Content > WATER_CONTENT_THRESHOLD
      nCloud_Layers = COUNT(Layer_Mask)
      IF ( nCloud_Layers == 0 ) CYCLE Cloud_loop

      ! ------------------------------------
      ! Loop over the current cloud's layers
      ! ------------------------------------
      Layer_Index(1:nCloud_Layers) = PACK((/(k, k=1,Atm%Cloud(n)%n_Layers)/), Layer_Mask)
      Cloud_Layer_loop: DO k = nCloud_Layers, 1, -1 !1, nCloud_Layers
        kc = Layer_Index(k)

        ! Initialise the individual
        ! cloud adjoint variables
        bs_AD = ZERO
        pcoeff_AD = ZERO
        ke_AD     = ZERO
        kb_AD     = ZERO
        w_AD      = ZERO

        ! Compute the adjoint of the
        ! phase matrix coefficients
        IF( n_Phase_Elements > 0 .and. CScat%Include_Scattering ) THEN
        ! Recompute the forward model volume scattering
        ! coefficient for the current cloud ONLY
        bs = Atm%Cloud(n)%Water_Content(kc) * CSV%ke(kc,n) * CSV%w(kc,n)

          DO m = 1, n_Phase_Elements
            DO l = 0, n_Legendre_Terms
              bs_AD = bs_AD + (CSV%pcoeff(l,m,kc,n) * CScat_AD%Phase_Coefficient(l,m,kc))
              pcoeff_AD(l,m) = pcoeff_AD(l,m) + (bs * CScat_AD%Phase_Coefficient(l,m,kc))
            END DO
          END DO
        ! NOTE: bs_AD is not reinitialised after this
        !       point since it is reinitialised at the
        !       start of the Cloud_Layer_loop.
        bs_AD = bs_AD + CScat_AD%Single_Scatter_Albedo(kc)
        w_AD  = w_AD  + (Atm%Cloud(n)%Water_Content(kc) * CSV%ke(kc,n)* bs_AD )
        END IF

        ! Compute the adjoint of the optical
        ! depth (absorption + scattering)
        Atm_AD%Cloud(n)%Water_Content(kc) = Atm_AD%Cloud(n)%Water_Content(kc) + &
                                            (CSV%ke(kc,n) * CScat_AD%Optical_Depth(kc))
        ke_AD = ke_AD + (Atm%Cloud(n)%Water_Content(kc) * CScat_AD%Optical_Depth(kc))

        IF ( SC(SensorIndex)%Is_Active_Sensor .and. CScat%Include_Scattering ) THEN
           Atm_AD%Cloud(n)%Water_Content(kc) = Atm_AD%Cloud(n)%Water_Content(kc) + &
                                               (CSV%kb(kc,n) * CScat_AD%Backscat_Coefficient(kc))
           kb_AD = kb_AD + (Atm%Cloud(n)%Water_Content(kc) * CScat_AD%Backscat_Coefficient(kc))
         END IF

        ! Compute the adjoint of the volume
        ! scattering coefficient.
        ke_AD = ke_AD + (Atm%Cloud(n)%Water_Content(kc) * bs_AD * CSV%w(kc,n) )
        Atm_AD%Cloud(n)%Water_Content(kc) = Atm_AD%Cloud(n)%Water_Content(kc) + &
                                              ( bs_AD * CSV%ke(kc,n) * CSV%w(kc,n) )

        ! interpolation quality control
        IF( CSV%w(kc,n) >= ONE ) THEN
          w_AD = ZERO
        END IF
        IF( CSV%ke(kc,n) <= ZERO ) THEN
          ke_AD = ZERO
          kb_AD = ZERO
          w_AD  = ZERO
        END IF

        IF( CSV%kb(kc,n) <= ZERO ) THEN
          kb_AD = ZERO
        END IF

        IF( CSV%w(kc,n) <= ZERO ) THEN
          w_AD = ZERO
          kb_AD = ZERO
          pcoeff_AD = ZERO
        END IF

        ! Call sensor specific routines
        IF ( Active_Cloud_Scheme == CRTM_EXP_CLOUDCOEFF ) THEN
          IF ( SpcCoeff_IsMicrowaveSensor(SC(SensorIndex)) ) THEN
            CALL Get_Cloud_Opt_MW_Exp_AD( CScat_AD, Atm%Cloud(n)%Type, CSV%ke(kc,n), CSV%w(kc,n), &
                                          ke_AD, kb_AD, w_AD, pcoeff_AD, &
                                          Atm_AD%Cloud(n)%Effective_Radius(kc), Atm_AD%Temperature(kc), &
                                          CSV%csi(kc,n) )
          ELSE
            ke_AD = ZERO ; kb_AD = ZERO ; w_AD = ZERO ; pcoeff_AD = ZERO
          END IF
        ELSE IF ( SpcCoeff_IsMicrowaveSensor(SC(SensorIndex)) ) THEN
          CALL Get_Cloud_Opt_MW_AD(CScat_AD                            , & ! Input
                                   Atm%Cloud(n)%Type                   , & ! Input
                                   CSV%ke(kc,n)                        , & ! Input
                                   CSV%kb(kc,n)                        , & ! Input
                                   CSV%w(kc,n)                         , & ! Input
                                   SC(SensorIndex)%Is_Active_Sensor    , & ! Input
                                   ke_AD                               , & ! AD  Input
                                   kb_AD                               , & ! AD  Input
                                   w_AD                                , & ! AD  Input
                                   pcoeff_AD                           , & ! AD  Input
                                   Atm_AD%Cloud(n)%Water_Density(kc)   , & ! AD  Output
                                   Atm_AD%Cloud(n)%Effective_Radius(kc), & ! AD  Output
                                   Atm_AD%Temperature(kc)              , & ! AD  Output
                                   CSV%csi(kc,n)                         ) ! Interpolation
        ELSE IF ( SpcCoeff_IsInfraredSensor(SC(SensorIndex)) .OR. &
                  SpcCoeff_IsUltravioletSensor(SC(SensorIndex)) .OR. &
                  SpcCoeff_IsVisibleSensor( SC(SensorIndex))      ) THEN
          CALL Get_Cloud_Opt_IR_AD(CScat_AD                            , & ! Input
                                   Atm%Cloud(n)%Type                   , & ! Input
                                   CSV%ke(kc,n)                        , & ! Input
                                   CSV%kb(kc,n)                        , & ! Input
                                   CSV%w(kc,n)                         , & ! Input
                                   SC(SensorIndex)%Is_Active_Sensor    , & ! Input
                                   ke_AD                               , & ! AD  Input
                                   kb_AD                               , & ! AD  Input
                                   w_AD                                , & ! AD  Input
                                   pcoeff_AD                           , & ! AD  Input
                                   Atm_AD%Cloud(n)%Effective_Radius(kc), & ! AD  Output
                                   CSV%csi(kc,n)                         ) ! Interpolation
        ELSE
          ke_AD     = ZERO
          kb_AD     = ZERO
          w_AD      = ZERO
          pcoeff_AD = ZERO
        END IF
      END DO Cloud_Layer_loop
    END DO Cloud_loop

  END FUNCTION CRTM_Compute_CloudScatter_AD



!################################################################################
!################################################################################
!##                                                                            ##
!##                        ## PRIVATE MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

  ! ------------------------------------------
  ! Subroutine to obtain the IR bulk
  ! optical properties of a cloud:
  !   extinction coefficient (ke),
  !   scattering coefficient (w)
  !   asymmetry factor (g), and
  !   spherical Legendre coefficients (pcoeff)
  ! ------------------------------------------
  SUBROUTINE Get_Cloud_Opt_IR( CloudScatter, &  ! Input  CloudScatter structure
                               Frequency   , &  ! Input  Frequency (cm^-1)
                               cloud_type  , &  ! Input  see CRTM_Cloud_Define.f90
                               Reff        , &  ! Input  effective radius (mm)
                               ke          , &  ! Output optical depth for 1 mm water content
                               kb          , &  ! Output optical depth for 1 mm water content
                               w           , &  ! Output single scattering albedo
                               pcoeff      , &  ! Output spherical Legendre coefficients
                               csi           )  ! Output interpolation data
    ! Arguments
    TYPE(CRTM_AtmOptics_type), INTENT(IN)     :: CloudScatter
    REAL(fp)                  , INTENT(IN)     :: Frequency
    INTEGER                   , INTENT(IN)     :: Cloud_Type
    REAL(fp)                  , INTENT(IN)     :: Reff
    REAL(fp)                  , INTENT(OUT)    :: ke, kb
    REAL(fp)                  , INTENT(OUT)    :: w
    REAL(fp)                  , INTENT(IN OUT) :: pcoeff(0:,:)
    TYPE(CSinterp_type)       , INTENT(IN OUT) :: csi
    ! Local variables
    INTEGER  :: k, l, m


    ! Find the frequency and effective
    ! radius indices for interpolation
    ! --------------------------------
    csi%f_int = MAX(MIN(CloudC%Frequency_IR(CloudC%n_IR_Frequencies),Frequency),CloudC%Frequency_IR(1))
    CALL find_index(CloudC%Frequency_IR, csi%f_int, csi%i1,csi%i2, csi%f_outbound)
    csi%f = CloudC%Frequency_IR(csi%i1:csi%i2)
    csi%r_int = MAX(MIN(CloudC%Reff_IR(CloudC%n_IR_Radii),Reff),CloudC%Reff_IR(1))
    CALL find_index(CloudC%Reff_IR, csi%r_int, csi%j1,csi%j2, csi%r_outbound)
    csi%r = CloudC%Reff_IR(csi%j1:csi%j2)

    ! Calculate the interpolating polynomials
    ! ---------------------------------------
    ! Frequency term
    CALL LPoly( csi%f, csi%f_int, &  ! Input
                csi%wlp           )  ! Output
    ! Effective radius term
    CALL LPoly( csi%r, csi%r_int, &  ! Input
                csi%xlp           )  ! Output

    ! Determine the density index, k, for the clouds
    ! based on CloudC LUT organisation
    ! ----------------------------------------------
    SELECT CASE (Cloud_Type)
      CASE(WATER_CLOUD)  ; k=1  ! Liquid
      CASE(ICE_CLOUD)    ; k=4  ! Solid
      CASE(RAIN_CLOUD)   ; k=1  ! Liquid
      CASE(SNOW_CLOUD)   ; k=2  ! Solid
      CASE(GRAUPEL_CLOUD); k=3  ! Solid
      CASE(HAIL_CLOUD)   ; k=4  ! Solid
      CASE DEFAULT       ; k=HUGE(k)   ! Case not found: Hopefully induce code to fail
    END SELECT

    ! Perform interpolation
    ! ---------------------
    CALL interp_2D( CloudC%ke_IR(csi%i1:csi%i2,csi%j1:csi%j2,k), csi%wlp, csi%xlp, ke )
    CALL interp_2D( CloudC%kb_IR(csi%i1:csi%i2,csi%j1:csi%j2,k), csi%wlp, csi%xlp, kb )
    CALL interp_2D( CloudC%w_IR(csi%i1:csi%i2,csi%j1:csi%j2,k) , csi%wlp, csi%xlp, w  )
    IF (CloudScatter%n_Phase_Elements > 0 .and. CloudScatter%Include_Scattering ) THEN
      pcoeff(0,1) = POINT_5
      DO m = 1, CloudScatter%n_Phase_Elements
      DO l = 1, CloudScatter%n_Legendre_Terms
        CALL interp_2D( CloudC%pcoeff_IR(csi%i1:csi%i2,csi%j1:csi%j2,k,l+CloudScatter%lOffset,m), &
                        csi%wlp, csi%xlp, pcoeff(l,m) )
      END DO
      END DO
    ELSE
      ! Absorption coefficient
      ke = ke *(ONE - w)
    END IF

  END SUBROUTINE Get_Cloud_Opt_IR


  ! ---------------------------------------------
  ! Subroutine to obtain the tangent-linear
  ! IR bulk optical properties of a cloud:
  !   extinction coefficient (ke_TL),
  !   scattereing coefficient (w_TL)
  !   spherical Legendre coefficients (pcoeff_TL)
  ! ---------------------------------------------
  SUBROUTINE Get_Cloud_Opt_IR_TL( CloudScatter_TL , &  ! Input      CloudScatter TL structure
                                  cloud_type      , &  ! Input      see CRTM_Cloud_Define.f90
                                  ke              , &  ! Input
                                  kb              , &  ! Input
                                  w               , &  ! Input
                                  Is_Active_Sensor, & ! Input
                                  Reff_TL         , &  ! TL  Input  effective radius (mm)
                                  ke_TL           , &  ! TL  Output extinction coefficient (=~ optical depth for 1 mm water content)
                                  kb_TL           , &  ! TL  Output backscat coefficient
                                  w_TL            , &  ! TL  Output single scattering albedo
                                  pcoeff_TL       , &  ! TL  Output spherical Legendre coefficients
                                  csi              )  ! Input interpolation data
    ! Arguments
    TYPE(CRTM_AtmOptics_type), INTENT(IN)     :: CloudScatter_TL
    INTEGER ,                   INTENT(IN)     :: Cloud_Type
    REAL(fp),                   INTENT(IN)     :: ke, kb, w, Reff_TL
    LOGICAL,                    INTENT(IN)     :: Is_Active_Sensor
    REAL(fp),                   INTENT(OUT)    :: ke_TL, kb_TL
    REAL(fp),                   INTENT(OUT)    :: w_TL
    REAL(fp),                   INTENT(IN OUT) :: pcoeff_TL(0:,:)
    TYPE(CSinterp_type),        INTENT(IN)     :: csi
    ! Local variables
    INTEGER  :: k, l, m
    REAL(fp) :: f_int_TL, r_int_TL
    REAL(fp) :: f_TL(NPTS), r_TL(NPTS)
    REAL(fp) :: z_TL(NPTS,NPTS)
    TYPE(LPoly_type) :: wlp_TL, xlp_TL
! JR    REAL(fp), POINTER :: z(:,:) => NULL()
    REAL(fp), POINTER :: z(:,:)

    NULLIFY(z)
    ! Setup
    ! -----
    ! No TL output when all dimensions
    ! are outside LUT bounds
    IF ( csi%f_outbound .AND. csi%r_outbound ) THEN
      ke_TL     = ZERO
      kb_TL     = ZERO
      w_TL      = ZERO
      pcoeff_TL = ZERO
      RETURN
    END IF
    ! The TL inputs
    f_int_TL = ZERO
    f_TL     = ZERO
    r_int_TL = Reff_TL
    r_TL     = ZERO
    z_TL = ZERO


    ! Calculate the TL interpolating polynomials
    ! ------------------------------------------
    ! Frequency term (always zero. This is a placeholder for testing)
    CALL LPoly_TL( csi%f, csi%f_int, & ! FWD Input
                   csi%wlp,          & ! FWD Input
                   f_TL, f_int_TL,   & ! TL  Input
                   wlp_TL            ) ! TL  Output
    ! Effective radius term
    CALL LPoly_TL( csi%r, csi%r_int, & ! FWD Input
                   csi%xlp,          & ! FWD Input
                   r_TL, r_int_TL,   & ! TL  Input
                   xlp_TL            ) ! TL  Output


    ! Determine the density index, k, for the clouds
    ! based on CloudC LUT organisation
    ! ----------------------------------------------
    SELECT CASE (Cloud_Type)
      CASE(WATER_CLOUD)  ; k=1  ! Liquid
      CASE(ICE_CLOUD)    ; k=4  ! Solid
      CASE(RAIN_CLOUD)   ; k=1  ! Liquid
      CASE(SNOW_CLOUD)   ; k=2  ! Solid
      CASE(GRAUPEL_CLOUD); k=3  ! Solid
      CASE(HAIL_CLOUD)   ; k=4  ! Solid
      CASE DEFAULT       ; k=HUGE(k)   ! Case not found: Hopefully induce code to fail
    END SELECT


    ! Perform interpolation based on cloud type
    ! -----------------------------------------
    ! Extinction coefficient
    z => CloudC%ke_IR(csi%i1:csi%i2,csi%j1:csi%j2,k)
    CALL interp_2D_TL( z   , csi%wlp, csi%xlp, &  ! FWD Input
                       z_TL, wlp_TL , xlp_TL , &  ! TL  Input
                       ke_TL                   )  ! TL  Output
    ! Extinction coefficient
    z => CloudC%kb_IR(csi%i1:csi%i2,csi%j1:csi%j2,k)
    CALL interp_2D_TL( z   , csi%wlp, csi%xlp, &  ! FWD Input
                       z_TL, wlp_TL , xlp_TL , &  ! TL  Input
                       kb_TL                   )  ! TL  Output
    ! Single scatter albedo
    z => CloudC%w_IR(csi%i1:csi%i2,csi%j1:csi%j2,k)
    CALL interp_2D_TL( z   , csi%wlp, csi%xlp, &  ! FWD Input
                       z_TL, wlp_TL , xlp_TL , &  ! TL  Input
                       w_TL                    )  ! TL  Output
    ! Phase matrix coefficients
    IF ( CloudScatter_TL%n_Phase_Elements > 0 .and. CloudScatter_TL%Include_Scattering ) THEN
      pcoeff_TL(0,1) = ZERO

      DO m = 1, CloudScatter_TL%n_Phase_Elements
      DO l = 1, CloudScatter_TL%n_Legendre_Terms
        z => CloudC%pcoeff_IR(csi%i1:csi%i2,csi%j1:csi%j2,k,l+CloudScatter_TL%lOffset,m)
        CALL interp_2D_TL( z   , csi%wlp, csi%xlp, &  ! FWD Input
                           z_TL, wlp_TL , xlp_TL , &  ! TL  Input
                           pcoeff_TL(l,m)          )  ! TL  Output
      END DO
      END DO
    ELSE
      ! Absorption coefficient
      IF( w < ONE ) THEN
        ke_TL = ke_TL * (ONE - w) - ke/(ONE -w) * w_TL
      ELSE
        ke_TL = ZERO
      END IF
    END IF
    NULLIFY(z)

  END SUBROUTINE Get_Cloud_Opt_IR_TL


  ! ---------------------------------------------
  ! Subroutine to obtain the adjoint of the
  ! IR bulk optical properties of a cloud:
  !   effective radius (Reff_AD)
  ! ---------------------------------------------
  SUBROUTINE Get_Cloud_Opt_IR_AD( CloudScatter_AD , &  ! Input      CloudScatter AD structure
                                  cloud_type      , &  ! Input      see CRTM_Cloud_Define.f90
                                  ke              , &  ! Input
                                  kb              , &  ! Input
                                  w               , &  ! Input
                                  Is_Active_Sensor, & ! Input
                                  ke_AD           , &  ! AD  Input  extinction coefficient (=~ optical depth for 1 mm water content)
                                  kb_AD           , &  ! AD  Input  Backscatering coefficient
                                  w_AD            , &  ! AD  Input  single scattering albedo
                                  pcoeff_AD       , &  ! AD  Input  spherical Legendre coefficients
                                  Reff_AD         , &  ! AD  Output effective radius (mm)
                                  csi              )  ! Input interpolation data
    ! Arguments
    TYPE(CRTM_AtmOptics_type), INTENT(IN)     :: CloudScatter_AD
    INTEGER ,                   INTENT(IN)     :: Cloud_Type
    REAL(fp),                   INTENT(IN)     :: ke, kb, w
    LOGICAL,                    INTENT(IN)     :: Is_Active_Sensor
    REAL(fp),                   INTENT(IN OUT) :: ke_AD, kb_AD    ! AD  Input
    REAL(fp),                   INTENT(IN OUT) :: w_AD            ! AD  Input
    REAL(fp),                   INTENT(IN OUT) :: pcoeff_AD(0:,:) ! AD  Input
    REAL(fp),                   INTENT(IN OUT) :: Reff_AD         ! AD  Output
    TYPE(CSinterp_type),        INTENT(IN)     :: csi
    ! Local variables
    INTEGER  :: k, l, m
    REAL(fp) :: f_int_AD, r_int_AD
    REAL(fp) :: f_AD(NPTS), r_AD(NPTS)
    REAL(fp) :: z_AD(NPTS,NPTS)
    TYPE(LPoly_type) :: wlp_AD, xlp_AD
!    REAL(fp), POINTER :: z(:,:) => NULL()
    REAL(fp), POINTER :: z(:,:)

    NULLIFY(z)
    ! Setup
    ! -----
    ! No AD output when all dimensions
    ! are outside LUT bounds
    IF ( csi%f_outbound .AND. csi%r_outbound ) THEN
      Reff_AD     = ZERO
      ke_AD     = ZERO
      kb_AD     = ZERO
      w_AD      = ZERO
      pcoeff_AD = ZERO
      RETURN
    END IF
    ! Initialise local adjoint variables
    f_int_AD = ZERO
    r_int_AD = ZERO
    f_AD = ZERO
    r_AD = ZERO
    z_AD = ZERO
    CALL Clear_LPoly(wlp_AD)
    CALL Clear_LPoly(xlp_AD)


    ! Determine the density index, k, for the clouds
    ! based on CloudC LUT organisation
    ! ----------------------------------------------
    SELECT CASE (Cloud_Type)
      CASE(WATER_CLOUD)  ; k=1  ! Liquid
      CASE(ICE_CLOUD)    ; k=4  ! Solid
      CASE(RAIN_CLOUD)   ; k=1  ! Liquid
      CASE(SNOW_CLOUD)   ; k=2  ! Solid
      CASE(GRAUPEL_CLOUD); k=3  ! Solid
      CASE(HAIL_CLOUD)   ; k=4  ! Solid
      CASE DEFAULT       ; k=HUGE(k)   ! Case not found: Hopefully induce code to fail
    END SELECT

    ! Perform interpolation
    ! ---------------------
    ! Phase matrix coefficients
    IF (CloudScatter_AD%n_Phase_Elements > 0 .and. CloudScatter_AD%Include_Scattering ) THEN
      DO m = 1, CloudScatter_AD%n_Phase_Elements
      DO l = 1, CloudScatter_AD%n_Legendre_Terms
        z => CloudC%pcoeff_IR(csi%i1:csi%i2,csi%j1:csi%j2,k,l+CloudScatter_AD%lOffset,m)
        CALL interp_2D_AD( z   , csi%wlp   , csi%xlp   , &  ! FWD Input
                           pcoeff_AD(l,m)      , &  ! AD  Input
                           z_AD, wlp_AD, xlp_AD  )  ! AD  Output
      END DO
      END DO
      pcoeff_AD(0,1) = ZERO
    ELSE
      ! Absorption coefficient
      IF( w < ONE ) THEN
        w_AD  = w_AD - ke/(ONE -w) * ke_AD
        ke_AD = ke_AD * (ONE - w)
      ELSE
        ke_AD = ZERO
      END IF
    END IF
    ! Single scatter albedo
    z => CloudC%w_IR(csi%i1:csi%i2,csi%j1:csi%j2,k)
    CALL interp_2D_AD( z   , csi%wlp, csi%xlp, &  ! FWD Input
                       w_AD                  , &  ! AD  Input
                       z_AD, wlp_AD, xlp_AD    )  ! AD  Output
    ! Extinction coefficient
    z => CloudC%ke_IR(csi%i1:csi%i2,csi%j1:csi%j2,k)
    CALL interp_2D_AD( z   , csi%wlp, csi%xlp, &  ! FWD Input
                       ke_AD                 , &  ! AD  Input
                       z_AD, wlp_AD, xlp_AD    )  ! AD  Output
    ! Backscatering coefficient
    z => CloudC%kb_IR(csi%i1:csi%i2,csi%j1:csi%j2,k)
    CALL interp_2D_AD( z   , csi%wlp, csi%xlp, &  ! FWD Input
                       kb_AD                 , &  ! AD  Input
                       z_AD, wlp_AD, xlp_AD    )  ! AD  Output
    NULLIFY(z)


    ! Compute the AD of the interpolating polynomials
    ! -----------------------------------------------
    ! Efective radius term
    CALL LPoly_AD( csi%r, csi%r_int, & ! FWD Input
                   csi%xlp,          & ! FWD Input
                   xlp_AD,           & ! AD  Input
                   r_AD, r_int_AD    ) ! AD  Output
    ! Frequency term (always zero. This is a placeholder for testing)
    CALL LPoly_AD( csi%f, csi%f_int, & ! FWD Input
                   csi%wlp,          & ! FWD Input
                   wlp_AD,           & ! AD  Input
                   f_AD, f_int_AD    ) ! AD  Output

    ! The AD outputs
    ! --------------
    Reff_AD = Reff_AD + r_int_AD

  END SUBROUTINE Get_Cloud_Opt_IR_AD


  ! ------------------------------------------
  ! Subroutine to obtain the MW bulk
  ! optical properties of a cloud:
  !   extinction coefficient (ke),
  !   scattereing coefficient (w)
  !   asymmetry factor (g), and
  !   spherical Legendre coefficients (pcoeff)
  ! ------------------------------------------
  SUBROUTINE Get_Cloud_Opt_MW( CloudScatter  , &  ! Input  CloudScatter structure
                               Frequency     , &  ! Input  Frequency (GHz)
                               cloud_type    , &  ! Input  see CRTM_Cloud_Define.f90
                               Water_Density , &  ! Input  water density (kg.m^-3)
                               Reff          , &  ! Input  effective radius (mm)
                               Temperature   , &  ! Input  cloudy temperature
                               ke            , &  ! Input optical depth for 1 mm water content
                               kb            , &  ! Input Backscatering coeff for 1 mm water content
                               w             , &  ! Input single scattering albedo
                               pcoeff        , &  ! Output spherical Legendre coefficients
                               csi             )  ! Output interpolation data
    ! Arguments
    TYPE(CRTM_AtmOptics_type) , INTENT(IN)     :: CloudScatter
    REAL(fp)                  , INTENT(IN)     :: Frequency
    INTEGER                   , INTENT(IN)     :: Cloud_Type
    REAL(fp)                  , INTENT(IN)     :: Water_Density
    REAL(fp)                  , INTENT(IN)     :: Reff
    REAL(fp)                  , INTENT(IN)     :: Temperature
    REAL(fp)                  , INTENT(OUT)    :: ke, kb
    REAL(fp)                  , INTENT(OUT)    :: w
    REAL(fp)                  , INTENT(IN OUT) :: pcoeff(0:,:)
    TYPE(CSinterp_type)       , INTENT(IN OUT) :: csi

    ! Local variables
    INTEGER :: cloud_state, cloud_loc
    INTEGER :: Error_Status
    INTEGER  :: j, k, l, m
    CHARACTER(ML) :: Message

    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Get_Cloud_Opt_MW'

    ! Initialise results that may
    ! not be interpolated
    ! ---------------------------
    w      = ZERO
    pcoeff = ZERO


    ! Find the frequency, effective radius (or water content)
    ! and temperature indices for interpolation
    ! -----------------------------------------
    csi%f_int = MAX(MIN(CloudC%Frequency_MW(CloudC%n_MW_Frequencies),Frequency),CloudC%Frequency_MW(1))
    CALL find_index(CloudC%Frequency_MW, csi%f_int, csi%i1,csi%i2, csi%f_outbound)
    csi%f = CloudC%Frequency_MW(csi%i1:csi%i2)

    ! If all Reff_MW existed in the CloudCoeff then  CloudC%Reff_MW should be greater than zero and
    ! will use Reff for interpolation otherwise will use the water content
    IF (ALL(CloudC%Reff_MW .GT. ZERO)) THEN
       csi%r_int = MAX(MIN(CloudC%Reff_MW(CloudC%n_MW_Radii),Reff),CloudC%Reff_MW(1))
       CALL find_index(CloudC%Reff_MW, csi%r_int, csi%j1,csi%j2, csi%r_outbound)
       csi%r = CloudC%Reff_MW(csi%j1:csi%j2)

       ! Find the index of the given cloud type (k) in CloudCoeff
       ! The array index starts from zero but findloc starts from 1
       cloud_loc = FINDLOC(CLOUD_TYPE_MIE_TAMU, Cloud_Type, DIM=1) - 1
       cloud_state = CLOUD_STATE_MIE_TAMU(cloud_loc)
       k = CLOUD_INDEX_MIE_TAMU(cloud_loc)
    ELSE
       csi%r_int = MAX(MIN(CloudC%Water_Density_MW(CloudC%n_MW_Radii),Water_Density),CloudC%Water_Density_MW(1))
       CALL find_index(CloudC%Water_Density_MW, csi%r_int, csi%j1,csi%j2, csi%r_outbound)
       csi%r = CloudC%Water_Density_MW(csi%j1:csi%j2)

       ! Find the index of the given cloud type (k) in CloudCoeff
       ! The array index starts from zero but findloc starts from 1
       cloud_loc = FINDLOC(CLOUD_TYPE_DDA_ARTS, Cloud_Type, DIM=1) - 1
       cloud_state = CLOUD_STATE_DDA_ARTS(cloud_loc)
       k = CLOUD_INDEX_DDA_ARTS(cloud_loc)
    ENDIF

    csi%t_int = MAX(MIN(CloudC%Temperature(CloudC%n_Temperatures),Temperature),CloudC%Temperature(1))
    CALL find_index(CloudC%Temperature, csi%t_int, csi%k1,csi%k2, csi%t_outbound)
    csi%t = CloudC%Temperature(csi%k1:csi%k2)

    ! Calculate the interpolating polynomials
    ! ---------------------------------------
    ! Frequency term
    CALL LPoly( csi%f, csi%f_int, &  ! Input
                csi%wlp           )  ! Output
    ! Effective radius term
    CALL LPoly( csi%r, csi%r_int, &  ! Input
                csi%xlp           )  ! Output
    ! Temperature term
    CALL LPoly( csi%t, csi%t_int, &  ! Input
                csi%ylp           )  ! Output

    ! Perform interpolation based on cloud type
    ! -----------------------------------------

    SELECT CASE (cloud_state)
      CASE (LIQUID)
        ! Only 2-D interpolation of extinction coefficient as a
        ! fn. of frequency and temperature for water cloud; i.e.
        ! we're only considering Rayleigh scattering here.
        IF(Cloud_Type .EQ. WATER_CLOUD) THEN
           j = 1
           CALL interp_2D( CloudC%ke_L_MW(csi%i1:csi%i2,j,csi%k1:csi%k2), csi%wlp, csi%ylp, ke )
           CALL interp_2D( CloudC%kb_L_MW(csi%i1:csi%i2,j,csi%k1:csi%k2), csi%wlp, csi%ylp, kb )
        ! All 3-D interpolations for rain cloud!
        ELSE
           CALL interp_3D( CloudC%ke_L_MW(csi%i1:csi%i2,csi%j1:csi%j2,csi%k1:csi%k2), csi%wlp, csi%xlp, csi%ylp, ke )
           CALL interp_3D( CloudC%kb_L_MW(csi%i1:csi%i2,csi%j1:csi%j2,csi%k1:csi%k2), csi%wlp, csi%xlp, csi%ylp, kb )
           CALL interp_3D( CloudC%w_L_MW(csi%i1:csi%i2,csi%j1:csi%j2,csi%k1:csi%k2) , csi%wlp, csi%xlp, csi%ylp, w  )
           IF ( CloudScatter%n_Phase_Elements > 0 .and. CloudScatter%Include_Scattering ) THEN
              pcoeff(0,1) = POINT_5
              DO m = 1, CloudScatter%n_Phase_Elements
                DO l = 1, CloudScatter%n_Legendre_Terms
                  CALL interp_3D( CloudC%pcoeff_L_MW(csi%i1:csi%i2,csi%j1:csi%j2,csi%k1:csi%k2,l+CloudScatter%lOffset,m), &
                                  csi%wlp, csi%xlp, csi%ylp, pcoeff(l,m) )
                END DO
              END DO
           ELSE
             ! Absorption coefficient
             ke = ke * (ONE - w)
           END IF ! Cloud scatter
        ENDIF !Cloud_Type
      CASE (FROZEN)
        IF (Cloud_Type .EQ. ICE_CLOUD) THEN
           j = 1
           CALL interp_1D( CloudC%ke_S_MW(csi%i1:csi%i2,j,k), csi%wlp, ke )
           CALL interp_1D( CloudC%kb_S_MW(csi%i1:csi%i2,j,k), csi%wlp, kb )
        ELSE
           ! Perform interpolation
           CALL interp_2D( CloudC%ke_S_MW(csi%i1:csi%i2,csi%j1:csi%j2,k), csi%wlp, csi%xlp, ke )
           CALL interp_2D( CloudC%kb_S_MW(csi%i1:csi%i2,csi%j1:csi%j2,k), csi%wlp, csi%xlp, kb )
           CALL interp_2D( CloudC%w_S_MW(csi%i1:csi%i2,csi%j1:csi%j2,k) , csi%wlp, csi%xlp, w  )
           IF (CloudScatter%n_Phase_Elements > 0 .and. CloudScatter%Include_Scattering ) THEN
             pcoeff(0,1) = POINT_5
             DO m = 1, CloudScatter%n_Phase_Elements
               DO l = 1, CloudScatter%n_Legendre_Terms
                 CALL interp_2D( CloudC%pcoeff_S_MW(csi%i1:csi%i2,csi%j1:csi%j2,k,l+CloudScatter%lOffset,m), &
                                 csi%wlp, csi%xlp, pcoeff(l,m) )
               END DO
             END DO

           ELSE
           ! Absorption coefficient
             ke = ke * (ONE - w)
           END IF
        END IF ! ICE_Cloud etc interpol
      CASE DEFAULT
          Error_Status = FAILURE
          WRITE(Message,'("Cloud state not defined")')
          CALL Display_Message( ROUTINE_NAME,Message,Error_Status )
          RETURN
    END SELECT

  END SUBROUTINE Get_Cloud_Opt_MW


  ! ---------------------------------------------
  ! Experimental ('CRTM-Exp') MW cloud optics.
  ! Interpolates CloudC_Exp over (Frequency, Dm, Temperature) for the given
  ! habit (cloud_type). The phase-function truncation order (Le) is taken from
  ! the LUT (n_Legendre_Eff) and is DECOUPLED from the RT stream count.
  ! v1: host Effective_Radius is used directly as Dm (identity map); n_Mu=1.
  ! pcoeff convention matches the legacy reader: pcoeff(l,m)=0.5*chi_l, with the
  ! (2l+1) carried in the normalized Legendre polynomials and pcoeff(0,1)=0.5.
  ! ---------------------------------------------
  SUBROUTINE Get_Cloud_Opt_MW_Exp( CloudScatter , &  ! Input
                                   Frequency    , &  ! Input  (GHz)
                                   cloud_type   , &  ! Input
                                   Dm_in        , &  ! Input  Dm proxy (microns)
                                   Temperature  , &  ! Input  (K)
                                   ke           , &  ! Output mass extinction (m^2/kg)
                                   kb           , &  ! Output mass backscatter (m^2/kg)
                                   w            , &  ! Output single-scatter albedo
                                   pcoeff       , &  ! Output phase coefficients
                                   Le           , &  ! Output effective truncation order
                                   csi            )  ! Interpolation
    TYPE(CRTM_AtmOptics_type), INTENT(IN)     :: CloudScatter
    REAL(fp),                  INTENT(IN)     :: Frequency, Dm_in, Temperature
    INTEGER,                   INTENT(IN)     :: cloud_type
    REAL(fp),                  INTENT(OUT)    :: ke, kb, w
    REAL(fp),                  INTENT(IN OUT) :: pcoeff(0:,:)
    INTEGER,                   INTENT(OUT)    :: Le
    TYPE(CSinterp_type),       INTENT(IN OUT) :: csi
    ! Local
    INTEGER  :: h, l, m, np, nfill
    REAL(fp) :: ka, tmp

    ke = ZERO; ka = ZERO; kb = ZERO; w = ZERO; Le = 0
    pcoeff = ZERO

    ! Habit index in the LUT (match the CRTM cloud-type integer)
    h = 0
    DO l = 1, INT(CloudC_Exp%n_Habit)
      IF ( INT(CloudC_Exp%Habit_Id(l)) == cloud_type ) THEN ; h = l ; EXIT ; END IF
    END DO
    IF ( h < 1 ) RETURN     ! habit not in this LUT -> no contribution

    ! Interpolation indices/polynomials: Frequency (wlp), Dm (xlp), Temperature (ylp)
    csi%f_int = MAX(MIN(CloudC_Exp%Frequency(CloudC_Exp%n_Frequency),Frequency),CloudC_Exp%Frequency(1))
    CALL find_index( CloudC_Exp%Frequency, csi%f_int, csi%i1, csi%i2, csi%f_outbound )
    csi%f = CloudC_Exp%Frequency(csi%i1:csi%i2)

    csi%r_int = MAX(MIN(CloudC_Exp%Dm(CloudC_Exp%n_Dm),Dm_in),CloudC_Exp%Dm(1))
    CALL find_index( CloudC_Exp%Dm, csi%r_int, csi%j1, csi%j2, csi%r_outbound )
    csi%r = CloudC_Exp%Dm(csi%j1:csi%j2)

    csi%t_int = MAX(MIN(CloudC_Exp%Temperature(CloudC_Exp%n_Temperature),Temperature),CloudC_Exp%Temperature(1))
    CALL find_index( CloudC_Exp%Temperature, csi%t_int, csi%k1, csi%k2, csi%t_outbound )
    csi%t = CloudC_Exp%Temperature(csi%k1:csi%k2)

    CALL LPoly( csi%f, csi%f_int, csi%wlp )   ! Frequency
    CALL LPoly( csi%r, csi%r_int, csi%xlp )   ! Dm
    CALL LPoly( csi%t, csi%t_int, csi%ylp )   ! Temperature

    ! Bulk optics (mu index = 1 in v1). Cubes are (Freq,T,Dm) -> (wlp,ylp,xlp).
    CALL interp_3D( CloudC_Exp%ke(csi%i1:csi%i2,csi%k1:csi%k2,1,csi%j1:csi%j2,h), csi%wlp, csi%ylp, csi%xlp, ke )
    CALL interp_3D( CloudC_Exp%ka(csi%i1:csi%i2,csi%k1:csi%k2,1,csi%j1:csi%j2,h), csi%wlp, csi%ylp, csi%xlp, ka )
    CALL interp_3D( CloudC_Exp%kb(csi%i1:csi%i2,csi%k1:csi%k2,1,csi%j1:csi%j2,h), csi%wlp, csi%ylp, csi%xlp, kb )
    IF ( ke > ZERO ) w = MAX( MIN( (ke-ka)/ke, ONE ), ZERO )

    IF ( CloudScatter%n_Phase_Elements > 0 .AND. CloudScatter%Include_Scattering ) THEN
      ! n_Legendre_Eff is a COUNT incl. order 0, so the max significant physical
      ! order is (count-1). Output pcoeff index = physical order; the LUT Fortran
      ! Legendre index = physical order + 1 (1-based vs 0-based). Cap so the LUT
      ! read (l+1) and the output index stay in bounds.
      nfill = MAXVAL( CloudC_Exp%n_Legendre_Eff(csi%i1:csi%i2,csi%k1:csi%k2,1,csi%j1:csi%j2,h) ) - 1
      nfill = MIN( nfill, INT(CloudC_Exp%n_Legendre)-1, MAX_N_LEGENDRE_TERMS-1 )
      nfill = MAX( nfill, 0 )
      np = MIN( CloudScatter%n_Phase_Elements, INT(CloudC_Exp%n_Phase_Elements) )
      pcoeff(0,1) = POINT_5
      DO m = 1, np
        DO l = 1, nfill
          CALL interp_3D( CloudC_Exp%pcoeff(m,l+1,csi%i1:csi%i2,csi%k1:csi%k2,1,csi%j1:csi%j2,h), &
                          csi%wlp, csi%ylp, csi%xlp, tmp )
          pcoeff(l,m) = POINT_5 * tmp
        END DO
      END DO
      Le = nfill + 1     ! n_Legendre_Terms: RT sums physical orders 0..nfill
    END IF

  END SUBROUTINE Get_Cloud_Opt_MW_Exp
!
!------------------------------------------------------------------------------
! Tangent-linear / adjoint of the experimental MW cloud optics interpolation.
! Perturbable inputs are Dm (= effective radius) and Temperature; Frequency is
! the (fixed) channel.  Mirrors Get_Cloud_Opt_MW_Exp; uses the saved csi.
!------------------------------------------------------------------------------
  SUBROUTINE Get_Cloud_Opt_MW_Exp_TL( CloudScatter, cloud_type, ke, w,   &  ! FWD Input
                                      Dm_TL, Temperature_TL,             &  ! TL  Input
                                      ke_TL, kb_TL, w_TL, pcoeff_TL,     &  ! TL  Output
                                      csi )                                 ! Interpolation
    TYPE(CRTM_AtmOptics_type), INTENT(IN)     :: CloudScatter
    INTEGER,                   INTENT(IN)     :: cloud_type
    REAL(fp),                  INTENT(IN)     :: ke, w, Dm_TL, Temperature_TL
    REAL(fp),                  INTENT(OUT)    :: ke_TL, kb_TL, w_TL
    REAL(fp),                  INTENT(IN OUT) :: pcoeff_TL(0:,:)
    TYPE(CSinterp_type),       INTENT(IN)     :: csi
    ! Local
    INTEGER  :: h, l, m, np, nfill
    REAL(fp) :: ka, ka_TL, tmp_TL, r_int_TL, t_int_TL
    REAL(fp) :: f_TL(NPTS), r_TL(NPTS), t_TL(NPTS), z3_TL(NPTS,NPTS,NPTS)
    TYPE(LPoly_type) :: wlp_TL, xlp_TL, ylp_TL

    ke_TL = ZERO; kb_TL = ZERO; w_TL = ZERO; pcoeff_TL = ZERO

    h = 0
    DO l = 1, INT(CloudC_Exp%n_Habit)
      IF ( INT(CloudC_Exp%Habit_Id(l)) == cloud_type ) THEN ; h = l ; EXIT ; END IF
    END DO
    IF ( h < 1 ) RETURN

    ! TL of the interpolation point.  Frequency is fixed (channel) -> zero TL; Dm and
    ! Temperature get zero TL when the value was clamped to the LUT bounds (outbound).
    f_TL = ZERO ; r_TL = ZERO ; t_TL = ZERO ; z3_TL = ZERO
    r_int_TL = MERGE( ZERO, Dm_TL,          csi%r_outbound )
    t_int_TL = MERGE( ZERO, Temperature_TL, csi%t_outbound )

    CALL LPoly_TL( csi%f, csi%f_int, csi%wlp, f_TL, ZERO,     wlp_TL )
    CALL LPoly_TL( csi%r, csi%r_int, csi%xlp, r_TL, r_int_TL, xlp_TL )
    CALL LPoly_TL( csi%t, csi%t_int, csi%ylp, t_TL, t_int_TL, ylp_TL )

    ! Cubes are (Freq,T,Dm) -> interp order (wlp,ylp,xlp), matching the forward.
    CALL interp_3D_TL( CloudC_Exp%ke(csi%i1:csi%i2,csi%k1:csi%k2,1,csi%j1:csi%j2,h), &
                       csi%wlp, csi%ylp, csi%xlp, z3_TL, wlp_TL, ylp_TL, xlp_TL, ke_TL )
    CALL interp_3D_TL( CloudC_Exp%ka(csi%i1:csi%i2,csi%k1:csi%k2,1,csi%j1:csi%j2,h), &
                       csi%wlp, csi%ylp, csi%xlp, z3_TL, wlp_TL, ylp_TL, xlp_TL, ka_TL )
    CALL interp_3D_TL( CloudC_Exp%kb(csi%i1:csi%i2,csi%k1:csi%k2,1,csi%j1:csi%j2,h), &
                       csi%wlp, csi%ylp, csi%xlp, z3_TL, wlp_TL, ylp_TL, xlp_TL, kb_TL )

    ! w = (ke-ka)/ke, clamped to [0,1] in the forward -> zero TL when clamped.
    IF ( ke > ZERO .AND. w > ZERO .AND. w < ONE ) THEN
      ka   = ke * ( ONE - w )
      w_TL = ka*ke_TL/(ke*ke) - ka_TL/ke
    END IF

    IF ( CloudScatter%n_Phase_Elements > 0 .AND. CloudScatter%Include_Scattering ) THEN
      nfill = MAXVAL( CloudC_Exp%n_Legendre_Eff(csi%i1:csi%i2,csi%k1:csi%k2,1,csi%j1:csi%j2,h) ) - 1
      nfill = MIN( nfill, INT(CloudC_Exp%n_Legendre)-1, MAX_N_LEGENDRE_TERMS-1 )
      nfill = MAX( nfill, 0 )
      np = MIN( CloudScatter%n_Phase_Elements, INT(CloudC_Exp%n_Phase_Elements) )
      DO m = 1, np
        DO l = 1, nfill
          CALL interp_3D_TL( CloudC_Exp%pcoeff(m,l+1,csi%i1:csi%i2,csi%k1:csi%k2,1,csi%j1:csi%j2,h), &
                             csi%wlp, csi%ylp, csi%xlp, z3_TL, wlp_TL, ylp_TL, xlp_TL, tmp_TL )
          pcoeff_TL(l,m) = POINT_5 * tmp_TL
        END DO
      END DO
    END IF
  END SUBROUTINE Get_Cloud_Opt_MW_Exp_TL
!
  SUBROUTINE Get_Cloud_Opt_MW_Exp_AD( CloudScatter, cloud_type, ke, w,   &  ! FWD Input
                                      ke_AD, kb_AD, w_AD, pcoeff_AD,     &  ! AD  Input (consumed)
                                      Dm_AD, Temperature_AD,             &  ! AD  Output (accumulated)
                                      csi )                                 ! Interpolation
    TYPE(CRTM_AtmOptics_type), INTENT(IN)     :: CloudScatter
    INTEGER,                   INTENT(IN)     :: cloud_type
    REAL(fp),                  INTENT(IN)     :: ke, w
    REAL(fp),                  INTENT(IN OUT) :: ke_AD, kb_AD, w_AD
    REAL(fp),                  INTENT(IN OUT) :: pcoeff_AD(0:,:)
    REAL(fp),                  INTENT(IN OUT) :: Dm_AD, Temperature_AD
    TYPE(CSinterp_type),       INTENT(IN)     :: csi
    ! Local
    INTEGER  :: h, l, m, np, nfill
    REAL(fp) :: ka, ka_AD, tmp_AD, r_int_AD, t_int_AD
    REAL(fp) :: f_AD(NPTS), r_AD(NPTS), t_AD(NPTS), z3_AD(NPTS,NPTS,NPTS)
    TYPE(LPoly_type) :: wlp_AD, xlp_AD, ylp_AD

    h = 0
    DO l = 1, INT(CloudC_Exp%n_Habit)
      IF ( INT(CloudC_Exp%Habit_Id(l)) == cloud_type ) THEN ; h = l ; EXIT ; END IF
    END DO
    IF ( h < 1 ) THEN     ! habit not in LUT -> consume input adjoints, no contribution
      ke_AD = ZERO; kb_AD = ZERO; w_AD = ZERO; pcoeff_AD = ZERO
      RETURN
    END IF

    ! Initialise local adjoints
    z3_AD = ZERO; f_AD = ZERO; r_AD = ZERO; t_AD = ZERO
    r_int_AD = ZERO; t_int_AD = ZERO; ka_AD = ZERO
    CALL Clear_LPoly( wlp_AD ); CALL Clear_LPoly( xlp_AD ); CALL Clear_LPoly( ylp_AD )

    ! Adjoint of the phase-coefficient interpolation
    IF ( CloudScatter%n_Phase_Elements > 0 .AND. CloudScatter%Include_Scattering ) THEN
      nfill = MAXVAL( CloudC_Exp%n_Legendre_Eff(csi%i1:csi%i2,csi%k1:csi%k2,1,csi%j1:csi%j2,h) ) - 1
      nfill = MIN( nfill, INT(CloudC_Exp%n_Legendre)-1, MAX_N_LEGENDRE_TERMS-1 )
      nfill = MAX( nfill, 0 )
      np = MIN( CloudScatter%n_Phase_Elements, INT(CloudC_Exp%n_Phase_Elements) )
      DO m = 1, np
        DO l = 1, nfill
          tmp_AD = POINT_5 * pcoeff_AD(l,m)
          CALL interp_3D_AD( CloudC_Exp%pcoeff(m,l+1,csi%i1:csi%i2,csi%k1:csi%k2,1,csi%j1:csi%j2,h), &
                             csi%wlp, csi%ylp, csi%xlp, tmp_AD, z3_AD, wlp_AD, ylp_AD, xlp_AD )
        END DO
      END DO
    END IF
    pcoeff_AD = ZERO

    ! Adjoint of w = (ke-ka)/ke (clamped) -> ke_AD, ka_AD
    IF ( ke > ZERO .AND. w > ZERO .AND. w < ONE ) THEN
      ka    = ke * ( ONE - w )
      ke_AD = ke_AD + ka*w_AD/(ke*ke)
      ka_AD = ka_AD - w_AD/ke
    END IF
    w_AD = ZERO

    ! Adjoint of the ke/ka/kb interpolation (accumulate into the LPoly adjoints)
    CALL interp_3D_AD( CloudC_Exp%kb(csi%i1:csi%i2,csi%k1:csi%k2,1,csi%j1:csi%j2,h), &
                       csi%wlp, csi%ylp, csi%xlp, kb_AD, z3_AD, wlp_AD, ylp_AD, xlp_AD )
    kb_AD = ZERO
    CALL interp_3D_AD( CloudC_Exp%ka(csi%i1:csi%i2,csi%k1:csi%k2,1,csi%j1:csi%j2,h), &
                       csi%wlp, csi%ylp, csi%xlp, ka_AD, z3_AD, wlp_AD, ylp_AD, xlp_AD )
    CALL interp_3D_AD( CloudC_Exp%ke(csi%i1:csi%i2,csi%k1:csi%k2,1,csi%j1:csi%j2,h), &
                       csi%wlp, csi%ylp, csi%xlp, ke_AD, z3_AD, wlp_AD, ylp_AD, xlp_AD )
    ke_AD = ZERO

    ! Adjoint of the interpolating polynomials -> interp-point adjoints
    CALL LPoly_AD( csi%t, csi%t_int, csi%ylp, ylp_AD, t_AD, t_int_AD )
    CALL LPoly_AD( csi%r, csi%r_int, csi%xlp, xlp_AD, r_AD, r_int_AD )
    ! (frequency is fixed -> wlp_AD discarded)

    ! Map interp-point adjoints to the cloud-state inputs (zero when clamped/outbound)
    IF ( .NOT. csi%t_outbound ) Temperature_AD = Temperature_AD + t_int_AD
    IF ( .NOT. csi%r_outbound ) Dm_AD          = Dm_AD          + r_int_AD
  END SUBROUTINE Get_Cloud_Opt_MW_Exp_AD


  ! ---------------------------------------------
  ! Subroutine to obtain the tangent-linear
  ! MW bulk optical properties of a cloud:
  !   extinction coefficient (ke_TL),
  !   scattereing coefficient (w_TL)
  !   spherical Legendre coefficients (pcoeff_TL)
  ! ---------------------------------------------
  SUBROUTINE Get_Cloud_Opt_MW_TL( CloudScatter_TL  , &  ! Input  CloudScatter TL structure
                                  cloud_type       , &  ! Input  see CRTM_Cloud_Define.f90
                                  ke               , &  ! Input
                                  kb               , &  ! Input
                                  w                , &  ! Input
                                  Is_Active_Sensor , & ! Input
                                  Water_Density_TL , &  ! TL  Input  effective radius (mm)
                                  Reff_TL          , &  ! TL  Input  effective radius (mm)
                                  Temperature_TL   , &  ! TL  Input  cloudy temperature
                                  ke_TL            , &  ! TL  Output extinction coefficient (=~ optical depth for 1 mm water content)
                                  kb_TL            , &  ! TL  Output Backscatering coefficient
                                  w_TL             , &  ! TL  Output single scattering albedo
                                  pcoeff_TL        , &  ! TL  Output spherical Legendre coefficients
                                  csi                )  ! Input interpolation data
    ! Arguments
    TYPE(CRTM_AtmOptics_type), INTENT(IN)     :: CloudScatter_TL
    INTEGER ,                   INTENT(IN)     :: Cloud_Type
    REAL(fp),                   INTENT(IN)     :: ke, kb, w, Reff_TL
    REAL(fp),                   INTENT(IN)     :: Water_Density_TL
    REAL(fp),                   INTENT(IN)     :: Temperature_TL
    LOGICAL,                    INTENT(IN)     :: Is_Active_Sensor
    REAL(fp),                   INTENT(OUT)    :: ke_TL, kb_TL
    REAL(fp),                   INTENT(OUT)    :: w_TL
    REAL(fp),                   INTENT(IN OUT) :: pcoeff_TL(0:,:)
    TYPE(CSinterp_type),        INTENT(IN)     :: csi
    ! Local variables
    INTEGER  :: j, k, l, m
    REAL(fp) :: f_int_TL, r_int_TL, t_int_TL
    REAL(fp) :: f_TL(NPTS), r_TL(NPTS), t_TL(NPTS)
    REAL(fp) :: z2_TL(NPTS,NPTS)
    REAL(fp) :: z3_TL(NPTS,NPTS,NPTS)
    TYPE(LPoly_type) :: wlp_TL, xlp_TL, ylp_TL

    ! Variables to find the cloud index and state
    INTEGER :: cloud_state, cloud_loc

!JR Static initialization means only 1 copy of the variable. OpenMP over profiles
!JR means $OPENMP_NUM_THREADS copies are needed. So change to run-time initialization
!JR    REAL(fp), POINTER :: z2(:,:) => NULL()
!JR    REAL(fp), POINTER :: z3(:,:,:) => NULL()
    REAL(fp), POINTER :: z2(:,:)
    REAL(fp), POINTER :: z3(:,:,:)

    NULLIFY(z2)
    NULLIFY(z3)
    ! Setup
    ! -----
    ! Initialise results that may
    ! not be interpolated
    w_TL      = ZERO
    pcoeff_TL = ZERO
    ! The local TL inputs
    f_int_TL = ZERO
    f_TL     = ZERO
    r_TL     = ZERO
    t_int_TL = Temperature_TL
    t_TL     = ZERO
    z2_TL = ZERO
    z3_TL = ZERO

     ! If all Reff_MW existed in the CloudCoeff then  CloudC%Reff_MW should be greater than zero and
    ! will use Reff for interpolation otherwise will use the water content
    IF (ALL(CloudC%Reff_MW .GT. ZERO)) THEN
       ! Find the index of the given cloud type (k) in CloudCoeff
       ! The array index starts from zero but findloc starts from 1
       cloud_loc = FINDLOC(CLOUD_TYPE_MIE_TAMU, Cloud_Type, DIM=1) - 1
       cloud_state = CLOUD_STATE_MIE_TAMU(cloud_loc)
       k = CLOUD_INDEX_MIE_TAMU(cloud_loc)
       r_int_TL = Reff_TL
    ELSE
       ! Find the index of the given cloud type (k) in CloudCoeff
       ! The array index starts from zero but findloc starts from 1
       cloud_loc = FINDLOC(CLOUD_TYPE_DDA_ARTS, Cloud_Type, DIM=1) - 1
       cloud_state = CLOUD_STATE_DDA_ARTS(cloud_loc)
       k = CLOUD_INDEX_DDA_ARTS(cloud_loc)
       r_int_TL = Water_Density_TL
    END IF


    ! Calculate the TL interpolating polynomials
    ! ------------------------------------------
    ! Frequency term (always zero. This is a placeholder for testing)
    CALL LPoly_TL( csi%f, csi%f_int, & ! FWD Input
                   csi%wlp,          & ! FWD Input
                   f_TL, f_int_TL,   & ! TL  Input
                   wlp_TL            ) ! TL  Output
    ! size term: effective radius or water content
    CALL LPoly_TL( csi%r, csi%r_int, & ! FWD Input
                   csi%xlp,          & ! FWD Input
                   r_TL, r_int_TL,   & ! TL  Input
                   xlp_TL            ) ! TL  Output
    ! Temperature term
    CALL LPoly_TL( csi%t, csi%t_int, & ! FWD Input
                   csi%ylp,          & ! FWD Input
                   t_TL, t_int_TL,   & ! TL  Input
                   ylp_TL            ) ! TL  Output

    ! Perform interpolation based on cloud type
    ! -----------------------------------------
    SELECT CASE (cloud_state)
     CASE (LIQUID)
        IF(Cloud_Type .EQ. WATER_CLOUD) THEN
          ! Only 2-D interpolation of extinction coefficient as a
          ! fn. of frequency and temperature for water cloud; i.e.
          ! we're only considering Rayleigh scattering here.
          ! No TL output when all dimensions
          ! are outside LUT bounds
          IF ( csi%f_outbound .AND. csi%t_outbound ) THEN
             ke_TL = ZERO
             kb_TL = ZERO
             RETURN
          END IF
          j = 1
          z2 => CloudC%ke_L_MW(csi%i1:csi%i2,j,csi%k1:csi%k2)
          CALL interp_2D_TL( z2   , csi%wlp, csi%ylp, &  ! FWD Input
                             z2_TL, wlp_TL , ylp_TL , &  ! TL  Input
                             ke_TL                    )  ! TL  Output
          z2 => CloudC%kb_L_MW(csi%i1:csi%i2,j,csi%k1:csi%k2)
          CALL interp_2D_TL( z2   , csi%wlp, csi%ylp, &  ! FWD Input
                             z2_TL, wlp_TL , ylp_TL , &  ! TL  Input
                             kb_TL                    )  ! TL  Output
        ELSE ! Rain_Cloud
           ! All 3-D interpolations for rain cloud!
           ! No TL output when all dimensions
           ! are outside LUT bounds
           IF ( csi%f_outbound .AND. csi%r_outbound .AND. csi%t_outbound ) THEN
              ke_TL = ZERO
              kb_TL = ZERO
              RETURN
           END IF
           ! Extinction coefficient
           z3 => CloudC%ke_L_MW(csi%i1:csi%i2,csi%j1:csi%j2,csi%k1:csi%k2)
           CALL interp_3D_TL( z3   , csi%wlp, csi%xlp, csi%ylp, &  ! FWD Input
                              z3_TL, wlp_TL , xlp_TL , ylp_TL , &  ! TL  Input
                              ke_TL                           )  ! TL  Output
           ! Backscatering coefficient
           IF ( Is_Active_Sensor .and. CloudScatter_TL%Include_Scattering ) THEN
              z3 => CloudC%kb_L_MW(csi%i1:csi%i2,csi%j1:csi%j2,csi%k1:csi%k2)
              CALL interp_3D_TL( z3   , csi%wlp, csi%xlp, csi%ylp, &  ! FWD Input
                                 z3_TL, wlp_TL , xlp_TL , ylp_TL , &  ! TL  Input
                                 kb_TL                           )  ! TL  Output
           END IF
           ! Single scatter albedo
           z3 => CloudC%w_L_MW(csi%i1:csi%i2,csi%j1:csi%j2,csi%k1:csi%k2)
           CALL interp_3D_TL( z3   , csi%wlp, csi%xlp, csi%ylp, &  ! FWD Input
                             z3_TL, wlp_TL , xlp_TL , ylp_TL , &  ! TL  Input
                             w_TL                              )  ! TL  Output
           ! Phase matrix coefficients
           IF ( CloudScatter_TL%n_Phase_Elements > 0 .and. CloudScatter_TL%Include_Scattering ) THEN
             pcoeff_TL(0,1) = ZERO
             DO m = 1, CloudScatter_TL%n_Phase_Elements
               DO l = 1, CloudScatter_TL%n_Legendre_Terms
                 z3 => CloudC%pcoeff_L_MW(csi%i1:csi%i2,csi%j1:csi%j2,csi%k1:csi%k2,l+CloudScatter_TL%lOffset,m)
                 CALL interp_3D_TL( z3   , csi%wlp, csi%xlp, csi%ylp, &  ! FWD Input
                                    z3_TL, wlp_TL , xlp_TL , ylp_TL , &  ! TL  Input
                                    pcoeff_TL(l,m)                    )  ! TL  Output
               END DO
             END DO
           ELSE
             ! Absorption coefficient
             IF( w < ONE ) THEN
                ke_TL = ke_TL * (ONE - w) - ke/(ONE -w) * w_TL
             ELSE
                ke_TL = ZERO
             END IF
           END IF
        END IF ! Cloud_Type
     CASE (FROZEN)
        IF (Cloud_Type .EQ. ICE_CLOUD) THEN
            ! No TL interpolation of extinction coefficient as it
            ! is only a fn. of frequency for ice cloud
            ke_TL = ZERO
            kb_TL = ZERO
        ELSE ! All other frozen habits
            ! The remaining cloud types have 2-D interpolation
            ! as a fn. of frequency and radius
            ! No TL output when all dimensions
            ! are outside LUT bounds
            IF ( csi%f_outbound .AND. csi%r_outbound ) THEN
              ke_TL = ZERO
              kb_TL = ZERO
              RETURN
            END IF

            ! Extinction coefficient
            z2 => CloudC%ke_S_MW(csi%i1:csi%i2,csi%j1:csi%j2,k)
            CALL interp_2D_TL( z2   , csi%wlp, csi%xlp, &  ! FWD Input
                               z2_TL, wlp_TL , xlp_TL , &  ! TL  Input
                               ke_TL                    )  ! TL  Output
            ! Backscatering coefficient
            IF (Is_Active_Sensor .and. CloudScatter_TL%Include_Scattering ) THEN
               z2 => CloudC%kb_S_MW(csi%i1:csi%i2,csi%j1:csi%j2,k)
               CALL interp_2D_TL( z2   , csi%wlp, csi%xlp, &  ! FWD Input
                                  z2_TL, wlp_TL , xlp_TL , &  ! TL  Input
                                  kb_TL                    )  ! TL  Output
            END IF
            ! Single scatter albedo
            z2 => CloudC%w_S_MW(csi%i1:csi%i2,csi%j1:csi%j2,k)
            CALL interp_2D_TL( z2   , csi%wlp, csi%xlp, &  ! FWD Input
                               z2_TL, wlp_TL , xlp_TL , &  ! TL  Input
                               w_TL                     )  ! TL  Output
            ! Phase matrix coefficients
            IF ( CloudScatter_TL%n_Phase_Elements > 0 .and. CloudScatter_TL%Include_Scattering ) THEN
              pcoeff_TL(0,1) = ZERO
              DO m = 1, CloudScatter_TL%n_Phase_Elements
                DO l = 1, CloudScatter_TL%n_Legendre_Terms
                  z2 => CloudC%pcoeff_S_MW(csi%i1:csi%i2,csi%j1:csi%j2,k,l+CloudScatter_TL%lOffset,m)
                  CALL interp_2D_TL( z2   , csi%wlp, csi%xlp, &  ! FWD Input
                                     z2_TL, wlp_TL , xlp_TL , &  ! TL  Input
                                     pcoeff_TL(l,m)           )  ! TL  Output
                END DO
              END DO
            ELSE
              ! Absorption coefficient
              IF( w < ONE ) THEN
                 ke_TL = ke_TL * (ONE - w) - ke/(ONE -w) * w_TL
              ELSE
                 ke_TL = ZERO
              END IF
            END IF
        END IF
    END SELECT

    NULLIFY(z2, z3)

  END SUBROUTINE Get_Cloud_Opt_MW_TL


  ! ---------------------------------------------
  ! Subroutine to obtain the adjoint of the
  ! MW bulk optical properties of a cloud:
  !   effective radius (Reff_AD),
  !   temperature (temperature_AD)
  ! ---------------------------------------------
  SUBROUTINE Get_Cloud_Opt_MW_AD(CloudScatter_AD  , &  ! Input      CloudScatter AD structure
                                 cloud_type       , &  ! Input      see CRTM_Cloud_Define.f90
                                 ke               , &  ! Input
                                 kb               , &  ! Input
                                 w                , &  ! Input
                                 Is_Active_Sensor , &  ! Input
                                 ke_AD            , &  ! AD  Input  extinction coefficient (=~ optical depth for 1 mm water content)
                                 kb_AD            , &  ! AD  Input  Backscatering coefficient
                                 w_AD             , &  ! AD  Input  single scattering albedo
                                 pcoeff_AD        , &  ! AD  Input  spherical Legendre coefficients
                                 Water_Density_AD , &  ! AD  Output water content (kg.m^-3)
                                 Reff_AD          , &  ! AD  Output effective radius (mm)
                                 Temperature_AD   , &  ! AD  Output temperature
                                 csi                 ) ! Input interpolation data
    ! Arguments
    TYPE(CRTM_AtmOptics_type), INTENT(IN)     :: CloudScatter_AD
    INTEGER ,                   INTENT(IN)     :: Cloud_Type
    REAL(fp),                   INTENT(IN)     :: ke, kb, w
    LOGICAL,                    INTENT(IN)     :: Is_Active_Sensor
    REAL(fp),                   INTENT(IN OUT) :: ke_AD, kb_AD      ! AD  Input
    REAL(fp),                   INTENT(IN OUT) :: w_AD              ! AD  Input
    REAL(fp),                   INTENT(IN OUT) :: pcoeff_AD(0:,:)   ! AD  Input
    REAL(fp),                   INTENT(IN OUT) :: Water_Density_AD  ! AD  Output
    REAL(fp),                   INTENT(IN OUT) :: Reff_AD           ! AD  Output
    REAL(fp),                   INTENT(IN OUT) :: Temperature_AD    ! AD  Output
    TYPE(CSinterp_type),        INTENT(IN)     :: csi
    ! Local variables
    INTEGER  :: j, k, l, m
    REAL(fp) :: f_int_AD, r_int_AD, t_int_AD
    REAL(fp) :: f_AD(NPTS), r_AD(NPTS), t_AD(NPTS)
    REAL(fp) :: z2_AD(NPTS,NPTS)
    REAL(fp) :: z3_AD(NPTS,NPTS,NPTS)
    TYPE(LPoly_type) :: wlp_AD, xlp_AD, ylp_AD
!JR Static initialization means only 1 copy of the variable. OpenMP over profiles
!JR means $OPENMP_NUM_THREADS copies are needed. So change to run-time initialization
!JR    REAL(fp), POINTER :: z2(:,:) => NULL()
!JR    REAL(fp), POINTER :: z3(:,:,:) => NULL()
    REAL(fp), POINTER :: z2(:,:)
    REAL(fp), POINTER :: z3(:,:,:)

    ! Cloud index and state in CloudCoeff
    INTEGER :: cloud_state, cloud_loc

    NULLIFY(z2)
    NULLIFY(z3)
    ! Setup
    ! -----
    ! Initialise local adjoint variables
    f_int_AD = ZERO
    f_AD     = ZERO
    r_int_AD = ZERO
    r_AD     = ZERO
    t_int_AD = ZERO
    t_AD     = ZERO
    z2_AD = ZERO
    z3_AD = ZERO
    CALL Clear_LPoly(wlp_AD)
    CALL Clear_LPoly(xlp_AD)
    CALL Clear_LPoly(ylp_AD)

    ! If all Reff_MW existed in the CloudCoeff then  CloudC%Reff_MW should be greater than zero and
    ! will use Reff for interpolation otherwise will use the water content
    IF (ALL(CloudC%Reff_MW .GT. ZERO)) THEN
       ! Find the index of the given cloud type (k) in CloudCoeff
       ! The array index starts from zero but findloc starts from 1
       cloud_loc = FINDLOC(CLOUD_TYPE_MIE_TAMU, Cloud_Type, DIM=1) - 1
       cloud_state = CLOUD_STATE_MIE_TAMU(cloud_loc)
       k = CLOUD_INDEX_MIE_TAMU(cloud_loc)
    ELSE
       ! Find the index of the given cloud type (k) in CloudCoeff
       ! The array index starts from zero but findloc starts from 1
       cloud_loc = FINDLOC(CLOUD_TYPE_DDA_ARTS, Cloud_Type, DIM=1) - 1
       cloud_state = CLOUD_STATE_DDA_ARTS(cloud_loc)
       k = CLOUD_INDEX_DDA_ARTS(cloud_loc)
    END IF


    ! Perform interpolation based on cloud type
    ! -----------------------------------------
    SELECT CASE (cloud_state)
     CASE (LIQUID)
        IF(Cloud_Type .EQ. WATER_CLOUD) THEN
           ! Only 2-D interpolation of extinction coefficient as a
           ! fn. of frequency and temperature for water cloud; i.e.
           ! we're only considering Rayleigh scattering here.
           ! ------------------------------------------------------
           ! No AD output when all dimensions
           ! are outside LUT bounds
           IF ( csi%f_outbound .AND. csi%t_outbound ) THEN
             ke_AD     = ZERO
             kb_AD     = ZERO
             w_AD      = ZERO
             pcoeff_AD = ZERO
             RETURN
           END IF
           ! Perform the AD interpolations
           j = 1
           z2 => CloudC%ke_L_MW(csi%i1:csi%i2,j,csi%k1:csi%k2)
           CALL interp_2D_AD( z2   , csi%wlp, csi%ylp, &  ! FWD Input
                              ke_AD                  , &  ! AD  Input
                              z2_AD, wlp_AD, ylp_AD    )  ! AD  Output

           IF (Is_Active_Sensor .and. CloudScatter_AD%Include_Scattering ) THEN
              z2 => CloudC%kb_L_MW(csi%i1:csi%i2,j,csi%k1:csi%k2)
              CALL interp_2D_AD( z2   , csi%wlp, csi%ylp, &  ! FWD Input
                                 kb_AD                  , &  ! AD  Input
                                 z2_AD, wlp_AD, ylp_AD    )  ! AD  Output
           END IF
           ! Compute the AD of the interpolating polynomials
           ! Temperature term
           CALL LPoly_AD( csi%t, csi%t_int, & ! FWD Input
                          csi%ylp,          & ! FWD Input
                          ylp_AD,           & ! AD  Input
                          t_AD, t_int_AD    ) ! AD  Output
           ! Frequency term (always zero. This is a placeholder for testing)
           CALL LPoly_AD( csi%f, csi%f_int, & ! FWD Input
                          csi%wlp,          & ! FWD Input
                          wlp_AD,           & ! AD  Input
                          f_AD, f_int_AD    ) ! AD  Output
           ! The AD outputs
           Temperature_AD = Temperature_AD + t_int_AD
        ! All 3-D interpolations for rain cloud!
        ! --------------------------------------
        ELSE
           ! No AD output when all dimensions
           ! are outside LUT bounds
           IF ( csi%f_outbound .AND. csi%r_outbound .AND. csi%t_outbound ) THEN
             ke_AD     = ZERO
             kb_AD     = ZERO
             w_AD      = ZERO
             pcoeff_AD = ZERO
             RETURN
           END IF
           ! Perform the AD interpolations
           ! Phase matrix coefficients
           IF (CloudScatter_AD%n_Phase_Elements > 0 .and. CloudScatter_AD%Include_Scattering ) THEN
             DO m = 1, CloudScatter_AD%n_Phase_Elements
               DO l = 1, CloudScatter_AD%n_Legendre_Terms
                 z3 => CloudC%pcoeff_L_MW(csi%i1:csi%i2,csi%j1:csi%j2,csi%k1:csi%k2,l+CloudScatter_AD%lOffset,m)
                 CALL interp_3D_AD( z3   , csi%wlp, csi%xlp, csi%ylp, &  ! FWD Input
                                    pcoeff_AD(l,m)                  , &  ! AD  Input
                                    z3_AD, wlp_AD , xlp_AD , ylp_AD   )  ! AD  Output
               END DO
             END DO
             pcoeff_AD(0,1) = ZERO

           ELSE
              ! Absorption coefficient
              IF( w < ONE ) THEN
                w_AD  = w_AD - ke/(ONE -w) * ke_AD
                ke_AD = ke_AD * (ONE - w)
              ELSE
                ke_AD = ZERO
              END IF
           END IF
           ! Single scatter albedo
           z3 => CloudC%w_L_MW(csi%i1:csi%i2,csi%j1:csi%j2,csi%k1:csi%k2)
           CALL interp_3D_AD( z3   , csi%wlp, csi%xlp, csi%ylp, &  ! FWD Input
                              w_AD                            , &  ! AD  Input
                              z3_AD, wlp_AD , xlp_AD , ylp_AD   )  ! AD  Output
           ! Extinction coefficient
           z3 => CloudC%ke_L_MW(csi%i1:csi%i2,csi%j1:csi%j2,csi%k1:csi%k2)
           CALL interp_3D_AD( z3   , csi%wlp, csi%xlp, csi%ylp, &  ! FWD Input
                              ke_AD                           , &  ! AD  Input
                              z3_AD, wlp_AD , xlp_AD , ylp_AD   )  ! AD  Output
           ! Backscatering coefficient
           IF (Is_Active_Sensor .and. CloudScatter_AD%Include_Scattering ) THEN
              z3 => CloudC%kb_L_MW(csi%i1:csi%i2,csi%j1:csi%j2,csi%k1:csi%k2)
              CALL interp_3D_AD( z3   , csi%wlp, csi%xlp, csi%ylp, &  ! FWD Input
                                 kb_AD                           , &  ! AD  Input
                                 z3_AD, wlp_AD , xlp_AD , ylp_AD   )  ! AD  Output
           END IF
           ! Compute the AD of the interpolating polynomials
           ! Temperature term
           CALL LPoly_AD( csi%t, csi%t_int, & ! FWD Input
                          csi%ylp,          & ! FWD Input
                          ylp_AD,           & ! AD  Input
                          t_AD, t_int_AD    ) ! AD  Output
           ! Effective radius term
           CALL LPoly_AD( csi%r, csi%r_int, & ! FWD Input
                          csi%xlp,          & ! FWD Input
                          xlp_AD,           & ! AD  Input
                          r_AD, r_int_AD    ) ! AD  Output
           ! Frequency term (always zero. This is a placeholder for testing)
           CALL LPoly_AD( csi%f, csi%f_int, & ! FWD Input
                          csi%wlp,          & ! FWD Input
                          wlp_AD,           & ! AD  Input
                          f_AD, f_int_AD    ) ! AD  Output
           ! The AD outputs
           Temperature_AD = Temperature_AD + t_int_AD
           IF (ALL(CloudC%Reff_MW .GT. ZERO)) THEN
              Reff_AD = Reff_AD + r_int_AD
           ELSE
              Water_Density_AD = Water_Density_AD + r_int_AD
           END IF
        END IF
     CASE (FROZEN)
        IF (Cloud_Type .EQ. ICE_CLOUD) THEN
           ! No AD interpolation as it is only a fn.
           ! of frequency for ice cloud
           ! ---------------------------------------
           ke_AD     = ZERO
           kb_AD     = ZERO
           w_AD      = ZERO
           pcoeff_AD = ZERO
        ELSE ! All other frozen habits
           ! The remaining cloud types have 2-D interpolation
           ! as a fn. of frequency and radius
           ! ------------------------------------------------
           ! No TL output when all dimensions
           ! are outside LUT bounds
           IF ( csi%f_outbound .AND. csi%r_outbound ) THEN
             ke_AD     = ZERO
             kb_AD     = ZERO
             w_AD      = ZERO
             pcoeff_AD = ZERO
             RETURN
           END IF
           ! Perform the AD interpolations
           ! Phase matrix coefficients
           IF (CloudScatter_AD%n_Phase_Elements > 0 .and. CloudScatter_AD%Include_Scattering ) THEN
             DO m = 1, CloudScatter_AD%n_Phase_Elements
               DO l = 1, CloudScatter_AD%n_Legendre_Terms
                 z2 => CloudC%pcoeff_S_MW(csi%i1:csi%i2,csi%j1:csi%j2,k,l+CloudScatter_AD%lOffset,m)
                 CALL interp_2D_AD( z2   , csi%wlp, csi%xlp, &  ! FWD Input
                                    pcoeff_AD(l,m)         , &  ! AD  Input
                                    z2_AD, wlp_AD , xlp_AD   )  ! AD  Output
               END DO
             END DO
             pcoeff_AD(0,1) = ZERO
           ELSE
             ! Absorption coefficient
             IF( w < ONE ) THEN
                w_AD  = w_AD - ke/(ONE -w) * ke_AD
                ke_AD = ke_AD * (ONE - w)
             ELSE
                 ke_AD = ZERO
             END IF
           END IF
           ! Single scatter albedo
           z2 => CloudC%w_S_MW(csi%i1:csi%i2,csi%j1:csi%j2,k)
           CALL interp_2D_AD( z2   , csi%wlp, csi%xlp, &  ! FWD Input
                              w_AD                   , &  ! AD  Input
                              z2_AD, wlp_AD , xlp_AD   )  ! AD  Output
           ! Extinction coefficient
           z2 => CloudC%ke_S_MW(csi%i1:csi%i2,csi%j1:csi%j2,k)
           CALL interp_2D_AD( z2   , csi%wlp, csi%xlp, &  ! FWD Input
                              ke_AD                  , &  ! AD  Input
                              z2_AD, wlp_AD , xlp_AD   )  ! AD  Output
           ! Backscatering coefficient
           IF ( Is_Active_Sensor .and. CloudScatter_AD%Include_Scattering ) THEN
              z2 => CloudC%kb_S_MW(csi%i1:csi%i2,csi%j1:csi%j2,k)
              CALL interp_2D_AD( z2   , csi%wlp, csi%xlp, &  ! FWD Input
                                 kb_AD                  , &  ! AD  Input
                                 z2_AD, wlp_AD , xlp_AD   )  ! AD  Output
           END IF
           ! Compute the AD of the interpolating polynomials
           ! Effective radius term
           CALL LPoly_AD( csi%r, csi%r_int, & ! FWD Input
                          csi%xlp,          & ! FWD Input
                          xlp_AD,           & ! AD  Input
                          r_AD, r_int_AD    ) ! AD  Output
           ! Frequency term (always zero. This is a placeholder for testing)
           CALL LPoly_AD( csi%f, csi%f_int, & ! FWD Input
                          csi%wlp,          & ! FWD Input
                          wlp_AD,           & ! AD  Input
                          f_AD, f_int_AD    ) ! AD  Output
           ! The AD outputs
           IF (ALL(CloudC%Reff_MW .GT. ZERO)) THEN
              Reff_AD = Reff_AD + r_int_AD
           ELSE
              Water_Density_AD = Water_Density_AD + r_int_AD
           END IF
        END IF ! Cloud_Type
    END SELECT
    NULLIFY(z2, z3)

  END SUBROUTINE Get_Cloud_Opt_MW_AD

END MODULE CRTM_CloudScatter
