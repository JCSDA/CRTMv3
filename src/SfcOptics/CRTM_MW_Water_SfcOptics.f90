!
! CRTM_MW_Water_SfcOptics
!
! Module to compute the surface optical properties for WATER surfaces at
! microwave frequencies required for determining the WATER surface
! contribution to the radiative transfer.
!
! This module is provided to allow developers to "wrap" their existing
! codes inside the provided functions to simplify integration into
! the main CRTM_SfcOptics module.
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 25-Jun-2005
!                       paul.vandelst@noaa.gov
!

MODULE CRTM_MW_Water_SfcOptics

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE Type_Kinds,               ONLY: fp
  USE Message_Handler,          ONLY: SUCCESS, FAILURE, WARNING, Display_Message
  USE CRTM_Parameters,          ONLY: SET, NOT_SET, &
                                      ZERO, ONE, &
                                      MAX_N_ANGLES, &
                                      N_STOKES => MAX_N_STOKES
  USE CRTM_SpcCoeff,            ONLY: SC
  USE CRTM_Surface_Define,      ONLY: CRTM_Surface_type
  USE CRTM_GeometryInfo_Define, ONLY: CRTM_GeometryInfo_type, &
                                      CRTM_GeometryInfo_GetValue
  USE CRTM_SfcOptics_Define,    ONLY: CRTM_SfcOptics_type
  USE CRTM_LowFrequency_MWSSEM, ONLY: LF_MWSSEM_type => iVar_type, &
                                      LowFrequency_MWSSEM, &
                                      LowFrequency_MWSSEM_TL, &
                                      LowFrequency_MWSSEM_AD
  USE CRTM_Fastem1,             ONLY: Fastem1
  USE CRTM_FastemX,             ONLY: FastemX_type => iVar_type, &
                                      Compute_FastemX,   &
                                      Compute_FastemX_TL,&
                                      Compute_FastemX_AD
  USE CRTM_MWwaterCoeff       , ONLY: MWwaterC
  USE CRTM_PARMIO,              ONLY: PARMIO_type => iVar_type, &
                                      Compute_PARMIO
  USE PARMIO_LUT_Interpolation, ONLY: PARMIO_LUT_Clamped_Axes
  USE CRTM_PARMIO_TL,           ONLY: Compute_PARMIO_TL
  USE CRTM_PARMIO_AD,           ONLY: Compute_PARMIO_AD
  USE CRTM_PARMIOCoeff,         ONLY: PARMIOC, CRTM_PARMIOCoeff_IsLoaded, &
                                      CRTM_PARMIOCoeff_Covers_Frequency
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
  PUBLIC :: Compute_MW_Water_SfcOptics
  PUBLIC :: Compute_MW_Water_SfcOptics_TL
  PUBLIC :: Compute_MW_Water_SfcOptics_AD
  ! Dispatch threshold (CRTM_LifeCycle reports the FASTEM fallback against it)
  PUBLIC :: PARMIO_FREQ_THRESHOLD
  ! Runtime control and the single predicate answering "will PARMIO be used
  ! at this frequency?". Everything that needs to know asks this rather than
  ! re-deriving it, so the policy and the table-coverage rule cannot drift
  ! apart between call sites.
  PUBLIC :: PARMIO_Set_Freq_Threshold
  PUBLIC :: PARMIO_Get_Freq_Threshold
  PUBLIC :: PARMIO_Is_Active_At


  ! -----------------
  ! Module parameters
  ! -----------------
  ! Low frequency model threshold
  REAL(fp), PARAMETER :: LOW_F_THRESHOLD = 20.0_fp ! GHz
  ! Finite-difference step (K) for the Fastem1 emissivity SST derivative. Fastem1
  ! returns only wind-speed derivatives, so d(emissivity)/d(Water_Temperature) is
  ! obtained by a central difference around the forward call (see below).
  REAL(fp), PARAMETER :: FASTEM1_DTS = 0.1_fp
  ! PARMIO LUT is the surface-emissivity backend at and above this frequency
  ! when the LUT has been loaded. Below this threshold the FASTEM/Stogryn
  ! legacy path is used.
  ! Set from obs-space validation against ATMS-NPP (see
  !   PARMIO vs FASTEM6 — ATMS-NPP obs-space validation report, 2026-05-13):
  ! FASTEM6 is tuned and competitive against real obs through the entire ATMS
  ! band (max 183.31 GHz); PARMIO's physical-reference advantage only shows
  ! cleanly where FASTEM6 extrapolates beyond its tuning band (e.g. 325 GHz
  ! synthetic-RT wind-roughness sign flip).
  REAL(fp), PARAMETER :: PARMIO_FREQ_THRESHOLD = 200.0_fp ! GHz

  ! The policy threshold above is only half the question. Being loaded and
  ! being above the threshold does not mean the table has data at a given
  ! frequency: the coefficient groups are gridded separately either side of
  ! the permittivity switch and their grids do not meet it. In the shipped
  ! production LUT sss_nominal_m stops at 183.31 GHz and sss_nominal_h starts
  ! at 229 GHz, so 183.31 to 229 selects a group with nothing in it and the
  ! interpolator clamps to the nearest grid edge without saying so. A
  ! 204.78 GHz channel was being evaluated at 229 GHz.
  !
  ! So the dispatch asks both questions: is PARMIO allowed here (policy), and
  ! does it have data here (the table). With the shipped LUT that makes the
  ! effective floor 229 GHz rather than 200, and it will become 200 on its own
  ! if the 183 to 229 band is ever filled in, with no constant to update.
  !
  ! The effective threshold is a variable rather than a parameter so users can
  ! test PARMIO outside its default band. It is written once at
  ! initialization, before any forward call, and read-only thereafter, so it
  ! carries the same threading characteristics as the coefficient data.
  REAL(fp), SAVE :: PARMIO_Threshold_InUse = PARMIO_FREQ_THRESHOLD

  ! Latch for the edge-clamp report, so it is stated once rather than on every
  ! angle of every channel of every profile. Re-armed whenever the dispatch
  ! threshold is set, which CRTM_Init always does, so each run reports afresh.
  LOGICAL, SAVE :: PARMIO_Clamp_Warn_Pending = .TRUE.


  ! --------------------------------------
  ! Structure definition to hold forward
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
  TYPE :: iVar_type
    PRIVATE
    ! FastemX model internal variable structure
    TYPE(FastemX_type), DIMENSION(MAX_N_ANGLES)   :: FastemX_Var
    ! Low frequency model internal variable structure
    TYPE(LF_MWSSEM_type), DIMENSION(MAX_N_ANGLES) :: LF_MWSSEM_Var
    ! PARMIO model internal variable structure
    TYPE(PARMIO_type), DIMENSION(MAX_N_ANGLES)    :: PARMIO_Var
    ! Fastem outputs
    REAL(fp), DIMENSION(MAX_N_ANGLES) :: dEH_dTs        = ZERO
    REAL(fp), DIMENSION(MAX_N_ANGLES) :: dEH_dWindSpeed = ZERO
    REAL(fp), DIMENSION(MAX_N_ANGLES) :: dEV_dTs        = ZERO
    REAL(fp), DIMENSION(MAX_N_ANGLES) :: dEV_dWindSpeed = ZERO
  END TYPE iVar_type

CONTAINS


!--------------------------------------------------------------------------------
!
! PARMIO dispatch control.
!
!   PARMIO_Set_Freq_Threshold  set the frequency at and above which PARMIO is
!                              permitted, for users who want to exercise it
!                              outside its default band. Call before any
!                              forward call; CRTM_Init does this from its
!                              optional PARMIO_Freq_Threshold argument.
!   PARMIO_Get_Freq_Threshold  report the value in force.
!   PARMIO_Is_Active_At        the one predicate answering whether PARMIO will
!                              actually serve a given frequency.
!
! Lowering the threshold does not let a caller into a band the table does not
! cover: coverage is a hard requirement, not a preference, because the
! alternative is a silently clamped result from the wrong frequency. If the
! band is genuinely wanted, the fix is LUT data, not a looser gate.
!
!--------------------------------------------------------------------------------

  SUBROUTINE PARMIO_Set_Freq_Threshold( Frequency_GHz )
    REAL(fp), INTENT(IN) :: Frequency_GHz
    PARMIO_Threshold_InUse    = Frequency_GHz
    PARMIO_Clamp_Warn_Pending = .TRUE.
  END SUBROUTINE PARMIO_Set_Freq_Threshold


  PURE FUNCTION PARMIO_Get_Freq_Threshold() RESULT( Frequency_GHz )
    REAL(fp) :: Frequency_GHz
    Frequency_GHz = PARMIO_Threshold_InUse
  END FUNCTION PARMIO_Get_Freq_Threshold


  PURE FUNCTION PARMIO_Is_Active_At( Frequency ) RESULT( Active )
    REAL(fp), INTENT(IN) :: Frequency
    LOGICAL :: Active
    Active = CRTM_PARMIOCoeff_IsLoaded()
    IF ( Active ) Active = ( Frequency >= PARMIO_Threshold_InUse )
    IF ( Active ) Active = CRTM_PARMIOCoeff_Covers_Frequency( Frequency )
  END FUNCTION PARMIO_Is_Active_At



!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       Compute_MW_Water_SfcOptics
!
! PURPOSE:
!       Function to compute the surface emissivity and reflectivity at microwave
!       frequencies over a water surface.
!
!       This function is a wrapper for third party code.
!
! CALLING SEQUENCE:
!       Error_Status = Compute_MW_Water_SfcOptics( Surface     , &  ! Input
!                                                  GeometryInfo, &  ! Input
!                                                  SensorIndex , &  ! Input
!                                                  ChannelIndex, &  ! Input
!                                                  SfcOptics   , &  ! Output
!                                                  iVar          )  ! Internal variable output
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
!       iVar:            Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        The contents of this structure are NOT accessible
!                        outside of the CRTM_MW_Water_SfcOptics module.
!                        UNITS:      N/A
!                        TYPE:       iVar_type
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
! COMMENTS:
!       Note the INTENT on the output SfcOptics argument is IN OUT rather
!       than just OUT as it is assumed to contain some data upon input.
!
!:sdoc-:
!----------------------------------------------------------------------------------

  FUNCTION Compute_MW_Water_SfcOptics( &
    Surface     , &  ! Input
    GeometryInfo, &  ! Input
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    SfcOptics   , &  ! Output
    iVar        ) &  ! Internal variable output
  RESULT( err_stat )
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
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Compute_MW_Water_SfcOptics'
    CHARACTER(64) :: Clamped_Axes
    ! Local variables
    INTEGER  :: i, j
    REAL(fp) :: Frequency
    REAL(fp) :: Source_Azimuth_Angle, Sensor_Azimuth_Angle
    REAL(fp) :: Reflectivity(N_STOKES)
    ! Fastem1 SST-derivative finite-difference scratch (V=1, H=2)
    REAL(fp) :: emis_pTs(2), emis_mTs(2), dwind_h, dwind_v


    ! Set up
    err_stat = SUCCESS
    SfcOptics%Reflectivity = ZERO
    ! ...Retrieve data from structures
    Frequency = SC(SensorIndex)%Frequency(ChannelIndex)
    CALL CRTM_GeometryInfo_GetValue( &
           GeometryInfo, &
           Source_Azimuth_Angle = Source_Azimuth_Angle, &
           Sensor_Azimuth_Angle = Sensor_Azimuth_Angle  )


    ! ------------------------------------------------------------------
    ! Polarimetric azimuth convention (authoritative definition)
    ! ------------------------------------------------------------------
    ! Every microwave water backend below takes its relative azimuth from
    ! the single expression
    !
    !     phi = Surface%Wind_Direction - Sensor_Azimuth_Angle   [degrees]
    !
    ! and each applies it as phi_radians = phi * DEGREES_TO_RADIANS with no
    ! further reflection or offset. The two terms are defined by CRTM as
    !
    !   Wind_Direction       the direction the wind blows TOWARD, clockwise
    !                        from North. Zero is a wind blowing toward the
    !                        north, i.e. a southerly. This is the opposite
    !                        of the meteorological convention
    !                        (CRTM_Surface_Define.f90, DEFAULT_WIND_DIRECTION).
    !   Sensor_Azimuth_Angle the azimuth of the horizontal projection of the
    !                        line from the satellite to the FOV, clockwise
    !                        from North (CRTM_Geometry_Define.f90:312-316).
    !
    ! so phi = 0 means the wind blows toward the same compass azimuth as the
    ! satellite-to-FOV horizontal projection.
    !
    ! The azimuthal emissivity is expanded as (Liu et al., FASTEM-4
    ! validation, NWPSAF-MO-VS-045, equations 2a-2d)
    !
    !     e_V = ... + SUM_m c_m cos(m phi)      e_U   = SUM_m e_m sin(m phi)
    !     e_H = ... + SUM_m d_m cos(m phi)      e_V4  = SUM_m g_m sin(m phi)
    !
    ! Cosine for V and H, sine for the third and fourth Stokes components.
    ! All three backends implement exactly this: Azimuth_Emissivity_Module
    ! (FASTEM4/5), Azimuth_Emissivity_F6_Module (FASTEM6) and
    ! PARMIO_Azimuth_Module. The third Stokes component follows the standard
    ! radiometric definition U = T(+45) - T(-45), as used by WindSat, whose
    ! measurements the FASTEM azimuth coefficients were fitted to, and by
    ! RTTOV.
    !
    ! Consequences worth knowing:
    !   * V and H are EVEN in phi and U and V4 are ODD. A sign error in the
    !     azimuth convention is therefore invisible in I and Q and shows up
    !     only in U and V4. test_VectorRT_SurfaceFrame pins that parity.
    !   * FASTEM6, the CRTM default, parameterises V and H only and returns
    !     the third and fourth Stokes components as identically zero. A
    !     polarimetric run has a real surface U and V4 only on FASTEM4 or
    !     PARMIO. Note that only FASTEM4 and FASTEM6 can actually be loaded;
    !     Azimuth_Emissivity_Module serves FASTEM4 and FASTEM5, but FASTEM5
    !     is not a selectable scheme and asking for it is a hard error.
    !
    ! See docs/design/polarimetric_conventions.md for the full statement,
    ! the literature basis, and the one item still open (whether this phi
    ! origin matches the one the FASTEM coefficients were fitted under).
    ! ------------------------------------------------------------------
    !
    ! Compute the surface optical parameters
    ! PARMIO dispatch is gated by frequency: at and above
    ! PARMIO_FREQ_THRESHOLD (and provided the LUT was loaded at CRTM_Init
    ! time) PARMIO is used as the MW-water emissivity backend. Below the
    ! threshold the FASTEM/Stogryn legacy path is used. If no LUT was
    ! loaded the path is byte-identical to a pre-PARMIO build at every
    ! frequency.
    IF( PARMIO_Is_Active_At( Frequency ) ) THEN

      ! PARMIO_MWSSEM (LUT-driven, replaces FASTEM at runtime)
      SfcOptics%Azimuth_Angle = Surface%Wind_Direction - Sensor_Azimuth_Angle
      DO i = 1, SfcOptics%n_Angles
        CALL Compute_PARMIO( &
               PARMIOC                                , &  ! Input PARMIO LUT coefficients
               Frequency                              , &  ! Input
               SfcOptics%n_Angles                     , &  ! Input
               SfcOptics%Angle(i)                     , &  ! Input
               Surface%Water_Temperature              , &  ! Input
               Surface%Salinity                       , &  ! Input
               Surface%Wind_Speed                     , &  ! Input
               iVar%PARMIO_Var(i)                     , &  ! Internal variable output
               SfcOptics%Emissivity(i,:)              , &  ! Output
               Reflectivity                           , &  ! Output
               Azimuth_Angle = SfcOptics%Azimuth_Angle, &  ! Optional input
               Transmittance = SfcOptics%Transmittance  )  ! Optional input
        DO j = 1, N_STOKES
          SfcOptics%Reflectivity(i,j,i,j) = Reflectivity(j)
        END DO
        ! Report an edge-clamped lookup once. The interpolator pins an
        ! out-of-range query to the nearest grid node and returns a confident
        ! number computed somewhere other than where it was asked, which is
        ! defensible as a fallback and not as a silent one. Frequency cannot
        ! clamp here because PARMIO_Is_Active_At already required coverage,
        ! but the state axes can and do: the table spans zenith 0 to 65 deg,
        ! wind 1 to 25 m/s and SST -2 to 30 C, all of which real scenes exceed.
        ! Latched, because this sits in the per-angle loop of every channel of
        ! every profile.
        IF ( PARMIO_Clamp_Warn_Pending ) THEN
          Clamped_Axes = PARMIO_LUT_Clamped_Axes( iVar%PARMIO_Var(i)%LUT_Var )
          IF ( LEN_TRIM(Clamped_Axes) > 0 ) THEN
            PARMIO_Clamp_Warn_Pending = .FALSE.
            CALL Display_Message( ROUTINE_NAME, &
              'PARMIO lookup clamped to the table edge on: '//TRIM(Clamped_Axes)//&
              '. Results there are evaluated at the nearest grid node, not at '//&
              'the requested value.', WARNING )
          END IF
        END IF
      END DO

    ELSE IF ( SfcOptics%Use_New_MWSSEM ) THEN

      ! FastemX model
      SfcOptics%Azimuth_Angle = Surface%Wind_Direction - Sensor_Azimuth_Angle
      DO i = 1, SfcOptics%n_Angles
        CALL Compute_FastemX( &
               MWwaterC                               , &  ! Input model coefficients
               Frequency                              , &  ! Input
               SfcOptics%n_Angles                     , &  ! Input
               SfcOptics%Angle(i)                     , &  ! Input
               Surface%Water_Temperature              , &  ! Input
               Surface%Salinity                       , &  ! Input
               Surface%Wind_Speed                     , &  ! Input
               iVar%FastemX_Var(i)                    , &  ! Internal variable output
               SfcOptics%Emissivity(i,:)              , &  ! Output
               Reflectivity                           , &  ! Output
               Azimuth_Angle = SfcOptics%Azimuth_Angle, &  ! Optional input
               Transmittance = SfcOptics%Transmittance  )  ! Optional input
        DO j = 1, N_STOKES
          SfcOptics%Reflectivity(i,j,i,j) = Reflectivity(j)
        END DO
      END DO

    ELSE

      ! Low frequency model coupled with Fastem1
      IF( Frequency < LOW_F_THRESHOLD ) THEN
        ! Call the low frequency model
        DO i = 1, SfcOptics%n_Angles
          CALL LowFrequency_MWSSEM( &
                 Frequency                , &  ! Input
                 SfcOptics%Angle(i)       , &  ! Input
                 Surface%Water_Temperature, &  ! Input
                 Surface%Salinity         , &  ! Input
                 Surface%Wind_Speed       , &  ! Input
                 SfcOptics%Emissivity(i,:), &  ! Output
                 iVar%LF_MWSSEM_Var(i)      )  ! Internal variable output
          SfcOptics%Reflectivity(i,1,i,1) = ONE-SfcOptics%Emissivity(i,1)
          SfcOptics%Reflectivity(i,2,i,2) = ONE-SfcOptics%Emissivity(i,2)
        END DO
      ELSE
        ! Call Fastem1
        DO i = 1, SfcOptics%n_Angles
          CALL Fastem1( Frequency                , & ! Input
                        SfcOptics%Angle(i)       , & ! Input
                        Surface%Water_Temperature, & ! Input
                        Surface%Wind_Speed       , & ! Input
                        SfcOptics%Emissivity(i,:), & ! Output
                        iVar%dEH_dWindSpeed(i)   , & ! Output
                        iVar%dEV_dWindSpeed(i)     ) ! Output
          ! Fastem1 returns no SST derivative; obtain d(emissivity)/d(Water_Temperature)
          ! by a central finite difference around the forward call so the TL/AD SST
          ! Jacobian (iVar%dE?_dTs, read below) is not silently zero.
          CALL Fastem1( Frequency, SfcOptics%Angle(i), Surface%Water_Temperature+FASTEM1_DTS, &
                        Surface%Wind_Speed, emis_pTs, dwind_h, dwind_v )
          CALL Fastem1( Frequency, SfcOptics%Angle(i), Surface%Water_Temperature-FASTEM1_DTS, &
                        Surface%Wind_Speed, emis_mTs, dwind_h, dwind_v )
          iVar%dEV_dTs(i) = (emis_pTs(1) - emis_mTs(1))/(2.0_fp*FASTEM1_DTS)  ! V (index 1)
          iVar%dEH_dTs(i) = (emis_pTs(2) - emis_mTs(2))/(2.0_fp*FASTEM1_DTS)  ! H (index 2)
          SfcOptics%Reflectivity(i,1,i,1) = ONE-SfcOptics%Emissivity(i,1)
          SfcOptics%Reflectivity(i,2,i,2) = ONE-SfcOptics%Emissivity(i,2)
        END DO
      END IF

    END IF

  END FUNCTION Compute_MW_Water_SfcOptics


!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       Compute_MW_Water_SfcOptics_TL
!
! PURPOSE:
!       Function to compute the tangent-linear surface emissivity and
!       reflectivity at microwave frequencies over a water surface.
!
!       This function is a wrapper for third party code.
!
! CALLING SEQUENCE:
!       Error_Status = Compute_MW_Water_SfcOptics_TL( Surface     , &  ! Input
!                                                     SfcOptics   , &  ! Input
!                                                     Surface_TL  , &  ! Input
!                                                     GeometryInfo, &  ! Input
!                                                     SensorIndex , &  ! Input
!                                                     ChannelIndex, &  ! Output
!                                                     SfcOptics_TL, &  ! Output
!                                                     iVar          )  ! Internal variable input
!
! INPUTS:
!       Surface:         CRTM_Surface structure containing the surface state
!                        data.
!                        UNITS:      N/A
!                        TYPE:       CRTM_Surface_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Surface_TL:      CRTM_Surface structure containing the tangent-linear
!                        surface state data.
!                        UNITS:      N/A
!                        TYPE:       CRTM_Surface_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       SfcOptics:       CRTM_SfcOptics structure containing the surface
!                        optical properties required for the radiative
!                        transfer calculation.
!                        UNITS:      N/A
!                        TYPE:       CRTM_SfcOptics_type
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
!       iVar:            Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        The contents of this structure are NOT accessible
!                        outside of the CRTM_MW_Water_SfcOptics module.
!                        UNITS:      N/A
!                        TYPE:       iVar_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       SfcOptics_TL:    CRTM_SfcOptics structure containing the tangent-linear
!                        surface optical properties required for the tangent-
!                        linear radiative transfer calculation.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_SfcOptics_type)
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
!       upon input. To prevent memory leaks, the IN OUT INTENT is a must.
!
!:sdoc-:
!----------------------------------------------------------------------------------

  FUNCTION Compute_MW_Water_SfcOptics_TL( &
    SfcOptics   , &  ! Input
    Surface_TL  , &  ! Input
    GeometryInfo, &  ! Input
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    SfcOptics_TL, &  ! Output
    iVar        ) &  ! Internal variable input
  RESULT( err_stat )
    ! Arguments
    TYPE(CRTM_Surface_type),      INTENT(IN)     :: Surface_TL
    TYPE(CRTM_SfcOptics_type),    INTENT(IN)     :: SfcOptics
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)     :: GeometryInfo
    INTEGER,                      INTENT(IN)     :: SensorIndex
    INTEGER,                      INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_SfcOptics_type),    INTENT(IN OUT) :: SfcOptics_TL
    TYPE(iVar_type),              INTENT(IN)     :: iVar
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Compute_MW_Water_SfcOptics_TL'
    ! Local variables
    INTEGER :: i, j
    REAL(fp) :: Frequency
    REAL(fp) :: Source_Azimuth_Angle, Sensor_Azimuth_Angle
    REAL(fp) :: Reflectivity_TL(N_STOKES)


    ! Set up
    err_stat = SUCCESS
    SfcOptics_TL%Reflectivity = ZERO
    ! ...Retrieve data from structures
    Frequency = SC(SensorIndex)%Frequency(ChannelIndex)
    CALL CRTM_GeometryInfo_GetValue( &
           GeometryInfo, &
           Source_Azimuth_Angle = Source_Azimuth_Angle, &
           Sensor_Azimuth_Angle = Sensor_Azimuth_Angle  )


    ! Compute the tangent-linear surface optical parameters
    ! Dispatch matches the Forward path: PARMIO at and above
    ! PARMIO_FREQ_THRESHOLD when the LUT is loaded.
    IF( PARMIO_Is_Active_At( Frequency ) ) THEN

      ! PARMIO_MWSSEM (LUT-driven)
      DO i = 1, SfcOptics%n_Angles
        CALL Compute_PARMIO_TL( &
               PARMIOC                                     , &  ! Input PARMIO LUT coefficients
               Surface_TL%Water_Temperature                , &  ! TL Input
               Surface_TL%Salinity                         , &  ! TL Input
               Surface_TL%Wind_Speed                       , &  ! TL Input
               iVar%PARMIO_Var(i)                          , &  ! Internal variable input
               SfcOptics_TL%Emissivity(i,:)                , &  ! TL Output
               Reflectivity_TL                             , &  ! TL Output
               Azimuth_Angle_TL = Surface_TL%Wind_Direction, &  ! Optional TL input
               Transmittance_TL = SfcOptics_TL%Transmittance )  ! Optional TL input
        DO j = 1, N_STOKES
          SfcOptics_TL%Reflectivity(i,j,i,j) = Reflectivity_TL(j)
        END DO
      END DO

    ELSE IF( SfcOptics%Use_New_MWSSEM ) THEN

      ! FastemX model
      DO i = 1, SfcOptics%n_Angles
        CALL Compute_FastemX_TL( &
               MWwaterC                                    , &  ! Input model coefficients
               Surface_TL%Water_Temperature                , &  ! TL Input
               Surface_TL%Salinity                         , &  ! TL Input
               Surface_TL%Wind_Speed                       , &  ! TL Input
               iVar%FastemX_Var(i)                         , &  ! Internal variable input
               SfcOptics_TL%Emissivity(i,:)                , &  ! TL Output
               Reflectivity_TL                             , &  ! TL Output
               Azimuth_Angle_TL = Surface_TL%Wind_Direction, &  ! Optional TL input
               Transmittance_TL = SfcOptics_TL%Transmittance )  ! Optional TL input
        DO j = 1, N_STOKES
          !we probably need further low-level check and modifications
          !SfcOptics_TL%Reflectivity(i,j,i,j) = -Reflectivity_TL(j)
          SfcOptics_TL%Reflectivity(i,j,i,j) = Reflectivity_TL(j)
        END DO
      END DO

    ELSE

      ! Low frequency model coupled with Fastem1
      IF( Frequency < LOW_F_THRESHOLD ) THEN
        ! Call the low frequency model
        DO i = 1, SfcOptics%n_Angles
          CALL LowFrequency_MWSSEM_TL( &
                 Surface_TL%Water_Temperature, &  ! TL  Input
                 Surface_TL%Salinity         , &  ! TL  Input
                 Surface_TL%Wind_Speed       , &  ! TL  Input
                 SfcOptics_TL%Emissivity(i,:), &  ! TL  Output
                 iVar%LF_MWSSEM_Var(i)         )  ! Internal variable input
          SfcOptics_TL%Reflectivity(i,1,i,1) = -SfcOptics_TL%Emissivity(i,1)
          SfcOptics_TL%Reflectivity(i,2,i,2) = -SfcOptics_TL%Emissivity(i,2)
        END DO
      ELSE
        ! Call Fastem1
        DO i = 1, SfcOptics%n_Angles
          SfcOptics_TL%Emissivity(i,2) = (iVar%dEH_dTs(i)*Surface_TL%Water_Temperature) + &
                                         (iVar%dEH_dWindSpeed(i)*Surface_TL%Wind_Speed)
          SfcOptics_TL%Emissivity(i,1) = (iVar%dEV_dTs(i)*Surface_TL%Water_Temperature) + &
                                         (iVar%dEV_dWindSpeed(i)*Surface_TL%Wind_Speed)
          SfcOptics_TL%Reflectivity(i,1,i,1) = -SfcOptics_TL%Emissivity(i,1)
          SfcOptics_TL%Reflectivity(i,2,i,2) = -SfcOptics_TL%Emissivity(i,2)
        END DO
      END IF
    END IF

  END FUNCTION Compute_MW_Water_SfcOptics_TL


!----------------------------------------------------------------------------------
!
! NAME:
!       Compute_MW_Water_SfcOptics_AD
!
! PURPOSE:
!       Function to compute the adjoint surface emissivity and
!       reflectivity at microwave frequencies over a water surface.
!
!       This function is a wrapper for third party code.
!
! CALLING SEQUENCE:
!       Error_Status = Compute_MW_Water_SfcOptics_AD( Surface     , &  ! Input
!                                                     SfcOptics   , &  ! Input
!                                                     SfcOptics_AD, &  ! Input
!                                                     GeometryInfo, &  ! Input
!                                                     SensorIndex , &  ! Input
!                                                     ChannelIndex, &  ! Output
!                                                     Surface_AD  , &  ! Output
!                                                     iVar          )  ! Internal variable input
!
! INPUT ARGUMENTS:
!       Surface:         CRTM_Surface structure containing the surface state
!                        data.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_Surface_type)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       SfcOptics:       CRTM_SfcOptics structure containing the surface
!                        optical properties required for the radiative
!                        transfer calculation.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_SfcOptics_type)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       SfcOptics_AD:    CRTM_SfcOptics structure containing the adjoint
!                        surface optical properties required for the adjoint
!                        radiative transfer calculation.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_SfcOptics_type)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
!       GeometryInfo:    CRTM_GeometryInfo structure containing the
!                        view geometry information.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_GeometryInfo_type)
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
!       iVar:            Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        The contents of this structure are NOT accessible
!                        outside of the CRTM_MW_Water_SfcOptics module.
!                        UNITS:      N/A
!                        TYPE:       iVar_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       Surface_AD:      CRTM_Surface structure containing the adjoint
!                        surface state data.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_Surface_type)
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
!       Note the INTENT on the input SfcOptics_AD argument is IN OUT rather
!       than just OUT. This is necessary because components of this argument
!       may need to be zeroed out upon output.
!
!       Note the INTENT on the output Surface_AD argument is IN OUT rather
!       than just OUT. This is necessary because the argument may be defined
!       upon input. To prevent memory leaks, the IN OUT INTENT is a must.
!
!----------------------------------------------------------------------------------

  FUNCTION Compute_MW_Water_SfcOptics_AD( &
    SfcOptics   , &  ! Input
    SfcOptics_AD, &  ! Input
    GeometryInfo, &  ! Input
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    Surface_AD  , &  ! Output
    iVar        ) &  ! Internal variable input
  RESULT( err_stat )
    ! Arguments
    TYPE(CRTM_SfcOptics_type),    INTENT(IN)     :: SfcOptics
    TYPE(CRTM_SfcOptics_type),    INTENT(IN OUT) :: SfcOptics_AD
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)     :: GeometryInfo
    INTEGER,                      INTENT(IN)     :: SensorIndex
    INTEGER,                      INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_Surface_type),      INTENT(IN OUT) :: Surface_AD
    TYPE(iVar_type),              INTENT(IN)     :: iVar
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Compute_MW_Water_SfcOptics_AD'
    ! Local variables
    INTEGER :: i, j
    REAL(fp) :: Frequency
    REAL(fp) :: Source_Azimuth_Angle, Sensor_Azimuth_Angle
    REAL(fp) :: Reflectivity_AD(N_STOKES)
    REAL(fp) :: Azimuth_Angle_AD


    ! Set up
    err_stat = SUCCESS
    ! ...Retrieve data from structures
    Frequency = SC(SensorIndex)%Frequency(ChannelIndex)
    CALL CRTM_GeometryInfo_GetValue( &
           GeometryInfo, &
           Source_Azimuth_Angle = Source_Azimuth_Angle, &
           Sensor_Azimuth_Angle = Sensor_Azimuth_Angle  )


    ! Compute the adjoint surface optical parameters
    ! Dispatch matches the Forward path: PARMIO at and above
    ! PARMIO_FREQ_THRESHOLD when the LUT is loaded.
    IF( PARMIO_Is_Active_At( Frequency ) ) THEN

      ! PARMIO_MWSSEM (LUT-driven)
      Azimuth_Angle_AD = ZERO
      DO i = 1, SfcOptics%n_Angles
        DO j = 1, N_STOKES
          Reflectivity_AD(j) = SfcOptics_AD%Reflectivity(i,j,i,j)
        END DO
        CALL Compute_PARMIO_AD( &
               PARMIOC                                      , &  ! Input PARMIO LUT coefficients
               SfcOptics_AD%Emissivity(i,:)                 , &  ! AD Input
               Reflectivity_AD                              , &  ! AD Input
               iVar%PARMIO_Var(i)                           , &  ! Internal variable input
               Surface_AD%Water_Temperature                 , &  ! AD Output
               Surface_AD%Salinity                          , &  ! AD Output
               Surface_AD%Wind_Speed                        , &  ! AD Output
               Azimuth_Angle_AD = Azimuth_Angle_AD          , &  ! Optional AD Output
               Transmittance_AD = SfcOptics_AD%Transmittance  )  ! Optional AD Output
      END DO
      Surface_AD%Wind_Direction = Surface_AD%Wind_Direction + Azimuth_Angle_AD

    ELSE IF( SfcOptics%Use_New_MWSSEM ) THEN

      ! FastemX model
      Azimuth_Angle_AD = ZERO
      DO i = 1, SfcOptics%n_Angles
        DO j = 1, N_STOKES
          Reflectivity_AD(j) = SfcOptics_AD%Reflectivity(i,j,i,j)
        END DO
        CALL Compute_FastemX_AD( &
               MWwaterC                                     , &  ! Input model coefficients
               SfcOptics_AD%Emissivity(i,:)                 , &  ! AD Input
               Reflectivity_ad                              , &  ! AD Input
               iVar%FastemX_Var(i)                          , &  ! Internal variable input
               Surface_AD%Water_Temperature                 , &  ! AD Output
               Surface_AD%Salinity                          , &  ! AD Output
               Surface_AD%Wind_Speed                        , &  ! AD Output
               Azimuth_Angle_AD = Azimuth_Angle_AD          , &  ! Optional AD Output
               Transmittance_AD = SfcOptics_AD%Transmittance  )  ! Optional AD Output
      END DO
      Surface_AD%Wind_Direction = Surface_AD%Wind_Direction + Azimuth_Angle_AD

    ELSE

      ! Low frequency model coupled with Fastem1
      IF( Frequency < LOW_F_THRESHOLD ) THEN
        ! Call the low frequency model
        DO i = 1, SfcOptics%n_Angles
          SfcOptics_AD%Emissivity(i,1) = SfcOptics_AD%Emissivity(i,1)-SfcOptics_AD%Reflectivity(i,1,i,1)
          SfcOptics_AD%Emissivity(i,2) = SfcOptics_AD%Emissivity(i,2)-SfcOptics_AD%Reflectivity(i,2,i,2)
          CALL LowFrequency_MWSSEM_AD( &
                 SfcOptics_AD%Emissivity(i,:), &  ! AD  Input
                 Surface_AD%Water_Temperature, &  ! AD  Output
                 Surface_AD%Salinity         , &  ! AD  Output
                 Surface_AD%Wind_Speed       , &  ! AD  Output
                 iVar%LF_MWSSEM_Var(i)         )  ! Internal variable input
        END DO
      ELSE
        ! Call Fastem1
        DO i = SfcOptics%n_Angles, 1, -1
          DO j = 1, 2
            SfcOptics_AD%Emissivity(i,j) = SfcOptics_AD%Emissivity(i,j) - &
                                           SfcOptics_AD%Reflectivity(i,j,i,j)
            SfcOptics_AD%Reflectivity(i,j,i,j) = ZERO
          END DO
          ! Vertical polarisation component
          Surface_AD%Water_Temperature  = Surface_AD%Water_Temperature + &
                                          (iVar%dEV_dTs(i)*SfcOptics_AD%Emissivity(i,1))
          Surface_AD%Wind_Speed         = Surface_AD%Wind_Speed + &
                                          (iVar%dEV_dWindSpeed(i)*SfcOptics_AD%Emissivity(i,1))
          SfcOptics_AD%Emissivity(i,1)  = ZERO
          ! Horizontal polarization component
          Surface_AD%Water_Temperature  = Surface_AD%Water_Temperature + &
                                          (iVar%dEH_dTs(i)*SfcOptics_AD%Emissivity(i,2))
          Surface_AD%Wind_Speed         = Surface_AD%Wind_Speed + &
                                          (iVar%dEH_dWindSpeed(i)*SfcOptics_AD%Emissivity(i,2))
          SfcOptics_AD%Emissivity(i,2)  = ZERO
        END DO
      END IF
    END IF

    SfcOptics_AD%Reflectivity = ZERO

  END FUNCTION Compute_MW_Water_SfcOptics_AD

END MODULE CRTM_MW_Water_SfcOptics
