!
! CRTM_VISsnowRF.f90
!
! Module containing functions to invoke the CRTM Visible/Near-IR
! Snow Reflectance Model (SNICAR-based LUT).
!
! The reflectance LUT is 5-dimensional:
!   Reflectance(Angle, Frequency, Grain_Size, Depth, Density)
!
! Three model functions are provided following the standard CRTM
! forward / tangent-linear / adjoint (FWD/TL/AD) pattern:
!
!   CRTM_Compute_VISsnowRefl     – forward model
!   CRTM_Compute_VISsnowRefl_TL  – tangent-linear model
!   CRTM_Compute_VISsnowRefl_AD  – adjoint model
!
! The TL and AD functions must be called *after* the forward function
! because they reuse the internal variable structure (iVar) populated
! by the forward call.
!
!
! CREATION HISTORY:
!       Written by:   Cheng Dang, Jun-2026
!                     dangch@ucar.edu
!

MODULE CRTM_VISsnowRF

  ! -----------------
  ! Environment setup
  ! -----------------
  USE Type_Kinds,          ONLY: fp
  USE Message_Handler,     ONLY: SUCCESS, FAILURE, Display_Message
  USE CRTM_Parameters,     ONLY: ZERO, ONE, DEGREES_TO_RADIANS
  USE CRTM_Interpolation,  ONLY: NPTS,          &
                                 LPoly,         &
                                 LPoly_type,    &
                                 Clear_LPoly,   &
                                 Find_Index,    &
                                 Interp_1D,     &
                                 Interp_4D,     &
                                 Interp_1D_TL,  &
                                 Interp_4D_TL,  &
                                 LPoly_TL,      &
                                 Interp_1D_AD,  &
                                 Interp_4D_AD,  &
                                 LPoly_AD
  USE VISsnowCoeff_Define, ONLY: VISsnowCoeff_type
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  ! Derived type
  PUBLIC :: iVar_type
  ! Procedures
  PUBLIC :: CRTM_Compute_VISsnowRF
  PUBLIC :: CRTM_Compute_VISsnowRF_TL
  PUBLIC :: CRTM_Compute_VISsnowRF_AD


  ! -----------------
  ! Module parameters
  ! -----------------
  INTEGER, PARAMETER :: ML = 256


  ! -------------------------------------------------------
  ! Einterp_type
  !
  ! Internal interpolation variable structure.  Holds all
  ! Lagrange polynomial objects, LUT sub-array indices, and
  ! scalar inputs needed by the TL and AD calls.
  !
  ! Dimension key (matching VISsnowCoeff_type):
  !   I  – Angle
  !   L  – Frequency
  !   G  – Grain_Size
  !   T  – Depth
  !   J  – Density
  ! -------------------------------------------------------
  TYPE :: Einterp_type
    ! Scalar dimensions
    INTEGER :: n_Angles = 0
    INTEGER :: n_Pts    = 0
    ! Allocation flag
    LOGICAL :: Is_Allocated = .FALSE.
    ! Interpolating polynomials
    TYPE(LPoly_type), ALLOCATABLE :: wlp(:)  ! Angle        (I) – one per angle
    TYPE(LPoly_type)              :: xlp     ! Frequency    (L)
    TYPE(LPoly_type)              :: ylp     ! Grain_Size   (G)
    TYPE(LPoly_type)              :: tlp     ! Depth        (T)
    TYPE(LPoly_type)              :: zlp     ! Density      (J)
    ! LUT interpolation indices
    INTEGER, ALLOCATABLE :: i1(:), i2(:)     ! Angle
    INTEGER              :: j1,    j2        ! Frequency
    INTEGER              :: k1,    k2        ! Grain_Size
    INTEGER              :: l1,    l2        ! Depth
    INTEGER              :: m1,    m2        ! Density
    ! Out-of-bounds flags
    LOGICAL, ALLOCATABLE :: a_outbound(:)    ! Angle
    LOGICAL              :: f_outbound       ! Frequency
    LOGICAL              :: r_outbound       ! Grain_Size
    LOGICAL              :: d_outbound       ! Depth
    LOGICAL              :: rho_outbound     ! Density
    ! Scalar interpolation inputs (saved for TL/AD)
    REAL(fp), ALLOCATABLE :: a_int(:)        ! Angle
    REAL(fp)              :: f_int           ! Frequency
    REAL(fp)              :: r_int           ! Grain_Size
    REAL(fp)              :: d_int           ! Depth
    REAL(fp)              :: rho_int         ! Density
    ! LUT sub-array data (saved for TL/AD)
    REAL(fp), ALLOCATABLE :: a(:,:)          ! Angle         (NPTS, n_Angles)
    REAL(fp)              :: f(NPTS)         ! Frequency     (NPTS)
    REAL(fp)              :: r(NPTS)         ! Grain_Size    (NPTS)
    REAL(fp)              :: d(NPTS)         ! Depth         (NPTS)
    REAL(fp)              :: rho(NPTS)       ! Density       (NPTS)
  END TYPE Einterp_type

  ! Outer iVar_type wrapping Einterp (keeps TL/AD interface stable)
  TYPE :: iVar_type
    PRIVATE
    TYPE(Einterp_type) :: ei
  END TYPE iVar_type


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
!       CRTM_Compute_VISsnowRF
!
! PURPOSE:
!       Function to compute the CRTM visible/near-IR snow surface reflectance
!       for input grain size, depth, density, frequency, and angles.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_Compute_VISsnowRefl( VISsnowCoeff    , &
!                                                Snow_Grain_Size , &
!                                                Snow_Depth      , &
!                                                Snow_Density    , &
!                                                Frequency       , &
!                                                Angle           , &
!                                                iVar            , &
!                                                Reflectance       )
!
! INPUTS:
!       VISsnowCoeff:     Visible/near-IR snow reflectance coefficient object.
!                         UNITS:      N/A
!                         TYPE:       VISsnowCoeff_type
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN)
!
!     Snow_Grain_Size:    Snow grain effective radius.
!                         UNITS:      microns (um)
!                         TYPE:       REAL(fp)
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN)
!
!         Snow_Depth:     Snow layer depth.
!                         UNITS:      metres (m)
!                         TYPE:       REAL(fp)
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN)
!
!       Snow_Density:     Snow bulk density.
!                         UNITS:      kg m^-3
!                         TYPE:       REAL(fp)
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN)
!
!         Frequency:      Solar/visible channel frequency.
!                         UNITS:      inverse centimetres (cm^-1)
!                         TYPE:       REAL(fp)
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN)
!
!             Angle:      Surface zenith angles.
!                         UNITS:      Degrees
!                         TYPE:       REAL(fp)
!                         DIMENSION:  Rank-1 (n_Angles)
!                         ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!              iVar:      Structure containing internal variables required for
!                         subsequent tangent-linear or adjoint model calls.
!                         The contents of this structure are NOT accessible
!                         outside of this module.
!                         UNITS:      N/A
!                         TYPE:       iVar_type
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(OUT)
!
!       Reflectance:      Snow surface reflectances for the requested grain
!                         size, depth, density, frequency, and angles.
!                         UNITS:      N/A
!                         TYPE:       REAL(fp)
!                         DIMENSION:  Same as input Angle argument.
!                         ATTRIBUTES: INTENT(OUT)
!
! FUNCTION RESULT:
!       Error_Status:     The return value is an integer defining the error status.
!                         If == SUCCESS the computation was successful.
!                            == FAILURE an unrecoverable error occurred.
!                         UNITS:      N/A
!                         TYPE:       INTEGER
!                         DIMENSION:  Scalar
!
!:sdoc-:
!--------------------------------------------------------------------------------

  FUNCTION CRTM_Compute_VISsnowRF( &
    VISsnowCoeff    , &  ! Input
    Snow_Grain_Size , &  ! Input
    Snow_Depth      , &  ! Input
    Snow_Density    , &  ! Input
    Frequency       , &  ! Input
    Angle           , &  ! Input
    iVar            , &  ! Internal variable output
    Reflectance     ) &  ! Output
  RESULT( err_stat )
    ! Arguments
    TYPE(VISsnowCoeff_type), INTENT(IN)  :: VISsnowCoeff
    REAL(fp)               , INTENT(IN)  :: Snow_Grain_Size   ! r
    REAL(fp)               , INTENT(IN)  :: Snow_Depth        ! d
    REAL(fp)               , INTENT(IN)  :: Snow_Density      ! rho
    REAL(fp)               , INTENT(IN)  :: Frequency         ! f
    REAL(fp)               , INTENT(IN)  :: Angle(:)          ! a
    TYPE(iVar_type)        , INTENT(OUT) :: iVar
    REAL(fp)               , INTENT(OUT) :: Reflectance(:)
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Compute_VISsnowRF'
    ! Local variables
    CHARACTER(ML) :: msg
    INTEGER :: n_Angles, i

    ! Set up
    err_stat = SUCCESS
    ! ...Check dimensions
    n_Angles = SIZE(Angle)
    IF ( SIZE(Reflectance) /= n_Angles ) THEN
      err_stat = FAILURE
      msg = 'Input Angle and output Reflectance array dimensions inconsistent.'
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
      RETURN
    END IF
    ! ...Allocate interpolation variable structure
    CALL Einterp_Create( iVar%ei, NPTS, n_Angles )
    IF ( .NOT. Einterp_Associated( iVar%ei ) ) THEN
      err_stat = FAILURE
      msg = 'Error allocating interpolation variable structure.'
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
      RETURN
    END IF


    ! Compute the density interpolating polynomial (outermost / 5th dim)
    iVar%ei%rho_int = Snow_Density
    CALL Find_Index( VISsnowCoeff%Density, &
                     iVar%ei%rho_int, iVar%ei%m1, iVar%ei%m2, iVar%ei%rho_outbound )
    iVar%ei%rho = VISsnowCoeff%Density( iVar%ei%m1:iVar%ei%m2 )
    CALL LPoly( iVar%ei%rho    , &
                iVar%ei%rho_int, &
                iVar%ei%zlp      )


    ! Compute the depth interpolating polynomial (4th dim)
    iVar%ei%d_int = Snow_Depth
    CALL Find_Index( VISsnowCoeff%Depth, &
                     iVar%ei%d_int, iVar%ei%l1, iVar%ei%l2, iVar%ei%d_outbound )
    iVar%ei%d = VISsnowCoeff%Depth( iVar%ei%l1:iVar%ei%l2 )
    CALL LPoly( iVar%ei%d    , &
                iVar%ei%d_int, &
                iVar%ei%tlp    )


    ! Compute the grain size interpolating polynomial (3rd dim)
    iVar%ei%r_int = Snow_Grain_Size
    CALL Find_Index( VISsnowCoeff%Grain_Size, &
                     iVar%ei%r_int, iVar%ei%k1, iVar%ei%k2, iVar%ei%r_outbound )
    iVar%ei%r = VISsnowCoeff%Grain_Size( iVar%ei%k1:iVar%ei%k2 )
    CALL LPoly( iVar%ei%r    , &
                iVar%ei%r_int, &
                iVar%ei%ylp    )


    ! Compute the frequency interpolating polynomial (2nd dim)
    iVar%ei%f_int = Frequency
    CALL Find_Index( VISsnowCoeff%Frequency, &
                     iVar%ei%f_int, iVar%ei%j1, iVar%ei%j2, iVar%ei%f_outbound )
    iVar%ei%f = VISsnowCoeff%Frequency( iVar%ei%j1:iVar%ei%j2 )
    CALL LPoly( iVar%ei%f    , &
                iVar%ei%f_int, &
                iVar%ei%xlp    )


    ! Loop over angles (1st / innermost dim of LUT)
    DO i = 1, n_Angles

      ! Find index and compute angle polynomial
      iVar%ei%a_int(i) = ABS( Angle(i) )
      CALL Find_Index( VISsnowCoeff%Angle, &
                       iVar%ei%a_int(i), iVar%ei%i1(i), iVar%ei%i2(i), iVar%ei%a_outbound(i) )
      iVar%ei%a(:,i) = VISsnowCoeff%Angle( iVar%ei%i1(i):iVar%ei%i2(i) )
      CALL LPoly( iVar%ei%a(:,i)  , &
                  iVar%ei%a_int(i), &
                  iVar%ei%wlp(i)    )

      ! 5-D interpolation decomposed as 4-D over inner dims + 1-D over density
      CALL Interp_5D( VISsnowCoeff%Reflectance( iVar%ei%i1(i):iVar%ei%i2(i) , &
                                                iVar%ei%j1   :iVar%ei%j2   , &
                                                iVar%ei%k1   :iVar%ei%k2   , &
                                                iVar%ei%l1   :iVar%ei%l2   , &
                                                iVar%ei%m1   :iVar%ei%m2  ), &
                      iVar%ei%wlp(i), &   ! Angle polynomial
                      iVar%ei%xlp   , &   ! Frequency polynomial
                      iVar%ei%ylp   , &   ! Grain_Size polynomial
                      iVar%ei%tlp   , &   ! Depth polynomial
                      iVar%ei%zlp   , &   ! Density polynomial
                      Reflectance(i)  )

    END DO

  END FUNCTION CRTM_Compute_VISsnowRF


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Compute_VISsnowRefl_TL
!
! PURPOSE:
!       Function to compute the tangent-linear CRTM visible/near-IR snow
!       reflectance for input grain size, depth, density, frequency, and
!       angles.
!
!       This function must be called *after* the forward model function,
!       CRTM_Compute_VISsnowRefl, has been called.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_Compute_VISsnowRefl_TL( VISsnowCoeff        , &
!                                                   Snow_Grain_Size_TL  , &
!                                                   Snow_Depth_TL       , &
!                                                   Snow_Density_TL     , &
!                                                   iVar                , &
!                                                   Reflectance_TL        )
!
! INPUTS:
!       VISsnowCoeff:         Visible/near-IR snow reflectance coefficient object.
!                             UNITS:      N/A
!                             TYPE:       VISsnowCoeff_type
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(IN)
!
!   Snow_Grain_Size_TL:       Tangent-linear snow grain size.
!                             UNITS:      microns (um)
!                             TYPE:       REAL(fp)
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(IN)
!
!       Snow_Depth_TL:        Tangent-linear snow depth.
!                             UNITS:      metres (m)
!                             TYPE:       REAL(fp)
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(IN)
!
!     Snow_Density_TL:        Tangent-linear snow density.
!                             UNITS:      kg m^-3
!                             TYPE:       REAL(fp)
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(IN)
!
!                iVar:        Structure containing internal variables from
!                             the forward model call.
!                             UNITS:      N/A
!                             TYPE:       iVar_type
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       Reflectance_TL:       Tangent-linear snow surface reflectance.
!                             UNITS:      N/A
!                             TYPE:       REAL(fp)
!                             DIMENSION:  Rank-1 (n_Angles)
!                             ATTRIBUTES: INTENT(OUT)
!
! FUNCTION RESULT:
!       Error_Status:         If == SUCCESS the computation was successful.
!                             If == FAILURE an unrecoverable error occurred.
!
!:sdoc-:
!--------------------------------------------------------------------------------

  FUNCTION CRTM_Compute_VISsnowRF_TL( &
    VISsnowCoeff       , &  ! Input
    Snow_Grain_Size_TL , &  ! Input
    Snow_Depth_TL      , &  ! Input
    Snow_Density_TL    , &  ! Input
    iVar               , &  ! Internal variable input
    Reflectance_TL     ) &  ! Output
  RESULT( err_stat )
    ! Arguments
    TYPE(VISsnowCoeff_type), INTENT(IN)  :: VISsnowCoeff
    REAL(fp)               , INTENT(IN)  :: Snow_Grain_Size_TL
    REAL(fp)               , INTENT(IN)  :: Snow_Depth_TL
    REAL(fp)               , INTENT(IN)  :: Snow_Density_TL
    TYPE(iVar_type)        , INTENT(IN)  :: iVar
    REAL(fp)               , INTENT(OUT) :: Reflectance_TL(:)
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Compute_VISsnowRefl_TL'
    ! Local variables
    CHARACTER(ML) :: msg
    INTEGER  :: i
    REAL(fp) :: r_TL(NPTS), d_TL(NPTS), rho_TL(NPTS)
    REAL(fp) :: e_TL_5D(NPTS,NPTS,NPTS,NPTS,NPTS)
    TYPE(LPoly_type) :: wlp_TL, xlp_TL, ylp_TL, tlp_TL, zlp_TL

    ! Set up
    err_stat = SUCCESS
    ! ...Check internal variable allocation
    IF ( .NOT. Einterp_Associated( iVar%ei ) ) THEN
      err_stat = FAILURE
      msg = 'Internal structure ei is not allocated'
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
      RETURN
    END IF
    ! ...Check dimensions
    IF ( SIZE( Reflectance_TL ) /= iVar%ei%n_Angles ) THEN
      err_stat = FAILURE
      msg = 'Reflectance_TL array dimensions inconsistent with number of angles.'
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
      RETURN
    END IF
    ! ...No TL if any input is out of LUT bounds
    IF ( iVar%ei%r_outbound .OR. iVar%ei%d_outbound .OR. iVar%ei%rho_outbound ) THEN
      Reflectance_TL = ZERO
      RETURN
    END IF
    ! ...Initialise local TL variables
    r_TL   = ZERO
    d_TL   = ZERO
    rho_TL = ZERO
    e_TL_5D = ZERO
    CALL Clear_LPoly(wlp_TL)
    CALL Clear_LPoly(xlp_TL)
    CALL Clear_LPoly(ylp_TL)
    CALL Clear_LPoly(tlp_TL)
    CALL Clear_LPoly(zlp_TL)


    ! TL interpolating polynomial for density (5th dim)
    CALL LPoly_TL( iVar%ei%rho, iVar%ei%rho_int, &
                   iVar%ei%zlp,                  &
                   rho_TL, Snow_Density_TL,      &
                   zlp_TL                        )

    ! TL interpolating polynomial for depth (4th dim)
    CALL LPoly_TL( iVar%ei%d, iVar%ei%d_int, &
                   iVar%ei%tlp,              &
                   d_TL, Snow_Depth_TL,      &
                   tlp_TL                    )

    ! TL interpolating polynomial for grain size (3rd dim)
    CALL LPoly_TL( iVar%ei%r, iVar%ei%r_int, &
                   iVar%ei%ylp,              &
                   r_TL, Snow_Grain_Size_TL, &
                   ylp_TL                    )


    ! Loop over angles
    DO i = 1, iVar%ei%n_Angles

      CALL Interp_5D_TL( VISsnowCoeff%Reflectance( iVar%ei%i1(i):iVar%ei%i2(i), &
                                                   iVar%ei%j1   :iVar%ei%j2   , &
                                                   iVar%ei%k1   :iVar%ei%k2   , &
                                                   iVar%ei%l1   :iVar%ei%l2   , &
                                                   iVar%ei%m1   :iVar%ei%m2  ), &
                         iVar%ei%wlp(i), &   ! FWD polynomials
                         iVar%ei%xlp   , &
                         iVar%ei%ylp   , &
                         iVar%ei%tlp   , &
                         iVar%ei%zlp   , &
                         e_TL_5D, wlp_TL, xlp_TL, ylp_TL, tlp_TL, zlp_TL, & ! TL input
                         Reflectance_TL(i) )  ! TL output

    END DO

  END FUNCTION CRTM_Compute_VISsnowRF_TL


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Compute_VISsnowRF_AD
!
! PURPOSE:
!       Function to compute the adjoint of the CRTM visible/near-IR snow
!       reflectance for input grain size, depth, density, frequency, and
!       angles.
!
!       This function must be called *after* the forward model function,
!       CRTM_Compute_VISsnowRF, has been called.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_Compute_VISsnowRF_AD( VISsnowCoeff        , &
!                                                 Reflectance_AD      , &
!                                                 iVar                , &
!                                                 Snow_Grain_Size_AD  , &
!                                                 Snow_Depth_AD       , &
!                                                 Snow_Density_AD       )
!
! INPUTS:
!       VISsnowCoeff:         Visible/near-IR snow reflectance coefficient object.
!                             UNITS:      N/A
!                             TYPE:       VISsnowCoeff_type
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(IN)
!
!     Reflectance_AD:         Adjoint snow surface reflectance.
!                             *** SET TO ZERO ON EXIT ***
!                             UNITS:      N/A
!                             TYPE:       REAL(fp)
!                             DIMENSION:  Rank-1 (n_Angles)
!                             ATTRIBUTES: INTENT(IN OUT)
!
!               iVar:         Structure containing internal variables from
!                             the forward model call.
!                             UNITS:      N/A
!                             TYPE:       iVar_type
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!   Snow_Grain_Size_AD:       Adjoint snow grain size.
!                             *** MUST HAVE VALUE ON ENTRY ***
!                             UNITS:      per micron (um^-1)
!                             TYPE:       REAL(fp)
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(IN OUT)
!
!       Snow_Depth_AD:        Adjoint snow depth.
!                             *** MUST HAVE VALUE ON ENTRY ***
!                             UNITS:      per metre (m^-1)
!                             TYPE:       REAL(fp)
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(IN OUT)
!
!     Snow_Density_AD:        Adjoint snow density.
!                             *** MUST HAVE VALUE ON ENTRY ***
!                             UNITS:      per (kg m^-3)^-1
!                             TYPE:       REAL(fp)
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(IN OUT)
!
! FUNCTION RESULT:
!       Error_Status:         If == SUCCESS the computation was successful.
!                             If == FAILURE an unrecoverable error occurred.
!
!:sdoc-:
!--------------------------------------------------------------------------------

  FUNCTION CRTM_Compute_VISsnowRF_AD( &
    VISsnowCoeff       , &  ! Input
    Reflectance_AD     , &  ! Input  (zeroed on exit)
    iVar               , &  ! Internal variable input
    Snow_Grain_Size_AD , &  ! Output (accumulates)
    Snow_Depth_AD      , &  ! Output (accumulates)
    Snow_Density_AD    ) &  ! Output (accumulates)
  RESULT( err_stat )
    ! Arguments
    TYPE(VISsnowCoeff_type), INTENT(IN)     :: VISsnowCoeff
    REAL(fp)               , INTENT(IN OUT) :: Reflectance_AD(:)
    TYPE(iVar_type)        , INTENT(IN)     :: iVar
    REAL(fp)               , INTENT(IN OUT) :: Snow_Grain_Size_AD
    REAL(fp)               , INTENT(IN OUT) :: Snow_Depth_AD
    REAL(fp)               , INTENT(IN OUT) :: Snow_Density_AD
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Compute_VISsnowRF_AD'
    ! Local variables
    CHARACTER(ML) :: msg
    INTEGER  :: i
    REAL(fp) :: e_AD_5D(NPTS,NPTS,NPTS,NPTS,NPTS)
    REAL(fp) :: r_AD(NPTS), d_AD(NPTS), rho_AD(NPTS)
    TYPE(LPoly_type) :: wlp_AD, xlp_AD, ylp_AD, tlp_AD, zlp_AD

    ! Set up
    err_stat = SUCCESS
    e_AD_5D = ZERO
    r_AD    = ZERO
    d_AD    = ZERO
    rho_AD  = ZERO
    ! ...Check internal variable allocation
    IF ( .NOT. Einterp_Associated( iVar%ei ) ) THEN
      err_stat = FAILURE
      msg = 'Internal structure ei is not allocated'
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
      RETURN
    END IF
    ! ...Check dimensions
    IF ( SIZE(Reflectance_AD) /= iVar%ei%n_Angles ) THEN
      err_stat = FAILURE
      msg = 'Reflectance_AD array dimensions inconsistent with number of angles.'
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
      RETURN
    END IF
    ! ...No AD if any input was out of LUT bounds during forward call
    IF ( iVar%ei%r_outbound .OR. iVar%ei%d_outbound .OR. iVar%ei%rho_outbound ) RETURN
    ! ...Initialise local AD polynomial structures
    CALL Clear_LPoly(wlp_AD)
    CALL Clear_LPoly(xlp_AD)
    CALL Clear_LPoly(ylp_AD)
    CALL Clear_LPoly(tlp_AD)
    CALL Clear_LPoly(zlp_AD)


    ! Loop over angles – accumulate adjoint polynomials
    DO i = 1, iVar%ei%n_Angles

      CALL Interp_5D_AD( VISsnowCoeff%Reflectance( iVar%ei%i1(i):iVar%ei%i2(i), &
                                                   iVar%ei%j1   :iVar%ei%j2   , &
                                                   iVar%ei%k1   :iVar%ei%k2   , &
                                                   iVar%ei%l1   :iVar%ei%l2   , &
                                                   iVar%ei%m1   :iVar%ei%m2  ), &
                         iVar%ei%wlp(i), &   ! FWD polynomials
                         iVar%ei%xlp   , &
                         iVar%ei%ylp   , &
                         iVar%ei%tlp   , &
                         iVar%ei%zlp   , &
                         Reflectance_AD(i),        &   ! AD input
                         e_AD_5D, wlp_AD, xlp_AD, ylp_AD, tlp_AD, zlp_AD )  ! AD output

      ! Zero the adjoint reflectance for this angle
      Reflectance_AD(i) = ZERO

    END DO


    ! AD of grain size polynomial
    CALL LPoly_AD( iVar%ei%r    , &
                   iVar%ei%r_int, &
                   iVar%ei%ylp  , &
                   ylp_AD       , &
                   r_AD         , &
                   Snow_Grain_Size_AD )

    ! AD of depth polynomial
    CALL LPoly_AD( iVar%ei%d    , &
                   iVar%ei%d_int, &
                   iVar%ei%tlp  , &
                   tlp_AD       , &
                   d_AD         , &
                   Snow_Depth_AD  )

    ! AD of density polynomial
    CALL LPoly_AD( iVar%ei%rho    , &
                   iVar%ei%rho_int, &
                   iVar%ei%zlp    , &
                   zlp_AD         , &
                   rho_AD         , &
                   Snow_Density_AD  )

  END FUNCTION CRTM_Compute_VISsnowRF_AD


!################################################################################
!################################################################################
!##                                                                            ##
!##                        ## PRIVATE MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

  ! --------------------------------------------------
  ! Einterp allocation helpers
  ! --------------------------------------------------
  ELEMENTAL FUNCTION Einterp_Associated( ei ) RESULT( Status )
    TYPE(Einterp_type), INTENT(IN) :: ei
    LOGICAL :: Status
    Status = ei%Is_Allocated
  END FUNCTION Einterp_Associated

  ELEMENTAL SUBROUTINE Einterp_Create( ei, n_Pts, n_Angles )
    TYPE(Einterp_type), INTENT(OUT) :: ei
    INTEGER,            INTENT(IN)  :: n_Pts
    INTEGER,            INTENT(IN)  :: n_Angles
    INTEGER :: alloc_stat
    IF ( n_Pts < 1 .OR. n_Angles < 1 ) RETURN
    ALLOCATE( ei%wlp(n_Angles)        , &
              ei%i1(n_Angles)         , &
              ei%i2(n_Angles)         , &
              ei%a_outbound(n_Angles) , &
              ei%a_int(n_Angles)      , &
              ei%a(n_Pts, n_Angles)   , &
              STAT = alloc_stat )
    IF ( alloc_stat /= 0 ) RETURN
    ei%n_Angles     = n_Angles
    ei%n_Pts        = n_Pts
    ei%Is_Allocated = .TRUE.
  END SUBROUTINE Einterp_Create


  ! --------------------------------------------------
  ! 5-D forward interpolation
  !
  !   z(u,v,w,x,y)  with polynomials ulp,vlp,wlp,xlp,ylp
  !
  ! Implemented by looping over the 5th (y) dimension and
  ! calling Interp_4D, then collapsing with Interp_1D.
  ! --------------------------------------------------
  SUBROUTINE Interp_5D( z, ulp, vlp, wlp, xlp, ylp, &
                        z_int )
    REAL(fp)        ,  INTENT(IN)     :: z(:,:,:,:,:)
    TYPE(LPoly_type),  INTENT(IN)     :: ulp, vlp, wlp, xlp, ylp
    REAL(fp)        ,  INTENT(IN OUT) :: z_int  ! INTENT(IN OUT) to preclude reinitialisation
    ! Local variables
    INTEGER  :: i
    REAL(fp) :: a(NPTS)
    ! Interpolate in u,v,w,x for each y slice
    DO i = 1, NPTS
      CALL Interp_4D( z(:,:,:,:,i), ulp, vlp, wlp, xlp, a(i) )
    END DO
    ! Interpolate the resulting vector in y
    CALL Interp_1D( a, ylp, z_int )
  END SUBROUTINE Interp_5D


  ! --------------------------------------------------
  ! 5-D tangent-linear interpolation
  ! --------------------------------------------------
  SUBROUTINE Interp_5D_TL( z   , ulp   , vlp   , wlp   , xlp   , ylp   , &
                           z_TL, ulp_TL, vlp_TL, wlp_TL, xlp_TL, ylp_TL, &
                           z_int_TL )
    REAL(fp)        ,  INTENT(IN)     :: z(:,:,:,:,:)
    TYPE(LPoly_type),  INTENT(IN)     :: ulp, vlp, wlp, xlp, ylp
    REAL(fp)        ,  INTENT(IN)     :: z_TL(:,:,:,:,:)
    TYPE(LPoly_type),  INTENT(IN)     :: ulp_TL, vlp_TL, wlp_TL, xlp_TL, ylp_TL
    REAL(fp)        ,  INTENT(IN OUT) :: z_int_TL  ! INTENT(IN OUT) to preclude reinitialisation
    ! Local variables
    INTEGER  :: i
    REAL(fp) :: a(NPTS), a_TL(NPTS)
    ! Forward and TL in u,v,w,x for each y slice
    DO i = 1, NPTS
      CALL Interp_4D   ( z(:,:,:,:,i), ulp, vlp, wlp, xlp, a(i) )
      CALL Interp_4D_TL( z(:,:,:,:,i), ulp, vlp, wlp, xlp, &
                         z_TL(:,:,:,:,i), ulp_TL, vlp_TL, wlp_TL, xlp_TL, a_TL(i) )
    END DO
    ! TL collapse in y
    CALL Interp_1D_TL( a, ylp, a_TL, ylp_TL, z_int_TL )
  END SUBROUTINE Interp_5D_TL


  ! --------------------------------------------------
  ! 5-D adjoint interpolation
  ! --------------------------------------------------
  SUBROUTINE Interp_5D_AD( z      , ulp   , vlp   , wlp   , xlp   , ylp   , &
                           z_int_AD,                                          &
                           z_AD   , ulp_AD, vlp_AD, wlp_AD, xlp_AD, ylp_AD  )
    REAL(fp)        ,  INTENT(IN)     :: z(:,:,:,:,:)
    TYPE(LPoly_type),  INTENT(IN)     :: ulp, vlp, wlp, xlp, ylp
    REAL(fp)        ,  INTENT(IN OUT) :: z_int_AD
    REAL(fp)        ,  INTENT(IN OUT) :: z_AD(:,:,:,:,:)
    TYPE(LPoly_type),  INTENT(IN OUT) :: ulp_AD, vlp_AD, wlp_AD, xlp_AD, ylp_AD
    ! Local variables
    INTEGER  :: i
    REAL(fp) :: a(NPTS), a_AD(NPTS)

    ! Forward pass: build a(i) = Interp_4D over y slices
    DO i = 1, NPTS
      CALL Interp_4D( z(:,:,:,:,i), ulp, vlp, wlp, xlp, a(i) )
    END DO

    ! Adjoint initialisation
    a_AD = ZERO

    ! Adjoint collapse in y
    CALL Interp_1D_AD( a, ylp, z_int_AD, a_AD, ylp_AD )

    ! Adjoint of inner 4-D interpolation for each y slice
    DO i = 1, NPTS
      CALL Interp_4D_AD( z(:,:,:,:,i), ulp, vlp, wlp, xlp, &
                         a_AD(i), &
                         z_AD(:,:,:,:,i), ulp_AD, vlp_AD, wlp_AD, xlp_AD )
    END DO

  END SUBROUTINE Interp_5D_AD

END MODULE CRTM_VISsnowRF
