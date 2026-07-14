!
! CRTM_Surface_Define
!
! Module defining the CRTM Surface structure and containing routines
! to manipulate it.
!
!
! CREATION HISTORY:
!       Written by:  Yong Han,       yong.han@noaa.gov
!                    Quanhua Liu,    quanhua.liu@noaa.gov
!                    Paul van Delst, paul.vandelst@noaa.gov
!                    07-May-2004
!

MODULE CRTM_Surface_Define

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Intrinsic modules
  USE ISO_Fortran_Env       , ONLY: OUTPUT_UNIT
  ! Module use
  USE netcdf
  USE Type_Kinds            , ONLY: fp
  USE Message_Handler       , ONLY: SUCCESS, FAILURE, WARNING, INFORMATION, Display_Message
  USE Compare_Float_Numbers , ONLY: DEFAULT_N_SIGFIG, &
                                    OPERATOR(.EqualTo.), &
                                    Compares_Within_Tolerance
  USE File_Utility          , ONLY: File_Open, File_Exists
  USE Binary_File_Utility   , ONLY: Open_Binary_File      , &
                                    WriteGAtts_Binary_File, &
                                    ReadGAtts_Binary_File
  USE CRTM_SensorData_Define, ONLY: CRTM_SensorData_type, &
                                    OPERATOR(==), &
                                    OPERATOR(+), &
                                    OPERATOR(-), &
                                    CRTM_SensorData_Associated, &
                                    CRTM_SensorData_Destroy, &
                                    CRTM_SensorData_Create, &
                                    CRTM_SensorData_Zero, &
                                    CRTM_SensorData_IsValid, &
                                    CRTM_SensorData_Inspect, &
                                    CRTM_SensorData_Compare, &
                                    CRTM_SensorData_ReadFile, &
                                    CRTM_SensorData_WriteFile
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Operators
  PUBLIC :: OPERATOR(==)
  PUBLIC :: OPERATOR(+)
  PUBLIC :: OPERATOR(-)
  ! SensorData enitities
  ! ...Structures
  PUBLIC :: CRTM_SensorData_type
  ! ...Procedures
  PUBLIC :: CRTM_SensorData_Associated
  PUBLIC :: CRTM_SensorData_Destroy
  PUBLIC :: CRTM_SensorData_Create
  PUBLIC :: CRTM_SensorData_Zero
  PUBLIC :: CRTM_SensorData_IsValid
  PUBLIC :: CRTM_SensorData_Inspect
  PUBLIC :: CRTM_SensorData_Compare
  ! Surface entities
  ! ...Gross surface parameters
  PUBLIC :: INVALID_SURFACE
  PUBLIC :: LAND_SURFACE
  PUBLIC :: WATER_SURFACE
  PUBLIC :: SNOW_SURFACE
  PUBLIC :: ICE_SURFACE
  PUBLIC :: N_VALID_SURFACE_TYPES
  PUBLIC :: SURFACE_TYPE_NAME
  ! ...Structures
  PUBLIC :: CRTM_Surface_type
  ! ...Procedures
  PUBLIC :: CRTM_Surface_Associated
  PUBLIC :: CRTM_Surface_Destroy
  PUBLIC :: CRTM_Surface_Create
  PUBLIC :: CRTM_Surface_NonVariableCopy
  PUBLIC :: CRTM_Surface_Zero
  PUBLIC :: CRTM_Surface_IsValid
  PUBLIC :: CRTM_Surface_Inspect
  PUBLIC :: CRTM_Surface_IsCoverageValid
  PUBLIC :: CRTM_Surface_CoverageType
  PUBLIC :: CRTM_Surface_Compare
  PUBLIC :: CRTM_Surface_InquireFile
  PUBLIC :: CRTM_Surface_ReadFile
  PUBLIC :: CRTM_Surface_WriteFile


  ! ---------------------
  ! Procedure overloading
  ! ---------------------
  INTERFACE OPERATOR(==)
    MODULE PROCEDURE CRTM_Surface_Equal
  END INTERFACE OPERATOR(==)

  INTERFACE OPERATOR(+)
    MODULE PROCEDURE CRTM_Surface_Add
  END INTERFACE OPERATOR(+)

  INTERFACE OPERATOR(-)
    MODULE PROCEDURE CRTM_Surface_Subtract
  END INTERFACE OPERATOR(-)

  INTERFACE CRTM_Surface_Inspect
    MODULE PROCEDURE Scalar_Inspect
    MODULE PROCEDURE Rank1_Inspect
    MODULE PROCEDURE Rank2_Inspect
  END INTERFACE CRTM_Surface_Inspect

  INTERFACE CRTM_Surface_ReadFile
    MODULE PROCEDURE Read_Surface_Rank1
    MODULE PROCEDURE Read_Surface_Rank2
  END INTERFACE CRTM_Surface_ReadFile

  INTERFACE CRTM_Surface_WriteFile
    MODULE PROCEDURE Write_Surface_Rank1
    MODULE PROCEDURE Write_Surface_Rank2
  END INTERFACE CRTM_Surface_WriteFile


  ! -----------------
  ! Module parameters
  ! -----------------
  ! Literal constants
  REAL(fp), PARAMETER :: ZERO = 0.0_fp
  REAL(fp), PARAMETER :: ONE  = 1.0_fp
  ! Message string length
  INTEGER, PARAMETER :: ML = 256
  ! File status on close after write error
  CHARACTER(*), PARAMETER :: WRITE_ERROR_STATUS = 'DELETE'

  ! ---------------------------------------------------------------------------
  ! netCDF I/O schema (used by the *_NetCDF file workers)
  !
  ! The Surface object is all-scalar-per-element, so every field is packed into
  ! a single rank-3 REAL variable Surface_Data(n_Channels, n_Profiles, n_Fields)
  ! with a fixed field ordering given by the IDX_* parameters below. INTEGER
  ! surface fields are stored as REAL(fp) and recovered with NINT on read (the
  ! values are small type codes, so the round-trip is exact). Because these
  ! files are written and read by the same build (ctest baselines), the only
  ! requirement is an exact round-trip; the schema is otherwise free.
  ! ---------------------------------------------------------------------------
  ! ...Dimension names
  CHARACTER(*), PARAMETER :: SFC_CHANNEL_DIMNAME = 'n_Channels'
  CHARACTER(*), PARAMETER :: SFC_PROFILE_DIMNAME = 'n_Profiles'
  CHARACTER(*), PARAMETER :: SFC_FIELD_DIMNAME   = 'n_Surface_Fields'
  ! ...Global attribute holding the true n_Channels. The channel dimension is
  !    MAX(n_Channels,1) so that a profile-only (rank-1) Surface, which has
  !    n_Channels==0, is representable (NF90_DEF_DIM treats 0 as UNLIMITED).
  CHARACTER(*), PARAMETER :: SFC_NCHANNELS_GATTNAME = 'n_Channels'
  ! ...Variable names
  CHARACTER(*), PARAMETER :: SFC_DATA_VARNAME    = 'Surface_Data'
  ! ...netCDF storage type for REAL(fp) data (fp is double; see Type_Kinds)
  INTEGER, PARAMETER :: SFC_FLOAT_TYPE = NF90_DOUBLE
  ! ...Packed field ordering
  INTEGER, PARAMETER :: IDX_LAND_COVERAGE          =  1
  INTEGER, PARAMETER :: IDX_WATER_COVERAGE         =  2
  INTEGER, PARAMETER :: IDX_SNOW_COVERAGE          =  3
  INTEGER, PARAMETER :: IDX_ICE_COVERAGE           =  4
  INTEGER, PARAMETER :: IDX_WIND_SPEED             =  5
  INTEGER, PARAMETER :: IDX_LAND_TEMPERATURE       =  6
  INTEGER, PARAMETER :: IDX_SOIL_MOISTURE_CONTENT  =  7
  INTEGER, PARAMETER :: IDX_CANOPY_WATER_CONTENT   =  8
  INTEGER, PARAMETER :: IDX_VEGETATION_FRACTION    =  9
  INTEGER, PARAMETER :: IDX_SOIL_TEMPERATURE       = 10
  INTEGER, PARAMETER :: IDX_LAI                    = 11
  INTEGER, PARAMETER :: IDX_WATER_TEMPERATURE      = 12
  INTEGER, PARAMETER :: IDX_WIND_DIRECTION         = 13
  INTEGER, PARAMETER :: IDX_SALINITY               = 14
  INTEGER, PARAMETER :: IDX_SNOW_TEMPERATURE       = 15
  INTEGER, PARAMETER :: IDX_SNOW_DEPTH             = 16
  INTEGER, PARAMETER :: IDX_SNOW_DENSITY           = 17
  INTEGER, PARAMETER :: IDX_SNOW_GRAIN_SIZE        = 18
  INTEGER, PARAMETER :: IDX_ICE_TEMPERATURE        = 19
  INTEGER, PARAMETER :: IDX_ICE_THICKNESS          = 20
  INTEGER, PARAMETER :: IDX_ICE_DENSITY            = 21
  INTEGER, PARAMETER :: IDX_ICE_ROUGHNESS          = 22
  INTEGER, PARAMETER :: IDX_LAND_TYPE              = 23
  INTEGER, PARAMETER :: IDX_SOIL_TYPE              = 24
  INTEGER, PARAMETER :: IDX_VEGETATION_TYPE        = 25
  INTEGER, PARAMETER :: IDX_WATER_TYPE             = 26
  INTEGER, PARAMETER :: IDX_SNOW_TYPE              = 27
  INTEGER, PARAMETER :: IDX_ICE_TYPE               = 28
  INTEGER, PARAMETER :: IDX_SENSORDATA_N_CHANNELS  = 29
  INTEGER, PARAMETER :: N_SURFACE_FIELDS           = 29

  ! The gross surface types. These are used for
  ! cross-checking with the coverage fractions
  ! of each gross surface types.
  INTEGER, PARAMETER :: INVALID_SURFACE = 0
  INTEGER, PARAMETER :: LAND_SURFACE    = 1
  INTEGER, PARAMETER :: WATER_SURFACE   = 2
  INTEGER, PARAMETER :: SNOW_SURFACE    = 4
  INTEGER, PARAMETER :: ICE_SURFACE     = 8
  INTEGER, PARAMETER :: N_VALID_SURFACE_TYPES = LAND_SURFACE  + &
                                                WATER_SURFACE + &
                                                SNOW_SURFACE  + &
                                                ICE_SURFACE
  CHARACTER(*), PARAMETER, DIMENSION( 0:N_VALID_SURFACE_TYPES ) :: &
    SURFACE_TYPE_NAME = (/ 'Invalid surface type     ', &
                           'Land                     ', &
                           'Water                    ', &
                           'Land + water             ', &
                           'Snow                     ', &
                           'Land + snow              ', &
                           'Water + snow             ', &
                           'Land + water + snow      ', &
                           'Ice                      ', &
                           'Land + ice               ', &
                           'Water + ice              ', &
                           'Land + water + ice       ', &
                           'Snow + ice               ', &
                           'Land + snow + ice        ', &
                           'Water + snow + ice       ', &
                           'Land + water + snow + ice' /)

  ! Default value parameters
  ! ...Land surface type data
  INTEGER,  PARAMETER :: DEFAULT_LAND_TYPE             = 1         ! First item in list
  REAL(fp), PARAMETER :: DEFAULT_LAND_TEMPERATURE      = 283.0_fp  ! K
  REAL(fp), PARAMETER :: DEFAULT_SOIL_MOISTURE_CONTENT = 0.05_fp   ! g/cm^3
  REAL(fp), PARAMETER :: DEFAULT_CANOPY_WATER_CONTENT  = 0.05_fp   ! g/cm^3
  REAL(fp), PARAMETER :: DEFAULT_VEGETATION_FRACTION   = 0.3_fp    ! 30%
  REAL(fp), PARAMETER :: DEFAULT_SOIL_TEMPERATURE      = 283.0_fp  ! K
  REAL(fp), PARAMETER :: DEFAULT_LAI                   = 3.5
  INTEGER,  PARAMETER :: DEFAULT_SOIL_TYPE             = 1         ! First item in list
  INTEGER,  PARAMETER :: DEFAULT_VEGETATION_TYPE       = 1         ! First item in list
  ! ...Water type data
  INTEGER,  PARAMETER :: DEFAULT_WATER_TYPE        = 1          ! First item in list
  REAL(fp), PARAMETER :: DEFAULT_WATER_TEMPERATURE = 283.0_fp   ! K
  REAL(fp), PARAMETER :: DEFAULT_WIND_SPEED        = 5.0_fp     ! m/s
  REAL(fp), PARAMETER :: DEFAULT_WIND_DIRECTION    = 0.0_fp     ! Southerly wind, i.e. FROM the south. Opposite from met. defn.
  REAL(fp), PARAMETER :: DEFAULT_SALINITY          = 33.0_fp    ! ppmv
  ! ...Snow surface type data
  INTEGER,  PARAMETER :: DEFAULT_SNOW_TYPE        = 1          ! First item in list
  REAL(fp), PARAMETER :: DEFAULT_SNOW_TEMPERATURE = 263.0_fp   ! K
  REAL(fp), PARAMETER :: DEFAULT_SNOW_DEPTH       = 50.0_fp    ! mm
  REAL(fp), PARAMETER :: DEFAULT_SNOW_DENSITY     = 0.2_fp     ! g/cm^3
  REAL(fp), PARAMETER :: DEFAULT_SNOW_GRAIN_SIZE  = 2.0e3_fp   ! um (microns)
  ! ...Ice surface type data
  INTEGER,  PARAMETER :: DEFAULT_ICE_TYPE        = 1         ! First item in list
  REAL(fp), PARAMETER :: DEFAULT_ICE_TEMPERATURE = 263.0_fp  ! K
  REAL(fp), PARAMETER :: DEFAULT_ICE_THICKNESS   = 10.0_fp   ! mm
  REAL(fp), PARAMETER :: DEFAULT_ICE_DENSITY     = 0.9_fp    ! g/cm^3
  REAL(fp), PARAMETER :: DEFAULT_ICE_ROUGHNESS   = ZERO


  ! ----------------------------
  ! Surface structure definition
  ! ----------------------------
  !:tdoc+:
  TYPE :: CRTM_Surface_type
    ! Allocation indicator
    LOGICAL :: Is_Allocated = .TRUE.  ! Placeholder for future expansion
    ! Dimension values
    ! ...None yet
    ! Gross type of surface determined by coverage
    REAL(fp) :: Land_Coverage  = ZERO
    REAL(fp) :: Water_Coverage = ZERO
    REAL(fp) :: Snow_Coverage  = ZERO
    REAL(fp) :: Ice_Coverage   = ZERO
    ! Land surface type data
    INTEGER  :: Land_Type             = DEFAULT_LAND_TYPE
    REAL(fp) :: Land_Temperature      = DEFAULT_LAND_TEMPERATURE
    REAL(fp) :: Soil_Moisture_Content = DEFAULT_SOIL_MOISTURE_CONTENT
    REAL(fp) :: Canopy_Water_Content  = DEFAULT_CANOPY_WATER_CONTENT
    REAL(fp) :: Vegetation_Fraction   = DEFAULT_VEGETATION_FRACTION
    REAL(fp) :: Soil_Temperature      = DEFAULT_SOIL_TEMPERATURE
    REAL(fp) :: LAI                   = DEFAULT_LAI
    INTEGER  :: Soil_Type             = DEFAULT_SOIL_TYPE
    INTEGER  :: Vegetation_Type       = DEFAULT_VEGETATION_TYPE
    ! Water type data
    INTEGER  :: Water_Type        = DEFAULT_WATER_TYPE
    REAL(fp) :: Water_Temperature = DEFAULT_WATER_TEMPERATURE
    REAL(fp) :: Wind_Speed        = DEFAULT_WIND_SPEED
    REAL(fp) :: Wind_Direction    = DEFAULT_WIND_DIRECTION
    REAL(fp) :: Salinity          = DEFAULT_SALINITY
    ! Snow surface type data
    INTEGER  :: Snow_Type        = DEFAULT_SNOW_TYPE
    REAL(fp) :: Snow_Temperature = DEFAULT_SNOW_TEMPERATURE
    REAL(fp) :: Snow_Depth       = DEFAULT_SNOW_DEPTH
    REAL(fp) :: Snow_Density     = DEFAULT_SNOW_DENSITY
    REAL(fp) :: Snow_Grain_Size  = DEFAULT_SNOW_GRAIN_SIZE
    ! Ice surface type data
    INTEGER  :: Ice_Type        = DEFAULT_ICE_TYPE
    REAL(fp) :: Ice_Temperature = DEFAULT_ICE_TEMPERATURE
    REAL(fp) :: Ice_Thickness   = DEFAULT_ICE_THICKNESS
    REAL(fp) :: Ice_Density     = DEFAULT_ICE_DENSITY
    REAL(fp) :: Ice_Roughness   = DEFAULT_ICE_ROUGHNESS
    ! SensorData containing channel brightness temperatures
    TYPE(CRTM_SensorData_type) :: SensorData
  END TYPE CRTM_Surface_type
  !:tdoc-:


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
!       CRTM_Surface_Associated
!
! PURPOSE:
!       Elemental function to test the status of the allocatable components
!       of a CRTM Surface object.
!
! CALLING SEQUENCE:
!       Status = CRTM_Surface_Associated( Sfc )
!
! OBJECTS:
!       Sfc:       Surface structure which is to have its member's
!                  status tested.
!                  UNITS:      N/A
!                  TYPE:       CRTM_Surface_type
!                  DIMENSION:  Scalar or any rank
!                  ATTRIBUTES: INTENT(IN)
!
! FUNCTION RESULT:
!       Status:    The return value is a logical value indicating the
!                  status of the Surface members.
!                    .TRUE.  - if the array components are allocated.
!                    .FALSE. - if the array components are not allocated.
!                  UNITS:      N/A
!                  TYPE:       LOGICAL
!                  DIMENSION:  Same as input
!
!:sdoc-:
!--------------------------------------------------------------------------------

  ELEMENTAL FUNCTION CRTM_Surface_Associated( Sfc ) RESULT( Status )
    TYPE(CRTM_Surface_type), INTENT(IN) :: Sfc
    LOGICAL :: Status

    Status = Sfc%Is_Allocated
    ! ...SensorData
    Status = Status .AND. CRTM_SensorData_Associated(Sfc%SensorData)

  END FUNCTION CRTM_Surface_Associated


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Surface_Destroy
!
! PURPOSE:
!       Elemental subroutine to re-initialize CRTM Surface objects.
!
! CALLING SEQUENCE:
!       CALL CRTM_Surface_Destroy( Sfc )
!
! OBJECTS:
!       Sfc:          Re-initialized Surface structure.
!                     UNITS:      N/A
!                     TYPE:       CRTM_Surface_type
!                     DIMENSION:  Scalar or any rank
!                     ATTRIBUTES: INTENT(OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE CRTM_Surface_Destroy( Sfc )
    TYPE(CRTM_Surface_type), INTENT(OUT) :: Sfc
    Sfc%Is_Allocated = .TRUE.  ! Placeholder for future expansion
  END SUBROUTINE CRTM_Surface_Destroy


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Surface_Create
!
! PURPOSE:
!       Elemental subroutine to create an instance of the CRTM Surface object.
!
! CALLING SEQUENCE:
!       CALL CRTM_Surface_Create( Sfc       , &
!                                 n_Channels  )
!
! OBJECTS:
!       Sfc:          Surface structure.
!                     UNITS:      N/A
!                     TYPE:       CRTM_Surface_type
!                     DIMENSION:  Scalar or any rank
!                     ATTRIBUTES: INTENT(OUT)
!
! INPUT ARGUMENTS:
!       n_Channels:   Number of channels dimension of SensorData
!                     substructure
!                     ** Note: Can be = 0 (i.e. no sensor data). **
!                     UNITS:      N/A
!                     TYPE:       INTEGER
!                     DIMENSION:  Same as Surface object
!                     ATTRIBUTES: INTENT(IN)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE CRTM_Surface_Create( &
    Sfc       , &  ! Output
    n_Channels  )  ! Input
    ! Arguments
    TYPE(CRTM_Surface_type), INTENT(OUT) :: Sfc
    INTEGER                , INTENT(IN)  :: n_Channels

    ! Check input
    IF ( n_Channels < 0 ) RETURN

    ! Perform the substructure allocation
    ! ...SensorData
    IF ( n_Channels > 0 ) CALL CRTM_SensorData_Create( Sfc%SensorData, n_Channels )

    ! Set allocation indicator
    Sfc%Is_Allocated = .TRUE.

  END SUBROUTINE CRTM_Surface_Create


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Surface_NonVariableCopy
!
! PURPOSE:
!       Elemental utility subroutine to copy the "non-variable" data (coverages
!       and surface types) from one instance of a CRTM Surface object to another
!       (usually a TL or AD one).
!
!       NOTE: No error checking is performed in this procedure.
!
! CALLING SEQUENCE:
!       CALL CRTM_Surface_NonVariableCopy( sfc, modified_sfc )
!
! OBJECTS:
!       sfc:             Surface object from which to copy.
!                        UNITS:      N/A
!                        TYPE:       CRTM_Surface_type
!                        DIMENSION:  Scalar or any rank
!                        ATTRIBUTES: INTENT(IN)
!
! IN/OUTPUTS:
!       modified_sfc:    Existing Surface object to be modified.
!                        UNITS:      N/A
!                        TYPE:       CRTM_Surface_type
!                        DIMENSION:  Conformable with sfc input
!                        ATTRIBUTES: INTENT(IN OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE CRTM_Surface_NonVariableCopy( sfc, modified_sfc )
    TYPE(CRTM_Surface_type), INTENT(IN)     :: sfc
    TYPE(CRTM_Surface_type), INTENT(IN OUT) :: modified_sfc

    modified_sfc%Land_Coverage  = sfc%Land_Coverage
    modified_sfc%Water_Coverage = sfc%Water_Coverage
    modified_sfc%Snow_Coverage  = sfc%Snow_Coverage
    modified_sfc%Ice_Coverage   = sfc%Ice_Coverage

    modified_sfc%Land_Type  = sfc%Land_Type
    modified_sfc%Water_Type = sfc%Water_Type
    modified_sfc%Snow_Type  = sfc%Snow_Type
    modified_sfc%Ice_Type   = sfc%Ice_Type

  END SUBROUTINE CRTM_Surface_NonVariableCopy


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Surface_Zero
!
! PURPOSE:
!       Elemental subroutine to zero out the data arrays
!       in a CRTM Surface object.
!
! CALLING SEQUENCE:
!       CALL CRTM_Surface_Zero( Sfc )
!
! OUTPUT ARGUMENTS:
!       Sfc:          CRTM Surface structure in which the data arrays
!                     are to be zeroed out.
!                     UNITS:      N/A
!                     TYPE:       CRTM_Surface_type
!                     DIMENSION:  Scalar or any rank
!                     ATTRIBUTES: INTENT(IN OUT)
!
! COMMENTS:
!       - The various surface type indicator flags are
!         *NOT* reset in this routine.
!
!:sdoc-:
!--------------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE CRTM_Surface_Zero( Sfc )
    TYPE(CRTM_Surface_type), INTENT(IN OUT) :: Sfc

    ! Zero the components
    ! ...Coverage fractions
    Sfc%Land_Coverage  = ZERO
    Sfc%Water_Coverage = ZERO
    Sfc%Snow_Coverage  = ZERO
    Sfc%Ice_Coverage   = ZERO
    ! ...The various surface types
    CALL CRTM_LandSurface_Zero(sfc)
    CALL CRTM_WaterSurface_Zero(sfc)
    CALL CRTM_SnowSurface_Zero(sfc)
    CALL CRTM_IceSurface_Zero(sfc)

    ! Reset the structure components
    IF ( CRTM_SensorData_Associated(Sfc%SensorData) ) CALL CRTM_SensorData_Zero(Sfc%SensorData)

  END SUBROUTINE CRTM_Surface_Zero


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Surface_IsValid
!
! PURPOSE:
!       Non-pure function to perform some simple validity checks on a
!       CRTM Surface object.
!
!       If invalid data is found, a message is printed to stdout.
!
! CALLING SEQUENCE:
!       result = CRTM_Surface_IsValid( Sfc )
!
!         or
!
!       IF ( CRTM_Surface_IsValid( Sfc ) ) THEN....
!
! OBJECTS:
!       Sfc:       CRTM Surface object which is to have its
!                  contents checked.
!                  UNITS:      N/A
!                  TYPE:       CRTM_Surface_type
!                  DIMENSION:  Scalar
!                  ATTRIBUTES: INTENT(IN)
!
! FUNCTION RESULT:
!       result:    Logical variable indicating whether or not the input
!                  passed the check.
!                  If == .FALSE., Surface object is unused or contains
!                                 invalid data.
!                     == .TRUE.,  Surface object can be used in CRTM.
!                  UNITS:      N/A
!                  TYPE:       LOGICAL
!                  DIMENSION:  Scalar
!
!:sdoc-:
!--------------------------------------------------------------------------------

  FUNCTION CRTM_Surface_IsValid( Sfc ) RESULT( IsValid )
    TYPE(CRTM_Surface_type), INTENT(IN) :: Sfc
    LOGICAL :: IsValid
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Surface_IsValid'
    CHARACTER(ML) :: msg

    ! Check the gross surface type indicators
    IsValid = CRTM_Surface_IsCoverageValid(sfc)
    IF ( .NOT. IsValid ) THEN
      msg = 'Invalid surface coverage fraction(s) found'
      CALL Display_Message( ROUTINE_NAME, TRIM(msg), INFORMATION )
    ENDIF

    ! Check the various surface types
    IF ( Sfc%Land_Coverage  > ZERO ) IsValid = CRTM_LandSurface_IsValid(sfc)  .AND. IsValid
    IF ( Sfc%Water_Coverage > ZERO ) IsValid = CRTM_WaterSurface_IsValid(sfc) .AND. IsValid
    IF ( Sfc%Snow_Coverage  > ZERO ) IsValid = CRTM_SnowSurface_IsValid(sfc)  .AND. IsValid
    IF ( Sfc%Ice_Coverage   > ZERO ) IsValid = CRTM_IceSurface_IsValid(sfc)   .AND. IsValid

    ! Structure components
    IF ( CRTM_SensorData_Associated(Sfc%SensorData) ) &
      IsValid = CRTM_SensorData_IsValid( Sfc%SensorData ) .AND. IsValid

  END FUNCTION CRTM_Surface_IsValid



!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Surface_Inspect
!
! PURPOSE:
!       Subroutine to print the contents of a CRTM Surface object to stdout.
!
! CALLING SEQUENCE:
!       CALL CRTM_Surface_Inspect( Sfc, Unit=unit )
!
! INPUTS:
!       Sfc:   CRTM Surface object to display.
!              UNITS:      N/A
!              TYPE:       CRTM_Surface_type
!              DIMENSION:  Scalar
!              ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUTS:
!       Unit:  Unit number for an already open file to which the output
!              will be written.
!              If the argument is specified and the file unit is not
!              connected, the output goes to stdout.
!              UNITS:      N/A
!              TYPE:       INTEGER
!              DIMENSION:  Scalar
!              ATTRIBUTES: INTENT(IN), OPTIONAL
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE Scalar_Inspect( Sfc, Unit )
    ! Arguments
    TYPE(CRTM_Surface_type), INTENT(IN) :: Sfc
    INTEGER,       OPTIONAL, INTENT(IN) :: Unit
    ! Local variables
    INTEGER :: fid

    ! Setup
    fid = OUTPUT_UNIT
    IF ( PRESENT(Unit) ) THEN
      IF ( File_Open(Unit) ) fid = Unit
    END IF


    WRITE(fid,'(1x,"Surface OBJECT")')
    ! Surface coverage
    WRITE(fid,'(3x,"Land Coverage :",1x,f6.3)') Sfc%Land_Coverage
    WRITE(fid,'(3x,"Water Coverage:",1x,f6.3)') Sfc%Water_Coverage
    WRITE(fid,'(3x,"Snow Coverage :",1x,f6.3)') Sfc%Snow_Coverage
    WRITE(fid,'(3x,"Ice Coverage  :",1x,f6.3)') Sfc%Ice_Coverage
    ! The various surface types
    IF ( sfc%Land_Coverage  > ZERO ) CALL CRTM_LandSurface_Inspect(sfc, Unit=Unit)
    IF ( sfc%Water_Coverage > ZERO ) CALL CRTM_WaterSurface_Inspect(sfc, Unit=Unit)
    IF ( sfc%Snow_Coverage  > ZERO ) CALL CRTM_SnowSurface_Inspect(sfc, Unit=Unit)
    IF ( sfc%Ice_Coverage   > ZERO ) CALL CRTM_IceSurface_Inspect(sfc, Unit=Unit)
    ! SensorData information
    IF ( CRTM_SensorData_Associated(Sfc%SensorData) ) &
      CALL CRTM_SensorData_Inspect(Sfc%SensorData, Unit=Unit)

  END SUBROUTINE Scalar_Inspect

  SUBROUTINE Rank1_Inspect( Sfc, Unit )
    TYPE(CRTM_Surface_type), INTENT(IN) :: Sfc(:)
    INTEGER,          OPTIONAL, INTENT(IN) :: Unit
    INTEGER :: fid
    INTEGER :: i
    fid = OUTPUT_UNIT
    IF ( PRESENT(Unit) ) THEN
      IF ( File_Open(Unit) ) fid = Unit
    END IF
    DO i = 1, SIZE(Sfc)
      WRITE(fid, FMT='(1x,"RANK-1 INDEX:",i0," - ")', ADVANCE='NO') i
      CALL Scalar_Inspect(Sfc(i), Unit=Unit)
    END DO
  END SUBROUTINE Rank1_Inspect

  SUBROUTINE Rank2_Inspect( Sfc, Unit )
    TYPE(CRTM_Surface_type), INTENT(IN) :: Sfc(:,:)
    INTEGER,          OPTIONAL, INTENT(IN) :: Unit
    INTEGER :: fid
    INTEGER :: i, j
    fid = OUTPUT_UNIT
    IF ( PRESENT(Unit) ) THEN
      IF ( File_Open(Unit) ) fid = Unit
    END IF
    DO j = 1, SIZE(Sfc,2)
      DO i = 1, SIZE(Sfc,1)
        WRITE(fid, FMT='(1x,"RANK-2 INDEX:",i0,",",i0," - ")', ADVANCE='NO') i,j
        CALL Scalar_Inspect(Sfc(i,j), Unit=Unit)
      END DO
    END DO
  END SUBROUTINE Rank2_Inspect

!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Surface_IsCoverageValid
!
! PURPOSE:
!       Function to determine if the coverage fractions are valid
!       for a CRTM Surface object.
!
! CALLING SEQUENCE:
!       result = CRTM_Surface_IsCoverageValid( Sfc )
!
! OBJECTS:
!       Sfc:       CRTM Surface object which is to have its
!                  coverage fractions checked.
!                  UNITS:      N/A
!                  TYPE:       CRTM_Surface_type
!                  DIMENSION:  Scalar
!                  ATTRIBUTES: INTENT(IN)
!
! FUNCTION RESULT:
!       result:    Logical variable indicating whether or not the input
!                  passed the check.
!                  If == .FALSE., Surface object coverage fractions are invalid.
!                     == .TRUE.,  Surface object coverage fractions are valid.
!                  UNITS:      N/A
!                  TYPE:       LOGICAL
!                  DIMENSION:  Scalar
!
!:sdoc-:
!--------------------------------------------------------------------------------

  FUNCTION CRTM_Surface_IsCoverageValid( Sfc ) RESULT( IsValid )
    TYPE(CRTM_Surface_type), INTENT(IN) :: Sfc
    LOGICAL :: IsValid
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Surface_IsCoverageValid'
    REAL(fp)    , PARAMETER :: TOLERANCE = 1.0e-6_fp
    CHARACTER(ML) :: msg
    REAL(fp) :: Total_Coverage

    ! Compute the total coverage
    Total_Coverage = Sfc%Land_Coverage + Sfc%Water_Coverage + &
                     Sfc%Snow_Coverage + Sfc%Ice_Coverage

    ! Check coverage fractions for < 0 and > 1
    IsValid = IsCoverageValid(Sfc%Land_Coverage, 'Land')
    IsValid = IsValid .AND. IsCoverageValid(Sfc%Water_Coverage, 'Water')
    IsValid = IsValid .AND. IsCoverageValid(Sfc%Snow_Coverage, 'Snow')
    IsValid = IsValid .AND. IsCoverageValid(Sfc%Ice_Coverage, 'Ice')

    ! Check total coverage sums to 1
    IF ( ABS(Total_Coverage-ONE) > TOLERANCE ) THEN
      WRITE( msg,'("Total coverage fraction does not sum to 1 +/- ",es22.15,G12.4)' ) TOLERANCE, Total_Coverage
      CALL Display_Message( ROUTINE_NAME,msg,INFORMATION )
      IsValid = .FALSE.
    END IF

  CONTAINS

    FUNCTION IsCoverageValid( Coverage, Name ) RESULT( IsValid )
      REAL(fp)    , INTENT(IN) :: Coverage
      CHARACTER(*), INTENT(IN) :: Name
      LOGICAL :: IsValid

      IsValid = .TRUE.

      ! Check for coverage < -TOLERANCE
      IF ( Coverage < -TOLERANCE ) THEN
        WRITE( msg,'(a," coverage fraction is < ",es22.15)' ) TRIM(Name), -TOLERANCE
        CALL Display_Message( ROUTINE_NAME,msg,INFORMATION )
        IsValid = .FALSE.
      END IF

      ! Check for coverage > 1+TOLERANCE
      IF ( Coverage > ONE+TOLERANCE ) THEN
        WRITE( msg,'(a," coverage fraction is > 1 +",es22.15)' ) TRIM(Name), TOLERANCE
        CALL Display_Message( ROUTINE_NAME,msg,INFORMATION )
        IsValid = .FALSE.
      END IF

    END FUNCTION IsCoverageValid

  END FUNCTION CRTM_Surface_IsCoverageValid



!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Surface_CoverageType
!
! PURPOSE:
!       Elemental function to return the gross surface type based on coverage.
!
! CALLING SEQUENCE:
!       type = CRTM_Surface_CoverageType( sfc )
!
! INPUTS:
!       Sfc:  CRTM Surface object for which the gross surface type is required.
!             UNITS:      N/A
!             TYPE:       CRTM_Surface_type
!             DIMENSION:  Scalar or any rank
!             ATTRIBUTES: INTENT(IN)
!
! FUNCTION:
!       type: Surface type indicator for the passed CRTM Surface object.
!             UNITS:      N/A
!             TYPE:       INTEGER
!             DIMENSION:  Same as input
!
! COMMENTS:
!       For a scalar Surface object, this function result can be used to
!       determine what gross surface types are included by using it to
!       index the SURFACE_TYPE_NAME parameter arrays, e.g.
!
!         WRITE(*,*) SURFACE_TYPE_NAME(CRTM_Surface_CoverageType(sfc))
!:sdoc-:
!--------------------------------------------------------------------------------

  ELEMENTAL FUNCTION CRTM_Surface_CoverageType( sfc ) RESULT( Coverage_Type )
    TYPE(CRTM_Surface_type), INTENT(IN) :: sfc
    INTEGER :: Coverage_Type
    Coverage_Type = 0
    IF ( sfc%Land_Coverage  > ZERO ) Coverage_Type = LAND_SURFACE
    IF ( sfc%Water_Coverage > ZERO ) Coverage_Type = Coverage_Type + WATER_SURFACE
    IF ( sfc%Snow_Coverage  > ZERO ) Coverage_Type = Coverage_Type + SNOW_SURFACE
    IF ( sfc%Ice_Coverage   > ZERO ) Coverage_Type = Coverage_Type + ICE_SURFACE
  END FUNCTION CRTM_Surface_CoverageType


!------------------------------------------------------------------------------
!:sdoc+:
! NAME:
!       CRTM_Surface_Compare
!
! PURPOSE:
!       Elemental function to compare two CRTM_Surface objects to within
!       a user specified number of significant figures.
!
! CALLING SEQUENCE:
!       is_comparable = CRTM_Surface_Compare( x, y, n_SigFig=n_SigFig )
!
! OBJECTS:
!       x, y:          Two CRTM Surface objects to be compared.
!                      UNITS:      N/A
!                      TYPE:       CRTM_Surface_type
!                      DIMENSION:  Scalar or any rank
!                      ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUTS:
!       n_SigFig:      Number of significant figure to compare floating point
!                      components.
!                      UNITS:      N/A
!                      TYPE:       INTEGER
!                      DIMENSION:  Scalar or same as input
!                      ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       is_equal:      Logical value indicating whether the inputs are equal.
!                      UNITS:      N/A
!                      TYPE:       LOGICAL
!                      DIMENSION:  Same as inputs.
!:sdoc-:
!------------------------------------------------------------------------------

  ELEMENTAL FUNCTION CRTM_Surface_Compare( &
    x, &
    y, &
    n_SigFig ) &
  RESULT( is_comparable )
    TYPE(CRTM_Surface_type), INTENT(IN) :: x, y
    INTEGER,       OPTIONAL, INTENT(IN) :: n_SigFig
    LOGICAL :: is_comparable
    ! Variables
    INTEGER :: n

    ! Set up
    is_comparable = .FALSE.
    IF ( PRESENT(n_SigFig) ) THEN
      n = ABS(n_SigFig)
    ELSE
      n = DEFAULT_N_SIGFIG
    END IF

    ! Compare gross surface type coverage
    IF ( (.NOT. Compares_Within_Tolerance(x%Land_Coverage ,y%Land_Coverage ,n)) .OR. &
         (.NOT. Compares_Within_Tolerance(x%Water_Coverage,y%Water_Coverage,n)) .OR. &
         (.NOT. Compares_Within_Tolerance(x%Snow_Coverage ,y%Snow_Coverage ,n)) .OR. &
         (.NOT. Compares_Within_Tolerance(x%Ice_Coverage  ,y%Ice_Coverage  ,n)) ) RETURN

    ! Compare the land surface type data
    IF ( .NOT. CRTM_LandSurface_Compare(x,y,n_SigFig=n) ) RETURN

    ! Compare the water surface type data
    IF ( .NOT. CRTM_WaterSurface_Compare(x,y,n_SigFig=n) ) RETURN

    ! Compare the snow surface type data
    IF ( .NOT. CRTM_SnowSurface_Compare(x,y,n_SigFig=n) ) RETURN

    ! Compare the ice surface type data
    IF ( .NOT. CRTM_IceSurface_Compare(x,y,n_SigFig=n) ) RETURN

    ! Check the SensorData
    IF ( CRTM_SensorData_Associated(x%SensorData) .AND. &
         CRTM_SensorData_Associated(y%SensorData) ) THEN
      IF ( .NOT. CRTM_SensorData_Compare(x%SensorData,y%SensorData,n_SigFig=n) ) RETURN
    END IF

    ! If we get here, the structures are comparable
    is_comparable = .TRUE.

  END FUNCTION CRTM_Surface_Compare


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Surface_InquireFile
!
! PURPOSE:
!       Function to inquire CRTM Surface object files.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_Surface_InquireFile( Filename               , &
!                                                n_Channels = n_Channels, &
!                                                n_Profiles = n_Profiles  )
!
! INPUTS:
!       Filename:       Character string specifying the name of a
!                       CRTM Surface data file to read.
!                       UNITS:      N/A
!                       TYPE:       CHARACTER(*)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
! OPTIONAL OUTPUTS:
!       n_Channels:     The number of spectral channels for which there is
!                       data in the file. Note that this value will always
!                       be 0 for a profile-only dataset-- it only has meaning
!                       for K-matrix data.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: OPTIONAL, INTENT(OUT)
!
!       n_Profiles:     The number of profiles in the data file.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: OPTIONAL, INTENT(OUT)
!
! FUNCTION RESULT:
!       Error_Status:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS, the file inquire was successful
!                          == FAILURE, an unrecoverable error occurred.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_Surface_InquireFile( &
    Filename   , &  ! Input
    n_Channels , &  ! Optional output
    n_Profiles , &  ! Optional output
    NetCDF     ) &  ! Optional input
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),           INTENT(IN)  :: Filename
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Channels
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Profiles
    LOGICAL     , OPTIONAL, INTENT(IN)  :: NetCDF
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Surface_InquireFile'
    ! Function variables
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: io_msg
    INTEGER :: io_stat
    INTEGER :: fid
    INTEGER :: l, m
    LOGICAL :: binary

    ! Set up
    err_stat = SUCCESS
    ! ...Check output format
    binary = .TRUE.
    IF ( PRESENT(NetCDF) ) binary = .NOT. NetCDF
    ! ...Dispatch to the netCDF reader if requested
    IF ( .NOT. binary ) THEN
      err_stat = CRTM_Surface_InquireFile_NetCDF( Filename, &
                   n_Channels = n_Channels, &
                   n_Profiles = n_Profiles  )
      RETURN
    END IF
    ! Check that the file exists
    IF ( .NOT. File_Exists( TRIM(Filename) ) ) THEN
      msg = 'File '//TRIM(Filename)//' not found.'
      CALL Inquire_Cleanup(); RETURN
    END IF

    ! Open the file
    err_stat = Open_Binary_File( Filename, fid )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error opening '//TRIM(Filename)
      CALL Inquire_Cleanup(); RETURN
    END IF

    ! Read the number of channels,profiles
    READ( fid, IOSTAT=io_stat,IOMSG=io_msg ) l, m
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading dimensions from '//TRIM(Filename)//' - '//TRIM(io_msg)
      CALL Inquire_Cleanup(); RETURN
    END IF

    ! Close the file
    CLOSE( fid, IOSTAT=io_stat,IOMSG=io_msg )
    IF ( io_stat /= 0 ) THEN
      msg = 'Error closing '//TRIM(Filename)//' - '//TRIM(io_msg)
      CALL Inquire_Cleanup(); RETURN
    END IF

    ! Set the return arguments
    IF ( PRESENT(n_Channels) ) n_Channels = l
    IF ( PRESENT(n_Profiles) ) n_Profiles = m

  CONTAINS

    SUBROUTINE Inquire_CleanUp()
      IF ( File_Open(fid) ) THEN
        CLOSE( fid,IOSTAT=io_stat,IOMSG=io_msg )
        IF ( io_stat /= SUCCESS ) &
          msg = TRIM(msg)//'; Error closing input file during error cleanup - '//TRIM(io_msg)
      END IF
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE Inquire_CleanUp

  END FUNCTION CRTM_Surface_InquireFile


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Surface_ReadFile
!
! PURPOSE:
!       Function to read CRTM Surface object files.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_Surface_ReadFile( Filename               , &
!                                             Surface                , &
!                                             Quiet      = Quiet     , &
!                                             n_Channels = n_Channels, &
!                                             n_Profiles = n_Profiles  )
!
! INPUTS:
!       Filename:     Character string specifying the name of an
!                     Surface format data file to read.
!                     UNITS:      N/A
!                     TYPE:       CHARACTER(*)
!                     DIMENSION:  Scalar
!                     ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       Surface:      CRTM Surface object array containing the Surface
!                     data. Note the following meanings attributed to the
!                     dimensions of the object array:
!                     Rank-1: Only profile data are to be read in. The file
!                             does not contain channel information. The
!                             dimension of the structure is understood to
!                             be the PROFILE dimension.
!                     Rank-2: Channel and profile data are to be read in.
!                             The file contains both channel and profile
!                             information. The first dimension of the
!                             structure is the CHANNEL dimension, the second
!                             is the PROFILE dimension. This is to allow
!                             K-matrix structures to be read in with the
!                             same function.
!                     UNITS:      N/A
!                     TYPE:       CRTM_Surface_type
!                     DIMENSION:  Rank-1 or Rank-2
!                     ATTRIBUTES: INTENT(OUT), ALLOCATABLE
!
! OPTIONAL INPUTS:
!       Quiet:        Set this logical argument to suppress INFORMATION
!                     messages being printed to stdout
!                     If == .FALSE., INFORMATION messages are OUTPUT [DEFAULT].
!                        == .TRUE.,  INFORMATION messages are SUPPRESSED.
!                     If not specified, default is .FALSE.
!                     UNITS:      N/A
!                     TYPE:       LOGICAL
!                     DIMENSION:  Scalar
!                     ATTRIBUTES: INTENT(IN), OPTIONAL
!
! OPTIONAL OUTPUTS:
!       n_Channels:   The number of channels for which data was read. Note that
!                     this value will always be 0 for a profile-only dataset--
!                     it only has meaning for K-matrix data.
!                     UNITS:      N/A
!                     TYPE:       INTEGER
!                     DIMENSION:  Scalar
!                     ATTRIBUTES: OPTIONAL, INTENT(OUT)
!
!       n_Profiles:   The number of profiles for which data was read.
!                     UNITS:      N/A
!                     TYPE:       INTEGER
!                     DIMENSION:  Scalar
!                     ATTRIBUTES: OPTIONAL, INTENT(OUT)
!
!
! FUNCTION RESULT:
!       Error_Status: The return value is an integer defining the error status.
!                     The error codes are defined in the Message_Handler module.
!                     If == SUCCESS, the file read was successful
!                        == FAILURE, an unrecoverable error occurred.
!                     UNITS:      N/A
!                     TYPE:       INTEGER
!                     DIMENSION:  Scalar
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION Read_Surface_Rank1( &
    Filename  , &  ! Input
    Surface   , &  ! Output
    NetCDF    , &  ! Optional input
    Quiet     , &  ! Optional input
    n_Channels, &  ! Optional output
    n_Profiles, &  ! Optional output
    Debug     ) &  ! Optional input (Debug output control)
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),                         INTENT(IN)  :: Filename
    TYPE(CRTM_Surface_type), ALLOCATABLE, INTENT(OUT) :: Surface(:)  ! M
    LOGICAL,       OPTIONAL,              INTENT(IN)  :: NetCDF
    LOGICAL,       OPTIONAL,              INTENT(IN)  :: Quiet
    INTEGER,       OPTIONAL,              INTENT(OUT) :: n_Channels
    INTEGER,       OPTIONAL,              INTENT(OUT) :: n_Profiles
    LOGICAL,       OPTIONAL,              INTENT(IN)  :: Debug
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Surface_ReadFile(M)'
    ! Function variables
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: io_msg
    CHARACTER(ML) :: alloc_msg
    INTEGER :: io_stat
    INTEGER :: alloc_stat
    LOGICAL :: noisy
    LOGICAL :: binary
    INTEGER :: fid
    INTEGER :: n_input_channels
    INTEGER :: m, n_input_profiles
    INTEGER :: nch, nprof
    TYPE(CRTM_Surface_type), ALLOCATABLE :: tmp2(:,:)


    ! Set up
    err_stat = SUCCESS
    ! ...Check Quiet argument
    noisy = .TRUE.
    IF ( PRESENT(Quiet) ) noisy = .NOT. Quiet
    ! ...Override Quiet settings if debug set.
    IF ( PRESENT(Debug) ) noisy = Debug
    ! ...Profile-only (rank-1) netCDF: read the n_Channels==0 file via the rank-2
    !    reader (returns a 1 x M array) and collapse the channel axis.
    binary = .TRUE.
    IF ( PRESENT(NetCDF) ) binary = .NOT. NetCDF
    IF ( .NOT. binary ) THEN
      err_stat = Read_Surface_Rank2_NetCDF( Filename, tmp2, noisy, &
                   n_Channels = nch, n_Profiles = nprof )
      IF ( err_stat == SUCCESS ) THEN
        ! Parity with the binary rank-1 path: a profile-only file must carry the
        ! true n_Channels == 0 (stored as a global attribute; the channel
        ! dimension is forced to 1). Reject a rank-2 (K-matrix) file handed to
        ! the rank-1 reader rather than silently returning its channel-1 slice.
        IF ( nch /= 0 ) THEN
          err_stat = FAILURE
          CALL Display_Message( ROUTINE_NAME, &
            'n_Channels in '//TRIM(Filename)//' is not zero for a rank-1 '//&
            '(profiles only) Surface read.', err_stat )
          RETURN
        END IF
        Surface = tmp2(1,:)   ! auto-allocates Surface(M)
        IF ( PRESENT(n_Channels) ) n_Channels = nch
        IF ( PRESENT(n_Profiles) ) n_Profiles = nprof
      END IF
      RETURN
    END IF


    ! Open the file
    err_stat = Open_Binary_File( Filename, fid )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error opening '//TRIM(Filename)
      CALL Read_Cleanup(); RETURN
    END IF


    ! Read the dimensions
    READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) n_input_channels, n_input_profiles
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading dimensions from '//TRIM(Filename)//' - '//TRIM(io_msg)
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...Check that n_Channels is zero
    IF ( n_input_channels /= 0 ) THEN
      msg = 'n_Channels dimensions in '//TRIM(Filename)//' is not zero for a rank-1 '//&
            '(i.e. profiles only) Surface read.'
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...Allocate the return structure array
   !ALLOCATE(Surface(n_input_profiles), STAT=alloc_stat, ERRMSG=alloc_msg)
    ALLOCATE(Surface(n_input_profiles), STAT=alloc_stat)
    IF ( alloc_stat /= 0 ) THEN
      msg = 'Error allocating Surface array - '//TRIM(alloc_msg)
      CALL Read_Cleanup(); RETURN
    END IF


    ! Loop over all the profiles
    Profile_Loop: DO m = 1, n_input_profiles
      err_stat = Read_Record( fid, Surface(m), &
                              Quiet = Quiet, &
                              Debug = Debug  )
      IF ( err_stat /= SUCCESS ) THEN
        WRITE( msg,'("Error reading Surface element (",i0,") from ",a)' ) m, TRIM(Filename)
        CALL Read_Cleanup(); RETURN
      END IF
    END DO Profile_Loop


    ! Close the file
    CLOSE( fid,IOSTAT=io_stat,IOMSG=io_msg )
    IF ( io_stat /= 0 ) THEN
      msg = 'Error closing '//TRIM(Filename)//' - '//TRIM(io_msg)
      CALL Read_Cleanup(); RETURN
    END IF


    ! Set the return values
    IF ( PRESENT(n_Channels) ) n_Channels = 0
    IF ( PRESENT(n_Profiles) ) n_Profiles = n_input_profiles


    ! Output an info message
    IF ( noisy ) THEN
      WRITE( msg,'("Number of profiles read from ",a,": ",i0)' ) TRIM(Filename), n_Input_Profiles
      CALL Display_Message( ROUTINE_NAME, msg, INFORMATION )
    END IF

  CONTAINS

    SUBROUTINE Read_CleanUp()
      IF ( File_Open( Filename ) ) THEN
        CLOSE( fid,IOSTAT=io_stat,IOMSG=io_msg )
        IF ( io_stat /= 0 ) &
          msg = TRIM(msg)//'; Error closing input file during error cleanup - '//TRIM(io_msg)
      END IF
      IF ( ALLOCATED(Surface) ) THEN
       !DEALLOCATE(Surface, STAT=alloc_stat, ERRMSG=alloc_msg)
        DEALLOCATE(Surface, STAT=alloc_stat)
        IF ( alloc_stat /= 0 ) &
          msg = TRIM(msg)//'; Error deallocating Surface array during error cleanup - '//&
                TRIM(alloc_msg)
      END IF
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE Read_CleanUp

  END FUNCTION Read_Surface_Rank1


  FUNCTION Read_Surface_Rank2( &
    Filename  , &  ! Input
    Surface   , &  ! Output
    NetCDF    , &  ! Optional input
    Quiet     , &  ! Optional input
    n_Channels, &  ! Optional output
    n_Profiles, &  ! Optional output
    Debug     ) &  ! Optional input (Debug output control)
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),                         INTENT(IN)  :: Filename
    TYPE(CRTM_Surface_type), ALLOCATABLE, INTENT(OUT) :: Surface(:,:)  ! L x M
    LOGICAL,       OPTIONAL,              INTENT(IN)  :: NetCDF
    LOGICAL,       OPTIONAL,              INTENT(IN)  :: Quiet
    INTEGER,       OPTIONAL,              INTENT(OUT) :: n_Channels
    INTEGER,       OPTIONAL,              INTENT(OUT) :: n_Profiles
    LOGICAL,       OPTIONAL,              INTENT(IN)  :: Debug
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Surface_ReadFile(L x M)'
    ! Function variables
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: io_msg
    CHARACTER(ML) :: alloc_msg
    INTEGER :: io_stat
    INTEGER :: alloc_stat
    LOGICAL :: noisy
    LOGICAL :: binary
    INTEGER :: fid
    INTEGER :: l, n_input_channels
    INTEGER :: m, n_input_profiles


    ! Set up
    err_stat = SUCCESS
    ! ...Check Quiet argument
    noisy = .TRUE.
    IF ( PRESENT(Quiet) ) noisy = .NOT. Quiet
    ! ...Override Quiet settings if debug set.
    IF ( PRESENT(Debug) ) noisy = Debug
    ! ...Check output format and dispatch to the netCDF reader if requested
    binary = .TRUE.
    IF ( PRESENT(NetCDF) ) binary = .NOT. NetCDF
    IF ( .NOT. binary ) THEN
      err_stat = Read_Surface_Rank2_NetCDF( Filename, Surface, noisy, &
                   n_Channels = n_Channels, &
                   n_Profiles = n_Profiles  )
      RETURN
    END IF


    ! Open the file
    err_stat = Open_Binary_File( Filename, fid )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error opening '//TRIM(Filename)
      CALL Read_Cleanup(); RETURN
    END IF


    ! Read the dimensions
    READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) n_input_channels, n_input_profiles
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading dimensions from '//TRIM(Filename)//' - '//TRIM(io_msg)
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...Allocate the return structure array
   !ALLOCATE(Surface(n_input_channels, n_input_profiles), STAT=alloc_stat, ERRMSG=alloc_msg)
    ALLOCATE(Surface(n_input_channels, n_input_profiles), STAT=alloc_stat)
    IF ( alloc_stat /= 0 ) THEN
      msg = 'Error allocating Surface array - '//TRIM(alloc_msg)
      CALL Read_Cleanup(); RETURN
    END IF


    ! Loop over all the profiles and channels
    Profile_Loop: DO m = 1, n_input_profiles
      Channel_Loop: DO l = 1, n_input_channels
        err_stat = Read_Record( fid, Surface(l,m), &
                                Quiet = Quiet, &
                                Debug = Debug )
        IF ( err_stat /= SUCCESS ) THEN
          WRITE( msg,'("Error reading Surface element (",i0,",",i0,") from ",a)' ) &
                 l, m, TRIM(Filename)
          CALL Read_Cleanup(); RETURN
        END IF
      END DO Channel_Loop
    END DO Profile_Loop


    ! Close the file
    CLOSE( fid,IOSTAT=io_stat,IOMSG=io_msg )
    IF ( io_stat /= 0 ) THEN
      msg = 'Error closing '//TRIM(Filename)//' - '//TRIM(io_msg)
      CALL Read_Cleanup(); RETURN
    END IF


    ! Set the return values
    IF ( PRESENT(n_Channels) ) n_Channels = n_input_channels
    IF ( PRESENT(n_Profiles) ) n_Profiles = n_input_profiles


    ! Output an info message
    IF ( noisy ) THEN
      WRITE( msg,'("Number of channels and profiles read from ",a,": ",i0,1x,i0)' ) &
             TRIM(Filename), n_Input_Channels, n_Input_Profiles
      CALL Display_Message( ROUTINE_NAME, msg, INFORMATION )
    END IF

  CONTAINS

    SUBROUTINE Read_CleanUp()
      IF ( File_Open( Filename ) ) THEN
        CLOSE( fid,IOSTAT=io_stat,IOMSG=io_msg )
        IF ( io_stat /= 0 ) &
          msg = TRIM(msg)//'; Error closing input file during error cleanup - '//TRIM(io_msg)
      END IF
      IF ( ALLOCATED(Surface) ) THEN
       !DEALLOCATE(Surface, STAT=alloc_stat, ERRMSG=alloc_msg)
        DEALLOCATE(Surface, STAT=alloc_stat)
        IF ( alloc_stat /= 0 ) &
          msg = TRIM(msg)//'; Error deallocating Surface array during error cleanup - '//&
                TRIM(alloc_msg)
      END IF
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE Read_CleanUp

  END FUNCTION Read_Surface_Rank2


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Surface_WriteFile
!
! PURPOSE:
!       Function to write CRTM Surface object files.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_Surface_WriteFile( Filename     , &
!                                              Surface      , &
!                                              Quiet = Quiet  )
!
! INPUTS:
!       Filename:     Character string specifying the name of the
!                     Surface format data file to write.
!                     UNITS:      N/A
!                     TYPE:       CHARACTER(*)
!                     DIMENSION:  Scalar
!                     ATTRIBUTES: INTENT(IN)
!
!       Surface:      CRTM Surface object array containing the Surface
!                     data. Note the following meanings attributed to the
!                     dimensions of the Surface array:
!                     Rank-1: M profiles.
!                             Only profile data are to be read in. The file
!                             does not contain channel information. The
!                             dimension of the array is understood to
!                             be the PROFILE dimension.
!                     Rank-2: L channels  x  M profiles
!                             Channel and profile data are to be read in.
!                             The file contains both channel and profile
!                             information. The first dimension of the
!                             array is the CHANNEL dimension, the second
!                             is the PROFILE dimension. This is to allow
!                             K-matrix structures to be read in with the
!                             same function.
!                     UNITS:      N/A
!                     TYPE:       CRTM_Surface_type
!                     DIMENSION:  Rank-1 (M) or Rank-2 (L x M)
!                     ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUTS:
!       Quiet:        Set this logical argument to suppress INFORMATION
!                     messages being printed to stdout
!                     If == .FALSE., INFORMATION messages are OUTPUT [DEFAULT].
!                        == .TRUE.,  INFORMATION messages are SUPPRESSED.
!                     If not specified, default is .FALSE.
!                     UNITS:      N/A
!                     TYPE:       LOGICAL
!                     DIMENSION:  Scalar
!                     ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status: The return value is an integer defining the error status.
!                     The error codes are defined in the Message_Handler module.
!                     If == SUCCESS, the file write was successful
!                        == FAILURE, an unrecoverable error occurred.
!                     UNITS:      N/A
!                     TYPE:       INTEGER
!                     DIMENSION:  Scalar
!
! SIDE EFFECTS:
!       - If the output file already exists, it is overwritten.
!       - If an error occurs during *writing*, the output file is deleted before
!         returning to the calling routine.
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION Write_Surface_Rank1( &
    Filename, &  ! Input
    Surface , &  ! Input
    NetCDF  , &  ! Optional input
    Quiet   , &  ! Optional input
    Debug   ) &  ! Optional input (Debug output control)
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),            INTENT(IN) :: Filename
    TYPE(CRTM_Surface_type), INTENT(IN) :: Surface(:)  ! M
    LOGICAL,       OPTIONAL, INTENT(IN) :: NetCDF
    LOGICAL,       OPTIONAL, INTENT(IN) :: Quiet
    LOGICAL,       OPTIONAL, INTENT(IN) :: Debug
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Surface_WriteFile(M)'
    ! Function variables
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: io_msg
    LOGICAL :: noisy
    LOGICAL :: binary
    INTEGER :: io_stat
    INTEGER :: fid
    INTEGER :: m, n_Output_Profiles

    ! Setup
    err_stat = SUCCESS
    ! ...Check Quiet argument
    noisy = .TRUE.
    IF ( PRESENT(Quiet) ) noisy = .NOT. Quiet
    ! ...Override Quiet settings if debug set.
    IF ( PRESENT(Debug) ) THEN
      IF ( Debug ) noisy = .TRUE.
    END IF
    ! ...Profile-only (rank-1) netCDF: store as an n_Channels==0 file by reusing
    !    the rank-2 writer with a 1 x M view (stored n_Channels attribute = 0).
    binary = .TRUE.
    IF ( PRESENT(NetCDF) ) binary = .NOT. NetCDF
    IF ( .NOT. binary ) THEN
      err_stat = Write_Surface_Rank2_NetCDF( Filename, &
                   RESHAPE( Surface, (/ 1, SIZE(Surface) /) ), 0, noisy )
      RETURN
    END IF
    ! Dimensions
    n_Output_Profiles = SIZE(Surface)


    ! Open the file
    err_stat = Open_Binary_File( Filename, fid, For_Output = .TRUE. )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error opening '//TRIM(Filename)
      CALL Write_Cleanup(); RETURN
    END IF


    ! Write the dimensions
    WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) 0, n_Output_Profiles
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing dimensions to '//TRIM(Filename)//'- '//TRIM(io_msg)
      CALL Write_Cleanup(); RETURN
    END IF


    ! Write the data
    Profile_Loop: DO m = 1, n_Output_Profiles
      err_stat = Write_Record( fid, Surface(m), &
                               Quiet = Quiet, &
                               Debug = Debug )
      IF ( err_stat /= SUCCESS ) THEN
        WRITE( msg,'("Error writing Surface element (",i0,") to ",a)' ) m, TRIM(Filename)
        CALL Write_Cleanup(); RETURN
      END IF
    END DO Profile_Loop


    ! Close the file (if error, no delete)
    CLOSE( fid,STATUS='KEEP',IOSTAT=io_stat,IOMSG=io_msg )
    IF ( io_stat /= 0 ) THEN
      msg = 'Error closing '//TRIM(Filename)//'- '//TRIM(io_msg)
      CALL Write_Cleanup(); RETURN
    END IF


    ! Output an info message
    IF ( noisy ) THEN
      WRITE( msg,'("Number of profiles written to ",a,": ",i0)' ) &
             TRIM(Filename), n_Output_Profiles
      CALL Display_Message( ROUTINE_NAME, msg, INFORMATION )
    END IF

  CONTAINS

    SUBROUTINE Write_CleanUp()
      IF ( File_Open( Filename ) ) THEN
        CLOSE( fid,STATUS=WRITE_ERROR_STATUS,IOSTAT=io_stat,IOMSG=io_msg )
        IF ( io_stat /= 0 ) &
          msg = TRIM(msg)//'; Error deleting output file during error cleanup - '//TRIM(io_msg)
      END IF
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE Write_CleanUp

  END FUNCTION Write_Surface_Rank1


  FUNCTION Write_Surface_Rank2( &
    Filename, &  ! Input
    Surface , &  ! Input
    NetCDF  , &  ! Optional input
    Quiet   , &  ! Optional input
    Debug   ) &  ! Optional input (Debug output control)
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),            INTENT(IN)  :: Filename
    TYPE(CRTM_Surface_type), INTENT(IN)  :: Surface(:,:)  ! L x M
    LOGICAL,       OPTIONAL, INTENT(IN)  :: NetCDF
    LOGICAL,       OPTIONAL, INTENT(IN)  :: Quiet
    LOGICAL,       OPTIONAL, INTENT(IN)  :: Debug
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Surface_WriteFile(L x M)'
    ! Function variables
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: io_msg
    LOGICAL :: noisy
    LOGICAL :: binary
    INTEGER :: io_stat
    INTEGER :: fid
    INTEGER :: l, n_Output_Channels
    INTEGER :: m, n_Output_Profiles

    ! Set up
    err_stat = SUCCESS
    ! ...Check Quiet argument
    noisy = .TRUE.
    IF ( PRESENT(Quiet) ) noisy = .NOT. Quiet
    ! ...Override Quiet settings if debug set.
    IF ( PRESENT(Debug) ) THEN
      IF ( Debug ) noisy = .TRUE.
    END IF
    ! ...Check output format and dispatch to the netCDF writer if requested
    binary = .TRUE.
    IF ( PRESENT(NetCDF) ) binary = .NOT. NetCDF
    IF ( .NOT. binary ) THEN
      err_stat = Write_Surface_Rank2_NetCDF( Filename, Surface, SIZE(Surface,DIM=1), noisy )
      RETURN
    END IF
    ! Dimensions
    n_Output_Channels = SIZE(Surface,DIM=1)
    n_Output_Profiles = SIZE(Surface,DIM=2)


    ! Open the file
    err_stat = Open_Binary_File( Filename, fid, For_Output = .TRUE. )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error opening '//TRIM(Filename)
      CALL Write_Cleanup(); RETURN
    END IF


    ! Write the dimensions
    WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) n_Output_Channels, n_Output_Profiles
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing dimensions to '//TRIM(Filename)//'- '//TRIM(io_msg)
      CALL Write_Cleanup(); RETURN
    END IF


    ! Write the data
    Profile_Loop: DO m = 1, n_Output_Profiles
      Channel_Loop: DO l = 1, n_Output_Channels
        err_stat = Write_Record( fid, Surface(l,m), &
                                 Quiet = Quiet, &
                                 Debug = Debug )
        IF ( err_stat /= SUCCESS ) THEN
          WRITE( msg,'("Error writing Surface element (",i0,",",i0,") to ",a)' ) &
                 l, m, TRIM(Filename)
          CALL Write_Cleanup(); RETURN
        END IF
      END DO Channel_Loop
    END DO Profile_Loop


    ! Close the file (if error, no delete)
    CLOSE( fid,STATUS='KEEP',IOSTAT=io_stat,IOMSG=io_msg )
    IF ( io_stat /= 0 ) THEN
      msg = 'Error closing '//TRIM(Filename)//'- '//TRIM(io_msg)
      CALL Write_Cleanup(); RETURN
    END IF


    ! Output an info message
    IF ( noisy ) THEN
      WRITE( msg,'("Number of channels and profiles written to ",a,": ",i0,1x,i0 )' ) &
             TRIM(Filename), n_Output_Channels, n_Output_Profiles
      CALL Display_Message( ROUTINE_NAME, msg, INFORMATION )
    END IF

  CONTAINS

    SUBROUTINE Write_CleanUp()
      IF ( File_Open( Filename ) ) THEN
        CLOSE( fid,STATUS=WRITE_ERROR_STATUS,IOSTAT=io_stat,IOMSG=io_msg )
        IF ( io_stat /= 0 ) &
          msg = TRIM(msg)//'; Error deleting output file during error cleanup - '//TRIM(io_msg)
      END IF
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE Write_CleanUp

  END FUNCTION Write_Surface_Rank2



!##################################################################################
!##################################################################################
!##                                                                              ##
!##                          ## PRIVATE MODULE ROUTINES ##                       ##
!##                                                                              ##
!##################################################################################
!##################################################################################

!--------------------------------------------------------------------------------
!
! NAME:
!       CRTM_Surface_Equal
!
! PURPOSE:
!       Elemental function to test the equality of two CRTM_Surface objects.
!       Used in OPERATOR(==) interface block.
!
! CALLING SEQUENCE:
!       is_equal = CRTM_Surface_Equal( x, y )
!
!         or
!
!       IF ( x == y ) THEN
!         ...
!       END IF
!
! OBJECTS:
!       x, y:          Two CRTM Surface objects to be compared.
!                      UNITS:      N/A
!                      TYPE:       CRTM_Surface_type
!                      DIMENSION:  Scalar or any rank
!                      ATTRIBUTES: INTENT(IN)
!
! FUNCTION RESULT:
!       is_equal:      Logical value indicating whether the inputs are equal.
!                      UNITS:      N/A
!                      TYPE:       LOGICAL
!                      DIMENSION:  Same as inputs.
!
!--------------------------------------------------------------------------------

  ELEMENTAL FUNCTION CRTM_Surface_Equal( x, y ) RESULT( is_equal )
    TYPE(CRTM_Surface_type) , INTENT(IN)  :: x, y
    LOGICAL :: is_equal

    ! Check the gross surface type coverage
    is_equal = ( (x%Land_Coverage  .EqualTo. y%Land_Coverage ) .AND. &
                 (x%Water_Coverage .EqualTo. y%Water_Coverage) .AND. &
                 (x%Snow_Coverage  .EqualTo. y%Snow_Coverage ) .AND. &
                 (x%Ice_Coverage   .EqualTo. y%Ice_Coverage  )       )
    IF ( .NOT. is_equal ) RETURN

    ! Check the land surface type data
    is_equal = is_equal .AND. CRTM_LandSurface_Equal(x,y)
    IF ( .NOT. is_equal ) RETURN

    ! Check the water surface type data
    is_equal = is_equal .AND. CRTM_WaterSurface_Equal(x,y)
    IF ( .NOT. is_equal ) RETURN

    ! Check the snow surface type data
    is_equal = is_equal .AND. CRTM_SnowSurface_Equal(x,y)
    IF ( .NOT. is_equal ) RETURN

    ! Check the ice surface type data
    is_equal = is_equal .AND. CRTM_IceSurface_Equal(x,y)
    IF ( .NOT. is_equal ) RETURN

    ! Check the SensorData
    IF ( CRTM_SensorData_Associated(x%SensorData) .AND. &
         CRTM_SensorData_Associated(y%SensorData) ) THEN
      is_equal = is_equal .AND. (x%SensorData == y%SensorData)
    END IF

  END FUNCTION CRTM_Surface_Equal


!--------------------------------------------------------------------------------
!
! NAME:
!       CRTM_Surface_Add
!
! PURPOSE:
!       Pure function to add two CRTM Surface objects.
!       Used in OPERATOR(+) interface block.
!
! CALLING SEQUENCE:
!       sfcsum = CRTM_Surface_Add( sfc1, sfc2 )
!
!         or
!
!       sfcsum = sfc1 + sfc2
!
!
! INPUTS:
!       sfc1, sfc2: The Surface objects to add.
!                   UNITS:      N/A
!                   TYPE:       CRTM_Surface_type
!                   DIMENSION:  Scalar
!                   ATTRIBUTES: INTENT(IN OUT)
!
! RESULT:
!       sfcsum:     Surface structure containing the added components.
!                   UNITS:      N/A
!                   TYPE:       CRTM_Surface_type
!                   DIMENSION:  Scalar
!
!--------------------------------------------------------------------------------

  ELEMENTAL FUNCTION CRTM_Surface_Add( sfc1, sfc2 ) RESULT( sfcsum )
    TYPE(CRTM_Surface_type), INTENT(IN) :: sfc1, sfc2
    TYPE(CRTM_Surface_type) :: sfcsum

    ! Copy the first structure
    sfcsum = sfc1

    ! And add its components to the second one
    sfcsum%Land_Temperature      = sfcsum%Land_Temperature      + sfc2%Land_Temperature
    sfcsum%Soil_Moisture_Content = sfcsum%Soil_Moisture_Content + sfc2%Soil_Moisture_Content
    sfcsum%Canopy_Water_Content  = sfcsum%Canopy_Water_Content  + sfc2%Canopy_Water_Content
    sfcsum%Vegetation_Fraction   = sfcsum%Vegetation_Fraction   + sfc2%Vegetation_Fraction
    sfcsum%Soil_Temperature      = sfcsum%Soil_Temperature      + sfc2%Soil_Temperature
    sfcsum%LAI                   = sfcsum%LAI                   + sfc2%LAI
    sfcsum%Water_Temperature     = sfcsum%Water_Temperature     + sfc2%Water_Temperature
    sfcsum%Wind_Speed            = sfcsum%Wind_Speed            + sfc2%Wind_Speed
    sfcsum%Wind_Direction        = sfcsum%Wind_Direction        + sfc2%Wind_Direction
    sfcsum%Salinity              = sfcsum%Salinity              + sfc2%Salinity
    sfcsum%Snow_Temperature      = sfcsum%Snow_Temperature      + sfc2%Snow_Temperature
    sfcsum%Snow_Depth            = sfcsum%Snow_Depth            + sfc2%Snow_Depth
    sfcsum%Snow_Density          = sfcsum%Snow_Density          + sfc2%Snow_Density
    sfcsum%Snow_Grain_Size       = sfcsum%Snow_Grain_Size       + sfc2%Snow_Grain_Size
    sfcsum%Ice_Temperature       = sfcsum%Ice_Temperature       + sfc2%Ice_Temperature
    sfcsum%Ice_Thickness         = sfcsum%Ice_Thickness         + sfc2%Ice_Thickness
    sfcsum%Ice_Density           = sfcsum%Ice_Density           + sfc2%Ice_Density
    sfcsum%Ice_Roughness         = sfcsum%Ice_Roughness         + sfc2%Ice_Roughness
    ! ...SensorData component
    IF ( CRTM_SensorData_Associated(sfc1%SensorData) .AND. &
         CRTM_SensorData_Associated(sfc2%SensorData)       ) THEN
      sfcsum%SensorData = sfcsum%SensorData + sfc2%SensorData
    END IF

  END FUNCTION CRTM_Surface_Add


!--------------------------------------------------------------------------------
!
! NAME:
!       CRTM_Surface_Subtract
!
! PURPOSE:
!       Pure function to subtract two CRTM Surface objects.
!       Used in OPERATOR(-) interface block.
!
! CALLING SEQUENCE:
!       sfcdiff = CRTM_Surface_Subtract( sfc1, sfc2 )
!
!         or
!
!       sfcdiff = sfc1 - sfc2
!
!
! INPUTS:
!       sfc1, sfc2: The Surface objects to subtract.
!                   UNITS:      N/A
!                   TYPE:       CRTM_Surface_type
!                   DIMENSION:  Scalar
!                   ATTRIBUTES: INTENT(IN OUT)
!
! RESULT:
!       sfcdiff:    Surface structure containing the differenced components.
!                   UNITS:      N/A
!                   TYPE:       CRTM_Surface_type
!                   DIMENSION:  Scalar
!
!--------------------------------------------------------------------------------

  ELEMENTAL FUNCTION CRTM_Surface_Subtract( sfc1, sfc2 ) RESULT( sfcdiff )
    TYPE(CRTM_Surface_type), INTENT(IN) :: sfc1, sfc2
    TYPE(CRTM_Surface_type) :: sfcdiff

    ! Copy the first structure
    sfcdiff = sfc1

    ! And subtract the second one's components from it.
    sfcdiff%Land_Temperature      = sfcdiff%Land_Temperature      - sfc2%Land_Temperature
    sfcdiff%Soil_Moisture_Content = sfcdiff%Soil_Moisture_Content - sfc2%Soil_Moisture_Content
    sfcdiff%Canopy_Water_Content  = sfcdiff%Canopy_Water_Content  - sfc2%Canopy_Water_Content
    sfcdiff%Vegetation_Fraction   = sfcdiff%Vegetation_Fraction   - sfc2%Vegetation_Fraction
    sfcdiff%Soil_Temperature      = sfcdiff%Soil_Temperature      - sfc2%Soil_Temperature
    sfcdiff%LAI                   = sfcdiff%LAI                   - sfc2%LAI
    sfcdiff%Water_Temperature     = sfcdiff%Water_Temperature     - sfc2%Water_Temperature
    sfcdiff%Wind_Speed            = sfcdiff%Wind_Speed            - sfc2%Wind_Speed
    sfcdiff%Wind_Direction        = sfcdiff%Wind_Direction        - sfc2%Wind_Direction
    sfcdiff%Salinity              = sfcdiff%Salinity              - sfc2%Salinity
    sfcdiff%Snow_Temperature      = sfcdiff%Snow_Temperature      - sfc2%Snow_Temperature
    sfcdiff%Snow_Depth            = sfcdiff%Snow_Depth            - sfc2%Snow_Depth
    sfcdiff%Snow_Density          = sfcdiff%Snow_Density          - sfc2%Snow_Density
    sfcdiff%Snow_Grain_Size       = sfcdiff%Snow_Grain_Size       - sfc2%Snow_Grain_Size
    sfcdiff%Ice_Temperature       = sfcdiff%Ice_Temperature       - sfc2%Ice_Temperature
    sfcdiff%Ice_Thickness         = sfcdiff%Ice_Thickness         - sfc2%Ice_Thickness
    sfcdiff%Ice_Density           = sfcdiff%Ice_Density           - sfc2%Ice_Density
    sfcdiff%Ice_Roughness         = sfcdiff%Ice_Roughness         - sfc2%Ice_Roughness
    ! ...SensorData component
    IF ( CRTM_SensorData_Associated(sfc1%SensorData) .AND. &
         CRTM_SensorData_Associated(sfc2%SensorData)       ) THEN
      sfcdiff%SensorData = sfcdiff%SensorData - sfc2%SensorData
    END IF

  END FUNCTION CRTM_Surface_Subtract



!##################################################################################
!##################################################################################
!##                                                                              ##
!##     ## PROCEDURES BELOW WILL EVENTUALLY BE MOVED TO THEIR OWN MODULE ##      ##
!##                                                                              ##
!##################################################################################
!##################################################################################

! =============================
! LAND TYPE SPECIFIC PROCEDURES
! =============================
  ELEMENTAL SUBROUTINE CRTM_LandSurface_Zero( Sfc )
    TYPE(CRTM_Surface_type), INTENT(IN OUT) :: Sfc
    ! Zero land surface type data
    Sfc%Land_Temperature      = ZERO
    Sfc%Soil_Moisture_Content = ZERO
    Sfc%Canopy_Water_Content  = ZERO
    Sfc%Vegetation_Fraction   = ZERO
    Sfc%Soil_Temperature      = ZERO
    Sfc%LAI                   = ZERO
  END SUBROUTINE CRTM_LandSurface_Zero


  FUNCTION CRTM_LandSurface_IsValid( Sfc ) RESULT( IsValid )
    TYPE(CRTM_Surface_type), INTENT(IN) :: Sfc
    LOGICAL :: IsValid
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_LandSurface_IsValid'
    CHARACTER(ML) :: msg

    ! Setup
    IsValid = .TRUE.

    ! Check the data
    IF ( Sfc%Land_Type < 1 ) THEN
      msg = 'Invalid Land Surface type'
      CALL Display_Message( ROUTINE_NAME, TRIM(msg), INFORMATION )
      IsValid = .FALSE.
    ENDIF
    IF ( Sfc%Land_Temperature      < ZERO .OR. &
         Sfc%Soil_Moisture_Content < ZERO .OR. &
         Sfc%Canopy_Water_Content  < ZERO .OR. &
         Sfc%Vegetation_Fraction   < ZERO .OR. &
         Sfc%Soil_Temperature      < ZERO .OR. &
         Sfc%LAI                   < ZERO      ) THEN
      msg = 'Invalid Land Surface data'
      CALL Display_Message( ROUTINE_NAME, TRIM(msg), INFORMATION )
      IsValid = .FALSE.
    ENDIF

  END FUNCTION CRTM_LandSurface_IsValid


  SUBROUTINE CRTM_LandSurface_Inspect( Sfc, Unit )
    TYPE(CRTM_Surface_type), INTENT(IN) :: Sfc
    INTEGER,       OPTIONAL, INTENT(IN) :: Unit
    INTEGER :: fid
    fid = OUTPUT_UNIT
    IF ( PRESENT(Unit) ) THEN
      IF ( File_Open(Unit) ) fid = Unit
    END IF
    WRITE(fid,'(3x,"Land type index      :",1x,i0)') Sfc%Land_Type
    WRITE(fid,'(3x,"Land Temperature     :",1x,es22.15)') Sfc%Land_Temperature
    WRITE(fid,'(3x,"Soil Moisture Content:",1x,es22.15)') Sfc%Soil_Moisture_Content
    WRITE(fid,'(3x,"Canopy Water Content :",1x,es22.15)') Sfc%Canopy_Water_Content
    WRITE(fid,'(3x,"Vegetation Fraction  :",1x,es22.15)') Sfc%Vegetation_Fraction
    WRITE(fid,'(3x,"Soil Temperature     :",1x,es22.15)') Sfc%Soil_Temperature
    WRITE(fid,'(3x,"Leaf Area Index      :",1x,es22.15)') Sfc%LAI
    WRITE(fid,'(3x,"Soil type index      :",1x,i0)') Sfc%Soil_Type
    WRITE(fid,'(3x,"Vegetation type index:",1x,i0)') Sfc%Vegetation_Type
  END SUBROUTINE CRTM_LandSurface_Inspect


  ELEMENTAL FUNCTION CRTM_LandSurface_Compare( x, y, n_SigFig ) RESULT( is_comparable )
    TYPE(CRTM_Surface_type), INTENT(IN) :: x, y
    INTEGER,       OPTIONAL, INTENT(IN) :: n_SigFig
    LOGICAL :: is_comparable
    ! Variables
    INTEGER :: n

    ! Set up
    is_comparable = .FALSE.
    IF ( PRESENT(n_SigFig) ) THEN
      n = ABS(n_SigFig)
    ELSE
      n = DEFAULT_N_SIGFIG
    END IF

    ! Check integers
    IF ( x%Land_Type       /= y%Land_Type .OR. &
         x%Soil_Type       /= y%Soil_Type .OR. &
         x%Vegetation_Type /= y%Vegetation_Type ) RETURN

    ! Check floats
    IF ( (.NOT. Compares_Within_Tolerance(x%Land_Temperature     ,y%Land_Temperature     ,n)) .OR. &
         (.NOT. Compares_Within_Tolerance(x%Soil_Moisture_Content,y%Soil_Moisture_Content,n)) .OR. &
         (.NOT. Compares_Within_Tolerance(x%Canopy_Water_Content ,y%Canopy_Water_Content ,n)) .OR. &
         (.NOT. Compares_Within_Tolerance(x%Vegetation_Fraction  ,y%Vegetation_Fraction  ,n)) .OR. &
         (.NOT. Compares_Within_Tolerance(x%Soil_Temperature     ,y%Soil_Temperature     ,n)) .OR. &
         (.NOT. Compares_Within_Tolerance(x%LAI                  ,y%LAI                  ,n)) ) RETURN

    ! If we get here, the structures are comparable
    is_comparable = .TRUE.

  END FUNCTION CRTM_LandSurface_Compare


  ELEMENTAL FUNCTION CRTM_LandSurface_Equal( x, y ) RESULT( is_equal )
    TYPE(CRTM_Surface_type) , INTENT(IN)  :: x, y
    LOGICAL :: is_equal
    is_equal = ( (x%Land_Type                ==     y%Land_Type            ) .AND. &
                 (x%Land_Temperature      .EqualTo. y%Land_Temperature     ) .AND. &
                 (x%Soil_Moisture_Content .EqualTo. y%Soil_Moisture_Content) .AND. &
                 (x%Canopy_Water_Content  .EqualTo. y%Canopy_Water_Content ) .AND. &
                 (x%Vegetation_Fraction   .EqualTo. y%Vegetation_Fraction  ) .AND. &
                 (x%Soil_Temperature      .EqualTo. y%Soil_Temperature     ) .AND. &
                 (x%LAI                   .EqualTo. y%LAI                  ) .AND. &
                 (x%Soil_Type                ==     y%Soil_Type            ) .AND. &
                 (x%Vegetation_Type          ==     y%Vegetation_Type      )       )
  END FUNCTION CRTM_LandSurface_Equal


! ==============================
! WATER TYPE SPECIFIC PROCEDURES
! ==============================
  ELEMENTAL SUBROUTINE CRTM_WaterSurface_Zero( Sfc )
    TYPE(CRTM_Surface_type), INTENT(IN OUT) :: Sfc
    ! Zero the water surface type data
    Sfc%Water_Temperature = ZERO
    Sfc%Wind_Speed        = ZERO
    Sfc%Wind_Direction    = ZERO
    Sfc%Salinity          = ZERO
  END SUBROUTINE CRTM_WaterSurface_Zero


  FUNCTION CRTM_WaterSurface_IsValid( Sfc ) RESULT( IsValid )
    TYPE(CRTM_Surface_type), INTENT(IN) :: Sfc
    LOGICAL :: IsValid
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_WaterSurface_IsValid'
    CHARACTER(ML) :: msg

    ! Setup
    IsValid = .TRUE.

    ! Check the data
    IF ( Sfc%Water_Type < 1 ) THEN
      msg = 'Invalid Water Surface type'
      CALL Display_Message( ROUTINE_NAME, TRIM(msg), INFORMATION )
      IsValid = .FALSE.
    ENDIF
    IF ( Sfc%Water_Temperature < ZERO .OR. &
         Sfc%Wind_Speed        < ZERO .OR. &
         Sfc%Wind_Direction    < ZERO .OR. &
         Sfc%Salinity          < ZERO      ) THEN
      msg = 'Invalid Water Surface data'
      CALL Display_Message( ROUTINE_NAME, TRIM(msg), INFORMATION )
      IsValid = .FALSE.
    END IF

  END FUNCTION CRTM_WaterSurface_IsValid


  SUBROUTINE CRTM_WaterSurface_Inspect( Sfc, Unit )
    TYPE(CRTM_Surface_type), INTENT(IN) :: Sfc
    INTEGER,       OPTIONAL, INTENT(IN) :: Unit
    INTEGER :: fid
    fid = OUTPUT_UNIT
    IF ( PRESENT(Unit) ) THEN
      IF ( File_Open(Unit) ) fid = Unit
    END IF
    WRITE(fid,'(3x,"Water Type index :",1x,i0)') Sfc%Water_Type
    WRITE(fid,'(3x,"Water Temperature:",1x,es22.15)') Sfc%Water_Temperature
    WRITE(fid,'(3x,"Wind Speed       :",1x,es22.15)') Sfc%Wind_Speed
    WRITE(fid,'(3x,"Wind Direction   :",1x,es22.15)') Sfc%Wind_Direction
    WRITE(fid,'(3x,"Salinity         :",1x,es22.15)') Sfc%Salinity
  END SUBROUTINE CRTM_WaterSurface_Inspect


  ELEMENTAL FUNCTION CRTM_WaterSurface_Compare( x, y, n_SigFig ) RESULT( is_comparable )
    TYPE(CRTM_Surface_type), INTENT(IN) :: x, y
    INTEGER,       OPTIONAL, INTENT(IN) :: n_SigFig
    LOGICAL :: is_comparable
    ! Variables
    INTEGER :: n

    ! Set up
    is_comparable = .FALSE.
    IF ( PRESENT(n_SigFig) ) THEN
      n = ABS(n_SigFig)
    ELSE
      n = DEFAULT_N_SIGFIG
    END IF

    ! Check integers
    IF ( x%Water_Type /= y%Water_Type ) RETURN

    ! Check floats
    IF ( (.NOT. Compares_Within_Tolerance(x%Water_Temperature,y%Water_Temperature,n)) .OR. &
         (.NOT. Compares_Within_Tolerance(x%Wind_Speed       ,y%Wind_Speed       ,n)) .OR. &
         (.NOT. Compares_Within_Tolerance(x%Wind_Direction   ,y%Wind_Direction   ,n)) .OR. &
         (.NOT. Compares_Within_Tolerance(x%Salinity         ,y%Salinity         ,n)) ) RETURN

    ! If we get here, the structures are comparable
    is_comparable = .TRUE.

  END FUNCTION CRTM_WaterSurface_Compare


  ELEMENTAL FUNCTION CRTM_WaterSurface_Equal( x, y ) RESULT( is_equal )
    TYPE(CRTM_Surface_type) , INTENT(IN)  :: x, y
    LOGICAL :: is_equal
    is_equal = ( (x%Water_Type           ==     y%Water_Type       ) .AND. &
                 (x%Water_Temperature .EqualTo. y%Water_Temperature) .AND. &
                 (x%Wind_Speed        .EqualTo. y%Wind_Speed       ) .AND. &
                 (x%Wind_Direction    .EqualTo. y%Wind_Direction   ) .AND. &
                 (x%Salinity          .EqualTo. y%Salinity         )       )
  END FUNCTION CRTM_WaterSurface_Equal


! =============================
! SNOW TYPE SPECIFIC PROCEDURES
! =============================
  ELEMENTAL SUBROUTINE CRTM_SnowSurface_Zero( Sfc )
    TYPE(CRTM_Surface_type), INTENT(IN OUT) :: Sfc
    ! Zero the snow surface type data
    Sfc%Snow_Temperature = ZERO
    Sfc%Snow_Depth       = ZERO
    Sfc%Snow_Density     = ZERO
    Sfc%Snow_Grain_Size  = ZERO
  END SUBROUTINE CRTM_SnowSurface_Zero


  FUNCTION CRTM_SnowSurface_IsValid( Sfc ) RESULT( IsValid )
    TYPE(CRTM_Surface_type), INTENT(IN) :: Sfc
    LOGICAL :: IsValid
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_SnowSurface_IsValid'
    CHARACTER(ML) :: msg

    ! Setup
    IsValid = .TRUE.

    ! Check the data
    IF ( Sfc%Snow_Type < 1 ) THEN
      msg = 'Invalid Snow Surface type'
      CALL Display_Message( ROUTINE_NAME, TRIM(msg), INFORMATION )
      IsValid = .FALSE.
    ENDIF
    IF ( Sfc%Snow_Temperature < ZERO .OR. &
         Sfc%Snow_Depth       < ZERO .OR. &
         Sfc%Snow_Density     < ZERO .OR. &
         Sfc%Snow_Grain_Size  < ZERO      ) THEN
      msg = 'Invalid Snow Surface data'
      CALL Display_Message( ROUTINE_NAME, TRIM(msg), INFORMATION )
      IsValid = .FALSE.
    END IF

  END FUNCTION CRTM_SnowSurface_IsValid


  SUBROUTINE CRTM_SnowSurface_Inspect( Sfc, Unit )
    TYPE(CRTM_Surface_type), INTENT(IN) :: Sfc
    INTEGER,       OPTIONAL, INTENT(IN) :: Unit
    INTEGER :: fid
    fid = OUTPUT_UNIT
    IF ( PRESENT(Unit) ) THEN
      IF ( File_Open(Unit) ) fid = Unit
    END IF
    WRITE(fid,'(3x,"Snow Type index :",1x,i0)') Sfc%Snow_Type
    WRITE(fid,'(3x,"Snow Temperature:",1x,es22.15)') Sfc%Snow_Temperature
    WRITE(fid,'(3x,"Snow Depth      :",1x,es22.15)') Sfc%Snow_Depth
    WRITE(fid,'(3x,"Snow Density    :",1x,es22.15)') Sfc%Snow_Density
    WRITE(fid,'(3x,"Snow Grain_Size :",1x,es22.15)') Sfc%Snow_Grain_Size
  END SUBROUTINE CRTM_SnowSurface_Inspect


  ELEMENTAL FUNCTION CRTM_SnowSurface_Compare( x, y, n_SigFig ) RESULT( is_comparable )
    TYPE(CRTM_Surface_type), INTENT(IN) :: x, y
    INTEGER,       OPTIONAL, INTENT(IN) :: n_SigFig
    LOGICAL :: is_comparable
    ! Variables
    INTEGER :: n

    ! Set up
    is_comparable = .FALSE.
    IF ( PRESENT(n_SigFig) ) THEN
      n = ABS(n_SigFig)
    ELSE
      n = DEFAULT_N_SIGFIG
    END IF

    ! Check integers
    IF ( x%Snow_Type /= y%Snow_Type ) RETURN

    ! Check floats
    IF ( (.NOT. Compares_Within_Tolerance(x%Snow_Temperature,y%Snow_Temperature,n)) .OR. &
         (.NOT. Compares_Within_Tolerance(x%Snow_Depth      ,y%Snow_Depth      ,n)) .OR. &
         (.NOT. Compares_Within_Tolerance(x%Snow_Density    ,y%Snow_Density    ,n)) .OR. &
         (.NOT. Compares_Within_Tolerance(x%Snow_Grain_Size ,y%Snow_Grain_Size ,n)) ) RETURN

    ! If we get here, the structures are comparable
    is_comparable = .TRUE.

  END FUNCTION CRTM_SnowSurface_Compare


  ELEMENTAL FUNCTION CRTM_SnowSurface_Equal( x, y ) RESULT( is_equal )
    TYPE(CRTM_Surface_type) , INTENT(IN)  :: x, y
    LOGICAL :: is_equal
    is_equal = ( (x%Snow_Type           ==     y%Snow_Type       ) .AND. &
                 (x%Snow_Temperature .EqualTo. y%Snow_Temperature) .AND. &
                 (x%Snow_Depth       .EqualTo. y%Snow_Depth      ) .AND. &
                 (x%Snow_Density     .EqualTo. y%Snow_Density    ) .AND. &
                 (x%Snow_Grain_Size  .EqualTo. y%Snow_Grain_Size )       )
  END FUNCTION CRTM_SnowSurface_Equal


! ============================
! ICE TYPE SPECIFIC PROCEDURES
! ============================
  ELEMENTAL SUBROUTINE CRTM_IceSurface_Zero( Sfc )
    TYPE(CRTM_Surface_type), INTENT(IN OUT) :: Sfc
    ! Zero the ice surface type data
    Sfc%Ice_Temperature = ZERO
    Sfc%Ice_Thickness   = ZERO
    Sfc%Ice_Density     = ZERO
    Sfc%Ice_Roughness   = ZERO
  END SUBROUTINE CRTM_IceSurface_Zero


  FUNCTION CRTM_IceSurface_IsValid( Sfc ) RESULT( IsValid )
    TYPE(CRTM_Surface_type), INTENT(IN) :: Sfc
    LOGICAL :: IsValid
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_IceSurface_IsValid'
    CHARACTER(ML) :: msg

    ! Setup
    IsValid = .TRUE.

    ! Check the data
    IF ( Sfc%Ice_Type < 1 ) THEN
      msg = 'Invalid Ice Surface type'
      CALL Display_Message( ROUTINE_NAME, TRIM(msg), INFORMATION )
      IsValid = .FALSE.
    ENDIF
    IF ( Sfc%Ice_Temperature < ZERO .OR. &
         Sfc%Ice_Thickness   < ZERO .OR. &
         Sfc%Ice_Density     < ZERO .OR. &
         Sfc%Ice_Roughness   < ZERO      ) THEN
      msg = 'Invalid Ice Surface data'
      CALL Display_Message( ROUTINE_NAME, TRIM(msg), INFORMATION )
      IsValid = .FALSE.
    END IF

  END FUNCTION CRTM_IceSurface_IsValid


  SUBROUTINE CRTM_IceSurface_Inspect( Sfc, Unit )
    TYPE(CRTM_Surface_type), INTENT(IN) :: Sfc
    INTEGER,       OPTIONAL, INTENT(IN) :: Unit
    INTEGER :: fid
    fid = OUTPUT_UNIT
    IF ( PRESENT(Unit) ) THEN
      IF ( File_Open(Unit) ) fid = Unit
    END IF
    WRITE(fid,'(3x,"Ice Type index :",1x,i0)') Sfc%Ice_Type
    WRITE(fid,'(3x,"Ice Temperature:",1x,es22.15)') Sfc%Ice_Temperature
    WRITE(fid,'(3x,"Ice Thickness  :",1x,es22.15)') Sfc%Ice_Thickness
    WRITE(fid,'(3x,"Ice Density    :",1x,es22.15)') Sfc%Ice_Density
    WRITE(fid,'(3x,"Ice Roughness  :",1x,es22.15)') Sfc%Ice_Roughness
  END SUBROUTINE CRTM_IceSurface_Inspect


  ELEMENTAL FUNCTION CRTM_IceSurface_Compare( x, y, n_SigFig ) RESULT( is_comparable )
    TYPE(CRTM_Surface_type), INTENT(IN) :: x, y
    INTEGER,       OPTIONAL, INTENT(IN) :: n_SigFig
    LOGICAL :: is_comparable
    ! Variables
    INTEGER :: n

    ! Set up
    is_comparable = .FALSE.
    IF ( PRESENT(n_SigFig) ) THEN
      n = ABS(n_SigFig)
    ELSE
      n = DEFAULT_N_SIGFIG
    END IF

    ! Check integers
    IF ( x%Ice_Type /= y%Ice_Type ) RETURN

    ! Check floats
    IF ( (.NOT. Compares_Within_Tolerance(x%Ice_Temperature,y%Ice_Temperature,n)) .OR. &
         (.NOT. Compares_Within_Tolerance(x%Ice_Thickness  ,y%Ice_Thickness  ,n)) .OR. &
         (.NOT. Compares_Within_Tolerance(x%Ice_Density    ,y%Ice_Density    ,n)) .OR. &
         (.NOT. Compares_Within_Tolerance(x%Ice_Roughness  ,y%Ice_Roughness  ,n)) ) RETURN

    ! If we get here, the structures are comparable
    is_comparable = .TRUE.

  END FUNCTION CRTM_IceSurface_Compare


  ELEMENTAL FUNCTION CRTM_IceSurface_Equal( x, y ) RESULT( is_equal )
    TYPE(CRTM_Surface_type) , INTENT(IN)  :: x, y
    LOGICAL :: is_equal
    is_equal = ( (x%Ice_Type           ==     y%Ice_Type       ) .AND. &
                 (x%Ice_Temperature .EqualTo. y%Ice_Temperature) .AND. &
                 (x%Ice_Thickness   .EqualTo. y%Ice_Thickness  ) .AND. &
                 (x%Ice_Density     .EqualTo. y%Ice_Density    ) .AND. &
                 (x%Ice_Roughness   .EqualTo. y%Ice_Roughness  )       )
  END FUNCTION CRTM_IceSurface_Equal


!
! NAME:
!       Read_Record
!
! PURPOSE:
!       Utility function to read a single surface data record
!

  FUNCTION Read_Record( &
    fid        , &  ! Input
    sfc        , &  ! Output
    Quiet      , &  ! Optional input
    Debug      ) &  ! Optional input (Debug output control)
  RESULT( err_stat )
    ! Arguments
    INTEGER,                 INTENT(IN)  :: fid
    TYPE(CRTM_Surface_type), INTENT(OUT) :: sfc
    LOGICAL,       OPTIONAL, INTENT(IN)  :: Quiet
    LOGICAL,       OPTIONAL, INTENT(IN)  :: Debug
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Surface_ReadFile(Record)'
    ! Function variables
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: io_msg
    LOGICAL :: noisy
    INTEGER :: io_stat
    INTEGER :: Coverage_Type
    INTEGER :: n_Channels

    ! Set up
    err_stat = SUCCESS
    ! ...Check Quiet argument
    noisy = .TRUE.
    IF ( PRESENT(Quiet) ) noisy = .NOT. Quiet
    ! ...Override Quiet settings if debug set.
    IF ( PRESENT(Debug) ) THEN
      IF ( Debug ) noisy = .TRUE.
    END IF


    ! Read the gross surface type coverage
    READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) &
      Coverage_Type, &
      sfc%Land_Coverage, &
      sfc%Water_Coverage, &
      sfc%Snow_Coverage, &
      sfc%Ice_Coverage
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading gross surface type data - '//TRIM(io_msg)
      CALL Read_Record_Cleanup(); RETURN
    END IF
    ! ...Check the coverage fractions
    IF ( .NOT. CRTM_Surface_IsCoverageValid(sfc) ) THEN
      msg = 'Invalid surface coverage fraction(s) found'
      CALL Read_Record_Cleanup(); RETURN
    END IF
    ! ...Check the coverge surface type
    IF ( CRTM_Surface_CoverageType( sfc ) /= Coverage_Type ) THEN
      msg = 'Coverage surface type, '//&
            TRIM(SURFACE_TYPE_NAME(CRTM_Surface_CoverageType(sfc)))//&
            ', inconsistent with that specified in file.'
      CALL Read_Record_Cleanup(); RETURN
    END IF


    ! Read the surface type independent data
    READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) sfc%Wind_Speed
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading surface type independent data - '//TRIM(io_msg)
      CALL Read_Record_Cleanup(); RETURN
    END IF


    ! Read the land surface type data
    READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) &
      sfc%Land_Type, &
      sfc%Land_Temperature, &
      sfc%Soil_Moisture_Content, &
      sfc%Canopy_Water_Content , &
      sfc%Vegetation_Fraction, &
      sfc%Soil_Temperature, &
      sfc%Lai
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading land surface type data - '//TRIM(io_msg)
      CALL Read_Record_Cleanup(); RETURN
    END IF


    ! Read the water surface type data
    READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) &
      sfc%Water_Type, &
      sfc%Water_Temperature, &
      sfc%Wind_Direction, &
      sfc%Salinity
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading water surface type data - '//TRIM(io_msg)
      CALL Read_Record_Cleanup(); RETURN
    END IF


    ! Read the snow surface type data
    READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) &
      sfc%Snow_Type, &
      sfc%Snow_Temperature, &
      sfc%Snow_Depth, &
      sfc%Snow_Density, &
      sfc%Snow_Grain_Size
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading snow surface type data - '//TRIM(io_msg)
      CALL Read_Record_Cleanup(); RETURN
    END IF


    ! Read the ice surface type data
    READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) &
      sfc%Ice_Type, &
      sfc%Ice_Temperature, &
      sfc%Ice_Thickness, &
      sfc%Ice_Density, &
      sfc%Ice_Roughness
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading ice surface type data - '//TRIM(io_msg)
      CALL Read_Record_Cleanup(); RETURN
    END IF


    ! Read the SensorData
    ! ...The dimensions
    READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) n_Channels
    IF ( io_stat /= 0 ) THEN
      msg = 'Error reading SensorData dimensions - '//TRIM(io_msg)
      CALL Read_Record_Cleanup(); RETURN
    END IF
    ! ...The data
    IF ( n_Channels > 0 ) THEN
      CALL CRTM_SensorData_Create(sfc%SensorData, n_Channels )
      IF ( .NOT. CRTM_SensorData_Associated(sfc%SensorData) ) THEN
        msg = 'Error creating SensorData object.'
        CALL Read_Record_Cleanup(); RETURN
      END IF
      READ( fid,IOSTAT=io_stat,IOMSG=io_msg ) &
        sfc%SensorData%Sensor_ID       , &
        sfc%SensorData%WMO_Satellite_ID, &
        sfc%SensorData%WMO_Sensor_ID   , &
        sfc%SensorData%Sensor_Channel  , &
        sfc%SensorData%Tb
      IF ( io_stat /= 0 ) THEN
        msg = 'Error reading SensorData  - '//TRIM(io_msg)
        CALL Read_Record_Cleanup(); RETURN
      END IF
    END IF

  CONTAINS

    SUBROUTINE Read_Record_Cleanup()
      CALL CRTM_Surface_Destroy( sfc )
      CLOSE( fid,IOSTAT=io_stat,IOMSG=io_msg )
      IF ( io_stat /= SUCCESS ) &
        msg = TRIM(msg)//'; Error closing file during error cleanup - '//TRIM(io_msg)
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE Read_Record_Cleanup

  END FUNCTION Read_Record


!
! NAME:
!       Write_Record
!
! PURPOSE:
!       Utility function to write a single surface data record
!

  FUNCTION Write_Record( &
    fid  , &  ! Input
    sfc  , &  ! Input
    Quiet, &  ! Optional input
    Debug) &  ! Optional input (Debug output control)
  RESULT( err_stat )
    ! Arguments
    INTEGER,                 INTENT(IN) :: fid
    TYPE(CRTM_Surface_type), INTENT(IN) :: sfc
    LOGICAL,       OPTIONAL, INTENT(IN) :: Quiet
    LOGICAL,       OPTIONAL, INTENT(IN) :: Debug
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Surface_WriteFile(Record)'
    ! Function variables
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: io_msg
    LOGICAL :: noisy
    INTEGER :: io_stat


    ! Set up
    err_stat = SUCCESS
    ! ...Check Quiet argument
    noisy = .TRUE.
    IF ( PRESENT(Quiet) ) noisy = .NOT. Quiet
    ! ...Override Quiet settings if debug set.
    IF ( PRESENT(Debug) ) THEN
      IF ( Debug ) noisy = .TRUE.
    END IF


    ! Write the gross surface type coverage
    WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) &
      CRTM_Surface_CoverageType(sfc), &
      sfc%Land_Coverage, &
      sfc%Water_Coverage, &
      sfc%Snow_Coverage, &
      sfc%Ice_Coverage
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing gross surface type data - '//TRIM(io_msg)
      CALL Write_Record_Cleanup(); RETURN
    END IF


    ! Write the surface type independent data
    WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) sfc%Wind_Speed
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing surface type independent data - '//TRIM(io_msg)
      CALL Write_Record_Cleanup(); RETURN
    END IF


    ! Write the land surface type data
    WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) &
      sfc%Land_Type, &
      sfc%Land_Temperature, &
      sfc%Soil_Moisture_Content, &
      sfc%Canopy_Water_Content, &
      sfc%Vegetation_Fraction, &
      sfc%Soil_Temperature, &
      sfc%Lai
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing land surface type data - '//TRIM(io_msg)
      CALL Write_Record_Cleanup(); RETURN
    END IF


    ! Write the water surface type data
    WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) &
      sfc%Water_Type, &
      sfc%Water_Temperature, &
      sfc%Wind_Direction, &
      sfc%Salinity
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing water surface type data - '//TRIM(io_msg)
      CALL Write_Record_Cleanup(); RETURN
    END IF


    ! Write the snow surface type data
    WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) &
      sfc%Snow_Type, &
      sfc%Snow_Temperature, &
      sfc%Snow_Depth, &
      sfc%Snow_Density, &
      sfc%Snow_Grain_Size
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing snow surface type data - '//TRIM(io_msg)
      CALL Write_Record_Cleanup(); RETURN
    END IF


    ! Write the ice surface type data
    WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) &
      sfc%Ice_Type, &
      sfc%Ice_Temperature, &
      sfc%Ice_Thickness, &
      sfc%Ice_Density, &
      sfc%Ice_Roughness
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing ice surface type data - '//TRIM(io_msg)
      CALL Write_Record_Cleanup(); RETURN
    END IF


    ! Write the SensorData object
    ! ...The dimensions
    WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) sfc%SensorData%n_Channels
    IF ( io_stat /= 0 ) THEN
      msg = 'Error writing SensorData dimensions - '//TRIM(io_msg)
      CALL Write_Record_Cleanup(); RETURN
    END IF
    ! ...The data
    IF ( sfc%SensorData%n_Channels > 0 ) THEN
      WRITE( fid,IOSTAT=io_stat,IOMSG=io_msg ) &
        sfc%SensorData%Sensor_ID       , &
        sfc%SensorData%WMO_Satellite_ID, &
        sfc%SensorData%WMO_Sensor_ID   , &
        sfc%SensorData%Sensor_Channel  , &
        sfc%SensorData%Tb
      IF ( io_stat /= 0 ) THEN
        msg = 'Error writing SensorData - '//TRIM(io_msg)
        CALL Write_Record_Cleanup(); RETURN
      END IF
    END IF

  CONTAINS

    SUBROUTINE Write_Record_Cleanup()
      CLOSE( fid,STATUS=WRITE_ERROR_STATUS,IOSTAT=io_stat,IOMSG=io_msg )
      IF ( io_stat /= SUCCESS ) &
        msg = TRIM(msg)//'; Error closing file during error cleanup - '//TRIM(io_msg)
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE Write_Record_Cleanup

  END FUNCTION Write_Record


!##############################################################################
!##############################################################################
!##                                                                          ##
!##                       ## netCDF I/O WORKER ROUTINES ##                    ##
!##                                                                          ##
!##############################################################################
!##############################################################################

!------------------------------------------------------------------------------
!
! NAME:
!       CRTM_Surface_InquireFile_NetCDF
!
! PURPOSE:
!       Function to inquire the dimensions of a netCDF CRTM Surface file.
!       n_Channels is returned as 0 if the file has no channel dimension
!       (i.e. a profile-only dataset).
!
!------------------------------------------------------------------------------

  FUNCTION CRTM_Surface_InquireFile_NetCDF( &
    Filename   , &  ! Input
    n_Channels , &  ! Optional output
    n_Profiles ) &  ! Optional output
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),           INTENT(IN)  :: Filename
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Channels
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Profiles
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Surface_InquireFile_NetCDF'
    ! Function variables
    CHARACTER(ML) :: msg
    LOGICAL :: Close_File
    INTEGER :: NF90_Status
    INTEGER :: FileId, DimId
    INTEGER :: l, m

    ! Set up
    err_stat = SUCCESS
    Close_File = .FALSE.
    ! ...Check that the file exists
    IF ( .NOT. File_Exists( TRIM(Filename) ) ) THEN
      msg = 'File '//TRIM(Filename)//' not found.'
      CALL Inquire_Cleanup(); RETURN
    END IF

    ! Open the file
    NF90_Status = NF90_OPEN( Filename,NF90_NOWRITE,FileId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error opening '//TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_Cleanup(); RETURN
    END IF
    ! ...Close the file if any error from here on
    Close_File = .TRUE.

    ! Get the number of profiles (always present)
    NF90_Status = NF90_INQ_DIMID( FileId,SFC_PROFILE_DIMNAME,DimId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring dimension ID for '//SFC_PROFILE_DIMNAME//' - '// &
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_INQUIRE_DIMENSION( FileId,DimId,Len=m )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error reading dimension value for '//SFC_PROFILE_DIMNAME//' - '// &
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_Cleanup(); RETURN
    END IF

    ! Get the true number of channels from the global attribute (0 for a
    ! profile-only/rank-1 Surface; the channel dimension is MAX(n_Channels,1))
    NF90_Status = NF90_GET_ATT( FileId,NF90_GLOBAL,SFC_NCHANNELS_GATTNAME,l )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error reading global attribute '//SFC_NCHANNELS_GATTNAME//' - '// &
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_Cleanup(); RETURN
    END IF

    ! Close the file
    NF90_Status = NF90_CLOSE( FileId ); Close_File = .FALSE.
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error closing '//TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_Cleanup(); RETURN
    END IF

    ! Set the return arguments
    IF ( PRESENT(n_Channels) ) n_Channels = l
    IF ( PRESENT(n_Profiles) ) n_Profiles = m

  CONTAINS

    SUBROUTINE Inquire_CleanUp()
      IF ( Close_File ) THEN
        NF90_Status = NF90_CLOSE( FileId )
        IF ( NF90_Status /= NF90_NOERR ) &
          msg = TRIM(msg)//'; Error closing input file during error cleanup - '//&
                TRIM(NF90_STRERROR( NF90_Status ))
      END IF
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE Inquire_CleanUp

  END FUNCTION CRTM_Surface_InquireFile_NetCDF


!------------------------------------------------------------------------------
!
! NAME:
!       CreateFile_Surface_netCDF
!
! PURPOSE:
!       Utility function to create a netCDF Surface file: defines the
!       dimensions and the single packed Surface_Data variable, leaving the
!       file open (out of define mode) for the caller to populate.
!
!------------------------------------------------------------------------------

  FUNCTION CreateFile_Surface_netCDF( &
    Filename  , &  ! Input
    n_Channels, &  ! Input
    n_Profiles, &  ! Input
    FileId    ) &  ! Output
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*), INTENT(IN)  :: Filename
    INTEGER     , INTENT(IN)  :: n_Channels
    INTEGER     , INTENT(IN)  :: n_Profiles
    INTEGER     , INTENT(OUT) :: FileId
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Surface_WriteFile(netCDF)'
    ! Local variables
    CHARACTER(ML) :: msg
    LOGICAL :: Close_File
    INTEGER :: NF90_Status
    INTEGER :: n_Channels_DimID
    INTEGER :: n_Profiles_DimID
    INTEGER :: n_Fields_DimID
    INTEGER :: VarID

    ! Setup
    err_stat = SUCCESS
    Close_File = .FALSE.

    ! Create the data file
    NF90_Status = NF90_CREATE( Filename,NF90_CLOBBER,FileId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error creating '//TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Create_Cleanup(); RETURN
    END IF
    ! ...Close the file if any error from here on
    Close_File = .TRUE.

    ! Define the dimensions (channel dim is MAX(n_Channels,1); see GATT below)
    NF90_Status = NF90_DEF_DIM( FileID,SFC_CHANNEL_DIMNAME,MAX(n_Channels,1),n_Channels_DimID )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error defining '//SFC_CHANNEL_DIMNAME//' dimension in '//&
            TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Create_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_DEF_DIM( FileID,SFC_PROFILE_DIMNAME,n_Profiles,n_Profiles_DimID )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error defining '//SFC_PROFILE_DIMNAME//' dimension in '//&
            TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Create_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_DEF_DIM( FileID,SFC_FIELD_DIMNAME,N_SURFACE_FIELDS,n_Fields_DimID )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error defining '//SFC_FIELD_DIMNAME//' dimension in '//&
            TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Create_Cleanup(); RETURN
    END IF

    ! Write the true n_Channels (0 for a profile-only/rank-1 Surface)
    NF90_Status = NF90_PUT_ATT( FileId,NF90_GLOBAL,SFC_NCHANNELS_GATTNAME,n_Channels )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error setting '//SFC_NCHANNELS_GATTNAME//' attribute in '//&
            TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Create_Cleanup(); RETURN
    END IF

    ! Define the packed data variable
    NF90_Status = NF90_DEF_VAR( FileID, &
      SFC_DATA_VARNAME, &
      SFC_FLOAT_TYPE, &
      dimIDs=(/n_Channels_DimID, n_Profiles_DimID, n_Fields_DimID/), &
      varID=VarID )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error defining '//SFC_DATA_VARNAME//' variable in '//&
            TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Create_Cleanup(); RETURN
    END IF

    ! Take the file out of define mode
    NF90_Status = NF90_ENDDEF( FileId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error taking file '//TRIM(Filename)// &
            ' out of define mode - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Create_Cleanup(); RETURN
    END IF

  CONTAINS

    SUBROUTINE Create_CleanUp()
      IF ( Close_File ) THEN
        NF90_Status = NF90_CLOSE( FileID )
        IF ( NF90_Status /= NF90_NOERR ) &
          msg = TRIM(msg)//'; Error closing file during error cleanup - '//&
                TRIM(NF90_STRERROR( NF90_Status ))
      END IF
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME,msg,err_stat )
    END SUBROUTINE Create_CleanUp

  END FUNCTION CreateFile_Surface_netCDF


!------------------------------------------------------------------------------
!
! NAME:
!       Write_Surface_Rank2_NetCDF
!
! PURPOSE:
!       Utility function to write a rank-2 (L x M) Surface array to a netCDF
!       file. Populated SensorData (n_Channels > 0) is not supported and is
!       rejected (no driver/baseline populates it).
!
!------------------------------------------------------------------------------

  FUNCTION Write_Surface_Rank2_NetCDF( &
    Filename         , &  ! Input
    Surface          , &  ! Input
    n_Channels_stored, &  ! Input (true n_Channels; 0 for profile-only rank-1)
    noisy            ) &  ! Input
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),            INTENT(IN) :: Filename
    TYPE(CRTM_Surface_type), INTENT(IN) :: Surface(:,:)  ! L x M
    INTEGER,                 INTENT(IN) :: n_Channels_stored
    LOGICAL,                 INTENT(IN) :: noisy
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Surface_WriteFile_netCDF'
    ! Function variables
    CHARACTER(ML) :: msg
    LOGICAL :: Close_File
    INTEGER :: NF90_Status, FileId, VarId
    INTEGER :: l, m, n_Channels, n_Profiles, alloc_stat
    REAL(fp), ALLOCATABLE :: Surface_Data(:,:,:)

    ! Set up
    err_stat = SUCCESS
    Close_File = .FALSE.
    n_Channels = SIZE(Surface,DIM=1)
    n_Profiles = SIZE(Surface,DIM=2)

    ! Reject populated SensorData (unsupported by the packed schema)
    IF ( ANY( Surface%SensorData%n_Channels > 0 ) ) THEN
      msg = 'Populated SensorData (n_Channels > 0) is not supported by the '//&
            'Surface netCDF writer.'
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
      RETURN
    END IF

    ! Pack the per-element data
    ALLOCATE( Surface_Data( n_Channels, n_Profiles, N_SURFACE_FIELDS ), STAT=alloc_stat )
    IF ( alloc_stat /= 0 ) THEN
      msg = 'Error allocating Surface_Data array'
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
      RETURN
    END IF
    DO m = 1, n_Profiles
      DO l = 1, n_Channels
        Surface_Data(l,m,IDX_LAND_COVERAGE)         = Surface(l,m)%Land_Coverage
        Surface_Data(l,m,IDX_WATER_COVERAGE)        = Surface(l,m)%Water_Coverage
        Surface_Data(l,m,IDX_SNOW_COVERAGE)         = Surface(l,m)%Snow_Coverage
        Surface_Data(l,m,IDX_ICE_COVERAGE)          = Surface(l,m)%Ice_Coverage
        Surface_Data(l,m,IDX_WIND_SPEED)            = Surface(l,m)%Wind_Speed
        Surface_Data(l,m,IDX_LAND_TEMPERATURE)      = Surface(l,m)%Land_Temperature
        Surface_Data(l,m,IDX_SOIL_MOISTURE_CONTENT) = Surface(l,m)%Soil_Moisture_Content
        Surface_Data(l,m,IDX_CANOPY_WATER_CONTENT)  = Surface(l,m)%Canopy_Water_Content
        Surface_Data(l,m,IDX_VEGETATION_FRACTION)   = Surface(l,m)%Vegetation_Fraction
        Surface_Data(l,m,IDX_SOIL_TEMPERATURE)      = Surface(l,m)%Soil_Temperature
        Surface_Data(l,m,IDX_LAI)                   = Surface(l,m)%LAI
        Surface_Data(l,m,IDX_WATER_TEMPERATURE)     = Surface(l,m)%Water_Temperature
        Surface_Data(l,m,IDX_WIND_DIRECTION)        = Surface(l,m)%Wind_Direction
        Surface_Data(l,m,IDX_SALINITY)              = Surface(l,m)%Salinity
        Surface_Data(l,m,IDX_SNOW_TEMPERATURE)      = Surface(l,m)%Snow_Temperature
        Surface_Data(l,m,IDX_SNOW_DEPTH)            = Surface(l,m)%Snow_Depth
        Surface_Data(l,m,IDX_SNOW_DENSITY)          = Surface(l,m)%Snow_Density
        Surface_Data(l,m,IDX_SNOW_GRAIN_SIZE)       = Surface(l,m)%Snow_Grain_Size
        Surface_Data(l,m,IDX_ICE_TEMPERATURE)       = Surface(l,m)%Ice_Temperature
        Surface_Data(l,m,IDX_ICE_THICKNESS)         = Surface(l,m)%Ice_Thickness
        Surface_Data(l,m,IDX_ICE_DENSITY)           = Surface(l,m)%Ice_Density
        Surface_Data(l,m,IDX_ICE_ROUGHNESS)         = Surface(l,m)%Ice_Roughness
        Surface_Data(l,m,IDX_LAND_TYPE)             = REAL(Surface(l,m)%Land_Type      , fp)
        Surface_Data(l,m,IDX_SOIL_TYPE)             = REAL(Surface(l,m)%Soil_Type      , fp)
        Surface_Data(l,m,IDX_VEGETATION_TYPE)       = REAL(Surface(l,m)%Vegetation_Type, fp)
        Surface_Data(l,m,IDX_WATER_TYPE)            = REAL(Surface(l,m)%Water_Type     , fp)
        Surface_Data(l,m,IDX_SNOW_TYPE)             = REAL(Surface(l,m)%Snow_Type      , fp)
        Surface_Data(l,m,IDX_ICE_TYPE)              = REAL(Surface(l,m)%Ice_Type       , fp)
        Surface_Data(l,m,IDX_SENSORDATA_N_CHANNELS) = REAL(Surface(l,m)%SensorData%n_Channels, fp)
      END DO
    END DO

    ! Create the output file (defines dims + variable). The stored n_Channels
    ! is the true value (0 for rank-1); the data array uses SIZE = MAX(.,1).
    err_stat = CreateFile_Surface_netCDF( Filename, n_Channels_stored, n_Profiles, FileId )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error creating output file '//TRIM(Filename)
      CALL Write_Cleanup(); RETURN
    END IF
    ! ...Close the file if any error from here on
    Close_File = .TRUE.

    ! Write the packed data
    NF90_Status = NF90_INQ_VARID( FileId,SFC_DATA_VARNAME,VarId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring '//TRIM(Filename)//' for '//SFC_DATA_VARNAME//&
            ' variable ID - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Write_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_PUT_VAR( FileId,VarId,Surface_Data )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error writing '//SFC_DATA_VARNAME//' to '//TRIM(Filename)//&
            ' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Write_Cleanup(); RETURN
    END IF

    ! Close the file
    NF90_Status = NF90_CLOSE( FileId ); Close_File = .FALSE.
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error closing output file - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Write_Cleanup(); RETURN
    END IF

    ! Output an info message
    IF ( noisy ) THEN
      WRITE( msg,'("Number of channels and profiles written to ",a,": ",i0,1x,i0 )' ) &
             TRIM(Filename), n_Channels, n_Profiles
      CALL Display_Message( ROUTINE_NAME, msg, INFORMATION )
    END IF

  CONTAINS

    SUBROUTINE Write_CleanUp()
      IF ( Close_File ) THEN
        NF90_Status = NF90_CLOSE( FileId )
        IF ( NF90_Status /= NF90_NOERR ) &
          msg = TRIM(msg)//'; Error closing output file during error cleanup - '//&
                TRIM(NF90_STRERROR( NF90_Status ))
      END IF
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE Write_CleanUp

  END FUNCTION Write_Surface_Rank2_NetCDF


!------------------------------------------------------------------------------
!
! NAME:
!       Read_Surface_Rank2_NetCDF
!
! PURPOSE:
!       Utility function to read a rank-2 (L x M) Surface array from a netCDF
!       file written by Write_Surface_Rank2_NetCDF.
!
!------------------------------------------------------------------------------

  FUNCTION Read_Surface_Rank2_NetCDF( &
    Filename  , &  ! Input
    Surface   , &  ! Output
    noisy     , &  ! Input
    n_Channels, &  ! Optional output
    n_Profiles) &  ! Optional output
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),                         INTENT(IN)  :: Filename
    TYPE(CRTM_Surface_type), ALLOCATABLE, INTENT(OUT) :: Surface(:,:)  ! L x M
    LOGICAL,                              INTENT(IN)  :: noisy
    INTEGER,       OPTIONAL,              INTENT(OUT) :: n_Channels
    INTEGER,       OPTIONAL,              INTENT(OUT) :: n_Profiles
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Surface_ReadFile_netCDF'
    ! Function variables
    CHARACTER(ML) :: msg
    LOGICAL :: Close_File
    INTEGER :: NF90_Status, FileId, VarId
    INTEGER :: l, m, n_File_Channels, n_File_Profiles, n_Buf_Channels, alloc_stat
    REAL(fp), ALLOCATABLE :: Surface_Data(:,:,:)

    ! Set up
    err_stat = SUCCESS
    Close_File = .FALSE.
    ! ...Check that the file exists
    IF ( .NOT. File_Exists( TRIM(Filename) ) ) THEN
      msg = 'File '//TRIM(Filename)//' not found.'
      CALL Read_Cleanup(); RETURN
    END IF

    ! Inquire the file for its dimensions. n_File_Channels is the true value
    ! (0 for a profile-only/rank-1 file); the stored array has MAX(.,1) channels.
    err_stat = CRTM_Surface_InquireFile_NetCDF( Filename, &
                 n_Channels = n_File_Channels, &
                 n_Profiles = n_File_Profiles  )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error obtaining Surface dimensions from '//TRIM(Filename)
      CALL Read_Cleanup(); RETURN
    END IF
    n_Buf_Channels = MAX( n_File_Channels, 1 )

    ! Allocate the return structure and the read buffer
    ALLOCATE( Surface( n_Buf_Channels, n_File_Profiles ), &
              Surface_Data( n_Buf_Channels, n_File_Profiles, N_SURFACE_FIELDS ), &
              STAT = alloc_stat )
    IF ( alloc_stat /= 0 ) THEN
      msg = 'Error allocating Surface/Surface_Data arrays'
      CALL Read_Cleanup(); RETURN
    END IF

    ! Open the file for reading
    NF90_Status = NF90_OPEN( Filename,NF90_NOWRITE,FileId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error opening '//TRIM(Filename)//' for read access - '//&
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...Close the file if any error from here on
    Close_File = .TRUE.

    ! Read the packed data
    NF90_Status = NF90_INQ_VARID( FileId,SFC_DATA_VARNAME,VarId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring '//TRIM(Filename)//' for '//SFC_DATA_VARNAME//&
            ' variable ID - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_GET_VAR( FileId,VarId,Surface_Data )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error reading '//SFC_DATA_VARNAME//' from '//TRIM(Filename)//&
            ' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF

    ! Close the file
    NF90_Status = NF90_CLOSE( FileId ); Close_File = .FALSE.
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error closing input file - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF

    ! Unpack into the return structure (INTEGER fields recovered with NINT)
    DO m = 1, n_File_Profiles
      DO l = 1, n_Buf_Channels
        Surface(l,m)%Land_Coverage         = Surface_Data(l,m,IDX_LAND_COVERAGE)
        Surface(l,m)%Water_Coverage        = Surface_Data(l,m,IDX_WATER_COVERAGE)
        Surface(l,m)%Snow_Coverage         = Surface_Data(l,m,IDX_SNOW_COVERAGE)
        Surface(l,m)%Ice_Coverage          = Surface_Data(l,m,IDX_ICE_COVERAGE)
        Surface(l,m)%Wind_Speed            = Surface_Data(l,m,IDX_WIND_SPEED)
        Surface(l,m)%Land_Temperature      = Surface_Data(l,m,IDX_LAND_TEMPERATURE)
        Surface(l,m)%Soil_Moisture_Content = Surface_Data(l,m,IDX_SOIL_MOISTURE_CONTENT)
        Surface(l,m)%Canopy_Water_Content  = Surface_Data(l,m,IDX_CANOPY_WATER_CONTENT)
        Surface(l,m)%Vegetation_Fraction   = Surface_Data(l,m,IDX_VEGETATION_FRACTION)
        Surface(l,m)%Soil_Temperature      = Surface_Data(l,m,IDX_SOIL_TEMPERATURE)
        Surface(l,m)%LAI                   = Surface_Data(l,m,IDX_LAI)
        Surface(l,m)%Water_Temperature     = Surface_Data(l,m,IDX_WATER_TEMPERATURE)
        Surface(l,m)%Wind_Direction        = Surface_Data(l,m,IDX_WIND_DIRECTION)
        Surface(l,m)%Salinity              = Surface_Data(l,m,IDX_SALINITY)
        Surface(l,m)%Snow_Temperature      = Surface_Data(l,m,IDX_SNOW_TEMPERATURE)
        Surface(l,m)%Snow_Depth            = Surface_Data(l,m,IDX_SNOW_DEPTH)
        Surface(l,m)%Snow_Density          = Surface_Data(l,m,IDX_SNOW_DENSITY)
        Surface(l,m)%Snow_Grain_Size       = Surface_Data(l,m,IDX_SNOW_GRAIN_SIZE)
        Surface(l,m)%Ice_Temperature       = Surface_Data(l,m,IDX_ICE_TEMPERATURE)
        Surface(l,m)%Ice_Thickness         = Surface_Data(l,m,IDX_ICE_THICKNESS)
        Surface(l,m)%Ice_Density           = Surface_Data(l,m,IDX_ICE_DENSITY)
        Surface(l,m)%Ice_Roughness         = Surface_Data(l,m,IDX_ICE_ROUGHNESS)
        Surface(l,m)%Land_Type             = NINT(Surface_Data(l,m,IDX_LAND_TYPE))
        Surface(l,m)%Soil_Type             = NINT(Surface_Data(l,m,IDX_SOIL_TYPE))
        Surface(l,m)%Vegetation_Type       = NINT(Surface_Data(l,m,IDX_VEGETATION_TYPE))
        Surface(l,m)%Water_Type            = NINT(Surface_Data(l,m,IDX_WATER_TYPE))
        Surface(l,m)%Snow_Type             = NINT(Surface_Data(l,m,IDX_SNOW_TYPE))
        Surface(l,m)%Ice_Type              = NINT(Surface_Data(l,m,IDX_ICE_TYPE))
        ! SensorData is guaranteed unpopulated (n_Channels == 0); leave default.
      END DO
    END DO

    ! Set the return values
    IF ( PRESENT(n_Channels) ) n_Channels = n_File_Channels
    IF ( PRESENT(n_Profiles) ) n_Profiles = n_File_Profiles

    ! Output an info message
    IF ( noisy ) THEN
      WRITE( msg,'("Number of channels and profiles read from ",a,": ",i0,1x,i0)' ) &
             TRIM(Filename), n_File_Channels, n_File_Profiles
      CALL Display_Message( ROUTINE_NAME, msg, INFORMATION )
    END IF

  CONTAINS

    SUBROUTINE Read_CleanUp()
      IF ( Close_File ) THEN
        NF90_Status = NF90_CLOSE( FileId )
        IF ( NF90_Status /= NF90_NOERR ) &
          msg = TRIM(msg)//'; Error closing input file during error cleanup - '//&
                TRIM(NF90_STRERROR( NF90_Status ))
      END IF
      IF ( ALLOCATED(Surface) ) DEALLOCATE(Surface, STAT=alloc_stat)
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE Read_CleanUp

  END FUNCTION Read_Surface_Rank2_NetCDF

END MODULE CRTM_Surface_Define
