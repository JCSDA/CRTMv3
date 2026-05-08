!
! PARMIOCoeff_netCDF_IO
!
! Reader for PARMIO LUT netCDF coefficient files produced by
! parmio/scripts/write_parmio_lut_netcdf.py. Loads the three coefficient
! groups (sss_dependent, sss_nominal_m, sss_nominal_h) into a
! PARMIOCoeff_type instance.

MODULE PARMIOCoeff_netCDF_IO

  USE Type_Kinds      , ONLY: fp, Long, Double
  USE Message_Handler , ONLY: SUCCESS, FAILURE, INFORMATION, Display_Message
  USE File_Utility    , ONLY: File_Exists
  USE PARMIOCoeff_Define, ONLY: &
    PARMIOCoeff_type,           &
    PARMIOCoeff_Group_type,     &
    PARMIOCoeff_RC_Group_type,  &
    PARMIOCoeff_Destroy,        &
    PARMIOCoeff_Create,         &
    PARMIOCoeff_Associated,     &
    PARMIOCoeff_Info,           &
    N_PARMIO_HARMONIC_TERMS,    &
    PARMIO_N_RC_POLARIZATIONS,  &
    PARMIO_RC_V_POL,            &
    PARMIO_RC_H_POL,            &
    PARMIO_GROUP_SSS_DEPENDENT, &
    PARMIO_GROUP_SSS_NOMINAL_M, &
    PARMIO_GROUP_SSS_NOMINAL_H, &
    PARMIO_N_GROUPS
  USE netcdf
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: PARMIOCoeff_netCDF_ReadFile

  ! Group names as written by the Python LUT generator
  CHARACTER(*), PARAMETER :: GROUP_NAME(PARMIO_N_GROUPS) = (/ &
    'sss_dependent  ', 'sss_nominal_m  ', 'sss_nominal_h  ' /)

  ! Variable / dimension names
  CHARACTER(*), PARAMETER :: DIM_FREQUENCY  = 'frequency'
  CHARACTER(*), PARAMETER :: DIM_THETA      = 'theta'
  CHARACTER(*), PARAMETER :: DIM_WIND_SPEED = 'wind_speed'
  CHARACTER(*), PARAMETER :: DIM_SST        = 'sst'
  CHARACTER(*), PARAMETER :: DIM_SSS        = 'sss'
  CHARACTER(*), PARAMETER :: DIM_FOAM_STATE = 'foam_state'
  CHARACTER(*), PARAMETER :: DIM_TRANSMITTANCE = 'transmittance'
  CHARACTER(*), PARAMETER :: VAR_FOAM       = 'Foam'
  CHARACTER(*), PARAMETER :: VAR_RDOWN_V    = 'Rdown_v'
  CHARACTER(*), PARAMETER :: VAR_RDOWN_H    = 'Rdown_h'

  ! 14 harmonic-coefficient variable names — order MUST match the harmonic
  ! index convention in PARMIOCoeff_Define.
  CHARACTER(8), PARAMETER :: HARMONIC_VARNAMES(N_PARMIO_HARMONIC_TERMS) = (/ &
    'evN     ', 'ehN     ',                                                  &
    'ev0     ', 'eh0     ',                                                  &
    'ev1     ', 'eh1     ', 'eU1     ', 'eV1     ',                          &
    'ev2     ', 'eh2     ', 'eU2     ', 'eV2     ',                          &
    'edv_MR  ', 'edh_MR  ' /)

  INTEGER, PARAMETER :: ML = 256

CONTAINS

  FUNCTION PARMIOCoeff_netCDF_ReadFile( &
      PARMIOCoeff, Filename, Quiet) RESULT(err_stat)
    TYPE(PARMIOCoeff_type), INTENT(OUT) :: PARMIOCoeff
    CHARACTER(*),           INTENT(IN)  :: Filename
    LOGICAL, OPTIONAL,      INTENT(IN)  :: Quiet
    INTEGER :: err_stat
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'PARMIOCoeff_netCDF_ReadFile'
    CHARACTER(ML) :: msg
    LOGICAL :: Noisy, Close_File
    INTEGER :: NF90_Status, FileId, g, gid, rcid, rcgid

    err_stat = SUCCESS
    Close_File = .FALSE.
    Noisy = .TRUE.; IF (PRESENT(Quiet)) Noisy = .NOT. Quiet

    IF (.NOT. File_Exists(Filename)) THEN
      msg = 'File '//TRIM(Filename)//' not found.'
      CALL Read_Cleanup(); RETURN
    END IF

    NF90_Status = NF90_OPEN(Filename, NF90_NOWRITE, FileId)
    IF (NF90_Status /= NF90_NOERR) THEN
      msg = 'Error opening '//TRIM(Filename)//' for read - '// &
            TRIM(NF90_STRERROR(NF90_Status))
      CALL Read_Cleanup(); RETURN
    END IF
    Close_File = .TRUE.

    ! Allocate master record
    CALL PARMIOCoeff_Create(PARMIOCoeff)

    ! Read root global attributes (best-effort)
    CALL Get_Global_Att_String(FileId, 'grid_name',           PARMIOCoeff%Grid_Name)
    CALL Get_Global_Att_String(FileId, 'source_rows_csv',     PARMIOCoeff%Source_Rows_CSV)
    CALL Get_Global_Att_String(FileId, 'permittivity_policy', PARMIOCoeff%Permittivity_Policy)
    CALL Get_Global_Att_String(FileId, 'foam_policy',         PARMIOCoeff%Foam_Policy)
    CALL Get_Global_Att_String(FileId, 'coefficient_units',   PARMIOCoeff%Coefficient_Units)

    ! Load each of the three groups
    DO g = 1, PARMIO_N_GROUPS
      NF90_Status = NF90_INQ_NCID(FileId, TRIM(GROUP_NAME(g)), gid)
      IF (NF90_Status /= NF90_NOERR) THEN
        msg = 'Group "'//TRIM(GROUP_NAME(g))//'" not found in '// &
              TRIM(Filename)//' - '//TRIM(NF90_STRERROR(NF90_Status))
        CALL Read_Cleanup(); RETURN
      END IF
      CALL Read_Group(gid, PARMIOCoeff%Group(g), msg, err_stat)
      IF (err_stat /= SUCCESS) THEN
        msg = 'In group "'//TRIM(GROUP_NAME(g))//'": '//TRIM(msg)
        CALL Read_Cleanup(); RETURN
      END IF
    END DO

    ! Optional PARMIO-native reflection-correction groups. Legacy LUT files
    ! without these groups remain valid and continue to fall back to the
    ! FASTEM reflection correction in the runtime module.
    NF90_Status = NF90_INQ_NCID(FileId, 'reflection_correction', rcid)
    IF (NF90_Status == NF90_NOERR) THEN
      DO g = 1, PARMIO_N_GROUPS
        NF90_Status = NF90_INQ_NCID(rcid, TRIM(GROUP_NAME(g)), rcgid)
        IF (NF90_Status /= NF90_NOERR) CYCLE
        CALL Read_RC_Group(rcgid, PARMIOCoeff%RC_Group(g), msg, err_stat)
        IF (err_stat /= SUCCESS) THEN
          msg = 'In reflection_correction/'//TRIM(GROUP_NAME(g))//': '//TRIM(msg)
          CALL Read_Cleanup(); RETURN
        END IF
      END DO
    END IF

    NF90_Status = NF90_CLOSE(FileId); Close_File = .FALSE.
    IF (NF90_Status /= NF90_NOERR) THEN
      msg = 'Error closing '//TRIM(Filename)//' - '// &
            TRIM(NF90_STRERROR(NF90_Status))
      CALL Read_Cleanup(); RETURN
    END IF

    IF (Noisy) THEN
      CALL PARMIOCoeff_Info(PARMIOCoeff, msg)
      CALL Display_Message(ROUTINE_NAME, &
                           'FILE: '//TRIM(Filename)//'; '//TRIM(msg), &
                           INFORMATION)
    END IF

  CONTAINS

    SUBROUTINE Read_Cleanup()
      IF (Close_File) THEN
        NF90_Status = NF90_CLOSE(FileId)
      END IF
      CALL PARMIOCoeff_Destroy(PARMIOCoeff)
      err_stat = FAILURE
      CALL Display_Message(ROUTINE_NAME, msg, err_stat)
    END SUBROUTINE Read_Cleanup

  END FUNCTION PARMIOCoeff_netCDF_ReadFile


  !-----------------------------------------------------------------
  !  Read one group into a PARMIOCoeff_Group_type instance.
  !-----------------------------------------------------------------
  SUBROUTINE Read_Group(gid, gcoeff, msg, err_stat)
    INTEGER,                       INTENT(IN)  :: gid
    TYPE(PARMIOCoeff_Group_type),  INTENT(OUT) :: gcoeff
    CHARACTER(*),                  INTENT(OUT) :: msg
    INTEGER,                       INTENT(OUT) :: err_stat
    INTEGER :: NF90_Status
    INTEGER :: nF, nT, nU, nS, nQ, nFoam
    LOGICAL :: has_sss
    INTEGER :: VarId, dimid_sss
    CHARACTER(16) :: sss_attr
    INTEGER :: k
    REAL(Double), ALLOCATABLE :: buf6(:,:,:,:,:,:)
    REAL(Double), ALLOCATABLE :: buf5(:,:,:,:,:)

    err_stat = FAILURE
    msg = ''

    ! Detect SSS axis presence by trying to look up the dim.
    NF90_Status = NF90_INQ_DIMID(gid, DIM_SSS, dimid_sss)
    has_sss = (NF90_Status == NF90_NOERR)

    ! Resolve dimension lengths
    IF (.NOT. Get_Dim_Len(gid, DIM_FREQUENCY,  nF, msg)) RETURN
    IF (.NOT. Get_Dim_Len(gid, DIM_THETA,      nT, msg)) RETURN
    IF (.NOT. Get_Dim_Len(gid, DIM_WIND_SPEED, nU, msg)) RETURN
    IF (.NOT. Get_Dim_Len(gid, DIM_SST,        nS, msg)) RETURN
    IF (.NOT. Get_Dim_Len(gid, DIM_FOAM_STATE, nFoam, msg)) RETURN
    IF (has_sss) THEN
      IF (.NOT. Get_Dim_Len(gid, DIM_SSS, nQ, msg)) RETURN
    ELSE
      nQ = 1
    END IF

    ! Allocate the group's storage
    CALL Group_Create_Wrapper(gcoeff, nF, nT, nU, nS, nQ, nFoam, has_sss)
    IF (.NOT. gcoeff%Is_Allocated) THEN
      msg = 'Failed to allocate group storage'
      RETURN
    END IF

    ! Read axis vectors
    IF (.NOT. Get_Var_1D(gid, DIM_FREQUENCY,  gcoeff%Frequency,  msg)) RETURN
    IF (.NOT. Get_Var_1D(gid, DIM_THETA,      gcoeff%Theta,      msg)) RETURN
    IF (.NOT. Get_Var_1D(gid, DIM_WIND_SPEED, gcoeff%Wind_Speed, msg)) RETURN
    IF (.NOT. Get_Var_1D(gid, DIM_SST,        gcoeff%SST,        msg)) RETURN
    IF (has_sss) THEN
      IF (.NOT. Get_Var_1D(gid, DIM_SSS,      gcoeff%SSS,        msg)) RETURN
    ELSE
      gcoeff%SSS = 35.0_Double
    END IF

    ! Read SSS-axis-active group attribute (string "true" / "false")
    sss_attr = ''
    NF90_Status = NF90_GET_ATT(gid, NF90_GLOBAL, 'sss_axis_active', sss_attr)
    IF (NF90_Status == NF90_NOERR) THEN
      gcoeff%SSS_Axis_Active = (TRIM(sss_attr) == 'true')
    ELSE
      gcoeff%SSS_Axis_Active = has_sss
    END IF

    ! Read 14 harmonic coefficient variables. Each is shape
    ! (freq, theta, U10, sst, [sss], foam_state) on disk; Fortran sees
    ! it reversed: (foam_state, [sss], sst, U10, theta, freq).
    IF (has_sss) THEN
      ALLOCATE(buf6(nFoam, nQ, nS, nU, nT, nF))
    ELSE
      ALLOCATE(buf5(nFoam, nS, nU, nT, nF))
    END IF

    DO k = 1, N_PARMIO_HARMONIC_TERMS
      NF90_Status = NF90_INQ_VARID(gid, TRIM(HARMONIC_VARNAMES(k)), VarId)
      IF (NF90_Status /= NF90_NOERR) THEN
        msg = 'Variable "'//TRIM(HARMONIC_VARNAMES(k))//'" not found - '// &
              TRIM(NF90_STRERROR(NF90_Status))
        RETURN
      END IF
      IF (has_sss) THEN
        NF90_Status = NF90_GET_VAR(gid, VarId, buf6)
        IF (NF90_Status /= NF90_NOERR) THEN
          msg = 'Error reading "'//TRIM(HARMONIC_VARNAMES(k))//'" - '// &
                TRIM(NF90_STRERROR(NF90_Status))
          RETURN
        END IF
        gcoeff%Coefficients(k, :, :, :, :, :, :) = buf6
      ELSE
        NF90_Status = NF90_GET_VAR(gid, VarId, buf5)
        IF (NF90_Status /= NF90_NOERR) THEN
          msg = 'Error reading "'//TRIM(HARMONIC_VARNAMES(k))//'" - '// &
                TRIM(NF90_STRERROR(NF90_Status))
          RETURN
        END IF
        ! Inject an SSS dim of length 1 by reshape-style assignment.
        gcoeff%Coefficients(k, :, 1, :, :, :, :) = buf5
      END IF
    END DO

    ! Read Foam fraction
    NF90_Status = NF90_INQ_VARID(gid, VAR_FOAM, VarId)
    IF (NF90_Status /= NF90_NOERR) THEN
      msg = 'Variable "Foam" not found - '//TRIM(NF90_STRERROR(NF90_Status))
      RETURN
    END IF
    IF (has_sss) THEN
      NF90_Status = NF90_GET_VAR(gid, VarId, buf6)
      IF (NF90_Status /= NF90_NOERR) THEN
        msg = 'Error reading "Foam" - '//TRIM(NF90_STRERROR(NF90_Status))
        RETURN
      END IF
      gcoeff%Foam = buf6
    ELSE
      NF90_Status = NF90_GET_VAR(gid, VarId, buf5)
      IF (NF90_Status /= NF90_NOERR) THEN
        msg = 'Error reading "Foam" - '//TRIM(NF90_STRERROR(NF90_Status))
        RETURN
      END IF
      gcoeff%Foam(:, 1, :, :, :, :) = buf5
    END IF

    IF (ALLOCATED(buf6)) DEALLOCATE(buf6)
    IF (ALLOCATED(buf5)) DEALLOCATE(buf5)
    err_stat = SUCCESS
  END SUBROUTINE Read_Group


  !-----------------------------------------------------------------
  !  Read one optional reflection-correction group.
  !-----------------------------------------------------------------
  SUBROUTINE Read_RC_Group(gid, rcoeff, msg, err_stat)
    INTEGER,                          INTENT(IN)  :: gid
    TYPE(PARMIOCoeff_RC_Group_type),  INTENT(OUT) :: rcoeff
    CHARACTER(*),                     INTENT(OUT) :: msg
    INTEGER,                          INTENT(OUT) :: err_stat
    INTEGER :: NF90_Status
    INTEGER :: nF, nT, nU, nS, nQ, nFoam, nTau
    LOGICAL :: has_sss
    INTEGER :: dimid_sss
    REAL(Double), ALLOCATABLE :: buf7(:,:,:,:,:,:,:)
    REAL(Double), ALLOCATABLE :: buf6(:,:,:,:,:,:)

    err_stat = FAILURE
    msg = ''

    NF90_Status = NF90_INQ_DIMID(gid, DIM_SSS, dimid_sss)
    has_sss = (NF90_Status == NF90_NOERR)

    IF (.NOT. Get_Dim_Len(gid, DIM_FREQUENCY,     nF,    msg)) RETURN
    IF (.NOT. Get_Dim_Len(gid, DIM_THETA,         nT,    msg)) RETURN
    IF (.NOT. Get_Dim_Len(gid, DIM_WIND_SPEED,    nU,    msg)) RETURN
    IF (.NOT. Get_Dim_Len(gid, DIM_SST,           nS,    msg)) RETURN
    IF (.NOT. Get_Dim_Len(gid, DIM_FOAM_STATE,    nFoam, msg)) RETURN
    IF (.NOT. Get_Dim_Len(gid, DIM_TRANSMITTANCE, nTau,  msg)) RETURN
    IF (has_sss) THEN
      IF (.NOT. Get_Dim_Len(gid, DIM_SSS, nQ, msg)) RETURN
    ELSE
      nQ = 1
    END IF

    CALL RC_Group_Create_Wrapper(rcoeff, nF, nT, nU, nS, nQ, nFoam, nTau, has_sss)
    IF (.NOT. rcoeff%Is_Allocated) THEN
      msg = 'Failed to allocate reflection-correction group storage'
      RETURN
    END IF

    IF (.NOT. Get_Var_1D(gid, DIM_FREQUENCY,     rcoeff%Frequency,     msg)) RETURN
    IF (.NOT. Get_Var_1D(gid, DIM_THETA,         rcoeff%Theta,         msg)) RETURN
    IF (.NOT. Get_Var_1D(gid, DIM_WIND_SPEED,    rcoeff%Wind_Speed,    msg)) RETURN
    IF (.NOT. Get_Var_1D(gid, DIM_SST,           rcoeff%SST,           msg)) RETURN
    IF (.NOT. Get_Var_1D(gid, DIM_TRANSMITTANCE, rcoeff%Transmittance, msg)) RETURN
    IF (has_sss) THEN
      IF (.NOT. Get_Var_1D(gid, DIM_SSS,         rcoeff%SSS,           msg)) RETURN
    ELSE
      rcoeff%SSS = 35.0_Double
    END IF

    IF (has_sss) THEN
      ALLOCATE(buf7(nTau, nFoam, nQ, nS, nU, nT, nF))
      IF (.NOT. Get_Var_7D(gid, VAR_RDOWN_V, buf7, msg)) RETURN
      rcoeff%Rdown(PARMIO_RC_V_POL, :, :, :, :, :, :, :) = buf7
      IF (.NOT. Get_Var_7D(gid, VAR_RDOWN_H, buf7, msg)) RETURN
      rcoeff%Rdown(PARMIO_RC_H_POL, :, :, :, :, :, :, :) = buf7
    ELSE
      ALLOCATE(buf6(nTau, nFoam, nS, nU, nT, nF))
      IF (.NOT. Get_Var_6D(gid, VAR_RDOWN_V, buf6, msg)) RETURN
      rcoeff%Rdown(PARMIO_RC_V_POL, :, :, 1, :, :, :, :) = buf6
      IF (.NOT. Get_Var_6D(gid, VAR_RDOWN_H, buf6, msg)) RETURN
      rcoeff%Rdown(PARMIO_RC_H_POL, :, :, 1, :, :, :, :) = buf6
    END IF

    IF (ALLOCATED(buf7)) DEALLOCATE(buf7)
    IF (ALLOCATED(buf6)) DEALLOCATE(buf6)
    err_stat = SUCCESS
  END SUBROUTINE Read_RC_Group


  !-----------------------------------------------------------------
  !  Helper: dimension-length lookup
  !-----------------------------------------------------------------
  LOGICAL FUNCTION Get_Dim_Len(ncid, dimname, n, msg) RESULT(ok)
    INTEGER,      INTENT(IN)  :: ncid
    CHARACTER(*), INTENT(IN)  :: dimname
    INTEGER,      INTENT(OUT) :: n
    CHARACTER(*), INTENT(OUT) :: msg
    INTEGER :: dimid, status
    ok = .FALSE.
    status = NF90_INQ_DIMID(ncid, dimname, dimid)
    IF (status /= NF90_NOERR) THEN
      msg = 'Dimension "'//TRIM(dimname)//'" not found - '// &
            TRIM(NF90_STRERROR(status))
      RETURN
    END IF
    status = NF90_INQUIRE_DIMENSION(ncid, dimid, len=n)
    IF (status /= NF90_NOERR) THEN
      msg = 'Error inquiring dimension "'//TRIM(dimname)//'" - '// &
            TRIM(NF90_STRERROR(status))
      RETURN
    END IF
    ok = .TRUE.
  END FUNCTION Get_Dim_Len


  !-----------------------------------------------------------------
  !  Helper: 1-D real(double) variable read
  !-----------------------------------------------------------------
  LOGICAL FUNCTION Get_Var_1D(ncid, varname, arr, msg) RESULT(ok)
    INTEGER,      INTENT(IN)  :: ncid
    CHARACTER(*), INTENT(IN)  :: varname
    REAL(Double), INTENT(OUT) :: arr(:)
    CHARACTER(*), INTENT(OUT) :: msg
    INTEGER :: varid, status
    ok = .FALSE.
    status = NF90_INQ_VARID(ncid, varname, varid)
    IF (status /= NF90_NOERR) THEN
      msg = 'Variable "'//TRIM(varname)//'" not found - '// &
            TRIM(NF90_STRERROR(status))
      RETURN
    END IF
    status = NF90_GET_VAR(ncid, varid, arr)
    IF (status /= NF90_NOERR) THEN
      msg = 'Error reading variable "'//TRIM(varname)//'" - '// &
            TRIM(NF90_STRERROR(status))
      RETURN
    END IF
    ok = .TRUE.
  END FUNCTION Get_Var_1D


  LOGICAL FUNCTION Get_Var_6D(ncid, varname, arr, msg) RESULT(ok)
    INTEGER,      INTENT(IN)  :: ncid
    CHARACTER(*), INTENT(IN)  :: varname
    REAL(Double), INTENT(OUT) :: arr(:,:,:,:,:,:)
    CHARACTER(*), INTENT(OUT) :: msg
    INTEGER :: varid, status
    ok = .FALSE.
    status = NF90_INQ_VARID(ncid, varname, varid)
    IF (status /= NF90_NOERR) THEN
      msg = 'Variable "'//TRIM(varname)//'" not found - '//TRIM(NF90_STRERROR(status))
      RETURN
    END IF
    status = NF90_GET_VAR(ncid, varid, arr)
    IF (status /= NF90_NOERR) THEN
      msg = 'Error reading variable "'//TRIM(varname)//'" - '//TRIM(NF90_STRERROR(status))
      RETURN
    END IF
    ok = .TRUE.
  END FUNCTION Get_Var_6D


  LOGICAL FUNCTION Get_Var_7D(ncid, varname, arr, msg) RESULT(ok)
    INTEGER,      INTENT(IN)  :: ncid
    CHARACTER(*), INTENT(IN)  :: varname
    REAL(Double), INTENT(OUT) :: arr(:,:,:,:,:,:,:)
    CHARACTER(*), INTENT(OUT) :: msg
    INTEGER :: varid, status
    ok = .FALSE.
    status = NF90_INQ_VARID(ncid, varname, varid)
    IF (status /= NF90_NOERR) THEN
      msg = 'Variable "'//TRIM(varname)//'" not found - '//TRIM(NF90_STRERROR(status))
      RETURN
    END IF
    status = NF90_GET_VAR(ncid, varid, arr)
    IF (status /= NF90_NOERR) THEN
      msg = 'Error reading variable "'//TRIM(varname)//'" - '//TRIM(NF90_STRERROR(status))
      RETURN
    END IF
    ok = .TRUE.
  END FUNCTION Get_Var_7D


  !-----------------------------------------------------------------
  !  Helper: best-effort global-attribute string read
  !-----------------------------------------------------------------
  SUBROUTINE Get_Global_Att_String(ncid, attname, value)
    INTEGER,      INTENT(IN)    :: ncid
    CHARACTER(*), INTENT(IN)    :: attname
    CHARACTER(*), INTENT(INOUT) :: value
    INTEGER :: status
    status = NF90_GET_ATT(ncid, NF90_GLOBAL, attname, value)
    IF (status /= NF90_NOERR) value = ''
  END SUBROUTINE Get_Global_Att_String


  !-----------------------------------------------------------------
  !  Wrapper for the (private) Group_Create routine in the Define
  !  module. Pulled out so Read_Group can call it without exporting
  !  the create routine to the outside world.
  !-----------------------------------------------------------------
  SUBROUTINE Group_Create_Wrapper(self, nF, nT, nU, nS, nQ, nFoam, has_sss)
    TYPE(PARMIOCoeff_Group_type), INTENT(OUT) :: self
    INTEGER, INTENT(IN) :: nF, nT, nU, nS, nQ, nFoam
    LOGICAL, INTENT(IN) :: has_sss
    ! Inline allocation matching PARMIOCoeff_Group_Create in the Define
    ! module. The Define-module routine is PURE; we can also call it
    ! directly via the alias above, but doing so requires that routine
    ! to be PUBLIC. Replicating the allocation here keeps the Define
    ! module's API surface minimal.
    INTEGER :: alloc_stat
    IF (nF < 1 .OR. nT < 1 .OR. nU < 1 .OR. nS < 1 .OR. nQ < 1 .OR. nFoam < 1) RETURN
    ALLOCATE( &
      self%Frequency(nF),                                  &
      self%Theta(nT),                                      &
      self%Wind_Speed(nU),                                 &
      self%SST(nS),                                        &
      self%SSS(nQ),                                        &
      self%Confidence_Label(nF),                           &
      self%Coefficients(N_PARMIO_HARMONIC_TERMS, nFoam,    &
                        nQ, nS, nU, nT, nF),               &
      self%Foam(nFoam, nQ, nS, nU, nT, nF),                &
      STAT=alloc_stat)
    IF (alloc_stat /= 0) RETURN
    self%n_Frequencies   = nF
    self%n_Angles        = nT
    self%n_Wind_Speeds   = nU
    self%n_SSTs          = nS
    self%n_SSSs          = nQ
    self%n_Foam_States   = nFoam
    self%SSS_Axis_Active = has_sss
    self%Frequency        = 0.0_Double
    self%Theta            = 0.0_Double
    self%Wind_Speed       = 0.0_Double
    self%SST              = 0.0_Double
    self%SSS              = 0.0_Double
    self%Confidence_Label = ''
    self%Coefficients     = 0.0_Double
    self%Foam             = 0.0_Double
    self%Is_Allocated     = .TRUE.
  END SUBROUTINE Group_Create_Wrapper


  SUBROUTINE RC_Group_Create_Wrapper(self, nF, nT, nU, nS, nQ, nFoam, nTau, has_sss)
    TYPE(PARMIOCoeff_RC_Group_type), INTENT(OUT) :: self
    INTEGER, INTENT(IN) :: nF, nT, nU, nS, nQ, nFoam, nTau
    LOGICAL, INTENT(IN) :: has_sss
    INTEGER :: alloc_stat
    IF (nF < 1 .OR. nT < 1 .OR. nU < 1 .OR. nS < 1 .OR. &
        nQ < 1 .OR. nFoam < 1 .OR. nTau < 1) RETURN
    ALLOCATE( &
      self%Frequency(nF),                                      &
      self%Theta(nT),                                          &
      self%Wind_Speed(nU),                                     &
      self%SST(nS),                                            &
      self%SSS(nQ),                                            &
      self%Transmittance(nTau),                                &
      self%Rdown(PARMIO_N_RC_POLARIZATIONS, nTau, nFoam,       &
                 nQ, nS, nU, nT, nF),                          &
      STAT=alloc_stat)
    IF (alloc_stat /= 0) RETURN
    self%n_Frequencies    = nF
    self%n_Angles         = nT
    self%n_Wind_Speeds    = nU
    self%n_SSTs           = nS
    self%n_SSSs           = nQ
    self%n_Foam_States    = nFoam
    self%n_Transmittances = nTau
    self%SSS_Axis_Active  = has_sss
    self%Frequency        = 0.0_Double
    self%Theta            = 0.0_Double
    self%Wind_Speed       = 0.0_Double
    self%SST              = 0.0_Double
    self%SSS              = 0.0_Double
    self%Transmittance    = 0.0_Double
    self%Rdown            = 0.0_Double
    self%Is_Allocated     = .TRUE.
  END SUBROUTINE RC_Group_Create_Wrapper

END MODULE PARMIOCoeff_netCDF_IO
