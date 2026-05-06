!
! PARMIOCoeff_Define
!
! Module defining the PARMIOCoeff object that holds the lookup-table
! coefficients used by the PARMIO_MWSSEM ocean surface emissivity model.
!
! The LUT is built offline from the PARMIO reference radiative transfer
! model (Dinnat 2023). It stores PARMIO's full azimuthal-harmonic
! decomposition of brightness temperature, divided by SST_K, on a 5-D
! grid (frequency, zenith angle, 10-m wind speed, SST, SSS). To respect
! PARMIO's discontinuous physics three coefficient *groups* are kept:
!
!   sss_dependent    : f <= 10.65 GHz, Meissner permittivity, SSS axis
!   sss_nominal_m    : 10.65 < f < 28.8 GHz, Meissner, SSS = 35 only
!   sss_nominal_h    : f >= 28.8 GHz, high-freq tabulated dielectric
!
! At runtime the PARMIO_MWSSEM module looks up the 14 harmonic terms in
! the appropriate group, recombines them through cos/sin in azimuth, and
! returns the 4-Stokes (V, H, U, V) emissivity / reflectivity at FastemX
! call sites.

MODULE PARMIOCoeff_Define

  USE Type_Kinds           , ONLY: fp, Long, Double
  USE Message_Handler      , ONLY: SUCCESS, FAILURE, INFORMATION, Display_Message
  USE Compare_Float_Numbers, ONLY: OPERATOR(.EqualTo.)

  IMPLICIT NONE

  PRIVATE

  ! Datatypes
  PUBLIC :: PARMIOCoeff_Group_type
  PUBLIC :: PARMIOCoeff_RC_Group_type
  PUBLIC :: PARMIOCoeff_type
  ! Operators
  PUBLIC :: OPERATOR(==)
  ! Procedures
  PUBLIC :: PARMIOCoeff_Associated
  PUBLIC :: PARMIOCoeff_Destroy
  PUBLIC :: PARMIOCoeff_Create
  PUBLIC :: PARMIOCoeff_Inspect
  PUBLIC :: PARMIOCoeff_ValidRelease
  PUBLIC :: PARMIOCoeff_Info
  ! Group-name helpers
  PUBLIC :: PARMIOCoeff_GroupName_For_Frequency

  ! ---------------------
  ! Procedure overloading
  ! ---------------------
  INTERFACE OPERATOR(==)
    MODULE PROCEDURE PARMIOCoeff_Equal
  END INTERFACE OPERATOR(==)

  ! -----------------
  ! Module parameters
  ! -----------------
  INTEGER, PARAMETER :: PARMIOCOEFF_RELEASE = 1
  INTEGER, PARAMETER :: PARMIOCOEFF_VERSION = 1
  INTEGER, PARAMETER :: ML = 256
  INTEGER, PARAMETER :: SL =  80
  REAL(fp), PARAMETER :: ZERO = 0.0_fp

  ! Frequency thresholds (GHz) for group selection. Must match the values
  ! used by the offline LUT generator (parmio/scripts/parmio_lut_grid.py).
  REAL(fp), PARAMETER :: SSS_CUTOFF_GHZ          = 10.65_fp
  REAL(fp), PARAMETER :: PERMITTIVITY_SWITCH_GHZ = 28.8_fp

  ! Number of harmonic coefficient terms stored per (axis, foam_state) cell.
  ! Order matches Tb.f output columns 7-18, 26-27 divided by SST_K:
  !   1: TvN     (V-pol specular emissivity)
  !   2: ThN     (H-pol specular emissivity)
  !   3: Tv0     (V-pol roughness 0th harmonic)
  !   4: Th0     (H-pol roughness 0th harmonic)
  !   5: Tv1     (V-pol 1st harmonic, cos phi)
  !   6: Th1     (H-pol 1st harmonic, cos phi)
  !   7: U1      (3rd Stokes 1st harmonic, sin phi)
  !   8: V1      (4th Stokes 1st harmonic, sin phi)
  !   9: Tv2     (V-pol 2nd harmonic, cos 2phi)
  !  10: Th2     (H-pol 2nd harmonic, cos 2phi)
  !  11: U2      (3rd Stokes 2nd harmonic, sin 2phi)
  !  12: V2      (4th Stokes 2nd harmonic, sin 2phi)
  !  13: dTV_MR  (V-pol multi-reflection correction)
  !  14: dTh_MR  (H-pol multi-reflection correction)
  INTEGER, PARAMETER, PUBLIC :: N_PARMIO_HARMONIC_TERMS = 14

  ! Foam state index: 1 = foam-off, 2 = foam-on. Mirrors the offline writer's
  ! foam_state(0=nofoam, 1=foam) zero-based ordering.
  INTEGER, PARAMETER, PUBLIC :: PARMIO_FOAM_OFF = 1
  INTEGER, PARAMETER, PUBLIC :: PARMIO_FOAM_ON  = 2

  ! Group identifiers
  INTEGER, PARAMETER, PUBLIC :: PARMIO_GROUP_SSS_DEPENDENT  = 1
  INTEGER, PARAMETER, PUBLIC :: PARMIO_GROUP_SSS_NOMINAL_M  = 2
  INTEGER, PARAMETER, PUBLIC :: PARMIO_GROUP_SSS_NOMINAL_H  = 3
  INTEGER, PARAMETER, PUBLIC :: PARMIO_N_GROUPS             = 3
  INTEGER, PARAMETER, PUBLIC :: PARMIO_RC_V_POL             = 1
  INTEGER, PARAMETER, PUBLIC :: PARMIO_RC_H_POL             = 2
  INTEGER, PARAMETER, PUBLIC :: PARMIO_N_RC_POLARIZATIONS   = 2

  ! ----------------------------------
  ! PARMIO LUT group data type
  ! ----------------------------------
  ! One PARMIOCoeff_Group_type per coefficient group. The SSS axis is
  ! optional: groups without an active SSS axis have n_SSS = 1 and a
  ! single-element SSS coordinate at the nominal value.
  TYPE :: PARMIOCoeff_Group_type
    LOGICAL       :: Is_Allocated   = .FALSE.
    LOGICAL       :: SSS_Axis_Active = .FALSE.
    ! Dimensions
    INTEGER(Long) :: n_Frequencies  = 0
    INTEGER(Long) :: n_Angles       = 0
    INTEGER(Long) :: n_Wind_Speeds  = 0
    INTEGER(Long) :: n_SSTs         = 0
    INTEGER(Long) :: n_SSSs         = 0
    INTEGER(Long) :: n_Foam_States  = 0  ! always 2 in current LUTs
    ! Axis vectors
    REAL(Double), ALLOCATABLE :: Frequency(:)         ! GHz, n_Frequencies
    REAL(Double), ALLOCATABLE :: Theta(:)             ! deg, n_Angles
    REAL(Double), ALLOCATABLE :: Wind_Speed(:)        ! m/s, n_Wind_Speeds
    REAL(Double), ALLOCATABLE :: SST(:)               ! deg C, n_SSTs
    REAL(Double), ALLOCATABLE :: SSS(:)               ! psu, n_SSSs (=1 if SSS not active)
    ! Per-frequency confidence label (validated/extrapolated-defensible/...)
    CHARACTER(SL), ALLOCATABLE :: Confidence_Label(:) ! n_Frequencies
    ! Harmonic coefficient table — dimensionless (PARMIO Tb / SST_K).
    ! Dimension order is the natural Fortran layout for direct reads from
    ! the netCDF file (which writes with C-order dims
    ! frequency,theta,wind_speed,sst,sss,foam_state):
    !   Coefficients(harmonic_idx, foam_state, sss, sst, wind_speed, theta, freq)
    ! Harmonic index varies fastest so that the 14 terms at a given
    ! (foam_state, sss, sst, U10, theta, freq) lattice point are contiguous.
    REAL(Double), ALLOCATABLE :: Coefficients(:,:,:,:,:,:,:)
    ! Foam fraction (percent), same axis order (no harmonic dim):
    !   Foam(foam_state, sss, sst, wind_speed, theta, freq)
    REAL(Double), ALLOCATABLE :: Foam(:,:,:,:,:,:)
  END TYPE PARMIOCoeff_Group_type

  ! ----------------------------------
  ! PARMIO reflection-correction group
  ! ----------------------------------
  ! Optional group that carries PARMIO-native effective reflectivity of
  ! downwelling atmospheric radiation. It is intentionally separate from the
  ! emissivity harmonic coefficients so legacy coefficient files without this
  ! group remain readable.
  TYPE :: PARMIOCoeff_RC_Group_type
    LOGICAL       :: Is_Allocated    = .FALSE.
    LOGICAL       :: SSS_Axis_Active = .FALSE.
    INTEGER(Long) :: n_Frequencies   = 0
    INTEGER(Long) :: n_Angles        = 0
    INTEGER(Long) :: n_Wind_Speeds   = 0
    INTEGER(Long) :: n_SSTs          = 0
    INTEGER(Long) :: n_SSSs          = 0
    INTEGER(Long) :: n_Foam_States   = 0
    INTEGER(Long) :: n_Transmittances = 0
    REAL(Double), ALLOCATABLE :: Frequency(:)      ! GHz
    REAL(Double), ALLOCATABLE :: Theta(:)          ! deg
    REAL(Double), ALLOCATABLE :: Wind_Speed(:)     ! m/s
    REAL(Double), ALLOCATABLE :: SST(:)            ! deg C
    REAL(Double), ALLOCATABLE :: SSS(:)            ! psu
    REAL(Double), ALLOCATABLE :: Transmittance(:)  ! 1
    ! Rdown(pol, transmittance, foam_state, sss, sst, U10, theta, freq)
    REAL(Double), ALLOCATABLE :: Rdown(:,:,:,:,:,:,:,:)
  END TYPE PARMIOCoeff_RC_Group_type

  ! ----------------------------------
  ! PARMIOCoeff master data type
  ! ----------------------------------
  TYPE :: PARMIOCoeff_type
    LOGICAL       :: Is_Allocated = .FALSE.
    INTEGER(Long) :: Release      = PARMIOCOEFF_RELEASE
    INTEGER(Long) :: Version      = PARMIOCOEFF_VERSION
    ! Threshold metadata (loaded from netCDF global attrs; copied here
    ! so client code does not need the strings).
    REAL(Double)  :: SSS_Cutoff_GHz          = SSS_CUTOFF_GHZ
    REAL(Double)  :: Permittivity_Switch_GHz = PERMITTIVITY_SWITCH_GHZ
    ! Provenance / metadata
    CHARACTER(SL) :: Grid_Name        = ''
    CHARACTER(ML) :: Source_Rows_CSV  = ''
    CHARACTER(ML) :: Permittivity_Policy = ''
    CHARACTER(ML) :: Foam_Policy         = ''
    CHARACTER(ML) :: Coefficient_Units   = ''
    ! The three coefficient groups
    TYPE(PARMIOCoeff_Group_type) :: Group(PARMIO_N_GROUPS)
    TYPE(PARMIOCoeff_RC_Group_type) :: RC_Group(PARMIO_N_GROUPS)
  END TYPE PARMIOCoeff_type


CONTAINS


  !-----------------------------------------------------------------
  !  Group-name helper. Frequency-driven choice of PARMIO_GROUP_*.
  !-----------------------------------------------------------------
  PURE FUNCTION PARMIOCoeff_GroupName_For_Frequency( &
      Frequency_GHz, &
      SSS_Cutoff_GHz_Override, &
      Permittivity_Switch_GHz_Override) RESULT(group_id)
    REAL(fp),           INTENT(IN) :: Frequency_GHz
    REAL(fp), OPTIONAL, INTENT(IN) :: SSS_Cutoff_GHz_Override
    REAL(fp), OPTIONAL, INTENT(IN) :: Permittivity_Switch_GHz_Override
    INTEGER :: group_id
    REAL(fp) :: f_sss, f_perm
    f_sss  = SSS_CUTOFF_GHZ
    f_perm = PERMITTIVITY_SWITCH_GHZ
    IF (PRESENT(SSS_Cutoff_GHz_Override))          f_sss  = SSS_Cutoff_GHz_Override
    IF (PRESENT(Permittivity_Switch_GHz_Override)) f_perm = Permittivity_Switch_GHz_Override
    IF (Frequency_GHz <= f_sss) THEN
      group_id = PARMIO_GROUP_SSS_DEPENDENT
    ELSE IF (Frequency_GHz < f_perm) THEN
      group_id = PARMIO_GROUP_SSS_NOMINAL_M
    ELSE
      group_id = PARMIO_GROUP_SSS_NOMINAL_H
    END IF
  END FUNCTION PARMIOCoeff_GroupName_For_Frequency


  !-----------------------------------------------------------------
  !  Group-level allocation / deallocation
  !-----------------------------------------------------------------
  PURE FUNCTION PARMIOCoeff_Group_Associated(self) RESULT(Status)
    TYPE(PARMIOCoeff_Group_type), INTENT(IN) :: self
    LOGICAL :: Status
    Status = self%Is_Allocated
  END FUNCTION PARMIOCoeff_Group_Associated


  PURE SUBROUTINE PARMIOCoeff_Group_Destroy(self)
    TYPE(PARMIOCoeff_Group_type), INTENT(IN OUT) :: self
    self%Is_Allocated   = .FALSE.
    self%SSS_Axis_Active = .FALSE.
    self%n_Frequencies  = 0
    self%n_Angles       = 0
    self%n_Wind_Speeds  = 0
    self%n_SSTs         = 0
    self%n_SSSs         = 0
    self%n_Foam_States  = 0
    IF (ALLOCATED(self%Frequency))        DEALLOCATE(self%Frequency)
    IF (ALLOCATED(self%Theta))            DEALLOCATE(self%Theta)
    IF (ALLOCATED(self%Wind_Speed))       DEALLOCATE(self%Wind_Speed)
    IF (ALLOCATED(self%SST))              DEALLOCATE(self%SST)
    IF (ALLOCATED(self%SSS))              DEALLOCATE(self%SSS)
    IF (ALLOCATED(self%Confidence_Label)) DEALLOCATE(self%Confidence_Label)
    IF (ALLOCATED(self%Coefficients))     DEALLOCATE(self%Coefficients)
    IF (ALLOCATED(self%Foam))             DEALLOCATE(self%Foam)
  END SUBROUTINE PARMIOCoeff_Group_Destroy


  PURE SUBROUTINE PARMIOCoeff_RC_Group_Destroy(self)
    TYPE(PARMIOCoeff_RC_Group_type), INTENT(IN OUT) :: self
    self%Is_Allocated     = .FALSE.
    self%SSS_Axis_Active  = .FALSE.
    self%n_Frequencies    = 0
    self%n_Angles         = 0
    self%n_Wind_Speeds    = 0
    self%n_SSTs           = 0
    self%n_SSSs           = 0
    self%n_Foam_States    = 0
    self%n_Transmittances = 0
    IF (ALLOCATED(self%Frequency))     DEALLOCATE(self%Frequency)
    IF (ALLOCATED(self%Theta))         DEALLOCATE(self%Theta)
    IF (ALLOCATED(self%Wind_Speed))    DEALLOCATE(self%Wind_Speed)
    IF (ALLOCATED(self%SST))           DEALLOCATE(self%SST)
    IF (ALLOCATED(self%SSS))           DEALLOCATE(self%SSS)
    IF (ALLOCATED(self%Transmittance)) DEALLOCATE(self%Transmittance)
    IF (ALLOCATED(self%Rdown))         DEALLOCATE(self%Rdown)
  END SUBROUTINE PARMIOCoeff_RC_Group_Destroy


  PURE SUBROUTINE PARMIOCoeff_Group_Create( &
      self, n_Frequencies, n_Angles, n_Wind_Speeds, n_SSTs, n_SSSs, &
      n_Foam_States, SSS_Axis_Active)
    TYPE(PARMIOCoeff_Group_type), INTENT(OUT) :: self
    INTEGER, INTENT(IN) :: n_Frequencies, n_Angles, n_Wind_Speeds
    INTEGER, INTENT(IN) :: n_SSTs, n_SSSs, n_Foam_States
    LOGICAL, INTENT(IN) :: SSS_Axis_Active
    INTEGER :: alloc_stat
    IF (n_Frequencies < 1 .OR. n_Angles < 1 .OR. n_Wind_Speeds < 1 .OR. &
        n_SSTs < 1 .OR. n_SSSs < 1 .OR. n_Foam_States < 1) RETURN
    ALLOCATE( &
      self%Frequency       (n_Frequencies),        &
      self%Theta           (n_Angles),             &
      self%Wind_Speed      (n_Wind_Speeds),        &
      self%SST             (n_SSTs),               &
      self%SSS             (n_SSSs),               &
      self%Confidence_Label(n_Frequencies),        &
      ! Natural-Fortran layout matching netCDF C-order on disk
      self%Coefficients(N_PARMIO_HARMONIC_TERMS,   &
                        n_Foam_States, n_SSSs,     &
                        n_SSTs, n_Wind_Speeds,     &
                        n_Angles, n_Frequencies),  &
      self%Foam(n_Foam_States, n_SSSs, n_SSTs,     &
                n_Wind_Speeds, n_Angles,           &
                n_Frequencies),                    &
      STAT=alloc_stat)
    IF (alloc_stat /= 0) RETURN
    self%n_Frequencies   = n_Frequencies
    self%n_Angles        = n_Angles
    self%n_Wind_Speeds   = n_Wind_Speeds
    self%n_SSTs          = n_SSTs
    self%n_SSSs          = n_SSSs
    self%n_Foam_States   = n_Foam_States
    self%SSS_Axis_Active = SSS_Axis_Active
    self%Frequency        = ZERO
    self%Theta            = ZERO
    self%Wind_Speed       = ZERO
    self%SST              = ZERO
    self%SSS              = ZERO
    self%Confidence_Label = ''
    self%Coefficients     = ZERO
    self%Foam             = ZERO
    self%Is_Allocated     = .TRUE.
  END SUBROUTINE PARMIOCoeff_Group_Create


  !-----------------------------------------------------------------
  !  Master-level passthroughs over the three groups
  !-----------------------------------------------------------------
  PURE FUNCTION PARMIOCoeff_Associated(self) RESULT(Status)
    TYPE(PARMIOCoeff_type), INTENT(IN) :: self
    LOGICAL :: Status
    INTEGER :: g
    Status = self%Is_Allocated
    IF (.NOT. Status) RETURN
    DO g = 1, PARMIO_N_GROUPS
      IF (.NOT. PARMIOCoeff_Group_Associated(self%Group(g))) THEN
        Status = .FALSE.
        RETURN
      END IF
    END DO
  END FUNCTION PARMIOCoeff_Associated


  PURE SUBROUTINE PARMIOCoeff_Destroy(self)
    TYPE(PARMIOCoeff_type), INTENT(IN OUT) :: self
    INTEGER :: g
    DO g = 1, PARMIO_N_GROUPS
      CALL PARMIOCoeff_Group_Destroy(self%Group(g))
      CALL PARMIOCoeff_RC_Group_Destroy(self%RC_Group(g))
    END DO
    self%Is_Allocated        = .FALSE.
    self%Grid_Name           = ''
    self%Source_Rows_CSV     = ''
    self%Permittivity_Policy = ''
    self%Foam_Policy         = ''
    self%Coefficient_Units   = ''
  END SUBROUTINE PARMIOCoeff_Destroy


  SUBROUTINE PARMIOCoeff_Create(self)
    TYPE(PARMIOCoeff_type), INTENT(IN OUT) :: self
    self%Is_Allocated = .TRUE.
  END SUBROUTINE PARMIOCoeff_Create


  !-----------------------------------------------------------------
  !  Inspect / Info
  !-----------------------------------------------------------------
  SUBROUTINE PARMIOCoeff_Inspect(self)
    TYPE(PARMIOCoeff_type), INTENT(IN) :: self
    INTEGER :: g
    CHARACTER(*), PARAMETER :: GROUP_NAME(PARMIO_N_GROUPS) = (/ &
      'sss_dependent  ', 'sss_nominal_m  ', 'sss_nominal_h  ' /)
    WRITE(*,'(/," PARMIOCoeff RELEASE.VERSION: ",i0,".",i0)') self%Release, self%Version
    WRITE(*,'(  "   Grid_Name              : ",a)')  TRIM(self%Grid_Name)
    WRITE(*,'(  "   SSS cutoff (GHz)       : ",f8.3)') self%SSS_Cutoff_GHz
    WRITE(*,'(  "   Permittivity switch GHz: ",f8.3)') self%Permittivity_Switch_GHz
    WRITE(*,'(  "   Permittivity policy    : ",a)') TRIM(self%Permittivity_Policy)
    WRITE(*,'(  "   Foam policy            : ",a)') TRIM(self%Foam_Policy)
    DO g = 1, PARMIO_N_GROUPS
      WRITE(*,'(/,"  -- Group ",i1," : ",a," --")') g, TRIM(GROUP_NAME(g))
      IF (.NOT. PARMIOCoeff_Group_Associated(self%Group(g))) THEN
        WRITE(*,'("   (not allocated)")')
        CYCLE
      END IF
      WRITE(*,'("   SSS axis active: ",l1)') self%Group(g)%SSS_Axis_Active
      WRITE(*,'("   Dimensions     : freq=",i0," theta=",i0," U10=",i0, &
              & " SST=",i0," SSS=",i0," foam=",i0)')                    &
            self%Group(g)%n_Frequencies, self%Group(g)%n_Angles,        &
            self%Group(g)%n_Wind_Speeds, self%Group(g)%n_SSTs,          &
            self%Group(g)%n_SSSs,        self%Group(g)%n_Foam_States
      IF (self%Group(g)%n_Frequencies > 0) THEN
        WRITE(*,'("   Frequency range: ",f8.3,"  -  ",f8.3," GHz")')    &
              self%Group(g)%Frequency(1),                               &
              self%Group(g)%Frequency(self%Group(g)%n_Frequencies)
      END IF
      IF (self%RC_Group(g)%Is_Allocated) THEN
        WRITE(*,'("   RC dimensions  : freq=",i0," theta=",i0," U10=",i0, &
                & " SST=",i0," SSS=",i0," foam=",i0," trans=",i0)')      &
              self%RC_Group(g)%n_Frequencies, self%RC_Group(g)%n_Angles, &
              self%RC_Group(g)%n_Wind_Speeds, self%RC_Group(g)%n_SSTs,   &
              self%RC_Group(g)%n_SSSs, self%RC_Group(g)%n_Foam_States,   &
              self%RC_Group(g)%n_Transmittances
      END IF
    END DO
  END SUBROUTINE PARMIOCoeff_Inspect


  SUBROUTINE PARMIOCoeff_Info(self, Info)
    TYPE(PARMIOCoeff_type), INTENT(IN)  :: self
    CHARACTER(*),           INTENT(OUT) :: Info
    CHARACTER(8) :: rel, ver
    WRITE(rel, '(i0)') self%Release
    WRITE(ver, '(i0)') self%Version
    Info = 'PARMIOCoeff RELEASE.VERSION: '// &
           TRIM(rel)//'.'//TRIM(ver)//'; grid='//TRIM(self%Grid_Name)
  END SUBROUTINE PARMIOCoeff_Info


  PURE FUNCTION PARMIOCoeff_ValidRelease(self) RESULT(Is_Valid)
    TYPE(PARMIOCoeff_type), INTENT(IN) :: self
    LOGICAL :: Is_Valid
    Is_Valid = (self%Release == PARMIOCOEFF_RELEASE)
  END FUNCTION PARMIOCoeff_ValidRelease


  !-----------------------------------------------------------------
  !  Equality (axis-vector + harmonic-coefficient comparison)
  !-----------------------------------------------------------------
  FUNCTION PARMIOCoeff_Equal(x, y) RESULT(is_equal)
    TYPE(PARMIOCoeff_type), INTENT(IN) :: x, y
    LOGICAL :: is_equal
    INTEGER :: g
    is_equal = .FALSE.
    IF (x%Is_Allocated .NEQV. y%Is_Allocated) RETURN
    IF (x%Release /= y%Release) RETURN
    DO g = 1, PARMIO_N_GROUPS
      IF (.NOT. Group_Equal(x%Group(g), y%Group(g))) RETURN
      IF (.NOT. RC_Group_Equal(x%RC_Group(g), y%RC_Group(g))) RETURN
    END DO
    is_equal = .TRUE.
  END FUNCTION PARMIOCoeff_Equal


  PURE FUNCTION Group_Equal(x, y) RESULT(is_equal)
    TYPE(PARMIOCoeff_Group_type), INTENT(IN) :: x, y
    LOGICAL :: is_equal
    is_equal = .FALSE.
    IF (x%Is_Allocated .NEQV. y%Is_Allocated) RETURN
    IF (.NOT. x%Is_Allocated) THEN
      is_equal = .TRUE.
      RETURN
    END IF
    IF (x%n_Frequencies /= y%n_Frequencies) RETURN
    IF (x%n_Angles      /= y%n_Angles)      RETURN
    IF (x%n_Wind_Speeds /= y%n_Wind_Speeds) RETURN
    IF (x%n_SSTs        /= y%n_SSTs)        RETURN
    IF (x%n_SSSs        /= y%n_SSSs)        RETURN
    IF (x%n_Foam_States /= y%n_Foam_States) RETURN
    IF (.NOT. ALL(x%Frequency  .EqualTo. y%Frequency))  RETURN
    IF (.NOT. ALL(x%Theta      .EqualTo. y%Theta))      RETURN
    IF (.NOT. ALL(x%Wind_Speed .EqualTo. y%Wind_Speed)) RETURN
    IF (.NOT. ALL(x%SST        .EqualTo. y%SST))        RETURN
    IF (.NOT. ALL(x%SSS        .EqualTo. y%SSS))        RETURN
    IF (.NOT. ALL(x%Coefficients .EqualTo. y%Coefficients)) RETURN
    IF (.NOT. ALL(x%Foam         .EqualTo. y%Foam))         RETURN
    is_equal = .TRUE.
  END FUNCTION Group_Equal


  PURE FUNCTION RC_Group_Equal(x, y) RESULT(is_equal)
    TYPE(PARMIOCoeff_RC_Group_type), INTENT(IN) :: x, y
    LOGICAL :: is_equal
    is_equal = .FALSE.
    IF (x%Is_Allocated .NEQV. y%Is_Allocated) RETURN
    IF (.NOT. x%Is_Allocated) THEN
      is_equal = .TRUE.
      RETURN
    END IF
    IF (x%n_Frequencies    /= y%n_Frequencies)    RETURN
    IF (x%n_Angles         /= y%n_Angles)         RETURN
    IF (x%n_Wind_Speeds    /= y%n_Wind_Speeds)    RETURN
    IF (x%n_SSTs           /= y%n_SSTs)           RETURN
    IF (x%n_SSSs           /= y%n_SSSs)           RETURN
    IF (x%n_Foam_States    /= y%n_Foam_States)    RETURN
    IF (x%n_Transmittances /= y%n_Transmittances) RETURN
    IF (.NOT. ALL(x%Frequency     .EqualTo. y%Frequency))     RETURN
    IF (.NOT. ALL(x%Theta         .EqualTo. y%Theta))         RETURN
    IF (.NOT. ALL(x%Wind_Speed    .EqualTo. y%Wind_Speed))    RETURN
    IF (.NOT. ALL(x%SST           .EqualTo. y%SST))           RETURN
    IF (.NOT. ALL(x%SSS           .EqualTo. y%SSS))           RETURN
    IF (.NOT. ALL(x%Transmittance .EqualTo. y%Transmittance)) RETURN
    IF (.NOT. ALL(x%Rdown         .EqualTo. y%Rdown))         RETURN
    is_equal = .TRUE.
  END FUNCTION RC_Group_Equal

END MODULE PARMIOCoeff_Define
