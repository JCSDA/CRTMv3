!
! CloudCoeff_Exp_Define
!
! Module defining the EXPERIMENTAL CloudCoeff data structure (v1) and routines to
! manipulate it. This is the opt-in 'CRTM-Exp' cloud optics scheme; the legacy
! CloudCoeff_Define / CloudCoeff_type is unchanged and remains the default.
!
! Format (see CloudCoeff_Experimental_Schema_v1.md):
!   - explicit habit axis + (Dm, mu) PSD-moment axes + temperature (all phases)
!   - one full / variable-length GSF expansion per entry with a per-entry effective
!     truncation order (n_Legendre_Eff) -- DECOUPLED from the RT stream count
!   - 6 phase elements (alpha1..alpha4, beta1, beta2)
!   - bulk optics ke, ka, kb (+ g); w = (ke-ka)/ke derived at runtime
!
! CREATION HISTORY:
!       Written for the experimental cloud-optics redesign, 2026-06-01
!

MODULE CloudCoeff_Exp_Define

  ! ------------------
  ! Environment set up
  ! ------------------
  USE Type_Kinds,       ONLY: Long, Double
  USE Message_Handler,  ONLY: SUCCESS, FAILURE, INFORMATION, Display_Message
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  PUBLIC :: CloudCoeff_Exp_type
  PUBLIC :: CloudCoeff_Exp_Associated
  PUBLIC :: CloudCoeff_Exp_Destroy
  PUBLIC :: CloudCoeff_Exp_Create
  PUBLIC :: CloudCoeff_Exp_Inspect
  PUBLIC :: CloudCoeff_Exp_Info
  PUBLIC :: CloudCoeff_Exp_DefineVersion
  PUBLIC :: INVALID_CLOUDCOEFF_EXP
  PUBLIC :: CRTM_EXP_CLOUDCOEFF


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: CloudCoeff_Exp_Define.f90  v1  2026-06-01 $'
  REAL(Double), PARAMETER :: ZERO = 0.0_Double
  INTEGER, PARAMETER :: ML = 256
  INTEGER, PARAMETER :: SL = 32                       ! habit-name string length
  ! Release/scheme identifiers
  INTEGER, PARAMETER :: CLOUDCOEFF_EXP_RELEASE = 1
  INTEGER, PARAMETER :: INVALID_CLOUDCOEFF_EXP = 0
  INTEGER, PARAMETER :: CRTM_EXP_CLOUDCOEFF    = 3    ! distinct from MIE_TAMU(1)/DDA_ARTS(2)


  ! ---------------------------------------
  ! Experimental CloudCoeff data definition
  ! ---------------------------------------
  TYPE :: CloudCoeff_Exp_type
    ! Release/version + allocation state
    INTEGER(Long) :: Release = CLOUDCOEFF_EXP_RELEASE
    INTEGER(Long) :: Version = 1
    LOGICAL       :: Is_Allocated = .FALSE.
    CHARACTER(SL) :: Scheme = 'CRTM-Exp'
    ! Dimensions
    INTEGER(Long) :: n_Frequency      = 0
    INTEGER(Long) :: n_Temperature    = 0
    INTEGER(Long) :: n_Mu             = 0
    INTEGER(Long) :: n_Dm             = 0
    INTEGER(Long) :: n_Habit          = 0
    INTEGER(Long) :: n_Legendre       = 0   ! L_max (full expansion length)
    INTEGER(Long) :: n_Phase_Elements = 0   ! up to 6
    ! Axes
    REAL(Double),  ALLOCATABLE :: Frequency(:)        ! n_Frequency  (GHz)
    REAL(Double),  ALLOCATABLE :: Temperature(:)      ! n_Temperature (K)
    REAL(Double),  ALLOCATABLE :: Mu(:)               ! n_Mu
    REAL(Double),  ALLOCATABLE :: Dm(:)               ! n_Dm (microns)
    ! Habit metadata
    INTEGER(Long), ALLOCATABLE :: Habit_Id(:)         ! n_Habit (CRTM cloud-type integer)
    INTEGER(Long), ALLOCATABLE :: Habit_Phase(:)      ! n_Habit (0=liquid,1=frozen)
    REAL(Double),  ALLOCATABLE :: mD_a(:)             ! n_Habit (mass-dimension prefactor)
    REAL(Double),  ALLOCATABLE :: mD_b(:)             ! n_Habit (mass-dimension exponent)
    CHARACTER(SL), ALLOCATABLE :: Habit_Name(:)       ! n_Habit
    ! Bulk optics : (Frequency, Temperature, Mu, Dm, Habit)
    REAL(Double),  ALLOCATABLE :: ke(:,:,:,:,:)
    REAL(Double),  ALLOCATABLE :: ka(:,:,:,:,:)
    REAL(Double),  ALLOCATABLE :: g (:,:,:,:,:)
    REAL(Double),  ALLOCATABLE :: kb(:,:,:,:,:)
    ! Per-entry effective truncation order : (Frequency, Temperature, Mu, Dm, Habit)
    INTEGER(Long), ALLOCATABLE :: n_Legendre_Eff(:,:,:,:,:)
    ! GSF expansion : (Phase_Elements, Legendre, Frequency, Temperature, Mu, Dm, Habit)
    REAL(Double),  ALLOCATABLE :: pcoeff(:,:,:,:,:,:,:)
  END TYPE CloudCoeff_Exp_type


CONTAINS


  ! Is the structure allocated?
  ELEMENTAL FUNCTION CloudCoeff_Exp_Associated( self ) RESULT( status )
    TYPE(CloudCoeff_Exp_type), INTENT(IN) :: self
    LOGICAL :: status
    status = self%Is_Allocated
  END FUNCTION CloudCoeff_Exp_Associated


  ! Deallocate
  ELEMENTAL SUBROUTINE CloudCoeff_Exp_Destroy( self )
    TYPE(CloudCoeff_Exp_type), INTENT(OUT) :: self
    self%Is_Allocated = .FALSE.
    self%n_Frequency=0; self%n_Temperature=0; self%n_Mu=0; self%n_Dm=0
    self%n_Habit=0; self%n_Legendre=0; self%n_Phase_Elements=0
  END SUBROUTINE CloudCoeff_Exp_Destroy


  ! Allocate
  SUBROUTINE CloudCoeff_Exp_Create( self, n_Frequency, n_Temperature, n_Mu, &
                                    n_Dm, n_Habit, n_Legendre, n_Phase_Elements )
    TYPE(CloudCoeff_Exp_type), INTENT(OUT) :: self
    INTEGER, INTENT(IN) :: n_Frequency, n_Temperature, n_Mu, n_Dm, &
                           n_Habit, n_Legendre, n_Phase_Elements
    INTEGER :: alloc_stat

    IF ( n_Frequency < 1 .OR. n_Temperature < 1 .OR. n_Mu < 1 .OR. n_Dm < 1 .OR. &
         n_Habit < 1 .OR. n_Legendre < 1 .OR. n_Phase_Elements < 1 ) RETURN

    ALLOCATE( self%Frequency(n_Frequency), &
              self%Temperature(n_Temperature), &
              self%Mu(n_Mu), &
              self%Dm(n_Dm), &
              self%Habit_Id(n_Habit), &
              self%Habit_Phase(n_Habit), &
              self%mD_a(n_Habit), &
              self%mD_b(n_Habit), &
              self%Habit_Name(n_Habit), &
              self%ke(n_Frequency,n_Temperature,n_Mu,n_Dm,n_Habit), &
              self%ka(n_Frequency,n_Temperature,n_Mu,n_Dm,n_Habit), &
              self%g (n_Frequency,n_Temperature,n_Mu,n_Dm,n_Habit), &
              self%kb(n_Frequency,n_Temperature,n_Mu,n_Dm,n_Habit), &
              self%n_Legendre_Eff(n_Frequency,n_Temperature,n_Mu,n_Dm,n_Habit), &
              self%pcoeff(n_Phase_Elements,n_Legendre,n_Frequency,n_Temperature,n_Mu,n_Dm,n_Habit), &
              STAT = alloc_stat )
    IF ( alloc_stat /= 0 ) RETURN

    self%n_Frequency      = n_Frequency
    self%n_Temperature    = n_Temperature
    self%n_Mu             = n_Mu
    self%n_Dm             = n_Dm
    self%n_Habit          = n_Habit
    self%n_Legendre       = n_Legendre
    self%n_Phase_Elements = n_Phase_Elements

    self%Frequency = ZERO; self%Temperature = ZERO; self%Mu = ZERO; self%Dm = ZERO
    self%Habit_Id = 0; self%Habit_Phase = 0; self%mD_a = ZERO; self%mD_b = ZERO
    self%Habit_Name = ' '
    self%ke = ZERO; self%ka = ZERO; self%g = ZERO; self%kb = ZERO
    self%n_Legendre_Eff = 0; self%pcoeff = ZERO

    self%Is_Allocated = .TRUE.
  END SUBROUTINE CloudCoeff_Exp_Create


  ! Inspect
  SUBROUTINE CloudCoeff_Exp_Inspect( self )
    TYPE(CloudCoeff_Exp_type), INTENT(IN) :: self
    INTEGER :: i
    WRITE(*,'(1x,"CloudCoeff_Exp OBJECT  scheme=",a)') TRIM(self%Scheme)
    WRITE(*,'(3x,"Release.Version :",i0,".",i0)') self%Release, self%Version
    IF ( .NOT. self%Is_Allocated ) THEN
      WRITE(*,'(3x,"(not allocated)")'); RETURN
    END IF
    WRITE(*,'(3x,"n_Frequency      :",i0)') self%n_Frequency
    WRITE(*,'(3x,"n_Temperature    :",i0)') self%n_Temperature
    WRITE(*,'(3x,"n_Mu             :",i0)') self%n_Mu
    WRITE(*,'(3x,"n_Dm             :",i0)') self%n_Dm
    WRITE(*,'(3x,"n_Habit          :",i0)') self%n_Habit
    WRITE(*,'(3x,"n_Legendre (max) :",i0)') self%n_Legendre
    WRITE(*,'(3x,"n_Phase_Elements :",i0)') self%n_Phase_Elements
    WRITE(*,'(3x,"Habits:")')
    DO i = 1, self%n_Habit
      WRITE(*,'(5x,i3,1x,a,"  phase=",i0,"  m-D a,b=",es9.2,1x,f5.2)') &
        self%Habit_Id(i), TRIM(self%Habit_Name(i)), self%Habit_Phase(i), self%mD_a(i), self%mD_b(i)
    END DO
  END SUBROUTINE CloudCoeff_Exp_Inspect


  ! One-line info string
  SUBROUTINE CloudCoeff_Exp_Info( self, Info )
    TYPE(CloudCoeff_Exp_type), INTENT(IN)  :: self
    CHARACTER(*),              INTENT(OUT) :: Info
    WRITE(Info,'("CloudCoeff_Exp R.V=",i0,".",i0,&
               &" nFreq=",i0," nT=",i0," nMu=",i0," nDm=",i0," nHabit=",i0,&
               &" Lmax=",i0," nPhase=",i0)') &
      self%Release, self%Version, self%n_Frequency, self%n_Temperature, self%n_Mu, &
      self%n_Dm, self%n_Habit, self%n_Legendre, self%n_Phase_Elements
  END SUBROUTINE CloudCoeff_Exp_Info


  SUBROUTINE CloudCoeff_Exp_DefineVersion( Id )
    CHARACTER(*), INTENT(OUT) :: Id
    Id = MODULE_VERSION_ID
  END SUBROUTINE CloudCoeff_Exp_DefineVersion

END MODULE CloudCoeff_Exp_Define
