!
! Type_Kinds
!
! Module to hold specification kinds for variable declaration, as well as 
! associated descriptors.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 12-Jun-2000
!                       paul.vandelst@noaa.gov
!

MODULE CSEM_Type_Kinds

  ! ---------------------------
  ! Disable all implicit typing
  ! ---------------------------
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything is private by default
  PRIVATE
  ! The integer types
  PUBLIC :: CSEM_Byte   
  PUBLIC :: CSEM_Short  
  PUBLIC :: CSEM_Long  
  PUBLIC :: CSEM_LLong  
  PUBLIC :: CSEM_IP_Kind  ! Default integer set by IIP
  PUBLIC :: CSEM_IP       ! Aliases for IP_Kind
  ! The floating point types
  PUBLIC :: CSEM_Single 
  PUBLIC :: CSEM_Double 
  PUBLIC :: CSEM_Quad    
  PUBLIC :: CSEM_FP_Kind  ! Default integer set by IFP
  PUBLIC :: CSEM_FP       ! Aliases for FP_Kind


  ! -------------------------------------------------------------------
  ! THE DEFAULT INTEGER INDEX. Change the value of IIP for the required
  ! integer kind. The following chart details the correspondence:
  !
  !    IIP        INTEGER(IP)
  !  ==============================
  !     1        Byte 
  !     2       Short (2 bytes)
  !     3        Long (4 bytes)
  !     4       LLong (8 bytes)  **IF AVAILABLE, Long OTHERWISE**
  !
  ! -------------------------------------------------------------------
  INTEGER, PARAMETER :: IIP = 3  ! 1=Byte, 2=Short, 3=Long, 4=LLong


  ! -------------------------------------------------------------------
  ! THE DEFAULT FLOATING POINT INDEX. Change the value of IFP for the
  ! required floating point kind. The following chart details the
  ! correspondence:
  !
  !    IFP          REAL(FP)
  !  ==============================
  !     1       Single (4  bytes)
  !     2       Double (8  bytes)
  !     3       Quad   (16 bytes)  **IF AVAILABLE, Double OTHERWISE**
  !
  ! -------------------------------------------------------------------
  INTEGER, PARAMETER :: IFP = 2  ! 1=Single, 2=Double, 3=Quad


  ! -------------------
  ! Integer definitions
  ! -------------------
  ! Integer types
  INTEGER, PARAMETER :: Byte    = SELECTED_INT_KIND(1)   ! Byte  integer
  INTEGER, PARAMETER :: Short   = SELECTED_INT_KIND(4)   ! Short integer
  INTEGER, PARAMETER :: Long    = SELECTED_INT_KIND(8)   ! Long  integer
  INTEGER, PARAMETER :: LLong   = SELECTED_INT_KIND(16)  ! LLong integer

  ! Define arrays for default definition
  INTEGER, PARAMETER :: N_IP = 4
  INTEGER, PARAMETER, DIMENSION(N_IP) :: IP_KIND_TYPES = (/ Byte,  &
                                                            Short, &
                                                            Long,  &
                                                            LLong  /) 
  ! Default values
  INTEGER, PARAMETER :: IP_Kind   =   IP_KIND_TYPES(IIP)
  INTEGER, PARAMETER :: IP        =   IP_Kind

  INTEGER, PARAMETER :: CSEM_Byte     =   Byte 
  INTEGER, PARAMETER :: CSEM_Short    =   Short 
  INTEGER, PARAMETER :: CSEM_Long     =   Long
  INTEGER, PARAMETER :: CSEM_LLong    =   LLong
  INTEGER, PARAMETER :: CSEM_IP_Kind  =   IP_Kind
  INTEGER, PARAMETER :: CSEM_IP       =   IP

  ! --------------------------
  ! Floating point definitions
  ! --------------------------
  ! Floating point types
  INTEGER, PARAMETER :: Single = SELECTED_REAL_KIND(6)  ! Single precision
  INTEGER, PARAMETER :: Double = SELECTED_REAL_KIND(15) ! Double precision
  INTEGER, PARAMETER :: Quad   = SELECTED_REAL_KIND(20) ! Quad precision

  ! Define arrays for default definition
  INTEGER, PARAMETER :: N_FP = 3
  INTEGER, PARAMETER, DIMENSION(N_FP) :: FP_KIND_TYPES = (/ Single, &
                                                            Double, &
                                                            Quad    /) 
  ! Default values
  INTEGER, PARAMETER :: FP_Kind   = FP_KIND_TYPES(IFP)
  INTEGER, PARAMETER :: FP        = FP_Kind
  
  INTEGER, PARAMETER :: CSEM_Single   = Single
  INTEGER, PARAMETER :: CSEM_Double   = Double
  INTEGER, PARAMETER :: CSEM_Quad     = Quad  
  INTEGER, PARAMETER :: CSEM_FP_Kind  = FP_Kind ! Default integer set by IFP
  INTEGER, PARAMETER :: CSEM_FP       = FP      ! Aliases for FP_Kind

  
END MODULE CSEM_Type_Kinds
