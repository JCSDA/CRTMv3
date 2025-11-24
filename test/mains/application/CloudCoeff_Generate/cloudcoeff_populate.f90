! Fix: put the OpenMP clauses on a single line and ensure the directive
! is immediately followed by the DO nest. Some compilers choke on
! continued directive lines. Also, no executable statements between the
! directive and the first DO.

module cloudcoeff_populate
  use Type_Kinds,        only: Double
  use CloudCoeff_Define , only: CloudCoeff_type
  use mie_cloudcoeff_kernels, only: compute_optics_mw, compute_optics_ir, MAXLEG
  implicit none
  private
  public :: fill_cloudcoeff_all
  integer, parameter :: wp = Double
contains

  subroutine fill_cloudcoeff_all(CC, mixflag)
    type(CloudCoeff_type), intent(inout) :: CC
    integer,               intent(in)    :: mixflag

    integer :: iF, iR, iT, iD, nL, nleg
    real(wp) :: kev, wv, gv
    real(wp) :: p1(MAXLEG+1), p2(MAXLEG+1), p3(MAXLEG+1), p4(MAXLEG+1)

    nL = min(CC%n_Legendre_Terms, MAXLEG)

    !========================
    ! MW: SOLID  I1 x I2 x I6
    !========================
!$omp parallel do 
    do iF = 1, CC%n_MW_Frequencies
      print *, 'S_mw:', iF, CC%n_MW_Frequencies
      do iR = 1, CC%n_MW_Radii
        do iD = 1, CC%n_Densities
          call compute_optics_mw( CC%Frequency_MW(iF), 273.15_wp, .false., mixflag, &
                                  CC%Reff_MW(iR), CC%Density(iD),                    &
                                  nleg, kev, wv, gv, p1, p2, p3, p4 )
          CC%ke_S_MW(iF,iR,iD) = kev
          CC%w_S_MW (iF,iR,iD) = wv
          CC%g_S_MW (iF,iR,iD) = gv
          CC%pcoeff_S_MW(iF,iR,iD,1:nL,1) = p1(1:nL)/4.0_wp
        end do
      end do
    end do
!$omp end parallel do

    !========================
    ! MW: LIQUID I1 x I2 x I5
    !========================
!$omp parallel do 
    do iF = 1, CC%n_MW_Frequencies
      print *, 'L_mw:', iF, CC%n_MW_Frequencies
      do iR = 1, CC%n_MW_Radii
        do iT = 1, CC%n_Temperatures
          call compute_optics_mw( CC%Frequency_MW(iF), CC%Temperature(iT), .true., mixflag, &
                                  CC%Reff_MW(iR), 1000.0_wp,                                 &
                                  nleg, kev, wv, gv, p1, p2, p3, p4 )
          CC%ke_L_MW(iF,iR,iT) = kev
          CC%w_L_MW (iF,iR,iT) = wv
          CC%g_L_MW (iF,iR,iT) = gv
          CC%pcoeff_L_MW(iF,iR,iT,1:nL,1) = p1(1:nL)/4.0_wp
        end do
      end do
    end do
!$omp end parallel do

    !========================
    ! IR: 0th density = liquid; 1..n_Densities = solid
    !========================
!$omp parallel do 
    do iF = 1, CC%n_IR_Frequencies
      print *, 'IR:', iF, CC%n_IR_Frequencies, CC%Frequency_IR(1:CC%n_IR_Frequencies)
      do iR = 1, CC%n_IR_Radii
        call compute_optics_ir( CC%Frequency_IR(iF), 273.15_wp, .true., mixflag,  &
                                CC%Reff_IR(iR), 1000.0_wp,                          &
                                nleg, kev, wv, gv, p1, p2, p3, p4 )
        CC%ke_IR(iF,iR,0) = kev
        CC%w_IR (iF,iR,0) = wv
        CC%g_IR (iF,iR,0) = gv
        CC%pcoeff_IR(iF,iR,0,1:nL) = p1(1:nL)/4.0_wp

        do iD = 1, CC%n_Densities
          call compute_optics_ir( CC%Frequency_IR(iF), 263.15_wp, .false., mixflag, &
                                  CC%Reff_IR(iR), CC%Density(iD),                    &
                                  nleg, kev, wv, gv, p1, p2, p3, p4 )
          CC%ke_IR(iF,iR,iD) = kev
          CC%w_IR (iF,iR,iD) = wv
          CC%g_IR (iF,iR,iD) = gv
          CC%pcoeff_IR(iF,iR,iD,1:nL) = p1(1:nL)/4.0_wp
        end do
      end do
    end do
!$omp end parallel do

  end subroutine fill_cloudcoeff_all
end module cloudcoeff_populate

!!$
!!$
!!$! File: cloudcoeff_populate.f90
!!$!
!!$! Conforms to your CloudCoeff_type:
!!$!   ke_L_MW      (I1 x I2 x I5)
!!$!   ke_S_MW      (I1 x I2 x I6)
!!$!   ke_IR        (I3 x I4 x 0:I6)  ! 0th density = liquid
!!$! Phase-function packing: use coef1 only, /4 scaling, truncate to n_Legendre_Terms.
!!$
!!$module cloudcoeff_populate
!!$  use Type_Kinds,        only: Double
!!$  use CloudCoeff_Define , only: CloudCoeff_type
!!$  use mie_cloudcoeff_kernels, only: compute_optics_mw, compute_optics_ir, MAXLEG
!!$  implicit none
!!$  private
!!$  public :: fill_cloudcoeff_all
!!$
!!$  integer, parameter :: wp = Double
!!$
!!$contains
!!$
!!$  subroutine fill_cloudcoeff_all(CC, mixflag)
!!$    type(CloudCoeff_type), intent(inout) :: CC
!!$    integer,               intent(in)    :: mixflag
!!$
!!$    integer :: iF, iR, iT, iD, nL, nleg
!!$    real(wp) :: kev, wv, gv
!!$    real(wp), allocatable :: p1(:), p2(:), p3(:), p4(:)
!!$
!!$    nL = min(CC%n_Legendre_Terms, MAXLEG)
!!$    allocate(p1(MAXLEG+1), p2(MAXLEG+1), p3(MAXLEG+1), p4(MAXLEG+1))
!!$
!!$    !------------------
!!$    ! MW: SOLID  I1 x I2 x I6
!!$    !------------------
!!$    do iF = 1, CC%n_MW_Frequencies
!!$       print *, 'S_mw:', iF, CC%n_MW_Frequencies
!!$      do iR = 1, CC%n_MW_Radii
!!$        do iD = 1, CC%n_Densities
!!$          call compute_optics_mw( CC%Frequency_MW(iF), 273.15_wp, .false., mixflag, &
!!$                                  CC%Reff_MW(iR), CC%Density(iD),                    &
!!$                                  nleg, kev, wv, gv, p1, p2, p3, p4 )
!!$          CC%ke_S_MW(iF,iR,iD) = kev
!!$          CC%w_S_MW (iF,iR,iD) = wv
!!$          CC%g_S_MW (iF,iR,iD) = gv
!!$          call pack_leg_phase(CC%pcoeff_S_MW(iF,iR,iD,1:nL,1:CC%n_Phase_Elements), p1, nL)
!!$        end do
!!$      end do
!!$    end do
!!$
!!$    !------------------
!!$    ! MW: LIQUID I1 x I2 x I5
!!$    !------------------
!!$    do iF = 1, CC%n_MW_Frequencies
!!$       print *, 'L_mw:', iF, CC%n_MW_Frequencies
!!$      do iR = 1, CC%n_MW_Radii
!!$        do iT = 1, CC%n_Temperatures
!!$          call compute_optics_mw( CC%Frequency_MW(iF), CC%Temperature(iT), .true., mixflag, &
!!$                                  CC%Reff_MW(iR), 1000.0_wp,                                 &
!!$                                  nleg, kev, wv, gv, p1, p2, p3, p4 )
!!$          CC%ke_L_MW(iF,iR,iT) = kev
!!$          CC%w_L_MW (iF,iR,iT) = wv
!!$          CC%g_L_MW (iF,iR,iT) = gv
!!$          call pack_leg_phase(CC%pcoeff_L_MW(iF,iR,iT,1:nL,1:CC%n_Phase_Elements), p1, nL)
!!$        end do
!!$      end do
!!$    end do
!!$
!!$    !------------------
!!$    ! IR: 0th density = LIQUID; 1..n_Densities = SOLID
!!$    ! ke/w/g: I3 x I4 x 0:I6
!!$    ! pcoeff : I3 x I4 x 0:I6 x I7
!!$    !------------------
!!$    do iF = 1, CC%n_IR_Frequencies
!!$      print *, 'IR:', iF, CC%n_MW_Frequencies
!!$      do iR = 1, CC%n_IR_Radii
!!$
!!$        ! liquid at density index 0
!!$        call compute_optics_ir( CC%Frequency_IR(iF), 273.15_wp, .true., mixflag, &
!!$                                CC%Reff_IR(iR), 1000.0_wp,                         &
!!$                                nleg, kev, wv, gv, p1, p2, p3, p4 )
!!$        CC%ke_IR(iF,iR,0) = kev
!!$        CC%w_IR (iF,iR,0) = wv
!!$        CC%g_IR (iF,iR,0) = gv
!!$        call pack_leg(CC%pcoeff_IR(iF,iR,0,1:nL), p1, nL)
!!$
!!$        ! solids 1..n_Densities
!!$        do iD = 1, CC%n_Densities
!!$          call compute_optics_ir( CC%Frequency_IR(iF), 263.15_wp, .false., mixflag, &
!!$                                  CC%Reff_IR(iR), CC%Density(iD),                    &
!!$                                  nleg, kev, wv, gv, p1, p2, p3, p4 )
!!$          CC%ke_IR(iF,iR,iD) = kev
!!$          CC%w_IR (iF,iR,iD) = wv
!!$          CC%g_IR (iF,iR,iD) = gv
!!$          call pack_leg(CC%pcoeff_IR(iF,iR,iD,1:nL), p1, nL)
!!$        end do
!!$
!!$      end do
!!$    end do
!!$
!!$    deallocate(p1,p2,p3,p4)
!!$
!!$  contains
!!$
!!$    subroutine pack_leg_phase(dst_leg_phase, p1, nL)
!!$      real(wp), intent(inout) :: dst_leg_phase(:,:)
!!$      real(wp), intent(in)    :: p1(:)
!!$      integer , intent(in)    :: nL
!!$      integer :: n
!!$      do n = 1, nL
!!$        dst_leg_phase(n,1) = p1(n)/4.0_wp
!!$      end do
!!$    end subroutine pack_leg_phase
!!$
!!$    subroutine pack_leg(dst_leg, p1, nL)
!!$      real(wp), intent(inout) :: dst_leg(:)
!!$      real(wp), intent(in)    :: p1(:)
!!$      integer , intent(in)    :: nL
!!$      integer :: n
!!$      do n = 1, nL
!!$        dst_leg(n) = p1(n)/4.0_wp
!!$      end do
!!$    end subroutine pack_leg
!!$
!!$  end subroutine fill_cloudcoeff_all
!!$
!!$end module cloudcoeff_populate
