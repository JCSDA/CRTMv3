! File: miecalc_safe.f90
module miecalc_safe
  use Type_Kinds, only: FP
  implicit none
  private
  public :: Dmiecalc
contains

  subroutine Dmiecalc(nterms, x, mn, a, b, istat)
    implicit none
    integer,      intent(inout) :: nterms
    real(FP),     intent(in)    :: x
    complex(FP),  intent(in)    :: mn
    complex(FP),  intent(out)   :: a(:), b(:)
    integer,      intent(out), optional :: istat

    integer, parameter :: maxterms_cap = 1000000
    integer :: ncap, nn, n
    real(FP)    :: psin, psim, chin, chim, tmp
    complex(FP) :: m, y, ctmp, xin, xim
    complex(FP), allocatable :: D(:)

    if (present(istat)) istat = 0
    if (x <= 0.0_FP) then
      nterms = 0
      if (size(a) > 0) a = (0.0_FP, 0.0_FP)
      if (size(b) > 0) b = (0.0_FP, 0.0_FP)
      if (present(istat)) istat = 1
      return
    end if

    if (nterms <= 0) nterms = int(x + 4.0_FP * x**(1.0_FP/3.0_FP) + 2.0_FP)
    ncap = min(nterms, size(a), size(b), maxterms_cap)
    nterms = ncap
    if (ncap <= 0) then
      if (present(istat)) istat = 2
      return
    end if

    allocate(D(ncap + 15))
    D = (0.0_FP, 0.0_FP)

    m  = conjg(mn)
    y  = m * x
    nn = ncap + 15
    D(nn) = (0.0_FP, 0.0_FP)
    do n = nn, 2, -1
      D(n-1) = real(n,FP)/y - 1.0_FP / ( D(n) + real(n,FP)/y )
    end do

    psim = cos(x);  psin = sin(x)
    chim = -sin(x); chin = cos(x)

    a = (0.0_FP, 0.0_FP)
    b = (0.0_FP, 0.0_FP)

    do n = 1, ncap
      tmp  = psin
      psin = ( real(2*n-1,FP)/x )*psin - psim
      psim = tmp
      tmp  = chin
      chin = ( real(2*n-1,FP)/x )*chin - chim
      chim = tmp

      xin  = cmplx(psin, -chin, kind=FP)
      xim  = cmplx(psim, -chim, kind=FP)

      ctmp = D(n)/m + real(n,FP)/x
      a(n) = (ctmp*psin - psim) / (ctmp*xin - xim)

      ctmp = m*D(n)  + real(n,FP)/x
      b(n) = (ctmp*psin - psim) / (ctmp*xin - xim)
    end do

    deallocate(D)
  end subroutine Dmiecalc

end module miecalc_safe
