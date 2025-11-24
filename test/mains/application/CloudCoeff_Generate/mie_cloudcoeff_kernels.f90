! File: mie_cloudcoeff_kernels.f90
!
! Uses Type_Kinds::Double (no dp alias). Explicit interfaces live in the
! specification part (above CONTAINS). MAXLEG = 39 to match authoritative file.

module mie_cloudcoeff_kernels
  use Type_Kinds, only: Double
  implicit none
  private
  public :: compute_optics_mw, compute_optics_ir, MAXLEG

  integer, parameter :: wp = Double
  integer, parameter :: MAXLEG = 39

  interface
    subroutine Dmie(nshp, wavelength, rhoavx, mspherex, mcore, corefrac, &
                    numDx, MAXLEG_in, ad, bd, alpha, gamma, IWCout,      &
                    volext, mext, volscat, albedo, volback, mback,       &
                    volasym, nlegen, coef1, coef2, coef3, coef4)
      use Type_Kinds, only: Double
      implicit none
      integer,      intent(in)  :: nshp, numDx, MAXLEG_in
      real(Double), intent(in)  :: wavelength, rhoavx, corefrac, ad, bd, alpha, gamma
      complex(Double), intent(in) :: mspherex, mcore
      real(Double), intent(out) :: IWCout, volext, mext, volscat, albedo, volback, mback, volasym
      integer,      intent(out) :: nlegen
      real(Double), intent(out) :: coef1(MAXLEG_in+1), coef2(MAXLEG_in+1), coef3(MAXLEG_in+1), coef4(MAXLEG_in+1)
    end subroutine Dmie

    subroutine cdmx(f1,f2,tc,nu,mflag,Kfactor,nav)
      use Type_Kinds, only: Double
      implicit none
      integer,      intent(in)  :: mflag
      real(Double), intent(in)  :: f1, f2, tc, nu
      real(Double), intent(out) :: Kfactor
      complex(Double), intent(out) :: nav
    end subroutine cdmx
  end interface

contains

  subroutine compute_optics_mw(freq_GHz, temp_K, phase_is_liquid, mixflag,  &
                               reff_um, rho_bulk_kgm3,                       &
                               nlegen, ke, w, g, p1, p2, p3, p4)
    implicit none
    real(wp),    intent(in)  :: freq_GHz, temp_K, reff_um, rho_bulk_kgm3
    logical,     intent(in)  :: phase_is_liquid
    integer,     intent(in)  :: mixflag
    integer,     intent(out) :: nlegen
    real(wp),    intent(out) :: ke, w, g
    real(wp),    intent(out) :: p1(MAXLEG+1), p2(MAXLEG+1), p3(MAXLEG+1), p4(MAXLEG+1)

    integer                 :: nshp, numDx
    real(Double)            :: wavelength_m, ad, bd, alpha_par, gamma_par
    real(Double)            :: IWCout, volext, mext, volscat, albedo, volback, mback, volasym
    complex(Double)         :: msphere, mcore
    real(wp)                :: temp_C, corefrac, rhoav_gcm3
    real(wp)                :: f1, f2

    call psd_params_from_reff(reff_um, rho_bulk_kgm3, ad, bd, alpha_par, gamma_par)

    temp_C = temp_K - 273.15_wp
    f1     = merge(1.0_wp, 0.0_wp, .not. phase_is_liquid)   ! ice fraction
    f2     = merge(1.0_wp, 0.0_wp,        phase_is_liquid)  ! water fraction
    call cdmx(f1, f2, temp_C, freq_GHz, mixflag, Kfactor = IWCout, nav = msphere)

    mcore    = msphere
    corefrac = 0.0_wp

    wavelength_m = 299792458.0_wp / (freq_GHz * 1.0e9_wp)   ! [m]
    rhoav_gcm3   = max(1.0e-6_wp, rho_bulk_kgm3 * 1.0e-3_wp)

    nshp  = 1
    numDx = 200

    p1 = 0.0_wp; p2 = 0.0_wp; p3 = 0.0_wp; p4 = 0.0_wp
    nlegen = 0; ke = 0.0_wp; w = 0.0_wp; g = 0.0_wp

    call Dmie(nshp, wavelength_m, rhoav_gcm3, msphere, mcore, corefrac, &
              numDx, MAXLEG, ad, bd, alpha_par, gamma_par,              &
              IWCout,  volext, mext, volscat, albedo, volback, mback,   &
              volasym, nlegen, p1, p2, p3, p4)

    ke = mext
    w  = albedo
    if (volscat > 0.0_wp) then
      g = volasym / volscat
    else
      g = 0.0_wp
    end if
  end subroutine compute_optics_mw


  subroutine compute_optics_ir(wavenumber_cm, temp_K, phase_is_liquid, mixflag,  &
                               reff_um, rho_bulk_kgm3,                           &
                               nlegen, ke, w, g, p1, p2, p3, p4)
    implicit none
    real(wp),    intent(in)  :: wavenumber_cm, temp_K, reff_um, rho_bulk_kgm3
    logical,     intent(in)  :: phase_is_liquid
    integer,     intent(in)  :: mixflag
    integer,     intent(out) :: nlegen
    real(wp),    intent(out) :: ke, w, g
    real(wp),    intent(out) :: p1(MAXLEG+1), p2(MAXLEG+1), p3(MAXLEG+1), p4(MAXLEG+1)

    integer                 :: nshp, numDx
    real(Double)            :: wavelength_m, ad, bd, alpha_par, gamma_par
    real(Double)            :: IWCout, volext, mext, volscat, albedo, volback, mback, volasym
    complex(Double)         :: msphere, mcore
    real(wp)                :: temp_C, corefrac, rhoav_gcm3
    real(wp)                :: f1, f2

    call psd_params_from_reff(reff_um, rho_bulk_kgm3, ad, bd, alpha_par, gamma_par)

    temp_C = temp_K - 273.15_wp
    f1     = merge(1.0_wp, 0.0_wp, .not. phase_is_liquid)
    f2     = merge(1.0_wp, 0.0_wp,        phase_is_liquid)

    call cdmx(f1, f2, temp_C, 2.99792458e-3_wp * wavenumber_cm, mixflag, Kfactor = IWCout, nav = msphere)

    mcore    = msphere
    corefrac = 0.0_wp

    wavelength_m = 1.0_wp / (wavenumber_cm * 100.0_wp)      ! [m]
    rhoav_gcm3   = max(1.0e-6_wp, rho_bulk_kgm3 * 1.0e-3_wp)

    nshp  = 1
    numDx = 200

    p1 = 0.0_wp; p2 = 0.0_wp; p3 = 0.0_wp; p4 = 0.0_wp
    nlegen = 0; ke = 0.0_wp; w = 0.0_wp; g = 0.0_wp

    call Dmie(nshp, wavelength_m, rhoav_gcm3, msphere, mcore, corefrac, &
              numDx, MAXLEG, ad, bd, alpha_par, gamma_par,              &
              IWCout,  volext, mext, volscat, albedo, volback, mback,   &
              volasym, nlegen, p1, p2, p3, p4)

    ke = mext
    w  = albedo
    if (volscat > 0.0_wp) then
      g = volasym / volscat
    else
      g = 0.0_wp
    end if
  end subroutine compute_optics_ir


  subroutine psd_params_from_reff(reff_um, rho_bulk_kgm3, ad, bd, alpha_par, gamma_par)
    implicit none
    real(wp),    intent(in)  :: reff_um, rho_bulk_kgm3
    real(Double),intent(out) :: ad, bd, alpha_par, gamma_par
    alpha_par = 2.0_Double
    gamma_par = 1.0_Double
    bd        = max(1.0e-12_Double, 3.0_Double / max(1.0e-6_Double, real(reff_um,Double)))
    ad        = 1.0_Double
  end subroutine psd_params_from_reff

end module mie_cloudcoeff_kernels
