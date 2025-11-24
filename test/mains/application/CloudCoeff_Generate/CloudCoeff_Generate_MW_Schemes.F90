! CloudCoeff_Generate_MW_Schemes.F90
! Single-source file: all helpers are internal; no external module deps.
! Variable orders in the NetCDF file use "reversed" axis order (F,R,...) per request.

module util_kinds
  use iso_fortran_env, only: real64
  implicit none
  integer, parameter :: rk = real64
end module util_kinds

module scheme_mass
  use util_kinds
  implicit none
contains
  pure real(rk) function mass_g_density1(d_um) result(mg)
    real(rk), intent(in) :: d_um
    real(rk) :: d_cm, vcirc_cm3
    d_cm = d_um*1.0e-4_rk
    vcirc_cm3 = (acos(-1.0_rk)/6.0_rk)*d_cm**3
    mg = 0.1_rk*vcirc_cm3
  end function mass_g_density1

  pure real(rk) function mass_g_density2(d_um) result(mg)
    real(rk), intent(in) :: d_um
    real(rk) :: d_mm
    d_mm = d_um*1.0e-3_rk
    mg = 6.9e-5_rk*d_mm**2.0_rk
  end function mass_g_density2

  pure real(rk) function mass_g_density3(d_um) result(mg)
    real(rk), intent(in) :: d_um
    real(rk) :: d_mm
    d_mm = d_um*1.0e-3_rk
    mg = 5.436e-5_rk*d_mm**2.05_rk
  end function mass_g_density3

  pure real(rk) function mass_g_density4(d_um) result(mg)
    real(rk), intent(in) :: d_um
    real(rk) :: d_mm
    d_mm = d_um*1.0e-3_rk
    mg = 8.90e-5_rk*d_mm**2.1_rk
  end function mass_g_density4

  pure real(rk) function reff_um_from_mass_g(mg) result(r_um)
    real(rk), intent(in) :: mg
    real(rk), parameter :: rho_ice_g_cm3 = 0.917_rk
    real(rk) :: r_cm
    r_cm = (3.0_rk*mg/(4.0_rk*acos(-1.0_rk)*rho_ice_g_cm3))**(1.0_rk/3.0_rk)
    r_um = r_cm*1.0e4_rk
  end function reff_um_from_mass_g
end module scheme_mass

module text_io
  use util_kinds
  implicit none
contains
  subroutine read_vector_file(filename, x, nread, allow_comments)
    character(*), intent(in)  :: filename
    real(rk),     intent(out) :: x(:)
    integer,      intent(out) :: nread
    logical,      intent(in),  optional :: allow_comments
    logical :: allowc
    integer :: u, ios, nmax
    character(len=8192) :: line
    real(rk) :: v
    allowc = .false.; if (present(allow_comments)) allowc = allow_comments
    nread = 0; nmax = size(x)
    open(newunit=u, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'failed to open file'
    do
      read(u,'(A)',iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      if (allowc) then
        if (line(1:1) == '#' .or. line(1:1) == '!') cycle
      end if
      read(line,*,iostat=ios) v
      if (ios /= 0) cycle
      nread = nread + 1
      if (nread <= nmax) then
        x(nread) = v
      else
        exit
      end if
    end do
    close(u)
  end subroutine read_vector_file
end module text_io

module legendre_tools
  use util_kinds
  implicit none
contains
  subroutine trapezoid_mu_weights(theta, w)
    real(rk), intent(in)  :: theta(:)
    real(rk), intent(out) :: w(:)
    integer :: n, i
    real(rk) :: dthm, dthp
    n = size(theta)
    w = 0.0_rk
    do i = 1, n
      if (i == 1) then
        dthp = theta(i+1) - theta(i)
        w(i) = sin(theta(i)) * dthp * 0.5_rk
      else if (i == n) then
        dthm = theta(i) - theta(i-1)
        w(i) = sin(theta(i)) * dthm * 0.5_rk
      else
        dthp = theta(i+1) - theta(i)
        dthm = theta(i)   - theta(i-1)
        w(i) = sin(theta(i)) * 0.5_rk * (dthp + dthm)
      end if
    end do
  end subroutine trapezoid_mu_weights

  subroutine legendre_P_l(mu, l, Pl)
    real(rk), intent(in)  :: mu(:)
    integer,  intent(in)  :: l
    real(rk), intent(out) :: Pl(:)
    integer :: n, k
    real(rk), allocatable :: P0(:), P1(:), Pn(:)
    n = size(mu)
    if (l == 0) then
      Pl = 1.0_rk
      return
    else if (l == 1) then
      Pl = mu
      return
    end if
    allocate(P0(n), P1(n), Pn(n))
    P0 = 1.0_rk
    P1 = mu
    do k = 2, l
      Pn = ((2.0_rk*k-1.0_rk)*mu*P1 - (k-1.0_rk)*P0)/real(k, rk)
      P0 = P1
      P1 = Pn
    end do
    Pl = P1
  end subroutine legendre_P_l

  subroutine p11_to_coef(mu, w, P11, nLeg, coef)
    real(rk), intent(in)  :: mu(:), w(:), P11(:)
    integer,  intent(in)  :: nLeg
    real(rk), intent(out) :: coef(0:nLeg-1)
    integer :: l
    real(rk) :: Pl(size(mu)), norm, twoPi
    twoPi = 2.0_rk*acos(-1.0_rk)
    norm  = twoPi * sum(P11 * w)
    if (norm <= 0.0_rk) then
      coef = 0.0_rk
      coef(0) = 1.0_rk
      return
    end if
    do l = 0, nLeg-1
      call legendre_P_l(mu, l, Pl)
      coef(l) = (2.0_rk*l+1.0_rk)/norm * twoPi * sum(P11 * Pl * w)
    end do
    if (coef(0) /= 0.0_rk) coef = coef/coef(0)
  end subroutine p11_to_coef
end module legendre_tools

module attr_tools
  use netcdf
  implicit none
contains
  subroutine put_text_att(ncid, name, val)
    integer,      intent(in) :: ncid
    character(*), intent(in) :: name, val
    integer :: ierr
    character(len=:), allocatable :: sval
    sval = trim(val)
    ierr = nf90_put_att(ncid, NF90_GLOBAL, trim(name), sval)
    if (ierr /= nf90_noerr) stop 'nf90_put_att text failed'
  end subroutine put_text_att

  subroutine put_int_att(ncid, name, ival)
    integer,      intent(in) :: ncid
    character(*), intent(in) :: name
    integer,      intent(in) :: ival
    integer :: ierr
    ierr = nf90_put_att(ncid, NF90_GLOBAL, trim(name), ival)
    if (ierr /= nf90_noerr) stop 'nf90_put_att int failed'
  end subroutine put_int_att
end module attr_tools

module mw_ir_netcdf_writer
  use util_kinds
  use netcdf
  use attr_tools, only: put_text_att, put_int_att
  implicit none
contains
  subroutine write_cloudcoeff_mw_ir(outfile, &
                                    freq_mw, reff_mw, temp_k, density_val, &
                                    ke_S_mw, w_S_mw, g_S_mw, pS_mw_PLTRF, &
                                    freq_ir, reff_ir, &
                                    ke_ir_solid_RxF, w_ir_solid_RxF, g_ir_solid_RxF, pcoeff_ir_solid_LxRxF)
    character(*), intent(in) :: outfile
    real(rk),     intent(in) :: freq_mw(:), reff_mw(:), temp_k(:), density_val
    real(rk),     intent(in) :: ke_S_mw(:,:), w_S_mw(:,:), g_S_mw(:,:)              ! (R,F)
    real(rk),     intent(in) :: pS_mw_PLTRF(:,:,:,:,:)                               ! (P,L,T,R,F)
    real(rk),     intent(in) :: freq_ir(:), reff_ir(:)
    real(rk),     intent(in) :: ke_ir_solid_RxF(:,:), w_ir_solid_RxF(:,:), g_ir_solid_RxF(:,:) ! (R,F)
    real(rk),     intent(in) :: pcoeff_ir_solid_LxRxF(:,:,:)                         ! (L,R,F)

    integer :: ncid, ierr
    integer :: did_t, did_fmw, did_rmw, did_fir, did_rir, did_d, did_dir, did_l, did_p
    integer :: nT, nFmw, nRmw, nFir, nRir, nLeg, nPhase, nDen_mw, nDen_ir
    integer :: vid_T, vid_fmw, vid_fir, vid_rmw, vid_rir, vid_D
    integer :: vid_keL_mw, vid_wL_mw, vid_gL_mw, vid_pL_mw
    integer :: vid_keS_mw, vid_wS_mw, vid_gS_mw, vid_pS_mw
    integer :: vid_ke_ir, vid_w_ir, vid_g_ir, vid_p_ir

    real(rk), allocatable :: density(:)

    ! temp arrays in requested reversed on-disk order
    real(rk), allocatable :: zero_f_r_t(:,:,:)
    real(rk), allocatable :: zero_pL_f_r_t_l_p(:,:,:,:,:)
    real(rk), allocatable :: keS_f_r_d(:,:,:), wS_f_r_d(:,:,:), gS_f_r_d(:,:,:)
    real(rk), allocatable :: pS_f_r_d_l_p(:,:,:,:,:)
    real(rk), allocatable :: keIR_f_r_d(:,:,:), wIR_f_r_d(:,:,:), gIR_f_r_d(:,:,:)
    real(rk), allocatable :: pIR_f_r_d_l(:,:,:,:)

    integer :: f, r, t, l

    nT = size(temp_k); nFmw = size(freq_mw); nRmw = size(reff_mw)
    nFir = size(freq_ir); nRir = size(reff_ir)
    nLeg = size(pcoeff_ir_solid_LxRxF,1)
    nPhase = 1
    nDen_mw = 1
    nDen_ir = 2

    allocate(density(nDen_mw)); density = density_val

    ierr = nf90_create(outfile, NF90_NETCDF4, ncid); if (ierr/=nf90_noerr) stop 'create failed'

    ierr = nf90_def_dim(ncid,'n_Temperatures',   nT,   did_t)
    ierr = nf90_def_dim(ncid,'n_MW_Frequencies', nFmw, did_fmw)
    ierr = nf90_def_dim(ncid,'n_MW_Radii',       nRmw, did_rmw)
    ierr = nf90_def_dim(ncid,'n_IR_Frequencies', nFir, did_fir)
    ierr = nf90_def_dim(ncid,'n_IR_Radii',       nRir, did_rir)
    ierr = nf90_def_dim(ncid,'n_Densities',      nDen_mw, did_d)
    ierr = nf90_def_dim(ncid,'n_IR_Densities',   nDen_ir, did_dir)
    ierr = nf90_def_dim(ncid,'n_Legendre_Terms', nLeg, did_l)
    ierr = nf90_def_dim(ncid,'n_Phase_Elements', nPhase, did_p)

    ierr = nf90_def_var(ncid,'Temperature', NF90_DOUBLE, (/did_t/), vid_T)
    ierr = nf90_put_att(ncid, vid_T, 'description', 'Temperature LUT dimension vector')
    ierr = nf90_put_att(ncid, vid_T, 'long_name',   'Temperature')
    ierr = nf90_put_att(ncid, vid_T, 'units',       'Kelvin (K)')
    ierr = nf90_put_att(ncid, vid_T, '_FillValue',  0.0_rk)

    ierr = nf90_def_var(ncid,'Frequency_MW', NF90_DOUBLE, (/did_fmw/), vid_fmw)
    ierr = nf90_put_att(ncid, vid_fmw, '_FillValue',  0.0_rk)
    ierr = nf90_put_att(ncid, vid_fmw, 'long_name',   'Frequency')
    ierr = nf90_put_att(ncid, vid_fmw, 'description','Microwave frequency LUT dimension vector')
    ierr = nf90_put_att(ncid, vid_fmw, 'units',      'GigaHertz (GHz)')

    ierr = nf90_def_var(ncid,'Frequency_IR', NF90_DOUBLE, (/did_fir/), vid_fir)
    ierr = nf90_put_att(ncid, vid_fir, '_FillValue',  0.0_rk)
    ierr = nf90_put_att(ncid, vid_fir, 'long_name',   'Frequency')
    ierr = nf90_put_att(ncid, vid_fir, 'description','Infrared frequency LUT dimension vector')
    ierr = nf90_put_att(ncid, vid_fir, 'units',      'Inverse centimetres (cm^-1)')

    ierr = nf90_def_var(ncid,'Reff_MW', NF90_DOUBLE, (/did_rmw/), vid_rmw)
    ierr = nf90_put_att(ncid, vid_rmw, '_FillValue',  0.0_rk)
    ierr = nf90_put_att(ncid, vid_rmw, 'long_name',   'Effective radius')
    ierr = nf90_put_att(ncid, vid_rmw, 'description','Microwave effective radius LUT dimension vector')
    ierr = nf90_put_att(ncid, vid_rmw, 'units',      'Microns (um)')

    ierr = nf90_def_var(ncid,'Reff_IR', NF90_DOUBLE, (/did_rir/), vid_rir)
    ierr = nf90_put_att(ncid, vid_rir, '_FillValue',  0.0_rk)
    ierr = nf90_put_att(ncid, vid_rir, 'long_name',   'Effective radius')
    ierr = nf90_put_att(ncid, vid_rir, 'description','Infrared effective radius LUT dimension vector')
    ierr = nf90_put_att(ncid, vid_rir, 'units',      'Microns (um)')

    ierr = nf90_def_var(ncid,'Density', NF90_DOUBLE, (/did_d/), vid_D)
    ierr = nf90_put_att(ncid, vid_D, '_FillValue',  0.0_rk)
    ierr = nf90_put_att(ncid, vid_D, 'long_name',   'Density')
    ierr = nf90_put_att(ncid, vid_D, 'description','Density LUT dimension vector')
    ierr = nf90_put_att(ncid, vid_D, 'units',      'Kilograms per cubic metre (kg.m^-3)')

    ! --------- Reversed orders for data variables ----------
    ! MW liquid: (F,R,T)
    ierr = nf90_def_var(ncid,'ke_L_MW', NF90_DOUBLE, (/did_fmw,did_rmw,did_t/), vid_keL_mw)
    ierr = nf90_put_att(ncid, vid_keL_mw, '_FillValue', 0.0_rk)
    ierr = nf90_put_att(ncid, vid_keL_mw, 'long_name',  'Microwave ke(L)')
    ierr = nf90_put_att(ncid, vid_keL_mw, 'description','Mass extinction coefficient for liquid phase microwave scatterers')
    ierr = nf90_put_att(ncid, vid_keL_mw, 'units',      'Metres squared per kilogram (m^2.kg^-1)')

    ierr = nf90_def_var(ncid,'w_L_MW', NF90_DOUBLE, (/did_fmw,did_rmw,did_t/), vid_wL_mw)
    ierr = nf90_put_att(ncid, vid_wL_mw, '_FillValue', 0.0_rk)
    ierr = nf90_put_att(ncid, vid_wL_mw, 'long_name',  'Microwave w(L)')
    ierr = nf90_put_att(ncid, vid_wL_mw, 'description','Single scatter albedo for liquid phase microwave scatterers')
    ierr = nf90_put_att(ncid, vid_wL_mw, 'units',      'N/A')

    ierr = nf90_def_var(ncid,'g_L_MW', NF90_DOUBLE, (/did_fmw,did_rmw,did_t/), vid_gL_mw)
    ierr = nf90_put_att(ncid, vid_gL_mw, '_FillValue', 0.0_rk)
    ierr = nf90_put_att(ncid, vid_gL_mw, 'long_name',  'Microwave g(L)')
    ierr = nf90_put_att(ncid, vid_gL_mw, 'description','Asymmetry parameter for liquid phase microwave scatterers')
    ierr = nf90_put_att(ncid, vid_gL_mw, 'units',      'N/A')

    ! MW liquid pcoeff: (F,R,T,L,P)
    ierr = nf90_def_var(ncid,'pcoeff_L_MW', NF90_DOUBLE, (/did_fmw,did_rmw,did_t,did_l,did_p/), vid_pL_mw)
    ierr = nf90_put_att(ncid, vid_pL_mw, '_FillValue', 0.0_rk)
    ierr = nf90_put_att(ncid, vid_pL_mw, 'long_name',  'Microwave pcoeff(L)')
    ierr = nf90_put_att(ncid, vid_pL_mw, 'description','Phase coefficients for liquid phase microwave scatterers')
    ierr = nf90_put_att(ncid, vid_pL_mw, 'units',      'N/A')

    ! MW solid: (F,R,D) and pcoeff (F,R,D,L,P)
    ierr = nf90_def_var(ncid,'ke_S_MW', NF90_DOUBLE, (/did_fmw,did_rmw,did_d/), vid_keS_mw)
    ierr = nf90_put_att(ncid, vid_keS_mw, '_FillValue', 0.0_rk)
    ierr = nf90_put_att(ncid, vid_keS_mw, 'long_name',  'Microwave ke(S)')
    ierr = nf90_put_att(ncid, vid_keS_mw, 'description','Mass extinction coefficient for solid phase microwave scatterers')
    ierr = nf90_put_att(ncid, vid_keS_mw, 'units',      'Metres squared per kilogram (m^2.kg^-1)')

    ierr = nf90_def_var(ncid,'w_S_MW', NF90_DOUBLE, (/did_fmw,did_rmw,did_d/), vid_wS_mw)
    ierr = nf90_put_att(ncid, vid_wS_mw, '_FillValue', 0.0_rk)
    ierr = nf90_put_att(ncid, vid_wS_mw, 'long_name',  'Microwave w(S)')
    ierr = nf90_put_att(ncid, vid_wS_mw, 'description','Single scatter albedo for solid phase microwave scatterers')
    ierr = nf90_put_att(ncid, vid_wS_mw, 'units',      'N/A')

    ierr = nf90_def_var(ncid,'g_S_MW', NF90_DOUBLE, (/did_fmw,did_rmw,did_d/), vid_gS_mw)
    ierr = nf90_put_att(ncid, vid_gS_mw, '_FillValue', 0.0_rk)
    ierr = nf90_put_att(ncid, vid_gS_mw, 'long_name',  'Microwave g(S)')
    ierr = nf90_put_att(ncid, vid_gS_mw, 'description','Asymmetry parameter for solid phase microwave scatterers')
    ierr = nf90_put_att(ncid, vid_gS_mw, 'units',      'N/A')

    ierr = nf90_def_var(ncid,'pcoeff_S_MW', NF90_DOUBLE, (/did_fmw,did_rmw,did_d,did_l,did_p/), vid_pS_mw)
    ierr = nf90_put_att(ncid, vid_pS_mw, '_FillValue', 0.0_rk)
    ierr = nf90_put_att(ncid, vid_pS_mw, 'long_name',  'Microwave pcoeff(S)')
    ierr = nf90_put_att(ncid, vid_pS_mw, 'description','Phase coefficients for solid phase microwave scatterers')
    ierr = nf90_put_att(ncid, vid_pS_mw, 'units',      'N/A')

    ! IR: ke/w/g (F,R,D), pcoeff (F,R,D,L)
    ierr = nf90_def_var(ncid,'ke_IR', NF90_DOUBLE, (/did_fir,did_rir,did_dir/), vid_ke_ir)
    ierr = nf90_put_att(ncid, vid_ke_ir, '_FillValue', 0.0_rk)
    ierr = nf90_put_att(ncid, vid_ke_ir, 'long_name',  'Infrared ke')
    ierr = nf90_put_att(ncid, vid_ke_ir, 'description','Mass extinction coefficient for infrared scatterers')
    ierr = nf90_put_att(ncid, vid_ke_ir, 'units',      'Metres squared per kilogram (m^2.kg^-1)')

    ierr = nf90_def_var(ncid,'w_IR', NF90_DOUBLE, (/did_fir,did_rir,did_dir/), vid_w_ir)
    ierr = nf90_put_att(ncid, vid_w_ir, '_FillValue', 0.0_rk)
    ierr = nf90_put_att(ncid, vid_w_ir, 'long_name',  'Infrared w')
    ierr = nf90_put_att(ncid, vid_w_ir, 'description','Single scatter albedo for infrared scatterers')
    ierr = nf90_put_att(ncid, vid_w_ir, 'units',      'N/A')

    ierr = nf90_def_var(ncid,'g_IR', NF90_DOUBLE, (/did_fir,did_rir,did_dir/), vid_g_ir)
    ierr = nf90_put_att(ncid, vid_g_ir, '_FillValue', 0.0_rk)
    ierr = nf90_put_att(ncid, vid_g_ir, 'long_name',  'Infrared g')
    ierr = nf90_put_att(ncid, vid_g_ir, 'description','Asymmetry parameter for infrared scatterers')
    ierr = nf90_put_att(ncid, vid_g_ir, 'units',      'N/A')

    ierr = nf90_def_var(ncid,'pcoeff_IR', NF90_DOUBLE, (/did_fir,did_rir,did_dir,did_l/), vid_p_ir)
    ierr = nf90_put_att(ncid, vid_p_ir, '_FillValue', 0.0_rk)
    ierr = nf90_put_att(ncid, vid_p_ir, 'long_name',  'Infrared pcoeff')
    ierr = nf90_put_att(ncid, vid_p_ir, 'description','Phase coefficients for infrared scatterers')
    ierr = nf90_put_att(ncid, vid_p_ir, 'units',      'N/A')

    ierr = nf90_enddef(ncid); if (ierr/=nf90_noerr) stop 'enddef failed'

    ierr = nf90_put_var(ncid, vid_T,    temp_k)
    ierr = nf90_put_var(ncid, vid_fmw,  freq_mw)
    ierr = nf90_put_var(ncid, vid_fir,  freq_ir)
    ierr = nf90_put_var(ncid, vid_rmw,  reff_mw)
    ierr = nf90_put_var(ncid, vid_rir,  reff_ir)
    ierr = nf90_put_var(ncid, vid_D,    density)

    ! ---- MW liquid zeros in (F,R,T) and (F,R,T,L,P)
    allocate(zero_f_r_t(nFmw,nRmw,nT)); zero_f_r_t = 0.0_rk
    allocate(zero_pL_f_r_t_l_p(nFmw,nRmw,nT,nLeg,1)); zero_pL_f_r_t_l_p = 0.0_rk
    ierr = nf90_put_var(ncid, vid_keL_mw, zero_f_r_t)
    ierr = nf90_put_var(ncid, vid_wL_mw,  zero_f_r_t)
    ierr = nf90_put_var(ncid, vid_gL_mw,  zero_f_r_t)
    ierr = nf90_put_var(ncid, vid_pL_mw,  zero_pL_f_r_t_l_p)

    ! ---- MW solid: map (R,F) -> (F,R,D=1)
    allocate(keS_f_r_d(nFmw,nRmw,1), wS_f_r_d(nFmw,nRmw,1), gS_f_r_d(nFmw,nRmw,1))
    do f=1,nFmw; do r=1,nRmw
      keS_f_r_d(f,r,1) = ke_S_mw(r,f)
      wS_f_r_d (f,r,1) = w_S_mw (r,f)
      gS_f_r_d (f,r,1) = g_S_mw (r,f)
    end do; end do
    ierr = nf90_put_var(ncid, vid_keS_mw, keS_f_r_d)
    ierr = nf90_put_var(ncid, vid_wS_mw,  wS_f_r_d)
    ierr = nf90_put_var(ncid, vid_gS_mw,  gS_f_r_d)

    ! ---- MW solid pcoeff: (P,L,T,R,F) -> (F,R,D=1,L,P)
    allocate(pS_f_r_d_l_p(nFmw,nRmw,1,nLeg,1)); pS_f_r_d_l_p = 0.0_rk
    do f=1,nFmw
      do r=1,nRmw
        do l=1,nLeg
          ! pS_mw_PLTRF(P=1,L=l,T=1..nT,R=r,F=f) reduces across T? Keep T=1 (as in your MW P11 at T190)
          pS_f_r_d_l_p(f,r,1,l,1) = pS_mw_PLTRF(1,l,1,r,f)
        end do
      end do
    end do
    ierr = nf90_put_var(ncid, vid_pS_mw,  pS_f_r_d_l_p)

    ! ---- IR solid: map (R,F) -> (F,R,D=2); zeros for D=1 (liquid)
    allocate(keIR_f_r_d(nFir,nRir,nDen_ir), wIR_f_r_d(nFir,nRir,nDen_ir), gIR_f_r_d(nFir,nRir,nDen_ir))
    keIR_f_r_d = 0.0_rk; wIR_f_r_d = 0.0_rk; gIR_f_r_d = 0.0_rk
    do f=1,nFir; do r=1,nRir
      keIR_f_r_d(f,r,2) = ke_ir_solid_RxF(r,f)
      wIR_f_r_d (f,r,2) = w_ir_solid_RxF (r,f)
      gIR_f_r_d (f,r,2) = g_ir_solid_RxF (r,f)
    end do; end do
    ierr = nf90_put_var(ncid, vid_ke_ir, keIR_f_r_d)
    ierr = nf90_put_var(ncid, vid_w_ir,  wIR_f_r_d)
    ierr = nf90_put_var(ncid, vid_g_ir,  gIR_f_r_d)

    ! ---- IR pcoeff: (L,R,F) -> (F,R,D=2,L)
    allocate(pIR_f_r_d_l(nFir,nRir,nDen_ir,nLeg)); pIR_f_r_d_l = 0.0_rk
    do f=1,nFir; do r=1,nRir; do l=1,nLeg
      pIR_f_r_d_l(f,r,2,l) = pcoeff_ir_solid_LxRxF(l,r,f)
    end do; end do; end do
    ierr = nf90_put_var(ncid, vid_p_ir, pIR_f_r_d_l)

    call put_text_att(ncid, 'Title',   'CRTM Cloud Optical Properties in Infrared/Visible and Microwave Ranges')
    call put_text_att(ncid, 'History', 'history text')
    call put_text_att(ncid, 'Comment', 'Comment text')
    call put_int_att(ncid,  'Release', 3)
    call put_int_att(ncid,  'Version', 4)

    ierr = nf90_close(ncid); if (ierr/=nf90_noerr) stop 'close failed'
  end subroutine write_cloudcoeff_mw_ir
end module mw_ir_netcdf_writer

program CloudCoeff_Generate_MW_Schemes
  use util_kinds
  use scheme_mass
  use text_io
  use legendre_tools
  use mw_ir_netcdf_writer, only: write_cloudcoeff_mw_ir
  implicit none

  integer, parameter :: nFmw = 70, nFir = 470, nR = 131, nT = 5, nLeg = 39, nPhase = 1
  character(len=*), parameter :: base_mw  = 'SnowflakeModel_MW'
  character(len=*), parameter :: base_ir  = 'SnowflakeModel_UVIR'
  character(len=*), parameter :: temps(nT) = (/'T190','T210','T230','T250','T270'/)
  real(rk) :: freq_mw(nFmw), freq_ir(nFir)
  real(rk) :: temp_k(nT)
  integer :: i
  character(len=3) :: tstr
  integer :: tval

  do i=1,nT
    tstr = temps(i)(2:4)
    read(tstr,'(I3)') tval
    temp_k(i) = real(tval, rk)
  end do

  call read_freq_mw(freq_mw)
  call do_scheme_combined(1, freq_mw, temp_k)
  call do_scheme_combined(2, freq_mw, temp_k)
  call do_scheme_combined(3, freq_mw, temp_k)
  call do_scheme_combined(4, freq_mw, temp_k)

contains

  subroutine read_freq_mw(freq_mw)
    real(rk), intent(out) :: freq_mw(:)
    integer :: nread
    call read_vector_file('MW_freq.txt', freq_mw, nread, .true.)
    if (nread /= size(freq_mw)) stop 'MW_freq.txt length mismatch'
  end subroutine read_freq_mw

  subroutine do_scheme_combined(sid, freq_mw, temp_k)
    integer,  intent(in) :: sid
    real(rk), intent(in) :: freq_mw(:), temp_k(:)

    character(len=512) :: ddir_mw, ddir_ir, tdir, iscapath_mw, iscapath_ir, p11path_ir, outfile, p11path_mw
    real(rk) :: fchk_mw(nFmw), dmax_um_mw(nR)
    real(rk) :: fchk_ir(nFir), dmax_um_ir(nR)
    real(rk) :: masses_g(nR), reff_um(nR), reff_ir(nR)
    real(rk) :: ke_sf(nR,nFmw), w_sf(nR,nFmw), g_sf(nR,nFmw)
    real(rk) :: ke_ir_s(nR,nFir), w_ir_s(nR,nFir), g_ir_s(nR,nFir)
    real(rk), allocatable :: pS_mw_PLTRF(:,:,:,:,:)         ! (P,L,T,R,F)
    real(rk), allocatable :: pcoeff_ir_leg_r_f(:,:,:)       ! (L,R,F)
    integer :: t, s

    write(ddir_mw,'(A,"/Density",I1)') trim(base_mw), sid
    write(ddir_ir,'(A,"/Density",I1)') trim(base_ir), sid

    write(tdir,'(A,"/",A)') trim(ddir_mw), 'T190'
    write(iscapath_mw,'(A,"/isca_freq.dat")') trim(tdir)
    call read_isca_block_header_local(iscapath_mw, nFmw, nR, fchk_mw, dmax_um_mw)
    call assert_match(freq_mw, fchk_mw, 1.0e-3_rk)

    write(iscapath_ir,'(A,"/isca_wn.dat")') trim(ddir_ir)
    call read_isca_wn_header_local(iscapath_ir, nFir, nR, fchk_ir, dmax_um_ir)

    do s=1,nR
      select case(sid)
      case(1); masses_g(s) = mass_g_density1(dmax_um_mw(s))
      case(2); masses_g(s) = mass_g_density2(dmax_um_mw(s))
      case(3); masses_g(s) = mass_g_density3(dmax_um_mw(s))
      case(4); masses_g(s) = mass_g_density4(dmax_um_mw(s))
      end select
      reff_um(s) = reff_um_from_mass_g(masses_g(s))
    end do
    reff_ir = reff_um
    call assert_monotonic_increasing(reff_um)

    ke_sf = 0.0_rk; w_sf = 0.0_rk; g_sf = 0.0_rk
    do t=1,nT
      write(tdir,'(A,"/",A)') trim(ddir_mw), temps(t)
      write(iscapath_mw,'(A,"/isca_freq.dat")') trim(tdir)
      call stream_isca_fill_local(iscapath_mw, nFmw, nR, freq_mw, masses_g, ke_sf, w_sf, g_sf)
    end do

    call stream_isca_wn_fill_local(iscapath_ir, nFir, nR, fchk_ir, masses_g, ke_ir_s, w_ir_s, g_ir_s)

    allocate(pS_mw_PLTRF(1,nLeg,nT,nR,nFmw)); pS_mw_PLTRF = 0.0_rk
    write(p11path_mw,'(A,"/T190/P11.dat")') trim(ddir_mw)
    call compute_pcoeff_mw_from_P11(p11path_mw, pS_mw_PLTRF)

    allocate(pcoeff_ir_leg_r_f(nLeg,nR,nFir)); pcoeff_ir_leg_r_f = 0.0_rk
    write(p11path_ir,'(A,"/P11.dat")') trim(ddir_ir)
    call compute_pcoeff_ir_from_P11_local(p11path_ir, nLeg, nR, nFir, pcoeff_ir_leg_r_f)

    write(outfile,'("CloudCoeff_MW_Snowflake_Density",I1,".nc4")') sid
    call write_cloudcoeff_mw_ir(outfile, &
                                freq_mw, reff_um, temp_k, density_value_for(sid), &
                                ke_sf, w_sf, g_sf, pS_mw_PLTRF, &
                                fchk_ir, reff_ir, &
                                ke_ir_s, w_ir_s, g_ir_s, pcoeff_ir_leg_r_f)
  end subroutine do_scheme_combined

  subroutine compute_pcoeff_mw_from_P11(p11file, pS_PLTRF)
    integer, parameter :: nFmw_loc = nFmw, nR_loc = nR, nLeg_loc = nLeg
    character(*), intent(in)  :: p11file
    real(rk),     intent(out) :: pS_PLTRF(1,nLeg_loc,nT,nR_loc,nFmw_loc)   ! (P,L,T,R,F)
    integer, parameter :: nAng = 498
    real(rk) :: theta_deg(nAng), theta(nAng), mu(nAng), w(nAng)
    real(rk) :: P11(nAng), coef1(0:nLeg_loc-1)
    integer :: u, ios, f, r
    open(newunit=u, file=trim(p11file), status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'open MW P11.dat failed'
    read(u,*,iostat=ios) theta_deg
    if (ios /= 0) stop 'failed reading MW P11 angle header'
    theta = theta_deg*(acos(-1.0_rk)/180.0_rk)
    mu    = cos(theta)
    call trapezoid_mu_weights(theta, w)
    do f=1,nFmw_loc
      do r=1,nR_loc
        read(u,*,iostat=ios) P11
        if (ios /= 0) stop 'bad MW P11 line'
        call p11_to_coef(mu, w, P11, nLeg_loc, coef1)
        call pack_stream_windows_mw(nLeg_loc, coef1, pS_PLTRF(1,1,1,r,f))
      end do
    end do
    close(u)
  end subroutine compute_pcoeff_mw_from_P11

  subroutine pack_stream_windows_mw(nLeg, coef1, pvec)
    integer,  intent(in) :: nLeg
    real(rk), intent(in) :: coef1(0:nLeg-1)
    real(rk), intent(out):: pvec(nLeg)     ! L = 1..nLeg
    integer, parameter :: DEF_N_STREAM_SETS = 5
    integer, parameter :: LOWS(DEF_N_STREAM_SETS)  = (/0, 0, 5, 12, 21/)
    integer, parameter :: HIGHS(DEF_N_STREAM_SETS) = (/2, 4,11, 20, 37/)
    integer :: i, low, high
    pvec = 0.0_rk
    do i=1,DEF_N_STREAM_SETS
      low  = LOWS(i); high = HIGHS(i)
      if (high >= low) pvec(low+1:high+1) = 0.5_rk*coef1(low:high)
    end do
  end subroutine pack_stream_windows_mw

  real(rk) function density_value_for(sid) result(rho)
    integer, intent(in) :: sid
    select case(sid)
    case(1)
      rho = 100.0_rk
    case default
      rho = 0.0_rk
    end select
  end function density_value_for

  subroutine assert_match(a, b, tol)
    real(rk), intent(in) :: a(:), b(:), tol
    integer :: i
    if (size(a) /= size(b)) stop 'vector length mismatch'
    do i=1,size(a)
      if (abs(a(i)-b(i)) > tol) stop 'value mismatch'
    end do
  end subroutine assert_match

  subroutine assert_monotonic_increasing(x)
    real(rk), intent(in) :: x(:)
    integer :: i
    do i=1,size(x)-1
      if (x(i+1) <= x(i)) stop 'Reff not strictly increasing'
    end do
  end subroutine assert_monotonic_increasing

  subroutine read_isca_block_header_local(path, nfreq, nsize, freq_out, dmax_um)
    character(*), intent(in)  :: path
    integer,      intent(in)  :: nfreq, nsize
    real(rk),     intent(out) :: freq_out(nfreq)
    real(rk),     intent(out) :: dmax_um(nsize)
    integer :: u, ios, f, s
    real(rk) :: fv, du, vv, av, qv, wv, gv
    open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'failed to open isca_freq.dat'
    do f=1,nfreq
      do s=1,nsize
        read(u,*,iostat=ios) fv, du, vv, av, qv, wv, gv
        if (ios /= 0) stop 'bad isca_freq.dat format'
        if (s == 1) freq_out(f) = fv
        if (f == 1) dmax_um(s)  = du
      end do
    end do
    close(u)
  end subroutine read_isca_block_header_local

  subroutine stream_isca_fill_local(path, nfreq, nsize, freq_ref, masses_g, ke_sf, w_sf, g_sf)
    character(*), intent(in)  :: path
    integer,      intent(in)  :: nfreq, nsize
    real(rk),     intent(in)  :: freq_ref(nfreq)
    real(rk),     intent(in)  :: masses_g(nsize)
    real(rk),     intent(out) :: ke_sf(nsize,nfreq)
    real(rk),     intent(out) :: w_sf(nsize,nfreq)
    real(rk),     intent(out) :: g_sf(nsize,nfreq)
    integer :: u, ios, f, s
    real(rk) :: fv, du, vv, av, qv, wv, gv
    real(rk) :: area_m2, m_kg
    open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'open isca failed'
    do f=1,nfreq
      do s=1,nsize
        read(u,*,iostat=ios) fv, du, vv, av, qv, wv, gv
        if (ios /= 0) stop 'bad isca record'
        if (abs(fv - freq_ref(f)) > 1.0e-3_rk) stop 'frequency mismatch in isca_freq.dat'
        area_m2     = av*1.0e-12_rk
        m_kg        = masses_g(s)*1.0e-3_rk
        ke_sf(s,f)  = qv*area_m2 / max(m_kg, 1.0e-300_rk)
        w_sf(s,f)   = min(max(wv, 0.0_rk), 1.0_rk)
        g_sf(s,f)   = max(min(gv, 1.0_rk), -1.0_rk)
      end do
    end do
    close(u)
  end subroutine stream_isca_fill_local

  subroutine read_isca_wn_header_local(path, nwn, nsize, wn_out, dmax_um)
    character(*), intent(in)  :: path
    integer,      intent(in)  :: nwn, nsize
    real(rk),     intent(out) :: wn_out(nwn)
    real(rk),     intent(out) :: dmax_um(nsize)
    integer :: u, ios, f, s
    real(rk) :: wavnum, du, vv, av, qv, wv, gv
    open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'failed to open isca_wn.dat'
    do f=1,nwn
      do s=1,nsize
        read(u,*,iostat=ios) wavnum, du, vv, av, qv, wv, gv
        if (ios /= 0) stop 'bad isca_wn.dat format'
        if (s == 1) wn_out(f) = wavnum
        if (f == 1) dmax_um(s) = du
      end do
    end do
    close(u)
  end subroutine read_isca_wn_header_local

  subroutine stream_isca_wn_fill_local(path, nwn, nsize, wn_ref, masses_g, ke_ir, w_ir, g_ir)
    character(*), intent(in)  :: path
    integer,      intent(in)  :: nwn, nsize
    real(rk),     intent(in)  :: wn_ref(nwn)
    real(rk),     intent(in)  :: masses_g(nsize)
    real(rk),     intent(out) :: ke_ir(nsize,nwn)
    real(rk),     intent(out) :: w_ir(nsize,nwn)
    real(rk),     intent(out) :: g_ir(nsize,nwn)
    integer :: u, ios, f, s
    real(rk) :: wavnum, du, vv, av, qv, wv, gv
    real(rk) :: area_m2, m_kg
    open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'open isca_wn failed'
    do f=1,nwn
      do s=1,nsize
        read(u,*,iostat=ios) wavnum, du, vv, av, qv, wv, gv
        if (ios /= 0) stop 'bad isca_wn record'
        if (abs(wavnum - wn_ref(f)) > 1.0e-8_rk) stop 'wavenumber mismatch in isca_wn.dat'
        area_m2     = av*1.0e-12_rk
        m_kg        = masses_g(s)*1.0e-3_rk
        ke_ir(s,f)  = qv*area_m2 / max(m_kg, 1.0e-300_rk)
        w_ir(s,f)   = min(max(wv, 0.0_rk), 1.0_rk)
        g_ir(s,f)   = max(min(gv, 1.0_rk), -1.0_rk)
      end do
    end do
    close(u)
  end subroutine stream_isca_wn_fill_local

  subroutine compute_pcoeff_ir_from_P11_local(p11file, nLeg, nR, nF, pcoeff_leg_rf)
    character(*), intent(in) :: p11file
    integer,      intent(in) :: nLeg, nR, nF
    real(rk),     intent(out):: pcoeff_leg_rf(nLeg,nR,nF)
    integer, parameter :: nAng = 498
    real(rk) :: theta_deg(nAng), theta(nAng), mu(nAng), w(nAng)
    real(rk) :: P11(nAng), coef1(0:nLeg-1)
    integer :: u, ios, f, r
    open(newunit=u, file=trim(p11file), status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'open IR P11.dat failed'
    read(u,*,iostat=ios) theta_deg
    if (ios /= 0) stop 'failed reading IR P11 angle header'
    theta = theta_deg*(acos(-1.0_rk)/180.0_rk)
    mu    = cos(theta)
    call trapezoid_mu_weights(theta, w)
    do f=1,nF
      do r=1,nR
        read(u,*,iostat=ios) P11
        if (ios /= 0) stop 'bad IR P11 line'
        call p11_to_coef(mu, w, P11, nLeg, coef1)
        call pack_stream_windows_ir_local(nLeg, coef1, pcoeff_leg_rf(1,r,f))
      end do
    end do
    close(u)
  end subroutine compute_pcoeff_ir_from_P11_local

  subroutine pack_stream_windows_ir_local(nLeg, coef1, pvec)
    integer,  intent(in) :: nLeg
    real(rk), intent(in) :: coef1(0:nLeg-1)
    real(rk), intent(out):: pvec(nLeg)
    integer, parameter :: DEF_N_STREAM_SETS = 5
    integer, parameter :: LOWS(DEF_N_STREAM_SETS)  = (/0, 0, 5, 12, 21/)
    integer, parameter :: HIGHS(DEF_N_STREAM_SETS) = (/2, 4,11, 20, 37/)
    integer :: i, low, high
    pvec = 0.0_rk
    do i=1,DEF_N_STREAM_SETS
      low  = LOWS(i); high = HIGHS(i)
      if (high >= low) pvec(low+1:high+1) = 0.5_rk*coef1(low:high)
    end do
  end subroutine pack_stream_windows_ir_local

end program CloudCoeff_Generate_MW_Schemes
