!==============================================================
! freq_io.F90 (updated)
! Adds readers for MW and UVIR effective radius lists (microns)
! Files expected:
!   MW_freq.txt      -> Frequency_MW (GHz)
!   UVIR_WN.txt      -> Frequency_IR (cm^-1)
!   Reff_MW.txt      -> Reff_MW (microns)
!   Reff_UVIR.txt    -> Reff_IR (microns)
!==============================================================

module freq_io
  use iso_fortran_env, only: real64
  implicit none
contains

  subroutine read_frequencies_mw(filename, Frequency_MW, n_expected)
    character(len=*), intent(in)  :: filename
    real(real64),     intent(out) :: Frequency_MW(:)
    integer,          intent(in),  optional :: n_expected
    call read_vector_file(filename, Frequency_MW, "GHz", n_expected)
  end subroutine read_frequencies_mw

  subroutine read_wavenumbers_ir(filename, Frequency_IR, n_expected)
    character(len=*), intent(in)  :: filename
    real(real64),     intent(out) :: Frequency_IR(:)
    integer,          intent(in),  optional :: n_expected
    call read_vector_file(filename, Frequency_IR, "cm^-1", n_expected)
  end subroutine read_wavenumbers_ir

  subroutine read_reff_mw(filename, Reff_MW, n_expected)
    character(len=*), intent(in)  :: filename
    real(real64),     intent(out) :: Reff_MW(:)
    integer,          intent(in),  optional :: n_expected
    call read_vector_file(filename, Reff_MW, "micron", n_expected)
  end subroutine read_reff_mw

  subroutine read_reff_uvir(filename, Reff_IR, n_expected)
    character(len=*), intent(in)  :: filename
    real(real64),     intent(out) :: Reff_IR(:)
    integer,          intent(in),  optional :: n_expected
    call read_vector_file(filename, Reff_IR, "micron", n_expected)
  end subroutine read_reff_uvir

  subroutine read_vector_file(filename, vec, units_label, n_expected)
    character(len=*), intent(in)  :: filename
    real(real64),     intent(out) :: vec(:)
    character(len=*), intent(in)  :: units_label
    integer,          intent(in),  optional :: n_expected

    integer :: u, ios, i, count_read
    character(len=4096) :: line
    real(real64) :: value

    count_read = 0
    open(newunit=u, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) then
      write(*,*) "ERROR: cannot open ", trim(filename)
      stop 1
    end if

    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (is_blank_or_comment(line)) cycle
      value = read_first_real_token(line, ios)
      if (ios == 0) then
        count_read = count_read + 1
        if (count_read <= size(vec)) then
          vec(count_read) = value
        else
          write(*,*) "ERROR: more entries in ", trim(filename), " than size(vec)"
          close(u)
          stop 1
        end if
      end if
    end do
    close(u)

    if (present(n_expected)) then
      if (n_expected /= count_read) then
        write(*,*) "ERROR: count mismatch for ", trim(filename), " expected=", n_expected, " read=", count_read
        stop 1
      end if
    end if

    if (count_read < size(vec)) then
      do i = count_read+1, size(vec)
        vec(i) = 0.0_real64
      end do
    end if
  end subroutine read_vector_file

  logical function is_blank_or_comment(line)
    character(len=*), intent(in) :: line
    integer :: k, n
    character :: c
    n = len_trim(line)
    is_blank_or_comment = .true.
    if (n == 0) return
    do k = 1, n
      c = line(k:k)
      if (c /= ' ' .and. c /= char(9)) then
        if (c == '#' .or. c == '!') then
          is_blank_or_comment = .true.
        else
          is_blank_or_comment = .false.
        end if
        return
      end if
    end do
  end function is_blank_or_comment

  real(real64) function read_first_real_token(line, ios) result(val)
    character(len=*), intent(in) :: line
    integer,          intent(out) :: ios
    integer :: p1, p2, n
    character(len=:), allocatable :: token
    character :: c

    n = len_trim(line)
    p1 = 0
    do p2 = 1, n
      c = line(p2:p2)
      if (c /= ' ' .and. c /= char(9) .and. c /= ',' ) then
        p1 = p2
        exit
      end if
    end do
    if (p1 == 0) then
      ios = 1
      val = 0.0_real64
      return
    end if

    do p2 = p1, n
      c = line(p2:p2)
      if (c == ' ' .or. c == char(9) .or. c == ',' ) exit
    end do
    if (p2 <= n) then
      token = line(p1:p2-1)
    else
      token = line(p1:n)
    end if

    read(token, *, iostat=ios) val
    if (ios /= 0) val = 0.0_real64
  end function read_first_real_token

end module freq_io
