module csv_data_io
    use csv_data_define, only: csv_data
    implicit none

    private
    public :: count_csv_columns
    public :: count_csv_rows
    public :: read_csv

contains
    
function count_csv_columns(filename) result(ncols)
    implicit none
    character(len=*), intent(in) :: filename
    integer :: ncols
    integer :: unit, ios, i, lenline
    character(len=4096) :: line
    CHARACTER(len=255) :: cwd

    ncols = -1   ! default: error indicator

    ! Open the file
    !inquire(file=filename, number=unit)
    !if (unit == 0) unit = 21

    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) "CSV file not found!"
        CALL getcwd(cwd)
        WRITE(*,*) TRIM(cwd), unit
        return
    end if

    ! Read until we find a non-empty line
    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) then
            close(unit)
            write(*,*) "Error reading CSV file lines!"
            return
        end if
        lenline = len_trim(line)
        if (lenline > 0) exit
    end do

    ! Count commas in the line
    ncols = 1
    do i = 1, len_trim(line)
        if (line(i:i) == ',') ncols = ncols + 1
    end do

    close(unit)
end function count_csv_columns

function count_csv_rows(filename) result(nrows)
    implicit none
    character(len=*), intent(in) :: filename
    integer :: nrows
    integer :: unit, ios
    character(len=4096) :: line

    nrows = -1   ! default (error)

    ! Get or assign a unit number
    !inquire(file=filename, number=unit)
    !if (unit == 0) unit = 21

    ! Open the file
    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) return

    nrows = 0

    ! Read line-by-line
    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit

        if (len_trim(line) > 0) then
            nrows = nrows + 1
        end if
    end do

    close(unit)
end function count_csv_rows

function read_csv(filename, csv) result(status)
        !! Reads numeric CSV file into csv%data (already allocated)
        !!
        !! Returns 0 on success, nonzero on error.
        implicit none

        character(len=*), intent(in) :: filename
        type(csv_data), intent(inout) :: csv
        integer :: status
        integer :: i, j, ios
        character(len=4096) :: line
        character(len=64) :: token
        integer :: pos, start, col

        status = 0

        open(unit=444, file=filename, status='old', action='read', &
             iostat=ios)
        if (ios /= 0) then
            status = 1
            return
        end if

        do i = 1, csv%rows
            read(444, '(A)', iostat=ios) line
            if (ios /= 0) then
                status = 2
                close(444)
                return
            end if

            ! Parse the line manually
            pos = 1
            col = 1

            do while (col <= csv%columns)
                call parse_next_token(line, pos, token)
                read(token, *, iostat=ios) csv%data(i, col)
                if (ios /= 0) then
                    status = 3
                    close(444)
                    return
                end if
                col = col + 1
            end do
        end do

        close(444)
    end function read_csv
    
    subroutine parse_next_token(line, pos, token)
        !! Extracts the next comma-separated token from `line`
        !! starting at character index `pos`.
        implicit none
        character(len=*), intent(in)    :: line
        integer, intent(inout)          :: pos
        character(len=*), intent(out)   :: token

        integer :: lenline, start, endpos

        lenline = len_trim(line)
        token = ""

        ! Skip leading spaces
        do while (pos <= lenline .and. line(pos:pos) == ' ')
            pos = pos + 1
        end do

        start = pos

        ! Find next comma or end of line
        do while (pos <= lenline .and. line(pos:pos) /= ',')
            pos = pos + 1
        end do

        endpos = pos - 1

        ! Copy token
        if (start <= endpos) then
            token = adjustl(line(start:endpos))
        else
            token = ""
        end if

        ! Skip comma
        if (pos <= lenline .and. line(pos:pos) == ',') pos = pos + 1
    end subroutine parse_next_token

end module csv_data_io
