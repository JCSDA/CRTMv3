module csv_data_define
    
    implicit none
    
    private

    type, public :: csv_data
        integer :: columns
        integer :: rows
        logical :: is_allocated = .false.
        character(12) :: description
        real, allocatable :: data(:,:)
    contains
        procedure, public, pass :: csv_data_create
        final :: csv_data_destroy 
    end type csv_data
    
contains

  !**************************************
  !
  ! CONSTRUCTOR
  !
  !**************************************
  subroutine csv_data_create(this,c,r)
      integer, intent(in) :: c
      integer, intent(in) :: r
      class(csv_data), intent(inout) :: this
      if(.not. this%is_allocated) then
        allocate(this%data(r,c))
        this%is_allocated = .true.
      end if
  end subroutine csv_data_create

  !*********************************
  !
  ! FINALIZER
  !
  !*********************************
  subroutine csv_data_destroy(this)
      type(csv_data), intent(inout) :: this
      if(this%is_allocated) then
        deallocate(this%data)
      end if
      this%is_allocated = .false.
  end subroutine csv_data_destroy

end module csv_data_define
