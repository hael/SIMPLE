module simple_stack
!$ use omp_lib
include 'simple_lib.f08'
implicit none
private
public :: stack
#include "simple_local_flags.inc"

type :: stack
    private
    integer           :: max_size   !< maximum size of the stack, default to 100
    real, allocatable :: val_arr(:) !< contains the values of the stack
    integer           :: last_ind   !< index of the last element of the stack
contains
procedure :: new
procedure :: push
procedure :: pop
procedure :: contains
! DESTRUCTOR
!procedure :: kill
end type stack

interface stack
    module procedure constructor
end interface stack

contains
    ! CONSTRUCTORS

    !>  \brief  is a constructor
    function constructor( max_sz, vals ) result( self ) !(FAILS W PRESENT GFORTRAN)
        integer, optional, intent(in) :: max_sz
        real   , optional, intent(in) :: vals(:)
        type(stack) :: self
        call self%new( max_sz, vals )
    end function constructor

    !>  \brief  Constructor for simple_image class
    subroutine new( self, max_sz, vals )
        class(stack),      intent(inout) :: self
        integer, optional, intent(in)    :: max_sz
        real   , optional, intent(in)    :: vals(:)
        integer :: k
        if( present(max_sz) )then
            self%max_size = max_sz
        else
            self%max_size = 100 
        endif
        allocate(self%val_arr(self%max_size), source=0.)
        if( present(vals) )then
            self%last_ind = size(vals) + 1
            do k = 1, size(vals)
                self%val_arr(k) = vals(k)
            enddo
        else
            self%last_ind = 1
        endif
    end subroutine new

    subroutine push(self, obj)
        class(stack), intent(inout) :: self
        integer     , intent(in)    :: obj
        self%val_arr(self%last_ind) = obj
        self%last_ind               = self%last_ind + 1
    end subroutine push

    function pop(self) result(ret)
        class(stack), intent(inout) :: self
        integer :: ret
        if( self%last_ind > 1 )then
            ret           = self%val_arr(self%last_ind - 1)
            self%last_ind = self%last_ind - 1
        else
            ret = -1
        endif
    end function pop

    function contains(self, obj_val) result(ret)
        class(stack), intent(inout) :: self
        integer     , intent(in)    :: obj_val
        logical :: ret
        integer :: k
        ret = .false.
        do k = 1, self%last_ind-1
            if( self%val_arr(k) == obj_val )then
                ret = .true.
            endif
        enddo
    end function contains
end module simple_stack