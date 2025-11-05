module simple_stack
!$ use omp_lib
include 'simple_lib.f08'
implicit none
private
public :: stack, pop, contains, is_empty, remove
#include "simple_local_flags.inc"

type :: stack
    private
    integer           :: max_size   !< maximum size of the stack, default to 100
    real, allocatable :: val_arr(:) !< contains the values of the stack
    integer           :: last_ind   !< index of the last element of the stack
contains
procedure :: new
procedure, private :: push_obj
procedure, private :: push_arr
generic   :: push => push_obj, push_arr
procedure :: pop
procedure :: get_at
procedure :: is_empty
procedure :: contains
procedure :: remove
procedure :: size_of
procedure :: clear
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

    subroutine clear(self)
        class(stack), intent(inout) :: self
        self%last_ind = 1
    end subroutine clear

    subroutine push_obj(self, obj)
        class(stack), intent(inout) :: self
        integer     , intent(in)    :: obj
        self%val_arr(self%last_ind) = obj
        self%last_ind               = self%last_ind + 1
    end subroutine push_obj

    subroutine push_arr(self, arr)
        class(stack), intent(inout) :: self
        integer     , intent(in)    :: arr(:)
        integer :: k
        do k = 1, size(arr)
            self%val_arr(self%last_ind) = arr(k)
            self%last_ind               = self%last_ind + 1
        enddo
    end subroutine push_arr

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

    subroutine remove(self, val)
        class(stack), intent(inout) :: self
        real    :: val
        integer :: k, l, N
        N = self%last_ind
        do k = 1, N
            if( self%val_arr(k) == val )then
                do l = k, N-1
                    self%val_arr(l) = self%val_arr(l+1)
                enddo
                self%last_ind = self%last_ind - 1
            endif
        enddo
    end subroutine remove
    
    function get_at(self, ind) result(ret)
        class(stack), intent(in) :: self
        integer     , intent(in) :: ind
        integer :: ret
        if( ind < self%last_ind )then
            ret = self%val_arr(ind)
        else
            ret = -1
        endif
    end function get_at

    function contains(self, obj_val) result(ret)
        class(stack), intent(in) :: self
        integer     , intent(in) :: obj_val
        logical :: ret
        integer :: k
        ret = .false.
        do k = 1, self%last_ind-1
            if( self%val_arr(k) == obj_val )then
                ret = .true.
            endif
        enddo
    end function contains

    function is_empty(self) result(ret)
        class(stack), intent(in) :: self
        logical :: ret
        ret = .false.
        if( self%last_ind == 1 ) ret = .true.
    end function is_empty
    
    function size_of(self) result(ret)
        class(stack), intent(in) :: self
        integer :: ret
        ret = self%last_ind - 1
    end function size_of
end module simple_stack