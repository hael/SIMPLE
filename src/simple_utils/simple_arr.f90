!------------------------------------------------------------------------------!
! SIMPLE v3.0 Elmlund & Elmlund Lab simplecryoem.com !
!------------------------------------------------------------------------------!
!> Simple array class
!
! The code is distributed with the hope that it will be useful, but WITHOUT ANY
! WARRANTY. Redistribution and modification is regulated by the GNU General
! Public License.
!-----------------------------------------------------------------------------!
module simple_arr
implicit none

public :: arr
private
!> General array type
type :: arr
    private
    integer, allocatable :: iarr(:)
    real, allocatable    :: rarr(:)
  contains
    procedure, private :: new_1
    procedure, private :: new_2
    generic   :: assignment(=) => new_1, new_2
    procedure :: iget
    procedure :: rget
    procedure :: display
    procedure :: kill
end type arr

contains

    subroutine new_1( self, iarr )
        class(arr), intent(inout) :: self
        integer, intent(in) :: iarr(:)        !< input integer array
        if( allocated(self%iarr) ) deallocate( self%iarr )
        allocate( self%iarr(size(iarr)), source=iarr )
    end subroutine

    subroutine new_2( self, rarr )
        class(arr), intent(inout) :: self
        real, intent(in) :: rarr(:)           !< input float array
        if( allocated(self%rarr) ) deallocate( self%rarr )
        allocate( self%rarr(size(rarr)), source=rarr )
    end subroutine

    function iget( self ) result( iarr )
        class(arr), intent(in) :: self
        integer, allocatable   :: iarr(:)  !< output integer array
        if( allocated(self%iarr) )then
            allocate( iarr(size(self%iarr)), source=self%iarr )
        else
            stop 'no info in iarr; get_1; simple_arr'
        endif
    end function

    function rget( self ) result( rarr )
        class(arr), intent(in) :: self
        real, allocatable      :: rarr(:) !< output float array
        if( allocated(self%rarr) )then
            allocate( rarr(size(self%rarr)), source=self%rarr )
        else
            stop 'no info in rarr; get_2; simple_arr'
        endif
    end function

    subroutine display( self )
        class(arr), intent(in) :: self
        if( allocated(self%iarr) ) print *, self%iarr
        if( allocated(self%rarr) ) print *, self%rarr
    end subroutine display

    subroutine kill( self )
        class(arr), intent(inout) :: self
        if( allocated(self%iarr) ) deallocate( self%iarr )
        if( allocated(self%rarr) ) deallocate( self%rarr )
    end subroutine

end module simple_arr
