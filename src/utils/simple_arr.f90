! array class (container class for the singly linked list)
module simple_arr
use simple_defs
use simple_error, only: allocchk, simple_exception
implicit none

public :: arr
private
#include "simple_local_flags.inc"

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
        if( allocated(self%iarr) )then
            deallocate( self%iarr,STAT=alloc_stat)
            if(alloc_stat /= 0) call allocchk(" In simple_arr::new_1  dealloc iarr ", alloc_stat )
        end if
        allocate( self%iarr(size(iarr)), source=iarr, stat=alloc_stat)
        if(alloc_stat /= 0) call allocchk(" In simple_arr::new_1  iarr " , alloc_stat)
    end subroutine new_1

    subroutine new_2( self, rarr )
        class(arr), intent(inout) :: self
        real, intent(in) :: rarr(:)           !< input float array
        if( allocated(self%rarr) ) then
            deallocate( self%rarr ,STAT=alloc_stat)
            if(alloc_stat /= 0) call allocchk(" In simple_arr::new_2 dealloc rarr" , alloc_stat)
          end if
        allocate( self%rarr(size(rarr)), source=rarr,STAT=alloc_stat )
        if(alloc_stat /= 0) call allocchk(" In simple_arr::new_2 rarr ", alloc_stat)
    end subroutine

    function iget( self ) result( iarr )
        class(arr), intent(in) :: self
        integer, allocatable   :: iarr(:)  !< output integer array
        if( allocated(self%iarr) )then
            allocate( iarr(size(self%iarr)), source=self%iarr ,stat=alloc_stat)
            if(alloc_stat /= 0) call allocchk(" In simple_arr::iget  iarr dealloc ", alloc_stat )
        else
            THROW_HARD('no info in iarr; iget')
        endif
    end function

    function rget( self ) result( rarr )
        class(arr), intent(in) :: self
        real, allocatable      :: rarr(:) !< output float array
        if( allocated(self%rarr) )then
            allocate( rarr(size(self%rarr)), source=self%rarr,STAT=alloc_stat)
            if(alloc_stat /= 0) call allocchk(" In simple_arr::rget  allocation fault ", alloc_stat )
        else
            THROW_HARD('no info in rarr; rget')
        endif
    end function

    subroutine display( self )
        class(arr), intent(in) :: self
        if( allocated(self%iarr) ) print *, self%iarr
        if( allocated(self%rarr) ) print *, self%rarr
    end subroutine display

    subroutine kill( self )
        class(arr), intent(inout) :: self
        if( allocated(self%iarr) ) then
            deallocate( self%iarr ,stat=alloc_stat)
            if(alloc_stat /= 0) call allocchk(" In simple_arr::kill  deallocate iarr fault ", alloc_stat)
        end if
          if( allocated(self%rarr) )then
              deallocate( self%rarr, stat=alloc_stat)
              if(alloc_stat /= 0) call allocchk(" In simple_arr::kill  deallocate rarr fault ", alloc_stat)
          end if
    end subroutine

end module simple_arr
