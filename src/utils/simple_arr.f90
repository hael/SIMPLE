! array class (container class for the singly linked list)
module simple_arr
    use simple_defs
    use simple_syslib, only: alloc_errchk
implicit none

public :: arr
private

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

#include "simple_local_flags.inc"


contains
    subroutine new_1( self, iarr )
        class(arr), intent(inout) :: self
        integer, intent(in) :: iarr(:)        !< input integer array
        integer              :: err
        character(len=STDLEN):: io_msg
        if( allocated(self%iarr) )then
            deallocate( self%iarr,STAT=err,ERRMSG=io_msg)
            if(alloc_stat/=0)call alloc_errchk(" In simple_arr::new_1  deallocation fault "//trim(io_msg),err)
        end if
        allocate( self%iarr(size(iarr)), source=iarr, stat=err,errmsg=io_msg)
        if(alloc_stat/=0)call alloc_errchk(" In simple_arr::new_1  allocation fault "//trim(io_msg),err)
    end subroutine new_1

    subroutine new_2( self, rarr )
        class(arr), intent(inout) :: self
        real, intent(in) :: rarr(:)           !< input float array
        integer              :: err
        character(len=STDLEN):: io_msg
        if( allocated(self%rarr) ) then
            deallocate( self%rarr ,STAT=err,ERRMSG=io_msg)
            if(alloc_stat/=0)call alloc_errchk(" In simple_arr::new_2 deallocation fault"//trim(io_msg),err)
          end if
        allocate( self%rarr(size(rarr)), source=rarr,STAT=err,ERRMSG=io_msg )
        if(alloc_stat/=0)call alloc_errchk(" In simple_arr::new_2 allocation fault "//trim(io_msg),err)
    end subroutine

    function iget( self ) result( iarr )
        class(arr), intent(in) :: self
        integer, allocatable   :: iarr(:)  !< output integer array
        integer              :: err
        character(len=STDLEN):: io_msg
        if( allocated(self%iarr) )then
            allocate( iarr(size(self%iarr)), source=self%iarr ,stat=err,errmsg=io_msg)
            if(alloc_stat/=0)call alloc_errchk(" In simple_arr::iget  allocation fault "//trim(io_msg),err)
        else
            stop 'no info in iarr; get_1; simple_arr'
        endif
    end function

    function rget( self ) result( rarr )
        class(arr), intent(in) :: self
        real, allocatable      :: rarr(:) !< output float array
        integer              :: err
        character(len=STDLEN):: io_msg
        if( allocated(self%rarr) )then
            allocate( rarr(size(self%rarr)), source=self%rarr,STAT=err,ERRMSG=io_msg)
            if(alloc_stat/=0)call alloc_errchk(" In simple_arr::rget  allocation fault "//trim(io_msg),err)
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
        integer              :: err
        character(len=STDLEN):: io_msg
        if( allocated(self%iarr) ) then
            deallocate( self%iarr ,stat=err,errmsg=io_msg)
            if(alloc_stat/=0)call alloc_errchk(" In simple_arr::kill  deallocation fault "//trim(io_msg),err)
        end if
          if( allocated(self%rarr) )then
              deallocate( self%rarr, stat=err,errmsg=io_msg)
              if(alloc_stat/=0)call alloc_errchk(" In simple_arr::kill  deallocation fault "//trim(io_msg),err)
          end if
         
    end subroutine

end module simple_arr
