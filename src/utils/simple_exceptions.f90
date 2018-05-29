module simple_exceptions
use simple_defs
use simple_user_interface, only: simple_program, get_prg_ptr
use simple_cmdline,        only: cmdline
implicit none

public :: check_required_keys
private

contains

    subroutine check_required_keys( cline, prg )
        class(cmdline),   intent(in) :: cline
        character(len=*), intent(in) :: prg
        type(simple_program), pointer :: ptr2prg => null()
        type(str4arr),    allocatable :: keys_required(:), keys_in_cline(:)
        integer, allocatable :: matches(:)
        integer :: nreq, i, j, nkeys
        ! get pointer to program user interface
        call get_prg_ptr(prg, ptr2prg)
        ! get required keys
        keys_required = ptr2prg%get_required_keys()
        if( .not. allocated(keys_required) ) return
        nreq = size(keys_required)
        ! get keys in cline
        keys_in_cline = cline%get_keys()
        if( allocated(keys_in_cline) )then
            nkeys = size(keys_in_cline)
        else
            write(*,*) 'ERROR! no keys avaliable in command line'
            write(*,*) 'Program: ', trim(prg), ' requires # keys: ', nreq
            stop 'simple_exceptions :: check_required_keys'
        endif
        ! find matches
        allocate(matches(nreq), source=0)
        do i=1,nreq
            do j=1,nkeys
                if( keys_required(i)%str .eq. keys_in_cline(j)%str )then
                    matches(i) = j
                    exit
                endif
            end do
        end do
        if( any(matches == 0) )then
            ! raise exception
            write(*,*) 'ERROR! missing required inputs'
            write(*,*) 'Program: ', trim(prg), ' requires the following input keys:'
            do i=1,nreq
                if( matches(i) == 0 ) write(*,*) keys_required(i)%str
            end do
            stop 'simple_exceptions :: check_required_keys'
        endif
    end subroutine check_required_keys

end module simple_exceptions
