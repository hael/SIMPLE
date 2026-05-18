!@descr: small workflow-policy handoffs for nonuniform filter refinement
module simple_nu_filter_policy
use simple_core_module_api
implicit none

private
public :: nu_align_lp_fname, write_nu_align_lp_for_state, read_nu_align_lp_for_state

contains

    function nu_align_lp_fname( state ) result( fname )
        integer, intent(in) :: state
        type(string) :: fname
        fname = 'nu_align_lp_state'//int2str_pad(state,2)//'.txt'
    end function nu_align_lp_fname

    subroutine write_nu_align_lp_for_state( state, lp )
        integer, intent(in) :: state
        real,    intent(in) :: lp
        type(string) :: fname
        integer :: funit, io_stat
        fname = nu_align_lp_fname(state)
        call fopen(funit, status='REPLACE', action='WRITE', file=fname, iostat=io_stat)
        if( io_stat == 0 )then
            write(funit,'(F12.6)') lp
            call fclose(funit)
        else
            write(logfhandle,'(A,1X,A)') '>>> WARNING: failed to write NU alignment LP file:', fname%to_char()
        endif
        call fname%kill
    end subroutine write_nu_align_lp_for_state

    subroutine read_nu_align_lp_for_state( state, lp, found )
        integer, intent(in)  :: state
        real,    intent(out) :: lp
        logical, intent(out) :: found
        type(string) :: fname
        integer :: funit, io_stat
        lp    = 0.
        found = .false.
        fname = nu_align_lp_fname(state)
        if( .not.file_exists(fname) )then
            call fname%kill
            return
        endif
        call fopen(funit, status='OLD', action='READ', file=fname, iostat=io_stat)
        if( io_stat == 0 )then
            read(funit, *, iostat=io_stat) lp
            call fclose(funit)
        endif
        found = io_stat == 0 .and. lp > TINY
        call fname%kill
    end subroutine read_nu_align_lp_for_state

end module simple_nu_filter_policy
