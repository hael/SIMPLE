module simple_pick_iter
use simple_picker
implicit none

public :: pick_iter
private

type :: pick_iter
  contains
    procedure :: iterate
end type pick_iter

contains

    subroutine iterate( self, cline, p, movie_counter, moviename_intg )
        use simple_params,       only: params
        use simple_oris,         only: oris
        use simple_cmdline,      only: cmdline
        use simple_filehandling, only: file_exists
        class(pick_iter),   intent(inout) :: self
        class(cmdline),     intent(in)    :: cline
        class(params),      intent(inout) :: p
        integer,            intent(inout) :: movie_counter
        character(len=*),   intent(in)    :: moviename_intg
        movie_counter = movie_counter + 1
        if( .not. file_exists(moviename_intg) )then
            write(*,*) 'inputted micrograph does not exist: ', trim(adjustl(moviename_intg))
        endif
        write(*,'(a,1x,a)') '>>> PICKING MICROGRAPH:', trim(adjustl(moviename_intg))
        if( cline%defined('thres') )then
            call init_picker(moviename_intg, p%refs, p%smpd, lp_in=p%lp, distthr_in=p%thres)
        else
            call init_picker(moviename_intg, p%refs, p%smpd, lp_in=p%lp)
        endif
        call exec_picker 
        call kill_picker
    end subroutine iterate

end module simple_pick_iter