! particle picker iterator
module simple_picker_iter
include 'simple_lib.f08'
use simple_picker
implicit none

public :: picker_iter
private

type :: picker_iter
  contains
    procedure :: iterate
end type picker_iter

contains

    subroutine iterate( self, cline, p, moviename_intg, boxfile, nptcls_out, dir_out )
        use simple_params,  only: params
        use simple_cmdline, only: cmdline
        class(picker_iter),    intent(inout)   :: self
        class(cmdline),        intent(in)      :: cline
        class(params),         intent(inout)   :: p
        character(len=*),      intent(in)      :: moviename_intg
        character(len=LONGSTRLEN), intent(out) :: boxfile
        integer,               intent(out)     :: nptcls_out
        character(len=*),      intent(in)      :: dir_out
        if( .not. file_exists(moviename_intg) )then
            write(*,*) 'inputted micrograph does not exist: ', trim(adjustl(moviename_intg))
        endif
        write(*,'(a,1x,a)') '>>> PICKING MICROGRAPH:', trim(adjustl(moviename_intg))
        if( cline%defined('thres') )then
            call init_picker(moviename_intg, p%refs, p%smpd, lp_in=p%lp,&
                &distthr_in=p%thres, ndev_in=p%ndev, dir_out=dir_out)
        else
            call init_picker(moviename_intg, p%refs, p%smpd, lp_in=p%lp, ndev_in=p%ndev, dir_out=dir_out)
        endif
        call exec_picker(boxfile, nptcls_out)
        call kill_picker
    end subroutine iterate

end module simple_picker_iter
