! particle picker iterator
module simple_picker_iter
include 'simple_lib.f08'
use simple_picker
use simple_parameters, only: params_glob
implicit none

public :: picker_iter

#include "simple_local_flags.inc"

private

type :: picker_iter
  contains
    procedure :: iterate
end type picker_iter

contains

    subroutine iterate( self, cline, smpd, moviename_intg, boxfile, nptcls_out, dir_out )
        use simple_cmdline, only: cmdline
        use simple_parameters, only : parameters
        class(picker_iter),        intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        real,                      intent(in)    :: smpd
        character(len=*),          intent(in)    :: moviename_intg
        character(len=LONGSTRLEN), intent(out)   :: boxfile
        integer,                   intent(out)   :: nptcls_out
        character(len=*),          intent(in)    :: dir_out
        if( .not. file_exists(moviename_intg) ) write(logfhandle,*) 'inputted micrograph does not exist: ', trim(adjustl(moviename_intg))
        write(logfhandle,'(a,1x,a)') '>>> PICKING MICROGRAPH:', trim(adjustl(moviename_intg))
        if( cline%defined('refs') .or. cline%defined('vol1') )then
            ! all good
        else
            THROW_HARD('Need references or volume!')
        endif
        if( cline%defined('thres') )then
            call init_picker(moviename_intg, params_glob%refs, smpd, lp_in=params_glob%lp,&
                &distthr_in=params_glob%thres, ndev_in=params_glob%ndev, dir_out=dir_out)
        else
            call init_picker(moviename_intg, params_glob%refs, smpd, lp_in=params_glob%lp, &
                &ndev_in=params_glob%ndev, dir_out=dir_out)
        endif
        call exec_picker(boxfile, nptcls_out)
        call kill_picker
    end subroutine iterate

end module simple_picker_iter
