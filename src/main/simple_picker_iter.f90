! particle picker iterator
module simple_picker_iter
include 'simple_lib.f08'
use simple_picker
use simple_parameters
use simple_picker_utils, only: picker_utils
use simple_cmdline,      only: cmdline
use simple_image,        only: image
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
        class(picker_iter),        intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        real,                      intent(in)    :: smpd
        character(len=*),          intent(in)    :: moviename_intg
        character(len=LONGSTRLEN), intent(out)   :: boxfile
        integer,                   intent(out)   :: nptcls_out
        character(len=*),          intent(in)    :: dir_out
        type(image)        :: micimg
        type(picker_utils) :: putils
        integer            :: ldim(3), ifoo
        if( .not. file_exists(moviename_intg) ) write(logfhandle,*) 'inputted micrograph does not exist: ', trim(adjustl(moviename_intg))
        write(logfhandle,'(a,1x,a)') '>>> PICKING MICROGRAPH:', trim(adjustl(moviename_intg))
        if( cline%defined('refs') .or. cline%defined('vol1') )then
            if( cline%defined('thres') )then
                call init_picker(moviename_intg, params_glob%refs, smpd, lp_in=params_glob%lp,&
                    &distthr_in=params_glob%thres, ndev_in=params_glob%ndev, dir_out=dir_out)
            else
                call init_picker(moviename_intg, params_glob%refs, smpd, lp_in=params_glob%lp, &
                    &ndev_in=params_glob%ndev, dir_out=dir_out)
            endif
            call exec_picker(boxfile, nptcls_out)
            call kill_picker
        else if( cline%defined('moldiam') )then
            call find_ldim_nptcls(moviename_intg, ldim, ifoo)
            call micimg%new(ldim, smpd)
            call micimg%read(moviename_intg)
            call putils%exec_gaupicker(micimg, smpd, params_glob%moldiam, params_glob%pcontrast, moviename_intg, boxfile, nptcls_out, dir_out=dir_out)
            call micimg%kill
        else
            THROW_HARD('Need 2D references or volume or modiam!')
        endif
    end subroutine iterate

end module simple_picker_iter
