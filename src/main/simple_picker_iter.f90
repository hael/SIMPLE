! particle picker iterator
module simple_picker_iter
include 'simple_lib.f08'
use simple_picker
use simple_segpicker,  only: segpicker
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

subroutine iterate( self, cline, moviename_intg, boxfile, nptcls_out, dir_out )
    use simple_cmdline, only: cmdline
    use simple_parameters, only : parameters
    class(picker_iter),        intent(inout) :: self
    class(cmdline),            intent(inout) :: cline
    character(len=*),          intent(in)    :: moviename_intg
    character(len=LONGSTRLEN), intent(out)   :: boxfile
    integer,                   intent(out)   :: nptcls_out
    character(len=*),          intent(in)    :: dir_out
    type(segpicker)  :: seg_picker
    type(parameters) :: params
    if( .not. file_exists(moviename_intg) )then
        write(logfhandle,*) 'inputted micrograph does not exist: ', trim(adjustl(moviename_intg))
    endif
    write(logfhandle,'(a,1x,a)') '>>> PICKING MICROGRAPH:', trim(adjustl(moviename_intg))
    if( cline%defined('refs') )then
        ! template based picking
        if( cline%defined('thres') )then
            call init_picker(moviename_intg, params_glob%refs, params_glob%smpd, lp_in=params_glob%lp,&
                &distthr_in=params_glob%thres, ndev_in=params_glob%ndev, dir_out=dir_out)
        else
            call init_picker(moviename_intg, params_glob%refs, params_glob%smpd, lp_in=params_glob%lp, ndev_in=params_glob%ndev, dir_out=dir_out)
        endif
        call exec_picker(boxfile, nptcls_out)
        call kill_picker
    else
        if( .not. cline%defined('min_rad') .or.  .not. cline%defined('max_rad'))then
            THROW_HARD('ERROR! min_rad and max_rad need to be present; exec_segpick')
        endif
        call params%new(cline)
        if( cline%defined('thres') .and. params%detector .ne. 'sobel')then
            THROW_HARD('ERROR! thres is compatible only with sobel detector; exec_segpick')
        endif
        if(cline%defined('draw_color')) then
            call seg_picker%new(moviename_intg, params_glob%min_rad, params_glob%max_rad,&
                &params_glob%smpd, params_glob%draw_color)
        else
            call seg_picker%new(moviename_intg, params_glob%min_rad, params_glob%max_rad,&
                &params_glob%smpd)
        endif
        if(cline%defined('detector')) then
            call seg_picker%preprocess_mic(params%detector)
        else
            call seg_picker%preprocess_mic(params_glob%detector)
        endif
        call seg_picker%identify_particle_positions()
        call seg_picker%output_identified_particle_positions()
        call seg_picker%write_boxfile()
        ! call segpick%print_info()
        call seg_picker%kill
    endif
end subroutine iterate

end module simple_picker_iter
