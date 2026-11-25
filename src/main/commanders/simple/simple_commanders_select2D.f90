!@descr: select2D: interactive class-average selection in a GTK 3 window, saved as a new project in <n>_selection/
module simple_commanders_select2D
use simple_commander_base, only: commander_base
use simple_cmdline,        only: cmdline
use simple_error,          only: simple_exception
implicit none

public :: commander_select2D
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_select2D
  contains
    procedure :: execute => exec_select2D
end type commander_select2D

contains

    ! mkdir=no: nothing is written until Save selection, which makes <n>_selection/ in the input's workspace root
    ! with the selected project, as prg=selection oritype=cls2D writes it (doc/implementation_notes/planned/select2D_viewer.md)
    subroutine exec_select2D( self, cline )
        use simple_parameters,      only: parameters
        use simple_sp_project,      only: sp_project
        use simple_chash,           only: chash
        use simple_string,          only: string
        use simple_image,           only: image
        use simple_imgarr_utils,    only: read_cavgs_into_imgarr, dealloc_imgarr
        use simple_stack_display,   only: stack_display
        use simple_cls2D_selection, only: cls2D_selection
        use simple_gtk_utils,       only: gtkui_display_available
        use simple_select2D_window, only: run_select2D_window
        use simple_jiffys,          only: simple_end
        class(commander_select2D), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(parameters)               :: params
        type(sp_project),      target  :: spproj
        type(cls2D_selection), target  :: sel
        type(stack_display),   target  :: disp
        type(chash)                    :: job_descr
        type(image), allocatable       :: imgs(:)
        call cline%set('mkdir', 'no')
        ! parse accepts every key SIMPLE knows: refuse append= rather than ignore it (decision 14)
        if( cline%defined('append') ) THROW_HARD('select2D has no append=: the class selection alone sets the particle states, as prg=selection append=no')
        call cline%set('oritype', 'cls2D')
        ! before params%new, so that a run without a display changes nothing
        if( .not. gtkui_display_available() ) THROW_HARD('select2D needs a desktop session: run it where a display is available (locally or over ssh -X)')
        call params%new(cline)
        ! the window needs the classes and their stack; the particles are read only when a selection is saved
        call spproj%read_segment('cls2D', params%projfile)
        call spproj%read_segment('out',   params%projfile)
        if( spproj%os_cls2D%get_noris() == 0 ) THROW_HARD('the project holds no class averages (its cls2D segment is empty)')
        call cline%gen_job_descr(job_descr, string('select2D'))
        call sel%new(spproj, params%projfile, job_descr)
        imgs = read_cavgs_into_imgarr(spproj)
        if( size(imgs) /= sel%get_ncls() ) THROW_HARD('the class-average stack and the cls2D segment hold different numbers of classes')
        call disp%new(imgs, sel%get_visible_mask())
        call dealloc_imgarr(imgs)
        call run_select2D_window(sel, disp, params%projname%to_char(), params%projfile%to_char())
        call disp%kill
        call sel%kill
        call job_descr%kill
        call spproj%kill
        call simple_end('**** SIMPLE_SELECT2D NORMAL STOP ****')
    end subroutine exec_select2D

end module simple_commanders_select2D
