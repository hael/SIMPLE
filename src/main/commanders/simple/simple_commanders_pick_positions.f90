!@descr: pick_positions: load an MRC frame for native GTK picking of initial nanoparticle tracking boxes
module simple_commanders_pick_positions
use simple_commander_base, only: commander_base
use simple_cmdline,        only: cmdline
use simple_error,          only: simple_exception
implicit none
private
public :: commander_pick_positions
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_pick_positions
  contains
    procedure :: execute => exec_pick_positions
end type commander_pick_positions

contains

    subroutine exec_pick_positions( self, cline )
        use simple_defs,                  only: dp
        use simple_defs_gtk,              only: GTKUI_PICK_BOX_START
        use simple_parameters,            only: parameters
        use simple_image,                 only: image
        use simple_imghead,               only: find_ldim_nptcls, check_stack_dims
        use simple_stack_display,         only: stack_display
        use simple_manual_picker,         only: manual_picker
        use simple_gtk_utils,             only: gtkui_display_available
        use simple_pick_positions_window, only: run_pick_positions_window
        use simple_jiffys,                only: simple_end
        class(commander_pick_positions), intent(inout) :: self
        class(cmdline), intent(inout) :: cline
        class(parameters), allocatable :: params
        type(image) :: imgs(1)
        type(stack_display), target :: disp
        type(manual_picker), target :: picks
        integer :: dims(3), nframes, frame
        if( cline%defined('box') ) THROW_HARD('pick_positions: choose the box size with the slider, not box=')
        call cline%set('mkdir', 'no')
        if( .not. cline%defined('fromf') ) call cline%set('fromf', 0)
        if( .not. cline%defined('outfile') ) call cline%set('outfile', 'positions_all.box')
        if( .not. gtkui_display_available() ) THROW_HARD('pick_positions needs a desktop session (locally or over ssh -X)')
        allocate(params)
        call params%new(cline)
        call find_ldim_nptcls(params%stk, dims, nframes)
        call check_stack_dims(nframes, dims, 'pick_positions')
        frame = params%fromf
        if( frame == 0 ) frame = nframes/2+1
        if( frame < 1 .or. frame > nframes ) THROW_HARD('fromf must identify an existing frame (one-based)')
        call picks%new(dims(:2), min(GTKUI_PICK_BOX_START, minval(dims(:2))))
        call imgs(1)%new([dims(1),dims(2),1], 1.)
        call imgs(1)%read(params%stk, frame)
        call disp%new(imgs, [.true.], auto_percentiles=[1._dp,99._dp])
        call imgs(1)%kill
        if( .not. disp%is_usable() ) THROW_HARD('the selected frame has no finite pixels')
        call run_pick_positions_window(picks, disp, params%stk%to_char(), frame, params%outfile%to_char())
        call picks%kill
        call disp%kill
        call simple_end('**** SIMPLE_PICK_POSITIONS NORMAL STOP ****')
    end subroutine exec_pick_positions

end module simple_commanders_pick_positions
