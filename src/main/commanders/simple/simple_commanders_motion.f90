!@descr: reference-based motion-model refinement commander
module simple_commanders_motion
use simple_commander_base, only: commander_base
implicit none

public :: commander_refine_motion_model
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_refine_motion_model
  contains
    procedure :: execute       => exec_refine_motion_model
end type commander_refine_motion_model

contains

    subroutine exec_refine_motion_model( self, cline )
        use simple_core_module_api,                    only: simple_end
        use simple_refine_motion_model_strategy,      only: refine_motion_model_strategy, &
                                                               create_refine_motion_model_strategy
        use simple_cmdline,                            only: cmdline
        use simple_parameters,                         only: parameters
        class(commander_refine_motion_model), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        class(refine_motion_model_strategy), allocatable :: strategy
        type(parameters)                                   :: params
        call cline%set('prg', 'refine_motion_model')
        strategy = create_refine_motion_model_strategy(cline)
        call strategy%apply_defaults(cline)
        call strategy%initialize(params, cline)
        call strategy%execute(params, cline)
        call strategy%finalize_run(params, cline)
        call strategy%cleanup(params, cline)
        call simple_end(strategy%end_message())
        if( allocated(strategy) ) deallocate(strategy)
    end subroutine exec_refine_motion_model

end module simple_commanders_motion
