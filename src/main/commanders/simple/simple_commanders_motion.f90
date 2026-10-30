!@descr: refine_motion_model commander: saved-motion fraction re-extraction and interval reconstruction
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
        if( .not. cline%defined('oritype')      ) call cline%set('oritype',     'ptcl3D')
        if( .not. cline%defined('mkdir')        ) call cline%set('mkdir',          'yes')
        if( .not. cline%defined('pcontrast')    ) call cline%set('pcontrast',    'black')
        if( .not. cline%defined('wfloat16')     ) call cline%set('wfloat16',        'no')
        if( .not. cline%defined('fromf')        ) call cline%set('fromf',              1)
        if( .not. cline%defined('tof')          ) call cline%set('tof',                0)
        ! These reconstructions use the saved poses with no noise-model weighting or refinement.
        call cline%set('objfun', 'cc')
        call cline%set('ml_reg', 'no')
        call cline%set('cc_emit_sigma', 'no')
        call cline%set('rec_backend', 'gridding')
        call cline%set('filt_mode', 'none')
        call cline%set('automsk', 'no')
        call cline%set('postprocess', 'no')
        call cline%set('m_estimator', 'no')
        call cline%set('update_frac', 1.)
        call cline%set('trail_rec', 'no')
        call cline%set('trs', 5.)
        call cline%delete('refine')
        strategy = create_refine_motion_model_strategy(cline)
        call strategy%initialize(params, cline)
        call strategy%execute(params, cline)
        call strategy%finalize_run(params, cline)
        call strategy%cleanup(params, cline)
        call simple_end(strategy%end_message())
        if( allocated(strategy) ) deallocate(strategy)
    end subroutine exec_refine_motion_model

end module simple_commanders_motion
