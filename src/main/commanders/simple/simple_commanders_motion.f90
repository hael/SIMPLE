!@descr: fractionate_model commander: saved-motion fraction re-extraction and interval reconstruction
module simple_commanders_motion
use simple_commander_base, only: commander_base
implicit none

public :: commander_fractionate_model
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_fractionate_model
  contains
    procedure :: execute       => exec_fractionate_model
end type commander_fractionate_model

contains

    subroutine exec_fractionate_model( self, cline )
        use simple_error,                              only: simple_exception
        use simple_core_module_api,                    only: simple_end
        use simple_fractionate_model_strategy,         only: fractionate_model_strategy, &
                                                               create_fractionate_model_strategy
        use simple_string,                             only: string
        use simple_cmdline,                            only: cmdline
        use simple_parameters,                         only: parameters
        class(commander_fractionate_model), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        class(fractionate_model_strategy), allocatable :: strategy
        type(parameters) :: params
        type(string)     :: oritype
        call cline%set('prg', 'fractionate_model')
        call cline%set('wfloat16', 'yes')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        oritype = cline%get_carg('oritype')
        select case(oritype%to_char())
            case('ptcl2D')
                ! User averaging options are retained for make_cavgs.
            case('ptcl3D')
                if( .not.cline%defined('pgrp') ) THROW_HARD('pgrp is required for oritype=ptcl3D')
            case default
                THROW_HARD('fractionate_model requires oritype=ptcl2D or ptcl3D')
        end select
        if( .not. cline%defined('mkdir')        ) call cline%set('mkdir',          'yes')
        if( .not. cline%defined('pcontrast')    ) call cline%set('pcontrast',    'black')
        if( .not. cline%defined('fromf')        ) call cline%set('fromf',              1)
        if( .not. cline%defined('tof')          ) call cline%set('tof',                0)
        if( oritype == 'ptcl3D' )then
            ! These reconstructions use the saved poses with no noise-model weighting or refinement.
            call cline%set('objfun', 'cc')
            call cline%set('ml_reg', 'no')
            call cline%set('cc_emit_sigma', 'no')
            call cline%set('rec_backend', 'gridding')
            call cline%set('filt_mode', 'none')
            call cline%set('automsk', 'no')
            call cline%set('postprocess', 'yes')
            call cline%set('bfac', 0.)
            call cline%set('mirr', 'no')
            call cline%set('m_estimator', 'no')
            call cline%set('update_frac', 1.)
            call cline%set('trail_rec', 'no')
            call cline%set('trs', 5.)
            call cline%delete('refine')
        endif
        strategy = create_fractionate_model_strategy(cline)
        call strategy%initialize(params, cline)
        call strategy%execute(params, cline)
        call strategy%finalize_run(params, cline)
        call strategy%cleanup(params, cline)
        call simple_end(strategy%end_message())
        if( allocated(strategy) ) deallocate(strategy)
    end subroutine exec_fractionate_model

end module simple_commanders_motion
