!@descr: for pre-processing (motion correction, CTF estimation etc.)
module simple_commanders_preprocess
use simple_commanders_api
use simple_motion_correct_utils, only: flip_gain
use simple_mini_stream_utils,    only: segdiampick_preprocess
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_preprocess
  contains
    procedure :: execute      => exec_preprocess
end type commander_preprocess

type, extends(commander_base) :: commander_motion_correct
  contains
    procedure :: execute      => exec_motion_correct
end type commander_motion_correct

type, extends(commander_base) :: commander_gen_pspecs_and_thumbs
  contains
    procedure :: execute      => exec_gen_pspecs_and_thumbs
end type commander_gen_pspecs_and_thumbs

type, extends(commander_base) :: commander_ctf_estimate
  contains
    procedure :: execute      => exec_ctf_estimate
end type commander_ctf_estimate

contains

    ! Unified preprocess workflow (shared-memory/worker + distributed master)
    ! driven by runtime polymorphism (strategy pattern).
    subroutine exec_preprocess( self, cline )
        use simple_core_module_api,     only: simple_end
        use simple_preprocess_strategy, only: preprocess_strategy, create_preprocess_strategy
        use simple_cmdline,             only: cmdline
        use simple_parameters,          only: parameters
        class(commander_preprocess), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        class(preprocess_strategy), allocatable :: strategy
        type(parameters) :: params
        ! Ensure the distributed scripts see the correct program name.
        call cline%set('prg', 'preprocess')
        strategy = create_preprocess_strategy(cline)
        call strategy%initialize(params, cline)
        call strategy%execute(params, cline)
        call strategy%finalize_run(params, cline)
        call strategy%cleanup(params, cline)
        call simple_end(strategy%end_message())
        if( allocated(strategy) ) deallocate(strategy)
    end subroutine exec_preprocess

    ! Unified motion_correct workflow (worker/shared-memory + distributed master)
    ! driven by runtime polymorphism (strategy pattern).
    subroutine exec_motion_correct( self, cline )
        use simple_core_module_api,         only: simple_end
        use simple_motion_correct_strategy, only: motion_correct_strategy, create_motion_correct_strategy
        use simple_cmdline,                 only: cmdline
        use simple_parameters,              only: parameters
        class(commander_motion_correct), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        class(motion_correct_strategy), allocatable :: strategy
        type(parameters) :: params
        ! Helps distributed job script generation if it relies on 'prg'
        call cline%set('prg', 'motion_correct')
        strategy = create_motion_correct_strategy(cline)
        call strategy%initialize(params, cline)
        call strategy%execute(params, cline)
        call strategy%finalize_run(params, cline)
        call strategy%cleanup(params, cline)
        call simple_end(strategy%end_message())
        if( allocated(strategy) ) deallocate(strategy)
    end subroutine exec_motion_correct

    ! Unified gen_pspecs_and_thumbs workflow (worker/shared-memory + distributed master)
    ! driven by runtime polymorphism (strategy pattern).
    subroutine exec_gen_pspecs_and_thumbs( self, cline )
        use simple_core_module_api,                only: simple_end
        use simple_gen_pspecs_and_thumbs_strategy, only: gen_pspecs_and_thumbs_strategy, create_gen_pspecs_and_thumbs_strategy
        use simple_cmdline,                        only: cmdline
        use simple_parameters,                     only: parameters
        class(commander_gen_pspecs_and_thumbs), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline
        class(gen_pspecs_and_thumbs_strategy), allocatable :: strategy
        type(parameters) :: params
        ! Helps distributed job script generation if it relies on 'prg'
        call cline%set('prg', 'gen_pspecs_and_thumbs')
        strategy = create_gen_pspecs_and_thumbs_strategy(cline)
        call strategy%initialize(params, cline)
        call strategy%execute(params, cline)
        call strategy%finalize_run(params, cline)
        call strategy%cleanup(params, cline)
        call simple_end(strategy%end_message())
        if( allocated(strategy) ) deallocate(strategy)
    end subroutine exec_gen_pspecs_and_thumbs

    ! Unified ctf_estimate workflow (worker/shared-memory + distributed master)
    ! driven by runtime polymorphism (strategy pattern).
    subroutine exec_ctf_estimate( self, cline )
        use simple_core_module_api,      only: simple_end
        use simple_ctf_estimate_strategy, only: ctf_estimate_strategy, create_ctf_estimate_strategy
        use simple_cmdline,              only: cmdline
        use simple_parameters,           only: parameters
        class(commander_ctf_estimate), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        class(ctf_estimate_strategy), allocatable :: strategy
        type(parameters) :: params
        ! Helps distributed job script generation if it relies on 'prg'
        call cline%set('prg', 'ctf_estimate')
        strategy = create_ctf_estimate_strategy(cline)
        call strategy%initialize(params, cline)
        call strategy%execute(params, cline)
        call strategy%finalize_run(params, cline)
        call strategy%cleanup(params, cline)
        call simple_end(strategy%end_message())
        if( allocated(strategy) ) deallocate(strategy)
    end subroutine exec_ctf_estimate

end module simple_commanders_preprocess
