!@descr: simultaneous 2D alignment and clustering of single-particle images
module simple_commanders_refine2D
use simple_commanders_api
use simple_pftc_srch_api
use simple_classaverager
use simple_commanders_cavgs,   only: commander_rank_cavgs
use simple_commanders_mkcavgs, only: commander_make_cavgs, commander_make_cavgs_distr, commander_cavgassemble
use simple_gui_utils,          only: mrc2jpeg_tiled
use simple_procimgstk,         only: selection_from_tseries_imgfile, random_selection_from_imgfile, copy_imgfile, noise_imgfile
use simple_progress,           only: progressfile_init, progressfile_update
use simple_commanders_imgops,  only: commander_scale
use simple_matcher_smpl_and_lplims, only: refine2D_requires_full_assignment
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_refine2D
  contains
    procedure :: execute      => exec_refine2D
end type commander_refine2D

type, extends(commander_base) :: commander_refine2D_distr_worker
  contains
    procedure :: execute      => exec_refine2D_distr_worker
end type commander_refine2D_distr_worker

contains

    !> Single entrypoint (shared-memory OR distributed master), driven by a strategy.
    subroutine exec_refine2D( self, cline )
        use simple_refine2D_strategy
        use simple_eul_prob_tab_utils, only: prob_beta_at
        use simple_defs,               only: PROB_BETA_TAIL
        class(commander_refine2D), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        class(refine2D_strategy), allocatable     :: strategy
        type(parameters)       :: params
        type(builder)          :: build
        logical                :: converged, l_anneal
        integer                :: niters, it_tail
        real                   :: beta_start
        ! local defaults (kept consistent with previous distributed master)
        call cline%set('prg', 'refine2D')
        call set_refine2D_defaults(cline)
        ! Select execution strategy (shared-memory vs distributed master)
        strategy = create_refine2D_strategy(cline)
        call strategy%initialize(params, build, cline)
        if( params%l_nonuniform ) THROW_HARD('nonuniform filtering is not available in 2D; exec_refine2D')
        ! the iteration jobs do not report to NICE
        if( cline%defined("niceserver") ) call cline%delete('niceserver')
        if( cline%defined("niceprocid") ) call cline%delete('niceprocid')
        if( trim(params%stream2d).eq.'no' ) call progressfile_init()
        ! Main loop counter semantics:
        !   - params%maxits is the *number of iterations to run* in this invocation.
        !   - params%which_iter starts at params%startit.
        niters            = 0
        params%which_iter = params%startit - 1
        if( cline%defined('extr_iter') )then
            params%extr_iter = params%extr_iter - 1
        else
            params%extr_iter = params%startit   - 1
        endif
        ! the probability tables' posterior draws, as exec_refine3D: the iteration's tempering factor on the command
        ! line for the workers, and the annealing tail to prob_beta_end once the class overlap converges
        ! (check_conv2D starts it)
        beta_start = params%prob_beta
        l_anneal   = params%l_prob_align_mode .and. abs(params%prob_beta_end - beta_start) > TINY
        do
            niters            = niters + 1
            params%which_iter = params%which_iter + 1
            params%extr_iter  = params%extr_iter  + 1
            if( params%l_prob_align_mode )then
                if( l_anneal )then
                    it_tail = strategy%anneal_tail_start()
                    if( it_tail > 0 )then
                        params%prob_beta = prob_beta_at(params%which_iter - it_tail + 1, PROB_BETA_TAIL,&
                            &beta_start, params%prob_beta_end)
                    else
                        params%prob_beta = beta_start
                    endif
                endif
                call cline%set('prob_beta', params%prob_beta)
            endif
            ! Strategy handles everything: alignment + cavgs + convergence
            call strategy%execute_iteration(params, build, cline, converged)
            call strategy%finalize_iteration(params, build)
            if( converged ) exit
            if( niters >= params%maxits .and. .not. refine2D_requires_full_assignment(params) ) exit
        end do
        ! Cleanup
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params)
        call qsys_cleanup(params)
        if( allocated(strategy) ) deallocate(strategy)
        ! Global teardown (strategies may have built different toolboxes)
        call build%kill_general_tbox()
        call build%kill_strategy2D_tbox()
        call simple_touch(REFINE2D_FINISHED)
        call simple_end('**** SIMPLE_REFINE2D NORMAL STOP ****')
    end subroutine exec_refine2D

    !> Distributed worker (single-iteration execution). This should be the command
    !> invoked by the scheduler for each partition.
    subroutine exec_refine2D_distr_worker( self, cline )
        use simple_strategy2D_matcher, only: refine2D_exec
        class(commander_refine2D_distr_worker), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        logical          :: l_frac
        call set_refine2D_defaults(cline)
        ! Flags required for worker execution
        if( .not. cline%defined('part')    ) THROW_HARD('PART must be defined for distributed worker execution')
        if( .not. cline%defined('outfile') ) THROW_HARD('OUTFILE must be defined for distributed worker execution')
        ! Worker needs the alignment toolboxes
        call build%init_params_and_build_strategy2D_tbox(cline, params, wthreads=.true.)
        if( params%l_nonuniform ) THROW_HARD('nonuniform filtering is not available in 2D; exec_refine2D_distr_worker')
        if( cline%defined('which_iter') )then
            params%which_iter = max(1, params%which_iter)
        else
            params%which_iter = max(1, params%startit)
        endif
        if( .not. cline%defined('extr_iter') ) params%extr_iter = params%which_iter
        call cline%set('which_iter', int2str(params%which_iter))
        call refine2D_exec(params, build, cline, params%which_iter, l_frac)
        call build%kill_strategy2D_tbox
        call build%kill_general_tbox
    end subroutine exec_refine2D_distr_worker

end module simple_commanders_refine2D
