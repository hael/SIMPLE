!@descr: supporting 3D orientation search
module simple_commanders_refine3D
use simple_commanders_api
use simple_pftc_srch_api
use simple_commanders_euclid
use simple_commanders_volops, only: commander_postprocess
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: nspace_commander
 contains
   procedure :: execute      => exec_nspace
end type nspace_commander

type, extends(commander_base) :: commander_refine3D_auto
  contains
    procedure :: execute      => exec_refine3D_auto
end type commander_refine3D_auto

type, extends(commander_base) :: commander_refine3D
  contains
    procedure :: execute      => exec_refine3D
end type commander_refine3D

type, extends(commander_base) :: commander_refine3D_distr_worker
  contains
    procedure :: execute      => exec_refine3D_distr_worker
end type commander_refine3D_distr_worker

contains

    subroutine exec_nspace(self,cline)
        class(nspace_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters) :: params
        type(oris)       :: o
        real             :: ares
        integer          :: i
        call params%new(cline)
        do i=500,5000,500
            o = oris(i, is_ptcl=.false.)
            call o%spiral
            ares = o%find_angres()
            write(logfhandle,'(A,1X,I7,1X,A,1X,F5.2)') 'NR OF PROJDIRS:', i, 'RESOLUTION:', resang(ares, params%moldiam)
        end do
        call simple_end('**** SIMPLE_NSPACE NORMAL STOP ****')
    end subroutine exec_nspace

    subroutine exec_refine3D_auto( self, cline )
        use simple_commanders_rec, only: commander_rec3D
        class(commander_refine3D_auto), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(cmdline)               :: cline_rec3D
        type(commander_postprocess) :: xpostprocess
        type(parameters)            :: params
        real,    parameter :: LP2SMPD_TARGET   = 1./3.
        real,    parameter :: SMPD_TARGET_MIN  = 1.3
        logical, parameter :: DEBUG            = .true.
        integer, parameter :: MINBOX           = 256
        integer, parameter :: NPDIRS4BAL       = 300
        type(string) :: str_state
        real         :: smpd_target, smpd_crop, scale, trslim
        integer      :: box_crop, maxits_phase1, maxits_phase2, iter
        logical      :: l_autoscale
        ! commanders
        type(commander_rec3D)    :: xrec3D
        type(commander_refine3D) :: xrefine3D
        ! hard defaults
        call cline%set('balance',         'no') ! balanced particle sampling based on available 3D solution
        call cline%set('greedy_sampling', 'no') ! stochastic within-class selection without consideration to objective function value
        call cline%set('trail_rec',      'yes') ! trailing average 3D reconstruction
        call cline%set('refine',  'prob_neigh') ! greedy multi-neighborhood 3D refinement 
        call cline%set('ml_reg',         'yes') ! ML regularization is on
        call cline%set('automsk',        'yes') ! envelope masking for background flattening
        call cline%set('overlap',         0.99) ! convergence if overlap > 99%
        call cline%set('nstates',            1) ! only single-state refinement is supported
        call cline%set('objfun',      'euclid') ! the objective function is noise-normalized Euclidean distance
        call cline%set('envfsc',         'yes') ! we use the envelope mask when calculating an FSC plot
        call cline%set('lplim_crit',     0.143) ! we use the 0.143 criterion for low-pass limitation
        call cline%set('incrreslim',      'no') ! if anything 'yes' makes it slightly worse, but no real difference right now
        ! overridable defaults
        if( .not. cline%defined('mkdir')       ) call cline%set('mkdir',        'yes')
        if( .not. cline%defined('center')      ) call cline%set('center',        'no') ! 4 now, probably fine
        if( .not. cline%defined('sigma_est')   ) call cline%set('sigma_est', 'global') ! 4 now, probably fine
        if( .not. cline%defined('combine_eo')  ) call cline%set('combine_eo',    'no') ! 4 now, to allow more rapid testing
        if( .not. cline%defined('prob_inpl')   ) call cline%set('prob_inpl',    'yes') ! no difference at this stage, so prefer 'yes'
        if( .not. cline%defined('update_frac') ) call cline%set('update_frac',    0.1) ! 4 now, needs testing/different logic (nsample?)
        if( .not. cline%defined('ml_reg')      ) call cline%set('ml_reg',       'yes') ! better map with ml_reg='yes'
        if( .not. cline%defined('filt_mode')   ) call cline%set('filt_mode',   'none')
        ! works, should be considered if the defaults are not satisfactory
        if( .not. cline%defined('maxits')      ) call cline%set('maxits',          20) ! ~2 passes over particles, sufficient ?
        if( .not. cline%defined('keepvol')     ) call cline%set('keepvol',       'no') ! we do not keep volumes for each iteration by deafult
        call params%new(cline)
        call cline%set('maxits_glob', params%maxits) ! needed for correct lambda annealing
        call cline%set('mkdir', 'no') ! to avoid nested directory structure
        if( mod(params%maxits,2) /= 0) THROW_HARD('Maximum number of iterations (MAXITS) need to be divisible with 2')
        maxits_phase1 = params%maxits / 2
        maxits_phase2 = maxits_phase1
        if( params%box <= MINBOX )then
            smpd_target = params%smpd
            smpd_crop   = params%smpd
            box_crop    = params%box
            scale       = 1.0
            l_autoscale = .false.
        else
            smpd_target = max(SMPD_TARGET_MIN, params%res_target * LP2SMPD_TARGET)
            call autoscale(params%box, params%smpd, smpd_target, box_crop, smpd_crop, scale, minbox=MINBOX)
            l_autoscale = box_crop < params%box
        endif
        trslim = min(8.,max(2.0, AHELIX_WIDTH / smpd_crop))
        if( DEBUG )then
            print *, 'smpd_target: ', smpd_target
            print *, 'box:         ', params%box
            print *, 'box_crop:    ', box_crop
            print *, 'smpd:        ', params%smpd
            print *, 'smpd_crop:   ', smpd_crop
            print *, 'scale:       ', scale
            print *, 'trslim:      ', trslim
            print *, 'l_autoscale: ', l_autoscale
        endif
        call cline%set('trs', trslim)
        if( l_autoscale )then
            call cline%set('box_crop',  box_crop)
            call cline%set('smpd_crop', smpd_crop)
        endif
        ! generate an initial 3D reconstruction
        cline_rec3D = cline
        call cline_rec3D%set('prg', 'reconstruct3D') ! required for distributed call
        call cline_rec3D%delete('trail_rec')
        call cline_rec3D%delete('nspace_next')
        call cline_rec3D%delete('objfun')
        call cline_rec3D%delete('sigma_est')
        call cline_rec3D%delete('update_frac')
        call cline_rec3D%set('objfun', 'cc') ! ugly, but this is how it works in parameters
        call xrec3D%execute(cline_rec3D)
        ! 3D refinement, phase1
        str_state = int2str_pad(1,2)
        call cline%set('vol1', string(VOL_FBODY)//str_state//MRC_EXT)
        call cline%set('prg',                'refine3D')
        call cline%set('ufrac_trec', params%update_frac)
        call cline%set('maxits',          maxits_phase1)
        call cline%set('filt_mode',           'uniform')
        call xrefine3D%execute(cline)
        ! iteration number bookkeeping
        iter = 0
        if( cline%defined('endit') )then
            iter = cline%get_iarg('endit')
            call cline%delete('endit')
        endif
        iter = iter + 1
        ! re-reconstruct from all particle images
        call xrec3D%execute(cline_rec3D)
        ! 3D refinement, phase2
        call cline%set('vol1', string(VOL_FBODY)//str_state//MRC_EXT)
        call cline%set('maxits',   maxits_phase1)
        call cline%set('filt_mode', params%filt_mode)
        call cline%set('startit',           iter)
        call cline%set('which_iter',        iter)
        call xrefine3D%execute(cline)
        ! re-reconstruct from all particle images
        call xrec3D%execute(cline_rec3D)
        ! postprocess
        call cline%set('prg', 'postprocess')
        call cline%set('mkdir', 'yes')
        call xpostprocess%execute(cline)        
    end subroutine exec_refine3D_auto

    !> Single entrypoint (shared-memory OR distributed master), driven by a strategy.
    subroutine exec_refine3D( self, cline )
        use simple_core_module_api
        use simple_refine3D_strategy
        class(commander_refine3D), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        class(refine3D_strategy), allocatable :: strategy
        type(parameters) :: params
        type(builder)    :: build
        logical          :: converged
        integer          :: niters
        ! local defaults (kept consistent with previous distributed master)
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('cenlp')   ) call cline%set('cenlp',        30.)
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call cline%set('prg', 'refine3D')
        ! Select execution strategy (shared-memory vs distributed master)
        strategy = create_refine3D_strategy(cline)
        call strategy%initialize(params, build, cline)
        ! Main loop counter semantics:
        !   - params%maxits is the *number of iterations to run* in this invocation.
        !   - params%which_iter starts at params%startit.
        niters            = 0
        params%which_iter = params%startit - 1
        if( cline%defined('extr_iter') )then
            params%extr_iter = params%extr_iter - 1
        else
            params%extr_iter = params%startit - 1
        endif
        do
            niters            = niters + 1
            params%which_iter = params%which_iter + 1
            params%extr_iter  = params%extr_iter  + 1
            call strategy%execute_iteration(params, build, cline, converged)
            call strategy%finalize_iteration(params, build)
            if( converged .or. niters >= params%maxits ) exit
        end do
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params)
        if( allocated(strategy) ) deallocate(strategy)
        ! Global teardown (strategies may have built different toolboxes)
        call build%kill_strategy3D_tbox
        call build%kill_general_tbox
        call build%pftc%kill
        call simple_end('**** SIMPLE_REFINE3D NORMAL STOP ****')
    end subroutine exec_refine3D

    !> Distributed worker (single-iteration execution). This should be the command
    !> invoked by the scheduler for each partition.
    subroutine exec_refine3D_distr_worker( self, cline )
        use simple_core_module_api
        use simple_strategy3D_matcher, only: refine3D_exec
        class(commander_refine3D_distr_worker), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        logical          :: converged
        logical          :: l_write_partial_recs
        ! Flags required for worker execution
        if( .not. cline%defined('part')    ) THROW_HARD('PART must be defined for distributed worker execution')
        if( .not. cline%defined('outfile') ) THROW_HARD('OUTFILE must be defined for distributed worker execution')
        ! Worker needs the alignment toolboxes
        call build%init_params_and_build_strategy3D_tbox(cline, params)
        if( params%which_iter < 1 ) params%which_iter = max(1, params%startit)
        if( .not. cline%defined('extr_iter') ) params%extr_iter = params%which_iter
        call cline%set('which_iter', int2str(params%which_iter))
        l_write_partial_recs = trim(params%volrec) .eq. 'yes' .or. params%l_polar
        call refine3D_exec(params, build, cline, params%which_iter, converged, l_write_partial_recs)
        call build%kill_strategy3D_tbox
        call build%kill_general_tbox
    end subroutine exec_refine3D_distr_worker

end module simple_commanders_refine3D
