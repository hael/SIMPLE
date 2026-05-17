!@descr: supporting 3D orientation search
module simple_commanders_refine3D
use simple_commanders_api
use simple_pftc_srch_api
use simple_refine3D_fnames,   only: refine3D_state_vol_fname
implicit none
#include "simple_local_flags.inc"

integer, parameter :: NSAMPLE_REFINE3D_AUTO = 25000

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
        type(parameters)            :: params
        type(sp_project)            :: spproj
        type(string)                :: init_vol
        real,    parameter :: SMPD_TARGET_DEFAULT = 1.3
        logical, parameter :: DEBUG  = .true.
        integer, parameter :: MINBOX = 256
        integer, parameter :: MINITS_REFINE3D_AUTO = 10
        real    :: smpd_target, smpd_crop, scale, trslim, init_smpd, update_frac_auto
        integer :: box_crop, init_box, nptcls_eff, nsample_target
        logical :: l_autoscale, l_have_init_vol, l_maxits_defined
        ! commanders
        type(commander_rec3D)    :: xrec3D
        type(commander_refine3D) :: xrefine3D
        ! hard defaults
        call cline%set('balance',         'no') ! no balancing based on 2D clustering
        call cline%set('greedy_sampling', 'no') ! only active when balance is 'yes'`
        call cline%set('trail_rec',      'yes') ! trailing average 3D reconstruction
        call cline%set('refine',  'prob_neigh') ! probabilistioc neighborhood 3D refinement 
        call cline%set('ml_reg',         'yes') ! ML regularization is on
        call cline%set('automsk',        'yes') ! envelope masking for background flattening
        call cline%set('overlap',         0.99) ! convergence if overlap > 99%
        call cline%set('nstates',            1) ! only single-state refinement is supported
        call cline%set('objfun',      'euclid') ! the objective function is noise-normalized Euclidean distance
        call cline%set('envfsc',         'yes') ! we use the envelope mask when calculating the FSC
        call cline%set('lplim_crit',     0.143) ! we use the 0.143 criterion for low-pass limitation
        call cline%set('incrreslim',      'no') ! if anything 'yes' makes it slightly worse, but no real difference right now
        ! overridable defaults
        if( .not. cline%defined('mkdir')       ) call cline%set('mkdir',            'yes')
        if( .not. cline%defined('center')      ) call cline%set('center',            'no') ! 4 now, probably fine
        if( .not. cline%defined('sigma_est')   ) call cline%set('sigma_est',     'global') ! 4 now, probably fine
        if( .not. cline%defined('combine_eo')  ) call cline%set('combine_eo',        'no') ! 4 now, to allow more rapid testing
        if( .not. cline%defined('prob_inpl')   ) call cline%set('prob_inpl',        'yes') ! no difference at this stage, so prefer 'yes'
        if( .not. cline%defined('nsample')     ) call cline%set('nsample', NSAMPLE_REFINE3D_AUTO)
        if( .not. cline%defined('ml_reg')      ) call cline%set('ml_reg',           'yes') ! better map with ml_reg='yes'
        if( .not. cline%defined('filt_mode')   ) call cline%set('filt_mode', 'nonuniform') ! obvioulsy
        if( .not. cline%defined('nu_refine')   ) call cline%set('nu_refine',        'yes') ! allow conservative NU resolution-bank expansion
        l_maxits_defined = cline%defined('maxits')
        if( cline%defined('minits') )then
            call cline%set('minits', max(MINITS_REFINE3D_AUTO, cline%get_iarg('minits')))
        else
            call cline%set('minits', MINITS_REFINE3D_AUTO)
        endif
        if( .not. cline%defined('keepvol')     ) call cline%set('keepvol',           'no') ! we do not keep volumes for each iteration by deafult
        call params%new(cline)
        call cline%set('mkdir', 'no') ! to avoid nested directory structure
        call set_refine3D_auto_sampling()
        l_have_init_vol = .false.
        init_box        = 0
        init_smpd       = 0.
        if( cline%defined('vol1') )then
            init_vol = cline%get_carg('vol1')
            if( .not. file_exists(init_vol) )then
                THROW_HARD('File: '//init_vol%to_char()//' does not exist! refine3D_auto')
            endif
            l_have_init_vol = .true.
            write(logfhandle,'(A,1X,A)') '>>> REFINE3D_AUTO USING INPUT VOLUME:', init_vol%to_char()
        else
            call spproj%read_segment('out', params%projfile)
            if( spproj%isthere_in_osout('vol', 1) )then
                call spproj%get_vol('vol', 1, init_vol, init_smpd, init_box)
                if( file_exists(init_vol) )then
                    if( project_init_vol_compatible() )then
                        l_have_init_vol = .true.
                        write(logfhandle,'(A,1X,A)') '>>> REFINE3D_AUTO USING PROJECT VOLUME:', init_vol%to_char()
                        write(logfhandle,'(A,I0,A,F8.4)') '>>> PROJECT VOLUME BOX/SMPD: ', init_box, '/', init_smpd
                    else
                        write(logfhandle,'(A,1X,A)') '>>> REFINE3D_AUTO PROJECT VOLUME SAMPLING MISMATCH, RECONSTRUCTING:', init_vol%to_char()
                        write(logfhandle,'(A,I0,A,F8.4)') '>>> PROJECT VOLUME BOX/SMPD: ', init_box, '/', init_smpd
                        write(logfhandle,'(A,I0,A,F8.4)') '>>> CURRENT RUN BOX/SMPD:    ', params%box, '/', params%smpd
                    endif
                else
                    write(logfhandle,'(A,1X,A)') '>>> REFINE3D_AUTO PROJECT VOLUME MISSING, RECONSTRUCTING:', init_vol%to_char()
                endif
            endif
            call spproj%kill
        endif
        smpd_target = SMPD_TARGET_DEFAULT
        if( params%box <= MINBOX )then
            smpd_target = params%smpd
            smpd_crop   = params%smpd
            box_crop    = params%box
            scale       = 1.0
            l_autoscale = .false.
        else
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
        call cline_rec3D%delete('objfun')
        call cline_rec3D%delete('sigma_est')
        call cline_rec3D%delete('update_frac')
        call cline_rec3D%set('objfun', 'cc') ! ugly, but this is how it works in parameters
        if( l_have_init_vol )then
            call cline%set('vol1', init_vol)
        else
            call xrec3D%execute(cline_rec3D)
            call cline%set('vol1', refine3D_state_vol_fname(1))
        endif
        ! 3D refinement iterations
        call cline%set('prg',                   'refine3D')
        call cline%set('ufrac_trec',    params%update_frac)
        call cline%set('maxits',             params%maxits)
        call xrefine3D%execute(cline)
        ! re-reconstruct from all particle images
        call cline_rec3D%set('postprocess', 'yes')
        call xrec3D%execute(cline_rec3D)       
        call init_vol%kill

    contains

        logical function project_init_vol_compatible() result( l_compatible )
            l_compatible = init_box == params%box .and. init_smpd > TINY .and. &
                &abs(init_smpd - params%smpd) <= 1.e-6
        end function project_init_vol_compatible

        subroutine set_refine3D_auto_sampling()
            type(sp_project) :: sampling_proj
            integer :: maxits_auto, nptcls_per_iter
            nsample_target = params%nsample
            if( nsample_target < 1 ) THROW_HARD('nsample must be >= 1 for refine3D_auto')
            call sampling_proj%read(params%projfile)
            nptcls_eff = sampling_proj%count_state_gt_zero()
            call sampling_proj%kill
            if( nptcls_eff < 1 ) THROW_HARD('no active particles available for refine3D_auto')
            nptcls_per_iter = min(nptcls_eff, nsample_target)
            if( nptcls_eff <= nsample_target )then
                params%update_frac   = 1.0
                params%l_update_frac = .false.
                params%l_trail_rec   = .false.
                call cline%delete('update_frac')
                write(logfhandle,'(A,I0,A,I0,A)') '>>> REFINE3D_AUTO ACTIVE PARTICLES/SAMPLE TARGET: ', &
                    nptcls_eff, '/', nsample_target, ' -> FULL UPDATE'
            else
                update_frac_auto = real(nsample_target) / real(nptcls_eff)
                if( update_frac_auto <= 0.99 )then
                    params%update_frac   = update_frac_auto
                    params%l_update_frac = .true.
                    params%l_trail_rec   = trim(params%trail_rec).eq.'yes'
                    call cline%set('update_frac', update_frac_auto)
                    write(logfhandle,'(A,I0,A,I0,A,F8.4)') '>>> REFINE3D_AUTO ACTIVE PARTICLES/SAMPLE TARGET/UPDATE_FRAC: ', &
                        nptcls_eff, '/', nsample_target, '/', update_frac_auto
                else
                    params%update_frac   = 1.0
                    params%l_update_frac = .false.
                    params%l_trail_rec   = .false.
                    call cline%delete('update_frac')
                    write(logfhandle,'(A,I0,A,I0,A)') '>>> REFINE3D_AUTO ACTIVE PARTICLES/SAMPLE TARGET: ', &
                        nptcls_eff, '/', nsample_target, ' -> FULL UPDATE'
                endif
            endif
            if( .not. l_maxits_defined )then
                maxits_auto = ceiling((2.0 * real(nptcls_eff)) / real(nptcls_per_iter))
                maxits_auto = max(params%minits, min(50, max(2, maxits_auto)))
                params%maxits = maxits_auto
                call cline%set('maxits', params%maxits)
                write(logfhandle,'(A,I0,A,I0,A)') '>>> REFINE3D_AUTO MAXITS FOR ~2 UPDATES/PARTICLE: ', &
                    &params%maxits, ' (MINIMUM: ', params%minits, ')'
            else if( params%maxits < params%minits )then
                params%maxits = params%minits
                call cline%set('maxits', params%maxits)
                write(logfhandle,'(A,I0)') '>>> REFINE3D_AUTO MAXITS RAISED TO MINIMUM ITERATIONS: ', params%maxits
            endif
        end subroutine set_refine3D_auto_sampling

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
        ! sanity check: multiple input volumes require nstates > 1
        if( cline%defined('vol2') )then
            if( .not. cline%defined('nstates') )then
                THROW_HARD('Multiple volumes (vol1, vol2, ...) provided on command line but NSTATES is not set; set NSTATES to the number of states')
            else if( cline%get_iarg('nstates') <= 1 )then
                THROW_HARD('Multiple volumes (vol1, vol2, ...) provided on command line but NSTATES <= 1; set NSTATES to the number of states')
            endif
        endif
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
        class(cmdline),                         intent(inout) :: cline
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
        l_write_partial_recs = trim(params%volrec) .eq. 'yes'
        call refine3D_exec(params, build, cline, params%which_iter, converged, l_write_partial_recs)
        call build%kill_strategy3D_tbox
        call build%kill_general_tbox
    end subroutine exec_refine3D_distr_worker

end module simple_commanders_refine3D
