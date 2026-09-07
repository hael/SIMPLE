!@descr: 3D reconstruction and associated things
module simple_commanders_rec
use simple_commanders_api
use simple_matcher_2Dprep
use simple_matcher_3Drec, only: calc_3Drec, calc_projdir3Drec
use simple_refine3D_fnames, only: refine3D_fsc_fname, refine3D_state_halfvol_fname, refine3D_state_vol_fname
use simple_sigma2_files, only: load_sigma2_groups
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_rec3D
  contains
    procedure :: execute => exec_rec3D
end type commander_rec3D

type, extends(commander_base) :: commander_bootstrap_rec3D
  contains
    procedure :: execute => exec_bootstrap_rec3D
end type commander_bootstrap_rec3D

type, extends(commander_base) :: commander_rec3D_worker
  contains
    procedure :: execute      => exec_rec3D_distr_worker
end type commander_rec3D_worker

type, extends(commander_base) :: random_rec_commander
  contains
    procedure :: execute      => exec_random_rec
end type random_rec_commander

contains

    subroutine exec_rec3D( self, cline )
        use simple_rec3D_strategy, only: rec3D_strategy, create_rec3D_strategy
        use simple_parameters,     only: parameters
        use simple_builder,        only: builder
        class(commander_rec3D), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        class(rec3D_strategy), allocatable :: strategy
        type(parameters) :: params
        type(builder)    :: build
        type(string)     :: rec_backend
        ! Commander-level defaults (apply to both modes)
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('trs')     ) call cline%set('trs', 5.)     ! to assure that shifts are being used
        if( .not. cline%defined('rec_backend') ) call cline%set('rec_backend', 'gridding')
        rec_backend = cline%get_carg('rec_backend')
        call cline%set('oritype', 'ptcl3D')
        call cline%delete('refine')
        ! Select and run strategy
        strategy = create_rec3D_strategy(cline)
        call strategy%initialize(params, build, cline)
        call strategy%execute(params, build, cline)
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params, build, cline)
        call rec_backend%kill
        ! End gracefully (single unified termination)
        call simple_end('**** SIMPLE_RECONSTRUCT3D NORMAL STOP ****', print_simple=.false.)
        if( allocated(strategy) ) deallocate(strategy)
    end subroutine exec_rec3D

    subroutine exec_bootstrap_rec3D( self, cline )
        use simple_commanders_euclid, only: commander_calc_pspec
        class(commander_bootstrap_rec3D), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(commander_rec3D) :: xrec3D
        type(commander_calc_pspec) :: xcalc_pspec
        type(cmdline)         :: cline_reg, cline_pspec
        type(parameters)      :: params
        integer               :: state, which_iter
        if( .not. cline%defined('mkdir')       ) call cline%set('mkdir',       'yes')
        call cline%set('oritype', 'ptcl3D')
        if( .not. cline%defined('nstates')     ) call cline%set('nstates',          1)
        call warn_for_forced_bootstrap_overrides(cline)
        call cline%set('sigma_est', 'global')
        if( .not. cline%defined('which_iter')  ) call cline%set('which_iter',       1)
        if( .not. cline%defined('postprocess') ) call cline%set('postprocess',  'yes')
        if( .not. cline%defined('combine_eo')  ) call cline%set('combine_eo',    'no')
        if( .not. cline%defined('envfsc')      ) call cline%set('envfsc',        'no')
        call cline%delete('objfun')
        call cline%delete('ml_reg')
        call params%new(cline)
        which_iter = max(1, params%which_iter)
        call cline%set('which_iter', which_iter)
        call cline%set('mkdir', 'no') ! child reconstruct3D calls must not create nested run directories
        ! One sigma2 basis for every bootstrap (2026-09-06): the particle
        ! power spectra, exactly what a fresh refinement seeds from, written
        ! as the grouped STAR of which_iter (legacy store) or the registered
        ! canonical state. The former half-map power estimator sat on a
        ! different basis than the residual sigmas a refinement then computes
        ! and conditioned the euclid system markedly worse (bgal residual
        ! 0.23 vs 0.08, refine3D_auto startup record). With the seed in hand
        ! the reconstruction is a single euclid ML-regularized pass; callers
        ! that ship a final map upgrade the seed with a residual pass
        ! (simple_sigma2_bootstrap).
        cline_pspec = cline
        call cline_pspec%set('prg',       'calc_pspec')
        call cline_pspec%set('mkdir',              'no')
        call cline_pspec%set('objfun',       'euclid')
        call cline_pspec%set('sigma_est',    'global')
        call cline_pspec%set('which_iter',  which_iter)
        call cline_pspec%delete('postprocess')
        call cline_pspec%delete('combine_eo')
        call cline_pspec%delete('rec_backend')
        call cline_pspec%delete('maxits_pcg')
        call cline_pspec%delete('rtol')
        call cline_pspec%delete('trail_seed')
        call cline_pspec%delete('outfile')
        write(logfhandle,'(A,I0)') '>>> BOOTSTRAP_REC3D SIGMA2 FROM PARTICLE POWER SPECTRA, ITERATION ', which_iter
        call xcalc_pspec%execute(cline_pspec)
        call cline_pspec%kill
        cline_reg = cline
        call prepare_bootstrap_rec_cline(cline_reg, l_regularized=.true.)
        write(logfhandle,'(A)') '>>> BOOTSTRAP_REC3D: EUCLID ML-REGULARIZED RECONSTRUCTION'
        call xrec3D%execute(cline_reg)
        call register_bootstrap_rec_outputs()
        do state = 1, params%nstates
            call cline%set('vol'//int2str(state), refine3D_state_vol_fname(state))
        enddo
        call cline_reg%kill
        call simple_end('**** SIMPLE_BOOTSTRAP_REC3D NORMAL STOP ****', print_simple=.false.)

    contains

        subroutine prepare_bootstrap_rec_cline( cline_rec, l_regularized )
            class(cmdline), intent(inout) :: cline_rec
            logical,        intent(in)    :: l_regularized
            integer :: state
            call cline_rec%set('prg',       'reconstruct3D')
            call cline_rec%set('mkdir',              'no')
            call cline_rec%set('oritype', params%oritype)
            call cline_rec%set('nstates', params%nstates)
            call cline_rec%set('sigma_est',     'global')
            call cline_rec%set('which_iter', which_iter)
            call cline_rec%set('trail_rec',        'no')
            call cline_rec%set('combine_eo',       'no')
            call cline_rec%delete('refine')
            call cline_rec%delete('update_frac')
            call cline_rec%delete('fillin')
            call cline_rec%delete('objfun_den')
            call cline_rec%delete('objfun_den_w')
            call cline_rec%delete('ufrac_trec')
            call cline_rec%delete('endit')
            call cline_rec%delete('vol_even')
            call cline_rec%delete('vol_odd')
            call cline_rec%delete('refs')
            call cline_rec%delete('refs_even')
            call cline_rec%delete('refs_odd')
            do state = 1, params%nstates
                call cline_rec%delete('vol'//int2str(state))
            enddo
            if( l_regularized )then
                call cline_rec%set('objfun', 'euclid')
                call cline_rec%set('ml_reg',    'yes')
            else
                call cline_rec%set('objfun',       'cc')
                call cline_rec%set('ml_reg',       'no')
                call cline_rec%set('postprocess',  'no')
                call cline_rec%set('filt_mode',    'none')
                call cline_rec%set('automsk',      'no')
            endif
        end subroutine prepare_bootstrap_rec_cline


        subroutine register_bootstrap_rec_outputs()
            type(sp_project) :: spproj
            type(string)     :: volname, fscname
            integer          :: state, pop
            character(len=16) :: imgkind
            call spproj%read_segment('out', params%projfile)
            call spproj%read_segment(params%oritype, params%projfile)
            select case(trim(params%oritype))
                case('cls3D')
                    imgkind = 'vol_cavg'
                case DEFAULT
                    imgkind = 'vol'
            end select
            do state = 1, params%nstates
                select case(trim(params%oritype))
                    case('cls3D')
                        pop = spproj%os_cls3D%get_pop(state, 'state')
                    case DEFAULT
                        pop = spproj%os_ptcl3D%get_pop(state, 'state')
                end select
                if( pop == 0 )cycle
                volname = refine3D_state_vol_fname(state)
                if( .not. file_exists(volname) )then
                    call volname%kill
                    cycle
                endif
                fscname = refine3D_fsc_fname(state)
                ! params%box_crop/smpd_crop are the effective reconstruction
                ! sampling resolved by params%new or explicitly pinned by
                ! callers that must avoid staged downsampling leakage.
                call spproj%add_vol2os_out(volname, params%smpd_crop, state, trim(imgkind), pop=pop)
                if( file_exists(fscname) ) call spproj%add_fsc2os_out(fscname, state, params%box_crop)
                call volname%kill
                call fscname%kill
            enddo
            call spproj%write_segment_inside('out', params%projfile)
            call spproj%kill
        end subroutine register_bootstrap_rec_outputs


        subroutine warn_for_forced_bootstrap_overrides( cline_in )
            class(cmdline), intent(inout) :: cline_in
            type(string) :: val
            if( cline_in%defined('sigma_est') )then
                val = cline_in%get_carg('sigma_est')
                if( val%to_char().ne.'global' )then
                    THROW_WARN('bootstrap_rec3D enforces sigma_est=global; ignoring input sigma_est='//val%to_char())
                endif
                call val%kill
            endif
            if( cline_in%defined('objfun') ) THROW_WARN('bootstrap_rec3D controls objfun internally; ignoring input objfun')
            if( cline_in%defined('ml_reg') ) THROW_WARN('bootstrap_rec3D controls ml_reg internally; ignoring input ml_reg')
        end subroutine warn_for_forced_bootstrap_overrides
    end subroutine exec_bootstrap_rec3D

    subroutine exec_rec3D_distr_worker( self, cline )
        use simple_rec3D_pcg_strategy, only: execute_rec3D_pcg_worker
        class(commander_rec3D_worker), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters)     :: params
        type(builder)        :: build
        integer, allocatable :: pinds(:)
        integer              :: nptcls2update
        logical              :: l_sigma_loaded
        call build%init_params_and_build_general_tbox(cline, params)
        call build%build_strategy3D_tbox(params)
        if( params%l_update_frac .and. build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls2update, pinds)
        else
            ! we sample all state > 0 and updatecnt > 0
            call build%spproj_field%sample4rec([params%fromp,params%top], nptcls2update, pinds)
        endif
        if( trim(params%rec_backend) == 'pcg' )then
            call execute_rec3D_pcg_worker(params, build, cline, pinds)
        else
            if( params%cc_objfun == OBJFUN_EUCLID )then
                call load_sigma2_groups(params, build%pftc, build%esig, build%spproj, build%spproj_field, &
                    &cline, l_sigma_loaded)
                if( .not. l_sigma_loaded ) THROW_HARD('gridding objfun=euclid requires sigma2 files')
            endif
            if( trim(params%projrec) == 'yes' )then
                call calc_projdir3Drec(params, build, cline, nptcls2update, pinds)
            else
                call calc_3Drec(params, build, cline, nptcls2update, pinds)
            endif
        endif
        ! cleanup
        call build%esig%kill
        call build%kill_strategy3D_tbox
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_rec :: exec_rec3D'))
    end subroutine exec_rec3D_distr_worker

    subroutine exec_random_rec( self, cline )
        class(random_rec_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(commander_rec3D) :: xrec3D
        type(parameters)      :: params
        type(builder)         :: build
        call cline%set('oritype', 'ptcl3D')
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'yes'   )
        call build%init_params_and_build_spproj(cline, params)
        call build%spproj%os_ptcl3D%rnd_oris
        call build%spproj%write_segment_inside('ptcl3D', params%projfile)
        call cline%set('mkdir', 'no') ! to avoid nested dirs
        call cline%set('prg',   'rec3D')
        call xrec3D%execute(cline)
        call build%spproj_field%kill
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_RANDOM_REC NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_random_rec

end module simple_commanders_rec
