module simple_commanders_prob
use simple_commanders_api
use simple_pftc_srch_api
implicit none
#include "simple_local_flags.inc"

type :: prob_bench_state
    integer(timer_int_kind) :: t_tot        = 0
    integer(timer_int_kind) :: t_phase      = 0
    real(timer_int_kind)    :: rt_init      = 0.
    real(timer_int_kind)    :: rt_sample    = 0.
    real(timer_int_kind)    :: rt_refprep   = 0.
    real(timer_int_kind)    :: rt_sigmas    = 0.
    real(timer_int_kind)    :: rt_alloc     = 0.
    real(timer_int_kind)    :: rt_memo_refs = 0.
    real(timer_int_kind)    :: rt_ptcls     = 0.
    real(timer_int_kind)    :: rt_tab_init  = 0.
    real(timer_int_kind)    :: rt_fill_tab  = 0.
    real(timer_int_kind)    :: rt_tab_io    = 0.
    real(timer_int_kind)    :: rt_assign    = 0.
    real(timer_int_kind)    :: rt_cleanup   = 0.
    real(timer_int_kind)    :: rt_tot       = 0.
    real(timer_int_kind)    :: rt_assign_normalize = 0.
    real(timer_int_kind)    :: rt_assign_sort      = 0.
    real(timer_int_kind)    :: rt_assign_graph     = 0.
    real(timer_int_kind)    :: rt_assign_loop      = 0.
    real(timer_int_kind)    :: rt_assign_fallback  = 0.
    real(timer_int_kind)    :: rt_fill_shift_seed   = 0.
    real(timer_int_kind)    :: rt_fill_ref_sweep    = 0.
    real(timer_int_kind)    :: rt_fill_select       = 0.
    real(timer_int_kind)    :: rt_fill_neigh_eval   = 0.
    real(timer_int_kind)    :: rt_fill_shift_refine = 0.
    integer                 :: nptcls       = 0
    integer                 :: nrefs        = 0
    integer                 :: nstates      = 0
    integer                 :: nspace       = 0
    integer                 :: nrots        = 0
    integer                 :: nparts       = 0
    integer                 :: part         = 0
    integer                 :: nthr         = 0
    integer                 :: fill_projs_ns            = 0
    integer                 :: fill_ref_evals           = 0
    integer                 :: fill_nsubs               = 0
    integer                 :: fill_neigh_evals         = 0
    integer                 :: fill_shift_seed_trials   = 0
    integer                 :: fill_shift_refine_trials = 0
    logical                 :: l_doshift    = .false.
end type prob_bench_state

type, extends(commander_base) :: commander_prob_tab
  contains
    procedure :: execute      => exec_prob_tab
end type commander_prob_tab

type, extends(commander_base) :: commander_prob_tab_neigh
    contains
        procedure :: execute      => exec_prob_tab_neigh
end type commander_prob_tab_neigh

type, extends(commander_base) :: commander_prob_align
  contains
    procedure :: execute      => exec_prob_align
end type commander_prob_align

type, extends(commander_base) :: commander_prob_align_neigh
    contains
        procedure :: execute      => exec_prob_align_neigh
end type commander_prob_align_neigh

type, extends(commander_base) :: commander_prob_tab2D
  contains
    procedure :: execute      => exec_prob_tab2D
end type commander_prob_tab2D

type, extends(commander_base) :: commander_prob_align2D
  contains
    procedure :: execute      => exec_prob_align2D
end type commander_prob_align2D

contains

    subroutine prob_bench_start( bench )
        type(prob_bench_state), intent(inout) :: bench
        if( .not. L_BENCH_GLOB ) return
        bench = prob_bench_state()
        bench%t_tot   = tic()
        bench%t_phase = bench%t_tot
    end subroutine prob_bench_start

    subroutine prob_bench_next( bench, elapsed )
        type(prob_bench_state), intent(inout) :: bench
        real(timer_int_kind),   intent(out)   :: elapsed
        elapsed       = 0.
        if( .not. L_BENCH_GLOB ) return
        elapsed       = toc(bench%t_phase)
        bench%t_phase = tic()
    end subroutine prob_bench_next

    subroutine prob_bench_finish( bench )
        type(prob_bench_state), intent(inout) :: bench
        if( .not. L_BENCH_GLOB ) return
        bench%rt_tot = toc(bench%t_tot)
    end subroutine prob_bench_finish

    subroutine prob_bench_capture_context( bench, params, nptcls, nrefs, nrots )
        type(prob_bench_state), intent(inout) :: bench
        type(parameters),       intent(in)    :: params
        integer,                intent(in)    :: nptcls, nrefs, nrots
        if( .not. L_BENCH_GLOB ) return
        bench%nptcls    = nptcls
        bench%nrefs     = nrefs
        bench%nstates   = params%nstates
        bench%nspace    = params%nspace
        bench%nrots     = nrots
        bench%nparts    = params%nparts
        bench%part      = params%part
        bench%nthr      = params%nthr
        bench%l_doshift = params%l_doshift
    end subroutine prob_bench_capture_context

    subroutine prob_bench_capture_assign( bench, rt_normalize, rt_sort, rt_graph, rt_loop, rt_fallback )
        type(prob_bench_state), intent(inout) :: bench
        real(timer_int_kind),   intent(in)    :: rt_normalize, rt_sort, rt_graph, rt_loop, rt_fallback
        if( .not. L_BENCH_GLOB ) return
        bench%rt_assign_normalize = rt_normalize
        bench%rt_assign_sort      = rt_sort
        bench%rt_assign_graph     = rt_graph
        bench%rt_assign_loop      = rt_loop
        bench%rt_assign_fallback  = rt_fallback
    end subroutine prob_bench_capture_assign

    subroutine write_prob_bench_report( params, bench, bench_kind, refine_mode, part_suffix )
        type(parameters),       intent(in) :: params
        type(prob_bench_state), intent(in) :: bench
        character(len=*),       intent(in) :: bench_kind, refine_mode
        logical, optional,      intent(in) :: part_suffix
        type(string) :: benchfname
        integer :: fnr
        logical :: l_part_suffix
        if( .not. L_BENCH_GLOB ) return
        l_part_suffix = .false.
        if( present(part_suffix) ) l_part_suffix = part_suffix
        if( l_part_suffix .and. params%part /= 1 ) return
        benchfname = string(uppercase(trim(bench_kind)))//'_BENCH_ITER'//int2str_pad(params%which_iter,3)
        benchfname = benchfname//TXT_EXT
        call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
        write(fnr,'(a)')    '*** BENCHMARK CONTEXT ***'
        write(fnr,'(a,a)')  'prob refine mode                    : ', trim(refine_mode)
        write(fnr,'(a,a)')  'prob benchmark kind                 : ', trim(bench_kind)
        write(fnr,'(a,i0)') 'prob iteration                      : ', params%which_iter
        write(fnr,'(a,i0)') 'prob part                           : ', bench%part
        write(fnr,'(a,i0)') 'prob nparts                         : ', bench%nparts
        write(fnr,'(a,i0)') 'prob nthr                           : ', bench%nthr
        write(fnr,'(a,i0)') 'prob nptcls sampled                 : ', bench%nptcls
        write(fnr,'(a,i0)') 'prob nrefs active                   : ', bench%nrefs
        write(fnr,'(a,i0)') 'prob nstates                        : ', bench%nstates
        write(fnr,'(a,i0)') 'prob nspace                         : ', bench%nspace
        write(fnr,'(a,i0)') 'prob nrots                          : ', bench%nrots
        write(fnr,'(a,l1)') 'prob l_doshift                      : ', bench%l_doshift
        write(fnr,'(a)') ''
        write(fnr,'(a)') '*** TIMINGS (s) ***'
        if( l_part_suffix )then
            write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' reference preparation      : ', bench%rt_refprep
            write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' particle allocation        : ', bench%rt_alloc
            write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' particle polar prep        : ', bench%rt_ptcls
            write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' table init                 : ', bench%rt_tab_init
            write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' table fill                 : ', bench%rt_fill_tab
            if( (bench%rt_fill_shift_seed + bench%rt_fill_ref_sweep + bench%rt_fill_select + &
                &bench%rt_fill_neigh_eval + bench%rt_fill_shift_refine) > 0. )then
                write(fnr,'(a)') ''
                write(fnr,'(a)') '*** TABLE FILL DETAILS (thread-s) ***'
                if( index(trim(bench_kind),'neigh') > 0 )then
                    write(fnr,'(a,i0)')       trim(bench_kind)//' fill nsubs                 : ', bench%fill_nsubs
                    write(fnr,'(a,i0)')       trim(bench_kind)//' fill peak subspaces/state  : ', bench%fill_projs_ns
                    write(fnr,'(a,i0)')       trim(bench_kind)//' fill coarse evals          : ', bench%fill_ref_evals
                    write(fnr,'(a,i0)')       trim(bench_kind)//' fill neighborhood evals    : ', bench%fill_neigh_evals
                    write(fnr,'(a,i0)')       trim(bench_kind)//' fill shift seed trials     : ', bench%fill_shift_seed_trials
                    write(fnr,'(a,i0)')       trim(bench_kind)//' fill shift refine trials   : ', bench%fill_shift_refine_trials
                    write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' fill shift seed search     : ', bench%rt_fill_shift_seed
                    write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' fill coarse subspace sweep : ', bench%rt_fill_ref_sweep
                    write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' fill neighborhood select   : ', bench%rt_fill_select
                    write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' fill neighborhood eval     : ', bench%rt_fill_neigh_eval
                    write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' fill shift refine          : ', bench%rt_fill_shift_refine
                else
                    write(fnr,'(a,i0)')       trim(bench_kind)//' fill projs_ns              : ', bench%fill_projs_ns
                    write(fnr,'(a,i0)')       trim(bench_kind)//' fill ref evals             : ', bench%fill_ref_evals
                    write(fnr,'(a,i0)')       trim(bench_kind)//' fill shift seed trials     : ', bench%fill_shift_seed_trials
                    write(fnr,'(a,i0)')       trim(bench_kind)//' fill shift refine trials   : ', bench%fill_shift_refine_trials
                    write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' fill shift seed search     : ', bench%rt_fill_shift_seed
                    write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' fill all-reference sweep   : ', bench%rt_fill_ref_sweep
                    write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' fill neighborhood select   : ', bench%rt_fill_select
                    write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' fill shift refine          : ', bench%rt_fill_shift_refine
                endif
                write(fnr,'(a)') ''
                write(fnr,'(a)') '*** TIMINGS (s), continued ***'
            endif
            write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' table write                : ', bench%rt_tab_io
        else
            write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' table generation           : ', bench%rt_fill_tab
            write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' table read/write           : ', bench%rt_tab_io
            write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' assignment                 : ', bench%rt_assign
            if( (bench%rt_assign_normalize + bench%rt_assign_sort + bench%rt_assign_graph + &
                &bench%rt_assign_loop + bench%rt_assign_fallback) > 0. )then
                write(fnr,'(a)') ''
                write(fnr,'(a)') '*** ASSIGNMENT DETAILS (s) ***'
                write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' assign normalize           : ', bench%rt_assign_normalize
                write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' assign sort                : ', bench%rt_assign_sort
                write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' assign graph/setup         : ', bench%rt_assign_graph
                write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' assign loop                : ', bench%rt_assign_loop
                write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' assign fallback            : ', bench%rt_assign_fallback
                write(fnr,'(a)') ''
                write(fnr,'(a)') '*** TIMINGS (s), continued ***'
            endif
        endif
        write(fnr,'(a,t52,f9.2)') trim(bench_kind)//' total time                 : ', bench%rt_tot
        call fclose(fnr)
        call benchfname%kill
    end subroutine write_prob_bench_report

    subroutine exec_prob_tab( self, cline )
        use simple_matcher_2Dprep
        use simple_matcher_refvol_utils,    only: read_reprojection_model
        use simple_matcher_ptcl_batch,      only: prep_sigmas_objfun, alloc_ptcl_imgs, build_batch_particles3D, clean_batch_particles3D
        use simple_eul_prob_tab,            only: eul_prob_tab
        class(commander_prob_tab), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        integer,     allocatable :: pinds(:)
        type(image), allocatable :: tmp_imgs(:), tmp_imgs_pad(:)
        type(string)             :: fname
        type(builder)            :: build
        type(parameters)         :: params
        type(eul_prob_tab)       :: eulprob_obj_part
        integer :: nptcls
        type(prob_bench_state) :: bench
        call prob_bench_start(bench)
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.true.)
        call prob_bench_next(bench, bench%rt_init)
        ! The policy here ought to be that nothing is done with regards to sampling other than reproducing
        ! what was generated in the driver (prob_align, below). Sampling is delegated to prob_align (below)
        ! and merely reproduced here
        if( build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls, pinds)
        else
            THROW_HARD('exec_prob_tab requires prior particle sampling (in exec_prob_align)')
        endif
        call prob_bench_next(bench, bench%rt_sample)
        ! PREPARE REFERENCES, SIGMAS, POLAR_CORRCALC, PTCLS
        call read_reprojection_model(params, build, nptcls, nmany_refs=params%nspace*params%nstates)
        call prob_bench_next(bench, bench%rt_refprep)
        call prep_sigmas_objfun(params, build, .false.)
        call prob_bench_next(bench, bench%rt_sigmas)
        call alloc_ptcl_imgs( params, build, tmp_imgs, tmp_imgs_pad, nptcls )
        call prob_bench_next(bench, bench%rt_alloc)
        call build%pftc%memoize_refs(eulspace=build%eulspace)
        call prob_bench_next(bench, bench%rt_memo_refs)
        ! Build polar particle images
        call build_batch_particles3D(params, build, nptcls, pinds, tmp_imgs, tmp_imgs_pad)
        call prob_bench_next(bench, bench%rt_ptcls)
        ! Filling prob table in eul_prob_tab
        call eulprob_obj_part%new(params, build, pinds)
        call prob_bench_capture_context(bench, params, nptcls, eulprob_obj_part%nrefs, build%pftc%get_nrots())
        call prob_bench_next(bench, bench%rt_tab_init)
        fname = string(DIST_FBODY)//int2str_pad(params%part,params%numlen)//'.dat'
        if( str_has_substr(params%refine, 'prob_state') )then
            call eulprob_obj_part%fill_tab_state_only
            call prob_bench_next(bench, bench%rt_fill_tab)
            call eulprob_obj_part%write_state_tab(fname)
        else
            call eulprob_obj_part%fill_tab
            bench%fill_projs_ns            = eulprob_obj_part%bench_fill_projs_ns
            bench%fill_ref_evals           = eulprob_obj_part%bench_fill_ref_evals
            bench%fill_shift_seed_trials   = eulprob_obj_part%bench_fill_shift_seed_trials
            bench%fill_shift_refine_trials = eulprob_obj_part%bench_fill_shift_refine_trials
            bench%rt_fill_shift_seed       = eulprob_obj_part%bench_fill_shift_seed
            bench%rt_fill_ref_sweep        = eulprob_obj_part%bench_fill_ref_sweep
            bench%rt_fill_select           = eulprob_obj_part%bench_fill_select
            bench%rt_fill_shift_refine     = eulprob_obj_part%bench_fill_shift_refine
            call prob_bench_next(bench, bench%rt_fill_tab)
            call eulprob_obj_part%write_tab(fname)
        endif
        call prob_bench_next(bench, bench%rt_tab_io)
        call eulprob_obj_part%kill
        call clean_batch_particles3D(build, tmp_imgs, tmp_imgs_pad)
        call build%pftc%kill
        call build%kill_general_tbox
        call prob_bench_next(bench, bench%rt_cleanup)
        call prob_bench_finish(bench)
        call write_prob_bench_report(params, bench, 'prob_tab', params%refine, part_suffix=.true.)
        call qsys_job_finished(params, string('simple_commanders_refine3D :: exec_prob_tab'))
        call simple_end('**** SIMPLE_PROB_TAB NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_tab

    subroutine exec_prob_tab_neigh( self, cline )
        use simple_matcher_2Dprep
        use simple_matcher_refvol_utils,    only: read_reprojection_model
        use simple_matcher_ptcl_batch,      only: prep_sigmas_objfun, alloc_ptcl_imgs, build_batch_particles3D, clean_batch_particles3D
        use simple_eul_prob_tab_neigh,      only: eul_prob_tab_neigh
        class(commander_prob_tab_neigh), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        integer,     allocatable :: pinds(:)
        type(image), allocatable :: tmp_imgs(:), tmp_imgs_pad(:)
        type(string)             :: fname
        type(builder)            :: build
        type(parameters)         :: params
        type(eul_prob_tab_neigh) :: eulprob_obj_part_neigh
        integer :: nptcls
        type(prob_bench_state) :: bench
        call prob_bench_start(bench)
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.true.)
        call prob_bench_next(bench, bench%rt_init)
        ! Sampling policy mirrors exec_prob_tab: only reproduce already sampled particles.
        if( build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls, pinds)
        else
            THROW_HARD('exec_prob_tab_neigh requires prior particle sampling (in exec_prob_align)')
        endif
        call prob_bench_next(bench, bench%rt_sample)
        ! PREPARE REFERENCES, SIGMAS, POLAR_CORRCALC, PTCLS
        call read_reprojection_model(params, build, nptcls)
        call prob_bench_next(bench, bench%rt_refprep)
        call prep_sigmas_objfun(params, build, .false.)
        call prob_bench_next(bench, bench%rt_sigmas)
        call alloc_ptcl_imgs( params, build, tmp_imgs, tmp_imgs_pad, nptcls )
        call prob_bench_next(bench, bench%rt_alloc)
        call build%pftc%memoize_refs(eulspace=build%eulspace)
        call prob_bench_next(bench, bench%rt_memo_refs)
        ! Build polar particle images
        call build_batch_particles3D(params, build, nptcls, pinds, tmp_imgs, tmp_imgs_pad)
        call prob_bench_next(bench, bench%rt_ptcls)
        call eulprob_obj_part_neigh%new_neigh(params, build, pinds)
        call prob_bench_capture_context(bench, params, nptcls, eulprob_obj_part_neigh%nrefs, build%pftc%get_nrots())
        call prob_bench_next(bench, bench%rt_tab_init)
        call eulprob_obj_part_neigh%fill_tab
        bench%fill_projs_ns            = eulprob_obj_part_neigh%bench_fill_projs_ns
        bench%fill_ref_evals           = eulprob_obj_part_neigh%bench_fill_ref_evals
        bench%fill_nsubs               = eulprob_obj_part_neigh%bench_fill_nsubs
        bench%fill_neigh_evals         = eulprob_obj_part_neigh%bench_fill_neigh_evals
        bench%fill_shift_seed_trials   = eulprob_obj_part_neigh%bench_fill_shift_seed_trials
        bench%fill_shift_refine_trials = eulprob_obj_part_neigh%bench_fill_shift_refine_trials
        bench%rt_fill_shift_seed       = eulprob_obj_part_neigh%bench_fill_shift_seed
        bench%rt_fill_ref_sweep        = eulprob_obj_part_neigh%bench_fill_ref_sweep
        bench%rt_fill_select           = eulprob_obj_part_neigh%bench_fill_select
        bench%rt_fill_neigh_eval       = eulprob_obj_part_neigh%bench_fill_neigh_eval
        bench%rt_fill_shift_refine     = eulprob_obj_part_neigh%bench_fill_shift_refine
        call prob_bench_next(bench, bench%rt_fill_tab)
        fname = string(DIST_FBODY)//'_neigh_'//int2str_pad(params%part,params%numlen)//'.dat'
        call eulprob_obj_part_neigh%write_tab(fname)
        call prob_bench_next(bench, bench%rt_tab_io)
        call eulprob_obj_part_neigh%kill
        call clean_batch_particles3D(build, tmp_imgs, tmp_imgs_pad)
        call build%pftc%kill
        call build%kill_general_tbox
        call prob_bench_next(bench, bench%rt_cleanup)
        call prob_bench_finish(bench)
        call write_prob_bench_report(params, bench, 'prob_tab_neigh', params%refine, part_suffix=.true.)
        call qsys_job_finished(params, string('simple_commanders_refine3D :: exec_prob_tab_neigh'))
        call simple_end('**** SIMPLE_PROB_TAB_NEIGH NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_tab_neigh

    subroutine exec_prob_align( self, cline )
        use simple_eul_prob_tab,            only: eul_prob_tab
        use simple_matcher_smpl_and_lplims, only: sample_ptcls4fillin, sample_ptcls4update3D
        use simple_builder,                 only: builder
        class(commander_prob_align), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        integer,     allocatable :: pinds(:)
        type(string)             :: fname
        type(builder)            :: build
        type(parameters)         :: params
        type(commander_prob_tab) :: xprob_tab
        type(eul_prob_tab)       :: eulprob_obj_glob
        type(cmdline)            :: cline_prob_tab
        type(qsys_env)           :: qenv
        type(chash)              :: job_descr
        integer :: nptcls, ipart
        real(timer_int_kind) :: rt_tmp
        type(prob_bench_state) :: bench
        call prob_bench_start(bench)
        call cline%set('mkdir',  'no')
        call cline%set('stream', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        call prob_bench_next(bench, bench%rt_init)
        if( params%startit == 1 ) call build%spproj_field%clean_entry('updatecnt', 'sampled')
        ! sampled incremented
        if( params%l_fillin .and. mod(params%startit,5) == 0 )then
            call sample_ptcls4fillin(params, build, [1,params%nptcls], .true., nptcls, pinds)
        else
            call sample_ptcls4update3D(params, build, [1,params%nptcls], .true., nptcls, pinds)
        endif
        call prob_bench_next(bench, bench%rt_sample)
        ! communicate to project file
        call build%spproj%write_segment_inside(params%oritype)
        ! more prep
        call eulprob_obj_glob%new(params, build, pinds)
        call prob_bench_capture_context(bench, params, nptcls, eulprob_obj_glob%nrefs, build%pftc%get_nrots())
        call prob_bench_next(bench, bench%rt_tab_init)
        ! generating all corrs on all parts
        cline_prob_tab = cline
        call cline_prob_tab%set('prg', 'prob_tab' ) ! required for distributed call
        ! execution
        if( .not.cline_prob_tab%defined('nparts') )then
            call xprob_tab%execute(cline_prob_tab)
        else
            ! setup the environment for distributed execution
            call qenv%new(params, params%nparts, nptcls=params%nptcls)
            call cline_prob_tab%gen_job_descr(job_descr)
            ! schedule
            call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        endif
        call prob_bench_next(bench, bench%rt_fill_tab)
        ! reading corrs from all parts
        if( str_has_substr(params%refine, 'prob_state') )then
            do ipart = 1, params%nparts
                fname = string(DIST_FBODY)//int2str_pad(ipart,params%numlen)//'.dat'
                call eulprob_obj_glob%read_state_tab(fname)
            enddo
            call prob_bench_next(bench, bench%rt_tab_io)
            call eulprob_obj_glob%state_assign
        else
            do ipart = 1, params%nparts
                fname = string(DIST_FBODY)//int2str_pad(ipart,params%numlen)//'.dat'
                call eulprob_obj_glob%read_tab_to_glob(fname)
            enddo
            call prob_bench_next(bench, bench%rt_tab_io)
            call eulprob_obj_glob%ref_assign
            call prob_bench_capture_assign(bench, eulprob_obj_glob%bench_assign_normalize, eulprob_obj_glob%bench_assign_sort,&
                &eulprob_obj_glob%bench_assign_graph, eulprob_obj_glob%bench_assign_loop, eulprob_obj_glob%bench_assign_fallback)
        endif
        call prob_bench_next(bench, bench%rt_assign)
        ! write the iptcl->(iref,istate) assignment
        fname = string(ASSIGNMENT_FBODY)//'.dat'
        call eulprob_obj_glob%write_assignment(fname)
        call prob_bench_next(bench, rt_tmp)
        bench%rt_tab_io = bench%rt_tab_io + rt_tmp
        ! cleanup
        call eulprob_obj_glob%kill
        call cline_prob_tab%kill
        call qenv%kill
        call job_descr%kill
        call build%kill_general_tbox
        call prob_bench_next(bench, bench%rt_cleanup)
        call prob_bench_finish(bench)
        call write_prob_bench_report(params, bench, 'prob_align', params%refine)
        call qsys_job_finished(params, string('simple_commanders_refine3D :: exec_prob_align'))
        call qsys_cleanup(params)
        call simple_end('**** SIMPLE_PROB_ALIGN NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_align

    subroutine exec_prob_align_neigh( self, cline )
        use simple_eul_prob_tab_neigh,      only: eul_prob_tab_neigh
        use simple_matcher_smpl_and_lplims, only: sample_ptcls4fillin, sample_ptcls4update3D
        use simple_builder,                 only: builder
        class(commander_prob_align_neigh), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        integer,           allocatable :: pinds(:)
        type(string)                   :: fname
        type(builder)                  :: build
        type(parameters)               :: params
        type(commander_prob_tab_neigh) :: xprob_tab_neigh
        type(eul_prob_tab_neigh)       :: eulprob_obj_glob_neigh
        type(cmdline)                  :: cline_prob_tab
        type(qsys_env)                 :: qenv
        type(chash)                    :: job_descr
        integer :: nptcls
        real(timer_int_kind) :: rt_tmp
        type(prob_bench_state) :: bench
        call prob_bench_start(bench)
        call cline%set('mkdir',  'no')
        call cline%set('stream', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        call prob_bench_next(bench, bench%rt_init)
        if( params%startit == 1 ) call build%spproj_field%clean_entry('updatecnt', 'sampled')
        if( params%l_fillin .and. mod(params%startit,5) == 0 )then
            call sample_ptcls4fillin(params, build, [1,params%nptcls], .true., nptcls, pinds)
        else
            call sample_ptcls4update3D(params, build, [1,params%nptcls], .true., nptcls, pinds)
        endif
        call prob_bench_next(bench, bench%rt_sample)
        call build%spproj%write_segment_inside(params%oritype)
        ! Global object only needs sampled-set maps before reading partition sparse tables.
        ! The neighborhood scoring itself is performed in each prob_tab_neigh partition job.
        call eulprob_obj_glob_neigh%new_neigh(params, build, pinds)
        call prob_bench_capture_context(bench, params, nptcls, eulprob_obj_glob_neigh%nrefs, build%pftc%get_nrots())
        call prob_bench_next(bench, bench%rt_tab_init)
        cline_prob_tab = cline
        call cline_prob_tab%set('prg', 'prob_tab_neigh')
        if( .not. cline_prob_tab%defined('nparts') )then
            call xprob_tab_neigh%execute(cline_prob_tab)
        else
            call qenv%new(params, params%nparts, nptcls=params%nptcls)
            call cline_prob_tab%gen_job_descr(job_descr)
            call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        endif
        call prob_bench_next(bench, bench%rt_fill_tab)
        call eulprob_obj_glob_neigh%read_tabs_to_glob(string(DIST_FBODY)//'_neigh_', params%nparts, params%numlen)
        call prob_bench_next(bench, bench%rt_tab_io)
        call eulprob_obj_glob_neigh%ref_assign
        call prob_bench_capture_assign(bench, eulprob_obj_glob_neigh%bench_assign_normalize, eulprob_obj_glob_neigh%bench_assign_sort,&
            &eulprob_obj_glob_neigh%bench_assign_graph, eulprob_obj_glob_neigh%bench_assign_loop, eulprob_obj_glob_neigh%bench_assign_fallback)
        call prob_bench_next(bench, bench%rt_assign)
        ! write the iptcl->(iref,istate) assignment
        fname = string(ASSIGNMENT_FBODY)//'.dat'
        call eulprob_obj_glob_neigh%write_assignment(fname)
        call prob_bench_next(bench, rt_tmp)
        bench%rt_tab_io = bench%rt_tab_io + rt_tmp
        ! cleanup
        call eulprob_obj_glob_neigh%kill
        call cline_prob_tab%kill
        call qenv%kill
        call job_descr%kill
        call build%kill_general_tbox
        call prob_bench_next(bench, bench%rt_cleanup)
        call prob_bench_finish(bench)
        call write_prob_bench_report(params, bench, 'prob_align_neigh', params%refine)
        call qsys_job_finished(params, string('simple_commanders_refine3D :: exec_prob_align_neigh'))
        call qsys_cleanup(params)
        call simple_end('**** SIMPLE_PROB_ALIGN_NEIGH NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_align_neigh

    subroutine exec_prob_tab2D( self, cline )
        use simple_matcher_smpl_and_lplims, only: set_bp_range2D
        use simple_strategy2D_matcher,  only: set_b_p_ptrs2D, &
                                              ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad
        use simple_matcher_pftc_prep,      only: prep_pftc4align2D
        use simple_matcher_ptcl_batch,  only: alloc_ptcl_imgs, build_batch_particles2D, clean_batch_particles2D
        use simple_imgarr_utils,        only: alloc_imgarr
        use simple_classaverager,       only: cavger_new, cavger_read_all, cavger_kill
        use simple_eul_prob_tab2D,      only: eul_prob_tab2D
        class(commander_prob_tab2D), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        integer,     allocatable :: pinds(:)
        type(string)             :: fname
        type(builder)            :: build
        type(parameters)         :: params
        type(eul_prob_tab2D)     :: eulprob_obj_part
        real    :: frac_srch_space
        integer :: nptcls
        logical :: l_alloc_read_cavgs
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        if( build%spproj_field%get_nevenodd() == 0 )then
            call build%spproj_field%partition_eo
            call build%spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        frac_srch_space  = build%spproj_field%get_avg('frac')
        call set_bp_range2D(params, build, cline, params%which_iter, frac_srch_space)
        ! reproduce particle sampling from exec_prob_align2D
        if( build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls, pinds)
        else
            THROW_HARD('exec_prob_tab2D requires prior particle sampling (in exec_prob_align2D)')
        endif
        call set_b_p_ptrs2D(params, build)
        call alloc_ptcl_imgs(params, build, ptcl_match_imgs, ptcl_match_imgs_pad, nptcls)
        call alloc_imgarr(nptcls, [params%box, params%box, 1], params%smpd, ptcl_imgs)
        ! mirror cluster2D_exec reference setup
        l_alloc_read_cavgs = l_distr_worker_glob .or. (params%which_iter==1)
        call cavger_new(params, build, alloccavgs=l_alloc_read_cavgs)
        if( .not. cline%defined('refs') ) THROW_HARD('exec_prob_tab2D requires refs on the command line')
        call cavger_read_all
        call prep_pftc4align2D(params, build, ptcl_match_imgs_pad, nptcls, params%which_iter, .false.)
        ! build polar particle images
        call build_batch_particles2D(params, build, nptcls, pinds, ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad)
        ! fill and write the 2D probability table
        call eulprob_obj_part%new(params, build, pinds)
        call eulprob_obj_part%fill_tab
        fname = string(DIST_FBODY)//int2str_pad(params%part,params%numlen)//'.dat'
        call eulprob_obj_part%write_tab(fname)
        call eulprob_obj_part%kill
        call clean_batch_particles2D(build, ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad)
        if( l_distr_worker_glob ) call cavger_kill
        call build%pftc%kill
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_prob :: exec_prob_tab2D'))
        call simple_end('**** SIMPLE_PROB_TAB2D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_tab2D

    subroutine exec_prob_align2D( self, cline )
        use simple_eul_prob_tab2D,          only: eul_prob_tab2D
        use simple_strategy2D_matcher,      only: set_b_p_ptrs2D
        use simple_matcher_smpl_and_lplims, only: sample_ptcls4update2D
        use simple_builder,                 only: builder
        class(commander_prob_align2D), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        integer,       allocatable :: pinds(:)
        type(string)               :: fname
        type(builder)              :: build
        type(parameters)           :: params
        type(commander_prob_tab2D) :: xprob_tab2D
        type(eul_prob_tab2D)       :: eulprob_obj_glob
        type(cmdline)              :: cline_prob_tab
        type(qsys_env)             :: qenv
        type(chash)                :: job_descr
        integer :: nptcls, ipart
        real(timer_int_kind) :: rt_tmp
        type(prob_bench_state) :: bench
        call prob_bench_start(bench)
        call cline%set('mkdir',  'no')
        call cline%set('stream', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        call set_b_p_ptrs2D(params, build)
        call prob_bench_next(bench, bench%rt_init)
        if( build%spproj_field%get_nevenodd() == 0 )then
            call build%spproj_field%partition_eo
            call build%spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        if( params%startit == 1 .and. params%which_iter == params%startit )then
            call build%spproj_field%clean_entry('updatecnt', 'sampled')
        endif
        ! Mirror the 3D workflow: sampled-update is active from the first stage onward.
        ! In probabilistic mode the sampled subset is reused within the current iteration
        ! by prob_tab2D/cluster2D_exec, but it is redrawn on later iterations.
        call sample_ptcls4update2D(params, build, [params%fromp,params%top], params%l_update_frac, nptcls, pinds)
        call prob_bench_next(bench, bench%rt_sample)
        ! write sampling to project
        call build%spproj%write_segment_inside(params%oritype)
        ! build the global prob table (nclasses x nptcls)
        call eulprob_obj_glob%new(params, build, pinds)
        call prob_bench_capture_context(bench, params, nptcls, eulprob_obj_glob%nclasses, 0)
        call prob_bench_next(bench, bench%rt_tab_init)
        ! generate partition-wise dist tables
        cline_prob_tab = cline
        call cline_prob_tab%set('prg', 'prob_tab2D')
        if( .not. cline_prob_tab%defined('nparts') )then
            call xprob_tab2D%execute(cline_prob_tab)
        else
            call qenv%new(params, params%nparts, nptcls=params%nptcls)
            call cline_prob_tab%gen_job_descr(job_descr)
            call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        endif
        call prob_bench_next(bench, bench%rt_fill_tab)
        ! merge all partition tables into global
        do ipart = 1, params%nparts
            fname = string(DIST_FBODY)//int2str_pad(ipart,params%numlen)//'.dat'
            call eulprob_obj_glob%read_tab_to_glob(fname)
        end do
        call prob_bench_next(bench, bench%rt_tab_io)
        ! global probabilistic class assignment
        call eulprob_obj_glob%ref_assign
        call prob_bench_capture_assign(bench, eulprob_obj_glob%bench_assign_normalize, eulprob_obj_glob%bench_assign_sort,&
            &eulprob_obj_glob%bench_assign_graph, eulprob_obj_glob%bench_assign_loop, eulprob_obj_glob%bench_assign_fallback)
        call prob_bench_next(bench, bench%rt_assign)
        ! write assignment to file
        fname = string(ASSIGNMENT_FBODY)//'.dat'
        call eulprob_obj_glob%write_assignment(fname)
        call prob_bench_next(bench, rt_tmp)
        bench%rt_tab_io = bench%rt_tab_io + rt_tmp
        ! cleanup
        call eulprob_obj_glob%kill
        call cline_prob_tab%kill
        call qenv%kill
        call job_descr%kill
        call build%kill_general_tbox
        call prob_bench_next(bench, bench%rt_cleanup)
        call prob_bench_finish(bench)
        call write_prob_bench_report(params, bench, 'prob_align2D', params%refine)
        call qsys_job_finished(params, string('simple_commanders_prob :: exec_prob_align2D'))
        call qsys_cleanup(params)
        call simple_end('**** SIMPLE_PROB_ALIGN2D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_align2D

end module simple_commanders_prob
