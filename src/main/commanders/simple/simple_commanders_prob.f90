!@descr: probability-table commanders of the 2D and 3D searches (prob_tab, prob_align and their neighbourhood variants)
module simple_commanders_prob
use simple_commanders_api
use simple_pftc_srch_api
use simple_eo_partition, only: ensure_field_halves
implicit none
private :: write_empty_prob_tab, write_prob_tab_bench, write_prob_align_bench, run_prob_align, run_prob_align_neigh
#include "simple_local_flags.inc"

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
    procedure :: execute              => exec_prob_align
    procedure :: execute_from_project => exec_prob_align_from_project
end type commander_prob_align

type, extends(commander_base) :: commander_prob_align_neigh
    contains
        procedure :: execute              => exec_prob_align_neigh
        procedure :: execute_from_project => exec_prob_align_neigh_from_project
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

    subroutine exec_prob_tab( self, cline )
        use simple_matcher_2Dprep
        use simple_matcher_refvol_utils, only: read_reprojection_model
        use simple_matcher_ptcl_batch,   only: prep_sigmas_objfun, alloc_ptcl_imgs, build_batch_particles3D, clean_batch_particles3D
        use simple_matcher_ptcl_io,      only: ptcl_read_counts, ptcl_read_counts_now, ptcl_read_counts_since, &
            &ptcl_read_session
        use simple_eul_prob_tab,         only: eul_prob_tab
        use simple_eul_prob_tab_utils,   only: prob_fill_profile
        class(commander_prob_tab), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        integer,     allocatable :: pinds(:)
        type(image), allocatable :: tmp_imgs(:), tmp_imgs_pad(:)
        type(string)             :: fname
        type(builder)            :: build
        type(parameters)         :: params
        type(eul_prob_tab)       :: eulprob_obj_part
        type(prob_fill_profile)  :: fill_prof
        type(ptcl_read_counts)   :: reads0
        type(ptcl_read_session)  :: reads
        integer :: nptcls, batchsz_max, nbatches, ibatch, batch_start, batch_end, batchsz
        integer, allocatable :: batches(:,:)
        integer(timer_int_kind) :: t_tot, t
        real(timer_int_kind)    :: rt_setup, rt_refs, rt_batch, rt_fill, rt_write
        if( L_BENCH_GLOB )then
            t_tot  = tic()
            reads0 = ptcl_read_counts_now()
        endif
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        ! The policy here ought to be that nothing is done with regards to sampling other than reproducing
        ! what was generated in the driver (prob_align, below). Sampling is delegated to prob_align (below)
        ! and merely reproduced here
        if( build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls, pinds, &
                &allow_empty=params%l_distr_worker)
        else
            THROW_HARD('exec_prob_tab requires particle sampling from exec_prob_align')
        endif
        fname = string(DIST_FBODY)//int2str_pad(params%part,params%numlen)//'.dat'
        if( nptcls < 1 )then
            ! a partition without sampled particles (state-0 rows only, e.g. the
            ! frozen rows of solve3D_addon) contributes an empty table
            call write_empty_prob_tab(params, build, pinds, fname, eulprob_obj_part)
            call build%kill_general_tbox
            call qsys_declare_part_finished(params, string('simple_commanders_refine3D :: exec_prob_tab'))
            call simple_end('**** SIMPLE_PROB_TAB NORMAL STOP ****', print_simple=.false.)
            return
        endif
        batchsz_max = min(nptcls, params%nthr * BATCHTHRSZ)
        nbatches    = ceiling(real(nptcls) / real(batchsz_max))
        batches     = split_nobjs_even(nptcls, nbatches)
        batchsz_max = maxval(batches(:,2) - batches(:,1) + 1)
        if( L_BENCH_GLOB )then
            rt_setup = toc(t_tot)
            t        = tic()
        endif
        ! PREPARE REFERENCES, SIGMAS, POLAR_CORRCALC, PTCLS
        call read_reprojection_model(params, build, batchsz_max)
        call prep_sigmas_objfun(params, build)
        call alloc_ptcl_imgs( params, build, tmp_imgs, tmp_imgs_pad, batchsz_max )
        call build%pftc%memoize_refs(eulspace=build%eulspace)
        ! Fill the partition table in matcher-sized batches to cap particle PFT memo memory.
        call eulprob_obj_part%new_worker(params,build,pinds)
        call eulprob_obj_part%begin_write(fname)
        if( L_BENCH_GLOB )then
            rt_refs  = toc(t)
            rt_batch = 0.
            rt_fill  = 0.
        endif
        call reads%new(params)
        call reads%expect(pinds(1:nptcls))
        do ibatch = 1, nbatches
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            batchsz     = batch_end - batch_start + 1
            if( L_BENCH_GLOB ) t = tic()
            call build_batch_particles3D(params, build, batchsz, pinds(batch_start:batch_end), tmp_imgs, tmp_imgs_pad, &
                &reads)
            if( L_BENCH_GLOB )then
                rt_batch = rt_batch + toc(t)
                t        = tic()
            endif
            call eulprob_obj_part%fill_tab_range(batch_start, batch_end)
            if( L_BENCH_GLOB ) rt_fill = rt_fill + toc(t)
        end do
        call reads%kill
        if( L_BENCH_GLOB ) t = tic()
        call eulprob_obj_part%write_tab(fname)
        if( L_BENCH_GLOB ) rt_write = toc(t)
        fill_prof = eulprob_obj_part%get_fill_profile()
        call eulprob_obj_part%kill
        if( allocated(batches) ) deallocate(batches)
        call build%pftc%kill
        call clean_batch_particles3D(build, tmp_imgs, tmp_imgs_pad)
        if( L_BENCH_GLOB ) call write_prob_tab_bench(params, nptcls, ptcl_read_counts_since(reads0), &
            &rt_setup, rt_refs, rt_batch, rt_fill, rt_write, toc(t_tot), fill=fill_prof)
        call build%kill_general_tbox
        call qsys_declare_part_finished(params, string('simple_commanders_refine3D :: exec_prob_tab'))
        call simple_end('**** SIMPLE_PROB_TAB NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_tab

    subroutine exec_prob_tab_neigh( self, cline )
        use simple_matcher_2Dprep
        use simple_matcher_refvol_utils, only: read_reprojection_model
        use simple_matcher_ptcl_batch,   only: prep_sigmas_objfun, alloc_ptcl_imgs, build_batch_particles3D, clean_batch_particles3D
        use simple_matcher_ptcl_io,      only: ptcl_read_counts, ptcl_read_counts_now, ptcl_read_counts_since, &
            &ptcl_read_session
        use simple_eul_prob_tab_neigh,   only: eul_prob_tab_neigh
        use simple_eul_prob_tab_utils,   only: prob_fill_profile
        class(commander_prob_tab_neigh), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        integer,     allocatable :: pinds(:)
        type(image), allocatable :: tmp_imgs(:), tmp_imgs_pad(:)
        type(string)             :: fname
        type(builder)            :: build
        type(parameters)         :: params
        type(eul_prob_tab_neigh) :: eulprob_obj_part_neigh
        type(prob_fill_profile)  :: fill_prof
        type(ptcl_read_counts)   :: reads0
        type(ptcl_read_session)  :: reads
        integer :: nptcls
        integer(timer_int_kind) :: t_tot, t
        real(timer_int_kind)    :: rt_setup, rt_refs, rt_batch, rt_fill, rt_write
        if( L_BENCH_GLOB )then
            t_tot  = tic()
            reads0 = ptcl_read_counts_now()
        endif
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.true.)
        ! Sampling policy mirrors exec_prob_tab: only reproduce already sampled particles.
        if( build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls, pinds, &
                &allow_empty=params%l_distr_worker)
        else
            THROW_HARD('exec_prob_tab_neigh requires particle sampling from exec_prob_align')
        endif
        fname = string(DIST_FBODY)//'_neigh_'//int2str_pad(params%part,params%numlen)//'.dat'
        if( nptcls < 1 )then
            ! a partition without sampled particles contributes an empty table
            call write_empty_prob_tab(params, build, pinds, fname, eulprob_obj_part_neigh)
        else
            ! All neighborhood modes can fill the table in matcher-sized batches; this
            ! caps particle-image/PFT memo memory without changing assignment ownership.
            call run_prob_tab_neigh_batch(fname)
        endif
        call fname%kill
        call build%kill_general_tbox
        call qsys_declare_part_finished(params, string('simple_commanders_refine3D :: exec_prob_tab_neigh'))
        call simple_end('**** SIMPLE_PROB_TAB_NEIGH NORMAL STOP ****', print_simple=.false.)

    contains

        subroutine prepare_prob_neigh_workspace(batchsz_here)
            integer, intent(in) :: batchsz_here
            call read_reprojection_model(params, build, batchsz_here)
            call prep_sigmas_objfun(params, build)
            call alloc_ptcl_imgs(params, build, tmp_imgs, tmp_imgs_pad, batchsz_here)
            call build%pftc%memoize_refs(eulspace=build%eulspace)
        end subroutine prepare_prob_neigh_workspace

        subroutine cleanup_prob_neigh_workspace
            call build%pftc%kill
            call clean_batch_particles3D(build, tmp_imgs, tmp_imgs_pad)
        end subroutine cleanup_prob_neigh_workspace

        subroutine run_prob_tab_neigh_batch(outfname)
            class(string), intent(in)  :: outfname
            integer, allocatable :: batches(:,:)
            integer :: ibatch, batch_start, batch_end, batchsz, batchsz_max, nbatches
            batchsz_max = min(nptcls, max(1, params%nthr * BATCHTHRSZ))
            nbatches    = ceiling(real(nptcls) / real(batchsz_max))
            batches     = split_nobjs_even(nptcls, nbatches)
            batchsz_max = maxval(batches(:,2) - batches(:,1) + 1)
            if( L_BENCH_GLOB )then
                rt_setup = toc(t_tot)
                t        = tic()
            endif
            call prepare_prob_neigh_workspace(batchsz_max)
            call eulprob_obj_part_neigh%new_neigh(params,build,pinds)
            call eulprob_obj_part_neigh%begin_write(outfname)
            if( L_BENCH_GLOB )then
                rt_refs  = toc(t)
                rt_batch = 0.
                rt_fill  = 0.
            endif
            call reads%new(params)
            call reads%expect(pinds(1:nptcls))
            do ibatch = 1, nbatches
                batch_start = batches(ibatch,1)
                batch_end   = batches(ibatch,2)
                batchsz     = batch_end - batch_start + 1
                if( L_BENCH_GLOB ) t = tic()
                call build_batch_particles3D(params, build, batchsz, pinds(batch_start:batch_end), tmp_imgs, tmp_imgs_pad, &
                &reads)
                if( L_BENCH_GLOB )then
                    rt_batch = rt_batch + toc(t)
                    t        = tic()
                endif
                call eulprob_obj_part_neigh%fill_tab_range(batch_start, batch_end)
                if( L_BENCH_GLOB ) rt_fill = rt_fill + toc(t)
            enddo
            call reads%kill
            if( L_BENCH_GLOB ) t = tic()
            call eulprob_obj_part_neigh%write_tab(outfname)
            if( L_BENCH_GLOB ) rt_write = toc(t)
            fill_prof = eulprob_obj_part_neigh%get_fill_profile()
            call eulprob_obj_part_neigh%kill
            if( allocated(batches) ) deallocate(batches)
            call cleanup_prob_neigh_workspace
            if( L_BENCH_GLOB ) call write_prob_tab_bench(params, nptcls, ptcl_read_counts_since(reads0), &
                &rt_setup, rt_refs, rt_batch, rt_fill, rt_write, toc(t_tot), fill=fill_prof)
        end subroutine run_prob_tab_neigh_batch

    end subroutine exec_prob_tab_neigh


    subroutine exec_prob_align( self, cline )
        class(commander_prob_align), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        call run_prob_align(cline)
    end subroutine exec_prob_align

    !> prob_align with the caller's project of projfile in memory, read and unchanged since: copied, not read again
    subroutine exec_prob_align_from_project( self, cline, spproj )
        class(commander_prob_align), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        class(sp_project),           intent(in)    :: spproj
        call run_prob_align(cline, spproj)
    end subroutine exec_prob_align_from_project

    subroutine run_prob_align( cline, spproj_src )
        use simple_eul_prob_tab,            only: eul_prob_tab
        use simple_matcher_smpl_and_lplims, only: sample_ptcls4fillin, sample_ptcls4update3D
        use simple_builder,                 only: builder
        class(cmdline),              intent(inout) :: cline
        class(sp_project), optional, intent(in)    :: spproj_src
        integer,     allocatable :: pinds(:)
        type(string)             :: fname
        type(builder)            :: build
        type(parameters)         :: params
        type(commander_prob_tab) :: xprob_tab
        type(eul_prob_tab)       :: eulprob_obj_glob
        type(cmdline)            :: cline_prob_tab
        type(qsys_env)           :: qenv
        type(chash)              :: job_descr
        integer :: nptcls
        integer(timer_int_kind) :: t_tot, t
        real(timer_int_kind)    :: rt_setup, rt_sample, rt_sched, rt_read, rt_assign, rt_write
        if( L_BENCH_GLOB ) t_tot = tic()
        call cline%set('mkdir',  'no')
        call cline%set('stream', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false., spproj_src=spproj_src)
        if( L_BENCH_GLOB )then
            rt_setup = toc(t_tot)
            t        = tic()
        endif
        ! Reset sampling history only on the very first iteration of the run.
        ! Guarding on startit alone re-cleans every iteration when startit==1,
        ! which prevents updatecnt/sampled from accumulating and makes the
        ! fractional-update subset non-progressing (see exec_prob_align2D).
        if( params%startit == 1 .and. params%which_iter == params%startit )then
            call build%spproj_field%clean_entry('updatecnt', 'sampled')
        endif
        ! sampled incremented
        if( params%l_fillin .and. mod(params%which_iter,5) == 0 )then
            call sample_ptcls4fillin(params, build, [1,params%nptcls], .true., nptcls, pinds)
        else
            call sample_ptcls4update3D(params, build, [1,params%nptcls], .true., nptcls, pinds)
        endif
        ! communicate to project file
        call build%spproj%write_segment_inside(params%oritype)
        call cleanup_prob_align_outputs(params, .false.)
        if( L_BENCH_GLOB )then
            rt_sample = toc(t)
            t         = tic()
        endif
        ! generating all corrs on all parts
        cline_prob_tab = cline
        call cline_prob_tab%set('prg', 'prob_tab' ) ! required for distributed call
        ! execution
        if( .not.cline_prob_tab%defined('nparts') )then
            call xprob_tab%execute(cline_prob_tab)
        else
            ! setup the environment for distributed execution: partitions balance the particles with state > 0
            call qenv%new(params, params%nparts, nptcls=params%nptcls, l_active=build%spproj_field%included())
            call cline_prob_tab%gen_job_descr(job_descr)
            ! schedule
            call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        endif
        if( L_BENCH_GLOB )then
            rt_sched = toc(t)
            t        = tic()
        endif
        ! Build the global table only after worker tables are complete.  Keeping it
        ! live while workers build dense partition tables roughly doubles peak RSS.
        call eulprob_obj_glob%new(params,build,pinds)
        call eulprob_obj_glob%read_tabs_to_glob(string(DIST_FBODY), params%nparts, params%numlen)
        if( L_BENCH_GLOB )then
            rt_read = toc(t)
            t       = tic()
        endif
        call eulprob_obj_glob%ref_assign
        if( L_BENCH_GLOB )then
            rt_assign = toc(t)
            t         = tic()
        endif
        ! write the iptcl->(iref,istate) assignment
        fname = string(ASSIGNMENT_FBODY)//'.dat'
        call eulprob_obj_glob%write_assignment(fname)
        if( L_BENCH_GLOB ) rt_write = toc(t)
        call eulprob_obj_glob%kill
        ! cleanup
        call cline_prob_tab%kill
        call qenv%kill
        call job_descr%kill
        if( L_BENCH_GLOB ) call write_prob_align_bench(params, 'PROB_ALIGN', nptcls, rt_setup, rt_sample, rt_sched, &
            &rt_read, rt_assign, rt_write, toc(t_tot))
        call build%kill_general_tbox
        call qsys_declare_part_finished(params, string('simple_commanders_refine3D :: exec_prob_align'))
        call qsys_cleanup(params)
        call simple_end('**** SIMPLE_PROB_ALIGN NORMAL STOP ****', print_simple=.false.)
    end subroutine run_prob_align

    subroutine exec_prob_align_neigh( self, cline )
        class(commander_prob_align_neigh), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        call run_prob_align_neigh(cline)
    end subroutine exec_prob_align_neigh

    !> prob_align_neigh with the caller's project of projfile in memory, read and unchanged since: copied
    subroutine exec_prob_align_neigh_from_project( self, cline, spproj )
        class(commander_prob_align_neigh), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        class(sp_project),                 intent(in)    :: spproj
        call run_prob_align_neigh(cline, spproj)
    end subroutine exec_prob_align_neigh_from_project

    subroutine run_prob_align_neigh( cline, spproj_src )
        use simple_eul_prob_tab_neigh,      only: eul_prob_tab_neigh
        use simple_matcher_smpl_and_lplims, only: sample_ptcls4fillin, sample_ptcls4update3D
        use simple_builder,                 only: builder
        class(cmdline),              intent(inout) :: cline
        class(sp_project), optional, intent(in)    :: spproj_src
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
        integer(timer_int_kind) :: t_tot, t
        real(timer_int_kind)    :: rt_setup, rt_sample, rt_sched, rt_read, rt_assign, rt_write
        if( L_BENCH_GLOB ) t_tot = tic()
        call cline%set('mkdir',  'no')
        call cline%set('stream', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true., spproj_src=spproj_src)
        if( L_BENCH_GLOB )then
            rt_setup = toc(t_tot)
            t        = tic()
        endif
        ! Reset sampling history only on the very first iteration of the run.
        ! Guarding on startit alone re-cleans every iteration when startit==1,
        ! which prevents updatecnt/sampled from accumulating and makes the
        ! fractional-update subset non-progressing (see exec_prob_align2D).
        if( params%startit == 1 .and. params%which_iter == params%startit )then
            call build%spproj_field%clean_entry('updatecnt', 'sampled')
        endif
        if( params%l_fillin .and. mod(params%which_iter,5) == 0 )then
            call sample_ptcls4fillin(params, build, [1,params%nptcls], .true., nptcls, pinds)
        else
            call sample_ptcls4update3D(params, build, [1,params%nptcls], .true., nptcls, pinds)
        endif
        call build%spproj%write_segment_inside(params%oritype)
        call cleanup_prob_align_outputs(params, .true.)
        if( L_BENCH_GLOB )then
            rt_sample = toc(t)
            t         = tic()
        endif
        cline_prob_tab = cline
        call cline_prob_tab%set('prg', 'prob_tab_neigh')
        if( .not. cline_prob_tab%defined('nparts') )then
            call xprob_tab_neigh%execute(cline_prob_tab)
        else
            call qenv%new(params, params%nparts, nptcls=params%nptcls, l_active=build%spproj_field%included())
            call cline_prob_tab%gen_job_descr(job_descr)
            call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        endif
        if( L_BENCH_GLOB )then
            rt_sched = toc(t)
            t        = tic()
        endif
        ! Construct global storage only after worker scoring has released its table.
        call eulprob_obj_glob_neigh%new_neigh_global(params,build,pinds)
        call eulprob_obj_glob_neigh%read_tabs_to_glob(string(DIST_FBODY)//'_neigh_', params%nparts, params%numlen)
        if( L_BENCH_GLOB )then
            rt_read = toc(t)
            t       = tic()
        endif
        call eulprob_obj_glob_neigh%ref_assign
        if( L_BENCH_GLOB )then
            rt_assign = toc(t)
            t         = tic()
        endif
        ! write the iptcl->(iref,istate) assignment
        fname = string(ASSIGNMENT_FBODY)//'.dat'
        call eulprob_obj_glob_neigh%write_assignment(fname)
        if( L_BENCH_GLOB ) rt_write = toc(t)
        call eulprob_obj_glob_neigh%kill
        ! cleanup
        call cline_prob_tab%kill
        call qenv%kill
        call job_descr%kill
        if( L_BENCH_GLOB ) call write_prob_align_bench(params, 'PROB_ALIGN_NEIGH', nptcls, rt_setup, rt_sample, &
            &rt_sched, rt_read, rt_assign, rt_write, toc(t_tot))
        call build%kill_general_tbox
        call qsys_declare_part_finished(params, string('simple_commanders_refine3D :: exec_prob_align_neigh'))
        call qsys_cleanup(params)
        call simple_end('**** SIMPLE_PROB_ALIGN_NEIGH NORMAL STOP ****', print_simple=.false.)
    end subroutine run_prob_align_neigh

    subroutine exec_prob_tab2D( self, cline )
        use simple_matcher_smpl_and_lplims, only: set_bp_range2D
        use simple_strategy2D_matcher,      only: set_b_p_ptrs2D, &
            &ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad
        use simple_matcher_pftc_prep,       only: prep_pftc4align2D
        use simple_matcher_ptcl_batch,      only: alloc_ptcl_imgs, build_batch_particles2D, clean_batch_particles2D
        use simple_ptcl_cache,              only: ptcl_cache_in_use, ptcl_cache_assert_ready
        use simple_classaverager,           only: cavger_new, cavger_read_all, cavger_kill
        use simple_eul_prob_tab2D,          only: eul_prob_tab2D
        use simple_eul_prob_tab_utils,      only: prob_fill_profile
        use simple_matcher_ptcl_io,         only: ptcl_read_counts, ptcl_read_counts_now, ptcl_read_counts_since, &
            &ptcl_read_session, ptcl_buffer_budget, fill_pass_buffer
        class(commander_prob_tab2D), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        integer,     allocatable :: pinds(:)
        type(string)             :: fname
        type(builder)            :: build
        type(parameters)         :: params
        type(eul_prob_tab2D)     :: eulprob_obj_part
        type(prob_fill_profile)  :: fill_prof
        type(ptcl_read_counts)   :: reads0
        type(ptcl_read_session)  :: reads
        real    :: frac_srch_space
        integer :: nptcls, batchsz_max, nbatches, ibatch, batch_start, batch_end, batchsz, buffer_status
        integer, allocatable :: batches(:,:)
        logical :: l_cached
        integer(timer_int_kind) :: t_tot, t
        real(timer_int_kind)    :: rt_setup, rt_refs, rt_batch, rt_fill, rt_write
        if( L_BENCH_GLOB )then
            t_tot  = tic()
            reads0 = ptcl_read_counts_now()
        endif
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        call ensure_field_halves(build%spproj, params%oritype, params%projfile)
        frac_srch_space  = build%spproj_field%get_avg('frac')
        call set_bp_range2D(params, build, cline, params%which_iter, frac_srch_space)
        ! reproduce particle sampling from exec_prob_align2D
        if( build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls, pinds, &
                &allow_empty=params%l_distr_worker)
        else
            THROW_HARD('exec_prob_tab2D requires particle sampling from exec_prob_align2D')
        endif
        fname = string(DIST_FBODY)//int2str_pad(params%part,params%numlen)//'.dat'
        if( nptcls < 1 )then
            ! a partition without sampled particles contributes an empty table
            call eulprob_obj_part%new_worker(params, build, pinds)
            call eulprob_obj_part%begin_write(fname)
            call eulprob_obj_part%write_tab(fname)
            call eulprob_obj_part%kill
            call build%kill_general_tbox
            call qsys_declare_part_finished(params, string('simple_commanders_prob :: exec_prob_tab2D'))
            call simple_end('**** SIMPLE_PROB_TAB2D NORMAL STOP ****', print_simple=.false.)
            return
        endif
        batchsz_max = min(nptcls, params%nthr * BATCHTHRSZ)
        nbatches    = ceiling(real(nptcls) / real(batchsz_max))
        batches     = split_nobjs_even(nptcls, nbatches)
        batchsz_max = maxval(batches(:,2) - batches(:,1) + 1)
        if( L_BENCH_GLOB )then
            rt_setup = toc(t_tot)
            t        = tic()
        endif
        call set_b_p_ptrs2D(params, build)
        ! no ptcl_imgs here: prob_tab2D scores against the references and never
        ! restores class averages, so the raw full-size images are not needed and the
        ! batch can be served from the downscaled cache
        call ptcl_cache_assert_ready(params, build)
        l_cached = ptcl_cache_in_use(params, build)
        if( l_cached )then
            call alloc_ptcl_imgs(params, build, ptcl_match_imgs, ptcl_match_imgs_pad, batchsz_max, &
                &imgbatch_box=params%box_crop, imgbatch_smpd=params%smpd_crop)
        else
            call alloc_ptcl_imgs(params, build, ptcl_match_imgs, ptcl_match_imgs_pad, batchsz_max)
        endif
        ! mirror refine2D_exec reference setup
        call cavger_new(params, build)
        if( .not. cline%defined('refs') ) THROW_HARD('exec_prob_tab2D requires refs on the command line')
        call cavger_read_all
        call prep_pftc4align2D(params, build, ptcl_match_imgs_pad, batchsz_max, params%which_iter)
        ! Fill the partition table in matcher-sized batches to cap polar FT memo memory.
        call eulprob_obj_part%new_worker(params,build,pinds)
        call eulprob_obj_part%begin_write(fname)
        if( L_BENCH_GLOB )then
            rt_refs  = toc(t)
            rt_batch = 0.
            rt_fill  = 0.
        endif
        call reads%new(params)
        call reads%expect(pinds(1:nptcls))
        ! in shared memory the matcher of the same iteration and process searches the same particles: the table
        ! reads them once into the particle buffer, and the matcher's batches are copied from it (it releases it)
        if( .not. params%l_distr_worker .and. .not. l_cached )then
            buffer_status = fill_pass_buffer(params, build, reads, nptcls, pinds, ptcl_buffer_budget(params, build))
        endif
        do ibatch = 1, nbatches
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            batchsz     = batch_end - batch_start + 1
            if( L_BENCH_GLOB ) t = tic()
            call build_batch_particles2D(params, build, batchsz, pinds(batch_start:batch_end), reads, &
                &ptcl_match_imgs=ptcl_match_imgs, ptcl_match_imgs_pad=ptcl_match_imgs_pad)
            if( L_BENCH_GLOB )then
                rt_batch = rt_batch + toc(t)
                t        = tic()
            endif
            call eulprob_obj_part%fill_tab_range(batch_start, batch_end)
            if( L_BENCH_GLOB ) rt_fill = rt_fill + toc(t)
        end do
        call reads%kill
        ! write the 2D probability table
        if( L_BENCH_GLOB ) t = tic()
        call eulprob_obj_part%write_tab(fname)
        if( L_BENCH_GLOB ) rt_write = toc(t)
        fill_prof = eulprob_obj_part%get_fill_profile()
        call eulprob_obj_part%kill
        if( allocated(batches) ) deallocate(batches)
        call clean_batch_particles2D(build, ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad)
        if( L_BENCH_GLOB ) call write_prob_tab_bench(params, nptcls, ptcl_read_counts_since(reads0), &
            &rt_setup, rt_refs, rt_batch, rt_fill, rt_write, toc(t_tot), fill=fill_prof)
        call cavger_kill
        call build%pftc%kill
        call build%kill_general_tbox
        call qsys_declare_part_finished(params, string('simple_commanders_prob :: exec_prob_tab2D'))
        call simple_end('**** SIMPLE_PROB_TAB2D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_tab2D

    subroutine exec_prob_align2D( self, cline )
        use simple_eul_prob_tab2D,          only: eul_prob_tab2D
        use simple_strategy2D_matcher,      only: set_b_p_ptrs2D
        use simple_matcher_smpl_and_lplims, only: sample_ptcls4update2D
        use simple_ptcl_cache,              only: ptcl_cache_ensure
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
        integer :: nptcls
        integer(timer_int_kind) :: t_tot, t
        real(timer_int_kind)    :: rt_setup, rt_sample, rt_sched, rt_read, rt_assign, rt_write
        if( L_BENCH_GLOB ) t_tot = tic()
        call cline%set('mkdir',  'no')
        call cline%set('stream', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        call set_b_p_ptrs2D(params, build)
        call ensure_field_halves(build%spproj, params%oritype, params%projfile)
        if( L_BENCH_GLOB )then
            rt_setup = toc(t_tot)
            t        = tic()
        endif
        if( params%startit == 1 .and. params%which_iter == params%startit )then
            call build%spproj_field%clean_entry('updatecnt', 'sampled')
        endif
        ! Mirror the 3D workflow: sampled-update is active from the first stage onward.
        ! In probabilistic mode the sampled subset is reused within the current iteration
        ! by prob_tab2D/refine2D_exec, but it is redrawn on later iterations.
        call sample_ptcls4update2D(params, build, [params%fromp,params%top], params%l_update_frac, nptcls, pinds)
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> PROB_ALIGN2D: sampled ', nptcls, ' particles over ', params%nparts, ' part(s)'
        call flush(logfhandle)
        ! write sampling to project
        call build%spproj%write_segment_inside(params%oritype)
        ! Build the downscaled particle cache before the workers start, so they all
        ! find it ready. A no-op once it exists, so the per-iteration cost is a file
        ! and key check; this is also the first point in the workflow where the same
        ! particles get read twice per iteration, once here and once in refine2D.
        ! Must precede the cline copy below so a fallback to cache=no reaches the
        ! prob_tab2D worker command lines.
        call ptcl_cache_ensure(params, build, cline)
        if( L_BENCH_GLOB )then
            rt_sample = toc(t)
            t         = tic()
        endif
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
        write(logfhandle,'(A)') '>>> PROB_ALIGN2D: prob_tab2D workers completed; merging partition tables'
        call flush(logfhandle)
        if( L_BENCH_GLOB )then
            rt_sched = toc(t)
            t        = tic()
        endif
        ! Allocate the global table only after worker scoring has released its batch-local storage.
        call eulprob_obj_glob%new(params, build, pinds)
        call eulprob_obj_glob%read_tabs_to_glob(string(DIST_FBODY),params%nparts,params%numlen)
        if( L_BENCH_GLOB )then
            rt_read = toc(t)
            t       = tic()
        endif
        ! global probabilistic class assignment
        write(logfhandle,'(A)') '>>> PROB_ALIGN2D: running global probabilistic assignment'
        call flush(logfhandle)
        call eulprob_obj_glob%ref_assign
        if( L_BENCH_GLOB )then
            rt_assign = toc(t)
            t         = tic()
        endif
        ! write assignment to file
        fname = string(ASSIGNMENT_FBODY)//'.dat'
        write(logfhandle,'(A,A)') '>>> PROB_ALIGN2D: writing assignment ', fname%to_char()
        call flush(logfhandle)
        call eulprob_obj_glob%write_assignment(fname)
        write(logfhandle,'(A)') '>>> PROB_ALIGN2D: assignment written'
        call flush(logfhandle)
        if( L_BENCH_GLOB ) rt_write = toc(t)
        ! cleanup
        call eulprob_obj_glob%kill
        call cline_prob_tab%kill
        call qenv%kill
        call job_descr%kill
        if( L_BENCH_GLOB ) call write_prob_align_bench(params, 'PROB_ALIGN2D', nptcls, rt_setup, rt_sample, rt_sched, &
            &rt_read, rt_assign, rt_write, toc(t_tot))
        call build%kill_general_tbox
        call qsys_declare_part_finished(params, string('simple_commanders_prob :: exec_prob_align2D'))
        call qsys_cleanup(params)
        call simple_end('**** SIMPLE_PROB_ALIGN2D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_align2D

    subroutine cleanup_prob_align_outputs( params, neigh )
        class(parameters), intent(in) :: params
        logical,           intent(in) :: neigh
        type(string) :: fname
        integer :: ipart
        if( neigh )then
            do ipart = 1, params%nparts
                fname = string(DIST_FBODY)//'_neigh_'//int2str_pad(ipart,params%numlen)//'.dat'
                call del_file(fname)
            end do
        else
            call del_files(DIST_FBODY, params%nparts, ext='.dat', numlen=params%numlen)
        endif
        call del_file(string(ASSIGNMENT_FBODY)//'.dat')
        call fname%kill
    end subroutine cleanup_prob_align_outputs

    !> The bench record of a probability-table part, in its partition's file of the prob_tab, prob_tab_neigh
    !! or prob_tab2D family; rt_batch is the batch builds, particle reading included, which is reported apart
    subroutine write_prob_tab_bench( params, nptcls, reads, rt_setup, rt_refs, rt_batch, rt_fill, rt_write, rt_tot, fill )
        use simple_matcher_ptcl_io,    only: ptcl_read_counts, write_ptcl_read_bench_context, write_mem_bench_context
        use simple_bench_io,           only: open_bench_report
        use simple_eul_prob_tab_utils, only: prob_fill_profile
        class(parameters),      intent(in) :: params
        integer,                intent(in) :: nptcls
        type(ptcl_read_counts), intent(in) :: reads
        real(timer_int_kind),   intent(in) :: rt_setup, rt_refs, rt_batch, rt_fill, rt_write, rt_tot
        type(prob_fill_profile), optional, intent(in) :: fill   !< cost record of the table fill
        real(timer_int_kind) :: rt_read
        real(dp) :: ts_sum
        integer :: fnr
        rt_read = reads%rt_read + reads%rt_cache
        call open_bench_report(fnr, uppercase(params%prg%to_char()), iter=params%which_iter, part=params%part, &
            &numlen=params%numlen)
        write(fnr,'(a)') '*** BENCHMARK CONTEXT ***'
        write(fnr,'(a,a)')  'prob_tab program                    : ', params%prg%to_char()
        write(fnr,'(a,i0)') 'prob_tab process partition          : ', params%part
        write(fnr,'(a,i0)') 'prob_tab process pid                : ', params%pid
        write(fnr,'(a,i0)') 'prob_tab nparts                     : ', params%nparts
        write(fnr,'(a,i0)') 'prob_tab worker threads             : ', params%nthr
        write(fnr,'(a,i0)') 'prob_tab box                        : ', params%box
        write(fnr,'(a,i0)') 'prob_tab box_crop                   : ', params%box_crop
        write(fnr,'(a,i0)') 'prob_tab sampled particles          : ', nptcls
        if( present(fill) )then
            if( fill%nscored > 0 )then
                ! a neighbourhood table (and the prob_snhc 2D table) scores a different set of references per particle
                write(fnr,'(a,i0)') 'prob_tab scored references (all)    : ', fill%nscored
                if( fill%nrefined > 0 ) write(fnr,'(a,i0)') 'prob_tab refined refs (all)         : ', fill%nrefined
            else
                write(fnr,'(a,i0)') 'prob_tab references                 : ', fill%nrefs
                write(fnr,'(a,i0)') 'prob_tab refined refs per particle  : ', fill%projs_ns
                write(fnr,'(a,i0)') 'prob_tab refined refs (all)         : ', fill%nrefined
            endif
            write(fnr,'(a,i0)') 'prob_tab shift searches             : ', fill%minimize_calls
            write(fnr,'(a,i0)') 'prob_tab alternations               : ', fill%alternations
            write(fnr,'(a,i0)') 'prob_tab grid evaluations           : ', fill%grid_points
            write(fnr,'(a,i0)') 'prob_tab objective evaluations      : ', fill%objective_evals
            write(fnr,'(a,i0)') 'prob_tab gradient evaluations (E)   : ', fill%gradient_evals
            if( params%nstates > 1 )then
                write(fnr,'(a,i0)') 'prob_tab grid fits kept             : ', fill%fit_counts(0)
                write(fnr,'(a,i0)') 'prob_tab fit fallback, edge         : ', fill%fit_counts(1)
                write(fnr,'(a,i0)') 'prob_tab fit fallback, curvature    : ', fill%fit_counts(2)
                write(fnr,'(a,i0)') 'prob_tab fit fallback, step         : ', fill%fit_counts(3)
                write(fnr,'(a,i0)') 'prob_tab fit fallback, search box   : ', fill%fit_counts(4)
                write(fnr,'(a,i0)') 'prob_tab fit fallback, value        : ', fill%fit_counts(5)
            endif
            write(fnr,'(a,i0)') 'prob_tab shifted particle memos     : ', fill%bank_builds
            if( fill%minimize_calls > 0 )then
                write(fnr,'(a,1x,f0.2)') 'prob_tab alternations per search      :', &
                    &real(fill%alternations,dp) / real(fill%minimize_calls,dp)
                write(fnr,'(a,1x,f0.2)') 'prob_tab gradient evals per search    :', &
                    &real(fill%gradient_evals,dp) / real(fill%minimize_calls,dp)
            endif
        endif
        call write_mem_bench_context(fnr, 'prob_tab')
        call write_ptcl_read_bench_context(fnr, 'prob_tab', reads)
        write(fnr,'(a)') ''
        write(fnr,'(a)') '*** TIMINGS (s) ***'
        write(fnr,'(a,1x,f0.2)') 'prob_tab startup/setup              :', rt_setup
        write(fnr,'(a,1x,f0.2)') 'prob_tab reference preparation      :', rt_refs
        write(fnr,'(a,1x,f0.2)') 'prob_tab particle reading           :', rt_read
        write(fnr,'(a,1x,f0.2)') 'prob_tab particle preparation       :', rt_batch - rt_read
        write(fnr,'(a,1x,f0.2)') 'prob_tab table fill                 :', rt_fill
        write(fnr,'(a,1x,f0.2)') 'prob_tab table write                :', rt_write
        write(fnr,'(a,1x,f0.2)') 'prob_tab total time                 :', rt_tot
        if( rt_tot > 0. )then
            write(fnr,'(a,1x,f0.2)') 'prob_tab % accounted for            :', &
                &((rt_setup + rt_refs + rt_batch + rt_fill + rt_write) / rt_tot) * 100.
        endif
        if( present(fill) )then
            ! thread-seconds are never part of '% accounted for'
            ts_sum = fill%ts_seed + fill%ts_score + fill%ts_refine
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** THREAD-SECONDS (summed over threads; not wall time) ***'
            write(fnr,'(a,1x,f0.2)') 'prob_tab seed estimation            :', fill%ts_seed
            write(fnr,'(a,1x,f0.2)') 'prob_tab score loop                 :', fill%ts_score
            write(fnr,'(a,1x,f0.2)') 'prob_tab refinement                 :', fill%ts_refine
            write(fnr,'(a,1x,f0.2)') 'prob_tab thread-seconds, sum        :', ts_sum
            write(fnr,'(a,1x,f0.2)') 'prob_tab fill region wall time      :', fill%wall_region
            if( fill%wall_region > 0.d0 .and. fill%nthr > 0 )then
                write(fnr,'(a,1x,f0.3)') 'prob_tab load balance               :', &
                    &ts_sum / (fill%wall_region * real(fill%nthr,dp))
            endif
            if( fill%nscored > 0 )then
                write(fnr,'(a,1x,es10.3)') 'prob_tab score loop s/(particle*ref) :', &
                    &fill%ts_score / real(fill%nscored,dp)
            else if( fill%nptcls > 0 .and. fill%nrefs > 0 )then
                write(fnr,'(a,1x,es10.3)') 'prob_tab score loop s/(particle*ref) :', &
                    &fill%ts_score / (real(fill%nptcls,dp) * real(fill%nrefs,dp))
            endif
        endif
        call fclose(fnr)
    end subroutine write_prob_tab_bench

    !> The bench record of the prob_align driver, all wall time, in family PROB_ALIGN (prob), PROB_ALIGN_NEIGH
    !! (prob_neigh) or PROB_ALIGN2D; rt_sched runs from the launch of the parts to the last part finished (the
    !! in-process table fill in shared memory).
    subroutine write_prob_align_bench( params, family, nptcls, rt_setup, rt_sample, rt_sched, rt_read, rt_assign, &
            &rt_write, rt_tot )
        use simple_bench_io, only: open_bench_report
        class(parameters),    intent(in) :: params
        character(len=*),     intent(in) :: family
        integer,              intent(in) :: nptcls
        real(timer_int_kind), intent(in) :: rt_setup, rt_sample, rt_sched, rt_read, rt_assign, rt_write, rt_tot
        integer :: fnr
        call open_bench_report(fnr, family, iter=params%which_iter)
        write(fnr,'(a)') '*** BENCHMARK CONTEXT ***'
        write(fnr,'(a,i0)') 'prob_align iteration                : ', params%which_iter
        write(fnr,'(a,i0)') 'prob_align nparts                   : ', params%nparts
        write(fnr,'(a,i0)') 'prob_align worker threads           : ', params%nthr
        if( family == 'PROB_ALIGN2D' )then
            write(fnr,'(a,i0)') 'prob_align classes                  : ', params%ncls
        else
            write(fnr,'(a,i0)') 'prob_align states                   : ', params%nstates
            write(fnr,'(a,i0)') 'prob_align projection directions    : ', params%nspace
        endif
        write(fnr,'(a,i0)') 'prob_align box_crop                 : ', params%box_crop
        write(fnr,'(a,i0)') 'prob_align sampled particles        : ', nptcls
        write(fnr,'(a)') ''
        write(fnr,'(a)') '*** TIMINGS (s) ***'
        write(fnr,'(a,1x,f0.2)') 'prob_align startup/setup            :', rt_setup
        write(fnr,'(a,1x,f0.2)') 'prob_align sampling, project write  :', rt_sample
        write(fnr,'(a,1x,f0.2)') 'prob_align scheduling (parts)       :', rt_sched
        write(fnr,'(a,1x,f0.2)') 'prob_align table reads              :', rt_read
        write(fnr,'(a,1x,f0.2)') 'prob_align ref_assign               :', rt_assign
        write(fnr,'(a,1x,f0.2)') 'prob_align assignment write         :', rt_write
        write(fnr,'(a,1x,f0.2)') 'prob_align total time               :', rt_tot
        if( rt_tot > 0. )then
            write(fnr,'(a,1x,f0.2)') 'prob_align % accounted for          :', &
                &((rt_setup + rt_sample + rt_sched + rt_read + rt_assign + rt_write) / rt_tot) * 100.
        endif
        call fclose(fnr)
    end subroutine write_prob_align_bench

    !> The table of a partition without sampled particles: the header, no
    !! candidates, so prob_align reads one table per partition. Written by the
    !! table type of the caller (dense or neighbourhood), which owns the format.
    subroutine write_empty_prob_tab( params, build, pinds, fname, tab )
        use simple_eul_prob_tab, only: eul_prob_tab
        class(parameters),   intent(in)    :: params
        class(builder),      intent(in)    :: build
        integer,             intent(in)    :: pinds(:)
        class(string),       intent(in)    :: fname
        class(eul_prob_tab), intent(inout) :: tab
        call tab%new_worker(params, build, pinds)
        call tab%begin_write(fname)
        call tab%write_tab(fname)
        call tab%kill
    end subroutine write_empty_prob_tab

end module simple_commanders_prob
