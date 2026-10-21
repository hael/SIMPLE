!@descr: probability alignment of refine3D and refine2D: the tables' sampling, assignment and its output
module simple_prob_align_strategy
use simple_core_module_api
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_parameters,     only: parameters
use simple_commander_base, only: commander_base
use simple_qsys_env,       only: qsys_env
use simple_qsys_funs,      only: qsys_declare_part_finished, qsys_cleanup
!$ use omp_lib,            only: omp_set_num_threads
implicit none
private
public :: run_prob_align, run_prob_align_neigh, prob_align_exec, prob_align_neigh_exec, run_prob_align2D, prob_align2D_exec
public :: write_prob_align_bench
#include "simple_local_flags.inc"

contains

    !> prob_align on its own parameters and builder, built from cline: samples the particles to update, fills the dense
    !! probability table (in process through the caller's table commander xprob_tab, or on nparts partitions), assigns
    !! and writes the assignment
    subroutine run_prob_align( cline, xprob_tab )
        class(cmdline),                  intent(inout) :: cline
        class(commander_base), optional, intent(inout) :: xprob_tab
        type(builder)           :: build
        type(parameters)        :: params
        integer(timer_int_kind) :: t_setup
        real(timer_int_kind)    :: rt_setup
        rt_setup = 0.
        if( L_BENCH_GLOB ) t_setup = tic()
        call cline%set('mkdir',  'no')
        call cline%set('stream', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        if( L_BENCH_GLOB ) rt_setup = toc(t_setup)
        call prob_align_exec(params, build, cline, params%nthr, xprob_tab=xprob_tab, rt_setup=rt_setup)
        call build%kill_general_tbox
        call qsys_declare_part_finished(params, string('simple_commanders_refine3D :: exec_prob_align'))
        call simple_end('**** SIMPLE_PROB_ALIGN NORMAL STOP ****', print_simple=.false.)
    end subroutine run_prob_align

    !> the dense prob_align on the caller's parameters and builder: nthr threads (the workers', as the command line
    !! gives them; the caller restores its own), the sample drawn in the caller's particle field and written to the
    !! project, the table filled (in process through xprob_tab, or on nparts partitions), assigned and written;
    !! rt_setup: the caller's set-up time, for the bench record
    subroutine prob_align_exec( params, build, cline, nthr, xprob_tab, rt_setup )
        use simple_eul_prob_tab,            only: eul_prob_tab
        use simple_matcher_smpl_and_lplims, only: sample_ptcls4fillin, sample_ptcls4update3D
        class(parameters), target,       intent(inout) :: params
        class(builder),    target,       intent(inout) :: build
        class(cmdline),                  intent(inout) :: cline
        integer,                         intent(in)    :: nthr
        class(commander_base), optional, intent(inout) :: xprob_tab
        real(timer_int_kind),  optional, intent(in)    :: rt_setup
        integer,     allocatable :: pinds(:)
        type(string)             :: fname
        type(eul_prob_tab)       :: eulprob_obj_glob
        type(cmdline)            :: cline_prob_tab
        type(qsys_env)           :: qenv
        type(chash)              :: job_descr
        integer :: nptcls
        integer(timer_int_kind) :: t_tot, t
        real(timer_int_kind)    :: rt_sample, rt_sched, rt_read, rt_assign, rt_write
        real(timer_int_kind)    :: rt_setup_here
        rt_setup_here = 0.
        if( present(rt_setup) ) rt_setup_here = rt_setup
        if( L_BENCH_GLOB )then
            t_tot = tic()
            t     = t_tot
        endif
        !$ call omp_set_num_threads(nthr)
        nthr_glob = nthr
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
            if( .not. present(xprob_tab) ) THROW_HARD('prob_align in shared memory needs the table commander')
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
        if( L_BENCH_GLOB ) call write_prob_align_bench(params, 'PROB_ALIGN', nptcls, rt_setup_here, rt_sample, rt_sched, &
            &rt_read, rt_assign, rt_write, toc(t_tot))
        call qsys_cleanup(params)
    end subroutine prob_align_exec

    !> prob_align_neigh: as run_prob_align with the neighbourhood table (xprob_tab_neigh fills it in process)
    subroutine run_prob_align_neigh( cline, xprob_tab_neigh )
        class(cmdline),                  intent(inout) :: cline
        class(commander_base), optional, intent(inout) :: xprob_tab_neigh
        type(builder)           :: build
        type(parameters)        :: params
        integer(timer_int_kind) :: t_setup
        real(timer_int_kind)    :: rt_setup
        rt_setup = 0.
        if( L_BENCH_GLOB ) t_setup = tic()
        call cline%set('mkdir',  'no')
        call cline%set('stream', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        if( L_BENCH_GLOB ) rt_setup = toc(t_setup)
        call prob_align_neigh_exec(params, build, cline, params%nthr, xprob_tab_neigh=xprob_tab_neigh, rt_setup=rt_setup)
        call build%kill_general_tbox
        call qsys_declare_part_finished(params, string('simple_commanders_refine3D :: exec_prob_align_neigh'))
        call simple_end('**** SIMPLE_PROB_ALIGN_NEIGH NORMAL STOP ****', print_simple=.false.)
    end subroutine run_prob_align_neigh

    !> the neighbourhood prob_align on the caller's parameters and builder: nthr threads (the workers', as the command line
    !! gives them; the caller restores its own), the sample drawn in the caller's particle field and written to the
    !! project, the table filled (in process through xprob_tab_neigh, or on nparts partitions), assigned and written;
    !! rt_setup: the caller's set-up time, for the bench record
    subroutine prob_align_neigh_exec( params, build, cline, nthr, xprob_tab_neigh, rt_setup )
        use simple_eul_prob_tab_neigh,      only: eul_prob_tab_neigh
        use simple_matcher_smpl_and_lplims, only: sample_ptcls4fillin, sample_ptcls4update3D
        class(parameters), target,       intent(inout) :: params
        class(builder),    target,       intent(inout) :: build
        class(cmdline),                  intent(inout) :: cline
        integer,                         intent(in)    :: nthr
        class(commander_base), optional, intent(inout) :: xprob_tab_neigh
        real(timer_int_kind),  optional, intent(in)    :: rt_setup
        integer,           allocatable :: pinds(:)
        type(string)                   :: fname
        type(eul_prob_tab_neigh)       :: eulprob_obj_glob_neigh
        type(cmdline)                  :: cline_prob_tab
        type(qsys_env)                 :: qenv
        type(chash)                    :: job_descr
        integer :: nptcls
        integer(timer_int_kind) :: t_tot, t
        real(timer_int_kind)    :: rt_sample, rt_sched, rt_read, rt_assign, rt_write
        real(timer_int_kind)    :: rt_setup_here
        rt_setup_here = 0.
        if( present(rt_setup) ) rt_setup_here = rt_setup
        if( L_BENCH_GLOB )then
            t_tot = tic()
            t     = t_tot
        endif
        !$ call omp_set_num_threads(nthr)
        nthr_glob = nthr
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
            if( .not. present(xprob_tab_neigh) ) THROW_HARD('prob_align_neigh in shared memory needs the table commander')
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
        if( L_BENCH_GLOB ) call write_prob_align_bench(params, 'PROB_ALIGN_NEIGH', nptcls, rt_setup_here, rt_sample, &
            &rt_sched, rt_read, rt_assign, rt_write, toc(t_tot))
        call qsys_cleanup(params)
    end subroutine prob_align_neigh_exec


    !> prob_align2D on its own parameters and builder, built from cline: samples the particles to update, fills the 2D
    !! probability table (in process through the caller's table commander xprob_tab2D, or on nparts partitions), assigns
    !! the classes and writes the assignment
    subroutine run_prob_align2D( cline, xprob_tab2D )
        class(cmdline),                  intent(inout) :: cline
        class(commander_base), optional, intent(inout) :: xprob_tab2D
        type(builder)           :: build
        type(parameters)        :: params
        integer(timer_int_kind) :: t_setup
        real(timer_int_kind)    :: rt_setup
        rt_setup = 0.
        if( L_BENCH_GLOB ) t_setup = tic()
        call cline%set('mkdir',  'no')
        call cline%set('stream', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        if( L_BENCH_GLOB ) rt_setup = toc(t_setup)
        call prob_align2D_exec(params, build, cline, params%nthr, xprob_tab2D, rt_setup)
        call build%kill_general_tbox
        call qsys_declare_part_finished(params, string('simple_commanders_prob :: exec_prob_align2D'))
        call simple_end('**** SIMPLE_PROB_ALIGN2D NORMAL STOP ****', print_simple=.false.)
    end subroutine run_prob_align2D

    !> prob_align2D on the caller's parameters and builder: nthr threads (the workers', as the command line gives them;
    !! the caller restores its own), the halves ensured, the sample drawn in the caller's particle field and written to
    !! the project, the particle cache made ready, the table filled (in process through xprob_tab2D, or on nparts
    !! partitions), assigned and written; rt_setup: the caller's set-up time, for the bench record
    subroutine prob_align2D_exec( params, build, cline, nthr, xprob_tab2D, rt_setup )
        use simple_eul_prob_tab2D,          only: eul_prob_tab2D, prob_assign_profile2D
        use simple_strategy2D_matcher,      only: set_b_p_ptrs2D
        use simple_matcher_smpl_and_lplims, only: sample_ptcls4update2D
        use simple_ptcl_cache,              only: ptcl_cache_ensure
        use simple_eo_partition,            only: ensure_field_halves
        class(parameters), target,       intent(inout) :: params
        class(builder),    target,       intent(inout) :: build
        class(cmdline),                  intent(inout) :: cline
        integer,                         intent(in)    :: nthr
        class(commander_base), optional, intent(inout) :: xprob_tab2D
        real(timer_int_kind),  optional, intent(in)    :: rt_setup
        integer,       allocatable :: pinds(:)
        type(string)               :: fname
        type(eul_prob_tab2D)       :: eulprob_obj_glob
        type(cmdline)              :: cline_prob_tab
        type(qsys_env)             :: qenv
        type(chash)                :: job_descr
        type(prob_assign_profile2D) :: assign_profile
        integer :: nptcls
        integer(timer_int_kind) :: t_tot, t
        real(timer_int_kind)    :: rt_setup_here, rt_sample, rt_sched, rt_read, rt_assign, rt_write
        rt_setup_here = 0.
        if( present(rt_setup) ) rt_setup_here = rt_setup
        if( L_BENCH_GLOB )then
            t_tot = tic()
            t     = t_tot
        endif
        !$ call omp_set_num_threads(nthr)
        nthr_glob = nthr
        call set_b_p_ptrs2D(params, build)
        call ensure_field_halves(build%spproj, params%oritype, params%projfile)
        if( params%startit == 1 .and. params%which_iter == params%startit )then
            call build%spproj_field%clean_entry('updatecnt', 'sampled')
        endif
        ! Mirror the 3D workflow: sampled-update is active from the first stage onward.
        ! In probabilistic mode the sampled subset is reused within the current iteration
        ! by prob_tab2D/refine2D_exec; later iterations redraw it, except on a fresh start
        ! without fillin (startit=1), which keeps the subset it drew first.
        call sample_ptcls4update2D(params, build, [params%fromp,params%top], params%l_update_frac, nptcls, pinds)
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
            if( .not. present(xprob_tab2D) ) THROW_HARD('prob_align2D in shared memory needs the table commander')
            call xprob_tab2D%execute(cline_prob_tab)
        else
            call qenv%new(params, params%nparts, nptcls=params%nptcls)
            call cline_prob_tab%gen_job_descr(job_descr)
            call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        endif
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
        call eulprob_obj_glob%ref_assign
        if( L_BENCH_GLOB )then
            rt_assign = toc(t)
            t         = tic()
        endif
        ! write assignment to file
        fname = string(ASSIGNMENT_FBODY)//'.dat'
        call eulprob_obj_glob%write_assignment(fname)
        write(logfhandle,'(A,I0,A,I0,A)') '>>> PROB_ALIGN2D: ', nptcls, ' particles over ', params%nparts,&
            &' part(s) assigned'
        call flush(logfhandle)
        if( L_BENCH_GLOB ) rt_write = toc(t)
        assign_profile = eulprob_obj_glob%get_assign_profile()
        ! cleanup
        call eulprob_obj_glob%kill
        call cline_prob_tab%kill
        call qenv%kill
        call job_descr%kill
        if( L_BENCH_GLOB ) call write_prob_align_bench(params, 'PROB_ALIGN2D', nptcls, rt_setup_here, rt_sample, rt_sched, &
            &rt_read, rt_assign, rt_write, toc(t_tot), assign2D=assign_profile)
        call qsys_cleanup(params)
    end subroutine prob_align2D_exec

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

    !> The bench record of the prob_align driver, all wall time, in family PROB_ALIGN (prob), PROB_ALIGN_NEIGH
    !! (prob_neigh) or PROB_ALIGN2D; rt_sched runs from the launch of the parts to the last part finished (the
    !! in-process table fill in shared memory); assign2D splits the 2D reads and assignment.
    subroutine write_prob_align_bench( params, family, nptcls, rt_setup, rt_sample, rt_sched, rt_read, rt_assign, &
            &rt_write, rt_tot, assign2D )
        use simple_bench_io,       only: open_bench_report
        use simple_eul_prob_tab2D, only: prob_assign_profile2D
        class(parameters),                     intent(in) :: params
        character(len=*),                      intent(in) :: family
        integer,                               intent(in) :: nptcls
        real(timer_int_kind),                  intent(in) :: rt_setup, rt_sample, rt_sched, rt_read, rt_assign, rt_write, rt_tot
        type(prob_assign_profile2D), optional, intent(in) :: assign2D
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
        if( present(assign2D) )then
            ! parts of the table reads and of ref_assign, not added to '% accounted for'
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** READ AND ASSIGNMENT DETAIL TIMINGS (s) ***'
            write(fnr,'(a,1x,f0.3)') 'prob_align2D read, candidate counts  :', assign2D%read_count
            write(fnr,'(a,1x,f0.3)') 'prob_align2D read, candidates        :', assign2D%read_fill
            write(fnr,'(a,1x,f0.3)') 'prob_align2D class lists, counts     :', assign2D%lists_count
            write(fnr,'(a,1x,f0.3)') 'prob_align2D class lists, fill       :', assign2D%lists_fill
            write(fnr,'(a,1x,f0.3)') 'prob_align2D class lists, sort       :', assign2D%lists_sort
            write(fnr,'(a,1x,f0.3)') 'prob_align2D frontier loop           :', assign2D%frontier
            write(fnr,'(a,1x,f0.3)') 'prob_align2D leftover particles      :', assign2D%leftover
        endif
        call fclose(fnr)
    end subroutine write_prob_align_bench

end module simple_prob_align_strategy
