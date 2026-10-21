!@descr: canonical sigma2 store steps shared by the refine2D and refine3D strategies
module simple_strategy_sigma2
use, intrinsic :: iso_fortran_env, only: int64
use simple_defs,              only: STDLEN
use simple_error,             only: simple_exception
use simple_string,            only: string
use simple_parameters,        only: parameters
use simple_builder,           only: builder
use simple_sigma2_state,      only: sigma2_state_prepare_update
use simple_sigma2_state_file, only: sigma2_state_header, sigma2_state_read_header
implicit none

public :: prepare_canonical_sigma_update, read_sigma2_bench_counts, run_commit_sigma2, commit_sigma2_exec
private
#include "simple_local_flags.inc"

contains

    !> before the parts write the next generation: the committed store validated, nothing left by another
    !! transaction of that generation
    subroutine prepare_canonical_sigma_update( params, build )
        type(parameters), intent(in)    :: params
        type(builder),    intent(inout) :: build
        type(string)          :: state_path
        integer               :: status
        logical               :: found
        character(len=STDLEN) :: message
        call build%spproj%get_sigma2_state_path(state_path, found)
        if( .not. found ) THROW_HARD('particle project has no canonical sigma2 state path')
        call sigma2_state_prepare_update(state_path%to_char(), params%nparts, params%numlen, status, message)
        if( status /= 0 ) THROW_HARD(trim(message))
        call state_path%kill
    end subroutine prepare_canonical_sigma_update

    !> the canonical store's size for a bench record: particles, shells, bytes and the rows updated this iteration;
    !! unchanged when the project has no store
    subroutine read_sigma2_bench_counts( build, nptcls, nshells, nbytes, nrows_updated )
        type(builder),  intent(inout) :: build
        integer,        intent(inout) :: nptcls, nshells, nrows_updated
        integer(int64), intent(inout) :: nbytes
        type(sigma2_state_header) :: header
        type(string)              :: state_path
        integer, allocatable      :: rep_rows(:), smp_rows(:), rec_rows(:)
        integer                   :: status
        logical                   :: found
        character(len=STDLEN)     :: message
        call build%spproj%get_sigma2_state_path(state_path, found)
        if( .not. found ) return
        call sigma2_state_read_header(state_path%to_char(), header, status, message)
        if( status /= 0 ) THROW_HARD(trim(message))
        nptcls  = header%nptcls
        nshells = header%kto - header%kfrom + 1
        nbytes  = header%file_bytes
        call build%spproj_field%get_update_rows(rep_rows, smp_rows, rec_rows)
        nrows_updated = size(smp_rows)
        deallocate(rep_rows, smp_rows, rec_rows)
        call state_path%kill
    end subroutine read_sigma2_bench_counts

    !> the sigma2 commit on its own parameters and builder, built from cline: the workers' ranges of the next generation
    !! assembled into the canonical store and published
    subroutine run_commit_sigma2( cline )
        use simple_core_module_api
        use simple_cmdline, only: cmdline
        class(cmdline), intent(inout) :: cline
        type(parameters)      :: params
        type(builder)         :: build
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',  'no')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        call commit_sigma2_exec(params, build, params%which_iter, params%nthr)
        ! cleanup
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_COMMIT_SIGMA2 NORMAL STOP ****', print_simple=.false.)
    end subroutine run_commit_sigma2

    !> the sigma2 commit on the caller's parameters and builder (the particles' states and halves and the canonical
    !! sigma2 path of its project): nthr threads (the caller restores its own); which_iter: the iteration the commit
    !! reports (refine2D passes the iteration the committed store feeds, refine3D the matcher iteration)
    subroutine commit_sigma2_exec( params, build, which_iter, nthr )
        use simple_core_module_api
        use simple_sigma2_state,      only: sigma2_state_candidate_path, sigma2_state_assemble_candidate, &
            &sigma2_state_range_path, sigma2_state_next_generation
        use simple_sigma2_state_file, only: sigma2_state_publish
        !$ use omp_lib,               only: omp_set_num_threads
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: which_iter, nthr
        !$ call omp_set_num_threads(nthr)
        nthr_glob = nthr
        call consolidate_canonical
        call simple_touch('COMMIT_SIGMA2_FINISHED')
    contains

        !> the streaming commit: a fresh candidate assembled from the workers' ranges in one pass that
        !! also forms the global spectra, then published; the ranges are deleted only after publication
        subroutine consolidate_canonical
            type(string) :: state_path, candidate_path
            type(string), allocatable :: range_paths(:)
            type(sigma2_state_header) :: header
            logical, allocatable :: scheduled(:), active(:)
            integer, allocatable :: eo_ids(:)
            integer(int64) :: next_gen
            integer(timer_int_kind) :: t_assemble, t_publish, t_tot
            real(timer_int_kind)    :: rt_assemble, rt_publish, rt_tot
            integer :: ipart, iptcl, status
            logical :: found
            character(len=STDLEN) :: message
            if( L_BENCH_GLOB ) t_tot = tic()
            call build%spproj%get_sigma2_state_path(state_path, found)
            if( .not. found ) THROW_HARD('particle project has no canonical sigma2 state path')
            ! the transaction the workers wrote their ranges for: the generation this commit makes
            call sigma2_state_next_generation(state_path%to_char(), next_gen, status, message)
            if( status /= 0 ) THROW_HARD(trim(message))
            candidate_path = sigma2_state_candidate_path(state_path%to_char(), next_gen)
            allocate(range_paths(params%nparts))
            do ipart = 1, params%nparts
                range_paths(ipart) = sigma2_state_range_path(state_path%to_char(), next_gen, ipart, params%numlen)
            enddo
            allocate(scheduled(params%nptcls), source=.true.)
            allocate(active(params%nptcls), eo_ids(params%nptcls))
            do iptcl = 1, params%nptcls
                active(iptcl) = build%spproj_field%get_state(iptcl) > 0
                eo_ids(iptcl) = build%spproj_field%get_eo(iptcl)
            enddo
            if( L_BENCH_GLOB ) t_assemble = tic()
            call sigma2_state_assemble_candidate(state_path%to_char(), candidate_path%to_char(), range_paths, &
                &scheduled, active, eo_ids, status, message)
            if( status /= 0 ) THROW_HARD(trim(message))
            if( L_BENCH_GLOB )then
                rt_assemble = toc(t_assemble)
                t_publish   = tic()
            endif
            call sigma2_state_publish(candidate_path%to_char(), state_path%to_char(), status, message)
            if( status /= 0 ) THROW_HARD(trim(message))
            if( L_BENCH_GLOB ) rt_publish = toc(t_publish)
            do ipart = 1, params%nparts
                call del_file(range_paths(ipart))
            enddo
            if( L_BENCH_GLOB )then
                rt_tot = toc(t_tot)
                call sigma2_state_read_header(state_path%to_char(), header, status, message)
                if( status /= 0 ) THROW_HARD(trim(message))
                call write_commit_bench(header, count(scheduled), rt_assemble, rt_publish, rt_tot)
            endif
            call range_paths(:)%kill
            deallocate(range_paths, scheduled, active, eo_ids)
            call state_path%kill
            call candidate_path%kill
        end subroutine consolidate_canonical

        !> The SIGMA2_COMMIT bench record, at the which_iter the commit is given (refine2D passes the
        !! iteration the committed store feeds, refine3D the matcher iteration)
        subroutine write_commit_bench( header, nrows_merged, rt_assemble, rt_publish, rt_tot )
            use simple_bench_io, only: open_bench_report
            type(sigma2_state_header), intent(in) :: header
            integer,                   intent(in) :: nrows_merged
            real(timer_int_kind),      intent(in) :: rt_assemble, rt_publish, rt_tot
            real(timer_int_kind) :: rt_accounted
            integer :: fnr
            rt_accounted = 0.
            if( rt_tot > TINY ) rt_accounted = 100. * (rt_assemble + rt_publish) / rt_tot
            call open_bench_report(fnr, 'SIGMA2_COMMIT', iter=which_iter)
            write(fnr,'(a)') '*** BENCHMARK CONTEXT ***'
            write(fnr,'(a,a)')  'sigma2 commit program               : ', params%prg%to_char()
            write(fnr,'(a,i0)') 'sigma2 commit which_iter            : ', which_iter
            write(fnr,'(a,i0)') 'sigma2 commit nparts (range files)  : ', params%nparts
            write(fnr,'(a,i0)') 'sigma2 store particles              : ', header%nptcls
            write(fnr,'(a,i0)') 'sigma2 store shells                 : ', header%kto - header%kfrom + 1
            write(fnr,'(a,i0)') 'sigma2 store bytes                  : ', header%file_bytes
            write(fnr,'(a,i0)') 'sigma2 rows merged from ranges      : ', nrows_merged
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f0.3)') 'sigma2 commit streaming assembly    :', rt_assemble
            write(fnr,'(a,1x,f0.3)') 'sigma2 commit publication           :', rt_publish
            write(fnr,'(a,1x,f0.3)') 'sigma2 commit total time            :', rt_tot
            write(fnr,'(a,1x,f0.2)') 'sigma2 commit % accounted for       :', rt_accounted
            call fclose(fnr)
        end subroutine write_commit_bench

    end subroutine commit_sigma2_exec

end module simple_strategy_sigma2
