!@descr: for sigma2 calculations in objfun=euclid 2D and 3D refinement
module simple_commanders_euclid
use, intrinsic :: iso_fortran_env, only: int64
use simple_commanders_api
use simple_sigma2_state, only: sigma2_state_candidate_path, sigma2_state_assemble_candidate, &
    &sigma2_state_range_path, sigma2_state_next_generation
use simple_sigma2_state_file, only: sigma2_state_header, sigma2_state_read_header, sigma2_state_publish
implicit none
private :: run_commit_sigma2
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_calc_pspec
  contains
    procedure :: execute      => exec_calc_pspec
end type commander_calc_pspec

type, extends(commander_base) :: commander_commit_sigma2
  contains
    procedure :: execute              => exec_commit_sigma2
    procedure :: execute_from_project => exec_commit_sigma2_from_project
end type commander_commit_sigma2

contains

    subroutine exec_calc_pspec( self,cline )
        use simple_calc_pspec_strategy, only: calc_pspec_strategy, create_calc_pspec_strategy
        class(commander_calc_pspec), intent(inout) :: self
        class(cmdline), intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        class(calc_pspec_strategy), allocatable :: strategy
        call cline%set('stream','no')
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',  'yes')
        if( .not. cline%defined('oritype') ) call cline%set('oritype','ptcl3D')
        ! Create and run strategy (shared-memory/worker vs distributed master)
        strategy = create_calc_pspec_strategy(cline)
        call strategy%initialize(params, build, cline)
        call strategy%execute(params, build, cline)
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params, build, cline)
        ! Unified termination point
        call simple_end(strategy%end_msg, print_simple=strategy%end_print_simple)
        if( allocated(strategy) ) deallocate(strategy)
    end subroutine exec_calc_pspec

    subroutine exec_commit_sigma2( self, cline )
        class(commander_commit_sigma2), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        call run_commit_sigma2(cline)
    end subroutine exec_commit_sigma2

    !> the commit with the caller's project of projfile in memory, as the file holds it: copied, not read again
    subroutine exec_commit_sigma2_from_project( self, cline, spproj )
        class(commander_commit_sigma2), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        class(sp_project),              intent(in)    :: spproj
        call run_commit_sigma2(cline, spproj)
    end subroutine exec_commit_sigma2_from_project

    subroutine run_commit_sigma2( cline, spproj_src )
        class(cmdline),              intent(inout) :: cline
        class(sp_project), optional, intent(in)    :: spproj_src
        type(parameters)      :: params
        type(builder)         :: build
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',  'no')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.,spproj_src=spproj_src)
        call consolidate_canonical
        ! cleanup
        call build%kill_general_tbox
        call simple_touch('COMMIT_SIGMA2_FINISHED')
        call simple_end('**** SIMPLE_COMMIT_SIGMA2 NORMAL STOP ****', print_simple=.false.)
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
            call open_bench_report(fnr, 'SIGMA2_COMMIT', iter=params%which_iter)
            write(fnr,'(a)') '*** BENCHMARK CONTEXT ***'
            write(fnr,'(a,a)')  'sigma2 commit program               : ', params%prg%to_char()
            write(fnr,'(a,i0)') 'sigma2 commit which_iter            : ', params%which_iter
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

    end subroutine run_commit_sigma2

end module simple_commanders_euclid
