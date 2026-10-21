!@descr: class-average assembly of refine2D: the workers' partial class sums reduced, restored and written
module simple_cavgassemble_strategy
use simple_core_module_api
use simple_builder,         only: builder
use simple_cmdline,         only: cmdline
use simple_parameters,      only: parameters
use simple_classaverager,   only: cavger_new, cavger_assemble_sums_from_parts, cavger_gen2Dclassdoc, cavger_write_all, &
    &cavger_kill
use simple_refine2D_fnames, only: refine2D_cavgs_fname, refine2D_cavgs_half_fname, refine2D_cavgs_iter
implicit none
private
public :: run_cavgassemble, cavgassemble_exec
#include "simple_local_flags.inc"

contains

    !> class-average assembly on its own parameters and builder, built from cline
    subroutine run_cavgassemble( cline )
        class(cmdline), intent(inout) :: cline
        type(parameters)        :: params
        type(builder)           :: build
        integer(timer_int_kind) :: t_init
        real(timer_int_kind)    :: rt_init
        rt_init = 0.
        if( L_BENCH_GLOB ) t_init = tic()
        call cline%set('oritype', 'ptcl2D')
        call build%init_params_and_build_strategy2D_tbox(cline, params, wthreads=.true.)
        if( L_BENCH_GLOB ) rt_init = toc(t_init)
        call cavgassemble_exec(params, build, cline, rt_init, params%nthr)
        ! multiple fields updated, do a full write
        call build%spproj%write(params%projfile)
        ! end gracefully
        call build%spproj%kill
        call build%kill_general_tbox
        call build%kill_strategy2D_tbox
        call simple_end('**** SIMPLE_CAVGASSEMBLE NORMAL STOP ****', print_simple=.false.)
        ! indicate completion (when run in a qsys env)
        call simple_touch('CAVGASSEMBLE_FINISHED')
    end subroutine run_cavgassemble

    !> class-average assembly on the caller's parameters and builder at nthr threads (the caller restores its own): the
    !! class averages named by which_iter (or the start references) reduced from the workers' partial sums, restored and
    !! written, the class document and the cls2D star file; the segments it changes (cls2D, cls3D, out and the particles'
    !! res) updated in the caller's project in memory, for the caller to write; rt_init: the caller's set-up time, for
    !! the bench record; prev_frcs: restoration starts from the previous iteration's class FRCs (params%frcs) when the
    !! file exists, as shared-memory refine2D does, instead of from none
    subroutine cavgassemble_exec( params, build, cline, rt_init, nthr, prev_frcs )
        use simple_bench_io,     only: open_bench_report
        use simple_starproject,  only: starproject
        use simple_stream_utils, only: terminate_stream
        !$ use omp_lib,          only: omp_set_num_threads
        class(parameters),    intent(inout) :: params
        class(builder),       intent(inout) :: build
        class(cmdline),       intent(inout) :: cline
        real(timer_int_kind), intent(in)    :: rt_init
        integer,              intent(in)    :: nthr
        logical, optional,    intent(in)    :: prev_frcs
        type(starproject)  :: starproj
        real, allocatable  :: states(:)
        integer            :: iter, io_stat, fnr
        integer(timer_int_kind) :: t_tot, t_phase
        real(timer_int_kind)    :: rt_reduce_partials, rt_classdoc, rt_write_cavgs
        real(timer_int_kind)    :: rt_export_star, rt_project_update, rt_cleanup, rt_tot
        if( L_BENCH_GLOB )then
            rt_reduce_partials= 0.
            rt_classdoc       = 0.
            rt_write_cavgs    = 0.
            rt_export_star    = 0.
            rt_project_update = 0.
            rt_cleanup        = 0.
            t_tot             = tic()
        endif
        !$ call omp_set_num_threads(nthr)
        nthr_glob = nthr
        if( cline%defined('which_iter') )then
            params%refs      = refine2D_cavgs_fname(params%which_iter, params%ext%to_char())
            params%refs_even = refine2D_cavgs_half_fname(params%which_iter, 'even', params%ext%to_char())
            params%refs_odd  = refine2D_cavgs_half_fname(params%which_iter, 'odd', params%ext%to_char())
        else if( .not. cline%defined('refs') )then
            params%refs      = 'start2Drefs'//params%ext%to_char()
            params%refs_even = 'start2Drefs_even'//params%ext%to_char()
            params%refs_odd  = 'start2Drefs_odd'//params%ext%to_char()
        endif
        ! restoration reads the class FRCs for the even/odd low-resolution averaging before it stores this assembly's:
        ! refine2D's previous iteration's (as its shared-memory matcher reads them), otherwise none
        call build%clsfrcs%new(params%ncls, params%box_crop, params%smpd_crop, params%nstates)
        if( present(prev_frcs) )then
            if( prev_frcs )then
                if( file_exists(params%frcs) ) call build%clsfrcs%read(params%frcs)
            endif
        endif
        if( L_BENCH_GLOB ) t_phase = tic()
        call cavger_new(params, build)
        call cavger_assemble_sums_from_parts
        if( L_BENCH_GLOB ) rt_reduce_partials = toc(t_phase)
        if( L_BENCH_GLOB ) t_phase = tic()
        call cavger_gen2Dclassdoc
        if( L_BENCH_GLOB ) rt_classdoc = toc(t_phase)
        call terminate_stream(params, 'SIMPLE_CAVGASSEMBLE HARD STOP')
        if( L_BENCH_GLOB ) t_phase = tic()
        call cavger_write_all(params%refs, params%refs_even, params%refs_odd)
        call cavger_kill
        if( L_BENCH_GLOB ) rt_write_cavgs = toc(t_phase)
        ! get iteration from which_iter else from refs filename and write cavgs starfile
        if( cline%defined('which_iter') ) then
            if( L_BENCH_GLOB ) t_phase = tic()
            call starproj%export_cls2D(build%spproj, params%which_iter)
            if( L_BENCH_GLOB ) rt_export_star = toc(t_phase)
            iter = params%which_iter
        else if( cline%defined('refs') .and. params%refs%substr_ind(CAVGS_ITER_FBODY) > 0 ) then
            iter = refine2D_cavgs_iter(params%refs, params%ext, io_stat)
            if( L_BENCH_GLOB ) t_phase = tic()
            call starproj%export_cls2D(build%spproj, iter)
            if( L_BENCH_GLOB ) rt_export_star = toc(t_phase)
        else
            iter = max(1, params%which_iter)
        end if
        ! updates project
        ! cls2D and state congruent cls3D
        if( L_BENCH_GLOB ) t_phase = tic()
        call build%spproj%os_cls3D%new(params%ncls, is_ptcl=.false.)
        states = build%spproj%os_cls2D%get_all('state')
        call build%spproj%os_cls3D%set_all('state',states)
        deallocate(states)
        call build%spproj%add_frcs2os_out( string(FRCS_FILE), 'frc2D')
        call build%spproj%add_cavgs2os_out(params%refs, build%spproj%get_smpd(), imgkind='cavg')
        if( L_BENCH_GLOB ) rt_project_update = toc(t_phase)
        call starproj%kill
        if( L_BENCH_GLOB )then
            rt_tot = toc(t_tot) + rt_init
            call open_bench_report(fnr, 'CAVGASSEMBLE', iter=iter)
            write(fnr,'(a)') '*** BENCHMARK CONTEXT ***'
            write(fnr,'(a,a)')  'cavgassemble assembly mode          : class'
            write(fnr,'(a,a)')  'cavgassemble refs                   : ', params%refs%to_char()
            write(fnr,'(a,i0)') 'cavgassemble nclasses               : ', params%ncls
            write(fnr,'(a,i0)') 'cavgassemble nparts                 : ', params%nparts
            write(fnr,'(a,i0)') 'cavgassemble kfrom                  : ', params%kfromto(1)
            write(fnr,'(a,i0)') 'cavgassemble kto                    : ', params%kfromto(2)
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f0.2)') 'cavgassemble setup/init             :', rt_init
            write(fnr,'(a,1x,f0.2)') 'cavgassemble reduce partial sums    :', rt_reduce_partials
            write(fnr,'(a,1x,f0.2)') 'cavgassemble class document         :', rt_classdoc
            write(fnr,'(a,1x,f0.2)') 'cavgassemble write class averages   :', rt_write_cavgs
            write(fnr,'(a,1x,f0.2)') 'cavgassemble export cls2D star      :', rt_export_star
            write(fnr,'(a,1x,f0.2)') 'cavgassemble project update         :', rt_project_update
            write(fnr,'(a,1x,f0.2)') 'cavgassemble cleanup                :', rt_cleanup
            write(fnr,'(a,1x,f0.2)') 'cavgassemble total time             :', rt_tot
            write(fnr,'(a,1x,f0.2)') 'cavgassemble % accounted for        :', &
                &((rt_init + rt_reduce_partials + rt_classdoc + rt_write_cavgs + &
                &  rt_export_star + rt_project_update + rt_cleanup) / rt_tot) * 100.
            call fclose(fnr)
        endif
    end subroutine cavgassemble_exec

end module simple_cavgassemble_strategy
