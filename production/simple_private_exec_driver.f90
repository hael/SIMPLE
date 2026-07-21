!@descr: callable driver for simple_private_exec private commander dispatch
module simple_private_exec_driver
use simple_private_exec_api
use simple_syslib, only: redirect_stdout_stderr, restore_stdout_stderr
use simple_memory_monitor, only: mem_monitor_init, mem_monitor_finish
implicit none
private
public :: run_private_exec_from_command_line, run_private_exec_line, run_coarray_direct
#include "simple_local_flags.inc"

logical, save :: private_exec_ui_ready = .false.

contains

    subroutine run_private_exec_from_command_line
        character(len=STDLEN)   :: xarg, prg
        type(cmdline)           :: cline
        integer                 :: cmdstat, cmdlen, pos
        integer(timer_int_kind) :: t0
        real(timer_int_kind)    :: rt_exec
        logical                 :: l_silent
        t0 = tic()
        call get_command_argument(1, xarg, cmdlen, cmdstat)
        if( cmdstat /= 0 ) call cmdline_err(cmdstat, cmdlen, xarg, 0)
        if( trim(xarg) == '--coarray' ) THROW_HARD('run_private_exec_from_command_line cannot dispatch --coarray')
        pos = index(xarg, '=')
        call cmdline_err(cmdstat, cmdlen, xarg, pos)
        prg = xarg(pos+1:)
        call ensure_private_exec_ui
        call cline%parse_private
        call mem_monitor_init(cline, 'simple_private_exec:'//trim(prg))
        call print_slurm_env
        call dispatch_private_prg(prg, cline, l_silent)
        rt_exec = toc(t0)
        if( .not. l_silent ) call simple_print_timer(rt_exec)
        call mem_monitor_finish
        call cline%kill
    end subroutine run_private_exec_from_command_line

    subroutine run_private_exec_line( job_arg_line, output_file )
        character(len=*),           intent(in) :: job_arg_line
        character(len=*), optional, intent(in) :: output_file
        if( present(output_file) ) call redirect_stdout_stderr(output_file)
        call run_private_exec_line_unredirected(job_arg_line)
        if( present(output_file) ) call restore_stdout_stderr
    end subroutine run_private_exec_line

    subroutine run_private_exec_line_unredirected( job_arg_line )
        character(len=*), intent(in) :: job_arg_line
        character(len=STDLEN)        :: prg
        type(cmdline)                :: cline
        integer(timer_int_kind)      :: t0
        real(timer_int_kind)         :: rt_exec
        logical                      :: l_silent
        t0 = tic()
        call extract_prg_from_line(job_arg_line, prg)
        call ensure_private_exec_ui
        call cline%parse_private_line(job_arg_line)
        call mem_monitor_init(cline, 'simple_private_exec:'//trim(prg))
        call print_slurm_env
        call dispatch_private_prg(prg, cline, l_silent)
        rt_exec = toc(t0)
        if( .not. l_silent ) call simple_print_timer(rt_exec)
        call mem_monitor_finish
        call cline%kill
    end subroutine run_private_exec_line_unredirected

    !> Coarray backend entry point. Each image executes its assigned private
    !! command lines in-process while redirecting stdout/stderr per partition.
    subroutine run_coarray_direct
#ifdef USE_COARRAYS
        character(len=XLONGSTRLEN)              :: part_log
        character(len=XLONGSTRLEN), allocatable :: part_job_args(:)
        logical, allocatable                    :: part_finished(:)
        integer :: nargs, nparts_run, from_part, to_part, part_log_width
        integer :: iarg, ipart, part_job_ind, argstat, arglen
        nargs = command_argument_count()
        if( nargs < 4 )then
            THROW_HARD('Usage: simple_private_exec --coarray <from_part> <to_part> <part_job_args...>')
        endif
        call parse_coarray_int_arg(2, from_part, 'from_part')
        call parse_coarray_int_arg(3, to_part,   'to_part')
        if( from_part < 1 .or. to_part < from_part )then
            THROW_HARD('Invalid coarray partition range')
        endif
        nparts_run = to_part - from_part + 1
        if( nargs /= 3 + nparts_run )then
            THROW_HARD('simple_private_exec --coarray received wrong number of partition job arguments')
        endif
        part_log_width = len(int2str(to_part))
        allocate(part_job_args(nparts_run), part_finished(from_part:to_part))
        part_finished = .false.
        do iarg = 1, nparts_run
            call get_command_argument(3 + iarg, part_job_args(iarg), arglen, argstat)
            if( argstat /= 0 )then
                THROW_HARD('could not read coarray partition job argument for part '//int2str(from_part + iarg - 1))
            endif
            if( len_trim(part_job_args(iarg)) == 0 )then
                THROW_HARD('empty coarray partition job argument for part '//int2str(from_part + iarg - 1))
            endif
        end do
        do ipart = from_part + this_image() - 1, to_part, num_images()
            part_job_ind = ipart - from_part + 1
            part_log     = 'simple_private_exec_coarray_part_'//int2str_pad(ipart, part_log_width)//'.out'
            call run_private_exec_line(trim(part_job_args(part_job_ind)), trim(part_log))
            call declare_coarray_part_finished(ipart, from_part, to_part, part_finished)
        end do
        call sync_coarray_part_completion(from_part, to_part, part_finished)
#else
        THROW_HARD('simple_private_exec --coarray requires a USE_COARRAYS build')
#endif
    end subroutine run_coarray_direct

#ifdef USE_COARRAYS
    subroutine declare_coarray_part_finished( ipart, from_part, to_part, part_finished )
        integer, intent(in)    :: ipart, from_part, to_part
        logical, intent(inout) :: part_finished(from_part:to_part)
        if( ipart < from_part .or. ipart > to_part )then
            THROW_HARD('coarray completion declaration outside partition range')
        endif
        part_finished(ipart) = .true.
    end subroutine declare_coarray_part_finished

    subroutine sync_coarray_part_completion( from_part, to_part, part_finished )
        integer, intent(in) :: from_part, to_part
        logical, intent(in) :: part_finished(from_part:to_part)
        integer :: ipart, syncstat
        do ipart = from_part + this_image() - 1, to_part, num_images()
            if( .not. part_finished(ipart) )then
                THROW_HARD('coarray image did not declare completion for part '//int2str(ipart))
            endif
        end do
        sync all(stat=syncstat)
        if( syncstat /= 0 ) THROW_HARD('simple_private_exec coarray sync failed with stat='//int2str(syncstat))
    end subroutine sync_coarray_part_completion

    subroutine parse_coarray_int_arg( iarg, val, name )
        integer,          intent(in)  :: iarg
        integer,          intent(out) :: val
        character(len=*), intent(in)  :: name
        character(len=XLONGSTRLEN)    :: raw
        integer                       :: parse_stat
        call get_command_argument(iarg, raw)
        read(raw,*,iostat=parse_stat) val
        if( parse_stat /= 0 )then
            THROW_HARD('Invalid coarray '//trim(name)//' argument: '//trim(raw))
        endif
    end subroutine parse_coarray_int_arg

#endif

    subroutine ensure_private_exec_ui
        if( private_exec_ui_ready ) return
        call make_ui
        call make_private_ui
        private_exec_ui_ready = .true.
    end subroutine ensure_private_exec_ui

    subroutine extract_prg_from_line( job_arg_line, prg )
        character(len=*),     intent(in)  :: job_arg_line
        character(len=STDLEN), intent(out) :: prg
        character(len=XLONGSTRLEN) :: line, first_arg
        integer :: pos_space, pos
        line = adjustl(trim(job_arg_line))
        if( len_trim(line) == 0 ) THROW_HARD('empty private command line')
        pos_space = index(line, ' ')
        if( pos_space > 0 )then
            first_arg = line(:pos_space-1)
        else
            first_arg = line
        endif
        pos = index(first_arg, '=')
        call cmdline_err(0, len_trim(first_arg), first_arg, pos)
        if( trim(first_arg(:pos-1)) /= 'prg' ) THROW_HARD('private command line must start with prg=')
        prg = first_arg(pos+1:)
    end subroutine extract_prg_from_line

    subroutine dispatch_private_prg( prg, cline, l_silent )
        character(len=*), intent(in)    :: prg
        class(cmdline),   intent(inout) :: cline
        logical,          intent(out)   :: l_silent
        ! PRE-PROCESSING PROGRAMS
        type(commander_preprocess)              :: xpreprocess
        type(commander_extract)                 :: xextract
        type(commander_reextract)               :: xreextract
        type(commander_motion_correct)          :: xmotion_correct
        type(commander_gen_pspecs_and_thumbs)   :: xgen_pspecs_and_thumbs
        type(commander_ctf_estimate)            :: xctf_estimate
        type(commander_pick_extract)            :: xpick_extract
        type(commander_pick)                    :: xpick
        type(commander_shape_rank_cavgs)        :: xshape_rank_cavgs
        type(commander_make_pickrefs)           :: xmake_pickrefs
        ! CLUSTER2D PROGRAMS
        type(commander_make_cavgs)              :: xmake_cavgs
        type(commander_cluster2D_distr_worker)  :: xcluster2D_worker
        type(commander_cluster2D)               :: xcluster2D
        type(commander_cavgassemble)            :: xcavgassemble
        type(commander_rank_cavgs)              :: xrank_cavgs
        type(commander_export_cavgs)            :: xexport_cavgs
        type(commander_cls_split)               :: xcls_split
        type(commander_denoise_project)         :: xdenoise_project
        ! REFINE3D PROGRAMS
        type(commander_refine3D_distr_worker)   :: xrefine3D_worker
        type(commander_calc_pspec)              :: xcalc_pspec
        type(commander_calc_pspec_assemble)     :: xcalc_pspec_assemble
        type(commander_calc_group_sigmas)       :: xcalc_group_sigmas
        type(commander_prob_tab)                :: xprob_tab
        type(commander_prob_tab_neigh)          :: xprob_tab_neigh
        type(commander_prob_align)              :: xprob_align
        type(commander_prob_align_neigh)        :: xprob_align_neigh
        type(commander_prob_tab2D)              :: xprob_tab2D
        type(commander_prob_align2D)            :: xprob_align2D
        type(commander_flex_eigenvol_reconstruct) :: xflex_eigenvol_reconstruct
        ! RECONSTRUCTION PROGRAMS
        type(commander_volassemble)             :: xvolassemble
        type(commander_rec3D_worker)            :: xrec3D
        ! CHECKER PROGRAMS
        type(commander_check_box)               :: xcheck_box
        type(commander_check_nptcls)            :: xcheck_nptcls
        type(commander_check_stoch_update)      :: xcheck_stoch_update
        type(commander_check_update_frac)       :: xcheck_update_frac
        ! VOLOPS PROGRAMS
        type(commander_postprocess)             :: xpostprocess
        type(commander_automask)                :: xautomask
        ! GENERAL IMAGE PROCESSING PROGRAMS
        type(commander_scale)                   :: xscale
        type(commander_binarize)                :: xbinarize
        type(commander_ppca_denoise)            :: xppca_denoise
        ! MISCELLANOUS PROGRAMS
        type(commander_aggregate_chunks)        :: xaggregate_chunks
        type(commander_fractionate_movies)      :: xfractionate_movies
        type(commander_kstest)                  :: xkstst
        type(commander_pearsn)                  :: xpearsn
        ! ORIENTATION DATA MANAGEMENT PROGRAMS
        type(commander_rotmats2oris)            :: xrotmats2oris
        type(commander_print_project_vals)      :: xprint_project_vals
        ! DATA MANAGEMENT PROGRAMS
        type(commander_prune_project)           :: xprune_project
        type(commander_scale_project)           :: xscale_project
        ! TIME-SERIES ANALYSIS PROGRAMS
        type(commander_track_particles)         :: xtrack_particles
        type(commander_tseries_motion_correct)  :: xtseries_mcorr
        ! PARALLEL PROCESSING PROGRAMS
        type(commander_split)                   :: xsplit
        l_silent = .false.
        select case(trim(prg))
            ! PRIVATE UTILITY PROGRAMS
            case( 'print_ui_json' )
                call print_ui_json
                l_silent = .true.
            case( 'write_ui_json' )
                call write_ui_json
            case( 'print_sym_subgroups' )
                call print_subgroups
                l_silent = .true.
            case( 'print_ui_stream' )
                call print_stream_ui_json
                l_silent = .true.
            ! PRE-PROCESSING PROGRAMS
            case( 'preprocess' )
                call xpreprocess%execute(cline)
            case( 'extract' )
                call xextract%execute(cline)
            case( 'reextract' )
                call xreextract%execute(cline)
            case( 'motion_correct' )
                call xmotion_correct%execute(cline)
            case( 'gen_pspecs_and_thumbs' )
                call xgen_pspecs_and_thumbs%execute(cline)
            case( 'ctf_estimate' )
                call xctf_estimate%execute(cline)
            case( 'pick_extract' )
                call xpick_extract%execute(cline)
            case( 'pick' )
                call xpick%execute(cline)
            case( 'shape_rank_cavgs' )
                call xshape_rank_cavgs%execute(cline)
            case( 'make_pickrefs' )
                call xmake_pickrefs%execute(cline)
            ! CLUSTER2D PROGRAMS
            case( 'make_cavgs' )
                call xmake_cavgs%execute(cline)
            case( 'cluster2D' )
                call xcluster2D_worker%execute(cline)
            case( 'cluster2D_distr' )
                call xcluster2D%execute(cline)
            case( 'cavgassemble' )
                call xcavgassemble%execute(cline)
            case( 'rank_cavgs' )
                call xrank_cavgs%execute(cline)
            case( 'export_cavgs' )
                call xexport_cavgs%execute(cline)
            case( 'cls_split' )
                call xcls_split%execute(cline)
            case( 'denoise_project' )
                call xdenoise_project%execute(cline)
            ! REFINE3D PROGRAMS
            case( 'refine3D' )
                call xrefine3D_worker%execute(cline)
            case( 'calc_pspec' )
                call xcalc_pspec%execute(cline)
            case( 'calc_pspec_assemble' )
                call xcalc_pspec_assemble%execute(cline)
            case( 'calc_group_sigmas' )
                call xcalc_group_sigmas%execute(cline)
            case( 'prob_tab' )
                call xprob_tab%execute(cline)
            case( 'prob_tab_neigh' )
                call xprob_tab_neigh%execute(cline)
            case( 'prob_align' )
                call xprob_align%execute(cline)
            case( 'prob_align_neigh' )
                call xprob_align_neigh%execute(cline)
            case( 'prob_tab2D' )
                call xprob_tab2D%execute(cline)
            case( 'prob_align2D' )
                call xprob_align2D%execute(cline)
            case( 'flex_eigenvol_reconstruct' )
                call xflex_eigenvol_reconstruct%execute(cline)
            ! RECONSTRUCTION PROGRAMS
            case( 'reconstruct3D' )
                call xrec3D%execute(cline)
            case( 'volassemble' )
                call xvolassemble%execute(cline)
            ! CHECKER PROGRAMS
            case( 'check_box' )
                call xcheck_box%execute(cline)
            case( 'check_nptcls' )
                call xcheck_nptcls%execute(cline)
            case( 'check_stoch_update' )
                call xcheck_stoch_update%execute(cline)
            case( 'check_update_frac' )
                call xcheck_update_frac%execute(cline)
            ! VOLOPS PROGRAMS
            case( 'postprocess' )
                call xpostprocess%execute(cline)
            case( 'automask' )
                call xautomask%execute(cline)
            ! GENERAL IMAGE PROCESSING PROGRAMS
            case( 'scale' )
                call xscale%execute(cline)
            case( 'binarize' )
                call xbinarize%execute(cline)
            case( 'ppca_denoise' )
                call xppca_denoise%execute(cline)
            ! MISCELLANOUS PROGRAMS
            case( 'aggregate_chunks' )
                call xaggregate_chunks%execute(cline)
            case( 'fractionate_movies' )
                call xfractionate_movies%execute(cline)
            case( 'kstest' )
                call xkstst%execute(cline)
            case( 'pearsn' )
                call xpearsn%execute(cline)
            ! ORIENTATION DATA MANAGEMENT PROGRAMS
            case( 'rotmats2oris' )
                call xrotmats2oris%execute(cline)
            case( 'print_project_vals' )
                call xprint_project_vals%execute(cline)
                l_silent = .true.
            ! DATA MANAGEMENT PROGRAMS
            case( 'prune_project' )
                call xprune_project%execute(cline)
            case( 'scale_project' )
                call xscale_project%execute(cline)
            ! TIME-SERIES ANALYSIS PROGRAMS
            case( 'tseries_motion_correct' )
                call xtseries_mcorr%execute(cline)
            case( 'track_particles' )
                call xtrack_particles%execute(cline)
            ! PARALLEL PROCESSING PROGRAMS
            case( 'split' )
                call xsplit%execute(cline)
            case DEFAULT
                THROW_HARD('prg='//trim(prg)//' is unsupported')
        end select
    end subroutine dispatch_private_prg

end module simple_private_exec_driver
