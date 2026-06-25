!@descr: executes shared-memory parallelized programs executed by distributed commanders
program simple_private_exec
use simple_private_exec_api
implicit none
#include "simple_local_flags.inc"

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
type(commander_diffmap_denoise_project) :: xdiffmap_denoise_project

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

! RECONSTRUCTION PROGRAMS
type(commander_volassemble)   :: xvolassemble
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

! ORIENTATION DATA MANAGEMENT PROGRAMS
type(commander_prune_project)           :: xprune_project
type(commander_scale_project)     :: xscale_project

! TIME-SERIES ANALYSIS PROGRAMS
type(commander_track_particles)         :: xtrack_particles
type(commander_tseries_motion_correct)  :: xtseries_mcorr

! PARALLEL PROCESSING PROGRAMS
type(commander_split)                   :: xsplit

! OTHER DECLARATIONS
character(len=STDLEN)   :: xarg, prg
type(cmdline)           :: cline
integer                 :: cmdstat, cmdlen, pos
integer(timer_int_kind) :: t0
real(timer_int_kind)    :: rt_exec
logical                 :: l_silent

! start timer
t0 = tic()
! parse command-line
call get_command_argument(1, xarg, cmdlen, cmdstat)
if( cmdstat /= 0 ) call cmdline_err( cmdstat, cmdlen, xarg, 0 )
if( trim(xarg) == '--coarray' )then
    call run_coarray_direct
    stop
endif
pos = index(xarg, '=') ! position of '='
call cmdline_err( cmdstat, cmdlen, xarg, pos )
prg = xarg(pos+1:)     ! this is the program name
! make UIs
call make_ui
call make_private_ui
! this parses all key=value pairs on the command line
call cline%parse_private
call print_slurm_env
l_silent = .false.
select case(prg)

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
    case( 'diffmap_denoise_project' )
        call xdiffmap_denoise_project%execute(cline)

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
        ! for convenience
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
! end timer and print
rt_exec = toc(t0)
if( .not. l_silent ) call simple_print_timer(rt_exec)
! cleanup
call cline%kill

contains

    !> Coarray backend entry point.  The parent invocation is:
    !!   cafrun -np <nimages> simple_private_exec --coarray <first> <last> <numlen> <job_args...>
    !! Each image runs the same simple_private_exec payload used by script-based
    !! backends for its assigned partitions, but no distr_simple_script_* files
    !! are written.
    subroutine run_coarray_direct
#ifdef USE_COARRAYS
        character(len=*), parameter :: UNSET_MPI_ENV_LOCAL = 'env -u OMPI_COMM_WORLD_SIZE -u OMPI_COMM_WORLD_RANK '//&
            &'-u OMPI_COMM_WORLD_LOCAL_SIZE -u OMPI_COMM_WORLD_LOCAL_RANK -u OMPI_COMM_WORLD_NODE_RANK '//&
            &'-u PMIX_NAMESPACE -u PMIX_RANK -u PMIX_SERVER_URI2 -u PMIX_SERVER_URI3 '//&
            &'-u PMI_RANK -u PMI_SIZE -u PMI_FD -u PMI_PORT -u PMI_ID -u PMI_KVSNAME -u PMI_PROCESS_MAPPING '//&
            &'-u PMI_APPNUM -u PMI_SPAWNED -u PMI_CONTROL_PORT -u PMI_CONTROL_FD '//&
            &'-u HYDRA_PROXY_ID -u HYDRA_PROXY_PORT -u HYDRA_CONTROL_FD -u HYDI_CONTROL_FD '//&
            &'-u MPI_LOCALRANKID -u MPI_LOCALNRANKS '//&
            &'-u MV2_COMM_WORLD_RANK -u MV2_COMM_WORLD_SIZE -u MV2_COMM_WORLD_LOCAL_RANK -u MV2_COMM_WORLD_LOCAL_SIZE'
        character(len=XLONGSTRLEN)              :: exec_binary, cmdmsg, job_log
        character(len=XLONGSTRLEN), allocatable :: job_args(:)
        character(len=:), allocatable           :: cmd, errmsg
        integer :: nargs, njobs, first_part, last_part, numlen_local
        integer :: iarg, ipart, job_ind, cmdstat_local, exitstat, syncstat, argstat, arglen
        nargs = command_argument_count()
        if( nargs < 5 )then
            THROW_HARD('Usage: simple_private_exec --coarray <first_part> <last_part> <numlen> <job_args...>')
        endif
        call parse_coarray_int_arg(2, first_part, 'first_part')
        call parse_coarray_int_arg(3, last_part,  'last_part')
        call parse_coarray_int_arg(4, numlen_local, 'numlen')
        if( first_part < 1 .or. last_part < first_part .or. numlen_local < 1 )then
            THROW_HARD('Invalid coarray partition range or numlen')
        endif
        njobs = last_part - first_part + 1
        if( nargs /= 4 + njobs )then
            THROW_HARD('simple_private_exec --coarray received wrong number of job arguments')
        endif
        allocate(job_args(njobs))
        do iarg = 1, njobs
            call get_command_argument(4 + iarg, job_args(iarg), arglen, argstat)
            if( argstat /= 0 ) THROW_HARD('could not read coarray job argument for part '//int2str(first_part + iarg - 1))
            if( len_trim(job_args(iarg)) == 0 ) THROW_HARD('empty coarray job argument for part '//int2str(first_part + iarg - 1))
        end do
        call get_command_argument(0, exec_binary)
        if( len_trim(exec_binary) == 0 ) exec_binary = 'simple_private_exec'
        do ipart = first_part + this_image() - 1, last_part, num_images()
            job_ind  = ipart - first_part + 1
            job_log  = 'simple_private_exec_coarray_part_'//int2str_pad(ipart, numlen_local)//'.out'
            cmd      = UNSET_MPI_ENV_LOCAL//' SIMPLE_COARRAY_IMAGE='//int2str(this_image())//' SIMPLE_COARRAY_IMAGES='//&
                &int2str(num_images())//' '//shell_quote(trim(exec_binary))//' '//trim(job_args(job_ind))//' >> '//&
                &shell_quote(trim(job_log))//' 2>&1'
            cmdstat_local = 0
            exitstat      = 0
            cmdmsg        = ''
            call execute_command_line(trim(cmd), wait=.true., exitstat=exitstat, cmdstat=cmdstat_local, cmdmsg=cmdmsg)
            if( cmdstat_local /= 0 .or. exitstat /= 0 )then
                if( len_trim(cmdmsg) > 0 ) THROW_WARN(trim(cmdmsg))
                errmsg = 'simple_private_exec coarray image '//int2str(this_image())//' failed part '//int2str(ipart)
                errmsg = errmsg//' exitstat='//int2str(exitstat)//' cmdstat='//int2str(cmdstat_local)
                errmsg = errmsg//' log='//trim(job_log)
                THROW_HARD(errmsg)
            endif
        end do
        sync all(stat=syncstat)
        if( syncstat /= 0 ) THROW_HARD('simple_private_exec coarray sync failed with stat='//int2str(syncstat))
#else
        THROW_HARD('simple_private_exec --coarray requires a USE_COARRAYS build')
#endif
    end subroutine run_coarray_direct

#ifdef USE_COARRAYS
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

    function shell_quote( raw ) result( quoted )
        character(len=*), intent(in) :: raw
        character(len=:), allocatable :: quoted
        integer :: i
        quoted = achar(39)
        do i = 1, len_trim(raw)
            if( raw(i:i) == achar(39) )then
                quoted = quoted//achar(39)//achar(92)//achar(39)//achar(39)
            else
                quoted = quoted//raw(i:i)
            endif
        end do
        quoted = quoted//achar(39)
    end function shell_quote
#endif

end program simple_private_exec
