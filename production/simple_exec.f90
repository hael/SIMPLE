!@descr: executes SIMPLE workflows
program simple_exec
use simple_exec_api
use simple_exec_project,    only: exec_project_commander
use simple_exec_preproc,    only: exec_preproc_commander
use simple_exec_cluster2D,  only: exec_cluster2D_commander
use simple_exec_cavgproc,   only: exec_cavgproc_commander
use simple_exec_abinitio3D, only: exec_abinitio3D_commander
use simple_exec_refine3D,   only: exec_refine3D_commander
implicit none
#include "simple_local_flags.inc"

! DENOISING
type(commander_icm2D)                       :: xicm2D
type(commander_icm3D)                       :: xicm3D
type(commander_ppca_denoise)                :: xppca_denoise
type(commander_ppca_denoise_classes)        :: xppca_denoise_classes

! FILTERING
type(commander_filter)                      :: xfilter
type(commander_uniform_filter2D)            :: xuniform_filter2D
type(commander_uniform_filter3D)            :: xuniform_filter3D

! GENERAL IMAGE PROCESSING
type(commander_binarize)                    :: xbinarize
type(commander_convert)                     :: xconvert
type(commander_ctf_phaseflip)               :: xctf_phaseflip
type(commander_ctfops)                      :: xctfops
type(commander_normalize)                   :: xnormalize
type(commander_scale)                       :: xscale
type(commander_stack)                       :: xstack
type(commander_stackops)                    :: xstackops

! MASKING
type(commander_auto_spher_mask)             :: xauto_spher_mask
type(commander_automask2D)                  :: xautomask2D
type(commander_mask)                        :: xmask

! ORIENTATION PROCESSING
type(commander_make_oris)                   :: xmake_oris
type(commander_orisops)                     :: xorisops
type(commander_oristats)                    :: xoristats
type(commander_vizoris)                     :: xvizoris

! PARALLEL UTILITIES
type(commander_split)                       :: xsplit

! PRINT INFO
type(commander_info_image)                  :: xinfo_image
type(commander_info_stktab)                 :: xinfo_stktab
type(commander_print_dose_weights)          :: xprint_dose_weights
type(commander_print_fsc)                   :: xprint_fsc
type(commander_print_magic_boxes)           :: xprint_magic_boxes

! RESOLUTION ESTIMATON
type(commander_clin_fsc)                    :: xclin_fsc
type(commander_fsc)                         :: xfsc

! SIMULATION
type(commander_pdb2mrc)                     :: xpdb2mrc
type(commander_simulate_movie)              :: xsimulate_movie
type(commander_simulate_noise)              :: xsimulate_noise
type(commander_simulate_particles)          :: xsimulate_particles

! STREAM VALIDATION
type(commander_check_refpick)               :: xcheck_refpick
type(commander_mini_stream)                 :: xmini_stream

! SYMMETRY
type(commander_symaxis_search)              :: xsymsrch
type(commander_symmetrize_map)              :: xsymmetrize_map
type(commander_symmetry_test)               :: xsymtst

! VOLUME DOCKING
type(commander_dock_volpair)                :: xdock_volpair
type(commander_volanalyze)                  :: xvolanalyze

! VOLUME PROCESSING
type(commander_centervol)                   :: xcenter
type(commander_reproject)                   :: xreproject
type(commander_volops)                      :: xvolops

! MODEL ANALYSIS
type(commander_model_validation)            :: xmodel_validation

! OTHER DECLARATIONS
character(len=STDLEN)      :: xarg, prg
character(len=XLONGSTRLEN) :: entire_line
type(cmdline)              :: cline
integer                    :: cmdstat, cmdlen, pos
integer(timer_int_kind)    :: t0
real(timer_int_kind)       :: rt_exec
logical                    :: l_silent, l_did_execute

! start timer
t0 = tic()
! parse command-line
call get_command_argument(1, xarg, cmdlen, cmdstat)
call get_command(entire_line)
pos = index(xarg, '=') ! position of '='
call cmdline_err( cmdstat, cmdlen, xarg, pos )
prg = xarg(pos+1:)     ! this is the program name
! make UI
call make_user_interface
if( str_has_substr(entire_line, 'prg=list') )then
    call list_simple_prgs_in_ui
    stop
endif
! parse command line into cline object
call cline%parse
! generate script for queue submission?
call script_exec(cline, string(trim(prg)), string('simple_exec'))
l_silent      = .false.
l_did_execute = .false.

call exec_project_commander(trim(prg),    cline, l_silent, l_did_execute)
call exec_preproc_commander(trim(prg),    cline, l_silent, l_did_execute)
call exec_cluster2D_commander(trim(prg),  cline, l_silent, l_did_execute)
call exec_cavgproc_commander(trim(prg),   cline, l_silent, l_did_execute)
call exec_abinitio3D_commander(trim(prg), cline, l_silent, l_did_execute)
call exec_refine3D_commander(trim(prg),   cline, l_silent, l_did_execute)

select case(trim(prg))


    !====================================================================
    ! DENOISING
    !====================================================================
    case( 'icm2D' )
        call xicm2D%execute(cline)
    case( 'icm3D' )
        call xicm3D%execute(cline)
    case( 'ppca_denoise' )
        call xppca_denoise%execute(cline)
    case( 'ppca_denoise_classes' )
        call xppca_denoise_classes%execute(cline)

    !====================================================================
    ! FILTERING
    !====================================================================
    case( 'filter' )
        call xfilter%execute(cline)
    case( 'uniform_filter2D' )
        call xuniform_filter2D%execute(cline)
    case( 'uniform_filter3D' )
        call xuniform_filter3D%execute(cline)

    !====================================================================
    ! GENERAL IMAGE PROCESSING
    !====================================================================
    case( 'binarize' )
        call xbinarize%execute(cline)
    case( 'convert' )
        call xconvert%execute(cline)
    case( 'ctf_phaseflip' )
        call xctf_phaseflip%execute(cline)
    case( 'ctfops' )
        call xctfops%execute(cline)
    case( 'normalize' )
        call xnormalize%execute(cline)
    case( 'scale' )
        call xscale%execute(cline)
    case( 'stack' )
        call xstack%execute(cline)
    case( 'stackops' )
        call xstackops%execute(cline)

    !====================================================================
    ! MASKING
    !====================================================================
    case( 'auto_spher_mask' )
        call xauto_spher_mask%execute(cline)
    case( 'automask2D' )
        call xautomask2D%execute(cline)
    case( 'mask' )
        call xmask%execute(cline)

    !====================================================================
    ! ORIENTATION PROCESSING
    !====================================================================
    case( 'make_oris' )
        call xmake_oris%execute(cline)
    case( 'orisops' )
        call xorisops%execute(cline)
    case( 'oristats' )
        call xoristats%execute(cline)
    case( 'vizoris' )
        call xvizoris%execute(cline)

    !====================================================================
    ! PARALLEL UTILITIES
    !====================================================================
    case( 'split' )
        call xsplit%execute(cline)

    !====================================================================
    ! PRINT INFO
    !====================================================================
    case( 'info_image' )
        call xinfo_image%execute(cline)
    case( 'info_stktab' )
        call xinfo_stktab%execute(cline)
    case( 'print_dose_weights' )
        call xprint_dose_weights%execute(cline)
        l_silent = .true.
    case( 'print_fsc' )
        call xprint_fsc%execute(cline)
        l_silent = .true.
    case( 'print_magic_boxes' )
        call xprint_magic_boxes%execute(cline)
        l_silent = .true.

    !====================================================================
    ! RESOLUTION ESTIMATION
    !====================================================================
    case( 'clin_fsc' )
        call xclin_fsc%execute(cline)
    case( 'fsc' )
        call xfsc%execute(cline)

    !====================================================================
    ! SIMULATION
    !====================================================================
    case( 'pdb2mrc' )
        call xpdb2mrc%execute(cline)
    case( 'simulate_movie' )
        call xsimulate_movie%execute(cline)
    case( 'simulate_noise' )
        call xsimulate_noise%execute(cline)
    case( 'simulate_particles' )
        call xsimulate_particles%execute(cline)

    !====================================================================
    ! STREAM VALIDATION
    !====================================================================
    case( 'check_refpick' )
        call xcheck_refpick%execute(cline)
    case( 'mini_stream' )
        call xmini_stream%execute(cline)

    !====================================================================
    ! SYMMETRY
    !====================================================================
    case( 'symaxis_search' )
        call xsymsrch%execute(cline)
    case( 'symmetrize_map' )
        call xsymmetrize_map%execute(cline)
    case( 'symmetry_test' )
        call xsymtst%execute(cline)

    !====================================================================
    ! VOLUME DOCKING
    !====================================================================
    case( 'dock_volpair' )
        call xdock_volpair%execute(cline)
    case( 'volanalyze' )
        call xvolanalyze%execute(cline)

    !====================================================================
    ! VOLUME PROCESSING
    !====================================================================
    case( 'center' )
        call xcenter%execute(cline)
    case( 'reproject' )
        call xreproject%execute(cline)
    case( 'volops' )
        call xvolops%execute(cline)

    !====================================================================
    ! MODEL ANALYSIS
    !====================================================================
    case( 'model_validation' )
        call xmodel_validation%execute(cline)

    case default
        THROW_HARD('prg='//trim(prg)//' is unsupported')

end select

call update_job_descriptions_in_project( cline )
! close log file
if( logfhandle .ne. OUTPUT_UNIT )then
    if( is_open(logfhandle) ) call fclose(logfhandle)
endif
if( .not. l_silent )then
    call simple_print_git_version('9961f231')
    ! end timer and print
    rt_exec = toc(t0)
    call simple_print_timer(rt_exec)
endif
end program simple_exec
