! shared-memory parallelised programs executed by distributed commanders
program simple_private_exec
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline, cmdline_err
use simple_user_interface, only: make_user_interface, print_ui_json, write_ui_json, print_stream_ui_json
use simple_private_prgs,   only: make_private_user_interface
use simple_symanalyzer,    only: print_subgroups
use simple_commanders_project
use simple_commanders_checks
use simple_commanders_distr
use simple_commanders_misc
use simple_commanders_mask
use simple_commanders_imgproc
use simple_commanders_oris
use simple_commanders_preprocess
use simple_commanders_cluster2D
use simple_commanders_cavgs
use simple_commanders_refine3D
use simple_commanders_euclid
use simple_commanders_rec
use simple_commanders_sim
use simple_commanders_volops
use simple_commanders_tseries
use simple_commanders_resolest
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
type(commander_fractionate_movies)      :: xfractionate_movies

! CLUSTER2D PROGRAMS
type(commander_make_cavgs)              :: xmake_cavgs
type(commander_cluster2D)               :: xcluster2D
type(commander_cluster2D_distr)         :: xcluster2D_distr
type(commander_cavgassemble)            :: xcavgassemble
type(commander_rank_cavgs)              :: xrank_cavgs
type(commander_export_cavgs)            :: xexport_cavgs
type(commander_prob_tab2D)              :: xprob_tab2D

! REFINE3D PROGRAMS
type(commander_refine3D)                :: xrefine3D
type(commander_calc_pspec_distr)        :: xcalc_pspec_distr
type(commander_calc_pspec)              :: xcalc_pspec
type(commander_calc_pspec_assemble)     :: xcalc_pspec_assemble
type(commander_check_3Dconv)            :: xcheck_3Dconv
type(commander_calc_group_sigmas)       :: xcalc_group_sigmas
type(commander_prob_tab)                :: xprob_tab
type(commander_prob_align)              :: xprob_align

! RECONSTRUCTION PROGRAMS
type(commander_volassemble)             :: xvolassemble
type(commander_reconstruct3D)           :: xreconstruct3D

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
type(commander_edge_detect)             :: xdetector

! MISCELLANOUS PROGRAMS
type(commander_kstest)                  :: xkstst
type(commander_pearsn)                  :: xpearsn

! ORIENTATION DATA MANAGEMENT PROGRAMS
type(commander_rotmats2oris)            :: xrotmats2oris
type(commander_print_project_vals)      :: xprint_project_vals

! ORIENTATION DATA MANAGEMENT PROGRAMS
type(commander_prune_project)           :: xprune_project
type(commander_scale_project_distr)     :: xscale_project_distr

! TIME-SERIES ANALYSIS PROGRAMS
type(commander_tseries_track_particles) :: xtseries_track_particles
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
pos = index(xarg, '=') ! position of '='
call cmdline_err( cmdstat, cmdlen, xarg, pos )
prg = xarg(pos+1:)     ! this is the program name
! make UIs
call make_user_interface
call make_private_user_interface
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
    case( 'fractionate_movies' )
        call xfractionate_movies%execute(cline)

    ! CLUSTER2D PROGRAMS
    case( 'make_cavgs' )
        call xmake_cavgs%execute(cline)
    case( 'cluster2D' )
        call xcluster2D%execute(cline)
    case( 'cluster2D_distr' )
        call xcluster2D_distr%execute(cline)
    case( 'cavgassemble' )
        call xcavgassemble%execute(cline)
    case( 'rank_cavgs' )
        call xrank_cavgs%execute(cline)
    case( 'export_cavgs' )
        call xexport_cavgs%execute(cline)
    case( 'prob_tab2D' )
        call xprob_tab2D%execute(cline)

    ! REFINE3D PROGRAMS
    case( 'refine3D' )
        call xrefine3D%execute(cline)
    case( 'calc_pspec_distr' )
        call xcalc_pspec_distr%execute(cline)
    case( 'calc_pspec' )
        call xcalc_pspec%execute(cline)
    case( 'calc_pspec_assemble' )
        call xcalc_pspec_assemble%execute(cline)
    case( 'check_3Dconv' )
        call xcheck_3Dconv%execute(cline)
    case( 'calc_group_sigmas' )
        call xcalc_group_sigmas%execute(cline)
    case( 'prob_tab' )
        call xprob_tab%execute(cline)
    case( 'prob_align' )
        call xprob_align%execute(cline)

    ! RECONSTRUCTION PROGRAMS
    case( 'reconstruct3D' )
        call xreconstruct3D%execute(cline)
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
    case('edge_detect')
        call xdetector%execute(cline)

    ! MISCELLANOUS PROGRAMS
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
    case( 'scale_project_distr' )
        ! for convenience
        call xscale_project_distr%execute(cline)

    ! TIME-SERIES ANALYSIS PROGRAMS
    case( 'tseries_motion_correct' )
        call xtseries_mcorr%execute(cline)
    case( 'tseries_track_particles' )
        call xtseries_track_particles%execute(cline)

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
end program simple_private_exec
