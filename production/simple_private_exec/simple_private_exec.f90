! shared-memory parallelised programs executed by distributed commanders
program simple_private_exec
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline, cmdline_err
use simple_user_interface, only: make_user_interface, write_ui_json, print_ui_latex
use simple_private_prgs,   only: make_private_user_interface
use simple_symanalyzer,    only: print_subgroups
use simple_commander_project
use simple_commander_checks
use simple_commander_distr
use simple_commander_misc
use simple_commander_mask
use simple_commander_imgproc
use simple_commander_oris
use simple_commander_preprocess
use simple_commander_cluster2D
use simple_commander_refine3D
use simple_commander_rec
use simple_commander_sim
use simple_commander_volops
use simple_commander_tseries
implicit none
#include "simple_local_flags.inc"
! include 'git_version.inc'

! PRE-PROCESSING PROGRAMS
type(preprocess_commander)            :: xpreprocess
type(extract_commander)               :: xextract
type(reextract_commander)             :: xreextract
type(motion_correct_commander)        :: xmotion_correct
type(gen_pspecs_and_thumbs_commander) :: xgen_pspecs_and_thumbs
type(ctf_estimate_commander)          :: xctf_estimate
type(pick_extract_commander)          :: xpick_extract
type(pick_commander)                  :: xpick
type(make_pickrefs_commander)         :: xmake_pickrefs

! CLUSTER2D PROGRAMS
type(make_cavgs_commander)            :: xmake_cavgs
type(cluster2D_commander)             :: xcluster2D
type(cluster2D_commander_distr)       :: xcluster2D_distr
type(cavgassemble_commander)          :: xcavgassemble
type(rank_cavgs_commander)            :: xrank_cavgs
type(export_cavgs_commander)          :: xexport_cavgs

! REFINE3D PROGRAMS
type(refine3D_commander)              :: xrefine3D
type(calc_pspec_commander)            :: xcalc_pspec
type(calc_pspec_assemble_commander)   :: xcalc_pspec_assemble
type(check_3Dconv_commander)          :: xcheck_3Dconv
type(calc_group_sigmas_commander)     :: xcalc_group_sigmas

! RECONSTRUCTION PROGRAMS
type(volassemble_commander)           :: xvolassemble
type(reconstruct3D_commander)         :: xreconstruct3D

! CHECKER PROGRAMS
type(check_box_commander)             :: xcheck_box
type(check_nptcls_commander)          :: xcheck_nptcls

! VOLOPS PROGRAMS
type(postprocess_commander)           :: xpostprocess
type(automask_commander)              :: xautomask

! GENERAL IMAGE PROCESSING PROGRAMS
type(scale_commander)                 :: xscale
type(binarize_commander)              :: xbinarize
type(edge_detect_commander)           :: xdetector

! MISCELLANOUS PROGRAMS
type(masscen_commander)               :: xmasscen
type(stk_corr_commander)              :: xstk_corr
type(kstest_commander)                :: xkstst

! ORIENTATION DATA MANAGEMENT PROGRAMS
type(rotmats2oris_commander)          :: xrotmats2oris
type(print_project_vals_commander)    :: xprint_project_vals

! ORIENTATION DATA MANAGEMENT PROGRAMS
type(prune_project_commander)         :: xprune_project
type(scale_project_commander_distr)   :: xscale_project_distr

! TIME-SERIES ANALYSIS PROGRAMS
type(tseries_track_particles_commander) :: xtseries_track_particles
type(tseries_motion_correct_commander)  :: xtseries_mcorr
type(pspec_int_rank_commander)          :: xpspec_int_rank

! PARALLEL PROCESSING PROGRAMS
type(split_commander)                   :: xsplit

! OTHER DECLARATIONS
character(len=STDLEN) :: xarg, prg
type(cmdline)         :: cline
integer               :: cmdstat, cmdlen, pos
integer(timer_int_kind)                     :: t0
real(timer_int_kind)                        :: rt_exec

! start timer
t0 = tic()  

! print git version
! call simple_print_git_version(GIT_HASH)
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

select case(prg)

    ! PRIVATE UTILITY PROGRAMS
    case( 'write_ui_json' )
        call write_ui_json
    case( 'print_ui_latex' )
        call print_ui_latex
    case( 'print_sym_subgroups' )
        call print_subgroups

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
    case( 'make_pickrefs' )
        call xmake_pickrefs%execute(cline)

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

    ! REFINE3D PROGRAMS
    case( 'refine3D' )
        call xrefine3D%execute(cline)
    case( 'calc_pspec' )
        call xcalc_pspec%execute(cline)
    case( 'calc_pspec_assemble' )
        call xcalc_pspec_assemble%execute(cline)
    case( 'check_3Dconv' )
        call xcheck_3Dconv%execute(cline)
    case( 'calc_group_sigmas' )
        call xcalc_group_sigmas%execute(cline)

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
    case( 'masscen' )
        call xmasscen%execute(cline)
    case( 'stk_corr' )
        call xstk_corr%execute(cline)
    case( 'kstest' )
        call xkstst%execute(cline)

    ! ORIENTATION DATA MANAGEMENT PROGRAMS
    case( 'rotmats2oris' )
        call xrotmats2oris%execute(cline)
    case( 'print_project_vals' )
        call xprint_project_vals%execute(cline)

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
    case( 'pspec_int_rank' )
        call xpspec_int_rank%execute(cline)

    ! PARALLEL PROCESSING PROGRAMS
    case( 'split' )
        call xsplit%execute(cline)
    case DEFAULT
        THROW_HARD('prg='//trim(prg)//' is unsupported')
    end select
!end timer and print
rt_exec = toc(t0)
! call simple_print_timer(rt_exec)
end program simple_private_exec
