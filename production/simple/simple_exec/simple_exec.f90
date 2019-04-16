! executes the shared-memory parallelised programs in SIMPLE
program simple_exec
include 'simple_lib.f08'
use simple_user_interface, only: make_user_interface,list_shmem_prgs_in_ui
use simple_cmdline,        only: cmdline, cmdline_err
use simple_spproj_hlev
use simple_commander_project
use simple_commander_checks
use simple_commander_distr
use simple_commander_imgproc
use simple_commander_mask
use simple_commander_misc
use simple_commander_oris
use simple_commander_preprocess
use simple_commander_cluster2D
use simple_commander_refine3D
use simple_commander_rec
use simple_commander_sim
use simple_commander_star
use simple_commander_volops
use simple_commander_tseries
use simple_commander_resolest
use simple_projection_frcs
implicit none
#include "simple_local_flags.inc"

! PROJECT MANAGEMENT PROGRAMS
type(new_project_commander)         :: xnew_project
type(update_project_commander)      :: xupdate_project
type(print_project_info_commander)  :: xprint_project_info
type(print_project_field_commander) :: xprint_project_field
type(import_movies_commander)       :: ximport_movies
type(import_boxes_commander)        :: ximport_boxes
type(import_particles_commander)    :: ximport_particles
type(import_cavgs_commander)        :: ximport_cavgs
type(prune_project_commander)       :: xprune_project
type(export_starproject_commander)  :: xexport_starproject
type(import_starproject_commander)  :: ximport_starproject
type(report_selection_commander)    :: xreport_selection

! SINGLE-PARTICLE WORKFLOW PROGRAMS
type(cluster_cavgs_commander)  :: xcluster_cavgs
type(symaxis_search_commander) :: xsymsrch
type(symmetry_test_commander)  :: xsymtst
type(postprocess_commander)    :: xpostprocess

! IMAGE PROCESSING PROGRAMS
type(pspec_stats_commander) :: xpspecstats
type(mask_commander)        :: xmask
type(fsc_commander)         :: xfsc
type(local_res_commander)   :: xlocal_res
type(centervol_commander)   :: xcenter
type(reproject_commander)   :: xreproject
type(volops_commander)      :: xvolops
type(convert_commander)     :: xconvert
type(ctfops_commander)      :: xctfops
type(filter_commander)      :: xfilter
type(normalize_commander)   :: xnormalize
type(scale_commander)       :: xscale
type(stack_commander)       :: xstack
type(stackops_commander)    :: xstackops
type(shift_commander)       :: xshift

! NANOPARTICLES PROGRAMS
type(detect_atoms_commander)  :: xdetectatoms  !I put it in simple_commander_preprocess not to write another file
type(compare_nano_commander)  :: xcomparenano


! ORIENTATION PROCESSING PROGRAMS
type(make_oris_commander) :: xmake_oris
type(orisops_commander)   :: xorisops
type(oristats_commander)  :: xoristats
type(vizoris_commander)   :: xvizoris

! PRINT INFO PROGRAMS
type(info_image_commander)        :: xinfo_image
type(info_stktab_commander)       :: xinfo_stktab
type(print_fsc_commander)         :: xprint_fsc
type(print_magic_boxes_commander) :: xprint_magic_boxes

! SIMULATOR PROGRAMS
type(simulate_noise_commander)       :: xsimulate_noise
type(simulate_particles_commander)   :: xsimulate_particles
type(simulate_movie_commander)       :: xsimulate_movie
type(simulate_subtomogram_commander) :: xsimulate_subtomogram

! TIME-SERIES PROGRAMS
type(tseries_import_commander)       :: xtseries_import
type(tseries_ctf_estimate_commander) :: xtseries_ctf_estimate

! SYSTEM INTERACTION PROGRAMS
type(mkdir_commander) :: xmkdir

! OTHER DECLARATIONS
character(len=STDLEN) :: xarg, prg, entire_line
type(cmdline)         :: cline
integer               :: cmdstat, cmdlen, pos

! parse command-line
call get_command_argument(1, xarg, cmdlen, cmdstat)
call get_command(entire_line)
pos = index(xarg, '=') ! position of '='
call cmdline_err( cmdstat, cmdlen, xarg, pos )
prg = xarg(pos+1:)     ! this is the program name
! make UI
call make_user_interface
if( str_has_substr(entire_line, 'prg=list') )then
    call list_shmem_prgs_in_ui
    stop
endif
! parse command line into cline object
call cline%parse

select case(prg)

    ! PROJECT MANAGEMENT PROGRAMS

    case( 'new_project' )
        call xnew_project%execute(cline)
    case( 'update_project' )
        call xupdate_project%execute(cline)
    case( 'print_project_info' )
        call xprint_project_info%execute(cline)
    case( 'print_project_field' )
        call xprint_project_field%execute(cline)
    case( 'import_movies' )
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('ctf')   ) call cline%set('ctf',   'yes')
        call ximport_movies%execute(cline)
    case( 'import_boxes' )
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call ximport_boxes%execute(cline)
    case( 'import_particles' )
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('ctf')   ) call cline%set('ctf',   'yes')
        call ximport_particles%execute(cline)
    case( 'import_cavgs' )
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call ximport_cavgs%execute(cline)
    case( 'prune_project' )
        call cline%set('mkdir', 'yes')
        call xprune_project%execute(cline)
    case( 'export_starproject' )
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',     'yes')
        if( .not. cline%defined('starfile')) call cline%set('starfile', 'NONE')
        call xexport_starproject%execute(cline)
    case( 'import_starproject' )
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',     'yes')
        if( .not. cline%defined('starfile')) call cline%set('starfile', 'NONE')
        call ximport_starproject%execute(cline)
    case( 'report_selection' )
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call xreport_selection%execute(cline)

    ! SINGLE-PARTICLE WORKFLOW PROGRAMS

    case('cluster_cavgs')
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir', 'yes')
        call xcluster_cavgs%execute(cline)
    case( 'symaxis_search' )
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir',  'yes')
        if( .not. cline%defined('cenlp')  ) call cline%set('cenlp',    30.)
        if( .not. cline%defined('center') ) call cline%set('center', 'yes')
        call xsymsrch%execute( cline )
    case( 'symmetry_test' )
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir',  'yes')
        if( .not. cline%defined('cenlp')  ) call cline%set('cenlp',    20.)
        if( .not. cline%defined('center') ) call cline%set('center', 'yes')
        call xsymtst%execute( cline )
    case( 'postprocess' )
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir', 'yes')
        call xpostprocess%execute(cline)


    ! IMAGE PROCESSING PROGRAMS
    case('pspec_stats')
        call xpspecstats%execute(cline)
    case( 'mask' )
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call xmask%execute(cline)
    case( 'fsc' )
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call xfsc%execute(cline)
    case( 'local_resolution' )
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call xlocal_res%execute(cline)
    case( 'center' )
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('cenlp') ) call cline%set('cenlp',   30.)
        call xcenter%execute(cline)
    case( 'reproject' )
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('wfun')  ) call cline%set('wfun',   'kb')
        if( .not. cline%defined('winsz') ) call cline%set('winsz',   1.5)
        if( .not. cline%defined('alpha') ) call cline%set('alpha',    2.)
        call xreproject%execute(cline)
    case( 'volops' )
        call xvolops%execute(cline)
    case( 'convert' )
        call xconvert%execute(cline)
    case( 'ctfops' )
        if( .not. cline%defined('stk')   ) call cline%set('box',    256.)
        call xctfops%execute(cline)
    case( 'filter' )
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call xfilter%execute(cline)
    case( 'normalize' )
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call xnormalize%execute(cline)
    case( 'scale' )
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call xscale%execute(cline)
    case( 'stack' )
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call xstack%execute(cline)
    case( 'stackops' )
        if( .not. cline%defined('outfile') ) call cline%set('outfile', 'outfile.txt')
        call xstackops%execute(cline)
    case( 'shift' )
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call xshift%execute(cline)

    ! NANOPARTICLES PROGRAMS
    case('detect_atoms')
        call xdetectatoms%execute(cline)
    case('compare_nano')
        call xcomparenano%execute(cline)

    ! ORIENTATION PROCESSING PROGRAMS

    case( 'make_oris' )
        call xmake_oris%execute(cline)
    case( 'orisops' )
        call xorisops%execute(cline)
    case( 'oristats' )
        call xoristats%execute(cline)
    case( 'vizoris' )
        call xvizoris%execute(cline)

    ! PRINT INFO PROGRAMS

    case( 'info_image' )
        call xinfo_image%execute(cline)
    case( 'info_stktab' )
        call xinfo_stktab%execute(cline)
    case( 'print_fsc' )
        call xprint_fsc%execute(cline)
    case( 'print_magic_boxes' )
        call xprint_magic_boxes%execute(cline)

    ! SIMULATOR PROGRAMS

    case( 'simulate_noise' )
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call xsimulate_noise%execute(cline)
    case( 'simulate_particles' )
        if( .not. cline%defined('mkdir')    ) call cline%set('mkdir', 'yes')
        call cline%set('nspace', cline%get_rarg('nptcls'))
        if( .not. cline%defined('sherr') .and. .not. cline%defined('oritab') ) call cline%set('sherr', 2.)
        if( .not. cline%defined('ctf')      ) call cline%set('ctf',    'yes')
        if( .not. cline%defined('dferr')    ) call cline%set('dferr',    1.5)
        if( .not. cline%defined('astigerr') ) call cline%set('astigerr', 0.5)
        if( .not. cline%defined('bfacerr')  ) call cline%set('bfacerr',  0.0)
        call cline%set('wfun', 'kb')
        call cline%set('winsz', 1.5)
        call cline%set('alpha',  2.)
        call cline%set('eo', '  no')
        call xsimulate_particles%execute(cline)
    case( 'simulate_movie' )
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('trs')     ) call cline%set('trs',      3.)
        if( .not. cline%defined('ctf')     ) call cline%set('ctf',   'yes')
        if( .not. cline%defined('bfac')    ) call cline%set('bfac',   200.)
        if( .not. cline%defined('nframes') ) call cline%set('nframes', 30.)
        call cline%set('wfun', 'kb')
        call cline%set('winsz', 1.5)
        call cline%set('alpha',  2.)
        call cline%set('eo',   'no')
        call xsimulate_movie%execute(cline)
    case( 'simulate_subtomogram' )
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call xsimulate_subtomogram%execute(cline)

    ! TIME-SERIES PROGRAM

    case( 'tseries_import' )
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( cline%defined('ctf') )then
            select case(trim(cline%get_carg('ctf')))
            case('no','flip')
                ! all good
            case DEFAULT
                THROW_HARD('Only CTF=NO|FLIP are allowed!')
            end select
        endif
        call xtseries_import%execute(cline)
    case( 'tseries_ctf_estimate' )
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        ! set defaults
        if( .not. cline%defined('hp') ) call cline%set('hp', 10.)
        if( .not. cline%defined('lp') ) call cline%set('lp', 2.5)
        call xtseries_ctf_estimate%execute(cline)

    ! SYSTEM INTERACTION PROGRAMS

    case( 'mkdir' )
        call cline%set('mkdir', 'yes')
        call xmkdir%execute(cline)
    case DEFAULT
        THROW_HARD('prg='//trim(prg)//' is unsupported')
end select
call update_job_descriptions_in_project( cline )
end program simple_exec
