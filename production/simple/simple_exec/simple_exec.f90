! executes the shared-memory parallelised programs in SIMPLE
program simple_exec
include 'simple_lib.f08'
use simple_user_interface
use simple_cmdline, only: cmdline, cmdline_err
use simple_commander_project
use simple_commander_checks
use simple_commander_comlin
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
use simple_commander_volops
use simple_commander_tseries
implicit none

! PROJECT MANAGEMENT
type(new_project_commander)          :: xnew_project
type(update_project_commander)       :: xupdate_project
type(print_project_info_commander)   :: xprint_project_info
type(import_movies_commander)        :: ximport_movies
type(import_particles_commander)     :: ximport_particles
type(import_cavgs_commander)         :: ximport_cavgs

! PART OF SP WORKFLOW
type(make_pickrefs_commander)        :: xmake_pickrefs
type(extract_commander)              :: xextract
type(cluster_cavgs_commander)        :: xcluster_cavgs
type(postprocess_commander)          :: xpostprocess

! IMAGE PROCESSING
type(mask_commander)                 :: xmask
type(fsc_commander)                  :: xfsc
type(center_commander)               :: xcenter
type(project_commander)              :: xproject
type(volops_commander)               :: xvolops
type(convert_commander)              :: xconvert
type(ctfops_commander)               :: xctfops
type(filter_commander)               :: xfilter
type(normalize_commander)            :: xnormalize
type(scale_commander)                :: xscale
type(stack_commander)                :: xstack
type(stackops_commander)             :: xstackops
type(shift_commander)                :: xshift

! ORIENTATION PROCESSING
type(make_oris_commander)            :: xmake_oris
type(orisops_commander)              :: xorisops
type(oristats_commander)             :: xoristats
type(vizoris_commander)              :: xvizoris

! PRINT INFO
type(info_image_commander)           :: xinfo_image
type(info_stktab_commander)          :: xinfo_stktab
type(print_fsc_commander)            :: xprint_fsc
type(print_magic_boxes_commander)    :: xprint_magic_boxes

! SIMULATORS
type(simulate_noise_commander)       :: xsimulate_noise
type(simulate_particles_commander)   :: xsimulate_particles
type(simulate_movie_commander)       :: xsimulate_movie
type(simulate_subtomogram_commander) :: xsimulate_subtomogram

! OTHER DECLARATIONS
!character(len=KEYLEN) :: keys_required(MAXNKEYS)='', keys_optional(MAXNKEYS)=''
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
if( str_has_substr(entire_line, 'prg=list') ) call list_shmem_prgs_in_ui

select case(prg)

    ! PROJECT MANAGEMENT

    case( 'new_project' )
        call cline%parse()
        call xnew_project%execute(cline)
    case( 'update_project' )
        call cline%parse()
        call xupdate_project%execute(cline)
    case( 'print_project_info' )
        call cline%parse()
        call xprint_project_info%execute(cline)
    case( 'import_movies' )
        call cline%parse()
        call ximport_movies%execute(cline)
    case( 'import_particles' )
        call cline%parse()
        call ximport_particles%execute(cline)
    case( 'import_cavgs' )
        call cline%parse()
        call ximport_cavgs%execute(cline)

    ! PART OF SP WORKFLOW

    case( 'make_pickrefs' )
        call cline%parse()
        if( .not. cline%defined('pcontrast')) call cline%set('pcontrast', 'black')
        if( .not. cline%defined('pgrp')     ) call cline%set('pgrp',      'd1'   )
        call xmake_pickrefs%execute(cline)
    case( 'extract' )
        call cline%parse()
        if( .not. cline%defined('pcontrast') )call cline%set('pcontrast', 'black')
        call cline%set('mkdir','yes')
        call xextract%execute(cline)
    case('cluster_cavgs')
        call cline%parse()
        call xcluster_cavgs%execute(cline)
    case( 'postprocess' )
        call cline%parse()
        call xpostprocess%execute(cline)

    ! IMAGE PROCESSING

    case( 'mask' )
        call cline%parse()
        call xmask%execute(cline)
    case( 'fsc' )
        call cline%parse()
        call xfsc%execute(cline)
    case( 'center' )
        call cline%parse()
        if( .not. cline%defined('cenlp') ) call cline%set('cenlp', 30.)
        call xcenter%execute(cline)
    case( 'project' )
        call cline%parse()
        if( .not. cline%defined('wfun')  ) call cline%set('wfun', 'kb')
        if( .not. cline%defined('winsz') ) call cline%set('winsz', 1.5)
        if( .not. cline%defined('alpha') ) call cline%set('alpha', 2.)
        call xproject%execute(cline)
    case( 'volops' )
        call cline%parse()
        call xvolops%execute(cline)
    case( 'convert' )
        call cline%parse()
        call xconvert%execute(cline)
    case( 'ctfops' )
        call cline%parse()
        if( .not. cline%defined('stk') ) call cline%set('box', 256.)
        call xctfops%execute(cline)
    case( 'filter' )
        call cline%parse()
        call xfilter%execute(cline)
    case( 'normalize' )
        call cline%parse()
        call xnormalize%execute(cline)
    case( 'scale' )
        call cline%parse()
        call xscale%execute(cline)
    case( 'stack' )
        call cline%parse()
        call xstack%execute(cline)
    case( 'stackops' )
        call cline%parse()
        if( .not. cline%defined('outfile') ) call cline%set('outfile', 'outfile.txt')
        call xstackops%execute(cline)
    case( 'shift' )
        call cline%parse()
        call xshift%execute(cline)

    ! ORIENTATION PROCESSING

    case( 'make_oris' )
        call cline%parse()
        call xmake_oris%execute(cline)
    case( 'orisops' )
        call cline%parse()
        call xorisops%execute(cline)
    case( 'oristats' )
        call cline%parse()
        call xoristats%execute(cline)
    case( 'vizoris' )
        call cline%parse()
        call xvizoris%execute(cline)

    ! PRINT INFO

    case( 'info_image' )
        call cline%parse()
        call xinfo_image%execute(cline)
    case( 'info_stktab' )
        call cline%parse()
        call xinfo_stktab%execute(cline)
    case( 'print_fsc' )
        call cline%parse()
        call xprint_fsc%execute(cline)
    case( 'print_magic_boxes' )
        call cline%parse()
        call xprint_magic_boxes%execute(cline)

    ! SIMULATORS

    case( 'simulate_noise' )
        call cline%parse()
        call xsimulate_noise%execute(cline)
    case( 'simulate_particles' )
        call cline%parse()
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
        call cline%parse()
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
        call cline%parse()
        call xsimulate_subtomogram%execute(cline)
    case DEFAULT
        write(*,'(a,a)') 'program key (prg) is: ', trim(prg)
        stop 'unsupported program'
    end select
end program simple_exec
