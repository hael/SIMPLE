! executes the shared-memory parallelised programs in SIMPLE
program simple_exec
include 'simple_lib.f08'
use simple_user_interface
use simple_cmdline, only: cmdline, cmdline_err
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

type(simulate_noise_commander)       :: xsimulate_noise
type(simulate_particles_commander)   :: xsimulate_particles
type(simulate_movie_commander)       :: xsimulate_movie
type(simulate_subtomogram_commander) :: xsimulate_subtomogram
type(select_commander)               :: xselect
type(make_pickrefs_commander)        :: xmake_pickrefs
type(extract_commander)              :: xextract
type(cluster_cavgs_commander)        :: xcluster_cavgs
type(mask_commander)                 :: xmask
type(info_image_commander)           :: xinfo_image
type(info_stktab_commander)          :: xinfo_stktab
type(fsc_commander)                  :: xfsc
type(center_commander)               :: xcenter
type(postprocess_commander)          :: xpostprocess
type(project_commander)              :: xproject
type(volops_commander)               :: xvolops
type(convert_commander)              :: xconvert
type(ctfops_commander)               :: xctfops
type(filter_commander)               :: xfilter
type(normalize_commander)            :: xnormalize
type(scale_commander)                :: xscale
type(stack_commander)                :: xstack
type(stackops_commander)             :: xstackops
type(print_fsc_commander)            :: xprint_fsc
type(print_magic_boxes_commander)    :: xprint_magic_boxes
type(shift_commander)                :: xshift
type(make_oris_commander)            :: xmake_oris
type(orisops_commander)              :: xorisops
type(oristats_commander)             :: xoristats
type(vizoris_commander)              :: xvizoris

! OTHER DECLARATIONS
character(len=KEYLEN) :: keys_required(MAXNKEYS)='', keys_optional(MAXNKEYS)=''
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
    case( 'select' )
        call cline%parse()
        if( .not. cline%defined('outfile') )  call cline%set('outfile', 'selected_lines.txt')
        call xselect%execute(cline)
    case( 'make_pickrefs' )
        call cline%parse()
        if( .not. cline%defined('pcontrast')) call cline%set('pcontrast', 'black')
        if( .not. cline%defined('pgrp')     ) call cline%set('pgrp',      'd1'   )
        call xmake_pickrefs%execute(cline)
    case( 'extract' )
        call cline%parse()
        if( .not. cline%defined('pcontrast') )call cline%set('pcontrast', 'black')
        call xextract%execute(cline)
    case('cluster_cavgs')
        call cline%parse()
        call xcluster_cavgs%execute(cline)
    case( 'mask' )
        call cline%parse()
        call xmask%execute(cline)
    case( 'info_image' )
        call cline%parse()
        call xinfo_image%execute(cline)
    case( 'info_stktab' )
        call cline%parse()
        call xinfo_stktab%execute(cline)
    case( 'fsc' )
        call cline%parse()
        call xfsc%execute(cline)
    case( 'center' )
        call cline%parse()
        if( .not. cline%defined('cenlp') ) call cline%set('cenlp', 30.)
        call xcenter%execute(cline)
    case( 'postprocess' )
        call cline%parse()
        call xpostprocess%execute(cline)
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
    case( 'print_fsc' )
        call cline%parse()
        call xprint_fsc%execute(cline)
    case( 'print_magic_boxes' )
        call cline%parse()
        call xprint_magic_boxes%execute(cline)
    case( 'shift' )
        call cline%parse()
        call xshift%execute(cline)
    case( 'make_oris' )
        call cline%parse()
        call xmake_oris%execute(cline)
    case( 'orisops' )
        !==Program orisops
        !
        ! <orisops/begin>is a program for analyzing SIMPLE orientation/parameter files (text files containing input parameters
        ! and/or parameters estimated by cluster2D or refine3D). If only oritab is inputted, there
        ! are a few options available. If errify=yes, then the program introduces uniform random angular errors
        ! .in.[-angerr,angerr], and uniform origin shift errors
        ! .in.[-sherr,sherr], and uniform random defocus errors .in.[-dferr,dferr]. If nstates > 1
        ! then random states are assigned .in.[1,nstates]. If mirr=2d, then the Euler angles in oritab
        ! are mirrored according to the relation e1=e1, e2=180.+e2, e3=-e3. If mirr=3d, then the Euler angles in
        ! oritab are mirrored according to the relation R=M(M*R), where R is the rotation matrix calculated from
        ! the Euler angle triplet and M is a 3D reflection matrix (like a unit matrix but with the 3,3-component sign swapped).
        ! If e1, e2, or e3 is inputted, the orientations in oritab are rotated correspondingly. If
        ! you input state as well, you rotate only the orientations assigned to state state. If mul
        ! is defined, you multiply the origin shifts with mul. If zero=yes, then the shifts are zeroed. If none of
        ! the above described parameter are defined, and oritab is still defined, the program projects the 3D orientation into
        ! the xy-plane and plots the resulting vector (this is useful for checking orientation coverage)<orisops/end>
        !
        ! set optional keys
        keys_optional(1)  = 'oritab'
        keys_optional(2)  = 'nptcls'
        keys_optional(3)  = 'outfile'
        keys_optional(4)  = 'e1'
        keys_optional(5)  = 'e2'
        keys_optional(6)  = 'e3'
        keys_optional(7)  = 'trs'
        keys_optional(8)  = 'nstates'
        keys_optional(9)  = 'pgrp'
        keys_optional(10) = 'defocus'
        keys_optional(11) = 'deftab'
        keys_optional(12) = 'angerr'
        keys_optional(13) = 'sherr'
        keys_optional(14) = 'dferr'
        keys_optional(15) = 'zero'
        keys_optional(16) = 'discrete'
        keys_optional(17) = 'ndiscrete'
        keys_optional(18) = 'state'
        keys_optional(19) = 'errify'
        keys_optional(20) = 'mul'
        keys_optional(21) = 'mirr'
        keys_optional(22) = 'xsh'
        keys_optional(23) = 'ysh'
        keys_optional(24) = 'zsh'
        keys_optional(25) = 'ctfreslim'
        keys_optional(26) = 'oritype'
        ! parse command line
        call cline%parse_oldschool(keys_optional=keys_optional(:26))
        ! execute
        call xorisops%execute(cline)
    case( 'oristats' )
        !==Program oristats
        !
        ! <oristats/begin>is a program for analyzing SIMPLE orientation/parameter files (text files containing input
        ! parameters and/or parameters estimated by cluster2D or refine3D). If two orientation
        ! tables (oritab and oritab2) are inputted, the program provides statistics of the distances between the orientations
        ! in the two documents. These statistics include the sum of angular distances between the orientations, the average
        ! angular distance between the orientations, the standard deviation of angular distances, the minimum angular
        ! distance, and the maximum angular distance<oristats/end>
        !
        ! Required keys
        keys_required(1)  = 'oritab'
        ! set optional keys
        keys_optional(1)  = 'nptcls'
        keys_optional(2)  = 'oritab2'
        keys_optional(3)  = 'outfile'
        keys_optional(4)  = 'nstates'
        keys_optional(5)  = 'state'
        keys_optional(6)  = 'ctfstats'
        keys_optional(7)  = 'trsstats'
        keys_optional(8)  = 'ncls'
        keys_optional(9)  = 'minp'
        keys_optional(10) = 'thres'
        keys_optional(11) = 'projstats'
        keys_optional(12) = 'nspace'
        keys_optional(13) = 'pgrp'
        keys_optional(14) = 'ndiscrete'
        keys_optional(15) = 'weights2D'
        keys_optional(16) = 'weights3D'
        keys_optional(17) = 'classtats'
        keys_optional(18) = 'oritype'
        ! parse command line
        call cline%parse_oldschool( keys_required(:1), keys_optional(:18) )
        ! set defaults
        if( .not. cline%defined('ndiscrete') ) call cline%set('ndiscrete', 100.)
        ! execute
        call xoristats%execute(cline)
    case( 'vizoris' )
        !==Program vizoris
        !
        ! <vizoris/begin>extract projection direction from an orientation direction
        ! for visualization in UCSF Chimera<vizoris/end>
        !
        ! Required keys
        keys_required(1) = 'oritab'
        ! set optional keys
        keys_optional(1) = 'nspace'
        keys_optional(2) = 'pgrp'
        keys_optional(3) = 'tseries'
        keys_optional(4) = 'oritype'
        ! parse command line
        call cline%parse_oldschool( keys_required(:1), keys_optional(:4) )
        ! execute
        call xvizoris%execute(cline)

    case DEFAULT
        write(*,'(a,a)') 'program key (prg) is: ', trim(prg)
        stop 'unsupported program'
    end select

end program simple_exec
