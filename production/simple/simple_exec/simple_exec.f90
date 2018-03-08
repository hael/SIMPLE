! executes the shared-memory parallelised programs in SIMPLE
program simple_exec
include 'simple_lib.f08'

use simple_user_interface, only: make_user_interface
use simple_cmdline,        only: cmdline, cmdline_err
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
type(print_cmd_dict_commander)       :: xprint_cmd_dict
type(print_fsc_commander)            :: xprint_fsc
type(print_magic_boxes_commander)    :: xprint_magic_boxes
type(shift_commander)                :: xshift
type(make_deftab_commander)          :: xmake_deftab
type(make_oris_commander)            :: xmake_oris
type(map2ptcls_commander)            :: xmap2ptcls
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
        ! set defaults
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
        ! set defaults
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
        ! set defaults
        if( .not. cline%defined('outfile') )  call cline%set('outfile', 'selected_lines.txt')
        call xselect%execute(cline)
    case( 'make_pickrefs' )
        call cline%parse()
        if( .not. cline%defined('pcontrast')) call cline%set('pcontrast', 'black')
        if( .not. cline%defined('pgrp')     ) call cline%set('pgrp',      'd1'   )
        ! execute
        call xmake_pickrefs%execute(cline)
    case( 'extract' )
        call cline%parse()
        ! set defaults
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
        ! set defaults
        if( .not. cline%defined('cenlp') ) call cline%set('cenlp', 30.)
        call xcenter%execute(cline)
    case( 'postprocess' )
        call cline%parse()
        ! execute
        call xpostprocess%execute(cline)

!!!!!!!!!!!!!!!!!!!!!!

    case( 'project' )
        !==Program project
        !
        ! <project/begin>is a program for projecting a volume using interpolation in Fourier space. Input is a SPIDER or
        ! MRC volume. Output is a stack of projection images of the same format as the inputted volume. Projections
        ! are generated by extraction of central sections from the Fourier volume and back transformation of the 2D FTs.
        ! nspace controls the number of projection images generated with quasi-even projection directions. The
        ! oritab parameter allows you to input the orientations that you wish to have your volume projected in. If
        ! rnd=yes, random rather than quasi-even projections are generated, trs then controls the halfwidth of
        ! the random origin shift. Less commonly used parameters are pgrp, which controls the point-group symmetry
        ! c (rotational), d (dihedral), t (tetrahedral), o (octahedral) or i (icosahedral). The point-group symmetry is
        ! used to restrict the set of projections to within the asymmetric unit.
        ! neg inverts the contrast of the projections. <project/end>
        !
        ! set required keys
        keys_required(1)  = 'vol1'
        keys_required(2)  = 'smpd'
        ! set optional keys
        keys_optional(1)  = 'nspace'
        keys_optional(2)  = 'outstk'
        keys_optional(3)  = 'oritab'
        keys_optional(4)  = 'nthr'
        keys_optional(5)  = 'rnd'
        keys_optional(6)  = 'trs'
        keys_optional(7)  = 'pgrp'
        keys_optional(8)  = 'neg'
        keys_optional(9)  = 'top'
        keys_optional(10) = 'msk'
        ! parse command line
        call cline%parse_oldschool(keys_required(:2), keys_optional(:10))
        ! set defaults
        if( .not. cline%defined('wfun')  ) call cline%set('wfun', 'kb')
        if( .not. cline%defined('winsz') ) call cline%set('winsz', 1.5)
        if( .not. cline%defined('alpha') ) call cline%set('alpha', 2.)
        ! execute
        call xproject%execute(cline)
    case( 'volops' )
        !==Program volops
        !
        ! <volops/begin>provides standard single-particle image processing routines that are applied to MRC or SPIDER volumes
        ! <volops/end>
        !
        ! set optional keys
        keys_optional(1)  = 'vol1'
        keys_optional(2)  = 'nthr'
        keys_optional(3)  = 'guinier'
        keys_optional(4)  = 'smpd'
        keys_optional(5)  = 'hp'
        keys_optional(6)  = 'lp'
        keys_optional(7)  = 'msk'
        keys_optional(8)  = 'neg'
        keys_optional(9)  = 'snr'
        keys_optional(10) = 'mirr'
        keys_optional(11) = 'bfac'
        keys_optional(12) = 'e1'
        keys_optional(13) = 'e2'
        keys_optional(14) = 'e3'
        keys_optional(15) = 'xsh'
        keys_optional(16) = 'ysh'
        keys_optional(17) = 'zsh'
        keys_optional(18) = 'outvol'
        keys_optional(19) = 'vollist'
        keys_optional(20) = 'outfile'
        ! parse command line
        call cline%parse_oldschool(keys_optional=keys_optional(:20))
        if( .not.cline%defined('vol1') .and. .not.cline%defined('vollist') )&
            &stop 'Input volume required!'
        ! execute
        call xvolops%execute(cline)

    ! GENERAL IMAGE PROCESSING PROGRAMS

    case( 'convert' )
        !==Program convert
        !
        ! <convert/begin>is a program for converting between SPIDER and MRC formats
        ! <convert/end>
        !
        ! set optional keys
        keys_optional(1) = 'stk'
        keys_optional(2) = 'vol1'
        keys_optional(3) = 'outstk'
        keys_optional(4) = 'outvol'
        ! parse command line
        call cline%parse_oldschool(keys_optional=keys_optional(:4))
        ! execute
        call xconvert%execute(cline)
    case( 'ctfops' )
        !==Program ctfops
        !
        ! <ctfops/begin>is a program for applying CTF to stacked images
        ! <ctfops/end>
        !
        ! Required keys
        keys_required(1) = 'smpd'
        keys_required(2) = 'ctf'
        ! Optional keys
        keys_optional(1) = 'stk'
        keys_optional(2) = 'outstk'
        keys_optional(3) = 'neg'
        keys_optional(4) = 'oritab'
        keys_optional(5) = 'deftab'
        keys_optional(6) = 'bfac'
        ! parse command line
        call cline%parse_oldschool(keys_required(:2), keys_optional(:6))
        ! set defaults
        if( .not. cline%defined('stk') ) call cline%set('box', 256.)
        ! execute
        call xctfops%execute(cline)
    case( 'filter' )
        !==Program filter
        !
        ! <filter/begin>is a program for filtering stacks/volumes
        ! <filter/end>
        !
        ! Required keys
        keys_required(1)  = 'smpd'
        ! Optional keys
        keys_optional(1)  = 'stk'
        keys_optional(2)  = 'vol1'
        keys_optional(3)  = 'outstk'
        keys_optional(4)  = 'outvol'
        keys_optional(5)  = 'lp'
        keys_optional(6)  = 'hp'
        keys_optional(7)  = 'phrand'
        keys_optional(8)  = 'bfac'
        keys_optional(9)  = 'winsz'
        keys_optional(10) = 'real_filter'
        keys_optional(11) = 'width'
        ! parse command line
        call cline%parse_oldschool(keys_required(:1), keys_optional(:11))
        ! execute
        call xfilter%execute(cline)
    case( 'normalize' )
        !==Program normalize
        !
        ! <normalize/begin>is a program for normalization of MRC or SPIDER stacks and volumes. If you want to
        ! normalize your images inputted with stk, set norm=yes. If you want to perform noise normalisation
        ! of the images set noise_norm=yes given a mask radius msk (pixels). If you want to normalize your
        ! images or volume (vol1) with respect to their power spectrum set shell_norm=yes
        ! <normalize/end>
        !
        ! set required keys
        keys_required(1) = 'smpd'
        ! set optional keys
        keys_optional(1) = 'stk'
        keys_optional(2) = 'msk'
        keys_optional(3) = 'norm'
        keys_optional(4) = 'noise_norm'
        keys_optional(5) = 'shell_norm'
        keys_optional(6) = 'nthr'
        ! parse command line
        call cline%parse_oldschool(keys_required(:1), keys_optional(:6))
        ! execute
        call xnormalize%execute(cline)
    case( 'scale' )
        call cline%parse()
        call xscale%execute(cline)
    case( 'stack' )
        !==Program stack
        !
        ! <stack/begin>is a program for stacking individual images or multiple stacks into one<stack/end>
        !
        ! set required keys
        keys_required(1) = 'filetab'
        keys_required(2) = 'outstk'
        keys_required(3) = 'smpd'
        ! set optional keys
        keys_optional(1) = 'clip'
        keys_optional(2) = 'nframes'
        keys_optional(3) = 'fbody'
        keys_optional(4) = 'numlen'
        keys_optional(5) = 'xdim'
        keys_optional(6) = 'ydim'
        keys_optional(7) = 'endian'
        ! parse command line
        call cline%parse_oldschool(keys_required(:3), keys_optional(:6))
        ! execute
        call xstack%execute(cline)
    case( 'stackops' )
        !==Program stackops
        !
        ! <stackops/begin>is a program that provides standard single-particle image processing routines that are applied to MRC or SPIDER
        ! stacks. If you want to extract a particular state, give an alignment document (oritab) and set state
        ! to the state that you want to extract. If you want to select the fraction of best particles (according to the goal function), input
        ! an alignment doc (oritab) and set frac. You can combine the state and frac options. If you
        ! want to apply noise to images, give the desired signal-to-noise ratio via snr. If you want to calculate the autocorrelation
        ! function of your images set acf=yes. If you want to extract a contiguous subset of particle images from the stack, set
        ! fromp and top. If you want to fish out a number of particle images from your stack at random, set nran to
        ! some nonzero integer number less than nptcls. With avg=yes the global average of the inputted stack is calculated.
        ! If you define nframesgrp to some integer number larger than one averages with chunk sizes of nframesgrp are produced,
        ! which may be useful for analysis of dose-fractionated image series. neg inverts the contrast of the images<stackops/end>
        !
        ! Required keys
        keys_required(1)  = 'smpd'
        ! set optional keys
        keys_optional(1)  = 'oritab'
        keys_optional(2)  = 'outstk'
        keys_optional(3)  = 'mirr'
        keys_optional(4)  = 'nran'
        keys_optional(5)  = 'frac'
        keys_optional(6)  = 'state'
        keys_optional(7)  = 'class'
        keys_optional(8)  = 'neg'
        keys_optional(9)  = 'acf'
        keys_optional(10) = 'avg'
        keys_optional(11) = 'nframesgrp'
        keys_optional(12) = 'vis'
        keys_optional(13) = 'snr'
        keys_optional(14) = 'fromp'
        keys_optional(15) = 'top'
        keys_optional(16) = 'nptcls'
        keys_optional(17) = 'order'
        keys_optional(18) = 'bfac'
        keys_optional(19) = 'outfile'
        keys_optional(20) = 'ctfreslim'
        keys_optional(21) = 'dfclose'
        keys_optional(22) = 'dffar'
        keys_optional(23) = 'stats'
        keys_optional(24) = 'stk'
        keys_optional(25) = 'stktab'
        ! parse command line
        call cline%parse_oldschool( keys_required(:1),keys_optional(:25) )
        ! sanity check
        if( cline%defined('stk') .or. cline%defined('stktab') )then
            ! all ok
        else
            stop 'stk or stktab need to be part of command line!'
        endif
        if( .not. cline%defined('outfile') ) call cline%set('outfile', 'outfile.txt')
        ! execute
        call xstackops%execute(cline)

    ! MISCELLANOUS PROGRAMS

    case( 'print_cmd_dict' )
        !==Program print_cmd_dict
        !
        ! <print_cmd_dict/begin>is a program for printing the command line key dictonary<print_cmd_dict/end>
        !
        ! set optional keys
        keys_optional(1) = 'outfile'
        ! parse command line
        call cline%parse_oldschool(keys_optional=keys_optional(:1))
        ! execute
        call xprint_cmd_dict%execute(cline)
    case( 'print_fsc' )
        !==Program print_fsc
        !
        ! <print_fsc/begin>is a program for printing the binary FSC files produced by PRIME3D<print_fsc/end>
        !
        ! set required keys
        keys_required(1)  = 'smpd'
        keys_required(2)  = 'box'
        keys_required(3)  = 'fsc'
        ! parse command line
        call cline%parse_oldschool(keys_required(:3))
        ! execute
        call xprint_fsc%execute(cline)
    case( 'print_magic_boxes' )
        !==Program print_magic_boxes
        !
        ! <print_magic_boxes/begin>is a program for printing magic box sizes (fast FFT)<print_magic_boxes/end>
        !
        ! set required keys
        keys_required(1) = 'smpd'
        keys_required(2) = 'moldiam'
        ! parse command line
        call cline%parse_oldschool(keys_required(:2))
        ! execute
        call xprint_magic_boxes%execute(cline)
    case( 'shift' )
        !==Program shift
        !
        ! <shift/begin>is a program for shifting a stack according to shifts in oritab<shift/end>
        !
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'smpd'
        keys_required(3) = 'oritab'
        ! set optional keys
        keys_optional(1) = 'outstk'
        keys_optional(2) = 'mul'
        keys_optional(3) = 'oritype'
        ! parse command line
        call cline%parse_oldschool(keys_required(:3), keys_optional(:3))
        ! execute
        call xshift%execute(cline)

    ! ORIENTATION DATA MANAGEMENT PROGRAMS

    case( 'make_deftab' )
        !==Program make_deftab
        !
        ! <make_deftab/begin>is a program for creating a SIMPLE conformant file of CTF parameter values (deftab).
        ! Input is either an earlier SIMPLE deftab/oritab. The purpose is to get the kv, cs, and fraca parameters
        ! as part of the CTF input doc as that is the new convention. The other alternative is to input a plain text
        ! file with CTF parameters dfx, dfy, angast, phshift according to the Frealign convention. Unit conversions are dealt
        ! with using optional variables. The units refer to the units in the inputted document<make_deftab/end>
        !
        ! Required keys
        keys_required(1) = 'smpd'
        keys_required(2) = 'kv'
        keys_required(3) = 'cs'
        keys_required(4) = 'fraca'
        keys_required(5) = 'outfile'
        ! set optional keys
        keys_optional(1) = 'plaintexttab'
        keys_optional(2) = 'oritab'
        keys_optional(3) = 'deftab'
        keys_optional(4) = 'dfunit'
        keys_optional(5) = 'angastunit'
        keys_optional(6) = 'phaseplate'
        keys_optional(7) = 'phshiftunit'
        keys_optional(8) = 'oritype'
        ! parse command line
        call cline%parse_oldschool(keys_required(:5),keys_optional(:8))
        ! execute
        call xmake_deftab%execute(cline)
    case( 'make_oris' )
        !==Program make_oris
        !
        ! <make_oris/begin>is a program for making SIMPLE orientation/parameter files (text files containing input parameters and/or
        ! parameters estimated by cluster2D or refine3D). The program generates random
        ! Euler angles e1.in.[0,360], e2.in.[0,180], and e3.in.[0,360] and random origin
        ! shifts x.in.[-trs,yrs] and y.in.[-trs,yrs]. If ndiscrete is set to an integer number > 0, the
        ! orientations produced are randomly sampled from the set of ndiscrete quasi-even projection directions, and the in-plane
        ! parameters are assigned randomly. If even=yes, then all nptcls orientations are assigned
        ! quasi-even projection directions,and random in-plane parameters. If nstates is set to some integer number > 0, then
        ! states are assigned randomly .in.[1,nstates]. If zero=yes in this mode of execution, the projection
        ! directions are zeroed and only the in-plane parameters are kept intact. If errify=yes and astigerr is defined,
        ! then uniform random astigmatism errors are introduced .in.[-astigerr,astigerr]<make_oris/end>
        !
        ! Required keys
        keys_required(1)  = 'nptcls'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'minp'
        keys_optional(3)  = 'ncls'
        keys_optional(4)  = 'outfile'
        keys_optional(5)  = 'trs'
        keys_optional(6)  = 'nstates'
        keys_optional(7)  = 'pgrp'
        keys_optional(8)  = 'ctf'
        keys_optional(9)  = 'defocus'
        keys_optional(10) = 'angerr'
        keys_optional(11) = 'sherr'
        keys_optional(12) = 'dferr'
        keys_optional(13) = 'even'
        keys_optional(14) = 'zero'
        keys_optional(15) = 'discrete'
        keys_optional(16) = 'ndiscrete'
        keys_optional(17) = 'diverse'
        keys_optional(18) = 'state'
        keys_optional(19) = 'nspace'
        keys_optional(20) = 'iares'
        keys_optional(21) = 'oritype'
        ! parse command line
        call cline%parse_oldschool(keys_required(:1),keys_optional(:21))
        ! execute
        call xmake_oris%execute(cline)
    case( 'map2ptcls' )
        !==Program map2ptcls
        !
        ! <map2ptcls/begin>is a program for mapping parameters that have been obtained using class averages to
        ! the individual particle images<map2ptcls/end>
        !
        ! set required keys
        keys_required(1) = 'stk2'
        keys_required(2) = 'stk3'
        keys_required(3) = 'oritab'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'oritab3D'
        keys_optional(3) = 'deftab'
        keys_optional(4) = 'outfile'
        keys_optional(5) = 'mul'
        keys_optional(6) = 'stk'
        keys_optional(7) = 'stktab'
        ! parse command line
        call cline%parse_oldschool(keys_required(:3), keys_optional(:7))
        ! set defaults
        if( .not. cline%defined('outfile') ) call cline%set('outfile', 'mapped_ptcls_params.txt')
        ! execute
        call xmap2ptcls%execute(cline)
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
