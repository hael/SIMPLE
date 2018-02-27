! executes the shared-memory parallelised programs in SIMPLE
program simple_exec
#include "simple_lib.f08"
use simple_cmdline, only: cmdline, cmdline_err
use simple_gen_doc
use simple_commander_checks
use simple_commander_comlin
use simple_commander_distr
use simple_commander_imgproc
use simple_commander_mask
use simple_commander_misc
use simple_commander_oris
use simple_commander_preprocess
use simple_commander_prime2D
use simple_commander_prime3D
use simple_commander_rec
use simple_commander_sim
use simple_commander_volops
use simple_commander_tseries
implicit none

! SIMULATOR PROGRAMS
type(simulate_noise_commander)       :: xsimulate_noise
type(simulate_particles_commander)   :: xsimulate_particles
type(simulate_movie_commander)       :: xsimulate_movie
type(simulate_subtomogram_commander) :: xsimulate_subtomogram

! PRE-PROCESSING PROGRAMS
type(select_commander)               :: xselect
type(make_pickrefs_commander)        :: xmake_pickrefs
type(extract_commander)              :: xextract

! CLASS AVERAGE ANALYSIS
type(cluster_cavgs_commander)        :: xcluster_cavgs

! MASK PROGRAMS
type(mask_commander)                 :: xmask

! CHECKER PROGRAMS
type(info_image_commander)           :: xinfo_image
type(info_stktab_commander)          :: xinfo_stktab

! VOLOPS PROGRAMS
type(fsc_commander)                  :: xfsc
type(center_commander)               :: xcenter
type(postprocess_commander)          :: xpostprocess
type(project_commander)              :: xproject
type(volops_commander)               :: xvolops

! GENERAL IMAGE PROCESSING PROGRAMS
type(convert_commander)              :: xconvert
type(ctfops_commander)               :: xctfops
type(filter_commander)               :: xfilter
type(normalize_commander)            :: xnormalize
type(scale_commander)                :: xscale
type(stack_commander)                :: xstack
type(stackops_commander)             :: xstackops

! MISCELLANOUS PROGRAMS
type(print_cmd_dict_commander)       :: xprint_cmd_dict
type(print_fsc_commander)            :: xprint_fsc
type(print_magic_boxes_commander)    :: xprint_magic_boxes
type(shift_commander)                :: xshift

! ORIENTATION DATA MANAGEMENT PROGRAMS
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
logical               :: describe

! parse command-line
call get_command_argument(1, xarg, cmdlen, cmdstat)
call get_command(entire_line)
if( str_has_substr(entire_line, 'prg=list') ) call list_all_simple_programs
describe = str_has_substr(entire_line, 'describe=yes')
pos = index(xarg, '=') ! position of '='
call cmdline_err( cmdstat, cmdlen, xarg, pos )
prg = xarg(pos+1:)     ! this is the program name
if( str_has_substr(prg, 'simple') ) stop 'giving program names with simple* prefix is depreciated'
select case(prg)

    ! SIMULATOR PROGRAMS

    case( 'simulate_noise' )
        !==Program simulate_noise
        !
        ! <simulate_noise/begin>is a program for generating pure noise images<simulate_noise/end>
        !
          ! set required keys
        keys_required(1) = 'box'
        keys_required(2) = 'nptcls'
        ! parse command line
        if( describe ) call print_doc_simulate_noise
        call cline%parse(keys_required(:2))
        ! execute
        call xsimulate_noise%execute(cline)
    case( 'simulate_particles' )
        !==Program simulate_particles
        !
        ! <simulate_particles/begin>is a program for simulating single-particle cryo-EM images. It is not a very
        ! sophisticated simulator, but it is nevertheless useful for testing purposes. It does not do any multi-
        ! slice simulation and it cannot be used for simulating molecules containing heavy atoms. It does not even
        ! accept a PDB file as an input. Input is a cryo-EM map, which we usually generate from a PDB file using
        ! EMANs program pdb2mrc. The volume is projected using Fourier interpolation, 20% of the total noise is
        ! added to the images (pink noise), they are then Fourier transformed and multiplied with astigmatic CTF and
        ! B-factor. Next, the they are inverse FTed before the remaining 80% of the noise (white noise) is added
        ! <simulate_particles/end>
        !
        ! set required keys
        keys_required(1)  = 'vol1'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'msk'
        keys_required(4)  = 'nptcls'
        keys_required(5)  = 'snr'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'xsh'
        keys_optional(3)  = 'ysh'
        keys_optional(4)  = 'sherr'
        keys_optional(5)  = 'oritab'
        keys_optional(6)  = 'outfile'
        keys_optional(7)  = 'outstk'
        keys_optional(8)  = 'single'
        keys_optional(9)  = 'ndiscrete'
        keys_optional(10) = 'diverse'
        keys_optional(11) = 'pgrp'
        ! set optional CTF-related keys
        keys_optional(12) = 'ctf'
        keys_optional(13) = 'kv'
        keys_optional(14) = 'cs'
        keys_optional(15) = 'fraca'
        keys_optional(16) = 'deftab'
        keys_optional(17) = 'defocus'
        keys_optional(18) = 'dferr'
        keys_optional(19) = 'astigerr'
        keys_optional(20) = 'bfac'
        keys_optional(21) = 'bfacerr'
        ! parse command line
        if( describe ) call print_doc_simulate_particles
        call cline%parse(keys_required(:5), keys_optional(:21))
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
        ! execute
        call xsimulate_particles%execute(cline)
    case( 'simulate_movie' )
        !==Program simulate_movie
        !
        ! <simulate_movie/begin>is a program for crude simulation of a DDD movie. Input is a set of projection images to place.
        ! Movie frames are then generated related by randomly shifting the base image and applying noise<simulate_movie/end>
        !
        ! set required keys
        keys_required(1)  = 'stk'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'msk'
        keys_required(4)  = 'xdim'
        keys_required(5)  = 'ydim'
        keys_required(6)  = 'snr'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'nframes'
        keys_optional(3)  = 'trs'
        keys_optional(4)  = 'vis'
        ! set optional CTF-related keys
        keys_optional(5)  = 'kv'
        keys_optional(6)  = 'cs'
        keys_optional(7)  = 'fraca'
        keys_optional(8)  = 'deftab'
        keys_optional(9)  = 'defocus'
        keys_optional(10) = 'bfac'
        ! parse command line
        if( describe ) call print_doc_simulate_movie
        call cline%parse(keys_required(:6), keys_optional(:10))
        ! set defaults
        if( .not. cline%defined('trs')     ) call cline%set('trs',      3.)
        if( .not. cline%defined('ctf')     ) call cline%set('ctf',   'yes')
        if( .not. cline%defined('bfac')    ) call cline%set('bfac',   200.)
        if( .not. cline%defined('nframes') ) call cline%set('nframes', 30.)
        call cline%set('wfun', 'kb')
        call cline%set('winsz', 1.5)
        call cline%set('alpha',  2.)
        call cline%set('eo',   'no')
        ! execute
        call xsimulate_movie%execute(cline)
    case( 'simulate_subtomogram' )
        !==Program simulate_subtomogram
        !
        ! <simulate_subtomogram/begin>is a program for crude simulation of a subtomogram<simulate_subtomogram/end>
        !
        ! set required keys
        keys_required(1) = 'vol1'
        keys_required(2) = 'smpd'
        keys_required(3) = 'nptcls'
        keys_required(4) = 'snr'
        ! set optional keys
        keys_optional(1) = 'nthr'
        ! parse command line
        if( describe ) call print_doc_simulate_subtomogram
        call cline%parse(keys_required(:4), keys_optional(:1))
        ! execute
        call xsimulate_subtomogram%execute(cline)
    case( 'select' )
        !==Program select
        !
        ! <select/begin>is a program for selecting files based on image correlation matching<select/end>
        !
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'stk2'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'stk3'
        keys_optional(3) = 'filetab'
        keys_optional(4) = 'outfile'
        keys_optional(5) = 'outstk'
        keys_optional(6) = 'dir_select'
        keys_optional(7) = 'dir_reject'
        ! parse command line
        if( describe ) call print_doc_select
        call cline%parse(keys_required(:2), keys_optional(:7))
        ! set defaults
        if( .not. cline%defined('outfile') )  call cline%set('outfile', 'selected_lines.txt')
        ! execute
        call xselect%execute(cline)
    case( 'make_pickrefs' )
        !==Program pickrefs
        !
        ! <make_pickrefs/begin>is a program for generating references for template-based particle picking<make_pickrefs/end>
        !
        ! set required keys
        keys_required(1) = 'pgrp'
        keys_required(2) = 'smpd'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'vol1'
        keys_optional(3) = 'stk'
        keys_optional(4) = 'pcontrast'
        ! parse command line
        if( describe ) call print_doc_make_pickrefs
        call cline%parse(keys_required(:1), keys_optional(:4))
        ! set defaults
        if( .not. cline%defined('pcontrast')) call cline%set('pcontrast', 'black')
        if( .not. cline%defined('pgrp')     ) call cline%set('pgrp',      'd1'   )
        ! execute
        call xmake_pickrefs%execute(cline)
    case( 'extract' )
        !==Program extract
        !
        ! <extract/begin>is a program that extracts particle images from integrated movies.
        ! Boxfiles are assumed to be in EMAN format but we provide a conversion script
        ! (relion2emanbox.pl) for *.star files containing particle coordinates obtained
        ! with Relion. In addition to one single-particle image stack per micrograph, the
        ! program produces a parameter files that should be concatenated for use in
        ! conjunction with other SIMPLE programs.
        ! <extract/end>
        !
        ! set required keys
        keys_required(1) = 'smpd'
        ! set optional keys
        keys_optional(1) = 'unidoc'
        keys_optional(2) = 'filetab'
        keys_optional(3) = 'boxtab'
        keys_optional(4) = 'ctffind_doc'
        keys_optional(5) = 'pcontrast'
        keys_optional(6) = 'box'
        keys_optional(7) = 'outside'
        keys_optional(8) = 'dir'
        ! parse command line
        if( describe ) call print_doc_extract
        call cline%parse(keys_required(:1),keys_optional(:8))
        ! parse command line
        if( .not. cline%defined('pcontrast') )call cline%set('pcontrast', 'black')
        ! execute
        call xextract%execute(cline)

    ! PRIME2D PROGRAMS

    case('cluster_cavgs')
        !==Program cluster_cavgs
        !
        ! <cluster_cavgs/begin>is a program for analyzing class averages with affinity propagation,
        ! in order to get a better understanding of the view distribution. The balance flag is used
        ! to apply a balancing restraint (on the class population). Adjust balance until you are
        ! satisfied with the shape of the histogram. <cluster_cavgs/end>
        !
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'smpd'
        keys_required(3) = 'msk'
        keys_required(4) = 'lp'
        keys_required(5) = 'classdoc'
        ! set optional keys
        keys_optional(1) = 'hp'
        keys_optional(2) = 'nthr'
        keys_optional(3) = 'balance'
        keys_optional(4) = 'objfun'
        ! parse command line
        ! if( describe ) call print_doc_cluster_cavgs
        call cline%parse(keys_required(:5), keys_optional(:4))
        ! execute
        call xcluster_cavgs%execute(cline)

    ! MASK PROGRAMS

    case( 'mask' )
        !==Program mask
        !
        ! <mask/begin>is a program for masking images and volumes.
        ! If you want to mask your images with a spherical mask with a soft falloff, set msk
        ! to the radius in pixels<mask/end>
        !
        ! set required keys
        keys_required(1)  = 'smpd'
        ! set optional keys
        keys_optional(1)  = 'stk'
        keys_optional(2)  = 'vol1'
        keys_optional(3)  = 'msktype'
        keys_optional(4)  = 'inner'
        keys_optional(5)  = 'width'
        keys_optional(6)  = 'outer'
        keys_optional(7)  = 'nthr'
        keys_optional(8)  = 'mw'
        keys_optional(9)  = 'edge'
        keys_optional(10) = 'amsklp'
        keys_optional(11) = 'automsk'
        keys_optional(12) = 'smpd'
        keys_optional(13) = 'outstk'
        keys_optional(14) = 'outvol'
        keys_optional(15) = 'mskfile'
        keys_optional(16) = 'taper_edges'
        keys_optional(17) = 'msk'
        keys_optional(18) = 'pdbfile'
        keys_optional(19) = 'oritab'
        keys_optional(20) = 'outfile'
        keys_optional(21) = 'center'
        ! parse command line
        if( describe ) call print_doc_mask
        call cline%parse( keys_required(:1), keys_optional(:21))
        ! execute
        call xmask%execute(cline)
    case( 'info_image' )
        !==Program info_image
        !
        ! <info_image/begin>is a program for printing header information in MRC and SPIDER stacks
        ! and volumes<info_image/end>
        !
        ! set required keys
        keys_required(1) = 'fname'
        ! set optional keys
        keys_optional(1) = 'box'
        keys_optional(2) = 'smpd'
        keys_optional(3) = 'stats'
        keys_optional(4) = 'endian'
        ! parse command line
        if( describe ) call print_doc_info_image
        call cline%parse(keys_required(:1), keys_optional(:4))
        ! execute
        call xinfo_image%execute(cline)
    case( 'info_stktab' )
        !==Program info_stktab
        !
        ! <info_stktab/begin>is a program for for printing information about stktab <info_stktab/end>
        !
        ! set required keys
        keys_required(1) = 'stktab'
        ! parse command line
        ! if( describe ) call print_doc_info_stktab
        call cline%parse(keys_required(:1))
        ! execute
        call xinfo_stktab%execute(cline)

    ! VOLOPS PROGRAMS

    case( 'fsc' )
        !==Program fsc
        !
        ! <fsc/begin>is a program for calculating the FSC between the two input volumes<fsc/end>
        !
        ! set required keys
        keys_required(1) = 'smpd'
        keys_required(2) = 'vol1'
        keys_required(3) = 'vol2'
        keys_required(4) = 'msk'
        ! set optional keys
        keys_optional(1) = 'mskfile'
        ! parse command line
        if( describe ) call print_doc_fsc
        call cline%parse(keys_required(:4), keys_optional(:1))
        ! execute
        call xfsc%execute(cline)
    case( 'center' )
        !==Program center
        !
        ! <center/begin>is a program for centering a volume and mapping the shift parameters
        ! back to the particle images<center/end>
        !
        ! set required keys
        keys_required(1) = 'vol1'
        keys_required(2) = 'smpd'
        ! set optional keys
        keys_optional(1) = 'oritab'
        keys_optional(2) = 'outfile'
        keys_optional(3) = 'cenlp'
        ! parse command line
        if( describe ) call print_doc_center
        call cline%parse(keys_required(:2), keys_optional(:3))
        ! set defaults
        if( .not. cline%defined('cenlp') ) call cline%set('cenlp', 30.)
        ! execute
        call xcenter%execute(cline)
    case( 'postprocess' )
        !==Program postprocess
        !
        ! <postprocess/begin>is a program for post-processing of volumes.
        ! Use program volops to estimate the B-factor with the Guinier plot
        ! <postprocess/end>
        !
        ! set required keys
        keys_required(1)  = 'vol1'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'msk'
        ! set optional keys
        keys_optional(1)  = 'fsc'
        keys_optional(2)  = 'lp'
        keys_optional(3)  = 'mw'
        keys_optional(4)  = 'bfac'
        keys_optional(5)  = 'automsk'
        keys_optional(6)  = 'amsklp'
        keys_optional(7)  = 'edge'
        keys_optional(8)  = 'binwidth'
        keys_optional(9)  = 'thres'
        keys_optional(10) = 'mskfile'
        keys_optional(11) = 'vol_filt'
        keys_optional(12) = 'inner'
        keys_optional(13) = 'mirr'
        ! parse command line
        if( describe ) call print_doc_postprocess
        call cline%parse(keys_required(:3), keys_optional(:13))
        ! execute
        call xpostprocess%execute(cline)
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
        if( describe ) call print_doc_project
        call cline%parse(keys_required(:2), keys_optional(:10))
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
        if( describe ) call print_doc_volops
        call cline%parse(keys_optional=keys_optional(:20))
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
        if( describe ) call print_doc_convert
        call cline%parse(keys_optional=keys_optional(:4))
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
        if( describe ) call print_doc_ctfops
        call cline%parse(keys_required(:2), keys_optional(:6))
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
        if( describe ) call print_doc_filter
        call cline%parse(keys_required(:1), keys_optional(:11))
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
        ! if( describe ) call print_doc_normalize
        call cline%parse(keys_required(:1), keys_optional(:6))
        ! execute
        call xnormalize%execute(cline)
    case( 'scale' )
        !==Program scale
        !
        ! <scale/begin>is a program that provides re-scaling and clipping routines for MRC or SPIDER stacks
        ! and volumes<scale/end>
        !
        ! set required keys
        keys_required(1)  = 'smpd'
        ! set optional keys
        keys_optional(1)  = 'stk'
        keys_optional(2)  = 'vol1'
        keys_optional(3)  = 'filetab'
        keys_optional(4)  = 'msk'
        keys_optional(5)  = 'scale'
        keys_optional(6)  = 'scale2'
        keys_optional(7)  = 'newbox'
        keys_optional(8)  = 'clip'
        keys_optional(9)  = 'outvol'
        keys_optional(10) = 'outstk'
        keys_optional(11) = 'outstk2'
        ! parse command line
        if( describe ) call print_doc_scale
        call cline%parse(keys_required(:1),keys_optional(:11))
        ! execute
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
        if( describe ) call print_doc_stack
        call cline%parse(keys_required(:3), keys_optional(:6))
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
        if( describe ) call print_doc_stackops
        call cline%parse( keys_required(:1),keys_optional(:25) )
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
        if( describe ) call print_doc_print_cmd_dict
        call cline%parse(keys_optional=keys_optional(:1))
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
        if( describe ) call print_doc_print_fsc
        call cline%parse(keys_required(:3))
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
        if( describe ) call print_doc_print_magic_boxes
        call cline%parse(keys_required(:2))
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
        if( describe ) call print_doc_shift
        call cline%parse(keys_required(:3), keys_optional(:3))
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
        if( describe ) call print_doc_make_deftab
        call cline%parse(keys_required(:5),keys_optional(:8))
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
        if( describe ) call print_doc_make_oris
        call cline%parse(keys_required(:1),keys_optional(:21))
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
        if( describe ) call print_doc_map2ptcls
        call cline%parse(keys_required(:3), keys_optional(:7))
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
        if( describe ) call print_doc_orisops
        call cline%parse(keys_optional=keys_optional(:26))
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
        if( describe ) call print_doc_oristats
        call cline%parse( keys_required(:1), keys_optional(:18) )
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
        if( describe ) call print_doc_vizoris
        call cline%parse( keys_required(:1), keys_optional(:4) )
        ! execute
        call xvizoris%execute(cline)
    case DEFAULT
        write(*,'(a,a)') 'program key (prg) is: ', trim(prg)
        stop 'unsupported program'
    end select

end program simple_exec
