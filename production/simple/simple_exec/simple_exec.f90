!==Program simple_exec
!
! <simple_exec/begin> executes every individual program in the SIMPLE library. To list all programs type
! simple_exec prg=list. To execute program prime3D type simple_exec prg=prime3D and instructions are given. 
! simple_exec is the client of the commanders that create concrete commander objects (xsimimgs etc.) and 
! prepares their reciever (cmdline). In this case, simple_exec is also the invoker as it executes the request. 
! <simple_exec/end>
!
! The code is distributed with the hope that it will be useful, but WITHOUT ANY WARRANTY.
! Redistribution and modification is regulated by the GNU General Public License.
! Authors: Cyril Reboul & Hans Elmlund 2016
!
program simple_exec
use simple_defs
use simple_cmdline, only: cmdline
use simple_strings, only: str_has_substr
use simple_jiffys,  only: cmdline_err
use simple_gen_doc
use simple_commander_checks
use simple_commander_comlin
use simple_commander_distr
use simple_commander_imgproc
use simple_commander_mask
use simple_commander_misc
use simple_commander_oris
use simple_commander_preproc
use simple_commander_prime2D
use simple_commander_prime3D
use simple_commander_rec
use simple_commander_sim
use simple_commander_volops
use simple_commander_tseries
implicit none

! SIMULATOR PROGRAMS
type(noiseimgs_commander)          :: xnoiseimgs
type(simimgs_commander)            :: xsimimgs
type(simmovie_commander)           :: xsimmovie
type(simsubtomo_commander)         :: xsimsubtomo

! PRE-PROCESSING PROGRAMS
type(preproc_commander)            :: xpreproc
type(select_frames_commander)      :: xselect_frames
type(boxconvs_commander)           :: xboxconvs
type(integrate_movies_commander)   :: xintegrate_movies
type(powerspecs_commander)         :: xpowerspecs
type(unblur_commander)             :: xunblur
type(ctffind_commander)            :: xctffind
type(select_commander)             :: xselect
type(makepickrefs_commander)       :: xmakepickrefs
type(pick_commander)               :: xpick
type(extract_commander)            :: xextract

! PRIME2D PROGRAMS
type(prime2D_init_commander)       :: xprime2D_init
type(prime2D_commander)            :: xprime2D
type(cavgassemble_commander)       :: xcavgassemble
type(check2D_conv_commander)       :: xcheck2D_conv
type(rank_cavgs_commander)         :: xrank_cavgs

! PRIME3D PROGRAMS
type(resrange_commander)           :: xresrange
type(npeaks_commander)             :: xnpeaks
type(nspace_commander)             :: xnspace
type(shellweight3D_commander)      :: xshellweight3D
type(prime3D_init_commander)       :: xprime3D_init
type(het_init_commander)           :: xhet_init
type(multiptcl_init_commander)     :: xmultiptcl_init
type(prime3D_commander)            :: xprime3D
type(cont3D_commander)             :: xcont3D
type(check3D_conv_commander)       :: xcheck3D_conv
    
! COMMON-LINES PROGRAMS
type(comlin_smat_commander)        :: xcomlin_smat
type(symsrch_commander)            :: xsymsrch
    
! MASK PROGRAMS
type(mask_commander)               :: xmask
type(automask2D_commander)         :: xautomask2D
type(automask3D_commander)         :: xautomask3D
    
! RECONSTRUCTION PROGRAMS
type(eo_volassemble_commander)     :: xeo_volassemble
type(recvol_commander)             :: xrecvol
type(volassemble_commander)        :: xvolassemble

! CHECKER PROGRAMS    
type(check_box_commander)          :: xcheck_box
type(check_nptcls_commander)       :: xcheck_nptcls
type(iminfo_commander)             :: ximinfo

! VOLOPS PROGRAMS
type(cenvol_commander)             :: xcenvol
type(postproc_vol_commander)       :: xpostproc_vol
type(projvol_commander)            :: xprojvol
type(volaverager_commander)        :: xvolaverager
type(volops_commander)             :: xvolops
type(volume_smat_commander)        :: xvolume_smat

! GENERAL IMAGE PROCESSING PROGRAMS
type(binarise_commander)           :: xbinarise
type(convert_commander)            :: xconvert
type(corrcompare_commander)        :: xcorrcompare
type(ctfops_commander)             :: xctfops
type(filter_commander)             :: xfilter
type(image_smat_commander)         :: ximage_smat
type(norm_commander)               :: xnorm
type(respimg_commander)            :: xrespimg
type(scale_commander)              :: xscale
type(stack_commander)              :: xstack
type(stackops_commander)           :: xstackops
type(tseries_extract_commander)    :: xtseries_extract
    
! MISCELLANOUS PROGRAMS
type(cluster_smat_commander)       :: xcluster_smat
type(find_nnimgs_commander)        :: xfind_nnimgs
type(masscen_commander)            :: xmasscen
type(print_cmd_dict_commander)     :: xprint_cmd_dict
type(print_dose_weights_commander) :: xprint_dose_weights
type(print_fsc_commander)          :: xprint_fsc
type(res_commander)                :: xres
type(shift_commander)              :: xshift

! ORIENTATION DATA MANAGEMENT PROGRAMS
type(cluster_oris_commander)       :: xcluster_oris
type(makedeftab_commander)         :: xmakedeftab
type(makeoris_commander)           :: xmakeoris
type(map2ptcls_commander)          :: xmap2ptcls
type(orisops_commander)            :: xorisops
type(oristats_commander)           :: xoristats
type(rotmats2oris_commander)       :: xrotmats2oris
    
! PARALLEL PROCESSING PROGRAMS
type(merge_algndocs_commander)     :: xmerge_algndocs
type(merge_nnmat_commander)        :: xmerge_nnmat
type(merge_shellweights_commander) :: xmerge_shellweights
type(merge_similarities_commander) :: xmerge_similarities  
type(split_pairs_commander)        :: xsplit_pairs
type(split_commander)              :: xsplit
    
! OTHER DECLARATIONS
integer, parameter    :: MAXNKEYS=100, KEYLEN=32
character(len=KEYLEN) :: keys_required(MAXNKEYS)='', keys_optional(MAXNKEYS)=''
character(len=STDLEN) :: arg, prg, entire_line
type(cmdline)         :: cline
integer               :: cmdstat, cmdlen, pos
logical               :: describe
call get_command_argument(1, arg, cmdlen, cmdstat)
call get_command(entire_line)
if( str_has_substr(entire_line, 'prg=list') ) call list_all_simple_programs
describe = str_has_substr(entire_line, 'describe=yes')
pos = index(arg, '=') ! position of '='
call cmdline_err( cmdstat, cmdlen, arg, pos )
prg = arg(pos+1:)     ! this is the program name
if( str_has_substr(prg, 'simple') ) stop 'giving program names with simple* prefix is depreciated'
select case(prg)

   ! SIMULATOR PROGRAMS

    case( 'noiseimgs' )
        !==Program noiseimgs
        !
        ! <noiseimgs/begin>is a program for generating pure noise images<noiseimgs/end>
        !
        ! set required keys
        keys_required(1) = 'box'
        keys_required(2) = 'nptcls'
        ! parse command line
        if( describe ) call print_doc_noiseimgs
        call cline%parse(keys_required(:2))
        ! execute
        call xnoiseimgs%execute(cline)
    case( 'simimgs' )
        !==Program simimgs
        !
        ! <simimgs/begin>is a program for simulating cryo-EM images. It is not a very sophisticated simulator, but it is 
        ! nevertheless useful for testing purposes. It does not do any multi-slice simulation and it cannot be used for
        ! simulating molecules containing heavy atoms. It does not even accept a PDB file as an input. Input is a cryo-EM 
        ! map, which we usually generate from a PDB file using EMANs program pdb2mrc. simimgs then projects the 
        ! volume using Fourier interpolation, applies 20% of the total noise to the images (pink noise), Fourier transforms 
        ! them, and multiplies them with astigmatic CTF and B-factor. The images are inverse FTed before the remaining 80% 
        ! of the noise (white noise) is added<simimgs/end>
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
        if( describe ) call print_doc_simimgs
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
        call xsimimgs%execute(cline)
    case( 'simmovie' )
        !==Program simmovie
        !
        ! <simmovie/begin>is a program for crude simulation of a DDD movie. Input is a set of projection images to place.
        ! Movie frames are then generated related by randomly shifting the base image and applying two different noise
        ! sources: shot and detector noise<simmovie/end>
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
        if( describe ) call print_doc_simmovie
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
        call xsimmovie%execute(cline)
    case( 'simsubtomo' )
        !==Program simsubtomo
        !
        ! <simsubtomo/begin>is a program for crude simulation of a subtomograms<simsubtomo/end>
        !
        ! set required keys
        keys_required(1) = 'vol1'
        keys_required(2) = 'smpd'
        keys_required(3) = 'nptcls'
        keys_required(4) = 'snr'
        ! set optional keys
        keys_optional(1) = 'nthr'
        ! parse command line
        if( describe ) call print_doc_simsubtomo
        call cline%parse(keys_required(:4), keys_optional(:1))
        ! execute
        call xsimsubtomo%execute(cline)

    ! PRE-PROCESSING PROGRAMS

    case( 'preproc' )
        !==Program preproc
        !
        ! <preproc/begin>is a program that executes unblur, ctffind & pick in sequence
        ! <preproc/end>
        !
        ! set required keys
        keys_required(1)   = 'filetab'
        keys_required(2)   = 'smpd'
        keys_required(3)   = 'kv'
        keys_required(4)   = 'cs'
        keys_required(5)   = 'fraca'
        keys_required(6)   = 'refs'
        ! set optional keys
        keys_optional(1)   = 'nthr'
        keys_optional(2)   = 'fbody'
        keys_optional(3)   = 'lpstart'
        keys_optional(4)   = 'lpstop'
        keys_optional(5)   = 'trs'
        keys_optional(6)   = 'exp_time'
        keys_optional(7)   = 'dose_rate'
        keys_optional(8)   = 'pspecsz_unblur'
        keys_optional(9)   = 'pspecsz_ctffind'
        keys_optional(10)  = 'numlen'
        keys_optional(11)  = 'startit'
        keys_optional(12)  = 'scale'
        keys_optional(13)  = 'frameavg'
        keys_optional(14)  = 'tomo'
        keys_optional(15)  = 'hp_ctffind'
        keys_optional(16)  = 'lp_ctffind'
        keys_optional(17)  = 'lp_pick'
        keys_optional(18)  = 'dfmin'
        keys_optional(19)  = 'dfmax'
        keys_optional(20)  = 'astigstep'
        keys_optional(21)  = 'expastig'
        keys_optional(22)  = 'phaseplate'
        keys_optional(23)  = 'thres'
        ! parse command line
        if( describe ) call print_doc_preproc
        call cline%parse(keys_required(:6), keys_optional(:23))
        ! set defaults
        if( .not. cline%defined('trs')             ) call cline%set('trs',        5.)
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',   15.)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',     8.)
        if( .not. cline%defined('pspecsz_unblur')  ) call cline%set('pspecsz',  512.)
        if( .not. cline%defined('pspecsz_ctffind') ) call cline%set('pspecsz', 1024.)
        if( .not. cline%defined('hp_ctffind')      ) call cline%set('hp',        30.)
        if( .not. cline%defined('lp_ctffind')      ) call cline%set('lp',         5.)
        if( .not. cline%defined('lp_pick')         ) call cline%set('lp_pick',   20.)
        call xpreproc%execute(cline)
    case( 'select_frames' )
        !==Program select_frames
        !
        ! <select_frames/begin>is a program for selecting contiguous segments of frames from DDD movies
        ! <select_frames/end>
        !
        ! set required keys
        keys_required(1) = 'filetab'
        keys_required(2) = 'fbody'
        keys_required(3) = 'fromf'
        keys_required(4) = 'tof'
        keys_required(5) = 'smpd'
        ! set optional keys
        keys_optional(1) = 'startit'
        ! parse command line
        if( describe ) call print_doc_select_frames
        call cline%parse(keys_required(:5), keys_optional(:1))
        ! execute
        call xselect_frames%execute(cline)
    case( 'boxconvs' )
        !==Program boxconvs
        !
        ! <boxconvs/begin>is a program for averaging overlapping boxes across a micrograph
        ! in order to check if gain correction was appropriately done<boxconvs/end>
        !
        ! set required keys
        keys_required(1) = 'fbody'
        ! set optional keys
        keys_optional(1) = 'stk'
        keys_optional(2) = 'filetab'
        keys_optional(3) = 'boxconvsz'
        keys_optional(4) = 'startit'
        ! parse command line
        if( describe ) call print_doc_boxconvs
        call cline%parse(keys_required(:1), keys_optional(:4))
        ! set defaults
        if( .not. cline%defined('boxconvsz') ) call cline%set('boxconvsz', 512.)
        ! execute
        call xboxconvs%execute(cline)
    case( 'integrate_movies' )
        !==Program integrate_movies
        !
        ! <integrate_movies/begin>is a program for integrating DDD movies<integrate_movies/end>
        !
        ! set required keys
        keys_required(1)  = 'filetab'
        keys_required(2)  = 'fbody'
        keys_required(3)  = 'smpd'
        ! set optional keys
        keys_optional(1)  = 'pspecsz'
        keys_optional(2)  = 'scale'
        ! parse command line
        if( describe ) call print_doc_integrate_movies
        call cline%parse(keys_required(:3), keys_optional(:2))
        ! execute
        call xintegrate_movies%execute(cline)
    case( 'powerspecs' )
        !==Program powerspecs
        !
        ! <powerspecs/begin>is a program for generating powerspectra from a stack or filetable<powerspecs/end>
        !
        ! set required keys
        keys_required(1)  = 'smpd'
        keys_required(2)  = 'fbody'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'stk'
        keys_optional(3)  = 'filetab'
        keys_optional(4)  = 'pspecsz'
        keys_optional(5)  = 'speckind'
        keys_optional(6)  = 'startit'
        keys_optional(7)  = 'lp'
        keys_optional(8)  = 'clip'
        ! parse command line
        if( describe ) call print_doc_powerspecs
        call cline%parse(keys_required(:2), keys_optional(:8))
        ! set defaults
        if( .not. cline%defined('pspecsz') ) call cline%set('pspecsz', 512.)
        if( .not. cline%defined('lp')      ) call cline%set('lp',        6.)
        if( .not. cline%defined('clip')    ) call cline%set('clip',    256.)
        ! execute
        call xpowerspecs%execute(cline)
    case( 'unblur' )
        !==Program unblur
        !
        ! <unblur/begin>is a program for movie alignment or unblurring.
        ! Input is a textfile with absolute paths to movie files in addition to a few obvious input
        ! parameters<unblur/end>
        !
        ! set required keys
        keys_required(1)  = 'filetab'
        keys_required(2)  = 'smpd'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'fbody'
        keys_optional(3)  = 'lpstart'
        keys_optional(4)  = 'lpstop'
        keys_optional(5)  = 'trs'
        keys_optional(6)  = 'exp_time'
        keys_optional(7)  = 'dose_rate'
        keys_optional(8)  = 'kv'
        keys_optional(9)  = 'pspecsz'
        keys_optional(10) = 'numlen'
        keys_optional(11) = 'startit'
        keys_optional(12) = 'scale'
        keys_optional(13) = 'frameavg'
        keys_optional(14) = 'tomo'
        ! parse command line
        if( describe ) call print_doc_unblur
        call cline%parse(keys_required(:2), keys_optional(:14))
        ! set defaults
        if( .not. cline%defined('trs')     ) call cline%set('trs',      5.)
        if( .not. cline%defined('lpstart') ) call cline%set('lpstart', 15.)
        if( .not. cline%defined('lpstop')  ) call cline%set('lpstop',   8.)
        ! execute
        call xunblur%execute(cline)
    case( 'ctffind' )
        !==Program ctffind
        !
        ! <ctffind/begin>is a wrapper program for CTFFIND4 (Grigorieff lab)<ctffind/end> 
        !
        ! set required keys
        keys_required(1) = 'filetab'
        keys_required(2) = 'smpd'
        keys_required(3) = 'kv'
        keys_required(4) = 'cs'
        keys_required(5) = 'fraca'
        ! set optional keys
        keys_optional(1) = 'pspecsz'
        keys_optional(2) = 'hp'
        keys_optional(3) = 'lp'
        keys_optional(4) = 'dfmin'
        keys_optional(5) = 'dfmax'
        keys_optional(6) = 'astigstep'
        keys_optional(7) = 'expastig'
        keys_optional(8) = 'phaseplate'
        ! parse command line
        if( describe ) call print_doc_ctffind
        call cline%parse(keys_required(:5), keys_optional(:8))
        ! set defaults
        if( .not. cline%defined('pspecsz') ) call cline%set('pspecsz', 1024.)
        if( .not. cline%defined('hp')      ) call cline%set('hp',        30.)
        if( .not. cline%defined('lp')      ) call cline%set('lp',         5.)
        ! execute
        call xctffind%execute(cline)
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
    case( 'makepickrefs' )
        !==Program pickrefs
        !
        ! <pickrefs/begin>is a program for generating references for template-based particle picking<pickrefs/end> 
        !
        ! set required keys
        keys_required(1) = 'pgrp'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'vol1'
        keys_optional(3) = 'stk'
        keys_optional(4) = 'smpd'
        ! parse command line
        ! if( describe ) call print_doc_makepickrefs
        call cline%parse(keys_required(:1), keys_optional(:4))
        ! execute
        call xmakepickrefs%execute(cline)
    case( 'pick' )
        !==Program pick
        !
        ! <pick/begin>is a template-based picker program<pick/end> 
        !
        ! set required keys
        keys_required(1) = 'filetab'
        keys_required(2) = 'refs'
        keys_required(3) = 'smpd'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'lp'
        keys_optional(3) = 'thres'
        ! parse command line
        if( describe ) call print_doc_pick
        call cline%parse(keys_required(:3), keys_optional(:3))
        ! execute
        call xpick%execute(cline)
    case( 'extract' )
        !==Program extract
        !
        ! <extract/begin>is a program that extracts particle images from DDD movies
        ! or integrated movies. Boxfiles are assumed to be in EMAN format but we provide
        ! a conversion script (relion2emanbox.pl) for *.star files containing
        ! particle coordinates obtained with Relion. The program creates one stack per movie
        ! frame as well as a stack of corrected framesums. In addition to single-particle
        ! image stacks, the program produces a parameter file extract_params.txt
        ! that can be used in conjunction with other SIMPLE programs. We obtain CTF parameters 
        ! with CTFFIND4<extract/end> 
        !
        ! set required keys
        keys_required(1) = 'filetab'
        keys_required(2) = 'boxtab'
        keys_required(3) = 'smpd'
        ! set optional keys
        keys_optional(1) = 'msk'
        keys_optional(2) = 'ctffind_doc'
        keys_optional(3) = 'neg'
        keys_optional(4) = 'box'
        keys_optional(5) = 'noise_norm'
        keys_optional(6) = 'outside'
        ! parse command line
        if( describe ) call print_doc_extract
        call cline%parse( keys_required(:3), keys_optional(:6))
        ! parse command line
        if( .not. cline%defined('neg') )call cline%set('neg', 'yes')
        ! execute
        call xextract%execute(cline)
        
    ! PRIME2D PROGRAMS

    case( 'prime2D_init' )
        !==Program prime2D_init
        !
        ! <prime2D_init/begin>is used  to produce the initial random references for prime2D execution.
        ! The random clustering and in-plane alignment is printed in the file prime2D_startdoc.txt 
        ! produced by the program. This file is used together with the initial references
        ! (startcavgs.ext) to execute prime2D<prime2D_init/end> 
        !
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'smpd'
        keys_required(3) = 'ncls'
        keys_required(4) = 'ctf'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'deftab'
        keys_optional(3) = 'oritab'
        keys_optional(4) = 'filwidth'
        keys_optional(5) = 'mul'
        keys_optional(6) = 'srch_inpl'
        keys_optional(7) = 'tseries'
        ! parse command line
        if( describe ) call print_doc_prime2D_init
        call cline%parse(keys_required(:4), keys_optional(:7))
        ! set defaults
        if( .not. cline%defined('eo') ) call cline%set('eo', 'no')
        ! execute
        call xprime2D_init%execute(cline)
    case( 'prime2D' )
        !==Program prime2D
        !
        ! <prime2D/begin>is a reference-free 2D alignment/clustering algorithm adopted from the prime3D 
        ! probabilistic ab initio 3D reconstruction algorithm<prime2D/end>
        !
        ! set required keys
        keys_required(1)  = 'stk'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'msk'
        keys_required(4)  = 'ncls'
        keys_required(5)  = 'refs'
        keys_required(6)  = 'ctf'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'deftab'
        keys_optional(3)  = 'oritab'
        keys_optional(4)  = 'hp'
        keys_optional(5)  = 'lp'
        keys_optional(6)  = 'cenlp'
        keys_optional(7)  = 'trs'
        keys_optional(8)  = 'automsk'
        keys_optional(9)  = 'amsklp'
        keys_optional(10) = 'inner'
        keys_optional(11) = 'width'
        keys_optional(12) = 'startit'
        keys_optional(13) = 'maxits'
        keys_optional(14) = 'srch_inpl'
        ! parse command line
        if( describe ) call print_doc_prime2D
        call cline%parse(keys_required(:6), keys_optional(:14))
        ! set defaults
        if( .not. cline%defined('lp')     ) call cline%set('lp',     20.)
        if( .not. cline%defined('amsklp') ) call cline%set('amsklp', 25.)
        if( .not. cline%defined('cenlp')  ) call cline%set('cenlp',  30.)
        if( .not. cline%defined('edge')   ) call cline%set('edge',   20.)
        if( .not. cline%defined('eo')     ) call cline%set('eo',    'no')
        ! execute
        call xprime2D%execute(cline)
    case( 'cavgassemble' )
        !==Program cavgassemble
        !
        ! <cavgassemble/begin>is a program that assembles class averages when the clustering 
        ! program (prime2D) has been executed in distributed mode<cavgassemble/end>
        !
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'smpd'
        keys_required(3) = 'ncls'
        keys_required(4) = 'oritab'
        keys_required(5) = 'nparts'
        keys_required(6) = 'ctf'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'deftab'
        keys_optional(3) = 'inner'
        keys_optional(4) = 'width'
        keys_optional(5) = 'which_iter'        
        ! parse command line
        if( describe ) call print_doc_cavgassemble
        call cline%parse(keys_required(:6), keys_optional(:5))
        ! execute
        call xcavgassemble%execute(cline)
    case( 'check2D_conv' )
        !==Program check2D_conv
        !
        ! <check2D_conv/begin>is a program for checking if a PRIME2D run has converged.
        ! The statistics outputted include (1) the overlap between the distribution of parameters
        ! for succesive runs. (2) The percentage of search space scanned, i.e. how many reference
        ! images are evaluated on average. (3) The average correlation between the images and 
        ! their corresponding best matching reference section. If convergence to a local optimum
        ! is achieved, the fraction increases. Convergence is achieved if the parameter distribution
        ! overlap is larger than 0.95 and more than 99% of the reference sections need to be 
        ! searched to find an improving solution<check2D_conv/end>
        !
        ! set required keys
        keys_required(1) = 'smpd'
        keys_required(2) = 'box'
        keys_required(3) = 'oritab'
        keys_required(4) = 'nptcls'
        ! set optional keys
        keys_optional(1) = 'lp'
        ! parse command line
        if( describe ) call print_doc_check2D_conv
        call cline%parse(keys_required(:4), keys_optional(:1))
        ! set defaults
        if( .not. cline%defined('lp') ) call cline%set('lp', 20.)
        ! execute
        call xcheck2D_conv%execute(cline)
    case( 'rank_cavgs' )
        !==Program rank_cavgs
        !
        ! <rank_cavgs/begin>is a program for ranking class averages by decreasing population, given the
        ! stack of class averages (stk argument) and the 2D orientations document (oritab)
        ! generated by prime2D<rank_cavgs/end>
        !
        ! set required keys
        keys_required(1)  = 'stk'
        keys_required(2)  = 'oritab'
        ! set optional keys
        keys_optional(1)  = 'outstk'
        ! parse command line
        if( describe ) call print_doc_rank_cavgs
        call cline%parse(keys_required(:2), keys_optional(:1))
        ! execute
        call xrank_cavgs%execute(cline)
        
    ! PRIME3D PROGRAMS
    
    case( 'resrange' )
        !==Program resrange
        !
        ! <resrange/begin>is a program for estimating the resolution range used in the heuristic
        ! resolution-stepping scheme in the PRIME3D initial model production procedure. The initial
        ! low-pass limit is set so that each image receives ten nonzero orientation weights. When 
        ! quasi-convergence has been reached, the limit is updated one Fourier index at the time,
        ! until PRIME reaches the condition where six nonzero orientation weights are assigned to
        ! each image. FSC-based filtering is unfortunately not possible to do in the ab initio
        ! 3D reconstruction step, because when the orientations are mostly random, the FSC overestimates
        ! the resolution. This program is used internally when executing PRIME in distributed mode. We 
        ! advise you to check the starting and stopping low-pass limits before executing PRIME3D using 
        ! this program. The resolution range estimate depends on the molecular diameter,
        ! which is estimated based on the box size. If you want to override this estimate, set moldiam 
        ! to the desired value (in A). This may be necessary if your images have a lot of background
        ! padding. However, for starting model generation it is probably better to clip the images snugly
        ! around the particle, because smaller images equal less computation<resrange/end>
        !
        ! set required keys
        keys_required(1) = 'smpd'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'nspace'
        keys_optional(3) = 'pgrp'
        keys_optional(4) = 'box'
        keys_optional(5) = 'moldiam'
        ! parse command line
        if( describe ) call print_doc_resrange
        call cline%parse(keys_required(:1), keys_optional(:5))
        ! set defaults
        if( .not. cline%defined('nspace') ) call cline%set('nspace', 1000.)
        ! execute
        call xresrange%execute(cline)
    case( 'npeaks' )
        !==Program npeaks
        !
        ! <npeaks/begin>is a program for checking the number of nonzero orientation weights (number of correlation peaks 
        ! included in the weighted reconstruction)<npeaks/end>
        !
        ! set required keys
        keys_required(1)  = 'smpd'
        keys_required(2)  = 'box'
        keys_required(3)  = 'lp'
        ! set optional keys
        keys_optional(1)  = 'nspace'
        keys_optional(2)  = 'moldiam'
        keys_optional(3)  = 'pgrp'
        ! parse command line
        if( describe ) call print_doc_npeaks
        call cline%parse(keys_required(:3), keys_optional(:3))
        ! set defaults
        if( .not. cline%defined('nspace') ) call cline%set('nspace', 1000.)
        if( .not. cline%defined('lp') )     call cline%set('lp',       20.)
        ! execute
        call xnpeaks%execute(cline)
    case( 'nspace' )
        !==Program nspace
        !
        ! <nspace/begin>is a program for calculating the expected resolution obtainable with different values of nspace
        ! (number of discrete projection directions used for discrete search)<nspace/end>
        !
        ! set required keys
        keys_required(1)  = 'moldiam'
        ! parse command line
        if( describe ) call print_doc_nspace
        call cline%parse(keys_required=keys_required(:1))
        ! execute
        call xnspace%execute(cline)
    case( 'shellweight3D' )
        !==Program shellweight3D
        !
        ! <shellweight3D/begin>is a program for calculating the shell-by-shell resolution weights in a global sense, so that 
        ! particles that do contribute with higher resolution information (as measured by the FRC) are given the appropriate 
        ! weight<shellweight3D/end>
        !
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'vol1'
        keys_required(3) = 'smpd'
        keys_required(4) = 'msk' 
        keys_required(5) = 'oritab'
        keys_required(6) = 'ctf'
        ! set optional CTF-related keys
        keys_optional(1) = 'deftab'
        keys_optional(2) = 'automsk'
        keys_optional(3) = 'mw'
        keys_optional(4) = 'amsklp'
        keys_optional(5) = 'edge'
        keys_optional(6) = 'binwidth'
        keys_optional(7) = 'inner'
        keys_optional(8) = 'width'        
        ! parse command line
        if( describe ) call print_doc_shellweight3D
        call cline%parse(keys_required(:6), keys_optional(:8))
        ! execute
        call xshellweight3D%execute(cline)
    case( 'prime3D_init' )
        !==Program prime3D_init
        !
        ! <prime3D_init/begin>is a program for generating a random initial model for initialisation of PRIME3D.
        ! If the data set is large (>5000 images), generating a random model can be slow. To speedup, set 
        ! nran to some smaller number, resulting in nran images selected randomly for 
        ! reconstruction<prime3D_init/end> 
        !
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'smpd'
        keys_required(3) = 'msk'
        keys_required(4) = 'ctf'
        keys_required(5) = 'pgrp'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'deftab'
        keys_optional(3) = 'lp'
        keys_optional(4) = 'inner'
        keys_optional(5) = 'width'
        keys_optional(6) = 'nspace'
        keys_optional(7) = 'nran'
        keys_optional(8) = 'npeaks'
        keys_optional(9) = 'xfel'      
        ! parse command line
        if( describe ) call print_doc_prime3D_init
        call cline%parse(keys_required(:5), keys_optional(:9))
        ! set defaults
        if( .not. cline%defined('eo')     ) call cline%set('eo',      'no')
        if( .not. cline%defined('nspace') ) call cline%set('nspace', 1000.)
        ! execute
        call xprime3D_init%execute(cline)
    case( 'multiptcl_init' )
        !==Program multiptcl_init
        !
        ! <multiptcl_init/begin>is a program for generating random initial models for initialisation of PRIME3D
        ! when run in multiparticle mode<multiptcl_init/end> 
        !
        ! set required keys
        keys_required(1)  = 'stk'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'oritab'
        keys_required(4)  = 'ctf'
        keys_required(5)  = 'nstates'
        keys_required(6)  = 'msk'
        ! set optionnal keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'deftab'
        keys_optional(3)  = 'inner'
        keys_optional(4)  = 'width'
        keys_optional(5)  = 'lp'
        keys_optional(6)  = 'eo'
        keys_optional(7)  = 'frac'
        keys_optional(8)  = 'state2split'
        keys_optional(9)  = 'norec'
        keys_optional(10) = 'mul'
        keys_optional(11) = 'zero'
        keys_optional(12) = 'tseries'
        ! parse command line
        if( describe ) call print_doc_multiptcl_init
        call cline%parse(keys_required(:6), keys_optional(:12))
        ! set defaults
        if( .not. cline%defined('trs') ) call cline%set('trs', 3.) ! to assure that shifts are being used
        !execute
        call xmultiptcl_init%execute(cline)     
    case( 'het_init' )
        !==Program het_init
        !
        ! <het_init/begin> <het_init/end> 
        !
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'smpd'
        keys_required(3) = 'oritab'
        keys_required(4) = 'nstates'
        keys_required(5) = 'ctf'
        ! set optionnal keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'deftab'
        keys_optional(3) = 'msk'
        keys_optional(4) = 'inner'
        keys_optional(5) = 'width'
        keys_optional(6) = 'lp'
        keys_optional(7) = 'eo'
        keys_optional(8) = 'frac'        
        ! parse command line
        call cline%parse(keys_required(:5), keys_optional(:8))
        ! set defaults
        if( .not. cline%defined('eo') ) call cline%set('eo', 'no')
        !execute
        call xhet_init%execute(cline)    
    case( 'prime3D' )
        !==Program prime3D
        !
        ! <prime3D/begin>is an ab inito reconstruction/refinement program based on probabilistic
        ! projection matching. PRIME is short for PRobabilistic Initial 3D Model generation for Single-
        ! particle cryo-Electron microscopy. Do not search the origin shifts initially, when the model is 
        ! of very low quality. If your images are far off centre, use stackops with option
        ! shalgn=yes instead to shiftalign the images beforehand (the algorithm implemented is the 
        ! same as EMANs cenalignint program). We recommend running the first round of PRIME with 
        ! the default dynamic resolution stepping dynlp=yes. The dynlp option implements 
        ! a heuristic resolution weighting/update scheme. The initial low-pass limit is set so that each
        ! image receives ten nonzero orientation weights. When quasi-convergence has been reached, the limit 
        ! is updated one Fourier index at the time until PRIME reaches the condition where six nonzero 
        ! orientation weights are assigned to each image. FSC-based filtering is unfortunately not possible
        ! to do in the ab initio reconstruction step, because when the orientations are mostly random, the 
        ! FSC overestimates the resolution. Once the initial model has converged, we recommend start searching 
        ! the shifts (by setting trs to some nonzero value) and applying the FSC for resolution-
        ! weighting (by setting eo=yes)<prime3D/end>
        !
        ! set required keys
        keys_required(1)  = 'stk'
        keys_required(2)  = 'vol1'
        keys_required(3)  = 'smpd'
        keys_required(4)  = 'msk'
        keys_required(5)  = 'ctf'
        keys_required(6)  = 'pgrp'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'vol2'
        keys_optional(3)  = 'oritab'
        keys_optional(4)  = 'deftab'
        keys_optional(5)  = 'trs'
        keys_optional(6)  = 'hp'
        keys_optional(7)  = 'lp'
        keys_optional(8)  = 'cenlp'
        keys_optional(9)  = 'dynlp'
        keys_optional(10) = 'lpstart'
        keys_optional(11) = 'lpstop'
        keys_optional(12) = 'eo'
        keys_optional(13) = 'refine'
        keys_optional(14) = 'frac'
        keys_optional(15) = 'automsk'
        keys_optional(16) = 'mw'
        keys_optional(17) = 'amsklp'
        keys_optional(18) = 'edge'
        keys_optional(19) = 'binwidth'
        keys_optional(20) = 'inner'
        keys_optional(21) = 'width'
        keys_optional(22) = 'nspace'
        keys_optional(23) = 'nstates'
        keys_optional(24) = 'npeaks'
        keys_optional(25) = 'startit'
        keys_optional(26) = 'maxits'
        keys_optional(27) = 'shbarrier'
        keys_optional(28) = 'noise'
        keys_optional(29) = 'xfel'
        keys_optional(30) = 'nnn'
        keys_optional(31) = 'shellw'
        ! parse command line
        if( describe ) call print_doc_prime3D
        call cline%parse(keys_required(:6), keys_optional(:31))
        ! set defaults
        if( .not. cline%defined('nspace')                  ) call cline%set('nspace', 1000.)
        if( cline%defined('lp') .or. cline%defined('find') ) call cline%set('dynlp',   'no')
        if( .not. cline%defined('cenlp')  )                  call cline%set('cenlp',    30.)
        if( .not. cline%defined('refine')                  ) call cline%set('refine',  'no')
        ! execute
        call xprime3D%execute(cline)
    case( 'cont3D' )
        !==Program cont3D
        !
        ! <cont3D/begin>is a continuous refinement code under development<cont3D/end>
        !
        ! set required keys
        keys_required(1)  = 'stk'
        keys_required(2)  = 'vol1'
        keys_required(3)  = 'smpd'
        keys_required(4)  = 'msk'
        keys_required(5)  = 'oritab'
        keys_required(6)  = 'trs'
        keys_required(7)  = 'ctf'
        keys_required(8)  = 'pgrp'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'deftab'
        keys_optional(3)  = 'frac'
        keys_optional(4)  = 'automsk'
        keys_optional(5)  = 'mw'
        keys_optional(6)  = 'amsklp'
        keys_optional(7)  = 'edge'
        keys_optional(8)  = 'inner'
        keys_optional(9)  = 'width'
        keys_optional(10) = 'hp'
        keys_optional(11) = 'lpstart'
        keys_optional(12) = 'lpstop'
        keys_optional(13) = 'startit'
        keys_optional(14) = 'maxits'
        keys_optional(15) = 'xfel'
        keys_optional(16) = 'shellw'
        ! parse command line
        if( describe ) call print_doc_cont3D
        call cline%parse(keys_required(:8), keys_optional(:16))
        ! set defaults
        call cline%set('eo',     'yes')
        !call cline%set('refine', 'yes')
        call cline%set('dynlp',  'no')
        if( .not.cline%defined('shellw') )call cline%set('shellw','no')
        ! execute
        call xcont3D%execute(cline)        
    case( 'check3D_conv' )
        !==Program check3D_conv
        !
        ! <check3D_conv/begin>is a program for checking if a PRIME3D run has converged. The statistics
        ! outputted include (1) angle of feasible region, which is proportional to the angular
        ! resolution of the set of discrete projection directions being searched. (2) The average angular
        ! distance between orientations in the present and previous iteration. In the early iterations,
        ! the distance is large because a diverse set of orientations is explored. If convergence to a
        ! local optimum is achieved, the distance decreases. (3) The percentage of search space scanned,
        ! i.e. how many reference images are evaluated on average. (4) The average correlation between
        ! the images and their corresponding best matching reference sections. (5) The average standard
        ! deviation of the Euler angles. Convergence is achieved if the angular distance between the
        ! orientations in successive iterations falls significantly below the angular resolution of the
        ! search space and more than 99% of the reference sections need to be matched on average 
        ! <check3D_conv/end>
        !
        ! set required keys
        keys_required(1) = 'smpd'
        keys_required(2) = 'box'
        keys_required(3) = 'oritab'
        keys_required(4) = 'nptcls'
        keys_required(5) = 'pgrp'
        ! set optional keys
        keys_optional(1) = 'lp'
        keys_optional(2) = 'nstates'
        keys_optional(3) = 'eo'
        keys_optional(4) = 'nspace'
        keys_optional(5) = 'find'
        keys_optional(6) = 'refine'
        ! parse command line
        if( describe ) call print_doc_check3D_conv
        call cline%parse(keys_required(:5), keys_optional(:6))
        ! set defaults
        if( .not. cline%defined('lp') .and. .not.cline%defined('lp') )call cline%set('lp', 20.)
        if( .not. cline%defined('nspace') ) call cline%set('nspace', 1000.)
        ! execute
        call xcheck3D_conv%execute(cline)
    
    ! COMMON-LINES PROGRAMS
    
    case( 'comlin_smat' )
        !==Program comlin_smat
        !
        ! <comlin_smat/begin>is a program for creating a similarity matrix based on common
        ! line correlation. The idea being that it should be possible to cluster images based
        ! on their 3D similarity witout having a 3D model by only operating on class averages
        ! and find averages that fit well together in 3D<comlin_smat/end>
        !
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'smpd'
        keys_required(3) = 'lp'
        keys_required(4) = 'msk'
        ! set optional keys
        keys_optional(1) = 'hp'
        keys_optional(2) = 'trs'
        ! parse command line
        if( describe ) call print_doc_comlin_smat
        call cline%parse(keys_required(:4), keys_optional(:2))
        ! execute
        call xcomlin_smat%execute(cline)
    case( 'symsrch' )
        !==Program symsrch
        !
        ! <symsrch/begin>is a program for searching for the principal symmetry axis of a volume 
        ! reconstructed without assuming any point-group symmetry or assessing the degree of symmetry
        ! of class averages or individual particles of higher pointgroups (dn,t,o,i). The program takes 
        ! as input an asymmetrical reconstruction or stack of class averages/individual particles.
        ! For volumes, the alignment document for all the particle images that have gone into the 3D 
        ! reconstruction and the desired point-group symmetry needs to be inputted. The 3D reconstruction 
        ! is then projected in 150 (default option) even directions, common lines-based optimisation is 
        ! used to identify the principal symmetry axis, the rotational transformation is applied to the 
        ! inputted orientations, and a new alignment document is produced. Input this document to 
        ! recvol together with the images and the point-group symmetry to generate a symmetrised map. 
        ! If you are unsure about the point-group, you should use the compare=yes mode and input the highest 
        ! conceviable point-group. The program then calculates probabilities for all lower groups inclusive.
        ! The class average/particle option operates in an equivalent fashion but with individual images. 
        ! The output is then a per-image correlation value that informs about how well the image conforms 
        ! to to inputted point-group. The state parameter allows you to apply symmetry for the given state.
        ! <symsrch/end>
        !        
        ! set required keys
        keys_required(1) = 'smpd'
        keys_required(2) = 'msk'
        keys_required(3) = 'pgrp'
        keys_required(4) = 'outfile'
        keys_required(5) = 'lp'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'vol1'
        keys_optional(3) = 'stk'
        keys_optional(4) = 'oritab'
        keys_optional(5) = 'cenlp'
        keys_optional(6) = 'hp'
        keys_optional(7) = 'nspace'
        keys_optional(8) = 'compare'
        ! parse command line
        if( describe ) call print_doc_symsrch
        call cline%parse(keys_required(:5), keys_optional(:8))
        ! set defaults
        if( .not. cline%defined('nspace') )then
            if( cline%defined('vol1') )then
                call cline%set('nptcls', 150.) ! 50 projections 4 symsrch
                call cline%set('nspace', 150.) ! 50 projections 4 symsrch
            else
                call cline%set('nptcls', 1.) ! single-particle search
                call cline%set('nspace', 1.) ! single-particle search
            endif
        else
            call cline%set('nptcls', cline%get_rarg('nspace'))
        endif
        if( .not. cline%defined('cenlp')   ) call cline%set('cenlp',    50.)
        if( .not. cline%defined('compare') ) call cline%set('compare', 'no')
        ! execute
        call xsymsrch%execute(cline)
        
    ! MASK PROGRAMS

    case( 'mask' )
       !==Program mask
       !
       ! <mask/begin>is a program for masking images and volumes.
       ! If you want to mask your images with a spherical mask with a soft falloff, set msk
       ! to the radius in pixels<mask/end> 
       !
       ! set required keys
       keys_required(1)  = 'msk'
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
       keys_optional(10) = 'binwidth'
       keys_optional(11) = 'amsklp'
       keys_optional(12) = 'automsk'
       keys_optional(13) = 'smpd'
       ! parse command line
       if( describe ) call print_doc_mask
       call cline%parse( keys_required(:1), keys_optional(:13))
       ! execute
       call xmask%execute(cline)
    case( 'automask2D' )
        !==Program automask2D
        !
        ! <automask2D/begin>is a program for solvent flattening of class averages. 
        ! The algorithm for background removal is based on low-pass filtering and binarization. 
        ! First, the class averages are low-pass filtered to amsklp. Binary representatives 
        ! are then generated by assigning foreground pixels using sortmeans. A cosine function
        ! softens the edge of the binary mask before it is  multiplied with the unmasked input 
        ! averages to accomplish flattening<automask2D/end>
        !
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'smpd'
        keys_required(3) = 'msk'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'amsklp'
        keys_optional(3) = 'edge'
        ! parse command line
        if( describe ) call print_doc_automask2D
        call cline%parse(keys_required(:3), keys_optional(:3))
        ! set defaults
        call cline%set('automsk', 'yes')
        if( .not. cline%defined('amsklp') ) call cline%set('amsklp', 25.)
        if( .not. cline%defined('edge')   ) call cline%set('edge', 20.)
        ! execute
        call xautomask2D%execute(cline)
    case( 'automask3D' )
        !==Program automask3D
        !
        ! <automask3D/begin>is a program for solvent flattening of a volume. The algorithm 
        ! for background removal is based on low-pass filtering and binarization. First, the volume is low-pass
        ! filtered to amsklp. A binary volume is then generated by assigning foreground pixels 
        ! (=1) based on the volume calculated from the molecular weight. A cosine function softens the edge 
        ! of the binary mask before it is  multiplied with the unmasked input to generate the flattened map
        ! <automask3D/end>
        !
        ! set required keys
        keys_required(1) = 'vol1'
        keys_required(2) = 'smpd'
        keys_required(3) = 'msk'
        keys_required(4) = 'mw'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'vol2'
        keys_optional(3) = 'amsklp'
        keys_optional(4) = 'edge'
        keys_optional(5) = 'binwidth'
        ! parse command line
        if( describe ) call print_doc_automask3D
        call cline%parse(keys_required(:4), keys_optional(:5))
        ! set defaults
        call cline%set('automsk', 'yes')
        ! execute
        call xautomask3D%execute(cline)

    ! RECONSTRUCTION PROGRAMS

    case( 'recvol' )
        !==Program recvol
        !
        ! <recvol/begin>is a program for reconstructing volumes from MRC and SPIDER stacks, given input 
        ! orientations and state assignments. The algorithm is based on direct Fourier inversion with a 
        ! Kaiser-Bessel (KB) interpolation kernel. This window function reduces the real-space ripple 
        ! artifacts associated with direct moving windowed-sinc interpolation. The feature sought when 
        ! implementing this algorithm was to enable quick, reliable reconstruction from aligned individual 
        ! particle images. mul is used to scale the origin shifts if down-sampled 
        ! were used for alignment and the original images are used for reconstruction. ctf=yes or ctf=flip 
        ! turns on the Wiener restoration. If the images were phase-flipped set ctf=flip. amsklp, mw, and edge 
        ! control the solvent mask: the low-pass limit used to generate the envelope; the molecular weight of the 
        ! molecule (protein assumed but it works reasonably well also for RNA; slight modification of mw 
        ! might be needed). The inner parameter controls the radius of the soft-edged mask used to remove 
        ! the unordered DNA/RNA core of spherical icosahedral viruses<recvol/end>
        !
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'smpd'
        keys_required(3) = 'oritab'
        keys_required(4) = 'msk'
        keys_required(5) = 'ctf'
        keys_required(6) = 'pgrp'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'eo'
        keys_optional(3) = 'deftab'
        keys_optional(4) = 'frac'
        keys_optional(5) = 'mw'
        keys_optional(6) = 'mul'
        keys_optional(7) = 'state'
        ! parse command line
        if( describe ) call print_doc_recvol
        call cline%parse(keys_required(:6), keys_optional(:7))
        ! set defaults
        if( .not. cline%defined('trs') ) call cline%set('trs', 5.) ! to assure that shifts are being used
        if( .not. cline%defined('eo')  ) call cline%set('eo', 'no')
        ! execute
        call xrecvol%execute(cline)
    case( 'eo_volassemble' )
        !==Program eo_volassemble
        !
        ! <eo_volassemble/begin>is a program that assembles volume(s) when the reconstruction
        ! program (recvol with eo=yes) has been executed in distributed mode. inner applies a soft-edged 
        ! inner mask. An inner mask is used for icosahedral virus reconstruction, because the
        ! DNA or RNA core is often unordered and  if not removed it may negatively impact the
        ! alignment. The width parameter controls the fall-off of the edge of the 
        ! inner mask<eo_volassemble/end>
        !
        ! set required keys
        keys_required(1)  = 'stk'
        keys_required(2)  = 'nparts'
        keys_required(3)  = 'smpd'
        keys_required(4)  = 'msk'
        keys_required(5)  = 'oritab'
        keys_required(6)  = 'ctf'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'deftab'
        keys_optional(3)  = 'mul'
        keys_optional(4)  = 'state'
        keys_optional(5)  = 'nstates'
        ! parse command line
        if( describe ) call print_doc_eo_volassemble
        call cline%parse(keys_required(:6), keys_optional(:5))
        ! set defaults
        if( cline%defined('state') ) call cline%set('nstates', 1.) ! to assure that shifts are being used
        ! execute
        call xeo_volassemble%execute(cline)
    case( 'volassemble' )
        !==Program volassemble
        !
        ! <volassemble/begin>is a program that assembles volume(s) when the reconstruction program 
        ! (recvol) has been executed in distributed mode. odd is used to assemble the odd reconstruction, 
        ! even is used to assemble the even reconstruction, eo is used to assemble both the even and the 
        ! odd reconstruction and state is used to assemble the inputted state. Normally, you do not fiddle with 
        ! these parameters. They are used internally<volassemble/end>
        ! 
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'nparts'
        keys_required(3) = 'smpd'
        keys_required(4) = 'oritab'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'even'
        keys_optional(3) = 'odd'
        keys_optional(4) = 'eo'
        keys_optional(5) = 'state'
        keys_optional(6) = 'xfel'
        ! parse command line
        if( describe ) call print_doc_volassemble
        call cline%parse(keys_required(:4), keys_optional(:6))
        ! set defaults
        if( cline%defined('state')    ) call cline%set('nstates', 1.)
        if( .not. cline%defined('eo') ) call cline%set('eo', 'no')
        ! execute
        call xvolassemble%execute(cline)

    ! CHECKER PROGRAMS
        
    case( 'check_box' )
        !==Program check_box
        !
        ! <check_box/begin>is a program for checking the image dimensions of MRC and SPIDER
        !  stacks and volumes<check_box/end>
        !
        ! set optional keys
        keys_optional(1) = 'stk'
        keys_optional(2) = 'vol1'
        ! parse command line
        if( describe ) call print_doc_check_box
        call cline%parse( keys_optional=keys_optional(:2))
        ! execute
        call xcheck_box%execute(cline)
    case( 'check_nptcls' )
        !==Program check_nptcls
        !
        ! <check_nptcls/begin>is a program for checking the number of images in MRC and SPIDER
        ! stacks<check_nptcls/end>
        !
        ! set optional keys
        keys_required(1) = 'stk'
        ! parse command line
        if( describe ) call print_doc_check_nptcls
        call cline%parse(keys_required(:1))
        ! execute
        call xcheck_nptcls%execute(cline) 
    case( 'iminfo' )
        !==Program iminfo
        !
        ! <iminfo/begin>is a program for printing header information in MRC and SPIDER stacks 
        ! and volumes<iminfo/end>
        !
        ! set required keys
        keys_required(1) = 'fname'
        ! set optional keys
        keys_optional(1) = 'box'
        keys_optional(2) = 'smpd'
        keys_optional(3) = 'stats'
        keys_optional(4) = 'endian'
        ! parse command line
        if( describe ) call print_doc_iminfo
        call cline%parse(keys_required(:1), keys_optional(:4))
        ! execute
        call ximinfo%execute(cline)
        
    ! VOLOPS PROGRAMS
    
    case( 'cenvol' )
        !==Program cenvol
        !
        ! <cenvol/begin>is a program for centering a volume and mapping the shift parameters
        ! back to the particle images<cenvol/end>
        !
        ! set required keys
        keys_required(1) = 'vol1'
        keys_required(2) = 'smpd'
        ! set optional keys
        keys_optional(1) = 'oritab'
        keys_optional(2) = 'outfile'
        keys_optional(3) = 'cenlp'
        ! parse command line
        if( describe ) call print_doc_cenvol
        call cline%parse(keys_required(:2), keys_optional(:3))
        ! set defaults
        if( .not. cline%defined('cenlp') ) call cline%set('cenlp', 50.)
        ! execute
        call xcenvol%execute(cline)
    case( 'postproc_vol' )
        !==Program postproc_vol
        !
        ! <postproc_vol/begin>is a program for post-processing of volumes<postproc_vol/end>
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
        ! parse command line
        if( describe ) call print_doc_postproc_vol
        call cline%parse(keys_required(:3), keys_optional(:8))
        ! execute
        call xpostproc_vol%execute(cline)
    case( 'projvol' )
        !==Program projvol
        !
        ! <projvol/begin>is a program for projecting a volume using interpolation in Fourier space. Input is a SPIDER or
        ! MRC volume. Output is a stack of projection images of the same format as the inputted volume. Projections
        ! are generated by extraction of central sections from the Fourier volume and back transformation of the 2D FTs.
        ! nspace controls the number of projection images generated with quasi-even projection directions. The
        ! oritab parameter allows you to input the orientations that you wish to have your volume projected in. If
        ! rnd=yes, random rather than quasi-even projections are generated, trs then controls the halfwidth of
        ! the random origin shift. Less commonly used parameters are pgrp, which controls the point-group symmetry
        ! c (rotational), d (dihedral), t (tetrahedral), o (octahedral) or i (icosahedral). The point-group symmetry is 
        ! used to restrict the set of projections to within the asymmetric unit. 
        ! neg inverts the contrast of the projections. mirr=yes mirrors the projection by modifying the Euler 
        ! angles. If mirr=x or mirr=y the projection is physically mirrored after it has been generated<projvol/end>
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
        keys_optional(9)  = 'mirr'
        keys_optional(10) = 'top'
        keys_optional(11) = 'xfel'
        ! parse command line
        if( describe ) call print_doc_projvol
        call cline%parse(keys_required(:2), keys_optional(:11))
        ! set defaults
        if( .not. cline%defined('wfun')  ) call cline%set('wfun', 'kb')
        if( .not. cline%defined('winsz') ) call cline%set('winsz', 1.5)
        if( .not. cline%defined('alpha') ) call cline%set('alpha', 2.)
        ! execute
        call xprojvol%execute(cline)
    case( 'volaverager' )
        !==Program volaverager
        !
        ! <volaverager/begin>is a program for averaging volumes according to state label in oritab
        ! <volaverager/end> 
        !
        ! set required keys
        keys_required(1) = 'vollist'
        keys_required(2) = 'oritab'
        ! set optional keys
        keys_optional(1) = 'nthr'
        ! parse command line
        if( describe ) call print_doc_volaverager
        call cline%parse(keys_required(:2), keys_optional(:1))
        ! execute
        call xvolaverager%execute(cline)
    case( 'volops' )
        !==Program volops
        !
        ! <volops/begin>provides standard single-particle image processing routines that are applied to MRC or SPIDER volumes
        ! <volops/end>
        !
        ! set required keys
        keys_required(1) = 'vol1'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'guinier'
        keys_optional(3) = 'smpd'
        keys_optional(4) = 'hp'
        keys_optional(5) = 'lp'
        keys_optional(6) = 'neg'
        keys_optional(7) = 'snr'
        keys_optional(8) = 'mirr'
        keys_optional(9) = 'outvol'
        ! parse command line
        if( describe ) call print_doc_volops
        call cline%parse(keys_required(:1), keys_optional(:9))
        ! execute
        call xvolops%execute(cline)
    case( 'volume_smat' )
        !==Program volume_smat
        !
        ! <volume_smat/begin>is a program for creating a similarity matrix based on volume2volume
        ! correlation<olume_smat/end>
        !
        ! set required keys
        keys_required(1) = 'vollist'
        keys_required(2) = 'smpd'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'hp'
        keys_optional(3) = 'lp'
        keys_optional(4) = 'msk'
        ! parse command line
        if( describe ) call print_doc_volume_smat
        call cline%parse(keys_required(:2), keys_optional(:4))
        ! execute
        call xvolume_smat%execute(cline)
        
    ! GENERAL IMAGE PROCESSING PROGRAMS
    
    case( 'binarise' )
        !==Program binarise
        !
        ! <binarise/begin>is a program for binarisation of stacks and volumes<binarise/end> 
        !
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'stk'
        keys_optional(3)  = 'vol1'
        keys_optional(4)  = 'thres'
        keys_optional(5)  = 'npix'
        keys_optional(6)  = 'grow'
        keys_optional(7)  = 'edge'
        keys_optional(8)  = 'neg'
        keys_optional(9)  = 'outvol'
        keys_optional(10) = 'outstk'
        ! parse command line
        if( describe ) call print_doc_binarise
        call cline%parse(keys_optional=keys_optional(:10))
        ! execute
        call xbinarise%execute(cline)
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
    case( 'corrcompare' )
        !==Program corrcompare
        !
        ! <corrcompare/begin>is a program for comparing stacked images using real-space and Fourier-based approaches
        ! <corrcompare/end> 
        !
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'stk2'
        ! set optional keys
        keys_optional(1) = 'msk'
        keys_optional(2) = 'stats'
        keys_optional(3) = 'lp'
        keys_optional(4) = 'smpd'
        ! parse command line
        if( describe ) call print_doc_corrcompare
        call cline%parse(keys_required(:2), keys_optional(:4))
        ! execute
        call xcorrcompare%execute(cline)
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
        keys_required(1) = 'smpd'
        ! Optional keys
        keys_optional(1) = 'stk'
        keys_optional(2) = 'vol1'
        keys_optional(3) = 'outstk'
        keys_optional(4) = 'outvol'
        keys_optional(5) = 'lp'
        keys_optional(6) = 'hp'
        keys_optional(7) = 'phrand'
        keys_optional(8) = 'bfac'
        ! parse command line
        if( describe ) call print_doc_filter
        call cline%parse(keys_required(:1), keys_optional(:8))
        ! execute
        call xfilter%execute(cline)
    case( 'image_smat' )
        !==Program image_smat
        !
        ! <image_smat/begin>is a program for creating a similarity matrix based on common line correlation. The idea
        ! being that it should be possible to cluster images based on their 3D similarity witout having a 3D model
        ! by only operating on class averages and find averages that fit well together in 3D<image_smat/end>
        !
        ! set required keys
        keys_required(1)  = 'stk'
        keys_required(2)  = 'smpd'
        ! set optional keys
        keys_optional(1)  = 'lp'
        keys_optional(2)  = 'msk'
        keys_optional(3)  = 'hp'
        keys_optional(4)  = 'nthr'
        ! parse command line
        if( describe ) call print_doc_image_smat
        call cline%parse(keys_required(:2), keys_optional(:4))
        ! execute
        call ximage_smat%execute(cline)
    case( 'norm' )
        !==Program norm
        !
        ! <norm/begin>is a program for normalization of MRC or SPIDER stacks and volumes. If you want to
        ! normalise your images inputted with stk, set norm=yes. hfun (e.g. hfun=sigm) controls
        ! the normalisation function. If you want to perform noise normalisation of the images set 
        ! noise_norm=yes given a mask radius msk (pixels). If you want to normalise your
        ! images or volume (vol1) with respect to their power spectrum set shell_norm=yes
        ! <norm/end> 
        !
        ! set optional keys
        keys_optional(1)  = 'stk'
        keys_optional(2)  = 'msk'
        keys_optional(3)  = 'norm'
        keys_optional(4)  = 'noise_norm'
        keys_optional(5)  = 'shell_norm'
        keys_optional(6)  = 'hfun'
        keys_optional(7)  = 'nthr'
        ! parse command line
        if( describe ) call print_doc_norm
        call cline%parse(keys_optional=keys_optional(:7))
        ! execute
        call xnorm%execute(cline)
    case( 'respimg' )
        !==Program respimg
        !
        ! <respimg/begin> 
        ! <respimg/end> 
        !
        ! set required keys
        keys_required(1)  = 'stk'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'msk'
        keys_required(4)  = 'oritab'
        keys_required(5)  = 'ctf'
        ! set optional keys
        keys_optional(1)  = 'outstk'
        keys_optional(2)  = 'deftab'
        keys_optional(3)  = 'nthr'
        ! parse command line
        if( describe ) call print_doc_respimg
        call cline%parse( keys_required(:5), keys_optional(:3))
        ! execute
        call xrespimg%execute(cline)
    case( 'scale' )
        !==Program scale
        !
        ! <scale/begin>is a program that provides re-scaling and clipping routines for MRC or SPIDER stacks
        ! and volumes<scale/end> 
        !
        keys_required(1) = 'smpd'
        ! set optional keys
        keys_optional(1) = 'stk'
        keys_optional(2) = 'vol1'
        keys_optional(3) = 'filetab'
        keys_optional(4) = 'msk'
        keys_optional(5) = 'newbox'
        keys_optional(6) = 'scale'
        keys_optional(7) = 'clip'
        keys_optional(8) = 'outvol'
        keys_optional(9) = 'outstk'
        ! parse command line
        if( describe ) call print_doc_scale
        call cline%parse(keys_required(:1),keys_optional(:9))
        ! set defaults
        if( .not. cline%defined('scale') ) call cline%set('scale', 1.)
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
        ! set optional keys
        keys_optional(1) = 'lp'
        keys_optional(2) = 'smpd'
        keys_optional(3) = 'clip'
        keys_optional(4) = 'nframes'
        keys_optional(5) = 'fbody'
        keys_optional(6) = 'numlen'
        keys_optional(7) = 'xdim'
        keys_optional(8) = 'ydim'
        keys_optional(9) = 'endian'
        ! parse command line
        if( describe ) call print_doc_stack
        call cline%parse(keys_required(:2), keys_optional(:9))
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
        ! If you define frameavg to some integer number larger than one averages with chunk sizes of frameavg are produced,
        ! which may be useful for analysis of dose-fractionated image series. neg inverts the contrast of the images<stackops/end>
        !
        ! Required keys
        keys_required(1) = 'stk'
        ! set optional keys
        keys_optional(1)  = 'smpd'
        keys_optional(2)  = 'oritab'
        keys_optional(3)  = 'outstk'
        keys_optional(4)  = 'mirr'
        keys_optional(5)  = 'nran'
        keys_optional(6)  = 'frac'
        keys_optional(7)  = 'state'
        keys_optional(8)  = 'class'
        keys_optional(9)  = 'neg'
        keys_optional(10) = 'acf'
        keys_optional(11) = 'avg'
        keys_optional(12) = 'frameavg'
        keys_optional(13) = 'vis'
        keys_optional(14) = 'snr'
        keys_optional(15) = 'fromp'
        keys_optional(16) = 'top'
        keys_optional(17) = 'stk2'
        keys_optional(18) = 'nptcls'
        keys_optional(19) = 'append'
        keys_optional(20) = 'order'
        keys_optional(21) = 'outfile'
        ! parse command line
        if( describe ) call print_doc_stackops
        call cline%parse( keys_required(:1),keys_optional(:21) )
        ! execute
        call xstackops%execute(cline)
    case( 'tseries_extract' )
        !==Program tseries_extract
        !
        ! <tseries_extract/begin>is a program for creating overlapping chunks of frameavg frames from time-series data
        ! <tseries_extract/end> 
        !
        ! set required keys
        keys_required(1) = 'filetab'
        keys_required(2) = 'smpd'
        keys_required(3) = 'frameavg'
        ! set optional keys
        keys_optional(1)  = 'fbody'
        ! parse command line
        ! if( describe ) call print_doc_tseries_extract
        call cline%parse(keys_required(:3), keys_optional(:1))
        ! execute
        call xtseries_extract%execute(cline)
        
    ! MISCELLANOUS PROGRAMS
        
    case( 'cluster_smat' )
        !==Program cluster_smat
        !
        ! <cluster_smat/begin>is a program for clustering a similarity matrix and use
        ! an combined cluster validation index to assess the quality of the clustering
        ! based on the number of clusters<cluster_smat/end>
        !
        ! set required keys
        keys_required(1) = 'nptcls'
        keys_required(2) = 'fname'
        keys_required(3) = 'ncls'
        keys_required(4) = 'label'
        ! set optional keys
        keys_optional(1) = 'nthr'
        ! parse command line
        if( describe ) call print_doc_cluster_smat
        call cline%parse(keys_required(:4), keys_optional(:1))
        ! execute
        call xcluster_smat%execute(cline)
    case( 'find_nnimgs' )
        !==Program find_nnimgs
        !
        ! <find_nnimgs/begin>is a program for cidentifying the nnn nearest neighbor
        ! images for each image in the inputted stack<find_nnimgs/end>
        !
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'smpd'
        keys_required(3) = 'msk'
        keys_required(4) = 'nnn'
        ! set optional keys
        keys_optional(1) = 'lp'
        keys_optional(2) = 'hp'
        ! parse command line
        if( describe ) call print_doc_find_nnimgs
        call cline%parse(keys_required(:4), keys_optional(:2))
        ! execute
        call xfind_nnimgs%execute(cline)
    case( 'masscen' )
        !==Program masscen
        !
        ! <masscen/begin>is a program for centering images acccording to their
        ! centre of mass<masscen/end>
        !
        ! set required keys
        keys_required(1) = 'stk' 
        keys_required(2) = 'smpd'
        keys_required(3) = 'lp'
        ! set optional keys
        keys_optional(1) = 'msk'
        keys_optional(2) = 'neg'
        keys_optional(3) = 'thres'
        keys_optional(4) = 'outstk'
        ! parse command line
        if( describe ) call print_doc_masscen
        call cline%parse(keys_required(:3), keys_optional(:4))
        ! execute
        call xmasscen%execute(cline)
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
    case( 'print_dose_weights' )
        !==Program print_dose_weights
        !
        ! <print_dose_weights/begin>is a program for printing the dose weights applied to individual frames<print_dose_weights/end>
        !
        ! set required keys
        keys_required(1) = 'nframes' 
        keys_required(2) = 'exp_time'
        keys_required(3) = 'dose_rate'
        keys_required(4) = 'box'
        keys_required(5) = 'smpd'
        ! set optional keys
        keys_optional(1) = 'kv'
        ! parse command line
        if( describe ) call print_doc_print_dose_weights
        call cline%parse(keys_required(:5),keys_optional(:1))
        ! execute
        call xprint_dose_weights%execute(cline)
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
    case( 'res' )
        !==Program res
        !
        ! <res/begin>is a program for checking the low-pass resolution limit for a given Fourier index<res/end>
        !
        !set required keys
        keys_required(1)  = 'smpd'
        keys_required(2)  = 'find'
        keys_required(3)  = 'box'
        ! parse command line
        if( describe ) call print_doc_res
        call cline%parse(keys_required(:3))
        !execute
        call xres%execute(cline)
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
        ! parse command line
        if( describe ) call print_doc_shift
        call cline%parse(keys_required(:3), keys_optional(:2))
        ! execute
        call xshift%execute(cline)

    ! ORIENTATION DATA MANAGEMENT PROGRAMS

    case( 'cluster_oris' )
        !==Program cluster_oris
        !
        ! <cluster_oris/begin>is a program for clustering orientations based on geodesic distance<cluster_oris/end>
        !
        ! Required keys
        keys_required(1)  = 'oritab'
        keys_required(2)  = 'ncls'
        ! set optional keys
        keys_optional(1) = 'nthr'
        ! parse command line
        if( describe ) call print_doc_cluster_oris
        call cline%parse(keys_required(:2), keys_optional(:1))
        ! execute
        call xcluster_oris%execute(cline)
    case( 'makedeftab' )
        !==Program makedeftab
        !
        ! <makedeftab/begin>is a program for creating a SIMPLE conformant file of CTF parameter values (deftab).
        ! Input is either an earlier SIMPLE deftab/oritab. The purpose is to get the kv, cs, and fraca parameters
        ! as part of the CTF input doc as that is the new convention. The other alternative is to input a plain text
        ! file with CTF parameters dfx, dfy, angast according to the Frealign convention. Unit conversions are dealt 
        ! with using optional variables. The units refer to the units in the inputted document<makedeftab/end>
        !
        ! Required keys
        keys_required(1) = 'kv'
        keys_required(2) = 'cs'
        keys_required(3) = 'fraca'
        keys_required(4) = 'outfile'
        ! set optional keys
        keys_optional(1) = 'plaintexttab'
        keys_optional(2) = 'oritab'
        keys_optional(3) = 'deftab'
        keys_optional(4) = 'dfunit'
        keys_optional(5) = 'angastunit'
        ! parse command line
        if( describe ) call print_doc_makedeftab
        call cline%parse(keys_required(:4),keys_optional(:5))
        ! execute
        call xmakedeftab%execute(cline)
    case( 'makeoris' )
        !==Program makeoris
        !
        ! <makeoris/begin>is a program for making SIMPLE orientation/parameter files (text files containing input parameters and/or 
        ! parameters estimated by prime2D or prime3D). The program generates random
        ! Euler angles e1.in.[0,360], e2.in.[0,180], and e3.in.[0,360] and random origin 
        ! shifts x.in.[-trs,yrs] and y.in.[-trs,yrs]. If ndiscrete is set to an integer number > 0, the 
        ! orientations produced are randomly sampled from the set of ndiscrete quasi-even projection directions, and the in-plane 
        ! parameters are assigned randomly. If even=yes, then all nptcls orientations are assigned 
        ! quasi-even projection directions,and random in-plane parameters. If nstates is set to some integer number > 0, then 
        ! states are assigned randomly .in.[1,nstates]. If zero=yes in this mode of execution, the projection 
        ! directions are zeroed and only the in-plane parameters are kept intact. If errify=yes and astigerr is defined, 
        ! then uniform random astigmatism errors are introduced .in.[-astigerr,astigerr]<makeoris/end>
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
        keys_optional(8)  = 'defocus'
        keys_optional(9)  = 'angerr'
        keys_optional(10) = 'sherr'
        keys_optional(11) = 'dferr'
        keys_optional(12) = 'even'
        keys_optional(13) = 'zero'
        keys_optional(14) = 'discrete'
        keys_optional(15) = 'ndiscrete'
        keys_optional(16) = 'diverse'
        keys_optional(17) = 'state'
        keys_optional(18) = 'nspace'
        keys_optional(19) = 'iares'
        ! parse command line
        if( describe ) call print_doc_makeoris
        call cline%parse(keys_required(:1),keys_optional(:19))
        ! execute
        call xmakeoris%execute(cline)
    case( 'map2ptcls' )
        !==Program map2ptcls
       !
       ! <map2ptcls/begin>is a program for mapping parameters that have been obtained using class averages to 
       ! the individual particle images<map2ptcls/end> 
       !
       ! set required keys
       keys_required(1)  = 'stk'
       keys_required(2)  = 'stk2'
       keys_required(3)  = 'stk3'
       keys_required(4)  = 'oritab'
       ! set optional keys
       keys_optional(1)  = 'oritab2'
       keys_optional(2)  = 'comlindoc'
       keys_optional(3)  = 'doclist'
       keys_optional(4)  = 'deftab'
       keys_optional(5)  = 'outfile'
       keys_optional(6)  = 'mul'
       keys_optional(7)  = 'nthr'
       ! parse command line
       if( describe ) call print_doc_map2ptcls
       call cline%parse(keys_required(:4), keys_optional(:7))
       ! set defaults
       if( .not. cline%defined('outfile') ) call cline%set('outfile', 'mapped_ptcls_params.txt')
       ! execute
       call xmap2ptcls%execute(cline)
    case( 'orisops' )
        !==Program orisops
        !
        ! <orisops/begin>is a program for analyzing SIMPLE orientation/parameter files (text files containing input parameters 
        ! and/or parameters estimated by prime2D or prime3D). If only oritab is inputted, there 
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
        keys_optional(3)  = 'plot'
        keys_optional(4)  = 'outfile'
        keys_optional(5)  = 'e1'
        keys_optional(6)  = 'e2'
        keys_optional(7)  = 'e3'
        keys_optional(8)  = 'trs'
        keys_optional(9)  = 'nstates'
        keys_optional(10)  = 'pgrp'
        keys_optional(11) = 'defocus'
        keys_optional(12) = 'deftab'
        keys_optional(13) = 'angerr'
        keys_optional(14) = 'sherr'
        keys_optional(15) = 'dferr'
        keys_optional(16) = 'zero'
        keys_optional(17) = 'discrete'
        keys_optional(18) = 'ndiscrete'
        keys_optional(19) = 'state'
        keys_optional(20) = 'errify'
        keys_optional(21) = 'mul'
        keys_optional(22) = 'mirr'
        keys_optional(23) = 'xsh'
        keys_optional(24) = 'ysh'
        keys_optional(25) = 'zsh'
        ! parse command line
        if( describe ) call print_doc_orisops
        call cline%parse(keys_optional=keys_optional(:25))
        ! execute
        call xorisops%execute(cline)
    case( 'oristats' )
        !==Program oristats
        !
        ! <oristats/begin>is a program for analyzing SIMPLE orientation/parameter files (text files containing input 
        ! parameters and/or parameters estimated by prime2D or prime3D). If two orientation 
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
        keys_optional(8)  = 'hist'
        keys_optional(9)  = 'ncls'
        keys_optional(10) = 'minp'
        keys_optional(11) = 'clustvalid'
        keys_optional(12) = 'thres'
        ! parse command line
        if( describe ) call print_doc_oristats
        call cline%parse( keys_required(:1), keys_optional(:12) )
        ! execute
        call xoristats%execute(cline)
    case( 'rotmats2oris' )
        !==Program rotmats2oris
        !
        ! <rotmats2oris/begin>converts a text file (9 records per line) describing 
        ! rotation matrices into a SIMPLE oritab<rotmats2oris/end>
        !
        ! Required keys
        keys_required(1)  = 'infile'
        ! set optional keys
        keys_optional(1)  = 'outfile'
        ! parse command line
        if( describe ) call print_doc_rotmats2oris
        call cline%parse( keys_required(:1), keys_optional(:1) )
        ! execute
        call xrotmats2oris%execute(cline)
        
    ! PARALLEL PROCESSING PROGRAMS
    
    case( 'merge_algndocs' )
        !==Program merge_algndocs
        !
        ! <merge_algndocs/begin>is a program for merging alignment documents from SIMPLE
        ! runs in distributed mode<merge_algndocs/end>
        !
        ! set required keys
        keys_required(1) = 'fbody'
        keys_required(2) = 'nptcls'
        keys_required(3) = 'ndocs'
        keys_required(4) = 'outfile'
        ! parse command line
        if( describe ) call print_doc_merge_algndocs
        call cline%parse(keys_required(:4))
        ! execute
        call xmerge_algndocs%execute(cline)
    case( 'merge_nnmat' )
        !==Program merge_nnmat
        !
        ! <merge_nnmat/begin>is a program for merging partial nearest neighbour matrices calculated 
        ! in distributed mode<merge_nnmat/end>
        !
        ! set required keys
        keys_required(1) = 'nptcls'
        keys_required(2) = 'nparts'
        keys_required(3) = 'nnn'
        ! parse command line
        if( describe ) call print_doc_merge_nnmat
        call cline%parse( keys_required(:3) )
        ! execute
        call xmerge_nnmat%execute(cline)
    case( 'merge_shellweights' )
        !==Program merge_shellweights
        !
        ! <merge_shellweights/begin>is a program for merging partial shellweight matrices calculated 
        ! in distributed mode<merge_shellweights/end>
        !
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'nparts'
        ! parse command line
        if( describe ) call print_doc_merge_shellweights
        call cline%parse(keys_required(:2))
        ! set defaults
        call cline%set('outfile', 'shellweight3D_doc.txt')
        ! execute
        call xmerge_shellweights%execute(cline)
    case( 'merge_similarities' )
        !==Program merge_similarities
        !
        ! <merge_similarities/begin>is a program for merging similarities calculated between pairs of objects 
        ! into a similarity matrix that can be inputted to cluster_smat<merge_similarities/end>
        !
        ! set required keys
        keys_required(1)  = 'nptcls'
        ! set optional keys
        keys_optional(1)  = 'nparts'
        ! parse command line
        if( describe ) call print_doc_merge_similarities
        call cline%parse(keys_required(:1), keys_optional(:1))
        ! execute
        call xmerge_similarities%execute(cline)    
    case( 'split_pairs' )
        !==Program split_pairs
        !
        ! <split_pairs/begin>is a program for splitting calculations between pairs of objects 
        ! into balanced partitions<split_pairs/end>
        !
        ! set required keys
        keys_required(1) = 'nptcls'
        keys_required(2) = 'nparts'
        ! parse command line
        if( describe ) call print_doc_split_pairs
        call cline%parse(keys_required(:2))
        ! execute
        call xsplit_pairs%execute(cline)
    case( 'split' )
        !==Program split
        !
        ! <split/begin>is a program for splitting of image stacks into partitions for parallel execution.
        ! This is done to reduce I/O latency<split/end>
        !
        ! set required keys
        keys_required(1) = 'stk'
        ! set optional keys
        keys_optional(1) = 'split_mode'
        keys_optional(2) = 'nparts'
        keys_optional(3) = 'neg'
        keys_optional(4) = 'ncls'
        keys_optional(5) = 'nspace'
        ! parse command line
        if( describe ) call print_doc_split
        call cline%parse(keys_required(:1), keys_optional(:5))
        ! execute
        call xsplit%execute(cline)
    case DEFAULT
        write(*,'(a,a)') 'program key (prg) is: ', trim(prg)
        stop 'unsupported program'
    end select
end program simple_exec
