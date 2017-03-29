!==Program simple_distr_exec
!
! <simple_distr_exec/begin> executes the parallel (or distributed workflows) of SIMPLE <simple_distr_exec/end>
!
! The code is distributed with the hope that it will be useful, but WITHOUT ANY WARRANTY.
! Redistribution and modification is regulated by the GNU General Public License.
! Authors: Cyril Reboul & Hans Elmlund 2016
!
program simple_distr_exec
use simple_defs  
use simple_cmdline, only: cmdline
use simple_strings, only: str_has_substr
use simple_jiffys,  only: cmdline_err
use simple_gen_doc
use simple_restart
use simple_commander_stream_wflows
use simple_commander_distr_wflows
use simple_commander_hlev_wflows
implicit none

! DISTRIBUTED COMMANDERS
! pre-processing
type(preproc_stream_commander)            :: xpreproc_stream
type(unblur_distr_commander)              :: xunblur_distr
type(unblur_tomo_movies_distr_commander)  :: xunblur_tomo_distr
type(ctffind_distr_commander)             :: xctffind_distr
type(pick_distr_commander)                :: xpick_distr
! PRIME2D
type(prime2D_init_distr_commander)        :: xprime2D_init_distr
type(prime2D_distr_commander)             :: xprime2D_distr
type(find_nnimgs_distr_commander)         :: xfind_nnimgs_distr
! PRIME3D
type(prime3D_init_distr_commander)        :: xprime3D_init_distr
type(prime3D_distr_commander)             :: xprime3D_distr
type(cont3D_distr_commander)             ::  xcont3D_distr
type(shellweight3D_distr_commander)       :: xshellweight3D_distr
type(recvol_distr_commander)              :: xrecvol_distr
! time-series workflows
type(tseries_track_distr_commander)       :: xtseries_track_distr
! high-level workflows
type(iterated_spectral_weights_commander) :: xisw_distr
type(ini3D_from_cavgs_commander)          :: xini3D_from_cavgs
type(het_ensemble_commander)              :: xhet_ensemble

! OTHER DECLARATIONS
integer, parameter    :: MAXNKEYS=100, KEYLEN=32
character(len=KEYLEN) :: keys_required(MAXNKEYS)='', keys_optional(MAXNKEYS)=''
character(len=STDLEN) :: arg, prg, entire_line
type(cmdline)         :: cline
integer               :: cmdstat, cmdlen, pos
logical               :: describe, is_restart
call get_command_argument(1, arg, cmdlen, cmdstat)
call get_command(entire_line)
if( str_has_substr(entire_line, 'prg=list') ) call list_all_simple_distr_programs
describe = str_has_substr(entire_line, 'describe=yes')
pos = index(arg, '=') ! position of '='
call cmdline_err( cmdstat, cmdlen, arg, pos )
prg = arg(pos+1:) ! this is the program name
if( str_has_substr(prg, 'simple_') ) stop 'giving program names with simple_* prefix is depreciated'
select case(prg)

    ! PRE-PROCESSING STREAM, LINKING UNBLUR + CTFFIND + PICK

    case( 'preproc' )
        !==Program preproc
        !
        ! <preproc/begin>is a program that executes unblur, ctffind & pick in sequence
        ! and in streaming mode as the microscope collects the data <preproc/end>
        !
        ! set required keys
        keys_required(1)   = 'smpd'
        keys_required(2)   = 'kv'
        keys_required(3)   = 'cs'
        keys_required(4)   = 'fraca'
        keys_required(5)   = 'refs'
        keys_required(6)   = 'dir_movies'
        keys_required(7)   = 'dir_target'
        keys_required(8)   = 'ncunits'
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
        keys_optional(24)  = 'rm_outliers'
        ! parse command line
        if( describe ) call print_doc_preproc
        call cline%parse(keys_required(:8), keys_optional(:24))
        ! set defaults
        if( .not. cline%defined('trs')             ) call cline%set('trs',        5.)
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',   15.)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',     8.)
        if( .not. cline%defined('pspecsz_unblur')  ) call cline%set('pspecsz',  512.)
        if( .not. cline%defined('pspecsz_ctffind') ) call cline%set('pspecsz', 1024.)
        if( .not. cline%defined('hp_ctffind')      ) call cline%set('hp',        30.)
        if( .not. cline%defined('lp_ctffind')      ) call cline%set('lp',         5.)
        if( .not. cline%defined('lp_pick')         ) call cline%set('lp_pick',   20.)
        call xpreproc_stream%execute(cline)

    ! UNBLUR_MOVIES

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
        keys_required(3)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'ncunits'
        keys_optional(3)  = 'fbody'
        keys_optional(4)  = 'lpstart'
        keys_optional(5)  = 'lpstop'
        keys_optional(6)  = 'trs'
        keys_optional(7)  = 'exp_time'
        keys_optional(8)  = 'dose_rate'
        keys_optional(9)  = 'kv'
        keys_optional(10) = 'pspecsz'
        keys_optional(11) = 'numlen'
        keys_optional(12) = 'startit'
        keys_optional(13) = 'scale'
        keys_optional(14) = 'frameavg'
        keys_optional(15) = 'fromf'
        keys_optional(16) = 'tof'
        ! parse command line
        if( describe ) call print_doc_unblur
        call cline%parse(keys_required(:3), keys_optional(:16))
        ! set defaults
        if( .not. cline%defined('trs')     ) call cline%set('trs',      5.)
        if( .not. cline%defined('lpstart') ) call cline%set('lpstart', 15.)
        if( .not. cline%defined('lpstop')  ) call cline%set('lpstop',   8.)
        ! execute
        call xunblur_distr%execute(cline)
    case( 'unblur_tomo' )
        !==Program unblur
        !
        ! <unblur_tomo/begin>is a program for movie alignment or unblurring of tomographic movies.
        ! Input is a textfile with absolute paths to movie files in addition to a few obvious input
        ! parameters<unblur/end>
        !
        ! set required keys
        keys_required(1)  = 'tomoseries'
        keys_required(2)  = 'exp_doc'
        keys_required(3)  = 'smpd'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'ncunits'
        keys_optional(3)  = 'lpstart'
        keys_optional(4)  = 'lpstop'
        keys_optional(5)  = 'trs'
        keys_optional(6)  = 'kv'
        keys_optional(7)  = 'pspecsz'
        keys_optional(8)  = 'numlen'
        keys_optional(9)  = 'startit'
        keys_optional(10) = 'scale'
        keys_optional(11) = 'frameavg'
        ! parse command line
        if( describe ) call print_doc_unblur_tomo
        call cline%parse(keys_required(:3), keys_optional(:11))
        ! set defaults
        if( .not. cline%defined('trs')     ) call cline%set('trs',      5.)
        if( .not. cline%defined('lpstart') ) call cline%set('lpstart', 15.)
        if( .not. cline%defined('lpstop')  ) call cline%set('lpstop',   8.)
        if( .not. cline%defined('tomo')    ) call cline%set('tomo',  'yes')
        ! execute
        call xunblur_tomo_distr%execute(cline)

    ! CTFFIND

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
        keys_required(6) = 'nparts'
        ! set optional keys
        keys_optional(1) = 'ncunits'
        keys_optional(2) = 'pspecsz'
        keys_optional(3) = 'hp'
        keys_optional(4) = 'lp'
        keys_optional(5) = 'dfmin'
        keys_optional(6) = 'dfmax'
        keys_optional(7) = 'astigstep'
        keys_optional(8) = 'expastig'
        keys_optional(9) = 'phaseplate'
        ! parse command line
        if( describe ) call print_doc_ctffind
        call cline%parse(keys_required(:6), keys_optional(:9))
        ! set defaults
        call cline%set('nthr', 1.0)
        if( .not. cline%defined('pspecsz') ) call cline%set('pspecsz', 1024.)
        if( .not. cline%defined('hp')      ) call cline%set('hp',        30.)
        if( .not. cline%defined('lp')      ) call cline%set('lp',         5.)
        ! execute
        call xctffind_distr%execute(cline)

    ! PARTICLE PICKER

    case( 'pick' )
        !==Program pick
        !
        ! <pick/begin>is a template-based picker program<pick/end> 
        !
        ! set required keys
        keys_required(1) = 'filetab'
        keys_required(2) = 'refs'
        keys_required(3) = 'smpd'
        keys_required(4) = 'nparts'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'lp'
        keys_optional(3) = 'thres'
        keys_optional(4) = 'rm_outliers'
        ! parse command line
        if( describe ) call print_doc_pick
        call cline%parse(keys_required(:4), keys_optional(:4))
        ! execute
        call xpick_distr%execute(cline)

    ! PRIME2D

    case( 'prime2D_init' )
        !==Program prime2D_init
        !
        ! <prime2D_init/begin>is used to produce the initial random references for prime2D execution.
        ! The random clustering and in-plane alignment is printed in the file
        ! prime2D_startdoc.txt produced by the program. This file is used together with the initial references
        ! (startcavgs.ext) to execute prime2D<prime2D_init/end> 
        !
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'smpd'
        keys_required(3) = 'ncls'
        keys_required(4) = 'ctf'
        keys_required(5) = 'nparts'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'ncunits'
        keys_optional(3) = 'deftab'
        keys_optional(4) = 'oritab'
        keys_optional(5) = 'filwidth'
        keys_optional(6) = 'mul'   
        ! parse command line
        if( describe ) call print_doc_prime2D_init
        call cline%parse(keys_required(:5), keys_optional(:6))
        ! execute
        call xprime2D_init_distr%execute(cline)
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
        keys_required(5)  = 'ctf'
        keys_required(6)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'ncunits'
        keys_optional(3)  = 'deftab'
        keys_optional(4)  = 'refine'
        keys_optional(5)  = 'refs'
        keys_optional(6)  = 'oritab'
        keys_optional(7)  = 'hp'
        keys_optional(8)  = 'lp'
        keys_optional(9)  = 'cenlp'
        keys_optional(10)  = 'trs'
        keys_optional(11) = 'automsk'
        keys_optional(12) = 'amsklp'
        keys_optional(13) = 'inner'
        keys_optional(14) = 'width'
        keys_optional(15) = 'startit'
        keys_optional(16) = 'maxits'
        keys_optional(17) = 'filwidth'
        keys_optional(18) = 'nnn'
        keys_optional(19) = 'minp'
        keys_optional(20) = 'center'     
        ! documentation
        if( describe ) call print_doc_prime2D
        ! parse command line
        ! call check_restart( entire_line, is_restart )
        ! if( is_restart )then
        !     call parse_restart('prime2D', entire_line, cline, keys_required(:7), keys_optional(:19))
        ! else
        !     call cline%parse( keys_required(:7), keys_optional(:19) )
        ! endif
        call cline%parse( keys_required(:6), keys_optional(:19) )
        ! set defaults
        if( .not. cline%defined('lp')     ) call cline%set('lp',     20.)
        if( .not. cline%defined('eo')     ) call cline%set('eo',    'no')
        if( .not. cline%defined('amsklp') ) call cline%set('amsklp', 25.)
        if( .not. cline%defined('cenlp')  ) call cline%set('cenlp',  30.)
        if( .not. cline%defined('edge')   ) call cline%set('edge',   20.)
        if( .not. cline%defined('center') ) call cline%set('center', 'yes')
        ! execute
        call xprime2D_distr%execute(cline)
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
        keys_required(4) = 'nparts'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'ncunits'
        keys_optional(3) = 'nnn'
        keys_optional(4) = 'lp'
        keys_optional(5) = 'hp'
        ! parse command line
        if( describe ) call print_doc_find_nnimgs
        call cline%parse(keys_required(:4), keys_optional(:5))
        ! execute
        call xfind_nnimgs_distr%execute(cline)

    ! PRIME3D

    case('prime3D_init')
        !==Program prime3D_init
        !
        ! <prime3D_init/begin>is a program for generating a random initial model for initialisation of PRIME3D
        ! <prime3D_init/end> 
        !
        ! set required keys
        keys_required(1)  = 'stk'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'msk'
        keys_required(4)  = 'ctf'
        keys_required(5)  = 'pgrp' 
        keys_required(6)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'ncunits'
        keys_optional(3)  = 'deftab'
        keys_optional(4)  = 'lp'
        keys_optional(5)  = 'inner'
        keys_optional(6)  = 'width'
        keys_optional(7)  = 'nspace'
        keys_optional(8)  = 'nran'
        keys_optional(9)  = 'npeaks'
        keys_optional(10) = 'xfel'        
        ! parse command line
        if( describe ) call print_doc_prime3D_init
        call cline%parse(keys_required(:6), keys_optional(:10))
        ! set defaults
        if( .not. cline%defined('nspace') ) call cline%set('nspace', 1000.)
        ! execute
        call xprime3D_init_distr%execute( cline )
    case('prime3D')
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
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'msk'
        keys_required(4)  = 'ctf'
        keys_required(5)  = 'pgrp'
        keys_required(6)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'ncunits'
        keys_optional(3)  = 'deftab'
        keys_optional(4)  = 'vol1'
        keys_optional(5)  = 'oritab'
        keys_optional(6)  = 'trs'
        keys_optional(7)  = 'hp'
        keys_optional(8)  = 'lp'
        keys_optional(9)  = 'cenlp'
        keys_optional(10) = 'dynlp'
        keys_optional(11) = 'lpstart'
        keys_optional(12) = 'lpstop'
        keys_optional(13) = 'eo'
        keys_optional(14) = 'refine'
        keys_optional(15) = 'frac'
        keys_optional(16) = 'automsk'
        keys_optional(17) = 'mw'
        keys_optional(18) = 'amsklp'
        keys_optional(19) = 'edge'
        keys_optional(20) = 'binwidth'
        keys_optional(21) = 'inner'
        keys_optional(22) = 'width'
        keys_optional(23) = 'nspace'
        keys_optional(24) = 'nstates'
        keys_optional(25) = 'npeaks'
        keys_optional(26) = 'startit'
        keys_optional(27) = 'maxits'
        keys_optional(28) = 'shbarrier'
        keys_optional(29) = 'noise'
        keys_optional(30) = 'xfel'
        keys_optional(31) = 'nnn'
        keys_optional(32) = 'shellw'
        keys_optional(33) = 'rrate'
        keys_optional(34) = 'norec'
        ! documentation
        if( describe ) call print_doc_prime3D
        ! parse command line
        call check_restart( entire_line, is_restart )
        if( is_restart )then
            call parse_restart('prime3D', entire_line, cline, keys_required(:6), keys_optional(:34))
        else
            call cline%parse( keys_required(:6), keys_optional(:34) )
        endif
        ! set defaults
        if( .not. cline%defined('nspace')                  ) call cline%set('nspace', 1000.)
        if( cline%defined('lp') .or. cline%defined('find') ) call cline%set('dynlp',   'no')
        if( .not. cline%defined('cenlp')                   ) call cline%set('cenlp',    30.)
        if( .not. cline%defined('refine')                  ) call cline%set('refine',  'no')
        if( .not. cline%defined('eo') )then
            call cline%set('eo', 'no')
        else
            if( cline%get_carg('eo').eq.'yes' )call cline%set('dynlp','no')
        endif
        ! execute
        call xprime3D_distr%execute(cline)
    case('cont3D')
        !==Program prime3D
        !
        ! <cont3D/begin><cont3D/end>
        !
        ! set required keys
        keys_required(1)  = 'stk'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'msk'
        keys_required(4)  = 'ctf'
        keys_required(5)  = 'pgrp'
        keys_required(6)  = 'oritab'
        keys_required(7)  = 'trs'
        keys_required(8)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'deftab'
        keys_optional(3)  = 'vol1'
        keys_optional(4)  = 'hp'
        keys_optional(5)  = 'lpstart'
        keys_optional(6)  = 'lpstop'
        keys_optional(7)  = 'frac'
        keys_optional(8)  = 'automsk'
        keys_optional(9)  = 'mw'
        keys_optional(10) = 'amsklp'
        keys_optional(11) = 'edge'
        keys_optional(12) = 'inner'
        keys_optional(13) = 'width'
        keys_optional(14) = 'startit'
        keys_optional(15) = 'maxits'
        keys_optional(16) = 'xfel'
        keys_optional(17) = 'shellw'
        ! documentation
        if( describe ) call print_doc_cont3D
        ! parse command line
        call check_restart( entire_line, is_restart )
        if( is_restart )then
            call parse_restart('cont3D', entire_line, cline, keys_required(:8), keys_optional(:17))
        else
            call cline%parse( keys_required(:8), keys_optional(:17) )
        endif
        ! set defaults
        if( .not. cline%defined('eo') )then
            call cline%set('eo', 'no')
        endif
        call cline%set('dynlp', 'no')
        if(.not.cline%defined('nspace'))call cline%set('nspace', 1000.)
        if(.not.cline%defined('shellw'))call cline%set('shellw', 'no')
        if(.not.cline%defined('refine'))call cline%set('refine', 'no')
        ! execute
        call xcont3D_distr%execute(cline)
    case('shellweight3D')
        !==Program shellweight3D
        !
        ! <shellweight3D/begin>is a program for calculating the shell-by-shell resolution weights in a global sense, so that 
        ! particles that do contribute with higher resolution information (as measure by the FRC) are given the appropriate 
        ! weight<shellweight3D/end>
        !
        ! set required keys     
        keys_required(1)  = 'stk'
        keys_required(2)  = 'vol1'
        keys_required(3)  = 'smpd'
        keys_required(4)  = 'msk' 
        keys_required(5)  = 'oritab'
        keys_required(6)  = 'ctf'
        keys_required(7)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'ncunits'
        keys_optional(3)  = 'deftab'
        keys_optional(4)  = 'automsk'
        keys_optional(5)  = 'mw'
        keys_optional(6)  = 'amsklp'
        keys_optional(7)  = 'edge'
        keys_optional(8)  = 'binwidth'
        keys_optional(9)  = 'inner'
        keys_optional(10) = 'width'
        keys_optional(11) = 'refine'       
        ! parse command line
        if( describe ) call print_doc_shellweight3D
        call cline%parse(keys_required(:7), keys_optional(:11))
        ! execute
        call xshellweight3D_distr%execute(cline)
        ! set defaults
        call cline%set('outfile', 'shellweight3D_doc.txt')
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
        keys_required(1)  = 'stk'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'oritab'
        keys_required(4)  = 'msk'
        keys_required(5)  = 'ctf'
        keys_required(6)  = 'pgrp'
        keys_required(7)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'ncunits'
        keys_optional(3)  = 'eo'
        keys_optional(4)  = 'deftab'
        keys_optional(5)  = 'frac'
        keys_optional(6)  = 'mw'
        keys_optional(7)  = 'mul'
        keys_optional(8)  = 'state'
        keys_optional(9)  = 'shellw'
        keys_optional(10) = 'vol1'
        keys_optional(11) = 'npeaks'
        ! parse command line
        if( describe ) call print_doc_recvol
        call cline%parse(keys_required(:7), keys_optional(:11))
        ! set defaults
        if( .not. cline%defined('trs') ) call cline%set('trs', 5.) ! to assure that shifts are being used
        if( .not. cline%defined('eo')  ) call cline%set('eo', 'no')
        ! execute
        call xrecvol_distr%execute( cline )

    ! TIME-SERIES DISTRIBUTED WORKFLOWS

    case( 'tseries_track' )
        !==Program tseries_track
        !
        ! <tseries_extract/begin>is a program for particle tracking in time-series data
        ! <tseries_extract/end> 
        !
        ! set required keys
        keys_required(1) = 'filetab'
        keys_required(2) = 'fbody'
        keys_required(3) = 'smpd'
        keys_required(4) = 'boxfile'
        keys_required(5) = 'ncunits'
        ! set optional keys
        keys_optional(1) = 'lp'
        keys_optional(2) = 'offset'
        ! parse command line
        ! if( describe ) call print_doc_tseries_track
        call cline%parse(keys_required(:5), keys_optional(:2))
        ! set defaults
        call cline%set('nthr', 1.0)
        if( .not. cline%defined('neg') ) call cline%set('neg', 'yes')
        if( .not. cline%defined('lp')  ) call cline%set('lp',    2.0)
        ! execute
        call xtseries_track_distr%execute( cline )

    ! HIGH-LEVEL DISTRIBUTED WORKFLOWS

    case( 'isw' )
        !==Program isw
        !
        ! <isw/begin>is a program for multi-particle 3D reconstruction by iterated spectral weights,
        ! a blind deconvolution approach for weighted state sorting<isw/end> 
        !
        ! set required keys
        keys_required(1)  = 'stk'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'oritab'
        keys_required(4)  = 'msk'
        keys_required(5)  = 'ctf'
        keys_required(6)  = 'pgrp'
        keys_required(7)  = 'nstates'
        keys_required(8)  = 'npeaks'
        keys_required(9)  = 'lp'
        keys_required(10) = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'ncunits'
        keys_optional(3)  = 'deftab'
        keys_optional(4)  = 'frac'
        keys_optional(5)  = 'maxits'
        ! parse command line
        ! if( describe ) call print_doc_isw
        call cline%parse(keys_required(:10), keys_optional(:5))
        ! set defaults
        call cline%set('refine', 'isw')
        ! execute
        call xisw_distr%execute( cline )
    case( 'ini3D_from_cavgs' )
        !==Program ini3D_from_cavgs
        !
        ! <ini3D_from_cavgs/begin>is a program for generating an initial 3D model from class averages
        ! obtained with prime2D<ini3D_from_cavgs/end> 
        !
        ! set required keys
        keys_required(1)  = 'stk'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'msk'
        keys_required(4)  = 'pgrp'
        keys_required(5)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'nthr_master'
        keys_optional(3)  = 'ncunits'
        keys_optional(4)  = 'hp'
        keys_optional(5)  = 'lp'
        keys_optional(6)  = 'frac'
        keys_optional(7)  = 'automsk'
        keys_optional(8)  = 'mw'
        keys_optional(9)  = 'amsklp'
        keys_optional(10) = 'edge'
        keys_optional(11) = 'binwidth'
        keys_optional(12) = 'inner'
        keys_optional(13) = 'width'
        keys_optional(14) = 'nspace'
        keys_optional(15) = 'shbarrier'
        ! parse command line
        if( describe ) call print_doc_ini3D_from_cavgs
        call cline%parse(keys_required(:5), keys_optional(:15))
        ! execute
        call xini3D_from_cavgs%execute( cline )
    case( 'het_ensemble' )
        !==Program het_ensemble
        !
        ! <het_ensemble/begin><het_ensemble/end> 
        !
        ! set required keys
        keys_required(1)  = 'stk'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'oritab'
        keys_required(4)  = 'msk'
        keys_required(5)  = 'pgrp'
        keys_required(6)  = 'ctf'
        keys_required(7)  = 'nstates'
        keys_required(8)  = 'lp'
        keys_required(9)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'frac'
        keys_optional(3)  = 'automsk'
        keys_optional(4)  = 'mw'
        keys_optional(5)  = 'amsklp'
        keys_optional(6)  = 'edge'
        keys_optional(7)  = 'binwidth'
        keys_optional(8)  = 'inner'
        keys_optional(9)  = 'width'
        keys_optional(10) = 'nspace'
        ! parse command line
        ! if( describe ) call print_doc_het_ensemble
        call cline%parse(keys_required(:9), keys_optional(:10))
        if( .not. cline%defined('shellw')  ) call cline%set('shellw', 'no')
        ! execute
        call xhet_ensemble%execute( cline )
    case DEFAULT
        write(*,'(a,a)') 'program key (prg) is: ', trim(prg)
        stop 'unsupported program'
    end select
end program simple_distr_exec
