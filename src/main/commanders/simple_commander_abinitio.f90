! concrete commander: high-level workflows
module simple_commander_abinitio
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_parameters,         only: parameters, params_glob
use simple_sp_project,         only: sp_project
use simple_stack_io,           only: stack_io
use simple_qsys_env,           only: qsys_env
use simple_commander_base,     only: commander_base
use simple_commander_volops,   only: reproject_commander, symaxis_search_commander, postprocess_commander
use simple_commander_rec,      only: reconstruct3D_commander, reconstruct3D_commander_distr
use simple_commander_refine3D, only: refine3D_commander, refine3D_commander_distr
use simple_procimgstk,         only: shift_imgfile
use simple_image,              only: image
use simple_builder,            only: builder
use simple_class_frcs,         only: class_frcs
use simple_commander_euclid
use simple_euclid_sigma2
use simple_qsys_funs
implicit none

public :: initial_3Dmodel_commander, abinitio_3Dmodel_autolp_commander, abinitio_3Dmodel_commander
public :: batch_abinitio_3Dmodel_commander
public :: abinitio_3Dmodel2_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: initial_3Dmodel_commander
    contains
    procedure :: execute => exec_initial_3Dmodel
end type initial_3Dmodel_commander

type, extends(commander_base) :: abinitio_3Dmodel_autolp_commander
    contains
    procedure :: execute => exec_abinitio_3Dmodel_autolp
end type abinitio_3Dmodel_autolp_commander

type, extends(commander_base) :: abinitio_3Dmodel_commander
    contains
    procedure :: execute => exec_abinitio_3Dmodel
end type abinitio_3Dmodel_commander

type, extends(commander_base) :: abinitio_3Dmodel2_commander
    contains
    procedure :: execute => exec_abinitio_3Dmodel2
end type abinitio_3Dmodel2_commander

type, extends(commander_base) :: batch_abinitio_3Dmodel_commander
    contains
    procedure :: execute => exec_batch_abinitio_3Dmodel
end type batch_abinitio_3Dmodel_commander

character(len=STDLEN), parameter :: REC_FBODY            = 'rec_final'
character(len=STDLEN), parameter :: REC_PPROC_FBODY      = trim(REC_FBODY)//trim(PPROC_SUFFIX)
character(len=STDLEN), parameter :: REC_PPROC_MIRR_FBODY = trim(REC_PPROC_FBODY)//trim(MIRR_SUFFIX)

contains

    ! !> for generation of an initial 3d model from class averages
    ! subroutine exec_initial_3Dmodel( self, cline )
    !     class(initial_3Dmodel_commander), intent(inout) :: self
    !     class(cmdline),                   intent(inout) :: cline
    !     ! constants
    !     real,                  parameter :: CENLP_DEFAULT   = 30.
    !     real,                  parameter :: STARTLP_DEFAULT = 20.
    !     real,                  parameter :: LP_SYMSRCH_LB   = 12.
    !     integer,               parameter :: MAXITS1         = 50
    !     integer,               parameter :: MAXITS2         = 30
    !     integer,               parameter :: NSPACE1         = 500
    !     integer,               parameter :: NSPACE2         = 1000
    !     integer,               parameter :: MINBOX          = 88
    !     character(len=STDLEN), parameter :: work_projfile   = 'initial_3Dmodel_tmpproj.simple'
    !     logical,               parameter :: L_FSC_STAGE2    = .true.
    !     logical,               parameter :: L_EXPGATE_MODE  = .false.
    !     ! distributed commanders
    !     type(calc_pspec_commander_distr) :: xcalc_pspec_distr
    !     ! shared-mem commanders
    !     type(refine3D_commander)         :: xrefine3D
    !     type(reconstruct3D_commander)    :: xreconstruct3D
    !     type(symaxis_search_commander)   :: xsymsrch
    !     type(reproject_commander)        :: xreproject
    !     type(postprocess_commander)      :: xpostprocess
    !     ! command lines
    !     type(cmdline) :: cline_refine3D_1, cline_refine3D_2
    !     type(cmdline) :: cline_symsrch
    !     type(cmdline) :: cline_reconstruct3D, cline_postprocess
    !     type(cmdline) :: cline_reproject, cline_calc_pspec
    !     ! other
    !     character(len=:), allocatable :: stk, stkpath, orig_stk, frcs_fname, shifted_stk, stk_even, stk_odd, ext
    !     integer,          allocatable :: states(:)
    !     real,             allocatable :: frcs_avg(:)
    !     character(len=2)      :: str_state
    !     type(ori)             :: o, o_even, o_odd
    !     type(qsys_env)        :: qenv
    !     type(parameters)      :: params
    !     type(ctfparams)       :: ctfvars
    !     type(sp_project)      :: spproj, work_proj
    !     type(sym)             :: se1,se2
    !     type(image)           :: img, vol
    !     type(stack_io)        :: stkio_r, stkio_r2, stkio_w
    !     type(class_frcs)      :: clsfrcs
    !     character(len=STDLEN) :: vol_iter, pgrp_init, pgrp_refine, vol_iter_pproc, vol_iter_pproc_mirr
    !     character(len=STDLEN) :: sigma2_fname, sigma2_fname_sc, orig_objfun, frckind
    !     real                  :: scale_factor1, scale_factor2, trslim, smpd_target, lplims(2), lp_est, lp, lp_sym
    !     integer               :: icls, ncavgs, cnt, iter, ipart, even_ind, odd_ind, state, find_start, find, filtsz
    !     logical               :: srch4symaxis, do_autoscale, symran_before_refine
    !     call cline%set('objfun',    'euclid') ! use noise normalized Euclidean distances from the start
    !     call cline%set('sigma_est', 'global') ! obviously
    !     call cline%set('oritype',      'out') ! because cavgs are part of out segment
    !     call cline%set('bfac',            0.) ! because initial models should not be sharpened
    !     call cline%set('ml_reg',        'no') ! not trusting FSC-based SSNR estimation on class averages by default
    !     if( .not. cline%defined('mkdir')       ) call cline%set('mkdir',      'yes')
    !     if( .not. cline%defined('autoscale')   ) call cline%set('autoscale',  'yes')
    !     if( .not. cline%defined('overlap')     ) call cline%set('overlap',     0.99) ! needed to prevent premature convergence
    !     if( L_EXPGATE_MODE )then
    !     if( .not. cline%defined('prob_athres') ) call cline%set('prob_athres',  10.) ! what Cong uses for exp gate
    !     if( .not. cline%defined('center')      ) call cline%set('center',      'no') ! what Cong uses for exp gate
    !     else
    !     if( .not. cline%defined('prob_athres') ) call cline%set('prob_athres',  90.) ! reduces # failed runs on trpv1 from 4->2/10
    !     endif
    !     if( .not. cline%defined('cenlp')       ) call cline%set('cenlp', CENLP_DEFAULT)
    !     if( .not. cline%defined('imgkind')     ) call cline%set('imgkind',   'cavg') ! whether to use classes generated from 2D/3D
    !     ! make master parameters
    !     call params%new(cline)
    !     call cline%delete('autoscale')
    !     call cline%set('lp_auto',    'yes')    ! automated frequency marching is the default method of choice
    !     call cline%set('mkdir',       'no')    ! to avoid nested directory structure
    !     call cline%set('oritype', 'ptcl3D')    ! from now on we are in the ptcl3D segment, final report is in the cls3D segment
    !     ! state string
    !     str_state    = int2str_pad(1,2)
    !     ! decide wether to search for the symmetry axis
    !     pgrp_init    = trim(params%pgrp_start)
    !     pgrp_refine  = trim(params%pgrp)
    !     srch4symaxis = trim(pgrp_refine) .ne. trim(pgrp_init)
    !     symran_before_refine = .false.
    !     if( pgrp_init.ne.'c1' .or. pgrp_refine.ne.'c1' )then
    !         se1 = sym(pgrp_init)
    !         se2 = sym(pgrp_refine)
    !         if(se1%get_nsym() > se2%get_nsym())then
    !             ! ensure se2 is a subgroup of se1
    !             if( .not. se1%has_subgrp(pgrp_refine) ) THROW_HARD('Incompatible symmetry groups; simple_commander_abinitio2')
    !             ! set flag for symmetry randomisation before refinmement
    !             ! in case we are moving from a higher to lower group
    !             symran_before_refine = .true.
    !         else if( se2%get_nsym() > se1%get_nsym() )then
    !             ! ensure se1 is a subgroup of se2
    !             if( .not. se2%has_subgrp(pgrp_init) ) THROW_HARD('Incompatible symmetry groups; simple_commander_abinitio2')
    !         endif
    !     endif
    !     ! read project & update sampling distance
    !     call spproj%read(params%projfile)
    !     call spproj%update_projinfo(cline)
    !     call spproj%write_segment_inside('projinfo', params%projfile)
    !     ! whether to use classes generated from 2D or 3D
    !     select case(trim(params%imgkind))
    !         case('cavg')
    !             frckind = 'frc2D'
    !             states  = nint(spproj%os_cls2D%get_all('state'))
    !         case('cavg3D')
    !             frckind = 'frc3D'
    !             states  = nint(spproj%os_cls3D%get_all('state'))
    !         case DEFAULT
    !             THROW_HARD('Unsupported IMGKIND!')
    !     end select
    !     call cline%delete('imgkind') ! no interference down the line
    !     ! retrieve cavgs stack info
    !     call spproj%get_cavgs_stk(stk, ncavgs, params%smpd, imgkind=params%imgkind, stkpath=stkpath)
    !     if(.not. file_exists(stk)) stk = trim(stkpath) // '/' // trim(stk)
    !     if(.not. file_exists(stk)) THROW_HARD('cavgs stk does not exist; simple_commander_abinitio')
    !     orig_stk        = stk
    !     ext             = '.'//fname2ext(stk)
    !     stk_even        = add2fbody(trim(stk), trim(ext), '_even')
    !     stk_odd         = add2fbody(trim(stk), trim(ext), '_odd')
    !     if( .not. file_exists(stk_even) ) THROW_HARD('Even cavgs stk: '//trim(stk_even)//' does not exist!')
    !     if( .not. file_exists(stk_odd)  ) THROW_HARD('Odd cavgs stk: '//trim(stk_odd)//' does not exist!')
    !     ctfvars%ctfflag = CTFFLAG_NO
    !     ctfvars%smpd    = params%smpd
    !     shifted_stk     = basename(add2fbody(stk, ext, '_shifted'))
    !     if( count(states==0) .eq. ncavgs )then
    !         THROW_HARD('no class averages detected in project file: '//trim(params%projfile)//'; initial_3Dmodel')
    !     endif
    !     ! retrieve FRC info
    !     call spproj%get_frcs(frcs_fname, frckind, fail=.false.)
    !     if( .not.file_exists(frcs_fname) )then
    !         ! 08/24 This is a backwards compatibility patch to account for error in metadata
    !         ! on exit of streaming related to GUI directory structure (now fixed and cf above get_cavgs_stk).
    !         ! Will need to harmonize (move to absolute path?).
    !         frcs_fname = trim(stkpath)//'/'//trim(frcs_fname)
    !         if( .not.file_exists(frcs_fname) )then
    !             THROW_HARD('the project file does not contain an FRCs file, which is required')
    !         endif
    !     endif
    !     params%frcs = trim(frcs_fname)
    !     call clsfrcs%read(frcs_fname)
    !     filtsz = clsfrcs%get_filtsz()
    !     allocate(frcs_avg(filtsz), source=0.)
    !     call clsfrcs%avg_frc_getter(frcs_avg, states)
    !     find = get_find_at_corr(frcs_avg, 0.9)
    !     lp   = calc_lowpass_lim(find, params%box, params%smpd)
    !     print *, 'lp@0.9   ', lp
    !     find = get_find_at_corr(frcs_avg, 0.8)
    !     lp   = calc_lowpass_lim(find, params%box, params%smpd)
    !     print *, 'lp@0.8   ', lp
    !     find = get_find_at_corr(frcs_avg, 0.7)
    !     lp   = calc_lowpass_lim(find, params%box, params%smpd)
    !     print *, 'lp@0.7   ', lp
    !     find = get_find_at_corr(frcs_avg, 0.6)
    !     lp   = calc_lowpass_lim(find, params%box, params%smpd)
    !     print *, 'lp@0.6   ', lp
    !     find = get_find_at_corr(frcs_avg, 0.5)
    !     lp   = calc_lowpass_lim(find, params%box, params%smpd)
    !     print *, 'lp@0.5   ', lp
    !     find = get_find_at_corr(frcs_avg, 0.4)
    !     lp   = calc_lowpass_lim(find, params%box, params%smpd)
    !     print *, 'lp@0.4   ', lp
    !     find = get_find_at_corr(frcs_avg, 0.3)
    !     lp   = calc_lowpass_lim(find, params%box, params%smpd)
    !     print *, 'lp@0.3   ', lp
    !     find = get_find_at_corr(frcs_avg, 0.143)
    !     lp   = calc_lowpass_lim(find, params%box, params%smpd)
    !     print *, 'lp@0.143 ', lp
    !     ! prepare a temporary project file
    !     work_proj%projinfo = spproj%projinfo
    !     work_proj%compenv  = spproj%compenv
    !     if( spproj%jobproc%get_noris()  > 0 ) work_proj%jobproc = spproj%jobproc
    !     ! name change
    !     call work_proj%projinfo%delete_entry('projname')
    !     call work_proj%projinfo%delete_entry('projfile')
    !     call cline%set('projfile', trim(work_projfile))
    !     call cline%set('projname', trim(get_fbody(trim(work_projfile),trim('simple'))))
    !     call work_proj%update_projinfo(cline)
    !     ! add stks to temporary project
    !     call work_proj%add_stk(stk_even, ctfvars)
    !     call work_proj%add_stk(stk_odd,  ctfvars)
    !     ! update orientations parameters
    !     do icls=1,ncavgs
    !         even_ind = icls
    !         odd_ind  = ncavgs + icls
    !         call work_proj%os_ptcl3D%get_ori(icls, o)
    !         call o%set('class', real(icls))
    !         call o%set('state', real(states(icls)))
    !         ! even
    !         o_even = o
    !         call o_even%set('eo', 0.)
    !         call o_even%set('stkind', work_proj%os_ptcl3D%get(even_ind,'stkind'))
    !         call work_proj%os_ptcl3D%set_ori(even_ind, o_even)
    !         ! odd
    !         o_odd = o
    !         call o_odd%set('eo', 1.)
    !         call o_odd%set('stkind', work_proj%os_ptcl3D%get(odd_ind,'stkind'))
    !         call work_proj%os_ptcl3D%set_ori(odd_ind, o_odd)
    !     enddo
    !     params_glob%nptcls = work_proj%get_nptcls()
    !     call work_proj%write()
    !     ! set lplims
    !     if( cline%defined('lpstart') )then
    !         lplims(1) = params%lpstart
    !     else if( any(frcs_avg > 0.8) )then
    !         find      = get_find_at_corr(frcs_avg, 0.8)
    !         lplims(1) = calc_lowpass_lim(find, params%box, params%smpd)
    !     else
    !         lplims(1) = STARTLP_DEFAULT
    !     endif
    !     if( any(frcs_avg > 0.5) )then
    !         find      = get_find_at_corr(frcs_avg, 0.5)
    !         lplims(2) = calc_lowpass_lim(find, params%box, params%smpd)
    !     else
    !         lplims(2) = STARTLP_DEFAULT
    !     endif
    !     write(logfhandle,'(A,F5.1)') '>>> DID SET STARTING  LOW-PASS LIMIT (IN A) TO: ', lplims(1)
    !     write(logfhandle,'(A,F5.1)') '>>> DID SET HARD      LOW-PASS LIMIT (IN A) TO: ', lplims(2)
    !     write(logfhandle,'(A,F5.1)') '>>> DID SET CENTERING LOW-PASS LIMIT (IN A) TO: ', params_glob%cenlp
    !     ! initial phase scaling
    !     smpd_target = max(params%smpd, lplims(2)*LP2SMPDFAC)
    !     call downscale(smpd_target, scale_factor1)
    !     ! prepare command lines from prototype
    !     ! projects names are subject to change depending on scaling and are updated individually
    !     call cline%delete('projname')
    !     call cline%delete('projfile')
    !     cline_reconstruct3D  = cline
    !     cline_refine3D_1     = cline
    !     cline_refine3D_2     = cline
    !     cline_reproject      = cline
    !     cline_symsrch        = cline
    !     ! initialise command line parameters
    !     ! (1) PROBABILISTIC AB INITIO STEP
    !     call cline_refine3D_1%set('prg',                 'refine3D')
    !     call cline_refine3D_1%set('projfile',   trim(work_projfile))
    !     call cline_refine3D_1%set('box_crop', real(params%box_crop))
    !     call cline_refine3D_1%set('smpd_crop',     params%smpd_crop)
    !     call cline_refine3D_1%set('lpstart',              lplims(1))
    !     call cline_refine3D_1%set('lpstop',               lplims(2))
    !     call cline_refine3D_1%set('nspace',           real(NSPACE1))
    !     call cline_refine3D_1%set('maxits',           real(MAXITS1))
    !     call cline_refine3D_1%set('silence_fsc',              'yes') ! no FSC plot printing in prob phase
    !     if( L_EXPGATE_MODE )then
    !     call cline_refine3D_1%set('refine',                  'prob')
    !     call cline_refine3D_1%set('lp_auto',                  'yes')
    !     call cline_refine3D_1%set('sh_first',                  'no')
    !     call cline_refine3D_1%set('prob_sh',                   'no')
    !     call cline_refine3D_1%set('ml_reg',                   'yes')
    !     call cline_refine3D_1%set('trs',                        0.0)
    !     else
    !     call cline_refine3D_1%set('refine',              'shc_smpl') ! best refine mode identified for class averages
    !     call cline_refine3D_1%set('lp_auto',                  'yes')
    !     call cline_refine3D_1%set('sh_first',                 'yes')
    !     call cline_refine3D_1%set('prob_sh',                   'no')
    !     call cline_refine3D_1%set('snr_noise_reg',              2.0)
    !     call cline_refine3D_1%set('ml_reg',                    'no')
    !     call cline_refine3D_1%set('trs',                     trslim) ! trslim set in call downscale(smpd_target, scale_factor1) above
    !     endif
    !     ! (2) SYMMETRY AXIS SEARCH
    !     if( srch4symaxis )then
    !         ! need to replace original point-group flag with c1/pgrp_start
    !         call cline_refine3D_1%set('pgrp', trim(pgrp_init))
    !         ! symsrch
    !         call qenv%new(1, exec_bin='simple_exec')
    !         call cline_symsrch%set('prg',         'symaxis_search') ! needed for cluster exec
    !         call cline_symsrch%set('pgrp',       trim(pgrp_refine))
    !         call cline_symsrch%set('smpd',        params%smpd_crop)
    !         call cline_symsrch%set('box',    real(params%box_crop))
    !         call cline_symsrch%set('projfile', trim(work_projfile))
    !         if( .not. cline_symsrch%defined('cenlp') ) call cline_symsrch%set('cenlp', CENLP_DEFAULT)
    !         call cline_symsrch%set('hp',                 params%hp)
    !         call cline_symsrch%set('oritype',             'ptcl3D')
    !         call cline_symsrch%delete('lp_auto')
    !     endif
    !     ! (3)  REFINEMENT
    !     call cline_refine3D_2%set('prg',               'refine3D')
    !     call cline_refine3D_2%set('projfile', trim(work_projfile))
    !     ! box_crop & smpd_crop set after downscaling, below
    !     ! lpstart & lpstop set after low-pass limit estimation, below
    !     call cline_refine3D_2%set('pgrp',       trim(pgrp_refine))
    !     call cline_refine3D_2%set('nspace',         real(NSPACE2))
    !     call cline_refine3D_2%set('maxits',         real(MAXITS2))
    !     call cline_refine3D_2%set('silence_fsc',             'no')
    !     call cline_refine3D_2%set('refine',                'prob')
    !     call cline_refine3D_2%set('sh_first',               'yes')
    !     call cline_refine3D_2%set('prob_sh',                'yes')
    !     if( L_FSC_STAGE2 )then
    !         call cline_refine3D_2%set('lp_auto',            'fsc')
    !         call cline_refine3D_2%set('lplim_crit',         0.143)
    !         call cline_refine3D_2%set('ml_reg',             'yes')
    !     else
    !         call cline_refine3D_2%set('lp_auto',            'yes') ! lpstart/lpstop set below
    !         call cline_refine3D_2%set('snr_noise_reg',        3.0)
    !         call cline_refine3D_2%set('ml_reg',              'no')
    !     endif
    !     ! trslim set after downscaling, below
    !     ! (4) RE-CONSTRUCT & RE-PROJECT VOLUME
    !     call cline_reconstruct3D%set('prg',       'reconstruct3D')
    !     call cline_reconstruct3D%set('box',      real(params%box))
    !     call cline_reconstruct3D%set('projfile',    work_projfile)
    !     call cline_reconstruct3D%set('needs_sigma',         'yes')
    !     call cline_postprocess%set('prg',           'postprocess')
    !     call cline_postprocess%set('projfile',      work_projfile)
    !     call cline_postprocess%set('mkdir',                  'no')
    !     call cline_postprocess%delete('bfac') ! sharpen final map
    !     call cline_reproject%set('prg',               'reproject')
    !     call cline_reproject%set('pgrp',        trim(pgrp_refine))
    !     call cline_reproject%set('outstk',         'reprojs'//ext)
    !     call cline_reproject%set('smpd',              params%smpd)
    !     call cline_reproject%set('box',          real(params%box))
    !     ! execute commanders
    !     write(logfhandle,'(A)') '>>>'
    !     write(logfhandle,'(A)') '>>> BAYESIAN 3D AB INITIO'
    !     write(logfhandle,'(A)') '>>>'
    !     call rndstart(cline_refine3D_1)
    !     call xrefine3D%execute_shmem(cline_refine3D_1)
    !     iter       = nint(cline_refine3D_1%get_rarg('endit'))
    !     vol_iter   = trim(VOL_FBODY)//trim(str_state)//ext
    !     call work_proj%read_segment('ptcl3D', work_projfile)
    !     lp_est     = work_proj%os_ptcl3D%get_avg('lp')
    !     find_start = calc_fourier_index(lp_est, params%box_crop, params%smpd_crop) - 2
    !     lplims(1)  = calc_lowpass_lim(find_start, params%box_crop, params%smpd_crop)
    !     lplims(2)  = calc_lplim_stage2(3) ! low-pass limit is median of three best (as in 2D)
    !     call cline_refine3D_2%set('lpstart', lplims(1))
    !     call cline_refine3D_2%set('lpstop',  lplims(2))
    !     write(logfhandle,'(A,F5.1)') '>>> ESTIMATED         LOW-PASS LIMIT (IN A) TO: ', lp_est
    !     write(logfhandle,'(A,F5.1)') '>>> DID SET STARTING  LOW-PASS LIMIT (IN A) TO: ', lplims(1)
    !     write(logfhandle,'(A,F5.1)') '>>> DID SET HARD      LOW-PASS LIMIT (IN A) TO: ', lplims(2)
    !     if( .not. file_exists(vol_iter) ) THROW_HARD('input volume to symmetry axis search does not exist')
    !     if( symran_before_refine )then
    !         call se1%symrandomize(work_proj%os_ptcl3D)
    !         call work_proj%write_segment_inside('ptcl3D', work_projfile)
    !     endif
    !     if( srch4symaxis )then
    !         lp_sym = max(LP_SYMSRCH_LB,lp_est)
    !         write(logfhandle,'(A,F5.1)') '>>> DID SET SYMSEARCH LOW-PASS LIMIT (IN A) TO: ', lp_sym
    !         call cline_symsrch%set('lp', lp_sym)
    !         write(logfhandle,'(A)') '>>>'
    !         write(logfhandle,'(A)') '>>> SYMMETRY AXIS SEARCH'
    !         write(logfhandle,'(A)') '>>>'
    !         call cline_symsrch%set('vol1', trim(vol_iter))
    !         call xsymsrch%execute_shmem(cline_symsrch)
    !         call del_file('SYMAXIS_SEARCH_FINISHED')
    !     endif
    !     ! refinement scaling
    !     smpd_target = max(params%smpd, lplims(2)*LP2SMPDFAC)
    !     call downscale(smpd_target, scale_factor2)
    !     ! refinement stage
    !     write(logfhandle,'(A)') '>>>'
    !     write(logfhandle,'(A)') '>>> PROBABLISTIC 3D REFINEMENT'
    !     write(logfhandle,'(A)') '>>>'
    !     iter = iter + 1
    !     call cline_refine3D_2%set('startit',             real(iter))
    !     call cline_refine3D_2%set('box_crop', real(params%box_crop))
    !     call cline_refine3D_2%set('smpd_crop',     params%smpd_crop)
    !     call cline_refine3D_2%set('trs',                     trslim) ! trslim set in call downscale(smpd_target, scale_factor2) above
    !     call cline_refine3D_2%set('which_iter',          real(iter))
    !     ! sigma2 at original sampling
    !     cline_calc_pspec = cline
    !     call cline_calc_pspec%delete('lp_auto')
    !     call cline_calc_pspec%set('prg',      'calc_pspec' )
    !     call cline_calc_pspec%set('projfile', work_projfile)
    !     call cline_calc_pspec%set('box',   real(params%box))
    !     call cline_calc_pspec%set('smpd',       params%smpd)
    !     call cline_calc_pspec%set('which_iter',  real(iter))
    !     call xcalc_pspec_distr%execute_shmem(cline_calc_pspec)
    !     call rec(cline_refine3D_2)
    !     call xrefine3D%execute_shmem(cline_refine3D_2)
    !     iter = nint(cline_refine3D_2%get_rarg('endit'))
    !     call cline_reconstruct3D%set('which_iter',real(iter))
    !     ! deals with final volume
    !     if( do_autoscale )then
    !         write(logfhandle,'(A)') '>>>'
    !         write(logfhandle,'(A)') '>>> RECONSTRUCTION AT ORIGINAL SAMPLING'
    !         write(logfhandle,'(A)') '>>>'
    !         call cline_reconstruct3D%delete('lp_auto')
    !         call cline_reconstruct3D%set('box',  real(params%box))
    !         call cline_reconstruct3D%set('smpd', params%smpd)
    !         ! reconstruction
    !         call xreconstruct3D%execute_shmem(cline_reconstruct3D)
    !         vol_iter = trim(VOL_FBODY)//trim(str_state)//ext
    !         ! because postprocess only updates project file when mkdir=yes
    !         call work_proj%read_segment('out', work_projfile)
    !         call work_proj%add_vol2os_out(vol_iter, params%smpd, 1, 'vol')
    !         call work_proj%add_fsc2os_out(FSC_FBODY//str_state//trim(BIN_EXT), 1, params%box)
    !         call work_proj%write_segment_inside('out', work_projfile)
    !         call xpostprocess%execute(cline_postprocess)
    !     else
    !         vol_iter = trim(VOL_FBODY)//trim(str_state)//ext
    !         call vol%new([params%box,params%box,params%box],params%smpd)
    !         call vol%read(vol_iter)
    !         call vol%mirror('x')
    !         call vol%write(add2fbody(vol_iter,ext,trim(PPROC_SUFFIX)//trim(MIRR_SUFFIX)))
    !         call vol%kill
    !     endif
    !     vol_iter_pproc      = add2fbody(vol_iter,ext,PPROC_SUFFIX)
    !     vol_iter_pproc_mirr = add2fbody(vol_iter,ext,trim(PPROC_SUFFIX)//trim(MIRR_SUFFIX))
    !     if( file_exists(vol_iter)            ) call simple_rename(vol_iter,            trim(REC_FBODY)//ext)
    !     if( file_exists(vol_iter_pproc)      ) call simple_rename(vol_iter_pproc,      trim(REC_PPROC_FBODY)//ext)
    !     if( file_exists(vol_iter_pproc_mirr) ) call simple_rename(vol_iter_pproc_mirr, trim(REC_PPROC_MIRR_FBODY)//ext)
    !     ! updates original cls3D segment
    !     call work_proj%read_segment('ptcl3D', work_projfile)
    !     call work_proj%os_ptcl3D%delete_entry('stkind')
    !     call work_proj%os_ptcl3D%delete_entry('eo')
    !     params_glob%nptcls = ncavgs
    !     call spproj%os_cls3D%new(ncavgs, is_ptcl=.false.)
    !     do icls=1,ncavgs
    !         call spproj%os_cls3D%transfer_ori(icls, work_proj%os_ptcl3D, icls)
    !     enddo
    !     call conv_eo(work_proj%os_ptcl3D)
    !     ! revert splitting
    !     call spproj%os_cls3D%set_all2single('stkind', 1.)
    !     ! map the orientation parameters obtained for the clusters back to the particles
    !     select case(trim(params%imgkind))
    !     case('cavg')
    !         call spproj%map2ptcls
    !     case('cavg3D')
    !         ! TODO
    !     end select
    !     ! add rec_final to os_out
    !     call spproj%add_vol2os_out(trim(REC_FBODY)//ext, params%smpd, 1, 'vol_cavg')
    !     ! reprojections
    !     call spproj%os_cls3D%write('final_oris.txt')
    !     write(logfhandle,'(A)') '>>>'
    !     write(logfhandle,'(A)') '>>> RE-PROJECTION OF THE FINAL VOLUME'
    !     write(logfhandle,'(A)') '>>>'
    !     call cline_reproject%delete('lp_auto')
    !     call cline_reproject%set('vol1',   trim(REC_PPROC_FBODY)//ext)
    !     call cline_reproject%set('oritab', 'final_oris.txt')
    !     call xreproject%execute(cline_reproject)
    !     ! write alternated stack
    !     call img%new([params%box,params%box,1], params%smpd)
    !     call stkio_r%open(orig_stk,            params%smpd, 'read',                                 bufsz=500)
    !     call stkio_r2%open('reprojs.mrc',      params%smpd, 'read',                                 bufsz=500)
    !     call stkio_w%open('cavgs_reprojs.mrc', params%smpd, 'write', box=params%box, is_ft=.false., bufsz=500)
    !     cnt = -1
    !     do icls=1,ncavgs
    !         cnt = cnt + 2
    !         call stkio_r%read(icls, img)
    !         call img%norm
    !         call stkio_w%write(cnt, img)
    !         call stkio_r2%read(icls, img)
    !         call img%norm
    !         call stkio_w%write(cnt + 1, img)
    !     enddo
    !     call stkio_r%close
    !     call stkio_r2%close
    !     call stkio_w%close
    !     ! produce shifted stack
    !     call shift_imgfile(orig_stk, shifted_stk, spproj%os_cls3D, params%smpd)
    !     ! add shifted stack to project
    !     call spproj%add_cavgs2os_out(simple_abspath(shifted_stk), params%smpd, 'cavg_shifted')
    !     ! write results (this needs to be a full write as multiple segments are updated)
    !     call spproj%write()
    !     ! end gracefully
    !     call se1%kill
    !     call se2%kill
    !     call img%kill
    !     call spproj%kill
    !     call o%kill
    !     call o_even%kill
    !     call o_odd%kill
    !     call clsfrcs%kill
    !     call work_proj%kill
    !     call del_file(work_projfile)
    !     call simple_rmdir(STKPARTSDIR)
    !     call simple_end('**** SIMPLE_INITIAL_3DMODEL NORMAL STOP ****')

    !     contains
            
    !         subroutine downscale( smpd_target, scale_factor )
    !             real, intent(in)    :: smpd_target
    !             real, intent(inout) :: scale_factor
    !             do_autoscale     = .false.
    !             scale_factor     = 1.0
    !             params%box       = work_proj%get_box()
    !             params%smpd_crop = params%smpd
    !             params%box_crop  = params%box
    !             params%msk_crop  = params%msk
    !             if( params%l_autoscale )then
    !                 call autoscale(params%box, params%smpd, smpd_target, params%box_crop, params%smpd_crop, scale_factor, minbox=MINBOX)
    !                 do_autoscale = params%box_crop < params%box
    !                 if( do_autoscale )then
    !                     params%msk_crop = params%msk * scale_factor1
    !                     write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',params%box_crop
    !                     trslim = max(2.0, AHELIX_WIDTH / params%smpd_crop / 2.0)
    !                 else
    !                     trslim = max(2.0, AHELIX_WIDTH / params%smpd / 2.0)
    !                 endif
    !             endif
    !         end subroutine downscale

    !         function calc_lplim_stage2( nbest ) result( lplim )
    !             integer, intent(in)  :: nbest
    !             real,    allocatable :: res(:), tmp_rarr(:)
    !             integer, allocatable :: states(:), tmp_iarr(:)
    !             real :: lplim
    !             tmp_rarr  = spproj%os_cls2D%get_all('res')
    !             tmp_iarr  = nint(spproj%os_cls2D%get_all('state'))
    !             res       = pack(tmp_rarr, mask=(tmp_iarr>0))
    !             call hpsort(res)
    !             lplim = median_nocopy(res(:nbest))
    !             deallocate(tmp_rarr, tmp_iarr, res)
    !         end function calc_lplim_stage2

    !         subroutine rndstart( cline )
    !             class(cmdline), intent(inout) :: cline
    !             call work_proj%os_ptcl3D%rnd_oris
    !             call work_proj%os_ptcl3D%zero_shifts
    !             call work_proj%write_segment_inside('ptcl3D', work_projfile)
    !             call cline%set('mkdir', 'no') ! to avoid nested dirs
    !             call cline%set('objfun', 'cc')
    !             call xreconstruct3D%execute_shmem(cline)
    !             call cline%set('objfun', trim(params%objfun))
    !             call simple_copy_file('recvol_state01_even.mrc', 'startvol_even_unfil.mrc')
    !             call simple_copy_file('recvol_state01_odd.mrc',  'startvol_odd_unfil.mrc')
    !             call simple_rename(   'recvol_state01_even.mrc', 'startvol_even.mrc')
    !             call simple_rename(   'recvol_state01_odd.mrc',  'startvol_odd.mrc')
    !             call simple_rename(   'recvol_state01.mrc',      'startvol.mrc')
    !             call cline%set('vol1', 'startvol.mrc')
    !         end subroutine rndstart

    !         subroutine rec( cline )
    !             class(cmdline), intent(inout) :: cline
    !             type(parameters) :: params
    !             type(builder)    :: build
    !             call build%init_params_and_build_spproj(cline, params)
    !             call cline%set('mkdir', 'no') ! to avoid nested dirs
    !             call xreconstruct3D%execute_shmem(cline)
    !             call build%spproj_field%kill
    !             call simple_copy_file('recvol_state01_even.mrc', 'startvol_even_unfil.mrc')
    !             call simple_copy_file('recvol_state01_odd.mrc',  'startvol_odd_unfil.mrc')
    !             call simple_rename(   'recvol_state01_even.mrc', 'startvol_even.mrc')
    !             call simple_rename(   'recvol_state01_odd.mrc',  'startvol_odd.mrc')
    !             call simple_rename(   'recvol_state01.mrc',      'startvol.mrc')
    !             call cline%set('vol1', 'startvol.mrc')
    !         end subroutine rec

    !         subroutine conv_eo( os )
    !             class(oris), intent(inout) :: os
    !             type(sym) :: se
    !             type(ori) :: o_odd, o_even
    !             real      :: avg_euldist, euldist
    !             integer   :: icls, ncls
    !             call se%new(pgrp_refine)
    !             avg_euldist = 0.
    !             ncls = 0
    !             do icls=1,os%get_noris()/2
    !                 call os%get_ori(icls, o_even)
    !                 if( o_even%get_state() == 0 )cycle
    !                 ncls    = ncls + 1
    !                 call os%get_ori(ncavgs+icls, o_odd)
    !                 euldist = rad2deg(o_odd.euldist.o_even)
    !                 if( se%get_nsym() > 1 )then
    !                     call o_odd%mirror2d
    !                     call se%rot_to_asym(o_odd)
    !                     euldist = min(rad2deg(o_odd.euldist.o_even), euldist)
    !                 endif
    !                 avg_euldist = avg_euldist + euldist
    !             enddo
    !             avg_euldist = avg_euldist/real(ncls)
    !             write(logfhandle,'(A)')'>>>'
    !             write(logfhandle,'(A,F6.1)')'>>> EVEN/ODD AVERAGE ANGULAR DISTANCE: ', avg_euldist
    !         end subroutine conv_eo

    ! end subroutine exec_initial_3Dmodel

    !> for generation of an initial 3d model from class averages
    subroutine exec_initial_3Dmodel( self, cline )
        class(initial_3Dmodel_commander), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        ! constants
        real,                  parameter :: LPSTART_LB=10., LPSTART_DEFAULT=20., LPSTOP_LB=6.
        real,                  parameter :: CENLP_DEFAULT   = 30.
        real,                  parameter :: LPSYMSRCH_LB    = 12.
        integer,               parameter :: NSTAGES         = 8
        integer,               parameter :: PHASES(3)       = [2,6,8]
        integer,               parameter :: MAXITS(3)       = [20,15,10]
        integer,               parameter :: MAXITS_GLOB(3)  = [2*20,4*15,2*10]
        integer,               parameter :: NSPACE(3)       = [500,1000,2500]
        integer,               parameter :: SYMSRCH_STAGE   = 3
        character(len=STDLEN), parameter :: work_projfile   = 'initial_3Dmodel_tmpproj.simple'
        ! distributed commanders
        type(calc_pspec_commander_distr) :: xcalc_pspec_distr
        ! shared-mem commanders
        type(refine3D_commander)         :: xrefine3D
        type(reconstruct3D_commander)    :: xreconstruct3D
        type(symaxis_search_commander)   :: xsymsrch
        type(reproject_commander)        :: xreproject
        type(postprocess_commander)      :: xpostprocess
        ! command lines
        type(cmdline) :: cline_refine3D
        type(cmdline) :: cline_symsrch
        type(cmdline) :: cline_reconstruct3D, cline_postprocess
        type(cmdline) :: cline_reproject, cline_calc_pspec
        ! other
        character(len=:),  allocatable :: stk, stkpath, orig_stk, frcs_fname, shifted_stk, stk_even, stk_odd, ext
        integer,           allocatable :: states(:)
        real,              allocatable :: frcs_avg(:)
        type(lp_crop_inf), allocatable :: lpinfo(:)
        character(len=2)      :: str_state
        type(ori)             :: o, o_even, o_odd
        type(qsys_env)        :: qenv
        type(parameters)      :: params
        type(ctfparams)       :: ctfvars
        type(sp_project)      :: spproj, work_proj
        type(sym)             :: se1,se2
        type(image)           :: img, vol
        type(stack_io)        :: stkio_r, stkio_r2, stkio_w
        type(class_frcs)      :: clsfrcs
        character(len=STDLEN) :: vol_iter, pgrp_init, pgrp_refine, vol_iter_pproc, vol_iter_pproc_mirr
        character(len=STDLEN) :: sigma2_fname, sigma2_fname_sc, orig_objfun, frckind
        real                  :: lpsym, lpfinal
        integer               :: icls, ncavgs, cnt, iter, ipart, even_ind, odd_ind, state, filtsz, istage
        logical               :: srch4symaxis, symran_before_refine
        call cline%set('objfun',    'euclid') ! use noise normalized Euclidean distances from the start
        call cline%set('sigma_est', 'global') ! obviously
        call cline%set('oritype',      'out') ! because cavgs are part of out segment
        call cline%set('bfac',            0.) ! because initial models should not be sharpened
        if( .not. cline%defined('mkdir')       ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('overlap')     ) call cline%set('overlap',     0.95) ! needed to prevent premature convergence
        if( .not. cline%defined('prob_athres') ) call cline%set('prob_athres',  90.) ! reduces # failed runs on trpv1 from 4->2/10
        if( .not. cline%defined('cenlp')       ) call cline%set('cenlp', CENLP_DEFAULT)
        if( .not. cline%defined('imgkind')     ) call cline%set('imgkind',   'cavg') ! whether to use classes generated from 2D/3D
        ! make master parameters
        call params%new(cline)
        call cline%set('mkdir',       'no')   ! to avoid nested directory structure
        call cline%set('oritype', 'ptcl3D')   ! from now on we are in the ptcl3D segment, final report is in the cls3D segment
        ! state string
        str_state    = int2str_pad(1,2)
        ! decide wether to search for the symmetry axis
        pgrp_init    = trim(params%pgrp_start)
        pgrp_refine  = trim(params%pgrp)
        srch4symaxis = trim(pgrp_refine) .ne. trim(pgrp_init)
        symran_before_refine = .false.
        if( pgrp_init.ne.'c1' .or. pgrp_refine.ne.'c1' )then
            se1 = sym(pgrp_init)
            se2 = sym(pgrp_refine)
            if(se1%get_nsym() > se2%get_nsym())then
                ! ensure se2 is a subgroup of se1
                if( .not. se1%has_subgrp(pgrp_refine) ) THROW_HARD('Incompatible symmetry groups; simple_commander_abinitio2')
                ! set flag for symmetry randomisation before refinmement
                ! in case we are moving from a higher to lower group
                symran_before_refine = .true.
            else if( se2%get_nsym() > se1%get_nsym() )then
                ! ensure se1 is a subgroup of se2
                if( .not. se2%has_subgrp(pgrp_init) ) THROW_HARD('Incompatible symmetry groups; simple_commander_abinitio2')
            endif
        endif
        ! read project & update sampling distance
        call spproj%read(params%projfile)
        call spproj%update_projinfo(cline)
        call spproj%write_segment_inside('projinfo', params%projfile)
        ! whether to use classes generated from 2D or 3D
        select case(trim(params%imgkind))
            case('cavg')
                frckind = 'frc2D'
                states  = nint(spproj%os_cls2D%get_all('state'))
            case('cavg3D')
                frckind = 'frc3D'
                states  = nint(spproj%os_cls3D%get_all('state'))
            case DEFAULT
                THROW_HARD('Unsupported IMGKIND!')
        end select
        call cline%delete('imgkind') ! no interference down the line
        ! retrieve cavgs stack info
        call spproj%get_cavgs_stk(stk, ncavgs, params%smpd, imgkind=params%imgkind, stkpath=stkpath)
        if(.not. file_exists(stk)) stk = trim(stkpath) // '/' // trim(stk)
        if(.not. file_exists(stk)) THROW_HARD('cavgs stk does not exist; simple_commander_abinitio')
        orig_stk        = stk
        ext             = '.'//fname2ext(stk)
        stk_even        = add2fbody(trim(stk), trim(ext), '_even')
        stk_odd         = add2fbody(trim(stk), trim(ext), '_odd')
        if( .not. file_exists(stk_even) ) THROW_HARD('Even cavgs stk: '//trim(stk_even)//' does not exist!')
        if( .not. file_exists(stk_odd)  ) THROW_HARD('Odd cavgs stk: '//trim(stk_odd)//' does not exist!')
        ctfvars%ctfflag = CTFFLAG_NO
        ctfvars%smpd    = params%smpd
        shifted_stk     = basename(add2fbody(stk, ext, '_shifted'))
        if( count(states==0) .eq. ncavgs )then
            THROW_HARD('no class averages detected in project file: '//trim(params%projfile)//'; initial_3Dmodel')
        endif
        ! retrieve FRC info
        call spproj%get_frcs(frcs_fname, frckind, fail=.false.)
        if( .not.file_exists(frcs_fname) )then
            ! 08/24 This is a backwards compatibility patch to account for error in metadata
            ! on exit of streaming related to GUI directory structure (now fixed and cf above get_cavgs_stk).
            ! Will need to harmonize (move to absolute path?).
            frcs_fname = trim(stkpath)//'/'//trim(frcs_fname)
            if( .not.file_exists(frcs_fname) )then
                THROW_HARD('the project file does not contain an FRCs file, which is required')
            endif
        endif
        ! work out low-pass limits and downscaling parameters
        params%frcs = trim(frcs_fname)
        call clsfrcs%read(frcs_fname)
        filtsz = clsfrcs%get_filtsz()
        allocate(frcs_avg(filtsz), source=0.)
        call clsfrcs%avg_frc_getter(frcs_avg, states)
        allocate(lpinfo(NSTAGES))
        lpfinal = max(LPSTOP_LB,calc_lplim_final_stage(3))
        call lpstages(params%box, NSTAGES, frcs_avg, params%smpd, LPSTART_LB, LPSTART_DEFAULT, lpfinal, lpinfo, verbose=.true.)
        ! prepare a temporary project file
        work_proj%projinfo = spproj%projinfo
        work_proj%compenv  = spproj%compenv
        if( spproj%jobproc%get_noris()  > 0 ) work_proj%jobproc = spproj%jobproc
        ! name change
        call work_proj%projinfo%delete_entry('projname')
        call work_proj%projinfo%delete_entry('projfile')
        call cline%set('projfile', trim(work_projfile))
        call cline%set('projname', trim(get_fbody(trim(work_projfile),trim('simple'))))
        call work_proj%update_projinfo(cline)
        ! add stks to temporary project
        call work_proj%add_stk(stk_even, ctfvars)
        call work_proj%add_stk(stk_odd,  ctfvars)
        ! update orientations parameters
        do icls=1,ncavgs
            even_ind = icls
            odd_ind  = ncavgs + icls
            call work_proj%os_ptcl3D%get_ori(icls, o)
            call o%set('class', real(icls))
            call o%set('state', real(states(icls)))
            ! even
            o_even = o
            call o_even%set('eo', 0.)
            call o_even%set('stkind', work_proj%os_ptcl3D%get(even_ind,'stkind'))
            call work_proj%os_ptcl3D%set_ori(even_ind, o_even)
            ! odd
            o_odd = o
            call o_odd%set('eo', 1.)
            call o_odd%set('stkind', work_proj%os_ptcl3D%get(odd_ind,'stkind'))
            call work_proj%os_ptcl3D%set_ori(odd_ind, o_odd)
        enddo
        params_glob%nptcls = work_proj%get_nptcls()
        call work_proj%write()
        ! prepare command lines from prototype
        cline_reconstruct3D = cline
        cline_refine3D      = cline
        cline_reproject     = cline
        cline_symsrch       = cline
        ! symmetry axis search
        if( srch4symaxis )then
            call qenv%new(1, exec_bin='simple_exec')
            call cline_symsrch%set('prg',         'symaxis_search') ! needed for cluster exec
            call cline_symsrch%set('pgrp',       trim(pgrp_refine))
            call cline_symsrch%set('smpd',        params%smpd_crop)
            call cline_symsrch%set('box',    real(params%box_crop))
            call cline_symsrch%set('projfile', trim(work_projfile))
            if( .not. cline_symsrch%defined('cenlp') ) call cline_symsrch%set('cenlp', CENLP_DEFAULT)
            call cline_symsrch%set('hp',                 params%hp)
            call cline_symsrch%set('oritype',             'ptcl3D')
        endif
        ! re-reconstruct & re-project volume
        call cline_reconstruct3D%set('prg',       'reconstruct3D')
        call cline_reconstruct3D%set('box',      real(params%box))
        call cline_reconstruct3D%set('projfile',    work_projfile)
        call cline_reconstruct3D%set('needs_sigma',         'yes')
        call cline_postprocess%set('prg',           'postprocess')
        call cline_postprocess%set('projfile',      work_projfile)
        call cline_postprocess%set('mkdir',                  'no')
        call cline_postprocess%delete('bfac') ! sharpen final map
        call cline_reproject%set('prg',               'reproject')
        call cline_reproject%set('pgrp',        trim(pgrp_refine))
        call cline_reproject%set('outstk',         'reprojs'//ext)
        call cline_reproject%set('smpd',              params%smpd)
        call cline_reproject%set('box',          real(params%box))
        ! Frequency marching
        do istage = 1, NSTAGES
            write(logfhandle,'(A)')'>>>'
            write(logfhandle,'(A,I3,A9,F5.1)')'>>> STAGE ', istage,' WITH LP =', lpinfo(istage)%lp
            if( istage == 1 )then
                call rndstart(cline_refine3D)
            else
                call rec(cline_refine3D)
            endif
            call set_cline_refine3D(istage)
            if( lpinfo(istage)%l_autoscale )then
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',lpinfo(istage)%box_crop
            endif
            call xrefine3D%execute_shmem(cline_refine3D)
            call work_proj%read_segment('ptcl3D', work_projfile)
            if( istage == SYMSRCH_STAGE )then
                if( symran_before_refine )then
                    call se1%symrandomize(work_proj%os_ptcl3D)
                    call work_proj%write_segment_inside('ptcl3D', work_projfile)
                endif
                if( srch4symaxis )then
                    vol_iter = trim(VOL_FBODY)//trim(str_state)//ext
                    if( .not. file_exists(vol_iter) ) THROW_HARD('input volume to symmetry axis search does not exist')
                    lpsym = max(LPSYMSRCH_LB,lpinfo(SYMSRCH_STAGE)%lp)
                    write(logfhandle,'(A,F5.1)') '>>> DID SET SYMSEARCH LOW-PASS LIMIT (IN A) TO: ', lpsym
                    call cline_symsrch%set('lp', lpsym)
                    write(logfhandle,'(A)') '>>>'
                    write(logfhandle,'(A)') '>>> SYMMETRY AXIS SEARCH'
                    write(logfhandle,'(A)') '>>>'
                    call cline_symsrch%set('vol1', trim(vol_iter))
                    call xsymsrch%execute_shmem(cline_symsrch)
                    call del_file('SYMAXIS_SEARCH_FINISHED')
                endif
            endif
        end do
        ! sigma2 at original sampling
        cline_calc_pspec = cline
        call cline_calc_pspec%set('prg',      'calc_pspec' )
        call cline_calc_pspec%set('projfile', work_projfile)
        call cline_calc_pspec%set('box',   real(params%box))
        call cline_calc_pspec%set('smpd',       params%smpd)
        call cline_calc_pspec%set('which_iter',  real(iter))
        call xcalc_pspec_distr%execute_shmem(cline_calc_pspec)
        call rec(cline_refine3D)
        iter = nint(cline_refine3D%get_rarg('endit'))
        call cline_reconstruct3D%set('which_iter',real(iter))
        ! deal with final volume
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> RECONSTRUCTION AT ORIGINAL SAMPLING'
        write(logfhandle,'(A)') '>>>'
        call cline_reconstruct3D%delete('lp_auto')
        call cline_reconstruct3D%set('box',  real(params%box))
        call cline_reconstruct3D%set('smpd', params%smpd)
        ! reconstruction
        call xreconstruct3D%execute_shmem(cline_reconstruct3D)
        vol_iter = trim(VOL_FBODY)//trim(str_state)//ext
        ! because postprocess only updates project file when mkdir=yes
        call work_proj%read_segment('out', work_projfile)
        call work_proj%add_vol2os_out(vol_iter, params%smpd, 1, 'vol')
        call work_proj%add_fsc2os_out(FSC_FBODY//str_state//trim(BIN_EXT), 1, params%box)
        call work_proj%write_segment_inside('out', work_projfile)
        call xpostprocess%execute(cline_postprocess)
        vol_iter_pproc      = add2fbody(vol_iter,ext,PPROC_SUFFIX)
        vol_iter_pproc_mirr = add2fbody(vol_iter,ext,trim(PPROC_SUFFIX)//trim(MIRR_SUFFIX))
        if( file_exists(vol_iter)            ) call simple_rename(vol_iter,            trim(REC_FBODY)//ext)
        if( file_exists(vol_iter_pproc)      ) call simple_rename(vol_iter_pproc,      trim(REC_PPROC_FBODY)//ext)
        if( file_exists(vol_iter_pproc_mirr) ) call simple_rename(vol_iter_pproc_mirr, trim(REC_PPROC_MIRR_FBODY)//ext)
        ! updates original cls3D segment
        call work_proj%read_segment('ptcl3D', work_projfile)
        call work_proj%os_ptcl3D%delete_entry('stkind')
        call work_proj%os_ptcl3D%delete_entry('eo')
        params_glob%nptcls = ncavgs
        call spproj%os_cls3D%new(ncavgs, is_ptcl=.false.)
        do icls=1,ncavgs
            call spproj%os_cls3D%transfer_ori(icls, work_proj%os_ptcl3D, icls)
        enddo
        call conv_eo(work_proj%os_ptcl3D)
        ! revert splitting
        call spproj%os_cls3D%set_all2single('stkind', 1.)
        ! map the orientation parameters obtained for the clusters back to the particles
        select case(trim(params%imgkind))
        case('cavg')
            call spproj%map2ptcls
        case('cavg3D')
            ! TODO
        end select
        ! add rec_final to os_out
        call spproj%add_vol2os_out(trim(REC_FBODY)//ext, params%smpd, 1, 'vol_cavg')
        ! reprojections
        call spproj%os_cls3D%write('final_oris.txt')
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> RE-PROJECTION OF THE FINAL VOLUME'
        write(logfhandle,'(A)') '>>>'
        call cline_reproject%delete('lp_auto')
        call cline_reproject%set('vol1',   trim(REC_PPROC_FBODY)//ext)
        call cline_reproject%set('oritab', 'final_oris.txt')
        call xreproject%execute(cline_reproject)
        ! write alternated stack
        call img%new([params%box,params%box,1], params%smpd)
        call stkio_r%open(orig_stk,            params%smpd, 'read',                                 bufsz=500)
        call stkio_r2%open('reprojs.mrc',      params%smpd, 'read',                                 bufsz=500)
        call stkio_w%open('cavgs_reprojs.mrc', params%smpd, 'write', box=params%box, is_ft=.false., bufsz=500)
        cnt = -1
        do icls=1,ncavgs
            cnt = cnt + 2
            call stkio_r%read(icls, img)
            call img%norm
            call stkio_w%write(cnt, img)
            call stkio_r2%read(icls, img)
            call img%norm
            call stkio_w%write(cnt + 1, img)
        enddo
        call stkio_r%close
        call stkio_r2%close
        call stkio_w%close
        ! produce shifted stack
        call shift_imgfile(orig_stk, shifted_stk, spproj%os_cls3D, params%smpd)
        ! add shifted stack to project
        call spproj%add_cavgs2os_out(simple_abspath(shifted_stk), params%smpd, 'cavg_shifted')
        ! write results (this needs to be a full write as multiple segments are updated)
        call spproj%write()
        ! end gracefully
        call se1%kill
        call se2%kill
        call img%kill
        call spproj%kill
        call o%kill
        call o_even%kill
        call o_odd%kill
        call clsfrcs%kill
        call work_proj%kill
        call del_file(work_projfile)
        call simple_rmdir(STKPARTSDIR)
        call simple_end('**** SIMPLE_INITIAL_3DMODEL NORMAL STOP ****')

        contains

           function calc_lplim_final_stage( nbest ) result( lplim )
                integer, intent(in)  :: nbest
                real,    allocatable :: res(:), tmp_rarr(:)
                integer, allocatable :: states(:), tmp_iarr(:)
                real :: lplim
                tmp_rarr  = spproj%os_cls2D%get_all('res')
                tmp_iarr  = nint(spproj%os_cls2D%get_all('state'))
                res       = pack(tmp_rarr, mask=(tmp_iarr>0))
                call hpsort(res)
                lplim = median_nocopy(res(:nbest))
                deallocate(tmp_rarr, tmp_iarr, res)
            end function calc_lplim_final_stage

            subroutine rndstart( cline )
                class(cmdline), intent(inout) :: cline
                call work_proj%os_ptcl3D%rnd_oris
                call work_proj%os_ptcl3D%zero_shifts
                call work_proj%write_segment_inside('ptcl3D', work_projfile)
                call cline%set('mkdir', 'no') ! to avoid nested dirs
                call cline%set('objfun', 'cc')
                call cline%set('silence_fsc', 'yes')
                call xreconstruct3D%execute_shmem(cline)
                call cline%set('objfun', trim(params%objfun))
                call simple_copy_file('recvol_state01_even.mrc', 'startvol_even_unfil.mrc')
                call simple_copy_file('recvol_state01_odd.mrc',  'startvol_odd_unfil.mrc')
                call simple_rename(   'recvol_state01_even.mrc', 'startvol_even.mrc')
                call simple_rename(   'recvol_state01_odd.mrc',  'startvol_odd.mrc')
                call simple_rename(   'recvol_state01.mrc',      'startvol.mrc')
                call cline%set('vol1', 'startvol.mrc')
            end subroutine rndstart

            subroutine rec( cline )
                class(cmdline), intent(inout) :: cline
                type(parameters) :: params
                type(builder)    :: build
                call build%init_params_and_build_spproj(cline, params)
                call cline%set('mkdir', 'no') ! to avoid nested dirs
                call xreconstruct3D%execute_shmem(cline)
                call build%spproj_field%kill
                call simple_copy_file('recvol_state01_even.mrc', 'startvol_even_unfil.mrc')
                call simple_copy_file('recvol_state01_odd.mrc',  'startvol_odd_unfil.mrc')
                call simple_rename(   'recvol_state01_even.mrc', 'startvol_even.mrc')
                call simple_rename(   'recvol_state01_odd.mrc',  'startvol_odd.mrc')
                call simple_rename(   'recvol_state01.mrc',      'startvol.mrc')
                call cline%set('vol1', 'startvol.mrc')
            end subroutine rec

            subroutine set_cline_refine3D( istage )
                integer, intent(in) :: istage
                character(len=:), allocatable :: silence_fsc, sh_first, prob_sh, ml_reg, refine, icm
                integer :: iphase
                real    :: trs, rnspace, rmaxits, rmaxits_glob, riter, snr_noise_reg
                ! iteration number bookkeeping
                if( cline_refine3D%defined('endit') )then
                    riter = cline_refine3D%get_rarg('endit')
                else
                    riter = 0.
                endif
                riter = riter + 1.0
                ! phase logics
                if(      istage <= PHASES(1) )then
                    iphase = 1
                else if( istage <= PHASES(2) )then
                    iphase = 2
                else if( istage <= PHASES(3) )then
                    iphase = 3
                else 
                    THROW_HARD('Invalid istage index')
                endif
                ! phase control parameters
                select case(iphase)
                    case(1)
                        refine        = 'shc_smpl'
                        rnspace       = real(NSPACE(1))
                        rmaxits       = real(MAXITS(1))
                        rmaxits_glob  = real(MAXITS_GLOB(1))
                        silence_fsc   = 'yes'
                        trs           = 0.
                        snr_noise_reg = 2.0
                        sh_first      = 'no'
                        prob_sh       = 'no'
                        ml_reg        = 'no'
                        icm           = 'no'
                    case(2)
                        refine        = 'shc_smpl'
                        rnspace       = real(NSPACE(2))
                        rmaxits       = real(MAXITS(2))
                        rmaxits_glob  = real(MAXITS_GLOB(2))
                        silence_fsc   = 'yes'
                        trs           = lpinfo(istage)%trslim
                        snr_noise_reg = 4.0
                        sh_first      = 'yes'
                        prob_sh       = 'no'
                        ml_reg        = 'yes'
                        icm           = 'no'
                    case(3)
                        refine        = 'prob'
                        rnspace       = real(NSPACE(3))
                        rmaxits       = real(MAXITS(3))
                        rmaxits_glob  = real(MAXITS_GLOB(3))
                        silence_fsc   = 'no'
                        trs           = lpinfo(istage)%trslim
                        snr_noise_reg = 6.0
                        sh_first      = 'yes'
                        prob_sh       = 'yes'
                        ml_reg        = 'yes'
                        icm           = 'yes'
                end select
                ! symmetry
                if( srch4symaxis )then
                    if( istage <= SYMSRCH_STAGE )then
                        ! need to replace original point-group flag with c1/pgrp_start
                        call cline_refine3D%set('pgrp', trim(pgrp_init))
                    else
                        call cline_refine3D%set('pgrp', trim(pgrp_refine))
                    endif
                endif
                ! command line update
                call cline_refine3D%set('prg',                     'refine3D')
                call cline_refine3D%set('startit',                      riter)
                call cline_refine3D%set('which_iter',                   riter)
                call cline_refine3D%set('refine',                      refine)
                call cline_refine3D%set('lp',               lpinfo(istage)%lp)
                call cline_refine3D%set('smpd_crop', lpinfo(istage)%smpd_crop)
                call cline_refine3D%set('box_crop',   lpinfo(istage)%box_crop)
                call cline_refine3D%set('nspace',                     rnspace)
                call cline_refine3D%set('maxits',                     rmaxits)
                call cline_refine3D%set('maxits_glob',           rmaxits_glob)
                call cline_refine3D%set('silence_fsc',            silence_fsc)
                call cline_refine3D%set('trs',                            trs)
                call cline_refine3D%set('snr_noise_reg',        snr_noise_reg)
                call cline_refine3D%set('sh_first',                  sh_first)
                call cline_refine3D%set('prob_sh',                    prob_sh)
                call cline_refine3D%set('ml_reg',                      ml_reg)
                call cline_refine3D%set('icm',                            icm)

                call cline_refine3D%printline

            end subroutine set_cline_refine3D
    
            subroutine conv_eo( os )
                class(oris), intent(inout) :: os
                type(sym) :: se
                type(ori) :: o_odd, o_even
                real      :: avg_euldist, euldist
                integer   :: icls, ncls
                call se%new(pgrp_refine)
                avg_euldist = 0.
                ncls = 0
                do icls=1,os%get_noris()/2
                    call os%get_ori(icls, o_even)
                    if( o_even%get_state() == 0 )cycle
                    ncls    = ncls + 1
                    call os%get_ori(ncavgs+icls, o_odd)
                    euldist = rad2deg(o_odd.euldist.o_even)
                    if( se%get_nsym() > 1 )then
                        call o_odd%mirror2d
                        call se%rot_to_asym(o_odd)
                        euldist = min(rad2deg(o_odd.euldist.o_even), euldist)
                    endif
                    avg_euldist = avg_euldist + euldist
                enddo
                avg_euldist = avg_euldist/real(ncls)
                write(logfhandle,'(A)')'>>>'
                write(logfhandle,'(A,F6.1)')'>>> EVEN/ODD AVERAGE ANGULAR DISTANCE: ', avg_euldist
            end subroutine conv_eo

    end subroutine exec_initial_3Dmodel

    !> for generation of an initial 3d model from particles
    subroutine exec_abinitio_3Dmodel_autolp( self, cline )
        class(abinitio_3Dmodel_autolp_commander), intent(inout) :: self
        class(cmdline),                           intent(inout) :: cline
        real,    parameter :: CENLP_DEFAULT   = 30.
        real,    parameter :: STARTLP_DEFAULT = 20.
        real,    parameter :: LP_SYMSRCH_LB   = 12.
        integer, parameter :: MAXITS1         = 20
        integer, parameter :: MAXITS2         = 10
        integer, parameter :: MAXITS3         = 30
        integer, parameter :: NSPACE1         = 500
        integer, parameter :: NSPACE2         = 1000
        integer, parameter :: NSPACE3         = 2500
        integer, parameter :: MINBOX          = 88
        ! commanders
        type(refine3D_commander_distr)      :: xrefine3D_distr
        type(reconstruct3D_commander_distr) :: xreconstruct3D_distr
        type(symaxis_search_commander)      :: xsymsrch
        type(postprocess_commander)         :: xpostprocess
        ! command lines
        type(cmdline) :: cline_refine3D_1, cline_refine3D_2, cline_refine3D_3
        type(cmdline) :: cline_reconstruct3D, cline_reconstruct3D_mlreg
        type(cmdline) :: cline_postprocess, cline_symsrch
        ! other
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(sym)                     :: se1, se2
        type(class_frcs)              :: clsfrcs
        type(image)                   :: final_vol, reprojs
        character(len=:), allocatable :: vol_type, str_state, vol, vol_pproc, vol_pproc_mirr, frcs_fname, vol_iter
        character(len=LONGSTRLEN)     :: vol_str
        integer,          allocatable :: states_cavg(:)
        real,             allocatable :: frcs_avg(:)
        integer :: find, state, filtsz, find_start, iter
        real    :: smpd_target, scale_factor1, scale_factor2, trslim, lplims(2), lp, lp_est, lp_sym
        logical :: l_autoscale, l_err, l_srch4symaxis, l_symran, l_sym
        call cline%set('objfun',    'euclid') ! use noise normalized Euclidean distances from the start
        call cline%set('sigma_est', 'global') ! obviously
        call cline%set('oritype',   'ptcl3D') ! obviously
        if( .not. cline%defined('mkdir')        ) call cline%set('mkdir',        'yes')
        if( .not. cline%defined('autoscale')    ) call cline%set('autoscale',    'yes')
        if( .not. cline%defined('center')       ) call cline%set('center',        'no')
        if( .not. cline%defined('pgrp')         ) call cline%set('pgrp',          'c1')
        if( .not. cline%defined('pgrp_start')   ) call cline%set('pgrp_start',    'c1')
        ! make master parameters
        call params%new(cline)
        call cline%delete('autoscale')
        call cline%set('mkdir', 'no') ! to avoid nested directory structure
        ! state string
        str_state      = int2str_pad(1,2)
        ! decide wether to search for the symmetry axis
        l_srch4symaxis = trim(params%pgrp) .ne. trim(params%pgrp_start)
        l_symran       = .false.
        l_sym          = l_srch4symaxis
        if( params%pgrp_start.ne.'c1' .or. params%pgrp.ne.'c1' )then
            se1 = sym(params%pgrp_start)
            se2 = sym(params%pgrp)
            if(se1%get_nsym() > se2%get_nsym())then
                ! ensure se2 is a subgroup of se1
                if( .not. se1%has_subgrp(params%pgrp) )THROW_HARD('Incompatible symmetry groups; exec_abinitio_3Dmodel')
                ! set flag for symmetry randomisation
                ! in case we are moving from a higher to lower group
                l_symran = .true.
            else if( se2%get_nsym() > se1%get_nsym() )then
                ! ensure se1 is a subgroup of se2
                if( .not. se2%has_subgrp(params%pgrp_start) )THROW_HARD('Incompatible symmetry groups; exec_abinitio_3Dmodel')
            endif
        endif
        ! read project
        call spproj%read(params%projfile)
        call spproj%update_projinfo(cline)
        call spproj%write_segment_inside('projinfo', params%projfile)
        ! retrieve FRC info
        states_cavg = nint(spproj%os_cls2D%get_all('state'))
        if( all(states_cavg == 0) )then
            THROW_HARD('no class averages detected in project file: '//trim(params%projfile))
        endif
        call spproj%get_frcs(frcs_fname, 'frc2D', fail=.false.)
        if( .not.file_exists(frcs_fname) )then
            THROW_HARD('the project file does not contain an FRCs file, which is required')
        endif
        params%frcs = trim(frcs_fname)
        call clsfrcs%read(frcs_fname)
        filtsz = clsfrcs%get_filtsz()
        allocate(frcs_avg(filtsz), source=0.)
        call clsfrcs%avg_frc_getter(frcs_avg, states_cavg)
        find = get_find_at_corr(frcs_avg, 0.9)
        lp   = calc_lowpass_lim(find, params%box, params%smpd)
        print *, 'lp@0.9   ', lp
        find = get_find_at_corr(frcs_avg, 0.8)
        lp   = calc_lowpass_lim(find, params%box, params%smpd)
        print *, 'lp@0.8   ', lp
        find = get_find_at_corr(frcs_avg, 0.7)
        lp   = calc_lowpass_lim(find, params%box, params%smpd)
        print *, 'lp@0.7   ', lp
        find = get_find_at_corr(frcs_avg, 0.6)
        lp   = calc_lowpass_lim(find, params%box, params%smpd)
        print *, 'lp@0.6   ', lp
        find = get_find_at_corr(frcs_avg, 0.5)
        lp   = calc_lowpass_lim(find, params%box, params%smpd)
        print *, 'lp@0.5   ', lp
        find = get_find_at_corr(frcs_avg, 0.4)
        lp   = calc_lowpass_lim(find, params%box, params%smpd)
        print *, 'lp@0.4   ', lp
        find = get_find_at_corr(frcs_avg, 0.3)
        lp   = calc_lowpass_lim(find, params%box, params%smpd)
        print *, 'lp@0.3   ', lp
        find = get_find_at_corr(frcs_avg, 0.143)
        lp   = calc_lowpass_lim(find, params%box, params%smpd)
        print *, 'lp@0.143 ', lp
        ! set lplims
        if( cline%defined('lpstart') )then
            lplims(1) = params%lpstart
        else if( any(frcs_avg > 0.8) )then
            find      = get_find_at_corr(frcs_avg, 0.8)
            lplims(1) = calc_lowpass_lim(find, params%box, params%smpd)
        else
            lplims(1) = STARTLP_DEFAULT
        endif
        if( any(frcs_avg > 0.5) )then
            find      = get_find_at_corr(frcs_avg, 0.5)
            lplims(2) = calc_lowpass_lim(find, params%box, params%smpd)
        else
            lplims(2) = STARTLP_DEFAULT
        endif
        if( .not. cline%defined('vol1') )then
            ! randomize projection directions
            select case(trim(params%oritype))
                case('ptcl3D')
                    if( spproj%os_ptcl3D%get_noris() < 1 )then
                        THROW_HARD('Particles could not be found in the project')
                    endif
                    vol_type = 'vol'
                    call spproj%os_ptcl3D%rnd_oris
                    call spproj%os_ptcl3D%set_all2single('w',1.)
                case DEFAULT
                    THROW_HARD('Unsupported ORITYPE')
            end select
            call spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        write(logfhandle,'(A,F5.1)') '>>> DID SET STARTING  LOW-PASS LIMIT (IN A) TO: ', lplims(1)
        write(logfhandle,'(A,F5.1)') '>>> DID SET HARD      LOW-PASS LIMIT (IN A) TO: ', lplims(2)
        write(logfhandle,'(A,F5.1)') '>>> DID SET CENTERING LOW-PASS LIMIT (IN A) TO: ', params_glob%cenlp
        ! initial phase scaling
        smpd_target = max(params%smpd, lplims(2)*LP2SMPDFAC)
        call downscale(smpd_target, scale_factor1)
        ! command-lines
        cline_refine3D_1          = cline
        cline_refine3D_2          = cline
        cline_refine3D_3          = cline
        cline_reconstruct3D       = cline
        cline_reconstruct3D_mlreg = cline
        cline_postprocess         = cline
        cline_symsrch             = cline
        ! initialise command line parameters
        ! (1)
        call cline_refine3D_1%set('prg',                 'refine3D')
        call cline_refine3D_1%set('projfile', trim(params%projfile))
        call cline_refine3D_1%set('box_crop', real(params%box_crop))
        call cline_refine3D_1%set('smpd_crop',     params%smpd_crop)
        call cline_refine3D_1%set('lpstart',              lplims(1))
        call cline_refine3D_1%set('lpstop',               lplims(2))
        call cline_refine3D_1%set('nspace',           real(NSPACE1))
        call cline_refine3D_1%set('maxits',           real(MAXITS1))
        call cline_refine3D_1%set('pgrp',         params%pgrp_start)
        call cline_refine3D_1%set('silence_fsc',              'no')
        call cline_refine3D_1%set('refine',              'shc_smpl') ! this works for class averages but may not work on particles
        call cline_refine3D_1%set('lp_auto',                  'yes')
        call cline_refine3D_1%set('sh_first',                  'no')
        call cline_refine3D_1%set('prob_sh',                   'no')
        call cline_refine3D_1%set('ml_reg',                    'no') ! No ML-reg in the start, we don't trust the FSC yet
        call cline_refine3D_1%set('snr_noise_reg',              4.0)
        call cline_refine3D_1%set('trs',                        0.0) ! no shifts in phase 1
        ! (2)
        call cline_refine3D_2%set('prg',                 'refine3D')
        call cline_refine3D_2%set('projfile', trim(params%projfile))
        call cline_refine3D_2%set('box_crop', real(params%box_crop))
        call cline_refine3D_2%set('smpd_crop',     params%smpd_crop)
        call cline_refine3D_2%set('nspace',           real(NSPACE2)) ! # projection directions are increased
        call cline_refine3D_2%set('maxits',           real(MAXITS2))
        call cline_refine3D_2%set('pgrp',         params%pgrp_start)
        call cline_refine3D_2%set('silence_fsc',               'no')
        call cline_refine3D_2%set('refine',                  'prob') ! changing to prob refinement
        call cline_refine3D_2%set('prob_athres',                90.) ! following the same logic as in initial_3Dmodel, phase 2
        call cline_refine3D_2%set('lp_auto',                  'yes')
        call cline_refine3D_2%set('sh_first',                 'yes') ! first shift logic is turned on
        call cline_refine3D_2%set('prob_sh',                  'yes') ! probabilistic shift search is turned on
        call cline_refine3D_2%set('ml_reg',                    'no') ! No ML-reg in the start, we don't trust the FSC yet
        call cline_refine3D_2%set('trs',                     trslim) ! shifts are turned on, trslim set by downscale above
        ! (3)
        call cline_refine3D_3%set('prg',                 'refine3D')
        call cline_refine3D_3%set('projfile', trim(params%projfile))
        call cline_refine3D_3%set('nspace',           real(NSPACE3)) ! # projection directions are increased
        call cline_refine3D_3%set('maxits',           real(MAXITS3))
        call cline_refine3D_3%set('pgrp',               params%pgrp) ! highest point-group symmetry is used
        call cline_refine3D_3%set('silence_fsc',               'no')
        if( cline%defined('refine') )then
            call cline_refine3D_3%set('refine', trim(params%refine)) ! possible to input desired refinement mode for final phase
        else
            call cline_refine3D_3%set('refine',              'prob')
        endif
        call cline_refine3D_2%set('prob_athres',                10.) ! following what have worked in the past on export gate
        call cline_refine3D_3%set('refine',                  'prob')
        call cline_refine3D_3%set('lp_auto',                  'fsc') ! changing to FSC-based lp estimation 
        call cline_refine3D_3%set('sh_first',                 'yes')
        call cline_refine3D_3%set('prob_sh',                  'yes')
        call cline_refine3D_3%set('ml_reg',                   'yes') ! turning on ML regularization to prevent overfitting
        ! rec3D & postproc
        call cline_reconstruct3D%set('prg',         'reconstruct3D')
        call cline_reconstruct3D%set('box',        real(params%box))
        call cline_reconstruct3D%set('projfile',    params%projfile)
        call cline_reconstruct3D%set('ml_reg',                 'no')
        call cline_reconstruct3D%set('needs_sigma',            'no')
        call cline_reconstruct3D%set('objfun',                 'cc')
        call cline_reconstruct3D%set('pgrp',            params%pgrp)
        call cline_reconstruct3D_mlreg%set('prg',   'reconstruct3D')
        call cline_reconstruct3D_mlreg%set('objfun',       'euclid')
        call cline_reconstruct3D_mlreg%set('needs_sigma',     'yes')
        call cline_reconstruct3D_mlreg%set('sigma_est', params%sigma_est)
        call cline_reconstruct3D_mlreg%set('ml_reg',          'yes')
        call cline_postprocess%set('prg',             'postprocess')
        call cline_postprocess%set('projfile',      params%projfile)
        call cline_postprocess%set('imgkind',              vol_type)
        ! SYMMETRY AXIS SEARCH
        if( l_srch4symaxis )then
            call cline_symsrch%set('prg',      'symaxis_search') ! needed for cluster exec
            call cline_symsrch%set('pgrp',          params%pgrp)
            call cline_symsrch%set('smpd',     params%smpd_crop)
            call cline_symsrch%set('box', real(params%box_crop))
            call cline_symsrch%set('projfile',  params%projfile)
            if( .not. cline_symsrch%defined('cenlp') ) call cline_symsrch%set('cenlp', CENLP_DEFAULT)
            call cline_symsrch%set('hp',              params%hp)
            call cline_symsrch%set('oritype',          'ptcl3D')
            call cline_symsrch%delete('lp_auto')
        endif
        ! execute commanders
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> BAYESIAN 3D AB INITIO'
        write(logfhandle,'(A)') '>>>'
        call xrefine3D_distr%execute(cline_refine3D_1)
        iter       = nint(cline_refine3D_1%get_rarg('endit'))
        call spproj%read_segment('ptcl3D', params%projfile)
        lp_est     = spproj%os_ptcl3D%get_avg('lp')
        find_start = calc_fourier_index(lp_est, params%box_crop, params%smpd_crop) - 2
        lplims(1)  = calc_lowpass_lim(find_start, params%box_crop, params%smpd_crop)
        lplims(2)  = calc_lplim_stage2(3) ! low-pass limit is median of three best (as in 2D)
        iter = iter + 1
        call cline_refine3D_2%set('startit', real(iter))
        call cline_refine3D_2%set('lpstart', lplims(1))
        call cline_refine3D_2%set('lpstop',  lplims(2))
        write(logfhandle,'(A,F5.1)') '>>> ESTIMATED         LOW-PASS LIMIT (IN A) TO: ', lp_est
        write(logfhandle,'(A,F5.1)') '>>> DID SET STARTING  LOW-PASS LIMIT (IN A) TO: ', lplims(1)
        write(logfhandle,'(A,F5.1)') '>>> DID SET HARD      LOW-PASS LIMIT (IN A) TO: ', lplims(2)
        call xrefine3D_distr%execute(cline_refine3D_2)
        iter     = nint(cline_refine3D_2%get_rarg('endit'))
        call spproj%read_segment('ptcl3D', params%projfile)
        lp_est   = spproj%os_ptcl3D%get_avg('lp')
        vol_iter = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
        if( l_srch4symaxis )then
            lp_sym = max(LP_SYMSRCH_LB,lp_est)
            write(logfhandle,'(A,F5.1)') '>>> DID SET SYMSEARCH LOW-PASS LIMIT (IN A) TO: ', lp_sym
            call cline_symsrch%set('lp', lp_sym)
            write(logfhandle,'(A)') '>>>'
            write(logfhandle,'(A)') '>>> SYMMETRY AXIS SEARCH'
            write(logfhandle,'(A)') '>>>'
            call cline_symsrch%set('vol1', trim(vol_iter))
            call xsymsrch%execute_shmem(cline_symsrch)
            call del_file('SYMAXIS_SEARCH_FINISHED')
        endif
        ! final phase scaling
        smpd_target = max(params%smpd, lplims(2)*LP2SMPDFAC)
        call downscale(smpd_target, scale_factor2)
        iter = iter + 1
        call cline_refine3D_3%set('startit',             real(iter))
        call cline_refine3D_3%set('box_crop', real(params%box_crop))
        call cline_refine3D_3%set('smpd_crop',     params%smpd_crop)
        call cline_refine3D_3%set('trs',                     trslim) ! trslim set by downscale above
        call xrefine3D_distr%execute(cline_refine3D_3)
        ! for visualization
        do state = 1, params%nstates
            str_state = int2str_pad(state,2)
            call final_vol%new([params%box_crop,params%box_crop,params%box_crop],params%smpd_crop)
            call final_vol%read(trim(VOL_FBODY)//trim(str_state)//trim(params%ext))
            call final_vol%generate_orthogonal_reprojs(reprojs)
            call reprojs%write_jpg('orthogonal_reprojs_state'//trim(str_state)//'.jpg')
            call final_vol%kill
            call reprojs%kill
        enddo
        ! Final reconstruction at original scale
        if( l_autoscale )then
            write(logfhandle,'(A)') '>>>'
            write(logfhandle,'(A)') '>>> RECONSTRUCTION AT ORIGINAL SAMPLING'
            write(logfhandle,'(A)') '>>>'
            ! no ML-filtering
            call cline_reconstruct3D%set('ml_reg',      'no')
            call cline_reconstruct3D%set('needs_sigma', 'no')
            call cline_reconstruct3D%set('objfun',      'cc')
            ! no fractional or stochastic updates
            call cline_reconstruct3D%delete('update_frac')
            call cline_reconstruct3D%delete('stoch_update')
            ! reconstruction
            call xreconstruct3D_distr%execute_shmem(cline_reconstruct3D)
            call spproj%read_segment('out',params%projfile)
            do state = 1, params%nstates
                str_state = int2str_pad(state,2)
                vol       = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
                call spproj%add_vol2os_out(vol, params%smpd, state, vol_type)
                if( trim(params%oritype).eq.'ptcl3D' )then
                    call spproj%add_fsc2os_out(FSC_FBODY//str_state//trim(BIN_EXT), state, params%box)
                endif
            enddo
            call spproj%write_segment_inside('out',params%projfile)
            ! post-processing
            do state = 1, params%nstates
                call cline_postprocess%delete('lp') ! to obtain optimal filtration
                call cline_postprocess%set('state', real(state))
                call xpostprocess%execute(cline_postprocess)
            enddo
        else
            do state = 1, params%nstates
                str_state = int2str_pad(state,2)
                vol       = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
                call final_vol%new([params%box,params%box,params%box],params%smpd)
                call final_vol%read(vol)
                call final_vol%mirror('x')
                call final_vol%write(add2fbody(vol,params%ext,trim(PPROC_SUFFIX)//trim(MIRR_SUFFIX)))
                call final_vol%kill
            enddo
        endif
        do state = 1, params%nstates
            str_state      = int2str_pad(state,2)
            vol            = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
            vol_pproc      = add2fbody(vol,params%ext,PPROC_SUFFIX)
            vol_pproc_mirr = add2fbody(vol,params%ext,trim(PPROC_SUFFIX)//trim(MIRR_SUFFIX))
            if( file_exists(vol)            ) call simple_rename(vol,            trim(REC_FBODY)           //trim(str_state)//trim(params%ext))
            if( file_exists(vol_pproc)      ) call simple_rename(vol_pproc,      trim(REC_PPROC_FBODY)     //trim(str_state)//trim(params%ext))
            if( file_exists(vol_pproc_mirr) ) call simple_rename(vol_pproc_mirr, trim(REC_PPROC_MIRR_FBODY)//trim(str_state)//trim(params%ext))
        enddo
        ! cleanup
        call spproj%kill
        call qsys_cleanup
        call simple_end('**** SIMPLE_ABINITIO_3DMODEL NORMAL STOP ****')

        contains

            subroutine downscale( smpd_target, scale_factor )
                real, intent(in)    :: smpd_target
                real, intent(inout) :: scale_factor
                l_autoscale      = .false.
                scale_factor     = 1.0
                params%box       = spproj%get_box()
                params%smpd_crop = params%smpd
                params%box_crop  = params%box
                params%msk_crop  = params%msk
                if( params%l_autoscale )then
                    call autoscale(params%box, params%smpd, smpd_target, params%box_crop, params%smpd_crop, scale_factor, minbox=MINBOX)
                    l_autoscale = params%box_crop < params%box
                    if( l_autoscale )then
                        params%msk_crop = params%msk * scale_factor1
                        write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',params%box_crop
                        trslim = max(2.0, AHELIX_WIDTH / params%smpd_crop / 2.0)
                    else
                        trslim = max(2.0, AHELIX_WIDTH / params%smpd / 2.0)
                    endif
                endif
            end subroutine downscale

            function calc_lplim_stage2( nbest ) result( lplim )
                integer, intent(in)  :: nbest
                real,    allocatable :: res(:), tmp_rarr(:)
                integer, allocatable :: states(:), tmp_iarr(:)
                real :: lplim
                tmp_rarr  = spproj%os_cls2D%get_all('res')
                tmp_iarr  = nint(spproj%os_cls2D%get_all('state'))
                res       = pack(tmp_rarr, mask=(tmp_iarr>0))
                call hpsort(res)
                lplim = median_nocopy(res(:nbest))
                deallocate(tmp_rarr, tmp_iarr, res)
            end function calc_lplim_stage2

    end subroutine exec_abinitio_3Dmodel_autolp

    !> for generation of an initial 3d model from particles
    subroutine exec_abinitio_3Dmodel( self, cline )
        use simple_convergence, only: convergence
        class(abinitio_3Dmodel_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        real,    parameter :: SCALEFAC        = 0.667
        real,    parameter :: CENLP_DEFAULT   = 30.
        real,    parameter :: LP_DEFAULT      = 6.
        real,    parameter :: LPSTART_DEFAULT = 30., LPSTOP_DEFAULT=LP_DEFAULT
        integer, parameter :: MINBOX  = 64
        integer, parameter :: NSTAGES_DEFAULT = 5
        integer, parameter :: MAXITS2 = 30, MAXITS_SHORT1 = 15, MAXITS_SHORT2 = 25
        integer, parameter :: NSPACE1 = 500, NSPACE2 = 1000, NSPACE3 = 2000
        integer, parameter :: SYMSEARCH_DEFAULT   = 3               ! in [1;NSTAGES]
        integer, parameter :: MLREG_ITER          = 1               ! in [1;NSTAGES+1]
        integer, parameter :: SHIFT_STAGE_DEFAULT = 4
        integer, parameter :: ICM_STAGE_DEFAULT   = NSTAGES_DEFAULT ! in [1;NSTAGES+1]
        ! commanders
        type(refine3D_commander_distr)      :: xrefine3D_distr
        type(reconstruct3D_commander_distr) :: xreconstruct3D_distr
        type(postprocess_commander)         :: xpostprocess
        type(symaxis_search_commander)      :: xsymsrch
        ! command lines
        type(cmdline)                 :: cline_refine3D, cline_reconstruct3D, cline_reconstruct3D_mlreg
        type(cmdline)                 :: cline_postprocess, cline_symsrch
        ! other
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(convergence)             :: conv
        type(sym)                     :: se1, se2
        type(class_frcs)              :: clsfrcs
        type(image)                   :: final_vol, reprojs
        character(len=:), allocatable :: vol_type, str_state, vol, vol_pproc, vol_pproc_mirr, frcs_fname
        character(len=LONGSTRLEN)     :: vol_str
        real    :: smpd_target, lp_target, scale, trslim, cenlp, symlp, dummy
        integer :: iter, it, prev_box_crop, maxits, state, frac_maxits_incr, maxits_glob
        integer :: nstages, symsearch_iter
        logical :: l_autoscale, l_lpset, l_err, l_srch4symaxis, l_symran, l_sym, l_lpstop_set
        logical :: l_lpstart_set
        if( .not. cline%defined('mkdir')        ) call cline%set('mkdir',        'yes')
        if( .not. cline%defined('refine')       ) call cline%set('refine',      'prob')
        if( .not. cline%defined('autoscale')    ) call cline%set('autoscale',    'yes')
        if( .not. cline%defined('ml_reg')       ) call cline%set('ml_reg',       'yes')
        if( .not. cline%defined('sigma_est')    ) call cline%set('sigma_est', 'global')
        if( .not. cline%defined('prob_sh')      ) call cline%set('prob_sh',      'yes')
        if( .not. cline%defined('sh_first')     ) call cline%set('sh_first',     'yes')
        if( .not. cline%defined('prob_athres')  ) call cline%set('prob_athres',    10.)
        ! if( .not. cline%defined('stoch_update') ) call cline%set('stoch_update', 'yes') ! off 4 now
        call cline%set('stoch_update', 'no')
        if( .not. cline%defined('center')       ) call cline%set('center',        'no')
        if( .not. cline%defined('objfun')       ) call cline%set('objfun',    'euclid')
        if( .not. cline%defined('oritype')      ) call cline%set('oritype',   'ptcl3D')
        if( .not. cline%defined('pgrp')         ) call cline%set('pgrp',          'c1')
        if( .not. cline%defined('pgrp_start')   ) call cline%set('pgrp_start',    'c1')
        if( .not. cline%defined('shift_stage')  ) call cline%set('shift_stage', SHIFT_STAGE_DEFAULT)
        if( .not. cline%defined('icm')          ) call cline%set('icm',           'no')
        if( .not. cline%defined('icm_stage')    ) call cline%set('icm_stage', ICM_STAGE_DEFAULT)
        if( .not. cline%defined('ptclw')        ) call cline%set('ptclw',         'no')
        ! resolution limit strategy
        l_lpset       = .false.
        l_lpstop_set  = cline%defined('lpstop')
        l_lpstart_set = cline%defined('lpstart')
        if( cline%defined('lp') )then
            if( l_lpstart_set .or. l_lpstop_set )then
                THROW_HARD('One of LP or LPSTART & LPSTOP must be defined!')
            endif
            l_lpset = .true.
        else
            if( .not.l_lpstart_set ) call cline%set('lpstart',LPSTART_DEFAULT)
            if( .not.l_lpstop_set  ) call cline%set('lpstop', LPSTOP_DEFAULT)
        endif
        ! make master parameters
        if( cline%defined('update_frac') ) call cline%delete('stoch_update')
        call params%new(cline)
        call cline%set('mkdir', 'no')
        call cline%delete('autoscale')
        call cline%delete('lp')
        call cline%delete('shift_stage')
        call cline%delete('icm_stage')
        call cline%delete('icm')
        ! stages specific parameters
        nstages          = NSTAGES_DEFAULT
        symsearch_iter   = SYMSEARCH_DEFAULT
        frac_maxits_incr = 0
        maxits_glob      = MAXITS_SHORT1 * (nstages - 2) + MAXITS_SHORT2 * 2 + MAXITS2
        if( params%l_frac_update .and. (.not.params%l_stoch_update) )then
            ! adjusting number ot iterations for frac_update alone
            frac_maxits_incr = max(2,min(10,nint(1./params%update_frac)))
            maxits_glob      = maxits_glob + (nstages+1)*frac_maxits_incr
        endif
        if( l_lpset )then
            params%lpstop  = params%lp
            params%lpstart = params%lp
        endif
        if( params%shift_stage < 1 .or. params%shift_stage > nstages+1 )then
            params%shift_stage = min(nstages+1,max(1,params%shift_stage))
            THROW_WARN('SHIFT_STAGE out of range, defaulting to: '//int2str(params%shift_stage))
        endif
        if( params%l_icm .and. (params%icm_stage < 1 .or. params%icm_stage > nstages+1) )then
            params%icm_stage = min(nstages+1,max(1,params%icm_stage))
            THROW_WARN('ICM_STAGE out of range, defaulting to: '//int2str(params%icm_stage))
        endif
        ! read project
        call spproj%read(params%projfile)
        call spproj%update_projinfo(cline)
        call spproj%write_segment_inside('projinfo', params%projfile)
        if( .not. cline%defined('vol1') )then
            ! randomize projection directions
            select case(trim(params%oritype))
                case('cls3D')
                    if( spproj%os_cls3D%get_noris() < 1 )then
                        THROW_HARD('Class averages could not be found in the project')
                    endif
                    vol_type = 'vol_cavg'
                    call spproj%os_cls3D%rnd_oris
                    call spproj%os_cls3D%set_all2single('stkind', 1.)
                    call spproj%os_cls3D%set_all2single('w',1.)
                case('ptcl3D')
                    if( spproj%os_ptcl3D%get_noris() < 1 )then
                        THROW_HARD('Particles could not be found in the project')
                    endif
                    vol_type = 'vol'
                    call spproj%os_ptcl3D%rnd_oris
                    call spproj%os_ptcl3D%set_all2single('w',1.)
                case DEFAULT
                    THROW_HARD('Unsupported ORITYPE; exec_abinitio_3Dmodel')
            end select
            call spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        write(logfhandle,'(A,F5.1)') '>>> STARTING RESOLUTION LIMIT (IN A): ', params%lpstart
        write(logfhandle,'(A,F5.1)') '>>> HARD     RESOLUTION LIMIT (IN A): ', params%lpstop
        if( trim(params%center).eq.'yes' )then
            write(logfhandle,'(A,F5.1)') '>>> CENTERING  LOW-PASS LIMIT (IN A): ', params%cenlp
        endif
        ! centering & symmetry resolution limit
        call mskdiam2lplimits(params%mskdiam, symlp, dummy, cenlp)
        if( .not. cline%defined('cenlp') )then
            params%cenlp = cenlp
            call cline%set('cenlp', params%cenlp)
        endif
        ! symmetry
        if( l_lpset )then
            ! from mskdiam2lplimits lpstart above
        else
            symlp = (params%lpstart+params%lpstop)/2.
        endif
        l_srch4symaxis = trim(params%pgrp) .ne. trim(params%pgrp_start)
        l_symran       = .false.
        l_sym          = l_srch4symaxis
        if( params%pgrp_start.ne.'c1' .or. params%pgrp.ne.'c1' )then
            se1 = sym(params%pgrp_start)
            se2 = sym(params%pgrp)
            if(se1%get_nsym() > se2%get_nsym())then
                ! ensure se2 is a subgroup of se1
                if( .not. se1%has_subgrp(params%pgrp) )THROW_HARD('Incompatible symmetry groups; exec_abinitio_3Dmodel')
                ! set flag for symmetry randomisation
                ! in case we are moving from a higher to lower group
                l_symran = .true.
            else if( se2%get_nsym() > se1%get_nsym() )then
                ! ensure se1 is a subgroup of se2
                if( .not. se2%has_subgrp(params%pgrp_start) )THROW_HARD('Incompatible symmetry groups; exec_abinitio_3Dmodel')
            endif
        endif
        ! dimensions defaults
        params%box       = spproj%get_box()
        params%smpd_crop = params%smpd
        params%box_crop  = params%box
        l_autoscale      = .false.
        call spproj%kill ! not needed anymore
        ! command-lines
        cline_refine3D            = cline
        cline_reconstruct3D       = cline
        cline_postprocess         = cline
        cline_symsrch             = cline
        cline_reconstruct3D_mlreg = cline_reconstruct3D
        call cline_refine3D%set('prg',                'refine3D')
        call cline_refine3D%set('projfile',      params%projfile)
        call cline_refine3D%set('pgrp',        params%pgrp_start)
        if( params%l_lpauto )then
            call cline_refine3D%set('lpstart', params%lpstart)
            call cline_refine3D%set('lpstop',  params%lpstop)
        endif
        call cline_reconstruct3D%set('prg',      'reconstruct3D')
        call cline_reconstruct3D%set('box',     real(params%box))
        call cline_reconstruct3D%set('projfile', params%projfile)
        call cline_reconstruct3D%set('ml_reg',              'no')
        call cline_reconstruct3D%set('needs_sigma',         'no')
        call cline_reconstruct3D%set('objfun',              'cc')
        call cline_reconstruct3D%set('pgrp',         params%pgrp)
        call cline_postprocess%set('prg',          'postprocess')
        call cline_postprocess%set('projfile',   params%projfile)
        call cline_postprocess%set('imgkind',           vol_type)
        if( l_srch4symaxis )then
            call cline_symsrch%set('prg',     'symaxis_search') ! needed for cluster exec
            call cline_symsrch%set('pgrp',     params%pgrp)
            call cline_symsrch%set('projfile', params%projfile)
            call cline_symsrch%set('hp',       params%hp)
            call cline_symsrch%set('center',   'yes')
        endif
        call cline_reconstruct3D_mlreg%set('prg',         'reconstruct3D')
        call cline_reconstruct3D_mlreg%set('objfun',      'euclid')
        call cline_reconstruct3D_mlreg%set('needs_sigma', 'yes')
        call cline_reconstruct3D_mlreg%set('sigma_est',   params%sigma_est)
        call cline_reconstruct3D_mlreg%set('ml_reg',      'yes')
        ! Frequency marching
        iter = 0
        do it = 1,nstages
            write(logfhandle,'(A)')'>>>'
            prev_box_crop = params%box_crop
            ! resolution limit
            if( it == 1 )then
                params%lp = params%lpstart
            else
                params%lp = max(params%lpstop, params%lpstop+(params%lp-params%lpstop)/2.)
            endif
            write(logfhandle,'(A,I3,A9,F5.1)')'>>> STAGE ',it,' WITH LP =',params%lp
            ! dimensions
            prev_box_crop = params%box_crop
            if( params%l_autoscale )then
                l_autoscale = .false.
                lp_target   = params%lp * SCALEFAC
                smpd_target = max(params%smpd, lp_target/2.)
                call autoscale(params%box, params%smpd, smpd_target, params%box_crop, params%smpd_crop, scale, minbox=MINBOX)
                trslim      = max(2.0, AHELIX_WIDTH / params%smpd_crop / 2.0)
                l_autoscale = params%box_crop < params%box
            else
                trslim      = max(2.0, AHELIX_WIDTH / params%smpd / 2.0)
            endif
            if( cline%defined('trs') ) trslim = params%trs
            if( l_autoscale )then
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',params%box_crop
            endif
            call cline_refine3D%set('center', params%center)
            ! stage updates
            call cline_refine3D%set('box_crop', params%box_crop)
            call cline_refine3D%set('lp',       params%lp)
            call cline_refine3D%set('startit',  iter+1)
            ! # of iterations
            maxits = MAXITS_SHORT1
            if( it >= nstages-1 ) maxits = MAXITS_SHORT2
            ! projection directions & shift
            if( it < params%shift_stage )then
                call cline_refine3D%set('nspace', NSPACE1)
                call cline_refine3D%set('trs',    0.)
            else
                call cline_refine3D%set('nspace', NSPACE2)
                call cline_refine3D%set('trs',    trslim)
            end if
            maxits = maxits + frac_maxits_incr
            call cline_refine3D%set('maxits',      maxits)
            call cline_refine3D%set('maxits_glob', maxits_glob)
            ! ICM filtering
            if( params%l_icm .and. (it >= params%icm_stage) )then
                call cline_refine3D%set('icm',   'yes')
                call cline_refine3D%set('lambda', params%lambda)
            endif
            ! Volume: ML regularization, symmetry
            call cline_refine3D%set('ml_reg', 'no')
            if( it > 1 ) call cline_refine3D%set('continue',  'yes')
            ! reconstruction is required after symmetry search
            if( it == symsearch_iter+1 .and. (l_srch4symaxis .or. l_symran) )then
                call cline_refine3D%delete('continue')
                do state = 1,params%nstates
                    vol_str = 'vol'//trim(int2str(state))
                    call cline_refine3D%delete(vol_str)
                enddo
                call cline_refine3D%set('pgrp', params%pgrp)
                l_srch4symaxis = .false.
                l_symran       = .false.
            endif
            ! reconstruction is required when switching to MLreg
            if( params%l_ml_reg .and. it >= MLREG_ITER )then
                call cline_refine3D%set('ml_reg', 'yes')
                if( it == MLREG_ITER .and. it > 1 )then
                    call cline_reconstruct3D_mlreg%set('which_iter', iter)
                    call cline_reconstruct3D_mlreg%set('box_crop',   params%box_crop)
                    call cline_reconstruct3D_mlreg%set('smpd_crop',  params%smpd_crop)
                    ! nspace must match refine3D when projrec=yes
                    call cline_reconstruct3D_mlreg%set('nspace',     cline_refine3D%get_rarg('nspace'))
                    ! reconstruction
                    call xreconstruct3D_distr%execute_shmem(cline_reconstruct3D_mlreg)
                    call cline_refine3D%set('continue',  'yes')
                    do state = 1,params%nstates
                        vol_str = 'vol'//trim(int2str(state))
                        call cline_refine3D%delete(vol_str)
                    enddo
                endif
            endif
            ! Execution
            call exec_refine3D(iter)
            ! Symmetrization
            if( it == SYMSEARCH_ITER ) call symmetrize
        enddo
        ! Final stage
        it = nstages+1
        write(logfhandle,'(A)')'>>>'
        write(logfhandle,'(A,F6.2)')'>>> FINAL STAGE WITH LP=',params%lpstop
        prev_box_crop = params%box_crop
        if( params%l_autoscale )then
            l_autoscale = .false.
            lp_target   = params%lpstop * SCALEFAC
            smpd_target = max(params%smpd, lp_target/2.)
            call autoscale(params%box, params%smpd, smpd_target, params%box_crop, params%smpd_crop, scale, minbox=MINBOX)
            trslim      = max(MINSHIFT, AHELIX_WIDTH / params%smpd_crop)
            l_autoscale = params%box_crop < params%box
        else
            trslim      = max(MINSHIFT, AHELIX_WIDTH / params%smpd)
        endif
        if( l_autoscale )then
            write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',params%box_crop
        endif
        call cline_refine3D%set('continue', 'yes')
        call cline_refine3D%set('ml_reg',   'no')
        if( params%l_ml_reg )then
            call cline_refine3D%set('ml_reg', 'yes')
            if( it == MLREG_ITER )then
                ! reconstruction required
                call cline_reconstruct3D_mlreg%set('which_iter', iter)
                call cline_reconstruct3D_mlreg%set('box_crop',   params%box_crop)
                call cline_reconstruct3D_mlreg%set('smpd_crop',  params%smpd_crop)
                call xreconstruct3D_distr%execute_shmem(cline_reconstruct3D_mlreg)
                do state = 1,params%nstates
                    vol_str = 'vol'//trim(int2str(state))
                    call cline_refine3D%delete(vol_str)
                enddo
            endif
        endif
        if( params%l_icm .and. (it >= params%icm_stage) )then
            call cline_refine3D%set('icm',   'yes')
            call cline_refine3D%set('lambda', params%lambda)
        endif
        ! Final step
        call cline_refine3D%set('box_crop',   params%box_crop)
        call cline_refine3D%set('lp',         params%lpstop)
        call cline_refine3D%set('trs',        trslim)
        call cline_refine3D%set('startit',    iter+1)
        call cline_refine3D%set('maxits',     MAXITS2+frac_maxits_incr)
        call cline_refine3D%set('nspace',     NSPACE3)
        ! execution
        call exec_refine3D(iter)
        ! for visualization
        do state = 1, params%nstates
            str_state = int2str_pad(state,2)
            call final_vol%new([params%box_crop,params%box_crop,params%box_crop],params%smpd_crop)
            call final_vol%read(trim(VOL_FBODY)//trim(str_state)//trim(params%ext))
            call final_vol%generate_orthogonal_reprojs(reprojs)
            call reprojs%write_jpg('orthogonal_reprojs_state'//trim(str_state)//'.jpg')
            call final_vol%kill
            call reprojs%kill
        enddo
        ! Final reconstruction at original scale
        if( l_autoscale )then
            write(logfhandle,'(A)') '>>>'
            write(logfhandle,'(A)') '>>> RECONSTRUCTION AT ORIGINAL SAMPLING'
            write(logfhandle,'(A)') '>>>'
            ! no ML-filtering
            call cline_reconstruct3D%set('ml_reg',      'no')
            call cline_reconstruct3D%set('needs_sigma', 'no')
            call cline_reconstruct3D%set('objfun',      'cc')
            ! no fractional or stochastic updates
            call cline_reconstruct3D%delete('update_frac')
            call cline_reconstruct3D%delete('stoch_update')
            ! individual particles reconstruction
            call cline_reconstruct3D%set('projrec', 'no')
            ! reconstruction
            call xreconstruct3D_distr%execute_shmem(cline_reconstruct3D)
            call spproj%read_segment('out',params%projfile)
            do state = 1, params%nstates
                str_state = int2str_pad(state,2)
                vol       = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
                call spproj%add_vol2os_out(vol, params%smpd, state, vol_type)
                if( trim(params%oritype).eq.'ptcl3D' )then
                    call spproj%add_fsc2os_out(FSC_FBODY//str_state//trim(BIN_EXT), state, params%box)
                endif
            enddo
            call spproj%write_segment_inside('out',params%projfile)
            ! post-processing
            do state = 1, params%nstates
                call cline_postprocess%delete('lp') ! so as to obtain optimal filtration
                call cline_postprocess%set('state', real(state))
                call xpostprocess%execute(cline_postprocess)
            enddo
        else
            do state = 1, params%nstates
                str_state = int2str_pad(state,2)
                vol       = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
                call final_vol%new([params%box,params%box,params%box],params%smpd)
                call final_vol%read(vol)
                call final_vol%mirror('x')
                call final_vol%write(add2fbody(vol,params%ext,trim(PPROC_SUFFIX)//trim(MIRR_SUFFIX)))
                call final_vol%kill
            enddo
        endif
        do state = 1, params%nstates
            str_state      = int2str_pad(state,2)
            vol            = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
            vol_pproc      = add2fbody(vol,params%ext,PPROC_SUFFIX)
            vol_pproc_mirr = add2fbody(vol,params%ext,trim(PPROC_SUFFIX)//trim(MIRR_SUFFIX))
            if( file_exists(vol)            ) call simple_rename(vol,            trim(REC_FBODY)           //trim(str_state)//trim(params%ext))
            if( file_exists(vol_pproc)      ) call simple_rename(vol_pproc,      trim(REC_PPROC_FBODY)     //trim(str_state)//trim(params%ext))
            if( file_exists(vol_pproc_mirr) ) call simple_rename(vol_pproc_mirr, trim(REC_PPROC_MIRR_FBODY)//trim(str_state)//trim(params%ext))
        enddo
        ! transfer cls3D parameters to particles
        if( trim(params%oritype) .eq. 'cls3D' )then
            call spproj%read_segment('cls3D', params%projfile)
            call spproj%read_segment('ptcl3D',params%projfile)
            call spproj%map2ptcls
            call spproj%write_segment_inside('ptcl3D',params%projfile)
            call spproj%os_cls3D%delete_entry('stkind')
            call spproj%write_segment_inside('cls3D',params%projfile)
        endif
        ! cleanup
        call spproj%kill
        call qsys_cleanup
        call simple_end('**** SIMPLE_ABINITIO_3DMODEL NORMAL STOP ****')
        contains

            subroutine exec_refine3D( iter )
                integer,          intent(out) :: iter
                character(len=:), allocatable :: stage
                call cline_refine3D%delete('endit')
                call xrefine3D_distr%execute_shmem(cline_refine3D)
                call conv%read(l_err)
                iter = nint(conv%get('iter'))
                call del_files(DIST_FBODY,      params_glob%nparts,ext='.dat')
                call del_files(ASSIGNMENT_FBODY,params_glob%nparts,ext='.dat')
                call del_file(trim(DIST_FBODY)      //'.dat')
                call del_file(trim(ASSIGNMENT_FBODY)//'.dat')
                if( it <= nstages )then
                    stage = '_stage_'//int2str(it)
                    do state = 1, params%nstates
                        str_state = int2str_pad(state,2)
                        vol       = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
                        vol_pproc = add2fbody(vol,params%ext,PPROC_SUFFIX)
                        if( file_exists(vol)      ) call simple_copy_file(vol,       add2fbody(vol,      params%ext,stage))
                        if( file_exists(vol_pproc)) call simple_copy_file(vol_pproc, add2fbody(vol_pproc,params%ext,stage))
                    enddo
                endif
            end subroutine exec_refine3D

            subroutine symmetrize()
                if( l_symran )then
                    call spproj%read_segment(params%oritype, params%projfile)
                    select case(trim(params%oritype))
                    case('cls3D')
                        call se1%symrandomize(spproj%os_cls3D)
                    case('ptcl3D')
                        call se1%symrandomize(spproj%os_ptcl3D)
                    end select
                    call spproj%write_segment_inside(params%oritype, params%projfile)
                endif
                if( l_srch4symaxis )then
                    write(logfhandle,'(A)') '>>>'
                    write(logfhandle,'(A)') '>>> SYMMETRY AXIS SEARCH'
                    write(logfhandle,'(A)') '>>>'
                    symlp = max(symlp, params%lp)
                    call cline_symsrch%set('lp',       symlp)
                    call cline_symsrch%set('box_crop', params%box_crop)
                    do state = 1,params%nstates
                        vol_str   = 'vol'//trim(int2str(state))
                        str_state = int2str_pad(state,2)
                        call cline_symsrch%set(vol_str, trim(VOL_FBODY)//trim(str_state)//trim(params%ext))
                    enddo
                    call xsymsrch%execute_shmem(cline_symsrch)
                    call del_file('SYMAXIS_SEARCH_FINISHED')
                endif
            end subroutine symmetrize

    end subroutine exec_abinitio_3Dmodel

    !> for generation of an initial 3d model from particles
    subroutine exec_batch_abinitio_3Dmodel( self, cline )
        use simple_convergence, only: convergence
        class(batch_abinitio_3Dmodel_commander), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        real,    parameter :: SCALEFAC      = 0.667
        real,    parameter :: CENLP_DEFAULT = 30.
        real,    parameter :: LP_DEFAULT    = 6.
        real,    parameter :: LPSTART_DEFAULT = 30., LPSTOP_DEFAULT=LP_DEFAULT
        real,    parameter :: BATCHFRAC_DEFAULT = 0.2
        integer, parameter :: MINBOX  = 64
        integer, parameter :: NSTAGES = 11
        integer, parameter :: NSPACE1=500, NSPACE2=1000, NSPACE3=2000
        integer, parameter :: SHIFT_STAGE_DEFAULT = NSTAGES-2
        integer, parameter :: ICM_STAGE_DEFAULT   = 2       ! in [1;NSTAGES]
        integer, parameter :: MIN_NPTCLS = 20000
        integer, parameter :: MAX_NPTCLS = 500000
        integer, parameter :: MAXITS_BETWEEN_DEFAULT = 4
        integer, parameter :: MLREG_ITER     = 1
        integer, parameter :: SYMSEARCH_ITER = 5
        ! commanders
        type(refine3D_commander_distr)      :: xrefine3D_distr
        type(reconstruct3D_commander_distr) :: xreconstruct3D_distr
        type(postprocess_commander)         :: xpostprocess
        type(calc_pspec_commander_distr)    :: xcalc_pspec_distr
        type(symaxis_search_commander)      :: xsymsrch
        ! command lines
        type(cmdline)                 :: cline_calc_pspec_distr
        type(cmdline)                 :: cline_refine3D, cline_reconstruct3D
        type(cmdline)                 :: cline_postprocess, cline_symsrch
        ! other
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(convergence)             :: conv
        type(image)                   :: final_vol, reprojs
        type(ran_tabu)                :: rt
        type(sym)                     :: se1, se2
        character(len=:), allocatable :: vol_type, str_state, vol, vol_pproc, vol_pproc_mirr
        integer,          allocatable :: vec(:), states(:), batches(:), pinds(:), sampled(:), cnts(:)
        logical,          allocatable :: mask(:)
        character(len=LONGSTRLEN)     :: vol_str
        real    :: lps(NSTAGES), smpds(NSTAGES), trs(NSTAGES)
        integer :: boxs(NSTAGES), nspaces(NSTAGES)
        real    :: smpd_target, lp_target, scale, trslim, cenlp, symlp, dummy
        integer :: iter, it, prev_box_crop, maxits, state, i,j, final_nptcls, nptcls_sel
        integer :: ini_nptcls, nptcls, iters_per_stage
        logical :: l_lpset, l_err, l_lpstop_set, l_lpstart_set, l_srch4symaxis, l_symran, l_sym
        call cline%set('stoch_update', 'no')
        call cline%set('oritype',      'ptcl3D')
        call cline%set('sigma_est',    'global')
        call cline%set('center',       'no')
        if( .not. cline%defined('mkdir')        ) call cline%set('mkdir',          'yes')
        if( .not. cline%defined('ml_reg')       )  call cline%set('ml_reg',        'yes')
        if( .not. cline%defined('refine')       )  call cline%set('refine',        'prob')
        if( .not. cline%defined('prob_sh')      )  call cline%set('prob_sh',       'yes')
        if( .not. cline%defined('prob_athres')  )  call cline%set('prob_athres',   10.)
        if( .not. cline%defined('objfun')       )  call cline%set('objfun',        'euclid')
        if( .not. cline%defined('shift_stage')  )  call cline%set('shift_stage',   SHIFT_STAGE_DEFAULT)
        if( .not. cline%defined('batchfrac')    )  call cline%set('batchfrac',     BATCHFRAC_DEFAULT)
        if( .not. cline%defined('icm')          )  call cline%set('icm',           'no')
        if( .not. cline%defined('icm_stage')    )  call cline%set('icm_stage',     ICM_STAGE_DEFAULT)
        if( .not. cline%defined('maxits_between')) call cline%set('maxits_between',MAXITS_BETWEEN_DEFAULT)
        if( .not. cline%defined('pgrp')         )  call cline%set('pgrp',          'c1')
        if( .not. cline%defined('pgrp_start')   )  call cline%set('pgrp_start',    'c1')
        ! resolution limit strategy
        l_lpset       = .false.
        l_lpstop_set  = cline%defined('lpstop')
        l_lpstart_set = cline%defined('lpstart')
        if( cline%defined('lp') )then
            if( l_lpstart_set .or. l_lpstop_set )then
                THROW_HARD('One of LP or LPSTART & LPSTOP must be defined!')
            endif
            l_lpset = .true.
        else
            if( .not.l_lpstart_set ) call cline%set('lpstart',LPSTART_DEFAULT)
            if( .not.l_lpstop_set  ) call cline%set('lpstop', LPSTOP_DEFAULT)
        endif
        ! make master parameters
        call params%new(cline)
        call cline%set('mkdir', 'no')
        call cline%delete('lp')
        call cline%delete('shift_stage')
        call cline%delete('icm_stage')
        call cline%delete('icm')
        call cline%delete('lambda')
        if( params%l_icm .and. (params%icm_stage < 1 .or. params%icm_stage > NSTAGES) )then
            params%icm_stage = min(NSTAGES,max(1,params%icm_stage))
            THROW_WARN('ICM_STAGE out of range, defaulting to: '//int2str(params%icm_stage))
        endif
        if( l_lpset )then
            params%lpstop  = params%lp
            params%lpstart = params%lp
        endif
        if( params%shift_stage < 1 .or. params%shift_stage > NSTAGES )then
            params%shift_stage = min(NSTAGES,max(1,params%shift_stage))
            THROW_WARN('SHIFT_STAGE out of range, defaulting to: '//int2str(params%shift_stage))
        endif
        trslim = max(MINSHIFT, AHELIX_WIDTH / params%smpd)
        if( cline%defined('trs') )then
            trslim = params%trs
            call cline%delete('trs')
        endif
        ! read project
        call spproj%read(params%projfile)
        call spproj%update_projinfo(cline)
        call spproj%write_segment_inside('projinfo', params%projfile)
        if( .not. cline%defined('vol1') )then
            ! randomize projection directions
            select case(trim(params%oritype))
                case('ptcl3D')
                    if( spproj%os_ptcl3D%get_noris() < 1 )then
                        THROW_HARD('Particles could not be found in the project')
                    endif
                    vol_type = 'vol'
                    call spproj%os_ptcl3D%rnd_oris
                    call spproj%os_ptcl3D%clean_updatecnt_sampled
                    call spproj%os_ptcl3D%set_all2single('batch',0.)
                    if( spproj%os_ptcl3D%get_nevenodd() == 0 ) call spproj%os_ptcl3D%partition_eo
                    call spproj%write_segment_inside(params%oritype, params%projfile)
                case DEFAULT
                    THROW_HARD('Unsupported ORITYPE; exec_abinitio_3Dmodel')
            end select
        endif
        write(logfhandle,'(A,F5.1)') '>>> STARTING RESOLUTION LIMIT (IN A): ', params%lpstart
        write(logfhandle,'(A,F5.1)') '>>> HARD     RESOLUTION LIMIT (IN A): ', params%lpstop
        if( trim(params%center).eq.'yes' )then
            write(logfhandle,'(A,F5.1)') '>>> CENTERING  LOW-PASS LIMIT (IN A): ', params%cenlp
        endif
        ! centering & symmetry resolution limit
        call mskdiam2lplimits(params%mskdiam, symlp, dummy, cenlp)
        if( .not. cline%defined('cenlp') )then
            params%cenlp = cenlp
            call cline%set('cenlp', params%cenlp)
        endif
        ! symmetry
        l_srch4symaxis = trim(params%pgrp) .ne. trim(params%pgrp_start)
        l_symran       = .false.
        l_sym          = l_srch4symaxis
        if( params%pgrp_start.ne.'c1' .or. params%pgrp.ne.'c1' )then
            se1 = sym(params%pgrp_start)
            se2 = sym(params%pgrp)
            if(se1%get_nsym() > se2%get_nsym())then
                ! ensure se2 is a subgroup of se1
                if( .not. se1%has_subgrp(params%pgrp) )THROW_HARD('Incompatible symmetry groups; exec_abinitio_3Dmodel')
                ! set flag for symmetry randomisation
                ! in case we are moving from a higher to lower group
                l_symran = .true.
            else if( se2%get_nsym() > se1%get_nsym() )then
                ! ensure se1 is a subgroup of se2
                if( .not. se2%has_subgrp(params%pgrp_start) )THROW_HARD('Incompatible symmetry groups; exec_abinitio_3Dmodel')
            endif
        endif
        ! dimensions defaults
        params%box       = spproj%get_box()
        params%smpd_crop = params%smpd
        params%box_crop  = params%box
        ! command-lines
        cline_refine3D            = cline
        cline_reconstruct3D       = cline
        cline_postprocess         = cline
        cline_calc_pspec_distr    = cline
        cline_symsrch             = cline
        call cline_refine3D%set('prg',                'refine3D')
        call cline_refine3D%set('projfile',      params%projfile)
        call cline_refine3D%set('pgrp',        params%pgrp_start)
        call cline_refine3D%set('ml_reg',                   'no')
        call cline_reconstruct3D%set('prg',      'reconstruct3D')
        call cline_reconstruct3D%set('box',     real(params%box))
        call cline_reconstruct3D%set('projfile', params%projfile)
        call cline_reconstruct3D%set('ml_reg',              'no')
        call cline_reconstruct3D%set('needs_sigma',         'no')
        call cline_reconstruct3D%set('objfun',              'cc')
        call cline_reconstruct3D%set('pgrp',   params%pgrp_start)
        call cline_reconstruct3D%delete('batchfrac')
        call cline_postprocess%set('prg',          'postprocess')
        call cline_postprocess%set('projfile',   params%projfile)
        call cline_postprocess%set('imgkind',           vol_type)
        call cline_calc_pspec_distr%set('prg',      'calc_pspec')
        if( l_srch4symaxis )then
            call cline_symsrch%set('prg',     'symaxis_search') ! needed for cluster exec
            call cline_symsrch%set('pgrp',     params%pgrp)
            call cline_symsrch%set('projfile', params%projfile)
            call cline_symsrch%set('hp',       params%hp)
            call cline_symsrch%set('center',   'yes')
        endif
        ! Frequency marching plan
        lps(1)       = params%lpstart
        lps(NSTAGES) = params%lpstop
        do it = 2,NSTAGES-1
            lps(it) = params%lpstop + (params%lpstart-params%lpstop) * real(NSTAGES-it) / real(NSTAGES)
        enddo
        if( l_lpset )then
            ! from mskdiam2lplimits lpstart above
        else
            symlp = (params%lpstart+params%lpstop)/2.
        endif
        ! dimensions
        do it = 1,NSTAGES
            lp_target   = lps(it) * SCALEFAC
            smpd_target = max(params%smpd, lp_target/2.)
            call autoscale(params%box, params%smpd, smpd_target, boxs(it), smpds(it), scale, minbox=MINBOX)
        enddo
        ! shift search & number of projections directions
        do it = 1,NSTAGES
            trs(it)     = 0.
            nspaces(it) = NSPACE1
            if( it >= params%shift_stage )then
                trs(it)     = trslim * params%smpd / smpds(it)
                nspaces(it) = NSPACE2
            endif
        enddo
        nspaces(NSTAGES) = NSPACE3
        ! Batch plan
        iters_per_stage  = ceiling(1./params%batchfrac)
        params%batchfrac = 1./real(iters_per_stage)
        states           = nint(spproj%os_ptcl3D%get_all('state'))
        nptcls_sel       = count(states==1)
        select case(trim(params%algorithm))
        case('fast')
            ! Number of particles increases throughout optimization
            final_nptcls = min(nptcls_sel, MAX_NPTCLS)
            ini_nptcls   = min(nptcls_sel, MIN_NPTCLS)
            allocate(vec(nptcls_sel),source=0)
            do it = NSTAGES,1,-1
                nptcls = ini_nptcls + nint(real((it-1)*(final_nptcls-ini_nptcls)) / real(NSTAGES-1))
                vec(1:nptcls) = it
            enddo
            call seed_rnd
            rt = ran_tabu(nptcls_sel)
            call rt%shuffle(vec)
            call rt%kill
            batches = states
            j = 0
            do i = 1,params%nptcls
                if( batches(i) == 1 )then
                    j = j + 1
                    batches(i) = vec(j)
                endif
            enddo
            deallocate(vec)
            ! Main loop
            iter = 0
            do it = 1,NSTAGES
                write(logfhandle,'(A)') '>>>'
                ! Dimensions & resolution limit
                params%box_crop  = boxs(it)
                params%smpd_crop = smpds(it)
                params%lp        = lps(it)
                ! Particles sampling
                allocate(vec(params%nptcls),source=0)
                where( (batches>0) .and. (batches<=it) ) vec = 1
                nptcls = count(vec==1)
                call spproj%read_segment('ptcl3D',     params%projfile)
                call spproj%os_ptcl3D%set_all('state', real(vec))
                if( it == 1 ) call spproj%os_ptcl3D%set_all('updatecnt', real(vec))
                call spproj%os_ptcl3D%set_all2single('sampled', 0.)
                call spproj%write_segment_inside('ptcl3D', params%projfile)
                deallocate(vec)
                write(logfhandle,'(A,I3,A12,F6.1)')'>>> STAGE ',it,' WITH LP  : ',params%lp
                write(logfhandle,'(A,I6)'         )'>>> MINIBATCH SIZE: ', nint(real(nptcls) / real(iters_per_stage))
                write(logfhandle,'(A,I3,A1,I3)'   )'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',params%box_crop
                ! Starting reconstruction
                if( it == 1 )then
                    call cline_reconstruct3D%set('smpd_crop', params%smpd_crop)
                    call cline_reconstruct3D%set('box_crop',  params%box_crop)
                    call xreconstruct3D_distr%execute_shmem(cline_reconstruct3D)
                endif
                ! Minibatch size
                call cline_refine3D%set('batchfrac', params%batchfrac)
                ! Iterations per minibatch
                if( it == 1 )then
                    call cline_refine3D%set('maxits', (params%maxits_between+1)*iters_per_stage)
                else if( it == NSTAGES )then
                    call cline_refine3D%set('maxits', (params%maxits_between+1)*iters_per_stage)
                else
                    call cline_refine3D%set('maxits', params%maxits_between*iters_per_stage)
                endif
                ! Shift search & projection directions
                call cline_refine3D%set('trs',    trs(it))
                call cline_refine3D%set('nspace', nspaces(it))
                ! ICM filter
                if( params%l_icm .and. (it >= params%icm_stage) )then
                    call cline_refine3D%set('icm',   'yes')
                    call cline_refine3D%set('lambda', params%lambda)
                endif
                ! Execution
                call exec_refine3D(iter)
                ! Reconstruction
                if( it < NSTAGES )then
                    if( params%box_crop /= boxs(it+1) )then
                        call cline_reconstruct3D%set('smpd_crop',  smpds(it+1))
                        call cline_reconstruct3D%set('box_crop',   boxs(it+1))
                        call cline_reconstruct3D%set('which_iter', iter)
                        call xreconstruct3D_distr%execute_shmem(cline_reconstruct3D)
                    endif
                endif
            enddo
            deallocate(batches)
        case DEFAULT
            ! Constant number of particles
            deallocate(states)
            ! Main loop
            iter = 0
            do it = 1,NSTAGES
                write(logfhandle,'(A)') '>>>'
                ! Dimensions & resolution limit
                params%box_crop  = boxs(it)
                params%smpd_crop = smpds(it)
                params%lp        = lps(it)
                write(logfhandle,'(A,I3,A12,F6.1)')'>>> STAGE ',it,' WITH LP  : ',params%lp
                write(logfhandle,'(A,I6)'         )'>>> MINIBATCH SIZE: ', nint(real(nptcls_sel) / real(iters_per_stage))
                write(logfhandle,'(A,I3,A1,I3)'   )'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',params%box_crop
                ! Starting reconstruction
                if( it == 1 )then
                    call cline_reconstruct3D%set('smpd_crop', params%smpd_crop)
                    call cline_reconstruct3D%set('box_crop',  params%box_crop)
                    if( params%l_ml_reg .and. it == MLREG_ITER )then
                        call xcalc_pspec_distr%execute_shmem(cline_calc_pspec_distr)
                        call cline_reconstruct3D%set('which_iter', 1)
                        call cline_reconstruct3D%set('ml_reg',     'yes')
                        call cline_reconstruct3D%set('needs_sigma','yes')
                        call cline_reconstruct3D%set('objfun',     'euclid')
                        call cline_refine3D%set('continue', 'yes')
                    else
                        call cline_refine3D%set('vol1', 'recvol_state01.mrc')
                    endif
                    call xreconstruct3D_distr%execute_shmem(cline_reconstruct3D)
                    call spproj%read_segment('ptcl3D', params%projfile)
                    call spproj%os_ptcl3D%set_all2single('updatecnt',0.)
                    call spproj%write_segment_inside('ptcl3D', params%projfile)
                    call spproj%read_segment('out', params%projfile)
                    call spproj%add_vol2os_out('recvol_state01.mrc', params%smpd_crop, 1, 'vol')
                    call spproj%write_segment_inside('out', params%projfile)
                else
                    call spproj%read_segment('ptcl3D', params%projfile)
                    call spproj%os_ptcl3D%set_all2single('sampled',0.)
                    call spproj%write_segment_inside('ptcl3D', params%projfile)
                endif
                ! Minibatch size
                call cline_refine3D%set('batchfrac', params%batchfrac)
                ! Iterations per minibatch
                if( (it == 1) .or. (it == NSTAGES) )then
                    call cline_refine3D%set('maxits', (params%maxits_between+1)*iters_per_stage)
                else
                    call cline_refine3D%set('maxits', params%maxits_between*iters_per_stage)
                endif
                ! Shift search & projection directions
                call cline_refine3D%set('trs',    trs(it))
                call cline_refine3D%set('nspace', nspaces(it))
                ! ICM filter
                if( params%l_icm .and. (it >= params%icm_stage) )then
                    call cline_refine3D%set('icm',   'yes')
                    call cline_refine3D%set('lambda', params%lambda)
                endif
                ! ML regularization
                if( params%l_ml_reg .and. it >= MLREG_ITER ) call cline_refine3D%set('ml_reg', 'yes')
                ! Execution
                call exec_refine3D(iter)
                ! symmetrization
                if( it == SYMSEARCH_ITER-1 )then
                    call symmetrize
                    call cline_refine3D%set('pgrp', params%pgrp)
                    call cline_reconstruct3D%set('pgrp', params%pgrp)
                    l_srch4symaxis = .false.
                    l_symran       = .false.
                endif
                ! Reconstruction
                if( it < NSTAGES )then
                    call cline_reconstruct3D%set('smpd_crop',  smpds(it+1))
                    call cline_reconstruct3D%set('box_crop',   boxs(it+1))
                    call cline_reconstruct3D%set('which_iter', iter)
                    call cline_reconstruct3D%set('batchfrac',  params%batchfrac)
                    call xreconstruct3D_distr%execute_shmem(cline_reconstruct3D)
                    call spproj%read_segment('out', params%projfile)
                    call spproj%add_vol2os_out('recvol_state01.mrc', params%smpd_crop, 1, 'vol')
                    call spproj%write_segment_inside('out', params%projfile)
                    if( params%l_ml_reg .and. it >= MLREG_ITER )then
                        call cline_refine3D%set('continue', 'yes')
                    else
                        call cline_refine3D%set('vol1', 'recvol_state01.mrc')
                    endif
                endif
            enddo
        end select
        ! for visualization
        do state = 1, params%nstates
            str_state = int2str_pad(state,2)
            call final_vol%new([params%box_crop,params%box_crop,params%box_crop],params%smpd_crop)
            call final_vol%read(trim(VOL_FBODY)//trim(str_state)//trim(PPROC_SUFFIX)//trim(params%ext))
            call final_vol%generate_orthogonal_reprojs(reprojs)
            call reprojs%write_jpg('orthogonal_reprojs_state'//trim(str_state)//'.jpg')
            call final_vol%kill
            call reprojs%kill
        enddo
        ! Final reconstruction at original scale
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> RECONSTRUCTION AT ORIGINAL SAMPLING'
        write(logfhandle,'(A)') '>>>'
        call cline_reconstruct3D%set('ml_reg',      'no')
        call cline_reconstruct3D%set('needs_sigma', 'no')
        call cline_reconstruct3D%set('objfun',      'cc')
        call cline_reconstruct3D%delete('smpd_crop')
        call cline_reconstruct3D%delete('box_crop')
        call xreconstruct3D_distr%execute_shmem(cline_reconstruct3D)
        call spproj%read_segment('out',params%projfile)
        do state = 1, params%nstates
            str_state = int2str_pad(state,2)
            vol       = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
            call spproj%add_vol2os_out(vol, params%smpd, state, vol_type)
            if( trim(params%oritype).eq.'ptcl3D' )then
                call spproj%add_fsc2os_out(FSC_FBODY//str_state//trim(BIN_EXT), state, params%box)
            endif
        enddo
        call spproj%write_segment_inside('out',params%projfile)
        ! post-processing
        do state = 1, params%nstates
            call cline_postprocess%delete('lp')
            call cline_postprocess%set('state', real(state))
            call xpostprocess%execute(cline_postprocess)
        enddo
        do state = 1, params%nstates
            str_state      = int2str_pad(state,2)
            vol            = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
            vol_pproc      = add2fbody(vol,params%ext,PPROC_SUFFIX)
            vol_pproc_mirr = add2fbody(vol,params%ext,trim(PPROC_SUFFIX)//trim(MIRR_SUFFIX))
            if( file_exists(vol)            ) call simple_rename(vol,            trim(REC_FBODY)           //trim(str_state)//trim(params%ext))
            if( file_exists(vol_pproc)      ) call simple_rename(vol_pproc,      trim(REC_PPROC_FBODY)     //trim(str_state)//trim(params%ext))
            if( file_exists(vol_pproc_mirr) ) call simple_rename(vol_pproc_mirr, trim(REC_PPROC_MIRR_FBODY)//trim(str_state)//trim(params%ext))
        enddo
        if( trim(params%algorithm) == 'fast' )then
            ! restore states
            call spproj%read_segment('ptcl3D', params%projfile)
            call spproj%os_ptcl3D%set_all('state', real(states))
            call spproj%write_segment_inside('ptcl3D', params%projfile)
            deallocate(states)
        endif
        ! cleanup
        call spproj%kill
        call qsys_cleanup
        call simple_end('**** SIMPLE_BATCH_ABINITIO_3DMODEL NORMAL STOP ****')
        contains

            subroutine exec_refine3D( iter )
                integer,          intent(inout) :: iter
                character(len=:), allocatable   :: stage
                integer :: i, enditer
                call cline_refine3D%delete('endit')
                call cline_refine3D%set('smpd_crop', params%smpd_crop)
                call cline_refine3D%set('box_crop',  params%box_crop)
                call cline_refine3D%set('lp',        params%lp)
                call cline_refine3D%set('startit',   iter+1)
                call xrefine3D_distr%execute_shmem(cline_refine3D)
                call conv%read(l_err)
                enditer = nint(conv%get('iter'))
                ! probability table
                if( (it < NSTAGES) )then
                    if( nspaces(it+1) /= nspaces(it) )call del_file(trim(DIST_FBODY)//'.dat')
                endif
                call del_files(DIST_FBODY,      params_glob%nparts,ext='.dat')
                call del_files(ASSIGNMENT_FBODY,params_glob%nparts,ext='.dat')
                call del_file(trim(ASSIGNMENT_FBODY)//'.dat')
                if( it < NSTAGES )then
                    ! stash volumes
                    stage = '_stage_'//int2str(it)
                    do state = 1, params%nstates
                        str_state = int2str_pad(state,2)
                        vol       = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
                        vol_pproc = add2fbody(vol,params%ext,PPROC_SUFFIX)
                        if( file_exists(vol)      ) call simple_copy_file(vol,       add2fbody(vol,      params%ext,stage))
                        if( file_exists(vol_pproc)) call simple_copy_file(vol_pproc, add2fbody(vol_pproc,params%ext,stage))
                    enddo
                endif
                do i = iter+1,enditer-1
                    call del_file(trim(FSC_FBODY)//int2str_pad(1,2)//'_iter'//int2str_pad(i,3)//'.pdf')
                    call del_file('RESOLUTION_STATE'//int2str_pad(1,2)//'_ITER'//int2str_pad(i,3))
                enddo
                iter = enditer
            end subroutine exec_refine3D

            subroutine symmetrize()
                if( l_symran )then
                    call spproj%read_segment(params%oritype, params%projfile)
                    select case(trim(params%oritype))
                    case('cls3D')
                        call se1%symrandomize(spproj%os_cls3D)
                    case('ptcl3D')
                        call se1%symrandomize(spproj%os_ptcl3D)
                    end select
                    call spproj%write_segment_inside(params%oritype, params%projfile)
                endif
                if( l_srch4symaxis )then
                    write(logfhandle,'(A)') '>>>'
                    write(logfhandle,'(A)') '>>> SYMMETRY AXIS SEARCH'
                    write(logfhandle,'(A)') '>>>'
                    symlp = max(symlp, params%lp)
                    call cline_symsrch%set('lp',       symlp)
                    call cline_symsrch%set('box_crop', params%box_crop)
                    do state = 1,params%nstates
                        vol_str   = 'vol'//trim(int2str(state))
                        str_state = int2str_pad(state,2)
                        call cline_symsrch%set(vol_str, trim(VOL_FBODY)//trim(str_state)//trim(params%ext))
                    enddo
                    call xsymsrch%execute_shmem(cline_symsrch)
                    call del_file('SYMAXIS_SEARCH_FINISHED')
                endif
            end subroutine symmetrize

    end subroutine exec_batch_abinitio_3Dmodel

    !> for generation of an initial 3d model from particles
    subroutine exec_abinitio_3Dmodel2( self, cline )
        use simple_convergence, only: convergence
        use simple_fsc,         only: plot_fsc
        class(abinitio_3Dmodel2_commander), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        real,    parameter :: SCALEFAC        = 0.667
        real,    parameter :: CENLP_DEFAULT   = 30.
        real,    parameter :: LP_DEFAULT      = 6.
        real,    parameter :: LPSTART_DEFAULT = 30., LPSTOP_DEFAULT=LP_DEFAULT
        integer, parameter :: NPARTS  = 4
        integer, parameter :: MINBOX  = 64
        integer, parameter :: NSTAGES_DEFAULT = 22
        integer, parameter :: MAXITS_SHORT = 5
        integer, parameter :: NSPACE1 = 500, NSPACE2 = 1000, NSPACE3 = 2000
        integer, parameter :: SYMSEARCH_DEFAULT = 5
        integer, parameter :: MLREG_ITER        = 1
        integer, parameter :: SHIFT_STAGE_DEFAULT = NSTAGES_DEFAULT-5 ! in [1;NSTAGES+1]
        ! commanders
        type(refine3D_commander_distr)      :: xrefine3D_distr
        type(reconstruct3D_commander_distr) :: xreconstruct3D_distr
        type(postprocess_commander)         :: xpostprocess
        type(symaxis_search_commander)      :: xsymsrch
        type(calc_pspec_commander_distr)    :: xcalc_pspec_distr
        ! command lines
        type(cmdline)                 :: cline_refine3D, cline_reconstruct3D, cline_reconstruct3D_mlreg
        type(cmdline)                 :: cline_postprocess, cline_symsrch, cline_calc_pspec_distr
        ! other
        type(parameters)              :: params
        type(sp_project)              :: spproj, spproj_part
        type(convergence)             :: conv
        type(sym)                     :: se1, se2
        type(class_frcs)              :: clsfrcs
        type(image)                   :: vol_even, vol_odd, reprojs, tmpvol, vol
        type(qsys_env)                :: qenv
        real,             allocatable :: fsc(:), res(:)
        character(len=:), allocatable :: str_state, vol_pproc, vol_pproc_mirr
        character(len=:), allocatable :: stack_name, dir, fsc_fname
        integer,          allocatable :: states(:), tmp(:), iters(:), prev_iters(:)
        character(len=STDLEN), allocatable :: completion_fnames(:)
        character(len=LONGSTRLEN)     :: vol_str
        real    :: lps(NSTAGES_DEFAULT), smpds(NSTAGES_DEFAULT), trs(NSTAGES_DEFAULT)
        integer :: boxs(NSTAGES_DEFAULT)
        real    :: smpd_target, lp_target, scale, trslim, cenlp, symlp, dummy, msk
        integer :: it, prev_box_crop, maxits, nptcls_sel, filtsz
        integer :: nstages, symsearch_iter, istk, part, iter, nptcls_part, i,j, cnt
        logical :: l_autoscale, l_lpset, l_err, l_srch4symaxis, l_symran, l_sym, l_lpstop_set
        logical :: l_lpstart_set
        call cline%set('oritype',      'ptcl3D')
        call cline%set('ml_reg',       'yes')
        call cline%set('stoch_update', 'no')
        call cline%set('icm',          'no')
        if( .not. cline%defined('mkdir')        ) call cline%set('mkdir',        'yes')
        if( .not. cline%defined('refine')       ) call cline%set('refine',      'prob')
        if( .not. cline%defined('autoscale')    ) call cline%set('autoscale',    'yes')
        if( .not. cline%defined('sigma_est')    ) call cline%set('sigma_est', 'global')
        if( .not. cline%defined('prob_sh')      ) call cline%set('prob_sh',      'yes')
        if( .not. cline%defined('prob_athres')  ) call cline%set('prob_athres',    10.)
        if( .not. cline%defined('center')       ) call cline%set('center',        'no')
        if( .not. cline%defined('objfun')       ) call cline%set('objfun',    'euclid')
        if( .not. cline%defined('oritype')      ) call cline%set('oritype',   'ptcl3D')
        if( .not. cline%defined('pgrp')         ) call cline%set('pgrp',          'c1')
        if( .not. cline%defined('pgrp_start')   ) call cline%set('pgrp_start',    'c1')
        if( .not. cline%defined('shift_stage')  ) call cline%set('shift_stage', SHIFT_STAGE_DEFAULT)
        if( .not. cline%defined('ptclw')        ) call cline%set('ptclw',         'no')
        if( .not. cline%defined('nparts')       ) call cline%set('nparts',      NPARTS)
        ! resolution limit strategy
        l_lpset       = .false.
        l_lpstop_set  = cline%defined('lpstop')
        l_lpstart_set = cline%defined('lpstart')
        if( cline%defined('lp') )then
            if( l_lpstart_set .or. l_lpstop_set )then
                THROW_HARD('One of LP or LPSTART & LPSTOP must be defined!')
            endif
            l_lpset = .true.
        else
            if( .not.l_lpstart_set ) call cline%set('lpstart',LPSTART_DEFAULT)
            if( .not.l_lpstop_set  ) call cline%set('lpstop', LPSTOP_DEFAULT)
        endif
        ! make master parameters
        if( cline%defined('update_frac') ) call cline%delete('stoch_update')
        call params%new(cline)
        call cline%set('mkdir', 'no')
        call cline%delete('autoscale')
        call cline%delete('lpstart')
        call cline%delete('lpstop')
        call cline%delete('lp')
        call cline%delete('shift_stage')
        allocate(completion_fnames(params%nparts),iters(params%nparts),prev_iters(params%nparts))
        call qenv%new(1)
        str_state = int2str_pad(1,2)
        ! stages specific parameters
        nstages        = NSTAGES_DEFAULT
        symsearch_iter = SYMSEARCH_DEFAULT
        if( l_lpset )then
            params%lpstop  = params%lp
            params%lpstart = params%lp
        endif
        if( params%shift_stage < 1 .or. params%shift_stage > nstages+1 )then
            params%shift_stage = min(nstages+1,max(1,params%shift_stage))
            THROW_WARN('SHIFT_STAGE out of range, defaulting to: '//int2str(params%shift_stage))
        endif
        ! read project
        call spproj%read(params%projfile)
        call spproj%update_projinfo(cline)
        call spproj%write_segment_inside('projinfo', params%projfile)
        if( .not. cline%defined('vol1') )then
            ! randomize projection directions
            if( spproj%os_ptcl3D%get_noris() < 1 )then
                THROW_HARD('Particles could not be found in the project')
            endif
            call spproj%os_ptcl3D%rnd_oris
            call spproj%os_ptcl3D%set_all2single('w',1.)
            states = nint(spproj%os_ptcl3D%get_all('state'))
            call spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        write(logfhandle,'(A,F5.1)') '>>> STARTING RESOLUTION LIMIT (IN A): ', params%lpstart
        write(logfhandle,'(A,F5.1)') '>>> HARD     RESOLUTION LIMIT (IN A): ', params%lpstop
        if( trim(params%center).eq.'yes' )then
            write(logfhandle,'(A,F5.1)') '>>> CENTERING  LOW-PASS LIMIT (IN A): ', params%cenlp
        endif
        ! centering & symmetry resolution limit
        call mskdiam2lplimits(params%mskdiam, symlp, dummy, cenlp)
        if( .not. cline%defined('cenlp') )then
            params%cenlp = cenlp
            call cline%set('cenlp', params%cenlp)
        endif
        ! symmetry
        if( l_lpset )then
            ! from mskdiam2lplimits lpstart above
        else
            symlp = (params%lpstart+params%lpstop)/2.
        endif
        l_srch4symaxis = trim(params%pgrp) .ne. trim(params%pgrp_start)
        l_symran       = .false.
        l_sym          = l_srch4symaxis
        if( params%pgrp_start.ne.'c1' .or. params%pgrp.ne.'c1' )then
            se1 = sym(params%pgrp_start)
            se2 = sym(params%pgrp)
            if(se1%get_nsym() > se2%get_nsym())then
                ! ensure se2 is a subgroup of se1
                if( .not. se1%has_subgrp(params%pgrp) )THROW_HARD('Incompatible symmetry groups; exec_abinitio_3Dmodel')
                ! set flag for symmetry randomisation
                ! in case we are moving from a higher to lower group
                l_symran = .true.
            else if( se2%get_nsym() > se1%get_nsym() )then
                ! ensure se1 is a subgroup of se2
                if( .not. se2%has_subgrp(params%pgrp_start) )THROW_HARD('Incompatible symmetry groups; exec_abinitio_3Dmodel')
            endif
        endif
        ! dimensions defaults
        params%box       = spproj%get_box()
        params%smpd_crop = params%smpd
        params%box_crop  = params%box
        l_autoscale      = .false.
        ! command-lines
        cline_refine3D            = cline
        cline_reconstruct3D       = cline
        cline_postprocess         = cline
        cline_symsrch             = cline
        cline_reconstruct3D_mlreg = cline_reconstruct3D
        cline_calc_pspec_distr    = cline
        call cline_refine3D%set('prg',                'refine3D')
        call cline_refine3D%set('projfile',      params%projfile)
        call cline_refine3D%set('pgrp',        params%pgrp_start)
        call cline_refine3D%set('maxits',                    999)
        call cline_refine3D%delete('nparts')
        call cline_reconstruct3D%set('prg',      'reconstruct3D')
        call cline_reconstruct3D%set('box',     real(params%box))
        call cline_reconstruct3D%set('projfile', params%projfile)
        call cline_reconstruct3D%set('ml_reg',              'no')
        call cline_reconstruct3D%set('needs_sigma',         'no')
        call cline_reconstruct3D%set('objfun',              'cc')
        call cline_reconstruct3D%set('pgrp',   params%pgrp_start)
        call cline_postprocess%set('prg',          'postprocess')
        call cline_postprocess%set('projfile',   params%projfile)
        call cline_postprocess%set('imgkind',              'vol')
        if( l_srch4symaxis )then
            call cline_symsrch%set('prg',     'symaxis_search') ! needed for cluster exec
            call cline_symsrch%set('pgrp',     params%pgrp)
            call cline_symsrch%set('projfile', params%projfile)
            call cline_symsrch%set('hp',       params%hp)
            call cline_symsrch%set('center',   'yes')
        endif
        call cline_reconstruct3D_mlreg%set('prg',         'reconstruct3D')
        call cline_reconstruct3D_mlreg%set('objfun',      'euclid')
        call cline_reconstruct3D_mlreg%set('needs_sigma', 'yes')
        call cline_reconstruct3D_mlreg%set('sigma_est',   params%sigma_est)
        call cline_reconstruct3D_mlreg%set('ml_reg',      'yes')
        call cline_calc_pspec_distr%set('prg',      'calc_pspec')
        ! Frequency marching plan
        lps(1) = params%lpstart
        do it = 2,nstages-1
            lps(it) = params%lpstop + (params%lpstart-params%lpstop) * real(nstages-it) / real(nstages)
        enddo
        lps(nstages) = params%lpstop
        if( l_lpset )then
            ! from mskdiam2lplimits lpstart above
        else
            symlp = (params%lpstart+params%lpstop)/2.
        endif
        ! dimensions
        do it = 1,nstages
            lp_target   = lps(it) * SCALEFAC
            smpd_target = max(params%smpd, lp_target/2.)
            call autoscale(params%box, params%smpd, smpd_target, boxs(it), smpds(it), scale, minbox=MINBOX)
            if( it < params%shift_stage )then
                trs(it) = 0.
            else
                trs(it) = max(2.0, AHELIX_WIDTH / smpds(it) / 2.0)
            endif
        enddo
        ! random reconstruction
        params%smpd_crop = smpds(1)
        params%box_crop  = boxs(1)
        call cline_reconstruct3D%set('smpd_crop', params%smpd_crop)
        call cline_reconstruct3D%set('box_crop',  params%box_crop)
        if( params%l_ml_reg .and. MLREG_ITER==1 )then
            call xcalc_pspec_distr%execute_shmem(cline_calc_pspec_distr)
            call cline_reconstruct3D%set('which_iter', 1)
            call cline_reconstruct3D%set('ml_reg',     'yes')
            call cline_reconstruct3D%set('needs_sigma','yes')
            call cline_reconstruct3D%set('objfun',     'euclid')
        endif
        call xreconstruct3D_distr%execute_shmem(cline_reconstruct3D)
        call spproj%read_segment('ptcl3D', params%projfile)
        call spproj%os_ptcl3D%set_all2single('updatecnt',0.)
        call spproj%write_segment_inside('ptcl3D', params%projfile)
        call spproj%read_segment('out', params%projfile)
        call spproj%add_vol2os_out('recvol_state01.mrc', params%smpd_crop, 1, 'vol')
        call spproj%write_segment_inside('out', params%projfile)
        ! updating stack names to absolute path
        call spproj%read_segment('stk', params%projfile)
        do istk = 1,spproj%os_stk%get_noris()
            stack_name = trim(spproj%get_stkname(istk))
            stack_name = simple_abspath(stack_name, check_exists=.false.)
            call spproj%os_stk%set(istk, 'stk', stack_name)
        enddo
        call spproj%write_segment_inside('stk', params%projfile)
        call spproj%read_segment('ptcl2D', params%projfile)
        call spproj%read_segment('ptcl3D', params%projfile)
        ! directory structure
        do part = 1,params%nparts
            dir = int2str(part)//'/'
            call simple_mkdir(dir)
        enddo
        ! Parts partitioning
        call cline_refine3D%delete('projfile')
        call cline_refine3D%set('projname', get_fbody(basename(params%projfile), 'simple'))
        call cline_refine3D%set('projfile', basename(params%projfile))
        nptcls_sel  = count(states==1)
        nptcls_part = ceiling(real(nptcls_sel)/real(params%nparts))
        tmp = states
        j   = 0
        do part = 1,params%nparts
            spproj_part%os_stk    = spproj%os_stk
            spproj_part%os_ptcl2D = spproj%os_ptcl2D
            spproj_part%os_ptcl3D = spproj%os_ptcl3D
            spproj_part%projinfo  = spproj%projinfo
            spproj_part%compenv   = spproj%compenv
            tmp = states
            if( j > 0 ) tmp(1:j) = 0
            cnt = 0
            do i = j+1,params%nptcls
                if( states(i) == 1)then
                    cnt       = cnt+1
                    states(i) = part
                    if( cnt == nptcls_part )then
                        j = i
                        exit
                    endif
                endif
            enddo
            if( part < params%nparts ) tmp(j+1:) = 0
            call spproj_part%os_ptcl2D%set_all('state', real(tmp))
            call spproj_part%os_ptcl3D%set_all('state', real(tmp))
            call spproj_part%prune_particles
            call chdir(int2str(part)//'/')
            call spproj_part%update_projinfo(cline_refine3D)
            call spproj_part%write
            call chdir('../')
            completion_fnames(part) = int2str(part)//'/'//trim(JOB_FINISHED_FBODY)
            call spproj_part%kill
        enddo
        deallocate(tmp)
        ! Stages loop
        iters(:) = 0
        do it = 1,nstages
            params%smpd_crop = smpds(it)
            params%box_crop  = boxs(it)
            params%lp        = lps(it)
            params%trs       = trs(it)
            write(logfhandle,'(A,I3,A9,F5.1)')'>>> STAGE ',it,' WITH LP =',params%lp
            write(logfhandle,'(A,I3)')        '>>> CROPPED IMAGE SIZE: ',params%box_crop
            call cline_refine3D%set('smpd_crop', params%smpd_crop)
            call cline_refine3D%set('box_crop',  params%box_crop)
            call cline_refine3D%set('lp',        params%lp)
            if( it == 1 )then
                call cline_refine3D%set('vol1', '../recvol_state01.mrc')
            else
                call cline_refine3D%set('vol1', '../recvol_state01_stage'//int2str_pad(it-1,2)//'.mrc')
            endif
            ! # of iterations
            call cline_refine3D%set('maxits', MAXITS_SHORT)
            ! projection directions & shift
            if( it < params%shift_stage )then
                call cline_refine3D%set('nspace', NSPACE1)
            else
                call cline_refine3D%set('nspace', NSPACE2)
            end if
            if( it >= nstages-2 ) call cline_refine3D%set('nspace', NSPACE3)
            call cline_refine3D%set('trs', params%trs)
            ! execution
            do part = 1,params%nparts
                call exec_refine3D(part)
            enddo
            ! waiting
            call qsys_watcher(completion_fnames)
            ! convergence, volume averaging, padding & cleanup
            prev_iters = iters
            call vol%new([params%box_crop, params%box_crop,params%box_crop],params%smpd_crop)
            call vol_even%new([params%box_crop, params%box_crop,params%box_crop],params%smpd_crop)
            call vol_odd%new([params%box_crop, params%box_crop,params%box_crop],params%smpd_crop)
            call tmpvol%new([params%box_crop, params%box_crop,params%box_crop],params%smpd_crop)
            do part = 1,params%nparts
                call chdir(int2str(part)//'/')
                ! convergence parameters
                call conv%read(l_err)
                iters(part) = nint(conv%get('iter'))
                write(logfhandle,'(A,I3,A,F7.3,A,F7.3)')'>>> PART ',part,'; PROJ OVERLAP: ',&
                    &conv%get('mi_proj'),'; SCORE: ',conv%get('score')
                ! volumes
                if( part == 1)then
                    call vol%read('recvol_state01_iter'//int2str_pad(iters(part),3)//'.mrc')
                    call vol_even%read('recvol_state01_iter'//int2str_pad(iters(part),3)//'_even.mrc')
                    call vol_odd%read('recvol_state01_iter'//int2str_pad(iters(part),3)//'_odd.mrc')
                else
                    call tmpvol%read('recvol_state01_iter'//int2str_pad(iters(part),3)//'.mrc')
                    call vol%add(tmpvol)
                    call tmpvol%read('recvol_state01_iter'//int2str_pad(iters(part),3)//'_even.mrc')
                    call vol_even%add(tmpvol)
                    call tmpvol%read('recvol_state01_iter'//int2str_pad(iters(part),3)//'_odd.mrc')
                    call vol_odd%add(tmpvol)
                endif
                ! cleanup
                call qsys_cleanup
                call del_files(DIST_FBODY,      1,ext='.dat')
                call del_files(ASSIGNMENT_FBODY,1,ext='.dat')
                call del_file(trim(DIST_FBODY)      //'.dat')
                call del_file(trim(ASSIGNMENT_FBODY)//'.dat')
                call del_file(JOB_FINISHED_FBODY)
                call del_file(trim(FSC_FBODY)//int2str_pad(1,2)//BIN_EXT)
                do i = prev_iters(part)+1,iters(part)-1
                    call del_file(trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(i,3)//PPROC_SUFFIX//params%ext)
                    call del_file(trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(i,3)//LP_SUFFIX//params%ext)
                    call del_file(trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(i,3)//'_even'//params%ext)
                    call del_file(trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(i,3)//'_odd'//params%ext)
                    call del_file(trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(i,3)//'_even_unfil'//params%ext)
                    call del_file(trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(i,3)//'_odd_unfil'//params%ext)
                    call del_file(trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(i,3)//params%ext)
                enddo
                call chdir('../')
            enddo
            ! averaging & fsc
            if( it < NSTAGES_DEFAULT )then
                ! Volume & FSC will be padded on the fly at the next refine3D run
                call vol%div(real(params%nparts))
                call vol_even%div(real(params%nparts))
                call vol_odd%div(real(params%nparts))
                call vol%write(trim(VOL_FBODY)//trim(str_state)//'_stage'//int2str_pad(it,2)//'.mrc')
                filtsz = fdim(params%box_crop) - 1
                msk    = real(params%box_crop / 2) - COSMSKHALFWIDTH - 1.
                allocate(fsc(filtsz),source=0.)
                call vol_even%mask(msk, 'soft', backgr=0.)
                call vol_odd%mask(msk, 'soft', backgr=0.)
                call vol_even%fft()
                call vol_odd%fft()
                call vol_even%fsc(vol_odd, fsc)
                fsc_fname = trim(FSC_FBODY)//int2str_pad(1,2)//BIN_EXT
                call arr2file(fsc, fsc_fname)
                call cline_refine3D%set('fsc', '../'//trim(fsc_fname))
                res = get_resarr(params%box_crop, params%smpd_crop)
                call plot_fsc(size(fsc), fsc, res, params%smpd_crop, 'fsc_stage_'//int2str_pad(it,2))
                deallocate(fsc,res)
                call vol%kill
                call vol_even%kill
                call vol_odd%kill
                call tmpvol%kill
            endif
            ! symmetrization
            if( it == SYMSEARCH_ITER-1 )then
                call consolidate_alnparms
                call cline_symsrch%set('vol1', trim(VOL_FBODY)//trim(str_state)//'_stage'//int2str_pad(it,2)//'.mrc')
                call symmetrize
                call cline_refine3D%set('pgrp', params%pgrp)
                call cline_reconstruct3D%set('pgrp', params%pgrp)
                l_srch4symaxis = .false.
                l_symran       = .false.
                ! transfer symmetrized parameters
                call spproj%read_segment('ptcl3D',params%projfile)
                do part = 1,params%nparts
                    call spproj_part%read_segment('ptcl3D', int2str(part)//'/'//basename(params%projfile))
                    j = 0
                    do i = 1,params%nptcls
                        if( states(i) /= part ) cycle
                        j = j+1
                        call spproj_part%os_ptcl3D%transfer_3Dparams(j, spproj%os_ptcl3D, i)
                    enddo
                    call spproj_part%write_segment_inside('ptcl3D',int2str(part)//'/'//basename(params%projfile))
                enddo
                call spproj_part%kill
            endif
        enddo
        ! gathering alignment parameters
        call consolidate_alnparms
        ! final reconstruction
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> RECONSTRUCTION AT ORIGINAL SAMPLING'
        write(logfhandle,'(A)') '>>>'
        ! no ML-filtering
        call cline_reconstruct3D%set('ml_reg',      'no')
        call cline_reconstruct3D%set('needs_sigma', 'no')
        call cline_reconstruct3D%set('objfun',      'cc')
        call cline_reconstruct3D%delete('smpd_crop')
        call cline_reconstruct3D%delete('box_crop')
        call xreconstruct3D_distr%execute_shmem(cline_reconstruct3D)
        vol_str = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
        call spproj%read_segment('out',params%projfile)
        call spproj%add_vol2os_out(vol_str, params%smpd, 1, 'vol')
        call spproj%add_fsc2os_out(FSC_FBODY//str_state//trim(BIN_EXT), 1, params%box)
        call spproj%write_segment_inside('out',params%projfile)
        ! post-processing
        call cline_postprocess%delete('lp')
        call cline_postprocess%set('state', 1)
        call xpostprocess%execute(cline_postprocess)
        vol_pproc      = add2fbody(vol_str,params%ext,PPROC_SUFFIX)
        vol_pproc_mirr = add2fbody(vol_str,params%ext,trim(PPROC_SUFFIX)//trim(MIRR_SUFFIX))
        call simple_rename(vol_str, trim(REC_FBODY)//trim(str_state)//trim(params%ext))
        if(file_exists(vol_pproc)     ) call simple_rename(vol_pproc,      trim(REC_PPROC_FBODY)     //trim(str_state)//trim(params%ext))
        if(file_exists(vol_pproc_mirr)) call simple_rename(vol_pproc_mirr, trim(REC_PPROC_MIRR_FBODY)//trim(str_state)//trim(params%ext))
        ! cleanup
        call spproj%kill
        call qsys_cleanup
        call simple_end('**** SIMPLE_ABINITIO_3DMODEL2 NORMAL STOP ****')
        contains

            subroutine exec_refine3D( part )
                integer,          intent(in)  :: part
                character(len=XLONGSTRLEN) :: cwd
                call cline_refine3D%set('startit', iters(part)+1)
                dir = int2str(part)//'/'
                call chdir(dir)
                call simple_getcwd(cwd)
                cwd_glob = trim(cwd)
                call qenv%new(1)
                call qenv%exec_simple_prg_in_queue_async(cline_refine3D, './refine3D', 'log_refine3D')
                call chdir('../')
                call simple_getcwd(cwd_glob)
            end subroutine exec_refine3D

            subroutine symmetrize()
                use simple_projector_hlev, only: rotvol_slim
                use simple_projector,      only: projector
                type(projector) :: vol_pad
                type(image)     :: rovol_pad, rovol
                type(ori)       :: o
                real    :: symaxis_rmat(3,3), symop_rmat(3,3), rmat(3,3)
                integer :: ldim_pd(3), boxpd,isym, nsym
                if( l_symran )then
                    call spproj%read_segment(params%oritype, params%projfile)
                    call se1%symrandomize(spproj%os_ptcl3D)
                    call spproj%write_segment_inside(params%oritype, params%projfile)
                endif
                if( l_srch4symaxis )then
                    write(logfhandle,'(A)') '>>>'
                    write(logfhandle,'(A)') '>>> SYMMETRY AXIS SEARCH'
                    write(logfhandle,'(A)') '>>>'
                    symlp = max(symlp, params%lp)
                    call cline_symsrch%set('lp',       symlp)
                    call cline_symsrch%set('box_crop', params%box_crop)
                    call xsymsrch%execute_shmem(cline_symsrch)
                    call del_file('SYMAXIS_SEARCH_FINISHED')
                    ! symmetrize volume
                    call vol%new([params%box_crop, params%box_crop,params%box_crop],params%smpd_crop)
                    call rovol%new([params%box_crop, params%box_crop,params%box_crop],params%smpd_crop)
                    boxpd   = 2 * round2even(KBALPHA * real(params%box_crop))
                    ldim_pd = [boxpd,boxpd,boxpd]
                    call rovol_pad%new(ldim_pd, params%smpd_crop)
                    call vol_pad%new(ldim_pd, params%smpd_crop)
                    call vol%read('vol_aligned2_'//trim(params%pgrp)//'axis'//params%ext)
                    call vol%pad(vol_pad)
                    call vol_pad%fft
                    call vol_pad%expand_cmat(KBALPHA)
                    nsym = se2%get_nsym()
                    do isym =2,nsym
                        call se2%get_symori(isym, o)
                        call rotvol_slim(vol_pad, rovol_pad, rovol, o)
                        call vol%add_workshare(rovol)
                    end do
                    call vol%div(real(nsym))
                    call vol%write(trim(VOL_FBODY)//trim(str_state)//'_stage'//int2str_pad(it,2)//'.mrc')
                    call o%kill
                    call rovol%kill
                    call rovol_pad%kill
                    call vol%kill
                    call vol_pad%kill
                endif
            end subroutine symmetrize

            subroutine consolidate_alnparms
                integer :: i,j,part
                do part = 1,params%nparts
                    call spproj_part%read_segment('ptcl3D', int2str(part)//'/'//basename(params%projfile))
                    j = 0
                    do i = 1,params%nptcls
                        if( states(i) /= part ) cycle
                        j = j+1
                        call spproj%os_ptcl3D%transfer_3Dparams(i, spproj_part%os_ptcl3D, j)
                    enddo
                enddo
                call spproj_part%kill
                call spproj%write_segment_inside('ptcl3D',params%projfile)
            end subroutine consolidate_alnparms

    end subroutine exec_abinitio_3Dmodel2

end module simple_commander_abinitio
