! concrete commander: high-level workflows
module simple_commanders_abinitio
include 'simple_lib.f08'
use simple_stack_io,            only: stack_io
use simple_sp_project,          only: sp_project
use simple_qsys_env,            only: qsys_env
use simple_procimgstk,          only: shift_imgfile
use simple_parameters,          only: parameters, params_glob
use simple_image,               only: image
use simple_convergence,         only: convergence
use simple_commanders_volops,   only: commander_reproject, commander_symaxis_search, commander_postprocess, commander_symmetrize_map
use simple_commanders_refine3D, only: commander_refine3D, commander_refine3D_distr
use simple_commanders_rec,      only: commander_reconstruct3D, commander_reconstruct3D_distr
use simple_commanders_mask,     only: commander_automask
use simple_commander_base,      only: commander_base
use simple_cmdline,             only: cmdline
use simple_cluster_seed,        only: gen_labelling
use simple_class_frcs,          only: class_frcs
use simple_qsys_funs
use simple_nice
use simple_euclid_sigma2
use simple_decay_funs
use simple_commanders_euclid
implicit none

public :: commander_abinitio3D_cavgs, commander_abinitio3D_cavgs_fast, commander_abinitio3D, commander_multivol_assign
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_abinitio3D_cavgs
    contains
    procedure :: execute => exec_abinitio3D_cavgs
end type commander_abinitio3D_cavgs

type, extends(commander_base) :: commander_abinitio3D_cavgs_fast
    contains
    procedure :: execute => exec_abinitio3D_cavgs_fast
end type commander_abinitio3D_cavgs_fast

type, extends(commander_base) :: commander_abinitio3D
    contains
    procedure :: execute => exec_abinitio3D
end type commander_abinitio3D

type, extends(commander_base) :: commander_multivol_assign
    contains
    procedure :: execute => exec_multivol_assign
end type commander_multivol_assign

! class constants
character(len=*), parameter :: REC_FBODY             = 'rec_final_state'
character(len=*), parameter :: STR_STATE_GLOB        = '01'
real,             parameter :: LPSTOP_BOUNDS(2)      = [4.5,6.0]
real,             parameter :: LPSTART_BOUNDS(2)     = [10.,20.]
real,             parameter :: CENLP_DEFAULT         = 30.
real,             parameter :: LPSYMSRCH_LB          = 12.
real,             parameter :: UPDATE_FRAC_MAX       = 0.9               !< to ensure fractional update is always on
real,             parameter :: UPDATE_FRAC_MIN       = 0.1               !< 10% of the particles updated each iteration
integer,          parameter :: NSTAGES               = 8
integer,          parameter :: NSTAGES_INI3D         = 4 ! # of ini3D stages used for initialization
integer,          parameter :: NSTAGES_INI3D_MAX     = 7
integer,          parameter :: PHASES(3)             = [2,6,NSTAGES]
integer,          parameter :: MAXITS(8)             = [20,20,17,17,17,17,15,30]
integer,          parameter :: MAXITS_GLOB           = SUM(MAXITS(1:7))  ! the last 30 iterations are not included in this estimate since the sampling method changes
integer,          parameter :: NSPACE(3)             = [500,1000,2500]
integer,          parameter :: SYMSRCH_STAGE         = 3
integer,          parameter :: PROBREFINE_STAGE      = 5
integer,          parameter :: ICM_STAGE             = PROBREFINE_STAGE     ! we switch from ML regularization when prob is switched on
integer,          parameter :: STOCH_SAMPL_STAGE     = PROBREFINE_STAGE     ! we switch from greedy to stochastic blanced class sampling when prob is switched on
integer,          parameter :: TRAILREC_STAGE_SINGLE = STOCH_SAMPL_STAGE    ! we start trailing when we start sampling particles randomly
integer,          parameter :: TRAILREC_STAGE_MULTI  = NSTAGES              ! we start trailing in the last stage
integer,          parameter :: LPAUTO_STAGE          = NSTAGES - 1          ! cannot be switched on too early
integer,          parameter :: RECALC_STARTREC_STAGE = LPAUTO_STAGE         ! re-estimate starting volume for optimal LP and AUTOMSK
integer,          parameter :: AUTOMSK_STAGE         = LPAUTO_STAGE         ! swith on automasking when lpauto is switched on
integer,          parameter :: HET_DOCKED_STAGE      = NSTAGES              ! stage at which state splitting is done when multivol_mode==docked
integer,          parameter :: STREAM_ANALYSIS_STAGE = 5                    ! when streaming on some analysis will be performed
integer,          parameter :: CAVGWEIGHTS_STAGE     = 3                    ! when to activate optional cavg weighing in abinitio3D_cavgs/cavgs_fast
integer,          parameter :: GAUREF_LAST_STAGE     = PHASES(1)            ! When to stop using gaussian filtering of the references with polar=yes
integer,          parameter :: NSPACE_PHASE_POLAR(3) = [  2, PROBREFINE_STAGE, NSTAGES]
integer,          parameter :: NSPACE_POLAR(3)       = [500,             1000,    1500]
integer,          parameter :: MAXITS_BETWEEN        = 10                   ! Development

! class variables
type(lp_crop_inf), allocatable :: lpinfo(:)
logical          :: l_srch4symaxis=.false., l_symran=.false., l_sym=.false., l_update_frac_dyn=.false., l_polar=.false.
logical          :: l_ini3D=.false., l_lpauto=.false., l_nsample_given=.false., l_nsample_stop_given=.false., l_automsk=.false.
type(sym)        :: se1, se2
type(cmdline)    :: cline_refine3D, cline_symmap, cline_reconstruct3D, cline_postprocess, cline_reproject
real             :: update_frac  = 1.0, update_frac_dyn  = 1.0
integer          :: nstates_glob = 1, nptcls_eff = 0, nsample_minmax(2), maxits_dyn=0

contains

    !> for generation of an initial 3D model from class averages
    subroutine exec_abinitio3D_cavgs( self, cline )
        use simple_estimate_ssnr, only: lpstages_fast
        class(commander_abinitio3D_cavgs), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        ! shared-mem commanders
        type(commander_refine3D)      :: xrefine3D
        type(commander_reconstruct3D) :: xreconstruct3D
        type(commander_reproject)     :: xreproject
        ! other
        type(string)                  :: stk, orig_stk, shifted_stk, stk_even, stk_odd, ext
        integer,          allocatable :: states(:)
        type(ori)                     :: o, o_even, o_odd
        type(parameters)              :: params
        type(ctfparams)               :: ctfvars
        type(sp_project)              :: spproj, work_proj
        type(image)                   :: img
        type(stack_io)                :: stkio_r, stkio_r2, stkio_w
        type(string)                  :: final_vol, work_projfile
        integer                       :: icls, ncavgs, cnt, even_ind, odd_ind, istage, nstages_ini3D, s
        if( cline%defined('nparts') ) THROW_HARD('abinitio3D_cavgs does not support distributed execution, remove nparts from command line')
        call cline%set('sigma_est', 'global') ! obviously
        call cline%set('oritype',      'out') ! because cavgs are part of out segment
        call cline%set('bfac',            0.) ! because initial models should not be sharpened
        if( .not. cline%defined('mkdir')       ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('objfun')      ) call cline%set('objfun',  'euclid') ! use noise normalized Euclidean distances from the start
        if( .not. cline%defined('overlap')     ) call cline%set('overlap',     0.95)
        if( .not. cline%defined('prob_athres') ) call cline%set('prob_athres',  90.) ! reduces # failed runs on trpv1 from 4->2/10
        if( .not. cline%defined('cenlp')       ) call cline%set('cenlp', CENLP_DEFAULT)
        if( .not. cline%defined('imgkind')     ) call cline%set('imgkind',   'cavg')
        if( .not. cline%defined('lp_auto')     ) call cline%set('lp_auto',    'yes')
        if( .not. cline%defined('noise_norm')  ) call cline%set('noise_norm',  'no')
        if( .not. cline%defined('cavgw')       ) call cline%set('cavgw',       'no')
        if( .not. cline%defined('lpstart')     ) call cline%set('lpstart',      20.)
        if( .not. cline%defined('lpstop')      ) call cline%set('lpstop',        8.)
        ! make master parameters
        call params%new(cline)
        call cline%set('mkdir',       'no')   ! to avoid nested directory structure
        call cline%set('oritype', 'ptcl3D')   ! from now on we are in the ptcl3D segment, final report is in the cls3D segment
        ! set work projfile
        work_projfile = 'abinitio3D_cavgs_tmpproj.simple'
        ! set class global lp_auto flag for low-pass limit estimation
        l_lpauto = .true.
        if( cline%defined('lp_auto') ) l_lpauto = params%l_lpauto
        ! set nstages_ini3D
        nstages_ini3D = NSTAGES_INI3D_MAX
        if( cline%defined('nstages') )then
            nstages_ini3D = min(NSTAGES_INI3D_MAX,params%nstages)   
        endif
        ! prepare class command lines
        call prep_class_command_lines(cline, work_projfile)
        ! set symmetry class variables
        call set_symmetry_class_vars
        ! read project
        call spproj%read(params%projfile)
        ! set low-pass limits and downscaling info from FRCs
        if( cline%defined('lpstart_ini3D').or.cline%defined('lpstop_ini3D') )then
            ! overrides resolution limits scheme based on frcs
            if( cline%defined('lpstart_ini3D').and.cline%defined('lpstop_ini3D') )then
                allocate(lpinfo(nstages_ini3D))
                call lpstages_fast(params%box, nstages_ini3D, params%smpd, params%lpstart_ini3D, params%lpstop_ini3D, lpinfo)
            else
                THROW_HARD('Both lpstart_ini3D & lpstop_ini3D must be inputted')
            endif
            call cline%delete('lpstart_ini3D')
            call cline%delete('lpstop_ini3D')
        else
            if( cline%defined('lpstart') .and. cline%defined('lpstop') )then
                call set_lplims_from_frcs(spproj, l_cavgs=.true., lpstart=params%lpstart, lpstop=params%lpstop)
            else if( cline%defined('lpstart') )then
                call set_lplims_from_frcs(spproj, l_cavgs=.true., lpstart=params%lpstart)
            else if( cline%defined('lpstop') )then
                call set_lplims_from_frcs(spproj, l_cavgs=.true., lpstop=params%lpstop)
            else
                call set_lplims_from_frcs(spproj, l_cavgs=.true.)
            endif
        endif
        ! whether to use classes generated from 2D or 3D
        select case(trim(params%imgkind))
            case('cavg')
                states  = nint(spproj%os_cls2D%get_all('state'))
            case('cavg3D')
                states  = nint(spproj%os_cls3D%get_all('state'))
            case DEFAULT
                THROW_HARD('Unsupported IMGKIND!')
        end select
        ! retrieve cavgs stack info
        call spproj%get_cavgs_stk(stk, ncavgs, params%smpd, imgkind=params%imgkind)
        if(.not. file_exists(stk)) THROW_HARD('cavgs stk does not exist; simple_commanders_abinitio')
        states          = nint(spproj%os_cls2D%get_all('state'))
        orig_stk        = stk
        ext             = string('.')//fname2ext(stk)
        stk_even        = add2fbody(stk, ext, '_even')
        stk_odd         = add2fbody(stk, ext, '_odd')
        if( .not. file_exists(stk_even) ) THROW_HARD('Even cavgs stk: '//stk_even%to_char()//' does not exist!')
        if( .not. file_exists(stk_odd)  ) THROW_HARD('Odd cavgs stk: '//stk_odd%to_char()//' does not exist!')
        ctfvars%ctfflag = CTFFLAG_NO
        ctfvars%smpd    = params%smpd
        shifted_stk     = basename(add2fbody(stk, ext, '_shifted'))
        if( count(states==0) .eq. ncavgs )then
            THROW_HARD('no class averages detected in project file: '//params%projfile%to_char()//'; abinitio3D_cavgs')
        endif
        ! prepare a temporary project file
        work_proj%projinfo = spproj%projinfo
        work_proj%compenv  = spproj%compenv
        if( spproj%jobproc%get_noris() > 0 ) work_proj%jobproc = spproj%jobproc
        ! name change
        call work_proj%projinfo%delete_entry('projname')
        call work_proj%projinfo%delete_entry('projfile')
        call cline%set('projfile', work_projfile)
        call cline%set('projname', get_fbody(work_projfile,'simple'))
        call work_proj%update_projinfo(cline)
        ! add stks to temporary project
        call work_proj%add_stk(stk_even, ctfvars)
        call work_proj%add_stk(stk_odd,  ctfvars)
        ! update orientations parameters
        do icls=1,ncavgs
            even_ind = icls
            odd_ind  = ncavgs + icls
            call work_proj%os_ptcl3D%get_ori(icls, o)
            call o%set('class', icls)
            call o%set('state', states(icls))
            ! even
            o_even = o
            call o_even%set('eo', 0)
            call o_even%set('stkind', work_proj%os_ptcl3D%get(even_ind,'stkind'))
            call work_proj%os_ptcl3D%set_ori(even_ind, o_even)
            ! odd
            o_odd = o
            call o_odd%set('eo', 1)
            call o_odd%set('stkind', work_proj%os_ptcl3D%get(odd_ind,'stkind'))
            call work_proj%os_ptcl3D%set_ori(odd_ind, o_odd)
        enddo
        params_glob%nptcls = work_proj%get_nptcls()
        call work_proj%write()
        ! Frequency marching
        call set_cline_refine3D(1, l_cavgs=.true.)
        call rndstart(cline_refine3D)
        do istage = 1, nstages_ini3D
            write(logfhandle,'(A)')'>>>'
            write(logfhandle,'(A,I3,A9,F5.1)')'>>> STAGE ', istage,' WITH LP =', lpinfo(istage)%lp
            ! Preparation of command line for probabilistic search
            call set_cline_refine3D(istage, l_cavgs=.true.)
            if( lpinfo(istage)%l_autoscale )then
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',lpinfo(istage)%box_crop
            endif
            ! Probabilistic search
            call exec_refine3D(istage, xrefine3D)
            ! Symmetrization
            if( istage == SYMSRCH_STAGE )then
                call symmetrize(istage, work_proj, work_projfile)
            endif
        end do
        ! update original cls3D segment
        call work_proj%read_segment('ptcl3D', work_projfile)
        call work_proj%read_segment('out',    work_projfile)
        call work_proj%os_ptcl3D%delete_entry('stkind')
        call work_proj%os_ptcl3D%delete_entry('eo')
        params_glob%nptcls = ncavgs
        call spproj%os_cls3D%new(ncavgs, is_ptcl=.false.)
        do icls=1,ncavgs
            if( work_proj%os_ptcl3D%get_state(icls) == 0 )then
                call spproj%os_cls3D%set_state(icls, 0)
            else
                ! e/o orientation with best score is selected
                if( work_proj%os_ptcl3D%get(icls, 'corr') > work_proj%os_ptcl3D%get(ncavgs+icls, 'corr') )then
                    cnt = icls
                else
                    cnt = ncavgs+icls
                endif
                ! alignment parameters
                call spproj%os_cls3D%set(icls, 'corr', work_proj%os_ptcl3D%get(cnt, 'corr'))
                call spproj%os_cls3D%set(icls, 'proj', work_proj%os_ptcl3D%get(cnt, 'proj'))
                call spproj%os_cls3D%set(icls, 'w',    work_proj%os_ptcl3D%get(cnt, 'w'))
                call spproj%os_cls3D%set_euler(icls, work_proj%os_ptcl3D%get_euler(cnt))
                call spproj%os_cls3D%set_shift(icls, work_proj%os_ptcl3D%get_2Dshift(cnt))
                call spproj%os_cls3D%set_state(icls, work_proj%os_ptcl3D%get_state(cnt))
            endif
        enddo
        call spproj%os_cls3D%set_all2single('stkind', 1)    ! revert splitting
        ! map the orientation parameters obtained for the clusters back to the particles
        call spproj%map2ptcls
        if( nstages_ini3D == NSTAGES_INI3D_MAX )then ! produce validation info
            ! check even odd convergence
            if( params%nstates > 1 ) call conv_eo_states(work_proj%os_ptcl3D)
            call conv_eo(work_proj%os_ptcl3D)
            ! calculate 3D reconstruction at original sampling
            call calc_final_rec(work_proj, work_projfile, xreconstruct3D)
            ! postprocess final 3D reconstruction
            call postprocess_final_rec(work_proj)
            ! add rec_final to os_out
            do s = 1,params%nstates
                if( .not.work_proj%isthere_in_osout('vol', s) )cycle
                final_vol = REC_FBODY//int2str_pad(s,2)//params%ext%to_char()
                if( file_exists(final_vol) )then
                    call spproj%add_vol2os_out(final_vol, params%smpd, s, 'vol_cavg')
                endif
            enddo
            ! reprojections
            call spproj%os_cls3D%write(string('final_oris.txt'))
            write(logfhandle,'(A)') '>>>'
            write(logfhandle,'(A)') '>>> RE-PROJECTION OF THE FINAL VOLUME'
            write(logfhandle,'(A)') '>>>'
            do s = 1,params%nstates
                if( .not.work_proj%isthere_in_osout('vol', s) )cycle
                call cline_reproject%set('vol'//int2str(s), REC_FBODY//int2str_pad(s,2)//PPROC_SUFFIX//params_glob%ext%to_char())
            enddo
            call xreproject%execute_safe(cline_reproject)
            ! write alternated stack
            call img%new([params%box,params%box,1],         params%smpd)
            call stkio_r%open(orig_stk,                     params%smpd, 'read',                                 bufsz=500)
            call stkio_r2%open(string('reprojs.mrc'),       params%smpd, 'read',                                 bufsz=500)
            call stkio_w%open(string('cavgs_reprojs.mrc'),  params%smpd, 'write', box=params%box, is_ft=.false., bufsz=500)
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
        endif
        ! write results (this needs to be a full write as multiple segments are updated)
        call spproj%write()
        ! rank classes based on agreement to volume (after writing)
        if( nstages_ini3D == NSTAGES_INI3D_MAX )then
            if( trim(params%rank_cavgs).eq.'yes' ) call rank_cavgs
        endif
        ! end gracefully
        call img%kill
        call spproj%kill
        call o%kill
        call o_even%kill
        call o_odd%kill
        call work_proj%kill
        call del_file(work_projfile)
        call simple_rmdir(string(STKPARTSDIR))
        call simple_end('**** SIMPLE_ABINITIO3D_CAVGS NORMAL STOP ****')
        contains

            subroutine rndstart( cline )
                class(cmdline), intent(inout) :: cline
                type(string) :: src, dest
                character(len=:), allocatable :: state
                integer :: s
                call work_proj%os_ptcl3D%rnd_oris
                call work_proj%os_ptcl3D%zero_shifts
                if( params%nstates > 1 )then
                    call gen_labelling(work_proj%os_ptcl3D, params%nstates, 'uniform')
                endif
                call work_proj%write_segment_inside('ptcl3D', work_projfile)
                call cline%set('mkdir', 'no') ! to avoid nested dirs
                call cline%set('objfun', 'cc')
                call xreconstruct3D%execute_safe(cline)
                call cline%set('objfun', trim(params%objfun))
                do s = 1,params%nstates
                    state = int2str_pad(s,2)
                    src   = VOL_FBODY//state//'.mrc'
                    dest  = STARTVOL_FBODY//state//'.mrc'
                    call simple_rename(src, dest)
                    call cline%set('vol'//int2str(s), dest)
                    src   = VOL_FBODY//state//'_even.mrc'
                    dest  = STARTVOL_FBODY//state//'_even_unfil.mrc'
                    call simple_copy_file(src, dest)
                    dest  = STARTVOL_FBODY//state//'_even.mrc'
                    call simple_rename(src, dest)
                    src   = VOL_FBODY//state//'_odd.mrc'
                    dest  = STARTVOL_FBODY//state//'_odd_unfil.mrc'
                    call simple_copy_file(src, dest)
                    dest  = STARTVOL_FBODY//state//'_odd.mrc'
                    call simple_rename(src, dest)
                enddo
            end subroutine rndstart
    
            subroutine conv_eo( os )
                class(oris), intent(in) :: os
                type(sym) :: se
                type(ori) :: o_odd, o_even
                real      :: avg_euldist, euldist
                integer   :: icls, ncls
                call se%new(params%pgrp)
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

            subroutine conv_eo_states( os )
                class(oris), intent(in) :: os
                real      :: score
                integer   :: icls, nsame_state, se, so
                nsame_state = 0
                do icls = 1,os%get_noris()/2
                    se = os%get_state(icls)
                    so = os%get_state(icls+ncavgs)
                    if( se == so ) nsame_state = nsame_state + 1
                enddo
                score = 100.0 * real(nsame_state) / real(ncavgs)
                write(logfhandle,'(A)')'>>>'
                write(logfhandle,'(A,F6.1,A1)')'>>> EVEN/ODD STATES OVERLAP: ', score,'%'
            end subroutine conv_eo_states

            subroutine rank_cavgs
                use simple_commanders_cavgs, only: commander_rank_cavgs
                type(commander_rank_cavgs) :: xrank_cavgs
                type(cmdline)              :: cline_rank_cavgs
                call cline_rank_cavgs%set('prg',      'rank_cavgs')
                call cline_rank_cavgs%set('projfile', params%projfile)
                call cline_rank_cavgs%set('flag',     'corr') ! rank by cavg vs. reproj agreement
                call cline_rank_cavgs%set('oritype',  'cls3D')
                call cline_rank_cavgs%set('stk',      orig_stk)
                call cline_rank_cavgs%set('outstk',   basename(add2fbody(stk, ext, '_sorted')))
                call xrank_cavgs%execute_safe(cline_rank_cavgs)
                call cline_rank_cavgs%kill
            end subroutine rank_cavgs

    end subroutine exec_abinitio3D_cavgs

    !> for crude generation of an initial 3D model from class averages
    subroutine exec_abinitio3D_cavgs_fast( self, cline )
        class(commander_abinitio3D_cavgs_fast), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline
        type(commander_abinitio3D_cavgs) :: xabinitio3D_cavgs
        type(parameters) :: params
        real             :: lpstart, lpstop
        if( .not.cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        ! resolution limits: lpstart in [12;20], lpstop in [6.;8.]
        lpstart = max(min(params%mskdiam/10., 20.), 12.)
        lpstop  = min(max(params%mskdiam/30.,  6.),  8.)
        if( cline%defined('lpstart') ) lpstart = params%lpstart
        if( cline%defined('lpstop')  ) lpstop  = params%lpstop
        if( lpstop > lpstart ) lpstop = lpstart
        call cline%delete('lpstart')
        call cline%delete('lpstop')
        ! command-line updates
        if( cline%defined('nstates') )then
            if( cline%get_iarg('nstates') > 1 )then
                call cline%set('multivol_mode', 'independent')
            endif
        endif
        call cline%set('mkdir',         'no')
        call cline%set('lp_auto',       'no')
        call cline%set('lpstart_ini3D', lpstart)
        call cline%set('lpstop_ini3D',  lpstop)
        call cline%set('lp_auto',       'no')
        call cline%set('nspace_max',    1500)
        call cline%set('nstages',       NSTAGES_INI3D_MAX)
        call cline%set('rank_cavgs',    'yes')
        ! prune junk
        call prune_junk_classes
        call cline%delete('prune')
        ! execution
        call xabinitio3D_cavgs%execute_safe( cline )
        ! end
        call simple_end('**** SIMPLE_ABINITIO3D_CAVGS_FAST NORMAL STOP ****')
      contains

        subroutine prune_junk_classes
            use simple_strategy2D_utils, only: flag_non_junk_cavgs
            use simple_imgarr_utils,     only: read_cavgs_into_imgarr
            type(sp_project)              :: spproj
            type(image),      allocatable :: cavg_imgs(:)
            logical,          allocatable :: l_non_junk(:)
            integer,          allocatable :: states(:)
            type(string) :: fname
            integer :: i, ncls, j
            if( trim(params%prune).eq.'yes' )then
                call spproj%read(params%projfile)
                cavg_imgs = read_cavgs_into_imgarr(spproj)
                call flag_non_junk_cavgs(cavg_imgs, 20.0, params%msk, l_non_junk, spproj%os_cls2D)
                if( .not.all(l_non_junk) )then
                    ncls = size(cavg_imgs)
                    allocate(states(ncls),source=1)
                    fname = 'cavgs_junk.mrc'
                    j = 0
                    do i = 1, ncls
                        if( .not. l_non_junk(i) )then
                            j = j + 1
                            call cavg_imgs(i)%write(fname, j)
                            call spproj%os_cls2D%set_state(i, 0)
                            states(i) = 0
                        endif
                        call cavg_imgs(i)%kill
                    enddo
                    deallocate(cavg_imgs)
                    call spproj%map_cavgs_selection(states)
                    call spproj%write(params%projfile)
                    write(logfhandle,'(A,I5)') '>>> # classes left after junk rejection ', count(l_non_junk)
                    call fname%kill
                endif
                call spproj%kill
            endif
        end subroutine prune_junk_classes


    end subroutine exec_abinitio3D_cavgs_fast

    !> for generation of an initial 3d model from particles
    subroutine exec_abinitio3D( self, cline )
        class(commander_abinitio3D), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        ! commanders
        type(commander_refine3D_distr)         :: xrefine3D
        type(commander_reconstruct3D_distr)    :: xreconstruct3D_distr
        ! other
        real,               allocatable :: rstates(:)
        integer,            allocatable :: tmpinds(:), clsinds(:), pinds(:)
        type(class_sample), allocatable :: clssmp(:)
        type(parameters)                :: params
        type(sp_project)                :: spproj
        type(simple_nice_communicator)  :: nice_communicator
        integer :: istage, icls, start_stage, nptcls2update, noris, nstates_on_cline, nstates_in_project, split_stage
        logical :: l_stream
        call cline%set('objfun',    'euclid') ! use noise normalized Euclidean distances from the start
        call cline%set('sigma_est', 'global') ! obviously
        call cline%set('bfac',            0.) ! because initial models should not be sharpened
        if( .not. cline%defined('mkdir')               ) call cline%set('mkdir',                              'yes')
        if( .not. cline%defined('overlap')             ) call cline%set('overlap',                             0.95)
        if( .not. cline%defined('prob_athres')         ) call cline%set('prob_athres',                          10.)
        if( .not. cline%defined('center')              ) call cline%set('center',                              'no')
        if( .not. cline%defined('cenlp')               ) call cline%set('cenlp',                      CENLP_DEFAULT)
        if( .not. cline%defined('oritype')             ) call cline%set('oritype',                         'ptcl3D')
        if( .not. cline%defined('pgrp')                ) call cline%set('pgrp',                                'c1')
        if( .not. cline%defined('pgrp_start')          ) call cline%set('pgrp_start',                          'c1')
        if( .not. cline%defined('ptclw')               ) call cline%set('ptclw',                               'no')
        if( .not. cline%defined('projrec')             ) call cline%set('projrec',                            'yes')
        if( .not. cline%defined('lp_auto')             ) call cline%set('lp_auto',                            'yes')
        if( .not. cline%defined('first_sigmas')        ) call cline%set('first_sigmas',                        'no')
        if( .not. cline%defined('ref_type')            ) call cline%set('ref_type',                 'comlin_noself')
        if( .not. cline%defined('gauref_last_stage')   ) call cline%set('gauref_last_stage',      GAUREF_LAST_STAGE)
        if( .not. cline%defined('inivol')              ) call cline%set('inivol',                          'sphere')
        if( .not. cline%defined('maxits_between')      ) call cline%set('maxits_between',            MAXITS_BETWEEN)
        ! adjust cartesian/polar options
        if( cline%get_carg('polar')=='yes' )then
        if( .not. cline%defined('gauref')              ) call cline%set('gauref',                             'yes')
        else
        if( .not. cline%defined('gauref')              ) call cline%set('gauref',                              'no')
        endif
        ! splitting stage
        split_stage = HET_DOCKED_STAGE
        if( cline%defined('split_stage') ) split_stage = cline%get_iarg('split_stage')
        ! adjust default multivol_mode unless given on command line
        if( cline%defined('nstates') )then
            nstates_on_cline = cline%get_iarg('nstates')
            if( nstates_on_cline > 1 .and. .not. cline%defined('multivol_mode') )then
                call cline%set('multivol_mode', 'independent')
            endif
        endif
        ! make master parameters
        l_stream = .false.
        if( cline%defined('stream') ) l_stream = cline%get_carg('stream').eq.'yes'
        call cline%delete('stream')
        call params%new(cline)
        call cline%set('mkdir', 'no')
        call cline%delete('algorithm')
        call cline%delete('maxits_between')
        ! Multiple states
        nstates_glob = params%nstates
        select case(trim(params%multivol_mode))
            case('single')
                if( nstates_glob /= 1 ) THROW_HARD('nstates /= 1 incompatible with multivol_mode:' //trim(params%multivol_mode))
            case('independent', 'docked', 'input_oris_start', 'input_oris_fixed')
                if( nstates_glob == 1 ) THROW_HARD('nstates == 1 incompatible with multivol_mode: '//trim(params%multivol_mode))
            case DEFAULT
                THROW_HARD('Unsupported multivol_mode: '//trim(params%multivol_mode))
        end select
        if( trim(params%multivol_mode).eq.'docked' )then
            params%nstates = 1
            call cline%delete('nstates')
        endif
        ! Polar representation
        if( trim(params%polar).eq.'yes' )then
            if( trim(params%multivol_mode).ne.'single' )then
                THROW_HARD('POLAR=YES not compatible with MULTIVOL_MODE='//trim(params%multivol_mode))
            endif
            if( trim(params%lp_auto).eq.'yes' )then
                THROW_WARN('POLAR=YES not compatible LP_AUTO=YES; reverting to LP_AUTO=NO')
            endif
            params%lp_auto = 'no'; params%l_lpauto = .false.; l_lpauto=.false.
            call cline%set('lp_auto', 'no')
            l_polar = .true. ! global parameter
        else
            call cline%delete('ref_type')
        endif
        ! nice communicator init
        call nice_communicator%init(params%niceprocid, params%niceserver)
        call nice_communicator%cycle()
        ! read project
        call spproj%read(params%projfile)
        ! provide initialization of 3D alignment using class averages?
        start_stage = 1
        l_ini3D     = .false.
        if( trim(params%cavg_ini).eq.'yes' )then
            if( str_has_substr(params%multivol_mode,'input_oris') ) THROW_HARD('Ini3D on cavgs not allowed for multivol_mode=input_oris*')
            ! nice
            nice_communicator%stat_root%stage = "initialising 3D volume from class averages"
            call nice_communicator%cycle()
            call ini3D_from_cavgs(cline)
            ! re-read the project file to update info in spproj
            call spproj%read(params%projfile)
            start_stage = NSTAGES_INI3D - 1 ! compute reduced to two overlapping stages
            l_ini3D     = .true.
            ! symmetry dealt with by ini3D
        else
            if( trim(params%multivol_mode).eq.'independent' )then
                ! turn off symmetry axis search and put the symmetry in from the start
                params%pgrp_start = params%pgrp
            endif
        endif
        ! nice
        nice_communicator%stat_root%stage = "preparing workflow"
        call nice_communicator%cycle()
        ! initialization on class averages done outside this workflow (externally)?
        if( trim(params%cavg_ini_ext).eq.'yes' )then 
            ! check that ptcl3D field is not virgin
            if( spproj%is_virgin_field('ptcl3D') )then
                THROW_HARD('Prior 3D alignment required for abinitio workflow when cavg_ini_ext is set to yes')
            endif
            start_stage = NSTAGES_INI3D - 1 ! compute reduced to two overlapping stages
            l_ini3D     = .true.
        endif
        ! set class global lp_auto flag for low-pass limit estimation
        l_lpauto = .true.
        if( cline%defined('lp_auto') ) l_lpauto = params%l_lpauto
        ! set class global automasking flag
        l_automsk = .false.
        if( cline%defined('automsk') )then
            if( trim(params%automsk).eq.'yes' )then
                if( trim(params%multivol_mode).eq.'single' )then
                    l_automsk = .true.
                else
                    THROW_WARN('automasking not supported for modes other than multivol_mod.eq.single, turning automasking off')
                    l_automsk = .false.
                endif
            endif
        endif
        ! prepare class command lines
        call prep_class_command_lines(cline, params%projfile)
        ! set symmetry class variables
        call set_symmetry_class_vars
        ! fall over if there are no particles
        if( spproj%os_ptcl3D%get_noris() < 1 ) THROW_HARD('Particles could not be found in the project')
        ! take care of class-biased particle sampling
        if( spproj%is_virgin_field('ptcl2D') )then
            THROW_HARD('Prior 2D clustering required for abinitio workflow')
        else
            l_update_frac_dyn    = .false.
            l_nsample_given      = .false.
            l_nsample_stop_given = .false.
            update_frac          = 1.0
            nptcls_eff           = spproj%count_state_gt_zero()
            if( cline%defined('nsample') )then
                update_frac = real(params%nsample * params%nstates) / real(nptcls_eff)
                l_nsample_given = .true.
            else if( cline%defined('update_frac') )then
                update_frac = params%update_frac
                l_nsample_given = .true.
            else if( cline%defined('nsample_start') )then
                if( params%nsample_start > nptcls_eff ) THROW_HARD('nsample_start > effective # ptcls, decrease!')
                nsample_minmax(1) = params%nsample_start
                if( cline%defined('nsample_stop') )then
                    nsample_minmax(2)    = min(nptcls_eff,params%nsample_stop)
                    l_nsample_stop_given = .true.
                else
                    nsample_minmax(2) = nptcls_eff
                endif
                update_frac       = calc_update_frac_dyn(nptcls_eff, params%nstates, nsample_minmax, 1, MAXITS_GLOB)
                l_update_frac_dyn = .true.
            else
                if( cline%defined('nsample_max') )then
                    update_frac = calc_update_frac(nptcls_eff, params%nstates, [NSAMPLE_MINMAX_DEFAULT(1),params%nsample_max])
                else
                    update_frac = calc_update_frac(nptcls_eff, params%nstates, NSAMPLE_MINMAX_DEFAULT)
                endif
            endif
            update_frac = min(UPDATE_FRAC_MAX, update_frac) ! to ensure fractional update is always on
            ! generate a data structure for class sampling on disk
            rstates = spproj%os_cls2D%get_all('state')
            if( trim(params%partition).eq.'yes' )then
                tmpinds = nint(spproj%os_cls2D%get_all('cluster'))
                where( rstates < 0.5 ) tmpinds = 0
                clsinds = (/(icls,icls=1,maxval(tmpinds))/)
                do icls = 1,size(clsinds)
                    if(count(tmpinds==icls) == 0) clsinds(icls) = 0
                enddo
                clsinds = pack(clsinds, mask=clsinds>0)
                call spproj%os_ptcl2D%get_class_sample_stats(clsinds, clssmp, label='cluster')
                deallocate(tmpinds)
            else
                clsinds = spproj%get_selected_clsinds()
                call spproj%os_ptcl2D%get_class_sample_stats(clsinds, clssmp)
            endif
            call write_class_samples(clssmp, string(CLASS_SAMPLING_FILE))
            deallocate(rstates, clsinds)
            if( spproj%os_ptcl3D%has_been_sampled() )then
                ! the ptcl3D field should be clean of sampling at this stage
                call spproj%os_ptcl3D%clean_entry('sampled')
                ! call spproj%os_ptcl3D%clean_entry('sampled', 'updatecnt')
                call spproj%write_segment_inside('ptcl3D', params%projfile)
            endif
        endif
        ! set low-pass limits and downscaling info from FRCs
         if( cline%defined('lpstart') .and. cline%defined('lpstop') )then
            call set_lplims_from_frcs(spproj, l_cavgs=.false., lpstart=params%lpstart, lpstop=params%lpstop)
        else if( cline%defined('lpstart') )then
            call set_lplims_from_frcs(spproj, l_cavgs=.false., lpstart=params%lpstart)
        else if( cline%defined('lpstop') )then
            call set_lplims_from_frcs(spproj, l_cavgs=.false., lpstop=params%lpstop)
        else
            call set_lplims_from_frcs(spproj, l_cavgs=.false.)
        endif
        ! starting volume logics
        if( str_has_substr(params%multivol_mode,'input_oris') )then
            ! check that ptcl3D field is not virgin
            if( spproj%is_virgin_field('ptcl3D') )then
                THROW_HARD('Prior 3D alignment is lacking for multi-volume assignment')
            endif
            ! create an initial sampling of all updated ptcls for 3D reconstruction
            noris = spproj%os_ptcl3D%get_noris()
            call spproj%os_ptcl3D%sample4update_updated([1,noris], nptcls2update, pinds, .true.)
            call spproj%os_ptcl3D%set_updatecnt(1, pinds) ! set all sampled updatecnts to 1 & the rest to zero
            deallocate(pinds) ! these are not needed
            ! start at the same stage as for multivol_mode==docked
            start_stage = split_stage
            ! create state labelling
            nstates_in_project = spproj%os_ptcl3D%get_n('state')
            if( nstates_in_project == params%nstates )then
                THROW_WARN('exec_abinitio3D: prior nstates equal to given nstates. No randomization!')
            elseif( nstates_in_project == 1 )then
                THROW_WARN('No previous state assignment detected in project. Randomizing states!')
                call gen_labelling(spproj%os_ptcl3D, params%nstates, 'squared_uniform')
            else
                THROW_HARD('Previous state assignment inconsistent with given number of states!')
            endif
            ! write updated project file
            call spproj%write_segment_inside(params%oritype, params%projfile)
            ! calc recs
            call calc_start_rec(params%projfile, xreconstruct3D_distr, start_stage)
        else if( .not. l_ini3D )then
            ! the ptcl3D field should be clean of updates at this stage
            call spproj%os_ptcl3D%clean_entry('updatecnt')
            ! randomize projection directions
            select case(trim(params%oritype))
                case('ptcl3D')
                    call spproj%os_ptcl3D%rnd_oris
                case DEFAULT
                    THROW_HARD('Unsupported ORITYPE; exec_abinitio3D')
            end select
            ! randomize states
            if( trim(params%multivol_mode).eq.'independent' )then
                call gen_labelling(spproj%os_ptcl3D, params%nstates, 'squared_uniform')
            endif
            call spproj%write_segment_inside(params%oritype, params%projfile)
            ! create noise starting volume(s)
            call generate_random_volumes( lpinfo(1)%box_crop, lpinfo(1)%smpd_crop, cline_refine3D )
        else
            ! check that ptcl3D field is not virgin
            if( spproj%is_virgin_field('ptcl3D') )then
                THROW_HARD('Prior 3D alignment is lacking for starting volume generation')
            endif
            ! randomize states
            if( trim(params%multivol_mode).eq.'independent' )then
                call gen_labelling(spproj%os_ptcl3D, params%nstates, 'squared_uniform')
            endif
            ! create an initial balanced greedy sampling
            noris = spproj%os_ptcl3D%get_noris()
            call spproj%os_ptcl3D%sample4update_class(clssmp, [1,noris], update_frac, nptcls2update, pinds, .true., .true.)
            call spproj%os_ptcl3D%set_updatecnt(1, pinds) ! set all sampled updatecnts to 1 & the rest to zero
            deallocate(pinds)                             ! these are not needed
            call deallocate_class_samples(clssmp)         ! done with this one
            ! write updated project file
            call spproj%write_segment_inside(params%oritype, params%projfile)
            ! create starting volume(s)
            call calc_start_rec(params%projfile, xreconstruct3D_distr, start_stage)
        endif
        ! Frequency marching
        maxits_dyn = 0
        if( start_stage < NSTAGES )then
            maxits_dyn = sum(MAXITS(start_stage:NSTAGES - 1)) ! the last stage is omitted in this estimate since the sampling method changes
        endif
        ! nice
        nice_communicator%stat_root%stage = "starting workflow"
        call nice_communicator%cycle()
        do istage = start_stage, NSTAGES
            ! nice
             if( nice_communicator%stop )then
                ! termination
                write(logfhandle,'(A)')'>>> USER COMMANDED STOP'
                call spproj%kill
                call qsys_cleanup
                call nice_communicator%terminate(stop=.true.)
                call simple_end('**** SIMPLE_ABINITIO3D USER STOP ****')
                call EXIT(0)
            endif
            nice_communicator%stat_root%stage = "running workflow"
            call nice_communicator%update_ini3D(stage=istage, number_states=nstates_glob, lp=lpinfo(istage)%lp) 
            call nice_communicator%cycle()
            write(logfhandle,'(A)')'>>>'
            write(logfhandle,'(A,I3,A9,F5.1)')'>>> STAGE ', istage,' WITH LP =', lpinfo(istage)%lp
            ! At the splitting stage of docked mode: reset the nstates in params
            if( params%multivol_mode.eq.'docked' .and. istage == split_stage )then
                params_glob%nstates = nstates_glob
                update_frac         = min(update_frac * nstates_glob, UPDATE_FRAC_MAX)
            endif
            ! Preparation of command line for refinement
            call set_cline_refine3D(istage, l_cavgs=.false.)
            ! Need to be here since rec cline depends on refine3D cline
            if( params%multivol_mode.eq.'docked' .and. istage == split_stage )then
                call randomize_states(spproj, params%projfile, xreconstruct3D_distr, istage=split_stage)
            else if( istage >= RECALC_STARTREC_STAGE )then
                if( .not.l_polar ) call calc_start_rec(params%projfile, xreconstruct3D_distr, istage)
            endif
            if( lpinfo(istage)%l_autoscale )then
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',lpinfo(istage)%box_crop
            endif
            ! Reconstruction for polar representation
            if( l_polar ) call calc_rec4polar( xreconstruct3D_distr, istage )
            ! Executing the refinement with the above settings
            call exec_refine3D(istage, xrefine3D)
            ! Symmetrization
            if( istage == SYMSRCH_STAGE )then
                call symmetrize(istage, spproj, params%projfile, xreconstruct3D_distr)
            endif
            ! nice
            call nice_communicator%update_ini3D(last_stage_completed=.true.) 
            call nice_communicator%cycle()
        enddo
        ! calculate 3D reconstruction at original sampling
        call calc_final_rec(spproj, params%projfile, xreconstruct3D_distr)
        ! for visualization
        call gen_ortho_reprojs4viz(spproj)
        ! postprocess final 3D reconstruction
        call postprocess_final_rec(spproj)
        ! termination
        nice_communicator%stat_root%stage = "terminating"
        call nice_communicator%cycle()
        ! cleanup
        call nice_communicator%terminate(export_project=spproj)
        call spproj%kill
        call qsys_cleanup
        if( l_stream ) call simple_touch(ABINITIO3D_FINISHED)
        call simple_end('**** SIMPLE_ABINITIO3D NORMAL STOP ****')
    end subroutine exec_abinitio3D

    subroutine exec_multivol_assign( self, cline )
        class(commander_multivol_assign), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        type(commander_abinitio3D) :: xabini3D
        type(string) :: srch_oris
        call cline%set('center',   'no')
        call cline%set('cavg_ini', 'no')
        call cline%set('prg',      'multivol_assign')
        if( .not. cline%defined('nstates')  ) THROW_HARD('nstates required on command line')
        srch_oris = 'yes'
        if( cline%defined('srch_oris') )then
            srch_oris = cline%get_carg('srch_oris')
        endif
        select case(srch_oris%to_char())
            case('yes')
                call cline%set('multivol_mode', 'input_oris_start')
            case('no')
                call cline%set('multivol_mode', 'input_oris_fixed')
            case DEFAULT
                THROW_HARD('Unsupported srch_oris flag')
        end select
        call xabini3D%execute(cline)
    end subroutine exec_multivol_assign

    ! private helper routines

    subroutine prep_class_command_lines( cline, projfile )
        class(cmdline), intent(in) :: cline
        class(string),  intent(in) :: projfile
        cline_refine3D      = cline
        cline_symmap        = cline
        cline_reconstruct3D = cline
        cline_postprocess   = cline
        cline_reproject     = cline
        ! refine3D
        call cline_refine3D%set('prg',                         'refine3D')
        call cline_refine3D%set('pgrp',                  params_glob%pgrp)
        call cline_refine3D%set('projfile',                      projfile)
        ! symmetrization
        call cline_symmap%set('prg',                     'symmetrize_map')
        call cline_symmap%set('pgrp',                    params_glob%pgrp)
        call cline_symmap%set('projfile',                        projfile)
        call cline_symmap%set('center',                             'yes')
        if( .not. cline_symmap%defined('cenlp') )then
        call cline_symmap%set('cenlp',                      CENLP_DEFAULT)
        endif
        call cline_symmap%set('hp',                        params_glob%hp)
        ! re-reconstruct volume
        call cline_reconstruct3D%set('prg',               'reconstruct3D')
        call cline_reconstruct3D%set('box',               params_glob%box)
        call cline_reconstruct3D%set('smpd',             params_glob%smpd)
        call cline_reconstruct3D%set('projfile',                 projfile)
        call cline_reconstruct3D%set('pgrp',             params_glob%pgrp)
        call cline_reconstruct3D%set('ml_reg',                       'no')
        call cline_reconstruct3D%set('needs_sigma',                  'no')
        call cline_reconstruct3D%set('objfun',                       'cc')
        call cline_reconstruct3D%delete('polar')
        ! no fractional update
        call cline_reconstruct3D%delete('update_frac')
        ! individual particles reconstruction
        call cline_reconstruct3D%set('projrec', 'no')
        ! postprocess volume
        call cline_postprocess%set('prg',                   'postprocess')
        call cline_postprocess%set('projfile',                   projfile)
        call cline_postprocess%set('mkdir',                          'no')
        call cline_postprocess%set('imgkind',                       'vol')
        call cline_postprocess%delete('lp')   ! to obtain optimal filtration
        ! re-project volume, only with cavgs
        call cline_reproject%set('prg',                       'reproject')
        call cline_reproject%set('pgrp',                 params_glob%pgrp)
        call cline_reproject%set('outstk', 'reprojs'//params_glob%ext%to_char())
        call cline_reproject%set('smpd',                 params_glob%smpd)
        call cline_reproject%set('box',                   params_glob%box)
        call cline_reproject%set('oritab',               'final_oris.txt')
        call cline_reproject%set('nstates',           params_glob%nstates)
        call cline_reproject%delete('projfile')
    end subroutine prep_class_command_lines

    subroutine set_symmetry_class_vars
        type(string) :: pgrp, pgrp_start
        pgrp           = lowercase(trim(params_glob%pgrp))
        pgrp_start     = lowercase(trim(params_glob%pgrp_start))
        l_srch4symaxis = pgrp .ne. pgrp_start
        l_symran       = .false.
        l_sym          = l_srch4symaxis
        if( pgrp_start.ne.'c1' .or. pgrp.ne.'c1' )then
            se1 = sym(pgrp_start%to_char())
            se2 = sym(pgrp%to_char())
            if(se1%get_nsym() > se2%get_nsym())then
                ! ensure se2 is a subgroup of se1
                if( .not. se1%has_subgrp(pgrp%to_char()) )THROW_HARD('Incompatible symmetry groups; exec_abinitio3D')
                ! set flag for symmetry randomisation
                ! in case we are moving from a higher to lower group
                l_symran = .true.
            else if( se2%get_nsym() > se1%get_nsym() )then
                ! ensure se1 is a subgroup of se2
                if( .not. se2%has_subgrp(pgrp_start%to_char()) )THROW_HARD('Incompatible symmetry groups; exec_abinitio3D')
            endif
        endif
    end subroutine set_symmetry_class_vars

    subroutine set_lplims_from_frcs( spproj, l_cavgs, lpstart, lpstop )
        class(sp_project), intent(inout) :: spproj
        logical,           intent(in)    :: l_cavgs
        real, optional,    intent(in)    :: lpstart, lpstop
        real,         allocatable :: frcs_avg(:)
        integer,      allocatable :: states(:)
        type(string)     :: frcs_fname
        type(class_frcs) :: clsfrcs
        real             :: lpfinal
        integer          :: filtsz
        ! retrieve FRC info
        call spproj%get_frcs(frcs_fname, 'frc2D', fail=.false.)
        ! work out low-pass limits and downscaling parameters
        params_glob%frcs = frcs_fname
        call clsfrcs%read(frcs_fname)
        filtsz = clsfrcs%get_filtsz()
        allocate(frcs_avg(filtsz), source=0.)
        states = nint(spproj%os_cls2D%get_all('state'))
        call clsfrcs%avg_frc_getter(frcs_avg, states)
        if( allocated(lpinfo) ) deallocate(lpinfo)
        allocate(lpinfo(NSTAGES))
        lpfinal = max(LPSTOP_BOUNDS(1),calc_lplim_final_stage(3))
        lpfinal = min(LPSTOP_BOUNDS(2),lpfinal)
        if( present(lpstop) ) lpfinal = max(lpstop,lpfinal)
        if( present(lpstart) )then
            call lpstages(params_glob%box, NSTAGES, frcs_avg, params_glob%smpd,&
            &lpstart, lpstart, lpfinal, lpinfo, l_cavgs, verbose=.not.l_polar)
        else
            call lpstages(params_glob%box, NSTAGES, frcs_avg, params_glob%smpd,&
            &LPSTART_BOUNDS(1), LPSTART_BOUNDS(2), lpfinal, lpinfo, l_cavgs, verbose=.not.l_polar)
        endif
        ! Stepped increase of dimensions with polar representation
        call edit_lpstages4polar( NSPACE_PHASE_POLAR, NSTAGES, lpinfo )
        ! cleanup
        call clsfrcs%kill
        contains

            function calc_lplim_final_stage( nbest ) result( lplim )
                integer, intent(in)  :: nbest
                real,    allocatable :: res(:), tmp_rarr(:)
                integer, allocatable :: tmp_iarr(:)
                real :: lplim
                tmp_rarr  = spproj%os_cls2D%get_all('res')
                tmp_iarr  = nint(spproj%os_cls2D%get_all('state'))
                res       = pack(tmp_rarr, mask=(tmp_iarr>0))
                call hpsort(res)
                lplim = median_nocopy(res(:nbest))
                deallocate(tmp_rarr, tmp_iarr, res)
            end function calc_lplim_final_stage

    end subroutine set_lplims_from_frcs

    subroutine ini3D_from_cavgs( cline )
        class(cmdline),    intent(inout) :: cline
        type(commander_abinitio3D_cavgs) :: xini3D
        type(cmdline)                    :: cline_ini3D
        type(string),    allocatable     :: files_that_stay(:)
        character(len=*), parameter      :: INI3D_DIR='abinitio3D_cavgs/'
        real,             parameter      :: LPSTART_INI3D = 20.
        real,             parameter      :: LPSTOP_INI3D  = 6.
        cline_ini3D = cline
        call cline_ini3D%set('nstages', NSTAGES_INI3D)
        if( .not. cline_ini3D%defined('lpstart_ini3D') ) call cline_ini3D%set('lpstart_ini3D', LPSTART_INI3D)
        if( .not. cline_ini3D%defined('lpstop_ini3D')  ) call cline_ini3D%set('lpstop_ini3D',  LPSTOP_INI3D)
        if( cline%defined('lpstart_ini3D') )then
            call cline_ini3D%set('lpstart', params_glob%lpstart_ini3D)
            call cline_ini3D%delete('lpstart_ini3D')
        endif
        if( cline%defined('lpstop_ini3D') )then
            call cline_ini3D%set('lpstop', params_glob%lpstop_ini3D)
            call cline_ini3D%delete('lpstop_ini3D')
        endif
        if( cline%defined('nthr_ini3D') )then
            call cline_ini3D%set('nthr', params_glob%nthr_ini3D)
            call cline_ini3D%delete('nthr_ini3D')
        endif
        call cline_ini3D%delete('nstates') ! cavg_ini under the assumption of one state
        call cline_ini3D%delete('nparts')
        call cline_ini3D%delete('projrec')
        call cline_ini3D%delete('oritype')
        call cline_ini3D%delete('imgkind')
        call cline_ini3D%delete('prob_athres')
        call xini3D%execute_safe(cline_ini3D)
        ! update point-group symmetry
        call cline%set('pgrp_start', params_glob%pgrp)
        params_glob%pgrp_start = params_glob%pgrp
        ! stash away files
        ! identfy files that stay
        allocate(files_that_stay(7))
        files_that_stay(1) = basename(params_glob%projfile)
        files_that_stay(2) = 'cavgs'
        files_that_stay(3) = 'nice'
        files_that_stay(4) = 'frcs'
        files_that_stay(5) = 'ABINITIO3D'
        files_that_stay(6) = 'execscript' ! only with streaming
        files_that_stay(7) = 'execlog'    ! only with streaming
        ! make the move
        call move_files_in_cwd(string(INI3D_DIR), files_that_stay)
    end subroutine ini3D_from_cavgs

    subroutine set_cline_refine3D( istage, l_cavgs )
        integer,          intent(in)  :: istage
        logical,          intent(in)  :: l_cavgs
        type(string) :: sh_first, prob_sh, ml_reg, fillin, cavgw, ref_type
        type(string) :: refine, icm, trail_rec, pgrp, balance, lp_auto, automsk
        integer :: iphase, iter, inspace, imaxits, nsample_dyn, nspace_phase
        real    :: trs, frac_best, overlap, fracsrch, lpstart, lpstop, snr_noise_reg, gaufreq
        ! iteration number bookkeeping
        iter = 0
        if( cline_refine3D%defined('endit') )then
            iter = cline_refine3D%get_iarg('endit')
        endif
        iter = iter + 1
        ! dynamic update frac
        if( istage == NSTAGES )then
            fillin = 'yes'
            if( params_glob%nstates > 1 ) fillin = 'no' ! fill-in doesn't work with multi-state
            if( l_nsample_stop_given )then
                update_frac_dyn = real(nsample_minmax(2)) / real(nptcls_eff)
            else if( l_nsample_given )then
                update_frac_dyn = update_frac
            else
                ! we change the sampling method for the last stage (accelerated refinement)
                nsample_dyn     = nint(UPDATE_FRAC_MIN * real(nptcls_eff) / real(params_glob%nstates))
                nsample_dyn     = max(NSAMPLE_MINMAX_DEFAULT(1), min(NSAMPLE_MINMAX_DEFAULT(2), nsample_dyn))
                update_frac_dyn = real(nsample_dyn * params_glob%nstates) / real(nptcls_eff)
            endif
        else
            fillin = 'no'
            update_frac_dyn = calc_update_frac_dyn(nptcls_eff, params_glob%nstates, nsample_minmax, iter, maxits_dyn)
        endif
        update_frac_dyn = min(UPDATE_FRAC_MAX, update_frac_dyn) ! to ensure fractional update is always on
        ! symmetry
        pgrp = trim(params_glob%pgrp)
        if( l_srch4symaxis )then
            if( istage <= SYMSRCH_STAGE )then
                ! need to replace original point-group flag with c1/pgrp_start
                pgrp = trim(params_glob%pgrp_start)
            endif
        endif
        ! refinement mode
        refine  = 'shc_smpl'
        prob_sh = 'no'
        if( istage >= PROBREFINE_STAGE )then
            refine  = 'prob'
            prob_sh = 'yes'
        endif
        if( trim(params_glob%multivol_mode).eq.'input_oris_fixed' )then ! only state sorting, no 3D ori refinement
            refine = 'prob_state'
        endif
        ! ICM regularization
        icm = 'no'
        if( istage >= ICM_STAGE ) icm = 'yes'
        ! forcing icm to no if user sets icm to no
        if( cline_refine3D%defined('icm') )then
            if( trim(params_glob%icm).eq.'no')then
                icm = 'no'
            endif
        endif
        ! balance
        balance = 'yes'
        ! Gaussian filtering of reference volume
        gaufreq = -1.
        if( (trim(params_glob%gauref)=='yes') .and. (istage <= params_glob%gauref_last_stage) )then
            gaufreq = lpinfo(istage)%lp
        endif
        ! trailing reconstruction
        trail_rec = 'no'
        select case(trim(params_glob%multivol_mode))
            case('single')
                if( istage >= TRAILREC_STAGE_SINGLE ) trail_rec = 'yes'
            case('independent')
                if( istage >= TRAILREC_STAGE_MULTI  ) trail_rec = 'yes'
            case('docked')
                if( istage == NSTAGES )then
                    trail_rec = 'no'
                else if( istage >= TRAILREC_STAGE_SINGLE )then
                    trail_rec = 'yes'
                endif
            case('input_oris_fixed')
                trail_rec = 'no'
            case('input_oris_start')
                trail_rec = 'no'
            case DEFAULT
                trail_rec = 'no' ! defaults to 'no' for safety
        end select
        ! automatic low-pass limit estimation
        lp_auto = 'no'
        if( istage >= LPAUTO_STAGE .and. l_lpauto )then
            lp_auto = trim(params_glob%lp_auto)
            lpstart = lpinfo(istage - 1)%lp
            if( istage == NSTAGES )then
                lpstop = lpinfo(istage)%smpd_crop * 2. ! Nyqvist limit
            else
                lpstop = lpinfo(istage + 1)%lp
            endif
        endif
        ! automasking
        automsk = 'no'
        if( (.not.l_cavgs) .and. (.not.l_polar) )then
            if( istage >= AUTOMSK_STAGE .and. l_automsk )then
                automsk = 'yes'
            endif
        endif
        ! cavgs weights, not supported for particles
        cavgw = 'no'
        if( l_cavgs )then
            if( (trim(params_glob%cavgw).eq.'yes') .and. (istage>=CAVGWEIGHTS_STAGE))then
                cavgw = 'yes'
            endif
        endif
        ! phase logics
        iphase = 0
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
                inspace       = NSPACE(1)
                imaxits       = MAXITS(istage)
                trs           = 0.
                sh_first      = 'no'
                ml_reg        = 'no'
                frac_best     = 1.0 ! means it does not control sampling, greedy selection
                overlap       = 0.90
                fracsrch      = 90.
                snr_noise_reg = 2.0
            case(2)
                inspace       = NSPACE(2)
                imaxits       = MAXITS(istage)
                trs           = lpinfo(istage)%trslim
                sh_first      = 'yes'
                ml_reg        = 'yes'
                if( istage >= STOCH_SAMPL_STAGE )then
                frac_best     = 0.5 ! means sampling is done from top-ranking 50% particles in class
                else
                frac_best     = 1.0 ! means it does not control sampling, greedy selection
                endif
                overlap       = 0.90
                fracsrch      = 90.
                snr_noise_reg = 4.0
            case(3)
                inspace       = NSPACE(3)
                imaxits       = MAXITS(istage)
                trs           = lpinfo(istage)%trslim
                sh_first      = 'yes'
                ml_reg        = 'yes'
                if( params_glob%nstates > 1 )then
                ! turn off balancing
                frac_best     = 0.98 ! max out balanced sampling
                else
                frac_best     = 0.85 ! means sampling is done from top-ranking 85% particles in class
                endif
                if( istage == NSTAGES )then
                overlap       = 0.95
                fracsrch      = 95.
                else
                overlap       = 0.90
                fracsrch      = 90.
                endif
                snr_noise_reg = 6.0
        end select
        ! Specific options for polar representation
        if( l_polar )then
            ! # of projection directions
            nspace_phase = 3
            if( istage <= NSPACE_PHASE_POLAR(1) )then
                nspace_phase = 1
            else if( istage <= NSPACE_PHASE_POLAR(2) )then
                nspace_phase = 2
            endif
            inspace  = NSPACE_POLAR(nspace_phase)
            ! volume filtering
            icm = 'no'
            if( trim(params_glob%gauref).eq.'no' ) gaufreq = -1.
            if( ml_reg.eq.'yes' )                  gaufreq = -1.
            ! CL-based approach
            ref_type = trim(params_glob%ref_type)
            if( (istage >= params_glob%switch_reftype_stage) .and.&
            &(params_glob%switch_reftype_stage > 0 )) ref_type = 'comlin'
        endif
        ! turn off ML-regularization when icm is on
        if( icm.eq.'yes' ) ml_reg = 'no'
        ! projection directions
        if( cline_refine3D%defined('nspace_max') )then
            inspace = min(inspace, params_glob%nspace_max)
        endif
        ! Development
        if( trim(params_glob%algorithm).eq.'mimic2D' )then
            imaxits  = params_glob%maxits_between   ! # of iterations
            inspace  = min(inspace, 2000)           ! nspace
            refine   = 'snhc_smpl'                  ! refinement
            sh_first = 'no'                         ! not implemented
            prob_sh  = 'no'
            lp_auto  = 'no'
            icm      = 'no'                         ! icm filter off
            if( istage <= params_glob%gauref_last_stage )then
             gaufreq = lpinfo(istage)%lp            ! gaussian on, ML off
             ml_reg  = 'no'
            else
             gaufreq = -1.0                         ! gaussian off, ML on
             ml_reg  = 'yes'
            endif
            call cline_refine3D%set('extr_lim', (NSTAGES-1)*imaxits)
            if( istage <= 1 )then
                call cline_refine3D%set('extr_iter', 1)
            elseif( istage < NSTAGES )then
                call cline_refine3D%set('extr_iter', (istage-1)*imaxits+1)
            else
                call cline_refine3D%delete('extr_iter')
                call cline_refine3D%delete('extr_lim')
            endif
        endif
        ! command line update
        call cline_refine3D%set('prg',                     'refine3D')
        ! class global control parameters
        if( l_update_frac_dyn .or. istage == NSTAGES )then
        call cline_refine3D%set('update_frac',        update_frac_dyn)
        call cline_refine3D%set('fillin',                      fillin)
        else
        call cline_refine3D%set('update_frac',            update_frac)
        call cline_refine3D%delete('fillin')
        endif
        call cline_refine3D%set('lp',             lpinfo(istage  )%lp)
        call cline_refine3D%set('smpd_crop', lpinfo(istage)%smpd_crop)
        call cline_refine3D%set('box_crop',   lpinfo(istage)%box_crop)
        ! iteration number
        call cline_refine3D%set('startit',                       iter)
        call cline_refine3D%set('which_iter',                    iter)
        ! dynamic control parameters
        call cline_refine3D%set('pgrp',                          pgrp)
        call cline_refine3D%set('refine',                      refine)
        call cline_refine3D%set('balance',                    balance)
        call cline_refine3D%set('trail_rec',                trail_rec)
        call cline_refine3D%set('lp_auto',                    lp_auto)
        if( lp_auto.eq.'yes' )then
        call cline_refine3D%set('lpstart',                    lpstart)
        call cline_refine3D%set('lpstop',                      lpstop)
        else
        call cline_refine3D%delete('lpstart')
        call cline_refine3D%delete('lpstop')
        endif
        call cline_refine3D%set('automsk',                    automsk)
        if( l_cavgs )then
        call cline_refine3D%set('cavgw',                        cavgw)
        endif
        ! phase control parameters
        call cline_refine3D%set('nspace',                     inspace)
        call cline_refine3D%set('maxits',                     imaxits)
        call cline_refine3D%set('trs',                            trs)
        call cline_refine3D%set('sh_first',                  sh_first)
        call cline_refine3D%set('prob_sh',                    prob_sh)
        call cline_refine3D%set('ml_reg',                      ml_reg)
        call cline_refine3D%set('icm',                            icm)
        call cline_refine3D%set('frac_best',                frac_best)
        call cline_refine3D%set('overlap',                    overlap)
        call cline_refine3D%set('fracsrch',                  fracsrch)
        if( l_cavgs )then
        call cline_refine3D%set('snr_noise_reg',        snr_noise_reg)
        call cline_refine3D%delete('update_frac') ! never on cavgs
        else
        call cline_refine3D%delete('snr_noise_reg')
        endif
        if( gaufreq > 0.)then
        call cline_refine3D%set('gauref',                      'yes')
        call cline_refine3D%set('gaufreq',                   gaufreq)
        else
        call cline_refine3D%delete('gauref')
        call cline_refine3D%delete('gaufreq')
        endif
        if( l_polar )then
            call cline_refine3D%set('center_type',           'params')
            call cline_refine3D%set('ref_type',              ref_type)
        endif
    end subroutine set_cline_refine3D

    subroutine exec_refine3D( istage, xrefine3D )
        integer,               intent(in)    :: istage
        class(commander_base), intent(inout) :: xrefine3D
        type(string) :: stage, str_state, vol_name, vol_pproc
        integer      :: state
        call cline_refine3D%delete('endit')
        call xrefine3D%execute_safe(cline_refine3D)
        call del_files(DIST_FBODY,      params_glob%nparts,ext='.dat')
        call del_files(ASSIGNMENT_FBODY,params_glob%nparts,ext='.dat')
        call del_file(DIST_FBODY//'.dat')
        call del_file(ASSIGNMENT_FBODY//'.dat')
        stage = '_stage_'//int2str(istage)
        if( l_polar )then
            vol_name = CAVGS_ITER_FBODY//int2str_pad(cline_refine3D%get_iarg('endit'),3)//params_glob%ext%to_char()
            call simple_copy_file(vol_name, string('cavgs')//stage//params_glob%ext)
        else
            do state = 1, params_glob%nstates
                str_state = int2str_pad(state,2)
                vol_name  = string(VOL_FBODY)//str_state//params_glob%ext
                vol_pproc = add2fbody(vol_name, params_glob%ext, PPROC_SUFFIX)
                if( file_exists(vol_name) ) call simple_copy_file(vol_name,  add2fbody(vol_name, params_glob%ext,stage))
                if( file_exists(vol_pproc)) call simple_copy_file(vol_pproc, add2fbody(vol_pproc,params_glob%ext,stage))
            enddo
        endif
    end subroutine exec_refine3D

    subroutine symmetrize( istage, spproj, projfile, xreconstruct3D )
        integer,                         intent(in)    :: istage
        class(sp_project),               intent(inout) :: spproj
        class(string),                   intent(in)    :: projfile
        class(commander_base), optional, intent(inout) :: xreconstruct3D
        type(commander_symmetrize_map) :: xsymmap
        type(cmdline)                  :: cline_asymrec, cline_symrec
        type(string) :: vol_iter, vol_sym
        real :: lpsym
        if( l_symran )then
            call se1%symrandomize(spproj%os_ptcl3D)
            call spproj%write_segment_inside('ptcl3D', projfile)
        endif
        if( l_srch4symaxis )then
            ! asymmetric/low symmetry reconstruction
            if( l_polar )then
                if( .not.present(xreconstruct3D) )then
                    THROW_HARD('Reconstructor required with polar=yes')
                endif
                cline_asymrec = cline_refine3D
                call cline_asymrec%set('prg',        'reconstruct3D')
                call cline_asymrec%set('mkdir',      'no')
                call cline_asymrec%set('projfile',   projfile)
                call cline_asymrec%set('pgrp',       params_glob%pgrp_start)
                call cline_asymrec%set('ml_reg',     'no') ! no ml reg for now
                call cline_asymrec%set('objfun',     'cc')
                call cline_asymrec%set('needs_sigma','no')
                call cline_asymrec%delete('which_iter')
                call cline_asymrec%delete('endit')
                call xreconstruct3D%execute_safe(cline_asymrec)
                vol_iter = 'asymmetric_map'//params_glob%ext%to_char()
                call simple_copy_file(string(VOL_FBODY)//int2str_pad(1,2)//params_glob%ext%to_char(), vol_iter)
                call cline_asymrec%kill
            else
                ! Volume from previous stage
                vol_iter = string(VOL_FBODY)//STR_STATE_GLOB//params_glob%ext%to_char()
            endif
            ! symmetry determination & map symmetrization
            if( .not. file_exists(vol_iter) ) THROW_HARD('input volume to map symmetrization does not exist')
            call cline_symmap%set('vol1', vol_iter)
            call cline_symmap%set('smpd', lpinfo(istage)%smpd_crop)
            call cline_symmap%set('box',  lpinfo(istage)%box_crop)
            vol_sym = 'symmetrized_map'//params_glob%ext%to_char()
            call cline_symmap%set('outvol', vol_sym)
            lpsym = max(LPSYMSRCH_LB,lpinfo(SYMSRCH_STAGE)%lp)
            call cline_symmap%set('lp', lpsym)
            write(logfhandle,'(A,F5.1)') '>>> DID SET MAP SYMMETRIZATION LOW-PASS LIMIT (IN A) TO: ', lpsym
            write(logfhandle,'(A)') '>>>'
            write(logfhandle,'(A)') '>>> MAP SYMMETRIZATION'
            write(logfhandle,'(A)') '>>>'
            call xsymmap%execute_safe(cline_symmap)
            call del_file('SYMAXIS_SEARCH_FINISHED')
            if( present(xreconstruct3D) )then
                ! symmetric reconstruction
                cline_symrec = cline_refine3D
                call cline_symrec%set('prg',        'reconstruct3D')
                call cline_symrec%set('mkdir',      'no')
                call cline_symrec%set('projfile',   projfile)
                call cline_symrec%set('pgrp',       params_glob%pgrp)
                call cline_symrec%set('which_iter', cline_refine3D%get_iarg('endit'))
                call cline_symrec%delete('endit')
                call xreconstruct3D%execute_safe(cline_symrec)
                vol_sym = VOL_FBODY//int2str_pad(1,2)//params_glob%ext%to_char()
                call simple_copy_file(vol_sym, string('symmetric_map')//params_glob%ext)
                call cline_symrec%kill
            endif
            call cline_refine3D%set('vol1', vol_sym)
        endif
    end subroutine symmetrize

    subroutine calc_start_rec( projfile, xreconstruct3D, istage )
        class(string),         intent(in)    :: projfile
        class(commander_base), intent(inout) :: xreconstruct3D
        integer,               intent(in)    :: istage
        type(commander_automask)       :: xautomask
        type(string)  :: str_state, vol, str, vol_even, vol_odd
        type(cmdline) :: cline_startrec, cline_automask
        integer       :: state
        cline_startrec = cline_refine3D
        call cline_startrec%set('prg',         'reconstruct3D')
        call cline_startrec%set('mkdir',       'no')
        call cline_startrec%set('projfile',    projfile)
        call cline_startrec%set('pgrp',        params_glob%pgrp)
        call cline_startrec%set('objfun',      'cc') ! ugly, but this is how it works in parameters 
        call cline_startrec%set('box_crop',    lpinfo(istage)%box_crop)
        call cline_startrec%set('projrec',     'no')
        call cline_startrec%delete('update_frac')    ! use all particles that have been updated
        call cline_startrec%delete('which_iter')
        call cline_startrec%delete('endit')
        call cline_startrec%delete('needs_sigma')
        call cline_startrec%delete('sigma_est')
        call cline_startrec%delete('automsk') ! no automask generated
        call cline_startrec%delete('mskfile') ! no masked FSC
        ! endif
        call xreconstruct3D%execute_safe(cline_startrec)
        do state = 1,params_glob%nstates
            ! rename volumes and update cline
            str_state = int2str_pad(state,2)
            vol       = string(VOL_FBODY)//str_state//params_glob%ext
            str       = string(STARTVOL_FBODY)//str_state//params_glob%ext
            call      simple_rename(vol, str)
            params_glob%vols(state) = str
            vol       = 'vol'//int2str(state)
            call      cline_refine3D%set(vol, str)
            vol_even  = string(VOL_FBODY)//str_state//'_even'//params_glob%ext%to_char()
            str       = string(STARTVOL_FBODY)//str_state//'_even_unfil'//params_glob%ext%to_char()
            call      simple_copy_file(vol_even, str)
            str       = string(STARTVOL_FBODY)//str_state//'_even'//params_glob%ext%to_char()
            call      simple_rename(vol_even, str)
            vol_odd   = string(VOL_FBODY)//str_state//'_odd' //params_glob%ext%to_char()
            str       = string(STARTVOL_FBODY)//str_state//'_odd_unfil'//params_glob%ext%to_char()
            call      simple_copy_file(vol_odd, str)
            str       = string(STARTVOL_FBODY)//str_state//'_odd'//params_glob%ext%to_char()
            call      simple_rename(vol_odd, str)
        enddo
        if( istage >= AUTOMSK_STAGE .and. l_automsk )then
            str_state = int2str_pad(1,2)
            vol_even  = string(STARTVOL_FBODY)//str_state//'_even'//params_glob%ext%to_char()
            vol_odd   = string(STARTVOL_FBODY)//str_state//'_odd'//params_glob%ext%to_char()
            call cline_automask%set('vol1', vol_odd)
            call cline_automask%set('vol2', vol_even)
            call cline_automask%set('smpd', lpinfo(istage)%smpd_crop)
            call cline_automask%set('amsklp', params_glob%amsklp)
            call cline_automask%set('automsk', 'yes')
            call cline_automask%set('mkdir',    'no')
            call cline_automask%set('nthr', params_glob%nthr)
            call xautomask%execute_safe(cline_automask)
            params_glob%mskfile = MSKVOL_FILE
            call cline_refine3D%set('mskfile', MSKVOL_FILE)
        endif
        call cline_startrec%kill
    end subroutine calc_start_rec

    ! Performs reconstruction at some set stages when polar=yes
    subroutine calc_rec4polar( xreconstruct3D, istage )
        use simple_class_frcs, only: class_frcs
        class(commander_base), intent(inout) :: xreconstruct3D
        integer,               intent(in)    :: istage
        real, allocatable :: fsc(:)
        type(cmdline)     :: cline_rec
        type(class_frcs)  :: frcs
        type(string)      :: src, dest, sstate, sstage, ext, pgrp
        integer :: i, inspace
        if( trim(params_glob%polar) /= 'yes' ) return
        if( (istage /= NSPACE_PHASE_POLAR(1)+1) .and. (istage /= NSPACE_PHASE_POLAR(2)+1) ) return
        ! Reconstruction
        pgrp = trim(params_glob%pgrp)
        if( istage <= SYMSRCH_STAGE ) pgrp = trim(params_glob%pgrp_start)
        cline_rec = cline_refine3D
        call cline_rec%set('prg',       'reconstruct3D')
        call cline_rec%set('mkdir',     'no')
        call cline_rec%set('projfile',  params_glob%projfile)
        call cline_rec%set('pgrp',      pgrp)
        call cline_rec%set('box_crop',  lpinfo(istage)%box_crop)
        call cline_rec%set('projrec',   'yes')
        call cline_rec%set('trail_rec', 'no')
        if( cline_rec%get_carg('ml_reg').ne.'yes' )then
            call cline_rec%set('objfun',    'cc')
        endif
        call cline_rec%delete('update_frac')
        call cline_rec%delete('endit')
        call cline_rec%delete('needs_sigma')
        call cline_rec%delete('automsk')
        call cline_rec%delete('mskfile')
        call xreconstruct3D%execute_safe(cline_rec)
        ! preserve volume (e/o will be overwritten in next iteration)
        sstate = int2str_pad(1,2)
        sstage = int2str_pad(istage-1,2)
        ext    = params_glob%ext
        src    = string(VOL_FBODY)//sstate//ext
        dest   = string(VOL_FBODY)//sstate//'_stage_'//sstage//ext
        call simple_rename(src, dest)
        ! Update refine3D command line
        call cline_refine3D%set('vol1', dest)
        ! Update FRCS
        fsc     = file2rarr(string(FSC_FBODY)//int2str_pad(1,2)//BIN_EXT)
        inspace = cline_refine3D%get_iarg('nspace')
        call frcs%new(inspace, lpinfo(istage)%box_crop, lpinfo(istage)%smpd_crop, 1)
        do i = 1,inspace
            call frcs%set_frc(i,fsc)
        enddo
        call frcs%write(string(FRCS_FILE))
        ! cleanup
        deallocate(fsc)
        call frcs%kill
        call cline_rec%kill
    end subroutine calc_rec4polar

    subroutine randomize_states( spproj, projfile, xreconstruct3D, istage )
        class(sp_project),     intent(inout) :: spproj
        class(string),         intent(in)    :: projfile
        class(commander_base), intent(inout) :: xreconstruct3D
        integer,               intent(in)    :: istage
        call spproj%read_segment('ptcl3D', projfile)
        call gen_labelling(spproj%os_ptcl3D, params_glob%nstates, 'squared_uniform')
        call spproj%write_segment_inside(params_glob%oritype, projfile)
        call cline_refine3D%set(     'nstates', params_glob%nstates)
        call cline_reconstruct3D%set('nstates', params_glob%nstates)
        call cline_postprocess%set(  'nstates', params_glob%nstates)
        call cline_reproject%set(    'nstates', params_glob%nstates)
        call calc_start_rec(projfile, xreconstruct3D, istage=istage)
    end subroutine randomize_states

    subroutine gen_ortho_reprojs4viz( spproj )
        type(sp_project), intent(inout) :: spproj
        type(string) :: str_state, fname
        type(image)  :: final_vol, reprojs
        integer      :: state, ifoo, ldim(3)
        real         :: smpd
        call spproj%read_segment('out', params_glob%projfile)
        do state = 1, params_glob%nstates
            if( .not.spproj%isthere_in_osout('vol', state) )cycle   ! empty-state case
            str_state = int2str_pad(state,2)
            fname = string(VOL_FBODY)//str_state//params_glob%ext
            if( .not. file_exists(fname) )cycle
            exit
        enddo
        call find_ldim_nptcls(fname, ldim, ifoo, smpd)
        call final_vol%new(ldim, smpd)
        do state = 1, params_glob%nstates
            str_state = int2str_pad(state,2)
            if( spproj%isthere_in_osout('vol', state) )then
                str_state = int2str_pad(state,2)
                fname     = VOL_FBODY//str_state%to_char()//params_glob%ext%to_char()
                if( .not. file_exists(fname) )cycle
                call final_vol%read(fname)
                call final_vol%generate_orthogonal_reprojs(reprojs)
                call reprojs%write_jpg(string('orthogonal_reprojs_state')//str_state//'.jpg')
                call reprojs%kill
            endif
        enddo
        call final_vol%kill
    end subroutine gen_ortho_reprojs4viz

    subroutine calc_final_rec( spproj, projfile, xreconstruct3D )
        class(sp_project),     intent(inout) :: spproj
        class(string),         intent(in)    :: projfile
        class(commander_base), intent(inout) :: xreconstruct3D
        type(string) :: str_state, vol_name
        integer      :: state, pop
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> RECONSTRUCTION AT ORIGINAL SAMPLING'
        write(logfhandle,'(A)') '>>>'
        call xreconstruct3D%execute_safe(cline_reconstruct3D)
        call spproj%read_segment('out', projfile)
        call spproj%read_segment('ptcl3D', projfile)
        do state = 1, params_glob%nstates
            pop = spproj%os_ptcl3D%get_pop(state, 'state')
            if( pop == 0 )cycle     ! empty-state case
            str_state = int2str_pad(state,2)
            vol_name  = string(VOL_FBODY)//str_state//params_glob%ext
            if( .not. file_exists(vol_name) )cycle
            call spproj%add_vol2os_out(vol_name, params_glob%smpd, state, 'vol', pop=pop)
            call spproj%add_fsc2os_out(string(FSC_FBODY)//str_state//BIN_EXT, state, params_glob%box)
        enddo
        call spproj%write_segment_inside('out', projfile)
    end subroutine calc_final_rec

    subroutine postprocess_final_rec( spproj )
        class(sp_project), intent(in) :: spproj
        type(commander_postprocess)   :: xpostprocess
        type(string) :: str_state, vol_name, vol_pproc, vol_pproc_mirr, vol_final
        integer :: state
        do state = 1, params_glob%nstates
            if( .not.spproj%isthere_in_osout('vol', state) )cycle ! empty-state case
            str_state      = int2str_pad(state,2)
            vol_name       = string(VOL_FBODY)//str_state//params_glob%ext  ! reconstruction from particles stored in project
            if( .not. file_exists(vol_name) )cycle
            call cline_postprocess%set('state', state)
            call xpostprocess%execute_safe(cline_postprocess)
            vol_pproc      = add2fbody(vol_name,params_glob%ext, PPROC_SUFFIX)
            vol_pproc_mirr = add2fbody(vol_name,params_glob%ext, PPROC_SUFFIX//MIRR_SUFFIX)
            vol_final      = string(REC_FBODY)//str_state//params_glob%ext
            if( file_exists(vol_name)       ) call simple_copy_file(vol_name,    vol_final)
            if( file_exists(vol_pproc)      ) call simple_rename(vol_pproc,      add2fbody(vol_final,params_glob%ext, PPROC_SUFFIX))
            if( file_exists(vol_pproc_mirr) ) call simple_rename(vol_pproc_mirr, add2fbody(vol_final,params_glob%ext, PPROC_SUFFIX//MIRR_SUFFIX))
        enddo
    end subroutine postprocess_final_rec

    ! create starting noise volume(s)
    subroutine generate_random_volumes( box, smpd, cline )
        integer,        intent(in)    :: box
        real,           intent(in)    :: smpd
        type(cmdline),  intent(inout) :: cline
        type(string) :: vol_name
        type(image)  :: noisevol, signal
        real         :: b
        integer      :: s
        call noisevol%new([box,box,box], smpd)
        select case(trim(params_glob%inivol))
            case('rand')
                ! random uniform distribution [0,1]
                ! resulting std dev in projections is ~box/10
                do s = 1, params_glob%nstates
                    call noisevol%ran
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'.mrc'
                    call cline%set('vol'//int2str(s), vol_name)
                    params_glob%vols(s) = vol_name
                    call noisevol%write(vol_name)
                    call noisevol%ran
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_even.mrc'
                    call noisevol%write(vol_name)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_even_unfil.mrc'
                    call noisevol%write(vol_name)
                    call noisevol%ran
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_odd.mrc'
                    call noisevol%write(vol_name)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_odd_unfil.mrc'
                    call noisevol%write(vol_name)
                end do
            case('rand_scaled')
                ! random uniform distribution [0,5/box]
                ! resulting std dev in projections is ~0.2
                b = 5.0/real(box)
                do s = 1, params_glob%nstates
                    call noisevol%ran(b=b)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'.mrc'
                    call cline%set('vol'//int2str(s), vol_name)
                    params_glob%vols(s) = vol_name
                    call noisevol%write(vol_name)
                    call noisevol%ran(b=b)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_even.mrc'
                    call noisevol%write(vol_name)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_even_unfil.mrc'
                    call noisevol%write(vol_name)
                    call noisevol%ran(b=b)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_odd.mrc'
                    call noisevol%write(vol_name)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_odd_unfil.mrc'
                    call noisevol%write(vol_name)
                end do
            case('sphere')
                ! random normal sphere N(0,5/box) + normal background N(0,5/box)
                b = 5.0/real(box)
                call signal%new([box,box,box], smpd)
                call signal%gauran(0.0, b)
                call signal%mask(0.25*real(box), 'soft', backgr=0.)
                do s = 1, params_glob%nstates
                    call noisevol%gauran(0., b)
                    call noisevol%add(signal)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'.mrc'
                    call cline%set('vol'//int2str(s), vol_name)
                    params_glob%vols(s) = vol_name
                    call noisevol%write(vol_name)
                    call noisevol%gauran(0., b)
                    call noisevol%add(signal)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_even.mrc'
                    call noisevol%write(vol_name)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_even_unfil.mrc'
                    call noisevol%write(vol_name)
                    call noisevol%gauran(0., b)
                    call noisevol%add(signal)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_odd.mrc'
                    call noisevol%write(vol_name)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_odd_unfil.mrc'
                    call noisevol%write(vol_name)
                end do
                call signal%kill
        end select
        call noisevol%kill
    end subroutine generate_random_volumes

end module simple_commanders_abinitio
