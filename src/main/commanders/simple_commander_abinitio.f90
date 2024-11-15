! concrete commander: high-level workflows
module simple_commander_abinitio
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_parameters,         only: parameters, params_glob
use simple_sp_project,         only: sp_project
use simple_stack_io,           only: stack_io
use simple_qsys_env,           only: qsys_env
use simple_commander_base,     only: commander_base
use simple_commander_volops,   only: reproject_commander, symaxis_search_commander, postprocess_commander, symmetrize_map_commander
use simple_commander_rec,      only: reconstruct3D_commander, reconstruct3D_commander_distr
use simple_commander_refine3D, only: refine3D_commander, refine3D_commander_distr
use simple_procimgstk,         only: shift_imgfile
use simple_image,              only: image
use simple_builder,            only: builder
use simple_class_frcs,         only: class_frcs
use simple_convergence,        only: convergence
use simple_cluster_seed,       only: gen_labelling
use simple_commander_euclid
use simple_euclid_sigma2
use simple_qsys_funs
use simple_decay_funs
implicit none

public :: abinitio3D_cavgs_commander, abinitio3D_commander, abinitio3D_parts_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: abinitio3D_cavgs_commander
    contains
    procedure :: execute => exec_abinitio3D_cavgs
end type abinitio3D_cavgs_commander

type, extends(commander_base) :: abinitio3D_commander
    contains
    procedure :: execute => exec_abinitio3D
end type abinitio3D_commander

type, extends(commander_base) :: abinitio3D_parts_commander
    contains
    procedure :: execute => exec_abinitio3D_parts
end type abinitio3D_parts_commander

! class constants
character(len=*), parameter :: REC_FBODY             = 'rec_final_state'
character(len=*), parameter :: STR_STATE_GLOB        = '01'
real,             parameter :: LPSTART_LB            = 10.
real,             parameter :: LPSTART_DEFAULT       = 20.
real,             parameter :: LPSTOP_LB             = 6.
real,             parameter :: CENLP_DEFAULT         = 30.
real,             parameter :: LPSYMSRCH_LB          = 12.
integer,          parameter :: NSTAGES               = 8
integer,          parameter :: NSTAGES_INI3D         = 4 ! # of ini3D stages used for initialization
integer,          parameter :: PHASES(3)             = [2,6,8]
integer,          parameter :: MAXITS(3)             = [20,17,15]
integer,          parameter :: MAXITS_GLOB           = 2*20 + 4*17 + 2*15
integer,          parameter :: NSPACE(3)             = [500,1000,2500]
integer,          parameter :: SYMSRCH_STAGE         = 3
integer,          parameter :: INCR_GREEDINESS_STAGE = 4 ! INCR_GREEDINESS_STAGE < PROBREFINE_STAGE must be true for it to have an effect
integer,          parameter :: PROBREFINE_STAGE      = 5
integer,          parameter :: ICM_STAGE             = PROBREFINE_STAGE
integer,          parameter :: TRAILREC_STAGE        = 7
! class variables
type(lp_crop_inf), allocatable :: lpinfo(:)
logical          :: l_srch4symaxis=.false., l_symran=.false., l_sym=.false., l_update_frac=.false.
logical          :: l_ml_reg=.true., l_icm_reg=.true., l_ini3D=.false., l_multistates = .false.
type(sym)        :: se1, se2
type(cmdline)    :: cline_refine3D, cline_symmap, cline_reconstruct3D, cline_postprocess, cline_reproject
real             :: update_frac = 1.0

contains

    !> for generation of an initial 3D model from class averages
    subroutine exec_abinitio3D_cavgs( self, cline )
        class(abinitio3D_cavgs_commander), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        character(len=*),      parameter :: work_projfile = 'abinitio3D_cavgs_tmpproj.simple'
        ! shared-mem commanders
        type(refine3D_commander)         :: xrefine3D
        type(reconstruct3D_commander)    :: xreconstruct3D
        type(reproject_commander)        :: xreproject
        ! other
        character(len=:),    allocatable :: stk, stkpath, orig_stk, shifted_stk, stk_even, stk_odd, ext
        integer,             allocatable :: states(:)
        type(ori)                        :: o, o_even, o_odd
        type(parameters)                 :: params
        type(ctfparams)                  :: ctfvars
        type(sp_project)                 :: spproj, work_proj
        type(image)                      :: img
        type(stack_io)                   :: stkio_r, stkio_r2, stkio_w
        character(len=STDLEN)            :: final_vol
        integer                          :: icls, ncavgs, cnt, even_ind, odd_ind, istage, nstages_ini3D
        call cline%set('objfun',    'euclid') ! use noise normalized Euclidean distances from the start
        call cline%set('sigma_est', 'global') ! obviously
        call cline%set('oritype',      'out') ! because cavgs are part of out segment
        call cline%set('bfac',            0.) ! because initial models should not be sharpened
        if( .not. cline%defined('mkdir')       ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('overlap')     ) call cline%set('overlap',     0.95)
        if( .not. cline%defined('prob_athres') ) call cline%set('prob_athres',  90.) ! reduces # failed runs on trpv1 from 4->2/10
        if( .not. cline%defined('cenlp')       ) call cline%set('cenlp', CENLP_DEFAULT)
        if( .not. cline%defined('imgkind')     ) call cline%set('imgkind',   'cavg') ! whether to use classes generated from 2D/3D
        ! make master parameters
        call params%new(cline)
        call cline%set('mkdir',       'no')   ! to avoid nested directory structure
        call cline%set('oritype', 'ptcl3D')   ! from now on we are in the ptcl3D segment, final report is in the cls3D segment
        ! set class global ML regularization flag
        l_ml_reg = .true.
        if( cline%defined('icm') )then
            l_ml_reg = params%l_ml_reg   
        endif
        ! set class global ICM regularization flag
        l_icm_reg = .true.
        if( cline%defined('icm') )then
            l_icm_reg = params%l_icm   
        endif
        ! set nstages_ini3D
        nstages_ini3D = NSTAGES
        if( cline%defined('nstages') )then
            nstages_ini3D = params%nstages   
        endif
        ! prepare class command lines
        call prep_class_command_lines(cline, work_projfile)
        ! set symmetry class variables
        call set_symmetry_class_vars
        ! read project
        call spproj%read(params%projfile)
        call spproj%update_projinfo(cline)
        call spproj%write_segment_inside('projinfo', params%projfile)
        ! set low-pass limits and downscaling info from FRCs
        call set_lplims_from_frcs(spproj)
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
        call spproj%get_cavgs_stk(stk, ncavgs, params%smpd, imgkind=params%imgkind, stkpath=stkpath)
        if(.not. file_exists(stk)) stk = trim(stkpath) // '/' // trim(stk)
        if(.not. file_exists(stk)) THROW_HARD('cavgs stk does not exist; simple_commander_abinitio')
        states          = nint(spproj%os_cls2D%get_all('state'))
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
            THROW_HARD('no class averages detected in project file: '//trim(params%projfile)//'; abinitio3D_cavgs')
        endif
        ! prepare a temporary project file
        work_proj%projinfo = spproj%projinfo
        work_proj%compenv  = spproj%compenv
        if( spproj%jobproc%get_noris() > 0 ) work_proj%jobproc = spproj%jobproc
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
        call work_proj%os_ptcl3D%delete_entry('stkind')
        call work_proj%os_ptcl3D%delete_entry('eo')
        params_glob%nptcls = ncavgs
        call spproj%os_cls3D%new(ncavgs, is_ptcl=.false.)
        do icls=1,ncavgs
            call spproj%os_cls3D%transfer_ori(icls, work_proj%os_ptcl3D, icls)
        enddo
        ! revert splitting
        call spproj%os_cls3D%set_all2single('stkind', 1)
        ! map the orientation parameters obtained for the clusters back to the particles
        call spproj%map2ptcls
        if( nstages_ini3D == NSTAGES )then ! produce validation info
            ! check even odd convergence
            call conv_eo(work_proj%os_ptcl3D)
            ! for visualization
            call gen_ortho_reprojs4viz
            ! calculate 3D reconstruction at original sampling
            call calc_final_rec(work_proj, work_projfile, xreconstruct3D)
            ! postprocess final 3D reconstruction
            call postprocess_final_rec
            ! add rec_final to os_out, one state assumed
            final_vol = trim(REC_FBODY)//STR_STATE_GLOB//params%ext
            call spproj%add_vol2os_out(final_vol, params%smpd, 1, 'vol_cavg')
            ! reprojections
            call spproj%os_cls3D%write('final_oris.txt')
            write(logfhandle,'(A)') '>>>'
            write(logfhandle,'(A)') '>>> RE-PROJECTION OF THE FINAL VOLUME'
            write(logfhandle,'(A)') '>>>'
            call xreproject%execute_safe(cline_reproject)
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
        endif
        ! write results (this needs to be a full write as multiple segments are updated)
        call spproj%write()
        ! end gracefully
        call img%kill
        call spproj%kill
        call o%kill
        call o_even%kill
        call o_odd%kill
        call work_proj%kill
        call del_file(work_projfile)
        call simple_rmdir(STKPARTSDIR)
        call simple_end('**** SIMPLE_ABINITIO3D_CAVGS NORMAL STOP ****')

        contains

            subroutine rndstart( cline )
                class(cmdline), intent(inout) :: cline
                call work_proj%os_ptcl3D%rnd_oris
                call work_proj%os_ptcl3D%zero_shifts
                call work_proj%write_segment_inside('ptcl3D', work_projfile)
                call cline%set('mkdir', 'no') ! to avoid nested dirs
                call cline%set('objfun', 'cc')
                call cline%set('silence_fsc', 'yes')
                call xreconstruct3D%execute_safe(cline)
                call cline%set('objfun', trim(params%objfun))
                call simple_copy_file('recvol_state01_even.mrc', 'startvol_even_unfil.mrc')
                call simple_copy_file('recvol_state01_odd.mrc',  'startvol_odd_unfil.mrc')
                call simple_rename(   'recvol_state01_even.mrc', 'startvol_even.mrc')
                call simple_rename(   'recvol_state01_odd.mrc',  'startvol_odd.mrc')
                call simple_rename(   'recvol_state01.mrc',      'startvol.mrc')
                call cline%set('vol1', 'startvol.mrc')
            end subroutine rndstart
    
            subroutine conv_eo( os )
                class(oris), intent(inout) :: os
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

    end subroutine exec_abinitio3D_cavgs

    !> for generation of an initial 3d model from particles
    subroutine exec_abinitio3D( self, cline )
        class(abinitio3D_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        ! commanders
        type(refine3D_commander_distr)         :: xrefine3D
        type(reconstruct3D_commander_distr)    :: xreconstruct3D_distr
        ! other
        real,                      parameter   :: UPDATE_FRAC_MAX = 0.9 !< to ensure fractional update is always on
        character(len=:),          allocatable :: vol_name
        real,                      allocatable :: rstates(:)
        integer,                   allocatable :: tmpinds(:), clsinds(:)
        type(class_sample),        allocatable :: clssmp(:)
        type(parameters)                       :: params
        type(sp_project)                       :: spproj
        type(image)                            :: noisevol
        integer :: istage, s, ncls, icls, nptcls_eff, i
        call cline%set('objfun',    'euclid') ! use noise normalized Euclidean distances from the start
        call cline%set('sigma_est', 'global') ! obviously
        if( .not. cline%defined('mkdir')       ) call cline%set('mkdir',         'yes')
        if( .not. cline%defined('overlap')     ) call cline%set('overlap',        0.95)
        if( .not. cline%defined('prob_athres') ) call cline%set('prob_athres',     10.)
        if( .not. cline%defined('center')      ) call cline%set('center',         'no')
        if( .not. cline%defined('cenlp')       ) call cline%set('cenlp', CENLP_DEFAULT)
        if( .not. cline%defined('oritype')     ) call cline%set('oritype',    'ptcl3D')
        if( .not. cline%defined('pgrp')        ) call cline%set('pgrp',           'c1')
        if( .not. cline%defined('pgrp_start')  ) call cline%set('pgrp_start',     'c1')
        if( .not. cline%defined('ptclw')       ) call cline%set('ptclw',          'no')
        if( .not. cline%defined('projrec')     ) call cline%set('projrec',       'yes')
        ! make master parameters
        call params%new(cline)
        call cline%set('mkdir', 'no')
        l_multistates = .false.
        if( cline%defined('nstates') )then
            if( params%nstates > 1  )then
                l_multistates = .true.
                call cline%set('projrec', 'no') ! not yet supported for multi-state
            endif
        endif
        ! read project
        call spproj%read(params%projfile)
        call spproj%update_projinfo(cline)
        call spproj%write_segment_inside('projinfo', params%projfile)
        ! provide initialization of 3D alignment using class averages?
        l_ini3D = .false.
        if( trim(params%cavg_ini).eq.'yes' )then
            call ini3D_from_cavgs(cline)
        endif
        ! initialization on class averages done outside this workflow (externally)?
        if( trim(params%cavg_ini_ext).eq.'yes' )then 
            ! check that ptcl3D field is not virgin
            if( spproj%is_virgin_field('ptcl3D') )then
                THROW_HARD('Prior 3D alignment required for abinitio workflow when cavg_ini_ext is set to yes')
            endif
            l_ini3D = .true.
        endif
        ! set class global ML regularization flag
        l_ml_reg = .true.
        if( cline%defined('icm') )then
            l_ml_reg = params%l_ml_reg   
        endif
        ! set class global ICM regularization flag
        l_icm_reg = .true.
        if( cline%defined('icm') )then
            l_icm_reg = params%l_icm
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
            l_update_frac = .true.
            update_frac   = 1.0
            nptcls_eff    = spproj%count_state_gt_zero()
            if( cline%defined('nsample') )then
                update_frac = real(params%nsample) / real(nptcls_eff)
            else if( cline%defined('update_frac') )then
                update_frac = params%update_frac
            else
                if( cline%defined('nsample_max') )then
                    update_frac = calc_update_frac(nptcls_eff, [NSAMPLE_MINMAX_DEFAULT(1),params%nsample_max])
                else
                    update_frac = calc_update_frac(nptcls_eff, NSAMPLE_MINMAX_DEFAULT)
                endif
            endif
            update_frac = min(UPDATE_FRAC_MAX, update_frac) ! to ensure fractional update is always on      
            ! generate a data structure for class sampling on disk
            ncls    = spproj%os_cls2D%get_noris()
            tmpinds = (/(icls,icls=1,ncls)/)
            rstates = spproj%os_cls2D%get_all('state')
            clsinds = pack(tmpinds, mask=rstates > 0.5)
            call spproj%os_ptcl2D%get_class_sample_stats(clsinds, clssmp)
            call write_class_samples(clssmp, CLASS_SAMPLING_FILE)
            deallocate(rstates, tmpinds, clsinds)
            call deallocate_class_samples(clssmp)
        endif
        ! set low-pass limits and downscaling info from FRCs
        call set_lplims_from_frcs(spproj)
        ! starting volume logics
        if( .not. cline%defined('vol1') )then
            if( .not. l_ini3D )then
                ! randomize projection directions
                select case(trim(params%oritype))
                    case('ptcl3D')
                        call spproj%os_ptcl3D%rnd_oris
                    case DEFAULT
                        THROW_HARD('Unsupported ORITYPE; exec_abinitio3D')
                end select
                ! randomize states
                if( l_multistates )then
                    call gen_labelling(spproj%os_ptcl3D, params%nstates, 'squared_uniform')
                endif
                call spproj%write_segment_inside(params%oritype, params%projfile)
                ! create noise starting volume(s)
                call noisevol%new([lpinfo(1)%box_crop,lpinfo(1)%box_crop,lpinfo(1)%box_crop], lpinfo(1)%smpd_crop)
                do s = 1, params%nstates
                    call noisevol%ran()
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'.mrc'
                    call cline_refine3D%set('vol'//int2str(s), vol_name)
                    params%vols(s) = vol_name
                    call noisevol%write(vol_name)
                    call noisevol%ran()
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_even.mrc'
                    call noisevol%write(vol_name)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_even_unfil.mrc'
                    call noisevol%write(vol_name)
                    call noisevol%ran()
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_odd.mrc'
                    call noisevol%write(vol_name)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_odd_unfil.mrc'
                    call noisevol%write(vol_name)
                end do
                call noisevol%kill
            else
                ! randomize states
                if( l_multistates )then
                    call gen_labelling(spproj%os_ptcl3D, params%nstates, 'squared_uniform')
                endif
                call spproj%write_segment_inside(params%oritype, params%projfile)
                ! create starting volume(s)
                call calc_start_rec(params%projfile, xreconstruct3D_distr)
            endif
        endif
        ! Frequency marching
        do istage = 1, NSTAGES
            write(logfhandle,'(A)')'>>>'
            write(logfhandle,'(A,I3,A9,F5.1)')'>>> STAGE ', istage,' WITH LP =', lpinfo(istage)%lp
            ! Preparation of command line for probabilistic search
            call set_cline_refine3D(istage, l_cavgs=.false.)
            if( lpinfo(istage)%l_autoscale )then
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',lpinfo(istage)%box_crop
            endif
            ! Probabilistic search
            call exec_refine3D(istage, xrefine3D)
            ! Symmetrization
            if( istage == SYMSRCH_STAGE )then
                call symmetrize(istage, spproj, params%projfile, xreconstruct3D_distr)
            endif
        enddo
        ! for visualization
        call gen_ortho_reprojs4viz
        ! calculate 3D reconstruction at original sampling
        call calc_final_rec(spproj, params%projfile, xreconstruct3D_distr)
        ! postprocess final 3D reconstruction
        call postprocess_final_rec
        ! cleanup
        call spproj%kill
        call qsys_cleanup
        call simple_end('**** SIMPLE_ABINITIO3D NORMAL STOP ****')
    end subroutine exec_abinitio3D
    
    subroutine exec_abinitio3D_parts( self, cline )
        use simple_commander_project,   only: new_project_commander, selection_commander
        use simple_exec_helpers,        only: gen_exec_cmd, async_exec
        use simple_commander_cluster2D, only: make_cavgs_commander_distr
        class(abinitio3D_parts_commander), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        ! commanders
        type(selection_commander)          :: xsel
        type(new_project_commander)        :: xnew_project
        type(make_cavgs_commander_distr)   :: xmake_cavgs_distr
        ! command lines
        type(cmdline)                      :: cline_split_bal, cline_new_proj, cline_mk_cavgs, cline_abinitio3D
        ! other vars
        character(len=:),      allocatable :: cmd
        character(len=STDLEN), allocatable :: projnames(:), projfnames(:)
        character(len=LONGSTRLEN)          :: cwd
        type(parameters)                   :: params
        type(sp_project)                   :: spproj
        integer                            :: iproj
        logical, parameter                 :: L_USE_CAVGS = .false.
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        ! make master parameters
        call params%new(cline)
        call cline%set('mkdir', 'no')
        ! provide initialization of 3D alignment using class averages?
        l_ini3D = .false.
        if( trim(params%cavg_ini).eq.'yes' )then
            call ini3D_from_cavgs(cline)
            call cline%set('cavg_ini_ext', 'yes')
        else
            call cline%set('cavg_ini_ext', 'no')
        endif
        ! split stack so it does not happen downstream
        call spproj%read(params%projfile)
        call spproj%split_stk(max(params%nparts,params%nparts_per_part*params%nparts))
        call spproj%kill
        ! conduct balanced split
        cline_split_bal = cline
        call cline_split_bal%set('balance', 'yes')
        call cline_split_bal%set('oritype', 'cls2D')
        call xsel%execute(cline_split_bal)
        ! make projects
        allocate(projfnames(params%nparts), projnames(params%nparts))
        do iproj = 1, params%nparts
            projnames(iproj)  = BALPROJPARTFBODY//int2str(iproj)
            projfnames(iproj) = BALPROJPARTFBODY//int2str(iproj)//'.simple'
            call cline_new_proj%set('projname', trim(projnames(iproj)))
            call cline_new_proj%set('projfile', trim(projfnames(iproj)))
            call xnew_project%execute_safe(cline_new_proj)
            call chdir('../')
            call del_file(trim(projfnames(iproj)))
        end do
        if( L_USE_CAVGS )then
            ! make class averages
            do iproj = 1, params%nparts
                call chdir('./'//trim(projnames(iproj)))
                call simple_getcwd(cwd)
                write(logfhandle,'(A)') 'CWD: '//trim(cwd)
                call cline_mk_cavgs%set('prg',      'make_cavgs')
                call cline_mk_cavgs%set('projfile',  trim(projfnames(iproj)))
                call cline_mk_cavgs%set('mkdir',    'no') ! to avoid nested directory structure
                call cline_mk_cavgs%set('refs',     'cavgs_'//BALPROJPARTFBODY//int2str(iproj)//'.mrc')
                call cline_mk_cavgs%set('nparts',    params%nparts_per_part)
                call cline_mk_cavgs%set('nthr',      params%nthr)
                ! cmd = gen_exec_cmd(cline_mk_cavgs, 'simple_exec')
                ! call exec_cmdline(cmd)
                call xmake_cavgs_distr%execute_safe(cline_mk_cavgs)
                call chdir('../')
            end do
            ! execute independent jobs for cross validation asynchronously
            do iproj = 1, params%nparts
                call chdir('./'//trim(projnames(iproj)))
                call simple_getcwd(cwd)
                write(logfhandle,'(A)') 'CWD: '//trim(cwd)
                cline_abinitio3D = cline
                call cline_abinitio3D%delete('nparts_per_part')
                call cline_abinitio3D%delete('nparts')
                call cline_abinitio3D%set('prg',        'abinitio3D_cavgs')
                call cline_abinitio3D%set('projfile',    trim(projfnames(iproj)))
                call cline_abinitio3D%set('mkdir',      'no') ! to avoid nested directory structure
                ! cmd = gen_exec_cmd(cline_abinitio3D, 'simple_exec', 'ABINITIO3D_CAVGS')
                call async_exec(cline_abinitio3D, 'simple_exec', 'ABINITIO3D_CAVGS')
                call chdir('../')
            end do
        else
            ! execute independent jobs for cross validation asynchronously
            do iproj = 1, params%nparts
                call chdir('./'//trim(projnames(iproj)))
                call simple_getcwd(cwd)
                write(logfhandle,'(A)') 'CWD: '//trim(cwd)
                cline_abinitio3D = cline
                call cline_abinitio3D%delete('nparts_per_part')
                call cline_abinitio3D%delete('nthr_ini3D')
                call cline_abinitio3D%set('cavg_ini',   'no')
                call cline_abinitio3D%set('prg',        'abinitio3D')
                call cline_abinitio3D%set('projfile',    trim(projfnames(iproj)))
                call cline_abinitio3D%set('mkdir',      'no') ! to avoid nested directory structure
                call cline_abinitio3D%set('nparts',      params%nparts_per_part)
                call cline_abinitio3D%set('update_frac', 1.0) ! maximal nsample
                ! cmd = gen_exec_cmd(cline_abinitio3D, 'simple_exec', 'ABINITIO3D')
                call async_exec(cline_abinitio3D, 'simple_exec', 'ABINITIO3D')
                call chdir('../')
            end do
        endif
        call simple_end('**** SIMPLE_ABINITIO3D_PARTS NORMAL STOP ****')
    end subroutine exec_abinitio3D_parts

    ! private helper routines

    subroutine prep_class_command_lines( cline, projfile )
        class(cmdline),   intent(in) :: cline
        character(len=*), intent(in) :: projfile
        cline_refine3D      = cline
        cline_symmap        = cline
        cline_reconstruct3D = cline
        cline_postprocess   = cline
        cline_reproject     = cline
        ! refine3D
        call cline_refine3D%set('prg',                         'refine3D')
        call cline_refine3D%set('pgrp',            trim(params_glob%pgrp))
        call cline_refine3D%set('projfile',                trim(projfile))
        ! symmetrization
        call cline_symmap%set('prg',                     'symmetrize_map')
        call cline_symmap%set('pgrp',              trim(params_glob%pgrp))
        call cline_symmap%set('projfile',                  trim(projfile))
        if( .not. cline_symmap%defined('cenlp') )then
        call cline_symmap%set('cenlp',                      CENLP_DEFAULT)
        endif
        call cline_symmap%set('hp',                        params_glob%hp)
        ! re-reconstruct volume
        call cline_reconstruct3D%set('prg',               'reconstruct3D')
        call cline_reconstruct3D%set('box',               params_glob%box)
        call cline_reconstruct3D%set('smpd',             params_glob%smpd)
        call cline_reconstruct3D%set('projfile',           trim(projfile))
        call cline_reconstruct3D%set('pgrp',       trim(params_glob%pgrp))
        call cline_reconstruct3D%set('ml_reg',                       'no')
        call cline_reconstruct3D%set('needs_sigma',                  'no')
        call cline_reconstruct3D%set('objfun',                       'cc')
        ! no fractional update
        call cline_reconstruct3D%delete('update_frac')
        ! individual particles reconstruction
        call cline_reconstruct3D%set('projrec', 'no')
        ! postprocess volume
        call cline_postprocess%set('prg',                   'postprocess')
        call cline_postprocess%set('projfile',             trim(projfile))
        call cline_postprocess%set('mkdir',                          'no')
        call cline_postprocess%set('imgkind',                       'vol')
        call cline_postprocess%delete('bfac') ! sharpen final map
        call cline_postprocess%delete('lp')   ! to obtain optimal filtration
        ! re-project volume
        call cline_reproject%set('prg',                       'reproject')
        call cline_reproject%set('vol1', REC_FBODY//STR_STATE_GLOB//PPROC_SUFFIX//params_glob%ext)
        call cline_reproject%set('pgrp',           trim(params_glob%pgrp))
        call cline_reproject%set('outstk',     'reprojs'//params_glob%ext)
        call cline_reproject%set('smpd',                 params_glob%smpd)
        call cline_reproject%set('box',                   params_glob%box)
        call cline_reproject%set('oritab',               'final_oris.txt')
        call cline_reproject%delete('projfile')
    end subroutine prep_class_command_lines

    subroutine set_symmetry_class_vars
        l_srch4symaxis = trim(params_glob%pgrp) .ne. trim(params_glob%pgrp_start)
        l_symran       = .false.
        l_sym          = l_srch4symaxis
        if( params_glob%pgrp_start.ne.'c1' .or. params_glob%pgrp.ne.'c1' )then
            se1 = sym(params_glob%pgrp_start)
            se2 = sym(params_glob%pgrp)
            if(se1%get_nsym() > se2%get_nsym())then
                ! ensure se2 is a subgroup of se1
                if( .not. se1%has_subgrp(params_glob%pgrp) )THROW_HARD('Incompatible symmetry groups; exec_abinitio3D')
                ! set flag for symmetry randomisation
                ! in case we are moving from a higher to lower group
                l_symran = .true.
            else if( se2%get_nsym() > se1%get_nsym() )then
                ! ensure se1 is a subgroup of se2
                if( .not. se2%has_subgrp(params_glob%pgrp_start) )THROW_HARD('Incompatible symmetry groups; exec_abinitio3D')
            endif
        endif
    end subroutine set_symmetry_class_vars

    subroutine set_lplims_from_frcs( spproj )
        class(sp_project), intent(inout) :: spproj
        character(len=:),  allocatable   :: frcs_fname, stk, imgkind, stkpath
        real,              allocatable   :: frcs_avg(:)
        integer,           allocatable   :: states(:)
        type(class_frcs) :: clsfrcs
        real             :: smpd, lpfinal
        integer          :: filtsz, ncavgs
        ! retrieve FRC info
        call spproj%get_frcs(frcs_fname, 'frc2D', fail=.false.)
        ! work out low-pass limits and downscaling parameters
        params_glob%frcs = trim(frcs_fname)
        call clsfrcs%read(frcs_fname)
        filtsz = clsfrcs%get_filtsz()
        allocate(frcs_avg(filtsz), source=0.)
        states = nint(spproj%os_cls2D%get_all('state'))
        if( params_glob%frc_weight .eq. 'yes' )then
            call clsfrcs%avg_frc_getter(frcs_avg, states, cur_oris=spproj%os_ptcl2D)
        else
            call clsfrcs%avg_frc_getter(frcs_avg, states)
        endif
        if( allocated(lpinfo) ) deallocate(lpinfo)
        allocate(lpinfo(NSTAGES))
        lpfinal = max(LPSTOP_LB,calc_lplim_final_stage(3))
        call lpstages(params_glob%box, NSTAGES, frcs_avg, params_glob%smpd, LPSTART_LB, LPSTART_DEFAULT, lpfinal, lpinfo, verbose=.true.)
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
        type(abinitio3D_cavgs_commander) :: xini3D
        type(cmdline)                    :: cline_ini3D
        character(len=:), allocatable    :: cavgs_stk
        type(str4arr),    allocatable    :: files_that_stay(:)
        character(len=*), parameter      :: INI3D_DIR='abinitio3D_cavgs/'
        integer :: ncavgs
        cline_ini3D = cline
        call cline_ini3D%set('nstages', NSTAGES_INI3D)
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
        allocate(files_that_stay(5))
        files_that_stay(1)%str = basename(trim(params_glob%projfile))
        files_that_stay(2)%str = 'cavgs'
        files_that_stay(3)%str = 'nice_'
        files_that_stay(4)%str = 'frcs'
        files_that_stay(5)%str = 'ABINITIO3D'
        ! make the move
        call move_files_in_cwd(INI3D_DIR, files_that_stay)
        ! flag completion
        l_ini3D = .true.
    end subroutine ini3D_from_cavgs

    subroutine set_cline_refine3D( istage, l_cavgs )
        integer,          intent(in)  :: istage
        logical,          intent(in)  :: l_cavgs
        character(len=:), allocatable :: silence_fsc, sh_first, prob_sh, ml_reg
        character(len=:), allocatable :: refine, icm, trail_rec, pgrp, balance
        integer :: iphase, iter, inspace, imaxits
        real    :: trs, snr_noise_reg, greediness, frac_best, overlap, fracsrch
        ! iteration number bookkeeping
        if( cline_refine3D%defined('endit') )then
            iter = cline_refine3D%get_iarg('endit')
        else
            iter = 0
        endif
        iter = iter + 1
        ! symmetry
        pgrp = trim(params_glob%pgrp)
        if( l_srch4symaxis )then
            if( istage <= SYMSRCH_STAGE )then
                ! need to replace original point-group flag with c1/pgrp_start
                pgrp = trim(params_glob%pgrp_start)
            endif
        endif
        ! refinement mode
        if( istage < PROBREFINE_STAGE )then
            refine    = 'shc_smpl'
            prob_sh   = 'no'
        else
            refine    = 'prob'
            prob_sh   = 'yes'
        endif
        ! ICM regularization
        if( istage < ICM_STAGE )then
            icm       = 'no'
        else
            icm       = 'yes'
        endif
        ! balance
        if( l_update_frac )then
            balance   = 'yes'
        else
            balance   = 'no'
        endif
        ! trailing reconstruction
        if( istage >= TRAILREC_STAGE .and. l_update_frac )then
            trail_rec = 'yes'
        else
            trail_rec = 'no'
        endif
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
                inspace       = NSPACE(1)
                imaxits       = MAXITS(1)
                silence_fsc   = 'yes'
                trs           = 0.
                sh_first      = 'no'
                ml_reg        = 'no'
                greediness    = 2.0 ! completely greedy balanced sampling based on objective function value
                frac_best     = 1.0 ! means it does not control sampling
                overlap       = 0.90
                fracsrch      = 90.
                snr_noise_reg = 2.0
            case(2)
                inspace       = NSPACE(2)
                imaxits       = MAXITS(2)
                silence_fsc   = 'yes'
                trs           = lpinfo(istage)%trslim
                sh_first      = 'yes'
                ml_reg        = 'yes'
                if( istage >= INCR_GREEDINESS_STAGE )then
                greediness    = 1.0 ! sample first half of each class as the best ones and the rest randomly
                frac_best     = 0.5 ! means sampling is done from top-ranking 50% particles in class
                else
                greediness    = 2.0 ! completely greedy balanced sampling based on objective function value
                frac_best     = 1.0 ! means it does not control sampling
                endif
                overlap       = 0.95
                fracsrch      = 95.
                snr_noise_reg = 4.0
            case(3)
                inspace       = NSPACE(3)
                imaxits       = MAXITS(3)
                if( l_cavgs )then
                silence_fsc   = 'yes'
                else
                silence_fsc   = 'no'
                endif
                trs           = lpinfo(istage)%trslim
                sh_first      = 'yes'
                ml_reg        = 'yes'
                greediness    = 0.0  ! completely random balanced sampling (only class assignment matters)
                frac_best     = 0.85 ! means sampling is done from top-ranking 85% particles in class
                overlap       = 0.99
                fracsrch      = 99.
                snr_noise_reg = 6.0
        end select
        ! overrride regularization parameters
        if( .not. l_ml_reg  ) ml_reg = 'no'
        if( .not. l_icm_reg ) icm    = 'no'
        ! command line update
        call cline_refine3D%set('prg',                     'refine3D')
        ! class global control parameters
        call cline_refine3D%set('update_frac',            update_frac)
        call cline_refine3D%set('lp',               lpinfo(istage)%lp)
        call cline_refine3D%set('smpd_crop', lpinfo(istage)%smpd_crop)
        call cline_refine3D%set('box_crop',   lpinfo(istage)%box_crop)
        call cline_refine3D%set('maxits_glob',            MAXITS_GLOB)
        ! iteration number
        call cline_refine3D%set('startit',                       iter)
        call cline_refine3D%set('which_iter',                    iter)
        ! dynamic control parameters
        call cline_refine3D%set('pgrp',                          pgrp)
        call cline_refine3D%set('refine',                      refine)
        call cline_refine3D%set('balance',                    balance)
        call cline_refine3D%set('trail_rec',                trail_rec)
        ! phase control parameters
        call cline_refine3D%set('nspace',                     inspace)
        call cline_refine3D%set('maxits',                     imaxits)
        call cline_refine3D%set('silence_fsc',            silence_fsc)
        call cline_refine3D%set('trs',                            trs)
        call cline_refine3D%set('sh_first',                  sh_first)
        call cline_refine3D%set('prob_sh',                    prob_sh)
        call cline_refine3D%set('ml_reg',                      ml_reg)
        call cline_refine3D%set('icm',                            icm)
        call cline_refine3D%set('greediness',              greediness)
        call cline_refine3D%set('frac_best',                frac_best)
        call cline_refine3D%set('overlap',                    overlap)
        call cline_refine3D%set('fracsrch',                  fracsrch)
        if( l_cavgs )then
            call cline_refine3D%set('snr_noise_reg',    snr_noise_reg)
        endif
    end subroutine set_cline_refine3D

    subroutine exec_refine3D( istage, xrefine3D )
        integer,               intent(in)    :: istage
        class(commander_base), intent(inout) :: xrefine3D
        character(len=:),      allocatable   :: stage, str_state, vol_name, vol_pproc
        integer :: state
        call cline_refine3D%delete('endit')
        call xrefine3D%execute_safe(cline_refine3D)
        call del_files(DIST_FBODY,      params_glob%nparts,ext='.dat')
        call del_files(ASSIGNMENT_FBODY,params_glob%nparts,ext='.dat')
        call del_file(DIST_FBODY      //'.dat')
        call del_file(ASSIGNMENT_FBODY//'.dat')
        stage = '_stage_'//int2str(istage)
        do state = 1, params_glob%nstates
            str_state = int2str_pad(state,2)
            vol_name  = VOL_FBODY//str_state//params_glob%ext
            vol_pproc = add2fbody(vol_name, params_glob%ext, PPROC_SUFFIX)
            if( file_exists(vol_name) ) call simple_copy_file(vol_name,  add2fbody(vol_name, params_glob%ext,stage))
            if( file_exists(vol_pproc)) call simple_copy_file(vol_pproc, add2fbody(vol_pproc,params_glob%ext,stage))
        enddo
    end subroutine exec_refine3D

    subroutine symmetrize( istage, spproj, projfile, xreconstruct3D )
        integer,                         intent(in)    :: istage
        class(sp_project),               intent(inout) :: spproj
        character(len=*),                intent(in)    :: projfile
        class(commander_base), optional, intent(inout) :: xreconstruct3D
        type(symmetrize_map_commander) :: xsymmap
        type(cmdline)                  :: cline_symrec
        character(len=:),  allocatable :: vol_iter, vol_sym
        real :: lpsym
        if( l_symran )then
            call se1%symrandomize(spproj%os_ptcl3D)
            call spproj%write_segment_inside('ptcl3D', projfile)
        endif
        if( l_srch4symaxis )then
            ! symmetry determination & map symmetrization
            vol_iter = VOL_FBODY//STR_STATE_GLOB//params_glob%ext
            if( .not. file_exists(vol_iter) ) THROW_HARD('input volume to map symmetrization does not exist')
            call cline_symmap%set('vol1', vol_iter)
            call cline_symmap%set('smpd', lpinfo(istage)%smpd_crop)
            call cline_symmap%set('box',  lpinfo(istage)%box_crop)
            vol_sym = 'symmetrized_map'//params_glob%ext
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
                vol_sym = VOL_FBODY//int2str_pad(1,2)//params_glob%ext
                call simple_copy_file(vol_sym,  'symmetric_map'//params_glob%ext)
                call cline_symrec%kill
            endif
            call cline_refine3D%set('vol1', vol_sym)
        endif
    end subroutine symmetrize

    subroutine calc_start_rec( projfile, xreconstruct3D )
        character(len=*),      intent(in)    :: projfile
        class(commander_base), intent(inout) :: xreconstruct3D
        character(len=:),  allocatable :: str_state, vol, str, vol_even, vol_odd
        type(cmdline) :: cline_startrec
        integer       :: state
        cline_startrec = cline_refine3D
        call cline_startrec%set('prg',         'reconstruct3D')
        call cline_startrec%set('mkdir',       'no')
        call cline_startrec%set('projfile',    projfile)
        call cline_startrec%set('pgrp',        params_glob%pgrp)
        call cline_startrec%set('objfun',      'cc') ! ugly, but this is how it works in parameters 
        call cline_startrec%set('silence_fsc', 'yes')
        call cline_startrec%delete('which_iter')
        call cline_startrec%delete('endit')
        call cline_startrec%delete('needs_sigma')
        call cline_startrec%delete('sigma_est')
        call xreconstruct3D%execute_safe(cline_startrec)
        do state = 1,params_glob%nstates
            ! rename volumes and update cline
            str_state = int2str_pad(state,2)
            vol       = trim(VOL_FBODY)//trim(str_state)//params_glob%ext
            str       = trim(STARTVOL_FBODY)//trim(str_state)//params_glob%ext
            call      simple_rename( trim(vol), trim(str) )
            params_glob%vols(state) = trim(str)
            vol       = 'vol'//trim(int2str(state))
            call      cline_refine3D%set( trim(vol), trim(str) )
            vol_even  = trim(VOL_FBODY)//trim(str_state)//'_even'//params_glob%ext
            str       = trim(STARTVOL_FBODY)//trim(str_state)//'_even_unfil'//params_glob%ext
            call      simple_copy_file( trim(vol_even), trim(str) )
            str       = trim(STARTVOL_FBODY)//trim(str_state)//'_even'//params_glob%ext
            call      simple_rename( trim(vol_even), trim(str) )
            vol_odd   = trim(VOL_FBODY)//trim(str_state)//'_odd' //params_glob%ext
            str       = trim(STARTVOL_FBODY)//trim(str_state)//'_odd_unfil'//params_glob%ext
            call      simple_copy_file( trim(vol_odd), trim(str) )
            str       = trim(STARTVOL_FBODY)//trim(str_state)//'_odd'//params_glob%ext
            call      simple_rename( trim(vol_odd), trim(str) )
        enddo
        call cline_startrec%kill
    end subroutine calc_start_rec

    subroutine gen_ortho_reprojs4viz
        character(len=:), allocatable :: str_state
        type(image) :: final_vol, reprojs
        integer     :: state, box_crop
        box_crop = lpinfo(NSTAGES)%box_crop 
        call final_vol%new([box_crop,box_crop,box_crop],lpinfo(NSTAGES)%smpd_crop)
        do state = 1, params_glob%nstates
            str_state = int2str_pad(state,2)
            call final_vol%read(VOL_FBODY//str_state//params_glob%ext)
            call final_vol%generate_orthogonal_reprojs(reprojs)
            call reprojs%write_jpg('orthogonal_reprojs_state'//str_state//'.jpg')
            call reprojs%kill
        enddo
        call final_vol%kill
    end subroutine gen_ortho_reprojs4viz

    subroutine calc_final_rec( spproj, projfile, xreconstruct3D )
        class(sp_project),     intent(inout) :: spproj
        character(len=*),      intent(in)    :: projfile
        class(commander_base), intent(inout) :: xreconstruct3D
        character(len=:),      allocatable   :: str_state, vol_name
        integer :: state
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> RECONSTRUCTION AT ORIGINAL SAMPLING'
        write(logfhandle,'(A)') '>>>'
        call xreconstruct3D%execute_safe(cline_reconstruct3D)
        call spproj%read_segment('out', projfile)
        do state = 1, params_glob%nstates
            str_state = int2str_pad(state,2)
            vol_name  = VOL_FBODY//str_state//params_glob%ext
            call spproj%add_vol2os_out(vol_name, params_glob%smpd, state, 'vol')
            call spproj%add_fsc2os_out(FSC_FBODY//str_state//BIN_EXT, state, params_glob%box)
        enddo
        call spproj%write_segment_inside('out', projfile)
    end subroutine calc_final_rec

    subroutine postprocess_final_rec
        type(postprocess_commander)   :: xpostprocess
        character(len=:), allocatable :: str_state, vol_name, vol_pproc, vol_pproc_mirr, vol_final
        integer :: state
        do state = 1, params_glob%nstates
            call cline_postprocess%set('state', state)
            call xpostprocess%execute_safe(cline_postprocess)
        enddo
        do state = 1, params_glob%nstates
            str_state      = int2str_pad(state,2)
            vol_name       = VOL_FBODY//str_state//params_glob%ext ! reconstruction from particles stored in project
            vol_pproc      = add2fbody(vol_name,params_glob%ext,PPROC_SUFFIX)
            vol_pproc_mirr = add2fbody(vol_name,params_glob%ext,PPROC_SUFFIX//MIRR_SUFFIX)
            vol_final      = REC_FBODY//str_state//params_glob%ext
            if( file_exists(vol_name)       ) call simple_copy_file(vol_name,    vol_final)
            if( file_exists(vol_pproc)      ) call simple_rename(vol_pproc,      add2fbody(vol_final,params_glob%ext,PPROC_SUFFIX))
            if( file_exists(vol_pproc_mirr) ) call simple_rename(vol_pproc_mirr, add2fbody(vol_final,params_glob%ext,PPROC_SUFFIX//MIRR_SUFFIX))
        enddo
    end subroutine postprocess_final_rec

end module simple_commander_abinitio
