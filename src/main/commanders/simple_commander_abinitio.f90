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
use simple_commander_euclid
use simple_euclid_sigma2
use simple_qsys_funs
use simple_decay_funs
implicit none

public :: initial_3Dmodel_commander, abinitio_3Dmodel_commander, abinitio_3Dmodel_parts_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: initial_3Dmodel_commander
    contains
    procedure :: execute => exec_initial_3Dmodel
end type initial_3Dmodel_commander

type, extends(commander_base) :: abinitio_3Dmodel_commander
    contains
    procedure :: execute => exec_abinitio_3Dmodel
end type abinitio_3Dmodel_commander

type, extends(commander_base) :: abinitio_3Dmodel_parts_commander
    contains
    procedure :: execute => exec_abinitio_3Dmodel_parts
end type abinitio_3Dmodel_parts_commander

type, extends(commander_base) :: abinitio_3Dmodel2_commander
    contains
    procedure :: execute => exec_abinitio_3Dmodel2
end type abinitio_3Dmodel2_commander

! class constants
character(len=*), parameter :: REC_FBODY             = 'rec_final_state'
character(len=*), parameter :: STR_STATE_GLOB        = '01'
real,             parameter :: LPSTART_LB            = 10.
real,             parameter :: LPSTART_DEFAULT       = 20.
real,             parameter :: LPSTOP_LB             = 6.
real,             parameter :: CENLP_DEFAULT         = 30.
real,             parameter :: LPSYMSRCH_LB          = 12.
integer,          parameter :: NSTAGES               = 8
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
logical          :: l_icm_reg=.true., l_ini3D=.false., l_greediness_given = .false.
type(sym)        :: se1,se2
type(cmdline)    :: cline_refine3D, cline_symmap, cline_reconstruct3D, cline_postprocess, cline_reproject
real             :: update_frac   = 1.0

contains

    !> for generation of an initial 3D model from class averages
    subroutine exec_initial_3Dmodel( self, cline )
        class(initial_3Dmodel_commander), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        character(len=*),      parameter :: work_projfile = 'initial_3Dmodel_tmpproj.simple'
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
            THROW_HARD('no class averages detected in project file: '//trim(params%projfile)//'; initial_3Dmodel')
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
        call simple_end('**** SIMPLE_INITIAL_3DMODEL NORMAL STOP ****')

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

    end subroutine exec_initial_3Dmodel

    !> for generation of an initial 3d model from particles
    subroutine exec_abinitio_3Dmodel( self, cline )
        class(abinitio_3Dmodel_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        ! commanders
        type(initial_3Dmodel_commander)        :: xini3D
        type(refine3D_commander_distr)         :: xrefine3D
        type(reconstruct3D_commander_distr)    :: xreconstruct3D_distr
        ! command lines
        type(cmdline)                          :: cline_ini3D
        ! other
        real,                      parameter   :: UPDATE_FRAC_MAX = 0.9 !< to ensure fractional update is always on
        integer,                   parameter   :: NSTAGES_INI3D   = 4   !< # of ini3D stages
        character(len=*),          parameter   :: INI3D_DIR='ini3D_on_cavgs/'
        character(len=:),          allocatable :: vol_name
        real,                      allocatable :: rstates(:)
        integer,                   allocatable :: tmpinds(:), clsinds(:)
        type(class_sample),        allocatable :: clssmp(:)
        character(len=LONGSTRLEN), allocatable :: file_list(:)
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
        ! provide initialization of 3D alignment using class averages?
        l_ini3D = .false.
        if( trim(params%cavg_ini).eq.'yes' )then 
            cline_ini3D = cline
            call cline_ini3D%set('nstages', NSTAGES_INI3D)
            if( cline%defined('nthr_ini3D') )then
                call cline_ini3D%set('nthr', params%nthr_ini3D)
                call cline_ini3D%delete('nthr_ini3D')
            endif
            call cline_ini3D%delete('nparts')
            call cline_ini3D%delete('projrec')
            call cline_ini3D%delete('oritype')
            call cline_ini3D%delete('imgkind')
            call cline_ini3D%delete('prob_athres')
            call xini3D%execute_safe(cline_ini3D)
            ! update point-group symmetry
            call cline%set('pgrp_start', params%pgrp)
            params%pgrp_start = params%pgrp
            ! stash away files
            call simple_list_files('*', file_list)
            call simple_mkdir(INI3D_DIR)
            do i = 1, size(file_list)
                if( str_has_substr(trim(file_list(i)), basename(trim(params%projfile))) )then
                    ! this stays
                else if( str_has_substr(trim(file_list(i)), 'nice_') )then
                    ! this stays
                else
                    call simple_rename(trim(file_list(i)), INI3D_DIR//trim(file_list(i)))
                endif
            end do
            deallocate(file_list)
            ! flag completion
            l_ini3D = .true.
        endif
        ! set class global ICM regularization flag
        l_icm_reg = .true.
        if( cline%defined('icm') )then
            l_icm_reg = params%l_icm
        endif
        ! set greediness flag
        l_greediness_given = cline%defined('greediness')
        ! prepare class command lines
        call prep_class_command_lines(cline, params%projfile)
        ! set symmetry class variables
        call set_symmetry_class_vars
        ! read project
        call spproj%read(params%projfile)
        call spproj%update_projinfo(cline)
        call spproj%write_segment_inside('projinfo', params%projfile)
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
                        THROW_HARD('Unsupported ORITYPE; exec_abinitio_3Dmodel')
                end select
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
        call simple_end('**** SIMPLE_ABINITIO_3DMODEL NORMAL STOP ****')
    end subroutine exec_abinitio_3Dmodel

    subroutine exec_abinitio_3Dmodel_parts( self, cline )
        use simple_commander_project, only: new_project_commander, selection_commander
        class(abinitio_3Dmodel_parts_commander), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        ! commanders
        type(selection_commander)   :: xsel
        type(new_project_commander) :: xnew_project
        ! command lines
        type(cmdline)               :: cline_split_bal, cline_new_proj      
        ! other vars
        type(parameters) :: params
        integer          :: iproj
        character(len=STDLEN), allocatable :: projnames(:), projfnames(:)
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        ! command-line inputs
        ! REQUIRED: nparts    
        ! OPTIONAL: nptcls_per_part

        ! make master parameters
        call params%new(cline)
        call cline%set('mkdir', 'no')
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
            call xnew_project%execute(cline_new_proj)
            call chdir('../')
            call del_file(trim(projfnames(iproj)))
        end do
        
        ! simple_exec prg=make_cavgs nparts=1 nthr=4 refs=cavgs_balanced_part1.mrc

    end subroutine exec_abinitio_3Dmodel_parts

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
                if( .not. se1%has_subgrp(params_glob%pgrp) )THROW_HARD('Incompatible symmetry groups; exec_abinitio_3Dmodel')
                ! set flag for symmetry randomisation
                ! in case we are moving from a higher to lower group
                l_symran = .true.
            else if( se2%get_nsym() > se1%get_nsym() )then
                ! ensure se1 is a subgroup of se2
                if( .not. se2%has_subgrp(params_glob%pgrp_start) )THROW_HARD('Incompatible symmetry groups; exec_abinitio_3Dmodel')
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
        if( .not.file_exists(frcs_fname) )then
            ! 08/24 This is a backwards compatibility patch to account for error in metadata
            ! on exit of streaming related to GUI directory structure (now fixed and cf above get_cavgs_stk).
            ! Will need to harmonize (move to absolute path?).
            call spproj%get_cavgs_stk(stk, ncavgs, smpd, imgkind=imgkind, stkpath=stkpath)
            frcs_fname = trim(stkpath)//'/'//trim(frcs_fname)
            if( .not.file_exists(frcs_fname) )then
                THROW_HARD('the project file does not contain an FRCs file, which is required')
            endif
        endif
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
        ! override phased greediness
        if( l_greediness_given ) greediness = params_glob%greediness
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
        ! type(class_frcs)              :: clsfrcs
        type(image)                   :: vol_even, vol_odd, tmpvol, vol
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
        call cline_reconstruct3D%set('box',           params%box)
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
            call xcalc_pspec_distr%execute_safe(cline_calc_pspec_distr)
            call cline_reconstruct3D%set('which_iter', 1)
            call cline_reconstruct3D%set('ml_reg',     'yes')
            call cline_reconstruct3D%set('needs_sigma','yes')
            call cline_reconstruct3D%set('objfun',     'euclid')
        endif
        call xreconstruct3D_distr%execute_safe(cline_reconstruct3D)
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
        call xreconstruct3D_distr%execute_safe(cline_reconstruct3D)
        vol_str = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
        call spproj%read_segment('out',params%projfile)
        call spproj%add_vol2os_out(vol_str, params%smpd, 1, 'vol')
        call spproj%add_fsc2os_out(FSC_FBODY//str_state//trim(BIN_EXT), 1, params%box)
        call spproj%write_segment_inside('out',params%projfile)
        ! post-processing
        call cline_postprocess%delete('lp')
        call cline_postprocess%set('state', 1)
        call xpostprocess%execute_safe(cline_postprocess)
        vol_pproc      = add2fbody(vol_str,params%ext,PPROC_SUFFIX)
        vol_pproc_mirr = add2fbody(vol_str,params%ext,trim(PPROC_SUFFIX)//trim(MIRR_SUFFIX))
        call simple_rename(vol_str, trim(REC_FBODY)//trim(str_state)//trim(params%ext))
        if(file_exists(vol_pproc)     ) call simple_rename(vol_pproc,      add2fbody(vol_str,params%ext,PPROC_SUFFIX))
        if(file_exists(vol_pproc_mirr)) call simple_rename(vol_pproc_mirr, add2fbody(vol_str,params%ext,PPROC_SUFFIX//MIRR_SUFFIX))
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
                    call xsymsrch%execute_safe(cline_symsrch)
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
