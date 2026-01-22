!@descr: commanders for abinitio 3D reconstruction in single- and multi-particle mode
module simple_commanders_abinitio
use simple_commander_module_api
use simple_abinitio_config
use simple_abinitio_utils
use simple_procimgstk,          only: shift_imgfile
use simple_commanders_volops,   only: commander_reproject
use simple_commanders_refine3D, only: commander_refine3D, commander_refine3D_distr
use simple_commanders_rec,      only: commander_reconstruct3D, commander_reconstruct3D_distr
use simple_cluster_seed,        only: gen_labelling
use simple_decay_funs,          only: calc_update_frac_dyn, calc_update_frac
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
        if( .not. cline%defined('mkdir')            ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('objfun')           ) call cline%set('objfun',  'euclid') ! use noise normalized Euclidean distances from the start
        if( .not. cline%defined('overlap')          ) call cline%set('overlap',     0.95)
        if( .not. cline%defined('prob_athres')      ) call cline%set('prob_athres',  90.) ! reduces # failed runs on trpv1 from 4->2/10
        if( .not. cline%defined('cenlp')            ) call cline%set('cenlp', CENLP_DEFAULT)
        if( .not. cline%defined('imgkind')          ) call cline%set('imgkind',   'cavg')
        if( .not. cline%defined('lp_auto')          ) call cline%set('lp_auto',    'yes')
        if( .not. cline%defined('noise_norm')       ) call cline%set('noise_norm',  'no')
        if( .not. cline%defined('cavgw')            ) call cline%set('cavgw',       'no')
        if( .not. cline%defined('lpstart')          ) call cline%set('lpstart',      20.)
        if( .not. cline%defined('lpstop')           ) call cline%set('lpstop',        8.)
        if( .not. cline%defined('ref_type')         ) call cline%set('ref_type', 'comlin_noself')
        if( .not. cline%defined('gauref_last_stage')) call cline%set('gauref_last_stage', GAUREF_LAST_STAGE)
        if( .not. cline%defined('gauref')           ) call cline%set('gauref',     'yes')
        ! make master parameters
        call params%new(cline)
        call cline%set('mkdir',       'no')   ! to avoid nested directory structure
        call cline%set('oritype', 'ptcl3D')   ! from now on we are in the ptcl3D segment, final report is in the cls3D segment
        ! set work projfile
        work_projfile = 'abinitio3D_cavgs_tmpproj.simple'
        ! set class global lp_auto flag for low-pass limit estimation
        l_lpauto = .true.
        if( cline%defined('lp_auto') ) l_lpauto = params%l_lpauto
        ! Polar representation
        if( params%l_polar )then
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
            ! Reconstruction for polar representation
            if( l_polar )then
                call calc_rec4polar( xreconstruct3D, istage, work_projfile )
            endif
            ! Probabilistic search
            call exec_refine3D(istage, xrefine3D) 
            ! Symmetrization
            if( istage == SYMSRCH_STAGE )then
                call symmetrize(istage, work_proj, work_projfile, xreconstruct3D)
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
        if( .not. cline%defined('inivol')              ) call cline%set('inivol',                          'sphere')
        if( .not. cline%defined('maxits_between')      ) call cline%set('maxits_between',            MAXITS_BETWEEN)
        if( .not. cline%defined('gauref_last_stage')   ) call cline%set('gauref_last_stage',      GAUREF_LAST_STAGE)
        if( .not. cline%defined('gauref')              ) call cline%set('gauref',                             'yes')
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
        if( params%l_polar )then
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
                    THROW_WARN('automasking not supported for modes other than multivol_mode.eq.single, turning automasking off')
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
            call spproj%os_ptcl3D%delete_3Dalignment(keepshifts=.true.)
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
            if( l_polar )then
                if( l_ini3D .and. (istage==start_stage) )then
                    ! reconstruction has been performed above
                else
                    call calc_rec4polar( xreconstruct3D_distr, istage )
                endif
            endif
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

    contains

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

end module simple_commanders_abinitio
