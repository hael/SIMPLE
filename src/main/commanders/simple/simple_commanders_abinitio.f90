!@descr: abinitio 3D reconstruction in single- and multi-particle mode
module simple_commanders_abinitio
use simple_commanders_api
use simple_abinitio_utils
use simple_procimgstk,           only: shift_imgfile
use simple_commanders_reproject, only: commander_reproject
use simple_commanders_refine3D,  only: commander_refine3D, commander_refine3D
use simple_commanders_rec,       only: commander_rec3D, commander_rec3D
use simple_cluster_seed,         only: gen_labelling
use simple_refine3D_fnames,      only: refine3D_startvol_fname, refine3D_startvol_half_fname, &
    &refine3D_state_vol_fname, refine3D_state_halfvol_fname
implicit none

public :: commander_abinitio3D_cavgs, commander_abinitio3D_cavg_sort, commander_abinitio3D, commander_multivol_assign
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_abinitio3D_cavgs
    contains
    procedure :: execute => exec_abinitio3D_cavgs
end type commander_abinitio3D_cavgs

type, extends(commander_base) :: commander_abinitio3D_cavg_sort
    contains
    procedure :: execute => exec_abinitio3D_cavg_sort
end type commander_abinitio3D_cavg_sort

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
        type(commander_refine3D)  :: xrefine3D
        type(commander_rec3D)     :: xrec3D
        type(commander_reproject) :: xreproject
        ! other
        type(string)              :: stk, orig_stk, shifted_stk, stk_even, stk_odd, ext
        integer, allocatable      :: states(:)
        type(ori)                 :: o, o_even, o_odd
        type(parameters)          :: params
        type(ctfparams)           :: ctfvars
        type(sp_project)          :: spproj, work_proj
        type(image)               :: img
        type(stack_io)            :: stkio_r, stkio_r2, stkio_w
        type(string)              :: final_vol, work_projfile
        integer                   :: icls, ncavgs, cnt, even_ind, odd_ind, istage, nstages_ini3D, s
        integer                   :: cavg_ldim(3), cavg_nimgs
        real                      :: cavg_smpd
        if( cline%defined('part') )then
            THROW_HARD('abinitio3D_cavgs distributed execution is master-only; remove part from command line')
        endif
        call cline%set('sigma_est', 'global') ! obviously
        call cline%set('oritype',      'out') ! because cavgs are part of out segment
        call cline%set('bfac',            0.) ! because initial models should not be sharpened
        call cline%set('filt_mode',   'none') ! no fancy filtering for cavgs route
        call cline%set('automsk',       'no') ! no envelope masking for cavgs route
        call cline%set('nu_refine',     'no') ! no nonuniform refinement for cavgs route
        if( .not. cline%defined('mkdir')            ) call cline%set('mkdir',                      'yes')
        if( .not. cline%defined('objfun')           ) call cline%set('objfun',                  'euclid') ! noise normalized Euclidean distances from the start
        if( .not. cline%defined('overlap')          ) call cline%set('overlap',                     0.95)
        if( .not. cline%defined('prob_athres')      ) call cline%set('prob_athres',                  90.) ! reduces # failed runs on trpv1 from 4->2/10
        if( .not. cline%defined('cenlp')            ) call cline%set('cenlp',   abinitio_cenlp_default())
        if( .not. cline%defined('imgkind')          ) call cline%set('imgkind',                   'cavg')
        if( .not. cline%defined('filt_mode')        ) call cline%set('filt_mode',                 'none')
        if( .not. cline%defined('noise_norm')       ) call cline%set('noise_norm',                  'no')
        if( .not. cline%defined('lpstart')          ) call cline%set('lpstart', abinitio_lpstart_ini3D())
        if( .not. cline%defined('lpstop')           ) call cline%set('lpstop',   abinitio_lpstop_ini3D())
        if( .not. cline%defined('gauref')           ) call cline%set('gauref',                     'yes')
        ! make master parameters
        call params%new(cline)
        call cline%set('mkdir',       'no')   ! to avoid nested directory structure
        call cline%set('oritype', 'ptcl3D')   ! from now on we are in the ptcl3D segment, final report is in the cls3D segment
        ! set work projfile
        work_projfile = 'abinitio3D_cavgs_tmpproj.simple'
        ! set class global filtering flags for staged refine3D policy
        l_nonuniform = .false.
        ! set nstages_ini3D
        nstages_ini3D = abinitio_nstages_ini3D_max()
        if( cline%defined('nstages') )then
            nstages_ini3D = min(abinitio_nstages_ini3D_max(),params%nstages)
        endif
        nstages_refine3D = nstages_ini3D
        ! prepare class command lines
        call prep_class_command_lines(params, cline, work_projfile)
        ! set symmetry class variables
        call set_symmetry_class_vars(params)
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
                call set_lplims_from_frcs(params, spproj, l_cavgs=.true., lpstart=params%lpstart, lpstop=params%lpstop)
            else if( cline%defined('lpstart') )then
                call set_lplims_from_frcs(params, spproj, l_cavgs=.true., lpstart=params%lpstart)
            else if( cline%defined('lpstop') )then
                call set_lplims_from_frcs(params, spproj, l_cavgs=.true., lpstop=params%lpstop)
            else
                call set_lplims_from_frcs(params, spproj, l_cavgs=.true.)
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
        params%nptcls = 2 * ncavgs
        call configure_cavgs_distributed_clines
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
        params%nptcls = work_proj%get_nptcls()
        call configure_cavgs_distributed_clines
        call work_proj%write()
        ! Frequency marching
        call set_cline_refine3D(params, 1, l_cavgs=.true.)
        call rndstart(cline_refine3D)
        do istage = 1, nstages_ini3D
            write(logfhandle,'(A)')'>>>'
            write(logfhandle,'(A,I3,A9,F5.1)')'>>> STAGE ', istage,' WITH LP =', lpinfo(istage)%lp
            ! Preparation of command line for probabilistic search
            call set_cline_refine3D(params, istage, l_cavgs=.true.)
            if( lpinfo(istage)%l_autoscale )then
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',lpinfo(istage)%box_crop
            endif
            ! Probabilistic search
            call exec_refine3D(params, istage, xrefine3D)
            ! Symmetrization
            if( istage == abinitio_symsrch_stage() )then
                call symmetrize(params, istage, work_proj, work_projfile, xrec3D)
            endif
        end do
        ! update original cls3D segment
        call work_proj%read_segment('ptcl3D', work_projfile)
        call work_proj%read_segment('out',    work_projfile)
        call work_proj%os_ptcl3D%delete_entry('stkind')
        call work_proj%os_ptcl3D%delete_entry('eo')
        params%nptcls = ncavgs
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
                call spproj%os_cls3D%set_euler(icls, work_proj%os_ptcl3D%get_euler(cnt))
                call spproj%os_cls3D%set_shift(icls, work_proj%os_ptcl3D%get_2Dshift(cnt))
                call spproj%os_cls3D%set_state(icls, work_proj%os_ptcl3D%get_state(cnt))
            endif
        enddo
        call spproj%os_cls3D%set_all2single('stkind', 1)    ! revert splitting
        ! map the orientation parameters obtained for the clusters back to the particles
        call spproj%map2ptcls
        if( nstages_ini3D == abinitio_nstages_ini3D_max() )then ! produce validation info
            call find_ldim_nptcls(orig_stk, cavg_ldim, cavg_nimgs)
            cavg_smpd = params%smpd
            if( cavg_nimgs < ncavgs ) THROW_HARD('fewer images in cavgs stack than expected; abinitio3D_cavgs')
            ! check even odd convergence
            if( params%nstates > 1 ) call conv_eo_states(work_proj%os_ptcl3D)
            call conv_eo(work_proj%os_ptcl3D)
            ! calculate 3D reconstruction at original sampling
            call calc_final_rec(params, work_proj, work_projfile, xrec3D, l_postprocess=.false.)
            ! final raw and low-pass diagnostic 3D reconstruction outputs
            call write_final_rec_outputs(params, work_proj, lpinfo(nstages_ini3D)%lp)
            ! add rec_final to os_out
            do s = 1,params%nstates
                if( .not.work_proj%isthere_in_osout('vol', s) )cycle
                final_vol = abinitio_rec_fbody()//int2str_pad(s,2)//MRC_EXT
                if( file_exists(final_vol) )then
                    call spproj%add_vol2os_out(final_vol, cavg_smpd, s, 'vol_cavg')
                endif
            enddo
            ! reprojections
            call spproj%os_cls3D%write(string('final_oris.txt'))
            write(logfhandle,'(A)') '>>>'
            write(logfhandle,'(A)') '>>> RE-PROJECTION OF THE FINAL VOLUME'
            write(logfhandle,'(A)') '>>>'
            do s = 1,params%nstates
                if( .not.work_proj%isthere_in_osout('vol', s) )cycle
                call cline_reproject%set('vol'//int2str(s), abinitio_rec_fbody()//int2str_pad(s,2)//LP_SUFFIX//MRC_EXT)
            enddo
            call cline_reproject%set('box',  cavg_ldim(1))
            call cline_reproject%set('smpd', cavg_smpd)
            call cline_reproject%delete('box_crop')
            call cline_reproject%delete('smpd_crop')
            call xreproject%execute(cline_reproject)
            ! write alternated stack
            call img%new([cavg_ldim(1),cavg_ldim(1),1],     cavg_smpd)
            call stkio_r%open(orig_stk,                     cavg_smpd, 'read',                                    bufsz=500)
            call stkio_r2%open(string('reprojs.mrc'),       cavg_smpd, 'read',                                    bufsz=500)
            call stkio_w%open(string('cavgs_reprojs.mrc'),  cavg_smpd, 'write', box=cavg_ldim(1), is_ft=.false., bufsz=500)
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
            call shift_imgfile(orig_stk, shifted_stk, spproj%os_cls3D, cavg_smpd)
            ! add shifted stack to project
            call spproj%add_cavgs2os_out(simple_abspath(shifted_stk), cavg_smpd, 'cavg_shifted')
        endif
        ! write results (this needs to be a full write as multiple segments are updated)
        call spproj%write()
        ! rank classes based on agreement to volume (after writing)
        if( nstages_ini3D == abinitio_nstages_ini3D_max() )then
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
        call simple_end('**** SIMPLE_ABINITIO3D_CAVGS NORMAL STOP ****', &
            verbose_exit=trim(params%verbose_exit).eq.'yes', verbose_exit_fname=params%verbose_exit_fname)
        contains

            subroutine rndstart( cline )
                class(cmdline), intent(inout) :: cline
                type(string)  :: src, dest
                type(cmdline) :: local_cline_rec
                integer :: s
                call work_proj%os_ptcl3D%rnd_oris
                call work_proj%os_ptcl3D%zero_shifts
                if( params%nstates > 1 )then
                    call gen_labelling(work_proj%os_ptcl3D, params%nstates, 'uniform')
                endif
                call work_proj%write_segment_inside('ptcl3D', work_projfile)
                local_cline_rec = cline
                ! Distributed rec3D schedules workers from PRG, so do not inherit refine3D here.
                call local_cline_rec%set('prg',   'reconstruct3D')
                call local_cline_rec%set('mkdir', 'no') ! to avoid nested dirs
                call local_cline_rec%set('objfun', 'cc')
                call xrec3D%execute(local_cline_rec)
                do s = 1,params%nstates
                    src   = refine3D_state_vol_fname(s)
                    dest  = refine3D_startvol_fname(s)
                    call simple_rename(src, dest)
                    ! updates refine3D command line with new volume
                    call cline%set('vol'//int2str(s), dest)
                    src   = refine3D_state_halfvol_fname(s, 'even')
                    dest  = refine3D_startvol_half_fname(s, 'even', unfil=.true.)
                    call simple_copy_file(src, dest)
                    dest  = refine3D_startvol_half_fname(s, 'even')
                    call simple_rename(src, dest)
                    src   = refine3D_state_halfvol_fname(s, 'odd')
                    dest  = refine3D_startvol_half_fname(s, 'odd', unfil=.true.)
                    call simple_copy_file(src, dest)
                    dest  = refine3D_startvol_half_fname(s, 'odd')
                    call simple_rename(src, dest)
                enddo
                call local_cline_rec%kill
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
                call xrank_cavgs%execute(cline_rank_cavgs)
                call cline_rank_cavgs%kill
            end subroutine rank_cavgs

            subroutine configure_cavgs_distributed_clines
                integer :: nparts_eff
                if( .not. cline%defined('nparts') ) return
                nparts_eff = min(params%nparts, max(1, params%nptcls))
                if( nparts_eff < params%nparts )then
                    write(logfhandle,'(A,I0,A,I0)') '>>> REDUCING NPARTS FROM ', params%nparts, &
                        ' TO THE NUMBER OF EVEN/ODD CLASS AVERAGE ENTRIES: ', nparts_eff
                endif
                params%nparts = nparts_eff
                params%numlen = len(int2str(params%nparts))
                call cline%set('nparts', params%nparts)
                call cline%set('numlen', params%numlen)
                ! Only refinement/reconstruction are distributed in this workflow.
                call sync_distributed_child(cline_refine3D)
                call sync_distributed_child(cline_reconstruct3D)
                call strip_distributed_child(cline_symmap)
                call strip_distributed_child(cline_reproject)
            end subroutine configure_cavgs_distributed_clines

            subroutine sync_distributed_child( child_cline )
                type(cmdline), intent(inout) :: child_cline
                call child_cline%set('nparts', params%nparts)
                call child_cline%set('numlen', params%numlen)
            end subroutine sync_distributed_child

            subroutine strip_distributed_child( child_cline )
                type(cmdline), intent(inout) :: child_cline
                call child_cline%delete('nparts')
                call child_cline%delete('numlen')
            end subroutine strip_distributed_child

    end subroutine exec_abinitio3D_cavgs

    !> sort class averages by consensus over restarted two-state abinitio3D_cavgs runs
    subroutine exec_abinitio3D_cavg_sort( self, cline )
        use simple_cavg_quality_analysis, only: evaluate_cavg_quality, write_cavg_quality_feature_table
        use simple_cavg_quality_model,    only: CAVG_QUALITY_MODEL_CHUNK_DEFAULT, cavg_quality_model
        use simple_cavg_quality_types,    only: cavg_quality_result
        use simple_imgarr_utils,          only: read_cavgs_into_imgarr, dealloc_imgarr
        class(commander_abinitio3D_cavg_sort), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline
        character(len=*), parameter :: RESTART_DIR_FBODY = 'abinitio3D_cavg_sort_restart_'
        character(len=*), parameter :: RESTART_DONE      = 'ABINITIO3D_CAVG_SORT_FINISHED'
        character(len=*), parameter :: CONSENSUS_REPORT  = 'abinitio3D_cavg_sort_consensus.txt'
        integer,          parameter :: SORT_NSTATES      = 2
        integer,          parameter :: SORT_NSTAGES      = 2
        type(parameters)          :: params
        type(qsys_env)            :: qenv
        type(sp_project)          :: spproj, restart_proj
        type(cavg_quality_model)  :: model
        type(cavg_quality_result) :: quality
        type(image), allocatable  :: cavg_imgs(:)
        type(string)              :: cwd, cwd_run, projbase, folder, restart_projfile, done_file
        type(string), allocatable :: restart_projfiles(:), done_files(:)
        integer, allocatable      :: restart_labels(:,:), mapped_labels(:,:), consensus(:), final_states(:)
        integer, allocatable      :: original_states(:), quality_auto_states(:), votes(:,:)
        real,    allocatable      :: quality_scores(:)
        integer                   :: ncls, irestart
        if( cline%defined('part') )then
            THROW_HARD('abinitio3D_cavg_sort is master-only; remove part from command line')
        endif
        if( cline%defined('nstates') .and. cline%get_iarg('nstates') /= SORT_NSTATES )then
            THROW_HARD('abinitio3D_cavg_sort requires nstates=2')
        endif
        if( cline%defined('nstages') .and. cline%get_iarg('nstages') /= SORT_NSTAGES )then
            THROW_HARD('abinitio3D_cavg_sort always runs abinitio3D_cavgs through nstages=2')
        endif
        call cline%set('nstates',  SORT_NSTATES)
        call cline%set('nstages',  SORT_NSTAGES)
        if( .not.cline%defined('nrestarts')     ) call cline%set('nrestarts', 3)
        if( .not.cline%defined('mkdir')         ) call cline%set('mkdir', 'yes')
        if( .not.cline%defined('quality_model') ) call cline%set('quality_model', CAVG_QUALITY_MODEL_CHUNK_DEFAULT)
        call params%new(cline)
        if( params%nrestarts < 1 ) THROW_HARD('abinitio3D_cavg_sort requires nrestarts >= 1')
        call spproj%read(params%projfile)
        ncls = spproj%os_cls2D%get_noris()
        if( ncls == 0 ) THROW_HARD('abinitio3D_cavg_sort: project has no cls2D entries')
        original_states = spproj%os_cls2D%get_all_asint('state')
        if( size(original_states) /= ncls ) THROW_HARD('abinitio3D_cavg_sort: invalid cls2D state array')
        allocate(restart_labels(params%nrestarts,ncls), source=0)
        allocate(mapped_labels(params%nrestarts,ncls),  source=0)
        allocate(votes(SORT_NSTATES,ncls),              source=0)
        allocate(consensus(ncls),                       source=0)
        allocate(final_states(ncls),                    source=0)
        allocate(restart_projfiles(params%nrestarts))
        allocate(done_files(params%nrestarts))
        call init_quality_model
        call evaluate_quality
        call submit_restarts
        call read_restart_labels
        call map_state_correspondence
        call build_consensus
        call assign_good_bad_states
        call spproj%map_cavgs_selection(final_states)
        call annotate_project
        call spproj%write(params%projfile)
        call write_consensus_report
        call cleanup
        call simple_end('**** SIMPLE_ABINITIO3D_CAVG_SORT NORMAL STOP ****', &
            verbose_exit=trim(params%verbose_exit).eq.'yes', verbose_exit_fname=params%verbose_exit_fname)

    contains

        subroutine init_quality_model
            if( trim(params%quality_model) == '' )then
                call model%init_preset(CAVG_QUALITY_MODEL_CHUNK_DEFAULT)
            else
                call model%init_preset(params%quality_model)
            endif
            if( cline%defined('infile') .and. trim(params%infile%to_char()) /= '' ) call model%read(params%infile%to_char())
        end subroutine init_quality_model

        subroutine evaluate_quality
            cavg_imgs = read_cavgs_into_imgarr(spproj)
            if( size(cavg_imgs) /= ncls ) THROW_HARD('abinitio3D_cavg_sort: # cavgs /= # cls2D entries')
            call evaluate_cavg_quality(cavg_imgs, spproj%os_cls2D, params%mskdiam, quality, model)
            quality_scores      = quality%scores
            quality_auto_states = quality%states
            write(logfhandle,'(A,A)') '>>> CAVG SORT QUALITY MODEL         : ', trim(model%name)
            write(logfhandle,'(A,A)') '>>> CAVG SORT QUALITY MODEL CONTEXT : ', trim(model%context)
        end subroutine evaluate_quality

        subroutine submit_restarts
            type(cmdline) :: cline_restart
            call simple_getcwd(cwd)
            CWD_GLOB = cwd%to_char()
            projbase = basename(params%projfile)
            call qenv%new(params, 1, exec_bin=string('simple_exec'))
            do irestart = 1, params%nrestarts
                folder           = RESTART_DIR_FBODY//int2str_pad(irestart, 3)
                restart_projfile = folder//'/'//projbase
                done_file        = folder//'/'//RESTART_DONE
                call simple_mkdir(folder)
                call simple_copy_file(params%projfile, restart_projfile)
                if( file_exists(done_file) ) call del_file(done_file)
                restart_projfiles(irestart) = simple_abspath(restart_projfile, check_exists=.false.)
                done_files(irestart)        = simple_abspath(done_file,        check_exists=.false.)
                cline_restart = cline
                call cline_restart%set('prg',                'abinitio3D_cavgs')
                call cline_restart%set('projfile',           basename(restart_projfile))
                call cline_restart%set('mkdir',              'no')
                call cline_restart%set('nstates',            SORT_NSTATES)
                call cline_restart%set('nstages',            SORT_NSTAGES)
                call cline_restart%set('verbose_exit',       'yes')
                call cline_restart%set('verbose_exit_fname', RESTART_DONE)
                call cline_restart%delete('nrestarts')
                call cline_restart%delete('nparts')
                call cline_restart%delete('numlen')
                call cline_restart%delete('dir_exec')
                call cline_restart%delete('outdir')
                call cline_restart%delete('quality_mode')
                call cline_restart%delete('quality_model')
                call cline_restart%delete('filetab')
                call cline_restart%delete('fname')
                call cline_restart%delete('infile')
                call simple_chdir(folder)
                call simple_getcwd(cwd_run)
                CWD_GLOB = cwd_run%to_char()
                call qenv%exec_simple_prg_in_queue_async(cline_restart, &
                    string('abinitio3D_cavg_sort_script'), string('abinitio3D_cavg_sort.log'))
                call simple_chdir(cwd)
                CWD_GLOB = cwd%to_char()
                call cline_restart%kill
                write(logfhandle,'(A,I0,A,A)') '>>> SUBMITTED ABINITIO3D_CAVGS RESTART ', irestart, &
                    ' IN ', folder%to_char()
            enddo
            call qsys_watcher(done_files)
            do irestart = 1, params%nrestarts
                if( .not.file_exists(done_files(irestart)) )then
                    write(logfhandle,'(A,A)') '>>> MISSING RESTART COMPLETION MARKER: ', done_files(irestart)%to_char()
                    THROW_HARD('abinitio3D_cavg_sort: one or more restarts did not finish')
                endif
            enddo
        end subroutine submit_restarts

        subroutine read_restart_labels
            do irestart = 1, params%nrestarts
                call restart_proj%kill
                call restart_proj%read_segment('cls3D', restart_projfiles(irestart))
                if( restart_proj%os_cls3D%get_noris() /= ncls )then
                    write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> RESTART ', irestart, ' #CLS3D: ', &
                        restart_proj%os_cls3D%get_noris(), ' EXPECTED: ', ncls
                    THROW_HARD('abinitio3D_cavg_sort: restart cls3D count mismatch')
                endif
                restart_labels(irestart,:) = restart_proj%os_cls3D%get_all_asint('state')
                call sanitize_restart_labels(irestart)
            enddo
            call restart_proj%kill
        end subroutine read_restart_labels

        subroutine sanitize_restart_labels( restart_ind )
            integer, intent(in) :: restart_ind
            integer :: icls
            do icls = 1, ncls
                if( original_states(icls) <= 0 )then
                    restart_labels(restart_ind,icls) = 0
                else if( restart_labels(restart_ind,icls) < 1 .or. restart_labels(restart_ind,icls) > SORT_NSTATES )then
                    restart_labels(restart_ind,icls) = 0
                endif
            enddo
        end subroutine sanitize_restart_labels

        subroutine map_state_correspondence
            integer :: same, swapped, icls
            logical :: flip
            mapped_labels(1,:) = restart_labels(1,:)
            do irestart = 2, params%nrestarts
                same    = 0
                swapped = 0
                do icls = 1, ncls
                    if( restart_labels(1,icls) < 1 .or. restart_labels(irestart,icls) < 1 ) cycle
                    if( restart_labels(1,icls) == restart_labels(irestart,icls) ) same = same + 1
                    if( restart_labels(1,icls) == swapped_label(restart_labels(irestart,icls)) ) swapped = swapped + 1
                enddo
                flip = swapped > same
                do icls = 1, ncls
                    if( restart_labels(irestart,icls) < 1 )then
                        mapped_labels(irestart,icls) = 0
                    else if( flip )then
                        mapped_labels(irestart,icls) = swapped_label(restart_labels(irestart,icls))
                    else
                        mapped_labels(irestart,icls) = restart_labels(irestart,icls)
                    endif
                enddo
                write(logfhandle,'(A,I0,A,I0,A,I0,A,L1)') '>>> RESTART ', irestart, &
                    ' STATE-LABEL AGREEMENT SAME/SWAPPED: ', same, ' / ', swapped, ' FLIPPED=', flip
            enddo
        end subroutine map_state_correspondence

        integer function swapped_label( label )
            integer, intent(in) :: label
            if( label == 1 )then
                swapped_label = 2
            else if( label == 2 )then
                swapped_label = 1
            else
                swapped_label = 0
            endif
        end function swapped_label

        subroutine build_consensus
            integer :: icls, label, best_label, best_votes
            votes = 0
            do icls = 1, ncls
                if( original_states(icls) <= 0 ) cycle
                do irestart = 1, params%nrestarts
                    label = mapped_labels(irestart,icls)
                    if( label >= 1 .and. label <= SORT_NSTATES ) votes(label,icls) = votes(label,icls) + 1
                enddo
                best_label = 0
                best_votes = -1
                do label = 1, SORT_NSTATES
                    if( votes(label,icls) > best_votes )then
                        best_label = label
                        best_votes = votes(label,icls)
                    endif
                enddo
                if( votes(1,icls) == votes(2,icls) .and. mapped_labels(1,icls) > 0 ) best_label = mapped_labels(1,icls)
                consensus(icls) = best_label
            enddo
        end subroutine build_consensus

        subroutine assign_good_bad_states
            integer :: good_consensus, bad_consensus, icls, pop1, pop2
            real    :: mean1, mean2
            mean1 = mean_quality_for_consensus(1, pop1)
            mean2 = mean_quality_for_consensus(2, pop2)
            if( pop1 == 0 .and. pop2 == 0 ) THROW_HARD('abinitio3D_cavg_sort: no active consensus classes')
            if( pop2 == 0 .or. (pop1 > 0 .and. mean1 >= mean2) )then
                good_consensus = 1
                bad_consensus  = 2
            else
                good_consensus = 2
                bad_consensus  = 1
            endif
            do icls = 1, ncls
                select case(consensus(icls))
                    case(1,2)
                        if( consensus(icls) == good_consensus )then
                            final_states(icls) = 1
                        else if( consensus(icls) == bad_consensus )then
                            final_states(icls) = 2
                        endif
                    case default
                        final_states(icls) = 0
                end select
            enddo
            write(logfhandle,'(A,I0,A,F8.3,A,I0)') '>>> CONSENSUS STATE 1 QUALITY MEAN / POP: ', 1, ' ', mean1, ' / ', pop1
            write(logfhandle,'(A,I0,A,F8.3,A,I0)') '>>> CONSENSUS STATE 2 QUALITY MEAN / POP: ', 2, ' ', mean2, ' / ', pop2
            write(logfhandle,'(A,I0,A,I0)') '>>> CAVG SORT GOOD/BAD CONSENSUS STATES: ', good_consensus, ' / ', bad_consensus
            write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> CAVG SORT FINAL GOOD/BAD/INACTIVE: ', &
                count(final_states == 1), ' / ', count(final_states == 2), ' / ', count(final_states == 0)
        end subroutine assign_good_bad_states

        real function mean_quality_for_consensus( label, pop )
            integer, intent(in)  :: label
            integer, intent(out) :: pop
            integer :: icls
            mean_quality_for_consensus = 0.0
            pop = 0
            do icls = 1, ncls
                if( consensus(icls) /= label ) cycle
                if( original_states(icls) <= 0 ) cycle
                pop = pop + 1
                mean_quality_for_consensus = mean_quality_for_consensus + quality_scores(icls)
            enddo
            if( pop > 0 ) mean_quality_for_consensus = mean_quality_for_consensus / real(pop)
        end function mean_quality_for_consensus

        subroutine annotate_project
            integer :: icls, ncls3d, max_vote
            ncls3d = spproj%os_cls3D%get_noris()
            do icls = 1, ncls
                max_vote = max(votes(1,icls), votes(2,icls))
                call spproj%os_cls2D%set(icls, 'quality',             quality_scores(icls))
                call spproj%os_cls2D%set(icls, 'accept',              merge(1, 0, final_states(icls) == 1))
                call spproj%os_cls2D%set(icls, 'quality_cluster',     quality%labels(icls))
                call spproj%os_cls2D%set(icls, 'cavg_sort_consensus', consensus(icls))
                call spproj%os_cls2D%set(icls, 'cavg_sort_votes',     max_vote)
                if( ncls3d == ncls )then
                    call spproj%os_cls3D%set(icls, 'quality',             quality_scores(icls))
                    call spproj%os_cls3D%set(icls, 'accept',              merge(1, 0, final_states(icls) == 1))
                    call spproj%os_cls3D%set(icls, 'quality_cluster',     quality%labels(icls))
                    call spproj%os_cls3D%set(icls, 'cavg_sort_consensus', consensus(icls))
                    call spproj%os_cls3D%set(icls, 'cavg_sort_votes',     max_vote)
                endif
            enddo
        end subroutine annotate_project

        subroutine write_consensus_report
            integer :: funit, icls
            open(newunit=funit, file=CONSENSUS_REPORT, status='replace', action='write')
            write(funit,'(A)') '# abinitio3D_cavg_sort consensus report'
            write(funit,'(A,I0)') '# nrestarts=', params%nrestarts
            write(funit,'(A,I0)') '# nstages=', SORT_NSTAGES
            write(funit,'(A,A)') '# quality_model=', trim(model%name)
            write(funit,'(A)', advance='no') 'class,original_state,consensus_state,final_state,votes_state1,votes_state2,quality_score,quality_auto_state'
            do irestart = 1, params%nrestarts
                write(funit,'(A,I0)', advance='no') ',restart_raw_', irestart
            enddo
            do irestart = 1, params%nrestarts
                write(funit,'(A,I0)', advance='no') ',restart_mapped_', irestart
            enddo
            write(funit,*)
            do icls = 1, ncls
                write(funit,'(I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,ES14.6,A,I0)', advance='no') &
                    icls, ',', original_states(icls), ',', consensus(icls), ',', final_states(icls), ',', &
                    votes(1,icls), ',', votes(2,icls), ',', quality_scores(icls), ',', quality_auto_states(icls)
                do irestart = 1, params%nrestarts
                    write(funit,'(A,I0)', advance='no') ',', restart_labels(irestart,icls)
                enddo
                do irestart = 1, params%nrestarts
                    write(funit,'(A,I0)', advance='no') ',', mapped_labels(irestart,icls)
                enddo
                write(funit,*)
            enddo
            close(funit)
            write(logfhandle,'(A,A)') '>>> WROTE ', CONSENSUS_REPORT
            call write_cavg_quality_feature_table(quality, model, 'abinitio3D_cavg_sort_quality_features.txt', &
                params%projfile%to_char())
        end subroutine write_consensus_report

        subroutine cleanup
            call spproj%kill
            call restart_proj%kill
            call quality%kill
            call dealloc_imgarr(cavg_imgs)
            if( allocated(restart_labels)      ) deallocate(restart_labels)
            if( allocated(mapped_labels)       ) deallocate(mapped_labels)
            if( allocated(consensus)           ) deallocate(consensus)
            if( allocated(final_states)        ) deallocate(final_states)
            if( allocated(original_states)     ) deallocate(original_states)
            if( allocated(quality_auto_states) ) deallocate(quality_auto_states)
            if( allocated(quality_scores)      ) deallocate(quality_scores)
            if( allocated(votes)               ) deallocate(votes)
            if( allocated(restart_projfiles)   ) deallocate(restart_projfiles)
            if( allocated(done_files)          ) deallocate(done_files)
        end subroutine cleanup

    end subroutine exec_abinitio3D_cavg_sort

    !> for generation of an initial 3d model from particles
    subroutine exec_abinitio3D( self, cline )
        class(commander_abinitio3D), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        ! commanders
        type(commander_refine3D)        :: xrefine3D
        type(commander_rec3D)           :: xrec3D
        ! other
        real,               allocatable :: rstates(:)
        integer,            allocatable :: tmpinds(:), clsinds(:), pinds(:)
        type(class_sample), allocatable :: clssmp(:)
        type(parameters)                :: params
        type(sp_project)                :: spproj
        type(simple_nice_comm)          :: nice_comm
        real    :: lprange(2)
        integer :: state, istage, icls, start_stage, nptcls2update, noris, nstates_on_cline
        integer :: nstates_in_project, split_stage, last_stage
        logical :: l_cavg_ini_ext, l_vol_ini_ext
        call cline%set('objfun',    'euclid') ! use noise normalized Euclidean distances from the start
        call cline%set('sigma_est', 'global') ! obviously
        call cline%set('bfac',            0.) ! because initial models should not be sharpened
        call cline%set('nu_refine',     'no') ! no nonuniform refinement
        if( .not. cline%defined('mkdir')               ) call cline%set('mkdir',                    'yes')
        if( .not. cline%defined('overlap')             ) call cline%set('overlap',                   0.95)
        if( .not. cline%defined('prob_athres')         ) call cline%set('prob_athres',                10.)
        if( .not. cline%defined('center')              ) call cline%set('center',                    'no')
        if( .not. cline%defined('cenlp')               ) call cline%set('cenlp', abinitio_cenlp_default())
        if( .not. cline%defined('oritype')             ) call cline%set('oritype',               'ptcl3D')
        if( .not. cline%defined('pgrp')                ) call cline%set('pgrp',                      'c1')
        if( .not. cline%defined('pgrp_start')          ) call cline%set('pgrp_start',                'c1')
        if( .not. cline%defined('filt_mode')           ) call cline%set('filt_mode',         'nonuniform')
        if( .not. cline%defined('automsk')             ) call cline%set('automsk',                   'no')
        if( .not. cline%defined('gauref')              ) call cline%set('gauref',                   'yes')
        if( .not. cline%defined('nsample_start')       ) call cline%set('nsample_start', abinitio_nsample_start_default())
        ! splitting stage
        split_stage = abinitio_het_docked_stage()
        if( cline%defined('split_stage') ) split_stage = cline%get_iarg('split_stage')
        ! adjust default multivol_mode unless given on command line
        if( cline%defined('nstates') )then
            nstates_on_cline = cline%get_iarg('nstates')
            if( nstates_on_cline > 1 .and. .not. cline%defined('multivol_mode') )then
                call cline%set('multivol_mode', 'independent')
            endif
        endif
        ! make master parameters
        call params%new(cline)
        select case(trim(params%filt_mode))
            case('uniform','fsc')
                THROW_HARD('abinitio3D no longer supports automatic low-pass filt_mode=uniform|fsc; &
                    &use none|nonuniform')
        end select
        call cline%set('mkdir', 'no')
        call cline%delete('algorithm')
        ! optional early stop stage, matching the abinitio3D_cavgs nstages policy
        last_stage = abinitio_nstages()
        if( cline%defined('nstages') )then
            if( params%nstages < 1 ) THROW_HARD('nstages must be >= 1 for abinitio3D')
            last_stage = min(abinitio_nstages(), params%nstages)
        endif
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
        ! nice communicator init
        call nice_comm%init(params%niceprocid, params%niceserver)
        call nice_comm%cycle()
        ! read project
        call spproj%read(params%projfile)
        ! provide initialization of 3D alignment using class averages?
        start_stage = 1
        l_ini3D     = .false.
        if( trim(params%cavg_ini).eq.'yes' )then
            if( str_has_substr(params%multivol_mode,'input_oris') ) THROW_HARD('Ini3D on cavgs not allowed for multivol_mode=input_oris*')
            if( last_stage < abinitio_nstages_ini3D() - 1 ) THROW_HARD('nstages must be >= first executable abinitio3D stage')
            ! nice
            nice_comm%stat_root%stage = "initialising 3D volume from class averages"
            call nice_comm%cycle()
            ! execution
            call ini3D_from_cavgs(cline)
            ! re-read the project file to update info in spproj
            call spproj%read(params%projfile)
            start_stage = abinitio_nstages_ini3D() - 1 ! compute reduced to two overlapping stages
            l_ini3D     = .true.
            ! symmetry dealt with by ini3D
        else
            if( trim(params%multivol_mode).eq.'independent' )then
                ! turn off symmetry axis search and put the symmetry in from the start
                params%pgrp_start = params%pgrp
            endif
        endif
        ! nice
        nice_comm%stat_root%stage = "preparing workflow"
        call nice_comm%cycle()
        ! initialization on class averages done outside this workflow (externally)?
        l_cavg_ini_ext = trim(params%cavg_ini_ext).eq.'yes'
        if( l_cavg_ini_ext )then
            if( last_stage < abinitio_symsrch_stage() + 1 ) THROW_HARD('nstages must be >= first executable abinitio3D stage')
            ! check that ptcl3D field is not virgin
            if( spproj%is_virgin_field('ptcl3D') )then
                THROW_HARD('Prior 3D alignment required for abinitio workflow when cavg_ini_ext is set to yes')
            endif
            call validate_cavg_ini_ext_states
            ! symmetry axis search is skipped: input orientations are assumed already symmetrized
            call cline%set('pgrp_start', params%pgrp)
            params%pgrp_start = params%pgrp
            start_stage = abinitio_symsrch_stage() + 1 ! start after the symmetry search stage
            l_ini3D     = .true.
        endif
        ! initialization of input volumes originating from outside the workflow
        l_vol_ini_ext = cline%defined('vol1')
        if( l_vol_ini_ext )then
            ! sanity checks, it is also assumed no 2D clustering info has been performed
            ! resolution limits have to be defined
            select case(trim(params%multivol_mode))
            case('single','independent','docked')
                ! volume input only allowed for these modes
                if( (params%nstates > 1)  )then
                    ! making sure all volumes are present (for 'docked', nstates==1 here)
                    do state = 2, params%nstates
                        if( .not. cline%defined('vol'//int2str(state)) )then
                            THROW_HARD('vol'//int2str(state)//' must be defined for state s='//int2str(state))
                        endif
                    enddo
                endif
            case DEFAULT
                THROW_HARD('Unsupported volume input and multivol_mode: '//trim(params%multivol_mode))
            end select
            if( l_ini3D ) THROW_HARD('Cannot have both class initialization and an input volume')
            if( trim(params%partition).eq.'yes' ) THROW_HARD('Volume input not currently supported with partition=yes')
            ! setting up random classes for particles sampling
            call spproj%os_ptcl2D%rnd_cls(100)
            call spproj%write_segment_inside('ptcl2D', params%projfile)
            call spproj%os_cls2D%new(100, is_ptcl=.false.)
            call spproj%os_cls2D%set_all2single('state', 1)
        endif
        ! set class global filtering flags for staged refine3D policy
        l_nonuniform = params%l_nonuniform
        nstages_refine3D = last_stage
        if( str_has_substr(params%multivol_mode,'input_oris') ) start_stage = split_stage
        if( nstages_refine3D < start_stage )then
            THROW_HARD('nstages must be >= first executable abinitio3D stage')
        endif
        ! set class global automasking flag (now supported for all multivol modes via state-specific masks)
        l_automsk = (cline%defined('automsk') .and. trim(params%automsk).ne.'no')
        ! prepare class command lines
        call prep_class_command_lines(params, cline, params%projfile)
        ! set symmetry class variables
        call set_symmetry_class_vars(params)
        ! fall over if there are no particles
        if( spproj%os_ptcl3D%get_noris() < 1 ) THROW_HARD('Particles could not be found in the project')
        ! take care of class-biased particle sampling
        if( spproj%is_virgin_field('ptcl2D') )then
            THROW_HARD('Prior 2D clustering required for abinitio workflow')
        else
            update_frac = 1.0
            nptcls_eff  = spproj%count_state_gt_zero()
            if( .not. cline%defined('nsample') ) params%nsample = abinitio_nsample_default()
            if( params%nsample < 1 ) THROW_HARD('nsample must be >= 1 for abinitio3D sampled update')
            if( cline%defined('nsample_start') )then
                if( params%nsample_start < 1 ) THROW_HARD('nsample_start must be >= 1 for abinitio3D sampled update')
                if( params%nsample_start > params%nsample )then
                    THROW_HARD('nsample_start must be <= nsample for abinitio3D sampled update ramp')
                endif
                write(logfhandle,'(A,I0,A,I0,A,I0)') &
                    &'>>> ABINITIO3D NSAMPLE RAMP: ', params%nsample_start, ' -> ', params%nsample, &
                    &' BY STAGE ', abinitio_symsrch_stage() + 2
            endif
            update_frac = real(params%nsample * params%nstates) / real(nptcls_eff)
            update_frac = min(abinitio_update_frac_max(), update_frac) ! to ensure fractional update is always on
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
        if( l_vol_ini_ext )then
            ! limits based on dimensions or input
            call mskdiam2lplimits( params%mskdiam, lprange(1), lprange(2), params%cenlp )
            if( .not.cline%defined('lpstart') ) params%lpstart = lprange(1)
            if( .not.cline%defined('lpstop')  )then
                params%lpstop = lprange(2)
                lprange       = abinitio_lpstop_bounds()
                params%lpstop = min(params%lpstop, lprange(1))
            endif
            call set_lplims_from_input(params, spproj, params%lpstart, params%lpstop)
        else
            if( cline%defined('lpstart') .and. cline%defined('lpstop') )then
                call set_lplims_from_frcs(params, spproj, l_cavgs=.false., lpstart=params%lpstart, lpstop=params%lpstop)
            else if( cline%defined('lpstart') )then
                call set_lplims_from_frcs(params, spproj, l_cavgs=.false., lpstart=params%lpstart)
            else if( cline%defined('lpstop') )then
                call set_lplims_from_frcs(params, spproj, l_cavgs=.false., lpstop=params%lpstop)
            else
                call set_lplims_from_frcs(params, spproj, l_cavgs=.false.)
            endif
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
            call calc_rec(params, params%projfile, xrec3D, start_stage)
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
            if( trim(params%multivol_mode).eq.'independent' .and. .not.l_cavg_ini_ext )then
                call gen_labelling(spproj%os_ptcl3D, params%nstates, 'squared_uniform')
            endif
            call spproj%write_segment_inside(params%oritype, params%projfile)
            if( l_vol_ini_ext )then
                ! user provided input volumes
                call normalize_input_volumes(params, cline_refine3D)
            else
                ! create noise starting volume(s)
                call generate_random_volumes(params, lpinfo(1)%box_crop, lpinfo(1)%smpd_crop, cline_refine3D)
            endif
        else
            ! check that ptcl3D field is not virgin
            if( spproj%is_virgin_field('ptcl3D') )then
                THROW_HARD('Prior 3D alignment is lacking for starting volume generation')
            endif
            ! randomize states
            if( trim(params%multivol_mode).eq.'independent' .and. .not.l_cavg_ini_ext )then
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
            call calc_rec(params, params%projfile, xrec3D, start_stage)
        endif
        if( cline%defined('nstages') )then
            write(logfhandle,'(A,I0,A,I0)')'>>> ABINITIO3D STAGE RANGE: ', start_stage, ' -> ', nstages_refine3D
            if( nstages_refine3D < abinitio_nstages() )then
                write(logfhandle,'(A)')'>>> ABINITIO3D EARLY STAGE STOP: SKIPPING FINAL ALL-PARTICLE RECONSTRUCTION'
            endif
        endif
        ! Frequency marching
        call print_states(params, 0)
        ! nice
        nice_comm%stat_root%stage = "starting workflow"
        call nice_comm%cycle()
        do istage = start_stage, nstages_refine3D
            ! nice
             if( nice_comm%stop )then
                ! termination
                write(logfhandle,'(A)')'>>> USER COMMANDED STOP'
                call spproj%kill
                call qsys_cleanup(params)
                call nice_comm%terminate(stop=.true.)
                call simple_end('**** SIMPLE_ABINITIO3D USER STOP ****')
                call EXIT(0)
            endif
            nice_comm%stat_root%stage = "running workflow"
            call nice_comm%update_ini3D(stage=istage, number_states=nstates_glob, lp=lpinfo(istage)%lp) 
            call nice_comm%cycle()
            ! At the splitting stage of docked mode: reset the nstates in params
            if( params%multivol_mode.eq.'docked' .and. istage == split_stage )then
                params%nstates = nstates_glob
                update_frac    = min(update_frac * nstates_glob, abinitio_update_frac_max())
            endif
            ! Preparation of command line for refinement
            call set_cline_refine3D(params, istage, l_cavgs=.false.)
            write(logfhandle,'(A)')'>>>'
            if( cline_refine3D%defined('lp') )then
                write(logfhandle,'(A,I3,A9,F5.1)')'>>> STAGE ', istage,' WITH LP =', cline_refine3D%get_rarg('lp')
            else
                write(logfhandle,'(A,I3,A)')'>>> STAGE ', istage,' WITH NU-SELECTED MATCHING LP'
            endif
            ! Need to be here since rec cline depends on refine3D cline
            if( params%multivol_mode.eq.'docked' .and. istage == split_stage )then
                call randomize_states(params, spproj, params%projfile, xrec3D, split_stage)
            endif
            if( lpinfo(istage)%l_autoscale )then
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',lpinfo(istage)%box_crop
            endif
            ! Executing the refinement with the above settings
            write(logfhandle,'(A,I0)')'>>> ABINITIO3D ENTERING REFINE3D STAGE ', istage
            call flush(logfhandle)
            call exec_refine3D(params, istage, xrefine3D)
            write(logfhandle,'(A,I0)')'>>> ABINITIO3D RETURNED FROM REFINE3D STAGE ', istage
            call flush(logfhandle)
            call print_states(params, istage)
            ! Symmetrization
            if( istage == abinitio_symsrch_stage() )then
                call symmetrize(params, istage, spproj, params%projfile, xrec3D)
            endif
            ! nice
            call nice_comm%update_ini3D(last_stage_completed=.true.) 
            call nice_comm%cycle()
        enddo
        if( nstages_refine3D == abinitio_nstages() )then
            ! calculate 3D reconstruction at original sampling
            call calc_final_rec(params, spproj, params%projfile, xrec3D, l_postprocess=.true.)
            ! for visualization
            call gen_ortho_reprojs4viz(params, spproj)
            ! final raw and low-pass diagnostic 3D reconstruction outputs
            call write_final_rec_outputs(params, spproj, lpinfo(nstages_refine3D)%lp)
        else
            write(logfhandle,'(A,I0)')'>>> ABINITIO3D EARLY STOP AFTER STAGE ', nstages_refine3D
            write(logfhandle,'(A)')'>>> FINAL ALL-PARTICLE RECONSTRUCTION SKIPPED'
        endif
        ! termination
        nice_comm%stat_root%stage = "terminating"
        call nice_comm%cycle()
        ! cleanup
        call nice_comm%terminate(export_project=spproj)
        call spproj%kill
        call qsys_cleanup(params)
        call simple_end('**** SIMPLE_ABINITIO3D NORMAL STOP ****')

    contains

        subroutine ini3D_from_cavgs( cline )
            class(cmdline),    intent(inout) :: cline
            type(commander_abinitio3D_cavgs) :: xini3D
            type(cmdline)                    :: cline_ini3D
            type(string),    allocatable     :: files_that_stay(:)
            character(len=*), parameter      :: INI3D_DIR='abinitio3D_cavgs/'
            cline_ini3D = cline
            call cline_ini3D%set('nstages', abinitio_nstages_ini3D())
            ! Resolution limits
            if( .not. cline_ini3D%defined('lpstart_ini3D') ) call cline_ini3D%set('lpstart_ini3D', abinitio_lpstart_ini3D())
            if( .not. cline_ini3D%defined('lpstop_ini3D')  ) call cline_ini3D%set('lpstop_ini3D',  abinitio_lpstop_ini3D())
            if( cline%defined('lpstart_ini3D') )then
                call cline_ini3D%set('lpstart', params%lpstart_ini3D)
                call cline_ini3D%delete('lpstart_ini3D')
            endif
            if( cline%defined('lpstop_ini3D') )then
                call cline_ini3D%set('lpstop', params%lpstop_ini3D)
                call cline_ini3D%delete('lpstop_ini3D')
            endif
            ! Compute
            if( cline%defined('nthr_ini3D') )then
                call cline_ini3D%set('nthr', params%nthr_ini3D)
                call cline_ini3D%delete('nthr_ini3D')
            endif
            call cline_ini3D%delete('nstates') ! cavg_ini under the assumption of one state
            call cline_ini3D%delete('oritype')
            call cline_ini3D%delete('imgkind')
            call cline_ini3D%delete('prob_athres')
            call xini3D%execute(cline_ini3D)
            ! update point-group symmetry
            call cline%set('pgrp_start', params%pgrp)
            params%pgrp_start = params%pgrp
            call prep_class_command_lines(params, cline, params%projfile)
            ! stash away files
            ! identfy files that stay
            allocate(files_that_stay(7))
            files_that_stay(1) = basename(params%projfile)
            files_that_stay(2) = 'cavgs'
            files_that_stay(3) = 'nice'
            files_that_stay(4) = 'frcs'
            files_that_stay(5) = 'ABINITIO3D'
            files_that_stay(6) = 'execscript' ! only with streaming
            files_that_stay(7) = 'execlog'    ! only with streaming
            ! make the move
            call move_files_in_cwd(string(INI3D_DIR), files_that_stay)
        end subroutine ini3D_from_cavgs

        subroutine validate_cavg_ini_ext_states
            integer :: state, pop
            if( params%nstates <= 1 ) return
            nstates_in_project = spproj%os_ptcl3D%get_n('state')
            if( nstates_in_project /= params%nstates )then
                write(logfhandle,*) 'requested nstates, project ptcl3D state bins: ', params%nstates, nstates_in_project
                THROW_HARD('cavg_ini_ext=yes with nstates>1 requires matching existing ptcl3D state assignments')
            endif
            do state = 1,params%nstates
                pop = spproj%os_ptcl3D%get_pop(state, 'state')
                if( pop < 1 )then
                    write(logfhandle,*) 'empty ptcl3D state for cavg_ini_ext: ', state
                    THROW_HARD('cavg_ini_ext=yes requires every requested state to be populated')
                endif
            enddo
        end subroutine validate_cavg_ini_ext_states

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
