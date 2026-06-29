!@descr: abinitio 3D reconstruction in single- and multi-particle mode
module simple_commanders_abinitio
use simple_commanders_api
use simple_abinitio_utils
use simple_qsys_funs,               only: qsys_watcher_diag
use simple_procimgstk,              only: shift_imgfile
use simple_commanders_project_core, only: commander_selection
use simple_commanders_reproject,    only: commander_reproject
use simple_commanders_refine3D,     only: commander_refine3D, commander_refine3D
use simple_commanders_rec,          only: commander_rec3D, commander_rec3D
use simple_cluster_seed,            only: gen_labelling
use simple_refine3D_fnames,         only: refine3D_startvol_fname, refine3D_startvol_half_fname, &
    &refine3D_state_vol_fname, refine3D_state_halfvol_fname
implicit none

public :: commander_abinitio3D_cavgs, commander_abinitio3D_cavgs_reject, commander_abinitio3D
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_abinitio3D_cavgs
    contains
    procedure :: execute => exec_abinitio3D_cavgs
end type commander_abinitio3D_cavgs

type, extends(commander_base) :: commander_abinitio3D_cavgs_reject
    contains
    procedure :: execute => exec_abinitio3D_cavgs_reject
end type commander_abinitio3D_cavgs_reject

type, extends(commander_base) :: commander_abinitio3D
    contains
    procedure :: execute => exec_abinitio3D
end type commander_abinitio3D

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
        l_state_continue_mode = .false.
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

    !> reject class averages by consensus over restarted multi-state abinitio3D_cavgs runs
    subroutine exec_abinitio3D_cavgs_reject( self, cline )
        use simple_cavg_quality_analysis, only: evaluate_cavg_quality, write_cavg_quality_feature_table
        use simple_cavg_quality_model,    only: CAVG_QUALITY_MODEL_CHUNK_DEFAULT, cavg_quality_model
        use simple_cavg_quality_types,    only: cavg_quality_result
        use simple_imgarr_utils,          only: read_cavgs_into_imgarr, dealloc_imgarr
        class(commander_abinitio3D_cavgs_reject), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline
        character(len=*), parameter :: RESTART_DIR_FBODY    = 'abinitio3D_cavgs_reject_restart_'
        character(len=*), parameter :: RESTART_DONE         = 'ABINITIO3D_CAVGS_REJECT_FINISHED'
        character(len=*), parameter :: DOCK_DIR_FBODY       = 'abinitio3D_cavgs_reject_dock_'
        character(len=*), parameter :: DOCK_DONE            = 'ABINITIO3D_CAVGS_REJECT_DOCK_FINISHED'
        character(len=*), parameter :: CONSENSUS_VOL        = 'abinitio3D_cavgs_reject_consensus_vol.mrc'
        character(len=*), parameter :: CONSENSUS_REPORT     = 'abinitio3D_cavgs_reject_consensus.txt'
        character(len=*), parameter :: CONSENSUS_VOL_REPORT = 'abinitio3D_cavgs_reject_consensus_volume.txt'
        integer,          parameter :: DEFAULT_SORT_NSTATES = 2
        integer,          parameter :: MAX_SORT_NSTATES     = 3
        integer,          parameter :: SORT_NSTAGES         = 3
        real,             parameter :: DEFAULT_DOCK_HP      = 100.0
        real,             parameter :: DEFAULT_DOCK_LP      = 10.0
        type(parameters)          :: params
        type(qsys_env)            :: qenv
        type(sp_project)          :: spproj, restart_proj
        type(cavg_quality_model)  :: model
        type(cavg_quality_result) :: quality
        type(image), allocatable  :: cavg_imgs(:)
        type(string)              :: cwd, cwd_run, projbase, folder, restart_projfile, done_file
        type(string), allocatable :: restart_projfiles(:), done_files(:), selected_vols(:), docked_vols(:), dock_reports(:)
        integer, allocatable      :: restart_labels(:,:), mapped_labels(:,:), state_maps(:,:)
        integer, allocatable      :: consensus(:), final_states(:), selection_states(:)
        integer, allocatable      :: original_states(:), quality_auto_states(:), votes(:,:)
        integer, allocatable      :: restart_good_state(:), restart_good_pop(:)
        real,    allocatable      :: quality_scores(:), restart_good_scores(:)
        integer                   :: ncls, irestart, sort_nstates, good_consensus, ref_restart, consensus_ldim(3)
        real                      :: consensus_smpd, dock_lp
        if( cline%defined('part') )then
            THROW_HARD('abinitio3D_cavgs_reject is master-only; remove part from command line')
        endif
        if( cline%defined('nstates') )then
            if( cline%get_iarg('nstates') < 2 .or. cline%get_iarg('nstates') > MAX_SORT_NSTATES )then
                THROW_HARD('abinitio3D_cavgs_reject supports nstates=2 or nstates=3')
            endif
        else
            call cline%set('nstates', DEFAULT_SORT_NSTATES)
        endif
        if( cline%defined('nstages') .and. cline%get_iarg('nstages') /= SORT_NSTAGES )then
            THROW_HARD('abinitio3D_cavgs_reject always runs abinitio3D_cavgs through nstages='//int2str(SORT_NSTAGES))
        endif
        call cline%set('nstages',  SORT_NSTAGES)
        call cline%set('pgrp',     'c1')
        call cline%delete('pgrp_start')
        if( .not.cline%defined('lpstart')       ) call cline%set('lpstart', abinitio_lpstart_ini3D())
        if( .not.cline%defined('lpstop')        ) call cline%set('lpstop',   abinitio_lpstop_ini3D())
        if( .not.cline%defined('nrestarts')     ) call cline%set('nrestarts', 3)
        if( .not.cline%defined('mkdir')         ) call cline%set('mkdir', 'yes')
        if( .not.cline%defined('quality_model') ) call cline%set('quality_model', CAVG_QUALITY_MODEL_CHUNK_DEFAULT)
        if( .not.cline%defined('prune')         ) call cline%set('prune', 'no')
        call params%new(cline)
        if( params%nrestarts < 1 ) THROW_HARD('abinitio3D_cavgs_reject requires nrestarts >= 1')
        sort_nstates  = params%nstates
        good_consensus = 0
        ref_restart    = 0
        consensus_ldim = 0
        consensus_smpd = 0.0
        dock_lp = DEFAULT_DOCK_LP
        call spproj%read(params%projfile)
        ncls = spproj%os_cls2D%get_noris()
        if( ncls == 0 ) THROW_HARD('abinitio3D_cavgs_reject: project has no cls2D entries')
        original_states = spproj%os_cls2D%get_all_asint('state')
        if( size(original_states) /= ncls ) THROW_HARD('abinitio3D_cavgs_reject: invalid cls2D state array')
        allocate(restart_labels(params%nrestarts,ncls),     source=0)
        allocate(mapped_labels(params%nrestarts,ncls),      source=0)
        allocate(state_maps(params%nrestarts,sort_nstates), source=0)
        allocate(votes(sort_nstates,ncls),                  source=0)
        allocate(consensus(ncls),                           source=0)
        allocate(final_states(ncls),                        source=0)
        allocate(selection_states(ncls),                    source=0)
        allocate(restart_good_state(params%nrestarts),      source=0)
        allocate(restart_good_pop(params%nrestarts),        source=0)
        allocate(restart_good_scores(params%nrestarts),     source=0.0)
        allocate(restart_projfiles(params%nrestarts))
        allocate(done_files(params%nrestarts))
        allocate(selected_vols(params%nrestarts))
        allocate(docked_vols(params%nrestarts))
        allocate(dock_reports(params%nrestarts))
        call init_quality_model
        call evaluate_quality
        call submit_restarts
        call read_restart_labels
        call map_state_correspondence
        call build_consensus
        call assign_good_bad_states
        call build_selection_states
        call write_selection_outputs
        call spproj%map_cavgs_selection(selection_states)
        call annotate_project
        call build_consensus_volume
        if( trim(params%prune) == 'yes' ) call spproj%prune_particles
        call spproj%write(params%projfile)
        call write_consensus_report
        call write_consensus_volume_report
        call cleanup
        call simple_end('**** SIMPLE_ABINITIO3D_CAVGS_REJECT NORMAL STOP ****', &
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
            if( size(cavg_imgs) /= ncls ) THROW_HARD('abinitio3D_cavgs_reject: # cavgs /= # cls2D entries')
            call evaluate_cavg_quality(cavg_imgs, spproj%os_cls2D, params%mskdiam, quality, model)
            quality_scores      = quality%scores
            quality_auto_states = quality%states
            write(logfhandle,'(A,A)') '>>> CAVGS REJECT QUALITY MODEL         : ', trim(model%name)
            write(logfhandle,'(A,A)') '>>> CAVGS REJECT QUALITY MODEL CONTEXT : ', trim(model%context)
            write(logfhandle,'(A,F8.3,A,F8.3,A,F8.3)') '>>> CAVGS REJECT QUALITY THRESHOLD RAW / OFFSET / EFFECTIVE : ', &
                quality%raw_threshold, ' / ', quality%threshold_offset, ' / ', quality%threshold
            write(logfhandle,'(A,F8.3,A,L1)') '>>> CAVGS REJECT QUALITY SEPARATION / USED THRESHOLD : ', &
                quality%separation, ' USED=', quality%used_threshold
        end subroutine evaluate_quality

        subroutine submit_restarts
            type(cmdline) :: cline_restart
            call simple_getcwd(cwd)
            CWD_GLOB = cwd%to_char()
            projbase = basename(params%projfile)
            call qenv%new(params, 1, exec_bin=string('simple_exec'))
            do irestart = 1, params%nrestarts
                folder           = RESTART_DIR_FBODY//int2str_pad(irestart, 3)
                restart_projfile = filepath(folder, projbase)
                done_file        = filepath(folder, string(RESTART_DONE))
                call simple_mkdir(folder)
                call simple_copy_file(params%projfile, restart_projfile)
                if( file_exists(done_file) ) call del_file(done_file)
                restart_projfiles(irestart) = filepath(cwd, folder, projbase)
                done_files(irestart)        = filepath(cwd, folder, string(RESTART_DONE))
                cline_restart = cline
                call cline_restart%set('prg',                'abinitio3D_cavgs')
                call cline_restart%set('projfile',           projbase)
                call cline_restart%set('mkdir',              'no')
                call cline_restart%set('nstates',            sort_nstates)
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
                    string('abinitio3D_cavgs_reject_script'), string('abinitio3D_cavgs_reject.log'))
                call simple_chdir(cwd)
                CWD_GLOB = cwd%to_char()
                call cline_restart%kill
                write(logfhandle,'(A,I0,A,A)') '>>> SUBMITTED ABINITIO3D_CAVGS RESTART ', irestart, &
                    ' IN ', folder%to_char()
            enddo
            call qsys_watcher_diag(done_files)
            do irestart = 1, params%nrestarts
                if( .not.file_exists(done_files(irestart)) )then
                    write(logfhandle,'(A,A)') '>>> MISSING RESTART COMPLETION MARKER: ', done_files(irestart)%to_char()
                    THROW_HARD('abinitio3D_cavgs_reject: one or more restarts did not finish')
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
                    THROW_HARD('abinitio3D_cavgs_reject: restart cls3D count mismatch')
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
                else if( restart_labels(restart_ind,icls) < 1 .or. restart_labels(restart_ind,icls) > sort_nstates )then
                    restart_labels(restart_ind,icls) = 0
                endif
            enddo
        end subroutine sanitize_restart_labels

        subroutine map_state_correspondence
            integer :: best_map(MAX_SORT_NSTATES), best_score, icls, label
            mapped_labels(1,:) = restart_labels(1,:)
            do label = 1, sort_nstates
                state_maps(1,label) = label
            enddo
            do irestart = 2, params%nrestarts
                call best_state_mapping(irestart, best_map, best_score)
                state_maps(irestart,:) = best_map(1:sort_nstates)
                do icls = 1, ncls
                    label = restart_labels(irestart,icls)
                    if( label < 1 )then
                        mapped_labels(irestart,icls) = 0
                    else
                        mapped_labels(irestart,icls) = best_map(label)
                    endif
                enddo
                write(logfhandle,'(A,I0,A,I0,A,3(I0,1X))') '>>> RESTART ', irestart, &
                    ' BEST STATE-LABEL AGREEMENT: ', best_score, ' MAP RAW->CONSENSUS: ', best_map
            enddo
        end subroutine map_state_correspondence

        subroutine best_state_mapping( restart_ind, best_map, best_score )
            integer, intent(in)  :: restart_ind
            integer, intent(out) :: best_map(MAX_SORT_NSTATES), best_score
            integer :: candidate(MAX_SORT_NSTATES), p1, p2, p3, score
            best_map   = 0
            best_score = -1
            select case(sort_nstates)
                case(2)
                    do p1 = 1, 2
                        do p2 = 1, 2
                            if( p2 == p1 ) cycle
                            candidate    = 0
                            candidate(1) = p1
                            candidate(2) = p2
                            score = mapping_score(restart_ind, candidate)
                            if( score > best_score )then
                                best_score = score
                                best_map   = candidate
                            endif
                        enddo
                    enddo
                case(3)
                    do p1 = 1, 3
                        do p2 = 1, 3
                            if( p2 == p1 ) cycle
                            do p3 = 1, 3
                                if( p3 == p1 .or. p3 == p2 ) cycle
                                candidate    = 0
                                candidate(1) = p1
                                candidate(2) = p2
                                candidate(3) = p3
                                score = mapping_score(restart_ind, candidate)
                                if( score > best_score )then
                                    best_score = score
                                    best_map   = candidate
                                endif
                            enddo
                        enddo
                    enddo
            end select
        end subroutine best_state_mapping

        integer function mapping_score( restart_ind, candidate_map )
            integer, intent(in) :: restart_ind, candidate_map(MAX_SORT_NSTATES)
            integer :: icls, ref_label, raw_label
            mapping_score = 0
            do icls = 1, ncls
                ref_label = restart_labels(1,icls)
                raw_label = restart_labels(restart_ind,icls)
                if( ref_label < 1 .or. raw_label < 1 ) cycle
                if( ref_label == candidate_map(raw_label) ) mapping_score = mapping_score + 1
            enddo
        end function mapping_score

        subroutine build_consensus
            integer :: icls, label, best_label, best_votes
            votes = 0
            do icls = 1, ncls
                if( original_states(icls) <= 0 ) cycle
                do irestart = 1, params%nrestarts
                    label = mapped_labels(irestart,icls)
                    if( label >= 1 .and. label <= sort_nstates ) votes(label,icls) = votes(label,icls) + 1
                enddo
                best_label = 0
                best_votes = -1
                do label = 1, sort_nstates
                    if( votes(label,icls) > best_votes )then
                        best_label = label
                        best_votes = votes(label,icls)
                    endif
                enddo
                label = mapped_labels(1,icls)
                if( label > 0 .and. votes(label,icls) == best_votes ) best_label = label
                consensus(icls) = best_label
            enddo
        end subroutine build_consensus

        subroutine assign_good_bad_states
            integer :: icls, label, pops(MAX_SORT_NSTATES)
            real    :: means(MAX_SORT_NSTATES), best_mean
            pops           = 0
            means          = 0.0
            good_consensus = 0
            best_mean      = -huge(0.0)
            do label = 1, sort_nstates
                means(label) = mean_quality_for_consensus(label, pops(label))
                write(logfhandle,'(A,I0,A,F8.3,A,I0)') '>>> CONSENSUS STATE QUALITY MEAN / POP: ', &
                    label, ' ', means(label), ' / ', pops(label)
                if( pops(label) > 0 .and. means(label) > best_mean )then
                    best_mean      = means(label)
                    good_consensus = label
                endif
            enddo
            if( good_consensus == 0 ) THROW_HARD('abinitio3D_cavgs_reject: no active consensus classes')
            do icls = 1, ncls
                if( consensus(icls) == 0 )then
                    final_states(icls) = 0
                else if( consensus(icls) == good_consensus )then
                    final_states(icls) = 1
                else
                    final_states(icls) = 2
                endif
            enddo
            write(logfhandle,'(A,I0)') '>>> CAVGS REJECT GOOD CONSENSUS STATE: ', good_consensus
            write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> CAVGS REJECT FINAL GOOD/BAD/INACTIVE: ', &
                count(final_states == 1), ' / ', count(final_states == 2), ' / ', count(final_states == 0)
        end subroutine assign_good_bad_states

        subroutine build_selection_states
            selection_states = merge(1, 0, final_states == 1)
            write(logfhandle,'(A,I6,A,I6)') '>>> CAVGS REJECT SELECTED / REJECTED : ', &
                count(selection_states > 0), ' / ', count(selection_states == 0)
        end subroutine build_selection_states

        subroutine build_consensus_volume
            call select_restart_consensus_volumes
            call dock_consensus_volumes
            call average_consensus_volume
        end subroutine build_consensus_volume

        subroutine select_restart_consensus_volumes
            type(image)  :: vol_probe
            type(string) :: restart_dir, vol_name
            integer :: nptcls
            real    :: best_score
            integer :: raw_state
            call simple_getcwd(cwd)
            CWD_GLOB = cwd%to_char()
            best_score = -huge(0.0)
            do irestart = 1, params%nrestarts
                raw_state = raw_state_for_consensus(irestart, good_consensus)
                if( raw_state == 0 )then
                    write(logfhandle,'(A,I0,A,I0)') '>>> RESTART ', irestart, &
                        ' HAS NO RAW STATE FOR CONSENSUS STATE ', good_consensus
                    THROW_HARD('abinitio3D_cavgs_reject: cannot identify restart consensus volume')
                endif
                restart_good_state(irestart)  = raw_state
                restart_good_scores(irestart) = mean_quality_for_restart_consensus(irestart, &
                    good_consensus, restart_good_pop(irestart))
                restart_dir = RESTART_DIR_FBODY//int2str_pad(irestart, 3)
                vol_name    = filepath(cwd, restart_dir, refine3D_state_vol_fname(raw_state))
                if( .not.file_exists(vol_name) )then
                    write(logfhandle,'(A,A)') '>>> MISSING RESTART CONSENSUS VOLUME: ', vol_name%to_char()
                    THROW_HARD('abinitio3D_cavgs_reject: selected restart volume missing')
                endif
                selected_vols(irestart) = vol_name
                if( restart_good_pop(irestart) > 0 .and. &
                    (ref_restart == 0 .or. restart_good_scores(irestart) > best_score) )then
                    ref_restart = irestart
                    best_score  = restart_good_scores(irestart)
                endif
                write(logfhandle,'(A,I0,A,I0,A,I0,A,F8.3,A,I0)') '>>> RESTART CONSENSUS VOLUME: ', &
                    irestart, ' RAW_STATE=', raw_state, ' CONSENSUS_STATE=', good_consensus, &
                    ' QUALITY_MEAN=', restart_good_scores(irestart), ' POP=', restart_good_pop(irestart)
            enddo
            if( ref_restart == 0 ) THROW_HARD('abinitio3D_cavgs_reject: no restart consensus volume selected')
            call find_ldim_nptcls(selected_vols(ref_restart), consensus_ldim, nptcls)
            consensus_smpd = find_img_smpd(selected_vols(ref_restart))
            call vol_probe%new(consensus_ldim, consensus_smpd)
            call vol_probe%read(selected_vols(ref_restart))
            consensus_ldim = vol_probe%get_ldim()
            consensus_smpd = vol_probe%get_smpd()
            if( consensus_smpd <= 0.0 ) consensus_smpd = params%smpd
            call vol_probe%kill
            write(logfhandle,'(A,I0,A,F8.3,A,I0)') &
                '>>> CAVGS REJECT CONSENSUS VOLUME REFERENCE RESTART: ', ref_restart, &
                ' QUALITY_MEAN=', best_score, ' RAW_STATE=', restart_good_state(ref_restart)
        end subroutine select_restart_consensus_volumes

        subroutine dock_consensus_volumes
            type(cmdline) :: cline_dock
            type(string), allocatable :: dock_wait_files(:)
            type(string) :: dock_out, dock_report, dock_done_path
            integer :: ndock, idock
            ndock = params%nrestarts - 1
            if( ndock > 0 ) allocate(dock_wait_files(ndock))
            idock = 0
            do irestart = 1, params%nrestarts
                if( irestart == ref_restart )then
                    docked_vols(irestart)  = selected_vols(irestart)
                    dock_reports(irestart) = 'reference'
                    cycle
                endif
                folder      = DOCK_DIR_FBODY//int2str_pad(irestart, 3)
                dock_out    = 'consensus_docked_restart_'//int2str_pad(irestart, 3)//MRC_EXT
                dock_report = 'consensus_dock_restart_'//int2str_pad(irestart, 3)//TXT_EXT
                dock_done_path = filepath(folder, string(DOCK_DONE))
                call simple_mkdir(folder)
                if( file_exists(dock_done_path) ) call del_file(dock_done_path)
                docked_vols(irestart)  = filepath(cwd, folder, dock_out)
                dock_reports(irestart) = filepath(cwd, folder, dock_report)
                if( file_exists(docked_vols(irestart)) ) call del_file(docked_vols(irestart))
                idock = idock + 1
                dock_wait_files(idock) = filepath(cwd, folder, string(DOCK_DONE))
                call cline_dock%kill
                call cline_dock%set('prg',                'dock_volpair')
                call cline_dock%set('vol1',               selected_vols(ref_restart))
                call cline_dock%set('vol2',               selected_vols(irestart))
                call cline_dock%set('outvol',             dock_out)
                call cline_dock%set('outfile',            dock_report)
                call cline_dock%set('smpd',               consensus_smpd)
                call cline_dock%set('hp',                 DEFAULT_DOCK_HP)
                call cline_dock%set('lp',                 dock_lp)
                call cline_dock%set('mskdiam',            params%mskdiam)
                call cline_dock%set('nthr',               params%nthr)
                call cline_dock%set('mkdir',              'no')
                call cline_dock%set('verbose_exit',       'yes')
                call cline_dock%set('verbose_exit_fname', DOCK_DONE)
                call simple_chdir(folder)
                call simple_getcwd(cwd_run)
                CWD_GLOB = cwd_run%to_char()
                call qenv%exec_simple_prg_in_queue_async(cline_dock, &
                    string('abinitio3D_cavgs_reject_dock_script'), string('abinitio3D_cavgs_reject_dock.log'))
                call simple_chdir(cwd)
                CWD_GLOB = cwd%to_char()
                write(logfhandle,'(A,I0,A,A,A,F8.3,A,F8.1)') '>>> SUBMITTED CONSENSUS VOLUME DOCK ', &
                    irestart, ' IN ', folder%to_char(), ' LP=', dock_lp, ' HP=', DEFAULT_DOCK_HP
            enddo
            call cline_dock%kill
            if( ndock > 0 )then
                call qsys_watcher_diag(dock_wait_files)
                do idock = 1, ndock
                    if( .not.file_exists(dock_wait_files(idock)) )then
                        write(logfhandle,'(A,A)') '>>> MISSING DOCK COMPLETION MARKER: ', &
                            dock_wait_files(idock)%to_char()
                        THROW_HARD('abinitio3D_cavgs_reject: one or more volume docks did not finish')
                    endif
                enddo
            endif
            do irestart = 1, params%nrestarts
                if( .not.file_exists(docked_vols(irestart)) )then
                    write(logfhandle,'(A,A)') '>>> MISSING DOCKED CONSENSUS VOLUME: ', &
                        docked_vols(irestart)%to_char()
                    THROW_HARD('abinitio3D_cavgs_reject: docked consensus volume missing')
                endif
            enddo
            if( allocated(dock_wait_files) ) deallocate(dock_wait_files)
        end subroutine dock_consensus_volumes

        subroutine average_consensus_volume
            type(image) :: avg_vol, add_vol
            integer :: ldim_here(3)
            real    :: smpd_here
            call avg_vol%new(consensus_ldim, consensus_smpd)
            call avg_vol%read(docked_vols(ref_restart))
            do irestart = 1, params%nrestarts
                if( irestart == ref_restart ) cycle
                call add_vol%kill
                call add_vol%new(consensus_ldim, consensus_smpd)
                call add_vol%read(docked_vols(irestart))
                ldim_here = add_vol%get_ldim()
                if( any(ldim_here /= consensus_ldim) )then
                    write(logfhandle,'(A,A)') '>>> NONCONFORMING DOCKED VOLUME: ', &
                        docked_vols(irestart)%to_char()
                    THROW_HARD('abinitio3D_cavgs_reject: docked volumes have inconsistent dimensions')
                endif
                smpd_here = add_vol%get_smpd()
                if( smpd_here > 0.0 .and. abs(smpd_here - consensus_smpd) > 1.0e-4 )then
                    write(logfhandle,'(A,A)') '>>> NONCONFORMING DOCKED VOLUME SAMPLING: ', &
                        docked_vols(irestart)%to_char()
                    THROW_HARD('abinitio3D_cavgs_reject: docked volumes have inconsistent sampling')
                endif
                call avg_vol%add(add_vol)
            enddo
            call avg_vol%div(real(params%nrestarts))
            call avg_vol%write(string(CONSENSUS_VOL), del_if_exists=.true.)
            call spproj%add_vol2os_out(string(CONSENSUS_VOL), consensus_smpd, 1, 'vol_cavg', &
                box=consensus_ldim(1))
            write(logfhandle,'(A,A)') '>>> WROTE CAVGS REJECT CONSENSUS VOLUME: ', CONSENSUS_VOL
            call avg_vol%kill
            call add_vol%kill
        end subroutine average_consensus_volume

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

        real function mean_quality_for_restart_consensus( restart_ind, label, pop )
            integer, intent(in)  :: restart_ind, label
            integer, intent(out) :: pop
            integer :: icls
            mean_quality_for_restart_consensus = 0.0
            pop = 0
            do icls = 1, ncls
                if( mapped_labels(restart_ind,icls) /= label ) cycle
                if( original_states(icls) <= 0 ) cycle
                pop = pop + 1
                mean_quality_for_restart_consensus = mean_quality_for_restart_consensus + quality_scores(icls)
            enddo
            if( pop > 0 ) mean_quality_for_restart_consensus = &
                mean_quality_for_restart_consensus / real(pop)
        end function mean_quality_for_restart_consensus

        integer function raw_state_for_consensus( restart_ind, label )
            integer, intent(in) :: restart_ind, label
            integer :: raw_state
            raw_state_for_consensus = 0
            do raw_state = 1, sort_nstates
                if( state_maps(restart_ind,raw_state) == label )then
                    raw_state_for_consensus = raw_state
                    return
                endif
            enddo
        end function raw_state_for_consensus

        subroutine annotate_project
            integer :: icls, ncls3d, max_vote
            ncls3d = spproj%os_cls3D%get_noris()
            do icls = 1, ncls
                max_vote = maxval(votes(1:sort_nstates,icls))
                call spproj%os_cls2D%set(icls, 'quality',             quality_scores(icls))
                call spproj%os_cls2D%set(icls, 'accept',              selection_states(icls))
                call spproj%os_cls2D%set(icls, 'quality_cluster',     quality%labels(icls))
                call spproj%os_cls2D%set(icls, 'cavgs_reject_consensus', consensus(icls))
                call spproj%os_cls2D%set(icls, 'cavgs_reject_votes',     max_vote)
                if( ncls3d == ncls )then
                    call spproj%os_cls3D%set(icls, 'quality',             quality_scores(icls))
                    call spproj%os_cls3D%set(icls, 'accept',              selection_states(icls))
                    call spproj%os_cls3D%set(icls, 'quality_cluster',     quality%labels(icls))
                    call spproj%os_cls3D%set(icls, 'cavgs_reject_consensus', consensus(icls))
                    call spproj%os_cls3D%set(icls, 'cavgs_reject_votes',     max_vote)
                endif
            enddo
        end subroutine annotate_project

        subroutine write_consensus_report
            integer :: funit, icls, label
            open(newunit=funit, file=CONSENSUS_REPORT, status='replace', action='write')
            write(funit,'(A)') '# abinitio3D_cavgs_reject consensus report'
            write(funit,'(A,I0)') '# nrestarts=', params%nrestarts
            write(funit,'(A,I0)') '# nstates=', sort_nstates
            write(funit,'(A,I0)') '# nstages=', SORT_NSTAGES
            write(funit,'(A,A)') '# quality_model=', trim(model%name)
            write(funit,'(A)', advance='no') 'class,original_state,consensus_state,final_state,selection_state'
            do label = 1, sort_nstates
                write(funit,'(A,I0)', advance='no') ',votes_state_', label
            enddo
            write(funit,'(A)', advance='no') ',quality_score,quality_auto_state'
            do irestart = 1, params%nrestarts
                write(funit,'(A,I0)', advance='no') ',restart_raw_', irestart
            enddo
            do irestart = 1, params%nrestarts
                write(funit,'(A,I0)', advance='no') ',restart_mapped_', irestart
            enddo
            write(funit,*)
            do icls = 1, ncls
                write(funit,'(I0,A,I0,A,I0,A,I0,A,I0)', advance='no') &
                    icls, ',', original_states(icls), ',', consensus(icls), ',', final_states(icls), ',', selection_states(icls)
                do label = 1, sort_nstates
                    write(funit,'(A,I0)', advance='no') ',', votes(label,icls)
                enddo
                write(funit,'(A,ES14.6,A,I0)', advance='no') ',', quality_scores(icls), ',', quality_auto_states(icls)
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
        end subroutine write_consensus_report

        subroutine write_consensus_volume_report
            integer :: funit
            open(newunit=funit, file=CONSENSUS_VOL_REPORT, status='replace', action='write')
            write(funit,'(A)') '# abinitio3D_cavgs_reject consensus volume report'
            write(funit,'(A,A)') '# consensus_volume=', CONSENSUS_VOL
            write(funit,'(A,I0)') '# reference_restart=', ref_restart
            write(funit,'(A,I0)') '# good_consensus_state=', good_consensus
            write(funit,'(A,I0)') '# box=', consensus_ldim(1)
            write(funit,'(A,F10.5)') '# smpd=', consensus_smpd
            write(funit,'(A,F8.3)') '# dock_lp=', dock_lp
            write(funit,'(A,F8.3)') '# dock_hp=', DEFAULT_DOCK_HP
            write(funit,'(A,A)') &
                'restart,raw_state,mapped_consensus_state,state_pop,state_quality_mean,role,', &
                'selected_volume,docked_volume,dock_report'
            do irestart = 1, params%nrestarts
                write(funit,'(I0,A,I0,A,I0,A,I0,A,ES14.6,A)', advance='no') &
                    irestart, ',', restart_good_state(irestart), ',', good_consensus, ',', &
                    restart_good_pop(irestart), ',', restart_good_scores(irestart), ','
                if( irestart == ref_restart )then
                    write(funit,'(A)', advance='no') 'reference'
                else
                    write(funit,'(A)', advance='no') 'docked'
                endif
                write(funit,'(A,A,A,A,A,A)') ',', selected_vols(irestart)%to_char(), ',', &
                    docked_vols(irestart)%to_char(), ',', dock_reports(irestart)%to_char()
            enddo
            close(funit)
            write(logfhandle,'(A,A)') '>>> WROTE ', CONSENSUS_VOL_REPORT
        end subroutine write_consensus_volume_report

        subroutine write_selection_outputs
            call write_cavg_quality_feature_table(quality, model, 'cavgs_quality_features.txt', &
                params%projfile%to_char(), manual_states=selection_states)
            call write_rejection_stack(string('quality_selected_cavgs'//MRC_EXT),  selected=.true.)
            call write_rejection_stack(string('quality_rejected_cavgs'//MRC_EXT), selected=.false.)
        end subroutine write_selection_outputs

        subroutine write_rejection_stack( fname, selected )
            type(string), intent(in) :: fname
            logical,      intent(in) :: selected
            integer :: icls, istk
            if( file_exists(fname) ) call del_file(fname)
            istk = 0
            do icls = 1, ncls
                if( selected .eqv. (selection_states(icls) > 0) )then
                    istk = istk + 1
                    call cavg_imgs(icls)%write(fname, istk)
                endif
            enddo
            write(logfhandle,'(A,A,A,I6)') '>>> WROTE ', fname%to_char(), ' #CAVGS: ', istk
        end subroutine write_rejection_stack

        subroutine cleanup
            call spproj%kill
            call restart_proj%kill
            call quality%kill
            call dealloc_imgarr(cavg_imgs)
            if( allocated(restart_labels)      ) deallocate(restart_labels)
            if( allocated(mapped_labels)       ) deallocate(mapped_labels)
            if( allocated(state_maps)          ) deallocate(state_maps)
            if( allocated(consensus)           ) deallocate(consensus)
            if( allocated(final_states)        ) deallocate(final_states)
            if( allocated(selection_states)    ) deallocate(selection_states)
            if( allocated(original_states)     ) deallocate(original_states)
            if( allocated(quality_auto_states) ) deallocate(quality_auto_states)
            if( allocated(quality_scores)      ) deallocate(quality_scores)
            if( allocated(votes)               ) deallocate(votes)
            if( allocated(restart_good_state)  ) deallocate(restart_good_state)
            if( allocated(restart_good_pop)    ) deallocate(restart_good_pop)
            if( allocated(restart_good_scores) ) deallocate(restart_good_scores)
            if( allocated(restart_projfiles)   ) deallocate(restart_projfiles)
            if( allocated(done_files)          ) deallocate(done_files)
            if( allocated(selected_vols)       ) deallocate(selected_vols)
            if( allocated(docked_vols)         ) deallocate(docked_vols)
            if( allocated(dock_reports)        ) deallocate(dock_reports)
        end subroutine cleanup

    end subroutine exec_abinitio3D_cavgs_reject

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
        logical :: l_cavg_ini_ext, l_vol_ini_ext, l_user_nstages, l_user_lpstop, l_run_final_rec
        logical :: l_state_continue
        logical :: l_force_full_sampling
        real    :: sampled_active_frac
        l_state_continue = cline%defined('state')
        l_force_full_sampling = .false.
        sampled_active_frac   = 0.
        l_state_continue_mode = .false.
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
        if( cline%defined('nsample_start') .or. cline%defined('nsample_stop') )then
            THROW_HARD('nsample_start/nsample_stop are no longer supported for abinitio3D; set nsample instead')
        endif
        if( l_state_continue )then
            if( cline%defined('multivol_mode') )then
                if( cline%get_carg('multivol_mode').ne.'single' )then
                    THROW_HARD('abinitio3D state continuation requires multivol_mode=single')
                endif
            endif
            call cline%set('multivol_mode', 'single')
            call cline%set('filt_mode',     'nonuniform')
        endif
        l_user_nstages = cline%defined('nstages')
        l_user_lpstop  = cline%defined('lpstop')
        ! splitting stage
        split_stage = abinitio_het_docked_stage()
        if( cline%defined('split_stage') ) split_stage = cline%get_iarg('split_stage')
        if( split_stage < 2 .or. split_stage > abinitio_nstages() )then
            THROW_HARD('split_stage must be between 2 and '//int2str(abinitio_nstages())//' for abinitio3D')
        endif
        call cline%set('split_stage', split_stage)
        ! adjust default multivol_mode unless given on command line
        if( cline%defined('nstates') )then
            nstates_on_cline = cline%get_iarg('nstates')
            if( nstates_on_cline > 1 .and. .not. cline%defined('multivol_mode') )then
                call cline%set('multivol_mode', 'independent')
            endif
        endif
        if( cline%defined('multivol_mode') )then
            if( cline%get_carg('multivol_mode').eq.'independent' )then
                ! Stop independent multi-state starts before prob_neigh/NU by default.
                if( .not. l_user_nstages ) call cline%set('nstages', abinitio_independent_nstages_default())
                if( .not. l_user_lpstop  ) call cline%set('lpstop',  abinitio_independent_lpstop_default())
            endif
        endif
        ! make master parameters
        call params%new(cline)
        write(logfhandle,'(A,A)') '>>> ABINITIO3D PARTICLE SOURCE: ', trim(params%ptcl_src)
        l_state_continue_mode = l_state_continue
        if( trim(params%multivol_mode).eq.'independent' )then
            if( .not. l_user_nstages ) write(logfhandle,'(A,I0)') &
                &'>>> ABINITIO3D INDEPENDENT MULTI-STATE DEFAULT NSTAGES: ', params%nstages
            if( .not. l_user_lpstop ) write(logfhandle,'(A,F4.1,A)') &
                &'>>> ABINITIO3D INDEPENDENT MULTI-STATE DEFAULT LPSTOP: ', params%lpstop, ' A'
        endif
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
            case('independent', 'docked')
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
        if( l_state_continue )then
            if( trim(params%cavg_ini).eq.'yes' .or. trim(params%cavg_ini_ext).eq.'yes' )then
                THROW_HARD('abinitio3D state continuation cannot be combined with cavg_ini/cavg_ini_ext')
            endif
            if( cline%defined('vol1') )then
                THROW_HARD('abinitio3D state continuation uses the selected project state; remove vol1')
            endif
            call prepare_state_continue_project
            call cline%set('pgrp_start', params%pgrp)
            params%pgrp_start = params%pgrp
            start_stage = abinitio_independent_nstages_default()
            l_ini3D     = .true.
            write(logfhandle,'(A,I0,A)') &
                &'>>> ABINITIO3D STATE CONTINUATION STARTING FROM STAGE ', start_stage, ' WITH NONUNIFORM FILTERING'
        endif
        if( trim(params%cavg_ini).eq.'yes' )then
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
            ! input volumes are assumed aligned to the target symmetry axis
            call cline%set('pgrp_start', params%pgrp)
            params%pgrp_start = params%pgrp
            ! setting up random classes for particles sampling
            call spproj%os_ptcl2D%rnd_cls(100)
            call spproj%write_segment_inside('ptcl2D', params%projfile)
            call spproj%os_cls2D%new(100, is_ptcl=.false.)
            call spproj%os_cls2D%set_all2single('state', 1)
        endif
        ! set class global filtering flags for staged refine3D policy
        l_nonuniform = params%l_nonuniform
        nstages_refine3D = last_stage
        if( nstages_refine3D < start_stage )then
            THROW_HARD('nstages must be >= first executable abinitio3D stage')
        endif
        l_run_final_rec = nstages_refine3D == abinitio_nstages() .or. trim(params%multivol_mode).eq.'independent'
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
            if( nptcls_eff < 1 ) THROW_HARD('No active particles selected in ptcl2D for abinitio3D')
            if( .not. cline%defined('nsample') ) params%nsample = abinitio_nsample_default()
            if( params%nsample < 1 ) THROW_HARD('nsample must be >= 1 for abinitio3D sampled update')
            sampled_active_frac   = real(params%nsample) / real(nptcls_eff)
            l_force_full_sampling = sampled_active_frac > abinitio_full_sample_switch_frac()
            if( l_force_full_sampling )then
                update_frac = 1.0
                write(logfhandle,'(A,F8.4,A,F8.4,A)') &
                    &'>>> ABINITIO3D NSAMPLE/ACTIVE FRACTION ', sampled_active_frac, ' > ', &
                    &abinitio_full_sample_switch_frac(), ' -> FORCING FULL ACTIVE SAMPLING (NO FRACTIONAL OR TRAILING UPDATE)'
            else
                update_frac = real(params%nsample * params%nstates) / real(nptcls_eff)
                update_frac = min(abinitio_update_frac_max(), update_frac) ! keep fractional update on below the switch threshold
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
            endif
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
        if( .not. l_ini3D )then
            call reset_ptcl3D_from_ptcl2D_selection
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
            if( l_force_full_sampling )then
                call spproj%os_ptcl3D%sample4update_all([1,noris], nptcls2update, pinds, .true.)
            else
                call spproj%os_ptcl3D%sample4update_class(clssmp, [1,noris], update_frac, nptcls2update, pinds, .true., .true.)
            endif
            call spproj%os_ptcl3D%set_updatecnt(1, pinds) ! set all sampled updatecnts to 1 & the rest to zero
            deallocate(pinds)                             ! these are not needed
            if( allocated(clssmp) ) call deallocate_class_samples(clssmp) ! done with this one
            ! write updated project file
            call spproj%write_segment_inside(params%oritype, params%projfile)
            ! create starting volume(s)
            call calc_rec(params, params%projfile, xrec3D, start_stage)
        endif
        if( cline%defined('nstages') )then
            write(logfhandle,'(A,I0,A,I0)')'>>> ABINITIO3D STAGE RANGE: ', start_stage, ' -> ', nstages_refine3D
            if( nstages_refine3D < abinitio_nstages() )then
                if( l_run_final_rec )then
                    write(logfhandle,'(A)')'>>> ABINITIO3D EARLY STAGE STOP: FINAL ALL-PARTICLE RECONSTRUCTION ENABLED'
                else
                    write(logfhandle,'(A)')'>>> ABINITIO3D EARLY STAGE STOP: SKIPPING FINAL ALL-PARTICLE RECONSTRUCTION'
                endif
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
            ! Splitting stage of docked mode
            if( params%multivol_mode.eq.'docked' .and. istage == split_stage )then
                ! map all particles to a projection direction
                call ensure_docked_multistate_particle_assignments
                ! reset the nstates in params, update sampling
                params%nstates = nstates_glob
                if( l_force_full_sampling )then
                    update_frac = 1.0
                else
                    update_frac = real(params%nsample * params%nstates) / real(nptcls_eff)
                    update_frac = min(abinitio_update_frac_max(), update_frac)
                endif
                write(logfhandle,'(A,I0,A,I0,A,F8.4)') &
                    &'>>> ABINITIO3D DOCKED SPLIT STAGE/NSTATES/POSTSPLIT_UPDATE_FRAC: ', &
                    &split_stage, '/', params%nstates, '/', update_frac
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
        if( l_run_final_rec )then
            select case(trim(params%multivol_mode))
                case('independent','docked')
                    call ensure_multistate_particle_assignments
            end select
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

        subroutine clean_ptcl3D_sampling
            call spproj%os_ptcl3D%clean_entry('updatecnt', 'sampled')
        end subroutine clean_ptcl3D_sampling

        subroutine prepare_state_continue_project
            type(commander_selection) :: xselection
            type(cmdline)             :: cline_selection
            type(string)              :: src_projfile, work_projfile, work_projname
            integer                   :: nselected
            if( params%state < 1 ) THROW_HARD('abinitio3D state continuation requires state >= 1')
            nselected = spproj%get_n_insegment_state('ptcl3D', params%state)
            if( nselected < 1 )then
                THROW_HARD('requested abinitio3D continuation state is absent from ptcl3D')
            endif
            if( spproj%is_virgin_field('ptcl3D') )then
                THROW_HARD('abinitio3D state continuation requires existing ptcl3D orientations')
            endif
            src_projfile  = params%projfile
            work_projfile = 'abinitio3D_state'//int2str_pad(params%state,2)//'_tmpproj.simple'
            work_projname = get_fbody(work_projfile,'simple')
            if( file_exists(work_projfile) ) call del_file(work_projfile)
            call simple_copy_file(src_projfile, work_projfile)
            cline_selection = cline
            call cline_selection%set('prg',      'selection')
            call cline_selection%set('projfile', work_projfile)
            call cline_selection%set('projname', work_projname)
            call cline_selection%set('oritype',  'ptcl3D')
            call cline_selection%set('state',    params%state)
            call cline_selection%set('prune',    'yes')
            call cline_selection%set('append',   'no')
            call cline_selection%set('mkdir',    'no')
            call xselection%execute(cline_selection)
            call cline%set('projfile', work_projfile)
            call cline%set('projname', work_projname)
            params%projfile = work_projfile
            params%projname = work_projname
            call spproj%read(params%projfile)
            call spproj%update_projinfo(params%projfile)
            call spproj%write(params%projfile)
            write(logfhandle,'(A,I0,A,I0,A,A)') &
                &'>>> ABINITIO3D STATE CONTINUATION STATE/PARTICLES: ', params%state, &
                &' / ', nselected, ' TEMP PROJECT: ', params%projfile%to_char()
            call cline_selection%kill
            call src_projfile%kill
            call work_projfile%kill
            call work_projname%kill
        end subroutine prepare_state_continue_project

        subroutine reset_ptcl3D_from_ptcl2D_selection
            integer :: iptcl, nptcls2D, nptcls3D, state2D, nactive
            nptcls2D = spproj%os_ptcl2D%get_noris()
            nptcls3D = spproj%os_ptcl3D%get_noris()
            if( nptcls2D /= nptcls3D )then
                THROW_HARD('Inconsistent number of particles in PTCL2D/PTCL3D segments; abinitio3D')
            endif
            if( .not. spproj%os_ptcl2D%isthere('state') )then
                THROW_HARD('state flag missing from ptcl2D; abinitio3D')
            endif
            call clean_ptcl3D_sampling
            call spproj%os_ptcl3D%delete_3Dalignment(keepshifts=.true.)
            call spproj%os_ptcl3D%transfer_2Dshifts(spproj%os_ptcl2D)
            nactive = 0
            do iptcl = 1,nptcls3D
                state2D = spproj%os_ptcl2D%get_state(iptcl)
                if( state2D > 0 )then
                    call spproj%os_ptcl3D%set_state(iptcl, 1)
                    nactive = nactive + 1
                else
                    call spproj%os_ptcl3D%set_state(iptcl, 0)
                endif
            enddo
            if( nactive < 1 ) THROW_HARD('No active particles selected in ptcl2D for abinitio3D')
        end subroutine reset_ptcl3D_from_ptcl2D_selection

        subroutine ensure_multistate_particle_assignments
            integer :: nactive, nupdated, nmissing
            call read_multistate_assignment_coverage(nactive, nupdated, nmissing)
            if( nactive < 1 )then
                THROW_HARD('multistate abinitio3D has no active particles after staged refinement')
            endif
            if( nmissing > 0 )then
                call run_multistate_missing_update(nmissing, nactive)
                call read_multistate_assignment_coverage(nactive, nupdated, nmissing)
                if( nmissing > 0 )then
                    THROW_HARD('multistate abinitio3D final missing-update pass failed to update every active particle')
                endif
            endif
        end subroutine ensure_multistate_particle_assignments

        subroutine ensure_docked_multistate_particle_assignments
            type(cmdline) :: cline_missing
            integer       :: nactive, nupdated, nmissing, iter_missing
            call read_multistate_assignment_coverage(nactive, nupdated, nmissing)
            if( nactive < 1 )then
                THROW_HARD('multistate abinitio3D has no active particles after staged refinement')
            endif
            if( nmissing > 0 )then
                iter_missing = next_refine3D_iteration()
                write(logfhandle,'(A,A,I0,A,I0,A,I0)') &
                &'>>> ABINITIO3D DOCKED MULTISTATE MISSING-UPDATE ASSIGNMENT', &
                &' MISSING/ACTIVE/ITER: ', nmissing, '/', nactive, '/', iter_missing
                call flush(logfhandle)
                cline_missing = cline_refine3D
                call cline_missing%set('prg',          'refine3D')
                call cline_missing%set('mkdir',              'no')
                call cline_missing%set('refine',           'prob')
                call cline_missing%set('balance',            'no')
                call cline_missing%set('frac_best',           1.0)
                call cline_missing%set('fillin',             'no')
                call cline_missing%set('update_frac',         1.0)
                call cline_missing%set('trail_rec',          'no')
                call cline_missing%set('volrec',             'no')
                call cline_missing%set('maxits',                1)
                call cline_missing%set('startit',    iter_missing)
                call cline_missing%set('which_iter', iter_missing)
                call cline_missing%set('extr_iter',  iter_missing)
                call cline_missing%delete('endit')
                call cline_missing%delete('greedy_sampling')
                call xrefine3D%execute(cline_missing)
                call del_files(DIST_FBODY,      params%nparts, ext='.dat')
                call del_files(ASSIGNMENT_FBODY,params%nparts, ext='.dat')
                call del_file(DIST_FBODY//'.dat')
                call del_file(ASSIGNMENT_FBODY//'.dat')
                call cline_missing%kill
                call read_multistate_assignment_coverage(nactive, nupdated, nmissing)
                if( nmissing > 0 )then
                    THROW_HARD('multistate abinitio3D final missing-update pass failed to update every active particle')
                endif
            endif
        end subroutine ensure_docked_multistate_particle_assignments

        subroutine read_multistate_assignment_coverage( nactive, nupdated, nmissing )
            integer, intent(out) :: nactive, nupdated, nmissing
            integer, allocatable :: states(:), updatecnts(:)
            call spproj%read_segment('ptcl3D', params%projfile)
            if( .not. spproj%os_ptcl3D%isthere('updatecnt') )then
                THROW_HARD('multistate abinitio3D requires post-label particle assignments before final reconstruction')
            endif
            states     = spproj%os_ptcl3D%get_all_asint('state')
            updatecnts = spproj%os_ptcl3D%get_all_asint('updatecnt')
            nactive    = count(states > 0)
            nupdated   = count(states > 0 .and. updatecnts > 0)
            nmissing   = nactive - nupdated
            write(logfhandle,'(A,A,A,I0,A,I0,A,I0)') &
                &'>>> ABINITIO3D MULTISTATE ASSIGNMENT COVERAGE MODE=', trim(params%multivol_mode), &
                &' UPDATED/ACTIVE/MISSING: ', nupdated, '/', nactive, '/', nmissing
            if( allocated(states)     ) deallocate(states)
            if( allocated(updatecnts) ) deallocate(updatecnts)
        end subroutine read_multistate_assignment_coverage

        subroutine run_multistate_missing_update( nmissing, nactive )
            integer, intent(in) :: nmissing, nactive
            type(cmdline) :: cline_missing
            integer       :: iter_missing
            iter_missing = next_refine3D_iteration()
            write(logfhandle,'(A,A,A,I0,A,I0,A,I0)') &
                &'>>> ABINITIO3D MULTISTATE FINAL MISSING-UPDATE GREEDY ASSIGNMENT MODE=', trim(params%multivol_mode), &
                &' MISSING/ACTIVE/ITER: ', nmissing, '/', nactive, '/', iter_missing
            call flush(logfhandle)
            cline_missing = cline_refine3D
            call cline_missing%set('prg',             'refine3D')
            call cline_missing%set('mkdir',                 'no')
            call cline_missing%set('refine',            'greedy')
            call cline_missing%set('balance',               'no')
            call cline_missing%set('greedy_sampling',      'yes')
            call cline_missing%set('frac_best',             1.0)
            call cline_missing%set('fillin',               'no')
            call cline_missing%set('update_missing',       'yes')
            call cline_missing%set('update_frac',           1.0)
            call cline_missing%set('trail_rec',             'no')
            call cline_missing%set('volrec',                'no')
            call cline_missing%set('maxits',                   1)
            call cline_missing%set('startit',       iter_missing)
            call cline_missing%set('which_iter',    iter_missing)
            call cline_missing%set('extr_iter',     iter_missing)
            call cline_missing%delete('endit')
            call xrefine3D%execute(cline_missing)
            call del_files(DIST_FBODY,      params%nparts, ext='.dat')
            call del_files(ASSIGNMENT_FBODY,params%nparts, ext='.dat')
            call del_file(DIST_FBODY//'.dat')
            call del_file(ASSIGNMENT_FBODY//'.dat')
            call cline_missing%kill
        end subroutine run_multistate_missing_update

        integer function next_refine3D_iteration() result(iter)
            iter = 1
            if( cline_refine3D%defined('endit') )then
                iter = cline_refine3D%get_iarg('endit') + 1
            else if( cline_refine3D%defined('which_iter') )then
                iter = cline_refine3D%get_iarg('which_iter') + 1
            endif
            iter = max(1, iter)
        end function next_refine3D_iteration

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

end module simple_commanders_abinitio
