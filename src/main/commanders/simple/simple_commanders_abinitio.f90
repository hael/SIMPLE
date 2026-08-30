!@descr: abinitio 3D reconstruction in single- and multi-particle mode
module simple_commanders_abinitio
use simple_commanders_api
use simple_abinitio_utils
use simple_procimgstk,              only: shift_imgfile
use simple_commanders_project_core, only: commander_selection
use simple_commanders_reproject,    only: commander_reproject
use simple_commanders_refine3D,     only: commander_refine3D, commander_refine3D
use simple_commanders_rec,          only: commander_rec3D, commander_rec3D
use simple_cluster_seed,            only: gen_labelling
use simple_refine3D_fnames,         only: refine3D_startvol_fname, refine3D_startvol_half_fname, &
    &refine3D_state_vol_fname, refine3D_state_halfvol_fname
implicit none

public :: commander_abinitio3D_cavgs, commander_abinitio3D
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_abinitio3D_cavgs
    contains
    procedure :: execute => exec_abinitio3D_cavgs
end type commander_abinitio3D_cavgs

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
        integer                   :: nstates_on_cline, nstates_target, split_stage
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
        call cline%set('objfun', 'euclid') ! noise normalized Euclidean distances from the start
        if( .not. cline%defined('overlap')          ) call cline%set('overlap',                     0.95)
        if( .not. cline%defined('prob_athres')      ) call cline%set('prob_athres',                  90.) ! reduces # failed runs on trpv1 from 4->2/10
        if( .not. cline%defined('cenlp')            ) call cline%set('cenlp',   abinitio_cenlp_default())
        if( .not. cline%defined('imgkind')          ) call cline%set('imgkind',                   'cavg')
        if( .not. cline%defined('filt_mode')        ) call cline%set('filt_mode',                 'none')
        if( .not. cline%defined('noise_norm')       ) call cline%set('noise_norm',                  'no')
        if( .not. cline%defined('lpstart')          ) call cline%set('lpstart', abinitio_lpstart_ini3D())
        if( .not. cline%defined('lpstop')           ) call cline%set('lpstop',   abinitio_lpstop_ini3D())
        if( .not. cline%defined('gauref')           ) call cline%set('gauref',                     'yes')
        ! splitting stage
        split_stage = abinitio_het_docked_stage()
        if( cline%defined('split_stage') ) split_stage = cline%get_iarg('split_stage')
        if( split_stage < 2 .or. split_stage > abinitio_nstages_ini3D_max() )then
            THROW_HARD('split_stage must be between 2 and '//int2str(abinitio_nstages_ini3D_max())//' for abinitio3D_cavgs')
        endif
        call cline%set('split_stage', split_stage)
        ! adjust default multivol_mode unless given on command line
        if( cline%defined('nstates') )then
            nstates_on_cline = cline%get_iarg('nstates')
            if( nstates_on_cline > 1 .and. .not. cline%defined('multivol_mode') )then
                call cline%set('multivol_mode', 'independent')
            endif
        endif
        ! make master parameters
        call params%new(cline)
        nstates_target = params%nstates
        nstates_glob   = nstates_target
        select case(trim(params%multivol_mode))
            case('single')
                if( nstates_target /= 1 ) THROW_HARD('nstates /= 1 incompatible with multivol_mode:' //trim(params%multivol_mode))
            case('independent', 'docked')
                if( nstates_target == 1 ) THROW_HARD('nstates == 1 incompatible with multivol_mode: '//trim(params%multivol_mode))
            case DEFAULT
                THROW_HARD('Unsupported multivol_mode: '//trim(params%multivol_mode))
        end select
        if( trim(params%multivol_mode).eq.'docked' )then
            params%nstates = 1
            call cline%delete('nstates')
        endif
        call cline%set('mkdir',       'no')   ! to avoid nested directory structure
        call cline%set('oritype', 'ptcl3D')   ! from now on we are in the ptcl3D segment, final report is in the cls3D segment
        params%oritype = 'ptcl3D'
        ! set work projfile
        work_projfile = 'abinitio3D_cavgs_tmpproj.simple'
        ! set class global filtering flags for staged refine3D policy
        l_nonuniform = .false.
        ! set nstages_ini3D
        nstages_ini3D = abinitio_nstages_ini3D_max()
        if( cline%defined('nstages') )then
            nstages_ini3D = min(abinitio_nstages_ini3D_max(),params%nstages)
        endif
        if( trim(params%multivol_mode).eq.'docked' .and. nstages_ini3D < split_stage )then
            THROW_HARD('multivol_mode=docked requires nstages >= split_stage for abinitio3D_cavgs')
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
                l_cavgs_mode = .true.
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
        if( trim(params%multivol_mode).eq.'docked' )then
            where( states > 0 ) states = 1
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
            ! Splitting stage of docked mode
            if( trim(params%multivol_mode).eq.'docked' )then
                if( istage == split_stage-1 )then
                    write(logfhandle,'(A,I0,A,I0)') &
                        &'>>> ABINITIO3D_CAVGS DOCKED PRE-SPLIT STAGE/NSTATES: ', istage, '/', params%nstates
                else if( istage == split_stage )then
                    params%nstates = nstates_target
                    write(logfhandle,'(A,I0,A,I0)') &
                        &'>>> ABINITIO3D_CAVGS DOCKED SPLIT STAGE/NSTATES: ', istage, '/', params%nstates
                endif
            endif
            ! Preparation of command line for probabilistic search
            call set_cline_refine3D(params, istage, l_cavgs=.true.)
            if( trim(params%multivol_mode).eq.'docked' .and. istage == split_stage )then
                call randomize_states(params, work_proj, work_projfile, xrec3D, split_stage)
            endif
            if( cline_refine3D%get_iarg('box_crop') < params%box )then
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',&
                    &cline_refine3D%get_iarg('box_crop')
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
        ! remove postprocessed (pproc) volumes; with bfac=0 they add nothing in cavgs mode
        call del_pproc_vols
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

            subroutine del_pproc_vols
                type(string), allocatable :: pproc_list(:)
                ! covers per-state, per-iteration and mirrored pproc volumes
                call simple_list_files(VOL_FBODY//'*'//PPROC_SUFFIX//'*'//MRC_EXT, pproc_list)
                if( allocated(pproc_list) ) call del_files(pproc_list)
            end subroutine del_pproc_vols

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
                call local_cline_rec%delete('objfun_den')
                call local_cline_rec%delete('objfun_den_w')
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

    !> for generation of an initial 3d model from particles
    subroutine exec_abinitio3D( self, cline )
        class(commander_abinitio3D), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        ! commanders
        type(commander_refine3D)        :: xrefine3D
        type(commander_rec3D)           :: xrec3D
        ! other
        integer,            allocatable :: tmpinds(:), clsinds(:), pinds(:), cls_states(:)
        type(class_sample), allocatable :: clssmp(:)
        type(parameters)                :: params
        type(sp_project)                :: spproj
        type(simple_nice_comm)          :: nice_comm
        real    :: lprange(2)
        integer :: state, istage, icls, start_stage, nptcls2update, noris, nstates_on_cline
        integer :: nstates_in_project, split_stage, last_stage, nptcls_cap
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
        ! nu_refine stays off by default: the stage ladder owns frequency
        ! marching in abinitio3D. An explicit nu_refine=yes is honored only on
        ! the pcg backend, where it enables the Q_NU evidence-bank shell
        ! extension (pcg_priors.md Stage 6.6); the gridding-path challenger
        ! remains unsupported here, and an explicit pcg_nu_lambda_rel=0 would
        ! disable the replay the extension rides on
        if( cline%defined('nu_refine') .and. cline%get_carg('nu_refine').eq.'yes' )then
            if( .not. (cline%defined('rec_backend') .and. cline%get_carg('rec_backend').eq.'pcg') )&
                &THROW_HARD('nu_refine=yes in abinitio3D requires rec_backend=pcg')
            if( cline%defined('pcg_nu_lambda_rel') )then
                if( .not. (cline%get_rarg('pcg_nu_lambda_rel') > 0.) )&
                    &THROW_HARD('nu_refine=yes in abinitio3D requires the Q_NU replay; pcg_nu_lambda_rel=0 disables it')
            endif
        else
            call cline%set('nu_refine', 'no')
        endif
        if( .not. cline%defined('mkdir')       ) call cline%set('mkdir',                    'yes')
        if( .not. cline%defined('overlap')     ) call cline%set('overlap',                   0.95)
        if( .not. cline%defined('prob_athres') ) call cline%set('prob_athres',                10.)
        if( .not. cline%defined('center')      ) call cline%set('center',                    'no')
        if( .not. cline%defined('cenlp')       ) call cline%set('cenlp', abinitio_cenlp_default())
        call cline%set('oritype', 'ptcl3D')
        if( .not. cline%defined('pgrp')        ) call cline%set('pgrp',                      'c1')
        if( .not. cline%defined('pgrp_start')  ) call cline%set('pgrp_start',                'c1')
        if( .not. cline%defined('filt_mode')   ) call cline%set('filt_mode',         'nonuniform')
        ! the pcg-backend automasking veto must be read BEFORE the default
        ! below is injected into the cline; after injection every run looks
        ! like an explicit automsk=no and the stage policy can never engage
        l_automsk_off = (cline%defined('automsk') .and. cline%get_carg('automsk') .eq. 'no')
        if( .not. cline%defined('automsk')     ) call cline%set('automsk',                   'no')
        if( .not. cline%defined('gauref')      ) call cline%set('gauref',                   'yes')
        if( .not. cline%defined('partition')   ) call cline%set('partition',                 'no')
        if( .not. cline%defined('envfsc')      ) call cline%set('envfsc',                    'no')
        if( .not. cline%defined('envmsklp')    ) call cline%set('envmsklp',      ENVMSKLP_DEFAULT)
        if( .not. cline%defined('binwidth')    ) call cline%set('binwidth',   ENVMSKWIDTH_DEFAULT)
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
        if( trim(params%multivol_mode).eq.'docked' .and. last_stage < split_stage )then
            THROW_HARD('multivol_mode=docked requires nstages >= split_stage unless running an explicit pre-split diagnostic')
        endif
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
            start_stage = abinitio_symsrch_stage() + 1
            l_ini3D     = .true.
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
        l_automsk     = (cline%defined('automsk') .and. trim(params%automsk).ne.'no')
        ! l_automsk_off (the EXPLICIT automsk=no veto of the pcg-backend
        ! automasking default) is set where the workflow defaults are
        ! injected, BEFORE automsk=no lands on the cline as a default
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
                if( trim(params%partition).eq.'yes' )then
                    if( .not. spproj%os_cls2D%isthere('cluster') )then
                        THROW_HARD('Missing CLUSTER metadata in CLS2D field needed for PARTITION=YES')
                    endif
                    cls_states = nint(spproj%os_cls2D%get_all('state'))
                    tmpinds    = nint(spproj%os_cls2D%get_all('cluster'))
                    where( cls_states == 0 ) tmpinds = 0
                    clsinds = (/(icls,icls=1,maxval(tmpinds))/)
                    do icls = 1,size(clsinds)
                        if(count(tmpinds==icls) == 0) clsinds(icls) = 0
                    enddo
                    clsinds = pack(clsinds, mask=clsinds>0)
                    call spproj%os_ptcl2D%get_class_sample_stats(clsinds, clssmp, label='cluster')
                    deallocate(cls_states,tmpinds)
                else
                    clsinds = spproj%get_selected_clsinds()
                    call spproj%os_ptcl2D%get_class_sample_stats(clsinds, clssmp)
                endif
                call write_class_samples(clssmp, string(CLASS_SAMPLING_FILE))
                deallocate(clsinds)
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
                call gen_labelling(spproj%os_ptcl3D, params%nstates, 'uniform')
            endif
            call spproj%write_segment_inside(params%oritype, params%projfile)
            if( l_vol_ini_ext )then
                ! user provided input volumes
                call normalize_input_volumes(params, cline_refine3D)
            else
                ! create noise starting volume(s)
                call generate_random_volumes(params, abinitio_stage_box_crop(params, 1), &
                    &abinitio_stage_smpd_crop(params, 1), cline_refine3D)
            endif
        else
            ! check that ptcl3D field is not virgin
            if( spproj%is_virgin_field('ptcl3D') )then
                THROW_HARD('Prior 3D alignment is lacking for starting volume generation')
            endif
            ! randomize states
            if( trim(params%multivol_mode).eq.'independent' .and. .not.l_cavg_ini_ext )then
                call gen_labelling(spproj%os_ptcl3D, params%nstates, 'uniform')
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
            ! This reconstruction feeds start_stage but runs before the stage
            ! loop. Emit the same controller policy here so a requested PCG
            ! backend never executes with the unprocessed outer command line.
            call set_cline_refine3D(params, start_stage, l_cavgs=.false.)
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
            if( params%multivol_mode.eq.'docked' )then
                if( istage == split_stage-1 )then
                    ! update pre-split sampling
                    if( l_force_full_sampling )then
                        update_frac = 1.0
                    else
                        update_frac = real(nstates_glob * params%nsample) / real(nptcls_eff)
                        update_frac = min(abinitio_update_frac_max(), update_frac)
                    endif
                    write(logfhandle,'(A,I0,A,F8.4)') &
                        &'>>> ABINITIO3D DOCKED SPLIT STAGE/PRE-SPLIT_UPDATE_FRAC: ',istage, '/',update_frac
                else if( istage == split_stage )then
                    if( l_force_full_sampling )then
                        update_frac = 1.0
                        call prepare_docked_particle_cohort(nptcls_eff, update_frac)
                    else
                        call calc_docked_multistate_max_sampling(params, nptcls_eff, nptcls_cap, update_frac)
                        call prepare_docked_particle_cohort(nptcls_cap, update_frac)
                        update_frac = real(nstates_glob * params%nsample) / real(nptcls_eff)
                        update_frac = min(abinitio_update_frac_max(), update_frac)
                    endif
                    params%nstates = nstates_glob
                    write(logfhandle,'(A,I0,A,I0,A,F8.4)') &
                        &'>>> ABINITIO3D DOCKED SPLIT STAGE/NSTATES/POSTSPLIT_UPDATE_FRAC: ', &
                        &istage, '/', params%nstates, '/', update_frac
                endif
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
                call randomize_states(params, spproj, params%projfile, xrec3D, split_stage,&
                    &clean_sampling=.false., reconstruct_states=.false.)
                if( l_force_full_sampling )then
                    call calc_rec(params, params%projfile, xrec3D, split_stage)
                else
                    call select_docked_split_reconstruction_sample
                    call calc_rec(params, params%projfile, xrec3D, split_stage, current_sample_only=.true.)
                endif
            endif
            if( cline_refine3D%get_iarg('box_crop') < params%box )then
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',&
                    &cline_refine3D%get_iarg('box_crop')
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
            call strip_pcg_backend_keys(cline_selection)
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

        subroutine prepare_docked_particle_cohort( nsample, ufrac )
            integer, intent(in) :: nsample
            real ,   intent(in) ::ufrac
            integer,  parameter :: DOCKED_NITERS_MISSING = 1
            type(cmdline)       :: cline_missing
            integer             :: nupdated, nmissing, iter_missing, nactive
            call spproj%read_segment('ptcl3D', params%projfile)
            call clean_ptcl3D_sampling
            call spproj%write_segment_inside(params%oritype, params%projfile)
            iter_missing = next_refine3D_iteration()
            write(logfhandle,'(A,A,I0,A,I0,A,I0)') &
                &'>>> ABINITIO3D DOCKED PRE-SPLIT COHORT ASSIGNMENT', &
                &' TARGET/ACTIVE/ITER: ', nsample, '/', nptcls_eff, '/', iter_missing
            cline_missing = cline_refine3D
            call cline_missing%set('prg',               'refine3D')
            call cline_missing%set('mkdir',                   'no')
            call cline_missing%set('refine',                'prob')
            call cline_missing%set('balance',                'yes')
            call cline_missing%set('nsample',              nsample)
            call cline_missing%set('frac_best',                1.0)
            call cline_missing%set('fillin',                  'no')
            call cline_missing%set('update_frac',            ufrac)
            call cline_missing%set('trail_rec',               'no')
            call cline_missing%set('volrec',                  'no')
            call cline_missing%set('sticky_class_sampling',   'no')
            call cline_missing%set('maxits', DOCKED_NITERS_MISSING)
            call cline_missing%set('startit',         iter_missing)
            call cline_missing%set('which_iter',      iter_missing)
            call cline_missing%set('extr_iter',       iter_missing)
            call cline_missing%delete('endit')
            call cline_missing%delete('greedy_sampling')
            call xrefine3D%execute(cline_missing)
            call del_files(DIST_FBODY,      params%nparts, ext='.dat')
            call del_files(ASSIGNMENT_FBODY,params%nparts, ext='.dat')
            call del_file(DIST_FBODY//'.dat')
            call del_file(ASSIGNMENT_FBODY//'.dat')
            call read_multistate_assignment_coverage(nactive, nupdated, nmissing)
            if( nupdated < nsample )then
                THROW_HARD('docked pre-split cohort assignment updated too few particles')
            endif
            call cline_missing%kill
        end subroutine prepare_docked_particle_cohort

        subroutine select_docked_split_reconstruction_sample
            type(class_sample), allocatable :: split_clssmp(:)
            integer, allocatable :: split_pinds(:), sampled(:), states(:)
            integer :: noris, nrequested, nselected, sample_ind, split_state, state_pop
            call spproj%read_segment('ptcl3D', params%projfile)
            noris   = spproj%os_ptcl3D%get_noris()
            if( nptcls_eff < 1 ) THROW_HARD('docked split reconstruction has no active particles')
            call read_class_samples(split_clssmp, string(CLASS_SAMPLING_FILE))
            call spproj%os_ptcl3D%sample4update_class(split_clssmp, [1,noris], update_frac,&
                &nselected, split_pinds, .true., .false., sampled_only=.true.)
            nrequested = nint(update_frac * real(nptcls_eff))
            if( nselected /= nrequested )then
                THROW_HARD('docked split reconstruction failed to select the exact requested sample')
            endif
            call spproj%write_segment_inside(params%oritype, params%projfile)
            sampled   = spproj%os_ptcl3D%get_all_asint('sampled')
            states    = spproj%os_ptcl3D%get_all_asint('state')
            sample_ind = maxval(sampled)
            do split_state = 1, params%nstates
                state_pop = count(states == split_state .and. sampled == sample_ind)
                if( state_pop <= 5 )then
                    THROW_HARD('docked split reconstruction sample has insufficient state population')
                endif
            enddo
            write(logfhandle,'(A,I0,A,I0,A,F8.4)') &
                &'>>> ABINITIO3D DOCKED SPLIT RECONSTRUCTION SAMPLE SELECTED/COHORT/FRACTION: ',&
                &nselected, '/', count(sampled > 0), '/', real(nselected) / real(count(sampled > 0))
            if( allocated(sampled)    ) deallocate(sampled)
            if( allocated(states)     ) deallocate(states)
            if( allocated(split_pinds)) deallocate(split_pinds)
            if( allocated(split_clssmp) ) call deallocate_class_samples(split_clssmp)
        end subroutine select_docked_split_reconstruction_sample

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
            call cline_missing%set('frac_best',              1.0)
            call cline_missing%set('fillin',                'no')
            call cline_missing%set('update_missing',       'yes')
            call cline_missing%set('update_frac',            1.0)
            call cline_missing%set('trail_rec',             'no')
            call cline_missing%set('volrec',                'no')
            call cline_missing%set('sticky_class_sampling', 'no')
            call cline_missing%set('maxits',                   1)
            call cline_missing%set('startit',       iter_missing)
            call cline_missing%set('which_iter',    iter_missing)
            call cline_missing%set('extr_iter',     iter_missing)
            call cline_missing%delete('endit')
            call cline_missing%delete('partition')
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
            ! Particle-stage PCG policy belongs to the outer abinitio3D run;
            ! class-average initialization retains its gridding workflow
            call strip_pcg_backend_keys(cline_ini3D)
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
            call cline_ini3D%delete('projrec') ! compact projection sums are for particle refinement stages
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
