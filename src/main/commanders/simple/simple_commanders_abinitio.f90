!@descr: abinitio 3D reconstruction in single- and multi-particle mode
module simple_commanders_abinitio
use simple_commanders_api
use simple_abinitio_utils
use simple_abinitio3D_split_checkpoint,            only: build_abinitio3D_split_checkpoint
use simple_final_rec,                              only: calc_final_rec
use simple_external_reference_pose_initialization, only: initialize_poses_against_external_references
use simple_procimgstk,                             only: shift_imgfile
use simple_commanders_project_core,                only: commander_selection
use simple_commanders_reproject,                   only: commander_reproject
use simple_commanders_refine3D,                    only: commander_refine3D, commander_refine3D_states, commander_bootstrap_rec3D
use simple_commanders_rec,                         only: commander_rec3D
use simple_cluster_seed,                           only: gen_labelling
use simple_refine3D_fnames,                        only: refine3D_startvol_fname, refine3D_startvol_half_fname, &
    &refine3D_state_vol_fname, refine3D_state_halfvol_fname, refine3D_frozen_context_fname, refine3D_fsc_fname
use simple_halfmap_diagnostics,                    only: rename_support_provenance
use simple_gui_communicator,                       only: gui_communicator
use simple_abinitio3D_manifest,                    only: abinitio3D_manifest, abinitio3D_stage_record, MANIFEST_FNAME, &
    &manifest_records_input
use simple_project_superset,                       only: project_superset, COHORT_WARN_FRAC
use simple_frozen_accum,                           only: frozen_accum
use simple_abinitio3D_addon_report,                only: abinitio3D_addon_report, ADDON_REPORT_FNAME
implicit none

public :: commander_abinitio3D_cavgs, commander_abinitio3D_cavgs_conditional_restarts
public :: commander_abinitio3D, commander_abinitio3D_addon
private
#include "simple_local_flags.inc"

character(len=*), parameter :: CAVGS_SIGMA2_STATE_FNAME = 'abinitio3D_cavgs_sigma2_state.bin'

type, extends(commander_base) :: commander_abinitio3D_cavgs
    contains
    procedure :: execute => exec_abinitio3D_cavgs
end type commander_abinitio3D_cavgs

type, extends(commander_base) :: commander_abinitio3D_cavgs_conditional_restarts
    contains
    procedure :: execute => exec_abinitio3D_cavgs_conditional_restarts
end type commander_abinitio3D_cavgs_conditional_restarts

type, extends(commander_base) :: commander_abinitio3D
    contains
    procedure :: execute => exec_abinitio3D
end type commander_abinitio3D

!> abinitio3D_addon: a thin wrapper that validates both projects and the base
!! run's manifest before any write, turns the manifest and the add-on's own
!! inputs into a fresh command line and enters exec_abinitio3D through the
!! internal addon_manifest handshake, which no command line can set
type, extends(commander_base) :: commander_abinitio3D_addon
    contains
    procedure :: execute => exec_abinitio3D_addon
end type commander_abinitio3D_addon


contains

    !> for generation of an initial 3D model from class averages
    subroutine exec_abinitio3D_cavgs( self, cline )
        use simple_estimate_ssnr, only: lpstages_fast
        class(commander_abinitio3D_cavgs), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        ! shared-mem commanders
        type(commander_refine3D)  :: xrefine3D
        type(commander_rec3D)     :: xrec3D
        type(commander_bootstrap_rec3D) :: xbootstrap_rec3D
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
        integer                   :: nstates_on_cline, nstates_target, split_stage, pop
        integer                   :: cavg_ldim(3), cavg_nimgs, final_nstates
        real                      :: cavg_smpd
        if( cline%defined('part') )then
            THROW_HARD('abinitio3D_cavgs distributed execution is master-only; remove part from command line')
        endif
        ! the parser accepts every vocabulary key for every program
        if( cline%defined('projfile_frozen') .or. cline%defined('addon_diag') )then
            THROW_HARD('projfile_frozen and addon_diag belong to abinitio3D_addon, not abinitio3D_cavgs')
        endif
        l_state_continue_mode = .false.
        call cline%set('sigma_est', 'global') ! obviously
        call cline%set('oritype',      'out') ! because cavgs are part of out segment
        call cline%set('bfac',            0.) ! because initial models should not be sharpened
        call cline%set('filt_mode',   'none') ! no fancy filtering for cavgs route
        call cline%set('automsk',       'no') ! no envelope masking for cavgs route
        call cline%set('objfun', 'euclid') ! noise normalized Euclidean distances from the start
        if( .not. cline%defined('mkdir')            ) call cline%set('mkdir',                      'yes')
        if( .not. cline%defined('overlap')          ) call cline%set('overlap',                     0.95)
        if( .not. cline%defined('prob_athres')      ) call cline%set('prob_athres',                  90.) ! reduces # failed runs on trpv1 from 4->2/10
        if( .not. cline%defined('cenlp')            ) call cline%set('cenlp',   abinitio_cenlp_default())
        if( .not. cline%defined('imgkind')          ) call cline%set('imgkind',                   'cavg')
        if( .not. cline%defined('filt_mode')        ) call cline%set('filt_mode',                 'none')
        if( .not. cline%defined('noise_norm')       ) call cline%set('noise_norm',                  'no')
        if( .not. cline%defined('lpstart')          ) call cline%set('lpstart', abinitio_lpstart_ini3D())
        if( .not. cline%defined('lpstop')           ) call cline%set('lpstop',   abinitio_lpstop_ini3D())
        if( .not. cline%defined('gauref')           ) call cline%set('gauref',                     'yes')
        if( .not. cline%defined('exit_collapse')    ) call cline%set('exit_collapse',               'no')
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
                if( allocated(lpinfo) ) deallocate(lpinfo)
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
        ! The temporary particle lineage must never inherit the input
        ! project's canonical state. Its class-average rows have a different
        ! identity and the temporary project is deleted when this workflow
        ! completes.
        call work_proj%projinfo%delete_entry('sigma2_state')
        ! name change
        call work_proj%projinfo%delete_entry('projname')
        call work_proj%projinfo%delete_entry('projfile')
        call cline%set('projfile', work_projfile)
        call cline%set('projname', get_fbody(work_projfile,'simple'))
        call work_proj%update_projinfo(cline)
        call del_file(CAVGS_SIGMA2_STATE_FNAME)
        call work_proj%set_sigma2_state_path(string(CAVGS_SIGMA2_STATE_FNAME))
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
            write(logfhandle,'(A,I3,A,F5.1,A)')'>>> STAGE ', istage,' WITH LP ', lpinfo(istage)%lp, ' A'
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
            ! Early exit on state collapse
            if( nstates_target > 1 .and. params%exit_collapse .eq. 'yes' ) then
                write(logfhandle,'(A,A)')'>>> CHECKING FOR STATE COLLAPSE...', work_projfile%to_char()
                call work_proj%read_segment('ptcl3D', work_projfile)
                states = nint(work_proj%os_ptcl3D%get_all('state'))
                final_nstates = 0
                do s = 1, nstates_target
                    pop = count(states == s)
                    if( pop > 0 ) final_nstates = final_nstates + 1
                    write(logfhandle, '(A,I0,A,I0)') '>>> FINAL POPULATION STATE ', s, ': ', pop
                enddo
                if( final_nstates < 2 )then
                    write(logfhandle,'(A)')'>>> EARLY EXIT DUE TO STATE COLLAPSE'
                    exit
                endif
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
            call calc_final_rec(params, work_proj, work_projfile, cline_refine3D, xrec3D, xbootstrap_rec3D, &
                &l_postprocess=.false., lp_snapshot=lpinfo(nstages_ini3D)%lp)
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
        ! Message for conditional restarts
        final_nstates = 1
        if( nstates_target > 1 )then
            states = nint(spproj%os_cls3D%get_all('state'))
            final_nstates = 0
            do s = 1, nstates_target
                pop = count(states == s)
                if( pop > 0 ) final_nstates = final_nstates + 1
                write(logfhandle, '(A,I0,A,I0)') '>>> FINAL POPULATION STATE ', s, ': ', pop
            enddo
        endif
        call cline%set('final_nstates', final_nstates)
        ! remove postprocessed (pproc) volumes; with bfac=0 they add nothing in cavgs mode
        call del_pproc_vols
        ! end gracefully
        call img%kill
        call spproj%kill
        call o%kill
        call o_even%kill
        call o_odd%kill
        call work_proj%kill
        call del_file(CAVGS_SIGMA2_STATE_FNAME)
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

    subroutine exec_abinitio3D_cavgs_conditional_restarts( self, cline )
        class(commander_abinitio3D_cavgs_conditional_restarts), intent(inout) :: self
        class(cmdline),                                         intent(inout) :: cline
        type(commander_abinitio3D_cavgs) :: xcommander_abinitio3D_cavgs
        type(cmdline) :: cline_backup
        integer       :: nrestarts, irestart, final_nstates, input_nstates, nstates_collapse
        logical       :: state_collapse, l_mkdir
        if( .not. cline%defined('nrestarts_collapse') )then
            THROW_HARD('nrestarts_collapse needs to be on the command line for abinitio3D_cavgs state collapse conditional restarts')
        endif
        if( cline%defined('nrestarts') )then
            THROW_HARD('nrestarts is not compatible with abinitio3D_cavgs state collapse conditional restarts')
        endif
        if( .not. cline%defined('nstates') )then
            THROW_HARD('nstates needs to be defined on command line for abinitio3D_cavgs state collapse conditional restarts')
        endif
        input_nstates  = cline%get_iarg('nstates')
        if( input_nstates == 1 )then
            THROW_HARD('nstates needs to be greater than 1 for abinitio3D_cavgs state collapse conditional restarts')
        endif
        if( .not.cline%defined('projfile') )then
            THROW_HARD('projfile needs to be defined on command line for abinitio3D_cavgs state collapse conditional restarts')
        endif
        if( cline%defined('mkdir') )then
            l_mkdir = cline%get_carg('mkdir')=='yes'
        else
            call cline%set('mkdir', 'yes')
            l_mkdir = .true.
        endif
        if( .not.l_mkdir ) THROW_HARD('MKDIR must be YES for abinitio3D_cavgs state collapse conditional restarts')
        cline_backup     = cline
        nrestarts        = cline%get_iarg('nrestarts_collapse')
        nstates_collapse = 0
        write(logfhandle,'(A,I0,A,I0)') '>>> INITIAL RUN WITH NSTATES=', input_nstates
        write(logfhandle,'(A)') '>>>'
        do irestart = 1, nrestarts
            cline = cline_backup
            call cline%delete('nrestarts_collapse')
            call cline%set('exit_collapse', 'yes')
            if (irestart == nrestarts ) call cline%delete('exit_collapse')
            call xcommander_abinitio3D_cavgs%execute(cline)
            if( l_mkdir ) call chdir('..')
            final_nstates = cline%get_iarg('final_nstates')
            call cline%delete('final_nstates')
            state_collapse = (final_nstates == 1)
            if( state_collapse )then
                nstates_collapse = nstates_collapse + 1
                write(logfhandle,'(A,I0)') '>>>'
                write(logfhandle,'(A,I0,A,I0)') '>>> STATE COLLAPSE DURING RUN ', irestart,&
                    &' - FINAL NSTATES=', final_nstates
                write(logfhandle,'(A)')    '>>>'
            else
                exit
            endif
        end do
        ! cleanup
        call cline_backup%kill
        call simple_touch(TASK_FINISHED)
    end subroutine exec_abinitio3D_cavgs_conditional_restarts

    !> Validate, translate and hand over: nothing is written before both
    !! projects, their shared particle index space and the base run's manifest
    !! have been validated.
    subroutine exec_abinitio3D_addon( self, cline )
        use simple_ui,         only: get_prg_ptr
        use simple_ui_program, only: ui_program
        class(commander_abinitio3D_addon), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(commander_abinitio3D)    :: xabinitio3D
        type(abinitio3D_manifest)     :: man
        type(abinitio3D_stage_record) :: stage
        type(sp_project)              :: spproj_cur, spproj_frz
        type(project_superset)        :: superset
        type(ui_program), pointer     :: ui_addon => null()
        type(cmdline)                 :: cline_run
        type(string), allocatable     :: keys(:)
        type(string)                  :: projfile, projfile_frz, sigma_path
        character(len=STDLEN)         :: msg
        character(len=:), allocatable :: key
        integer :: i, status
        logical :: found
        ! The command line is the program's UI contract: its declared inputs and
        ! the execution environment. A setting of the frozen solution comes
        ! from its manifest; refusal is by key, so an inherited value cannot be
        ! confirmed silently.
        call get_prg_ptr(string('abinitio3D_addon'), ui_addon)
        if( .not. associated(ui_addon) ) THROW_HARD('the abinitio3D_addon user interface is not registered')
        keys = cline%get_keys()
        do i = 1, size(keys)
            key = trim(keys(i)%to_char())
            select case(key)
                case('prg', 'projfile', 'mkdir')
                    cycle
            end select
            if( ui_addon%accepts(key) ) cycle
            if( manifest_records_input(key) )then
                THROW_HARD(key//' is set by the frozen run (recorded in its manifest); remove it from the command line')
            else
                THROW_HARD(key//' is not an input of abinitio3D_addon')
            endif
        enddo
        if( .not. cline%defined('projfile') )        THROW_HARD('abinitio3D_addon requires projfile')
        if( .not. cline%defined('projfile_frozen') ) THROW_HARD('abinitio3D_addon requires projfile_frozen')
        ! both project paths are normalised before any change of directory
        projfile     = simple_abspath(cline%get_carg('projfile'))
        projfile_frz = simple_abspath(cline%get_carg('projfile_frozen'))
        if( projfile%to_char() == projfile_frz%to_char() )then
            THROW_HARD('projfile and projfile_frozen are the same file (or aliases of it)')
        endif
        call spproj_frz%read(projfile_frz)
        call spproj_cur%read(projfile)
        ! the manifest is the only route into the add-on
        call man%read_registered(spproj_frz, projfile_frz, status, msg)
        if( status /= 0 ) THROW_HARD('projfile_frozen: '//trim(msg))
        call man%validate_frozen(spproj_frz, status, msg)
        if( status /= 0 ) THROW_HARD(trim(msg))
        call man%get_artifact('sigma2_state', 0, sigma_path, found)
        if( .not. found ) THROW_HARD('the frozen run recorded no committed residual sigma2 state')
        if( .not. man%matches_artifact('sigma2_state', 0, sigma_path) )then
            THROW_HARD('the frozen run''s sigma2 state is missing or changed: '//sigma_path%to_char())
        endif
        if( spproj_cur%get_box() /= man%get_box() .or. abs(spproj_cur%get_smpd() - man%get_smpd()) > 1.e-4*man%get_smpd() )then
            THROW_HARD('the current project''s native box or sampling differs from the frozen solution''s')
        endif
        ! every consuming box gets a frozen set from reconstruct3D, which never
        ! upsamples: a ladder whose stage box exceeds the native box (small
        ! boxes, where the crop rounds up to a larger magic box) is refused
        do i = abinitio_symsrch_stage(), man%get_last_stage()
            stage = man%get_stage(i)
            if( stage%box_crop > man%get_box() )then
                THROW_HARD('the base run''s ladder upsamples (stage '//int2str(i)//' box '//int2str(stage%box_crop)//' > native '//int2str(man%get_box())//'); unsupported by abinitio3D_addon')
            endif
        enddo
        ! one particle index space, the superset relation, the membership
        call superset%new(spproj_cur, spproj_frz, man%get_nstates(), man%get_ptcl_src(), status, msg)
        if( status /= 0 ) THROW_HARD(trim(msg))
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> ABINITIO3D_ADDON FROZEN/COHORT/NEVER-UPDATED/STATES: ', &
            &superset%get_nfrozen(), '/', superset%get_ncohort(), '/', superset%get_nnever_updated(), '/', man%get_nstates()
        if( superset%is_small_cohort() ) THROW_WARN('the cohort is below '//int2str(nint(100.*COHORT_WARN_FRAC))//'% of the frozen population')
        call superset%kill
        call spproj_cur%kill
        call spproj_frz%kill
        ! a fresh, sparse command line from the allowlisted manifest records
        call cline_run%set('prg',             'abinitio3D_addon')
        call cline_run%set('projfile',        projfile)
        call cline_run%set('projfile_frozen', projfile_frz)
        call cline_run%set('mkdir',           'yes')
        if( cline%defined('mkdir') ) call cline_run%set('mkdir', cline%get_carg('mkdir'))
        ! the base run's replayed inputs, stage-line shape and solution
        call man%replay(cline_run)
        call cline_run%set('center', 'no')  ! a centring shift would never reach the frozen term
        ! the add-on's own inputs and the execution environment, as given
        do i = 1, size(keys)
            key = trim(keys(i)%to_char())
            select case(key)
                case('prg', 'projfile', 'projfile_frozen', 'mkdir')
                    cycle
            end select
            call cline_run%copy_arg(cline, key)
        enddo
        call cline_run%set('addon_manifest', man%get_fname())
        call man%kill
        call xabinitio3D%execute(cline_run)
        ! all done: the finished project replaces the original project file,
        ! and simple_exec's job record goes to it
        call publish_result(cline_run%get_carg('projfile'))
        call cline%set('projfile', projfile)
        call cline_run%kill
        call simple_touch(TASK_FINISHED)

    contains

        !> The finished working project (frozen rows restored, the cohort
        !! refined) replaces the original project file: written beside it under
        !! a temporary name, then renamed over it, so a failed run never touches
        !! the original. The run manifest is registered by absolute path, since
        !! a bare name resolves against the project file's own directory. Under
        !! mkdir=no the working project is the original, already in place.
        subroutine publish_result( working )
            class(string), intent(in) :: working
            type(abinitio3D_manifest) :: man_run
            type(sp_project)          :: spproj_out
            type(string)              :: working_abs, tmpname, manifest_abs
            character(len=STDLEN)     :: msg_here
            integer :: status_here
            working_abs = simple_abspath(working)
            if( working_abs%to_char() == projfile%to_char() ) return
            call spproj_out%read(working_abs)
            call man_run%read_registered(spproj_out, working_abs, status_here, msg_here)
            if( status_here == 0 )then
                manifest_abs = man_run%get_fname()
                call man_run%register(spproj_out, manifest_abs%to_char())
            else
                THROW_WARN('abinitio3D_addon: the published project registers no run manifest: '//trim(msg_here))
            endif
            call spproj_out%projinfo%set(1, 'projfile', projfile%to_char())
            tmpname = stemname(projfile)//'/abinitio3D_addon_publish_'//basename(projfile)
            call spproj_out%write(tmpname)
            call simple_rename(tmpname, projfile, overwrite=.true.)
            write(logfhandle,'(A,A)') '>>> ABINITIO3D_ADDON: PUBLISHED THE FINISHED PROJECT TO ', projfile%to_char()
            call man_run%kill
            call spproj_out%kill
            call working_abs%kill
            call tmpname%kill
            call manifest_abs%kill
        end subroutine publish_result

    end subroutine exec_abinitio3D_addon

    !> for generation of an initial 3d model from particles
    subroutine exec_abinitio3D( self, cline )
        class(commander_abinitio3D), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        ! commanders
        type(commander_refine3D)        :: xrefine3D
        type(commander_refine3D_states) :: xrefine3D_states
        type(commander_rec3D)           :: xrec3D
        type(commander_bootstrap_rec3D) :: xbootstrap_rec3D
        ! other
        integer,            allocatable :: tmpinds(:), clsinds(:), pinds(:), cls_states(:)
        type(class_sample), allocatable :: clssmp(:)
        type(string),       allocatable :: external_refs(:), external_checkpoint(:)
        type(parameters)                :: params
        type(sp_project)                :: spproj
        type(gui_communicator)          :: gui_comm
        real    :: lprange(2)
        integer :: state, istage, icls, start_stage, nptcls2update, noris, nstates_on_cline
        integer :: nstates_in_project, split_stage, last_stage, pose_init_iter
        logical :: l_cavg_ini_ext, l_vol_ini_ext, l_user_nstages, l_user_lpstop, l_run_final_rec
        logical :: l_state_continue
        logical :: l_force_full_sampling
        logical :: l_states_handoff_complete
        real    :: sampled_active_frac
        real    :: update_frac_post_split
        ! run manifest: the command line as given and the emitted stage ladder
        type(cmdline)     :: cline_entry
        real, allocatable :: emitted_lp(:), emitted_lpstop(:)
        character(len=64) :: run_id
        ! abinitio3D_addon route, selected only by the internal addon_manifest
        ! handshake that commander_abinitio3D_addon sets (outside the argument
        ! vocabulary, so no command line can select it)
        type(abinitio3D_manifest)  :: man_addon
        type(project_superset)     :: superset
        type(sp_project)           :: spproj_frz
        type(abinitio3D_addon_ctx) :: addon_ctx
        type(string)               :: frozen_copy, frozen_ctx_fname, frozen_sigma_copy
        logical :: l_addon, l_cohort_chain_seeded
        character(len=*), parameter :: ADDON_DIAG_DIR  = 'addon_diag' !< addon_diag=yes: the cohort-only map
        character(len=*), parameter :: FROZEN_COPY_DIR = 'frozen'     !< the frozen project's copy
        cline_entry = cline
        l_addon = cline%defined('addon_manifest')
        l_cohort_chain_seeded = .false.
        if( l_addon )then
            run_id = new_run_id('abinitio3D_addon')
        else
            run_id = new_run_id('abinitio3D')
            if( cline%defined('projfile_frozen') .or. cline%defined('addon_diag') )then
                THROW_HARD('projfile_frozen and addon_diag belong to abinitio3D_addon, not abinitio3D')
            endif
        endif
        l_state_continue = cline%defined('state')
        l_force_full_sampling = .false.
        l_states_handoff_complete = .false.
        sampled_active_frac   = 0.
        l_state_continue_mode = .false.
        ! Particle caching is a 2D-only feature; reject rather than silently ignore
        if( cline%defined('cache') .or. cline%defined('cache_dir') )then
            THROW_HARD('cache=yes (downscaled particle cache) is not supported for 3D workflows')
        endif
        call cline%set('objfun',    'euclid') ! use noise normalized Euclidean distances from the start
        call cline%set('sigma_est', 'global') ! obviously
        call cline%set('bfac',            0.) ! because initial models should not be sharpened
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
        if( cline%defined('multivol_mode') .and. .not. l_addon )then
            if( cline%get_carg('multivol_mode').eq.'independent' )then
                ! Stop independent multi-state starts before prob_neigh/NU by default.
                if( .not. l_user_nstages ) call cline%set('nstages', abinitio_independent_nstages_default())
                if( .not. l_user_lpstop  ) call cline%set('lpstop',  abinitio_independent_lpstop_default())
            endif
        endif
        ! make master parameters
        call params%new(cline)
        call gui_comm%new(params)
        if( l_addon ) call addon_parse
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
        ! read project
        call spproj%read(params%projfile)
        ! A fresh abinitio3D never continues another run's sigma2 estimate: a
        ! canonical registration inherited with the project is dropped so the
        ! first euclid stage seeds in this run's own directory (2026-09-07)
        if( spproj%projinfo%get_noris() == 1 )then
            if( spproj%projinfo%isthere(1, 'sigma2_state') )then
                call spproj%projinfo%delete_entry('sigma2_state')
                call spproj%write_segment_inside('projinfo', params%projfile)
                write(logfhandle,'(A)') '>>> ABINITIO3D: dropped an inherited canonical sigma2 registration; sigmas are seeded here'
            endif
        endif
        ! add-on prologue: frozen copy, physical identity, frozen-row mask; the
        ! mask precedes every count below, so the established sampling
        ! initialisation sees the cohort alone
        if( l_addon ) call addon_prologue
        ! provide initialization of 3D alignment using class averages?
        start_stage = 1
        l_ini3D     = .false.
        ! abinitio3D_addon has one entry route: stage 3 (prob) with trusted
        ! frozen references, no symmetry search (pgrp_start = pgrp) and no CC
        ! pose initialisation
        if( l_addon ) start_stage = abinitio_symsrch_stage()
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
            ! execution
            call ini3D_from_cavgs(cline)
            ! re-read the project file to update info in spproj
            call spproj%read(params%projfile)
            start_stage = abinitio_nstages_ini3D() - 1 ! compute reduced to two overlapping stages
            l_ini3D     = .true.
            ! symmetry dealt with by ini3D
        endif
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
            ! CC pose initialization below owns the external-reference route.
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
        if( l_addon )then
            ! the base run's ladder at the limits it emitted; never planned from class FRCs
            call set_lplims_from_manifest(man_addon)
        else if( l_vol_ini_ext )then
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
        if( l_user_lpstop ) write(logfhandle,'(A,F8.3,A)') &
            &'>>> ABINITIO3D COMMAND-LINE LPSTOP CEILING: ', params%lpstop, ' A'
        ! the limits the controller actually emitted, per stage the loop runs
        ! (0 = not on the stage line); -1 marks a stage the run never ran
        allocate(emitted_lp(size(lpinfo)), emitted_lpstop(size(lpinfo)))
        emitted_lp     = -1.
        emitted_lpstop = -1.
        ! starting volume logics
        if( l_addon )then
            ! random cohort poses and labels, the per-box frozen sets and the
            ! native frozen references as the stage-3 vol1..volN
            call addon_starting_state
        else if( .not. l_ini3D )then
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
                call set_cline_refine3D(params, start_stage, l_cavgs=.false.)
                pose_init_iter = max(1, cline_refine3D%get_iarg('startit') - 1)
                allocate(external_refs(params%nstates), external_checkpoint(params%nstates))
                external_refs = params%vols(1:params%nstates)
                write(logfhandle,'(A)') &
                    &'>>> ABINITIO3D EXTERNAL REFERENCES ARE UNTRUSTED UNTIL CC POSE INITIALIZATION COMPLETES'
                call initialize_poses_against_external_references(params, cline_refine3D, xrefine3D, xrec3D, &
                    &nptcls_eff, external_refs, external_checkpoint, pose_init_iter)
                params%vols(1:params%nstates) = external_checkpoint
                do state = 1,params%nstates
                    call cline_refine3D%set('vol'//int2str(state), params%vols(state))
                enddo
                call cline_refine3D%set('endit', pose_init_iter)
                deallocate(external_refs, external_checkpoint)
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
        do istage = start_stage, nstages_refine3D
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
                endif
            endif
            ! Preparation of command line for refinement
            if( params%multivol_mode.eq.'docked' .and. istage == split_stage )then
                ! A local receives the intent(out) update fraction: the checkpoint
                ! routine reads the module-level update_frac through
                ! set_cline_refine3D, so passing the module variable itself would
                ! alias an intent(out) dummy with host-associated reads (F2018
                ! 15.5.2.13).
                call build_abinitio3D_split_checkpoint(params, spproj, xrefine3D, xrec3D, split_stage, &
                    &nptcls_eff, nstates_glob, l_force_full_sampling, update_frac_post_split)
                update_frac = update_frac_post_split
            else if( l_addon )then
                call set_cline_refine3D(params, istage, l_cavgs=.false., addon=addon_ctx)
            else
                call set_cline_refine3D(params, istage, l_cavgs=.false.)
            endif
            call record_emitted_limits(istage)
            if( l_addon ) call addon_stage_boundary(istage)
            write(logfhandle,'(A)')'>>>'
            if( cline_refine3D%defined('lp') )then
                if( l_refine3D_lp_override )then
                    write(logfhandle,'(A,I3,A,F5.1,A)')'>>> STAGE ', istage,' WITH LP ', &
                        &cline_refine3D%get_rarg('lp'), ' A (command line)'
                else if( abs(cline_refine3D%get_rarg('lp') - lpinfo(istage)%lp) > 1.e-3 )then
                    write(logfhandle,'(A,I3,A,F5.1,A,F5.1,A)')'>>> STAGE ', istage,' WITH LP ', &
                        &cline_refine3D%get_rarg('lp'), ' A (planned ', lpinfo(istage)%lp, ' A, promoted by FSC=0.5)'
                else
                    write(logfhandle,'(A,I3,A,F5.1,A)')'>>> STAGE ', istage,' WITH LP ', cline_refine3D%get_rarg('lp'), ' A'
                endif
            else
                write(logfhandle,'(A,I3,A)')'>>> STAGE ', istage,' WITH NU-SELECTED MATCHING LP'
            endif
            if( params%multivol_mode.eq.'docked' .and. istage == split_stage )then
                call handoff_split_checkpoint_to_refine3D_states
                l_states_handoff_complete = .true.
                exit
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
            ! update GUI
            call spproj%read_segment('cls3D',  params%projfile)
            call spproj%read_segment('ptcl3D', params%projfile)
            call spproj%read_segment('out',    params%projfile)
            call gen_ortho_reprojs4viz(params, spproj)
            call gui_comm%add_metadata(spproj, oritype='cls3D', stage=istage)
        enddo
        if( l_states_handoff_complete )then
            write(logfhandle,'(A)') &
                &'>>> ABINITIO3D POST-SPLIT REFINEMENT AND FINAL RECONSTRUCTION OWNED BY REFINE3D_STATES'
        else if( l_run_final_rec )then
            select case(trim(params%multivol_mode))
                case('independent','docked')
                    call ensure_multistate_particle_assignments
            end select
            ! the shared ending: final all-particle reconstruction at original
            ! sampling, project registration, final products and reprojections
            ! abinitio3D_addon: every particle is active again and the
            ! cohort-only sigma2 registration is gone, so the final
            ! reconstruction bootstraps the union's own state over every
            ! particle (bootstrap_rec3D) and ships the union map on it
            if( l_addon ) call addon_restore_union
            call calc_final_rec(params, spproj, params%projfile, cline_refine3D, xrec3D, xbootstrap_rec3D, &
                &l_postprocess=.true., lp_snapshot=lpinfo(nstages_refine3D)%lp)
        else
            write(logfhandle,'(A,I0)')'>>> ABINITIO3D EARLY STOP AFTER STAGE ', nstages_refine3D
            write(logfhandle,'(A)')'>>> FINAL ALL-PARTICLE RECONSTRUCTION SKIPPED'
        endif
        ! the run manifest, last and never fatal: a completed run (final
        ! reconstruction or refine3D_states handoff) is a candidate frozen input
        if( l_addon )then
            ! epilogue: union metadata and the report against the base
            ! solution; with its union sigma2 state the output is a frozen
            ! input for the next add-on
            if( .not. l_run_final_rec ) THROW_HARD('abinitio3D_addon inherited a ladder without a final reconstruction')
            call addon_epilogue
            call write_run_manifest(.true., 'abinitio3D_addon')
        else if( l_states_handoff_complete .or. l_run_final_rec )then
            call write_run_manifest(.true., 'abinitio3D')
        endif
        ! final update GUI
        call spproj%read_segment('cls2D',  params%projfile)
        call spproj%read_segment('ptcl3D', params%projfile)
        call spproj%read_segment('out',    params%projfile)
        call gui_comm%add_metadata(spproj, oritype='cls3D', stage=0, selection=.true.) ! stage=0 signifies final
        ! cleanup
        call spproj%kill
        call qsys_cleanup(params)
        call gui_comm%kill()
        call simple_end('**** SIMPLE_ABINITIO3D NORMAL STOP ****', &
            verbose_exit=trim(params%verbose_exit).eq.'yes', verbose_exit_fname=params%verbose_exit_fname)

    contains

        subroutine record_emitted_limits( istage_run )
            integer, intent(in) :: istage_run
            if( istage_run < 1 .or. istage_run > size(emitted_lp) ) return
            emitted_lp(istage_run)     = 0.
            emitted_lpstop(istage_run) = 0.
            if( cline_refine3D%defined('lp') )     emitted_lp(istage_run)     = cline_refine3D%get_rarg('lp')
            if( cline_refine3D%defined('lpstop') ) emitted_lpstop(istage_run) = cline_refine3D%get_rarg('lpstop')
        end subroutine record_emitted_limits

        ! ------------------------------------------------------------------
        ! abinitio3D_addon route
        ! ------------------------------------------------------------------

        !> typed add-on inputs: the manifest named by the handshake, the add-on
        !! context; the add-on keys leave the command line so no child sees them
        subroutine addon_parse
            character(len=STDLEN) :: msg
            type(string) :: mpath
            integer :: status
            mpath = cline%get_carg('addon_manifest')
            call man_addon%read(mpath, status, msg)
            if( status /= 0 ) THROW_HARD(trim(msg))
            if( .not. file_exists(params%projfile_frozen) ) THROW_HARD('abinitio3D_addon: projfile_frozen does not exist')
            if( params%nstates /= man_addon%get_nstates() ) THROW_HARD('abinitio3D_addon: state layout differs from the manifest')
            call cline%delete('projfile_frozen')
            call cline%delete('addon_diag')
            call cline%delete('addon_manifest')
            frozen_ctx_fname   = simple_abspath(refine3D_frozen_context_fname(), check_exists=.false.)
            ! the copy keeps the frozen project's file name, in a directory of its
            ! own: reading a project resets projname to its file name, and
            ! projname is the lineage of the sigma2 layout digest
            frozen_copy        = simple_abspath(string(FROZEN_COPY_DIR//'/')//basename(params%projfile_frozen), &
                &check_exists=.false.)
            frozen_sigma_copy  = simple_abspath(string(FROZEN_COPY_DIR//'/frozen_sigma2_state.bin'), check_exists=.false.)
            addon_ctx%active     = .true.
            addon_ctx%overlap    = params%overlap
            addon_ctx%frozen_rec = frozen_ctx_fname
            write(logfhandle,'(A,A)') '>>> ABINITIO3D_ADDON RUN: ', trim(run_id)
            write(logfhandle,'(A,A)') '>>> ABINITIO3D_ADDON BASE RUN: ', trim(man_addon%get_run_id())
            call mpath%kill
        end subroutine addon_parse

        !> The frozen project is copied, never written: the copy is its own
        !! project file under the frozen project's name in a directory of its
        !! own (reading a project resets projname, the sigma2 layout lineage, to
        !! the file name, and projinfo projfile, the target of every segment write
        !! without a file name, to the file itself, so nothing can reach the
        !! working copy or the frozen run), and it owns a copy of the frozen run's
        !! committed residual sigma2 state, registered by absolute path so that
        !! no working-directory convention enters its resolution. Then the
        !! identity of the run-directory copies is re-validated and the frozen
        !! rows of the working copy are masked.
        subroutine addon_prologue
            character(len=STDLEN) :: msg
            type(string) :: sigma_src
            integer :: status
            logical :: found
            call spproj_frz%read(params%projfile_frozen)
            call superset%new(spproj, spproj_frz, man_addon%get_nstates(), man_addon%get_ptcl_src(), status, msg)
            if( status /= 0 ) THROW_HARD(trim(msg))
            call man_addon%get_artifact('sigma2_state', 0, sigma_src, found)
            if( .not. found ) THROW_HARD('the frozen run recorded no committed residual sigma2 state')
            call simple_mkdir(FROZEN_COPY_DIR)
            call simple_copy_file(sigma_src, frozen_sigma_copy)
            call verify_frozen_sigma
            call spproj_frz%projinfo%set(1, 'projfile', frozen_copy%to_char())
            sigma_src = stemname(frozen_copy)
            call spproj_frz%projinfo%set(1, 'cwd',      sigma_src%to_char())
            call spproj_frz%projinfo%set(1, 'sigma2_state', frozen_sigma_copy%to_char())
            call spproj_frz%write(frozen_copy)
            call superset%mask(spproj)
            call spproj%write_segment_inside('ptcl2D', params%projfile)
            call spproj%write_segment_inside('ptcl3D', params%projfile)
            write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> ABINITIO3D_ADDON MASKED FROZEN/COHORT/NEVER-UPDATED: ', &
                &superset%get_nfrozen(), '/', superset%get_ncohort(), '/', superset%get_nnever_updated()
            write(logfhandle,'(A,A)') '>>> ABINITIO3D_ADDON FROZEN COPY: ', frozen_copy%to_char()
            call sigma_src%kill
        end subroutine addon_prologue

        !> The frozen term is weighted by the base run's committed residual
        !! sigma2 state, never by a re-estimate: its copy must stay byte-equal
        !! to the base run's through every frozen accumulation
        subroutine verify_frozen_sigma
            if( .not. man_addon%matches_artifact('sigma2_state', 0, frozen_sigma_copy) ) &
                &THROW_HARD('the frozen sigma2 state was not consumed as committed by the base run')
        end subroutine verify_frozen_sigma

        !> Cohort poses and labels from random, the frozen sets at every distinct
        !! consuming box of the inherited ladder plus the native box, and the
        !! native frozen-only maps as the stage-3 references (trusted: no CC
        !! pose-initialisation pass).
        subroutine addon_starting_state
            type(frozen_accum) :: store
            type(string) :: src, dest
            character(len=STDLEN) :: msg
            integer, allocatable :: boxes(:)
            integer :: s, b, status
            call reset_ptcl3D_from_ptcl2D_selection
            call spproj%os_ptcl3D%rnd_oris
            if( params%nstates > 1 ) call gen_labelling(spproj%os_ptcl3D, params%nstates, 'uniform')
            call superset%validate_cohort_states(spproj, status, msg)
            if( status /= 0 ) THROW_HARD(trim(msg))
            call spproj%write_segment_inside(params%oritype, params%projfile)
            ! the stage-3 controls define the frozen accumulations
            call set_cline_refine3D(params, start_stage, l_cavgs=.false., addon=addon_ctx)
            call store%new(run_id, trim(params%rec_backend), spproj%os_ptcl3D%get_noris(), &
                &spproj_frz%os_ptcl3D%get_noris(), [(superset%get_nfrozen_state(s), s=1,params%nstates)])
            call store%write(frozen_ctx_fname)
            allocate(boxes(0))
            do istage = start_stage, nstages_refine3D
                b = lpinfo(istage)%box_crop
                if( b == params%box ) cycle
                if( any(boxes == b) ) cycle
                boxes = [boxes, b]
            enddo
            do s = 1, size(boxes)
                call calc_frozen_rec(params, frozen_copy, xrec3D, boxes(s), frozen_ctx_fname)
                call verify_frozen_sigma
            enddo
            call calc_frozen_rec(params, frozen_copy, xrec3D, params%box, frozen_ctx_fname)
            call verify_frozen_sigma
            do s = 1, params%nstates
                src  = refine3D_state_vol_fname(s)
                dest = refine3D_startvol_fname(s)
                call simple_rename(src, dest)
                call rename_support_provenance(src, dest)
                call inject_refine3D_volume(params, s, dest)
                src  = refine3D_state_halfvol_fname(s, 'even')
                dest = refine3D_startvol_half_fname(s, 'even', unfil=.true.)
                call simple_copy_file(src, dest)
                dest = refine3D_startvol_half_fname(s, 'even')
                call simple_rename(src, dest)
                src  = refine3D_state_halfvol_fname(s, 'odd')
                dest = refine3D_startvol_half_fname(s, 'odd', unfil=.true.)
                call simple_copy_file(src, dest)
                dest = refine3D_startvol_half_fname(s, 'odd')
                call simple_rename(src, dest)
                call report_frozen_provenance(s)
            enddo
            call store%kill
            call src%kill
            call dest%kill
        end subroutine addon_starting_state

        !> The frozen-only native map against the base run's registered final
        !! map: agreement shows that poses, halves, sigmas and settings were
        !! reproduced (reported; no refusal in the first release)
        subroutine report_frozen_provenance( s )
            integer, intent(in) :: s
            type(image)  :: vol_frozen, vol_base
            type(string) :: fname
            integer :: box_here, ldim(3), nfoo
            real    :: smpd_here, corr
            if( .not. spproj_frz%isthere_in_osout('vol', s) ) return
            call spproj_frz%get_vol('vol', s, fname, smpd_here, box_here)
            if( .not. file_exists(fname) ) return
            call find_ldim_nptcls(fname, ldim, nfoo)
            if( ldim(1) /= params%box ) return
            call vol_frozen%new(ldim, params%smpd)
            call vol_base%new(ldim, params%smpd)
            call vol_frozen%read(refine3D_startvol_fname(s))
            call vol_base%read(fname)
            corr = vol_frozen%real_corr(vol_base)
            write(logfhandle,'(A,I0,A,F7.4)') '>>> ABINITIO3D_ADDON PROVENANCE: FROZEN-ONLY VS BASE FINAL MAP, STATE ', &
                &s, ', CORRELATION ', corr
            if( corr < 0.9 ) THROW_WARN('the frozen-only map does not reproduce the base run''s final map closely')
            call vol_frozen%kill
            call vol_base%kill
            call fname%kill
        end subroutine report_frozen_provenance

        !> Per stage: the emitted limits next to the base run's (the add-on
        !! plans and promotes by the legacy rule, so they may differ); before
        !! the first trailing stage the boundary reconstruction seeds the
        !! cohort-only chain (trail_seed) with the frozen term added to its maps
        subroutine addon_stage_boundary( istage_here )
            integer, intent(in) :: istage_here
            type(abinitio3D_stage_record) :: base
            if( istage_here <= man_addon%get_nstages() )then
                base = man_addon%get_stage(istage_here)
                if( base%lp_emitted < 0. )then
                    write(logfhandle,'(A,I0,2(A,F7.3),A)') '>>> ABINITIO3D_ADDON STAGE ', istage_here, ' LP ', &
                        &emitted_lp(istage_here), ' LPSTOP ', emitted_lpstop(istage_here), ' (NOT RUN BY THE BASE RUN)'
                else
                    write(logfhandle,'(A,I0,4(A,F7.3),A)') '>>> ABINITIO3D_ADDON STAGE ', istage_here, ' LP ', &
                        &emitted_lp(istage_here), ' LPSTOP ', emitted_lpstop(istage_here), ' (BASE RUN LP ', &
                        &base%lp_emitted, ' LPSTOP ', base%lpstop_emitted, ')'
                endif
            endif
            if( l_cohort_chain_seeded ) return
            if( .not. cline_refine3D%defined('trail_rec') ) return
            if( cline_refine3D%get_carg('trail_rec') /= 'yes' ) return
            write(logfhandle,'(A,I0)') '>>> ABINITIO3D_ADDON: SEEDING THE COHORT TRAILING CHAIN BEFORE STAGE ', istage_here
            call calc_rec(params, params%projfile, xrec3D, istage_here)
            l_cohort_chain_seeded = .true.
        end subroutine addon_stage_boundary

        !> The union before its final reconstruction: the frozen rows' 3D
        !! records and every saved ptcl2D state back from the frozen
        !! project. The cohort-only sigma2 registration is
        !! dropped and the frozen term leaves the reconstruction command line:
        !! the final reconstruction reads every particle and finds no
        !! consumable state, so it bootstraps the union's (image-power seed,
        !! bootstrap map, one residual pass) exactly as bootstrap_rec3D does
        !! for any project without consumable sigmas
        subroutine addon_restore_union
            call spproj%read(params%projfile)
            call superset%restore(spproj, spproj_frz)
            if( spproj%projinfo%isthere(1, 'sigma2_state') ) call spproj%projinfo%delete_entry('sigma2_state')
            call spproj%write(params%projfile)
            call cline_refine3D%delete('frozen_rec')
            write(logfhandle,'(A,I0,A)') '>>> ABINITIO3D_ADDON: RESTORED ', superset%get_nfrozen(), &
                &' FROZEN PARTICLES; THE FINAL RECONSTRUCTION BOOTSTRAPS THE UNION SIGMA2 STATE'
        end subroutine addon_restore_union

        !> The cohort-only diagnostic (addon_diag=yes), union-aware res/res05
        !! for every row and the validation against the base solution: FSC
        !! verdict, map correlation at the base resolution, the cohort-only map
        !! against the base map, and both runs' stage limits
        !! (abinitio3D_addon_report.txt; a regression is warned about, the
        !! result is published regardless)
        subroutine addon_epilogue
            type(abinitio3D_addon_report) :: report
            type(abinitio3D_stage_record) :: base_stage
            real, allocatable :: fsc(:), res(:), fsc_base(:)
            type(string) :: fsc_name, vol_base
            real    :: fsc05, fsc0143, lp_base, lpstop_base, smpd_vol
            integer :: s, i, box_fsc, box_vol
            logical :: l_union, l_base
            if( trim(params%addon_diag) == 'yes' ) call addon_cohort_diagnostic
            call spproj%read(params%projfile)
            call report%new(params%nstates, size(emitted_lp), params%smpd, params%mskdiam)
            do i = start_stage, size(emitted_lp)
                if( emitted_lp(i) < 0. ) cycle
                lp_base     = -1.
                lpstop_base = -1.
                if( i <= man_addon%get_nstages() )then
                    base_stage  = man_addon%get_stage(i)
                    lp_base     = base_stage%lp_emitted
                    lpstop_base = base_stage%lpstop_emitted
                endif
                call report%set_stage(i, emitted_lp(i), emitted_lpstop(i), lp_base, lpstop_base)
            enddo
            res = get_resarr(params%box, params%smpd)
            do s = 1, params%nstates
                l_union  = .false.
                fsc_name = refine3D_fsc_fname(s)
                if( file_exists(fsc_name) )then
                    fsc     = file2rarr(fsc_name)
                    l_union = size(fsc) == size(res)
                endif
                if( l_union )then
                    call get_resolution(fsc, res, fsc05, fsc0143)
                    do i = 1, spproj%os_ptcl3D%get_noris()
                        if( spproj%os_ptcl3D%get_state(i) /= s ) cycle
                        call spproj%os_ptcl3D%set(i, 'res',   fsc0143)
                        call spproj%os_ptcl3D%set(i, 'res05', fsc05)
                    enddo
                endif
                call report%set_populations(s, spproj%os_ptcl3D%get_pop(s, 'state'), superset%get_nfrozen_state(s))
                l_base = .false.
                if( spproj_frz%isthere_in_osout('fsc', s) )then
                    call spproj_frz%get_fsc(s, fsc_name, box_fsc)
                    if( file_exists(fsc_name) )then
                        fsc_base = file2rarr(fsc_name)
                        l_base   = size(fsc_base) == size(res)
                    endif
                endif
                if( l_union .and. l_base ) call report%compare_fsc(s, fsc_base, fsc, res)
                if( .not. spproj_frz%isthere_in_osout('vol', s) ) cycle
                call spproj_frz%get_vol('vol', s, vol_base, smpd_vol, box_vol)
                call report%compare_maps(s, vol_base, refine3D_state_vol_fname(s))
                if( trim(params%addon_diag) == 'yes' ) call report%compare_cohort(s, vol_base, &
                    &string(ADDON_DIAG_DIR//'/')//refine3D_state_vol_fname(s))
            enddo
            call spproj%write(params%projfile)
            call report%print
            call report%write(string(ADDON_REPORT_FNAME))
            write(logfhandle,'(A,A)') '>>> ABINITIO3D_ADDON REPORT WRITTEN: ', ADDON_REPORT_FNAME
            if( report%any_regressed() ) THROW_WARN('abinitio3D_addon: a state regressed against the base solution')
            call report%kill
            call fsc_name%kill
            call vol_base%kill
        end subroutine addon_epilogue

        !> addon_diag=yes: the cohort alone at the native box, without the
        !! frozen term, on a copy with the frozen rows masked again and the
        !! union's sigma2 state, in its own directory so no output of the run
        !! is touched
        subroutine addon_cohort_diagnostic
            type(sp_project) :: spproj_diag
            type(cmdline)    :: cline_diag
            type(string)     :: diag_proj, sigma_path, cwd_run
            logical :: found
            integer :: status, s
            call simple_mkdir(ADDON_DIAG_DIR)
            call simple_getcwd(cwd_run)
            ! the working copy's file name keeps the sigma2 layout lineage
            diag_proj = simple_abspath(string(ADDON_DIAG_DIR//'/')//basename(params%projfile), check_exists=.false.)
            call spproj_diag%read(params%projfile)
            call superset%mask(spproj_diag)
            call spproj_diag%get_sigma2_state_path(sigma_path, found)
            if( found )then
                sigma_path = simple_abspath(sigma_path)
                call spproj_diag%projinfo%set(1, 'sigma2_state', sigma_path%to_char())
            endif
            call spproj_diag%projinfo%set(1, 'projfile', diag_proj%to_char())
            call spproj_diag%write(diag_proj)
            call spproj_diag%kill
            cline_diag = cline_reconstruct3D
            call apply_refine3D_reconstruction_controls(cline_diag)
            call cline_diag%delete('frozen_rec')
            call cline_diag%set('prg',       'reconstruct3D')
            call cline_diag%set('mkdir',     'no')
            call cline_diag%set('projfile',  diag_proj)
            call cline_diag%set('pgrp',      params%pgrp)
            call cline_diag%set('trail_rec', 'no')
            call cline_diag%delete('box_crop')
            call cline_diag%delete('update_frac')
            call cline_diag%delete('trail_seed')
            do s = 1, params%nstates
                call cline_diag%delete('vol'//int2str(s))
            enddo
            call strip_refine3D_planning_keys(cline_diag)
            call simple_chdir(string(ADDON_DIAG_DIR), status)
            if( status /= 0 ) THROW_HARD('cannot enter the add-on diagnostic directory')
            write(logfhandle,'(A)') '>>> ABINITIO3D_ADDON: COHORT-ONLY DIAGNOSTIC RECONSTRUCTION (addon_diag/)'
            call xrec3D%execute(cline_diag)
            call simple_chdir(cwd_run, status)
            if( status /= 0 ) THROW_HARD('cannot leave the add-on diagnostic directory')
            call cline_diag%kill
            call diag_proj%kill
            call sigma_path%kill
            call cwd_run%kill
        end subroutine addon_cohort_diagnostic

        !> Describe the completed run in its manifest and register it in the
        !! project. Every failure is reported and leaves the completed run as
        !! it is: the manifest is published last, atomically, or not at all.
        subroutine write_run_manifest( l_eligible, program_name )
            logical,          intent(in) :: l_eligible
            character(len=*), intent(in) :: program_name
            type(abinitio3D_manifest)                  :: man
            type(abinitio3D_stage_record), allocatable :: stages(:)
            type(sp_project)      :: spproj_man
            type(string)          :: fname
            character(len=STDLEN) :: msg
            integer :: i, status
            call spproj_man%read(params%projfile)
            call man%new(run_id, program_name, l_eligible, spproj_man, trim(params%ptcl_src))
            call man%set_solution(nstates_glob, params%pgrp, params%box, params%smpd, params%mskdiam, &
                &params%multivol_mode, split_stage)
            call man%set_sampling(params%nsample, nptcls_eff, update_frac, l_force_full_sampling)
            call man%set_stage_line(l_refine3D_lp_override, params%lp, l_refine3D_lpstop_override, params%lpstop)
            allocate(stages(size(lpinfo)))
            do i = 1, size(lpinfo)
                stages(i)%lp_planned     = lpinfo(i)%lp
                stages(i)%lp_emitted     = emitted_lp(i)
                stages(i)%lpstop_emitted = emitted_lpstop(i)
                stages(i)%box_crop       = lpinfo(i)%box_crop
                stages(i)%smpd_crop      = lpinfo(i)%smpd_crop
                stages(i)%scale          = lpinfo(i)%scale
                stages(i)%trslim         = lpinfo(i)%trslim
                stages(i)%frc_crit       = lpinfo(i)%frc_crit
                stages(i)%l_autoscale    = lpinfo(i)%l_autoscale
                stages(i)%l_lpset        = lpinfo(i)%l_lpset
            enddo
            call man%set_ladder(start_stage, nstages_refine3D, stages)
            call man%record_inputs(cline_entry)
            call man%record_artifacts(spproj_man)
            fname = MANIFEST_FNAME
            call man%write(fname, status, msg)
            if( status == 0 )then
                call man%register(spproj_man, MANIFEST_FNAME)
                call spproj_man%write_segment_inside('projinfo', params%projfile)
                write(logfhandle,'(A,A)') '>>> ABINITIO3D RUN MANIFEST WRITTEN: ', trim(run_id)
            else
                THROW_WARN('abinitio3D run manifest not written: '//trim(msg))
            endif
            call man%kill
            call spproj_man%kill
            call fname%kill
        end subroutine write_run_manifest

        subroutine handoff_split_checkpoint_to_refine3D_states
            type(cmdline) :: cline_states
            integer       :: state, first_iter, remaining_niters, nsample_handoff
            cline_states = cline_refine3D
            first_iter = next_refine3D_iteration()
            remaining_niters = abinitio_remaining_niters(split_stage, nstages_refine3D)
            if( remaining_niters < 1 )then
                THROW_HARD('abinitio3D split checkpoint has no remaining refine3D_states iterations')
            endif
            nsample_handoff = min(nptcls_eff, max(1, nint(update_frac * real(nptcls_eff))))
            call cline_states%set('prg',          'refine3D_states')
            call cline_states%set('mkdir',                       'no')
            call cline_states%set('pose_policy',              'local')
            call cline_states%set('flex',                        'no') ! the checkpoint already carries the states
            call cline_states%set('nstates',            nstates_glob)
            call cline_states%set('nsample',          nsample_handoff)
            call cline_states%set('maxits',          remaining_niters)
            call cline_states%set('lpstart',   lpinfo(split_stage)%lp)
            call cline_states%set('lpstop', lpinfo(nstages_refine3D)%lp)
            call cline_states%set('startit',              first_iter)
            call cline_states%set('which_iter',           first_iter)
            call cline_states%set('extr_iter',            first_iter)
            call cline_states%set('filt_mode',    'nonuniform_lpset')
            if( l_force_full_sampling )then
                call cline_states%set('sticky_class_sampling', 'no')
            else
                call cline_states%set('sticky_class_sampling', 'yes')
            endif
            call cline_states%delete('multivol_mode')
            call cline_states%delete('prob_neigh_mode')
            call cline_states%delete('refine')
            call cline_states%delete('nspace')
            call cline_states%delete('nspace_sub')
            call cline_states%delete('lp')
            call cline_states%delete('minits')
            call cline_states%delete('endit')
            ! refine3D_states owns the state-overlap convergence policy; the
            ! inherited abinitio stage target must not override its default
            call cline_states%delete('overlap')
            ! refine3D_states rejects vol1..volN inputs and takes its starting
            ! state maps from the project out segment, where calc_rec registered
            ! the split-checkpoint reconstructions
            call spproj%read_segment('out', params%projfile)
            do state = 1,nstates_glob
                if( .not. spproj%isthere_in_osout('vol', state) )then
                    THROW_HARD('abinitio3D split checkpoint did not register every state volume in the project')
                endif
                call cline_states%delete('vol'//int2str(state))
            enddo
            write(logfhandle,'(A,I0,A,I0,A,I0,A,F7.2,A,F7.2)') &
                &'>>> ABINITIO3D -> REFINE3D_STATES FIRST_ITER/MAXITS/NSAMPLE/LPSTART/LPSTOP: ', &
                &first_iter, '/', remaining_niters, '/', nsample_handoff, '/', &
                &lpinfo(split_stage)%lp, '/', lpinfo(nstages_refine3D)%lp
            call xrefine3D_states%execute(cline_states)
            call spproj%read(params%projfile)
            call cline_states%kill
        end subroutine handoff_split_checkpoint_to_refine3D_states

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
            ! the stage-LP promotion reads res05: no resolution left by an
            ! earlier refinement of this project may promote the first stage
            call spproj%os_ptcl3D%delete_entry('res')
            call spproj%os_ptcl3D%delete_entry('res05')
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

    !> a run identifier unique to this process and moment, blank-free
    function new_run_id( prefix ) result( id )
        character(len=*), intent(in) :: prefix
        character(len=64) :: id
        integer :: v(8)
        call date_and_time(values=v)
        write(id,'(A,A,I4.4,2I2.2,A,3I2.2,A,I3.3,A,I0)') trim(prefix), '_', v(1), v(2), v(3), 'T', v(5), v(6), v(7), &
            &'.', v(8), '_', get_process_id()
    end function new_run_id

end module simple_commanders_abinitio
