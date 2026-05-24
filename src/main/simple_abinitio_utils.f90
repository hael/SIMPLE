!@descr: utilities for ab initio 3D reconstruction used by commanders_abinitio
module simple_abinitio_utils
use simple_commanders_api
use simple_commanders_rec, only: commander_bootstrap_rec3D
use simple_commanders_volops, only: commander_symmetrize_map
use simple_cluster_seed,      only: gen_labelling
use simple_class_frcs,        only: class_frcs
use simple_euclid_sigma2,     only: sigma2_star_from_iter
use simple_matcher_refvol_utils, only: remove_ref_section_files
use simple_parameters,        only: parameters
use simple_refine3D_fnames,   only: refine3D_fsc_fname, refine3D_startvol_fbody, &
    &refine3D_startvol_fname, refine3D_startvol_half_fname, &
    &refine3D_state_halfvol_fname, refine3D_state_vol_fbody, refine3D_state_vol_fname
implicit none
#include "simple_local_flags.inc"

! singleton constants
character(len=*), parameter :: REC_FBODY             = 'rec_final_state'
character(len=*), parameter :: STR_STATE_GLOB        = '01'
real,             parameter :: LPSTOP_BOUNDS(2)      = [4.5,6.0]
real,             parameter :: LPSTART_BOUNDS(2)     = [10.,20.]
real,             parameter :: CENLP_DEFAULT         = 30.
real,             parameter :: LPSYMSRCH_LB          = 12.
real,             parameter :: UPDATE_FRAC_MAX       = 0.9                  ! to ensure fractional update is always on
real,             parameter :: LPSTART_INI3D         = 20.                  ! Default lpstart for abinitio3D_cavgs/cavgs_ini
real,             parameter :: LPSTOP_INI3D          = 8.                   ! Default lpstop for abinitio3D_cavgs/cavgs_ini
integer,          parameter :: NSTAGES               = 8
integer,          parameter :: NSTAGES_INI3D         = 4 ! # of ini3D stages used for initialization
integer,          parameter :: NSTAGES_INI3D_MAX     = 7
integer,          parameter :: MAXITS(8)             = [20,20,17,17,17,15,15,30]
integer,          parameter :: MAXITS_GLOB           = SUM(MAXITS(1:7))     ! the last stage is omitted in this estimate since the sampling method changes
integer,          parameter :: SYMSRCH_STAGE         = 3
integer,          parameter :: TRAILREC_STAGE_SINGLE = 5                    ! first stage where trail_rec behavior changes
integer,          parameter :: TRAILREC_STAGE_MULTI  = NSTAGES              ! first stage where trail_rec is enabled for independent mode
integer,          parameter :: AUTOMSK_STAGE         = 6                    ! switch on automasking when lpauto is switched on
integer,          parameter :: HET_DOCKED_STAGE      = NSTAGES              ! stage at which state splitting is done when multivol_mode==docked
integer,          parameter :: STREAM_ANALYSIS_STAGE = 5                    ! when streaming on some analysis will be performed
integer,          parameter :: GAUREF_LAST_STAGE     = 2                    ! When to stop gaussian filtering in early stages
integer,          parameter :: MAXITS_BETWEEN        = 10                   ! Development
integer,          parameter :: NSAMPLE_ABINITIO3D_DEFAULT = 10000           ! default # particles sampled per state

! singleton variables
type(lp_crop_inf), allocatable :: lpinfo(:)
logical          :: l_srch4symaxis    = .false., l_symran        = .false.
logical          :: l_ini3D           = .false.
logical          :: l_lpauto          = .false.
logical          :: l_automsk         = .false.
logical          :: l_nonuniform      = .false.
logical          :: l_staged_nonuniform_mode = .false.
type(sym)        :: se1, se2
type(cmdline)    :: cline_refine3D, cline_symmap, cline_reconstruct3D, cline_reproject
real             :: update_frac  = 1.0
integer          :: nstates_glob = 1, nptcls_eff = 0
integer          :: nstages_refine3D = NSTAGES

! In submodule: simple_abinitio_controller.f90
interface
    module subroutine set_cline_refine3D( params, istage, l_cavgs )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: istage
        logical,           intent(in) :: l_cavgs
    end subroutine set_cline_refine3D
end interface

contains

    subroutine prep_class_command_lines( params, cline, projfile )
        class(parameters), intent(in) :: params
        class(cmdline),    intent(in) :: cline
        class(string),     intent(in) :: projfile
        cline_refine3D      = cline
        cline_symmap        = cline
        cline_reconstruct3D = cline
        cline_reproject     = cline
        ! refine3D
        call cline_refine3D%set('prg',                    'refine3D')
        call cline_refine3D%set('pgrp',                  params%pgrp)
        call cline_refine3D%set('projfile',                 projfile)
        call cline_refine3D%delete('box')
        call cline_refine3D%delete('smpd')
        call cline_refine3D%delete('smpd_crop')
        ! symmetrization
        call cline_symmap%set('prg',                'symmetrize_map')
        call cline_symmap%set('pgrp',                    params%pgrp)
        call cline_symmap%set('projfile',                   projfile)
        call cline_symmap%set('center',                        'yes')
        if( .not. cline_symmap%defined('cenlp') )then
        call cline_symmap%set('cenlp',                 CENLP_DEFAULT)
        endif
        call cline_symmap%set('hp',                        params%hp)
        ! re-reconstruct volume
        call cline_reconstruct3D%set('prg',          'reconstruct3D')
        call cline_reconstruct3D%set('projfile',            projfile)
        call cline_reconstruct3D%set('pgrp',             params%pgrp)
        call cline_reconstruct3D%set('ml_reg',                  'no')
        call cline_reconstruct3D%set('objfun',                  'cc')
        call cline_reconstruct3D%delete('box')
        call cline_reconstruct3D%delete('smpd')
        call cline_reconstruct3D%delete('smpd_crop')
        ! no fractional update
        call cline_reconstruct3D%delete('update_frac')
        call cline_reconstruct3D%delete('refs')
        call cline_reconstruct3D%delete('refs_even')
        call cline_reconstruct3D%delete('refs_odd')
        ! re-project volume, only with cavgs
        call cline_reproject%set('prg',                  'reproject')
        call cline_reproject%set('pgrp',                 params%pgrp)
        call cline_reproject%set('outstk',        'reprojs'//MRC_EXT)
        call cline_reproject%set('smpd',                 params%smpd)
        call cline_reproject%set('box',                   params%box)
        call cline_reproject%set('oritab',          'final_oris.txt')
        call cline_reproject%set('nstates',           params%nstates)
        call cline_reproject%delete('projfile')
    end subroutine prep_class_command_lines

    subroutine strip_refine3D_planning_keys( child_cline, delete_which_iter )
        class(cmdline), intent(inout) :: child_cline
        logical, optional, intent(in) :: delete_which_iter
        logical :: l_delete_which_iter
        l_delete_which_iter = .false.
        if( present(delete_which_iter) ) l_delete_which_iter = delete_which_iter
        if( l_delete_which_iter ) call child_cline%delete('which_iter')
        call child_cline%delete('endit')
        call child_cline%delete('automsk')
        call child_cline%delete('filt_mode')
        call child_cline%delete('refs')
        call child_cline%delete('refs_even')
        call child_cline%delete('refs_odd')
        call child_cline%delete('box')
        call child_cline%delete('smpd')
        call child_cline%delete('smpd_crop')
    end subroutine strip_refine3D_planning_keys

    subroutine inject_refine3D_volume( params, state, vol )
        class(parameters), intent(inout) :: params
        integer,           intent(in)    :: state
        class(string),     intent(in)    :: vol
        type(string) :: vol_key
        vol_key = 'vol'//int2str(state)
        call cline_refine3D%set(vol_key%to_char(), vol)
        params%vols(state) = vol
        call remove_ref_section_files
        call vol_key%kill
    end subroutine inject_refine3D_volume

    subroutine register_stage_volume( params, state, vol_name, projfile )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: state
        class(string),     intent(in) :: vol_name
        class(string),     intent(in), optional :: projfile
        type(sp_project)            :: spproj
        type(string) :: fsc_name
        integer :: ldim(3), nptcls
        real    :: smpd
        if( .not. file_exists(vol_name) ) return
        call find_ldim_nptcls(vol_name, ldim, nptcls)
        smpd = params%smpd_crop
        if( present(projfile) )then
            call spproj%read_segment('out', projfile)
        else
            call spproj%read_segment('out', params%projfile)
        endif
        call spproj%add_vol2os_out(vol_name, smpd, state, 'vol')
        fsc_name = refine3D_fsc_fname(state)
        if( file_exists(fsc_name) ) call spproj%add_fsc2os_out(fsc_name, state, ldim(1))
        if( present(projfile) )then
            call spproj%write_segment_inside('out', projfile)
        else
            call spproj%write_segment_inside('out', params%projfile)
        endif
        call spproj%kill
        call fsc_name%kill
    end subroutine register_stage_volume

    subroutine write_abinitio_lowpass_snapshot( vol_in, lp, vol_out, smpd )
        class(string), intent(in) :: vol_in, vol_out
        real,          intent(in) :: lp, smpd
        type(image) :: vol_lp
        integer :: ldim(3), nptcls
        real    :: lp_eff
        if( .not. file_exists(vol_in) ) return
        call find_ldim_nptcls(vol_in, ldim, nptcls)
        lp_eff = max(2.0 * smpd, lp)
        call vol_lp%new(ldim, smpd)
        call vol_lp%read(vol_in)
        call vol_lp%fft()
        call vol_lp%bp(0., lp_eff)
        call vol_lp%ifft()
        call vol_lp%write(vol_out, del_if_exists=.true.)
        call vol_lp%kill
    end subroutine write_abinitio_lowpass_snapshot

    real function abinitio_state_fsc_lowpass( state, box, smpd, fallback_lp ) result( lp )
        integer, intent(in) :: state, box
        real,    intent(in) :: smpd, fallback_lp
        type(string) :: fsc_name
        real, allocatable :: fsc(:), res(:)
        real :: fsc05, fsc0143
        lp = fallback_lp
        fsc_name = refine3D_fsc_fname(state)
        if( file_exists(fsc_name) )then
            fsc = file2rarr(fsc_name)
            res = get_resarr(box, smpd)
            call get_resolution(fsc, res, fsc05, fsc0143)
            if( fsc0143 > 0. .and. fsc0143 == fsc0143 ) lp = fsc0143
        endif
        if( allocated(fsc) ) deallocate(fsc)
        if( allocated(res) ) deallocate(res)
        call fsc_name%kill
    end function abinitio_state_fsc_lowpass

    subroutine set_symmetry_class_vars( params )
        class(parameters), intent(in) :: params
        type(string) :: pgrp, pgrp_start
        pgrp           = lowercase(trim(params%pgrp))
        pgrp_start     = lowercase(trim(params%pgrp_start))
        l_srch4symaxis = pgrp .ne. pgrp_start
        l_symran       = .false.
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

    subroutine set_lplims_from_frcs( params, spproj, l_cavgs, lpstart, lpstop )
        class(parameters), intent(inout) :: params
        class(sp_project), intent(in)    :: spproj
        logical,           intent(in)    :: l_cavgs
        real, optional,    intent(in)    :: lpstart, lpstop
        real,         allocatable :: frcs_avg(:)
        integer,      allocatable :: states(:)
        type(string)     :: frcs_fname
        type(class_frcs) :: clsfrcs
        real             :: lpfinal
        integer          :: filtsz
        if( trim(params%force_lp_range).eq.'yes' )then
            if( .not.(present(lpstart) .and. present(lpstop)) )then
                THROW_HARD('force_lp_range=yes requires both lpstart and lpstop')
            endif
            if( allocated(lpinfo) ) deallocate(lpinfo)
            allocate(lpinfo(NSTAGES))
            call lpstages_fast(params%box, NSTAGES, params%smpd, lpstart, lpstop, lpinfo)
            return
        endif
        ! retrieve FRC info
        call spproj%get_frcs(frcs_fname, 'frc2D', fail=.false.)
        ! work out low-pass limits and downscaling parameters
        params%frcs = frcs_fname
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
            call lpstages(params%box, NSTAGES, frcs_avg, params%smpd,&
            &lpstart, lpstart, lpfinal, lpinfo, l_cavgs, verbose=.true.)
        else
            call lpstages(params%box, NSTAGES, frcs_avg, params%smpd,&
            &LPSTART_BOUNDS(1), LPSTART_BOUNDS(2), lpfinal, lpinfo, l_cavgs, verbose=.true.)
        endif
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

    subroutine exec_refine3D( params, istage, xrefine3D )
        class(parameters),     intent(inout) :: params
        integer,               intent(in)    :: istage
        class(commander_base), intent(inout) :: xrefine3D
        type(string) :: stage, vol_name, vol_stage, vol_lp_stage
        integer      :: state
        real         :: lp_snapshot
        call cline_refine3D%delete('endit')
        call xrefine3D%execute(cline_refine3D)
        call del_files(DIST_FBODY,      params%nparts,ext='.dat')
        call del_files(ASSIGNMENT_FBODY,params%nparts,ext='.dat')
        call del_file(DIST_FBODY//'.dat')
        call del_file(ASSIGNMENT_FBODY//'.dat')
        stage = '_stage'//int2str_pad(istage,2)
        do state = 1, params%nstates
            vol_name  = refine3D_state_vol_fname(state)
            vol_stage = add2fbody(vol_name, string(MRC_EXT),stage)
            vol_lp_stage = add2fbody(vol_stage, MRC_EXT, LP_SUFFIX)
            if( file_exists(vol_name) )then
                lp_snapshot = abinitio_state_fsc_lowpass(state, lpinfo(istage)%box_crop, &
                    &lpinfo(istage)%smpd_crop, lpinfo(istage)%lp)
                call write_abinitio_lowpass_snapshot(vol_name, lp_snapshot, vol_lp_stage, lpinfo(istage)%smpd_crop)
            endif
        enddo
        call vol_stage%kill
        call vol_lp_stage%kill
    end subroutine exec_refine3D

    subroutine symmetrize( params, istage, spproj, projfile, xrec3D )
        class(parameters),               intent(inout) :: params
        integer,                         intent(in)    :: istage
        class(sp_project),               intent(inout) :: spproj
        class(string),                   intent(in)    :: projfile
        class(commander_base), optional, intent(inout) :: xrec3D
        type(commander_symmetrize_map) :: xsymmap
        type(cmdline)                  :: cline_symrec
        type(string) :: vol_iter, vol_sym, stage, vol_stage, vol_lp_stage
        integer :: state
        real :: lp_snapshot
        real :: lpsym
        if( l_symran )then
            call se1%symrandomize(spproj%os_ptcl3D)
            call spproj%write_segment_inside('ptcl3D', projfile)
        endif
        if( l_srch4symaxis )then
            lpsym = max(LPSYMSRCH_LB,lpinfo(SYMSRCH_STAGE)%lp)
            write(logfhandle,'(A,F5.1)') '>>> DID SET MAP SYMMETRIZATION LOW-PASS LIMIT (IN A) TO: ', lpsym
            write(logfhandle,'(A)') '>>>'
            if( params%nstates > 1 )then
                write(logfhandle,'(A)') '>>> STATE-WISE MAP SYMMETRIZATION'
            else
                write(logfhandle,'(A)') '>>> MAP SYMMETRIZATION'
            endif
            write(logfhandle,'(A)') '>>>'
            call cline_symmap%set('smpd', lpinfo(istage)%smpd_crop)
            call cline_symmap%set('box',  lpinfo(istage)%box_crop)
            call cline_symmap%set('lp', lpsym)
            do state = 1,params%nstates
                vol_iter = refine3D_state_vol_fname(state)
                if( .not. file_exists(vol_iter) )then
                    THROW_HARD('input volume to map symmetrization does not exist for state '//int2str(state))
                endif
                call cline_symmap%set('vol1', vol_iter)
                if( params%nstates > 1 )then
                    vol_sym = 'symmetrized_map_state'//int2str_pad(state,2)//MRC_EXT
                    call cline_symmap%set('state', state)
                    write(logfhandle,'(A,I0)') '>>> MAP SYMMETRIZATION STATE ', state
                else
                    vol_sym = 'symmetrized_map'//MRC_EXT
                    call cline_symmap%delete('state')
                endif
                call cline_symmap%set('outvol', vol_sym)
                call xsymmap%execute(cline_symmap)
                call del_file('SYMAXIS_SEARCH_FINISHED')
            enddo
            if( present(xrec3D) )then
                ! symmetric reconstruction
                cline_symrec = cline_refine3D
                call cline_symrec%set('prg',        'reconstruct3D')
                call cline_symrec%set('mkdir',      'no')
                call cline_symrec%set('projfile',   projfile)
                call cline_symrec%set('pgrp',       params%pgrp)
                call cline_symrec%set('which_iter', cline_refine3D%get_iarg('endit'))
                call strip_refine3D_planning_keys(cline_symrec)
                call xrec3D%execute(cline_symrec)
                do state = 1,params%nstates
                    vol_sym = refine3D_state_vol_fname(state)
                    if( params%nstates > 1 )then
                        call simple_copy_file(vol_sym, string('symmetric_map_state')//int2str_pad(state,2)//MRC_EXT)
                    else
                        call simple_copy_file(vol_sym, string('symmetric_map')//MRC_EXT)
                    endif
                enddo
                call cline_symrec%kill
            endif
            stage        = '_stage'//int2str_pad(istage,2)
            do state = 1,params%nstates
                vol_sym      = refine3D_state_vol_fname(state)
                vol_stage    = add2fbody(vol_sym, string(MRC_EXT), stage)
                vol_lp_stage = add2fbody(vol_stage, MRC_EXT, LP_SUFFIX)
                lp_snapshot  = abinitio_state_fsc_lowpass(state, lpinfo(istage)%box_crop, &
                    &lpinfo(istage)%smpd_crop, lpinfo(istage)%lp)
                call write_abinitio_lowpass_snapshot(vol_sym, lp_snapshot, vol_lp_stage, lpinfo(istage)%smpd_crop)
                call inject_refine3D_volume(params, state, vol_sym)
                call vol_stage%kill
                call vol_lp_stage%kill
            enddo
        endif
    end subroutine symmetrize

    ! Performs reconstruction at selected stage boundaries.
    subroutine calc_rec( params, projfile, xrec3D, istage )
        class(parameters),       intent(inout) :: params
        class(string),           intent(in)    :: projfile
        class(commander_base),   intent(inout) :: xrec3D
        integer,                 intent(in)    :: istage
        type(string)      :: vol_even, vol_odd, tmpl, src, dest, dest_main, dest_even, dest_odd, sstate, sstage, pgrp, vol_diag
        type(cmdline)     :: cline_rec
        integer           :: state
        real              :: lp_snapshot
        logical           :: have_even_stage, have_odd_stage
        ! Reconstruction
        pgrp = trim(params%pgrp)
        if( istage <= SYMSRCH_STAGE ) pgrp = trim(params%pgrp_start)
        cline_rec = cline_refine3D
        call cline_rec%set('prg',       'reconstruct3D')
        call cline_rec%set('mkdir',     'no')
        call cline_rec%set('projfile',  projfile)
        call cline_rec%set('pgrp',      pgrp)
        call cline_rec%set('box_crop',  lpinfo(istage)%box_crop)
        call cline_rec%set('trail_rec', 'no')
        if( cline_rec%get_carg('ml_reg').ne.'yes' ) call cline_rec%set('objfun','cc')
        do state = 1,params%nstates
            call cline_rec%delete('vol'//int2str(state))
        enddo
        call cline_rec%delete('vol_even')
        call cline_rec%delete('vol_odd')
        call cline_rec%delete('update_frac')
        call strip_refine3D_planning_keys(cline_rec)
        call xrec3D%execute(cline_rec)
        ! Rename volumes, update cline & project
        sstage  = int2str_pad(istage-1,2)
        do state = 1,params%nstates
            sstate = int2str_pad(state,2)
            ! Rename volumes
            if( istage == 1 )then
                tmpl = refine3D_startvol_fbody(state)
            else
                tmpl = refine3D_state_vol_fbody(state)//'_stage'//sstage
            endif
            have_even_stage = .false.
            have_odd_stage  = .false.
            src  = refine3D_state_vol_fname(state)
            dest_main = tmpl//MRC_EXT
            call simple_rename(src, dest_main)
            vol_diag = add2fbody(dest_main, MRC_EXT, LP_SUFFIX)
            lp_snapshot = abinitio_state_fsc_lowpass(state, lpinfo(istage)%box_crop, &
                &lpinfo(istage)%smpd_crop, lpinfo(istage)%lp)
            call write_abinitio_lowpass_snapshot(dest_main, lp_snapshot, vol_diag, lpinfo(istage)%smpd_crop)
            vol_even = refine3D_state_halfvol_fname(state, 'even')
            if( file_exists(vol_even) )then
                dest = tmpl//'_even_unfil'//MRC_EXT
                call simple_copy_file(vol_even, dest)
                dest = tmpl//'_even'//MRC_EXT
                call simple_rename(vol_even, dest)
                dest_even = dest
                have_even_stage = .true.
            endif
            vol_odd  = refine3D_state_halfvol_fname(state, 'odd')
            if( file_exists(vol_odd) )then
                dest = tmpl//'_odd_unfil'//MRC_EXT
                call simple_copy_file(vol_odd, dest)
                dest = tmpl//'_odd'//MRC_EXT
                call simple_rename(vol_odd, dest)
                dest_odd = dest
                have_odd_stage = .true.
            endif
            ! Update refine3D command line
            call inject_refine3D_volume(params, state, dest_main)
            if( have_even_stage .and. have_odd_stage )then
                params%vols_even(state) = dest_even
                params%vols_odd(state)  = dest_odd
            endif
            ! Update project
            call register_stage_volume(params, state, dest_main, projfile)
        enddo
        call vol_diag%kill
        call cline_rec%kill
    end subroutine calc_rec

    subroutine randomize_states( params, spproj, projfile, xrec3D, istage )
        use simple_commanders_euclid,  only: commander_calc_group_sigmas
        class(parameters),     intent(inout) :: params
        class(sp_project),     intent(inout) :: spproj
        class(string),         intent(in)    :: projfile
        class(commander_base), intent(inout) :: xrec3D
        integer,               intent(in)    :: istage
        type(commander_calc_group_sigmas) :: xcalc_group_sigmas
        type(cmdline)                     :: cline_calc_group_sigmas
        call spproj%read_segment('ptcl3D', projfile)
        call gen_labelling(spproj%os_ptcl3D, params%nstates, 'squared_uniform')
        call spproj%write_segment_inside(params%oritype, projfile)
        call cline_refine3D%set(     'nstates', params%nstates)
        call cline_reconstruct3D%set('nstates', params%nstates)
        call cline_reproject%set(    'nstates', params%nstates)
        ! Sigma2 need be to consolidated because we bypass the general
        ! path to reconstruct new volumes
        if( cline_refine3D%get_carg('ml_reg').eq.'yes' )then
            cline_calc_group_sigmas = cline_refine3D
            call cline_calc_group_sigmas%set('prg', 'calc_group_sigmas')
            call xcalc_group_sigmas%execute(cline_calc_group_sigmas)
            call cline_calc_group_sigmas%kill
        endif
        ! Multi-state reconstruction
        call calc_rec(params, projfile, xrec3D, istage)
    end subroutine randomize_states

    subroutine gen_ortho_reprojs4viz( params, spproj )
        class(parameters), intent(in)    :: params
        type(sp_project),  intent(inout) :: spproj
        type(string) :: str_state, fname
        type(image)  :: final_vol, reprojs
        integer      :: state, ifoo, ldim(3)
        real         :: smpd
        call spproj%read_segment('out', params%projfile)
        do state = 1, params%nstates
            if( .not.spproj%isthere_in_osout('vol', state) )cycle   ! empty-state case
            str_state = int2str_pad(state,2)
            fname = refine3D_state_vol_fname(state)
            if( .not. file_exists(fname) )cycle
            exit
        enddo
        call find_ldim_nptcls(fname, ldim, ifoo)
        smpd = params%smpd
        call final_vol%new(ldim, smpd)
        do state = 1, params%nstates
            str_state = int2str_pad(state,2)
            if( spproj%isthere_in_osout('vol', state) )then
                fname = refine3D_state_vol_fname(state)
                if( .not. file_exists(fname) )cycle
                call final_vol%read(fname)
                call final_vol%generate_orthogonal_reprojs(reprojs)
                call reprojs%write_jpg(string('orthogonal_reprojs_state')//str_state//'.jpg')
                call reprojs%kill
            endif
        enddo
        call final_vol%kill
    end subroutine gen_ortho_reprojs4viz

    subroutine calc_final_rec( params, spproj, projfile, xrec3D, l_postprocess )
        class(parameters),     intent(in)    :: params
        class(sp_project),     intent(inout) :: spproj
        class(string),         intent(in)    :: projfile
        class(commander_base), intent(inout) :: xrec3D
        logical,               intent(in)    :: l_postprocess
        type(string) :: str_state, vol_name, stkname, vol_pproc, vol_mirr, sigma_star
        type(commander_bootstrap_rec3D) :: xbootstrap_rec3D
        integer      :: state, pop, stkind, ind_in_stk, nptcls, ldim(3), sigma_iter, bootstrap_sigma_iter
        real         :: smpd
        logical      :: l_bootstrap_sigmas
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> RECONSTRUCTION AT ORIGINAL SAMPLING'
        write(logfhandle,'(A)') '>>>'
        call spproj%read(projfile) ! ensure we have the latest project info
        call spproj%map_ptcl_ind2stk_ind('ptcl3D', 1, stkind, ind_in_stk)
        stkname = spproj%os_stk%get_str(stkind, 'stk')
        call find_ldim_nptcls(stkname, ldim, nptcls)
        smpd = spproj%os_stk%get(stkind, 'smpd')
        write(logfhandle,'(A,I0,A,F8.4)') '>>> FINAL RECONSTRUCTION SAMPLING: box=', ldim(1), ' smpd=', smpd
        call prep_final_rec_cline(cline_reconstruct3D, 'reconstruct3D')
        sigma_iter = final_rec_sigma_iter()
        l_bootstrap_sigmas = final_rec_needs_bootstrap_sigmas(sigma_iter)
        if( sigma_iter > 0 .and. .not. l_bootstrap_sigmas )then
            call cline_reconstruct3D%set('which_iter', sigma_iter)
            write(logfhandle,'(A,I0)') '>>> FINAL RECONSTRUCTION SIGMA ITERATION: ', sigma_iter
        endif
        if( l_bootstrap_sigmas )then
            bootstrap_sigma_iter = final_rec_bootstrap_sigma_iter(sigma_iter)
            call prep_final_rec_cline(cline_reconstruct3D, 'bootstrap_rec3D')
            call cline_reconstruct3D%set('which_iter', bootstrap_sigma_iter)
            write(logfhandle,'(A,I0)') '>>> FINAL RECONSTRUCTION BOOTSTRAP SIGMA ITERATION: ', bootstrap_sigma_iter
            call xbootstrap_rec3D%execute(cline_reconstruct3D)
        else
            call xrec3D%execute(cline_reconstruct3D)
        endif
        if( .not. l_postprocess )then
            do state = 1, params%nstates
                vol_name  = refine3D_state_vol_fname(state)
                vol_pproc = add2fbody(vol_name, MRC_EXT, PPROC_SUFFIX)
                if( file_exists(vol_pproc) ) call del_file(vol_pproc)
                vol_mirr = add2fbody(vol_pproc, MRC_EXT, MIRR_SUFFIX)
                if( file_exists(vol_mirr) ) call del_file(vol_mirr)
                call vol_name%kill
                call vol_pproc%kill
                call vol_mirr%kill
            enddo
        endif
        call spproj%read_segment('out', projfile)
        call spproj%read_segment('ptcl3D', projfile)
        do state = 1, params%nstates
            pop = spproj%os_ptcl3D%get_pop(state, 'state')
            if( pop == 0 )cycle     ! empty-state case
            str_state = int2str_pad(state,2)
            vol_name  = refine3D_state_vol_fname(state)
            if( .not. file_exists(vol_name) )cycle
            call spproj%add_vol2os_out(vol_name, smpd, state, 'vol', pop=pop)
            call spproj%add_fsc2os_out(refine3D_fsc_fname(state), state, ldim(1))
        enddo
        call spproj%write_segment_inside('out', projfile)
        call stkname%kill
        call sigma_star%kill

        contains

            integer function final_rec_sigma_iter() result(iter)
                integer :: candidates(4), i
                iter = 0
                if( .not. final_stage_uses_ml_reg() ) return
                candidates = 0
                if( cline_refine3D%defined('endit') )then
                    candidates(1) = cline_refine3D%get_iarg('endit') + 1
                    candidates(3) = cline_refine3D%get_iarg('endit')
                endif
                if( cline_refine3D%defined('which_iter') )then
                    candidates(2) = cline_refine3D%get_iarg('which_iter') + 1
                    candidates(4) = cline_refine3D%get_iarg('which_iter')
                endif
                do i = 1,size(candidates)
                    if( candidates(i) <= 0 )cycle
                    sigma_star = sigma2_star_from_iter(candidates(i))
                    if( file_exists(sigma_star) )then
                        iter = candidates(i)
                        return
                    endif
                enddo
            end function final_rec_sigma_iter

            logical function final_rec_needs_bootstrap_sigmas( sigma_iter ) result( l_bootstrap )
                integer, intent(in) :: sigma_iter
                integer :: reg_box
                l_bootstrap = .false.
                if( .not. final_stage_uses_ml_reg() ) return
                if( sigma_iter <= 0 )then
                    l_bootstrap = .true.
                    write(logfhandle,'(A)') '>>> FINAL RECONSTRUCTION: no compatible sigma file found; bootstrapping sigmas'
                    return
                endif
                reg_box  = params%box_crop
                if( cline_refine3D%defined('box_crop')  ) reg_box  = cline_refine3D%get_iarg('box_crop')
                if( reg_box > 0 .and. reg_box /= ldim(1) )then
                    l_bootstrap = .true.
                    write(logfhandle,'(A,I0,A,I0)') &
                        &'>>> FINAL RECONSTRUCTION: registration/final boxes differ; bootstrapping sigmas: ', &
                        &reg_box, ' -> ', ldim(1)
                endif
            end function final_rec_needs_bootstrap_sigmas

            integer function final_rec_bootstrap_sigma_iter( sigma_iter ) result( iter )
                integer, intent(in) :: sigma_iter
                iter = 1
                if( sigma_iter > 0 )then
                    ! Write bootstrap sigmas to the next index so an existing
                    ! compatible star is never overwritten. In the common
                    ! refine3D-finalized case this advances endit+1 to endit+2.
                    iter = sigma_iter + 1
                else if( cline_refine3D%defined('endit') )then
                    iter = cline_refine3D%get_iarg('endit') + 2
                else if( cline_refine3D%defined('which_iter') )then
                    iter = cline_refine3D%get_iarg('which_iter') + 2
                endif
                iter = max(1, iter)
            end function final_rec_bootstrap_sigma_iter

            subroutine prep_final_rec_cline( child_cline, prg )
                class(cmdline), intent(inout) :: child_cline
                character(len=*), intent(in)  :: prg
                call child_cline%kill
                call child_cline%set('prg',      prg)
                call child_cline%set('mkdir',    'no')
                call child_cline%set('projfile', projfile)
                ! volassemble appends _STATENN and writes the extension-less
                ! resolution document next to rec_final_stateNN.mrc.
                call child_cline%set('outfile', 'RESOLUTION_FINAL.txt')
                call child_cline%set('pgrp',    params%pgrp)
                if( params%nthr    > 1  ) call child_cline%set('nthr',    params%nthr)
                if( params%mskdiam > 0. ) call child_cline%set('mskdiam', params%mskdiam)
                if( params%nparts  > 1  ) call child_cline%set('nparts',  params%nparts)
                if( params%nstates > 1  ) call child_cline%set('nstates', params%nstates)
                if( .not. l_postprocess )then
                    call child_cline%set('postprocess', 'no')
                endif
                if( prg.eq.'reconstruct3D' .and. .not. final_stage_uses_ml_reg() )then
                    call child_cline%set('objfun', 'cc')
                    call child_cline%set('ml_reg', 'no')
                endif
            end subroutine prep_final_rec_cline

            logical function final_stage_uses_ml_reg() result( l_ml_reg )
                l_ml_reg = .false.
                if( .not. cline_refine3D%defined('ml_reg') ) return
                if( cline_refine3D%get_carg('ml_reg').ne.'yes' ) return
                if( cline_refine3D%defined('objfun') )then
                    l_ml_reg = cline_refine3D%get_carg('objfun').eq.'euclid'
                else
                    l_ml_reg = .true.
                endif
            end function final_stage_uses_ml_reg

    end subroutine calc_final_rec

    subroutine write_final_rec_outputs( params, spproj, lp )
        class(parameters), intent(in) :: params
        class(sp_project), intent(in) :: spproj
        real,              intent(in) :: lp
        type(string) :: str_state, vol_name, vol_final, vol_final_lp
        type(string) :: vol_pproc, vol_final_pproc, vol_mirr, vol_final_mirr
        integer :: state
        real    :: lp_snapshot
        do state = 1, params%nstates
            if( .not.spproj%isthere_in_osout('vol', state) )cycle ! empty-state case
            str_state      = int2str_pad(state,2)
            vol_name       = refine3D_state_vol_fname(state) ! reconstruction from particles stored in project
            if( .not. file_exists(vol_name) )cycle
            vol_final      = string(REC_FBODY)//str_state//MRC_EXT
            call simple_copy_file(vol_name, vol_final)
            vol_final_lp = add2fbody(vol_final, MRC_EXT, LP_SUFFIX)
            lp_snapshot = abinitio_state_fsc_lowpass(state, params%box, params%smpd, lp)
            call write_abinitio_lowpass_snapshot(vol_final, lp_snapshot, vol_final_lp, params%smpd)
            vol_pproc = add2fbody(vol_name, MRC_EXT, PPROC_SUFFIX)
            vol_final_pproc = add2fbody(vol_final, MRC_EXT, PPROC_SUFFIX)
            if( file_exists(vol_pproc) )then
                call simple_copy_file(vol_pproc, vol_final_pproc)
                vol_mirr = add2fbody(vol_pproc, MRC_EXT, MIRR_SUFFIX)
                if( file_exists(vol_mirr) )then
                    vol_final_mirr = add2fbody(vol_final_pproc, MRC_EXT, MIRR_SUFFIX)
                    call simple_copy_file(vol_mirr, vol_final_mirr)
                endif
            else
                if( file_exists(vol_final_pproc) ) call del_file(vol_final_pproc)
                vol_final_mirr = add2fbody(vol_final_pproc, MRC_EXT, MIRR_SUFFIX)
                if( file_exists(vol_final_mirr) ) call del_file(vol_final_mirr)
            endif
        enddo
        call vol_final_lp%kill
    end subroutine write_final_rec_outputs

    ! create starting noise volume(s)
    subroutine generate_random_volumes( params, box, smpd, cline )
        class(parameters), intent(inout) :: params
        integer,           intent(in)    :: box
        real,              intent(in)    :: smpd
        type(cmdline),     intent(inout) :: cline
        type(string) :: vol_name
        type(image)  :: noisevol, signal
        real         :: b
        integer      :: s
        ! The starting volume is noisy sphere whose values are scaled
        ! to suit the euclid/sigma2 alignment scheme
        ! random normal sphere N(0,5/box) + normal background N(0,5/box)
        call noisevol%new([box,box,box], smpd)
        b = 5.0/real(box)
        call signal%new([box,box,box], smpd)
        call signal%gauran(0.0, b)
        call signal%mask3D_soft(0.25*real(box), backgr=0.)
        do s = 1, params%nstates
            call noisevol%gauran(0., b)
            call noisevol%add(signal)
            vol_name = refine3D_startvol_fname(s)
            call cline%set('vol'//int2str(s), vol_name)
            params%vols(s) = vol_name
            call noisevol%write(vol_name)
            call noisevol%gauran(0., b)
            call noisevol%add(signal)
            vol_name = refine3D_startvol_half_fname(s, 'even')
            call noisevol%write(vol_name)
            vol_name = refine3D_startvol_half_fname(s, 'even', unfil=.true.)
            call noisevol%write(vol_name)
            call noisevol%gauran(0., b)
            call noisevol%add(signal)
            vol_name = refine3D_startvol_half_fname(s, 'odd')
            call noisevol%write(vol_name)
            vol_name = refine3D_startvol_half_fname(s, 'odd', unfil=.true.)
            call noisevol%write(vol_name)
        end do
        call signal%kill
        call noisevol%kill
    end subroutine generate_random_volumes

    subroutine print_states( params, istage )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: istage
        type(sp_project)     :: spproj
        integer, allocatable :: states(:)
        if( nstates_glob == 1 )return
        if( params%print_states /= 'yes' )return
        call spproj%read_segment('ptcl3D', params%projfile)
        states = nint(spproj%os_ptcl3D%get_all('state'))
        call arr2txtfile(states, string('states_stage'//int2str(istage)//'.txt'))
        deallocate(states)
        call spproj%kill
    end subroutine

end module simple_abinitio_utils
