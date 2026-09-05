!@descr: utilities for ab initio 3D reconstruction used by commanders_abinitio
module simple_abinitio_utils
use, intrinsic :: iso_fortran_env, only: int64
use simple_commanders_api
use simple_commanders_rec,       only: commander_bootstrap_rec3D
use simple_commanders_volops,    only: commander_symmetrize_map
use simple_cluster_seed,         only: gen_labelling
use simple_class_frcs,           only: class_frcs
use simple_euclid_sigma2,        only: sigma2_star_from_iter
use simple_matcher_refvol_utils, only: remove_ref_section_files
use simple_parameters,           only: parameters
use simple_refine3D_fnames,      only: refine3D_fsc_fname, refine3D_startvol_fbody, &
    &refine3D_startvol_fname, refine3D_startvol_half_fname, &
    &refine3D_state_halfvol_fname, refine3D_state_vol_fbody, refine3D_state_vol_fname
use simple_sigma2_state,         only: sigma2_state_project_layout_digest, sigma2_state_validate_identity
use simple_sigma2_state_file,    only: sigma2_state_validate_file, SIGMA2_GROUP_GLOBAL, &
    &SIGMA2_GROUP_STACK, SIGMA2_STATE_COMMITTED
use simple_vol_pproc_policy,     only: state_mask_is_compatible
implicit none
#include "simple_local_flags.inc"

! singleton variables
type(lp_crop_inf), allocatable :: lpinfo(:)
logical          :: l_srch4symaxis    = .false., l_symran = .false.
logical          :: l_ini3D           = .false.
logical          :: l_automsk         = .false.
logical          :: l_automsk_off     = .false.
logical          :: l_nonuniform      = .false.
logical          :: l_state_continue_mode = .false.
logical          :: l_refine3D_mode_override   = .false.
logical          :: l_refine3D_lp_override     = .false.
logical          :: l_refine3D_lpstop_override = .false.
logical          :: l_cavgs_mode = .false.
type(sym)        :: se1, se2
type(cmdline)    :: cline_refine3D, cline_symmap, cline_reconstruct3D, cline_reproject
type(string)     :: refine3D_mode_override
real             :: update_frac  = 1.0
integer          :: nstates_glob = 1, nptcls_eff = 0
integer          :: nstages_refine3D = 0
integer, parameter :: FINAL_PCG_MAXITS_FLOOR = 5

! In submodule: simple_abinitio_controller.f90
interface
    module function abinitio_rec_fbody() result(fbody)
        character(len=15) :: fbody
    end function abinitio_rec_fbody

    module function abinitio_lpstop_bounds() result(bounds)
        real :: bounds(2)
    end function abinitio_lpstop_bounds

    module function abinitio_lpstart_bounds() result(bounds)
        real :: bounds(2)
    end function abinitio_lpstart_bounds

    module function abinitio_cenlp_default() result(cenlp)
        real :: cenlp
    end function abinitio_cenlp_default

    module function abinitio_lpsymsrch_lb() result(lp)
        real :: lp
    end function abinitio_lpsymsrch_lb

    module function abinitio_update_frac_max() result(frac)
        real :: frac
    end function abinitio_update_frac_max

    module function abinitio_full_sample_switch_frac() result(frac)
        real :: frac
    end function abinitio_full_sample_switch_frac

    module function abinitio_lpstart_ini3D() result(lp)
        real :: lp
    end function abinitio_lpstart_ini3D

    module function abinitio_lpstop_ini3D() result(lp)
        real :: lp
    end function abinitio_lpstop_ini3D

    module function abinitio_independent_lpstop_default() result(lp)
        real :: lp
    end function abinitio_independent_lpstop_default

    module function abinitio_nstages() result(nstages_out)
        integer :: nstages_out
    end function abinitio_nstages

    module function abinitio_independent_nstages_default() result(nstages_out)
        integer :: nstages_out
    end function abinitio_independent_nstages_default

    module function abinitio_nstages_ini3D() result(nstages_out)
        integer :: nstages_out
    end function abinitio_nstages_ini3D

    module function abinitio_nstages_ini3D_max() result(nstages_out)
        integer :: nstages_out
    end function abinitio_nstages_ini3D_max

    module function abinitio_cavgs_early_nstages() result(nstages_out)
        integer :: nstages_out
    end function abinitio_cavgs_early_nstages

    module function abinitio_symsrch_stage() result(istage)
        integer :: istage
    end function abinitio_symsrch_stage

    module function abinitio_ml_reg_start_stage() result(istage)
        integer :: istage
    end function abinitio_ml_reg_start_stage

    module function abinitio_het_docked_stage() result(istage)
        integer :: istage
    end function abinitio_het_docked_stage

    module function abinitio_docked_cohort_active( params, istage ) result(l_active)
        class(parameters), intent(in) :: params
        integer,           intent(in) :: istage
        logical :: l_active
    end function abinitio_docked_cohort_active

    module function abinitio_stoch_sampl_stage(params) result(istage)
        class(parameters), intent(in) :: params
        integer :: istage
    end function abinitio_stoch_sampl_stage

    module function abinitio_nsample_default() result(nsample)
        integer :: nsample
    end function abinitio_nsample_default

    module function abinitio_remaining_niters(first_stage, last_stage) result(niters)
        integer, intent(in) :: first_stage, last_stage
        integer :: niters
    end function abinitio_remaining_niters

    module subroutine set_cline_refine3D( params, istage, l_cavgs )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: istage
        logical,           intent(in) :: l_cavgs
    end subroutine set_cline_refine3D

    module subroutine calc_docked_multistate_max_sampling( params, nptcls, nptcls_cap, ufrac_cap )
        class(parameters), intent(in)  :: params
        integer,           intent(in)  :: nptcls
        integer,           intent(out) :: nptcls_cap
        real,              intent(out) :: ufrac_cap
    end subroutine calc_docked_multistate_max_sampling

end interface

contains

    integer function abinitio_stage_box_crop( params, istage ) result(box_crop)
        class(parameters), intent(in) :: params
        integer,           intent(in) :: istage
        box_crop = lpinfo(istage)%box_crop
    end function abinitio_stage_box_crop

    real function abinitio_stage_smpd_crop( params, istage ) result(smpd_crop)
        class(parameters), intent(in) :: params
        integer,           intent(in) :: istage
        smpd_crop = real(params%box) / real(abinitio_stage_box_crop(params, istage)) * params%smpd
    end function abinitio_stage_smpd_crop

    subroutine prep_class_command_lines( params, cline, projfile )
        class(parameters), intent(in) :: params
        class(cmdline),    intent(in) :: cline
        class(string),     intent(in) :: projfile
        cline_refine3D      = cline
        cline_symmap        = cline
        cline_reconstruct3D = cline
        cline_reproject     = cline
        l_refine3D_mode_override   = cline%defined('refine')
        l_refine3D_lp_override     = cline%defined('lp')
        l_refine3D_lpstop_override = cline%defined('lpstop')
        if( l_refine3D_mode_override ) refine3D_mode_override = trim(params%refine)
        if( l_refine3D_mode_override ) write(logfhandle,'(A,A)') &
            &'>>> ABINITIO3D REFINE MODE OVERRIDE: ', refine3D_mode_override%to_char()
        ! refine3D
        call cline_refine3D%set('prg',                    'refine3D')
        call cline_refine3D%set('pgrp',                  params%pgrp)
        call cline_refine3D%set('projfile',                 projfile)
        call cline_refine3D%delete('box')
        call cline_refine3D%delete('smpd')
        call cline_refine3D%delete('smpd_crop')
        ! symmetrization
        call cline_symmap%set('prg',                'symmetrize_map')
        call strip_pcg_backend_keys(cline_symmap)
        call cline_symmap%set('pgrp',                    params%pgrp)
        call cline_symmap%set('projfile',                   projfile)
        call cline_symmap%set('center',                        'yes')
        if( .not. cline_symmap%defined('cenlp') )then
        call cline_symmap%set('cenlp',                 abinitio_cenlp_default())
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
        call strip_pcg_backend_keys(cline_reproject)
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
        call child_cline%delete('nu_refine')
        call child_cline%delete('refs')
        call child_cline%delete('refs_even')
        call child_cline%delete('refs_odd')
        call child_cline%delete('box')
        call child_cline%delete('smpd')
        call child_cline%delete('smpd_crop')
    end subroutine strip_refine3D_planning_keys

    ! Remove the PCG backend key together with every PCG-only control. Any
    ! child command line that leaves the PCG backend must go through here:
    ! the reconstruct3D commander hard-errors on a positive pcg_nu_lambda_rel
    ! without rec_backend=pcg (explicit-activation contract, pcg_priors.md
    ! S6.2). This list is the strip-side dual of the PCG subset copied by
    ! apply_refine3D_reconstruction_controls below -- extend both together.
    subroutine strip_pcg_backend_keys( child_cline )
        class(cmdline), intent(inout) :: child_cline
        call child_cline%delete('rec_backend')
        call child_cline%delete('pcgop')
        call child_cline%delete('maxits_pcg')
        call child_cline%delete('rtol')
        call child_cline%delete('pcg_nu_lambda_rel')
        call child_cline%delete('pcg_nu_supp_target')
    end subroutine strip_pcg_backend_keys

    ! Copy only controls that genuinely define the reconstruction performed at
    ! the current refinement stage. refine3D maxits is an outer alignment count
    ! and is deliberately not part of this interface.
    subroutine apply_refine3D_reconstruction_controls( child_cline )
        class(cmdline), intent(inout) :: child_cline
        if( cline_refine3D%defined('rec_backend') )then
            call child_cline%set('rec_backend', cline_refine3D%get_carg('rec_backend'))
        endif
        if( cline_refine3D%defined('maxits_pcg') )then
            call child_cline%set('maxits_pcg', cline_refine3D%get_iarg('maxits_pcg'))
        endif
        if( cline_refine3D%defined('rtol') )then
            call child_cline%set('rtol', cline_refine3D%get_rarg('rtol'))
        endif
        if( cline_refine3D%defined('pcg_nu_lambda_rel') )then
            call child_cline%set('pcg_nu_lambda_rel', cline_refine3D%get_rarg('pcg_nu_lambda_rel'))
        endif
        if( cline_refine3D%defined('pcg_nu_supp_target') )then
            call child_cline%set('pcg_nu_supp_target', cline_refine3D%get_rarg('pcg_nu_supp_target'))
        endif
        if( cline_refine3D%defined('ml_reg') )then
            call child_cline%set('ml_reg', cline_refine3D%get_carg('ml_reg'))
        endif
        if( cline_refine3D%defined('objfun') )then
            call child_cline%set('objfun', cline_refine3D%get_carg('objfun'))
        endif
        if( cline_refine3D%defined('envfsc') )then
            call child_cline%set('envfsc', cline_refine3D%get_carg('envfsc'))
        endif
        if( cline_refine3D%defined('envmsklp') )then
            call child_cline%set('envmsklp', cline_refine3D%get_rarg('envmsklp'))
        endif
        if( cline_refine3D%defined('conical_fsc') )then
            call child_cline%set('conical_fsc', cline_refine3D%get_carg('conical_fsc'))
        endif
        if( cline_refine3D%defined('ptcl_src') )then
            call child_cline%set('ptcl_src', cline_refine3D%get_carg('ptcl_src'))
        endif
        if( cline_refine3D%defined('which_iter') )then
            call child_cline%set('which_iter', cline_refine3D%get_iarg('which_iter'))
        endif
    end subroutine apply_refine3D_reconstruction_controls

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
        smpd = find_img_smpd(vol_name)
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

    !> smpd/box describe the sampling the caller expects; the volume on disk
    !! may be at another box (e.g. the full box after the symmetry search while
    !! the stage runs cropped), in which case the pixel size is rescaled with
    !! the box*smpd invariant so the header stays physically correct.
    subroutine write_abinitio_lowpass_snapshot( vol_in, lp, vol_out, smpd, box )
        class(string),     intent(in) :: vol_in, vol_out
        real,              intent(in) :: lp, smpd
        integer, optional, intent(in) :: box
        type(image) :: vol_lp
        integer :: ldim(3), nptcls
        real    :: lp_eff, smpd_eff
        if( .not. file_exists(vol_in) ) return
        call find_ldim_nptcls(vol_in, ldim, nptcls)
        smpd_eff = smpd
        if( present(box) )then
            if( box > 0 .and. ldim(1) > 0 .and. ldim(1) /= box ) smpd_eff = smpd * real(box) / real(ldim(1))
        endif
        lp_eff = max(2.0 * smpd_eff, lp)
        call vol_lp%new(ldim, smpd_eff)
        call vol_lp%read(vol_in)
        call vol_lp%fft()
        call vol_lp%bp(0., lp_eff)
        call vol_lp%ifft()
        call vol_lp%write(vol_out, del_if_exists=.true.)
        call vol_lp%kill
    end subroutine write_abinitio_lowpass_snapshot

    real function abinitio_state_fsc_lowpass( state, box, smpd, fallback_lp, istage ) result( lp )
        integer,           intent(in) :: state, box
        real,              intent(in) :: smpd, fallback_lp
        integer, optional, intent(in) :: istage
        type(string) :: fsc_name
        real, allocatable :: fsc(:), res(:)
        real :: fsc05, fsc0143
        lp = fallback_lp
        ! In cavgs mode the even/odd class-average FSC is only unreliable for the
        ! first stages (ml_reg off). Once ml_reg is active the FSC resolution is
        ! trusted for filtering; final outputs (no istage) always trust it.
        if( l_cavgs_mode )then
            if( present(istage) )then
                if( istage < abinitio_ml_reg_start_stage() ) return
            endif
        endif
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
        real             :: lpstart_bounds(2), lpstop_bounds(2)
        integer          :: filtsz, nstages
        l_cavgs_mode = l_cavgs
        nstages = active_lp_schedule_nstages()
        if( trim(params%force_lp_range).eq.'yes' )then
            if( .not.(present(lpstart) .and. present(lpstop)) )then
                THROW_HARD('force_lp_range=yes requires both lpstart and lpstop')
            endif
            if( allocated(lpinfo) ) deallocate(lpinfo)
            allocate(lpinfo(nstages))
            call lpstages_fast(params%box, nstages, params%smpd, lpstart, lpstop, lpinfo)
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
        allocate(lpinfo(nstages))
        lpstop_bounds  = abinitio_lpstop_bounds()
        lpstart_bounds = abinitio_lpstart_bounds()
        lpfinal = max(lpstop_bounds(1),calc_lplim_final_stage(3))
        lpfinal = min(lpstop_bounds(2),lpfinal)
        if( present(lpstop) ) lpfinal = max(lpstop,lpfinal)
        if( present(lpstart) )then
            call lpstages(params%box, nstages, frcs_avg, params%smpd,&
            &lpstart, lpstart, lpfinal, lpinfo, l_cavgs, verbose=.true.)
        else
            call lpstages(params%box, nstages, frcs_avg, params%smpd,&
            &lpstart_bounds(1), lpstart_bounds(2), lpfinal, lpinfo, l_cavgs, verbose=.true.)
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

    subroutine set_lplims_from_input( params, spproj, lpstart, lpstop )
        class(parameters), intent(inout) :: params
        class(sp_project), intent(in)    :: spproj
        real,              intent(in)    :: lpstart, lpstop
        integer :: nstages
        l_cavgs_mode = .false.
        nstages = active_lp_schedule_nstages()
        if( allocated(lpinfo) ) deallocate(lpinfo)
        allocate(lpinfo(nstages))
        call lpstages_setlims(params%box, nstages, params%smpd, lpstart, lpstop, lpinfo)
    end subroutine set_lplims_from_input

    integer function active_lp_schedule_nstages() result(nstages)
        nstages = abinitio_nstages()
        if( nstages_refine3D > 0 ) nstages = min(nstages_refine3D, nstages)
        nstages = max(1, nstages)
    end function active_lp_schedule_nstages

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
                lp_snapshot = abinitio_state_fsc_lowpass(state, abinitio_stage_box_crop(params, istage), &
                    &abinitio_stage_smpd_crop(params, istage), lpinfo(istage)%lp, istage)
                call write_abinitio_lowpass_snapshot(vol_name, lp_snapshot, vol_lp_stage, &
                    &abinitio_stage_smpd_crop(params, istage), box=abinitio_stage_box_crop(params, istage))
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
            lpsym = max(abinitio_lpsymsrch_lb(),lpinfo(abinitio_symsrch_stage())%lp)
            write(logfhandle,'(A,F5.1)') '>>> DID SET MAP SYMMETRIZATION LOW-PASS LIMIT (IN A) TO: ', lpsym
            write(logfhandle,'(A)') '>>>'
            if( params%nstates > 1 )then
                write(logfhandle,'(A)') '>>> STATE-WISE MAP SYMMETRIZATION'
            else
                write(logfhandle,'(A)') '>>> MAP SYMMETRIZATION'
            endif
            write(logfhandle,'(A)') '>>>'
            call cline_symmap%set('smpd', abinitio_stage_smpd_crop(params, istage))
            call cline_symmap%set('box',  abinitio_stage_box_crop(params, istage))
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
                cline_symrec = cline_reconstruct3D
                call apply_refine3D_reconstruction_controls(cline_symrec)
                ! the symmetric reconstruction belongs to the stage and runs at
                ! its cropped sampling, like the stage-boundary reconstruction;
                ! the base reconstruct3D command line carries no box keys
                call cline_symrec%set('box_crop',   abinitio_stage_box_crop(params, istage))
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
                lp_snapshot  = abinitio_state_fsc_lowpass(state, abinitio_stage_box_crop(params, istage), &
                    &abinitio_stage_smpd_crop(params, istage), lpinfo(istage)%lp, istage)
                call write_abinitio_lowpass_snapshot(vol_sym, lp_snapshot, vol_lp_stage, &
                    &abinitio_stage_smpd_crop(params, istage), box=abinitio_stage_box_crop(params, istage))
                call inject_refine3D_volume(params, state, vol_sym)
                call vol_stage%kill
                call vol_lp_stage%kill
            enddo
        endif
    end subroutine symmetrize

    ! Performs reconstruction at selected stage boundaries.
    subroutine calc_rec( params, projfile, xrec3D, istage, current_sample_only )
        class(parameters),       intent(inout) :: params
        class(string),           intent(in)    :: projfile
        class(commander_base),   intent(inout) :: xrec3D
        integer,                 intent(in)    :: istage
        logical, optional,       intent(in)    :: current_sample_only
        type(string)      :: vol_even, vol_odd, vol_even_unfil, vol_odd_unfil
        type(string)      :: tmpl, src, dest, dest_main, dest_even, dest_odd, sstate, sstage, pgrp, vol_diag
        type(cmdline)     :: cline_rec
        integer           :: state
        real              :: lp_snapshot
        logical           :: have_even_stage, have_odd_stage, l_current_sample_only, l_seed_trail_chain, l_pcg_rec
        l_current_sample_only = .false.
        if( present(current_sample_only) ) l_current_sample_only = current_sample_only
        ! Seed the trailing accumulator chain only when the stage this boundary
        ! feeds actually trails (its refine3D cline is already configured);
        ! seeding earlier would park stale full-weight alignments in the chain.
        l_seed_trail_chain = .false.
        if( cline_refine3D%defined('trail_rec') )then
            l_seed_trail_chain = cline_refine3D%get_carg('trail_rec') .eq. 'yes'
        endif
        ! Reconstruction
        pgrp = trim(params%pgrp)
        if( istage <= abinitio_symsrch_stage() ) pgrp = trim(params%pgrp_start)
        cline_rec = cline_reconstruct3D
        call apply_refine3D_reconstruction_controls(cline_rec)
        call cline_rec%set('prg',       'reconstruct3D')
        call cline_rec%set('mkdir',     'no')
        call cline_rec%set('projfile',  projfile)
        call cline_rec%set('pgrp',      pgrp)
        call cline_rec%set('box_crop',  abinitio_stage_box_crop(params, istage))
        call cline_rec%set('trail_rec', 'no')
        call cline_rec%delete('sticky_class_sampling')
        if( cline_rec%get_carg('ml_reg').ne.'yes' ) call cline_rec%set('objfun','cc')
        l_pcg_rec = cline_rec%defined('rec_backend')
        if( l_pcg_rec ) l_pcg_rec = cline_rec%get_carg('rec_backend').eq.'pcg'
        do state = 1,params%nstates
            call cline_rec%delete('vol'//int2str(state))
        enddo
        call cline_rec%delete('vol_even')
        call cline_rec%delete('vol_odd')
        if( l_current_sample_only )then
            if( .not. cline_refine3D%defined('update_frac') )then
                THROW_HARD('current-sample reconstruction requires update_frac')
            endif
            call cline_rec%set('update_frac', cline_refine3D%get_rarg('update_frac'))
        else
            call cline_rec%delete('update_frac')
            ! A full stage-boundary reconstruction is the producer of the
            ! trailing accumulator chain: seed it at full-dataset weight so the
            ! consuming trail_rec stage starts from complete blended statistics.
            if( l_seed_trail_chain ) call cline_rec%set('trail_seed', 'yes')
        endif
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
            lp_snapshot = abinitio_state_fsc_lowpass(state, abinitio_stage_box_crop(params, istage), &
                &abinitio_stage_smpd_crop(params, istage), lpinfo(istage)%lp, istage)
            call write_abinitio_lowpass_snapshot(dest_main, lp_snapshot, vol_diag, &
                &abinitio_stage_smpd_crop(params, istage), box=abinitio_stage_box_crop(params, istage))
            vol_even = refine3D_state_halfvol_fname(state, 'even')
            if( file_exists(vol_even) )then
                dest = tmpl//'_even_unfil'//MRC_EXT
                vol_even_unfil = refine3D_state_halfvol_fname(state, 'even', unfil=.true.)
                if( l_pcg_rec .and. file_exists(vol_even_unfil) )then
                    call simple_rename(vol_even_unfil, dest)
                else
                    call simple_copy_file(vol_even, dest)
                endif
                dest = tmpl//'_even'//MRC_EXT
                call simple_rename(vol_even, dest)
                dest_even = dest
                have_even_stage = .true.
            endif
            vol_odd  = refine3D_state_halfvol_fname(state, 'odd')
            if( file_exists(vol_odd) )then
                dest = tmpl//'_odd_unfil'//MRC_EXT
                vol_odd_unfil = refine3D_state_halfvol_fname(state, 'odd', unfil=.true.)
                if( l_pcg_rec .and. file_exists(vol_odd_unfil) )then
                    call simple_rename(vol_odd_unfil, dest)
                else
                    call simple_copy_file(vol_odd, dest)
                endif
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
        call vol_even_unfil%kill
        call vol_odd_unfil%kill
        call cline_rec%kill
    end subroutine calc_rec

    subroutine randomize_states( params, spproj, projfile, xrec3D, istage, clean_sampling, reconstruct_states )
        use simple_commanders_euclid,  only: commander_calc_group_sigmas
        class(parameters),     intent(inout) :: params
        class(sp_project),     intent(inout) :: spproj
        class(string),         intent(in)    :: projfile
        class(commander_base), intent(inout) :: xrec3D
        integer,               intent(in)    :: istage
        logical, optional,     intent(in)    :: clean_sampling, reconstruct_states
        integer, parameter :: MIN_SPLIT_STATE_POP = 5
        type(commander_calc_group_sigmas) :: xcalc_group_sigmas
        type(cmdline)                     :: cline_calc_group_sigmas
        integer :: pop, state
        logical :: l_clean_sampling, l_reconstruct_states
        l_clean_sampling     = .true.
        l_reconstruct_states = .true.
        if( present(clean_sampling)     ) l_clean_sampling     = clean_sampling
        if( present(reconstruct_states) ) l_reconstruct_states = reconstruct_states
        call spproj%read_segment('ptcl3D', projfile)
        if( l_clean_sampling ) call spproj%os_ptcl3D%clean_entry('updatecnt', 'sampled')
        call gen_labelling(spproj%os_ptcl3D, params%nstates, 'uniform')
        do state = 1, params%nstates
            pop = spproj%os_ptcl3D%get_pop(state, 'state')
            if( pop <= MIN_SPLIT_STATE_POP )then
                THROW_HARD('docked split generated insufficient state population; reduce nstates or provide more active particles/classes')
            endif
        enddo
        call spproj%write_segment_inside(params%oritype, projfile)
        call cline_refine3D%set(     'nstates', params%nstates)
        call cline_reconstruct3D%set('nstates', params%nstates)
        call cline_reproject%set(    'nstates', params%nstates)
        ! Legacy reconstruction selects an iteration STAR, so materialize it
        ! before this matcher-bypassing reconstruction. Canonical state is
        ! already committed and state relabelling does not change its row
        ! identity or global/stack grouping; reconstruct3D consumes it directly.
        if( cline_refine3D%get_carg('ml_reg').eq.'yes' )then
            if( params%l_sigma_canonical )then
                write(logfhandle,'(A)') '>>> DOCKED SPLIT: reusing committed canonical sigmas'
            else
                cline_calc_group_sigmas = cline_refine3D
                call cline_calc_group_sigmas%set('prg', 'calc_group_sigmas')
                call strip_pcg_backend_keys(cline_calc_group_sigmas)
                call xcalc_group_sigmas%execute(cline_calc_group_sigmas)
                call cline_calc_group_sigmas%kill
            endif
        endif
        ! Multi-state reconstruction
        if( l_reconstruct_states ) call calc_rec(params, projfile, xrec3D, istage)
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

    subroutine configure_final_pcg_solve_budget( source_cline, final_cline )
        class(cmdline), intent(in)    :: source_cline
        class(cmdline), intent(inout) :: final_cline
        integer :: maxits_final
        maxits_final = FINAL_PCG_MAXITS_FLOOR
        if( source_cline%defined('maxits_pcg') ) &
            &maxits_final = max(FINAL_PCG_MAXITS_FLOOR, source_cline%get_iarg('maxits_pcg'))
        call final_cline%set('maxits_pcg', maxits_final)
    end subroutine configure_final_pcg_solve_budget

    subroutine calc_final_rec( params, spproj, projfile, xrec3D, l_postprocess )
        use simple_commanders_euclid, only: commander_calc_pspec
        class(parameters),     intent(in)    :: params
        class(sp_project),     intent(inout) :: spproj
        class(string),         intent(in)    :: projfile
        class(commander_base), intent(inout) :: xrec3D
        logical,               intent(in)    :: l_postprocess
        type(string) :: str_state, vol_name, stkname, vol_pproc, vol_mirr, sigma_star, vol_envmsk
        type(commander_bootstrap_rec3D) :: xbootstrap_rec3D
        type(commander_calc_pspec) :: xcalc_pspec
        type(cmdline) :: cline_calc_pspec
        integer      :: ldim(3), state, pop, stkind, ind_in_stk, nptcls, sigma_iter, bootstrap_sigma_iter
        real         :: smpd
        logical      :: l_bootstrap_sigmas, l_mask_exists, l_mask_compatible
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
        if( params%l_sigma_canonical .and. final_stage_uses_ml_reg() )then
            if( canonical_final_rec_needs_bootstrap() )then
                cline_calc_pspec = cline_reconstruct3D
                call cline_calc_pspec%set('prg', 'calc_pspec')
                write(logfhandle,'(A)') &
                    &'>>> FINAL RECONSTRUCTION: rebuilding canonical sigmas at original sampling'
                call xcalc_pspec%execute(cline_calc_pspec)
                call spproj%read(projfile)
                call cline_calc_pspec%kill
            else
                write(logfhandle,'(A)') &
                    &'>>> FINAL RECONSTRUCTION: reusing committed canonical sigmas'
            endif
            sigma_iter        = 0
            l_bootstrap_sigmas = .false.
        else
            sigma_iter = final_rec_sigma_iter()
            l_bootstrap_sigmas = final_rec_needs_bootstrap_sigmas(sigma_iter)
        endif
        if( sigma_iter > 0 .and. .not. l_bootstrap_sigmas )then
            call cline_reconstruct3D%set('which_iter', sigma_iter)
            write(logfhandle,'(A,I0)') '>>> FINAL RECONSTRUCTION SIGMA ITERATION: ', sigma_iter
        endif
        if( l_bootstrap_sigmas )then
            bootstrap_sigma_iter = final_rec_bootstrap_sigma_iter(sigma_iter)
            call prep_final_rec_cline(cline_reconstruct3D, 'bootstrap_rec3D')
            call cline_reconstruct3D%set('which_iter', bootstrap_sigma_iter)
            write(logfhandle,'(A,I0)') '>>> FINAL RECONSTRUCTION BOOTSTRAP SIGMA ITERATION: ', bootstrap_sigma_iter
            if( trim(params%rec_backend) == 'pcg' ) write(logfhandle,'(A,I0)') &
                &'>>> FINAL PCG COLD-SOLVE ITERATION BUDGET: ', cline_reconstruct3D%get_iarg('maxits_pcg')
            call xbootstrap_rec3D%execute(cline_reconstruct3D)
        else
            if( trim(params%rec_backend) == 'pcg' ) write(logfhandle,'(A,I0)') &
                &'>>> FINAL PCG COLD-SOLVE ITERATION BUDGET: ', cline_reconstruct3D%get_iarg('maxits_pcg')
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
            if( params%l_envfsc )then
                vol_envmsk = AUTOMASK_FBODY//trim(str_state%to_char())//MRC_EXT
                call state_mask_is_compatible(vol_envmsk, ldim(1), smpd, l_mask_exists, l_mask_compatible)
                if( l_mask_compatible )then
                    call spproj%add_vol2os_out(vol_envmsk, smpd, state, 'vol_msk', ldim(1))
                else if( l_mask_exists )then
                    THROW_WARN('>>> FINAL RECONSTRUCTION: envfsc mask has incompatible dimensions or sampling for state '//str_state%to_char())
                else
                    THROW_WARN('>>> FINAL RECONSTRUCTION: expected envfsc mask file does not exist for state '//str_state%to_char())
                endif
            endif
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

            logical function canonical_final_rec_needs_bootstrap() result( l_bootstrap )
                type(string) :: state_path
                integer(int64) :: layout_digest
                integer :: expected_grouping, expected_ngroups, iptcl, nprojptcls, status
                logical :: found
                character(len=STDLEN) :: message
                l_bootstrap = .true.
                call spproj%get_sigma2_state_path(state_path, found)
                if( .not. found )then
                    write(logfhandle,'(A)') &
                        &'>>> FINAL RECONSTRUCTION: canonical sigma state is not registered'
                    return
                endif
                call sigma2_state_validate_file(state_path%to_char(), status, message, deep=.true.)
                if( status /= 0 )then
                    write(logfhandle,'(A)') &
                        &'>>> FINAL RECONSTRUCTION: rebuilding canonical sigmas: '//trim(message)
                    call state_path%kill
                    return
                endif
                nprojptcls   = spproj%os_ptcl3D%get_noris()
                layout_digest = sigma2_state_project_layout_digest(spproj, spproj%os_ptcl3D)
                if( params%l_sigma_glob )then
                    expected_grouping = SIGMA2_GROUP_GLOBAL
                    expected_ngroups  = 1
                else
                    expected_grouping = SIGMA2_GROUP_STACK
                    expected_ngroups  = 0
                    do iptcl = 1, nprojptcls
                        if( spproj%os_ptcl3D%get_state(iptcl) <= 0 ) cycle
                        expected_ngroups = max(expected_ngroups, spproj%os_ptcl3D%get_int(iptcl, 'stkind'))
                    enddo
                endif
                call sigma2_state_validate_identity(state_path%to_char(), ldim(1), smpd, 1, &
                    &fdim(ldim(1))-1, nprojptcls, layout_digest, status, message, &
                    &expected_state=SIGMA2_STATE_COMMITTED, expected_grouping=expected_grouping, &
                    &expected_ngroups=expected_ngroups)
                l_bootstrap = status /= 0
                if( l_bootstrap ) write(logfhandle,'(A)') &
                    &'>>> FINAL RECONSTRUCTION: rebuilding canonical sigmas: '//trim(message)
                call state_path%kill
            end function canonical_final_rec_needs_bootstrap

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
                call child_cline%set('sigma_store', params%sigma_store)
                call child_cline%set('sigma_est',   params%sigma_est)
                ! volassemble appends _STATENN and writes the extension-less
                ! resolution document next to rec_final_stateNN.mrc.
                call child_cline%set('outfile', 'RESOLUTION_FINAL.txt')
                call child_cline%set('pgrp',    params%pgrp)
                call child_cline%set('ptcl_src', params%ptcl_src)
                call child_cline%set('envfsc',   params%envfsc)
                call child_cline%set('envmsklp', params%envmsklp)
                call child_cline%set('binwidth', params%binwidth)
                if( params%nthr    > 1  ) call child_cline%set('nthr',    params%nthr)
                if( params%mskdiam > 0. ) call child_cline%set('mskdiam', params%mskdiam)
                if( params%nparts  > 1  ) call child_cline%set('nparts',  params%nparts)
                if( params%nstates > 1  ) call child_cline%set('nstates', params%nstates)
                if( final_stage_uses_ml_reg() ) call child_cline%set('conical_fsc', params%conical_fsc)
                if( .not. l_postprocess )then
                    call child_cline%set('postprocess', 'no')
                endif
                if( prg.eq.'reconstruct3D' .and. .not. final_stage_uses_ml_reg() )then
                    call child_cline%set('objfun', 'cc')
                    call child_cline%set('ml_reg', 'no')
                endif
                ! the final reconstruction runs on the refinement's backend;
                ! with the euclid ML replay active the Q_NU prior regularizes
                ! it in-solve. When bootstrap_rec3D changes the reconstruction
                ! grid, it retains the learned suppression target but measures
                ! and corrects the Q_NU strength on that grid before shipping
                ! the final map. Explicit controls remain pinned.
                if( trim(params%rec_backend) == 'pcg' )then
                    call child_cline%set('rec_backend', 'pcg')
                    call configure_final_pcg_solve_budget(cline_refine3D, child_cline)
                    if( cline_refine3D%defined('rtol') )&
                        &call child_cline%set('rtol', cline_refine3D%get_rarg('rtol'))
                    if( final_stage_uses_ml_reg() )then
                        call child_cline%set('filt_mode', 'nonuniform')
                        if( cline_refine3D%defined('pcg_nu_lambda_rel') )&
                            &call child_cline%set('pcg_nu_lambda_rel', cline_refine3D%get_rarg('pcg_nu_lambda_rel'))
                        if( cline_refine3D%defined('pcg_nu_supp_target') )&
                            &call child_cline%set('pcg_nu_supp_target', cline_refine3D%get_rarg('pcg_nu_supp_target'))
                    endif
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
            vol_final      = string(abinitio_rec_fbody())//str_state//MRC_EXT
            call simple_copy_file(vol_name, vol_final)
            vol_final_lp = add2fbody(vol_final, MRC_EXT, LP_SUFFIX)
            lp_snapshot = abinitio_state_fsc_lowpass(state, params%box, params%smpd, lp)
            call write_abinitio_lowpass_snapshot(vol_final, lp_snapshot, vol_final_lp, params%smpd, box=params%box)
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
        ! to suit the euclid/sigma2 alignment scheme. References are no
        ! longer multiplied by the original box on projection (the legacy
        ! division/multiplication pair is retired), so the volume itself
        ! must carry the data-quotient scale: the legacy N(0,5/box) noise
        ! times the retired projection factor (the original box size)
        ! random normal sphere + normal background, N(0, 5*box_orig/box)
        call noisevol%new([box,box,box], smpd)
        b = 5.0*real(params%box)/real(box)
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

    subroutine normalize_input_volumes( params, cline )
        class(parameters), intent(in)    :: params
        type(cmdline),     intent(inout) :: cline
        type(string) :: vol_name
        type(image)  :: vol
        real         :: smpd, msk, ave, stdev, maxv, minv, v
        integer      :: ldim(3), s, n
        do s = 1, params%nstates
            if( .not. file_exists(params%vols(s)) )then
                THROW_HARD('Input volume for state '//int2str(s)//' does not exist: '//params%vols(s)%to_char())
            endif
            call find_ldim_nptcls(params%vols(s), ldim, n)
            smpd = find_img_smpd(params%vols(s))
            msk  = 0.5 * params%mskdiam / smpd
            call vol%new(ldim, smpd)
            call vol%read(params%vols(s))
            ! normalization of inner mask region to the data-quotient
            ! reference scale (was 1/box when projection multiplied by box)
            call vol%stats('foreground', ave, stdev, maxv, minv, msk=msk)
            v = stdev*real(ldim(1))/real(params%box)
            call vol%norm_ext(ave, v)
            vol_name = refine3D_startvol_fname(s)
            call vol%write(vol_name)
            call cline%set('vol'//int2str(s), vol_name)
        enddo
        call vol%kill
        call vol_name%kill
    end subroutine normalize_input_volumes

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
