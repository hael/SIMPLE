!@descr: shared helpers for reading, masking, filtering and reprojecting reference volumes
module simple_matcher_refvol_utils
use simple_core_module_api
use simple_builder,    only: builder
use simple_parameters, only: parameters
use simple_cmdline,    only: cmdline
use simple_image,      only: image
use simple_pftc_srch_api, only: polarft_dims_from_file_header
use simple_matcher_smpl_and_lplims, only: enforce_3D_pftc_k_range
implicit none

public :: read_mask_filter_reproject_refvols
public :: pick_lp_est_state, estimate_lp_from_refs
public :: any_volume_source_defined, complete_volume_source_defined
public :: complete_volume_source_available
public :: polar_ref_sections_available
public :: ensure_polar_refs_on_disk, materialize_polar_refs_from_volume_source
public :: remove_polar_ref_section_files, write_polar_refs_from_current_pftc
private
#include "simple_local_flags.inc"

contains

    logical function any_volume_source_defined( cline, nstates ) result( l_defined )
        class(cmdline), intent(in) :: cline
        integer,        intent(in) :: nstates
        integer :: state
        l_defined = .false.
        do state = 1,nstates
            if( cline%defined('vol'//int2str(state)) )then
                l_defined = .true.
                return
            endif
        enddo
    end function any_volume_source_defined

    logical function complete_volume_source_defined( cline, nstates ) result( l_defined )
        class(cmdline), intent(in) :: cline
        integer,        intent(in) :: nstates
        integer :: state
        l_defined = .true.
        do state = 1,nstates
            if( .not. cline%defined('vol'//int2str(state)) )then
                l_defined = .false.
                return
            endif
        enddo
    end function complete_volume_source_defined

    logical function complete_volume_source_available( params ) result( l_available )
        class(parameters), intent(in) :: params
        integer :: state
        l_available = .true.
        do state = 1,params%nstates
            if( .not. file_exists(params%vols(state)) )then
                l_available = .false.
                return
            endif
        enddo
    end function complete_volume_source_available

    logical function polar_ref_sections_available( params )
        class(parameters), intent(in) :: params
        type(string) :: refs, refs_even, refs_odd
        logical      :: have_refs, have_even, have_odd
        refs      = POLAR_REFS_FBODY//BIN_EXT
        refs_even = POLAR_REFS_FBODY//'_even'//BIN_EXT
        refs_odd  = POLAR_REFS_FBODY//'_odd'//BIN_EXT
        have_refs = file_exists(refs)
        have_even = file_exists(refs_even)
        have_odd  = file_exists(refs_odd)
        if( have_even .neqv. have_odd )then
            polar_ref_sections_available = .false.
        elseif( have_refs .and. have_even )then
            polar_ref_sections_available = polar_ref_file_compatible(refs) .and. &
                &polar_ref_file_compatible(refs_even) .and. &
                &polar_ref_file_compatible(refs_odd)
        elseif( have_even )then
            polar_ref_sections_available = polar_ref_file_compatible(refs_even) .and. &
                &polar_ref_file_compatible(refs_odd)
        elseif( have_refs )then
            polar_ref_sections_available = polar_ref_file_compatible(refs)
        else
            polar_ref_sections_available = .false.
        endif
        call refs%kill
        call refs_even%kill
        call refs_odd%kill

    contains

        logical function polar_ref_file_compatible( fname )
            type(string), intent(in) :: fname
            integer :: pftsz_here, kfromto_here(2), nrefs_here
            integer :: expected_nrefs, expected_interpklim
            if( .not. file_exists(fname) )then
                polar_ref_file_compatible = .false.
                return
            endif
            call polarft_dims_from_file_header(fname, pftsz_here, kfromto_here, nrefs_here)
            expected_nrefs      = params%nspace * params%nstates
            expected_interpklim = fdim(params%box_crop) - 1
            ! Header-only gate for stale POLAR_REFS*. Lower-k range
            ! differences are accepted because the reader owns overlap
            ! transfer and k-range zero padding.
            polar_ref_file_compatible = nrefs_here == expected_nrefs
            polar_ref_file_compatible = polar_ref_file_compatible .and. pftsz_here == params%pftsz
            polar_ref_file_compatible = polar_ref_file_compatible .and. kfromto_here(2) <= expected_interpklim
        end function polar_ref_file_compatible
    end function polar_ref_sections_available

    subroutine remove_polar_ref_section_files
        if( file_exists(POLAR_REFS_FBODY//BIN_EXT) ) call del_file(POLAR_REFS_FBODY//BIN_EXT)
        if( file_exists(POLAR_REFS_FBODY//'_even'//BIN_EXT) ) call del_file(POLAR_REFS_FBODY//'_even'//BIN_EXT)
        if( file_exists(POLAR_REFS_FBODY//'_odd'//BIN_EXT) ) call del_file(POLAR_REFS_FBODY//'_odd'//BIN_EXT)
    end subroutine remove_polar_ref_section_files

    subroutine write_polar_refs_from_current_pftc( params, build, nspace_refs, skip_distr_nonroot )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer, optional, intent(in)    :: nspace_refs
        logical, optional, intent(in)    :: skip_distr_nonroot
        integer :: nspace_write, nrefs_write
        if( present(skip_distr_nonroot) )then
            if( skip_distr_nonroot .and. params%l_distr_worker .and. params%part /= 1 ) return
        endif
        nspace_write = params%nspace
        if( present(nspace_refs) ) nspace_write = nspace_refs
        nrefs_write = nspace_write * params%nstates
        call remove_polar_ref_section_files
        call build%pftc%polar_cavger_new(.true., nrefs=nrefs_write)
        call build%pftc%polar_cavger_write_eo_pftcrefs(string(POLAR_REFS_FBODY))
    end subroutine write_polar_refs_from_current_pftc

    subroutine materialize_polar_refs_from_volume_source( params, build, cline, batchsz, cleanup_pftc )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(in)    :: cline
        integer,           intent(in)    :: batchsz
        logical, optional, intent(in)    :: cleanup_pftc
        logical :: l_cleanup
        if( .not. complete_volume_source_available(params) )then
            THROW_HARD('cannot materialize POLAR_REFS without complete volume source files')
        endif
        call read_mask_filter_reproject_refvols(params, build, cline, batchsz)
        call write_polar_refs_from_current_pftc(params, build)
        l_cleanup = .false.
        if( present(cleanup_pftc) ) l_cleanup = cleanup_pftc
        if( l_cleanup )then
            call build%pftc%polar_cavger_kill
            call build%pftc%kill
        endif
    end subroutine materialize_polar_refs_from_volume_source

    subroutine ensure_polar_refs_on_disk( params, build, cline, batchsz, caller )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(in)    :: cline
        integer,           intent(in)    :: batchsz
        character(len=*),  intent(in)    :: caller
        if( any_volume_source_defined(cline, params%nstates) &
            &.and. (.not. complete_volume_source_defined(cline, params%nstates)) )then
            THROW_HARD('incomplete multi-state volume source; provide vol1..volN or use POLAR_REFS; '//caller)
        endif
        if( complete_volume_source_defined(cline, params%nstates) )then
            call materialize_polar_refs_from_volume_source(params, build, cline, batchsz, cleanup_pftc=.true.)
            return
        endif
        if( polar_ref_sections_available(params) ) return
        if( complete_volume_source_available(params) )then
            call materialize_polar_refs_from_volume_source(params, build, cline, batchsz, cleanup_pftc=.true.)
            return
        endif
        THROW_HARD('polar reference sections are missing and no complete volume source is available; '//caller)
    end subroutine ensure_polar_refs_on_disk

    !>  \brief  determines the reference volume shift and map shifts back to particles
    !>          reference volume shifting is performed in shift_and_mask_refvol
    subroutine calcrefvolshift_and_mapshifts2ptcls( params, build, s, volfname, do_center, xyz, map_shift )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: s
        class(string),     intent(in)    :: volfname
        logical,           intent(out)   :: do_center
        real,              intent(out)   :: xyz(3)
        logical,           intent(in)    :: map_shift
        real    :: crop_factor
        logical :: has_been_searched
        do_center = .true.
        ! centering
        if( params%center .eq. 'no' .or. params%nstates > 1 .or. &
            .not. params%l_doshift .or. params%pgrp(:1) .ne. 'c' .or. &
            params%l_update_frac )then
            do_center = .false.
            xyz       = 0.
            return
        endif
        ! taking care of volume dimensions
        call build%vol%read_and_crop(volfname, params%smpd, params%box_crop, params%smpd_crop)
        ! offset
        xyz = build%vol%calc_shiftcen(params%cenlp,params%msk_crop)
        if( params%pgrp .ne. 'c1' ) xyz(1:2) = 0.     ! shifts only along z-axis for C2 and above
        if( arg(xyz) <= CENTHRESH )then
            do_center = .false.
            xyz       = 0.
            return
        endif
        if( map_shift )then
            ! map back to particle oritentations
            has_been_searched = .not.build%spproj%is_virgin_field(params%oritype)
            if( has_been_searched )then
                crop_factor = real(params%box) / real(params%box_crop)
                call build%spproj_field%map3dshift22d(-xyz(:)*crop_factor, state=s)
            endif
        endif
    end subroutine calcrefvolshift_and_mapshifts2ptcls

    subroutine read_mask_filter_refvols( params, build, s )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: s
        type(string)         :: vol_even, vol_odd, vol_avg
        real    :: cur_fil(params%box_crop)
        integer :: filtsz
        logical :: have_even, have_odd, l_nonuniform_mode
        l_nonuniform_mode = trim(params%filt_mode).eq.'nonuniform'
        vol_avg = params%vols(s)
        ! READ: try nonuniform refs first if requested, then regular refs, then average
        if( l_nonuniform_mode )then
            vol_even = add2fbody(params%vols_even(s), params%ext, NUFILT_SUFFIX)
            vol_odd  = add2fbody(params%vols_odd(s),  params%ext, NUFILT_SUFFIX)
        else
            vol_even = params%vols_even(s)
            vol_odd  = params%vols_odd(s)
        endif
        have_even = file_exists(vol_even)
        have_odd  = file_exists(vol_odd)
        ! fall back to regular refs if nonuniform versions don't exist yet (will be created later in volassemble)
        if( .not. have_even .or. .not. have_odd )then
            if( l_nonuniform_mode )then
                vol_even = params%vols_even(s)
                vol_odd  = params%vols_odd(s)
                have_even = file_exists(vol_even)
                have_odd  = file_exists(vol_odd)
            endif
        endif
        ! read even/odd pair or fall back to average volume
        if( have_even .and. have_odd )then
            call build%vol%read_and_crop(   vol_even, params%smpd, params%box_crop, params%smpd_crop)
            call build%vol_odd%read_and_crop(vol_odd,  params%smpd, params%box_crop, params%smpd_crop)
        else
            vol_avg = params%vols(s)
            if( .not. file_exists(vol_avg) )then
                THROW_HARD('No usable reference volume inputs for state='//int2str(s)//'; need volN or vol_even/vol_odd')
            endif
            call build%vol%read_and_crop(vol_avg, params%smpd, params%box_crop, params%smpd_crop)
            call build%vol_odd%copy_fast(build%vol)
        endif
        ! noise regularization
        if( params%l_noise_reg )then
            call build%vol%add_gauran(params%eps)
            call build%vol_odd%add_gauran(params%eps)
        endif
        ! MASK: use circular masking (volassemble handles automask if needed)
        call build%vol%mask3D_soft(params%msk_crop, backgr=0.0)
        call build%vol_odd%mask3D_soft(params%msk_crop, backgr=0.0)
        ! FILTER
        if( params%l_icm )then
            call build%vol%ICM3D_eo(build%vol_odd, params%lambda)
            if( params%l_lpset .or. trim(params%combine_eo).eq.'yes' )then ! no independent volume registration, so average eo pairs
                call build%vol%add(build%vol_odd)
                call build%vol%mul(0.5)
                call build%vol_odd%copy_fast(build%vol)
            endif
            ! FT
            call build%vol%fft
            call build%vol_odd%fft
        else if( params%l_lpset .and. .not. l_nonuniform_mode )then
            ! read average volume that will occupy both even and odd
            if( .not. file_exists(vol_avg) )then
                THROW_HARD('Missing average reference for state='//int2str(s)//' : '//vol_avg%to_char())
            endif
            call build%vol%read_and_crop(vol_avg, params%smpd, params%box_crop, params%smpd_crop)
            ! noise regularization
            if( params%l_noise_reg )then
                call build%vol%add_gauran(params%eps)
            endif
            ! mask with circular mask (volassemble handles automask if needed)
            call build%vol%mask3D_soft(params%msk_crop, backgr=0.0)
            ! FT & odd <- even
            call build%vol%fft
            call build%vol_odd%copy_fast(build%vol)
        else
            ! FT
            call build%vol%fft
            call build%vol_odd%fft
        endif
        if( params%l_ml_reg .or. l_nonuniform_mode )then
            ! filtering done when volumes are assembled
        else if( params%l_icm )then
            ! filtering done above
        else if( params%l_lpset )then
            if( trim(params%gauref)=='yes' )then
                ! Gaussian filter
                call build%vol%bpgau3D(0., params%gaufreq)
                call build%vol_odd%bpgau3D(0., params%gaufreq)
            else
                ! Cosine low-pass filter, works best for nanoparticles
                call build%vol%bp(0., params%lp)
                call build%vol_odd%bp(0., params%lp)
            endif
        else
            ! Optimal filter
            filtsz = build%vol%get_filtsz()
            if( any(build%fsc(s,:) > 0.143) )then
                call fsc2optlp_sub(filtsz,build%fsc(s,:),cur_fil)
                call build%vol%apply_filter(cur_fil)
                call build%vol_odd%apply_filter(cur_fil)
            endif
        endif
    end subroutine read_mask_filter_refvols

    subroutine pick_lp_est_state( params, build, state )
        class(parameters), intent(in)  :: params
        class(builder),    intent(in)  :: build
        integer,           intent(out) :: state
        integer :: s
        real    :: res_avg(params%nstates), res_val
        integer :: n_particles
        logical, allocatable :: mask(:)
        integer, allocatable :: states(:)
        res_avg = huge(res_avg(1))
        states = build%spproj_field%get_all_asint('state')
        do s = 1, params%nstates
            ! Create mask for particles in this state
            mask = states == s
            n_particles = count(mask)
            if( n_particles > 0 )then
                ! Get average res for this state
                res_avg(s) = sum(build%spproj_field%get_all('res'), mask=mask) / real(n_particles)
            endif
        end do
        deallocate(mask,states)
        ! Select state with best (minimum) resolution
        state = minloc(res_avg, dim=1)
    end subroutine pick_lp_est_state

    subroutine estimate_lp_from_refs( params, build, cline, lpstart, lpstop, state )
        use simple_opt_filter,   only: estimate_lplim3D
        use simple_polarft_calc, only: polarft_estimate_lplim3D
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(in)    :: cline
        real,              intent(in)    :: lpstart, lpstop
        integer,           intent(in)    :: state
        type(string)       :: vol_even, vol_odd
        type(image) :: mskvol
        integer     :: npix
        real        :: lpopt
        if( lpstart < lpstop )then
            THROW_HARD('Invalid low-pass range ordering: lpstart must be >= lpstop')
        endif
        if( params%l_polar .and. (.not. complete_volume_source_defined(cline, params%nstates)) )then
            call polarft_estimate_lplim3D(params%box_crop, params%smpd_crop, lpstart, lpstop, lpopt)
        else
            ! Use circular mask for low-pass estimation (volassemble handles automasking)
            call mskvol%disc([params%box_crop,  params%box_crop, params%box_crop], params%smpd_crop,&
                &params%msk_crop, npix )
            vol_even = params%vols_even(state)
            vol_odd  = params%vols_odd(state)
            call build%vol%read_and_crop(    vol_even, params%smpd, params%box_crop, params%smpd_crop)
            call build%vol_odd%read_and_crop(vol_odd,  params%smpd, params%box_crop, params%smpd_crop)
            call estimate_lplim3D(build%vol_odd, build%vol, mskvol, lpstart, lpstop, lpopt)
        endif
        call build%spproj_field%set_all2single('lp_est', lpopt)
        if( params%l_lpauto )then
            params%lp = lpopt
            params%kfromto(2) = calc_fourier_index(params%lp, params%box_crop, params%smpd_crop)
            call build%spproj_field%set_all2single('lp',params%lp)
        endif
        call mskvol%kill
    end subroutine estimate_lp_from_refs

    !> Producer-side reference-section reprojection used by Cartesian assembly
    !> and polar bootstrap. It cannot assume matcher setup already reconciled
    !> the requested band limit with the cropped PFTC interpolation range.
    subroutine read_mask_filter_reproject_refvols( params, build, cline, batchsz )
        use simple_polarft_calc, only: vol_pad2ref_pfts
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(in)    :: cline
        integer,           intent(in)    :: batchsz
        real      :: xyz(3)
        integer   :: s, nrefs, state
        logical   :: do_center
        call pick_lp_est_state(params, build, state)
        if( cline%defined('lpstart') .and. cline%defined('lpstop') )then
            call estimate_lp_from_refs(params, build, cline, params%lpstart, params%lpstop, state)
        endif
        if( build%eulspace%get_noris() /= params%nspace )then
            call build%eulspace%kill
            call build%eulspace%new(params%nspace, is_ptcl=.false.)
            call build%pgrpsyms%build_refspiral(build%eulspace)
        endif
        call enforce_3D_pftc_k_range(params, build, 'read_mask_filter_reproject_refvols')
        nrefs = params%nspace * params%nstates
        call build%pftc%new(params, nrefs, [1, 1], params%kfromto)
        do s = 1, params%nstates
            if( params%l_prob_align_mode )then
                call calcrefvolshift_and_mapshifts2ptcls(params, build, s, params%vols(s), &
                    & do_center, xyz, map_shift=l_distr_worker_glob)
            else
                call calcrefvolshift_and_mapshifts2ptcls(params, build, s, params%vols(s), &
                    & do_center, xyz, map_shift=.true.)
            endif
            call read_mask_filter_refvols(params, build, s)
            call build%vol_pad%new([params%box_croppd, params%box_croppd, params%box_croppd], &
                & params%smpd_crop, wthreads=.true.)
            if( do_center )then
                call build%vol%fft()
                call build%vol%shift(xyz)
            endif
            call build%vol%ifft()
            call build%vol%pad_fft(build%vol_pad)
            call build%vol_pad%expand_cmat(params%box)
            call vol_pad2ref_pfts(build%pftc, build%vol_pad, build%eulspace, s, .true.)
            call build%vol_pad%kill
            call build%vol_pad%kill_expanded
            call build%vol_odd_pad%new([params%box_croppd, params%box_croppd, params%box_croppd], &
                & params%smpd_crop, wthreads=.true.)
            if( do_center )then
                call build%vol_odd%fft()
                call build%vol_odd%shift(xyz)
            endif
            call build%vol_odd%ifft()
            call build%vol_odd%pad_fft(build%vol_odd_pad)
            call build%vol_odd_pad%expand_cmat(params%box)
            call vol_pad2ref_pfts(build%pftc, build%vol_odd_pad, build%eulspace, s, .false.)
            call build%vol_odd_pad%kill
            call build%vol_odd_pad%kill_expanded
        end do
    end subroutine read_mask_filter_reproject_refvols

end module simple_matcher_refvol_utils
