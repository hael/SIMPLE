!@descr: shared helpers for reading, masking and filtering reference volumes for matcher reprojection
module simple_matcher_refvol_utils
use simple_core_module_api
use simple_builder,    only: builder
use simple_parameters, only: parameters
use simple_cmdline,    only: cmdline
use simple_image,      only: image
implicit none

public :: read_mask_filter_reproject_refvols
public :: report_resolution, estimate_lp_from_refs
private
#include "simple_local_flags.inc"

contains

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
                THROW_HARD('No usable reference volume inputs for state='//int2str(s)//'; need vol1 or vol_even/vol_odd')
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

    subroutine report_resolution( params, build, state )
        class(parameters), intent(in)  :: params
        class(builder),    intent(in)  :: build
        integer,           intent(out) :: state
        real :: res_fsc05, res_fsc0143
        real, allocatable  :: res(:)
        integer :: s
        real    :: lpest(params%nstates)
        res = get_resarr(params%box_crop, params%smpd_crop)
        do s = 1, params%nstates
            call get_resolution_at_fsc(build%fsc(s,:), res, 0.5, lpest(s))
        end do
        state = minloc(lpest,dim=1)
        call get_resolution_at_fsc(build%fsc(state,:), res, 0.50,  res_fsc05)
        call get_resolution_at_fsc(build%fsc(state,:), res, 0.143, res_fsc0143)
        call build%spproj_field%set_all2single('res', res_fsc0143)
    end subroutine report_resolution

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
        if( params%l_polar .and. (.not.cline%defined('vol1')) )then
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

    subroutine read_mask_filter_reproject_refvols( params, build, cline, batchsz )
        use simple_polarft_calc, only: vol_pad2ref_pfts
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(in)    :: cline
        integer,           intent(in)    :: batchsz
        real      :: xyz(3)
        integer   :: s, nrefs, state
        logical   :: do_center
        call report_resolution(params, build, state)
        if( cline%defined('lpstart') .and. cline%defined('lpstop') )then
            call estimate_lp_from_refs(params, build, cline, params%lpstart, params%lpstop, state)
        endif
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
            call vol_pad2ref_pfts(build%pftc, build%vol_pad, build%eulspace, s, .true., build%l_resmsk)
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
            call vol_pad2ref_pfts(build%pftc, build%vol_odd_pad, build%eulspace, s, .false., build%l_resmsk)
            call build%vol_odd_pad%kill
            call build%vol_odd_pad%kill_expanded
        end do
    end subroutine read_mask_filter_reproject_refvols

end module simple_matcher_refvol_utils