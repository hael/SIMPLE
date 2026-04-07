!@descr: shared helpers for reading, masking and filtering reference volumes for reprojection
module simple_reproj_refvol_utils
use simple_core_module_api
use simple_builder,    only: builder
use simple_parameters, only: parameters
implicit none

public :: read_mask_filter_refvols, calcrefvolshift_and_mapshifts2ptcls
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
            params%l_filemsk .or. params%l_update_frac )then
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
        logical, allocatable :: l_msk(:,:,:)
        real    :: cur_fil(params%box_crop)
        integer :: filtsz
        logical :: have_even, have_odd
        ! READ
        vol_even = params%vols_even(s)
        vol_odd  = params%vols_odd(s)
        vol_avg  = params%vols(s)
        have_even = file_exists(vol_even)
        have_odd  = file_exists(vol_odd)
        if( have_even .and. have_odd )then
            call build%vol%read_and_crop(   vol_even, params%smpd, params%box_crop, params%smpd_crop)
            call build%vol_odd%read_and_crop(vol_odd,  params%smpd, params%box_crop, params%smpd_crop)
        else
            if( .not. file_exists(vol_avg) )then
                THROW_HARD('No usable reference volume inputs for state='//int2str(s)//'; need vol1 or vol_even/vol_odd')
            endif
            call build%vol%read_and_crop(vol_avg, params%smpd, params%box_crop, params%smpd_crop)
            call build%vol_odd%copy_fast(build%vol)
        endif
        if( s == 1 .and. params%l_filemsk )then
            ! read 3D envelope mask
            call build%mskvol%read_and_crop(params%mskfile, params%smpd, params%box_crop, params%smpd_crop)
        endif
        ! noise regularization
        if( params%l_noise_reg )then
            call build%vol%add_gauran(params%eps)
            call build%vol_odd%add_gauran(params%eps)
        endif
        ! MASK
        if( params%l_filemsk )then
            ! envelope masking
            call build%vol%zero_env_background(build%mskvol)
            call build%vol_odd%zero_env_background(build%mskvol)
            call build%vol%mul(build%mskvol)
            call build%vol_odd%mul(build%mskvol)
        else
            ! circular masking
            call build%vol%mask3D_soft(params%msk_crop, backgr=0.0)
            call build%vol_odd%mask3D_soft(params%msk_crop, backgr=0.0)
        endif
        ! FILTER
        if( params%l_icm )then
            if( params%l_filemsk )then
                l_msk = build%mskvol%bin2logical()
                call build%vol%ICM3D_eo(build%vol_odd, params%lambda, l_msk)
            else
                call build%vol%ICM3D_eo(build%vol_odd, params%lambda)
            endif
            if( params%l_lpset .or. trim(params%combine_eo).eq.'yes' )then ! no independent volume registration, so average eo pairs
                call build%vol%add(build%vol_odd)
                call build%vol%mul(0.5)
                call build%vol_odd%copy_fast(build%vol)
            endif
            ! FT
            call build%vol%fft
            call build%vol_odd%fft
        else if( params%l_lpset )then
            ! read average volume that will occupy both even and odd
            call build%vol%read_and_crop(vol_avg, params%smpd, params%box_crop, params%smpd_crop)
            ! noise regularization
            if( params%l_noise_reg )then
                call build%vol%add_gauran(params%eps)
            endif
            ! mask again, BP filter performed below
            if( params%l_filemsk )then
                ! envelope masking
                call build%vol%zero_env_background(build%mskvol)
                call build%vol%mul(build%mskvol)
            else
                ! circular masking
                call build%vol%mask3D_soft(params%msk_crop, backgr=0.0)
            endif
            ! FT & odd <- even
            call build%vol%fft
            call build%vol_odd%copy_fast(build%vol)
        else
            ! FT
            call build%vol%fft
            call build%vol_odd%fft
        endif
        if( params%l_ml_reg )then
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

end module simple_reproj_refvol_utils