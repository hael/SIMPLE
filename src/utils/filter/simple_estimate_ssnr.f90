!@descr: spectral signal-to-noise ratio estimation routines
module simple_estimate_ssnr
!$ use omp_lib
!$ use omp_lib_kinds
use simple_defs
use simple_type_defs
use simple_defs_stream
use simple_defs_conv
use simple_srch_sort_loc
use simple_syslib
use simple_math_ft
use simple_string_utils
use simple_fileio
implicit none

public :: fsc2ssnr, fsc2wiener_regularizer, fsc2optlp, fsc2optlp_sub, fsc2cref, ssnr2fsc, ssnr2optlp
public :: lowpass_from_klim, gaussian_filter, bfac_cap_filter
public :: mskdiam2lplimits, mskdiam2streamresthreshold
public :: calc_dose_weights, get_resolution, get_resolution_at_fsc
public :: lpstages_fast, lpstages_setlims, stage_trslim
private
#include "simple_local_flags.inc"
real,    parameter :: LPSTAGES_STAGE1_LP_FLOOR = 20.
integer, parameter :: STAGE_MINBOX             = 88 ! smallest cropped box of a stage

contains

    !> \brief  converts the FSC to SSNR (the 2.* is because of the division of the data)
    function fsc2ssnr( corrs ) result( ssnr )
        real, intent(in)  :: corrs(:) !< FSC
        real, allocatable :: ssnr(:)  !< SSNR
        integer :: nyq, k
        real    :: fsc
        nyq = size(corrs)
        allocate( ssnr(nyq) )
        do k=1,nyq
            fsc     = min(abs(corrs(k)), 0.999)
            ssnr(k) = (2. * fsc) / (1. - fsc)
        end do
    end function fsc2ssnr

    !> Convert a merged-map FSC to the per-observation N/S term used by a
    !! Wiener CTF correction. nobs is the number of independent original
    !! observations represented by the merged FSC, not the number per half-map.
    function fsc2wiener_regularizer( corrs, nobs ) result( noise_to_signal )
        real,    intent(in) :: corrs(:)
        integer, intent(in) :: nobs
        real, allocatable :: noise_to_signal(:), ssnr(:)
        if( nobs < 1 ) THROW_HARD('nobs must be positive; fsc2wiener_regularizer')
        ssnr = fsc2ssnr(corrs)
        allocate(noise_to_signal(size(ssnr)), source=1. / TINY)
        where( ssnr > TINY ) noise_to_signal = real(nobs) / ssnr
    end function fsc2wiener_regularizer

    !> \brief  converts the FSC to the optimal low-pass filter
    function fsc2optlp( corrs ) result( filt )
        real, intent(in)  :: corrs(:) !< fsc plot (correlations)
        real, allocatable :: filt(:)  !< output filter coefficients
        integer :: nyq
        nyq = size(corrs)
        allocate( filt(nyq) )
        filt = 0.
        where( corrs > 0. )     filt = 2. * corrs / (corrs + 1.)
        where( filt  > 0.99999 ) filt = 0.99999
    end function fsc2optlp

    !> \brief  converts the FSC to the Rosenthal & Henderson (2003) figure-of-merit
    !!         weight C_ref = sqrt(2FSC/(1+FSC)), with RELION's truncation
    !!         (postprocessing.cpp applyFscWeighting): zero where FSC <= 0 and
    !!         for every shell from the first one where FSC < 1e-4 onwards.
    !!         The square root of fsc2optlp: the correlation-optimal weight of a
    !!         merged unregularized map, milder than the Wiener form (0.5 against
    !!         0.25 at FSC=0.143)
    function fsc2cref( corrs ) result( filt )
        real, intent(in)  :: corrs(:) !< fsc plot (correlations)
        real, allocatable :: filt(:)  !< output weights
        integer :: k
        filt = sqrt(fsc2optlp(corrs))
        do k = 1,size(corrs)
            if( corrs(k) < 1.e-4 )then
                filt(k:) = 0.
                exit
            endif
        end do
    end function fsc2cref


    !> \brief  converts the FSC to the optimal low-pass filter
    subroutine fsc2optlp_sub( filtsz, corrs, filt, merged )
        integer, intent(in)  :: filtsz        !< sz of filter
        real,    intent(in)  :: corrs(filtsz) !< fsc plot (correlations)
        real,    intent(out) :: filt(filtsz)  !< output filter coefficients
        logical, optional, intent(in) :: merged
        logical :: l_merged
        l_merged = .false.
        if( present(merged) ) l_merged = merged
        filt = 0.
        if( l_merged )then
            where( corrs > 0. ) filt = corrs
        else
            where( corrs > 0. ) filt = 2. * corrs / (corrs + 1.)
        endif
        where( filt  > 0.99999 ) filt = 0.99999
    end subroutine fsc2optlp_sub

    !> \brief  converts the SSNR to FSC
    function ssnr2fsc( ssnr ) result( corrs )
        real, intent(in)  :: ssnr(:)  !< input SSNR array
        real, allocatable :: corrs(:) !< output FSC result
        integer :: nyq, k
        nyq = size(ssnr)
        allocate( corrs(nyq) )
        do k=1,nyq
            corrs(k) = ssnr(k) / (ssnr(k) + 1.)
        end do
    end function ssnr2fsc

    ! !> \brief  converts the SSNR 2 the optimal low-pass filter
    function ssnr2optlp( ssnr ) result( w )
        real, intent(in)  :: ssnr(:) !<  instrument SSNR
        real, allocatable :: w(:) !<  FIR low-pass filter
        integer :: nyq, k
        nyq = size(ssnr)
        allocate( w(nyq) )
        do k=1,nyq
            w(k) = ssnr(k) / (ssnr(k) + 1.)
        end do
    end function ssnr2optlp

    subroutine lowpass_from_klim( klim, nyq, filter, width )
        integer,        intent(in)    :: klim, nyq
        real,           intent(inout) :: filter(nyq)
        real, optional, intent(in)    :: width
        real    :: freq, lplim_freq, wwidth
        integer :: k
        wwidth = 10.
        if( present(width) ) wwidth = width
        lplim_freq = real(klim)
        do k = 1,nyq
            freq = real(k)
            if( k > klim )then
                filter(k) = 0.
            else if(k .ge. klim - wwidth)then
                filter(k) = (cos(((freq-(lplim_freq-wwidth))/wwidth)*pi)+1.)/2.
            else
                filter(k) = 1.
            endif
        end do
    end subroutine lowpass_from_klim

    subroutine gaussian_filter( freq, smpd, box, filter )
        real,    intent(in)    :: freq, smpd
        integer, intent(in)    :: box
        real,    intent(inout) :: filter(:)
        real    :: A
        integer :: k
        A = freq / smpd / real(box)
        A = 0.5 * (PI * 2.0 * A / 2.35482)**2
        do k = 1,size(filter)
            filter(k) = real(exp(-real(k**2) * A))
        enddo
    end subroutine gaussian_filter

    !> \brief  per-shell factor that caps a B-factor sharpening at shell kcut:
    !!         1 up to kcut and exp((B/4)(s_k^2 - s_kcut^2)) beyond, s_k = k/(box smpd),
    !!         so that a volume sharpened by apply_bfac(B) (exp(-B s^2/4), B < 0) and then
    !!         filtered with it carries the cutoff shell's gain exp(-B s_kcut^2/4) on every
    !!         shell beyond kcut instead of a gain that keeps growing to Nyquist. Identity
    !!         for B >= 0 (nothing to cap). Shell indices as apply_filter reads them
    function bfac_cap_filter( filtsz, box, smpd, bfac, kcut ) result( filt )
        integer, intent(in) :: filtsz, box, kcut
        real,    intent(in) :: smpd, bfac
        real    :: filt(filtsz)
        real    :: s2, s2cut
        integer :: k
        filt = 1.
        if( bfac >= 0. ) return
        s2cut = (real(max(0,kcut)) / (real(box) * smpd))**2
        do k = max(1,kcut+1), filtsz
            s2      = (real(k) / (real(box) * smpd))**2
            filt(k) = exp((bfac / 4.) * (s2 - s2cut))
        end do
    end function bfac_cap_filter

    subroutine mskdiam2lplimits( mskdiam, lpstart,lpstop, lpcen )
        real, intent(in)    :: mskdiam
        real, intent(inout) :: lpstart,lpstop, lpcen
        lpstart = max(min(mskdiam/12., 15.),  8.)
        lpstop  = min(max(mskdiam/22.,  5.),  8.)
        lpcen   = min(max(mskdiam/6.,  20.), 30.)
    end subroutine mskdiam2lplimits

    subroutine mskdiam2streamresthreshold( mskdiam, res_threshold )
        real, intent(in)  :: mskdiam
        real, intent(out) :: res_threshold
        res_threshold = max(min(mskdiam/6.,150.),STREAM_RES_THRESHOLD)
    end subroutine mskdiam2streamresthreshold

    ! Following Grant & Grigorieff; eLife 2015;4:e06980
    subroutine calc_dose_weights( nframes, box, smpd, kV, total_dose, weights )
        integer,           intent(in)    :: nframes, box
        real,              intent(in)    :: smpd, kV, total_dose
        real, allocatable, intent(inout) :: weights(:,:)
        real, parameter :: A=0.245, B=-1.665, C=2.81
        real            :: acc_doses(nframes), spaFreq
        real            :: twoNe, limksq, dose_per_frame
        integer         :: filtsz, iframe, k
        filtsz = fdim(box) - 1
        ! accumulated doses
        dose_per_frame = total_dose / real(nframes) ! e-/Angs2/frame
        do iframe=1,nframes
            acc_doses(iframe) = real(iframe) * dose_per_frame ! e-/Angs2
        end do
        ! voltage scaling
        if( is_equal(kV,200.) )then
            acc_doses = acc_doses / 0.8
        else if( is_equal(kV,100.) )then
            acc_doses = acc_doses / 0.64
        endif
        if( allocated(weights) ) deallocate( weights )
        allocate( weights(nframes,filtsz), source=0.)
        ! dose normalization
        limksq = real(box*smpd)**2.
        do k = 1,filtsz
            spaFreq      = sqrt(real(k*k)/limksq)
            twoNe        = 2.*(A * spaFreq**B + C)
            weights(:,k) = exp(-acc_doses / twoNe)
            weights(:,k) = weights(:,k) / sqrt(sum(weights(:,k) * weights(:,k)))
        enddo
    end subroutine calc_dose_weights

    !>   calculates the resolution values given corrs and res params
    !! \param corrs Fourier shell correlations
    !! \param res resolution value
    subroutine get_resolution( corrs, res, fsc05, fsc0143 )
        real, intent(in)  :: corrs(:), res(:) !<  corrs Fourier shell correlation
        real, intent(out) :: fsc05, fsc0143   !<  fsc05 resolution at FSC=0.5,  fsc0143 resolution at FSC=0.143
        integer           :: n, ires0143, ires05
        n        = size(corrs)
        ! Start the scan from the first shell that is above the 0.143 threshold.
        ! For a normal monotonically-decreasing FSC this is always shell 1.
        ! For a non-monotonic FSC (e.g. spike in low shells) the true crossing
        ! can begin at shell 2 or 3; bounds are guarded by the n >= k checks.
        ires0143 = 1
        if( corrs(1) < 0.143 ) then
            if( n >= 2 .and. corrs(2) > 0.143 ) then
                ires0143 = 2
            else if( n >= 3 .and. corrs(3) > 0.143 ) then
                ires0143 = 3
            end if
        end if
        do while( ires0143 <= n )
            if( corrs(ires0143) >= 0.143 )then
                ires0143 = ires0143 + 1
                cycle
            else
                exit
            endif
        end do
        ires0143 = ires0143 - 1
        if( ires0143 == 0 )then
            fsc0143 = 0.
        else
            fsc0143 = res(ires0143)
        endif
        ! Same non-monotonic FSC logic as for ires0143 above, applied at the 0.5 threshold.
        ires05 = 1
        if( corrs(1) < 0.5 ) then
            if( n >= 2 .and. corrs(2) > 0.5 ) then
                ires05 = 2
            else if( n >= 3 .and. corrs(3) > 0.5 ) then
                ires05 = 3
            end if
        end if
        do while( ires05 <= n )
            if( corrs(ires05) >= 0.5 )then
                ires05 = ires05+1
                cycle
            else
                exit
            endif
        end do
        ires05 = ires05 - 1
        if( ires05 == 0 )then
            fsc05 = 0.
        else
            fsc05 = res(ires05)
        endif
    end subroutine get_resolution

    subroutine get_resolution_at_fsc( corrs, res, fsc_val, res_at_fsc )
        real, intent(in)  :: corrs(:), res(:), fsc_val
        real, intent(out) :: res_at_fsc
        integer           :: n, ires
        n    = size(corrs)
        ires = 1
        do while( ires <= n )
            if( corrs(ires) >= fsc_val )then
                ires = ires + 1
                cycle
            else
                exit
            endif
        end do
        ires = ires - 1
        if( ires == 0 )then
            res_at_fsc = 0.
        else
            res_at_fsc = res(ires)
        endif 
    end subroutine get_resolution_at_fsc

    !> stage @p inf (its lp set) of a box @p box at @p smpd autoscaled towards @p smpd_target: the cropped
    !! box, sampling and scale, the translation limit and whether it crops
    subroutine stage_autoscale( box, smpd, smpd_target, inf )
        use simple_magic_boxes, only: autoscale
        integer,           intent(in)    :: box
        real,              intent(in)    :: smpd, smpd_target
        type(lp_crop_inf), intent(inout) :: inf
        call autoscale(box, smpd, smpd_target, inf%box_crop, inf%smpd_crop, inf%scale, minbox=STAGE_MINBOX)
        inf%trslim      = stage_trslim(inf%smpd_crop)
        inf%l_autoscale = inf%box_crop < box
    end subroutine stage_autoscale

    !> stage @p inf cropped to the given box @p box_crop of a box @p box at @p smpd
    subroutine stage_box( box, smpd, box_crop, inf )
        integer,           intent(in)    :: box, box_crop
        real,              intent(in)    :: smpd
        type(lp_crop_inf), intent(inout) :: inf
        inf%box_crop    = box_crop
        inf%scale       = real(box_crop) / real(box)
        inf%smpd_crop   = smpd / inf%scale
        inf%trslim      = stage_trslim(inf%smpd_crop)
        inf%l_autoscale = box_crop < box
    end subroutine stage_box

    !> a stage's translation limit (pixels): an alpha-helix width at the cropped sampling, within [2,8]
    elemental real function stage_trslim( smpd_crop )
        real, intent(in) :: smpd_crop
        stage_trslim = min(8., max(2.0, AHELIX_WIDTH / smpd_crop))
    end function stage_trslim

    subroutine lpstages_fast( box, nstages, smpd, lpstart, lpstop, lpinfo, force_lpstart, verbose, trslim_override )
        integer,           intent(in)  :: box, nstages
        real,              intent(in)  :: smpd, lpstart, lpstop
        type(lp_crop_inf), intent(out) :: lpinfo(nstages)
        logical, optional, intent(in)  :: force_lpstart, verbose
        real,    optional, intent(in)  :: trslim_override
        real,    parameter :: LP2SMPD_TARGET   = 0.4
        real,    parameter :: SMPD_TARGET_MIN  = 2.5
        integer :: i
        real    :: lpstart_eff
        logical :: l_verbose
        l_verbose = .true.
        if( present(verbose) ) l_verbose = verbose
        if( nstages < 1 ) THROW_HARD('nstages must be >= 1 in lpstages_fast')
        lpstart_eff = max(lpstart, LPSTAGES_STAGE1_LP_FLOOR)
        if( present(force_lpstart) )then
            if( force_lpstart ) lpstart_eff = lpstart
        endif
        lpinfo(:)%l_lpset = .true.
        do i = 1,nstages
            if( nstages == 1 )then
                lpinfo(i)%lp = max(lpstop, lpstart_eff)
            else
                lpinfo(i)%lp = lpstart_eff - real(i-1)*(lpstart_eff-lpstop)/real(nstages-1)
            endif
            call stage_autoscale(box, smpd, max(SMPD_TARGET_MIN, lpinfo(i)%lp * LP2SMPD_TARGET), lpinfo(i))
        enddo
        if( present(trslim_override) )then
            if( trslim_override > 0. ) lpinfo(:)%trslim = trslim_override
        endif
        if( l_verbose )then
            do i = 1,nstages
                print *, 'lpset lp box_crop smpd_crop trslim ', lpinfo(i)%l_lpset, lpinfo(i)%lp, lpinfo(i)%box_crop, &
                    &lpinfo(i)%smpd_crop, lpinfo(i)%trslim
            enddo
        endif
    end subroutine lpstages_fast

    !> a linear ladder from lpstart to lpstop; the end stages cropped to their limits, a linear box ramp between
    subroutine lpstages_setlims( box, nstages, smpd, lpstart, lpstop, lpinfo, verbose, first_stage, trslim_override )
        use simple_magic_boxes
        integer,           intent(in)  :: box, nstages
        real,              intent(in)  :: smpd, lpstart, lpstop
        type(lp_crop_inf), intent(out) :: lpinfo(nstages)
        logical, optional, intent(in)  :: verbose
        integer, optional, intent(in)  :: first_stage
        real,    optional, intent(in)  :: trslim_override
        real,    parameter :: LP2SMPD_TARGET   = 1./3.
        real,    parameter :: SMPD_TARGET_MIN  = 2.0
        integer :: i, box_trial, ifirst
        real    :: rbox_stepsz
        logical :: l_verbose
        l_verbose = .true.
        if( present(verbose) ) l_verbose = verbose
        if( nstages < 1 ) THROW_HARD('nstages must be >= 1 in lpstages_setlims')
        ifirst = 1
        if( present(first_stage) )then
            if( first_stage < 1 .or. first_stage > nstages ) THROW_HARD('first_stage must be within the lpstages_setlims schedule')
            ifirst = first_stage
        endif
        lpinfo(:)%l_lpset = .true.
        do i = 1,nstages
            if( present(first_stage) )then
                if( i <= ifirst )then
                    lpinfo(i)%lp = lpstart
                else
                    lpinfo(i)%lp = lpstart - real(i-ifirst)*(lpstart-lpstop)/real(nstages-ifirst)
                endif
            else if( nstages == 1 )then
                lpinfo(i)%lp = lpstop
            else
                lpinfo(i)%lp = lpstart - real(i-1)*(lpstart-lpstop)/real(nstages-1)
            endif
            call stage_autoscale(box, smpd, max(SMPD_TARGET_MIN, lpinfo(i)%lp * LP2SMPD_TARGET), lpinfo(i))
        enddo
        if( nstages > ifirst )then
            rbox_stepsz = real(lpinfo(nstages)%box_crop - lpinfo(ifirst)%box_crop)/real(nstages - ifirst)
            do i = ifirst + 1, nstages - 1 ! linear box_crop scheme
                box_trial = nint(real(lpinfo(ifirst)%box_crop) + real(i-ifirst)*rbox_stepsz)
                call stage_box(box, smpd, min(lpinfo(nstages)%box_crop, find_magic_box(box_trial)), lpinfo(i))
            end do
        endif
        if( present(trslim_override) )then
            if( trslim_override > 0. ) lpinfo(:)%trslim = trslim_override
        endif
        if( l_verbose )then
            do i = ifirst, nstages
                write(logfhandle,*) 'lpset lp box_crop smpd_crop trslim ', lpinfo(i)%l_lpset, lpinfo(i)%lp, lpinfo(i)%box_crop, &
                    &lpinfo(i)%smpd_crop, lpinfo(i)%trslim
            end do
        endif
    end subroutine lpstages_setlims

end module simple_estimate_ssnr
