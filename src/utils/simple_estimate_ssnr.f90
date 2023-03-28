! spectral signal-to-noise ratio estimation routines
module simple_estimate_ssnr
!$ use omp_lib
!$ use omp_lib_kinds
use simple_defs
use simple_srch_sort_loc
use simple_syslib
use simple_math_ft
use simple_strings
use simple_fileio
implicit none

public :: fsc2ssnr, fsc2optlp, fsc2optlp_sub, ssnr2fsc, ssnr2optlp
public :: lowpass_from_klim, mskdiam2lplimits, calc_dose_weights, get_resolution
private
#include "simple_local_flags.inc"

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

    !> \brief  converts the FSC to the optimal low-pass filter
    subroutine fsc2optlp_sub( filtsz, corrs, filt )
        integer, intent(in)  :: filtsz        !< sz of filter
        real,    intent(in)  :: corrs(filtsz) !< fsc plot (correlations)
        real,    intent(out) :: filt(filtsz)  !< output filter coefficients
        filt = 0.
        where( corrs > 0. )     filt = 2. * corrs / (corrs + 1.)
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

    subroutine mskdiam2lplimits( mskdiam, lpstart,lpstop, lpcen )
        real, intent(in)    :: mskdiam
        real, intent(inout) :: lpstart,lpstop, lpcen
        lpstart = max(min(mskdiam/12., 15.),  8.)
        lpstop  = min(max(mskdiam/22.,  5.),  8.)
        lpcen   = min(max(mskdiam/6.,  20.), 30.)
    end subroutine mskdiam2lplimits


    ! Following Grant & Grigorieff; eLife 2015;4:e06980
    subroutine calc_dose_weights( nframes, box, smpd, kV, total_dose, weights )
        integer,           intent(in)    :: nframes, box
        real,              intent(in)    :: smpd, kV, total_dose
        real, allocatable, intent(inout) :: weights(:,:)
        real, parameter :: A=0.245, B=-1.665, C=2.81
        real            :: acc_doses(nframes), spaFreq, current_time
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
        n = size(corrs)
        ires0143 = 1
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
        ires05 = 1
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

end module simple_estimate_ssnr
