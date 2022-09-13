! spectral signal-to-noise ratio estimation routines
module simple_estimate_ssnr
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
implicit none

public :: fsc2ssnr, fsc2TVfilt, fsc2TVfilt_fast, fsc2optlp, fsc2optlp_sub, ssnr2fsc, ssnr2optlp, subsample_optlp
public :: plot_fsc, lowpass_from_klim, mskdiam2lplimits, calc_dose_weights
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

    !> \brief  converts the FSC to TV-Filter SSNR/(SSNR + 1)
    subroutine fsc2TVfilt( fsc, flims, TVfilt )
        real,    intent(in)    :: fsc(:)            !< FSC correlations (depnds on k <=> resolution)
        integer, intent(in)    :: flims(3,2)
        real,    intent(inout) :: TVfilt(size(fsc)) !< TV-Filter derived from FSC
        integer :: nfcomps(size(fsc)), nyq, h, k, l, sh
        real :: nr, snr,  one_o_sqrt_nr
        nyq     = size(fsc)
        nfcomps = 0
        do h = flims(1,1),flims(1,2)
            do k = flims(2,1),flims(2,2)
                do l = flims(3,1),flims(3,2)
                    sh = nint(hyp(real(h),real(k),real(l)))
                    if (sh < 1 .or. sh > nyq) cycle
                    nfcomps(sh) = nfcomps(sh) + 1
                end do
            end do
        end do
        do k = 1,nyq
            nr = real(nfcomps(k))
            one_o_sqrt_nr = 1./sqrt(nr)
            if( fsc(k) > 1. )then ! round-off error
                TVfilt(k) = 1.
            else if( fsc(k) < one_o_sqrt_nr )then
                TVfilt(k) = 0.
            else
                snr = sqrt(1./nr - (one_o_sqrt_nr - fsc(k))/(1. - fsc(k))) - one_o_sqrt_nr
                TVfilt(k) = snr / (snr + 1.)
            endif
        enddo
    end subroutine fsc2TVfilt

    subroutine fsc2TVfilt_fast( fsc, nfcomps, TVfilt )
        real,    intent(in)    :: fsc(:)            !< FSC correlations (depnds on k <=> resolution)
        integer, intent(in)    :: nfcomps(size(fsc))
        real,    intent(inout) :: TVfilt(size(fsc)) !< TV-Filter derived from FSC
        integer :: nyq, h, k, l, sh
        real :: nr, snr, one_o_sqrt_nr
        nyq = size(fsc)
        do k = 1,nyq
            nr = real(nfcomps(k))
            one_o_sqrt_nr = 1./sqrt(nr)
            if( fsc(k) > 1. )then ! round-off error
                TVfilt(k) = 1.
            else if( fsc(k) < one_o_sqrt_nr )then
                TVfilt(k) = 0.
            else
                snr = sqrt(1./nr - (one_o_sqrt_nr - fsc(k))/(1. - fsc(k))) - one_o_sqrt_nr
                TVfilt(k) = snr / (snr + 1.)
            endif
        enddo
    end subroutine fsc2TVfilt_fast

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

    subroutine subsample_optlp( filtsz, subfiltsz, filt, subfilt )
        integer, intent(in)  :: filtsz, subfiltsz   !< sz of filters
        real,    intent(in)  :: filt(filtsz)        !< filter coefficients
        real,    intent(out) :: subfilt(subfiltsz)  !< output filter coefficients
        real    :: x, fracx, step
        integer :: i, floorx
        if( filtsz < subfiltsz )then
            THROW_HARD('Invalid filter sizes!')
        else if( filtsz == subfiltsz )then
            subfilt = filt
        else
            x    = 1.
            step = real(filtsz-1) / real(subfiltsz-1)
            do i = 2,subfiltsz-1
                x          = x+step
                floorx     = floor(x)
                fracx      = x-real(floorx)
                subfilt(i) = (1.-fracx)*filt(floorx) + fracx*filt(ceiling(x))
                subfilt(i) = max(min(subfilt(i),1.),0.) ! always in [0.,1.]
            enddo
            subfilt(1)         = max(min(filt(1),1.),0.)
            subfilt(subfiltsz) = max(min(filt(filtsz),1.),0.)
        endif
    end subroutine subsample_optlp

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
    subroutine calc_dose_weights( nframes, box, smpd, kV, exp_time, dose_rate, weights )
        integer,           intent(in)    :: nframes, box
        real,              intent(in)    :: smpd, kV, exp_time, dose_rate
        real, allocatable, intent(inout) :: weights(:,:)
        real, parameter :: A=0.245, B=-1.665, C=2.81
        real            :: frame_dose(nframes), acc_doses(nframes), spaFreq, current_time
        real            :: twoNe, limksq, time_per_frame
        integer         :: filtsz, iframe, k
        time_per_frame = exp_time/real(nframes)               ! unit: s
        do iframe=1,nframes
            current_time      = real(iframe) * time_per_frame ! unit: s
            acc_doses(iframe) = dose_rate * current_time      ! unit: e/A2/s * s = e/A2
        end do
        filtsz = fdim(box) - 1
        ! doses
        limksq = real(box*smpd)**2.
        do iframe = 1,nframes
            frame_dose(iframe) = acc_doses(iframe)
            if( is_equal(kV,200.) )then
                frame_dose(iframe) = frame_dose(iframe) / 0.8
            else if( is_equal(kV,100.) )then
                frame_dose(iframe) = frame_dose(iframe) / 0.64
            endif
        enddo
        if( allocated(weights) ) deallocate( weights )
        allocate( weights(nframes,filtsz), source=0.)
        ! dose normalization
        do k = 1,filtsz
            spaFreq      = sqrt(real(k*k)/limksq)
            twoNe        = 2.*(A * spaFreq**B + C)
            weights(:,k) = exp(-frame_dose / twoNe)
            weights(:,k) = weights(:,k) / sqrt(sum(weights(:,k) * weights(:,k)))
        enddo
    end subroutine calc_dose_weights

    subroutine plot_fsc( n, fsc, res, smpd, tmpl_fname )
        use CPlot2D_wrapper_module
        integer,           intent(in) :: n
        real,              intent(in) :: fsc(n), res(n), smpd
        character(len=*),  intent(in) :: tmpl_fname
        real, parameter           :: SCALE = 40.
        type(str4arr)             :: title
        type(CPlot2D_type)        :: plot2D
        type(CDataSet_type)       :: dataSet
        character(len=LONGSTRLEN) :: ps2pdf_cmd, fname_pdf, fname_eps
        integer  :: k,iostat
        if( n == 0 ) THROW_HARD('Empty FSC vector; plot_fsc')
        fname_eps  = trim(tmpl_fname)//'.eps'
        fname_pdf  = trim(tmpl_fname)//'.pdf'
        call CPlot2D__new(plot2D, trim(tmpl_fname)//C_NULL_CHAR)
        call CPlot2D__SetXAxisSize(plot2D, 400.d0)
        call CPlot2D__SetYAxisSize(plot2D, 400.d0)
        call CPlot2D__SetDrawLegend(plot2D, C_FALSE)
        call CPlot2D__SetFlipY(plot2D, C_FALSE)
        call CDataSet__new(dataSet)
        call CDataSet__SetDrawMarker(dataSet, C_FALSE)
        call CDataSet__SetDatasetColor(dataSet, 0.d0,0.d0,1.d0)
        do k = 1,n
            call CDataSet_addpoint(dataSet, 1.0/res(k), fsc(k))
        end do
        call CPlot2D__AddDataSet(plot2D, dataset)
        call CDataSet__delete(dataset)
        title%str = 'Resolution (Angstroms^-1)'//C_NULL_CHAR
        call CPlot2D__SetXAxisTitle(plot2D, title%str)
        title%str = 'Fourier Shell Correlations'//C_NULL_CHAR
        call CPlot2D__SetYAxisTitle(plot2D, title%str)
        call CPlot2D__OutputPostScriptPlot(plot2D, trim(fname_eps)//C_NULL_CHAR)
        call CPlot2D__delete(plot2D)
        ! conversion to PDF
        ps2pdf_cmd = 'gs -q -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dDEVICEWIDTHPOINTS=600 -dDEVICEHEIGHTPOINTS=600 -sOutputFile='&
            &//trim(fname_pdf)//' '//trim(fname_eps)
        call exec_cmdline(trim(adjustl(ps2pdf_cmd)), suppress_errors=.true., exitstat=iostat)
        if( iostat == 0 ) call del_file(fname_eps)
    end subroutine plot_fsc

end module simple_estimate_ssnr
