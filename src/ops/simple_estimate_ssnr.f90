! spectral signal-to-noise ratio estimation routines
module simple_estimate_ssnr
use simple_defs
use simple_error, only: allocchk, simple_exception
use simple_math,  only: find
implicit none

public :: fsc2ssnr, fsc2optlp, fsc2optlp_sub, ssnr2fsc, ssnr2optlp, acc_dose2filter, dose_weight
private
#include "simple_local_flags.inc"

contains

    !> \brief  converts the FSC to SSNR (the 2.* is because of the division of the data)
    function fsc2ssnr( corrs ) result( ssnr )
        real, intent(in)  :: corrs(:) !<  instrument FSC
        real, allocatable :: ssnr(:) !<  instrument SSNR
        integer :: nyq, k
        real    :: fsc
        nyq = size(corrs)
        allocate( ssnr(nyq),stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("in simple_estimate_ssnr::fsc2ssnr ssnr ",alloc_stat)
        do k=1,nyq
            fsc = min(abs(corrs(k)),0.999)
            ssnr(k) = (2.*fsc)/(1.-fsc)
        end do
    end function fsc2ssnr

    !> \brief  converts the FSC to the optimal low-pass filter
    function fsc2optlp( corrs ) result( filt )
        real, intent(in)  :: corrs(:) !< fsc plot (correlations)
        real, allocatable :: filt(:)  !< output filter coefficients
        integer :: nyq
        nyq = size(corrs)
        allocate( filt(nyq),stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("in simple_estimate_ssnr::fsc2optlp filt ",alloc_stat)
        filt = 0.
        where( corrs > 0. )     filt = sqrt( 2. * corrs / (corrs + 1.) )
        where( filt  > 0.9999 ) filt = 0.99999
    end function fsc2optlp

    !> \brief  converts the FSC to the optimal low-pass filter
    subroutine fsc2optlp_sub( filtsz, corrs, filt )
        integer, intent(in)  :: filtsz        !< sz of filter
        real,    intent(in)  :: corrs(filtsz) !< fsc plot (correlations)
        real,    intent(out) :: filt(filtsz)  !< output filter coefficients
        filt = 0.
        where( corrs > 0. )     filt = sqrt( 2. * corrs / (corrs + 1.) )
        where( filt  > 0.9999 ) filt = 0.99999
    end subroutine fsc2optlp_sub

    !> \brief  converts the SSNR to FSC
    function ssnr2fsc( ssnr ) result( corrs )
        real, intent(in)  :: ssnr(:)  !< input SSNR array
        real, allocatable :: corrs(:) !< output FSC result
        integer :: nyq, k
        nyq = size(ssnr)
        allocate( corrs(nyq),stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("in simple_estimate_ssnr::ssnr2fsc corrs ",alloc_stat)
        do k=1,nyq
            corrs(k) = ssnr(k)/(ssnr(k)+1.)
        end do
    end function ssnr2fsc

    ! !> \brief  converts the SSNR 2 the optimal low-pass filter
    function ssnr2optlp( ssnr ) result( w )
        real, intent(in)  :: ssnr(:) !<  instrument SSNR
        real, allocatable :: w(:) !<  FIR low-pass filter
        integer :: nyq, k
        nyq = size(ssnr)
        allocate( w(nyq),stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("in simple_estimate_ssnr::ssnr2optlp w ",alloc_stat)
        do k=1,nyq
            w(k) = ssnr(k)/(ssnr(k)+1.)
        end do
    end function ssnr2optlp

    !> DOSE FILTERING (Grant, Grigorieff eLife 2015)
    !! input is template image, accumulative dose (in e/A2) and acceleration voltage
    !!         output is filter coefficients
    !! \f$  \mathrm{dose}_\mathrm{acc} = \int^{N}_{1} \mathrm{dose\_weight}(a,F,V),\ n_\mathrm{e}/\si{\angstrom\squared}  \f$
    !! \param acc_dose accumulative dose (in \f$n_\mathrm{e}^- per \si{\angstrom\squared}\f$)
    function acc_dose2filter( img, acc_dose, kV ) result( filter )
        use simple_image,   only: image
        type(image), intent(in) :: img           !< input image
        real,        intent(in) :: acc_dose, kV  !< acceleration voltage
        real, allocatable       :: filter(:)
        integer :: find, sz
        sz = img%get_filtsz()
        allocate(filter(sz),stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("simple_estimate_ssnr::acc_dose2filter ",alloc_stat)
        do find=1,sz
            filter(find) = dose_weight(acc_dose, img%get_spat_freq(find), kV)
        end do
    end function acc_dose2filter

    !>  \brief Calculate dose weight. Input is accumulative dose (in e/A2) and
    !>  spatial frequency (in 1/A)
    !!         output is resolution dependent weight applied to individual frames
    !!         before correlation search and averaging
    !!
    !! \f$  \mathrm{dose\_weight}(a,F,V) = \exp\left(- \frac{A_\mathrm{dose}}{k\times (2.0\times A\times f^B + C)} \right), \f$
    !! where \f$k\f$ is 0.75 for \f$V<200\f$ kV, 1.0 for \f$200 \leqslant  V \leqslant  300\f$
    real function dose_weight( acc_dose, spat_freq, kV )
        real, intent(in) :: acc_dose                !< accumulative dose (in e/A2)
        real, intent(in) :: spat_freq               !< spatial frequency (in 1/A)
        real, intent(in) :: kV                      !< accelleration voltage
        real, parameter  :: A=0.245, B=-1.665, C=2.81, kV_factor=0.75
        real             :: critical_exp !< critical exposure (only depends on spatial frequency)
        critical_exp = A*(spat_freq**B)+C
        if( abs(kV-300.) < 0.001 )then
            ! critical exposure does not need modification
        else if( abs(kV-200.) < 0.001 )then
            ! critical exposure at 200 kV expected to be ~25% lower
            critical_exp = critical_exp*kV_factor
        else
            THROW_HARD('unsupported kV (acceleration voltage); dose_weight')
        endif
        dose_weight = exp(-acc_dose/(2.0*critical_exp))
    end function dose_weight

    subroutine local_res( even, odd, mskimg, corr_thres, locres_finds, half_winsz )
        use simple_image, only: image
        class(image),         intent(inout) :: even, odd, mskimg
        real,                 intent(in)    :: corr_thres
        integer, allocatable, intent(out)   :: locres_finds(:,:,:)
        integer, optional,    intent(in)    :: half_winsz
        logical, allocatable :: l_mask(:,:,:)
        real,    allocatable :: rmat_even(:,:,:), rmat_odd(:,:,:)
        integer     :: hwinsz, filtsz, ldim(3), i, j, k, kind
        logical     :: is2d
        type(image) :: ecopy, ocopy
        real        :: smpd, cc
        ! set parameters
        hwinsz = 3
        if( present(half_winsz) ) hwinsz = half_winsz
        ldim   = even%get_ldim()
        is2d   = ldim(3) == 1
        filtsz = even%get_filtsz()
        smpd   = even%get_smpd()
        l_mask = mskimg%bin2logical()
        ! to avoid allocation in the loop
        call ecopy%new(ldim, smpd)
        call ocopy%new(ldim, smpd)
        if( allocated(locres_finds) ) deallocate(locres_finds)
        allocate(locres_finds(ldim(1),ldim(2),ldim(3)), source=0)
        allocate(rmat_even(ldim(1),ldim(2),ldim(3)), rmat_odd(ldim(1),ldim(2),ldim(3)))
        call even%fft
        call odd%fft
        ! loop over resolution shells
        do kind=1,filtsz
            call ecopy%copy(even)
            call ocopy%copy(odd)
            call ecopy%tophat(k)
            call ocopy%tophat(k)
            call ecopy%ifft
            call ocopy%ifft
            call ecopy%get_rmat_sub(rmat_even)
            call ocopy%get_rmat_sub(rmat_odd)
            do k=1,ldim(3)
                do j=1,ldim(2)
                    do i=1,ldim(1)
                        if( l_mask(i,j,k) )then
                            cc = neigh_cc([i,j,k])
                            if( cc >= corr_thres ) locres_finds(i,j,k) = kind
                        endif
                    end do
                end do
            end do
        end do
        ! make sure no 0 elements
        where( l_mask )
            where(locres_finds == 0) locres_finds = filtsz
        end where
        ! return e/o images in real-space
        call even%ifft
        call odd%ifft
        ! destruct
        deallocate(l_mask, rmat_even, rmat_odd)
        call ecopy%kill
        call ocopy%kill

        contains

            function neigh_cc( loc ) result( cc )
                integer, intent(in) :: loc(3)
                integer :: lb(3), ub(3), npix
                real    :: ae, ao, diffe, diffo, see, soo, seo, cc
                ! set bounds
                lb = loc - hwinsz
                where(lb < 1) lb = 1
                ub = loc + hwinsz
                where(ub > ldim ) ub = ldim
                if( is2d )then
                    lb(3) = 1
                    ub(3) = 1
                endif
                ! calc avgs
                npix = product(ub - lb + 1)
                ae   = sum(rmat_even(lb(1):ub(1),lb(2):ub(2),lb(2):ub(3))) / real(npix)
                ao   = sum(rmat_odd(lb(1):ub(1),lb(2):ub(2),lb(2):ub(3))) / real(npix)
                see  = 0.
                soo  = 0.
                seo  = 0.
                ! calc corr
                do k=lb(3),ub(3)
                    do j=lb(2),ub(2)
                        do i=lb(1),ub(1)
                            diffe = rmat_even(i,j,k) - ae
                            diffo = rmat_odd(i,j,k)  - ao
                            see   = see + diffe * diffe
                            soo   = soo + diffo * diffo
                            seo   = seo + diffe * diffo
                        end do
                    end do
                end do
                if( see > 0. .and. soo > 0. )then
                    cc = seo / sqrt(see * soo)
                else
                    cc = 0.
                endif
            end function neigh_cc

    end subroutine local_res

end module simple_estimate_ssnr
