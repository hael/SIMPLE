! spectral signal-to-noise ratio estimation routines
module simple_estimate_ssnr
include 'simple_lib.f08'
implicit none

public :: fsc2ssnr, fsc2optlp, fsc2optlp_sub, ssnr2fsc, ssnr2optlp, acc_dose2filter, dose_weight
public :: local_res
private
#include "simple_local_flags.inc"

contains

    !> \brief  converts the FSC to SSNR (the 2.* is because of the division of the data)
    function fsc2ssnr( corrs ) result( ssnr )
        real, intent(in)  :: corrs(:) !< instrument FSC
        real, allocatable :: ssnr(:)  !< instrument SSNR
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

    ! even odd images assumed to be appropriately masked before calling this routine
    ! mskimg should be a hard mask (real-space)
    subroutine local_res( even, odd, mskimg, corr_thres, locres_finds, half_winsz )
        use simple_image, only: image
        class(image),         intent(inout) :: even, odd, mskimg
        real,                 intent(in)    :: corr_thres
        integer, allocatable, intent(out)   :: locres_finds(:,:,:)
        integer, optional,    intent(in)    :: half_winsz
        logical,  allocatable :: l_mask(:,:,:)
        real,     allocatable :: fsc(:), res(:)
        real(dp), allocatable :: vec1(:), vec2(:)
        integer     :: hwinsz, filtsz, ldim(3), i, j, k, kind, vecsz, find_hres, find_lres
        logical     :: is2d
        type(image) :: ecopy, ocopy
        real        :: smpd, cc, ccavg, res_fsc05, res_fsc0143
        ! check input
        if( .not. even%is_ft() ) THROW_HARD('even vol not FTed; local_res')
        if( .not. odd%is_ft()  ) THROW_HARD('odd vol not FTed; local_res')
        ! set parameters
        hwinsz = 3
        if( present(half_winsz) ) hwinsz = half_winsz
        ldim   = even%get_ldim()
        is2d   = ldim(3) == 1
        if( is2d )then
            vecsz = (2*hwinsz + 1)**2
        else
            vecsz = (2*hwinsz + 1)**3
        endif
        filtsz = even%get_filtsz()
        smpd   = even%get_smpd()
        l_mask = mskimg%bin2logical()
        res    = even%get_res()
        ! to avoid allocation in the loop
        call ecopy%new(ldim, smpd)
        call ocopy%new(ldim, smpd)
        if( allocated(locres_finds) ) deallocate(locres_finds)
        allocate(locres_finds(ldim(1),ldim(2),ldim(3)), vec1(vecsz), vec2(vecsz), fsc(filtsz))
        locres_finds = 0
        ! loop over resolution shells
        do kind=1,filtsz
            ! call progress(kind, filtsz)
            call ecopy%copy(even)
            call ocopy%copy(odd)
            call ecopy%tophat(kind)
            call ocopy%tophat(kind)
            call ecopy%ifft
            call ocopy%ifft
            ccavg = 0.
            ! parallel section needs to start here as the above steps are threaded as well
            !$omp parallel do collapse(3) default(shared) private(i,j,k,cc)&
            !$omp schedule(static) proc_bind(close) reduction(+:ccavg)
            do k=1,ldim(3)
                do j=1,ldim(2)
                    do i=1,ldim(1)
                        if( .not. l_mask(i,j,k) ) cycle
                        cc = neigh_cc([i,j,k])
                        ccavg = ccavg + cc
                        if( cc >= corr_thres ) locres_finds(i,j,k) = kind
                    end do
                end do
            end do
            !$omp end parallel do
            ccavg     = ccavg / real(count(l_mask))
            fsc(kind) = ccavg
            write(*,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(kind), '>>> CORRELATION:', fsc(kind)
        end do
        call get_resolution(fsc, res, res_fsc05, res_fsc0143)
        write(*,'(A,1X,F6.2)') '>>> GLOBAL RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
        write(*,'(A,1X,F6.2)') '>>> GLOBAL RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
        find_hres = maxval(locres_finds, l_mask)
        find_lres = minval(locres_finds, l_mask .and. locres_finds > 0)
        write(*,'(A,1X,F6.2)') '>>> HIGHEST LOCAL RESOLUTION DETERMINED TO:', even%get_lp(find_hres)
        write(*,'(A,1X,F6.2)') '>>> LOWEST  LOCAL RESOLUTION DETERMINED TO:', even%get_lp(find_lres)
        ! make sure no 0 elements
        where( l_mask )
            where(locres_finds == 0) locres_finds = filtsz
        end where
        ! destruct
        deallocate(l_mask)
        call ecopy%kill
        call ocopy%kill

        contains

            function neigh_cc( loc ) result( cc )
                integer, intent(in) :: loc(3)
                integer  :: lb(3), ub(3), i, j, k, cnt
                real     :: cc
                ! set bounds
                lb = loc - hwinsz
                ub = loc + hwinsz
                if( is2d )then
                    lb(3) = 1
                    ub(3) = 1
                endif
                ! extract components
                cnt = 0
                do k=lb(3),ub(3)
                    if( k < 1 .or. k > ldim(3) ) cycle
                    do j=lb(2),ub(2)
                        if( j < 1 .or. j > ldim(2) ) cycle
                        do i=lb(1),ub(1)
                            if( i < 1 .or. i > ldim(1) ) cycle
                            cnt = cnt + 1
                            vec1(cnt) = dble(ecopy%get_rmat_at(i,j,k))
                            vec2(cnt) = dble(ocopy%get_rmat_at(i,j,k))
                        enddo
                    enddo
                enddo
                cc = pearsn_serial_8(vec1(:cnt), vec2(:cnt))
            end function neigh_cc

    end subroutine local_res

end module simple_estimate_ssnr
