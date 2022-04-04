! spectral signal-to-noise ratio estimation routines
module simple_estimate_ssnr
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
implicit none

public :: fsc2ssnr, fsc2TVfilt, fsc2TVfilt_fast, fsc2optlp, fsc2optlp_sub, ssnr2fsc, ssnr2optlp, subsample_optlp
public :: nonuniform_phase_ran, nonuniform_fscTVfilt, local_res_lp
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
        real, intent(in)  :: mskdiam
        real, intent(out) :: lpstart,lpstop, lpcen
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

    subroutine nonuniform_phase_ran( even, odd, mskimg )
        use simple_image, only: image
        class(image), intent(inout) :: even, odd, mskimg
        integer,      parameter     :: SUBBOX=32, LSHIFT=15, RSHIFT=16, CPIX=LSHIFT + 1, CHUNKSZ=20
        real,         parameter     :: SUBMSK=real(SUBBOX)/2. - COSMSKHALFWIDTH - 1.
        type(image),  allocatable   :: subvols_even(:) ! one per thread
        type(image),  allocatable   :: subvols_odd(:)  ! one per thread
        logical,      allocatable   :: l_mask(:,:,:)
        integer     :: ldim(3), i, j, k, ithr, cnt, npix, lb(3), ub(3)
        real        :: smpd
        ! check input
        if( even%is_ft()   ) THROW_HARD('even vol FTed; nonuniform_phase_ran')
        if( odd%is_ft()    ) THROW_HARD('odd  vol FTed; nonuniform_phase_ran')
        if( mskimg%is_ft() ) THROW_HARD('msk  vol FTed; nonuniform_phase_ran')
        ! set parameters
        ldim    = even%get_ldim()
        if( ldim(3) == 1 ) THROW_HARD('not intended for 2D images; nonuniform_phase_ran')
        smpd    = even%get_smpd()
        l_mask  = mskimg%bin2logical()
        npix    = count(l_mask)
        ! construct
        allocate(subvols_even(nthr_glob), subvols_odd(nthr_glob))
        do ithr = 1, nthr_glob
            call subvols_even(ithr)%new([SUBBOX,SUBBOX,SUBBOX], smpd, wthreads=.false.)
            call subvols_odd(ithr)%new( [SUBBOX,SUBBOX,SUBBOX], smpd, wthreads=.false.)
        end do
        ! determine loop bounds for better load balancing in the following parallel loop
        call bounds_from_mask3D( l_mask, lb, ub )
        ! loop over pixels
        !$omp parallel do collapse(3) default(shared) private(i,j,k,ithr) schedule(dynamic,CHUNKSZ) proc_bind(close)
        do k = lb(3), ub(3)
            do j = lb(2), ub(2)
                do i = lb(1), ub(1)
                    if( .not.l_mask(i,j,k) ) cycle
                    ithr = omp_get_thread_num() + 1
                    call set_subvols_msk_fft_phran_ifft([i,j,k], ithr)
                    call even%set_rmat_at(i,j,k, subvols_even(ithr)%get_rmat_at(CPIX,CPIX,CPIX))
                    call odd%set_rmat_at(i,j,k, subvols_odd(ithr)%get_rmat_at(CPIX,CPIX,CPIX))
                end do
            end do
        end do
        !$omp end parallel do
        ! destruct
        do ithr = 1, nthr_glob
            call subvols_even(ithr)%kill
            call subvols_odd(ithr)%kill
        end do
        deallocate(l_mask)

        contains

            subroutine set_subvols_msk_fft_phran_ifft( loc, ithr )
                integer, intent(in) :: loc(3), ithr
                integer  :: lb(3), ub(3), i, j, k, ii, jj, kk, isub, jsub, ksub
                call subvols_even(ithr)%zero_and_unflag_ft
                call subvols_odd(ithr)%zero_and_unflag_ft
                ! set bounds
                lb = loc - LSHIFT
                ub = loc + RSHIFT
                ! extract components
                kk  = -LSHIFT
                do k = lb(3), ub(3)
                    if( k < 1 .and. k > ldim(3) ) cycle
                    jj = -LSHIFT
                    do j = lb(2), ub(2)
                        if( j < 1 .and. j > ldim(2) ) cycle
                        ii = -LSHIFT
                        do i = lb(1), ub(1)
                            if( i < 1 .and. i > ldim(1) ) cycle
                            isub = ii + LSHIFT + 1
                            jsub = jj + LSHIFT + 1
                            ksub = kk + LSHIFT + 1
                            call subvols_even(ithr)%set_rmat_at(isub,jsub,ksub, even%get_rmat_at(i,j,k))
                            call subvols_odd(ithr)%set_rmat_at( isub,jsub,ksub, odd%get_rmat_at(i,j,k))
                            ii = ii + 1
                        enddo
                        jj = jj + 1
                    enddo
                    kk = kk + 1
                enddo
                ! mask
                call subvols_even(ithr)%mask(SUBMSK, 'soft')
                call subvols_odd(ithr)%mask( SUBMSK, 'soft')
                ! fft
                call subvols_even(ithr)%fft
                call subvols_odd(ithr)%fft
                ! randomize phases below noise power
                call subvols_even(ithr)%ran_phases_below_noise_power(subvols_odd(ithr))
                ! back to real space
                call subvols_even(ithr)%ifft
                call subvols_odd(ithr)%ifft
            end subroutine set_subvols_msk_fft_phran_ifft

    end subroutine nonuniform_phase_ran

    subroutine nonuniform_fscTVfilt( even, odd, mskimg, l_match_filt, map2filt, phran )
        use simple_image, only: image
        class(image),           intent(inout) :: even, odd, mskimg
        logical,                intent(in)    :: l_match_filt
        class(image), optional, intent(inout) :: map2filt
        logical,      optional, intent(in)    :: phran
        integer,      parameter   :: SUBBOX=32, LSHIFT=15, RSHIFT=16, CPIX=LSHIFT + 1, CHUNKSZ=20
        real,         parameter   :: SUBMSK=real(SUBBOX)/2. - COSMSKHALFWIDTH - 1.
        type(image),  allocatable :: subvols_even(:)  ! one per thread
        type(image),  allocatable :: subvols_odd(:)   ! one per thread
        type(image),  allocatable :: subvols_2filt(:) ! one per thread
        type(image)               :: tmpvol, tmpvol_odd
        logical,      allocatable :: l_mask(:,:,:)
        real,         allocatable :: subvol_res(:)
        integer,      allocatable :: nfcomps(:)
        integer :: filtsz, ldim(3), i, j, h, k, l, funit, io_stat, npix, flims(3,2), sh
        integer :: find_hres, find_lres, ithr, lb(3), ub(3), min_find
        real    :: smpd
        logical :: map2filt_present, pphran
        ! check input
        if( even%is_ft()   ) THROW_HARD('even vol FTed; nonuniform_fscTVfilt')
        if( odd%is_ft()    ) THROW_HARD('odd  vol FTed; nonuniform_fscTVfilt')
        if( mskimg%is_ft() ) THROW_HARD('msk  vol FTed; nonuniform_fscTVfilt')
        map2filt_present = present(map2filt)
        if( map2filt_present )then
            if( map2filt%is_ft() ) THROW_HARD('map2filt vol FTed; nonuniform_fscTVfilt')
        endif
        pphran = .true.
        if( present(phran) ) pphran = phran
        ! set parameters
        ldim    = even%get_ldim()
        if( ldim(3) == 1 ) THROW_HARD('not intended for 2D images; nonuniform_fscTVfilt')
        smpd    = even%get_smpd()
        l_mask  = mskimg%bin2logical()
        npix    = count(l_mask)
        filtsz  = fdim(SUBBOX) - 1
        ! construct
        allocate(subvols_even(nthr_glob), subvols_odd(nthr_glob), subvols_2filt(nthr_glob), nfcomps(filtsz))
        do ithr = 1, nthr_glob
            call subvols_even(ithr)%new(  [SUBBOX,SUBBOX,SUBBOX], smpd, wthreads=.false.)
            call subvols_odd(ithr)%new(   [SUBBOX,SUBBOX,SUBBOX], smpd, wthreads=.false.)
            if( map2filt_present ) call subvols_2filt(ithr)%new( [SUBBOX,SUBBOX,SUBBOX], smpd, wthreads=.false.)
        end do
        subvol_res = subvols_even(1)%get_res()
        ! set # fcomps per shell
        flims   = subvols_even(1)%loop_lims(2)
        nfcomps = 0
        do h = flims(1,1),flims(1,2)
            do k = flims(2,1),flims(2,2)
                do l = flims(3,1),flims(3,2)
                    sh = nint(hyp(real(h),real(k),real(l)))
                    if( sh < 1 .or. sh > filtsz ) cycle
                    nfcomps(sh) = nfcomps(sh) + 1
                end do
            end do
        end do
        if( map2filt_present )then
            call tmpvol%new(ldim, smpd)
        else
            call tmpvol%new(ldim, smpd)
            call tmpvol_odd%new(ldim, smpd)
        endif
        ! determine loop bounds for better load balancing in the following parallel loop
        call bounds_from_mask3D( l_mask, lb, ub )
        ! loop over pixels
        !$omp parallel do collapse(3) default(shared) private(i,j,k,ithr) schedule(dynamic,CHUNKSZ) proc_bind(close)
        do k = lb(3), ub(3)
            do j = lb(2), ub(2)
                do i = lb(1), ub(1)
                    if( .not.l_mask(i,j,k) ) cycle
                    ithr = omp_get_thread_num() + 1
                    call set_subvols_msk_fft_fscTVfilt_ifft([i,j,k], ithr)
                    if( map2filt_present )then
                        call tmpvol%set_rmat_at(i,j,k, subvols_2filt(ithr)%get_rmat_at(CPIX,CPIX,CPIX))
                    else
                        call tmpvol%set_rmat_at(i,j,k, subvols_even(ithr)%get_rmat_at(CPIX,CPIX,CPIX))
                        call tmpvol_odd%set_rmat_at(i,j,k, subvols_odd(ithr)%get_rmat_at(CPIX,CPIX,CPIX))
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        if( map2filt_present )then
            call map2filt%copy(tmpvol)
        else
            call even%copy(tmpvol)
            call odd%copy(tmpvol_odd)
        endif
        ! destruct
        do ithr = 1, nthr_glob
            call subvols_even(ithr)%kill
            call subvols_odd(ithr)%kill
            if( map2filt_present ) call subvols_2filt(ithr)%kill
        end do
        if( map2filt_present )then
            call tmpvol%kill
        else
            call tmpvol%kill
            call tmpvol_odd%kill
        endif
        deallocate(l_mask)

    contains

        subroutine set_subvols_msk_fft_fscTVfilt_ifft( loc, ithr )
            integer, intent(in) :: loc(3), ithr
            integer :: lb(3), ub(3), i, j, k, ii, jj, kk, isub, jsub, ksub
            real    :: corrs(filtsz), filt(filtsz)
            call subvols_even(ithr)%zero_and_unflag_ft
            call subvols_odd(ithr)%zero_and_unflag_ft
            if( map2filt_present ) call subvols_2filt(ithr)%zero_and_unflag_ft
            ! set bounds
            lb = loc - LSHIFT
            ub = loc + RSHIFT
            ! extract components
            kk  = -LSHIFT
            do k = lb(3), ub(3)
                if( k < 1 .and. k > ldim(3) ) cycle
                jj = -LSHIFT
                do j = lb(2), ub(2)
                    if( j < 1 .and. j > ldim(2) ) cycle
                    ii = -LSHIFT
                    do i = lb(1), ub(1)
                        if( i < 1 .and. i > ldim(1) ) cycle
                        isub = ii + LSHIFT + 1
                        jsub = jj + LSHIFT + 1
                        ksub = kk + LSHIFT + 1
                        call subvols_even(ithr)%set_rmat_at(  isub,jsub,ksub, even%get_rmat_at(i,j,k))
                        call subvols_odd(ithr)%set_rmat_at(   isub,jsub,ksub, odd%get_rmat_at(i,j,k))
                        if( map2filt_present ) call subvols_2filt(ithr)%set_rmat_at( isub,jsub,ksub, map2filt%get_rmat_at(i,j,k))
                        ii = ii + 1
                    enddo
                    jj = jj + 1
                enddo
                kk = kk + 1
            enddo
            ! mask
            call subvols_even(ithr)%mask(  SUBMSK, 'soft')
            call subvols_odd(ithr)%mask(   SUBMSK, 'soft')
            if( map2filt_present ) call subvols_2filt(ithr)%mask( SUBMSK, 'soft')
            ! fft
            call subvols_even(ithr)%fft
            call subvols_odd(ithr)%fft
            if( map2filt_present ) call subvols_2filt(ithr)%fft
            ! fsc
            call subvols_even(ithr)%fsc(subvols_odd(ithr), corrs)
            ! calculate TV filter
            call fsc2TVfilt_fast(corrs, nfcomps, filt)
            ! randomize phases below noise power
            if( pphran ) call subvols_even(ithr)%ran_phases_below_noise_power(subvols_odd(ithr), 0.5)
            ! apply TV filter
            if( map2filt_present )then
                if( l_match_filt )then
                    call subvols_2filt(ithr)%shellnorm_and_apply_filter(filt)
                else
                    call subvols_2filt(ithr)%apply_filter(filt)
                endif
                ! back to real space
                call subvols_even(ithr)%zero_and_unflag_ft
                call subvols_odd(ithr)%zero_and_unflag_ft
                call subvols_2filt(ithr)%ifft
            else
                if( l_match_filt )then
                    call subvols_even(ithr)%shellnorm_and_apply_filter(filt)
                    call subvols_odd(ithr)%shellnorm_and_apply_filter(filt)
                else
                    call subvols_even(ithr)%apply_filter(filt)
                    call subvols_odd(ithr)%apply_filter(filt)
                endif
                call subvols_even(ithr)%ifft
                call subvols_odd(ithr)%ifft
            endif
        end subroutine set_subvols_msk_fft_fscTVfilt_ifft

    end subroutine nonuniform_fscTVfilt

    ! filtering of the map according to local resolution estimates
    subroutine local_res_lp( locres_finds, img2filter )
        use simple_image, only: image
        integer,      intent(in)    :: locres_finds(:,:,:)
        class(image), intent(inout) :: img2filter
        real, allocatable :: rmat_filt(:,:,:), rmat_lp(:,:,:)
        type(image) :: resimg
        integer     :: ldim_finds(3), ldim(3), k, kstart, kstop
        real        :: smpd, lp
        ! sanity checks
        if( any(locres_finds == 0 ) ) THROW_HARD('zero Fourier indices not allowed in locres_finds; local_res_lp')
        ldim_finds(1) = size(locres_finds,1)
        ldim_finds(2) = size(locres_finds,2)
        ldim_finds(3) = size(locres_finds,3)
        ldim = img2filter%get_ldim()
        if( .not. all(ldim_finds == ldim) ) THROW_HARD('nonconforming dims of img2filter and locres_finds; local_res_lp')
        ! fetch params
        smpd   = img2filter%get_smpd()
        kstart = minval(locres_finds)
        kstop  = maxval(locres_finds)
        ! to avoid allocation in the loop
        call resimg%new(ldim, smpd)
        allocate(rmat_filt(ldim(1),ldim(2),ldim(3)), rmat_lp(ldim(1),ldim(2),ldim(3)), source=0.)
        ! loop over shells
        do k=kstart,kstop
            if( any(locres_finds == k ) )then
                lp = calc_lowpass_lim(k, ldim(1), smpd)
                call resimg%copy(img2filter)
                call resimg%bp(0.,lp, width=6.0)
                call resimg%get_rmat_sub(rmat_lp)
                where( locres_finds == k ) rmat_filt = rmat_lp
            endif
        enddo
        call img2filter%set_rmat(rmat_filt,.false.)
        call resimg%kill
    end subroutine local_res_lp

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
