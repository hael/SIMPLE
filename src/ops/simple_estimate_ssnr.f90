! spectral signal-to-noise ratio estimation routines
module simple_estimate_ssnr
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
implicit none

public :: fsc2ssnr, fsc2optlp, fsc2optlp_sub, ssnr2fsc, ssnr2optlp, subsample_optlp
public :: acc_dose2filter, dose_weight, nonuniform_lp, local_res, local_res_lp
public :: plot_fsc
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
        allocate( filt(nyq),stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("in simple_estimate_ssnr::fsc2optlp filt ",alloc_stat)
        filt = 0.
        where( corrs > 0. )     filt = 2. * corrs / (corrs + 1.)
        where( filt  > 0.9999 ) filt = 0.99999
    end function fsc2optlp

    !> \brief  converts the FSC to the optimal low-pass filter
    subroutine fsc2optlp_sub( filtsz, corrs, filt )
        integer, intent(in)  :: filtsz        !< sz of filter
        real,    intent(in)  :: corrs(filtsz) !< fsc plot (correlations)
        real,    intent(out) :: filt(filtsz)  !< output filter coefficients
        filt = 0.
        where( corrs > 0. )     filt = 2. * corrs / (corrs + 1.)
        where( filt  > 0.9999 ) filt = 0.99999
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
        allocate( corrs(nyq),stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("in simple_estimate_ssnr::ssnr2fsc corrs ",alloc_stat)
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
        allocate( w(nyq),stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("in simple_estimate_ssnr::ssnr2optlp w ",alloc_stat)
        do k=1,nyq
            w(k) = ssnr(k) / (ssnr(k) + 1.)
        end do
    end function ssnr2optlp

    !> DOSE FILTERING (Grant, Grigorieff eLife 2015)
    !! input is template image, accumulative dose (in e/A2) and acceleration voltage
    !!         output is filter coefficients
    !! \f$  \mathrm{dose}_\mathrm{acc} = \int^{N}_{1} \mathrm{dose\_weight}(a,F,V),\ n_\mathrm{e}/\si{\angstrom\squared}  \f$
    !! \param acc_dose accumulative dose (in \f$n_\mathrm{e}^- per \si{\angstrom\squared}\f$)
    function acc_dose2filter( img, acc_dose, kV, filtsz ) result( filter )
        use simple_image,   only: image
        type(image), intent(in) :: img           !< input image
        real,        intent(in) :: acc_dose, kV  !< acceleration voltage
        integer,     intent(in) :: filtsz
        real, allocatable       :: filter(:)
        integer :: find
        allocate(filter(filtsz),stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("simple_estimate_ssnr::acc_dose2filter ",alloc_stat)
        do find=1,filtsz
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
            critical_exp = critical_exp * kV_factor
        else
            THROW_HARD('unsupported kV (acceleration voltage); dose_weight')
        endif
        dose_weight = exp(-acc_dose/(2.0*critical_exp))
    end function dose_weight

    subroutine nonuniform_lp( even, odd, mskimg )
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
        if( even%is_ft()   ) THROW_HARD('even vol FTed; nonuniform_lp')
        if( odd%is_ft()    ) THROW_HARD('odd  vol FTed; nonuniform_lp')
        if( mskimg%is_ft() ) THROW_HARD('msk  vol FTed; nonuniform_lp')
        ! set parameters
        ldim    = even%get_ldim()
        if( ldim(3) == 1 ) THROW_HARD('not intended for 2D images; nonuniform_lp')
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
                    call set_subvols_msk_fft_filter_ifft([i,j,k], ithr)
                    call even%set_rmat_at(i,j,k, subvols_even(ithr)%get_rmat_at(CPIX,CPIX,CPIX))
                    call odd%set_rmat_at(i,j,k, subvols_odd(ithr)%get_rmat_at(CPIX,CPIX,CPIX))
                end do
            end do
        end do
        !$omp end parallel do
        ! destruct
        deallocate(l_mask)

        contains

            subroutine set_subvols_msk_fft_filter_ifft( loc, ithr )
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
                ! filter
                call subvols_even(ithr)%zero_fcomps_below_noise_power(subvols_odd(ithr))
                ! back to real space
                call subvols_even(ithr)%ifft
                call subvols_odd(ithr)%ifft
            end subroutine set_subvols_msk_fft_filter_ifft

    end subroutine nonuniform_lp

    subroutine local_res( even, odd, mskimg, fsc_crit, locres_finds )
        use simple_image, only: image
        class(image),         intent(inout) :: even, odd, mskimg
        real,                 intent(in)    :: fsc_crit
        integer, allocatable, intent(inout) :: locres_finds(:,:,:)
        integer,      parameter   :: SUBBOX=32, LSHIFT=15, RSHIFT=16, CPIX=LSHIFT + 1, CHUNKSZ=20
        real,         parameter   :: SUBMSK=real(SUBBOX)/2. - COSMSKHALFWIDTH - 1.
        type(image),  allocatable :: subvols_even(:) ! one per thread
        type(image),  allocatable :: subvols_odd(:)  ! one per thread
        logical,      allocatable :: l_mask(:,:,:)
        integer :: filtsz, ldim(3), i, j, k, funit, io_stat, npix
        integer :: find_hres, find_lres, ithr, lb(3), ub(3)
        real    :: smpd
        ! check input
        if( even%is_ft()   ) THROW_HARD('even vol FTed; local_res')
        if( odd%is_ft()    ) THROW_HARD('odd  vol FTed; local_res')
        if( mskimg%is_ft() ) THROW_HARD('msk  vol FTed; local_res')
        ! set parameters
        ldim    = even%get_ldim()
        if( ldim(3) == 1 ) THROW_HARD('not intended for 2D images; local_res')
        smpd    = even%get_smpd()
        l_mask  = mskimg%bin2logical()
        npix    = count(l_mask)
        filtsz  = fdim(SUBBOX) - 1
        ! allocate res indices
        if( allocated(locres_finds) ) deallocate(locres_finds)
        allocate(locres_finds(ldim(1),ldim(2),ldim(3)))
        locres_finds = 0
        ! determine loop bounds for better load balancing in the following parallel loop
        call bounds_from_mask3D( l_mask, lb, ub )
        ! loop over pixels
        !$omp parallel do collapse(3) default(shared) private(i,j,k,ithr) schedule(dynamic,CHUNKSZ) proc_bind(close)
        do k = lb(3), ub(3)
            do j = lb(2), ub(2)
                do i = lb(1), ub(1)
                    if( .not.l_mask(i,j,k) ) cycle
                    ithr = omp_get_thread_num() + 1
                    call set_subvols_msk_fft_filter_fsc([i,j,k], ithr)
                end do
            end do
        end do
        !$omp end parallel do
        find_hres = maxval(locres_finds, l_mask)
        find_lres = max(2, minval(locres_finds, l_mask .and. locres_finds > 0))
        write(logfhandle,'(A,1X,F6.2)') '>>> HIGHEST LOCAL RESOLUTION DETERMINED TO:', even%get_lp(find_hres)
        write(logfhandle,'(A,1X,F6.2)') '>>> LOWEST  LOCAL RESOLUTION DETERMINED TO:', even%get_lp(find_lres)
        ! make sure no 0 elements within the mask
        where( l_mask .and. locres_finds == 0 ) locres_finds = filtsz
        ! make sure background at lowest estimated local resolution
        where( locres_finds == 0 ) locres_finds = find_lres
        ! destruct
        deallocate(l_mask)
        ! write output
        call fopen(funit, LOCRESMAP3D_FILE, access='STREAM', action='WRITE',&
            &status='REPLACE', form='UNFORMATTED', iostat=io_stat)
        call fileiochk('local_res, file: '//LOCRESMAP3D_FILE, io_stat)
        write(unit=funit,pos=1) locres_finds
        call fclose(funit)

    contains

        subroutine set_subvols_msk_fft_filter_fsc( loc, ithr )
            integer, intent(in) :: loc(3), ithr
            integer :: lb(3), ub(3), i, j, k, ii, jj, kk, isub, jsub, ksub
            real    :: corrs(filtsz)
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
            ! fsc
            call subvols_even(ithr)%fsc(subvols_odd(ithr), corrs)
            ! Fourier index at fsc_crit
            call get_find_at_crit(filtsz, corrs, fsc_crit, locres_finds(loc(1),loc(2),loc(3)))
            ! back to real space
            call subvols_even(ithr)%zero_and_unflag_ft
            call subvols_odd(ithr)%zero_and_unflag_ft
        end subroutine set_subvols_msk_fft_filter_fsc

    end subroutine local_res

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
        call img2filter%set_rmat(rmat_filt)
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
