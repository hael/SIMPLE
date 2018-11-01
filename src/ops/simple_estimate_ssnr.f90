! spectral signal-to-noise ratio estimation routines
module simple_estimate_ssnr
include 'simple_lib.f08'
implicit none

public :: fsc2ssnr, fsc2optlp, fsc2optlp_sub, ssnr2fsc, ssnr2optlp, acc_dose2filter, dose_weight
public :: local_res2D, local_res2D_lp, local_res, local_res_lp
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

    ! local resolution estimation for class averages
    subroutine local_res2D( even_imgs, odd_imgs, mskrad, corr_thres, locres_finds, half_winsz, fbody )
        use simple_image, only: image
        class(image),               intent(inout) :: even_imgs(:), odd_imgs(:)
        real,                       intent(in)    :: mskrad, corr_thres
        integer, allocatable,       intent(out)   :: locres_finds(:,:,:)
        integer,          optional, intent(in)    :: half_winsz
        character(len=*), optional, intent(in)    :: fbody
        character(len=:), allocatable :: fname_finds
        logical,     allocatable :: l_mask(:,:,:)
        real,        allocatable :: frc(:), res(:)
        complex,     allocatable :: cmat(:,:,:)
        real(dp),    allocatable :: vec1(:), vec2(:)
        type(image), allocatable :: eimgs_copy(:), oimgs_copy(:)
        type(image) :: mskimg
        integer :: hwinsz, filtsz, ldim(3), i, j, k, kind, funit, io_stat, iptcl, cmdims(3)
        integer :: vecsz, find_hres, find_lres, cnt, npix, hwinszsq, nptcls, find_min, find_max
        integer ::  lb(2), ub(2), m, ivec, kk, mm
        real    :: smpd, cc, ccavg, res_fsc05, res_fsc0143, frac_nccs
        nptcls = size(even_imgs)
        if( nptcls /= size(odd_imgs) ) THROW_HARD('Nonconforming image array sizes; local_res2D')
        do iptcl=1,nptcls
            if( .not. even_imgs(iptcl)%is_ft() ) THROW_HARD('even imgs not FTed; local_res2D')
            if( .not. odd_imgs(iptcl)%is_ft()  ) THROW_HARD('odd imgs not FTed; local_res2D')
        enddo
        ! set parameters
        hwinsz   = 3
        if( present(half_winsz) ) hwinsz = half_winsz
        hwinszsq = hwinsz * hwinsz
        ldim     = even_imgs(1)%get_ldim()
        vecsz    = (2 * hwinsz + 1)**2
        filtsz   = even_imgs(1)%get_filtsz()
        smpd     = even_imgs(1)%get_smpd()
        call mskimg%disc(ldim, smpd, mskrad, l_mask)
        call mskimg%kill
        npix     = count(l_mask)
        res      = even_imgs(1)%get_res()
        cmdims   = even_imgs(1)%get_array_shape()
        ! to avoid allocation in the loop
        if( allocated(locres_finds) ) deallocate(locres_finds)
        allocate(locres_finds(nptcls,ldim(1),ldim(2)), vec1(vecsz), vec2(vecsz),&
            &eimgs_copy(nptcls), oimgs_copy(nptcls), cmat(cmdims(1),cmdims(2),cmdims(3)))
        locres_finds = 0
        do iptcl=1,nptcls
            call eimgs_copy(iptcl)%new(ldim, smpd, wthreads=.false.)
            call oimgs_copy(iptcl)%new(ldim, smpd, wthreads=.false.)
        end do
        ! loop over resolution shells
        do kind=1,filtsz
            ccavg = 0.
            cnt   = 0
            !$omp parallel do default(shared) private(iptcl,cmat,j,i,lb,ub,ivec,m,k,mm,kk,vec1,vec2,cc) proc_bind(close) schedule(static) reduction(+:ccavg,cnt)
            do iptcl=1,nptcls
                ! tophat filter out the shell
                call even_imgs(iptcl)%get_cmat_sub(cmat)
                call eimgs_copy(iptcl)%set_cmat(cmat)
                call eimgs_copy(iptcl)%tophat(kind)
                call eimgs_copy(iptcl)%ifft
                call odd_imgs(iptcl)%get_cmat_sub(cmat)
                call oimgs_copy(iptcl)%set_cmat(cmat)
                call oimgs_copy(iptcl)%tophat(kind)
                call oimgs_copy(iptcl)%ifft
                ! neighborhood correlation analysis
                do j=1,ldim(2)
                    do i=1,ldim(1)
                        if( .not. l_mask(i,j,1) ) cycle
                        ! set bounds
                        lb   = [i,j] - hwinsz
                        ub   = [i,j] + hwinsz
                        ! extract components within circle
                        ivec = 0
                        mm   = -hwinsz
                        do m=lb(2),ub(2)
                            if( m >= 1 .and. m <= ldim(2) )then
                                kk = -hwinsz
                                do k=lb(1),ub(1)
                                    if( k >= 1 .and. k <= ldim(1) )then
                                        if( mm*mm + kk*kk <= hwinszsq )then
                                            ivec       = ivec + 1
                                            vec1(ivec) = dble(eimgs_copy(iptcl)%get_rmat_at(k,m,1))
                                            vec2(ivec) = dble(oimgs_copy(iptcl)%get_rmat_at(k,m,1))
                                        endif
                                    endif
                                    kk = kk + 1
                                enddo
                            endif
                            mm = mm + 1
                        enddo
                        ! correlate
                        cc    = pearsn_serial_8(ivec, vec1(:ivec), vec2(:ivec))
                        ccavg = ccavg + cc
                        if( cc >= corr_thres ) locres_finds(iptcl,i,j) = kind
                        if( cc >= 0.5 ) cnt = cnt + 1
                    enddo
                enddo
            end do
            !$omp end parallel do
            ccavg     = ccavg / real(npix * nptcls)
            frac_nccs = real(cnt) / real(npix * nptcls) * 100.
            write(*,'(A,1X,F6.2,1X,A,1X,F7.3,1X,A,1X,F6.2)') '>>> RESOLUTION:', res(kind), '>>> CORRELATION:', ccavg,&
            &'>>> % NEIGH_CCS >= 0.5', frac_nccs
            if( ccavg < 0.01 ) exit
        end do
        ! identify highest and lowest local resolution
        do iptcl=1,nptcls
            find_max = maxval(locres_finds(iptcl,:,:), l_mask(:,:,1))
            if( iptcl == 1 )then
                find_hres = find_max
            else if( find_max > find_hres )then
                find_hres = find_max
            endif
            if( count(l_mask(:,:,1) .and. locres_finds(iptcl,:,:) > 0) > 0 )then
                find_min = minval(locres_finds(iptcl,:,:), l_mask(:,:,1) .and. locres_finds(iptcl,:,:) > 0)
            else
                find_min = 2
            endif
            if( iptcl == 1 )then
                find_lres = max(find_min,2)
            else if( find_min < find_lres )then
                find_lres = max(find_min,2)
            endif
            ! make sure no 0 elements within the mask
            where( l_mask(:,:,1)       .and. locres_finds(iptcl,:,:) == 0 ) locres_finds(iptcl,:,:) = filtsz
            ! make sure background at highest estimated local resolution
            where( .not. l_mask(:,:,1) .and. locres_finds(iptcl,:,:) == 0 ) locres_finds(iptcl,:,:) = max(2,find_max)
        end do
        write(*,'(A,1X,F6.2)') '>>> HIGHEST LOCAL RESOLUTION DETERMINED TO:', even_imgs(1)%get_lp(find_hres)
        write(*,'(A,1X,F6.2)') '>>> LOWEST  LOCAL RESOLUTION DETERMINED TO:', even_imgs(1)%get_lp(find_lres)
        ! destruct
        do iptcl=1,nptcls
            call eimgs_copy(iptcl)%kill
            call oimgs_copy(iptcl)%kill
        end do
        deallocate(l_mask, eimgs_copy, oimgs_copy)
        ! write output
        if( present(fbody) )then
            fname_finds  = trim(fbody)//'_finds.bin'
        else
            fname_finds  = 'locresmap2D_finds.bin'
        endif
        call fopen(funit, fname_finds, access='STREAM', action='WRITE',&
            &status='REPLACE', form='UNFORMATTED', iostat=io_stat)
        call fileiochk('local_res, file: '//fname_finds, io_stat)
        write(unit=funit,pos=1) locres_finds
        call fclose(funit)
    end subroutine local_res2D

    ! filtering of the class averages according to local resolution estimates
    subroutine local_res2D_lp( locres_finds, imgs2filter, shellnorm )
        use simple_image, only: image
        integer,           intent(in)    :: locres_finds(:,:,:)
        class(image),      intent(inout) :: imgs2filter(:)
        logical, optional, intent(in)    :: shellnorm
        real,         allocatable   :: rmats_filt(:,:,:,:), rmats_lp(:,:,:,:), rmats(:,:,:,:)
        type(image),  allocatable   :: resimgs(:)
        integer :: ldim_finds(3), ldim(3), k, kstart, kstop, nptcls, iptcl
        real    :: smpd, lp
        logical :: sshellnorm
        sshellnorm = .false.
        if( present(shellnorm) ) sshellnorm = shellnorm
        ! sanity checks
        if( any(locres_finds == 0 ) ) THROW_HARD('zero Fourier indices not allowed in locres_finds; local_res2D_lp')
        ldim_finds(1) = size(locres_finds,1)
        ldim_finds(2) = size(locres_finds,2)
        ldim_finds(3) = size(locres_finds,3)
        nptcls = size(imgs2filter)
        if( nptcls /= ldim_finds(1) ) THROW_HARD('nonconforming dim 1 of imgs2filter and locres_finds; local_res2D_lp')
        ldim = imgs2filter(1)%get_ldim()
        if( .not. all(ldim_finds(2:3) == ldim(1:2)) ) THROW_HARD('nonconforming dims 2-3 locres_finds vs imgs2filter; local_res2D_lp')
        ! fetch params
        smpd   = imgs2filter(1)%get_smpd()
        kstart = minval(locres_finds)
        kstop  = maxval(locres_finds)
        ! to avoid allocation in the loop
        allocate(resimgs(nptcls), rmats_filt(nptcls,ldim(1),ldim(2),1),&
            &rmats_lp(nptcls,ldim(1),ldim(2),1), rmats(nptcls,ldim(1),ldim(2),1))
        do iptcl=1,nptcls
            call resimgs(iptcl)%new(ldim, smpd, wthreads=.false.)
        end do
        rmats_filt = 0.
        rmats_lp    = 0.
        ! loop over shells
        do k=kstart,kstop
            call progress(k,kstop)
            !$omp parallel do default(shared) private(iptcl,lp) proc_bind(close) schedule(static)
            do iptcl=1,nptcls
                if( any(locres_finds(iptcl,:,:) == k ) )then
                    lp = calc_lowpass_lim(k, ldim(1), smpd)
                    if( sshellnorm ) call imgs2filter(iptcl)%shellnorm
                    call imgs2filter(iptcl)%get_rmat_sub(rmats(iptcl,:,:,:))
                    call resimgs(iptcl)%set_rmat(rmats(iptcl,:,:,:))
                    call resimgs(iptcl)%bp(0., lp, width=6.0)
                    call resimgs(iptcl)%get_rmat_sub(rmats_lp(iptcl,:,:,:))
                    where( locres_finds(iptcl,:,:) == k ) rmats_filt(iptcl,:,:,1) = rmats_lp(iptcl,:,:,1)
                endif
            end do
            !$omp end parallel do
        enddo
        ! set the filtered images and destroy resimgs
        do iptcl=1,nptcls
            call imgs2filter(iptcl)%set_rmat(rmats_filt(iptcl,:,:,:))
            call resimgs(iptcl)%kill
        end do
    end subroutine local_res2D_lp

    ! even odd images assumed to be appropriately masked before calling this routine
    ! mskimg should be a hard mask (real-space)
    subroutine local_res( even, odd, mskimg, corr_thres, locres_finds, half_winsz, fbody )
        use simple_image, only: image
        class(image),               intent(inout) :: even, odd, mskimg
        real,                       intent(in)    :: corr_thres
        integer, allocatable,       intent(out)   :: locres_finds(:,:,:)
        integer,          optional, intent(in)    :: half_winsz
        character(len=*), optional, intent(in)    :: fbody
        character(len=:), allocatable :: fname_finds
        logical,  allocatable :: l_mask(:,:,:)
        real,     allocatable :: fsc(:), res(:)
        integer     :: hwinsz, filtsz, ldim(3), i, j, k, kind, funit, io_stat
        integer     :: vecsz, find_hres, find_lres, cnt, npix, hwinszsq
        type(image) :: ecopy, ocopy
        real        :: smpd, cc, ccavg, res_fsc05, res_fsc0143, frac_nccs
        ! check input
        if( .not. even%is_ft() ) THROW_HARD('even vol not FTed; local_res')
        if( .not. odd%is_ft()  ) THROW_HARD('odd vol not FTed; local_res')
        ! set parameters
        hwinsz   = 3
        if( present(half_winsz) ) hwinsz = half_winsz
        hwinszsq = hwinsz * hwinsz
        ldim     = even%get_ldim()
        if( ldim(3) == 1 )  THROW_HARD('not intended for 2D images; local_res')
        vecsz = (2 * hwinsz + 1)**3
        filtsz = even%get_filtsz()
        smpd   = even%get_smpd()
        l_mask = mskimg%bin2logical()
        npix   = count(l_mask)
        res    = even%get_res()
        ! to avoid allocation in the loop
        call ecopy%new(ldim, smpd)
        call ocopy%new(ldim, smpd)
        if( allocated(locres_finds) ) deallocate(locres_finds)
        allocate(locres_finds(ldim(1),ldim(2),ldim(3)), fsc(filtsz))
        locres_finds = 0
        fsc          = 0.
        ! loop over resolution shells
        do kind=1,filtsz
            ! tophat filter out the shell
            call ecopy%copy(even)
            call ocopy%copy(odd)
            call ecopy%tophat(kind)
            call ocopy%tophat(kind)
            call ecopy%ifft
            call ocopy%ifft
            ! neighborhood correlation analysis
            ccavg = 0.
            cnt   = 0
            ! parallel section needs to start here as the above steps are threaded as well
            !$omp parallel do collapse(3) default(shared) private(i,j,k,cc)&
            !$omp schedule(static) proc_bind(close) reduction(+:ccavg,cnt)
            do k=1,ldim(3)
                do j=1,ldim(2)
                    do i=1,ldim(1)
                        if( .not. l_mask(i,j,k) ) cycle
                        cc    = neigh_cc([i,j,k])
                        ccavg = ccavg + cc
                        if( cc >= corr_thres ) locres_finds(i,j,k) = kind
                        if( cc >= 0.5 ) cnt = cnt + 1
                    end do
                end do
            end do
            !$omp end parallel do
            ccavg     = ccavg / real(npix)
            fsc(kind) = ccavg
            frac_nccs = real(cnt) / real(npix) * 100.
            write(*,'(A,1X,F6.2,1X,A,1X,F7.3,1X,A,1X,F6.2)') '>>> RESOLUTION:', res(kind), '>>> CORRELATION:', fsc(kind),&
            &'>>> % NEIGH_CCS >= 0.5', frac_nccs
            if( frac_nccs < 1. .or. ccavg < 0.01 ) exit ! save the compute
        end do
        call get_resolution(fsc, res, res_fsc05, res_fsc0143)
        write(*,'(A,1X,F6.2)') '>>> GLOBAL RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
        write(*,'(A,1X,F6.2)') '>>> GLOBAL RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
        find_hres = maxval(locres_finds, l_mask)
        find_lres = max(2, minval(locres_finds, l_mask .and. locres_finds > 0))
        write(*,'(A,1X,F6.2)') '>>> HIGHEST LOCAL RESOLUTION DETERMINED TO:', even%get_lp(find_hres)
        write(*,'(A,1X,F6.2)') '>>> LOWEST  LOCAL RESOLUTION DETERMINED TO:', even%get_lp(find_lres)
        ! make sure no 0 elements within the mask
        where( l_mask .and. locres_finds == 0 ) locres_finds = filtsz
        ! make sure background at lowest estimated local resolution
        where( locres_finds == 0 ) locres_finds = find_lres
        ! destruct
        deallocate(l_mask)
        call ecopy%kill
        call ocopy%kill
        ! write output
        if( present(fbody) )then
            fname_finds  = trim(fbody)//'_finds.bin'
        else
            fname_finds  = 'locresmap_finds.bin'
        endif
        call fopen(funit, fname_finds, access='STREAM', action='WRITE',&
            &status='REPLACE', form='UNFORMATTED', iostat=io_stat)
        call fileiochk('local_res, file: '//fname_finds, io_stat)
        write(unit=funit,pos=1) locres_finds
        call fclose(funit)

        contains

            function neigh_cc( loc ) result( cc )
                integer, intent(in) :: loc(3)
                integer  :: lb(3), ub(3), i, j, k, cnt, ii, jj, kk
                real     :: cc
                real(dp) :: vec1(vecsz), vec2(vecsz)
                ! set bounds
                lb = loc - hwinsz
                ub = loc + hwinsz
                ! extract components wihin sphere
                cnt = 0
                kk  = -hwinsz
                do k=lb(3),ub(3)
                    if( k >= 1 .and. k <= ldim(3) )then
                        jj = -hwinsz
                        do j=lb(2),ub(2)
                            if( j >= 1 .and. j <= ldim(2) )then
                                ii = -hwinsz
                                do i=lb(1),ub(1)
                                    if( i >= 1 .and. i <= ldim(1) )then
                                        if( kk*kk + jj*jj + ii*ii <= hwinszsq )then
                                            cnt       = cnt + 1
                                            vec1(cnt) = dble(ecopy%get_rmat_at(i,j,k))
                                            vec2(cnt) = dble(ocopy%get_rmat_at(i,j,k))
                                        endif
                                    endif
                                    ii = ii + 1
                                enddo
                            endif
                            jj = jj + 1
                        enddo
                    endif
                    kk = kk + 1
                enddo
                ! correlate
                cc = pearsn_serial_8(cnt, vec1(:cnt), vec2(:cnt))
            end function neigh_cc

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

end module simple_estimate_ssnr
