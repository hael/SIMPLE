! optimization(search)-based filter (uniform/nonuniform)
module simple_opt_filter
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,      only: image, image_ptr
use simple_parameters, only: params_glob
implicit none
#include "simple_local_flags.inc"

public :: nonuni_filt2D_sub, nonuni_filt2D, nonuni_filt3D, butterworth_filter, uni_filt2D, uni_filt2D_sub, uni_filt3D, exponential_reg
public :: test_filt2D
private

interface butterworth_filter
    module procedure butterworth_filter_1, butterworth_filter_2, butterworth_filter_3, butterworth_filter_4
end interface butterworth_filter

type optfilt2Dvars
    real, allocatable :: cur_fil(:)
    type(image)       :: odd_copy_rmat, even_copy_rmat
    type(image)       :: odd_copy_cmat, even_copy_cmat
    type(image)       :: odd_filt, even_filt, diff_img_odd, diff_img_even, diff_img_opt_odd, diff_img_opt_even
    real, allocatable :: frc(:)
    integer           :: best_ind_odd, best_ind_even ! for uniform filter
    integer           :: lplim_hres
    logical           :: have_mask = .false.
end type optfilt2Dvars

logical :: TEST = .false.

contains

    subroutine nonuni_filt2D_sub( even, odd, mask )
        use simple_class_frcs, only: class_frcs
        class(image),           intent(inout) :: even(:), odd(:)
        class(image), optional, intent(inout) :: mask(:)
        character(len=:),       allocatable   :: frcs_fname
        type(optfilt2Dvars),    allocatable   :: optf2Dvars(:)
        real,                   allocatable   :: frc(:)
        type(class_frcs) :: clsfrcs
        type(image)      :: weights_img
        type(kbinterpol) :: kb
        real             :: smpd, lpstart, lp, val, rk, dk, scale, s
        integer          :: iptcl, box, box_pd, filtsz, ldim(3), ldim_pd(3), smooth_ext
        integer          :: nptcls, hpind_fsc, find, m, n, k, lk, uk
        logical          :: lpstart_fallback, l_phaseplate, have_mask
        write(logfhandle,'(A)') '>>> 2D NONUNIFORM FILTERING'
        ! init
        filtsz       = even(1)%get_filtsz()
        smooth_ext   = params_glob%smooth_ext
        ldim         = even(1)%get_ldim()
        ldim(3)      = 1 ! because we operate on stacks
        box          = ldim(1)
        box_pd       = find_larger_magic_box(box + 2 * smooth_ext)
        ldim_pd      = [box_pd, box_pd, 1]
        frcs_fname   = trim(params_glob%frcs)
        smpd         = even(1)%get_smpd()
        nptcls       = size(even)
        lpstart      = params_glob%lpstart
        hpind_fsc    = params_glob%hpind_fsc
        l_phaseplate = params_glob%l_phaseplate
        have_mask    = present(mask)
        ! retrieve FRCs
        lpstart_fallback = .false.
        if( file_exists(frcs_fname) )then
            call clsfrcs%read(frcs_fname)
            if( clsfrcs%get_filtsz().ne.even(1)%get_filtsz() )then
                write(logfhandle,*) 'img filtsz:  ', even(1)%get_filtsz()
                write(logfhandle,*) 'frcs filtsz: ', clsfrcs%get_filtsz()
                THROW_HARD('Inconsistent filter dimensions; nonuni_filt2D_sub')
            endif
        else
            THROW_WARN('Class average FRCs file '//frcs_fname//' does not exist, falling back on lpstart: '//real2str(lpstart))
            lpstart_fallback = .true.
        endif
        filtsz = clsfrcs%get_filtsz()
        ! allocate
        allocate(optf2Dvars(nptcls), frc(filtsz))
        frc = 0.
        ! calculate high-res low-pass limits
        if( lpstart_fallback )then
            optf2Dvars(:)%lplim_hres = calc_fourier_index(lpstart, box_pd, smpd)
        else
            ! the below required to retrieve the right Fourier index limit when we are padding
            do iptcl = 1, nptcls
                call clsfrcs%frc_getter(iptcl, hpind_fsc, l_phaseplate, frc)
                find = get_lplim_at_corr(frc, 0.1)          ! little overshoot, resolution is limited in polarft_corrcalc anyway
                lp   = calc_lowpass_lim(find, box, smpd)    ! resolution limit in original image
                optf2Dvars(iptcl)%lplim_hres = calc_fourier_index(lp, box_pd, smpd) ! this is the Fourier index limit for the padded image
                if( TEST )then
                    ! frc interpolation to padded dimensions
                    allocate(optf2Dvars(iptcl)%frc(1:box_pd/2),source=0.)
                    scale = real(box) / real(box_pd)
                    do k = 1,box_pd/2
                        rk = scale * real(k)
                        lk = floor(rk)
                        dk = rk - real(lk)
                        uk = lk+1
                        if( lk == 0 )then
                            optf2Dvars(iptcl)%frc(k) = (1.0-dk)
                        else
                            optf2Dvars(iptcl)%frc(k) = (1.0-dk)*frc(lk)
                        endif
                        if( uk < box/2 ) optf2Dvars(iptcl)%frc(k) = optf2Dvars(iptcl)%frc(k) + dk*frc(uk)
                    enddo
                endif

            end do
        endif
        ! fill up optf2Dvars struct
        !$omp parallel do default(shared) private(iptcl) schedule(static) proc_bind(close)
        do iptcl = 1, nptcls
            call even(iptcl)%pad_mirr(ldim_pd)
            call odd(iptcl)%pad_mirr(ldim_pd)
            call optf2Dvars(iptcl)%even_filt%copy(even(iptcl))
            call optf2Dvars(iptcl)%odd_filt %copy( odd(iptcl))
            call optf2Dvars(iptcl)%diff_img_odd    %new(ldim_pd, smpd, .false.)
            call optf2Dvars(iptcl)%diff_img_even   %new(ldim_pd, smpd, .false.)
            call optf2Dvars(iptcl)%diff_img_opt_odd %new(ldim_pd, smpd, .false.)
            call optf2Dvars(iptcl)%diff_img_opt_even%new(ldim_pd, smpd, .false.)
            call optf2Dvars(iptcl)%odd_copy_rmat%copy(odd(iptcl))
            call optf2Dvars(iptcl)%odd_copy_cmat%copy(odd(iptcl))
            call optf2Dvars(iptcl)%odd_copy_cmat%fft
            call optf2Dvars(iptcl)%even_copy_rmat%copy(even(iptcl))
            call optf2Dvars(iptcl)%even_copy_cmat%copy(even(iptcl))
            call optf2Dvars(iptcl)%even_copy_cmat%fft
            optf2Dvars(iptcl)%have_mask = .false.
            allocate(optf2Dvars(iptcl)%cur_fil(box_pd), source=0.)
            if( have_mask )then
                if( mask(iptcl)%nforeground() > 0 )then
                    optf2Dvars(iptcl)%have_mask = .true.
                    call mask(iptcl)%pad_inplace(ldim_pd)
                endif
            endif
        end do
        !$omp end parallel do
        ! make weight image for diff convolution
        call weights_img%new(ldim_pd, smpd, .false.)
        call weights_img%zero_and_unflag_ft()
        if( TEST )then
            kb = kbinterpol(1.5,1.)
            s  = 0.0
            do m = -min(2*params_glob%smooth_ext,box_pd/2-1), min(2*params_glob%smooth_ext,box_pd/2-1)
                do n = -min(2*params_glob%smooth_ext,box_pd/2-1), min(2*params_glob%smooth_ext,box_pd/2-1)
                    val = kb%apod(hyp(m,n)/real(params_glob%smooth_ext)) / kb%apod(0.)
                    s   = s + val
                    call weights_img%set_rmat_at(box_pd/2+m+1, box_pd/2+n+1, 1, val)
                enddo
            enddo
            call weights_img%div(s)
        else
            do m = -params_glob%smooth_ext, params_glob%smooth_ext
                do n = -params_glob%smooth_ext, params_glob%smooth_ext
                    val = -hyp(m,n) / (params_glob%smooth_ext + 1) + 1.
                    if( val > 0 ) call weights_img%set_rmat_at(box_pd/2+m+1, box_pd/2+n+1, 1, val)
                enddo
            enddo
        endif
        call weights_img%fft()
        ! filter
        if( have_mask )then
            !$omp parallel do default(shared) private(iptcl) schedule(static) proc_bind(close)
            do iptcl = 1, nptcls
                call nonuni_filt2D_masked(odd(iptcl), even(iptcl), mask(iptcl), weights_img, optf2Dvars(iptcl))
            enddo
            !$omp end parallel do
        else
            !$omp parallel do default(shared) private(iptcl) schedule(static) proc_bind(close)
            do iptcl = 1, nptcls
                if( TEST )then
                    call test_filt2D(odd(iptcl), even(iptcl), weights_img, optf2Dvars(iptcl))
                else
                    call nonuni_filt2D(odd(iptcl), even(iptcl), weights_img, optf2Dvars(iptcl))
                endif
            enddo
            !$omp end parallel do
        endif
        ! destruct
        call clsfrcs%kill
        call weights_img%kill
        do iptcl = 1, nptcls
            call optf2Dvars(iptcl)%odd_copy_rmat%kill
            call optf2Dvars(iptcl)%even_copy_rmat%kill
            call optf2Dvars(iptcl)%odd_copy_cmat%kill
            call optf2Dvars(iptcl)%even_copy_cmat%kill
            call optf2Dvars(iptcl)%even_filt%kill
            call optf2Dvars(iptcl)%odd_filt%kill
            call optf2Dvars(iptcl)%diff_img_odd%kill
            call optf2Dvars(iptcl)%diff_img_even%kill
            call optf2Dvars(iptcl)%diff_img_opt_odd%kill
            call optf2Dvars(iptcl)%diff_img_opt_even%kill
            if( allocated(optf2Dvars(iptcl)%frc) )     deallocate(optf2Dvars(iptcl)%frc)
            if( allocated(optf2Dvars(iptcl)%cur_fil) ) deallocate(optf2Dvars(iptcl)%cur_fil)
        enddo
        do iptcl = 1, nptcls
            call even(iptcl)%clip_inplace(ldim)
            call odd( iptcl)%clip_inplace(ldim)
            if( have_mask )then
                if( optf2Dvars(iptcl)%have_mask )then
                    call mask(iptcl)%clip_inplace(ldim)
                endif
            endif
        enddo
    end subroutine nonuni_filt2D_sub

    ! Compute the value of the Butterworth transfer function of order n(th)
    ! at a given frequency s, with the cut-off frequency fc
    ! SOURCE :
    ! https://en.wikipedia.org/wiki/Butterworth_filter
    pure function butterworth(s, n, fc) result(val)
        real   , intent(in)  :: s
        integer, intent(in)  :: n
        real   , intent(in)  :: fc
        real                 :: val
        real,    parameter :: AN(11,10) = reshape((/ 1., 1.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    , 0.    , 0.,&
                                                    &1., 1.4142,  1.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    , 0.    , 0.,&
                                                    &1., 2.    ,  2.    ,  1.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    , 0.    , 0.,&
                                                    &1., 2.6131,  3.4142,  2.6131,  1.    ,  0.    ,  0.    ,  0.    ,  0.    , 0.    , 0.,&
                                                    &1., 3.2361,  5.2361,  5.2361,  3.2361,  1.    ,  0.    ,  0.    ,  0.    , 0.    , 0.,&
                                                    &1., 3.8637,  7.4641,  9.1416,  7.4641,  3.8637,  1.    ,  0.    ,  0.    , 0.    , 0.,&
                                                    &1., 4.4940, 10.0978, 14.5918, 14.5918, 10.0978,  4.4940,  1.    ,  0.    , 0.    , 0.,&
                                                    &1., 5.1258, 13.1371, 21.8462, 25.6884, 21.8462, 13.1371,  5.1258,  1.    , 0.    , 0.,&
                                                    &1., 5.7588, 16.5817, 31.1634, 41.9864, 41.9864, 31.1634, 16.5817,  5.7588, 1.    , 0.,&
                                                    &1., 6.3925, 20.4317, 42.8021, 64.8824, 74.2334, 64.8824, 42.8021, 20.4317, 6.3925, 1. /),&
                                                    &(/11,10/))
        complex, parameter :: J = (0, 1) ! Complex identity: j = sqrt(-1)
        complex :: Bn, Kn                ! Normalized Butterworth polynomial, its derivative and its reciprocal
        complex :: js                    ! frequency is multiplied by the complex identity j
        integer :: k
        Bn  = (0., 0.)
        if (s/fc < 100) then
            js  = J*s/fc
            do k = 0, n
                Bn  = Bn + AN(k+1,n)*js**k
            end do
            Kn  = 1/Bn
            val = cabs(Kn)
        else
            val = epsilon(val)
        endif
    end function butterworth

    subroutine butterworth_filter_1(img, cutoff_find, cur_fil)
        class(image), intent(inout) :: img
        integer,      intent(in)    :: cutoff_find
        real,         intent(inout) :: cur_fil(:)
        integer,      parameter     :: BW_ORDER = 8
        integer :: freq_val
        do freq_val = 1, size(cur_fil)
            cur_fil(freq_val) = butterworth(real(freq_val-1), BW_ORDER, real(cutoff_find))
        enddo
        call img%apply_filter(cur_fil)
    end subroutine butterworth_filter_1

    subroutine butterworth_filter_2(cutoff_find, cur_fil)
        integer,      intent(in)    :: cutoff_find
        real,         intent(inout) :: cur_fil(:)
        integer,      parameter     :: BW_ORDER = 8
        integer :: freq_val
        do freq_val = 1, size(cur_fil)
            cur_fil(freq_val) = butterworth(real(freq_val-1), BW_ORDER, real(cutoff_find))
        enddo
    end subroutine butterworth_filter_2

    !< notch filter squeezed between two cut_off frequencies
    subroutine butterworth_filter_3(c1, c2, cur_fil)
        integer,      intent(in)    :: c1, c2           ! two cut-off frequencies
        real,         intent(inout) :: cur_fil(:)
        integer,      parameter     :: BW_ORDER = 8
        integer :: freq_val
        do freq_val = 1, size(cur_fil)
            cur_fil(freq_val) =        butterworth(real(freq_val-1), BW_ORDER, real(c2)) * &
                                &( 1 - butterworth(real(freq_val-1), BW_ORDER, real(c1)) )
        enddo
    end subroutine butterworth_filter_3

    subroutine butterworth_filter_4(cutoff_find, kfromto, cur_fil)
        integer,      intent(in)    :: cutoff_find
        integer,      intent(in)    :: kfromto(2)
        real,         intent(inout) :: cur_fil(kfromto(1):kfromto(2))
        integer,      parameter     :: BW_ORDER = 8
        integer :: freq_val
        do freq_val = kfromto(1), kfromto(2)
            cur_fil(freq_val) = butterworth(real(freq_val-1), BW_ORDER, real(cutoff_find))
        enddo
    end subroutine butterworth_filter_4

    subroutine test_filt2D(odd, even, weights_img, optf2Dvars)
        class(image),        intent(inout) :: odd, even, weights_img
        type(optfilt2Dvars), intent(inout) :: optf2Dvars
        real,             parameter :: alpha_max = 6.0
        real,             parameter :: beta      = 0.002
        type(image_ptr)             :: pdiff_opt_odd, pdiff_opt_even, pdiff_odd, pdiff_even, pweights
        real                        :: frc_opt(odd%get_filtsz())
        real(kind=c_float), pointer :: rmat_odd(:,:,:), rmat_even(:,:,:), rmat_odd_filt(:,:,:), rmat_even_filt(:,:,:)
        real            :: a, alpha
        integer         :: k,l, box, ldim(3), it, ext, lb(3), ub(3), n, hnsearch, find
        ! init
        ldim = odd%get_ldim()
        box  = ldim(1)
        if( ldim(3) > 1 ) THROW_HARD('This test_filt2D is strictly for 2D!')
        ext  = params_glob%smooth_ext
        lb   = (/ ext+1  , ext+1  , 1/)
        ub   = (/ box-ext, box-ext, 1/)
        n    = (box-2*ext)**2
        hnsearch = merge(params_glob%nsearch/2, (params_glob%nsearch+1)/2, is_even(params_glob%nsearch))
        call weights_img%get_mat_ptrs(pweights)
        call optf2Dvars%diff_img_odd %get_mat_ptrs(pdiff_odd)
        call optf2Dvars%diff_img_even%get_mat_ptrs(pdiff_even)
        call optf2Dvars%diff_img_opt_odd %get_mat_ptrs(pdiff_opt_odd)
        call optf2Dvars%diff_img_opt_even%get_mat_ptrs(pdiff_opt_even)
        pdiff_opt_odd %rmat = huge(alpha)
        pdiff_opt_even%rmat = huge(alpha)
        call odd%get_rmat_ptr(rmat_odd)
        call even%get_rmat_ptr(rmat_even)
        call optf2Dvars%odd_filt%get_rmat_ptr(rmat_odd_filt)
        call optf2Dvars%even_filt%get_rmat_ptr(rmat_even_filt)
        rmat_odd  = 0.
        rmat_even = 0.
        ! bounding frc
        where( optf2Dvars%frc > 0.99999 ) optf2Dvars%frc = 0.99999
        where( optf2Dvars%frc < 0.00001 ) optf2Dvars%frc = 0.0
        ! main loop 
        do it = -hnsearch,hnsearch
            ! scaling factor
            if( it == 0 )then
                alpha = 1.0
            else if( it > 0 )then
                alpha =         alpha_max * real(it) / real(hnsearch)
            else
                alpha = -1.0 / (alpha_max * real(it) / real(hnsearch))
            endif
            optf2Dvars%cur_fil = 2.0*optf2Dvars%frc / (1.0-optf2Dvars%frc)      ! SSNR
            optf2Dvars%cur_fil = alpha * optf2Dvars%cur_fil                     ! SSNR scaling
            frc_opt            = optf2Dvars%cur_fil / (optf2Dvars%cur_fil+2.0)  ! scaled FRC
            find               = get_lplim_at_corr(frc_opt,0.143)               ! scaled resolution limit
            optf2Dvars%cur_fil = optf2Dvars%cur_fil / (optf2Dvars%cur_fil+1.0)  ! modulated optimal filter
            ! butterworth-like falloff
            do k = find+1,odd%get_filtsz(),1
                optf2Dvars%cur_fil(k) = optf2Dvars%cur_fil(k) * butterworth(real(k),8,real(find)) / butterworth(real(find),8,real(find))
            enddo
            ! filtering
            call optf2Dvars%odd_filt%copy_fast(optf2Dvars%odd_copy_cmat)
            call optf2Dvars%even_filt%copy_fast(optf2Dvars%even_copy_cmat)
            call optf2Dvars%odd_filt%apply_filter_serial(optf2Dvars%cur_fil)
            call optf2Dvars%even_filt%apply_filter_serial(optf2Dvars%cur_fil)
            ! square difference
            call optf2Dvars%even_filt%ifft
            call optf2Dvars%odd_filt%ifft
            call optf2Dvars% odd_filt%sqeuclid_matrix(optf2Dvars%even_copy_rmat, optf2Dvars%diff_img_odd)
            call optf2Dvars%even_filt%sqeuclid_matrix(optf2Dvars% odd_copy_rmat, optf2Dvars%diff_img_even)
            ! neighbourhood smoothing
            call optf2Dvars%diff_img_odd%fft
            call optf2Dvars%diff_img_even%fft
            pdiff_odd %cmat = pdiff_odd %cmat * pweights%cmat
            pdiff_even%cmat = pdiff_even%cmat * pweights%cmat
            call optf2Dvars%diff_img_odd%ifft
            call optf2Dvars%diff_img_even%ifft
            ! deviation penalty
            if( alpha > 1. )then
                a = beta * (alpha-1.)**2 / real(n)
                call optf2Dvars%diff_img_odd%add( a )
                call optf2Dvars%diff_img_even%add( a )
            endif
            ! pixel-wise optimization
            do l = lb(2),ub(2)
                do k = lb(1),ub(1)
                    if( pdiff_odd%rmat(k,l,1) < pdiff_opt_odd%rmat(k,l,1) )then
                        rmat_odd(          k,l,1) = rmat_odd_filt( k,l,1)
                        pdiff_opt_odd%rmat(k,l,1) = pdiff_odd%rmat(k,l,1)
                    endif
                    if( pdiff_even%rmat(k,l,1) < pdiff_opt_even%rmat(k,l,1) )then
                        rmat_even(          k,l,1) = rmat_even_filt(k,l,1)
                        pdiff_opt_even%rmat(k,l,1) = pdiff_even%rmat(k,l,1)
                    endif
                enddo
            enddo
        enddo
    end subroutine test_filt2D

    subroutine nonuni_filt2D(odd, even, weights_img, optf2Dvars)
        class(image),        intent(inout) :: odd, even, weights_img
        type(optfilt2Dvars), intent(inout) :: optf2Dvars
        real(kind=c_float),        pointer :: rmat_odd(:,:,:), rmat_even(:,:,:), rmat_odd_filt(:,:,:), rmat_even_filt(:,:,:)
        integer         :: k,l, box, ldim(3), find_start, find_stop, iter_no, ext, cutoff_find, lb(3), ub(3)
        real            :: find_stepsz, val
        type(image_ptr) :: pdiff_opt_odd, pdiff_opt_even, pdiff_odd, pdiff_even, pweights
        ! init
        ldim        = odd%get_ldim()
        box         = ldim(1)
        if( ldim(3) > 1 ) THROW_HARD('This nonuni_filt2D is strictly for 2D!')
        ext         = params_glob%smooth_ext
        find_stop   = optf2Dvars%lplim_hres
        find_start  = calc_fourier_index(params_glob%lpstart_nonuni, box, even%get_smpd())
        find_stepsz = real(find_stop - find_start)/(params_glob%nsearch - 1)
        lb          = (/ ext+1  , ext+1  , 1/)
        ub          = (/ box-ext, box-ext, 1/)
        ! searching for the best fourier index from here and generating the optimized filter
        call weights_img%get_mat_ptrs(pweights)
        call optf2Dvars%diff_img_odd %get_mat_ptrs(pdiff_odd)
        call optf2Dvars%diff_img_even%get_mat_ptrs(pdiff_even)
        call optf2Dvars%diff_img_opt_odd %get_mat_ptrs(pdiff_opt_odd)
        call optf2Dvars%diff_img_opt_even%get_mat_ptrs(pdiff_opt_even)
        pdiff_opt_odd %rmat = huge(val)
        pdiff_opt_even%rmat = huge(val)
        call odd%get_rmat_ptr(rmat_odd)
        call even%get_rmat_ptr(rmat_even)
        call optf2Dvars%odd_filt%get_rmat_ptr(rmat_odd_filt)
        call optf2Dvars%even_filt%get_rmat_ptr(rmat_even_filt)
        rmat_odd  = 0.
        rmat_even = 0.
        do iter_no = 1, params_glob%nsearch
            cutoff_find = nint(find_start + (iter_no - 1)*find_stepsz)
            ! filtering odd/even
            call optf2Dvars%odd_filt%copy_fast(optf2Dvars%odd_copy_cmat)
            call optf2Dvars%even_filt%copy_fast(optf2Dvars%even_copy_cmat)
            call butterworth_filter(optf2Dvars%odd_filt,  cutoff_find, optf2Dvars%cur_fil)
            call butterworth_filter(optf2Dvars%even_filt, cutoff_find, optf2Dvars%cur_fil)
            call optf2Dvars%even_filt%ifft
            call optf2Dvars%odd_filt%ifft
            call optf2Dvars% odd_filt%sqeuclid_matrix(optf2Dvars%even_copy_rmat, optf2Dvars%diff_img_odd)
            call optf2Dvars%even_filt%sqeuclid_matrix(optf2Dvars% odd_copy_rmat, optf2Dvars%diff_img_even)
            ! do the non-uniform, i.e. optimizing at each voxel=
            call optf2Dvars%diff_img_odd%fft
            call optf2Dvars%diff_img_even%fft
            pdiff_odd %cmat  = pdiff_odd %cmat * pweights%cmat
            pdiff_even%cmat  = pdiff_even%cmat * pweights%cmat
            call optf2Dvars%diff_img_odd%ifft
            call optf2Dvars%diff_img_even%ifft
            do l = lb(2),ub(2)
                do k = lb(1),ub(1)
                    if( pdiff_odd%rmat(k,l,1) < pdiff_opt_odd%rmat(k,l,1) )then
                        rmat_odd(          k,l,1) = rmat_odd_filt( k,l,1)
                        pdiff_opt_odd%rmat(k,l,1) = pdiff_odd%rmat(k,l,1)
                    endif
                    if( pdiff_even%rmat(k,l,1) < pdiff_opt_even%rmat(k,l,1) )then
                        rmat_even(          k,l,1) = rmat_even_filt(k,l,1)
                        pdiff_opt_even%rmat(k,l,1) = pdiff_even%rmat(k,l,1)
                    endif
                enddo
            enddo
        enddo
    end subroutine nonuni_filt2D

    subroutine nonuni_filt2D_masked(odd, even, mask, weights_img, optf2Dvars)
        class(image),        intent(inout) :: odd, even, mask, weights_img
        type(optfilt2Dvars), intent(inout) :: optf2Dvars
        real(kind=c_float),        pointer :: rmat_odd(:,:,:), rmat_even(:,:,:), rmat_odd_filt(:,:,:), rmat_even_filt(:,:,:)
        real(kind=c_float),        pointer :: rmat_mask(:,:,:), rmat_odd_lowres(:,:,:), rmat_even_lowres(:,:,:)
        type(image)     :: odd_filt_lowres, even_filt_lowres
        integer         :: k, l, box, ldim(3), find_start, find_stop, iter_no, ext, cutoff_find, lb(3), ub(3)
        real            :: find_stepsz, val, m
        type(image_ptr) :: pdiff_opt_odd, pdiff_opt_even, pdiff_odd, pdiff_even, pweights
        if( .not. optf2Dvars%have_mask )then
            call nonuni_filt2D(odd, even, weights_img, optf2Dvars)
            return
        endif
        ! init
        ldim        = odd%get_ldim()
        box         = ldim(1)
        if( ldim(3) > 1 ) THROW_HARD('This nonuni_filt2D is strictly for 2D!')
        ext         = params_glob%smooth_ext
        find_stop   = optf2Dvars%lplim_hres
        find_start  = calc_fourier_index(params_glob%lpstart_nonuni, box, even%get_smpd())
        find_stepsz = real(find_stop - find_start)/(params_glob%nsearch - 1)
        lb          = (/ ext+1  , ext+1  , 1/)
        ub          = (/ box-ext, box-ext, 1/)
        ! searching for the best fourier index from here and generating the optimized filter
        call weights_img%get_mat_ptrs(pweights)
        call optf2Dvars%diff_img_odd %get_mat_ptrs(pdiff_odd)
        call optf2Dvars%diff_img_even%get_mat_ptrs(pdiff_even)
        call optf2Dvars%diff_img_opt_odd %get_mat_ptrs(pdiff_opt_odd)
        call optf2Dvars%diff_img_opt_even%get_mat_ptrs(pdiff_opt_even)
        pdiff_opt_odd %rmat = huge(val)
        pdiff_opt_even%rmat = huge(val)
        call odd%get_rmat_ptr(rmat_odd)
        call even%get_rmat_ptr(rmat_even)
        call optf2Dvars%odd_filt%get_rmat_ptr(rmat_odd_filt)
        call optf2Dvars%even_filt%get_rmat_ptr(rmat_even_filt)
        rmat_odd  = 0.
        rmat_even = 0.
        ! generate the e/o images with the lowest resolution cutoff frequency applied
        call odd_filt_lowres%copy(optf2Dvars%odd_copy_cmat)
        call even_filt_lowres%copy(optf2Dvars%even_copy_cmat)
        call butterworth_filter(odd_filt_lowres,  find_start, optf2Dvars%cur_fil)
        call butterworth_filter(even_filt_lowres, find_start, optf2Dvars%cur_fil)
        call even_filt_lowres%ifft
        call odd_filt_lowres%ifft
        ! set pointers to e/o lowres and mask
        call odd_filt_lowres%get_rmat_ptr(rmat_odd_lowres)
        call even_filt_lowres%get_rmat_ptr(rmat_even_lowres)
        call mask%get_rmat_ptr(rmat_mask)
        do iter_no = 1, params_glob%nsearch
            cutoff_find = nint(find_start + (iter_no - 1)*find_stepsz)
            ! filtering odd/even
            call optf2Dvars%odd_filt%copy_fast(optf2Dvars%odd_copy_cmat)
            call optf2Dvars%even_filt%copy_fast(optf2Dvars%even_copy_cmat)
            call butterworth_filter(optf2Dvars%odd_filt,  cutoff_find, optf2Dvars%cur_fil)
            call butterworth_filter(optf2Dvars%even_filt, cutoff_find, optf2Dvars%cur_fil)
            call optf2Dvars%even_filt%ifft
            call optf2Dvars%odd_filt%ifft
            call optf2Dvars% odd_filt%sqeuclid_matrix(optf2Dvars%even_copy_rmat, optf2Dvars%diff_img_odd)
            call optf2Dvars%even_filt%sqeuclid_matrix(optf2Dvars% odd_copy_rmat, optf2Dvars%diff_img_even)
            ! do the non-uniform, i.e. optimizing at each voxel
            call optf2Dvars%diff_img_odd%fft
            call optf2Dvars%diff_img_even%fft
            pdiff_odd %cmat  = pdiff_odd %cmat * pweights%cmat
            pdiff_even%cmat  = pdiff_even%cmat * pweights%cmat
            call optf2Dvars%diff_img_odd%ifft
            call optf2Dvars%diff_img_even%ifft
            do l = lb(2),ub(2)
                do k = lb(1),ub(1)
                    if( pdiff_odd%rmat(k,l,1) < pdiff_opt_odd%rmat(k,l,1) )then
                        m = rmat_mask(k,l,1)
                        if( m > 0.99 )then
                            rmat_odd( k,l,1) = rmat_odd_filt( k,l,1)
                        else if( m < 0.01 )then
                            rmat_odd( k,l,1) = rmat_odd_lowres( k,l,1)
                        else
                            rmat_odd( k,l,1) = m * rmat_odd_filt( k,l,1) + (1. - m) * rmat_odd_lowres( k,l,1)
                        endif
                        pdiff_opt_odd%rmat(k,l,1) = pdiff_odd%rmat(k,l,1)
                    endif
                    if( pdiff_even%rmat(k,l,1) < pdiff_opt_even%rmat(k,l,1) )then
                        m = rmat_mask(k,l,1)
                        if( m > 0.99 )then
                            rmat_even(k,l,1) = rmat_even_filt(k,l,1)
                        else if( m < 0.01 )then
                            rmat_even(k,l,1) = rmat_even_lowres(k,l,1)
                        else
                            rmat_even(k,l,1) = m * rmat_even_filt(k,l,1) + (1. - m) * rmat_even_lowres(k,l,1)
                        endif
                        pdiff_opt_even%rmat(k,l,1) = pdiff_even%rmat(k,l,1)
                    endif
                enddo
            enddo
        enddo
        call odd_filt_lowres%kill
        call even_filt_lowres%kill
    end subroutine nonuni_filt2D_masked

    ! 3D optimization(search)-based nonuniform filter, paralellized version
    subroutine nonuni_filt3D(odd, even, mskimg, lpstop)
        class(image),           intent(inout) :: odd, even
        class(image), optional, intent(inout) :: mskimg
        real,         optional, intent(in)    :: lpstop
        type(image)                   ::  odd_copy_rmat, odd_copy_cmat, even_copy_rmat, even_copy_cmat, weights_img,&
                                        &diff_img_opt_odd, diff_img_opt_even, diff_img_odd, diff_img_even, odd_filt, even_filt
        integer                       :: k,l,m, box, ldim(3), find_start, find_stop, iter_no
        integer                       :: filtsz, cutoff_find, lb(3), ub(3), smooth_ext
        real                          :: rad, find_stepsz, val, smpd
        type(image_ptr)               :: pdiff_odd, pdiff_even, pdiff_opt_odd, pdiff_opt_even, pweights
        integer,          parameter   :: CHUNKSZ = 20, N_IMGS = 2
        real,             pointer     :: rmat_odd(:,:,:), rmat_even(:,:,:), rmat_odd_filt(:,:,:), rmat_even_filt(:,:,:)
        real,             allocatable :: fsc(:), cur_fil(:)
        character(len=:), allocatable :: fsc_fname
        real(   kind=c_float),         pointer ::  in(:,:,:,:)
        complex(kind=c_float_complex), pointer :: out(:,:,:,:)
        ldim       = odd%get_ldim()
        filtsz     = odd%get_filtsz()
        smooth_ext = params_glob%smooth_ext
        box        = ldim(1)
        fsc_fname  = trim(params_glob%fsc)
        smpd       = even%get_smpd()
        find_start = calc_fourier_index(params_glob%lpstart_nonuni, box, even%get_smpd())
        if( present(lpstop) )then
            find_stop = calc_fourier_index(lpstop, box, smpd)
        else
            ! calculate Fourier index low-pass limit for search based on FSC
            if( .not.file_exists(fsc_fname) ) then
                THROW_HARD('FSC file: '//fsc_fname//' not found')
            else
                ! retrieve FSC and calculate low-pass limit
                fsc       = file2rarr(fsc_fname)
                find_stop = min(get_lplim_at_corr(fsc, 0.1),calc_fourier_index(params_glob%lpstop, box, smpd)) ! little overshoot, filter function anyway applied in polarft_corrcalc
            endif
        endif
        find_stepsz = real(find_stop - find_start)/(params_glob%nsearch - 1)
        if( find_start >= find_stop ) THROW_HARD('nonuni_filt3D: starting Fourier index is larger than ending Fourier index!')
        allocate(in(ldim(1),ldim(2),ldim(3),2), out(ldim(1),ldim(2),ldim(3),2))
        call       weights_img%new(ldim, smpd)
        call      diff_img_odd%new(ldim, smpd)
        call     diff_img_even%new(ldim, smpd)
        call  diff_img_opt_odd%new(ldim, smpd)
        call diff_img_opt_even%new(ldim, smpd)
        call       weights_img%get_mat_ptrs(pweights)
        call      diff_img_odd%get_mat_ptrs(pdiff_odd)
        call     diff_img_even%get_mat_ptrs(pdiff_even)
        call  diff_img_opt_odd%get_mat_ptrs(pdiff_opt_odd)
        call diff_img_opt_even%get_mat_ptrs(pdiff_opt_even)
        call  odd_copy_rmat%copy(odd)
        call even_copy_rmat%copy(even)
        call  odd_copy_cmat%copy(odd)
        call even_copy_cmat%copy(even)
        call       odd_filt%copy(odd)
        call      even_filt%copy(even)
        call even_copy_cmat%fft
        call odd_copy_cmat%fft
        allocate(cur_fil(box), source=0.)
        if( present(mskimg) )then
            call bounds_from_mask3D(mskimg%bin2logical(), lb, ub)
        else
            lb = (/ 1, 1, 1/)
            ub = (/ box, box, box /)
        endif
        do k = 1, 3
            if( lb(k) < smooth_ext + 1 )   lb(k) = smooth_ext+1
            if( ub(k) > box - smooth_ext ) ub(k) = box - smooth_ext
        enddo
        call weights_img%zero_and_unflag_ft()
        do k = -smooth_ext, smooth_ext
            do l = -smooth_ext, smooth_ext
                do m = -smooth_ext, smooth_ext
                    rad = hyp(k,l,m)
                    val = -rad/(smooth_ext + 1) + 1.
                    if( val > 0 ) call weights_img%set_rmat_at(box/2+k+1, box/2+l+1, box/2+m+1, val)
                enddo
            enddo
        enddo
        call weights_img%fft()
        ! searching for the best fourier index from here and generating the optimized filter
        pdiff_opt_odd%rmat  = huge(val)
        pdiff_opt_even%rmat = huge(val)
        call  odd%get_rmat_ptr(rmat_odd)
        call even%get_rmat_ptr(rmat_even)
        rmat_odd  = 0.
        rmat_even = 0.
        call  odd_filt%get_rmat_ptr( rmat_odd_filt)
        call even_filt%get_rmat_ptr(rmat_even_filt)
        do iter_no = 1, params_glob%nsearch
            cutoff_find = nint(find_start + (iter_no - 1)*find_stepsz)
            ! filtering odd/even
            call  odd_filt%copy_fast( odd_copy_cmat)
            call even_filt%copy_fast(even_copy_cmat)
            call butterworth_filter( odd_filt, cutoff_find, cur_fil)
            call butterworth_filter(even_filt, cutoff_find, cur_fil)
            call even_filt%ifft
            call odd_filt%ifft
            call  odd_filt%sqeuclid_matrix(even_copy_rmat, diff_img_odd)
            call even_filt%sqeuclid_matrix( odd_copy_rmat, diff_img_even)
            ! do the non-uniform, i.e. optimizing at each voxel
            call diff_img_even%fft
            call diff_img_odd%fft
            !$omp parallel workshare
            pdiff_odd% cmat = pdiff_odd %cmat * pweights%cmat
            pdiff_even%cmat = pdiff_even%cmat * pweights%cmat
            !$omp end parallel workshare
            call diff_img_even%ifft
            call diff_img_odd%ifft
            !$omp parallel do collapse(3) default(shared) private(k,l,m) schedule(dynamic,CHUNKSZ) proc_bind(close)
            do m = lb(3),ub(3)
                do l = lb(2),ub(2)
                    do k = lb(1),ub(1)
                        if( pdiff_odd%rmat(k,l,m) < pdiff_opt_odd%rmat(k,l,m) )then
                            rmat_odd(          k,l,m) = rmat_odd_filt( k,l,m)
                            pdiff_opt_odd%rmat(k,l,m) = pdiff_odd%rmat(k,l,m)
                        endif
                        if( pdiff_even%rmat(k,l,m) < pdiff_opt_even%rmat(k,l,m) )then
                            rmat_even(          k,l,m) = rmat_even_filt( k,l,m)
                            pdiff_opt_even%rmat(k,l,m) = pdiff_even%rmat(k,l,m)
                        endif
                    enddo
                enddo
            enddo
            !$omp end parallel do
        enddo
        deallocate(cur_fil)
        call odd_copy_rmat%kill
        call odd_copy_cmat%kill
        call even_copy_rmat%kill
        call even_copy_cmat%kill
        call odd_filt%kill
        call even_filt%kill
        call weights_img%kill
        call diff_img_odd%kill
        call diff_img_even%kill
        call diff_img_opt_odd%kill
        call diff_img_opt_even%kill
    end subroutine nonuni_filt3D

    ! 3D optimization(search)-based uniform filter, paralellized version
    subroutine uni_filt3D(odd, even, mskimg, cutoff_finds_eo)
        class(image),                   intent(inout) :: odd, even
        class(image),         optional, intent(inout) :: mskimg
        integer,              optional, intent(inout) :: cutoff_finds_eo(2)
        real(   kind=c_float),          pointer       ::  in(:,:,:,:)
        complex(kind=c_float_complex),  pointer       :: out(:,:,:,:)
        integer,                        parameter     :: CHUNKSZ = 20, N_IMGS = 2
        real,                           allocatable   :: fsc(:), cur_fil(:)
        character(len=:),               allocatable   :: fsc_fname
        type(image)     ::  odd_copy_rmat, odd_copy_cmat, even_copy_rmat, even_copy_cmat,&
                            &diff_img_odd, diff_img_even, odd_filt, even_filt
        type(image_ptr) :: pdiff_odd, pdiff_even, pweights, pmask
        integer         :: box, ldim(3), find_start, find_stop, iter_no
        integer         :: cutoff_find, best_ind_even, best_ind_odd
        real            :: find_stepsz, smpd, cur_cost_odd, cur_cost_even
        real            :: best_cost_odd, best_cost_even
        logical         :: present_cuofindeo, present_mskimg
        ldim              = odd%get_ldim()
        box               = ldim(1)
        fsc_fname         = trim(params_glob%fsc)
        smpd              = even%get_smpd()
        present_cuofindeo = present(cutoff_finds_eo)
        present_mskimg    = present(mskimg)
        ! calculate Fourier index limits for search
        if( .not.file_exists(fsc_fname) ) then
            THROW_WARN('FSC file: '//fsc_fname//' not found, falling back on lpstart: '//real2str(params_glob%lpstart))
            find_stop = calc_fourier_index(params_glob%lpstart, box, smpd)
        else
            ! retrieve FSC and calculate optimal filter
            fsc       = file2rarr(fsc_fname)
            find_stop = min(get_lplim_at_corr(fsc, 0.1),calc_fourier_index(params_glob%lpstop, box, smpd)) ! little overshoot, filter function anyway applied in polarft_corrcalc
        endif
        find_start  = calc_fourier_index(params_glob%lpstart_nonuni, box, smpd)
        find_stepsz = real(find_stop - find_start)/(params_glob%nsearch - 1)
        if( find_start >= find_stop ) THROW_HARD('uni_filt3D: starting Fourier index is larger than ending Fourier index!')
        allocate(in(ldim(1),ldim(2),ldim(3),2), out(ldim(1),ldim(2),ldim(3),2))
        call   diff_img_odd%new(ldim, smpd)
        call  diff_img_even%new(ldim, smpd)
        call   diff_img_odd%get_mat_ptrs(pdiff_odd)
        call  diff_img_even%get_mat_ptrs(pdiff_even)
        call  odd_copy_rmat%copy(odd)
        call even_copy_rmat%copy(even)
        call  odd_copy_cmat%copy(odd)
        call even_copy_cmat%copy(even)
        call       odd_filt%copy(odd)
        call      even_filt%copy(even)
        call even_copy_cmat%fft
        call odd_copy_cmat%fft
        allocate(cur_fil(box), source=0.)
        if( present_mskimg )then
            call mskimg%get_mat_ptrs(pmask)
        endif
        ! searching for the best fourier index from here and generating the optimized filter
        best_ind_odd   = find_start
        best_ind_even  = find_start
        best_cost_even = huge(best_cost_even)
        best_cost_odd  = huge(best_cost_odd)
        if( present_mskimg )then
            do iter_no = 1, params_glob%nsearch
                cutoff_find = nint(find_start + (iter_no - 1)*find_stepsz)
                ! filtering odd/even
                call  odd_filt%copy_fast( odd_copy_cmat)
                call even_filt%copy_fast(even_copy_cmat)
                call butterworth_filter( odd_filt, cutoff_find, cur_fil)
                call butterworth_filter(even_filt, cutoff_find, cur_fil)
                call even_filt%ifft
                call odd_filt%ifft
                call odd_filt%sqeuclid_matrix(even_copy_rmat, diff_img_odd)
                call even_filt%sqeuclid_matrix( odd_copy_rmat, diff_img_even)
                call diff_img_even%fft
                call diff_img_odd%fft
                !$omp parallel workshare
                pdiff_odd% cmat = pdiff_odd %cmat * pweights%cmat
                pdiff_even%cmat = pdiff_even%cmat * pweights%cmat
                !$omp end parallel workshare
                call diff_img_even%ifft
                call diff_img_odd%ifft
                !$omp parallel workshare
                cur_cost_odd  = sum(pdiff_odd %rmat * pmask%rmat)
                cur_cost_even = sum(pdiff_even%rmat * pmask%rmat)
                !$omp end parallel workshare
                if( cur_cost_odd < best_cost_odd )then
                    best_cost_odd = cur_cost_odd
                    best_ind_odd  = cutoff_find
                endif
                if( cur_cost_even < best_cost_even )then
                    best_cost_even = cur_cost_even
                    best_ind_even  = cutoff_find
                endif
            enddo
        else
            do iter_no = 1, params_glob%nsearch
                cutoff_find = nint(find_start + (iter_no - 1)*find_stepsz)
                ! filtering odd/even
                call  odd_filt%copy_fast( odd_copy_cmat)
                call even_filt%copy_fast(even_copy_cmat)
                call butterworth_filter( odd_filt, cutoff_find, cur_fil)
                call butterworth_filter(even_filt, cutoff_find, cur_fil)
                call even_filt%ifft
                call odd_filt%ifft
                call  odd_filt%sqeuclid_matrix(even_copy_rmat, diff_img_odd)
                call even_filt%sqeuclid_matrix( odd_copy_rmat, diff_img_even)
                call diff_img_even%fft
                call diff_img_odd%fft
                !$omp parallel workshare
                pdiff_odd% cmat = pdiff_odd %cmat * pweights%cmat
                pdiff_even%cmat = pdiff_even%cmat * pweights%cmat
                !$omp end parallel workshare
                call diff_img_even%ifft
                call diff_img_odd%ifft
                !$omp parallel workshare
                cur_cost_odd  = sum(pdiff_odd %rmat(1:ldim(1),1:ldim(2),1:ldim(3)))
                cur_cost_even = sum(pdiff_even%rmat(1:ldim(1),1:ldim(2),1:ldim(3)))
                !$omp end parallel workshare
                if( cur_cost_odd < best_cost_odd )then
                    best_cost_odd = cur_cost_odd
                    best_ind_odd  = cutoff_find
                endif
                if( cur_cost_even < best_cost_even )then
                    best_cost_even = cur_cost_even
                    best_ind_even  = cutoff_find
                endif
            enddo
        endif
        call odd%fft
        call even%fft
        call butterworth_filter(odd,  best_ind_odd,  cur_fil)
        call butterworth_filter(even, best_ind_even, cur_fil)
        call odd%ifft
        call even%ifft
        print *, 'best resolution cut_off = ', calc_lowpass_lim(best_ind_odd, box, smpd)
        if( present_cuofindeo )then
            cutoff_finds_eo(1) = best_ind_odd
            cutoff_finds_eo(2) = best_ind_even
        endif
        deallocate(cur_fil)
        call odd_copy_rmat%kill
        call odd_copy_cmat%kill
        call even_copy_rmat%kill
        call even_copy_cmat%kill
        call odd_filt%kill
        call even_filt%kill
        call diff_img_odd%kill
        call diff_img_even%kill
    end subroutine uni_filt3D

    ! searching for the best index of the cost function |sum(filter(img) - img)|
    ! also return the filtered img at this best index
    subroutine uni_filt2D( odd, even, mask, optf2Dvars )
        class(image),        intent(inout) :: odd, even, mask
        type(optfilt2Dvars), intent(inout) :: optf2Dvars
        real(kind=c_float),        pointer :: rmat_odd(:,:,:), rmat_even(:,:,:), rmat_odd_filt(:,:,:), rmat_even_filt(:,:,:), rmat_mask(:,:,:)
        integer :: box, dim3, ldim(3), find_start, find_stop, iter_no, cutoff_find
        real    :: find_stepsz, best_cost_odd, best_cost_even, cur_cost_odd, cur_cost_even
        ! init
        ldim        = odd%get_ldim()
        box         = ldim(1)
        dim3        = ldim(3)
        if( dim3 > 1 ) THROW_HARD('This nonuni_filt2D is strictly for 2D!')
        find_stop   = optf2Dvars%lplim_hres
        find_start  = calc_fourier_index(params_glob%lpstart_nonuni, box, even%get_smpd())
        find_stepsz = real(find_stop - find_start)/(params_glob%nsearch - 1)
        call odd%get_rmat_ptr(rmat_odd)
        call even%get_rmat_ptr(rmat_even)
        call optf2Dvars%odd_filt%get_rmat_ptr(rmat_odd_filt)
        call optf2Dvars%even_filt%get_rmat_ptr(rmat_even_filt)
        call mask%get_rmat_ptr(rmat_mask)
        best_cost_odd  = huge(best_cost_odd)
        best_cost_even = huge(best_cost_even)
        optf2Dvars%best_ind_odd  = find_start
        optf2Dvars%best_ind_even = find_start
        do iter_no = 1, params_glob%nsearch
            cutoff_find = nint(find_start + (iter_no - 1)*find_stepsz)
            ! filtering odd/even
            call optf2Dvars%odd_filt%copy_fast(optf2Dvars%odd_copy_cmat)
            call optf2Dvars%even_filt%copy_fast(optf2Dvars%even_copy_cmat)
            call butterworth_filter(optf2Dvars%odd_filt,  cutoff_find, optf2Dvars%cur_fil)
            call butterworth_filter(optf2Dvars%even_filt, cutoff_find, optf2Dvars%cur_fil)
            call optf2Dvars%even_filt%ifft
            call optf2Dvars%odd_filt%ifft
            cur_cost_odd  = sum(( rmat_odd_filt - rmat_even)**2 * rmat_mask) ! within the mask
            cur_cost_even = sum((rmat_even_filt - rmat_odd )**2 * rmat_mask) ! within the mask
            if( cur_cost_odd < best_cost_odd )then
                best_cost_odd           = cur_cost_odd
                optf2Dvars%best_ind_odd = cutoff_find
            endif
            if( cur_cost_even < best_cost_even )then
                best_cost_even           = cur_cost_even
                optf2Dvars%best_ind_even = cutoff_find
            endif
        enddo
        call odd%fft
        call even%fft
        call butterworth_filter(odd,  optf2Dvars%best_ind_odd,  optf2Dvars%cur_fil)
        call butterworth_filter(even, optf2Dvars%best_ind_even, optf2Dvars%cur_fil)
        call odd%ifft
        call even%ifft
    end subroutine uni_filt2D

    subroutine uni_filt2D_sub( even, odd, mask )
        use simple_class_frcs, only: class_frcs
        class(image),        intent(inout) :: even(:), odd(:)
        class(image),        intent(inout) :: mask(:)
        character(len=:),    allocatable   :: frcs_fname
        type(optfilt2Dvars), allocatable   :: optf2Dvars(:)
        real,                allocatable   :: frc(:)
        type(class_frcs) :: clsfrcs
        real             :: smpd, lpstart
        integer          :: iptcl, box, filtsz, ldim(3)
        integer          :: nptcls, hpind_fsc, find
        logical          :: lpstart_fallback, l_phaseplate
        write(logfhandle,'(A)') '>>> 2D UNIFORM FILTERING'
        ! init
        ldim         = even(1)%get_ldim()
        filtsz       = even(1)%get_filtsz()
        ldim(3)      = 1 ! because we operate on stacks
        box          = ldim(1)
        frcs_fname   = trim(params_glob%frcs)
        smpd         = even(1)%get_smpd()
        nptcls       = size(even)
        lpstart      = params_glob%lpstart
        hpind_fsc    = params_glob%hpind_fsc
        l_phaseplate = params_glob%l_phaseplate
        ! retrieve FRCs
        lpstart_fallback = .false.
        if( file_exists(frcs_fname) )then
            call clsfrcs%read(frcs_fname)
            if( clsfrcs%get_filtsz().ne.even(1)%get_filtsz() )then
                write(logfhandle,*) 'img filtsz:  ', even(1)%get_filtsz()
                write(logfhandle,*) 'frcs filtsz: ', clsfrcs%get_filtsz()
                THROW_HARD('Inconsistent filter dimensions; nonuni_filt2D_sub')
            endif
        else
            THROW_WARN('Class average FRCs file '//frcs_fname//' does not exist, falling back on lpstart: '//real2str(lpstart))
            lpstart_fallback = .true.
        endif
        filtsz = clsfrcs%get_filtsz()
        ! allocate
        allocate(optf2Dvars(nptcls), frc(filtsz))
        frc = 0.
        ! calculate high-res low-pass limits
        if( lpstart_fallback )then
            optf2Dvars(:)%lplim_hres = calc_fourier_index(lpstart, box, smpd)
        else
            do iptcl = 1, nptcls
                call clsfrcs%frc_getter(iptcl, hpind_fsc, l_phaseplate, frc)
                optf2Dvars(iptcl)%lplim_hres = get_lplim_at_corr(frc, 0.1)
            end do
        endif
        ! fill up optf2Dvars struct
        !$omp parallel do default(shared) private(iptcl) schedule(static) proc_bind(close)
        do iptcl = 1, nptcls
            call optf2Dvars(iptcl)%even_filt%copy(even(iptcl))
            call optf2Dvars(iptcl)%odd_filt %copy( odd(iptcl))
            call optf2Dvars(iptcl)%odd_copy_cmat%copy(odd(iptcl))
            call optf2Dvars(iptcl)%odd_copy_cmat%fft
            call optf2Dvars(iptcl)%even_copy_cmat%copy(even(iptcl))
            call optf2Dvars(iptcl)%even_copy_cmat%fft
            allocate(optf2Dvars(iptcl)%cur_fil(box), source=0.)
        end do
        !$omp end parallel do
        ! filter
        !$omp parallel do default(shared) private(iptcl) schedule(static) proc_bind(close)
        do iptcl = 1, nptcls
            call uni_filt2D(odd(iptcl), even(iptcl), mask(iptcl), optf2Dvars(iptcl))
        enddo
        !$omp end parallel do
        do iptcl = 1, nptcls
            call optf2Dvars(iptcl)%odd_copy_cmat%kill
            call optf2Dvars(iptcl)%even_copy_cmat%kill
            call optf2Dvars(iptcl)%even_filt%kill
            call optf2Dvars(iptcl)%odd_filt%kill
        enddo
        call clsfrcs%kill
    end subroutine uni_filt2D_sub

    ! generating and applying the exponential regularization
    subroutine exponential_reg( vol, lambda, eps )
        class(image),           intent(inout) :: vol
        real,         optional, intent(in)    :: lambda, eps
        type(image_ptr) :: pvol
        integer         :: k, l, m, box, ldim(3)
        real            :: lambda_here, eps_here
        ldim = vol%get_ldim()
        box  = ldim(1)
        call vol%get_mat_ptrs(pvol)
        lambda_here = 1.
        eps_here    = 0.01
        if( present(lambda) ) lambda_here = lambda
        if( present(eps   ) ) eps_here    = eps
        !$omp parallel do collapse(3) default(shared) private(k,l,m) schedule(static) proc_bind(close)
        do m = 1, ldim(3)
            do l = 1, ldim(2)
                do k = 1,ldim(1)
                    pvol%rmat(k,l,m) = pvol%rmat(k,l,m) - eps_here * lambda_here**2 * exp( - lambda_here * pvol%rmat(k,l,m) )
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine exponential_reg

end module simple_opt_filter
