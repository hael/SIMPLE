! optimization(search)-based filter (uniform/nonuniform)
module simple_opt_filter
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_defs
use simple_fftw3
use simple_image,      only: image, image_ptr
use simple_parameters, only: params_glob
implicit none
#include "simple_local_flags.inc"

public :: opt_2D_filter_sub, opt_filter_2D, opt_filter_3D, butterworth_filter
private

type fft_vars_type
    type(c_ptr)                            :: plan_fwd, plan_bwd
    real(   kind=c_float),         pointer ::  in(:,:,:)
    complex(kind=c_float_complex), pointer :: out(:,:,:)
end type fft_vars_type

contains

    subroutine batch_fft_2D( even, odd, fft_vars )
        class(image),        intent(inout) :: even, odd
        type(fft_vars_type), intent(in)    :: fft_vars
        integer            :: ldim(3), k, l
        type(image_ptr)    :: peven, podd
        ldim = even%get_ldim()
        call even%set_ft(.false.)
        call  odd%set_ft(.false.)
        call even%get_mat_ptrs(peven)
        call  odd%get_mat_ptrs(podd)
        !$omp parallel do collapse(2) default(shared) private(k,l) schedule(static) proc_bind(close)
        do l = 1, ldim(2)
            do k = 1, ldim(1)
                fft_vars%in(k,l,1) = peven%rmat(k,l,1)
                fft_vars%in(k,l,2) =  podd%rmat(k,l,1)
            enddo
        enddo
        !$omp end parallel do
        call fftwf_execute_dft_r2c(fft_vars%plan_fwd, fft_vars%in, fft_vars%out)
        call even%set_ft(.true.)
        call  odd%set_ft(.true.)
        !$omp parallel do collapse(2) default(shared) private(k,l) schedule(static) proc_bind(close)
        do l = 1, ldim(2)
            do k = 1, ldim(1)/2+1
                if( mod(k+l,2) == 1 )then
                    peven%cmat(k,l,1) = -fft_vars%out(k,l,1)/product(ldim(1:2))
                    podd %cmat(k,l,1) = -fft_vars%out(k,l,2)/product(ldim(1:2))
                else
                    peven%cmat(k,l,1) =  fft_vars%out(k,l,1)/product(ldim(1:2))
                    podd %cmat(k,l,1) =  fft_vars%out(k,l,2)/product(ldim(1:2))
                endif
            enddo
        enddo
        !$omp end parallel do
    end subroutine batch_fft_2D

    subroutine batch_ifft_2D( even, odd, fft_vars )
        class(image),        intent(inout) :: even, odd
        type(fft_vars_type), intent(in)    :: fft_vars
        integer            :: ldim(3), k, l
        type(image_ptr)    :: peven, podd
        ldim = even%get_ldim()
        call even%set_ft(.true.)
        call  odd%set_ft(.true.)
        call even%get_mat_ptrs(peven)
        call  odd%get_mat_ptrs(podd)
        !$omp parallel do collapse(2) default(shared) private(k,l) schedule(static) proc_bind(close)
        do l = 1, ldim(2)
            do k = 1, ldim(1)/2+1
                if( mod(k+l,2) == 1 )then
                    fft_vars%out(k,l,1) = -peven%cmat(k,l,1)
                    fft_vars%out(k,l,2) = - podd%cmat(k,l,1)
                else
                    fft_vars%out(k,l,1) = peven%cmat(k,l,1)
                    fft_vars%out(k,l,2) =  podd%cmat(k,l,1)
                endif
            enddo
        enddo
        !$omp end parallel do
        call fftwf_execute_dft_c2r(fft_vars%plan_bwd, fft_vars%out, fft_vars%in)
        call even%set_ft(.false.)
        call  odd%set_ft(.false.)
        !$omp parallel do collapse(2) default(shared) private(k,l) schedule(static) proc_bind(close)
        do l = 1, ldim(2)
            do k = 1, ldim(1)
                peven%rmat(k,l,1) = fft_vars%in(k,l,1)
                podd %rmat(k,l,1) = fft_vars%in(k,l,2)
            enddo
        enddo
        !$omp end parallel do
    end subroutine batch_ifft_2D

    subroutine batch_fft_3D( even, odd, in, out, plan_fwd)
        class(image),                           intent(inout) :: even, odd
        real(   kind=c_float),         pointer, intent(inout) ::  in(:,:,:,:)
        complex(kind=c_float_complex), pointer, intent(inout) :: out(:,:,:,:)
        type(c_ptr),                            intent(in)    :: plan_fwd
        integer            :: ldim(3), k, l, m
        type(image_ptr)    :: peven, podd
        ldim = even%get_ldim()
        call even%set_ft(.false.)
        call  odd%set_ft(.false.)
        call even%get_mat_ptrs(peven)
        call  odd%get_mat_ptrs(podd)
        !$omp parallel do collapse(3) default(shared) private(k,l,m) schedule(static) proc_bind(close)
        do m = 1, ldim(3)
            do l = 1, ldim(2)
                do k = 1, ldim(1)
                    in(k,l,m,1) = peven%rmat(k,l,m)
                    in(k,l,m,2) =  podd%rmat(k,l,m)
                enddo
            enddo
        enddo
        !$omp end parallel do
        call fftwf_execute_dft_r2c(plan_fwd, in, out)
        call even%set_ft(.true.)
        call  odd%set_ft(.true.)
        !$omp parallel do collapse(3) default(shared) private(k,l,m) schedule(static) proc_bind(close)
        do m = 1, ldim(3)
            do l = 1, ldim(2)
                do k = 1, ldim(1)/2+1
                    if( mod(k+l+m,2) == 0 )then
                        peven%cmat(k,l,m) = -out(k,l,m,1)/product(ldim)
                        podd %cmat(k,l,m) = -out(k,l,m,2)/product(ldim)
                    else
                        peven%cmat(k,l,m) =  out(k,l,m,1)/product(ldim)
                        podd %cmat(k,l,m) =  out(k,l,m,2)/product(ldim)
                    endif
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine batch_fft_3D

    subroutine batch_ifft_3D( even, odd, in, out, plan_bwd)
        class(image),                           intent(inout) :: even, odd
        real(   kind=c_float),         pointer, intent(inout) ::  in(:,:,:,:)
        complex(kind=c_float_complex), pointer, intent(inout) :: out(:,:,:,:)
        type(c_ptr),                            intent(in)    :: plan_bwd
        integer            :: ldim(3), k, l, m
        type(image_ptr)    :: peven, podd
        ldim = even%get_ldim()
        call even%set_ft(.true.)
        call  odd%set_ft(.true.)
        call even%get_mat_ptrs(peven)
        call  odd%get_mat_ptrs(podd)
        !$omp parallel do collapse(3) default(shared) private(k,l,m) schedule(static) proc_bind(close)
        do m = 1, ldim(3)
            do l = 1, ldim(2)
                do k = 1, ldim(3)/2+1
                    if( mod(k+l+m,2) == 0 )then
                        out(k,l,m,1) = -peven%cmat(k,l,m)
                        out(k,l,m,2) =  -podd%cmat(k,l,m)
                    else
                        out(k,l,m,1) = peven%cmat(k,l,m)
                        out(k,l,m,2) =  podd%cmat(k,l,m)
                    endif
                enddo
            enddo
        enddo
        !$omp end parallel do
        call fftwf_execute_dft_c2r(plan_bwd, out, in)
        call even%set_ft(.false.)
        call  odd%set_ft(.false.)
        !$omp parallel do collapse(3) default(shared) private(k,l,m) schedule(static) proc_bind(close)
        do m = 1, ldim(3)
            do l = 1, ldim(2)
                do k = 1, ldim(1)
                    peven%rmat(k,l,m) = in(k,l,m,1)
                    podd%rmat(k,l,m)  = in(k,l,m,2)
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine batch_ifft_3D

    subroutine opt_2D_filter_sub( even, odd )
        use simple_class_frcs, only: class_frcs
        class(image),      intent(inout) :: even(:), odd(:)
        character(len=:),    allocatable :: frcs_fname
        type(class_frcs)                 :: clsfrcs
        type(image)                      :: weights_img
        type(image),         allocatable :: diff_img(:), diff_img_opt(:), even_filt(:), odd_filt(:),&
                                             &odd_copy_rmat(:),  odd_copy_cmat(:),&
                                            &even_copy_rmat(:), even_copy_cmat(:)
        real,                allocatable :: weights_2D(:,:,:), frc(:), filt(:), cur_fil(:,:)
        integer,             allocatable :: lplims_hres(:)
        real                             :: smpd, lpstart, lp, val
        integer                          :: iptcl, box, filtsz, ldim(3), ldim_pd(3), smooth_ext
        integer                          :: nptcls, hpind_fsc, find, c_shape(3), m, n
        logical                          :: lpstart_fallback, l_phaseplate
        type(c_ptr)                      :: ptr
        integer,             parameter   :: N_IMGS = 2  ! for batch_fft (2 images batch)
        type(fft_vars_type), allocatable :: fft_vars(:)
        ! init
        ldim           = even(1)%get_ldim()
        filtsz         = even(1)%get_filtsz()
        ldim(3)        = 1 ! because we operate on stacks
        smooth_ext     = params_glob%smooth_ext
        ldim_pd        = ldim + 2 * smooth_ext
        ldim_pd(3)     = 1 ! because we operate on stacks
        box            = ldim_pd(1)
        frcs_fname     = trim(params_glob%frcs)
        smpd           = params_glob%smpd
        nptcls         = size(even)
        lpstart        = params_glob%lpstart
        hpind_fsc      = params_glob%hpind_fsc
        l_phaseplate   = params_glob%l_phaseplate
        ! retrieve FRCs
        call clsfrcs%new(nptcls, box, smpd, 1)
        lpstart_fallback = .false.
        if( file_exists(frcs_fname) )then
            call clsfrcs%read(frcs_fname)
            if( clsfrcs%get_filtsz().ne.even(1)%get_filtsz() )then
                write(logfhandle,*) 'img filtsz:  ', even(1)%get_filtsz()
                write(logfhandle,*) 'frcs filtsz: ', clsfrcs%get_filtsz()
                THROW_HARD('Inconsistent filter dimensions; opt_2D_filter_sub')
            endif
        else
            THROW_WARN('Class average FRCs file '//frcs_fname//' does not exist, falling back on lpstart: '//real2str(lpstart))
            lpstart_fallback = .true.
        endif
        filtsz = clsfrcs%get_filtsz()
        ! allocate
        allocate( odd_copy_rmat(nptcls), odd_copy_cmat(nptcls), even_copy_rmat(nptcls), even_copy_cmat(nptcls),&
                      &diff_img(nptcls),  diff_img_opt(nptcls),      even_filt(nptcls),       odd_filt(nptcls))
        allocate(weights_2D(smooth_ext*2+1,smooth_ext*2+1,nptcls), frc(filtsz), cur_fil(box,nptcls), source=0.)
        allocate(lplims_hres(nptcls))
        ! calculate high-res low-pass limits
        if( lpstart_fallback )then
            lplims_hres = calc_fourier_index(lpstart, box, smpd)
        else
            do iptcl = 1, nptcls
                call clsfrcs%frc_getter(iptcl, hpind_fsc, l_phaseplate, frc)
                ! the below required to retrieve the right Fouirer index limit when we are padding
                find = get_lplim_at_corr(frc, 0.1)                     ! little overshoot, filter function anyway applied in polarft_corrcalc
                lp   = calc_lowpass_lim(find, box, smpd)               ! box is the padded box size
                lplims_hres(iptcl) = calc_fourier_index(lp, box, smpd) ! this is the Fourier index limit for the padded images
            end do
        endif
        ! construct: all allocations and initialization
        call weights_img%new(ldim_pd, smpd, .false.)
        call weights_img%zero_and_unflag_ft()
        do m = -params_glob%smooth_ext, params_glob%smooth_ext
            do n = -params_glob%smooth_ext, params_glob%smooth_ext
                val = -hyp(real(m), real(n))/(params_glob%smooth_ext + 1) + 1.
                if( val > 0 ) call weights_img%set_rmat_at(box/2+m+1, box/2+n+1, 1, val)
            enddo
        enddo
        call weights_img%fft()
        do iptcl = 1, nptcls
            call even(iptcl)%pad_mirr(ldim_pd)
            call odd( iptcl)%pad_mirr(ldim_pd)
            call even_filt(iptcl)%copy(even(iptcl))
            call odd_filt( iptcl)%copy(odd( iptcl))
            call diff_img( iptcl)%new(ldim_pd, smpd, .false.)
            call diff_img_opt( iptcl)%new(ldim_pd, smpd, .false.)
            call odd_copy_rmat(iptcl)%copy(odd(iptcl))
            call odd_copy_cmat(iptcl)%copy(odd(iptcl))
            call odd_copy_cmat(iptcl)%fft
            call even_copy_rmat(iptcl)%copy(even(iptcl))
            call even_copy_cmat(iptcl)%copy(even(iptcl))
            call even_copy_cmat(iptcl)%fft
        enddo
        ! batch_fft stuffs
        allocate(fft_vars(nptcls))
        call fftwf_plan_with_nthreads(nthr_glob)
        c_shape = [ldim_pd(1), ldim_pd(2), N_IMGS]
        do iptcl = 1, nptcls
            ptr = fftwf_alloc_complex(int(product(c_shape),c_size_t))
            call c_f_pointer(ptr,fft_vars(iptcl)%out,c_shape)
            call c_f_pointer(ptr,fft_vars(iptcl)%in ,c_shape)
            !$omp critical
            fft_vars(iptcl)%plan_fwd = fftwf_plan_many_dft_r2c(2,    [ldim_pd(2), ldim_pd(1)], N_IMGS,&
                                                &fft_vars(iptcl)%in ,[ldim_pd(2), ldim_pd(1)], 1, product([ldim_pd(2), ldim_pd(1)]),&
                                                &fft_vars(iptcl)%out,[ldim_pd(2), ldim_pd(1)], 1, product([ldim_pd(2), ldim_pd(1)]),FFTW_ESTIMATE)
            fft_vars(iptcl)%plan_bwd = fftwf_plan_many_dft_c2r(2,    [ldim_pd(2), ldim_pd(1)], N_IMGS,&
                                                &fft_vars(iptcl)%out,[ldim_pd(2), ldim_pd(1)], 1, product([ldim_pd(2), ldim_pd(1)]),&
                                                &fft_vars(iptcl)%in ,[ldim_pd(2), ldim_pd(1)], 1, product([ldim_pd(2), ldim_pd(1)]),FFTW_ESTIMATE)
            !$omp end critical
        enddo
        call fftwf_plan_with_nthreads(1)
        ! filter
        !$omp parallel do default(shared) private(iptcl) schedule(static) proc_bind(close)
        do iptcl = 1, nptcls
            call opt_filter_2D(odd(iptcl), even(iptcl), odd_copy_rmat(iptcl), odd_copy_cmat(iptcl),&
            &even_copy_rmat(iptcl), even_copy_cmat(iptcl), odd_filt(iptcl), even_filt(iptcl),&
            &cur_fil(:,iptcl), weights_2D(:,:,iptcl), lplims_hres(iptcl), weights_img, diff_img(iptcl),&
            &diff_img_opt(iptcl), fft_vars(iptcl))
        enddo
        !$omp end parallel do
        ! destruct
        call weights_img%kill
        do iptcl = 1, nptcls
            call odd_copy_rmat( iptcl)%kill
            call even_copy_rmat(iptcl)%kill
            call odd_copy_cmat( iptcl)%kill
            call even_copy_cmat(iptcl)%kill
            call even_filt(     iptcl)%kill
            call odd_filt(      iptcl)%kill
            call diff_img(      iptcl)%kill
            call diff_img_opt(  iptcl)%kill
            call even(iptcl)%clip_inplace(ldim)
            call odd( iptcl)%clip_inplace(ldim)
            call fftwf_destroy_plan(fft_vars(iptcl)%plan_fwd)
            call fftwf_destroy_plan(fft_vars(iptcl)%plan_bwd)
        enddo
    end subroutine opt_2D_filter_sub

    ! Compute the value of the Butterworth transfer function of order n(th)
    ! at a given frequency s, with the cut-off frequency fc
    ! SOURCE :
    ! https://en.wikipedia.org/wiki/Butterworth_filter
    function butterworth(s, n, fc) result(val)
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

    ! Compute the Butterworth kernel of the order n-th of width w
    ! with the cut-off frequency fc
    ! https://en.wikipedia.org/wiki/Butterworth_filter
    subroutine butterworth_filter(ker, n, fc)
        real,    intent(inout) :: ker(:)
        integer, intent(in)    :: n
        real   , intent(in)    :: fc
        integer :: freq_val
        do freq_val = 1, size(ker)
            ker(freq_val) = butterworth(real(freq_val-1), n, fc)
        enddo
    end subroutine butterworth_filter

    subroutine apply_opt_filter(img, cur_ind, find_start, find_stop, cur_fil, use_cache, do_ifft)
        class(image), intent(inout) :: img
        integer,      intent(in)    :: cur_ind
        integer,      intent(in)    :: find_start
        integer,      intent(in)    :: find_stop
        real,         intent(inout) :: cur_fil(:)
        logical,      intent(in)    :: use_cache
        logical,      intent(in), optional :: do_ifft
        integer, parameter :: BW_ORDER = 8
        if( .not. use_cache ) call butterworth_filter(cur_fil, BW_ORDER, real(cur_ind))
        call img%apply_filter(cur_fil)
        if( (.not. present(do_ifft)) .or. do_ifft ) call img%ifft()
    end subroutine apply_opt_filter

    subroutine opt_filter_2D(odd, even, odd_copy_rmat, odd_copy_cmat, even_copy_rmat, even_copy_cmat,&
    &odd_filt, even_filt, cur_fil, weights_2D, kstop, weights_img, diff_img, diff_img_opt, fft_vars)
        class(image),         intent(inout) :: odd, even, odd_filt, even_filt
        class(image),         intent(in)    :: odd_copy_rmat, odd_copy_cmat, even_copy_rmat, even_copy_cmat
        real,                 intent(inout) :: cur_fil(:), weights_2D(:,:)
        integer,              intent(in)    :: kstop
        class(image),         intent(inout) :: weights_img, diff_img, diff_img_opt
        type(fft_vars_type) , intent(in)    :: fft_vars
        real(kind=c_float),   pointer       :: rmat_odd(:,:,:), rmat_even(:,:,:), rmat_odd_filt(:,:,:), rmat_even_filt(:,:,:)
        integer :: k,l, box, dim3, ldim(3), find_start, find_stop, iter_no, ext, best_ind, cur_ind, lb(3), ub(3) 
        real    :: rad, find_stepsz, val
        type(image_ptr) :: pdiff_opt, pdiff, pweights
        ! init
        ldim              = odd%get_ldim()
        box               = ldim(1)
        dim3              = ldim(3)
        if( dim3 > 1 ) THROW_HARD('This opt_filter_2D is strictly for 2D!')
        ext               = params_glob%smooth_ext
        find_stop         = kstop
        find_start        = calc_fourier_index(params_glob%lp_lowres,  box, params_glob%smpd)
        find_stepsz       = real(find_stop - find_start)/(params_glob%nsearch - 1)
        lb                = (/ ext+1  , ext+1  , 1/)
        ub                = (/ box-ext, box-ext, dim3 /)
        ! searching for the best fourier index from here and generating the optimized filter
        best_ind          = find_start
        call  weights_img%get_mat_ptrs(pweights)
        call     diff_img%get_mat_ptrs(pdiff)
        call diff_img_opt%get_mat_ptrs(pdiff_opt)
        pdiff_opt%rmat = huge(val)
        call  odd%get_rmat_ptr(rmat_odd)
        call even%get_rmat_ptr(rmat_even)
        call  odd_filt%get_rmat_ptr(rmat_odd_filt)
        call even_filt%get_rmat_ptr(rmat_even_filt)
        rmat_odd  = 0.
        rmat_even = 0.
        do iter_no = 1, params_glob%nsearch
            cur_ind = nint(find_start + (iter_no - 1)*find_stepsz)
            ! filtering odd/even
            call  odd_filt%copy_fast(odd_copy_cmat)
            call even_filt%copy_fast(even_copy_cmat)
            call apply_opt_filter( odd_filt, cur_ind, find_start, find_stop, cur_fil, .false., do_ifft=.false.)
            call apply_opt_filter(even_filt, cur_ind, find_start, find_stop, cur_fil, .true. , do_ifft=.false.)
            call batch_ifft_2D(even_filt, odd_filt, fft_vars)
            call even_filt%opt_filter_costfun(odd_copy_rmat, odd_filt, even_copy_rmat, diff_img)
            ! do the non-uniform, i.e. optimizing at each voxel
            call diff_img%fft
            pdiff%cmat  = pdiff%cmat * pweights%cmat
            call diff_img%ifft
            do l = lb(2),ub(2)
                do k = lb(1),ub(1)
                    if( pdiff%rmat(k,l,1) < pdiff_opt%rmat(k,l,1) )then
                        rmat_odd( k,l,1)      = rmat_odd_filt( k,l,1)
                        rmat_even(k,l,1)      = rmat_even_filt(k,l,1)
                        pdiff_opt%rmat(k,l,1) = pdiff%rmat(k,l,1)
                    endif
                enddo
            enddo
        enddo
    end subroutine opt_filter_2D

    ! 3D optimization(search)-based uniform/nonuniform filter, paralellized version
    subroutine opt_filter_3D(odd, even, mskimg)
        class(image),           intent(inout) :: odd, even
        class(image), optional, intent(inout) :: mskimg
        type(image)                   ::  odd_copy_rmat, odd_copy_cmat, even_copy_rmat, even_copy_cmat,&
                                        &weights_img, diff_img_opt, diff_img, odd_filt, even_filt
        integer                       :: k,l,m, box, ldim(3), find_start, find_stop, iter_no, fnr
        integer                       :: filtsz, best_ind, cur_ind, lb(3), ub(3), smooth_ext
        real                          :: rad, find_stepsz, val, smpd
        type(image_ptr)               :: pdiff, pdiff_opt, pweights
        type(c_ptr)                   :: plan_fwd, plan_bwd
        integer,          parameter   :: CHUNKSZ = 20, N_IMGS = 2
        real,             pointer     :: rmat_odd(:,:,:), rmat_even(:,:,:), rmat_odd_filt(:,:,:), rmat_even_filt(:,:,:)
        real,             allocatable :: weights_3D(:,:,:), fsc(:), cur_fil(:)
        character(len=:), allocatable :: fsc_fname
        real(   kind=c_float),         pointer ::  in(:,:,:,:)
        complex(kind=c_float_complex), pointer :: out(:,:,:,:)
        ldim         = odd%get_ldim()
        filtsz       = odd%get_filtsz()
        smooth_ext   = params_glob%smooth_ext
        box          = ldim(1)
        fsc_fname    = trim(params_glob%fsc)
        smpd         = params_glob%smpd
        ! retrieve FSC and calculate optimal filter
        if( .not.file_exists(fsc_fname) ) THROW_HARD('FSC file: '//fsc_fname//' not found')
        fsc   = file2rarr(fsc_fname)
        ! calculate Fourier index limits for search
        find_stop   = get_lplim_at_corr(fsc, 0.1) ! little overshoot, filter function anyway applied in polarft_corrcalc
        if( params_glob%lp_stopres > 0 ) find_stop = calc_fourier_index(params_glob%lp_stopres, box, smpd)
        find_start  = calc_fourier_index(params_glob%lp_lowres, box, smpd)
        find_stepsz = real(find_stop - find_start)/(params_glob%nsearch - 1)
        allocate( in(ldim(1),ldim(2),ldim(3),2), out(ldim(1),ldim(2),ldim(3),2))
        !$omp critical
        call fftwf_plan_with_nthreads(nthr_glob)
        plan_fwd = fftwf_plan_many_dft_r2c(3, ldim, N_IMGS, in , ldim, 1, product(ldim), out, ldim, 1, product(ldim), FFTW_ESTIMATE)
        plan_bwd = fftwf_plan_many_dft_c2r(3, ldim, N_IMGS, out, ldim, 1, product(ldim),  in, ldim, 1, product(ldim), FFTW_ESTIMATE)
        call fftwf_plan_with_nthreads(1)
        !$omp end critical
        call   weights_img%new(ldim, smpd)
        call      diff_img%new(ldim, smpd)
        call  diff_img_opt%new(ldim, smpd)
        call   weights_img%get_mat_ptrs(pweights)
        call      diff_img%get_mat_ptrs(pdiff)
        call  diff_img_opt%get_mat_ptrs(pdiff_opt)
        call  odd_copy_rmat%copy(odd)
        call even_copy_rmat%copy(even)
        call  odd_copy_cmat%copy(odd)
        call even_copy_cmat%copy(even)
        call       odd_filt%copy(odd)
        call      even_filt%copy(even)
        call batch_fft_3D(even_copy_cmat, odd_copy_cmat, in, out, plan_fwd)
        allocate(weights_3D(smooth_ext*2+1,smooth_ext*2+1, smooth_ext*2+1), cur_fil(box), source=0.)
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
                    rad = hyp(real(k), real(l), real(m))
                    val = -rad/(smooth_ext + 1) + 1.
                    if( val > 0 ) call weights_img%set_rmat_at(box/2+k+1, box/2+l+1, box/2+m+1, val)
                enddo
            enddo
        enddo
        call weights_img%fft()
        ! searching for the best fourier index from here and generating the optimized filter
        best_ind       = find_start
        pdiff_opt%rmat = huge(val)
        call  odd%get_rmat_ptr(rmat_odd)
        call even%get_rmat_ptr(rmat_even)
        call  odd_filt%get_rmat_ptr( rmat_odd_filt)
        call even_filt%get_rmat_ptr(rmat_even_filt)
        do iter_no = 1, params_glob%nsearch
            cur_ind = nint(find_start + (iter_no - 1)*find_stepsz)
            ! filtering odd/even
            call  odd_filt%copy_fast( odd_copy_cmat)
            call even_filt%copy_fast(even_copy_cmat)
            call apply_opt_filter( odd_filt, cur_ind, find_start, find_stop, cur_fil, .false., .false.)
            call apply_opt_filter(even_filt, cur_ind, find_start, find_stop, cur_fil, .true. , .false.)
            call batch_ifft_3D(even_filt, odd_filt, in, out, plan_bwd)
            call even_filt%opt_filter_costfun_workshare(odd_copy_rmat, odd_filt, even_copy_rmat, diff_img)
            ! do the non-uniform, i.e. optimizing at each voxel
            call diff_img%set_ft(.false.)
            call diff_img%fft
            !$omp parallel workshare
            pdiff%cmat = pdiff%cmat * pweights%cmat
            !$omp end parallel workshare
            call diff_img%ifft
            !$omp parallel do collapse(3) default(shared) private(k,l,m) schedule(dynamic,CHUNKSZ) proc_bind(close)
            do m = lb(3),ub(3)
                do l = lb(2),ub(2)
                    do k = lb(1),ub(1)
                        ! pdiff_opt keeps the minimized cost value at each voxel of the search
                        ! opt_odd  keeps the best voxel of the form B*odd
                        ! opt_even keeps the best voxel of the form B*even
                        if( pdiff%rmat(k,l,m) < pdiff_opt%rmat(k,l,m) )then
                            rmat_odd( k,l,m)      = rmat_odd_filt( k,l,m)
                            rmat_even(k,l,m)      = rmat_even_filt(k,l,m)
                            pdiff_opt%rmat(k,l,m) = pdiff%rmat(k,l,m)
                        endif
                    enddo
                enddo
            enddo
            !$omp end parallel do
        enddo
        deallocate(cur_fil, weights_3D)
        call odd_copy_rmat%kill
        call odd_copy_cmat%kill
        call even_copy_rmat%kill
        call even_copy_cmat%kill
        call odd_filt%kill
        call even_filt%kill
        call weights_img%kill
        call diff_img%kill
        call diff_img_opt%kill
        call fftwf_destroy_plan(plan_fwd)
        call fftwf_destroy_plan(plan_bwd)
    end subroutine opt_filter_3D

end module simple_opt_filter
