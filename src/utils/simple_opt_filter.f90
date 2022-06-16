! optimization(search)-based filter (uniform/nonuniform)
module simple_opt_filter
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_defs
use simple_image,      only: image
use simple_parameters, only: params_glob
implicit none
#include "simple_local_flags.inc"

type opt_vol
    real :: opt_val
    real :: opt_diff
    real :: opt_freq
end type opt_vol

contains

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
            val = sqrt(real(Kn)**2 + aimag(Kn)**2)
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

    subroutine apply_opt_filter(img, cur_ind, find_start, find_stop, cur_fil, use_cache, tvfilt_in)
        use simple_tvfilter, only: tvfilter
        class(image), intent(inout) :: img
        integer,      intent(in)    :: cur_ind
        integer,      intent(in)    :: find_start
        integer,      intent(in)    :: find_stop
        real,         intent(inout) :: cur_fil(:)
        logical,      intent(in)    :: use_cache
        type(tvfilter), optional, intent(inout) :: tvfilt_in
        integer, parameter :: BW_ORDER = 8
        real,    parameter :: LAMBDA_MIN = .5 , LAMBDA_MAX = 5.    ! for TV filter
        real               :: param
        type(tvfilter)     :: tvfilt_loc
        select case(params_glob%filt_enum)
            case(FILT_LP)
                call img%lp(cur_ind)
                call img%ifft()
            case(FILT_TV)
                param = LAMBDA_MIN + (cur_ind - find_start)*(LAMBDA_MAX - LAMBDA_MIN)/(find_stop - find_start)
                if( .not. present(tvfilt_in) )then
                    call tvfilt_loc%new
                    if( img%is_2d() )then
                        call tvfilt_loc%apply_filter(img, param)
                    else
                        call tvfilt_loc%apply_filter_3d(img, param)
                    endif
                    call tvfilt_loc%kill
                else
                    if( img%is_2d() )then
                        call tvfilt_in%apply_filter(img, param)
                    else
                        call tvfilt_in%apply_filter_3d(img, param)
                    endif
                endif
                call img%ifft()
            case(FILT_BW8)
                if( .not. use_cache ) call butterworth_filter(cur_fil, BW_ORDER, real(cur_ind))
                call img%apply_filter(cur_fil)
                call img%ifft()
            case DEFAULT
                THROW_HARD('unsupported filter type')
        end select
    end subroutine apply_opt_filter

    ! 2D optimization(search)-based uniform/nonuniform filter, serial (strictly non-paralellized) version
    subroutine opt_filter_2D(odd, even,&
                            &odd_copy_rmat,  odd_copy_cmat,  odd_copy_shellnorm,&
                            &even_copy_rmat, even_copy_cmat, even_copy_shellnorm,&
                            &tvfilt_in)
        use simple_tvfilter, only: tvfilter
        class(image),   intent(inout) :: odd
        class(image),   intent(inout) :: even
        class(image),   intent(in)    :: odd_copy_rmat,  odd_copy_cmat,  odd_copy_shellnorm,&
                                        &even_copy_rmat, even_copy_cmat, even_copy_shellnorm
        type(tvfilter), intent(inout) :: tvfilt_in
        integer              :: k,l,m,n, box, dim3, ldim(3), find_start, find_stop, iter_no, ext
        integer              :: best_ind, cur_ind, k1, l1, m1, lb(3), ub(3), mid_ext
        real                 :: min_sum_odd, min_sum_even, ref_diff_odd, ref_diff_even, rad, find_stepsz, val
        character(len=90)    :: file_tag
        real,    pointer     :: rmat_odd(:,:,:)=>null(), rmat_even(:,:,:)=>null()
        real,    allocatable :: cur_diff_odd( :,:,:), cur_diff_even(:,:,:)
        real,    allocatable :: cur_fil(:), weights_2D(:,:)
        type(opt_vol), allocatable :: opt_odd(:,:,:), opt_even(:,:,:)
        !$omp critical
        ldim        = odd%get_ldim()
        box         = ldim(1)
        dim3        = ldim(3)
        if( dim3 > 1 ) THROW_HARD('This opt_filter_2D is strictly for 2D case only!')
        find_stop   = calc_fourier_index(2. * params_glob%smpd , box, params_glob%smpd)
        find_start  = calc_fourier_index(     params_glob%lp_lb, box, params_glob%smpd)
        find_stepsz = real(find_stop - find_start)/(params_glob%nsearch - 1)
        lb          = (/ params_glob%smooth_ext+1  , params_glob%smooth_ext+1  , 1/)
        ub          = (/ box-params_glob%smooth_ext, box-params_glob%smooth_ext, dim3 /)
        allocate(cur_diff_odd( box,box,dim3), cur_diff_even(box,box,dim3),&
                &cur_fil(box),weights_2D(params_glob%smooth_ext*2+1, params_glob%smooth_ext*2+1), source=0.)
        allocate(opt_odd(box,box,dim3), opt_even(box,box,dim3))
        ! searching for the best fourier index from here and generating the optimized filter
        min_sum_odd       = huge(min_sum_odd)
        min_sum_even      = huge(min_sum_even)
        best_ind          = find_start
        opt_odd%opt_val   = 0.
        opt_odd%opt_diff  = 0.
        opt_odd%opt_freq  = 0.
        opt_even%opt_val  = 0.
        opt_even%opt_diff = 0.
        opt_even%opt_freq = 0.
        opt_odd( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3))%opt_diff = huge(min_sum_odd)
        opt_even(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3))%opt_diff = huge(min_sum_odd)
        do iter_no = 1, params_glob%nsearch
            cur_ind = find_start + (iter_no - 1)*find_stepsz
            if( L_VERBOSE_GLOB ) write(*,*) '('//int2str(iter_no)//'/'//int2str(params_glob%nsearch)//') current Fourier index = ', cur_ind
            ! filtering odd
            call odd%copy_fast(odd_copy_cmat)
            call apply_opt_filter(odd, cur_ind, find_start, find_stop, cur_fil, .false., tvfilt_in)
            call odd%sqeuclid_matrix(even_copy_rmat, cur_diff_odd)
            if( params_glob%l_match_filt )then
                call odd%copy_fast(odd_copy_shellnorm)
                call apply_opt_filter(odd, cur_ind, find_start, find_stop, cur_fil, .false., tvfilt_in)
            endif
            call odd%get_rmat_ptr(rmat_odd)
            ! filtering even
            call even%copy_fast(even_copy_cmat)
            call apply_opt_filter(even, cur_ind, find_start, find_stop, cur_fil, .true., tvfilt_in)
            call even%sqeuclid_matrix(odd_copy_rmat, cur_diff_even)
            if( params_glob%l_match_filt )then
                call even%copy_fast(even_copy_shellnorm)
                call apply_opt_filter(even, cur_ind, find_start, find_stop, cur_fil, .false., tvfilt_in)
            endif
            call even%get_rmat_ptr(rmat_even)
            ! do the non-uniform, i.e. optimizing at each voxel
            if( params_glob%l_nonuniform )then
                ! searching through the smoothing extension here
                do ext = params_glob%smooth_ext,params_glob%smooth_ext
                    ! setting up the 2D weights
                    weights_2D = 0.
                    mid_ext    = 1 + ext
                    do m = 1, 2*ext+1
                        do n = 1, 2*ext+1
                            rad = hyp(real(m-mid_ext), real(n-mid_ext))
                            weights_2D(m,n) = -rad/(ext + 1) + 1.  ! linear function: 1 at rad = 0 and 0 at rad = smooth_ext + 1
                            if (weights_2D(m,n) < 0.) then
                                weights_2D(m,n) = 0.
                            endif
                        enddo
                    enddo
                    weights_2D = weights_2D/sum(weights_2D) ! weights has energy of 1
                    do l = lb(2),ub(2)
                        do k = lb(1),ub(1)
                            ! applying the smoothing extension to the difference
                            k1 = k - ext
                            l1 = l - ext
                            ref_diff_odd  = sum(sum(cur_diff_odd(k1:k1+2*ext,&
                                                                &l1:l1+2*ext,1)*weights_2D(1:2*ext+1, 1:2*ext+1),&
                                                    &dim=2), dim=1)
                            ref_diff_even = sum(sum(cur_diff_even(k1:k1+2*ext,&
                                                                 &l1:l1+2*ext,1)*weights_2D(1:2*ext+1, 1:2*ext+1),&
                                                    &dim=2), dim=1)
                            ! opt_diff keeps the minimized cost value at each voxel of the search
                            ! opt_odd  keeps the best voxel of the form B*odd
                            ! opt_even keeps the best voxel of the form B*even
                            if (ref_diff_odd < opt_odd(k,l,1)%opt_diff) then
                                opt_odd(k,l,1)%opt_val  = rmat_odd(k,l,1)
                                opt_odd(k,l,1)%opt_diff = ref_diff_odd
                                opt_odd(k,l,1)%opt_freq = cur_ind
                            endif
                            if (ref_diff_even < opt_even(k,l,1)%opt_diff) then
                                opt_even(k,l,1)%opt_val  = rmat_even(k,l,1)
                                opt_even(k,l,1)%opt_diff = ref_diff_even
                                opt_even(k,l,1)%opt_freq = cur_ind
                            endif
                        enddo
                    enddo
                enddo
            else
                ! keep the theta which gives the lowest cost (over all voxels)
                if (sum(cur_diff_odd) < min_sum_odd) then
                    opt_odd(:,:,:)%opt_val  = rmat_odd(1:box, 1:box, 1:dim3)
                    opt_odd(:,:,:)%opt_freq = cur_ind
                    min_sum_odd  = sum(cur_diff_odd)
                    best_ind     = cur_ind
                endif
                if (sum(cur_diff_even) < min_sum_even) then
                    opt_even(:,:,:)%opt_val  = rmat_even(1:box, 1:box, 1:dim3)
                    opt_even(:,:,:)%opt_freq = cur_ind
                    min_sum_even  = sum(cur_diff_even)
                endif
            endif
            if( L_VERBOSE_GLOB ) write(*,*) 'current cost (odd) = ', sum(cur_diff_odd)
        enddo
        if( L_VERBOSE_GLOB )then
            if( .not. params_glob%l_nonuniform ) write(*,*) 'minimized cost at resolution = ', box*params_glob%smpd/best_ind
        endif
        call odd%set_rmat( opt_odd( :,:,:)%opt_val, .false.)
        call even%set_rmat(opt_even(:,:,:)%opt_val, .false.)
        deallocate(opt_odd, opt_even, cur_diff_odd, cur_diff_even, cur_fil, weights_2D)
        !$omp end critical
    end subroutine opt_filter_2D

    ! 3D optimization(search)-based uniform/nonuniform filter, paralellized version
    subroutine opt_filter_3D(odd, even, mskimg)
        class(image),           intent(inout) :: odd
        class(image),           intent(inout) :: even
        class(image), optional, intent(inout) :: mskimg
        type(image)          ::  odd_copy_rmat,  odd_copy_cmat,  odd_copy_shellnorm, freq_img,&
                               &even_copy_rmat, even_copy_cmat, even_copy_shellnorm
        integer              :: k,l,m, box, ldim(3), find_start, find_stop, iter_no, fnr
        integer              :: best_ind, cur_ind, k1, l1, m1, lb(3), ub(3), mid_ext
        real                 :: min_sum_odd, min_sum_even, ref_diff_odd, ref_diff_even, rad, find_stepsz
        logical              :: mskimg_present
        character(len=90)    :: file_tag
        integer, parameter   :: CHUNKSZ = 20
        real,    pointer     :: rmat_odd(:,:,:)=>null(), rmat_even(:,:,:)=>null()
        real,    allocatable :: cur_diff_odd( :,:,:), cur_diff_even(:,:,:)
        real,    allocatable :: cur_fil(:), weights_3D(:,:,:)
        logical, allocatable :: l_mask(:,:,:)
        type(opt_vol), allocatable :: opt_odd(:,:,:), opt_even(:,:,:)
        character(len=LONGSTRLEN)     :: benchfname
        integer(timer_int_kind)       ::  t_tot,  t_filter_odd,  t_filter_even,  t_search_opt, &
                                        & t_chop_copy,  t_chop_filter,  t_chop_sqeu
        real(timer_int_kind)          :: rt_tot, rt_filter_odd, rt_filter_even, rt_search_opt, &
                                        &rt_chop_copy, rt_chop_filter, rt_chop_sqeu
        mskimg_present   = present(mskimg)
        if( mskimg_present ) l_mask = mskimg%bin2logical()
        ldim        = odd%get_ldim()
        box         = ldim(1)
        mid_ext     = 1 + params_glob%smooth_ext
        find_stop   = calc_fourier_index(2. * params_glob%smpd , box, params_glob%smpd)
        find_start  = calc_fourier_index(     params_glob%lp_lb, box, params_glob%smpd)
        find_stepsz = real(find_stop - find_start)/(params_glob%nsearch - 1)
        call freq_img%new([box,box,box], params_glob%smpd)
        call odd_copy_rmat%copy(odd)
        call odd_copy_cmat%copy(odd)
        call odd_copy_cmat%fft
        call odd_copy_shellnorm%copy(odd)
        call odd_copy_shellnorm%shellnorm(return_ft=.true.)
        call even_copy_rmat%copy(even)
        call even_copy_cmat%copy(even)
        call even_copy_cmat%fft
        call even_copy_shellnorm%copy(even)
        call even_copy_shellnorm%shellnorm(return_ft=.true.)
        allocate(cur_diff_odd( box,box,box), cur_diff_even(box,box,box),&
                &cur_fil(box), weights_3D(params_glob%smooth_ext*2+1,&
                &params_glob%smooth_ext*2+1, params_glob%smooth_ext*2+1), source=0.)
        allocate(opt_odd(box,box,box), opt_even(box,box,box))
        ! assign the weights of the neighboring voxels
        ! 3D weights
        do k = 1, 2*params_glob%smooth_ext+1
            do l = 1, 2*params_glob%smooth_ext+1
                do m = 1, 2*params_glob%smooth_ext+1
                    rad = hyp(real(k-mid_ext), real(l-mid_ext), real(m-mid_ext))
                    weights_3D(k,l,m) = -rad/(params_glob%smooth_ext + 1) + 1.  ! linear function: 1 at rad = 0 and 0 at rad = smooth_ext + 1
                    if (weights_3D(k,l,m) < 0.) then
                        weights_3D(k,l,m) = 0.
                    endif
                enddo
            enddo
        enddo
        weights_3D = weights_3D/sum(weights_3D) ! weights has energy of 1
        if( mskimg_present )then
            call bounds_from_mask3D(l_mask, lb, ub)
        else
            lb = (/ 1, 1, 1/)
            ub = (/ box, box, box /)
        endif
        do k = 1, 3
            if( lb(k) < params_glob%smooth_ext + 1 )   lb(k) = params_glob%smooth_ext+1
            if( ub(k) > box - params_glob%smooth_ext ) ub(k) = box - params_glob%smooth_ext
        enddo
        ! searching for the best fourier index from here and generating the optimized filter
        min_sum_odd       = huge(min_sum_odd)
        min_sum_even      = huge(min_sum_even)
        best_ind          = find_start
        opt_odd%opt_val   = 0.
        opt_odd%opt_diff  = 0.
        opt_odd%opt_freq  = 0.
        opt_even%opt_val  = 0.
        opt_even%opt_diff = 0.
        opt_even%opt_freq = 0.
        opt_odd( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3))%opt_diff = huge(min_sum_odd)
        opt_even(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3))%opt_diff = huge(min_sum_odd)
        if( L_BENCH_GLOB )then
            t_tot          = tic()
            rt_filter_odd  = 0.
            rt_filter_even = 0.
            rt_search_opt  = 0.
            rt_chop_copy   = 0.
            rt_tot         = 0.
            rt_chop_filter = 0.
            rt_chop_sqeu   = 0.
        endif
        do iter_no = 1, params_glob%nsearch
            cur_ind = find_start + (iter_no - 1)*find_stepsz
            if( L_VERBOSE_GLOB ) write(*,*) '('//int2str(iter_no)//'/'//int2str(params_glob%nsearch)//') current Fourier index = ', cur_ind
            if( L_BENCH_GLOB )then
                t_filter_odd = tic()
                t_chop_copy  = tic()
            endif
            ! filtering odd
            call odd%copy_fast(odd_copy_cmat)
            if( L_BENCH_GLOB )then
                rt_chop_copy  = rt_chop_copy + toc(t_chop_copy)
                t_chop_filter = tic()
            endif
            call apply_opt_filter(odd, cur_ind, find_start, find_stop, cur_fil, .false.)
            if( L_BENCH_GLOB )then
                rt_chop_filter = rt_chop_filter + toc(t_chop_filter)
                t_chop_sqeu    = tic()
            endif
            call odd%sqeuclid_matrix(even_copy_rmat, cur_diff_odd)
            if( L_BENCH_GLOB )then
                rt_chop_sqeu = rt_chop_sqeu + toc(t_chop_sqeu)
            endif
            if( params_glob%l_match_filt )then
                call odd%copy_fast(odd_copy_shellnorm)
                if( L_BENCH_GLOB ) t_chop_filter = tic()
                call apply_opt_filter(odd, cur_ind, find_start, find_stop, cur_fil, .false.)
                if( L_BENCH_GLOB ) rt_chop_filter = rt_chop_filter + toc(t_chop_filter)
            endif
            call odd%get_rmat_ptr(rmat_odd)
            if( L_BENCH_GLOB )then
                rt_filter_odd = rt_filter_odd + toc(t_filter_odd)
                t_filter_even = tic()
            endif
            ! filtering even
            call even%copy_fast(even_copy_cmat)
            call apply_opt_filter(even, cur_ind, find_start, find_stop, cur_fil, .true.)
            call even%sqeuclid_matrix(odd_copy_rmat, cur_diff_even)
            if( params_glob%l_match_filt )then
                call even%copy_fast(even_copy_shellnorm)
                call apply_opt_filter(even, cur_ind, find_start, find_stop, cur_fil, .false.)
            endif
            call even%get_rmat_ptr(rmat_even)
            if( L_BENCH_GLOB )then
                rt_filter_even = rt_filter_even + toc(t_filter_even)
                t_search_opt   = tic()
            endif
            ! do the non-uniform, i.e. optimizing at each voxel
            if( params_glob%l_nonuniform )then
                !$omp parallel do collapse(3) default(shared) private(k,l,m,k1,l1,m1,ref_diff_odd,ref_diff_even) schedule(dynamic,CHUNKSZ) proc_bind(close)
                do m = lb(3),ub(3)
                    do l = lb(2),ub(2)
                        do k = lb(1),ub(1)
                            ! applying an average window to each diff (eq 7 in the nonuniform paper)
                            k1 = k - params_glob%smooth_ext
                            l1 = l - params_glob%smooth_ext
                            m1 = m - params_glob%smooth_ext
                            ref_diff_odd  = sum(sum(sum(cur_diff_odd( k1:k1+2*params_glob%smooth_ext,&
                                                                        &l1:l1+2*params_glob%smooth_ext,&
                                                                        &m1:m1+2*params_glob%smooth_ext)*weights_3D,&
                                                        &dim=3), dim=2), dim=1)
                            ref_diff_even = sum(sum(sum(cur_diff_even(k1:k1+2*params_glob%smooth_ext,&
                                                                        &l1:l1+2*params_glob%smooth_ext,&
                                                                        &m1:m1+2*params_glob%smooth_ext)*weights_3D,&
                                                        &dim=3), dim=2), dim=1)
                            ! opt_diff keeps the minimized cost value at each voxel of the search
                            ! opt_odd  keeps the best voxel of the form B*odd
                            ! opt_even keeps the best voxel of the form B*even
                            if (ref_diff_odd < opt_odd(k,l,m)%opt_diff) then
                                opt_odd(k,l,m)%opt_val  = rmat_odd(k,l,m)
                                opt_odd(k,l,m)%opt_diff = ref_diff_odd
                                opt_odd(k,l,m)%opt_freq = cur_ind
                            endif
                            if (ref_diff_even < opt_even(k,l,m)%opt_diff) then
                                opt_even(k,l,m)%opt_val  = rmat_even(k,l,m)
                                opt_even(k,l,m)%opt_diff = ref_diff_even
                                opt_even(k,l,m)%opt_freq = cur_ind
                            endif
                        enddo
                    enddo
                enddo
                !$omp end parallel do
            else
                ! keep the theta which gives the lowest cost (over all voxels)
                if (sum(cur_diff_odd) < min_sum_odd) then
                    opt_odd(:,:,:)%opt_val  = rmat_odd(1:box, 1:box, 1:box)
                    opt_odd(:,:,:)%opt_freq = cur_ind
                    min_sum_odd  = sum(cur_diff_odd)
                    best_ind     = cur_ind
                endif
                if (sum(cur_diff_even) < min_sum_even) then
                    opt_even(:,:,:)%opt_val  = rmat_even(1:box, 1:box, 1:box)
                    opt_even(:,:,:)%opt_freq = cur_ind
                    min_sum_even  = sum(cur_diff_even)
                endif
            endif
            if( L_VERBOSE_GLOB ) write(*,*) 'current cost (odd) = ', sum(cur_diff_odd)
            if( L_BENCH_GLOB )then
                rt_search_opt = rt_search_opt + toc(t_search_opt)
            endif
        enddo
        if( L_BENCH_GLOB )then
            rt_tot     = toc(t_tot)
            benchfname = 'OPT_FILTER_BENCH.txt'
            call fopen(fnr, FILE=trim(benchfname), STATUS='REPLACE', action='WRITE')
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f9.2)') 'copy_fast            : ', rt_chop_copy
            write(fnr,'(a,1x,f9.2)') 'lp_filter and ifft   : ', rt_chop_filter
            write(fnr,'(a,1x,f9.2)') 'sqeuclid_matrix      : ', rt_chop_sqeu
            write(fnr,'(a,1x,f9.2)') 'odd filtering        : ', rt_filter_odd
            write(fnr,'(a,1x,f9.2)') 'even filtering       : ', rt_filter_even
            write(fnr,'(a,1x,f9.2)') 'searching/optimizing : ', rt_search_opt
            write(fnr,'(a,1x,f9.2)') 'total time           : ', rt_tot
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
            write(fnr,'(a,1x,f9.2)') 'odd filtering        : ', (rt_filter_odd /rt_tot) * 100. 
            write(fnr,'(a,1x,f9.2)') 'even filtering       : ', (rt_filter_even/rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') 'searching/optimizing : ', (rt_search_opt /rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') '% accounted for      : ',&
            &((rt_filter_odd+rt_filter_even+rt_search_opt)/rt_tot) * 100.
            call fclose(fnr)
        endif
        if( L_VERBOSE_GLOB )then
            if( .not. params_glob%l_nonuniform ) write(*,*) 'minimized cost at resolution = ', box*params_glob%smpd/best_ind
        endif
        call odd%set_rmat( opt_odd( :,:,:)%opt_val, .false.)
        call even%set_rmat(opt_even(:,:,:)%opt_val, .false.)
        ! output the optimized frequency map to see the nonuniform parts
        if( params_glob%l_nonuniform )then
            file_tag = 'nonuniform_filter_'//trim(params_glob%filter)//'_ext_'//int2str(params_glob%smooth_ext)
        else
            file_tag = 'uniform_filter_'//trim(params_glob%filter)//'_ext_'//int2str(params_glob%smooth_ext)
        endif
        call freq_img%set_rmat(opt_odd(:,:,:)%opt_freq, .false.)
        call freq_img%write('opt_freq_odd_map_'//trim(file_tag)//'.mrc')
        call freq_img%set_rmat(box*params_glob%smpd/opt_odd(:,:,:)%opt_freq, .false.) ! resolution map
        call freq_img%write('opt_resolution_odd_map_'//trim(file_tag)//'.mrc')
        call freq_img%kill
        deallocate(opt_odd, opt_even, cur_diff_odd, cur_diff_even, cur_fil, weights_3D)
        call odd_copy_rmat%kill
        call odd_copy_cmat%kill
        call even_copy_rmat%kill
        call even_copy_cmat%kill
    end subroutine opt_filter_3D

end module simple_opt_filter
    