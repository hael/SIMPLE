! optimization(search)-based filter (uniform/nonuniform)
module simple_opt_filter
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_defs
use simple_image,      only: image
use simple_parameters, only: params_glob
implicit none

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

    subroutine apply_opt_filter(img, param, cur_fil, use_cache)
        use simple_tvfilter, only: tvfilter
        class(image),     intent(inout) :: img
        real,             intent(in)    :: param
        real,             intent(inout) :: cur_fil(:)
        logical,          intent(in)    :: use_cache
        integer                         :: bw_order, io_stat
        type(tvfilter)                  :: tvfilt
        select case(trim(params_glob%filter))
            case('lp')
                call img%lp(int(param))
            case('tv')
                call tvfilt%new
                call tvfilt%apply_filter_3d(img, param)
                call tvfilt%kill
            case DEFAULT ! default to butterworth8, even if wrong filter type is entered
                bw_order = 8
                ! extract butterworth order number
                if( str_has_substr(params_glob%filter, 'butterworth') )then
                    call str2int(params_glob%filter(12:len_trim(params_glob%filter)), io_stat, bw_order)
                    if (bw_order < 1 .or. bw_order > 10) bw_order = 8
                endif
                if( .not. use_cache ) call butterworth_filter(cur_fil, bw_order, param)
                call img%apply_filter(cur_fil)
        end select
        call img%ifft()
    end subroutine apply_opt_filter

    ! optimization(search)-based uniform/nonuniform filter, using the (low-pass/butterworth)
    subroutine opt_filter(odd, even, mskimg)
        class(image),           intent(inout) :: odd
        class(image),           intent(inout) :: even
        class(image), optional, intent(inout) :: mskimg
        type(image)          :: odd_copy_rmat, odd_copy_cmat, even_copy_rmat, even_copy_cmat, freq_img
        integer              :: k,l,m, box, dim3, ldim(3), find_start, find_stop, iter_no, fnr
        integer              :: best_ind, cur_ind, k1,l1,m1,k_ind,l_ind,m_ind, lb(3), ub(3), mid_ext
        real                 :: min_sum_odd, min_sum_even, ref_diff_odd, ref_diff_even, rad, param, find_stepsz
        logical              :: mskimg_present
        character(len=90)    :: file_tag
        integer, parameter   :: CHUNKSZ = 20
        real,    parameter   :: LAMBDA_MIN = .5, LAMBDA_MAX = 5.
        real,    pointer     :: rmat_odd(:,:,:)=>null(), rmat_even(:,:,:)=>null()
        real,    allocatable :: opt_odd( :,:,:), cur_diff_odd( :,:,:), opt_diff_odd( :,:,:), opt_freq_odd(:,:,:)
        real,    allocatable :: opt_even(:,:,:), cur_diff_even(:,:,:), opt_diff_even(:,:,:), opt_freq_even(:,:,:)
        real,    allocatable :: cur_fil(:), weights_3D(:,:,:), weights_2D(:,:)
        logical, allocatable :: l_mask(:,:,:)
        character(len=LONGSTRLEN)     :: benchfname
        integer(timer_int_kind)       ::  t_tot,  t_filter_odd,  t_filter_even,  t_search_opt
        real(timer_int_kind)          :: rt_tot, rt_filter_odd, rt_filter_even, rt_search_opt
        mskimg_present   = present(mskimg)
        if( mskimg_present ) l_mask = mskimg%bin2logical()
        ldim        = odd%get_ldim()
        box         = ldim(1)
        dim3        = ldim(3)
        mid_ext     = 1 + params_glob%smooth_ext
        find_stop   = calc_fourier_index(2. * params_glob%smpd, box, params_glob%smpd)
        find_start  = calc_fourier_index(params_glob%lp_lb, box, params_glob%smpd)
        find_stepsz = real(find_stop - find_start)/(params_glob%nsearch - 1)
        call freq_img%new([box,box,dim3], params_glob%smpd)
        call odd_copy_rmat%copy(odd)
        call odd_copy_cmat%copy(odd)
        call odd_copy_cmat%fft
        call even_copy_rmat%copy(even)
        call even_copy_cmat%copy(even)
        call even_copy_cmat%fft
        allocate(opt_odd( box,box,dim3), cur_diff_odd( box,box,dim3), opt_diff_odd( box,box,dim3), opt_freq_odd( box,box,dim3),&
                &opt_even(box,box,dim3), cur_diff_even(box,box,dim3), opt_diff_even(box,box,dim3), opt_freq_even(box,box,dim3),&
                &cur_fil(box),weights_2D(params_glob%smooth_ext*2+1, params_glob%smooth_ext*2+1), weights_3D(params_glob%smooth_ext*2+1, params_glob%smooth_ext*2+1, params_glob%smooth_ext*2+1), source=0.)
        ! assign the weights of the neighboring voxels
        ! 2D weights
        do k = 1, 2*params_glob%smooth_ext+1
            do l = 1, 2*params_glob%smooth_ext+1
                rad = hyp(real(k-mid_ext), real(l-mid_ext))
                weights_2D(k,l) = -rad/(params_glob%smooth_ext + 1) + 1.  ! linear function: 1 at rad = 0 and 0 at rad = smooth_ext + 1
                if (weights_2D(k,l) < 0.) then
                    weights_2D(k,l) = 0.
                endif
            enddo
        enddo
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
        weights_2D = weights_2D/sum(weights_2D) ! weights has energy of 1
        ! determine loop bounds for better load balancing in the following parallel loop
        if( mskimg_present )then
            call bounds_from_mask3D(l_mask, lb, ub)
        else
            lb = (/ 1, 1, 1/)
            ub = (/ box, box, box /)
        endif
        ! searching for the best fourier index from here and generating the optimized filter
        opt_diff_odd  = 0.
        opt_diff_even = 0.
        opt_freq_odd  = 0.   ! record the optimized cutoff frequency
        opt_freq_even = 0.   ! record the optimized cutoff frequency
        min_sum_odd   = huge(min_sum_odd)
        min_sum_even  = huge(min_sum_even)
        best_ind      = find_start
        opt_diff_odd( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) = huge(min_sum_odd)
        opt_diff_even(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) = huge(min_sum_odd)
        opt_freq_odd( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) = huge(min_sum_odd)
        opt_freq_even(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) = huge(min_sum_odd)
        if( L_BENCH_GLOB )then
            t_tot          = tic()
            rt_filter_odd  = 0.
            rt_filter_even = 0.
            rt_search_opt  = 0.
        endif
        do iter_no = 1, params_glob%nsearch
            cur_ind = find_start + (iter_no - 1)*find_stepsz
            if( params_glob%filter == 'tv' )then
                param = LAMBDA_MIN + (cur_ind - find_start)*(LAMBDA_MAX - LAMBDA_MIN)/(find_stop - find_start)
                write(*, *) '('//int2str(iter_no)//'/'//int2str(params_glob%nsearch)//') current lambda = ', param
            else
                param = real(cur_ind)
                write(*, *) '('//int2str(iter_no)//'/'//int2str(params_glob%nsearch)//') current Fourier index = ', param
            endif
            if( L_BENCH_GLOB ) t_filter_odd = tic()
            ! filtering odd
            call odd%copy_fast(odd_copy_cmat)
            call apply_opt_filter(odd, param, cur_fil, .false.)
            call odd%get_rmat_ptr(rmat_odd)
            call even_copy_rmat%get_rmat_ptr(rmat_even)
            call odd%sqeuclid_matrix(even_copy_rmat, cur_diff_odd)
            if( params_glob%l_match_filt )then
                call odd%set_ft(.true.)
                call odd%copy_fast(odd_copy_cmat)
                call odd%shellnorm
                call apply_opt_filter(odd, param, cur_fil, .false.)
            endif
            if( L_BENCH_GLOB )then
                rt_filter_odd = rt_filter_odd + toc(t_filter_odd)
                t_filter_even = tic()
            endif
            ! filtering even
            call even%copy_fast(even_copy_cmat)
            call apply_opt_filter(even, param, cur_fil, .true.)
            call even%get_rmat_ptr(rmat_even)
            call odd_copy_rmat%get_rmat_ptr(rmat_odd)
            call even%sqeuclid_matrix(odd_copy_rmat, cur_diff_even)
            if( params_glob%l_match_filt )then
                call even%set_ft(.true.)
                call even%copy_fast(even_copy_cmat)
                call even%shellnorm
                call apply_opt_filter(even, param, cur_fil, .false.)
            endif
            call odd%get_rmat_ptr(rmat_odd) ! point rmat_odd to the filtered odd, rmat_even should be filtered even
            if( L_BENCH_GLOB )then
                rt_filter_even = rt_filter_even + toc(t_filter_even)
                t_search_opt   = tic()
            endif
            ! do the non-uniform, i.e. optimizing at each voxel
            if( params_glob%l_nonuniform )then
                ! 2D vs 3D cases
                if( dim3 == 1 )then
                    !$omp parallel do collapse(2) default(shared) private(k,l,k1,l1,k_ind,l_ind,ref_diff_odd, ref_diff_even) schedule(dynamic,CHUNKSZ) proc_bind(close)
                    do k = lb(1),ub(1)
                        do l = lb(2),ub(2)
                            ref_diff_odd  = 0.
                            ref_diff_even = 0.
                            ! applying an average window to each diff (eq 7 in the nonuniform paper)
                            do k_ind = 1, 2*params_glob%smooth_ext+1
                                k1 = k - params_glob%smooth_ext + k_ind - 1
                                do l_ind = 1, 2*params_glob%smooth_ext+1
                                    l1 = l - params_glob%smooth_ext + l_ind - 1
                                    if ((k1 >= 1 .and. k1 <= box) .and. (l1 >= 1 .and. l1 <= box)) then
                                        ref_diff_odd  = ref_diff_odd  + cur_diff_odd( k1,l1,1)*weights_2D(k_ind,l_ind)
                                        ref_diff_even = ref_diff_even + cur_diff_even(k1,l1,1)*weights_2D(k_ind,l_ind)
                                    endif
                                enddo
                            enddo
                            ! opt_diff keeps the minimized cost value at each voxel of the search
                            ! opt_odd  keeps the best voxel of the form B*odd
                            ! opt_even keeps the best voxel of the form B*even
                            if (ref_diff_odd < opt_diff_odd(k,l,1)) then
                                opt_odd(k,l,1)      = rmat_odd(k,l,1)
                                opt_diff_odd(k,l,1) = ref_diff_odd
                                opt_freq_odd(k,l,1) = cur_ind
                            endif
                            if (ref_diff_even < opt_diff_even(k,l,1)) then
                                opt_even(k,l,1)      = rmat_even(k,l,1)
                                opt_diff_even(k,l,1) = ref_diff_even
                                opt_freq_even(k,l,1) = cur_ind
                            endif
                        enddo
                    enddo
                    !$omp end parallel do
                else
                    !omp parallel do collapse(3) default(shared) private(k,l,m,k1,l1,m1,k_ind,l_ind,m_ind,ref_diff_odd,ref_diff_even) schedule(dynamic,CHUNKSZ) proc_bind(close)
                    do k = lb(1),ub(1)
                        do l = lb(2),ub(2)
                            do m = lb(3),ub(3)
                                ref_diff_odd  = 0.
                                ref_diff_even = 0.
                                ! applying an average window to each diff (eq 7 in the nonuniform paper)
                                do k_ind = 1, 2*params_glob%smooth_ext+1
                                    k1 = k - params_glob%smooth_ext + k_ind - 1
                                    do l_ind = 1, 2*params_glob%smooth_ext+1
                                        l1 = l - params_glob%smooth_ext + l_ind - 1
                                        do m_ind = 1, 2*params_glob%smooth_ext+1
                                            m1 = m - params_glob%smooth_ext + m_ind - 1
                                            if ((k1 >= 1 .and. k1 <= box) .and. (l1 >= 1 .and. l1 <= box) .and. (m1 >= 1 .and. m1 <= dim3)) then
                                                ref_diff_odd  = ref_diff_odd  + cur_diff_odd( k1,l1,m1)*weights_3D(k_ind,l_ind,m_ind)
                                                ref_diff_even = ref_diff_even + cur_diff_even(k1,l1,m1)*weights_3D(k_ind,l_ind,m_ind)
                                            endif
                                        enddo
                                    enddo
                                enddo
                                ! opt_diff keeps the minimized cost value at each voxel of the search
                                ! opt_odd  keeps the best voxel of the form B*odd
                                ! opt_even keeps the best voxel of the form B*even
                                if (ref_diff_odd < opt_diff_odd(k,l,m)) then
                                    opt_odd(k,l,m)      = rmat_odd(k,l,m)
                                    opt_diff_odd(k,l,m) = ref_diff_odd
                                    opt_freq_odd(k,l,m) = cur_ind
                                endif
                                if (ref_diff_even < opt_diff_even(k,l,m)) then
                                    opt_even(k,l,m)      = rmat_even(k,l,m)
                                    opt_diff_even(k,l,m) = ref_diff_even
                                    opt_freq_even(k,l,m) = cur_ind
                                endif
                            enddo
                        enddo
                    enddo
                    !omp end parallel do
                endif
                min_sum_odd  = sum(opt_diff_odd) ! TODO
                min_sum_even = sum(opt_diff_even) ! TODO
            else
                ! keep the theta which gives the lowest cost (over all voxels)
                if (sum(cur_diff_odd) < min_sum_odd) then
                    opt_odd      = rmat_odd
                    min_sum_odd  = sum(cur_diff_odd)
                    best_ind     = cur_ind
                    opt_freq_odd = cur_ind
                endif
                if (sum(cur_diff_even) < min_sum_even) then
                    opt_even      = rmat_even
                    min_sum_even  = sum(cur_diff_even)
                    opt_freq_even = cur_ind
                endif
            endif
            write(*, *) 'min cost val (odd) = ', min_sum_odd, '; current cost (odd) = ', sum(cur_diff_odd)
            if( L_BENCH_GLOB )then
                rt_search_opt = rt_search_opt + toc(t_search_opt)
            endif
        enddo
        if( L_BENCH_GLOB )then
            rt_tot     = toc(t_tot)
            benchfname = 'OPT_FILTER_BENCH.txt'
            call fopen(fnr, FILE=trim(benchfname), STATUS='REPLACE', action='WRITE')
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f9.2)') 'odd filtering        : ', rt_filter_odd
            write(fnr,'(a,1x,f9.2)') 'even filtering       : ', rt_filter_even
            write(fnr,'(a,1x,f9.2)') 'searching/optimizing : ', rt_search_opt
            write(fnr,'(a,1x,f9.2)') 'total time           : ', rt_tot
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
            write(fnr,'(a,1x,f9.2)') 'odd filtering        : ', (rt_filter_odd /rt_tot) * 100. 
            write(fnr,'(a,1x,f9.2)') 'even filtering       : ', (rt_filter_even/rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') 'searching/optimizing : ', (rt_search_opt /rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') '% accounted for          : ',&
            &((rt_filter_odd+rt_filter_even+rt_search_opt)/rt_tot) * 100.
            call fclose(fnr)
        endif
        if( .not. params_glob%l_nonuniform ) write(*, *) 'minimized cost at resolution = ', box*params_glob%smpd/best_ind
        call odd%set_rmat(opt_odd,   .false.)
        call even%set_rmat(opt_even, .false.)
        ! output the optimized frequency map to see the nonuniform parts
        if( dim3 > 1 )then
            if( params_glob%l_nonuniform )then
                file_tag = 'nonuniform_filter_'//trim(params_glob%filter)//'_ext_'//int2str(params_glob%smooth_ext)
            else
                file_tag = 'uniform_filter_'//trim(params_glob%filter)//'_ext_'//int2str(params_glob%smooth_ext)
            endif
            call freq_img%set_rmat(opt_freq_odd, .false.)
            call freq_img%write('opt_freq_odd_map_'//trim(file_tag)//'.mrc')
            call freq_img%set_rmat(box*params_glob%smpd/opt_freq_odd, .false.) ! resolution map
            call freq_img%write('opt_resolution_odd_map_'//trim(file_tag)//'.mrc')
            call freq_img%kill
        endif
        deallocate(opt_odd, opt_even, cur_diff_odd, opt_diff_odd, cur_diff_even, opt_diff_even, cur_fil, weights_3D, weights_2D)
        call odd_copy_rmat%kill
        call odd_copy_cmat%kill
        call even_copy_rmat%kill
        call even_copy_cmat%kill
    end subroutine opt_filter

end module simple_opt_filter
    