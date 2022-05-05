! optimization(search)-based filter (uniform/nonuniform)
module simple_opt_filter
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_defs
use simple_image, only: image
use simple_math,  only: hyp
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

    subroutine squared_diff(odd, even, diff)
        real, intent(in)    :: odd(:,:,:)
        real, intent(in)    :: even(:,:,:)
        real, intent(inout) :: diff(:,:,:)
        diff = (odd - even)**2
    end subroutine squared_diff

    ! normalized to 1 then take the squared diff
    subroutine same_energy_squared_diff(odd, even, diff)
        real, intent(in)    :: odd(:,:,:)
        real, intent(in)    :: even(:,:,:)
        real, intent(inout) :: diff(:,:,:)
        call squared_diff(odd/sum(odd), even/sum(even), diff)
    end subroutine same_energy_squared_diff

    ! from https://stats.stackexchange.com/questions/136232/definition-of-normalized-euclidean-distance#:~:text=The%20normalized%20squared%20euclidean%20distance,not%20related%20to%20Mahalanobis%20distance
    subroutine normalized_squared_diff(odd, even, diff)
        real, intent(in)    :: odd(:,:,:)
        real, intent(in)    :: even(:,:,:)
        real, intent(inout) :: diff(:,:,:)
        real :: mean_odd, mean_even
        mean_odd  = sum(odd)/product(shape(odd))
        mean_even = sum(even)/product(shape(even))
        call squared_diff(odd-mean_odd, even-mean_even, diff)
        diff = diff/(sum(odd-mean_odd)**2 + sum(even-mean_even)**2)
    end subroutine normalized_squared_diff

    subroutine apply_opt_filter(img, filter_type, param, cur_fil, use_cache)
        use simple_tvfilter, only: tvfilter
        type(image),      intent(inout) :: img
        character(len=*), intent(in)    :: filter_type
        real,             intent(in)    :: param
        real,             intent(inout) :: cur_fil(:)
        logical,          intent(in)    :: use_cache
        integer                         :: bw_order, io_stat
        type(tvfilter)                  :: tvfilt
        call img%fft()
        if (filter_type == 'lp') then
            call img%lp(int(param))
        elseif (filter_type == 'tv') then
            call tvfilt%new
            call tvfilt%apply_filter_3d(img, param)
            call tvfilt%kill
        else    ! default to butterworth8, even if wrong filter type is entered
            ! extract butterworth order number
            bw_order = 8
            if (filter_type(1:11) == 'butterworth') then
                call str2int(filter_type(12:len_trim(filter_type)), io_stat, bw_order)
                if (bw_order < 1 .or. bw_order > 10) then
                    bw_order = 8
                endif
            endif
            if( .not. use_cache )then
                call butterworth_filter(cur_fil, bw_order, param)
            endif
            call img%apply_filter(cur_fil)
        endif
        call img%ifft()
    end subroutine apply_opt_filter

    ! optimization(search)-based uniform/nonuniform filter, using the (low-pass/butterworth)
    subroutine opt_filter(odd, even, smpd, is_uniform, smooth_ext, filter_type, max_res, nsearch, mskimg, map2filt)
        type(image),            intent(inout) :: odd
        type(image),            intent(inout) :: even
        real,                   intent(in)    :: smpd
        character(len=*),       intent(in)    :: is_uniform
        integer,                intent(in)    :: smooth_ext
        character(len=*),       intent(in)    :: filter_type
        real,                   intent(in)    :: max_res
        integer,                intent(in)    :: nsearch
        type(image),  optional, intent(inout) :: mskimg
        class(image), optional, intent(inout) :: map2filt
        type(image)          :: odd_copy, even_copy, map2filt_copy, freq_img
        integer              :: k,l,m, box, dim3, ldim(3), find_start, find_stop, iter_no
        integer              :: best_ind, cur_ind, k1,l1,m1,k_ind,l_ind,m_ind, lb(3), ub(3), mid_ext
        real                 :: cur_min_sum, ref_diff, rad, param, find_stepsz
        logical              :: map2filt_present, mskimg_present
        character(len=90)    :: file_tag
        integer, parameter   :: CHUNKSZ = 20
        real,    parameter   :: LAMBDA_MIN = .1, LAMBDA_MAX = 2.
        real,    pointer     :: rmat_odd(:,:,:)=>null(), rmat_even(:,:,:)=>null(), rmat_map2filt(:,:,:)=>null()
        real,    allocatable :: opt_odd(:,:,:), opt_even(:,:,:), cur_diff(:,:,:), opt_diff(:,:,:), opt_freq(:,:,:)
        real,    allocatable :: cur_fil(:), weights_3D(:,:,:), weights_2D(:,:), opt_map2filt(:,:,:)
        logical, allocatable :: l_mask(:,:,:)
        map2filt_present = present(map2filt)
        mskimg_present   = present(mskimg)
        if( mskimg_present )then
            l_mask = mskimg%bin2logical()
        endif
        ldim        = odd%get_ldim()
        box         = ldim(1)
        dim3        = ldim(3)
        mid_ext     = 1 + smooth_ext
        find_stop   = calc_fourier_index(2. * smpd, box, smpd)
        find_start  = calc_fourier_index(max_res, box, smpd)
        find_stepsz = real(find_stop - find_start)/(nsearch-1)
        call freq_img%new([box,box,dim3], smpd)
        call odd_copy%copy(odd)
        call even_copy%copy(even)
        if( map2filt_present )then
            allocate(opt_map2filt(box,box,dim3), source=0.)
            call map2filt_copy%copy(map2filt)
        endif
        allocate(opt_odd(box,box,dim3), opt_even(box,box,dim3), cur_diff(box,box,dim3), opt_diff(box,box,dim3),opt_freq(box,box,dim3),&
        &cur_fil(box),weights_2D(smooth_ext*2+1, smooth_ext*2+1), weights_3D(smooth_ext*2+1, smooth_ext*2+1, smooth_ext*2+1), source=0.)
        ! assign the weights of the neighboring voxels
        ! 2D weights
        do k = 1, 2*smooth_ext+1
            do l = 1, 2*smooth_ext+1
                rad = hyp(real(k-mid_ext), real(l-mid_ext))
                weights_2D(k,l) = -rad/(smooth_ext + 1) + 1.  ! linear function: 1 at rad = 0 and 0 at rad = smooth_ext + 1
                if (weights_2D(k,l) < 0.) then
                    weights_2D(k,l) = 0.
                endif
            enddo
        enddo
        ! 3D weights
        do k = 1, 2*smooth_ext+1
            do l = 1, 2*smooth_ext+1
                do m = 1, 2*smooth_ext+1
                    rad = hyp(real(k-mid_ext), real(l-mid_ext), real(m-mid_ext))
                    weights_3D(k,l,m) = -rad/(smooth_ext + 1) + 1.  ! linear function: 1 at rad = 0 and 0 at rad = smooth_ext + 1
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
        ! searching for the best fourier index from here
        opt_diff     = 0.
        opt_freq     = 0.   ! record the optimized cutoff frequency
        cur_min_sum  = huge(cur_min_sum)   
        best_ind     = find_start
        iter_no      = 0
        opt_diff(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) = huge(cur_min_sum)
        opt_freq(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) = huge(cur_min_sum)
        do iter_no = 1, nsearch
            cur_ind = find_start + (iter_no - 1)*find_stepsz
            if( filter_type == 'tv' )then
                param = LAMBDA_MIN + (cur_ind - find_start)*(LAMBDA_MAX - LAMBDA_MIN)/(find_stop - find_start)
                write(*, *) '('//int2str(iter_no)//'/'//int2str(nsearch)//') current lambda = ', param
            else
                param = real(cur_ind)
                write(*, *) '('//int2str(iter_no)//'/'//int2str(nsearch)//') current Fourier index = ', param
            endif
            ! filtering odd
            call odd%copy_fast(odd_copy)
            call apply_opt_filter(odd, filter_type, param, cur_fil, .false.)
            call odd%get_rmat_ptr(rmat_odd)
            call even%copy_fast(even_copy)
            call even%get_rmat_ptr(rmat_even)
            call squared_diff(rmat_odd, rmat_even, cur_diff)
            ! filtering even using the same filter
            call apply_opt_filter(even, filter_type, param, cur_fil, .true.)
            call even%get_rmat_ptr(rmat_even)
            if( map2filt_present )then
                call map2filt%copy_fast(map2filt_copy)
                call apply_opt_filter(map2filt, filter_type, param, cur_fil, .true.)
                call map2filt%get_rmat_ptr(rmat_map2filt)
            endif
            ! do the non-uniform, i.e. optimizing at each voxel
            if( is_uniform == 'no')then
                ! 2D vs 3D cases
                if( dim3 == 1 )then
                    !$omp parallel do collapse(2) default(shared) private(k,l,k1,l1,k_ind,l_ind,ref_diff) schedule(dynamic,CHUNKSZ) proc_bind(close)
                    do k = lb(1),ub(1)
                        do l = lb(2),ub(2)
                            ref_diff = 0.
                            ! applying an average window to each diff (eq 7 in the nonuniform paper)
                            do k_ind = 1, 2*smooth_ext+1
                                k1 = k - smooth_ext + k_ind - 1
                                do l_ind = 1, 2*smooth_ext+1
                                    l1 = l - smooth_ext + l_ind - 1
                                    if ((k1 >= 1 .and. k1 <= box) .and. (l1 >= 1 .and. l1 <= box)) then
                                        ref_diff = ref_diff + cur_diff(k1,l1,1)*weights_2D(k_ind,l_ind)
                                    endif
                                enddo
                            enddo
                            ! opt_diff keeps the minimized cost value at each voxel of the search
                            ! opt_odd  keeps the best voxel of the form B*odd
                            ! opt_even keeps the best voxel of the form B*even
                            if (ref_diff < opt_diff(k,l,1)) then
                                opt_odd(k,l,1)  = rmat_odd(k,l,1)
                                opt_even(k,l,1) = rmat_even(k,l,1)
                                opt_diff(k,l,1) = ref_diff
                                opt_freq(k,l,1) = cur_ind
                                if( map2filt_present ) opt_map2filt(k,l,1) = rmat_map2filt(k,l,1)
                            endif
                        enddo
                    enddo
                    !$omp end parallel do
                else
                    !$omp parallel do collapse(3) default(shared) private(k,l,m,k1,l1,m1,k_ind,l_ind,m_ind,ref_diff) schedule(dynamic,CHUNKSZ) proc_bind(close)
                    do k = lb(1),ub(1)
                        do l = lb(2),ub(2)
                            do m = lb(3),ub(3)
                                ref_diff = 0.
                                ! applying an average window to each diff (eq 7 in the nonuniform paper)
                                do k_ind = 1, 2*smooth_ext+1
                                    k1 = k - smooth_ext + k_ind - 1
                                    do l_ind = 1, 2*smooth_ext+1
                                        l1 = l - smooth_ext + l_ind - 1
                                        do m_ind = 1, 2*smooth_ext+1
                                            m1 = m - smooth_ext + m_ind - 1
                                            if ((k1 >= 1 .and. k1 <= box) .and. (l1 >= 1 .and. l1 <= box) .and. (m1 >= 1 .and. m1 <= dim3)) then
                                                ref_diff = ref_diff + cur_diff(k1,l1,m1)*weights_3D(k_ind,l_ind,m_ind)
                                            endif
                                        enddo
                                    enddo
                                enddo
                                ! opt_diff keeps the minimized cost value at each voxel of the search
                                ! opt_odd  keeps the best voxel of the form B*odd
                                ! opt_even keeps the best voxel of the form B*even
                                if (ref_diff < opt_diff(k,l,m)) then
                                    opt_odd(k,l,m)  = rmat_odd(k,l,m)
                                    opt_even(k,l,m) = rmat_even(k,l,m)
                                    opt_diff(k,l,m) = ref_diff
                                    opt_freq(k,l,m) = cur_ind
                                    if(map2filt_present ) opt_map2filt(k,l,m) = rmat_map2filt(k,l,m)
                                endif
                            enddo
                        enddo
                    enddo
                    !$omp end parallel do
                endif
                cur_min_sum = sum(opt_diff)
            else
                ! keep the theta which gives the lowest cost (over all voxels)
                if (sum(cur_diff) < cur_min_sum) then
                    opt_odd      = rmat_odd
                    opt_even     = rmat_even
                    cur_min_sum  = sum(cur_diff)
                    best_ind     = cur_ind
                    opt_freq     = cur_ind
                    if( map2filt_present ) opt_map2filt = rmat_map2filt
                endif
            endif
            write(*, *) 'min cost val = ', cur_min_sum, '; current cost = ', sum(cur_diff)
        enddo
        if(is_uniform == 'yes' ) write(*, *) 'minimized cost at resolution = ', box*smpd/best_ind
        call odd%set_rmat(opt_odd,   .false.)
        call even%set_rmat(opt_even, .false.)
        ! output the optimized frequency map to see the nonuniform parts
        if (dim3 > 1) then
            file_tag = trim(is_uniform)//'_uniform_filter_'//trim(filter_type)//'_ext_'//int2str(smooth_ext)
            call freq_img%set_rmat(opt_freq, .false.)
            call freq_img%write('opt_freq_map_'//trim(file_tag)//'.mrc')
            call freq_img%set_rmat(box*smpd/opt_freq, .false.) ! resolution map
            call freq_img%write('opt_resolution_map_'//trim(file_tag)//'.mrc')
            call freq_img%kill
        endif
        if( map2filt_present ) call map2filt%set_rmat(opt_map2filt, .false.)
        deallocate(opt_odd, opt_even, cur_diff, opt_diff, cur_fil, weights_3D, weights_2D)
    end subroutine opt_filter
end module simple_opt_filter
    