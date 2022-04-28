module simple_butterworth
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
        real,    parameter :: AN(9) = (/ 1., 5.1258, 13.1371, 21.8462, 25.6884, 21.8462, 13.1371, 5.1258, 1./)
        complex, parameter :: J = (0, 1) ! Complex identity: j = sqrt(-1)
        complex :: Bn, Kn                ! Normalized Butterworth polynomial, its derivative and its reciprocal
        complex :: js                    ! frequency is multiplied by the complex identity j
        integer :: k
        Bn  = (0., 0.)
        if (s/fc < 100) then
            js  = J*s/fc
            do k = 0, n
                Bn  = Bn + AN(k+1)*js**k
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
        integer :: k, l, j, half_w, ldim3, freq_val
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
        real                :: mean_odd, mean_even
        mean_odd  =  sum(odd)/product(shape(odd))
        mean_even = sum(even)/product(shape(even))
        call squared_diff(odd-mean_odd, even-mean_even, diff)
        diff = diff/(sum(odd-mean_odd)**2 + sum(even-mean_even)**2)
    end subroutine normalized_squared_diff

    ! optimized voxelwise uniform/nonuniform filter, using the (low-pass/butterworth)
    subroutine opt_voxel_fil(odd, even, smpd, is_uniform, mskimg, map2filt)
        type(image),            intent(inout) :: odd
        type(image),            intent(inout) :: even
        real,                   intent(in)    :: smpd
        character(len=*),       intent(in)    :: is_uniform
        type(image),  optional, intent(inout) :: mskimg
        class(image), optional, intent(inout) :: map2filt
        type(image)          :: odd_copy, even_copy, map2filt_copy
        integer              :: k,l,m,max_lplim, box, dim3, ldim(3), find_start, find_stop, best_ind, cur_ind, k1,l1,m1,k_ind,l_ind,m_ind, lb(3), ub(3)
        real                 :: cur_min_sum, ref_diff, rad
        logical              :: map2filt_present, mskimg_present
        integer, parameter   :: CHUNKSZ=20, BW_ORDER=8, FIND_STEPSZ=2
        real,    parameter   :: LP_START = 30.                ! 30 A resolution
        integer, parameter   :: SPA_SUP = 0, MID = 1+SPA_SUP  ! support of the window function
        real,    pointer     :: rmat_odd(:,:,:)=>null(), rmat_even(:,:,:)=>null(), rmat_map2filt(:,:,:)=>null()
        real,    allocatable :: opt_odd(:,:,:), opt_even(:,:,:), cur_diff(:,:,:), opt_diff(:,:,:), but_fil(:), weights_3D(:,:,:), weights_2D(:,:), opt_map2filt(:,:,:)
        logical, allocatable :: l_mask(:,:,:)
        map2filt_present = present(map2filt)
        mskimg_present   = present(mskimg)
        if (mskimg_present) then
            l_mask  = mskimg%bin2logical()
        endif
        ldim    = odd%get_ldim()
        box     = ldim(1)
        dim3    = ldim(3)
        find_stop  = calc_fourier_index(2. * smpd, box, smpd)
        find_start = calc_fourier_index(LP_START, box, smpd)
        call odd_copy%copy(odd)
        call even_copy%copy(even)
        if (map2filt_present) then
            allocate(opt_map2filt(box,box,dim3), source=0.)
            call map2filt_copy%copy(map2filt)
        endif
        allocate(opt_odd(box,box,dim3), opt_even(box,box,dim3), cur_diff(box,box,dim3), opt_diff(box,box,dim3), but_fil(box), source=0.)
        ! assign the weights of the neighboring voxels
        allocate(weights_2D(SPA_SUP*2+1, SPA_SUP*2+1), weights_3D(SPA_SUP*2+1, SPA_SUP*2+1, SPA_SUP*2+1), source=0.)
        ! 2D weights
        do k = 1, 2*SPA_SUP+1
            do l = 1, 2*SPA_SUP+1
                rad = hyp(real(k-MID), real(l-MID))
                weights_2D(k,l) = -rad/(SPA_SUP + 1) + 1.  ! linear function: 1 at rad = 0 and 0 at rad = SPA_SUP + 1
                if (weights_2D(k,l) < 0.) then
                    weights_2D(k,l) = 0.
                endif
            enddo
        enddo
        ! 3D weights
        do k = 1, 2*SPA_SUP+1
            do l = 1, 2*SPA_SUP+1
                do m = 1, 2*SPA_SUP+1
                    rad = hyp(real(k-MID), real(l-MID), real(m-MID))
                    weights_3D(k,l,m) = -rad/(SPA_SUP + 1) + 1.  ! linear function: 1 at rad = 0 and 0 at rad = SPA_SUP + 1
                    if (weights_3D(k,l,m) < 0.) then
                        weights_3D(k,l,m) = 0.
                    endif
                enddo
            enddo
        enddo
        weights_3D = weights_3D/sum(weights_3D) ! weights has energy of 1
        weights_2D = weights_2D/sum(weights_2D) ! weights has energy of 1
        ! determine loop bounds for better load balancing in the following parallel loop
        if (mskimg_present) then
            call bounds_from_mask3D( l_mask, lb, ub )
        else
            lb = (/ 1, 1, 1/)
            ub = (/ box, box, box /)
        endif
        ! starting the searching for the best fourier index from here
        opt_diff     = 0.
        opt_diff(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) = huge(cur_min_sum)
        cur_min_sum  = huge(cur_min_sum)   
        best_ind     = find_start
        do cur_ind = find_start, find_stop, FIND_STEPSZ
            write(*, *) 'current Fourier index = ', cur_ind
            ! filtering odd
            call odd%copy(odd_copy)
            call odd%fft
            !call odd%lp(cur_ind)
            call butterworth_filter(but_fil, BW_ORDER, real(cur_ind))
            call odd%apply_filter(but_fil)
            call odd%ifft
            call odd%get_rmat_ptr(rmat_odd)
            call even%copy(even_copy)
            call even%get_rmat_ptr(rmat_even)
            call normalized_squared_diff(rmat_odd, rmat_even, cur_diff)
            ! filtering even using the same filter
            call even%fft
            call even%apply_filter(but_fil)
            call even%ifft
            call even%get_rmat_ptr(rmat_even)
            if (map2filt_present) then
                call map2filt%copy(map2filt_copy)
                call map2filt%fft
                call map2filt%apply_filter(but_fil)
                call map2filt%ifft
                call map2filt%get_rmat_ptr(rmat_map2filt)
            endif
            ! do the non-uniform, i.e. optimizing at each voxel
            if (is_uniform == 'no') then
                ! 2D vs 3D cases
                if (dim3 == 1) then
                    !$omp parallel do collapse(2) default(shared) private(k,l,k1,l1,k_ind,l_ind,ref_diff) schedule(dynamic,CHUNKSZ) proc_bind(close)
                    do k = lb(1),ub(1)
                        do l = lb(2),ub(2)
                            ref_diff = 0.
                            ! applying an average window to each diff (eq 7 in the nonuniform paper)
                            do k_ind = 1, 2*SPA_SUP+1
                                k1 = k - SPA_SUP + k_ind - 1
                                do l_ind = 1, 2*SPA_SUP+1
                                    l1 = l - SPA_SUP + l_ind - 1
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
                                if (map2filt_present) then
                                    opt_map2filt(k,l,1) = rmat_map2filt(k,l,1)
                                endif
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
                                do k_ind = 1, 2*SPA_SUP+1
                                    k1 = k - SPA_SUP + k_ind - 1
                                    do l_ind = 1, 2*SPA_SUP+1
                                        l1 = l - SPA_SUP + l_ind - 1
                                        do m_ind = 1, 2*SPA_SUP+1
                                            m1 = m - SPA_SUP + m_ind - 1
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
                                    if (map2filt_present) then
                                        opt_map2filt(k,l,m) = rmat_map2filt(k,l,m)
                                    endif
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
                    if (map2filt_present) then
                        opt_map2filt = rmat_map2filt
                    endif
                endif
            endif
            write(*, *) 'min cost val = ', cur_min_sum, '; current cost = ', sum(cur_diff)
        enddo
        if (is_uniform == 'yes') then
            write(*, *) 'minimized cost at index = ', best_ind
        endif
        call odd%set_rmat(opt_odd,   .false.)
        call even%set_rmat(opt_even, .false.)
        if (map2filt_present) then
            call map2filt%set_rmat(opt_map2filt, .false.)
        endif
        deallocate(l_mask, opt_odd, opt_even, cur_diff, opt_diff, but_fil, weights_3D, weights_2D)
    end subroutine opt_voxel_fil
end module simple_butterworth
    