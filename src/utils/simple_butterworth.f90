
module simple_butterworth
!$ use omp_lib
!$ use omp_lib_kinds
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
        real                 :: val(2)
        real,    parameter :: AN(9) = (/ 1., 5.1258, 13.1371, 21.8462, 25.6884, 21.8462, 13.1371, 5.1258, 1./)
        complex, parameter :: J = (0, 1) ! Complex identity: j = sqrt(-1)
        complex :: Bn, dBn, Kn, dKn      ! Normalized Butterworth polynomial, its derivative and its reciprocal
        complex :: js                    ! frequency is multiplied by the complex identity j
        integer :: k
        val = [0., 0.]
        Bn  = (0., 0.)
        dBn = (0., 0.)
        if (s/fc < 100) then
            js  = j*s/fc
            do k = 0, n
                Bn  = Bn  +   AN(k+1)*js**k
                dBn = dBn + k*AN(k+1)*js**k
            end do
            dBn = -dBn/fc
            Kn  = 1/Bn
            dKn = -dBn/Bn/Bn
            val(1) = sqrt(real(Kn)**2 + aimag(Kn)**2)
            val(2) = real( Kn*conjg(dKn) )/val(1)
        else
            val(1) = epsilon(val(1))
            val(2) = epsilon(val(2))
        endif
    end function butterworth

    ! Compute the Butterworth kernel of the order n-th of width w
    ! with the cut-off frequency fc
    ! https://en.wikipedia.org/wiki/Butterworth_filter
    subroutine butterworth_kernel(ker, ker_der, w, n, fc)
        real,    intent(inout) :: ker    (:, :, :)
        real,    intent(inout) :: ker_der(:, :, :)
        integer, intent(in)    :: w
        integer, intent(in)    :: n
        real   , intent(in)    :: fc
        integer :: k, l, j, half_w, ldim3
        real    :: freq_val, val(2) ! current frequency value
        freq_val = 0
        half_w   = int(w/2)
        ldim3    = size(ker, 3)
        ! loop over pixels
        if( ldim3 == 1 )then
            !$omp parallel do collapse(2) default(shared) private(l,k,freq_val,val) schedule(static) proc_bind(close)
            do k = 1, w
                do l = 1, w
                    freq_val = hyp(real(k-half_w), real(l-half_w))
                    ! compute the value of Butterworth transfer function at current frequency value
                    val            = butterworth(freq_val, n, fc)
                    ker(k,l,1)     = val(1)
                    ker_der(k,l,1) = val(2)
                end do
            end do
            !$omp end parallel do
        else
            !$omp parallel do collapse(3) default(shared) private(l,j,k,freq_val,val) schedule(static) proc_bind(close)
            do k = 1, w
                do l = 1, w
                    do j = 1, ldim3
                        freq_val = hyp(real(k-half_w), real(l-half_w), real(j-half_w))
                        ! compute the value of Butterworth transfer function at current frequency value
                        val            = butterworth(freq_val, n, fc)
                        ker(k,l,j)     = val(1)
                        ker_der(k,l,j) = val(2)
                    end do
                end do
            end do
            !$omp end parallel do
        endif
    end subroutine butterworth_kernel

    subroutine find_butterworth_opt(cmat_odd, rmat_even, cur_mat, ldim, smpd, is_uniform)
        real,                          intent(in)    :: rmat_even(:,:,:)
        complex(kind=c_float_complex), intent(in)    :: cmat_odd(:,:,:)
        real,                          intent(inout) :: cur_mat(:,:,:)
        integer,                       intent(in)    :: ldim(3)
        real,                          intent(in)    :: smpd
        character(len=*),              intent(in)    :: is_uniform
        integer              :: k,l,m,n,k1,l1,m1,k_ind,l_ind,m_ind, max_sup, box, dim3
        real                 :: rad, cur_min_sum, ref_diff, sup, theta
        real   , parameter   :: A = 47.27, B = -0.1781, C = 7.69, D = -0.02228  ! Fitting constants (constructed in MATLAB) of theta(FT_support) = a*exp(b*x) + c*exp(d*x)
        real   , parameter   :: MIN_SUP = 0.5, RES_LB = 30.                     ! lower bound of resolution is 30 Angstrom, upper bound is nyquist, hence .5
        integer, parameter   :: N_SUP = 20                                      ! number of intervals between MIN_SUP and MAX_SUP
        integer, parameter   :: SPA_SUP = 2, MID = 1+SPA_SUP                    ! support of the window function
        integer, parameter   :: CHUNKSZ=20, BW_ORDER=8
        real,    allocatable :: rmat_ker(:,:,:), orig_ker(:,:,:), orig_ker_der(:,:,:), prev_diff(:,:,:), cur_diff(:,:,:), weights_3D(:,:,:), weights_2D(:,:)
        type(image)          :: ker_odd_img
        complex(kind=c_float_complex), pointer :: cmat_conv(:,:,:)=>null()
        box  = ldim(1)
        dim3 = ldim(3)
        allocate(rmat_ker(box,box,dim3), orig_ker(box,box,dim3), orig_ker_der( box,box,dim3),&
        &prev_diff(box,box,dim3), cur_diff(box,box,dim3), source=0.)
        call ker_odd_img%new([box,box,dim3], smpd)
        max_sup      = int(RES_LB/smpd)*2 ! multiplication factor depending on the definition of support, set to 2 for now
        prev_diff    = huge(rad)
        cur_min_sum  = huge(rad)     
        ! assign the weights of the neighboring voxels
        allocate(weights_2D(SPA_SUP*2+1, SPA_SUP*2+1))
        allocate(weights_3D(SPA_SUP*2+1, SPA_SUP*2+1, SPA_SUP*2+1))
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
        do n = 1, N_SUP
            sup   = MIN_SUP + (n-1.)*(max_sup-MIN_SUP)/(N_SUP-1.)
            theta = A*exp(B*sup) + C*exp(D*sup)
            write(*, *) 'support = ', sup, '; theta = ', theta
            call butterworth_kernel(orig_ker, orig_ker_der, box, BW_ORDER, theta)
            ! computing B_kernel 'convolve' odd
            call ker_odd_img%set_rmat(orig_ker, .false.)
            call ker_odd_img%fft()
            call ker_odd_img%get_cmat_ptr(cmat_conv)
            cmat_conv = cmat_conv * cmat_odd
            call ker_odd_img%ifft()
            call ker_odd_img%get_rmat_sub(rmat_ker)
            ! Normalize to energy of 1, so B*odd is comparable with even in the cost function
            rmat_ker = rmat_ker/sum(rmat_ker)
            cur_diff = (rmat_ker - rmat_even)**2
            ! do the non-uniform, i.e. optimizing at each voxel
            if (is_uniform == 'no') then
                ! 2D vs 3D cases
                if (dim3 == 1) then
                    !$omp parallel do collapse(2) default(shared) private(k,l,k1,l1,k_ind,l_ind,ref_diff) schedule(dynamic,CHUNKSZ) proc_bind(close)
                    do k = 1,box
                        do l = 1,box
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
                            ! prev_diff keeps the lowest cost value at each voxel of the search
                            ! cur_mat   keeps the best voxel of the form B*odd
                            if (ref_diff < prev_diff(k,l,1)) then
                                cur_mat(k,l,1)   = rmat_ker(k,l,1)
                                prev_diff(k,l,1) = ref_diff
                            endif
                        enddo
                    enddo
                    !$omp end parallel do
                else
                    !$omp parallel do collapse(3) default(shared) private(k,l,m,k1,l1,m1,k_ind,l_ind,m_ind,ref_diff) schedule(dynamic,CHUNKSZ) proc_bind(close)
                    do k = 1,box
                        do l = 1,box
                            do m = 1,box
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
                                ! prev_diff keeps the lowest cost value at each voxel of the search
                                ! cur_mat   keeps the best voxel of the form B*odd
                                if (ref_diff < prev_diff(k,l,m)) then
                                    cur_mat(k,l,m)   = rmat_ker(k,l,m)
                                    prev_diff(k,l,m) = ref_diff
                                endif
                            enddo
                        enddo
                    enddo
                    !$omp end parallel do
                endif
            else
                ! keep the theta which gives the lowest cost (over all voxels)
                if (sum(cur_diff) < cur_min_sum) then
                    cur_mat     = rmat_ker
                    cur_min_sum = sum(cur_diff)
                endif
            endif
        enddo
    end subroutine find_butterworth_opt
end module simple_butterworth
    