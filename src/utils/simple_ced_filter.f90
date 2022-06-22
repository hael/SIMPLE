! coherence-enhancing diffusion (CED) filter
module simple_ced_filter
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_defs
use simple_image,      only: image
use simple_parameters, only: params_glob
implicit none
#include "simple_local_flags.inc"

contains
    subroutine ced_filter_2D(img, sigma, niter_in, step_size_in, rho_in, c1_in, c2_in)
        class(image),      intent(inout) :: img
        real,              intent(in)    :: sigma
        integer, optional, intent(in)    :: niter_in
        real,    optional, intent(in)    :: step_size_in, rho_in, c1_in, c2_in
        type(image)          :: img_ker, ker, J11, J12, J22
        integer              :: ldim(3), n, k, l, gaussian_ext, niter
        real                 :: smpd, step_size, rho, c1, c2, sh
        real,    parameter   :: RHO_DEF = 4., C1_DEF = 0.001, C2_DEF = 1., DT_DEF = 0.15
        integer, parameter   :: NITER_DEF = 10
        integer, allocatable :: pos_ind_1(:), pos_ind_2(:), neg_ind_1(:), neg_ind_2(:)
        real,    allocatable :: D1(:,:,:), D2(:,:,:), eig_val(:,:,:,:), eig_vec(:,:,:,:),&
                                &lambda(:,:,:,:), a(:,:,:), b(:,:,:), c(:,:,:), discrete_table(:,:,:,:,:)
        real,    pointer     :: J11_rmat(:,:,:)=>null(), J12_rmat(:,:,:)=>null(), J22_rmat(:,:,:)=>null(),&
                                &cur_img_rmat(:,:,:)=>null(), ker_rmat(:,:,:)=>null()
        step_size = DT_DEF
        rho       = RHO_DEF
        c1        = C1_DEF
        c2        = C2_DEF
        niter     = NITER_DEF
        if( present(step_size_in) ) step_size = step_size_in
        if( present(rho_in) )       rho       = rho_in
        if( present(c1_in) )        c1        = c1_in
        if( present(c2_in) )        c2        = c2_in
        if( present(niter_in) )     niter     = niter_in
        ldim = img%get_ldim()
        smpd = img%get_smpd()
        if( ldim(3) > 1 ) THROW_HARD('This ced_filter_2D is strictly for 2D case only!')
        call ker%new(ldim, smpd)
        call J11%new(ldim, smpd)
        call J12%new(ldim, smpd)
        call J22%new(ldim, smpd)
        call img_ker%new(ldim, smpd)
        allocate(D1(ldim(1),ldim(2),1), D2(ldim(1),ldim(2),1), &
                &eig_val(ldim(1),ldim(2),1,2),eig_vec(ldim(1),ldim(2),1,2),&
                &lambda(ldim(1),ldim(2),1,2),discrete_table(ldim(1),ldim(2),1,3,3),&
                &a(ldim(1),ldim(2),1), b(ldim(1),ldim(2),1), c(ldim(1),ldim(2),1),&
                &source=0.)
        allocate(neg_ind_1(ldim(1)), pos_ind_1(ldim(1)), neg_ind_2(ldim(2)), pos_ind_2(ldim(2)))
        neg_ind_1 = [ldim(1), (k, k = 1,ldim(1)-1)]
        neg_ind_2 = [ldim(2), (k, k = 1,ldim(2)-1)]
        pos_ind_1 = [(k, k = 2,ldim(1)), 1]
        pos_ind_2 = [(k, k = 2,ldim(2)), 1]
        call img%get_rmat_ptr(cur_img_rmat)
        do n = 1, niter
            call img_ker%copy_fast(img)
            call ker%zero_and_unflag_ft()
            ! build the Gaussian kernel with sigma
            gaussian_ext = ceiling(2*sigma)
            do k = -gaussian_ext, gaussian_ext
            do l = -gaussian_ext, gaussian_ext
                sh = hyp(real(k),real(l))
                call ker%set_rmat_at(ldim(1)/2 + k, ldim(2)/2 + l, 1, exp(-(sh**2/(2*sigma**2)))) 
            enddo
            enddo
            call ker%get_rmat_ptr(ker_rmat)
            ker_rmat = ker_rmat/sum(ker_rmat)
            ! convolving img with the kernel
            call img_ker%fft()
            call ker%fft()
            img_ker = img_ker*ker
            call img_ker%ifft()
            ! computing the gradient using central difference scheme
            call img_ker%gradient(D1, D2)
            ! scalling D1, D2 due to fft down-scaling
            D1 = D1*product(ldim)
            D2 = D2*product(ldim)
            ! build the Gaussian kernel with rho
            gaussian_ext = ceiling(3*rho)
            call ker%zero_and_unflag_ft()
            do k = -gaussian_ext, gaussian_ext
            do l = -gaussian_ext, gaussian_ext
                sh = hyp(real(k),real(l))
                call ker%set_rmat_at(ldim(1)/2 + k, ldim(2)/2 + l, 1, exp(-(sh**2/(2*rho**2)))) 
            enddo
            enddo
            call ker%get_rmat_ptr(ker_rmat)
            ker_rmat = ker_rmat/sum(ker_rmat)
            ! construct the structure tensor
            call J11%set_rmat(D1**2, .false.)
            call J12%set_rmat(D1*D2, .false.)
            call J22%set_rmat(D2**2, .false.)
            call J11%fft()
            call J12%fft()
            call J22%fft()
            call ker%fft()
            J11 = J11*ker
            call J11%ifft()
            call J11%get_rmat_ptr(J11_rmat)
            J11_rmat = J11_rmat*product(ldim)
            J12 = J12*ker
            call J12%ifft()
            call J12%get_rmat_ptr(J12_rmat)
            J12_rmat = J12_rmat*product(ldim)
            J22 = J22*ker
            call J22%ifft()
            call J22%get_rmat_ptr(J22_rmat)
            J22_rmat = J22_rmat*product(ldim)
            ! computing eigenvalues/eigenvectors of the structure tensor
            eig_val(:,:,:,1) = (J11_rmat + J22_rmat + sqrt((J11_rmat - J22_rmat)**2 + 4*J12_rmat**2))/2.
            eig_val(:,:,:,2) = (J11_rmat + J22_rmat - sqrt((J11_rmat - J22_rmat)**2 + 4*J12_rmat**2))/2.
            lambda( :,:,:,1) = c1
            lambda( :,:,:,2) = c1 + (1-c1)*exp(-c2/(eig_val(:,:,:,1) - eig_val(:,:,:,2))**2)
            eig_vec(:,:,:,1) = 2*J12_rmat
            eig_vec(:,:,:,2) = J22_rmat - J11_rmat + sqrt((J11_rmat - J22_rmat)**2 + 4*J12_rmat**2)
            ! normalize the eigenvectors (critical, but not mentioned in the CED paper)
            do k = 1,ldim(1)
            do l = 1,ldim(2)
                sh = sqrt(eig_vec(k,l,1,1)**2 + eig_vec(k,l,1,2)**2)
                if( sh > epsilon(sh) ) eig_vec(k,l,1,:) = eig_vec(k,l,1,:)/sh
            enddo
            enddo
            a =  lambda(:,:,:,1)*eig_vec(:,:,:,1)**2 + lambda(:,:,:,2)*eig_vec(:,:,:,2)**2
            c =  lambda(:,:,:,1)*eig_vec(:,:,:,2)**2 + lambda(:,:,:,2)*eig_vec(:,:,:,1)**2
            b = (lambda(:,:,:,1) - lambda(:,:,:,2))*eig_vec(:,:,:,1)*eig_vec(:,:,:,2)
            ! solving the diffusion equations
            discrete_table(:,:,1,1,1) = (cur_img_rmat(neg_ind_1, pos_ind_2, 1) - cur_img_rmat(:,:,1))*&
                                    &( abs(b(neg_ind_1, pos_ind_2, 1)) - b(neg_ind_1, pos_ind_2, 1)+&
                                    &  abs(b(:,:,1)) - b(:,:,1) )/4.
            discrete_table(:,:,1,1,2) = (cur_img_rmat(:, pos_ind_2, 1) - cur_img_rmat(:,:,1))*&
                                    &( c(:, pos_ind_2, 1) + c(:,:,1) - abs(b(:,pos_ind_2,1)) - abs(b(:,:,1)) )/2.
            discrete_table(:,:,1,1,3) = (cur_img_rmat(pos_ind_1, pos_ind_2, 1) - cur_img_rmat(:,:,1))*&
                                    &( abs(b(pos_ind_1, pos_ind_2, 1)) + b(pos_ind_1, pos_ind_2, 1)+&
                                    &  abs(b(:,:,1)) + b(:,:,1) )/4.
            discrete_table(:,:,1,2,1) = (cur_img_rmat(neg_ind_1, :, 1) - cur_img_rmat(:,:,1))*&
                                    &( a(neg_ind_1, :, 1) + a(:,:,1) - abs(b(neg_ind_1,:,1)) - abs(b(:,:,1)) )/2.
            discrete_table(:,:,1,2,3) = (cur_img_rmat(pos_ind_1, :, 1) - cur_img_rmat(:,:,1))*&
                                    &( a(pos_ind_1, :, 1) + a(:,:,1) - abs(b(pos_ind_1,:,1)) - abs(b(:,:,1)) )/2.
            discrete_table(:,:,1,3,1) = (cur_img_rmat(neg_ind_1, neg_ind_2, 1) - cur_img_rmat(:,:,1))*&
                                    &( abs(b(neg_ind_1, neg_ind_2, 1)) + b(neg_ind_1, neg_ind_2, 1)+&
                                    &  abs(b(:,:,1)) + b(:,:,1) )/4.
            discrete_table(:,:,1,3,2) = (cur_img_rmat(:, neg_ind_2, 1) - cur_img_rmat(:,:,1))*&
                                    &( c(:, neg_ind_2, 1) + c(:,:,1) - abs(b(:,neg_ind_2,1)) - abs(b(:,:,1)) )/2.
            discrete_table(:,:,1,3,3) = (cur_img_rmat(pos_ind_1, neg_ind_2, 1) - cur_img_rmat(:,:,1))*&
                                    &( abs(b(pos_ind_1, neg_ind_2, 1)) - b(pos_ind_1, neg_ind_2, 1)+&
                                    &  abs(b(:,:,1)) - b(:,:,1) )/4.
            cur_img_rmat = cur_img_rmat + step_size*sum(sum(discrete_table,5),4)
        enddo
    end subroutine ced_filter_2D
end module simple_ced_filter