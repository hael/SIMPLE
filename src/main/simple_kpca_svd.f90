! kPCA using 'Learning to Find Pre-Images', using svd for eigvals/eigvecs
module simple_kpca_svd
include 'simple_lib.f08'
use simple_defs
use simple_pca, only: pca
implicit none

public :: kpca_svd
private

type, extends(pca) :: kpca_svd
    private
    real, allocatable :: E_zn(:,:)  !< expectations (feature vecs)
    real, allocatable :: data(:,:)  !< projected data on feature vecs
    logical           :: existence=.false.
    contains
    ! CONSTRUCTOR
    procedure :: new      => new_kpca
    ! GETTERS
    procedure :: get_feat => get_feat_kpca
    procedure :: generate => generate_kpca
    ! CALCULATORS
    procedure :: master   => master_kpca
    ! DESTRUCTOR
    procedure :: kill     => kill_kpca
    procedure :: kernel_center
    procedure :: rbf_kernel_1, rbf_kernel_2
    generic   :: rbf_kernel => rbf_kernel_1, rbf_kernel_2
    procedure :: compute_alpha
    procedure :: test_all, test_kernel_center, test_rbf_kernel, test_alpha
end type

logical :: L_PRINT = .false.

contains

    ! CONSTRUCTORS

    !>  \brief  is a constructor
    subroutine new_kpca( self, N, D, Q )
        class(kpca_svd), intent(inout) :: self
        integer,         intent(in)    :: N, D, Q
        call self%kill
        self%N = N
        self%D = D
        self%Q = Q
        ! allocate principal subspace and feature vectors
        allocate( self%E_zn(self%Q,self%N), self%data(self%D,self%N), source=0.)
        self%existence = .true.
    end subroutine new_kpca

    ! GETTERS

    pure integer function get_N( self )
        class(kpca_svd), intent(in) :: self
        get_N = self%N
    end function get_N

    pure integer function get_D( self )
        class(kpca_svd), intent(in) :: self
        get_D = self%D
    end function get_D

    pure integer function get_Q( self )
        class(kpca_svd), intent(in) :: self
        get_Q = self%Q
    end function get_Q

    !>  \brief  is for getting a feature vector
    function get_feat_kpca( self, i ) result( feat )
        class(kpca_svd), intent(inout) :: self
        integer,    intent(in)    :: i
        real,       allocatable   :: feat(:)
        allocate(feat(self%Q), source=self%E_zn(:,i))
    end function get_feat_kpca

    !>  \brief  is for sampling the generative model at a given image index
    subroutine generate_kpca( self, i, avg, dat )
        class(kpca_svd), intent(inout) :: self
        integer,         intent(in)    :: i
        real,            intent(in)    :: avg(self%D)
        real,            intent(inout) :: dat(self%D)
        dat = avg + self%data(:,i)
    end subroutine generate_kpca

    ! CALCULATORS

    subroutine master_kpca( self, pcavecs, maxpcaits )
        class(kpca_svd),   intent(inout) :: self
        real,              intent(in)    :: pcavecs(self%D,self%N)
        integer, optional, intent(in)    :: maxpcaits ! redundant
    end subroutine master_kpca

    subroutine kernel_center( self, ker, cen_ker )
        class(kpca_svd), intent(inout) :: self
        real,            intent(in)    :: ker(self%N,self%N)
        real,            intent(inout) :: cen_ker(self%N,self%N)
        real :: ones(self%N,self%N)
        ones = 1. / real(self%N)
        ! Appendix D.2.2 Centering in Feature Space from Schoelkopf, Bernhard, Support vector learning, 1997
        cen_ker  = cen_ker - matmul(ones, ker) - matmul(cen_ker, ones) + matmul(matmul(ones, ker), ones)
    end subroutine kernel_center

    subroutine rbf_kernel_1( self, mat, c, ker )
        class(kpca_svd), intent(inout) :: self
        real,            intent(in)    :: mat(self%N,self%D)
        real,            intent(in)    :: c
        real,            intent(out)   :: ker(self%N,self%N)
        integer :: i, j
        ! squared euclidean distance between pairs of rows
        do i = 1,self%N
            do j = 1,self%N
                ker(i,j) = euclid(mat(i,:), mat(j,:))**2
            enddo
        enddo
        ! normalization and centering
        ker = exp(-ker/real(self%Q)/c)
        call self%kernel_center(ker, ker)
    end subroutine rbf_kernel_1

    ! rbf kernel from the previous kernel
    subroutine rbf_kernel_2( self, mat_test, mat_train, c, ker_train, ker )
        class(kpca_svd), intent(inout) :: self
        real,            intent(in)    :: mat_test( self%N,self%D)
        real,            intent(in)    :: mat_train(self%N,self%D)
        real,            intent(in)    :: c
        real,            intent(in)    :: ker_train(self%N,self%N)
        real,            intent(out)   :: ker(self%N,self%N)
        integer :: i, j
        ! squared euclidean distance between pairs of rows
        do i = 1,self%N
            do j = 1,self%N
                ker(i,j) = euclid(mat_test(i,:), mat_train(j,:))**2
            enddo
        enddo
        ! normalization and centering
        ker = exp(-ker/real(self%Q)/c)
        call self%kernel_center(ker_train, ker)
    end subroutine rbf_kernel_2

    subroutine compute_alpha( self, ker, alpha )
        class(kpca_svd), intent(inout) :: self
        real,            intent(in)    :: ker(self%N,self%N)
        real,            intent(out)   :: alpha(self%N,self%Q)
        real    :: eig_vecs(self%N,self%N), eig_vals(self%N), tmp(self%N,self%N), tmp_ker(self%N,self%N)
        integer :: i, inds(self%N)
        tmp_ker = ker
        do i = 1, self%N
            tmp_ker(:,i) = tmp_ker(:,i) - sum(tmp_ker(:,i))/real(self%N)
        enddo
        ! computing eigvals/eigvecs
        eig_vecs = tmp_ker
        call svdcmp(eig_vecs, eig_vals, tmp)
        eig_vals  = eig_vals**2 / real(self%N)
        inds      = (/(i,i=1,self%N)/)
        call hpsort(eig_vals, inds)
        call reverse(eig_vals)
        call reverse(inds)
        tmp = tmp(:, inds)
        do i = 1, self%N
            alpha(i,:) = tmp(i,1:self%Q) / sqrt(eig_vals(1:self%Q))
        enddo
        ! change the sorted order of alpha
        alpha = alpha(:,self%Q:1:-1)
    end subroutine compute_alpha

    subroutine compute_gamma( self, mat, gamma )
        class(kpca_svd), intent(inout) :: self
        real,            intent(in)    :: mat(self%N,self%D)
        real,            intent(out)   :: gamma(self%N,self%Q)
    end subroutine compute_gamma

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill_kpca( self )
        class(kpca_svd), intent(inout) :: self
        if( self%existence )then
            deallocate( self%E_zn, self%data )
            self%existence = .false.
        endif
    end subroutine kill_kpca

    ! UNIT TESTS
    subroutine test_all( self )
        class(kpca_svd), intent(inout) :: self
        call self%test_kernel_center
        call self%test_rbf_kernel
        call self%test_alpha
    end subroutine test_all

    subroutine test_kernel_center( self )
        class(kpca_svd), intent(inout) :: self
        real :: ker(2,2), truth(2,2)
        self%N   = 2
        self%D   = 2
        self%Q   = 2
        ker(1,:) = [ 1., 2.]
        ker(2,:) = [ 3., 1.]
        call self%kernel_center(ker, ker)
        truth(1,:) = [ -0.75,  0.75 ]
        truth(2,:) = [  0.75, -0.75 ]
        if( maxval(abs(ker - truth)) < TINY )then
            print *, 'Unit test of kernel_center: PASSED'
        else
            print *, 'Unit test of kernel_center: FAILED'
            print *, 'ker = ', ker
        endif
    end subroutine test_kernel_center

    subroutine test_rbf_kernel( self )
        class(kpca_svd), intent(inout) :: self
        real :: mat(2,2), truth(2,2), ker(2,2)
        self%N   = 2
        self%D   = 2
        self%Q   = 2
        mat(1,:) = [ 0., 2.]
        mat(2,:) = [ 1., 1.]
        call self%rbf_kernel(mat, 1., ker)
        truth(1,:) = [  0.31606028, -0.31606028 ]
        truth(2,:) = [ -0.31606028,  0.31606028 ]
        if( maxval(abs(ker - truth)) < TINY )then
            print *, 'Unit test of rbf_kernel: PASSED'
        else
            print *, 'Unit test of rbf_kernel: FAILED'
            print *, 'ker = ', ker
        endif
    end subroutine test_rbf_kernel

    subroutine test_alpha( self )
        class(kpca_svd), intent(inout) :: self
        real :: mat(4,3), truth(4,2), alpha(4,2), ker(4,4), beta(4,2), gamma(4,4)
        real :: z_prev(3), z_cur(3), zcoeff(4), s
        integer :: i, ind
        self%N   = 4
        self%D   = 3
        self%Q   = 2
        mat(1,:) = [  0.,  1.,  7.]
        mat(2,:) = [  2.,  5.,  3.]
        mat(3,:) = [  1., -2.,  0.]
        mat(4,:) = [ -2.,  4.,  5.]
        call self%rbf_kernel(mat, 1., ker)
        call self%compute_alpha(ker, alpha)
        beta  = matmul(ker, alpha)
        gamma = matmul(beta, transpose(alpha))
        truth(1,1) =  0.31606028
        truth(2,1) = -0.31606028
        if( maxval(abs(alpha - truth)) < TINY )then
            print *, 'Unit test of compute_alpha: PASSED'
        else
            print *, 'Unit test of compute_alpha: FAILED'
            print *, 'ker = ', ker(:,1)
            print *, 'alpha = ', alpha(:,1)
            print *, 'alpha = ', alpha(:,2)
            print *, 'beta = ', beta(:,1)
            print *, 'gamma = ', gamma(:,1)
        endif
        ! testing pre-imaging iterations
        do ind = 1, self%N
            z_cur  = mat(ind,:)
            z_prev = 0.
            do while( euclid(z_cur,z_prev) > 0.01 )
                z_prev = z_cur
                do i = 1,self%N
                    zcoeff(i) = euclid(z_prev, mat(i,:))**2
                enddo
                zcoeff = exp(-zcoeff/real(self%Q)/1.)    ! c = 1.
                zcoeff = gamma(:,ind) * zcoeff
                s      = sum(zcoeff)
                do i = 1,self%D
                    z_cur(i) = sum(mat(:,i) * zcoeff)
                enddo
                z_cur = z_cur/s
            enddo
            self%data(:,ind) = z_cur
            print *, z_cur
        enddo
    end subroutine test_alpha

end module simple_kpca_svd