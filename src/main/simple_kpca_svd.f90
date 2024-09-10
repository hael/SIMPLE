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
    procedure :: kernel_center
    procedure :: rbf_kernel
    procedure :: cosine_kernel
    procedure :: compute_alpha
    ! DESTRUCTOR
    procedure :: kill     => kill_kpca
end type

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
        real, parameter   :: TOL = 0.0001, MAX_ITS = 100
        real    :: mat(self%N,self%D), alpha(self%N,self%Q), ker(self%N,self%N), gamma_t(self%N,self%N)
        real    :: prev_data(self%D), cur_data(self%D), coeffs(self%N), s, denom
        integer :: i, ind, iter
        mat = transpose(pcavecs)
        call self%cosine_kernel(pcavecs, pcavecs, ker)
        call self%compute_alpha(ker, alpha)
        gamma_t = matmul(matmul(ker, alpha), transpose(alpha))
        ! pre-imaging iterations
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(ind,cur_data,prev_data,iter,i,coeffs,s,denom)
        do ind = 1, self%N
            cur_data  = pcavecs(:,ind)
            prev_data = 0.
            iter      = 1
            do while( euclid(cur_data,prev_data) > TOL .and. iter < MAX_ITS )
                prev_data = cur_data
                do i = 1,self%N
                    denom     = sqrt(sum(prev_data**2) * sum(pcavecs(:,i)**2))
                    coeffs(i) = 0.
                    if( denom > TINY ) coeffs(i) = sum(prev_data * pcavecs(:,i)) / denom
                enddo
                coeffs = gamma_t(:,ind) * coeffs
                s      = sum(coeffs)        ! CHECK FOR s = 0 ?
                do i = 1,self%D
                    cur_data(i) = sum(mat(:,i) * coeffs)
                enddo
                cur_data = cur_data/s
                iter     = iter + 1
            enddo
            self%data(:,ind) = cur_data
        enddo
        !$omp end parallel do
    end subroutine master_kpca

    subroutine kernel_center( self, ker, cen_ker )
        class(kpca_svd), intent(inout) :: self
        real,            intent(in)    :: ker(self%N,self%N)
        real,            intent(inout) :: cen_ker(self%N,self%N)
        real :: ones(self%N,self%N)
        ones = 1. / real(self%N)
        ! Appendix D.2.2 Centering in Feature Space from Schoelkopf, Bernhard, Support vector learning, 1997
        cen_ker = cen_ker - matmul(ones, ker) - matmul(cen_ker, ones) + matmul(matmul(ones, ker), ones)
    end subroutine kernel_center

    ! rbf kernel from the previous kernel
    subroutine rbf_kernel( self, mat_test, mat_train, c, ker_train, ker )
        class(kpca_svd), intent(inout) :: self
        real,            intent(in)    :: mat_test( self%D,self%N)
        real,            intent(in)    :: mat_train(self%D,self%N)
        real,            intent(in)    :: c
        real,            intent(inout) :: ker_train(self%N,self%N)
        real,            intent(inout) :: ker(self%N,self%N)
        integer :: i, j
        ! squared euclidean distance between pairs of rows
        do i = 1,self%N
            do j = 1,self%N
                ker(i,j) = euclid(mat_test(:,i), mat_train(:,j))**2
            enddo
        enddo
        ! normalization and centering
        ker = exp(-ker/real(self%Q)/c)
        call self%kernel_center(ker_train, ker)
    end subroutine rbf_kernel

    subroutine cosine_kernel( self, mat_test, mat_train, ker )
        class(kpca_svd), intent(inout) :: self
        real,            intent(in)    :: mat_test( self%D,self%N)
        real,            intent(in)    :: mat_train(self%D,self%N)
        real,            intent(out)   :: ker(self%N,self%N)
        integer :: i, j
        real    :: denom
        ! squared cosine similarity between pairs of rows
        do i = 1,self%N
            do j = 1,self%N
                denom    = sqrt(sum(mat_test(:,i)**2) * sum(mat_train(:,j)**2))
                ker(i,j) = 0.
                if( denom > TINY ) ker(i,j) = sum(mat_test(:,i) * mat_train(:,j)) / denom
            enddo
        enddo
        call self%kernel_center(ker, ker)
    end subroutine cosine_kernel

    subroutine compute_alpha( self, ker, alpha )
        class(kpca_svd), intent(inout) :: self
        real,            intent(in)    :: ker(self%N,self%N)
        real,            intent(inout) :: alpha(self%N,self%Q)
        real    :: eig_vecs(self%N,self%N), eig_vals(self%N), tmp(self%N,self%N), tmp_ker(self%N,self%N)
        integer :: i, inds(self%N)
        tmp_ker = ker
        do i = 1, self%N
            tmp_ker(:,i) = tmp_ker(:,i) - sum(tmp_ker(:,i))/real(self%N)
        enddo
        ! computing eigvals/eigvecs
        eig_vecs = tmp_ker
        call svdcmp(eig_vecs, eig_vals, tmp)
        eig_vals = eig_vals**2 / real(self%N)
        inds     = (/(i,i=1,self%N)/)
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

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill_kpca( self )
        class(kpca_svd), intent(inout) :: self
        if( self%existence )then
            deallocate( self%E_zn, self%data )
            self%existence = .false.
        endif
    end subroutine kill_kpca

end module simple_kpca_svd