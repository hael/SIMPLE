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
    ! PRIVATE
    procedure, private :: kernel_center
    procedure, private :: cosine_kernel
    procedure, private :: rbf_kernel
    procedure, private :: compute_eigvecs
end type

real, parameter :: C_CONST = 0.4_dp  ! for rbf_kernel for testing

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
        use simple_parameters, only: params_glob
        class(kpca_svd),   intent(inout) :: self
        real,              intent(in)    :: pcavecs(self%D,self%N)
        integer, optional, intent(in)    :: maxpcaits
        integer,  parameter :: MAX_ITS = 100
        real(dp), parameter :: TOL = 0.01_dp
        real(dp) :: ker(self%N,self%N), prev_data(self%D), cur_data(self%D), proj_data(self%N), s, denom,&
                    &sum_vecs(self%N), ker_weight(self%N,self%N), mat(self%N,self%D), eig_vecs(self%N,self%Q)
        integer  :: i, ind, iter, its
        mat = real(transpose(pcavecs),dp)
        ! compute the kernel
        select case(trim(params_glob%kpca_ker))
            case('rbf')
                call self%rbf_kernel(   real(pcavecs,dp), ker)
            case('cosine')
                call self%cosine_kernel(real(pcavecs,dp), ker)
        end select
        ! compute the sorted principle components of the kernel above
        call self%compute_eigvecs(ker, eig_vecs)
        ! pre-imaging:
        ! 1. projecting each image on kernel space
        ! 2. applying the principle components to the projected vector
        ! 3. computing the pre-image (of the image in step 1) using the result in step 2
        its = MAX_ITS
        if( present(maxpcaits) ) its = maxpcaits
        ker_weight = matmul(matmul(ker, eig_vecs), transpose(eig_vecs))
        select case(trim(params_glob%kpca_ker))
            case('rbf')
                !$omp parallel do default(shared) proc_bind(close) schedule(static) private(ind,cur_data,prev_data,iter,i,proj_data,s)
                do ind = 1, self%N
                    cur_data  = real(pcavecs(:,ind), dp)
                    prev_data = 0._dp
                    iter      = 1
                    do while( euclid(cur_data,prev_data) > TOL .and. iter < its )
                        prev_data = cur_data
                        ! 1. projecting each image on kernel space
                        do i = 1,self%N
                            proj_data(i) = euclid(prev_data, real(pcavecs(:,i),dp))**2
                        enddo
                        ! 2. applying the principle components to the projected vector
                        proj_data = dexp(-proj_data/real(self%Q,dp)/C_CONST) * ker_weight(:,ind)
                        ! 3. computing the pre-image (of the image in step 1) using the result in step 2
                        s = sum(proj_data)
                        do i = 1,self%D
                            cur_data(i) = sum(proj_data * mat(:,i))
                        enddo
                        if( s > DTINY ) cur_data = cur_data/s
                        iter = iter + 1
                    enddo
                    self%data(:,ind) = cur_data
                enddo
                !$omp end parallel do
            case('cosine')
                sum_vecs = 0._dp
                !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i)
                do i = 1,self%N
                    sum_vecs(i) = sum(real(pcavecs(:,i),dp)**2)
                enddo
                !$omp end parallel do
                !$omp parallel do default(shared) proc_bind(close) schedule(static) private(ind,cur_data,prev_data,iter,i,proj_data,s,denom)
                do ind = 1, self%N
                    cur_data  = real(pcavecs(:,ind), dp)
                    prev_data = 0._dp
                    iter      = 1
                    do while( euclid(cur_data,prev_data) > TOL .and. iter < its )
                        prev_data = cur_data
                        ! 1. projecting each image on kernel space
                        do i = 1,self%N
                            denom        = dsqrt(sum(prev_data**2) * sum_vecs(i))
                            proj_data(i) = 0._dp
                            if( denom > DTINY ) proj_data(i) = sum(prev_data * real(pcavecs(:,i),dp)) / denom
                        enddo
                        ! 2. applying the principle components to the projected vector
                        proj_data = proj_data * ker_weight(:,ind)
                        ! 3. computing the pre-image (of the image in step 1) using the result in step 2
                        s = sum(proj_data)
                        do i = 1,self%D
                            cur_data(i) = sum(proj_data * mat(:,i))
                        enddo
                        if( s > TINY ) cur_data = cur_data/s
                        iter = iter + 1
                    enddo
                    self%data(:,ind) = cur_data
                enddo
                !$omp end parallel do
        end select
    end subroutine master_kpca

    subroutine kernel_center( self, ker )
        class(kpca_svd), intent(inout) :: self
        real(dp),        intent(inout) :: ker(self%N,self%N)
        real(dp) :: ones(self%N,self%N), ones_ker(self%N,self%N)
        ones     = 1._dp / real(self%N,dp)
        ones_ker = matmul(ones, ker)
        ! Appendix D.2.2 Centering in Feature Space from Schoelkopf, Bernhard, Support vector learning, 1997
        ker = ker - ones_ker - transpose(ones_ker) + matmul(ones_ker, ones)
    end subroutine kernel_center

    subroutine cosine_kernel( self, mat, ker )
        class(kpca_svd), intent(inout) :: self
        real(dp),        intent(in)    :: mat(self%D,self%N)
        real(dp),        intent(out)   :: ker(self%N,self%N)
        integer  :: i, j
        real(dp) :: denom, row_sums2(self%N)
        ! squared cosine similarity between pairs of rows
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i)
        do i = 1,self%N
            row_sums2(i) = sum(mat(:,i)**2)
        enddo
        !$omp end parallel do
        !$omp parallel do collapse(2) default(shared) proc_bind(close) schedule(static) private(i,j,denom)
        do i = 1,self%N
            do j = 1,self%N
                denom    = dsqrt(row_sums2(i) * row_sums2(j))
                ker(i,j) = 0._dp
                if( denom > TINY ) ker(i,j) = sum(mat(:,i) * mat(:,j)) / denom
            enddo
        enddo
        !$omp end parallel do
        call self%kernel_center(ker)
    end subroutine cosine_kernel

    subroutine rbf_kernel( self, mat, ker )
        class(kpca_svd), intent(inout) :: self
        real(dp),        intent(in)    :: mat(self%D,self%N)
        real(dp),        intent(out)   :: ker(self%N,self%N)
        integer :: i, j
        ! squared euclidean distance between pairs of rows
        !$omp parallel do collapse(2) default(shared) proc_bind(close) schedule(static) private(i,j)
        do i = 1,self%N
            do j = 1,self%N
                ker(i,j) = euclid(mat(:,i), mat(:,j))**2
            enddo
        enddo
        !$omp end parallel do
        ! normalization and centering
        ker = dexp(-ker/real(self%Q,dp)/C_CONST)
        call self%kernel_center(ker)
    end subroutine rbf_kernel

    subroutine compute_eigvecs( self, ker, eig_vecs )
        class(kpca_svd), intent(inout) :: self
        real(dp),        intent(in)    :: ker(self%N,self%N)
        real(dp),        intent(inout) :: eig_vecs(self%N,self%Q)
        real(dp) :: eig_vals(self%Q), tmp_ker(self%N,self%N), eig_vecs_ori(self%N,self%Q)
        integer  :: i
        tmp_ker = ker
        do i = 1, self%N
            tmp_ker(:,i) = tmp_ker(:,i) - sum(tmp_ker(:,i))/real(self%N,dp)
        enddo
        ! computing eigvals/eigvecs
        call eigh(self%N, tmp_ker, self%Q, eig_vals, eig_vecs_ori)
        eig_vals = eig_vals**2 / real(self%N)
        do i = 1, self%Q
            eig_vecs(:,i) = eig_vecs_ori(:,i) / dsqrt(eig_vals(i))
        enddo
        ! reverse the sorted order of eig_vecs
        eig_vecs = eig_vecs(:,self%Q:1:-1)
    end subroutine compute_eigvecs

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