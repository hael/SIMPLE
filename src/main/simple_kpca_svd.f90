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
    procedure :: cosine_kernel
    procedure :: rbf_kernel
    procedure :: compute_eigvecs
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
        use simple_parameters, only: params_glob
        class(kpca_svd),   intent(inout) :: self
        real,              intent(in)    :: pcavecs(self%D,self%N)
        integer, optional, intent(in)    :: maxpcaits
        real, parameter   :: TOL = 0.0001, MAX_ITS = 100, C_CONST = 0.4 ! for rbf_kernel for testing
        real    :: mat(self%N,self%D), eig_vecs(self%N,self%Q), ker(self%N,self%N), ker_weight(self%N,self%N)
        real    :: prev_data(self%D), cur_data(self%D), proj_data(self%N), s, denom
        integer :: i, ind, iter, its
        mat = transpose(pcavecs)
        ! compute the kernel
        select case(trim(params_glob%kpca_ker))
            case('rbf')
                call self%rbf_kernel(pcavecs, ker)
            case('cosine')
                call self%cosine_kernel(pcavecs, ker)
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
                    cur_data  = pcavecs(:,ind)
                    prev_data = 0.
                    iter      = 1
                    do while( euclid(cur_data,prev_data) > TOL .and. iter < its )
                        prev_data = cur_data
                        ! 1. projecting each image on kernel space
                        do i = 1,self%N
                            proj_data(i) = euclid(prev_data, pcavecs(:,i))**2
                        enddo
                        ! 2. applying the principle components to the projected vector
                        proj_data = exp(-proj_data/real(self%Q)/C_CONST)
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
            case('cosine')
                !$omp parallel do default(shared) proc_bind(close) schedule(static) private(ind,cur_data,prev_data,iter,i,proj_data,s,denom)
                do ind = 1, self%N
                    cur_data  = pcavecs(:,ind)
                    prev_data = 0.
                    iter      = 1
                    do while( euclid(cur_data,prev_data) > TOL .and. iter < its )
                        prev_data = cur_data
                        ! 1. projecting each image on kernel space
                        do i = 1,self%N
                            denom        = sqrt(sum(prev_data**2) * sum(pcavecs(:,i)**2))
                            proj_data(i) = 0.
                            if( denom > TINY ) proj_data(i) = sum(prev_data * pcavecs(:,i)) / denom
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
        real,            intent(inout) :: ker(self%N,self%N)
        real :: ones(self%N,self%N), ones_ker(self%N,self%N)
        ones     = 1. / real(self%N)
        ones_ker = matmul(ones, ker)
        ! Appendix D.2.2 Centering in Feature Space from Schoelkopf, Bernhard, Support vector learning, 1997
        ker = ker - ones_ker - transpose(ones_ker) + matmul(ones_ker, ones)
    end subroutine kernel_center

    subroutine cosine_kernel( self, mat, ker )
        class(kpca_svd), intent(inout) :: self
        real,            intent(in)    :: mat(self%D,self%N)
        real,            intent(out)   :: ker(self%N,self%N)
        integer :: i, j
        real    :: denom, row_sums2(self%N)
        ! squared cosine similarity between pairs of rows
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i)
        do i = 1,self%N
            row_sums2(i) = sum(mat(:,i)**2)
        enddo
        !$omp end parallel do
        !$omp parallel do collapse(2) default(shared) proc_bind(close) schedule(static) private(i,j,denom)
        do i = 1,self%N
            do j = 1,self%N
                denom    = sqrt(row_sums2(i) * row_sums2(j))
                ker(i,j) = 0.
                if( denom > TINY ) ker(i,j) = sum(mat(:,i) * mat(:,j)) / denom
            enddo
        enddo
        !$omp end parallel do
        call self%kernel_center(ker)
    end subroutine cosine_kernel

    subroutine rbf_kernel( self, mat, ker )
        class(kpca_svd), intent(inout) :: self
        real,            intent(in)    :: mat(self%D,self%N)
        real,            intent(out)   :: ker(self%N,self%N)
        real,            parameter     :: C_CONST = 0.4 ! for rbf_kernel for testing
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
        ker = exp(-ker/real(self%Q)/C_CONST)
        call self%kernel_center(ker)
    end subroutine rbf_kernel

    subroutine compute_eigvecs( self, ker, eig_vecs )
        class(kpca_svd), intent(inout) :: self
        real,            intent(in)    :: ker(self%N,self%N)
        real,            intent(inout) :: eig_vecs(self%N,self%Q)
        real    :: eig_vecs_all(self%N,self%N), eig_vals(self%N), tmp(self%N,self%N), tmp_ker(self%N,self%N)
        integer :: i, inds(self%N)
        tmp_ker = ker
        do i = 1, self%N
            tmp_ker(:,i) = tmp_ker(:,i) - sum(tmp_ker(:,i))/real(self%N)
        enddo
        ! computing eigvals/eigvecs
        eig_vecs_all = tmp_ker
        call svdcmp(eig_vecs_all, eig_vals, tmp)
        eig_vals = eig_vals**2 / real(self%N)
        inds     = (/(i,i=1,self%N)/)
        call hpsort(eig_vals, inds)
        call reverse(eig_vals)
        call reverse(inds)
        eig_vecs_all = eig_vecs_all(:,inds)
        do i = 1, self%N
            eig_vecs(i,:) = eig_vecs_all(i,1:self%Q) / sqrt(eig_vals(1:self%Q))
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