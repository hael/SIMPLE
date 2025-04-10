! kPCA using 'Learning to Find Pre-Images', using svd for eigvals/eigvecs
module simple_kpca_svd
!$ use omp_lib
!$ use omp_lib_kinds
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

real, parameter :: C_CONST = 0.4  ! for rbf_kernel for testing

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
        integer,         intent(in)    :: i
        real,            allocatable   :: feat(:)
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
        integer, parameter :: MAX_ITS = 500
        real,    parameter :: TOL     = 0.0001
        logical, parameter :: DEBUG   = .false.
        integer(int64)     :: start_time, end_time
        real(real64)       :: rate
        real(dp) :: denom
        real     :: ker(self%N,self%N), eig_vecs(self%N,self%Q), norm_pcavecs(self%D,self%N),norm_pcavecs_t(self%N,self%D),&
                   &proj_data(self%N,params_glob%nthr), norm_prev(self%D,params_glob%nthr), norm_data(self%D,params_glob%nthr)
        integer  :: i, ind, iter, its, ithr
        ! compute the kernel
        select case(trim(params_glob%kpca_ker))
            case('rbf')
                call self%rbf_kernel(   pcavecs, ker)
            case('cosine')
                call self%cosine_kernel(pcavecs, ker)
        end select
        ! compute the sorted principle components of the kernel above
        if( DEBUG ) call system_clock(start_time, rate)
        call self%compute_eigvecs(ker, eig_vecs)
        if( DEBUG )then
            call system_clock(end_time)
            print *, "Eigh time: ", real(end_time-start_time)/real(rate), " seconds"
        endif
        ! pre-imaging:
        ! 1. projecting each image on kernel space
        ! 2. applying the principle components to the projected vector
        ! 3. computing the pre-image (of the image in step 1) using the result in step 2
        its = MAX_ITS
        if( present(maxpcaits) ) its = maxpcaits
        if( DEBUG ) call system_clock(start_time, rate)
        ker = matmul(matmul(ker, eig_vecs), transpose(eig_vecs))
        select case(trim(params_glob%kpca_ker))
            case('rbf')
                !$omp parallel do default(shared) proc_bind(close) schedule(static) private(ind,iter,ithr,i,denom)
                do ind = 1, self%N
                    self%data(:,ind)  = pcavecs(:,ind)
                    ithr              = omp_get_thread_num() + 1
                    norm_prev(:,ithr) = 0.
                    iter              = 1
                    do while( euclid(self%data(:,ind),norm_prev(:,ithr)) > TOL .and. iter < its )
                        norm_prev(:,ithr) = self%data(:,ind)
                        ! 1. projecting each image on kernel space
                        do i = 1,self%N
                            proj_data(i,ithr) = euclid(norm_prev(:,ithr), pcavecs(:,i))**2
                        enddo
                        ! 2. applying the principle components to the projected vector
                        proj_data(:,ithr) = exp(-proj_data(:,ithr)/real(self%Q)/C_CONST) * ker(:,ind)
                        ! 3. computing the pre-image (of the image in step 1) using the result in step 2
                        denom             = sum(real(proj_data(:,ithr),dp))
                        self%data(:,ind)  = matmul(pcavecs, proj_data(:,ithr))
                        if( denom > DTINY ) self%data(:,ind) = self%data(:,ind)/real(denom)
                        iter = iter + 1
                    enddo
                enddo
                !$omp end parallel do
            case('cosine')
                norm_pcavecs = pcavecs
                !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,denom)
                do i = 1,self%N
                    denom = dsqrt(sum(real(pcavecs(:,i),dp)**2))
                    if( denom > DTINY ) norm_pcavecs(:,i) = pcavecs(:,i) / real(denom)
                enddo
                !$omp end parallel do
                norm_pcavecs_t = transpose(norm_pcavecs)
                if( trim(params_glob%kpca_target) .eq. 'ptcl' )then
                    !$omp parallel do default(shared) proc_bind(close) schedule(static) private(ind,ithr,iter,denom)
                    do ind = 1, self%N
                        ithr              = omp_get_thread_num() + 1
                        self%data(:,ind)  = pcavecs(:,ind)
                        norm_prev(:,ithr) = 1. / sqrt(real(self%D))
                        norm_data(:,ithr) = norm_pcavecs(:,ind)
                        iter              = 1
                        do while( abs(sum(norm_data(:,ithr) * norm_prev(:,ithr)) - 1.) > TOL .and. iter < its )
                            norm_prev(:,ithr) = norm_data(:,ithr)
                            ! 1. projecting each image on kernel space
                            proj_data(:,ithr) = matmul(norm_pcavecs_t, norm_prev(:,ithr)) * ker(:,ind)
                            ! 2. computing the pre-image
                            denom = sum(abs(real(proj_data(:,ithr),dp)))
                            if( denom < DTINY ) exit
                            self%data(:,ind)  = matmul(pcavecs, proj_data(:,ithr)) / real(denom)
                            denom             = dsqrt(sum(real(self%data(:,ind),dp)**2))
                            if( denom < DTINY ) exit
                            norm_data(:,ithr) = self%data(:,ind) / real(denom)
                            iter = iter + 1
                        enddo
                    enddo
                    !$omp end parallel do
                else
                    !$omp parallel do default(shared) proc_bind(close) schedule(static) private(ind,ithr,iter,denom)
                    do ind = 1, self%N
                        ithr              = omp_get_thread_num() + 1
                        self%data(:,ind)  = pcavecs(:,ind)
                        norm_prev(:,ithr) = 1. / sqrt(real(self%D))
                        norm_data(:,ithr) = norm_pcavecs(:,ind)
                        iter              = 1
                        do while( abs(sum(norm_data(:,ithr) * norm_prev(:,ithr)) - 1.) > TOL .and. iter < its )
                            norm_prev(:,ithr) = norm_data(:,ithr)
                            ! 1. projecting each image on kernel space
                            proj_data(:,ithr) = matmul(norm_pcavecs_t, norm_prev(:,ithr)) * ker(:,ind)
                            ! 2. computing the pre-image
                            self%data(:,ind)  = matmul(norm_pcavecs, proj_data(:,ithr))
                            denom             = dsqrt(sum(real(self%data(:,ind),dp)**2))
                            if( denom < DTINY ) exit
                            norm_data(:,ithr) = self%data(:,ind) / real(denom)
                            iter = iter + 1
                        enddo
                    enddo
                    !$omp end parallel do
                endif
        end select
        if( DEBUG )then
            call system_clock(end_time)
            print *, "Pre-imaging time: ", real(end_time-start_time)/real(rate), " seconds"
        endif
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
        integer  :: i, j
        real(dp) :: denom
        real     :: norm_mat(self%D,self%N)
        ! squared cosine similarity between pairs of rows
        norm_mat = mat
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,denom)
        do i = 1,self%N
            denom = dsqrt(sum(real(mat(:,i),dp)**2))
            if( denom > DTINY ) norm_mat(:,i) = mat(:,i) / real(denom)
        enddo
        !$omp end parallel do
        !$omp parallel do collapse(2) default(shared) proc_bind(close) schedule(static) private(i,j)
        do j = 1,self%N
            do i = 1,self%N
                ker(i,j) = sum(norm_mat(:,i) * norm_mat(:,j))
            enddo
        enddo
        !$omp end parallel do
        call self%kernel_center(ker)
    end subroutine cosine_kernel

    subroutine rbf_kernel( self, mat, ker )
        class(kpca_svd), intent(inout) :: self
        real,            intent(in)    :: mat(self%D,self%N)
        real,            intent(out)   :: ker(self%N,self%N)
        integer :: i, j
        ! squared euclidean distance between pairs of rows
        !$omp parallel do collapse(2) default(shared) proc_bind(close) schedule(static) private(i,j)
        do j = 1,self%N
            do i = 1,self%N
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
        real    :: eig_vals(self%Q), tmp_ker(self%N,self%N)
        integer :: i
        tmp_ker = ker
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i)
        do i = 1, self%N
            tmp_ker(:,i) = tmp_ker(:,i) - sum(tmp_ker(:,i))/real(self%N)
        enddo
        !$omp end parallel do
        ! computing eigvals/eigvecs
        call eigh(self%N, tmp_ker, self%Q, eig_vals, eig_vecs)
        eig_vals = eig_vals**2 / real(self%N)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i)
        do i = 1, self%Q
            eig_vecs(:,i) = eig_vecs(:,i) / sqrt(eig_vals(i))
        enddo
        !$omp end parallel do
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