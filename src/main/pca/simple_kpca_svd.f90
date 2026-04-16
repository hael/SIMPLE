!@descr: kPCA using 'Learning to Find Pre-Images', using svd for eigvals/eigvecs
module simple_kpca_svd
use simple_core_module_api
use simple_pca, only: pca
implicit none

public :: kpca_svd
private
#include "simple_local_flags.inc"

type, extends(pca) :: kpca_svd
    private
    real, allocatable :: E_zn(:,:)  !< expectations (feature vecs)
    real, allocatable :: data(:,:)  !< projected data on feature vecs
    integer           :: nthr                          !< number of threads
    character(len=16) :: kpca_backend = ''            !< backend ('exact' or 'nystrom')
    character(len=16) :: kpca_ker = ''                !< kernel type ('rbf' or 'cosine')
    character(len=16) :: kpca_target = ''             !< target type ('ptcl' or other)
    integer           :: kpca_nystrom_npts = 0        !< nr of Nystrom landmarks (0 => auto)
    logical           :: existence=.false.
    contains
    ! CONSTRUCTOR
    procedure :: new        => new_kpca
    ! SETTERS
    procedure :: set_params => set_params_kpca
    ! GETTERS
    procedure :: get_feat   => get_feat_kpca
    procedure :: generate   => generate_kpca
    ! CALCULATORS
    procedure :: master     => master_kpca
    ! DESTRUCTOR
    procedure :: kill       => kill_kpca
    ! PRIVATE
    procedure, private :: kernel_center
    procedure, private :: cosine_kernel
    procedure, private :: rbf_kernel
    procedure, private :: master_nystrom
    procedure, private :: compute_eigvecs
    procedure, private :: projected_kernel_col
    procedure, private :: select_nystrom_inds
end type

real, parameter :: C_CONST = 0.4  ! for rbf_kernel for testing

contains

    ! CONSTRUCTORS

    !>  \brief  is a constructor
    subroutine new_kpca( self, N, D, Q )
        class(kpca_svd), intent(inout) :: self
        integer,         intent(in)    :: N, D, Q
        call self%kill
        self%N           = N
        self%D           = D
        self%Q           = Q
        ! Initialize with defaults (use set_params() to override)
        self%nthr = 1
        self%kpca_backend = 'exact'
        self%kpca_ker = 'cosine'
        self%kpca_target = 'ptcl'
        self%kpca_nystrom_npts = 0
        ! allocate principal subspace and feature vectors
        allocate( self%E_zn(self%Q,self%N), self%data(self%D,self%N), source=0.)
        self%existence = .true.
    end subroutine new_kpca

    ! SETTERS

    !>  \brief  setter for runtime parameters
    subroutine set_params_kpca( self, nthr, kpca_ker, kpca_target, kpca_backend, kpca_nystrom_npts )
        class(kpca_svd),            intent(inout) :: self
        integer,          optional, intent(in)    :: nthr
        character(len=*), optional, intent(in)    :: kpca_ker
        character(len=*), optional, intent(in)    :: kpca_target
        character(len=*), optional, intent(in)    :: kpca_backend
        integer,          optional, intent(in)    :: kpca_nystrom_npts
        if( present(nthr) )then
            self%nthr = nthr
        endif
        if( present(kpca_ker) )then
            self%kpca_ker = trim(kpca_ker)
        endif
        if( present(kpca_target) )then
            self%kpca_target = trim(kpca_target)
        endif
        if( present(kpca_backend) )then
            self%kpca_backend = trim(kpca_backend)
        endif
        if( present(kpca_nystrom_npts) )then
            self%kpca_nystrom_npts = kpca_nystrom_npts
        endif
    end subroutine set_params_kpca

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
                   &proj_data(self%N,self%nthr), norm_prev(self%D,self%nthr), norm_data(self%D,self%nthr)
        integer  :: i, ind, iter, its, ithr
        if( trim(self%kpca_backend) .eq. 'nystrom' )then
            call self%master_nystrom(pcavecs, maxpcaits)
            return
        endif
        ! compute the kernel
        select case(trim(self%kpca_ker))
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
        select case(trim(self%kpca_ker))
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
                if( trim(self%kpca_target) .eq. 'ptcl' )then
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

    subroutine master_nystrom( self, pcavecs, maxpcaits )
        class(kpca_svd),   intent(inout) :: self
        real,              intent(in)    :: pcavecs(self%D,self%N)
        integer, optional, intent(in)    :: maxpcaits
        integer, parameter :: MAX_ITS = 500
        real,    parameter :: TOL     = 0.0001
        integer  :: m, q_used, r, its, ind, iter, ithr, i, j
        integer  :: landmark_inds(self%N)
        real(dp) :: denom, eig_tol
        real     :: eig_q(self%Q), alpha(self%N,self%Q), ker_col(self%N,self%nthr), proj_data(self%N,self%nthr),&
                   &norm_prev(self%D,self%nthr), norm_data(self%D,self%nthr), norm_pcavecs(self%D,self%N), norm_pcavecs_t(self%N,self%D)
        real, allocatable :: ker_nm(:,:), ker_mm(:,:), feat(:,:), feat_center(:), eig_w(:), eigvec_w(:,:), tmp_ker_mm(:,:),&
                            &tmp_gram(:,:), gram_eigvecs(:,:), landmark_mat(:,:), gram_small(:,:), gram_eigvecs_small(:,:)
        m = self%kpca_nystrom_npts
        if( m <= 0 ) m = min(self%N, max(4*self%Q, 32))
        m = min(self%N, max(self%Q, m))
        if( m < 1 )then
            THROW_HARD('kpca Nystrom backend requires at least one landmark')
        endif
        allocate(ker_nm(self%N,m), ker_mm(m,m), feat(self%N,m), feat_center(m), eig_w(m), eigvec_w(m,m), tmp_ker_mm(m,m),&
                 &tmp_gram(m,m), gram_eigvecs(m,self%Q), landmark_mat(self%D,m), source=0.)
        call self%select_nystrom_inds(m, landmark_inds(1:m))
        landmark_mat(:,1:m) = pcavecs(:,landmark_inds(1:m))
        select case(trim(self%kpca_ker))
            case('cosine')
                norm_pcavecs = pcavecs
                !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,denom)
                do i = 1,self%N
                    denom = dsqrt(sum(real(pcavecs(:,i),dp)**2))
                    if( denom > DTINY ) norm_pcavecs(:,i) = pcavecs(:,i) / real(denom)
                enddo
                !$omp end parallel do
                ker_nm(1:self%N,1:m) = matmul(transpose(norm_pcavecs), norm_pcavecs(:,landmark_inds(1:m)))
                ker_mm(1:m,1:m)      = matmul(transpose(norm_pcavecs(:,landmark_inds(1:m))), norm_pcavecs(:,landmark_inds(1:m)))
            case('rbf')
                !$omp parallel do default(shared) proc_bind(close) schedule(static) private(j,i)
                do j = 1,m
                    do i = 1,self%N
                        ker_nm(i,j) = euclid(pcavecs(:,i), landmark_mat(:,j))**2
                    enddo
                enddo
                !$omp end parallel do
                ker_nm(1:self%N,1:m) = exp(-ker_nm(1:self%N,1:m)/real(self%Q)/C_CONST)
                !$omp parallel do default(shared) proc_bind(close) schedule(static) private(j,i)
                do j = 1,m
                    ker_mm(j,j) = 1.
                    do i = 1,j-1
                        ker_mm(i,j) = euclid(landmark_mat(:,i), landmark_mat(:,j))**2
                        ker_mm(j,i) = ker_mm(i,j)
                    enddo
                enddo
                !$omp end parallel do
                do j = 1,m
                    do i = 1,j-1
                        ker_mm(i,j) = exp(-ker_mm(i,j)/real(self%Q)/C_CONST)
                        ker_mm(j,i) = ker_mm(i,j)
                    enddo
                enddo
            case default
                THROW_HARD('Unsupported kPCA kernel backend: '//trim(self%kpca_ker))
        end select
        tmp_ker_mm = ker_mm(1:m,1:m)
        call eigh(m, tmp_ker_mm, m, eig_w(1:m), eigvec_w)
        eig_w(1:m)    = eig_w(m:1:-1)
        eigvec_w(:,:) = eigvec_w(:,m:1:-1)
        eig_tol = max(real(DTINY,dp), 1.e-6_dp * max(real(maxval(eig_w(1:m)),dp), 1._dp))
        r = count(real(eig_w(1:m),dp) > eig_tol)
        if( r < 1 )then
            THROW_HARD('Nystrom kernel approximation is rank deficient')
        endif
        feat(1:self%N,1:r) = matmul(ker_nm(1:self%N,1:m), eigvec_w(:,1:r))
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i)
        do i = 1,r
            feat(:,i) = feat(:,i) / sqrt(max(eig_w(i), real(DTINY)))
        enddo
        !$omp end parallel do
        feat_center(1:r) = sum(feat(1:self%N,1:r), dim=1) / real(self%N)
        feat(1:self%N,1:r) = feat(1:self%N,1:r) - spread(feat_center(1:r), dim=1, ncopies=self%N)
        tmp_gram(1:r,1:r) = matmul(transpose(feat(1:self%N,1:r)), feat(1:self%N,1:r))
        q_used = min(self%Q, r)
        allocate(gram_small(r,r), gram_eigvecs_small(r,q_used))
        gram_small = tmp_gram(1:r,1:r)
        call eigh(r, gram_small, q_used, eig_q(1:q_used), gram_eigvecs_small)
        eig_q(1:q_used)                  = eig_q(q_used:1:-1)
        gram_eigvecs(1:r,1:q_used)       = gram_eigvecs_small(:,q_used:1:-1)
        alpha(:,1:q_used)                = matmul(feat(1:self%N,1:r), gram_eigvecs(1:r,1:q_used))
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i)
        do i = 1,q_used
            if( eig_q(i) > real(DTINY) ) alpha(:,i) = alpha(:,i) / sqrt(eig_q(i))
        enddo
        !$omp end parallel do
        if( q_used < self%Q )then
            alpha(:,q_used+1:self%Q) = 0.
            eig_q(q_used+1:self%Q)   = 0.
        endif
        its = MAX_ITS
        if( present(maxpcaits) ) its = maxpcaits
        select case(trim(self%kpca_ker))
            case('rbf')
                !$omp parallel do default(shared) proc_bind(close) schedule(static) private(ind,iter,ithr,i,denom)
                do ind = 1, self%N
                    ithr              = omp_get_thread_num() + 1
                    call self%projected_kernel_col(alpha, eig_q, q_used, ind, ker_col(:,ithr))
                    self%data(:,ind)  = pcavecs(:,ind)
                    norm_prev(:,ithr) = 0.
                    iter              = 1
                    do while( euclid(self%data(:,ind),norm_prev(:,ithr)) > TOL .and. iter < its )
                        norm_prev(:,ithr) = self%data(:,ind)
                        do i = 1,self%N
                            proj_data(i,ithr) = euclid(norm_prev(:,ithr), pcavecs(:,i))**2
                        enddo
                        proj_data(:,ithr) = exp(-proj_data(:,ithr)/real(self%Q)/C_CONST) * ker_col(:,ithr)
                        denom             = sum(real(proj_data(:,ithr),dp))
                        self%data(:,ind)  = matmul(pcavecs, proj_data(:,ithr))
                        if( denom > DTINY ) self%data(:,ind) = self%data(:,ind) / real(denom)
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
                if( trim(self%kpca_target) .eq. 'ptcl' )then
                    !$omp parallel do default(shared) proc_bind(close) schedule(static) private(ind,ithr,iter,denom)
                    do ind = 1, self%N
                        ithr              = omp_get_thread_num() + 1
                        call self%projected_kernel_col(alpha, eig_q, q_used, ind, ker_col(:,ithr))
                        self%data(:,ind)  = pcavecs(:,ind)
                        norm_prev(:,ithr) = 1. / sqrt(real(self%D))
                        norm_data(:,ithr) = norm_pcavecs(:,ind)
                        iter              = 1
                        do while( abs(sum(norm_data(:,ithr) * norm_prev(:,ithr)) - 1.) > TOL .and. iter < its )
                            norm_prev(:,ithr) = norm_data(:,ithr)
                            proj_data(:,ithr) = matmul(norm_pcavecs_t, norm_prev(:,ithr)) * ker_col(:,ithr)
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
                        call self%projected_kernel_col(alpha, eig_q, q_used, ind, ker_col(:,ithr))
                        self%data(:,ind)  = pcavecs(:,ind)
                        norm_prev(:,ithr) = 1. / sqrt(real(self%D))
                        norm_data(:,ithr) = norm_pcavecs(:,ind)
                        iter              = 1
                        do while( abs(sum(norm_data(:,ithr) * norm_prev(:,ithr)) - 1.) > TOL .and. iter < its )
                            norm_prev(:,ithr) = norm_data(:,ithr)
                            proj_data(:,ithr) = matmul(norm_pcavecs_t, norm_prev(:,ithr)) * ker_col(:,ithr)
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
        deallocate(ker_nm, ker_mm, feat, feat_center, eig_w, eigvec_w, tmp_ker_mm, tmp_gram, gram_eigvecs, landmark_mat, &
                   &gram_small, gram_eigvecs_small)
    end subroutine master_nystrom

    subroutine kernel_center( self, ker )
        class(kpca_svd), intent(inout) :: self
        real,            intent(inout) :: ker(self%N,self%N)
        real :: row_mean(self%N), col_mean(self%N), total_mean
        row_mean   = sum(ker, dim=2) / real(self%N)
        col_mean   = sum(ker, dim=1) / real(self%N)
        total_mean = sum(row_mean) / real(self%N)
        ! Appendix D.2.2 Centering in Feature Space from Schoelkopf, Bernhard, Support vector learning, 1997
        ker = ker - spread(row_mean, dim=2, ncopies=self%N) - spread(col_mean, dim=1, ncopies=self%N) + total_mean
    end subroutine kernel_center

    subroutine cosine_kernel( self, mat, ker )
        class(kpca_svd), intent(inout) :: self
        real,            intent(in)    :: mat(self%D,self%N)
        real,            intent(out)   :: ker(self%N,self%N)
        integer  :: i
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
        ker = matmul(transpose(norm_mat), norm_mat)
        call self%kernel_center(ker)
    end subroutine cosine_kernel

    subroutine rbf_kernel( self, mat, ker )
        class(kpca_svd), intent(inout) :: self
        real,            intent(in)    :: mat(self%D,self%N)
        real,            intent(out)   :: ker(self%N,self%N)
        integer :: i, j
        ! squared euclidean distance between pairs of rows
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,j)
        do j = 1,self%N
            ker(j,j) = 0.
            do i = 1,j-1
                ker(i,j) = euclid(mat(:,i), mat(:,j))**2
                ker(j,i) = ker(i,j)
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

    subroutine projected_kernel_col( self, eig_vecs, eig_vals, q_used, ind, ker_col )
        class(kpca_svd), intent(inout) :: self
        real,            intent(in)    :: eig_vecs(self%N,self%Q), eig_vals(self%Q)
        integer,         intent(in)    :: q_used, ind
        real,            intent(out)   :: ker_col(self%N)
        real :: weights(self%Q)
        weights = 0.
        if( q_used > 0 ) weights(1:q_used) = eig_vals(1:q_used) * eig_vecs(ind,1:q_used)
        ker_col = matmul(eig_vecs, weights)
    end subroutine projected_kernel_col

    subroutine select_nystrom_inds( self, m, inds )
        class(kpca_svd), intent(inout) :: self
        integer,         intent(in)    :: m
        integer,         intent(out)   :: inds(m)
        integer :: i
        if( m == 1 )then
            inds(1) = 1
            return
        endif
        do i = 1,m
            inds(i) = 1 + ((i-1) * (self%N-1)) / (m-1)
        enddo
    end subroutine select_nystrom_inds

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
