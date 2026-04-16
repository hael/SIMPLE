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
    real              :: kpca_rbf_gamma = 0.          !< RBF gamma (0 => auto)
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
    procedure, private :: get_rbf_gamma
    procedure, private :: projected_kernel_col
    procedure, private :: select_nystrom_inds
    procedure, private :: dense_tmm
    procedure, private :: dense_mm
    procedure, private :: gram_symmetric
    procedure, private :: center_columns
    procedure, private :: partial_eigh_sym
    procedure, private :: orthonormalize_cols
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
        self%kpca_rbf_gamma = 0.
        ! allocate principal subspace and feature vectors
        allocate( self%E_zn(self%Q,self%N), self%data(self%D,self%N), source=0.)
        self%existence = .true.
    end subroutine new_kpca

    ! SETTERS

    !>  \brief  setter for runtime parameters
    subroutine set_params_kpca( self, nthr, kpca_ker, kpca_target, kpca_backend, kpca_nystrom_npts, kpca_rbf_gamma )
        class(kpca_svd),            intent(inout) :: self
        integer,          optional, intent(in)    :: nthr
        character(len=*), optional, intent(in)    :: kpca_ker
        character(len=*), optional, intent(in)    :: kpca_target
        character(len=*), optional, intent(in)    :: kpca_backend
        integer,          optional, intent(in)    :: kpca_nystrom_npts
        real,             optional, intent(in)    :: kpca_rbf_gamma
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
        if( present(kpca_rbf_gamma) )then
            self%kpca_rbf_gamma = kpca_rbf_gamma
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
        real(dp) :: denom, rbf_gamma
        real     :: ker(self%N,self%N), eig_vecs(self%N,self%Q), eig_vals(self%Q), norm_pcavecs(self%D,self%N),&
                   &norm_pcavecs_t(self%N,self%D)
        real, allocatable :: ker_col(:,:), proj_data(:,:), norm_prev(:,:), norm_data(:,:)
        integer  :: i, ind, iter, its, ithr, nworkthr
        if( trim(self%kpca_backend) .eq. 'nystrom' )then
            call self%master_nystrom(pcavecs, maxpcaits)
            return
        endif
        nworkthr = max(1, max(self%nthr, omp_get_max_threads()))
        allocate(ker_col(self%N,nworkthr), proj_data(self%N,nworkthr), norm_prev(self%D,nworkthr), norm_data(self%D,nworkthr), source=0.)
        rbf_gamma = 0._dp
        if( trim(self%kpca_ker) .eq. 'rbf' ) rbf_gamma = self%get_rbf_gamma(pcavecs)
        ! compute the kernel
        select case(trim(self%kpca_ker))
            case('rbf')
                call self%rbf_kernel(   pcavecs, ker)
            case('cosine')
                call self%cosine_kernel(pcavecs, ker)
        end select
        ! compute the sorted principle components of the kernel above
        if( DEBUG ) call system_clock(start_time, rate)
        call self%compute_eigvecs(ker, eig_vecs, eig_vals)
        self%E_zn = transpose(matmul(ker, eig_vecs))
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
        select case(trim(self%kpca_ker))
            case('rbf')
                !$omp parallel do default(shared) proc_bind(close) schedule(static) private(ind,iter,ithr,i,denom)
                do ind = 1, self%N
                    ithr              = omp_get_thread_num() + 1
                    call self%projected_kernel_col(eig_vecs, eig_vals, self%Q, ind, ker_col(:,ithr))
                    self%data(:,ind)  = pcavecs(:,ind)
                    norm_prev(:,ithr) = 0.
                    iter              = 1
                    do while( euclid(self%data(:,ind),norm_prev(:,ithr)) > TOL .and. iter < its )
                        norm_prev(:,ithr) = self%data(:,ind)
                        ! 1. projecting each image on kernel space
                        do i = 1,self%N
                            proj_data(i,ithr) = euclid(norm_prev(:,ithr), pcavecs(:,i))**2
                        enddo
                        ! 2. applying the principle components to the projected vector
                        proj_data(:,ithr) = max(0., ker_col(:,ithr)) * exp(-real(rbf_gamma) * proj_data(:,ithr))
                        ! 3. computing the pre-image (of the image in step 1) using the result in step 2
                        denom             = sum(real(proj_data(:,ithr),dp))
                        if( denom < DTINY ) exit
                        self%data(:,ind)  = matmul(pcavecs, proj_data(:,ithr)) / real(denom)
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
                        call self%projected_kernel_col(eig_vecs, eig_vals, self%Q, ind, ker_col(:,ithr))
                        self%data(:,ind)  = pcavecs(:,ind)
                        norm_prev(:,ithr) = 1. / sqrt(real(self%D))
                        norm_data(:,ithr) = norm_pcavecs(:,ind)
                        iter              = 1
                        do while( abs(sum(norm_data(:,ithr) * norm_prev(:,ithr)) - 1.) > TOL .and. iter < its )
                            norm_prev(:,ithr) = norm_data(:,ithr)
                            ! 1. projecting each image on kernel space
                            proj_data(:,ithr) = matmul(norm_pcavecs_t, norm_prev(:,ithr)) * ker_col(:,ithr)
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
                        call self%projected_kernel_col(eig_vecs, eig_vals, self%Q, ind, ker_col(:,ithr))
                        self%data(:,ind)  = pcavecs(:,ind)
                        norm_prev(:,ithr) = 1. / sqrt(real(self%D))
                        norm_data(:,ithr) = norm_pcavecs(:,ind)
                        iter              = 1
                        do while( abs(sum(norm_data(:,ithr) * norm_prev(:,ithr)) - 1.) > TOL .and. iter < its )
                            norm_prev(:,ithr) = norm_data(:,ithr)
                            ! 1. projecting each image on kernel space
                            proj_data(:,ithr) = matmul(norm_pcavecs_t, norm_prev(:,ithr)) * ker_col(:,ithr)
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
        deallocate(ker_col, proj_data, norm_prev, norm_data)
    end subroutine master_kpca

    subroutine master_nystrom( self, pcavecs, maxpcaits )
        class(kpca_svd),   intent(inout) :: self
        real,              intent(in)    :: pcavecs(self%D,self%N)
        integer, optional, intent(in)    :: maxpcaits
        integer, parameter :: MAX_ITS = 500
        real,    parameter :: TOL     = 0.0001
        logical, parameter :: PROFILE = .true.
        integer  :: m, q_used, r, r_keep, its, ind, iter, ithr, i, j, nworkthr
        integer  :: landmark_inds(self%N)
        real(dp) :: denom, eig_tol, rbf_gamma
        integer(int64) :: t0, t1
        real(real64)   :: trate
        real     :: eig_q(self%Q), alpha(self%N,self%Q), norm_pcavecs(self%D,self%N), norm_pcavecs_t(self%N,self%D)
        real, allocatable :: ker_col(:,:), proj_data(:,:), norm_prev(:,:), norm_data(:,:)
        real, allocatable :: ker_nm(:,:), ker_mm(:,:), feat(:,:), feat_center(:), eig_w(:), eigvec_w(:,:), tmp_ker_mm(:,:),&
                            &tmp_gram(:,:), gram_eigvecs(:,:), landmark_mat(:,:), gram_small(:,:), gram_eigvecs_small(:,:)
        m = self%kpca_nystrom_npts
        if( m <= 0 ) m = min(self%N, max(4*self%Q, 32))
        m = min(self%N, max(self%Q, m))
        if( m < 1 )then
            THROW_HARD('kpca Nystrom backend requires at least one landmark')
        endif
        if( PROFILE ) call system_clock(t0, trate)
        nworkthr = max(1, max(self%nthr, omp_get_max_threads()))
        rbf_gamma = 0._dp
        if( trim(self%kpca_ker) .eq. 'rbf' ) rbf_gamma = self%get_rbf_gamma(pcavecs)
        r_keep = min(m, max(self%Q, min(2*self%Q, self%Q + 32)))
        allocate(ker_nm(self%N,m), ker_mm(m,m), feat(self%N,m), feat_center(m), eig_w(m), eigvec_w(m,m), tmp_ker_mm(m,m),&
                 &tmp_gram(m,m), gram_eigvecs(m,self%Q), landmark_mat(self%D,m), ker_col(self%N,nworkthr), proj_data(self%N,nworkthr),&
                 &norm_prev(self%D,nworkthr), norm_data(self%D,nworkthr), source=0.)
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
                call self%dense_tmm(norm_pcavecs, norm_pcavecs(:,landmark_inds(1:m)), ker_nm(1:self%N,1:m))
                call self%dense_tmm(norm_pcavecs(:,landmark_inds(1:m)), norm_pcavecs(:,landmark_inds(1:m)), ker_mm(1:m,1:m))
            case('rbf')
                !$omp parallel do default(shared) proc_bind(close) schedule(static) private(j,i)
                do j = 1,m
                    do i = 1,self%N
                        ker_nm(i,j) = euclid(pcavecs(:,i), landmark_mat(:,j))**2
                    enddo
                enddo
                !$omp end parallel do
                ker_nm(1:self%N,1:m) = exp(-real(rbf_gamma) * ker_nm(1:self%N,1:m))
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
                        ker_mm(i,j) = exp(-real(rbf_gamma) * ker_mm(i,j))
                        ker_mm(j,i) = ker_mm(i,j)
                    enddo
                enddo
            case default
                THROW_HARD('Unsupported kPCA kernel backend: '//trim(self%kpca_ker))
        end select
        if( PROFILE )then
            call system_clock(t1)
            write(logfhandle,'(A,F8.3,A,I8,A,I8)') 'kPCA Nyström kernel/features setup: ', real(t1-t0)/real(trate), ' s; N=', self%N, ' m=', m
            call system_clock(t0)
        endif
        tmp_ker_mm = ker_mm(1:m,1:m)
        call self%partial_eigh_sym(tmp_ker_mm, r_keep, eig_w(1:r_keep), eigvec_w(:,1:r_keep))
        if( PROFILE )then
            call system_clock(t1)
            write(logfhandle,'(A,F8.3,A,I8)') 'kPCA Nyström landmark eigensolve: ', real(t1-t0)/real(trate), ' s; r_keep=', r_keep
            call system_clock(t0)
        endif
        eig_tol = max(real(DTINY,dp), 1.e-6_dp * max(real(maxval(eig_w(1:r_keep)),dp), 1._dp))
        r = count(real(eig_w(1:r_keep),dp) > eig_tol)
        if( r < 1 )then
            THROW_HARD('Nystrom kernel approximation is rank deficient')
        endif
        call self%dense_mm(ker_nm(1:self%N,1:m), eigvec_w(:,1:r), feat(1:self%N,1:r))
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i)
        do i = 1,r
            feat(:,i) = feat(:,i) / sqrt(max(eig_w(i), real(DTINY)))
        enddo
        !$omp end parallel do
        call self%center_columns(feat(1:self%N,1:r), feat_center(1:r))
        call self%gram_symmetric(feat(1:self%N,1:r), tmp_gram(1:r,1:r))
        q_used = min(self%Q, r)
        allocate(gram_small(r,r), gram_eigvecs_small(r,q_used))
        gram_small = tmp_gram(1:r,1:r)
        if( PROFILE )then
            call system_clock(t1)
            write(logfhandle,'(A,F8.3,A,I8)') 'kPCA Nyström feature/gram build: ', real(t1-t0)/real(trate), ' s; r=', r
            call system_clock(t0)
        endif
        call self%partial_eigh_sym(gram_small, q_used, eig_q(1:q_used), gram_eigvecs_small)
        gram_eigvecs(1:r,1:q_used)       = gram_eigvecs_small(:,1:q_used)
        call self%dense_mm(feat(1:self%N,1:r), gram_eigvecs(1:r,1:q_used), alpha(:,1:q_used))
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i)
        do i = 1,q_used
            if( eig_q(i) > real(DTINY) ) alpha(:,i) = alpha(:,i) / sqrt(eig_q(i))
        enddo
        !$omp end parallel do
        if( PROFILE )then
            call system_clock(t1)
            write(logfhandle,'(A,F8.3,A,I8)') 'kPCA Nyström reduced eigensolve/proj: ', real(t1-t0)/real(trate), ' s; q=', q_used
            call system_clock(t0)
        endif
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
                        proj_data(:,ithr) = max(0., ker_col(:,ithr)) * exp(-real(rbf_gamma) * proj_data(:,ithr))
                        denom             = sum(real(proj_data(:,ithr),dp))
                        if( denom < DTINY ) exit
                        self%data(:,ind)  = matmul(pcavecs, proj_data(:,ithr)) / real(denom)
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
        if( PROFILE )then
            call system_clock(t1)
            write(logfhandle,'(A,F8.3,A)') 'kPCA Nyström pre-image/reconstruct: ', real(t1-t0)/real(trate), ' s'
        endif
        deallocate(ker_nm, ker_mm, feat, feat_center, eig_w, eigvec_w, tmp_ker_mm, tmp_gram, gram_eigvecs, landmark_mat, &
                   &gram_small, gram_eigvecs_small, ker_col, proj_data, norm_prev, norm_data)
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
        integer  :: i, j
        real(dp) :: gamma
        gamma = self%get_rbf_gamma(mat)
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
        ker = exp(-real(gamma) * ker)
        call self%kernel_center(ker)
    end subroutine rbf_kernel

    subroutine compute_eigvecs( self, ker, eig_vecs, eig_vals )
        class(kpca_svd), intent(inout) :: self
        real,            intent(in)    :: ker(self%N,self%N)
        real,            intent(out)   :: eig_vals(self%Q)
        real,            intent(inout) :: eig_vecs(self%N,self%Q)
        real    :: tmp_ker(self%N,self%N)
        integer :: i
        tmp_ker = ker
        ! computing eigvals/eigvecs
        call eigh(self%N, tmp_ker, self%Q, eig_vals, eig_vecs)
        eig_vals = eig_vals(self%Q:1:-1)
        eig_vecs = eig_vecs(:,self%Q:1:-1)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i)
        do i = 1, self%Q
            if( eig_vals(i) > real(DTINY) )then
                eig_vecs(:,i) = eig_vecs(:,i) / sqrt(eig_vals(i))
            else
                eig_vecs(:,i) = 0.
                eig_vals(i)   = 0.
            endif
        enddo
        !$omp end parallel do
    end subroutine compute_eigvecs

    real(dp) function get_rbf_gamma( self, mat ) result( gamma )
        class(kpca_svd), intent(in) :: self
        real,            intent(in) :: mat(self%D,self%N)
        integer  :: i, j, npairs
        real(dp) :: mean_sqdist
        if( self%kpca_rbf_gamma > 0. )then
            gamma = real(self%kpca_rbf_gamma, dp)
            return
        endif
        mean_sqdist = 0._dp
        npairs      = 0
        do j = 1,self%N
            do i = 1,j-1
                mean_sqdist = mean_sqdist + euclid(real(mat(:,i),dp), real(mat(:,j),dp))**2
                npairs = npairs + 1
            enddo
        enddo
        if( npairs < 1 )then
            gamma = 1._dp
        else
            mean_sqdist = mean_sqdist / real(npairs, dp)
            gamma = 1._dp / max(mean_sqdist, DTINY)
        endif
    end function get_rbf_gamma

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

    subroutine dense_tmm( self, left, right, out )
        class(kpca_svd), intent(in)  :: self
        real,            intent(in)  :: left(:,:), right(:,:)
        real,            intent(out) :: out(:,:)
        integer  :: i, j, k, nrow, ninner, ncol
        real(dp) :: acc
        out = 0.
        nrow   = size(left, 2)
        ninner = size(left, 1)
        ncol   = size(right, 2)
        !$omp parallel do collapse(2) default(shared) proc_bind(close) schedule(static) private(i,j,k,acc)
        do j = 1, ncol
            do i = 1, nrow
                acc = 0._dp
                do k = 1, ninner
                    acc = acc + real(left(k,i),dp) * real(right(k,j),dp)
                enddo
                out(i,j) = real(acc)
            enddo
        enddo
        !$omp end parallel do
    end subroutine dense_tmm

    subroutine dense_mm( self, left, right, out )
        class(kpca_svd), intent(in)  :: self
        real,            intent(in)  :: left(:,:), right(:,:)
        real,            intent(out) :: out(:,:)
        integer  :: i, j, k, nrow, ninner, ncol
        real(dp) :: acc
        out = 0.
        nrow   = size(left, 1)
        ninner = size(left, 2)
        ncol   = size(right, 2)
        !$omp parallel do collapse(2) default(shared) proc_bind(close) schedule(static) private(i,j,k,acc)
        do j = 1, ncol
            do i = 1, nrow
                acc = 0._dp
                do k = 1, ninner
                    acc = acc + real(left(i,k),dp) * real(right(k,j),dp)
                enddo
                out(i,j) = real(acc)
            enddo
        enddo
        !$omp end parallel do
    end subroutine dense_mm

    subroutine gram_symmetric( self, feat, gram )
        class(kpca_svd), intent(in)  :: self
        real,            intent(in)  :: feat(:,:)
        real,            intent(out) :: gram(:,:)
        integer  :: i, j, k, nrow, ncol
        real(dp) :: acc
        nrow = size(feat, 1)
        ncol = size(feat, 2)
        gram = 0.
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,j,k,acc)
        do j = 1, ncol
            do i = 1, j
                acc = 0._dp
                do k = 1, nrow
                    acc = acc + real(feat(k,i),dp) * real(feat(k,j),dp)
                enddo
                gram(i,j) = real(acc)
                gram(j,i) = gram(i,j)
            enddo
        enddo
        !$omp end parallel do
    end subroutine gram_symmetric

    subroutine center_columns( self, mat, col_mean )
        class(kpca_svd), intent(in)    :: self
        real,            intent(inout) :: mat(:,:)
        real,            intent(out)   :: col_mean(:)
        integer  :: i, j, nrow, ncol
        real(dp) :: acc
        nrow = size(mat, 1)
        ncol = size(mat, 2)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,j,acc)
        do j = 1, ncol
            acc = 0._dp
            do i = 1, nrow
                acc = acc + real(mat(i,j),dp)
            enddo
            col_mean(j) = real(acc / real(nrow, dp))
            do i = 1, nrow
                mat(i,j) = mat(i,j) - col_mean(j)
            enddo
        enddo
        !$omp end parallel do
    end subroutine center_columns

    subroutine partial_eigh_sym( self, mat, neigs, eigvals, eigvecs )
        class(kpca_svd), intent(in)    :: self
        real,            intent(in)    :: mat(:,:)
        integer,         intent(in)    :: neigs
        real,            intent(out)   :: eigvals(neigs)
        real,            intent(out)   :: eigvecs(size(mat,1),neigs)
        integer, parameter :: MAX_IT = 10
        integer :: n, ksub, iter, j
        real, allocatable :: basis(:,:), z(:,:), small(:,:), small_vecs(:,:)
        n = size(mat,1)
        if( neigs >= n .or. n <= max(64, 2*neigs) )then
            allocate(small(n,n))
            small = mat
            call eigh(n, small, neigs, eigvals, eigvecs)
            eigvals = eigvals(neigs:1:-1)
            eigvecs = eigvecs(:,neigs:1:-1)
            deallocate(small)
            return
        endif
        ksub = min(n, max(neigs + 8, min(2*neigs, neigs + 32)))
        allocate(basis(n,ksub), z(n,ksub), small(ksub,ksub), small_vecs(ksub,neigs), source=0.)
        do j = 1, ksub
            basis(:,j) = mat(:,j)
            if( sum(abs(real(basis(:,j),dp))) <= DTINY ) basis(j,j) = 1.
        enddo
        call self%orthonormalize_cols(basis)
        do iter = 1, MAX_IT
            call self%dense_mm(mat, basis, z)
            call self%orthonormalize_cols(z)
            basis = z
        enddo
        call self%dense_mm(mat, basis, z)
        call self%dense_tmm(basis, z, small)
        call eigh(ksub, small, neigs, eigvals, small_vecs)
        eigvals    = eigvals(neigs:1:-1)
        small_vecs = small_vecs(:,neigs:1:-1)
        call self%dense_mm(basis, small_vecs, eigvecs)
        deallocate(basis, z, small, small_vecs)
    end subroutine partial_eigh_sym

    subroutine orthonormalize_cols( self, mat )
        class(kpca_svd), intent(in)    :: self
        real,            intent(inout) :: mat(:,:)
        integer  :: i, j, nrow, ncol
        real(dp) :: proj, nrm
        nrow = size(mat,1)
        ncol = size(mat,2)
        do j = 1, ncol
            do i = 1, j-1
                proj = sum(real(mat(:,i),dp) * real(mat(:,j),dp))
                mat(:,j) = mat(:,j) - real(proj) * mat(:,i)
            enddo
            nrm = sqrt(sum(real(mat(:,j),dp) * real(mat(:,j),dp)))
            if( nrm <= DTINY )then
                mat(:,j) = 0.
                mat(1 + mod(j-1, nrow), j) = 1.
            else
                mat(:,j) = mat(:,j) / real(nrm)
            endif
        enddo
    end subroutine orthonormalize_cols

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
