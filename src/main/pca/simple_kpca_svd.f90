!@descr: kPCA using 'Learning to Find Pre-Images', using svd for eigvals/eigvecs
module simple_kpca_svd
use simple_core_module_api
use simple_pca, only: pca
implicit none

public :: kpca_svd, suggest_kpca_nystrom_neigs
private
#include "simple_local_flags.inc"

integer, parameter :: NYSTROM_AUTO_NPTS = 128

type, extends(pca) :: kpca_svd
    private
    real, allocatable :: E_zn(:,:)                   !< expectations (feature vecs)
    real, allocatable :: data(:,:)                   !< projected data on feature vecs
    integer           :: nthr                        !< number of threads
    character(len=16) :: kpca_backend      = ''      !< backend ('exact' or 'nystrom')
    character(len=16) :: kpca_ker          = ''      !< kernel type ('rbf' or 'cosine')
    character(len=16) :: kpca_target       = ''      !< target type ('ptcl' or other)
    integer           :: kpca_nystrom_npts = 512       !< nr of Nystrom landmarks
    integer           :: kpca_nystrom_local_nbrs = 128 !< max extra local support for Nyström reconstruction
    real              :: kpca_rbf_gamma    = 0.      !< RBF gamma (0 => auto)
    logical           :: existence         = .false.
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
    procedure, private :: select_local_support_inds
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
        ! Initialize with defaults (use set_params() to override)
        self%nthr              = 1
        self%kpca_backend      = 'nystrom'
        self%kpca_ker          = 'cosine'
        self%kpca_target       = 'ptcl'
        self%kpca_nystrom_npts = 512
        self%kpca_nystrom_local_nbrs = 128
        self%kpca_rbf_gamma    = 0.
        ! allocate principal subspace and feature vectors
        allocate( self%E_zn(self%Q,self%N), self%data(self%D,self%N), source=0.)
        self%existence = .true.
    end subroutine new_kpca

    ! SETTERS

    !>  \brief  setter for runtime parameters
    subroutine set_params_kpca( self, nthr, kpca_ker, kpca_target, kpca_backend, kpca_nystrom_npts, kpca_rbf_gamma, kpca_nystrom_local_nbrs )
        class(kpca_svd),            intent(inout) :: self
        integer,          optional, intent(in)    :: nthr
        character(len=*), optional, intent(in)    :: kpca_ker
        character(len=*), optional, intent(in)    :: kpca_target
        character(len=*), optional, intent(in)    :: kpca_backend
        integer,          optional, intent(in)    :: kpca_nystrom_npts
        real,             optional, intent(in)    :: kpca_rbf_gamma
        integer,          optional, intent(in)    :: kpca_nystrom_local_nbrs
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
        if( present(kpca_nystrom_local_nbrs) )then
            self%kpca_nystrom_local_nbrs = kpca_nystrom_local_nbrs
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
        logical, parameter :: DEBUG   = .true.
        integer(int64)     :: start_time, end_time
        real(real64)       :: rate
        real(dp) :: denom, rbf_gamma
        integer  :: nthr_use
        real     :: ker(self%N,self%N), eig_vecs(self%N,self%Q), eig_vals(self%Q)
        real, allocatable :: ker_col(:,:), proj_data(:,:), norm_prev(:,:), norm_data(:,:)
        real, allocatable :: norm_pcavecs(:,:), norm_pcavecs_t(:,:)
        integer  :: i, ind, iter, its, ithr
        if( DEBUG )then
            write(logfhandle,'(A,A,A,A,A,A,A,I8,A,I8,A,I8)') 'kPCA master entry: backend=', trim(self%kpca_backend), &
                '; kernel=', trim(self%kpca_ker), '; target=', trim(self%kpca_target), '; N=', self%N, ' D=', self%D, ' Q=', self%Q
            call flush(logfhandle)
        endif
        if( trim(self%kpca_backend) .eq. 'nystrom' )then
            if( DEBUG )then
                write(logfhandle,'(A)') 'kPCA master dispatch: entering Nyström backend'
                call flush(logfhandle)
            endif
            call self%master_nystrom(pcavecs, maxpcaits)
            return
        endif
        nthr_use = max(1, max(self%nthr, omp_get_max_threads()))
        allocate(ker_col(self%N,nthr_use), proj_data(self%N,nthr_use), norm_prev(self%D,nthr_use), norm_data(self%D,nthr_use), source=0.)
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
                allocate(norm_pcavecs(self%D,self%N), source=pcavecs)
                norm_pcavecs = pcavecs
                !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,denom)
                do i = 1,self%N
                    denom = dsqrt(sum(real(pcavecs(:,i),dp)**2))
                    if( denom > DTINY ) norm_pcavecs(:,i) = pcavecs(:,i) / real(denom)
                enddo
                !$omp end parallel do
                allocate(norm_pcavecs_t(self%N,self%D))
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
        if( allocated(norm_pcavecs)   ) deallocate(norm_pcavecs)
        if( allocated(norm_pcavecs_t) ) deallocate(norm_pcavecs_t)
        deallocate(ker_col, proj_data, norm_prev, norm_data)
    end subroutine master_kpca

    subroutine master_nystrom( self, pcavecs, maxpcaits )
        class(kpca_svd),   intent(inout) :: self
        real,              intent(in)    :: pcavecs(self%D,self%N)
        integer, optional, intent(in)    :: maxpcaits
        integer, parameter :: MAX_ITS = 500
        integer, parameter :: EARLY_STOP_PATIENCE = 3
        real,    parameter :: TOL     = 0.0001
        real(dp), parameter :: EARLY_STOP_REL = 1.e-5_dp
        logical, parameter :: PROFILE = .true.
        integer  :: m, q_used, r, r_keep, its, ind, iter, ithr, i, j, k, nthr_use, local_nbrs
        integer  :: stagn_count
        integer  :: landmark_inds(self%N)
        integer  :: support_count
        real(dp) :: denom, eig_tol, rbf_gamma, err_prev, err_curr, score
        integer(int64) :: t0, t1
        real(real64)   :: trate
        real, allocatable :: eig_q(:), alpha(:,:)
        real, allocatable :: ker_col(:,:), proj_data(:,:), norm_prev(:,:), norm_data(:,:)
        real, allocatable :: norm_pcavecs(:,:), norm_landmarks(:,:), norm_landmarks_t(:,:)
        real, allocatable :: ker_nm(:,:), ker_mm(:,:), feat(:,:), feat_center(:), eig_w(:), eigvec_w(:,:), tmp_ker_mm(:,:),&
                            &tmp_gram(:,:), gram_eigvecs(:,:), landmark_mat(:,:), gram_small(:,:), gram_eigvecs_small(:,:),&
                            &landmark_eigvals(:), landmark_eigvecs(:,:), gram_eigvals(:)
        real, allocatable :: support_scores(:,:)
        integer, allocatable :: support_inds(:,:), support_counts(:)
        logical, allocatable :: is_landmark(:)
        m = resolve_nystrom_npts(self%N, self%Q, self%kpca_nystrom_npts)
        if( m < 1 )then
            THROW_HARD('kpca Nystrom backend requires at least one landmark')
        endif
        if( PROFILE )then
            call system_clock(t0, trate)
            write(logfhandle,'(A)') 'kPCA Nyström entered'
            call flush(logfhandle)
        endif
        if( PROFILE )then
            write(logfhandle,'(A,A,A,I8,A,I8,A,I8,A,I8)') 'kPCA Nyström start: kernel=', trim(self%kpca_ker), &
                '; N=', self%N, ' D=', self%D, ' Q=', self%Q, ' m=', m
            call flush(logfhandle)
        endif
        rbf_gamma = 0._dp
        if( trim(self%kpca_ker) .eq. 'rbf' ) rbf_gamma = self%get_rbf_gamma(pcavecs)
        r_keep = min(m, max(self%Q, min(2*self%Q, self%Q + 32)))
        if( PROFILE ) call system_clock(t0)
        allocate(eig_q(self%Q), alpha(self%N,self%Q), source=0.)
        if( PROFILE )then
            call system_clock(t1)
            write(logfhandle,'(A,F8.3,A,I8,A,I8)') 'kPCA Nyström alpha/eig alloc: ', real(t1-t0)/real(trate), ' s; N=', self%N, ' Q=', self%Q
            call flush(logfhandle)
        endif
        if( PROFILE ) call system_clock(t0)
        nthr_use = max(1, max(self%nthr, omp_get_max_threads()))
        local_nbrs = 0
        local_nbrs = min(max(0, self%kpca_nystrom_local_nbrs), max(0, self%N - m))
        allocate(ker_nm(self%N,m), ker_mm(m,m), feat(self%N,m), feat_center(m), eig_w(m), eigvec_w(m,m), tmp_ker_mm(m,m),&
                 &tmp_gram(m,m), gram_eigvecs(m,self%Q), landmark_mat(self%D,m), ker_col(self%N,nthr_use), proj_data(m,nthr_use),&
                 &norm_prev(self%D,nthr_use), norm_data(self%D,nthr_use), source=0.)
        if( local_nbrs > 0 )then
            allocate(support_scores(local_nbrs,nthr_use), support_inds(local_nbrs,nthr_use), support_counts(nthr_use), is_landmark(self%N))
            support_scores = 0.
            support_inds   = 0
            support_counts = 0
            is_landmark = .false.
        endif
        if( PROFILE )then
            call system_clock(t1)
            write(logfhandle,'(A,F8.3,A)') 'kPCA Nyström work alloc: ', real(t1-t0)/real(trate), ' s'
            call flush(logfhandle)
            call system_clock(t0)
        endif
        call self%select_nystrom_inds(m, landmark_inds(1:m), pcavecs)
        if( local_nbrs > 0 ) is_landmark(landmark_inds(1:m)) = .true.
        if( PROFILE )then
            call system_clock(t1)
            write(logfhandle,'(A,F8.3,A)') 'kPCA Nyström landmark index selection: ', real(t1-t0)/real(trate), ' s'
            call flush(logfhandle)
            call system_clock(t0)
        endif
        landmark_mat(:,1:m) = pcavecs(:,landmark_inds(1:m))
        if( PROFILE )then
            call system_clock(t1)
            write(logfhandle,'(A,F8.3,A)') 'kPCA Nyström landmark copy: ', real(t1-t0)/real(trate), ' s'
            call flush(logfhandle)
            call system_clock(t0)
        endif
        if( PROFILE )then
            call system_clock(t1)
            write(logfhandle,'(A,F8.3,A,I8,A,I8,A,I8)') 'kPCA Nyström alloc/init: ', real(t1-t0)/real(trate), ' s; N=', self%N, ' D=', self%D, ' m=', m
            call flush(logfhandle)
            call system_clock(t0)
        endif
        select case(trim(self%kpca_ker))
            case('cosine')
                allocate(norm_pcavecs(self%D,self%N), source=pcavecs)
                norm_pcavecs = pcavecs
                !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,denom)
                do i = 1,self%N
                    denom = dsqrt(sum(real(pcavecs(:,i),dp)**2))
                    if( denom > DTINY ) norm_pcavecs(:,i) = pcavecs(:,i) / real(denom)
                enddo
                !$omp end parallel do
                if( PROFILE )then
                    call system_clock(t1)
                    write(logfhandle,'(A,F8.3,A)') 'kPCA Nyström cosine normalize: ', real(t1-t0)/real(trate), ' s'
                    call flush(logfhandle)
                    call system_clock(t0)
                endif
                allocate(norm_landmarks(self%D,m), norm_landmarks_t(m,self%D), source=0.)
                norm_landmarks(:,1:m) = norm_pcavecs(:,landmark_inds(1:m))
                norm_landmarks_t      = transpose(norm_landmarks)
                call self%dense_tmm(norm_pcavecs, norm_pcavecs(:,landmark_inds(1:m)), ker_nm(1:self%N,1:m))
                if( PROFILE )then
                    call system_clock(t1)
                    write(logfhandle,'(A,F8.3,A,I8,A,I8)') 'kPCA Nyström cosine K_nm: ', real(t1-t0)/real(trate), ' s; N=', self%N, ' m=', m
                    call flush(logfhandle)
                    call system_clock(t0)
                endif
                call self%dense_tmm(norm_pcavecs(:,landmark_inds(1:m)), norm_pcavecs(:,landmark_inds(1:m)), ker_mm(1:m,1:m))
                if( PROFILE )then
                    call system_clock(t1)
                    write(logfhandle,'(A,F8.3,A,I8)') 'kPCA Nyström cosine K_mm: ', real(t1-t0)/real(trate), ' s; m=', m
                    call flush(logfhandle)
                    call system_clock(t0)
                endif
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
            write(logfhandle,'(A,F8.3,A,I8,A,I8)') 'kPCA Nyström kernel setup: ', real(t1-t0)/real(trate), ' s; N=', self%N, ' m=', m
            call flush(logfhandle)
            call system_clock(t0)
        endif
        allocate(landmark_eigvals(r_keep), landmark_eigvecs(m,r_keep))
        tmp_ker_mm = ker_mm(1:m,1:m)
        call self%partial_eigh_sym(tmp_ker_mm, r_keep, landmark_eigvals, landmark_eigvecs)
        eig_w(1:r_keep)        = landmark_eigvals
        eigvec_w(:,1:r_keep)   = landmark_eigvecs
        deallocate(landmark_eigvals, landmark_eigvecs)
        if( PROFILE )then
            call system_clock(t1)
            write(logfhandle,'(A,F8.3,A,I8)') 'kPCA Nyström landmark eigensolve: ', real(t1-t0)/real(trate), ' s; r_keep=', r_keep
            call flush(logfhandle)
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
        allocate(gram_small(r,r), gram_eigvecs_small(r,q_used), gram_eigvals(q_used))
        gram_small = tmp_gram(1:r,1:r)
        if( PROFILE )then
            call system_clock(t1)
            write(logfhandle,'(A,F8.3,A,I8)') 'kPCA Nyström feature/gram build: ', real(t1-t0)/real(trate), ' s; r=', r
            call flush(logfhandle)
            call system_clock(t0)
        endif
        call self%partial_eigh_sym(gram_small, q_used, gram_eigvals, gram_eigvecs_small)
        eig_q(1:q_used)                       = gram_eigvals
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
            call flush(logfhandle)
            call system_clock(t0)
        endif
        if( q_used < self%Q )then
            alpha(:,q_used+1:self%Q) = 0.
            eig_q(q_used+1:self%Q)   = 0.
        endif
        its = MAX_ITS
        if( present(maxpcaits) ) its = maxpcaits
        if( PROFILE )then
            write(logfhandle,'(A,I8,A,I8)') 'kPCA Nyström pre-image start: N=', self%N, ' max_its=', its
            call flush(logfhandle)
            call system_clock(t0)
        endif
        select case(trim(self%kpca_ker))
            case('rbf')
                !$omp parallel do default(shared) proc_bind(close) schedule(static) private(ind,iter,ithr,i,denom,stagn_count,err_prev,err_curr,support_count,k,score)
                do ind = 1, self%N
                    ithr              = omp_get_thread_num() + 1
                    call self%projected_kernel_col(alpha, eig_q, q_used, ind, ker_col(:,ithr))
                    if( local_nbrs > 0 )then
                        call self%select_local_support_inds(ker_col(:,ithr), is_landmark, ind, local_nbrs, support_count, &
                            support_inds(:,ithr), support_scores(:,ithr))
                        support_counts(ithr) = support_count
                    else
                        support_counts(ithr) = 0
                    endif
                    self%data(:,ind)  = pcavecs(:,ind)
                    norm_prev(:,ithr) = 0.
                    err_prev          = huge(1._dp)
                    stagn_count       = 0
                    iter              = 1
                    do while( euclid(self%data(:,ind),norm_prev(:,ithr)) > TOL .and. iter < its )
                        norm_prev(:,ithr) = self%data(:,ind)
                        do i = 1,m
                            proj_data(i,ithr) = euclid(norm_prev(:,ithr), landmark_mat(:,i))**2
                        enddo
                        proj_data(1:m,ithr) = max(0., ker_col(landmark_inds(1:m),ithr)) * exp(-real(rbf_gamma) * proj_data(1:m,ithr))
                        denom             = sum(real(proj_data(1:m,ithr),dp))
                        if( support_counts(ithr) > 0 )then
                            do k = 1, support_counts(ithr)
                                score = real(support_scores(k,ithr), dp) * exp(-real(rbf_gamma) * euclid(norm_prev(:,ithr), pcavecs(:,support_inds(k,ithr)))**2)
                                denom = denom + score
                            enddo
                        endif
                        if( denom < DTINY ) exit
                        self%data(:,ind)  = matmul(landmark_mat, proj_data(1:m,ithr)) / real(denom)
                        if( support_counts(ithr) > 0 )then
                            do k = 1, support_counts(ithr)
                                score = real(support_scores(k,ithr), dp) * exp(-real(rbf_gamma) * euclid(norm_prev(:,ithr), pcavecs(:,support_inds(k,ithr)))**2)
                                self%data(:,ind) = self%data(:,ind) + real(score / denom) * pcavecs(:,support_inds(k,ithr))
                            enddo
                        endif
                        err_curr = euclid(self%data(:,ind), norm_prev(:,ithr))
                        if( abs(err_prev - err_curr) <= EARLY_STOP_REL * max(1._dp, err_prev) )then
                            stagn_count = stagn_count + 1
                            if( stagn_count >= EARLY_STOP_PATIENCE ) exit
                        else
                            stagn_count = 0
                        endif
                        err_prev = err_curr
                        iter = iter + 1
                    enddo
                    if( PROFILE .and. mod(ind, 256) == 0 )then
                        !$omp critical(kpca_nystrom_preimg_progress)
                        call system_clock(t1)
                        write(logfhandle,'(A,I8,A,I8,A,F8.3,A)') 'kPCA Nyström pre-image progress: ', ind, '/', self%N, &
                            '; elapsed=', real(t1-t0)/real(trate), ' s'
                        call flush(logfhandle)
                        !$omp end critical(kpca_nystrom_preimg_progress)
                    endif
                enddo
                !$omp end parallel do
            case('cosine')
                if( trim(self%kpca_target) .eq. 'ptcl' )then
                    !$omp parallel do default(shared) proc_bind(close) schedule(static) private(ind,ithr,iter,denom,stagn_count,err_prev,err_curr,support_count,k,score)
                    do ind = 1, self%N
                        ithr              = omp_get_thread_num() + 1
                        call self%projected_kernel_col(alpha, eig_q, q_used, ind, ker_col(:,ithr))
                        if( local_nbrs > 0 )then
                            call self%select_local_support_inds(ker_col(:,ithr), is_landmark, ind, local_nbrs, support_count, &
                                support_inds(:,ithr), support_scores(:,ithr))
                            support_counts(ithr) = support_count
                        else
                            support_counts(ithr) = 0
                        endif
                        self%data(:,ind)  = pcavecs(:,ind)
                        norm_prev(:,ithr) = 1. / sqrt(real(self%D))
                        norm_data(:,ithr) = norm_pcavecs(:,ind)
                        err_prev          = huge(1._dp)
                        stagn_count       = 0
                        iter              = 1
                        do while( abs(sum(norm_data(:,ithr) * norm_prev(:,ithr)) - 1.) > TOL .and. iter < its )
                            norm_prev(:,ithr) = norm_data(:,ithr)
                            proj_data(1:m,ithr) = matmul(norm_landmarks_t, norm_prev(:,ithr)) * ker_col(landmark_inds(1:m),ithr)
                            denom = sum(abs(real(proj_data(1:m,ithr),dp)))
                            self%data(:,ind)  = matmul(landmark_mat, proj_data(1:m,ithr))
                            if( support_counts(ithr) > 0 )then
                                do k = 1, support_counts(ithr)
                                    score = real(support_scores(k,ithr), dp) * sum(real(norm_pcavecs(:,support_inds(k,ithr)),dp) * real(norm_prev(:,ithr),dp))
                                    self%data(:,ind) = self%data(:,ind) + real(score) * pcavecs(:,support_inds(k,ithr))
                                    denom = denom + abs(score)
                                enddo
                            endif
                            if( denom < DTINY ) exit
                            self%data(:,ind)  = self%data(:,ind) / real(denom)
                            denom             = dsqrt(sum(real(self%data(:,ind),dp)**2))
                            if( denom < DTINY ) exit
                            norm_data(:,ithr) = self%data(:,ind) / real(denom)
                            err_curr = abs(sum(real(norm_data(:,ithr),dp) * real(norm_prev(:,ithr),dp)) - 1._dp)
                            if( abs(err_prev - err_curr) <= EARLY_STOP_REL * max(1._dp, err_prev) )then
                                stagn_count = stagn_count + 1
                                if( stagn_count >= EARLY_STOP_PATIENCE ) exit
                            else
                                stagn_count = 0
                            endif
                            err_prev = err_curr
                            iter = iter + 1
                        enddo
                        if( PROFILE .and. mod(ind, 256) == 0 )then
                            !$omp critical(kpca_nystrom_preimg_progress)
                            call system_clock(t1)
                            write(logfhandle,'(A,I8,A,I8,A,F8.3,A)') 'kPCA Nyström pre-image progress: ', ind, '/', self%N, &
                                '; elapsed=', real(t1-t0)/real(trate), ' s'
                            call flush(logfhandle)
                            !$omp end critical(kpca_nystrom_preimg_progress)
                        endif
                    enddo
                    !$omp end parallel do
                else
                    !$omp parallel do default(shared) proc_bind(close) schedule(static) private(ind,ithr,iter,denom,stagn_count,err_prev,err_curr,support_count,k,score)
                    do ind = 1, self%N
                        ithr              = omp_get_thread_num() + 1
                        call self%projected_kernel_col(alpha, eig_q, q_used, ind, ker_col(:,ithr))
                        if( local_nbrs > 0 )then
                            call self%select_local_support_inds(ker_col(:,ithr), is_landmark, ind, local_nbrs, support_count, &
                                support_inds(:,ithr), support_scores(:,ithr))
                            support_counts(ithr) = support_count
                        else
                            support_counts(ithr) = 0
                        endif
                        self%data(:,ind)  = pcavecs(:,ind)
                        norm_prev(:,ithr) = 1. / sqrt(real(self%D))
                        norm_data(:,ithr) = norm_pcavecs(:,ind)
                        err_prev          = huge(1._dp)
                        stagn_count       = 0
                        iter              = 1
                        do while( abs(sum(norm_data(:,ithr) * norm_prev(:,ithr)) - 1.) > TOL .and. iter < its )
                            norm_prev(:,ithr) = norm_data(:,ithr)
                            proj_data(1:m,ithr) = matmul(norm_landmarks_t, norm_prev(:,ithr)) * ker_col(landmark_inds(1:m),ithr)
                            self%data(:,ind)  = matmul(norm_landmarks, proj_data(1:m,ithr))
                            if( support_counts(ithr) > 0 )then
                                do k = 1, support_counts(ithr)
                                    score = real(support_scores(k,ithr), dp) * sum(real(norm_pcavecs(:,support_inds(k,ithr)),dp) * real(norm_prev(:,ithr),dp))
                                    self%data(:,ind) = self%data(:,ind) + real(score) * norm_pcavecs(:,support_inds(k,ithr))
                                enddo
                            endif
                            denom             = dsqrt(sum(real(self%data(:,ind),dp)**2))
                            if( denom < DTINY ) exit
                            norm_data(:,ithr) = self%data(:,ind) / real(denom)
                            err_curr = abs(sum(real(norm_data(:,ithr),dp) * real(norm_prev(:,ithr),dp)) - 1._dp)
                            if( abs(err_prev - err_curr) <= EARLY_STOP_REL * max(1._dp, err_prev) )then
                                stagn_count = stagn_count + 1
                                if( stagn_count >= EARLY_STOP_PATIENCE ) exit
                            else
                                stagn_count = 0
                            endif
                            err_prev = err_curr
                            iter = iter + 1
                        enddo
                        if( PROFILE .and. mod(ind, 256) == 0 )then
                            !$omp critical(kpca_nystrom_preimg_progress)
                            call system_clock(t1)
                            write(logfhandle,'(A,I8,A,I8,A,F8.3,A)') 'kPCA Nyström pre-image progress: ', ind, '/', self%N, &
                                '; elapsed=', real(t1-t0)/real(trate), ' s'
                            call flush(logfhandle)
                            !$omp end critical(kpca_nystrom_preimg_progress)
                        endif
                    enddo
                    !$omp end parallel do
                endif
        end select
        if( PROFILE )then
            call system_clock(t1)
            write(logfhandle,'(A,F8.3,A)') 'kPCA Nyström pre-image/reconstruct: ', real(t1-t0)/real(trate), ' s'
            call flush(logfhandle)
        endif
        if( allocated(norm_pcavecs)   ) deallocate(norm_pcavecs)
        if( allocated(norm_landmarks) ) deallocate(norm_landmarks)
        if( allocated(norm_landmarks_t) ) deallocate(norm_landmarks_t)
        if( allocated(support_scores) ) deallocate(support_scores)
        if( allocated(support_inds) ) deallocate(support_inds)
        if( allocated(support_counts) ) deallocate(support_counts)
        if( allocated(is_landmark) ) deallocate(is_landmark)
        deallocate(eig_q, alpha)
        deallocate(ker_nm, ker_mm, feat, feat_center, eig_w, eigvec_w, tmp_ker_mm, tmp_gram, gram_eigvecs, landmark_mat, &
                   &gram_small, gram_eigvecs_small, gram_eigvals, ker_col, proj_data, norm_prev, norm_data)
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
        real, allocatable :: eig_vals_raw(:), eig_vecs_raw(:,:)
        integer :: i
        tmp_ker = ker
        ! computing eigvals/eigvecs
        allocate(eig_vals_raw(self%Q), eig_vecs_raw(self%N,self%Q))
        call eigh(self%N, tmp_ker, self%Q, eig_vals_raw, eig_vecs_raw)
        eig_vals = eig_vals_raw(self%Q:1:-1)
        eig_vecs = eig_vecs_raw(:,self%Q:1:-1)
        deallocate(eig_vals_raw, eig_vecs_raw)
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

    integer function resolve_nystrom_npts( n_samples, q_hint, npts_user ) result( m )
        integer, intent(in) :: n_samples, q_hint, npts_user
        if( npts_user > 0 )then
            m = npts_user
        else
            m = max(NYSTROM_AUTO_NPTS, 2 * max(q_hint, 1))
        endif
        m = min(n_samples, min(512, max(max(q_hint, 1), m)))
    end function resolve_nystrom_npts

    real(dp) function auto_rbf_gamma_from_data( mat ) result( gamma )
        real, intent(in) :: mat(:,:)
        integer  :: i, j, npairs, n
        real(dp) :: mean_sqdist
        n = size(mat, 2)
        mean_sqdist = 0._dp
        npairs      = 0
        do j = 1,n
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
    end function auto_rbf_gamma_from_data

    integer function choose_auto_q_from_eigs( eigvals, energy_frac ) result( q_auto )
        real,    intent(in) :: eigvals(:)
        real(dp), intent(in) :: energy_frac
        integer  :: i, npos
        real(dp) :: total_energy, running_energy
        npos = count(real(eigvals,dp) > real(DTINY,dp))
        if( npos < 1 )then
            q_auto = 1
            return
        endif
        total_energy = sum(real(eigvals(1:npos),dp))
        if( total_energy <= real(DTINY,dp) )then
            q_auto = 1
            return
        endif
        running_energy = 0._dp
        q_auto = npos
        do i = 1,npos
            running_energy = running_energy + real(eigvals(i),dp)
            if( running_energy / total_energy >= energy_frac )then
                q_auto = i
                exit
            endif
        enddo
        q_auto = max(1, q_auto)
    end function choose_auto_q_from_eigs

    integer function suggest_kpca_nystrom_neigs( pcavecs, kpca_ker, kpca_nystrom_npts, kpca_rbf_gamma ) result( q_auto )
        real,             intent(in) :: pcavecs(:,:)
        character(len=*), intent(in) :: kpca_ker
        integer, optional, intent(in) :: kpca_nystrom_npts
        real,    optional, intent(in) :: kpca_rbf_gamma
        real(dp), parameter :: ENERGY_FRAC = 0.99_dp
        integer  :: n, d, m, r_keep, r, i, j
        integer, allocatable :: landmark_inds(:)
        real(dp) :: gamma, denom, eig_tol
        real, allocatable :: ker_nm(:,:), ker_mm(:,:), feat(:,:), feat_center(:), eig_w(:), eigvec_w(:,:), tmp_ker_mm(:,:)
        real, allocatable :: tmp_gram(:,:), landmark_eigvals(:), landmark_eigvecs(:,:), gram_small(:,:), gram_eigvals(:), gram_eigvecs(:,:)
        real, allocatable :: norm_pcavecs(:,:), norm_landmarks(:,:)
        n = size(pcavecs, 2)
        d = size(pcavecs, 1)
        if( n <= 1 )then
            q_auto = 1
            return
        endif
        if( present(kpca_nystrom_npts) )then
            m = resolve_nystrom_npts(n, 0, kpca_nystrom_npts)
        else
            m = resolve_nystrom_npts(n, 0, 0)
        endif
        r_keep = min(m, max(16, min(m, 64)))
        allocate(landmark_inds(m), source=0)
        call choose_nystrom_inds_from_data(m, n, d, landmark_inds, pcavecs)
        allocate(ker_nm(n,m), ker_mm(m,m), feat(n,m), feat_center(m), eig_w(m), eigvec_w(m,m), tmp_ker_mm(m,m), tmp_gram(m,m), source=0.)
        select case(trim(kpca_ker))
            case('cosine')
                allocate(norm_pcavecs(d,n), norm_landmarks(d,m), source=0.)
                norm_pcavecs = pcavecs
                do i = 1,n
                    denom = dsqrt(sum(real(pcavecs(:,i),dp)**2))
                    if( denom > DTINY ) norm_pcavecs(:,i) = pcavecs(:,i) / real(denom)
                enddo
                norm_landmarks = norm_pcavecs(:,landmark_inds)
                ker_nm = matmul(transpose(norm_pcavecs), norm_landmarks)
                ker_mm = matmul(transpose(norm_landmarks), norm_landmarks)
                deallocate(norm_pcavecs, norm_landmarks)
            case('rbf')
                if( present(kpca_rbf_gamma) )then
                    gamma = real(kpca_rbf_gamma, dp)
                else
                    gamma = 0._dp
                endif
                if( gamma <= 0._dp ) gamma = auto_rbf_gamma_from_data(pcavecs)
                do j = 1,m
                    do i = 1,n
                        ker_nm(i,j) = exp(-real(gamma) * euclid(pcavecs(:,i), pcavecs(:,landmark_inds(j)))**2)
                    enddo
                enddo
                do j = 1,m
                    ker_mm(j,j) = 1.
                    do i = 1,j-1
                        ker_mm(i,j) = exp(-real(gamma) * euclid(pcavecs(:,landmark_inds(i)), pcavecs(:,landmark_inds(j)))**2)
                        ker_mm(j,i) = ker_mm(i,j)
                    enddo
                enddo
            case default
                q_auto = min(n-1, max(8, min(64, m/2)))
                deallocate(landmark_inds, ker_nm, ker_mm, feat, feat_center, eig_w, eigvec_w, tmp_ker_mm, tmp_gram)
                return
        end select
        allocate(landmark_eigvals(r_keep), landmark_eigvecs(m,r_keep))
        tmp_ker_mm = ker_mm
        call eigh(m, tmp_ker_mm, r_keep, landmark_eigvals, landmark_eigvecs)
        eig_w(1:r_keep)      = landmark_eigvals(r_keep:1:-1)
        eigvec_w(:,1:r_keep) = landmark_eigvecs(:,r_keep:1:-1)
        deallocate(landmark_eigvals, landmark_eigvecs)
        eig_tol = max(real(DTINY,dp), 1.e-6_dp * max(real(maxval(eig_w(1:r_keep)),dp), 1._dp))
        r = count(real(eig_w(1:r_keep),dp) > eig_tol)
        if( r < 1 )then
            q_auto = 1
            deallocate(landmark_inds, ker_nm, ker_mm, feat, feat_center, eig_w, eigvec_w, tmp_ker_mm, tmp_gram)
            return
        endif
        feat(:,1:r) = matmul(ker_nm, eigvec_w(:,1:r))
        do i = 1,r
            feat(:,i) = feat(:,i) / sqrt(max(eig_w(i), real(DTINY)))
        enddo
        feat_center(1:r) = sum(feat(:,1:r), dim=1) / real(n)
        do i = 1,r
            feat(:,i) = feat(:,i) - feat_center(i)
        enddo
        allocate(gram_small(r,r), gram_eigvals(r), gram_eigvecs(r,r))
        gram_small = matmul(transpose(feat(:,1:r)), feat(:,1:r))
        call eigh(r, gram_small, r, gram_eigvals, gram_eigvecs)
        gram_eigvals = gram_eigvals(r:1:-1)
        q_auto = choose_auto_q_from_eigs(gram_eigvals, ENERGY_FRAC)
        q_auto = min(max(q_auto, min(8, r)), r)
        deallocate(landmark_inds, ker_nm, ker_mm, feat, feat_center, eig_w, eigvec_w, tmp_ker_mm, tmp_gram, gram_small, gram_eigvals, gram_eigvecs)
    end function suggest_kpca_nystrom_neigs

    real(dp) function get_rbf_gamma( self, mat ) result( gamma )
        class(kpca_svd), intent(in) :: self
        real,            intent(in) :: mat(self%D,self%N)
        if( self%kpca_rbf_gamma > 0. )then
            gamma = real(self%kpca_rbf_gamma, dp)
            return
        endif
        gamma = auto_rbf_gamma_from_data(mat)
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

    subroutine select_nystrom_inds( self, m, inds, pcavecs )
        class(kpca_svd), intent(inout) :: self
        integer,         intent(in)    :: m
        integer,         intent(out)   :: inds(m)
        real, optional,  intent(in)    :: pcavecs(self%D,self%N)
        if( present(pcavecs) )then
            call choose_nystrom_inds_from_data(m, self%N, self%D, inds, pcavecs)
        else
            call choose_nystrom_inds_from_data(m, self%N, 0, inds)
        endif
    end subroutine select_nystrom_inds

    subroutine choose_nystrom_inds_from_data( m, n, d, inds, pcavecs )
        integer,         intent(in)  :: m, n, d
        integer,         intent(out) :: inds(m)
        real, optional,  intent(in)  :: pcavecs(d,n)
        integer, parameter :: SKETCH_DIM = 32
        integer :: i, j, k, ncand, nsketch, ibeg, iend, best_i, best_pos, cand_idx
        integer, allocatable :: candidates(:), sketch_pos(:)
        logical, allocatable :: chosen(:)
        real(dp), allocatable :: cand_norms(:), min_dist(:), sketch(:,:)
        real(dp) :: cur_norm, best_norm, dist_ij, delta, nrm
        if( m <= 1 .or. n <= 1 )then
            inds(1) = 1
            return
        endif
        if( .not. present(pcavecs) .or. d <= 0 .or. m >= n )then
            do i = 1,m
                inds(i) = 1 + ((i-1) * (n-1)) / max(1, m-1)
            enddo
            return
        endif
        ncand = min(n, max(m, min(4*m, n)))
        allocate(candidates(ncand), cand_norms(ncand), chosen(ncand))
        candidates = 0
        cand_norms = 0._dp
        chosen     = .false.
        do i = 1,ncand
            ibeg = 1 + ((i-1) * n) / ncand
            iend = max(ibeg, (i * n) / ncand)
            best_i = ibeg
            best_norm = -1._dp
            do k = ibeg, iend
                cur_norm = sum(real(pcavecs(:,k),dp) * real(pcavecs(:,k),dp))
                if( cur_norm > best_norm )then
                    best_norm = cur_norm
                    best_i    = k
                endif
            enddo
            candidates(i) = best_i
            cand_norms(i) = max(best_norm, DTINY)
        enddo
        nsketch = min(SKETCH_DIM, d)
        allocate(sketch_pos(nsketch), sketch(nsketch,ncand), min_dist(ncand))
        sketch_pos = 0
        sketch     = 0._dp
        min_dist   = 0._dp
        if( nsketch == 1 )then
            sketch_pos(1) = 1
        else
            do i = 1,nsketch
                sketch_pos(i) = 1 + ((i-1) * (d-1)) / (nsketch-1)
            enddo
        endif
        do i = 1,ncand
            cand_idx = candidates(i)
            nrm = sqrt(cand_norms(i))
            do j = 1,nsketch
                sketch(j,i) = real(pcavecs(sketch_pos(j), cand_idx), dp) / nrm
            enddo
        enddo
        best_pos = 1
        best_norm = -1._dp
        do i = 1,ncand
            cur_norm = sum(sketch(:,i) * sketch(:,i))
            if( cur_norm > best_norm )then
                best_norm = cur_norm
                best_pos  = i
            endif
        enddo
        inds(1) = candidates(best_pos)
        chosen(best_pos) = .true.
        do i = 1,ncand
            dist_ij = 0._dp
            do j = 1,nsketch
                delta = sketch(j,i) - sketch(j,best_pos)
                dist_ij = dist_ij + delta * delta
            enddo
            min_dist(i) = dist_ij
        enddo
        min_dist(best_pos) = -1._dp
        do k = 2,m
            best_pos = 1
            best_norm = -1._dp
            do i = 1,ncand
                if( .not. chosen(i) .and. min_dist(i) > best_norm )then
                    best_norm = min_dist(i)
                    best_pos  = i
                endif
            enddo
            inds(k) = candidates(best_pos)
            chosen(best_pos) = .true.
            do i = 1,ncand
                if( .not. chosen(i) )then
                    dist_ij = 0._dp
                    do j = 1,nsketch
                        delta = sketch(j,i) - sketch(j,best_pos)
                        dist_ij = dist_ij + delta * delta
                    enddo
                    min_dist(i) = min(min_dist(i), dist_ij)
                endif
            enddo
        enddo
        deallocate(candidates, cand_norms, chosen, sketch_pos, sketch, min_dist)
    end subroutine choose_nystrom_inds_from_data

    subroutine select_local_support_inds( self, weights, is_landmark, self_ind, max_keep, n_keep, inds, vals )
        class(kpca_svd), intent(in)  :: self
        real,            intent(in)  :: weights(:)
        logical,         intent(in)  :: is_landmark(:)
        integer,         intent(in)  :: self_ind, max_keep
        integer,         intent(out) :: n_keep, inds(max_keep)
        real,            intent(out) :: vals(max_keep)
        real(dp), parameter :: KEEP_FRAC = 0.90_dp
        integer :: i, j, pos
        real    :: w, tmpv
        integer :: tmpi
        real(dp) :: total_pos, keep_target, running
        inds = 0
        vals = 0.
        n_keep = 0
        if( max_keep <= 0 ) return
        total_pos = 0._dp
        do i = 1, size(weights)
            if( i == self_ind ) cycle
            if( is_landmark(i) ) cycle
            w = weights(i)
            if( w <= 0. ) cycle
            total_pos = total_pos + real(w, dp)
            if( n_keep < max_keep )then
                n_keep = n_keep + 1
                inds(n_keep) = i
                vals(n_keep) = w
            else
                pos = 1
                do j = 2, max_keep
                    if( vals(j) < vals(pos) ) pos = j
                enddo
                if( w > vals(pos) )then
                    inds(pos) = i
                    vals(pos) = w
                endif
            endif
        enddo
        if( n_keep <= 1 ) return
        do i = 1, n_keep-1
            do j = i+1, n_keep
                if( vals(j) > vals(i) )then
                    tmpv    = vals(i)
                    vals(i) = vals(j)
                    vals(j) = tmpv
                    tmpi    = inds(i)
                    inds(i) = inds(j)
                    inds(j) = tmpi
                endif
            enddo
        enddo
        keep_target = KEEP_FRAC * total_pos
        running = 0._dp
        do i = 1, n_keep
            running = running + real(vals(i), dp)
            if( running >= keep_target )then
                n_keep = i
                exit
            endif
        enddo
        if( n_keep < max_keep )then
            inds(n_keep+1:max_keep) = 0
            vals(n_keep+1:max_keep) = 0.
        endif
    end subroutine select_local_support_inds

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
