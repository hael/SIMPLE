!@descr: complex/Hermitian probabilistic PCA for Fourier-domain vectors
!
! Key ideas / intended use:
! - Streaming comes first: this type is designed to accumulate sufficient
!   statistics over enormous numbers of complex vectors without storing the
!   full data matrix in memory.
! - Complex (Hermitian) formulation: we work directly with complex Fourier
!   vectors and Hermitian covariance matrices instead of splitting into
!   real/imaginary parts.
! - Probabilistic PCA view: after fitting the Hermitian covariance spectrum,
!   we interpret the retained subspace with an isotropic residual noise model.
! - Wiener-style denoising: denoising is done by MMSE shrinkage in the complex
!   eigenbasis rather than by bare hard projection alone.
! - Fourier-domain geometry: this preserves phase-aware linear structure and is
!   a better match to polar Fourier line data than ordinary real PCA.
! - Scalable path for line denoising: the intended workflow for cryo-EM line
!   data is stream -> accumulate mean/covariance -> finalize fit -> denoise.
module simple_complex_ppca
use simple_core_module_api
implicit none

public :: complex_ppca
private

type :: complex_ppca
    private
    integer :: N = 0
    integer :: D = 0
    integer :: Q = 0
    complex(dp), allocatable :: mu(:)       !< complex mean
    complex(dp), allocatable :: U(:,:)      !< Hermitian eigenvectors, descending eigenvalues
    complex(dp), allocatable :: W(:,:)      !< PPCA loading matrix
    complex(dp), allocatable :: E_zn(:,:)   !< posterior latent means
    complex(dp), allocatable :: data(:,:)   !< reconstructed centered data
    real(dp),    allocatable :: eigvals(:)
    real(dp),    allocatable :: shrink(:)   !< MMSE clean-signal shrink factors in eigenbasis
    complex(dp), allocatable :: sum_x(:)    !< streaming first moment accumulator
    complex(dp), allocatable :: sum_xx(:,:) !< streaming second moment accumulator
    real(dp) :: sigma2 = 1._dp
    integer(kind=8) :: nobs = 0_8
    real(dp) :: sum_w = 0._dp
    logical  :: stream_ready = .false.
    logical  :: existence = .false.
  contains
    procedure :: new                  => new_complex_ppca
    procedure :: master               => master_complex_ppca
    procedure :: stream_reset         => stream_reset_complex_ppca
    procedure :: stream_update        => stream_update_complex_ppca
    procedure :: stream_finalize      => stream_finalize_complex_ppca
    procedure :: fit_covariance       => fit_covariance_complex_ppca
    procedure :: denoise              => denoise_complex_ppca
    procedure :: reconstruct_external => reconstruct_external_complex_ppca
    procedure :: get_feat             => get_feat_complex_ppca
    procedure :: generate             => generate_complex_ppca
    procedure :: get_mean             => get_mean_complex_ppca
    procedure :: get_eigvals          => get_eigvals_complex_ppca
    procedure :: get_sigma2           => get_sigma2_complex_ppca
    procedure :: kill                 => kill_complex_ppca
end type complex_ppca

contains

    subroutine new_complex_ppca( self, N, D, Q )
        class(complex_ppca), intent(inout) :: self
        integer,             intent(in)    :: N, D, Q
        call self%kill
        self%N = N
        self%D = D
        self%Q = min(max(Q, 1), D)
        allocate(self%mu(self%D), self%U(self%D,self%D), self%W(self%D,self%Q), &
                 self%E_zn(self%Q,max(1,self%N)), self%data(self%D,max(1,self%N)), &
                 source=cmplx(0._dp, 0._dp, kind=dp))
        allocate(self%eigvals(self%D), self%shrink(self%Q), source=0._dp)
        allocate(self%sum_x(self%D), self%sum_xx(self%D,self%D), source=cmplx(0._dp, 0._dp, kind=dp))
        self%sigma2 = 1._dp
        self%nobs = 0_8
        self%sum_w = 0._dp
        self%stream_ready = .false.
        self%existence = .true.
    end subroutine new_complex_ppca

    subroutine master_complex_ppca( self, cvecs )
        class(complex_ppca), intent(inout) :: self
        complex(dp),         intent(in)    :: cvecs(self%D,self%N)
        integer :: i
        if( .not. self%existence )then
            call simple_exception('complex_ppca not allocated; master', 'simple_complex_ppca.f90', __LINE__)
        endif
        call self%stream_reset()
        do i = 1,self%N
            call self%stream_update(cvecs(:,i))
        enddo
        call self%stream_finalize()
        do i = 1,self%N
            call self%reconstruct_external(cvecs(:,i) - self%mu, self%data(:,i), self%E_zn(:,i))
        enddo
    end subroutine master_complex_ppca

    subroutine stream_reset_complex_ppca( self )
        class(complex_ppca), intent(inout) :: self
        if( .not. self%existence )then
            call simple_exception('complex_ppca not allocated; stream_reset', 'simple_complex_ppca.f90', __LINE__)
        endif
        self%sum_x = cmplx(0._dp, 0._dp, kind=dp)
        self%sum_xx = cmplx(0._dp, 0._dp, kind=dp)
        self%mu = cmplx(0._dp, 0._dp, kind=dp)
        self%nobs = 0_8
        self%sum_w = 0._dp
        self%stream_ready = .false.
    end subroutine stream_reset_complex_ppca

    subroutine stream_update_complex_ppca( self, x, weight )
        class(complex_ppca), intent(inout) :: self
        complex(dp),         intent(in)    :: x(self%D)
        real(dp), optional,  intent(in)    :: weight
        real(dp) :: w
        if( .not. self%existence )then
            call simple_exception('complex_ppca not allocated; stream_update', 'simple_complex_ppca.f90', __LINE__)
        endif
        w = 1._dp
        if( present(weight) ) w = weight
        if( w <= 0._dp ) return
        self%sum_x  = self%sum_x  + cmplx(w, 0._dp, kind=dp) * x
        self%sum_xx = self%sum_xx + cmplx(w, 0._dp, kind=dp) * outer_hermitian(x)
        self%nobs   = self%nobs + 1_8
        self%sum_w  = self%sum_w + w
    end subroutine stream_update_complex_ppca

    subroutine stream_finalize_complex_ppca( self, sigma2_in )
        class(complex_ppca), intent(inout) :: self
        real(dp), optional,  intent(in)    :: sigma2_in
        complex(dp) :: mean_here(self%D), cov(self%D,self%D)
        real(dp) :: inv_n
        if( .not. self%existence )then
            call simple_exception('complex_ppca not allocated; stream_finalize', 'simple_complex_ppca.f90', __LINE__)
        endif
        if( self%nobs < 1_8 .or. self%sum_w <= 0._dp )then
            call simple_exception('cannot finalize empty complex_ppca stream', 'simple_complex_ppca.f90', __LINE__)
        endif
        inv_n = 1._dp / self%sum_w
        mean_here = self%sum_x * cmplx(inv_n, 0._dp, kind=dp)
        cov = self%sum_xx * cmplx(inv_n, 0._dp, kind=dp) - outer_hermitian(mean_here)
        if( present(sigma2_in) )then
            call self%fit_covariance(mean_here, cov, sigma2_in)
        else
            call self%fit_covariance(mean_here, cov)
        endif
        self%stream_ready = .true.
    end subroutine stream_finalize_complex_ppca

    subroutine fit_covariance_complex_ppca( self, mu, cov, sigma2_in )
        class(complex_ppca), intent(inout) :: self
        complex(dp),         intent(in)    :: mu(self%D)
        complex(dp),         intent(in)    :: cov(self%D,self%D)
        real(dp), optional,  intent(in)    :: sigma2_in
        complex(dp) :: hcov(self%D,self%D)
        real(dp)    :: signal_var
        integer :: j, ntail
        if( .not. self%existence )then
            call simple_exception('complex_ppca not allocated; fit_covariance', 'simple_complex_ppca.f90', __LINE__)
        endif
        self%mu = mu
        hcov = 0.5_dp * (cov + transpose(conjg(cov)))
        call hermitian_jacobi(hcov, self%eigvals, self%U)
        call sort_eigs_desc(self%eigvals, self%U)
        if( present(sigma2_in) )then
            self%sigma2 = max(sigma2_in, real(DTINY,dp))
        else if( self%Q < self%D )then
            ntail = self%D - self%Q
            self%sigma2 = max(sum(self%eigvals(self%Q+1:self%D)) / real(ntail,dp), real(DTINY,dp))
        else
            self%sigma2 = max(minval(self%eigvals), real(DTINY,dp))
        endif
        do j = 1,self%Q
            signal_var = max(self%eigvals(j) - self%sigma2, 0._dp)
            self%W(:,j) = self%U(:,j) * cmplx(sqrt(signal_var), 0._dp, kind=dp)
            if( self%eigvals(j) > real(DTINY,dp) )then
                self%shrink(j) = signal_var / self%eigvals(j)
            else
                self%shrink(j) = 0._dp
            endif
        enddo
    end subroutine fit_covariance_complex_ppca

    subroutine denoise_complex_ppca( self, x, xhat )
        class(complex_ppca), intent(in)  :: self
        complex(dp),         intent(in)  :: x(self%D)
        complex(dp),         intent(out) :: xhat(self%D)
        complex(dp) :: centered(self%D), recon(self%D)
        if( .not. self%existence )then
            call simple_exception('complex_ppca not allocated; denoise', 'simple_complex_ppca.f90', __LINE__)
        endif
        if( .not. self%stream_ready )then
            call simple_exception('complex_ppca fit not finalized; denoise', 'simple_complex_ppca.f90', __LINE__)
        endif
        centered = x - self%mu
        call self%reconstruct_external(centered, recon)
        xhat = self%mu + recon
    end subroutine denoise_complex_ppca

    subroutine reconstruct_external_complex_ppca( self, centered_in, centered_out, feat )
        class(complex_ppca),           intent(in)  :: self
        complex(dp),                   intent(in)  :: centered_in(self%D)
        complex(dp),                   intent(out) :: centered_out(self%D)
        complex(dp), optional, target, intent(out) :: feat(self%Q)
        complex(dp) :: zloc(self%Q)
        integer :: j
        if( .not. self%existence )then
            call simple_exception('complex_ppca not allocated; reconstruct_external', 'simple_complex_ppca.f90', __LINE__)
        endif
        if( .not. self%stream_ready )then
            call simple_exception('complex_ppca fit not finalized; reconstruct_external', 'simple_complex_ppca.f90', __LINE__)
        endif
        do j = 1,self%Q
            zloc(j) = dot_product(self%U(:,j), centered_in)
            zloc(j) = cmplx(self%shrink(j), 0._dp, kind=dp) * zloc(j)
        enddo
        centered_out = matmul(self%U(:,1:self%Q), zloc)
        if( present(feat) ) feat = zloc
    end subroutine reconstruct_external_complex_ppca

    function get_feat_complex_ppca( self, i ) result( feat )
        class(complex_ppca), intent(in) :: self
        integer,             intent(in) :: i
        complex(dp), allocatable :: feat(:)
        if( .not. self%existence )then
            call simple_exception('complex_ppca not allocated; get_feat', 'simple_complex_ppca.f90', __LINE__)
        endif
        allocate(feat(self%Q), source=self%E_zn(:,i))
    end function get_feat_complex_ppca

    subroutine generate_complex_ppca( self, i, avg, dat )
        class(complex_ppca), intent(in)    :: self
        integer,             intent(in)    :: i
        complex(dp),         intent(in)    :: avg(self%D)
        complex(dp),         intent(inout) :: dat(self%D)
        if( .not. self%existence )then
            call simple_exception('complex_ppca not allocated; generate', 'simple_complex_ppca.f90', __LINE__)
        endif
        dat = avg + self%data(:,i)
    end subroutine generate_complex_ppca

    function get_mean_complex_ppca( self ) result( mu )
        class(complex_ppca), intent(in) :: self
        complex(dp), allocatable :: mu(:)
        allocate(mu(self%D), source=self%mu)
    end function get_mean_complex_ppca

    function get_eigvals_complex_ppca( self ) result( eigvals )
        class(complex_ppca), intent(in) :: self
        real(dp), allocatable :: eigvals(:)
        allocate(eigvals(self%D), source=self%eigvals)
    end function get_eigvals_complex_ppca

    real(dp) function get_sigma2_complex_ppca( self ) result( sigma2 )
        class(complex_ppca), intent(in) :: self
        sigma2 = self%sigma2
    end function get_sigma2_complex_ppca

    subroutine kill_complex_ppca( self )
        class(complex_ppca), intent(inout) :: self
        if( allocated(self%mu)      ) deallocate(self%mu)
        if( allocated(self%U)       ) deallocate(self%U)
        if( allocated(self%W)       ) deallocate(self%W)
        if( allocated(self%E_zn)    ) deallocate(self%E_zn)
        if( allocated(self%data)    ) deallocate(self%data)
        if( allocated(self%eigvals) ) deallocate(self%eigvals)
        if( allocated(self%shrink)  ) deallocate(self%shrink)
        if( allocated(self%sum_x)   ) deallocate(self%sum_x)
        if( allocated(self%sum_xx)  ) deallocate(self%sum_xx)
        self%N = 0
        self%D = 0
        self%Q = 0
        self%sigma2 = 1._dp
        self%nobs = 0_8
        self%sum_w = 0._dp
        self%stream_ready = .false.
        self%existence = .false.
    end subroutine kill_complex_ppca

    pure function outer_hermitian( x ) result( cov )
        complex(dp), intent(in) :: x(:)
        complex(dp) :: cov(size(x),size(x))
        integer :: i, j
        do j = 1,size(x)
            do i = 1,size(x)
                cov(i,j) = x(i) * conjg(x(j))
            enddo
        enddo
    end function outer_hermitian

    subroutine hermitian_jacobi( A, evals, evecs )
        complex(dp), intent(in)  :: A(:,:)
        real(dp),    intent(out) :: evals(size(A,1))
        complex(dp), intent(out) :: evecs(size(A,1),size(A,1))
        integer, parameter :: MAX_SWEEPS = 80
        real(dp), parameter :: TOL = 1.e-12_dp
        complex(dp) :: B(size(A,1),size(A,2))
        complex(dp) :: cphase, vip, viq, apq
        real(dp) :: c, s, tau, t, app, aqq, offmax, abs_apq
        integer :: n, p, q, i, sweep
        n = size(A,1)
        B = A
        evecs = cmplx(0._dp, 0._dp, kind=dp)
        do i = 1,n
            evecs(i,i) = cmplx(1._dp, 0._dp, kind=dp)
        enddo
        do sweep = 1,MAX_SWEEPS
            offmax = 0._dp
            do q = 2,n
                do p = 1,q-1
                    offmax = max(offmax, abs(B(p,q)))
                enddo
            enddo
            if( offmax <= TOL * max(1._dp, maxval(abs([(real(B(i,i),dp), i=1,n)]))) ) exit
            do q = 2,n
                do p = 1,q-1
                    apq = B(p,q)
                    abs_apq = abs(apq)
                    if( abs_apq <= TOL ) cycle
                    app = real(B(p,p),dp)
                    aqq = real(B(q,q),dp)
                    tau = (aqq - app) / (2._dp * abs_apq)
                    t = sign(1._dp, tau) / (abs(tau) + sqrt(1._dp + tau*tau))
                    if( tau == 0._dp ) t = 1._dp
                    c = 1._dp / sqrt(1._dp + t*t)
                    s = c * t
                    cphase = apq / abs_apq
                    do i = 1,n
                        if( i /= p .and. i /= q )then
                            vip = B(i,p)
                            viq = B(i,q)
                            B(i,p) = c * vip - s * cphase * viq
                            B(p,i) = conjg(B(i,p))
                            B(i,q) = c * viq + s * conjg(cphase) * vip
                            B(q,i) = conjg(B(i,q))
                        endif
                    enddo
                    B(p,p) = cmplx(c*c*app + s*s*aqq - 2._dp*c*s*abs_apq, 0._dp, kind=dp)
                    B(q,q) = cmplx(s*s*app + c*c*aqq + 2._dp*c*s*abs_apq, 0._dp, kind=dp)
                    B(p,q) = cmplx(0._dp, 0._dp, kind=dp)
                    B(q,p) = cmplx(0._dp, 0._dp, kind=dp)
                    do i = 1,n
                        vip = evecs(i,p)
                        viq = evecs(i,q)
                        evecs(i,p) = c * vip - s * cphase * viq
                        evecs(i,q) = c * viq + s * conjg(cphase) * vip
                    enddo
                enddo
            enddo
        enddo
        do i = 1,n
            evals(i) = real(B(i,i),dp)
        enddo
    end subroutine hermitian_jacobi

    subroutine sort_eigs_desc( evals, evecs )
        real(dp),    intent(inout) :: evals(:)
        complex(dp), intent(inout) :: evecs(:,:)
        real(dp) :: tmpval
        complex(dp) :: tmpvec(size(evecs,1))
        integer :: i, j, imax, n
        n = size(evals)
        do i = 1,n-1
            imax = i
            do j = i+1,n
                if( evals(j) > evals(imax) ) imax = j
            enddo
            if( imax /= i )then
                tmpval = evals(i)
                evals(i) = evals(imax)
                evals(imax) = tmpval
                tmpvec = evecs(:,i)
                evecs(:,i) = evecs(:,imax)
                evecs(:,imax) = tmpvec
            endif
        enddo
    end subroutine sort_eigs_desc

end module simple_complex_ppca
