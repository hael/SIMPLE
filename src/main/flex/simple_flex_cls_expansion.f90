!@descr: cls_expansion backend: per-class low-rank covariance model on Fourier coefficients (fit, placement, weighted sub-averages, reproducibility); arrays in, arrays out
module simple_flex_cls_expansion
use simple_defs,               only: sp, dp, PI, DPI, CMPLX_ZERO, DCMPLX_ZERO, logfhandle
use simple_error,              only: simple_exception
use simple_rnd,                only: ran3
use simple_flex_pca_targets,   only: kmeans_latent_targets
use simple_flex_pca_gmm,       only: gmm_state_weights
use simple_flex_pca_posterior, only: spd_inv_dp
use simple_stat,               only: kish_ess
implicit none

public :: flex_cls_model, flex_cls_fit, flex_cls_place_states, flex_cls_restore_states, flex_cls_shell_noise, flex_cls_half_reproducibility
public :: flex_cls_signal_power, flex_cls_pose_tangents, flex_cls_weighted_mean, flex_cls_fit_crossed
private
#include "simple_local_flags.inc"

integer,  parameter :: FLEX_CLS_N_ALS      = 4      !< prior-free probe iterations
integer,  parameter :: FLEX_CLS_N_PPCA     = 8      !< PPCA iterations (upper bound, convergence stops earlier)
integer,  parameter :: FLEX_CLS_N_RESTARTS = 3      !< random starts of the ALS probe
integer,  parameter :: FLEX_CLS_NFOLD      = 5      !< cross-fitting folds (basis on 4/5 of the members)
real(dp), parameter :: FLEX_CLS_RIDGE_REL  = 1.d-3  !< Tikhonov ridge on the basis block, relative to G(1,1)
real(dp), parameter :: FLEX_CLS_CONV_TOL   = 1.d-3  !< relative latent change that stops the PPCA loop
real(dp), parameter :: FLEX_CLS_WIENER_EPS = 0.1d0  !< sub-average denominator ctf^2 + eps, per member
real(dp), parameter :: FLEX_CLS_NUIS_RIDGE = 1.d-6  !< ridge on the free nuisance block of the E-step, relative to max A(q,q)
real(dp), parameter :: FLEX_CLS_INV_RIDGE  = 1.d-12 !< ridge before inverting a small SPD matrix, relative to max |A|
real(sp), parameter :: FLEX_CLS_W_CUTOFF   = 1.e-3  !< kernel weights below this are zero
integer,  parameter :: FLEX_CLS_MIN_LEAF   = 8      !< a leaf is cut only when both halves can hold this many members

type flex_cls_model
    integer :: nfit = 0, ncomp = 0, nptcls = 0
    integer,     allocatable :: fitidx(:)   !< (nfit) indices of the fitted coefficients in the caller's arrays
    complex(dp), allocatable :: mu(:)       !< (nfit) mean
    complex(dp), allocatable :: U(:,:)      !< (nfit,ncomp) basis
    real(dp),    allocatable :: z(:,:)      !< (nptcls,ncomp) MAP latent of every member
    real(dp),    allocatable :: prec(:,:,:) !< (ncomp,ncomp,nptcls) posterior precision
    integer                  :: nnuis = 0   !< fixed nuisance columns (pose-residual tangents), fitted jointly, never placed on
    complex(dp), allocatable :: N(:,:)      !< (nfit,nnuis) orthonormalised nuisance columns
    real(dp),    allocatable :: nu(:,:)     !< (nptcls,nnuis) nuisance coefficients of every member
    real(dp),    allocatable :: resid(:)    !< (nptcls) weighted fit residual per coefficient (1 = noise level)
    real(dp),    allocatable :: cpow(:)     !< (nptcls) log of the member's weighted CTF power over the fit band
    logical :: exists = .false.
  contains
    procedure :: new  => model_new
    procedure :: kill => model_kill
end type flex_cls_model

contains

    subroutine model_new( self, nfit, ncomp, nptcls )
        class(flex_cls_model), intent(inout) :: self
        integer,               intent(in)    :: nfit, ncomp, nptcls
        call self%kill
        if( nfit < 1 .or. ncomp < 1 .or. nptcls < 1 ) THROW_HARD('invalid flex_cls_model dimensions')
        self%nfit   = nfit
        self%ncomp  = ncomp
        self%nptcls = nptcls
        allocate(self%fitidx(nfit), source=0)
        allocate(self%mu(nfit), source=DCMPLX_ZERO)
        allocate(self%U(nfit,ncomp), source=DCMPLX_ZERO)
        allocate(self%z(nptcls,ncomp), source=0.d0)
        allocate(self%prec(ncomp,ncomp,nptcls), source=0.d0)
        allocate(self%resid(nptcls), source=1.d0)
        allocate(self%cpow(nptcls), source=0.d0)
        self%exists = .true.
    end subroutine model_new

    subroutine model_kill( self )
        class(flex_cls_model), intent(inout) :: self
        if( allocated(self%fitidx) ) deallocate(self%fitidx)
        if( allocated(self%mu)     ) deallocate(self%mu)
        if( allocated(self%U)      ) deallocate(self%U)
        if( allocated(self%z)      ) deallocate(self%z)
        if( allocated(self%prec)   ) deallocate(self%prec)
        if( allocated(self%N)      ) deallocate(self%N)
        if( allocated(self%nu)     ) deallocate(self%nu)
        if( allocated(self%resid)  ) deallocate(self%resid)
        if( allocated(self%cpow)   ) deallocate(self%cpow)
        self%nfit = 0; self%ncomp = 0; self%nptcls = 0; self%nnuis = 0
        self%exists = .false.
    end subroutine model_kill

    !> fits mean and basis on the coefficients in fitmask; leaves every member's latent and precision in the model
    subroutine flex_cls_fit( model, y, c, w, wq, fitmask, ncomp, verbose, nuis )
        type(flex_cls_model),  intent(inout) :: model
        complex(sp),           intent(in)    :: y(:,:)      !< (ncoeff,nptcls) coefficients in the class frame
        real(sp),              intent(in)    :: c(:,:)      !< (ncoeff,nptcls) CTF in the class frame
        real(sp),              intent(in)    :: w(:,:)      !< (ncoeff,nptcls) inverse noise variance
        real(sp),              intent(in)    :: wq(:)       !< (ncoeff) half-plane quadrature weight
        logical,               intent(in)    :: fitmask(:)  !< (ncoeff) coefficient enters the fit
        integer,               intent(in)    :: ncomp
        logical,     optional, intent(in)    :: verbose
        complex(sp), optional, intent(in)    :: nuis(:,:)   !< (ncoeff,nnuis) nuisance columns, fitted jointly, never placed on
        complex(sp), allocatable :: yf(:,:)
        real(sp),    allocatable :: cf(:,:), wf(:,:), wqf(:)
        real(dp),    allocatable :: ea(:,:), eaa(:,:,:), zprev(:,:), ubest(:,:), zbest(:,:)
        complex(dp), allocatable :: mubest(:)
        real(dp) :: dz, znorm, sq, gsum, resid, resid_best
        integer  :: nptcls, ncoeff, nfit, i, j, q, it, nals, nppca, nit, nrs, irs
        logical  :: l_prior, l_verb
        nals  = FLEX_CLS_N_ALS
        nppca = FLEX_CLS_N_PPCA
        nrs   = FLEX_CLS_N_RESTARTS   ! random starts of the ALS probe; the smallest weighted residual continues
        l_verb = .false.
        if( present(verbose) ) l_verb = verbose
        ncoeff = size(y,1)
        nptcls = size(y,2)
        if( size(c,1) /= ncoeff .or. size(c,2) /= nptcls ) THROW_HARD('flex_cls_fit: c shape')
        if( size(w,1) /= ncoeff .or. size(w,2) /= nptcls ) THROW_HARD('flex_cls_fit: w shape')
        if( size(wq) /= ncoeff .or. size(fitmask) /= ncoeff ) THROW_HARD('flex_cls_fit: wq/fitmask shape')
        nfit = count(fitmask)
        if( nfit < ncomp + 1 ) THROW_HARD('flex_cls_fit: fewer fitted coefficients than components')
        if( nptcls < ncomp + 1 ) THROW_HARD('flex_cls_fit: fewer members than components')
        call model%new(nfit, ncomp, nptcls)
        model%fitidx = pack([(j, j=1,ncoeff)], fitmask)
        allocate(yf(nfit,nptcls), cf(nfit,nptcls), wf(nfit,nptcls), wqf(nfit))
        do i = 1, nptcls
            yf(:,i) = y(model%fitidx,i)
            cf(:,i) = c(model%fitidx,i)
            wf(:,i) = w(model%fitidx,i)
        end do
        wqf = wq(model%fitidx)
        if( present(nuis) )then
            if( size(nuis,1) /= ncoeff ) THROW_HARD('flex_cls_fit: nuis shape')
            model%nnuis = size(nuis,2)
            allocate(model%N(nfit,model%nnuis), model%nu(nptcls,model%nnuis))
            model%N  = cmplx(nuis(model%fitidx,:), kind=dp)
            model%nu = 0.d0
            call orthonormalise(model%N, wqf, nfit, model%nnuis)
        endif
        allocate(ea(0:ncomp,nptcls), eaa(0:ncomp,0:ncomp,nptcls), zprev(nptcls,ncomp), source=0.d0)
        call weighted_average(yf, cf, wf, nfit, nptcls, FLEX_CLS_WIENER_EPS, model%mu)
        ! ALS probe from nrs random starts; the start with the smallest weighted residual continues
        if( nrs > 1 .and. nals > 0 )then
            allocate(ubest(nfit,2*ncomp), zbest(nptcls,ncomp), mubest(nfit))
            resid_best = huge(1.d0)
            do irs = 1, nrs
                call init_basis(model, yf, cf, wf, wqf, nfit, nptcls, ncomp)
                do it = 1, nals
                    call estep(model, yf, cf, wf, wqf, nfit, nptcls, ncomp, .false., ea, eaa)
                    call mstep(model, denuisanced(model, yf, cf, nfit, nptcls), cf, wf, nfit, nptcls, ncomp, ea, eaa)
                    call orthonormalise(model%U, wqf, nfit, ncomp)
                end do
                call estep(model, yf, cf, wf, wqf, nfit, nptcls, ncomp, .false., ea, eaa)
                resid = weighted_residual(model, denuisanced(model, yf, cf, nfit, nptcls), cf, wf, wqf, nfit, nptcls, ncomp)
                if( l_verb ) write(logfhandle,'(A,I0,A,ES12.5)') '>>> FLEX_CLS probe start ', irs, ' weighted residual=', resid
                if( resid < resid_best )then
                    resid_best = resid
                    ubest(:,1:ncomp)         = real(model%U, dp)
                    ubest(:,ncomp+1:2*ncomp) = aimag(model%U)
                    zbest  = model%z
                    mubest = model%mu
                endif
            end do
            model%U  = cmplx(ubest(:,1:ncomp), ubest(:,ncomp+1:2*ncomp), kind=dp)
            model%z  = zbest
            model%mu = mubest
            deallocate(ubest, zbest, mubest)
            nals = 0   ! the probe is done; the loop below runs the PPCA iterations only
            ! scale the basis so the latent has unit second moment
            do q = 1, ncomp
                sq = sqrt(max(sum(model%z(:,q)**2) / real(nptcls,dp), 1.d-30))
                model%U(:,q) = model%U(:,q) * sq
            end do
        else
            call init_basis(model, yf, cf, wf, wqf, nfit, nptcls, ncomp)
        endif
        nit = nals + nppca
        do it = 1, nit
            l_prior = it > nals
            if( it == nals + 1 .and. nals > 0 )then
                ! scale the basis so the latent has unit second moment
                do q = 1, ncomp
                    sq = sqrt(max(sum(model%z(:,q)**2) / real(nptcls,dp), 1.d-30))
                    model%U(:,q) = model%U(:,q) * sq
                end do
            endif
            zprev = model%z
            call estep(model, yf, cf, wf, wqf, nfit, nptcls, ncomp, l_prior, ea, eaa)
            call mstep(model, denuisanced(model, yf, cf, nfit, nptcls), cf, wf, nfit, nptcls, ncomp, ea, eaa)
            if( .not. l_prior ) call orthonormalise(model%U, wqf, nfit, ncomp)
            dz    = sqrt(sum((model%z - zprev)**2))
            znorm = sqrt(max(sum(model%z**2), 1.d-30))
            if( l_verb )then
                gsum = 0.d0
                do q = 1, ncomp
                    gsum = gsum + sum(wqf * real(conjg(model%U(:,q))*model%U(:,q),dp))
                end do
                write(logfhandle,'(A,I3,A,L1,A,ES10.3,A,ES10.3)') '>>> FLEX_CLS iter ', it, ' prior=', l_prior, &
                    &' basis power=', gsum, ' latent change=', dz/znorm
            endif
            if( l_prior .and. dz < FLEX_CLS_CONV_TOL * znorm ) exit
        end do
        ! final E-step under the prior
        call estep(model, yf, cf, wf, wqf, nfit, nptcls, ncomp, .true., ea, eaa)
        do i = 1, nptcls
            model%resid(i) = member_residual(model, yf(:,i), cf(:,i), wf(:,i), wqf, nfit, model%z(i,:), i)
            model%cpow(i)  = log(max(sum(real(wqf,dp) * real(wf(:,i),dp) * real(cf(:,i),dp)**2) / max(sum(real(wqf,dp)),1.d-30), 1.d-30))
        end do
        deallocate(yf, cf, wf, wqf, ea, eaa, zprev)
    end subroutine flex_cls_fit

    !> cross-fitted embedding: FLEX_CLS_NFOLD folds, every member embedded with a basis fitted on the other
    !! folds and mapped into the full fit's frame; the model holds the full basis and every member's latent, precision and nuisance
    subroutine flex_cls_fit_crossed( model, y, c, w, wq, fitmask, ncomp, verbose, nuis )
        type(flex_cls_model),  intent(inout) :: model
        complex(sp),           intent(in)    :: y(:,:)
        real(sp),              intent(in)    :: c(:,:), w(:,:), wq(:)
        logical,               intent(in)    :: fitmask(:)
        integer,               intent(in)    :: ncomp
        logical,     optional, intent(in)    :: verbose
        complex(sp), optional, intent(in)    :: nuis(:,:)
        type(flex_cls_model) :: fold
        complex(sp), allocatable :: yt(:,:)
        real(sp),    allocatable :: ct(:,:), wt(:,:)
        real(dp),    allocatable :: ea(:,:), eaa(:,:,:), zfull(:,:), T(:,:), G(:,:), Bm(:,:), cov(:,:), tmp(:,:)
        integer,     allocatable :: fold_of(:)
        real(dp) :: nudummy(1)
        integer :: nptcls, ncoeff, nfold, f, i, k, ntrain, nfit, q, r
        logical :: l_verb
        l_verb = .false.
        if( present(verbose) ) l_verb = verbose
        ncoeff = size(y,1)
        nptcls = size(y,2)
        nfold  = min(FLEX_CLS_NFOLD, max(2, nptcls / (3 * (ncomp + 1))))
        allocate(fold_of(nptcls))
        do i = 1, nptcls
            fold_of(i) = 1 + mod(i - 1, nfold)
        end do
        if( present(nuis) )then
            call flex_cls_fit(model, y, c, w, wq, fitmask, ncomp, verbose=l_verb, nuis=nuis)
        else
            call flex_cls_fit(model, y, c, w, wq, fitmask, ncomp, verbose=l_verb)
        endif
        nfit = model%nfit
        allocate(zfull(nptcls,ncomp), source=model%z)
        allocate(ea(0:ncomp,1), eaa(0:ncomp,0:ncomp,1))
        allocate(T(ncomp,ncomp), G(ncomp,ncomp), Bm(ncomp,ncomp), cov(ncomp,ncomp), tmp(ncomp,ncomp))
        ! Gram matrix of the full basis in the fit-band inner product
        do q = 1, ncomp
            do r = 1, ncomp
                G(q,r) = sum(real(wq(model%fitidx),dp) * real(conjg(model%U(:,q)) * model%U(:,r)))
            end do
        end do
        call spd_inverse(G, ncomp)
        do f = 1, nfold
            ntrain = count(fold_of /= f)
            allocate(yt(ncoeff,ntrain), ct(ncoeff,ntrain), wt(ncoeff,ntrain))
            k = 0
            do i = 1, nptcls
                if( fold_of(i) == f ) cycle
                k = k + 1
                yt(:,k) = y(:,i); ct(:,k) = c(:,i); wt(:,k) = w(:,i)
            end do
            if( present(nuis) )then
                call flex_cls_fit(fold, yt, ct, wt, wq, fitmask, ncomp, verbose=.false., nuis=nuis)
            else
                call flex_cls_fit(fold, yt, ct, wt, wq, fitmask, ncomp, verbose=.false.)
            endif
            ! z_full = G^-1 (U_full^H W U_fold) z_fold
            do q = 1, ncomp
                do r = 1, ncomp
                    Bm(q,r) = sum(real(wq(model%fitidx),dp) * real(conjg(model%U(:,q)) * fold%U(:,r)))
                end do
            end do
            T = matmul(G, Bm)
            do i = 1, nptcls
                if( fold_of(i) /= f ) cycle
                if( model%nnuis > 0 )then
                    call embed_one(fold, y(:,i), c(:,i), w(:,i), wq, nfit, ncomp, model%z(i,:), model%prec(:,:,i), &
                        &model%nu(i,:), ea, eaa)
                    model%resid(i) = member_residual(fold, y(fold%fitidx,i), c(fold%fitidx,i), w(fold%fitidx,i), &
                        &wq(fold%fitidx), nfit, model%z(i,:), 0, model%nu(i,:))
                else
                    call embed_one(fold, y(:,i), c(:,i), w(:,i), wq, nfit, ncomp, model%z(i,:), model%prec(:,:,i), &
                        &nudummy, ea, eaa)
                    model%resid(i) = member_residual(fold, y(fold%fitidx,i), c(fold%fitidx,i), w(fold%fitidx,i), &
                        &wq(fold%fitidx), nfit, model%z(i,:), 0)
                endif
                model%z(i,:) = matmul(T, model%z(i,:))
                cov = model%prec(:,:,i)
                call spd_inverse(cov, ncomp)
                tmp = matmul(T, matmul(cov, transpose(T)))
                call spd_inverse(tmp, ncomp)
                model%prec(:,:,i) = tmp
            end do
            call fold%kill
            deallocate(yt, ct, wt)
        end do
        deallocate(T, G, Bm, cov, tmp)
        if( l_verb ) call flex_cls_calibration(model, zfull)
        deallocate(fold_of, zfull, ea, eaa)
    end subroutine flex_cls_fit_crossed

    !> one member's latent, precision and nuisance under a given basis
    subroutine embed_one( basis, yi, ci, wi, wq, nfit, ncomp, z, prec, nu, ea, eaa )
        type(flex_cls_model), intent(inout) :: basis
        complex(sp),          intent(in)    :: yi(:)
        real(sp),             intent(in)    :: ci(:), wi(:), wq(:)
        integer,              intent(in)    :: nfit, ncomp
        real(dp),             intent(out)   :: z(ncomp), prec(ncomp,ncomp), nu(:)
        real(dp),             intent(inout) :: ea(0:ncomp,1), eaa(0:ncomp,0:ncomp,1)
        complex(sp) :: yf(nfit,1)
        real(sp)    :: cf(nfit,1), wf(nfit,1), wqf(nfit)
        type(flex_cls_model) :: one
        yf(:,1) = yi(basis%fitidx); cf(:,1) = ci(basis%fitidx); wf(:,1) = wi(basis%fitidx); wqf = wq(basis%fitidx)
        call one%new(nfit, ncomp, 1)
        one%fitidx = basis%fitidx; one%mu = basis%mu; one%U = basis%U
        if( basis%nnuis > 0 )then
            one%nnuis = basis%nnuis
            allocate(one%N(nfit,one%nnuis), source=basis%N)
            allocate(one%nu(1,one%nnuis), source=0.d0)
        endif
        call estep(one, yf, cf, wf, wqf, nfit, 1, ncomp, .true., ea, eaa)
        z    = one%z(1,:)
        prec = one%prec(:,:,1)
        if( one%nnuis > 0 ) nu = one%nu(1,:)
        call one%kill
    end subroutine embed_one

    !> one member's weighted fit residual per coefficient (1 at the noise level); nuisance from row inu or nu
    function member_residual( basis, yf, cf, wf, wqf, nfit, z, inu, nu ) result( r )
        type(flex_cls_model), intent(in) :: basis
        integer,              intent(in) :: nfit, inu
        complex(sp),          intent(in) :: yf(nfit)
        real(sp),             intent(in) :: cf(nfit), wf(nfit), wqf(nfit)
        real(dp),             intent(in) :: z(:)
        real(dp), optional,   intent(in) :: nu(:)
        real(dp) :: r
        complex(dp) :: pred(nfit)
        integer :: q
        pred = basis%mu
        do q = 1, basis%ncomp
            pred = pred + basis%U(:,q) * z(q)
        end do
        if( basis%nnuis > 0 )then
            do q = 1, basis%nnuis
                if( present(nu) )then
                    pred = pred + basis%N(:,q) * nu(q)
                else if( inu > 0 )then
                    pred = pred + basis%N(:,q) * basis%nu(inu,q)
                endif
            end do
        endif
        pred = cmplx(yf, kind=dp) - real(cf,dp) * pred
        r = sum(real(wqf,dp) * real(wf,dp) * real(conjg(pred) * pred)) / max(sum(real(wqf,dp)), 1.d-30)
    end function member_residual

    !> in-place inverse of a small SPD matrix
    subroutine spd_inverse( A, k )
        integer,  intent(in)    :: k
        real(dp), intent(inout) :: A(k,k)
        real(dp) :: Ainv(k,k)
        integer  :: q
        A = 0.5d0 * (A + transpose(A))
        do q = 1, k
            A(q,q) = A(q,q) + FLEX_CLS_INV_RIDGE * max(maxval(abs(A)), 1.d-300)
        end do
        call spd_inv_dp(A, Ainv, k)
        A = Ainv
    end subroutine spd_inverse

    !> log: cross-fitted and in-sample latent scatter over the posterior width, per component
    subroutine flex_cls_calibration( model, zfull )
        type(flex_cls_model), intent(in) :: model
        real(dp),             intent(in) :: zfull(:,:)
        real(dp) :: var_x, var_f, pw
        integer  :: q, n, i
        n = model%nptcls
        do q = 1, model%ncomp
            var_x = sum((model%z(:,q) - sum(model%z(:,q))/n)**2) / max(n-1,1)
            var_f = sum((zfull(:,q)   - sum(zfull(:,q))/n)**2)   / max(n-1,1)
            pw    = sum(1.d0 / max([(model%prec(q,q,i), i=1,n)], 1.d-30)) / n
            write(logfhandle,'(A,I2,A,F8.2,A,F8.2)') '>>> FLEX_CLS calibration comp ', q, &
                &': cross-fit scatter/posterior width=', sqrt(var_x / max(pw,1.d-30)), &
                &'  in-sample=', sqrt(var_f / max(pw,1.d-30))
        end do
    end subroutine flex_cls_calibration

    !> exactly ncls subclasses: divisive bisection of the standardised latent (plus the log fit residual) scored by
    !! Ashman's D, a tied-covariance GMM for the labels, posterior-precision kernel weights, Kish neff
    subroutine flex_cls_place_states( model, ncls, weights, labels, neff, separation )
        type(flex_cls_model), intent(in)  :: model
        integer,              intent(in)  :: ncls
        real(sp),             intent(out) :: weights(:,:)  !< (nptcls,ncls)
        integer,              intent(out) :: labels(:)     !< (nptcls)
        real(sp),             intent(out) :: neff(:)       !< (ncls)
        real(sp), optional,   intent(out) :: separation    !< Ashman's D of the first cut (standardised latent)
        real(dp), allocatable :: tcen(:,:), wcomp(:), zstd(:,:), cen2(:,:,:), dsep_of(:), zaug(:,:), paug(:,:,:)
        integer,  allocatable :: memb(:), nmemb(:)
        real(sp), allocatable :: bw(:)
        real(dp) :: var, dfirst, rmean, rvar
        integer  :: i, s, q, nleaf, isplit, nsub, k, nd
        if( .not. model%exists ) THROW_HARD('flex_cls_place_states: no model')
        if( ncls < 2 ) THROW_HARD('flex_cls_place_states: ncls must be >= 2')
        if( size(weights,1) /= model%nptcls .or. size(weights,2) /= ncls ) THROW_HARD('flex_cls_place_states: weights shape')
        if( size(labels) /= model%nptcls .or. size(neff) /= ncls ) THROW_HARD('flex_cls_place_states: labels/neff shape')
        ! placement latent = model latent + log fit residual (precision = population scatter)
        nd = model%ncomp + 1
        allocate(zaug(model%nptcls,nd), paug(nd,nd,model%nptcls), source=0.d0)
        zaug(:,1:model%ncomp) = model%z
        zaug(:,nd) = log(max(model%resid, 1.d-30))
        call detrend(zaug(:,nd), model%cpow, model%nptcls)
        rmean = sum(zaug(:,nd)) / real(model%nptcls,dp)
        rvar  = sum((zaug(:,nd) - rmean)**2) / real(max(model%nptcls-1,1),dp)
        do i = 1, model%nptcls
            paug(1:model%ncomp,1:model%ncomp,i) = model%prec(:,:,i)
            paug(nd,nd,i) = 1.d0 / max(rvar, 1.d-30)
        end do
        allocate(tcen(nd,ncls), wcomp(nd), bw(ncls), zstd(model%nptcls,nd))
        allocate(cen2(nd,2,ncls), dsep_of(ncls), memb(model%nptcls), nmemb(ncls))
        do q = 1, nd
            var = sum((zaug(:,q) - sum(zaug(:,q))/real(model%nptcls,dp))**2) / real(max(model%nptcls-1,1),dp)
            zstd(:,q) = (zaug(:,q) - sum(zaug(:,q))/real(model%nptcls,dp)) / sqrt(max(var, 1.d-30))
        end do
        memb  = 1
        nleaf = 1
        call best_bisection(zstd, model%nptcls, nd, memb, 1, cen2(:,:,1), dsep_of(1))
        dfirst = dsep_of(1)
        do while( nleaf < ncls )
            isplit = maxloc(dsep_of(1:nleaf), dim=1)
            if( dsep_of(isplit) < 0.d0 ) exit    ! nothing left to cut (leaves too small)
            nleaf = nleaf + 1
            do i = 1, model%nptcls
                if( memb(i) /= isplit ) cycle
                if( sum((zstd(i,:) - cen2(:,2,isplit))**2) < sum((zstd(i,:) - cen2(:,1,isplit))**2) ) memb(i) = nleaf
            end do
            call best_bisection(zstd, model%nptcls, nd, memb, isplit, cen2(:,:,isplit), dsep_of(isplit))
            call best_bisection(zstd, model%nptcls, nd, memb, nleaf,  cen2(:,:,nleaf),  dsep_of(nleaf))
        end do
        if( nleaf < ncls )then
            ! too few leaves: halve the largest
            do k = nleaf + 1, ncls
                isplit = maxloc([(count(memb == s), s=1,k-1)], dim=1)
                nsub = 0
                do i = 1, model%nptcls
                    if( memb(i) == isplit )then
                        nsub = nsub + 1
                        if( mod(nsub,2) == 0 ) memb(i) = k
                    endif
                end do
            end do
            nleaf = ncls
        endif
        do s = 1, ncls
            nmemb(s) = count(memb == s)
            do q = 1, nd
                if( nmemb(s) > 0 )then
                    tcen(q,s) = sum(zaug(:,q), mask=(memb == s)) / real(nmemb(s),dp)
                else
                    tcen(q,s) = 0.d0
                endif
            end do
        end do
        write(logfhandle,'(A,I0,A,F7.3,A,*(1X,I0))') '>>> FLEX_CLS divisive split into ', ncls, &
            &' subclasses; first cut D=', dfirst, ' leaf populations:', nmemb
        if( present(separation) ) separation = real(dfirst, sp)
        wcomp = 1.d0
        weights = 0.
        labels  = 0
        neff    = 0.
        bw      = 0.
        call gmm_state_weights(zaug, model%nptcls, nd, nd, ncls, tcen, wcomp, &
            &weights, neff, bw, labels, respawn=.false.)
        ! singular GMM covariance: the divisive leaves, one-hot
        if( any(sum(weights, dim=2) < 0.5) )then
            write(logfhandle,'(A)') '>>> FLEX_CLS GMM responsibilities unavailable; hard divisive subclasses'
            weights = 0.
            do i = 1, model%nptcls
                labels(i) = memb(i)
                weights(i,memb(i)) = 1.
            end do
        endif
        call posterior_kernel_weights(zaug, paug, model%nptcls, nd, ncls, labels, weights, neff)
        deallocate(tcen, wcomp, bw, zstd, cen2, dsep_of, memb, nmemb, zaug, paug)
    end subroutine flex_cls_place_states

    !> w_is = exp(-1/2 (z_i - c_s)^T A_i (z_i - c_s)), c_s the mean latent of the members labelled s; neff = Kish
    subroutine posterior_kernel_weights( z, prec, nptcls, nd, ncls, labels, weights, neff )
        integer,  intent(in)  :: nptcls, nd, ncls, labels(:)
        real(dp), intent(in)  :: z(nptcls,nd), prec(nd,nd,nptcls)
        real(sp), intent(out) :: weights(:,:), neff(:)
        real(dp) :: cen(nd), d(nd), q
        integer  :: s, i, n
        do s = 1, ncls
            n = count(labels == s)
            if( n > 0 )then
                do i = 1, nd
                    cen(i) = sum(z(:,i), mask=(labels == s)) / real(n,dp)
                end do
            else
                cen = 0.d0
            endif
            do i = 1, nptcls
                d = z(i,:) - cen
                q = dot_product(d, matmul(prec(:,:,i), d))
                weights(i,s) = real(exp(-0.5d0 * q), sp)
            end do
            where( weights(:,s) < FLEX_CLS_W_CUTOFF ) weights(:,s) = 0.
            where( labels == s ) weights(:,s) = max(weights(:,s), 1.)
            neff(s) = real(kish_ess(real(weights(:,s),dp)), sp)
        end do
    end subroutine posterior_kernel_weights

    !> the best two-way cut of one leaf: 2-means on the full latent and on each coordinate, scored by Ashman's D; dsep = -1 when too small
    subroutine best_bisection( zstd, nptcls, ncomp, memb, leaf, cen, dsep )
        integer,  intent(in)  :: nptcls, ncomp, memb(nptcls), leaf
        real(dp), intent(in)  :: zstd(nptcls,ncomp)
        real(dp), intent(out) :: cen(ncomp,2), dsep
        real(dp), allocatable :: zsub(:,:), ccand(:,:), wcomp(:)
        real(dp) :: d
        integer  :: n, i, icand
        n = count(memb == leaf)
        cen  = 0.d0
        dsep = -1.d0
        if( n < 2 * FLEX_CLS_MIN_LEAF ) return
        allocate(zsub(n,ncomp), ccand(ncomp,2), wcomp(ncomp))
        n = 0
        do i = 1, nptcls
            if( memb(i) == leaf )then
                n = n + 1
                zsub(n,:) = zstd(i,:)
            endif
        end do
        ! candidates: the full latent (without the residual coordinate), then each single coordinate
        do icand = 0, ncomp
            wcomp = 0.d0
            if( icand == 0 )then
                wcomp = 1.d0
                if( ncomp > 1 ) wcomp(ncomp) = 0.d0
            else
                wcomp(icand) = 1.d0
            endif
            call kmeans_latent_targets(zsub, n, ncomp, 2, wcomp, ccand)
            d = partition_separation(zsub, n, ncomp, 2, wcomp, ccand)
            if( d > dsep )then
                dsep = d
                cen  = ccand
            endif
        end do
        deallocate(zsub, ccand, wcomp)
    end subroutine best_bisection

    !> Ashman's D of the two best-separated clusters along their centre-to-centre axis
    function partition_separation( z, nptcls, ncomp, ncls, wcomp, cen ) result( dsep )
        integer,  intent(in) :: nptcls, ncomp, ncls
        real(dp), intent(in) :: z(nptcls,ncomp), wcomp(ncomp), cen(ncomp,ncls)
        real(dp) :: dsep
        real(dp) :: axis(ncomp), pm(ncls), ps(ncls), d2, best, proj, d
        integer  :: memb(nptcls), cnt(ncls), i, s, t, ibest
        do i = 1, nptcls
            best = huge(1.d0); ibest = 1
            do s = 1, ncls
                d2 = sum(wcomp * (z(i,:) - cen(:,s))**2)
                if( d2 < best )then
                    best = d2; ibest = s
                endif
            end do
            memb(i) = ibest
        end do
        dsep = 0.d0
        do s = 1, ncls - 1
            do t = s + 1, ncls
                axis = wcomp * (cen(:,t) - cen(:,s))
                if( sum(axis**2) < 1.d-30 ) cycle
                axis = axis / sqrt(sum(axis**2))
                pm = 0.d0; ps = 0.d0; cnt = 0
                do i = 1, nptcls
                    if( memb(i) /= s .and. memb(i) /= t ) cycle
                    proj = sum(axis * z(i,:))
                    pm(memb(i))  = pm(memb(i)) + proj
                    ps(memb(i))  = ps(memb(i)) + proj * proj
                    cnt(memb(i)) = cnt(memb(i)) + 1
                end do
                if( cnt(s) < 2 .or. cnt(t) < 2 ) cycle
                pm(s) = pm(s) / real(cnt(s),dp); pm(t) = pm(t) / real(cnt(t),dp)
                ps(s) = max(ps(s) / real(cnt(s),dp) - pm(s)**2, 1.d-30)
                ps(t) = max(ps(t) / real(cnt(t),dp) - pm(t)**2, 1.d-30)
                d = abs(pm(s) - pm(t)) / sqrt(0.5d0 * (ps(s) + ps(t)))
                dsep = max(dsep, d)
            end do
        end do
    end function partition_separation

    !> the CTF-weighted, Wiener-regularised mean
    subroutine flex_cls_weighted_mean( y, c, w, mu )
        complex(sp), intent(in)  :: y(:,:)
        real(sp),    intent(in)  :: c(:,:), w(:,:)
        complex(dp), intent(out) :: mu(:)
        call weighted_average(y, c, w, size(y,1), size(y,2), FLEX_CLS_WIENER_EPS, mu)
    end subroutine flex_cls_weighted_mean

    !> in-plane pose tangents of a class mean on the half-plane lattice: d/dx, d/dy (phase ramps), d/dtheta (central differences)
    subroutine flex_cls_pose_tangents( mu, hidx, kidx, box, nuis )
        complex(dp), intent(in)  :: mu(:)          !< (ncoeff)
        integer,     intent(in)  :: hidx(:), kidx(:), box
        complex(sp), intent(out) :: nuis(:,:)      !< (ncoeff,3)
        integer, allocatable :: lookup(:,:)
        complex(dp) :: dh, dk, mp, mm
        integer :: ncoeff, j, hmax, kmin, kmax
        ncoeff = size(mu)
        hmax = maxval(hidx); kmin = minval(kidx); kmax = maxval(kidx)
        allocate(lookup(-hmax:hmax, kmin-1:kmax+1), source=0)
        do j = 1, ncoeff
            lookup(hidx(j),kidx(j)) = j
            if( -hidx(j) >= -hmax .and. -kidx(j) >= kmin-1 .and. -kidx(j) <= kmax+1 ) lookup(-hidx(j),-kidx(j)) = -j
        end do
        do j = 1, ncoeff
            nuis(j,1) = cmplx(cmplx(0.d0, 2.d0*DPI*real(hidx(j),dp)/real(box,dp), kind=dp) * mu(j), kind=sp)
            nuis(j,2) = cmplx(cmplx(0.d0, 2.d0*DPI*real(kidx(j),dp)/real(box,dp), kind=dp) * mu(j), kind=sp)
            mp = value_at(hidx(j)+1, kidx(j)); mm = value_at(hidx(j)-1, kidx(j)); dh = 0.5d0 * (mp - mm)
            mp = value_at(hidx(j), kidx(j)+1); mm = value_at(hidx(j), kidx(j)-1); dk = 0.5d0 * (mp - mm)
            nuis(j,3) = cmplx(-real(kidx(j),dp) * dh + real(hidx(j),dp) * dk, kind=sp)
        end do
        deallocate(lookup)
    contains
        complex(dp) function value_at( h, k )
            integer, intent(in) :: h, k
            integer :: idx
            value_at = DCMPLX_ZERO
            if( abs(h) > hmax .or. k < kmin-1 .or. k > kmax+1 ) return
            idx = lookup(h,k)
            if( idx > 0 )then
                value_at = mu(idx)
            else if( idx < 0 )then
                value_at = conjg(mu(-idx))
            endif
        end function value_at
    end subroutine flex_cls_pose_tangents

    !> per-class noise spectrum when no canonical sigma2 is available: residual power to the CTF-weighted mean, per shell
    subroutine flex_cls_shell_noise( y, c, shell, nsh, s2 )
        complex(sp), intent(in)  :: y(:,:)     !< (ncoeff,nptcls)
        real(sp),    intent(in)  :: c(:,:)     !< (ncoeff,nptcls)
        integer,     intent(in)  :: shell(:)   !< (ncoeff) shell index in 1..nsh
        integer,     intent(in)  :: nsh
        real(sp),    intent(out) :: s2(nsh)
        complex(dp), allocatable :: mu(:)
        real(sp),    allocatable :: w1(:,:)
        real(dp) :: pow(nsh)
        integer  :: cnt(nsh), ncoeff, nptcls, i, j
        ncoeff = size(y,1)
        nptcls = size(y,2)
        allocate(mu(ncoeff), w1(ncoeff,nptcls))
        w1 = 1.
        call weighted_average(y, c, w1, ncoeff, nptcls, FLEX_CLS_WIENER_EPS, mu)
        pow = 0.d0; cnt = 0
        do i = 1, nptcls
            do j = 1, ncoeff
                if( shell(j) < 1 .or. shell(j) > nsh ) cycle
                pow(shell(j)) = pow(shell(j)) + abs(cmplx(y(j,i),kind=dp) - real(c(j,i),dp) * mu(j))**2
                cnt(shell(j)) = cnt(shell(j)) + 1
            end do
        end do
        do j = 1, nsh
            if( cnt(j) > 0 )then
                s2(j) = real(max(pow(j) / real(cnt(j),dp), 1.d-30), sp)
            else
                s2(j) = 1.
            endif
        end do
        deallocate(mu, w1)
    end subroutine flex_cls_shell_noise

    !> CTF-corrected weighted sub-averages; without shell a constant eps regularises, with shell and tau2 den += 1/tau2,
    !! with shell alone unregularised (shell_den returns the mean den per shell)
    subroutine flex_cls_restore_states( y, c, w, weights, avgs, shell, tau2, shell_den )
        complex(sp),        intent(in)  :: y(:,:)         !< (ncoeff,nptcls)
        real(sp),           intent(in)  :: c(:,:)         !< (ncoeff,nptcls)
        real(sp),           intent(in)  :: w(:,:)         !< (ncoeff,nptcls)
        real(sp),           intent(in)  :: weights(:,:)   !< (nptcls,ncls)
        complex(sp),        intent(out) :: avgs(:,:)      !< (ncoeff,ncls)
        integer,  optional, intent(in)  :: shell(:)       !< (ncoeff) shell index of every coefficient
        real(dp), optional, intent(in)  :: tau2(:,:)      !< (nsh,ncls) prior signal power per shell: den += 1/tau2
        real(dp), optional, intent(out) :: shell_den(:,:) !< (nsh,ncls) mean of den over the shell (1/noise variance of the average)
        real(dp),    allocatable :: den0(:), rsum(:)
        complex(dp), allocatable :: num0(:)
        integer,     allocatable :: cnt(:)
        real(dp)    :: eps, r, wc, den
        complex(dp) :: num
        integer  :: ncoeff, nptcls, ncls, i, j, s, nsh, sh
        logical  :: l_wiener, l_shell
        l_shell  = present(shell)
        l_wiener = l_shell .and. present(tau2)
        eps = FLEX_CLS_WIENER_EPS
        if( l_shell ) eps = 0.d0
        ncoeff = size(y,1)
        nptcls = size(y,2)
        ncls   = size(weights,2)
        if( size(weights,1) /= nptcls ) THROW_HARD('flex_cls_restore_states: weights shape')
        if( size(avgs,1) /= ncoeff .or. size(avgs,2) /= ncls ) THROW_HARD('flex_cls_restore_states: avgs shape')
        allocate(num0(ncoeff), den0(ncoeff))
        nsh = 0
        if( l_shell )then
            nsh = maxval(shell)
            if( present(tau2) )      nsh = size(tau2,1)
            if( present(shell_den) ) nsh = size(shell_den,1)
        endif
        allocate(rsum(max(nsh,1)), cnt(max(nsh,1)))
        do s = 1, ncls
            !$omp parallel do default(shared) private(i,j,r,wc,num,den) schedule(static)
            do j = 1, ncoeff
                num = DCMPLX_ZERO
                den = 0.d0
                do i = 1, nptcls
                    r = real(weights(i,s),dp)
                    if( r <= 0.d0 ) cycle
                    wc  = r * real(w(j,i),dp)
                    num = num + wc * real(c(j,i),dp) * cmplx(y(j,i),kind=dp)
                    den = den + wc * (real(c(j,i),dp)**2 + eps)
                end do
                num0(j) = num; den0(j) = den
            end do
            !$omp end parallel do
            if( l_shell )then
                rsum = 0.d0; cnt = 0
                do j = 1, ncoeff
                    sh = shell(j)
                    if( sh < 1 .or. sh > nsh ) cycle
                    rsum(sh) = rsum(sh) + den0(j); cnt(sh) = cnt(sh) + 1
                end do
                if( present(shell_den) )then
                    do sh = 1, nsh
                        shell_den(sh,s) = 0.d0
                        if( cnt(sh) > 0 ) shell_den(sh,s) = rsum(sh) / real(cnt(sh),dp)
                    end do
                endif
            endif
            if( l_wiener )then
                do j = 1, ncoeff
                    sh = shell(j)
                    if( sh >= 1 .and. sh <= nsh ) den0(j) = den0(j) + 1.d0 / max(tau2(sh,s), 1.d-30)
                end do
            endif
            do j = 1, ncoeff
                if( den0(j) > 1.d-30 )then
                    avgs(j,s) = cmplx(num0(j) / den0(j), kind=sp)
                else
                    avgs(j,s) = CMPLX_ZERO
                endif
            end do
        end do
        deallocate(num0, den0, rsum, cnt)
    end subroutine flex_cls_restore_states

    !> prior signal power per shell and subclass from the even/odd sub-averages: ssnr * noise variance of the half-average
    subroutine flex_cls_signal_power( ae, ao, den_e, den_o, shell, nsh, tau2 )
        complex(sp), intent(in)  :: ae(:,:), ao(:,:)         !< (ncoeff,ncls)
        real(dp),    intent(in)  :: den_e(:,:), den_o(:,:)   !< (nsh,ncls)
        integer,     intent(in)  :: shell(:), nsh
        real(dp),    intent(out) :: tau2(nsh,size(ae,2))
        real(sp) :: frc(nsh,size(ae,2))
        real(dp) :: cc, ssnr, sig2
        integer  :: s, sh
        call flex_cls_shell_frc(ae, ao, shell, nsh, frc)
        do s = 1, size(ae,2)
            do sh = 1, nsh
                cc   = min(0.999d0, max(0.001d0, real(frc(sh,s),dp)))
                ssnr = cc / (1.d0 - cc)
                sig2 = 1.d0 / max(0.5d0 * (den_e(sh,s) + den_o(sh,s)), 1.d-30)
                tau2(sh,s) = max(ssnr * sig2, 1.d-30)
            end do
        end do
    end subroutine flex_cls_signal_power

    !> per-shell Fourier ring correlation between two sets of sub-averages
    subroutine flex_cls_shell_frc( ae, ao, shell, nsh, frc )
        complex(sp), intent(in)  :: ae(:,:), ao(:,:)   !< (ncoeff,ncls)
        integer,     intent(in)  :: shell(:), nsh
        real(sp),    intent(out) :: frc(nsh,size(ae,2))
        real(dp) :: num(nsh), pe(nsh), po(nsh)
        integer  :: j, s, sh
        do s = 1, size(ae,2)
            num = 0.d0; pe = 0.d0; po = 0.d0
            do j = 1, size(ae,1)
                sh = shell(j)
                if( sh < 1 .or. sh > nsh ) cycle
                num(sh) = num(sh) + real(cmplx(ae(j,s),kind=dp) * conjg(cmplx(ao(j,s),kind=dp)))
                pe(sh)  = pe(sh)  + real(cmplx(ae(j,s),kind=dp) * conjg(cmplx(ae(j,s),kind=dp)))
                po(sh)  = po(sh)  + real(cmplx(ao(j,s),kind=dp) * conjg(cmplx(ao(j,s),kind=dp)))
            end do
            do sh = 1, nsh
                frc(sh,s) = real(num(sh) / sqrt(max(pe(sh) * po(sh), 1.d-60)), sp)
            end do
        end do
    end subroutine flex_cls_shell_frc

    !> even/odd sub-averages and, per subclass, the best cross-half correlation of its difference image to a sibling
    subroutine flex_cls_half_reproducibility( y, c, w, wq, fitmask, weights, repro, avgs_even, avgs_odd, shell, tau2, den_even, den_odd )
        complex(sp),        intent(in)  :: y(:,:)
        real(sp),           intent(in)  :: c(:,:), w(:,:), wq(:)
        logical,            intent(in)  :: fitmask(:)
        real(sp),           intent(in)  :: weights(:,:)   !< (nptcls,ncls)
        real(sp),           intent(out) :: repro(:)       !< (ncls)
        complex(sp),        intent(out) :: avgs_even(:,:), avgs_odd(:,:)   !< (ncoeff,ncls)
        integer,  optional, intent(in)  :: shell(:)
        real(dp), optional, intent(in)  :: tau2(:,:)
        real(dp), optional, intent(out) :: den_even(:,:), den_odd(:,:)   !< (nsh,ncls) mean den per shell of each half
        real(sp), allocatable :: wh(:,:)
        real(dp), allocatable :: om(:)
        complex(dp) :: de, do_
        real(dp) :: num, dee, doo, r
        integer  :: ncoeff, nptcls, ncls, i, j, s, t
        ncoeff = size(y,1); nptcls = size(y,2); ncls = size(weights,2)
        allocate(wh(nptcls,ncls), om(ncoeff))
        wh = weights
        do i = 1, nptcls
            if( mod(i,2) == 1 ) wh(i,:) = 0.
        end do
        if( present(shell) .and. present(tau2) )then
            call flex_cls_restore_states(y, c, w, wh, avgs_even, shell=shell, tau2=tau2)
        else if( present(shell) .and. present(den_even) )then
            call flex_cls_restore_states(y, c, w, wh, avgs_even, shell=shell, shell_den=den_even)
        else
            call flex_cls_restore_states(y, c, w, wh, avgs_even)
        endif
        wh = weights
        do i = 1, nptcls
            if( mod(i,2) == 0 ) wh(i,:) = 0.
        end do
        if( present(shell) .and. present(tau2) )then
            call flex_cls_restore_states(y, c, w, wh, avgs_odd, shell=shell, tau2=tau2)
        else if( present(shell) .and. present(den_odd) )then
            call flex_cls_restore_states(y, c, w, wh, avgs_odd, shell=shell, shell_den=den_odd)
        else
            call flex_cls_restore_states(y, c, w, wh, avgs_odd)
        endif
        do j = 1, ncoeff
            om(j) = 0.d0
            if( fitmask(j) ) om(j) = real(wq(j),dp) * sum(real(w(j,:),dp)) / real(nptcls,dp)
        end do
        repro = 0.
        do s = 1, ncls
            do t = 1, ncls
                if( t == s ) cycle
                num = 0.d0; dee = 0.d0; doo = 0.d0
                do j = 1, ncoeff
                    if( om(j) <= 0.d0 ) cycle
                    de  = cmplx(avgs_even(j,s) - avgs_even(j,t), kind=dp)
                    do_ = cmplx(avgs_odd(j,s)  - avgs_odd(j,t),  kind=dp)
                    num = num + om(j) * real(de * conjg(do_))
                    dee = dee + om(j) * real(de * conjg(de))
                    doo = doo + om(j) * real(do_ * conjg(do_))
                end do
                r = num / sqrt(max(dee * doo, 1.d-60))
                repro(s) = max(repro(s), real(r, sp))
            end do
        end do
        deallocate(wh, om)
    end subroutine flex_cls_half_reproducibility

    ! ===== private

    !> remove the least-squares linear trend of v on x (in place)
    subroutine detrend( v, x, n )
        integer,  intent(in)    :: n
        real(dp), intent(inout) :: v(n)
        real(dp), intent(in)    :: x(n)
        real(dp) :: xm, vm, sxx, sxv, b
        xm  = sum(x) / real(n,dp); vm = sum(v) / real(n,dp)
        sxx = sum((x - xm)**2);    sxv = sum((x - xm) * (v - vm))
        if( sxx > 1.d-30 )then
            b = sxv / sxx
            v = v - b * (x - xm)
        endif
    end subroutine detrend

    subroutine weighted_average( yf, cf, wf, nfit, nptcls, eps, mu )
        integer,     intent(in)  :: nfit, nptcls
        complex(sp), intent(in)  :: yf(nfit,nptcls)
        real(sp),    intent(in)  :: cf(nfit,nptcls), wf(nfit,nptcls)
        real(dp),    intent(in)  :: eps
        complex(dp), intent(out) :: mu(nfit)
        real(dp) :: den, wc
        complex(dp) :: num
        integer :: i, j
        !$omp parallel do default(shared) private(i,j,num,den,wc) schedule(static)
        do j = 1, nfit
            num = DCMPLX_ZERO
            den = 0.d0
            do i = 1, nptcls
                wc  = real(wf(j,i),dp)
                num = num + wc * real(cf(j,i),dp) * cmplx(yf(j,i),kind=dp)
                den = den + wc * (real(cf(j,i),dp)**2 + eps)
            end do
            if( den > 1.d-30 )then
                mu(j) = num / den
            else
                mu(j) = DCMPLX_ZERO
            endif
        end do
        !$omp end parallel do
    end subroutine weighted_average

    !> random start: random projections of the CTF-weighted residuals, orthonormalised
    subroutine init_basis( model, yf, cf, wf, wqf, nfit, nptcls, ncomp )
        type(flex_cls_model), intent(inout) :: model
        integer,              intent(in)    :: nfit, nptcls, ncomp
        complex(sp),          intent(in)    :: yf(nfit,nptcls)
        real(sp),             intent(in)    :: cf(nfit,nptcls), wf(nfit,nptcls), wqf(nfit)
        real(dp), allocatable :: g(:,:)
        real(dp) :: u1, u2
        integer  :: i, j, q
        allocate(g(nptcls,ncomp))
        do q = 1, ncomp
            do i = 1, nptcls
                u1 = max(real(ran3(),dp), 1.d-12)
                u2 = real(ran3(),dp)
                g(i,q) = sqrt(-2.d0*log(u1)) * cos(2.d0*PI*u2)
            end do
        end do
        model%U = DCMPLX_ZERO
        !$omp parallel do default(shared) private(j,q,i) schedule(static)
        do j = 1, nfit
            do q = 1, ncomp
                do i = 1, nptcls
                    model%U(j,q) = model%U(j,q) + g(i,q) * real(wf(j,i),dp) * real(cf(j,i),dp) * &
                        &(cmplx(yf(j,i),kind=dp) - real(cf(j,i),dp) * model%mu(j))
                end do
            end do
        end do
        !$omp end parallel do
        call orthonormalise(model%U, wqf, nfit, ncomp)
        deallocate(g)
    end subroutine init_basis

    !> sum_i sum_j wq_j w_ij |y_ij - c_ij (mu_j + U_j z_i)|^2 over the fitted coefficients
    function weighted_residual( model, yf, cf, wf, wqf, nfit, nptcls, ncomp ) result( resid )
        type(flex_cls_model), intent(in) :: model
        integer,              intent(in) :: nfit, nptcls, ncomp
        complex(sp),          intent(in) :: yf(nfit,nptcls)
        real(sp),             intent(in) :: cf(nfit,nptcls), wf(nfit,nptcls), wqf(nfit)
        real(dp) :: resid
        complex(dp) :: pred
        integer :: i, j
        resid = 0.d0
        !$omp parallel do default(shared) private(i,j,pred) schedule(static) reduction(+:resid)
        do i = 1, nptcls
            do j = 1, nfit
                pred  = real(cf(j,i),dp) * (model%mu(j) + sum(model%U(j,:) * model%z(i,:)))
                resid = resid + real(wqf(j),dp) * real(wf(j,i),dp) * abs(cmplx(yf(j,i),kind=dp) - pred)**2
            end do
        end do
        !$omp end parallel do
    end function weighted_residual

    !> modified Gram-Schmidt under the real half-plane inner product <u,v> = sum_j wq_j Re(u_j* v_j)
    subroutine orthonormalise( U, wqf, nfit, ncomp )
        integer,     intent(in)    :: nfit, ncomp
        complex(dp), intent(inout) :: U(nfit,ncomp)
        real(sp),    intent(in)    :: wqf(nfit)
        real(dp) :: nrm, dotp
        integer  :: q, r
        do q = 1, ncomp
            do r = 1, q - 1
                dotp   = sum(real(wqf,dp) * real(conjg(U(:,r)) * U(:,q), dp))
                U(:,q) = U(:,q) - dotp * U(:,r)
            end do
            nrm = sqrt(sum(real(wqf,dp) * real(conjg(U(:,q)) * U(:,q), dp)))
            if( nrm > 1.d-30 )then
                U(:,q) = U(:,q) / nrm
            else
                U(:,q) = DCMPLX_ZERO
            endif
        end do
    end subroutine orthonormalise

    !> posterior of every member's latent and the moments E[a], E[a a^T] of a = [1; z] for the M-step
    subroutine estep( model, yf, cf, wf, wqf, nfit, nptcls, ncomp, l_prior, ea, eaa )
        type(flex_cls_model), intent(inout) :: model
        integer,              intent(in)    :: nfit, nptcls, ncomp
        complex(sp),          intent(in)    :: yf(nfit,nptcls)
        real(sp),             intent(in)    :: cf(nfit,nptcls), wf(nfit,nptcls), wqf(nfit)
        logical,              intent(in)    :: l_prior
        real(dp),             intent(out)   :: ea(0:ncomp,nptcls), eaa(0:ncomp,0:ncomp,nptcls)
        real(dp),    allocatable :: A(:,:), S(:,:), b(:), m(:)
        complex(dp), allocatable :: col(:,:), mui(:)
        real(dp)    :: wc, cc, ridge
        complex(dp) :: res
        integer :: i, j, q, r, nn, ntot
        logical :: ok
        nn   = model%nnuis
        ntot = ncomp + nn
        !$omp parallel do default(shared) private(i,j,q,r,A,S,b,m,col,mui,wc,cc,res,ok,ridge) schedule(dynamic,8)
        do i = 1, nptcls
            allocate(A(ntot,ntot), S(ntot,ntot), b(ntot), m(ntot), col(nfit,ntot), mui(nfit))
            col(:,1:ncomp) = model%U
            mui = model%mu
            if( nn > 0 ) col(:,ncomp+1:ntot) = model%N
            A = 0.d0
            b = 0.d0
            do j = 1, nfit
                wc  = real(wqf(j),dp) * real(wf(j,i),dp)
                cc  = real(cf(j,i),dp)
                res = cmplx(yf(j,i),kind=dp) - cc * mui(j)
                do q = 1, ntot
                    b(q) = b(q) + wc * cc * real(conjg(col(j,q)) * res, dp)
                    do r = 1, q
                        A(q,r) = A(q,r) + wc * cc * cc * real(conjg(col(j,q)) * col(j,r), dp)
                    end do
                end do
            end do
            ridge = FLEX_CLS_NUIS_RIDGE * maxval([(A(q,q), q=1,ntot)])
            do q = 1, ntot
                do r = q + 1, ntot
                    A(q,r) = A(r,q)
                end do
                if( l_prior .and. q <= ncomp ) A(q,q) = A(q,q) + 1.d0
                if( q > ncomp ) A(q,q) = A(q,q) + ridge
            end do
            call spd_solve_inverse(A, ntot, b, m, S, ok)
            if( .not. ok )then
                m = 0.d0
                S = 0.d0
                do q = 1, ntot
                    S(q,q) = 1.d0
                    A(q,q) = A(q,q) + 1.d0
                end do
            endif
            model%z(i,:)      = m(1:ncomp)
            model%prec(:,:,i) = A(1:ncomp,1:ncomp)
            if( nn > 0 ) model%nu(i,:) = m(ncomp+1:ntot)
            ea(0,i)   = 1.d0
            ea(1:,i)  = m
            eaa(0,0,i) = 1.d0
            do q = 1, ncomp
                eaa(q,0,i) = m(q)
                eaa(0,q,i) = m(q)
                do r = 1, ncomp
                    eaa(q,r,i) = S(q,r) + m(q) * m(r)
                end do
            end do
            deallocate(A, S, b, m, col, mui)
        end do
        !$omp end parallel do
    end subroutine estep

    !> the data with the nuisance part removed: y_ij - c_ij sum_k N_jk nu_ik (identity without nuisance)
    function denuisanced( model, yf, cf, nfit, nptcls ) result( yr )
        type(flex_cls_model), intent(in) :: model
        integer,              intent(in) :: nfit, nptcls
        complex(sp),          intent(in) :: yf(nfit,nptcls)
        real(sp),             intent(in) :: cf(nfit,nptcls)
        complex(sp) :: yr(nfit,nptcls)
        integer :: i, j
        if( model%nnuis < 1 )then
            yr = yf
            return
        endif
        !$omp parallel do default(shared) private(i,j) schedule(static)
        do i = 1, nptcls
            do j = 1, nfit
                yr(j,i) = yf(j,i) - cmplx(real(cf(j,i),dp) * sum(model%N(j,:) * model%nu(i,:)), kind=sp)
            end do
        end do
        !$omp end parallel do
    end function denuisanced

    !> per coefficient: [mu_j; U_j(:)] = (G_j + ridge)^-1 h_j
    subroutine mstep( model, yf, cf, wf, nfit, nptcls, ncomp, ea, eaa )
        type(flex_cls_model), intent(inout) :: model
        integer,              intent(in)    :: nfit, nptcls, ncomp
        complex(sp),          intent(in)    :: yf(nfit,nptcls)
        real(sp),             intent(in)    :: cf(nfit,nptcls), wf(nfit,nptcls)
        real(dp),             intent(in)    :: ea(0:ncomp,nptcls), eaa(0:ncomp,0:ncomp,nptcls)
        real(dp)    :: G(0:ncomp,0:ncomp), Ginv(0:ncomp,0:ncomp), xr(0:ncomp), xi(0:ncomp), wc, cc, ridge
        complex(dp) :: h(0:ncomp)
        integer :: i, j, q
        logical :: ok
        !$omp parallel do default(shared) private(i,j,q,G,Ginv,h,xr,xi,wc,cc,ridge,ok) schedule(static)
        do j = 1, nfit
            G = 0.d0
            h = DCMPLX_ZERO
            do i = 1, nptcls
                cc = real(cf(j,i),dp)
                wc = real(wf(j,i),dp) * cc
                G  = G + wc * cc * eaa(:,:,i)
                h  = h + wc * cmplx(yf(j,i),kind=dp) * ea(:,i)
            end do
            ridge = FLEX_CLS_RIDGE_REL * max(G(0,0), 1.d-30)
            do q = 1, ncomp
                G(q,q) = G(q,q) + ridge
            end do
            call spd_solve_inverse(G, ncomp+1, real(h,dp), xr, Ginv, ok)
            if( ok )then
                xi = matmul(Ginv, aimag(h))
                model%mu(j)  = cmplx(xr(0), xi(0), kind=dp)
                model%U(j,:) = cmplx(xr(1:), xi(1:), kind=dp)
            else
                model%mu(j)  = DCMPLX_ZERO
                model%U(j,:) = DCMPLX_ZERO
            endif
        end do
        !$omp end parallel do
    end subroutine mstep

    !> solve A x = b and return A^-1 for a small SPD A; ok is false when the result is not finite
    subroutine spd_solve_inverse( A, n, b, x, Ainv, ok )
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: A(n,n), b(n)
        real(dp), intent(out) :: x(n), Ainv(n,n)
        logical,  intent(out) :: ok
        real(dp) :: Acopy(n,n)
        Acopy = A
        call spd_inv_dp(Acopy, Ainv, n)
        x  = matmul(Ainv, b)
        ok = all(abs(x) < huge(1.d0)) .and. all(abs(Ainv) < huge(1.d0)) .and. .not. any(x /= x)
    end subroutine spd_solve_inverse

end module simple_flex_cls_expansion
