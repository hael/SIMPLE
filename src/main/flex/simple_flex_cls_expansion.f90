!@descr: cls_expansion backend: per-class low-rank covariance model on Fourier coefficients (fit, basis extension, band selection, CTF means and FRCs) for the kernel sub-class averages; arrays in, arrays out
! Also the subclass budget: how many subclasses each class gets out of a total.
module simple_flex_cls_expansion
use simple_defs,               only: sp, dp, PI, DPI, CMPLX_ZERO, DCMPLX_ZERO, logfhandle
use simple_error,              only: simple_exception
use simple_rnd,                only: ran3
use simple_flex_pca_posterior, only: spd_inv_dp
implicit none

public :: flex_cls_model, flex_cls_fit, flex_cls_shell_noise, flex_cls_pose_tangents, flex_cls_weighted_mean
public :: flex_cls_ctf_mean, flex_cls_extend_basis, flex_cls_shell_frc, flex_cls_half_frc
public :: flex_cls_frc_cutoff, flex_cls_allocate_subclasses, flex_cls_split_floor
private
#include "simple_local_flags.inc"

integer,  parameter :: FLEX_CLS_N_ALS      = 4      !< prior-free probe iterations
integer,  parameter :: FLEX_CLS_N_PPCA     = 8      !< PPCA iterations (upper bound, convergence stops earlier)
integer,  parameter :: FLEX_CLS_N_RESTARTS = 3      !< random starts of the ALS probe
real(dp), parameter :: FLEX_CLS_RIDGE_REL  = 1.d-3  !< Tikhonov ridge on the basis block, relative to G(1,1)
real(dp), parameter :: FLEX_CLS_CONV_TOL   = 1.d-3  !< relative latent change that stops the PPCA loop
real(dp), parameter :: FLEX_CLS_WIENER_EPS = 0.1d0  !< sub-average denominator ctf^2 + eps, per member
real(dp), parameter :: FLEX_CLS_NUIS_RIDGE = 1.d-6  !< ridge on the free nuisance block of the E-step, relative to max A(q,q)
real(dp), parameter :: FLEX_CLS_INV_RIDGE  = 1.d-12 !< ridge before inverting a small SPD matrix, relative to max |A|
real(dp), parameter :: FLEX_CLS_GROUP_ALPHA = 0.5d0 !< a group's subclass total grows as its population**alpha

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

    !> Members a class needs before the covariance model can split it into nsub subclasses (ncomp = neigs)
    pure integer function flex_cls_split_floor( nsub, ncomp ) result( nmin )
        integer, intent(in) :: nsub, ncomp
        nmin = max(2 * nsub, 3 * (ncomp + 1))
    end function flex_cls_split_floor

    !> Subclasses per class out of ntot: group totals by highest averages (Webster) on pop**alpha from one per class,
    !! then per class by highest averages on pop; a class of at least minsplit members holds up to pop/minpop, any
    !! other class one. ngot < ntot when the capacities run out, ngot > ntot when ntot is below the class count.
    pure subroutine flex_cls_allocate_subclasses( pops, groups, ntot, minpop, minsplit, nsub, ngot )
        integer, intent(in)  :: pops(:)    !< class populations
        integer, intent(in)  :: groups(:)  !< group label of every class (any positive integers)
        integer, intent(in)  :: ntot, minpop, minsplit
        integer, intent(out) :: nsub(size(pops)), ngot
        integer,  allocatable :: cap(:), gind(:), glab(:), gk(:), gcap(:), gtot(:)
        real(dp), allocatable :: gpow(:)
        real(dp) :: q, qbest
        integer  :: m, ng, c, g, k, ibest, iextra
        m = size(pops)
        nsub = 1
        ngot = m
        if( m < 1 ) return
        allocate(cap(m), gind(m), glab(m), source=0)
        do c = 1, m
            cap(c) = 1
            if( pops(c) >= minsplit ) cap(c) = max(1, pops(c) / max(1, minpop))
        end do
        ng = 0
        do c = 1, m
            k = 0
            if( ng > 0 ) k = findloc(glab(1:ng), groups(c), dim=1)
            if( k == 0 )then
                ng       = ng + 1
                glab(ng) = groups(c)
                k        = ng
            endif
            gind(c) = k
        end do
        allocate(gk(ng), gcap(ng), gtot(ng), source=0)
        allocate(gpow(ng), source=0.d0)
        do g = 1, ng
            gk(g)   = count(gind == g)
            gcap(g) = sum(cap, mask=(gind == g))
            gpow(g) = real(sum(pops, mask=(gind == g)), dp)**FLEX_CLS_GROUP_ALPHA
        end do
        gtot = gk
        do while( ngot < ntot )
            ibest = 0
            qbest = -1.d0
            do g = 1, ng
                if( gtot(g) >= gcap(g) ) cycle
                q = gpow(g) / (real(gtot(g),dp) + 0.5d0)
                if( q > qbest )then
                    ibest = g
                    qbest = q
                endif
            end do
            if( ibest == 0 ) exit
            gtot(ibest) = gtot(ibest) + 1
            ngot        = ngot + 1
        end do
        do g = 1, ng
            do iextra = 1, gtot(g) - gk(g)
                ibest = 0
                qbest = -1.d0
                do c = 1, m
                    if( gind(c) /= g .or. nsub(c) >= cap(c) ) cycle
                    q = real(pops(c),dp) / (real(nsub(c),dp) + 0.5d0)
                    if( q > qbest )then
                        ibest = c
                        qbest = q
                    endif
                end do
                nsub(ibest) = nsub(ibest) + 1
            end do
        end do
    end subroutine flex_cls_allocate_subclasses

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
        integer,     parameter   :: CHUNK = 2048  !< coefficients per thread block
        complex(dp), allocatable :: mu(:)
        real(dp),    allocatable :: den(:)
        real(dp) :: pow(nsh)
        integer  :: cnt(nsh), ncoeff, nptcls, i, j, ic, j0, j1
        ncoeff = size(y,1)
        nptcls = size(y,2)
        allocate(mu(ncoeff), den(ncoeff))
        ! the CTF-weighted mean with unit noise weights, blockwise
        !$omp parallel do default(shared) private(ic,j0,j1,i,j) schedule(static)
        do ic = 1, (ncoeff + CHUNK - 1) / CHUNK
            j0 = (ic - 1) * CHUNK + 1
            j1 = min(ncoeff, ic * CHUNK)
            mu(j0:j1) = DCMPLX_ZERO; den(j0:j1) = 0.d0
            do i = 1, nptcls
                do j = j0, j1
                    mu(j)  = mu(j)  + real(c(j,i),dp) * cmplx(y(j,i),kind=dp)
                    den(j) = den(j) + real(c(j,i),dp)**2 + FLEX_CLS_WIENER_EPS
                end do
            end do
            do j = j0, j1
                if( den(j) > 1.d-30 )then
                    mu(j) = mu(j) / den(j)
                else
                    mu(j) = DCMPLX_ZERO
                endif
            end do
        end do
        !$omp end parallel do
        pow = 0.d0; cnt = 0
        !$omp parallel do default(shared) private(i,j) schedule(static) reduction(+:pow,cnt)
        do i = 1, nptcls
            do j = 1, ncoeff
                if( shell(j) < 1 .or. shell(j) > nsh ) cycle
                pow(shell(j)) = pow(shell(j)) + abs(cmplx(y(j,i),kind=dp) - real(c(j,i),dp) * mu(j))**2
                cnt(shell(j)) = cnt(shell(j)) + 1
            end do
        end do
        !$omp end parallel do
        do j = 1, nsh
            if( cnt(j) > 0 )then
                s2(j) = real(max(pow(j) / real(cnt(j),dp), 1.d-30), sp)
            else
                s2(j) = 1.
            endif
        end do
        deallocate(mu, den)
    end subroutine flex_cls_shell_noise

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

    !> the CTF-corrected mean of the members, sum_i w_i c_i y_i / sum_i w_i c_i^2 per coefficient (no regularisation:
    !! the averages are filtered afterwards by the FRC of their independent halves), blockwise
    subroutine flex_cls_ctf_mean( y, c, w, avg )
        complex(sp), intent(in)  :: y(:,:)        !< (ncoeff,nptcls)
        real(sp),    intent(in)  :: c(:,:)        !< (ncoeff,nptcls)
        real(sp),    intent(in)  :: w(:,:)        !< (ncoeff,nptcls)
        complex(sp), intent(out) :: avg(:)        !< (ncoeff)
        integer,     parameter   :: CHUNK = 2048  !< coefficients per thread block
        real(dp)    :: den, wc
        complex(dp) :: num
        integer     :: ncoeff, nptcls, i, j, ic, j0, j1
        ncoeff = size(y,1)
        nptcls = size(y,2)
        if( size(avg) /= ncoeff ) THROW_HARD('flex_cls_ctf_mean: avg shape')
        !$omp parallel do default(shared) private(ic,j0,j1,i,j,wc,num,den) schedule(static)
        do ic = 1, (ncoeff + CHUNK - 1) / CHUNK
            j0 = (ic - 1) * CHUNK + 1
            j1 = min(ncoeff, ic * CHUNK)
            do j = j0, j1
                num = DCMPLX_ZERO
                den = 0.d0
                do i = 1, nptcls
                    wc  = real(w(j,i),dp)
                    num = num + wc * real(c(j,i),dp) * cmplx(y(j,i),kind=dp)
                    den = den + wc * real(c(j,i),dp)**2
                end do
                avg(j) = CMPLX_ZERO
                if( den > 1.d-30 ) avg(j) = cmplx(num / den, kind=sp)
            end do
        end do
        !$omp end parallel do
    end subroutine flex_cls_ctf_mean

    !> The basis on every coefficient: the PPCA M-step with the members' latents and nuisances fixed, x = [1, z, nu],
    !! on all, the even and the odd members; above the fit band each shell is shrunk by 2F/(1+F), F the FRC of the
    !! half solutions, and the fitted coefficients keep the fitted basis. even: the halves (default: alternately)
    subroutine flex_cls_extend_basis( model, y, c, w, shell, uext, frc_basis, even )
        type(flex_cls_model), intent(in)  :: model
        complex(sp),          intent(in)  :: y(:,:)        !< (ncoeff,nptcls)
        real(sp),             intent(in)  :: c(:,:), w(:,:) !< (ncoeff,nptcls)
        integer,              intent(in)  :: shell(:)      !< (ncoeff)
        complex(dp),          intent(out) :: uext(:,:)     !< (ncoeff,ncomp)
        real(dp),    optional, intent(out) :: frc_basis(:) !< (nsh) even/odd FRC of the columns on every shell
        logical,     optional, intent(in)  :: even(:)      !< (nptcls) the even half
        integer,     allocatable :: half(:)
        integer,  parameter :: CHUNK = 512
        real(dp),    allocatable :: cov(:,:), wsh(:), num(:), dee(:), doo(:)
        real(dp),    allocatable :: xm(:,:,:), exm(:,:,:), am(:,:), bre(:,:), bim(:,:), lpk(:,:,:), rre(:,:,:), rim(:,:,:)
        complex(dp), allocatable :: ue(:,:), uo(:,:)
        integer,     allocatable :: hidx(:,:), pos(:)
        logical,     allocatable :: l_fit(:)
        real(dp)    :: x(1+model%ncomp+model%nnuis), ex(1+model%ncomp+model%nnuis,1+model%ncomp+model%nnuis)
        real(dp)    :: a(1+model%ncomp+model%nnuis,1+model%ncomp+model%nnuis), cw, f
        complex(dp) :: b(1+model%ncomp+model%nnuis,3), th(1+model%ncomp+model%nnuis,3)
        integer     :: ncoeff, nptcls, r, nn, q, npk, nsh, nchunk, ic, j0, j1, nj, j, i, h, sh, p, k, m, nh(2)
        ncoeff = size(y,1); nptcls = size(y,2); r = model%ncomp; nn = model%nnuis; q = 1 + r + nn; npk = q * (q + 1) / 2
        if( nptcls /= model%nptcls ) THROW_HARD('flex_cls_extend_basis: member count')
        if( size(uext,1) /= ncoeff .or. size(uext,2) /= r ) THROW_HARD('flex_cls_extend_basis: uext shape')
        nsh = maxval(shell)
        allocate(half(nptcls))
        do i = 1, nptcls
            half(i) = 1 + mod(i,2)                            ! 1: even, 2: odd
            if( present(even) )then
                half(i) = 2
                if( even(i) ) half(i) = 1
            endif
        end do
        ! per half: the members' regressors x = [1, z, nu] and their second moments E[x x^T] (the latent block plus
        ! its posterior covariance), packed upper triangle by columns
        nh(1) = count(half == 1); nh(2) = nptcls - nh(1)
        allocate(hidx(maxval(nh),2), pos(nptcls), source=0)
        allocate(xm(maxval(nh),q,2), exm(maxval(nh),npk,2), source=0.d0)
        nh = 0
        do i = 1, nptcls
            nh(half(i)) = nh(half(i)) + 1
            hidx(nh(half(i)),half(i)) = i
            pos(i) = nh(half(i))
        end do
        !$omp parallel do default(shared) private(i,h,m,x,ex,cov,p,k) schedule(static)
        do i = 1, nptcls
            h = half(i)
            m = pos(i)
            allocate(cov(r,r))
            cov = model%prec(:,:,i)
            call spd_inverse(cov, r)
            x(1) = 1.d0; x(2:1+r) = model%z(i,:)
            if( nn > 0 ) x(2+r:q) = model%nu(i,:)
            do p = 1, q
                ex(:,p) = x * x(p)
            end do
            ex(2:1+r,2:1+r) = ex(2:1+r,2:1+r) + cov
            xm(m,:,h) = x
            do p = 1, q
                k = p * (p - 1) / 2
                exm(m,k+1:k+p,h) = ex(1:p,p)
            end do
            deallocate(cov)
        end do
        !$omp end parallel do
        allocate(l_fit(ncoeff), source=.false.)
        l_fit(model%fitidx) = .true.
        allocate(ue(ncoeff,r), uo(ncoeff,r), source=DCMPLX_ZERO)
        uext = DCMPLX_ZERO
        nchunk = (ncoeff + CHUNK - 1) / CHUNK
        !$omp parallel do default(shared) private(ic,j0,j1,nj,am,bre,bim,lpk,rre,rim,h,m,i,cw,j,a,b,th,p,k) schedule(dynamic)
        do ic = 1, nchunk
            j0 = (ic - 1) * CHUNK + 1
            j1 = min(ncoeff, ic * CHUNK)
            nj = j1 - j0 + 1
            allocate(lpk(nj,npk,2), rre(nj,q,2), rim(nj,q,2))
            do h = 1, 2                                       ! 1: even, 2: odd (the halves)
                ! sum_i c^2 w E[x x^T] and sum_i c w y x over the half's members, as matrix products
                allocate(am(nj,nh(h)), bre(nj,nh(h)), bim(nj,nh(h)))
                do m = 1, nh(h)
                    i = hidx(m,h)
                    do j = j0, j1
                        cw = real(c(j,i),dp) * real(w(j,i),dp)
                        am(j-j0+1,m)  = cw * real(c(j,i),dp)
                        bre(j-j0+1,m) = cw * real(real(y(j,i)),dp)
                        bim(j-j0+1,m) = cw * real(aimag(y(j,i)),dp)
                    end do
                end do
                lpk(:,:,h) = matmul(am,  exm(1:nh(h),:,h))
                rre(:,:,h) = matmul(bre, xm(1:nh(h),:,h))
                rim(:,:,h) = matmul(bim, xm(1:nh(h),:,h))
                deallocate(am, bre, bim)
            end do
            do j = j0, j1
                ! full = even + odd; solve the three systems
                call unpack_lhs(lpk(j-j0+1,:,1) + lpk(j-j0+1,:,2), a)
                b(:,1) = cmplx(rre(j-j0+1,:,1) + rre(j-j0+1,:,2), rim(j-j0+1,:,1) + rim(j-j0+1,:,2), kind=dp)
                call solve3(a, b(:,1:1), th(:,1:1))
                uext(j,:) = th(2:1+r,1)
                call unpack_lhs(lpk(j-j0+1,:,1), a); b(:,2) = cmplx(rre(j-j0+1,:,1), rim(j-j0+1,:,1), kind=dp)
                call solve3(a, b(:,2:2), th(:,2:2))
                ue(j,:) = th(2:1+r,2)
                call unpack_lhs(lpk(j-j0+1,:,2), a); b(:,3) = cmplx(rre(j-j0+1,:,2), rim(j-j0+1,:,2), kind=dp)
                call solve3(a, b(:,3:3), th(:,3:3))
                uo(j,:) = th(2:1+r,3)
            end do
            deallocate(lpk, rre, rim)
        end do
        !$omp end parallel do
        ! per-shell shrinkage of the extension from the FRC of the halves: the columns pooled above the fit band
        allocate(wsh(nsh), num(nsh), dee(nsh), doo(nsh), source=0.d0)
        do j = 1, ncoeff
            sh = shell(j)
            if( sh < 1 .or. l_fit(j) ) cycle
            num(sh) = num(sh) + sum(real(ue(j,:) * conjg(uo(j,:)), dp))
            dee(sh) = dee(sh) + sum(abs(ue(j,:))**2)
            doo(sh) = doo(sh) + sum(abs(uo(j,:))**2)
        end do
        call frc_weights(wsh)
        do j = 1, ncoeff
            sh = shell(j)
            if( sh < 1 )then
                uext(j,:) = DCMPLX_ZERO
                cycle
            endif
            if( .not. l_fit(j) ) uext(j,:) = uext(j,:) * wsh(sh)
        end do
        ! the fitted coefficients keep the fitted basis
        do j = 1, model%nfit
            uext(model%fitidx(j),:) = model%U(j,:)
        end do
        if( present(frc_basis) )then
            ! the columns' even/odd agreement on every shell, the fit band included
            num = 0.d0; dee = 0.d0; doo = 0.d0
            do j = 1, ncoeff
                sh = shell(j)
                if( sh < 1 ) cycle
                num(sh) = num(sh) + sum(real(ue(j,:) * conjg(uo(j,:)), dp))
                dee(sh) = dee(sh) + sum(abs(ue(j,:))**2)
                doo(sh) = doo(sh) + sum(abs(uo(j,:))**2)
            end do
            frc_basis = 0.d0
            do sh = 1, min(nsh, size(frc_basis))
                if( dee(sh) * doo(sh) > 0.d0 ) frc_basis(sh) = num(sh) / sqrt(dee(sh) * doo(sh))
            end do
        endif
        deallocate(xm, exm, hidx, pos, l_fit, ue, uo, wsh, num, dee, doo, half)

      contains

        !> 2F/(1+F) per shell from num/dee/doo (0 where F <= 0)
        subroutine frc_weights( wout )
            real(dp), intent(out) :: wout(:)
            integer :: k
            do k = 1, size(wout)
                wout(k) = 0.d0
                if( dee(k) * doo(k) > 0.d0 )then
                    f = num(k) / sqrt(dee(k) * doo(k))
                    if( f > 0.d0 ) wout(k) = 2.d0 * f / (1.d0 + f)
                endif
            end do
        end subroutine frc_weights

        !> the symmetric q x q system from its packed upper triangle
        subroutine unpack_lhs( pk, amat )
            real(dp), intent(in)  :: pk(:)
            real(dp), intent(out) :: amat(:,:)
            integer :: pp, kk
            do pp = 1, size(amat,1)
                kk = pp * (pp - 1) / 2
                amat(1:pp,pp) = pk(kk+1:kk+pp)
                amat(pp,1:pp-1) = pk(kk+1:kk+pp-1)
            end do
        end subroutine unpack_lhs

        !> a q x q SPD solve with one complex right-hand side (real and imaginary parts)
        subroutine solve3( amat, bvec, tvec )
            real(dp),    intent(inout) :: amat(:,:)
            complex(dp), intent(in)    :: bvec(:,:)
            complex(dp), intent(out)   :: tvec(:,:)
            real(dp) :: ainv(size(amat,1),size(amat,1))
            ainv = amat
            call spd_inverse(ainv, size(amat,1))
            tvec(:,1) = cmplx(matmul(ainv, real(bvec(:,1),dp)), matmul(ainv, aimag(bvec(:,1))), kind=dp)
        end subroutine solve3

    end subroutine flex_cls_extend_basis

    !> even/odd FRC per shell of the class's CTF-weighted mean, blockwise
    subroutine flex_cls_half_frc( y, c, w, shell, nsh, frc, even )
        complex(sp), intent(in)  :: y(:,:)          !< (ncoeff,nptcls)
        real(sp),    intent(in)  :: c(:,:), w(:,:)  !< (ncoeff,nptcls)
        integer,     intent(in)  :: shell(:), nsh
        real(dp),    intent(out) :: frc(nsh)
        logical,     intent(in)  :: even(:)         !< (nptcls) the even half
        integer,     parameter   :: CHUNK = 2048
        complex(dp), allocatable :: me(:), mo(:)
        real(dp),    allocatable :: de(:), do_(:)
        real(dp) :: num(nsh), pe(nsh), po(nsh), wc
        integer  :: ncoeff, nptcls, ic, j0, j1, i, j, sh
        ncoeff = size(y,1); nptcls = size(y,2)
        allocate(me(ncoeff), mo(ncoeff), de(ncoeff), do_(ncoeff))
        !$omp parallel do default(shared) private(ic,j0,j1,i,j,wc) schedule(static)
        do ic = 1, (ncoeff + CHUNK - 1) / CHUNK
            j0 = (ic - 1) * CHUNK + 1
            j1 = min(ncoeff, ic * CHUNK)
            me(j0:j1) = DCMPLX_ZERO; mo(j0:j1) = DCMPLX_ZERO; de(j0:j1) = 0.d0; do_(j0:j1) = 0.d0
            do i = 1, nptcls
                do j = j0, j1
                    wc = real(w(j,i),dp) * real(c(j,i),dp)
                    if( even(i) )then
                        me(j) = me(j) + wc * cmplx(y(j,i),kind=dp); de(j) = de(j) + wc * real(c(j,i),dp)
                    else
                        mo(j) = mo(j) + wc * cmplx(y(j,i),kind=dp); do_(j) = do_(j) + wc * real(c(j,i),dp)
                    endif
                end do
            end do
            do j = j0, j1
                me(j) = me(j) / max(de(j), 1.d-30); mo(j) = mo(j) / max(do_(j), 1.d-30)
            end do
        end do
        !$omp end parallel do
        num = 0.d0; pe = 0.d0; po = 0.d0
        do j = 1, ncoeff
            sh = shell(j)
            if( sh < 1 .or. sh > nsh ) cycle
            num(sh) = num(sh) + real(me(j) * conjg(mo(j)), dp)
            pe(sh)  = pe(sh)  + abs(me(j))**2
            po(sh)  = po(sh)  + abs(mo(j))**2
        end do
        frc = 0.d0
        do sh = 1, nsh
            if( pe(sh) * po(sh) > 0.d0 ) frc(sh) = num(sh) / sqrt(pe(sh) * po(sh))
        end do
        deallocate(me, mo, de, do_)
    end subroutine flex_cls_half_frc

    !> the last shell before an FRC curve first falls below thr (scanned from shell 3: the first shells hold few
    !! coefficients); the last shell when it never does
    integer function flex_cls_frc_cutoff( frc, thr ) result( k )
        real(dp), intent(in) :: frc(:), thr
        integer :: sh
        k = size(frc)
        do sh = 3, size(frc)
            if( frc(sh) < thr )then
                k = sh - 1
                return
            endif
        end do
    end function flex_cls_frc_cutoff

    ! ===== private

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
            ea(1:,i)  = m(1:ncomp)
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
