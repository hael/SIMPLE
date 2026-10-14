!@descr: cls_expansion kernel route: every member's latent pre-image over its class, rank and bandwidth by coefficient-split cross-validation, the budget's representatives and their Stein-corrected averages
! A member's pre-image weighs the members by their likelihoods of being at its latent: ML latents of the class
! model, noise-whitened, in the leading directions of scatter over noise. Method, readouts and tests:
! doc/implementation_notes/completed/cls_expansion_kernel_preimages.md
module simple_flex_cls_kernel
use simple_defs,               only: sp, dp, DCMPLX_ZERO, logfhandle
use simple_error,              only: simple_exception
use simple_flex_cls_expansion, only: flex_cls_model, flex_cls_fit
use simple_linalg,             only: cholesky, chol_forward, jacobi, spd_inverse, spd_logdet
use simple_stat,               only: kish_ess
use simple_rnd,                only: seed_rnd_fixed
implicit none

public :: flex_cls_kernel, run_all_flex_cls_kernel_tests
private
#include "simple_local_flags.inc"

real(dp), parameter :: TAUS(8)  = [0.d0, 0.25d0, 0.5d0, 1.d0, 2.d0, 4.d0, 8.d0, -1.d0] !< bandwidths in noise units; < 0: the class mean
integer,  parameter :: RANKS(5) = [1, 2, 4, 8, 16]  !< latent ranks the cross-validation tries
integer,  parameter :: BLOCK    = 4                 !< side of the Fourier blocks of the coefficient split
integer,  parameter :: NSAMPLE  = 1000              !< validation members and candidate representatives, at most
real(dp), parameter :: RIDGE    = 1.d-6             !< ridge of the ML latent system, relative to its largest diagonal
real(dp), parameter :: WIENER   = 0.1d0             !< CV prediction denominator + this times the mean member density
real(dp), parameter :: LWMIN    = -46.d0            !< log weights this far below the largest are zero

!> One class: its members' latents in the directions the cross-validation chose, the pre-image bandwidth and the
!! representatives' weights. new -> select -> average, preimage, labels -> kill
type flex_cls_kernel
    private
    integer  :: n = 0, nk = 0, nn = 0, nsig = 0, nrep = 0, ncv = 0
    real(dp) :: tau_cv = -1.d0, tau = -1.d0, zscore = 0.d0, overlap = 0.d0
    logical  :: l_pre  = .false.                         !< the representatives are pre-images (else the class mean)
    logical  :: exists = .false.
    real(dp),    allocatable :: cv(:), cv_rank(:)        !< CV error per bandwidth at the chosen rank; best per rank
    logical,     allocatable :: info(:)                  !< (n) the members with a latent
    real(dp),    allocatable :: cov(:,:,:), snz(:,:,:)   !< full-model latent covariance and nuisance-latent block
    real(dp),    allocatable :: wsig(:,:)                !< (nk,nsig) the chosen directions, zr = z wsig
    real(dp),    allocatable :: zr(:,:), covr(:,:,:)     !< latents and covariances in them
    real(dp),    allocatable :: covc(:,:)                !< a pre-image centre's covariance, the members' mean
    integer,     allocatable :: fitidx(:)                !< the model's band, basis and nuisance columns (the Stein term)
    complex(dp), allocatable :: basis(:,:), nuisance(:,:)
    integer,     allocatable :: imem(:), lab(:)          !< (nrep) whose pre-image; (n) hard labels
    real(dp),    allocatable :: kw(:,:)                  !< (n,nrep) weights, a representative's largest 1
  contains
    procedure :: new
    procedure :: select
    procedure :: average
    procedure :: preimage
    procedure :: labels
    procedure :: neff
    procedure :: get_nrep
    procedure :: get_member
    procedure :: get_zscore
    procedure :: report
    procedure :: kill
    ! bound for the self-test submodule (external linkage)
    procedure, private :: cross_validate
    procedure, private :: logweights
    procedure, private :: kernel_mean
    procedure, private, nopass :: ml_latents
    procedure, private, nopass :: coeff_sets
end type flex_cls_kernel

interface
    !> white-box tests (submodule simple_flex_cls_kernel_tester, built only with BUILD_TESTS=ON)
    module subroutine run_all_flex_cls_kernel_tests()
    end subroutine run_all_flex_cls_kernel_tests
end interface

contains

    !> the class model's latents, the two coefficient-set models (seeded seed and seed + 1) and the
    !! cross-validation of rank and bandwidth; without it (band too narrow for the split, too few latents) select
    !! gives the class mean
    subroutine new( self, model, nuis, y, c, w, wq, hidx, kidx, seed )
        class(flex_cls_kernel), intent(inout) :: self
        type(flex_cls_model),   intent(in)    :: model
        complex(sp),            intent(in)    :: nuis(:,:)       !< (ncoeff,nnuis) the model's nuisance columns
        complex(sp),            intent(in)    :: y(:,:)          !< (ncoeff,nptcls) coefficients in the class frame
        real(sp),               intent(in)    :: c(:,:), w(:,:)  !< (ncoeff,nptcls) CTF, inverse noise variance
        real(sp),               intent(in)    :: wq(:)           !< (ncoeff) half-plane quadrature weight
        integer,                intent(in)    :: hidx(:), kidx(:), seed
        type(flex_cls_model)  :: moda, modb
        real(dp), allocatable :: z(:,:), wall(:,:)
        logical,  allocatable :: fitmask(:), seta(:), setb(:)
        integer :: ncoeff, kcv
        logical :: ok
        call self%kill
        if( .not. model%exists ) THROW_HARD('flex_cls_kernel%new: no model')
        ncoeff  = size(y,1)
        self%n  = size(y,2)
        self%nk = model%ncomp
        self%nn = model%nnuis
        self%fitidx = model%fitidx
        self%basis  = model%U
        if( self%nn > 0 ) self%nuisance = model%N
        allocate(self%cv(size(TAUS)), self%cv_rank(size(RANKS)), source=-1.d0)
        allocate(z(self%n,self%nk), self%cov(self%nk,self%nk,self%n), self%snz(max(self%nn,1),self%nk,self%n))
        call ml_latents(model, y, c, w, wq, z, self%cov, self%info, self%snz)
        call signal_space(z, self%cov, self%info, wall, ok)
        allocate(fitmask(ncoeff), seta(ncoeff), setb(ncoeff), source=.false.)
        fitmask(model%fitidx) = .true.
        call coeff_sets(hidx, kidx, fitmask, seta, setb)
        ! the set models at the class rank, lower when a set holds too few coefficients for it
        kcv = min(self%nk, min(count(seta), count(setb)) / 3 - self%nn)
        if( ok .and. kcv >= 1 )then
            call seed_rnd_fixed(seed)
            call flex_cls_fit(moda, y, c, w, wq, seta, kcv, verbose=.false., nuis=nuis)
            call seed_rnd_fixed(seed + 1)
            call flex_cls_fit(modb, y, c, w, wq, setb, kcv, verbose=.false., nuis=nuis)
            call self%cross_validate(moda, modb, y, c, w, wq, seta, setb)
            call moda%kill
            call modb%kill
        endif
        if( self%nsig > 0 )then
            self%wsig = wall(:,1:self%nsig)
            call project(z, self%cov, self%wsig, self%zr, self%covr)
            self%covc = mean_cov(self%covr, self%info)
        endif
        self%exists = .true.
        call whole(self)
    end subroutine new

    !> the representatives: the pre-images of min(nsub, candidates) members, greedily by coverage of the members'
    !! likelihoods (facility location), first among those holding minpop effective members; with nsub < 2 or
    !! without a cross-validation, one: the class mean
    subroutine select( self, nsub, minpop )
        class(flex_cls_kernel), intent(inout) :: self
        integer,                intent(in)    :: nsub, minpop
        real(dp), allocatable :: prec(:,:,:), logdet(:), kc(:,:), ne(:), kn(:), lmax(:), sc(:), cover(:), gain(:)
        real(dp), allocatable :: lw(:,:)
        integer,  allocatable :: cinds(:), pick(:)
        logical,  allocatable :: okp(:), l_tier(:), l_used(:), okl(:)
        real(dp) :: lwj(self%n), d(max(self%nsig,1)), lwmax
        integer  :: n, nc, ntake, t, i, j, r
        if( .not. self%exists ) THROW_HARD('flex_cls_kernel%select: no kernel')
        call whole(self)
        if( nsub < 2 .or. self%nsig < 1 ) return
        n = self%n
        allocate(prec(self%nsig,self%nsig,n), logdet(n))
        call precisions(self%covr, self%tau, self%covc, prec, logdet, okp)
        okp   = okp .and. self%info
        cinds = pack([(j, j=1,n)], okp)
        cinds = cinds(1:size(cinds):max(1, (size(cinds) + NSAMPLE - 1) / NSAMPLE))
        nc    = size(cinds)
        ntake = min(nsub, nc)
        if( ntake < 2 ) return
        ! every candidate's pre-image: the members' likelihoods of being at its latent
        allocate(kc(n,nc), ne(nc), kn(nc), lmax(nc))
        !$omp parallel do default(shared) private(t,i,j,d,lwj,lwmax) schedule(dynamic,4)
        do t = 1, nc
            i     = cinds(t)
            lwmax = -huge(1.d0)
            do j = 1, n
                if( .not. okp(j) ) cycle
                d      = self%zr(j,:) - self%zr(i,:)
                lwj(j) = -0.5d0 * (dot_product(d, matmul(prec(:,:,j), d)) + logdet(j))
                lwmax  = max(lwmax, lwj(j))
            end do
            lmax(t) = lwmax
            kc(:,t) = 0.d0
            do j = 1, n
                if( okp(j) .and. lwj(j) - lwmax > LWMIN ) kc(j,t) = exp(lwj(j) - lwmax)
            end do
            ne(t) = kish_ess(kc(:,t))
            kn(t) = sqrt(sum(kc(:,t)**2))
        end do
        !$omp end parallel do
        ! greedy coverage of L_tj = sc(t) kc(j,t); the candidates below minpop once the others are used up
        sc = exp(max(lmax - maxval(lmax), -700.d0))
        allocate(cover(n), source=0.d0)
        allocate(gain(nc), pick(ntake))
        allocate(l_used(nc), source=.false.)
        l_tier = ne >= real(minpop,dp)
        do r = 1, ntake
            if( .not. any(l_tier .and. .not. l_used) ) l_tier = .not. l_used
            !$omp parallel do default(shared) private(t) schedule(static)
            do t = 1, nc
                gain(t) = -1.d0
                if( l_tier(t) .and. .not. l_used(t) ) gain(t) = sum(max(0.d0, sc(t) * kc(:,t) - cover))
            end do
            !$omp end parallel do
            t         = maxloc(gain, dim=1)
            pick(r)   = t
            l_used(t) = .true.
            cover     = max(cover, sc(t) * kc(:,t))
        end do
        self%nrep  = ntake
        self%l_pre = .true.
        self%imem  = cinds(pick)
        self%kw    = kc(:,pick)
        do r = 2, ntake
            do t = 1, r - 1
                self%overlap = max(self%overlap, &
                    &sum(kc(:,pick(r)) * kc(:,pick(t))) / max(kn(pick(r)) * kn(pick(t)), 1.d-300))
            end do
        end do
        ! a member's hard label: the representative of its largest likelihood
        allocate(lw(n,ntake))
        call self%logweights(transpose(self%zr(self%imem,:)), lw, okl)
        do j = 1, n
            self%lab(j) = maxloc(lw(j,:), dim=1)
        end do
    end subroutine select

    !> representative irep's kernel mean over the members in members
    subroutine average( self, irep, y, c, w, members, avg )
        class(flex_cls_kernel), intent(in)  :: self
        integer,                intent(in)  :: irep
        complex(sp),            intent(in)  :: y(:,:)          !< (ncoeff,nptcls)
        real(sp),               intent(in)  :: c(:,:), w(:,:)  !< (ncoeff,nptcls)
        logical,                intent(in)  :: members(:)      !< (nptcls)
        complex(sp),            intent(out) :: avg(:)          !< (ncoeff)
        real(dp) :: cen(max(self%nsig,1))
        if( irep < 1 .or. irep > self%nrep ) THROW_HARD('flex_cls_kernel%average: representative index')
        cen = 0.d0
        if( self%l_pre ) cen = self%zr(self%imem(irep),:)
        call self%kernel_mean(cen, self%kw(:,irep), y, c, w, members, avg, self%l_pre)
    end subroutine average

    !> the generative model: member i's pre-image over the members in members (the class mean without pre-images)
    subroutine preimage( self, i, y, c, w, members, avg )
        class(flex_cls_kernel), intent(in)  :: self
        integer,                intent(in)  :: i
        complex(sp),            intent(in)  :: y(:,:)
        real(sp),               intent(in)  :: c(:,:), w(:,:)
        logical,                intent(in)  :: members(:)
        complex(sp),            intent(out) :: avg(:)
        real(dp), allocatable :: lw(:,:), wts(:)
        logical,  allocatable :: okl(:)
        real(dp) :: cen(max(self%nsig,1)), lwmax
        if( i < 1 .or. i > self%n ) THROW_HARD('flex_cls_kernel%preimage: member index')
        allocate(wts(self%n), source=1.d0)
        cen = 0.d0
        if( self%l_pre )then
            cen = self%zr(i,:)
            allocate(lw(self%n,1))
            call self%logweights(reshape(cen, [self%nsig,1]), lw, okl)
            lwmax = maxval(lw(:,1), mask=okl)
            wts   = 0.d0
            where( okl .and. lw(:,1) - lwmax > LWMIN ) wts = exp(lw(:,1) - lwmax)
        endif
        call self%kernel_mean(cen, wts, y, c, w, members, avg, self%l_pre)
    end subroutine preimage

    !> every member's hard label
    function labels( self ) result( lab )
        class(flex_cls_kernel), intent(in) :: self
        integer, allocatable :: lab(:)
        lab = self%lab
    end function labels

    !> Kish effective population of representative irep's weights, over the members in sel when present
    real(dp) function neff( self, irep, sel )
        class(flex_cls_kernel), intent(in) :: self
        integer,                intent(in) :: irep
        logical, optional,      intent(in) :: sel(:)
        if( present(sel) )then
            neff = kish_ess(pack(self%kw(:,irep), sel))
        else
            neff = kish_ess(self%kw(:,irep))
        endif
    end function neff

    integer function get_nrep( self )
        class(flex_cls_kernel), intent(in) :: self
        get_nrep = self%nrep
    end function get_nrep

    !> the member whose pre-image representative irep is (0: the class mean)
    integer function get_member( self, irep )
        class(flex_cls_kernel), intent(in) :: self
        integer,                intent(in) :: irep
        get_member = self%imem(irep)
    end function get_member

    !> the z-score of the pre-images' CV gain over the class mean (below about 3: no variation the held-out data confirm)
    real(dp) function get_zscore( self )
        class(flex_cls_kernel), intent(in) :: self
        get_zscore = self%zscore
    end function get_zscore

    !> the cross-validation curves and the representatives, in the log
    subroutine report( self, cls_id )
        class(flex_cls_kernel), intent(in) :: self
        integer,                intent(in) :: cls_id
        if( self%nsig < 1 )then
            write(logfhandle,'(A,I8,A)') 'Cls expansion kernel: class=', cls_id, &
                &' kept whole: no cross-validation (fit band too narrow for the coefficient split, or too few latents)'
            return
        endif
        write(logfhandle,'(A,I8,A,*(1X,F8.4))') 'Cls expansion kernel CV: class=', cls_id, &
            &' whitened error per held-out coefficient by latent rank 1 2 4 8 16 (best bandwidth; -1 not tried):', self%cv_rank
        write(logfhandle,'(A,I8,A,I3,A,*(1X,F8.4))') 'Cls expansion kernel CV: class=', cls_id, ' at rank', self%nsig, &
            &', tau = 0 0.25 0.5 1 2 4 8 (noise units) and the class mean:', self%cv
        write(logfhandle,'(A,I8,A,I8,A,I3,A,F6.2,A,F6.2,A,F8.2,A,I4,A,F6.3)') 'Cls expansion kernel: class=', cls_id, &
            &' members=', self%n, ' kernel rank=', self%nsig, ' CV tau=', self%tau_cv, ' pre-image tau=', self%tau, &
            &' CV z (pre-images over the class mean)=', self%zscore, ' representatives=', self%nrep, &
            &' max overlap=', self%overlap
    end subroutine report

    subroutine kill( self )
        class(flex_cls_kernel), intent(inout) :: self
        if( allocated(self%cv)      ) deallocate(self%cv)
        if( allocated(self%cv_rank) ) deallocate(self%cv_rank)
        if( allocated(self%info)    ) deallocate(self%info)
        if( allocated(self%cov)     ) deallocate(self%cov)
        if( allocated(self%snz)     ) deallocate(self%snz)
        if( allocated(self%wsig)    ) deallocate(self%wsig)
        if( allocated(self%zr)      ) deallocate(self%zr)
        if( allocated(self%covr)    ) deallocate(self%covr)
        if( allocated(self%covc)    ) deallocate(self%covc)
        if( allocated(self%fitidx)  ) deallocate(self%fitidx)
        if( allocated(self%basis)   ) deallocate(self%basis)
        if( allocated(self%nuisance)) deallocate(self%nuisance)
        if( allocated(self%imem)    ) deallocate(self%imem)
        if( allocated(self%lab)     ) deallocate(self%lab)
        if( allocated(self%kw)      ) deallocate(self%kw)
        self%n = 0; self%nk = 0; self%nn = 0; self%nsig = 0; self%nrep = 0; self%ncv = 0
        self%tau_cv = -1.d0; self%tau = -1.d0; self%zscore = 0.d0; self%overlap = 0.d0
        self%l_pre  = .false.
        self%exists = .false.
    end subroutine kill

    ! ===== private

    !> Two-way coefficient-split CV: latents selected on one set predict the other set's coefficients, over the
    !! leading RANKS directions and the bandwidths TAUS. Sets the rank (smallest error at a finite bandwidth), the CV
    !! optimum, the best finite bandwidth (the pre-images') and the z-score of its paired gain over the class mean
    subroutine cross_validate( self, moda, modb, y, c, w, wq, seta, setb )
        class(flex_cls_kernel), intent(inout) :: self
        type(flex_cls_model),   intent(in)    :: moda, modb
        complex(sp),            intent(in)    :: y(:,:)
        real(sp),               intent(in)    :: c(:,:), w(:,:), wq(:)
        logical,                intent(in)    :: seta(:), setb(:)
        real(dp), allocatable :: za(:,:), cova(:,:,:), zb(:,:), covb(:,:,:), wa(:,:), wb(:,:), zra(:,:), covra(:,:,:)
        real(dp), allocatable :: zrb(:,:), covrb(:,:,:), erra(:,:), errb(:,:), err(:,:), errbest(:,:), emean(:), cvr(:)
        real(dp), allocatable :: g(:)
        logical,  allocatable :: infa(:), infb(:)
        integer,  allocatable :: vinds(:)
        real(dp) :: wnorm, best, gmean, gsd
        integer  :: n, nv, ntau, ir, r, itau, i
        logical  :: oka, okb
        n    = self%n
        ntau = size(TAUS)
        allocate(za(n,moda%ncomp), cova(moda%ncomp,moda%ncomp,n), zb(n,modb%ncomp), covb(modb%ncomp,modb%ncomp,n))
        call ml_latents(moda, y, c, w, wq, za, cova, infa)
        call ml_latents(modb, y, c, w, wq, zb, covb, infb)
        call signal_space(za, cova, infa, wa, oka)
        call signal_space(zb, covb, infb, wb, okb)
        if( .not. (oka .and. okb) ) return
        vinds = [(i, i=1,n,max(1, (n + NSAMPLE - 1) / NSAMPLE))]
        nv    = size(vinds)
        wnorm = real(nv,dp) * (sum(real(wq,dp), mask=seta) + sum(real(wq,dp), mask=setb))
        allocate(erra(nv,ntau), errb(nv,ntau), errbest(nv,ntau))
        best = huge(1.d0)
        do ir = 1, size(RANKS)
            r = RANKS(ir)
            if( r > min(self%nk, moda%ncomp, modb%ncomp) ) exit
            call project(za, cova, wa(:,1:r), zra, covra)
            call project(zb, covb, wb(:,1:r), zrb, covrb)
            ! the class mean does not depend on the latents: once
            call cv_errors(zra, covra, infa, y, c, w, wq, setb, vinds, erra, ir == 1)
            call cv_errors(zrb, covrb, infb, y, c, w, wq, seta, vinds, errb, ir == 1)
            if( ir == 1 ) emean = erra(:,ntau) + errb(:,ntau)
            err = erra + errb
            err(:,ntau) = emean
            cvr = sum(err, dim=1) / max(wnorm, 1.d-30)
            self%cv_rank(ir) = minval(cvr(1:ntau-1))
            if( self%cv_rank(ir) < best )then
                best      = self%cv_rank(ir)
                self%nsig = r
                self%cv   = cvr
                errbest   = err
            endif
        end do
        self%ncv    = nv
        self%tau_cv = TAUS(minloc(self%cv, dim=1))
        itau        = minloc(self%cv(1:ntau-1), dim=1)
        self%tau    = TAUS(itau)
        if( nv > 1 )then
            g     = errbest(:,ntau) - errbest(:,itau)
            gmean = sum(g) / real(nv,dp)
            gsd   = sqrt(sum((g - gmean)**2) / real(nv - 1,dp))
            if( gsd > 0.d0 ) self%zscore = gmean / (gsd / sqrt(real(nv,dp)))
        endif
    end subroutine cross_validate

    !> log N(zr_j; cen_l, covr_j + covc + tau^2 I) of every member for the locations cen, the constant dropped;
    !! ok is false for a member without a latent or a positive-definite kernel covariance
    subroutine logweights( self, cen, lw, ok )
        class(flex_cls_kernel), intent(in)  :: self
        real(dp),               intent(in)  :: cen(:,:)    !< (nsig,nloc)
        real(dp),               intent(out) :: lw(:,:)     !< (n,nloc)
        logical, allocatable,   intent(out) :: ok(:)
        real(dp), allocatable :: prec(:,:,:), logdet(:)
        real(dp) :: d(self%nsig)
        integer  :: j, l
        allocate(prec(self%nsig,self%nsig,self%n), logdet(self%n))
        call precisions(self%covr, self%tau, self%covc, prec, logdet, ok)
        ok = ok .and. self%info
        lw = 0.d0
        !$omp parallel do default(shared) private(j,l,d) schedule(static)
        do j = 1, self%n
            if( .not. ok(j) ) cycle
            do l = 1, size(cen,2)
                d       = self%zr(j,:) - cen(:,l)
                lw(j,l) = -0.5d0 * (dot_product(d, matmul(prec(:,:,j), d)) + logdet(j))
            end do
        end do
        !$omp end parallel do
    end subroutine logweights

    !> sum_j w_j c_j W_j y_j / sum_j w_j c_j^2 W_j per coefficient; with l_stein, on the fit band, Stein's correction
    !! of the selection bias, -c w (U . cov g + N . snz g), g = wsig P_j (zr_j - cen), P_j the kernel's precision
    !! (there the weights depend on the member's own noise)
    subroutine kernel_mean( self, cen, wts, y, c, w, members, avg, l_stein )
        class(flex_cls_kernel), intent(in)  :: self
        real(dp),               intent(in)  :: cen(:), wts(:)
        complex(sp),            intent(in)  :: y(:,:)
        real(sp),               intent(in)  :: c(:,:), w(:,:)
        logical,                intent(in)  :: members(:)
        complex(sp),            intent(out) :: avg(:)
        logical,                intent(in)  :: l_stein
        integer,     parameter   :: CHUNK = 2048
        real(dp),    allocatable :: wk(:), ga(:,:), gb(:,:), prec(:,:,:), logdet(:), sa(:), sb(:)
        logical,     allocatable :: okp(:)
        complex(dp), allocatable :: num(:)
        real(dp)    :: den, cw, g(self%nk)
        integer     :: ncoeff, n, nk, nn, ic, j0, j1, i, j, jj
        ncoeff = size(y,1)
        n      = self%n
        nk     = self%nk
        nn     = self%nn
        if( size(y,2) /= n .or. size(members) /= n .or. size(wts) /= n ) THROW_HARD('flex_cls_kernel: member count')
        if( size(avg) /= ncoeff ) THROW_HARD('flex_cls_kernel: average shape')
        allocate(wk(n), source=0.d0)
        where( members ) wk = wts
        allocate(num(ncoeff))
        !$omp parallel do default(shared) private(ic,j0,j1,i,j,den,cw) schedule(static)
        do ic = 1, (ncoeff + CHUNK - 1) / CHUNK
            j0 = (ic - 1) * CHUNK + 1
            j1 = min(ncoeff, ic * CHUNK)
            do j = j0, j1
                num(j) = DCMPLX_ZERO
                den    = 0.d0
                do i = 1, n
                    if( wk(i) <= 0.d0 ) cycle
                    cw     = wk(i) * real(c(j,i),dp) * real(w(j,i),dp)
                    num(j) = num(j) + cw * cmplx(y(j,i), kind=dp)
                    den    = den + cw * real(c(j,i),dp)
                end do
                if( den > 1.d-30 )then
                    num(j) = num(j) / den
                else
                    num(j) = DCMPLX_ZERO
                endif
            end do
        end do
        !$omp end parallel do
        if( l_stein )then
            allocate(prec(self%nsig,self%nsig,n), logdet(n), ga(nk,n), gb(max(nn,1),n), source=0.d0)
            call precisions(self%covr, self%tau, self%covc, prec, logdet, okp)
            do i = 1, n
                if( wk(i) <= 0.d0 .or. .not. okp(i) ) cycle
                g       = matmul(self%wsig, matmul(prec(:,:,i), self%zr(i,:) - cen))
                ga(:,i) = matmul(self%cov(:,:,i), g)
                if( nn > 0 ) gb(:,i) = matmul(self%snz(:,:,i), g)
            end do
            !$omp parallel do default(shared) private(jj,j,i,cw,den,sa,sb) schedule(static)
            do jj = 1, size(self%fitidx)
                j = self%fitidx(jj)
                allocate(sa(nk), sb(max(nn,1)), source=0.d0)
                den = 0.d0
                do i = 1, n
                    if( wk(i) <= 0.d0 ) cycle
                    cw  = wk(i) * real(c(j,i),dp)**2 * real(w(j,i),dp)
                    den = den + cw
                    sa  = sa + cw * ga(:,i)
                    if( nn > 0 ) sb = sb + cw * gb(:,i)
                end do
                if( den > 1.d-30 )then
                    num(j) = num(j) + sum(self%basis(jj,:) * sa) / den
                    if( nn > 0 ) num(j) = num(j) + sum(self%nuisance(jj,:) * sb(1:nn)) / den
                endif
                deallocate(sa, sb)
            end do
            !$omp end parallel do
        endif
        avg = cmplx(num, kind=sp)
    end subroutine kernel_mean

    !> every member's ML latent under a model (no latent prior, the nuisance profiled out), its covariance and the
    !! nuisance-latent block of the inverse; info is false where the joint system is not positive definite
    subroutine ml_latents( model, y, c, w, wq, z, cov, info, snz )
        type(flex_cls_model), intent(in)  :: model
        complex(sp),          intent(in)  :: y(:,:)
        real(sp),             intent(in)  :: c(:,:), w(:,:), wq(:)
        real(dp),             intent(out) :: z(:,:)          !< (nptcls,ncomp)
        real(dp),             intent(out) :: cov(:,:,:)      !< (ncomp,ncomp,nptcls)
        logical, allocatable, intent(out) :: info(:)         !< (nptcls)
        real(dp), optional,   intent(out) :: snz(:,:,:)      !< (nnuis,ncomp,nptcls)
        real(dp),    allocatable :: A(:,:), S(:,:), b(:)
        complex(dp), allocatable :: col(:)
        complex(dp) :: res
        real(dp)    :: wc, cc, ridge_abs
        integer     :: nptcls, nk, nt, i, jj, j, q, r
        logical     :: ok
        nptcls = size(y,2)
        nk = model%ncomp
        nt = nk + model%nnuis
        allocate(info(nptcls))
        z   = 0.d0
        cov = 0.d0
        if( present(snz) ) snz = 0.d0
        !$omp parallel do default(shared) private(i,jj,j,q,r,A,S,b,col,res,wc,cc,ridge_abs,ok) schedule(dynamic,8)
        do i = 1, nptcls
            allocate(A(nt,nt), S(nt,nt), b(nt), source=0.d0)
            allocate(col(nt))
            do jj = 1, model%nfit
                j   = model%fitidx(jj)
                wc  = real(wq(j),dp) * real(w(j,i),dp)
                cc  = real(c(j,i),dp)
                res = cmplx(y(j,i), kind=dp) - cc * model%mu(jj)
                col(1:nk) = model%U(jj,:)
                if( nt > nk ) col(nk+1:nt) = model%N(jj,:)
                do q = 1, nt
                    b(q) = b(q) + wc * cc * real(conjg(col(q)) * res, dp)
                    do r = 1, q
                        A(q,r) = A(q,r) + wc * cc * cc * real(conjg(col(q)) * col(r), dp)
                    end do
                end do
            end do
            do q = 1, nt
                A(q,q+1:nt) = A(q+1:nt,q)
            end do
            ridge_abs = RIDGE * max(maxval([(A(q,q), q=1,nt)]), 1.d-300)
            do q = 1, nt
                A(q,q) = A(q,q) + ridge_abs
            end do
            call spd_inverse(A, S, nt, ok)
            info(i) = ok
            if( ok )then
                z(i,:)     = matmul(S(1:nk,:), b)
                cov(:,:,i) = S(1:nk,1:nk)
                if( present(snz) .and. nt > nk ) snz(:,:,i) = S(nk+1:nt,1:nk)
            endif
            deallocate(A, S, b, col)
        end do
        !$omp end parallel do
    end subroutine ml_latents

    !> the fit band in two interleaved coefficient sets, BLOCK-pixel Fourier blocks in a checkerboard (blocks,
    !! because the class-frame mask correlates neighbouring coefficients)
    pure subroutine coeff_sets( hidx, kidx, fitmask, seta, setb )
        integer, intent(in)  :: hidx(:), kidx(:)
        logical, intent(in)  :: fitmask(:)
        logical, intent(out) :: seta(:), setb(:)
        integer :: j
        logical :: l_a
        do j = 1, size(fitmask)
            l_a     = modulo(floor(real(hidx(j)) / BLOCK) + floor(real(kidx(j)) / BLOCK), 2) == 0
            seta(j) = fitmask(j) .and. l_a
            setb(j) = fitmask(j) .and. .not. l_a
        end do
    end subroutine coeff_sets

    !> the latent directions by scatter over noise: generalised eigenvectors of the informative members' latent
    !! scatter against their mean covariance, normalised to whiten it, by decreasing eigenvalue
    subroutine signal_space( z, cov, info, wall, ok )
        real(dp),              intent(in)  :: z(:,:), cov(:,:,:)
        logical,               intent(in)  :: info(:)
        real(dp), allocatable, intent(out) :: wall(:,:)    !< (k,k) the directions, by column
        logical,               intent(out) :: ok
        real(dp), allocatable :: cbar(:,:), scat(:,:), l(:,:), linv(:,:), m(:,:), v(:,:), d(:), zbar(:), e(:)
        integer :: k, j, q, nrot, ninfo
        k     = size(z,2)
        ninfo = count(info)
        allocate(wall(k,k), source=0.d0)
        ok = .false.
        if( ninfo <= k + 1 ) return
        cbar = mean_cov(cov, info)
        allocate(zbar(k), scat(k,k), source=0.d0)
        do j = 1, size(z,1)
            if( info(j) ) zbar = zbar + z(j,:) / real(ninfo,dp)
        end do
        do j = 1, size(z,1)
            if( .not. info(j) ) cycle
            do q = 1, k
                scat(:,q) = scat(:,q) + (z(j,:) - zbar) * (z(j,q) - zbar(q)) / real(ninfo - 1,dp)
            end do
        end do
        ! cbar = L L^T: the scatter in whitened coordinates is L^-1 scat L^-T, the directions L^-T v
        allocate(l(k,k), linv(k,k), e(k), d(k), v(k,k))
        call cholesky(cbar, l, k, ok)
        if( .not. ok ) return
        do q = 1, k
            e    = 0.d0
            e(q) = 1.d0
            call chol_forward(l, e, linv(:,q), k)
        end do
        m = matmul(linv, matmul(scat, transpose(linv)))
        m = 0.5d0 * (m + transpose(m))
        call jacobi(m, k, k, d, v, nrot)
        do q = 1, k
            wall(:,q) = matmul(transpose(linv), v(:,k - q + 1))   ! dsyev: ascending
        end do
    end subroutine signal_space

    !> the latents and their covariances in the directions wsig
    subroutine project( z, cov, wsig, zr, covr )
        real(dp),              intent(in)  :: z(:,:), cov(:,:,:), wsig(:,:)
        real(dp), allocatable, intent(out) :: zr(:,:), covr(:,:,:)
        integer :: j
        zr = matmul(z, wsig)
        allocate(covr(size(wsig,2),size(wsig,2),size(z,1)))
        !$omp parallel do default(shared) private(j) schedule(static)
        do j = 1, size(z,1)
            covr(:,:,j) = matmul(transpose(wsig), matmul(cov(:,:,j), wsig))
        end do
        !$omp end parallel do
    end subroutine project

    !> One direction of the CV: the latents zs (selected on the other set) predict every validation member's
    !! coefficients in valset, leaving it out; err(t,itau) the weighted squared error. l_mean: also the class mean
    subroutine cv_errors( zs, covs, info, y, c, w, wq, valset, vinds, err, l_mean )
        real(dp),    intent(in)  :: zs(:,:), covs(:,:,:)
        logical,     intent(in)  :: info(:)
        complex(sp), intent(in)  :: y(:,:)
        real(sp),    intent(in)  :: c(:,:), w(:,:), wq(:)
        logical,     intent(in)  :: valset(:)
        integer,     intent(in)  :: vinds(:)
        real(dp),    intent(out) :: err(:,:)        !< (size(vinds),size(TAUS))
        logical,     intent(in)  :: l_mean
        real(dp), allocatable :: are(:,:), aim(:,:), aden(:,:), eps(:), prec(:,:,:), logdet(:), lw(:), kw(:)
        real(dp), allocatable :: nre(:), nim(:), nden(:), covbar(:,:)
        integer,  allocatable :: qidx(:)
        logical,  allocatable :: okp(:)
        real(dp)    :: d(size(zs,2)), tau, lwmax, e
        complex(dp) :: xhat
        integer     :: n, nq, qq, q, j, t, i, itau
        n    = size(zs,1)
        qidx = pack([(q, q=1,size(valset))], valset)
        nq   = size(qidx)
        allocate(are(nq,n), aim(nq,n), aden(nq,n), prec(size(zs,2),size(zs,2),n), logdet(n))
        !$omp parallel do default(shared) private(j,qq,q) schedule(static)
        do j = 1, n
            do qq = 1, nq
                q = qidx(qq)
                are(qq,j)  = real(c(q,j),dp) * real(w(q,j),dp) * real(real(y(q,j)),dp)
                aim(qq,j)  = real(c(q,j),dp) * real(w(q,j),dp) * real(aimag(y(q,j)),dp)
                aden(qq,j) = real(c(q,j),dp)**2 * real(w(q,j),dp)
            end do
        end do
        !$omp end parallel do
        eps    = WIENER * sum(aden, dim=2) / real(max(n,1),dp)
        covbar = mean_cov(covs, info)
        err    = -1.d0
        do itau = 1, size(TAUS)
            tau = TAUS(itau)
            if( tau < 0.d0 .and. .not. l_mean ) cycle
            if( tau >= 0.d0 )then
                call precisions(covs, tau, covbar, prec, logdet, okp)
                okp = okp .and. info
            endif
            !$omp parallel do default(shared) private(t,i,j,d,lw,kw,lwmax,nre,nim,nden,qq,q,xhat,e) schedule(dynamic,4)
            do t = 1, size(vinds)
                i = vinds(t)
                allocate(lw(n), kw(n), source=0.d0)
                allocate(nre(nq), nim(nq), nden(nq), source=0.d0)
                if( tau < 0.d0 )then
                    kw = 1.d0
                else
                    lwmax = -huge(1.d0)
                    do j = 1, n
                        if( .not. okp(j) .or. j == i ) cycle
                        d     = zs(i,:) - zs(j,:)
                        lw(j) = -0.5d0 * (dot_product(d, matmul(prec(:,:,j), d)) + logdet(j))
                        lwmax = max(lwmax, lw(j))
                    end do
                    do j = 1, n
                        if( okp(j) .and. j /= i .and. lw(j) - lwmax > LWMIN ) kw(j) = exp(lw(j) - lwmax)
                    end do
                endif
                kw(i) = 0.d0
                do j = 1, n
                    if( kw(j) <= 0.d0 ) cycle
                    nre  = nre  + kw(j) * are(:,j)
                    nim  = nim  + kw(j) * aim(:,j)
                    nden = nden + kw(j) * aden(:,j)
                end do
                e = 0.d0
                do qq = 1, nq
                    q    = qidx(qq)
                    xhat = cmplx(nre(qq), nim(qq), kind=dp) / (nden(qq) + eps(qq))
                    e    = e + real(wq(q),dp) * real(w(q,i),dp) * abs(cmplx(y(q,i),kind=dp) - real(c(q,i),dp) * xhat)**2
                end do
                err(t,itau) = e
                deallocate(lw, kw, nre, nim, nden)
            end do
            !$omp end parallel do
        end do
    end subroutine cv_errors

    !> one representative of unit weights, the class mean
    subroutine whole( self )
        class(flex_cls_kernel), intent(inout) :: self
        if( allocated(self%imem) ) deallocate(self%imem)
        if( allocated(self%lab)  ) deallocate(self%lab)
        if( allocated(self%kw)   ) deallocate(self%kw)
        allocate(self%imem(1), source=0)
        allocate(self%lab(self%n), source=1)
        allocate(self%kw(self%n,1), source=1.d0)
        self%nrep    = 1
        self%l_pre   = .false.
        self%overlap = 0.d0
    end subroutine whole

    !> (cov_j + add + tau^2 I)^-1 and its log determinant for every member; ok is false where not positive definite
    subroutine precisions( cov, tau, add, prec, logdet, ok )
        real(dp),             intent(in)  :: cov(:,:,:), tau, add(:,:)
        real(dp),             intent(out) :: prec(:,:,:), logdet(:)
        logical, allocatable, intent(out) :: ok(:)
        real(dp) :: m(size(cov,1),size(cov,1))
        integer  :: j, p, k
        k = size(cov,1)
        allocate(ok(size(cov,3)))
        !$omp parallel do default(shared) private(j,m,p) schedule(static)
        do j = 1, size(cov,3)
            m = cov(:,:,j) + add
            do p = 1, k
                m(p,p) = m(p,p) + max(tau, 0.d0)**2
            end do
            call spd_inverse(m, prec(:,:,j), k, ok(j))
            if( ok(j) ) call spd_logdet(m, k, logdet(j), ok(j))
        end do
        !$omp end parallel do
    end subroutine precisions

    !> the mean covariance of the members in sel
    pure function mean_cov( covs, sel ) result( cbar )
        real(dp), intent(in) :: covs(:,:,:)
        logical,  intent(in) :: sel(:)
        real(dp) :: cbar(size(covs,1),size(covs,2))
        integer  :: j
        cbar = 0.d0
        do j = 1, size(covs,3)
            if( sel(j) ) cbar = cbar + covs(:,:,j)
        end do
        cbar = cbar / real(max(count(sel),1),dp)
    end function mean_cov

end module simple_flex_cls_kernel
