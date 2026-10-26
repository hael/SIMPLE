!@descr: trajectory_denoise pca_mode=flex: every frame the kernel pre-image of the frames that look like it in the latent of a flex class model and lie near it in time, rank and bandwidths by coefficient-split cross-validation; arrays in, arrays out
! Frames as given (no alignment, no interpolation). Without CTF and with one noise spectrum every frame's ML latent
! has the same covariance, so the kernel is Euclidean in noise-whitened latents and the CV and the pre-images are
! matrix products. The pre-images use the weights the CV scored: the fit band's two coefficient sets each averaged
! with weights from the other set's model (no selection bias), every other coefficient with the mean of the two.
! Design, toy and validation: doc/implementation_notes/planned/trajectory_denoise_flex.md
module simple_flex_traj_denoise
use simple_defs,               only: sp, dp, DCMPLX_ZERO, logfhandle
use simple_error,              only: simple_exception
use simple_flex_cls_expansion, only: flex_cls_model, flex_cls_fit, flex_cls_weighted_mean
use simple_linalg,             only: cholesky, chol_forward, jacobi, spd_inverse
use simple_stat,               only: kish_ess, median
use simple_rnd,                only: seed_rnd_fixed
implicit none

public :: flex_traj_denoise, flex_traj_noise_spectrum, flex_traj_fit_band
private
#include "simple_local_flags.inc"

real(dp), parameter :: TAUS(10)  = [0.d0, 0.25d0, 0.5d0, 1.d0, 2.d0, 4.d0, 8.d0, 16.d0, 32.d0, -1.d0] !< noise units; < 0: none
integer,  parameter :: RANKS(8)  = [1, 2, 4, 8, 16, 32, 64, 128]   !< latent ranks tried, up to the models'
real(dp), parameter :: HS(8)     = [2.d0, 4.d0, 8.d0, 16.d0, 32.d0, 64.d0, 128.d0, -1.d0] !< frames; < 0: none
integer,  parameter :: BLOCK     = 4        !< side of the Fourier blocks of the coefficient split
integer,  parameter :: BLOCKSZ   = 10       !< frames per block when the sets are windows
integer,  parameter :: NVMAX     = 1000     !< validation frames per set, at most
integer,  parameter :: NQMAX     = 1000     !< held-out coefficients per coefficient set, at most
integer,  parameter :: NMIN      = 12       !< fewest frames a set needs for the models
real(dp), parameter :: RIDGE     = 1.d-6    !< ridge of the ML latent system, relative to its largest diagonal
real(dp), parameter :: WIENER    = 0.1d0    !< CV prediction denominator + this
real(dp), parameter :: LWMIN     = -46.d0   !< log weights this far below the largest are zero
real(dp), parameter :: SNR_BAND  = 0.1d0    !< derived fit band: up to the last shell whose spectral SNR reaches this
integer,  parameter :: KBAND_MIN = 6        !< derived fit band: at least this many shells
integer,  parameter :: CHUNK     = 1024     !< coefficients per matrix product of the pre-images

!> a set's latents: ML latents of a model, noise-whitened (unit noise per coordinate), by decreasing scatter
type whitened_latents
    integer :: k = 0
    logical :: ok = .false.
    real(dp), allocatable :: zeta(:,:)   !< (m,k)
end type whitened_latents

contains

    !> y: (ncoeff,n) the frames' half-plane coefficients (soft-masked frames, as given). Returns every frame's
    !! pre-image yden and the Kish population of its weights, plus its CV choice (rank, tau, h; rank 0: the temporal
    !! kernel alone; tau, h < 0: none). halfwin <= 0: one set, the whole trajectory, every frame a candidate of every
    !! frame; halfwin > 0: blocks of BLOCKSZ frames, each with the frames within halfwin of its centre as its set.
    !! l_eo: a frame pools only frames of its own row parity (the even/odd halves stay independent)
    subroutine flex_traj_denoise( y, hidx, kidx, shell, wq, ncomp, kfit, halfwin, l_eo, seed, yden, ess, choice, kband )
        complex(sp), intent(in)  :: y(:,:)
        integer,     intent(in)  :: hidx(:), kidx(:), shell(:)
        real(sp),    intent(in)  :: wq(:)
        integer,     intent(in)  :: ncomp, kfit, halfwin, seed
        logical,     intent(in)  :: l_eo
        complex(sp), intent(out) :: yden(:,:)
        real(dp),    intent(out) :: ess(:)
        real(dp),    intent(out) :: choice(:,:)    !< (3,n) rank, tau, h per frame
        integer,     intent(out) :: kband
        real(sp), allocatable :: s2(:), w(:)
        logical,  allocatable :: fitmask(:), seta(:), setb(:)
        integer,  allocatable :: times(:)
        integer :: n, ncoeff, nsh, b, nb, c0, c1, ctr, lo, hi, t
        n      = size(y,2)
        ncoeff = size(y,1)
        nsh    = maxval(shell)
        if( size(hidx) /= ncoeff .or. size(kidx) /= ncoeff .or. size(shell) /= ncoeff .or. size(wq) /= ncoeff ) &
            &THROW_HARD('flex_traj_denoise: lattice sizes')
        if( size(yden,1) /= ncoeff .or. size(yden,2) /= n .or. size(ess) /= n ) THROW_HARD('flex_traj_denoise: outputs')
        if( size(choice,1) /= 3 .or. size(choice,2) /= n ) THROW_HARD('flex_traj_denoise: choice shape')
        ! the noise spectrum and the fit band of the whole trajectory
        call flex_traj_noise_spectrum(y, shell, nsh, s2)
        allocate(w(ncoeff))
        w = 1. / s2(shell)
        kband = flex_traj_fit_band(y, shell, nsh, s2, kfit)
        allocate(fitmask(ncoeff), seta(ncoeff), setb(ncoeff))
        fitmask = shell >= 1 .and. shell <= kband
        call coeff_sets(hidx, kidx, fitmask, seta, setb)
        write(logfhandle,'(A,I7,A,I4,A,I4,A,I7,A,I7,A,I7)') 'Trajectory flex: frames=', n, ' fit band shells 1-', kband, &
            &' of ', nsh, ' coefficients=', count(fitmask), ' set A=', count(seta), ' set B=', count(setb)
        write(logfhandle,'(A,*(1X,ES10.3))') 'Trajectory flex: noise spectrum (shells 1, 1/4, 1/2, 3/4, last):', &
            &s2(1), s2(max(1,nsh/4)), s2(max(1,nsh/2)), s2(max(1,(3*nsh)/4)), s2(nsh)
        call flush(logfhandle)
        times = [(t, t=1,n)]
        if( halfwin <= 0 .or. 2 * halfwin + 1 >= n )then
            call denoise_set(y, w, wq, fitmask, seta, setb, times, ncomp, l_eo, seed, times, yden, ess, choice, 0, .true.)
        else
            nb = (n + BLOCKSZ - 1) / BLOCKSZ
            do b = 1, nb
                c0  = (b - 1) * BLOCKSZ + 1
                c1  = min(n, b * BLOCKSZ)
                ctr = (c0 + c1) / 2
                lo  = max(1, ctr - halfwin)
                hi  = min(n, ctr + halfwin)
                call denoise_set(y(:,lo:hi), w, wq, fitmask, seta, setb, times(lo:hi), ncomp, l_eo, seed + 10 * b, &
                    &[(t - lo + 1, t=c0,c1)], yden(:,c0:c1), ess(c0:c1), choice(:,c0:c1), b, .false.)
            end do
        endif
    end subroutine flex_traj_denoise

    !> per shell, the median over consecutive frame pairs of |y_t - y_t+1|^2 / 2 averaged over the shell's
    !! coefficients: model-free and robust to events; it includes the signal change between consecutive frames
    subroutine flex_traj_noise_spectrum( y, shell, nsh, s2 )
        complex(sp),           intent(in)  :: y(:,:)
        integer,               intent(in)  :: shell(:), nsh
        real(sp), allocatable, intent(out) :: s2(:)
        real(dp), allocatable :: pw(:,:)
        integer,  allocatable :: cnt(:)
        integer :: n, t, j, sh
        n = size(y,2)
        allocate(s2(nsh), source=1.)
        if( n < 2 ) return
        allocate(cnt(nsh), source=0)
        do j = 1, size(shell)
            if( shell(j) >= 1 .and. shell(j) <= nsh ) cnt(shell(j)) = cnt(shell(j)) + 1
        end do
        allocate(pw(nsh,n-1), source=0.d0)
        !$omp parallel do default(shared) private(t,j,sh) schedule(static)
        do t = 1, n - 1
            do j = 1, size(shell)
                sh = shell(j)
                if( sh < 1 .or. sh > nsh ) cycle
                pw(sh,t) = pw(sh,t) + 0.5d0 * abs(cmplx(y(j,t+1), kind=dp) - cmplx(y(j,t), kind=dp))**2
            end do
        end do
        !$omp end parallel do
        do sh = 1, nsh
            if( cnt(sh) > 0 ) s2(sh) = real(max(median(pw(sh,:)) / real(cnt(sh),dp), 1.d-30), sp)
        end do
    end subroutine flex_traj_noise_spectrum

    !> kfit > 0: shells 1..kfit. Otherwise up to the last shell whose spectral SNR (mean frame power over the noise
    !! spectrum, minus one) reaches SNR_BAND, gaps allowed (the lattice fringes sit beyond shells of weak signal),
    !! at least KBAND_MIN shells
    integer function flex_traj_fit_band( y, shell, nsh, s2, kfit ) result( kband )
        complex(sp), intent(in) :: y(:,:)
        integer,     intent(in) :: shell(:), nsh, kfit
        real(sp),    intent(in) :: s2(:)
        real(dp), allocatable :: pw(:)
        integer,  allocatable :: cnt(:)
        integer :: j, sh
        if( kfit > 0 )then
            kband = min(kfit, nsh)
            return
        endif
        allocate(pw(nsh), source=0.d0)
        allocate(cnt(nsh), source=0)
        do j = 1, size(shell)
            sh = shell(j)
            if( sh < 1 .or. sh > nsh ) cycle
            pw(sh)  = pw(sh) + sum(abs(cmplx(y(j,:), kind=dp))**2) / real(size(y,2),dp)
            cnt(sh) = cnt(sh) + 1
        end do
        kband = 0
        do sh = 1, nsh
            if( cnt(sh) == 0 ) cycle
            if( pw(sh) / real(cnt(sh),dp) / real(s2(sh),dp) - 1.d0 >= SNR_BAND ) kband = sh
        end do
        kband = max(kband, min(KBAND_MIN, nsh))
    end function flex_traj_fit_band

    ! ===== private

    !> one set of frames (the columns of y; their frame indices in times), the candidates of its targets: a model
    !! on each coefficient set of the fit band with the set mean as the one nuisance column (contrast; the pose stays
    !! in the latent), noise-whitened latents, the CV over rank x tau x h, the targets' pre-images
    subroutine denoise_set( y, w, wq, fitmask, seta, setb, times, ncomp, l_eo, seed, targets, yden, ess, choice, iset, l_verbose )
        complex(sp), intent(in)  :: y(:,:)
        real(sp),    intent(in)  :: w(:), wq(:)
        logical,     intent(in)  :: fitmask(:), seta(:), setb(:)
        integer,     intent(in)  :: times(:), ncomp, seed, targets(:), iset
        logical,     intent(in)  :: l_eo, l_verbose
        complex(sp), intent(out) :: yden(:,:)
        real(dp),    intent(out) :: ess(:)
        real(dp),    intent(out) :: choice(:,:)
        type(flex_cls_model)     :: moda, modb
        type(whitened_latents)   :: lata, latb
        real(sp),    allocatable :: cm(:,:), wm(:,:)
        complex(sp), allocatable :: nuis(:,:)
        complex(dp), allocatable :: mu(:)
        logical,     allocatable :: even(:)
        real(dp)    :: cv_best, cv_mean, zscore, tau, h
        integer     :: m, ncoeff, kset, r, i
        logical     :: l_lat
        m      = size(y,2)
        ncoeff = size(y,1)
        allocate(even(m))
        if( l_eo )then
            even = mod(times, 2) == 0
        else
            even = .true.
        endif
        if( m < NMIN )then
            ! too few frames for any model: every target is its own frame
            do i = 1, size(targets)
                yden(:,i) = y(:,targets(i))
            end do
            ess    = 1.d0
            choice = 0.d0
            return
        endif
        allocate(cm(ncoeff,m), source=1.)
        allocate(wm(ncoeff,m))
        do i = 1, m
            wm(:,i) = w
        end do
        allocate(mu(ncoeff), nuis(ncoeff,1))
        call flex_cls_weighted_mean(y, cm, wm, mu)
        nuis(:,1) = cmplx(mu, kind=sp)
        kset  = min(ncomp, m - 2, min(count(seta), count(setb)) / 3 - 1)
        l_lat = kset >= 1
        if( l_lat )then
            call seed_rnd_fixed(seed)
            call flex_cls_fit(moda, y, cm, wm, wq, seta, kset, verbose=.false., nuis=nuis)
            call seed_rnd_fixed(seed + 1)
            call flex_cls_fit(modb, y, cm, wm, wq, setb, kset, verbose=.false., nuis=nuis)
            call latents(moda, y, w, wq, lata)
            call latents(modb, y, w, wq, latb)
            l_lat = lata%ok .and. latb%ok
            call moda%kill
            call modb%kill
        endif
        deallocate(cm, wm)
        call cross_validate(y, w, wq, fitmask, seta, setb, times, even, l_lat, lata, latb, kset, &
            &r, tau, h, cv_best, cv_mean, zscore, l_verbose)
        call preimages(y, seta, setb, times, even, targets, l_lat, lata, latb, r, tau, h, yden, ess)
        do i = 1, size(targets)
            choice(:,i) = [real(r,dp), tau, h]
        end do
        if( l_verbose .or. iset > 0 )then
            write(logfhandle,'(A,I6,A,I6,A,I6,A,I4,A,I4,A,F6.2,A,F7.1,A,F9.5,A,F9.5,A,F8.2,A,F8.1)') &
                &'Trajectory flex set=', iset, ' frames ', times(1), '-', times(m), ' model rank=', kset, &
                &' chosen rank=', r, ' tau=', tau, ' h=', h, ' CV=', cv_best, ' (parity mean ', cv_mean, ') z=', zscore, &
                &' mean neff=', sum(ess) / real(size(ess),dp)
            call flush(logfhandle)
        endif
    end subroutine denoise_set

    !> a model's ML latents (no latent prior, the nuisance profiled out), whitened by their common covariance (c and
    !! w are the same for every frame) and rotated to decreasing scatter over the set
    subroutine latents( model, y, w, wq, lat )
        type(flex_cls_model),   intent(in)    :: model
        complex(sp),            intent(in)    :: y(:,:)
        real(sp),               intent(in)    :: w(:), wq(:)
        type(whitened_latents), intent(inout) :: lat
        real(dp),    allocatable :: A(:,:), S(:,:), bvec(:,:), z(:,:), l(:,:), linv(:,:), u(:,:), scat(:,:), v(:,:), d(:)
        real(dp),    allocatable :: e(:), ubar(:)
        complex(dp), allocatable :: col(:,:)
        real(dp)    :: wc, ridge_abs
        integer     :: m, nk, nn, nt, jj, j, q, r, i, nrot
        logical     :: ok
        m  = size(y,2)
        nk = model%ncomp
        nn = model%nnuis
        nt = nk + nn
        lat%k  = nk
        lat%ok = .false.
        ! the columns on the fit band and the common normal matrix
        allocate(col(model%nfit,nt))
        col(:,1:nk) = model%U
        if( nn > 0 ) col(:,nk+1:nt) = model%N
        allocate(A(nt,nt), source=0.d0)
        do jj = 1, model%nfit
            j  = model%fitidx(jj)
            wc = real(wq(j),dp) * real(w(j),dp)
            do q = 1, nt
                do r = 1, q
                    A(q,r) = A(q,r) + wc * real(conjg(col(jj,q)) * col(jj,r), dp)
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
        allocate(S(nt,nt))
        call spd_inverse(A, S, nt, ok)
        if( .not. ok ) return
        ! the right-hand sides of every frame, then the latents
        allocate(bvec(nt,m), source=0.d0)
        !$omp parallel do default(shared) private(i,jj,j,wc,q) schedule(static)
        do i = 1, m
            do jj = 1, model%nfit
                j  = model%fitidx(jj)
                wc = real(wq(j),dp) * real(w(j),dp)
                do q = 1, nt
                    bvec(q,i) = bvec(q,i) + wc * real(conjg(col(jj,q)) * (cmplx(y(j,i), kind=dp) - model%mu(jj)), dp)
                end do
            end do
        end do
        !$omp end parallel do
        z = matmul(S(1:nk,:), bvec)                         ! (nk,m)
        ! whiten by the common covariance S_zz = L L^T, then the directions of decreasing scatter
        allocate(l(nk,nk), linv(nk,nk), e(nk))
        call cholesky(S(1:nk,1:nk), l, nk, ok)
        if( .not. ok ) return
        do q = 1, nk
            e    = 0.d0
            e(q) = 1.d0
            call chol_forward(l, e, linv(:,q), nk)
        end do
        u    = matmul(linv, z)                              ! (nk,m), unit noise
        ubar = sum(u, dim=2) / real(m,dp)
        allocate(scat(nk,nk), source=0.d0)
        do i = 1, m
            do q = 1, nk
                scat(:,q) = scat(:,q) + (u(:,i) - ubar) * (u(q,i) - ubar(q)) / real(max(m - 1,1),dp)
            end do
        end do
        scat = 0.5d0 * (scat + transpose(scat))
        allocate(v(nk,nk), d(nk))
        call jacobi(scat, nk, nk, d, v, nrot)
        v = v(:,nk:1:-1)                                    ! decreasing scatter
        lat%zeta = transpose(matmul(transpose(v), u))       ! (m,nk)
        lat%ok   = .true.
    end subroutine latents

    !> Two-way coefficient-split CV over RANKS x TAUS x HS: the latents of one set's model weigh the other frames
    !! (Euclidean in whitened latents, squared distance over 2 + tau^2, times the temporal kernel), leaving the
    !! validation frame out and pooling only its parity, and predict its coefficients of the other set; without
    !! latents the temporal kernel alone predicts every fit-band coefficient. Sets the optimum (r = 0: the temporal
    !! kernel alone), the CV error there and of the parity mean (tau, h < 0), and the z-score of the paired gain
    subroutine cross_validate( y, w, wq, fitmask, seta, setb, times, even, l_lat, lata, latb, kset, &
                               r, tau, h, cv_best, cv_mean, zscore, l_verbose )
        complex(sp),            intent(in)  :: y(:,:)
        real(sp),               intent(in)  :: w(:), wq(:)
        logical,                intent(in)  :: fitmask(:), seta(:), setb(:), even(:)
        integer,                intent(in)  :: times(:), kset
        logical,                intent(in)  :: l_lat, l_verbose
        type(whitened_latents), intent(in)  :: lata, latb
        integer,                intent(out) :: r
        real(dp),               intent(out) :: tau, h, cv_best, cv_mean, zscore
        real(dp), allocatable :: ea(:,:,:,:), eb(:,:,:,:), e(:), eref(:), ebest(:), cvh(:), cvr(:)
        integer,  allocatable :: vinds(:)
        real(dp) :: wnorm, cvv, gmean, gsd
        integer  :: m, nv, nr, nt, nh, ir, it, ih, rmax, irb, itb, ihb, i
        m  = size(y,2)
        nt = size(TAUS)
        nh = size(HS)
        vinds = [(i, i=1,m,max(1, (m + NVMAX - 1) / NVMAX))]
        nv    = size(vinds)
        nr    = 0
        if( l_lat )then
            rmax = min(lata%k, latb%k)
            do ir = 1, size(RANKS)
                if( RANKS(ir) <= rmax ) nr = ir
            end do
        endif
        ! ea(v,it,ih,ir): held-out error of validation frame v; it = nt is the temporal kernel alone, at ir = 1
        allocate(ea(nv,nt,nh,max(nr,1)), eb(nv,nt,nh,max(nr,1)), source=0.d0)
        if( nr > 0 )then
            call cv_direction(lata%zeta, y, w, wq, setb, times, even, vinds, nr, ea)
            call cv_direction(latb%zeta, y, w, wq, seta, times, even, vinds, nr, eb)
            wnorm = real(nv,dp) * (sum(real(wq,dp), mask=pick(setb)) + sum(real(wq,dp), mask=pick(seta)))
        else
            call cv_direction(reshape([0.d0], [1,1]), y, w, wq, fitmask, times, even, vinds, 0, ea)
            wnorm = real(nv,dp) * sum(real(wq,dp), mask=pick(fitmask))
        endif
        ea = ea + eb
        ! the optimum over the grid
        cv_best = huge(1.d0)
        irb = 1; itb = nt; ihb = nh
        allocate(cvh(nh), source=huge(1.d0))
        allocate(cvr(size(RANKS)), source=-1.d0)
        do ih = 1, nh
            do it = 1, nt
                do ir = 1, max(nr,1)
                    if( it == nt .and. ir > 1 ) exit
                    if( it < nt .and. nr == 0 ) exit
                    cvv = sum(ea(:,it,ih,ir)) / max(wnorm, 1.d-30)
                    cvh(ih) = min(cvh(ih), cvv)
                    if( it < nt )then
                        if( cvr(ir) < 0.d0 .or. cvv < cvr(ir) ) cvr(ir) = cvv
                    endif
                    if( cvv < cv_best )then
                        cv_best = cvv
                        irb = ir; itb = it; ihb = ih
                    endif
                end do
            end do
        end do
        cv_mean = sum(ea(:,nt,nh,1)) / max(wnorm, 1.d-30)
        r = 0
        if( itb < nt ) r = RANKS(irb)
        tau = TAUS(itb)
        h   = HS(ihb)
        zscore = 0.d0
        if( nv > 1 )then
            eref  = ea(:,nt,nh,1)
            ebest = ea(:,itb,ihb,irb)
            e     = eref - ebest
            gmean = sum(e) / real(nv,dp)
            gsd   = sqrt(sum((e - gmean)**2) / real(nv - 1,dp))
            if( gsd > 0.d0 ) zscore = gmean / (gsd / sqrt(real(nv,dp)))
        endif
        if( l_verbose )then
            write(logfhandle,'(A,I6,A,I4,A)') 'Trajectory flex CV: validation frames=', nv, ' model rank=', kset, &
                &' (whitened error per held-out coefficient; 1 is the noise floor)'
            write(logfhandle,'(A,*(1X,F9.5))') 'Trajectory flex CV by temporal bandwidth h = 2 4 8 16 32 64 128 frames, none:', cvh
            write(logfhandle,'(A,*(1X,F9.5))') 'Trajectory flex CV by latent rank 1 2 4 8 16 32 64 128 (-1 not tried):', cvr
            write(logfhandle,'(A,*(1X,F9.5))') 'Trajectory flex CV by tau = 0 .25 .5 1 2 4 8 16 32 (noise units), none, at the '//&
                &'chosen rank and h:', (sum(ea(:,it,ihb,merge(irb, 1, it < nt))) / max(wnorm, 1.d-30), it=1,nt)
            call flush(logfhandle)
        endif
    end subroutine cross_validate

    !> one direction of the CV: latents zeta (a set model's, or none when nr = 0) predict the validation frames'
    !! coefficients in valset (at most NQMAX of them) from the other frames; err(v,it,ih,ir) the weighted squared error
    subroutine cv_direction( zeta, y, w, wq, valset, times, even, vinds, nr, err )
        real(dp),    intent(in)    :: zeta(:,:)
        complex(sp), intent(in)    :: y(:,:)
        real(sp),    intent(in)    :: w(:), wq(:)
        logical,     intent(in)    :: valset(:), even(:)
        integer,     intent(in)    :: times(:), vinds(:), nr
        real(dp),    intent(inout) :: err(:,:,:,:)
        real(sp), allocatable :: rm(:,:), kmat(:,:), pred(:,:)
        real(dp), allocatable :: d2(:,:), lw(:,:), ksum(:)
        integer,  allocatable :: qidx(:)
        logical,  allocatable :: allowed(:,:)
        real(dp) :: tau, h, lwmax, sc
        integer  :: m, nv, nq, nt, nh, ir, it, ih, v, j, q, rprev, r, i
        m  = size(y,2)
        nv = size(vinds)
        nt = size(TAUS)
        nh = size(HS)
        ! the held-out coefficients, weighted so that the error is a plain sum of squares: [Re | Im] * sqrt(wq w)
        qidx = pack([(q, q=1,size(valset))], valset)
        if( size(qidx) > NQMAX ) qidx = qidx(1:size(qidx):(size(qidx) + NQMAX - 1) / NQMAX)
        nq = size(qidx)
        allocate(rm(m,2*nq))
        !$omp parallel do default(shared) private(j,q,sc) schedule(static)
        do q = 1, nq
            sc = sqrt(real(wq(qidx(q)),dp) * real(w(qidx(q)),dp))
            do j = 1, m
                rm(j,q)    = real(sc * real(y(qidx(q),j)), sp)
                rm(j,nq+q) = real(sc * aimag(y(qidx(q),j)), sp)
            end do
        end do
        !$omp end parallel do
        ! the pairs that may pool: same parity, not the validation frame itself
        allocate(allowed(nv,m))
        do j = 1, m
            do v = 1, nv
                allowed(v,j) = (even(j) .eqv. even(vinds(v))) .and. j /= vinds(v)
            end do
        end do
        allocate(d2(nv,m), source=0.d0)
        allocate(lw(nv,m), kmat(nv,m), pred(nv,2*nq), ksum(nv))
        rprev = 0
        do ir = 0, nr
            ! ir = 0: the temporal kernel alone; ir > 0: the latent squared distances over the leading RANKS(ir)
            if( ir > 0 )then
                r = RANKS(ir)
                !$omp parallel do default(shared) private(j,v,q) schedule(static)
                do j = 1, m
                    do v = 1, nv
                        do q = rprev + 1, r
                            d2(v,j) = d2(v,j) + (zeta(vinds(v),q) - zeta(j,q))**2
                        end do
                    end do
                end do
                !$omp end parallel do
                rprev = r
            endif
            do it = 1, nt
                tau = TAUS(it)
                if( ir == 0 .and. tau >= 0.d0 ) cycle
                if( ir >  0 .and. tau <  0.d0 ) cycle
                do ih = 1, nh
                    h = HS(ih)
                    !$omp parallel do default(shared) private(j,v) schedule(static)
                    do j = 1, m
                        do v = 1, nv
                            lw(v,j) = 0.d0
                            if( tau >= 0.d0 ) lw(v,j) = -0.5d0 * d2(v,j) / (2.d0 + tau * tau)
                            if( h > 0.d0 ) lw(v,j) = lw(v,j) - 0.5d0 * (real(times(vinds(v)) - times(j),dp) / h)**2
                        end do
                    end do
                    !$omp end parallel do
                    !$omp parallel do default(shared) private(v,j,lwmax) schedule(static)
                    do v = 1, nv
                        lwmax = -huge(1.d0)
                        do j = 1, m
                            if( allowed(v,j) ) lwmax = max(lwmax, lw(v,j))
                        end do
                        ksum(v) = 0.d0
                        do j = 1, m
                            kmat(v,j) = 0.
                            if( allowed(v,j) .and. lw(v,j) - lwmax > LWMIN ) kmat(v,j) = real(exp(lw(v,j) - lwmax), sp)
                            ksum(v) = ksum(v) + real(kmat(v,j),dp)
                        end do
                    end do
                    !$omp end parallel do
                    call sgemm('N', 'N', nv, 2*nq, m, 1.0, kmat, nv, rm, m, 0.0, pred, nv)
                    !$omp parallel do default(shared) private(v,i) schedule(static)
                    do v = 1, nv
                        i = max(ir, 1)
                        err(v,it,ih,i) = sum((real(rm(vinds(v),:),dp) - real(pred(v,:),dp) / (ksum(v) + WIENER))**2)
                    end do
                    !$omp end parallel do
                end do
            end do
        end do
    end subroutine cv_direction

    !> the targets' pre-images over the frames of their parity at (r, tau, h), with the weights the CV scored: the
    !! fit band's set-A coefficients weighted by the set-B model's latents and vice versa, every other coefficient by
    !! the mean of the two normalised weights; a target's own weight is the expected weight of a frame with its clean
    !! image (at distance 0 its own noise would dominate). ess: the Kish population of those mean weights
    subroutine preimages( y, seta, setb, times, even, targets, l_lat, lata, latb, r, tau, h, yden, ess )
        complex(sp),            intent(in)  :: y(:,:)
        logical,                intent(in)  :: seta(:), setb(:), even(:), l_lat
        integer,                intent(in)  :: times(:), targets(:), r
        type(whitened_latents), intent(in)  :: lata, latb
        real(dp),               intent(in)  :: tau, h
        complex(sp),            intent(out) :: yden(:,:)
        real(dp),               intent(out) :: ess(:)
        real(sp), allocatable :: ka(:,:), kb(:,:), ko(:,:)
        real(dp), allocatable :: wa(:), wb(:)
        integer,  allocatable :: ia(:), ib(:), io(:)
        integer :: m, nt, ncoeff, ti, q
        logical :: l_use
        m      = size(y,2)
        nt     = size(targets)
        ncoeff = size(y,1)
        l_use  = l_lat .and. r > 0 .and. tau >= 0.d0
        allocate(ka(nt,m), kb(nt,m), ko(nt,m), source=0.)
        !$omp parallel do default(shared) private(ti,wa,wb) schedule(dynamic,8)
        do ti = 1, nt
            call kernel_row(latb, targets(ti), wa)          ! for the set-A coefficients
            call kernel_row(lata, targets(ti), wb)          ! for the set-B coefficients
            ka(ti,:) = real(wa, sp)
            kb(ti,:) = real(wb, sp)
            ko(ti,:) = real(0.5d0 * (wa + wb), sp)
            ess(ti)  = kish_ess(0.5d0 * (wa + wb))
        end do
        !$omp end parallel do
        ia = pack([(q, q=1,ncoeff)], seta)
        ib = pack([(q, q=1,ncoeff)], setb)
        io = pack([(q, q=1,ncoeff)], .not. (seta .or. setb))
        call kernel_means(ka, ia)
        call kernel_means(kb, ib)
        call kernel_means(ko, io)

      contains

        !> one target's normalised weights over the frames, from the latents lat
        subroutine kernel_row( lat, i, wk )
            type(whitened_latents), intent(in)  :: lat
            integer,                intent(in)  :: i
            real(dp), allocatable,  intent(out) :: wk(:)
            real(dp) :: lw(m), lwmax, wsum
            integer  :: j
            allocate(wk(m), source=0.d0)
            lw    = 0.d0
            lwmax = -huge(1.d0)
            do j = 1, m
                if( even(j) .neqv. even(i) ) cycle
                if( l_use ) lw(j) = -0.5d0 * sum((lat%zeta(j,1:r) - lat%zeta(i,1:r))**2) / (2.d0 + tau * tau)
                if( h > 0.d0 ) lw(j) = lw(j) - 0.5d0 * (real(times(j) - times(i),dp) / h)**2
                ! the frame itself: its latent differs from a clean twin's by two frames' noise, variance 2 per direction
                if( j == i .and. l_use ) lw(j) = 0.5d0 * real(r,dp) * log((2.d0 + tau * tau) / (4.d0 + tau * tau))
                lwmax = max(lwmax, lw(j))
            end do
            do j = 1, m
                if( (even(j) .eqv. even(i)) .and. lw(j) - lwmax > LWMIN ) wk(j) = exp(lw(j) - lwmax)
            end do
            wsum = sum(wk)
            if( wsum <= 0.d0 )then
                wk(i) = 1.d0
                wsum  = 1.d0
            endif
            wk = wk / wsum
        end subroutine kernel_row

        !> yden(idx,:) = the weights kmat times the frames' coefficients idx, a chunk of coefficients at a time
        subroutine kernel_means( kmat, idx )
            real(sp), intent(in) :: kmat(:,:)
            integer,  intent(in) :: idx(:)
            real(sp), allocatable :: ych(:,:), out(:,:)
            integer :: c0, c1, nc, j, t
            do c0 = 1, size(idx), CHUNK
                c1 = min(size(idx), c0 + CHUNK - 1)
                nc = c1 - c0 + 1
                allocate(ych(m,2*nc), out(nt,2*nc))
                !$omp parallel do default(shared) private(j) schedule(static)
                do j = 1, m
                    ych(j,1:nc)      = real(y(idx(c0:c1),j))
                    ych(j,nc+1:2*nc) = aimag(y(idx(c0:c1),j))
                end do
                !$omp end parallel do
                call sgemm('N', 'N', nt, 2*nc, m, 1.0, kmat, nt, ych, m, 0.0, out, nt)
                do t = 1, nt
                    yden(idx(c0:c1),t) = cmplx(out(t,1:nc), out(t,nc+1:2*nc), kind=sp)
                end do
                deallocate(ych, out)
            end do
        end subroutine kernel_means

    end subroutine preimages

    !> the fit band in two interleaved coefficient sets, BLOCK-pixel Fourier blocks in a checkerboard (blocks,
    !! because the soft mask correlates neighbouring coefficients); as simple_flex_cls_kernel's
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

    !> the held-out coefficients cv_direction uses of a set (at most NQMAX, evenly strided)
    function pick( set ) result( mask )
        logical, intent(in) :: set(:)
        logical :: mask(size(set))
        integer, allocatable :: idx(:)
        integer :: q
        idx = pack([(q, q=1,size(set))], set)
        if( size(idx) > NQMAX ) idx = idx(1:size(idx):(size(idx) + NQMAX - 1) / NQMAX)
        mask = .false.
        mask(idx) = .true.
    end function pick

end module simple_flex_traj_denoise
