!@descr: trajectory_denoise pca_mode=flex: every frame the kernel pre-image of the frames that look like it in the latent of a flex class model and lie near it in time, rank and bandwidths by coefficient-split cross-validation; arrays in, arrays out
! Frames as given (no alignment, no interpolation). Without CTF and with one noise spectrum every frame's ML latent
! has the same covariance, so the kernel is Euclidean in noise-whitened latents and the CV and the pre-images are
! matrix products. The pre-images use the weights the CV scored: the fit band's two coefficient sets each averaged
! with weights from the other set's model (no selection bias), every other coefficient with the mean of the two.
! Two copies of the frames: soft-masked ones (ym) for the similarity (models, latents, CV) and the frames themselves
! (yu) for the averages, so the output is the whole box.
! Design, toy and validation: doc/implementation_notes/planned/trajectory_denoise_flex.md
module simple_flex_traj_denoise
use simple_defs,               only: sp, dp, DCMPLX_ZERO, logfhandle
use simple_error,              only: simple_exception
use simple_flex_cls_expansion, only: flex_cls_model, flex_cls_fit, flex_cls_weighted_mean
use simple_linalg,             only: cholesky, chol_forward, jacobi, spd_inverse
use simple_stat,               only: kish_ess, median
use simple_rnd,                only: seed_rnd_fixed
use simple_timer,              only: tic, toc, timer_int_kind
!$ use omp_lib,                only: omp_get_max_threads, omp_in_parallel
implicit none

public :: flex_traj_denoise, flex_traj_noise_spectrum, flex_traj_fit_band
private
#include "simple_local_flags.inc"

real(dp), parameter :: TAUS(10)  = [0.d0, 0.25d0, 0.5d0, 1.d0, 2.d0, 4.d0, 8.d0, 16.d0, 32.d0, -1.d0] !< noise units; < 0: none
integer,  parameter :: RANKS(8)  = [1, 2, 4, 8, 16, 32, 64, 128]   !< latent ranks tried, up to the models'
real(dp), parameter :: HS(8)     = [2.d0, 4.d0, 8.d0, 16.d0, 32.d0, 64.d0, 128.d0, -1.d0] !< frames; < 0: none
integer,  parameter :: BLOCKSZ   = 10       !< frames per block when the sets are windows
integer,  parameter :: NVMAX     = 1000     !< validation frames per set, at most
integer,  parameter :: NQMAX     = 1000     !< held-out coefficients per coefficient set, at most
integer,  parameter :: NMIN      = 12       !< fewest frames a set needs for the models
real(dp), parameter :: RIDGE     = 1.d-6    !< ridge of the ML latent system, relative to its largest diagonal
real(dp), parameter :: WIENER    = 0.1d0    !< CV prediction denominator + this
real(dp), parameter :: LWMIN     = -46.d0   !< log weights this far below the largest are zero
real(dp), parameter :: SNR_BAND  = 0.1d0    !< derived fit band: up to the last shell whose spectral SNR reaches this
integer,  parameter :: KBAND_MIN = 6        !< derived fit band: at least this many shells
integer,  parameter :: CHUNK     = 128      !< coefficients per matrix product of the pre-images (one per thread)
integer,  parameter :: NVB_MIN   = 16       !< validation frames per CV block, at least (one block per thread)
integer,  parameter :: NVB_MAX   = 64       !< and at most (bounds a thread's kernel block)
integer,  parameter :: LIBCOLS   = 1024     !< threaded BLAS: columns of a CV product (validation frames x bandwidths)
integer,  parameter :: LIBCHUNK  = 512      !< threaded BLAS: coefficients per pre-image product

! whether the linked BLAS threads one product (probe_blas, once per run). A serial BLAS is called inside the OpenMP
! team, one product per thread; a threaded one (an OpenMP OpenBLAS, MKL, a pthreads OpenBLAS) from the serial region
! with wide products, the kernels around it in parallel; inside a team (halfwin > 0) the compiler's matmul instead,
! since a pthreads BLAS would start its threads in every team thread. No environment variable either way
logical :: l_blas_probed   = .false.
logical :: l_blas_threaded = .false.
! the build's BLAS stays on the calling thread inside an OpenMP team (CMake: an OpenMP OpenBLAS, MKL, FlexiBLAS with
! an OpenMP backend), so a team may call it even when it threads serial-region products
#ifdef SIMPLE_BLAS_TEAM_SAFE
logical, parameter :: L_BLAS_TEAM_SAFE = .true.
#else
logical, parameter :: L_BLAS_TEAM_SAFE = .false.
#endif

!> a set's latents: ML latents of a model, noise-whitened (unit noise per coordinate), by decreasing scatter
type whitened_latents
    integer :: k = 0
    logical :: ok = .false.
    real(dp), allocatable :: zeta(:,:)   !< (m,k)
end type whitened_latents

contains

    !> ym, yu: (ncoeff,n) the half-plane coefficients of the soft-masked frames (the similarity) and of the frames
    !! themselves (the averages), frames as given. Returns every frame's pre-image yden and the Kish population of its
    !! weights, plus its CV choice (rank, tau, h; rank 0: the temporal kernel alone; tau, h < 0: none). halfwin <= 0:
    !! one set, the whole trajectory, every frame a candidate of every frame; halfwin > 0: blocks of BLOCKSZ frames,
    !! each with the frames within halfwin of its centre as its set. l_eo: a frame pools only frames of its own row
    !! parity (the even/odd halves stay independent). block: side of the Fourier blocks of the coefficient split
    !! (the tighter the mask, the more it correlates neighbouring coefficients, the larger the blocks)
    subroutine flex_traj_denoise( ym, yu, hidx, kidx, shell, wq, ncomp, kfit, halfwin, l_eo, seed, block, &
                                  &yden, ess, choice, kband )
        complex(sp), intent(in)  :: ym(:,:), yu(:,:)
        integer,     intent(in)  :: hidx(:), kidx(:), shell(:)
        real(sp),    intent(in)  :: wq(:)
        integer,     intent(in)  :: ncomp, kfit, halfwin, seed, block
        logical,     intent(in)  :: l_eo
        complex(sp), intent(out) :: yden(:,:)
        real(dp),    intent(out) :: ess(:)
        real(dp),    intent(out) :: choice(:,:)    !< (3,n) rank, tau, h per frame
        integer,     intent(out) :: kband
        real(sp), allocatable :: s2(:), w(:)
        logical,  allocatable :: fitmask(:), seta(:), setb(:)
        integer,  allocatable :: times(:)
        integer(timer_int_kind) :: t0
        integer :: n, ncoeff, nsh, b, nb, c0, c1, ctr, lo, hi, t
        n      = size(ym,2)
        ncoeff = size(ym,1)
        nsh    = maxval(shell)
        if( size(yu,1) /= ncoeff .or. size(yu,2) /= n ) THROW_HARD('flex_traj_denoise: ym and yu differ in shape')
        if( size(hidx) /= ncoeff .or. size(kidx) /= ncoeff .or. size(shell) /= ncoeff .or. size(wq) /= ncoeff ) &
            &THROW_HARD('flex_traj_denoise: lattice sizes')
        if( size(yden,1) /= ncoeff .or. size(yden,2) /= n .or. size(ess) /= n ) THROW_HARD('flex_traj_denoise: outputs')
        if( size(choice,1) /= 3 .or. size(choice,2) /= n ) THROW_HARD('flex_traj_denoise: choice shape')
        call probe_blas
        ! the noise spectrum and the fit band of the whole trajectory
        t0 = tic()
        call flex_traj_noise_spectrum(ym, shell, nsh, s2)
        allocate(w(ncoeff))
        w = 1. / s2(shell)
        kband = flex_traj_fit_band(ym, shell, nsh, s2, kfit)
        allocate(fitmask(ncoeff), seta(ncoeff), setb(ncoeff))
        fitmask = shell >= 1 .and. shell <= kband
        call coeff_sets(hidx, kidx, fitmask, max(block,1), seta, setb)
        write(logfhandle,'(A,I7,A,I4,A,I4,A,I7,A,I7,A,I7,A,I3)') 'Trajectory flex: frames=', n, ' fit band shells 1-', kband, &
            &' of ', nsh, ' coefficients=', count(fitmask), ' set A=', count(seta), ' set B=', count(setb), ' split block=', block
        write(logfhandle,'(A,*(1X,ES10.3))') 'Trajectory flex: noise spectrum (shells 1, 1/4, 1/2, 3/4, last):', &
            &s2(1), s2(max(1,nsh/4)), s2(max(1,nsh/2)), s2(max(1,(3*nsh)/4)), s2(nsh)
        write(logfhandle,'(A,F9.2)') 'Trajectory flex time: noise spectrum and fit band, s=', toc(t0)
        call flush(logfhandle)
        times = [(t, t=1,n)]
        if( halfwin <= 0 .or. 2 * halfwin + 1 >= n )then
            call denoise_set(ym, yu, w, wq, fitmask, seta, setb, times, ncomp, l_eo, seed, times, yden, ess, choice, 0, .true.)
        else
            ! the blocks are independent: one per thread (the work inside each runs on that thread)
            t0 = tic()
            nb = (n + BLOCKSZ - 1) / BLOCKSZ
            !$omp parallel do default(shared) private(b,c0,c1,ctr,lo,hi,t) schedule(dynamic,1)
            do b = 1, nb
                c0  = (b - 1) * BLOCKSZ + 1
                c1  = min(n, b * BLOCKSZ)
                ctr = (c0 + c1) / 2
                lo  = max(1, ctr - halfwin)
                hi  = min(n, ctr + halfwin)
                call denoise_set(ym(:,lo:hi), yu(:,lo:hi), w, wq, fitmask, seta, setb, times(lo:hi), ncomp, l_eo, seed + 10 * b, &
                    &[(t - lo + 1, t=c0,c1)], yden(:,c0:c1), ess(c0:c1), choice(:,c0:c1), b, .false.)
            end do
            !$omp end parallel do
            write(logfhandle,'(A,I6,A,F9.2)') 'Trajectory flex time: ', nb, ' window blocks, s=', toc(t0)
            call flush(logfhandle)
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
        real(dp) :: pw(nsh)
        integer  :: cnt(nsh), j, sh
        if( kfit > 0 )then
            kband = min(kfit, nsh)
            return
        endif
        pw  = 0.d0
        cnt = 0
        !$omp parallel do default(shared) private(j,sh) schedule(static) reduction(+:pw,cnt)
        do j = 1, size(shell)
            sh = shell(j)
            if( sh < 1 .or. sh > nsh ) cycle
            pw(sh)  = pw(sh) + sum(abs(cmplx(y(j,:), kind=dp))**2) / real(size(y,2),dp)
            cnt(sh) = cnt(sh) + 1
        end do
        !$omp end parallel do
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
    subroutine denoise_set( y, yu, w, wq, fitmask, seta, setb, times, ncomp, l_eo, seed, targets, yden, ess, choice, iset, &
                            &l_verbose )
        complex(sp), intent(in)  :: y(:,:), yu(:,:)   !< the masked frames (similarity), the frames (averages)
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
        integer(timer_int_kind) :: t0
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
                yden(:,i) = yu(:,targets(i))
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
        t0    = tic()
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
        if( l_verbose ) call log_time('the two set models and their latents', t0)
        call cross_validate(y, w, wq, fitmask, seta, setb, times, even, l_lat, lata, latb, kset, &
            &r, tau, h, cv_best, cv_mean, zscore, l_verbose)
        if( l_verbose ) call log_time('cross-validation', t0)
        call preimages(yu, seta, setb, times, even, targets, l_lat, lata, latb, r, tau, h, yden, ess)
        if( l_verbose ) call log_time('pre-images', t0)
        do i = 1, size(targets)
            choice(:,i) = [real(r,dp), tau, h]
        end do
        if( l_verbose .or. iset > 0 )then
            !$omp critical (traj_flex_log)
            write(logfhandle,'(A,I6,A,I6,A,I6,A,I4,A,I4,A,F6.2,A,F7.1,A,F9.5,A,F9.5,A,F8.2,A,F8.1)') &
                &'Trajectory flex set=', iset, ' frames ', times(1), '-', times(m), ' model rank=', kset, &
                &' chosen rank=', r, ' tau=', tau, ' h=', h, ' CV=', cv_best, ' (parity mean ', cv_mean, ') z=', zscore, &
                &' mean neff=', sum(ess) / real(size(ess),dp)
            call flush(logfhandle)
            !$omp end critical (traj_flex_log)
        endif
    end subroutine denoise_set

    !> one timing line, then the clock restarts
    subroutine log_time( step, t0 )
        character(len=*),        intent(in)    :: step
        integer(timer_int_kind), intent(inout) :: t0
        write(logfhandle,'(A,A,A,F9.2)') 'Trajectory flex time: ', step, ', s=', toc(t0)
        call flush(logfhandle)
        t0 = tic()
    end subroutine log_time

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
    !! coefficients in valset (at most NQMAX of them) from the other frames; err(v,it,ih,ir) the weighted squared error.
    !! The validation frames go in blocks, each through the whole grid, every temporal bandwidth of a (rank, tau) in one
    !! product: with a serial BLAS one or more blocks per thread, each its own products; with a threaded BLAS wide
    !! blocks one after the other, the products from the serial region (probe_blas)
    subroutine cv_direction( zeta, y, w, wq, valset, times, even, vinds, nr, err )
        real(dp),    intent(in)    :: zeta(:,:)
        complex(sp), intent(in)    :: y(:,:)
        real(sp),    intent(in)    :: w(:), wq(:)
        logical,     intent(in)    :: valset(:), even(:)
        integer,     intent(in)    :: times(:), vinds(:), nr
        real(dp),    intent(inout) :: err(:,:,:,:)
        real(sp), allocatable :: rmt(:,:)
        integer,  allocatable :: qidx(:)
        real(dp) :: sc
        integer  :: m, nv, nq, nt, nh, j, q, nthr, nvb, nblk, ib, v0
        m  = size(y,2)
        nv = size(vinds)
        nt = size(TAUS)
        nh = size(HS)
        qidx = pack([(q, q=1,size(valset))], valset)
        if( size(qidx) > NQMAX ) qidx = qidx(1:size(qidx):(size(qidx) + NQMAX - 1) / NQMAX)
        nq = size(qidx)
        ! the held-out coefficients of every frame, [Re; Im] * sqrt(wq w): the error is a plain sum of squares
        allocate(rmt(2*nq,m))
        !$omp parallel do default(shared) private(j,q,sc) schedule(static)
        do j = 1, m
            do q = 1, nq
                sc = sqrt(real(wq(qidx(q)),dp) * real(w(qidx(q)),dp))
                rmt(q,j)    = real(sc * real(y(qidx(q),j)), sp)
                rmt(nq+q,j) = real(sc * aimag(y(qidx(q),j)), sp)
            end do
        end do
        !$omp end parallel do
        nthr = 1
        !$ if( .not. omp_in_parallel() ) nthr = omp_get_max_threads()
        if( l_blas_threaded .or. nthr == 1 )then
            nvb = max(1, min(nv, LIBCOLS / nh))
            do v0 = 0, nv - 1, nvb
                call cv_block(v0, min(nvb, nv - v0), nthr > 1)
            end do
        else
            nvb  = min(NVB_MAX, max(NVB_MIN, (nv + nthr - 1) / nthr))
            nblk = (nv + nvb - 1) / nvb
            !$omp parallel do default(shared) private(ib) schedule(dynamic,1) proc_bind(close)
            do ib = 1, nblk
                call cv_block((ib - 1) * nvb, min(nvb, nv - (ib - 1) * nvb), .false.)
            end do
            !$omp end parallel do
        endif

      contains

        !> validation frames v0+1..v0+nb through the whole grid; l_par: their kernels and errors in parallel
        subroutine cv_block( v0, nb, l_par )
            integer, intent(in) :: v0, nb
            logical, intent(in) :: l_par
            real(sp), allocatable :: kloc(:,:), ploc(:,:)
            real(dp), allocatable :: d2(:,:), ksum(:)
            real(dp) :: tau
            integer  :: ir, ie, it, v, i, q, rprev, col
            allocate(d2(m,nb), source=0.d0)
            allocate(kloc(m,nb*nh), ploc(2*nq,nb*nh), ksum(nb*nh))
            rprev = 0
            do ir = 0, nr
                ! ir = 0: the temporal kernel alone; ir > 0: the latent squared distances over the leading RANKS(ir)
                if( ir > 0 )then
                    !$omp parallel do default(shared) private(v,i,q) schedule(static) if(l_par)
                    do v = 1, nb
                        i = vinds(v0 + v)
                        do q = rprev + 1, RANKS(ir)
                            d2(:,v) = d2(:,v) + (zeta(:,q) - zeta(i,q))**2
                        end do
                    end do
                    !$omp end parallel do
                    rprev = RANKS(ir)
                endif
                ie = max(ir, 1)
                do it = 1, nt
                    tau = TAUS(it)
                    if( ir == 0 .and. tau >= 0.d0 ) cycle
                    if( ir >  0 .and. tau <  0.d0 ) cycle
                    ! column (ih - 1) * nb + v: validation frame v0 + v at temporal bandwidth HS(ih)
                    !$omp parallel do default(shared) private(col,v) schedule(static) if(l_par)
                    do col = 1, nb * nh
                        v = mod(col - 1, nb) + 1
                        call kernel_col(v0 + v, (col - 1) / nb + 1, tau, d2(:,v), kloc(:,col), ksum(col))
                    end do
                    !$omp end parallel do
                    call gemm_nn(rmt, kloc, ploc)
                    !$omp parallel do default(shared) private(col,v) schedule(static) if(l_par)
                    do col = 1, nb * nh
                        v = v0 + mod(col - 1, nb) + 1
                        err(v,it,(col - 1) / nb + 1,ie) = sum((real(rmt(:,vinds(v)),dp) - &
                            &real(ploc(:,col),dp) / (ksum(col) + WIENER))**2)
                    end do
                    !$omp end parallel do
                end do
            end do
            deallocate(d2, kloc, ploc, ksum)
        end subroutine cv_block

        !> validation frame v's kernel at temporal bandwidth HS(ih) over the frames of its parity, itself left out
        subroutine kernel_col( v, ih, tau, d2v, kcol, ks )
            integer,  intent(in)  :: v, ih
            real(dp), intent(in)  :: tau, d2v(:)
            real(sp), intent(out) :: kcol(:)
            real(dp), intent(out) :: ks
            real(dp), allocatable :: lw(:)
            real(dp) :: lwmax
            integer  :: i, j
            logical  :: l_ok
            allocate(lw(m))
            i     = vinds(v)
            lwmax = -huge(1.d0)
            do j = 1, m
                lw(j) = 0.d0
                if( tau >= 0.d0 ) lw(j) = -0.5d0 * d2v(j) / (2.d0 + tau * tau)
                if( HS(ih) > 0.d0 ) lw(j) = lw(j) - 0.5d0 * (real(times(i) - times(j),dp) / HS(ih))**2
                if( (even(j) .eqv. even(i)) .and. j /= i ) lwmax = max(lwmax, lw(j))
            end do
            ks = 0.d0
            do j = 1, m
                kcol(j) = 0.
                l_ok = (even(j) .eqv. even(i)) .and. j /= i
                if( l_ok )then
                    if( lw(j) - lwmax > LWMIN ) kcol(j) = real(exp(lw(j) - lwmax), sp)
                endif
                ks = ks + real(kcol(j),dp)
            end do
            deallocate(lw)
        end subroutine kernel_col

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
        ! the weights, a column per target
        allocate(ka(m,nt), kb(m,nt), ko(m,nt), source=0.)
        !$omp parallel do default(shared) private(ti,wa,wb) schedule(dynamic,8)
        do ti = 1, nt
            call kernel_row(latb, targets(ti), wa)          ! for the set-A coefficients
            call kernel_row(lata, targets(ti), wb)          ! for the set-B coefficients
            ka(:,ti) = real(wa, sp)
            kb(:,ti) = real(wb, sp)
            ko(:,ti) = real(0.5d0 * (wa + wb), sp)
            ess(ti)  = kish_ess(0.5d0 * (wa + wb))
        end do
        !$omp end parallel do
        ia = pack([(q, q=1,ncoeff)], seta)
        ib = pack([(q, q=1,ncoeff)], setb)
        io = pack([(q, q=1,ncoeff)], .not. (seta .or. setb))
        call kernel_means(y, ka, ia, yden)
        call kernel_means(y, kb, ib, yden)
        call kernel_means(y, ko, io, yden)

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

    end subroutine preimages

    !> yden(idx,:) = the frames' coefficients idx (rows of y) times the weights kmat (frames x targets), in chunks of
    !! coefficients: with a serial BLAS one or more chunks per thread, each its own product; with a threaded BLAS
    !! wider chunks one after the other, the products from the serial region (probe_blas)
    subroutine kernel_means( y, kmat, idx, yden )
        complex(sp),          intent(in)    :: y(:,:)
        real(sp), contiguous, intent(in)    :: kmat(:,:)
        integer,              intent(in)    :: idx(:)
        complex(sp),          intent(inout) :: yden(:,:)
        integer :: m, nt, nthr, nchunk, csz, ic, c0
        if( size(idx) == 0 ) return
        m    = size(y,2)
        nt   = size(kmat,2)
        nthr = 1
        !$ if( .not. omp_in_parallel() ) nthr = omp_get_max_threads()
        if( l_blas_threaded .or. nthr == 1 )then
            csz = min(size(idx), LIBCHUNK)
            do c0 = 1, size(idx), csz
                call mean_chunk(c0, min(size(idx), c0 + csz - 1), nthr > 1)
            end do
        else
            ! chunks of at most CHUNK coefficients, their number a multiple of the threads (even shares)
            nchunk = nthr * max(1, (size(idx) + nthr * CHUNK - 1) / (nthr * CHUNK))
            csz    = (size(idx) + nchunk - 1) / nchunk
            nchunk = (size(idx) + csz - 1) / csz
            !$omp parallel do default(shared) private(ic) schedule(dynamic,1) proc_bind(close)
            do ic = 1, nchunk
                call mean_chunk((ic - 1) * csz + 1, min(size(idx), ic * csz), .false.)
            end do
            !$omp end parallel do
        endif

      contains

        !> the coefficients idx(c0:c1); l_par: their gather and scatter in parallel
        subroutine mean_chunk( c0, c1, l_par )
            integer, intent(in) :: c0, c1
            logical, intent(in) :: l_par
            real(sp), allocatable :: ycht(:,:), outt(:,:)
            integer :: nc, j, t
            nc = c1 - c0 + 1
            allocate(ycht(2*nc,m), outt(2*nc,nt))
            !$omp parallel do default(shared) private(j) schedule(static) if(l_par)
            do j = 1, m
                ycht(1:nc,j)      = real(y(idx(c0:c1),j))
                ycht(nc+1:2*nc,j) = aimag(y(idx(c0:c1),j))
            end do
            !$omp end parallel do
            call gemm_nn(ycht, kmat, outt)
            !$omp parallel do default(shared) private(t) schedule(static) if(l_par)
            do t = 1, nt
                yden(idx(c0:c1),t) = cmplx(outt(1:nc,t), outt(nc+1:2*nc,t), kind=sp)
            end do
            !$omp end parallel do
            deallocate(ycht, outt)
        end subroutine mean_chunk

    end subroutine kernel_means

    !> once per run: does the linked BLAS run one product on several threads? The process CPU time of serial-region
    !! products against their wall time, after a warm-up long enough for the OpenMP threads of earlier regions to stop
    !! spinning (their spinning would count as BLAS threads)
    subroutine probe_blas
        integer,  parameter :: NP = 512, NREP = 8
        real,     parameter :: WARMUP = 0.1     !< s
        real(sp), allocatable :: a(:,:), b(:,:), c(:,:)
        integer(timer_int_kind) :: t0
        real    :: c0, c1, wall, ratio
        integer :: irep
        if( l_blas_probed ) return
        allocate(a(NP,NP), b(NP,NP), c(NP,NP))
        a  = 1.e-3
        b  = 1.e-3
        t0 = tic()
        do
            call sgemm('N', 'N', NP, NP, NP, 1.0, a, NP, b, NP, 0.0, c, NP)
            if( real(toc(t0)) >= WARMUP ) exit
        end do
        call cpu_time(c0)
        t0 = tic()
        do irep = 1, NREP
            call sgemm('N', 'N', NP, NP, NP, 1.0, a, NP, b, NP, 0.0, c, NP)
        end do
        wall = max(real(toc(t0)), 1.e-6)
        call cpu_time(c1)
        ratio           = (c1 - c0) / wall
        l_blas_threaded = ratio > 1.8
        l_blas_probed   = .true.
        if( l_blas_threaded )then
            write(logfhandle,'(A,F6.2,A)') 'Trajectory flex: BLAS cpu/wall ', ratio, &
                &' (threaded): wide products from the serial region'
        else
            write(logfhandle,'(A,F6.2,A)') 'Trajectory flex: BLAS cpu/wall ', ratio, &
                &' (serial): one BLAS product per thread'
        endif
        call flush(logfhandle)
    end subroutine probe_blas

    !> c = a b, single precision: the BLAS, except inside an OpenMP team when the BLAS threads (probe_blas) and the
    !! build does not know it stays on the calling thread there (L_BLAS_TEAM_SAFE)
    subroutine gemm_nn( a, b, c )
        real(sp), contiguous, intent(in)  :: a(:,:), b(:,:)
        real(sp), contiguous, intent(out) :: c(:,:)
        logical :: l_team
        l_team = .false.
        !$ l_team = omp_in_parallel()
        if( l_team .and. l_blas_threaded .and. .not. L_BLAS_TEAM_SAFE )then
            c = matmul(a, b)
        else
            call sgemm('N', 'N', size(a,1), size(b,2), size(a,2), 1.0, a, size(a,1), b, size(b,1), 0.0, c, size(c,1))
        endif
    end subroutine gemm_nn

    !> the fit band in two interleaved coefficient sets, block-pixel Fourier blocks in a checkerboard (blocks,
    !! because the soft mask correlates neighbouring coefficients); as simple_flex_cls_kernel's
    pure subroutine coeff_sets( hidx, kidx, fitmask, block, seta, setb )
        integer, intent(in)  :: hidx(:), kidx(:), block
        logical, intent(in)  :: fitmask(:)
        logical, intent(out) :: seta(:), setb(:)
        integer :: j
        logical :: l_a
        do j = 1, size(fitmask)
            l_a     = modulo(floor(real(hidx(j)) / block) + floor(real(kidx(j)) / block), 2) == 0
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
