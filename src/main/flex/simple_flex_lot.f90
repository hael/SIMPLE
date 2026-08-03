!@descr: Linearised optimal-transport (LOT) coordinates for a set of conformational volumes.
!
!        MOTIVATION. flex_pca describes conformational variability in DENSITY space,
!        where a domain motion of amplitude R*Phi costs ~2*R*Phi/d principal components at
!        resolution d -- so a large-amplitude motion is spread over many components and no
!        single one carries it. Optimal transport measures the same variability by how far
!        mass MOVES, which costs O(1) coordinates regardless of amplitude. Measured on 20
!        IgG-RL state maps (box 64, lp 15): density-space PCA needs 13 components for 90% of
!        the variance and its best single component correlates 0.556 with the true Fab swing;
!        the LOT embedding needs 3 and reaches 0.927. Ordering the states along the transport
!        mode gives corr(order, true swing) = 0.936 against 0.685 for the delivered
!        density-space ordering (oracle 0.988).
!
!        WHAT THIS MODULE DOES AND DOES NOT DO. It produces COORDINATES, never density. No
!        volume is ever warped: the transport map is used to describe and order the state maps
!        that weighted backprojection already produced, so every delivered voxel remains a
!        weighted sum of measured data and FSC/local-resolution/sharpening all keep their usual
!        meaning. Warping a map to produce output would break every one of those (a warp is a
!        spatially varying resampling, so it low-passes by a position-dependent amount and makes
!        local resolution a function of the Jacobian rather than of the signal).
!
!        SCOPE LIMIT, MEASURED. LOT coordinates cannot be transferred to individual particles by
!        fitting psi = A*z on the states: cross-validated, no linear function of z beats the best
!        single latent component (rho 0.7903 vs 0.7971), because for any smooth reparameterisation
!        the between/within variance ratio is invariant to first order. Per-particle transport
!        coordinates would have to be fitted from the images directly. This module therefore
!        operates on VOLUMES only.
!
!        WHICH MODE IS THE MOTION IS NOT SOLVED. Five ground-truth-free selectors were tested and
!        all failed to pick the swing mode: half-set reproducibility (~1.000 for every mode, since
!        both half-set state sets sit at the same latent targets), Helmholtz curl fraction,
!        displacement localisation, rigid-domain fit, and the induced-density relevance score that
!        works for density eigenvolumes. `lot_mode` therefore selects the mode by hand, which is
!        tractable because only ~3 modes carry variance; the induced-density maps written here are
!        meant to be looked at.
!
!        See doc/for_developers/ideas/flex_pca_state_and_directions.md Part 9.
module simple_flex_lot
use simple_core_module_api
use simple_image,         only: image
use simple_parameters,    only: parameters
use simple_linalg,        only: jacobi, eigsrt
use simple_srch_sort_loc, only: hpsort
implicit none
private
#include "simple_local_flags.inc"

public :: lot_extract_cloud, lot_sinkhorn_map, lot_embed_volumes, lot_modes_from_disp
public :: lot_pullback_metric

! Sinkhorn defaults. EPS is in (voxel)^2 and is the entropic regularisation; it also acts as a
! transport length scale, since exp(-C/eps) underflows past |d| ~ sqrt(700*eps) voxels. At the
! default 6.0 that is ~65 voxels, comfortably beyond any real domain motion at box 64.
real(dp), parameter :: LOT_EPS_DEFAULT  = 6.0_dp
integer,  parameter :: LOT_NITER_DEFAULT = 120
integer,  parameter :: LOT_NPT_MAX       = 6000   ! npt^2 doubles must stay in memory

contains

    !>  Masked, density-weighted point cloud of a volume: the npt heaviest non-negative voxels
    !!  inside mskrad. Coordinates are returned in VOXELS relative to the box centre and the
    !!  weights are normalised to sum to one, so clouds of different volumes are directly
    !!  comparable as probability measures.
    !!
    !!  Taking the heaviest voxels rather than thresholding keeps the cloud size fixed across
    !!  volumes, which is what lets every state's displacement field live in one common vector
    !!  space indexed by the reference points.
    subroutine lot_extract_cloud( rmat, ldim, mskrad, npt, pts, wts )
        real,     intent(in)  :: rmat(:,:,:)
        integer,  intent(in)  :: ldim(3), npt
        real,     intent(in)  :: mskrad
        real(dp), allocatable, intent(out) :: pts(:,:)     ! (npt,3) voxel coords, box-centred
        real(dp), allocatable, intent(out) :: wts(:)       ! (npt) normalised
        real,    allocatable :: vals(:)
        integer, allocatable :: idx(:)
        integer  :: cx, cy, cz, i, j, k, n, m
        real     :: r2, mr2, v
        real(dp) :: wsum
        if( npt < 8 )          THROW_HARD('lot_extract_cloud: npt too small')
        if( npt > LOT_NPT_MAX ) THROW_HARD('lot_extract_cloud: npt exceeds LOT_NPT_MAX; the Sinkhorn &
            &kernel is npt^2 doubles')
        cx = ldim(1)/2 + 1; cy = ldim(2)/2 + 1; cz = ldim(3)/2 + 1
        mr2 = mskrad*mskrad
        ! count admissible voxels first so the gather is exact
        n = 0
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    r2 = real((i-cx)**2 + (j-cy)**2 + (k-cz)**2)
                    if( r2 > mr2 ) cycle
                    if( rmat(i,j,k) <= 0. ) cycle
                    n = n + 1
                end do
            end do
        end do
        if( n < npt ) THROW_HARD('lot_extract_cloud: fewer positive voxels inside the mask than npt')
        allocate(vals(n), idx(n))
        n = 0
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    r2 = real((i-cx)**2 + (j-cy)**2 + (k-cz)**2)
                    if( r2 > mr2 ) cycle
                    v = rmat(i,j,k)
                    if( v <= 0. ) cycle
                    n = n + 1
                    vals(n) = -v                       ! negate: hpsort is ascending, we want the top
                    idx(n)  = (k-1)*ldim(2)*ldim(1) + (j-1)*ldim(1) + i
                end do
            end do
        end do
        call hpsort(vals, idx)
        allocate(pts(npt,3), wts(npt))
        wsum = 0._dp
        do m = 1, npt
            n = idx(m) - 1
            i = mod(n, ldim(1)) + 1
            j = mod(n/ldim(1), ldim(2)) + 1
            k = n/(ldim(1)*ldim(2)) + 1
            pts(m,1) = real(i-cx, dp)
            pts(m,2) = real(j-cy, dp)
            pts(m,3) = real(k-cz, dp)
            wts(m)   = real(-vals(m), dp)
            wsum     = wsum + wts(m)
        end do
        if( wsum <= 0._dp ) THROW_HARD('lot_extract_cloud: zero total mass')
        wts = wts / wsum
        deallocate(vals, idx)
    end subroutine lot_extract_cloud

    !>  Entropic optimal transport between two weighted point clouds by Sinkhorn scaling, and the
    !!  barycentric projection of the resulting plan, T(x_i) = sum_j P_ij y_j / sum_j P_ij.
    !!
    !!  Plain (non-log-domain) Sinkhorn on purpose: the whole iteration is then two DGEMVs against
    !!  a precomputed kernel, which is orders of magnitude faster than a log-domain version and is
    !!  numerically fine at the eps used here. Rows of the kernel that underflow entirely are
    !!  guarded, and the achieved marginal error is returned so a caller can check convergence
    !!  rather than assume it.
    subroutine lot_sinkhorn_map( xs, a, ns, xt, b, nt, eps, niter, tmap, marg_err )
        integer,  intent(in)  :: ns, nt, niter
        real(dp), intent(in)  :: xs(ns,3), a(ns), xt(nt,3), b(nt), eps
        real(dp), intent(out) :: tmap(ns,3)
        real(dp), intent(out) :: marg_err
        real(dp), allocatable :: kmat(:,:), u(:), v(:), kv(:), ktu(:), rs(:)
        real(dp) :: c, d1, d2, d3
        integer  :: i, j, it
        ! Sinkhorn scalings legitimately reach ~1e-300 in the tails of the kernel; guarding with
        ! DTINY (1e-10) instead clamps them and stalls the iteration -- measured, it left the
        ! marginal error at 1e-2 where the true floor gives 1e-7.
        real(dp), parameter :: SK_FLOOR = 1.0e-300_dp
        allocate(kmat(ns,nt), u(ns), v(nt), kv(ns), ktu(nt), rs(ns))
        !$omp parallel do default(shared) private(i,j,d1,d2,d3,c) schedule(static) collapse(2)
        do j = 1, nt
            do i = 1, ns
                d1 = xs(i,1) - xt(j,1); d2 = xs(i,2) - xt(j,2); d3 = xs(i,3) - xt(j,3)
                c  = d1*d1 + d2*d2 + d3*d3
                kmat(i,j) = exp(-c/eps)
            end do
        end do
        !$omp end parallel do
        u = 1._dp; v = 1._dp
        do it = 1, niter
            call dgemv('N', ns, nt, 1._dp, kmat, ns, v, 1, 0._dp, kv, 1)
            do i = 1, ns
                u(i) = a(i) / max(kv(i), SK_FLOOR)
            end do
            call dgemv('T', ns, nt, 1._dp, kmat, ns, u, 1, 0._dp, ktu, 1)
            do j = 1, nt
                v(j) = b(j) / max(ktu(j), SK_FLOOR)
            end do
        end do
        ! row sums of the plan, and the barycentric projection
        call dgemv('N', ns, nt, 1._dp, kmat, ns, v, 1, 0._dp, kv, 1)
        do i = 1, ns
            rs(i) = u(i)*kv(i)
        end do
        marg_err = sum(abs(rs - a))
        tmap = 0._dp
        !$omp parallel do default(shared) private(i,j,c) schedule(static)
        do i = 1, ns
            do j = 1, nt
                c = u(i)*kmat(i,j)*v(j)
                tmap(i,1) = tmap(i,1) + c*xt(j,1)
                tmap(i,2) = tmap(i,2) + c*xt(j,2)
                tmap(i,3) = tmap(i,3) + c*xt(j,3)
            end do
            tmap(i,:) = tmap(i,:) / max(rs(i), SK_FLOOR)
        end do
        !$omp end parallel do
        deallocate(kmat, u, v, kv, ktu, rs)
    end subroutine lot_sinkhorn_map

    !>  LOT embedding of nvol volumes: reference = their arithmetic mean, and each volume's
    !!  coordinate is the displacement field u_s = T_s - id sampled at the reference points.
    !!  All displacement fields are indexed by the same reference cloud, so they live in one
    !!  vector space and Euclidean distance between them approximates the 2-Wasserstein distance.
    subroutine lot_embed_volumes( vols, nvol, ldim, mskrad, npt, eps, niter, disp, refpts, refwts )
        type(image), intent(inout) :: vols(nvol)
        integer,     intent(in)    :: nvol, ldim(3), npt, niter
        real,        intent(in)    :: mskrad
        real(dp),    intent(in)    :: eps
        real(dp), allocatable, intent(out) :: disp(:,:)     ! (nvol, 3*npt)
        real(dp), allocatable, intent(out) :: refpts(:,:)   ! (npt,3)
        real(dp), allocatable, intent(out) :: refwts(:)
        real,     allocatable :: rmean(:,:,:), rtmp(:,:,:)
        real(dp), allocatable :: xt(:,:), bt(:), tmap(:,:)
        real(dp) :: merr, umean, umax, un
        integer  :: s, i
        if( nvol < 3 ) THROW_HARD('lot_embed_volumes: need at least 3 volumes')
        allocate(rmean(ldim(1),ldim(2),ldim(3)), source=0.)
        do s = 1, nvol
            rtmp  = vols(s)%get_rmat()
            rmean = rmean + rtmp(:ldim(1),:ldim(2),:ldim(3))
            deallocate(rtmp)
        end do
        rmean = rmean / real(nvol)
        call lot_extract_cloud(rmean, ldim, mskrad, npt, refpts, refwts)
        deallocate(rmean)
        allocate(disp(nvol, 3*npt), tmap(npt,3))
        write(logfhandle,'(A,I0,A,I0,A,F6.2)') '>>> FLEX_LOT transport: ', nvol, ' volumes, ', &
            &npt, ' points, eps=', eps
        call flush(logfhandle)
        do s = 1, nvol
            rtmp = vols(s)%get_rmat()
            call lot_extract_cloud(rtmp(:ldim(1),:ldim(2),:ldim(3)), ldim, mskrad, npt, xt, bt)
            deallocate(rtmp)
            call lot_sinkhorn_map(refpts, refwts, npt, xt, bt, npt, eps, niter, tmap, merr)
            umean = 0._dp; umax = 0._dp
            do i = 1, npt
                disp(s, 3*i-2) = tmap(i,1) - refpts(i,1)
                disp(s, 3*i-1) = tmap(i,2) - refpts(i,2)
                disp(s, 3*i  ) = tmap(i,3) - refpts(i,3)
                un = sqrt(sum((tmap(i,:) - refpts(i,:))**2))
                umean = umean + un
                umax  = max(umax, un)
            end do
            umean = umean / real(npt, dp)
            write(logfhandle,'(A,I3,A,F7.2,A,F7.2,A,ES9.2)') '>>>   vol ', s, '  mean |u| ', &
                &umean, '  max ', umax, ' voxels  marginal err ', merr
            deallocate(xt, bt)
        end do
        call flush(logfhandle)
        deallocate(tmap)
    end subroutine lot_embed_volumes

    !>  PCA of the displacement fields. nd = 3*npt is far larger than nvol, so the modes are
    !!  obtained from the nvol x nvol Gram matrix rather than by a full SVD -- the same idiom
    !!  orthonormalize_representatives uses in the covariance path.
    !!
    !!  MEAN-CENTRING IS LOAD-BEARING, not cosmetic. The mean displacement field is dominated by
    !!  an inward contraction (measured: 71 % of points move inward, |radial|/|total| = 0.613),
    !!  because the arithmetic mean volume has larger support than any individual state and
    !!  transporting it to any one of them must contract. That contraction is NOT conformational
    !!  and it is exactly what sinks a Wasserstein-barycenter reference. Removing the mean takes
    !!  the net radial component to 0.000 and the radial share to 0.438, leaving the variation
    !!  that does carry the motion.
    subroutine lot_modes_from_disp( disp, nvol, nd, scores, modes, varfrac )
        integer,  intent(in)  :: nvol, nd
        real(dp), intent(in)  :: disp(nvol, nd)
        real(dp), allocatable, intent(out) :: scores(:,:)   ! (nvol, nvol) volume coords per mode
        real(dp), allocatable, intent(out) :: modes(:,:)    ! (nvol, nd)   unit displacement fields
        real(dp), allocatable, intent(out) :: varfrac(:)    ! (nvol)
        real(dp), allocatable :: dc(:,:), gram(:,:), evec(:,:), eval(:), mu(:)
        real(dp) :: tot, nrm
        integer  :: s, t, k, nrot
        allocate(dc(nvol,nd), mu(nd))
        do k = 1, nd
            mu(k) = sum(disp(:,k)) / real(nvol, dp)
        end do
        do s = 1, nvol
            dc(s,:) = disp(s,:) - mu
        end do
        allocate(gram(nvol,nvol), evec(nvol,nvol), eval(nvol))
        call dgemm('N', 'T', nvol, nvol, nd, 1._dp, dc, nvol, dc, nvol, 0._dp, gram, nvol)
        nrot = 0
        call jacobi(gram, nvol, nvol, eval, evec, nrot)
        call eigsrt(eval, evec, nvol, nvol)                 ! descending
        allocate(scores(nvol,nvol), modes(nvol,nd), varfrac(nvol))
        tot = max(sum(max(eval, 0._dp)), DTINY)
        do k = 1, nvol
            varfrac(k) = max(eval(k), 0._dp) / tot
            ! score of volume s on mode k is sqrt(lambda_k) * evec(s,k)
            nrm = sqrt(max(eval(k), 0._dp))
            do s = 1, nvol
                scores(s,k) = nrm * evec(s,k)
            end do
            ! mode k as a displacement field: sum_s evec(s,k) * dc(s,:) / sqrt(lambda_k)
            modes(k,:) = 0._dp
            do s = 1, nvol
                modes(k,:) = modes(k,:) + evec(s,k)*dc(s,:)
            end do
            if( nrm > DTINY )then
                modes(k,:) = modes(k,:) / nrm
            else
                modes(k,:) = 0._dp
            endif
        end do
        deallocate(dc, mu, gram, evec, eval)
    end subroutine lot_modes_from_disp

    !>  Pullback of the transport geometry onto the latent space, as a metric for the state kernel.
    !!
    !!  The transport coordinates live on VOLUMES, so they cannot weight particles directly. But the
    !!  states have known latent targets, so a ridge-regularised fit `psi ~ J' z` over the (target,
    !!  transport score) pairs gives a linear map, and `M = J J'` is its pullback: a metric on the
    !!  latent space that measures distance by how much DENSITY MOVES rather than by how far the
    !!  coordinates differ. Two particles separated only along a latent direction that shifts no
    !!  mass describe the same conformation and should pool; the plain metric separates them.
    !!
    !!  This is NOT the refuted `psi = A z` reparameterisation. That replaced the coordinate, was
    !!  fitted with 21 parameters on 20 states, and moved the k-means TARGETS. This keeps the
    !!  coordinate and the targets and changes only the metric, at rank nmodes, blended with the
    !!  identity so the fit cannot dominate. `blend = 0` recovers the current behaviour exactly.
    !!
    !!  Measured offline on run 49 with a simplified kernel: mean weighted ground-truth spread per
    !!  state 4.049 -> 3.619 (blend 1.0) -> 3.367 (blend 0.5) at comparable effective support, i.e.
    !!  each state pools a more conformationally homogeneous particle set. Blend 0.5 is the default
    !!  because the pure pullback is rank-nmodes and would collapse the remaining directions.
    subroutine lot_pullback_metric( scores, nvol, nmodes, targets, zsd, ncomp, ridge, blend, metric )
        integer,  intent(in)  :: nvol, nmodes, ncomp
        real(dp), intent(in)  :: scores(nvol,nvol)      ! transport scores per state
        real(dp), intent(in)  :: targets(ncomp,nvol)    ! latent target of each state
        real(dp), intent(in)  :: zsd(ncomp)             ! per-component latent sd, for standardisation
        real(dp), intent(in)  :: ridge, blend
        real(dp), allocatable, intent(out) :: metric(:,:)
        real(dp), allocatable :: ts(:,:), gram(:,:), rhs(:,:), jmat(:,:), mm(:,:)
        real(dp) :: tr, sc
        integer  :: q, r, s, k, nm, info
        integer,  allocatable :: ipiv(:)
        nm = max(1, min(nmodes, nvol))
        allocate(ts(nvol,ncomp))
        do s = 1, nvol
            do q = 1, ncomp
                ts(s,q) = targets(q,s) / max(zsd(q), DTINY)
            end do
        end do
        allocate(gram(ncomp,ncomp), rhs(ncomp,nm), jmat(ncomp,nm), ipiv(ncomp))
        call dgemm('T', 'N', ncomp, ncomp, nvol, 1._dp, ts, nvol, ts, nvol, 0._dp, gram, ncomp)
        do q = 1, ncomp
            gram(q,q) = gram(q,q) + ridge
        end do
        do k = 1, nm
            do q = 1, ncomp
                rhs(q,k) = 0._dp
                do s = 1, nvol
                    rhs(q,k) = rhs(q,k) + ts(s,q)*scores(s,k)
                end do
            end do
        end do
        jmat = rhs
        call dposv('U', ncomp, nm, gram, ncomp, jmat, ncomp, info)
        if( info /= 0 ) THROW_HARD('lot_pullback_metric: ridge system not positive definite')
        allocate(mm(ncomp,ncomp), metric(ncomp,ncomp))
        call dgemm('N', 'T', ncomp, ncomp, nm, 1._dp, jmat, ncomp, jmat, ncomp, 0._dp, mm, ncomp)
        tr = 0._dp
        do q = 1, ncomp
            tr = tr + mm(q,q)
        end do
        sc = real(ncomp, dp) / max(tr, DTINY)          ! normalise so the metric does not rescale distances
        mm = sc*mm
        do r = 1, ncomp
            do q = 1, ncomp
                metric(q,r) = blend*mm(q,r)
            end do
            metric(r,r) = metric(r,r) + (1._dp - blend)
        end do
        write(logfhandle,'(A,I0,A,F5.2,A,ES10.3)') '>>> FLEX_LOT pullback metric: rank ', nm, &
            &', blend ', blend, ', ridge ', ridge
        call flush(logfhandle)
        deallocate(ts, gram, rhs, jmat, mm, ipiv)
    end subroutine lot_pullback_metric

end module simple_flex_lot
