!@descr: UMAP projection of a high-dimensional embedding, for plotting it
module simple_umap
use simple_core_module_api
use simple_kd_tree, only: kd_tree, knn_table
implicit none
private
#include "simple_local_flags.inc"

public :: umap_embed, umap_subsample

integer, parameter :: NDIM         = 2      !< output dimensionality
integer, parameter :: NNN          = 15     !< neighbourhood size
integer, parameter :: SMALL_N      = 10000  !< below this, more epochs (UMAP's own heuristic)
integer, parameter :: NEPOCH_SMALL = 500, NEPOCH_LARGE = 200
integer, parameter :: NBISECT      = 64     !< bandwidth bisection steps
real,    parameter :: NEG_RATE     = 5.     !< negative samples per positive
real,    parameter :: CLIP         = 4.     !< SGD displacement clamp; without it coincident pairs diverge
real,    parameter :: INIT_SCALE   = 10.    !< coordinate spread of the initialisation
!> a and b of the 1/(1 + a d^2b) output kernel, least-squares fitted to UMAP's defaults
!! min_dist=0.1, spread=1.0. Constants because nothing here varies them.
real,    parameter :: KERN_A = 1.577, KERN_B = 0.895

contains

    !>  UMAP of X(nsamples, ndims) into Y(nsamples, NDIM). Columns are z-scored first, so a
    !!  direction's influence comes from its structure and not from its variance. The kNN is the
    !!  repo's EXACT k-d tree, so only the SGD's negative samples carry randomness, and those come
    !!  from a stream local to this call -- a plot never perturbs the global one.
    subroutine umap_embed( X, seed, Y )
        real,              intent(in)    :: X(:,:)
        integer,           intent(in)    :: seed
        real, allocatable, intent(inout) :: Y(:,:)
        real,    allocatable :: W(:,:), feats(:,:), w_e(:), eps(:), next_smp(:), next_neg(:)
        integer, allocatable :: head(:), tail(:)
        type(kd_tree)   :: tree
        type(knn_table) :: knn
        integer(kind=longer) :: rng
        real    :: mu, sdev, wmax, alpha, d2, coeff, delta, ycur
        integer :: n, d, i, j, k, e, t, kk, nits, knn_k, nedges, nneg, epoch, other
        n = size(X,1)
        d = size(X,2)
        if( n < 4 .or. d < NDIM + 1 ) THROW_HARD('UMAP needs >3 samples in >2 dimensions')
        allocate(W(n,d), source=X)
        do j = 1, d
            mu     = sum(W(:,j))/real(n)
            W(:,j) = W(:,j) - mu
            sdev   = sqrt(max(sum(W(:,j)**2)/real(n - 1), 0.))
            if( sdev > TINY ) W(:,j) = W(:,j)/sdev
        end do
        nits = NEPOCH_SMALL
        if( n > SMALL_N ) nits = NEPOCH_LARGE
        knn_k = max(2, min(NNN, n - 1))
        allocate(feats(d,n))
        feats = transpose(W)                        ! the k-d tree wants one point per COLUMN
        call tree%build(feats)
        call tree%query_all(feats, knn_k, knn)
        call tree%kill
        deallocate(feats)
        call fuzzy_graph(knn, n, head, tail, w_e, nedges)
        call knn%kill
        if( nedges < 1 ) THROW_HARD('UMAP produced an empty neighbourhood graph')
        call pca_init(W, Y)
        rng = seed_rng(seed)
        do i = 1, n                                 ! duplicate rows coincide, and the gradient is
            do kk = 1, NDIM                         ! singular there, so break the ties
                Y(i,kk) = Y(i,kk) + 1.e-2*(2.*rnd_uni(rng) - 1.)
            end do
        end do
        ! an edge is visited a number of times proportional to its membership strength; one too weak
        ! for a single visit over the whole schedule is dropped, not rounded up to one
        wmax = maxval(w_e(1:nedges))
        allocate(eps(nedges), next_smp(nedges), next_neg(nedges))
        do e = 1, nedges
            eps(e) = -1.
            if( w_e(e)*real(nits) >= wmax ) eps(e) = wmax/max(w_e(e), TINY)
            next_smp(e) = eps(e)
            next_neg(e) = eps(e)/NEG_RATE
        end do
        ! serial by choice: the reference implementation races on Y, and reproducibility is worth
        ! more here than the wall clock a hogwild sweep would save
        do epoch = 1, nits
            alpha = 1. - real(epoch - 1)/real(nits)
            do e = 1, nedges
                if( eps(e) <= 0. .or. next_smp(e) > real(epoch) ) cycle
                i = head(e); j = tail(e)
                d2    = sum((Y(i,:) - Y(j,:))**2)
                coeff = 0.                          ! attraction: d/dy of log(1 + a d^2b)
                if( d2 > 0. ) coeff = -2.*KERN_A*KERN_B*(d2**(KERN_B - 1.))/(1. + KERN_A*(d2**KERN_B))
                do kk = 1, NDIM
                    delta   = max(-CLIP, min(CLIP, coeff*(Y(i,kk) - Y(j,kk))))
                    ycur    = Y(i,kk)
                    Y(i,kk) = ycur       + alpha*delta
                    Y(j,kk) = Y(j,kk)    - alpha*delta
                end do
                next_smp(e) = next_smp(e) + eps(e)
                nneg = int((real(epoch) - next_neg(e))*NEG_RATE/eps(e))
                do t = 1, nneg
                    other = rnd_below(rng, n)
                    if( other == i ) cycle
                    d2    = sum((Y(i,:) - Y(other,:))**2)
                    coeff = 0.                      ! repulsion
                    if( d2 > 0. ) coeff = 2.*KERN_B/((0.001 + d2)*(1. + KERN_A*(d2**KERN_B)))
                    do kk = 1, NDIM
                        delta = CLIP
                        if( coeff > 0. ) delta = max(-CLIP, min(CLIP, coeff*(Y(i,kk) - Y(other,kk))))
                        Y(i,kk) = Y(i,kk) + alpha*delta
                    end do
                end do
                next_neg(e) = next_neg(e) + real(nneg)*eps(e)/NEG_RATE
            end do
        end do
        do k = 1, NDIM
            Y(:,k) = Y(:,k) - sum(Y(:,k))/real(n)
        end do
        write(logfhandle,'(A,I0,A,I0)') '>>> UMAP epochs=', nits, ' edges=', nedges
        deallocate(W, w_e, eps, next_smp, next_neg, head, tail)
    end subroutine umap_embed

    !>  Deterministic size-capped row subset, ascending. Depends only on (n, nsub, seed).
    subroutine umap_subsample( n, nsub, seed, inds )
        integer,              intent(in)  :: n, nsub, seed
        integer, allocatable, intent(out) :: inds(:)
        integer, allocatable :: perm(:)
        integer(kind=longer) :: rng
        integer :: i, j, tmp
        if( n < 1 ) THROW_HARD('umap_subsample needs a positive population')
        if( nsub < 1 .or. n <= nsub )then
            allocate(inds(n))
            inds = [(i, i=1,n)]
            return
        endif
        allocate(perm(n))
        perm = [(i, i=1,n)]
        rng  = seed_rng(seed)
        do i = 1, nsub                              ! partial Fisher-Yates: only nsub draws needed
            j = i - 1 + rnd_below(rng, n - i + 1)
            tmp = perm(i); perm(i) = perm(j); perm(j) = tmp
        end do
        allocate(inds(nsub), source=perm(1:nsub))
        call hpsort(inds)
        deallocate(perm)
    end subroutine umap_subsample

    !>  Local kernels, then the probabilistic union of each pair's two directed memberships. The
    !!  reverse membership is looked up directly in the neighbour's own list -- O(k) with k a
    !!  handful -- so the edge set needs no sort and no deduplication pass.
    subroutine fuzzy_graph( knn, n, head, tail, w_e, nedges )
        type(knn_table),      intent(in)  :: knn
        integer,              intent(in)  :: n
        integer, allocatable, intent(out) :: head(:), tail(:)
        real,    allocatable, intent(out) :: w_e(:)
        integer,              intent(out) :: nedges
        real,    allocatable :: dists(:,:), wdir(:,:)
        real    :: rho, sigma, lo, hi, tgt, s, dj, wij, wji
        integer :: k, i, j, it, nb, pos
        k = knn%k
        allocate(dists(k,n), wdir(k,n), source=0.)
        dists = sqrt(max(knn%distance2, 0.))
        tgt   = log(real(k))/log(2.)
        !$omp parallel do default(shared) private(i,j,it,rho,sigma,lo,hi,s,dj) &
        !$omp schedule(static) proc_bind(close)
        do i = 1, n
            rho = 0.                                ! local connectivity 1: full membership to the
            do j = 1, k                             ! nearest distinct neighbour, which is what keeps
                if( dists(j,i) > TINY )then         ! the graph connected in sparse regions
                    rho = dists(j,i)
                    exit
                endif
            end do
            lo = 0.; hi = huge(0.); sigma = 1.
            do it = 1, NBISECT                      ! bisect so the row's memberships sum to log2(k)
                s = 0.
                do j = 1, k
                    s = s + exp(-max(dists(j,i) - rho, 0.)/sigma)
                end do
                if( abs(s - tgt) < 1.e-5 ) exit
                if( s > tgt )then
                    hi    = sigma
                    sigma = 0.5*(lo + hi)
                else
                    lo = sigma
                    if( hi > 0.5*huge(0.) )then
                        sigma = 2.*sigma
                    else
                        sigma = 0.5*(lo + hi)
                    endif
                endif
            end do
            do j = 1, k
                wdir(j,i) = exp(-max(dists(j,i) - rho, 0.)/sigma)
                if( knn%neighbor(j,i) == i ) wdir(j,i) = 0.
            end do
        end do
        !$omp end parallel do
        deallocate(dists)
        allocate(head(n*k), tail(n*k), w_e(n*k))
        nedges = 0
        do i = 1, n
            do j = 1, k
                nb  = knn%neighbor(j,i)
                if( nb == i ) cycle
                pos = nbr_pos(knn%neighbor(:,nb), i)
                ! emit each undirected pair once: from the lower endpoint when both directions are
                ! present, otherwise from whichever end holds the only one
                if( nb < i .and. pos > 0 ) cycle
                wji = 0.
                if( pos > 0 ) wji = wdir(pos,nb)
                wij = wdir(j,i)
                if( wij + wji <= TINY ) cycle
                nedges = nedges + 1
                head(nedges) = i
                tail(nedges) = nb
                w_e(nedges)  = wij + wji - wij*wji
            end do
        end do
        deallocate(wdir)
    end subroutine fuzzy_graph

    pure integer function nbr_pos( nbrs, want ) result( pos )
        integer, intent(in) :: nbrs(:), want
        integer :: j
        pos = 0
        do j = 1, size(nbrs)
            if( nbrs(j) == want )then
                pos = j
                return
            endif
        end do
    end function nbr_pos

    !>  PCA initialisation, scaled so the widest axis has standard deviation INIT_SCALE. Better
    !!  global arrangement than a random start, and reproducible. The d x d covariance route avoids
    !!  ever forming the n x n matrix; d is the latent rank, so it is small.
    subroutine pca_init( W, Y )
        real,              intent(in)    :: W(:,:)
        real, allocatable, intent(inout) :: Y(:,:)
        real(dp), allocatable :: C(:,:), V(:,:), eigs(:)
        real    :: sdev, smax
        integer :: n, d, i, j, q, nrot
        n = size(W,1)
        d = size(W,2)
        if( allocated(Y) ) deallocate(Y)
        allocate(Y(n,NDIM), source=0.)
        allocate(C(d,d), V(d,d), source=0.d0)
        allocate(eigs(d), source=0.d0)
        do j = 1, d
            do i = 1, j
                C(i,j) = sum(real(W(:,i), dp)*real(W(:,j), dp))/real(n - 1, dp)
                C(j,i) = C(i,j)
            end do
        end do
        nrot = 0
        call jacobi(C, d, d, eigs, V, nrot)
        call eigsrt(eigs, V, d, d)
        do q = 1, NDIM
            do i = 1, n
                Y(i,q) = real(sum(real(W(i,:), dp)*V(:,q)), sp)
            end do
        end do
        smax = 0.
        do q = 1, NDIM
            sdev = sqrt(max(sum(Y(:,q)**2)/real(n), 0.))
            smax = max(smax, sdev)
        end do
        if( smax > TINY ) Y = Y*(INIT_SCALE/smax)
        deallocate(C, V, eigs)
    end subroutine pca_init

    ! ---- xorshift64, kept local so a plot never perturbs the global random stream

    integer(kind=longer) function seed_rng( seed ) result( state )
        integer, intent(in) :: seed
        integer :: i
        state = ieor(88172645463325252_longer, ishft(int(seed, longer), 21))
        state = ieor(state, int(seed, longer))
        if( state == 0_longer ) state = 88172645463325252_longer
        do i = 1, 8                                 ! a low-entropy seed shows through the first words
            state = ieor(state, ishft(state,  13))
            state = ieor(state, ishft(state,  -7))
            state = ieor(state, ishft(state,  17))
        end do
    end function seed_rng

    real function rnd_uni( state ) result( r )
        integer(kind=longer), intent(inout) :: state
        state = ieor(state, ishft(state,  13))
        state = ieor(state, ishft(state,  -7))
        state = ieor(state, ishft(state,  17))
        r = real(iand(ishft(state, -24), 16777215_longer), sp)/16777216.0_sp
    end function rnd_uni

    integer function rnd_below( state, n ) result( i )
        integer(kind=longer), intent(inout) :: state
        integer,              intent(in)    :: n
        i = max(1, min(n, 1 + int(rnd_uni(state)*real(n))))
    end function rnd_below

end module simple_umap
