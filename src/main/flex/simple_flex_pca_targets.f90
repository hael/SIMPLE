!@descr: flex_pca latent targets: k-means, diffusion k-centre, FINCH, path and reliability-path placement; basis rotations
module simple_flex_pca_targets
use simple_core_module_api
use simple_flex_pca_em, only: cov_env_int_pub
use simple_image, only: image
use simple_parameters, only: parameters
use simple_srch_sort_loc, only: hpsort
use simple_finch, only: finch_hierarchy, fit_finch, finch_representatives, select_finch_level, refine_finch_level
use simple_kd_tree, only: kd_tree, knn_table
use simple_linalg, only: jacobi, eigsrt, matinv
implicit none
private
#include "simple_local_flags.inc"

public :: diffusion_kcenter_targets, kmeans_latent_targets, path_latent_targets, reliability_path_targets, component_reliability_proxy, invert_lower, orth_block, sort_block_desc


contains

    !> Manifold-covering state targets by diffusion-map k-center. k-means allocates by DENSITY and misses
    !! sparse states; a 1-D path merges states on a BRANCHED manifold. A diffusion embedding makes geodesic
    !! structure Euclidean, so greedy farthest-point covers any shape. Reliability weighting is essential:
    !! the largest-eigenvalue component is typically the worst measured.
    !! Measurements: doc/implementation_notes/flex_pca_state_placement_measurements.md
    subroutine diffusion_kcenter_targets( z, nptcls, ncomp, nstates, wcomp, rho, centroids, ok )
        integer,  intent(in)  :: nptcls, ncomp, nstates
        real(dp), intent(in)  :: z(nptcls,ncomp), wcomp(ncomp), rho(ncomp)
        real(dp), intent(out) :: centroids(ncomp,nstates)
        logical,  intent(out) :: ok
        integer,  parameter :: NNODE_MAX = 20000   ! graph nodes; kd-tree kNN keeps this affordable
        integer,  parameter :: KNN       = 24      ! neighbours per node in the affinity graph
        integer,  parameter :: NDIFF     = 2       ! non-trivial diffusion coordinates retained
        integer,  parameter :: NPOW      = 300     ! orthogonal-iteration sweeps
        real(dp), parameter :: RHO_FLOOR = 0.1d0
        type(kd_tree)   :: tree
        type(knn_table) :: knntab
        real,     allocatable :: feats(:,:)
        integer,  allocatable :: nodes(:), er(:), ec(:), sel(:), cell(:), ccnt(:)
        real(dp), allocatable :: ev(:), sig(:), qdeg(:), ddeg(:), V(:,:), W(:,:), lam(:)
        real(dp), allocatable :: psi(:,:), dmin(:), rw(:), zbar(:), sdv(:)
        real(dp) :: d2, s, nrm, best, wk1
        integer  :: nnode, i, j, q, e, nedge, it, m, c, ibest, ni
        ok = .false.
        if( nstates < 2 .or. nptcls < 100 ) return
        nnode = min(nptcls, NNODE_MAX)
        if( nnode <= KNN + 2 ) return
        m = NDIFF + 1                                  ! + the trivial eigenvector
        allocate(nodes(nnode), rw(ncomp), zbar(ncomp), sdv(ncomp))
        ! deterministic stride subsample: it already carries the data's density, and is reproducible run to run
        do i = 1, nnode
            nodes(i) = 1 + int(real(i-1,dp)*real(nptcls-1,dp)/real(max(1,nnode-1),dp))
        end do
        do q = 1, ncomp
            zbar(q) = sum(z(:,q)) / real(nptcls,dp)
            sdv(q)  = sqrt(max(sum((z(:,q)-zbar(q))**2)/real(nptcls,dp), DTINY))
        end do
        s = maxval(rho)
        if( s <= DTINY ) s = 1.d0
        do q = 1, ncomp
            ! wcomp is 1/var per component, so sqrt(wcomp) ALREADY carries the 1/sd standardisation.
            ! Dividing by sdv again standardises twice and flattens the geometry the graph reads.
            rw(q) = sqrt(max(wcomp(q),0.d0)) * max(rho(q)/s, RHO_FLOOR)
        end do
        allocate(feats(ncomp,nnode))
        do i = 1, nnode
            do q = 1, ncomp
                feats(q,i) = real((z(nodes(i),q) - zbar(q)) * rw(q))
            end do
        end do
        call tree%build(feats)
        call tree%query_all(feats, KNN, knntab)
        ! self-tuning bandwidth: sigma_i is the K-th neighbour distance, so affinity adapts to local density
        allocate(sig(nnode))
        do i = 1, nnode
            sig(i) = sqrt(max(real(knntab%distance2(KNN,i),dp), DTINY))
        end do
        nedge = 2*nnode*KNN
        allocate(er(nedge), ec(nedge), ev(nedge))
        e = 0
        do i = 1, nnode
            do j = 1, KNN
                ni = knntab%neighbor(j,i)
                if( ni < 1 .or. ni > nnode .or. ni == i ) cycle
                d2 = real(knntab%distance2(j,i),dp)
                s  = exp(-d2/max(sig(i)*sig(ni), DTINY))
                e = e + 1; er(e) = i;  ec(e) = ni; ev(e) = s
                e = e + 1; er(e) = ni; ec(e) = i;  ev(e) = s
            end do
        end do
        nedge = e
        if( nedge < nnode ) return
        ! alpha = 1 (Laplace-Beltrami): divide out the sampling density so the embedding reflects manifold
        ! GEOMETRY, not how heavily each region is populated -- the exact failure mode of k-means here.
        allocate(qdeg(nnode), source=0.d0)
        do e = 1, nedge
            qdeg(er(e)) = qdeg(er(e)) + ev(e)
        end do
        do e = 1, nedge
            ev(e) = ev(e) / max(qdeg(er(e))*qdeg(ec(e)), DTINY)
        end do
        allocate(ddeg(nnode), source=0.d0)
        do e = 1, nedge
            ddeg(er(e)) = ddeg(er(e)) + ev(e)
        end do
        do i = 1, nnode
            ddeg(i) = 1.d0 / sqrt(max(ddeg(i), DTINY))
        end do
        do e = 1, nedge
            ev(e) = ev(e) * ddeg(er(e)) * ddeg(ec(e))       ! S = D^-1/2 W D^-1/2, symmetric
        end do
        ! leading eigenvectors of S by orthogonal iteration (deterministic start)
        allocate(V(nnode,m), W(nnode,m), lam(m))
        do j = 1, m
            do i = 1, nnode
                V(i,j) = sin(real(i*j,dp)*0.7717d0) + 0.1d0*real(j,dp)
            end do
        end do
        call orth_block(V, nnode, m)
        do it = 1, NPOW
            W = 0.d0
            do e = 1, nedge
                do j = 1, m
                    W(er(e),j) = W(er(e),j) + ev(e)*V(ec(e),j)
                end do
            end do
            V = W
            call orth_block(V, nnode, m)
        end do
        W = 0.d0
        do e = 1, nedge
            do j = 1, m
                W(er(e),j) = W(er(e),j) + ev(e)*V(ec(e),j)
            end do
        end do
        do j = 1, m
            lam(j) = sum(V(:,j)*W(:,j))
        end do
        call sort_block_desc(V, lam, nnode, m)
        ! psi = D^-1/2 V, dropping the trivial leading eigenvector; the 1/sqrt(1-lambda) commute-time
        ! scaling puts Euclidean distance on the diffusion metric, which is what makes k-center meaningful.
        allocate(psi(nnode,NDIFF))
        wk1 = 1.d0 / sqrt(max(1.d0 - min(lam(2), 1.d0-1.d-9), 1.d-9))
        do j = 1, NDIFF
            s = 1.d0 / sqrt(max(1.d0 - min(lam(j+1), 1.d0-1.d-9), 1.d-9))
            nrm = sqrt(max(sum((V(:,j+1)*ddeg)**2), DTINY))
            do i = 1, nnode
                psi(i,j) = V(i,j+1)*ddeg(i)/nrm * (s/wk1)
            end do
        end do
        ! greedy k-center: seed farthest from the centroid, then farthest from everything chosen -- coverage
        allocate(sel(nstates), dmin(nnode))
        do j = 1, NDIFF
            s = sum(psi(:,j))/real(nnode,dp)
            psi(:,j) = psi(:,j) - s
        end do
        best = -1.d0; ibest = 1
        do i = 1, nnode
            d2 = sum(psi(i,:)**2)
            if( d2 > best )then
                best = d2; ibest = i
            endif
        end do
        sel(1) = ibest
        do i = 1, nnode
            dmin(i) = sum((psi(i,:)-psi(ibest,:))**2)
        end do
        do c = 2, nstates
            best = -1.d0; ibest = 1
            do i = 1, nnode
                if( dmin(i) > best )then
                    best = dmin(i); ibest = i
                endif
            end do
            sel(c) = ibest
            do i = 1, nnode
                dmin(i) = min(dmin(i), sum((psi(i,:)-psi(ibest,:))**2))
            end do
        end do
        ! each node takes the nearest selected node's CELL; the target is the cell's mean RAW latent, which
        ! regresses the RIM nodes k-center picks by design back onto the state they cover.
        allocate(cell(nnode), ccnt(nstates), source=0)
        do i = 1, nnode
            best = huge(1.d0); ibest = 1
            do c = 1, nstates
                d2 = sum((psi(i,:) - psi(sel(c),:))**2)
                if( d2 < best )then
                    best = d2; ibest = c
                endif
            end do
            cell(i) = ibest
        end do
        centroids = 0.d0
        do i = 1, nnode
            c = cell(i)
            ccnt(c) = ccnt(c) + 1
            centroids(:,c) = centroids(:,c) + z(nodes(i),:)
        end do
        do c = 1, nstates
            if( ccnt(c) > 0 )then
                centroids(:,c) = centroids(:,c) / real(ccnt(c),dp)
            else
                centroids(:,c) = z(nodes(sel(c)),:)     ! degenerate cell: keep the node itself
            endif
        end do
        write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA diffusion cell occupancy: min=',minval(ccnt), &
            &' max=',maxval(ccnt)
        write(logfhandle,'(A,I0,A,I0,A,F7.4,A,F7.4)') '>>> FLEX_PCA diffusion k-center: nodes=',nnode, &
            &' knn=',KNN,' lambda2=',lam(2),' lambda3=',lam(3)
        ok = .true.
        call tree%kill; call knntab%kill
        deallocate(feats, nodes, er, ec, ev, sig, qdeg, ddeg, V, W, lam, psi, dmin, sel, rw, zbar, sdv, cell, ccnt)
    end subroutine diffusion_kcenter_targets

    !> Deterministic k-means in the SAME wcomp metric the kernel uses, so placement and weighting agree.
    !! Seeded farthest-point from the particle nearest the latent mean, so no RNG is involved.
    subroutine kmeans_latent_targets( z, nptcls, ncomp, nstates, wcomp, centroids )
        integer,  intent(in)  :: nptcls, ncomp, nstates
        real(dp), intent(in)  :: z(nptcls,ncomp), wcomp(ncomp)
        real(dp), intent(out) :: centroids(ncomp,nstates)
        integer, parameter    :: MAXIT = 50
        real(dp), allocatable :: mind(:), csum(:,:), zbar(:)
        integer,  allocatable :: cnt(:), memb(:)
        real(dp) :: d2, best, dmax
        integer  :: i, q, s, it, ibest, iseed, nchanged
        logical  :: l_reseed
        allocate(mind(nptcls), csum(ncomp,nstates), zbar(ncomp), cnt(nstates), memb(nptcls))
        do q = 1, ncomp
            zbar(q) = sum(z(:,q)) / real(nptcls,dp)
        end do
        best = huge(1.d0); iseed = 1
        do i = 1, nptcls
            d2 = 0.d0
            do q = 1, ncomp
                d2 = d2 + wcomp(q)*(z(i,q)-zbar(q))**2
            end do
            if( d2 < best )then
                best  = d2
                iseed = i
            endif
        end do
        centroids(:,1) = z(iseed,:)
        mind = huge(1.d0)
        do s = 2, nstates
            dmax = -1.d0; iseed = 1
            do i = 1, nptcls
                d2 = 0.d0
                do q = 1, ncomp
                    d2 = d2 + wcomp(q)*(z(i,q)-centroids(q,s-1))**2
                end do
                mind(i) = min(mind(i), d2)
                if( mind(i) > dmax )then
                    dmax  = mind(i)
                    iseed = i
                endif
            end do
            centroids(:,s) = z(iseed,:)
        end do
        memb = 0
        do it = 1, MAXIT
            nchanged = 0
            !$omp parallel do default(shared) private(i,q,s,d2,best,ibest) schedule(static) reduction(+:nchanged)
            do i = 1, nptcls
                best = huge(1.d0); ibest = 1
                do s = 1, nstates
                    d2 = 0.d0
                    do q = 1, ncomp
                        d2 = d2 + wcomp(q)*(z(i,q)-centroids(q,s))**2
                    end do
                    if( d2 < best )then
                        best  = d2
                        ibest = s
                    endif
                end do
                if( memb(i) /= ibest ) nchanged = nchanged + 1
                memb(i) = ibest
                mind(i)   = best
            end do
            !$omp end parallel do
            csum = 0.d0; cnt = 0
            do i = 1, nptcls
                cnt(memb(i))    = cnt(memb(i)) + 1
                csum(:,memb(i)) = csum(:,memb(i)) + z(i,:)
            end do
            l_reseed = .false.
            do s = 1, nstates
                if( cnt(s) > 0 )then
                    centroids(:,s) = csum(:,s) / real(cnt(s),dp)
                else
                    ! empty cluster: reseed on the worst-fitted particle
                    iseed          = maxloc(mind, dim=1)
                    centroids(:,s) = z(iseed,:)
                    mind(iseed)    = -1.d0
                    l_reseed       = .true.
                endif
            end do
            if( nchanged == 0 .and. .not. l_reseed ) exit
        end do
        write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA k-means iterations=',min(it,MAXIT), &
            &'  reassigned on last pass=',nchanged
        deallocate(mind, csum, zbar, cnt, memb)
    end subroutine kmeans_latent_targets

    subroutine path_latent_targets( z, nptcls, ncomp, nstates, wcomp, centroids )
        integer,  intent(in)  :: nptcls, ncomp, nstates
        real(dp), intent(in)  :: z(nptcls,ncomp), wcomp(ncomp)
        real(dp), intent(out) :: centroids(ncomp,nstates)
        real(dp), allocatable :: zbar(:), pa(:), pb(:), dirv(:), proj(:)
        real,     allocatable :: sproj(:)
        integer,  allocatable :: cnt(:)
        real(dp) :: d2, dmax, dnorm, lo, hi, t
        integer  :: i, q, s, ia, ib, islot
        allocate(zbar(ncomp), pa(ncomp), pb(ncomp), dirv(ncomp), proj(nptcls), &
            &sproj(nptcls), cnt(nstates))
        do q = 1, ncomp
            zbar(q) = sum(z(:,q)) / real(nptcls,dp)
        end do
        dmax = -1.d0; ia = 1
        do i = 1, nptcls
            d2 = 0.d0
            do q = 1, ncomp
                d2 = d2 + wcomp(q)*(z(i,q)-zbar(q))**2
            end do
            if( d2 > dmax )then
                dmax = d2; ia = i
            endif
        end do
        pa = z(ia,:)
        dmax = -1.d0; ib = 1
        do i = 1, nptcls
            d2 = 0.d0
            do q = 1, ncomp
                d2 = d2 + wcomp(q)*(z(i,q)-pa(q))**2
            end do
            if( d2 > dmax )then
                dmax = d2; ib = i
            endif
        end do
        pb = z(ib,:)
        dirv = pb - pa
        dnorm = sum(wcomp*dirv*dirv)
        if( dnorm <= DTINY )then
            do s = 1, nstates
                centroids(:,s) = zbar
            end do
            deallocate(zbar, pa, pb, dirv, proj, sproj, cnt)
            return
        endif
        ! project every particle onto the segment, in the same metric the kernel uses
        !$omp parallel do default(shared) private(i,q,d2) schedule(static)
        do i = 1, nptcls
            d2 = 0.d0
            do q = 1, ncomp
                d2 = d2 + wcomp(q)*(z(i,q)-pa(q))*dirv(q)
            end do
            proj(i) = d2 / dnorm
        end do
        !$omp end parallel do
        ! equal-WIDTH slices between the 0.1 % and 99.9 % projections, each target the slice mean
        sproj = real(proj)
        call hpsort(sproj)
        centroids = 0.d0; cnt = 0
        do i = 1, nptcls
            lo = real(sproj(max(1,min(nptcls,1+int(real(nptcls,dp)*0.001d0)))),dp)
            hi = real(sproj(max(1,min(nptcls,nptcls-int(real(nptcls,dp)*0.001d0)))),dp)
            if( hi-lo <= DTINY )then
                islot = 1
            else
                t     = (proj(i)-lo)/(hi-lo)
                islot = 1 + int(t*real(nstates,dp))
                islot = max(1, min(nstates, islot))
            endif
            cnt(islot) = cnt(islot) + 1
            do q = 1, ncomp
                centroids(q,islot) = centroids(q,islot) + z(i,q)
            end do
        end do
        do s = 1, nstates
            if( cnt(s) > 0 )then
                centroids(:,s) = centroids(:,s) / real(cnt(s),dp)
            else
                t = real(s-1,dp)/real(max(1,nstates-1),dp)
                centroids(:,s) = pa + t*dirv
            endif
        end do
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA path endpoints from particles ',ia,' and ',ib, &
            &'; slice occupancies min=',minval(cnt)
        deallocate(zbar, pa, pb, dirv, proj, sproj, cnt)
    end subroutine path_latent_targets

    !> Equal-occupancy targets along a reliability-ordered principal direction. Two departures from
    !! k-means: the direction is the RELIABILITY-weighted leading eigenvector, so a high-variance but
    !! poorly measured nuisance mode cannot set it; and slices carry equal PARTICLE COUNTS, not width.
    subroutine reliability_path_targets( z, nptcls, ncomp, nstates, wcomp, rho, centroids, proj_out, tproj_out )
        integer,  intent(in)  :: nptcls, ncomp, nstates
        real(dp), intent(in)  :: z(nptcls,ncomp), wcomp(ncomp), rho(ncomp)
        real(dp), intent(out) :: centroids(ncomp,nstates)
        ! coordinate ALONG the path. A full-rank Mahalanobis kernel also measures the off-path directions, so
        ! an on-path but noisy particle falls outside every support -- that is what strands the dataset.
        real(dp), optional, intent(out) :: proj_out(nptcls), tproj_out(nstates)
        integer,  parameter   :: NPOWER = 128
        real(dp), parameter   :: RHO_PATH_FLOOR = 0.1d0
        real(dp), allocatable :: zbar(:), u(:), unew(:), cov(:,:), proj(:), rw(:), zc(:,:), pedge(:)
        real,     allocatable :: sproj(:)
        integer,  allocatable :: cnt(:)
        real(dp) :: rmax, nrm, d
        integer  :: i, q, s, it, ilo, islot, slo, shi
        allocate(zbar(ncomp), u(ncomp), unew(ncomp), cov(ncomp,ncomp), proj(nptcls), &
            &sproj(nptcls), cnt(nstates), rw(ncomp), zc(nptcls,ncomp), pedge(max(1,nstates-1)))
        do q = 1, ncomp
            zbar(q) = sum(z(:,q)) / real(nptcls,dp)
        end do
        ! reliability RELATIVE to the best-measured component, floored so no direction is removed outright
        rmax = maxval(rho)
        if( rmax <= DTINY ) rmax = 1.d0
        do q = 1, ncomp
            rw(q) = sqrt(max(wcomp(q),0.d0)) * max(rho(q)/rmax, RHO_PATH_FLOOR)
        end do
        !$omp parallel do default(shared) private(i,q) schedule(static)
        do i = 1, nptcls
            do q = 1, ncomp
                zc(i,q) = rw(q)*(z(i,q)-zbar(q))
            end do
        end do
        !$omp end parallel do
        cov = matmul(transpose(zc), zc)
        ! leading eigenvector by power iteration: deterministic, and ncomp is at most a few dozen
        u = 1.d0 / sqrt(real(ncomp,dp))
        do it = 1, NPOWER
            unew = matmul(cov, u)
            nrm  = sqrt(sum(unew*unew))
            if( nrm <= DTINY ) exit
            unew = unew / nrm
            if( sum(abs(unew-u)) < 1.d-12 )then
                u = unew
                exit
            endif
            u = unew
        end do
        !$omp parallel do default(shared) private(i,q,d) schedule(static)
        do i = 1, nptcls
            d = 0.d0
            do q = 1, ncomp
                d = d + zc(i,q)*u(q)
            end do
            proj(i) = d
        end do
        !$omp end parallel do
        sproj = real(proj)
        call hpsort(sproj)
        do s = 1, nstates-1
            ilo = max(1, min(nptcls, nint(real(s,dp)/real(nstates,dp)*real(nptcls,dp))))
            pedge(s) = real(sproj(ilo),dp)
        end do
        centroids = 0.d0; cnt = 0
        do i = 1, nptcls
            islot = nstates
            do s = 1, nstates-1
                if( proj(i) < pedge(s) )then
                    islot = s
                    exit
                endif
            end do
            cnt(islot) = cnt(islot) + 1
            do q = 1, ncomp
                centroids(q,islot) = centroids(q,islot) + z(i,q)
            end do
        end do
        do s = 1, nstates
            if( cnt(s) > 0 )then
                centroids(:,s) = centroids(:,s) / real(cnt(s),dp)
            endif
        end do
        ! ties can leave a slice empty; interpolate so the polyline stays ordered and no state hits the mean
        do s = 1, nstates
            if( cnt(s) > 0 ) cycle
            slo = 0; shi = 0
            do i = s-1, 1, -1
                if( cnt(i) > 0 )then
                    slo = i; exit
                endif
            end do
            do i = s+1, nstates
                if( cnt(i) > 0 )then
                    shi = i; exit
                endif
            end do
            if( slo > 0 .and. shi > 0 )then
                d = real(s-slo,dp)/real(shi-slo,dp)
                centroids(:,s) = (1.d0-d)*centroids(:,slo) + d*centroids(:,shi)
            else if( slo > 0 )then
                centroids(:,s) = centroids(:,slo)
            else if( shi > 0 )then
                centroids(:,s) = centroids(:,shi)
            else
                centroids(:,s) = zbar
            endif
        end do
        if( present(proj_out) ) proj_out = proj
        if( present(tproj_out) )then
            tproj_out = 0.d0
            do i = 1, nptcls
                islot = nstates
                do s = 1, nstates-1
                    if( proj(i) < pedge(s) )then
                        islot = s
                        exit
                    endif
                end do
                tproj_out(islot) = tproj_out(islot) + proj(i)
            end do
            do s = 1, nstates
                if( cnt(s) > 0 )then
                    tproj_out(s) = tproj_out(s) / real(cnt(s),dp)
                else if( s > 1 )then
                    tproj_out(s) = tproj_out(s-1)
                endif
            end do
        endif
        write(logfhandle,'(A,I0,A,I0,A,F6.3)') '>>> FLEX_PCA path ordering direction over ',ncomp, &
            &' components; slice occupancies min=',minval(cnt),'  leading |u| on z1=',abs(u(1))
        deallocate(zbar, u, unew, cov, proj, sproj, cnt, rw, zc, pedge)
    end subroutine reliability_path_targets

    !> Reliability proxy from a cached embedding: observed spread over mean posterior variance, mapped
    !! through r/(1+r) onto the split-half rho scale. Posterior variances are stride-sampled.
    subroutine component_reliability_proxy( z, precision, nptcls, ncomp, rho )
        integer,  intent(in)  :: nptcls, ncomp
        real(dp), intent(in)  :: z(nptcls,ncomp), precision(ncomp,ncomp,nptcls)
        real(dp), intent(out) :: rho(ncomp)
        integer,  parameter   :: NSAMPLE = 2000
        real(dp), allocatable :: cfull(:,:), pvar(:)
        real(dp) :: zbar, spread, ratio
        integer  :: i, q, errflg, stride, nused
        allocate(cfull(ncomp,ncomp), pvar(ncomp), source=0.d0)
        stride = max(1, nptcls/NSAMPLE)
        nused  = 0
        do i = 1, nptcls, stride
            call matinv(precision(:,:,i), cfull, ncomp, errflg)
            if( errflg /= 0 ) cycle
            do q = 1, ncomp
                pvar(q) = pvar(q) + max(cfull(q,q), 0.d0)
            end do
            nused = nused + 1
        end do
        if( nused > 0 ) pvar = pvar / real(nused,dp)
        do q = 1, ncomp
            zbar   = sum(z(:,q)) / real(nptcls,dp)
            spread = sum((z(:,q)-zbar)**2) / real(nptcls,dp)
            ratio  = spread / max(pvar(q), DTINY)
            rho(q) = ratio / (1.d0 + ratio)
        end do
        deallocate(cfull, pvar)
    end subroutine component_reliability_proxy


    subroutine invert_lower( L, n )
        integer,  intent(in)    :: n
        real(dp), intent(inout) :: L(n,n)
        real(dp), allocatable :: X(:,:)
        integer  :: i, j, k
        real(dp) :: s
        allocate(X(n,n), source=0.d0)
        do j = 1, n
            X(j,j) = 1.d0 / L(j,j)
            do i = j+1, n
                s = 0.d0
                do k = j, i-1
                    s = s + L(i,k)*X(k,j)
                end do
                X(i,j) = -s / L(i,i)
            end do
        end do
        L = X
        deallocate(X)
    end subroutine invert_lower

    subroutine orth_block( V, n, m )
        integer,  intent(in)    :: n, m
        real(dp), intent(inout) :: V(n,m)
        real(dp) :: dot, nrm
        integer  :: j, k
        do j = 1, m
            do k = 1, j-1
                dot = sum(V(:,k)*V(:,j))
                V(:,j) = V(:,j) - dot*V(:,k)
            end do
            nrm = sqrt(max(sum(V(:,j)**2), DTINY))
            V(:,j) = V(:,j) / nrm
        end do
    end subroutine orth_block

    subroutine sort_block_desc( V, lam, n, m )
        integer,  intent(in)    :: n, m
        real(dp), intent(inout) :: V(n,m), lam(m)
        real(dp) :: tl
        real(dp), allocatable :: tv(:)
        integer  :: i, j, imax
        allocate(tv(n))
        do i = 1, m-1
            imax = i
            do j = i+1, m
                if( lam(j) > lam(imax) ) imax = j
            end do
            if( imax /= i )then
                tl = lam(i); lam(i) = lam(imax); lam(imax) = tl
                tv = V(:,i); V(:,i) = V(:,imax); V(:,imax) = tv
            endif
        end do
        deallocate(tv)
    end subroutine sort_block_desc

end module simple_flex_pca_targets
