!@descr: CSR graph construction helpers for diffusion-map class splitting
module simple_diff_map_graphs
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
!$ use omp_lib
use simple_ori,                 only: ori
use simple_parameters,          only: parameters
use simple_srch_sort_loc,       only: hpsort
use simple_sp_project,          only: sp_project
implicit none
#include "simple_local_flags.inc"

private

public :: diffmap_graph
public :: build_cls_split_graph
public :: build_euclidean_knn_graph
public :: build_gated_euclidean_knn_graph
public :: find_gated_euclidean_neighbors_rows
public :: build_gated_euclidean_graph_from_neighbors
public :: build_orientation_knn_graph
public :: graph_matvec
public :: graph_directed_edges

type :: diffmap_graph
    integer :: n   = 0
    integer :: nnz = 0
    integer :: k_nn = 0
    character(len=16) :: metric   = 'euc'
    integer, allocatable :: rowptr(:)
    integer, allocatable :: colind(:)
    real,    allocatable :: w(:)
    real,    allocatable :: wnorm(:)
contains
    procedure :: kill      => kill_diffmap_graph
    procedure :: normalize => normalize_diffmap_graph
    procedure :: degree    => diffmap_graph_degree
end type diffmap_graph

contains

    subroutine build_cls_split_graph( params, spproj, pinds, pcavecs, graph )
        type(parameters),           intent(in)    :: params
        type(sp_project),           intent(inout) :: spproj
        integer,                    intent(in)    :: pinds(:)
        real, optional,             intent(in)    :: pcavecs(:,:)
        type(diffmap_graph),        intent(out)   :: graph
        character(len=16) :: metric
        metric = lowercase(trim(params%graph))
        if( trim(metric) /= 'euc' .and. trim(metric) /= 'ori' )then
            THROW_HARD('graph must be euc or ori in build_cls_split_graph')
        endif
        select case(trim(metric))
            case('ori')
                call build_orientation_knn_graph(params, spproj, pinds, max(2, params%k_nn), graph, &
                    &bandwidth_mode=params%bandwidth_mode, bandwidth_tune=params%bandwidth_tune)
            case DEFAULT
                if( .not. present(pcavecs) ) THROW_HARD('Euclidean graph requires pcavecs')
                call build_euclidean_knn_graph(pcavecs, max(2, params%k_nn), graph, &
                    &bandwidth_mode=params%bandwidth_mode, bandwidth_tune=params%bandwidth_tune)
        end select
    end subroutine build_cls_split_graph

    subroutine build_euclidean_knn_graph( pcavecs, k_nn, graph, bandwidth_mode, bandwidth_tune )
        real,                 intent(in)  :: pcavecs(:,:)
        integer,              intent(in)  :: k_nn
        type(diffmap_graph),  intent(out) :: graph
        character(len=*), optional, intent(in) :: bandwidth_mode
        real, optional, intent(in) :: bandwidth_tune
        integer, allocatable :: nbrs(:,:)
        real,    allocatable :: d2s(:,:), kth_d2(:)
        integer :: n, k_used
        n = size(pcavecs, 2)
        if( n < 1 ) THROW_HARD('empty Euclidean graph')
        if( n == 1 )then
            call make_singleton_graph('euc', graph)
            return
        endif
        k_used = min(max(1, k_nn), n - 1)
        allocate(nbrs(k_used,n), source=0)
        allocate(d2s(k_used,n), kth_d2(n), source=0.)
        call find_euclidean_neighbors(pcavecs, k_used, nbrs, d2s)
        kth_d2 = d2s(k_used,:)
        call pack_scalar_knn_to_csr(n, k_used, nbrs, d2s, kth_d2, 'euc', graph, bandwidth_mode, bandwidth_tune)
        deallocate(nbrs, d2s, kth_d2)
    end subroutine build_euclidean_knn_graph

    !> Build a Euclidean kNN graph without an all-pairs particle comparison.
    !! Projection directions define candidate bins. Neighboring bins are
    !! visited in angular order until nang_nbrs particle candidates have been
    !! examined, and only the k_nn closest registered-residuals are retained.
    subroutine build_gated_euclidean_knn_graph( features, proj_ids, proj_dirs, k_nn, nang_nbrs, graph, &
        &ncandidates_min, ncandidates_max, ncandidates_mean, bandwidth_mode, bandwidth_tune )
        real,                intent(in)  :: features(:,:)
        integer,             intent(in)  :: proj_ids(:)
        real,                intent(in)  :: proj_dirs(:,:)
        integer,             intent(in)  :: k_nn, nang_nbrs
        type(diffmap_graph), intent(out) :: graph
        integer, optional,   intent(out) :: ncandidates_min, ncandidates_max
        real, optional,      intent(out) :: ncandidates_mean
        character(len=*), optional, intent(in) :: bandwidth_mode
        real, optional, intent(in) :: bandwidth_tune
        integer, allocatable :: nbrs(:,:),ncandidates(:),rows(:)
        real, allocatable :: d2s(:,:)
        integer :: n, ndim, nproj, k_used, cap_used
        integer :: i
        n     = size(features,2)
        ndim  = size(features,1)
        nproj = size(proj_dirs,2)
        if( n < 1 .or. ndim < 1 ) THROW_HARD('empty gated Euclidean graph')
        if( size(proj_ids) /= n ) THROW_HARD('projection-id size mismatch in gated Euclidean graph')
        if( size(proj_dirs,1) /= 3 .or. nproj < 1 ) THROW_HARD('invalid projection directions in gated Euclidean graph')
        if( any(proj_ids < 1) .or. any(proj_ids > nproj) ) THROW_HARD('projection id outside direction table')
        if( n == 1 )then
            call make_singleton_graph('euc_gated', graph)
            if( present(ncandidates_min)  ) ncandidates_min  = 0
            if( present(ncandidates_max)  ) ncandidates_max  = 0
            if( present(ncandidates_mean) ) ncandidates_mean = 0.
            return
        endif
        k_used=min(max(1,k_nn),n-1); cap_used=min(max(k_used,nang_nbrs),n-1)
        allocate(rows(n)); rows=[(i,i=1,n)]
        call find_gated_euclidean_neighbors_rows(features,proj_ids,proj_dirs,k_used,cap_used,rows,nbrs,d2s,ncandidates)
        call build_gated_euclidean_graph_from_neighbors(n,nbrs,d2s,ncandidates,graph,bandwidth_mode,bandwidth_tune)
        if( present(ncandidates_min)  ) ncandidates_min  = minval(ncandidates)
        if( present(ncandidates_max)  ) ncandidates_max  = maxval(ncandidates)
        if( present(ncandidates_mean) ) ncandidates_mean = real(sum(int(ncandidates,kind=8)),kind=sp) / real(n,kind=sp)
        deallocate(rows,nbrs,d2s,ncandidates)
    end subroutine build_gated_euclidean_knn_graph

    subroutine find_gated_euclidean_neighbors_rows( features, proj_ids, proj_dirs, k_nn, nang_nbrs, rows, &
        &nbrs, d2s, ncandidates )
        real, intent(in) :: features(:,:),proj_dirs(:,:)
        integer, intent(in) :: proj_ids(:),k_nn,nang_nbrs,rows(:)
        integer, allocatable, intent(out) :: nbrs(:,:),ncandidates(:)
        real, allocatable, intent(out) :: d2s(:,:)
        integer, allocatable :: counts(:),offsets(:),cursor(:),members(:)
        integer, allocatable :: row_counts(:),row_offsets(:),row_cursor(:),row_members(:)
        real(dp) :: d2
        real :: dotdir
        integer :: n,ndim,nproj,k_used,cap_used,ir,i,j,k,m,p,q,f,nseen,pos
        n=size(features,2); ndim=size(features,1); nproj=size(proj_dirs,2)
        if( n<2 .or. ndim<1 .or. size(proj_ids)/=n ) THROW_HARD('invalid distributed gated graph inputs')
        if( any(rows<1).or.any(rows>n) ) THROW_HARD('gated graph row assignment outside feature table')
        if( any(proj_ids<1).or.any(proj_ids>nproj) ) THROW_HARD('projection id outside direction table')
        k_used=min(max(1,k_nn),n-1); cap_used=min(max(k_used,nang_nbrs),n-1)
        allocate(counts(nproj),offsets(nproj+1),cursor(nproj),source=0)
        do i=1,n; counts(proj_ids(i))=counts(proj_ids(i))+1; end do
        offsets(1)=1
        do p=1,nproj; offsets(p+1)=offsets(p)+counts(p); end do
        cursor=offsets(1:nproj); allocate(members(n),source=0)
        do i=1,n
            p=proj_ids(i); members(cursor(p))=i; cursor(p)=cursor(p)+1
        end do
        allocate(row_counts(nproj),row_offsets(nproj+1),row_cursor(nproj),source=0)
        do ir=1,size(rows); row_counts(proj_ids(rows(ir)))=row_counts(proj_ids(rows(ir)))+1; end do
        row_offsets(1)=1
        do p=1,nproj; row_offsets(p+1)=row_offsets(p)+row_counts(p); end do
        row_cursor=row_offsets(1:nproj); allocate(row_members(size(rows)),source=0)
        do ir=1,size(rows)
            p=proj_ids(rows(ir)); row_members(row_cursor(p))=ir; row_cursor(p)=row_cursor(p)+1
        end do
        allocate(nbrs(k_used,size(rows)),source=0)
        allocate(d2s(k_used,size(rows)),source=huge(1.))
        allocate(ncandidates(size(rows)),source=0)
        !$omp parallel do default(shared) private(q,dotdir,pos,ir,i,nseen,m,k,j,d2,f) schedule(dynamic)
        do p=1,nproj
            block
                integer :: proj_order(nproj)
                real :: angular_key(nproj)
                if( row_counts(p)==0 ) cycle
                proj_order=[(q,q=1,nproj)]
                do q=1,nproj
                    dotdir=max(-1.,min(1.,dot_product(proj_dirs(:,p),proj_dirs(:,q))))
                    angular_key(q)=1.-dotdir
                end do
                call hpsort(angular_key,proj_order)
                do pos=row_offsets(p),row_offsets(p+1)-1
                    ir=row_members(pos); i=rows(ir); nseen=0
                    direction_loop: do m=1,nproj
                        q=proj_order(m)
                        if( counts(q)==0 ) cycle
                        do k=offsets(q),offsets(q+1)-1
                            j=members(k)
                            if( j==i ) cycle
                            d2=0._dp
                            do f=1,ndim
                                d2=d2+real(features(f,i)-features(f,j),dp)**2
                            end do
                            call insert_neighbor(j,real(d2),nbrs(:,ir),d2s(:,ir))
                            nseen=nseen+1
                            if( nseen>=cap_used ) exit direction_loop
                        end do
                    end do direction_loop
                    ncandidates(ir)=nseen
                end do
            end block
        end do
        !$omp end parallel do
        if( any(nbrs(k_used,:)<1) ) THROW_HARD('too few angular candidates for requested flex k_nn')
        deallocate(counts,offsets,cursor,members,row_counts,row_offsets,row_cursor,row_members)
    end subroutine find_gated_euclidean_neighbors_rows

    subroutine build_gated_euclidean_graph_from_neighbors( n, nbrs, d2s, ncandidates, graph, &
        &bandwidth_mode, bandwidth_tune )
        integer, intent(in) :: n,nbrs(:,:),ncandidates(:)
        real, intent(in) :: d2s(:,:)
        type(diffmap_graph), intent(out) :: graph
        character(len=*), optional, intent(in) :: bandwidth_mode
        real, optional, intent(in) :: bandwidth_tune
        real, allocatable :: kth_d2(:)
        integer :: k_used
        if( n<2 .or. size(nbrs,2)/=n ) THROW_HARD('invalid gated neighbor table assembly')
        k_used=size(nbrs,1)
        if( k_used<1 .or. size(d2s,1)/=k_used .or. size(d2s,2)/=n ) &
            &THROW_HARD('invalid gated neighbor distance table assembly')
        if( size(ncandidates)/=n .or. any(nbrs<1).or.any(nbrs>n) ) THROW_HARD('incomplete gated neighbor table')
        allocate(kth_d2(n))
        kth_d2=d2s(k_used,:)
        call pack_scalar_knn_to_csr(n,k_used,nbrs,d2s,kth_d2,'euc_gated',graph,bandwidth_mode,bandwidth_tune)
        deallocate(kth_d2)
    end subroutine build_gated_euclidean_graph_from_neighbors

    subroutine build_orientation_knn_graph( params, spproj, pinds, k_nn, graph, bandwidth_mode, bandwidth_tune )
        type(parameters),      intent(in)    :: params
        type(sp_project),      intent(inout) :: spproj
        integer,               intent(in)    :: pinds(:), k_nn
        type(diffmap_graph),   intent(out)   :: graph
        character(len=*), optional, intent(in) :: bandwidth_mode
        real, optional, intent(in) :: bandwidth_tune
        integer, allocatable :: nbrs(:,:)
        real,    allocatable :: d2s(:,:), kth_d2(:)
        integer :: n, k_used
        n = size(pinds)
        if( n < 1 ) THROW_HARD('empty orientation graph')
        if( n == 1 )then
            call make_singleton_graph('ori', graph)
            return
        endif
        k_used = min(max(1, k_nn), n - 1)
        allocate(nbrs(k_used,n), source=0)
        allocate(d2s(k_used,n), kth_d2(n), source=0.)
        call find_orientation_neighbors(spproj, pinds, k_used, nbrs, d2s)
        kth_d2 = d2s(k_used,:)
        call pack_scalar_knn_to_csr(n, k_used, nbrs, d2s, kth_d2, 'ori', graph, bandwidth_mode, bandwidth_tune)
        deallocate(nbrs, d2s, kth_d2)
    end subroutine build_orientation_knn_graph

    subroutine find_euclidean_neighbors( pcavecs, k_used, nbrs, d2s )
        real,    intent(in)  :: pcavecs(:,:)
        integer, intent(in)  :: k_used
        integer, intent(out) :: nbrs(:,:)
        real,    intent(out) :: d2s(:,:)
        integer, allocatable :: local_nbrs(:,:,:)
        real,    allocatable :: local_d2s(:,:,:)
        real :: d2
        integer :: n, ndim, nthreads, tid, i, j, k, m
        n    = size(pcavecs, 2)
        ndim = size(pcavecs, 1)
        nthreads = 1
        !$ nthreads = omp_get_max_threads()
        allocate(local_nbrs(k_used,n,nthreads), source=0)
        allocate(local_d2s(k_used,n,nthreads))
        local_d2s = huge(1.)
        nbrs = 0
        d2s  = huge(1.)
        !$omp parallel default(shared) private(tid,i,j,k,d2)
        tid = 1
        !$ tid = omp_get_thread_num() + 1
        !$omp do schedule(dynamic)
        do i = 1,n - 1
            do j = i + 1,n
                d2 = 0.
                do k = 1,ndim
                    d2 = d2 + (pcavecs(k,i) - pcavecs(k,j))**2
                end do
                call insert_neighbor(j, d2, local_nbrs(:,i,tid), local_d2s(:,i,tid))
                call insert_neighbor(i, d2, local_nbrs(:,j,tid), local_d2s(:,j,tid))
            end do
        end do
        !$omp end do
        !$omp end parallel
        do tid = 1,nthreads
            do i = 1,n
                do m = 1,k_used
                    if( local_nbrs(m,i,tid) < 1 ) cycle
                    call insert_neighbor(local_nbrs(m,i,tid), local_d2s(m,i,tid), nbrs(:,i), d2s(:,i))
                end do
            end do
        end do
        deallocate(local_nbrs, local_d2s)
    end subroutine find_euclidean_neighbors

    subroutine find_orientation_neighbors( spproj, pinds, k_used, nbrs, d2s )
        type(sp_project), intent(inout) :: spproj
        integer,          intent(in)    :: pinds(:), k_used
        integer,          intent(out)   :: nbrs(:,:)
        real,             intent(out)   :: d2s(:,:)
        type(ori) :: oi, oj
        real :: dist
        integer :: n, i, j
        n = size(pinds)
        nbrs = 0
        d2s  = huge(1.)
        do i = 1,n
            call spproj%os_ptcl3D%get_ori(pinds(i), oi)
            do j = 1,n
                if( i == j ) cycle
                call spproj%os_ptcl3D%get_ori(pinds(j), oj)
                dist = oi .euldist. oj
                call insert_neighbor(j, dist**2, nbrs(:,i), d2s(:,i))
                call oj%kill
            end do
            call oi%kill
        end do
    end subroutine find_orientation_neighbors

    subroutine pack_scalar_knn_to_csr( n, k_used, nbrs, d2s, kth_d2, metric, graph, bandwidth_mode, bandwidth_tune )
        integer,              intent(in)  :: n, k_used, nbrs(:,:)
        real,                 intent(in)  :: d2s(:,:), kth_d2(:)
        character(len=*),     intent(in)  :: metric
        type(diffmap_graph),  intent(out) :: graph
        character(len=*), optional, intent(in) :: bandwidth_mode
        real, optional, intent(in) :: bandwidth_tune
        integer, allocatable :: counts(:), cursor(:)
        real :: eps, w
        integer :: i, m, j, p, nnz
        logical :: mutual
        logical :: ok
        character(len=STDLEN) :: bw_mode
        real :: bw_tune
        allocate(counts(n), cursor(n), source=0)
        do i = 1,n
            do m = 1,k_used
                j = nbrs(m,i)
                if( j < 1 ) cycle
                counts(i) = counts(i) + 1
                mutual = neighbor_contains(nbrs(:,j), i)
                if( .not. mutual ) counts(j) = counts(j) + 1
            end do
        end do
        nnz = sum(counts)
        graph%n        = n
        graph%nnz      = nnz
        graph%metric   = metric
        graph%k_nn     = k_used
        allocate(graph%rowptr(n+1), graph%colind(nnz), graph%w(nnz))
        graph%rowptr(1) = 1
        do i = 1,n
            graph%rowptr(i+1) = graph%rowptr(i) + counts(i)
        end do
        cursor = graph%rowptr(1:n)
        bw_mode = 'median'
        if( present(bandwidth_mode) ) bw_mode = lowercase(trim(bandwidth_mode))
        select case(trim(bw_mode))
        case('median')
            eps = median_positive(kth_d2)
            if( eps < DTINY ) eps = max(sum(kth_d2) / real(max(size(kth_d2),1)), 1.e-6)
        case('ferguson')
            bw_tune = 3.0
            if( present(bandwidth_tune) ) bw_tune = bandwidth_tune
            call estimate_ferguson_bandwidth(d2s, bw_tune, eps, ok)
            if( .not. ok )then
                eps = median_positive(kth_d2)
                if( eps < DTINY ) eps = max(sum(kth_d2) / real(max(size(kth_d2),1)), 1.e-6)
                write(logfhandle,'(A)') '>>> DIFFMAP bandwidth_mode=ferguson fallback -> median kth-NN scale'
            endif
        case DEFAULT
            THROW_HARD('invalid bandwidth_mode in diffusion graph; use median|ferguson')
        end select
        do i = 1,n
            do m = 1,k_used
                j = nbrs(m,i)
                if( j < 1 ) cycle
                w = exp(-max(d2s(m,i), 0.) / eps)
                p = cursor(i)
                graph%colind(p) = j
                graph%w(p)      = w
                cursor(i)       = p + 1
                mutual = neighbor_contains(nbrs(:,j), i)
                if( .not. mutual )then
                    p = cursor(j)
                    graph%colind(p) = i
                    graph%w(p)      = w
                    cursor(j)       = p + 1
                endif
            end do
        end do
        call graph%normalize()
        deallocate(counts, cursor)
    contains
        logical function neighbor_contains( row_nbrs, target ) result(found)
            integer, intent(in) :: row_nbrs(:), target
            integer :: q
            found = .false.
            do q = 1,size(row_nbrs)
                if( row_nbrs(q) /= target ) cycle
                found = .true.
                return
            end do
        end function neighbor_contains
    end subroutine pack_scalar_knn_to_csr

    subroutine estimate_ferguson_bandwidth( d2s, tune, eps, ok )
        real, intent(in) :: d2s(:,:), tune
        real, intent(out) :: eps
        logical, intent(out) :: ok
        integer, parameter :: NSCAN=121, NA=28, NB=33
        real(dp), allocatable :: d2(:), x(:), y(:), xz(:)
        real(dp) :: d2med, eps_min, eps_max, logeps_min, logeps_max
        real(dp) :: xmu, xsd, step, le, epsi, sumw, arg
        real(dp) :: a, b, t, s1, st, stt, sy, sty, det, cfit, dfit, pred, res, sse
        real(dp) :: best_sse, best_a, best_b, a0, b0, logeps_star, sigma2, tune_eff
        integer :: i, q, nvalid, ia, ib
        eps = 1.e-6
        ok = .false.
        allocate(d2(size(d2s)))
        nvalid = 0
        do i = 1,size(d2s,2)
            do q = 1,size(d2s,1)
                if( d2s(q,i) <= 0. .or. .not. ieee_is_finite(d2s(q,i)) ) cycle
                nvalid = nvalid + 1
                d2(nvalid) = real(d2s(q,i),dp)
            end do
        end do
        if( nvalid < 16 )then
            deallocate(d2)
            return
        endif
        d2med = real(median_positive(real(d2(1:nvalid))),dp)
        if( d2med <= real(DTINY,dp) )then
            deallocate(d2)
            return
        endif
        eps_min = max(d2med * 1.0d-3, 1.0d-12)
        eps_max = max(d2med * 1.0d+3, eps_min * 10.0d0)
        logeps_min = log(eps_min)
        logeps_max = log(eps_max)
        step = (logeps_max - logeps_min) / real(max(NSCAN-1,1),dp)
        allocate(x(NSCAN), y(NSCAN), xz(NSCAN))
        do i = 1,NSCAN
            le = logeps_min + real(i-1,dp) * step
            x(i) = le
            epsi = exp(le)
            sumw = 0.d0
            do q = 1,nvalid
                arg = -d2(q) / max(epsi, real(DTINY,dp))
                if( arg > -80.d0 ) sumw = sumw + exp(arg)
            end do
            y(i) = log(max(sumw, real(DTINY,dp)))
        end do
        xmu = sum(x) / real(NSCAN,dp)
        xsd = sqrt(sum((x - xmu)**2) / real(max(NSCAN-1,1),dp))
        if( xsd <= real(DTINY,dp) ) xsd = 1.d0
        xz = (x - xmu) / xsd
        best_sse = huge(1.d0)
        best_a = 0.d0
        best_b = 0.d0
        do ia = 1,NA
            a = 0.2d0 + real(ia-1,dp) * (3.8d0 / real(max(NA-1,1),dp))
            do ib = 1,NB
                b = -4.0d0 + real(ib-1,dp) * (8.0d0 / real(max(NB-1,1),dp))
                s1 = 0.d0; st = 0.d0; stt = 0.d0; sy = 0.d0; sty = 0.d0
                do i = 1,NSCAN
                    t = tanh(a * xz(i) + b)
                    s1 = s1 + 1.d0
                    st = st + t
                    stt = stt + t * t
                    sy = sy + y(i)
                    sty = sty + t * y(i)
                end do
                det = stt * s1 - st * st
                if( abs(det) <= 1.d-14 ) cycle
                cfit = (sty * s1 - sy * st) / det
                if( cfit <= 0.d0 ) cycle
                dfit = (sy * stt - sty * st) / det
                sse = 0.d0
                do i = 1,NSCAN
                    pred = dfit + cfit * tanh(a * xz(i) + b)
                    res = y(i) - pred
                    sse = sse + res * res
                end do
                if( sse < best_sse )then
                    best_sse = sse
                    best_a = a
                    best_b = b
                endif
            end do
        end do
        if( .not. ieee_is_finite(best_sse) .or. best_sse >= huge(1.d0) )then
            deallocate(d2, x, y, xz)
            return
        endif
        a0 = best_a / xsd
        b0 = best_b - best_a * xmu / xsd
        if( abs(a0) <= 1.d-10 )then
            deallocate(d2, x, y, xz)
            return
        endif
        logeps_star = -b0 / a0
        tune_eff = max(real(tune,dp), 0.d0)
        sigma2 = tune_eff * tune_eff * 2.d0 * exp(logeps_star)
        if( .not. ieee_is_finite(sigma2) .or. sigma2 <= real(DTINY,dp) )then
            deallocate(d2, x, y, xz)
            return
        endif
        eps = real(max(sigma2, 1.d-12), kind=kind(eps))
        ok = .true.
        deallocate(d2, x, y, xz)
    end subroutine estimate_ferguson_bandwidth

    subroutine make_singleton_graph( metric, graph )
        character(len=*),    intent(in)  :: metric
        type(diffmap_graph), intent(out) :: graph
        graph%n      = 1
        graph%nnz    = 1
        graph%k_nn   = 0
        graph%metric = metric
        allocate(graph%rowptr(2), graph%colind(1), graph%w(1), graph%wnorm(1))
        graph%rowptr = [1,2]
        graph%colind = [1]
        graph%w      = [1.]
        graph%wnorm  = [1.]
    end subroutine make_singleton_graph

    subroutine normalize_diffmap_graph( self )
        class(diffmap_graph), intent(inout) :: self
        real, allocatable :: deg(:)
        integer :: i, j, p
        if( self%n < 1 ) return
        allocate(deg(self%n), source=0.)
        call self%degree(deg, normalized=.false.)
        if( allocated(self%wnorm) ) deallocate(self%wnorm)
        allocate(self%wnorm(self%nnz), source=0.)
        do i = 1,self%n
            do p = self%rowptr(i), self%rowptr(i+1) - 1
                j = self%colind(p)
                self%wnorm(p) = self%w(p) / sqrt(max(deg(i), DTINY) * max(deg(j), DTINY))
            end do
        end do
        deallocate(deg)
    end subroutine normalize_diffmap_graph

    subroutine diffmap_graph_degree( self, deg, normalized )
        class(diffmap_graph), intent(in)  :: self
        real,                 intent(out) :: deg(:)
        logical, optional,    intent(in)  :: normalized
        logical :: use_norm
        integer :: i, p
        if( size(deg) /= self%n ) THROW_HARD('graph degree shape mismatch')
        use_norm = .false.
        if( present(normalized) ) use_norm = normalized
        deg = 0.
        do i = 1,self%n
            do p = self%rowptr(i), self%rowptr(i+1) - 1
                if( use_norm .and. allocated(self%wnorm) )then
                    deg(i) = deg(i) + self%wnorm(p)
                else
                    deg(i) = deg(i) + self%w(p)
                endif
            end do
        end do
    end subroutine diffmap_graph_degree

    subroutine graph_matvec( ctx, x, y )
        class(*), intent(in)  :: ctx
        real,     intent(in)  :: x(:)
        real,     intent(out) :: y(:)
        integer :: i, p
        select type(graph => ctx)
        type is (diffmap_graph)
            if( size(x) /= graph%n .or. size(y) /= graph%n ) THROW_HARD('sparse graph matvec shape mismatch')
            y = 0.
            do i = 1,graph%n
                do p = graph%rowptr(i), graph%rowptr(i+1) - 1
                    if( allocated(graph%wnorm) )then
                        y(i) = y(i) + graph%wnorm(p) * x(graph%colind(p))
                    else
                        y(i) = y(i) + graph%w(p) * x(graph%colind(p))
                    endif
                end do
            end do
        class default
            THROW_HARD('invalid graph context for matvec')
        end select
    end subroutine graph_matvec

    integer function graph_directed_edges( graph ) result(nedges)
        type(diffmap_graph), intent(in) :: graph
        nedges = graph%nnz
    end function graph_directed_edges

    subroutine kill_diffmap_graph( self )
        class(diffmap_graph), intent(inout) :: self
        if( allocated(self%rowptr) ) deallocate(self%rowptr)
        if( allocated(self%colind) ) deallocate(self%colind)
        if( allocated(self%w)      ) deallocate(self%w)
        if( allocated(self%wnorm)  ) deallocate(self%wnorm)
        self%n      = 0
        self%nnz    = 0
        self%k_nn   = 0
        self%metric = 'euc'
    end subroutine kill_diffmap_graph

    subroutine insert_neighbor( ind, dist, inds, dists )
        integer, intent(in)    :: ind
        real,    intent(in)    :: dist
        integer, intent(inout) :: inds(:)
        real,    intent(inout) :: dists(:)
        integer :: pos, k
        pos = 0
        do k = 1,size(dists)
            if( dist < dists(k) )then
                pos = k
                exit
            endif
        end do
        if( pos == 0 ) return
        do k = size(dists),pos+1,-1
            dists(k) = dists(k-1)
            inds(k)  = inds(k-1)
        end do
        dists(pos) = dist
        inds(pos)  = ind
    end subroutine insert_neighbor

    real function median_positive( vals ) result(med)
        real, intent(in) :: vals(:)
        real, allocatable :: work(:)
        integer :: i, n
        allocate(work(size(vals)))
        n = 0
        do i = 1,size(vals)
            if( vals(i) <= 0. .or. .not. ieee_is_finite(vals(i)) ) cycle
            n = n + 1
            work(n) = vals(i)
        end do
        if( n < 1 )then
            med = 1.e-6
        else
            call hpsort(work(1:n))
            if( mod(n,2) == 0 )then
                med = 0.5 * (work(n/2) + work(n/2 + 1))
            else
                med = work((n + 1) / 2)
            endif
        endif
        deallocate(work)
    end function median_positive

end module simple_diff_map_graphs
