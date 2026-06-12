!@descr: CSR graph construction helpers for diffusion-map class splitting
module simple_diff_map_graphs
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
!$ use omp_lib
use simple_builder,             only: builder
use simple_cmdline,             only: cmdline
use simple_default_clines,      only: set_cluster2D_defaults
use simple_image,               only: image
use simple_imgarr_utils,        only: dealloc_imgarr, write_imgarr
use simple_ori,                 only: ori
use simple_parameters,          only: parameters
use simple_srch_sort_loc,       only: hpsort
use simple_sp_project,          only: sp_project
implicit none
#include "simple_local_flags.inc"

private
integer, parameter :: GRAPH_GAUGE_CLUSTER2D_MAXITS = 5

public :: diffmap_graph
public :: diffmap_gauge
public :: build_cls_split_graph
public :: build_euclidean_knn_graph
public :: build_orientation_knn_graph
public :: run_single_class_cluster2d_gauge
public :: graph_matvec
public :: estimate_graph_shift_scale
public :: graph_directed_edges

type :: diffmap_graph
    integer :: n   = 0
    integer :: nnz = 0
    integer :: k_nn = 0
    character(len=16) :: metric   = 'euclidean'
    character(len=16) :: steering = 'none'
    integer, allocatable :: rowptr(:)
    integer, allocatable :: colind(:)
    real,    allocatable :: w(:)
    real,    allocatable :: wnorm(:)
    real,    allocatable :: theta(:)
    real,    allocatable :: shift_x(:)
    real,    allocatable :: shift_y(:)
contains
    procedure :: kill      => kill_diffmap_graph
    procedure :: normalize => normalize_diffmap_graph
    procedure :: degree    => diffmap_graph_degree
    procedure :: has_theta => diffmap_graph_has_theta
    procedure :: has_shift => diffmap_graph_has_shift
end type diffmap_graph

type :: diffmap_gauge
    real, allocatable :: angle(:)
    real, allocatable :: shift_x(:)
    real, allocatable :: shift_y(:)
    real, allocatable :: img_dist(:)
contains
    procedure :: kill => kill_diffmap_gauge
end type diffmap_gauge

contains

    subroutine build_cls_split_graph( params, spproj, pinds, pcavecs, imgs, graph )
        type(parameters),           intent(in)    :: params
        type(sp_project),           intent(inout) :: spproj
        integer,                    intent(in)    :: pinds(:)
        real, optional,             intent(in)    :: pcavecs(:,:)
        type(image), optional,      intent(inout) :: imgs(:)
        type(diffmap_graph),        intent(out)   :: graph
        type(diffmap_gauge) :: gauge
        character(len=16) :: metric, steering
        logical :: need_gauge
        metric = lowercase(trim(params%graph))
        select case(trim(metric))
            case('', 'auto')
                metric = 'euclidean'
                if( trim(params%oritype) == 'ptcl3D' ) metric = 'orientation'
            case('euc', 'euclidean')
                metric = 'euclidean'
            case('ori', 'orientation')
                metric = 'orientation'
            case DEFAULT
                THROW_HARD('graph must be euc or ori in build_cls_split_graph')
        end select
        steering = lowercase(trim(params%steering))
        if( trim(steering) == '' .or. trim(steering) == 'auto' ) steering = 'none'
        select case(trim(steering))
            case('none', 'so2', 'se2')
            case DEFAULT
                THROW_HARD('unsupported graph steering in build_cls_split_graph')
        end select
        if( trim(metric) == 'euclidean' .and. trim(steering) /= 'none' )then
            THROW_HARD('graph=euc supports steering=none only in class splitting')
        endif
        need_gauge = trim(metric) == 'orientation' .and. trim(steering) /= 'none'
        if( need_gauge )then
            if( .not. present(imgs) ) THROW_HARD('steered orientation graph requires images')
            call run_single_class_cluster2d_gauge(params, pinds, imgs, gauge, steering)
        endif
        select case(trim(metric))
            case('orientation')
                if( need_gauge )then
                    call build_orientation_knn_graph(params, spproj, pinds, max(2, params%k_nn), steering, graph, gauge)
                else
                    call build_orientation_knn_graph(params, spproj, pinds, max(2, params%k_nn), steering, graph)
                endif
            case DEFAULT
                if( .not. present(pcavecs) ) THROW_HARD('Euclidean graph requires pcavecs')
                call build_euclidean_knn_graph(pcavecs, max(2, params%k_nn), 'none', graph)
        end select
        if( need_gauge ) call gauge%kill
    end subroutine build_cls_split_graph

    subroutine build_euclidean_knn_graph( pcavecs, k_nn, steering, graph, gauge )
        real,                 intent(in)  :: pcavecs(:,:)
        integer,              intent(in)  :: k_nn
        character(len=*),     intent(in)  :: steering
        type(diffmap_graph),  intent(out) :: graph
        type(diffmap_gauge), optional, intent(in) :: gauge
        integer, allocatable :: nbrs(:,:)
        real,    allocatable :: d2s(:,:), kth_d2(:)
        integer :: n, k_used
        n = size(pcavecs, 2)
        if( n < 1 ) THROW_HARD('empty Euclidean graph')
        if( n == 1 )then
            call make_singleton_graph('euclidean', steering, graph)
            return
        endif
        k_used = min(max(1, k_nn), n - 1)
        allocate(nbrs(k_used,n), source=0)
        allocate(d2s(k_used,n), kth_d2(n), source=0.)
        call find_euclidean_neighbors(pcavecs, k_used, nbrs, d2s)
        kth_d2 = d2s(k_used,:)
        call pack_knn_to_csr(n, k_used, nbrs, d2s, kth_d2, 'euclidean', steering, graph, gauge)
        deallocate(nbrs, d2s, kth_d2)
    end subroutine build_euclidean_knn_graph

    subroutine build_orientation_knn_graph( params, spproj, pinds, k_nn, steering, graph, gauge )
        type(parameters),      intent(in)    :: params
        type(sp_project),      intent(inout) :: spproj
        integer,               intent(in)    :: pinds(:), k_nn
        character(len=*),      intent(in)    :: steering
        type(diffmap_graph),   intent(out)   :: graph
        type(diffmap_gauge), optional, intent(in) :: gauge
        integer, allocatable :: nbrs(:,:)
        real,    allocatable :: d2s(:,:), kth_d2(:)
        integer :: n, k_used
        n = size(pinds)
        if( n < 1 ) THROW_HARD('empty orientation graph')
        if( n == 1 )then
            call make_singleton_graph('orientation', steering, graph)
            return
        endif
        k_used = min(max(1, k_nn), n - 1)
        allocate(nbrs(k_used,n), source=0)
        allocate(d2s(k_used,n), kth_d2(n), source=0.)
        call find_orientation_neighbors(spproj, pinds, k_used, nbrs, d2s)
        kth_d2 = d2s(k_used,:)
        call pack_knn_to_csr(n, k_used, nbrs, d2s, kth_d2, 'orientation', steering, graph, gauge)
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

    subroutine pack_knn_to_csr( n, k_used, nbrs, d2s, kth_d2, metric, steering, graph, gauge )
        integer,              intent(in)  :: n, k_used, nbrs(:,:)
        real,                 intent(in)  :: d2s(:,:), kth_d2(:)
        character(len=*),     intent(in)  :: metric, steering
        type(diffmap_graph),  intent(out) :: graph
        type(diffmap_gauge), optional, intent(in) :: gauge
        integer, allocatable :: rows(:), cols(:)
        real,    allocatable :: weights(:), theta(:), sx(:), sy(:)
        real :: eps, w, th, dx, dy
        integer :: max_edges, pos, i, m, j
        logical :: use_theta, use_shift
        use_theta = present(gauge) .and. (trim(steering) == 'so2' .or. trim(steering) == 'se2')
        use_shift = present(gauge) .and. trim(steering) == 'se2'
        if( .not. use_theta )then
            call pack_scalar_knn_to_csr(n, k_used, nbrs, d2s, kth_d2, metric, steering, graph)
            return
        endif
        max_edges = 2 * n * k_used
        allocate(rows(max_edges), cols(max_edges), source=0)
        allocate(weights(max_edges), source=0.)
        if( use_theta ) allocate(theta(max_edges), source=0.)
        if( use_shift ) allocate(sx(max_edges), sy(max_edges), source=0.)
        eps = median_positive(kth_d2)
        if( eps < DTINY ) eps = max(sum(kth_d2) / real(max(size(kth_d2),1)), 1.e-6)
        pos = 0
        do i = 1,n
            do m = 1,k_used
                j = nbrs(m,i)
                if( j < 1 ) cycle
                w = exp(-max(d2s(m,i), 0.) / eps)
                call add_edge(i, j, w, 1.)
                call add_edge(j, i, w, -1.)
            end do
        end do
        if( use_theta .and. use_shift )then
            call coalesce_coo_to_csr(n, rows(:pos), cols(:pos), weights(:pos), metric, steering, graph, &
                                     theta(:pos), sx(:pos), sy(:pos))
        else if( use_theta )then
            call coalesce_coo_to_csr(n, rows(:pos), cols(:pos), weights(:pos), metric, steering, graph, theta(:pos))
        else
            call coalesce_coo_to_csr(n, rows(:pos), cols(:pos), weights(:pos), metric, steering, graph)
        endif
        graph%k_nn = k_used
        call graph%normalize()
        deallocate(rows, cols, weights)
        if( allocated(theta) ) deallocate(theta)
        if( allocated(sx)    ) deallocate(sx)
        if( allocated(sy)    ) deallocate(sy)
    contains
        subroutine add_edge( ii, jj, ww, sign )
            integer, intent(in) :: ii, jj
            real,    intent(in) :: ww, sign
            pos = pos + 1
            rows(pos) = ii
            cols(pos) = jj
            weights(pos) = ww
            if( use_theta )then
                th = deg2rad(wrap_angle_pm180(gauge%angle(jj) - gauge%angle(ii)))
                theta(pos) = sign * th
            endif
            if( use_shift )then
                dx = gauge%shift_x(jj) - gauge%shift_x(ii)
                dy = gauge%shift_y(jj) - gauge%shift_y(ii)
                sx(pos) = sign * dx
                sy(pos) = sign * dy
            endif
        end subroutine add_edge
    end subroutine pack_knn_to_csr

    subroutine pack_scalar_knn_to_csr( n, k_used, nbrs, d2s, kth_d2, metric, steering, graph )
        integer,              intent(in)  :: n, k_used, nbrs(:,:)
        real,                 intent(in)  :: d2s(:,:), kth_d2(:)
        character(len=*),     intent(in)  :: metric, steering
        type(diffmap_graph),  intent(out) :: graph
        integer, allocatable :: counts(:), cursor(:)
        real :: eps, w
        integer :: i, m, j, p, nnz
        logical :: mutual
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
        graph%steering = steering
        graph%k_nn     = k_used
        allocate(graph%rowptr(n+1), graph%colind(nnz), graph%w(nnz))
        graph%rowptr(1) = 1
        do i = 1,n
            graph%rowptr(i+1) = graph%rowptr(i) + counts(i)
        end do
        cursor = graph%rowptr(1:n)
        eps = median_positive(kth_d2)
        if( eps < DTINY ) eps = max(sum(kth_d2) / real(max(size(kth_d2),1)), 1.e-6)
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

    subroutine coalesce_coo_to_csr( n, rows, cols, weights, metric, steering, graph, theta, sx, sy )
        integer,              intent(in)  :: n, rows(:), cols(:)
        real,                 intent(in)  :: weights(:)
        character(len=*),     intent(in)  :: metric, steering
        type(diffmap_graph),  intent(out) :: graph
        real, optional,       intent(in)  :: theta(:), sx(:), sy(:)
        integer, allocatable :: out_rows(:), out_cols(:), counts(:), cursor(:)
        real,    allocatable :: out_w(:), out_theta(:), out_sx(:), out_sy(:)
        integer :: e, f, nnz, p, i
        logical :: found, with_theta, with_shift
        if( size(rows) /= size(cols) .or. size(rows) /= size(weights) ) THROW_HARD('COO graph size mismatch')
        with_theta = present(theta)
        with_shift = present(sx) .and. present(sy)
        allocate(out_rows(size(rows)), out_cols(size(cols)), source=0)
        allocate(out_w(size(weights)), source=0.)
        if( with_theta ) allocate(out_theta(size(rows)), source=0.)
        if( with_shift ) allocate(out_sx(size(rows)), out_sy(size(rows)), source=0.)
        nnz = 0
        do e = 1,size(rows)
            if( rows(e) < 1 .or. rows(e) > n .or. cols(e) < 1 .or. cols(e) > n ) THROW_HARD('COO graph index out of range')
            if( weights(e) <= DTINY .or. .not. ieee_is_finite(weights(e)) ) cycle
            found = .false.
            do f = 1,nnz
                if( out_rows(f) /= rows(e) .or. out_cols(f) /= cols(e) ) cycle
                found = .true.
                if( weights(e) > out_w(f) )then
                    out_w(f) = weights(e)
                    if( with_theta ) out_theta(f) = theta(e)
                    if( with_shift )then
                        out_sx(f) = sx(e)
                        out_sy(f) = sy(e)
                    endif
                endif
                exit
            end do
            if( .not. found )then
                nnz = nnz + 1
                out_rows(nnz) = rows(e)
                out_cols(nnz) = cols(e)
                out_w(nnz)    = weights(e)
                if( with_theta ) out_theta(nnz) = theta(e)
                if( with_shift )then
                    out_sx(nnz) = sx(e)
                    out_sy(nnz) = sy(e)
                endif
            endif
        end do
        graph%n        = n
        graph%nnz      = nnz
        graph%metric   = metric
        graph%steering = steering
        allocate(graph%rowptr(n+1), graph%colind(nnz), graph%w(nnz))
        if( with_theta ) allocate(graph%theta(nnz), source=0.)
        if( with_shift ) allocate(graph%shift_x(nnz), graph%shift_y(nnz), source=0.)
        allocate(counts(n), cursor(n), source=0)
        do e = 1,nnz
            counts(out_rows(e)) = counts(out_rows(e)) + 1
        end do
        graph%rowptr(1) = 1
        do i = 1,n
            graph%rowptr(i+1) = graph%rowptr(i) + counts(i)
        end do
        cursor = graph%rowptr(1:n)
        do e = 1,nnz
            p = cursor(out_rows(e))
            graph%colind(p) = out_cols(e)
            graph%w(p)      = out_w(e)
            if( with_theta ) graph%theta(p) = out_theta(e)
            if( with_shift )then
                graph%shift_x(p) = out_sx(e)
                graph%shift_y(p) = out_sy(e)
            endif
            cursor(out_rows(e)) = cursor(out_rows(e)) + 1
        end do
        deallocate(out_rows, out_cols, out_w, counts, cursor)
        if( allocated(out_theta) ) deallocate(out_theta)
        if( allocated(out_sx)    ) deallocate(out_sx)
        if( allocated(out_sy)    ) deallocate(out_sy)
    end subroutine coalesce_coo_to_csr

    subroutine make_singleton_graph( metric, steering, graph )
        character(len=*),    intent(in)  :: metric, steering
        type(diffmap_graph), intent(out) :: graph
        graph%n        = 1
        graph%nnz      = 1
        graph%k_nn     = 0
        graph%metric   = metric
        graph%steering = steering
        allocate(graph%rowptr(2), graph%colind(1), graph%w(1), graph%wnorm(1))
        graph%rowptr = [1,2]
        graph%colind = [1]
        graph%w      = [1.]
        graph%wnorm  = [1.]
        if( trim(steering) == 'so2' .or. trim(steering) == 'se2' ) allocate(graph%theta(1), source=0.)
        if( trim(steering) == 'se2' ) allocate(graph%shift_x(1), graph%shift_y(1), source=0.)
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

    logical function diffmap_graph_has_theta( self ) result(l_has)
        class(diffmap_graph), intent(in) :: self
        l_has = allocated(self%theta)
    end function diffmap_graph_has_theta

    logical function diffmap_graph_has_shift( self ) result(l_has)
        class(diffmap_graph), intent(in) :: self
        l_has = allocated(self%shift_x) .and. allocated(self%shift_y)
    end function diffmap_graph_has_shift

    subroutine kill_diffmap_graph( self )
        class(diffmap_graph), intent(inout) :: self
        if( allocated(self%rowptr)  ) deallocate(self%rowptr)
        if( allocated(self%colind)  ) deallocate(self%colind)
        if( allocated(self%w)       ) deallocate(self%w)
        if( allocated(self%wnorm)   ) deallocate(self%wnorm)
        if( allocated(self%theta)   ) deallocate(self%theta)
        if( allocated(self%shift_x) ) deallocate(self%shift_x)
        if( allocated(self%shift_y) ) deallocate(self%shift_y)
        self%n = 0
        self%nnz = 0
        self%k_nn = 0
        self%metric = 'euclidean'
        self%steering = 'none'
    end subroutine kill_diffmap_graph

    subroutine kill_diffmap_gauge( self )
        class(diffmap_gauge), intent(inout) :: self
        if( allocated(self%angle)    ) deallocate(self%angle)
        if( allocated(self%shift_x)  ) deallocate(self%shift_x)
        if( allocated(self%shift_y)  ) deallocate(self%shift_y)
        if( allocated(self%img_dist) ) deallocate(self%img_dist)
    end subroutine kill_diffmap_gauge

    subroutine run_single_class_cluster2d_gauge( params, pinds, imgs, gauge, steering )
        use simple_strategy2D_matcher, only: cluster2D_exec
        type(parameters),    intent(in)    :: params
        integer,             intent(in)    :: pinds(:)
        type(image),         intent(inout) :: imgs(:)
        type(diffmap_gauge), intent(out)   :: gauge
        character(len=*), optional, intent(in) :: steering
        type(parameters) :: c2d_params
        type(builder)    :: c2d_build
        type(cmdline)    :: c2d_cline
        type(sp_project) :: local_proj
        type(oris)       :: local_os
        type(ori)        :: o
        type(ctfparams)  :: ctfvars
        type(image), allocatable :: ref_imgs(:)
        type(string) :: cwd_before, cwd_now, workdir, stk_fname, projfile, refs_fname, refs_even_fname, refs_odd_fname, dot
        integer :: n, i, iter, ldim(3), istat, next_dir, nlocal_maxits, first_label
        logical :: converged
        real :: smpd, mskdiam_local, local_trs
        character(len=16) :: local_steering
        n = size(imgs)
        if( n /= size(pinds) ) THROW_HARD('gauge particle/image count mismatch')
        allocate(gauge%angle(n), gauge%shift_x(n), gauge%shift_y(n), gauge%img_dist(n), source=0.)
        if( n < 1 ) return
        ldim = imgs(1)%get_ldim()
        smpd = imgs(1)%get_smpd()
        call simple_getcwd(cwd_before)
        dot      = '.'
        next_dir = find_next_int_dir_prefix(dot)
        first_label = merge(pinds(1), 1, size(pinds) > 0)
        workdir  = int2str(next_dir)//'_diffmap_gauge_cluster2d_'//int2str(first_label)
        call simple_mkdir(workdir, verbose=.false.)
        call simple_chdir(workdir, status=istat)
        if( istat /= 0 ) THROW_HARD('failed to enter diffusion-map gauge cluster2D directory')
        call simple_getcwd(cwd_now)
        CWD_GLOB = cwd_now%to_char()

        stk_fname = 'particles.mrcs'
        projfile  = 'gauge_cluster2d.simple'
        refs_fname      = 'start2Drefs.mrc'
        refs_even_fname = 'start2Drefs_even.mrc'
        refs_odd_fname  = 'start2Drefs_odd.mrc'
        call write_imgarr(imgs, stk_fname)
        allocate(ref_imgs(1))
        call ref_imgs(1)%new(ldim, smpd)
        call ref_imgs(1)%zero_and_unflag_ft
        do i = 1,n
            call ref_imgs(1)%add(imgs(i))
        end do
        call ref_imgs(1)%div(real(max(n,1)))
        call write_imgarr(ref_imgs, refs_fname)
        call write_imgarr(ref_imgs, refs_even_fname)
        call write_imgarr(ref_imgs, refs_odd_fname)
        call local_proj%update_projinfo(projfile)
        call local_os%new(n, is_ptcl=.true.)
        do i = 1,n
            call o%new(is_ptcl=.true.)
            call o%set('state', 1.)
            call o%set('w',     1.)
            call o%set('x',     0.)
            call o%set('y',     0.)
            call o%set('class', 1.)
            call o%e3set(0.)
            call local_os%set_ori(i, o)
            call o%kill
        end do
        ctfvars%ctfflag = CTFFLAG_NO
        ctfvars%smpd    = smpd
        ctfvars%kv      = max(params%kv, 300.)
        ctfvars%cs      = max(params%cs, 2.7)
        ctfvars%fraca   = max(params%fraca, 0.1)
        call local_proj%add_single_stk(stk_fname, ctfvars, local_os)
        call local_proj%write(projfile)

        mskdiam_local = params%mskdiam
        if( mskdiam_local <= 0. ) mskdiam_local = real(ldim(1)) * smpd
        local_steering = lowercase(trim(params%steering))
        if( present(steering) ) local_steering = lowercase(trim(steering))
        local_trs = 0.
        if( trim(local_steering) == 'se2' ) local_trs = max(params%trs, 3.)
        nlocal_maxits = GRAPH_GAUGE_CLUSTER2D_MAXITS
        call c2d_cline%set('prg',           'cluster2D')
        call c2d_cline%set('projfile',      projfile%to_char())
        call c2d_cline%set('refs',          refs_fname%to_char())
        call c2d_cline%set('refs_even',     refs_even_fname%to_char())
        call c2d_cline%set('refs_odd',      refs_odd_fname%to_char())
        call c2d_cline%set('outfile',       ALGN_FBODY//int2str(1)//METADATA_EXT)
        call c2d_cline%set('ncls',          1)
        call c2d_cline%set('mskdiam',       mskdiam_local)
        call c2d_cline%set('smpd',          smpd)
        call c2d_cline%set('box',           ldim(1))
        call c2d_cline%set('ctf',           'no')
        call c2d_cline%set('objfun',        'cc')
        call c2d_cline%set('refine',        'inpl')
        call c2d_cline%set('ml_reg',        'no')
        call c2d_cline%set('maxits',        nlocal_maxits)
        call c2d_cline%set('minits',        nlocal_maxits)
        call c2d_cline%set('converge',      'no')
        call c2d_cline%set('startit',       1)
        call c2d_cline%set('nthr',          max(1, params%nthr))
        call c2d_cline%set('trs',           local_trs)
        call c2d_cline%set('mkdir',         'no')
        call c2d_cline%set('restore_cavgs', 'yes')
        call set_cluster2D_defaults(c2d_cline)
        call c2d_build%init_params_and_build_strategy2D_tbox(c2d_cline, c2d_params, wthreads=.true.)
        if( c2d_build%spproj_field%get_nevenodd() == 0 )then
            call c2d_build%spproj_field%partition_eo
            call c2d_build%spproj%write_segment_inside(c2d_params%oritype, c2d_params%projfile)
        endif
        c2d_params%which_iter = c2d_params%startit - 1
        if( c2d_cline%defined('extr_iter') )then
            c2d_params%extr_iter = c2d_params%extr_iter - 1
        else
            c2d_params%extr_iter = c2d_params%startit - 1
        endif
        do iter = 1,nlocal_maxits
            c2d_params%which_iter = c2d_params%which_iter + 1
            c2d_params%extr_iter  = c2d_params%extr_iter  + 1
            call c2d_cline%set('which_iter', c2d_params%which_iter)
            call c2d_cline%set('extr_iter',  c2d_params%extr_iter)
            call c2d_cline%set('outfile', ALGN_FBODY//int2str_pad(c2d_params%part,c2d_params%numlen)//METADATA_EXT)
            c2d_params%outfile = ALGN_FBODY//int2str_pad(c2d_params%part,c2d_params%numlen)//METADATA_EXT
            call cluster2D_exec(c2d_params, c2d_build, c2d_cline, c2d_params%which_iter, converged)
        end do
        do i = 1,n
            gauge%angle(i) = c2d_build%spproj%os_ptcl2D%e3get(i)
            if( c2d_build%spproj%os_ptcl2D%isthere('x') ) gauge%shift_x(i) = c2d_build%spproj%os_ptcl2D%get(i, 'x')
            if( c2d_build%spproj%os_ptcl2D%isthere('y') ) gauge%shift_y(i) = c2d_build%spproj%os_ptcl2D%get(i, 'y')
        end do
        call local_proj%kill
        call local_os%kill
        if( allocated(ref_imgs) ) call dealloc_imgarr(ref_imgs)
        call c2d_build%kill_general_tbox()
        call c2d_build%kill_strategy2D_tbox()
        call c2d_cline%kill()
        call simple_chdir(cwd_before, status=istat)
        if( istat /= 0 ) THROW_HARD('failed to restore cwd after diffusion-map gauge cluster2D')
        call simple_getcwd(cwd_now)
        CWD_GLOB = cwd_now%to_char()
    end subroutine run_single_class_cluster2d_gauge

    real function estimate_graph_shift_scale( graph ) result(scale)
        type(diffmap_graph), intent(in) :: graph
        real :: acc, r
        integer :: p, n
        scale = 1.
        if( .not. graph%has_shift() ) return
        acc = 0.
        n   = 0
        do p = 1,graph%nnz
            r = sqrt(graph%shift_x(p)**2 + graph%shift_y(p)**2)
            if( r <= DTINY .or. .not. ieee_is_finite(r) ) cycle
            acc = acc + r
            n = n + 1
        end do
        if( n > 0 ) scale = max(acc / real(n), 1.e-3)
    end function estimate_graph_shift_scale

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

    real function wrap_angle_pm180( angle_deg ) result(wrapped)
        real, intent(in) :: angle_deg
        wrapped = modulo(angle_deg + 180., 360.) - 180.
    end function wrap_angle_pm180

end module simple_diff_map_graphs
