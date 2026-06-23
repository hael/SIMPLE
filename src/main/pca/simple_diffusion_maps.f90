!@descr: sparse diffusion-map embeddings for class-split graphs
module simple_diffusion_maps
use simple_core_module_api
use simple_diff_map_graphs,  only: diffmap_graph, build_euclidean_knn_graph, graph_matvec, estimate_graph_shift_scale
use simple_linalg,           only: sparse_eigh
use simple_stat,             only: pearsn
implicit none
#include "simple_local_flags.inc"

private
public :: diffusion_map_embedder
public :: steerable_diffusion_map_embedder
public :: embed_graph
public :: embed_so2_graph
public :: embed_se2_graph
public :: normalize_coords
public :: auto_ndiff_from_eigengap

integer, parameter :: SE2_NTRANS_MODES = 2

type :: diffusion_map_embedder
    integer :: ndiff = 4
    integer :: k_nn  = 10
contains
    procedure :: set_params
    procedure :: embed
end type diffusion_map_embedder

type :: steerable_diffusion_map_embedder
    integer :: ndiff  = 4
    integer :: k_nn   = 10
    integer :: nmodes = 4
contains
    procedure :: set_params => steerable_set_params
    procedure :: embed      => steerable_embed
end type steerable_diffusion_map_embedder

type :: connection_context
    type(diffmap_graph), pointer :: graph => null()
    integer :: mode = 0
    real    :: kx = 0.
    real    :: ky = 0.
    real    :: trans_factor = 0.
    logical :: use_shift = .false.
end type connection_context

contains

    subroutine set_params(self, ndiff, k_nn)
        class(diffusion_map_embedder), intent(inout) :: self
        integer,                       intent(in)    :: ndiff, k_nn
        self%ndiff = max(0, ndiff)
        self%k_nn  = max(2, k_nn)
    end subroutine set_params

    subroutine embed(self, pcavecs, coords, eigvals)
        class(diffusion_map_embedder), intent(inout) :: self
        real,                          intent(in)    :: pcavecs(:,:)
        real, allocatable,             intent(out)   :: coords(:,:)
        real, allocatable, optional,   intent(out)   :: eigvals(:)
        type(diffmap_graph), target :: graph
        real, allocatable :: eigs(:)
        integer :: nptcls
        nptcls = size(pcavecs,2)
        if( nptcls < 3 )then
            allocate(coords(1,nptcls), source=0.)
            if( present(eigvals) ) allocate(eigvals(1), source=0.)
            return
        endif
        call build_euclidean_knn_graph(pcavecs, min(max(2,self%k_nn), nptcls-1), 'none', graph)
        call embed_graph(graph, self%ndiff, coords, eigs)
        if( present(eigvals) ) allocate(eigvals(size(eigs)), source=eigs)
        call graph%kill()
        if( allocated(eigs) ) deallocate(eigs)
    end subroutine embed

    subroutine steerable_set_params(self, ndiff, k_nn, nmodes)
        class(steerable_diffusion_map_embedder), intent(inout) :: self
        integer,                                 intent(in)    :: ndiff, k_nn, nmodes
        self%ndiff  = max(0, ndiff)
        self%k_nn   = max(2, k_nn)
        self%nmodes = max(0, nmodes)
    end subroutine steerable_set_params

    subroutine steerable_embed(self, graph, coords, eigvals, shift_scale)
        class(steerable_diffusion_map_embedder), intent(inout) :: self
        type(diffmap_graph), target,             intent(in)    :: graph
        real, allocatable,                       intent(out)   :: coords(:,:)
        real, allocatable, optional,             intent(out)   :: eigvals(:)
        real, optional,                          intent(in)    :: shift_scale
        real, allocatable :: eigs(:)
        real :: se2_shift_scale
        if( graph%n < 1 ) THROW_HARD('steerable diffusion map embedding requires a prebuilt graph')
        select case(lowercase(trim(graph%steering)))
            case('so2')
                if( .not. graph%has_theta() ) THROW_HARD('SO2 steerable embedding requires theta payloads')
                if( graph%n < 3 )then
                    allocate(coords(1,graph%n), source=0.)
                    if( present(eigvals) ) allocate(eigvals(1), source=0.)
                    return
                endif
                call embed_so2_graph(graph, self%ndiff, self%nmodes, coords, eigs)
            case('se2')
                if( .not. graph%has_theta() .or. .not. graph%has_shift() )then
                    THROW_HARD('SE2 steerable embedding requires theta/shift payloads')
                endif
                if( graph%n < 3 )then
                    allocate(coords(1,graph%n), source=0.)
                    if( present(eigvals) ) allocate(eigvals(1), source=0.)
                    return
                endif
                if( present(shift_scale) )then
                    se2_shift_scale = shift_scale
                else
                    se2_shift_scale = estimate_graph_shift_scale(graph)
                endif
                call embed_se2_graph(graph, self%ndiff, self%nmodes, se2_shift_scale, coords, eigs)
            case DEFAULT
                THROW_HARD('steerable diffusion map embedding requires graph%steering=so2|se2')
        end select
        if( present(eigvals) ) allocate(eigvals(size(eigs)), source=eigs)
        if( allocated(eigs) ) deallocate(eigs)
    end subroutine steerable_embed

    subroutine embed_graph( graph, ndiff_req, coords, eigvals )
        type(diffmap_graph), target, intent(in)  :: graph
        integer,                   intent(in)    :: ndiff_req
        real, allocatable,         intent(out)   :: coords(:,:), eigvals(:)
        real, allocatable :: evals(:), evecs(:,:), diff_evals(:)
        integer :: n, ndiff_scan, ndiff_used, nev, eig_info, max_basis, i, k, j
        n = graph%n
        if( n < 3 )then
            allocate(coords(1,n), eigvals(1), source=0.)
            return
        endif
        if( ndiff_req <= 0 )then
            ndiff_scan = min(16, max(1, n - 2))
        else
            ndiff_scan = min(max(1, ndiff_req), max(1, n - 2))
        endif
        nev = ndiff_scan + 1
        allocate(evals(nev), evecs(n,nev))
        max_basis = min(n, max(160, 8 * nev + 80))
        call sparse_eigh(graph_matvec, graph, n, nev, evals, evecs, tol=1.e-5, max_basis=max_basis, info=eig_info)
        allocate(diff_evals(ndiff_scan), source=0.)
        do k = 1,ndiff_scan
            diff_evals(k) = evals(nev - k)
        end do
        if( ndiff_req <= 0 )then
            ndiff_used = auto_ndiff_from_eigengap(diff_evals)
        else
            ndiff_used = ndiff_scan
        endif
        allocate(coords(ndiff_used,n), eigvals(ndiff_scan), source=0.)
        eigvals = diff_evals
        do k = 1,ndiff_used
            j = nev - k
            do i = 1,n
                coords(k,i) = evals(j) * evecs(i,j)
            end do
        end do
        call normalize_coords(coords)
        deallocate(evals, evecs, diff_evals)
    end subroutine embed_graph

    subroutine embed_so2_graph( graph, ndiff_req, nmodes_req, coords, eigvals )
        type(diffmap_graph), target, intent(in)  :: graph
        integer,                   intent(in)    :: ndiff_req, nmodes_req
        real, allocatable,         intent(out)   :: coords(:,:), eigvals(:)
        real, allocatable :: cand_coords(:,:), cand_eigvals(:), out_eigvals(:)
        integer :: n, ndiff_scan, nmodes, max_cand, ncand
        n = graph%n
        if( n < 3 )then
            allocate(coords(1,n), eigvals(1), source=0.)
            return
        endif
        if( .not. graph%has_theta() ) THROW_HARD('SO2 graph embedding requires theta payloads')
        if( ndiff_req <= 0 )then
            ndiff_scan = min(16, max(1, n - 2))
        else
            ndiff_scan = min(max(1, ndiff_req), max(1, n - 2))
        endif
        nmodes = min(max(0, nmodes_req), max(0, n - 1))
        max_cand = max(1, ndiff_scan * (nmodes + 1))
        allocate(cand_coords(n,max_cand), cand_eigvals(max_cand), source=0.)
        call collect_so2_candidates(graph, nmodes, ndiff_scan, cand_coords, cand_eigvals, ncand)
        call select_output_candidates(cand_coords, cand_eigvals, ncand, ndiff_scan, coords, out_eigvals)
        call normalize_coords(coords)
        allocate(eigvals(size(out_eigvals)), source=out_eigvals)
        deallocate(cand_coords, cand_eigvals, out_eigvals)
    end subroutine embed_so2_graph

    subroutine embed_se2_graph( graph, ndiff_req, nmodes_req, shift_scale, coords, eigvals )
        type(diffmap_graph), target, intent(in)  :: graph
        integer,                   intent(in)    :: ndiff_req, nmodes_req
        real,                      intent(in)    :: shift_scale
        real, allocatable,         intent(out)   :: coords(:,:), eigvals(:)
        real, allocatable :: cand_coords(:,:), cand_eigvals(:), out_eigvals(:)
        integer :: n, ndiff_scan, nmodes, max_cand, ncand
        n = graph%n
        if( n < 3 )then
            allocate(coords(1,n), eigvals(1), source=0.)
            return
        endif
        if( .not. graph%has_theta() .or. .not. graph%has_shift() ) THROW_HARD('SE2 graph embedding requires theta/shift payloads')
        if( ndiff_req <= 0 )then
            ndiff_scan = min(16, max(1, n - 2))
        else
            ndiff_scan = min(max(1, ndiff_req), max(1, n - 2))
        endif
        nmodes = min(max(0, nmodes_req), max(0, n - 1))
        max_cand = ndiff_scan * (1 + nmodes + (nmodes + 1) * SE2_NTRANS_MODES)
        allocate(cand_coords(n,max(1,max_cand)), cand_eigvals(max(1,max_cand)), source=0.)
        call collect_se2_candidates(graph, nmodes, ndiff_scan, max(shift_scale, 1.e-6), cand_coords, cand_eigvals, ncand)
        call select_output_candidates(cand_coords, cand_eigvals, ncand, ndiff_scan, coords, out_eigvals)
        call normalize_coords(coords)
        allocate(eigvals(size(out_eigvals)), source=out_eigvals)
        deallocate(cand_coords, cand_eigvals, out_eigvals)
    end subroutine embed_se2_graph

    subroutine collect_so2_candidates(graph, nmodes, ndiff_scan, cand_coords, cand_eigvals, ncand)
        type(diffmap_graph), target, intent(in) :: graph
        integer, intent(in) :: nmodes, ndiff_scan
        real, intent(inout) :: cand_coords(:,:), cand_eigvals(:)
        integer, intent(out) :: ncand
        type(connection_context) :: ctx
        real, allocatable :: evals(:), evecs(:,:), feat(:)
        integer :: n, n2, nev, mode, idx, k, info, max_basis
        n = graph%n
        n2 = 2 * n
        ncand = 0
        allocate(feat(n))
        nev = min(ndiff_scan + 1, n)
        allocate(evals(nev), evecs(n,nev))
        max_basis = min(n, max(160, 8 * nev + 80))
        call sparse_eigh(graph_matvec, graph, n, nev, evals, evecs, tol=1.e-5, max_basis=max_basis, info=info)
        do k = 1,min(ndiff_scan, nev - 1)
            idx = nev - k
            feat = evals(idx) * evecs(:,idx)
            call add_steerable_candidate(feat, evals(idx), cand_coords, cand_eigvals, ncand)
        end do
        deallocate(evals, evecs)
        ctx%graph => graph
        do mode = 1,nmodes
            if( ncand >= size(cand_eigvals) ) exit
            ctx%mode = mode
            ctx%use_shift = .false.
            nev = min(2 * ndiff_scan + 2, n2)
            allocate(evals(nev), evecs(n2,nev))
            max_basis = min(n2, max(160, 8 * nev + 80))
            call sparse_eigh(connection_matvec, ctx, n2, nev, evals, evecs, tol=1.e-5, max_basis=max_basis, info=info)
            do idx = nev,1,-1
                feat = evals(idx) * sqrt(evecs(1:n,idx)**2 + evecs(n+1:n2,idx)**2)
                call add_steerable_candidate(feat, evals(idx), cand_coords, cand_eigvals, ncand)
                if( ncand >= size(cand_eigvals) ) exit
            end do
            deallocate(evals, evecs)
        end do
        deallocate(feat)
    end subroutine collect_so2_candidates

    subroutine collect_se2_candidates(graph, nmodes, ndiff_scan, shift_scale, cand_coords, cand_eigvals, ncand)
        type(diffmap_graph), target, intent(in) :: graph
        integer, intent(in) :: nmodes, ndiff_scan
        real,    intent(in) :: shift_scale
        real,    intent(inout) :: cand_coords(:,:), cand_eigvals(:)
        integer, intent(out) :: ncand
        type(connection_context) :: ctx
        real, allocatable :: evals(:), evecs(:,:), feat(:)
        real :: kx(SE2_NTRANS_MODES), ky(SE2_NTRANS_MODES), trans_factor
        integer :: n, n2, nev, mode, tmode, idx, info, max_basis
        call collect_so2_candidates(graph, nmodes, ndiff_scan, cand_coords, cand_eigvals, ncand)
        if( ncand >= size(cand_eigvals) ) return
        n = graph%n
        n2 = 2 * n
        kx = [1., 0.]
        ky = [0., 1.]
        trans_factor = TWOPI / max(2. * shift_scale, 1.e-6)
        allocate(feat(n))
        ctx%graph => graph
        ctx%use_shift = .true.
        ctx%trans_factor = trans_factor
        do mode = 0,nmodes
            ctx%mode = mode
            do tmode = 1,SE2_NTRANS_MODES
                if( ncand >= size(cand_eigvals) ) exit
                ctx%kx = kx(tmode)
                ctx%ky = ky(tmode)
                nev = min(2 * ndiff_scan + 2, n2)
                allocate(evals(nev), evecs(n2,nev))
                max_basis = min(n2, max(160, 8 * nev + 80))
                call sparse_eigh(connection_matvec, ctx, n2, nev, evals, evecs, tol=1.e-5, max_basis=max_basis, info=info)
                do idx = nev,1,-1
                    feat = evals(idx) * sqrt(evecs(1:n,idx)**2 + evecs(n+1:n2,idx)**2)
                    call add_steerable_candidate(feat, evals(idx), cand_coords, cand_eigvals, ncand)
                    if( ncand >= size(cand_eigvals) ) exit
                end do
                deallocate(evals, evecs)
            end do
            if( ncand >= size(cand_eigvals) ) exit
        end do
        deallocate(feat)
    end subroutine collect_se2_candidates

    subroutine connection_matvec(ctx_any, x, y)
        class(*), intent(in)  :: ctx_any
        real,     intent(in)  :: x(:)
        real,     intent(out) :: y(:)
        type(diffmap_graph), pointer :: graph
        real :: w, phase, c, s
        integer :: n, i, j, p, mode
        real :: kx, ky, trans_factor
        logical :: use_shift
        select type(ctx => ctx_any)
        type is (connection_context)
            graph => ctx%graph
            mode = ctx%mode
            kx = ctx%kx
            ky = ctx%ky
            trans_factor = ctx%trans_factor
            use_shift = ctx%use_shift
        class default
            THROW_HARD('invalid connection matvec context')
        end select
        n = graph%n
        if( size(x) /= 2*n .or. size(y) /= 2*n ) THROW_HARD('connection matvec shape mismatch')
        y = 0.
        do i = 1,n
            do p = graph%rowptr(i), graph%rowptr(i+1) - 1
                j = graph%colind(p)
                if( allocated(graph%wnorm) )then
                    w = graph%wnorm(p)
                else
                    w = graph%w(p)
                endif
                phase = real(mode) * graph%theta(p)
                if( use_shift ) phase = phase + trans_factor * (kx * graph%shift_x(p) + ky * graph%shift_y(p))
                c = cos(phase)
                s = sin(phase)
                y(i)   = y(i)   + w * ( c * x(j) - s * x(n+j) )
                y(n+i) = y(n+i) + w * ( s * x(j) + c * x(n+j) )
            end do
        end do
    end subroutine connection_matvec

    subroutine add_steerable_candidate(feat, eigval, cand_coords, cand_eigvals, ncand)
        real,    intent(in)    :: feat(:), eigval
        real,    intent(inout) :: cand_coords(:,:), cand_eigvals(:)
        integer, intent(inout) :: ncand
        real :: mu, sigma, corr
        integer :: i
        if( ncand >= size(cand_eigvals) ) return
        mu = sum(feat) / real(size(feat))
        sigma = sqrt(sum((feat - mu)**2) / real(max(size(feat)-1, 1)))
        if( sigma < 1.e-6 ) return
        do i = 1,ncand
            corr = pearsn(feat, cand_coords(:,i))
            if( abs(corr) > 0.999 ) return
        end do
        ncand = ncand + 1
        cand_coords(:,ncand) = feat
        cand_eigvals(ncand)  = eigval
    end subroutine add_steerable_candidate

    subroutine select_output_candidates(cand_coords, cand_eigvals, ncand, ndiff_scan, coords, eigvals)
        real,              intent(in)  :: cand_coords(:,:), cand_eigvals(:)
        integer,           intent(in)  :: ncand, ndiff_scan
        real, allocatable, intent(out) :: coords(:,:), eigvals(:)
        logical, allocatable :: used(:)
        real :: best
        integer :: n, nkeep, k, i, best_i
        n = size(cand_coords,1)
        if( ncand < 1 )then
            allocate(coords(1,n), eigvals(1), source=0.)
            return
        endif
        nkeep = min(ndiff_scan, ncand)
        allocate(coords(nkeep,n), eigvals(nkeep), source=0.)
        allocate(used(ncand), source=.false.)
        do k = 1,nkeep
            best = -huge(best)
            best_i = 0
            do i = 1,ncand
                if( used(i) ) cycle
                if( cand_eigvals(i) > best )then
                    best = cand_eigvals(i)
                    best_i = i
                endif
            end do
            if( best_i == 0 ) exit
            used(best_i) = .true.
            coords(k,:) = cand_coords(:,best_i)
            eigvals(k)  = cand_eigvals(best_i)
        end do
        deallocate(used)
    end subroutine select_output_candidates

    integer function auto_ndiff_from_eigengap(eigvals) result(ndiff)
        real, intent(in) :: eigvals(:)
        real :: gap, best_gap
        integer :: k, n
        n = size(eigvals)
        ndiff = min(4, max(1, n))
        if( n <= 1 ) return
        best_gap = -huge(1.)
        do k = 1,n-1
            gap = abs(eigvals(k) - eigvals(k+1))
            if( gap > best_gap )then
                best_gap = gap
                ndiff = k
            endif
        end do
        ndiff = min(max(1, ndiff), n)
    end function auto_ndiff_from_eigengap

    subroutine normalize_coords(coords, coord_mu, coord_sigma)
        real,              intent(inout) :: coords(:,:)
        real, allocatable, optional, intent(out) :: coord_mu(:), coord_sigma(:)
        real, allocatable :: mu(:), sigma(:)
        integer :: i, n
        n = size(coords,2)
        if( n < 1 ) return
        allocate(mu(size(coords,1)), sigma(size(coords,1)), source=0.)
        do i = 1,size(coords,1)
            mu(i) = sum(coords(i,:)) / real(n)
            sigma(i) = sqrt(sum((coords(i,:) - mu(i))**2) / real(max(n - 1, 1)))
            if( sigma(i) < 1.e-6 ) sigma(i) = 1.
            coords(i,:) = (coords(i,:) - mu(i)) / sigma(i)
        end do
        if( present(coord_mu)    ) allocate(coord_mu(size(mu)), source=mu)
        if( present(coord_sigma) ) allocate(coord_sigma(size(sigma)), source=sigma)
        deallocate(mu, sigma)
    end subroutine normalize_coords

end module simple_diffusion_maps
