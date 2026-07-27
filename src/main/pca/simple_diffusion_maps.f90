!@descr: sparse diffusion-map embeddings for class-split graphs
module simple_diffusion_maps
use simple_core_module_api
use simple_diff_map_graphs,  only: diffmap_graph, build_euclidean_knn_graph, graph_matvec
use simple_linalg,           only: sparse_eigh
implicit none
#include "simple_local_flags.inc"

private
public :: diffusion_map_embedder
public :: embed_graph
public :: normalize_coords
public :: auto_ndiff_from_eigengap

type :: diffusion_map_embedder
    integer :: ndiff = 4
    integer :: k_nn  = 5
contains
    procedure :: set_params
    procedure :: embed
end type diffusion_map_embedder

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
        call build_euclidean_knn_graph(pcavecs, min(max(2,self%k_nn), nptcls-1), graph)
        call embed_graph(graph, self%ndiff, coords, eigs)
        if( present(eigvals) ) allocate(eigvals(size(eigs)), source=eigs)
        call graph%kill()
        if( allocated(eigs) ) deallocate(eigs)
    end subroutine embed



    subroutine embed_graph( graph, ndiff_req, coords, eigvals, raw_coords, eigenfunctions, nystrom_coords )
        type(diffmap_graph), target, intent(in)  :: graph
        integer,                   intent(in)    :: ndiff_req
        real, allocatable,         intent(out)   :: coords(:,:), eigvals(:)
        real, allocatable, optional, intent(out) :: raw_coords(:,:)
        real, allocatable, optional, intent(out) :: eigenfunctions(:,:), nystrom_coords(:,:)
        real, allocatable :: evals(:), evecs(:,:), diff_evals(:)
        real :: coeff
        integer :: n, ndiff_scan, ndiff_used, nev, eig_info, max_basis, i, k, j, p, nbr
        n = graph%n
        if( n < 3 )then
            allocate(coords(1,n), eigvals(1), source=0.)
            if( present(raw_coords) ) allocate(raw_coords(1,n), source=0.)
            if( present(eigenfunctions) ) allocate(eigenfunctions(1,n), source=0.)
            if( present(nystrom_coords) ) allocate(nystrom_coords(1,n), source=0.)
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
        if( eig_info /= 0 ) THROW_HARD('sparse eigensolve failed in diffusion-map embedding')
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
        if( present(raw_coords) ) allocate(raw_coords(ndiff_used,n), source=coords)
        if( present(eigenfunctions) )then
            allocate(eigenfunctions(ndiff_used,n), source=0.)
            do k = 1,ndiff_used
                eigenfunctions(k,:) = evecs(:,nev-k)
            end do
        endif
        if( present(nystrom_coords) )then
            allocate(nystrom_coords(ndiff_used,n), source=0.)
            do k = 1,ndiff_used
                j = nev-k
                if( abs(evals(j)) <= real(DTINY) ) cycle
                do i = 1,n
                    coeff = 0.
                    do p = graph%rowptr(i),graph%rowptr(i+1)-1
                        nbr = graph%colind(p)
                        if( nbr < 1 .or. nbr > n ) cycle
                        if( allocated(graph%wnorm) )then
                            coeff = coeff + graph%wnorm(p)*evecs(nbr,j)
                        else
                            coeff = coeff + graph%w(p)*evecs(nbr,j)
                        endif
                    end do
                    nystrom_coords(k,i) = coeff/evals(j)
                end do
            end do
        endif
        call normalize_coords(coords)
        deallocate(evals, evecs, diff_evals)
    end subroutine embed_graph





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
