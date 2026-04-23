!@descr: diffusion maps embedding for nonlinear class splitting
module simple_diffusion_maps
use simple_core_module_api
use simple_linalg, only: eigh
implicit none
#include "simple_local_flags.inc"

private
public :: diffusion_map_embedder

type :: diffusion_map_embedder
    integer :: ndiff = 4  !< 0 requests automatic dimensionality selection from the diffusion spectrum
    integer :: k_nn  = 10
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
        real, allocatable, optional,    intent(out)   :: eigvals(:)
        real, allocatable :: d2(:,:), aff(:,:), evals(:), evecs(:,:), row(:), kth_d2(:), deg(:), norm_aff(:,:)
        real, allocatable :: diff_evals(:)
        real :: eps, scale
        integer :: nptcls, npix, i, j, k, ndiff_used, ndiff_scan, k_used, nev
        integer(int64) :: t0, t1, trate, t_total0
        npix   = size(pcavecs, 1)
        nptcls = size(pcavecs, 2)
        if( nptcls < 3 )then
            allocate(coords(1, nptcls), source=0.)
            return
        endif
        if( self%ndiff <= 0 )then
            ndiff_scan = min(16, max(1, nptcls - 2))
        else
            ndiff_scan = min(max(1, self%ndiff), max(1, nptcls - 2))
        endif
        ndiff_used = ndiff_scan
        k_used     = min(max(2, self%k_nn), nptcls - 1)
        call system_clock(t_total0, trate)
        write(logfhandle,'(A,I8,A,I8,A,I8,A,I8)') 'Diffusion maps start: N=', nptcls, ' D=', npix, ' ndiff=', ndiff_scan, ' k_nn=', k_used
        call flush(logfhandle)
        allocate(d2(nptcls,nptcls), aff(nptcls,nptcls), kth_d2(nptcls), deg(nptcls), source=0.)
        call system_clock(t0)
        do i = 1, nptcls - 1
            do j = i + 1, nptcls
                scale   = euclid(pcavecs(:,i), pcavecs(:,j))
                d2(i,j) = scale * scale
                d2(j,i) = d2(i,j)
            end do
        end do
        call system_clock(t1)
        write(logfhandle,'(A,F8.3,A)') 'Diffusion maps pairwise distances: ', real(t1-t0)/real(trate), ' s'
        call flush(logfhandle)
        call system_clock(t0)
        do i = 1, nptcls
            allocate(row(nptcls-1))
            k = 0
            do j = 1, nptcls
                if( j == i ) cycle
                k = k + 1
                row(k) = d2(i,j)
            end do
            call sort_reals_asc(row)
            kth_d2(i) = row(k_used)
            deallocate(row)
        end do
        eps = median_real(kth_d2)
        if( eps < DTINY ) eps = max(sum(kth_d2) / real(max(size(kth_d2),1)), 1.e-6)
        call system_clock(t1)
        write(logfhandle,'(A,F8.3,A,ES10.3)') 'Diffusion maps local scale selection: ', real(t1-t0)/real(trate), ' s; eps=', eps
        call flush(logfhandle)
        call system_clock(t0)
        do i = 1, nptcls
            do j = 1, nptcls
                if( i == j ) cycle
                if( d2(i,j) <= kth_d2(i) .or. d2(i,j) <= kth_d2(j) )then
                    aff(i,j) = exp(-d2(i,j) / eps)
                endif
            end do
        end do
        aff = 0.5 * (aff + transpose(aff))
        deg = sum(aff, dim=2)
        do i = 1, nptcls
            if( deg(i) > DTINY ) cycle
            j = nearest_neighbor_index(d2(i,:), i)
            aff(i,j) = 1.
            aff(j,i) = 1.
        end do
        deg = sum(aff, dim=2)
        allocate(norm_aff(nptcls,nptcls), source=0.)
        do i = 1, nptcls
            do j = 1, nptcls
                if( aff(i,j) <= DTINY ) cycle
                norm_aff(i,j) = aff(i,j) / sqrt(max(deg(i), DTINY) * max(deg(j), DTINY))
            end do
        end do
        call system_clock(t1)
        write(logfhandle,'(A,F8.3,A)') 'Diffusion maps graph/normalization: ', real(t1-t0)/real(trate), ' s'
        call flush(logfhandle)
        nev = ndiff_scan + 1
        allocate(evals(nev), evecs(nptcls,nev))
        call system_clock(t0)
        call eigh(nptcls, norm_aff, nev, evals, evecs)
        call system_clock(t1)
        write(logfhandle,'(A,F8.3,A,I8)') 'Diffusion maps eigensolve: ', real(t1-t0)/real(trate), ' s; nev=', nev
        call flush(logfhandle)
        call system_clock(t0)
        allocate(diff_evals(ndiff_scan), source=0.)
        do k = 1, ndiff_scan
            diff_evals(k) = evals(nev - k)
        end do
        if( self%ndiff <= 0 ) ndiff_used = auto_ndiff_from_eigengap(diff_evals)
        allocate(coords(ndiff_used, nptcls), source=0.)
        if( present(eigvals) ) allocate(eigvals(ndiff_scan), source=diff_evals)
        do k = 1, ndiff_used
            j = nev - k
            coords(k,:) = evals(j) * evecs(:,j)
        end do
        call normalize_coords(coords)
        call system_clock(t1)
        write(logfhandle,'(A,F8.3,A,I8)') 'Diffusion maps embedding build: ', real(t1-t0)/real(trate), ' s; dims=', ndiff_used
        call flush(logfhandle)
        write(logfhandle,'(A,F8.3,A,I8,A,ES10.3,A,ES10.3)') 'Diffusion maps total: ', real(t1-t_total0)/real(trate), &
            ' s; dims=', ndiff_used, ' lambda1=', diff_evals(1), ' lambda_last=', diff_evals(ndiff_scan)
        call flush(logfhandle)
        deallocate(d2, aff, kth_d2, deg, norm_aff, evals, evecs, diff_evals)
    end subroutine embed

    integer function auto_ndiff_from_eigengap(eigvals) result(ndiff)
        real, intent(in) :: eigvals(:)
        real :: best_gap, gap
        integer :: i, maxcand
        maxcand = size(eigvals)
        if( maxcand <= 1 )then
            ndiff = 1
            return
        endif
        best_gap = -huge(best_gap)
        ndiff = min(4, maxcand)
        do i = 1, maxcand - 1
            gap = eigvals(i) - eigvals(i+1)
            if( gap > best_gap )then
                best_gap = gap
                ndiff = i
            endif
        end do
        ! Avoid collapsing the embedding to one coordinate unless the spectrum is tiny.
        ndiff = max(min(2, maxcand), min(maxcand, ndiff))
    end function auto_ndiff_from_eigengap

    subroutine normalize_coords(coords)
        real, intent(inout) :: coords(:,:)
        real :: mu, sigma
        integer :: i
        do i = 1, size(coords,1)
            mu = sum(coords(i,:)) / real(size(coords,2))
            sigma = sqrt(sum((coords(i,:) - mu)**2) / real(max(size(coords,2)-1, 1)))
            if( sigma < 1.e-6 ) sigma = 1.0
            coords(i,:) = (coords(i,:) - mu) / sigma
        end do
    end subroutine normalize_coords

    integer function nearest_neighbor_index(row, self_idx) result(ind)
        real,    intent(in) :: row(:)
        integer, intent(in) :: self_idx
        real :: best
        integer :: i
        ind  = 1
        best = huge(best)
        do i = 1, size(row)
            if( i == self_idx ) cycle
            if( row(i) < best )then
                best = row(i)
                ind = i
            endif
        end do
    end function nearest_neighbor_index

    real function median_real(vals) result(med)
        real, intent(in) :: vals(:)
        real, allocatable :: work(:)
        integer :: n
        n = size(vals)
        allocate(work(n), source=vals)
        call sort_reals_asc(work)
        if( mod(n,2) == 0 )then
            med = 0.5 * (work(n/2) + work(n/2 + 1))
        else
            med = work((n + 1) / 2)
        endif
        deallocate(work)
    end function median_real

    subroutine sort_reals_asc(vals)
        real, intent(inout) :: vals(:)
        integer :: i, j
        real :: tmp
        do i = 1, size(vals) - 1
            do j = i + 1, size(vals)
                if( vals(j) < vals(i) )then
                    tmp = vals(i)
                    vals(i) = vals(j)
                    vals(j) = tmp
                endif
            end do
        end do
    end subroutine sort_reals_asc

end module simple_diffusion_maps
