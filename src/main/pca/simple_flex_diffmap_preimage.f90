!@descr: representative diffusion-manifold targets selected with SIMPLE k-medoids
module simple_flex_diffmap_preimage
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_clustering_utils, only: cluster_dmat
implicit none
private
#include "simple_local_flags.inc"

public :: select_flex_diffmap_preimages

contains

    subroutine select_flex_diffmap_preimages( pinds, coords, nstates_req, medoids, labels )
        integer,              intent(in)  :: pinds(:), nstates_req
        real,                 intent(in)  :: coords(:,:)
        integer, allocatable, intent(out) :: medoids(:), labels(:)
        real, allocatable :: dmat(:,:), coords_mode_major(:,:)
        real :: dval
        integer :: n, nstates, i, j, u
        n = size(coords,1)
        if( size(pinds) /= n .or. size(coords,2) < 1 ) THROW_HARD('invalid flex pre-image coordinate table')
        if( .not.all(ieee_is_finite(coords)) ) THROW_HARD('nonfinite flex pre-image coordinate')
        nstates = min(max(2,nstates_req),n)
        if( nstates /= nstates_req )then
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX PRE-IMAGE requested_states=',nstates_req, &
                &' clamped_states=',nstates
        endif
        allocate(dmat(n,n), source=0.)
        allocate(coords_mode_major(size(coords,2),n),source=transpose(coords))
        !$omp parallel do default(shared) private(i,j,dval) schedule(dynamic) proc_bind(close)
        do i = 1,n-1
            do j = i+1,n
                dval = euclid(coords_mode_major(:,i),coords_mode_major(:,j))
                dmat(i,j) = dval
                dmat(j,i) = dval
            end do
        end do
        !$omp end parallel do
        call cluster_dmat(dmat,'kmed',nstates,medoids,labels)
        deallocate(dmat,coords_mode_major)
        call del_file('flex_diffmap_preimages.txt')
        open(newunit=u,file='flex_diffmap_preimages.txt',status='replace',action='write')
        write(u,'(A)') '# state medoid_row particle hard_population'
        do i = 1,nstates
            write(u,'(I6,1X,I10,1X,I10,1X,I10)') i,medoids(i),pinds(medoids(i)),count(labels==i)
            write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX PRE-IMAGE state=',i, &
                &' exemplar_particle=',pinds(medoids(i)),' hard_population=',count(labels==i)
        end do
        close(u)
    end subroutine select_flex_diffmap_preimages

end module simple_flex_diffmap_preimage
