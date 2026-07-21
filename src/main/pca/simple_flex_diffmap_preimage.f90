!@descr: representative diffusion-manifold targets and soft pre-image weights
module simple_flex_diffmap_preimage
use simple_core_module_api
use simple_clustering_utils, only: cluster_coords_kmedoids_coverage
use simple_stat,             only: median
implicit none
private
#include "simple_local_flags.inc"

public :: select_flex_diffmap_preimages

contains

    subroutine select_flex_diffmap_preimages( pinds, raw_coords, nstates_req, medoids, labels, weights, neff, bandwidths )
        integer,              intent(in)  :: pinds(:),nstates_req
        real,                 intent(in)  :: raw_coords(:,:)
        integer, allocatable, intent(out) :: medoids(:),labels(:)
        real,    allocatable, intent(out) :: weights(:,:),neff(:),bandwidths(:)
        real, allocatable :: radii(:)
        real(dp), allocatable :: nearest2(:)
        real(dp) :: d2,sumw,sumw2,bw2
        integer :: n,nstates,i,j,k,nmembers,u
        n=size(raw_coords,1)
        if( size(pinds)/=n .or. size(raw_coords,2)<1 ) THROW_HARD('invalid flex pre-image coordinate table')
        nstates=min(max(2,nstates_req),n)
        if( nstates/=nstates_req )then
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX PRE-IMAGE requested_states=',nstates_req,' clamped_states=',nstates
        endif
        call cluster_coords_kmedoids_coverage(raw_coords,nstates,medoids,labels)
        allocate(weights(n,nstates),neff(nstates),bandwidths(nstates),source=0.)
        allocate(nearest2(n),source=0._dp)
        do k=1,nstates
            nmembers=count(labels==k)
            allocate(radii(nmembers))
            j=0
            do i=1,n
                d2=sum((real(raw_coords(i,:),dp)-real(raw_coords(medoids(k),:),dp))**2)
                nearest2(i)=d2
                if( labels(i)==k )then
                    j=j+1; radii(j)=sqrt(real(max(0._dp,d2),sp))
                endif
            end do
            bandwidths(k)=median(radii)
            if( bandwidths(k)<=real(DTINY) )then
                if( any(nearest2>real(DTINY,dp)) )then
                    bandwidths(k)=sqrt(real(minval(nearest2,mask=nearest2>real(DTINY,dp)),sp))
                else
                    THROW_HARD('degenerate flex diffusion coordinates')
                endif
            endif
            if( bandwidths(k)<=real(DTINY) ) THROW_HARD('degenerate flex pre-image bandwidth')
            bw2=real(bandwidths(k),dp)**2
            do i=1,n
                weights(i,k)=real(exp(-0.5_dp*nearest2(i)/bw2),sp)
            end do
            sumw=sum(real(weights(:,k),dp)); sumw2=sum(real(weights(:,k),dp)**2)
            if( sumw<=DTINY .or. sumw2<=DTINY ) THROW_HARD('empty flex pre-image kernel')
            neff(k)=real(sumw*sumw/sumw2,sp)
            deallocate(radii)
        end do
        call del_file('flex_diffmap_preimages.txt')
        open(newunit=u,file='flex_diffmap_preimages.txt',status='replace',action='write')
        write(u,'(A)') '# state medoid_row particle hard_population bandwidth effective_population'
        do k=1,nstates
            write(u,'(I6,1X,I10,1X,I10,1X,I10,1X,ES16.8,1X,ES16.8)') k,medoids(k),pinds(medoids(k)), &
                &count(labels==k),bandwidths(k),neff(k)
            write(logfhandle,'(A,I0,A,I0,A,I0,A,F10.2,A,F10.3)') '>>> FLEX PRE-IMAGE state=',k, &
                &' exemplar_particle=',pinds(medoids(k)),' hard_population=',count(labels==k), &
                &' effective_population=',neff(k),' bandwidth=',bandwidths(k)
        end do
        close(u)
        deallocate(nearest2)
    end subroutine select_flex_diffmap_preimages

end module simple_flex_diffmap_preimage
