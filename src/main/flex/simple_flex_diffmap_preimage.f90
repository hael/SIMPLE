!@descr: flex-analysis diffusion-manifold targets selected with SIMPLE k-medoids
module simple_flex_diffmap_preimage
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_clustering_utils, only: cluster_dmat
use simple_srch_sort_loc, only: hpsort
implicit none
private
#include "simple_local_flags.inc"

public :: select_flex_diffmap_preimages
public :: build_flex_preimage_kernel_weights

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

    !> Build soft state weights on the diffusion manifold using medoid
    !! descriptors as target points. Each column of `weights` is normalized to
    !! one and can be used directly as weighted reconstruction coefficients.
    subroutine build_flex_preimage_kernel_weights( coords, medoids, weights, bandwidths, neff )
        real,                 intent(in)  :: coords(:,:)
        integer,              intent(in)  :: medoids(:)
        real, allocatable,    intent(out) :: weights(:,:), bandwidths(:), neff(:)
        real(dp), allocatable :: d2(:,:), qdens(:), base_bandwidths(:)
        real(dp) :: dval, eps, wsum, qsum, qval, qfloor, neff_floor, uniform_pop, bw_scale
        integer :: n, nstates, i, j, s, pos_count, ibw, ibw_used
        logical :: has_zero
        integer, parameter :: BW_ADAPT_MAX_ITERS=10
        n=size(coords,1)
        nstates=size(medoids)
        if( n<1 .or. size(coords,2)<1 .or. nstates<1 ) THROW_HARD('invalid flex kernel pre-image table dimensions')
        if( .not.all(ieee_is_finite(coords)) ) THROW_HARD('nonfinite flex kernel pre-image coordinate')
        allocate(weights(n,nstates),bandwidths(nstates),neff(nstates),source=0.)
        allocate(d2(n,nstates),qdens(n),base_bandwidths(nstates),source=0.d0)
        do s=1,nstates
            if( medoids(s)<1 .or. medoids(s)>n ) THROW_HARD('flex kernel pre-image medoid outside coordinate table')
            do i=1,n
                dval=euclid(coords(i,:),coords(medoids(s),:))
                d2(i,s)=dval*dval
            end do
        end do
        ! Local kernel width from median positive distance per state.
        do s=1,nstates
            pos_count=0
            do i=1,n
                if( d2(i,s)>DTINY ) pos_count=pos_count+1
            end do
            if( pos_count>0 )then
                eps=median_positive_d2(d2(:,s))
            else
                eps=1.d0
            endif
            if( eps<=DTINY ) eps=max(sum(d2(:,s))/real(max(1,n),dp),1.d-6)
            base_bandwidths(s)=eps
        end do
        ! Enforce a minimum soft-neighborhood size with global bandwidth inflation.
        uniform_pop=real(n,dp)/real(max(1,nstates),dp)
        neff_floor=max(32.d0,min(256.d0,0.35d0*uniform_pop))
        bw_scale=1.d0
        do ibw=1,BW_ADAPT_MAX_ITERS
            bandwidths=real(base_bandwidths*bw_scale)
            ! Density correction proxy q(i)=sum_s exp(-d2/eps_s) to reduce oversampled clumps.
            qsum=0.d0
            do i=1,n
                qval=0.d0
                do s=1,nstates
                    qval=qval+exp(-d2(i,s)/max(real(bandwidths(s),dp),DTINY))
                end do
                qsum=qsum+qval
                qdens(i)=qval
            end do
            qfloor=max(DTINY,0.01d0*qsum/real(max(1,n),dp))
            do s=1,nstates
                wsum=0.d0
                do i=1,n
                    qval=max(qdens(i),qfloor)
                    weights(i,s)=real(exp(-d2(i,s)/max(real(bandwidths(s),dp),DTINY))/qval)
                    wsum=wsum+real(weights(i,s),dp)
                end do
                if( wsum<=DTINY ) THROW_HARD('degenerate flex kernel pre-image weights')
                weights(:,s)=weights(:,s)/real(wsum)
                neff(s)=1./sum(weights(:,s)*weights(:,s))
            end do
            if( minval(real(neff,dp))>=neff_floor ) exit
            bw_scale=bw_scale*1.5d0
        end do
        ibw_used=min(ibw,BW_ADAPT_MAX_ITERS)
        has_zero=.false.
        do s=1,nstates
            if( neff(s)<2. ) has_zero=.true.
        end do
        write(logfhandle,'(A,F8.2,A,F8.2,A,F8.2,A,I0,A,L1)') '>>> FLEX PRE-IMAGE kernel NEFF min=',minval(neff), &
            &' max=',maxval(neff),' floor=',real(neff_floor),' bw_iters=',ibw_used,' sparse_state_warning=',has_zero
        write(logfhandle,'(A,ES12.4,A,ES12.4)') '>>> FLEX PRE-IMAGE bandwidth scale=',bw_scale, &
            &' median=',median_default_real(bandwidths)
        deallocate(d2,qdens,base_bandwidths)
    contains
        function median_positive_d2( vals ) result( med )
            real(dp), intent(in) :: vals(:)
            real(dp) :: med
            real, allocatable :: work(:)
            integer :: i, m, k
            m=0
            do i=1,size(vals)
                if( vals(i)>DTINY ) m=m+1
            end do
            if( m<1 )then
                med=0.d0
                return
            endif
            allocate(work(m))
            k=0
            do i=1,size(vals)
                if( vals(i)>DTINY )then
                    k=k+1
                    work(k)=real(vals(i))
                endif
            end do
            call hpsort(work)
            if( mod(m,2)==1 )then
                med=real(work((m+1)/2),dp)
            else
                med=0.5d0*real(work(m/2)+work(m/2+1),dp)
            endif
            deallocate(work)
        end function median_positive_d2

        function median_default_real( vals ) result( med )
            real, intent(in) :: vals(:)
            real :: med
            real, allocatable :: work(:)
            integer :: nvals
            nvals=size(vals)
            if( nvals<1 )then
                med=0.
                return
            endif
            allocate(work(nvals),source=vals)
            call hpsort(work)
            if( mod(nvals,2)==1 )then
                med=work((nvals+1)/2)
            else
                med=0.5*(work(nvals/2)+work(nvals/2+1))
            endif
            deallocate(work)
        end function median_default_real
    end subroutine build_flex_preimage_kernel_weights

end module simple_flex_diffmap_preimage
