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
public :: test_flex_preimage_bandwidth_decoupling

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
        real(dp), allocatable :: d2(:,:), qdens(:), base_bandwidths(:), bw_scales(:)
        real(dp) :: dval, eps, wsum, qsum, qval, qfloor, neff_floor, uniform_pop, raw_neff
        integer :: n, nstates, i, j, s, pos_count, ibw
        logical :: has_zero
        integer, parameter :: BW_ADAPT_MAX_ITERS=10
        n=size(coords,1)
        nstates=size(medoids)
        if( n<1 .or. size(coords,2)<1 .or. nstates<1 ) THROW_HARD('invalid flex kernel pre-image table dimensions')
        if( .not.all(ieee_is_finite(coords)) ) THROW_HARD('nonfinite flex kernel pre-image coordinate')
        allocate(weights(n,nstates),bandwidths(nstates),neff(nstates),source=0.)
        allocate(d2(n,nstates),qdens(n),base_bandwidths(nstates),bw_scales(nstates),source=0.d0)
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
        ! Enforce a minimum soft-neighborhood size with *per-state* bandwidth
        ! inflation. A sparse state must not broaden the bandwidth of dense
        ! states: each state independently grows its own kernel width until its
        ! own single-state effective sample size reaches the floor (or the
        ! iteration cap). Selection uses the raw (pre-density-correction) kernel
        ! so states remain fully decoupled at the bandwidth-selection stage.
        uniform_pop=real(n,dp)/real(max(1,nstates),dp)
        neff_floor=max(32.d0,min(256.d0,0.35d0*uniform_pop))
        do s=1,nstates
            bw_scales(s)=1.d0
            do ibw=1,BW_ADAPT_MAX_ITERS
                raw_neff=raw_neff_state(d2(:,s),base_bandwidths(s)*bw_scales(s))
                if( raw_neff>=neff_floor ) exit
                bw_scales(s)=bw_scales(s)*1.5d0
            end do
            bandwidths(s)=real(base_bandwidths(s)*bw_scales(s))
        end do
        ! Density correction proxy q(i)=sum_s exp(-d2/eps_s) to reduce oversampled
        ! clumps, evaluated once with the final per-state bandwidths.
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
        has_zero=.false.
        do s=1,nstates
            if( neff(s)<2. ) has_zero=.true.
        end do
        write(logfhandle,'(A,F8.2,A,F8.2,A,F8.2,A,L1)') '>>> FLEX PRE-IMAGE kernel NEFF min=',minval(neff), &
            &' max=',maxval(neff),' floor=',real(neff_floor),' sparse_state_warning=',has_zero
        write(logfhandle,'(A,ES12.4,A,ES12.4,A,ES12.4)') '>>> FLEX PRE-IMAGE bandwidth scale min=',minval(bw_scales), &
            &' max=',maxval(bw_scales),' median=',median_default_real(bandwidths)
        deallocate(d2,qdens,base_bandwidths,bw_scales)
    contains
        !> Single-state effective sample size of the raw (unnormalized,
        !!
        !! non-density-corrected) Gaussian kernel. Used only to size each
        !! state's bandwidth independently of the other states.
        function raw_neff_state( d2col, bw ) result( neff_raw )
            real(dp), intent(in) :: d2col(:)
            real(dp), intent(in) :: bw
            real(dp) :: neff_raw, w, s1, s2, bw_eff
            integer  :: i
            bw_eff = max(bw,DTINY)
            s1 = 0.d0
            s2 = 0.d0
            do i=1,size(d2col)
                w  = exp(-d2col(i)/bw_eff)
                s1 = s1 + w
                s2 = s2 + w*w
            end do
            if( s2<=DTINY )then
                neff_raw = 0.d0
            else
                neff_raw = (s1*s1)/s2
            endif
        end function raw_neff_state

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

    !> Regression test for the per-state bandwidth decoupling in
    !! build_flex_preimage_kernel_weights.  A dense state's kernel bandwidth must
    !! be independent of how sparsely populated an unrelated state is: under the
    !! previous global bandwidth-inflation logic a sparse state broadened every
    !! state's bandwidth.  Two problems share an identical dense cluster (its
    !! medoid at the centre of a unit circle so every positive intra-cluster
    !! distance equals 1, making the median positive distance robust to the
    !! sparse population) and differ only in the number of sparse points.  The
    !! dense-state bandwidth must be identical across the two problems.
    subroutine test_flex_preimage_bandwidth_decoupling()
        integer, parameter :: NDENSE=90, NSPARSE_A=10, NSPARSE_B=60
        real :: bw_dense_a, bw_dense_b, rel_diff
        bw_dense_a = dense_bandwidth(NSPARSE_A)
        bw_dense_b = dense_bandwidth(NSPARSE_B)
        rel_diff   = abs(bw_dense_a - bw_dense_b) / max(abs(bw_dense_b), 1.e-30)
        write(logfhandle,'(A,ES12.4,A,ES12.4,A,ES12.4)') '>>> FLEX PRE-IMAGE bandwidth decoupling dense_bw(sparse=10)=', &
            &bw_dense_a,' dense_bw(sparse=60)=',bw_dense_b,' rel_diff=',rel_diff
        if( rel_diff > 1.e-6 ) &
            &THROW_HARD('dense-state kernel bandwidth depends on the sparse-state population (bandwidth decoupling regression)')
        write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE BANDWIDTH DECOUPLING TEST PASSED'

    contains

        !> Build a two-cluster problem with a fixed dense cluster (centre + unit
        !! circle) and nsparse far-away sparse points, then return the fitted
        !! bandwidth of the dense state.
        function dense_bandwidth( nsparse ) result( bw_dense )
            integer, intent(in) :: nsparse
            real    :: bw_dense
            real,    allocatable :: coords(:,:), weights(:,:), bandwidths(:), neff(:)
            integer, allocatable :: medoids(:)
            real(dp) :: twopi, ang, phi
            integer  :: n, i, k
            twopi = 2.d0*acos(-1.d0)
            n     = NDENSE + nsparse
            allocate(coords(n,2), source=0.)
            ! Dense cluster: index 1 is the medoid at the origin; indices 2..NDENSE
            ! lie on the unit circle so every positive distance to the medoid is 1.
            coords(1,1) = 0.
            coords(1,2) = 0.
            do i=2,NDENSE
                ang = twopi*real(i-2,dp)/real(NDENSE-1,dp)
                coords(i,1) = real(cos(ang))
                coords(i,2) = real(sin(ang))
            end do
            ! Sparse cluster: far from the dense cluster; index NDENSE+1 is its medoid.
            do k=1,nsparse
                i   = NDENSE + k
                phi = twopi*real(k-1,dp)/real(max(1,nsparse),dp)
                coords(i,1) = 100. + real(5.d0*cos(phi))
                coords(i,2) =        real(5.d0*sin(phi))
            end do
            allocate(medoids(2))
            medoids(1) = 1
            medoids(2) = NDENSE + 1
            call build_flex_preimage_kernel_weights(coords, medoids, weights, bandwidths, neff)
            bw_dense = bandwidths(1)
            deallocate(coords, medoids, weights, bandwidths, neff)
        end function dense_bandwidth

    end subroutine test_flex_preimage_bandwidth_decoupling

end module simple_flex_diffmap_preimage
