!@descr: clustering of a distance matrix by algorithm name, and the seeding shared by the clustering modules
module simple_clustering_utils
use simple_kmedoids,      only: kmedoids
use simple_aff_prop,      only: aff_prop
use simple_hac,           only: hac
use simple_stat,          only: calc_ap_pref
use simple_srch_sort_loc, only: hpsort
use simple_core_module_api
implicit none

public :: cluster_dmat, farthest_point_seeds, equal_mass_quantile_start
private
#include "simple_local_flags.inc"

contains

    subroutine cluster_dmat( dmat, algorithm, nclust, i_medoids, labels, ap_pref, nclust_max )
        real,                 intent(in)    :: dmat(:,:)
        character(len=*),     intent(in)    :: algorithm
        integer,              intent(inout) :: nclust
        integer, allocatable, intent(inout) :: i_medoids(:), labels(:)
        real,    optional,    intent(in)    :: ap_pref
        integer, optional,    intent(in)    :: nclust_max
        real,    allocatable :: smat(:,:)
        type(kmedoids)       :: kmed
        type(aff_prop)       :: aprop
        type(hac)            :: avglink
        integer :: n
        real    :: pref, simsum
        n = size(dmat, dim=1)
        if( allocated(i_medoids) ) deallocate(i_medoids)
        select case(trim(algorithm))
            case('aprop')
                if( allocated(labels) ) deallocate(labels)
                write(logfhandle,'(A)') '>>> CLUSTERING DISTANCE MATRIX WITH AFFINITY PROPAGATION'
                smat = dmat2smat(dmat)
                pref = calc_ap_pref(smat, 'median')
                if( present(ap_pref) ) pref = ap_pref
                write(logfhandle,'(A,ES14.6)') '>>> AFFINITY PROPAGATION PREFERENCE: ', pref
                call aprop%new(n, smat, pref=pref)
                call aprop%propagate(i_medoids, labels, simsum)
                call aprop%kill
                nclust = size(i_medoids)
                write(logfhandle,'(A,I3)') '>>> # CLUSTERS FOUND BY AFFINITY PROPAGATION (AP): ', nclust
                call merge_if_necessary
            case('avglink')
                write(logfhandle,'(A,I0,A)') '>>> CLUSTERING DISTANCE MATRIX WITH AVERAGE LINKAGE INTO ', nclust, ' CLUSTERS'
                call avglink%new(n, real(dmat,dp), 'average', nclust=nclust)
                call avglink%cluster(i_medoids, labels)
                call avglink%kill
            case('kmed')
                if( allocated(labels) ) deallocate(labels)
                write(logfhandle,'(A)') '>>> CLUSTERING DISTANCE MATRIX WITH K-MEDOIDS'
                if( nclust < 2 ) THROW_HARD('Invalid nclust input')
                call kmed%new(n, dmat, nclust)
                call kmed%init
                call kmed%cluster
                allocate(labels(n), i_medoids(nclust), source=0)
                call kmed%get_labels(labels)
                call kmed%get_medoids(i_medoids)
                call kmed%kill
            case DEFAULT
                THROW_HARD('unsupported algorithm flag: '//trim(algorithm))
        end select

        contains

            subroutine merge_if_necessary
                if( present(nclust_max) )then
                    if( nclust > nclust_max )then
                        call kmed%new(labels, dmat)
                        nclust = nclust_max
                        call kmed%merge(nclust)
                        if( allocated(i_medoids) ) deallocate(i_medoids)
                        if( allocated(labels)    ) deallocate(labels) 
                        allocate(i_medoids(nclust), labels(n), source=0)
                        call kmed%get_labels(labels)
                        call kmed%get_medoids(i_medoids)
                        call kmed%kill
                    endif
                endif
            end subroutine merge_if_necessary

    end subroutine cluster_dmat

    !>  \brief  farthest-point seeding of the rows of X(n,d): seeds(1) = first, then each next seed is the
    !!          point farthest from every seed chosen so far (squared distance, per-dimension weights wdim when
    !!          given), the first index on a tie
    subroutine farthest_point_seeds( X, k, first, seeds, wdim )
        real(dp),           intent(in)  :: X(:,:)
        integer,            intent(in)  :: k, first
        integer,            intent(out) :: seeds(k)
        real(dp), optional, intent(in)  :: wdim(:)
        real(dp), allocatable :: mind(:), w(:)
        real(dp) :: d2, dmax
        integer  :: n, d, i, q, s
        n = size(X,1)
        d = size(X,2)
        if( k < 1 .or. k > n )       THROW_HARD('k must be in [1,n]; farthest_point_seeds')
        if( first < 1 .or. first > n ) THROW_HARD('first seed out of range; farthest_point_seeds')
        allocate(w(d), source=1._dp)
        if( present(wdim) ) w = wdim
        allocate(mind(n), source=huge(1._dp))
        seeds(1) = first
        do s = 2, k
            dmax = -1._dp
            seeds(s) = 1
            do i = 1, n
                d2 = 0._dp
                do q = 1, d
                    d2 = d2 + w(q)*(X(i,q) - X(seeds(s-1),q))**2
                end do
                mind(i) = min(mind(i), d2)
                if( mind(i) > dmax )then
                    dmax     = mind(i)
                    seeds(s) = i
                endif
            end do
        end do
        deallocate(mind, w)
    end subroutine farthest_point_seeds

    !>  \brief  equal-mass quantile start of k means: the rows of X(n,d) sorted by key (in single precision,
    !!          hpsort) are cut into k consecutive blocks of (nearly) equal size, and means(:,c) is the mean of
    !!          block c, the lowest keys first
    subroutine equal_mass_quantile_start( X, key, k, means )
        real(dp), intent(in)  :: X(:,:)     !< (n,d)
        real(dp), intent(in)  :: key(:)     !< (n) sort key, e.g. the projection on a dominant axis
        integer,  intent(in)  :: k
        real(dp), intent(out) :: means(:,:) !< (d,k)
        real,    allocatable :: skey(:)
        integer, allocatable :: ord(:)
        integer :: n, i, c, i0, i1
        n = size(X,1)
        if( size(key) /= n )  THROW_HARD('key must have one entry per point; equal_mass_quantile_start')
        if( k < 1 .or. k > n ) THROW_HARD('k must be in [1,n]; equal_mass_quantile_start')
        allocate(skey(n), ord(n))
        do i = 1, n
            skey(i) = real(key(i))
            ord(i)  = i
        end do
        call hpsort(skey, ord)
        means = 0._dp
        do c = 1, k
            i0 = nint(real(c-1,dp)*real(n,dp)/real(k,dp)) + 1
            i1 = max(nint(real(c,dp)*real(n,dp)/real(k,dp)), i0)
            do i = i0, i1
                means(:,c) = means(:,c) + X(ord(i),:)
            end do
            means(:,c) = means(:,c) / real(i1 - i0 + 1,dp)
        end do
        deallocate(skey, ord)
    end subroutine equal_mass_quantile_start

end module simple_clustering_utils
