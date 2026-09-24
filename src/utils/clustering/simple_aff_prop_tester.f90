!@descr: unit tests for affinity propagation (simple_aff_prop): exemplars, partition, preference limits and determinism
! Replaces the in-module test_aff_prop, which printed its verdict without recording it. Three clusters
! of twelve points each (a 4x3 grid with a small third coordinate around (0,0), (7,0) and (0,7)),
! similarity minus the squared distance. A float32 emulation of the algorithm (numpy) gives the
! exemplars 7, 19 and 31 for every preference from -1 to -200, the cluster medoids, which the test
! finds by brute force; at a preference of 0, the largest similarity, every point is an exemplar.
module simple_aff_prop_tester
use simple_test_utils
use simple_defs
use simple_aff_prop,     only: aff_prop
use simple_string_utils, only: real2str
implicit none
private
public :: run_all_aff_prop_tests

integer, parameter :: NPER = 12, NTOT = 3 * NPER

contains

    subroutine run_all_aff_prop_tests()
        write(*,'(A)') '**** running all aff_prop tests ****'
        call test_three_clusters()
        call test_preference_limits()
        call test_determinism()
    end subroutine run_all_aff_prop_tests

    subroutine make_data( simmat, truth )
        real,    intent(out) :: simmat(NTOT,NTOT)
        integer, intent(out) :: truth(NTOT)
        real    :: datavecs(NTOT,5), centers_true(3,5)
        integer :: i, j, k
        centers_true        = 0.
        centers_true(2,1:2) = [7.0, 0.0]
        centers_true(3,1:2) = [0.0, 7.0]
        do k = 1,3
            do i = 1,NPER
                j = (k - 1) * NPER + i
                truth(j)      = k
                datavecs(j,:) = centers_true(k,:)
                datavecs(j,1) = datavecs(j,1) + 0.12 * real(mod(i - 1, 4))
                datavecs(j,2) = datavecs(j,2) + 0.10 * real((i - 1) / 4)
                datavecs(j,3) = 0.03 * real(mod(i, 3))
            end do
        end do
        do i = 1,NTOT
            do j = 1,NTOT
                simmat(i,j) = -sum((datavecs(i,:) - datavecs(j,:))**2)
            end do
        end do
    end subroutine make_data

    ! the point of each true cluster with the largest summed similarity to its members
    function medoids( simmat, truth ) result( meds )
        real,    intent(in) :: simmat(NTOT,NTOT)
        integer, intent(in) :: truth(NTOT)
        integer :: meds(3), k, i
        real    :: best, score
        do k = 1,3
            best = -huge(best)
            do i = 1,NTOT
                if( truth(i) /= k ) cycle
                score = sum(simmat(i,:), mask=(truth == k))
                if( score > best )then
                    best    = score
                    meds(k) = i
                endif
            end do
        end do
    end function medoids

    subroutine test_three_clusters()
        real, parameter      :: PREFS(2) = [-1., -30.]
        type(aff_prop)       :: ap
        real                 :: simmat(NTOT,NTOT), simmat_in(NTOT,NTOT), simsum, simsum_ref
        integer              :: truth(NTOT), meds(3), i, j, ip, nerr
        integer, allocatable :: centers(:), labels(:)
        character(len=:), allocatable :: tag
        write(*,'(A)') 'test_three_clusters'
        call make_data(simmat, truth)
        simmat_in = simmat
        meds      = medoids(simmat, truth)
        do ip = 1,size(PREFS)
            tag = 'pref '//trim(real2str(PREFS(ip)))//': '
            call ap%new(NTOT, simmat, pref=PREFS(ip), lam=0.7, maxits=1000)
            call ap%propagate(centers, labels, simsum)
            call assert_int(3, size(centers), tag//'three clusters')
            if( size(centers) == 3 )then
                call assert_true(all(centers == meds), tag//'the exemplars are the cluster medoids')
                call assert_true(all(labels(centers) == [1,2,3]), tag//'each exemplar labels its own cluster')
            endif
            nerr = 0
            do i = 1,NTOT - 1
                do j = i + 1,NTOT
                    if( (truth(i) == truth(j)) .neqv. (labels(i) == labels(j)) ) nerr = nerr + 1
                end do
            end do
            call assert_int(0, nerr, tag//'the labels reproduce the true partition')
            ! simsum: the mean over all points of the similarity of each non-exemplar to its exemplar
            simsum_ref = 0.
            do j = 1,NTOT
                if( j /= centers(labels(j)) ) simsum_ref = simsum_ref + simmat(centers(labels(j)),j)
            end do
            simsum_ref = simsum_ref / real(NTOT)
            call assert_real(simsum_ref, simsum, 1.e-5, tag//'simsum is the mean similarity to the exemplars')
            call ap%kill
        end do
        call assert_true(all(simmat == simmat_in), 'new() leaves the input similarity matrix untouched')
    end subroutine test_three_clusters

    ! a preference at the largest similarity (0) makes every point its own exemplar; one far below
    ! every similarity leaves a single cluster
    subroutine test_preference_limits()
        type(aff_prop)       :: ap
        real                 :: simmat(NTOT,NTOT), simsum
        integer              :: truth(NTOT), i
        integer, allocatable :: centers(:), labels(:)
        write(*,'(A)') 'test_preference_limits'
        call make_data(simmat, truth)
        call ap%new(NTOT, simmat, pref=0., lam=0.7, maxits=1000)
        call ap%propagate(centers, labels, simsum)
        call assert_int(NTOT, size(centers), 'pref 0: every point is an exemplar')
        if( size(centers) == NTOT )then
            call assert_true(all(centers == [(i, i=1,NTOT)]) .and. all(labels == [(i, i=1,NTOT)]), 'pref 0: every point labels itself')
        endif
        call assert_real(0., simsum, 0., 'pref 0: simsum 0 without non-exemplars')
        call ap%kill
        call ap%new(NTOT, simmat, pref=-1000., lam=0.7, maxits=1000)
        call ap%propagate(centers, labels, simsum)
        call assert_int(1, size(centers), 'pref -1000: a single cluster')
        call assert_true(all(labels == 1), 'pref -1000: every point in it')
        call ap%kill
    end subroutine test_preference_limits

    ! AP restarts must be deterministic (volcluster, cluster_cavgs)
    subroutine test_determinism()
        type(aff_prop)       :: ap
        real                 :: simmat(NTOT,NTOT), simsum, simsum2
        integer              :: truth(NTOT), rep
        integer, allocatable :: centers(:), labels(:), centers2(:), labels2(:)
        logical              :: same
        write(*,'(A)') 'test_determinism'
        call make_data(simmat, truth)
        call ap%new(NTOT, simmat, pref=-1., lam=0.7, maxits=1000)
        call ap%propagate(centers, labels, simsum)
        same = .true.
        do rep = 1,3
            call ap%new(NTOT, simmat, pref=-1., lam=0.7, maxits=1000)
            call ap%propagate(centers2, labels2, simsum2)
            if( size(centers2) /= size(centers) )then
                same = .false.
            else
                if( any(centers2 /= centers) .or. any(labels2 /= labels) .or. simsum2 /= simsum ) same = .false.
            endif
        end do
        call assert_true(same, 'three restarts reproduce the exemplars, labels and simsum exactly')
        call ap%kill
    end subroutine test_determinism

end module simple_aff_prop_tester
