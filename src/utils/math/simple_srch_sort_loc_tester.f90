!@descr: unit test routines for simple_srch_sort_loc (sorting, searching, locating, selection, ordering)
! Heap sorts in all their forms, table search (locate, find), n-largest/smallest locators, peak finding,
! reversal, selection, unique, orderings and reordering, each against a brute-force reference on small
! fixed arrays with ties and duplicates.
module simple_srch_sort_loc_tester
use simple_test_utils    ! assertions etc.
use simple_defs          ! dp
use simple_string_utils, only: int2str
use simple_srch_sort_loc
implicit none
private
public :: run_all_srch_sort_loc_tests

integer, parameter :: NARR = 12
! fixed unsorted data with a duplicate (2.5 twice) and a tie in the integers (7 twice)
real,    parameter :: RVALS(NARR) = [3.5, -1.0, 2.5, 9.25, 0.0, 2.5, -4.75, 6.0, 1.5, 8.0, -2.0, 7.25]
integer, parameter :: IVALS(NARR) = [7, -3, 12, 0, 7, 5, -8, 19, 2, 11, -1, 4]

contains

    subroutine run_all_srch_sort_loc_tests()
        write(*,'(A)') '**** running all search/sort/locate tests ****'
        call test_hpsort_reals()
        call test_hpsort_integers()
        call test_hpsort_with_comparator()
        call test_hpsort_edge_sizes()
        call test_locate()
        call test_find()
        call test_maxnloc_minnloc()
        call test_min3()
        call test_peakfinder()
        call test_reverse()
        call test_selec()
        call test_unique()
        call test_orderings()
        call test_mask2inds_and_reorder()
    end subroutine run_all_srch_sort_loc_tests

    !---------------- references ----------------

    ! insertion sort, ascending, as the reference
    pure function sorted_reals( arr ) result( out )
        real, intent(in) :: arr(:)
        real    :: out(size(arr)), tmp
        integer :: i, j
        out = arr
        do i = 2,size(out)
            tmp = out(i)
            j   = i - 1
            do while( j >= 1 )
                if( out(j) <= tmp ) exit
                out(j+1) = out(j)
                j = j - 1
            enddo
            out(j+1) = tmp
        enddo
    end function sorted_reals

    pure function sorted_ints( arr ) result( out )
        integer, intent(in) :: arr(:)
        integer :: out(size(arr)), tmp, i, j
        out = arr
        do i = 2,size(out)
            tmp = out(i)
            j   = i - 1
            do while( j >= 1 )
                if( out(j) <= tmp ) exit
                out(j+1) = out(j)
                j = j - 1
            enddo
            out(j+1) = tmp
        enddo
    end function sorted_ints

    pure logical function is_permutation( inds, n )
        integer, intent(in) :: inds(:), n
        logical :: seen(n)
        integer :: i
        seen = .false.
        is_permutation = size(inds) == n
        if( .not. is_permutation ) return
        do i = 1,n
            if( inds(i) < 1 .or. inds(i) > n )then
                is_permutation = .false.
                return
            endif
            if( seen(inds(i)) )then
                is_permutation = .false.
                return
            endif
            seen(inds(i)) = .true.
        enddo
    end function is_permutation

    function descending( p1, p2 ) result( val )
        integer, intent(in) :: p1, p2
        logical :: val
        val = p1 > p2
    end function descending

    !---------------- heap sorts ----------------

    subroutine test_hpsort_reals()
        real    :: arr(NARR), arr2(NARR), ref(NARR)
        integer :: inds(NARR), i
        logical :: consistent
        write(*,'(A)') 'test_hpsort_reals'
        ref = sorted_reals(RVALS)
        ! hpsort(rarr)
        arr = RVALS
        call hpsort(arr)
        call assert_true(all(arr == ref), 'hpsort(rarr): ascending, duplicates kept')
        ! hpsort(rarr, iarr): the index array follows the values
        arr  = RVALS
        inds = [(i, i=1,NARR)]
        call hpsort(arr, inds)
        call assert_true(all(arr == ref), 'hpsort(rarr, iarr): values ascending')
        call assert_true(is_permutation(inds, NARR), 'hpsort(rarr, iarr): index array is a permutation')
        consistent = .true.
        do i = 1,NARR
            if( RVALS(inds(i)) /= arr(i) ) consistent = .false.
        enddo
        call assert_true(consistent, 'hpsort(rarr, iarr): arr(i) == original(iarr(i))')
        ! hpsort(rarr, rarr2): the second real array follows the first
        arr  = RVALS
        arr2 = 10. * RVALS
        call hpsort(arr, arr2)
        call assert_true(all(arr == ref), 'hpsort(rarr, rarr2): values ascending')
        call assert_true(maxval(abs(arr2 - 10. * arr)) < 1.e-6, 'hpsort(rarr, rarr2): the carried array moved with the values')
    end subroutine test_hpsort_reals

    subroutine test_hpsort_integers()
        integer :: arr(NARR), inds(NARR), ref(NARR), i
        logical :: consistent
        write(*,'(A)') 'test_hpsort_integers'
        ref = sorted_ints(IVALS)
        arr = IVALS
        call hpsort(arr)
        call assert_true(all(arr == ref), 'hpsort(iarr): ascending with the tie kept')
        arr  = IVALS
        inds = [(i, i=1,NARR)]
        call hpsort(arr, inds)
        call assert_true(all(arr == ref), 'hpsort(iarr, inds): values ascending')
        call assert_true(is_permutation(inds, NARR), 'hpsort(iarr, inds): index array is a permutation')
        consistent = .true.
        do i = 1,NARR
            if( IVALS(inds(i)) /= arr(i) ) consistent = .false.
        enddo
        call assert_true(consistent, 'hpsort(iarr, inds): arr(i) == original(inds(i))')
    end subroutine test_hpsort_integers

    subroutine test_hpsort_with_comparator()
        integer :: arr(NARR), ref(NARR), i
        write(*,'(A)') 'test_hpsort_with_comparator'
        ref = sorted_ints(IVALS)
        arr = IVALS
        call hpsort(arr, descending)
        call assert_true(all(arr == [(ref(NARR+1-i), i=1,NARR)]), 'hpsort(iarr, p1_lt_p2): a "greater than" comparator sorts descending')
    end subroutine test_hpsort_with_comparator

    subroutine test_hpsort_edge_sizes()
        real    :: one(1), two(2), same(5)
        integer :: ione(1), inds(2)
        write(*,'(A)') 'test_hpsort_edge_sizes'
        one = [4.]
        call hpsort(one)
        call assert_real(4., one(1), 0., 'hpsort: a single element is untouched')
        ione = [4]
        call hpsort(ione)
        call assert_int(4, ione(1), 'hpsort: a single integer is untouched')
        two  = [2., 1.]
        inds = [1, 2]
        call hpsort(two, inds)
        call assert_true(all(two == [1., 2.]) .and. all(inds == [2, 1]), 'hpsort: two elements swap with their indices')
        same = 3.
        call hpsort(same)
        call assert_true(all(same == 3.), 'hpsort: all-equal input is unchanged')
    end subroutine test_hpsort_edge_sizes

    !---------------- table search ----------------

    subroutine test_locate()
        real,    parameter :: TAB(6)  = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0]
        integer, parameter :: ITAB(6) = [1, 2, 4, 8, 16, 32]
        real,    parameter :: DTAB(4) = [10., 5., 2., 1.] ! decreasing tables are allowed
        write(*,'(A)') 'test_locate'
        call assert_int(3, locate(TAB, 6, 5.0),   'locate: 5 lies between arr(3) = 4 and arr(4) = 8')
        call assert_int(1, locate(TAB, 6, 1.5),   'locate: 1.5 lies in the first interval')
        call assert_int(5, locate(TAB, 6, 20.0),  'locate: 20 lies in the last interval')
        call assert_int(0, locate(TAB, 6, 0.5),   'locate: below the table gives 0')
        call assert_int(6, locate(TAB, 6, 40.0),  'locate: above the table gives n')
        call assert_int(1, locate(TAB, 6, 1.0),   'locate: exactly the first entry gives 1')
        call assert_int(5, locate(TAB, 6, 32.0),  'locate: exactly the last entry gives n-1')
        call assert_int(3, locate(TAB, 6, 4.0),   'locate: exactly an interior entry gives its interval')
        call assert_int(2, locate(DTAB, 4, 3.0),  'locate: decreasing table, 3 lies between 5 and 2')
        call assert_int(3, locate(ITAB, 6, 5),    'locate: integer table')
        call assert_int(0, locate(ITAB, 6, -3),   'locate: integer below the table gives 0')
        call assert_int(1, locate(ITAB, 6, 1),    'locate: integer first entry gives 1')
        call assert_int(5, locate(ITAB, 6, 32),   'locate: integer last entry gives n-1')
    end subroutine test_locate

    subroutine test_find()
        real,    parameter :: TAB(6)  = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0]
        integer, parameter :: ITAB(6) = [1, 2, 4, 8, 16, 32]
        integer :: j, idist
        real    :: dist
        write(*,'(A)') 'test_find'
        call find(TAB, 6, 4.9, j, dist)
        call assert_int(3, j,           'find: 4.9 is closest to 4')
        call assert_real(-0.9, dist, 1.e-6, 'find: signed distance arr(j) - x')
        call find(TAB, 6, 7.0, j, dist)
        call assert_int(4, j,           'find: 7 is closest to 8')
        call assert_real(1.0, dist, 1.e-6, 'find: distance to 8')
        call find(TAB, 6, 6.0, j, dist)
        call assert_int(4, j,           'find: an exact tie (6 between 4 and 8) goes to the upper neighbour')
        call find(TAB, 6, 0.1, j, dist)
        call assert_int(1, j,           'find: below the table gives the first entry')
        call find(TAB, 6, 100., j, dist)
        call assert_int(6, j,           'find: above the table gives the last entry')
        call find(ITAB, 6, 13, j, idist)
        call assert_int(5, j,           'find: integer 13 is closest to 16')
        call assert_int(3, idist,       'find: integer distance arr(j) - x')
    end subroutine test_find

    !---------------- n largest / smallest ----------------

    subroutine test_maxnloc_minnloc()
        integer, parameter :: NSEL = 4
        real    :: ref(NARR), shuffled(30)
        integer :: loc(NSEL), i, loc1(1), locall(NARR), loc10(10)
        logical :: descending_vals, ascending_vals
        write(*,'(A)') 'test_maxnloc_minnloc'
        ref = sorted_reals(RVALS)
        loc = maxnloc(RVALS, NSEL)
        call assert_true(all(RVALS(loc) == [(ref(NARR+1-i), i=1,NSEL)]), 'maxnloc: the n largest values, largest first')
        call assert_true(all(loc >= 1 .and. loc <= NARR), 'maxnloc: indices in range')
        call assert_int(4, loc(1), 'maxnloc: the largest value (9.25) is at index 4')
        loc = minnloc(RVALS, NSEL)
        call assert_true(all(RVALS(loc) == ref(1:NSEL)), 'minnloc: the n smallest values, smallest first')
        call assert_int(7, loc(1), 'minnloc: the smallest value (-4.75) is at index 7')
        ! n equal to the array size: a full ordering
        locall = maxnloc(RVALS, NARR)
        call assert_true(is_permutation(locall, NARR), 'maxnloc with n = size: a permutation')
        call assert_true(all(RVALS(locall) == [(ref(NARR+1-i), i=1,NARR)]), 'maxnloc with n = size: full descending order')
        locall = minnloc(RVALS, NARR)
        call assert_true(all(RVALS(locall) == ref), 'minnloc with n = size: full ascending order')
        ! n = 1
        loc1 = maxnloc(RVALS, 1)
        call assert_int(4, loc1(1), 'maxnloc with n = 1 is maxloc')
        loc1 = minnloc(RVALS, 1)
        call assert_int(7, loc1(1), 'minnloc with n = 1 is minloc')
        ! the old test: a shuffled 1..30, the ten largest are 30 down to 21
        do i = 1,30
            shuffled(i) = real(mod(7 * i, 31)) ! 7 is coprime to 31: a permutation of 1..30
        enddo
        descending_vals = .true.
        ascending_vals  = .true.
        loc10 = maxnloc(shuffled, 10)
        do i = 1,10
            if( nint(shuffled(loc10(i))) /= 31 - i ) descending_vals = .false.
        enddo
        loc10 = minnloc(shuffled, 10)
        do i = 1,10
            if( nint(shuffled(loc10(i))) /= i ) ascending_vals = .false.
        enddo
        call assert_true(descending_vals, 'maxnloc: ten largest of a shuffled 1..30 are 30..21')
        call assert_true(ascending_vals,  'minnloc: ten smallest of a shuffled 1..30 are 1..10')
    end subroutine test_maxnloc_minnloc

    subroutine test_min3()
        real :: m3(3), ref(NARR)
        write(*,'(A)') 'test_min3'
        ref = sorted_reals(RVALS)
        m3  = min3(RVALS)
        call assert_true(all(m3 == ref(1:3)), 'min3: the three smallest, ascending')
        m3 = min3([5., 1., 3.])
        call assert_true(all(m3 == [5., 1., 3.]), 'min3: three or fewer elements are returned as they are')
        m3 = min3([2., 2., 2., 1.])
        call assert_true(all(m3 == [1., 2., 2.]), 'min3: a later smaller value displaces the ties')
    end subroutine test_min3

    !---------------- peaks, reversal, selection ----------------

    subroutine test_peakfinder()
        logical, allocatable :: peaks(:)
        write(*,'(A)') 'test_peakfinder'
        peaks = peakfinder([1., 3., 2., 2., 5., 4., 6.])
        call assert_true(all(peaks .eqv. [.false., .true., .false., .false., .true., .false., .true.]), &
            &'peakfinder: interior local maxima and a rising end are peaks')
        peaks = peakfinder([4., 1., 2., 1.])
        call assert_true(all(peaks .eqv. [.true., .false., .true., .false.]), 'peakfinder: a falling start is a peak')
        peaks = peakfinder([1., 2., 2., 1.])
        call assert_true(all(peaks .eqv. [.false., .true., .true., .false.]), 'peakfinder: a plateau counts on both sides (>=)')
    end subroutine test_peakfinder

    subroutine test_reverse()
        integer  :: iarr(5), ieven(4)
        real     :: rarr(5), reven(4)
        real(dp) :: darr(3), deven(4)
        write(*,'(A)') 'test_reverse'
        iarr = [1, 2, 3, 4, 5]
        call reverse(iarr)
        call assert_true(all(iarr == [5, 4, 3, 2, 1]), 'reverse: odd-length integer array')
        ieven = [1, 2, 3, 4]
        call reverse(ieven)
        call assert_true(all(ieven == [4, 3, 2, 1]), 'reverse: even-length integer array')
        rarr = [1., 2., 3., 4., 5.]
        call reverse(rarr)
        call assert_true(all(rarr == [5., 4., 3., 2., 1.]), 'reverse: odd-length real array')
        reven = [1., 2., 3., 4.]
        call reverse(reven)
        call assert_true(all(reven == [4., 3., 2., 1.]), 'reverse: even-length real array')
        darr = [1.d0, 2.d0, 3.d0]
        call reverse(darr)
        call assert_true(all(darr == [3.d0, 2.d0, 1.d0]), 'reverse: odd-length double precision array')
        deven = [1.d0, 2.d0, 3.d0, 4.d0]
        call reverse(deven)
        call assert_true(all(deven == [4.d0, 3.d0, 2.d0, 1.d0]), 'reverse: even-length double precision array')
        ! reverse_f mirrors about the Fourier origin at element 1: an even length keeps element 1 in place
        rarr = [1., 2., 3., 4., 5.]
        call reverse_f(rarr)
        call assert_true(all(rarr == [5., 4., 3., 2., 1.]), 'reverse_f: odd length is a plain reversal')
        reven = [1., 2., 3., 4.]
        call reverse_f(reven)
        call assert_true(all(reven == [1., 4., 3., 2.]), 'reverse_f: even length keeps element 1 (the origin) and reverses the rest')
    end subroutine test_reverse

    subroutine test_selec()
        real    :: arr(NARR), ref(NARR)
        integer :: k
        logical :: all_ok
        write(*,'(A)') 'test_selec'
        ref = sorted_reals(RVALS)
        all_ok = .true.
        do k = 1,NARR
            arr = RVALS ! selec reorders its input
            if( selec(k, NARR, arr) /= ref(k) ) all_ok = .false.
        enddo
        call assert_true(all_ok, 'selec(k, n, arr) is the k-th smallest value for every k (the header comment says largest)')
    end subroutine test_selec

    !---------------- unique, orderings, reordering ----------------

    subroutine test_unique()
        integer, allocatable :: u(:)
        write(*,'(A)') 'test_unique'
        call unique(IVALS, u)
        call assert_int(NARR - 1, size(u), 'unique: one duplicate removed')
        call assert_true(all(u == [-8, -3, -1, 0, 2, 4, 5, 7, 11, 12, 19]), 'unique: sorted distinct values')
        call unique([3, 3, 3], u)
        call assert_true(size(u) == 1 .and. u(1) == 3, 'unique: all-equal input gives one value')
    end subroutine test_unique

    subroutine test_orderings()
        integer, allocatable :: order(:)
        real    :: ref(NARR)
        integer :: i
        write(*,'(A)') 'test_orderings'
        ref   = sorted_reals(RVALS)
        order = scores2order(RVALS)
        call assert_true(is_permutation(order, NARR), 'scores2order: a permutation')
        call assert_true(all(RVALS(order) == [(ref(NARR+1-i), i=1,NARR)]), 'scores2order: best score first')
        order = dists2order(RVALS)
        call assert_true(is_permutation(order, NARR), 'dists2order: a permutation')
        call assert_true(all(RVALS(order) == ref), 'dists2order: shortest distance first')
    end subroutine test_orderings

    subroutine test_mask2inds_and_reorder()
        integer, allocatable :: inds(:)
        real    :: rarr(5)
        integer :: iarr(5)
        write(*,'(A)') 'test_mask2inds_and_reorder'
        inds = mask2inds([.true., .false., .false., .true., .true.])
        call assert_true(all(inds == [1, 4, 5]), 'mask2inds: indices of the true entries')
        inds = mask2inds([.false., .false.])
        call assert_int(0, size(inds), 'mask2inds: no true entries gives an empty array')
        rarr = [10., 20., 30., 40., 50.]
        call reorder(rarr, [5, 3, 1, 2, 4])
        call assert_true(all(rarr == [50., 30., 10., 20., 40.]), 'reorder: real array permuted by the order')
        iarr = [1, 2, 3, 4, 5]
        call reorder(iarr, [2, 1, 4, 3, 5])
        call assert_true(all(iarr == [2, 1, 4, 3, 5]), 'reorder: integer array permuted by the order')
    end subroutine test_mask2inds_and_reorder

end module simple_srch_sort_loc_tester
