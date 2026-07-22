program simple_test_rnd_shuffle
use simple_rnd,        only: seed_rnd, shuffle, partial_shuffle
use simple_math,       only: hpsort
use simple_test_utils, only: assert_true, report_summary, tests_failed
implicit none

integer, parameter :: N = 100, NPARTIAL = 17
integer            :: iarr(N), iarr_copy(N), expected(N), i
real               :: rarr(N), rarr_copy(N), rexpected(N)

call seed_rnd
expected  = [(i,i=1,N)]
rexpected = real(expected)

iarr = expected
call shuffle(iarr)
iarr_copy = iarr
call hpsort(iarr_copy)
call assert_true(all(iarr_copy == expected), 'integer full shuffle preserves the input permutation')

rarr = rexpected
call shuffle(rarr)
rarr_copy = rarr
call hpsort(rarr_copy)
call assert_true(all(rarr_copy == rexpected), 'real full shuffle preserves the input permutation')

! A zero-length partial shuffle is a no-op.
iarr = expected
call partial_shuffle(iarr, 0)
call assert_true(all(iarr == expected), 'zero-length integer partial shuffle is a no-op')

! A partial shuffle selects an ordered, duplicate-free prefix while retaining
! every original value in the complete array.
call partial_shuffle(iarr, NPARTIAL)
iarr_copy = iarr
call hpsort(iarr_copy)
call assert_true(all(iarr_copy == expected), 'integer partial shuffle preserves the input permutation')
iarr_copy(:NPARTIAL) = iarr(:NPARTIAL)
call hpsort(iarr_copy(:NPARTIAL))
call assert_true(all(iarr_copy(2:NPARTIAL) /= iarr_copy(1:NPARTIAL-1)), &
    &'integer partial shuffle prefix contains no duplicates')

rarr = rexpected
call partial_shuffle(rarr, NPARTIAL)
rarr_copy = rarr
call hpsort(rarr_copy)
call assert_true(all(rarr_copy == rexpected), 'real partial shuffle preserves the input permutation')

! Selecting the complete prefix is a full Fisher-Yates permutation.
iarr = expected
call partial_shuffle(iarr, N)
iarr_copy = iarr
call hpsort(iarr_copy)
call assert_true(all(iarr_copy == expected), 'complete partial shuffle preserves the input permutation')

call report_summary()
if( tests_failed > 0 ) error stop 1
end program simple_test_rnd_shuffle
