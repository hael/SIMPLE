! test_tridiag.f90 --
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
!     Experiment with test-driven development in a numerical
!     context: can we sensibly formulate tests before having
!     the implementation that needs to be tested?
!

! tridiag --
!     Simple module for solving systems of linear equations
!     of the form:
!         a(k) x(k-1)  +  b(k) x(k)  +   c(k) x(k+1)  =  d(k)
!
!     Note:
!     With the original routine (solve_original) all tests
!     fail - as required.
!
!     Note:
!     The subroutines change the arguments b and d
!
module tridiag
    implicit none

contains

subroutine solve_original( a, b, c, d, x )
    real, dimension(:) :: a, b, c, d, x

    x = 0.0 ! Just a stub for now

end subroutine solve_original

subroutine solve( a, b, c, d, x )
    real, dimension(:) :: a, b, c, d, x

    integer :: i
    integer :: n
    real    :: factor

    n = size(a)
    do i = 2,n
        factor = a(i) / b(i-1)
        b(i)   = b(i) - factor * c(i-1)
        d(i)   = d(i) - factor * d(i-1)
    enddo

    x(n) = d(n) / b(n)
    do i = n-1,1,-1
        x(i) = (d(i) - c(i) * x(i+1)) / b(i)
    enddo

end subroutine solve

end module tridiag

! test_tridiag --
!     Module containing all the tests for the tridiag module
!
!     Note: I found the need for the "boundary condition" in test_basic
!     only after checking the results for all tests. All answers were
!     as expected (modulo some truncation errors in the solution for the
!     largest system - test_diagonal_dom3) but this was not the case
!     for test_basic. Without the different value for d in the last
!     row, the solution is completely different than a set of ones.
!
!     This incident shows that it does pay off to formulate these
!     tests on the outset. It forces you to think about tests that
!     may or may not be relevant, but that at least can be used.
!
module test_tridiag
    use ftnunit
    use tridiag

    implicit none

    real, parameter :: margin = 1.0e-6

contains

subroutine test_all

    call test( test_trivial,       "Solve trivial system a=0, b=3, c=0, d=1" )
    call test( test_basic,         "Solve basic system a=0, b=6, c=-5, d=1" )
    call test( test_diagonal_dom1, "Solve diagonally dominant system - n=3")
    call test( test_diagonal_dom2, "Solve diagonally dominant system - n=10")
    call test( test_diagonal_dom3, "Solve diagonally dominant system - n=100")

end subroutine test_all

subroutine test_trivial

    integer, parameter :: rows = 10
    real, dimension(rows) :: a, b, c, d, x, y

    a = 0.0
    b = 3.0 ! Using 3 because of numerical rounding issues
    c = 0.0
    d = 1.0

    y = 1.0/3.0 ! Expected solution

    call solve( a, b, c, d, x )

    call assert_comparable( x, y, margin, "Solution is uniformly 1/3" )

end subroutine test_trivial

subroutine test_basic

    integer, parameter :: rows = 10
    real, dimension(rows) :: a, b, c, d, x, y

    a = 0.0
    b = 6.0
    c = -5.0
    d = 1.0

    d(rows) = 6.0 ! Boundary condition!

    !

    y = 1.0 ! Expected solution

    call solve( a, b, c, d, x )

    call assert_comparable( x, y, margin, "Solution is uniformly 1" )

end subroutine test_basic

subroutine test_diagonal_dom1

    integer, parameter :: rows = 3
    real, dimension(rows) :: a, b, c, d, x, y

    a = -1.0
    b = 2.0
    c = -1.0
    d(1) =              b(1) * 1.0 + c(1) * 1.0
    d(2) = a(2) * 1.0 + b(2) * 1.0 + c(2) * 1.0
    d(3) = a(3) * 1.0 + b(3) * 1.0

    y = 1.0 ! Expected solution

    call solve( a, b, c, d, x )

    call assert_comparable( x, y, margin, "Solution is uniformly 1" )

end subroutine test_diagonal_dom1

subroutine test_diagonal_dom2

    integer, parameter    :: rows = 10
    real, dimension(rows) :: a, b, c, d, x, y
    integer               :: i


    a = -1.0
    b = 2.0
    c = -1.0
    y = (/ (1.0/i ,i=1,rows) /) ! Expected solution

    d(2:rows-1) = a(2:rows-1) * y(1:rows-2) + b(2:rows-1) * y(2:rows-1) + &
                  c(2:rows-1) * y(3:rows)
    d(1)        =                             b(1) * y(1) + c(1) * y(2)
    d(rows)     = a(rows) * y(rows-1)        + b(rows) * y(rows)

    call solve( a, b, c, d, x )

    call assert_comparable( x, y, margin, "Solution is 1/k (10 rows)" )

end subroutine test_diagonal_dom2

subroutine test_diagonal_dom3

    integer, parameter    :: rows = 100
    real, dimension(rows) :: a, b, c, d, x, y
    integer               :: i

    a = -1.0
    b = 2.0
    c = -1.0
    y = (/ (1.0/i ,i=1,rows) /) ! Expected solution

    d(2:rows-1) = a(2:rows-1) * y(1:rows-2) + b(2:rows-1) * y(2:rows-1) + &
                  c(2:rows-1) * y(3:rows)
    d(1)        =                             b(1) * y(1) + c(1) * y(2)
    d(rows)     = a(rows) * y(rows-1)        + b(rows) * y(rows)

    call solve( a, b, c, d, x )

    call assert_comparable( x, y, margin, "Solution is 1/k (100 rows)" )

end subroutine test_diagonal_dom3

end module test_tridiag

! run_tridiag_tests
!     Program to run the tests
!
program run_tridiag_tests
    use ftnunit
    use test_tridiag

    implicit none

    call runtests_init
    call runtests( test_all )
    call runtests_final

end program run_tridiag_tests
