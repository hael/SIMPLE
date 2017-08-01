program simd_example1
!$ use omp_lib
!$ use omp_lib_kinds
implicit none
integer, parameter :: N=32
integer :: i
double precision :: a(N), b(N)
do i=1,N
    a(i) = i-1
    b(i) = N-(i-1)
end do 
call work(a,b,N)
do i=1,N
    print *, i, a(i)
end do
end program simd_example1

! The uniform(fact) clause indicates that the variable fact 
! is invariant across the SIMD lanes
function add1(a,b,fact) result(c)
!$omp declare simd(add1) uniform(fact)
    double precision :: a,b,fact, c
    c = a + b + fact
end function add1

! Here, a and b are included in the unform list because the 
! Fortran array references are constant. The i index used in 
! the add2 function is included in a linear clause with a 
! constant-linear-step of 1, to guarantee a unity increment 
! of the associated loop.
function add2(a,b,i, fact) result(c)
!$omp declare simd(add2) uniform(a,b,fact) linear(i:1)
    integer :: i
    double precision :: a(*),b(*),fact, c
    c = a(i) + b(i) + fact
end function add2

subroutine work(a, b, n )
    double precision :: a(n),b(n), tmp
    integer :: n, i
    double precision, external :: add1, add2
    !$omp simd private(tmp)
    do i = 1,n
        tmp = add1(a(i), b(i), 1.0d0)
        a(i) = add2(a, b, i, 1.0d0) + tmp
        a(i) = a(i) + b(i) + 1.0d0
    end do
end subroutine work
