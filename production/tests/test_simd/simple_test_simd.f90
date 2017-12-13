program simple_test_simd
use simple_timer
implicit none
integer, parameter :: N=10000000
real :: a(N), b(N), c(N), t1, t2
integer(timer_int_kind) :: t_loop, t_loop_simd
real(timer_int_kind)    :: rt_loop, rt_loop_simd


integer :: i
a = 0.
b = 0.
c = 0.

! t_loop = tic()
! do i=1,N
!     a(i) = b(i) + c(i)
! end do
! rt_loop = toc(t_loop)
! print *, 'time(loop): ', rt_loop

! t_loop_simd = tic()
! !$omp simd
! do i=1,N
!     a(i) = b(i) + c(i)
! end do
! !$omp end simd
! rt_loop_simd = toc(t_loop_simd)
! print *, 'time(loop_simd): ', rt_loop_simd
! print *, 'speedup with simd: ', rt_loop / rt_loop_simd

t_loop = tic()
do i=1,N
    t1 = func1(b(i), c(i))
    t2 = func2(b(i), c(i))
    a(i) = t1 + t2
end do
rt_loop = toc(t_loop)
print *, 'time(loop): ', rt_loop

t_loop_simd = tic()
!$omp simd private(t1,t2)
do i=1,N
    t1 = func1(b(i), c(i))
    t2 = func2(b(i), c(i))
    a(i) = t1 + t2
end do
!$omp end simd
rt_loop_simd = toc(t_loop_simd)
print *, 'time(loop_simd): ', rt_loop_simd
print *, 'speedup with simd: ', rt_loop / rt_loop_simd

contains

    real function func1( a, b )
        real, intent(in) :: a, b
        func1 = a * b
	end function

    real function func2( a, b )
        real, intent(in) :: a, b
        func2 = (a + b)**2.0
    end function

end program simple_test_simd
