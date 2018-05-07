program calc_pi_serial
implicit none

integer :: i, n
double precision, parameter :: pi25dt = 3.141592653589793238462643d0
double precision :: a,h,pi,sum,x

do
    print *, 'Enter the number of intervals: (0 quits)'
    read *, n
    if( n <= 0 ) exit
    h = 1.d0 / dble(n) ! stride
    ! calculation of the quadrature (summation)
    sum = 0.d0
    do i=1,n
        x = h * (dble(i) - 0.5d0)
        sum = sum + f(x)
    end do
    pi = h * sum
    print *, 'pi is approximately: ', pi, ' Error is: ', abs(pi-pi25dt)

end do

contains

    double precision function f( a )
        double precision, intent(in) :: a
        f = 4.d0 / (1.d0 + a*a)
    end function f

end program calc_pi_serial
