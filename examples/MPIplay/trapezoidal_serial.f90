! calculate definite integral using trapezoidal rule.
! The function f(x) is hardwired.
! Input: a, b, n.
! Output: estimate of integral from a to b of f(x) using n trapezoids.
program trapezoidal_serial
implicit none

real    :: integral, a, b, h, x
integer :: n, i

print *, 'Enter a, b, and n'
read *, a, b, n

h = (b-a)/real(n)
integral = (f(a) + f(b)) / 2.0
x = a
do i=1,n-1
    x = x + h
    integral = integral + f(x)
end do
integral = integral * h
print *, 'With n = ', n, ' trapezoids, our estimate of the integral from ', a, ' to ', b, ' = ', integral


contains

    real function f(x)
        real, intent(in) :: x
        f = x*x
    end function f


end program trapezoidal_serial
