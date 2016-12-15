! derivative_near_singular.f90 --
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
!     Simple program to determine the performance differences between
!     using array-valued functions and ordinary do-loops
!
!     Determine (automatically and manually) the derivative near a
!     mildly singular point
!
program derivative_near_singular

    use automatic_differentiation

    implicit none

    type(autoderiv) :: x
    type(autoderiv) :: result
    real(kind=wp)   :: approx
    integer         :: i

    real(kind=wp), dimension(10) :: value = &
        (/ 0.001, 0.0005, 0.0002, 0.0001, 0.00005, 0.00002, 0.00001, &
           0.000005, 0.000002, 0.000001 /)

    do i = 1,size(value)
        x%v  = 1.0_wp + value(i)
        x%dv = 1.0_wp

        result = f(x)
        approx = fprime(x%v)

        write(*,'(10f12.8)') x%v, result%v, result%dv, approx, &
            result%dv+0.5_wp+2.0*(x%v-1.0_wp)/3.0_wp, &
            approx   +0.5_wp+2.0*(x%v-1.0_wp)/3.0_wp
    enddo

contains

function f( x )
    type(autoderiv), intent(in) :: x
    type(autoderiv)             :: f

    f = log(x) / (x-1.0_wp)
end function f

real(kind=wp) function fprime( x )
    real(kind=wp) :: x

    real(kind=wp) :: eps

    eps    = x - 1.0_wp
    fprime = (-0.5_wp + eps / 6.0_wp - eps**2 / 12.0_wp + eps**3 / 20.0_wp) / x
end function fprime

end program derivative_near_singular
