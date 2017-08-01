! newton.f90 --
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
!     Straightforward implementation of the Newton-Raphson method
!     for finding roots of an equation
!
module newton_raphson

    implicit none

contains

subroutine find_root( f, xinit, tol, maxiter, result, success )

    interface
        real function f(x)
            real, intent(in) :: x
        end function f
    end interface

    real, intent(in)     :: xinit
    real, intent(in)     :: tol
    integer, intent(in)  :: maxiter
    real, intent(out)    :: result
    logical, intent(out) :: success

    real                 :: eps = 1.0e-4
    real                 :: fx1
    real                 :: fx2
    real                 :: fprime
    real                 :: x
    real                 :: xnew
    integer              :: i

    result  = 0.0
    success = .false.

    x = xinit
    do i = 1,max(1,maxiter)
        fx1    = f(x)
        fx2    = f(x+eps)
        write(*,*) i, fx1, fx2, eps
        fprime = (fx2 - fx1) / eps

        xnew   = x - fx1 / fprime

        if ( abs(xnew-x) <= tol ) then
            success = .true.
            result  = xnew
            exit
        endif

        x = xnew
        write(*,*) i, x
     enddo

end subroutine find_root

end module

module functions
    implicit none

contains
real function f1( x )
    real, intent(in) :: x

    f1 = exp(x) + x  ! exp(x)-x is also interesting
end function

real function f2( x )
    real, intent(in) :: x

    f2 = x**2 + 1.0
end function

real function f3( x )
    real, intent(in) :: x

    ! f3 = sqrt(abs(x)) + 1.0  - cycle: 12.4 -- -130.08
    ! f3 = sqrt(abs(x)) + 0.1  - slow divergence
    ! f3 = sqrt(abs(x)) + 0.001  - cycle: 1.003 -- 1.004 -- 1.003 -- -1.004
    f3 = sqrt(abs(x)) + 0.001

end function

real function f4( x )
    real, intent(in) :: x

    f4 = abs(x) ** 0.3333
    if ( x > 0.0 ) then
       f4 = f4 - 1.0
    else
       f4 = -f4 - 1.0
    endif

end function

real function fln( x )
    real, intent(in) :: x

    fln = log(x) - 1.0

end function

real function fparabola( x )
    real, intent(in) :: x

    fparabola = x ** 2

end function

real function fparabola2( x )
    real, intent(in) :: x

    fparabola2 = x ** 2 + 1.0

end function

real function fsqrt( x )
    real, intent(in) :: x

    fsqrt = sqrt(abs(x))

end function


end module

program test_newton
    use newton_raphson
    use functions

    implicit none

    real           :: x
    real           :: root
    integer        :: maxiter = 20
    logical        :: success

    write(*,*) 'fln: 0.1'
    x = 0.1
    call find_root( fln, x, 1.0e-5, maxiter, root, success )
    write(*,*) root, fln(root), success

    write(*,*) 'fln: 10.0'
    x = 10.0
    call find_root( fln, x, 1.0e-5, maxiter, root, success )
    write(*,*) root, fln(root), success

    write(*,*) 'fparabola: 10.0'
    x = 10.0
    call find_root( fparabola, x, 1.0e-5, maxiter, root, success )
    write(*,*) root, fparabola(root), success

    write(*,*) 'fparabola2: 10.0'
    x = 10.0
    call find_root( fparabola2, x, 1.0e-5, maxiter, root, success )
    write(*,*) root, fparabola2(root), success

    write(*,*) 'fsqrt: 1.0'
    x = 1.0
    call find_root( fsqrt, x, 1.0e-5, maxiter, root, success )
    write(*,*) root, fsqrt(root), success

!
!    x = -2.0   ! x = 2.0 works, -2.0 does not
!    call find_root( f4, x, 1.0e-5, maxiter, root, success )
!
!    write(*,*) root, f4(root), success
!
end program test_newton
