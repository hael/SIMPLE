! explicit_vs_assumed.f90 --
!     Mismatch in interfaces (explicit-shape versus assumed-shape)
!     leading to unexpected output
!
!     Adapted from a posting by Clive Page on comp.lang.fortran
!     (9 may 2011)
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
module mymod
    implicit none
contains

subroutine mysub(param, result)
    real, intent(in)  :: param(3)
    real, intent(out) :: result
    print *,'param=', param
    result = 0.0
end subroutine mysub

subroutine minim(param, subr, result)
    real, intent(in) :: param(:)
    interface
        subroutine subr(p, r)
            real, intent(in)  :: p(:)
            real, intent(out) :: r
        end subroutine subr
    end interface
    real, intent(out):: result

    call subr(param, result)
end subroutine minim
end module mymod

program main
    use mymod
    implicit none

    real :: param(3) = [1.0, 2.0, 3.0], result

    call minim(param, mysub, result)
end program main
