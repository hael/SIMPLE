! compfunction.f90 --
!     Compute a simple function - Ftcl package
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
! functions --
!     Module for computing a simple function
!
module functions
    use ftcl

    implicit none
contains

! compute_func --
!     Compute the function:
!         f(x) = exp(-x) * cos(ax)
!
! Arguments:
!     cmdname         Name of the Tcl command
!     noargs          Number of arguments
!
subroutine compute_func( cmdname, noargs )
    character(len=*) :: cmdname
    integer          :: noargs

    real             :: param
    real             :: x
    real             :: result

    if ( noargs .ne. 2 ) then
        call ftcl_set_error( "Usage: " // trim(cmdname) // " param x" )
        return
    else
        call ftcl_get_arg( 1, param )
        call ftcl_get_arg( 2, x )

        result = exp(-x) * cos(param * x)

        call ftcl_put_result( result )
    endif
end subroutine compute_func
end module functions

! package_init --
!     Register the commands and the package itself
!
! Arguments:
!     error            Zero if all is okay, otherwise an error
!                      occurred while initialising
!
subroutine package_init( error )
    use ftcl
    use functions

    implicit none
    integer :: error

    error = 0

    call ftcl_make_command( compute_func, "func" )

    call ftcl_provide_package( "functions", "1.0", error )

end subroutine package_init
