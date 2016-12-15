! prng_exponential_plugin.f90 --
!     Module for exponential PRNGs
!
!     Note: keep the creation routine, create_prng, outside the module
!     That makes retrieving the procedure pointer easier
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
module prng_exponential_plugin
    use prng_class

    implicit none

    type, extends(prng) :: prng_exponential
        real :: xmean
    contains
        procedure :: get => get_exponential
    end type

contains

subroutine create_prngX( params, p1, p2 )
!DEC$ attributes dllexport :: create_prng

    class(prng_exponential), pointer :: params
    real                             :: p1
    real                             :: p2

    allocate( params )

    params%xmean = p1

end subroutine create_prngX

subroutine delete_prng( params )
!DEC$ attributes dllexport :: delete_prng

    class(prng_exponential), pointer :: params

    deallocate( params )

end subroutine delete_prng


real function get_exponential( params )
!DEC$ attributes dllexport :: get_exponential

    class(prng_exponential) :: params

    real               :: r

    call random_number( r )

    get_exponential = -params%xmean * log( r )

end function get_exponential

end module prng_exponential_plugin

subroutine create_prng( params, p1, p2 )
!DEC$ attributes dllexport :: create_prng
    use prng_exponential_plugin

    class(prng_exponential), pointer :: params
    real                             :: p1
    real                             :: p2

    allocate( params )

    params%xmean = p1

end subroutine create_prng
