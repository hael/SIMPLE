! prng_uniform_plugin.f90 --
!     Module for uniform PRNGs
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
module prng_uniform_plugin
    use prng_class

    implicit none

    type, extends(prng) :: prng_uniform
        real :: xmin, xmax
    contains
        procedure :: get => get_uniform
    end type

contains
subroutine create_prngX( params, p1, p2 )
!DEC$ attributes dllexport :: create_prng

    class(prng_uniform), pointer :: params
    real                         :: p1
    real                         :: p2

    allocate( params )

    params%xmin = p1
    params%xmax = p2

end subroutine create_prngX

subroutine delete_prng( params )
!DEC$ attributes dllexport :: delete_prng

    class(prng_uniform), pointer :: params

    deallocate( params )

end subroutine delete_prng


real function get_uniform( params )

    class(prng_uniform) :: params

    real               :: r

    call random_number( r )

    get_uniform = params%xmin + (params%xmax-params%xmin) * r

end function get_uniform

end module prng_uniform_plugin

subroutine create_prng( params, p1, p2 )
!DEC$ attributes dllexport :: create_prng

    use prng_uniform_plugin

    class(prng_uniform), pointer :: params
    real                         :: p1
    real                         :: p2

    allocate( params )

    params%xmin = p1
    params%xmax = p2

end subroutine create_prng
