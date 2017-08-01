! prng_use_plugin.f90 --
!     Example of a simple "factory" for pseudo-random number
!     generators (PRNG):
!     The main program using the DLL
!
!     Note:
!     Use strings instead of parameters. If the library is extended,
!     we do not know what new parameters (actually integer values)
!     are available. With strings it is more "natural" to incorporate
!     the literal values.
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
program prng_use_plugin
    use prng_factory

    class(prng), pointer :: p
    integer              :: i

    p => prng_create( "uniform", 1.0, 2.0 )

    do i = 1,10
        write(*,*) p%get()
    enddo

    p => prng_create( "exponential", 10.0 )

    do i = 1,10
        write(*,*) p%get()
    enddo

end program prng_use_plugin
