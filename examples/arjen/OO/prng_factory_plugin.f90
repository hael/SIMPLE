! prng_factory_plugin.f90 --
!     Plugin handler for pseudo-random number
!     generators (PRNG)
!
!     Note:
!     For some reason the interface does not get resolved with
!     Intel Fortran 11.1 and I need to use get_procedure_global
!     instead
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
module prng_factory
    use prng_class
    use dynamic_libraries

    implicit none

    private
    public :: prng_create, prng

    type :: prng_creator
        character(len=20)                               :: name   = '?'
        procedure(creation_procedure), pointer, nopass  :: create => null()
    endtype prng_creator

    abstract interface
        subroutine creation_procedure( params, p1, p2 )

            import               :: prng
            class(prng), pointer :: params
            real                 :: p1
            real                 :: p2
        end subroutine creation_procedure
    end interface


    type(prng_creator), dimension(2), save :: prng_creators
    logical, save                          :: initialised = .false.

contains

subroutine initialise_factory

    character(len=20), dimension(2) :: prng_name
    character(len=20), dimension(2) :: libname
    type(dynamic_library)           :: dynlib
    logical                         :: success

    integer                         :: i

    prng_name = (/ 'uniform             ', 'exponential         ' /)
    libname   = (/ 'prng_uniform.dll    ', 'prng_exponential.dll' /)

    do i = 1,size(libname)
        call load_library( dynlib, libname(i), success )
        if ( success ) then
            call get_procedure_global( dynlib, 'create_prng', prng_creators(i)%create, success )
            if ( .not. success ) then
                write(*,*) 'Could not load create_prng - ', libname(i)
            endif
            prng_creators(i)%name = prng_name(i)
        endif
    enddo

end subroutine initialise_factory

function prng_create( type, param1, param2 )
!DEC$ attributes dllexport :: prng_create

     character(len=*)     :: type
     real                 :: param1
     real, optional       :: param2

     class(prng), pointer :: prng_create
     real                 :: param2opt
     integer              :: i

     if ( .not. initialised ) then
         initialised = .true.
         call initialise_factory
     endif

     param2opt = 0.0
     if ( present(param2) ) then
         param2opt = param2
     endif

     prng_create => null()

     do i = 1,size(prng_creators)
         if ( prng_creators(i)%name == type ) then
             call prng_creators(i)%create( prng_create, param1, param2opt )
             exit
         endif
     enddo
end function prng_create

end module prng_factory
