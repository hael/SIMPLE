! linux_dynlib --
!     Implementation for Linux and OS/X of low-level routines
!     that deal with dynamic libraries
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
module system_dynamic_libraries
    use iso_c_binding
    implicit none

    integer(kind=c_int), parameter :: rtld_lazy = 1

    interface
        function c_load_library( cname, load_type ) bind(c, name='dlopen' )
            use iso_c_binding
            character(kind=c_char), dimension(*) :: cname
            integer(kind=c_int), value           :: load_type
            integer(kind=c_long)                 :: c_load_library
        end function c_load_library
    end interface

    interface
        function c_get_procedure( handle, cname ) bind(c, name='dlsym' )
            use iso_c_binding
            integer(kind=c_long), value          :: handle
            character(kind=c_char), dimension(*) :: cname
            type(c_funptr)                       :: c_get_procedure
        end function c_get_procedure
    end interface

    interface
        subroutine c_print_error( ) bind(c, name='print_error' )
            use iso_c_binding
        end subroutine c_print_error
    end interface


contains

! system_load_library --
!     Load the library
!
! Arguments:
!     handle         Handle to the library
!     cname          Null-terminated name of the library
!
! Returns:
!     Handle to the library
!
subroutine system_load_library( handle, cname )
    integer(kind=c_long)           :: handle
    character(len=1), dimension(*) :: cname

    integer(kind=c_int), parameter :: load_type = rtld_lazy

    handle = c_load_library( cname, load_type )
    call c_print_error
end subroutine system_load_library

! system_get_procedure --
!     Get the procedure
!
! Arguments:
!     handle         Handle to the library
!     cname          Null-terminated name of the procedure
!     cproc          C-style procedure pointer
!     success        Whether successful or not
!
! Returns:
!     Handle to the library
!
subroutine system_get_procedure( handle, cname, cproc, success )
    integer(kind=c_long), value    :: handle
    character(len=1), dimension(*) :: cname
    type(c_funptr)                 :: cproc
    logical                        :: success

    call c_print_error
    cproc = c_get_procedure( handle, cname )
    call c_print_error

    success = c_associated(cproc)

end subroutine system_get_procedure

end module system_dynamic_libraries
