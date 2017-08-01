! win_dynlib --
!     Implementation for (GNU Fortran on) Windows of low-level routines
!     that deal with dynamic libraries
!
!     Note:
!     You need gfortran 4.5.0 or later - version 4.4.x gave a preoblem
!     with the "LoadLibrary" function.
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

    interface
        function c_load_library( cname ) bind(c, name='LoadLibraryA')
            use iso_c_binding
            !GCC$ ATTRIBUTES STDCALL :: c_load_library
            character(kind=c_char), dimension(*) :: cname
            integer(kind=c_intptr_t)             :: c_load_library
        end function c_load_library
    end interface

    interface
        function c_get_procedure( handle, cname ) bind(c, name='GetProcAddress')
            use iso_c_binding
            !GCC$ ATTRIBUTES STDCALL :: c_get_procedure
            integer(kind=c_intptr_t), value      :: handle
            character(kind=c_char), dimension(*) :: cname
            type(c_funptr)                       :: c_get_procedure
        end function c_get_procedure
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

    handle = c_load_library( cname )
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
    integer(kind=c_long)           :: handle
    character(len=1), dimension(*) :: cname
    type(c_funptr)                 :: cproc
    logical                        :: success

    cproc = c_get_procedure( handle, cname )

    success = c_associated(cproc)

end subroutine system_get_procedure

end module system_dynamic_libraries
