! dynlib.f90 --
!     Module to handle procedures from dynamic libraries
!     Requires Fortran 2003
!
!     Note:
!     There are two implementations of the basic routines,
!     one for Windows and one for Linux/OSX, as the
!     underlying system routines have different APIs
!
!     TODO:
!     Register the contents of an internal library
!     and retrieving that
!
!     TODO:
!     Handle module names
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
module dynamic_libraries
    use iso_c_binding
    use system_dynamic_libraries

    implicit none

    private

    type dynamic_library
        private
        integer(kind=c_long) :: handle   = 0
        logical              :: loaded   = .false.
        logical              :: internal = .false.
        character(len=512)   :: name
    end type dynamic_library

    public :: dynamic_library
    public :: load_library
    public :: get_procedure
    public :: get_procedure_global
    !public :: register_internal_library
    !public :: register_procedure

    !public :: test_it

    interface get_procedure
        module procedure get_procedure_global
        module procedure get_procedure_module
    end interface

contains

! load_library --
!     Load a dynamic library (DLL/SO)
!
! Arguments:
!     dynlib            Variable holding the handle to the library
!     name              Name of the library to load
!     success           Whether it was loaded successfully or not
!
subroutine load_library( dynlib, name, success )

    type(dynamic_library), intent(inout) :: dynlib
    character(len=*), intent(in)           :: name
    logical, intent(out)                   :: success

    character(kind=c_char), dimension(512) :: cname
    integer                                :: i

    success = .false.

    dynlib%loaded   = .false.
    dynlib%internal = .false.
    dynlib%name     = name

    cname(1:len(name)+1) = (/ ( name(i:i), i = 1,len(name) ), char(0) /)

    call system_load_library( dynlib%handle, cname )

    if ( dynlib%handle /= 0_c_long ) then
        dynlib%loaded = .true.
        success = .true.
    endif
end subroutine load_library

! get_procedure_global --
!     Get a procedure pointer from a loaded dynamic library (DLL/SO)
!     (outside a module)
!
! Arguments:
!     dynlib            Variable holding the handle to the library
!     module            Name of the module that contains it (if any)
!     name              Name of the procedure to get
!     proc              Pointer to the procedure
!     success           Whether it was loaded successfully or not
!
subroutine get_procedure_global( dynlib, name, proc, success )

    type(dynamic_library), intent(inout)   :: dynlib
    character(len=*), intent(in)           :: name
    procedure(), pointer                   :: proc
    logical, intent(out)                   :: success

    call get_procedure_module( dynlib, '', name, proc, success )
end subroutine get_procedure_global

! get_procedure_module --
!     Get a procedure pointer from a loaded dynamic library (DLL/SO)
!
! Arguments:
!     dynlib            Variable holding the handle to the library
!     module            Name of the module that contains it (if any)
!     name              Name of the procedure to get
!     proc              Pointer to the procedure
!     success           Whether it was loaded successfully or not
!
subroutine get_procedure_module( dynlib, module, name, proc, success )

    type(dynamic_library), intent(inout)   :: dynlib
    character(len=*), intent(in)           :: module
    character(len=*), intent(in)           :: name
    procedure(), pointer                   :: proc
    logical, intent(out)                   :: success

    character(kind=c_char,len=512)         :: cname
    type(c_funptr)                         :: cproc
    integer                                :: i
    integer                                :: number_cases

    success = .false.
    proc    => null()

    number_cases = 4
    if ( module /= ' ' ) then
        number_cases = 2
    endif

    if ( .not. dynlib%loaded ) then
        return
    endif

    if ( .not. dynlib%internal ) then
        do i = 1,number_cases
            call convert_name( i, trim(module), trim(name), cname )
            call system_get_procedure( dynlib%handle, cname, cproc, success )

            if ( success  ) then
                call c_f_procpointer( cproc, proc )
                exit
            endif
        enddo
    else
        ! TODO
    endif
end subroutine get_procedure_module

! convert_name --
!     Convert the name of the procedure to a possible external name
!     (including the name of the module)
!
! Arguments:
!     type              Type of conversion
!     module            Name of the module that contains it (if any)
!     name              Name of the procedure to get
!     cname             Constructed name (suitable for passing to C)
!
subroutine convert_name( type, module, name, cname )
    integer, intent(in)          :: type
    character(len=*), intent(in) :: module
    character(len=*), intent(in) :: name
    character(kind=c_char,len=*) :: cname

    character(len=len(cname))    :: module_name
    integer                      :: i

    if ( module == ' ' ) then
        select case ( type )
            case ( 1 )
                cname = tolower(name) // char(0)

            case ( 2 )
                cname = tolower(name) // '_' // char(0)

            case ( 3 )
                cname = tolower(name) // '__' // char(0)

            case ( 4 )
                cname = toupper(name) // char(0)
        end select
    else
        select case ( type )
            case ( 1 )
                module_name = '__' // tolower(module) // '_MOD_' // tolower(name)

            case ( 2 )
                module_name = toupper(module) // '_mp_' // tolower(name)
        end select

        cname = trim(module_name) // char(0)
    endif
end subroutine convert_name

! tolower --
!     Return a string in lower case
!
! Arguments:
!     string         String to be converted
!
function tolower(string)
    character(len=*), intent(in) :: string
    character(len=len(string))   :: tolower

    integer                      :: i
    integer                      :: ch

    tolower = string
    do i = 1,len(string)
        ch = iachar(string(i:i))
        if ( ch >= 65 .and. ch <= 90 ) then
            tolower(i:i) = achar(ch + 32)
        endif
    enddo
end function tolower

! toupper --
!     Return a string in upper case
!
! Arguments:
!     string         String to be converted
!
function toupper(string)
    character(len=*), intent(in) :: string
    character(len=len(string))   :: toupper

    integer                      :: i
    integer                      :: ch

    toupper = string
    do i = 1,len(string)
        ch = iachar(string(i:i))
        if ( ch >= 97 .and. ch <= 122 ) then
            toupper(i:i) = achar(ch - 32)
        endif
    enddo
end function toupper

! Test code
!
! test_it --
!     Subroutine to test the various conversion cases
!
!subroutine test_it
!    character(kind=c_char,len=512) :: cname
!
!    integer               :: i
!    integer, dimension(1) :: pos
!
!    do i = 1,4
!        call convert_name( i, '', 'ProcedureName', cname )
!        pos = index( cname, char(0) )
!        write( *, '(60a)' ) cname(1:pos(1)-1)
!    enddo
!
!    do i = 1,2
!        call convert_name( i, 'ModuleName', 'ProcedureName', cname )
!        pos = index( cname, char(0) )
!        write( *, '(60a)' ) cname(1:pos(1)-1)
!    enddo
!end subroutine test_it
!
end module dynamic_libraries
!
!program test_dynlib
!    use dynamic_libraries
!
!    call test_it
!
!end program test_dynlib
