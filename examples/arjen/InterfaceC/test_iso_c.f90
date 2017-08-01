! test_iso_c.f90 --
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
!     Small program to test/illustrate the iso_c_binding module
!
!     Note:
!     The C routine that is called is a dummy routine
!     The most important aspect is: do the names on both sides match?
!
program test_iso_c
    use, intrinsic :: iso_c_binding

    implicit none

    interface
        function sqlite3_column_count( stmt ) bind(C)
            use, intrinsic     :: iso_c_binding
            type(c_ptr), value :: stmt
            integer(c_int)     :: sqlite3_column_count
        end function sqlite3_column_count
    end interface

    interface
        function sqlite3_do_c( handle, command, &
                             errmsg, len_errmsg ) bind(C, name = 'Sqlite3DoC')
            use, intrinsic                       :: iso_c_binding
            type(c_ptr), value                   :: handle
            character(kind=c_char), dimension(*) :: command
            character(kind=c_char), dimension(*) :: errmsg
            integer(kind=c_int), value           :: len_errmsg
            integer(kind=c_int)                  :: sqlite3_do_c
        end function sqlite3_do_c
    end interface

    type(c_ptr) :: stmt
    type(c_ptr) :: handle
    character(len=20) :: command
    character(len=20) :: errmsg
    integer           :: errcode

    write(*,*) 'Count: ', sqlite3_column_count(stmt)

    command = 'Dummy command' // char(0)
    errcode = sqlite3_do_c( handle, command, errmsg, len(errmsg) )
    write(*,*) 'Error message: ', errcode, errmsg
    write(*,*) 'Character at 14:', ichar(errmsg(14:14))

end program test_iso_c
