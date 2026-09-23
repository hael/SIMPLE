!@descr: extension type providing typed convenience accessors for ui_program
! Wraps the string->polymorphic vrefhash so that the UI program tables store and
! return type(ui_program) pointers without a select type at every call site.
! Only the overloads production uses are kept: set by character key
! (simple_ui_utils::add_ui_program) and get by string key (simple_ui).
module simple_ui_hash
use simple_vrefhash,   only: vrefhash    ! core polymorphic hash
use simple_string,     only: string      ! string helper
use simple_ui_program, only: ui_program  ! UI type
implicit none
private
public :: ui_hash

type, extends(vrefhash) :: ui_hash
contains
    procedure :: set_ui_program
    procedure :: get_ui_program
end type ui_hash

contains

    subroutine set_ui_program( self, key, obj )
        class(ui_hash),           intent(inout) :: self
        character(len=*),         intent(in)    :: key
        type(ui_program), target, intent(inout) :: obj
        class(*), pointer :: p
        character(:), allocatable :: k
        k = trim(adjustl(key))
        p => obj
        call self%set_ref(k, p)
    end subroutine set_ui_program

    subroutine get_ui_program( self, key, pobj, found )
        class(ui_hash),            intent(in)  :: self
        type(string),              intent(in)  :: key
        type(ui_program), pointer, intent(out) :: pobj
        logical, optional,         intent(out) :: found
        class(*), pointer :: ptmp
        logical :: lfound
        character(:), allocatable :: k
        nullify(pobj)
        k = adjustl(key%to_char())
        call self%get_ref(k, ptmp, lfound)
        if (present(found)) found = lfound
        if (.not. lfound) return
        select type (ptmp)
            type is (ui_program)
                pobj => ptmp
            class default
                nullify(pobj)
                if (present(found)) found = .false.
        end select
    end subroutine get_ui_program

end module simple_ui_hash
