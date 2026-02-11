!@descr: extension type providing typed convenience accessors for ui_param & ui_program
module simple_ui_hash
use simple_vrefhash,   only: vrefhash             ! core polymorphic hash
use simple_string,     only: string               ! string helper
use simple_ui_program, only: ui_param, ui_program ! UI types
implicit none
private
public :: ui_hash

!----------------------------------------------------------------
! Extended type: all convenience methods live inside.
!----------------------------------------------------------------
type, extends(vrefhash) :: ui_hash
contains
    procedure, private :: set_ref_ui_param_char
    procedure, private :: set_ref_ui_param_str
    generic            :: set_ui_param => set_ref_ui_param_char, set_ref_ui_param_str
    procedure, private :: get_ref_ui_param_char
    procedure, private :: get_ref_ui_param_str
    generic            :: get_ui_param => get_ref_ui_param_char, get_ref_ui_param_str
    procedure, private :: set_ref_ui_program_char
    procedure, private :: set_ref_ui_program_str
    generic            :: set_ui_program => set_ref_ui_program_char, set_ref_ui_program_str
    procedure, private :: get_ref_ui_program_char
    procedure, private :: get_ref_ui_program_str
    generic            :: get_ui_program => get_ref_ui_program_char, get_ref_ui_program_str
end type ui_hash

contains

    subroutine set_ref_ui_param_char(self, key, obj)
        class(ui_hash),         intent(inout) :: self
        character(len=*),       intent(in)    :: key
        type(ui_param), target, intent(inout) :: obj
        class(*), pointer :: p
        character(:), allocatable :: k
        k = trim(adjustl(key))
        p => obj
        ! Core API: set polymorphic pointer by string key
        call self%set_ref(k, p)
    end subroutine set_ref_ui_param_char

    subroutine set_ref_ui_param_str(self, key, obj)
        class(ui_hash),         intent(inout) :: self
        type(string),           intent(in)    :: key
        type(ui_param), target, intent(inout) :: obj
        class(*), pointer :: p
        character(:), allocatable :: k
        k = adjustl(key%to_char())
        p => obj
        call self%set_ref(k, p)
    end subroutine set_ref_ui_param_str

    subroutine get_ref_ui_param_char(self, key, pobj, found)
        class(ui_hash),          intent(in)  :: self
        character(len=*),        intent(in)  :: key
        type(ui_param), pointer, intent(out) :: pobj
        logical, optional,       intent(out) :: found
        class(*), pointer :: ptmp
        logical :: lfound
        character(:), allocatable :: k
        nullify(pobj)
        k = trim(adjustl(key))
        call self%get_ref(k, ptmp, lfound) ! core getter: returns polymorphic pointer and found flag
        if (present(found)) found = lfound
        if (.not. lfound) return
        select type (ptmp)
            type is (ui_param)
                pobj => ptmp
            class default
                ! mismatched dynamic type -> signal not found to caller
                nullify(pobj)
                if (present(found)) found = .false.
        end select
    end subroutine get_ref_ui_param_char

    subroutine get_ref_ui_param_str(self, key, pobj, found)
        class(ui_hash),          intent(in)  :: self
        type(string),            intent(in)  :: key
        type(ui_param), pointer, intent(out) :: pobj
        logical, optional,       intent(out) :: found
        class(*), pointer :: ptmp
        logical :: lfound
        character(:), allocatable :: k
        nullify(pobj)
        k = adjustl(key%to_char())
        call self%get_ref(k, ptmp, lfound)
        if (present(found)) found = lfound
        if (.not. lfound) return
        select type (ptmp)
            type is (ui_param)
                pobj => ptmp
          class default
                ! mismatched dynamic type -> signal not found to caller
                nullify(pobj)
                if (present(found)) found = .false.
        end select
    end subroutine get_ref_ui_param_str

    subroutine set_ref_ui_program_char(self, key, obj)
        class(ui_hash),           intent(inout) :: self
        character(len=*),         intent(in)    :: key
        type(ui_program), target, intent(inout) :: obj
        class(*), pointer :: p
        character(:), allocatable :: k
        k = trim(adjustl(key))
        p => obj
        call self%set_ref(k, p)
    end subroutine set_ref_ui_program_char

    subroutine set_ref_ui_program_str(self, key, obj)
        class(ui_hash),           intent(inout) :: self
        type(string),             intent(in)    :: key
        type(ui_program), target, intent(inout) :: obj
        class(*), pointer :: p
        character(:), allocatable :: k
        k = adjustl(key%to_char())
        p => obj
        call self%set_ref(k, p)
    end subroutine set_ref_ui_program_str

    subroutine get_ref_ui_program_char(self, key, pobj, found)
        class(ui_hash),            intent(in)  :: self
        character(len=*),          intent(in)  :: key
        type(ui_program), pointer, intent(out) :: pobj
        logical, optional,         intent(out) :: found
        class(*), pointer :: ptmp
        logical :: lfound
        character(:), allocatable :: k
        nullify(pobj)
        k = trim(adjustl(key))
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
    end subroutine get_ref_ui_program_char

    subroutine get_ref_ui_program_str(self, key, pobj, found)
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
    end subroutine get_ref_ui_program_str

end module simple_ui_hash
