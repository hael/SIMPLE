!@descr: extension type providing typed convenience accessors for ui_param & ui_program
module simple_ui_hash
use simple_vrefhash,   only: vrefhash             ! core polymorphic hash
use simple_string,     only: string               ! string helper
use simple_ui_program, only: ui_param, ui_program ! UI types
implicit none
private
public :: ui_hash, test_ui_hash

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

    subroutine set_ref_ui_param_char( self, key, obj )
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

    subroutine set_ref_ui_param_str( self, key, obj)
        class(ui_hash),         intent(inout) :: self
        type(string),           intent(in)    :: key
        type(ui_param), target, intent(inout) :: obj
        class(*), pointer :: p
        character(:), allocatable :: k
        k = adjustl(key%to_char())
        p => obj
        call self%set_ref(k, p)
    end subroutine set_ref_ui_param_str

    subroutine get_ref_ui_param_char( self, key, pobj, found)
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

    subroutine get_ref_ui_param_str( self, key, pobj, found )
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

    subroutine set_ref_ui_program_char( self, key, obj )
        class(ui_hash),           intent(inout) :: self
        character(len=*),         intent(in)    :: key
        type(ui_program), target, intent(inout) :: obj
        class(*), pointer :: p
        character(:), allocatable :: k
        k = trim(adjustl(key))
        p => obj
        call self%set_ref(k, p)
    end subroutine set_ref_ui_program_char

    subroutine set_ref_ui_program_str( self, key, obj )
        class(ui_hash),           intent(inout) :: self
        type(string),             intent(in)    :: key
        type(ui_program), target, intent(inout) :: obj
        class(*), pointer :: p
        character(:), allocatable :: k
        k = adjustl(key%to_char())
        p => obj
        call self%set_ref(k, p)
    end subroutine set_ref_ui_program_str

    subroutine get_ref_ui_program_char( self, key, pobj, found )
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

    subroutine get_ref_ui_program_str( self, key, pobj, found )
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

    !> Standalone test subroutine for ui_hash: comprehensive unit testing
    subroutine test_ui_hash()
        type(ui_hash)     :: hash
        type(ui_param)    :: param1, param2, param_new
        type(ui_program)  :: prog1, prog2
        type(ui_param), pointer    :: pparam
        type(ui_program), pointer  :: pprog
        type(string)      :: str_key
        logical           :: found
        integer           :: ntests, npassed
        ntests  = 0
        npassed = 0
        write(*,'(A)') '======================================'
        write(*,'(A)') '>>> UI_HASH UNIT TESTS'
        write(*,'(A)') '======================================'
        ! ===== TEST 1: Set and get ui_param with character key =====
        ntests = ntests + 1
        write(*,'(A)', advance='no') 'TEST 1: set/get ui_param (char key)... '
        call hash%set_ui_param('param1', param1)
        call hash%get_ui_param('param1', pparam, found)
        if (associated(pparam) .and. found) then
            npassed = npassed + 1
            write(*,'(A)') 'PASS'
        else
            write(*,'(A)') 'FAIL - pointer not associated or found=.false.'
        endif
        ! ===== TEST 2: Set and get ui_param with string key =====
        ntests = ntests + 1
        write(*,'(A)', advance='no') 'TEST 2: set/get ui_param (string key)... '
        str_key = string('param2')
        call hash%set_ui_param(str_key, param2)
        call hash%get_ui_param(str_key, pparam, found)
        if (associated(pparam) .and. found) then
            npassed = npassed + 1
            write(*,'(A)') 'PASS'
        else
            write(*,'(A)') 'FAIL - pointer not associated or found=.false.'
        endif
        ! ===== TEST 3: Set and get ui_program with character key =====
        ntests = ntests + 1
        write(*,'(A)', advance='no') 'TEST 3: set/get ui_program (char key)... '
        call hash%set_ui_program('prog1', prog1)
        call hash%get_ui_program('prog1', pprog, found)
        if (associated(pprog) .and. found) then
            npassed = npassed + 1
            write(*,'(A)') 'PASS'
        else
            write(*,'(A)') 'FAIL - pointer not associated or found=.false.'
        endif
        ! ===== TEST 4: Set and get ui_program with string key =====
        ntests = ntests + 1
        write(*,'(A)', advance='no') 'TEST 4: set/get ui_program (string key)... '
        str_key = string('prog2')
        call hash%set_ui_program(str_key, prog2)
        call hash%get_ui_program(str_key, pprog, found)
        if (associated(pprog) .and. found) then
            npassed = npassed + 1
            write(*,'(A)') 'PASS'
        else
            write(*,'(A)') 'FAIL - pointer not associated or found=.false.'
        endif
        ! ===== TEST 5: Get non-existent key returns found=.false. =====
        ntests = ntests + 1
        write(*,'(A)', advance='no') 'TEST 5: non-existent key returns found=.false.... '
        call hash%get_ui_param('nonexistent_key', pparam, found)
        if (.not. found .and. .not. associated(pparam)) then
            npassed = npassed + 1
            write(*,'(A)') 'PASS'
        else
            write(*,'(A)') 'FAIL - found should be .false. and pparam nullified'
        endif
        ! ===== TEST 6: Type mismatch handling (get program when param stored) =====
        ntests = ntests + 1
        write(*,'(A)', advance='no') 'TEST 6: type mismatch handling (prog vs param)... '
        call hash%get_ui_program('param1', pprog, found)
        if (.not. found .and. .not. associated(pprog)) then
            npassed = npassed + 1
            write(*,'(A)') 'PASS'
        else
            write(*,'(A)') 'FAIL - type mismatch should return found=.false.'
        endif
        ! ===== TEST 7: Type mismatch handling (get param when program stored) =====
        ntests = ntests + 1
        write(*,'(A)', advance='no') 'TEST 7: type mismatch handling (param vs prog)... '
        call hash%get_ui_param('prog1', pparam, found)
        if (.not. found .and. .not. associated(pparam)) then
            npassed = npassed + 1
            write(*,'(A)') 'PASS'
        else
            write(*,'(A)') 'FAIL - type mismatch should return found=.false.'
        endif
        ! ===== TEST 8: Overwriting existing key =====
        ntests = ntests + 1
        write(*,'(A)', advance='no') 'TEST 8: overwriting existing key... '
        call hash%set_ui_param('param1', param_new)
        call hash%get_ui_param('param1', pparam, found)
        if (associated(pparam) .and. found) then
            npassed = npassed + 1
            write(*,'(A)') 'PASS'
        else
            write(*,'(A)') 'FAIL - overwritten reference not found'
        endif
        ! ===== Summary =====
        write(*,'(A)') '======================================'
        write(*,'(A,I3,A,I3,A)') 'RESULTS: ', npassed, '/', ntests, ' tests passed'
        write(*,'(A)') '======================================'
        nullify(pparam)
        nullify(pprog)
    end subroutine test_ui_hash

end module simple_ui_hash
