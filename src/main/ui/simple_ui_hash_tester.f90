!@descr: unit test routines for the typed UI program table (simple_ui_hash)
! Pins the two accessors production relies on: set_ui_program by character key
! (simple_ui_utils::add_ui_program) and get_ui_program by string key (simple_ui):
! reference semantics, the typed miss (wrong dynamic type or absent key gives
! found=.false. and a null pointer), overwrite retargeting the pointer, key
! trimming and the inherited vrefhash bookkeeping.
module simple_ui_hash_tester
use simple_ui_hash,    only: ui_hash
use simple_ui_program, only: ui_program
use simple_string,     only: string
use simple_test_utils
implicit none
private
public :: run_all_ui_hash_tests

contains

    subroutine run_all_ui_hash_tests()
        write(*,'(A)') '**** running all ui_hash tests ****'
        call test_set_get_by_key()
        call test_missing_key_is_a_typed_miss()
        call test_wrong_dynamic_type_is_a_typed_miss()
        call test_overwrite_retargets_pointer()
        call test_keys_are_trimmed()
        call test_found_is_optional()
    end subroutine run_all_ui_hash_tests

    subroutine make_prg( prg, name )
        type(ui_program), intent(inout) :: prg
        character(len=*), intent(in)    :: name
        call prg%new(name, 'summary text long enough to pass the length check', &
            &'help text', 'simple_exec', .false.)
    end subroutine make_prg

    subroutine test_set_get_by_key()
        type(ui_hash)              :: tab
        type(ui_program), target   :: prg1, prg2
        type(ui_program), pointer  :: p
        logical :: found
        write(*,'(A)') 'test_set_get_by_key'
        call make_prg(prg1, 'abinitio2D')
        call make_prg(prg2, 'refine3D')
        call tab%set_ui_program('abinitio2D', prg1)
        call tab%set_ui_program('refine3D',   prg2)
        call assert_int(2, tab%count(), 'two programs stored')
        call assert_true(tab%has_key('abinitio2D'), 'has_key sees the stored program')
        call tab%get_ui_program(string('abinitio2D'), p, found)
        call assert_true(found, 'get by string key finds the program')
        call assert_true(associated(p, prg1), 'pointer targets the stored object, not a copy')
        call assert_string_eq('abinitio2D', p%name, 'the right program comes back')
        call tab%get_ui_program(string('refine3D'), p, found)
        call assert_true(associated(p, prg2), 'second key resolves to its own object')
        ! reference semantics: a change through the pointer is visible in the object
        p%visibility = 42
        call assert_int(42, prg2%visibility, 'the table stores a reference, not a copy')
        call tab%clear()
    end subroutine test_set_get_by_key

    subroutine test_missing_key_is_a_typed_miss()
        type(ui_hash)              :: tab
        type(ui_program), target   :: prg1
        type(ui_program), pointer  :: p
        logical :: found
        write(*,'(A)') 'test_missing_key_is_a_typed_miss'
        call make_prg(prg1, 'abinitio2D')
        call tab%set_ui_program('abinitio2D', prg1)
        p => prg1   ! must be reset by the call, intent(out)
        found = .true.
        call tab%get_ui_program(string('nonexistent'), p, found)
        call assert_false(found, 'absent key gives found=.false.')
        call assert_false(associated(p), 'absent key gives a null pointer')
        call tab%clear()
    end subroutine test_missing_key_is_a_typed_miss

    subroutine test_wrong_dynamic_type_is_a_typed_miss()
        type(ui_hash)              :: tab
        integer,          target   :: not_a_program
        type(ui_program), pointer  :: p
        logical :: found
        write(*,'(A)') 'test_wrong_dynamic_type_is_a_typed_miss'
        not_a_program = 7
        ! stored through the vrefhash base as a polymorphic reference of the wrong type
        call tab%set_ref('oddity', not_a_program)
        call assert_true(tab%has_key('oddity'), 'the base hash holds the key')
        found = .true.
        call tab%get_ui_program(string('oddity'), p, found)
        call assert_false(found, 'wrong dynamic type is reported as not found')
        call assert_false(associated(p), 'wrong dynamic type gives a null pointer')
        call tab%clear()
    end subroutine test_wrong_dynamic_type_is_a_typed_miss

    subroutine test_overwrite_retargets_pointer()
        type(ui_hash)              :: tab
        type(ui_program), target   :: old_prg, new_prg
        type(ui_program), pointer  :: p
        logical :: found
        write(*,'(A)') 'test_overwrite_retargets_pointer'
        call make_prg(old_prg, 'cluster2D')
        call make_prg(new_prg, 'cluster2D')
        new_prg%visibility = 99
        call tab%set_ui_program('cluster2D', old_prg)
        call tab%set_ui_program('cluster2D', new_prg)
        call assert_int(1, tab%count(), 'overwrite keeps one entry')
        call tab%get_ui_program(string('cluster2D'), p, found)
        call assert_true(found, 'overwritten key still found')
        call assert_true(associated(p, new_prg), 'pointer targets the new object')
        call assert_false(associated(p, old_prg), 'pointer no longer targets the old object')
        call assert_int(99, p%visibility, 'the new object is the one read back')
        call tab%clear()
    end subroutine test_overwrite_retargets_pointer

    subroutine test_keys_are_trimmed()
        type(ui_hash)              :: tab
        type(ui_program), target   :: prg1
        type(ui_program), pointer  :: p
        logical :: found
        write(*,'(A)') 'test_keys_are_trimmed'
        call make_prg(prg1, 'symmetrize_map')
        call tab%set_ui_program('   symmetrize_map  ', prg1)
        call assert_true(tab%has_key('symmetrize_map'), 'leading and trailing blanks are dropped on set')
        call tab%get_ui_program(string('  symmetrize_map'), p, found)
        call assert_true(found, 'leading blanks are dropped on get')
        call assert_true(associated(p, prg1), 'trimmed key resolves to the object')
        call tab%clear()
    end subroutine test_keys_are_trimmed

    subroutine test_found_is_optional()
        type(ui_hash)              :: tab
        type(ui_program), target   :: prg1
        type(ui_program), pointer  :: p
        write(*,'(A)') 'test_found_is_optional'
        call make_prg(prg1, 'postprocess')
        call tab%set_ui_program('postprocess', prg1)
        call tab%get_ui_program(string('postprocess'), p)
        call assert_true(associated(p, prg1), 'hit without the found argument')
        call tab%get_ui_program(string('missing'), p)
        call assert_false(associated(p), 'miss without the found argument gives a null pointer')
        call tab%clear()
    end subroutine test_found_is_optional

end module simple_ui_hash_tester
