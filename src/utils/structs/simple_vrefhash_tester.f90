!@descr: unit test routines for vrefhash (string->polymorphic reference hash)
module simple_vrefhash_tester
use simple_vrefhash
use simple_string
use simple_test_utils
implicit none
private
public :: run_all_vrefhash_tests

type :: my_cfg
    integer :: n = 0
    real    :: x = 0.0
end type my_cfg

type :: my_obj
    integer :: counter = 0
contains
    procedure :: inc
end type my_obj

contains

    subroutine inc(self, k)
        class(my_obj), intent(inout) :: self
        integer,       intent(in)    :: k
        self%counter = self%counter + k
    end subroutine inc

    subroutine run_all_vrefhash_tests()
        write(*,'(A)') '**** running all vrefhash tests ****'
        call test_init_clear_destroy_and_count()
        call test_set_get_has_key_with_char_and_string_keys()
        call test_reference_semantics_updates_visible()
        call test_replace_pointer_on_same_key()
        call test_store_multiple_dynamic_types()
        call test_delete()
        call test_missing_key_get_ref_returns_null()
        call test_collision_sanity_many_keys()
        call test_keys_returns_all_keys()
        call test_keys_sorted_returns_sorted_keys()
        ! call report_summary()
    end subroutine run_all_vrefhash_tests

    !---------------- lifecycle / count ----------------

    subroutine test_init_clear_destroy_and_count()
        type(vrefhash) :: h
        integer, target :: i
        class(*), pointer :: p
        logical :: ok
        write(*,'(A)') 'test_init_clear_destroy_and_count'
        call h%init(16)
        call assert_int(0, h%count(), 'count() after init == 0')
        i = 7
        call h%set_ref('i', i)
        call assert_int(1, h%count(), 'count() after one insert == 1')
        call h%clear()
        call assert_int(0, h%count(), 'count() after clear() == 0')
        call h%get_ref('i', p, ok)
        call assert_true(.not. ok, 'after clear, key not found')
        call h%destroy()
        call assert_int(0, h%count(), 'count() after destroy() == 0 (safe)')
    end subroutine test_init_clear_destroy_and_count

    !---------------- set/get/has_key with char + string keys ----------------

    subroutine test_set_get_has_key_with_char_and_string_keys()
        type(vrefhash) :: h
        integer, target :: i
        class(*), pointer :: p
        logical :: ok
        type(string) :: skey
        write(*,'(A)') 'test_set_get_has_key_with_char_and_string_keys'
        call h%init(32)
        i = 42
        call h%set_ref('answer', i)
        call assert_true(h%has_key('answer'), 'has_key(char) true')
        call assert_true(.not. h%has_key('missing'), 'has_key(char) false')
        call h%get_ref('answer', p, ok)
        call assert_true(ok, 'get_ref(char) found')
        call assert_true(associated(p), 'get_ref returns associated pointer')
        select type(t => p)
        type is (integer)
            call assert_int(42, t, 'stored integer value is 42')
        class default
            call assert_true(.false., 'wrong dynamic type for "answer"')
        end select
        skey = string('  answer  ')
        call assert_true(h%has_key(skey), 'has_key(string) true (string has spaces)')
        call h%get_ref(skey, p, ok)
        call assert_true(ok, 'get_ref(string) found')
        select type(t => p)
        type is (integer)
            call assert_int(42, t, 'get_ref(string) returns same integer')
        class default
            call assert_true(.false., 'wrong dynamic type for string key get')
        end select
        call h%destroy()
    end subroutine test_set_get_has_key_with_char_and_string_keys

    !---------------- reference semantics: updates visible ----------------

    subroutine test_reference_semantics_updates_visible()
        type(vrefhash) :: h
        type(my_cfg), target :: cfg
        class(*), pointer :: p
        logical :: ok
        write(*,'(A)') 'test_reference_semantics_updates_visible'
        call h%init(64)
        cfg%n = 1
        cfg%x = 1.5
        call h%set_ref('cfg', cfg)
        ! mutate original
        cfg%n = 99
        cfg%x = 3.25
        call h%get_ref('cfg', p, ok)
        call assert_true(ok, 'cfg found')
        select type(c => p)
        type is (my_cfg)
            call assert_int(99, c%n, 'mutated cfg%n visible via hash ref')
            call assert_real(3.25, c%x, 1.0e-6, 'mutated cfg%x visible via hash ref')
        class default
            call assert_true(.false., 'wrong dynamic type for cfg')
        end select

        call h%destroy()
    end subroutine test_reference_semantics_updates_visible

    !---------------- replace pointer on same key ----------------

    subroutine test_replace_pointer_on_same_key()
        type(vrefhash) :: h
        integer, target :: a, b
        class(*), pointer :: p
        logical :: ok
        write(*,'(A)') 'test_replace_pointer_on_same_key'
        call h%init(32)
        a = 1
        b = 2
        call h%set_ref('k', a)
        call assert_int(1, h%count(), 'count=1 after first set_ref')
        call h%get_ref('k', p, ok)
        call assert_true(ok, 'k found')
        select type(t => p)
        type is (integer)
            call assert_int(1, t, 'k points to a initially')
        class default
            call assert_true(.false., 'wrong type for k initially')
        end select
        ! replace reference under same key (count must not grow)
        call h%set_ref('k', b)
        call assert_int(1, h%count(), 'count unchanged after replacing same key')
        call h%get_ref('k', p, ok)
        call assert_true(ok, 'k found after replace')
        select type(t => p)
        type is (integer)
            call assert_int(2, t, 'k now points to b')
        class default
            call assert_true(.false., 'wrong type for k after replace')
        end select
        ! mutate b; should be visible via hash
        b = 77
        call h%get_ref('k', p, ok)
        select type(t => p)
        type is (integer)
            call assert_int(77, t, 'mutation of b visible via hash')
        class default
            call assert_true(.false., 'wrong type after mutation')
        end select
        call h%destroy()
    end subroutine test_replace_pointer_on_same_key

    !---------------- store multiple dynamic types ----------------

    subroutine test_store_multiple_dynamic_types()
        type(vrefhash) :: h
        integer,  target :: i
        real,     target :: r
        logical,  target :: l
        character(:), allocatable, target :: s
        type(my_obj), target :: o
        class(*), pointer :: p
        logical :: ok
        write(*,'(A)') 'test_store_multiple_dynamic_types'
        call h%init(128)
        i = 10
        r = 2.5
        l = .true.
        s = 'hello'
        o%counter = 5
        call h%set_ref('i', i)
        call h%set_ref('r', r)
        call h%set_ref('l', l)
        call h%set_ref('s', s)
        call h%set_ref('o', o)
        call assert_int(5, h%count(), 'count=5 after inserting five keys')
        call h%get_ref('i', p, ok)
        select type(t => p); type is (integer); call assert_int(10, t, 'i==10'); class default; call assert_true(.false., 'i type'); end select
        call h%get_ref('r', p, ok)
        select type(t => p); type is (real);    call assert_real(2.5, t, 1.0e-6, 'r==2.5'); class default; call assert_true(.false., 'r type'); end select
        call h%get_ref('l', p, ok)
        select type(t => p); type is (logical); call assert_true(t, 'l==true'); class default; call assert_true(.false., 'l type'); end select
        call h%get_ref('s', p, ok)
        select type(t => p)
        type is (character(*))
            call assert_true(t == 'hello', 's=="hello"')
        class default
            call assert_true(.false., 's type')
        end select
        ! derived type with method; ensure we can mutate through original and see it
        call o%inc(3)  ! o%counter = 8
        call h%get_ref('o', p, ok)
        select type(t => p)
        type is (my_obj)
            call assert_int(8, t%counter, 'derived object mutation visible')
        class default
            call assert_true(.false., 'o type')
        end select
        call h%destroy()
    end subroutine test_store_multiple_dynamic_types

    !---------------- delete ----------------

    subroutine test_delete()
        type(vrefhash) :: h
        integer, target :: a, b
        class(*), pointer :: p
        logical :: ok
        write(*,'(A)') 'test_delete'
        call h%init(32)
        a = 1; b = 2
        call h%set_ref('a', a)
        call h%set_ref('b', b)
        call assert_int(2, h%count(), 'count=2 before del')
        call h%del('a')
        call assert_int(1, h%count(), 'count=1 after del(a)')
        call assert_true(.not. h%has_key('a'), 'a removed')
        call assert_true(h%has_key('b'), 'b remains')
        call h%get_ref('a', p, ok)
        call assert_true(.not. ok, 'get_ref(a) not found after delete')
        call h%del('missing')  ! no-op
        call assert_int(1, h%count(), 'delete missing key is no-op')
        call h%destroy()
    end subroutine test_delete

    !---------------- missing key behavior ----------------

    subroutine test_missing_key_get_ref_returns_null()
        type(vrefhash) :: h
        class(*), pointer :: p
        logical :: ok
        write(*,'(A)') 'test_missing_key_get_ref_returns_null'
        call h%init(16)
        call h%get_ref('nope', p, ok)
        call assert_true(.not. ok, 'missing key returns ok=.false.')
        call assert_true(.not. associated(p), 'missing key returns null pointer')
        call h%destroy()
    end subroutine test_missing_key_get_ref_returns_null

    !---------------- collision sanity (many keys) ----------------

    subroutine test_collision_sanity_many_keys()
        type(vrefhash) :: h
        integer, parameter :: nkeys = 200
        integer :: i
        integer, target :: vals(nkeys)
        character(len=32) :: k
        class(*), pointer :: p
        logical :: ok
        write(*,'(A)') 'test_collision_sanity_many_keys'
        ! small bucket count forces collisions
        call h%init(8)
        do i = 1, nkeys
            write(k,'(A,I0)') 'key_', i
            vals(i) = i * 10
            call h%set_ref(trim(k), vals(i))
        end do
        call assert_int(nkeys, h%count(), 'count==nkeys after many inserts')
        ! spot check a few
        do i = 1, nkeys, 37
            write(k,'(A,I0)') 'key_', i
            call h%get_ref(trim(k), p, ok)
            call assert_true(ok, 'spot check found: '//trim(k))
            select type(t => p)
            type is (integer)
                call assert_int(i*10, t, 'spot check value matches for '//trim(k))
            class default
                call assert_true(.false., 'spot check wrong type for '//trim(k))
            end select
        end do
        call h%destroy()
    end subroutine test_collision_sanity_many_keys

    !---------------- keys() ----------------

    subroutine test_keys_returns_all_keys()
        type(vrefhash) :: h
        integer, target :: a, b, c
        type(string), allocatable :: ks(:)
        logical :: fa, fb, fc
        integer :: i
        write(*,'(A)') 'test_keys_returns_all_keys'
        call h%init(32)
        a = 1; b = 2; c = 3
        call h%set_ref('  alpha  ', a)          ! trimmed/adjustl in set_ref_char
        call h%set_ref(string('beta'), b)       ! adjustl(to_char()) in set_ref_str
        call h%set_ref('GAMMA', c)
        ks = h%keys()
        call assert_int(h%count(), size(ks), 'keys() size == count()')
        fa = .false.; fb = .false.; fc = .false.
        do i = 1, size(ks)
            if (ks(i) == 'alpha') fa = .true.
            if (ks(i) == 'beta')  fb = .true.
            if (ks(i) == 'GAMMA') fc = .true.
        end do
        call assert_true(fa, 'keys() contains "alpha"')
        call assert_true(fb, 'keys() contains "beta"')
        call assert_true(fc, 'keys() contains "GAMMA"')
        call h%destroy()
    end subroutine test_keys_returns_all_keys

    !---------------- keys_sorted() ----------------

    subroutine test_keys_sorted_returns_sorted_keys()
        type(vrefhash) :: h
        integer, target :: a, b, c
        type(string), allocatable :: ks(:)
        write(*,'(A)') 'test_keys_sorted_returns_sorted_keys'
        call h%init(32)
        a = 1; b = 2; c = 3
        ! Insert with mixed case/spaces; stored keys will be normalized by set_ref_*:
        !  - char keys: trim(adjustl())
        !  - string keys: adjustl(to_char()) (to_char trims trailing)
        call h%set_ref('gamma', a)
        call h%set_ref('  Alpha  ', b)
        call h%set_ref(string('beta'), c)
        ks = h%keys_sorted()
        call assert_int(3, size(ks), 'keys_sorted() size == 3')
        ! Expect case-insensitive lexicographic order: Alpha, beta, gamma
        call assert_string_eq('Alpha', ks(1), 'keys_sorted(1) == "Alpha"')
        call assert_string_eq('beta',  ks(2), 'keys_sorted(2) == "beta"')
        call assert_string_eq('gamma', ks(3), 'keys_sorted(3) == "gamma"')
        call h%destroy()
    end subroutine test_keys_sorted_returns_sorted_keys

end module simple_vrefhash_tester
