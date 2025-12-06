module simple_chash_tester
use simple_chash
use simple_string
use simple_defs
use simple_fileio
use simple_test_utils
implicit none
private
public :: run_all_chash_tests

contains

    subroutine run_all_chash_tests()
        write(*,'(A)') '**** running all chash tests ****'
        call test_construction_and_kill()
        call test_push_and_set_basic()
        call test_set_replaces_existing()
        call test_delete()
        call test_isthere_lookup_reverse()
        call test_get_by_key_and_index()
        call test_chash2str_and_strlen()
        call test_sort()
        ! call report_summary()
    end subroutine run_all_chash_tests

    !---------------- construction / destruction ----------------

    subroutine test_construction_and_kill()
        type(chash) :: h
        write(*,'(A)') 'test_construction_and_kill'
        h = chash()
        call assert_int(0, h%size_of(), 'new chash has size_of() == 0')
        call h%set('a', '1')
        call h%set('b', '2')
        call assert_int(2, h%size_of(), 'size_of() after 2 inserts')
        call h%kill()
        call assert_int(0, h%size_of(), 'size_of() after kill() == 0')
    end subroutine test_construction_and_kill

    !---------------- push / set semantics ----------------

    subroutine test_push_and_set_basic()
        type(chash) :: h
        write(*,'(A)') 'test_push_and_set_basic'
        h = chash()
        call h%push('alpha', '1')
        call h%push('beta',  string('2'))
        call h%push(string('gamma'), string('3'))
        call assert_int(3, h%size_of(), 'push_* increments size_of() to 3')
        call assert_true(h%isthere('alpha'), 'isthere(alpha)')
        call assert_true(h%isthere('beta' ), 'isthere(beta)')
        call assert_true(h%isthere('gamma'), 'isthere(gamma)')
        call assert_true(.not. h%isthere('delta'), 'isthere(delta) false')
        ! Check values via assertion helper
        call assert_string_eq('1', h%get('alpha'), 'get(alpha) == 1')
        call assert_string_eq('2', h%get('beta' ), 'get(beta) == 2')
        call assert_string_eq('3', h%get('gamma'), 'get(gamma) == 3')
    end subroutine test_push_and_set_basic

    subroutine test_set_replaces_existing()
        type(chash) :: h
        write(*,'(A)') 'test_set_replaces_existing'
        h = chash()
        call h%set('key', 'value1')
        call assert_int(1, h%size_of(), 'size_of() after first set(key)')
        call assert_string_eq('value1', h%get('key'), 'initial value for key')
        call h%set('key', 'value2')
        call assert_int(1, h%size_of(), 'size_of() unchanged after set(key) again')
        call assert_string_eq('value2', h%get('key'), 'value replaced for key')
        call h%set(string('key'), string('value3'))
        call assert_int(1, h%size_of(), 'size_of() still 1 after set(string,string)')
        call assert_string_eq('value3', h%get('key'), 'value3 after set(string,string)')
        call h%set('  key   ', 'value4')
        call assert_true(h%isthere('key'), 'isthere with trimmed key')
        call assert_string_eq('value4', h%get('key'), 'value updated when key has spaces')
    end subroutine test_set_replaces_existing

    !---------------- delete ----------------

    subroutine test_delete()
        type(chash) :: h
        write(*,'(A)') 'test_delete'
        h = chash()
        call h%set('k1','v1')
        call h%set('k2','v2')
        call h%set('k3','v3')
        call assert_int(3, h%size_of(), 'size_of() == 3 before delete')
        call h%delete('k2')
        call assert_int(2, h%size_of(), 'size_of() == 2 after delete(k2)')
        call assert_true(.not. h%isthere('k2'), 'k2 removed')
        call assert_true(h%isthere('k1'), 'k1 still present')
        call assert_true(h%isthere('k3'), 'k3 still present')
        call h%delete('k1')
        call assert_int(1, h%size_of(), 'size_of() == 1 after delete(k1)')
        call assert_true(.not. h%isthere('k1'), 'k1 removed')
        call assert_true(h%isthere('k3'), 'k3 remains')
        call h%delete('k3')
        call assert_int(0, h%size_of(), 'size_of() == 0 after delete(k3)')
        call assert_true(.not. h%isthere('k3'), 'k3 removed')
        call h%delete('does_not_exist')
        call assert_int(0, h%size_of(), 'delete of missing key leaves size_of() == 0')
    end subroutine test_delete

    !---------------- isthere / lookup / reverselookup ----------------

    subroutine test_isthere_lookup_reverse()
        type(chash) :: h
        integer     :: idx
        write(*,'(A)') 'test_isthere_lookup_reverse'
        h = chash()
        call h%set('a','1')
        call h%set('b','2')
        call h%set('c','1')
        call assert_true(h%isthere('a'), 'isthere(a)')
        call assert_true(.not. h%isthere('x'), 'isthere(x) false')
        idx = h%lookup('b')
        call assert_true(idx >= 1 .and. idx <= h%size_of(), 'lookup(b) returns valid index')
        idx = h%lookup('x')
        call assert_int(0, idx, 'lookup(x) returns 0 (not found)')
        idx = h%reverselookup('1')
        call assert_true(idx >= 1 .and. idx <= h%size_of(), 'reverselookup(1) returns some index')
        call assert_string_eq('1', h%get(idx), 'reverselookup(1) maps to value==1')
    end subroutine test_isthere_lookup_reverse

    !---------------- get(key) / get(index) / get_key(index) ----------------

    subroutine test_get_by_key_and_index()
        type(chash) :: h
        type(string) :: v, k
        write(*,'(A)') 'test_get_by_key_and_index'
        h = chash()
        call h%set('key1', 'val1')
        call h%set('key2', 'val2')
        v = h%get('key1')
        call assert_string_eq('val1', v, 'get(key1) == val1')
        v = h%get('key2')
        call assert_string_eq('val2', v, 'get(key2) == val2')
        v = h%get(2)
        call assert_string_eq('val2', v, 'get(2) == val2')
        k = h%get_key(1)
        call assert_string_eq('key1', k, 'get_key(1) == key1')
    end subroutine test_get_by_key_and_index

    !---------------- chash2str / chash_strlen ----------------

    subroutine test_chash2str_and_strlen()
        type(chash) :: h
        type(string) :: s
        integer :: n
        write(*,'(A)') 'test_chash2str_and_strlen'
        h = chash()
        n = h%chash_strlen()
        call assert_int(0, n, 'chash_strlen(empty) == 0')
        s = h%chash2str()
        call assert_true(s%to_char().eq.'', 'chash2str(empty) returns blank string')
        call h%set('x', '1')
        n = h%chash_strlen()
        call assert_int(len('x') + len('1') + 1, n, 'chash_strlen for 1 entry')
        s = h%chash2str()
        call assert_string_eq('x=1', s, 'chash2str single entry "x=1"')
        call h%set('y', '2')
        call h%set('z', '3')
        n = h%chash_strlen()
        call assert_int( &
            len('x')+len('1') + len('y')+len('2') + len('z')+len('3') + 3 + 2, &
            n, 'chash_strlen for 3 entries' )
        s = h%chash2str()
        call assert_string_eq('x=1 y=2 z=3', s, 'chash2str for 3 entries')
    end subroutine test_chash2str_and_strlen

    !---------------- sort ----------------

    subroutine test_sort()
        type(chash) :: h
        type(string) :: s
        write(*,'(A)') 'test_sort'
        h = chash()
        call h%set('gamma', '3')
        call h%set('alpha', '1')
        call h%set('beta',  '2')
        s = h%chash2str()
        call assert_string_eq('gamma=3 alpha=1 beta=2', s, 'insertion order before sort')
        call h%sort()
        s = h%chash2str()
        call assert_string_eq('alpha=1 beta=2 gamma=3', s, 'lexicographical order after sort')
        call assert_string_eq('1', h%get('alpha'), 'alpha still maps to 1')
        call assert_string_eq('2', h%get('beta' ), 'beta  still maps to 2')
        call assert_string_eq('3', h%get('gamma'), 'gamma still maps to 3')
    end subroutine test_sort

end module simple_chash_tester
