module simple_hash_tester
use simple_defs
use simple_hash
use simple_string
use simple_string_utils
use simple_test_utils
implicit none
private
public :: run_all_hash_tests

contains

    subroutine run_all_hash_tests
        write(*,'(A)') '**** running all hash tests ****'
        call test_constructor_and_kill()
        call test_push_and_set()
        call test_delete()
        call test_lookup_and_isthere()
        call test_getters_numeric()
        call test_get_keys_and_get_str()
        call test_hash2str_and_strlen()
        call test_copy_and_realloc()
        ! call report_summary()
    end subroutine run_all_hash_tests

    !---------------------------
    ! constructors / kill
    !---------------------------
    subroutine test_constructor_and_kill()
        type(hash) :: h
        ! constructor_1 via generic hash()
        h = hash()
        call assert_true(h%size_of() == 0, 'constructor_1: size_of == 0')
        call assert_int(0, h%size_of(), 'constructor_1: hash_index == 0')
        ! push something and then kill
        call h%set('a', 1.0_dp)
        call assert_int(1, h%size_of(), 'after set, size_of == 1')
        call h%kill()
        call assert_int(0, h%size_of(), 'kill resets size_of to 0')
    end subroutine test_constructor_and_kill

    !---------------------------
    ! push / set
    !---------------------------
    subroutine test_push_and_set()
        type(hash) :: h
        real(dp)   :: val
        h = hash(2)  ! constructor_2
        ! push_1: real(dp)
        call h%push('k1', 1.5_dp)
        call assert_int(1, h%size_of(), 'push_1: size_of == 1 after first push')
        val = h%get_dp('k1')
        call assert_double(1.5_dp, val, 'push_1: value for k1')
        ! push_2: integer
        call h%push('k2', 3)
        call assert_int(2, h%size_of(), 'push_2: size_of == 2 after second push')
        val = h%get_dp('k2')
        call assert_double(3.0_dp, val, 'push_2: value for k2')
        ! set_1: real
        call h%set('k1', 2.5_dp)
        val = h%get_dp('k1')
        call assert_double(2.5_dp, val, 'set_1: replace existing real')
        ! set_2: integer
        call h%set('k2', 10)
        val = h%get_dp('k2')
        call assert_double(10.0_dp, val, 'set_2: replace existing integer')
        ! set_3: real(dp)
        call h%set('k3', 4.25_dp)
        call assert_int(3, h%size_of(), 'set_3: new key increases size')
        val = h%get_dp('k3')
        call assert_double(4.25_dp, val, 'set_3: new dp value')
    end subroutine test_push_and_set

    !---------------------------
    ! delete
    !---------------------------
    subroutine test_delete()
        type(hash) :: h
        real(dp)   :: v
        h = hash()
        call h%set('a', 1.0_dp)
        call h%set('b', 2.0_dp)
        call h%set('c', 3.0_dp)
        call assert_int(3, h%size_of(), 'delete: initial size 3')
        ! delete middle element 'b'
        call h%delete('b')
        call assert_int(2, h%size_of(), 'delete: size reduced to 2')
        call assert_true(.not. h%isthere('b'), 'delete: key b removed')
        call assert_true(h%isthere('a'), 'delete: key a still present')
        call assert_true(h%isthere('c'), 'delete: key c still present')
        ! verify order shift: we don’t strictly depend on order semantics,
        ! but we can at least confirm values for remaining keys are ok
        v = h%get_dp('a')
        call assert_double(1.0_dp, v, 'delete: value for a after delete')
        v = h%get_dp('c')
        call assert_double(3.0_dp, v, 'delete: value for c after delete')
        ! deleting non-existent key is a no-op
        call h%delete('not_there')
        call assert_int(2, h%size_of(), 'delete: deleting non-existent key leaves size unchanged')
    end subroutine test_delete

    !---------------------------
    ! lookup / isthere
    !---------------------------
    subroutine test_lookup_and_isthere()
        type(hash) :: h
        integer    :: idx
        h = hash()
        call h%set('x', 10.0_dp)
        call h%set('y', 20.0_dp)
        call assert_true(h%isthere('x'), 'isthere: x present')
        call assert_true(h%isthere('y'), 'isthere: y present')
        call assert_true(.not. h%isthere('z'), 'isthere: z not present')
        idx = h%lookup('x')
        call assert_true(idx > 0, 'lookup: index for x > 0')
        idx = h%lookup('y')
        call assert_true(idx > 0, 'lookup: index for y > 0')
        idx = h%lookup('z')
        call assert_int(0, idx, 'lookup: index for z == 0')
    end subroutine test_lookup_and_isthere

    !---------------------------
    ! getters for numeric values
    !---------------------------
    subroutine test_getters_numeric()
        type(hash) :: h
        real(dp)   :: vdp
        real       :: vsp
        integer    :: ival
        h = hash()
        call h%set('a', 1.25_dp)
        call h%set('b', 5)        ! integer
        call h%set('c', 3.75_dp)
        ! get() returns real(sp)
        vsp = h%get('a')
        call assert_double(1.25_dp, real(vsp,dp), 'get (sp): a')
        ! get_dp()
        vdp = h%get_dp('c')
        call assert_double(3.75_dp, vdp, 'get_dp: c')
        ! getter_1: real(sp)
        vsp = 0.0
        call h%getter('b', vsp)
        call assert_double(5.0_dp, real(vsp,dp), 'getter_1 (sp): b')
        ! getter_2: real(dp)
        vdp = 0.0_dp
        call h%getter('b', vdp)
        call assert_double(5.0_dp, vdp, 'getter_2 (dp): b')
        ! getter_3: integer
        ival = 0
        call h%getter('b', ival)
        call assert_int(5, ival, 'getter_3 (int): b')
        ! get_value_at
        vdp = h%get_value_at(h%lookup('a'))
        call assert_double(1.25_dp, vdp, 'get_value_at: a')
    end subroutine test_getters_numeric

    !---------------------------
    ! get_keys / get_str
    !---------------------------
    subroutine test_get_keys_and_get_str()
        type(hash)   :: h
        type(string), allocatable :: ks(:)
        integer      :: n, idx
        type(string) :: s
        h = hash()
        call h%set('k1', 1.0_dp)
        call h%set('k2', 2.0_dp)
        call h%set('k3', 3.0_dp)
        ks = h%get_keys()
        n  = size(ks)
        call assert_int(3, n, 'get_keys: size 3')
        ! just check presence; order is implementation detail but expected to match insertion
        call assert_string_eq('k1', ks(1), 'get_keys: first key')
        call assert_string_eq('k2', ks(2), 'get_keys: second key')
        call assert_string_eq('k3', ks(3), 'get_keys: third key')
        ! get_str(keyindx)
        idx = h%lookup('k2')
        s = h%get_key(idx)
        call assert_string_eq('k2', s, 'get_key: key for index of k2')
        if (allocated(ks)) then
            call ks%kill
            deallocate(ks)
        end if
    end subroutine test_get_keys_and_get_str

    !---------------------------
    ! hash2str / hash_strlen
    !---------------------------
    subroutine test_hash2str_and_strlen()
        type(hash)   :: h
        type(string) :: sh
        integer      :: lenh
        h = hash()
        call h%set('a', 1.0_dp)
        call h%set('b', 2.0_dp)
        sh   = h%hash2str()
        lenh = h%hash_strlen()
        ! We mainly test that strlen is consistent with the actual string length
        call assert_int(len_trim(sh%to_char()), lenh, 'hash_strlen matches hash2str length')
        ! Very basic content check (ordering-dependent)
        ! Example expected pattern: "a=1 b=2" (but values are formatted by real2str)
        call assert_true(index(sh%to_char(), 'a=') > 0, 'hash2str contains a=')
        call assert_true(index(sh%to_char(), 'b=') > 0, 'hash2str contains b=')
        call sh%kill
    end subroutine test_hash2str_and_strlen

    !---------------------------
    ! copy and realloc
    !---------------------------
    subroutine test_copy_and_realloc()
        type(hash) :: h1, h2
        real(dp)   :: v
        h1 = hash(1)
        call h1%set('x', 10.0_dp)
        call h1%set('y', 20.0_dp)
        ! This should force realloc internally
        call h1%set('z', 30.0_dp)
        ! copy via assignment
        h2 = h1
        call assert_int(h1%size_of(), h2%size_of(), 'copy: size_of preserved')
        v = h2%get_dp('x')
        call assert_double(10.0_dp, v, 'copy: x')
        v = h2%get_dp('y')
        call assert_double(20.0_dp, v, 'copy: y')
        v = h2%get_dp('z')
        call assert_double(30.0_dp, v, 'copy: z')
        ! mutating h2 must not affect h1 if copy is deep
        call h2%set('x', 11.0_dp)
        v = h2%get_dp('x')
        call assert_double(11.0_dp, v, 'copy: h2 mutated x')
        v = h1%get_dp('x')
        call assert_double(10.0_dp, v, 'copy: h1 x unchanged after h2 mutation')
    end subroutine test_copy_and_realloc

end module simple_hash_tester
