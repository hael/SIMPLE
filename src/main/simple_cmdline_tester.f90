module simple_cmdline_tester
use simple_defs
use simple_cmdline
use simple_chash
use simple_string
use simple_string_utils
use simple_test_utils
implicit none
private
public :: run_all_cmdline_tests

contains

    subroutine run_all_cmdline_tests()
        write(*,'(A)') '**** running all cmdline tests ****'
        call test_set_and_get_numeric()
        call test_set_and_get_char_and_string()
        call test_delete_behavior()
        call test_defined_and_lookup()
        call test_get_keys()
        call test_copy_and_assignment()
        call test_read_parsing()
        call test_gen_job_descr()
        ! call report_summary()
    end subroutine run_all_cmdline_tests

    !-----------------------------------------
    ! 1. basic setters / numeric getters
    !-----------------------------------------
    subroutine test_set_and_get_numeric()
        type(cmdline) :: cl
        real(dp)      :: r
        integer       :: i
        write(*,'(A)') 'test_set_and_get_numeric'
        call cl%kill()   ! start from clean state
        ! set_1: real(dp)
        call cl%set('alpha', 1.5_dp)
        call assert_true(cl%defined('alpha'), 'set_1: alpha defined')
        r = cl%get_rarg('alpha')
        call assert_double(1.5_dp, real(r,dp), 'set_1 / get_rarg: alpha')
        ! set_2: real (sp)
        call cl%set('beta', real(2.0, kind=kind(0.0)))
        r = cl%get_rarg('beta')
        call assert_double(2.0_dp, real(r,dp), 'set_2 / get_rarg: beta')
        ! set_4: integer as real
        call cl%set('box', 10)
        i = cl%get_iarg('box')
        call assert_int(10, i, 'set_4 / get_iarg: box')
        ! overwrite existing numeric key
        call cl%set('box', 20)
        i = cl%get_iarg('box')
        call assert_int(20, i, 'set overwrite numeric: box')
        ! ensure size (argcnt) tracked properly
        call assert_int(3, cl%get_argcnt(), 'argcnt after numeric insertions')
    end subroutine test_set_and_get_numeric

    !-----------------------------------------
    ! 2. char / string setters & get_carg
    !-----------------------------------------
    subroutine test_set_and_get_char_and_string()
        type(cmdline) :: cl
        type(string)  :: sval
        write(*,'(A)') 'test_set_and_get_char_and_string'
        call cl%kill()
        ! set_3: char
        call cl%set('dir', 'hello')
        call assert_true(cl%defined('dir'), 'set_3: dir defined')
        sval = cl%get_carg('dir')
        call assert_string_eq('hello', sval, 'get_carg: dir')
        ! set_5: (key,char) -> string value
        call cl%set('fname', string('world'))
        sval = cl%get_carg('fname')
        call assert_string_eq('world', sval, 'set_5 / get_carg: fname')
        ! set_6: (string,string)
        call cl%set(string('infile'), string('xy'))
        sval = cl%get_carg('infile')
        call assert_string_eq('xy', sval, 'set_6 / get_carg: infile')
        ! overwrite character key
        call cl%set('dir', 'hello2')
        sval = cl%get_carg('dir')
        call assert_string_eq('hello2', sval, 'overwrite char value: dir')
    end subroutine test_set_and_get_char_and_string

    !-----------------------------------------
    ! 3. delete logic
    !-----------------------------------------
    subroutine test_delete_behavior()
        type(cmdline) :: cl
        write(*,'(A)') 'test_delete_behavior'
        call cl%kill()
        call cl%set('a', 1.0_dp)
        call cl%set('b', 'foo')
        call cl%set('c', 3)
        call assert_int(3, cl%get_argcnt(), 'delete: initial argcnt 3')
        ! delete middle element
        call cl%delete('b')
        call assert_int(2, cl%get_argcnt(), 'delete: argcnt after deleting b')
        call assert_true(.not. cl%defined('b'), 'delete: b no longer defined')
        ! a and c must still be defined
        call assert_true(cl%defined('a'),   'delete: a still defined')
        call assert_true(cl%defined('c'),   'delete: c still defined')
        ! deleting non-existent key is a no-op
        call cl%delete('not_there')
        call assert_int(2, cl%get_argcnt(), 'delete: argcnt unchanged when deleting non-existent key')
    end subroutine test_delete_behavior

    !-----------------------------------------
    ! 4. defined / lookup
    !-----------------------------------------
    subroutine test_defined_and_lookup()
        type(cmdline) :: cl
        integer       :: idx
        write(*,'(A)') 'test_defined_and_lookup'
        call cl%kill()
        call cl%set('x', 1.0_dp)
        call cl%set('y', 'abc')
        call assert_true(cl%defined('x'), 'defined: x')
        call assert_true(cl%defined('y'), 'defined: y')
        call assert_true(.not. cl%defined('z'), 'defined: z not present')
        idx = cl%lookup('x')
        call assert_true(idx > 0, 'lookup: x index > 0')
        idx = cl%lookup('y')
        call assert_true(idx > 0, 'lookup: y index > 0')
        idx = cl%lookup('z')
        call assert_int(0, idx,   'lookup: z index == 0')
    end subroutine test_defined_and_lookup

    !-----------------------------------------
    ! 5. get_keys
    !-----------------------------------------
    subroutine test_get_keys()
        type(cmdline)   :: cl
        type(string), allocatable :: ks(:)
        write(*,'(A)') 'test_get_keys'
        call cl%kill()
        call cl%set('k1', 1.0_dp)
        call cl%set('k2', 2.0_dp)
        call cl%set('k3', 'three')
        ks = cl%get_keys()
        call assert_true(allocated(ks),    'get_keys: allocated result')
        call assert_int(3, size(ks),       'get_keys: size == argcnt')
        call assert_string_eq('k1', ks(1), 'get_keys: first key')
        call assert_string_eq('k2', ks(2), 'get_keys: second key')
        call assert_string_eq('k3', ks(3), 'get_keys: third key')
        call ks%kill
        deallocate(ks)
    end subroutine test_get_keys

    !-----------------------------------------
    ! 6. copy / assignment
    !-----------------------------------------
    subroutine test_copy_and_assignment()
        type(cmdline) :: cl1, cl2
        real          :: r
        type(string)  :: sval
        write(*,'(A)') 'test_copy_and_assignment'
        call cl1%kill()
        call cl1%set('a', 1.0_dp)
        call cl1%set('b', 'foo')
        cl2 = cl1   ! uses assignment(=) => assign => copy
        call assert_int(cl1%get_argcnt(), cl2%get_argcnt(), 'copy: argcnt preserved')
        r = cl2%get_rarg('a')
        call assert_double(1.0_dp, real(r,dp),              'copy: value a')
        sval = cl2%get_carg('b')
        call assert_string_eq('foo', sval,                  'copy: value b')
        ! Mutate cl2 and ensure cl1 is unchanged (deep copy semantics)
        call cl2%set('a', 2.0_dp)
        r = cl2%get_rarg('a')
        call assert_double(2.0_dp, real(r,dp),              'copy: cl2 a mutated')
        r = cl1%get_rarg('a')
        call assert_double(1.0_dp, real(r,dp),              'copy: cl1 a unchanged')
    end subroutine test_copy_and_assignment

    !-----------------------------------------
    ! 7. read() – parsing from a single string
    !-----------------------------------------
    subroutine test_read_parsing()
        type(cmdline) :: cl
        real          :: r
        integer       :: i
        type(string)  :: sval
        write(*,'(A)') 'test_read_parsing'
        call cl%kill()
        ! This feeds a synthetic command line via read(), bypassing real argv
        call cl%read('alpha=1.25 box=3 dir=hello')
        call assert_true(cl%defined('alpha'), 'read: alpha defined')
        call assert_true(cl%defined('box'), 'read: box defined')
        call assert_true(cl%defined('dir'), 'read: dir defined')
        r = cl%get_rarg('alpha')
        call assert_double(1.25_dp, real(r,dp), 'read / get_rarg: alpha')
        i = cl%get_iarg('box')
        call assert_int(3, i,                   'read / get_iarg: box')
        sval = cl%get_carg('dir')
        call assert_string_eq('hello', sval,    'read / get_carg: dir')
    end subroutine test_read_parsing

    !-----------------------------------------
    ! 8. gen_job_descr -> chash
    !-----------------------------------------
    subroutine test_gen_job_descr()
        type(cmdline) :: cl
        type(chash)   :: h
        type(string)  :: v, s
        write(*,'(A)') 'test_gen_job_descr'
        call cl%kill()
        call cl%set('a', 1.0_dp)
        call cl%set('b', 'foo')
        call cl%gen_job_descr(h, prg=string('myprog'))
        ! hash should contain 'a', 'b', and 'prg'
        call assert_true(h%isthere('a'),   'gen_job_descr: a in chash')
        call assert_true(h%isthere('b'),   'gen_job_descr: b in chash')
        call assert_true(h%isthere('prg'), 'gen_job_descr: prg in chash')
        v = h%get('a')
        s = realdp2str(1.0_dp)
        call assert_string_eq(s%to_char(), v,    'gen_job_descr: a value (real2str)')  ! exact text depends on real2str
        call v%kill
        v = h%get('b')
        call assert_string_eq('foo', v,    'gen_job_descr: b value')
        call v%kill
        v = h%get('prg')
        call assert_string_eq('myprog', v, 'gen_job_descr: prg value')
        call v%kill
        call h%kill
        call s%kill
    end subroutine test_gen_job_descr

end module simple_cmdline_tester
