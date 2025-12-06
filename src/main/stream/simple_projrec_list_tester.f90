module simple_projrec_list_tester
use simple_projrec_list, only: projrec_list, projrec
use simple_string,       only: string
use simple_test_utils
implicit none

public :: run_all_projrec_list_tests
private

contains

    subroutine run_all_projrec_list_tests
        write(*,'(A)') '**** running all project_list tests ****'
        call test_constructor_and_new()
        call test_is_allocated_and_len()
        call test_getters_on_unallocated_and_zero_length()
        call test_push_and_set()
        call test_copy_from()
        call test_replace()
        call test_kill()
        call test_assignment()
        call test_append()
        call test_is_project()
        call test_totals()
        ! call report_summary()
    end subroutine run_all_projrec_list_tests

    ! Helper: int -> char for project names
    function to_char_int(i) result(ch)
        integer, intent(in) :: i
        character(len=16)   :: ch
        write(ch,'(I0)') i
    end function to_char_int

    ! Helper: build a record with predictable contents
    function make_projrec(i) result(rec)
        integer, intent(in) :: i
        type(projrec)       :: rec
        rec%projname   = 'proj_' // trim(adjustl(to_char_int(i)))
        rec%micind     = i
        rec%nptcls     = 10 * i
        rec%nptcls_sel = 5 * i
        rec%included   = (mod(i,2) == 0)
    end function make_projrec

    ! Helper: fill list with n records using set (grows via set/push)
    subroutine build_list(lst, n)
        type(projrec_list), intent(inout) :: lst
        integer,            intent(in)    :: n
        integer :: i
        type(projrec)      :: r
        call lst%kill()  ! start from unallocated
        do i = 1, n
            r = make_projrec(i)
            call lst%set(i, r)
        end do
    end subroutine build_list

    ! -----------------------------
    ! 1. constructor and new
    ! -----------------------------
    subroutine test_constructor_and_new()
        type(projrec_list) :: a, b
        integer :: n
        write(*,'(A)') 'test_constructor_and_new'
        ! constructor with explicit n
        n = 3
        a = projrec_list(n)
        call assert_true(a%is_allocated(), 'constructor(n): list should be allocated')
        call assert_int(n, a%len(),      'constructor(n): size should equal n')
        ! constructor without n -> zero-length but allocated
        b = projrec_list()
        call assert_true(b%is_allocated(), 'constructor(): list should be allocated (len 0)')
        call assert_int(0, b%len(),       'constructor(): size should be 0')
        ! new with explicit n
        call a%new(5)
        call assert_true(a%is_allocated(), 'new(n): list should be allocated')
        call assert_int(5, a%len(),       'new(n): size should equal n')
        ! new without n -> zero-length but allocated
        call a%new()
        call assert_true(a%is_allocated(), 'new(): list should be allocated')
        call assert_int(0, a%len(),       'new(): size should be 0')
    end subroutine test_constructor_and_new

    ! -----------------------------
    ! 2. is_allocated and size
    ! -----------------------------
    subroutine test_is_allocated_and_len()
        type(projrec_list) :: a
        write(*,'(A)') 'test_is_allocated_and_size'
        a = projrec_list()
        call assert_true(a%is_allocated(), 'is_allocated: freshly constructed list should be allocated')
        call assert_int(0, a%len(),       'size: freshly constructed empty list should be 0')
        call a%kill()
        call assert_true(.not. a%is_allocated(), 'is_allocated: after kill() should be .false.')
        call assert_int(0, a%len(),             'size: after kill() should be 0')
    end subroutine test_is_allocated_and_len

    ! ------------------------------------------------------------
    ! 3. getters on unallocated and zero-length list (edge cases)
    ! ------------------------------------------------------------
    subroutine test_getters_on_unallocated_and_zero_length()
        type(projrec_list) :: a
        type(string), allocatable :: pnames(:)
        integer, allocatable :: i_arr(:)
        integer :: val
        write(*,'(A)') 'test_getters_on_unallocated_and_zero_length'
        ! Start with unallocated
        a = projrec_list()
        call a%kill()
        call assert_true(.not. a%is_allocated(), 'getters: a should be unallocated before tests')
        call assert_string_eq('', a%get_projname(1), 'get_projname on unallocated should be empty')
        val = a%get_micind(1)
        call assert_int(0, val,                     'get_micind on unallocated should be 0')
        val = a%get_nptcls(1)
        call assert_int(0, val,                     'get_nptcls on unallocated should be 0')
        val = a%get_nptcls_sel(1)
        call assert_int(0, val,                     'get_nptcls_sel on unallocated should be 0')
        call assert_int(0, a%get_nptcls_tot(),      'get_nptcls_tot on unallocated should be 0')
        call assert_int(0, a%get_nptcls_sel_tot(),  'get_nptcls_sel_tot on unallocated should be 0')
        pnames = a%get_projname_arr()
        call assert_int(0, size(pnames), 'get_projname_arr on unallocated should be zero-length')
        i_arr = a%get_micind_arr()
        call assert_int(0, size(i_arr),  'get_micind_arr on unallocated should be zero-length')
        i_arr = a%get_nptcls_arr()
        call assert_int(0, size(i_arr),  'get_nptcls_arr on unallocated should be zero-length')
        i_arr = a%get_nptcls_sel_arr()
        call assert_int(0, size(i_arr),  'get_nptcls_sel_arr on unallocated should be zero-length')
        ! Now zero-length but allocated
        a = projrec_list()
        call assert_true(a%is_allocated(), 'zero-length allocated: is_allocated() should be true')
        call assert_int(0, a%len(),       'zero-length allocated: len() should be 0')
        pnames = a%get_projname_arr()
        call assert_int(0, size(pnames), 'get_projname_arr on zero-length should be zero-length')
        i_arr = a%get_micind_arr()
        call assert_int(0, size(i_arr),  'get_micind_arr on zero-length should be zero-length')
        i_arr = a%get_nptcls_arr()
        call assert_int(0, size(i_arr),  'get_nptcls_arr on zero-length should be zero-length')
        i_arr = a%get_nptcls_sel_arr()
        call assert_int(0, size(i_arr),  'get_nptcls_sel_arr on zero-length should be zero-length')
        call assert_int(0, a%get_nptcls_tot(),     'get_nptcls_tot on zero-length should be 0')
        call assert_int(0, a%get_nptcls_sel_tot(), 'get_nptcls_sel_tot on zero-length should be 0')
        call assert_string_eq('', a%get_projname(1), 'get_projname on zero-length should be empty')
        call assert_int(0, a%get_micind(1),          'get_micind on zero-length should be 0')
        call assert_int(0, a%get_nptcls(1),          'get_nptcls on zero-length should be 0')
        call assert_int(0, a%get_nptcls_sel(1),      'get_nptcls_sel on zero-length should be 0')
    end subroutine test_getters_on_unallocated_and_zero_length

    ! ------------------------------------------------
    ! 4. push and set (specifically exercising them)
    ! ------------------------------------------------
    subroutine test_push_and_set()
        type(projrec_list) :: a
        type(projrec)      :: r
        integer :: i
        write(*,'(A)') 'test_push_and_set'
        call a%kill()
        ! push into unallocated list: expect size 1
        r = make_projrec(1)
        call a%push(r)
        call assert_int(1, a%len(), 'push: size should be 1 after first push')
        call assert_string_eq('proj_1', a%get_projname(1), 'push: proj_1 at index 1')
        ! set at index size+1 should append
        r = make_projrec(2)
        call a%set(2, r)
        call assert_int(2, a%len(), 'set(size+1): should append to size 2')
        call assert_string_eq('proj_2', a%get_projname(2), 'set(size+1): proj_2 at index 2')
        ! overwrite middle element
        r = make_projrec(99)
        r%projname = 'override'
        call a%set(1, r)
        call assert_int(2, a%len(), 'set(existing): size unchanged')
        call assert_string_eq('override', a%get_projname(1), 'set(existing): projname overwritten')
        ! grow further using set i = size+1 multiple times
        do i = 3, 5
            r = make_projrec(i)
            call a%set(i,r)
        end do
        call assert_int(5, a%len(), 'set(size+1): repeated growth to 5')
        call assert_string_eq('proj_5', a%get_projname(5), 'set(size+1): proj_5 at index 5')
        ! NOTE: we intentionally do NOT test out-of-bounds error path
        ! because THROW_HARD is expected to abort the program.
    end subroutine test_push_and_set

    ! -----------------------------
    ! 5. copy_from
    ! -----------------------------
    subroutine test_copy_from()
        type(projrec_list) :: src, dst
        integer       :: i
        type(projrec) :: r
        write(*,'(A)') 'test_copy_from'
        ! Case 1: dst larger than src
        call build_list(src, 3)
        call dst%kill()
        ! make dst size 5 with known pattern "dst_i"
        do i = 1, 5
            r = make_projrec(-i)
            call dst%set(i, r)   ! initial dummy with negative indices
            r = make_projrec(100 + i)
            call dst%set(i, r)   ! overwrite with distinct pattern
            r = make_projrec(100 + i)
            call dst%set(i, r)   ! (still, values irrelevant; only size matters)
        end do
        call dst%copy_from(src)
        call assert_int(5, dst%len(), 'copy_from: dst size should remain 5 when src smaller')
        ! first 3 elements match src
        do i = 1, 3
            call assert_string_eq('proj_'//trim(adjustl(to_char_int(i))), dst%get_projname(i), &
                                  'copy_from: projname copied correctly (dst larger)')
            call assert_int( i,      dst%get_micind(i),     'copy_from: micind copied correctly (dst larger)')
            call assert_int(10*i,    dst%get_nptcls(i),     'copy_from: nptcls copied correctly (dst larger)')
            call assert_int(5*i,     dst%get_nptcls_sel(i), 'copy_from: nptcls_sel copied correctly (dst larger)')
        end do
        ! Case 2: src larger -> dst grows
        call build_list(src, 5)
        call build_list(dst, 3)     ! smaller dst
        call dst%copy_from(src)
        call assert_int(5, dst%len(), 'copy_from: dst size should grow to src size when src larger')
        do i = 1, 5
            call assert_string_eq('proj_'//trim(adjustl(to_char_int(i))), dst%get_projname(i), &
                                  'copy_from: projname copied correctly (src larger)')
            call assert_int( i,      dst%get_micind(i),     'copy_from: micind copied correctly (src larger)')
            call assert_int(10*i,    dst%get_nptcls(i),     'copy_from: nptcls copied correctly (src larger)')
            call assert_int(5*i,     dst%get_nptcls_sel(i), 'copy_from: nptcls_sel copied correctly (src larger)')
        end do
        ! Case 3: unallocated src => no-op
        call src%kill()
        call assert_true(.not. src%is_allocated(), 'copy_from: src killed (unallocated)')
        call dst%copy_from(src)
        call assert_true(dst%is_allocated(), 'copy_from: dst remains allocated when src unallocated')
        call assert_int(5, dst%len(),       'copy_from: dst size unchanged when src unallocated')
    end subroutine test_copy_from

    ! -----------------------------
    ! 6. replace (move semantics)
    ! -----------------------------
    subroutine test_replace()
        type(projrec_list) :: src, dst
        write(*,'(A)') 'test_replace'
        call build_list(src, 4)
        call build_list(dst, 1)
        call dst%replace(src)
        call assert_int(4, dst%len(), 'replace: dst size should become src size after move_alloc')
        call assert_true(.not. src%is_allocated(), 'replace: src should be unallocated after move_alloc')
        call assert_string_eq('proj_1', dst%get_projname(1), 'replace: proj_1 present')
        call assert_string_eq('proj_4', dst%get_projname(4), 'replace: proj_4 present')
        call assert_int(4, dst%get_micind(4),                'replace: micind(4) correct')
        call src%kill()
        call dst%replace(src) ! no-op
        call assert_int(4, dst%len(), 'replace: dst unchanged when src unallocated')
    end subroutine test_replace

    ! -----------------------------
    ! 7. kill
    ! -----------------------------
    subroutine test_kill()
        type(projrec_list) :: a
        write(*,'(A)') 'test_kill'
        call build_list(a, 3)
        call assert_true(a%is_allocated(), 'kill: list allocated before kill')
        call a%kill()
        call assert_true(.not. a%is_allocated(), 'kill: list unallocated after kill')
        call assert_int(0, a%len(),             'kill: len() is 0 after kill')
        call a%kill()
        call assert_true(.not. a%is_allocated(), 'kill: repeated kill keeps unallocated state')
        call assert_int(0, a%len(),             'kill: repeated kill keeps size 0')
    end subroutine test_kill

    ! -----------------------------
    ! 8. assignment (=)
    ! -----------------------------
    subroutine test_assignment()
        type(projrec_list) :: a, b
        type(string) :: pnam1, pnam2
        integer      :: i
        write(*,'(A)') 'test_assignment'
        call build_list(a, 3)
        b = projrec_list()
        b = a
        call assert_true(b%is_allocated(), 'assignment: lhs allocated after assignment')
        call assert_int(a%len(), b%len(),  'assignment: sizes should match')
        do i = 1, a%len()
            pnam1 = a%get_projname(i)
            pnam2 = b%get_projname(i)
            call assert_string_eq(pnam1%to_char(), pnam2,               'assignment: projname copy')
            call assert_int(a%get_micind(i),       b%get_micind(i),     'assignment: micind copy')
            call assert_int(a%get_nptcls(i),       b%get_nptcls(i),     'assignment: nptcls copy')
            call assert_int(a%get_nptcls_sel(i),   b%get_nptcls_sel(i), 'assignment: nptcls_sel copy')
        end do
        call a%kill()
        b = a
        call assert_true(.not. b%is_allocated(), 'assignment: lhs unallocated when rhs unallocated')
        call assert_int(0, b%len(),              'assignment: size 0 when rhs unallocated')
    end subroutine test_assignment

    ! -----------------------------
    ! 9. append (// operator)
    ! -----------------------------
    subroutine test_append()
        type(projrec_list) :: a, b, c
        integer :: i
        write(*,'(A)') 'test_append'
        call build_list(a, 2)
        call build_list(b, 3)
        c = a // b
        call assert_int(5, c%len(), 'append: size(a//b) should be 5')
        do i = 1, 2
            call assert_string_eq('proj_'//trim(adjustl(to_char_int(i))), c%get_projname(i), &
                                  'append: projname a-part correct')
        end do
        do i = 1, 3
            call assert_string_eq('proj_'//trim(adjustl(to_char_int(i))), c%get_projname(2+i), &
                                  'append: projname b-part correct')
        end do
        b = projrec_list()
        c = a // b
        call assert_int(a%len(), c%len(), 'append: lhs non-empty, rhs empty => size lhs')
        a = projrec_list()
        call build_list(b, 4)
        c = a // b
        call assert_int(b%len(), c%len(), 'append: lhs empty, rhs non-empty => size rhs')
        a = projrec_list()
        b = projrec_list()
        c = a // b
        call assert_int(0, c%len(), 'append: both empty => size 0')
        call assert_true(c%is_allocated(), 'append: result allocated even when empty')
    end subroutine test_append

    ! -----------------------------
    ! 10. is_project
    ! -----------------------------
    subroutine test_is_project()
        type(projrec_list) :: a
        logical :: res
        write(*,'(A)') 'test_is_project'
        call build_list(a, 3)
        res = a%is_project(1, 'proj_1')
        call assert_true(res, 'is_project: valid index, matching name => .true.')
        res = a%is_project(2, 'proj_1')
        call assert_true(.not. res, 'is_project: valid index, non-matching name => .false.')
        res = a%is_project(0, 'proj_1')
        call assert_true(.not. res, 'is_project: index 0 => .false.')
        res = a%is_project(-1, 'proj_1')
        call assert_true(.not. res, 'is_project: negative index => .false.')
        res = a%is_project(4, 'proj_4')
        call assert_true(.not. res, 'is_project: index > size => .false.')
        call a%kill()
        res = a%is_project(1, 'proj_1')
        call assert_true(.not. res, 'is_project: on unallocated => .false.')
    end subroutine test_is_project

    ! -----------------------------
    ! 11. total counters
    ! -----------------------------
    subroutine test_totals()
        type(projrec_list) :: a
        integer :: expected_n, expected_sel, i, n
        write(*,'(A)') 'test_totals'
        a = projrec_list()
        call a%kill()
        call assert_int(0, a%get_nptcls_tot(),     'totals: unallocated get_nptcls_tot==0')
        call assert_int(0, a%get_nptcls_sel_tot(), 'totals: unallocated get_nptcls_sel_tot==0')
        a = projrec_list()
        call assert_int(0, a%get_nptcls_tot(),     'totals: zero-length get_nptcls_tot==0')
        call assert_int(0, a%get_nptcls_sel_tot(), 'totals: zero-length get_nptcls_sel_tot==0')
        n = 4
        call build_list(a, n)
        expected_n   = 0
        expected_sel = 0
        do i = 1, n
            expected_n   = expected_n   + 10 * i
            expected_sel = expected_sel +  5 * i
        end do
        call assert_int(expected_n,   a%get_nptcls_tot(),     'totals: get_nptcls_tot correct')
        call assert_int(expected_sel, a%get_nptcls_sel_tot(), 'totals: get_nptcls_sel_tot correct')
    end subroutine test_totals

end module simple_projrec_list_tester
