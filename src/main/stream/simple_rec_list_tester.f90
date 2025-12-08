module simple_rec_list_tester
use simple_rec_list
use simple_defs
use simple_test_utils
use simple_string
implicit none

public :: run_all_rec_list_tests
private

contains

    subroutine run_all_rec_list_tests()
        write(*,'(A)') '**** running rec_list tests ****'
        call test_default_and_size()
        call test_push_back_and_at()
        call test_iterator_traversal()
        call test_replace_at()
        call test_iterator_replace()
        call test_included_mask()
        call test_flag_included()
        call test_particle_sums()
        call test_append_operator()
        call test_slice_subset()
        call test_assign_semantics()
        call test_kill_behavior()
    end subroutine run_all_rec_list_tests

    !----------------------------------------------------------
    ! Default & size semantics
    !----------------------------------------------------------
    subroutine test_default_and_size()
        type(rec_list) :: lst
        write(*,'(A)') 'test_default_and_size'
        call assert_int(0, lst%size(), 'default list size = 0')
    end subroutine test_default_and_size

    !----------------------------------------------------------
    ! Push_back + at() + polymorphic safety
    !----------------------------------------------------------
    subroutine test_push_back_and_at()
        type(rec_list) :: lst
        type(project_rec) :: pr
        class(rec), allocatable :: r
        write(*,'(A)') 'test_push_back_and_at'
        pr%id = 11
        pr%projname = string("alpha")
        call lst%push_back(pr)
        call assert_int(1, lst%size(), 'size after push_back = 1')
        call lst%at_rec(1, r)
        call assert_true(allocated(r), 'result of at() is allocated polymorphic')
        select type(r)
        type is(project_rec)
            call assert_int(11, r%id, 'at() retrieves same ID')
            call assert_char('alpha', r%projname%to_char(), 'stored name ok')
        class default
            call assert_true(.false., 'wrong dynamic type returned from at()')
        end select
    end subroutine test_push_back_and_at

    !----------------------------------------------------------
    ! Iterator traversal
    !----------------------------------------------------------
    subroutine test_iterator_traversal()
        type(rec_list) :: lst
        type(project_rec) :: a,b
        type(rec_iterator) :: it
        integer :: count
        write(*,'(A)') 'test_iterator_traversal'
        a%id = 1; a%projname = string("A")
        b%id = 2; b%projname = string("B")
        call lst%push_back(a)
        call lst%push_back(b)
        it = lst%begin()
        count = 0
        do while (it%valid())
            count = count + 1
            call it%next()
        end do
        call assert_int(2, count, 'iterator should visit both elements')
    end subroutine test_iterator_traversal

    !----------------------------------------------------------
    ! replace_at()
    !----------------------------------------------------------
    subroutine test_replace_at()
        type(rec_list) :: lst
        type(project_rec) :: pr1, pr2
        class(rec), allocatable :: r
        write(*,'(A)') 'test_replace_at'
        pr1%id = 10
        pr2%id = 20
        call lst%push_back(pr1)
        call lst%replace_at(1, pr2)
        call lst%at_rec(1, r)
        select type(r)
        type is(project_rec)
            call assert_int(20, r%id, 'replace_at stored new value')
        class default
            call assert_true(.false., 'wrong type in replace_at result')
        end select
    end subroutine test_replace_at

    !----------------------------------------------------------
    ! Iterator-based replacement
    !----------------------------------------------------------
    subroutine test_iterator_replace()
        type(rec_list) :: lst
        type(process_rec) :: p1, p2
        type(rec_iterator) :: it
        class(rec), allocatable :: r
        write(*,'(A)') 'test_iterator_replace'
        p1%id = 100
        p2%id = 200
        call lst%push_back(p1)
        it = lst%begin()
        call it%replace(lst, p2)
        call lst%at_rec(1, r)
        select type(r)
        type is(process_rec)
            call assert_int(200, r%id, 'iterator replace updated record')
        class default
            call assert_true(.false., 'iterator replace returned wrong type')
        end select
    end subroutine test_iterator_replace

    subroutine test_included_mask()
        type(rec_list) :: lst
        type(project_rec) :: r(3)
        logical, allocatable :: m(:)
        write(*,'(A)') 'test_included_mask'
        r(1)%id=1; r(1)%included=.false.
        r(2)%id=2; r(2)%included=.true.
        r(3)%id=3; r(3)%included=.false.
        call lst%push_back(r(1))
        call lst%push_back(r(2))
        call lst%push_back(r(3))
        m = lst%get_included_flags()
        call assert_true(all(m .eqv. [.false., .true., .false.]), 'included() logical array correct')
    end subroutine test_included_mask

    subroutine test_flag_included()
        type(rec_list) :: lst
        type(project_rec) :: r(3)
        logical, allocatable :: m(:)
        write(*,'(A)') 'test_flag_included'
        r(1)%included=.false.; r(2)%included=.false.; r(3)%included=.false.
        call lst%push_back(r(1))
        call lst%push_back(r(2))
        call lst%push_back(r(3))
        call lst%set_included_flags([2,3])
        m = lst%get_included_flags()
        call assert_true(all(m .eqv. [.false., .true., .true.]), 'range flagging sets correct subset')
    end subroutine test_flag_included

    subroutine test_particle_sums()
        type(rec_list) :: lst
        type(project_rec) :: a,b
        integer :: sum1, sum2
        write(*,'(A)') 'test_particle_sums'
        a%nptcls=100; a%nptcls_sel=25; a%included=.true.
        b%nptcls=200; b%nptcls_sel=50; b%included=.false.
        call lst%push_back(a)
        call lst%push_back(b)
        sum1 = lst%get_nptcls_tot()
        sum2 = lst%get_nptcls_tot(.true.)
        call assert_int(300, sum1, 'sum included+excluded')
        call assert_int(200, sum2, 'sum of non-included only')
        call assert_int(75, lst%get_nptcls_sel_tot(), 'selected total')
    end subroutine test_particle_sums

    !----------------------------------------------------------
    ! Operator // append
    !----------------------------------------------------------
    subroutine test_append_operator()
        type(rec_list) :: a,b,c
        type(chunk_rec) :: r1,r2
        write(*,'(A)') 'test_append_operator'
        r1%id = 1
        r2%id = 2
        call a%push_back(r1)
        call b%push_back(r2)
        c = a // b
        call assert_int(2, c%size(), 'append creates combined list')
    end subroutine test_append_operator

    !----------------------------------------------------------
    ! slice(), subset_to()
    !----------------------------------------------------------
    subroutine test_slice_subset()
        type(rec_list) :: lst, sub1, sub2
        type(project_rec) :: r(4)
        integer :: i
        write(*,'(A)') 'test_slice_subset_copy'
        do i=1,4
            r(i)%id = i
            call lst%push_back(r(i))
        end do
        call lst%slice(2,3,sub1)
        call assert_int(2, sub1%size(), 'slice extracts correct size')
        call sub2%subset_to(lst, 2, 4)
        call assert_int(3, sub2%size(), 'subset_to extracts correct size')
    end subroutine test_slice_subset

    !----------------------------------------------------------
    ! Assignment semantics
    !----------------------------------------------------------
    subroutine test_assign_semantics()
        type(rec_list) :: a,b
        type(chunk_rec) :: r
        write(*,'(A)') 'test_assign_semantics'
        r%id = 99
        call a%push_back(r)
        b = a
        call assert_int(1, b%size(), 'assignment copies content')
    end subroutine test_assign_semantics

    !----------------------------------------------------------
    ! kill()
    !----------------------------------------------------------
    subroutine test_kill_behavior()
        type(rec_list) :: lst
        type(project_rec) :: r
        write(*,'(A)') 'test_kill_behavior'
        r%id = 4
        call lst%push_back(r)
        call lst%kill()
        call assert_int(0, lst%size(), 'kill empties the list')
    end subroutine test_kill_behavior

end module simple_rec_list_tester
