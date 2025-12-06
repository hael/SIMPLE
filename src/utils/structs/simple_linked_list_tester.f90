module simple_linked_list_tester
use simple_linked_list
use simple_test_utils
implicit none

private
public :: run_all_list_tests

contains

    subroutine run_all_list_tests()
        write(*,'(A)') '**** running linked_list tests ****'
        call test_push_and_size()
        call test_pop_front()
        call test_front_back_at()
        call test_assign_and_copy_semantics()
        call test_append()
        call test_slice()
        call test_iteration()
        call test_iterator_index_and_advance()
        call test_end_iterator_behavior()
        call test_kill()
    end subroutine run_all_list_tests

    !------------------- push, size, empty -------------------

    subroutine test_push_and_size()
        type(linked_list) :: lst
        integer :: i
        write(*,'(A)') 'test_push_and_size'
        call assert_true(lst%is_empty(), 'initial list empty')
        call assert_int(0, lst%size(), 'initial size=0')
        do i = 1,5
            call lst%push_back(i)
        end do
        call assert_false(lst%is_empty(), 'list not empty after push')
        call assert_int(5, lst%size(), 'size after 5 pushes')
        call lst%kill()
    end subroutine test_push_and_size

    !------------------- pop_front -------------------

    subroutine test_pop_front()
        type(linked_list) :: lst
        class(*), allocatable :: x
        logical :: ok
        write(*,'(A)') 'test_pop_front'
        call lst%push_back(10)
        call lst%push_back(20)
        call lst%push_back(30)
        call lst%pop_front(x, ok)
        call assert_true(ok, 'pop_front success #1')
        call assert_int(10, transfer(x,0), 'pop value = 10')
        call lst%pop_front(x, ok)
        call assert_int(20, transfer(x,0), 'pop value = 20')
        call lst%pop_front(x, ok)
        call assert_int(30, transfer(x,0), 'pop value = 30')
        call lst%pop_front(x, ok)
        call assert_false(ok, 'pop empty returns false')
        call assert_true(lst%is_empty(), 'list empty after all pops')
        call lst%kill()
    end subroutine test_pop_front

    !------------------- accessors: front/back/at -------------------

    subroutine test_front_back_at()
        type(linked_list) :: lst
        class(*), allocatable :: x
        write(*,'(A)') 'test_front_back_at'
        call lst%push_back(5)
        call lst%push_back(10)
        call lst%push_back(15)
        call lst%front(x); call assert_int(5,  transfer(x,0), 'front = 5')
        call lst%back(x);  call assert_int(15, transfer(x,0), 'back = 15')
        call lst%at(2,x);  call assert_int(10, transfer(x,0), 'at(2) = 10')
        call lst%kill()
    end subroutine test_front_back_at

    !------------------- assign / deep copy -------------------

    subroutine test_assign_and_copy_semantics()
        type(linked_list) :: a, b
        class(*), allocatable :: x
        write(*,'(A)') 'test_assign_and_copy_semantics'
        call a%push_back(1)
        call a%push_back(2)
        call a%push_back(3)
        call b%assign(a)
        call assert_int(a%size(), b%size(), 'assign preserves size')
        call b%at(3, x)
        call assert_int(3, transfer(x,0), 'assign preserves contents deep')
        ! modify b and ensure a unchanged
        call b%push_back(4)
        call assert_int(3, a%size(), 'copy is independent')
        call a%kill(); call b%kill()
    end subroutine test_assign_and_copy_semantics

    !------------------- append -------------------

    subroutine test_append()
        type(linked_list) :: a, b, c
        class(*), allocatable :: x
        write(*,'(A)') 'test_append'
        call a%push_back(1)
        call a%push_back(2)
        call b%push_back(3)
        call b%push_back(4)
        c = a%append(b)
        call assert_int(4, c%size(), 'append size=4')
        call c%at(4,x); call assert_int(4, transfer(x,0), 'append final element')
        call a%kill(); call b%kill(); call c%kill()
    end subroutine test_append

    !------------------- slice -------------------

    subroutine test_slice()
        type(linked_list) :: lst, sl
        class(*), allocatable:: x
        write(*,'(A)') 'test_slice'
        call lst%push_back(10)
        call lst%push_back(20)
        call lst%push_back(30)
        call lst%slice(2,3,sl)
        call assert_int(2, sl%size(), 'slice size=2')
        call sl%front(x); call assert_int(20, transfer(x,0), 'slice first')
        call sl%back(x);  call assert_int(30, transfer(x,0), 'slice last')
        call lst%kill(); call sl%kill()
    end subroutine test_slice

    !------------------- iterator: manual loop -------------------

    subroutine test_iteration()
        type(linked_list) :: lst
        type(list_iterator) :: it
        class(*), allocatable :: x
        integer :: expected(3) = [10,20,30]
        integer :: i
        write(*,'(A)') 'test_iteration'
        call lst%push_back(10)
        call lst%push_back(20)
        call lst%push_back(30)
        it = lst%begin()
        i = 1
        do while (it%has_value())
            call it%get(x)
            call assert_int(expected(i), transfer(x,0), 'iterator value match')
            call it%next()
            i=i+1
        end do
        call lst%kill()
    end subroutine test_iteration

    !------------------- iterator index + advance -------------------

    subroutine test_iterator_index_and_advance()
        type(linked_list) :: lst
        type(list_iterator) :: it
        class(*), allocatable :: x
        write(*,'(A)') 'test_iterator_index_and_advance'
        call lst%push_back(100)
        call lst%push_back(200)
        call lst%push_back(300)
        it=lst%begin()
        call assert_int(1, it%index(lst), 'index(start)=1')
        call it%advance(2)
        call assert_int(3, it%index(lst), 'index(after advance 2)=3')
        call it%get(x)
        call assert_int(300, transfer(x,0), 'correct value after advance')
        call lst%kill()
    end subroutine test_iterator_index_and_advance


    !------------------- end() iterator semantics -------------------

    subroutine test_end_iterator_behavior()
        type(linked_list) :: lst
        type(list_iterator) :: it_end, it_begin
        write(*,'(A)') 'test_end_iterator_behavior'
        call lst%push_back(1)
        call lst%push_back(2)
        it_begin = lst%begin()
        it_end   = lst%end_iter()
        call assert_false(it_end%has_value(), 'end has no value')
        call assert_false(it_begin%equals(it_end), 'begin != end')
        call lst%kill()
    end subroutine test_end_iterator_behavior

    !------------------- kill -------------------

    subroutine test_kill()
        type(linked_list) :: lst
        write(*,'(A)') 'test_kill'
        call lst%push_back(999)
        call lst%kill()
        call assert_true(lst%is_empty(), 'kill leaves list empty')
        call assert_int(0, lst%size(), 'size after kill=0')
    end subroutine test_kill

end module simple_linked_list_tester
