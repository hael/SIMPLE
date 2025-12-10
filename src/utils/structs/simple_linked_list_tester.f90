module simple_linked_list_tester
use simple_linked_list
use simple_test_utils
implicit none

public :: run_all_list_tests
private

!-----------------------------------
! Derived type used in testing
!-----------------------------------
type :: particle
    integer :: id
    real    :: score
    character(:), allocatable :: label
end type particle

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
        ! call report_summary()
        write(*,'(A)') '**** running parametric linked_list tests ****'
        call test_intrinsic_integer()
        call test_real()
        call test_complex()
        call test_logical()
        call test_string_values()
        call test_derived_type_values()
        call test_derived_type_independence()
        call test_nested_allocations()
        call test_move_semantics()
        call test_type_safe_wrappers()
        call test_replace_iterator_stability()
        ! call report_summary()
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
        b = a
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
        c = a//b
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
            call it%getter(x)
            call assert_int(expected(i), transfer(x,0), 'iterator value match')
            ! Replace the middle element during iteration
            if (i == 2) call lst%replace_iterator(it, 777)
            call it%next()
            i=i+1
        end do
        call lst%at(2,x)
        call assert_int(777, transfer(x,0), 'replace_iterator modified correct node')
        deallocate(x)
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
        call it%getter(x)
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

    !===========================
    !  INTRINSIC TYPE TESTS
    !===========================

    subroutine test_intrinsic_integer()
        type(linked_list) :: lst
        class(*), allocatable :: x
        integer :: i
        write(*,'(A)') 'test_intrinsic_integer'
        do i=1,5
            call lst%push_back(i)
        end do
        do i=1,5
            call lst%at(i, x)
            if (.not. allocated(x)) stop "tmp not allocated"
            select type(t => x)
            type is (integer)
                call assert_int(i, t, 'integer match')
            class default
                stop "unexpected type"
            end select
            deallocate(x)
        end do
        call lst%kill()
    end subroutine test_intrinsic_integer

    subroutine test_real()
        type(linked_list) :: lst
        class(*), allocatable :: x
        write(*,'(A)') 'test_real'
        call lst%push_back(1.5)
        call lst%push_back(-2.25)
        call lst%at(1,x)
        select type(t=>x)
        type is (real)
            call assert_real(1.5, t, 1.0e-4, 'real value 1')
        class default
            stop "unexpected type"
        end select
        deallocate(x)
        call lst%at(2,x)
        select type(t=>x)
        type is (real)
            call assert_real(-2.25, t, 1.0e-4, 'real value 2')
        class default
            stop "unexpected type"
        end select
        deallocate(x)
        call lst%kill()
    end subroutine test_real

    subroutine test_complex()
        type(linked_list) :: lst
        class(*), allocatable :: x
        complex :: expected
        write(*,'(A)') 'test_complex'
        expected = (3.0, -4.0)
        call lst%push_back(expected)
        call lst%front(x)
        select type(t=>x)
        type is (complex)
            if (abs(t - expected) > 1.0e-12) stop "complex mismatch"
        class default
            stop "unexpected type"
        end select
        deallocate(x)
        call lst%kill()
    end subroutine test_complex

    subroutine test_logical()
        type(linked_list) :: lst
        class(*), allocatable :: x
        write(*,'(A)') 'test_logical'
        call lst%push_back(.true.)
        call lst%push_back(.false.)
        call lst%at(1,x)
        select type(t=>x)
        type is (logical)
            call assert_true(t, 'logical truth')
        class default
            stop "unexpected type"
        end select
        deallocate(x)
        call lst%at(2,x)
        select type(t=>x)
        type is (logical)
            call assert_false(t, 'logical false')
        class default
            stop "unexpected type"
        end select
        deallocate(x)
        call lst%kill()
    end subroutine test_logical

    !===========================
    ! STRING TESTS
    !===========================

    subroutine test_string_values()
        type(linked_list) :: lst
        class(*), allocatable :: x
        character(:), allocatable :: s
        write(*,'(A)') 'test_string_values'
        call lst%push_back('Alpha')
        call lst%push_back('Beta')
        call lst%push_back('Gamma')
        call lst%at(2,x)
        select type(t=>x)
        type is (character(*))
            s = t
        class default
            stop "unexpected type"
        end select
        call assert_char(s, 'Beta', 'string match')
        deallocate(x); deallocate(s)
        call lst%kill()
    end subroutine test_string_values

    !===========================
    ! DERIVED TYPE TESTS
    !===========================

    subroutine test_derived_type_values()
        type(linked_list) :: lst
        type(particle) :: p
        class(*), allocatable :: x
        type(particle) :: p2
        write(*,'(A)') 'test_derived_type_values'
        p%id=42; p%score=0.987; p%label='P42'
        call lst%push_back(p)
        call lst%front(x)
        select type(t => x)
        type is (particle)
            p2 = t
        class default
            stop "unexpected type"
        end select
        call assert_int(42, p2%id, 'particle id match')
        call assert_real(0.987, p2%score, 1.0e-4, 'particle score match')
        call assert_char(p2%label, 'P42', 'particle label match')
        deallocate(x)
        call lst%kill()
    end subroutine test_derived_type_values

    subroutine test_derived_type_independence()
        type(linked_list) :: lst
        type(particle) :: p
        class(*), allocatable :: x
        type(particle) :: copied
        write(*,'(A)') 'test_derived_type_independence'
        p%id=1; p%score=10; p%label='A'
        call lst%push_back(p)
        ! Modify original
        p%id=999
        p%label='CHANGED'
        ! Retrieve stored version
        call lst%front(x)
        select type(t=>x)
        type is (particle)
            copied = t
        class default
            stop "unexpected type"
        end select
        call assert_int(1, copied%id, 'stored object unchanged by caller mutation')
        deallocate(x)

        call lst%kill()
    end subroutine test_derived_type_independence

    !---------------------------------------
    ! Test deep copy of derived with nested allocatables
    !---------------------------------------

    subroutine test_nested_allocations()
        type(linked_list) :: lst
        type(particle) :: p
        class(*), allocatable :: x
        type(particle) :: restored
        write(*,'(A)') 'test_nested_allocations'
        p%id=77
        p%score=3.14
        p%label = 'NestedTest'
        call lst%push_back(p)
        call lst%at(1,x)
        select type(t=>x)
        type is (particle)
            restored = t
        class default
            stop "unexpected type"
        end select
        ! Modify original external variable
        p%label='corrupted'
        call assert_char(restored%label, 'NestedTest', 'copy preserved allocatable component')
        deallocate(x)
        call lst%kill()
    end subroutine test_nested_allocations

    !---------------------------------------
    ! Test move semantics
    !---------------------------------------

    subroutine test_move_semantics()
        type(linked_list) :: a, b
        class(*), allocatable :: x
        write(*,'(A)') 'test_move_semantics'
        call a%push_back(1)
        call a%push_back(2)
        ! move into b
        call b%replace_with(a)
        call assert_true(a%is_empty(), 'source cleared after move')
        call assert_int(0, a%size(), 'source empty after move')
        call assert_int(2, b%size(), 'dest has elements after move')
        call b%front(x)
        select type(t=>x)
        type is (integer)
            call assert_int(1, t, 'first element preserved after move')
        end select
        deallocate(x)
        call a%kill(); call b%kill()
    end subroutine test_move_semantics

    !---------------------------------------
    ! Test type-safe wrappers
    !---------------------------------------

    subroutine test_type_safe_wrappers()
        type(linked_list) :: lst
        character(:), allocatable :: s
        write(*,'(A)') 'test_type_safe_wrappers'
        call lst%push_back_int(101)
        call lst%push_back_char('hello')
        call lst%at_char(2, s)
        call assert_char(s, 'hello', 'at_char works')
        deallocate(s)
        call lst%kill()
    end subroutine test_type_safe_wrappers

    subroutine test_replace_iterator_stability()
        type(linked_list) :: lst
        type(list_iterator) :: it
        class(*), allocatable :: x
        write(*,'(A)') 'test_replace_iterator_stability'
        call lst%push_back(10)
        call lst%push_back(20)
        call lst%push_back(30)
        it = lst%begin()
        call it%next()  ! now points to second element (20)
        call lst%replace_iterator(it, 999)
        ! Iterator should still reference same logical position
        call it%getter(x)
        call assert_int(999, transfer(x,0), 'iterator kept position after replace')
        ! List structure intact
        call lst%at(1,x); call assert_int(10,  transfer(x,0), 'first ok')
        call lst%at(2,x); call assert_int(999, transfer(x,0), 'replacement ok')
        call lst%at(3,x); call assert_int(30,  transfer(x,0), 'third ok')
        deallocate(x)
        call lst%kill()
    end subroutine test_replace_iterator_stability

end module simple_linked_list_tester
