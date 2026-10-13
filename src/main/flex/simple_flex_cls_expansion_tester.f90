!@descr: unit tests for the cls_expansion subclass budget (simple_flex_cls_expansion)
! Expected allocations come from a Python emulation of the same highest-averages rule (sqrt group
! population between groups, population within a group, pop/minpop ceiling, split threshold).
module simple_flex_cls_expansion_tester
use simple_flex_cls_expansion, only: flex_cls_allocate_subclasses, flex_cls_split_floor
use simple_test_utils
implicit none
private
public :: run_all_flex_cls_expansion_tests

contains

    subroutine run_all_flex_cls_expansion_tests()
        write(*,'(A)') '**** running all flex cls_expansion tests ****'
        call test_split_floor()
        call test_allocation_views()
        call test_allocation_equal_views()
        call test_allocation_capacity()
        call test_allocation_one_per_class()
        call test_allocation_group_labels()
    end subroutine run_all_flex_cls_expansion_tests

    subroutine test_split_floor()
        write(*,'(A)') 'test_split_floor'
        call assert_int(27, flex_cls_split_floor(2, 8),  'model floor dominates at neigs=8')
        call assert_int(40, flex_cls_split_floor(20, 8), 'two members per subclass dominates')
    end subroutine test_split_floor

    !> a preferred view over six classes, a view in two classes, a single class and a small view
    subroutine test_allocation_views()
        integer, parameter :: POPS(12)   = [3200, 3100, 3000, 2900, 2900, 2800, 2500, 1800, 600, 400, 250, 120]
        integer, parameter :: GROUPS(12) = [1, 1, 1, 1, 1, 1, 2, 2, 3, 4, 4, 4]
        integer, parameter :: NSUB_REF(12) = [4, 3, 3, 3, 3, 3, 5, 4, 4, 2, 1, 1] ! alloc_ref.py
        integer :: nsub(12), ngot, c
        write(*,'(A)') 'test_allocation_views'
        call flex_cls_allocate_subclasses(POPS, GROUPS, 36, 50, 100, nsub, ngot)
        call assert_int(36, ngot,      'the total is met')
        call assert_int(36, sum(nsub), 'the counts add up to the total')
        do c = 1, 12
            call assert_int(NSUB_REF(c), nsub(c), 'subclasses of a class match the reference')
        end do
        call assert_int(1, nsub(12), 'the smallest class of its group gets no extra subclass')
    end subroutine test_allocation_views

    !> one class per view, equal populations: every class gets the same count
    subroutine test_allocation_equal_views()
        integer :: pops(60), groups(60), nsub(60), ngot, c
        write(*,'(A)') 'test_allocation_equal_views'
        pops   = 60
        groups = [(c, c=1,60)]
        call flex_cls_allocate_subclasses(pops, groups, 180, 20, 40, nsub, ngot)
        call assert_int(180, ngot, 'the total is met')
        call assert_true(all(nsub == 3), 'three subclasses per view')
    end subroutine test_allocation_equal_views

    !> a total above the capacities delivers the capacities
    subroutine test_allocation_capacity()
        integer :: nsub(3), ngot
        write(*,'(A)') 'test_allocation_capacity'
        call flex_cls_allocate_subclasses([130, 30, 260], [7, 7, 9], 100, 50, 100, nsub, ngot)
        call assert_int(8, ngot,     'the capacities are delivered')
        call assert_int(2, nsub(1),  'pop/minpop caps a class')
        call assert_int(1, nsub(2),  'a small class is kept whole')
        call assert_int(5, nsub(3),  'pop/minpop caps a class')
    end subroutine test_allocation_capacity

    !> a total equal to the class count splits nothing
    subroutine test_allocation_one_per_class()
        integer :: nsub(4), ngot
        write(*,'(A)') 'test_allocation_one_per_class'
        call flex_cls_allocate_subclasses([900, 800, 700, 600], [1, 1, 2, 3], 4, 25, 50, nsub, ngot)
        call assert_int(4, ngot, 'one per class')
        call assert_true(all(nsub == 1), 'every class kept whole')
    end subroutine test_allocation_one_per_class

    !> group labels need be neither consecutive nor ordered
    subroutine test_allocation_group_labels()
        integer :: nsub(3), ngot
        write(*,'(A)') 'test_allocation_group_labels'
        call flex_cls_allocate_subclasses([500, 500, 1000], [5, 2, 5], 8, 100, 200, nsub, ngot)
        call assert_int(8, ngot,    'the total is met')
        call assert_int(2, nsub(1), 'subclasses of a class match the reference')
        call assert_int(3, nsub(2), 'subclasses of a class match the reference')
        call assert_int(3, nsub(3), 'subclasses of a class match the reference')
    end subroutine test_allocation_group_labels

end module simple_flex_cls_expansion_tester
