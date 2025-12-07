module simple_poly_rec_list_tester
use simple_poly_rec_list
use simple_linked_list
use simple_string
implicit none

public :: run_poly_rec_list_test
private

contains

    subroutine run_poly_rec_list_test
        type(poly_rec_list)   :: lst, lst2, lst3, lst4
        type(project_rec)     :: p, p2
        type(process_rec)     :: pr, pr2
        type(chunk_rec)       :: c, c2
        type(list_iterator)   :: it
        class(*), allocatable :: tmp
        integer :: ierr, i
        logical :: failed
        failed = .false.
        print *, "========================================="
        print *, " UNIT TEST: poly_rec_list (allocatable)"
        print *, "========================================="
        ! ------------------------
        ! Test PUSH + SIZE
        ! ------------------------
        print *, "Test: push() and size()"
        p%projname ="projA"; p%micind=1; p%nptcls=120
        pr%str_id = "PROC_1"
        c%id = 101
        call lst%push(p)
        call lst%push(pr)
        call lst%push(c)
        if (lst%size() /= 3) then
            print *, "** FAIL: Size after push incorrect"
            failed = .true.
        endif
        ! ------------------------
        ! Test GETTERS (type safety)
        ! ------------------------
        print *, "Test: getters + type enforcement"
        call lst%get(1, p2)
        if (p2%projname /= "projA") then
            print *, "** FAIL: get returned wrong record"
            failed = .true.
        endif
        call lst%get(2, pr2)
        if (pr2%str_id%to_char() /= "PROC_1") then
            print *, "** FAIL: process_rec getter failed"
            failed = .true.
        endif
        call lst%get(3, c2)
        if (c2%id /= 101) then
            print *, "** FAIL: chunk_rec getter failed"
            failed = .true.
        endif
        ! ------------------------
        ! Test SET (replace existing + append)
        ! ------------------------
        print *, "Test: set() behavior"
        p2%projname="projB"
        call lst%set(1, p2)
        call lst%get(1, p)
        if( p%projname%to_char() /= "projB") then
            print *, "** FAIL: set() did not update value correctly"
            failed = .true.
        endif
        ! Append using set(index=n+1)
        pr2%str_id="PROC_2"
        call lst%set(lst%size()+1, pr2)
        if (lst%size() /= 4) then
            print *, "** FAIL: set(n+1) did not append"
            failed = .true.
        endif
        ! ------------------------
        ! Test COPY_FROM
        ! ------------------------
        print *, "Test: copy_from()"
        lst2 = lst        ! initial copy
        call lst2%copy_from(lst)
        if (lst2%size() /= lst%size()) then
            print *, "** FAIL: copy_from size mismatch"
            failed = .true.
        endif
        ! ------------------------
        ! Test concatenation operator //
        ! ------------------------
        print *, "Test: operator //"
        lst3 = lst // lst2
        if (lst3%size() /= lst%size()*2) then
            print *, "** FAIL: // append operator incorrect"
            failed = .true.
        endif
        ! ------------------------
        ! Test assignment (=) move semantics
        ! ------------------------
        print *, "Test: move assignment"
        lst4 = lst3
        if (lst4%size() /= lst3%size()) then
            print *, "** FAIL: assignment operator failed or shallow copy mismatch"
            failed = .true.
        endif
        ! ------------------------
        ! Test kill()
        ! ------------------------
        print *, "Test: kill() clear list"
        call lst4%kill()
        if (.not. lst4%is_empty()) then
            print *, "** FAIL: kill() did not empty list"
            failed = .true.
        endif
        ! ------------------------
        ! Summary
        ! ------------------------
        if (failed) then
            print *, "❌ TEST SUITE FAILED"
            stop 1
        else
            print *, "✅ All tests passed successfully."
        endif
    end subroutine run_poly_rec_list_test

end module simple_poly_rec_list_tester
