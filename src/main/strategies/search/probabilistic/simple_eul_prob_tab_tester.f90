!@descr: unit tests for the 3D probability table: the shift recorded for an assigned, unrefined candidate (simple_eul_prob_tab_utils)
! The seed table holds the particle-frame seed every candidate of a particle was scored with; an assigned, unrefined
! candidate records it as it is, whatever its in-plane angle. A refined candidate keeps the shift of its search.
module simple_eul_prob_tab_tester
use simple_core_module_api
use simple_eul_prob_tab_utils, only: record_seed_shift
use simple_type_defs,          only: ptcl_ref
use simple_test_utils
implicit none
private
public :: run_all_eul_prob_tab_tests

! eight in-plane rotations
integer, parameter :: NROTS = 8
real,    parameter :: SEED(2) = [1.5, -0.5]

contains

    subroutine run_all_eul_prob_tab_tests()
        write(*,'(A)') '**** running all 3D probability table tests ****'
        call test_record_seed_shift()
    end subroutine run_all_eul_prob_tab_tests

    subroutine test_record_seed_shift()
        type(ptcl_ref) :: assgn
        integer :: inpl
        write(*,'(A)') 'test_record_seed_shift'
        ! the particle-frame seed is recorded whatever the angle
        do inpl = 1, NROTS
            assgn = unrefined(inpl)
            call record_seed_shift(assgn, SEED, .true., .true.)
            call assert_true(assgn%has_sh .and. assgn%x == SEED(1) .and. assgn%y == SEED(2), &
                &'the seed is recorded as it is, at every in-plane angle')
        enddo
        ! a refined candidate keeps its own shift, a particle without seed or without shift search is left alone
        assgn        = unrefined(3)
        assgn%has_sh = .true.
        assgn%x      = 0.25
        assgn%y      = -0.75
        call record_seed_shift(assgn, SEED, .true., .true.)
        call assert_true(assgn%x == 0.25 .and. assgn%y == -0.75, 'a refined candidate keeps the shift of its search')
        assgn = unrefined(3)
        call record_seed_shift(assgn, SEED, .false., .true.)
        call assert_false(assgn%has_sh, 'no seed, nothing recorded')
        call record_seed_shift(assgn, SEED, .true., .false.)
        call assert_false(assgn%has_sh, 'no shift search, nothing recorded')
    end subroutine test_record_seed_shift

    function unrefined( inpl ) result( assgn )
        integer, intent(in) :: inpl
        type(ptcl_ref) :: assgn
        assgn%inpl = inpl
    end function unrefined

end module simple_eul_prob_tab_tester
