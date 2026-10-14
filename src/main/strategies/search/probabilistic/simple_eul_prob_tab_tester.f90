!@descr: unit tests for the 3D probability table: the shift recorded for an assigned, unrefined candidate (simple_eul_prob_tab)
! The seed table holds either a reference-frame seed, rotated by the candidate's in-plane angle when it is
! assigned, or with prob_fast (single state) the particle-frame seed every candidate was scored with, recorded as it is.
! A refined candidate keeps the shift of its search under both rules.
module simple_eul_prob_tab_tester
use simple_core_module_api
use simple_eul_prob_tab, only: record_seed_shift
use simple_type_defs,    only: ptcl_ref
use simple_test_utils
implicit none
private
public :: run_all_eul_prob_tab_tests

! eight rotations: index 3 is a quarter turn
integer, parameter :: NROTS = 8
real,    parameter :: SEED(2) = [1.5, -0.5]

contains

    subroutine run_all_eul_prob_tab_tests()
        write(*,'(A)') '**** running all 3D probability table tests ****'
        call test_record_seed_shift()
    end subroutine run_all_eul_prob_tab_tests

    subroutine test_record_seed_shift()
        type(ptcl_ref) :: assgn
        integer :: iseed
        logical :: l_seed_ptcl
        write(*,'(A)') 'test_record_seed_shift'
        ! default: the reference-frame seed rotated by the candidate's angle, 90 degrees at index 3
        assgn = unrefined(3)
        call record_seed_shift(assgn, SEED, .true., .true., NROTS, .false.)
        call assert_true(assgn%has_sh, 'default: the seed is recorded')
        call assert_true(abs(assgn%x - 0.5) <= 1.e-5 .and. abs(assgn%y - 1.5) <= 1.e-5, &
            &'default: the seed is rotated by the in-plane angle of the candidate')
        ! prob_fast (single state): the particle-frame seed is recorded whatever the angle
        do iseed = 1, NROTS
            assgn = unrefined(iseed)
            call record_seed_shift(assgn, SEED, .true., .true., NROTS, .true.)
            call assert_true(assgn%has_sh .and. assgn%x == SEED(1) .and. assgn%y == SEED(2), &
                &'prob_fast: the seed is recorded as it is, at every in-plane angle')
        enddo
        ! a refined candidate keeps its own shift, a particle without seed or without shift search is left alone
        do iseed = 0, 1
            l_seed_ptcl = iseed == 1
            assgn        = unrefined(3)
            assgn%has_sh = .true.
            assgn%x      = 0.25
            assgn%y      = -0.75
            call record_seed_shift(assgn, SEED, .true., .true., NROTS, l_seed_ptcl)
            call assert_true(assgn%x == 0.25 .and. assgn%y == -0.75, 'a refined candidate keeps the shift of its search')
            assgn = unrefined(3)
            call record_seed_shift(assgn, SEED, .false., .true., NROTS, l_seed_ptcl)
            call assert_false(assgn%has_sh, 'no seed, nothing recorded')
            call record_seed_shift(assgn, SEED, .true., .false., NROTS, l_seed_ptcl)
            call assert_false(assgn%has_sh, 'no shift search, nothing recorded')
        enddo
    end subroutine test_record_seed_shift

    function unrefined( inpl ) result( assgn )
        integer, intent(in) :: inpl
        type(ptcl_ref) :: assgn
        assgn%inpl = inpl
    end function unrefined

end module simple_eul_prob_tab_tester
