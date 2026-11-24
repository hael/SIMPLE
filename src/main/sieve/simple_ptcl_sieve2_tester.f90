!@descr: unit tests of particle_sieving_2's rules: the seed and batch cut with the seed's two limits, the retraining rounds, the junk labels, the greedy update fraction, the final states and the two stops
submodule (simple_ptcl_sieve2) simple_ptcl_sieve2_tester
use simple_test_utils,     only: assert_true, assert_false, assert_int, assert_real
use simple_cmdline,        only: cmdline
use simple_forked_process, only: forked_process, FORK_STATUS_FAILED, FORK_STATUS_STOPPED
implicit none

!> a check run in a child process: with l_violate the input breaks the rule and the child must stop, else it must
!! end normally
type, extends(forked_process) :: classes_left_run
    logical :: l_violate = .false.
contains
    procedure :: execute => run_classes_left
end type classes_left_run

type, extends(forked_process) :: none_pending_run
    logical :: l_violate = .false.
contains
    procedure :: execute => run_none_pending
end type none_pending_run

contains

    module subroutine run_all_ptcl_sieve2_tests()
        write(*,'(A)') '**** running all ptcl_sieve2 tests ****'
        call test_plan_batches()
        call test_plan_without_batches()
        call test_plan_seed_cap()
        call test_retrains_after()
        call test_class_labels()
        call test_greedy_frac()
        call test_final_states()
        call test_both_halves()
        call test_classes_left_stop()
        call test_none_pending_stop()
    end subroutine run_all_ptcl_sieve2_tests

    !> without the particle cap the seed takes the first nmics stacks with selected particles; batches close at
    !! batch_nptcls in project order, the tail is the last batch; a stack without selected particles is in none
    subroutine test_plan_batches()
        class(ptcl_sieve2), allocatable :: sieve
        integer, allocatable :: stk_batch(:)
        integer, parameter   :: NSEL(9) = [5, 0, 3, 4, 6, 2, 7, 0, 1]
        integer :: nbatches, b
        write(*,'(A)') 'test_plan_batches'
        allocate(sieve)
        call sieve%plan_batches(NSEL, 2, huge(1), 6, stk_batch, nbatches)
        call assert_int(3, nbatches, 'three batches, the tail the third')
        call assert_true(all(stk_batch == [0, -1, 0, 1, 1, 2, 2, -1, 3]), 'the seed, the batches and the empty stacks')
        do b = 1,nbatches - 1
            call assert_true(sum(NSEL, mask=stk_batch == b) >= 6, 'a batch but the tail holds at least the batch size')
        enddo
        call assert_int(1, sum(NSEL, mask=stk_batch == nbatches), 'the tail holds what is left')
        call assert_int(sum(NSEL), sum(NSEL, mask=stk_batch >= 0), 'every selected particle in the seed or a batch')
        deallocate(sieve)
    end subroutine test_plan_batches

    !> no more stacks than the seed takes: no batch
    subroutine test_plan_without_batches()
        class(ptcl_sieve2), allocatable :: sieve
        integer, allocatable :: stk_batch(:)
        integer :: nbatches
        write(*,'(A)') 'test_plan_without_batches'
        allocate(sieve)
        call sieve%plan_batches([4, 0, 2], 2, huge(1), 10, stk_batch, nbatches)
        call assert_int(0, nbatches, 'no batch')
        call assert_true(all(stk_batch == [0, -1, 0]), 'the seed and the empty stack')
        call sieve%plan_batches([4, 2, 9, 1], 2, huge(1), 3, stk_batch, nbatches)
        call assert_int(2, nbatches, 'a batch that reaches the batch size closes')
        call assert_true(all(stk_batch == [0, 0, 1, 2]), 'and the next stack starts another')
        deallocate(sieve)
    end subroutine test_plan_without_batches

    !> the seed stops at whichever limit comes first: once it holds seed_nptcls particles (the stack that reaches the
    !! cap is its last), or at nmics stacks; what follows goes to the batches
    subroutine test_plan_seed_cap()
        class(ptcl_sieve2), allocatable :: sieve
        integer, allocatable :: stk_batch(:)
        integer, parameter   :: NSEL(6) = [5, 0, 3, 4, 6, 2]
        integer :: nbatches
        write(*,'(A)') 'test_plan_seed_cap'
        allocate(sieve)
        ! the cap first: 5 + 3 reaches 7, so the third stack is the seed's last although nmics allows ten
        call sieve%plan_batches(NSEL, 10, 7, 6, stk_batch, nbatches)
        call assert_true(all(stk_batch == [0, -1, 0, 1, 1, 2]), 'the cap ends the seed with the stack that reaches it')
        call assert_int(8, sum(NSEL, mask=stk_batch == 0), 'the seed holds the cap and the rest of its last stack')
        call assert_int(2, nbatches, 'the stacks after the seed in batches')
        ! a cap met exactly: the next stack is not the seed's
        call sieve%plan_batches(NSEL, 10, 8, 6, stk_batch, nbatches)
        call assert_true(all(stk_batch == [0, -1, 0, 1, 1, 2]), 'a seed at the cap takes no further stack')
        ! nmics first: two stacks hold 8 of a cap of 100
        call sieve%plan_batches(NSEL, 2, 100, 6, stk_batch, nbatches)
        call assert_true(all(stk_batch == [0, -1, 0, 1, 1, 2]), 'nmics ends the seed before the cap')
        ! neither binds: five stacks of 20 particles in all, no batch
        call sieve%plan_batches(NSEL, 10, 100, 6, stk_batch, nbatches)
        call assert_int(0, nbatches, 'every stack in the seed when neither limit is reached')
        call assert_true(all(stk_batch == [0, -1, 0, 0, 0, 0]), 'the seed takes every stack with selected particles')
        deallocate(sieve)
    end subroutine test_plan_seed_cap

    !> the student is trained again after exactly the listed rounds; without a list after none
    subroutine test_retrains_after()
        class(ptcl_sieve2), allocatable :: sieve
        integer, allocatable :: rounds(:)
        write(*,'(A)') 'test_retrains_after'
        allocate(sieve)
        call assert_false(sieve%retrains_after(rounds, 1), 'no list: no retraining')
        rounds = [5, 2]
        call assert_true(sieve%retrains_after(rounds, 2) .and. sieve%retrains_after(rounds, 5), 'after each listed round')
        call assert_false(sieve%retrains_after(rounds, 1) .or. sieve%retrains_after(rounds, 3), 'after no other round')
        deallocate(sieve)
    end subroutine test_retrains_after

    !> junk: populated, not empty and p_junk above the threshold (strictly); a threshold above 1 marks none, one below
    !! 0 every populated class that is not empty
    subroutine test_class_labels()
        class(ptcl_sieve2), allocatable :: sieve
        logical, allocatable :: bad(:)
        real(dp), parameter :: PJUNK(6) = [0.95d0, 0.5d0, 0.91d0, 0.9d0, 0.99d0, 0.2d0]
        logical,  parameter :: EMPTY(6) = [.false., .false., .false., .false., .true., .false.]
        integer,  parameter :: POPS(6)  = [3, 4, 0, 2, 0, 5]
        write(*,'(A)') 'test_class_labels'
        allocate(sieve)
        call sieve%class_labels(PJUNK, EMPTY, POPS, 0.9d0, bad)
        call assert_true(all(bad .eqv. [.true., .false., .false., .false., .false., .false.]), &
            &'junk above 0.9 only when populated and not empty; 0.9 itself is not junk')
        call sieve%class_labels(PJUNK, EMPTY, POPS, 1.01d0, bad)
        call assert_false(any(bad), 'a threshold above 1 marks no class')
        call sieve%class_labels(PJUNK, EMPTY, POPS, -0.01d0, bad)
        call assert_true(all(bad .eqv. (POPS > 0 .and. .not. EMPTY)), 'a threshold below 0 marks every populated class')
        deallocate(sieve)
    end subroutine test_class_labels

    !> the new particles and as many earlier ones, at most every particle
    subroutine test_greedy_frac()
        class(ptcl_sieve2), allocatable :: sieve
        write(*,'(A)') 'test_greedy_frac'
        allocate(sieve)
        call assert_real(0.2,  sieve%greedy_frac(100, 1000), 1.e-6, 'twice the new particles over the ingested')
        call assert_real(1.0,  sieve%greedy_frac(600, 1000), 1.e-6, 'at most every particle')
        call assert_real(1.0,  sieve%greedy_frac(500, 1000), 1.e-6, 'half new: every particle')
        call assert_real(1.0,  sieve%greedy_frac(0, 0),      1.e-6, 'no particle: a full update')
        deallocate(sieve)
    end subroutine test_greedy_frac

    !> accepted: selected, ingested and of a class that is not junk; a particle without a class is rejected
    subroutine test_final_states()
        class(ptcl_sieve2), allocatable :: sieve
        integer, allocatable :: states(:)
        logical, parameter :: SEL(6)      = [.true., .true., .true., .false., .true., .true.]
        logical, parameter :: INGESTED(6) = [.true., .true., .true., .true.,  .true., .true.]
        integer, parameter :: CLASSES(6)  = [1, 2, 3, 1, 0, 4]
        logical, parameter :: BAD(3)      = [.false., .true., .false.]
        write(*,'(A)') 'test_final_states'
        allocate(sieve)
        call sieve%final_states(SEL, INGESTED, CLASSES, BAD, states)
        call assert_true(all(states == [1, 0, 1, 0, 0, 0]), &
            &'accepted in a good class; rejected in a junk class, unselected, without a class or beyond the classes')
        deallocate(sieve)
    end subroutine test_final_states

    !> the pending set is seeded only when it holds both halves (the extension takes each half's mean over it)
    subroutine test_both_halves()
        class(ptcl_sieve2), allocatable :: sieve
        write(*,'(A)') 'test_both_halves'
        allocate(sieve)
        call assert_true(sieve%both_halves([.true., .true., .false.], [0, 1, 1]), 'an even and an odd row')
        call assert_false(sieve%both_halves([.true., .true., .false.], [0, 0, 1]), 'the odd row outside the set')
        call assert_false(sieve%both_halves([.false., .false.], [0, 1]), 'an empty set')
        call assert_false(sieve%both_halves([.true., .true.], [-1, 1]), 'a row in no half counts in neither')
        deallocate(sieve)
    end subroutine test_both_halves

    !> every populated class junk stops the run; one good populated class does not
    subroutine test_classes_left_stop()
        type(classes_left_run) :: child
        write(*,'(A)') 'test_classes_left_stop'
        child%l_violate = .true.
        call child%start(name=string('SIEVE2_NO_CLASS_LEFT'))
        call child%await_final_status()
        call assert_int(FORK_STATUS_FAILED, child%status(), 'no class left stops')
        call child%destroy
        child%l_violate = .false.
        call child%start(name=string('SIEVE2_A_CLASS_LEFT'))
        call child%await_final_status()
        call assert_int(FORK_STATUS_STOPPED, child%status(), 'a good populated class does not')
        call child%destroy
    end subroutine test_classes_left_stop

    !> a selected particle never ingested stops the result; none pending does not
    subroutine test_none_pending_stop()
        type(none_pending_run) :: child
        write(*,'(A)') 'test_none_pending_stop'
        child%l_violate = .true.
        call child%start(name=string('SIEVE2_PENDING_ROW'))
        call child%await_final_status()
        call assert_int(FORK_STATUS_FAILED, child%status(), 'a pending row stops')
        call child%destroy
        child%l_violate = .false.
        call child%start(name=string('SIEVE2_NONE_PENDING'))
        call child%await_final_status()
        call assert_int(FORK_STATUS_STOPPED, child%status(), 'none pending does not')
        call child%destroy
    end subroutine test_none_pending_stop

    subroutine run_classes_left( self, cline )
        class(classes_left_run), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        class(ptcl_sieve2), allocatable :: sieve
        allocate(sieve)
        if( self%l_violate )then
            ! the empty class is no class left
            call sieve%check_classes_left([.true., .true., .false.], [3, 1, 0])
        else
            call sieve%check_classes_left([.true., .false., .false.], [3, 1, 0])
        endif
        deallocate(sieve)
    end subroutine run_classes_left

    subroutine run_none_pending( self, cline )
        class(none_pending_run), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        class(ptcl_sieve2), allocatable :: sieve
        allocate(sieve)
        if( self%l_violate )then
            call sieve%check_none_pending([.true., .true., .false.], [.true., .false., .false.])
        else
            ! an unselected particle need not be ingested
            call sieve%check_none_pending([.true., .true., .false.], [.true., .true., .false.])
        endif
        deallocate(sieve)
    end subroutine run_none_pending

end submodule simple_ptcl_sieve2_tester
