!@descr: unit tests of particle_sieving_2's rules: the seed and batch cut with the seed's two limits, the retraining rounds, the junk labels, the rounds' low-pass limit and agreement, the student's table by class, the freeze, the greedy update fraction, the final states and the two stops
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
        call test_round_lowpass()
        call test_round_agreement()
        call test_class_agreement()
        call test_kept_order()
        call test_class_pops()
        call test_freeze_rules()
        call test_final_states()
        call test_both_halves()
        call test_classes_left_stop()
        call test_none_pending_stop()
    end subroutine run_all_ptcl_sieve2_tests

    !> without the particle cap the seed takes the first nmics stacks with selected particles; batches close at
    !! batch_nptcls in project order, the tail is the last batch; a stack without selected particles is in none
    subroutine test_plan_batches()
        class(ptcl_sieve2), allocatable :: sieve
        integer,            allocatable :: stk_batch(:)
        integer,            parameter   :: NSEL(9) = [5, 0, 3, 4, 6, 2, 7, 0, 1]
        integer                         :: nbatches, b
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
        integer,            allocatable :: stk_batch(:)
        integer                         :: nbatches
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
        integer,            allocatable :: stk_batch(:)
        integer,            parameter   :: NSEL(6) = [5, 0, 3, 4, 6, 2]
        integer                         :: nbatches
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
        integer,            allocatable :: rounds(:)
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
        logical,            allocatable :: bad(:)
        real(dp),           parameter   :: PJUNK(6) = [0.95d0, 0.5d0, 0.91d0, 0.9d0, 0.99d0, 0.2d0]
        logical,            parameter   :: EMPTY(6) = [.false., .false., .false., .false., .true., .false.]
        integer,            parameter   :: POPS(6)  = [3, 4, 0, 2, 0, 5]
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

    !> the new particles and old_frac of the earlier ones, at least as many earlier as new, at most every particle;
    !! old_frac 0 is twice the new particles
    subroutine test_greedy_frac()
        class(ptcl_sieve2), allocatable :: sieve
        write(*,'(A)') 'test_greedy_frac'
        allocate(sieve)
        call assert_real(0.325, sieve%greedy_frac(100, 1000, 0.25), 1.e-6, 'the new and a quarter of the 900 earlier')
        call assert_real(0.6,   sieve%greedy_frac(300, 1000, 0.25), 1.e-6, 'never fewer earlier particles than new')
        call assert_real(0.2,   sieve%greedy_frac(100, 1000, 0.),   1.e-6, 'old_frac 0: twice the new particles')
        call assert_real(1.0,   sieve%greedy_frac(100, 1000, 1.),   1.e-6, 'old_frac 1: every particle')
        call assert_real(1.0,   sieve%greedy_frac(600, 1000, 0.25), 1.e-6, 'at most every particle')
        call assert_real(1.0,   sieve%greedy_frac(0, 0, 0.25),      1.e-6, 'no particle: a full update')
        deallocate(sieve)
    end subroutine test_greedy_frac

    !> a round's low-pass limit is the class FRCs' estimate, held at lpstop when the FRCs reach finer
    subroutine test_round_lowpass()
        class(ptcl_sieve2), allocatable :: sieve
        write(*,'(A)') 'test_round_lowpass'
        allocate(sieve)
        call assert_real(10., sieve%round_lowpass(6.4, 10.),  1.e-6, 'FRCs finer than lpstop: held at lpstop')
        call assert_real(14., sieve%round_lowpass(14., 10.),  1.e-6, 'FRCs coarser than lpstop: the FRCs'' limit')
        deallocate(sieve)
    end subroutine test_round_lowpass

    !> the round's agreement: new particles (no update before) that kept their class, also among those whose class was
    !! junk, earlier ones searched again that kept theirs; particles the round did not search and inactive ones count
    !! in neither
    subroutine test_round_agreement()
        class(ptcl_sieve2), allocatable :: sieve
        logical, parameter :: ACTIVE(7) = [.true., .true., .true., .true., .true., .true., .false.]
        integer, parameter :: CLS0(7)   = [3, 4, 5, 6, 7, 8, 9]
        integer, parameter :: CLS1(7)   = [3, 1, 5, 2, 7, 8, 1]
        integer, parameter :: CNT0(7)   = [0, 0, 1, 1, 2, 1, 0]
        integer, parameter :: CNT1(7)   = [1, 1, 2, 2, 2, 1, 1]
        logical, parameter :: BAD(9)    = [.false., .false., .false., .true., .false., .false., .false., .false., .true.]
        integer :: nnew, nnew_kept, njunk, njunk_kept, nold, nold_kept
        write(*,'(A)') 'test_round_agreement'
        allocate(sieve)
        call sieve%round_agreement(ACTIVE, CLS0, CLS1, CNT0, CNT1, BAD, nnew, nnew_kept, njunk, njunk_kept, nold, &
            &nold_kept)
        call assert_int(2, nnew,       'two new particles searched (the inactive one not counted)')
        call assert_int(1, nnew_kept,  'one of them kept the student''s class')
        call assert_int(1, njunk,      'one of them came with a junk class (the inactive one in class 9 not counted)')
        call assert_int(0, njunk_kept, 'and left it')
        call assert_int(2, nold,      'two earlier particles searched again (the unsearched ones not counted)')
        call assert_int(1, nold_kept, 'one of them kept its class')
        deallocate(sieve)
    end subroutine test_round_agreement

    !> the class counts: a new particle searched counts in its predicted class (kept when the round left it there) and
    !! in the class the round put it in; an earlier, an unsearched and an inactive particle in none
    subroutine test_class_agreement()
        class(ptcl_sieve2), allocatable :: sieve
        logical, parameter :: ACTIVE(5) = [.true., .true., .true., .true., .false.]
        integer, parameter :: CLS0(5)   = [2, 2, 3, 1, 2]
        integer, parameter :: CLS1(5)   = [2, 3, 3, 3, 2]
        integer, parameter :: CNT0(5)   = [0, 0, 1, 0, 0]
        integer, parameter :: CNT1(5)   = [1, 1, 2, 0, 1]
        integer :: npred(3), nkept(3), ninto(3)
        write(*,'(A)') 'test_class_agreement'
        allocate(sieve)
        npred = 0
        nkept = 0
        ninto = 0
        call sieve%class_agreement(ACTIVE, CLS0, CLS1, CNT0, CNT1, npred, nkept, ninto)
        call assert_true(all(npred == [0, 2, 0]), 'two new particles predicted into class 2 count')
        call assert_true(all(nkept == [0, 1, 0]), 'one of them kept')
        call assert_true(all(ninto == [0, 1, 1]), 'one placed in class 2, one in class 3')
        deallocate(sieve)
    end subroutine test_class_agreement

    !> the classes by the share kept, best first, ties in class order, those without a prediction last
    subroutine test_kept_order()
        class(ptcl_sieve2), allocatable :: sieve
        integer, allocatable :: order(:)
        write(*,'(A)') 'test_kept_order'
        allocate(sieve)
        order = sieve%kept_order([4, 0, 10, 2, 8], [1, 0, 9, 1, 4])
        call assert_true(all(order == [3, 4, 5, 1, 2]), '90 %, 50 % before 50 % in class order, 25 %, none last')
        deallocate(sieve)
    end subroutine test_kept_order

    !> the populations count the masked rows by class; classes outside 1..ncls count in none
    subroutine test_class_pops()
        class(ptcl_sieve2), allocatable :: sieve
        integer, allocatable :: pops(:)
        write(*,'(A)') 'test_class_pops'
        allocate(sieve)
        call sieve%class_pops([.true., .true., .false., .true., .true.], [2, 2, 1, 0, 3], 3, pops)
        call assert_true(all(pops == [0, 2, 1]), 'two in class 2, one in class 3, the masked and the unclassified not')
        deallocate(sieve)
    end subroutine test_class_pops

    !> the freeze: rounds at the threshold or above count up, one below starts again; a frozen round's fraction is the
    !! new particles over the ingested ones
    subroutine test_freeze_rules()
        class(ptcl_sieve2), allocatable :: sieve
        write(*,'(A)') 'test_freeze_rules'
        allocate(sieve)
        call assert_int(1, sieve%freeze_count(0, 0.995, 0.99), 'a round at the threshold or above counts')
        call assert_int(2, sieve%freeze_count(1, 0.99,  0.99), 'the threshold itself counts')
        call assert_int(0, sieve%freeze_count(2, 0.985, 0.99), 'a round below starts the count again')
        call assert_int(0, sieve%freeze_count(0, 1.,    1.01), 'a threshold above 1 never counts')
        call assert_real(0.25, sieve%newonly_frac(250, 1000), 1.e-6, 'the new particles over the ingested')
        call assert_real(1.,   sieve%newonly_frac(0, 0),      1.e-6, 'no particle: a full update')
        deallocate(sieve)
    end subroutine test_freeze_rules

    !> accepted: selected, ingested and of a class that is not junk; a particle without a class is rejected
    subroutine test_final_states()
        class(ptcl_sieve2), allocatable :: sieve
        integer,            allocatable :: states(:)
        logical,            parameter   :: SEL(6)      = [.true., .true., .true., .false., .true., .true.]
        logical,            parameter   :: INGESTED(6) = [.true., .true., .true., .true.,  .true., .true.]
        integer,            parameter   :: CLASSES(6)  = [1, 2, 3, 1, 0, 4]
        logical,            parameter   :: BAD(3)      = [.false., .true., .false.]
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
