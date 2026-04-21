!@descr: unit tests for the cluster2D_rejector type (simple_cluster2D_rejector module)
!
! Testable surface
! ----------------
! reject_pop and reject_res operate only on a cls oris object and the
! per-class rejection mask — no image pixel data is read. They are exercised
! with fully synthetic oris whose pop/res fields are set to known values, so
! outcomes are deterministic.
!
! reject_mask operates on the stored image array. An all-zero image produces
! an all-zero band-pass output, which Otsu-thresholds to all-zero, yielding
! no connected components → every class is rejected. This is the one
! image-based test included here.
!
! reject_local_variance is exercised with a deterministic assertion: all-zero
! images produce zero inside/outside variance scores, which are unconditionally
! rejected before z-scoring. All classes in an all-zero set are therefore
! guaranteed rejected regardless of the z-score distribution.
!
! Constants assumed from the production module (not re-exported, so tested
! implicitly through observable outcomes):
!   POP_PERCENT_THRESHOLD = 0.005   threshold = ceiling(sum(pop) * 0.005)
!   RES_THRESHOLD         = 40.0 Å  reject if res > 40.0

module simple_cluster2D_rejector_tester
  use simple_cluster2D_rejector, only: cluster2D_rejector
  use simple_test_utils
  use simple_image,              only: image
  use simple_oris,               only: oris
  implicit none

  public :: run_all_cluster2D_rejector_tests
  private

contains

  subroutine run_all_cluster2D_rejector_tests()
    write(*,'(A)') '**** running cluster2D_rejector tests ****'
    call test_kill_on_default_object()
    call test_kill_idempotent()
    call test_new_initial_state()
    call test_get_rejected_correct_size()
    call test_get_states_correct_size()
    call test_reject_pop_flags_underpopulated()
    call test_reject_pop_keeps_above_threshold()
    call test_reject_res_flags_poor_resolution()
    call test_reject_res_keeps_good_resolution()
    call test_reject_res_boundary_at_threshold()
    call test_rejection_accumulates_across_criteria()
    call test_get_states_mirrors_get_rejected()
    call test_reject_mask_all_zero_images_all_rejected()
    call test_reject_local_variance_all_zero_all_rejected()
  end subroutine run_all_cluster2D_rejector_tests

  ! ---------------------------------------------------------------------------
  ! Helpers
  ! ---------------------------------------------------------------------------

  ! Allocates n zero-filled images at box x box x 1 pixels with the given smpd.
  subroutine make_imgs( imgs, n, box, smpd )
    type(image), allocatable, intent(out) :: imgs(:)
    integer,                  intent(in)  :: n, box
    real,                     intent(in)  :: smpd
    integer :: i
    allocate(imgs(n))
    do i = 1, n
      call imgs(i)%new([box, box, 1], smpd)
    end do
  end subroutine make_imgs

  subroutine kill_imgs( imgs )
    type(image), allocatable, intent(inout) :: imgs(:)
    integer :: i
    if( .not. allocated(imgs) ) return
    do i = 1, size(imgs)
      call imgs(i)%kill()
    end do
    deallocate(imgs)
  end subroutine kill_imgs

  ! ---------------------------------------------------------------------------
  ! kill() safety
  ! ---------------------------------------------------------------------------

  subroutine test_kill_on_default_object()
    type(cluster2D_rejector) :: rjct
    write(*,'(A)') 'test_kill_on_default_object'
    call rjct%kill()
    call assert_true(.true., 'kill on default-declared object completes without error')
  end subroutine test_kill_on_default_object

  subroutine test_kill_idempotent()
    type(cluster2D_rejector) :: rjct
    type(image), allocatable :: imgs(:)
    write(*,'(A)') 'test_kill_idempotent'
    call make_imgs(imgs, 3, 32, 2.0)
    call rjct%new(imgs, 50.0)
    call rjct%kill()
    call rjct%kill()
    call assert_true(.true., 'double kill completes without error')
    call kill_imgs(imgs)
  end subroutine test_kill_idempotent

  ! ---------------------------------------------------------------------------
  ! Initial state after new()
  ! ---------------------------------------------------------------------------

  subroutine test_new_initial_state()
    type(cluster2D_rejector) :: rjct
    type(image), allocatable :: imgs(:)
    logical, allocatable     :: rejected(:)
    integer, allocatable     :: states(:)
    integer                  :: i
    write(*,'(A)') 'test_new_initial_state'
    call make_imgs(imgs, 4, 32, 2.0)
    call rjct%new(imgs, 50.0)
    rejected = rjct%get_rejected()
    states   = rjct%get_states()
    call assert_int(4, size(rejected), 'get_rejected size equals class count after new')
    call assert_int(4, size(states),   'get_states size equals class count after new')
    do i = 1, 4
      call assert_false(rejected(i), 'no class pre-rejected immediately after new')
      call assert_int(1, states(i),  'all states = 1 (kept) immediately after new')
    end do
    call rjct%kill()
    call kill_imgs(imgs)
  end subroutine test_new_initial_state

  subroutine test_get_rejected_correct_size()
    type(cluster2D_rejector) :: rjct
    type(image), allocatable :: imgs(:)
    logical, allocatable     :: rejected(:)
    write(*,'(A)') 'test_get_rejected_correct_size'
    call make_imgs(imgs, 6, 32, 2.0)
    call rjct%new(imgs, 50.0)
    rejected = rjct%get_rejected()
    call assert_int(6, size(rejected), 'get_rejected returns one entry per class')
    call rjct%kill()
    call kill_imgs(imgs)
  end subroutine test_get_rejected_correct_size

  subroutine test_get_states_correct_size()
    type(cluster2D_rejector) :: rjct
    type(image), allocatable :: imgs(:)
    integer, allocatable     :: states(:)
    write(*,'(A)') 'test_get_states_correct_size'
    call make_imgs(imgs, 6, 32, 2.0)
    call rjct%new(imgs, 50.0)
    states = rjct%get_states()
    call assert_int(6, size(states), 'get_states returns one entry per class')
    call rjct%kill()
    call kill_imgs(imgs)
  end subroutine test_get_states_correct_size

  ! ---------------------------------------------------------------------------
  ! reject_pop
  !
  ! threshold = ceiling(sum(pop) * POP_PERCENT_THRESHOLD)
  !           = ceiling(sum(pop) * 0.005)
  !
  ! Case A — 4 classes, pops = [0, 50, 50, 50]
  !   sum = 150, threshold = ceiling(0.75) = 1
  !   class 1 (pop = 0): 0 < 1  → rejected
  !   classes 2-4 (pop = 50): 50 >= 1 → kept
  ! ---------------------------------------------------------------------------

  subroutine test_reject_pop_flags_underpopulated()
    type(cluster2D_rejector) :: rjct
    type(image), allocatable :: imgs(:)
    type(oris)               :: cls_os
    logical, allocatable     :: rejected(:)
    write(*,'(A)') 'test_reject_pop_flags_underpopulated'
    call make_imgs(imgs, 4, 32, 2.0)
    call rjct%new(imgs, 50.0)
    call cls_os%new(4, is_ptcl=.false.)
    call cls_os%set(1, 'pop',  0.0)
    call cls_os%set(2, 'pop', 50.0)
    call cls_os%set(3, 'pop', 50.0)
    call cls_os%set(4, 'pop', 50.0)
    call rjct%reject_pop(cls_os)
    rejected = rjct%get_rejected()
    call assert_true(rejected(1),  'reject_pop flags class with pop = 0')
    call assert_false(rejected(2), 'reject_pop keeps class with pop = 50 (2)')
    call assert_false(rejected(3), 'reject_pop keeps class with pop = 50 (3)')
    call assert_false(rejected(4), 'reject_pop keeps class with pop = 50 (4)')
    call cls_os%kill()
    call rjct%kill()
    call kill_imgs(imgs)
  end subroutine test_reject_pop_flags_underpopulated

  ! Case B — 3 classes, pops = [20, 30, 50]
  !   sum = 100, threshold = ceiling(0.5) = 1
  !   all pops >= 1 → none rejected
  subroutine test_reject_pop_keeps_above_threshold()
    type(cluster2D_rejector) :: rjct
    type(image), allocatable :: imgs(:)
    type(oris)               :: cls_os
    logical, allocatable     :: rejected(:)
    integer                  :: i
    write(*,'(A)') 'test_reject_pop_keeps_above_threshold'
    call make_imgs(imgs, 3, 32, 2.0)
    call rjct%new(imgs, 50.0)
    call cls_os%new(3, is_ptcl=.false.)
    call cls_os%set(1, 'pop', 20.0)
    call cls_os%set(2, 'pop', 30.0)
    call cls_os%set(3, 'pop', 50.0)
    call rjct%reject_pop(cls_os)
    rejected = rjct%get_rejected()
    do i = 1, 3
      call assert_false(rejected(i), 'reject_pop keeps all classes when all above threshold')
    end do
    call cls_os%kill()
    call rjct%kill()
    call kill_imgs(imgs)
  end subroutine test_reject_pop_keeps_above_threshold

  ! ---------------------------------------------------------------------------
  ! reject_res  (RES_THRESHOLD = 40.0 Å; reject if res > 40.0)
  !
  ! Case A — 4 classes, res = [15.0, 30.0, 40.0, 55.0]
  !   class 4 (res = 55): 55 > 40 → rejected
  !   classes 1-3: <= 40 → kept (res = 40.0 is at the boundary, not > 40)
  ! ---------------------------------------------------------------------------

  subroutine test_reject_res_flags_poor_resolution()
    type(cluster2D_rejector) :: rjct
    type(image), allocatable :: imgs(:)
    type(oris)               :: cls_os
    logical, allocatable     :: rejected(:)
    write(*,'(A)') 'test_reject_res_flags_poor_resolution'
    call make_imgs(imgs, 4, 32, 2.0)
    call rjct%new(imgs, 50.0)
    call cls_os%new(4, is_ptcl=.false.)
    call cls_os%set(1, 'res', 15.0)
    call cls_os%set(2, 'res', 30.0)
    call cls_os%set(3, 'res', 40.0)
    call cls_os%set(4, 'res', 55.0)
    call rjct%reject_res(cls_os)
    rejected = rjct%get_rejected()
    call assert_false(rejected(1), 'reject_res keeps class with res = 15 Å')
    call assert_false(rejected(2), 'reject_res keeps class with res = 30 Å')
    call assert_false(rejected(3), 'reject_res keeps class with res = 40 Å (boundary, not >)')
    call assert_true(rejected(4),  'reject_res flags class with res = 55 Å')
    call cls_os%kill()
    call rjct%kill()
    call kill_imgs(imgs)
  end subroutine test_reject_res_flags_poor_resolution

  ! Case B — all classes at or below threshold → none rejected
  subroutine test_reject_res_keeps_good_resolution()
    type(cluster2D_rejector) :: rjct
    type(image), allocatable :: imgs(:)
    type(oris)               :: cls_os
    logical, allocatable     :: rejected(:)
    integer                  :: i
    write(*,'(A)') 'test_reject_res_keeps_good_resolution'
    call make_imgs(imgs, 3, 32, 2.0)
    call rjct%new(imgs, 50.0)
    call cls_os%new(3, is_ptcl=.false.)
    call cls_os%set(1, 'res', 10.0)
    call cls_os%set(2, 'res', 20.0)
    call cls_os%set(3, 'res', 38.0)
    call rjct%reject_res(cls_os)
    rejected = rjct%get_rejected()
    do i = 1, 3
      call assert_false(rejected(i), 'reject_res keeps all classes when all res <= 40')
    end do
    call cls_os%kill()
    call rjct%kill()
    call kill_imgs(imgs)
  end subroutine test_reject_res_keeps_good_resolution

  ! Case C — class at exactly RES_THRESHOLD (40.0) must not be rejected
  subroutine test_reject_res_boundary_at_threshold()
    type(cluster2D_rejector) :: rjct
    type(image), allocatable :: imgs(:)
    type(oris)               :: cls_os
    logical, allocatable     :: rejected(:)
    write(*,'(A)') 'test_reject_res_boundary_at_threshold'
    call make_imgs(imgs, 2, 32, 2.0)
    call rjct%new(imgs, 50.0)
    call cls_os%new(2, is_ptcl=.false.)
    call cls_os%set(1, 'res', 40.0)
    call cls_os%set(2, 'res', 40.001)   ! just above threshold
    call rjct%reject_res(cls_os)
    rejected = rjct%get_rejected()
    call assert_false(rejected(1), 'res = 40.0 is not > threshold, must not be rejected')
    call assert_true(rejected(2),  'res > 40.0 must be rejected')
    call cls_os%kill()
    call rjct%kill()
    call kill_imgs(imgs)
  end subroutine test_reject_res_boundary_at_threshold

  ! ---------------------------------------------------------------------------
  ! Accumulation: rejections from separate criteria stack additively.
  !
  ! pop pass:  pops = [0, 50, 50]  → class 1 rejected
  ! res pass:  res  = [10, 10, 55] → class 3 rejected
  ! Expected after both: rejected = [T, F, T]; class 2 remains kept.
  ! ---------------------------------------------------------------------------

  subroutine test_rejection_accumulates_across_criteria()
    type(cluster2D_rejector) :: rjct
    type(image), allocatable :: imgs(:)
    type(oris)               :: cls_os_pop, cls_os_res
    logical, allocatable     :: rejected(:)
    write(*,'(A)') 'test_rejection_accumulates_across_criteria'
    call make_imgs(imgs, 3, 32, 2.0)
    call rjct%new(imgs, 50.0)
    call cls_os_pop%new(3, is_ptcl=.false.)
    call cls_os_pop%set(1, 'pop',  0.0)
    call cls_os_pop%set(2, 'pop', 50.0)
    call cls_os_pop%set(3, 'pop', 50.0)
    call rjct%reject_pop(cls_os_pop)
    call cls_os_res%new(3, is_ptcl=.false.)
    call cls_os_res%set(1, 'res', 10.0)
    call cls_os_res%set(2, 'res', 10.0)
    call cls_os_res%set(3, 'res', 55.0)
    call rjct%reject_res(cls_os_res)
    rejected = rjct%get_rejected()
    call assert_true(rejected(1),  'class 1 remains rejected after both passes (pop criterion)')
    call assert_false(rejected(2), 'class 2 remains kept after both passes')
    call assert_true(rejected(3),  'class 3 newly rejected by res criterion')
    call cls_os_pop%kill()
    call cls_os_res%kill()
    call rjct%kill()
    call kill_imgs(imgs)
  end subroutine test_rejection_accumulates_across_criteria

  ! ---------------------------------------------------------------------------
  ! get_states consistency: states(i) = 0 iff rejected(i) = .true.
  ! ---------------------------------------------------------------------------

  subroutine test_get_states_mirrors_get_rejected()
    type(cluster2D_rejector) :: rjct
    type(image), allocatable :: imgs(:)
    type(oris)               :: cls_os
    logical, allocatable     :: rejected(:)
    integer, allocatable     :: states(:)
    integer                  :: i
    write(*,'(A)') 'test_get_states_mirrors_get_rejected'
    ! 4 classes, res = [10, 50, 10, 45] → classes 2 and 4 rejected
    call make_imgs(imgs, 4, 32, 2.0)
    call rjct%new(imgs, 50.0)
    call cls_os%new(4, is_ptcl=.false.)
    call cls_os%set(1, 'res', 10.0)
    call cls_os%set(2, 'res', 50.0)
    call cls_os%set(3, 'res', 10.0)
    call cls_os%set(4, 'res', 45.0)
    call rjct%reject_res(cls_os)
    rejected = rjct%get_rejected()
    states   = rjct%get_states()
    do i = 1, 4
      if( rejected(i) ) then
        call assert_int(0, states(i), 'state = 0 for a rejected class')
      else
        call assert_int(1, states(i), 'state = 1 for a kept class')
      end if
    end do
    call cls_os%kill()
    call rjct%kill()
    call kill_imgs(imgs)
  end subroutine test_get_states_mirrors_get_rejected

  ! ---------------------------------------------------------------------------
  ! reject_mask with all-zero images
  !
  ! An all-zero image → bp filter output is all-zero → Otsu threshold produces
  ! an all-zero binary image → find_ccs gives nccs = 0 → every class rejected.
  ! ---------------------------------------------------------------------------

  subroutine test_reject_mask_all_zero_images_all_rejected()
    type(cluster2D_rejector) :: rjct
    type(image), allocatable :: imgs(:)
    logical, allocatable     :: rejected(:)
    integer                  :: i
    write(*,'(A)') 'test_reject_mask_all_zero_images_all_rejected'
    call make_imgs(imgs, 3, 32, 2.0)
    call rjct%new(imgs, 50.0)
    call rjct%reject_mask()
    rejected = rjct%get_rejected()
    do i = 1, 3
      call assert_true(rejected(i), &
        'reject_mask rejects every zero-filled class (no surviving CCs)')
    end do
    call rjct%kill()
    call kill_imgs(imgs)
  end subroutine test_reject_mask_all_zero_images_all_rejected

  ! ---------------------------------------------------------------------------
  ! reject_local_variance with all-zero images
  !
  ! All-zero images produce scores_inside = 0 and scores_outside = 0 for every
  ! class. These are rejected unconditionally before z-scoring, so the outcome
  ! is deterministic: every class must be rejected.
  ! ---------------------------------------------------------------------------

  subroutine test_reject_local_variance_all_zero_all_rejected()
    type(cluster2D_rejector) :: rjct
    type(image), allocatable :: imgs(:)
    logical, allocatable     :: rejected(:)
    integer                  :: i
    write(*,'(A)') 'test_reject_local_variance_all_zero_all_rejected'
    call make_imgs(imgs, 3, 32, 2.0)
    call rjct%new(imgs, 50.0)
    call rjct%reject_local_variance()
    rejected = rjct%get_rejected()
    do i = 1, 3
      call assert_true(rejected(i), &
        'reject_local_variance rejects every zero-variance class (no signal in either region)')
    end do
    call rjct%kill()
    call kill_imgs(imgs)
  end subroutine test_reject_local_variance_all_zero_all_rejected

end module simple_cluster2D_rejector_tester
