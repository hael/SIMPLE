!@descr: unit tests for the microchunked2D type (simple_microchunked2D module)
!
! Testable surface
! ----------------
! Most behaviour of microchunked2D requires a live parameters object, a queue
! environment, and on-disk project files — these are integration concerns and
! are intentionally excluded here. The tests below exercise only the subset
! of the public interface that is self-contained:
!
!   • Query functions on a freshly declared or kill()ed object:
!       get_n_microchunks_pass_{1,2,match}, get_n_chunks_running,
!       get_n_pass_{1,2}_non_rejected_ptcls,
!       get_n_accepted_ptcls, get_n_rejected_ptcls, get_finished
!   • kill() safety: on a default-declared object and when called twice
!   • get_references:  .false. return contract and inout-array cleanup
!   • get_latest_match: same contract as get_references
!   • get_reference_selection: .false. return contract
!
! chunk2D is a module-private type so chunks cannot be constructed here;
! any test requiring actual chunk state belongs in an integration test suite.

module simple_microchunked2D_tester
  use simple_microchunked2D, only: microchunked2D
  use simple_test_utils
  use simple_string,         only: string
  implicit none

  public :: run_all_microchunked2D_tests
  private

contains

  subroutine run_all_microchunked2D_tests()
    write(*,'(A)') '**** running microchunked2D tests ****'
    call test_initial_chunk_counts()
    call test_initial_running_count()
    call test_initial_ptcl_counters()
    call test_initial_non_rejected_ptcl_counters()
    call test_get_finished_false_when_empty()
    call test_kill_on_default_object()
    call test_kill_idempotent()
    call test_get_references_false_on_fresh()
    call test_get_references_cleans_inout_arrays()
    call test_get_latest_match_false_on_fresh()
    call test_get_latest_match_cleans_inout_arrays()
    call test_get_reference_selection_false_on_fresh()
  end subroutine run_all_microchunked2D_tests

  ! ---------------------------------------------------------------------------
  ! Chunk count queries on a freshly killed object
  ! ---------------------------------------------------------------------------

  subroutine test_initial_chunk_counts()
    type(microchunked2D) :: obj
    write(*,'(A)') 'test_initial_chunk_counts'
    call obj%kill()
    call assert_int(0, obj%get_n_microchunks_pass_1(), &
      'get_n_microchunks_pass_1 = 0 after kill')
    call assert_int(0, obj%get_n_microchunks_pass_2(), &
      'get_n_microchunks_pass_2 = 0 after kill')
    call assert_int(0, obj%get_n_microchunks_match(), &
      'get_n_microchunks_match = 0 after kill')
  end subroutine test_initial_chunk_counts

  ! ---------------------------------------------------------------------------

  subroutine test_initial_running_count()
    type(microchunked2D) :: obj
    write(*,'(A)') 'test_initial_running_count'
    call obj%kill()
    call assert_int(0, obj%get_n_chunks_running(), &
      'get_n_chunks_running = 0 after kill')
  end subroutine test_initial_running_count

  ! ---------------------------------------------------------------------------

  subroutine test_initial_ptcl_counters()
    type(microchunked2D) :: obj
    write(*,'(A)') 'test_initial_ptcl_counters'
    call obj%kill()
    call assert_int(0, obj%get_n_accepted_ptcls(), &
      'get_n_accepted_ptcls = 0 after kill')
    call assert_int(0, obj%get_n_rejected_ptcls(), &
      'get_n_rejected_ptcls = 0 after kill')
  end subroutine test_initial_ptcl_counters

  ! ---------------------------------------------------------------------------

  subroutine test_initial_non_rejected_ptcl_counters()
    type(microchunked2D) :: obj
    write(*,'(A)') 'test_initial_non_rejected_ptcl_counters'
    call obj%kill()
    call assert_int(0, obj%get_n_pass_1_non_rejected_ptcls(), &
      'get_n_pass_1_non_rejected_ptcls = 0 after kill')
    call assert_int(0, obj%get_n_pass_2_non_rejected_ptcls(), &
      'get_n_pass_2_non_rejected_ptcls = 0 after kill')
  end subroutine test_initial_non_rejected_ptcl_counters

  ! ---------------------------------------------------------------------------
  ! get_finished: must return .false. when no pass-1 chunks exist
  ! ---------------------------------------------------------------------------

  subroutine test_get_finished_false_when_empty()
    type(microchunked2D) :: obj
    write(*,'(A)') 'test_get_finished_false_when_empty'
    call obj%kill()
    call assert_false(obj%get_finished(), &
      'get_finished = .false. when no chunks exist')
  end subroutine test_get_finished_false_when_empty

  ! ---------------------------------------------------------------------------
  ! kill() safety
  ! ---------------------------------------------------------------------------

  subroutine test_kill_on_default_object()
    type(microchunked2D) :: obj
    write(*,'(A)') 'test_kill_on_default_object'
    ! kill() on a never-initialised object must not crash; verify queries
    ! return defined values afterwards
    call obj%kill()
    call assert_int(0, obj%get_n_microchunks_pass_1(), &
      'pass_1 count = 0 after first kill on default object')
    call assert_int(0, obj%get_n_accepted_ptcls(), &
      'accepted count = 0 after first kill on default object')
  end subroutine test_kill_on_default_object

  ! ---------------------------------------------------------------------------

  subroutine test_kill_idempotent()
    type(microchunked2D) :: obj
    write(*,'(A)') 'test_kill_idempotent'
    call obj%kill()
    call obj%kill()  ! second kill must be safe
    call assert_int(0, obj%get_n_microchunks_pass_1(), &
      'pass_1 count = 0 after double kill')
    call assert_false(obj%get_finished(), &
      'get_finished = .false. after double kill')
  end subroutine test_kill_idempotent

  ! ---------------------------------------------------------------------------
  ! get_references: .false. return contract
  !   When the refchunk is not yet complete, get_references must:
  !     • return .false.
  !     • deallocate any pre-allocated inout arrays passed by the caller
  !     • zero the scalar intent(out) arguments
  ! ---------------------------------------------------------------------------

  subroutine test_get_references_false_on_fresh()
    type(microchunked2D) :: obj
    integer, allocatable  :: jpeg_inds(:), jpeg_pops(:)
    real,    allocatable  :: jpeg_res(:)
    type(string)          :: jpeg, stk
    integer               :: xtiles, ytiles
    logical               :: ok
    write(*,'(A)') 'test_get_references_false_on_fresh'
    call obj%kill()
    ok = obj%get_references(jpeg_inds, jpeg_pops, jpeg_res, jpeg, stk, xtiles, ytiles)
    call assert_false(ok,    'get_references returns .false. when refchunk not complete')
    call assert_int(0, xtiles, 'get_references zeros xtiles on .false. return')
    call assert_int(0, ytiles, 'get_references zeros ytiles on .false. return')
  end subroutine test_get_references_false_on_fresh

  ! ---------------------------------------------------------------------------

  subroutine test_get_references_cleans_inout_arrays()
    type(microchunked2D) :: obj
    integer, allocatable  :: jpeg_inds(:), jpeg_pops(:)
    real,    allocatable  :: jpeg_res(:)
    type(string)          :: jpeg, stk
    integer               :: xtiles, ytiles
    logical               :: ok
    write(*,'(A)') 'test_get_references_cleans_inout_arrays'
    call obj%kill()
    ! pre-allocate the inout arrays to verify they are released on .false. return
    allocate(jpeg_inds(10), jpeg_pops(10), jpeg_res(10))
    ok = obj%get_references(jpeg_inds, jpeg_pops, jpeg_res, jpeg, stk, xtiles, ytiles)
    call assert_false(ok, 'get_references .false. even with pre-allocated arrays')
    call assert_false(allocated(jpeg_inds), &
      'get_references deallocates jpeg_inds on .false. return')
    call assert_false(allocated(jpeg_pops), &
      'get_references deallocates jpeg_pops on .false. return')
    call assert_false(allocated(jpeg_res), &
      'get_references deallocates jpeg_res on .false. return')
  end subroutine test_get_references_cleans_inout_arrays

  ! ---------------------------------------------------------------------------
  ! get_latest_match: .false. return contract (mirrors get_references)
  ! ---------------------------------------------------------------------------

  subroutine test_get_latest_match_false_on_fresh()
    type(microchunked2D) :: obj
    integer, allocatable  :: jpeg_inds(:), jpeg_pops(:)
    real,    allocatable  :: jpeg_res(:)
    type(string)          :: jpeg, stk
    integer               :: xtiles, ytiles
    logical               :: ok
    write(*,'(A)') 'test_get_latest_match_false_on_fresh'
    call obj%kill()
    ok = obj%get_latest_match(jpeg_inds, jpeg_pops, jpeg_res, jpeg, stk, xtiles, ytiles)
    call assert_false(ok, 'get_latest_match returns .false. when no match available')
    call assert_int(0, xtiles, 'get_latest_match zeros xtiles on .false. return')
    call assert_int(0, ytiles, 'get_latest_match zeros ytiles on .false. return')
  end subroutine test_get_latest_match_false_on_fresh

  ! ---------------------------------------------------------------------------

  subroutine test_get_latest_match_cleans_inout_arrays()
    type(microchunked2D) :: obj
    integer, allocatable  :: jpeg_inds(:), jpeg_pops(:)
    real,    allocatable  :: jpeg_res(:)
    type(string)          :: jpeg, stk
    integer               :: xtiles, ytiles
    logical               :: ok
    write(*,'(A)') 'test_get_latest_match_cleans_inout_arrays'
    call obj%kill()
    allocate(jpeg_inds(8), jpeg_pops(8), jpeg_res(8))
    ok = obj%get_latest_match(jpeg_inds, jpeg_pops, jpeg_res, jpeg, stk, xtiles, ytiles)
    call assert_false(ok, 'get_latest_match .false. even with pre-allocated arrays')
    call assert_false(allocated(jpeg_inds), &
      'get_latest_match deallocates jpeg_inds on .false. return')
    call assert_false(allocated(jpeg_pops), &
      'get_latest_match deallocates jpeg_pops on .false. return')
    call assert_false(allocated(jpeg_res), &
      'get_latest_match deallocates jpeg_res on .false. return')
  end subroutine test_get_latest_match_cleans_inout_arrays

  ! ---------------------------------------------------------------------------
  ! get_reference_selection: .false. return and unallocated output on fresh object
  ! ---------------------------------------------------------------------------

  subroutine test_get_reference_selection_false_on_fresh()
    type(microchunked2D) :: obj
    integer, allocatable  :: selection(:)
    logical               :: ok
    write(*,'(A)') 'test_get_reference_selection_false_on_fresh'
    call obj%kill()
    ok = obj%get_reference_selection(selection)
    call assert_false(ok, 'get_reference_selection returns .false. before refchunk completes')
    call assert_false(allocated(selection), &
      'get_reference_selection leaves selection unallocated on .false. return')
  end subroutine test_get_reference_selection_false_on_fresh

end module simple_microchunked2D_tester
