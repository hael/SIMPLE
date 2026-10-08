!@descr: unit tests for the canonical sigma2 state files: transactions, grouping and recovery guards (simple_sigma2_state, simple_sigma2_state_file)
! The noise-power spectra of the likelihood objective live in one committed file per lineage,
! updated through a candidate that local ranges are merged into, groups reduced in and that is
! published by commit. Pinned: global and per-stack grouping end to end, the version-two
! file layout (no section after the particle spectra), rejection of missing and overlapping
! ranges and of a truncated candidate without touching the committed file, generation-scoped
! paths and the prepared update, and the reduction skipping a record with a non-positive shell.
module simple_sigma2_state_tester
use, intrinsic :: iso_fortran_env, only: int8, int32, int64, real32
use simple_string,            only: string
use simple_syslib,            only: del_file, file_exists
use simple_sigma2_state_file, only: sigma2_state_header, sigma2_state_init_header, &
    &sigma2_state_create_candidate, sigma2_state_write_local_range, sigma2_state_read_header, &
    &sigma2_state_read_local_range, sigma2_state_read_particles, sigma2_state_read_groups, &
    &sigma2_state_validate_file, SIGMA2_GROUP_GLOBAL, SIGMA2_GROUP_STACK, SIGMA2_PROV_PSPEC, &
    &SIGMA2_PROV_RESIDUAL, SIGMA2_STATE_COMMITTED
use simple_ptcl_layout,       only: ptcl_layout_digest
use simple_sigma2_state,      only: sigma2_state_merge_local_ranges, &
    &sigma2_state_reduce_groups, sigma2_state_validate_identity, sigma2_state_validate_science, &
    &sigma2_state_candidate_path, sigma2_state_commit, sigma2_state_prepare_update, &
    &sigma2_state_range_path, sigma2_state_next_generation
use simple_test_utils
implicit none
private
public :: run_all_sigma2_state_tests

contains

    subroutine run_all_sigma2_state_tests()
        write(*,'(A)') '**** running all sigma2 state tests ****'
        call test_policy('global', SIGMA2_GROUP_GLOBAL, 1, 4)
        call test_policy('group',  SIGMA2_GROUP_STACK,  2, 8)
        call test_particle_io_layout()
        call test_recovery_guards()
        call test_update_preparation()
        call test_invalid_record_skip()
        call test_estimate_available_on_disk()
    end subroutine run_all_sigma2_state_tests

    subroutine test_policy(prefix, grouping, ngroups, nptcls)
        character(len=*), intent(in) :: prefix
        integer(int32),   intent(in) :: grouping
        integer,          intent(in) :: ngroups, nptcls
        type(sigma2_state_header) :: header
        type(string) :: ranges(2), refs(2)
        real(real32), allocatable :: spectra(:,:)
        logical, allocatable :: scheduled(:), active(:)
        integer, allocatable :: eo(:), groups(:), stack_ids(:), stack_indices(:)
        character(len=128) :: candidate, committed, message
        integer(int64) :: digest, prefix_digest
        integer :: i, status, midpoint
        write(*,'(A)') 'test_policy '//trim(prefix)
        candidate = trim(prefix)//'_sigma2_state.next'
        committed = trim(prefix)//'_sigma2_state.bin'
        ranges(1) = trim(prefix)//'_range_1.bin'
        ranges(2) = trim(prefix)//'_range_2.bin'
        call cleanup(candidate, committed, ranges(1)%to_char(), ranges(2)%to_char())
        refs(1) = '/data/particles_a.mrcs'
        refs(2) = '/data/particles_b.mrcs'
        allocate(stack_ids(nptcls), stack_indices(nptcls))
        do i = 1, nptcls
            stack_ids(i) = merge(1, 2, i <= nptcls/2)
            stack_indices(i) = 1 + modulo(i-1, max(1,nptcls/2))
        enddo
        digest = ptcl_layout_digest('test-lineage', refs, stack_ids, stack_indices)
        call require(digest /= 0_int64, 'layout digest is nonzero')
        prefix_digest = ptcl_layout_digest('test-lineage', refs, stack_ids, stack_indices, nptcls-1)
        call require(prefix_digest /= 0_int64 .and. prefix_digest /= digest, 'prefix layout digest is distinct')
        call sigma2_state_init_header(header, 1, 3, nptcls, 8, 1.5, ngroups, grouping, &
            &1_int64, digest, SIGMA2_PROV_PSPEC)
        call sigma2_state_create_candidate(candidate, header, status, message)
        call require_ok(status, message)
        allocate(spectra(3,nptcls))
        do i = 1, nptcls
            spectra(:,i) = real([i, 2*i, 3*i], real32)
        enddo
        midpoint = nptcls/2
        call sigma2_state_write_local_range(ranges(1)%to_char(), 1_int64, digest, 1, &
            &spectra(:,:midpoint), 1, 3, status, message)
        call require_ok(status, message)
        call sigma2_state_write_local_range(ranges(2)%to_char(), 1_int64, digest, midpoint+1, &
            &spectra(:,midpoint+1:), 1, 3, status, message)
        call require_ok(status, message)
        allocate(scheduled(nptcls), active(nptcls), eo(nptcls), groups(nptcls))
        scheduled = .true.
        active = .true.
        do i = 1, nptcls
            eo(i) = modulo(i-1,2)
            if( grouping == SIGMA2_GROUP_GLOBAL )then
                groups(i) = 1
            else
                groups(i) = 1 + (i-1)/(nptcls/2)
            endif
        enddo
        call sigma2_state_merge_local_ranges(candidate, ranges, scheduled, status, message)
        call require_ok(status, message)
        call sigma2_state_reduce_groups(candidate, active, eo, groups, status, message)
        call require_ok(status, message)
        call sigma2_state_validate_science(candidate, active, eo, groups, status, message)
        call require_ok(status, message)
        call sigma2_state_commit(candidate, committed, active, eo, groups, status, message)
        call require_ok(status, message)
        call sigma2_state_validate_identity(committed, 8, 1.5, 1, 3, nptcls, digest, status, message, &
            &expected_state=SIGMA2_STATE_COMMITTED)
        call require_ok(status, message)
        call require(file_exists(committed), 'committed sigma2 state published')
        call require(.not. file_exists(candidate), 'candidate consumed by publication')
        call cleanup(candidate, committed, ranges(1)%to_char(), ranges(2)%to_char())
    end subroutine test_policy

    !> sigma2_estimate_available on a project file: the committed state of
    !! the project's own layout is consumable; another grid is not. The layout
    !! digest covers the stack table, so the project's stk segment must be
    !! read (it was not: every state was reported "layout digest undefined")
    subroutine test_estimate_available_on_disk()
        use simple_sp_project,       only: sp_project
        use simple_ptcl_layout,      only: ptcl_layout_digest
        use simple_sigma2_bootstrap, only: sigma2_estimate_available
        use simple_syslib,           only: simple_getcwd
        integer, parameter :: NP = 6, BOXT = 8
        real,    parameter :: SMPDT = 1.5
        character(len=*), parameter :: PROJ = 'tmp_sigma2_state_tester_proj.simple'
        character(len=*), parameter :: CAND = 'tmp_sigma2_state_tester.next'
        character(len=*), parameter :: COMM = 'tmp_sigma2_state_tester.bin'
        character(len=*), parameter :: RNG  = 'tmp_sigma2_state_tester_range.bin'
        type(sp_project) :: project
        type(sigma2_state_header) :: header
        type(string) :: ranges(1), cwd
        real(real32) :: spectra(BOXT/2, NP)
        logical :: scheduled(NP), active(NP)
        integer :: eo(NP), groups(NP), i, status
        integer(int64) :: digest
        character(len=128) :: message
        write(*,'(A)') 'test_estimate_available_on_disk'
        call cleanup(CAND, COMM, RNG, PROJ)
        call project%projinfo%new(1, is_ptcl=.false.)
        call project%projinfo%set(1, 'projname', 'sigma_lineage')
        call project%projinfo%set(1, 'projfile', PROJ)
        call simple_getcwd(cwd)
        call project%projinfo%set(1, 'cwd', cwd%to_char())
        call project%os_stk%new(1, is_ptcl=.false.)
        call project%os_stk%set(1, 'stk',   '/data/sigma_particles.mrcs')
        call project%os_stk%set(1, 'fromp', 1)
        call project%os_stk%set(1, 'top',   NP)
        call project%os_stk%set(1, 'nptcls_stk', NP)
        call project%os_stk%set(1, 'box',   BOXT)
        call project%os_stk%set(1, 'smpd',  SMPDT)
        call project%os_ptcl3D%new(NP, is_ptcl=.true.)
        do i = 1, NP
            call project%os_ptcl3D%set(i, 'stkind', 1)
            call project%os_ptcl3D%set(i, 'indstk', i)
            call project%os_ptcl3D%set_state(i, 1)
        enddo
        project%os_ptcl2D = project%os_ptcl3D
        digest = ptcl_layout_digest(project, project%os_ptcl3D)
        call require(digest /= 0_int64, 'the in-memory project has a layout digest')
        call sigma2_state_init_header(header, 1, BOXT/2, NP, BOXT, SMPDT, 1, SIGMA2_GROUP_GLOBAL, &
            &1_int64, digest, SIGMA2_PROV_RESIDUAL)
        call sigma2_state_create_candidate(CAND, header, status, message)
        call require_ok(status, message)
        do i = 1, NP
            spectra(:,i) = real(i, real32)
        enddo
        call sigma2_state_write_local_range(RNG, 1_int64, digest, 1, spectra, 1, BOXT/2, status, message)
        call require_ok(status, message)
        ranges(1) = RNG
        scheduled = .true.
        active    = .true.
        eo        = [(modulo(i-1,2), i = 1, NP)]
        groups    = 1
        call sigma2_state_merge_local_ranges(CAND, ranges, scheduled, status, message)
        call require_ok(status, message)
        call sigma2_state_reduce_groups(CAND, active, eo, groups, status, message)
        call require_ok(status, message)
        call sigma2_state_commit(CAND, COMM, active, eo, groups, status, message)
        call require_ok(status, message)
        call project%set_sigma2_state_path(string(COMM))
        call project%write(string(PROJ))
        call assert_true(sigma2_estimate_available(string(PROJ), BOXT, SMPDT, .true.), &
            &'a committed state of the project''s own layout is consumable from the project file')
        call assert_false(sigma2_estimate_available(string(PROJ), BOXT+2, SMPDT, .true.), &
            &'negative control: the same state is not consumable on another grid')
        call project%kill
        call cleanup(CAND, COMM, RNG, PROJ)
    end subroutine test_estimate_available_on_disk

    subroutine test_particle_io_layout()
        type(sigma2_state_header) :: header
        type(string) :: range_path
        real(real32), allocatable :: loaded(:,:)
        real(real32) :: spectra(3,4)
        logical :: scheduled(4), active(4)
        integer :: eo(4), groups(4), status
        integer :: first_row, last_row, kfrom, kto
        integer(int64), parameter :: DIGEST = 661_int64
        integer(int64), parameter :: RANGE_HEADER_BYTES = 256_int64
        integer(int64) :: generation, layout_digest, expected_bytes, actual_bytes
        character(len=128) :: message
        character(len=*), parameter :: CANDIDATE = 'layout_sigma2_state.next'
        character(len=*), parameter :: COMMITTED = 'layout_sigma2_state.bin'
        character(len=*), parameter :: RANGE = 'layout_sigma2_range.bin'
        write(*,'(A)') 'test_particle_io_layout'
        call del_file(CANDIDATE)
        call del_file(COMMITTED)
        call del_file(RANGE)
        spectra(:,1) = [1.0,2.0,3.0]
        spectra(:,2) = [2.0,3.0,4.0]
        spectra(:,3) = [3.0,4.0,5.0]
        spectra(:,4) = [4.0,5.0,6.0]
        scheduled = .true.; active = .true.; eo = [0,1,0,1]; groups = 1
        call sigma2_state_init_header(header, 1, 3, 4, 8, 1.5, 1, SIGMA2_GROUP_GLOBAL, &
            &1_int64, DIGEST, SIGMA2_PROV_PSPEC)
        call sigma2_state_create_candidate(CANDIDATE, header, status, message)
        call require_ok(status, message)
        range_path = RANGE
        call sigma2_state_write_local_range(RANGE, 1_int64, DIGEST, 1, spectra, 1, 3, status, message)
        call require_ok(status, message)
        ! a local range is its header and its spectra, nothing after them
        expected_bytes = RANGE_HEADER_BYTES + int(size(spectra),int64)*4_int64
        inquire(file=RANGE, size=actual_bytes, iostat=status)
        call require(status == 0 .and. actual_bytes == expected_bytes, &
            &'local range ends with its spectra (version-two layout)')
        call sigma2_state_read_local_range(RANGE, generation, layout_digest, first_row, last_row, &
            &kfrom, kto, loaded, status, message)
        call require_ok(status, message)
        call require(generation == 1_int64 .and. layout_digest == DIGEST, &
            &'local range retains its transaction identity')
        call require(first_row == 1 .and. last_row == 4 .and. kfrom == 1 .and. kto == 3, &
            &'local range retains its bounds')
        call require(all(loaded == spectra), 'local range round-trips its spectra')
        deallocate(loaded)

        call sigma2_state_merge_local_ranges(CANDIDATE, [range_path], scheduled, status, message)
        call require_ok(status, message)
        call sigma2_state_read_particles(CANDIDATE, 1, 4, loaded, status, message)
        call require_ok(status, message)
        call require(all(loaded == spectra), 'bulk particle write preserves every spectrum')
        deallocate(loaded)
        ! a state file is its header, its grouped section and its particle section
        call sigma2_state_read_header(CANDIDATE, header, status, message)
        call require_ok(status, message)
        inquire(file=CANDIDATE, size=actual_bytes, iostat=status)
        call require(status == 0 .and. actual_bytes == header%file_bytes, 'candidate size equals its recorded size')
        call require(header%file_bytes == header%particle_offset + 3_int64*4_int64*4_int64 - 1_int64, &
            &'the particle section is the last section of the state file')

        call sigma2_state_reduce_groups(CANDIDATE, active, eo, groups, status, message)
        call require_ok(status, message)
        call sigma2_state_validate_file(CANDIDATE, status, message, deep=.true.)
        call require_ok(status, message)
        call sigma2_state_commit(CANDIDATE, COMMITTED, active, eo, groups, status, message)
        call require_ok(status, message)
        call sigma2_state_read_particles(COMMITTED, 1, 4, loaded, status, message)
        call require_ok(status, message)
        call require(all(loaded == spectra), 'committed particle spectra round-trip')
        deallocate(loaded)
        call del_file(CANDIDATE)
        call del_file(COMMITTED)
        call del_file(RANGE)
        call range_path%kill
    end subroutine test_particle_io_layout

    subroutine test_recovery_guards()
        type(sigma2_state_header) :: header, committed_header
        type(string) :: ranges(2)
        real(real32) :: spectra(3,4)
        logical :: scheduled(4), active(4)
        integer :: eo(4), groups(4), status, funit
        integer(int8), allocatable :: committed_bytes(:)
        integer(int64), parameter :: DIGEST = 771_int64
        character(len=128) :: message
        character(len=*), parameter :: CANDIDATE='guard_sigma2_state.next'
        character(len=*), parameter :: COMMITTED='guard_sigma2_state.bin'
        write(*,'(A)') 'test_recovery_guards'
        ranges(1) = 'guard_range_1.bin'
        ranges(2) = 'guard_range_2.bin'
        call cleanup(CANDIDATE, COMMITTED, ranges(1)%to_char(), ranges(2)%to_char())
        spectra(:,1) = [1.0,2.0,3.0]
        spectra(:,2) = [2.0,3.0,4.0]
        spectra(:,3) = [3.0,4.0,5.0]
        spectra(:,4) = [4.0,5.0,6.0]
        scheduled = .true.; active = .true.; eo = [0,1,0,1]; groups = 1
        call sigma2_state_init_header(header, 1, 3, 4, 8, 1.5, 1, SIGMA2_GROUP_GLOBAL, &
            &1_int64, DIGEST, SIGMA2_PROV_PSPEC)
        call sigma2_state_create_candidate(CANDIDATE, header, status, message)
        call require_ok(status, message)
        call sigma2_state_write_local_range(ranges(1)%to_char(), 1_int64, DIGEST, 1, spectra, &
            &1, 3, status, message)
        call require_ok(status, message)
        call sigma2_state_merge_local_ranges(CANDIDATE, ranges(:1), scheduled, status, message)
        call require(status == 0, 'single complete range accepted')
        call sigma2_state_reduce_groups(CANDIDATE, active, eo, groups, status, message)
        call require_ok(status, message)
        call sigma2_state_commit(CANDIDATE, COMMITTED, active, eo, groups, status, message)
        call require_ok(status, message)
        call sigma2_state_read_header(COMMITTED, committed_header, status, message)
        call require_ok(status, message)
        call read_file_bytes(COMMITTED, committed_bytes)

        call del_file(ranges(1))
        header%generation = 2_int64
        call sigma2_state_create_candidate(CANDIDATE, header, status, message, source_path=COMMITTED)
        call require_ok(status, message)
        call sigma2_state_write_local_range(ranges(1)%to_char(), 2_int64, DIGEST, 1, spectra(:,:2), &
            &1, 3, status, message)
        call require_ok(status, message)
        call sigma2_state_merge_local_ranges(CANDIDATE, ranges(:1), scheduled, status, message)
        call require(status /= 0, 'missing scheduled range rejected')
        call assert_committed_generation(COMMITTED, committed_header%generation)
        call assert_file_unchanged(COMMITTED, committed_bytes)

        call del_file(CANDIDATE); call del_file(ranges(1))
        call sigma2_state_create_candidate(CANDIDATE, header, status, message, source_path=COMMITTED)
        call require_ok(status, message)
        call sigma2_state_write_local_range(ranges(1)%to_char(), 2_int64, DIGEST, 1, spectra(:,:3), &
            &1, 3, status, message)
        call require_ok(status, message)
        call sigma2_state_write_local_range(ranges(2)%to_char(), 2_int64, DIGEST, 3, spectra(:,3:), &
            &1, 3, status, message)
        call require_ok(status, message)
        call sigma2_state_merge_local_ranges(CANDIDATE, ranges, scheduled, status, message)
        call require(status /= 0, 'overlapping scheduled ranges rejected')
        call assert_committed_generation(COMMITTED, committed_header%generation)
        call assert_file_unchanged(COMMITTED, committed_bytes)

        call del_file(CANDIDATE)
        open(newunit=funit, file=CANDIDATE, access='stream', form='unformatted', status='replace')
        write(funit) 1
        close(funit)
        call sigma2_state_validate_file(CANDIDATE, status, message, deep=.true.)
        call require(status /= 0, 'truncated candidate rejected')
        call assert_committed_generation(COMMITTED, committed_header%generation)
        call assert_file_unchanged(COMMITTED, committed_bytes)
        call cleanup(CANDIDATE, COMMITTED, ranges(1)%to_char(), ranges(2)%to_char())
    end subroutine test_recovery_guards

    subroutine test_update_preparation()
        type(sigma2_state_header) :: header
        type(string) :: candidate_path, range_path
        real(real32) :: spectra(3,4)
        logical :: active(4)
        integer :: eo(4), groups(4), status
        integer(int64), parameter :: DIGEST = 991_int64
        integer(int64) :: next_gen
        character(len=128) :: message
        character(len=*), parameter :: COMMITTED = 'sigma2_state.bin'
        write(*,'(A)') 'test_update_preparation'
        candidate_path = sigma2_state_candidate_path(COMMITTED, 1_int64)
        range_path = sigma2_state_range_path(COMMITTED, 1_int64, 2, 3)
        call del_file(COMMITTED)
        call del_file(candidate_path)
        call del_file(range_path)
        call del_file(sigma2_state_candidate_path(COMMITTED, 2_int64))
        spectra(:,1) = [1.0,2.0,3.0]
        spectra(:,2) = [2.0,3.0,4.0]
        spectra(:,3) = [3.0,4.0,5.0]
        spectra(:,4) = [4.0,5.0,6.0]
        active = .true.; eo = [0,1,0,1]; groups = 1
        call sigma2_state_init_header(header, 1, 3, 4, 8, 1.5, 1, SIGMA2_GROUP_GLOBAL, &
            &1_int64, DIGEST, SIGMA2_PROV_PSPEC)
        call sigma2_state_create_candidate(candidate_path%to_char(), header, status, message)
        call require_ok(status, message)
        call sigma2_state_write_local_range(range_path%to_char(), 1_int64, DIGEST, 1, spectra, &
            &1, 3, status, message)
        call require_ok(status, message)
        call sigma2_state_merge_local_ranges(candidate_path%to_char(), [range_path], &
            &[.true.,.true.,.true.,.true.], status, message)
        call require_ok(status, message)
        call sigma2_state_reduce_groups(candidate_path%to_char(), active, eo, groups, status, message)
        call require_ok(status, message)
        call sigma2_state_commit(candidate_path%to_char(), COMMITTED, active, eo, groups, status, message)
        call require_ok(status, message)
        call require(index(candidate_path%to_char(), 'sigma2_state.g1.next') > 0, &
            &'candidate path is scoped to the generation it commits')
        call require(index(range_path%to_char(), 'sigma2_state.g1.part002.range') > 0, &
            &'canonical range path is scoped to the generation and includes the padded partition')
        ! the next transaction is named for the generation it will commit
        call sigma2_state_next_generation(COMMITTED, next_gen, status, message)
        call require_ok(status, message)
        call require(next_gen == 2_int64, 'next generation follows the committed one')
        candidate_path = sigma2_state_candidate_path(COMMITTED, next_gen)
        call require(index(candidate_path%to_char(), 'sigma2_state.g2.next') > 0, &
            &'next candidate path is scoped to the next generation')
        call sigma2_state_prepare_update(COMMITTED, candidate_path%to_char(), status, message)
        call require_ok(status, message)
        call sigma2_state_read_header(candidate_path%to_char(), header, status, message)
        call require_ok(status, message)
        call require(header%generation == 2_int64, 'prepared update advances the generation')
        call require(header%provenance == SIGMA2_PROV_RESIDUAL, 'prepared update records residual provenance')
        call assert_committed_generation(COMMITTED, 1_int64)
        call del_file(COMMITTED)
        call del_file(candidate_path)
        call del_file(range_path)
    end subroutine test_update_preparation

    !> an active record with a non-positive shell is skipped by the
    !! reduction (with a warning) instead of aborting; the grouped model
    !! is the mean of the valid records and commit-time validation agrees
    subroutine test_invalid_record_skip()
        type(sigma2_state_header) :: header
        type(string) :: candidate_path, range_path
        real(real32), allocatable :: stored(:,:,:)
        real(real32) :: spectra(3,4), expected_even(3), expected_odd(3)
        logical :: active(4)
        integer :: eo(4), groups(4), status
        integer(int64), parameter :: DIGEST = 1231_int64
        character(len=128) :: message
        character(len=*), parameter :: COMMITTED = 'skip_sigma2_state.bin'
        write(*,'(A)') 'test_invalid_record_skip'
        candidate_path = sigma2_state_candidate_path(COMMITTED, 1_int64)
        range_path     = sigma2_state_range_path(COMMITTED, 1_int64, 1, 1)
        call del_file(COMMITTED)
        call del_file(candidate_path)
        call del_file(range_path)
        spectra(:,1) = [1.0,2.0,3.0]
        spectra(:,2) = [2.0,3.0,4.0]
        spectra(:,3) = [3.0,0.0,5.0] ! invalid: a non-positive shell
        spectra(:,4) = [4.0,5.0,6.0]
        active = .true.; eo = [0,1,0,1]; groups = 1
        expected_even = spectra(:,1)
        expected_odd  = 0.5*(spectra(:,2)+spectra(:,4))
        call sigma2_state_init_header(header, 1, 3, 4, 8, 1.5, 1, SIGMA2_GROUP_GLOBAL, &
            &1_int64, DIGEST, SIGMA2_PROV_PSPEC)
        call sigma2_state_create_candidate(candidate_path%to_char(), header, status, message)
        call require_ok(status, message)
        call sigma2_state_write_local_range(range_path%to_char(), 1_int64, DIGEST, 1, spectra, &
            &1, 3, status, message)
        call require_ok(status, message)
        call sigma2_state_merge_local_ranges(candidate_path%to_char(), [range_path], &
            &[.true.,.true.,.true.,.true.], status, message)
        call require_ok(status, message)
        call sigma2_state_reduce_groups(candidate_path%to_char(), active, eo, groups, status, message)
        call require(status == 0, 'invalid record is skipped, not fatal')
        call sigma2_state_read_groups(candidate_path%to_char(), stored, status, message)
        call require_ok(status, message)
        call require(all(abs(stored(:,1,1)-expected_even) <= 1.e-6*expected_even), &
            &'even group mean excludes the invalid record')
        call require(all(abs(stored(:,2,1)-expected_odd) <= 1.e-6*expected_odd), &
            &'odd group mean is the mean of the valid records')
        call sigma2_state_validate_science(candidate_path%to_char(), active, eo, groups, status, message)
        call require(status == 0, 'commit-time validation applies the same skip rule')
        call sigma2_state_commit(candidate_path%to_char(), COMMITTED, active, eo, groups, status, message)
        call require_ok(status, message)
        call del_file(COMMITTED)
        call del_file(candidate_path)
        call del_file(range_path)
    end subroutine test_invalid_record_skip

    subroutine assert_committed_generation(path, expected)
        character(len=*), intent(in) :: path
        integer(int64),   intent(in) :: expected
        type(sigma2_state_header) :: header
        character(len=128) :: message
        integer :: status
        call sigma2_state_read_header(path, header, status, message)
        call require_ok(status, message)
        call require(header%generation == expected, 'candidate transaction preserves committed generation until publication')
    end subroutine assert_committed_generation

    subroutine assert_file_unchanged(path, expected)
        character(len=*), intent(in) :: path
        integer(int8),    intent(in) :: expected(:)
        integer(int8), allocatable :: actual(:)
        call read_file_bytes(path, actual)
        call require(size(actual) == size(expected), 'failed candidate preserves committed file size')
        call require(all(actual == expected), 'failed candidate preserves committed bytes')
    end subroutine assert_file_unchanged

    subroutine read_file_bytes(path, bytes)
        character(len=*), intent(in) :: path
        integer(int8), allocatable, intent(out) :: bytes(:)
        integer(int64) :: file_bytes
        integer :: funit, io_stat
        inquire(file=path, size=file_bytes, iostat=io_stat)
        call require(io_stat == 0 .and. file_bytes > 0_int64, 'committed file can be sized')
        allocate(bytes(file_bytes))
        open(newunit=funit, file=path, access='stream', form='unformatted', status='old', &
            &action='read', iostat=io_stat)
        call require(io_stat == 0, 'committed file can be opened')
        read(funit, pos=1, iostat=io_stat) bytes
        close(funit)
        call require(io_stat == 0, 'committed file can be read')
    end subroutine read_file_bytes

    subroutine cleanup(path1, path2, path3, path4)
        character(len=*), intent(in) :: path1, path2, path3, path4
        call del_file(path1); call del_file(path2); call del_file(path3); call del_file(path4)
    end subroutine cleanup

    !> a state operation returns status 0; on failure its own message names what went wrong
    subroutine require_ok(status, message)
        integer,          intent(in) :: status
        character(len=*), intent(in) :: message
        call assert_int(0, status, 'sigma2 state operation succeeds: '//trim(message))
    end subroutine require_ok

    subroutine require(condition, message)
        logical,          intent(in) :: condition
        character(len=*), intent(in) :: message
        call assert_true(condition, message)
    end subroutine require

end module simple_sigma2_state_tester
