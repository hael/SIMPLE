!@descr: validates canonical sigma2 transactions, grouping and recovery guards
program simple_test_sigma2_state
use, intrinsic :: iso_fortran_env, only: int8, int32, int64, real32
use simple_string,            only: string
use simple_syslib,            only: del_file, file_exists
use simple_sigma2_state_file, only: sigma2_state_header, sigma2_state_init_header, &
    &sigma2_state_create_candidate, sigma2_state_write_local_range, sigma2_state_read_header, &
    &sigma2_state_read_groups, sigma2_state_validate_file, SIGMA2_GROUP_GLOBAL, SIGMA2_GROUP_STACK, &
    &SIGMA2_PROV_PSPEC, SIGMA2_PROV_RESIDUAL, SIGMA2_STATE_COMMITTED
use simple_sigma2_state,      only: sigma2_state_layout_digest, sigma2_state_merge_local_ranges, &
    &sigma2_state_reduce_groups, sigma2_state_validate_identity, sigma2_state_validate_science, &
    &sigma2_state_candidate_path, sigma2_state_commit, sigma2_state_prepare_update, &
    &sigma2_state_range_path, sigma2_state_next_generation
implicit none

call exercise_policy('global', SIGMA2_GROUP_GLOBAL, 1, 4)
call exercise_policy('group',  SIGMA2_GROUP_STACK,  2, 8)
call exercise_recovery_guards()
call exercise_update_preparation()
call exercise_invalid_record_skip()
write(*,'(A)') 'SIMPLE_TEST_SIGMA2_STATE NORMAL STOP'

contains

    subroutine exercise_policy(prefix, grouping, ngroups, nptcls)
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
        digest = sigma2_state_layout_digest('test-lineage', refs, stack_ids, stack_indices)
        call require(digest /= 0_int64, 'layout digest is nonzero')
        prefix_digest = sigma2_state_layout_digest('test-lineage', refs, stack_ids, stack_indices, nptcls-1)
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
    end subroutine exercise_policy

    subroutine exercise_recovery_guards()
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
    end subroutine exercise_recovery_guards

    subroutine exercise_update_preparation()
        type(sigma2_state_header) :: header
        type(string) :: candidate_path, range_path
        real(real32) :: spectra(3,4)
        logical :: active(4)
        integer :: eo(4), groups(4), status
        integer(int64), parameter :: DIGEST = 991_int64
        integer(int64) :: next_gen
        character(len=128) :: message
        character(len=*), parameter :: COMMITTED = 'sigma2_state.bin'
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
        call del_file(COMMITTED)
        call del_file(candidate_path)
        call del_file(range_path)
    end subroutine exercise_update_preparation

    !> an active record with a non-positive shell is skipped by the
    !! reduction (with a warning) instead of aborting; the grouped model
    !! is the mean of the valid records and commit-time validation agrees
    subroutine exercise_invalid_record_skip()
        type(sigma2_state_header) :: header
        type(string) :: candidate_path, range_path
        real(real32), allocatable :: stored(:,:,:)
        real(real32) :: spectra(3,4), expected_even(3), expected_odd(3)
        logical :: active(4)
        integer :: eo(4), groups(4), status
        integer(int64), parameter :: DIGEST = 1231_int64
        character(len=128) :: message
        character(len=*), parameter :: COMMITTED = 'skip_sigma2_state.bin'
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
    end subroutine exercise_invalid_record_skip

    subroutine assert_committed_generation(path, expected)
        character(len=*), intent(in) :: path
        integer(int64),   intent(in) :: expected
        type(sigma2_state_header) :: header
        character(len=128) :: message
        integer :: status
        call sigma2_state_read_header(path, header, status, message)
        call require_ok(status, message)
        call require(header%generation == expected, 'failed candidate preserves committed generation')
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

    subroutine require_ok(status, message)
        integer,          intent(in) :: status
        character(len=*), intent(in) :: message
        if( status /= 0 )then
            write(*,'(A)') trim(message)
            error stop 'canonical sigma2 state operation failed'
        endif
    end subroutine require_ok

    subroutine require(condition, message)
        logical,          intent(in) :: condition
        character(len=*), intent(in) :: message
        if( .not. condition )then
            write(*,'(A)') trim(message)
            error stop 'canonical sigma2 state assertion failed'
        endif
    end subroutine require

end program simple_test_sigma2_state
