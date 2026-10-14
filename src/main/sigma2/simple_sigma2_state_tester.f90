!@descr: unit tests for the canonical sigma2 state files: transactions, the global spectra and recovery guards (simple_sigma2_state, simple_sigma2_state_file)
! The noise-power spectra of the likelihood objective live in one committed file per lineage. An update
! streams a fresh candidate from the workers' local ranges and the committed store, forming the two global
! spectra (one per even/odd half) as the rows pass, and publishes it. Pinned: the streaming commit end to
! end, its global spectra against a full reduction of the records, rows outside the schedule taken from the
! committed store, the refusal of the retired header words (stack grouping, STAR-seed provenance, direct
! range mode), the retired checksum words (zero in new files, ignored on read), the version-two file
! layout (no section after the particle spectra), refusal of uncovered
! and doubly covered rows and of a failed publication without touching the committed file or the ranges,
! generation-scoped paths and the prepared update, the block-copying concatenation, and the skipping of a
! record with a non-positive shell.
module simple_sigma2_state_tester
use, intrinsic :: iso_fortran_env, only: int8, int64, real32, real64
use simple_fileio,            only: simple_copy_file
use simple_string,            only: string
use simple_syslib,            only: del_file, file_exists, simple_mkdir, simple_rmdir
use simple_sigma2_state_file, only: sigma2_state_header, sigma2_state_init_header, &
    &sigma2_state_create_candidate, sigma2_state_write_local_range, sigma2_state_read_header, &
    &sigma2_state_read_local_range, sigma2_state_read_particles, sigma2_state_read_global, &
    &sigma2_state_validate_file, sigma2_state_publish, sigma2_state_write_particles, &
    &sigma2_state_write_global, SIGMA2_PROV_PSPEC, &
    &SIGMA2_PROV_RESIDUAL, SIGMA2_STATE_COMMITTED
use simple_ptcl_layout,       only: ptcl_layout_digest
use simple_sigma2_state,      only: sigma2_state_assemble_candidate, sigma2_state_concatenate, &
    &sigma2_state_reduce, sigma2_state_validate_identity, &
    &sigma2_state_candidate_path, sigma2_state_prepare_update, &
    &sigma2_state_range_path, sigma2_state_next_generation
use simple_test_utils
implicit none
private
public :: run_all_sigma2_state_tests, write_test_sigma2_store, patch_sigma2_header_word

contains

    subroutine run_all_sigma2_state_tests()
        write(*,'(A)') '**** running all sigma2 state tests ****'
        call test_global_transaction(4)
        call test_global_transaction(9000)
        call test_unscheduled_rows_from_committed()
        call test_retired_header_words()
        call test_particle_io_layout()
        call test_recovery_guards()
        call test_failed_publication()
        call test_update_preparation()
        call test_concatenation()
        call test_invalid_record_skip()
        call test_estimate_available_on_disk()
    end subroutine run_all_sigma2_state_tests

    !> The streaming commit end to end: two ranges covering every row of a committed generation-1 store
    !! assemble the generation-2 candidate, whose global spectra equal a full reduction of its records;
    !! publication replaces the store and leaves the ranges for the caller. 9000 rows span blocks.
    subroutine test_global_transaction(nptcls)
        integer,          intent(in) :: nptcls
        character(len=*), parameter  :: prefix = 'global'
        type(sigma2_state_header) :: header
        type(string) :: ranges(2), refs(2)
        real(real32), allocatable :: spectra(:,:), loaded(:,:)
        logical, allocatable :: scheduled(:), active(:)
        integer, allocatable :: eo(:), stack_ids(:), stack_indices(:)
        character(len=128) :: candidate, committed, message
        integer(int64) :: digest, prefix_digest
        integer :: i, status, midpoint
        write(*,'(A)') 'test_global_transaction'
        committed = trim(prefix)//'_sigma2_state.bin'
        candidate = trim(prefix)//'_sigma2_state.next'
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
        allocate(spectra(3,nptcls), scheduled(nptcls), active(nptcls), eo(nptcls))
        scheduled = .true.
        do i = 1, nptcls
            active(i) = modulo(i,7) /= 0
            eo(i)     = modulo(i-1,2)
        enddo
        call sigma2_state_init_header(header, 1, 3, nptcls, 8, 1.5, 1_int64, digest, SIGMA2_PROV_PSPEC)
        spectra = 1.0_real32
        call write_committed_store(committed, header, spectra, active, eo)
        do i = 1, nptcls
            spectra(:,i) = real([i, 2*i, 3*i], real32)
        enddo
        midpoint = nptcls/2
        call sigma2_state_write_local_range(ranges(1)%to_char(), 2_int64, digest, 1, &
            &spectra(:,:midpoint), 1, 3, status, message)
        call require_ok(status, message)
        call sigma2_state_write_local_range(ranges(2)%to_char(), 2_int64, digest, midpoint+1, &
            &spectra(:,midpoint+1:), 1, 3, status, message)
        call require_ok(status, message)
        call sigma2_state_assemble_candidate(committed, candidate, ranges, scheduled, active, eo, status, message)
        call require_ok(status, message)
        call sigma2_state_read_header(candidate, header, status, message)
        call require_ok(status, message)
        call require(header%generation == 2_int64 .and. header%provenance == SIGMA2_PROV_RESIDUAL, &
            &'the streamed candidate is the next generation, of residual provenance')
        call sigma2_state_read_particles(candidate, 1, nptcls, loaded, status, message)
        call require_ok(status, message)
        call require(all(loaded == spectra), 'the streamed candidate holds the rows of the ranges')
        call require_full_reduction(candidate, active, eo)
        call sigma2_state_publish(candidate, committed, status, message)
        call require_ok(status, message)
        call sigma2_state_validate_identity(committed, 8, 1.5, 1, 3, nptcls, digest, status, message, &
            &expected_state=SIGMA2_STATE_COMMITTED)
        call require_ok(status, message)
        call assert_committed_generation(committed, 2_int64)
        call require(.not. file_exists(candidate), 'candidate consumed by publication')
        call require(file_exists(ranges(1)%to_char()) .and. file_exists(ranges(2)%to_char()), &
            &'the ranges stay until the caller deletes them after publication')
        call cleanup(candidate, committed, ranges(1)%to_char(), ranges(2)%to_char())
    end subroutine test_global_transaction

    !> A row outside the schedule is taken from the committed store, and the global spectra still equal a
    !! full reduction of the candidate's records
    subroutine test_unscheduled_rows_from_committed()
        integer,          parameter :: NP = 6
        integer(int64),   parameter :: DIGEST = 551_int64
        character(len=*), parameter :: COMMITTED = 'partial_sigma2_state.bin'
        character(len=*), parameter :: RANGE = 'partial_sigma2_range.bin'
        type(sigma2_state_header) :: header
        real(real32), allocatable :: loaded(:,:)
        real(real32) :: old(3,NP), new(3,NP)
        logical :: scheduled(NP), active(NP)
        integer :: eo(NP), i, status
        character(len=128) :: candidate, message
        write(*,'(A)') 'test_unscheduled_rows_from_committed'
        candidate = 'partial_sigma2_state.next'
        call cleanup(candidate, COMMITTED, RANGE, RANGE)
        do i = 1, NP
            old(:,i) = real([10*i, 10*i+1, 10*i+2], real32)
            new(:,i) = real([i, i+1, i+2], real32)
        enddo
        active = .true.
        eo     = [(modulo(i-1,2), i = 1, NP)]
        call sigma2_state_init_header(header, 1, 3, NP, 8, 1.5, 1_int64, DIGEST, SIGMA2_PROV_PSPEC)
        call write_committed_store(COMMITTED, header, old, active, eo)
        ! rows 3:4 scheduled and covered; the others come from the committed store
        scheduled = .false.
        scheduled(3:4) = .true.
        call sigma2_state_write_local_range(RANGE, 2_int64, DIGEST, 3, new(:,3:4), 1, 3, status, message)
        call require_ok(status, message)
        call sigma2_state_assemble_candidate(COMMITTED, candidate, [string(RANGE)], scheduled, active, eo, &
            &status, message)
        call require_ok(status, message)
        call sigma2_state_read_particles(candidate, 1, NP, loaded, status, message)
        call require_ok(status, message)
        call require(all(loaded(:,3:4) == new(:,3:4)), 'scheduled rows come from their range')
        call require(all(loaded(:,1:2) == old(:,1:2)) .and. all(loaded(:,5:) == old(:,5:)), &
            &'unscheduled rows come from the committed store')
        call require_full_reduction(candidate, active, eo)
        call cleanup(candidate, COMMITTED, RANGE, RANGE)
    end subroutine test_unscheduled_rows_from_committed

    !> A store whose header carries a retired word is refused at read
    !! and validation: the stack grouping policy, two groups, the STAR-seed provenance (3) and the
    !! direct range mode (2). Each is patched into a valid committed store.
    subroutine test_retired_header_words()
        integer, parameter :: W_NGROUPS=8, W_GROUPING=9, W_PROVENANCE=11, W_RANGE_IO=13
        character(len=*), parameter :: COMMITTED = 'retired_sigma2_state.bin'
        type(sigma2_state_header) :: header
        type(string) :: candidate
        real(real32) :: spectra(3,4), globals(3,2)
        character(len=128) :: message
        integer :: status
        write(*,'(A)') 'test_retired_header_words'
        candidate = trim(COMMITTED)//'.next'
        call write_valid_store
        call sigma2_state_read_header(COMMITTED, header, status, message)
        call require_ok(status, message)
        call patch_sigma2_header_word(COMMITTED, W_GROUPING, 2_int64)
        call sigma2_state_read_header(COMMITTED, header, status, message)
        call require(status /= 0 .and. index(message, 'grouping') > 0, 'the stack grouping policy is refused')
        call write_valid_store
        call patch_sigma2_header_word(COMMITTED, W_NGROUPS, 2_int64)
        call sigma2_state_validate_file(COMMITTED, status, message)
        call require(status /= 0, 'more than one group is refused')
        call write_valid_store
        call patch_sigma2_header_word(COMMITTED, W_PROVENANCE, 3_int64)
        call sigma2_state_read_header(COMMITTED, header, status, message)
        call require(status /= 0 .and. index(message, 'provenance') > 0, 'the STAR-seed provenance is refused')
        call write_valid_store
        call patch_sigma2_header_word(COMMITTED, W_RANGE_IO, 2_int64)
        call sigma2_state_read_header(COMMITTED, header, status, message)
        call require(status /= 0 .and. index(message, 'range') > 0, 'the direct range mode is refused')
        call del_file(COMMITTED)
        call del_file(candidate)
        call candidate%kill

      contains

        subroutine write_valid_store
            call del_file(COMMITTED)
            call del_file(candidate)
            call sigma2_state_init_header(header, 1, 3, 4, 8, 1.5, 1_int64, 4242_int64, SIGMA2_PROV_RESIDUAL)
            call sigma2_state_create_candidate(candidate%to_char(), header, status, message)
            call require_ok(status, message)
            spectra = 1.0_real32
            call sigma2_state_write_particles(candidate%to_char(), 1, spectra, status, message)
            call require_ok(status, message)
            globals = 1.0_real32
            call sigma2_state_write_global(candidate%to_char(), globals, status, message)
            call require_ok(status, message)
            call sigma2_state_publish(candidate%to_char(), COMMITTED, status, message)
            call require_ok(status, message)
        end subroutine write_valid_store

    end subroutine test_retired_header_words

    !> sigma2_estimate_available on a project file: the committed state of
    !! the project's own layout is consumable; another grid is not. The layout
    !! digest covers the stack table, so the project's stk segment must be
    !! read (it was not: every state was reported "layout digest undefined")
    subroutine test_estimate_available_on_disk()
        use simple_sp_project,       only: sp_project
        use simple_ptcl_layout,      only: ptcl_layout_digest
        use simple_sigma2_bootstrap, only: sigma2_estimate_available
        integer, parameter :: NP = 6, BOXT = 8
        real,    parameter :: SMPDT = 1.5
        character(len=*), parameter :: PROJ = 'tmp_sigma2_state_tester_proj.simple'
        character(len=*), parameter :: COMM = 'tmp_sigma2_state_tester.bin'
        character(len=*), parameter :: COPY_DIR  = 'tmp_sigma2_state_tester_copy'
        character(len=*), parameter :: COPY_PROJ = COPY_DIR//'/test.simple'
        type(sp_project) :: project
        type(sigma2_state_header) :: header
        real(real32) :: spectra(BOXT/2, NP)
        logical :: active(NP)
        integer :: eo(NP), i, status
        integer(int64) :: digest
        write(*,'(A)') 'test_estimate_available_on_disk'
        call cleanup(COMM, trim(COMM)//'.next', PROJ, PROJ)
        ! as a project read from its file: named after the file, in the file's directory
        call project%update_projinfo(string(PROJ))
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
        call sigma2_state_init_header(header, 1, BOXT/2, NP, BOXT, SMPDT, &
            &1_int64, digest, SIGMA2_PROV_RESIDUAL)
        do i = 1, NP
            spectra(:,i) = real(i, real32)
        enddo
        active = .true.
        eo     = [(modulo(i-1,2), i = 1, NP)]
        call write_committed_store(COMM, header, spectra, active, eo)
        call project%set_sigma2_state_path(string(COMM))
        call project%write(string(PROJ))
        call assert_true(sigma2_estimate_available(string(PROJ), 'ptcl3D', BOXT, SMPDT), &
            &'a committed state of the project''s own layout is consumable from the project file')
        call assert_false(sigma2_estimate_available(string(PROJ), 'ptcl3D', BOXT+2, SMPDT), &
            &'negative control: the same state is not consumable on another grid')
        call simple_mkdir(COPY_DIR)
        call simple_copy_file(string(PROJ), string(COPY_PROJ))
        call assert_false(sigma2_estimate_available(string(COPY_PROJ), 'ptcl3D', BOXT, SMPDT), &
            &'a copied project cannot consume a relative state from its source directory')
        call del_file(COPY_PROJ)
        call simple_rmdir(COPY_DIR, status)
        call assert_int(0, status, 'copied project test directory is removed')
        call project%kill
        call cleanup(COMM, trim(COMM)//'.next', PROJ, PROJ)
    end subroutine test_estimate_available_on_disk

    subroutine test_particle_io_layout()
        type(sigma2_state_header) :: header
        real(real32), allocatable :: loaded(:,:), globals(:,:)
        real(real32) :: spectra(3,4)
        logical :: active(4)
        integer :: eo(4), status, funit
        integer(int64) :: range_words(16), state_words(32)
        integer :: first_row, last_row, kfrom, kto
        integer(int64), parameter :: DIGEST = 661_int64
        integer(int64), parameter :: RANGE_HEADER_BYTES = 256_int64
        integer(int64) :: generation, layout_digest, expected_bytes, actual_bytes
        character(len=128) :: message
        character(len=*), parameter :: CANDIDATE = 'layout_sigma2_state.next'
        character(len=*), parameter :: COMMITTED = 'layout_sigma2_state.bin'
        character(len=*), parameter :: RANGE = 'layout_sigma2_range.bin'
        write(*,'(A)') 'test_particle_io_layout'
        call cleanup(CANDIDATE, COMMITTED, RANGE, RANGE)
        spectra(:,1) = [1.0,2.0,3.0]
        spectra(:,2) = [2.0,3.0,4.0]
        spectra(:,3) = [3.0,4.0,5.0]
        spectra(:,4) = [4.0,5.0,6.0]
        active = .true.; eo = [0,1,0,1]
        call sigma2_state_init_header(header, 1, 3, 4, 8, 1.5, &
            &1_int64, DIGEST, SIGMA2_PROV_PSPEC)
        call sigma2_state_create_candidate(CANDIDATE, header, status, message)
        call require_ok(status, message)
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
        ! the retired range checksum word: zero when written, ignored when read
        open(newunit=funit, file=RANGE, access='stream', form='unformatted', status='old', action='readwrite')
        read(funit, pos=17) range_words
        call require(range_words(16) == 0_int64, 'new local range leaves the retired checksum word zero')
        range_words(16) = 12345_int64
        write(funit, pos=17) range_words
        close(funit)
        call sigma2_state_read_local_range(RANGE, generation, layout_digest, first_row, last_row, &
            &kfrom, kto, loaded, status, message)
        call require_ok(status, message)
        call require(all(loaded == spectra), 'local range accepts a historical nonzero checksum word')
        deallocate(loaded)
        open(newunit=funit, file=RANGE, access='stream', form='unformatted', status='old', action='readwrite')
        range_words(1) = 99_int64
        write(funit, pos=17) range_words
        close(funit)
        call sigma2_state_read_local_range(RANGE, generation, layout_digest, first_row, last_row, &
            &kfrom, kto, loaded, status, message)
        call require(status /= 0, 'local-range version corruption still fails structural validation')

        call sigma2_state_write_particles(CANDIDATE, 1, spectra, status, message)
        call require_ok(status, message)
        call sigma2_state_read_particles(CANDIDATE, 1, 4, loaded, status, message)
        call require_ok(status, message)
        call require(all(loaded == spectra), 'bulk particle write preserves every spectrum')
        deallocate(loaded)
        ! a state file is its header, its global spectra and its particle section
        call sigma2_state_read_header(CANDIDATE, header, status, message)
        call require_ok(status, message)
        inquire(file=CANDIDATE, size=actual_bytes, iostat=status)
        call require(status == 0 .and. actual_bytes == header%file_bytes, 'candidate size equals its recorded size')
        call require(header%file_bytes == header%particle_offset + 3_int64*4_int64*4_int64 - 1_int64, &
            &'the particle section is the last section of the state file')

        call sigma2_state_reduce(CANDIDATE, active, eo, status, message)
        call require_ok(status, message)
        call sigma2_state_validate_file(CANDIDATE, status, message)
        call require_ok(status, message)
        ! structural validation does not read the payload: an edited global spectrum still validates
        call sigma2_state_read_global(CANDIDATE, globals, status, message)
        call require_ok(status, message)
        open(newunit=funit, file=CANDIDATE, access='stream', form='unformatted', status='old', action='readwrite')
        write(funit, pos=header%global_offset) globals(1,1) + 10.0_real32
        close(funit)
        call sigma2_state_validate_file(CANDIDATE, status, message)
        call require_ok(status, message)
        open(newunit=funit, file=CANDIDATE, access='stream', form='unformatted', status='old', action='readwrite')
        write(funit, pos=header%global_offset) globals(1,1)
        close(funit)
        deallocate(globals)
        call sigma2_state_publish(CANDIDATE, COMMITTED, status, message)
        call require_ok(status, message)
        call sigma2_state_read_particles(COMMITTED, 1, 4, loaded, status, message)
        call require_ok(status, message)
        call require(all(loaded == spectra), 'committed particle spectra round-trip')
        deallocate(loaded)
        ! the retired state checksum words (18 and 32): zero when written, ignored when read
        open(newunit=funit, file=COMMITTED, access='stream', form='unformatted', status='old', action='read')
        read(funit, pos=17) state_words
        close(funit)
        call require(state_words(18) == 0_int64 .and. state_words(32) == 0_int64, &
            &'a new state leaves both retired checksum words zero')
        call patch_sigma2_header_word(COMMITTED, 32, 12345_int64)
        call patch_sigma2_header_word(COMMITTED, 18, 67890_int64)
        call sigma2_state_validate_file(COMMITTED, status, message)
        call require_ok(status, message)
        call sigma2_state_read_header(COMMITTED, header, status, message)
        call require_ok(status, message)
        call require(header%generation == 1_int64, 'historical checksum words are ignored on read')
        call patch_sigma2_header_word(COMMITTED, 1, 99_int64)
        call sigma2_state_validate_file(COMMITTED, status, message)
        call require(status /= 0, 'state version corruption still fails structural validation')
        call cleanup(CANDIDATE, COMMITTED, RANGE, RANGE)
    end subroutine test_particle_io_layout

    !> An uncovered row, a doubly covered row and a range of another transaction are refused before
    !! anything is published, the committed file untouched; a truncated candidate is refused
    subroutine test_recovery_guards()
        type(sigma2_state_header) :: header, committed_header
        type(string) :: ranges(2)
        real(real32) :: spectra(3,4)
        logical :: scheduled(4), active(4)
        integer :: eo(4), status, funit
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
        scheduled = .true.; active = .true.; eo = [0,1,0,1]
        call sigma2_state_init_header(header, 1, 3, 4, 8, 1.5, &
            &1_int64, DIGEST, SIGMA2_PROV_PSPEC)
        call write_committed_store(COMMITTED, header, spectra, active, eo)
        call sigma2_state_read_header(COMMITTED, committed_header, status, message)
        call require_ok(status, message)
        call read_file_bytes(COMMITTED, committed_bytes)

        ! rows 3:4 scheduled but not covered
        call sigma2_state_write_local_range(ranges(1)%to_char(), 2_int64, DIGEST, 1, spectra(:,:2), &
            &1, 3, status, message)
        call require_ok(status, message)
        call sigma2_state_assemble_candidate(COMMITTED, CANDIDATE, ranges(:1), scheduled, active, eo, status, message)
        call require(status /= 0, 'an uncovered scheduled row is refused')
        call assert_committed_generation(COMMITTED, committed_header%generation)
        call assert_file_unchanged(COMMITTED, committed_bytes)

        ! row 3 covered twice
        call del_file(CANDIDATE); call del_file(ranges(1))
        call sigma2_state_write_local_range(ranges(1)%to_char(), 2_int64, DIGEST, 1, spectra(:,:3), &
            &1, 3, status, message)
        call require_ok(status, message)
        call sigma2_state_write_local_range(ranges(2)%to_char(), 2_int64, DIGEST, 3, spectra(:,3:), &
            &1, 3, status, message)
        call require_ok(status, message)
        call sigma2_state_assemble_candidate(COMMITTED, CANDIDATE, ranges, scheduled, active, eo, status, message)
        call require(status /= 0, 'a doubly covered row is refused')
        call assert_committed_generation(COMMITTED, committed_header%generation)
        call assert_file_unchanged(COMMITTED, committed_bytes)

        ! a complete range of another transaction
        call del_file(CANDIDATE); call del_file(ranges(1)); call del_file(ranges(2))
        call sigma2_state_write_local_range(ranges(1)%to_char(), 3_int64, DIGEST, 1, spectra, &
            &1, 3, status, message)
        call require_ok(status, message)
        call sigma2_state_assemble_candidate(COMMITTED, CANDIDATE, ranges(:1), scheduled, active, eo, status, message)
        call require(status /= 0, 'a range of another generation is refused')
        call assert_file_unchanged(COMMITTED, committed_bytes)

        call del_file(CANDIDATE)
        open(newunit=funit, file=CANDIDATE, access='stream', form='unformatted', status='replace')
        write(funit) 1
        close(funit)
        call sigma2_state_validate_file(CANDIDATE, status, message)
        call require(status /= 0, 'truncated candidate rejected')
        call assert_committed_generation(COMMITTED, committed_header%generation)
        call assert_file_unchanged(COMMITTED, committed_bytes)
        call cleanup(CANDIDATE, COMMITTED, ranges(1)%to_char(), ranges(2)%to_char())
    end subroutine test_recovery_guards

    !> A publication that fails (the assembled candidate's header damaged before it is published) leaves the
    !! committed store and the range files in place, and the commit reruns from the same ranges
    subroutine test_failed_publication()
        integer(int64),   parameter :: DIGEST = 881_int64
        character(len=*), parameter :: CANDIDATE = 'failpub_sigma2_state.next'
        character(len=*), parameter :: COMMITTED = 'failpub_sigma2_state.bin'
        character(len=*), parameter :: RANGE     = 'failpub_sigma2_range.bin'
        type(sigma2_state_header) :: header
        real(real32) :: spectra(3,4)
        logical :: scheduled(4), active(4)
        integer :: eo(4), status
        integer(int8), allocatable :: committed_bytes(:)
        character(len=128) :: message
        write(*,'(A)') 'test_failed_publication'
        call cleanup(CANDIDATE, COMMITTED, RANGE, RANGE)
        spectra = 1.0_real32
        scheduled = .true.; active = .true.; eo = [0,1,0,1]
        call sigma2_state_init_header(header, 1, 3, 4, 8, 1.5, 1_int64, DIGEST, SIGMA2_PROV_PSPEC)
        call write_committed_store(COMMITTED, header, spectra, active, eo)
        call read_file_bytes(COMMITTED, committed_bytes)
        spectra = 2.0_real32
        call sigma2_state_write_local_range(RANGE, 2_int64, DIGEST, 1, spectra, 1, 3, status, message)
        call require_ok(status, message)
        call sigma2_state_assemble_candidate(COMMITTED, CANDIDATE, [string(RANGE)], scheduled, active, eo, &
            &status, message)
        call require_ok(status, message)
        ! a candidate whose header no longer validates (an unsupported version)
        call patch_sigma2_header_word(CANDIDATE, 1, 99_int64)
        call sigma2_state_publish(CANDIDATE, COMMITTED, status, message)
        call require(status /= 0, 'a damaged candidate is not published')
        call assert_file_unchanged(COMMITTED, committed_bytes)
        call require(file_exists(RANGE), 'a failed publication leaves the range files')
        ! the rerun
        call sigma2_state_assemble_candidate(COMMITTED, CANDIDATE, [string(RANGE)], scheduled, active, eo, &
            &status, message)
        call require_ok(status, message)
        call sigma2_state_publish(CANDIDATE, COMMITTED, status, message)
        call require_ok(status, message)
        call assert_committed_generation(COMMITTED, 2_int64)
        call cleanup(CANDIDATE, COMMITTED, RANGE, RANGE)
    end subroutine test_failed_publication

    subroutine test_update_preparation()
        type(sigma2_state_header) :: header
        type(string) :: candidate_path, range_path
        real(real32) :: spectra(3,4)
        logical :: active(4)
        integer :: eo(4), status, funit
        integer(int64), parameter :: DIGEST = 991_int64
        integer(int64) :: next_gen
        character(len=128) :: message
        character(len=*), parameter :: COMMITTED = 'sigma2_state.bin'
        write(*,'(A)') 'test_update_preparation'
        candidate_path = sigma2_state_candidate_path(COMMITTED, 1_int64)
        range_path = sigma2_state_range_path(COMMITTED, 1_int64, 2, 3)
        call require(index(candidate_path%to_char(), 'sigma2_state.g1.next') > 0, &
            &'candidate path is scoped to the generation it commits')
        call require(index(range_path%to_char(), 'sigma2_state.g1.part002.range') > 0, &
            &'canonical range path is scoped to the generation and includes the padded partition')
        call del_file(COMMITTED)
        call del_file(trim(COMMITTED)//'.next')
        spectra(:,1) = [1.0,2.0,3.0]
        spectra(:,2) = [2.0,3.0,4.0]
        spectra(:,3) = [3.0,4.0,5.0]
        spectra(:,4) = [4.0,5.0,6.0]
        active = .true.; eo = [0,1,0,1]
        call sigma2_state_init_header(header, 1, 3, 4, 8, 1.5, &
            &1_int64, DIGEST, SIGMA2_PROV_PSPEC)
        call write_committed_store(COMMITTED, header, spectra, active, eo)
        ! the next transaction is named for the generation it will commit
        call sigma2_state_next_generation(COMMITTED, next_gen, status, message)
        call require_ok(status, message)
        call require(next_gen == 2_int64, 'next generation follows the committed one')
        candidate_path = sigma2_state_candidate_path(COMMITTED, next_gen)
        range_path     = sigma2_state_range_path(COMMITTED, next_gen, 2, 3)
        call require(index(candidate_path%to_char(), 'sigma2_state.g2.next') > 0, &
            &'next candidate path is scoped to the next generation')
        ! what a crashed transaction of that generation left
        open(newunit=funit, file=candidate_path%to_char(), status='replace')
        close(funit)
        open(newunit=funit, file=range_path%to_char(), status='replace')
        close(funit)
        call sigma2_state_prepare_update(COMMITTED, 3, 3, status, message)
        call require_ok(status, message)
        call require(.not. file_exists(candidate_path%to_char()), 'the prepared update removes a stale candidate')
        call require(.not. file_exists(range_path%to_char()), 'the prepared update removes a stale range')
        call assert_committed_generation(COMMITTED, 1_int64)
        ! a store that is not committed is no update source
        call patch_sigma2_header_word(COMMITTED, 12, 1_int64)
        call sigma2_state_prepare_update(COMMITTED, 3, 3, status, message)
        call require(status /= 0, 'an update of a store that is not committed is refused')
        call del_file(COMMITTED)
    end subroutine test_update_preparation

    !> Two committed stores concatenated in blocks into one candidate: the rows follow one another and the
    !! global spectra equal a full reduction of the result
    subroutine test_concatenation()
        integer,          parameter :: NA = 5000, NB = 3
        character(len=*), parameter :: STORE_A = 'concat_a_sigma2_state.bin'
        character(len=*), parameter :: STORE_B = 'concat_b_sigma2_state.bin'
        character(len=*), parameter :: CANDIDATE = 'concat_sigma2_state.next'
        type(sigma2_state_header) :: header
        real(real32), allocatable :: loaded(:,:)
        real(real32) :: rows_a(3,NA), rows_b(3,NB)
        logical :: active(NA+NB)
        integer :: eo(NA+NB), i, status
        character(len=128) :: message
        write(*,'(A)') 'test_concatenation'
        call cleanup(STORE_A, STORE_B, CANDIDATE, CANDIDATE)
        do i = 1, NA
            rows_a(:,i) = real([i, i+1, i+2], real32)
        enddo
        do i = 1, NB
            rows_b(:,i) = real([7*i, 7*i+1, 7*i+2], real32)
        enddo
        active = .true.
        active(2) = .false.
        eo = [(modulo(i-1,2), i = 1, NA+NB)]
        call sigma2_state_init_header(header, 1, 3, NA, 8, 1.5, 4_int64, 111_int64, SIGMA2_PROV_RESIDUAL)
        call write_committed_store(STORE_A, header, rows_a, active(:NA), eo(:NA))
        call sigma2_state_init_header(header, 1, 3, NB, 8, 1.5, 2_int64, 222_int64, SIGMA2_PROV_RESIDUAL)
        call write_committed_store(STORE_B, header, rows_b, active(NA+1:), eo(NA+1:))
        call sigma2_state_init_header(header, 1, 3, NA+NB, 8, 1.5, 5_int64, 333_int64, SIGMA2_PROV_RESIDUAL)
        call sigma2_state_concatenate([string(STORE_A), string(STORE_B)], CANDIDATE, header, active, eo, &
            &status, message)
        call require_ok(status, message)
        call sigma2_state_read_particles(CANDIDATE, 1, NA+NB, loaded, status, message)
        call require_ok(status, message)
        call require(all(loaded(:,:NA) == rows_a) .and. all(loaded(:,NA+1:) == rows_b), &
            &'the concatenated rows follow one another')
        call require_full_reduction(CANDIDATE, active, eo)
        ! more rows than the candidate holds
        call sigma2_state_init_header(header, 1, 3, NA, 8, 1.5, 5_int64, 333_int64, SIGMA2_PROV_RESIDUAL)
        call sigma2_state_concatenate([string(STORE_A), string(STORE_B)], CANDIDATE, header, active(:NA), eo(:NA), &
            &status, message)
        call require(status /= 0, 'sources with more rows than the candidate are refused')
        call cleanup(STORE_A, STORE_B, CANDIDATE, CANDIDATE)
    end subroutine test_concatenation

    !> an active record with a non-positive shell is skipped by the streaming commit (with a warning)
    !! instead of aborting; the global spectra are the means of the valid records, as a full reduction
    !! forms them
    subroutine test_invalid_record_skip()
        type(sigma2_state_header) :: header
        type(string) :: range_path
        real(real32), allocatable :: stored(:,:)
        real(real32) :: spectra(3,4), expected_even(3), expected_odd(3)
        logical :: active(4)
        integer :: eo(4), status
        integer(int64), parameter :: DIGEST = 1231_int64
        character(len=128) :: message
        character(len=*), parameter :: COMMITTED = 'skip_sigma2_state.bin'
        character(len=*), parameter :: CANDIDATE = 'skip_sigma2_state.next'
        write(*,'(A)') 'test_invalid_record_skip'
        range_path = sigma2_state_range_path(COMMITTED, 2_int64, 1, 1)
        call cleanup(COMMITTED, CANDIDATE, range_path%to_char(), range_path%to_char())
        active = .true.; eo = [0,1,0,1]
        spectra = 1.0_real32
        call sigma2_state_init_header(header, 1, 3, 4, 8, 1.5, &
            &1_int64, DIGEST, SIGMA2_PROV_PSPEC)
        call write_committed_store(COMMITTED, header, spectra, active, eo)
        spectra(:,1) = [1.0,2.0,3.0]
        spectra(:,2) = [2.0,3.0,4.0]
        spectra(:,3) = [3.0,0.0,5.0] ! invalid: a non-positive shell
        spectra(:,4) = [4.0,5.0,6.0]
        expected_even = spectra(:,1)
        expected_odd  = 0.5*(spectra(:,2)+spectra(:,4))
        call sigma2_state_write_local_range(range_path%to_char(), 2_int64, DIGEST, 1, spectra, &
            &1, 3, status, message)
        call require_ok(status, message)
        call sigma2_state_assemble_candidate(COMMITTED, CANDIDATE, [range_path], [.true.,.true.,.true.,.true.], &
            &active, eo, status, message)
        call require(status == 0, 'invalid record is skipped, not fatal')
        call sigma2_state_read_global(CANDIDATE, stored, status, message)
        call require_ok(status, message)
        call require(all(abs(stored(:,1)-expected_even) <= 1.e-6*expected_even), &
            &'even global mean excludes the invalid record')
        call require(all(abs(stored(:,2)-expected_odd) <= 1.e-6*expected_odd), &
            &'odd global mean is the mean of the valid records')
        call require_full_reduction(CANDIDATE, active, eo)
        call sigma2_state_publish(CANDIDATE, COMMITTED, status, message)
        call require_ok(status, message)
        call cleanup(COMMITTED, CANDIDATE, range_path%to_char(), range_path%to_char())
    end subroutine test_invalid_record_skip

    !> a committed store at path of identity header and records spectra, its global spectra a full
    !! reduction (the fixture a transaction starts from)
    subroutine write_committed_store(path, header, spectra, active, eo)
        character(len=*),          intent(in)    :: path
        type(sigma2_state_header), intent(inout) :: header
        real(real32),              intent(in)    :: spectra(:,:)
        logical,                   intent(in)    :: active(:)
        integer,                   intent(in)    :: eo(:)
        character(len=128) :: message
        integer :: status
        call del_file(path)
        call sigma2_state_create_candidate(trim(path)//'.next', header, status, message)
        call require_ok(status, message)
        call sigma2_state_write_particles(trim(path)//'.next', 1, spectra, status, message)
        call require_ok(status, message)
        call sigma2_state_reduce(trim(path)//'.next', active, eo, status, message)
        call require_ok(status, message)
        call sigma2_state_publish(trim(path)//'.next', path, status, message)
        call require_ok(status, message)
    end subroutine write_committed_store

    !> The candidate's global spectra equal a full reduction of its records: the reduction re-reads every
    !! record and rewrites the global section, which must not move beyond rounding (the check the commit
    !! made by re-reading the store before the streaming commit)
    subroutine require_full_reduction(candidate, active, eo)
        character(len=*), intent(in) :: candidate
        logical,          intent(in) :: active(:)
        integer,          intent(in) :: eo(:)
        real(real32), allocatable :: streamed(:,:), reduced(:,:)
        character(len=128) :: message
        integer :: status
        call sigma2_state_read_global(candidate, streamed, status, message)
        call require_ok(status, message)
        call sigma2_state_reduce(candidate, active, eo, status, message)
        call require_ok(status, message)
        call sigma2_state_read_global(candidate, reduced, status, message)
        call require_ok(status, message)
        call require(all(abs(real(streamed,real64)-real(reduced,real64)) <= &
            &8.0_real64*epsilon(1.0_real32)*max(1.0_real64,abs(real(reduced,real64)))), &
            &'the streamed global spectra equal a full reduction of the records')
    end subroutine require_full_reduction

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

    !> overwrite header word iword of the store fname (the header carries no checksum)
    subroutine patch_sigma2_header_word( fname, iword, value )
        character(len=*), intent(in) :: fname
        integer,        intent(in) :: iword
        integer(int64), intent(in) :: value
        integer(int64) :: words(32)
        integer :: funit, io_stat
        open(newunit=funit, file=fname, access='stream', form='unformatted', status='old', &
            &action='readwrite', iostat=io_stat)
        call require(io_stat == 0, 'store can be opened for patching')
        read(funit, pos=17, iostat=io_stat) words
        words(iword) = value
        write(funit, pos=17, iostat=io_stat) words
        close(funit)
        call require(io_stat == 0, 'store header can be patched')
    end subroutine patch_sigma2_header_word

    !> A committed unit store of nptcls records over the shells 1..box/2 (fdim(box)-1): the fixture
    !! testers read their residual records from (read_part), as a matcher reads them
    subroutine write_test_sigma2_store( fname, nptcls, box )
        character(len=*), intent(in) :: fname
        integer,          intent(in) :: nptcls, box
        type(sigma2_state_header) :: header
        real(real32), allocatable :: records(:,:)
        real(real32) :: globals(box/2,2)
        character(len=128) :: candidate, message
        integer :: status
        candidate = trim(fname)//'.next'
        call del_file(fname)
        call del_file(candidate)
        call sigma2_state_init_header(header, 1, box/2, nptcls, box, 1.0, 1_int64, 4242_int64, SIGMA2_PROV_PSPEC)
        call sigma2_state_create_candidate(candidate, header, status, message)
        call require_ok(status, message)
        allocate(records(box/2,nptcls), source=1.0_real32)
        call sigma2_state_write_particles(candidate, 1, records, status, message)
        call require_ok(status, message)
        globals = 1.0_real32
        call sigma2_state_write_global(candidate, globals, status, message)
        call require_ok(status, message)
        call sigma2_state_publish(candidate, fname, status, message)
        call require_ok(status, message)
    end subroutine write_test_sigma2_store

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
