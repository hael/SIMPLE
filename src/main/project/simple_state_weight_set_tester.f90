!@descr: unit tests of state_weight_set: round trip, refusals, a simulated crash, a work project's derived set, withdrawal and projects of earlier releases
module simple_state_weight_set_tester
use, intrinsic :: iso_fortran_env, only: int8, int64, real32, real64
use simple_test_utils
use simple_oris,               only: oris
use simple_sp_project,         only: sp_project
use simple_string,             only: string
use simple_syslib,             only: simple_getcwd, simple_chdir, simple_mkdir, file_exists
use simple_state_weight_set,   only: state_weight_set, discard_state_weight_set, STATE_WEIGHTS_KIND_PARTITION, &
    &STATE_WEIGHTS_KIND_KERNEL
use simple_state_weights_file, only: state_weights_file_header, state_weights_fname, state_weights_write_state_file, &
    &STATE_WEIGHTS_MANIFEST_FNAME
implicit none
private
public :: run_all_state_weight_set_tests

integer,          parameter :: NP      = 12
integer,          parameter :: NSTATES = 3
character(len=*), parameter :: PROJ    = 'state_weight_set_test.simple'

contains

    subroutine run_all_state_weight_set_tests()
        type(string) :: cwd_saved, root
        integer      :: nfail0
        write(*,'(A)') '**** running all state weight set tests ****'
        nfail0 = tests_failed
        call enter_fixture('state_weight_set', cwd_saved, root)
        call test_round_trip_partition()
        call test_round_trip_kernel()
        call test_refusals()
        call test_crash_between_files_and_manifest()
        call test_work_state()
        call test_discard()
        call test_earlier_release_entries()
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine run_all_state_weight_set_tests

    !> rows 1..NP-2 selected; PARTITION responsibilities, label the argmax
    subroutine partition_table( pinds, weights, labels )
        integer, allocatable, intent(out) :: pinds(:), labels(:)
        real,    allocatable, intent(out) :: weights(:,:)
        integer :: i, nsel
        nsel = NP - 2
        allocate(pinds(nsel), weights(nsel,NSTATES), labels(nsel))
        do i = 1, nsel
            pinds(i)     = i
            weights(i,:) = [real(i), real(2*mod(i,3)+1), 3.5]
            weights(i,:) = weights(i,:) / sum(weights(i,:))
            labels(i)    = maxloc(weights(i,:), dim=1)
        enddo
    end subroutine partition_table

    subroutine test_round_trip_partition()
        class(sp_project), allocatable :: project
        type(state_weight_set)    :: wset, wread
        integer,      allocatable :: pinds(:), labels(:), labels_read(:)
        real,         allocatable :: weights(:,:)
        real(real32), allocatable :: column(:)
        logical,      allocatable :: selected(:)
        real(real64) :: mass
        integer      :: s, i, status
        logical      :: l_exact
        character(len=256) :: message
        write(*,'(A)') 'test_round_trip_partition'
        allocate(project)
        call make_project(project, 'swset_particles')
        call partition_table(pinds, weights, labels)
        call wset%publish(project, project%os_ptcl3D, string(PROJ), pinds, weights, labels, 'tester')
        call assert_int(STATE_WEIGHTS_KIND_PARTITION, wset%get_kind(), 'responsibilities publish as PARTITION')
        call assert_true(wset%get_generation() == 1_int64, 'the first publication is generation 1')
        call wset%kill
        ! a fresh read through the project file, as a consumer sees it
        call project%kill
        call project%read(string(PROJ))
        call wread%new(project, project%os_ptcl3D, status, message)
        call assert_int(0, status, 'the published set validates: '//trim(message))
        call assert_int(NSTATES, wread%get_nstates(), 'state count round trip')
        call assert_int(NP, wread%get_nptcls(), 'particle count round trip')
        call assert_char('tester', wread%get_producer(), 'producer round trip')
        l_exact = .true.
        do s = 1, NSTATES
            call wread%get_weights(s, column)
            do i = 1, NP
                if( i <= size(pinds) )then
                    l_exact = l_exact .and. column(i) == real(weights(i,s), real32)
                else
                    l_exact = l_exact .and. column(i) == 0._real32
                endif
            enddo
            mass = sum(real(real(weights(:,s), real32), real64))
            call assert_true(abs(wread%get_mass(s) - mass) <= 1.e-12_real64 * mass, 'manifest mass equals the rows')
            call assert_int(count(labels == s), wread%get_pop(s), 'manifest hard population equals the labels')
            call assert_true(abs(wread%applied_mass(s, 0., pinds) - mass) <= 1.e-12_real64 * mass, &
                &'applied mass over the selection equals the manifest mass')
        enddo
        call assert_true(l_exact, 'weights round trip bit for bit, zero outside the selection')
        call wread%get_labels(labels_read)
        call assert_true(all(labels_read(1:size(pinds)) == labels) .and. all(labels_read(size(pinds)+1:) == 0), &
            &'hard labels round trip')
        call wread%get_selection(0., selected)
        call assert_true(all(selected(1:size(pinds))) .and. .not. any(selected(size(pinds)+1:)), &
            &'the selection is the rows with positive weight')
        call wread%get_applied_weights(1, 0.2, column)
        call assert_true(all(column <= 0. .or. column > 0.2), 'applied weights are zero at or below the threshold')
        call wread%kill
        call project%kill
        deallocate(project)
    end subroutine test_round_trip_partition

    subroutine test_round_trip_kernel()
        class(sp_project), allocatable :: project
        type(state_weight_set) :: wset
        integer, allocatable   :: pinds(:), labels(:)
        real,    allocatable   :: weights(:,:)
        integer :: status
        character(len=256) :: message
        write(*,'(A)') 'test_round_trip_kernel'
        allocate(project)
        call make_project(project, 'swset_particles')
        call partition_table(pinds, weights, labels)
        weights = min(1., 1.5 * weights)
        call wset%publish(project, project%os_ptcl3D, string(PROJ), pinds, weights, labels, 'tester')
        call assert_int(STATE_WEIGHTS_KIND_KERNEL, wset%get_kind(), 'rows that do not sum to one publish as KERNEL')
        call wset%kill
        call wset%new(project, project%os_ptcl3D, status, message)
        call assert_int(0, status, 'the KERNEL set validates: '//trim(message))
        call wset%kill
        call project%kill
        deallocate(project)
    end subroutine test_round_trip_kernel

    subroutine test_refusals()
        class(sp_project), allocatable :: project, other
        type(state_weight_set) :: wset
        type(string)           :: manifest, fname
        integer, allocatable   :: pinds(:), labels(:)
        real,    allocatable   :: weights(:,:)
        integer(int8) :: byte
        integer       :: status, funit
        logical       :: found
        character(len=256) :: message
        write(*,'(A)') 'test_refusals'
        allocate(project, other)
        call make_project(project, 'swset_particles')
        call partition_table(pinds, weights, labels)
        call wset%publish(project, project%os_ptcl3D, string(PROJ), pinds, weights, labels, 'tester')
        fname = state_weights_fname(wset%get_generation(), 2)
        call wset%kill
        call project%get_state_weights(manifest, found)
        ! another particle layout (another stack) and another particle count
        call make_project(other, 'other_particles')
        call other%add_state_weights2os_out(manifest)
        call wset%new(other, other%os_ptcl3D, status, message)
        call assert_true(status /= 0, 'a set of another particle layout is refused')
        call other%kill
        call make_project(other, 'swset_particles', NP + 1)
        call other%add_state_weights2os_out(manifest)
        call wset%new(other, other%os_ptcl3D, status, message)
        call assert_true(status /= 0, 'a set of another particle count is refused')
        ! a byte flipped in place: same size, checksum mismatch
        open(newunit=funit, file=fname%to_char(), access='stream', form='unformatted', status='old', action='readwrite')
        read(funit, pos=101) byte
        write(funit, pos=101) not(byte)
        close(funit)
        call wset%new(project, project%os_ptcl3D, status, message)
        call assert_true(status /= 0, 'a state file with a wrong checksum is refused')
        ! restore, then append a byte: size mismatch
        open(newunit=funit, file=fname%to_char(), access='stream', form='unformatted', status='old', action='readwrite')
        write(funit, pos=101) byte
        close(funit)
        call wset%new(project, project%os_ptcl3D, status, message)
        call assert_int(0, status, 'the restored file validates again: '//trim(message))
        call wset%kill
        open(newunit=funit, file=fname%to_char(), access='stream', form='unformatted', status='old', &
            &action='write', position='append')
        write(funit) byte
        close(funit)
        call wset%new(project, project%os_ptcl3D, status, message)
        call assert_true(status /= 0, 'a state file of the wrong size is refused')
        call wset%kill
        call other%kill
        call project%kill
        deallocate(project, other)
    end subroutine test_refusals

    !> files of a new generation published without their manifest leave the previous generation selected
    subroutine test_crash_between_files_and_manifest()
        class(sp_project), allocatable :: project
        type(state_weight_set)          :: wset
        type(state_weights_file_header) :: header
        integer,       allocatable :: pinds(:), labels(:)
        real,          allocatable :: weights(:,:)
        real(real32),  allocatable :: column(:)
        integer(int8), allocatable :: flags(:)
        integer(int64) :: gen
        integer        :: status, s
        character(len=256) :: message
        write(*,'(A)') 'test_crash_between_files_and_manifest'
        allocate(project)
        call make_project(project, 'swset_particles')
        call partition_table(pinds, weights, labels)
        call wset%publish(project, project%os_ptcl3D, string(PROJ), pinds, weights, labels, 'tester')
        gen = wset%get_generation()
        call wset%kill
        ! the crashed publication: the next generation's files, other weights, no manifest
        allocate(column(NP), source=0.25_real32)
        allocate(flags(NP),  source=0_int8)
        do s = 1, NSTATES
            header%nptcls        = NP
            header%nstates       = NSTATES
            header%state         = s
            header%generation    = gen + 1_int64
            header%layout_digest = 7_int64
            header%kind          = STATE_WEIGHTS_KIND_KERNEL
            call state_weights_write_state_file(state_weights_fname(gen + 1_int64, s), header, column, flags, status, message)
        enddo
        call wset%new(project, project%os_ptcl3D, status, message)
        call assert_int(0, status, 'after the simulated crash the set still validates: '//trim(message))
        call assert_true(wset%get_generation() == gen, 'the previous generation stays selected')
        call wset%get_weights(1, column)
        call assert_true(column(1) == real(weights(1,1), real32), 'and its weights are read')
        call wset%kill
        ! a complete publication afterwards supersedes it and removes the older generation's files
        call wset%publish(project, project%os_ptcl3D, string(PROJ), pinds, weights, labels, 'tester')
        call assert_true(wset%get_generation() == gen + 1_int64, 'the next publication is the next generation')
        call assert_false(file_exists(state_weights_fname(gen, 1)), 'the superseded generation is deleted')
        call wset%kill
        call project%kill
        deallocate(project)
    end subroutine test_crash_between_files_and_manifest

    !> a project whose out segment carries flex_weights entries of an earlier release loads; the entries are ignored
    !> a withdrawn set leaves the project (in memory and on disk) and its manifest and state files are deleted
    subroutine test_discard()
        class(sp_project), allocatable :: project
        type(state_weight_set) :: wset
        type(string)           :: manifest
        integer, allocatable   :: pinds(:), labels(:)
        real,    allocatable   :: weights(:,:)
        integer(int64)         :: generation
        logical                :: found, l_files_gone
        integer                :: s
        write(*,'(A)') 'test_discard'
        allocate(project)
        call make_project(project, 'discard_particles')
        call partition_table(pinds, weights, labels)
        call wset%publish(project, project%os_ptcl3D, string(PROJ), pinds, weights, labels, 'tester')
        generation = wset%get_generation()
        call wset%kill
        call discard_state_weight_set(project, string(PROJ))
        call project%get_state_weights(manifest, found)
        call assert_false(found, 'a withdrawn set is no longer registered in memory')
        call project%kill
        call project%read(string(PROJ))
        call project%get_state_weights(manifest, found)
        call assert_false(found, 'a withdrawn set is no longer registered in the project file')
        l_files_gone = .not. file_exists(STATE_WEIGHTS_MANIFEST_FNAME)
        do s = 1, NSTATES
            if( file_exists(state_weights_fname(generation, s)) ) l_files_gone = .false.
        enddo
        call assert_true(l_files_gone, 'the manifest and the state files of a withdrawn set are deleted')
        call discard_state_weight_set(project, string(PROJ))
        call assert_false(file_exists(STATE_WEIGHTS_MANIFEST_FNAME), 'withdrawing from a project without a set is a no-op')
        call project%kill
        deallocate(project)
    end subroutine test_discard

    subroutine test_earlier_release_entries()
        class(sp_project), allocatable :: project
        type(string) :: manifest
        logical      :: found
        integer      :: n
        write(*,'(A)') 'test_earlier_release_entries'
        allocate(project)
        call make_project(project, 'old_particles')
        call project%os_out%new(1, is_ptcl=.false.)
        call project%os_out%set(1, 'flex_weights', '/somewhere/flex_weights_state_001.bin')
        call project%os_out%set(1, 'imgkind',      'flex_weights')
        call project%os_out%set(1, 'state',        1)
        call project%write(string('old_release.simple'))
        call project%kill
        call project%read(string('old_release.simple'))
        n = project%os_out%get_noris()
        call assert_int(1, n, 'the out segment of the earlier release reads')
        call project%get_state_weights(manifest, found)
        call assert_false(found, 'its flex_weights entry is not a state weight set')
        call project%kill
        deallocate(project)
    end subroutine test_earlier_release_entries

    !> One state of a PARTITION parent published for a work project of some of its rows (in its own
    !! directory, as refine3D_auto state=X does): one column, rows remapped, kind and identity kept
    subroutine test_work_state()
        integer, parameter :: STATE = 2
        class(sp_project), allocatable :: parent, work
        type(state_weight_set)    :: pset, wset, wread
        type(string)              :: root
        integer,      allocatable :: pinds(:), labels(:), labels_read(:), parent_rows(:)
        real,         allocatable :: weights(:,:)
        real(real32), allocatable :: column(:)
        integer(int64) :: pgen, pdigest
        integer        :: i, status, nwork
        character(len=256) :: message
        write(*,'(A)') 'test_work_state'
        allocate(parent, work)
        call make_project(parent, 'swset_particles')
        call partition_table(pinds, weights, labels)
        call pset%publish(parent, parent%os_ptcl3D, string(PROJ), pinds, weights, labels, 'tester')
        pgen    = pset%get_generation()
        pdigest = pset%get_layout_digest()
        parent_rows = [2, 4, 5, 7, 9, 10]
        nwork       = size(parent_rows)
        call simple_getcwd(root)
        call simple_mkdir('work_state')
        call simple_chdir('work_state')
        call make_project(work, 'swset_particles', n=nwork)
        do i = 1, nwork
            call work%os_ptcl3D%set(i, 'indstk', parent_rows(i))
        enddo
        call wset%publish_work_state(pset, STATE, parent_rows, work, work%os_ptcl3D, string(PROJ))
        call wset%kill
        call wread%new(work, work%os_ptcl3D, status, message)
        call assert_int(0, status, 'the derived single-state set validates: '//trim(message))
        call assert_int(1, wread%get_nstates(), 'the derived set holds one state')
        call assert_int(STATE_WEIGHTS_KIND_PARTITION, wread%get_kind(), 'the derived set keeps its parent''s kind')
        call assert_true(wread%get_parent_generation() == pgen .and. wread%get_parent_layout_digest() == pdigest, &
            &'the derived set records its parent''s identity')
        call assert_int(STATE, wread%get_parent_state(), 'the derived set records its parent''s state')
        call wread%get_weights(1, column)
        call assert_true(all(column == real(weights(parent_rows,STATE), real32)), &
            &'the derived weights are the parent''s column at the mapped rows, bit for bit')
        call wread%get_labels(labels_read)
        call assert_true(all(labels_read == merge(1, 0, labels(parent_rows) == STATE)), &
            &'the derived hard labels are the parent''s labels of that state')
        call assert_int(count(labels(parent_rows) == STATE), wread%get_pop(1), 'the derived hard population')
        call wread%kill
        call simple_chdir(root)
        call pset%kill
        call pset%new(parent, parent%os_ptcl3D, status, message)
        call assert_int(0, status, 'the parent set is untouched: '//trim(message))
        call assert_true(pset%get_generation() == pgen, 'the parent keeps its generation')
        call pset%kill
        call parent%kill
        call work%kill
        deallocate(parent, work)
    end subroutine test_work_state

    !> an in-memory project of n particle rows over one stack, written to PROJ (whose name is the lineage)
    subroutine make_project( project, stkname, n )
        class(sp_project), intent(inout) :: project
        character(len=*),  intent(in)    :: stkname
        integer, optional, intent(in)    :: n
        type(string) :: cwd
        integer :: i, np_here
        np_here = NP
        if( present(n) ) np_here = n
        call project%kill
        call project%projinfo%new(1, is_ptcl=.false.)
        call project%projinfo%set(1, 'projname', PROJ(1:len(PROJ)-7))
        call project%projinfo%set(1, 'projfile', PROJ)
        call simple_getcwd(cwd)
        call project%projinfo%set(1, 'cwd', cwd%to_char())
        call project%os_stk%new(1, is_ptcl=.false.)
        call project%os_stk%set(1, 'stk',        '/data/'//stkname//'.mrcs')
        call project%os_stk%set(1, 'fromp',      1)
        call project%os_stk%set(1, 'top',        np_here)
        call project%os_stk%set(1, 'nptcls_stk', np_here)
        call project%os_stk%set(1, 'box',        16)
        call project%os_stk%set(1, 'smpd',       1.5)
        call project%os_ptcl3D%new(np_here, is_ptcl=.true.)
        do i = 1, np_here
            call project%os_ptcl3D%set(i, 'stkind', 1)
            call project%os_ptcl3D%set(i, 'indstk', i)
            call project%os_ptcl3D%set_state(i, 1)
        enddo
        project%os_ptcl2D = project%os_ptcl3D
        call project%write(string(PROJ))
    end subroutine make_project

end module simple_state_weight_set_tester
