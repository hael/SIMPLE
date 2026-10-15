!@descr: unit tests of the project-owned even/odd partition: micrograph groups, routes, appended stacks, earlier releases, refusals, merge and the isolated random stream
module simple_eo_partition_tester
use, intrinsic :: iso_fortran_env, only: int64
use simple_test_utils
use simple_defs,           only: l_distr_worker_glob, STDLEN
use simple_oris,           only: oris
use simple_sp_project,     only: sp_project
use simple_string,         only: string
use simple_cmdline,        only: cmdline
use simple_rnd,            only: ran3
use simple_forked_process, only: forked_process, FORK_STATUS_FAILED
use simple_projfile_utils, only: merge_selected_project_files
use simple_eo_partition,   only: ensure_field_halves, eo_transaction, eo_partition_status, EO_COMPLETE, &
    &EO_MISSING, EO_PARTIAL, EO_INCONSISTENT
implicit none
private
public :: run_all_eo_partition_tests

integer,          parameter :: REFUSE_WORKER_MISSING = 1, REFUSE_MASTER_INCONSISTENT = 2
character(len=*), parameter :: EO_PROJ = 'eo_partition_test.simple'

!> a process that must stop: the project's halves are refused
type, extends(forked_process) :: refusing_partition
    integer :: what = REFUSE_WORKER_MISSING
contains
    procedure :: execute => refuse_partition
end type refusing_partition

contains

    subroutine run_all_eo_partition_tests()
        type(string) :: cwd_saved, root
        integer      :: nfail0
        write(*,'(A)') '**** running all eo partition tests ****'
        nfail0 = tests_failed
        call enter_fixture('eo_partition', cwd_saved, root)
        call test_micrograph_groups()
        call test_no_split_by_size()
        call test_too_few_groups_status()
        call test_per_index_route()
        call test_class_average_route()
        call test_appended_stacks()
        call test_earlier_release()
        call test_inconsistent()
        call test_refusals()
        call test_project_file_and_merge()
        call test_isolated_stream_and_seed()
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine run_all_eo_partition_tests

    !> stacks of the given sizes, one per micrograph when l_mics; every seventh particle inactive
    subroutine make_groups_project( proj, sizes, l_mics, tag )
        class(sp_project), intent(inout) :: proj
        integer,           intent(in)    :: sizes(:)
        logical,           intent(in)    :: l_mics
        character(len=*),  intent(in)    :: tag
        integer :: s, i, j, nptcls, fromp
        call proj%kill
        nptcls = sum(sizes)
        call proj%os_stk%new(size(sizes), is_ptcl=.false.)
        if( l_mics ) call proj%os_mic%new(size(sizes), is_ptcl=.false.)
        call proj%os_ptcl2D%new(nptcls, is_ptcl=.true.)
        call proj%os_ptcl3D%new(nptcls, is_ptcl=.true.)
        fromp = 1
        i     = 0
        do s = 1, size(sizes)
            call fill_stack_row(proj%os_stk, s, tag, sizes(s), fromp)
            if( l_mics )then
                call proj%os_mic%set(s, 'intg',   '/data/'//tag//'_mic'//char(48+mod(s,10))//'.mrc')
                call proj%os_mic%set(s, 'nptcls', sizes(s))
                call proj%os_mic%set_state(s, 1)
            endif
            do j = 1, sizes(s)
                i = i + 1
                call fill_particle_row(proj%os_ptcl2D, i, s, j)
                call fill_particle_row(proj%os_ptcl3D, i, s, j)
            end do
            fromp = fromp + sizes(s)
        end do
    end subroutine make_groups_project

    subroutine fill_stack_row( os, s, tag, n, fromp )
        class(oris),      intent(inout) :: os
        integer,          intent(in)    :: s, n, fromp
        character(len=*), intent(in)    :: tag
        call os%set(s, 'stk',        '/data/'//tag//'_stk'//char(48+mod(s,10))//'.mrcs')
        call os%set(s, 'ctf',        'yes')
        call os%set(s, 'smpd',       1.25)
        call os%set(s, 'kv',         300.)
        call os%set(s, 'cs',         2.7)
        call os%set(s, 'fraca',      0.1)
        call os%set(s, 'phshift',    0.)
        call os%set(s, 'box',        64)
        call os%set(s, 'nptcls',     n)
        call os%set(s, 'nptcls_stk', n)
        call os%set(s, 'fromp',      fromp)
        call os%set(s, 'top',        fromp + n - 1)
        call os%set_state(s, 1)
        call os%set_ogid(s, 1)
    end subroutine fill_stack_row

    subroutine fill_particle_row( os, i, s, j )
        class(oris), intent(inout) :: os
        integer,     intent(in)    :: i, s, j
        call os%set_stkind(i, s)
        call os%set(i, 'indstk', j)
        call os%set_state(i, merge(0, 1, mod(i,7) == 0))
        call os%set_dfx(i, 1.5)
        call os%set_dfy(i, 1.6)
        call os%set(i, 'angast', 10.)
        call os%set(i, 'phshift', 0.)
        call os%set_ogid(i, 1)
    end subroutine fill_particle_row

    !> active particles in the halves
    subroutine active_halves( proj, pops )
        class(sp_project), intent(inout) :: proj
        integer,           intent(out)   :: pops(2)
        integer :: i, eo
        pops = 0
        do i = 1, proj%os_ptcl3D%get_noris()
            if( proj%os_ptcl3D%get_state(i) == 0 ) cycle
            eo = proj%os_ptcl3D%get_eo(i)
            if( eo == 0 .or. eo == 1 ) pops(eo+1) = pops(eo+1) + 1
        end do
    end subroutine active_halves

    !> every row of every stack (inactive rows included) carries its stack's label, in both fields
    logical function groups_whole( proj )
        class(sp_project), intent(inout) :: proj
        integer :: i, s, label
        groups_whole = .true.
        do i = 1, proj%os_ptcl3D%get_noris()
            s     = proj%os_ptcl3D%get_int(i, 'stkind')
            label = proj%os_ptcl3D%get_eo(proj%os_stk%get_fromp(s))
            if( proj%os_ptcl3D%get_eo(i) /= label ) groups_whole = .false.
            if( proj%os_ptcl2D%get_eo(i) /= proj%os_ptcl3D%get_eo(i) ) groups_whole = .false.
        end do
    end function groups_whole

    subroutine test_micrograph_groups()
        class(sp_project), allocatable :: proj
        integer, allocatable :: sizes(:)
        character(len=STDLEN) :: msg
        integer :: s, pops(2)
        logical :: l_changed
        write(*,'(A)') 'test_micrograph_groups'
        allocate(proj)
        sizes = [(1 + mod(17*s, 23), s=1,40)]
        call make_groups_project(proj, sizes, .true., 'groups')
        call assert_int(EO_MISSING, eo_partition_status(proj, msg), 'a fresh project has no partition')
        call eo_transaction(proj, l_changed)
        call assert_true(l_changed, 'the transaction labels a fresh project')
        call assert_int(EO_COMPLETE, eo_partition_status(proj, msg), 'the partition is complete: '//trim(msg))
        call assert_true(groups_whole(proj), 'one label per micrograph group, the same in ptcl2D and ptcl3D')
        call active_halves(proj, pops)
        call assert_true(abs(pops(1) - pops(2)) <= maxval(sizes), 'the halves differ by at most the largest group')
        call assert_true(minval(pops) > 0, 'both halves are populated')
        do s = 1, size(sizes)
            call assert_true(proj%os_stk%isthere(s, 'eo_group'), 'every stack row records its group half')
        end do
        call eo_transaction(proj, l_changed)
        call assert_false(l_changed, 'a complete partition is left as it is')
        call proj%kill
    end subroutine test_micrograph_groups

    subroutine test_no_split_by_size()
        class(sp_project), allocatable :: proj
        integer, allocatable :: sizes(:)
        character(len=STDLEN) :: msg
        integer :: s, pops(2)
        logical :: l_changed
        write(*,'(A)') 'test_no_split_by_size'
        allocate(proj)
        sizes = [500, (10, s=1,25)]
        call make_groups_project(proj, sizes, .true., 'large')
        call eo_transaction(proj, l_changed)
        call assert_int(EO_COMPLETE, eo_partition_status(proj, msg), 'complete with one large group: '//trim(msg))
        call assert_true(groups_whole(proj), 'the large group is not split for its size')
        call active_halves(proj, pops)
        call assert_true(abs(pops(1) - pops(2)) <= maxval(sizes), 'the imbalance is bounded by the largest group')
        call proj%kill
    end subroutine test_no_split_by_size

    !> one micrograph group with active particles cannot give two halves; the run stops (checked in a child)
    subroutine test_too_few_groups_status()
        class(sp_project), allocatable :: proj
        character(len=STDLEN) :: msg
        write(*,'(A)') 'test_too_few_groups_status'
        allocate(proj)
        call make_groups_project(proj, [30], .true., 'single')
        call assert_int(EO_MISSING, eo_partition_status(proj, msg), 'a single group project starts unlabelled')
        call proj%kill
    end subroutine test_too_few_groups_status

    subroutine test_per_index_route()
        class(sp_project), allocatable :: proj
        type(oris) :: ref
        character(len=STDLEN) :: msg
        integer :: i, n
        logical :: l_changed, l_same
        write(*,'(A)') 'test_per_index_route'
        allocate(proj)
        ! imported stacks without micrograph rows, as the beta-galactosidase test set
        call make_groups_project(proj, [57, 13, 101], .false., 'imported')
        n = proj%os_ptcl3D%get_noris()
        call ref%new(n, is_ptcl=.true.)
        call ref%partition_eo
        call eo_transaction(proj, l_changed)
        call assert_int(EO_COMPLETE, eo_partition_status(proj, msg), 'the per-index partition is complete: '//trim(msg))
        l_same = .true.
        do i = 1, n
            if( proj%os_ptcl3D%get_eo(i) /= ref%get_eo(i) ) l_same = .false.
            if( proj%os_ptcl2D%get_eo(i) /= ref%get_eo(i) ) l_same = .false.
        end do
        call assert_true(l_same, 'without micrograph rows the labels are the per-index ones, in both fields')
        call assert_true(proj%os_stk%isthere(1, 'eo_byindex'), 'the stack rows record the per-index route')
        call assert_false(proj%os_stk%isthere(1, 'eo_group'), 'no group half is recorded on the per-index route')
        call ref%kill
        call proj%kill
    end subroutine test_per_index_route

    !> class-average fields keep the per-index routine, the project's particle partition is not involved
    subroutine test_class_average_route()
        class(sp_project), allocatable :: proj
        type(oris) :: ref
        integer :: i
        logical :: l_same
        write(*,'(A)') 'test_class_average_route'
        allocate(proj)
        call make_groups_project(proj, [20, 20, 20], .true., 'cavgs')
        call proj%os_cls3D%new(9, is_ptcl=.false.)
        do i = 1, 9
            call proj%os_cls3D%set_class(i, i)
            call proj%os_cls3D%set_state(i, 1)
        end do
        call proj%update_projinfo(string(EO_PROJ))
        call proj%write(string(EO_PROJ))
        call ensure_field_halves(proj, 'cls3D', string(EO_PROJ))
        call ref%new(9, is_ptcl=.false.)
        call ref%partition_eo
        l_same = .true.
        do i = 1, 9
            if( proj%os_cls3D%get_eo(i) /= ref%get_eo(i) ) l_same = .false.
        end do
        call assert_true(l_same, 'the class-average field takes the per-index routine')
        call assert_int(0, proj%os_ptcl3D%get_nodd(), 'the particle fields are left alone')
        call assert_false(proj%os_stk%isthere(1, 'eo_group'), 'no stack row is labelled for class averages')
        call ref%kill
        call proj%kill
    end subroutine test_class_average_route

    !> stacks appended to a partitioned project get labels by the same rule; existing labels stay
    subroutine test_appended_stacks()
        class(sp_project), allocatable :: proj, more
        integer, allocatable :: sizes(:), more_sizes(:), before(:)
        character(len=STDLEN) :: msg
        integer :: s, i, n0, nstk0
        logical :: l_changed, l_kept
        write(*,'(A)') 'test_appended_stacks'
        allocate(proj, more)
        sizes      = [(5 + mod(11*s, 13), s=1,24)]
        more_sizes = [(3 + mod(7*s, 9),   s=1,12)]
        call make_groups_project(proj, sizes, .true., 'first')
        call eo_transaction(proj, l_changed)
        n0    = proj%os_ptcl3D%get_noris()
        nstk0 = proj%os_stk%get_noris()
        before = [(proj%os_ptcl3D%get_eo(i), i=1,n0)]
        ! the appended rows, unlabelled, as a stream ingestion or an import adds them
        call make_groups_project(more, more_sizes, .true., 'second')
        call append_rows(proj, more, n0)
        call assert_int(EO_PARTIAL, eo_partition_status(proj, msg), 'an appended project is partly labelled')
        call eo_transaction(proj, l_changed)
        call assert_true(l_changed, 'the appended stacks are labelled')
        call assert_int(EO_COMPLETE, eo_partition_status(proj, msg), 'complete after the append: '//trim(msg))
        l_kept = .true.
        do i = 1, n0
            if( proj%os_ptcl3D%get_eo(i) /= before(i) ) l_kept = .false.
        end do
        call assert_true(l_kept, 'existing labels are untouched')
        call assert_true(groups_whole(proj), 'appended groups carry one label each, the same in both fields')
        do s = nstk0 + 1, proj%os_stk%get_noris()
            call assert_true(proj%os_stk%isthere(s, 'eo_group'), 'appended stack rows record their half')
        end do
        call more%kill
        call proj%kill
    end subroutine test_appended_stacks

    !> the rows of more appended to proj, the stack and particle indices shifted
    subroutine append_rows( proj, more, n0 )
        class(sp_project), intent(inout) :: proj, more
        integer,           intent(in)    :: n0
        integer :: nstk0, s, i
        nstk0 = proj%os_stk%get_noris()
        do s = 1, more%os_stk%get_noris()
            call more%os_stk%set(s, 'fromp', more%os_stk%get_fromp(s) + n0)
            call more%os_stk%set(s, 'top',   more%os_stk%get_top(s)   + n0)
        end do
        do i = 1, more%os_ptcl3D%get_noris()
            call more%os_ptcl2D%set_stkind(i, more%os_ptcl2D%get_int(i, 'stkind') + nstk0)
            call more%os_ptcl3D%set_stkind(i, more%os_ptcl3D%get_int(i, 'stkind') + nstk0)
        end do
        call proj%os_mic%append(more%os_mic)
        call proj%os_stk%append(more%os_stk)
        call proj%os_ptcl2D%append(more%os_ptcl2D)
        call proj%os_ptcl3D%append(more%os_ptcl3D)
    end subroutine append_rows

    !> a project of an earlier release: labels (in one field), no provenance; the labels are kept
    subroutine test_earlier_release()
        class(sp_project), allocatable :: proj
        integer, allocatable :: before(:)
        character(len=STDLEN) :: msg
        integer :: i, n, s
        logical :: l_changed, l_kept
        write(*,'(A)') 'test_earlier_release'
        allocate(proj)
        call make_groups_project(proj, [(9, s=1,30)], .true., 'legacy')
        n = proj%os_ptcl3D%get_noris()
        do i = 1, n
            call proj%os_ptcl3D%set(i, 'eo', real(mod(i,2)))
        end do
        before = [(proj%os_ptcl3D%get_eo(i), i=1,n)]
        call assert_int(EO_PARTIAL, eo_partition_status(proj, msg), 'labels without provenance are not complete')
        call eo_transaction(proj, l_changed)
        call assert_true(l_changed, 'the provenance is recorded')
        l_kept = .true.
        do i = 1, n
            if( proj%os_ptcl3D%get_eo(i) /= before(i) ) l_kept = .false.
            if( proj%os_ptcl2D%get_eo(i) /= before(i) ) l_kept = .false.
        end do
        call assert_true(l_kept, 'the earlier labels are kept and copied to the unlabelled field')
        call assert_true(proj%os_stk%isthere(1, 'eo_byindex'), 'the earlier labels are recorded as per-index')
        call assert_int(EO_COMPLETE, eo_partition_status(proj, msg), 'complete once recorded: '//trim(msg))
        call proj%kill
    end subroutine test_earlier_release

    subroutine test_inconsistent()
        class(sp_project), allocatable :: proj
        character(len=STDLEN) :: msg
        integer :: s, i
        logical :: l_changed
        write(*,'(A)') 'test_inconsistent'
        allocate(proj)
        call make_groups_project(proj, [(8, s=1,24)], .true., 'bad')
        call eo_transaction(proj, l_changed)
        ! an active particle moved to the other half of its group
        i = 1
        call proj%os_ptcl3D%set(i, 'eo', real(1 - proj%os_ptcl3D%get_eo(i)))
        call assert_int(EO_INCONSISTENT, eo_partition_status(proj, msg), 'a particle against its group half')
        call proj%os_ptcl3D%set(i, 'eo', real(proj%os_ptcl2D%get_eo(i)))
        call assert_int(EO_COMPLETE, eo_partition_status(proj, msg), 'consistent again: '//trim(msg))
        ! fields that are not congruent
        call proj%os_ptcl2D%new(10, is_ptcl=.true.)
        call assert_int(EO_INCONSISTENT, eo_partition_status(proj, msg), 'ptcl2D and ptcl3D of different rows')
        call proj%kill
    end subroutine test_inconsistent

    !> a worker without halves and a master with an inconsistent partition stop rather than make one
    subroutine test_refusals()
        class(sp_project), allocatable :: proj
        type(refusing_partition) :: child
        integer :: s
        logical :: l_changed
        write(*,'(A)') 'test_refusals'
        allocate(proj)
        call make_groups_project(proj, [(6, s=1,24)], .true., 'refuse')
        call proj%update_projinfo(string(EO_PROJ))
        call proj%write(string(EO_PROJ))
        child%what = REFUSE_WORKER_MISSING
        call child%start(name=string('EO_WORKER_MISSING'))
        call child%await_final_status()
        call assert_int(FORK_STATUS_FAILED, child%status(), 'a worker refuses a project without halves')
        call child%destroy
        call eo_transaction(proj, l_changed)
        call proj%os_ptcl2D%set(2, 'eo', real(1 - proj%os_ptcl2D%get_eo(2)))
        call proj%write(string(EO_PROJ))
        child%what = REFUSE_MASTER_INCONSISTENT
        call child%start(name=string('EO_MASTER_INCONSISTENT'))
        call child%await_final_status()
        call assert_int(FORK_STATUS_FAILED, child%status(), 'the master refuses an inconsistent partition')
        call child%destroy
        ! a single micrograph group cannot give two halves
        call make_groups_project(proj, [40], .true., 'one')
        call proj%update_projinfo(string(EO_PROJ))
        call proj%write(string(EO_PROJ))
        child%what = REFUSE_MASTER_INCONSISTENT
        call child%start(name=string('EO_ONE_GROUP'))
        call child%await_final_status()
        call assert_int(FORK_STATUS_FAILED, child%status(), 'one micrograph group gives no two halves')
        call child%destroy
        call proj%kill
    end subroutine test_refusals

    subroutine refuse_partition( self, cline )
        class(refusing_partition), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(sp_project), allocatable :: proj
        allocate(proj)
        call proj%read(string(EO_PROJ))
        select case(self%what)
            case(REFUSE_WORKER_MISSING)
                l_distr_worker_glob = .true.
                call ensure_field_halves(proj, 'ptcl3D', string(EO_PROJ))
            case(REFUSE_MASTER_INCONSISTENT)
                l_distr_worker_glob = .false.
                call ensure_field_halves(proj, 'ptcl3D', string(EO_PROJ))
        end select
        ! reached only when the partition was not refused: a clean exit fails the test
        call proj%kill
    end subroutine refuse_partition

    !> the transaction through the project file, a field loaded alone, and the labels after a merge
    subroutine test_project_file_and_merge()
        class(sp_project), allocatable :: proj, part2D, reread, merged
        type(string), allocatable :: files(:)
        integer, allocatable :: eo1(:), eo2(:)
        character(len=STDLEN) :: msg
        integer :: s, i, n1, n2
        logical :: l_same
        write(*,'(A)') 'test_project_file_and_merge'
        allocate(proj, part2D, reread, merged)
        ! a 2D run loads ptcl2D and stk only: the transaction reads the rest from the file
        call make_groups_project(proj, [(4 + mod(5*s, 7), s=1,22)], .true., 'filea')
        call proj%update_projinfo(string('eo_merge_a.simple'))
        call proj%write(string('eo_merge_a.simple'))
        call part2D%read_segment('ptcl2D', string('eo_merge_a.simple'))
        call part2D%read_segment('stk',    string('eo_merge_a.simple'))
        call ensure_field_halves(part2D, 'ptcl2D', string('eo_merge_a.simple'))
        call reread%read(string('eo_merge_a.simple'))
        call assert_int(EO_COMPLETE, eo_partition_status(reread, msg), 'the file holds a complete partition: '//trim(msg))
        call assert_true(groups_whole(reread), 'the file holds group labels in both fields')
        n1 = reread%os_ptcl3D%get_noris()
        l_same = .true.
        do i = 1, n1
            if( part2D%os_ptcl2D%get_eo(i) /= reread%os_ptcl2D%get_eo(i) ) l_same = .false.
        end do
        call assert_true(l_same, 'the loaded field receives the labels the file holds')
        call assert_true(part2D%os_stk%isthere(1, 'eo_group'), 'the loaded stack rows receive the provenance')
        eo1 = [(reread%os_ptcl3D%get_eo(i), i=1,n1)]
        ! a second partitioned project, then the merge
        call make_groups_project(proj, [(3 + mod(3*s, 8), s=1,21)], .true., 'fileb')
        call proj%update_projinfo(string('eo_merge_b.simple'))
        call proj%write(string('eo_merge_b.simple'))
        call ensure_field_halves(proj, 'ptcl3D', string('eo_merge_b.simple'))
        n2  = proj%os_ptcl3D%get_noris()
        eo2 = [(proj%os_ptcl3D%get_eo(i), i=1,n2)]
        allocate(files(2))
        files(1) = 'eo_merge_a.simple'
        files(2) = 'eo_merge_b.simple'
        call merge_selected_project_files(files, string('eo_merge_ab.simple'), merged, write_proj=.true.)
        call reread%kill
        call reread%read(string('eo_merge_ab.simple'))
        call assert_int(n1 + n2, reread%os_ptcl3D%get_noris(), 'the merged project holds both')
        l_same = .true.
        do i = 1, n1
            if( reread%os_ptcl3D%get_eo(i) /= eo1(i) .or. reread%os_ptcl2D%get_eo(i) /= eo1(i) ) l_same = .false.
        end do
        do i = 1, n2
            if( reread%os_ptcl3D%get_eo(n1+i) /= eo2(i) .or. reread%os_ptcl2D%get_eo(n1+i) /= eo2(i) ) l_same = .false.
        end do
        call assert_true(l_same, 'after the merge ptcl2D and ptcl3D carry the sources'' labels')
        call ensure_field_halves(reread, 'ptcl3D', string('eo_merge_ab.simple'))
        l_same = .true.
        do i = 1, n1
            if( reread%os_ptcl3D%get_eo(i) /= eo1(i) ) l_same = .false.
        end do
        do i = 1, n2
            if( reread%os_ptcl3D%get_eo(n1+i) /= eo2(i) ) l_same = .false.
        end do
        call assert_true(l_same, 'the merged partition is kept by the next run')
        call assert_int(EO_COMPLETE, eo_partition_status(reread, msg), 'the merged partition is complete: '//trim(msg))
        call proj%kill
        call part2D%kill
        call reread%kill
        call merged%kill
    end subroutine test_project_file_and_merge

    !> making the split draws nothing from the process-wide stream; a fixed seed gives the same split
    subroutine test_isolated_stream_and_seed()
        class(sp_project), allocatable :: proj
        integer, allocatable :: sizes(:), eo_a(:), eo_b(:), eo_c(:)
        real    :: draws_ref(5), draws(5)
        integer :: s, i, n, k
        logical :: l_changed
        write(*,'(A)') 'test_isolated_stream_and_seed'
        allocate(proj)
        sizes = [(2 + mod(13*s, 17), s=1,60)]
        call set_fixed_seed(4711)
        draws_ref = [(ran3(), k=1,5)]
        call make_groups_project(proj, sizes, .true., 'seed')
        call set_fixed_seed(4711)
        call eo_transaction(proj, l_changed)
        draws = [(ran3(), k=1,5)]
        call assert_true(all(draws == draws_ref), 'the split leaves the particle sampling stream unchanged')
        n    = proj%os_ptcl3D%get_noris()
        eo_a = [(proj%os_ptcl3D%get_eo(i), i=1,n)]
        call make_groups_project(proj, sizes, .true., 'seed')
        call eo_transaction(proj, l_changed)
        eo_b = [(proj%os_ptcl3D%get_eo(i), i=1,n)]
        call assert_true(all(eo_a == eo_b), 'the fixed seed gives the same split')
        call make_groups_project(proj, sizes, .true., 'seed')
        call eo_transaction(proj, l_changed, seed=7_int64)
        eo_c = [(proj%os_ptcl3D%get_eo(i), i=1,n)]
        call assert_true(any(eo_a /= eo_c), 'another seed gives another split')
        call proj%kill
    end subroutine test_isolated_stream_and_seed

end module simple_eo_partition_tester
