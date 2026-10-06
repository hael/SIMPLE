!@descr: unit tests for the stream's 2D pool type (simple_stream_pool2D) and its stateless helpers
! A pool is built without a queue: new, and for the tests that need its dimensions the state part of
! start (init_state) in a fresh fixture directory. No iteration is submitted; the iterations, the
! history and the dimension updates are covered by the high-level stream tests. The rows a set
! adds are checked on a plain project (append_project_sets, which the pool's append_sets calls).
! Sub-suite "pool 2D object" of unit_stream.
module simple_stream_pool2D_tester
use simple_test_utils
use simple_defs,                  only: COSMSKHALFWIDTH
use simple_string,                only: string
use simple_fileio,                only: file_exists
use simple_cmdline,               only: cmdline
use simple_parameters,            only: parameters
use simple_sp_project,            only: sp_project
use simple_stream_utils,          only: create_stream_project
use simple_stream_refine2D_utils, only: append_project_sets, draw_new_classes
use simple_stream_pool2D,         only: stream_pool2D, stream_pool2D_stats
implicit none
private
public :: run_all_stream_pool2D_tests

character(len=*), parameter :: TEST_PROJFILE = 'test_pool2D_object.simple'

contains

    subroutine run_all_stream_pool2D_tests()
        write(*,'(A)') '**** running all stream pool 2D object tests ****'
        call test_new_kill()
        call test_append_project_sets()
        call test_append_sets()
        call test_draw_new_classes()
        call test_mask_clamped_to_box()
        call test_restart_after_kill()
        call test_unclassified_publishes_nothing()
        call test_snapshot_without_frcs()
    end subroutine run_all_stream_pool2D_tests

    !> an empty pool answers as before any iteration, and kill leaves an empty pool, twice over
    subroutine test_new_kill()
        class(stream_pool2D), allocatable :: pool
        type(stream_pool2D_stats) :: st
        allocate(pool)
        write(*,'(A)') 'test_new_kill'
        call pool%new()
        st = pool%stats()
        call assert_int(0, pool%iteration(),         'an empty pool has run no iteration')
        call assert_int(0, st%last_complete_iter,    'and completed none')
        call assert_false(pool%available(),          'it is not available before it starts')
        call assert_false(pool%failed(),             'nor failed')
        call assert_int(0, st%nassigned,             'it holds no particle')
        call assert_real(999., st%resolution, 1.e-4, 'and has no resolution')
        call assert_int(0, st%cavgs_mrc%strlen(),    'no class averages')
        call assert_false(allocated(st%jpeg_map),    'and no sprite sheet')
        call pool%kill()
        call assert_int(0, pool%iteration(),         'a killed pool is empty')
        call pool%kill() ! idempotence
    end subroutine test_new_kill

    !> the sets are appended in order, with their stacks renumbered after the project's and their
    !! particles as new ones pointing at their stack; a set without micrographs adds nothing
    subroutine test_append_project_sets()
        type(sp_project) :: proj, sets(3)
        write(*,'(A)') 'test_append_project_sets'
        call make_set(sets(1), [3, 2], nrejected=1)
        call make_set(sets(2), [4])
        call append_project_sets(proj, sets(1:2))
        call assert_int(3, proj%os_mic%get_noris(),      'their micrographs are in the project')
        call assert_int(3, proj%os_stk%get_noris(),      'with a stack each')
        call assert_int(9, proj%os_ptcl2D%get_noris(),   'and every particle')
        call assert_int(4, proj%os_stk%get_fromp(2),     'the second stack follows the first')
        call assert_int(6, proj%os_stk%get_fromp(3),     'the second set''s stack follows the first set''s')
        call assert_int(9, proj%os_stk%get_top(3),       'to the last particle')
        call assert_int(3, proj%os_ptcl2D%get_int(7, 'stkind'),    'a particle points at its stack in the project')
        call assert_int(0, proj%os_ptcl2D%get_class(7),            'the particles come without classes')
        call assert_int(0, proj%os_ptcl2D%get_int(7, 'updatecnt'), 'as new particles')
        call assert_real(1.5, proj%os_ptcl2D%get(7, 'x'), 1.e-4,   'keeping their shifts')
        call assert_int(0, proj%os_ptcl2D%get_state(1),            'and their selection')
        call make_set(sets(3), [2])
        call append_project_sets(proj, sets(3:3))
        call assert_int(10, proj%os_stk%get_fromp(4),    'a later set comes after the particles already there')
        call assert_int(11, proj%os_ptcl2D%get_noris(),  'the project grows')
        call sets(3)%kill
        call append_project_sets(proj, sets(3:3))
        call assert_int(4, proj%os_mic%get_noris(),      'a set without micrographs adds nothing')
        call proj%kill
        call sets(1)%kill
        call sets(2)%kill
    end subroutine test_append_project_sets

    !> the pool takes the sets of an import and reports its micrographs and selected particles
    subroutine test_append_sets()
        class(stream_pool2D), allocatable :: pool
        type(sp_project) :: sets(2)
        integer          :: nmics, nsel
        allocate(pool)
        write(*,'(A)') 'test_append_sets'
        call pool%new()
        call make_set(sets(1), [3, 2], nrejected=1)
        call make_set(sets(2), [4])
        call pool%append_sets(sets, nmics, nsel)
        call assert_int(3, nmics, 'the pool''s micrographs')
        call assert_int(8, nsel,  'and its selected particles')
        call sets(1)%kill
        call make_set(sets(1), [2])
        call pool%append_sets(sets(1:1), nmics, nsel)
        call assert_int(4,  nmics, 'a later import adds to them')
        call assert_int(10, nsel,  'with its selected particles')
        call pool%kill()
        call sets(1)%kill
        call sets(2)%kill
    end subroutine test_append_sets

    !> new particles get a populated class, the same ones for the same iteration; particles
    !! already updated, and deselected ones, keep theirs
    subroutine test_draw_new_classes()
        type(sp_project) :: proj1, proj2
        integer :: i
        logical :: l_same, l_populated
        write(*,'(A)') 'test_draw_new_classes'
        call make_draw_project(proj1)
        call make_draw_project(proj2)
        call draw_new_classes(proj1, 12, 7, 4)
        call draw_new_classes(proj2, 12, 7, 4)
        l_same      = .true.
        l_populated = .true.
        do i = 1,12
            l_same = l_same .and. proj1%os_ptcl2D%get_class(i) == proj2%os_ptcl2D%get_class(i)
            if( i > 6 .and. i < 12 ) l_populated = l_populated .and. any(proj1%os_ptcl2D%get_class(i) == [2, 4])
        enddo
        call assert_true(l_same,      'the same iteration draws the same classes')
        call assert_true(l_populated, 'only populated classes are drawn')
        call assert_int(1, proj1%os_ptcl2D%get_class(1),  'an updated particle keeps its class')
        call assert_int(3, proj1%os_ptcl2D%get_class(12), 'a deselected particle keeps its class')
        call proj1%kill
        call proj2%kill

    contains

        ! 4 classes, populated 2 and 4; particles 1-6 updated (class 1), 7-12 new (class 3), 12 deselected
        subroutine make_draw_project( proj )
            type(sp_project), intent(inout) :: proj
            integer :: j
            call proj%os_cls2D%new(4, is_ptcl=.false.)
            call proj%os_cls2D%set_all('pop', real([0, 5, 0, 3]))
            call proj%os_ptcl2D%new(12, is_ptcl=.true.)
            do j = 1,12
                call proj%os_ptcl2D%set_state(j, 1)
                if( j <= 6 )then
                    call proj%os_ptcl2D%set(j, 'updatecnt', 1)
                    call proj%os_ptcl2D%set_class(j, 1)
                else
                    call proj%os_ptcl2D%set(j, 'updatecnt', 0)
                    call proj%os_ptcl2D%set_class(j, 3)
                endif
            enddo
            call proj%os_ptcl2D%set_state(12, 0)
        end subroutine make_draw_project

    end subroutine test_draw_new_classes

    !> a mask diameter beyond the pool's box reaches its command line clamped (D40), and the
    !! low-pass ramp follows the mask; a 64 px box is below CHUNK_MINBOXSZ, so the pool keeps its
    !! native dimensions
    subroutine test_mask_clamped_to_box()
        class(stream_pool2D), allocatable :: pool
        class(parameters),    allocatable :: params
        type(stream_pool2D_stats) :: st
        type(cmdline)    :: cline
        type(sp_project) :: spproj
        type(string)     :: cwd_saved, root
        real             :: mskdiam_max
        integer          :: nfail0
        allocate(pool, params)
        write(*,'(A)') 'test_mask_clamped_to_box'
        nfail0 = tests_failed
        call enter_fixture('pool2D_object_mask', cwd_saved, root)
        call make_started_pool(pool, params, cline, spproj, 100.)
        st = pool%stats()
        call assert_real(100., st%mskdiam, 1.e-4,  'a started pool has the mask it was given')
        call assert_true(pool%available(),         'and is available for an iteration')
        mskdiam_max = (64. - COSMSKHALFWIDTH) * 2.0
        call pool%set_mskdiam(1000)
        st = pool%stats()
        call assert_true(st%mskdiam  <= mskdiam_max + 1.e-3,          'the diameter is clamped to the box')
        call assert_true(st%msk_crop <= (64. - COSMSKHALFWIDTH) / 2., 'and the radius')
        call pool%set_mskdiam(100)
        st = pool%stats()
        call assert_real(100., st%mskdiam, 1.e-4,                      'a diameter within the box is kept')
        call pool%kill()
        call spproj%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_mask_clamped_to_box

    !> kill resets a started pool; it can be made and started again
    subroutine test_restart_after_kill()
        class(stream_pool2D), allocatable :: pool
        class(parameters),    allocatable :: params
        type(stream_pool2D_stats) :: st
        type(cmdline)    :: cline
        type(sp_project) :: spproj
        type(string)     :: cwd_saved, root
        integer          :: nfail0
        allocate(pool, params)
        write(*,'(A)') 'test_restart_after_kill'
        nfail0 = tests_failed
        call enter_fixture('pool2D_object_restart', cwd_saved, root)
        call make_started_pool(pool, params, cline, spproj, 100.)
        call pool%kill()
        st = pool%stats()
        call assert_false(pool%available(),         'a killed pool is not available')
        call assert_real(0., st%mskdiam, 1.e-4,     'and has no mask')
        call spproj%kill
        call cline%kill
        deallocate(params)
        allocate(params)
        call make_started_pool(pool, params, cline, spproj, 80.)
        st = pool%stats()
        call assert_true(pool%available(),          'it starts again')
        call assert_real(80., st%mskdiam, 1.e-4,    'with the new mask')
        call assert_int(0, pool%iteration(),        'and no iteration')
        call pool%kill()
        call spproj%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_restart_after_kill

    !> a pool whose particles have never been through an iteration publishes nothing for 3D
    subroutine test_unclassified_publishes_nothing()
        class(stream_pool2D), allocatable :: pool
        type(sp_project) :: sets(1)
        type(string)     :: cwd_saved, root
        integer          :: nfail0, nmics, nsel, nstks
        allocate(pool)
        write(*,'(A)') 'test_unclassified_publishes_nothing'
        nfail0 = tests_failed
        call enter_fixture('pool2D_object_publish', cwd_saved, root)
        call pool%new()
        call make_set(sets(1), [3])
        call pool%append_sets(sets, nmics, nsel)
        call pool%publish(string('00001.simple'), nstks, string(''), .false.)
        call assert_int(0, nstks,                                 'no stack is published')
        call assert_false(file_exists(string('00001.simple')),    'nor a project written')
        call pool%kill()
        call sets(1)%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_unclassified_publishes_nothing

    !> a snapshot of the current iteration before its frcs.bin exists is reported as not written
    !! (no particles, no file) instead of stopping the stage; an empty pool's current iteration is 0
    subroutine test_snapshot_without_frcs()
        class(stream_pool2D), allocatable :: pool
        type(string)         :: cwd_saved, root, jpeg, mrc
        integer, allocatable :: idx(:), pop(:)
        real,    allocatable :: res(:)
        integer              :: nfail0, nptcls, ntx, nty
        allocate(pool)
        write(*,'(A)') 'test_snapshot_without_frcs'
        nfail0 = tests_failed
        call enter_fixture('pool2D_object_snapshot', cwd_saved, root)
        call pool%new()
        call pool%write_snapshot(0, [1], string('snapshots/snap/snap.simple'), string('snapshots/snap/snap'), string(''), 0,&
            &nptcls, jpeg, mrc, ntx, nty, idx, pop, res)
        call assert_int(0, nptcls, 'nothing written')
        call assert_false(file_exists(string('snapshots/snap/snap.simple')), 'no project')
        call pool%kill()
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_snapshot_without_frcs

    ! ---- fixtures ------------------------------------------------------------

    ! a command line under a program name outside every UI table, with the refinement the pool
    ! accepts and a local queue system for the project's computing environment
    subroutine set_test_cline( cline )
        type(cmdline), intent(inout) :: cline
        call cline%set('prg',       'stream_pool2D_tester')
        call cline%set('mkdir',     'no')
        call cline%set('projfile',  TEST_PROJFILE)
        call cline%set('outdir',    '')
        call cline%set('nthr',      1)
        call cline%set('nparts',    1)
        call cline%set('ncls',      10)
        call cline%set('refine',    'snhc')
        call cline%set('qsys_name', 'local')
    end subroutine set_test_cline

    ! a new pool started without a queue (init_state) on 64 px particles at 2 A, with the mask
    ! diameter @p mskdiam (A), in the current directory
    subroutine make_started_pool( pool, params, cline, spproj, mskdiam )
        class(stream_pool2D), intent(inout) :: pool
        class(parameters),    intent(inout) :: params
        type(cmdline),        intent(inout) :: cline
        type(sp_project),     intent(inout) :: spproj
        real,                 intent(in)    :: mskdiam
        call set_test_cline(cline)
        call create_stream_project(spproj, cline, string('pool2D'))
        call params%new(cline)
        call pool%new()
        call pool%init_state(params, cline, spproj, 64, 2.0, mskdiam)
    end subroutine make_started_pool

    ! a sieve set in memory: one micrograph and one stack per entry of @p nptcls_mic, its particles
    ! selected, in class 3 with a shift, the first @p nrejected of them deselected
    subroutine make_set( proj, nptcls_mic, nrejected )
        type(sp_project),  intent(inout) :: proj
        integer,           intent(in)    :: nptcls_mic(:)
        integer, optional, intent(in)    :: nrejected
        integer :: imic, fromp, iptcl
        call proj%os_mic%new(size(nptcls_mic), is_ptcl=.false.)
        call proj%os_stk%new(size(nptcls_mic), is_ptcl=.false.)
        call proj%os_ptcl2D%new(sum(nptcls_mic), is_ptcl=.true.)
        fromp = 1
        do imic = 1,size(nptcls_mic)
            call proj%os_mic%set_state(imic, 1)
            call proj%os_mic%set(imic, 'nptcls', nptcls_mic(imic))
            call proj%os_stk%set(imic, 'nptcls', nptcls_mic(imic))
            call proj%os_stk%set(imic, 'fromp',  fromp)
            call proj%os_stk%set(imic, 'top',    fromp + nptcls_mic(imic) - 1)
            do iptcl = fromp,fromp + nptcls_mic(imic) - 1
                call proj%os_ptcl2D%set_stkind(iptcl, imic)
                call proj%os_ptcl2D%set_state(iptcl, 1)
                call proj%os_ptcl2D%set_class(iptcl, 3)
                call proj%os_ptcl2D%set(iptcl, 'x', 1.5)
            enddo
            fromp = fromp + nptcls_mic(imic)
        enddo
        if( present(nrejected) )then
            do iptcl = 1,nrejected
                call proj%os_ptcl2D%set_state(iptcl, 0)
            enddo
        endif
    end subroutine make_set

end module simple_stream_pool2D_tester
