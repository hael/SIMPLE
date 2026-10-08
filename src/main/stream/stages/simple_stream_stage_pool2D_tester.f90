!@descr: unit tests for the steps of the stream pool 2D stage (simple_stream_stage_pool2D)
! Each test assembles a stage in a fresh fixture directory from init_params and init_gui, with no
! waits and a settle time of -1. The command line names no registered program and carries
! qsys_name=local for the stage project's computing environment. The upstream is a sieve directory
! whose completed folder holds handed-off sets. The pool itself is not started (its iterations need
! a queue): sets are handed to the stage's empty pool and checked through the stage's counts; the
! pool object's own tests, rows included, are simple_stream_pool2D_tester (its own sub-suite), and
! the pool runs are covered by the high-level stream tests.
module simple_stream_stage_pool2D_tester
use, intrinsic :: iso_c_binding, only: c_int
use unix,                        only: c_pipe, c_close, c_fcntl, F_GETFL, F_SETFL, O_NONBLOCK
use simple_test_utils
use simple_defs_fname,                          only: TERM_STREAM, METADATA_EXT, USER_PARAMS2D, REFINE2D_FINISHED
use simple_defs_stream,                         only: DIR_STREAM_COMPLETED, POOL_EXIT_CODE, POOL_INPUT_PROJFILE, NPTCLS_FIRST3D
use simple_string,                              only: string
use simple_string_utils,                        only: int2str_pad
use simple_fileio,                              only: basename, del_file, file_exists, simple_getcwd, simple_touch, add2fbody
use simple_rec_list,                            only: rec_iterator, chunk_rec
use simple_syslib,                              only: dir_exists, simple_mkdir
use simple_cmdline,                             only: cmdline
use simple_sp_project,                          only: sp_project
use simple_gui_metadata_utils,                  only: max_metadata_size
use simple_gui_metadata_types,                  only: GUI_METADATA_STREAM_POOL2D_TYPE, GUI_METADATA_STREAM_POOL2D_SNAPSHOT_TYPE,&
                                                     &GUI_METADATA_STREAM_POOL2D_SNAPSHOT_CLS2D_TYPE, GUI_METADATA_STREAM_UPDATE_TYPE
use simple_gui_metadata_stream_pool2D,          only: gui_metadata_stream_pool2D
use simple_gui_metadata_stream_snapshot,        only: gui_metadata_stream_snapshot
use simple_gui_metadata_stream_update,          only: gui_metadata_stream_update
use simple_stream_pipe,                         only: stream_pipe
use simple_stream_stage_pool2D,                 only: stream_stage_pool2D
use simple_stream_refine2D_utils,              only: build_pool_publication, pool_publication_nselected,&
                                                     &build_sieve_publication, combine_sieve_classes
use simple_image,                               only: image
use simple_class_frcs,                          only: class_frcs
use simple_imghead,                             only: find_ldim_nptcls
use simple_optics_maps,                         only: publish_optics_map
implicit none
private
public :: run_all_stream_stage_pool2D_tests

character(len=*), parameter :: TEST_PROJFILE = 'test_pool2D.simple'
character(len=*), parameter :: UPSTREAM      = 'sieve'  ! the sieve stage directory (dir_target)
real,             parameter :: MSKDIAM_SIEVE = 180.     ! the mask diameter of the sieve's 2D

contains

    subroutine run_all_stream_stage_pool2D_tests()
        write(*,'(A)') '**** running all stream pool 2D stage tests ****'
        call test_init_params()
        call test_restart_cleans()
        call test_export_numbering()
        call test_publication_holds_classified_stacks()
        call test_sieve_publication()
        call test_attach_and_watch()
        call test_sets_taken_in_order()
        call test_transfer_sets()
        call test_transfer_stepwise()
        call test_sieve_final_set()
        call test_pause_rules()
        call test_gui_mskdiam_update()
        call test_given_mskdiam_kept()
        call test_send_status()
        call test_send_snapshot()
        call test_iterate_waits()
        call test_finished()
    end subroutine run_all_stream_stage_pool2D_tests

    subroutine test_init_params()
        class(stream_stage_pool2D), allocatable :: stage
        type(cmdline)             :: cline
        type(string)              :: cwd_saved, root, val
        integer                   :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_init_params'
        nfail0 = tests_failed
        call enter_fixture('p2_stage_init_params', cwd_saved, root)
        call set_test_cline(cline)
        call cline%set('nicedispid', 3)
        call cline%set('stepwise',   'yes')
        call cline%set('mkdir',      'yes') ! no execution directory: the program is not registered
        call make_test_stage(stage, cline)
        call assert_true(allocated(stage%params), 'initialization allocates the owned parameters')
        call assert_int(0, stage%spproj%os_mic%get_noris(),        'the project starts without micrographs')
        call assert_true(dir_exists(string(DIR_STREAM_COMPLETED)), 'the completed folder is made')
        call assert_false(stage%l_restart,                         'a fresh run is no restart')
        call assert_true(stage%l_stepwise,                         'stepwise is read from its parameter')
        call assert_int(1000, stage%optics_id_offset,              'optics ids are offset per GUI display')
        call assert_true(associated(stage%cline),                  'the stage keeps its command line')
        val = stage%cline%get_carg('mkdir')
        call assert_char('no', val%to_char(),                      'with mkdir=no for the pool')
        val = cline%get_carg('mkdir')
        call assert_char('yes', val%to_char(),                     'the caller''s command line is not changed')
        call assert_false(stage%l_pool_started,                    'no pool before the first import')
        call stage%kill
        call assert_false(allocated(stage%params), 'cleanup releases the owned parameters')
        call assert_false(associated(stage%cline),                 'kill releases the command line')
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_init_params

    !> restart: the previous pool's files go, its iteration's completion and exit status with them,
    !! and the stage's project and the snapshots stay
    subroutine test_restart_cleans()
        class(stream_stage_pool2D), allocatable :: stage
        type(cmdline)             :: cline
        type(string)              :: cwd_saved, root
        integer                   :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_restart_cleans'
        nfail0 = tests_failed
        call enter_fixture('p2_stage_restart', cwd_saved, root)
        call simple_mkdir('previous')
        call simple_mkdir('snapshots')
        call simple_touch(TERM_STREAM)
        call simple_touch(USER_PARAMS2D)
        call simple_touch('cavgs_iter003.jpg')
        call simple_touch(REFINE2D_FINISHED)
        call simple_touch(POOL_EXIT_CODE)
        call simple_touch(POOL_INPUT_PROJFILE)
        call set_test_cline(cline)
        call cline%set('outdir', 'previous')
        call make_test_stage(stage, cline)
        call assert_true(stage%l_restart, 'an existing output folder is a restart')
        call stage%clean_previous_run()
        call assert_false(file_exists(TERM_STREAM),            'the termination file is removed')
        call assert_false(file_exists(USER_PARAMS2D),          'the user parameters are removed')
        call assert_false(file_exists('cavgs_iter003.jpg'),    'the pool''s images are removed')
        call assert_false(file_exists(REFINE2D_FINISHED),      'the previous iteration''s completion is removed')
        call assert_false(file_exists(POOL_EXIT_CODE),         'and its exit status')
        call assert_false(file_exists(POOL_INPUT_PROJFILE),    'and its project as made')
        call assert_true(file_exists(TEST_PROJFILE),           'the stage''s project stays')
        call assert_true(dir_exists(string('snapshots')),      'the snapshots stay')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_restart_cleans

    !> exports for 3D continue after the highest one in the completed folder
    subroutine test_export_numbering()
        class(stream_stage_pool2D), allocatable :: stage
        type(cmdline)             :: cline
        type(string)              :: cwd_saved, root
        integer                   :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_export_numbering'
        nfail0 = tests_failed
        call enter_fixture('p2_stage_export_ids', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%restore_export_id()
        call assert_int(1, stage%last_export_id, 'no export yet: the first is number 1')
        call simple_touch(DIR_STREAM_COMPLETED//'00003'//METADATA_EXT)
        call simple_touch(DIR_STREAM_COMPLETED//'00007'//METADATA_EXT)
        call simple_touch(DIR_STREAM_COMPLETED//'00008.tmp')
        call stage%restore_export_id()
        call assert_int(8, stage%last_export_id, 'a restart continues after the highest export')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_export_numbering

    !> a publication for 3D holds only the stacks whose particles have been through an iteration
    !! (never particles just imported), renumbered from 1 with their particle ranges, stack indices
    !! and image indices, a 3D copy of the 2D particles, and the classes
    subroutine test_publication_holds_classified_stacks()
        integer, parameter :: NPTCLS_STK = 3
        type(sp_project) :: pool, pub
        type(string)     :: stk, mic
        integer          :: istk, iptcl, nstks
        write(*,'(A)') 'test_publication_holds_classified_stacks'
        ! three stacks of three particles (with nptcls_stk and indstk, as release 4 requires): stack 1
        ! classified, stack 2 just imported (never updated), stack 3 with one particle classified
        call pool%os_mic%new(3, is_ptcl=.false.)
        call pool%os_stk%new(3, is_ptcl=.false.)
        call pool%os_ptcl2D%new(3 * NPTCLS_STK, is_ptcl=.true.)
        call pool%os_cls2D%new(4, is_ptcl=.false.)
        do istk = 1,3
            call pool%os_mic%set(istk, 'intg', 'mic_'//int2str_pad(istk, 1)//'.mrc')
            call pool%os_stk%set(istk, 'stk',   'stk_'//int2str_pad(istk, 1)//'.mrcs')
            call pool%os_stk%set(istk, 'fromp', (istk - 1) * NPTCLS_STK + 1)
            call pool%os_stk%set(istk, 'top',   istk * NPTCLS_STK)
            call pool%os_stk%set(istk, 'nptcls_stk', NPTCLS_STK)
            do iptcl = (istk - 1) * NPTCLS_STK + 1,istk * NPTCLS_STK
                call pool%os_ptcl2D%set_stkind(iptcl, istk)
                call pool%os_ptcl2D%set(iptcl, 'indstk', iptcl - (istk - 1) * NPTCLS_STK)
                call pool%os_ptcl2D%set_state(iptcl, 1)
                call pool%os_ptcl2D%set(iptcl, 'x', real(iptcl))
            enddo
        enddo
        do iptcl = 1,NPTCLS_STK
            call pool%os_ptcl2D%set(iptcl, 'updatecnt', 1)
            call pool%os_ptcl2D%set_class(iptcl, 2)
        enddo
        call pool%os_ptcl2D%set(8, 'updatecnt', 1)
        call build_pool_publication(pool, pub, nstks)
        call assert_int(2, nstks,                        'the classified stacks are published')
        call assert_int(4, pool_publication_nselected(pool), 'the particles it selects are counted')
        call assert_int(pub%os_ptcl2D%count_state_gt_zero(), pool_publication_nselected(pool), 'as the publication selects them')
        call assert_int(2, pub%os_stk%get_noris(),       'stacks 1 and 3')
        call assert_int(2, pub%os_mic%get_noris(),       'with their micrographs')
        call assert_int(6, pub%os_ptcl2D%get_noris(),    'and particles; the just-imported stack is left out')
        call assert_int(6, pub%os_ptcl3D%get_noris(),    'the 3D particles are a copy of the 2D')
        stk = pub%os_stk%get_str(2, 'stk')
        mic = pub%os_mic%get_str(2, 'intg')
        call assert_char('stk_3.mrcs', stk%to_char(), 'the second published stack is stack 3')
        call assert_char('mic_3.mrc',  mic%to_char(), 'with its micrograph')
        call assert_int(4, pub%os_stk%get_fromp(2),      'its particle range is renumbered')
        call assert_int(6, pub%os_stk%get_top(2),        'to the published rows')
        call assert_int(2, pub%os_ptcl2D%get_int(4, 'stkind'), 'its particles point at it')
        call assert_int(1, pub%os_ptcl2D%get_int(4, 'indstk'), 'with their image index in the stack')
        call assert_real(7., pub%os_ptcl2D%get(4, 'x'), 1.e-4, 'and their parameters')
        call assert_int(2, pub%os_ptcl2D%get_class(1),   'the classes of classified particles')
        call assert_int(1, pub%os_ptcl2D%get_state(5),   'a classified particle of a published stack stays selected')
        call assert_int(0, pub%os_ptcl2D%get_state(4),   'one never updated is published deselected')
        call assert_int(0, pub%os_ptcl2D%get_state(6),   'like every never-updated particle of the stack')
        call assert_int(0, pub%os_ptcl3D%get_state(4),   'in both particle segments')
        call assert_int(4, pub%os_cls2D%get_noris(),     'and the pool''s class table')
        call assert_int(0, pub%os_ptcl3D%get_class(1),   'the 3D particles have no 2D class')
        call assert_real(1., pub%os_ptcl3D%get(1, 'x'), 1.e-4, 'and keep their shifts')
        call assert_int(0, pub%os_ptcl3D%get_int(1, 'updatecnt'), 'and no update count')
        call assert_int(0, pub%os_optics%get_noris(),    'without optics maps, no optics table')
        call pub%kill
        ! with an optics map: its groups and optics table
        call simple_mkdir('optics')
        call make_optics_map('optics', pool)
        call build_pool_publication(pool, pub, nstks, string('optics'))
        call assert_int(1, pub%os_optics%get_noris(),    'the publication carries the map''s optics table')
        call assert_int(7, pub%os_mic%get_int(1, 'ogid'), 'and its groups')
        call assert_int(7, pub%os_ptcl3D%get_int(1, 'ogid'), 'down to the particles')
        call pub%kill
        ! nothing classified yet: nothing to publish
        do iptcl = 1,3 * NPTCLS_STK
            call pool%os_ptcl2D%set(iptcl, 'updatecnt', 0)
        enddo
        call build_pool_publication(pool, pub, nstks)
        call assert_int(0, nstks,                        'a pool never classified publishes nothing')
        call pub%kill
        call pool%kill
    end subroutine test_publication_holds_classified_stacks

    !> the first publication from the sieve's own 2D (sieve_ini3D): the sets' stacks and particles
    !! with their class labels offset per set and the sieve's selection, the class tables
    !! concatenated with their states and the populations of the selected particles; and the sets'
    !! class averages, halves and FRCs concatenated in set order, refused when a half is missing or
    !! the boxes differ
    subroutine test_sieve_publication()
        integer, parameter :: BOX = 8
        real,    parameter :: SMPD = 2.0
        type(sp_project) :: sets(3), pub
        type(class_frcs) :: frcs_chk
        type(string)     :: cwd_saved, root, stks(3), frcs(3), stk_wide, frcs_wide
        integer          :: nfail0, nstks, ncls, ldim(3), n
        write(*,'(A)') 'test_sieve_publication'
        nfail0 = tests_failed
        call enter_fixture('p2_sieve_publication', cwd_saved, root)
        ! set A: a stack of three particles in classes 1, 2 and 2, its class 2 rejected by the sieve
        ! with its particles; set B: two particles, the second in a class outside its table; set C:
        ! the sieve's empty final set
        call make_sieve_set(sets(1), 'A', [1, 2, 2], [1, 0, 0], [1, 0])
        call make_sieve_set(sets(2), 'B', [1, 5],    [1, 1],    [1])
        call build_sieve_publication(sets, pub, nstks, ncls)
        call assert_int(2, nstks,                            'the sets'' stacks are published')
        call assert_int(3, ncls,                             'with their classes')
        call assert_int(5, pub%os_ptcl2D%get_noris(),        'and their particles')
        call assert_int(2, pub%os_ptcl2D%get_class(2),       'set A''s classes keep their labels')
        call assert_int(3, pub%os_ptcl2D%get_class(4),       'set B''s are offset by set A''s classes')
        call assert_int(0, pub%os_ptcl2D%get_state(2),       'the sieve''s selection is kept')
        call assert_int(0, pub%os_ptcl2D%get_state(5),       'a particle outside its set''s classes is deselected')
        call assert_int(4, pub%os_stk%get_fromp(2),          'the second stack''s range follows the first''s')
        call assert_int(1, pub%os_ptcl2D%get_int(4, 'indstk'), 'each particle with its image index in its stack')
        call assert_int(0, pub%os_cls2D%get_state(2),        'the sieve''s class states are kept')
        call assert_int(1, pub%os_cls2D%get_int(1, 'pop'),   'a population counts the selected particles')
        call assert_int(0, pub%os_cls2D%get_int(2, 'pop'),   'none in a rejected class')
        call assert_int(1, pub%os_cls2D%get_int(3, 'pop'),   'set B''s class counts its own')
        call assert_int(0, pub%os_ptcl3D%get_class(1),       'the 3D particles carry no 2D clustering')
        ! the class averages: set A's two, set B's one, with halves and FRCs
        call write_class_files('setA', 2, BOX, stks(1), frcs(1))
        call write_class_files('setB', 1, BOX, stks(2), frcs(2))
        stks(3) = ''
        frcs(3) = ''
        call assert_true(combine_sieve_classes(stks, frcs, [2, 1, 0], SMPD, string('pub_cavgs.mrcs'),&
            &string('pub_frcs.bin')), 'the class averages are combined')
        call find_ldim_nptcls(string('pub_cavgs.mrcs'), ldim, n)
        call assert_int(3, n, 'every class in set order')
        call find_ldim_nptcls(string('pub_cavgs_odd.mrcs'), ldim, n)
        call assert_int(3, n, 'with its halves')
        call frcs_chk%read(string('pub_frcs.bin'))
        call assert_int(3, frcs_chk%get_ncls(), 'and its FRCs')
        call frcs_chk%kill
        ! a set of another box
        call write_class_files('setW', 1, BOX + 2, stk_wide, frcs_wide)
        call assert_false(combine_sieve_classes([stks(1), stk_wide], [frcs(1), frcs_wide], [2, 1], SMPD,&
            &string('wide_cavgs.mrcs'), string('wide_frcs.bin')), 'sets of other boxes are refused')
        ! a set whose odd half is missing
        call del_file(add2fbody(stks(2), string('.mrc'), string('_odd')))
        call assert_false(combine_sieve_classes(stks, frcs, [2, 1, 0], SMPD, string('next_cavgs.mrcs'),&
            &string('next_frcs.bin')), 'a missing half is refused')
        call assert_false(file_exists(string('next_cavgs.mrcs')), 'and nothing is written')
        call pub%kill
        call leave_fixture(cwd_saved, root, nfail0)

    contains

        ! a sieve set: one micrograph and stack of particles with @p classes and @p states, and a
        ! class table of @p cls_states
        subroutine make_sieve_set( set, name, classes, states, cls_states )
            type(sp_project), intent(inout) :: set
            character(len=*), intent(in)    :: name
            integer,          intent(in)    :: classes(:), states(:), cls_states(:)
            integer :: iptcl, n, icls
            n = size(classes)
            call set%os_mic%new(1, is_ptcl=.false.)
            call set%os_stk%new(1, is_ptcl=.false.)
            call set%os_mic%set(1, 'intg', 'mic_'//name//'.mrc')
            call set%os_stk%set(1, 'stk',  'stk_'//name//'.mrcs')
            call set%os_stk%set(1, 'fromp', 1)
            call set%os_stk%set(1, 'top',   n)
            call set%os_stk%set(1, 'nptcls_stk', n)
            call set%os_ptcl2D%new(n, is_ptcl=.true.)
            do iptcl = 1,n
                call set%os_ptcl2D%set_stkind(iptcl, 1)
                call set%os_ptcl2D%set(iptcl, 'indstk', iptcl)
                call set%os_ptcl2D%set_class(iptcl, classes(iptcl))
                call set%os_ptcl2D%set_state(iptcl, states(iptcl))
            enddo
            call set%os_cls2D%new(size(cls_states), is_ptcl=.false.)
            do icls = 1,size(cls_states)
                call set%os_cls2D%set_state(icls, cls_states(icls))
            enddo
        end subroutine make_sieve_set

        ! @p n class averages of @p box pixels, their halves and their FRCs
        subroutine write_class_files( fbody, n, box, stk, frcs )
            character(len=*), intent(in)    :: fbody
            integer,          intent(in)    :: n, box
            type(string),     intent(inout) :: stk, frcs
            type(image)      :: img
            type(class_frcs) :: clsfrcs
            integer :: i
            stk  = fbody//'.mrc'
            frcs = fbody//'_frcs.bin'
            call img%new([box, box, 1], SMPD)
            do i = 1,n
                call img%write(stk, i)
                call img%write(string(fbody//'_even.mrc'), i)
                call img%write(string(fbody//'_odd.mrc'),  i)
            enddo
            call img%kill
            call clsfrcs%new(n, box, SMPD)
            call clsfrcs%write(frcs)
            call clsfrcs%kill
        end subroutine write_class_files

    end subroutine test_sieve_publication

    !> the sets present are taken in the order the sieve handed them off, by the number ending
    !! their name (10 after 2), with the final set last, whatever order the folder lists them in
    subroutine test_sets_taken_in_order()
        class(stream_stage_pool2D), allocatable :: stage
        type(cmdline)             :: cline
        type(string)              :: cwd_saved, root, cwd, set_any, names(4)
        type(rec_iterator)        :: it
        type(chunk_rec)           :: crec
        character(len=*), parameter :: EXPECTED(4) = [character(len=24) :: '00002.simple', '00003.simple',&
            &'00010.simple', 'sieve_final_c3_f3.simple']
        integer                   :: nfail0, i
        allocate(stage)
        write(*,'(A)') 'test_sets_taken_in_order'
        nfail0 = tests_failed
        call enter_fixture('p2_stage_order', cwd_saved, root)
        call simple_getcwd(cwd)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call make_upstream()
        call write_empty_final_set(cwd//'/'//UPSTREAM//'/'//DIR_STREAM_COMPLETED//'sieve_final_c3_f3.simple')
        set_any = write_sieved_set(10, [2])
        set_any = write_sieved_set(2,  [2])
        set_any = write_sieved_set(3,  [2])
        call stage%attach_upstream()
        call stage%watch_sets()
        call assert_int(4, stage%setslist%size(), 'every set present is recorded')
        if( stage%setslist%size() == 4 )then
            it = stage%setslist%begin()
            do i = 1,4
                call it%get(crec)
                names(i) = basename(crec%projfile)
                call it%next()
            enddo
            do i = 1,4
                call assert_char(trim(EXPECTED(i)), names(i)%to_char(), 'set '//int2str_pad(i, 1)//' in hand-off order')
            enddo
        endif
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_sets_taken_in_order

    !> the stage waits for the sieve's completed folder, then takes each set once; the first gives
    !! the sieve's mask diameter
    subroutine test_attach_and_watch()
        class(stream_stage_pool2D), allocatable :: stage
        type(cmdline)             :: cline
        type(string)              :: cwd_saved, root, set1
        integer                   :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_attach_and_watch'
        nfail0 = tests_failed
        call enter_fixture('p2_stage_watch', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%attach_upstream()
        call assert_false(stage%l_attached,      'no sieve folder: not attached')
        call assert_true(stage%l_waiting_logged, 'the wait is logged')
        call make_upstream()
        set1 = write_sieved_set(1, [3, 2])
        set1 = write_sieved_set(2, [4])
        call stage%attach_upstream()
        call assert_true(stage%l_attached, 'attached once the folder exists')
        call stage%watch_sets()
        call assert_int(2, stage%setslist%size(),                    'one record per set')
        call assert_true(all(.not. stage%setslist%get_included_flags()), 'none in the pool yet')
        call assert_true(stage%l_mskdiam_read,                         'the sieve''s mask diameter is read')
        call assert_real(MSKDIAM_SIEVE, stage%final_mskdiam, 1.e-4,    'from the set''s class averages')
        call stage%watch_sets()
        call assert_int(2, stage%setslist%size(), 'a set is taken once')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_attach_and_watch

    !> the sets go to the pool in order, each once, and a later set after the others; the stage
    !! counts the selected particles and the micrographs (the rows: simple_stream_pool2D_tester)
    subroutine test_transfer_sets()
        class(stream_stage_pool2D), allocatable :: stage
        type(cmdline)             :: cline
        type(string)              :: cwd_saved, root, set_file
        integer                   :: nfail0, nimported
        allocate(stage)
        write(*,'(A)') 'test_transfer_sets'
        nfail0 = tests_failed
        call enter_fixture('p2_stage_transfer', cwd_saved, root)
        call make_upstream()
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%attach_upstream()
        ! one watch per set, so the records keep the order the sets are written in
        set_file = write_sieved_set(1, [3, 2], nrejected=1)
        call stage%watch_sets()
        set_file = write_sieved_set(2, [4])
        call stage%watch_sets()
        call stage%transfer_sets(nimported)
        call assert_int(2, nimported,                    'both sets are imported')
        call assert_int(8, stage%nptcls_glob,            'the selected particles are counted')
        call assert_int(8, stage%nptcls_glob_state_1,    'so are those of the pool with state > 0')
        call assert_int(3, stage%nmics,                  'and the micrographs')
        call assert_int(3, stage%state_1_particle_rate,  'particles per micrograph, rounded up')
        call assert_true(all(stage%setslist%get_included_flags()), 'both sets are included')
        call stage%transfer_sets(nimported)
        call assert_int(0, nimported, 'a set is imported once')
        set_file = write_sieved_set(3, [2])
        call stage%watch_sets()
        call stage%transfer_sets(nimported)
        call assert_int(1,  nimported,                   'a later set is imported')
        call assert_int(4,  stage%nmics,                 'the pool grows by its micrograph')
        call assert_int(10, stage%nptcls_glob,           'and its particles')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_transfer_sets

    !> stepwise=yes takes only enough sets to reach the threshold; the rest come with the next import
    subroutine test_transfer_stepwise()
        class(stream_stage_pool2D), allocatable :: stage
        type(cmdline)             :: cline
        type(string)              :: cwd_saved, root, set_file
        integer                   :: nfail0, nimported
        allocate(stage)
        write(*,'(A)') 'test_transfer_stepwise'
        nfail0 = tests_failed
        call enter_fixture('p2_stage_stepwise', cwd_saved, root)
        call make_upstream()
        call set_test_cline(cline)
        call cline%set('stepwise', 'yes')
        call make_test_stage(stage, cline)
        stage%nptcls_threshold = 20
        call stage%attach_upstream()
        ! one watch per set, so the records keep the order the sets are written in
        set_file = write_sieved_set(1, [15])
        call stage%watch_sets()
        set_file = write_sieved_set(2, [10])
        call stage%watch_sets()
        set_file = write_sieved_set(3, [10])
        call stage%watch_sets()
        call stage%transfer_sets(nimported)
        call assert_int(2,  nimported,                 'the sets up to the threshold are imported')
        call assert_int(25, stage%nptcls_glob,         'their particles')
        call stage%transfer_sets(nimported)
        call assert_int(1,  nimported,                 'the deferred set comes with the next import')
        call assert_int(35, stage%nptcls_glob,         'after the others')
        ! past the threshold, an import still takes the sets its own particles need
        set_file = write_sieved_set(4, [10])
        call stage%watch_sets()
        set_file = write_sieved_set(5, [10])
        call stage%watch_sets()
        set_file = write_sieved_set(6, [10])
        call stage%watch_sets()
        call stage%transfer_sets(nimported)
        call assert_int(2,  nimported,                 'only this import''s particles count against the threshold')
        call assert_int(55, stage%nptcls_glob,         'two more sets')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_transfer_stepwise

    subroutine test_sieve_final_set()
        class(stream_stage_pool2D), allocatable :: stage
        type(cmdline)             :: cline
        type(string)              :: cwd_saved, root, set_file
        integer                   :: nfail0, nimported, nmics0
        allocate(stage)
        write(*,'(A)') 'test_sieve_final_set'
        nfail0 = tests_failed
        call enter_fixture('p2_stage_final', cwd_saved, root)
        call make_upstream()
        set_file = write_sieved_set(1, [3])
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%attach_upstream()
        call stage%watch_sets()
        call stage%transfer_sets(nimported)
        call assert_false(stage%l_sieve_final, 'an ordinary set')
        set_file = write_sieved_set(2, [3], l_final=.true.)
        call stage%watch_sets()
        call stage%transfer_sets(nimported)
        call assert_true(stage%l_sieve_final, 'the sieve''s final set is noted')
        ! a later set with particles that is not final: the sieve had more after all
        set_file = write_sieved_set(3, [3])
        call stage%watch_sets()
        call stage%transfer_sets(nimported)
        call assert_false(stage%l_sieve_final, 'a later ordinary set takes the final note back')
        ! the sieve's empty final set ends the intake and transfers nothing
        call write_empty_final_set(string(UPSTREAM//'/'//DIR_STREAM_COMPLETED//'sieve_final_c1_f1'//METADATA_EXT))
        nmics0 = stage%nmics
        call stage%watch_sets()
        call stage%transfer_sets(nimported)
        call assert_true(stage%l_sieve_final,          'an empty final set is noted')
        call assert_int(0, nimported,                  'and imports nothing')
        call assert_int(nmics0, stage%nmics,           'the pool is unchanged')
        call assert_int(0, count(.not. stage%setslist%get_included_flags()), 'the empty set is taken')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_sieve_final_set

    ! the sieve's empty final set (ptcl_sieve: hand_off_final_set)
    subroutine write_empty_final_set( fname )
        class(string), intent(in) :: fname
        type(sp_project) :: final_set
        call final_set%os_out%new(1, is_ptcl=.false.)
        call final_set%os_out%set(1, 'sieve_final', 'yes')
        call final_set%write(fname)
        call final_set%kill
    end subroutine write_empty_final_set

    !> the pause rule, the particle targets, the final run and the default mask diameter
    subroutine test_pause_rules()
        class(stream_stage_pool2D), allocatable :: stage
        allocate(stage)
        write(*,'(A)') 'test_pause_rules'
        call assert_int(0,   stage%pause_rate_factor(0,  -1), 'no pause before the first iteration')
        call assert_int(0,   stage%pause_rate_factor(1,  -1), 'nor after it')
        call assert_int(0,   stage%pause_rate_factor(2,   1), 'iterations 2-20: one iteration without import is no pause')
        call assert_int(50,  stage%pause_rate_factor(2,   0), 'two are')
        call assert_int(50,  stage%pause_rate_factor(20, 18), 'up to iteration 20')
        call assert_int(0,   stage%pause_rate_factor(21, 21), 'after iteration 20: an import this iteration is no pause')
        call assert_int(500, stage%pause_rate_factor(21, 20), 'one iteration without is')
        call assert_int(300, stage%target_nptcls(100, 10, 3,  50), 'at least 20 particles per class')
        call assert_int(600, stage%target_nptcls(100, 10, 10, 50), 'or the micrographs'' worth')
        call assert_true(stage%runs_to_final(24, .true.),  'the final set runs the pool uninterrupted')
        call assert_false(stage%runs_to_final(25, .true.), 'up to iteration 25')
        call assert_false(stage%runs_to_final(5, .false.), 'only with the final set')
        call assert_true(stage%publishes_final(25, .true.),  'the publication after iteration 25 with the final set is final')
        call assert_true(stage%publishes_final(31, .true.),  'as is one after a later iteration')
        call assert_false(stage%publishes_final(24, .true.), 'not before iteration 25')
        call assert_false(stage%publishes_final(30, .false.), 'nor without the final set')
        call assert_true(stage%exports_after(10, .true., 0),   'a fresh pool publishes after iteration 10')
        call assert_false(stage%exports_after(10, .false., 0), 'a pool with publications on disk does not')
        call assert_false(stage%exports_after(9, .true., NPTCLS_FIRST3D - 1), 'nor before it, with fewer selected than the first 3D takes')
        call assert_true(stage%exports_after(4, .true., NPTCLS_FIRST3D),  'but before it, once as many are selected')
        call assert_false(stage%exports_after(4, .false., NPTCLS_FIRST3D), 'only when it has published nothing yet')
        call assert_false(stage%exports_after(11, .true., NPTCLS_FIRST3D), 'nothing between 10 and 25')
        call assert_false(stage%exports_after(24, .true., 0),  'up to iteration 24')
        call assert_true(stage%exports_after(25, .false., 0),  'from iteration 25 every pool publishes')
        call assert_true(stage%exports_after(31, .true., 0),   'after each iteration')
        call assert_real(100., stage%default_mskdiam(64, 2.0), 1.e-4, 'the default mask diameter is in Angstroms')
    end subroutine test_pause_rules

    !> a new mask diameter from the GUI is taken, resumes a paused pool and replaces the sieve's;
    !! the same value again changes nothing; a snapshot request waits for the pool
    subroutine test_gui_mskdiam_update()
        class(stream_stage_pool2D), allocatable        :: stage
        type(cmdline)                    :: cline
        type(stream_pipe)                :: writer
        type(gui_metadata_stream_update) :: update
        type(string)                     :: cwd_saved, root
        integer(c_int)                   :: fds(2)
        integer                          :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_gui_mskdiam_update'
        nfail0 = tests_failed
        call enter_fixture('p2_stage_gui_update', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call open_loopback(fds)
        call stage%init_gui(int(fds(1)), -1)
        call writer%new(-1, int(fds(2)), max_metadata_size(), 'test writer')
        stage%l_pause       = .true.
        stage%final_mskdiam = MSKDIAM_SIEVE
        call update%new(GUI_METADATA_STREAM_UPDATE_TYPE)
        call update%set_mskdiam2D_update(150.)
        call writer%send_meta(update)
        call stage%apply_gui_updates()
        call assert_real(150., stage%mskdiam, 1.e-4, 'the mask diameter is updated')
        call assert_false(stage%l_pause,                    'a paused pool resumes')
        call assert_real(0., stage%final_mskdiam, 1.e-4,    'the sieve''s pending mask diameter is dropped')
        stage%l_pause = .true.
        call writer%send_meta(update)
        call stage%apply_gui_updates()
        call assert_true(stage%l_pause, 'the same mask diameter again changes nothing')
        call update%kill
        call update%new(GUI_METADATA_STREAM_UPDATE_TYPE)
        call update%set_snapshot2D_update(1, 3, [1, 2], string('snap.simple'))
        call writer%send_meta(update)
        call stage%apply_gui_updates()
        call assert_int(1, stage%last_snapshot_id, 'a snapshot before the pool runs is answered at once')
        call assert_int(0, stage%snapshot_nptcls,  'as not written')
        call update%kill
        call writer%kill
        call stage%kill
        call close_loopback(fds)
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_gui_mskdiam_update

    !> from iteration 10 the sieve's mask diameter replaces the pool's, unless one was given on the
    !! command line
    subroutine test_given_mskdiam_kept()
        class(stream_stage_pool2D), allocatable :: stage
        type(cmdline) :: cline
        type(string)  :: cwd_saved, root
        integer       :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_given_mskdiam_kept'
        nfail0 = tests_failed
        call enter_fixture('p2_stage_given_mskdiam', cwd_saved, root)
        call set_test_cline(cline)
        call cline%set('mskdiam', 150.)
        call make_test_stage(stage, cline)
        call assert_true(stage%l_mskdiam_given, 'a mask diameter on the command line is given')
        stage%final_mskdiam = MSKDIAM_SIEVE
        call stage%apply_final_mskdiam(10)
        call assert_real(150., stage%mskdiam, 1.e-4, 'it is kept at iteration 10')
        call assert_real(0., stage%final_mskdiam, 1.e-4, 'and the sieve''s is dropped')
        call stage%kill
        call cline%kill
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call assert_false(stage%l_mskdiam_given, 'none given')
        stage%final_mskdiam = MSKDIAM_SIEVE
        call stage%apply_final_mskdiam(10)
        call assert_real(MSKDIAM_SIEVE, stage%mskdiam, 1.e-4, 'the sieve''s applies at iteration 10')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_given_mskdiam_kept

    !> one status message per call, with the stage name, the particles imported and the mask
    subroutine test_send_status()
        class(stream_stage_pool2D), allocatable        :: stage
        type(cmdline)                    :: cline
        type(stream_pipe)                :: reader
        type(gui_metadata_stream_pool2D) :: status
        character(len=:), allocatable    :: buffer
        type(string)                     :: cwd_saved, root, stage_name
        integer(c_int)                   :: fds(2)
        integer                          :: nfail0, meta_type, iter, nimported, naccepted, nrejected, tlast, mskdiam
        real                             :: mskscale, res
        logical                          :: l_assigned, l_user_input
        allocate(stage)
        write(*,'(A)') 'test_send_status'
        nfail0 = tests_failed
        call enter_fixture('p2_stage_status', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call open_loopback(fds)
        call stage%init_gui(-1, int(fds(2)))
        call reader%new(int(fds(1)), -1, max_metadata_size(), 'test reader')
        stage%nptcls_glob    = 42
        stage%mskdiam = 160.
        call stage%send_status(string('waiting for sieved particles'))
        call assert_true(reader%receive(buffer), 'a status message is sent')
        if( allocated(buffer) )then
            meta_type = transfer(buffer, meta_type)
            call assert_int(GUI_METADATA_STREAM_POOL2D_TYPE, meta_type, 'it is a pool 2D status')
            if( meta_type == GUI_METADATA_STREAM_POOL2D_TYPE )then
                status     = transfer(buffer, status)
                l_assigned = status%get(stage_name, iter, nimported, naccepted, nrejected, tlast, l_user_input,&
                    &mskdiam, mskscale, res)
                call assert_char('waiting for sieved particles', stage_name%to_char(), 'the stage name')
                call assert_int(42,  nimported, 'the particles imported')
                call assert_int(160, mskdiam,   'the mask diameter')
            endif
        endif
        call assert_false(reader%receive(buffer), 'one message per call')
        call reader%kill
        call stage%kill
        call close_loopback(fds)
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_send_status

    !> a snapshot is announced with its project and particle count, then its selected classes
    !! follow as tiles; a snapshot that could not be written is announced with no particles, no
    !! file and no tiles
    subroutine test_send_snapshot()
        class(stream_stage_pool2D), allocatable                 :: stage
        type(cmdline)                             :: cline
        type(stream_pipe)                         :: reader
        type(gui_metadata_stream_snapshot) :: snapshot
        character(len=:), allocatable             :: buffer
        type(string)                              :: cwd_saved, root, fname
        integer(c_int)                            :: fds(2)
        integer                                   :: nfail0, meta_type, id, nptcls, stime, ncavgs
        logical                                   :: l_assigned
        allocate(stage)
        write(*,'(A)') 'test_send_snapshot'
        nfail0 = tests_failed
        call enter_fixture('p2_stage_snapshot', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call open_loopback(fds)
        call stage%init_gui(-1, int(fds(2)))
        call reader%new(int(fds(1)), -1, max_metadata_size(), 'test reader')
        stage%last_snapshot_id  = 3
        stage%snapshot_dir      = '/data/snapshots/snap'
        stage%snapshot_filename = 'snap.simple'
        stage%snapshot_jpeg     = '/data/snapshots/snap/thumb2D.jpeg'
        stage%snapshot_mrc      = '/data/snapshots/snap/cavgs.mrc'
        stage%snapshot_idx      = [2, 5, 7]
        stage%snapshot_pop      = [10, 20, 30]
        stage%snapshot_res      = [8., 9., 10.]
        stage%snapshot_ntilesx  = 2
        stage%snapshot_ntilesy  = 2
        stage%snapshot_nptcls   = 42
        call stage%send_snapshot()
        call assert_true(reader%receive(buffer), 'the snapshot is announced')
        if( allocated(buffer) )then
            meta_type = transfer(buffer, meta_type)
            call assert_int(GUI_METADATA_STREAM_POOL2D_SNAPSHOT_TYPE, meta_type, 'as a pool 2D snapshot')
            if( meta_type == GUI_METADATA_STREAM_POOL2D_SNAPSHOT_TYPE )then
                snapshot   = transfer(buffer, snapshot)
                l_assigned = snapshot%get(id, fname, nptcls, stime)
                call assert_int(3, id, 'with its id')
                call assert_char('/data/snapshots/snap/snap.simple', fname%to_char(), 'and its project')
                call assert_int(42, nptcls, 'and its particle count')
            endif
        endif
        ncavgs = 0
        do while( reader%receive(buffer) )
            meta_type = transfer(buffer, meta_type)
            if( meta_type == GUI_METADATA_STREAM_POOL2D_SNAPSHOT_CLS2D_TYPE ) ncavgs = ncavgs + 1
        enddo
        call assert_int(3, ncavgs, 'then one tile per selected class')
        ! a snapshot of an iteration the pool no longer keeps
        stage%last_snapshot_id = 4
        stage%snapshot_nptcls  = 0
        stage%snapshot_idx     = [integer ::]
        call stage%send_snapshot()
        call assert_true(reader%receive(buffer), 'a snapshot not written is announced')
        if( allocated(buffer) )then
            meta_type = transfer(buffer, meta_type)
            if( meta_type == GUI_METADATA_STREAM_POOL2D_SNAPSHOT_TYPE )then
                snapshot   = transfer(buffer, snapshot)
                l_assigned = snapshot%get(id, fname, nptcls, stime)
                call assert_int(4, id,     'with its id')
                call assert_int(0, nptcls, 'with no particles')
                call assert_int(0, fname%strlen_trim(), 'and no file')
            endif
        endif
        call assert_false(reader%receive(buffer), 'and no tiles')
        call reader%kill
        call stage%kill
        call close_loopback(fds)
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_send_snapshot

    !> the public loop waits for the sieve's folder, then attaches; no pool until a set arrives
    subroutine test_iterate_waits()
        class(stream_stage_pool2D), allocatable :: stage
        type(cmdline)             :: cline
        type(string)              :: cwd_saved, root
        integer                   :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_iterate_waits'
        nfail0 = tests_failed
        call enter_fixture('p2_stage_iterate', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%iterate()
        call assert_false(stage%l_attached, 'pass 1: no sieve folder, not attached')
        call make_upstream()
        call stage%iterate()
        call assert_true(stage%l_attached,      'pass 2: attached')
        call assert_false(stage%l_pool_started, 'pass 2: no set yet, no pool')
        call stage%kill
        call stage%kill ! idempotence
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_iterate_waits

    subroutine test_finished()
        class(stream_stage_pool2D), allocatable :: stage
        type(cmdline)             :: cline
        type(string)              :: cwd_saved, root
        integer                   :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_finished'
        nfail0 = tests_failed
        call enter_fixture('p2_stage_finished', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call assert_false(stage%finished(), 'a fresh stage is not finished')
        call simple_touch(TERM_STREAM)
        call assert_true(stage%finished(), 'the termination file finishes the stage')
        call del_file(TERM_STREAM)
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_finished

    ! ---- fixtures ------------------------------------------------------------

    ! the stage's command line under a program name outside every UI table, with a local queue
    ! system for the computing environment of the stage's project
    subroutine set_test_cline( cline )
        type(cmdline), intent(inout) :: cline
        call cline%set('prg',        'stream_pool2D_stage_tester')
        call cline%set('mkdir',      'no')
        call cline%set('projfile',   TEST_PROJFILE)
        call cline%set('outdir',     '')
        call cline%set('nthr',       1)
        call cline%set('nparts',     1)
        call cline%set('ncls',       10)
        call cline%set('dir_target', UPSTREAM)
        call cline%set('qsys_name',  'local')
    end subroutine set_test_cline

    ! a stage without a pipe, with no waits and a settle time that takes files written in the same
    ! second
    subroutine make_test_stage( stage, cline )
        class(stream_stage_pool2D), intent(inout) :: stage
        type(cmdline),             intent(inout) :: cline
        call stage%init_params(cline)
        call stage%init_gui(-1, -1)
        stage%settle_s = -1
        stage%wait_s   = 0
        stage%l_exists = .true.
    end subroutine make_test_stage

    ! the folders the sieve makes
    subroutine make_upstream()
        call simple_mkdir(UPSTREAM)
        call simple_mkdir(UPSTREAM//'/'//DIR_STREAM_COMPLETED)
    end subroutine make_upstream

    ! handed-off set number @p id: one micrograph and one stack per entry of @p nptcls_mic, its
    ! particles in class 3 with a shift, the first @p nrejected of them deselected, and the class
    ! averages' entry with the sieve's mask diameter (and sieve_final=yes when @p l_final); returns
    ! its absolute path
    function write_sieved_set( id, nptcls_mic, nrejected, l_final ) result( fname )
        integer,           intent(in) :: id, nptcls_mic(:)
        integer, optional, intent(in) :: nrejected
        logical, optional, intent(in) :: l_final
        type(string)     :: fname
        type(sp_project) :: proj
        type(string)     :: cwd
        integer          :: imic, fromp, iptcl, nptcls
        call simple_getcwd(cwd)
        fname  = cwd//'/'//UPSTREAM//'/'//DIR_STREAM_COMPLETED//int2str_pad(id, 5)//METADATA_EXT
        nptcls = sum(nptcls_mic)
        call proj%os_mic%new(size(nptcls_mic), is_ptcl=.false.)
        call proj%os_stk%new(size(nptcls_mic), is_ptcl=.false.)
        call proj%os_ptcl2D%new(nptcls, is_ptcl=.true.)
        fromp = 1
        do imic = 1,size(nptcls_mic)
            call proj%os_mic%set_state(imic, 1)
            call proj%os_mic%set(imic, 'nptcls',    nptcls_mic(imic))
            call proj%os_mic%set(imic, 'importind', (id - 1) * 10 + imic)
            call proj%os_stk%set(imic, 'nptcls',    nptcls_mic(imic))
            call proj%os_stk%set(imic, 'fromp',     fromp)
            call proj%os_stk%set(imic, 'top',       fromp + nptcls_mic(imic) - 1)
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
        call proj%os_out%new(1, is_ptcl=.false.)
        call proj%os_out%set(1, 'imgkind', 'cavg')
        call proj%os_out%set(1, 'stk',     'cavgs.mrc')
        call proj%os_out%set(1, 'nptcls',  10)
        call proj%os_out%set(1, 'mskdiam', MSKDIAM_SIEVE)
        if( present(l_final) )then
            if( l_final ) call proj%os_out%set(1, 'sieve_final', 'yes')
        endif
        call proj%update_projinfo(fname)
        call proj%write(fname)
        call proj%kill
    end function write_sieved_set

    ! optics map 1 in @p dir: every micrograph of @p proj (import index 1, 2, ...) in group 7, and
    ! an optics table of that one group
    subroutine make_optics_map( dir, proj )
        character(len=*), intent(in)    :: dir
        type(sp_project), intent(inout) :: proj
        type(sp_project) :: mapproj
        integer :: imic
        call mapproj%os_mic%new(proj%os_mic%get_noris(), is_ptcl=.false.)
        do imic = 1,proj%os_mic%get_noris()
            call mapproj%os_mic%set(imic, 'importind', imic)
            call mapproj%os_mic%set(imic, 'ogid',      7)
            call proj%os_mic%set(imic, 'importind', imic)
        enddo
        call mapproj%os_optics%new(1, is_ptcl=.false.)
        call mapproj%os_optics%set(1, 'ogid', 7)
        call mapproj%os_optics%set(1, 'pop',  proj%os_mic%get_noris())
        call publish_optics_map(mapproj, string(dir), 1, 5)
        call mapproj%kill
    end subroutine make_optics_map

    ! a pipe with a non-blocking read end
    subroutine open_loopback( fds )
        integer(c_int), intent(inout) :: fds(2)
        integer(c_int) :: rc, flags
        fds   = -1
        rc    = c_pipe(fds)
        call assert_int(0, int(rc), 'pipe() succeeds')
        flags = c_fcntl(fds(1), F_GETFL, 0_c_int)
        rc    = c_fcntl(fds(1), F_SETFL, ior(flags, O_NONBLOCK))
        call assert_int(0, int(rc), 'the read end is made non-blocking')
    end subroutine open_loopback

    subroutine close_loopback( fds )
        integer(c_int), intent(inout) :: fds(2)
        integer(c_int) :: rc
        rc = c_close(fds(1))
        call assert_int(0, int(rc), 'the read end closes')
        rc = c_close(fds(2))
        call assert_int(0, int(rc), 'the write end closes')
        fds = -1
    end subroutine close_loopback

end module simple_stream_stage_pool2D_tester
