!@descr: unit tests for the steps of the stream multistate 3D stage (simple_stream_stage_solve3D)
! Each test assembles a stage in a fresh fixture directory from init_params and init_gui: no queue
! environment, no waits, and a settle time of -1. The command line names no registered program
! and carries qsys_name=local for the stage project's computing environment. The upstream is a pool
! 2D directory whose completed folder holds exports. The class-average selection and the 3D jobs
! (which need a queue) are left to the high-level stream tests; the import is tested on sets made
! in memory, the class averages a publication brings on fixture stacks, and the volume messages
! on a fixture project with a volume and an FSC; a 3D snapshot is written from a fixture result
! project, on a request sent through the stage's update pipe.
module simple_stream_stage_solve3D_tester
use, intrinsic :: iso_c_binding, only: c_int
use unix,                        only: c_pipe, c_close, c_fcntl, F_GETFL, F_SETFL, O_NONBLOCK
use simple_test_utils
use simple_defs_fname,                                only: TERM_STREAM, METADATA_EXT
use simple_refine3D_fnames,                           only: refine3D_reprojs_fname
use simple_defs_stream,                               only: DIR_STREAM_COMPLETED
use simple_string,                                    only: string
use simple_string_utils,                              only: int2str, int2str_pad
use simple_fileio,                                    only: arr2file, del_file, file_exists, simple_getcwd, simple_touch
use simple_syslib,                                    only: simple_mkdir, dir_exists
use simple_math_ft,                                   only: get_resarr
use simple_cmdline,                                   only: cmdline
use simple_image,                                     only: image
use simple_sp_project,                                only: sp_project
use simple_rec_list,                                  only: rec_iterator, chunk_rec
use simple_gui_metadata_utils,                        only: max_metadata_size
use simple_gui_metadata_types,                        only: GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_TYPE, GUI_METADATA_VOL3D_TYPE,&
                                                           &GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_REPROJ_TYPE,&
                                                           &GUI_METADATA_STREAM_SOLVE3D_SNAPSHOT_TYPE, GUI_METADATA_STREAM_UPDATE_TYPE
use simple_gui_metadata_stream_solve3D_multistate, only: gui_metadata_stream_solve3D_multistate
use simple_gui_metadata_stream_snapshot,              only: gui_metadata_stream_snapshot
use simple_gui_metadata_stream_update,                only: gui_metadata_stream_update
use simple_gui_metadata_vol3D,                        only: gui_metadata_vol3D
use simple_stream_pipe,                               only: stream_pipe
use simple_stream_stage_solve3D,                   only: stream_stage_solve3D, PHASE_IMPORTING, PHASE_SOLVE3D,&
                                                           &PHASE_IDLE, PHASE_ADDON, JOB_NONE, JOB_SOLVE3D, JOB_ADDON
implicit none
private
public :: run_all_stream_stage_solve3D_tests

character(len=*), parameter :: TEST_PROJFILE = 'test_solve3D.simple'
character(len=*), parameter :: UPSTREAM      = 'pool2D' ! the pool 2D stage directory (dir_target)
integer,          parameter :: NSTATES       = 3
integer,          parameter :: VOL_BOX       = 16
real,             parameter :: VOL_SMPD      = 2.0

contains

    subroutine run_all_stream_stage_solve3D_tests()
        write(*,'(A)') '**** running all stream solve 3D stage tests ****'
        call test_init_params()
        call test_restart_removes_term_stream()
        call test_watch_order_and_mskdiam()
        call test_merge_publications()
        call test_take_cavgs()
        call test_first_set()
        call test_first_set_fallback()
        call test_mskdiam_from_each_publication()
        call test_rows_problem()
        call test_rules()
        call test_cohort()
        call test_retention()
        call test_send_status()
        call test_send_volumes()
        call test_snapshot3D()
        call test_iterate_waits()
        call test_finished()
    end subroutine run_all_stream_stage_solve3D_tests

    subroutine test_init_params()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)                 :: cline
        type(string)                  :: cwd_saved, root
        integer                       :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_init_params'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_init_params', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call assert_true(allocated(stage%params), 'initialization allocates the owned parameters')
        call assert_int(0, stage%spproj%os_mic%get_noris(), 'the project starts without micrographs')
        call assert_false(stage%l_restart,                  'a fresh run is no restart')
        call assert_int(PHASE_IMPORTING, stage%phase,       'the stage starts importing')
        call assert_int(NSTATES, size(stage%state_res),     'one resolution per state')
        call assert_int(0, size(stage%stk_names),           'no stack yet')
        call stage%kill
        call assert_false(allocated(stage%params), 'cleanup releases the owned parameters')
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_init_params

    !> a restart removes a leftover termination file, which would end the restarted stage at once
    subroutine test_restart_removes_term_stream()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)                 :: cline
        type(string)                  :: cwd_saved, root
        integer                       :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_restart_removes_term_stream'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_restart', cwd_saved, root)
        call simple_mkdir('previous')
        call simple_touch(TERM_STREAM)
        call set_test_cline(cline)
        call cline%set('outdir', 'previous')
        call make_test_stage(stage, cline)
        call assert_true(stage%l_restart,           'an existing output folder is a restart')
        call assert_false(file_exists(TERM_STREAM), 'the leftover termination file is removed')
        call assert_false(stage%finished(),         'so the restarted stage runs')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_restart_removes_term_stream

    !> exports are taken in export order whatever order they are listed in, each once
    subroutine test_watch_order_and_mskdiam()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)                 :: cline
        type(rec_iterator)            :: it
        type(chunk_rec)               :: crec
        type(string)                  :: cwd_saved, root
        integer                       :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_watch_order_and_mskdiam'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_watch', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%attach_upstream()
        call assert_false(stage%l_attached,      'no pool 2D folder: not attached')
        call assert_true(stage%l_waiting_logged, 'the wait is logged')
        call simple_mkdir(UPSTREAM)
        call simple_mkdir(UPSTREAM//'/'//DIR_STREAM_COMPLETED)
        call write_export(2, 170.)
        call write_export(1, 150.)
        call stage%attach_upstream()
        call assert_true(stage%l_attached, 'attached once the folder exists')
        call stage%watch_sets()
        call assert_int(2, stage%setslist%size(), 'one record per export')
        if( stage%setslist%size() == 2 )then
            it = stage%setslist%begin()
            call it%get(crec)
            call assert_true(crec%projfile%has_substr('00001'//METADATA_EXT), 'export 1 comes first')
            call it%next()
            call it%get(crec)
            call assert_true(crec%projfile%has_substr('00002'//METADATA_EXT), 'then export 2')
        endif
        call stage%watch_sets()
        call assert_int(2, stage%setslist%size(), 'an export is taken once')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_watch_order_and_mskdiam

    !> publications are merged into rows that only grow: a new stack is appended with both particle
    !! segments; a stack the stage holds is matched by name in any order, its particles take the
    !! publication's class and selection and keep their 3D parameters; a stack a publication lacks
    !! keeps its rows, deselected; the classes are the newest publication's
    subroutine test_merge_publications()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)                 :: cline
        type(sp_project)              :: set
        type(string)                  :: cwd_saved, root
        integer                       :: nfail0, i
        allocate(stage)
        write(*,'(A)') 'test_merge_publications'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_merge', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        ! publication 1: stacks A (3) and B (2), the first particle deselected, class 1
        call make_set(set, ['A', 'B'], [3, 2], 2, nrejected=1, icls=1)
        call stage%merge_publication(set, 1)
        call set%kill
        call assert_int(2, stage%spproj%os_mic%get_noris(),    'publication 1: its micrographs')
        call assert_int(5, stage%spproj%os_ptcl3D%get_noris(), 'its 3D particles')
        call assert_int(5, stage%spproj%os_ptcl2D%get_noris(), 'and 2D particles')
        call assert_int(4, stage%nptcls_selected,              'the selected particles are counted')
        ! a 3D result for the first particle, which later publications must keep, and a run's
        ! multistate label on the second
        call stage%spproj%os_ptcl3D%set(1, 'e1', 33.)
        call stage%spproj%os_ptcl3D%set_state(2, 2)
        ! publication 2, as a restarted pool lists it: B, A, then the new stack C (4), class 2
        call make_set(set, ['B', 'A', 'C'], [2, 3, 4], 5, icls=2)
        call stage%merge_publication(set, 2)
        call set%kill
        call assert_int(3, stage%spproj%os_mic%get_noris(),           'publication 2: one stack appended')
        call assert_int(9, stage%spproj%os_ptcl3D%get_noris(),        'with its particles; no row renumbered')
        call assert_int(6, stage%spproj%os_stk%get_fromp(3),          'the new stack follows the rows held')
        call assert_int(3, stage%spproj%os_ptcl3D%get_int(7, 'stkind'), 'its particles point at it')
        call assert_real(3., stage%spproj%os_ptcl3D%get(7, 'x'), 1.e-4, 'and carry its parameters')
        call assert_real(1., stage%spproj%os_ptcl3D%get(1, 'x'), 1.e-4, 'row 1 is still stack A''s first particle')
        call assert_int(2, stage%spproj%os_ptcl2D%get_class(1),       'which takes the new class')
        call assert_int(1, stage%spproj%os_ptcl3D%get_state(1),       'and the new selection')
        call assert_int(2, stage%spproj%os_ptcl3D%get_state(2),       'a selected particle keeps its 3D state label')
        call assert_real(33., stage%spproj%os_ptcl3D%get(1, 'e1'), 1.e-4, 'and keeps its 3D parameters')
        call assert_int(1, stage%spproj%os_ptcl2D%get_int(1, 'stkind'), 'and its stack')
        call assert_int(5, stage%spproj%os_cls2D%get_noris(),         'the classes are the newest publication''s')
        call assert_int(9, stage%nptcls_selected,                     'every particle selected')
        ! publication 3 lacks stack B (rows 4 and 5)
        call make_set(set, ['A', 'C'], [3, 4], 5, icls=2)
        call stage%merge_publication(set, 3)
        call set%kill
        call assert_int(9, stage%spproj%os_ptcl3D%get_noris(), 'publication 3: the rows are kept')
        call assert_int(0, stage%spproj%os_ptcl3D%get_state(4), 'a stack it lacks is deselected')
        call assert_int(0, stage%spproj%os_ptcl2D%get_state(5), 'in both particle segments')
        call assert_int(7, stage%nptcls_selected,               'and no longer counted')
        ! publication 4 lists stack A's images in reverse order: each updates its own row
        call make_set(set, ['A'], [3], 13, icls=1)
        do i = 1,3
            call set%os_ptcl2D%set(i, 'indstk', 4 - i)
            call set%os_ptcl2D%set_class(i, 14 - i)
        enddo
        call stage%merge_publication(set, 4)
        call set%kill
        call assert_int(11, stage%spproj%os_ptcl2D%get_class(1), 'a particle is matched by its image index, not its row')
        call assert_int(13, stage%spproj%os_ptcl2D%get_class(3), 'whatever order the publication lists them in')
        ! a publication with the optics map's table brings it
        call make_set(set, ['A'], [3], 13, icls=1)
        call set%os_optics%new(2, is_ptcl=.false.)
        call set%os_optics%set(1, 'ogid', 4)
        call set%os_optics%set(2, 'ogid', 9)
        call stage%merge_publication(set, 5)
        call set%kill
        call assert_int(2, stage%spproj%os_optics%get_noris(), 'the stage takes the publication''s optics table')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_merge_publications

    !> a publication's class averages and FRCs are copied into its quality folder and replace the
    !! stage's earlier ones, whose state volume stays; a publication without FRCs cannot be used
    subroutine test_take_cavgs()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)      :: cline
        type(sp_project)   :: set
        type(image)        :: img
        type(string)       :: cwd_saved, root, stk, frcs, problem
        real               :: smpd, mskdiam
        integer            :: nfail0, ncls
        integer, parameter :: NCLS_PUB = 4, NCLS_OLD = 2, BOX = 8
        allocate(stage)
        write(*,'(A)') 'test_take_cavgs'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_take_cavgs', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        ! the stage's earlier class averages and FRCs (a solve2D's), and a state volume
        call write_stack(string('old_cavgs.mrcs'), NCLS_OLD)
        call simple_touch('old_frcs.bin')
        call stage%spproj%add_cavgs2os_out(string('old_cavgs.mrcs'), VOL_SMPD, 'cavg')
        call stage%spproj%add_frcs2os_out(string('old_frcs.bin'), 'frc2D')
        call img%new([BOX, BOX, BOX], VOL_SMPD)
        call img%write(string('recvol_state01.mrc'))
        call img%kill
        call stage%spproj%add_vol2os_out(string('recvol_state01.mrc'), VOL_SMPD, 1, 'vol')
        ! a publication with more classes
        call simple_mkdir('pub')
        call write_stack(string('pub/00003_cavgs.mrcs'), NCLS_PUB)
        call set%add_cavgs2os_out(string('pub/00003_cavgs.mrcs'), VOL_SMPD, 'cavg', mskdiam=150.)
        problem = stage%publication_problem(set)
        call assert_char('its FRCs are missing', problem%to_char(), 'a publication without FRCs cannot be used')
        call simple_touch('pub/00003_frcs.bin')
        call set%add_frcs2os_out(string('pub/00003_frcs.bin'), 'frc2D')
        problem = stage%publication_problem(set)
        call assert_char('', problem%to_char(), 'one with them can')
        ! its classes, as merge_publication takes them, then its class averages and FRCs
        stage%spproj%os_cls2D = set%os_cls2D
        call stage%take_cavgs(set, string('00003'))
        call stage%spproj%get_cavgs_stk(stk, ncls, smpd)
        call assert_int(NCLS_PUB, ncls, 'the class averages are the publication''s')
        call assert_true(stk%has_substr('quality_selection/00003/00003_cavgs.mrcs'), 'copied into its quality folder')
        call assert_true(file_exists(stk), 'where the copy is')
        call stage%spproj%get_frcs(frcs, 'frc2D')
        call assert_true(frcs%has_substr('quality_selection/00003/00003_frcs.bin'), 'and so are its FRCs')
        call assert_true(file_exists(frcs), 'copied too')
        call stage%spproj%get_mskdiam('cavg', mskdiam)
        call assert_real(150., mskdiam, 1.e-4, 'with its mask diameter')
        call assert_int(NCLS_PUB, stage%spproj%os_cls2D%get_noris(), 'the classes stay the publication''s')
        call assert_true(stage%spproj%isthere_in_osout('vol', 1), 'the state volume stays')
        ! pool 2D removes its older publications
        call del_file('pub/00003_cavgs.mrcs')
        call del_file('pub/00003_frcs.bin')
        call assert_true(file_exists(stk),  'the class averages outlive the publication')
        call assert_true(file_exists(frcs), 'and the FRCs')
        call set%kill
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)

    contains

        subroutine write_stack( fname, n )
            class(string), intent(in) :: fname
            integer,       intent(in) :: n
            integer :: i
            call img%new([BOX, BOX, 1], VOL_SMPD)
            do i = 1,n
                call img%write(fname, i)
            enddo
            call img%kill
        end subroutine write_stack

    end subroutine test_take_cavgs

    !> the first publication into a stage without rows is the first set: the particles it selects
    !! are due for solve2D, and the pool model's selection is kept for a failed one; later
    !! publications change the first set's 2D parameters but not its selection, also when they
    !! lack its stack, and select every other row as before
    subroutine test_first_set()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)                 :: cline
        type(sp_project)              :: set
        type(string)                  :: cwd_saved, root
        integer                       :: nfail0, i
        allocate(stage)
        write(*,'(A)') 'test_first_set'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_first_set', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        ! publication 1: stacks A (3) and B (2), the first particle never updated (deselected);
        ! the pool model would keep particles 2 and 4
        call make_set(set, ['A', 'B'], [3, 2], 2, nrejected=1, icls=1)
        call stage%merge_first_set(set, 1, [.false., .true., .false., .true., .false.])
        call set%kill
        call assert_true(stage%l_solve2D_due,           'the first set is due for solve2D')
        call assert_int(4, count(stage%first_set),      'every particle the publication selects is in it')
        call assert_false(stage%in_first_set(1),        'a particle never updated is not')
        call assert_int(4, stage%nptcls_selected,       'the rows take the publication''s own selection')
        call assert_int(2, count(stage%first_fallback), 'the pool model''s selection is kept for a failed solve2D')
        ! as solve2D and its selection leave them: rows 2 and 4 selected, rows 3 and 5 rejected
        do i = 3,5,2
            call stage%spproj%os_ptcl2D%set_state(i, 0)
            call stage%spproj%os_ptcl3D%set_state(i, 0)
        enddo
        stage%l_solve2D_due = .false.
        ! publication 2 selects every particle, in class 2, with the new stack C (2)
        call make_set(set, ['A', 'B', 'C'], [3, 2, 2], 5, icls=2)
        call stage%merge_publication(set, 2)
        call set%kill
        call assert_int(0, stage%spproj%os_ptcl3D%get_state(3), 'a first-set row its selection rejected stays rejected')
        call assert_int(0, stage%spproj%os_ptcl2D%get_state(5), 'in both segments')
        call assert_int(2, stage%spproj%os_ptcl2D%get_class(3), 'and takes the publication''s 2D parameters')
        call assert_int(1, stage%spproj%os_ptcl3D%get_state(2), 'a first-set row it selected stays selected')
        call assert_int(1, stage%spproj%os_ptcl3D%get_state(1), 'a row outside the first set takes the publication''s selection')
        call assert_int(1, stage%spproj%os_ptcl3D%get_state(6), 'as do the new stack''s rows')
        call assert_false(stage%in_first_set(6),                'which are not in the first set')
        call assert_int(5, stage%nptcls_selected,               'the selected rows are counted')
        ! publication 3 lacks stack B (rows 4 and 5)
        call make_set(set, ['A', 'C'], [3, 2], 5, icls=2)
        call stage%merge_publication(set, 3)
        call set%kill
        call assert_int(1, stage%spproj%os_ptcl3D%get_state(4), 'a first-set row keeps its selection when its stack is missing')
        call assert_int(5, stage%nptcls_selected,               'so the count is unchanged')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_first_set

    !> a failed solve2D: the first set's rows take the pool model's selection of the first
    !! publication and are no longer the first set
    subroutine test_first_set_fallback()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)                 :: cline
        type(sp_project)              :: set
        type(string)                  :: cwd_saved, root
        integer                       :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_first_set_fallback'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_first_set_fallback', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call make_set(set, ['A', 'B'], [3, 2], 2, nrejected=1, icls=1)
        call stage%merge_first_set(set, 1, [.false., .true., .false., .true., .false.])
        call set%kill
        call stage%fallback_first_set()
        call assert_int(2, stage%nptcls_selected,               'the pool model''s selection applies')
        call assert_int(0, stage%spproj%os_ptcl3D%get_state(3), 'which deselects what it rejected')
        call assert_int(1, stage%spproj%os_ptcl2D%get_state(4), 'and keeps what it selected')
        call assert_false(allocated(stage%first_set),           'no row is the first set any more')
        call assert_false(stage%l_solve2D_due,                  'no solve2D is due')
        call assert_int(PHASE_IMPORTING, stage%phase,           'solve3D waits for its particles as before')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_first_set_fallback

    !> the mask diameter comes from every publication taken; a change is taken for the next run, and
    !! a publication without one leaves it
    subroutine test_mskdiam_from_each_publication()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)    :: cline
        type(sp_project) :: set
        type(string)     :: cwd_saved, root
        integer          :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_mskdiam_from_each_publication'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_mskdiam', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call add_cavgs_entry(set, 150.)
        call stage%take_mskdiam(set)
        call assert_real(150., stage%mskdiam, 1.e-4, 'the first publication''s mask diameter')
        call set%kill
        call add_cavgs_entry(set, 170.)
        call stage%take_mskdiam(set)
        call assert_real(170., stage%mskdiam, 1.e-4, 'a later publication''s new one')
        call set%kill
        call add_cavgs_entry(set, 0.)
        call stage%take_mskdiam(set)
        call assert_real(170., stage%mskdiam, 1.e-4, 'one without keeps the latest')
        call set%kill
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_mskdiam_from_each_publication

    !> a publication whose stacks disagree with the rows is found before anything is merged
    subroutine test_rows_problem()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)    :: cline
        type(sp_project) :: set
        type(string)     :: cwd_saved, root, problem
        integer          :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_rows_problem'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_rows_problem', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call make_set(set, ['A'], [3], 2, icls=1)
        call stage%merge_publication(set, 1)
        call set%kill
        call make_set(set, ['A', 'B'], [3, 2], 2, icls=1)
        problem = stage%rows_problem(set)
        call assert_int(0, problem%strlen(), 'a publication that agrees with the rows')
        call set%kill
        call make_set(set, ['A'], [4], 2, icls=1)
        problem = stage%rows_problem(set)
        call assert_true(problem%has_substr('changed size'), 'a stack that changed size is found')
        call assert_int(3, stage%spproj%os_ptcl3D%get_noris(), 'and the rows are untouched')
        call set%kill
        call make_set(set, ['A'], [3], 2, icls=1)
        call set%os_ptcl2D%set(1, 'indstk', 7)
        problem = stage%rows_problem(set)
        call assert_true(problem%has_substr('outside its stack'), 'an image index outside its stack is found')
        call set%kill
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_rows_problem

    !> the addon's cohort: selected rows not active in the frozen solution, appended or selected again
    subroutine test_cohort()
        class(stream_stage_solve3D), allocatable :: stage
        integer, parameter :: STATES(6) = [1, 1, 0, 1, 1, 0]
        integer :: i
        allocate(stage)
        allocate(stage%spproj) ! made by init_params in the stage; this test sets rows directly
        write(*,'(A)') 'test_cohort'
        call assert_int(0, stage%count_frozen(), 'before a run nothing is frozen')
        call stage%spproj%os_ptcl3D%new(6, is_ptcl=.true.)
        do i = 1,6
            call stage%spproj%os_ptcl3D%set_state(i, STATES(i))
        enddo
        call assert_int(4, stage%count_cohort(), 'before a run every selected particle is in the cohort')
        ! rows 1 to 4 were in the latest result, row 4 deselected there
        stage%frozen_active = [.true., .true., .true., .false.]
        call assert_int(3, stage%count_frozen(), 'the frozen particles are those active in the result')
        call assert_int(2, stage%count_cohort(), 'a particle selected again and one appended')
        call stage%spproj%kill
        deallocate(stage%spproj)
        deallocate(stage%frozen_active)
    end subroutine test_cohort

    !> the newest quality folders and those runs started from are kept, and so is the latest addon
    !! iteration; earlier iterations and folders set aside for unfinished jobs go
    subroutine test_retention()
        class(stream_stage_solve3D), allocatable :: stage
        type(string) :: cwd_saved, root
        integer      :: nfail0, i, funit
        allocate(stage)
        write(*,'(A)') 'test_retention'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_retention', cwd_saved, root)
        call simple_mkdir('quality_selection')
        do i = 1,6
            call simple_mkdir('quality_selection/'//int2str_pad(i, 5))
        enddo
        open(newunit=funit, file='quality_selection/runs.txt', action='write')
        write(funit,'(A)') '00002'
        close(funit)
        call stage%prune_quality_dirs()
        call assert_false(dir_exists('quality_selection/00001'), 'an old quality folder goes')
        call assert_true(dir_exists('quality_selection/00002'),  'one a run started from stays')
        call assert_false(dir_exists('quality_selection/00003'), 'like every other old one')
        call assert_true(dir_exists('quality_selection/00004'),  'the newest three stay')
        call assert_true(dir_exists('quality_selection/00006'),  'up to the newest')
        call simple_mkdir('solve3D_addon')
        do i = 1,3
            call simple_mkdir('solve3D_addon/it_'//int2str_pad(i, 1))
        enddo
        call simple_mkdir('solve3D_addon/it_2_unfinished1')
        call simple_mkdir('solve3D_unfinished1')
        stage%naddon_runs = 3
        call stage%prune_run_dirs()
        call assert_false(dir_exists('solve3D_addon/it_1'), 'an addon iteration before the frozen base goes')
        call assert_false(dir_exists('solve3D_addon/it_2'), 'every one')
        call assert_true(dir_exists('solve3D_addon/it_3'),  'the frozen base stays')
        call assert_false(dir_exists('solve3D_addon/it_2_unfinished1'), 'a folder set aside for an unfinished job goes')
        call assert_false(dir_exists('solve3D_unfinished1'), 'for solve3D too')
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_retention

    !> which job starts when, and the mask diameter that fits the class averages
    subroutine test_rules()
        class(stream_stage_solve3D), allocatable :: stage
        allocate(stage)
        write(*,'(A)') 'test_rules'
        ! next_job(phase, selected, cohort, frozen, nstates, cohort of the last refused addon run)
        call assert_int(JOB_NONE,    stage%next_job(PHASE_IMPORTING, 0,   0, 0, 3, -1), 'no particles: no job')
        call assert_int(JOB_NONE,    stage%next_job(PHASE_IMPORTING, 14,  0, 0, 3, -1), 'fewer than 5 per state: no solve3D')
        call assert_int(JOB_SOLVE3D, stage%next_job(PHASE_IMPORTING, 15,  0, 0, 3, -1), '5 per state start solve3D')
        call assert_int(JOB_NONE,    stage%next_job(PHASE_SOLVE3D,   100, 0, 0, 3, -1), 'nothing starts while solve3D runs')
        call assert_int(JOB_NONE,    stage%next_job(PHASE_IDLE, 100, 0,   100,  3, -1), 'no cohort: no addon run')
        call assert_int(JOB_NONE,    stage%next_job(PHASE_IDLE, 114, 14,  100,  3, -1), 'a cohort under 5 per state waits')
        call assert_int(JOB_ADDON,   stage%next_job(PHASE_IDLE, 115, 15,  100,  3, -1), '5 per state start an addon run')
        call assert_int(JOB_NONE,    stage%next_job(PHASE_IDLE, 1099, 99, 1000, 3, -1), 'under 10% of the frozen particles waits')
        call assert_int(JOB_ADDON,   stage%next_job(PHASE_IDLE, 1100, 100, 1000, 3, -1), '10% of them start it')
        call assert_int(JOB_NONE,    stage%next_job(PHASE_IDLE, 120, 20,  100,  3, 20), 'after a failed run, the same cohort waits')
        call assert_int(JOB_ADDON,   stage%next_job(PHASE_IDLE, 121, 21,  100,  3, 20), 'a larger one starts the next')
        call assert_int(JOB_NONE,    stage%next_job(PHASE_ADDON, 200, 100, 100, 3, -1), 'nothing starts while it runs')
        call assert_real(100., stage%fit_mskdiam(100.,  64, 2.0), 1.e-4, 'a mask diameter that fits is kept')
        call assert_real(116., stage%fit_mskdiam(400.,  64, 2.0), 1.e-4, 'one too large for the box gets the box default')
        call assert_real(116., stage%fit_mskdiam(0.,    64, 2.0), 1.e-4, 'none gets the box default')
        call assert_real(116., stage%fit_mskdiam(-7.8,  64, 2.0), 1.e-4, 'a negative one too')
        ! after a rolled-back addon run of 120 with 1000 frozen and 3 states, the next waits for
        ! the cadence step (max(15, 100)) beyond it
        call assert_int(219, stage%retry_cohort(120, 1000, 3), 'the cohort to exceed after a rollback')
        call assert_int(JOB_NONE,  stage%next_job(PHASE_IDLE, 1219, 219, 1000, 3, stage%retry_cohort(120, 1000, 3)),&
            &'a few more particles do not start a retry')
        call assert_int(JOB_ADDON, stage%next_job(PHASE_IDLE, 1220, 220, 1000, 3, stage%retry_cohort(120, 1000, 3)),&
            &'a cadence step more does')
        call check_final_publication()
    end subroutine test_rules

    !> a publication flagged pool_final=yes in its out segment is the pool's final one
    subroutine check_final_publication()
        type(sp_project) :: set
        class(stream_stage_solve3D), allocatable :: stage
        allocate(stage)
        call add_cavgs_entry(set, 150.)
        call assert_false(stage%is_final_publication(set), 'a publication without the flag is not final')
        call set%os_out%set(1, 'pool_final', 'yes')
        call assert_true(stage%is_final_publication(set), 'one with pool_final=yes is')
        call set%kill
    end subroutine check_final_publication

    !> one status message per call, with the stage name, the phase and the particle count
    subroutine test_send_status()
        class(stream_stage_solve3D), allocatable                   :: stage
        type(cmdline)                                   :: cline
        type(stream_pipe)                               :: reader
        type(gui_metadata_stream_solve3D_multistate) :: status
        character(len=:), allocatable                   :: buffer
        type(string)                                    :: cwd_saved, root, stage_name
        integer(c_int)                                  :: fds(2)
        integer :: nfail0, meta_type, solve3D_stage, refine_it, nstates_got, nimported, nlast, tlast
        logical :: l_assigned
        real    :: res
        allocate(stage)
        write(*,'(A)') 'test_send_status'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_status', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call open_loopback(fds)
        call stage%init_gui(-1, int(fds(2)))
        call reader%new(int(fds(1)), -1, max_metadata_size(), 'test reader')
        call stage%send_status()
        call assert_true(reader%receive(buffer), 'a status message is sent')
        if( allocated(buffer) )then
            meta_type = transfer(buffer, meta_type)
            call assert_int(GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_TYPE, meta_type, 'it is a multistate 3D status')
            if( meta_type == GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_TYPE )then
                status     = transfer(buffer, status)
                l_assigned = status%get(stage_name, solve3D_stage, refine_it, nstates_got, nimported, nlast, tlast, res)
                call assert_char('waiting on pool 2D', stage_name%to_char(), 'before pool 2D''s folder exists')
                call assert_int(0,       solve3D_stage, 'no 3D yet')
                call assert_int(NSTATES, nstates_got,      'the number of states')
            endif
        endif
        call assert_false(reader%receive(buffer), 'one message per call')
        call reader%kill
        call stage%kill
        call close_loopback(fds)
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_send_status

    !> after a run, each state with a volume is sent with its reprojection tiles, and its
    !! resolution is kept; a state without a volume is not sent, and a state without an FSC
    !! curve is sent without one (not with the previous state's)
    subroutine test_send_volumes()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)                 :: cline
        type(stream_pipe)             :: reader
        type(image)                   :: vol
        type(gui_metadata_vol3D)      :: meta_vol
        character(len=:), allocatable :: buffer
        real,             allocatable :: res(:), fsc(:), invres_got(:), fsc_got(:)
        type(string)                  :: cwd_saved, root, volfile, fscfile
        integer(c_int)                :: fds(2)
        integer                       :: nfail0, meta_type, nvols, ntiles, iptcl, k
        logical                       :: l_fsc_state(2)
        allocate(stage)
        write(*,'(A)') 'test_send_volumes'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_volumes', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call open_loopback(fds)
        call stage%init_gui(-1, int(fds(2)))
        call reader%new(int(fds(1)), -1, max_metadata_size(), 'test reader')
        ! state 1: a volume, its FSC and its reprojections; four particles in it
        volfile = 'recvol_state01.mrc'
        call vol%new([VOL_BOX, VOL_BOX, VOL_BOX], VOL_SMPD)
        call vol%write(volfile)
        call vol%kill
        res = get_resarr(VOL_BOX, VOL_SMPD)
        allocate(fsc(size(res)))
        do k = 1,size(fsc)
            fsc(k) = max(0., 1. - real(k - 1) / real(size(fsc) - 2))
        enddo
        fscfile = 'fsc_state01.bin'
        call arr2file(fsc, fscfile)
        call simple_touch(refine3D_reprojs_fname(1))
        call stage%spproj%add_vol2os_out(volfile, VOL_SMPD, 1, 'vol')
        call stage%spproj%add_fsc2os_out(fscfile, 1, VOL_BOX)
        ! state 2: a volume only; two particles in it
        volfile = 'recvol_state02.mrc'
        call vol%new([VOL_BOX, VOL_BOX, VOL_BOX], VOL_SMPD)
        call vol%write(volfile)
        call vol%kill
        call stage%spproj%add_vol2os_out(volfile, VOL_SMPD, 2, 'vol')
        call stage%spproj%os_ptcl3D%new(6, is_ptcl=.true.)
        do iptcl = 1,6
            call stage%spproj%os_ptcl3D%set_state(iptcl, merge(1, 2, iptcl <= 4))
        enddo
        call stage%send_volumes()
        nvols       = 0
        ntiles      = 0
        l_fsc_state = .false.
        do while( reader%receive(buffer) )
            meta_type = transfer(buffer, meta_type)
            if( meta_type == GUI_METADATA_VOL3D_TYPE )then
                nvols    = nvols + 1
                meta_vol = transfer(buffer, meta_vol)
                if( meta_vol%get_state() >= 1 .and. meta_vol%get_state() <= 2 )then
                    l_fsc_state(meta_vol%get_state()) = meta_vol%get_fsc(invres_got, fsc_got)
                endif
            endif
            if( meta_type == GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_REPROJ_TYPE ) ntiles = ntiles + 1
        enddo
        call assert_int(2, nvols,  'one volume message per state with a volume')
        call assert_int(3, ntiles, 'and the three reprojection tiles of state 1')
        call assert_true(l_fsc_state(1),                    'state 1 is sent with its FSC curve')
        call assert_false(l_fsc_state(2),                   'state 2, without one, is sent without a curve')
        call assert_true(stage%state_res(1) > 0.,          'the resolution of state 1 is kept for the status')
        call assert_real(0., stage%state_res(2), 1.e-6,     'a state without an FSC curve has none')
        call assert_real(0., stage%state_res(3), 1.e-6,     'nor does a state without a volume')
        call reader%kill
        call stage%kill
        call close_loopback(fds)
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_send_volumes

    !> a 3D snapshot request: before a result, answered as not written; from a result, the
    !! particles of the selected states merged into state 1, the state volumes out of the project
    !! and the selected ones copied beside it, answered with the file and the particle count; each
    !! request is written once
    subroutine test_snapshot3D()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)                      :: cline
        type(stream_pipe)                  :: writer, reader
        type(gui_metadata_stream_update)   :: update
        type(gui_metadata_stream_snapshot) :: report
        type(sp_project)                   :: proj
        type(image)                        :: vol
        character(len=:), allocatable      :: buffer
        integer,          allocatable      :: states(:)
        type(string)                       :: cwd_saved, root, fname, volfile
        integer(c_int)                     :: upd(2), gui(2)
        integer                            :: nfail0, id, nptcls, tstamp, istate, iptcl
        integer, parameter                 :: PTCL_STATES(6) = [1, 1, 2, 2, 3, 0]
        allocate(stage)
        write(*,'(A)') 'test_snapshot3D'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_snapshot3D', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call open_loopback(upd)
        call open_loopback(gui)
        call stage%init_gui(int(upd(1)), int(gui(2)))
        call writer%new(-1, int(upd(2)), max_metadata_size(), 'test writer')
        call reader%new(int(gui(1)), -1, max_metadata_size(), 'test reader')
        ! before a result: answered as not written
        call send_request(1, [1])
        call stage%apply_gui_updates()
        call assert_true(receive_report(), 'a request before a result is answered')
        call assert_int(1, id,                 'with its id')
        call assert_int(0, nptcls,             'and no particles')
        call assert_char('', fname%to_char(),  'and no file')
        ! a result: six particles in states 1, 1, 2, 2, 3 and none, and a volume per state
        call proj%os_ptcl2D%new(size(PTCL_STATES), is_ptcl=.true.)
        call proj%os_ptcl3D%new(size(PTCL_STATES), is_ptcl=.true.)
        do iptcl = 1,size(PTCL_STATES)
            call proj%os_ptcl2D%set_state(iptcl, PTCL_STATES(iptcl))
            call proj%os_ptcl3D%set_state(iptcl, PTCL_STATES(iptcl))
        enddo
        do istate = 1,NSTATES
            volfile = 'recvol_state'//int2str_pad(istate, 2)//'.mrc'
            call vol%new([VOL_BOX, VOL_BOX, VOL_BOX], VOL_SMPD)
            call vol%write(volfile)
            call vol%kill
            call proj%add_vol2os_out(volfile, VOL_SMPD, istate, 'vol')
        enddo
        call proj%write(string('result.simple'))
        call proj%kill
        stage%result_projfile = 'result.simple'
        ! states 1 and 3
        call send_request(2, [1, 3])
        call stage%apply_gui_updates()
        call assert_true(receive_report(), 'the request is answered')
        call assert_int(2, id,     'with its id')
        call assert_int(3, nptcls, 'and the particles of states 1 and 3')
        states = report%get_states()
        call assert_int(2, size(states), 'and its states')
        call assert_true(file_exists(fname), 'the snapshot project is written')
        if( file_exists(fname) )then
            call proj%read(fname)
            call assert_int(3, proj%os_ptcl3D%count_state_gt_zero(),           'with the selected particles')
            call assert_int(1, maxval(proj%os_ptcl3D%get_all_asint('state')), 'merged into state 1')
            call assert_int(3, proj%os_ptcl2D%count_state_gt_zero(),           'in both particle segments')
            call assert_false(proj%isthere_in_osout('vol', 1),                 'without the state volumes')
            call proj%kill
        endif
        call assert_true(file_exists(string('snapshots/snapshot_2/vol_state01.mrc')),  'the volume of state 1 is copied')
        call assert_true(file_exists(string('snapshots/snapshot_2/vol_state03.mrc')),  'and that of state 3')
        call assert_false(file_exists(string('snapshots/snapshot_2/vol_state02.mrc')), 'but not that of state 2')
        ! the same request again
        call send_request(2, [1, 3])
        call stage%apply_gui_updates()
        call assert_false(reader%receive(buffer), 'a request is written once')
        call update%kill
        call writer%kill
        call reader%kill
        call stage%kill
        call close_loopback(upd)
        call close_loopback(gui)
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)

    contains

        subroutine send_request( snapshot_id, selection )
            integer, intent(in) :: snapshot_id, selection(:)
            call update%kill
            call update%new(GUI_METADATA_STREAM_UPDATE_TYPE)
            call update%set_snapshot3D_update(snapshot_id, selection, string('snapshot_'//int2str(snapshot_id)//METADATA_EXT))
            call writer%send_meta(update)
        end subroutine send_request

        ! the latest 3D snapshot report the stage sent, into report, id, fname and nptcls
        logical function receive_report()
            integer :: meta_type
            receive_report = .false.
            do while( reader%receive(buffer) )
                meta_type = transfer(buffer, meta_type)
                if( meta_type /= GUI_METADATA_STREAM_SOLVE3D_SNAPSHOT_TYPE ) cycle
                report         = transfer(buffer, report)
                receive_report = report%get(id, fname, nptcls, tstamp)
            enddo
        end function receive_report

    end subroutine test_snapshot3D

    !> the public loop waits for pool 2D's folder, then attaches; no job without particles
    subroutine test_iterate_waits()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)                 :: cline
        type(string)                  :: cwd_saved, root
        integer                       :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_iterate_waits'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_iterate', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%iterate()
        call assert_false(stage%l_attached, 'pass 1: no pool 2D folder, not attached')
        call simple_mkdir(UPSTREAM)
        call simple_mkdir(UPSTREAM//'/'//DIR_STREAM_COMPLETED)
        call stage%iterate()
        call assert_true(stage%l_attached,            'pass 2: attached')
        call assert_int(PHASE_IMPORTING, stage%phase, 'pass 2: no particles, no job')
        call stage%kill
        call stage%kill ! idempotence
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_iterate_waits

    subroutine test_finished()
        class(stream_stage_solve3D), allocatable :: stage
        type(cmdline)                 :: cline
        type(string)                  :: cwd_saved, root
        integer                       :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_finished'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_finished', cwd_saved, root)
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
        call cline%set('prg',        'stream_solve3D_stage_tester')
        call cline%set('mkdir',      'no')
        call cline%set('projfile',   TEST_PROJFILE)
        call cline%set('outdir',     '')
        call cline%set('nthr',       1)
        call cline%set('nparts',     1)
        call cline%set('nstates',    NSTATES)
        call cline%set('dir_target', UPSTREAM)
        call cline%set('qsys_name',  'local')
    end subroutine set_test_cline

    ! a stage without a queue environment or a pipe, with no waits and a settle time that takes
    ! files written in the same second
    subroutine make_test_stage( stage, cline )
        class(stream_stage_solve3D), intent(inout) :: stage
        type(cmdline),                 intent(inout) :: cline
        call stage%init_params(cline)
        call stage%init_gui(-1, -1)
        stage%settle_s = -1
        stage%wait_s   = 0
        stage%l_exists = .true.
    end subroutine make_test_stage

    ! export number @p id of pool 2D, with only the class averages' entry and its mask diameter
    subroutine write_export( id, mskdiam )
        integer, intent(in) :: id
        real,    intent(in) :: mskdiam
        type(sp_project) :: proj
        type(string)     :: cwd, fname
        call simple_getcwd(cwd)
        fname = cwd//'/'//UPSTREAM//'/'//DIR_STREAM_COMPLETED//int2str_pad(id, 5)//METADATA_EXT
        call proj%os_out%new(1, is_ptcl=.false.)
        call proj%os_out%set(1, 'imgkind', 'cavg')
        call proj%os_out%set(1, 'stk',     'cavgs.mrc')
        call proj%os_out%set(1, 'nptcls',  10)
        call proj%os_out%set(1, 'mskdiam', mskdiam)
        call proj%update_projinfo(fname)
        call proj%write(fname)
        call proj%kill
    end subroutine write_export

    ! a publication in memory: one micrograph and stack per entry of @p stks (stack name
    ! 'stack_<name>.mrcs') with @p nptcls_stk particles each, whose 'x' is the stack's position in
    ! the alphabet; the first @p nrejected particles deselected (as the class-average selection
    ! leaves them, in both segments); every particle in class @p icls when given; @p ncls classes
    subroutine make_set( set, stks, nptcls_stk, ncls, nrejected, icls )
        type(sp_project),  intent(inout) :: set
        character(len=1),  intent(in)    :: stks(:)
        integer,           intent(in)    :: nptcls_stk(:), ncls
        integer, optional, intent(in)    :: nrejected, icls
        integer :: istk, iptcl, fromp, nptcls
        nptcls = sum(nptcls_stk)
        call set%os_mic%new(size(stks), is_ptcl=.false.)
        call set%os_stk%new(size(stks), is_ptcl=.false.)
        call set%os_ptcl2D%new(nptcls, is_ptcl=.true.)
        call set%os_ptcl3D%new(nptcls, is_ptcl=.true.)
        call set%os_cls2D%new(ncls, is_ptcl=.false.)
        fromp = 1
        do istk = 1,size(stks)
            call set%os_mic%set_state(istk, 1)
            call set%os_stk%set(istk, 'stk',    'stack_'//stks(istk)//'.mrcs')
            call set%os_stk%set(istk, 'nptcls', nptcls_stk(istk))
            call set%os_stk%set(istk, 'fromp',  fromp)
            call set%os_stk%set(istk, 'top',    fromp + nptcls_stk(istk) - 1)
            do iptcl = fromp,fromp + nptcls_stk(istk) - 1
                call set%os_ptcl2D%set_state(iptcl, 1)
                call set%os_ptcl3D%set_state(iptcl, 1)
                call set%os_ptcl3D%set(iptcl, 'x', real(iachar(stks(istk)) - iachar('A') + 1))
                if( present(icls) ) call set%os_ptcl2D%set_class(iptcl, icls)
            enddo
            fromp = fromp + nptcls_stk(istk)
        enddo
        if( present(nrejected) )then
            do iptcl = 1,nrejected
                call set%os_ptcl2D%set_state(iptcl, 0)
                call set%os_ptcl3D%set_state(iptcl, 0)
            enddo
        endif
    end subroutine make_set

    ! a publication's class averages entry with mask diameter @p mskdiam (none when 0)
    subroutine add_cavgs_entry( set, mskdiam )
        type(sp_project), intent(inout) :: set
        real,             intent(in)    :: mskdiam
        call set%os_out%new(1, is_ptcl=.false.)
        call set%os_out%set(1, 'imgkind', 'cavg')
        call set%os_out%set(1, 'stk',     'cavgs.mrc')
        call set%os_out%set(1, 'nptcls',  10)
        if( mskdiam > 0. ) call set%os_out%set(1, 'mskdiam', mskdiam)
    end subroutine add_cavgs_entry

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

end module simple_stream_stage_solve3D_tester
