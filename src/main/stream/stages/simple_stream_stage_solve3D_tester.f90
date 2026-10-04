!@descr: unit tests for the steps of the stream multistate 3D stage (simple_stream_stage_solve3D)
! Each test assembles a stage in a fresh fixture directory from init_params and init_gui: no queue
! environment, no waits, and a settle time of -1. The command line names no registered program
! and carries qsys_name=local for the stage project's computing environment. The upstream is a pool
! 2D directory whose completed folder holds exports. The class-average selection and the 3D jobs
! (which need a queue) are left to the high-level stream tests; the import is tested on sets made
! in memory, and the volume messages on a fixture project with a volume and an FSC.
module simple_stream_stage_solve3D_tester
use, intrinsic :: iso_c_binding, only: c_int
use unix,                        only: c_pipe, c_close, c_fcntl, F_GETFL, F_SETFL, O_NONBLOCK
use simple_test_utils
use simple_defs_fname,                                only: TERM_STREAM, METADATA_EXT, JPG_EXT
use simple_defs_stream,                               only: DIR_STREAM_COMPLETED
use simple_string,                                    only: string
use simple_string_utils,                              only: int2str_pad
use simple_fileio,                                    only: arr2file, del_file, file_exists, simple_getcwd, simple_touch
use simple_syslib,                                    only: simple_mkdir
use simple_math_ft,                                   only: get_resarr
use simple_cmdline,                                   only: cmdline
use simple_image,                                     only: image
use simple_sp_project,                                only: sp_project
use simple_rec_list,                                  only: rec_iterator, chunk_rec
use simple_gui_metadata_utils,                        only: max_metadata_size
use simple_gui_metadata_types,                        only: GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_TYPE, GUI_METADATA_VOL3D_TYPE,&
                                                           &GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_REPROJ_TYPE
use simple_gui_metadata_stream_solve3D_multistate, only: gui_metadata_stream_solve3D_multistate
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
        call test_rules()
        call test_send_status()
        call test_send_volumes()
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

    !> exports are taken in export order whatever order they are listed in, each once; the first
    !! gives pool 2D's mask diameter
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
        call assert_real(150., stage%mskdiam, 1.e-4, 'the mask diameter of the first export')
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
        integer                       :: nfail0
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
        ! a 3D result for the first particle, which later publications must keep
        call stage%spproj%os_ptcl3D%set(1, 'e1', 33.)
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
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_merge_publications

    !> which job starts when, and the mask diameter that fits the class averages
    subroutine test_rules()
        class(stream_stage_solve3D), allocatable :: stage
        allocate(stage)
        write(*,'(A)') 'test_rules'
        call assert_int(JOB_NONE,     stage%next_job(PHASE_IMPORTING, 0,   0),   'no particles: no job')
        call assert_int(JOB_SOLVE3D, stage%next_job(PHASE_IMPORTING, 10,  0),   'the first particles start solve3D')
        call assert_int(JOB_NONE,     stage%next_job(PHASE_SOLVE3D,  100, 0),   'nothing starts while solve3D runs')
        call assert_int(JOB_NONE,     stage%next_job(PHASE_IDLE,      100, 100), 'no growth: no addon run')
        call assert_int(JOB_ADDON,    stage%next_job(PHASE_IDLE,      120, 100), 'growth starts an addon run')
        call assert_int(JOB_NONE,     stage%next_job(PHASE_ADDON,     200, 100), 'nothing starts while it runs')
        call assert_real(100., stage%fit_mskdiam(100.,  64, 2.0), 1.e-4, 'a mask diameter that fits is kept')
        call assert_real(116., stage%fit_mskdiam(400.,  64, 2.0), 1.e-4, 'one too large for the box gets the box default')
        call assert_real(116., stage%fit_mskdiam(0.,    64, 2.0), 1.e-4, 'none gets the box default')
        call assert_real(116., stage%fit_mskdiam(-7.8,  64, 2.0), 1.e-4, 'a negative one too')
    end subroutine test_rules

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
        logical :: l_assigned, l_user_input
        real    :: res
        allocate(stage)
        write(*,'(A)') 'test_send_status'
        nfail0 = tests_failed
        call enter_fixture('a3_stage_status', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call open_loopback(fds)
        call stage%init_gui(int(fds(2)))
        call reader%new(int(fds(1)), -1, max_metadata_size(), 'test reader')
        call stage%send_status()
        call assert_true(reader%receive(buffer), 'a status message is sent')
        if( allocated(buffer) )then
            meta_type = transfer(buffer, meta_type)
            call assert_int(GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_TYPE, meta_type, 'it is a multistate 3D status')
            if( meta_type == GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_TYPE )then
                status     = transfer(buffer, status)
                l_assigned = status%get(stage_name, solve3D_stage, refine_it, nstates_got, nimported, nlast, tlast,&
                    &l_user_input, res)
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
        call stage%init_gui(int(fds(2)))
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
        call simple_touch('orthogonal_reprojs_state01'//JPG_EXT)
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
        call stage%init_gui(-1)
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
