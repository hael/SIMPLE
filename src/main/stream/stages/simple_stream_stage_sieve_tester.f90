!@descr: unit tests for the steps of the stream particle-sieving stage (simple_stream_stage_sieve)
! Each test assembles a stage in a fresh fixture directory from init_params and init_gui: the
! stage's own queue environment (which starts the persistent chunk workers) is not made, no waits,
! and a settle time of -1. The command line names no registered program and carries
! qsys_name=local for the stage project's computing environment, which the sieve's own queue
! environment reads. The upstream is a reference-picking directory with completed sets and the
! moldiam.txt that make_pickrefs writes. Chunk 2D runs and rejection are covered by the
! ptcl_sieve tests and the high-level stream tests.
module simple_stream_stage_sieve_tester
use, intrinsic :: iso_c_binding, only: c_int
use unix,                        only: c_pipe, c_close, c_fcntl, F_GETFL, F_SETFL, O_NONBLOCK
use simple_test_utils
use simple_defs_fname,                           only: TERM_STREAM, STREAM_MOLDIAM, METADATA_EXT
use simple_defs_stream,                          only: DIR_STREAM_COMPLETED, STREAM_IDLE_MARKER
use simple_string,                               only: string
use simple_string_utils,                         only: int2str_pad
use simple_fileio,                               only: del_file, file_exists, simple_getcwd, simple_touch
use simple_syslib,                               only: dir_exists, simple_mkdir
use simple_cmdline,                              only: cmdline
use simple_oris,                                 only: oris
use simple_sp_project,                           only: sp_project
use simple_gui_metadata_utils,                   only: max_metadata_size
use simple_gui_metadata_types,                   only: GUI_METADATA_STREAM_PARTICLE_SIEVING_TYPE
use simple_gui_metadata_stream_particle_sieving, only: gui_metadata_stream_particle_sieving
use simple_stream_pipe,                          only: stream_pipe
use simple_stream_stage_sieve,                   only: stream_stage_sieve
use simple_ptcl_sieve,                           only: CHUNKED_MICS
implicit none
private
public :: run_all_stream_stage_sieve_tests

character(len=*), parameter :: TEST_PROJFILE = 'test_sieve.simple'
character(len=*), parameter :: UPSTREAM      = 'refpick'   ! the reference-picking stage directory (dir_target)
real,             parameter :: MSKDIAM_REFS  = 136.        ! the mask diameter make_pickrefs decided

contains

    subroutine run_all_stream_stage_sieve_tests()
        write(*,'(A)') '**** running all stream particle sieving stage tests ****'
        call test_init_params()
        call test_restart_removes_term_stream()
        call test_restore_and_attach()
        call test_restart_resumes_sieve()
        call test_import_projects()
        call test_read_mask_diameter()
        call test_start_sieve()
        call test_final_ingestion_follows_marker()
        call test_send_status()
        call test_iterate_waits()
        call test_finished()
    end subroutine run_all_stream_stage_sieve_tests

    subroutine test_init_params()
        class(stream_stage_sieve), allocatable :: stage
        type(cmdline)            :: cline
        type(string)             :: cwd_saved, root
        integer                  :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_init_params'
        nfail0 = tests_failed
        call enter_fixture('sv_stage_init_params', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call assert_true(allocated(stage%params), 'initialization allocates the owned parameters')
        call assert_int(0, stage%spproj%os_mic%get_noris(),         'the project starts without micrographs')
        call assert_true(dir_exists(string(DIR_STREAM_COMPLETED)),  'the completed folder is made')
        call assert_false(stage%l_sieve_active,                     'no sieve before the first import')
        call assert_false(allocated(stage%sieve), 'initialization leaves the sieve unallocated')
        allocate(stage%sieve)
        call stage%kill
        call assert_false(allocated(stage%params), 'cleanup releases the owned parameters')
        call assert_false(allocated(stage%sieve), 'cleanup releases an inactive sieve')
        call stage%kill
        call assert_false(allocated(stage%sieve), 'inactive-sieve cleanup is idempotent')
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_init_params

    !> an existing output folder is a restart, and a leftover termination file is removed, so the
    !! restarted stage runs
    subroutine test_restart_removes_term_stream()
        class(stream_stage_sieve), allocatable :: stage
        type(cmdline)            :: cline
        type(string)             :: cwd_saved, root
        integer                  :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_restart_removes_term_stream'
        nfail0 = tests_failed
        call enter_fixture('sv_stage_restart', cwd_saved, root)
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

    !> restart: a set the sieve has chunked from (its chunked_mics.txt) is imported again with the
    !! chunked micrographs marked, so the rest of a partly chunked set is still sieved, and goes into
    !! the watcher history once the upstream folder is attached
    subroutine test_restore_and_attach()
        class(stream_stage_sieve), allocatable :: stage
        type(cmdline)            :: cline
        type(string)             :: cwd_saved, root, set1, set2
        integer                  :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_restore_and_attach'
        nfail0 = tests_failed
        call enter_fixture('sv_stage_restore', cwd_saved, root)
        call simple_mkdir(UPSTREAM)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%attach_upstream()
        call assert_false(stage%l_attached,      'no completed folder: not attached')
        call assert_true(stage%l_waiting_logged, 'the wait is logged')
        call make_upstream()
        set1 = write_completed_set(1, [10, 20])
        set2 = write_completed_set(2, [30])
        ! set 1 is partly chunked: its first micrograph only
        call write_chunked_mics(set1, [1])
        call stage%restore_imports()
        call assert_true(allocated(stage%restored_imports), 'the sets chunked from are read back')
        call assert_int(2,  stage%project_list%size(),  'every micrograph of a partly chunked set is imported again')
        call assert_int(20, stage%project_list%get_nptcls_tot(l_not_included=.true.),&
            &'only its unchunked micrograph (the second, 20 particles) is left for the sieve')
        call assert_int(2,  stage%n_mics_imported,      'the restored micrographs are counted')
        call assert_int(30, stage%n_ptcls_imported,     'and their particles')
        call stage%attach_upstream()
        call assert_true(stage%l_attached, 'attached once the folder exists')
        call assert_true(stage%project_buff%is_past(set1),  'an already chunked project is in the history')
        call assert_false(stage%project_buff%is_past(set2), 'a project not yet chunked is not')
        call assert_false(allocated(stage%restored_imports), 'the restart''s list is used once')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_restore_and_attach

    !> a restart with chunks to take up makes the sieve with nothing new to import
    subroutine test_restart_resumes_sieve()
        class(stream_stage_sieve), allocatable :: stage
        type(cmdline)            :: cline
        type(string)             :: cwd_saved, root
        integer                  :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_restart_resumes_sieve'
        nfail0 = tests_failed
        call enter_fixture('sv_stage_resume', cwd_saved, root)
        call make_upstream()
        call write_moldiam()
        call simple_mkdir('previous')
        call set_test_cline(cline)
        call cline%set('outdir', 'previous')
        call make_test_stage(stage, cline)
        call assert_true(stage%resumable(), 'a restart with the picking references'' mask diameter can resume')
        call stage%iterate()
        call assert_true(stage%l_attached,     'attached')
        call assert_int(0, stage%project_list%size(), 'nothing new to import')
        call assert_true(stage%l_sieve_active, 'the sieve is made all the same')
        call stage%finalize()
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_restart_resumes_sieve

    !> one record per micrograph of each newly completed set, each set once
    subroutine test_import_projects()
        class(stream_stage_sieve), allocatable :: stage
        type(cmdline)            :: cline
        type(string)             :: cwd_saved, root, set1, set2
        integer                  :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_import_projects'
        nfail0 = tests_failed
        call enter_fixture('sv_stage_import', cwd_saved, root)
        call make_upstream()
        set1 = write_completed_set(1, [10, 20])
        set2 = write_completed_set(2, [30])
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%attach_upstream()
        call stage%import_projects()
        call assert_int(3,  stage%project_list%size(), 'one record per micrograph')
        call assert_int(3,  stage%n_mics_imported,     'the micrographs are counted')
        call assert_int(60, stage%n_ptcls_imported,    'their particles are counted')
        call stage%import_projects()
        call assert_int(3,  stage%project_list%size(), 'a set is imported once')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_import_projects

    !> the mask diameter is the one make_pickrefs wrote beside the picking references
    subroutine test_read_mask_diameter()
        class(stream_stage_sieve), allocatable :: stage
        type(cmdline)            :: cline
        type(string)             :: cwd_saved, root
        integer                  :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_read_mask_diameter'
        nfail0 = tests_failed
        call enter_fixture('sv_stage_mskdiam', cwd_saved, root)
        call make_upstream()
        call write_moldiam()
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%read_mask_diameter()
        call assert_real(MSKDIAM_REFS, stage%mskdiam, 1.e-4, 'the mask diameter of the picking references')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_read_mask_diameter

    !> the sieve is made on the first import, with the picking references' mask diameter, and runs
    !! its warm-up cycles; too few particles for a chunk, so none is generated
    subroutine test_start_sieve()
        class(stream_stage_sieve), allocatable :: stage
        type(cmdline)            :: cline
        type(string)             :: cwd_saved, root, set1
        integer                  :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_start_sieve'
        nfail0 = tests_failed
        call enter_fixture('sv_stage_start', cwd_saved, root)
        call make_upstream()
        call write_moldiam()
        set1 = write_completed_set(1, [10, 20])
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%attach_upstream()
        call stage%import_projects()
        call assert_false(allocated(stage%sieve), 'import alone does not allocate the sieve')
        call stage%start_sieve()
        call assert_true(stage%l_sieve_active, 'the sieve is made')
        call assert_true(allocated(stage%sieve), 'starting the sieve allocates owned state')
        call assert_real(MSKDIAM_REFS, stage%mskdiam, 1.e-4, 'with the picking references'' mask diameter')
        call assert_int(0, stage%sieve%get_n_chunks_coarse(), 'too few particles for a chunk')
        call assert_true(dir_exists(string('chunks_coarse')), 'the sieve''s chunk folders are made')
        call stage%kill
        call assert_false(allocated(stage%sieve), 'cleanup releases an active sieve')
        call assert_false(stage%l_sieve_active, 'cleanup resets sieve activity')
        call stage%kill
        call assert_false(allocated(stage%sieve), 'active-sieve cleanup is idempotent')
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_start_sieve

    !> final ingestion is set once reference picking's idle marker has been seen and a later watch
    !! found nothing, and withdrawn when the marker goes
    subroutine test_final_ingestion_follows_marker()
        class(stream_stage_sieve), allocatable :: stage
        type(cmdline)            :: cline
        type(string)             :: cwd_saved, root, set1
        integer                  :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_final_ingestion_follows_marker'
        nfail0 = tests_failed
        call enter_fixture('sv_stage_final', cwd_saved, root)
        call make_upstream()
        call write_moldiam()
        set1 = write_completed_set(1, [10, 20])
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%attach_upstream()
        call stage%import_projects()
        call stage%start_sieve()
        call stage%update_final_ingestion()
        call assert_false(stage%l_final, 'no marker: no final ingestion')
        call simple_touch(UPSTREAM//'/'//STREAM_IDLE_MARKER)
        call stage%update_final_ingestion()
        call assert_false(stage%l_final, 'the marker first seen: wait for a later watch')
        call assert_true(stage%upstream_done_since > 0, 'the time it was seen is kept')
        stage%last_watch = stage%upstream_done_since + 1
        call stage%update_final_ingestion()
        call assert_true(stage%l_final, 'a later watch found nothing: final ingestion')
        call del_file(UPSTREAM//'/'//STREAM_IDLE_MARKER)
        call stage%update_final_ingestion()
        call assert_false(stage%l_final, 'the marker gone: final ingestion withdrawn')
        call assert_int(0, stage%upstream_done_since, 'and the wait starts again')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_final_ingestion_follows_marker

    !> one status message per call, with the particle count and the classes the sieve selected
    subroutine test_send_status()
        class(stream_stage_sieve), allocatable                   :: stage
        type(cmdline)                              :: cline
        type(stream_pipe)                          :: reader
        type(gui_metadata_stream_particle_sieving) :: status
        character(len=:), allocatable              :: buffer
        type(string)                               :: cwd_saved, root, stage_name
        integer(c_int)                             :: fds(2)
        integer                                    :: nfail0, meta_type, nimported, naccepted, nrejected, tlast
        logical                                    :: l_assigned
        allocate(stage)
        write(*,'(A)') 'test_send_status'
        nfail0 = tests_failed
        call enter_fixture('sv_stage_status', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call open_loopback(fds)
        call stage%init_gui(-1, int(fds(2)))
        call reader%new(int(fds(1)), -1, max_metadata_size(), 'test reader')
        stage%n_ptcls_imported = 50
        stage%latest_inds      = [3, 5, 7]
        stage%latest_selection = [1, 0, 1]
        call stage%send_status(string('importing and sieving particles'))
        call assert_false(allocated(stage%sieve), 'sending initial status does not allocate a sieve')
        call assert_true(reader%receive(buffer), 'a status message is sent')
        if( allocated(buffer) )then
            meta_type = transfer(buffer, meta_type)
            call assert_int(GUI_METADATA_STREAM_PARTICLE_SIEVING_TYPE, meta_type, 'it is a particle-sieving status')
            if( meta_type == GUI_METADATA_STREAM_PARTICLE_SIEVING_TYPE )then
                status     = transfer(buffer, status)
                l_assigned = status%get(stage_name, nimported, naccepted, nrejected, tlast)
                call assert_int(50, nimported, 'the particles imported')
                call assert_int(0, naccepted, 'initial status has no accepted particles')
                call assert_int(0, nrejected, 'initial status has no rejected particles')
                call assert_char('importing and sieving particles', stage_name%to_char(), 'the stage name')
            endif
        endif
        call assert_false(reader%receive(buffer), 'one message per call')
        call reader%kill
        call stage%kill
        call close_loopback(fds)
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_send_status

    !> the public loop waits for the upstream folder, then attaches; no sieve until a set arrives
    subroutine test_iterate_waits()
        class(stream_stage_sieve), allocatable :: stage
        type(cmdline)            :: cline
        type(string)             :: cwd_saved, root
        integer                  :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_iterate_waits'
        nfail0 = tests_failed
        call enter_fixture('sv_stage_iterate', cwd_saved, root)
        call simple_mkdir(UPSTREAM)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%iterate()
        call assert_false(stage%l_attached, 'pass 1: no upstream output, not attached')
        call make_upstream()
        call stage%iterate()
        call assert_true(stage%l_attached,      'pass 2: attached')
        call assert_false(stage%l_sieve_active, 'pass 2: no set yet, no sieve')
        call stage%kill
        call stage%kill ! idempotence
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_iterate_waits

    subroutine test_finished()
        class(stream_stage_sieve), allocatable :: stage
        type(cmdline)            :: cline
        type(string)             :: cwd_saved, root
        integer                  :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_finished'
        nfail0 = tests_failed
        call enter_fixture('sv_stage_finished', cwd_saved, root)
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
    ! system for the computing environment of the stage's project (the sieve's queue reads it)
    subroutine set_test_cline( cline )
        type(cmdline), intent(inout) :: cline
        call cline%set('prg',         'stream_sieve_stage_tester')
        call cline%set('mkdir',       'no')
        call cline%set('projfile',    TEST_PROJFILE)
        call cline%set('outdir',      '')
        call cline%set('nthr',        1)
        call cline%set('nparts',      1)
        call cline%set('nchunks',     1)
        call cline%set('walltime',    60)
        call cline%set('single_pass', 'yes')
        call cline%set('use_model',   'no')
        call cline%set('dir_target',  UPSTREAM)
        call cline%set('qsys_name',   'local')
    end subroutine set_test_cline

    ! a stage without its own queue environment or a pipe, with no waits and a settle time that
    ! takes files written in the same second
    subroutine make_test_stage( stage, cline )
        class(stream_stage_sieve), intent(inout) :: stage
        type(cmdline),            intent(inout) :: cline
        call stage%init_params(cline)
        call stage%init_gui(-1, -1)
        stage%settle_s = -1
        stage%wait_s   = 0
        stage%l_exists = .true.
    end subroutine make_test_stage

    ! the folders reference picking makes
    subroutine make_upstream()
        call simple_mkdir(UPSTREAM)
        call simple_mkdir(UPSTREAM//'/'//DIR_STREAM_COMPLETED)
    end subroutine make_upstream

    ! the moldiam.txt make_pickrefs writes in the reference-picking directory
    subroutine write_moldiam()
        type(oris) :: moldiam
        call moldiam%new(1, is_ptcl=.false.)
        call moldiam%set(1, 'mskdiam',         MSKDIAM_REFS)
        call moldiam%set(1, 'box_for_extract', 160)
        call moldiam%write(string(UPSTREAM//'/'//STREAM_MOLDIAM))
        call moldiam%kill
    end subroutine write_moldiam

    ! the sieve's record of the micrographs @p micinds of set @p set put in a chunk
    subroutine write_chunked_mics( set, micinds )
        type(string), intent(in) :: set
        integer,      intent(in) :: micinds(:)
        integer :: funit, i
        open(newunit=funit, file=CHUNKED_MICS, status='replace', action='write')
        do i = 1,size(micinds)
            write(funit,'(A,1X,I0)') set%to_char(), micinds(i)
        enddo
        close(funit)
    end subroutine write_chunked_mics

    ! completed reference-picking set number @p id: one micrograph and one stack per entry of
    ! @p nptcls_mic; returns its absolute path
    function write_completed_set( id, nptcls_mic ) result( fname )
        integer, intent(in) :: id, nptcls_mic(:)
        type(string)     :: fname
        type(sp_project) :: proj
        type(string)     :: cwd
        integer          :: imic, fromp
        call simple_getcwd(cwd)
        fname = cwd//'/'//UPSTREAM//'/'//DIR_STREAM_COMPLETED//int2str_pad(id, 5)//METADATA_EXT
        call proj%os_mic%new(size(nptcls_mic), is_ptcl=.false.)
        call proj%os_stk%new(size(nptcls_mic), is_ptcl=.false.)
        fromp = 1
        do imic = 1,size(nptcls_mic)
            call proj%os_mic%set_state(imic, 1)
            call proj%os_mic%set(imic, 'nptcls',    nptcls_mic(imic))
            call proj%os_mic%set(imic, 'importind', (id - 1) * 10 + imic)
            call proj%os_stk%set(imic, 'fromp',     fromp)
            call proj%os_stk%set(imic, 'top',       fromp + nptcls_mic(imic) - 1)
            fromp = fromp + nptcls_mic(imic)
        enddo
        call proj%update_projinfo(fname)
        call proj%write(fname)
        call proj%kill
    end function write_completed_set

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

end module simple_stream_stage_sieve_tester
