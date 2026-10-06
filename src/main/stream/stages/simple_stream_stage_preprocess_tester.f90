!@descr: unit tests for the steps of the stream preprocessing stage (simple_stream_stage_preprocess)
! Each test assembles a stage in a fresh fixture directory with only the steps it needs:
! parameters from a command line that names no registered program (so params%new neither
! requires a project nor makes an output directory), the job directories, GUI metadata on
! no pipe or on a loopback pipe, and zero waits. Nothing is submitted to a queue: job
! submission and collection, and gain resolution from movies, are left to the high-level
! stream_preproc test.
module simple_stream_stage_preprocess_tester
use, intrinsic :: iso_c_binding, only: c_int
use unix,                        only: c_pipe, c_close, c_fcntl, F_GETFL, F_SETFL, O_NONBLOCK
use simple_test_utils
use simple_defs_fname,                 only: TERM_STREAM, GAIN_THUMBNAIL
use simple_defs_stream,                only: DIR_STREAM_COMPLETED, STREAM_NMOVS_SET
use simple_string,                     only: string
use simple_string_utils,               only: int2str, int2str_pad
use simple_fileio,                     only: file_exists, del_file, simple_getcwd, simple_touch
use simple_syslib,                     only: simple_mkdir
use simple_cmdline,                    only: cmdline
use simple_image,                      only: image
use simple_oris,                       only: oris
use simple_sp_project,                 only: sp_project
use simple_stream_watcher,             only: stream_watcher
use simple_gui_metadata_utils,         only: max_metadata_size
use simple_gui_metadata_types,         only: GUI_METADATA_STREAM_UPDATE_TYPE, GUI_METADATA_STREAM_PREPROCESS_TYPE
use simple_gui_metadata_stream_update, only: gui_metadata_stream_update
use simple_stream_pipe,                only: stream_pipe
use simple_stream_stage_preprocess,    only: stream_stage_preprocess
implicit none
private
public :: run_all_stream_stage_preprocess_tests

character(len=*), parameter :: TEST_PROJFILE = 'test_preprocess_stage.simple'
real,             parameter :: SMPD = 1.3, CS = 2.7, KV = 300.0, FRACA = 0.1
real,             parameter :: CTFRES_THRES = 10.0, ICEFRAC_THRES = 1.0, ASTIG_THRES = 10.0

contains

    subroutine run_all_stream_stage_preprocess_tests()
        write(*,'(A)') '**** running all stream preprocessing stage tests ****'
        call test_init_params_split_mode()
        call test_import_completed()
        call test_process_imports()
        call test_import_previous_projects()
        call test_create_movies_set_project()
        call test_build_worker_cline()
        call test_resolve_gain_static_flip()
        call test_apply_gui_updates()
        call test_send_status()
        call test_finished()
    end subroutine run_all_stream_stage_preprocess_tests

    !> parameter derivation resets split_mode to 'even'; the stage must restore 'stream', or the
    !! queue gets one partition for ncunits computing units (out-of-bounds jobs_done)
    subroutine test_init_params_split_mode()
        class(stream_stage_preprocess), allocatable :: stage
        type(cmdline)                 :: cline
        type(string)                  :: cwd_saved, root
        integer                       :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_init_params_split_mode'
        nfail0 = tests_failed
        call enter_fixture('pp_stage_split_mode', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call assert_true(allocated(stage%params), 'initialization allocates the owned parameters')
        call assert_char('stream', trim(stage%params%split_mode), 'split_mode is stream after init_params')
        call assert_false(stage%l_restart, 'no output directory given: not a restart')
        call stage%kill
        call assert_false(allocated(stage%params), 'cleanup releases the owned parameters')
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_init_params_split_mode

    !> two finished jobs: the accepted micrographs reach the global project, the others count as
    !! failed, the thresholds are written back to the job project, and only a job with accepted
    !! micrographs moves to the completed folder
    subroutine test_import_completed()
        class(stream_stage_preprocess), allocatable :: stage
        type(cmdline)                 :: cline
        type(sp_project)              :: job
        type(string)                  :: cwd_saved, root, job_dir, done(2), partial(1)
        integer                       :: nfail0, n_imported
        allocate(stage)
        write(*,'(A)') 'test_import_completed'
        nfail0 = tests_failed
        call enter_fixture('pp_stage_import_completed', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        job_dir = stage%sets%get_job_dir()
        done(1) = job_dir//'/00001.simple'
        done(2) = job_dir//'/00002.simple'
        call write_mic_project(done(1), [1,0,1,1,1], [4.,4.,4.,15.,4.], 'job1_')
        call write_mic_project(done(2), [0,0,0,0,0], [4.,4.,4.,4.,4.],  'job2_')
        call stage%import_completed(done, n_imported)
        call assert_int(4, n_imported,                           'the four accepted micrographs are imported')
        call assert_int(4, stage%spproj_glob%os_mic%get_noris(), 'the global project holds them')
        call assert_int(6, stage%n_failed_jobs,                  'the other six count as failed')
        call assert_int(1, stage%spproj_glob%os_mic%get_state(3), 'thresholds are not applied to the global project here')
        call assert_true(file_exists(string(DIR_STREAM_COMPLETED//'00001.simple')), 'job 1 moves to the completed folder')
        call assert_true(file_exists(done(2)),                                      'job 2, with nothing accepted, stays')
        if( file_exists(string(DIR_STREAM_COMPLETED//'00001.simple')) )then
            call job%read_segment('mic', string(DIR_STREAM_COMPLETED//'00001.simple'))
            call assert_int(0, job%os_mic%get_state(4), 'the over-threshold micrograph is rejected in the job project')
            call assert_int(1, job%os_mic%get_state(1), 'an accepted micrograph stays accepted in the job project')
            call job%kill
        endif
        ! a partial set of three micrographs
        partial(1) = job_dir//'/00003.simple'
        call write_mic_project(partial(1), [1,1,0], [4.,4.,4.], 'job3_')
        call stage%import_completed(partial, n_imported)
        call assert_int(2, n_imported,                           'a partial set: its two accepted micrographs')
        call assert_int(6, stage%spproj_glob%os_mic%get_noris(), 'join the global project')
        call assert_int(7, stage%n_failed_jobs,                  'and its third counts as failed')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_import_completed

    !> after an import the thresholds reject micrographs of the global project, and below 1000
    !! micrographs the STAR file is written on every import
    subroutine test_process_imports()
        class(stream_stage_preprocess), allocatable :: stage
        type(cmdline)                 :: cline
        type(string)                  :: cwd_saved, root
        real,    parameter            :: CTFRES(5) = [4., 15., 4., 4., 4.]
        integer                       :: nfail0, imic
        allocate(stage)
        write(*,'(A)') 'test_process_imports'
        nfail0 = tests_failed
        call enter_fixture('pp_stage_process_imports', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%spproj_glob%os_mic%new(5, is_ptcl=.false.)
        do imic = 1,5
            call set_mic(stage%spproj_glob%os_mic, imic, 1, CTFRES(imic), string('movie_'//int2str(imic)//'.mrcs'))
        enddo
        call stage%process_imports()
        call assert_int(0, stage%spproj_glob%os_mic%get_state(2), 'ctfres above the threshold rejects')
        call assert_int(1, stage%spproj_glob%os_mic%get_state(1), 'ctfres below the threshold keeps')
        call assert_true(file_exists(string('micrographs.star')), 'the micrographs STAR file is written')
        call assert_true(stage%l_haschanged,                       'the import marks a change')
        call assert_int(0, stage%nmic_star,                        'below 1000 micrographs the snapshot count is not advanced')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_process_imports

    !> restart: the set counter continues after the highest completed set even when that set has
    !! no accepted micrograph; accepted micrographs come back; every movie processed, rejected ones
    !! too, goes into the history, and import indices continue after the highest given
    subroutine test_import_previous_projects()
        class(stream_stage_preprocess), allocatable :: stage
        type(cmdline)                 :: cline
        type(string)                  :: cwd_saved, root
        integer                       :: nfail0, imic
        allocate(stage)
        write(*,'(A)') 'test_import_previous_projects'
        nfail0 = tests_failed
        call enter_fixture('pp_stage_previous_projects', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        stage%movie_buff = stream_watcher(0, string('movies'))
        call write_mic_project(string(DIR_STREAM_COMPLETED//'00007.simple'), [1,1,1], [4.,4.,4.], 'set7_', importind0=1)
        call write_mic_project(string(DIR_STREAM_COMPLETED//'00009.simple'), [0,0],   [4.,4.],    'set9_', importind0=4)
        ! the watcher only records movies that exist
        do imic = 1,3
            call simple_touch(string('set7_'//int2str(imic)//'.mrcs'))
        enddo
        do imic = 1,2
            call simple_touch(string('set9_'//int2str(imic)//'.mrcs'))
        enddo
        call stage%import_previous_projects()
        call assert_int(9, stage%sets%get_counter(),              'the set counter continues after the highest set')
        call assert_int(3, stage%spproj_glob%os_mic%get_noris(),  'the accepted micrographs are re-imported')
        call assert_int(5, stage%import_counter,                  'the import counter continues after the highest index given')
        call assert_true(stage%movie_buff%is_past(string('set7_1.mrcs')), 'accepted movies are in the watcher history')
        call assert_true(stage%movie_buff%is_past(string('set9_1.mrcs')), 'so are rejected ones: none is processed again')
        call assert_true(file_exists(string(DIR_STREAM_COMPLETED//'00009.simple')), 'a set with nothing accepted is kept')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_import_previous_projects

    !> one set of movies becomes a project in the job directory: consecutive import indices, the
    !! XML metadata path without the _fractions suffix, and the worker command line pointed at it
    subroutine test_create_movies_set_project()
        class(stream_stage_preprocess), allocatable :: stage
        type(cmdline)                 :: cline
        type(sp_project)              :: set_proj
        type(string)                  :: cwd_saved, root, cwd, movies(STREAM_NMOVS_SET), meta, projfile
        integer                       :: nfail0, imov
        allocate(stage)
        write(*,'(A)') 'test_create_movies_set_project'
        nfail0 = tests_failed
        call enter_fixture('pp_stage_movies_set', cwd_saved, root)
        call simple_getcwd(cwd)
        call simple_mkdir('xml')
        call simple_mkdir('movies')
        do imov = 1,STREAM_NMOVS_SET
            movies(imov) = cwd//'/movies/mov_'//int2str_pad(imov, 3)//'_fractions.mrcs'
            call write_movie(movies(imov))
        enddo
        call set_test_cline(cline)
        call cline%set('dir_meta', 'xml')
        call make_test_stage(stage, cline)
        call stage%create_movies_set_project(movies)
        call assert_int(1, stage%sets%get_counter(), 'the first set is number 1')
        projfile = stage%sets%get_job_dir()//'/'//int2str_pad(1, 5)//'.simple'
        call assert_true(file_exists(projfile), 'the set project is written in the job directory')
        if( file_exists(projfile) )then
            call set_proj%read_segment('mic', projfile)
            call assert_int(STREAM_NMOVS_SET, set_proj%os_mic%get_noris(), 'the set project holds the movies')
            if( set_proj%os_mic%get_noris() == STREAM_NMOVS_SET )then
                call assert_int(1, set_proj%os_mic%get_int(1, 'importind'), 'the first import index is 1')
                call assert_int(STREAM_NMOVS_SET, set_proj%os_mic%get_int(STREAM_NMOVS_SET, 'importind'),&
                    &'import indices are consecutive')
                meta = set_proj%os_mic%get_str(1, 'meta')
                call assert_true(meta%has_substr('/mov_001.xml'), 'the XML path drops the _fractions suffix')
            endif
            call set_proj%kill
        endif
        projfile = stage%cline_exec%get_carg('projfile')
        call assert_char('00001.simple', projfile%to_char(), 'the worker command line names the set project')
        ! the movies short of a set, as a partial set
        call stage%create_movies_set_project(movies(1:3))
        projfile = stage%sets%get_job_dir()//'/'//int2str_pad(2, 5)//'.simple'
        call assert_true(file_exists(projfile), 'a partial set is written')
        if( file_exists(projfile) )then
            call set_proj%read_segment('mic', projfile)
            call assert_int(3, set_proj%os_mic%get_noris(), 'with its three movies')
            call set_proj%kill
        endif
        call assert_int(3, stage%cline_exec%get_iarg('top'), 'and its job processes three')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_create_movies_set_project

    !> the worker command line runs preprocess on one set in the job directory, with the gain
    !! reference as resolved and no further flipping
    subroutine test_build_worker_cline()
        class(stream_stage_preprocess), allocatable :: stage
        type(cmdline)                 :: cline
        type(string)                  :: val
        allocate(stage)
        write(*,'(A)') 'test_build_worker_cline'
        call cline%set('prg',      'preprocess_stream')
        call cline%set('flipgain', 'x')
        call cline%set('gainref',  'gain.mrc')
        call stage%build_worker_cline(cline)
        val = stage%cline_exec%get_carg('prg')
        call assert_char('preprocess', val%to_char(), 'workers run preprocess')
        val = stage%cline_exec%get_carg('mkdir')
        call assert_char('no',         val%to_char(), 'workers make no output directory')
        val = stage%cline_exec%get_carg('dir')
        call assert_char('../',        val%to_char(), 'workers write to the stage directory')
        call assert_int(1,                stage%cline_exec%get_iarg('fromp'), 'one set: first movie')
        call assert_int(STREAM_NMOVS_SET, stage%cline_exec%get_iarg('top'),   'one set: last movie')
        val = stage%cline_exec%get_carg('flipgain')
        call assert_char('no',         val%to_char(), 'workers do not flip the resolved gain again')
        val = stage%cline_exec%get_carg('gainref')
        call assert_char('gain.mrc',   val%to_char(), 'workers get the gain reference')
        call stage%cline_exec%kill
        call cline%kill
    end subroutine test_build_worker_cline

    !> a static flip is applied once, by the stage: a flipped copy of the gain reference is
    !! written and both the parameters and the command line point to it
    subroutine test_resolve_gain_static_flip()
        class(stream_stage_preprocess), allocatable :: stage
        type(cmdline)                 :: cline
        type(string)                  :: cwd_saved, root, gainref
        integer                       :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_resolve_gain_static_flip'
        nfail0 = tests_failed
        call enter_fixture('pp_stage_gain_flip', cwd_saved, root)
        call write_gainref(string('gainref.mrc'))
        call simple_touch(GAIN_THUMBNAIL) ! the thumbnail is not under test
        call set_test_cline(cline)
        call cline%set('gainref',  'gainref.mrc')
        call cline%set('flipgain', 'x')
        call make_test_stage(stage, cline)
        call stage%resolve_gain(cline)
        gainref = cline%get_carg('gainref')
        call assert_true(gainref%has_substr('_flipX'), 'the command line names the flipped gain reference')
        call assert_true(file_exists(gainref),                  'the flipped gain reference is written')
        call assert_char(gainref%to_char(), stage%gainref%to_char(), 'the stage names the same file')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_resolve_gain_static_flip

    !> two threshold updates queued before one drain are both applied, to the parameters and to
    !! the worker command line; a field left at 0 (unset) changes nothing
    subroutine test_apply_gui_updates()
        class(stream_stage_preprocess), allocatable    :: stage
        type(cmdline)                    :: cline
        type(stream_pipe)                :: writer
        type(gui_metadata_stream_update) :: update
        character(len=:), allocatable    :: buffer
        type(string)                     :: cwd_saved, root
        integer(c_int)                   :: fds(2)
        integer                          :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_apply_gui_updates'
        nfail0 = tests_failed
        call enter_fixture('pp_stage_gui_updates', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call open_loopback(fds)
        call stage%init_gui(int(fds(1)), -1)
        call writer%new(-1, int(fds(2)), max_metadata_size(), 'test writer')
        call update%new(GUI_METADATA_STREAM_UPDATE_TYPE)
        call update%set_ctfres_update(7.5)
        call update%serialise(buffer)
        call writer%send(buffer)
        call update%kill
        call update%new(GUI_METADATA_STREAM_UPDATE_TYPE)
        call update%set_astigmatism_update(3.0)
        call update%serialise(buffer)
        call writer%send(buffer)
        call update%kill
        call stage%apply_gui_updates()
        call assert_real(7.5,           stage%ctfres_thres,  1.e-6, 'the ctfres update is applied')
        call assert_real(3.0,           stage%astig_thres,   1.e-6, 'the queued astigmatism update is applied')
        call assert_real(ICEFRAC_THRES, stage%icefrac_thres, 1.e-6, 'an unset field changes nothing')
        call assert_real(7.5, stage%cline_exec%get_rarg('ctfresthreshold'), 1.e-6, 'workers get the new ctfres threshold')
        call assert_real(3.0, stage%cline_exec%get_rarg('astigthreshold'),  1.e-6, 'workers get the new astigmatism threshold')
        call writer%kill
        call stage%kill
        call close_loopback(fds)
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_apply_gui_updates

    !> one status message per call, of the preprocessing status type
    subroutine test_send_status()
        class(stream_stage_preprocess), allocatable :: stage
        type(cmdline)                 :: cline
        type(stream_pipe)             :: reader
        character(len=:), allocatable :: buffer
        type(string)                  :: cwd_saved, root
        integer(c_int)                :: fds(2)
        integer                       :: nfail0, imic, meta_type
        allocate(stage)
        write(*,'(A)') 'test_send_status'
        nfail0 = tests_failed
        call enter_fixture('pp_stage_status', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call open_loopback(fds)
        call stage%init_gui(-1, int(fds(2)))
        call reader%new(int(fds(1)), -1, max_metadata_size(), 'test reader')
        call stage%spproj_glob%os_mic%new(3, is_ptcl=.false.)
        do imic = 1,3
            call set_mic(stage%spproj_glob%os_mic, imic, merge(0, 1, imic == 2), 4., string('movie_'//int2str(imic)//'.mrcs'))
        enddo
        stage%n_failed_jobs = 2
        call stage%send_status()
        call assert_true(reader%receive(buffer), 'a status message is sent')
        if( allocated(buffer) )then
            meta_type = transfer(buffer, meta_type)
            call assert_int(GUI_METADATA_STREAM_PREPROCESS_TYPE, meta_type, 'it is a preprocessing status')
        endif
        call assert_false(reader%receive(buffer), 'exactly one message per call')
        call reader%kill
        call stage%kill
        call close_loopback(fds)
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_send_status

    !> the stage finishes on the stream's termination file or once enough micrographs are imported
    subroutine test_finished()
        class(stream_stage_preprocess), allocatable :: stage
        type(cmdline)                 :: cline
        type(string)                  :: cwd_saved, root
        integer                       :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_finished'
        nfail0 = tests_failed
        call enter_fixture('pp_stage_finished', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call assert_false(stage%finished(), 'a fresh stage is not finished')
        call simple_touch(TERM_STREAM)
        call assert_true(stage%finished(), 'the termination file finishes the stage')
        call del_file(TERM_STREAM)
        stage%l_nmics_reached = .true.
        call assert_true(stage%finished(), 'reaching the requested micrographs finishes the stage')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_finished

    ! ---- fixtures ------------------------------------------------------------

    ! the stage's command line, under a program name outside every UI table, so params%new
    ! neither requires a project nor makes an output directory
    subroutine set_test_cline( cline )
        type(cmdline), intent(inout) :: cline
        call cline%set('prg',              'stream_preprocess_stage_tester')
        call cline%set('mkdir',            'no')
        call cline%set('projfile',         TEST_PROJFILE)
        call cline%set('nthr',             1)
        call cline%set('nparts',           1)
        call cline%set('numlen',           5)
        call cline%set('stream',           'yes')
        call cline%set('dir_movies',       'movies')
        call cline%set('smpd',             SMPD)
        call cline%set('cs',               CS)
        call cline%set('kv',               KV)
        call cline%set('fraca',            FRACA)
        call cline%set('ctfresthreshold',  CTFRES_THRES)
        call cline%set('icefracthreshold', ICEFRAC_THRES)
        call cline%set('astigthreshold',   ASTIG_THRES)
    end subroutine set_test_cline

    ! a stage with parameters, job directories, GUI metadata on no pipe, and no waits
    subroutine make_test_stage( stage, cline )
        class(stream_stage_preprocess), intent(inout) :: stage
        type(cmdline),                 intent(inout) :: cline
        type(sp_project) :: proj
        call simple_mkdir('movies')
        ! an existing project file, so init_params does not build one from the environment
        call proj%update_projinfo(string(TEST_PROJFILE))
        call proj%write(string(TEST_PROJFILE))
        call proj%kill
        stage%settle_s     = 0
        stage%wait_s       = 0
        stage%sniff_wait_s = 0
        call stage%init_params(cline)
        call stage%init_job_dirs()
        call stage%init_gui(-1, -1)
        stage%cline_exec = cline
        stage%l_exists   = .true.
    end subroutine make_test_stage

    ! a project of micrographs with the given states and ctfres; movies are <prefix><i>.mrcs
    subroutine write_mic_project( fname, states, ctfres, prefix, importind0 )
        class(string),     intent(in) :: fname
        integer,           intent(in) :: states(:)
        real,              intent(in) :: ctfres(:)
        character(len=*),  intent(in) :: prefix
        integer, optional, intent(in) :: importind0 ! the first micrograph's import index
        type(sp_project) :: proj
        integer :: imic
        call proj%os_mic%new(size(states), is_ptcl=.false.)
        do imic = 1,size(states)
            call set_mic(proj%os_mic, imic, states(imic), ctfres(imic), string(prefix//int2str(imic)//'.mrcs'))
            if( present(importind0) ) call proj%os_mic%set(imic, 'importind', importind0 + imic - 1)
        enddo
        call proj%update_projinfo(fname)
        call proj%write(fname)
        call proj%kill
    end subroutine write_mic_project

    subroutine set_mic( os_mic, imic, state, ctfres, movie )
        class(oris),   intent(inout) :: os_mic
        integer,       intent(in)    :: imic, state
        real,          intent(in)    :: ctfres
        class(string), intent(in)    :: movie
        call os_mic%set_state(imic, state)
        call os_mic%set(imic, 'smpd',      SMPD)
        call os_mic%set(imic, 'cs',        CS)
        call os_mic%set(imic, 'kv',        KV)
        call os_mic%set(imic, 'fraca',     FRACA)
        call os_mic%set(imic, 'ctfres',    ctfres)
        call os_mic%set(imic, 'importind', real(imic))
        call os_mic%set(imic, 'movie',     movie)
    end subroutine set_mic

    ! a two-frame 8x8 movie
    subroutine write_movie( fname )
        class(string), intent(in) :: fname
        type(image) :: frame
        integer     :: iframe
        call frame%new([8,8,1], SMPD, wthreads=.false.)
        do iframe = 1,2
            frame = real(iframe)
            call frame%write(fname, iframe, del_if_exists=(iframe == 1))
        enddo
        call frame%kill
    end subroutine write_movie

    ! a 16x16 gain reference with a gradient, so a flip changes it
    subroutine write_gainref( fname )
        class(string), intent(in) :: fname
        type(image) :: gain
        real        :: rmat(16,16,1)
        integer     :: i, j
        do j = 1,16
            do i = 1,16
                rmat(i,j,1) = real(i + 10 * j)
            enddo
        enddo
        call gain%new([16,16,1], SMPD, wthreads=.false.)
        call gain%set_rmat(rmat, .false.)
        call gain%write(fname, del_if_exists=.true.)
        call gain%kill
    end subroutine write_gainref

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

end module simple_stream_stage_preprocess_tester
