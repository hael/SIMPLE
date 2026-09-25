!@descr: stream preprocessing workflow test: simulated movies through the stream's preprocessing stage and its worker jobs
! The only stream stage test that submits jobs (local queue), hence a workflow entry of its own
! (stream_preproc). The other stream stage tests of this module moved by the stream review (plan,
! section 9.7): optics assignment, picking references and pick and extract to lib_stream
! (simple_stream_tester), the sieve's collect-and-reject to the particle sieve tests (unit_project),
! the master heartbeat to the forked_process platform entry (simple_gui_assembler_tester);
! abinitio2D_stream was retired. doc/refactoring_notes/stream_area_tests_handover.md says what this
! test should compare with the simulation truth.
module simple_commanders_test_stream
use simple_commanders_api
use simple_commanders_sim, only: commander_simulate_movie
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_stream_preproc
  contains
    procedure :: execute => exec_test_stream_preproc
end type commander_test_stream_preproc

contains

subroutine exec_test_stream_preproc( self, cline )
    use simple_stream_p01_preprocess_new, only: stream_p01_preprocess
    use simple_ui,                        only: make_ui
    use, intrinsic :: ieee_arithmetic,    only: ieee_is_finite
    class(commander_test_stream_preproc), intent(inout) :: self
    class(cmdline),                intent(inout) :: cline
    character(len=*), parameter :: MOVIE_FILE      = 'simulate_movie.mrc'
    character(len=*), parameter :: OPTIMAL_FILE    = 'optimal_movie_average.mrc'
    character(len=*), parameter :: PARAMS_FILE     = 'simulate_movie_params.txt'
    character(len=*), parameter :: PARTICLE_STACK  = 'synthetic_particles.mrcs'
    character(len=*), parameter :: TEST_PROJFILE   = 'test_preproc.simple'
    character(len=*), parameter :: TEST_OUTDIR     = 'stream_preproc'
    real,             parameter :: SMPD             = 1.3
    real,             parameter :: CS               = 2.7
    real,             parameter :: KV               = 300.0
    real,             parameter :: FRACA            = 0.1
    integer,          parameter :: PARTICLE_BOX     = 64
    integer,          parameter :: NPARTICLES       = 12
    integer,          parameter :: MOVIE_DIM        = 512
    integer,          parameter :: NFRAMES          = 8
    type(stream_p01_preprocess)  :: xpreproc
    type(commander_simulate_movie) :: xsimov
    type(cmdline)                :: cline_preproc, cline_sim_mov
    type(image)                  :: particle, feature
    type(sp_project)             :: spproj
    type(string)                 :: cwd_root, fixture_root, movies_dir, particle_stack_path
    type(string)                 :: project_path, movie_name, value
    character(len=XLONGSTRLEN)   :: fixture_root_path
    integer                      :: i, status, ldim(3), nimages
    real                         :: dfx, dfy, ctfres

    ! simple_test_exec initializes only test-program metadata.  The movie
    ! simulator uses the production UI while constructing its parameters.
    call make_ui
    call simple_getcwd(cwd_root)
    fixture_root = 'test_preproc_'//int2str(get_process_id())
    if( dir_exists(fixture_root) ) THROW_HARD('TEST_PREPROC FAILED: fixture directory already exists')
    call simple_mkdir(fixture_root)
    call simple_chdir(fixture_root, status)
    if( status /= 0 ) THROW_HARD('TEST_PREPROC FAILED: could not enter fixture directory')
    call simple_getcwd(fixture_root)
    fixture_root_path = fixture_root%to_char()

    ! Build a small asymmetric particle stack.  The production movie simulator
    ! distributes these particles, applies known CTF parameters, adds motion and
    ! noise, and therefore provides realistic input to both preprocessing steps.
    call particle%new([PARTICLE_BOX, PARTICLE_BOX, 1], SMPD, wthreads=.false.)
    call particle%square(7)
    call feature%new([PARTICLE_BOX, PARTICLE_BOX, 1], SMPD, wthreads=.false.)
    call feature%square(4)
    call feature%shift([12., -8., 0.])
    call particle%add(feature)
    do i = 1, NPARTICLES
        call particle%write(string(PARTICLE_STACK), i, del_if_exists=(i == 1))
    enddo
    call particle%kill
    call feature%kill
    particle_stack_path = simple_abspath(string(PARTICLE_STACK))

    write(logfhandle,'(a,i0,a)') '>>> TEST_PREPROC: generating ', STREAM_NMOVS_SET, ' synthetic movies'
    do i = 1, STREAM_NMOVS_SET
        call cline_sim_mov%set('prg',       'simulate_movie')
        call cline_sim_mov%set('stk',       particle_stack_path)
        call cline_sim_mov%set('xdim',      MOVIE_DIM)
        call cline_sim_mov%set('ydim',      MOVIE_DIM)
        call cline_sim_mov%set('nframes',   NFRAMES)
        call cline_sim_mov%set('smpd',      SMPD)
        call cline_sim_mov%set('snr',       0.5)
        call cline_sim_mov%set('kv',        KV)
        call cline_sim_mov%set('cs',        CS)
        call cline_sim_mov%set('fraca',     FRACA)
        call cline_sim_mov%set('defocus',   1.5 + 0.25 * real(i - 1))
        call cline_sim_mov%set('trs',       1.0)
        call cline_sim_mov%set('nthr',      1)
        if( i == 1 )then
            call cline_sim_mov%set('mkdir',   'yes')
            call cline_sim_mov%set('dir_exec','simulate_movies')
        else
            call cline_sim_mov%set('mkdir',   'no')
        endif
        call xsimov%execute(cline_sim_mov)
        if( i == 1 ) call simple_getcwd(movies_dir)
        if( .not. file_exists(MOVIE_FILE) ) THROW_HARD('TEST_PREPROC FAILED: synthetic movie was not generated')
        movie_name = 'synthetic_movie_'//int2str_pad(i, 3)//MRC_EXT
        call simple_rename(MOVIE_FILE, movie_name)
        if( file_exists(OPTIMAL_FILE) ) call del_file(OPTIMAL_FILE)
        if( file_exists(PARAMS_FILE)  ) call del_file(PARAMS_FILE)
        call cline_sim_mov%kill
    enddo
    call simple_chdir(fixture_root_path, status)
    if( status /= 0 ) THROW_HARD('TEST_PREPROC FAILED: could not leave movie directory')

    project_path = filepath(fixture_root, TEST_PROJFILE)
    call cline_preproc%set('prg',                'preproc')
    call cline_preproc%set('projfile',           project_path)
    call cline_preproc%set('dir_movies',         movies_dir)
    call cline_preproc%set('outdir',             TEST_OUTDIR)
    call cline_preproc%set('qsys_name',          'local')
    call cline_preproc%set('nparts',             1)
    call cline_preproc%set('nthr',               1)
    call cline_preproc%set('nmics',              STREAM_NMOVS_SET)
    call cline_preproc%set('smpd',               SMPD)
    call cline_preproc%set('cs',                 CS)
    call cline_preproc%set('kv',                 KV)
    call cline_preproc%set('fraca',              FRACA)
    call cline_preproc%set('total_dose',         40.0)
    call cline_preproc%set('flipgain',           'none')
    call cline_preproc%set('algorithm',          'iso')
    call cline_preproc%set('mcpatch',            'no')
    call cline_preproc%set('ctfpatch',           'no')
    call cline_preproc%set('pspecsz',            256)
    call cline_preproc%set('ctfresthreshold',    100.0)
    call cline_preproc%set('icefracthreshold',   100.0)
    call cline_preproc%set('astigthreshold',     100.0)
    call cline_preproc%printline(unit=6)
    call xpreproc%execute(cline_preproc)

    project_path = cline_preproc%get_carg('projfile')
    call simple_chdir(fixture_root_path, status)
    if( status /= 0 ) THROW_HARD('TEST_PREPROC FAILED: could not leave preprocessing directory')
    if( .not. file_exists(project_path) ) THROW_HARD('TEST_PREPROC FAILED: output project was not created')
    call spproj%read_segment('mic', project_path)
    if( spproj%os_mic%get_noris() /= STREAM_NMOVS_SET )&
        &THROW_HARD('TEST_PREPROC FAILED: project does not contain five processed micrographs')

    do i = 1, STREAM_NMOVS_SET
        call spproj%os_mic%getter(i, 'imgkind', value)
        if( trim(value%to_char()) /= 'mic' )&
            &THROW_HARD('TEST_PREPROC FAILED: movie record was not converted to a micrograph')
        call assert_output_file(i, 'intg')
        call spproj%os_mic%getter(i, 'intg', value)
        call find_ldim_nptcls(value, ldim, nimages)
        if( any(ldim(1:2) /= MOVIE_DIM) .or. nimages /= 1 )&
            &THROW_HARD('TEST_PREPROC FAILED: integrated micrograph has unexpected dimensions')
        call assert_output_file(i, 'mc_starfile')
        call assert_output_file(i, 'ctfjpg')
        if( .not. spproj%os_mic%isthere(i, 'dfx') .or. .not. spproj%os_mic%isthere(i, 'dfy') .or. &
            &.not. spproj%os_mic%isthere(i, 'ctfres') )&
            &THROW_HARD('TEST_PREPROC FAILED: CTF parameters are missing')
        dfx    = spproj%os_mic%get(i, 'dfx')
        dfy    = spproj%os_mic%get(i, 'dfy')
        ctfres = spproj%os_mic%get(i, 'ctfres')
        if( .not. ieee_is_finite(dfx) .or. .not. ieee_is_finite(dfy) .or. &
            &.not. ieee_is_finite(ctfres) .or. dfx <= 0. .or. dfy <= 0. .or. ctfres <= 0. )&
            &THROW_HARD('TEST_PREPROC FAILED: CTF parameters are not physically valid')
        if( abs(spproj%os_mic%get(i, 'smpd') - SMPD) > 0.01 )&
            &THROW_HARD('TEST_PREPROC FAILED: micrograph sampling distance changed unexpectedly')
    enddo
    call spproj%kill
    call cline_preproc%kill
    call simple_chdir(cwd_root, status)
    if( status /= 0 ) THROW_HARD('TEST_PREPROC FAILED: could not restore the original directory')
    write(logfhandle,'(a,a)') '>>> TEST_PREPROC: validated motion and CTF outputs in ', fixture_root%to_char()
    call simple_end('**** SIMPLE_TEST_PREPROC NORMAL STOP ****')

  contains

    subroutine assert_output_file( imic, key )
        integer,          intent(in) :: imic
        character(len=*), intent(in) :: key
        type(string) :: fname
        if( .not. spproj%os_mic%isthere(imic, key) )&
            &THROW_HARD('TEST_PREPROC FAILED: missing '//key//' project field')
        call spproj%os_mic%getter(imic, key, fname)
        if( .not. file_exists(fname) ) THROW_HARD('TEST_PREPROC FAILED: missing '//key//' output file')
    end subroutine assert_output_file

end subroutine exec_test_stream_preproc

end module simple_commanders_test_stream
