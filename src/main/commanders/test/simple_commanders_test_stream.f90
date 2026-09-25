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
    use simple_oris,                      only: oris
    use simple_starfile_wrappers,         only: starfile_table_type, starfile_table__new, &
        &starfile_table__read, starfile_table__firstobject, starfile_table__nextobject, &
        &starfile_table__numberofobjects, starfile_table__getValue_double, &
        &starfile_table__delete, EMDL_MICROGRAPH_SHIFT_X, EMDL_MICROGRAPH_SHIFT_Y
    use simple_stream_p01_preprocess_new, only: stream_p01_preprocess
    use simple_ui,                        only: make_ui
    use, intrinsic :: iso_c_binding,      only: C_long, C_double
    use, intrinsic :: ieee_arithmetic,    only: ieee_is_finite
    class(commander_test_stream_preproc), intent(inout) :: self
    class(cmdline),                intent(inout) :: cline
    character(len=*), parameter :: MOVIE_FILE      = 'simulate_movie.mrc'
    character(len=*), parameter :: OPTIMAL_FILE    = 'optimal_movie_average.mrc'
    character(len=*), parameter :: PARAMS_FILE     = 'simulate_movie_params.txt'
    character(len=*), parameter :: PARTICLE_STACK  = 'synthetic_particles.mrcs'
    character(len=*), parameter :: TEST_PROJFILE   = 'test_preproc.simple'
    character(len=*), parameter :: TEST_OUTDIR     = 'stream_preproc'
    character(len=*), parameter :: TRUTH_DIRNAME   = 'simulation_truth'
    real,             parameter :: SMPD             = 1.3
    real,             parameter :: CS               = 2.7
    real,             parameter :: KV               = 300.0
    real,             parameter :: FRACA            = 0.1
    real,             parameter :: DEFOCUS_TOL      = 0.10
    real,             parameter :: MOTION_TOL       = 0.50
    real,             parameter :: CORR_FLOOR       = 0.50
    real,             parameter :: CORR_MARGIN      = 0.05
    integer,          parameter :: PARTICLE_BOX     = 64
    integer,          parameter :: NPARTICLES       = 12
    integer,          parameter :: MOVIE_DIM        = 512
    integer,          parameter :: NFRAMES          = 8
    type(stream_p01_preprocess)  :: xpreproc
    type(commander_simulate_movie) :: xsimov
    type(cmdline)                :: cline_preproc, cline_sim_mov
    type(image)                  :: particle, feature, integrated, optimal, wrong_optimal
    type(oris)                   :: simulation_truth
    type(sp_project)             :: spproj
    type(string)                 :: cwd_root, fixture_root, movies_dir, truth_dir, particle_stack_path
    type(string)                 :: project_path, movie_name, optimal_name, params_name, value
    type(string)                 :: movie_paths(STREAM_NMOVS_SET), optimal_paths(STREAM_NMOVS_SET)
    type(string)                 :: params_paths(STREAM_NMOVS_SET)
    character(len=XLONGSTRLEN)   :: fixture_root_path
    integer                      :: i, imovie, wrong_movie, status, ldim(3), nimages
    real                         :: dfx, dfy, ctfres, truth_dfx, truth_dfy, defocus_error
    real                         :: motion_error, image_corr, wrong_corr

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
    truth_dir = filepath(fixture_root, TRUTH_DIRNAME)
    call simple_mkdir(truth_dir)

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
        optimal_name = 'optimal_movie_average_'//int2str_pad(i, 3)//MRC_EXT
        params_name = 'simulate_movie_params_'//int2str_pad(i, 3)//TXT_EXT
        call simple_rename(MOVIE_FILE, movie_name)
        call simple_rename(OPTIMAL_FILE, filepath(truth_dir, optimal_name))
        call simple_rename(PARAMS_FILE, filepath(truth_dir, params_name))
        movie_paths(i)   = simple_abspath(movie_name)
        optimal_paths(i) = filepath(truth_dir, optimal_name)
        params_paths(i)  = filepath(truth_dir, params_name)
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
        call spproj%os_mic%getter(i, 'movie', value)
        imovie = find_simulation(value)
        if( imovie == 0 ) THROW_HARD('TEST_PREPROC FAILED: output micrograph cannot be matched to simulation truth')
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

        ! The simulator writes the actual defocus values after its deliberate
        ! random perturbation, so this compares against truth rather than the
        ! requested nominal defocus.
        call simulation_truth%new(1, is_ptcl=.false.)
        call simulation_truth%read(params_paths(imovie))
        truth_dfx = simulation_truth%get_dfx(1)
        truth_dfy = simulation_truth%get_dfy(1)
        ! The two astigmatic axes are physically unordered: an equivalent fit
        ! can exchange dfx/dfy and rotate the astigmatism angle by 90 degrees.
        defocus_error = min(max(abs(dfx - truth_dfx), abs(dfy - truth_dfy)), &
            &max(abs(dfx - truth_dfy), abs(dfy - truth_dfx)))
        write(logfhandle,'(a,i0,a,f7.3,a,f7.3)') '>>> TEST_PREPROC movie ', imovie, &
            &': defocus max error ', defocus_error, ' um; limit ', DEFOCUS_TOL
        if( defocus_error > DEFOCUS_TOL )&
            &THROW_HARD('TEST_PREPROC FAILED: recovered defocus differs from simulation truth by more than 0.10 um')

        call spproj%os_mic%getter(i, 'mc_starfile', value)
        motion_error = max_motion_error(value, simulation_truth)
        write(logfhandle,'(a,i0,a,f7.3,a,f7.3)') '>>> TEST_PREPROC movie ', imovie, &
            &': motion max error ', motion_error, ' pixels; limit ', MOTION_TOL
        if( motion_error > MOTION_TOL )&
            &THROW_HARD('TEST_PREPROC FAILED: recovered frame motion differs from simulation truth by more than 0.50 pixels')
        call simulation_truth%kill

        call spproj%os_mic%getter(i, 'intg', value)
        call integrated%new(ldim, SMPD, wthreads=.false.)
        call optimal%new(ldim, SMPD, wthreads=.false.)
        call wrong_optimal%new(ldim, SMPD, wthreads=.false.)
        call integrated%read(value)
        call optimal%read(optimal_paths(imovie))
        image_corr = integrated%real_corr(optimal)
        wrong_movie = modulo(imovie, STREAM_NMOVS_SET) + 1
        call wrong_optimal%read(optimal_paths(wrong_movie))
        wrong_corr = integrated%real_corr(wrong_optimal)
        write(logfhandle,'(a,i0,a,f7.3,a,f7.3,a,f7.3)') '>>> TEST_PREPROC movie ', imovie, &
            &': correlation with truth ', image_corr, '; control ', wrong_corr, '; floor ', CORR_FLOOR
        if( image_corr < CORR_FLOOR )&
            &THROW_HARD('TEST_PREPROC FAILED: integrated micrograph correlation with truth is below 0.50')
        if( image_corr < wrong_corr + CORR_MARGIN )&
            &THROW_HARD('TEST_PREPROC FAILED: integrated micrograph does not identify its own simulation truth')
        call integrated%kill
        call optimal%kill
        call wrong_optimal%kill
    enddo
    call spproj%kill
    call cline_preproc%kill
    call simple_chdir(cwd_root, status)
    if( status /= 0 ) THROW_HARD('TEST_PREPROC FAILED: could not restore the original directory')
    write(logfhandle,'(a,i0,a,a)') 'PASS: stream_preproc validated ', STREAM_NMOVS_SET, &
        &' movies against simulation truth in ', fixture_root%to_char()
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

    integer function find_simulation( movie_path ) result(ind)
        type(string), intent(in) :: movie_path
        type(string) :: movie_basename, truth_basename
        integer :: isim
        ind = 0
        movie_basename = basename(movie_path)
        do isim = 1, STREAM_NMOVS_SET
            truth_basename = basename(movie_paths(isim))
            if( trim(movie_basename%to_char()) == trim(truth_basename%to_char()) )then
                ind = isim
                return
            endif
        enddo
    end function find_simulation

    real function max_motion_error( star_path, truth ) result(max_error)
        type(string), intent(in) :: star_path
        type(oris),   intent(in) :: truth
        type(starfile_table_type) :: table
        integer(C_long) :: object_id, num_objects
        real(C_double)  :: shift_x_dp, shift_y_dp
        real :: estimated(NFRAMES,2), expected(NFRAMES,2)
        integer :: iframe, fixed_frame
        logical :: got_x, got_y

        call starfile_table__new(table)
        call starfile_table__read(table, star_path, 'global_shift')
        object_id  = starfile_table__firstobject(table)
        num_objects = starfile_table__numberofobjects(table)
        if( int(num_objects - object_id) /= NFRAMES )&
            &THROW_HARD('TEST_PREPROC FAILED: motion STAR file has an unexpected number of frames')
        iframe = 0
        do while( object_id < num_objects .and. object_id >= 0_C_long )
            iframe = iframe + 1
            got_x = starfile_table__getValue_double(table, EMDL_MICROGRAPH_SHIFT_X, shift_x_dp)
            got_y = starfile_table__getValue_double(table, EMDL_MICROGRAPH_SHIFT_Y, shift_y_dp)
            if( .not. got_x .or. .not. got_y )&
                &THROW_HARD('TEST_PREPROC FAILED: motion STAR file is missing a frame shift')
            estimated(iframe,:) = [real(shift_x_dp), real(shift_y_dp)]
            expected(iframe,1) = truth%get(1, 'x'//int2str(iframe))
            expected(iframe,2) = truth%get(1, 'y'//int2str(iframe))
            object_id = starfile_table__nextobject(table)
        enddo
        call starfile_table__delete(table)
        fixed_frame = nint(real(NFRAMES) / 2.)
        estimated = estimated - spread(estimated(fixed_frame,:), 1, NFRAMES)
        expected  = expected  - spread(expected(fixed_frame,:),  1, NFRAMES)
        max_error = maxval(abs(estimated - expected))
    end function max_motion_error

end subroutine exec_test_stream_preproc

end module simple_commanders_test_stream
