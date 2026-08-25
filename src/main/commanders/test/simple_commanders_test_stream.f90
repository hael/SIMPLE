!@descr: tests for SIMPLE_stream workflows
module simple_commanders_test_stream
use simple_commanders_api
use simple_commanders_sim, only: commander_simulate_movie
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_abinitio2D_stream
  contains
    procedure :: execute => exec_test_abinitio2D_stream
end type commander_test_abinitio2D_stream

type, extends(commander_base) :: commander_test_assign_optics
  contains
    procedure :: execute => exec_test_assign_optics
end type commander_test_assign_optics

type, extends(commander_base) :: commander_test_gen_pickrefs
  contains
    procedure :: execute => exec_test_gen_pickrefs
end type commander_test_gen_pickrefs

type, extends(commander_base) :: commander_test_master
  contains
    procedure :: execute => exec_test_master
end type commander_test_master

type, extends(commander_base) :: commander_test_pick_extract
  contains
    procedure :: execute => exec_test_pick_extract
end type commander_test_pick_extract

type, extends(commander_base) :: commander_test_preproc
  contains
    procedure :: execute => exec_test_preproc
end type commander_test_preproc

type, extends(commander_base) :: commander_test_sieve_cavgs
  contains
    procedure :: execute => exec_test_sieve_cavgs
end type commander_test_sieve_cavgs

contains

subroutine exec_test_abinitio2D_stream( self, cline )
    class(commander_test_abinitio2D_stream), intent(inout) :: self
    class(cmdline),                         intent(inout) :: cline
    write(logfhandle,'(a)') '>>> TEST_ABINITIO2D_STREAM: DUMMY'
end subroutine exec_test_abinitio2D_stream

subroutine exec_test_assign_optics( self, cline )
    use simple_stream_p02_assign_optics_new, only: stream_p02_assign_optics
    use simple_ui,                          only: make_ui
    class(commander_test_assign_optics), intent(inout) :: self
    class(cmdline),                      intent(inout) :: cline
    character(len=*), parameter :: TEST_PROJFILE = 'test_assign_optics.simple'
    character(len=*), parameter :: TEST_OUTDIR   = 'stream_assign_optics'
    real,             parameter :: SMPD           = 1.3
    real,             parameter :: CS             = 2.7
    real,             parameter :: KV             = 300.0
    real,             parameter :: FRACA          = 0.1
    real,             parameter :: SHIFT_X(STREAM_NMOVS_SET) = [0.0, 0.1, 5.0, 5.1, 4.9]
    real,             parameter :: SHIFT_Y(STREAM_NMOVS_SET) = [0.0,-0.1, 5.0, 4.9, 5.1]
    integer,          parameter :: TILT_GROUP(STREAM_NMOVS_SET) = [1, 1, 2, 2, 2]
    type(stream_p02_assign_optics) :: xassign_optics
    type(cmdline)                  :: cline_assign
    type(sp_project)               :: upstream, result
    type(string)                   :: cwd_root, fixture_root, preproc_root, completed_dir
    type(string)                   :: batch_path, project_path, output_dir, fname
    integer                        :: i, status, group_a, group_b

    call make_ui
    call simple_getcwd(cwd_root)
    fixture_root = 'test_assign_optics_'//int2str(get_process_id())
    if( dir_exists(fixture_root) ) THROW_HARD('TEST_ASSIGN_OPTICS FAILED: fixture directory already exists')
    call simple_mkdir(fixture_root)
    call simple_chdir(fixture_root, status)
    if( status /= 0 ) THROW_HARD('TEST_ASSIGN_OPTICS FAILED: could not enter fixture directory')
    call simple_getcwd(fixture_root)

    ! Reproduce the two directories and completed five-micrograph project
    ! emitted by Stream preprocessing.
    preproc_root = filepath(fixture_root, 'preproc')
    call simple_mkdir(preproc_root)
    call simple_mkdir(filepath(preproc_root, 'spprojs'))
    completed_dir = filepath(preproc_root, 'spprojs_completed')
    call simple_mkdir(completed_dir)
    batch_path = filepath(completed_dir, '00001.simple')
    call upstream%os_mic%new(STREAM_NMOVS_SET, is_ptcl=.false.)
    do i = 1, STREAM_NMOVS_SET
        fname = 'synthetic_'//int2str_pad(i, 3)
        call upstream%os_mic%set(i, 'movie',       filepath(preproc_root, fname//'.mrc'))
        call upstream%os_mic%set(i, 'intg',        filepath(preproc_root, fname//'_intg.mrc'))
        call upstream%os_mic%set(i, 'mc_starfile', filepath(preproc_root, fname//'.star'))
        call upstream%os_mic%set(i, 'ctfjpg',      filepath(preproc_root, fname//'_ctf.jpg'))
        call upstream%os_mic%set(i, 'boxfile',     filepath(preproc_root, fname//'.box'))
        call upstream%os_mic%set(i, 'imgkind',     'mic')
        call upstream%os_mic%set(i, 'state',       1.0)
        call upstream%os_mic%set(i, 'importind',   real(i))
        call upstream%os_mic%set(i, 'xdim',        512.0)
        call upstream%os_mic%set(i, 'ydim',        512.0)
        call upstream%os_mic%set(i, 'nframes',     8.0)
        call upstream%os_mic%set(i, 'smpd',        SMPD)
        call upstream%os_mic%set(i, 'cs',          CS)
        call upstream%os_mic%set(i, 'kv',          KV)
        call upstream%os_mic%set(i, 'fraca',       FRACA)
        call upstream%os_mic%set(i, 'dfx',         1.5 + 0.1 * real(i))
        call upstream%os_mic%set(i, 'dfy',         1.6 + 0.1 * real(i))
        call upstream%os_mic%set(i, 'angast',      10.0 * real(i))
        call upstream%os_mic%set(i, 'phshift',     0.0)
        call upstream%os_mic%set(i, 'ctfres',      4.0)
        call upstream%os_mic%set(i, 'icefrac',     0.1)
        call upstream%os_mic%set(i, 'astig',       0.1)
        call upstream%os_mic%set(i, 'shiftx',      SHIFT_X(i))
        call upstream%os_mic%set(i, 'shifty',      SHIFT_Y(i))
        call upstream%os_mic%set(i, 'tiltgrp',     real(TILT_GROUP(i)))
    enddo
    call upstream%update_projinfo(batch_path)
    call upstream%write(batch_path)
    call upstream%kill

    project_path = filepath(fixture_root, TEST_PROJFILE)
    call cline_assign%set('prg',        'assign_optics')
    call cline_assign%set('projfile',   project_path)
    call cline_assign%set('dir_target', preproc_root)
    call cline_assign%set('outdir',     TEST_OUTDIR)
    call cline_assign%set('nthr',       1)
    call cline_assign%set('nmics',      STREAM_NMOVS_SET)
    call cline_assign%set('beamtilt',   'yes')
    call cline_assign%set('tilt_thres', 0.5)
    call xassign_optics%execute(cline_assign)

    project_path = cline_assign%get_carg('projfile')
    call simple_getcwd(output_dir)
    call simple_chdir(fixture_root, status)
    if( status /= 0 ) THROW_HARD('TEST_ASSIGN_OPTICS FAILED: could not leave execution directory')
    if( .not. file_exists(project_path) ) THROW_HARD('TEST_ASSIGN_OPTICS FAILED: output project was not created')
    call result%read(project_path)
    if( result%os_mic%get_noris() /= STREAM_NMOVS_SET )&
        &THROW_HARD('TEST_ASSIGN_OPTICS FAILED: project does not contain five micrographs')
    if( result%os_optics%get_noris() /= 2 )&
        &THROW_HARD('TEST_ASSIGN_OPTICS FAILED: expected two optics groups')
    group_a = result%os_mic%get_int(1, 'ogid')
    group_b = result%os_mic%get_int(3, 'ogid')
    if( group_a == group_b .or. group_a < 1 .or. group_a > 2 .or. group_b < 1 .or. group_b > 2 )&
        &THROW_HARD('TEST_ASSIGN_OPTICS FAILED: controlled shift clusters were not separated')
    if( result%os_mic%get_int(2, 'ogid') /= group_a .or. &
        &any([result%os_mic%get_int(4, 'ogid'), result%os_mic%get_int(5, 'ogid')] /= group_b) )&
        &THROW_HARD('TEST_ASSIGN_OPTICS FAILED: micrograph optics assignments are inconsistent')
    if( nint(result%os_optics%get(group_a, 'pop')) /= 2 .or. &
        &nint(result%os_optics%get(group_b, 'pop')) /= 3 )&
        &THROW_HARD('TEST_ASSIGN_OPTICS FAILED: optics-group populations are incorrect')
    if( abs(result%os_optics%get(group_a, 'opcx') - 0.05) > 0.01 .or. &
        &abs(result%os_optics%get(group_a, 'opcy') + 0.05) > 0.01 .or. &
        &abs(result%os_optics%get(group_b, 'opcx') - 5.00) > 0.01 .or. &
        &abs(result%os_optics%get(group_b, 'opcy') - 5.00) > 0.01 )&
        &THROW_HARD('TEST_ASSIGN_OPTICS FAILED: optics-group centroids are incorrect')
    call result%kill

    if( .not. file_exists(filepath(output_dir, 'optics.star')) )&
        &THROW_HARD('TEST_ASSIGN_OPTICS FAILED: optics STAR file was not created')
    if( .not. file_exists(filepath(output_dir, 'micrographs.star')) )&
        &THROW_HARD('TEST_ASSIGN_OPTICS FAILED: micrographs STAR file was not created')
    if( .not. file_exists(filepath(output_dir, OPTICS_MAP_PREFIX//'1'//TXT_EXT)) )&
        &THROW_HARD('TEST_ASSIGN_OPTICS FAILED: optics map table was not created')
    if( .not. file_exists(filepath(output_dir, OPTICS_MAP_PREFIX//'1'//METADATA_EXT)) )&
        &THROW_HARD('TEST_ASSIGN_OPTICS FAILED: optics map project was not created')
    call cline_assign%kill
    call simple_chdir(cwd_root, status)
    if( status /= 0 ) THROW_HARD('TEST_ASSIGN_OPTICS FAILED: could not restore the original directory')
    write(*,'(a,a)') '>>> TEST_ASSIGN_OPTICS: validated two optics groups in ', fixture_root%to_char(); call flush(6)
    call simple_end('**** SIMPLE_TEST_ASSIGN_OPTICS NORMAL STOP ****')
end subroutine exec_test_assign_optics

subroutine exec_test_gen_pickrefs( self, cline )
    use simple_commanders_pick, only: commander_make_pickrefs
    use simple_ui,              only: make_ui
    class(commander_test_gen_pickrefs), intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    character(len=*), parameter :: SOURCE_STACK    = 'source_references.mrcs'
    character(len=*), parameter :: TEST_OUTDIR     = 'stream_gen_pickrefs'
    character(len=*), parameter :: OUTPUT_STACK    = PICKREFS_FBODY//MRC_EXT
    character(len=*), parameter :: SOURCE_PREVIEW  = 'pickrefs_source.jpeg'
    character(len=*), parameter :: FINISHED_FILE   = 'MAKE_PICKREFS_FINISHED'
    real,             parameter :: SMPD             = 2.0
    integer,          parameter :: SOURCE_BOX       = 64
    integer,          parameter :: NSOURCE_REFS     = 3
    integer,          parameter :: NROTS            = 4
    type(commander_make_pickrefs) :: xmake_pickrefs
    type(cmdline)                  :: cline_make_pickrefs
    type(image)                    :: reference, feature
    type(oris)                     :: moldiam_meta
    type(string)                   :: cwd_root, fixture_root, source_path, output_dir
    type(string)                   :: output_path, moldiam_path
    integer                        :: i, status, ldim(3), nrefs
    real                           :: output_smpd, diam_max, mskdiam
    integer                        :: box_for_pick, box_for_extract

    call make_ui
    call simple_getcwd(cwd_root)
    fixture_root = 'test_gen_pickrefs_'//int2str(get_process_id())
    if( dir_exists(fixture_root) ) THROW_HARD('TEST_GEN_PICKREFS FAILED: fixture directory already exists')
    call simple_mkdir(fixture_root)
    call simple_chdir(fixture_root, status)
    if( status /= 0 ) THROW_HARD('TEST_GEN_PICKREFS FAILED: could not enter fixture directory')
    call simple_getcwd(fixture_root)

    ! Construct three compact asymmetric class-average-like images.  The
    ! production make_pickrefs path automasks, centers, filters, rotates,
    ! mirrors, and clips these images into the final picking-reference stack.
    do i = 1, NSOURCE_REFS
        call reference%new([SOURCE_BOX, SOURCE_BOX, 1], SMPD, wthreads=.false.)
        call reference%square(5 + i)
        call feature%new([SOURCE_BOX, SOURCE_BOX, 1], SMPD, wthreads=.false.)
        call feature%square(2 + i)
        call feature%shift([real(7 + 2 * i), real(-5 + i), 0.])
        call reference%add(feature)
        call reference%write(string(SOURCE_STACK), i, del_if_exists=(i == 1))
        call feature%kill
        call reference%kill
    enddo
    source_path = simple_abspath(string(SOURCE_STACK))

    call cline_make_pickrefs%set('prg',          'make_pickrefs')
    call cline_make_pickrefs%set('pickrefs',     source_path)
    call cline_make_pickrefs%set('smpd',         SMPD)
    call cline_make_pickrefs%set('nrots',        NROTS)
    call cline_make_pickrefs%set('mirr',         'yes')
    call cline_make_pickrefs%set('trust_header', 'yes')
    call cline_make_pickrefs%set('ncls',         0)
    call cline_make_pickrefs%set('nthr',         1)
    call cline_make_pickrefs%set('mkdir',        'yes')
    call cline_make_pickrefs%set('dir_exec',     TEST_OUTDIR)
    call cline_make_pickrefs%printline(unit=6)
    call xmake_pickrefs%execute(cline_make_pickrefs)

    call simple_getcwd(output_dir)
    output_path  = filepath(output_dir, OUTPUT_STACK)
    moldiam_path = filepath(output_dir, STREAM_MOLDIAM)
    if( .not. file_exists(output_path) ) THROW_HARD('TEST_GEN_PICKREFS FAILED: output stack was not created')
    call find_ldim_nptcls(output_path, ldim, nrefs)
    if( nrefs /= NSOURCE_REFS * NROTS * 2 )&
        &THROW_HARD('TEST_GEN_PICKREFS FAILED: rotation and mirror expansion count is incorrect')
    if( ldim(1) <= 0 .or. ldim(1) /= ldim(2) .or. ldim(1) > SOURCE_BOX )&
        &THROW_HARD('TEST_GEN_PICKREFS FAILED: output reference dimensions are invalid')
    output_smpd = find_img_smpd(output_path)
    if( abs(output_smpd - SMPD) > 0.01 )&
        &THROW_HARD('TEST_GEN_PICKREFS FAILED: output sampling distance changed unexpectedly')
    if( .not. file_exists(moldiam_path) ) THROW_HARD('TEST_GEN_PICKREFS FAILED: diameter metadata was not created')
    call moldiam_meta%new(1, .false.)
    call moldiam_meta%read(moldiam_path)
    if( .not. moldiam_meta%isthere(1, 'diam_max') .or. .not. moldiam_meta%isthere(1, 'mskdiam') .or. &
        &.not. moldiam_meta%isthere(1, 'box_for_pick') .or. .not. moldiam_meta%isthere(1, 'box_for_extract') )&
        &THROW_HARD('TEST_GEN_PICKREFS FAILED: required diameter metadata fields are missing')
    diam_max        = moldiam_meta%get(1, 'diam_max')
    mskdiam         = moldiam_meta%get(1, 'mskdiam')
    box_for_pick    = moldiam_meta%get_int(1, 'box_for_pick')
    box_for_extract = moldiam_meta%get_int(1, 'box_for_extract')
    if( diam_max <= 0. .or. mskdiam <= 0. .or. box_for_pick /= ldim(1) .or. box_for_extract < box_for_pick )&
        &THROW_HARD('TEST_GEN_PICKREFS FAILED: diameter metadata values are invalid')
    call moldiam_meta%kill
    if( .not. file_exists(filepath(output_dir, SOURCE_PREVIEW)) )&
        &THROW_HARD('TEST_GEN_PICKREFS FAILED: source preview was not created')
    if( .not. file_exists(filepath(output_dir, FINISHED_FILE)) )&
        &THROW_HARD('TEST_GEN_PICKREFS FAILED: completion marker was not created')

    call cline_make_pickrefs%kill
    call simple_chdir(cwd_root, status)
    if( status /= 0 ) THROW_HARD('TEST_GEN_PICKREFS FAILED: could not restore the original directory')
    write(*,'(a,i0,a,a)') '>>> TEST_GEN_PICKREFS: validated ', nrefs, ' picking references in ', output_dir%to_char(); call flush(6)
    call simple_end('**** SIMPLE_TEST_GEN_PICKREFS NORMAL STOP ****')
end subroutine exec_test_gen_pickrefs

subroutine exec_test_master( self, cline )
    class(commander_test_master), intent(inout) :: self
    class(cmdline),                intent(inout) :: cline
    write(logfhandle,'(a)') '>>> TEST_MASTER: DUMMY'
end subroutine exec_test_master

subroutine exec_test_pick_extract( self, cline )
    class(commander_test_pick_extract), intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    write(logfhandle,'(a)') '>>> TEST_PICK_EXTRACT: DUMMY'
end subroutine exec_test_pick_extract

subroutine exec_test_preproc( self, cline )
    use simple_stream_p01_preprocess_new, only: stream_p01_preprocess
    use simple_ui,                        only: make_ui
    use, intrinsic :: ieee_arithmetic,    only: ieee_is_finite
    class(commander_test_preproc), intent(inout) :: self
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

end subroutine exec_test_preproc

subroutine exec_test_sieve_cavgs( self, cline )
    class(commander_test_sieve_cavgs), intent(inout) :: self
    class(cmdline),                    intent(inout) :: cline
    write(logfhandle,'(a)') '>>> TEST_SIEVE_CAVGS: DUMMY'
end subroutine exec_test_sieve_cavgs

end module simple_commanders_test_stream
