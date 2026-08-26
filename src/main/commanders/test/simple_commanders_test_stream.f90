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
    use simple_commanders_abinitio2D,   only: commander_abinitio2D
    use simple_commanders_project_ptcl, only: commander_import_particles
    use simple_defs_fname,              only: ABINITIO2D_FINISHED, METADATA_EXT, MRC_EXT
    use simple_ui,                      only: make_ui
    class(commander_test_abinitio2D_stream), intent(inout) :: self
    class(cmdline),                         intent(inout) :: cline
    character(len=*), parameter :: PROJNAME       = 'stream_abinitio2D_chunk'
    character(len=*), parameter :: PROJFILE       = PROJNAME//METADATA_EXT
    character(len=*), parameter :: PARTICLE_STACK = 'stream_particles'//MRC_EXT
    real,             parameter :: SMPD            = 2.0
    real,             parameter :: MSKDIAM         = 32.0
    integer,          parameter :: BOX             = 32
    integer,          parameter :: NCLASSES        = 2
    integer,          parameter :: NPARTICLES      = 24
    integer,          parameter :: NPER_CLASS      = NPARTICLES / NCLASSES
    type(commander_import_particles) :: ximport_particles
    type(commander_abinitio2D)       :: xabinitio2D
    type(cmdline)                     :: cline_import, cline_abinitio2D
    type(image)                       :: particle, feature
    type(sp_project)                  :: project, result
    type(string)                      :: cwd_root, fixture_root, project_path, stack_path, cavgs_path
    integer                           :: i, icls, status, ldim(3), nclasses_out, nimages
    real                              :: cavgs_smpd

    call make_ui
    call simple_getcwd(cwd_root)
    fixture_root = 'test_abinitio2D_stream_'//int2str(get_process_id())
    if( dir_exists(fixture_root) ) THROW_HARD('TEST_ABINITIO2D_STREAM FAILED: fixture directory already exists')
    call simple_mkdir(fixture_root)
    call simple_chdir(fixture_root, status)
    if( status /= 0 ) THROW_HARD('TEST_ABINITIO2D_STREAM FAILED: could not enter fixture directory')
    call simple_getcwd(fixture_root)

    ! Build two asymmetric particle families representative of one completed
    ! stream-extraction chunk.  Each family contains a centered component and
    ! a different off-axis feature so that the images contain real 2D signal.
    do i = 1, NPARTICLES
        icls = 1 + (i - 1) / NPER_CLASS
        call particle%new([BOX, BOX, 1], SMPD, wthreads=.false.)
        call particle%zero()
        call feature%new([BOX, BOX, 1], SMPD, wthreads=.false.)
        call feature%square(5)
        call feature%mul(2.)
        call particle%add(feature)
        call feature%kill()
        call feature%new([BOX, BOX, 1], SMPD, wthreads=.false.)
        if( icls == 1 )then
            call feature%square(2)
            call feature%shift([7., -4., 0.])
        else
            call feature%square(3)
            call feature%shift([-6., 5., 0.])
        endif
        call particle%add(feature)
        call particle%write(string(PARTICLE_STACK), i, del_if_exists=(i == 1))
        call feature%kill()
        call particle%kill()
    enddo
    stack_path   = simple_abspath(string(PARTICLE_STACK))
    project_path = filepath(fixture_root, PROJFILE)
    call project%update_projinfo(project_path)
    call project%write(project_path)
    call project%kill()

    ! Import exactly as a stream chunk project presents extracted particles to
    ! its finite abinitio2D child command.
    call cline_import%set('prg',      'import_particles')
    call cline_import%set('projfile', project_path)
    call cline_import%set('stk',      stack_path)
    call cline_import%set('smpd',     SMPD)
    call cline_import%set('ctf',      'no')
    call cline_import%set('mkdir',    'no')
    call cline_import%printline(unit=6)
    call ximport_particles%execute(cline_import)

    ! The top-level stream commander is a watcher and cannot terminate inside
    ! a unit test.  Its chunk analysis dispatches this production commander.
    ! One stage and one iteration keep the behavioral test compact.
    call cline_abinitio2D%set('prg',            'abinitio2D')
    call cline_abinitio2D%set('projfile',       project_path)
    call cline_abinitio2D%set('ncls',           NCLASSES)
    call cline_abinitio2D%set('mskdiam',        MSKDIAM)
    call cline_abinitio2D%set('nthr',           1)
    call cline_abinitio2D%set('mkdir',          'no')
    call cline_abinitio2D%set('autoscale',      'no')
    call cline_abinitio2D%set('center',         'no')
    call cline_abinitio2D%set('cls_init',       'ptcl')
    call cline_abinitio2D%set('refine',         'snhc_smpl')
    call cline_abinitio2D%set('inpl_cont',      'no')
    call cline_abinitio2D%set('sigma_est',      'global')
    call cline_abinitio2D%set('ml_reg',         'no')
    call cline_abinitio2D%set('lp',             20.0)
    call cline_abinitio2D%set('extr_lim',       4)
    call cline_abinitio2D%set('nits_per_stage', 1)
    call cline_abinitio2D%set('nstages',        1)
    call cline_abinitio2D%set('rank_cavgs',     'no')
    call cline_abinitio2D%printline(unit=6)
    call xabinitio2D%execute(cline_abinitio2D)

    if( .not. file_exists(filepath(fixture_root, ABINITIO2D_FINISHED)) )&
        &THROW_HARD('TEST_ABINITIO2D_STREAM FAILED: completion marker was not created')
    call result%read(project_path)
    if( result%get_nptcls() /= NPARTICLES )&
        &THROW_HARD('TEST_ABINITIO2D_STREAM FAILED: particle count changed during classification')
    if( result%os_cls2D%get_noris() /= NCLASSES )&
        &THROW_HARD('TEST_ABINITIO2D_STREAM FAILED: expected two output classes')
    do i = 1, NPARTICLES
        icls = result%os_ptcl2D%get_class(i)
        if( icls < 1 .or. icls > NCLASSES )&
            &THROW_HARD('TEST_ABINITIO2D_STREAM FAILED: particle lacks a valid class assignment')
    enddo
    if( nint(sum(result%os_cls2D%get_all('pop'))) /= NPARTICLES )&
        &THROW_HARD('TEST_ABINITIO2D_STREAM FAILED: class populations do not cover all particles')
    call result%get_cavgs_stk(cavgs_path, nclasses_out, cavgs_smpd, imgkind='cavg')
    if( nclasses_out /= NCLASSES .or. abs(cavgs_smpd - SMPD) > 0.01 )&
        &THROW_HARD('TEST_ABINITIO2D_STREAM FAILED: class-average metadata is incorrect')
    if( .not. file_exists(cavgs_path) )&
        &THROW_HARD('TEST_ABINITIO2D_STREAM FAILED: final class-average stack was not created')
    call find_ldim_nptcls(cavgs_path, ldim, nimages)
    if( nimages /= NCLASSES .or. any(ldim(1:2) /= [BOX, BOX]) )&
        &THROW_HARD('TEST_ABINITIO2D_STREAM FAILED: final class-average stack is invalid')
    call result%kill()

    call cline_import%kill()
    call cline_abinitio2D%kill()
    call simple_chdir(cwd_root, status)
    if( status /= 0 ) THROW_HARD('TEST_ABINITIO2D_STREAM FAILED: could not restore the original directory')
    write(*,'(a,a)') '>>> TEST_ABINITIO2D_STREAM: validated two class averages and all particle assignments in ', &
        &fixture_root%to_char()
    call flush(6)
    call simple_end('**** SIMPLE_TEST_ABINITIO2D_STREAM NORMAL STOP ****')
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
    use unix,                  only: c_usleep
    use simple_forked_process, only: forked_process, FORK_POLL_TIME
    use simple_gui_assembler,  only: gui_assembler
    use simple_gui_metadata_api, only: CK, json_core, json_value
    class(commander_test_master), intent(inout) :: self
    class(cmdline),                intent(inout) :: cline
    integer,          parameter :: TEST_JOB_ID = 42
    integer,          parameter :: NSTAGES     = 8
    character(len=24), parameter :: STAGE_NAMES(NSTAGES) = [character(len=24) :: &
        &'preprocessing', 'assign_optics', 'initial_picking', 'opening2D', &
        &'reference_picking', 'particle_sieving', 'pool2D', 'master']
    type(forked_process) :: fork_preprocess, fork_assign_optics, fork_opening2D
    type(forked_process) :: fork_reference_picking, fork_particle_sieving, fork_pool2D
    type(gui_assembler)  :: assembler
    type(string)         :: running_heartbeat, finished_heartbeat
    integer              :: rc

#if defined(_WIN32)
    write(*,'(a)') '>>> TEST_MASTER: skipped because forked processes are unavailable on Windows'
#else
    write(*,'(a,i0)') '>>> TEST_MASTER JOB ID: ', TEST_JOB_ID
    write(*,'(a,i0)') '>>> TEST_MASTER HEARTBEAT ENTRIES: ', NSTAGES

    ! The production master reports six forked workers as seven GUI stages:
    ! initial_picking and opening2D deliberately share the opening2D process.
    ! Use the finite default fork worker so the orchestration lifecycle can be
    ! exercised without starting the persistent Stream pipeline or NICE.
    call assembler%new(TEST_JOB_ID)
    call fork_preprocess%start(        name=string('TEST_MASTER_PREPROCESS'))
    call fork_assign_optics%start(      name=string('TEST_MASTER_ASSIGN_OPTICS'))
    call fork_opening2D%start(           name=string('TEST_MASTER_OPENING2D'))
    call fork_reference_picking%start(   name=string('TEST_MASTER_REFERENCE_PICKING'))
    call fork_particle_sieving%start(    name=string('TEST_MASTER_PARTICLE_SIEVING'))
    call fork_pool2D%start(               name=string('TEST_MASTER_POOL2D'))
    rc = c_usleep(FORK_POLL_TIME * 5)

    call assembler%assemble_stream_heartbeat(fork_preprocess, fork_assign_optics, fork_opening2D, &
        &fork_reference_picking, fork_particle_sieving, fork_pool2D)
    running_heartbeat = assembler%to_string()

    call fork_preprocess%terminate()
    call fork_assign_optics%terminate()
    call fork_opening2D%terminate()
    call fork_reference_picking%terminate()
    call fork_particle_sieving%terminate()
    call fork_pool2D%terminate()
    call fork_preprocess%await_final_status()
    call fork_assign_optics%await_final_status()
    call fork_opening2D%await_final_status()
    call fork_reference_picking%await_final_status()
    call fork_particle_sieving%await_final_status()
    call fork_pool2D%await_final_status()

    call assembler%set_stoptime()
    call assembler%assemble_stream_heartbeat(fork_preprocess, fork_assign_optics, fork_opening2D, &
        &fork_reference_picking, fork_particle_sieving, fork_pool2D)
    finished_heartbeat = assembler%to_string()
    call assembler%kill()

    if( .not. heartbeat_matches(running_heartbeat, 'running') )&
        &THROW_HARD('TEST_MASTER FAILED: running heartbeat is incomplete or incorrect')
    if( .not. heartbeat_matches(finished_heartbeat, 'finished') )&
        &THROW_HARD('TEST_MASTER FAILED: finished heartbeat is incomplete or incorrect')

    write(*,'(a)') '>>> TEST_MASTER: validated running and finished lifecycle heartbeats for all Stream stages'
#endif
    call flush(6)
    call simple_end('**** SIMPLE_TEST_MASTER NORMAL STOP ****')

contains

    function heartbeat_matches( payload, expected_status ) result( matches )
        type(string),     intent(in) :: payload
        character(len=*), intent(in) :: expected_status
        type(json_core)              :: json
        type(json_value), pointer    :: root
        character(kind=CK,len=:), allocatable :: actual_status
        character(len=:),         allocatable :: path
        integer :: i, job_id, pid, starttime, stoptime
        logical :: found, matches

        matches = .false.
        nullify(root)
        call json%initialize()
        call json%parse(root, payload%to_char())
        if( json%failed() .or. .not. associated(root) ) goto 100
        call json%get(root, 'jobid', job_id, found)
        if( .not. found ) goto 100
        if( job_id /= TEST_JOB_ID ) goto 100

        do i = 1, NSTAGES
            path = 'stream_heartbeat.'//trim(STAGE_NAMES(i))
            call json%get(root, path//'.status', actual_status, found)
            if( .not. found ) goto 100
            if( actual_status /= expected_status ) goto 100
            call json%get(root, path//'.pid', pid, found)
            if( .not. found ) goto 100
            if( pid <= 0 ) goto 100
            call json%get(root, path//'.starttime', starttime, found)
            if( .not. found ) goto 100
            if( starttime <= 0 ) goto 100
            call json%get(root, path//'.stoptime', stoptime, found)
            if( .not. found ) goto 100
            if( expected_status == 'running' .and. stoptime /= 0 ) goto 100
            if( expected_status == 'finished' .and. stoptime <= 0 ) goto 100
        enddo
        matches = .true.

100     if( associated(root) ) call json%destroy(root)
    end function heartbeat_matches
end subroutine exec_test_master

subroutine exec_test_pick_extract( self, cline )
    use simple_commanders_pick, only: commander_pick_extract
    use simple_ui,              only: make_ui
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    class(commander_test_pick_extract), intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    character(len=*), parameter :: MICROGRAPH_FILE = 'synthetic_micrograph.mrc'
    character(len=*), parameter :: PICKREFS_FILE   = 'synthetic_pickrefs.mrcs'
    character(len=*), parameter :: TEST_PROJFILE   = 'test_pick_extract.simple'
    real,             parameter :: SMPD             = 2.0
    integer,          parameter :: MICROGRAPH_BOX   = 256
    integer,          parameter :: PICKREF_BOX      = 64
    integer,          parameter :: EXTRACT_BOX      = 64
    integer,          parameter :: NPARTICLES       = 3
    integer,          parameter :: PARTICLE_COORDS(2, NPARTICLES) = reshape(&
        &[32, 32, 160, 40, 96, 152], [2, NPARTICLES])
    type(commander_pick_extract) :: xpick_extract
    type(cmdline)                :: cline_pick_extract
    type(image)                  :: micrograph, reference, feature, extracted
    type(sp_project)             :: fixture, result
    type(string)                 :: cwd_root, fixture_root, micrograph_path, pickrefs_path
    type(string)                 :: project_path, boxfile_path, thumbnail_path, stack_path
    integer                      :: i, status, npicked, nimages, ldim(3)
    real                         :: output_smpd, extracted_variance

    call make_ui
    call simple_getcwd(cwd_root)
    fixture_root = 'test_pick_extract_'//int2str(get_process_id())
    if( dir_exists(fixture_root) ) THROW_HARD('TEST_PICK_EXTRACT FAILED: fixture directory already exists')
    call simple_mkdir(fixture_root)
    call simple_chdir(fixture_root, status)
    if( status /= 0 ) THROW_HARD('TEST_PICK_EXTRACT FAILED: could not enter fixture directory')
    call simple_getcwd(fixture_root)

    ! Build an asymmetric reference and place three exact copies on a weak-noise
    ! micrograph.  This gives the reference picker controlled, unambiguous peaks.
    call reference%new([PICKREF_BOX, PICKREF_BOX, 1], SMPD, wthreads=.false.)
    call reference%square(7)
    call feature%new([PICKREF_BOX, PICKREF_BOX, 1], SMPD, wthreads=.false.)
    call feature%square(3)
    call feature%shift([11., -7., 0.])
    call reference%add(feature)
    call reference%mul(5.)
    call reference%write(string(PICKREFS_FILE), 1, del_if_exists=.true.)
    call micrograph%new([MICROGRAPH_BOX, MICROGRAPH_BOX, 1], SMPD, wthreads=.false.)
    call micrograph%gauran(0., 0.02)
    do i = 1, NPARTICLES
        call micrograph%add_window(reference, PARTICLE_COORDS(:, i))
    enddo
    call micrograph%write(string(MICROGRAPH_FILE), 1, del_if_exists=.true.)
    call feature%kill
    call reference%kill
    call micrograph%kill
    micrograph_path = simple_abspath(string(MICROGRAPH_FILE))
    pickrefs_path   = simple_abspath(string(PICKREFS_FILE))

    ! Reproduce the single-micrograph project handed to a Stream pick_extract job.
    project_path = filepath(fixture_root, TEST_PROJFILE)
    call fixture%os_mic%new(1, is_ptcl=.false.)
    call fixture%os_mic%set(1, 'movie',     micrograph_path)
    call fixture%os_mic%set(1, 'intg',      micrograph_path)
    call fixture%os_mic%set(1, 'imgkind',   'mic')
    call fixture%os_mic%set(1, 'state',     1.0)
    call fixture%os_mic%set(1, 'importind', 1.0)
    call fixture%os_mic%set(1, 'xdim',      real(MICROGRAPH_BOX))
    call fixture%os_mic%set(1, 'ydim',      real(MICROGRAPH_BOX))
    call fixture%os_mic%set(1, 'nframes',   1.0)
    call fixture%os_mic%set(1, 'smpd',      SMPD)
    call fixture%update_projinfo(project_path)
    call fixture%write(project_path)
    call fixture%kill

    call cline_pick_extract%set('prg',             'pick_extract')
    call cline_pick_extract%set('projfile',         project_path)
    call cline_pick_extract%set('pickrefs',         pickrefs_path)
    call cline_pick_extract%set('dir',              fixture_root)
    call cline_pick_extract%set('stream',           'yes')
    call cline_pick_extract%set('extract',          'yes')
    call cline_pick_extract%set('picker',           'new')
    call cline_pick_extract%set('pcontrast',        'white')
    call cline_pick_extract%set('pick_roi',         'no')
    call cline_pick_extract%set('backgr_subtr',     'no')
    call cline_pick_extract%set('refpick_backend',  'legacy')
    call cline_pick_extract%set('nboxes_max',       NPARTICLES)
    call cline_pick_extract%set('box_extract',      EXTRACT_BOX)
    call cline_pick_extract%set('fromp',            1)
    call cline_pick_extract%set('top',              1)
    call cline_pick_extract%set('nthr',             1)
    call cline_pick_extract%set('mkdir',            'no')
    call cline_pick_extract%printline(unit=6)
    call xpick_extract%execute(cline_pick_extract)

    if( .not. file_exists(project_path) ) THROW_HARD('TEST_PICK_EXTRACT FAILED: output project was not created')
    call result%read(project_path)
    if( result%os_mic%get_noris() /= 1 )&
        &THROW_HARD('TEST_PICK_EXTRACT FAILED: output project does not contain one micrograph')
    if( result%os_mic%get_state(1) /= 1 )&
        &THROW_HARD('TEST_PICK_EXTRACT FAILED: controlled micrograph was rejected')
    npicked = result%os_mic%get_int(1, 'nptcls')
    if( npicked /= NPARTICLES )&
        &THROW_HARD('TEST_PICK_EXTRACT FAILED: expected three picked particles')
    boxfile_path  = result%os_mic%get_str(1, 'boxfile')
    thumbnail_path = result%os_mic%get_str(1, 'thumb_den')
    if( .not. file_exists(boxfile_path) ) THROW_HARD('TEST_PICK_EXTRACT FAILED: box file was not created')
    if( .not. file_exists(thumbnail_path) ) THROW_HARD('TEST_PICK_EXTRACT FAILED: density thumbnail was not created')
    if( result%os_ptcl2D%get_noris() /= npicked .or. result%os_ptcl3D%get_noris() /= npicked )&
        &THROW_HARD('TEST_PICK_EXTRACT FAILED: particle metadata count does not match picking output')
    if( result%os_stk%get_noris() /= 1 ) THROW_HARD('TEST_PICK_EXTRACT FAILED: extraction stack metadata is missing')
    stack_path = result%os_stk%get_str(1, 'stk')
    if( .not. file_exists(stack_path) ) THROW_HARD('TEST_PICK_EXTRACT FAILED: extracted particle stack was not created')
    call find_ldim_nptcls(stack_path, ldim, nimages)
    if( nimages /= npicked ) THROW_HARD('TEST_PICK_EXTRACT FAILED: extracted image count does not match picking output')
    if( any(ldim(1:2) /= [EXTRACT_BOX, EXTRACT_BOX]) )&
        &THROW_HARD('TEST_PICK_EXTRACT FAILED: extracted particle dimensions are incorrect')
    output_smpd = find_img_smpd(stack_path)
    if( abs(output_smpd - SMPD) > 0.01 ) THROW_HARD('TEST_PICK_EXTRACT FAILED: extraction sampling distance changed')
    call extracted%new([EXTRACT_BOX, EXTRACT_BOX, 1], SMPD, wthreads=.false.)
    call extracted%read(stack_path, 1)
    extracted_variance = extracted%variance()
    if( .not. ieee_is_finite(extracted_variance) .or. extracted_variance <= TINY )&
        &THROW_HARD('TEST_PICK_EXTRACT FAILED: extracted particle data are invalid')
    call extracted%kill
    call result%kill

    call cline_pick_extract%kill
    call simple_chdir(cwd_root, status)
    if( status /= 0 ) THROW_HARD('TEST_PICK_EXTRACT FAILED: could not restore the original directory')
    write(*,'(a,i0,a,a)') '>>> TEST_PICK_EXTRACT: validated ', npicked, ' extracted particles in ', fixture_root%to_char()
    call flush(6)
    call simple_end('**** SIMPLE_TEST_PICK_EXTRACT NORMAL STOP ****')
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
    use simple_defs_fname, only: ABINITIO2D_FINISHED, JPG_EXT, METADATA_EXT, MRC_EXT
    use simple_fileio,     only: simple_touch
    use simple_ptcl_sieve, only: ptcl_sieve
    use simple_ui,         only: make_ui
    class(commander_test_sieve_cavgs), intent(inout) :: self
    class(cmdline),                    intent(inout) :: cline
    character(len=*), parameter :: CHUNK_TIER_DIR = 'chunks_coarse'
    character(len=*), parameter :: CHUNK_STEM     = 'chunk_coarse_1'
    character(len=*), parameter :: CAVG_STACK     = 'cavgs_iter001'//MRC_EXT
    character(len=*), parameter :: COMPLETED_DIR  = 'completed'
    real,             parameter :: SMPD            = 2.0
    integer,          parameter :: CAVG_BOX        = 64
    integer,          parameter :: NCLASSES        = 2
    integer,          parameter :: NPARTICLES      = 8
    integer,          parameter :: NPARTICLES_PER_CLASS = NPARTICLES / NCLASSES
    type(ptcl_sieve)              :: sieve
    type(parameters)              :: params_sieve
    type(cmdline)                 :: cline_sieve
    type(image)                   :: cavg_good, cavg_bad, feature
    type(sp_project)              :: chunk_project, result
    type(string)                  :: cwd_root, fixture_root, chunk_dir, completed_path
    type(string)                  :: chunk_projfile, cavg_path, completed_projfile
    type(string)                  :: latest_jpeg, latest_stk, rejection_reason
    type(string)                  :: selected_jpeg, rejected_jpeg, reasons_jpeg, reasons_key
    integer,          allocatable :: latest_inds(:), latest_pops(:), latest_selection(:)
    real,             allocatable :: latest_res(:)
    integer                       :: i, icls, status, ldim(3), nimages, xtiles, ytiles
    logical                       :: has_latest

    call make_ui
    call simple_getcwd(cwd_root)
    fixture_root = 'test_sieve_cavgs_'//int2str(get_process_id())
    if( dir_exists(fixture_root) ) THROW_HARD('TEST_SIEVE_CAVGS FAILED: fixture directory already exists')
    call simple_mkdir(fixture_root)
    call simple_chdir(fixture_root, status)
    if( status /= 0 ) THROW_HARD('TEST_SIEVE_CAVGS FAILED: could not enter fixture directory')
    call simple_getcwd(fixture_root)

    chunk_dir = filepath(fixture_root, CHUNK_TIER_DIR)
    call simple_mkdir(chunk_dir)
    chunk_dir = filepath(chunk_dir, CHUNK_STEM)
    call simple_mkdir(chunk_dir)
    completed_path = filepath(fixture_root, COMPLETED_DIR)
    call simple_mkdir(completed_path)
    chunk_projfile = filepath(chunk_dir, CHUNK_STEM//METADATA_EXT)
    cavg_path       = filepath(chunk_dir, CAVG_STACK)

    ! Class 1 contains a strong centered component and should survive the
    ! hard gates.  Class 2 is blank and must be rejected as NO_COMPONENT.
    call cavg_good%new([CAVG_BOX, CAVG_BOX, 1], SMPD, wthreads=.false.)
    call cavg_good%gauran(0., 0.02)
    call feature%new([CAVG_BOX, CAVG_BOX, 1], SMPD, wthreads=.false.)
    call feature%square(8)
    call feature%mul(5.)
    call cavg_good%add(feature)
    call cavg_good%write(cavg_path, 1, del_if_exists=.true.)
    call cavg_bad%new([CAVG_BOX, CAVG_BOX, 1], SMPD, wthreads=.false.)
    call cavg_bad%zero()
    call cavg_bad%write(cavg_path, 2)
    call feature%kill()
    call cavg_good%kill()
    call cavg_bad%kill()

    ! Reproduce a completed coarse abinitio2D chunk awaiting sieve rejection.
    call chunk_project%os_mic%new(1, is_ptcl=.false.)
    call chunk_project%os_mic%set_state(1, 1)
    call chunk_project%os_mic%set(1, 'imgkind', 'mic')
    call chunk_project%os_mic%set(1, 'nptcls', NPARTICLES)
    call chunk_project%os_mic%set(1, 'smpd', SMPD)
    call chunk_project%os_ptcl2D%new(NPARTICLES, is_ptcl=.true.)
    do i = 1, NPARTICLES
        icls = 1 + (i - 1) / NPARTICLES_PER_CLASS
        call chunk_project%os_ptcl2D%set_class(i, icls)
        call chunk_project%os_ptcl2D%set_state(i, 1)
        call chunk_project%os_ptcl2D%set_stkind(i, 1)
        call chunk_project%os_ptcl2D%set(i, 'indstk', i)
    enddo
    chunk_project%os_ptcl3D = chunk_project%os_ptcl2D
    call chunk_project%add_cavgs2os_out(cavg_path, SMPD, imgkind='cavg')
    do icls = 1, NCLASSES
        call chunk_project%os_cls2D%set_class(icls, icls)
        call chunk_project%os_cls2D%set_state(icls, 1)
        call chunk_project%os_cls2D%set(icls, 'pop', NPARTICLES_PER_CLASS)
        call chunk_project%os_cls2D%set(icls, 'res', 10.0)
        call chunk_project%os_cls2D%set(icls, 'corr', 0.9)
    enddo
    chunk_project%os_cls3D = chunk_project%os_cls2D
    call chunk_project%update_projinfo(chunk_projfile)

    ! Use the production ptcl_sieve collector in coarse-only mode.  Disabling
    ! the learned model isolates the deterministic sieve hard-gate contract.
    call cline_sieve%set('prg',               'sieve_cavgs')
    call cline_sieve%set('projfile',          chunk_projfile)
    call cline_sieve%set('dir_target',        fixture_root)
    call cline_sieve%set('ncls',              NCLASSES)
    call cline_sieve%set('nptcls_per_cls',    NPARTICLES_PER_CLASS)
    call cline_sieve%set('nchunksperset',     1)
    call cline_sieve%set('nchunks',           1)
    call cline_sieve%set('nparts',            1)
    call cline_sieve%set('nthr',              1)
    call cline_sieve%set('nptcls_coarse',     NPARTICLES)
    call cline_sieve%set('ncls_coarse',       NCLASSES)
    call cline_sieve%set('box_coarse',        CAVG_BOX)
    call cline_sieve%set('nsample_coarse',    NPARTICLES)
    call cline_sieve%set('lpstart',           20.0)
    call cline_sieve%set('lpstop_coarse',     15.0)
    call cline_sieve%set('mskdiam',           80.0)
    call cline_sieve%set('single_pass',       'yes')
    call cline_sieve%set('use_model',         'no')
    call cline_sieve%set('qsys_name',         'local')
    call cline_sieve%set('walltime',          60)
    call cline_sieve%set('mkdir',             'no')
    call cline_sieve%printline(unit=6)
    call chunk_project%update_compenv(cline_sieve)
    call chunk_project%write(chunk_projfile)
    call chunk_project%kill()
    call simple_touch(filepath(chunk_dir, ABINITIO2D_FINISHED))
    call params_sieve%new(cline_sieve)
    call sieve%new(params_sieve, completed_path)
    call sieve%collect_and_reject()

    if( sieve%get_n_chunks_coarse() /= 1 ) THROW_HARD('TEST_SIEVE_CAVGS FAILED: coarse chunk was not imported')
    if( sieve%get_n_coarse_accepted_ptcls() /= NPARTICLES_PER_CLASS )&
        &THROW_HARD('TEST_SIEVE_CAVGS FAILED: accepted-particle count is incorrect')
    if( sieve%get_n_coarse_rejected_ptcls() /= NPARTICLES_PER_CLASS )&
        &THROW_HARD('TEST_SIEVE_CAVGS FAILED: rejected-particle count is incorrect')
    if( sieve%get_n_accepted_ptcls() /= NPARTICLES_PER_CLASS .or. &
        &sieve%get_n_rejected_ptcls() /= NPARTICLES_PER_CLASS .or. &
        &sieve%get_n_total_particles() /= NPARTICLES )&
        &THROW_HARD('TEST_SIEVE_CAVGS FAILED: final sieve counters are inconsistent')
    if( sieve%get_n_accepted_micrographs() /= 1 )&
        &THROW_HARD('TEST_SIEVE_CAVGS FAILED: accepted-micrograph count is incorrect')
    if( .not. sieve%get_finished() ) THROW_HARD('TEST_SIEVE_CAVGS FAILED: coarse-only sieve did not finish')
    if( .not. file_exists(filepath(chunk_dir, 'REJECTION_FINISHED')) )&
        &THROW_HARD('TEST_SIEVE_CAVGS FAILED: rejection sentinel was not created')
    if( .not. file_exists(filepath(chunk_dir, 'COMPLETE')) )&
        &THROW_HARD('TEST_SIEVE_CAVGS FAILED: completion sentinel was not created')

    call result%read(chunk_projfile)
    if( any(result%os_cls2D%get_all_asint('state') /= [1, 0]) )&
        &THROW_HARD('TEST_SIEVE_CAVGS FAILED: class-average selection is incorrect')
    if( result%os_ptcl2D%count_state_gt_zero() /= NPARTICLES_PER_CLASS .or. &
        &result%os_ptcl3D%count_state_gt_zero() /= NPARTICLES_PER_CLASS )&
        &THROW_HARD('TEST_SIEVE_CAVGS FAILED: class selection was not mapped to particles')
    ! The persisted orientation reader tokenizes character values at blanks,
    ! so the stable round-tripped portion of the full reason is its tier prefix.
    rejection_reason = result%os_cls2D%get_str(2, 'rejection_reason')
    write(*,'(a,a)') '>>> TEST_SIEVE_CAVGS REJECTION REASON: ', rejection_reason%to_char()
    if( .not. rejection_reason%has_substr('coarse_reject') )&
        &THROW_HARD('TEST_SIEVE_CAVGS FAILED: rejected class does not record coarse rejection')
    call result%kill()

    completed_projfile = filepath(completed_path, CHUNK_STEM//METADATA_EXT)
    if( .not. file_exists(completed_projfile) )&
        &THROW_HARD('TEST_SIEVE_CAVGS FAILED: completed project was not exported')
    selected_jpeg = filepath(chunk_dir, CHUNK_STEM//'_selected'//JPG_EXT)
    rejected_jpeg = filepath(chunk_dir, CHUNK_STEM//'_rejected'//JPG_EXT)
    reasons_jpeg  = filepath(chunk_dir, CHUNK_STEM//'_all_reasons'//JPG_EXT)
    reasons_key   = reasons_jpeg//'.key.txt'
    if( .not. file_exists(selected_jpeg) .or. .not. file_exists(rejected_jpeg) )&
        &THROW_HARD('TEST_SIEVE_CAVGS FAILED: selected/rejected previews were not created')
    if( .not. file_exists(reasons_jpeg) .or. .not. file_exists(reasons_key) )&
        &THROW_HARD('TEST_SIEVE_CAVGS FAILED: rejection-reason report was not created')

    has_latest = sieve%get_latest(latest_inds, latest_pops, latest_res, latest_jpeg, latest_stk, &
        &xtiles, ytiles, latest_selection)
    if( .not. has_latest ) THROW_HARD('TEST_SIEVE_CAVGS FAILED: latest class-average metadata is missing')
    if( size(latest_inds) /= NCLASSES .or. any(latest_inds /= [1, 2]) .or. &
        &any(latest_pops /= [NPARTICLES_PER_CLASS, NPARTICLES_PER_CLASS]) .or. &
        &any(latest_selection /= [1, 0]) )&
        &THROW_HARD('TEST_SIEVE_CAVGS FAILED: latest class-average metadata is incorrect')
    if( any(abs(latest_res - 10.0) > 0.01) .or. xtiles * ytiles < NCLASSES )&
        &THROW_HARD('TEST_SIEVE_CAVGS FAILED: latest preview geometry or resolution is incorrect')
    if( .not. file_exists(latest_jpeg) .or. .not. file_exists(latest_stk) )&
        &THROW_HARD('TEST_SIEVE_CAVGS FAILED: latest preview files are missing')
    call find_ldim_nptcls(latest_stk, ldim, nimages)
    if( nimages /= NCLASSES .or. any(ldim(1:2) /= [CAVG_BOX, CAVG_BOX]) )&
        &THROW_HARD('TEST_SIEVE_CAVGS FAILED: retained class-average stack is invalid')

    if( allocated(latest_inds)      ) deallocate(latest_inds)
    if( allocated(latest_pops)      ) deallocate(latest_pops)
    if( allocated(latest_res)       ) deallocate(latest_res)
    if( allocated(latest_selection) ) deallocate(latest_selection)
    call sieve%kill()
    call cline_sieve%kill()
    call simple_chdir(cwd_root, status)
    if( status /= 0 ) THROW_HARD('TEST_SIEVE_CAVGS FAILED: could not restore the original directory')
    write(*,'(a,a)') '>>> TEST_SIEVE_CAVGS: validated one selected and one rejected class in ', fixture_root%to_char()
    call flush(6)
    call simple_end('**** SIMPLE_TEST_SIEVE_CAVGS NORMAL STOP ****')
end subroutine exec_test_sieve_cavgs

end module simple_commanders_test_stream
