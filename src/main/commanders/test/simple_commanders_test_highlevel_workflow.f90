!@descr: simulated end-to-end workflow test of simple_commanders_test_highlevel
submodule (simple_commanders_test_highlevel) simple_commanders_test_highlevel_workflow
implicit none
#include "simple_local_flags.inc"

contains

module subroutine exec_test_simulated_workflow( self, cline )
    use simple_atoms,              only: atoms
    use simple_test_truth_metrics, only: validate_reconstructed_volume
    use simple_molecule_data,      only: molecule_data, betagal_1jyx, sars_cov2_spkgp_6vxx
    use simple_refine3D_fnames,    only: refine3D_state_vol_fname
    use simple_string_utils,       only: lowercase
    use simple_ui,                 only: make_ui
    class(commander_test_simulated_workflow), intent(inout) :: self
    class(cmdline),                           intent(inout) :: cline
    character(len=*), parameter :: PROJNAME       = 'simulated_workflow'
    character(len=*), parameter :: PROJFILE       = PROJNAME//'.simple'
    character(len=*), parameter :: MOVIE_FILE     = 'simulate_movie.mrc'
    character(len=*), parameter :: SUBSET_FILE    = 'random_reprojections.mrcs'
    character(len=*), parameter :: OPTIMAL_FILE   = 'optimal_movie_average.mrc'
    character(len=*), parameter :: PARAMS_FILE    = 'simulate_movie_params.txt'
    character(len=*), parameter :: PICKREFS_FILE  = 'pickrefs.mrc'
    character(len=*), parameter :: FILETAB_FILE   = 'simulated_movies.txt'
    character(len=*), parameter :: VOL_DIR        = '0_pdb2mrc'
    character(len=*), parameter :: REPROJ_DIR     = '1_reproject'
    character(len=*), parameter :: IMPORT_DIR     = '1_import_movies'
    character(len=*), parameter :: MOTION_DIR     = '2_motion_correct'
    character(len=*), parameter :: CTF_DIR        = '3_ctf_estimate'
    character(len=*), parameter :: PICK_DIR       = '4_pick'
    character(len=*), parameter :: EXTRACT_DIR    = '5_extract'
    character(len=*), parameter :: ABINIT2D_DIR   = '6_solve2D'
    real,             parameter :: SMPD           = 1.3
    real,             parameter :: MSKDIAM        = 180.0
    real,             parameter :: CS             = 2.7
    real,             parameter :: KV             = 300.0
    real,             parameter :: FRACA          = 0.1
    integer,          parameter :: NPROJS         = 100
    integer,          parameter :: NPER_MOVIE     = 10
    integer,          parameter :: MOVIE_DIM      = 1536
    integer,          parameter :: NFRAMES        = 16
    integer,          parameter :: NMOVIES        = 10
    integer,          parameter :: EXTRACT_BOX    = 192
    integer,          parameter :: NTHR           = 4
    integer,          parameter :: TEST_SEED      = 20260923
    real,             parameter :: MIN_VOL_CORR   = 0.80
    real,             parameter :: MAX_FSC0143    = 40.0
    real,             parameter :: DOCK_HP        = 100.0
    real,             parameter :: DOCK_LP        = 20.0
    type(cmdline)                       :: cline_projection, cline_sim_mov, cline_new_project
    type(cmdline)                       :: cline_import_movies, cline_mot_corr, cline_ctf_est
    type(cmdline)                       :: cline_pick, cline_extract, cline_solve2D, cline_solve3D
    type(commander_new_project)         :: xnew_project
    type(commander_reproject)           :: xreproject
    type(commander_simulate_movie)      :: xsimov
    type(commander_motion_correct)      :: xmotcorr
    type(commander_ctf_estimate)        :: xctf_estimate
    type(commander_import_movies)       :: ximport_movies
    type(commander_pick)                :: xpick
    type(commander_extract)             :: xextract
    type(commander_solve2D)             :: xsolve2D
    type(commander_solve3D)             :: xsolve3D
    type(molecule_data)                 :: mol
    type(atoms)                         :: molecule
    type(image)                         :: projection
    type(sp_project)                    :: spproj
    type(string)                        :: cwd_root, workflow_root, project_path, reproj_path, subset_path, filetab_path
    type(string)                        :: system_name, workflow_picker, pgrp, test_workdir, vol_file, reproj_file
    type(string)                        :: movie_fname, subset_fname, optimal_fname, params_fname
    type(string)                        :: truth_volume, solve3D_dir, final_volume
    type(string)                        :: movie_files(NMOVIES)
    character(len=XLONGSTRLEN)          :: workflow_root_path
    integer                             :: i, j, proj_inds(NPER_MOVIE), projection_order(NPROJS), ldim(3), nprojs_stk
    integer                             :: reproj_box, npickrefs, nptcls, ncls, status
    real                                :: pickref_smpd, pickref_width
    real                                :: volume_corr, volume_fsc0143
    real                                :: dock_corr_direct, dock_corr_mirrored, dock_corr_selected
    integer                             :: rnd_defocus
    logical                             :: volume_ok

    ! The test executable initializes only the test UI.  Directly invoked SIMPLE
    ! commanders need the regular UI metadata to implement mkdir=yes correctly.
    call make_ui
    if( .not. cline%defined('suite') ) THROW_HARD('The suite keyword is required; use suite=6vxx or suite=1jxy')
    system_name = cline%get_carg('suite')
    system_name = lowercase(system_name%to_char())
    if( system_name == 'list' )then
        write(logfhandle,'(a)') 'Available suites for simulated_workflow:'
        write(logfhandle,'(a)') '  6vxx'
        write(logfhandle,'(a)') '  1jxy'
        return
    endif
    workflow_picker = 'segdiam'
    if( cline%defined('picker') ) workflow_picker = cline%get_carg('picker')
    workflow_picker = lowercase(workflow_picker%to_char())
    select case(workflow_picker%to_char())
        case('segdiam','new')
        case default
            THROW_HARD('Simulated workflow picker must be segdiam or new')
    end select
    select case(system_name%to_char())
        case('6vxx')
            vol_file    = '6VXX.mrc'
            reproj_file = 'reprojs_6VXX.mrcs'
            pgrp        = 'c3'
        case('1jxy')
            vol_file    = '1JXY.mrc'
            reproj_file = 'reprojs_1JXY.mrcs'
            pgrp        = 'd2'
        case default
            THROW_HARD('no sub-suite '//system_name%to_char()//' in simulated_workflow; use suite=list')
    end select
    call set_fixed_seed(TEST_SEED, propagate=.true.)
    write(logfhandle,'(a,i0)') '>>> Deterministic workflow seed: ', TEST_SEED
    test_workdir = 'test_simulated_workflow_'//system_name%to_char()
    call simple_getcwd(cwd_root)
    if( file_exists(test_workdir%to_char()) )then
        call simple_rmdir(test_workdir%to_char(), status)
        if( status /= 0 ) THROW_HARD('Could not reset '//test_workdir%to_char())
    endif
    call simple_mkdir(test_workdir%to_char())
    call simple_chdir(test_workdir%to_char(), status)
    if( status /= 0 ) THROW_HARD('Could not enter '//test_workdir%to_char())
    call simple_getcwd(workflow_root)
    workflow_root_path = workflow_root%to_char()

    ! Both coordinate sets are embedded in SIMPLE, so this test has no network dependency.
    write(logfhandle,'(a,a)') '>>> Step 1: create a volume from ', system_name%to_char()
    call simple_mkdir(VOL_DIR)
    call simple_chdir(VOL_DIR, status)
    if( status /= 0 ) THROW_HARD('Could not enter the volume-generation directory')
    select case(system_name%to_char())
        case('6vxx')
            mol = sars_cov2_spkgp_6vxx()
        case('1jxy')
            ! SIMPLE's embedded provider uses the underlying 1JYX PDB identifier.
            mol = betagal_1jyx()
    end select
    call molecule%pdb2mrc(volfile=vol_file, smpd=SMPD, mol=mol, center_pdb=.true., &
        &vol_dim=[EXTRACT_BOX, EXTRACT_BOX, EXTRACT_BOX])
    call molecule%kill()
    truth_volume = simple_abspath(vol_file)
    call simple_chdir(workflow_root, status)
    if( status /= 0 ) THROW_HARD('Could not leave the volume-generation directory')

    write(logfhandle,'(a)') '>>> Step 2: generate well-spaced spiral reprojections'
    call cline_projection%set('prg',                       'reproject')
    call cline_projection%set('mkdir',                           'yes')
    call cline_projection%set('vol1', VOL_DIR//'/'//vol_file%to_char())
    call cline_projection%set('outstk',                    reproj_file)
    call cline_projection%set('smpd',                             SMPD)
    call cline_projection%set('pgrp',                             pgrp)
    call cline_projection%set('mskdiam',                       MSKDIAM)
    call cline_projection%set('nspace',                         NPROJS)
    call cline_projection%set('nthr',                             NTHR)
    call xreproject%execute(cline_projection)
    call cline_projection%kill()
    call return_to_stage_root('reproject')
    reproj_path = simple_abspath(string(REPROJ_DIR//'/'//reproj_file%to_char()))
    call find_ldim_nptcls(reproj_path, ldim, nprojs_stk)
    if( nprojs_stk /= NPROJS ) THROW_HARD('Unexpected number of generated reprojections')
    reproj_box = ldim(1)
    ldim(3) = 1
    call projection%new(ldim, SMPD)

    write(logfhandle,'(a,i0,a,i0,a,i0,a)') '>>> Step 3: generate ', NMOVIES, &
        &' simulated movies with ', NFRAMES, ' frames and ', NPER_MOVIE, ' shuffled spiral projections each'
    if( NPROJS /= NMOVIES * NPER_MOVIE ) THROW_HARD('Projection count must divide equally among simulated movies')
    do i = 1, NPROJS
        projection_order(i) = i
    enddo
    call seed_rnd
    call shuffle(projection_order)
    do i = 1,NMOVIES
        proj_inds = projection_order((i - 1) * NPER_MOVIE + 1:i * NPER_MOVIE)
        if( file_exists(SUBSET_FILE) ) call del_file(SUBSET_FILE)
        do j = 1,NPER_MOVIE
            call projection%read(reproj_path, proj_inds(j))
            call projection%write(string(SUBSET_FILE), j)
        enddo
        subset_path = simple_abspath(string(SUBSET_FILE))
        call cline_sim_mov%set('prg',           'simulate_movie')
        if( i == 1 )then
            call cline_sim_mov%set('mkdir',                'yes')
            call cline_sim_mov%set('dir_exec', 'simulate_movies')
        else
            call cline_sim_mov%set('mkdir',                 'no')
        endif
        call cline_sim_mov%set('stk',                subset_path)
        call cline_sim_mov%set('xdim',                 MOVIE_DIM)
        call cline_sim_mov%set('ydim',                 MOVIE_DIM)
        call cline_sim_mov%set('nframes',                NFRAMES)
        call cline_sim_mov%set('smpd',                      SMPD)
        call cline_sim_mov%set('snr',                        0.2)
        call cline_sim_mov%set('kv',                          KV)
        call cline_sim_mov%set('cs',                          CS)
        call cline_sim_mov%set('fraca',                    FRACA)
        call seed_rnd
        rnd_defocus = irnd_uni(3)
        call cline_sim_mov%set('defocus',            rnd_defocus)
        call cline_sim_mov%set('trs',                        2.0)
        call cline_sim_mov%set('nthr',                      NTHR)
        call xsimov%execute(cline_sim_mov)
        if( .not. file_exists(MOVIE_FILE) ) THROW_HARD('Simulated movie was not generated')
        movie_fname = string('simulate_movie_')//int2str_pad(i,3)//MRC_EXT
        subset_fname = string('random_reprojections_')//int2str_pad(i,3)//'.mrcs'
        optimal_fname = string('optimal_movie_average_')//int2str_pad(i,3)//MRC_EXT
        params_fname = string('simulate_movie_params_')//int2str_pad(i,3)//TXT_EXT
        call simple_rename(MOVIE_FILE, movie_fname)
        call simple_rename(subset_path, subset_fname)
        call simple_rename(OPTIMAL_FILE, optimal_fname)
        call simple_rename(PARAMS_FILE, params_fname)
        movie_files(i) = simple_abspath(movie_fname)
        call cline_sim_mov%kill()
    enddo
    call projection%kill()
    call return_to_stage_root('simulate_movie')
    call write_filetable(string(FILETAB_FILE), movie_files)
    filetab_path = simple_abspath(string(FILETAB_FILE))

    write(logfhandle,'(a)') '>>> Step 4: create a project and import the movies'
    call cline_new_project%set('projname',              PROJNAME)
    call cline_new_project%set('qsys_name',              'local')
    call xnew_project%execute(cline_new_project)
    call cline_new_project%kill()
    call simple_getcwd(workflow_root)
    workflow_root_path = workflow_root%to_char()
    project_path       = simple_abspath(string(PROJFILE))
    call cline_import_movies%set('prg',          'import_movies')
    call cline_import_movies%set('mkdir',                  'yes')
    call cline_import_movies%set('projfile',        project_path)
    call cline_import_movies%set('filetab',         filetab_path)
    call cline_import_movies%set('cs',                        CS)
    call cline_import_movies%set('fraca',                  FRACA)
    call cline_import_movies%set('kv',                        KV)
    call cline_import_movies%set('smpd',                    SMPD)
    call cline_import_movies%set('ctf',                    'yes')
    call ximport_movies%execute(cline_import_movies)
    call cline_import_movies%kill()
    call update_project_path
    call return_to_stage_root('import_movies')

    write(logfhandle,'(a)') '>>> Step 5: motion correction'
    call cline_mot_corr%set('prg',              'motion_correct')
    call cline_mot_corr%set('projfile',             project_path)
    call cline_mot_corr%set('mkdir',                       'yes')
    call cline_mot_corr%set('nparts',                          1)
    call cline_mot_corr%set('nthr',                         NTHR)
    call xmotcorr%execute(cline_mot_corr)
    call cline_mot_corr%kill()
    call update_project_path
    call return_to_stage_root('motion_correct')

    write(logfhandle,'(a)') '>>> Step 6: CTF estimation'
    call cline_ctf_est%set('prg',                 'ctf_estimate')
    call cline_ctf_est%set('projfile',              project_path)
    call cline_ctf_est%set('mkdir',                        'yes')
    call cline_ctf_est%set('nparts',                           1)
    call cline_ctf_est%set('nthr',                          NTHR)
    call xctf_estimate%execute(cline_ctf_est)
    call cline_ctf_est%kill()
    call update_project_path
    call return_to_stage_root('ctf_estimate')

    write(logfhandle,'(a,a)') '>>> Step 7: particle picking with ', workflow_picker%to_char()
    call cline_pick%set('prg',                            'pick')
    call cline_pick%set('projfile',                 project_path)
    call cline_pick%set('mkdir',                           'yes')
    call cline_pick%set('pcontrast',                     'black')
    select case(workflow_picker%to_char())
        case('segdiam')
            call cline_pick%set('picker',              'segdiam')
            call cline_pick%set('moldiam_max',           MSKDIAM)
        case('new')
            call cline_pick%set('picker',                  'new')
            call cline_pick%set('pickrefs',          reproj_path)
            call cline_pick%set('moldiam',               MSKDIAM)
            call cline_pick%set('pick_roi',                 'no')
        case default
            THROW_HARD('Unsupported simulated-workflow picker')
    end select
    call cline_pick%set('nparts',                              1)
    call cline_pick%set('nthr',                             NTHR)
    call xpick%execute(cline_pick)
    call cline_pick%kill()
    if( workflow_picker%to_char() == 'new' )then
        if( .not. file_exists(PICKREFS_FILE) ) THROW_HARD('Picking references were not generated')
        call find_ldim_nptcls(string(PICKREFS_FILE), ldim, npickrefs)
        pickref_smpd = find_img_smpd(string(PICKREFS_FILE))
        pickref_width = real(ldim(1)) * pickref_smpd
        if( ldim(1) /= ldim(2) .or. ldim(1) > reproj_box )then
            THROW_HARD('Generated picking-reference dimensions are inconsistent with the reprojections')
        endif
        if( abs(pickref_smpd - SMPD) > 0.01 )then
            THROW_HARD('Generated picking-reference sampling distance is inconsistent with the micrographs')
        endif
        if( pickref_width < 0.75 * MSKDIAM .or. pickref_width > 1.25 * MSKDIAM )then
            THROW_HARD('Generated picking-reference physical size is inconsistent with the particle diameter')
        endif
        if( EXTRACT_BOX < ldim(1) )then
            THROW_HARD('Extraction box is smaller than the generated picking references')
        endif
        write(logfhandle,'(a,i0,a,f6.1,a)') '>>> VALIDATED PICKING REFERENCES: ', ldim(1), &
            &' pixels, ', pickref_width, ' A across'
    endif
    call update_project_path
    call return_to_stage_root('pick')

    write(logfhandle,'(a)') '>>> Step 8: particle extraction'
    call cline_extract%set('prg',                      'extract')
    call cline_extract%set('projfile',              project_path)
    call cline_extract%set('mkdir',                        'yes')
    call cline_extract%set('box',                    EXTRACT_BOX)
    call cline_extract%set('nparts',                           1)
    call cline_extract%set('nthr',                          NTHR)
    call xextract%execute(cline_extract)
    call cline_extract%kill()
    call update_project_path
    call return_to_stage_root('extract')

    call spproj%read(project_path)
    nptcls = spproj%get_nptcls()
    call spproj%kill()
    if( nptcls < 4 ) THROW_HARD('Too few particles were extracted for initial model tests')
    ncls = min(4, max(2, nptcls / 5))

    write(logfhandle,'(a)') '>>> Step 9: solve2D'
    call cline_solve2D%set('prg',                   'solve2D')
    call cline_solve2D%set('projfile',              project_path)
    call cline_solve2D%set('mkdir',                        'yes')
    call cline_solve2D%set('mskdiam',                    MSKDIAM)
    call cline_solve2D%set('ncls',                          ncls)
    call cline_solve2D%set('nthr',                          NTHR)
    call xsolve2D%execute(cline_solve2D)
    call cline_solve2D%kill()
    call update_project_path
    call return_to_stage_root('solve2D')

    call spproj%read(project_path)
    if( spproj%os_cls2D%get_noris() < 1 ) THROW_HARD('solve2D produced no classes')
    call spproj%kill()

    write(logfhandle,'(a)') '>>> Step 10: solve3D'
    call cline_solve3D%set('prg',                   'solve3D')
    call cline_solve3D%set('projfile',              project_path)
    call cline_solve3D%set('mkdir',                        'yes')
    call cline_solve3D%set('pgrp',                          pgrp)
    if( system_name == '1jxy' ) call cline_solve3D%set('pgrp_start', pgrp)
    call cline_solve3D%set('mskdiam',                    MSKDIAM)
    call cline_solve3D%set('nthr',                          NTHR)
    call xsolve3D%execute(cline_solve3D)
    call simple_getcwd(solve3D_dir)
    call cline_solve3D%kill()
    final_volume = filepath(solve3D_dir, refine3D_state_vol_fname(1))
    call validate_reconstructed_volume(truth_volume, final_volume, SMPD, EXTRACT_BOX, 0.01, MSKDIAM, &
        &DOCK_HP, DOCK_LP, MIN_VOL_CORR, MAX_FSC0143, volume_corr, volume_fsc0143, &
        &dock_corr_direct, dock_corr_mirrored, dock_corr_selected, volume_ok)

    call simple_chdir(cwd_root, status)
    if( status /= 0 ) THROW_HARD('Could not restore the original working directory')
    if( .not. volume_ok ) THROW_HARD('TEST_SIMULATED_WORKFLOW FAILED: final-volume validation failed')
    write(logfhandle,'(a,a,a,f7.4,a,f7.4,a,f7.4,a,f7.4,a,f7.2,a)') 'PASS: simulated_workflow ', &
        &system_name%to_char(), ' docking correlation direct=', dock_corr_direct, ', mirrored=', dock_corr_mirrored, &
        &', selected=', dock_corr_selected, ', whole-volume correlation=', volume_corr, &
        &', FSC=0.143 at ', volume_fsc0143, ' A'
    call simple_end('**** SIMPLE_TEST_SIMULATED_WORKFLOW NORMAL STOP ****')

  contains

    subroutine update_project_path
        if( file_exists(PROJFILE) ) project_path = simple_abspath(string(PROJFILE))
    end subroutine update_project_path

    subroutine return_to_stage_root( stage )
        character(len=*), intent(in) :: stage
        call simple_chdir(trim(workflow_root_path), status)
        if( status /= 0 ) THROW_HARD('Could not leave simulated workflow stage: '//stage)
    end subroutine return_to_stage_root

end subroutine exec_test_simulated_workflow

end submodule simple_commanders_test_highlevel_workflow
