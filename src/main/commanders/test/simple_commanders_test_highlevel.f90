!@descr: for all highlevel tests
module simple_commanders_test_highlevel
use simple_commanders_api
use simple_stream_api
use simple_commanders_project_core, only: commander_new_project, commander_selection
use simple_commanders_project_mov,  only: commander_import_movies
use simple_commanders_reproject,    only: commander_reproject
use simple_commanders_pick,         only: commander_pick, commander_extract
use simple_commanders_sim,          only: commander_simulate_noise, commander_simulate_particles, commander_simulate_movie
use simple_commanders_preprocess,   only: commander_ctf_estimate, commander_motion_correct, commander_preprocess
use simple_commanders_abinitio2D,   only: commander_abinitio2D
use simple_commanders_abinitio,     only: commander_abinitio3D
use simple_commanders_stkops,       only: commander_stack
use simple_micproc,                 only: sample_filetab
use simple_commanders_validate,     only: commander_mini_stream
use simple_qsys_env,                only: qsys_env
use simple_projfile_utils,          only: merge_chunk_projfiles
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_mini_stream
  contains
    procedure :: execute      => exec_test_mini_stream
end type commander_test_mini_stream

type, extends(commander_base) :: commander_test_simulate_particles
  contains
    procedure :: execute      => exec_test_simulate_particles
end type commander_test_simulate_particles

type, extends(commander_base) :: commander_test_reproject
  contains
    procedure :: execute      => exec_test_reproject
end type commander_test_reproject

type, extends(commander_base) :: commander_test_simulated_workflow
  contains
    procedure :: execute      => exec_test_simulated_workflow
end type commander_test_simulated_workflow

type, extends(commander_base) :: commander_test_subproject_distr
  contains
    procedure :: execute      => exec_test_subproject_distr
end type commander_test_subproject_distr

type, extends(commander_base) :: commander_test_ptcls_ppca_subproject_distr
    contains
        procedure :: execute      => exec_test_ptcls_ppca_subproject_distr
end type commander_test_ptcls_ppca_subproject_distr

type, extends(commander_base) :: commander_test_pcg_recon
  contains
    procedure :: execute      => exec_test_pcg_recon
end type commander_test_pcg_recon

contains

subroutine exec_test_mini_stream( self, cline )
    class(commander_test_mini_stream),  intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    real,         parameter       :: CTFRES_THRES = 8.0, ICE_THRES = 1.0, OVERSHOOT = 1.2
    type(string), allocatable     :: dataset_cmds(:)
    type(string), allocatable     :: micstab(:), filetab(:), movfnames(:)
    type(string)                  :: output_dir, imgkind
    integer,      allocatable     :: orimap(:)
    type(cmdline)                 :: cline_dataset, cline_new_project, cline_import_movies, cline_preprocess
    type(cmdline)                 :: cline_select, cline_mini_stream
    type(parameters)              :: params
    type(commander_new_project)   :: xnew_project
    type(commander_preprocess)    :: xpreprocess
    type(commander_import_movies) :: ximport_movies
    type(commander_mini_stream)   :: xmini_stream
    type(commander_selection)     :: xsel
    type(sp_project)              :: spproj
    type(stream_watcher)          :: movie_buff
    integer                       :: i, ndata_sets, n_nonzero, nmovf
    type(string)                  :: abspath, projfile
    character(len=*), parameter   :: filetab_file='filetab.txt'
    ! Parsing
    if( command_argument_count() < 1 )then
        write(logfhandle,'(a)') 'ERROR! Usage: simple_test_mini_stream fname=filetab.txt'
        call exit(-1)
    else 
        call cline%parse_oldschool
    endif
    call cline%checkvar('fname',        1)
    call cline%check()
    call params%new(cline)
    call read_filetable(params%fname, dataset_cmds)
    ndata_sets=size(dataset_cmds)
    ! projname=system_name smpd=1.3 cs=2.7 kv=300 fraca=0.1 total_dose=53 dir_movies=/usr/local/data/movies gainref=gainref.mrc nparts=4 nthr=16 moldiam_max=200 nram=100
    call simple_getcwd(abspath)
    output_dir=abspath
    do i = 1, ndata_sets
        call cline_dataset%read(dataset_cmds(i)%to_char())
        call cline_dataset%checkvar('projname',        1)
        call cline_dataset%checkvar('smpd',            2)
        call cline_dataset%checkvar('cs',              3)
        call cline_dataset%checkvar('kv',              4)
        call cline_dataset%checkvar('fraca',           5)
        call cline_dataset%checkvar('total_dose',      6)
        call cline_dataset%checkvar('dir_movies',      7)
        call cline_dataset%checkvar('gainref',         8)
        call cline_dataset%checkvar('nparts',          9)
        call cline_dataset%checkvar('nthr',           10)
        call cline_dataset%checkvar('moldiam_max',    11)
        call cline_dataset%checkvar('nran',           12)
        call cline_dataset%check()
        call params%new(cline_dataset)
        call cline_dataset%kill()
        ! project creation
        call cline_new_project%set('projname',  params%projname)
        call xnew_project%execute(cline_new_project)
        call cline_new_project%kill()
        projfile = params%projname//'.simple'
        ! create filetab with a subset of overshoot randomly selected movies
        movie_buff = stream_watcher(1,params%dir_movies)
        call movie_buff%watch(nmovf, movfnames)
        call movie_buff%kill
        filetab = sample_filetab(movfnames, ceiling(real(params%nran)*OVERSHOOT))
        call write_filetable(string(filetab_file), filetab)
        ! movie import
        call cline_import_movies%set('prg',                'import_movies')
        call cline_import_movies%set('mkdir',                        'yes')
        call cline_import_movies%set('cs',                       params%cs)
        call cline_import_movies%set('fraca',                 params%fraca)
        call cline_import_movies%set('kv',                       params%kv)
        call cline_import_movies%set('smpd',                   params%smpd)
        call cline_import_movies%set('filetab',               filetab_file)
        call cline_import_movies%set('ctf',                          'yes')
        call ximport_movies%execute(cline_import_movies)
        call cline_import_movies%kill()
        ! check either movies or micrographs
        call spproj%read(projfile)
        imgkind = spproj%get_mic_kind(1)
        if( imgkind.eq.'intg' )then
            ! nothing to do
        else
            ! preprocess
            call simple_chdir(output_dir//'/'//params%projname%to_char())
            call cline_preprocess%set('prg',                  'preprocess')
            call cline_preprocess%set('mkdir',                       'yes')
            call cline_preprocess%set('gainref',            params%gainref)
            call cline_preprocess%set('total_dose',      params%total_dose)
            call cline_preprocess%set('dfmin',               DFMIN_DEFAULT)
            call cline_preprocess%set('dfmax',               DFMAX_DEFAULT)
            call cline_preprocess%set('hp',                            30.)
            call cline_preprocess%set('lp',                             2.)
            call cline_preprocess%set('mcpatch',                      'no')
            call cline_preprocess%set('nparts',              params%nparts)
            call cline_preprocess%set('nthr',                  params%nthr)
            call cline_preprocess%check()
            call xpreprocess%execute(cline_preprocess)
            call cline_preprocess%kill()
        endif
        ! reject based on CTF resolution and ice score
        call simple_chdir(output_dir//'/'//params%projname%to_char())
        call cline_select%delete('nran')
        call cline_select%set('prg',                           'selection')
        call cline_select%set('mkdir',                               'yes')
        call cline_select%set('oritype',                             'mic')
        call cline_select%set('ctfresthreshold',              CTFRES_THRES)
        call cline_select%set('icefracthreshold',                ICE_THRES)
        call xsel%execute(cline_select)
        call cline_select%kill()
        ! state=0/1 should now be in project file on disk
        ! re-run for random selection
        call spproj%read(projfile)
        n_nonzero = spproj%get_n_insegment_state('mic', 1)
        if( n_nonzero > params%nran )then
            ! make random selection
            call simple_chdir(output_dir//'/'//params%projname%to_char())
            call cline_select%delete('ctfresthreshold')
            call cline_select%delete('icefracthreshold')
            call cline_select%set('prg',                       'selection')
            call cline_select%set('mkdir',                           'yes') 
            call cline_select%set('oritype',                         'mic')
            call cline_select%set('nran',                      params%nran)
            call xsel%execute(cline_select)
            call cline_select%kill()
        endif
        call spproj%read(projfile)
        call spproj%get_mics_table(micstab, orimap)
        call simple_chdir(output_dir//'/'//params%projname%to_char())
        call write_filetable(string('intgs.txt'),micstab)
        ! mini stream 
        call cline_mini_stream%set('prg',                    'mini_stream')
        call cline_mini_stream%set('mkdir',                          'yes')
        call cline_mini_stream%set('filetab',                  'intgs.txt')
        call cline_mini_stream%set('smpd',                     params%smpd)
        call cline_mini_stream%set('fraca',                   params%fraca)
        call cline_mini_stream%set('kv',                         params%kv)
        call cline_mini_stream%set('cs',                         params%cs)
        call cline_mini_stream%set('moldiam_max',       params%moldiam_max)
        call cline_mini_stream%set('nparts',                 params%nparts)
        call cline_mini_stream%set('nthr',                     params%nthr)
        call xmini_stream%execute(cline_mini_stream)
        call cline_dataset%kill()
        call cline_new_project%kill()
        call cline_import_movies%kill()
        call cline_preprocess%kill()
        call cline_mini_stream%kill()
        call cline_dataset%kill()
        call simple_chdir(output_dir)
    enddo
    call simple_end('**** SIMPLE_TEST_MINI_STREAM_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_mini_stream

subroutine exec_test_simulate_particles( self, cline )
    use simple_atoms,         only: atoms
    use simple_molecule_data, only: molecule_data, sars_cov2_spkgp_6vxx
    use simple_imghead,       only: find_ldim_nptcls, find_img_smpd
    class(commander_test_simulate_particles), intent(inout) :: self
    class(cmdline),                           intent(inout) :: cline
    type(cmdline)                       :: cline_sim
    type(parameters)                    :: params
    type(commander_simulate_particles)  :: xsim_ptcls
    type(atoms)                         :: molecule
    real, parameter                     :: smpd = 1.3
    type(molecule_data)                 :: mol
    integer, parameter                  :: NPTCLS_SIM = 200
    type(string)                        :: outstk, outfile, vol_file
    integer                             :: ldim(3), nptcls_stk, nlines_ori
    real                                :: smpd_stk
    logical                             :: all_ok
    mol = sars_cov2_spkgp_6vxx()
    call molecule%pdb2mrc(smpd=smpd, mol=mol, center_pdb=.true.)
    call params%new(cline)
    all_ok = .true.
    ! ---- simulate particles ----
    write(logfhandle,'(a)') '>>> TEST_SIMULATE_PARTICLES:'
    call cline_sim%set('prg',      'simulate_particles')
    call cline_sim%set('vol1',           'molecule.mrc')
    call cline_sim%set('smpd',                     smpd)
    call cline_sim%set('mskdiam',                   180)
    call cline_sim%set('nthr',                       16)
    call cline_sim%set('nptcls',             NPTCLS_SIM)
    call cline_sim%set('pgrp',                     'c1')
    call cline_sim%set('snr',                      0.01)
    call cline_sim%set('ctf',                     'yes')
    call cline_sim%set('sherr',                     0.0)
    call cline_sim%set('even',                     'on')
    call xsim_ptcls%execute(cline_sim)
    ! ---- define expected output file names ----
    vol_file = 'molecule.mrc'
    outstk   = 'simulated_particles.mrc'
    outfile  = 'simulated_oris'//trim(TXT_EXT)
    ! ---- check volume was generated ----
    write(logfhandle,'(a)') '>>> CHECK: volume file exists'
    if( .not. file_exists(vol_file) )then
        write(logfhandle,'(a)') '    FAIL: '//vol_file%to_char()//' not found'
        THROW_HARD('TEST_SIMULATE_PARTICLES FAILED: volume not generated')
    else
        call find_ldim_nptcls(vol_file, ldim, nptcls_stk)
        smpd_stk = find_img_smpd(vol_file)
        write(logfhandle,'(a,i4,a,i4,a,i4,a,f6.2)') '    PASS: volume dims = [', &
            ldim(1),',',ldim(2),',',ldim(3),' ], smpd = ', smpd_stk
        if( ldim(1) /= ldim(2) .or. ldim(1) < 1 )then
            write(logfhandle,'(a)') '    FAIL: volume has invalid dimensions'
            all_ok = .false.
        endif
        if( abs(smpd_stk - smpd) > 0.01 )then
            write(logfhandle,'(a,f6.2,a,f6.2)') '    FAIL: smpd mismatch, expected ', smpd, ' got ', smpd_stk
            all_ok = .false.
        endif
    endif
    ! ---- validate output stack ----
    write(logfhandle,'(a)') '>>> CHECK: output particle stack'
    if( .not. file_exists(outstk) )then
        write(logfhandle,'(a)') '    FAIL: '//outstk%to_char()//' not found'
        all_ok = .false.
    else
        call find_ldim_nptcls(outstk, ldim, nptcls_stk)
        smpd_stk = find_img_smpd(outstk)
        write(logfhandle,'(a,i6)')  '    particles in stack: ', nptcls_stk
        write(logfhandle,'(a,i4,a,i4)') '    box size:           ', ldim(1), ' x ', ldim(2)
        write(logfhandle,'(a,f6.2)')    '    smpd:               ', smpd_stk
        if( nptcls_stk /= NPTCLS_SIM )then
            write(logfhandle,'(a,i6,a,i6)') '    FAIL: expected ', NPTCLS_SIM, ' particles, got ', nptcls_stk
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: particle count matches'
        endif
        if( ldim(1) /= ldim(2) .or. ldim(1) < 1 )then
            write(logfhandle,'(a)') '    FAIL: invalid box dimensions'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: box dimensions valid'
        endif
        if( abs(smpd_stk - smpd) > 0.01 )then
            write(logfhandle,'(a,f6.2,a,f6.2)') '    FAIL: smpd mismatch, expected ', smpd, ' got ', smpd_stk
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: sampling distance matches'
        endif
    endif
    ! ---- validate orientations file ----
    write(logfhandle,'(a)') '>>> CHECK: orientations file'
    if( .not. file_exists(outfile) )then
        write(logfhandle,'(a)') '    FAIL: '//outfile%to_char()//' not found'
        all_ok = .false.
    else
        nlines_ori = nlines(outfile)
        write(logfhandle,'(a,i6)') '    orientation records: ', nlines_ori
        if( nlines_ori /= NPTCLS_SIM )then
            write(logfhandle,'(a,i6,a,i6)') '    FAIL: expected ', NPTCLS_SIM, ' records, got ', nlines_ori
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: orientation count matches'
        endif
    endif
    ! ---- final verdict ----
    if( all_ok )then
        call simple_end('**** SIMPLE_TEST_SIMULATE_PARTICLES NORMAL STOP ****')
    else
        THROW_HARD('TEST_SIMULATE_PARTICLES FAILED')
    endif
end subroutine exec_test_simulate_particles

subroutine exec_test_reproject( self, cline )
    use simple_atoms,         only: atoms
    use simple_molecule_data, only: molecule_data, sars_cov2_spkgp_6vxx
    use simple_imghead,       only: find_ldim_nptcls, find_img_smpd
    class(commander_test_reproject), intent(inout) :: self
    class(cmdline),                  intent(inout) :: cline
    integer, parameter              :: NSPACE = 100
    real,    parameter              :: SMPD   = 1.3
    type(cmdline)                   :: cline_reproj
    type(parameters)                :: params
    type(commander_reproject)       :: xreproject
    type(atoms)                     :: molecule
    type(molecule_data)             :: mol
    type(string)                    :: vol_file, outstk, outori
    integer                         :: ldim(3), nptcls_stk, nlines_ori
    real                            :: smpd_stk
    logical                         :: all_ok
        ! ---- define expected output file names ----
    vol_file = '6VXX.mrc'
    outstk   = 'reprojs.mrcs'
    outori   = 'reproject_oris'//trim(TXT_EXT)
    ! ---- generate 6VXX volume from built-in molecule data ----
    write(logfhandle,'(a)') '>>> TEST_REPROJECT: generating 6VXX.mrc volume'
    mol = sars_cov2_spkgp_6vxx()
    call molecule%pdb2mrc(smpd=SMPD, volfile=vol_file, mol=mol)
    call params%new(cline)
    all_ok = .true.
    ! ---- check volume was generated ----
    write(logfhandle,'(a)') '>>> CHECK: volume file exists'
    if( .not. file_exists(vol_file) )then
        write(logfhandle,'(a)') '    FAIL: '//vol_file%to_char()//' not found'
        THROW_HARD('TEST_REPROJECT FAILED: volume not generated')
    else
        call find_ldim_nptcls(vol_file, ldim, nptcls_stk)
        smpd_stk = find_img_smpd(vol_file)
        write(logfhandle,'(a,i4,a,i4,a,i4,a,f6.2)') '    PASS: volume dims = [', &
            ldim(1),',',ldim(2),',',ldim(3),' ], smpd = ', smpd_stk
        if( ldim(1) /= ldim(2) .or. ldim(1) < 1 )then
            write(logfhandle,'(a)') '    FAIL: volume has invalid dimensions'
            all_ok = .false.
        endif
        if( abs(smpd_stk - SMPD) > 0.01 )then
            write(logfhandle,'(a,f6.2,a,f6.2)') '    FAIL: smpd mismatch, expected ', SMPD, ' got ', smpd_stk
            all_ok = .false.
        endif
    endif
    ! ---- run reproject ----
    write(logfhandle,'(a)') '>>> TEST_REPROJECT: generating reprojections'
    call cline_reproj%set('prg',               'reproject')
    call cline_reproj%set('vol1',              '6VXX.mrc')
    call cline_reproj%set('smpd',                    SMPD)
    call cline_reproj%set('pgrp',                    'c1')
    call cline_reproj%set('mskdiam',                 180.)
    call cline_reproj%set('nspace',           real(NSPACE))
    call cline_reproj%set('nthr',                     16.)
    call xreproject%execute(cline_reproj)
    call cline_reproj%kill()
    ! ---- validate output stack ----
    write(logfhandle,'(a)') '>>> CHECK: output reprojection stack'
    if( .not. file_exists(outstk) )then
        write(logfhandle,'(a)') '    FAIL: '//outstk%to_char()//' not found'
        all_ok = .false.
    else
        call find_ldim_nptcls(outstk, ldim, nptcls_stk)
        smpd_stk = find_img_smpd(outstk)
        write(logfhandle,'(a,i6)')       '    projections in stack: ', nptcls_stk
        write(logfhandle,'(a,i4,a,i4)')  '    box size:             ', ldim(1), ' x ', ldim(2)
        write(logfhandle,'(a,f6.2)')     '    smpd:                 ', smpd_stk
        if( nptcls_stk /= NSPACE )then
            write(logfhandle,'(a,i6,a,i6)') '    FAIL: expected ', NSPACE, ' projections, got ', nptcls_stk
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: projection count matches'
        endif
        if( ldim(1) /= ldim(2) .or. ldim(1) < 1 )then
            write(logfhandle,'(a)') '    FAIL: invalid box dimensions'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: box dimensions valid'
        endif
        if( abs(smpd_stk - SMPD) > 0.01 )then
            write(logfhandle,'(a,f6.2,a,f6.2)') '    FAIL: smpd mismatch, expected ', SMPD, ' got ', smpd_stk
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: sampling distance matches'
        endif
    endif
    ! ---- validate orientations file ----
    write(logfhandle,'(a)') '>>> CHECK: orientations file'
    if( .not. file_exists(outori) )then
        write(logfhandle,'(a)') '    FAIL: '//outori%to_char()//' not found'
        all_ok = .false.
    else
        nlines_ori = nlines(outori)
        write(logfhandle,'(a,i6)') '    orientation records: ', nlines_ori
        if( nlines_ori /= NSPACE )then
            write(logfhandle,'(a,i6,a,i6)') '    FAIL: expected ', NSPACE, ' records, got ', nlines_ori
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: orientation count matches'
        endif
    endif
    ! ---- final verdict ----
    if( all_ok )then
        call simple_end('**** SIMPLE_TEST_REPROJECT NORMAL STOP ****')
    else
        THROW_HARD('TEST_REPROJECT FAILED')
    endif
end subroutine exec_test_reproject

subroutine exec_test_simulated_workflow( self, cline )
    use simple_atoms,         only: atoms
    use simple_molecule_data, only: molecule_data, betagal_1jyx, sars_cov2_spkgp_6vxx
    use simple_string_utils,  only: lowercase
    use simple_ui,            only: make_ui
    class(commander_test_simulated_workflow), intent(inout) :: self
    class(cmdline),                           intent(inout) :: cline
    character(len=*), parameter :: PROJNAME        = 'simulated_workflow'
    character(len=*), parameter :: PROJFILE        = PROJNAME//'.simple'
    character(len=*), parameter :: MOVIE_FILE      = 'simulate_movie.mrc'
    character(len=*), parameter :: SUBSET_FILE     = 'random_reprojections.mrcs'
    character(len=*), parameter :: OPTIMAL_FILE    = 'optimal_movie_average.mrc'
    character(len=*), parameter :: PARAMS_FILE     = 'simulate_movie_params.txt'
    character(len=*), parameter :: WORKFLOW_PICKER = 'segdiam'
    ! Reference-picker alternative: comment the line above and uncomment the line below.
    ! character(len=*), parameter :: WORKFLOW_PICKER = 'new'
    character(len=*), parameter :: PICKREFS_FILE  = 'pickrefs.mrc'
    character(len=*), parameter :: FILETAB_FILE   = 'simulated_movies.txt'
    character(len=*), parameter :: VOL_DIR        = '0_pdb2mrc'
    character(len=*), parameter :: REPROJ_DIR     = '1_reproject'
    character(len=*), parameter :: IMPORT_DIR     = '1_import_movies'
    character(len=*), parameter :: MOTION_DIR     = '2_motion_correct'
    character(len=*), parameter :: CTF_DIR        = '3_ctf_estimate'
    character(len=*), parameter :: PICK_DIR       = '4_pick'
    character(len=*), parameter :: EXTRACT_DIR    = '5_extract'
    character(len=*), parameter :: ABINIT2D_DIR   = '6_abinitio2D'
    real,             parameter :: SMPD           = 1.3
    real,             parameter :: MSKDIAM        = 180.0
    real,             parameter :: CS             = 2.7
    real,             parameter :: KV             = 300.0
    real,             parameter :: FRACA          = 0.1
    integer,          parameter :: NPROJS         = 100
    integer,          parameter :: NPER_MOVIE     = 12
    integer,          parameter :: MOVIE_DIM      = 1536
    integer,          parameter :: NFRAMES        = 16
    integer,          parameter :: NMOVIES        = 10
    integer,          parameter :: EXTRACT_BOX    = 192
    integer,          parameter :: NTHR           = 4
    type(cmdline)                       :: cline_projection, cline_sim_mov, cline_new_project
    type(cmdline)                       :: cline_import_movies, cline_mot_corr, cline_ctf_est
    type(cmdline)                       :: cline_pick, cline_extract, cline_abinitio2D, cline_abinitio3D
    type(commander_new_project)         :: xnew_project
    type(commander_reproject)           :: xreproject
    type(commander_simulate_movie)      :: xsimov
    type(commander_motion_correct)      :: xmotcorr
    type(commander_ctf_estimate)        :: xctf_estimate
    type(commander_import_movies)       :: ximport_movies
    type(commander_pick)                :: xpick
    type(commander_extract)             :: xextract
    type(commander_abinitio2D)          :: xabinitio2D
    type(commander_abinitio3D)          :: xabinitio3D
    type(molecule_data)                 :: mol
    type(atoms)                         :: molecule
    type(image)                         :: projection
    type(sp_project)                    :: spproj
    type(string)                        :: cwd_root, workflow_root, project_path, reproj_path, subset_path, filetab_path
    type(string)                        :: system_name, pgrp, test_workdir, vol_file, reproj_file
    type(string)                        :: movie_fname, subset_fname, optimal_fname, params_fname
    type(string)                        :: movie_files(NMOVIES)
    character(len=XLONGSTRLEN)          :: workflow_root_path
    integer                             :: i, j, iproj, proj_inds(NPER_MOVIE), ldim(3), nprojs_stk
    integer                             :: reproj_box, npickrefs, nptcls, ncls, status
    real                                :: pickref_smpd, pickref_width
    integer                             :: rnd_defocus

    ! The test executable initializes only the test UI.  Directly invoked SIMPLE
    ! commanders need the regular UI metadata to implement mkdir=yes correctly.
    call make_ui
    if( .not. cline%defined('system') ) THROW_HARD('The system keyword is required; expected 6vxx or 1jxy')
    system_name = cline%get_carg('system')
    system_name = lowercase(system_name%to_char())
    select case(system_name%to_char())
        case('6vxx')
            vol_file    = '6VXX.mrc'
            reproj_file = 'reprojs_6VXX.mrcs'
            pgrp        = 'c3'
        case('1jxy')
            vol_file    = '1JXY.mrc'
            reproj_file = 'reprojs_1JXY.mrcs'
            pgrp        = 'c1'
        case default
            THROW_HARD('Unsupported simulated-workflow system: '//system_name%to_char()//'; expected 6vxx or 1jxy')
    end select
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
    call molecule%pdb2mrc(volfile=vol_file, smpd=SMPD, mol=mol, center_pdb=.true.)
    call molecule%kill()
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
        &' simulated movies with ', NFRAMES, ' frames and ', NPER_MOVIE, ' random projections each'
    do i = 1,NMOVIES
        proj_inds = 0
        do j = 1,NPER_MOVIE
            do
                iproj = irnd_uni(NPROJS)
                if( j > 1 )then
                    if( any(proj_inds(:j-1) == iproj) ) cycle
                endif
                exit
            enddo
            proj_inds(j) = iproj
        enddo
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

    write(logfhandle,'(a,a)') '>>> Step 7: particle picking with ', WORKFLOW_PICKER
    call cline_pick%set('prg',                            'pick')
    call cline_pick%set('projfile',                 project_path)
    call cline_pick%set('mkdir',                           'yes')
    call cline_pick%set('pcontrast',                     'black')
    select case(WORKFLOW_PICKER)
        case('segdiam')
            call cline_pick%set('picker',              'segdiam')
            call cline_pick%set('moldiam_max',           MSKDIAM)
        case('new')
            call cline_pick%set('picker',                  'new')
            call cline_pick%set('pickrefs',          reproj_path)
            call cline_pick%set('moldiam',               MSKDIAM)
            call cline_pick%set('pick_roi',                 'no')
        case default
            THROW_HARD('Unsupported simulated-workflow picker: '//WORKFLOW_PICKER)
    end select
    call cline_pick%set('nparts',                              1)
    call cline_pick%set('nthr',                             NTHR)
    call xpick%execute(cline_pick)
    call cline_pick%kill()
    if( WORKFLOW_PICKER == 'new' )then
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

    write(logfhandle,'(a)') '>>> Step 9: ab initio 2D'
    call cline_abinitio2D%set('prg',                'abinitio2D')
    call cline_abinitio2D%set('projfile',           project_path)
    call cline_abinitio2D%set('mkdir',                     'yes')
    call cline_abinitio2D%set('mskdiam',                 MSKDIAM)
    call cline_abinitio2D%set('ncls',                       ncls)
    call cline_abinitio2D%set('nthr',                       NTHR)
    call xabinitio2D%execute(cline_abinitio2D)
    call cline_abinitio2D%kill()
    call update_project_path
    call return_to_stage_root('abinitio2D')

    call spproj%read(project_path)
    if( spproj%os_cls2D%get_noris() < 1 ) THROW_HARD('Ab initio 2D produced no classes')
    call spproj%kill()

    write(logfhandle,'(a)') '>>> Step 10: ab initio 3D'
    call cline_abinitio3D%set('prg',                'abinitio3D')
    call cline_abinitio3D%set('projfile',           project_path)
    call cline_abinitio3D%set('mkdir',                     'yes')
    call cline_abinitio3D%set('pgrp',                      pgrp)
    call cline_abinitio3D%set('mskdiam',                 MSKDIAM)
    call cline_abinitio3D%set('nthr',                       NTHR)
    call xabinitio3D%execute(cline_abinitio3D)
    call cline_abinitio3D%kill()

    call simple_chdir(cwd_root, status)
    if( status /= 0 ) THROW_HARD('Could not restore the original working directory')
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

!>  \brief  Integration test: split project into subprojects, run in parallel, merge back.
!  This test exercises the new generate_scripts_subprojects / gen_subproject_scripts_and_schedule
!  machinery end-to-end.
!    1. Generates a small synthetic particle stack and imports it into a project
!    2. Splits the project into deterministically balanced subprojects
!    3. Runs one cluster2D iteration on each subproject in parallel
!    4. Merges all subproject results with merge_chunk_projfiles
!    5. Validates that particle count is preserved
subroutine exec_test_subproject_distr( self, cline )
    use simple_commanders_project_ptcl, only: commander_import_particles
    class(commander_test_subproject_distr), intent(inout) :: self
    class(cmdline),                         intent(inout) :: cline
    integer,          parameter :: NPTCLS_SIM  = 32
    integer,          parameter :: NCLS_COARSE = 4
    integer,          parameter :: BOX_SIM     = 32
    integer,          parameter :: NTHR_JOB    = 1
    integer,          parameter :: MAXKEYS     = 20
    real,             parameter :: SMPD        = 1.3
    real,             parameter :: MSKDIAM     = 24.0
    character(len=*), parameter :: PROJNAME = 'test_subproj_distr'
    type(parameters)                    :: params
    type(sp_project)                    :: spproj, spproj_sub, spproj_merged
    type(commander_simulate_noise)      :: xsim_noise
    type(commander_new_project)         :: xnew_project
    type(commander_import_particles)    :: ximport_particles
    type(cmdline)                       :: cline_sim, cline_new_proj, cline_import, cline_sub
    type(qsys_env)                      :: qenv
    type(chash),  allocatable           :: jobs_descr(:)
    type(string), allocatable           :: subproj_fnames(:), subproj_dirs(:)
    integer,      allocatable           :: class_labels(:)
    real,         allocatable           :: state_labels(:)
    type(string)                        :: projname_sub, projfile_sub, cwd_root, sim_stk
    integer :: icls, iptcl, nsub, isub, nptcls_sub, nptcls_merged, nptcls_orig
    logical :: all_ok
    call cline%set('ncunits', NCLS_COARSE)
    call params%new(cline)
    call simple_getcwd(cwd_root)
    ! 1. Generate a deliberately small stack; molecular simulation and abinitio2D
    ! are unrelated to the subproject scheduling/merge behavior under test.
    write(logfhandle,'(a)') '>>> Step 1: generate compact synthetic particle stack'
    sim_stk = 'simulated_particles.mrc'
    call cline_sim%set('prg',        'simulate_noise')
    call cline_sim%set('mkdir',                  'no')
    call cline_sim%set('box',                 BOX_SIM)
    call cline_sim%set('smpd',                   SMPD)
    call cline_sim%set('nptcls',           NPTCLS_SIM)
    call cline_sim%set('outstk',              sim_stk)
    call xsim_noise%execute(cline_sim)
    call cline_sim%kill()
    ! 2. Create project and import synthetic particles without CTF metadata.
    write(logfhandle,'(a)') '>>> Step 2: create project & import particles'
    call cline_new_proj%set('projname', PROJNAME)
    call xnew_project%execute(cline_new_proj)
    call cline_new_proj%kill()
    ! import_particles (new_project changed cwd into project dir)
    call cline_import%set('prg',             'import_particles')
    call cline_import%set('mkdir',                         'no')
    call cline_import%set('projfile',       PROJNAME//'.simple')
    call cline_import%set('stk',       '../'//sim_stk%to_char())
    call cline_import%set('smpd',                          SMPD)
    call cline_import%set('ctf',                           'no')
    call ximport_particles%execute(cline_import)
    call cline_import%kill()
    write(logfhandle,'(a)') '    project populated with '//int2str(NPTCLS_SIM)//' particles'
    ! update cwd_root to project directory (new_project changed cwd)
    call simple_getcwd(cwd_root)
    ! 3. Split deterministically so setup cost and class populations are stable.
    write(logfhandle,'(a)') '>>> Step 3: split into balanced subprojects'
    call spproj%read(string(PROJNAME//'.simple'))
    nptcls_orig = spproj%get_nptcls()
    nsub        = NCLS_COARSE
    allocate(class_labels(nptcls_orig))
    do iptcl = 1, nptcls_orig
        class_labels(iptcl) = mod(iptcl - 1, nsub) + 1
    enddo
    allocate(jobs_descr(nsub), subproj_fnames(nsub), subproj_dirs(nsub))
    do icls = 1, nsub
        projname_sub = 'subproj_'//int2str_pad(icls, 2)
        projfile_sub = projname_sub%to_char()//'.simple'
        ! create subproject directory
        call simple_mkdir(projname_sub)
        subproj_dirs(icls) = cwd_root%to_char()//'/'//projname_sub%to_char()
        ! build subproject: copy parent, keep only particles in this class
        spproj_sub = spproj 
        ! set state labels: 1 for particles in this class, 0 otherwise
        allocate(state_labels(size(class_labels)))
        where(class_labels == icls)
            state_labels = 1.0
        elsewhere
            state_labels = 0.0
        endwhere
        call spproj_sub%os_ptcl2D%set_all('state', state_labels)
        call spproj_sub%os_ptcl3D%set_all('state', state_labels)
        call spproj_sub%prune_particles
        deallocate(state_labels)
        ! reset 2D clustering so cluster2D starts fresh (avoids early return
        ! in prep_strategy2D_glob that skips class_space_corrs allocation)
        call spproj_sub%os_ptcl2D%delete_2Dclustering
        call spproj_sub%os_cls2D%new(0, is_ptcl=.false.)
        nptcls_sub = spproj_sub%get_nptcls()
        write(logfhandle,'(a,i2,a,i6,a)') '    subproject ', icls, ': ', nptcls_sub, ' particles'
        ! update project info & write into subproject directory
        call cline_sub%set('projname', projname_sub)
        call spproj_sub%update_projinfo(cline_sub)
        call spproj_sub%write(string(subproj_dirs(icls)%to_char()//'/'//projfile_sub%to_char()))
        subproj_fnames(icls) = subproj_dirs(icls)%to_char()//'/'//projfile_sub%to_char()
        ! build job description for this subproject (no nparts => shared-memory execution)
        call jobs_descr(icls)%new(MAXKEYS)
        call jobs_descr(icls)%set('prg',                   'cluster2D_distr')
        call jobs_descr(icls)%set('projfile',         projfile_sub%to_char())
        call jobs_descr(icls)%set('ncls',     int2str(max(2, nptcls_sub/50)))
        call jobs_descr(icls)%set('nthr',                  int2str(NTHR_JOB))
        call jobs_descr(icls)%set('mskdiam',               real2str(MSKDIAM))
        call jobs_descr(icls)%set('maxits',                              '1')
        call jobs_descr(icls)%set('mkdir',                              'no')
        call jobs_descr(icls)%set('objfun',                             'cc')
        call spproj_sub%kill
    end do
    ! 4. Generate scripts and execute subprojects in parallel ----
    write(logfhandle,'(a)') '>>> Step 4: parallel execution of subprojects'
    params%projfile = PROJNAME//'.simple'
    params%nptcls   = nptcls_orig
    params%ncunits  = nsub
    params%nthr     = NTHR_JOB
    params%mskdiam  = MSKDIAM
    params%smpd     = SMPD
    call qenv%new(params, nsub)
    call qenv%gen_subproject_scripts_and_schedule(jobs_descr, subproj_dirs=subproj_dirs)
    write(logfhandle,'(a)') '    all subprojects completed'
    ! 5. Merge subproject results
    write(logfhandle,'(a)') '>>> Step 5: merging subproject results'
    call merge_chunk_projfiles(subproj_fnames, cwd_root, spproj_merged)
    nptcls_merged = spproj_merged%get_nptcls()
    ! 6. Validate
    write(logfhandle,'(a)') '>>> Step 6: validation'
    all_ok = .true.
    write(logfhandle,'(a,i6)') '    original particles:  ', nptcls_orig
    write(logfhandle,'(a,i6)') '    merged particles:    ', nptcls_merged
    if( nptcls_merged /= nptcls_orig )then
        write(logfhandle,'(a)') '    FAIL: particle count mismatch after merge!'
        all_ok = .false.
    else
        write(logfhandle,'(a)') '    PASS: particle count preserved'
    endif
    ! 7.Cleanup
    do isub = 1, nsub
        call jobs_descr(isub)%kill
    end do
    deallocate(jobs_descr, subproj_fnames, subproj_dirs)
    call spproj%kill
    call spproj_merged%kill
    call qenv%kill
    if( all_ok )then
        call simple_end('**** TEST_SUBPROJECT_DISTR NORMAL STOP ****')
    else
        THROW_HARD('TEST_SUBPROJECT_DISTR FAILED')
    endif
end subroutine exec_test_subproject_distr

subroutine exec_test_ptcls_ppca_subproject_distr( self, cline )
    class(commander_test_ptcls_ppca_subproject_distr), intent(inout) :: self
    class(cmdline),                                    intent(inout) :: cline
    integer,          parameter :: MAXKEYS = 20
    type(parameters)            :: params
    type(qsys_env)              :: qenv
    type(commander_stack)       :: xstack
    type(cmdline)               :: cline_stack
    type(chash),  allocatable   :: jobs_descr(:)
    type(string), allocatable   :: movie_fnames(:), chunk_fnames(:)
    type(string), allocatable   :: subproj_dirs(:), denoised_all(:)
    type(string)                :: cwd_root, chunk_filetab, chunk_stk, denoised_stk
    character(len=8)            :: subid
    integer                     :: nptcls, nsub, isub, i, fromp, top, nchunk
    integer                     :: nptcls_in, ldim(3)
    call params%new(cline)
    call read_filetable(params%filetab, movie_fnames)
    nptcls = size(movie_fnames)
    if( nptcls < 1 )then
        THROW_HARD('empty filetab; test_ptcls_ppca_subproject_distr')
    endif
    nsub = min(nptcls, params%nparts)
    call simple_getcwd(cwd_root)
    allocate(jobs_descr(nsub), subproj_dirs(nsub), denoised_all(nsub))
    write(logfhandle,'(a)') '>>> Step 1: split movie list into equal chunks and generate chunk stacks'
    do isub = 1, nsub
        fromp  = ((isub - 1) * nptcls) / nsub + 1
        top    = (isub * nptcls) / nsub
        nchunk = top - fromp + 1
        write(subid,'(I8.8)') isub
        subproj_dirs(isub) = cwd_root%to_char()//'/'//'ptcls_denoise_subproj_'//trim(adjustl(subid))
        call simple_mkdir(subproj_dirs(isub))
        allocate(chunk_fnames(nchunk))
        do i = 1, nchunk
            chunk_fnames(i) = movie_fnames(fromp + i - 1)
        enddo
        chunk_filetab = subproj_dirs(isub)%to_char()//'/'//'chunk_ptcls_'//trim(adjustl(subid))//'.txt'
        chunk_stk     = subproj_dirs(isub)%to_char()//'/'//'chunk_ptcls_'//trim(adjustl(subid))//'.mrcs'
        denoised_stk  = 'ppca_denoised_chunk_'//trim(adjustl(subid))//'.mrcs'
        call write_filetable(chunk_filetab, chunk_fnames)
        call cline_stack%set('prg',           'stack')
        call cline_stack%set('mkdir',            'no')
        call cline_stack%set('filetab', chunk_filetab)
        call cline_stack%set('outstk',      chunk_stk)
        call cline_stack%set('smpd',      params%smpd)
        call xstack%execute(cline_stack)
        call cline_stack%kill()
        call find_ldim_nptcls(chunk_stk, ldim, nptcls_in)
        write(logfhandle,'(a,i4,a,i8,a,i8,a,i8)') '    subproject ', isub,&
            &': chunk imgs=', nptcls_in, ' box=', ldim(1), 'x', ldim(2)
        call jobs_descr(isub)%new(MAXKEYS)
        call jobs_descr(isub)%set('prg',                                      'ppca_denoise')
        call jobs_descr(isub)%set('mkdir',                                              'no')
        call jobs_descr(isub)%set('stk',       'chunk_ptcls_'//trim(adjustl(subid))//'.mrcs')
        call jobs_descr(isub)%set('outstk',                                     denoised_stk)
        call jobs_descr(isub)%set('smpd',                              real2str(params%smpd))
        call jobs_descr(isub)%set('nthr',                               int2str(params%nthr))
        call jobs_descr(isub)%set('neigs',                                             '160')
        call jobs_descr(isub)%set('pca_mode',                                         'ppca')
        denoised_all(isub)        = subproj_dirs(isub)%to_char()//'/'//denoised_stk%to_char()
        deallocate(chunk_fnames)
    enddo
    write(logfhandle,'(a)') '>>> Step 2: run ppca_denoise in parallel across subprojects'
    call params%new(cline)
    call qenv%new(params, nsub)
    call qenv%gen_subproject_scripts_and_schedule(jobs_descr, subproj_dirs=subproj_dirs)
    write(logfhandle,'(a)') '    all denoising subprojects completed'
    write(logfhandle,'(a)') '>>> Step 3: merge denoised chunk stacks'
    call write_filetable(string('ppca_denoised_chunks.txt'), denoised_all)
    call cline_stack%set('prg',                        'stack')
    call cline_stack%set('mkdir',                         'no')
    call cline_stack%set('filetab', 'ppca_denoised_chunks.txt')
    call cline_stack%set('outstk',    'ppca_denoised_all.mrcs')
    call cline_stack%set('smpd',                   params%smpd)
    call xstack%execute(cline_stack)
    call cline_stack%kill()
    do isub = 1, nsub
        call jobs_descr(isub)%kill
    end do
    deallocate(jobs_descr, subproj_dirs, denoised_all)
    call qenv%kill
    call simple_end('**** TEST_PTCLS_PPCA_SUBPROJECT_DISTR NORMAL STOP ****')
end subroutine exec_test_ptcls_ppca_subproject_distr
!>  \brief  Single gate for the PCG reconstruction operator and solver,
!  see doc/policies/reconstruct3D_pcg_policy.md.
!
!  Replaces the former pcg_recon_ctf_free / _ctf_hetero / _kernel / _deapod
!  tests. Those shared a phantom, an RNG seed and three near-identical stages;
!  the CTF-free normal-operator and recovery stages were strict special cases
!  (T_i = 1) of the heterogeneous ones and are not repeated here. The CTF-free
!  ADJOINT check is kept, because it is the cheap isolation gate: if stage 2
!  fails while stage 1 passes, the fault is in build_transfer rather than in
!  the gather/scatter pair.
!
!  Stages run in increasing cost and fail fast, so the first failure is the
!  most local one:
!
!    1. adjoint identity, T = 1              -- gather/scatter pair alone
!    2. adjoint identity, T = C*S/sqrt(s2)   -- adds build_transfer
!    3. normal operator: symmetry + PSD      -- the two properties CG requires
!    4. synthetic recovery                   -- end-to-end solve
!    5. kernel vs matrix-free equivalence    -- the approximation's error budget
!    6. kernel invariants + preconditioner   -- |T|^2 drops the shift, keeps CTF
!    7. streaming vs monolithic accumulation -- begin_accum/accumulate_batch/
!                                               end_accum/solve_accum must
!                                               reproduce solve() exactly
!    8. deapodization                        -- the only stage without an
!                                               inverse crime (see below)
!    9. symmetry by coordinate replication   -- in-operator point-group
!                                               replication must equal a c1
!                                               solve of the expanded set
!   10. lambda/data-mass scaling              -- duplicating all weighted data
!                                               must not change the solution
!   11. crop invariants                       -- common-band raw B,D and the
!                                               cropped kernel operator
!   12. ML prior operator                     -- positive FSC/SSNR diagonal and
!                                               kernel/matrix-free parity
!
!  INVERSE CRIME, deliberately. Stages 1-7 generate observations with
!  forward_plane, so the operator's own KB envelope appears identically in the
!  data and the model and cancels. That is what makes them clean gates on the
!  operator ALGEBRA, and it is also why they cannot see an envelope bias: they
!  scored 0.9998 correlation while the forward model was still wrong for real
!  particles. Stage 8 is the honest one and exists for exactly that reason.
!
!  Does not touch reconstructor, reconstructor_eo, or volassemble.
subroutine exec_test_pcg_recon( self, cline )
    use simple_reconstructor_pcg, only: reconstructor_pcg, pcg_solver_outcome, &
        &PCG_OP_MATRIXFREE, PCG_OP_KERNEL
    use simple_sym,               only: sym
    class(commander_test_pcg_recon), intent(inout) :: self
    class(cmdline),                  intent(inout) :: cline
    integer,          parameter :: BOX = 24, CROP_BOX = 16, NPROJS = 40, NBLOBS = 4, NCTF = 5
    integer,          parameter :: BATCHSZ = 7          ! deliberately not a divisor of NPROJS
    real,             parameter :: SMPD = 1.5, LAMBDA = 1.0e-3
    real,             parameter :: CROP_SMPD = SMPD * real(BOX) / real(CROP_BOX)
    real,             parameter :: MASS_LAMBDA_REL = 1.0e-3
    real,             parameter :: ADJOINT_RELTOL = 1.0e-5, NORMAL_OP_RELTOL = 1.0e-4
    real,             parameter :: RECON_CORR_THRES = 0.9, MONOTONIC_SLACK = 1.0e-3
    ! kernel epsilon set from the single-precision roundoff the operator stages
    ! measure (~1e-7 relative), loosened for the interior to allow for the
    ! kernel's shift-invariance approximation. NOT tuned after inspecting a
    ! reconstruction.
    real,             parameter :: EPS_INTERIOR = 5.0e-2
    ! calibrate_kernel now applies the ANALYTIC factor padsc**2 instead of
    ! fitting one. The fit is retained in measure_kernel_scale purely so this
    ! test can assert the constant has not drifted: it returns 1.0 exactly when
    ! the analytic value is right. Observed fits were within 0.5% of it.
    real,             parameter :: KSCALE_TOL = 2.0e-2
    ! Fused and monolithic accumulation differ at single-precision roundoff
    ! (~1e-7 observed for both B and Khat). Gate those statistics directly and
    ! strictly. The preconditioner is a guarded reciprocal of D, and 20 Krylov
    ! recurrences amplify that harmless perturbation; give the final volume a
    ! separate 5e-4 bound so solver conditioning cannot hide an accumulator
    ! defect or manufacture a false failure.
    real,             parameter :: STREAM_ACCUM_RELTOL = 1.0e-6
    real,             parameter :: STREAM_SOLVE_RELTOL = 5.0e-4
    real,             parameter :: MASS_SCALE_RELTOL = 5.0e-6
    real,             parameter :: MASS_SOLVE_RELTOL = 5.0e-4
    real,             parameter :: CROP_RAW_RELTOL = 2.0e-5
    integer,          parameter :: STREAM_ITS    = 20
    integer,          parameter :: KERNEL_COMPARE_ITS = 8
    real,             parameter :: CTRS(3,NBLOBS) = reshape([&
        &-5.0,-3.0, 2.0,&
        &4.0, 5.0,-3.0,&
        &0.0,-6.0,-5.0,&
        &3.0,-2.0, 6.0], [3,NBLOBS])
    real,             parameter :: SIGMAS(NBLOBS) = [2.0, 2.5, 1.8, 2.2]
    real,             parameter :: AMPS(NBLOBS)   = [1.0, 0.8, 0.6, 0.5]
    real,             parameter :: KV = 300., CS = 2.7, FRACA = 0.1
    real,             parameter :: DFX_VALS(NCTF)    = [1.0, 1.5, 2.0, 2.5, 3.0]
    real,             parameter :: ASTIG_VALS(NCTF)  = [0.10, 0.15, 0.20, 0.12, 0.18]
    real,             parameter :: ANGAST_VALS(NCTF) = [0., 20., 40., 60., 80.]
    type(reconstructor_pcg) :: pcgop, pcg_reduce, pcg_crop, pcg_ml
    type(oris)              :: projdirs, projdirs_exp, projdirs_crop
    type(ori)               :: e, e_exp
    type(ctfparams)         :: ctfparms
    type(sym)               :: c1sym, c2sym
    type(string)            :: raw_part1, raw_part2, raw_ml
    type(pcg_solver_outcome) :: solver_outcome
    real,    allocatable    :: phantom(:,:,:), p_probe(:,:,:), q_probe(:,:,:)
    real,    allocatable    :: hp(:,:,:), hq(:,:,:), hm(:,:,:), hk(:,:,:)
    real,    allocatable    :: recon(:,:,:), recon_str(:,:,:), rel_res_hist(:)
    real,    allocatable    :: recon_mf(:,:,:), recon_kernel(:,:,:)
    real,    allocatable    :: recon_mass(:,:,:), recon_mass_dup(:,:,:)
    real,    allocatable    :: p_crop(:,:,:), hm_crop(:,:,:), hk_crop(:,:,:)
    real,    allocatable    :: hm_ml(:,:,:), hk_ml(:,:,:), ml_prior_diag(:,:,:)
    real,    allocatable    :: xdiv(:,:,:), recon_on(:,:,:), recon_off(:,:,:)
    real,    allocatable    :: env(:,:,:), invenv(:,:,:)
    real,    allocatable    :: khat_a(:,:,:), khat_b(:,:,:), b_mono(:,:,:), b_str(:,:,:)
    real,    allocatable    :: qplane_re(:,:), qplane_im(:,:), sig2arr(:), sig2_2d(:,:)
    real,    allocatable    :: sig2_crop(:,:), draw_full(:,:,:), draw_crop(:,:,:), fsc_prior(:)
    complex, allocatable    :: gx_plane(:,:), qplane(:,:), mplane(:,:), wplane(:,:), Ti(:,:)
    complex, allocatable    :: adj_out(:,:,:), y_planes(:,:,:), y_exp(:,:,:), y_crop(:,:,:)
    complex, allocatable    :: braw_full(:,:,:), braw_crop(:,:,:)
    integer :: lims2(2,2), lims3(3,2), i, j, k, b, g, c, niters, iseed_n
    integer :: R, margin, lo, hi, ifrom, nb, nsym, nraw, nraw_total, ml_prior_npositive
    integer :: lims2_crop(2,2), lims3_crop(3,2), raw_lim, hraw, kraw, mraw
    integer, allocatable :: iseed(:)
    real    :: ctr, dx, dy, dz, adjoint_err, corr, shift(2)
    real    :: err_all, err_int, err_max, den_all, den_int, kdiff, stream_err, rhs_err, kscale
    real    :: solution_err, solution_norm_ratio, energy_ratio
    real    :: corr_on, corr_off, env_ctr, env_edge, recip_err
    real    :: data_scale, data_scale_dup, lambda_eff, lambda_eff_dup
    real    :: mass_scale_err, mass_lambda_err, mass_solution_err
    real    :: crop_b_err, crop_d_err, crop_kernel_err, crop_factor
    real    :: ml_kernel_err, ml_prior_min, ml_prior_max
    real    :: ml_prior_positive_min, ml_prior_positive_max, ml_prior_to_khat_l1, ml_prior_to_khat_rms
    real(dp):: lhs, rhs, dp_p_hq, dp_hp_q, dp_p_hp
    real(dp):: crop_b_num, crop_b_den, crop_d_num, crop_d_den, ml_prior_energy
    logical :: all_ok
    all_ok = .true.

    ! ---- deterministic RNG seed: every probe below must be reproducible or a
    !      tolerance failure cannot be told from a different random draw ----
    call random_seed(size=iseed_n)
    allocate(iseed(iseed_n), source=42)
    call random_seed(put=iseed)

    ! ---- deterministic, asymmetric phantom (sum of off-centre Gaussian blobs).
    !      Asymmetric on purpose: a symmetric phantom hides orientation bugs. ----
    write(logfhandle,'(a)') '>>> TEST_PCG_RECON: building deterministic phantom'
    allocate(phantom(BOX,BOX,BOX), source=0.0)
    ctr = real(BOX)/2.0 + 0.5
    do k = 1,BOX
        do j = 1,BOX
            do i = 1,BOX
                do b = 1,NBLOBS
                    dx = real(i)-ctr-CTRS(1,b); dy = real(j)-ctr-CTRS(2,b); dz = real(k)-ctr-CTRS(3,b)
                    phantom(i,j,k) = phantom(i,j,k) + AMPS(b)*exp(-(dx*dx+dy*dy+dz*dz)/(2.0*SIGMAS(b)**2))
                end do
            end do
        end do
    end do

    call pcgop%new(BOX, SMPD, LAMBDA)
    ! Deapodization OFF for stages 1-7 -- see the inverse-crime note in the
    ! header. Stage 8 turns it back on and supplies envelope-free data.
    call pcgop%set_deapod(.false.)
    lims2 = pcgop%get_lims2()
    lims3 = pcgop%get_lims3()
    R     = lims2(1,2)
    allocate(sig2arr(0:R))
    do i = 0, R
        sig2arr(i) = 1.0 + 0.15*real(i)
    end do
    call e%new(.false.)

    ! ================= STAGE 1: adjoint identity, T = 1 =================
    write(logfhandle,'(a)') '>>> STAGE 1: adjoint dot-product identity (T = 1)'
    call e%set_euler([30.,55.,70.])
    allocate(gx_plane(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2)))
    call pcgop%set_volume(phantom)
    call pcgop%forward_plane(e, gx_plane)
    allocate(qplane_re(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2)))
    allocate(qplane_im(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2)))
    call random_number(qplane_re); call random_number(qplane_im)
    allocate(qplane(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2)))
    qplane = cmplx(qplane_re-0.5, qplane_im-0.5)
    lhs = sum(real(conjg(gx_plane)*qplane, dp))
    allocate(adj_out(lims3(1,1):lims3(1,2), lims3(2,1):lims3(2,2), lims3(3,1):lims3(3,2)), source=cmplx(0.,0.))
    call pcgop%adjoint_plane_add(qplane, e, adj_out)
    ! <x, G^dagger q> over the operator's own oversampled lattice; fourier_dot
    ! keeps that lattice and its Friedel packing an implementation detail
    rhs = pcgop%fourier_dot(adj_out)
    adjoint_err = real(abs(lhs-rhs) / max(1.0_dp, abs(lhs), abs(rhs)))
    write(logfhandle,'(a,es14.6,a,es14.6,a,es14.6)') '    <Gx,q>=', real(lhs), ' <x,G^Tq>=', real(rhs), &
        &' rel_err=', adjoint_err
    if( adjoint_err > ADJOINT_RELTOL )then
        write(logfhandle,'(a)') '    FAIL: adjoint dot-product identity violated'
        all_ok = .false.
    else
        write(logfhandle,'(a)') '    PASS: adjoint dot-product identity holds'
    endif

    ! ============ STAGE 2: adjoint identity, weighted transfer ============
    ! The T = 1 case above is trivially self-adjoint and cannot exercise
    ! build_transfer; this one carries a real astigmatic CTF, a nonzero shift
    ! and a sigma2 profile.
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 2: adjoint identity with T = C*S/sqrt(sigma2)'
        ctfparms%smpd = SMPD; ctfparms%kv = KV; ctfparms%cs = CS; ctfparms%fraca = FRACA
        ctfparms%dfx = DFX_VALS(1); ctfparms%dfy = DFX_VALS(1)+ASTIG_VALS(1)
        ctfparms%angast = ANGAST_VALS(1); ctfparms%phshift = 0.
        shift = [1.7, -2.3]
        allocate(Ti(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2)))
        Ti = pcgop%build_transfer(ctfparms, shift, sig2arr)
        allocate(mplane(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2)))
        mplane = Ti * gx_plane
        lhs = sum(real(conjg(mplane)*qplane, dp))
        allocate(wplane(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2)))
        wplane  = conjg(Ti) * qplane
        adj_out = cmplx(0.,0.)
        call pcgop%adjoint_plane_add(wplane, e, adj_out)
        rhs = pcgop%fourier_dot(adj_out)
        adjoint_err = real(abs(lhs-rhs) / max(1.0_dp, abs(lhs), abs(rhs)))
        write(logfhandle,'(a,es14.6,a,es14.6,a,es14.6)') '    <T*Gx,q>=', real(lhs), &
            &' <x,G^T(conjg(T)*q)>=', real(rhs), ' rel_err=', adjoint_err
        if( adjoint_err > ADJOINT_RELTOL )then
            write(logfhandle,'(a)') '    FAIL: weighted adjoint identity violated (build_transfer)'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: weighted adjoint identity holds'
        endif
    else
        write(logfhandle,'(a)') '>>> STAGE 2 SKIPPED: stage 1 failed'
    endif

    ! ---- heterogeneous selection used by stages 3-7 ----
    call projdirs%new(NPROJS, .false.)
    call projdirs%spiral
    do i = 1, NPROJS
        call projdirs%get_ori(i, e)
        g = mod(i-1, NCTF) + 1
        ctfparms%smpd = SMPD; ctfparms%kv = KV; ctfparms%cs = CS; ctfparms%fraca = FRACA
        ctfparms%dfx = DFX_VALS(g); ctfparms%dfy = DFX_VALS(g)+ASTIG_VALS(g)
        ctfparms%angast = ANGAST_VALS(g); ctfparms%phshift = 0.
        call e%set_ctfvars(ctfparms)
        call e%set_shift([1.0+0.3*real(mod(i-1,5)), -1.0+0.2*real(mod(i-1,7))])
        call projdirs%set_ori(i, e)
    end do
    allocate(sig2_2d(0:R,NPROJS))
    do i = 1, NPROJS
        sig2_2d(:,i) = sig2arr
    end do

    ! ========== STAGE 3: normal operator, symmetry and PSD ==========
    ! CG requires both. Symmetry is checked by the dot-product identity across
    ! two independent random probes; PSD by a single quadratic form.
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 3: normal-operator symmetry and positive-definiteness'
        call pcgop%prep_particles(projdirs, use_ctf=.true., sig2=sig2_2d)
        allocate(p_probe(BOX,BOX,BOX), q_probe(BOX,BOX,BOX))
        call random_number(p_probe); call random_number(q_probe)
        p_probe = p_probe - 0.5; q_probe = q_probe - 0.5
        hp = pcgop%apply_normal(p_probe)
        hq = pcgop%apply_normal(q_probe)
        dp_p_hq = pcgop%dot_real_volume(p_probe, hq)
        dp_hp_q = pcgop%dot_real_volume(hp, q_probe)
        dp_p_hp = pcgop%dot_real_volume(p_probe, hp)
        adjoint_err = real(abs(dp_p_hq-dp_hp_q) / max(1.0_dp, abs(dp_p_hq), abs(dp_hp_q)))
        write(logfhandle,'(a,es14.6,a,es14.6,a,es14.6)') '    dot(p,Hq)=', real(dp_p_hq), &
            &' dot(Hp,q)=', real(dp_hp_q), ' rel_err=', adjoint_err
        if( adjoint_err > NORMAL_OP_RELTOL )then
            write(logfhandle,'(a)') '    FAIL: normal operator is not symmetric'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: normal operator is symmetric'
        endif
        write(logfhandle,'(a,es14.6)') '    dot(p,Hp)=', real(dp_p_hp)
        if( dp_p_hp <= 0.0_dp )then
            write(logfhandle,'(a)') '    FAIL: normal operator is not positive-definite'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: normal operator is positive-definite'
        endif
    else
        write(logfhandle,'(a)') '>>> STAGE 3 SKIPPED: an earlier stage failed'
    endif

    ! ============== STAGE 4: heterogeneous synthetic recovery ==============
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 4: synthetic recovery through CTF + shift + sigma'
        allocate(y_planes(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), NPROJS))
        call pcgop%set_volume(phantom)
        do i = 1, NPROJS
            call projdirs%get_ori(i, e)
            call pcgop%forward_plane(e, gx_plane)
            Ti = pcgop%build_transfer(e%get_ctfvars(), e%get_2Dshift(), sig2arr)
            y_planes(:,:,i) = Ti * gx_plane
        end do
        ! build_operators here rather than a bare solve, so stage 7 can compare
        ! the two accumulation paths with identical preconditioner state
        call pcgop%build_operators(.false.)
        allocate(recon(BOX,BOX,BOX), source=0.0)
        call pcgop%solve(y_planes, recon, maxits=40, rtol=1.0e-3, &
            &rel_res_hist=rel_res_hist, niters=niters)
        write(logfhandle,'(a,i0,a)') '    PCG ran ', niters, ' iterations'
        do i = 2, niters
            if( rel_res_hist(i) > rel_res_hist(i-1) + MONOTONIC_SLACK )then
                write(logfhandle,'(a,i0)') '    FAIL: relative residual increased at iteration ', i
                all_ok = .false.
            endif
        end do
        corr = corr_of(recon, phantom)
        write(logfhandle,'(a,f8.5)') '    recon-vs-phantom correlation = ', corr
        if( corr < RECON_CORR_THRES )then
            write(logfhandle,'(a)') '    FAIL: recovered volume does not correlate with the known phantom'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: recovered volume correlates with the known phantom'
        endif
    else
        write(logfhandle,'(a)') '>>> STAGE 4 SKIPPED: an earlier stage failed'
    endif

    ! ===== STAGE 5: kernel vs matrix-free operator and solve baseline =====
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 5: kernelized vs matrix-free normal operator'
        call pcgop%build_kernel
        hm = pcgop%apply_normal_matrixfree(p_probe)
        hk = pcgop%apply_normal_kernel(p_probe)
        den_all = max(1.0, sqrt(sum(hm*hm)))
        err_all = sqrt(sum((hk-hm)**2)) / den_all
        ! interior margin at least the KB support, so the shift-invariance
        ! approximation is not judged on the boundary where it is known to differ
        margin  = 4
        lo      = margin + 1
        hi      = BOX - margin
        den_int = max(1.0, sqrt(sum(hm(lo:hi,lo:hi,lo:hi)**2)))
        err_int = sqrt(sum((hk(lo:hi,lo:hi,lo:hi)-hm(lo:hi,lo:hi,lo:hi))**2)) / den_int
        err_max = maxval(abs(hk-hm)) / max(1.0, maxval(abs(hm)))
        write(logfhandle,'(a,es14.6)') '    relative error, all voxels = ', err_all
        write(logfhandle,'(a,es14.6)') '    relative error, interior   = ', err_int
        write(logfhandle,'(a,es14.6)') '    relative maximum error     = ', err_max
        if( err_int > EPS_INTERIOR )then
            write(logfhandle,'(a)') '    FAIL: kernelized operator disagrees with the matrix-free reference'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: kernelized operator agrees with the reference in the interior'
        endif
        ! A scale error of a few percent would pass the interior check above but
        ! is exactly what a drift in the analytic constant would look like, so
        ! assert on it directly.
        kscale = pcgop%measure_kernel_scale()
        write(logfhandle,'(a,f10.6)') '    measured kernel scale (1.0 = analytic padsc**2 exact) = ', kscale
        if( abs(kscale - 1.0) > KSCALE_TOL )then
            write(logfhandle,'(a)') '    FAIL: kernel scale has drifted from the analytic constant'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: analytic kernel scale is correct'
        endif
        dp_p_hq = pcgop%dot_real_volume(p_probe, hm)
        dp_hp_q = pcgop%dot_real_volume(p_probe, hk)
        energy_ratio = real(dp_hp_q / max(abs(dp_p_hq), epsilon(1.0_dp)))
        write(logfhandle,'(a,es14.6)') '    kernel/matrix-free probe energy ratio = ', energy_ratio

        ! Establish the fixed-iteration reconstruction baseline without tuning
        ! an acceptance threshold before current Linux/macOS CI evidence exists.
        ! Same RHS, initial state and iteration count prevent convergence-stop
        ! jitter from masquerading as an operator difference.
        allocate(recon_mf(BOX,BOX,BOX), recon_kernel(BOX,BOX,BOX), source=0.0)
        call pcgop%set_op_mode(PCG_OP_MATRIXFREE)
        call pcgop%solve(y_planes, recon_mf, maxits=KERNEL_COMPARE_ITS, rtol=0.0)
        call pcgop%set_op_mode(PCG_OP_KERNEL)
        call pcgop%solve(y_planes, recon_kernel, maxits=KERNEL_COMPARE_ITS, rtol=0.0)
        solution_err = sqrt(sum((recon_kernel-recon_mf)**2)) / &
            &max(1.0, sqrt(sum(recon_mf*recon_mf)))
        solution_norm_ratio = sqrt(sum(recon_kernel*recon_kernel)) / &
            &max(1.0, sqrt(sum(recon_mf*recon_mf)))
        write(logfhandle,'(a,i0,a,es14.6)') '    fixed-', KERNEL_COMPARE_ITS, &
            &'-iteration solution rel_err = ', solution_err
        write(logfhandle,'(a,es14.6)') '    fixed-iteration solution norm ratio = ', solution_norm_ratio
        deallocate(recon_mf, recon_kernel)
    else
        write(logfhandle,'(a)') '>>> STAGE 5 SKIPPED: an earlier stage failed'
    endif

    ! ========= STAGE 6: kernel invariants and preconditioner =========
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 6: kernel invariants and preconditioner'
        khat_a = pcgop%apply_normal_kernel(p_probe)
        ! changing ONLY the shifts must leave the normal operator unchanged: the
        ! shift is a unit-modulus phase and cancels in |T|^2
        do i = 1, NPROJS
            call projdirs%get_ori(i, e)
            call e%set_shift([-2.0+0.11*real(i), 3.0-0.07*real(i)])
            call projdirs%set_ori(i, e)
        end do
        call pcgop%prep_particles(projdirs, use_ctf=.true., sig2=sig2_2d)
        call pcgop%build_kernel
        khat_b = pcgop%apply_normal_kernel(p_probe)
        kdiff  = sqrt(sum((khat_b-khat_a)**2)) / max(1.0, sqrt(sum(khat_a*khat_a)))
        write(logfhandle,'(a,es14.6)') '    relative change after shift-only edit = ', kdiff
        if( kdiff > 1.0e-5 )then
            write(logfhandle,'(a)') '    FAIL: kernel changed when only shifts changed'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: kernel is shift-invariant'
        endif
        ! changing a CTF MUST change it
        do i = 1, NPROJS
            call projdirs%get_ori(i, e)
            ctfparms = e%get_ctfvars()
            ctfparms%dfx = ctfparms%dfx + 0.7
            ctfparms%dfy = ctfparms%dfy + 0.7
            call e%set_ctfvars(ctfparms)
            call projdirs%set_ori(i, e)
        end do
        call pcgop%prep_particles(projdirs, use_ctf=.true., sig2=sig2_2d)
        call pcgop%build_kernel
        khat_b = pcgop%apply_normal_kernel(p_probe)
        kdiff  = sqrt(sum((khat_b-khat_a)**2)) / max(1.0, sqrt(sum(khat_a*khat_a)))
        write(logfhandle,'(a,es14.6)') '    relative change after CTF edit        = ', kdiff
        if( kdiff < 1.0e-4 )then
            write(logfhandle,'(a)') '    FAIL: kernel did not change when the CTF changed'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: kernel tracks the CTF'
        endif
        call pcgop%set_op_mode(PCG_OP_MATRIXFREE)
        call pcgop%build_precond
        hm = pcgop%apply_normal_matrixfree(p_probe)
        if( any(hm /= hm) )then
            write(logfhandle,'(a)') '    FAIL: operator produced non-finite values after build_precond'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: preconditioner built, operator still finite'
        endif
    else
        write(logfhandle,'(a)') '>>> STAGE 6 SKIPPED: an earlier stage failed'
    endif

    ! ================== STAGE 7: streaming accumulation ==================
    ! begin_accum/accumulate_batch/end_accum/solve_accum must reproduce what
    ! solve() produces from all planes at once. The accumulated B and D-derived
    ! kernel must agree near single-precision roundoff; the fixed-step solution
    ! has its own bound because the reciprocal preconditioner and CG recurrence
    ! amplify that input perturbation.
    ! BATCHSZ is deliberately not a divisor of NPROJS so the short final batch
    ! is exercised.
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 7: streaming accumulation vs monolithic solve'
        call pcgop%prep_particles(projdirs, use_ctf=.true., sig2=sig2_2d)
        ! regenerate observations: stage 6 edited the shifts and CTFs
        call pcgop%set_volume(phantom)
        do i = 1, NPROJS
            call projdirs%get_ori(i, e)
            call pcgop%forward_plane(e, gx_plane)
            Ti = pcgop%build_transfer(e%get_ctfvars(), e%get_2Dshift(), sig2arr)
            y_planes(:,:,i) = Ti * gx_plane
        end do
        ! BOTH paths run a FIXED number of iterations with the tolerance
        ! disabled. Comparing two tolerance-stopped solves is meaningless:
        ! round-off in the accumulators can put the stopping decision one
        ! iteration apart, and a single CG step near convergence moves the
        ! solution by O(rtol) -- orders of magnitude above anything worth
        ! asserting. Forcing identical steps leaves accumulator round-off as the
        ! only difference between the two paths, which is what is being tested.
        ! reference path: one shot
        call pcgop%build_operators(.false.)
        recon = 0.0
        call pcgop%solve(y_planes, recon, maxits=STREAM_ITS, rtol=0.0, &
            &niters=niters, outcome=solver_outcome)
        call pcgop%get_rhs(b_mono)
        write(logfhandle,'(a,i0,a)') '    monolithic: ', niters, ' iterations'
        if( trim(solver_outcome%stop_reason) /= 'fixed_iterations' .or. &
            &solver_outcome%iteration_count /= STREAM_ITS .or. solver_outcome%converged )then
            write(logfhandle,'(a)') '    FAIL: fixed-iteration solver outcome is inconsistent'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: fixed-iteration solver outcome is explicit'
        endif
        ! Two streaming runs, so a failure localizes itself rather than just
        ! reporting a number: ONE batch covering everything isolates the
        ! begin/end_accum machinery, and many batches adds the per-batch
        ! indexing on top. If single-batch agrees and multi-batch does not, the
        ! fault is in the batch range handling; if neither agrees, it is in the
        ! accumulate/fold path shared by both.
        allocate(recon_str(BOX,BOX,BOX))
        do j = 1, 2
            if( j == 1 )then
                nb = NPROJS                     ! single batch
            else
                nb = BATCHSZ                    ! many batches, short final one
            endif
            call pcgop%begin_accum
            do ifrom = 1, NPROJS, nb
                k = min(nb, NPROJS - ifrom + 1)
                call pcgop%accumulate_batch(y_planes(:,:,ifrom:ifrom+k-1), k, ifrom)
            end do
            call pcgop%end_accum(.false.)
            ! compare the RHS before the solve: if b already differs the fault is
            ! in fused B accumulation or the fold; if b agrees but the solutions
            ! do not, it is in the fused D accumulator feeding the preconditioner
            call pcgop%get_rhs(b_str)
            rhs_err = sqrt(sum((b_str-b_mono)**2)) / max(1.0, sqrt(sum(b_mono*b_mono)))
            recon_str = 0.0
            call pcgop%solve_accum(recon_str, maxits=STREAM_ITS, rtol=0.0, niters=niters)
            stream_err = sqrt(sum((recon_str-recon)**2)) / max(1.0, sqrt(sum(recon*recon)))
            write(logfhandle,'(a,i3,a,es14.6,a,es14.6)') '    batch size ', nb, &
                &': rel_err(b) = ', rhs_err, '   rel_err(x) = ', stream_err
            if( rhs_err > STREAM_ACCUM_RELTOL )then
                write(logfhandle,'(a)') '    FAIL: fused RHS does not reproduce monolithic accumulation'
                all_ok = .false.
            else if( stream_err > STREAM_SOLVE_RELTOL )then
                write(logfhandle,'(a)') '    FAIL: fused accumulation changes the fixed-step solve beyond tolerance'
                all_ok = .false.
            else
                write(logfhandle,'(a)') '    PASS: streaming accumulation reproduces the monolithic solve'
            endif
        end do
        ! The command's DEFAULT is pcgop=kernel, which reaches the kernel through
        ! end_accum(.true.) -- a branch the runs above never touch, since they
        ! pass .false. Deriving Khat from a streamed accumulator must match
        ! deriving it from a monolithic one, or the default path is untested.
        call pcgop%prep_particles(projdirs, use_ctf=.true., sig2=sig2_2d)
        call pcgop%build_kernel
        khat_a = pcgop%apply_normal_kernel(p_probe)
        do j = 1, 2
            if( j == 1 )then
                nb = NPROJS
            else
                nb = BATCHSZ
            endif
            call pcgop%begin_accum
            do ifrom = 1, NPROJS, nb
                k = min(nb, NPROJS - ifrom + 1)
                call pcgop%accumulate_batch(y_planes(:,:,ifrom:ifrom+k-1), k, ifrom)
            end do
            call pcgop%end_accum(.true.)
            khat_b = pcgop%apply_normal_kernel(p_probe)
            kdiff  = sqrt(sum((khat_b-khat_a)**2)) / max(1.0, sqrt(sum(khat_a*khat_a)))
            write(logfhandle,'(a,i3,a,es14.6)') '    kernel batch size ', nb, &
                &': streamed vs monolithic rel_err = ', kdiff
            if( kdiff > STREAM_ACCUM_RELTOL )then
                write(logfhandle,'(a)') '    FAIL: kernel built from a streamed accumulator differs'
                all_ok = .false.
            else
                write(logfhandle,'(a)') '    PASS: kernel is identical either way'
            endif
        end do
        ! Distributed contract: workers publish raw, unfolded B and D. The
        ! master adds artifacts in ascending part order and only then folds and
        ! finalizes. Split the same particle sequence into two artifacts so the
        ! gate covers serialization, fixed-order reduction and master-only
        ! finalization without involving a scheduler.
        raw_part1 = 'test_pcg_raw_part1.dat'
        raw_part2 = 'test_pcg_raw_part2.dat'
        call pcgop%begin_accum
        call pcgop%accumulate_batch(y_planes(:,:,1:NPROJS/2), NPROJS/2, 1)
        call pcgop%write_raw_accum(raw_part1, 1, 0, 1, 2, NPROJS/2, 'pcg_recon_test_v1')
        call pcgop%begin_accum
        call pcgop%accumulate_batch(y_planes(:,:,NPROJS/2+1:NPROJS), NPROJS/2, NPROJS/2+1)
        call pcgop%write_raw_accum(raw_part2, 1, 0, 2, 2, NPROJS/2, 'pcg_recon_test_v1')
        call pcgop%end_accum(.false.)
        call pcg_reduce%new(BOX, SMPD, LAMBDA)
        call pcg_reduce%set_deapod(.false.)
        call pcg_reduce%begin_reduction
        nraw_total = 0
        call pcg_reduce%add_raw_accum(raw_part1, 1, 0, 1, 2, 'pcg_recon_test_v1', nraw)
        nraw_total = nraw_total + nraw
        call pcg_reduce%add_raw_accum(raw_part2, 1, 0, 2, 2, 'pcg_recon_test_v1', nraw)
        nraw_total = nraw_total + nraw
        call pcg_reduce%end_accum(.true.)
        call pcg_reduce%get_rhs(b_str)
        khat_b = pcg_reduce%apply_normal_kernel(p_probe)
        rhs_err = sqrt(sum((b_str-b_mono)**2)) / max(1.0, sqrt(sum(b_mono*b_mono)))
        kdiff   = sqrt(sum((khat_b-khat_a)**2)) / max(1.0, sqrt(sum(khat_a*khat_a)))
        write(logfhandle,'(a,i0,a,es14.6,a,es14.6)') '    raw reduction particles ', nraw_total, &
            &': rel_err(b) = ', rhs_err, '   rel_err(kernel) = ', kdiff
        if( nraw_total /= NPROJS .or. rhs_err > STREAM_ACCUM_RELTOL .or. kdiff > STREAM_ACCUM_RELTOL )then
            write(logfhandle,'(a)') '    FAIL: raw fixed-order reduction differs from monolithic accumulation'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: raw fixed-order reduction reproduces monolithic accumulation'
        endif
        call pcg_reduce%kill
        call del_file(raw_part1)
        call del_file(raw_part2)
        call raw_part1%kill
        call raw_part2%kill
    else
        write(logfhandle,'(a)') '>>> STAGE 7 SKIPPED: an earlier stage failed'
    endif

    ! ============ STAGE 8: deapodization, the honest stage ============
    ! Run last because it flips deapod state and rebuilds the selection.
    ! forward_plane(x/env) = FT[env . (x/env)] = FT[x], i.e. the true
    ! envelope-free central section -- no inverse crime.
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 8: deapodization against envelope-free observations'
        env    = pcgop%get_env()
        invenv = pcgop%get_invenv()
        c         = BOX/2 + 1
        env_ctr   = env(c,c,c)
        env_edge  = env(1,c,c)
        recip_err = maxval(abs(env*invenv - 1.0))
        write(logfhandle,'(a,f10.6)')  '    envelope at centre       = ', env_ctr
        write(logfhandle,'(a,f10.6)')  '    envelope at edge (1,c,c) = ', env_edge
        write(logfhandle,'(a,es14.6)') '    max|env*invenv - 1|      = ', recip_err
        if( abs(env_ctr - 1.0) > 1.0e-5 )then
            write(logfhandle,'(a)') '    FAIL: envelope is not unity at the box centre'
            all_ok = .false.
        else if( env_edge >= env_ctr )then
            write(logfhandle,'(a)') '    FAIL: envelope does not decay away from the centre'
            all_ok = .false.
        else if( recip_err > 1.0e-4 )then
            write(logfhandle,'(a)') '    FAIL: invenv is not the reciprocal of env'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: envelope measured, normalized and invertible'
        endif
    endif
    if( all_ok )then
        allocate(xdiv(BOX,BOX,BOX))
        xdiv = phantom * invenv
        call projdirs%new(NPROJS, .false.)
        call projdirs%spiral
        call pcgop%set_volume(xdiv)
        do i = 1, NPROJS
            call projdirs%get_ori(i, e)
            call pcgop%forward_plane(e, y_planes(:,:,i))
        end do
        call pcgop%set_deapod(.true.)
        call pcgop%prep_particles(projdirs)
        allocate(recon_on(BOX,BOX,BOX), source=0.0)
        call pcgop%solve(y_planes, recon_on, maxits=40, rtol=1.0e-3, niters=niters)
        corr_on = corr_of(recon_on, phantom)
        write(logfhandle,'(a,i0,a,f9.5)') '    deapod ON : ', niters, ' iters, corr = ', corr_on
        call pcgop%set_deapod(.false.)
        call pcgop%prep_particles(projdirs)
        allocate(recon_off(BOX,BOX,BOX), source=0.0)
        call pcgop%solve(y_planes, recon_off, maxits=40, rtol=1.0e-3, niters=niters)
        corr_off = corr_of(recon_off, phantom)
        write(logfhandle,'(a,i0,a,f9.5)') '    deapod OFF: ', niters, ' iters, corr = ', corr_off
        if( corr_on <= corr_off )then
            write(logfhandle,'(a)') '    FAIL: deapodization did not improve recovery of the true volume'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: deapodization recovers the true volume better'
        endif
    else
        write(logfhandle,'(a)') '>>> STAGE 8 SKIPPED: an earlier stage failed'
    endif

    ! ============ STAGE 9: symmetry by coordinate replication ============
    ! Coordinate replication over a point group must produce the SAME normal
    ! system -- kernel Khat AND RHS b -- as reconstructing the symmetry-EXPANDED
    ! particle set at c1: the M symmetry mates are one measurement written M
    ! times (note section 2). Compared at the operator level, not through a
    ! solve: the kernelized operator is only approximately SPD, so a
    ! zero-tolerance solve driven to convergence loses positive-definiteness
    ! (that approximation is stages 5-7's business, not this stage's). c2 is
    ! lattice-exact, so the two builds agree to accumulator round-off plus the
    ! euler round-trip used to compose the expanded orientation set.
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 9: symmetry replication == particle-set expansion (c2, kernel)'
        call c1sym%new('c1')
        call c2sym%new('c2')
        nsym = c2sym%get_nsym()
        call pcgop%set_op_mode(PCG_OP_KERNEL)
        call pcgop%set_deapod(.true.)   ! stage 8 leaves it OFF; test the shipped config
        ! fresh, well-conditioned data: spiral projections of the phantom, T = 1
        call projdirs%new(NPROJS, .false.)
        call projdirs%spiral
        if( allocated(y_planes) ) deallocate(y_planes)
        allocate(y_planes(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), NPROJS))
        call pcgop%set_volume(phantom)
        do i = 1, NPROJS
            call projdirs%get_ori(i, e)
            call pcgop%forward_plane(e, y_planes(:,:,i))
        end do
        ! Path A: in-operator replication over c2 on NPROJS particles
        call pcgop%set_sym(c2sym)
        call pcgop%prep_particles(projdirs)
        call pcgop%build_operators(.true.)
        khat_a = pcgop%apply_normal_kernel(p_probe)
        b_mono = pcgop%apply_adjoint_all(y_planes)
        ! Path B: explicit c1 build of the c2-expanded particle set. The mate
        ! orientation is composed through sym%apply -- production's own symmetry
        ! adaptor, the one commander_reconstruct3D reaches through insert_fplane
        ! -- NOT by repeating the reconstructor's internal matmul(R_i, S_g). That
        ! is the point: composing both paths the same way would make this stage
        ! agree under a transposed or reversed composition too, and prove
        ! nothing. Going through sym%apply makes the stage an assertion that the
        ! operator's convention IS production's.
        call projdirs_exp%new(NPROJS*nsym, .false.)
        allocate(y_exp(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), NPROJS*nsym))
        c = 0
        do i = 1, NPROJS
            call projdirs%get_ori(i, e)
            do g = 1, nsym
                c = c + 1
                call c2sym%apply(e, g, e_exp)
                call projdirs_exp%set_ori(c, e_exp)
                y_exp(:,:,c) = y_planes(:,:,i)
            end do
        end do
        call pcgop%set_sym(c1sym)          ! back to no replication
        call pcgop%prep_particles(projdirs_exp)
        call pcgop%build_operators(.true.)
        khat_b = pcgop%apply_normal_kernel(p_probe)
        b_str  = pcgop%apply_adjoint_all(y_exp)
        kdiff      = sqrt(sum((khat_b-khat_a)**2)) / max(1.0, sqrt(sum(khat_a*khat_a)))
        stream_err = sqrt(sum((b_str -b_mono)**2)) / max(1.0, sqrt(sum(b_mono*b_mono)))
        write(logfhandle,'(a,es14.6,a,es14.6)') '    rel_err(Khat) = ', kdiff, '   rel_err(b) = ', stream_err
        if( kdiff > 5.0e-3 .or. stream_err > 5.0e-3 )then
            write(logfhandle,'(a)') '    FAIL: symmetry replication does not match the expanded-set build'
            all_ok = .false.
        else if( sqrt(sum(khat_a*khat_a)) <= TINY )then
            write(logfhandle,'(a)') '    FAIL: replicated kernel is trivially zero'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: symmetry replication equals the expanded-set build'
        endif
    else
        write(logfhandle,'(a)') '>>> STAGE 9 SKIPPED: an earlier stage failed'
    endif

    ! ============ STAGE 10: lambda scaling with effective data mass ============
    ! If every weighted observation is duplicated, B and H_data both double.
    ! A fixed absolute lambda would then become half as strong. Relative lambda
    ! must double with D so the complete normal system is scaled uniformly and
    ! a fixed-iteration PCG solve remains invariant.
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 10: lambda scaling with effective data mass'
        call pcgop%kill
        call pcgop%new(BOX, SMPD, 0.0)
        call pcgop%set_deapod(.false.)
        call pcgop%set_lambda_relative(MASS_LAMBDA_REL)
        call projdirs%new(NPROJS, .false.)
        call projdirs%spiral
        call pcgop%set_volume(phantom)
        do i = 1, NPROJS
            call projdirs%get_ori(i, e)
            call pcgop%forward_plane(e, y_planes(:,:,i))
        end do
        call pcgop%prep_particles(projdirs)
        call pcgop%begin_accum
        call pcgop%accumulate_batch(y_planes, NPROJS, 1)
        call pcgop%end_accum(.true.)
        call pcgop%set_op_mode(PCG_OP_KERNEL)
        data_scale = pcgop%get_data_scale()
        lambda_eff = pcgop%get_effective_lambda()
        allocate(recon_mass(BOX,BOX,BOX), source=0.0)
        call pcgop%solve_accum(recon_mass, maxits=2, rtol=0.0)

        call projdirs_exp%kill
        call projdirs_exp%new(2*NPROJS, .false.)
        do i = 1, NPROJS
            call projdirs%get_ori(i, e)
            call projdirs_exp%set_ori(i, e)
            call projdirs_exp%set_ori(NPROJS+i, e)
        end do
        if( allocated(y_exp) ) deallocate(y_exp)
        allocate(y_exp(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), 2*NPROJS))
        y_exp(:,:,1:NPROJS)          = y_planes
        y_exp(:,:,NPROJS+1:2*NPROJS) = y_planes
        ! Reach the doubled-data solve through the production distributed
        ! boundary: two raw worker artifacts, fixed-order master reduction,
        ! then relative-lambda derivation from the reduced D.
        raw_part1 = 'test_pcg_mass_part1.dat'
        raw_part2 = 'test_pcg_mass_part2.dat'
        call pcgop%kill
        call pcgop%new(BOX, SMPD, 0.0)
        call pcgop%set_deapod(.false.)
        call pcgop%prep_particles(projdirs_exp)
        call pcgop%begin_accum
        call pcgop%accumulate_batch(y_exp(:,:,1:NPROJS), NPROJS, 1)
        call pcgop%write_raw_accum(raw_part1, 1, 0, 1, 2, NPROJS, 'pcg_lambda_mass_v1')
        call pcgop%begin_accum
        call pcgop%accumulate_batch(y_exp(:,:,NPROJS+1:2*NPROJS), NPROJS, NPROJS+1)
        call pcgop%write_raw_accum(raw_part2, 1, 0, 2, 2, NPROJS, 'pcg_lambda_mass_v1')
        call pcgop%kill
        call pcg_reduce%new(BOX, SMPD, 0.0)
        call pcg_reduce%set_deapod(.false.)
        call pcg_reduce%set_lambda_relative(MASS_LAMBDA_REL)
        call pcg_reduce%begin_reduction
        call pcg_reduce%add_raw_accum(raw_part1, 1, 0, 1, 2, 'pcg_lambda_mass_v1', nraw)
        call pcg_reduce%add_raw_accum(raw_part2, 1, 0, 2, 2, 'pcg_lambda_mass_v1', nraw)
        call pcg_reduce%end_accum(.true.)
        call pcg_reduce%set_op_mode(PCG_OP_KERNEL)
        data_scale_dup = pcg_reduce%get_data_scale()
        lambda_eff_dup = pcg_reduce%get_effective_lambda()
        allocate(recon_mass_dup(BOX,BOX,BOX), source=0.0)
        call pcg_reduce%solve_accum(recon_mass_dup, maxits=2, rtol=0.0)

        mass_scale_err = abs(data_scale_dup - 2.0*data_scale) / max(TINY, 2.0*data_scale)
        mass_lambda_err = abs(lambda_eff_dup - 2.0*lambda_eff) / max(TINY, 2.0*lambda_eff)
        mass_solution_err = sqrt(sum((recon_mass_dup-recon_mass)**2)) / &
            &max(1.0, sqrt(sum(recon_mass*recon_mass)))
        write(logfhandle,'(a,es14.6,a,es14.6)') '    data scale: original = ', data_scale, &
            &' duplicated = ', data_scale_dup
        write(logfhandle,'(a,es14.6,a,es14.6)') '    lambda_eff: original = ', lambda_eff, &
            &' duplicated = ', lambda_eff_dup
        write(logfhandle,'(a,es14.6,a,es14.6,a,es14.6)') '    rel_err(scale x2) = ', mass_scale_err, &
            &' rel_err(lambda x2) = ', mass_lambda_err, ' rel_err(solution) = ', mass_solution_err
        if( mass_scale_err > MASS_SCALE_RELTOL .or. mass_lambda_err > MASS_SCALE_RELTOL )then
            write(logfhandle,'(a)') '    FAIL: effective lambda does not scale with duplicated data'
            all_ok = .false.
        else if( mass_solution_err > MASS_SOLVE_RELTOL )then
            write(logfhandle,'(a)') '    FAIL: relative lambda changes the solution under data duplication'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: relative lambda preserves the solution across data mass'
        endif
        call del_file(raw_part1)
        call del_file(raw_part2)
        call raw_part1%kill
        call raw_part2%kill
    else
        write(logfhandle,'(a)') '>>> STAGE 10 SKIPPED: an earlier stage failed'
    endif

    ! ============ STAGE 11: cropped accumulator and kernel invariants ============
    ! A cropped reconstruction covers the same physical extent with fewer
    ! Fourier samples. Native and cropped operators must therefore deposit the
    ! same raw sufficient statistics inside the common band, provided the crop
    ! boundary and the KB stencil are excluded. The cropped Khat must also
    ! retain the matrix-free equivalence already required at the native box.
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 11: full vs cropped raw accumulators and cropped kernel'
        call pcgop%kill
        call pcgop%new(BOX, SMPD, 0.0)
        call pcgop%set_deapod(.false.)
        call projdirs%kill
        call projdirs%new(NPROJS, .false.)
        call projdirs%spiral
        call projdirs_crop%new(NPROJS, .false.)
        crop_factor = real(CROP_BOX) / real(BOX)
        do i = 1, NPROJS
            call projdirs%get_ori(i, e)
            g = mod(i-1, NCTF) + 1
            ctfparms%smpd = SMPD; ctfparms%kv = KV; ctfparms%cs = CS; ctfparms%fraca = FRACA
            ctfparms%dfx = DFX_VALS(g); ctfparms%dfy = DFX_VALS(g)+ASTIG_VALS(g)
            ctfparms%angast = ANGAST_VALS(g); ctfparms%phshift = 0.
            shift = [1.0+0.3*real(mod(i-1,5)), -1.0+0.2*real(mod(i-1,7))]
            call e%set_ctfvars(ctfparms)
            call e%set_shift(shift)
            call projdirs%set_ori(i, e)
            ctfparms%smpd = CROP_SMPD
            call e%set_ctfvars(ctfparms)
            call e%set_shift(shift*crop_factor)
            call projdirs_crop%set_ori(i, e)
        enddo
        call pcgop%set_volume(phantom)
        do i = 1, NPROJS
            call projdirs%get_ori(i, e)
            call pcgop%forward_plane(e, gx_plane)
            Ti = pcgop%build_transfer(e%get_ctfvars(), e%get_2Dshift(), sig2arr)
            y_planes(:,:,i) = Ti * gx_plane
        enddo
        call pcgop%prep_particles(projdirs, use_ctf=.true., sig2=sig2_2d)
        call pcgop%begin_accum
        call pcgop%accumulate_batch(y_planes, NPROJS, 1)
        call pcgop%get_raw_accum(braw_full, draw_full)

        call pcg_crop%new(CROP_BOX, CROP_SMPD, 0.0)
        call pcg_crop%set_deapod(.true.)
        lims2_crop = pcg_crop%get_lims2()
        lims3_crop = pcg_crop%get_lims3()
        allocate(sig2_crop(0:CROP_BOX/2,NPROJS))
        do i = 1, NPROJS
            sig2_crop(:,i) = sig2arr(0:CROP_BOX/2)
        enddo
        allocate(y_crop(lims2_crop(1,1):lims2_crop(1,2), &
                       &lims2_crop(2,1):lims2_crop(2,2), NPROJS))
        y_crop = y_planes(lims2_crop(1,1):lims2_crop(1,2), &
            &lims2_crop(2,1):lims2_crop(2,2),:)
        call pcg_crop%prep_particles(projdirs_crop, use_ctf=.true., sig2=sig2_crop)
        call pcg_crop%begin_accum
        call pcg_crop%accumulate_batch(y_crop, NPROJS, 1)
        call pcg_crop%get_raw_accum(braw_crop, draw_crop)
        raw_ml = 'test_pcg_ml_prior.dat'
        call pcg_crop%write_raw_accum(raw_ml, 1, 0, 1, 1, NPROJS, 'pcg_ml_prior_v1')

        raw_lim    = CROP_BOX - 5
        crop_b_num = 0.0_dp
        crop_b_den = 0.0_dp
        crop_d_num = 0.0_dp
        crop_d_den = 0.0_dp
        do mraw = lims3_crop(3,1), lims3_crop(3,2)
            do kraw = lims3_crop(2,1), lims3_crop(2,2)
                do hraw = lims3_crop(1,1), lims3_crop(1,2)
                    if( hraw*hraw + kraw*kraw + mraw*mraw > raw_lim*raw_lim ) cycle
                    crop_b_num = crop_b_num + real(abs(braw_crop(hraw,kraw,mraw) - &
                        &braw_full(hraw,kraw,mraw))**2,dp)
                    crop_b_den = crop_b_den + real(abs(braw_full(hraw,kraw,mraw))**2,dp)
                    crop_d_num = crop_d_num + real((draw_crop(hraw,kraw,mraw) - &
                        &draw_full(hraw,kraw,mraw))**2,dp)
                    crop_d_den = crop_d_den + real(draw_full(hraw,kraw,mraw)**2,dp)
                enddo
            enddo
        enddo
        crop_b_err = real(sqrt(crop_b_num / max(1.0_dp,crop_b_den)))
        crop_d_err = real(sqrt(crop_d_num / max(1.0_dp,crop_d_den)))
        write(logfhandle,'(a,es14.6,a,es14.6)') '    common-band rel_err(B) = ', crop_b_err, &
            &' rel_err(D) = ', crop_d_err
        if( crop_b_err > CROP_RAW_RELTOL .or. crop_d_err > CROP_RAW_RELTOL )then
            write(logfhandle,'(a)') '    FAIL: cropped raw accumulators disagree with the full-box common band'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: cropped raw accumulators preserve the full-box common band'
        endif

        call pcg_crop%end_accum(.true.)
        allocate(p_crop(CROP_BOX,CROP_BOX,CROP_BOX))
        call random_number(p_crop)
        p_crop = p_crop - 0.5
        hm_crop = pcg_crop%apply_normal_matrixfree(p_crop)
        hk_crop = pcg_crop%apply_normal_kernel(p_crop)
        margin = 4
        lo = margin + 1
        hi = CROP_BOX - margin
        crop_kernel_err = sqrt(sum((hk_crop(lo:hi,lo:hi,lo:hi) - &
            &hm_crop(lo:hi,lo:hi,lo:hi))**2)) / &
            &max(1.0,sqrt(sum(hm_crop(lo:hi,lo:hi,lo:hi)**2)))
        write(logfhandle,'(a,es14.6)') '    cropped kernel interior rel_err = ', crop_kernel_err
        if( crop_kernel_err > EPS_INTERIOR )then
            write(logfhandle,'(a)') '    FAIL: cropped kernel disagrees with the matrix-free reference'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: cropped kernel agrees with the matrix-free reference'
        endif
    else
        write(logfhandle,'(a)') '>>> STAGE 11 SKIPPED: an earlier stage failed'
    endif

    ! ============ STAGE 12: FSC/SSNR ML prior operator ============
    ! Replay the cropped raw artifact so the prior is derived only after D is
    ! available, exactly as the production two-map path will do after the base
    ! independent-half FSC has been measured.
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 12: FSC/SSNR ML prior positivity and operator parity'
        call pcg_ml%new(CROP_BOX, CROP_SMPD, 0.0)
        call pcg_ml%set_deapod(.true.)
        call pcg_ml%prep_particles(projdirs_crop, use_ctf=.true., sig2=sig2_crop)
        call pcg_ml%begin_reduction
        call pcg_ml%add_raw_accum(raw_ml, 1, 0, 1, 1, 'pcg_ml_prior_v1', nraw)
        allocate(fsc_prior(CROP_BOX/2), source=0.5)
        do i = 1, size(fsc_prior)
            fsc_prior(i) = max(0.05, 0.9 - 0.1*real(i-1))
        enddo
        call pcg_ml%set_ml_prior(fsc_prior, 1.0, 100.0)
        call pcg_ml%end_accum(.true.)
        call pcg_ml%get_ml_prior(ml_prior_diag)
        call pcg_ml%get_ml_prior_stats(ml_prior_npositive, ml_prior_positive_min, ml_prior_positive_max, &
            &ml_prior_to_khat_l1, ml_prior_to_khat_rms)
        ml_prior_min = minval(ml_prior_diag)
        ml_prior_max = maxval(ml_prior_diag)
        write(logfhandle,'(a,es14.6,a,es14.6)') '    prior diagonal min = ', ml_prior_min, &
            &' max = ', ml_prior_max
        if( ml_prior_min < 0.0 .or. ml_prior_max <= 0.0 )then
            write(logfhandle,'(a)') '    FAIL: ML prior is not a nonzero positive-semidefinite diagonal'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: ML prior is a nonzero positive-semidefinite diagonal'
        endif
        write(logfhandle,'(a,i0,4(a,es14.6))') '    positive bins = ', ml_prior_npositive, &
            &' min+ = ', ml_prior_positive_min, ' max+ = ', ml_prior_positive_max, &
            &' prior/Khat L1 = ', ml_prior_to_khat_l1, ' RMS = ', ml_prior_to_khat_rms
        if( ml_prior_npositive /= count(ml_prior_diag > 0.0) .or. &
            &abs(ml_prior_positive_min-minval(ml_prior_diag, mask=ml_prior_diag > 0.0)) > 0.0 .or. &
            &abs(ml_prior_positive_max-ml_prior_max) > 0.0 .or. &
            &ml_prior_to_khat_l1 <= 0.0 .or. ml_prior_to_khat_rms <= 0.0 )then
            write(logfhandle,'(a)') '    FAIL: ML prior summary diagnostics are inconsistent'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: ML prior summary diagnostics are consistent'
        endif
        hm_ml = pcg_ml%apply_normal_matrixfree(p_crop)
        hk_ml = pcg_ml%apply_normal_kernel(p_crop)
        margin = 4
        lo = margin + 1
        hi = CROP_BOX - margin
        ml_kernel_err = sqrt(sum((hk_ml(lo:hi,lo:hi,lo:hi) - &
            &hm_ml(lo:hi,lo:hi,lo:hi))**2)) / &
            &max(1.0,sqrt(sum(hm_ml(lo:hi,lo:hi,lo:hi)**2)))
        ml_prior_energy = sum(real(p_crop,dp) * real(hk_ml-hk_crop,dp))
        write(logfhandle,'(a,es14.6,a,es14.6)') '    ML kernel interior rel_err = ', ml_kernel_err, &
            &' dot(p,P_tau*p) = ', real(ml_prior_energy)
        if( ml_kernel_err > EPS_INTERIOR )then
            write(logfhandle,'(a)') '    FAIL: ML kernel disagrees with the matrix-free reference'
            all_ok = .false.
        else if( ml_prior_energy <= 0.0_dp )then
            write(logfhandle,'(a)') '    FAIL: ML prior does not add positive quadratic energy'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: ML prior preserves operator parity and adds positive energy'
        endif
        call del_file(raw_ml)
        call raw_ml%kill
    else
        write(logfhandle,'(a)') '>>> STAGE 12 SKIPPED: an earlier stage failed'
    endif

    call pcgop%kill
    call pcg_reduce%kill
    call pcg_crop%kill
    call pcg_ml%kill
    call projdirs%kill
    call projdirs_exp%kill
    call projdirs_crop%kill
    call e%kill
    call e_exp%kill
    call c1sym%kill
    call c2sym%kill
    if( all_ok )then
        call simple_end('**** SIMPLE_TEST_PCG_RECON NORMAL STOP ****')
    else
        THROW_HARD('TEST_PCG_RECON FAILED')
    endif

  contains

    real function corr_of( a, b )
        real, intent(in) :: a(:,:,:), b(:,:,:)
        real :: ma, mb, num, da, db
        ma  = sum(a)/real(size(a))
        mb  = sum(b)/real(size(b))
        num = sum((a-ma)*(b-mb))
        da  = sum((a-ma)**2)
        db  = sum((b-mb)**2)
        corr_of = num / sqrt(max(da*db, TINY))
    end function corr_of

end subroutine exec_test_pcg_recon

end module simple_commanders_test_highlevel
