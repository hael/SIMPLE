!@descr: for all highlevel tests
module simple_commanders_test_highlevel
use simple_commanders_api
use simple_stream_api
use simple_commanders_project_core, only: commander_new_project, commander_selection
use simple_commanders_project_mov,  only: commander_import_movies
use simple_commanders_atoms,        only: commander_pdb2mrc
use simple_commanders_reproject,    only: commander_reproject
use simple_commanders_pick,         only: commander_make_pickrefs, commander_pick
use simple_commanders_sim,          only: commander_simulate_particles, commander_simulate_movie
use simple_commanders_preprocess,   only: commander_ctf_estimate, commander_motion_correct, commander_preprocess
use simple_commanders_abinitio2D,   only: commander_abinitio2D
use simple_commanders_cluster2D,    only: commander_cluster2D
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

type, extends(commander_base) :: commander_test_movie_ppca_subproject_distr
    contains
        procedure :: execute      => exec_test_movie_ppca_subproject_distr
end type commander_test_movie_ppca_subproject_distr

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
    use simple_imghead,       only: find_ldim_nptcls
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
        call find_ldim_nptcls(vol_file, ldim, nptcls_stk, smpd_stk)
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
        call find_ldim_nptcls(outstk, ldim, nptcls_stk, smpd_stk)
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
    use simple_imghead,       only: find_ldim_nptcls
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
    ! ---- generate 6VXX volume from built-in molecule data ----
    write(logfhandle,'(a)') '>>> TEST_REPROJECT: generating 6VXX.mrc volume'
    mol = sars_cov2_spkgp_6vxx()
    call molecule%pdb2mrc(smpd=SMPD, mol=mol)
    call params%new(cline)
    all_ok = .true.
    ! ---- define expected output file names ----
    vol_file = '6VXX.mrc'
    outstk   = 'reprojs.mrcs'
    outori   = 'reproject_oris'//trim(TXT_EXT)
    ! ---- check volume was generated ----
    write(logfhandle,'(a)') '>>> CHECK: volume file exists'
    if( .not. file_exists(vol_file) )then
        write(logfhandle,'(a)') '    FAIL: '//vol_file%to_char()//' not found'
        THROW_HARD('TEST_REPROJECT FAILED: volume not generated')
    else
        call find_ldim_nptcls(vol_file, ldim, nptcls_stk, smpd_stk)
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
        call find_ldim_nptcls(outstk, ldim, nptcls_stk, smpd_stk)
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
    class(commander_test_simulated_workflow), intent(inout) :: self
    class(cmdline),                           intent(inout) :: cline
    character(len=*), parameter         :: filetab_file='filetab.txt'
    type(cmdline)                       :: cline_pdb2mrc, cline_projection, cline_make_pick_refs, cline_sim_mov, cline_ctf_est, cline_mot_corr
    type(cmdline)                       :: cline_new_project, cline_import_movies, cline_segpick, cline_refpick, cline_sim_ptcls, cline_abinitio2D
    type(cmdline)                       :: cline_cluster2D 
    type(parameters)                    :: params
    type(commander_new_project)         :: xnew_project
    type(commander_pdb2mrc)             :: xpdb2mrc
    type(commander_reproject)           :: xreproject
    type(commander_make_pickrefs)       :: xmakepickrefs
    type(commander_simulate_movie)      :: xsimov
    type(commander_motion_correct)      :: xmotcorr
    type(commander_ctf_estimate)        :: xctf_estimate
    type(commander_import_movies)       :: ximport_movies
    type(commander_pick)                :: xsegpick, xrefpick
    type(commander_simulate_particles)  :: xsim_ptcls      
    type(commander_abinitio2D)          :: xabinitio2D
    type(commander_cluster2D)           :: xcluster2D
    integer                             :: rc
    type(string)                        :: cmd, projfile
    logical                             :: mrc_exists

    ! Download pdb & generate volume and reprojections
    inquire(file="1JYX.mrc", exist=mrc_exists)
    if( .not. mrc_exists )then
       write(*, *) 'Downloading the example dataset...'
       cmd = 'curl -s -o 1JYX.pdb https://files.rcsb.org/download/1JYX.pdb'
       call execute_command_line(cmd%to_char(), exitstat=rc)
       print *, 'Converting .pdb to .mrc...'
       call cline_pdb2mrc%set('smpd',            1.)
       call cline_pdb2mrc%set('pdbfile', '1JYX.pdb')
       call cline_pdb2mrc%checkvar('smpd',        1)
       call cline_pdb2mrc%checkvar('pdbfile',     2)
       call cline_pdb2mrc%check()
       call xpdb2mrc%execute(cline_pdb2mrc)
       call cline_pdb2mrc%kill()
       cmd = 'rm 1JYX.pdb'
       call execute_command_line(cmd%to_char(), exitstat=rc)
    endif

    ! Generate reprojections
    print *, 'Projecting 1JYX.mrc...'
    call cline_projection%set('vol1', '1JYX.mrc')
    call cline_projection%set('smpd',        1.3)
    call cline_projection%set('pgrp',       'c1')
    call cline_projection%set('mskdiam',    180.)
    call cline_projection%set('nspace',      10.)
    call cline_projection%set('nthr',        16.)
    call xreproject%execute(cline_projection)
    call cline_projection%kill()

    ! Make pick references
    call cline_make_pick_refs%set('pickrefs', 'reprojs.mrcs')
    call cline_make_pick_refs%set('nthr',                16.)
    call xmakepickrefs%execute(cline_make_pick_refs)
    call cline_make_pick_refs%kill()

    ! Simulate movie
    call cline_sim_mov%set('stk', 'reprojs.mrcs')
    call cline_sim_mov%set('xdim',          4096)
    call cline_sim_mov%set('ydim',          4096)
    call cline_sim_mov%set('nthr',           16.)
    call xsimov%execute(cline_sim_mov)
    call cline_sim_mov%kill()

    ! Project creation
    call cline_new_project%set('projname',  params%projname)
    call xnew_project%execute(cline_new_project)
    call cline_new_project%kill()
    projfile = params%projname//'.simple'

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

    ! Motion correction - algorithm iso
    call cline_mot_corr%set('stk', 'simulated_movies.mrcs')
    call cline_mot_corr%set('algorithm',    'iso')
    call cline_mot_corr%set('nparts',           4)
    call cline_mot_corr%set('nthr',            16)
    call xmotcorr%execute(cline_mot_corr)
    call cline_mot_corr%kill()
    ! Motion correction - algorithm patch
    call cline_mot_corr%set('stk',      'simulated_movies.mrcs')
    call cline_mot_corr%set('algorithm',                'patch')
    call cline_mot_corr%set('nparts',                         4)
    call cline_mot_corr%set('nthr',                          16)
    call xmotcorr%execute(cline_mot_corr)
    call cline_mot_corr%kill()

    ! CTF estimate - patch yes
    call cline_ctf_est%set('ctfpatch', 'yes')
    call cline_ctf_est%set('nparts',       4)
    call cline_ctf_est%set('nthr',        16)
    call xctf_estimate%execute(cline_ctf_est)
    call cline_ctf_est%kill()
    ! CTF estimate - patch yes
    call cline_ctf_est%set('prg',     'ctf_estimate')
    call cline_ctf_est%set('ctfpatch',          'no')
    call cline_ctf_est%set('nparts',               4)
    call cline_ctf_est%set('nthr',                16)
    call xctf_estimate%execute(cline_ctf_est)
    call cline_ctf_est%kill()

    ! Segmentation-based picking
    call cline_segpick%set('prg',       'pick')
    call cline_segpick%set('picker', 'segdiam')
    call cline_segpick%set('nparts',         4)
    call cline_segpick%set('nthr',          16)      
    call xsegpick%execute(cline_segpick)
    call cline_segpick%kill()  

    ! Reference-based picking
     call cline_refpick%set('prg',               'pick')
     call cline_refpick%set('pickrefs', 'pickrefs.mrcs')
     call cline_refpick%set('nparts',                 4)
     call cline_refpick%set('nthr',                  16)      
     call xrefpick%execute(cline_refpick)
     call cline_refpick%kill()  

    ! 2D analysis
    ! generate simulate particles
    call cline_sim_ptcls%set('prg', 'simulate_particles')
    call cline_sim_ptcls%set('vol1',          '1JXY.mrc')
    call cline_sim_ptcls%set('ctf',                'yes')
    call cline_sim_ptcls%set('nptcls',               500)
    call cline_sim_ptcls%set('smpd',                 1.3)
    call cline_sim_ptcls%set('snr',                  0.5)
    call cline_sim_ptcls%set('pgrp',                'c1')
    call cline_sim_ptcls%set('mskdiam',              180)
    call cline_sim_ptcls%set('nthr',                  16)
    call xsim_ptcls%execute(cline_sim_ptcls)
    call cline_sim_ptcls%kill()
    ! abinitio2D
     call cline_abinitio2D%set('prg', 'abinitio2D')
     call cline_abinitio2D%set('mskdiam',      180)
     call cline_abinitio2D%set('ncls',         100)
     call cline_abinitio2D%set('nthr',          16)
     call xabinitio2D%execute(cline_abinitio2D)
     call cline_abinitio2D%kill()

    ! cluster2D
    call cline_cluster2D%set('prg', 'cluster2D')
    call cline_cluster2D%set('mskdiam',     180)
    call cline_cluster2D%set('ncls',        100)
    call cline_cluster2D%set('nthr',         16)
    call xcluster2D%execute(cline_abinitio2D)
    call cline_cluster2D%kill()
    call simple_end('**** SIMPLE_TEST_SIMULATED_WORKFLOW_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_simulated_workflow

!>  \brief  Integration test: split project into subprojects, run in parallel, merge back.
!  This test exercises the new generate_scripts_subprojects / gen_subproject_scripts_and_schedule
!  machinery end-to-end.
!    1. Simulates particles and runs coarse abinitio2D to obtain class labels
!    2. Splits the project into per-class subprojects (extract_subproj pattern)
!    3. Runs a program on each subproject in parallel via generate_scripts_subprojects
!    4. Merges all subproject results with merge_chunk_projfiles
!    5. Validates that particle count is preserved
subroutine exec_test_subproject_distr( self, cline )
    use simple_atoms,                   only: atoms
    use simple_molecule_data,           only: molecule_data, sars_cov2_spkgp_6vxx
    use simple_imghead,                 only: find_ldim_nptcls
    use simple_commanders_project_ptcl, only: commander_import_particles
    class(commander_test_subproject_distr), intent(inout) :: self
    class(cmdline),                         intent(inout) :: cline
    integer,          parameter :: NPTCLS_SIM  = 500      ! particles to simulate
    integer,          parameter :: NCLS_COARSE = 4        ! coarse classes
    integer,          parameter :: MAXKEYS     = 20       ! chash capacity
    real,             parameter :: SMPD        = 1.3
    character(len=*), parameter :: PROJNAME = 'test_subproj_distr'
    type(parameters)                    :: params
    type(sp_project)                    :: spproj, spproj_sub, spproj_merged
    type(commander_simulate_particles)  :: xsim_ptcls
    type(commander_new_project)         :: xnew_project
    type(commander_import_particles)    :: ximport_particles
    type(commander_abinitio2D)          :: xabinitio2D
    type(cmdline)                       :: cline_sim, cline_new_proj, cline_import, cline_2D, cline_sub
    type(qsys_env)                      :: qenv
    type(chash),  allocatable           :: jobs_descr(:)
    type(string), allocatable           :: subproj_fnames(:), subproj_dirs(:)
    integer,      allocatable           :: class_labels(:)
    real,         allocatable           :: state_labels(:)
    type(string)                        :: projname_sub, projfile_sub, cwd_root, sim_stk, sim_oris
    type(molecule_data)                 :: mol
    type(atoms)                         :: molecule
    integer :: icls, nsub, isub, nptcls_sub, nptcls_merged, nptcls_orig
    logical :: all_ok
    call params%new(cline)
    call simple_getcwd(cwd_root)
    ! 1. Simulate particles and populate project
    write(logfhandle,'(a)') '>>> Step 1: simulate particles from 6VXX molecule data'
    mol = sars_cov2_spkgp_6vxx()
    call molecule%pdb2mrc(smpd=SMPD, mol=mol)
    call cline_sim%set('prg',      'simulate_particles')
    call cline_sim%set('vol1',           'molecule.mrc')
    call cline_sim%set('smpd',                     SMPD)
    call cline_sim%set('mskdiam',                   180)
    call cline_sim%set('nthr',                       16)
    call cline_sim%set('nptcls',             NPTCLS_SIM)
    call cline_sim%set('pgrp',                     'c1')
    call cline_sim%set('snr',                      0.01)
    call cline_sim%set('ctf',                     'yes')
    call params%new(cline_sim)
    call xsim_ptcls%execute(cline_sim)
    call cline_sim%kill()
    ! store paths to simulated output (relative to cwd_root)
    sim_stk  = 'simulated_particles.mrc'
    sim_oris = 'simulated_oris'//trim(TXT_EXT)
    ! 1b. Create project and import simulated particles
    write(logfhandle,'(a)') '>>> Step 1b: create project & import particles'
    call cline_new_proj%set('projname', PROJNAME)
    call xnew_project%execute(cline_new_proj)
    call cline_new_proj%kill()
    ! import_particles (new_project changed cwd into project dir)
    call cline_import%set('prg',             'import_particles')
    call cline_import%set('mkdir',                         'no')
    call cline_import%set('projfile',       PROJNAME//'.simple')
    call cline_import%set('stk',       '../'//sim_stk%to_char())
    call cline_import%set('deftab',   '../'//sim_oris%to_char())
    call cline_import%set('smpd',                          SMPD)
    call cline_import%set('kv',                           300.0)
    call cline_import%set('cs',                             2.7)
    call cline_import%set('fraca',                          0.1)
    call cline_import%set('ctf',                          'yes')
    call params%new(cline_import)
    call ximport_particles%execute(cline_import)
    call cline_import%kill()
    write(logfhandle,'(a)') '    project populated with '//int2str(NPTCLS_SIM)//' particles'
    ! update cwd_root to project directory (new_project changed cwd)
    call simple_getcwd(cwd_root)
    ! 2. Coarse abinitio2D to get class labels
    ! (cwd is now inside project dir after new_project)
    write(logfhandle,'(a)') '>>> Step 2: coarse abinitio2D'
    call cline_2D%set('prg',             'abinitio2D')
    call cline_2D%set('projfile', PROJNAME//'.simple')
    call cline_2D%set('ncls',             NCLS_COARSE)
    call cline_2D%set('mkdir',                   'no')
    call cline_2D%set('mskdiam',                  180)
    call cline_2D%set('smpd',                    SMPD)
    call cline_2D%set('nthr',                      16)
    call params%new(cline_2D)
    call xabinitio2D%execute(cline_2D)
    ! 3. Read project and split into subprojects by class 
    write(logfhandle,'(a)') '>>> Step 3: splitting into subprojects'
    call spproj%read(string(PROJNAME//'.simple'))
    nptcls_orig  = spproj%get_nptcls()
    class_labels = spproj%os_ptcl2D%get_all_asint('class')
    nsub         = NCLS_COARSE
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
        call jobs_descr(icls)%set('nthr',               int2str(params%nthr))
        call jobs_descr(icls)%set('mskdiam',        real2str(params%mskdiam))
        call jobs_descr(icls)%set('mkdir',                              'no')
        call jobs_descr(icls)%set('objfun',                             'cc')
        call spproj_sub%kill
    end do
    ! 4. Generate scripts and execute subprojects in parallel ----
    write(logfhandle,'(a)') '>>> Step 4: parallel execution of subprojects'
    call cline%set('projfile', PROJNAME//'.simple')
    call cline%set('mskdiam',                 180.)
    call cline%set('smpd',                    SMPD)
    call cline%set('nthr',                     16.)
    call cline%set('ncunits',                 nsub)
    call params%new(cline)
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

subroutine exec_test_movie_ppca_subproject_distr( self, cline )
    class(commander_test_movie_ppca_subproject_distr), intent(inout) :: self
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
    integer                     :: nmovies, nsub, isub, i, fromp, top, nchunk
    integer                     :: nptcls_in, ldim(3)
    ! if( .not. cline%defined('filetab') )then
    !     THROW_HARD('filetab is required; '//'test_movie_ppca_subproject_distr')
    ! endif
    ! if( .not. cline%defined('smpd') )then
    !     THROW_HARD('smpd is required; '//'test_movie_ppca_subproject_distr')
    ! endif
    !if( .not. cline%defined('nthr')    ) call cline%set('nthr', 1)
    call params%new(cline)
    call read_filetable(params%filetab, movie_fnames)
    nmovies = size(movie_fnames)
    if( nmovies < 1 )then
        THROW_HARD('empty filetab; test_movie_ppca_subproject_distr')
    endif
    nsub = min(nmovies, params%nparts)
    call simple_getcwd(cwd_root)
    allocate(jobs_descr(nsub), subproj_dirs(nsub), denoised_all(nsub))
    write(logfhandle,'(a)') '>>> Step 1: split movie list into equal chunks and generate chunk stacks'
    do isub = 1, nsub
        fromp = ((isub - 1) * nmovies) / nsub + 1
        top   = (isub * nmovies) / nsub
        nchunk = top - fromp + 1
        write(subid,'(I8.8)') isub
        subproj_dirs(isub) = cwd_root%to_char()//'/'//&
            &'movie_denoise_subproj_'//trim(adjustl(subid))
        call simple_mkdir(subproj_dirs(isub))
        allocate(chunk_fnames(nchunk))
        do i = 1, nchunk
            chunk_fnames(i) = movie_fnames(fromp + i - 1)
        enddo
        chunk_filetab = subproj_dirs(isub)%to_char()//'/'//'chunk_movies_'//trim(adjustl(subid))//'.txt'
        chunk_stk     = subproj_dirs(isub)%to_char()//'/'//'chunk_movies_'//trim(adjustl(subid))//'.mrcs'
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
        call jobs_descr(isub)%set('stk',      'chunk_movies_'//trim(adjustl(subid))//'.mrcs')
        call jobs_descr(isub)%set('outstk',                                     denoised_stk)
        call jobs_descr(isub)%set('smpd',                              real2str(params%smpd))
        call jobs_descr(isub)%set('nthr',                               int2str(params%nthr))
        call jobs_descr(isub)%set('neigs',                                             '160')
        call jobs_descr(isub)%set('pca_mode',                                        'mppca')
        call jobs_descr(isub)%set('mppca_k',                                             '4')
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
    call cline_stack%set('prg',     'stack')
    call cline_stack%set('mkdir',      'no')
    call cline_stack%set('filetab', 'ppca_denoised_chunks.txt')
    call cline_stack%set('outstk',  'ppca_denoised_all.mrcs')
    call cline_stack%set('smpd',    params%smpd)
    call xstack%execute(cline_stack)
    call cline_stack%kill()
    do isub = 1, nsub
        call jobs_descr(isub)%kill
    end do
    deallocate(jobs_descr, subproj_dirs, denoised_all)
    call qenv%kill
    call simple_end('**** TEST_MOVIE_PPCA_SUBPROJECT_DISTR NORMAL STOP ****')
end subroutine exec_test_movie_ppca_subproject_distr

end module simple_commanders_test_highlevel
