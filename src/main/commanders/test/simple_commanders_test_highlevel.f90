!@descr: for all highlevel tests
module simple_commanders_test_highlevel
use simple_commanders_api
use simple_stream_api
use simple_commanders_project_core, only: commander_new_project, commander_selection
use simple_commanders_project_mov,  only: commander_import_movies
use simple_commanders_reproject,    only: commander_reproject
use simple_commanders_pick,         only: commander_pick, commander_extract
use simple_commanders_sim,          only: commander_simulate_particles, commander_simulate_movie
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

type, extends(commander_base) :: commander_test_flex_preimage_identity
    contains
        procedure :: execute      => exec_test_flex_preimage_identity
end type commander_test_flex_preimage_identity

type, extends(commander_base) :: commander_test_flex_preimage_basis_ab
    contains
        procedure :: execute      => exec_test_flex_preimage_basis_ab
end type commander_test_flex_preimage_basis_ab

contains

subroutine exec_test_flex_preimage_identity( self, cline )
    use simple_builder,                 only: builder
    use simple_flex_diffmap_rec3D,      only: test_fake_preimage_against_reconstruct3D
    use simple_parameters,              only: parameters
    class(commander_test_flex_preimage_identity), intent(inout) :: self
    class(cmdline),                                  intent(inout) :: cline
    type(parameters) :: params
    type(builder)    :: build
    integer, allocatable :: pinds(:)
    integer :: nptcls
    if( .not.cline%defined('projfile') ) THROW_HARD('flex_preimage_identity requires projfile=flex_registered_particles.simple')
    if( .not.cline%defined('vol1') ) THROW_HARD('flex_preimage_identity requires vol1=<fixed mean volume>')
    if( .not.cline%defined('nspace') ) THROW_HARD('flex_preimage_identity requires nspace=<projection grid size>')
    if( .not.cline%defined('oritype') ) call cline%set('oritype','ptcl3D')
    if( .not.cline%defined('mkdir') ) call cline%set('mkdir','yes')
    call cline%set('ml_reg','no')
    call build%init_params_and_build_general_tbox(cline,params,do3d=.true.)
    call build%spproj_field%sample4rec([params%fromp,params%top],nptcls,pinds)
    if( nptcls<3 ) THROW_HARD('flex_preimage_identity requires at least three active ptcl3D particles')
    call test_fake_preimage_against_reconstruct3D(params,build,pinds)
    deallocate(pinds)
    call build%kill_general_tbox
    call simple_end('**** SIMPLE_TEST_FLEX_PREIMAGE_IDENTITY NORMAL STOP ****')
end subroutine exec_test_flex_preimage_identity

subroutine exec_test_flex_preimage_basis_ab( self, cline )
    use simple_builder,                only: builder
    use simple_flex_analysis_strategy, only: flex_analysis_strategy, create_flex_analysis_strategy, &
        &run_flex_preimage_basis_ab_test
    use simple_parameters,             only: parameters
    class(commander_test_flex_preimage_basis_ab), intent(inout) :: self
    class(cmdline),                                intent(inout) :: cline
    type(parameters) :: params
    type(builder) :: build
    class(flex_analysis_strategy), allocatable :: strategy
    if( .not.cline%defined('projfile') ) THROW_HARD('flex_preimage_basis_ab requires projfile=<input project>')
    if( .not.cline%defined('vol1') ) THROW_HARD('flex_preimage_basis_ab requires vol1=<fixed mean volume>')
    if( .not.cline%defined('nspace') ) THROW_HARD('flex_preimage_basis_ab requires nspace=<projection grid size>')
    if( .not.cline%defined('neigs') ) THROW_HARD('flex_preimage_basis_ab requires neigs=<diffusion rank>')
    if( .not.cline%defined('oritype') ) call cline%set('oritype','ptcl3D')
    if( .not.cline%defined('mkdir') ) call cline%set('mkdir','yes')
    call cline%set('nparts','1')
    strategy=create_flex_analysis_strategy(cline)
    call strategy%initialize(params,build,cline)
    call run_flex_preimage_basis_ab_test(params,build,cline)
    call strategy%finalize_run(params,build,cline)
    call strategy%cleanup(params)
    deallocate(strategy)
    call build%kill_general_tbox
    call simple_end('**** SIMPLE_TEST_FLEX_PREIMAGE_BASIS_AB NORMAL STOP ****')
end subroutine exec_test_flex_preimage_basis_ab

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
    integer,          parameter :: NPER_MOVIE     = 16
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
    type(string)                        :: system_name, test_workdir, vol_file, reproj_file
    type(string)                        :: movie_fname, subset_fname, optimal_fname, params_fname
    type(string)                        :: movie_files(NMOVIES)
    character(len=XLONGSTRLEN)          :: workflow_root_path
    integer                             :: i, j, iproj, proj_inds(NPER_MOVIE), ldim(3), nprojs_stk
    integer                             :: reproj_box, npickrefs, nptcls, ncls, status
    real                                :: pickref_smpd, pickref_width

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
        case('1jxy')
            vol_file    = '1JXY.mrc'
            reproj_file = 'reprojs_1JXY.mrcs'
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
    call cline_projection%set('pgrp',                             'c1')
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
        call cline_sim_mov%set('defocus',                    2.0)
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
    call cline_abinitio3D%set('pgrp',                       'c1')
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

end module simple_commanders_test_highlevel
