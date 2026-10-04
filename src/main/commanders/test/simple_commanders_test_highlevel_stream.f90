!@descr: end-to-end stream and particle-simulation tests of simple_commanders_test_highlevel
submodule (simple_commanders_test_highlevel) simple_commanders_test_highlevel_stream
implicit none
#include "simple_local_flags.inc"

contains

module subroutine exec_test_mini_stream_quantitative( self, cline )
    use simple_atoms,         only: atoms
    use simple_imghead,       only: find_ldim_nptcls
    use simple_molecule_data, only: molecule_data, betagal_1jyx, sars_cov2_spkgp_6vxx
    use simple_oris,          only: oris
    use simple_string_utils,  only: lowercase
    use simple_ui,            only: make_ui
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    class(commander_test_mini_stream), intent(inout) :: self
    class(cmdline),                    intent(inout) :: cline
    character(len=*), parameter :: MOVIE_FILE     = 'simulate_movie.mrc'
    character(len=*), parameter :: OPTIMAL_FILE   = 'optimal_movie_average.mrc'
    character(len=*), parameter :: PARAMS_FILE    = 'simulate_movie_params.txt'
    character(len=*), parameter :: FILETAB_FILE   = 'mini_stream_micrographs.txt'
    character(len=*), parameter :: PROJFILE       = 'mini_stream.simple'
    real,             parameter :: SMPD            = 1.3
    real,             parameter :: CS              = 2.7
    real,             parameter :: KV              = 300.0
    real,             parameter :: FRACA           = 0.1
    real,             parameter :: DEFOCUS_TOL     = 0.10
    ! First deterministic 6VXX run recovered 22/48 = 0.458; retain a
    ! 0.058 absolute margin while still detecting a substantial regression.
    real,             parameter :: PICK_RECALL_MIN = 0.40
    real,             parameter :: PICK_RATIO_MAX  = 1.50
    real,             parameter :: MAX_AREA_FRACTION = 0.30
    real,             parameter :: MSKDIAM         = 180.0
    integer,          parameter :: PARTICLE_BOX    = 192
    integer,          parameter :: MICROGRAPH_BOX_6VXX = 1024
    integer,          parameter :: MICROGRAPH_BOX_1JXY = 1280
    integer,          parameter :: NPARTICLES      = 8
    integer,          parameter :: NMICROGRAPHS    = 6
    integer,          parameter :: NFRAMES         = 8
    type(commander_reproject)      :: xreproject
    type(commander_simulate_movie) :: xsim_movie
    type(commander_mini_stream)    :: xmini_stream
    type(cmdline)                  :: cline_reproj, cline_sim, cline_mini
    type(atoms)                    :: molecule
    type(molecule_data)            :: mol
    type(image)                    :: cavg
    type(oris)                     :: truth
    type(sp_project)               :: result
    type(string)                   :: cwd_saved, fixture_root, mini_dir
    type(string)                   :: particle_stack_path, filetab_path, project_path
    type(string)                   :: movie_name, params_name, value, cavg_stack
    type(string)                   :: requested_suite, system_name, vol_file, pgrp
    type(string)                   :: movie_paths(NMICROGRAPHS), params_paths(NMICROGRAPHS)
    integer, allocatable           :: classes(:), populations(:), shape_ranks(:)
    integer                       :: i, imic, isuite, rank, status, ldim(3), nimages, micrograph_box
    integer                       :: npicked, nclasses, ncavgs, nranked, expected_particles
    real                          :: dfx, dfy, truth_dfx, truth_dfy, defocus_error
    real                          :: pick_ratio, cavg_smpd, cavg_variance, area_fraction

    call make_ui
    call set_fixed_seed(20260925)
    call simple_getcwd(cwd_saved)
    requested_suite = ''
    if( cline%defined('suite') )then
        requested_suite = cline%get_carg('suite')
        requested_suite = lowercase(requested_suite%to_char())
    endif
    if( requested_suite == 'list' )then
        write(logfhandle,'(a)') 'Available suites for mini_stream:'
        write(logfhandle,'(a)') '  6vxx'
        write(logfhandle,'(a)') '  1jxy'
        return
    endif
    if( requested_suite%strlen_trim() > 0 .and. requested_suite /= '6vxx' .and. requested_suite /= '1jxy' )&
        &THROW_HARD('no sub-suite '//requested_suite%to_char()//' in mini_stream; use suite=list')

    do isuite = 1, 2
        if( isuite == 1 )then
            system_name = '6vxx'
            vol_file    = '6VXX.mrc'
            pgrp        = 'c3'
            micrograph_box = MICROGRAPH_BOX_6VXX
        else
            system_name = '1jxy'
            vol_file    = '1JXY.mrc'
            pgrp        = 'c1'
            micrograph_box = MICROGRAPH_BOX_1JXY
        endif
        if( requested_suite%strlen_trim() > 0 .and. requested_suite /= system_name ) cycle
        write(logfhandle,'(a)') '---- TEST SUITE: mini_stream '//system_name%to_char()//' ----'
        call set_fixed_seed(20260925 + isuite)
        fixture_root = filepath(cwd_saved, 'test_mini_stream_'//system_name%to_char()//'_'//int2str(get_process_id()))
    if( dir_exists(fixture_root) ) THROW_HARD('TEST_MINI_STREAM FAILED: fixture directory already exists')
    call simple_mkdir(fixture_root)
    call simple_chdir(fixture_root, status)
    if( status /= 0 ) THROW_HARD('TEST_MINI_STREAM FAILED: could not enter fixture directory')

    select case(system_name%to_char())
        case('6vxx')
            mol = sars_cov2_spkgp_6vxx()
        case('1jxy')
            ! SIMPLE's embedded provider follows the underlying 1JYX PDB identifier.
            mol = betagal_1jyx()
    end select
    call molecule%pdb2mrc(volfile=vol_file, smpd=SMPD, mol=mol, center_pdb=.true., &
        &vol_dim=[PARTICLE_BOX, PARTICLE_BOX, PARTICLE_BOX])
    call molecule%kill()
    call cline_reproj%set('prg',      'reproject')
    call cline_reproj%set('mkdir',           'no')
    call cline_reproj%set('vol1',        vol_file)
    call cline_reproj%set('outstk', 'mini_stream_particles.mrcs')
    call cline_reproj%set('smpd',            SMPD)
    call cline_reproj%set('pgrp',            pgrp)
    call cline_reproj%set('mskdiam',      MSKDIAM)
    call cline_reproj%set('nspace',    NPARTICLES)
    call cline_reproj%set('nthr',               1)
    call xreproject%execute(cline_reproj)
    call cline_reproj%kill()
    particle_stack_path = simple_abspath(string('mini_stream_particles.mrcs'))
    if( .not. file_exists(particle_stack_path) )&
        &THROW_HARD('TEST_MINI_STREAM FAILED: molecular reprojections were not created')
    call find_ldim_nptcls(particle_stack_path, ldim, nimages)
    if( nimages /= NPARTICLES )&
        &THROW_HARD('TEST_MINI_STREAM FAILED: molecular reprojection count is incorrect')
    if( any(ldim(1:2) /= [PARTICLE_BOX, PARTICLE_BOX]) )&
        &THROW_HARD('TEST_MINI_STREAM FAILED: molecular reprojection dimensions are incorrect')
    area_fraction = real(NPARTICLES * ldim(1) * ldim(2)) / real(micrograph_box * micrograph_box)
    write(logfhandle,'(a,f7.3,a,f7.3)') '>>> TEST_MINI_STREAM particle area fraction ', area_fraction, &
        &'; maximum ', MAX_AREA_FRACTION
    if( area_fraction > MAX_AREA_FRACTION )&
        &THROW_HARD('TEST_MINI_STREAM FAILED: simulated micrograph is too densely occupied')

    write(logfhandle,'(a,a,a,i0,a)') '>>> TEST_MINI_STREAM ', system_name%to_char(), ': generating ', &
        &NMICROGRAPHS, ' synthetic micrographs'
    do i = 1, NMICROGRAPHS
        call cline_sim%set('prg',       'simulate_movie')
        call cline_sim%set('stk',       particle_stack_path)
        call cline_sim%set('xdim',      micrograph_box)
        call cline_sim%set('ydim',      micrograph_box)
        call cline_sim%set('nframes',   NFRAMES)
        call cline_sim%set('smpd',      SMPD)
        call cline_sim%set('snr',       0.5)
        call cline_sim%set('kv',        KV)
        call cline_sim%set('cs',        CS)
        call cline_sim%set('fraca',     FRACA)
        call cline_sim%set('defocus',   1.5 + 0.25 * real(i - 1))
        call cline_sim%set('trs',       1.0)
        call cline_sim%set('nthr',      1)
        if( i == 1 )then
            call cline_sim%set('mkdir',   'yes')
            call cline_sim%set('dir_exec','simulate_micrographs')
        else
            call cline_sim%set('mkdir',   'no')
        endif
        call xsim_movie%execute(cline_sim)
        movie_name  = 'synthetic_micrograph_'//int2str_pad(i, 3)//MRC_EXT
        params_name = 'synthetic_micrograph_truth_'//int2str_pad(i, 3)//TXT_EXT
        call simple_rename(OPTIMAL_FILE, movie_name)
        call simple_rename(PARAMS_FILE, params_name)
        if( file_exists(MOVIE_FILE) ) call del_file(MOVIE_FILE)
        movie_paths(i)  = simple_abspath(movie_name)
        params_paths(i) = simple_abspath(params_name)
        call cline_sim%kill
    enddo
    call simple_chdir(fixture_root, status)
    if( status /= 0 ) THROW_HARD('TEST_MINI_STREAM FAILED: could not leave simulation directory')
    filetab_path = filepath(fixture_root, FILETAB_FILE)
    call write_filetable(filetab_path, movie_paths)

    call cline_mini%set('prg',            'mini_stream')
    call cline_mini%set('mkdir',                  'yes')
    call cline_mini%set('filetab',        filetab_path)
    call cline_mini%set('smpd',                   SMPD)
    call cline_mini%set('fraca',                 FRACA)
    call cline_mini%set('kv',                       KV)
    call cline_mini%set('cs',                       CS)
    call cline_mini%set('moldiam_max',         MSKDIAM)
    call cline_mini%set('pcontrast',            'black')
    call cline_mini%set('pick_roi',                'no')
    call cline_mini%set('ncls',                       2)
    call cline_mini%set('nptcls_per_cls',             10)
    call cline_mini%set('pspecsz',                   256)
    call cline_mini%set('nparts',                       1)
    call cline_mini%set('nthr',                         1)
    call xmini_stream%execute(cline_mini)
    call simple_getcwd(mini_dir)
    project_path = filepath(mini_dir, PROJFILE)
    if( .not. file_exists(project_path) ) THROW_HARD('TEST_MINI_STREAM FAILED: project was not created')
    call result%read(project_path)

    if( result%os_mic%get_noris() /= NMICROGRAPHS )&
        &THROW_HARD('TEST_MINI_STREAM FAILED: project has an unexpected number of micrographs')
    do i = 1, result%os_mic%get_noris()
        call result%os_mic%getter(i, 'intg', value)
        imic = find_simulation(value)
        if( imic == 0 ) THROW_HARD('TEST_MINI_STREAM FAILED: micrograph cannot be matched to simulation truth')
        call truth%new(1, is_ptcl=.false.)
        call truth%read(params_paths(imic))
        dfx       = result%os_mic%get_dfx(i)
        dfy       = result%os_mic%get_dfy(i)
        truth_dfx = truth%get_dfx(1)
        truth_dfy = truth%get_dfy(1)
        defocus_error = min(max(abs(dfx - truth_dfx), abs(dfy - truth_dfy)), &
            &max(abs(dfx - truth_dfy), abs(dfy - truth_dfx)))
        write(logfhandle,'(a,i0,a,f7.3,a,f7.3)') '>>> TEST_MINI_STREAM micrograph ', imic, &
            &': defocus max error ', defocus_error, ' um; limit ', DEFOCUS_TOL
        if( defocus_error > DEFOCUS_TOL )&
            &THROW_HARD('TEST_MINI_STREAM FAILED: CTF defocus error exceeds 0.10 um')
        call truth%kill
    enddo

    expected_particles = NMICROGRAPHS * NPARTICLES
    npicked = result%os_ptcl2D%get_noris()
    pick_ratio = real(npicked) / real(expected_particles)
    write(logfhandle,'(a,i0,a,i0,a,f7.3)') '>>> TEST_MINI_STREAM particles picked ', npicked, &
        &'; placed ', expected_particles, '; ratio ', pick_ratio
    if( pick_ratio < PICK_RECALL_MIN ) THROW_HARD('TEST_MINI_STREAM FAILED: fewer than 40 percent of placed particles were picked')
    if( pick_ratio > PICK_RATIO_MAX ) THROW_HARD('TEST_MINI_STREAM FAILED: pick count exceeds 150 percent of placed particles')
    if( result%os_ptcl3D%get_noris() /= npicked )&
        &THROW_HARD('TEST_MINI_STREAM FAILED: 2D and 3D particle counts differ')

    nclasses = result%os_cls2D%get_noris()
    if( nclasses < 1 ) THROW_HARD('TEST_MINI_STREAM FAILED: solve2D produced no classes')
    classes = result%os_ptcl2D%get_all_asint('class')
    if( size(classes) /= npicked .or. any(classes < 1) .or. any(classes > nclasses) )&
        &THROW_HARD('TEST_MINI_STREAM FAILED: particle class assignments are incomplete or invalid')
    populations = result%os_cls2D%get_all_asint('pop')
    if( size(populations) /= nclasses .or. any(populations < 0) )&
        &THROW_HARD('TEST_MINI_STREAM FAILED: class populations are incomplete or invalid')
    do i = 1, nclasses
        if( populations(i) /= count(classes == i) )&
            &THROW_HARD('TEST_MINI_STREAM FAILED: class populations disagree with particle assignments')
    enddo
    if( sum(populations) /= npicked )&
        &THROW_HARD('TEST_MINI_STREAM FAILED: class populations do not account for all particles')
    if( .not. result%os_cls2D%isthere('shape_rank') )&
        &THROW_HARD('TEST_MINI_STREAM FAILED: shape-ranking output is absent')
    shape_ranks = result%os_cls2D%get_all_asint('shape_rank')
    if( any(shape_ranks < 0) .or. any(shape_ranks > nclasses) )&
        &THROW_HARD('TEST_MINI_STREAM FAILED: shape ranks are outside the valid range')
    if( any((shape_ranks > 0) .and. (populations == 0)) )&
        &THROW_HARD('TEST_MINI_STREAM FAILED: an empty class received a shape rank')
    nranked = count(shape_ranks > 0)
    do rank = 1, nranked
        if( count(shape_ranks == rank) /= 1 )&
            &THROW_HARD('TEST_MINI_STREAM FAILED: positive shape ranks are not contiguous and unique')
    enddo
    write(logfhandle,'(a,i0,a,i0)') '>>> TEST_MINI_STREAM quantitatively consistent classes ', nclasses, &
        &'; class averages passing shape-quality selection ', nranked

    call result%get_cavgs_stk(cavg_stack, ncavgs, cavg_smpd)
    if( ncavgs /= nclasses ) THROW_HARD('TEST_MINI_STREAM FAILED: class-average stack count differs from project classes')
    call find_ldim_nptcls(cavg_stack, ldim, nimages)
    if( nimages /= nclasses ) THROW_HARD('TEST_MINI_STREAM FAILED: class-average file has an unexpected image count')
    call cavg%new([ldim(1), ldim(2), 1], cavg_smpd, wthreads=.false.)
    do i = 1, nimages
        call cavg%read(cavg_stack, i)
        cavg_variance = cavg%variance()
        write(logfhandle,'(a,i0,a,i0,a,es12.4)') '>>> TEST_MINI_STREAM class ', i, &
            &' population ', populations(i), '; variance ', cavg_variance
        if( .not. ieee_is_finite(cavg_variance) )&
            &THROW_HARD('TEST_MINI_STREAM FAILED: class average contains non-finite data')
        if( populations(i) > 0 .and. cavg_variance <= TINY )&
            &THROW_HARD('TEST_MINI_STREAM FAILED: populated class average is constant')
    enddo
    call cavg%kill
    call result%kill
    call cline_mini%kill
    call simple_chdir(cwd_saved, status)
    if( status /= 0 ) THROW_HARD('TEST_MINI_STREAM FAILED: could not restore original directory')
    write(logfhandle,'(a,a,a,i0,a,i0,a,i0,a)') 'PASS: mini_stream ', system_name%to_char(), ' validated ', &
        &NMICROGRAPHS, ' micrographs, ', npicked, ' particles, and ', nclasses, ' classes'
    enddo
    call simple_end('**** SIMPLE_TEST_MINI_STREAM WORKFLOW NORMAL STOP ****')

  contains

    integer function find_simulation( micrograph_path ) result(ind)
        type(string), intent(in) :: micrograph_path
        type(string) :: observed_name, truth_name
        integer :: isim
        ind = 0
        observed_name = basename(micrograph_path)
        do isim = 1, NMICROGRAPHS
            truth_name = basename(movie_paths(isim))
            if( trim(observed_name%to_char()) == trim(truth_name%to_char()) )then
                ind = isim
                return
            endif
        enddo
    end function find_simulation

end subroutine exec_test_mini_stream_quantitative

!> Hermetic quantitative validation: one embedded 6VXX volume is reprojected
!! and turned into CTF-affected particles. Counts, dimensions, sampling, image
!! variance, orientation diversity, shifts, and CTF metadata are checked.
module subroutine exec_test_simulate_particles( self, cline )
    use simple_atoms,         only: atoms
    use simple_molecule_data, only: molecule_data, sars_cov2_spkgp_6vxx
    use simple_imghead,       only: find_ldim_nptcls, find_img_smpd
    use simple_oris,          only: oris
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    class(commander_test_simulate_particles), intent(inout) :: self
    class(cmdline),                           intent(inout) :: cline
    real,    parameter                  :: SMPD        = 1.3
    real,    parameter                  :: MSKDIAM     = 180.
    real,    parameter                  :: KV          = 300.0
    real,    parameter                  :: CS          = 2.7
    real,    parameter                  :: FRACA       = 0.1
    real,    parameter                  :: DEFOCUS     = 2.0
    real,    parameter                  :: DFERR       = 0.2
    real,    parameter                  :: ASTIGERR    = 0.05
    real,    parameter                  :: SHIFT_LIMIT = 2.0
    integer, parameter                  :: NSPACE      = 100
    integer, parameter                  :: NPTCLS_SIM  = 200
    type(cmdline)                       :: cline_reproj, cline_sim
    type(parameters)                    :: params
    type(commander_reproject)           :: xreproject
    type(commander_simulate_particles)  :: xsim_ptcls
    type(atoms)                         :: molecule
    type(molecule_data)                 :: mol
    type(image)                         :: volume
    type(string)                        :: vol_file, cwd_saved, fixture_root
    integer                             :: ldim(3), nsections, status
    real                                :: smpd_vol, volume_variance
    logical                             :: all_ok
    call set_fixed_seed(20260925)
    call simple_getcwd(cwd_saved)
    fixture_root = filepath(cwd_saved, 'test_simulate_particles_'//int2str(get_process_id()))
    if( dir_exists(fixture_root) ) call simple_rmdir(fixture_root)
    call simple_mkdir(fixture_root)
    call simple_chdir(fixture_root, status)
    if( status /= 0 ) THROW_HARD('TEST_SIMULATE_PARTICLES FAILED: could not enter fixture directory')
    call params%new(cline)
    all_ok   = .true.
    vol_file = '6VXX.mrc'
    ! ---- one volume from the embedded coordinates ----
    write(logfhandle,'(a)') '>>> TEST_SIMULATE_PARTICLES: generating '//vol_file%to_char()
    mol = sars_cov2_spkgp_6vxx()
    call molecule%pdb2mrc(smpd=SMPD, volfile=vol_file, mol=mol, center_pdb=.true.)
    call molecule%kill()
    if( .not. file_exists(vol_file) ) THROW_HARD('TEST_SIMULATE_PARTICLES FAILED: volume not generated')
    call find_ldim_nptcls(vol_file, ldim, nsections)
    smpd_vol = find_img_smpd(vol_file)
    write(logfhandle,'(a,i4,a,i4,a,i4,a,f6.2,a,i0)') '    volume dims = [', ldim(1), ',', ldim(2), ',', &
        &ldim(3), ' ], smpd = ', smpd_vol, ', sections = ', nsections
    if( ldim(1) /= ldim(2) .or. ldim(1) /= ldim(3) .or. ldim(1) < 1 )then
        write(logfhandle,'(a)') '    FAIL: volume is not a cube'
        all_ok = .false.
    endif
    if( abs(smpd_vol - SMPD) > 0.01 )then
        write(logfhandle,'(a,f6.2,a,f6.2)') '    FAIL: volume smpd mismatch, expected ', SMPD, ' got ', smpd_vol
        all_ok = .false.
    endif
    if( all(ldim > 0) )then
        call volume%new(ldim, smpd_vol, wthreads=.false.)
        call volume%read(vol_file)
        volume_variance = volume%variance()
        write(logfhandle,'(a,es12.4)') '    volume variance = ', volume_variance
        if( .not. ieee_is_finite(volume_variance) .or. volume_variance <= TINY )then
            write(logfhandle,'(a)') '    FAIL: volume contains invalid or constant density'
            all_ok = .false.
        endif
        call volume%kill
    endif
    ! ---- reproject ----
    write(logfhandle,'(a)') '>>> TEST_SIMULATE_PARTICLES: reproject'
    call cline_reproj%set('prg',      'reproject')
    call cline_reproj%set('vol1',      vol_file)
    call cline_reproj%set('smpd',      SMPD)
    call cline_reproj%set('pgrp',      'c1')
    call cline_reproj%set('mskdiam',   MSKDIAM)
    call cline_reproj%set('nspace',    NSPACE)
    call cline_reproj%set('nthr',      params%nthr)
    call xreproject%execute(cline_reproj)
    call cline_reproj%kill()
    call check_stack(string('reprojs.mrcs'), NSPACE, 'reprojection')
    call check_oris(string('reproject_oris'//trim(TXT_EXT)), NSPACE, 'reprojection')
    ! ---- simulate particles ----
    write(logfhandle,'(a)') '>>> TEST_SIMULATE_PARTICLES: simulate_particles'
    call cline_sim%set('prg',      'simulate_particles')
    call cline_sim%set('vol1',      vol_file)
    call cline_sim%set('smpd',      SMPD)
    call cline_sim%set('mskdiam',   MSKDIAM)
    call cline_sim%set('nthr',      params%nthr)
    call cline_sim%set('nptcls',    NPTCLS_SIM)
    call cline_sim%set('pgrp',      'c1')
    call cline_sim%set('snr',       0.01)
    call cline_sim%set('ctf',       'yes')
    call cline_sim%set('kv',         KV)
    call cline_sim%set('cs',         CS)
    call cline_sim%set('fraca',      FRACA)
    call cline_sim%set('defocus',    DEFOCUS)
    call cline_sim%set('dferr',      DFERR)
    call cline_sim%set('astigerr',   ASTIGERR)
    call cline_sim%set('trs',        SHIFT_LIMIT)
    call cline_sim%set('even',       'yes')
    call xsim_ptcls%execute(cline_sim)
    call cline_sim%kill()
    call check_stack(string('simulated_particles.mrc'), NPTCLS_SIM, 'particle')
    call check_oris(string('simulated_oris'//trim(TXT_EXT)), NPTCLS_SIM, 'particle', check_ctf=.true.)
    ! ---- final verdict ----
    call simple_chdir(cwd_saved, status)
    if( status /= 0 ) THROW_HARD('TEST_SIMULATE_PARTICLES FAILED: could not restore original directory')
    if( all_ok )then
        write(logfhandle,'(a,i0,a,i0,a)') 'PASS: simulate_particles validated ', NSPACE, &
            &' reprojections and ', NPTCLS_SIM, ' CTF-affected particles'
        call simple_end('**** SIMPLE_TEST_SIMULATE_PARTICLES NORMAL STOP ****')
    else
        THROW_HARD('TEST_SIMULATE_PARTICLES FAILED')
    endif

  contains

    !> the stack exists, holds nexpected square images at the volume's smpd
    subroutine check_stack( fname, nexpected, what )
        type(string),     intent(in) :: fname
        integer,          intent(in) :: nexpected
        character(len=*), intent(in) :: what
        type(image) :: img
        integer :: iimg, ldim_stk(3), nimgs
        real    :: img_variance, max_variance, min_variance, smpd_stk
        write(logfhandle,'(a)') '>>> CHECK: '//what//' stack '//fname%to_char()
        if( .not. file_exists(fname) )then
            write(logfhandle,'(a)') '    FAIL: '//fname%to_char()//' not found'
            all_ok = .false.
            return
        endif
        call find_ldim_nptcls(fname, ldim_stk, nimgs)
        smpd_stk = find_img_smpd(fname)
        write(logfhandle,'(a,i6,a,i4,a,i4,a,f6.2)') '    images: ', nimgs, ', box: ', &
            &ldim_stk(1), ' x ', ldim_stk(2), ', smpd: ', smpd_stk
        if( nimgs /= nexpected )then
            write(logfhandle,'(a,i6,a,i6)') '    FAIL: expected ', nexpected, ' images, got ', nimgs
            all_ok = .false.
        endif
        if( ldim_stk(1) /= ldim_stk(2) .or. ldim_stk(1) < 1 )then
            write(logfhandle,'(a)') '    FAIL: invalid box dimensions'
            all_ok = .false.
        endif
        if( any(ldim_stk(1:2) /= ldim(1:2)) )then
            write(logfhandle,'(a)') '    FAIL: stack box differs from the source volume'
            all_ok = .false.
        endif
        if( abs(smpd_stk - SMPD) > 0.01 )then
            write(logfhandle,'(a,f6.2,a,f6.2)') '    FAIL: smpd mismatch, expected ', SMPD, ' got ', smpd_stk
            all_ok = .false.
        endif
        if( nimgs < 1 .or. ldim_stk(1) < 1 .or. ldim_stk(2) < 1 ) return
        call img%new([ldim_stk(1), ldim_stk(2), 1], smpd_stk, wthreads=.false.)
        min_variance = huge(1.0)
        max_variance = 0.0
        do iimg = 1, nimgs
            call img%read(fname, iimg)
            img_variance = img%variance()
            if( .not. ieee_is_finite(img_variance) .or. img_variance <= TINY )then
                write(logfhandle,'(a,i0,a,es12.4)') '    FAIL: image ', iimg, &
                    &' has invalid or zero variance ', img_variance
                all_ok = .false.
            else
                min_variance = min(min_variance, img_variance)
                max_variance = max(max_variance, img_variance)
            endif
        enddo
        call img%kill
        write(logfhandle,'(a,es12.4,a,es12.4)') '    image variance range: ', min_variance, ' to ', max_variance
    end subroutine check_stack

    !> The orientation table has one valid record per image. Particle metadata
    !! additionally has bounded CTF parameters and shifts matching the fixture.
    subroutine check_oris( fname, nexpected, what, check_ctf )
        type(string),     intent(in) :: fname
        integer,          intent(in) :: nexpected
        character(len=*), intent(in) :: what
        logical, optional, intent(in) :: check_ctf
        type(oris) :: metadata
        integer, allocatable :: states(:)
        real,    allocatable :: dfx(:), dfy(:), e1(:), e2(:), e3(:), xs(:), ys(:)
        real,    allocatable :: kvs(:), css(:), fracas(:)
        integer :: nrecs
        logical :: inspect_ctf
        write(logfhandle,'(a)') '>>> CHECK: '//what//' orientations '//fname%to_char()
        if( .not. file_exists(fname) )then
            write(logfhandle,'(a)') '    FAIL: '//fname%to_char()//' not found'
            all_ok = .false.
            return
        endif
        nrecs = nlines(fname)
        write(logfhandle,'(a,i0)') '    orientation records: ', nrecs
        if( nrecs /= nexpected )then
            write(logfhandle,'(a,i6,a,i6)') '    FAIL: expected ', nexpected, ' records, got ', nrecs
            all_ok = .false.
        endif
        if( nrecs < 1 )then
            return
        endif
        call metadata%new(nrecs, is_ptcl=.true.)
        call metadata%read(fname)
        states = metadata%get_all_asint('state')
        if( any(states /= 1) )then
            write(logfhandle,'(a,i0)') '    FAIL: inactive orientation records: ', count(states /= 1)
            all_ok = .false.
        endif
        if( metadata%isthere('e1') .and. metadata%isthere('e2') .and. metadata%isthere('e3') )then
            e1 = metadata%get_all('e1')
            e2 = metadata%get_all('e2')
            e3 = metadata%get_all('e3')
            if( any(.not. ieee_is_finite(e1)) .or. any(.not. ieee_is_finite(e2)) .or. &
                &any(.not. ieee_is_finite(e3)) )then
                write(logfhandle,'(a)') '    FAIL: Euler angles contain non-finite values'
                all_ok = .false.
            endif
            if( maxval(e2) - minval(e2) < 10.0 )then
                write(logfhandle,'(a)') '    FAIL: projection directions lack angular diversity'
                all_ok = .false.
            endif
        else
            write(logfhandle,'(a)') '    FAIL: orientation table lacks Euler angles'
            all_ok = .false.
        endif
        inspect_ctf = .false.
        if( present(check_ctf) ) inspect_ctf = check_ctf
        if( inspect_ctf )then
            if( metadata%isthere('dfx') .and. metadata%isthere('dfy') )then
                dfx = metadata%get_all('dfx')
                dfy = metadata%get_all('dfy')
                write(logfhandle,'(a,f7.3,a,f7.3)') '    dfx range: ', minval(dfx), ' to ', maxval(dfx)
                if( any(.not. ieee_is_finite(dfx)) .or. any(abs(dfx - DEFOCUS) > DFERR + 1.e-4) )then
                    write(logfhandle,'(a)') '    FAIL: dfx lies outside the requested defocus interval'
                    all_ok = .false.
                endif
                if( maxval(dfx) - minval(dfx) <= TINY )then
                    write(logfhandle,'(a)') '    FAIL: requested defocus variation was not generated'
                    all_ok = .false.
                endif
                if( any(.not. ieee_is_finite(dfy)) .or. any(abs(dfy - dfx) > ASTIGERR + 1.e-4) )then
                    write(logfhandle,'(a)') '    FAIL: dfy lies outside the requested astigmatism interval'
                    all_ok = .false.
                endif
                if( maxval(abs(dfy - dfx)) <= TINY )then
                    write(logfhandle,'(a)') '    FAIL: requested astigmatism variation was not generated'
                    all_ok = .false.
                endif
            else
                write(logfhandle,'(a)') '    FAIL: particle orientations lack CTF defocus values'
                all_ok = .false.
            endif
            if( metadata%isthere('kv') .and. metadata%isthere('cs') .and. metadata%isthere('fraca') )then
                kvs    = metadata%get_all('kv')
                css    = metadata%get_all('cs')
                fracas = metadata%get_all('fraca')
                if( any(abs(kvs - KV) > 0.01) .or. any(abs(css - CS) > 0.01) .or. &
                    &any(abs(fracas - FRACA) > 1.e-4) )then
                    write(logfhandle,'(a)') '    FAIL: particle CTF constants differ from their inputs'
                    all_ok = .false.
                endif
            else
                write(logfhandle,'(a)') '    FAIL: particle orientations lack CTF constants'
                all_ok = .false.
            endif
            if( metadata%isthere('x') .and. metadata%isthere('y') )then
                xs = metadata%get_all('x')
                ys = metadata%get_all('y')
                if( any(.not. ieee_is_finite(xs)) .or. any(.not. ieee_is_finite(ys)) .or. &
                    &any(abs(xs) > SHIFT_LIMIT + 1.e-4) .or. any(abs(ys) > SHIFT_LIMIT + 1.e-4) )then
                    write(logfhandle,'(a)') '    FAIL: particle shifts exceed the requested limit'
                    all_ok = .false.
                endif
                if( max(maxval(abs(xs)), maxval(abs(ys))) <= TINY )then
                    write(logfhandle,'(a)') '    FAIL: requested particle shifts were not generated'
                    all_ok = .false.
                endif
            else
                write(logfhandle,'(a)') '    FAIL: particle orientations lack shifts'
                all_ok = .false.
            endif
        endif
        call metadata%kill
    end subroutine check_oris

end subroutine exec_test_simulate_particles

end submodule simple_commanders_test_highlevel_stream
