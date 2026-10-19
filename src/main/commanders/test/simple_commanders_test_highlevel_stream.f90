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
    use simple_test_gate,     only: test_gate
    use simple_bench_io,      only: bench_report_fname
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
    ! The floors by the floor rule, from the profile of 2026-10-18 (20 seed offsets at 8 threads,
    ! scratch/keep/profiles/phase5/after_e; the Phase 1 profile of 2026-10-16 gave 0.10 um and pick
    ! floors 0.16 / 0.31, which this wider sample does not admit):
    ! - defocus error per micrograph: medians 0.002-0.003 um, worst 0.125 (6vxx micrograph 3: a rare
    !   wrong fit), admissible from 0.127 (ceiling 0.13); the farthest other truth lies 0.3 um and more
    !   away;
    ! - pick ratio: 6vxx median 0.333, s 0.093, worst 0.188, admissible up to 0.055 (floor 0.05);
    !   1jxy median 0.490, s 0.077, worst 0.396, admissible up to 0.258 (floor 0.25). The simulator
    !   saves no particle positions, so there is no control; the former 0.40 came from one unseeded run.
    real,             parameter :: DEFOCUS_TOL     = 0.13
    real,             parameter :: PICK_RECALL_MIN_6VXX = 0.05
    real,             parameter :: PICK_RECALL_MIN_1JXY = 0.25
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
    type(test_gate)                :: gate
    type(string)                   :: cwd_saved, fixture_root, mini_dir
    type(string)                   :: particle_stack_path, filetab_path, project_path
    type(string)                   :: movie_name, params_name, value, cavg_stack
    type(string)                   :: requested_suite, system_name, vol_file, pgrp
    type(string)                   :: movie_paths(NMICROGRAPHS), params_paths(NMICROGRAPHS)
    integer, allocatable           :: classes(:), populations(:), shape_ranks(:)
    integer                       :: i, imic, isuite, rank, status, ldim(3), nimages, micrograph_box
    integer                       :: npicked, nclasses, ncavgs, nranked, expected_particles
    real                          :: dfx, dfy, truth_dfx, truth_dfy, defocus_error
    real                          :: pick_ratio, cavg_smpd, cavg_variance, area_fraction, pick_recall_min, wrong_error
    logical                       :: ok

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
    call gate%new(filepath(cwd_saved, 'metrics.tsv'))

    do isuite = 1, 2
        if( isuite == 1 )then
            system_name = '6vxx'
            vol_file    = '6VXX.mrc'
            pgrp        = 'c3'
            micrograph_box = MICROGRAPH_BOX_6VXX
            pick_recall_min = PICK_RECALL_MIN_6VXX
        else
            system_name = '1jxy'
            vol_file    = '1JXY.mrc'
            pgrp        = 'c1'
            micrograph_box = MICROGRAPH_BOX_1JXY
            pick_recall_min = PICK_RECALL_MIN_1JXY
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
    call gate%metric(system_name%to_char()//'_particle_area_fraction', area_fraction, MAX_AREA_FRACTION, &
        &area_fraction <= MAX_AREA_FRACTION, upper=.true.)
    if( .not. area_fraction <= MAX_AREA_FRACTION )&
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

    call require('_micrographs_in_project', result%os_mic%get_noris() == NMICROGRAPHS, &
        &'project has an unexpected number of micrographs')
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
        ! min and max may drop a non-finite argument
        if( .not. (ieee_is_finite(dfx) .and. ieee_is_finite(dfy)) ) dfx = huge(1.)
        defocus_error = min(max(abs(dfx - truth_dfx), abs(dfy - truth_dfy)), &
            &max(abs(dfx - truth_dfy), abs(dfy - truth_dfx)))
        write(logfhandle,'(a,i0,a,f7.3,a,f7.3)') '>>> TEST_MINI_STREAM micrograph ', imic, &
            &': defocus max error ', defocus_error, ' um; limit ', DEFOCUS_TOL
        call gate%metric(system_name%to_char()//'_defocus_error_um_mic'//int2str(imic), defocus_error, DEFOCUS_TOL, &
            &defocus_error <= DEFOCUS_TOL, upper=.true.)
        if( .not. defocus_error <= DEFOCUS_TOL )&
            &THROW_HARD('TEST_MINI_STREAM FAILED: CTF defocus error exceeds 0.13 um')
        call truth%kill
        ! the negative control: the fit against the truth of the other micrograph whose truth lies
        ! farthest from this one's (the neighbour's came within the ceiling in some runs)
        wrong_error = farthest_truth_error(imic, truth_dfx, truth_dfy, dfx, dfy)
        call gate%check(system_name%to_char()//'_control_defocus_mic'//int2str(imic)//'_against_the_farthest_truth_fails', &
            &.not. (wrong_error <= DEFOCUS_TOL))
    enddo

    expected_particles = NMICROGRAPHS * NPARTICLES
    npicked = result%os_ptcl2D%get_noris()
    pick_ratio = real(npicked) / real(expected_particles)
    write(logfhandle,'(a,i0,a,i0,a,f7.3)') '>>> TEST_MINI_STREAM particles picked ', npicked, &
        &'; placed ', expected_particles, '; ratio ', pick_ratio
    call gate%metric(system_name%to_char()//'_pick_ratio', pick_ratio, pick_recall_min, &
        &pick_ratio >= pick_recall_min)
    call gate%metric(system_name%to_char()//'_pick_ratio_against_ceiling', pick_ratio, PICK_RATIO_MAX, &
        &pick_ratio <= PICK_RATIO_MAX, upper=.true.)
    if( .not. pick_ratio >= pick_recall_min ) THROW_HARD('TEST_MINI_STREAM FAILED: too few placed particles were picked')
    if( .not. pick_ratio <= PICK_RATIO_MAX ) THROW_HARD('TEST_MINI_STREAM FAILED: pick count exceeds 150 percent of placed particles')
    call require('_particle_counts_2D_3D_agree', result%os_ptcl3D%get_noris() == npicked, '2D and 3D particle counts differ')

    nclasses = result%os_cls2D%get_noris()
    call gate%metric(system_name%to_char()//'_classes', real(nclasses), 1., nclasses >= 1)
    if( nclasses < 1 ) THROW_HARD('TEST_MINI_STREAM FAILED: solve2D produced no classes')
    classes = result%os_ptcl2D%get_all_asint('class')
    ok = size(classes) == npicked
    if( ok ) ok = all(classes >= 1 .and. classes <= nclasses)
    call require('_class_assignments_valid', ok, 'particle class assignments are incomplete or invalid')
    populations = result%os_cls2D%get_all_asint('pop')
    ok = size(populations) == nclasses
    if( ok ) ok = all(populations >= 0)
    call require('_class_populations_valid', ok, 'class populations are incomplete or invalid')
    ok = .true.
    do i = 1, nclasses
        ok = ok .and. populations(i) == count(classes == i)
    enddo
    call require('_class_populations_match_assignments', ok .and. sum(populations) == npicked, &
        &'class populations disagree with particle assignments')
    call require('_shape_ranks_written', result%os_cls2D%isthere('shape_rank'), 'shape-ranking output is absent')
    shape_ranks = result%os_cls2D%get_all_asint('shape_rank')
    ok = all(shape_ranks >= 0 .and. shape_ranks <= nclasses)
    if( ok ) ok = .not. any((shape_ranks > 0) .and. (populations == 0))
    nranked = count(shape_ranks > 0)
    do rank = 1, nranked
        ok = ok .and. count(shape_ranks == rank) == 1
    enddo
    call require('_shape_ranks_valid', ok, 'shape ranks are out of range, unique or contiguous, or rank an empty class')
    write(logfhandle,'(a,i0,a,i0)') '>>> TEST_MINI_STREAM quantitatively consistent classes ', nclasses, &
        &'; class averages passing shape-quality selection ', nranked
    call gate%report(system_name%to_char()//'_shape_ranked_classes', real(nranked))

    call result%get_cavgs_stk(cavg_stack, ncavgs, cavg_smpd)
    call find_ldim_nptcls(cavg_stack, ldim, nimages)
    call require('_class_average_counts', ncavgs == nclasses .and. nimages == nclasses, &
        &'class-average stack count differs from project classes')
    call cavg%new([ldim(1), ldim(2), 1], cavg_smpd, wthreads=.false.)
    ok = .true.
    do i = 1, nimages
        call cavg%read(cavg_stack, i)
        cavg_variance = cavg%variance()
        write(logfhandle,'(a,i0,a,i0,a,es12.4)') '>>> TEST_MINI_STREAM class ', i, &
            &' population ', populations(i), '; variance ', cavg_variance
        ok = ok .and. ieee_is_finite(cavg_variance)
        if( populations(i) > 0 ) ok = ok .and. cavg_variance > TINY
    enddo
    call require('_class_averages_finite_and_populated_ones_not_constant', ok, &
        &'a class average is not finite, or a populated one is constant')
    ! the path: solve2D ran and recorded its refine2D iterations (the strategies' records)
    call require('_path_solve2D_refine2D_iterations_recorded', file_exists(bench_report_fname('SOLVE2D')) .and. &
        &nrecords(bench_report_fname('REFINE2D_STRATEGY')) > 0, 'solve2D left no record of its refine2D iterations')
    call cavg%kill
    call result%kill
    call cline_mini%kill
    call simple_chdir(cwd_saved, status)
    if( status /= 0 ) THROW_HARD('TEST_MINI_STREAM FAILED: could not restore original directory')
    write(logfhandle,'(a,a,a,i0,a,i0,a,i0,a)') 'PASS: mini_stream ', system_name%to_char(), ' validated ', &
        &NMICROGRAPHS, ' micrographs, ', npicked, ' particles, and ', nclasses, ' classes'
    enddo
    call gate%kill
    call simple_end('**** SIMPLE_TEST_MINI_STREAM WORKFLOW NORMAL STOP ****')

  contains

    !> a consistency gate: its row, then the stop when it fails
    subroutine require( name, ok, what )
        character(len=*), intent(in) :: name, what
        logical,          intent(in) :: ok
        character(len=:), allocatable :: msg
        call gate%check(system_name%to_char()//name, ok)
        if( .not. ok )then
            msg = 'TEST_MINI_STREAM FAILED: '//what
            THROW_HARD(msg)
        endif
    end subroutine require

    !> the error of the fit (fdx, fdy) against the truth of the other micrograph whose truth lies
    !! farthest from micrograph imic's (tdx, tdy), the astigmatic axes unordered as in the error itself
    real function farthest_truth_error( imic, tdx, tdy, fdx, fdy ) result( err )
        integer, intent(in) :: imic
        real,    intent(in) :: tdx, tdy, fdx, fdy
        type(oris) :: other
        real    :: odx, ody, d, dmax
        integer :: j
        err  = 0.
        dmax = -1.
        do j = 1, NMICROGRAPHS
            if( j == imic ) cycle
            call other%new(1, is_ptcl=.false.)
            call other%read(params_paths(j))
            odx = other%get_dfx(1)
            ody = other%get_dfy(1)
            call other%kill
            d = 0.5 * (odx + ody) - 0.5 * (tdx + tdy)
            if( abs(d) > dmax )then
                dmax = abs(d)
                err  = min(max(abs(fdx - odx), abs(fdy - ody)), max(abs(fdx - ody), abs(fdy - odx)))
            endif
        enddo
    end function farthest_truth_error

    !> the number of records in a benchmark record file (0 when it is absent)
    integer function nrecords( fname )
        type(string), intent(in) :: fname
        character(len=256) :: line
        integer :: funit, ios
        nrecords = 0
        if( .not. file_exists(fname) ) return
        open(newunit=funit, file=fname%to_char(), status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        do
            read(funit, '(a)', iostat=ios) line
            if( ios /= 0 ) exit
            if( index(line, 'record iteration') == 1 ) nrecords = nrecords + 1
        enddo
        close(funit)
    end function nrecords

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
    use simple_test_gate,     only: test_gate
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
    type(test_gate)                     :: gate
    type(string)                        :: vol_file, cwd_saved, fixture_root
    integer                             :: ldim(3), nsections, status
    real                                :: smpd_vol, volume_variance
    logical                             :: all_ok, ok_before
    call set_fixed_seed(20260925)
    call simple_getcwd(cwd_saved)
    call gate%new(filepath(cwd_saved, 'metrics.tsv'))
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
    if( .not. (abs(smpd_vol - SMPD) <= 0.01) )then
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
    call gate%check('volume_valid', all_ok)
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
    call gate_check_stack(string('reprojs.mrcs'), NSPACE, 'reprojection')
    call gate_check_oris(string('reproject_oris'//trim(TXT_EXT)), NSPACE, 'reprojection', .false.)
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
    call gate_check_stack(string('simulated_particles.mrc'), NPTCLS_SIM, 'particle')
    call gate_check_oris(string('simulated_oris'//trim(TXT_EXT)), NPTCLS_SIM, 'particle', .true.)
    ! ---- the images against the reprojections at their recorded poses ----
    call check_images_against_projections
    ! ---- several states: exact counts, each image from its own state's volume ----
    call check_mixture('0.5,0.3,0.2',       90,  [45, 27, 18], 'mixture_532')
    call check_mixture('0.333,0.333,0.334', 100, [33, 33, 34], 'mixture_equal')
    all_ok = all_ok .and. gate%passed()
    ! ---- final verdict ----
    call gate%kill
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

    !> simulate_particles with three blob volumes (no CTF, no shift, high SNR) and the population
    !! fractions pops: the state counts are exact (each fraction times n rounded, the remainder to
    !! the largest fractions), and every image correlates best with the projection of its own
    !! state's volume at its true pose
    !> 50 particles of the 6VXX map without CTF or shift at SNR 2: every image correlates better with
    !! the reprojection at its own recorded pose than with the reprojection at another particle's pose
    !! (the negative control, whose mean correlation is reported)
    subroutine check_images_against_projections
        use simple_simple_volinterp, only: reproject
        use simple_rnd,              only: shuffle
        integer, parameter :: NIMG = 50
        class(cmdline), allocatable :: cl
        type(oris)    :: os
        type(image)   :: vol, img
        type(image), allocatable :: projs(:)
        integer, allocatable :: perm(:)
        integer :: i, nwrong, ldimv(3), nsec
        real    :: cc_own, cc_other, sum_own, sum_other
        write(logfhandle,'(a)') '>>> TEST_SIMULATE_PARTICLES: images against their reprojections'
        call simple_mkdir('images_vs_projections')
        call simple_chdir('images_vs_projections')
        allocate(cl)
        call cl%set('prg',     'simulate_particles')
        call cl%set('vol1',    filepath(fixture_root, vol_file))
        call cl%set('smpd',    SMPD)
        call cl%set('mskdiam', MSKDIAM)
        call cl%set('nptcls',  NIMG)
        call cl%set('pgrp',    'c1')
        call cl%set('snr',     2.)
        call cl%set('ctf',     'no')
        call cl%set('sherr',   0.)
        call cl%set('nthr',    params%nthr)
        call set_fixed_seed(20261018)
        call xsim_ptcls%execute(cl)
        call cl%kill
        deallocate(cl)
        call os%new(NIMG, is_ptcl=.true.)
        call os%read(string('simulated_oris'//trim(TXT_EXT)), [1,NIMG])
        call find_ldim_nptcls(filepath(fixture_root, vol_file), ldimv, nsec)
        call vol%new(ldimv, SMPD)
        call vol%read(filepath(fixture_root, vol_file))
        projs = reproject(vol, os)
        perm  = [(i, i = 1, NIMG)]
        call shuffle(perm)
        call img%new([ldimv(1), ldimv(2), 1], SMPD)
        nwrong    = 0
        sum_own   = 0.
        sum_other = 0.
        do i = 1, NIMG
            call img%read(string('simulated_particles.mrc'), i)
            cc_own   = img%real_corr(projs(i))
            ! a partner other than itself
            if( perm(i) == i ) perm(i) = modulo(i, NIMG) + 1
            cc_other = img%real_corr(projs(perm(i)))
            sum_own   = sum_own + cc_own
            sum_other = sum_other + cc_other
            if( .not. (cc_own > cc_other) ) nwrong = nwrong + 1
        enddo
        call gate%report('images_mean_corr_with_own_pose_reprojection', sum_own / real(NIMG))
        call gate%report('control_images_mean_corr_with_another_pose_reprojection', sum_other / real(NIMG))
        call gate%metric('images_closer_to_another_pose_reprojection', real(nwrong), 0., nwrong == 0, upper=.true.)
        if( nwrong > 0 ) all_ok = .false.
        do i = 1, NIMG
            call projs(i)%kill
        enddo
        call vol%kill
        call img%kill
        call os%kill
        call simple_chdir(fixture_root)
    end subroutine check_images_against_projections

    subroutine check_mixture( pops, n, counts, tag )
        use simple_flex_pca_application_tester, only: write_phantom_truths
        use simple_simple_volinterp,            only: reproject
        character(len=*), intent(in) :: pops, tag
        integer,          intent(in) :: n, counts(3)
        type(string)  :: truths(3), tmean, tdiff, dir
        class(cmdline), allocatable :: cl
        type(oris)    :: os
        type(image)   :: vols(3), img
        type(image), allocatable :: projs(:,:)
        integer :: k, i, nwrong, ldimv(3), nsec
        integer, allocatable :: states(:)
        real    :: cc(3)
        write(logfhandle,'(a)') '>>> TEST_SIMULATE_PARTICLES: mixture '//tag
        call simple_mkdir(tag)
        call simple_chdir(tag)
        call simple_getcwd(dir)
        call write_phantom_truths(dir, 3, truths, tmean, tdiff)
        allocate(cl)
        do k = 1, 3
            call cl%set('vol'//int2str(k), truths(k))
        enddo
        call cl%set('prg',        'simulate_particles')
        call cl%set('nstates',    3)
        call cl%set('state_pops', pops)
        call cl%set('smpd',       2.0)
        call cl%set('mskdiam',    100.)
        call cl%set('nptcls',     n)
        call cl%set('pgrp',       'c1')
        call cl%set('snr',        2.)
        call cl%set('ctf',        'no')
        call cl%set('sherr',      0.)
        call cl%set('nthr',       params%nthr)
        call set_fixed_seed(20261017)
        call xsim_ptcls%execute(cl)
        call cl%kill
        deallocate(cl)
        call os%new(n, is_ptcl=.true.)
        call os%read(string('simulated_oris'//trim(TXT_EXT)), [1,n])
        states = nint(os%get_all('state'))
        do k = 1, 3
            call gate%metric(tag//'_state'//int2str(k)//'_count', real(count(states == k)), real(counts(k)), &
                &count(states == k) == counts(k))
        enddo
        call find_ldim_nptcls(truths(1), ldimv, nsec)
        allocate(projs(n,3))
        do k = 1, 3
            call vols(k)%new(ldimv, 2.0)
            call vols(k)%read(truths(k))
            projs(:,k) = reproject(vols(k), os)
        enddo
        call img%new([ldimv(1), ldimv(2), 1], 2.0)
        nwrong = 0
        do i = 1, n
            call img%read(string('simulated_particles.mrc'), i)
            do k = 1, 3
                cc(k) = img%real_corr(projs(i,k))
            enddo
            ! maxloc passes over a non-finite correlation
            if( .not. all(ieee_is_finite(cc)) )then
                nwrong = nwrong + 1
            else if( maxloc(cc, dim=1) /= states(i) )then
                nwrong = nwrong + 1
            endif
        enddo
        call gate%metric(tag//'_images_closer_to_another_state', real(nwrong), 0., nwrong == 0, upper=.true.)
        if( any([(count(states == k) /= counts(k), k = 1, 3)]) .or. nwrong > 0 ) all_ok = .false.
        do k = 1, 3
            call vols(k)%kill
            do i = 1, n
                call projs(i,k)%kill
            enddo
        enddo
        call img%kill
        call os%kill
        call simple_chdir(fixture_root)
    end subroutine check_mixture

    !> check_stack, and its verdict as a row of the gate
    subroutine gate_check_stack( fname, nexpected, what )
        type(string),     intent(in) :: fname
        integer,          intent(in) :: nexpected
        character(len=*), intent(in) :: what
        ok_before = all_ok
        all_ok    = .true.
        call check_stack(fname, nexpected, what)
        call gate%check(what//'_stack_valid', all_ok)
        all_ok = ok_before .and. all_ok
    end subroutine gate_check_stack

    !> check_oris, and its verdict as a row of the gate
    subroutine gate_check_oris( fname, nexpected, what, check_ctf )
        type(string),     intent(in) :: fname
        integer,          intent(in) :: nexpected
        character(len=*), intent(in) :: what
        logical,          intent(in) :: check_ctf
        ok_before = all_ok
        all_ok    = .true.
        call check_oris(fname, nexpected, what, check_ctf)
        call gate%check(what//'_orientations_valid', all_ok)
        all_ok = ok_before .and. all_ok
    end subroutine gate_check_oris

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
        if( .not. (abs(smpd_stk - SMPD) <= 0.01) )then
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
            call gate%metric(what//'_e2_range_deg', maxval(e2) - minval(e2), 10.0, &
                &maxval(e2) - minval(e2) >= 10.0 .and. all(ieee_is_finite(e2)))
            if( .not. (maxval(e2) - minval(e2) >= 10.0) )then
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
                ! the simulator's contract: every draw inside its requested interval (the largest of
                ! many uniform draws lies just inside the bound by construction, so these are checks,
                ! not floors; the largest deviations are reported)
                call gate%report(what//'_max_defocus_deviation_um', maxval(abs(dfx - DEFOCUS)))
                call gate%check(what//'_defocus_within_requested_interval', all(abs(dfx - DEFOCUS) <= DFERR + 1.e-4))
                call gate%report(what//'_max_astigmatism_um', maxval(abs(dfy - dfx)))
                call gate%check(what//'_astigmatism_within_requested_interval', all(abs(dfy - dfx) <= ASTIGERR + 1.e-4))
                if( .not. all(abs(dfx - DEFOCUS) <= DFERR + 1.e-4) )then
                    write(logfhandle,'(a)') '    FAIL: dfx lies outside the requested defocus interval'
                    all_ok = .false.
                endif
                if( maxval(dfx) - minval(dfx) <= TINY )then
                    write(logfhandle,'(a)') '    FAIL: requested defocus variation was not generated'
                    all_ok = .false.
                endif
                if( .not. all(abs(dfy - dfx) <= ASTIGERR + 1.e-4) )then
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
                if( .not. (all(abs(kvs - KV) <= 0.01) .and. all(abs(css - CS) <= 0.01) .and. &
                    &all(abs(fracas - FRACA) <= 1.e-4)) )then
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
                call gate%report(what//'_max_shift_px', max(maxval(abs(xs)), maxval(abs(ys))))
                call gate%check(what//'_shifts_within_requested_limit', &
                    &all(abs(xs) <= SHIFT_LIMIT + 1.e-4) .and. all(abs(ys) <= SHIFT_LIMIT + 1.e-4))
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
