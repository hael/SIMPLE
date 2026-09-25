!@descr: SINGLE test commanders: the nanoparticle atoms pipeline and the SINGLE workflow
module simple_commanders_test_single
use simple_commanders_api
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_single_atoms_stats
  contains
    procedure :: execute      => exec_test_single_atoms_stats
end type commander_test_single_atoms_stats

type, extends(commander_base) :: commander_test_single_workflow
  contains
    procedure :: execute      => exec_test_single_workflow
end type commander_test_single_workflow

integer, parameter :: BOX          = 160
integer, parameter :: MOLDIAM      = 20

contains

subroutine exec_test_single_atoms_stats( self, cline )
    use simple_commanders_atoms, only: commander_detect_atoms
    use simple_commanders_sim,   only: commander_simulate_nanoparticle
    use simple_commanders_atoms, only: commander_atoms_stats
    use simple_atoms,            only: atoms
    use simple_imghead,          only: find_ldim_nptcls, find_img_smpd
    use simple_test_utils,       only: set_fixed_seed
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    class(commander_test_single_atoms_stats), intent(inout) :: self
    class(cmdline),                           intent(inout) :: cline
    character(len=*), parameter :: SIM_VOL       = 'outvol.mrc'
    character(len=*), parameter :: SIM_PDB       = 'simatms.pdb'
    character(len=*), parameter :: DET_CC_VOL    = 'outvol_CC.mrc'
    character(len=*), parameter :: DET_SIM_VOL   = 'outvol_SIM.mrc'
    character(len=*), parameter :: DET_PDB       = 'outvol_ATMS.pdb'
    character(len=*), parameter :: ATOM_STATS    = 'atoms_stats.csv'
    character(len=*), parameter :: NP_STATS      = 'nanoparticle_stats.csv'
    character(len=*), parameter :: CN_STATS      = 'cn_dependent_stats.csv'
    integer,          parameter :: CN_MIN         = 3
    integer,          parameter :: CN_MAX         = 13
    real,             parameter :: POSITION_TOL  = 1.0
    real,             parameter :: MIN_RECALL    = 0.90
    real,             parameter :: MIN_PRECISION = 0.90
    real,             parameter :: MAX_RMS_ERROR = 0.75
    real,             parameter :: MIN_DENSITY_CORR = 0.80
    type(cmdline)                         :: cline_sim, cline_detat, cline_atstats
    type(parameters)                      :: params
    type(commander_simulate_nanoparticle) :: xsim_nptcl
    type(commander_detect_atoms)          :: xdetat
    type(commander_atoms_stats)           :: xatstats
    type(atoms)                           :: simulated_atoms, detected_atoms
    type(string)                          :: cwd_saved, fixture_root
    integer                               :: status, nsimulated, ndetected
    real                                  :: recall, precision, rms_error, density_corr
    logical                               :: all_ok
    write(logfhandle,'(a)') '>>> TEST_SINGLE_ATOMS_STATS:'
    call set_fixed_seed(20260925)
    call simple_getcwd(cwd_saved)
    fixture_root = filepath(cwd_saved, 'test_single_atoms_stats_'//int2str(get_process_id()))
    if( dir_exists(fixture_root) ) call simple_rmdir(fixture_root)
    call simple_mkdir(fixture_root)
    call simple_chdir(fixture_root, status)
    if( status /= 0 ) THROW_HARD('TEST_SINGLE_ATOMS_STATS FAILED: could not enter fixture directory')
    if( .not. cline%defined('smpd') )    call cline%set('smpd', 0.358)
    if( .not. cline%defined('element') ) call cline%set('element', 'Pt')
    call params%new(cline)
    all_ok = .true.
    density_corr = 0.0
    call cline_sim%set('prg',      'simulate_nanoparticle')
    call cline_sim%set('box',                          BOX)
    call cline_sim%set('smpd',                 params%smpd)
    call cline_sim%set('moldiam',                  MOLDIAM)
    call cline_sim%set('element',           params%element)
    call cline_sim%set('nthr',                 params%nthr)
    call cline_sim%set('outvol',                    SIM_VOL)
    call cline_sim%set('pdbout',                    SIM_PDB)
    call xsim_nptcl%execute(cline_sim)
    call cline_sim%kill()
    call check_volume(SIM_VOL, 'simulated')
    if( file_exists(SIM_PDB) )then
        call simulated_atoms%new(string(SIM_PDB))
        nsimulated = simulated_atoms%get_n()
        write(logfhandle,'(a,i0)') '    simulated atoms: ', nsimulated
        if( nsimulated < 1 ) call fail('simulator produced no atoms')
    else
        nsimulated = 0
        call fail(SIM_PDB//' was not generated')
    endif
    call cline_detat%set('prg',             'detect_atoms')
    call cline_detat%set('vol1',                    SIM_VOL)
    call cline_detat%set('smpd',               params%smpd)
    call cline_detat%set('element',         params%element)
    call cline_detat%set('nthr',                      params%nthr)
    call xdetat%execute(cline_detat)
    call cline_detat%kill()
    call check_volume(DET_CC_VOL, 'connected-component')
    call check_volume(DET_SIM_VOL, 'detected-atom density')
    call compare_density_volumes
    if( file_exists(DET_PDB) )then
        call detected_atoms%new(string(DET_PDB))
        ndetected = detected_atoms%get_n()
        write(logfhandle,'(a,i0)') '    detected atoms:  ', ndetected
        if( ndetected < 1 ) call fail('atom detector produced no atoms')
    else
        ndetected = 0
        call fail(DET_PDB//' was not generated')
    endif
    if( nsimulated > 0 .and. ndetected > 0 ) call compare_atom_coordinates
    call cline_atstats%set('prg',             'atoms_stats')
    call cline_atstats%set('vol1',                  SIM_VOL)
    call cline_atstats%set('vol2',               DET_CC_VOL)
    call cline_atstats%set('pdbfile',               DET_PDB)
    call cline_atstats%set('smpd',              params%smpd)
    call cline_atstats%set('element',        params%element)
    call cline_atstats%set('nthr',                     params%nthr)
    call xatstats%execute(cline_atstats)
    call cline_atstats%kill()
    call check_statistics_files(ndetected)
    call simulated_atoms%kill
    call detected_atoms%kill
    call simple_chdir(cwd_saved, status)
    if( status /= 0 ) THROW_HARD('TEST_SINGLE_ATOMS_STATS FAILED: could not restore original directory')
    if( all_ok )then
        write(logfhandle,'(a,f6.3,a,f6.3,a,f6.3,a,f6.3)') 'PASS: single_atoms_stats recall=', recall, &
            &', precision=', precision, ', RMS position error=', rms_error, ' A, density correlation=', density_corr
        call simple_end('**** SIMPLE_TEST_SINGLE_ATOMS_STATS NORMAL STOP ****')
    else
        THROW_HARD('TEST_SINGLE_ATOMS_STATS FAILED')
    endif

  contains

    subroutine fail( message )
        character(len=*), intent(in) :: message
        write(logfhandle,'(a)') '    FAIL: '//trim(message)
        all_ok = .false.
    end subroutine fail

    subroutine check_volume( filename, description )
        character(len=*), intent(in) :: filename, description
        type(image) :: vol
        integer     :: ldim(3), nsections
        real        :: smpd, variance
        if( .not. file_exists(filename) )then
            call fail(trim(description)//' volume '//filename//' was not generated')
            return
        endif
        call find_ldim_nptcls(string(filename), ldim, nsections)
        smpd = find_img_smpd(string(filename))
        write(logfhandle,'(a,a,a,3(i0,1x),a,i0,a,f7.3)') '    ', trim(description), &
            &' volume: dimensions ', ldim, ', sections ', nsections, ', smpd ', smpd
        if( any(ldim /= [BOX, BOX, BOX]) ) call fail(trim(description)//' volume dimensions are incorrect')
        if( abs(smpd - params%smpd) > 0.001 ) call fail(trim(description)//' volume sampling is incorrect')
        if( all(ldim > 0) )then
            call vol%new(ldim, smpd, wthreads=.false.)
            call vol%read(string(filename))
            variance = vol%variance()
            call vol%kill
            write(logfhandle,'(a,es12.4)') '        density variance: ', variance
            if( .not. ieee_is_finite(variance) .or. variance <= TINY )then
                call fail(trim(description)//' volume has invalid or constant density')
            endif
        endif
    end subroutine check_volume

    subroutine compare_density_volumes
        type(image) :: simulated_density, detected_density
        integer     :: simulated_ldim(3), detected_ldim(3), nsections
        real        :: simulated_smpd, detected_smpd
        if( .not. file_exists(SIM_VOL) .or. .not. file_exists(DET_SIM_VOL) ) return
        call find_ldim_nptcls(string(SIM_VOL), simulated_ldim, nsections)
        call find_ldim_nptcls(string(DET_SIM_VOL), detected_ldim, nsections)
        simulated_smpd = find_img_smpd(string(SIM_VOL))
        detected_smpd  = find_img_smpd(string(DET_SIM_VOL))
        if( any(simulated_ldim /= detected_ldim) )then
            call fail('simulated and detected-atom density dimensions differ')
            return
        endif
        if( abs(simulated_smpd - detected_smpd) > 0.001 )then
            call fail('simulated and detected-atom density sampling differs')
            return
        endif
        call simulated_density%new(simulated_ldim, simulated_smpd, wthreads=.false.)
        call detected_density%new(detected_ldim, detected_smpd, wthreads=.false.)
        call simulated_density%read(string(SIM_VOL))
        call detected_density%read(string(DET_SIM_VOL))
        density_corr = simulated_density%real_corr(detected_density)
        call simulated_density%kill
        call detected_density%kill
        write(logfhandle,'(a,f7.4)') '    simulated/detected density correlation: ', density_corr
        if( .not. ieee_is_finite(density_corr) )then
            call fail('simulated/detected density correlation is not finite')
        elseif( density_corr < MIN_DENSITY_CORR )then
            call fail('simulated/detected density correlation is below 0.80')
        endif
    end subroutine compare_density_volumes

    subroutine compare_atom_coordinates
        integer :: i, j, nsim_matched, ndet_matched
        real    :: distance, min_distance, squared_error
        nsim_matched = 0
        ndet_matched = 0
        squared_error = 0.0
        do i = 1, nsimulated
            min_distance = huge(1.0)
            do j = 1, ndetected
                distance = norm2(simulated_atoms%get_coord(i) - detected_atoms%get_coord(j))
                min_distance = min(min_distance, distance)
            enddo
            if( min_distance <= POSITION_TOL )then
                nsim_matched = nsim_matched + 1
                squared_error = squared_error + min_distance**2
            endif
        enddo
        do i = 1, ndetected
            min_distance = huge(1.0)
            do j = 1, nsimulated
                distance = norm2(detected_atoms%get_coord(i) - simulated_atoms%get_coord(j))
                min_distance = min(min_distance, distance)
            enddo
            if( min_distance <= POSITION_TOL ) ndet_matched = ndet_matched + 1
        enddo
        recall    = real(nsim_matched) / real(nsimulated)
        precision = real(ndet_matched) / real(ndetected)
        if( nsim_matched > 0 )then
            rms_error = sqrt(squared_error / real(nsim_matched))
        else
            rms_error = huge(1.0)
        endif
        write(logfhandle,'(a,f7.4,a,f7.4,a,f7.4,a)') '    coordinate match: recall=', recall, &
            &', precision=', precision, ', RMS error=', rms_error, ' A'
        if( recall < MIN_RECALL ) call fail('atom-detection recall is below 0.90')
        if( precision < MIN_PRECISION ) call fail('atom-detection precision is below 0.90')
        if( rms_error > MAX_RMS_ERROR ) call fail('atom RMS position error exceeds 0.75 A')
    end subroutine compare_atom_coordinates

    subroutine check_statistics_files( expected_atoms )
        integer, intent(in) :: expected_atoms
        integer            :: funit, ios
        real               :: csv_natoms, csv_naniso, csv_diameter
        character(len=STDLEN) :: header
        call check_line_count(ATOM_STATS, expected_atoms + 1)
        call check_line_count(NP_STATS, 2)
        call check_coordination_statistics(expected_atoms)
        if( .not. file_exists(NP_STATS) ) return
        ios = 0
        call fopen(funit, file=string(NP_STATS), status='old', action='read', iostat=ios)
        if( ios /= 0 )then
            call fail('could not open '//NP_STATS)
            return
        endif
        read(funit,'(a)',iostat=ios) header
        if( ios == 0 ) read(funit,*,iostat=ios) csv_natoms, csv_naniso, csv_diameter
        call fclose(funit)
        if( ios /= 0 )then
            call fail(NP_STATS//' does not contain a readable statistics row')
            return
        endif
        write(logfhandle,'(a,f8.1,a,f8.1,a,f8.3)') '    CSV summary: atoms=', csv_natoms, &
            &', anisotropic atoms=', csv_naniso, ', diameter=', csv_diameter
        if( .not. ieee_is_finite(csv_natoms) )then
            call fail(NP_STATS//' contains an invalid atom count')
        elseif( nint(csv_natoms) /= expected_atoms )then
            call fail(NP_STATS//' atom count disagrees with detected PDB')
        endif
        if( .not. ieee_is_finite(csv_naniso) .or. csv_naniso < 0.0 .or. csv_naniso > csv_natoms )then
            call fail(NP_STATS//' contains an invalid anisotropic-atom count')
        endif
        if( .not. ieee_is_finite(csv_diameter) .or. csv_diameter < 0.75 * real(MOLDIAM) .or. &
            &csv_diameter > 1.25 * real(MOLDIAM) )then
            call fail(NP_STATS//' diameter is inconsistent with the simulated nanoparticle')
        endif
    end subroutine check_statistics_files

    subroutine check_coordination_statistics( expected_atoms )
        integer, intent(in) :: expected_atoms
        integer             :: cn_counts(CN_MIN:CN_MAX)
        integer             :: funit, ios, i, cn, expected_cn_rows, actual_cn_rows
        real                :: atom_index, atom_nvox, atom_cn
        real                :: csv_cn, csv_natoms, csv_naniso
        logical             :: cn_seen(CN_MIN:CN_MAX)
        character(len=STDLEN) :: header
        cn_counts = 0
        cn_seen   = .false.
        if( .not. file_exists(ATOM_STATS) ) return
        ios = 0
        call fopen(funit, file=string(ATOM_STATS), status='old', action='read', iostat=ios)
        if( ios /= 0 )then
            call fail('could not open '//ATOM_STATS)
            return
        endif
        read(funit,'(a)',iostat=ios) header
        do i = 1, expected_atoms
            if( ios /= 0 ) exit
            read(funit,*,iostat=ios) atom_index, atom_nvox, atom_cn
            if( ios /= 0 ) exit
            if( .not. ieee_is_finite(atom_index) .or. .not. ieee_is_finite(atom_nvox) .or. &
                &.not. ieee_is_finite(atom_cn) )then
                call fail(ATOM_STATS//' contains non-finite identifying values')
                cycle
            endif
            if( nint(atom_index) /= i ) call fail(ATOM_STATS//' contains a non-sequential atom index')
            if( atom_nvox < 3.0 ) call fail(ATOM_STATS//' contains an undersized atom component')
            cn = nint(atom_cn)
            if( abs(atom_cn - real(cn)) > 0.001 )then
                call fail(ATOM_STATS//' contains a non-integral coordination number')
            elseif( cn >= CN_MIN .and. cn <= CN_MAX )then
                cn_counts(cn) = cn_counts(cn) + 1
            endif
        enddo
        call fclose(funit)
        if( ios /= 0 )then
            call fail(ATOM_STATS//' ended before all detected atoms were read')
            return
        endif
        expected_cn_rows = count(cn_counts >= 2)
        call check_line_count(CN_STATS, expected_cn_rows + 1)
        if( .not. file_exists(CN_STATS) ) return
        actual_cn_rows = nlines(string(CN_STATS)) - 1
        ios = 0
        call fopen(funit, file=string(CN_STATS), status='old', action='read', iostat=ios)
        if( ios /= 0 )then
            call fail('could not open '//CN_STATS)
            return
        endif
        read(funit,'(a)',iostat=ios) header
        do i = 1, actual_cn_rows
            if( ios /= 0 ) exit
            read(funit,*,iostat=ios) csv_cn, csv_natoms, csv_naniso
            if( ios /= 0 ) exit
            if( .not. ieee_is_finite(csv_cn) .or. .not. ieee_is_finite(csv_natoms) .or. &
                &.not. ieee_is_finite(csv_naniso) )then
                call fail(CN_STATS//' contains non-finite identifying values')
                cycle
            endif
            cn = nint(csv_cn)
            if( abs(csv_cn - real(cn)) > 0.001 .or. cn < CN_MIN .or. cn > CN_MAX )then
                call fail(CN_STATS//' contains an invalid coordination number')
                cycle
            endif
            if( cn_seen(cn) ) call fail(CN_STATS//' contains a duplicate coordination number')
            cn_seen(cn) = .true.
            if( cn_counts(cn) < 2 ) call fail(CN_STATS//' contains a group with fewer than two atoms')
            if( nint(csv_natoms) /= cn_counts(cn) ) call fail(CN_STATS//' atom count disagrees with atom statistics')
            if( csv_naniso < 0.0 .or. csv_naniso > csv_natoms )then
                call fail(CN_STATS//' contains an invalid anisotropic-atom count')
            endif
        enddo
        call fclose(funit)
        if( ios /= 0 ) call fail(CN_STATS//' ended before all coordination groups were read')
        do cn = CN_MIN, CN_MAX
            if( cn_counts(cn) >= 2 .and. .not. cn_seen(cn) )then
                call fail(CN_STATS//' is missing an eligible coordination group')
            endif
        enddo
    end subroutine check_coordination_statistics

    subroutine check_line_count( filename, expected_lines )
        character(len=*), intent(in) :: filename
        integer,          intent(in) :: expected_lines
        integer :: actual_lines
        if( .not. file_exists(filename) )then
            call fail(filename//' was not generated')
            return
        endif
        actual_lines = nlines(string(filename))
        write(logfhandle,'(a,a,a,i0)') '    ', filename, ' rows including header: ', actual_lines
        if( actual_lines /= expected_lines ) call fail(filename//' has an incorrect row count')
    end subroutine check_line_count
end subroutine exec_test_single_atoms_stats

subroutine exec_test_single_workflow( self, cline )
    use single_commanders_nano2D,       only: commander_analysis2D_nano
    use simple_commanders_sim,          only: commander_simulate_nanoparticle
    use simple_commanders_reproject,    only: commander_reproject
    use simple_commanders_stkops,       only: commander_stackops
    use simple_dock_vols,               only: dock_vols
    use simple_refine3D_fnames,         only: refine3D_state_vol_fbody
    use single_commanders_trajectory,   only: commander_trajectory_denoise
    use simple_commanders_project_ptcl, only: commander_import_particles
    use simple_commanders_project_core, only: commander_new_project
    use single_commanders_nano3D,       only: commander_autorefine3D_nano
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    class(commander_test_single_workflow), intent(inout) :: self
    class(cmdline),                        intent(inout) :: cline
    type(cmdline)                         :: cline_sim, cline_reproject, cline_trajectory, cline_denoise
    type(cmdline)                         :: cline_nproj, cline_imptcls, cline_an2Dnano, cline_aref3Dnano
    type(parameters)                      :: params
    type(commander_simulate_nanoparticle) :: xsim_nptcl
    type(commander_reproject)             :: xreproject
    type(commander_stackops)              :: xtrajectory
    type(commander_trajectory_denoise)    :: xdenoise
    type(commander_new_project)           :: xnproj
    type(commander_import_particles)      :: ximptcls
    type(commander_analysis2D_nano)       :: xan2Dnano
    type(commander_autorefine3D_nano)     :: xaref3Dnano
    type(string)                          :: projname, projfile, project_dir, startvol
    type(string)                          :: simulated_vol, reprojections, trajectory, denoised_trajectory
    type(string)                          :: autorefine_dir, final_volume
    character(len=*), parameter           :: REPROJECTIONS_STK = 'reprojections.mrc'
    character(len=*), parameter           :: TRAJECTORY_STK    = 'simulated_trajectory.mrc'
    character(len=*), parameter           :: DENOISED_STK      = 'denoised_trajectory.mrc'
    character(len=*), parameter           :: TRAJECTORY_ORITAB = 'glc_trajectory_oris.txt'
    character(len=*), parameter           :: SIMULATION_DIR    = '1_simulate_nanoparticle'
    character(len=*), parameter           :: REPROJECTION_DIR  = '2_generate_reprojections'
    character(len=*), parameter           :: TRAJECTORY_DIR    = '3_generate_trajectory'
    character(len=*), parameter           :: DENOISE_DIR       = '4_trajectory_denoise'
    character(len=*), parameter           :: IMPORT_DIR        = '5_import_particles'
    character(len=*), parameter           :: ANALYSIS2D_DIR    = '6_analysis2D_nano'
    integer,          parameter           :: NREPROJS = 5000, MASKDIAM = 40, NREFINE_ITERS = 5
    integer,          parameter           :: NFRAMES_PER_GROUP = 50
    integer                               :: chdir_status
    real,             parameter           :: TRAJECTORY_SNR    = 0.2
    real,             parameter           :: MIN_VOL_CORR      = 0.30
    real,             parameter           :: MAX_FSC0143       = 5.0
    real,             parameter           :: DOCK_HP           = 100.0
    real,             parameter           :: DOCK_LP           = 5.0
    real                               :: volume_corr, volume_fsc0143
    logical                            :: volume_ok
    write(logfhandle,'(a)') '>>> TEST_SINGLE_WORKFLOW:'
    if( .not. cline%defined('smpd') )    call cline%set('smpd', 0.358)
    if( .not. cline%defined('element') ) call cline%set('element', 'Pt')
    projname = 'test_single_workflow'
    call params%new(cline)
    projfile = projname%to_char()//'.simple'
    call cline_nproj%set('prg',                       'new_project')
    call cline_nproj%set('projname',             projname%to_char())
    call xnproj%execute(cline_nproj)
    call simple_getcwd(project_dir)
    projfile = filepath(project_dir, projfile)
    simulated_vol       = filepath(filepath(project_dir, SIMULATION_DIR), 'outvol.mrc')
    reprojections       = filepath(filepath(project_dir, REPROJECTION_DIR), REPROJECTIONS_STK)
    trajectory          = filepath(filepath(project_dir, TRAJECTORY_DIR), TRAJECTORY_STK)
    denoised_trajectory = filepath(filepath(project_dir, DENOISE_DIR), DENOISED_STK)
    startvol            = filepath(filepath(project_dir, ANALYSIS2D_DIR), 'startvol.mrc')
    
    call enter_workflow_stage(SIMULATION_DIR, projfile)
    call cline_sim%set('prg',               'simulate_nanoparticle')
    call cline_sim%set('box',                                   BOX)
    call cline_sim%set('smpd',                          params%smpd)
    call cline_sim%set('moldiam',                           MOLDIAM)
    call cline_sim%set('element',                    params%element)
    call cline_sim%set('nthr',                          params%nthr)
    call xsim_nptcl%execute(cline_sim)
    call return_to_project_dir

    call enter_workflow_stage(REPROJECTION_DIR, projfile)
    call make_glc_trajectory_oris(TRAJECTORY_ORITAB, NREPROJS, NFRAMES_PER_GROUP)
    call cline_reproject%set('prg',                     'reproject')
    call cline_reproject%set('pgrp',                           'c1')
    call cline_reproject%set('vol1',        simulated_vol%to_char())
    call cline_reproject%set('smpd',                    params%smpd)
    call cline_reproject%set('oritab',            TRAJECTORY_ORITAB)
    call cline_reproject%set('mskdiam',                          20)
    call cline_reproject%set('outstk',            REPROJECTIONS_STK)
    call cline_reproject%set('nthr',                    params%nthr)
    call xreproject%execute(cline_reproject)
    call return_to_project_dir

    call enter_workflow_stage(TRAJECTORY_DIR, projfile)
    call cline_trajectory%set('prg',                     'stackops')
    call cline_trajectory%set('mkdir',                         'no')
    call cline_trajectory%set('stk',         reprojections%to_char())
    call cline_trajectory%set('outstk',              TRAJECTORY_STK)
    call cline_trajectory%set('smpd',                   params%smpd)
    call cline_trajectory%set('snr',                 TRAJECTORY_SNR)
    call cline_trajectory%set('nthr',                   params%nthr)
    call xtrajectory%execute(cline_trajectory)
    call return_to_project_dir

    call enter_workflow_stage(DENOISE_DIR, projfile)
    call cline_denoise%set('prg',              'trajectory_denoise')
    call cline_denoise%set('mkdir',                            'no')
    call cline_denoise%set('stk',              trajectory%to_char())
    call cline_denoise%set('outstk',                   DENOISED_STK)
    call cline_denoise%set('smpd',                      params%smpd)
    call cline_denoise%set('nthr',                      params%nthr)
    call xdenoise%execute(cline_denoise)
    call return_to_project_dir

    call enter_workflow_stage(IMPORT_DIR, projfile)
    call cline_imptcls%set('prg',                'import_particles')
    call cline_imptcls%set('mkdir',                            'no')
    call cline_imptcls%set('projfile',           projfile%to_char())
    call cline_imptcls%set('stk',     denoised_trajectory%to_char())
    call cline_imptcls%set('smpd',                      params%smpd)
    call cline_imptcls%set('ctf',                              'no')
    call ximptcls%execute(cline_imptcls)
    call return_to_project_dir

    call enter_workflow_stage(ANALYSIS2D_DIR, projfile)
    call cline_an2Dnano%set('prg',                'analysis2D_nano')
    call cline_an2Dnano%set('mkdir',                           'no')
    call cline_an2Dnano%set('projfile',          projfile%to_char())
    call cline_an2Dnano%set('element',               params%element)
    call cline_an2Dnano%set('nthr',                     params%nthr)
    call xan2Dnano%execute(cline_an2Dnano)
    if( .not. file_exists(startvol) ) THROW_HARD('analysis2D_nano did not generate '//startvol%to_char())
    call return_to_project_dir

    call cline_aref3Dnano%set('prg',            'autorefine3D_nano')
    call cline_aref3Dnano%set('projfile',        projfile%to_char())
    call cline_aref3Dnano%set('vol1',            startvol%to_char())
    call cline_aref3Dnano%set('smpd',                   params%smpd)
    call cline_aref3Dnano%set('element',             params%element)
    call cline_aref3Dnano%set('nthr',                   params%nthr)
    call cline_aref3Dnano%set('pgrp',                          'c1')
    call cline_aref3Dnano%set('lp',                             1.5)  
    call cline_aref3Dnano%set('mskdiam',                   MASKDIAM)
    call cline_aref3Dnano%set('maxits',               NREFINE_ITERS)
    call xaref3Dnano%execute(cline_aref3Dnano)
    call simple_getcwd(autorefine_dir)
    final_volume = filepath(filepath(autorefine_dir, 'final_results'), &
        &refine3D_state_vol_fbody(1)//'_iter'//int2str_pad(NREFINE_ITERS, 3)//MRC_EXT)
    call validate_reconstructed_volume(simulated_vol, final_volume, params%smpd, BOX, real(MOLDIAM), &
        &MIN_VOL_CORR, MAX_FSC0143, volume_corr, volume_fsc0143, volume_ok)
    call return_to_project_dir
    if( .not. volume_ok ) THROW_HARD('TEST_SINGLE_WORKFLOW FAILED: final-volume validation failed')
    write(logfhandle,'(a,f7.4,a,f7.2,a)') 'PASS: single_workflow whole-volume correlation=', volume_corr, &
        &', FSC=0.143 at ', volume_fsc0143, ' A'
    call simple_end('**** SIMPLE_TEST_SINGLE_WORKFLOW NORMAL STOP ****')

contains
    subroutine enter_workflow_stage( stage, stage_projfile )
        character(len=*), intent(in)    :: stage
        type(string),     intent(inout) :: stage_projfile
        type(string)                    :: stage_dir, previous_projfile
        stage_dir         = filepath(project_dir, stage)
        previous_projfile = stage_projfile
        stage_projfile    = filepath(stage_dir, basename(previous_projfile))
        call simple_mkdir(stage_dir)
        call simple_copy_file(previous_projfile, stage_projfile)
        call simple_chdir(stage_dir, chdir_status)
        if( chdir_status /= 0 ) THROW_HARD('Could not enter single_workflow stage')
    end subroutine enter_workflow_stage

    subroutine return_to_project_dir
        call simple_chdir(project_dir, chdir_status)
        if( chdir_status /= 0 ) THROW_HARD('Could not return to the single_workflow project directory')
    end subroutine return_to_project_dir

    subroutine validate_reconstructed_volume( truth_fname, reconstruction_fname, expected_smpd, expected_box, &
        &mask_diameter, min_corr, max_fsc0143, corr, fsc0143, passed )
        type(string), intent(in) :: truth_fname, reconstruction_fname
        real,         intent(in) :: expected_smpd, mask_diameter, min_corr, max_fsc0143
        integer,      intent(in) :: expected_box
        real,         intent(out) :: corr, fsc0143
        logical,      intent(out) :: passed
        character(len=*), parameter :: TRUTH_COMPARE = 'workflow_truth_compare.mrc'
        character(len=*), parameter :: RECON_MIRROR  = 'workflow_reconstruction_mirror.mrc'
        character(len=*), parameter :: DOCKED_DIRECT = 'workflow_reconstruction_docked.mrc'
        character(len=*), parameter :: DOCKED_MIRROR = 'workflow_reconstruction_mirror_docked.mrc'
        type(dock_vols) :: docker
        type(image)     :: truth, reconstruction
        type(string)    :: selected_reconstruction
        real, allocatable :: fsc(:), resolutions(:)
        integer :: truth_ldim(3), reconstruction_ldim(3), nsections, nyq
        real    :: truth_smpd, reconstruction_smpd, direct_cc, mirror_cc
        real    :: eulers(3), shifts(3), fsc05, mask_radius

        passed  = .false.
        corr    = 0.0
        fsc0143 = 0.0
        if( .not. file_exists(truth_fname) )then
            write(logfhandle,'(a)') '    FAIL: simulated truth volume was not generated'
            return
        endif
        if( .not. file_exists(reconstruction_fname) )then
            write(logfhandle,'(a,a)') '    FAIL: final reconstruction was not generated: ', &
                &reconstruction_fname%to_char()
            return
        endif
        call find_ldim_nptcls(truth_fname, truth_ldim, nsections)
        truth_smpd = find_img_smpd(truth_fname)
        call find_ldim_nptcls(reconstruction_fname, reconstruction_ldim, nsections)
        reconstruction_smpd = find_img_smpd(reconstruction_fname)
        write(logfhandle,'(a,3(i0,1x),a,f7.3)') '>>> Simulated truth dimensions/sampling: ', truth_ldim, &
            &' / ', truth_smpd
        write(logfhandle,'(a,3(i0,1x),a,f7.3)') '>>> Final volume dimensions/sampling:    ', reconstruction_ldim, &
            &' / ', reconstruction_smpd
        if( any(reconstruction_ldim /= [expected_box, expected_box, expected_box]) )then
            write(logfhandle,'(a,i0)') '    FAIL: final volume does not have the expected cubic box ', expected_box
            return
        endif
        if( abs(reconstruction_smpd - expected_smpd) > 0.001 )then
            write(logfhandle,'(a,f7.3)') '    FAIL: final volume has incorrect sampling; expected ', expected_smpd
            return
        endif
        if( any(truth_ldim /= reconstruction_ldim) )then
            write(logfhandle,'(a)') '    FAIL: simulated truth and final volume dimensions do not match'
            return
        endif
        if( abs(truth_smpd - reconstruction_smpd) > 0.001 )then
            write(logfhandle,'(a)') '    FAIL: simulated truth and final volume sampling do not match'
            return
        endif

        call truth%new(truth_ldim, truth_smpd, wthreads=.false.)
        call truth%read(truth_fname)
        call truth%write(string(TRUTH_COMPARE))
        call reconstruction%new(reconstruction_ldim, reconstruction_smpd, wthreads=.false.)
        call reconstruction%read(reconstruction_fname)
        call reconstruction%mirror('x')
        call reconstruction%write(string(RECON_MIRROR))
        call reconstruction%kill

        call docker%new(string(TRUTH_COMPARE), reconstruction_fname, reconstruction_smpd, &
            &DOCK_HP, DOCK_LP, mask_diameter)
        call docker%srch()
        call docker%get_dock_info(eulers, shifts, direct_cc)
        call docker%rotate_target(reconstruction_fname, string(DOCKED_DIRECT))
        call docker%kill()
        call docker%new(string(TRUTH_COMPARE), string(RECON_MIRROR), reconstruction_smpd, &
            &DOCK_HP, DOCK_LP, mask_diameter)
        call docker%srch()
        call docker%get_dock_info(eulers, shifts, mirror_cc)
        call docker%rotate_target(string(RECON_MIRROR), string(DOCKED_MIRROR))
        call docker%kill()
        if( direct_cc >= mirror_cc )then
            selected_reconstruction = DOCKED_DIRECT
        else
            selected_reconstruction = DOCKED_MIRROR
        endif
        write(logfhandle,'(a,f7.4,a,f7.4)') '>>> Docking correlation: direct=', direct_cc, ', mirrored=', mirror_cc

        call reconstruction%new(reconstruction_ldim, reconstruction_smpd, wthreads=.false.)
        call reconstruction%read(selected_reconstruction)
        corr = truth%real_corr(reconstruction)
        write(logfhandle,'(a,f7.4,a,f7.4)') '>>> Registered whole-volume Pearson correlation: ', corr, &
            &'; minimum ', min_corr
        if( .not. ieee_is_finite(corr) .or. corr < min_corr )then
            write(logfhandle,'(a)') '    FAIL: final-volume Pearson correlation is below the required minimum'
            call truth%kill
            call reconstruction%kill
            return
        endif

        mask_radius = 0.5 * mask_diameter / reconstruction_smpd
        call truth%mask3D_soft(mask_radius, backgr=0.0)
        call reconstruction%mask3D_soft(mask_radius, backgr=0.0)
        call truth%fft()
        call reconstruction%fft()
        nyq = truth%get_filtsz()
        allocate(fsc(nyq), source=0.0)
        call truth%fsc(reconstruction, fsc)
        resolutions = truth%get_res()
        if( any(.not. ieee_is_finite(fsc)) )then
            write(logfhandle,'(a)') '    FAIL: final-volume FSC contains non-finite values'
            call truth%kill
            call reconstruction%kill
            deallocate(fsc, resolutions)
            return
        endif
        call get_resolution(fsc, resolutions, fsc05, fsc0143)
        if( fsc05 > 0.0 )   fsc05   = max(fsc05,   2.0 * reconstruction_smpd)
        if( fsc0143 > 0.0 ) fsc0143 = max(fsc0143, 2.0 * reconstruction_smpd)
        write(logfhandle,'(a,f7.2,a,f7.2,a,f7.2,a)') '>>> Masked truth FSC: 0.500 at ', fsc05, &
            &' A; 0.143 at ', fsc0143, ' A; maximum ', max_fsc0143, ' A'
        passed = ieee_is_finite(fsc0143) .and. fsc0143 > 0.0 .and. fsc0143 <= max_fsc0143
        if( .not. passed ) write(logfhandle,'(a)') '    FAIL: final-volume FSC resolution is outside the accepted range'
        call truth%kill
        call reconstruction%kill
        deallocate(fsc, resolutions)
    end subroutine validate_reconstructed_volume
end subroutine exec_test_single_workflow

subroutine make_glc_trajectory_oris( oritab, nreprojs, nframes_per_group )
    character(len=*), intent(in) :: oritab
    integer,          intent(in) :: nreprojs, nframes_per_group
    real, parameter              :: ANGULAR_SPAN = 2.0
    type(oris)                   :: group_oris, trajectory_oris
    real                         :: base_euls(3), euls(3), frame_frac
    integer                      :: iframe, igroup, iproj, ngroups
    if( nframes_per_group < 2 ) THROW_HARD('GLC orientation groups require at least two frames')
    if( mod(nreprojs, nframes_per_group) /= 0 ) THROW_HARD('GLC reprojection count must divide into equal frame groups')
    ngroups = nreprojs / nframes_per_group
    call group_oris%new(ngroups, is_ptcl=.false.)
    call group_oris%spiral()
    call trajectory_oris%new(nreprojs, is_ptcl=.false.)
    do igroup = 1, ngroups
        base_euls = group_oris%get_euler(igroup)
        do iframe = 1, nframes_per_group
            iproj      = (igroup - 1) * nframes_per_group + iframe
            frame_frac = real(iframe - 1) / real(nframes_per_group - 1)
            euls       = base_euls + ANGULAR_SPAN * (frame_frac - 0.5) * [1.0, 0.5, 1.0]
            call trajectory_oris%set_euler(iproj, euls)
        enddo
    enddo
    call trajectory_oris%set_all2single('state', 1.0)
    call trajectory_oris%write(string(oritab), [1, nreprojs])
    call trajectory_oris%kill
    call group_oris%kill
end subroutine make_glc_trajectory_oris

end module simple_commanders_test_single
