!@descr: SINGLE test commanders: the nanoparticle atoms pipeline, species discovery and the SINGLE workflow
module simple_commanders_test_single
use simple_commanders_api
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_single_atoms_stats
  contains
    procedure :: execute      => exec_test_single_atoms_stats
end type commander_test_single_atoms_stats

type, extends(commander_base) :: commander_test_species_discovery
  contains
    procedure :: execute      => exec_test_species_discovery
end type commander_test_species_discovery

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

subroutine exec_test_species_discovery( self, cline )
    use simple_commanders_atoms, only: commander_detect_atoms
    use simple_commanders_sim,   only: commander_simulate_nanoparticle
    use simple_atoms,            only: atoms
    use simple_test_utils,       only: set_fixed_seed
    use simple_rnd,              only: gasdev, ran3
    class(commander_test_species_discovery), intent(inout) :: self
    class(cmdline),                                 intent(inout) :: cline
    ! floors of section 10.4 of doc/implementation_notes/planned/species_discovery.md, set before the first run
    real,    parameter :: MIN_WEAK_RECALL  = 0.90   ! weak class, over its atoms of predicted SNR >= SNR_RECALL
    real,    parameter :: SNR_RECALL       = 6.5
    integer, parameter :: MAX_FALSE        = 1      ! false atoms per particle
    real,    parameter :: RATIO_TOL        = 0.10   ! class intensity ratio, relative
    real,    parameter :: WIDTH_TOL_STRONG = 0.03   ! per-shell width of the strongest class, relative
    real,    parameter :: WIDTH_TOL_LIGHT  = 0.10   ! per-shell width of the light class, relative
    real,    parameter :: MATCH_NN_FRAC    = 0.3    ! found atom matches a generating atom within this many d_NN
    real,    parameter :: SIGMA_REF        = 0.4196 ! A, width of the detection template (B 13.9 A**2)
    ! generating models of section 8 on the Pt lattice of single_atoms_stats
    real,    parameter :: SMPD_FIX      = 0.358
    real,    parameter :: SIGMA_CORE    = 0.35      ! A, strongest class in the core
    real,    parameter :: SIGMA_SCATTER = 0.05      ! relative per-atom scatter of sigma
    real,    parameter :: Q_LIGHT       = 1. / 6.   ! intensity of the light class
    real,    parameter :: LIGHT_FRAC    = 0.25
    real,    parameter :: RADIAL_VAR    = 2.        ! radial cases: variance at the surface over the centre
    real,    parameter :: LIGHT_VAR_R6  = 1.3       ! R6: light over heavy variance at every radius
    integer, parameter :: NSHELL_MAX    = 5
    integer, parameter :: NSHELL_ATOMS  = 20
    integer, parameter :: NCASES        = 4
    character(len=5), parameter :: CASES(NCASES)       = ['case2', 'case9', 'R3   ', 'R6   ']
    real,             parameter :: NOISE(NCASES)       = [0.03, 0.05, 0.05, 0.03] ! sdev over the peak of a core atom
    logical,          parameter :: TWO_SPECIES(NCASES) = [.true., .false., .false., .true.]
    logical,          parameter :: RADIAL(NCASES)      = [.false., .false., .true., .true.]
    character(len=*), parameter :: PRODUCTS(5) = [character(len=12) :: 'map_ATMS.pdb', 'map_BIN.mrc', 'map_CC.mrc',&
        &'map_MSK.mrc', 'map_SIM.mrc']
    character(len=*), parameter :: RUNS(3) = [character(len=6) :: 'plain', 'disc', 'halves']
    type(commander_simulate_nanoparticle) :: xsim
    type(commander_detect_atoms)          :: xdet
    type(cmdline)                         :: cline_sim, cline_det
    type(atoms)                           :: lattice, model
    type(string)                          :: cwd_saved, fixture_root, case_dir
    real,    allocatable :: lat(:,:), gxyz(:,:), gsig(:), gq(:)
    integer, allocatable :: gcls(:)
    real    :: dnn_gen, cen(3), rmax
    integer :: status, nthr, nlat, icase, irun, i
    logical :: all_ok
    write(logfhandle,'(a)') '>>> TEST_SPECIES_DISCOVERY:'
    call set_fixed_seed(20261008)
    nthr = 8
    if( cline%defined('nthr') ) nthr = cline%get_iarg('nthr')
    all_ok = .true.
    call simple_getcwd(cwd_saved)
    fixture_root = filepath(cwd_saved, 'test_species_discovery_'//int2str(get_process_id()))
    if( dir_exists(fixture_root) ) call simple_rmdir(fixture_root)
    call simple_mkdir(fixture_root)
    call simple_chdir(fixture_root, status)
    if( status /= 0 ) THROW_HARD('TEST_SPECIES_DISCOVERY FAILED: could not enter fixture directory')
    ! the Pt lattice of phase 0: positions only
    call cline_sim%set('prg',     'simulate_nanoparticle')
    call cline_sim%set('element', 'Pt')
    call cline_sim%set('moldiam', real(MOLDIAM))
    call cline_sim%set('box',     BOX)
    call cline_sim%set('smpd',    SMPD_FIX)
    call cline_sim%set('nthr',    nthr)
    call cline_sim%set('outvol',  'lattice.mrc')
    call cline_sim%set('pdbout',  'lattice.pdb')
    call xsim%execute(cline_sim)
    call cline_sim%kill
    call del_file('lattice.mrc')
    call lattice%new(string('lattice.pdb'))
    nlat = lattice%get_n()
    allocate(lat(3,nlat))
    do i = 1,nlat
        lat(:,i) = lattice%get_coord(i)
    enddo
    call lattice%kill
    dnn_gen = median_nn(lat)
    cen     = sum(lat, dim=2) / real(nlat)
    rmax    = maxval(sqrt(sum((lat - spread(cen, 2, nlat))**2, dim=1)))
    write(logfhandle,'(a,i0,a,f7.4,a,f7.3,a)') '    lattice: ', nlat, ' atoms, d_NN ', dnn_gen, ' A, radius ', rmax, ' A'
    do icase = 1,NCASES
        case_dir = filepath(fixture_root, trim(CASES(icase)))
        call simple_mkdir(case_dir)
        call simple_chdir(case_dir, status)
        if( status /= 0 ) THROW_HARD('TEST_SPECIES_DISCOVERY FAILED: could not enter a case directory')
        write(logfhandle,'(a)') '>>> CASE '//trim(CASES(icase))
        call make_fixture(icase)
        do irun = 1,3
            call run_detect(irun)
        enddo
        call compare_products
        call evaluate(2)
        call evaluate(3)
        deallocate(gcls, gxyz, gsig, gq)
        call simple_chdir(fixture_root, status)
    enddo
    call simple_chdir(cwd_saved, status)
    if( status /= 0 ) THROW_HARD('TEST_SPECIES_DISCOVERY FAILED: could not restore original directory')
    if( all_ok )then
        write(logfhandle,'(a)') 'PASS: species_discovery'
        call simple_end('**** SIMPLE_TEST_SPECIES_DISCOVERY NORMAL STOP ****')
    else
        THROW_HARD('TEST_SPECIES_DISCOVERY FAILED')
    endif

  contains

    subroutine fail( message )
        character(len=*), intent(in) :: message
        write(logfhandle,'(a)') '    FAIL: '//trim(message)
        all_ok = .false.
    end subroutine fail

    ! median nearest-neighbour distance, the upper middle value for an even count
    real function median_nn( xyz )
        real, intent(in) :: xyz(:,:)
        real    :: nn(size(xyz,2))
        integer :: j, l
        nn = huge(1.)
        do j = 1,size(xyz,2)
            do l = 1,size(xyz,2)
                if( l /= j ) nn(j) = min(nn(j), norm2(xyz(:,j) - xyz(:,l)))
            enddo
        enddo
        call hpsort(nn)
        median_nn = nn(size(nn) / 2 + 1)
    end function median_nn

    ! pseudo-atom model of the case on the lattice, rendered, with noise; half maps whose average is the map
    subroutine make_fixture( ic )
        integer, intent(in) :: ic
        type(image)       :: clean, nimg, even, odd
        real, allocatable :: keys(:)
        integer, allocatable :: order(:)
        real    :: sig, r, bpeak, sdev_noise
        integer :: j, nlight
        allocate(gcls(nlat), source=1)
        if( TWO_SPECIES(ic) )then
            ! a random quarter of the atoms is light
            allocate(keys(nlat))
            do j = 1,nlat
                keys(j) = ran3()
            enddo
            order = [(j, j=1,nlat)]
            call hpsort(keys, order)
            nlight = nint(LIGHT_FRAC * real(nlat))
            gcls(order(:nlight)) = 2
        endif
        call model%new(nlat, dummy=.true.)
        do j = 1,nlat
            sig = SIGMA_CORE * (1. + SIGMA_SCATTER * gasdev(0., 1.))
            if( RADIAL(ic) )then
                r   = norm2(lat(:,j) - cen)
                sig = sig * sqrt(1. + (RADIAL_VAR - 1.) * (r / rmax)**2)
            endif
            if( trim(CASES(ic)) == 'R6' .and. gcls(j) == 2 ) sig = sig * sqrt(LIGHT_VAR_R6)
            if( gcls(j) == 1 )then
                call model%set_element(j, 'X1')
                call model%set_occupancy(j, 1.)
            else
                call model%set_element(j, 'X2')
                call model%set_occupancy(j, Q_LIGHT)
            endif
            call model%set_coord(j, lat(:,j))
            call model%set_beta(j, 8. * PI**2 * sig**2)
            call model%set_num(j, j)
            call model%set_resnum(j, j)
        enddo
        call model%writepdb(string('model.pdb'))
        call model%kill
        ! the generating model as written: two decimals of occupancy and B
        call model%new(string('model.pdb'))
        allocate(gxyz(3,nlat), gsig(nlat), gq(nlat))
        do j = 1,nlat
            gxyz(:,j) = model%get_coord(j)
            gq(j)     = model%get_occupancy(j)
            gsig(j)   = sqrt(model%get_beta(j) / (8. * PI**2))
        enddo
        call model%kill
        call cline_sim%set('prg',     'simulate_nanoparticle')
        call cline_sim%set('pdbfile', 'model.pdb')
        call cline_sim%set('box',     BOX)
        call cline_sim%set('smpd',    SMPD_FIX)
        call cline_sim%set('nthr',    nthr)
        call cline_sim%set('outvol',  'clean.mrc')
        call cline_sim%set('pdbout',  'model_out.pdb')
        call xsim%execute(cline_sim)
        call cline_sim%kill
        ! noise relative to the peak of a core atom of the strongest class
        bpeak      = 8. * PI**2 * SIGMA_CORE**2
        sdev_noise = NOISE(ic) * (4. * PI / bpeak)**1.5
        write(logfhandle,'(a,i0,a,i0,a,es12.4)') '    model: ', count(gcls == 1), ' strong, ', count(gcls == 2),&
            &' light atoms; noise sdev ', sdev_noise
        call clean%new([BOX,BOX,BOX], SMPD_FIX)
        call clean%read(string('clean.mrc'))
        call nimg%new([BOX,BOX,BOX], SMPD_FIX)
        call even%copy(clean)
        call nimg%gauran(0., sqrt(2.) * sdev_noise)
        call even%add(nimg)
        call odd%copy(clean)
        call nimg%gauran(0., sqrt(2.) * sdev_noise)
        call odd%add(nimg)
        call even%write(string('even.mrc'))
        call odd%write(string('odd.mrc'))
        call even%add(odd)
        call even%div(2.)
        call even%write(string('map.mrc'))
        call clean%kill
        call nimg%kill
        call even%kill
        call odd%kill
        call del_file('clean.mrc')
    end subroutine make_fixture

    ! detect_atoms without element in a directory of its own: plain, with discovery, with discovery and half maps.
    ! Full threads: the product comparison between the three runs also checks that detection is thread-invariant.
    subroutine run_detect( ir )
        integer, intent(in) :: ir
        integer :: st
        call simple_mkdir(trim(RUNS(ir)))
        call simple_chdir(trim(RUNS(ir)), st)
        call cline_det%set('prg',  'detect_atoms')
        call cline_det%set('vol1', filepath(case_dir, 'map.mrc'))
        call cline_det%set('smpd', SMPD_FIX)
        call cline_det%set('nthr', nthr)
        if( ir >= 2 ) call cline_det%set('discover_species', 'yes')
        if( ir == 3 )then
            call cline_det%set('vol_even', filepath(case_dir, 'even.mrc'))
            call cline_det%set('vol_odd',  filepath(case_dir, 'odd.mrc'))
        endif
        call xdet%execute(cline_det)
        call cline_det%kill
        call simple_chdir(case_dir, st)
    end subroutine run_detect

    ! the present products do not depend on discover_species
    subroutine compare_products
        integer :: ip, ir
        do ir = 2,3
            do ip = 1,size(PRODUCTS)
                if( .not. files_identical(string(trim(RUNS(1))//'/'//trim(PRODUCTS(ip))), string(trim(RUNS(ir))//'/'//trim(PRODUCTS(ip)))) )then
                    call fail(trim(CASES(icase))//': '//trim(PRODUCTS(ip))//' differs with '//trim(RUNS(ir)))
                endif
            enddo
        enddo
        write(logfhandle,'(a)') '    present products compared with and without discover_species'
    end subroutine compare_products

    ! value of key in the key = value report
    real function report_value( fname, key )
        character(len=*), intent(in) :: fname, key
        character(len=256) :: line
        integer :: u, ios, ieq
        report_value = -huge(1.)
        open(newunit=u, file=fname, status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        do
            read(u,'(a)',iostat=ios) line
            if( ios /= 0 ) exit
            ieq = index(line, '=')
            if( ieq < 2 ) cycle
            if( trim(adjustl(line(:ieq-1))) /= key ) cycle
            read(line(ieq+1:),*,iostat=ios) report_value
            exit
        enddo
        close(u)
    end function report_value

    ! text value of key in the key = value report
    function report_string( fname, key ) result( val )
        character(len=*), intent(in) :: fname, key
        character(len=:), allocatable :: val
        character(len=256) :: line
        integer :: u, ios, ieq
        val = ''
        open(newunit=u, file=fname, status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        do
            read(u,'(a)',iostat=ios) line
            if( ios /= 0 ) exit
            ieq = index(line, '=')
            if( ieq < 2 ) cycle
            if( trim(adjustl(line(:ieq-1))) /= key ) cycle
            val = trim(adjustl(line(ieq+1:)))
            exit
        enddo
        close(u)
    end function report_string

    ! a discovery run against the generating model
    subroutine evaluate( ir )
        integer, intent(in) :: ir
        integer, parameter :: NCOL = 21
        character(len=512) :: line
        character(len=:), allocatable :: tag, csv, txt, pdb
        real,    allocatable :: rows(:,:), dist(:), sig_fit(:), sig_gen(:), rad(:)
        integer, allocatable :: match(:), hits(:), order(:), sel(:)
        type(atoms) :: species_atms
        real    :: tol, s_a, ratio, ratio_gen, snr, recall, rms_fit, rms_gen, rel
        integer :: u, ios, nrow, j, l, nfalse, nwrong, nrec, k_found, nweak, nweak_found, kc, nsh, ish, m, lo, hi
        tag = trim(CASES(icase))//'/'//trim(RUNS(ir))
        csv = trim(RUNS(ir))//'/map_species.csv'
        txt = trim(RUNS(ir))//'/map_species.txt'
        if( .not. file_exists(csv) .or. .not. file_exists(txt) )then
            call fail(tag//': species files were not written')
            return
        endif
        nrow = nlines(string(csv)) - 1
        allocate(rows(NCOL,nrow))
        open(newunit=u, file=csv, status='old', action='read')
        read(u,'(a)') line
        do j = 1,nrow
            read(u,'(a)') line
            do l = 1,len_trim(line)
                if( line(l:l) == ',' ) line(l:l) = ' '
            enddo
            read(line,*,iostat=ios) rows(:,j)
            if( ios /= 0 ) call fail(tag//': unreadable species table row')
        enddo
        close(u)
        ! the species PDB carries one atom per table row, with the class in the element column
        pdb = trim(RUNS(ir))//'/map_species.pdb'
        if( .not. file_exists(pdb) )then
            call fail(tag//': species PDB was not written')
        else
            call species_atms%new(string(pdb))
            if( species_atms%get_n() /= nrow )then
                call fail(tag//': species PDB and table disagree on the atom count')
            else
                do j = 1,nrow
                    if( species_atms%get_element(j) /= 'X'//int2str(nint(rows(16,j))) )then
                        call fail(tag//': species PDB element column does not match the class of the table')
                        exit
                    endif
                enddo
            endif
            call species_atms%kill
        endif
        k_found = nint(report_value(txt, 'K'))
        s_a     = report_value(txt, 's_A')
        if( ir == 3 .neqv. report_string(txt, 'noise_source') == 'half_maps' )then
            call fail(tag//': the noise source is not the one given')
        endif
        ! match every found atom to its nearest generating atom
        tol = MATCH_NN_FRAC * dnn_gen
        allocate(match(nrow), hits(nlat), source=0)
        allocate(dist(nrow), source=0.)
        do j = 1,nrow
            match(j) = minloc(sum((gxyz - spread(rows(3:5,j), 2, nlat))**2, dim=1), dim=1)
            dist(j)  = norm2(gxyz(:,match(j)) - rows(3:5,j))
            if( dist(j) > tol )then
                match(j) = 0
            else
                hits(match(j)) = hits(match(j)) + 1
            endif
        enddo
        ! unmatched atoms and second atoms on one site are false
        nfalse = count(match == 0) + sum(max(hits - 1, 0))
        nrec   = count(nint(rows(2,:)) == 1)
        write(logfhandle,'(3a,i0,a,i0,a,i0,a,i0,a,es11.4)') '    ', tag, ': found ', nrow, ' (recovered ', nrec,&
            &'), K = ', k_found, ', false ', nfalse, ', s_A ', s_a
        if( nfalse > MAX_FALSE ) call fail(tag//': more than one false atom')
        if( .not. any(gcls == 2) )then
            if( k_found /= 1 ) call fail(tag//': a single species gave K /= 1')
            if( nrec /= 0 )    call fail(tag//': a single species gave recovered atoms')
            return
        endif
        if( k_found /= 2 ) call fail(tag//': two species gave K /= 2')
        ! labels of the found atoms
        nwrong = 0
        do j = 1,nrow
            if( match(j) == 0 ) cycle
            if( nint(rows(16,j)) /= gcls(match(j)) ) nwrong = nwrong + 1
        enddo
        write(logfhandle,'(a,i0)') '        wrong labels: ', nwrong
        if( nwrong > 0 ) call fail(tag//': a found atom has the wrong label')
        ! recall of the weak class where its predicted signal-to-noise is SNR_RECALL or more
        nweak       = 0
        nweak_found = 0
        do j = 1,nlat
            if( gcls(j) /= 2 ) cycle
            ! section 7: the amplitude over s_A, scaled from the template width as (sigma / sigma_ref)**1.5
            snr = gq(j) / (2. * PI * gsig(j)**2)**1.5 * (gsig(j) / SIGMA_REF)**1.5 / s_a
            if( snr < SNR_RECALL ) cycle
            nweak = nweak + 1
            if( hits(j) > 0 ) nweak_found = nweak_found + 1
        enddo
        recall = real(nweak_found) / real(max(nweak, 1))
        write(logfhandle,'(a,i0,a,i0,a,f7.4,a,i0)') '        weak recall at predicted SNR >= 6.5: ', nweak_found, ' of ', nweak,&
            &' = ', recall, '; weak atoms in total ', count(gcls == 2)
        if( nweak == 0 ) call fail(tag//': no weak atom reaches the predicted SNR of the recall floor')
        if( recall < MIN_WEAK_RECALL ) call fail(tag//': weak-class recall below 0.90')
        ! intensity ratio
        if( k_found == 2 )then
            ratio     = report_value(txt, 'class_ratio_2')
            ratio_gen = sum(gq, mask=gcls == 2) / real(count(gcls == 2)) / (sum(gq, mask=gcls == 1) / real(count(gcls == 1)))
            write(logfhandle,'(a,f8.4,a,f8.4)') '        intensity ratio ', ratio, ', generating ', ratio_gen
            if( abs(ratio / ratio_gen - 1.) > RATIO_TOL ) call fail(tag//': intensity ratio off by more than 10%')
        endif
        ! per-shell widths of the right-labelled found atoms against their generating widths
        do kc = 1,2
            sel = pack([(j, j=1,nrow)], match > 0)
            sel = pack(sel, nint(rows(16,sel)) == kc .and. gcls(match(sel)) == kc)
            m   = size(sel)
            if( m == 0 ) cycle
            rad   = rows(6,sel)
            order = [(j, j=1,m)]
            call hpsort(rad, order)
            sig_fit = sqrt(rows(13,sel(order)) / (8. * PI**2))
            sig_gen = gsig(match(sel(order)))
            nsh = max(1, min(NSHELL_MAX, m / NSHELL_ATOMS))
            do ish = 1,nsh
                lo = ((ish-1) * m) / nsh + 1
                hi = (ish * m) / nsh
                rms_fit = sqrt(sum(sig_fit(lo:hi)**2) / real(hi - lo + 1))
                rms_gen = sqrt(sum(sig_gen(lo:hi)**2) / real(hi - lo + 1))
                rel     = rms_fit / rms_gen - 1.
                write(logfhandle,'(a,i0,a,i0,a,i0,a,f7.4,a,f7.4,a,f7.4)') '        class ', kc, ' shell ', ish, ' (', hi - lo + 1,&
                    &' atoms): sigma ', rms_fit, ', generating ', rms_gen, ', relative ', rel
                if( kc == 1 .and. abs(rel) > WIDTH_TOL_STRONG ) call fail(tag//': strong-class shell width off by more than 3%')
                if( kc == 2 .and. abs(rel) > WIDTH_TOL_LIGHT )  call fail(tag//': light-class shell width off by more than 10%')
            enddo
        enddo
    end subroutine evaluate

end subroutine exec_test_species_discovery

subroutine exec_test_single_workflow( self, cline )
    use single_commanders_nano2D,       only: commander_analysis2D_nano
    use simple_commanders_sim,          only: commander_simulate_nanoparticle
    use simple_commanders_reproject,    only: commander_reproject
    use simple_commanders_stkops,       only: commander_stackops
    use simple_test_truth_metrics,      only: validate_reconstructed_volume
    use simple_refine3D_fnames,         only: refine3D_state_vol_fbody
    use single_commanders_trajectory,   only: commander_trajectory_denoise
    use simple_commanders_project_ptcl, only: commander_import_particles
    use simple_commanders_project_core, only: commander_new_project
    use single_commanders_nano3D,       only: commander_autorefine3D_nano, commander_refine3D_nano
    use simple_test_truth_metrics,      only: pair_pose_error
    use simple_test_utils,              only: set_fixed_seed
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
    type(commander_refine3D_nano)         :: xref3Dnano
    type(cmdline)                         :: cline_cont
    type(sp_project)                      :: run_proj
    type(oris)                            :: truth_oris
    type(string)                          :: cont_projfile
    integer, allocatable                  :: pinds(:)
    real                                  :: pair_polar, pair_cont, frac5, pair_nopolish, pair_polish
    type(string)                          :: aref_projfile, nano_projfile
    integer                               :: i, npairs_ptcls
    type(string)                          :: suite_name, projname, projfile, project_dir, startvol
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
    character(len=*), parameter           :: CONT_DIR          = '8_refine3D_nano_cont'
    character(len=*), parameter           :: NOPOLISH_DIR      = '9_refine3D_nano_nopolish'
    character(len=*), parameter           :: POLISH_DIR        = '10_refine3D_nano_polish'
    ! N29, Phase 9 (declared before the first run): refine3D_nano pose_cont=yes and pose_cont=no
    ! continue from the same autorefine3D_nano project; with the polish the pair pose error is no
    ! worse than without, within the same 0.5 deg
    ! N29 (pose_cont refactoring, Phase 8; ruling R5): refine3D_nano refine=cont continues from the
    ! autorefine3D_nano project; its poses against the trajectory truth (frame- and
    ! hand-independent pair metric, 20 000 pairs) are no worse than the polar result's, within
    ! 0.5 deg (declared before the first run)
    integer,          parameter           :: NCONT_ITERS = 2, NPAIRS = 20000, PAIR_SEED = 20261001
    real,             parameter           :: MAX_PAIR_LOSS = 0.5
    integer,          parameter           :: NREPROJS = 1000, MASKDIAM = 40, NREFINE_ITERS = 5, NTHR = 8
    integer,          parameter           :: TEST_SEED = 20260923
    integer,          parameter           :: NFRAMES_PER_GROUP = 10
    integer                               :: chdir_status
    real,             parameter           :: TRAJECTORY_SNR    = 0.2
    real,             parameter           :: MIN_VOL_CORR      = 0.90
    real,             parameter           :: MAX_FSC0143       = 5.0
    real,             parameter           :: DOCK_HP           = 100.0
    real,             parameter           :: DOCK_LP           = 5.0
    real                               :: volume_corr, volume_fsc0143
    real                               :: dock_corr_direct, dock_corr_mirrored, dock_corr_selected
    logical                            :: volume_ok
    write(logfhandle,'(a)') '>>> TEST_SINGLE_WORKFLOW:'
    if( .not. cline%defined('suite') ) THROW_HARD('The suite keyword is required; use suite=fcc or suite=wurtzite')
    suite_name = cline%get_carg('suite')
    suite_name = lowercase(suite_name%to_char())
    if( suite_name == 'list' )then
        write(logfhandle,'(a)') 'Available suites for single_workflow:'
        write(logfhandle,'(a)') '  fcc'
        write(logfhandle,'(a)') '  wurtzite'
        return
    endif
    select case(suite_name%to_char())
        case('fcc')
            call cline%set('element', 'Pt')
        case('wurtzite')
            call cline%set('element', 'CdSeW')
        case default
            THROW_HARD('no sub-suite '//suite_name%to_char()//' in single_workflow; use suite=list')
    end select
    call set_fixed_seed(TEST_SEED, propagate=.true.)
    write(logfhandle,'(a,i0)') '>>> Deterministic workflow seed: ', TEST_SEED
    if( .not. cline%defined('smpd') ) call cline%set('smpd', 0.358)
    call cline%set('nthr', NTHR)
    projname = 'test_single_workflow_'//suite_name%to_char()
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
    call validate_reconstructed_volume(simulated_vol, final_volume, params%smpd, BOX, 0.001, real(MOLDIAM), &
        &DOCK_HP, DOCK_LP, MIN_VOL_CORR, MAX_FSC0143, volume_corr, volume_fsc0143, &
        &dock_corr_direct, dock_corr_mirrored, dock_corr_selected, volume_ok, corr_lp=DOCK_LP)
    call return_to_project_dir
    if( .not. volume_ok ) THROW_HARD('TEST_SINGLE_WORKFLOW FAILED: final-volume validation failed')
    ! N29: refine3D_nano refine=cont on the autorefine3D_nano project (projfile, the last stage's
    ! copy, which autorefine3D_nano updated in place)
    call truth_oris%new(NREPROJS, is_ptcl=.false.)
    call truth_oris%read(filepath(filepath(project_dir, REPROJECTION_DIR), TRAJECTORY_ORITAB), [1,NREPROJS])
    cont_projfile = projfile
    call run_proj%read(cont_projfile)
    npairs_ptcls = min(NREPROJS, run_proj%os_ptcl3D%get_noris())
    if( npairs_ptcls < 2 ) THROW_HARD('TEST_SINGLE_WORKFLOW FAILED: the autorefine3D_nano project holds no particles')
    allocate(pinds(npairs_ptcls))
    pinds = [(i, i=1,npairs_ptcls)]
    pair_polar = pair_pose_error(run_proj%os_ptcl3D, truth_oris, pinds, pinds, NPAIRS, PAIR_SEED, frac5)
    call enter_workflow_stage(CONT_DIR, cont_projfile)
    call cline_cont%set('prg',      'refine3D_nano')
    call cline_cont%set('mkdir',    'no')
    call cline_cont%set('projfile', cont_projfile%to_char())
    call cline_cont%set('vol1',     final_volume%to_char())
    call cline_cont%set('refine',   'cont')
    call cline_cont%set('smpd',     params%smpd)
    call cline_cont%set('pgrp',     'c1')
    call cline_cont%set('lp',       1.5)
    call cline_cont%set('mskdiam',  MASKDIAM)
    call cline_cont%set('maxits',   NCONT_ITERS)
    call cline_cont%set('nthr',     params%nthr)
    call xref3Dnano%execute(cline_cont)
    call run_proj%read(cont_projfile)
    pair_cont = pair_pose_error(run_proj%os_ptcl3D, truth_oris, pinds, pinds, NPAIRS, PAIR_SEED, frac5)
    call return_to_project_dir
    write(logfhandle,'(a,f8.3,a,f8.3,a)') 'single_workflow pair pose error against the truth: autorefine3D_nano ', &
        &pair_polar, ' deg, then refine3D_nano refine=cont ', pair_cont, ' deg'
    if( pair_cont > pair_polar + MAX_PAIR_LOSS ) &
        &THROW_HARD('TEST_SINGLE_WORKFLOW FAILED: refine=cont worsened the poses of autorefine3D_nano')
    ! N29, Phase 9: refine3D_nano with and without the polish from the autorefine3D_nano project
    aref_projfile = projfile
    call run_nano_refine(NOPOLISH_DIR, 'no',  pair_nopolish)
    call run_nano_refine(POLISH_DIR,   'yes', pair_polish)
    write(logfhandle,'(a,f8.3,a,f8.3,a)') 'single_workflow pair pose error against the truth: refine3D_nano pose_cont=no ', &
        &pair_nopolish, ' deg, pose_cont=yes ', pair_polish, ' deg'
    if( pair_polish > pair_nopolish + MAX_PAIR_LOSS ) &
        &THROW_HARD('TEST_SINGLE_WORKFLOW FAILED: the polish worsened the poses of refine3D_nano')
    call run_proj%kill
    call truth_oris%kill
    write(logfhandle,'(a,a,a,f7.4,a,f7.4,a,f7.4,a,f7.4,a,f7.2,a)') &
        &'PASS: single_workflow ', suite_name%to_char(), ' docking correlation direct=', dock_corr_direct, &
        &', mirrored=', dock_corr_mirrored, &
        &', selected=', dock_corr_selected, ', soft-masked correlation to 5 A=', volume_corr, &
        &', FSC=0.143 at ', volume_fsc0143, ' A'
    call simple_end('**** SIMPLE_TEST_SINGLE_WORKFLOW NORMAL STOP ****')

contains
    !> NCONT_ITERS refine3D_nano iterations (refine=neigh, objfun=cc) from the autorefine3D_nano
    !! project with pose_cont, in their own stage directory; the pair pose error of the result
    subroutine run_nano_refine( dir, pose_cont, pair_err )
        character(len=*), intent(in)  :: dir, pose_cont
        real,             intent(out) :: pair_err
        type(cmdline) :: cline_nano
        nano_projfile = aref_projfile
        call enter_workflow_stage(dir, nano_projfile)
        call cline_nano%set('prg',       'refine3D_nano')
        call cline_nano%set('mkdir',     'no')
        call cline_nano%set('projfile',  nano_projfile%to_char())
        call cline_nano%set('vol1',      final_volume%to_char())
        call cline_nano%set('pose_cont', pose_cont)
        call cline_nano%set('smpd',      params%smpd)
        call cline_nano%set('pgrp',      'c1')
        call cline_nano%set('lp',        1.5)
        call cline_nano%set('mskdiam',   MASKDIAM)
        call cline_nano%set('maxits',    NCONT_ITERS)
        call cline_nano%set('nthr',      params%nthr)
        call xref3Dnano%execute(cline_nano)
        call cline_nano%kill
        call run_proj%read(nano_projfile)
        pair_err = pair_pose_error(run_proj%os_ptcl3D, truth_oris, pinds, pinds, NPAIRS, PAIR_SEED, frac5)
        call return_to_project_dir
    end subroutine run_nano_refine

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
