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
    real,             parameter :: INNER_FRAC    = 0.30 ! binary crystals: interior atoms, by distance from the centre
    real,             parameter :: MAX_A_ERR     = 0.02 ! binary crystals: fitted lattice constant against the table
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
    if( ndetected > 0 ) call check_binary_lattice
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

    ! binary crystals (rocksalt, zincblende, wurtzite): the interior atoms have the coordination of the first shell
    ! (6, 4, 4) and the fitted lattice constant is the table's within MAX_A_ERR
    subroutine check_binary_lattice
        use simple_defs_atoms,         only: get_lattice_params
        use simple_nanoparticle_utils, only: fit_lattice, binary_lattice
        use simple_srch_sort_loc,      only: hpsort
        character(len=10)     :: crystal_system
        character(len=5)      :: el_ucase
        character(len=4096)   :: header
        real,    allocatable  :: vals(:), cendist(:), cn(:), centers(:,:), sorted(:)
        real    :: a_tab(3), a_fit(3), rinner
        integer :: funit, ios, i, icn, icen, ncol, cn_expected, ninner, nwrong
        el_ucase = upperCase(trim(adjustl(params%element)))
        call get_lattice_params(el_ucase, crystal_system, a_tab)
        if( .not. binary_lattice(crystal_system) ) return
        cn_expected = 4
        if( trim(crystal_system) == 'rocksalt' ) cn_expected = 6
        ! fitted lattice constant of the detected atoms
        allocate(centers(3,ndetected))
        do i = 1,ndetected
            centers(:,i) = detected_atoms%get_coord(i)
        enddo
        call fit_lattice(el_ucase, centers, a_fit)
        write(logfhandle,'(a,a,a,3f9.4,a,3f9.4)') '    ', trim(crystal_system), ' lattice fitted ', a_fit, ', table ', a_tab
        if( abs(a_fit(1) - a_tab(1)) > MAX_A_ERR * a_tab(1) ) call fail('fitted lattice constant is more than 2% from the table')
        ! interior coordination from atoms_stats.csv, columns found by name
        if( .not. file_exists(ATOM_STATS) ) return
        call fopen(funit, file=string(ATOM_STATS), status='old', action='read', iostat=ios)
        if( ios /= 0 )then
            call fail('could not open '//ATOM_STATS)
            return
        endif
        read(funit,'(a)',iostat=ios) header
        call column_of(header, 'CN_STD',  icn,  ncol)
        call column_of(header, 'CENDIST', icen, ncol)
        if( icn == 0 .or. icen == 0 )then
            call fclose(funit)
            call fail(ATOM_STATS//' has no CN_STD or CENDIST column')
            return
        endif
        allocate(vals(ncol), cn(ndetected), cendist(ndetected))
        do i = 1,ndetected
            read(funit,*,iostat=ios) vals
            if( ios /= 0 ) exit
            cn(i)      = vals(icn)
            cendist(i) = vals(icen)
        enddo
        call fclose(funit)
        if( ios /= 0 )then
            call fail(ATOM_STATS//' ended before all detected atoms were read')
            return
        endif
        sorted = cendist
        call hpsort(sorted)
        ninner = max(1, nint(INNER_FRAC * real(ndetected)))
        rinner = sorted(ninner)
        nwrong = count(cendist <= rinner .and. nint(cn) /= cn_expected)
        write(logfhandle,'(a,i0,a,i0,a,i0,a,f6.2,a)') '    interior atoms: ', count(cendist <= rinner), ', with CN /= ', cn_expected,&
            &': ', nwrong, ' (CN range ', minval(cn, mask=cendist <= rinner), ')'
        if( nwrong > 0 ) call fail('an interior atom of the binary crystal has the wrong coordination number')
    end subroutine check_binary_lattice

    ! 1-based column of name in a comma-separated header, 0 if absent, and the number of columns
    subroutine column_of( header, name, icol, ncol )
        character(len=*), intent(in)  :: header, name
        integer,          intent(out) :: icol, ncol
        integer :: i, first
        icol  = 0
        ncol  = 1
        first = 1
        do i = 1,len_trim(header) + 1
            if( i > len_trim(header) .or. header(min(i,len(header)):min(i,len(header))) == ',' )then
                if( trim(adjustl(header(first:i-1))) == name ) icol = ncol
                if( i <= len_trim(header) ) ncol = ncol + 1
                first = i + 1
            endif
        enddo
    end subroutine column_of

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
    use simple_test_gate,        only: test_gate
    use simple_rnd,              only: ran3
    use simple_nano_species,     only: gauss_filter3D, fit_gauss_width, enclosed_fraction
    use simple_nanoparticle_utils, only: find_rMax
    class(commander_test_species_discovery), intent(inout) :: self
    class(cmdline),                          intent(inout) :: cline
    ! floors of section 10.6 of doc/implementation_notes/planned/species_discovery.md, set before the first run
    real,    parameter :: MIN_PT_RECALL    = 0.98   ! Pt over its eligible atoms: runs with discovery and the pure case
    real,    parameter :: MIN_ALL_RECALL   = 0.95   ! all eligible atoms of all species, runs with discovery
    integer, parameter :: MAX_FALSE        = 1      ! false atoms per run
    real,    parameter :: MIN_LIGHT_RECALL = 0.90   ! the second element over its eligible atoms, runs with the list
    real,    parameter :: MAX_PLAIN_LIGHT  = 0.90   ! light case: element=Pt finds fewer than this of the eligible Al
    real,    parameter :: MIN_HALF_AGREE   = 0.95   ! label agreement between the half maps
    real,    parameter :: MIN_ELIG_SNR     = 6.5    ! eligibility: filtered peak of the atom's own kernel over filtered noise
    real,    parameter :: MATCH_NN_FRAC    = 0.3    ! found atom matches a generating atom within this many d_NN
    real,    parameter :: SIGMA_REF        = 0.4196 ! A, width of the detection template (B 13.9 A**2)
    ! fixtures: the Pt lattice of single_atoms_stats with the element kernels, made to look like the real Pt map
    ! of the maintainer's ruling of 2026-10-09 (section 10.6): per-atom B, the measured signal transfer, noise
    ! shaped by the measured background spectrum and scaled to the measured core-peak signal-to-noise
    real,    parameter :: SMPD_FIX    = 0.358
    real,    parameter :: B_CORE      = 12.8    ! A**2, per-atom B = B_CORE + B_RISE (r / r_max)**2 on top of the blur
    real,    parameter :: B_RISE      = 6.
    real,    parameter :: PEAK_SD     = 17.     ! mean core-Pt peak of the filtered clean render over the noise sdev
    real,    parameter :: R_NOISE     = 18.     ! A, the noise sdev is measured outside this radius
    real,    parameter :: INNER_FRAC  = 0.3     ! core atoms: the inner 30% of the Pt atoms by radius
    real,    parameter :: SECOND_FRAC = 0.25    ! fraction of the atoms given the second element
    integer, parameter :: BOX_ATOM    = 64      ! ground truth: one atom of each element
    integer, parameter :: BOX_ELIG    = 64      ! eligibility: each generating atom alone
    ! signal transfer relative to 2-2.5 A at shell centres (1/A): 0.37 beyond 20 A, 0.50 at 10-20 A, 0.86 at 5-10 A,
    ! 1 from 5 A to 1.6 A, then a cosine fall to 0 at 1.1 A
    real,    parameter :: TF_S(4)     = [0.025, 0.075, 0.15, 0.2]
    real,    parameter :: TF_V(4)     = [0.37, 0.50, 0.86, 1.0]
    real,    parameter :: TF_FLAT     = 1. / 1.6
    real,    parameter :: TF_END      = 1. / 1.1
    ! background amplitude per resolution shell (bounds in A), relative; constant beyond the outer shell centres
    real,    parameter :: NZ_RES(12)  = [20., 10., 5., 3.3, 2.5, 2., 1.6, 1.3, 1.1, 0.95, 0.8, 0.72]
    real,    parameter :: NZ_V(11)    = [0.364, 0.339, 0.415, 0.322, 0.249, 0.172, 0.136, 0.119, 0.105, 0.096, 0.088]
    real,    parameter :: CUTOFF_ELIG = 12. * SMPD_FIX ! the kernel cutoff of simulate_nanoparticle pdbfile=
    integer, parameter :: NELEM = 3, NCASES = 4, NRUNS = 16
    character(len=2),  parameter :: ELEMS(NELEM)       = ['PT', 'NI', 'AL']
    character(len=5),  parameter :: CASES(NCASES)      = ['pure ', 'alloy', 'core ', 'light']
    character(len=2),  parameter :: SECOND(NCASES)     = ['  ', 'NI', 'NI', 'AL']
    integer,           parameter :: RUN_CASE(NRUNS)    = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4]
    character(len=11), parameter :: RUN_NAME(NRUNS)    = [character(len=11) :: 'pt', 'pt_disc', 'ptni', 'noel',&
        &'noel_disc', 'pt', 'ptni', 'ptni_halves', 'noel', 'noel_disc', 'pt', 'ptni', 'pt', 'ptal', 'noel', 'noel_disc']
    character(len=5),  parameter :: RUN_EL(NRUNS)      = [character(len=5) :: 'Pt', 'Pt', 'Pt,Ni', '', '', 'Pt',&
        &'Pt,Ni', 'Pt,Ni', '', '', 'Pt', 'Pt,Ni', 'Pt', 'Pt,Al', '', '']
    logical,           parameter :: RUN_KEY(NRUNS)     = [.false., .true., .false., .false., .true., .false., .false.,&
        &.false., .false., .true., .false., .false., .false., .false., .false., .true.]
    logical,           parameter :: RUN_HALVES(NRUNS)  = [.false., .false., .false., .false., .false., .false., .false.,&
        &.true., .false., .false., .false., .false., .false., .false., .false., .false.]
    character(len=*),  parameter :: PRODUCTS(5) = [character(len=12) :: 'map_ATMS.pdb', 'map_BIN.mrc', 'map_CC.mrc',&
        &'map_MSK.mrc', 'map_SIM.mrc']
    type(commander_simulate_nanoparticle) :: xsim
    type(commander_detect_atoms)          :: xdet
    type(cmdline)                         :: cline_sim, cline_det
    type(test_gate)                       :: gate
    type(atoms)                           :: lattice
    type(string)                          :: cwd_saved, fixture_root, case_dir
    real,             allocatable :: lat(:,:), gxyz(:,:), gb(:)
    character(len=2), allocatable :: gel(:)
    logical,          allocatable :: elig(:)
    real    :: dnn_gen, cen(3), rmax, int_el(NELEM), peak_el(NELEM), nz_s(11)
    integer :: status, nthr, nlat, icase, ir, i
    logical :: l_passed
    write(logfhandle,'(a)') '>>> TEST_SPECIES_DISCOVERY:'
    ! every commander run in process reseeds in parameters%new, from SIMPLE_SEED when it is set (CTest sets
    ! SIMPLE_SEED=20260923, the setting the floors were validated with) and from /dev/urandom otherwise: run by hand,
    ! export SIMPLE_SEED=20260923 first, or the fixtures differ from run to run
    call set_fixed_seed(20261008)
    nthr = 8
    if( cline%defined('nthr') ) nthr = cline%get_iarg('nthr')
    call simple_getcwd(cwd_saved)
    call gate%new(filepath(cwd_saved, 'metrics.tsv'))
    fixture_root = filepath(cwd_saved, 'test_species_discovery_'//int2str(get_process_id()))
    if( dir_exists(fixture_root) ) call simple_rmdir(fixture_root)
    call simple_mkdir(fixture_root)
    call simple_chdir(fixture_root, status)
    if( status /= 0 ) THROW_HARD('TEST_SPECIES_DISCOVERY FAILED: could not enter fixture directory')
    ! the Pt lattice of single_atoms_stats: positions only
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
    do i = 1,11
        nz_s(i) = 0.5 * (1. / NZ_RES(i) + 1. / NZ_RES(i+1))
    enddo
    call ground_truth
    do icase = 1,NCASES
        case_dir = filepath(fixture_root, trim(CASES(icase)))
        call simple_mkdir(case_dir)
        call simple_chdir(case_dir, status)
        if( status /= 0 ) THROW_HARD('TEST_SPECIES_DISCOVERY FAILED: could not enter a case directory')
        write(logfhandle,'(a)') '>>> CASE '//trim(CASES(icase))
        call make_fixture(icase)
        do ir = 1,NRUNS
            if( RUN_CASE(ir) == icase ) call run_detect(ir)
        enddo
        call compare_products(icase)
        do ir = 1,NRUNS
            if( RUN_CASE(ir) == icase ) call evaluate(ir)
        enddo
        call light_plain_recall(icase)
        deallocate(gxyz, gb, gel, elig)
        call simple_chdir(fixture_root, status)
    enddo
    call simple_chdir(cwd_saved, status)
    if( status /= 0 ) THROW_HARD('TEST_SPECIES_DISCOVERY FAILED: could not restore original directory')
    ! the fixture tree goes, pass or fail; metrics.tsv stays in the test's directory
    call simple_rmdir(fixture_root)
    l_passed = gate%passed()
    call gate%kill
    if( l_passed )then
        write(logfhandle,'(a)') 'PASS: species_discovery'
        call simple_end('**** SIMPLE_TEST_SPECIES_DISCOVERY NORMAL STOP ****')
    else
        THROW_HARD('TEST_SPECIES_DISCOVERY FAILED')
    endif

  contains

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

    integer function elem_index( el )
        character(len=2), intent(in) :: el
        elem_index = findloc(ELEMS, el, dim=1)
    end function elem_index

    ! one isolated atom of each element at the core B, rendered and filtered as the particles are; its aperture
    ! intensity, measured as the fitted atoms' are, is the generating intensity
    subroutine ground_truth
        type(atoms) :: one
        type(image) :: vol
        real, allocatable :: r(:,:,:)
        integer :: ie
        do ie = 1,NELEM
            call one%new(1, dummy=.true.)
            call one%set_name(1, ELEMS(ie)//'  ')
            call one%set_element(1, ELEMS(ie))
            call one%set_coord(1, real(BOX_ATOM / 2) * [SMPD_FIX, SMPD_FIX, SMPD_FIX])
            call one%set_occupancy(1, 1.)
            call one%set_beta(1, B_CORE)
            call one%set_num(1, 1)
            call one%set_resnum(1, 1)
            call one%writepdb(string('atom.pdb'))
            call one%kill
            call cline_sim%set('prg',      'simulate_nanoparticle')
            call cline_sim%set('pdbfile',  'atom.pdb')
            call cline_sim%set('pdb_bfac', 'yes')
            call cline_sim%set('box',      BOX_ATOM)
            call cline_sim%set('smpd',     SMPD_FIX)
            call cline_sim%set('nthr',     nthr)
            call cline_sim%set('outvol',   'atom.mrc')
            call cline_sim%set('pdbout',   'atom_out.pdb')
            call xsim%execute(cline_sim)
            call cline_sim%kill
            call vol%new([BOX_ATOM,BOX_ATOM,BOX_ATOM], SMPD_FIX)
            call vol%read(string('atom.mrc'))
            call radial_filter(vol, .true.)
            r = vol%get_rmat()
            call vol%kill
            int_el(ie)  = aperture_intensity(r, real(BOX_ATOM / 2 + 1))
            peak_el(ie) = maxval(r)
            call gate%report('single atom '//ELEMS(ie)//': aperture intensity', int_el(ie))
            call gate%report('single atom '//ELEMS(ie)//': peak', peak_el(ie))
            call del_file('atom.mrc')
        enddo
        do ie = 2,NELEM
            call gate%report('single atom ratio '//ELEMS(ie)//'/PT: aperture intensity', int_el(ie) / int_el(1))
            call gate%report('single atom ratio '//ELEMS(ie)//'/PT: peak', peak_el(ie) / peak_el(1))
        enddo
    end subroutine ground_truth

    ! the aperture intensity of discover_species for an atom alone at voxel c: the voxel sum within d_NN / 2 times
    ! smpd**3 over the enclosed fraction of the Gaussian fitted in the same sphere
    real function aperture_intensity( r, c )
        real, intent(in) :: r(:,:,:), c
        real, allocatable :: y(:), r2(:)
        real    :: rad, dd, amp, bfac
        integer :: i1, i2, i3
        rad = 0.5 * dnn_gen
        allocate(y(0), r2(0))
        do i3 = 1,size(r,3)
            do i2 = 1,size(r,2)
                do i1 = 1,size(r,1)
                    dd = sum((real([i1,i2,i3]) - c)**2) * SMPD_FIX**2
                    if( dd > rad**2 ) cycle
                    y  = [y, r(i1,i2,i3)]
                    r2 = [r2, dd]
                enddo
            enddo
        enddo
        ! a flat width prior: the isolated atom's own density determines its width
        call fit_gauss_width(y, r2, 1., log(13.9), 100., amp, bfac)
        aperture_intensity = SMPD_FIX**3 * sum(y) / enclosed_fraction(rad / sqrt(bfac / (8. * PI**2)))
    end function aperture_intensity

    ! the measured signal transfer (signal) or background amplitude (noise) as a function of spatial frequency s, 1/A
    real function profile( s, l_signal )
        real,    intent(in) :: s
        logical, intent(in) :: l_signal
        if( l_signal )then
            profile = interp(s, TF_S, TF_V)
            if( s >= TF_S(size(TF_S)) ) profile = 1.
            if( s > TF_FLAT ) profile = 0.5 * (1. + cos(PI * (s - TF_FLAT) / (TF_END - TF_FLAT)))
            if( s >= TF_END ) profile = 0.
        else
            profile = interp(s, nz_s, NZ_V)
        endif
    end function profile

    ! linear between the points, constant beyond the ends
    real function interp( s, xs, vs )
        real, intent(in) :: s, xs(:), vs(:)
        integer :: m
        interp = vs(size(vs))
        if( s <= xs(1) )then
            interp = vs(1)
            return
        endif
        do m = 2,size(xs)
            if( s <= xs(m) )then
                interp = vs(m-1) + (vs(m) - vs(m-1)) * (s - xs(m-1)) / (xs(m) - xs(m-1))
                return
            endif
        enddo
    end function interp

    ! multiply the Fourier transform of a cubic real-space volume by the signal or noise profile
    subroutine radial_filter( vol, l_signal )
        class(image), intent(inout) :: vol
        logical,      intent(in)    :: l_signal
        integer :: lims(3,2), h, k, l, phys(3), ld(3)
        real    :: s
        ld = vol%get_ldim()
        call vol%fft()
        lims = vol%loop_lims(2)
        do h = lims(1,1),lims(1,2)
            do k = lims(2,1),lims(2,2)
                do l = lims(3,1),lims(3,2)
                    s    = sqrt(real(h*h + k*k + l*l)) / (real(ld(1)) * SMPD_FIX)
                    phys = vol%comp_addr_phys(h,k,l)
                    call vol%set_fcomp([h,k,l], phys, vol%get_fcomp([h,k,l], phys) * profile(s, l_signal))
                enddo
            enddo
        enddo
        call vol%ifft()
    end subroutine radial_filter

    ! sdev of a volume outside R_NOISE from the lattice centre
    real function outside_sdev( a )
        real, intent(in) :: a(:,:,:)
        real(dp) :: s1, s2
        integer  :: i1, i2, i3, m
        s1 = 0._dp
        s2 = 0._dp
        m  = 0
        do i3 = 1,size(a,3)
            do i2 = 1,size(a,2)
                do i1 = 1,size(a,1)
                    if( norm2(real([i1,i2,i3] - 1) * SMPD_FIX - cen) <= R_NOISE ) cycle
                    s1 = s1 + a(i1,i2,i3)
                    s2 = s2 + a(i1,i2,i3)**2
                    m  = m + 1
                enddo
            enddo
        enddo
        outside_sdev = real(sqrt(s2 / m - (s1 / m)**2))
    end function outside_sdev

    ! the case's elements and B factors on the lattice, rendered with the element kernels and filtered by the signal
    ! transfer, with shaped noise; the eligible atoms from each atom's own kernel and the added noise, both filtered
    ! with the detection template
    subroutine make_fixture( ic )
        integer, intent(in) :: ic
        type(atoms)       :: model
        type(image)       :: clean, nimg, nimg2
        real, allocatable :: keys(:), rad2(:), cmat(:,:,:), nmat(:,:,:), fmat(:,:,:), rpt(:)
        integer, allocatable :: order(:)
        real    :: sd_filt, snr, rmin_snr, peak, rin, scale
        integer :: j, nsecond, ncore, ipt(3)
        allocate(gel(nlat))
        gel = ELEMS(1)
        nsecond = nint(SECOND_FRAC * real(nlat))
        select case(trim(CASES(ic)))
            case('alloy', 'light')
                ! a random quarter
                allocate(keys(nlat))
                do j = 1,nlat
                    keys(j) = ran3()
                enddo
                order = [(j, j=1,nlat)]
                call hpsort(keys, order)
                gel(order(:nsecond)) = SECOND(ic)
            case('core')
                ! the innermost quarter under a Pt skin; a light species in the surface layer is not a target
                rad2  = sum((lat - spread(cen, 2, nlat))**2, dim=1)
                order = [(j, j=1,nlat)]
                call hpsort(rad2, order)
                gel(order(:nsecond)) = SECOND(ic)
        end select
        call model%new(nlat, dummy=.true.)
        do j = 1,nlat
            call model%set_name(j, gel(j)//'  ')
            call model%set_element(j, gel(j))
            call model%set_coord(j, lat(:,j))
            call model%set_occupancy(j, 1.)
            call model%set_beta(j, B_CORE + B_RISE * sum((lat(:,j) - cen)**2) / rmax**2)
            call model%set_num(j, j)
            call model%set_resnum(j, j)
        enddo
        call model%writepdb(string('model.pdb'))
        call model%kill
        ! the generating model as written: two decimals of B
        call model%new(string('model.pdb'))
        allocate(gxyz(3,nlat), gb(nlat))
        do j = 1,nlat
            gxyz(:,j) = model%get_coord(j)
            gb(j)     = model%get_beta(j)
            gel(j)    = model%get_element(j)
        enddo
        call model%kill
        call cline_sim%set('prg',      'simulate_nanoparticle')
        call cline_sim%set('pdbfile',  'model.pdb')
        call cline_sim%set('pdb_bfac', 'yes')
        call cline_sim%set('box',      BOX)
        call cline_sim%set('smpd',     SMPD_FIX)
        call cline_sim%set('nthr',     nthr)
        call cline_sim%set('outvol',   'clean.mrc')
        call cline_sim%set('pdbout',   'model_out.pdb')
        call xsim%execute(cline_sim)
        call cline_sim%kill
        call clean%new([BOX,BOX,BOX], SMPD_FIX)
        call clean%read(string('clean.mrc'))
        call del_file('clean.mrc')
        call radial_filter(clean, .true.)
        cmat = clean%get_rmat()
        call clean%kill
        ! mean peak of the core Pt atoms (the inner INNER_FRAC of them by radius) on the filtered clean render
        rpt = pack(sqrt(sum((gxyz - spread(cen, 2, nlat))**2, dim=1)), gel == ELEMS(1))
        call hpsort(rpt)
        rin   = rpt(max(1, nint(INNER_FRAC * real(size(rpt)))))
        peak  = 0.
        ncore = 0
        do j = 1,nlat
            if( gel(j) /= ELEMS(1) .or. norm2(gxyz(:,j) - cen) > rin ) cycle
            ipt   = nint(gxyz(:,j) / SMPD_FIX) + 1
            peak  = peak + cmat(ipt(1),ipt(2),ipt(3))
            ncore = ncore + 1
        enddo
        peak = peak / real(ncore)
        ! shaped noise, scaled so that the mean core peak is PEAK_SD noise sdevs outside the particle
        call nimg%new([BOX,BOX,BOX], SMPD_FIX)
        call nimg%gauran(0., 1.)
        call radial_filter(nimg, .false.)
        nmat  = nimg%get_rmat()
        scale = peak / (PEAK_SD * outside_sdev(nmat))
        if( any(RUN_HALVES .and. RUN_CASE == ic) )then
            ! half maps with independent noise of sqrt(2) times the sdev, whose average is the map
            call nimg2%new([BOX,BOX,BOX], SMPD_FIX)
            call nimg2%gauran(0., 1.)
            call radial_filter(nimg2, .false.)
            call write_vol(cmat + sqrt(2.) * scale * nmat, 'even.mrc')
            call write_vol(cmat + sqrt(2.) * scale * nimg2%get_rmat(), 'odd.mrc')
            nmat = scale * (nmat + nimg2%get_rmat()) / sqrt(2.)
            call nimg2%kill
        else
            nmat = scale * nmat
        endif
        call nimg%kill
        call write_vol(cmat + nmat, 'map.mrc')
        write(logfhandle,'(a,i0,a,a,a,i0,a,f9.3,a,es12.4)') '    model: ', count(gel == ELEMS(1)), ' PT, ', SECOND(ic), ' ',&
            &count(gel /= ELEMS(1)), '; mean core peak ', peak, ', noise sdev outside the particle ', outside_sdev(nmat)
        call gate%report(trim(CASES(ic))//': mean core-PT peak of the filtered clean render', peak)
        call gate%report(trim(CASES(ic))//': noise sdev outside the particle', outside_sdev(nmat))
        ! the added noise filtered with the detection template
        allocate(fmat(BOX,BOX,BOX))
        call gauss_filter3D(nmat, SIGMA_REF / SMPD_FIX, fmat)
        sd_filt = sqrt(sum((fmat - sum(fmat) / real(size(fmat)))**2) / real(size(fmat)))
        deallocate(cmat, nmat, fmat)
        ! eligibility from each atom's own kernel, alone, at its sub-voxel offset
        allocate(elig(nlat))
        rmin_snr = huge(1.)
        do j = 1,nlat
            snr      = own_peak(gel(j), gb(j), gxyz(:,j)) / sd_filt
            elig(j)  = snr >= MIN_ELIG_SNR
            rmin_snr = min(rmin_snr, snr)
        enddo
        call prune_eligibility(ic)
        call gate%report(trim(CASES(ic))//': filtered noise sdev', sd_filt)
        call gate%report(trim(CASES(ic))//': lowest eligibility signal-to-noise', rmin_snr)
        call gate%report(trim(CASES(ic))//': eligible PT atoms', real(count(elig .and. gel == ELEMS(1))))
        if( trim(SECOND(ic)) /= '' )then
            call gate%report(trim(CASES(ic))//': eligible '//SECOND(ic)//' atoms', real(count(elig .and. gel == SECOND(ic))))
        endif
    end subroutine make_fixture

    ! the pruning of recovered atoms (section 3.9) on the generating model: in the zone beyond the 85% radius quantile
    ! of the Pt atoms, an atom with fewer Pt atoms within the contact cutoff than the threshold discard_atoms would
    ! derive for the Pt atoms is not eligible (ruling of 2026-10-09, second)
    subroutine prune_eligibility( ic )
        integer, intent(in) :: ic
        integer, parameter :: CS_CEIL = 12        ! fcc
        real,    allocatable :: rpt(:)
        integer, allocatable :: cs(:)
        real    :: cpt(3), rzone, rmax_c
        integer :: j, l, npt, cn, cthres, nexcl_pt, nexcl_2
        logical :: l_pt
        npt    = count(gel == ELEMS(1))
        cpt    = sum(gxyz, dim=2, mask=spread(gel == ELEMS(1), 1, 3)) / real(npt)
        rmax_c = find_rMax('Pt')
        allocate(cs(nlat), source=0)
        do j = 1,nlat
            do l = 1,nlat
                if( l == j .or. gel(l) /= ELEMS(1) ) cycle
                if( norm2(gxyz(:,j) - gxyz(:,l)) < rmax_c ) cs(j) = cs(j) + 1
            enddo
        enddo
        rpt = pack(sqrt(sum((gxyz - spread(cpt, 2, nlat))**2, dim=1)), gel == ELEMS(1))
        call hpsort(rpt)
        rzone  = rpt(nint(0.85 * real(npt)))
        cthres = CS_CEIL
        do cn = 1,CS_CEIL
            if( real(count(cs >= cn .and. gel == ELEMS(1))) / real(npt) * 100. <= 95. )then
                cthres = cn
                exit
            endif
        enddo
        if( cthres > CS_CEIL / 2 ) cthres = CS_CEIL / 2
        nexcl_pt = 0
        nexcl_2  = 0
        do j = 1,nlat
            if( norm2(gxyz(:,j) - cpt) <= rzone .or. cs(j) >= cthres ) cycle
            if( .not. elig(j) ) cycle
            elig(j) = .false.
            l_pt    = gel(j) == ELEMS(1)
            if( l_pt )then
                nexcl_pt = nexcl_pt + 1
            else
                nexcl_2 = nexcl_2 + 1
            endif
        enddo
        call gate%report(trim(CASES(ic))//': pruning-policy zone radius', rzone)
        call gate%report(trim(CASES(ic))//': pruning-policy contact threshold', real(cthres))
        call gate%report(trim(CASES(ic))//': PT atoms excluded by the pruning policy', real(nexcl_pt))
        if( trim(SECOND(ic)) /= '' )then
            call gate%report(trim(CASES(ic))//': '//SECOND(ic)//' atoms excluded by the pruning policy', real(nexcl_2))
        endif
    end subroutine prune_eligibility

    subroutine write_vol( arr, fname )
        real,             intent(in) :: arr(:,:,:)
        character(len=*), intent(in) :: fname
        type(image) :: vol
        call vol%new([BOX,BOX,BOX], SMPD_FIX)
        call vol%set_rmat(arr, .false.)
        call vol%write(string(fname))
        call vol%kill
    end subroutine write_vol

    ! peak of one atom's kernel (element, B) rendered alone, filtered by the signal transfer and the detection template
    real function own_peak( el, bfac, xyz )
        character(len=2), intent(in) :: el
        real,             intent(in) :: bfac, xyz(3)
        type(atoms) :: one
        type(image) :: vol
        real, allocatable :: r(:,:,:), f(:,:,:)
        real :: offset(3)
        offset = xyz / SMPD_FIX - real(floor(xyz / SMPD_FIX))
        call one%new(1, dummy=.true.)
        call one%set_element(1, el)
        call one%set_coord(1, (real(BOX_ELIG / 2) + offset) * SMPD_FIX)
        call one%set_occupancy(1, 1.)
        call one%set_beta(1, bfac)
        call vol%new([BOX_ELIG,BOX_ELIG,BOX_ELIG], SMPD_FIX, wthreads=.false.)
        call one%convolve(vol, CUTOFF_ELIG, bfac_pdb=.true.)
        call radial_filter(vol, .true.)
        r = vol%get_rmat()
        allocate(f(BOX_ELIG,BOX_ELIG,BOX_ELIG))
        call gauss_filter3D(r, SIGMA_REF / SMPD_FIX, f)
        own_peak = maxval(f)
        call vol%kill
        call one%kill
    end function own_peak

    ! detect_atoms on the case's map in a directory of its own, with the run's element value, key and half maps
    subroutine run_detect( ir )
        integer, intent(in) :: ir
        integer :: st
        call simple_mkdir(trim(RUN_NAME(ir)))
        call simple_chdir(trim(RUN_NAME(ir)), st)
        call cline_det%set('prg',  'detect_atoms')
        call cline_det%set('vol1', filepath(case_dir, 'map.mrc'))
        call cline_det%set('smpd', SMPD_FIX)
        call cline_det%set('nthr', nthr)
        if( len_trim(RUN_EL(ir)) > 0 ) call cline_det%set('element', trim(RUN_EL(ir)))
        if( RUN_KEY(ir) ) call cline_det%set('discover_species', 'yes')
        if( RUN_HALVES(ir) )then
            call cline_det%set('vol_even', filepath(case_dir, 'even.mrc'))
            call cline_det%set('vol_odd',  filepath(case_dir, 'odd.mrc'))
        endif
        call xdet%execute(cline_det)
        call cline_det%kill
        call simple_chdir(case_dir, st)
    end subroutine run_detect

    ! a run with discovery: the key, a list, or both
    logical function discovery( ir )
        integer, intent(in) :: ir
        discovery = RUN_KEY(ir) .or. index(RUN_EL(ir), ',') > 0
    end function discovery

    ! the five present products are identical across the runs of a case with the same first element (none, or Pt)
    subroutine compare_products( ic )
        integer, intent(in) :: ic
        integer :: ir, iref, ip
        logical :: same
        do ir = 1,NRUNS
            if( RUN_CASE(ir) /= ic ) cycle
            ! the first run of the case with the same level-1 element
            iref = findloc([(RUN_CASE(i) == ic .and. RUN_EL(i)(1:2) == RUN_EL(ir)(1:2), i=1,NRUNS)], .true., dim=1)
            if( iref == ir ) cycle
            same = .true.
            do ip = 1,size(PRODUCTS)
                if( .not. files_identical(string(trim(RUN_NAME(iref))//'/'//trim(PRODUCTS(ip))),&
                    &string(trim(RUN_NAME(ir))//'/'//trim(PRODUCTS(ip)))) )then
                    same = .false.
                    write(logfhandle,'(a)') '    differs: '//trim(PRODUCTS(ip))
                endif
            enddo
            call gate%check(trim(CASES(ic))//'/'//trim(RUN_NAME(ir))//': present products identical to '//trim(RUN_NAME(iref)), same)
        enddo
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

    ! the numeric rows of a comma-separated table with a header line; the header gives the column count
    subroutine read_table( fname, rows )
        character(len=*),  intent(in)  :: fname
        real, allocatable, intent(out) :: rows(:,:)
        character(len=1024) :: line
        integer :: u, ios, nrow, ncol, j, l
        nrow = nlines(string(fname)) - 1
        open(newunit=u, file=fname, status='old', action='read')
        read(u,'(a)') line
        ncol = count([(line(l:l) == ',', l=1,len_trim(line))]) + 1
        allocate(rows(ncol,max(nrow,0)), source=0.)
        do j = 1,nrow
            read(u,'(a)') line
            do l = 1,len_trim(line)
                if( line(l:l) == ',' ) line(l:l) = ' '
            enddo
            read(line,*,iostat=ios) rows(:,j)
            if( ios /= 0 ) call gate%check(fname//': row '//int2str(j)//' readable', .false.)
        enddo
        close(u)
    end subroutine read_table

    ! recall of a set of generating atoms: the fraction of them with a found atom on their site
    real function recall_of( sel, hits )
        logical, intent(in) :: sel(:)
        integer, intent(in) :: hits(:)
        recall_of = real(count(sel .and. hits > 0)) / real(max(count(sel), 1))
    end function recall_of

    ! a run against the generating model
    subroutine evaluate( ir )
        integer, intent(in) :: ir
        character(len=:), allocatable :: tag, csv, txt, pdb, rad_csv
        character(len=2), allocatable :: syms(:)
        real,    allocatable :: rows(:,:), xyz(:,:), radial(:,:), rsort(:)
        integer, allocatable :: match(:), hits(:), order(:), cls(:), sel(:)
        type(atoms) :: found_atms
        character(len=2) :: el2
        real    :: tol, dpt, ratio, ratio_gen, rise_fit, rise_gen, rin, rout
        integer, parameter :: MAXLEV = 8
        integer :: nrow, j, nfalse, nwrong, k_found, nrec, nsh, nin, nout, kc, j1, j2, ic, nfound(0:2,0:MAXLEV)
        logical :: l_list, l_disc
        ic     = RUN_CASE(ir)
        el2    = SECOND(ic)
        tag    = trim(CASES(ic))//'/'//trim(RUN_NAME(ir))
        l_list = index(RUN_EL(ir), ',') > 0
        l_disc = discovery(ir)
        csv    = trim(RUN_NAME(ir))//'/map_species.csv'
        txt    = trim(RUN_NAME(ir))//'/map_species.txt'
        pdb    = trim(RUN_NAME(ir))//'/map_species.pdb'
        ! found atoms: the species table with discovery, the present atoms otherwise
        if( l_disc )then
            call gate%check(tag//': species files written', file_exists(csv) .and. file_exists(txt) .and. file_exists(pdb))
            if( .not. (file_exists(csv) .and. file_exists(txt) .and. file_exists(pdb)) ) return
            call read_table(csv, rows)
            nrow = size(rows, 2)
            xyz  = rows(3:5,:)
            cls  = nint(rows(16,:))
        else
            call gate%check(tag//': atoms written', file_exists(trim(RUN_NAME(ir))//'/map_ATMS.pdb'))
            if( .not. file_exists(trim(RUN_NAME(ir))//'/map_ATMS.pdb') ) return
            call found_atms%new(string(trim(RUN_NAME(ir))//'/map_ATMS.pdb'))
            nrow = found_atms%get_n()
            allocate(xyz(3,nrow))
            do j = 1,nrow
                xyz(:,j) = found_atms%get_coord(j)
            enddo
            call found_atms%kill
        endif
        ! match every found atom to its nearest generating atom; unmatched atoms and second atoms on a site are false
        tol = MATCH_NN_FRAC * dnn_gen
        allocate(match(nrow), source=0)
        allocate(hits(nlat),  source=0)
        do j = 1,nrow
            match(j) = minloc(sum((gxyz - spread(xyz(:,j), 2, nlat))**2, dim=1), dim=1)
            if( norm2(gxyz(:,match(j)) - xyz(:,j)) > tol )then
                match(j) = 0
            else
                hits(match(j)) = hits(match(j)) + 1
            endif
        enddo
        nfalse = count(match == 0) + sum(max(hits - 1, 0))
        do j = 1,nrow
            if( match(j) /= 0 ) cycle
            ! a false atom: its position, its distance to the nearest Pt atom and its detection z
            dpt = minval(sqrt(sum((gxyz - spread(xyz(:,j), 2, nlat))**2, dim=1)), mask=gel == ELEMS(1))
            if( l_disc )then
                write(logfhandle,'(a,3f9.3,a,f7.3,a,i0,a,f8.3)') '    false atom at ', xyz(:,j), ' A, nearest PT ', dpt,&
                    &' A, stage ', nint(rows(8,j)), ', detection z ', rows(10,j)
            else
                write(logfhandle,'(a,3f9.3,a,f7.3,a)') '    false atom at ', xyz(:,j), ' A, nearest PT ', dpt, ' A'
            endif
        enddo
        call gate%report(tag//': found atoms', real(nrow))
        call gate%metric(tag//': false atoms', real(nfalse), real(MAX_FALSE), nfalse <= MAX_FALSE)
        if( l_disc .or. trim(el2) == '' )then
            call gate%metric(tag//': PT recall over eligible', recall_of(elig .and. gel == ELEMS(1), hits), MIN_PT_RECALL,&
                &recall_of(elig .and. gel == ELEMS(1), hits) >= MIN_PT_RECALL)
        else
            ! the present path on a mixed particle: reported (ruling of 2026-10-09, second)
            call gate%report(tag//': PT recall over eligible', recall_of(elig .and. gel == ELEMS(1), hits))
        endif
        if( trim(el2) /= '' )then
            call gate%report(tag//': '//el2//' recall over eligible', recall_of(elig .and. gel == el2, hits))
        endif
        if( .not. l_disc ) return
        call gate%metric(tag//': recall over all eligible', recall_of(elig, hits), MIN_ALL_RECALL,&
            &recall_of(elig, hits) >= MIN_ALL_RECALL)
        k_found = nint(report_value(txt, 'K'))
        nrec    = count(nint(rows(2,:)) == 1)
        call gate%report(tag//': K', real(k_found))
        call gate%report(tag//': recovered atoms', real(nrec))
        ! the species PDB: one atom per table row, the element column the class's symbol
        call found_atms%new(string(pdb))
        if( l_list )then
            syms = [character(len=2) :: 'PT', upperCase(RUN_EL(ir)(4:5))]
        else
            syms = [character(len=2) :: 'X1', 'X2', 'X3']
        endif
        nwrong = 0
        if( found_atms%get_n() == nrow )then
            do j = 1,nrow
                if( cls(j) < 1 .or. cls(j) > size(syms) )then
                    nwrong = nwrong + 1
                elseif( found_atms%get_element(j) /= syms(cls(j)) )then
                    nwrong = nwrong + 1
                endif
            enddo
        endif
        call gate%check(tag//': species PDB has one atom per table row, element column the class symbol',&
            &found_atms%get_n() == nrow .and. nwrong == 0)
        call found_atms%kill
        ! noise source
        call gate%check(tag//': noise source '//report_string(txt, 'noise_source'),&
            &RUN_HALVES(ir) .eqv. report_string(txt, 'noise_source') == 'half_maps')
        if( RUN_HALVES(ir) )then
            call gate%metric(tag//': half-map label agreement', report_value(txt, 'halfmap_label_agreement'),&
                &MIN_HALF_AGREE, report_value(txt, 'halfmap_label_agreement') >= MIN_HALF_AGREE)
        endif
        if( trim(el2) == '' )then
            ! a single species
            if( l_list )then
                call gate%check(tag//': no recovered atom', nrec == 0)
                call gate%check(tag//': two-class fit inadmissible', report_string(txt, 'admissible_K2') == 'no')
            else
                call gate%check(tag//': K = 1', k_found == 1)
                call gate%check(tag//': no recovered atom', nrec == 0)
            endif
        else
            ! two species: every matched atom carries its generating element (class 1 PT, class 2 the second element)
            nwrong = 0
            do j = 1,nrow
                if( match(j) == 0 ) cycle
                if( cls(j) == 1 .and. gel(match(j)) == ELEMS(1) ) cycle
                if( cls(j) == 2 .and. gel(match(j)) == el2 )   cycle
                nwrong = nwrong + 1
                write(logfhandle,'(a,3f9.3,a,a,a,i0,a,es12.4,a,f7.3,a,i0)') '    wrong element at ', xyz(:,j), ' A: ',&
                    &gel(match(j)), ' in class ', cls(j), ', aperture intensity ', rows(14,j), ', radius ', rows(6,j),&
                    &' A, stage ', nint(rows(8,j))
            enddo
            call gate%report(tag//': atoms with the wrong element', real(nwrong))
            ! the classes in the right order: class 1, the brighter, holds the heavier element
            call gate%check(tag//': classes in the right order (class 1 PT, the brighter)',&
                &count(match > 0 .and. cls == 1 .and. gel(max(match,1)) == ELEMS(1)) >&
                &count(match > 0 .and. cls == 1 .and. gel(max(match,1)) == el2) .and.&
                &count(match > 0 .and. cls == 2 .and. gel(max(match,1)) == el2) >&
                &count(match > 0 .and. cls == 2 .and. gel(max(match,1)) == ELEMS(1)))
            if( l_list )then
                if( trim(CASES(ic)) == 'core' )then
                    ! reported in phase 5, floored again once the atoms are fitted with the map-filtered kernel (ruling
                    ! of 2026-10-09, fifth): the aperture intensity's radial bias can mislabel a Ni atom at the centre
                    call gate%report(tag//': found atoms with the wrong element (reported)', real(nwrong))
                else
                    call gate%check(tag//': every found atom carries its element', nwrong == 0)
                endif
                call gate%check(tag//': two-class fit admissible', report_string(txt, 'admissible_K2') == 'yes')
                call gate%metric(tag//': '//el2//' recall over eligible', recall_of(elig .and. gel == el2, hits),&
                    &MIN_LIGHT_RECALL, recall_of(elig .and. gel == el2, hits) >= MIN_LIGHT_RECALL)
                ratio     = report_value(txt, 'class_ratio_2')
                ratio_gen = int_el(elem_index(el2)) / int_el(1)
                call gate%report(tag//': single-atom intensity ratio', ratio_gen)
                ! the ratio is reported until the atoms are fitted with the map-filtered kernel (ruling of 2026-10-09, second)
                call gate%report(tag//': fitted intensity ratio', ratio)
                call gate%report(tag//': fitted intensity ratio, relative error', abs(ratio / ratio_gen - 1.))
                ! the second element found by stage (0 level 1, 1 residual stage A, 2 stage B) and level
                nfound = 0
                do j = 1,nrow
                    if( match(j) == 0 ) cycle
                    if( gel(match(j)) /= el2 ) cycle
                    nfound(min(max(nint(rows(8,j)),0),2), min(max(nint(rows(9,j)),0),MAXLEV)) = &
                        &nfound(min(max(nint(rows(8,j)),0),2), min(max(nint(rows(9,j)),0),MAXLEV)) + 1
                enddo
                do j1 = 0,2
                    do j2 = 0,MAXLEV
                        if( nfound(j1,j2) == 0 ) cycle
                        call gate%report(tag//': '//el2//' found at stage '//int2str(j1)//', level '//int2str(j2),&
                            &real(nfound(j1,j2)))
                    enddo
                enddo
            else
                call gate%check(tag//': K = 2', k_found == 2)
                call gate%check(tag//': every label right, class 1 PT', nwrong == 0)
            endif
        endif
        ! widths, reported: the stage-2 B of the PT class (class 1) from its inner to its outer radial shell, and the
        ! generated rise over the generating atoms of the same shells. Not floored: on maps filtered as reconstructions
        ! are, the neighbours' negative halos narrow interior atoms, so the fitted rise carries a coordination artefact
        ! (section 7 of the plan)
        rad_csv = trim(RUN_NAME(ir))//'/map_species_radial.csv'
        call read_table(rad_csv, radial)
        sel = pack([(j, j=1,size(radial,2))], nint(radial(1,:)) == 1)
        nsh = size(sel)
        if( nsh < 2 )then
            call gate%check(tag//': the PT class has two radial shells or more', .false.)
            return
        endif
        nin  = nint(radial(3,sel(1)))
        nout = nint(radial(3,sel(nsh)))
        order = pack([(j, j=1,nrow)], cls == 1)
        rsort = rows(6,order)
        call hpsort(rsort, order)
        call gate%check(tag//': radial shells cover the PT class', nint(sum(radial(3,sel))) == size(order))
        rin  = 0.
        rout = 0.
        j1   = 0
        j2   = 0
        do kc = 1,size(order)
            j = order(kc)
            if( match(j) == 0 ) cycle
            if( kc <= nin )then
                rin = rin + sum((gxyz(:,match(j)) - cen)**2) / rmax**2
                j1  = j1 + 1
            elseif( kc > size(order) - nout )then
                rout = rout + sum((gxyz(:,match(j)) - cen)**2) / rmax**2
                j2   = j2 + 1
            endif
        enddo
        rise_gen = B_RISE * (rout / real(max(j2,1)) - rin / real(max(j1,1)))
        rise_fit = radial(7,sel(nsh)) - radial(7,sel(1))
        call gate%report(tag//': PT class B, inner shell', radial(7,sel(1)))
        call gate%report(tag//': PT class B, outer shell', radial(7,sel(nsh)))
        call gate%report(tag//': generated rise of B', rise_gen)
        call gate%report(tag//': fitted over generated rise of B', rise_fit / rise_gen)
    end subroutine evaluate

    ! light case: element=Pt must miss Al atoms that the list run finds, or the fixture does not test the recovery
    subroutine light_plain_recall( ic )
        integer, intent(in) :: ic
        type(atoms) :: found_atms
        integer, allocatable :: hits(:)
        real    :: xyz(3), rec
        integer :: j, m
        if( trim(CASES(ic)) /= 'light' ) return
        call found_atms%new(string('pt/map_ATMS.pdb'))
        allocate(hits(nlat), source=0)
        do j = 1,found_atms%get_n()
            xyz = found_atms%get_coord(j)
            m   = minloc(sum((gxyz - spread(xyz, 2, nlat))**2, dim=1), dim=1)
            if( norm2(gxyz(:,m) - xyz) <= MATCH_NN_FRAC * dnn_gen ) hits(m) = hits(m) + 1
        enddo
        call found_atms%kill
        rec = recall_of(elig .and. gel == SECOND(ic), hits)
        call gate%metric('light/pt: AL recall over eligible below the ceiling (the fixture tests the recovery)', rec,&
            &MAX_PLAIN_LIGHT, rec < MAX_PLAIN_LIGHT)
    end subroutine light_plain_recall

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
