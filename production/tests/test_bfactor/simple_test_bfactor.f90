! Tests anisotropic and isotropic bfactor calculations in atoms_stats
! Requires an input reconstructed volume, connected component binary volume
! and a PDB file.
program simple_test_bfactor
    include 'simple_lib.f08'
    use simple_cmdline,            only: cmdline
    !use simple_builder,            only: builder
    !use simple_parameters,         only: params_glob
    use simple_commander_atoms,    only: atoms_stats_commander
    !use simple_image,              only: image
    !use simple_binimage,           only: binimage
    use simple_atoms,              only: atoms
    implicit none

    ! Output atomic statistics in atoms_stats.csv file. 
    type :: stats
        real :: index, nvox, cnstd, nnbondl, cngen, diam, avgint, maxint
        real :: cendist, vcorr, displ, maxndispl 
        real :: bfac, iso_corr, aniso_corr, semiax(3) ! Relevant
        real :: anisoxyz(3)
        real :: azimuth, polar ! Relevant
        real :: rstrain
    end type

    type (cmdline)              :: cline_atoms_stats
    type(atoms)                 :: inatoms
    type(atoms_stats_commander) :: xatoms_stats
    type(stats), allocatable    :: ideal(:), fitiso(:), fitaniso(:)
    real, parameter             :: smpd=0.358, mskdiam=40., nthr=16.
    integer                     :: natoms, fu_ideal, fu_iso, fu_aniso, cc, io
    integer                     :: fu_test1, fu_test2
    character(len=256)          :: invol, incc, inpdb
    character(len=2)            :: element
    character(len=256)          :: junk(7)
    character(*), parameter     :: dir='simple_test_bfactor_out'
    character(*), parameter     :: isodir='test_isotropic_bfactor'
    character(*), parameter     :: anisodir='test_anisotropic_bfactor'
    character(*), parameter     :: fn_stats='atoms_stats.csv'
    character(*), parameter     :: fn_isofit='fit_isotropic.mrc'
    character(*), parameter     :: fn_anisofit='fit_anisotropic.mrc'
    character(*), parameter     :: fn_test1='results_simple_test_bfactor_iso_input.csv'
    character(*), parameter     :: fn_test2='results_simple_test_bfactor_aniso_input.csv'
    character(*), parameter :: HEADTEST1 = 'INDEX'//CSV_DELIM//'BFAC_IN'//CSV_DELIM&
        &//'BFAC_OUT'//CSV_DELIM//'MAJSAX_OUT'//CSV_DELIM//'MEDSAX_OUT'//CSV_DELIM//&
        &'MINSAX_OUT'//CSV_DELIM//'AZIMUTH_OUT'//CSV_DELIM//'POLAR_OUT'//CSV_DELIM//&
        &'ISOCORR'//CSV_DELIM//'ANISO_CORR'
    character(*), parameter :: HEADTEST2 = 'INDEX'//CSV_DELIM//'MAJSAX_IN'//&
        &CSV_DELIM//'MEDSAX_IN'//CSV_DELIM//'MINSAX_IN'//CSV_DELIM//'AZIMUTH_IN'//&
        &'POLAR_INT'//CSV_DELIM//'BFAC_OUT'//CSV_DELIM//'MAJSAX_OUT'//CSV_DELIM//&
        &'MEDSAX_OUT'//CSV_DELIM//'MINSAX_OUT'//CSV_DELIM//'AZIMUTH_OUT'//CSV_DELIM//&
        &'POLAR_OUT'//CSV_DELIM//'ISOCORR'//CSV_DELIM//'ANISO_CORR'


    if( command_argument_count() /= 3 )then
        write(logfhandle,'(a)') 'Usage: simple_test_bfactor vol.mrc cc.mrc &
                                &file.pdb'
        stop
    else
        call get_command_argument(1, invol)
        call get_command_argument(2, incc)
        call get_command_argument(3, inpdb)
    endif

    call inatoms%new(inpdb)
    element = inatoms%get_element(1)
    natoms = inatoms%get_n()

    ! Generate ideal isotropic and anisotropic Gaussian fits from atoms stats
    write(logfhandle,'(a)') '>>> PASSING INPUT FILES TO ATOMS_STATS'
    call simple_mkdir(dir)
    call simple_chdir(dir)
    call cline_atoms_stats%set('vol1', '../'//invol)
    call cline_atoms_stats%set('vol2', '../'//incc)
    call cline_atoms_stats%set('pdbfile', '../'//inpdb)
    call cline_atoms_stats%set('smpd', smpd)
    call cline_atoms_stats%set('element', element)
    call cline_atoms_stats%set('mskdiam', mskdiam)
    call cline_atoms_stats%set('nthr', nthr)
    call xatoms_stats%execute(cline_atoms_stats)

    ! Test 1: Use a density with ideal isotropic 3D Gaussians
    write(logfhandle,'(a)') '>>> TEST 1: IDEAL ISOTROPIC 3D GAUSSIANS'
    call simple_mkdir(isodir)
    call simple_chdir(isodir)
    call cline_atoms_stats%set('vol1', '../'//fn_isofit)
    call cline_atoms_stats%set('vol2', '../../'//incc)
    call cline_atoms_stats%set('pdbfile', '../../'//inpdb)
    call cline_atoms_stats%set('smpd', smpd)
    call cline_atoms_stats%set('element', element)
    call cline_atoms_stats%set('mskdiam', mskdiam)
    call cline_atoms_stats%set('nthr', nthr)
    call xatoms_stats%execute(cline_atoms_stats)
    call simple_chdir('../')

    ! Test 2: Use a density with ideal anisotropic 3D Gaussians
    write(logfhandle,'(a)') '>>> TEST 2: IDEAL ANISOTROPIC 3D GAUSSIANS'
    call simple_mkdir(anisodir)
    call simple_chdir(anisodir)
    call cline_atoms_stats%set('vol1', '../'//fn_anisofit)
    call cline_atoms_stats%set('vol2', '../../'//incc)
    call cline_atoms_stats%set('pdbfile', '../../'//inpdb)
    call cline_atoms_stats%set('smpd', smpd)
    call cline_atoms_stats%set('element', element)
    call cline_atoms_stats%set('mskdiam', mskdiam)
    call cline_atoms_stats%set('nthr', nthr)
    call xatoms_stats%execute(cline_atoms_stats)
    call simple_chdir('../')

    ! Get the ideal and fit stats: isotropic Bfactors, anisotropic Bfactors,
    ! and anisotropic major semi-axis orientations
    ! Last dimension of stats is test number.
    allocate(ideal(natoms), fitiso(natoms), fitaniso(natoms))
    call fopen(fu_ideal, FILE=trim(fn_stats), STATUS='OLD', action='READ')
    call fopen(fu_iso, FILE=trim(isodir//'/'//fn_stats), STATUS='OLD', action='READ')
    call fopen(fu_aniso, FILE=trim(anisodir//'/'//fn_stats), STATUS='OLD', action='READ')
    call fopen(fu_test1, FILE=trim(fn_test1), STATUS='REPLACE', action='WRITE', iostat=io)
    call fopen(fu_test2, FILE=trim(fn_test2), STATUS='REPLACE', action='WRITE')
    read (fu_ideal, *)      ! Skip header
    read (fu_iso, *)        ! Skip header
    read (fu_aniso, *)      ! Skip header
    write (fu_test1, '(A)') HEADTEST1
    write (fu_test2, '(A)') HEADTEST2
    do cc=1,natoms
        read (fu_ideal, *) ideal(cc)
        read (fu_iso, *) fitiso(cc)
        read(fu_aniso, *) fitaniso(cc)
        ! Write results for both tests
        call write_test(cc, fu_test1, ideal, fitiso, .true.)
        call write_test(cc, fu_test2, ideal, fitaniso, .false.)
    end do
    call fclose(fu_ideal)
    call fclose(fu_iso)
    call fclose(fu_aniso)
    call fclose(fu_test1)
    call fclose(fu_test2)
    deallocate(ideal, fitiso, fitaniso)
    write(logfhandle, '(a)') '***END SIMPLE_TEST_BFACTOR***'

contains

    subroutine write_test(cc, fu, stats_in, stats_out, iso_in)
        integer, intent(in) :: cc, fu
        type(stats), allocatable, intent(in) :: stats_in(:), stats_out(:)
        logical, intent(in) :: iso_in

        601 format(F8.4,A2)
        602 format(F8.4)
        ! various per-atom parameters
        write(fu,601,advance='no') real(cc),                    CSV_DELIM ! INDEX
        if (iso_in) then
            write(fu,601,advance='no') stats_in(cc)%bfac,       CSV_DELIM ! BFAC_IN
        else
            write(fu,601,advance='no') stats_in(cc)%semiax(1),  CSV_DELIM ! MAJSAX_IN
            write(fu,601,advance='no') stats_in(cc)%semiax(2),  CSV_DELIM ! MEDSAX_IN
            write(fu,601,advance='no') stats_in(cc)%semiax(3),  CSV_DELIM ! MINSAX_IN
            write(fu,601,advance='no') stats_in(cc)%azimuth,    CSV_DELIM ! AZMTH_IN
            write(fu,601,advance='no') stats_in(cc)%polar,      CSV_DELIM ! POLAR_IN
        end if
        write(fu,601,advance='no') stats_out(cc)%bfac,          CSV_DELIM ! BFAC_OUT
        write(fu,601,advance='no') stats_out(cc)%semiax(1),     CSV_DELIM ! MAJSAX_OUT
        write(fu,601,advance='no') stats_out(cc)%semiax(2),     CSV_DELIM ! MEDSAX_OUT
        write(fu,601,advance='no') stats_out(cc)%semiax(3),     CSV_DELIM ! MINSAX_OUT
        write(fu,601,advance='no') stats_out(cc)%azimuth,       CSV_DELIM ! AZIMUTH_OUT
        write(fu,601,advance='no') stats_out(cc)%polar,         CSV_DELIM ! POLAR_OUT
        write(fu,601,advance='no') stats_out(cc)%iso_corr,      CSV_DELIM ! ISO_CORR
        write(fu,602)              stats_out(cc)%aniso_corr               ! ANISO_CORR
    end subroutine write_test


end program simple_test_bfactor