!@descr: unit tests for atomic models (simple_atoms): access, geometry, PDB I/O, density simulation and per-atom validation
! Replaces the in-module self-test test_atoms of simple_atoms (its private assertions stopped at the
! first failure, the I/O and density parts only ran, and three PDB files were left behind); the
! routines only that self-test called were removed (plan, section 9.7, the open items, 2026-09-25).
! atom_validate is pinned against a map simulated from the model itself: every atom must correlate
! with its own simulated density (0.99 in a numpy emulation of convolve and the window mask), which
! the one-voxel offset of its window (fixed 2026-09-25) broke (0.35-0.44 in the same emulation).
module simple_atoms_tester
use simple_test_utils
use simple_defs
use simple_string, only: string
use simple_syslib, only: del_file
use simple_image,  only: image, unmemoize_mask_coords
use simple_atoms,  only: atoms
implicit none
private
public :: run_all_atoms_tests

character(len=*), parameter :: TMP_PDB   = 'tmp_atoms_tester.pdb'
character(len=*), parameter :: TMP_ANISO = 'tmp_atoms_tester_aniso.pdb'

contains

    subroutine run_all_atoms_tests()
        write(*,'(A)') '**** running all atoms tests ****'
        call test_access()
        call test_geometry()
        call test_pdb_roundtrip()
        call test_anisou()
        call test_validation()
    end subroutine run_all_atoms_tests

    ! three dummy atoms named C, O and N
    subroutine make_model( a )
        type(atoms), intent(inout) :: a
        call a%new(3, dummy=.true.)
        call a%set_name(1, ' C  ')
        call a%set_name(2, ' O  ')
        call a%set_name(3, ' N  ')
        call a%set_coord(1, [1.25, -2.50,  3.75])
        call a%set_coord(2, [0.50,  4.00, -1.00])
        call a%set_coord(3, [-3.0,  1.50,  0.25])
        call a%set_resnum(1, 7)
        call a%set_resnum(2, 7)
        call a%set_resnum(3, 8)
        call a%guess_element()
    end subroutine make_model

    subroutine test_access()
        type(atoms) :: a, b, c, one
        write(*,'(A)') 'test_access'
        call make_model(a)
        call assert_int(3, a%get_n(), 'new(3) holds three atoms')
        call assert_true(all(a%get_coord(1) == [1.25, -2.50, 3.75]), 'set_coord/get_coord round trip')
        call assert_int(8, a%get_resnum(3), 'set_resnum/get_resnum round trip')
        call a%set_beta(2, 12.5)
        call assert_real(12.5, a%get_beta(2), 0., 'set_beta/get_beta round trip')
        call assert_char(' O  ', a%get_name(2), 'set_name/get_name round trip')
        call assert_char('C ', a%get_element(1), 'guess_element: C from the name')
        call assert_char('O ', a%get_element(2), 'guess_element: O from the name')
        call assert_char('N ', a%get_element(3), 'guess_element: N from the name')
        call assert_int(6, a%get_atomicnumber(1), 'guess_element: Z of carbon')
        call assert_true(a%element_exists('C'),  'element_exists: C')
        call assert_true(a%element_exists('Pd'), 'element_exists: Pd')
        call b%copy(a)
        call assert_int(a%get_n(), b%get_n(), 'copy keeps the atom count')
        call assert_true(all(b%get_coord(3) == a%get_coord(3)), 'copy keeps the coordinates')
        c = a
        call assert_true(all(c%get_coord(2) == a%get_coord(2)), 'assignment keeps the coordinates')
        call a%extract_atom(one, 2)
        call assert_int(1, one%get_n(), 'extract_atom gives one atom')
        call assert_true(all(one%get_coord(1) == a%get_coord(2)), 'extract_atom keeps the coordinates')
        call assert_char('O ', one%get_element(1), 'extract_atom keeps the element')
        call a%kill
        call assert_int(0, a%get_n(), 'kill leaves no atoms')
        call b%kill
        call c%kill
        call one%kill
    end subroutine test_access

    subroutine test_geometry()
        real, parameter :: CEN(3) = [(1.25+0.50-3.0)/3., (-2.50+4.00+1.50)/3., (3.75-1.00+0.25)/3.]
        real, parameter :: SMPD = 1.5
        type(atoms) :: a
        write(*,'(A)') 'test_geometry'
        call make_model(a)
        call assert_true(all(abs(a%get_geom_center() - CEN) < 1.e-6), 'get_geom_center is the mean position')
        call a%translate([1., 2., 3.])
        call assert_true(all(abs(a%get_geom_center() - (CEN + [1., 2., 3.])) < 1.e-6), 'translate moves the centre by the shift')
        call a%center_pdbcoord([33,33,33], SMPD)
        call assert_true(all(abs(a%get_geom_center() - 16. * SMPD) < 1.e-5), 'center_pdbcoord puts the centre at (ldim-1)/2*smpd')
        call a%center_inbox(2, 20, SMPD)
        call assert_true(all(abs(a%get_coord(2) - 10. * SMPD) < 1.e-6), 'center_inbox puts the atom at box/2*smpd')
        call a%kill
    end subroutine test_geometry

    ! PDB coordinates have three decimals
    subroutine test_pdb_roundtrip()
        type(atoms) :: a, b
        integer     :: i
        logical     :: ok
        write(*,'(A)') 'test_pdb_roundtrip'
        call make_model(a)
        call del_file(TMP_PDB)
        call a%writepdb(string(TMP_PDB))
        call b%new(string(TMP_PDB))
        call assert_int(a%get_n(), b%get_n(), 'a PDB round trip keeps the atom count')
        ok = .true.
        do i = 1,a%get_n()
            if( any(abs(b%get_coord(i) - a%get_coord(i)) > 1.e-3) ) ok = .false.
            if( b%get_atomicnumber(i) /= a%get_atomicnumber(i) )    ok = .false.
            if( b%get_resnum(i)  /= a%get_resnum(i) )               ok = .false.
        end do
        call assert_true(ok, 'a PDB round trip keeps coordinates (to 1e-3 A), atomic numbers and residue numbers')
        call del_file(TMP_PDB)
        call a%kill
        call b%kill
    end subroutine test_pdb_roundtrip

    ! one ANISOU record per atom, U in units of 1e-4 A**2 in columns 29-70 (PDB format)
    subroutine test_anisou()
        type(atoms)       :: a
        real, allocatable :: aniso(:,:,:)
        character(len=6)  :: tag
        character(len=128):: line
        integer :: i, funit, ios, natom, naniso, u(6)
        logical :: ok
        write(*,'(A)') 'test_anisou'
        call make_model(a)
        allocate(aniso(3,3,3), source=0.)
        do i = 1,3
            aniso(1,1,i) = 1.e-3 * i
            aniso(2,2,i) = 2.e-3 * i
            aniso(3,3,i) = 3.e-3 * i
            aniso(1,2,i) = 5.e-4
        end do
        call del_file(TMP_ANISO)
        call a%writepdb_aniso(string(TMP_ANISO), aniso)
        open(newunit=funit, file=TMP_ANISO, status='old', action='read', iostat=ios)
        call assert_int(0, ios, 'the ANISOU file opens')
        natom  = 0
        naniso = 0
        ok     = .true.
        if( ios == 0 )then
            do
                read(funit, '(A)', iostat=ios) line
                if( ios /= 0 ) exit
                tag = line(1:6)
                if( tag == 'ATOM  ' .or. tag == 'HETATM' ) natom = natom + 1
                if( tag == 'ANISOU' )then
                    naniso = naniso + 1
                    read(line, '(28X,6I7)') u
                    if( any(u /= [10*naniso, 20*naniso, 30*naniso, 5, 0, 0]) ) ok = .false.
                endif
            end do
            close(funit)
        endif
        call assert_int(3, natom,  'writepdb_aniso writes one coordinate record per atom')
        call assert_int(3, naniso, 'writepdb_aniso writes one ANISOU record per atom')
        call assert_true(ok, 'ANISOU records hold U11 U22 U33 U12 U13 U23 in units of 1e-4 A**2')
        call del_file(TMP_ANISO)
        call a%kill
    end subroutine test_anisou

    ! three atoms 8 A apart on grid points of a 1 A map simulated from the model: map_validate of the
    ! map against itself gives 1 for every atom; atom_validate correlates each atom's simulated density
    ! with its window of the map, which holds that atom alone at the same place
    subroutine test_validation()
        integer, parameter :: B3 = 32
        real,    parameter :: SMPD = 1.0
        type(atoms) :: a
        type(image) :: vol, vol2
        integer     :: i
        real        :: beta_min
        logical     :: ok
        write(*,'(A)') 'test_validation'
        call make_model(a)
        call a%set_coord(1, [ 8., 16., 16.])
        call a%set_coord(2, [16., 16., 16.])
        call a%set_coord(3, [24., 16., 16.])
        call vol%new([B3,B3,B3], SMPD, wthreads=.false.)
        call a%convolve(vol, cutoff=8.*SMPD)
        call assert_true(vol%get_rmat_at(17,17,17) > vol%get_rmat_at(13,17,17), 'convolve: density peaks at the atom')
        call vol2%copy(vol)
        call a%map_validate(vol, vol2)
        ok = .true.
        do i = 1,a%get_n()
            if( abs(a%get_beta(i) - 1.) > 1.e-4 ) ok = .false.
        end do
        call assert_true(ok, 'map_validate of a map against itself scores every atom 1')
        call a%atom_validate(vol)
        beta_min = minval([(a%get_beta(i), i=1,a%get_n())])
        write(logfhandle,'(A,F8.4)') 'atom_validate, lowest atom correlation: ', beta_min
        call assert_true(beta_min > 0.95, 'atom_validate: every atom correlates with its own simulated density')
        call a%kill
        call vol%kill
        call vol2%kill
        call unmemoize_mask_coords ! atom_validate memoised the mask of its 4-voxel windows
    end subroutine test_validation

end module simple_atoms_tester
