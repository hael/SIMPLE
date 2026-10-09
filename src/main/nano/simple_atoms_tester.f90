!@descr: unit tests for atomic models (simple_atoms): access, geometry, PDB I/O, density simulation, pseudo-atoms and per-atom validation
! atom_validate is pinned against a map simulated from the model itself: every atom must correlate
! with its own simulated density. Guards a one-voxel offset of the atom_validate window.
module simple_atoms_tester
use simple_test_utils
use simple_defs
use simple_defs_atoms, only: Z_PSEUDO_FIRST, Z_PSEUDO_LAST
use simple_string, only: string
use simple_syslib, only: del_file
use simple_image,  only: image, unmemoize_mask_coords
use simple_atoms,  only: atoms
!$ use omp_lib,    only: omp_get_max_threads, omp_set_num_threads
implicit none
private
public :: run_all_atoms_tests

character(len=*), parameter :: TMP_PDB   = 'tmp_atoms_tester.pdb'
character(len=*), parameter :: TMP_ANISO = 'tmp_atoms_tester_aniso.pdb'
character(len=*), parameter :: TMP_PSEUDO = 'tmp_atoms_tester_pseudo.pdb'

contains

    subroutine run_all_atoms_tests()
        write(*,'(A)') '**** running all atoms tests ****'
        call test_access()
        call test_geometry()
        call test_pdb_roundtrip()
        call test_anisou()
        call test_validation()
        call test_pseudo_symbols()
        call test_pseudo_density()
        call test_pseudo_pdb_roundtrip()
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
        integer     :: i, nthr_saved
        real        :: beta_min
        logical     :: ok
        write(*,'(A)') 'test_validation'
        call make_model(a)
        call a%set_coord(1, [ 8., 16., 16.])
        call a%set_coord(2, [16., 16., 16.])
        call a%set_coord(3, [24., 16., 16.])
        call vol%new([B3,B3,B3], SMPD, wthreads=.false.)
        call vol2%new([B3,B3,B3], SMPD, wthreads=.false.)
        nthr_saved = 1
        !$ nthr_saved = omp_get_max_threads()
        !$ call omp_set_num_threads(1)
        call a%convolve(vol, cutoff=8.*SMPD)
        !$ call omp_set_num_threads(3)
        call a%convolve(vol2, cutoff=8.*SMPD)
        !$ call omp_set_num_threads(nthr_saved)
        call assert_true(all(vol%get_rmat() == vol2%get_rmat()), 'convolve gives the same voxels with one or three threads')
        call assert_true(vol%get_rmat_at(17,17,17) > vol%get_rmat_at(13,17,17), 'convolve: density peaks at the atom')
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

    ! X1..X3 are known symbols with distinct sentinel Z in the reserved range, clear of CDSE's 999
    subroutine test_pseudo_symbols()
        type(atoms) :: a
        integer     :: i, z(3)
        write(*,'(A)') 'test_pseudo_symbols'
        call assert_true(a%element_exists('X1') .and. a%element_exists('X2') .and. a%element_exists('X3'),&
            &'element_exists: X1, X2 and X3')
        call a%new(3, dummy=.true.)
        call a%set_element(1, 'X1')
        call a%set_element(2, 'x2')
        call a%set_name(3, 'X3  ')
        do i = 1,3
            z(i) = a%get_atomicnumber(i)
        enddo
        call assert_true(all(z >= Z_PSEUDO_FIRST .and. z <= Z_PSEUDO_LAST), 'pseudo-atoms have Z in the reserved range')
        call assert_true(z(1) /= z(2) .and. z(2) /= z(3) .and. z(1) /= z(3), 'X1, X2 and X3 have distinct Z')
        call assert_true(Z_PSEUDO_LAST < 999 .or. Z_PSEUDO_FIRST > 999, 'the reserved range excludes the CDSE sentinel 999')
        call assert_char('X2', a%get_element(2), 'set_element upper-cases the pseudo-atom symbol')
        call assert_char('X3', a%get_element(3), 'set_name sets the pseudo-atom symbol')
        call assert_real(1., a%get_radius(1), 0., 'pseudo-atom radius is 1 A')
        call a%kill
    end subroutine test_pseudo_symbols

    ! one pseudo-atom of integrated intensity Q and width B: peak Q (4 pi / B)**1.5 at the atom, voxel sum
    ! Q / smpd**3 within 1% at a cutoff of 4 sigma; two atoms add linearly; lp does not blur a pseudo-atom
    subroutine test_pseudo_density()
        integer, parameter :: BX = 48
        real,    parameter :: SMPD = 0.358, B1 = 13.9, Q1 = 2.0, B2 = 20.0, Q2 = 0.5
        type(atoms)        :: a, a1, a2
        type(image)        :: vol, vol1, vol2
        real, allocatable  :: r(:,:,:), r1(:,:,:), r2(:,:,:)
        real :: sigma, cutoff, peak
        write(*,'(A)') 'test_pseudo_density'
        sigma  = sqrt(B1 / (8. * PI**2))
        cutoff = 4. * sigma
        call vol%new([BX,BX,BX], SMPD, wthreads=.false.)
        call vol1%new([BX,BX,BX], SMPD, wthreads=.false.)
        call vol2%new([BX,BX,BX], SMPD, wthreads=.false.)
        ! on the voxel (25,25,25): voxel j is at (j - 1) * smpd
        call make_pseudo(a1, [24., 24., 24.] * SMPD, Q1, B1)
        call a1%convolve(vol1, cutoff)
        peak = Q1 * (4. * PI / B1)**1.5
        call assert_real(peak, vol1%get_rmat_at(25,25,25), 1.e-4 * peak, 'pseudo-atom peak is q (4 pi / B)**1.5')
        ! (29,25,25) is 1.43 A from the atom, (29,29,25) 2.03 A: inside and outside the 4 sigma = 1.68 A cutoff
        call assert_true(vol1%get_rmat_at(29,25,25) > 0., 'pseudo-atom density reaches the cutoff')
        call assert_real(0., vol1%get_rmat_at(29,29,25), 0., 'pseudo-atom density stops at the cutoff')
        ! off the grid
        call make_pseudo(a2, [23.3, 24.6, 25.2] * SMPD, Q1, B1)
        call a2%convolve(vol2, cutoff)
        call assert_real(Q1, sum(vol2%get_rmat()) * SMPD**3, 0.01 * Q1, 'pseudo-atom voxel sum is q / smpd**3 within 1% at 4 sigma')
        call a2%convolve(vol, cutoff, lp=1.5)
        call assert_true(all(vol%get_rmat() == vol2%get_rmat()), 'lp does not blur a pseudo-atom')
        ! two overlapping atoms 2 A apart
        call make_pseudo(a2, [25.0, 25.0, 25.0] * SMPD + [2., 0., 0.], Q2, B2)
        call a2%convolve(vol2, 6. * sigma)
        call a1%convolve(vol1, 6. * sigma)
        call a%new(2, dummy=.true.)
        call a%set_element(1, 'X1')
        call a%set_element(2, 'X2')
        call a%set_coord(1, a1%get_coord(1))
        call a%set_coord(2, a2%get_coord(1))
        call a%set_occupancy(1, Q1)
        call a%set_occupancy(2, Q2)
        call a%set_beta(1, B1)
        call a%set_beta(2, B2)
        call a%convolve(vol, 6. * sigma)
        r  = vol%get_rmat()
        r1 = vol1%get_rmat()
        r2 = vol2%get_rmat()
        call assert_true(all(abs(r - (r1 + r2)) <= 1.e-6 * maxval(r)), 'two pseudo-atoms add linearly')
        call assert_real(Q1 + Q2, sum(r) * SMPD**3, 0.01 * (Q1 + Q2), 'two pseudo-atoms integrate to q1 + q2')
        call a%kill
        call a1%kill
        call a2%kill
        call vol%kill
        call vol1%kill
        call vol2%kill
    end subroutine test_pseudo_density

    subroutine make_pseudo( a, xyz, q, bfac )
        type(atoms), intent(inout) :: a
        real,        intent(in)    :: xyz(3), q, bfac
        call a%new(1, dummy=.true.)
        call a%set_element(1, 'X1')
        call a%set_coord(1, xyz)
        call a%set_occupancy(1, q)
        call a%set_beta(1, bfac)
    end subroutine make_pseudo

    ! symbol, occupancy (q) and beta (B) survive writepdb and new; the PDB columns hold two decimals
    subroutine test_pseudo_pdb_roundtrip()
        character(len=2), parameter :: SYMS(3) = ['X1', 'X2', 'X3']
        real,             parameter :: QS(3) = [1.00, 0.17, 0.50], BS(3) = [13.90, 20.25, 7.50]
        type(atoms) :: a, b
        integer     :: i
        logical     :: ok
        write(*,'(A)') 'test_pseudo_pdb_roundtrip'
        call a%new(3, dummy=.true.)
        do i = 1,3
            call a%set_element(i, SYMS(i))
            call a%set_coord(i, [1.5 * i, -2.25, 3.125])
            call a%set_occupancy(i, QS(i))
            call a%set_beta(i, BS(i))
        enddo
        call del_file(TMP_PSEUDO)
        call a%writepdb(string(TMP_PSEUDO))
        call b%new(string(TMP_PSEUDO))
        call assert_int(3, b%get_n(), 'pseudo-atom PDB round trip keeps the atom count')
        ok = .true.
        do i = 1,3
            if( b%get_element(i) /= SYMS(i) )                     ok = .false.
            if( b%get_atomicnumber(i) /= a%get_atomicnumber(i) )  ok = .false.
            if( abs(b%get_occupancy(i) - QS(i)) > 0.005 )         ok = .false.
            if( abs(b%get_beta(i) - BS(i)) > 0.005 )              ok = .false.
            if( any(abs(b%get_coord(i) - a%get_coord(i)) > 1.e-3) ) ok = .false.
        enddo
        call assert_true(ok, 'pseudo-atom PDB round trip keeps symbol, Z, occupancy and B (to 0.005) and coordinates')
        call del_file(TMP_PSEUDO)
        call a%kill
        call b%kill
    end subroutine test_pseudo_pdb_roundtrip

end module simple_atoms_tester
