!@descr: unit tests for the C-alpha candidate search in density maps (simple_calpha_finder)
! Moved from the detect_calpha case of simple_test_exec (Ruben's, 2026-09-23) as it was: a 32^3 map
! of three Gaussian residues (CA, N, C sites, the CA brighter) along the diagonal, searched with a
! 180 degree step for ten candidates. The checks are the original ones and are weak (some candidate
! within 1.5 A of some residue centre); doc/refactoring_notes/single_area_tests_handover.md says
! what they should pin.
module simple_calpha_finder_tester
use simple_core_module_api
use simple_image,         only: image
use simple_atoms,         only: atoms
use simple_calpha_finder, only: calpha_finder
use simple_test_utils
implicit none
private
public :: run_all_calpha_finder_tests

contains

    subroutine run_all_calpha_finder_tests()
        write(*,'(A)') '**** running all C-alpha finder tests ****'
        call test_three_residues()
    end subroutine run_all_calpha_finder_tests

    subroutine test_three_residues()
        character(len=*), parameter :: PDB_FILE = 'test_calpha_candidates.pdb'
        character(len=*), parameter :: CSV_FILE = 'test_calpha_candidates.csv'
        character(len=*), parameter :: MRC_FILE = 'test_calpha_scores.mrc'
        integer,          parameter :: TEST_BOX = 32, NRES = 3
        type(image)         :: workvol
        type(atoms)         :: candidates
        type(calpha_finder) :: finder
        real(kind=c_float), pointer :: density(:,:,:)
        real    :: centers(3,NRES), atom_sites(3,3), amplitudes(3), rotation(3,3)
        real    :: xyz(3), site(3), delta(3), distance, closest
        integer :: ldim(3), ires, iatom, ix, iy, iz, ncand
        write(*,'(A)') 'test_three_residues'
        ldim            = [TEST_BOX,TEST_BOX,TEST_BOX]
        centers(:,1)    = [8.,8.,8.]
        centers(:,2)    = [16.,16.,16.]
        centers(:,3)    = [24.,24.,24.]
        atom_sites(:,1) = [0.,0.,0.]
        atom_sites(:,2) = [1.458*cos(111.2*PI/180.), 1.458*sin(111.2*PI/180.), 0.]
        atom_sites(:,3) = [1.525,0.,0.]
        amplitudes      = [1.25,1.0,1.0]
        rotation(:,1)   = [0.,1.,0.]
        rotation(:,2)   = [-0.5,0.,sqrt(0.75)]
        rotation(:,3)   = [sqrt(0.75),0.,0.5]
        call workvol%new(ldim, 1.0)
        call workvol%get_rmat_ptr(density)
        density = 0.
        do ires = 1, NRES
            do iatom = 1, size(atom_sites,2)
                site = centers(:,ires) + matmul(rotation, atom_sites(:,iatom))
                do iz = 1, TEST_BOX
                    do iy = 1, TEST_BOX
                        do ix = 1, TEST_BOX
                            xyz   = real([ix,iy,iz] - 1)
                            delta = xyz - site
                            density(ix,iy,iz) = density(ix,iy,iz) + amplitudes(iatom) * &
                                exp(-0.5 * sum(delta * delta) / 0.85**2)
                        enddo
                    enddo
                enddo
            enddo
        enddo
        call finder%new(1.0, 4.0)
        call finder%search(workvol, 180.0, 10, 0.5, string(PDB_FILE), string(MRC_FILE))
        ncand = nlines(string(PDB_FILE))
        call assert_true(ncand > 0, 'the finder writes candidates for three residues')
        if( ncand > 0 )then
            call candidates%new(string(PDB_FILE))
            closest = huge(1.)
            do iatom = 1, candidates%get_n()
                do ires = 1, NRES
                    distance = sqrt(sum((candidates%get_coord(iatom) - centers(:,ires))**2))
                    closest  = min(closest, distance)
                enddo
            enddo
            call assert_true(closest <= 1.5, 'a candidate lies within 1.5 A of a residue centre')
            write(*,'(A,F7.3,A)') '  closest candidate to a residue centre: ', closest, ' A'
            call candidates%kill()
        endif
        call finder%kill()
        call workvol%kill()
        if( file_exists(PDB_FILE) ) call del_file(PDB_FILE)
        if( file_exists(CSV_FILE) ) call del_file(CSV_FILE)
        if( file_exists(MRC_FILE) ) call del_file(MRC_FILE)
    end subroutine test_three_residues

end module simple_calpha_finder_tester
