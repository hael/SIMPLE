!@descr: library tests of pdb2mrc (simple_atoms): density maps from the built-in 6VXX and 1JYX models
! Each model is converted with one and three threads; the resulting MRC files must be byte-identical.
module simple_pdb2mrc_tester
use simple_string,        only: string
use simple_syslib,        only: del_file, file_exists
use simple_imghead,       only: find_img_smpd, find_ldim_nptcls
use simple_atoms,         only: atoms
use simple_molecule_data, only: molecule_data, betagal_1jyx, sars_cov2_spkgp_6vxx
use simple_fileio,        only: files_identical
use simple_test_utils,    only: assert_real, assert_true
!$ use omp_lib,           only: omp_get_max_threads, omp_set_num_threads
implicit none
private
public :: run_all_pdb2mrc_tests

real, parameter :: SMPD = 1.3

contains

    subroutine run_all_pdb2mrc_tests()
        write(*,'(A)') '**** running all pdb2mrc tests ****'
        call test_pdb2mrc_model('6VXX', sars_cov2_spkgp_6vxx())
        call test_pdb2mrc_model('1JYX', betagal_1jyx())
    end subroutine run_all_pdb2mrc_tests

    subroutine test_pdb2mrc_model( code, mol )
        character(len=*),    intent(in) :: code
        type(molecule_data), intent(in) :: mol
        type(atoms)  :: molecule
        type(string) :: pdb_file, vol_file
        integer      :: nthr_saved
        write(*,'(A)') 'test_pdb2mrc_model '//code
        call assert_true(mol%n > 0, code//': the built-in model has atoms')
        pdb_file = 'molecule.pdb'
        vol_file = code//'.mrc'
        call molecule%new(mol)
        call molecule%writepdb(pdb_file)
        call molecule%kill
        nthr_saved = 1
        !$ nthr_saved = omp_get_max_threads()
        !$ call omp_set_num_threads(1)
        call molecule%pdb2mrc(pdbfile=pdb_file, smpd=SMPD)
        call check_map(string('molecule.mrc'), code//' (one thread)')
        call molecule%kill
        !$ call omp_set_num_threads(3)
        call molecule%pdb2mrc(pdbfile=pdb_file, volfile=vol_file, smpd=SMPD)
        !$ call omp_set_num_threads(nthr_saved)
        call check_map(vol_file, code//' (three threads)')
        call assert_true(files_identical(string('molecule.mrc'), vol_file),&
            &code//': pdb2mrc output is byte-identical with one or three threads')
        call molecule%kill
        call del_file(string('molecule.pdb'))
        call del_file(string('molecule_centered.pdb'))
        call del_file(string('molecule.mrc'))
        call del_file(vol_file)
    end subroutine test_pdb2mrc_model

    subroutine check_map( vol_file, what )
        type(string),     intent(in) :: vol_file
        character(len=*), intent(in) :: what
        integer :: ldim(3), nptcls
        call assert_true(file_exists(vol_file), what//': the map is written')
        if( .not. file_exists(vol_file) ) return
        call find_ldim_nptcls(vol_file, ldim, nptcls)
        call assert_true(all(ldim >= 1), what//': the map has positive dimensions')
        ! 0.01 A: slack on the header's sampling, not derived
        call assert_real(SMPD, find_img_smpd(vol_file), 0.01, what//': the map has the requested sampling')
    end subroutine check_map

end module simple_pdb2mrc_tester
