!@descr: library tests of pdb2mrc (simple_atoms): density maps from the built-in 6VXX and 1JYX models
! Moved from the pdb2mrc test program by the utils review (plan, section 9.7), with the checks of
! its simple_test_exec twin as assertions. Each model is converted twice, with the default file
! names (molecule.pdb, molecule.mrc) and with explicit ones; the map has to exist, have positive
! dimensions and the requested sampling. Run nightly in lib_single (the full spike and
! beta-galactosidase at 1.3 A); the files are deleted afterwards.
module simple_pdb2mrc_tester
use simple_core_module_api
use simple_atoms,         only: atoms
use simple_molecule_data, only: molecule_data, betagal_1jyx, sars_cov2_spkgp_6vxx
use simple_test_utils
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
        write(*,'(A)') 'test_pdb2mrc_model '//code
        call assert_true(mol%n > 0, code//': the built-in model has atoms')
        ! default file names
        call molecule%pdb2mrc(smpd=SMPD, mol=mol)
        call check_map(string('molecule.mrc'), code//' (default file names)')
        call molecule%kill
        ! explicit file names
        pdb_file = code//'.pdb'
        vol_file = code//'.mrc'
        call molecule%pdb2mrc(pdbfile=pdb_file, volfile=vol_file, smpd=SMPD, mol=mol)
        call check_map(vol_file, code//' (explicit file names)')
        call molecule%kill
        call del_file(string('molecule.pdb'))
        call del_file(string('molecule.mrc'))
        call del_file(pdb_file)
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
        call assert_real(SMPD, find_img_smpd(vol_file), 0.01, what//': the map has the requested sampling')
    end subroutine check_map

end module simple_pdb2mrc_tester
