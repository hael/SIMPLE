program simple_test_pdb2mrc
use simple_core_module_api
use simple_image, only : image
use simple_atoms, only : atoms
use simple_molecule_data, only : molecule_data, betagal_1jyx, sars_cov2_spkgp_6vxx
implicit none
type(string)        :: pdb_file, vol_file
type(atoms)         :: molecule
type(molecule_data) :: mol
real, parameter     :: smpd = 1.3
mol = sars_cov2_spkgp_6vxx()
pdb_file = '6VXX.pdb'
vol_file = '6VXX.mrc'
call molecule%pdb2mrc( pdb_file, vol_file, smpd, mol=mol )
mol = betagal_1jyx()
pdb_file = '1JYX.pdb'
vol_file = '1JYX.mrc'
call molecule%pdb2mrc( pdb_file, vol_file, smpd, mol=mol )
end program simple_test_pdb2mrc
