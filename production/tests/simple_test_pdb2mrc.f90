program simple_test_pdb2mrc
use simple_core_module_api
use simple_image, only : image
use simple_atoms, only : atoms
use simple_molecule_data, only : molecule_data, betagal_1jyx_molecule
implicit none
type(string)        :: pdb_file, vol_file
type(atoms)         :: betagal
type(molecule_data) :: mol
real, parameter     :: smpd = 1.
mol = betagal_1jyx_molecule()
pdb_file = '1JYX.pdb'
vol_file = '1JYX.mrc'
call betagal%pdb2mrc( pdb_file, vol_file, smpd, mol=mol )
end program simple_test_pdb2mrc