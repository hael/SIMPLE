program simple_test_atomfit
use simple_core_module_api
use simple_image, only : image
use simple_atoms, only : atoms
implicit none
type(atoms)         :: molecule
real, parameter     :: smpd = 0.4
call molecule%pdb2mrc(pdbfile=string('CdSe_wurtzite_2nm_centered.pdb'),smpd=smpd)
call molecule%fit_5gauss( volfile='molecule.mrc', smpd_target=0.05 )
end program simple_test_atomfit
