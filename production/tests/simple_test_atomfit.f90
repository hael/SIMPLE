program simple_test_atomfit
use simple_core_module_api
use simple_image, only : image
use simple_atoms, only : atoms
implicit none
type(atoms)         :: molecule
real, parameter     :: smpd = 0.358
call molecule%pdb2mrc(pdbfile=string('CdSe_wurtzite_2nm_centered.pdb'),smpd=smpd)
call molecule%fit_bfactors( volfile='CdSe_wurtzite_2nm_centered.mrc', smpd_target=0.05 )
! call molecule%pdb2mrc(pdbfile=string('simatms.pdb'),smpd=smpd)
! call molecule%fit_bfactors( volfile='outvol.mrc', smpd_target=0.05 )
end program simple_test_atomfit
