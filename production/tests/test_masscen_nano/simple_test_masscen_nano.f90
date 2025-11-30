program simple_test_masscen_nano
! include 'simple_lib.f08'
! use simple_atoms
! implicit none
! character(len=*), parameter :: atmsin = '1.pdb'
! character(len=*), parameter :: atmsout = '2'
! real        :: m(3)
! type(atoms) :: atom_centers
! call atom_centers%new(atmsin)
! m = atom_centers%find_masscen()
! print *, 'center of mass: ', m
! !stop
! call atom_centers%writePDB(string(atmsout)//'.pdb')
! call atom_centers%kill
end program simple_test_masscen_nano
