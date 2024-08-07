program simple_test_val_per_atom
include 'simple_lib.f08' 
use simple_atoms, only: atoms
use simple_image
implicit none
!type(image)                   :: sim_vol, vol
character(len=STDLEN)         :: pdb_file, vol_file, sim_vol_file
!real, allocatable             :: xyz(:,:)
!real, parameter               :: smpd=0.358
real, parameter               :: smpd=0.6
logical                       :: pdb_exists
integer                       :: rc
character(len=:), allocatable :: cmd
type(atoms)                   :: molecule 
integer                       :: ifoo, ldim(3)
!integer                       :: i, natoms, ifoo, ldim(3)
!real                          :: shift(3)

if( command_argument_count() /= 1 )then
    write(logfhandle,'(a)') 'ERROR! Usage: simple_test_val_per_atom coords.pdb vol.mrc'
    write(logfhandle,'(a)') 'coords.pdb: PDB with the cartesian coordinates of the molecule'
    write(logfhandle,'(a)') 'vol.mrc: 3D recontrusted volume'
    write(logfhandle,'(a)') 'Example: https://www.rcsb.org/structure/1jyx with smpd=1. mskdiam=180'
    write(logfhandle,'(a)') 'DEFAULT TEST (example above) is running now...'
    inquire(file="1JYX.pdb", exist=pdb_exists)
    if( .not. pdb_exists )then
        write(*, *) 'Downloading the example dataset...'
        cmd = 'curl -s -o 1JYX.pdb https://files.rcsb.org/download/1JYX.pdb'
        call execute_command_line(cmd, exitstat=rc)
    endif
    write(*, *) 'Creating vol from pdb...'
    !cmd      = 'e2pdb2mrc.py 1JYX.pdb vol.mrc'
    cmd      = 'e2pdb2mrc.py --apix 0.358 1JYX.pdb vol.mrc'
    !cmd      = 'e2pdb2mrc.py 1JYX.pdb vol.mrc'
    call execute_command_line(cmd, exitstat=rc)
    pdb_file      = '1JYX.pdb'
else
    call get_command_argument(1, pdb_file)
endif
vol_file     = 'vol.mrc'
sim_vol_file = 'sim_vol.mrc'
! extract coordinates and elements from PDB file
call molecule%new(pdb_file)
call molecule%pdb2mrc(pdb_file, sim_vol_file, smpd)
call molecule%kill
call find_ldim_nptcls(vol_file, ldim, ifoo ) ! ifoo is not used but must be included in the function call
print *, "ldim", ldim
!stop
!!call molecule%translate_posqd
!natoms = molecule%get_n()
!print *, natoms
!allocate(xyz(natoms,3), source=0.)
!do i = 1, natoms
!    print *, molecule%get_element(i), molecule%get_coord(i)
!enddo
!! Load the experimental image
!call find_ldim_nptcls(vol_file, ldim, ifoo ) ! ifoo is not used but must be included in the function call
!! 0,0,0 in PDB space is map to the center of the volume
!shift = ( ldim/2. + 1 ) * smpd
!call molecule%translate(shift)
!call vol%new(ldim, smpd)
!call sim_vol%new(ldim, smpd)
!call molecule%convolve(sim_vol, cutoff = 8*smpd ) 
!    call sim_vol%write(sim_vol_file)
!
end program simple_test_val_per_atom 





