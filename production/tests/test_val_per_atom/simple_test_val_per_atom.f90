program simple_test_val_per_atom
include 'simple_lib.f08' 
use simple_atoms, only: atoms
use simple_image
implicit none
character(len=STDLEN)         :: pdb_file, vol_file
type(atoms)                   :: molecule 
type(image)                   :: vol 
real                          :: smpd=1.0
logical                       :: pdb_exists
integer                       :: rc
character(len=:), allocatable :: cmd
integer                       :: ifoo, ldim(3)
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
    !write(*, *) 'Creating vol from pdb...'
    !cmd      = 'e2pdb2mrc.py --apix 0.358 1JYX.pdb vol.mrc'
    !call execute_command_line(cmd, exitstat=rc)
    pdb_file      = '1JYX.pdb'
else
    call get_command_argument(1, pdb_file)
endif
smpd     = 1.000
vol_file = swap_suffix(pdb_file,'mrc','pdb')
call molecule%new(pdb_file)
call molecule%pdb2mrc(pdb_file, vol_file, smpd)
call find_ldim_nptcls(vol_file, ldim, ifoo)
call vol%new(ldim, smpd)
call vol%read(vol_file)
call molecule%atom_validation(vol, 'validation_correlation')
call molecule%kill
call vol%kill
end program simple_test_val_per_atom 





