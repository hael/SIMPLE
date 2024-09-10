program simple_test_val_map
include 'simple_lib.f08' 
use simple_atoms, only: atoms
use simple_image
implicit none
character(len=STDLEN)         :: pdb_file, pdb_out, exp_vol_file, sim_vol_file, even_vol_file, odd_vol_file
type(atoms)                   :: molecule 
type(image)                   :: exp_vol, sim_vol, even_vol, odd_vol 
real                          :: smpd=1.0
logical                       :: pdb_exists
integer                       :: rc
character(len=:), allocatable :: cmd
integer                       :: ifoo, ldim(3)
! if( command_argument_count() /= 1 )then
!     write(logfhandle,'(a)') 'ERROR! Usage: simple_test_val_per_atom coords.pdb vol.mrc'
!     write(logfhandle,'(a)') 'coords.pdb: PDB with the cartesian coordinates of the molecule'
!     write(logfhandle,'(a)') 'vol.mrc: 3D recontrusted volume'
!     write(logfhandle,'(a)') 'Example: https://www.rcsb.org/structure/1jyx with smpd=1. mskdiam=180'
!     write(logfhandle,'(a)') 'DEFAULT TEST (example above) is running now...'
!     inquire(file="1JYX.pdb", exist=pdb_exists)
!     if( .not. pdb_exists )then
!         write(*, *) 'Downloading the example dataset...'
!         cmd = 'curl -s -o 1JYX.pdb https://files.rcsb.org/download/1JYX.pdb'
!         call execute_command_line(cmd, exitstat=rc)
!     endif
!     !write(*, *) 'Creating vol from pdb...'
!     !cmd      = 'e2pdb2mrc.py --apix 0.358 1JYX.pdb vol.mrc'
!     !call execute_command_line(cmd, exitstat=rc)
!     pdb_file      = '1JYX.pdb'
! else
!     call get_command_argument(1, pdb_file)
! endif
smpd          = 1.000
pdb_file      = "1JYX.pdb"
exp_vol_file  = swap_suffix(pdb_file,'mrc','pdb')
sim_vol_file  = trim(get_fbody(pdb_file,'pdb'))//'_sim.mrc'
even_vol_file = trim(get_fbody(pdb_file,'pdb'))//'_even.mrc'
odd_vol_file  = trim(get_fbody(pdb_file,'pdb'))//'_odd.mrc'
pdb_out       = trim(get_fbody(pdb_file,'pdb'))//'_centered'

! the dimensions of all volumes need to be consisted
call find_ldim_nptcls(exp_vol_file, ldim, ifoo)
call molecule%new(pdb_file)
call molecule%pdb2mrc(pdb_file, sim_vol_file, smpd, pdb_out, vol_dim=ldim)
call sim_vol%new(ldim, smpd)
call sim_vol%read(sim_vol_file)
call molecule%atom_validation(sim_vol, 'sim_val_corr')
call sim_vol%kill()
call molecule%kill()

! experimental volume - even and odd

call molecule%new(pdb_file)
call molecule%center_pdbcoord(ldim, smpd)
call exp_vol%new(ldim, smpd)
call exp_vol%read(exp_vol_file)
call molecule%atom_validation(exp_vol, 'exp_val_corr')
call exp_vol%kill()
! even
call exp_vol%new(ldim, smpd)
call exp_vol%read(even_vol_file)
call molecule%atom_validation(exp_vol, 'even_val_corr')
call exp_vol%kill()
! odd
call exp_vol%new(ldim, smpd)
call exp_vol%read(odd_vol_file)
call molecule%atom_validation(exp_vol, 'odd_val_corr')
call exp_vol%kill()
call molecule%kill()

end program simple_test_val_map





