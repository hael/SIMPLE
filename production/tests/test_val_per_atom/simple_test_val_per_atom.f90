program simple_test_val_per_atom
! include 'simple_lib.f08' 
! use simple_atoms, only: atoms
! use simple_image
! implicit none
! #include "simple_local_flags.inc"
! character(len=STDLEN)         :: pdb_file, vol_file, pdb_out
! type(atoms)                   :: molecule 
! type(image)                   :: vol 
! real                          :: smpd
! logical                       :: pdb_exists
! integer                       :: rc, slen
! character(len=:), allocatable :: cmd
! integer                       :: ifoo, ldim(3)
! character(len=:), allocatable :: smpd_char
! if( command_argument_count() < 1 )then
!     write(logfhandle,'(a)') 'ERROR! Usage: simple_test_val_per_atom coords.pdb smpd'
!     write(logfhandle,'(a)') 'coords.pdb: PDB with the cartesian coordinates of the molecule'
!     write(logfhandle,'(a)') 'smpd               : SMPD value in Angstrom per voxel ' 
!     write(logfhandle,'(a)') 'Example: https://www.rcsb.org/structure/1jyx with smpd=0.358 mskdiam=180'
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
!     smpd          = 0.358
! else
!     call get_command_argument(1, pdb_file)
!     call get_command_argument(2, length=slen)
!     allocate(character(slen) :: smpd_char)
!     call get_command_argument(2, smpd_char)
!     read(smpd_char, *) smpd
! endif
! vol_file = swap_suffix(pdb_file,'mrc','pdb')
! pdb_out  = trim(get_fbody(pdb_file,'pdb'))//'_centered.pdb'
! call molecule%new(pdb_file)
! call molecule%pdb2mrc(pdb_file, vol_file, smpd, pdb_out=pdb_out)
! call find_ldim_nptcls(vol_file, ldim, ifoo)
! call vol%new(ldim, smpd)
! call vol%read(vol_file)
! call molecule%atom_validation(vol, 'validation_correlation')
! call molecule%kill
! call vol%kill
end program simple_test_val_per_atom 
