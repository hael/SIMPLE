program simple_test_pdb2mrc
! include 'simple_lib.f08' 
! use simple_atoms, only: atoms
! use simple_image
! implicit none
! character(len=STDLEN)         :: pdb_file, vol_file, str_name, pdb_out
! character(len=:), allocatable :: smpd_char
! real                          :: smpd = 1.
! logical                       :: pdb_exists
! integer                       :: rc, slen
! character(len=:), allocatable :: cmd
! type(atoms)                   :: molecule 

! if( command_argument_count() /= 2 )then
!     write(logfhandle,'(a)') 'ERROR! Usage: simple_test_pdb2mrc mol.pdb            smpd'
!     write(logfhandle,'(a)') '              simple test_pdb2mrc pdb_structure_name smpd'
!     write(logfhandle,'(a)') 'mol.pdb            : PDB with the cartesian coordinates of the molecule'
!     write(logfhandle,'(a)') 'pdb_structure_name : PDB name structure from PDB database'
!     write(logfhandle,'(a)') 'smpd               : SMPD value in Angstrom per voxel ' 
!     write(logfhandle,'(a)') 'Example: https://www.rcsb.org/structure/1jyx with smpd=1. mskdiam=180'
!     write(logfhandle,'(a)') 'DEFAULT TEST (example above) is running now...'
!     pdb_file = '1JYX.pdb'
!     inquire(file=pdb_file, exist=pdb_exists)
!     if( .not. pdb_exists )then
!         write(logfhandle, *) 'Downloading the example dataset...'
!         cmd = 'curl -s -o 1JYX.pdb https://files.rcsb.org/download/1JYX.pdb'
!         call execute_command_line(cmd, exitstat=rc)
!     endif
!     smpd = 1.
! else
!     call get_command_argument(1, pdb_file)
!     call get_command_argument(2, length=slen)
!     allocate(character(slen) :: smpd_char)
!     call get_command_argument(2, smpd_char)
!     read(smpd_char, *) smpd
!     if( fname2ext(pdb_file) .ne. 'pdb' )then
!         str_name = pdb_file 
!         pdb_file = trim(str_name)//'.pdb'
!         inquire(file=pdb_file, exist=pdb_exists)
!         if( .not. pdb_exists )then
!             write(logfhandle,'(A,A)') "Experimental PDB Structure: ",trim(str_name)
!             write(logfhandle, *) 'Downloading '//trim(str_name)//' from PDB database:'
!             cmd = 'curl -s -o '//trim(pdb_file)//' https://files.rcsb.org/download/'//trim(pdb_file)
!             write(logfhandle, *) ' https://files.rcsb.org/download/'//trim(pdb_file)
!             call execute_command_line(cmd, exitstat=rc)
!         else
!             write(logfhandle,'(A,A)') "Experimental PDB Structure: ",trim(str_name)
!         endif
!     else 
!         inquire(file=pdb_file, exist=pdb_exists)
!         if( .not. pdb_exists )then
!             str_name = get_fbody(pdb_file,'.pdb')
!             write(logfhandle,'(A,A)') "Experimental PDB Structure: ",trim(str_name)
!             write(logfhandle, *) 'Downloading '//trim(str_name)//' from PDB database:'
!             cmd = 'curl -s -o'//trim(pdb_file)//' https://files.rcsb.org/download/'//trim(pdb_file)
!             call execute_command_line(cmd, exitstat=rc)
!         else
!             write(logfhandle,'(A,A)') "Experimental PDB Structure: ",trim(str_name)
!         endif
!     endif 
! endif
! vol_file = swap_suffix(pdb_file,'mrc','pdb') 
! pdb_out  = trim(get_fbody(pdb_file,'pdb'))//'_centered'//'.pdb'
! call molecule%new(pdb_file)
! call molecule%pdb2mrc(pdb_file, vol_file, smpd, pdb_out=pdb_out)
! call molecule%kill
end program simple_test_pdb2mrc
