program simple_test_radial_cc
! include 'simple_lib.f08'
! use simple_image, only: image
! implicit none
! #include "simple_local_flags.inc"
! real,    parameter            :: smpd=0.358 
! type(image)                   :: img1, img2, img_w
! real,    allocatable          :: rad_corrs(:), rad_dists(:), filt(:)
! character(len=256)            :: fn_img1, fn_img2
! character(len=:), allocatable :: cmd
! integer                       :: ldim_refs(3), ifoo, ldim1(3), n_shells, i, rc
! logical                       :: pdb_exists
! if( command_argument_count() /= 2 )then
!     write(logfhandle,'(a)') 'ERROR! Usage: simple_test_radial_cc img1.mrc img2.mrc'
!     write(logfhandle,'(a)') 'Example: https://www.rcsb.org/structure/1jyx with smpd=1. mskdiam=180'
!     write(logfhandle,'(a)') 'DEFAULT TEST (example above) is running now...'
!     inquire(file="1JYX.pdb", exist=pdb_exists)
!     if( .not. pdb_exists )then
!         write(*, *) 'Downloading the example dataset...'
!         cmd = 'curl -s -o 1JYX.pdb https://files.rcsb.org/download/1JYX.pdb'
!     endif
!     call execute_command_line(cmd, exitstat=rc)
!     write(*, *) 'Creating two mrcs vols...'
!     cmd     = 'e2pdb2mrc.py 1JYX.pdb vol1.mrc'
!     call execute_command_line(cmd, exitstat=rc)
!     cmd     = 'cp vol1.mrc vol2.mrc'
!     call execute_command_line(cmd, exitstat=rc)
!     cmd     = 'rm 1JYX.pdb'
!     call execute_command_line(cmd, exitstat=rc)
!     fn_img1 = 'vol1.mrc'
!     fn_img2 = 'vol2.mrc'
! else
!     call get_command_argument(1, fn_img1)
!     call get_command_argument(2, fn_img2)
! endif
! call find_ldim_nptcls(fn_img1, ldim1, ifoo)
! ldim_refs = [ldim1(1), ldim1(2), ldim1(3)]
! n_shells  = nint(ldim_refs(1) / 2.)
! call img1%new(ldim_refs, smpd)
! call img2%new(ldim_refs, smpd)
! call img_w%new(ldim_refs, smpd)
! call img1%read(fn_img1)
! call img2%read(fn_img2)
! if( .not. (img1.eqdims.img2)     ) THROW_HARD('Nonconforming dimensions of images')
! if( .not. (img1%same_smpd(img2)) ) THROW_HARD('Nonconforming smpd of images')
! allocate(rad_corrs(n_shells), rad_dists(n_shells))
! call img1%radial_cc(img2, img_w, smpd, rad_corrs, rad_dists)
! write(logfhandle,*) '             Rad. dist.         Rad. corr. '
! do i = 1, n_shells
!     write(logfhandle,*) i, rad_dists(i), rad_corrs(i)
! enddo
! call img1%kill()
! call img2%kill()
! call img_w%kill()
! cmd  = 'rm vol1.mrc vol2.mrc .eman2log.txt'
! call execute_command_line(cmd, exitstat=rc)
end program simple_test_radial_cc


