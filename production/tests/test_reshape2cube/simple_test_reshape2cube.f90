program simple_test_reshape2cube
! include 'simple_lib.f08'
! use simple_image, only: image
! implicit none
! character(len=:), allocatable :: smpd_char
! integer               :: ifoo, ldim(3), slen
! real                  :: smpd
! type(image)           :: vol_in, vol_out
! character(len=STDLEN) :: fname
! call get_command_argument(1, fname)
! call get_command_argument(2, length=slen)
! allocate(character(slen) :: smpd_char)
! call get_command_argument(2, smpd_char)
! read(smpd_char, *) smpd
! call find_ldim_nptcls(fname, ldim, ifoo, smpd = smpd) ! ifoo is not used but must be included in the function call
! call vol_in%new(ldim, smpd)
! call vol_in%read(fname)
! call vol_in%reshape2cube(vol_out)
! call vol_out%write('cube.mrc') 
end program simple_test_reshape2cube
