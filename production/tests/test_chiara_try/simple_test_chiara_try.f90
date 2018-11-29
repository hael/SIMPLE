module simple_test_chiara_try_mod
  include 'simple_lib.f08'
  use simple_image, only : image
  implicit none
  public

  contains

end module simple_test_chiara_try_mod

program simple_test_chiara_try
  include 'simple_lib.f08'
  use simple_powerspec_analysis
  use gnufor2
  use simple_picker_chiara
  use simple_micops
  use simple_image
  use simple_stackops
  use simple_math
  use simple_test_chiara_try_mod
  use simple_segmentation
  use simple_parameters, only: parameters
  use simple_cmdline,    only: cmdline
  type(image)       :: img, img_out
  real, allocatable :: rmat(:,:,:), x(:)
  integer :: i, ldim(3), nptcls, npxls_at_mode
  real    ::  m(1), smpd
  real :: vec(10)
  real,    allocatable :: xhist(:)
  integer, allocatable :: yhist(:)

  smpd = 1.32
  call find_ldim_nptcls('/home/lenovoc30/Desktop/PickingResults/SomeExamples/0001_forctf.mrc', ldim, nptcls)
  ! call find_ldim_nptcls('/home/lenovoc30/Desktop/PickingResults/SomeExamples/NegStainingWorking/shrunken_hpassfiltered.mrc', ldim, nptcls)
  call img%new(ldim, smpd)
  call img%read('/home/lenovoc30/Desktop/PickingResults/SomeExamples/0001_forctf.mrc')
  ! call img%read('/home/lenovoc30/Desktop/PickingResults/SomeExamples/NegStainingWorking/shrunken_hpassfiltered.mrc')
  call img%hist_stretching(img_out)
  call img_out%write('HistogramStretched.mrc')
end program simple_test_chiara_try
! matrix = reshape(real([ 1,1,1,0,0,6,5, &
!                  & 1,1,0,0,6,6,6, &
!                  & 1,0,0,2,0,6,0, &
!                  & 0,0,2,2,0,0,4, &
!                  & 0,5,0,0,0,4,4, &
!                  & 0,5,5,5,0,0,0, &
!                  & 0,5,5,0,0,3,3]),[7,7,1])
