module simple_test_chiara_binarize_cavgs_mod

    include 'simple_lib.f08'
    use simple_image,         only : image

contains

! In this subroutine binarize stack of class averages
! Image processing steps: 1) Low pass filtering
!                         2) Binarization (Otsu)
subroutine binarize_cavg(stk_in, stk_out, smpd)
        use simple_segmentation, only : otsu_img
        character(len=*), intent(in) :: stk_in
        character(len=*), intent(in) :: stk_out
        real,             intent(in) :: smpd
        type(image)       :: img
        integer           :: nptcls, i, ldim(3)
        call find_ldim_nptcls(stk_in, ldim, nptcls)
        ldim(3)=1
        call img%new(ldim, smpd)
        do i = 1,nptcls
        call img%read(stk_in, i)
        ! 1) Low filtering
        call img%fft()
        call img%bp(0.,5.) !lp = 5
        call img%ifft()
        ! 2) Binarization
        call img%mask(30.,'hard')
        call otsu_img(img)
        call img%write(stk_out,i)
        enddo
end subroutine binarize_cavg
end module simple_test_chiara_binarize_cavgs_mod

program simple_test_chiara_binarize_cavgs
  include 'simple_lib.f08'
  use simple_image,         only : image
  use simple_parameters,    only: parameters
  use simple_cmdline,       only: cmdline
  use simple_test_chiara_binarize_cavgs_mod
  type(cmdline)      :: cline
  type(parameters)   :: params
  type(image)        :: img, img_cc, imgwin, img_pad
  real               :: smpd
  real               :: lp
  integer :: nptcls, i, ldim(3), px(2)
  real, allocatable :: rmat(:,:,:)
  logical :: outside
  if( command_argument_count() < 2 )then
      write(logfhandle,'(a)',advance='no') 'simple_test_chiara_try smpd=<sampling distance(in A)> [stk = stack file name]'
      stop
  endif
  call cline%parse_oldschool
  call cline%checkvar('stk', 1)
  call cline%checkvar('smpd',2)
  !Set defaults
  call params%new(cline)  !<read cline parameters
  call binarize_cavg(params%stk, 'OtsuBin.mrc', params%smpd)
end program simple_test_chiara_binarize_cavgs
