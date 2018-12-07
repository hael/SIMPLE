module simple_test_chiara_binarize_cavgs_mod

    include 'simple_lib.f08'
    use simple_image,         only : image

contains

! In this subroutine binarize stack of class averages
! Image processing steps: 1) Low pass filtering
!                         2) Binarization (Otsu)
!                         3) Median Filtering (don't know if necessary)
subroutine binarize_cavg(stk_in, stk_out, smpd)
        use simple_segmentation, only : otsu
        character(len=*), intent(in) :: stk_in
        character(len=*), intent(in) :: stk_out
        real,             intent(in) :: smpd
        type(image) :: img, img_bin
        integer :: nptcls, i, ldim(3)
        call find_ldim_nptcls(stk_in, ldim, nptcls)
        ldim(3)=1
        call img%new(ldim, smpd)
        call img_bin%new(ldim, smpd)
        do i = 1,nptcls
        call img%read(stk_in, i)
        ! 1) Low filtering
        call img%fft()
        call img%bp(0.,5.) !lp = 5
        call img%ifft()
        ! 2) Binarization
        call img%mask(30.,'hard')
        call otsu(img, img_bin)
        call img_bin%write('OtsuBin.mrc',i)
        call img_bin%real_space_filter(3,'median')
        call img_bin%write(stk_out,i)
        enddo
end subroutine binarize_cavg
end module simple_test_chiara_binarize_cavgs_mod

program simple_test_chiara_binarize_cavgs
  include 'simple_lib.f08'
  use simple_image,         only : image
  use simple_parameters,    only: parameters
  use simple_cmdline,       only: cmdline
  use simple_picker_chiara
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
  call binarize_cavg(params%stk, 'OtsuBinMedian.mrc', params%smpd)

  ! !Centering again
  ! call img%new(ldim, params%smpd)
  ! call img_cc%new(ldim, params%smpd)
  ! do i = 1,nptcls
  !         outside = .false.
  !         call img%read('OtsuBinMedian.mrc', i)
  !         ! CC selection (cuz of the masking there might be a white frame around the particle)
  !         call img%find_connected_comps(img_cc)
  !         call img_cc%write('CC.mrc',i)
  !         rmat = img_cc%get_rmat()
  !         if(maxval(rmat) > 1.) then
  !             where(rmat == 1.) rmat = 0. !delete the frame
  !         endif

          ! How to procede: are there still more than one cc?
          ! morpho_open the image until you have just one cc
          ! than calculate the central pixel of the cc (as usual)
          ! and center back the cavg

          !Are there still more than one cc??
          ! do while(maxval(rmat) > 1.)
          !     call img_out%morpho_opening()
          !     rmat = img%get_rmat()
          ! enddo
          ! rmat_masked = int(rmat)
          ! px = center_cc(rmat_masked)
          ! call img%pad(img_pad)
          ! call img%write('img_copy.mrc', i)
          ! call img_pad%window_slim(px-ldim(1)/8, ldim(1), imgwin, outside)
          ! if (.not. outside) call imgwin%write('Centered.mrc', i)
! enddo
end program simple_test_chiara_binarize_cavgs
