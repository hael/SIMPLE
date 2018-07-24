module simple_test_chiara_NLmean_mod
  include 'simple_lib.f08'
  use simple_image, only : image
  implicit none
  integer, private :: ldim(3,2)
  integer, public, parameter  :: dim_sw = 7
  contains

  subroutine set_ldim_pad(ldim_in, pad)
    integer, intent(in) :: ldim_in(3), pad
    ldim(1,1) = 1-pad
    ldim(1,2) = ldim_in(1)+pad
    ldim(2,1) = 1-pad
    ldim(2,2) = ldim_in(2)+pad
    ldim(3,:) = 1
  end subroutine set_ldim_pad

  subroutine print_mat1(matrix)
    real, intent(in)    :: matrix(:,:,:)
    integer             :: j, s(3)
    s = shape(matrix)
    if(s(3) /= 1) then
      print *, 'You can only print 2D matrices!'
      stop
    endif
    do j = 1, s(2)
      print *, matrix(:,j,1)
    enddo
  end subroutine print_mat1

  function similarity_window(rmat,px) result(sw)
    real,        intent(in) :: rmat(ldim(1,1):ldim(1,2),ldim(2,1):ldim(2,2),ldim(3,1):ldim(3,2))
    integer,     intent(in) :: px(3)
    real                    :: sw_mat(dim_sw,dim_sw,1), sw(dim_sw*dim_sw)

    sw_mat(:dim_sw,:dim_sw,1) = rmat(px(1):px(1)+dim_sw-1, px(2):px(2)+dim_sw-1,1)
    sw = reshape(sw_mat,[dim_sw*dim_sw])
  end function similarity_window

  subroutine build_ellipse(img_in,center, axis, rot)
      type(image), intent(inout) :: img_in
      real,        intent(in)    :: center(2), axis(2), rot
      real, allocatable          :: theta(:)
      integer                    :: ldim(3), i, j, k
      real :: rot_r
      if(rot < 0. .or. rot > 360. ) then
        print *, "Please insert an angle in the range [0,360]"
        stop
      endif
      rot_r = deg2rad(rot)
      ldim = img_in%get_ldim()
      if(ldim(3) /= 1) then
        print *, "The image has to be 2D!"
        stop
      endif
      theta = (/ (deg2rad(real(i)),i=1,360) /)
      do k = 1,size(theta)
          do i = 1,ldim(1)
              do j = 1,ldim(2)
                  if(abs(real(i)-center(1)-axis(1)*cos(theta(k))*cos(rot_r)+axis(2)*sin(theta(k))*sin(rot_r))<1 .and. &
                   & abs(real(j)-center(1)-axis(1)*cos(theta(k))*sin(rot_r)-axis(2)*sin(theta(k))*cos(rot_r))<1) then
                  call img_in%set([i,j,1], 1.)
                  end if
              enddo
          enddo
      enddo
      deallocate(theta)
  end subroutine build_ellipse
end module simple_test_chiara_NLmean_mod

program simple_test_chiara_NLmean
use simple_test_chiara_NLmean_mod
implicit none
real, allocatable  ::  NL(:,:,:), Z(:,:,:), rmat_pad(:,:,:), rmat(:,:,:)
real               :: sw(dim_sw*dim_sw), sw_px(dim_sw*dim_sw)
integer            :: i, j, m, n, ldim(3), pad
integer, parameter :: cfr_box = 10, box = 300
real, parameter    :: sigma = 5., h =4.*sigma !SIGMA, H..?
type(image)        :: img_in, img, img_pad
real               :: w
!pad = 2*cfr_box
pad = (dim_sw-1)/2
ldim = [box, box, 1]
call set_ldim_pad(ldim, pad)
call img_pad%new([ldim(1)+2*pad,ldim(2)+2*pad,1], 1.)  !you can also insert a control so that you pad only if necessary
allocate(rmat_pad(1-pad:ldim(1)+pad, 1-pad:ldim(2)+pad,1), source=0.)
!call img_in%new([3710,3838,1],1.)
!call img_in%read('/home/lenovoc30/Desktop/ANTERGOS/forctf/0001_forctf.mrc')
!call img_in%vis
!ldim = img_in%get_ldim()
call img_in%disc( ldim, 1.0 , 10.)
call build_ellipse(img_in,[10.,10.], [4.,4.], 0.)
call build_ellipse(img_in,[20.,20.], [4.,4.], 0.)
call build_ellipse(img_in,[90.,10.], [4.,4.], 0.)
call build_ellipse(img_in,[30.,30.], [4.,4.], 0.)
call build_ellipse(img_in,[40.,40.], [4.,4.], 0.)
call build_ellipse(img_in,[50.,50.], [4.,4.], 0.)
call build_ellipse(img_in,[60.,60.], [4.,4.], 0.)
call build_ellipse(img_in,[70.,70.], [4.,4.], 0.)
call build_ellipse(img_in,[80.,80.], [4.,4.], 0.)
call build_ellipse(img_in,[90.,90.], [4.,4.], 0.)
call build_ellipse(img_in,[100.,100.], [4.,4.], 0.)
call build_ellipse(img_in,[110.,110.], [4.,4.], 0.)
call build_ellipse(img_in,[120.,120.], [4.,4.], 0.)
call img_in%vis
allocate(NL(box,box,1),Z(box,box,1), source=0.)  !substitute box with ldim
call img_in%add_gauran(10.)
call img_in%vis
rmat = img_in%get_rmat()
rmat_pad(1:ldim(1),1:ldim(2),:) = rmat
!deallocate(rmat)
do m = 1,ldim(1) !I fix pixel (m,n)
  do n = 1,ldim(2)
    !sw_px = similarity_window(img_in,[m,n,1])
    do i = -cfr_box,cfr_box       !As suggested in the paper, I just consider a box 21x21
        do j = -cfr_box,cfr_box
          if(i/=0 .or. j/=0) then !Not to take pixel (m,n)
              if(m+i>0 .and. m+i<=ldim(1) .and. n+j>0 .and. n+j<=ldim(2)) then
                  !sw = similarity_window(img_in,[m+i,n+j,1])
                  Z(m,n,1) = Z(m,n,1)+exp(-((norm2(similarity_window(rmat_pad,[m,n,1])-similarity_window(rmat_pad,[m+i,n+j,1])))**2./h**2.))
              endif
          endif
        enddo
    enddo
  enddo
enddo

do m = 1,ldim(1) !I fix pixel (m,n)
  do n = 1,ldim(2)
    sw_px = similarity_window(rmat_pad,[m,n,1])
    do i = -cfr_box,cfr_box   !As suggested in the paper, I just consider a box 21x21
        do j = -cfr_box,cfr_box
          if(i/=0 .or. j/=0) then !Not to take pixel (m,n)
              if(m+i>0 .and. m+i<=ldim(1) .and. n+j>0 .and. n+j<=ldim(2)) then
                  sw = similarity_window(rmat_pad,[m+i,n+j,1])
                  w = (exp(-((norm2(sw_px-sw))**2./h**2.)))/Z(m,n,1)
                  NL(m,n,1) = NL(m,n,1)+w*img_in%get([m+i,n+j,1])
              endif
          endif
        enddo
    enddo
  enddo
enddo
call img_in%set_rmat(NL)
call img_in%vis
deallocate(rmat_pad,Z,NL)
!print *, "Finished Work"
end program simple_test_chiara_NLmean
