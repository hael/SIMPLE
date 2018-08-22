module simple_test_chiara_try_mod
include 'simple_lib.f08'
use simple_image, only : image
implicit none
#include "simple_local_flags.inc"

contains

  function center_template(self, part_radius) result(xy)
      class(image), intent(inout) :: self
      real,         intent(in)    :: part_radius
      integer      :: xy(2)         !coordinates of the new center
      type(image)  :: template, self_pad, self_compare
      integer      :: i, j          !for loops
      integer      :: ldim(3)       !size of the window self
      integer      :: cc_cont, cnt  !counters
      real,    allocatable :: mat_compare(:,:,:), rmat(:,:,:)
      logical, allocatable :: l_mask(:,:)
      real,    allocatable :: correlation_coefficients(:,:)
      integer, allocatable :: max_loc(:)
      logical, parameter   :: DEBUG = .true.
      !ldim = self%ldim
      ldim = self%get_ldim()
      if(ldim(3) /=1 ) stop 'simple_image:: center_template. This function is for 2D images!'
      allocate(mat_compare( ldim(1),ldim(2),1), source = 0. )
      call self_compare%new([ldim(1),ldim(2),1], 1.)  !to compare self and template
      !Creating template
      call template%disc([ldim(1),ldim(2),1], 1., part_radius)
      !Preparation and calculation of cc
      call self_pad%new([2*ldim(1),2*ldim(2),1],1.) !padding preparation
      call self%pad(self_pad)
      rmat = self_pad%get_rmat()
      !cc_cont = count(self_pad%rmat>0)
      cc_cont = count(rmat>0)
      if(DEBUG) print *, 'You will calculate # ', cc_cont, ' correlations'
      allocate(correlation_coefficients(3,cc_cont), source = 0.) !correlation and pixel position
      cnt = 0
      do i = 1, ldim(1)
        do j = 1, ldim(2)
          !if(self_pad%rmat(i,j,1) > 0.) then   !I am not sure about this, it's ok but first I have to close....
          if(rmat(i,j,1) > 0.) then
            cnt = cnt + 1
            mat_compare(:,:,:) = rmat( i-ldim(1)/2:i+ldim(1)/2-1, j-ldim(2)/2:j+ldim(2)/2-1,:)
            call self_compare%set_rmat(mat_compare)
            correlation_coefficients(:2,cnt) = [i,j]
            correlation_coefficients(3, cnt) = template%real_corr(self_compare)
          endif
        enddo
      enddo
      !Creating mask for maximum calculation along the first column
      allocate(l_mask(3, cc_cont), source = .false.)
      l_mask(3,:) = .true.
      max_loc = maxloc(correlation_coefficients, l_mask)
      if(DEBUG) print *, "Location (pixel) :", correlation_coefficients(:2,max_loc(2))
      xy = correlation_coefficients(:2,max_loc(2)) - ldim(:2)
  end function center_template
end module simple_test_chiara_try_mod

program simple_test_chiara_try
include 'simple_lib.f08'
use simple_micops
use simple_image
use simple_stackops
use simple_math
use simple_edge_detector
use simple_test_chiara_try_mod
type(image) :: img_in
real        :: part_radius



end program simple_test_chiara_try
! matrix = reshape(real([ 1,1,1,0,0,6,5, &
!                  & 1,1,0,0,6,6,6, &
!                  & 1,0,0,2,0,6,0, &
!                  & 0,0,2,2,0,0,4, &
!                  & 0,5,0,0,0,4,4, &
!                  & 0,5,5,5,0,0,0, &
!                  & 0,5,5,0,0,3,3]),[7,7,1])
! ldim = [box,box,1]
! call img_in%new(ldim,1.)
! call img_in%disc( ldim, 1.0 , 10.)
! call build_ellipse(img_in,[10.,10.], [4.,4.], 0.)
! call build_ellipse(img_in,[20.,20.], [4.,4.], 0.)
