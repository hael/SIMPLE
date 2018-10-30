module simple_test_chiara_thonrings
  include 'simple_lib.f08'
  use simple_image, only : image
  implicit none
  contains
      subroutine build_resolutions_vector(smpd, vector)
          real,    intent(in) :: smpd
          real, intent(out) :: vector(5,2)
          vector(:,1) = 0.
          vector(:,2) = [0., 1., 2., 3., 4.]
          !shells:   -) 0 - 30
          !          -) 30 - 20
          !          -) 20 - 15
          !          -) 15 - 10
          !          -) 10 - 8
          !          -) 8 - 5
          !          -) 5 - 4
          !          -) 4 - 3
          !          -) 3 - NyFr
      end subroutine build_resolutions_vector
end module simple_test_chiara_thonrings

program simple_test_chiara_try
  include 'simple_lib.f08'
  use simple_image
  use simple_math
  use simple_test_chiara_thonrings
  integer :: i, j, ldim(3)
  real :: matrix(7,7,1)
  integer            :: sh, h, k, ind(2)
  real :: vector(5,2), smpd
  integer ::  counter(5)
  logical :: mask(5,2)
  type(image) :: img
  real, allocatable :: rmat(:,:,:)
matrix = reshape(real([ 1,1,1,0,0,6,5, &
                 & 1,1,0,0,6,6,6, &
                 & 1,0,0,2,0,6,0, &
                 & 0,0,2,2,0,0,4, &
                 & 0,5,0,0,0,4,4, &
                 & 0,5,5,5,0,0,0, &
                 & 0,5,5,0,0,3,3]),[7,7,1])
call img%new([7,7,1],1.)
call img%set_rmat(matrix)
ldim = [7,7,1]
smpd = 1.
rmat = img%get_rmat()
call build_resolutions_vector(smpd, vector)
counter = 0
print *, 'matrix = '
call vis_mat(matrix)
mask(:,1) = .false.
mask(:,2) = .true.
do i = 1, ldim(1)
    do j = 1, ldim(2)
            h = -int(7/2) + i - 1
            k = -int(7/2) + j - 1
            sh = nint(hyp(real(h),real(k)))
            ind = minloc(abs(vector-sh),mask)
            counter(ind(1)) = counter(ind(1)) + 1   !!Number of pixels per shell, it is fixed, I could calculate apart
            if (rmat(i,j,1) > 0.5) then !binary image
                print *, 'pixel ', i,j, 'h k = ', h,k, 'sh = ', sh,  'ind', ind(1)
                vector(ind(1),1) = vector(ind(1),1) + 1. !update # of white pixel in the shell
        endif
    enddo
enddo
print *, 'counter ', counter
print *, 'vector '
call vis_mat(vector)
! normalise
where(counter > 0) vector(:,1) = vector(:,1) / counter(:)
print *, 'vector normalized'
call vis_mat(vector)
end program simple_test_chiara_try
