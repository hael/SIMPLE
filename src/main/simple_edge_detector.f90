module simple_edge_detector
include 'simple_lib.f08'
use simple_image, only: image
implicit none

public :: sobel, automatic_thresh_sobel, canny

private
#include "simple_local_flags.inc"
contains

  ! Classic Sobel edge detection routine.
  subroutine sobel(img_in,img_out,thresh)
      class(image), intent(inout) :: img_in,img_out   !image input and output
      real,         intent(in)    :: thresh(1)        !threshold for Sobel algorithm
      real,  allocatable :: grad(:,:,:), rmat(:,:,:)  !matrices
      integer :: ldim(3)
      ldim = img_in%get_ldim()
      call img_out%new(ldim,1.)          !reset if not empty
      allocate(rmat(ldim(1),ldim(2),1 ), source = 0.)
      call img_in%calc_gradient(grad)
      where(grad > thresh(1) ) rmat = 1.
      call img_out%set_rmat(rmat)
      deallocate(grad,rmat)
  end subroutine sobel

  ! This soubroutine performs sobel edge detection with an automatic iterative
  ! threshold selection in order to obtian a threshold which leads to
  ! a ratio between foreground and background pixels in the binary image
  ! close to 'goal'. This particular number is selected because it is what needed
  ! for Chiara's particle picking routine.
  subroutine automatic_thresh_sobel(img_in, bin_img, goal, thresh)
      class(image),    intent(inout)   :: img_in
      class(image),    intent(out)     :: bin_img
      real, optional,  intent(in)      :: goal
      real, optional,  intent(out)     :: thresh(1) !to keep track of the selected threshold
      real, allocatable  :: grad(:,:,:)
      real               :: tthresh(1), ggoal       !selected threshold and ratio nforeground/nbackground
      real               :: t(3)                    !to calculate 3 possible sobel thresholds for each iteration
      real               :: ratio(3)                !n_foreground/nbackground for 3 different values of the thresh
      real               :: diff(3)                 !difference between the ideal ratio (= goal) and the obtained 3 ratios
      integer, parameter :: N_MAXIT = 100           !max number of iterations
      integer            :: n_it
      logical            :: DEBUG
      real               :: r                       !random # which guarantees me not to get stuck
      DEBUG  = .true.
      ggoal = 2.
      if(present(goal)) ggoal = goal
      if(bin_img%exists()) call bin_img%kill
      call bin_img%new(img_in%get_ldim(), img_in%get_smpd())
      call img_in%calc_gradient(grad)
      tthresh(1) = (maxval(grad) + minval(grad))/2   !initial guessing
      diff = 1.     !initialization
      if(DEBUG) print*, 'Initial guessing: ', tthresh
      n_it = 0
      do while(abs(minval(diff(:))) > 0.005 .and. n_it < N_MAXIT)
          n_it = n_it + 1
          call random_seed()
          call random_number(r)                   !rand # in  [ 0,1 ), uniform distribution
          r = r*(maxval(grad) + minval(grad))/100 !rand # in  [ 0,(max(grad)-min(grad))/100 )
          t(1) = tthresh(1) + r
          t(2) = 3./2.*tthresh(1) + r
          t(3) = 3./4.*tthresh(1) + r
          if(t(2) > maxval(grad) .or. t(3) < minval(grad)) THROW_HARD('cannot find the threshold! automatic_thresh_sobel')
          call sobel(img_in, bin_img, t(1))
          if(abs(real(bin_img%nbackground())) > TINY) then !do not divide by 0
              ratio(1) = real(bin_img%nforeground())/real(bin_img%nbackground())  !obtained ratio with threshold = t(1)
          else
              ratio(1) = 0.
          endif
          call sobel(img_in, bin_img, t(2))
          if(abs(real(bin_img%nbackground())) > TINY) then !do not divide by 0
              ratio(2) = real(bin_img%nforeground())/real(bin_img%nbackground())  !obtained ratio with threshold = t(2)
          else
              ratio(2) = 0.
          endif
          call sobel(img_in, bin_img, t(3))
          if(abs(real(bin_img%nbackground())) > TINY) then !do not divide by 0
              ratio(3) = real(bin_img%nforeground())/real(bin_img%nbackground())  !obtained ratio with threshold = t(3)
          else
              ratio(3) = 0.
          endif
          diff(:) = abs(ggoal-ratio(:))
          tthresh(:) = t(minloc(diff))
          if(DEBUG) then
              print *, 'ITERATION ', n_it
              print *, 'ratios = ', ratio
              print *, 'thresholds = ', t
              print *, 'threshold selected = ', tthresh
          endif
      end do
      if(present(thresh)) thresh = tthresh
      call sobel(img_in, bin_img, tthresh)
      if(DEBUG) print *, 'Final threshold = ', tthresh
      deallocate(grad)
  end subroutine automatic_thresh_sobel

  ! Canny edge detection.
  ! Performs canny edge detection on img_in and saves the result in img_out.
  ! PARAMETERS:
  !     -) 'thresh': optional.
  !                 If present the algorithm performs canny
  !                 with thresh as double threshold. If thresh
  !                 isn't present, then the double threshold for
  !                 is automatically selected.
  ! Note: 'sigma'.
  !        It is meaningful when the threshold has to be automatically
  !        selected. It determines 'how much' has to be detected.
  !        It could be changed by the user, but in this version of
  !        the algorithm it is fixed to its default value (= 0.33)
  !        as suggested in the source.
  ! If it doesn't work as you would like you can also check "Auto Canny" in GIT directory.
  subroutine canny(img_in, img_out, thresh)
      type(image),    intent(inout) :: img_in, img_out
      real, optional, intent(in)    :: thresh(2)
      real                       :: tthresh(2),  sigma, m
      real, allocatable          :: grad(:,:,:), vector(:)
      integer                    :: ldim(3)
      if(present(thresh)) then
         call canny_edge(img_in, img_out, thresh)
         return
      endif
      sigma = 0.33
      call img_in%scale_pixels([0.,255.])
      ldim = img_in%get_ldim()
      call img_in%calc_gradient(grad)
      allocate(vector(size(grad)), source = pack(grad, .true.))
      m = median(vector) !Use the gradient, the source talks about the image itself
      !https://www.pyimagesearch.com/2015/04/06/zero-parameter-automatic-canny-edge-detection-with-python-and-opencv/
      tthresh(1) = max(0.   ,(1-sigma)*m) !lower
      tthresh(2) = min(255., (1+sigma)*m) !upper
      call canny_edge(img_in, img_out, tthresh)
      deallocate(vector)
  end subroutine canny

  ! NON_MAX_SUPPRESSION
  ! using the estimates of the Gx and Gy image gradients and the edge direction angle
  ! determines whether the magnitude of the gradient assumes a local  maximum in the gradient direction
  subroutine non_max_supp(mat_in,dir_mat)
      real, intent(inout) :: mat_in(:,:,:)
      real,  intent(in)   :: dir_mat(:,:,:)
      real, allocatable   :: rmat(:,:,:), temp_mat(:,:,:)
      integer             :: ldim(3), i, j
      ldim = shape(mat_in)
      if(ldim(3) /= 1) THROW_HARD("The matrix has to be 2D!")
      !if the rounded edge direction angle is   0 degrees, checks the north and south directions
      !if the rounded edge direction angle is  45 degrees, checks the northwest and southeast directions
      !if the rounded edge direction angle is  90 degrees, checks the east and west directions
      !if the rounded edge direction angle is 135 degrees, checks the east and west directions
      allocate(temp_mat(ldim(1),ldim(2),1), source = mat_in)
      do i = 1, ldim(1)
          do j = 1, ldim(2)
              if( (pi/4. + pi/8. < dir_mat(i,j,1) .and. dir_mat(i,j,1) <= pi/2.) &         !case(90), north-south
              &     .or.(-pi/2. <= dir_mat(i,j,1) .and. (dir_mat(i,j,1) < -pi/4.-pi/8.) )) then
                      if(j+1 > ldim(2)) then
                          if(mat_in(i,j,1) < mat_in(i,j-1,1)) temp_mat(i,j,1) = 0.
                      else if(j-1 < 1) then
                          if(mat_in(i,j,1) < mat_in(i,j+1,1)) temp_mat(i,j,1) = 0.
                      else
                          if(mat_in(i,j,1) < mat_in(i,j+1,1) .or. mat_in(i,j,1) < mat_in(i,j-1,1)) temp_mat(i,j,1) = 0.
                      end if
              else if( -pi/8. <= dir_mat(i,j,1) .and. dir_mat(i,j,1) <= pi/8.) then       !case(0), west-east
                      if(i+1 > ldim(1)) then
                          if(mat_in(i,j,1) < mat_in(i-1,j,1)) temp_mat(i,j,1) = 0.
                      else if(i-1 < 1) then
                          if(mat_in(i,j,1) < mat_in(i+1,j,1)) temp_mat(i,j,1) = 0.
                      else
                          if(mat_in(i,j,1) < mat_in(i+1,j,1) .or. mat_in(i,j,1) < mat_in(i-1,j,1)) temp_mat(i,j,1) = 0.
                      endif
              else if(-pi/4.-pi/8. <= dir_mat(i,j,1) .and. dir_mat(i,j,1)< -pi/8.) then  !case(135) northeast - southwest
                      if( (i-1 < 1 .and. j-1 >0) .or.(j+1 > ldim(2) .and. i+1 <= ldim(1)) ) then
                          if(mat_in(i,j,1) < mat_in(i+1,j-1,1)) temp_mat(i,j,1) = 0.       !just northeast
                      else if( (i-1 >0  .and. j-1 < 1) .or. ( i+1 > ldim(1) .and. j+1 <= ldim(2)) ) then
                          if(mat_in(i,j,1) < mat_in(i-1,j+1,1)) temp_mat(i,j,1) = 0.       !just southwest
                      else if(i-1 > 0 .and. j+1 <= ldim(2) .and. i+1 <= ldim(1) .and. j-1 > 0) then
                          if(mat_in(i,j,1) < mat_in(i-1,j+1,1) .or. mat_in(i,j,1) < mat_in(i+1,j-1,1)) temp_mat(i,j,1) = 0. !northeast- southwest
                          !In the angles I decide not to to anything
                      endif
              else if(pi/8. < dir_mat(i,j,1) .and. dir_mat(i,j,1) <= pi/4. + pi/8.) then !case(45), northwest- southeast
                      if((j-1 < 1 .and. i+1 <= ldim(1))  .or. (i-1 < 1 .and. j+1 <= ldim(2))) then
                          if(mat_in(i,j,1) < mat_in(i+1,j+1,1)) temp_mat(i,j,1) = 0.        !just southeast
                      elseif((i+1 > ldim(1) .and. j-1 > 0) .or. (j+1 > ldim(2) .and. i-1 > 0)) then
                          if(mat_in(i,j,1) < mat_in(i-1,j-1,1)) temp_mat(i,j,1) = 0.    !just northwest
                      else if(i-1 > 0 .and. j+1 <= ldim(2) .and. i+1 <= ldim(1) .and. j-1 > 0) then
                          if(mat_in(i,j,1) < mat_in(i-1,j-1,1) .or. mat_in(i,j,1) < mat_in(i+1,j+1,1)) temp_mat(i,j,1) = 0. !northwest - southeast
                          !In the angles I decide not to to anything
                      endif
              else !case default
                  print *, "There is an error in the direction matrix"
              end if
          enddo
      enddo
      mat_in = temp_mat
      deallocate(temp_mat)
   end subroutine non_max_supp

  ! Edges stronger than thresh(2) are mantained, edges weaker than thresh(1)
  ! are discarded, in between edges are marked as "weak edges"
  subroutine double_thresh(mat_in, thresh)
    real, intent(inout) :: mat_in(:,:,:)
    real, intent(in)    :: thresh(2)
    real, allocatable :: tmp(:,:,:)
    integer           :: ldim(3), i, j
    ldim = shape(mat_in)
    allocate(tmp(ldim(1),ldim(2),1), source = 0.)
    do i = 1, ldim(1)
      do j = 1, ldim(2)
            if(mat_in(i,j,1) > thresh(2))                                    tmp(i,j,1) = 1.  !strong edges
            if(mat_in(i,j,1) < thresh(1))                                    tmp(i,j,1) = 0.  !not edges
            if(thresh(1) <= mat_in(i,j,1) .and. mat_in(i,j,1) <= thresh(2))  tmp(i,j,1) = 0.5 !weak edges
      enddo
    enddo
    mat_in = tmp
    deallocate(tmp)
  end subroutine double_thresh

  ! Canny routine for edge detection, in its classical verison.
  ! It requires a double threshold.
  subroutine canny_edge(img_in,img_out, thresh)
        type(image), intent(inout) :: img_in, img_out                               !input and output image
        real,        intent(in)    :: thresh(2)                                     !low and high thresholds
        type(image) :: Gr !gradient image
        integer     :: ldim(3), s(3), r(3), k, i, j                        !just for implementation
        real,    allocatable :: dir_mat(:,:,:), Dc(:,:,:), Dr(:,:,:), grad(:,:,:)   !derivates, gradient
        integer, allocatable :: neigh_8(:,:,:)                                      !8-neighborhoods of a pixel
        ldim = img_in%get_ldim()
        if(ldim(3) /= 1) THROW_HARD("The image has to be 2D!")
        if(img_out%exists()) call img_out%kill !reset if not empty
        call img_out%new(ldim,1.)
        call      Gr%new(ldim,1.)
        !STEP 1: SMOOTHING
        !Apply a Gaussian filter to reduce noise

        !STEP 2: FINDING GRADIENTS
        allocate(dir_mat(ldim(1),ldim(2),1), source = 0.)
        call img_in%calc_gradient(grad, Dc, Dr)
        do i = 1, ldim(1)
          do j = 1, ldim(2)
            if(Dc(i,j,1) /= 0.) then                     !Do not divide by 0
              dir_mat(i,j,1) = atan(Dr(i,j,1)/Dc(i,j,1)) !Output of ATAN is in radians
            else
              dir_mat(i,j,1) = pi/2.                     !My choice (limit)
            endif
          enddo
        enddo
        deallocate(Dc,Dr)
        !STEP 3: NON-MAXIMUM SUPPRESSION
        call non_max_supp(grad,dir_mat)
        deallocate(dir_mat)
        !STEP 4: DOUBLE THRESHOLDING
        call double_thresh(grad,thresh)
        !STEP 5: EDGE TRACKING BY HYSTERESIS
        call Gr%set_rmat(grad)
        do i = 1, ldim(1) !now grad is in {0,0.5,1}
            do j = 1, ldim(2)
                if( abs( Gr%get([i,j,1]) - 0.5 ) < TINY ) then
                    call Gr%set([i,j,1],0.)                         !suppress this edge, later I might re-build it
                    call Gr%calc_neigh_8([i,j,1], neigh_8)
                    s = shape(neigh_8)
                    do k = 1, s(2)
                        r = neigh_8(1:3,k,1)                                !indexes of 8 neighbourhoods
                        if(abs( Gr%get([r(1),r(2),1])- 1. ) < TINY) then    !one of the 8 neigh is a strong edge
                           call Gr%set([i,j,1], 1.)                         !re-build edge
                           exit
                        endif
                    enddo
                endif
            enddo
        enddo
        img_out = Gr
        deallocate(grad)
  end subroutine canny_edge
end module simple_edge_detector
