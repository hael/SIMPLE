module simple_segmentation
include 'simple_lib.f08'
use simple_image, only: image
implicit none

logical, parameter :: IMPROVED = .true.
public :: sobel, automatic_thresh_sobel, canny, iterative_thresholding, otsu_img

private
#include "simple_local_flags.inc"
contains

    ! Classic Sobel edge detection routine.
    subroutine sobel(img_in,thresh)
        class(image), intent(inout) :: img_in           !image input and output
        real,         intent(in)    :: thresh(1)        !threshold for Sobel algorithm
        real,  allocatable :: grad(:,:,:), rmat(:,:,:)  !matrices
        integer :: ldim(3)
        ldim = img_in%get_ldim()
        allocate(grad(ldim(1),ldim(2),ldim(3)), rmat(ldim(1),ldim(2),ldim(3)), source = 0.)
        if(IMPROVED) then
            call img_in%calc_gradient_improved(grad)
        else
            call img_in%calc_gradient(grad)
        endif
        where(grad > thresh(1) ) rmat = 1.
        call img_in%set_rmat(rmat)
        deallocate(grad,rmat)
    end subroutine sobel

    ! This soubroutine performs sobel edge detection with an automatic iterative
    ! threshold selection in order to obtian a threshold which leads to
    ! a ratio between foreground and background pixels in the binary image
    ! close to 'goal'. This particular number is selected because it is what needed
    ! for Chiara's particle picking routine.
    subroutine automatic_thresh_sobel(img, goal, thresh)
        class(image),    intent(inout)   :: img
        real, optional,  intent(in)      :: goal
        real, optional,  intent(out)     :: thresh(1) !to keep track of the selected threshold
        real, allocatable  :: grad(:,:,:)
        real               :: tthresh(1), ggoal       !selected threshold and ratio nforeground/nbackground
        real               :: t(3)                    !to calculate 3 possible sobel thresholds for each iteration
        real               :: ratio(3)                !n_foreground/nbackground for 3 different values of the thresh
        real               :: diff(3)                 !difference between the ideal ratio (= goal) and the obtained 3 ratios
        integer, parameter :: N_MAXIT = 100           !max number of iterations
        integer            :: n_it, ldim(3)
        logical            :: DEBUG
        real               :: r                       !random # which guarantees me not to get stuck
        real               :: ave, sdev, minv, maxv
        DEBUG  = .false.
        ggoal = 2.
        ldim = img%get_ldim()
        allocate(grad(ldim(1),ldim(2),ldim(3)), source = 0.)
        if(present(goal)) ggoal = goal
        if(IMPROVED) then
            call img%calc_gradient_improved(grad)
        else
            call img%calc_gradient(grad)
        endif
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
            call sobel(img,t(1))
            if(abs(real(img%nbackground())) > TINY) then !do not divide by 0
                ratio(1) = real(img%nforeground())/real(img%nbackground())  !obtained ratio with threshold = t(1)
            else
                ratio(1) = 0.
            endif
            call sobel(img,t(2))
            if(abs(real(img%nbackground())) > TINY) then !do not divide by 0
                ratio(2) = real(img%nforeground())/real(img%nbackground())  !obtained ratio with threshold = t(2)
            else
                ratio(2) = 0.
            endif
            call sobel(img,t(3))
            if(abs(real(img%nbackground())) > TINY) then !do not divide by 0
                ratio(3) = real(img%nforeground())/real(img%nbackground())  !obtained ratio with threshold = t(3)
            else
                ratio(3) = 0.
            endif
            diff(:) = abs(ggoal-ratio(:))
            tthresh(:) = t(minloc(diff))
            if(DEBUG) then
                write(logfhandle,*) 'ITERATION ', n_it
                write(logfhandle,*) 'ratios = ', ratio
                write(logfhandle,*) 'thresholds = ', t
                write(logfhandle,*) 'threshold selected = ', tthresh
            endif
        end do
        if(present(thresh)) thresh = tthresh
        call sobel(img, tthresh)
        if(DEBUG) write(logfhandle,*) 'Final threshold = ', tthresh
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
  subroutine canny(img_in,img_out,scale_range,thresh,lp)
      type(image),           intent(inout) :: img_in
      type(image), optional, intent(out)   :: img_out
      real, optional, intent(in)    :: scale_range(2)
      real, optional, intent(in)    :: thresh(2)
      real, optional, intent(in)    :: lp(1)
      real              :: tthresh(2),  sigma, m
      real, allocatable :: grad(:,:,:), vector(:)
      integer           :: ldim(3)
      logical           :: debug
      real              :: sscale_range(2)
      real              :: smpd
      debug = .false.
      ldim  = img_in%get_ldim()
      smpd  = img_in%get_smpd()
      sscale_range = [1.,real(ldim(1))] !real(ldim(1))*2.+1., modification 25/07
      if(present(scale_range)) sscale_range = scale_range
      if(present(img_out)) then
          call img_out%new(ldim, smpd)
          call img_out%copy(img_in)
          call img_out%scale_pixels([sscale_range(1),sscale_range(2)])
          if(present(thresh)) then
                 call canny_edge(img_out,thresh)
                 return
          endif
          sigma = 0.33
          allocate(grad(ldim(1),ldim(2),ldim(3)), source = 0.)
          if(IMPROVED) then
              call img_out%calc_gradient_improved(grad)
          else
              call img_out%calc_gradient(grad)
          endif
          allocate(vector(size(grad)), source = 0.)
          vector(:) = pack(grad, .true.)
          m = median_nocopy(vector) !Use the gradient, the source talks about the image itself
          !https://www.pyimagesearch.com/2015/04/06/zero-parameter-automatic-canny-edge-detection-with-python-and-opencv/
          tthresh(1) = max(sscale_range(1), (1.-sigma)*m) !lower
          tthresh(2) = min(sscale_range(2), (1.+sigma)*m) !upper
          if(debug) write(logfhandle,*) 'Selected thresholds: ', tthresh
          if (present(lp)) then
              call canny_edge(img_out,tthresh,lp)
          else
              call canny_edge(img_out,tthresh)
          endif
      else
          call img_in%scale_pixels([sscale_range(1),sscale_range(2)]) !to be coherent with the other case
          if(present(thresh)) then
                 call canny_edge(img_in,thresh)
                 return
          endif
          sigma = 0.33
          allocate(grad(ldim(1),ldim(2),ldim(3)), source = 0.)
          if(IMPROVED) then
              call img_in%calc_gradient_improved(grad)
          else
              call img_in%calc_gradient(grad)
          endif
          allocate(vector(size(grad)), source = 0.)!= pack(grad, .true.))
          vector(:) = pack(grad, .true.)
          m = median_nocopy(vector) !Use the gradient, the source talks about the image itself
          !https://www.pyimagesearch.com/2015/04/06/zero-parameter-automatic-canny-edge-detection-with-python-and-opencv/
          tthresh(1) = max(sscale_range(1), (1.-sigma)*m) !lower
          tthresh(2) = min(sscale_range(2), (1.+sigma)*m) !upper
          if(debug) write(logfhandle,*) 'Selected thresholds: ', tthresh
          if (present(lp)) then
              call canny_edge(img_in,tthresh,lp)
          else
              call canny_edge(img_in,tthresh)
          endif
      endif
      deallocate(vector,grad)
  end subroutine canny

    ! NON_MAX_SUPPRESSION
    ! using the estimates of the Gx and Gy image gradients and the edge direction angle
    ! determines whether the magnitude of the gradient assumes a local  maximum in the gradient direction
    subroutine non_max_supp(mat_in,dir_mat)
        real,  intent(inout) :: mat_in(:,:,:)
        real,  intent(in)    :: dir_mat(:,:,:)
        real, allocatable    :: temp_mat(:,:,:)
        integer              :: ldim(3), i, j
        ldim = shape(mat_in)
        !if the rounded edge direction angle is   0 degrees, checks the north and south directions
        !if the rounded edge direction angle is  45 degrees, checks the northwest and southeast directions
        !if the rounded edge direction angle is  90 degrees, checks the east and west directions
        !if the rounded edge direction angle is 135 degrees, checks the northeast and southwest directions
        allocate(temp_mat(ldim(1),ldim(2),ldim(3)), source = mat_in)
        do i = 1, ldim(1)
            do j = 1, ldim(2)
                if( (pi/4. + pi/8. < dir_mat(i,j,1) .and. dir_mat(i,j,1) <= pi/2.) &  !case(90), north-south
                    & .or.(-pi/2. <= dir_mat(i,j,1) .and. (dir_mat(i,j,1) < -pi/4.-pi/8.) )) then
                    if(j+1 > ldim(2)) then
                        if(mat_in(i,j,1) < mat_in(i,j-1,1)) temp_mat(i,j,1) = 0.
                    else if(j-1 < 1) then
                        if(mat_in(i,j,1) < mat_in(i,j+1,1)) temp_mat(i,j,1) = 0.
                    else
                        if(mat_in(i,j,1) < mat_in(i,j+1,1) .or. mat_in(i,j,1) < mat_in(i,j-1,1)) temp_mat(i,j,1) = 0.
                    end if
                else if( -pi/8. <= dir_mat(i,j,1) .and. dir_mat(i,j,1) <= pi/8.) then !case(0), west-east
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
                    write(logfhandle,*) "There is an error in the direction matrix"
                end if
            enddo
        enddo
        mat_in(:ldim(1),:ldim(2),:ldim(3)) = temp_mat(:ldim(1),:ldim(2),:ldim(3))
        deallocate(temp_mat)
    end subroutine non_max_supp

    ! Edges stronger than thresh(2) are mantained, edges weaker than thresh(1)
    ! are discarded, in between edges are marked as "weak edges"
    subroutine double_thresh(mat_in, thresh)
        real, intent(inout) :: mat_in(:,:,:)
        real, intent(in)    :: thresh(2)
        real, allocatable   :: tmp(:,:,:)
        integer             :: ldim(3), i, j
        ldim = shape(mat_in)
        allocate(tmp(ldim(1),ldim(2),ldim(3)), source = 0.)
        do i = 1, ldim(1)
            do j = 1, ldim(2)
                if(mat_in(i,j,1) > thresh(2))                                    tmp(i,j,1) = 1.  !strong edges
                if(mat_in(i,j,1) < thresh(1))                                    tmp(i,j,1) = 0.  !not edges
                if(thresh(1) <= mat_in(i,j,1) .and. mat_in(i,j,1) <= thresh(2))  tmp(i,j,1) = 0.5 !weak edges
            enddo
        enddo
        mat_in(:ldim(1),:ldim(2),:ldim(3)) = tmp(:ldim(1),:ldim(2),:ldim(3))
        deallocate(tmp)
    end subroutine double_thresh

  ! Canny routine for edge detection, in its classical verison.
  ! It requires a double threshold.
  subroutine canny_edge(img_in, thresh, lp)
        type(image),    intent(inout) :: img_in     !input and output image
        real,           intent(in)    :: thresh(2)  !low and high thresholds
        real, optional, intent(in)    :: lp(1)        !lp filtering
        type(image) :: Gr !gradient image
        integer     :: ldim(3), s(3), r(3), k, i, j, nsz                   !just for implementation
        real,    allocatable :: dir_mat(:,:,:), Dc(:,:,:), Dr(:,:,:), grad(:,:,:)   !derivates, gradient
        integer :: neigh_8(3,8,1) !8-neighborhoods of a pixel
        real    :: smpd, llp(1)
        ldim = img_in%get_ldim()
        smpd = img_in%get_smpd()
        allocate(grad(ldim(1),ldim(2),ldim(3)), Dc(ldim(1),ldim(2),ldim(3)), Dr(ldim(1),ldim(2),ldim(3)), &
        & source = 0.)
        if(ldim(3) /= 1) THROW_HARD("The image has to be 2D!")
        call Gr%new(ldim,smpd)
        !STEP 1: SMOOTHING
        !Apply a lp to reduce noise (traditionally it would be Gaussian Filter)
        llp = 10.
        if(present(lp)) llp = lp
        call img_in%bp(0.,llp(1))   !threshold selected in order to facilitate particle picking
        call img_in%ifft()
        !STEP 2: FINDING GRADIENTS
        allocate(dir_mat(ldim(1),ldim(2),ldim(3)), source = 0.)
        if(IMPROVED) then
            call img_in%calc_gradient_improved(grad, Dc, Dr)
        else
            call img_in%calc_gradient(grad, Dc, Dr)
        endif
        !!!!!!!!!!!!!!!!!!
        do i = 1, ldim(1)
            do j = 1, ldim(2)
                if(.not. is_zero(Dc(i,j,1))) then                     !Do not divide by 0
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
                    call Gr%calc_neigh_8([i,j,1], neigh_8, nsz)
                    do k = 1, nsz
                        r = neigh_8(1:3,k,1)                                !indexes of 8 neighbourhoods
                        if(is_equal( Gr%get([r(1),r(2),1]), 1.) ) then      !one of the 8 neigh is a strong edge
                            call Gr%set([i,j,1], 1.)                         !re-build edge
                            exit
                        endif
                    enddo
                endif
            enddo
        enddo
        call img_in%copy(Gr)
        deallocate(grad)
    end subroutine canny_edge
  !
    !This subroutine binarises an image through iterative thresholding.
    !ALGORITHM: -) select initial threshold T
    !           -) split the image in 2 groups according to T
    !           -) calc the mean grey level value in each group M1, M2
    !           -) update T as the mean value of M1 and M2
    !           -) stop when 2 threshold in a row are 'similar'
    subroutine iterative_thresholding(img_in,img_out,thresh)
      class(image),   intent(inout) :: img_in
      type(image),    intent(out)   :: img_out
      real, optional, intent(out)   :: thresh  !selected threshold for binarisation
      real,    allocatable :: rmat(:,:,:)
      logical, allocatable :: mask(:,:,:)
      real, parameter :: TOLL = 0.05           !tollerance for threshold selection
      real            :: T_n, T
      integer         :: cnt
      call img_out%new(img_in%get_ldim(),img_in%get_smpd())
      rmat = img_in%get_rmat()
      allocate(mask(size(rmat, dim =1), size(rmat, dim =1), size(rmat, dim = 3)), source = .false.)
      T = (maxval(rmat) + minval(rmat))/2.    !initial guessing
      T_n = T + 1                             !to get into the while loop
      cnt = 0                                 !# of iterations
      do while(abs(T-T_n)> TOLL)
          cnt = cnt + 1
          if(cnt > 0)  T = T_n
          where(rmat <= T) mask = .true.                            !in order to distinguish the 2 groups
          T_n = ((sum(rmat,       mask)/real(count(     mask))) + &
               & (sum(rmat, .not. mask)/real(count(.not. mask))))/2.!new threshold
      enddo
      if(present(thresh)) thresh = T_n
      where(rmat > T_n)  !thresholding
          rmat = 1.
      elsewhere
          rmat = 0.
      endwhere
      call img_out%set_rmat(rmat)
      deallocate(rmat, mask)
    end subroutine iterative_thresholding

   !Otsu binarization for images
    subroutine otsu_img(img, thresh)
        use simple_math, only : otsu
        class(image),   intent(inout) :: img
        real, optional, intent(out)   :: thresh
        real, pointer     :: rmat(:,:,:)
        real, allocatable :: x(:)
        integer           :: ldim(3)
        real :: selected_t
        ldim = img%get_ldim()
        call img%get_rmat_ptr(rmat)
        x = pack(rmat(1:ldim(1),1:ldim(2),1:ldim(3)), .true.)
        call otsu(x, selected_t)
        if(present(thresh)) thresh= selected_t
        deallocate(x)
        where(rmat(1:ldim(1),1:ldim(2),1:ldim(3))>=selected_t)
            rmat(1:ldim(1),1:ldim(2),1:ldim(3)) = 1.
        elsewhere
            rmat(1:ldim(1),1:ldim(2),1:ldim(3)) = 0.
        endwhere
    end subroutine otsu_img
end module simple_segmentation
