module simple_segmentation
include 'simple_lib.f08'
use simple_image, only: image
implicit none

public :: sobel, automatic_thresh_sobel, canny, iterative_thresholding, otsu_img, otsu_robust_fast, otsu_img_robust
private
#include "simple_local_flags.inc"

logical, parameter :: IMPROVED   = .true.
logical, parameter :: DOPRINT    = .false.

contains

    ! Classic Sobel edge detection routine
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
        real, optional,  intent(inout)   :: thresh(1) !to keep track of the selected threshold
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
    subroutine canny(img_in,img_out,thresh,lp)
        type(image),           intent(inout) :: img_in
        type(image), optional, intent(inout) :: img_out
        real,        optional, intent(in)    :: thresh(2)
        real,        optional, intent(in)    :: lp(1)
        type(image)       :: img_copy
        real              :: tthresh(2),  sigma, m
        real, allocatable :: grad(:,:,:), vector(:)
        integer           :: ldim(3)
        logical           :: debug
        real              :: smpd
        debug = .false.
        ldim  = img_in%get_ldim()
        smpd  = img_in%get_smpd()
        if(present(img_out)) then
            call img_out%copy(img_in)
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
            tthresh(1) = max(minval(grad), (1.-sigma)*m) !lower
            tthresh(2) = min(maxval(grad), (1.+sigma)*m) !upper
            write(logfhandle,*) 'Selected thresholds: ', tthresh
            if (present(lp)) then
                call canny_edge(img_out,tthresh,lp)
            else
                call canny_edge(img_out,tthresh)
            endif
        else
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
          allocate(vector(size(grad)), source = 0.)
          vector(:) = pack(grad, .true.)
          m = median_nocopy(vector) !Use the gradient, the source talks about the image itself
          !https://www.pyimagesearch.com/2015/04/06/zero-parameter-automatic-canny-edge-detection-with-python-and-opencv/
          tthresh(1) = max(minval(grad), (1.-sigma)*m) !lower
          tthresh(2) = min(maxval(grad), (1.+sigma)*m) !upper
          write(logfhandle,*) 'Selected thresholds: ', tthresh
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
        integer :: neigh_8(3,8) !8-neighborhoods of a pixel
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
        do i = 1, ldim(1)
            do j = 1, ldim(2)
                if(.not. is_zero(Dc(i,j,1))) then              !Do not divide by 0
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
                        r = neigh_8(1:3,k)                                !indexes of 8 neighbourhoods
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

    !This subroutine binarises an image through iterative thresholding.
    !ALGORITHM: -) select initial threshold T
    !           -) split the image in 2 groups according to T
    !           -) calc the mean grey level value in each group M1, M2
    !           -) update T as the mean value of M1 and M2
    !           -) stop when 2 threshold in a row are 'similar'
    subroutine iterative_thresholding(img_in,img_out,thresh)
        class(image),   intent(inout) :: img_in
        type(image),    intent(inout) :: img_out
        real, optional, intent(inout) :: thresh  !selected threshold for binarisation
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
            where(rmat <= T) mask = .true.                          ! in order to distinguish the 2 groups
            T_n = ((sum(rmat,     mask)/real(count(     mask))) + &
               & (sum(rmat, .not. mask)/real(count(.not. mask))))/2.! new threshold
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

    ! otsu binarization for images, based on the implementation
    ! of otsu algo for 1D vectors
    subroutine otsu_img(img, thresh, positive)
        use simple_math, only : otsu
        class(image),      intent(inout) :: img
        real,    optional, intent(inout) :: thresh
        logical, optional, intent(in)    :: positive ! consider just the positive range (specific for nanoparticles)
        real, pointer     :: rmat(:,:,:)
        real, allocatable :: x(:)
        integer           :: ldim(3)
        real :: selected_t
        ldim = img%get_ldim()
        call img%get_rmat_ptr(rmat)
        if(present(positive) .and. positive .eqv. .true.) then
            x = pack(rmat(1:ldim(1),1:ldim(2),1:ldim(3)), rmat(1:ldim(1),1:ldim(2),1:ldim(3)) > 0.)
        else
            x = pack(rmat(1:ldim(1),1:ldim(2),1:ldim(3)), .true.)
        endif
        call otsu(x, selected_t)
        if(present(thresh)) thresh= selected_t
        deallocate(x)
        where(rmat(1:ldim(1),1:ldim(2),1:ldim(3))>=selected_t)
            rmat(1:ldim(1),1:ldim(2),1:ldim(3)) = 1.
        elsewhere
            rmat(1:ldim(1),1:ldim(2),1:ldim(3)) = 0.
        endwhere
    end subroutine otsu_img

    ! Source : An Equivalent 3D Otsu’s Thresholding Method,
    ! Puthipong Sthitpattanapongsa and Thitiwan Srinark
    ! Idea: use 3 times otsu1D on the original vol,
    ! the avg vol and the median filtered vol.
    ! Find 3 thresholds and 3 corresponding bin vols.
    ! For each voxel select the most frequent value between
    ! background and foreground in the bin volumes.
    ! It is very robust to salt & pepper noise. Much more than
    ! otsu_img_robust
    subroutine otsu_robust_fast(img, is2D, noneg, thresh)
        type(image), intent(inout) :: img
        logical,     intent(in)    :: is2D   ! is it a 2D image
        logical,     intent(in)    :: noneg ! is it a nanoparticle
        real,        intent(inout) :: thresh(3)
        type(image) :: img_copy         ! copy of the original img
        type(image) :: img_avg, img_med
        integer :: ldim(3)
        integer :: i,j,k
        integer :: count_back          ! counts how many times a px is selected as background
        real    :: smpd
        real, pointer :: rmat(:,:,:)   ! original matrix
        real, pointer :: rmat_t(:,:,:) ! for the thresholded img
        real, pointer :: rmat_avg(:,:,:),rmat_med(:,:,:)
        ldim = img%get_ldim()
        smpd = img%get_smpd()
        if(is2D .eqv. .true.) ldim(3) = 1
        ! Initialise
        call img_copy%copy(img)
        call img_avg%new(ldim, smpd)
        call img_med%new(ldim, smpd)
        ! Generate avg and median volumes
        call generate_avg_and_median_vols(img,ldim,img_avg,img_med)
        ! Apply otsu1D on each img
        if(noneg .eqv. .true. ) then
            call otsu_img(img_copy, thresh(1), positive = .true.)
            call otsu_img(img_avg,  thresh(2), positive = .true.)
            call otsu_img(img_med,  thresh(3), positive = .true.)
        else
            call otsu_img(img_copy, thresh(1), positive = .false.)
            call otsu_img(img_avg,  thresh(2), positive = .false.)
            call otsu_img(img_med,  thresh(3), positive = .false.)
        endif
        ! Fetch pointers
        call img%get_rmat_ptr(rmat_t)
        call img_copy%get_rmat_ptr(rmat)
        call img_avg%get_rmat_ptr(rmat_avg)
        call img_med%get_rmat_ptr(rmat_med)
        ! Find most frequent value in the binary version of each voxel
        do i = 1, ldim(1)
            do j = 1, ldim(2)
                do k = 1, ldim(3)
                    count_back = 0
                    if(rmat(i,j,k)     < 0.5) count_back = count_back + 1
                    if(rmat_avg(i,j,k) < 0.5) count_back = count_back + 1
                    if(rmat_med(i,j,k) < 0.5) count_back = count_back + 1
                    if(count_back < 2) then
                        rmat_t(i,j,k) = 1.
                    else
                        rmat_t(i,j,k) = 0.
                    endif
                enddo
            enddo
        enddo
        ! kill
        call img_copy%kill
        call img_avg%kill
        call img_med%kill
    contains
        ! This subroutine generates 2 images:
        ! 1) img_avg: image where in each pixel is contained
        ! the average value of the gray levels in the neighbours of the pixel.
        ! 2) img_med: image where in each pixel is contained
        ! the median value of the gray levels in the neighbours of the pixel.
        subroutine generate_avg_and_median_vols(img,ldim,img_avg,img_med)
            type(image), intent(inout) :: img, img_avg, img_med
            integer,     intent(in)    :: ldim(3)
            integer :: i, j, k
            integer :: nsz
            real, allocatable :: neigh_8(:)
            real, pointer     :: rmat_avg(:,:,:),rmat_med(:,:,:)
            call img_avg%get_rmat_ptr(rmat_avg)
            call img_med%get_rmat_ptr(rmat_med)
            if(ldim(3) == 1) then
                allocate(neigh_8(9))
            else
                allocate(neigh_8(27))
            endif
            do i = 1, ldim(1)
                do j = 1, ldim(2)
                    do k = 1, ldim(3)
                        if(ldim(3) == 1 ) then
                            call img%calc_neigh_8([i,j,k],neigh_8,nsz)
                        else
                            call img%calc3D_neigh_8([i,j,k],neigh_8,nsz)
                        endif
                        rmat_avg(i,j,k) = sum(neigh_8(1:nsz))/real(nsz)
                        rmat_med(i,j,k) = median(neigh_8(1:nsz))
                    enddo
                enddo
            enddo
            deallocate(neigh_8)
        end subroutine generate_avg_and_median_vols
    end subroutine otsu_robust_fast

    ! otsu binarization for images, based on the implementation
    ! of otsu algo for 2D arrays. It should be better in case
    ! of salt and pepper noise. It is computationally expensive.
    ! Don't use it unless necessary
    subroutine otsu_img_robust(img, thresh)
        type(image),    intent(inout) :: img
        real, optional, intent(inout) :: thresh
        integer, parameter   :: MIN_VAL = 0, MAX_VAL = 255
        integer, allocatable :: yhist_orig(:)
        real,    allocatable :: xhist_orig(:)
        real,    allocatable :: x_orig(:) ! vectorization of the img
        real,    pointer     :: rmat(:,:,:), rmat_avg(:,:,:)
        integer     :: i, j, k, h, ldim(3), S, T, ss, tt
        real        :: hists(MAX_VAL-MIN_VAL+1,MAX_VAL-MIN_VAL+1), smpd, tr, max
        real        :: sc, old_range(2) ! scale factor and old pixel range
        type(image) :: img_avg
        if( DOPRINT ) write(logfhandle,*) '*** initialising otsu'
        ldim = img%get_ldim()
        smpd = img%get_smpd()
        if(ldim(3) > 1) THROW_HARD('Not implemented for volumes! otsu_img_robust')
        call img%scale_pixels(real([MIN_VAL,MAX_VAL]), sc, old_range)
        call img%get_rmat_ptr(rmat)
        x_orig = pack(rmat(1:ldim(1),1:ldim(2),1), .true.)
        call create_hist_vector(x_orig, MAX_VAL-MIN_VAL+1, xhist_orig, yhist_orig)
        deallocate(yhist_orig)
        call img_avg%new(ldim, smpd)
        call generate_neigh_avg_mat(img, img_avg)
        call img_avg%get_rmat_ptr(rmat_avg)
        if( DOPRINT ) write(logfhandle,*)  '*** hists generation'
        hists = 0.
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,h,k)&
        !$omp proc_bind(close) reduction(+:hists)
        do i = 1, ldim(1)
            do j = 1, ldim(2)
                do k = MIN_VAL, MAX_VAL
                    do h = MIN_VAL, MAX_VAL
                        if(k == MAX_VAL .and. h == MAX_VAL) then
                            if(rmat(i,j,1)    >= xhist_orig(k) .and. &
                               rmat_avg(i,j,1)>= xhist_orig(h)) then
                                   hists(k+1,h+1) = hists(k+1,h+1) + 1.
                            endif
                            cycle
                        endif
                        if(k == MAX_VAL) then
                            if(rmat(i,j,1)    >= xhist_orig(k) .and. &
                               rmat_avg(i,j,1)>= xhist_orig(h) .and. rmat_avg(i,j,1) < xhist_orig(h+1)) then
                                   hists(k+1,h+1) = hists(k+1,h+1) + 1.
                             endif
                            cycle
                        endif
                        if(h == MAX_VAL) then
                            if(rmat(i,j,1)    >= xhist_orig(k) .and. rmat(i,j,1)     < xhist_orig(k+1) .and. &
                               rmat_avg(i,j,1)>= xhist_orig(h)) then
                                   hists(k+1,h+1) = hists(k+1,h+1) + 1.
                             endif
                            cycle
                        endif
                        if(rmat(i,j,1)    >= xhist_orig(k) .and. rmat(i,j,1)        < xhist_orig(k+1) .and. &
                           rmat_avg(i,j,1)>= xhist_orig(h) .and. rmat_avg(i,j,1)    < xhist_orig(h+1)) then
                               hists(k+1,h+1) = hists(k+1,h+1) + 1.
                        endif
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do
        ! hists is the joint probability, normalise
        hists = hists/real(ldim(1)*ldim(2))
        if( DOPRINT ) write(logfhandle,*) 'sum ', sum(hists)
        ! Threshold selection
        if( DOPRINT ) write(logfhandle,*) '*** selecting thresholds'
        call calc_tr(hists,MIN_VAL,MAX_VAL,threshold = tr)
        where(rmat >= tr)
            rmat = 1.
        elsewhere
            rmat = 0.
        endwhere
        if(present(thresh))  thresh =  tr/sc+old_range(1) !rescale threshold in the old range
        call img_avg%kill
        if( DOPRINT ) write(logfhandle,*) '*** binarization successfully completed'

        contains

            subroutine calc_tr(hists,MIN_VAL,MAX_VAL,threshold)
                real,    intent(inout) :: hists(:,:)
                integer, intent(in)    :: MIN_VAL,MAX_VAL
                real,    intent(inout) :: threshold
                real :: tr, maximum
                real :: p_0(MAX_VAL-MIN_VAL+1,MAX_VAL-MIN_VAL+1)
                integer :: i, j
                real :: mu_T(2)
                real :: mu_k(MAX_VAL-MIN_VAL+1,MAX_VAL-MIN_VAL+1), mu_h(MAX_VAL-MIN_VAL+1,MAX_VAL-MIN_VAL+1)
                mu_T(:) = 0.
                !$omp parallel do collapse(2) schedule(static) default(shared) private(k,h)&
                !$omp reduction(+:mu_T) proc_bind(close)
                do i = MIN_VAL+1, MAX_VAL+1
                    do j = MIN_VAL+1, MAX_VAL+1
                        mu_T(1) = mu_T(1) + (real(i)*hists(i,j))
                        mu_T(2) = mu_T(2) + (real(j)*hists(i,j))
                    enddo
                enddo
                !$omp end parallel do
                maximum   = 0.
                threshold = 0.
                p_0     = 0.
                mu_k    = 0.
                mu_h    = 0.
                do i = MIN_VAL+1, MAX_VAL + 1
                    do j = MIN_VAL+1, MAX_VAL + 1
                        if(j == 1) then
                            if(i == 1) then
                                p_0(1,1) = hists(1,1)
                            else
                                p_0(i,1)  = p_0(i-1,1)+hists(i,1)
                                mu_k(i,1) = mu_k(i-1,1)+real(i-1)*hists(i,1)
                                mu_h(i,1) = mu_h(i-1,1)
                            endif
                        elseif(i == 1) then
                            p_0(i,j)  = p_0(i,j-1)+hists(i,j)
                            mu_k(i,j) = mu_k(i,j-1)
                            mu_h(i,j) = mu_h(i,j-1)+real(j-1)*hists(i,j)
                        else
                             p_0(i,j)  = p_0(i,j-1)+p_0(i-1,j)-p_0(i-1,j-1)+hists(i,j)
                             mu_k(i,j) = mu_k(i,j-1)+mu_k(i-1,j)-mu_k(i-1,j-1)+real(i-1)*hists(i,j)
                             mu_h(i,j) = mu_h(i,j-1)+mu_h(i-1,j)-mu_h(i-1,j-1)+real(j-1)*hists(i,j)
                         endif
                         if(abs(p_0(i,j)) < TINY) cycle
                         if(p_0(i,j)-0.99 > 0.)  exit
                         if(i == 2 .and. j == 2) exit
                         tr = ((mu_k(i,j)-p_0(i,j)*mu_T(1))**2+(mu_h(i,j)-p_0(i,j)*mu_T(2))**2)/(p_0(i,j)*(1.-p_0(i,j)))
                         if(tr >= maximum) then
                             threshold = real(i)
                             maximum = tr
                         endif
                    enddo
                enddo
            end subroutine calc_tr

            ! This subroutine generates an image where in each pixel is contained
            ! the average value of the gray levels in the neighbours of the pixel.
            subroutine generate_neigh_avg_mat(img, img_avg)
                type(image), intent(inout) :: img, img_avg
                integer :: ldim(3) !shape of the matrix
                integer :: i, j
                integer :: nsz
                real    :: neigh_8(9)
                real, pointer :: rmat(:,:,:)
                call img_avg%get_rmat_ptr(rmat)
                ldim = img%get_ldim()
                do i = 1, ldim(1)
                    do j = 1, ldim(2)
                        call img%calc_neigh_8([i,j,1],neigh_8,nsz)
                        rmat(i,j,1) = sum(neigh_8(1:nsz))/real(nsz)
                    enddo
                enddo
            end subroutine generate_neigh_avg_mat

    end subroutine otsu_img_robust

end module simple_segmentation
