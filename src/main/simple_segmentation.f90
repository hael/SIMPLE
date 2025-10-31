module simple_segmentation
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,    only: image
use simple_binimage, only: binimage
use simple_neighs
implicit none
#include "simple_local_flags.inc"

interface detect_peak_thres
    module procedure detect_peak_thres_1, detect_peak_thres_2
end interface detect_peak_thres

contains

    subroutine detect_peak_thres_1( n, n_ub, level, x, t )
        integer, intent(in)    :: n, n_ub, level
        real,    intent(in)    :: x(n)
        real,    intent(inout) :: t
        real,    allocatable   :: arr(:)
        integer, allocatable   :: locn(:)
        real    :: ts(2)
        if( level < 0 .or. level > 2 )then
            call simple_exception('peak detection level out of range', 'simple_math.f90', __LINE__,l_stop=.true.)
        endif
        locn  = maxnloc(x, n_ub)
        ts(1) = minval(x(locn))
        t     = ts(1)
        if( level == 0 ) return
        arr   = pack(x, mask=x >= ts(1))
        call otsu(size(arr), arr, ts(2))
        t     = ts(2)
        if( level == 1 ) return
        ts(1) = ts(2)
        arr   = pack(x, mask=x >= ts(1))
        call otsu(size(arr), arr, ts(2))
        t     = ts(2)
    end subroutine detect_peak_thres_1

    subroutine detect_peak_thres_2( n, level, x, t )
        integer, intent(in)    :: n, level
        real,    intent(in)    :: x(n)
        real,    intent(inout) :: t
        real,    allocatable   :: arr(:)
        real    :: ts(2)
        if( level < 1 .or. level > 2 )then
            call simple_exception('peak detection level out of range', 'simple_math.f90', __LINE__,l_stop=.true.)
        endif
        arr   = pack(x, mask=.true.)
        call otsu(size(arr), arr, ts(1))
        t     = ts(1)
        if( level == 1 ) return
        arr   = pack(x, mask=x >= ts(1))
        call otsu(size(arr), arr, ts(2))
        t     = ts(2)
    end subroutine detect_peak_thres_2

    subroutine detect_peak_thres_for_npeaks( n, npeaks, x, t )
        integer, intent(in)    :: n, npeaks
        real,    intent(in)    :: x(n)
        real,    intent(inout) :: t 
        real, allocatable :: tmp(:)
        integer :: ind
        allocate(tmp(n), source=x)
        call hpsort(tmp)
        ind = max(1,n - npeaks + 1)
        t   = tmp(ind)
    end subroutine detect_peak_thres_for_npeaks

    subroutine detect_peak_thres_sortmeans( n, level, x, t_out )
        integer, intent(in)  :: n, level
        real,    intent(in)  :: x(n)
        real,    intent(out) :: t_out
        integer, parameter   :: NQUANTA = 10
        real,    allocatable :: arr(:), means(:), peak_ts(:), frac_peaks(:), arr1(:), arr2(:)
        integer, allocatable :: labels(:)
        integer :: iq, n_fg, nvals, cnt, iq_min, ind
        real    :: diff, diff_min, med1, med2
        ! quantize with sortmeans
        allocate( means(NQUANTA), labels(NQUANTA), frac_peaks(NQUANTA), peak_ts(NQUANTA) )
        means  = 0.
        labels = 0
        arr    = pack(x, mask=.true.)
        nvals  = size(arr)
        call sortmeans(arr, NQUANTA, means, labels)
        cnt    = 0
        do iq = 1, NQUANTA
            if( count(labels == iq) == 0 )cycle
            cnt          = cnt + 1
            peak_ts(cnt) = minval(arr, mask=labels == iq)
            n_fg         = count(arr >= peak_ts(cnt))
            if( n_fg == nvals )then
                frac_peaks(cnt) = 1.
            else if( n_fg == 0 )then
                frac_peaks(cnt) = 0.
            else
                frac_peaks(cnt) = real(n_fg) / real(nvals)
            endif
        end do
        ! binary clustering for dynamic threshold determination
        ! init: for iq=1, arr1=arr & arr2 is empty
        iq_min   = 1
        med1     = median_nocopy(arr)
        diff_min = sum(abs(arr - med1))
        do iq = 2, cnt
            arr1 = pack(arr, mask=arr >= peak_ts(iq))
            arr2 = pack(arr, mask=arr <  peak_ts(iq))
            med1 = median_nocopy(arr1)
            med2 = median_nocopy(arr2)
            diff = sum(abs(arr - med1), mask=arr >= peak_ts(iq)) + sum(abs(arr - med2), mask=arr < peak_ts(iq))
            if( diff < diff_min )then
                diff_min = diff
                iq_min   = iq
            endif
            deallocate(arr1,arr2)
        end do
        ! apply particle crowding level adjustement
        select case(level)
            case(1)
                ind = iq_min              ! fewer peaks
            case(2)
                ind = max(1,  iq_min - 1) ! optimal # peaks
            case(3)
                ind = max(1,  iq_min - 2) ! more peaks
            case DEFAULT
                ind = iq_min
        end select
        t_out = peak_ts(ind)
    end subroutine detect_peak_thres_sortmeans

    ! classic Sobel edge detection routine
    subroutine sobel( img_in, thresh )
        class(image), intent(inout) :: img_in    ! image input and output
        real,         intent(in)    :: thresh(1) ! threshold for Sobel algorithm
        real,  allocatable :: grad(:,:,:), rmat(:,:,:)
        integer :: ldim(3)
        ldim = img_in%get_ldim()
        allocate(grad(ldim(1),ldim(2),ldim(3)), rmat(ldim(1),ldim(2),ldim(3)), source=0.)
        call img_in%calc_gradient(grad)
        where(grad > thresh(1) ) rmat = 1.
        call img_in%set_rmat(rmat,.false.)
        deallocate(grad,rmat)
    end subroutine sobel

    ! This soubroutine performs sobel edge detection with an automatic iterative
    ! threshold selection in order to obtain a threshold which leads to
    ! a ratio between foreground and background pixels in the binary image
    ! close to 'goal'. This particular number is selected because it is what needed
    ! for Chiara's particle picking routine.
    subroutine automatic_thresh_sobel( img, goal, thresh )
        class(image), intent(inout)    :: img
        real, optional,  intent(in)    :: goal
        real, optional,  intent(inout) :: thresh(1) !to keep track of the selected threshold
        real, allocatable  :: grad(:,:,:)
        real               :: tthresh(1), ggoal       !selected threshold and ratio nforeground/nbackground
        real               :: t(3)                    !to calculate 3 possible sobel thresholds for each iteration
        real               :: ratio(3)                !n_foreground/nbackground for 3 different values of the thresh
        real               :: diff(3)                 !difference between the ideal ratio (= goal) and the obtained 3 ratios
        integer, parameter :: N_MAXIT = 100           !max number of iterations
        integer            :: n_it, ldim(3)
        logical            :: DEBUG
        real               :: r                       !random # which guarantees me not to get stuck
        DEBUG  = .false.
        ggoal = 2.
        ldim = img%get_ldim()
        allocate(grad(ldim(1),ldim(2),ldim(3)), source = 0.)
        if(present(goal)) ggoal = goal
        call img%calc_gradient(grad)
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
    subroutine canny( img_in, img_out, thresh, lp )
        class(image),           intent(inout) :: img_in
        class(image), optional, intent(inout) :: img_out
        real,         optional, intent(in)    :: thresh(2)
        real,         optional, intent(in)    :: lp(1)
        real              :: tthresh(2),  sigma, m
        real, allocatable :: grad(:,:,:), vector(:)
        integer           :: ldim(3)
        logical           :: debug
        real              :: smpd
        debug = .false.
        ldim  = img_in%get_ldim()
        smpd  = img_in%get_smpd()
        if( present(img_out) )then
            call img_out%copy(img_in)
            if(present(thresh)) then
                call canny_edge(img_out, thresh)
                return
            endif
            sigma = 0.33
            allocate(grad(ldim(1),ldim(2),ldim(3)), source = 0.)
            call img_out%calc_gradient(grad)
            allocate(vector(size(grad)), source = 0.)
            vector(:) = pack(grad, .true.)
            m = median_nocopy(vector) !Use the gradient, the source talks about the image itself
            !https://www.pyimagesearch.com/2015/04/06/zero-parameter-automatic-canny-edge-detection-with-python-and-opencv/
            tthresh(1) = max(minval(grad), (1.-sigma)*m) !lower
            tthresh(2) = min(maxval(grad), (1.+sigma)*m) !upper
            if( present(lp) ) then
                call canny_edge(img_out, tthresh, lp)
            else
                call canny_edge(img_out, tthresh)
            endif
        else
            if( present(thresh) )then
                call canny_edge(img_in, thresh)
                return
            endif
            sigma = 0.33
            allocate(grad(ldim(1),ldim(2),ldim(3)), source=0.)
            call img_in%calc_gradient(grad)
            allocate(vector(size(grad)), source=0.)
            vector(:) = pack(grad, .true.)
            m = median_nocopy(vector) !Use the gradient, the source talks about the image itself
            !https://www.pyimagesearch.com/2015/04/06/zero-parameter-automatic-canny-edge-detection-with-python-and-opencv/
            tthresh(1) = max(minval(grad), (1.-sigma)*m) !lower
            tthresh(2) = min(maxval(grad), (1.+sigma)*m) !upper
            write(logfhandle,*) 'Selected thresholds: ', tthresh
            if( present(lp) )then
                call canny_edge(img_in, tthresh, lp)
            else
                call canny_edge(img_in, tthresh)
            endif
        endif
        deallocate(vector,grad)
    end subroutine canny

    ! NON_MAX_SUPPRESSION
    ! using the estimates of the Gx and Gy image gradients and the edge direction angle
    ! determines whether the magnitude of the gradient assumes a local maximum in the gradient direction
    subroutine non_max_supp( mat_in, dir_mat )
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
                        !In the angles I decide not to to do anything
                    endif
                else if(pi/8. < dir_mat(i,j,1) .and. dir_mat(i,j,1) <= pi/4. + pi/8.) then !case(45), northwest- southeast
                    if((j-1 < 1 .and. i+1 <= ldim(1))  .or. (i-1 < 1 .and. j+1 <= ldim(2))) then
                        if(mat_in(i,j,1) < mat_in(i+1,j+1,1)) temp_mat(i,j,1) = 0.        !just southeast
                    elseif((i+1 > ldim(1) .and. j-1 > 0) .or. (j+1 > ldim(2) .and. i-1 > 0)) then
                        if(mat_in(i,j,1) < mat_in(i-1,j-1,1)) temp_mat(i,j,1) = 0.    !just northwest
                    else if(i-1 > 0 .and. j+1 <= ldim(2) .and. i+1 <= ldim(1) .and. j-1 > 0) then
                        if(mat_in(i,j,1) < mat_in(i-1,j-1,1) .or. mat_in(i,j,1) < mat_in(i+1,j+1,1)) temp_mat(i,j,1) = 0. !northwest - southeast
                        !In the angles I decide not to to do anything
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
    subroutine double_thresh( mat_in, thresh )
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

    ! Canny routine for edge detection, in its classical version.
    ! It requires a double threshold.
    subroutine canny_edge( img_in, thresh, lp )
        class(image), intent(inout) :: img_in    ! input and output image
        real,            intent(in) :: thresh(2) ! low and high thresholds
        real, optional,  intent(in) :: lp(1)     ! lp filtering
        type(binimage) :: Gr !gradient image
        integer        :: ldim(3), r(3), k, i, j, nsz                   !just for implementation
        real,    allocatable :: dir_mat(:,:,:), Dc(:,:,:), Dr(:,:,:), grad(:,:,:)   !derivates, gradient
        integer :: neigh_8_inds(3,8) ! 8-neighborhoods of a pixel
        real    :: smpd, llp(1)
        ldim = img_in%get_ldim()
        smpd = img_in%get_smpd()
        allocate(grad(ldim(1),ldim(2),ldim(3)), Dc(ldim(1),ldim(2),ldim(3)), Dr(ldim(1),ldim(2),ldim(3)), &
        & source = 0.)
        if(ldim(3) /= 1) THROW_HARD("The image has to be 2D!; canny_edge")
        call Gr%new(ldim,smpd)
        !STEP 1: SMOOTHING
        !Apply a lp to reduce noise (traditionally it would be Gaussian Filter)
        llp = 8.
        if(present(lp)) llp = lp
        call img_in%bp(0.,llp(1))   !threshold selected in order to facilitate particle picking
        call img_in%ifft()
        !STEP 2: FINDING GRADIENTS
        allocate(dir_mat(ldim(1),ldim(2),ldim(3)), source = 0.)
        call img_in%calc_gradient(grad, Dc, Dr)
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
        call Gr%set_rmat(grad,.false.)
        do i = 1, ldim(1) !now grad is in {0,0.5,1}
            do j = 1, ldim(2)
                if( abs( Gr%get([i,j,1]) - 0.5 ) < TINY ) then
                    call Gr%set([i,j,1],0.)                            ! suppress this edge, later I might re-build it
                    call neigh_8(ldim, [i,j,1], neigh_8_inds, nsz)
                    do k = 1, nsz
                        r = neigh_8_inds(1:3,k)                        ! indexes of 8 neighbourhoods
                        if(is_equal( Gr%get([r(1),r(2),1]), 1.) ) then ! one of the 8 neigh is a strong edge
                            call Gr%set([i,j,1], 1.)                   ! re-build edge
                            exit
                        endif
                    enddo
                endif
            enddo
        enddo
        call img_in%copy(Gr)
        deallocate(grad)
    end subroutine canny_edge

    ! line Hough transform for line detection
    subroutine hough_line( img_in, threshold )
        class(image),      intent(inout) :: img_in
        integer, optional, intent(in)    :: threshold
        type(image) :: img_edge
        integer :: thresh, ldim(3), r_size, t_size, i, x, y, t, r, curr_rho_r, draw
        real    :: diagonal, curr_rho
        real,    allocatable :: emat(:,:,:), rho(:), theta(:), sins(:), coss(:)
        integer, allocatable :: line_pos(:,:), accumulator(:,:)
        real,    parameter :: theta_step=PI/180.
        integer, parameter :: tthresh=10, rho_step=1
        ldim     = img_in%get_ldim()
        if( ldim(3) /= 1 ) THROW_HARD('not yet implemented for 3D; hough_line')
        if( .not.present(threshold) )then 
            thresh = threshold
        else
            thresh = tthresh
        endif 
        ! pre-processing
        call canny(img_in, img_edge)
        emat     = img_edge%get_rmat()
        diagonal = ceiling(sqrt(real(ldim(1))**2. + real(ldim(2))**2.))
        r_size   = 2*diagonal
        t_size   = nint(PI/theta_step) + 1
        allocate(rho(r_size), theta(t_size))
        rho      = [(-diagonal + (i - 1)*rho_step, i = 1, r_size)]
        theta    = [(-PI + (i - 1)*theta_step, i = 1, t_size)]
        allocate(sins(t_size), coss(t_size), source=0.)
        allocate(accumulator(r_size, t_size), source=0)
        sins(:)  = sin(theta(:))
        coss(:)  = cos(theta(:))
        do x = 1, ldim(1)
            do y = 1, ldim(2)
                if( emat(x, y, 1) > 0. )then
                    do t = 1, t_size
                        curr_rho                  = x*coss(t) + y*sins(t)
                        curr_rho_r                = minloc(abs(rho - curr_rho), 1)
                        accumulator(curr_rho_r,t) = accumulator(curr_rho_r,t) + 1                
                    enddo
                endif 
            enddo 
        enddo  
        allocate(line_pos(ldim(1), ldim(2)), source=0)
        ! finding local maxima
        do t = 1, t_size
            do r = 1, r_size
                if( accumulator(r, t) > thresh .and. theta(t) > PI/2. - 0.01 .and. theta(t) < PI/2. + 0.01 )then 
                    do draw = 0, accumulator(r,t)
                        line_pos(1, nint(rho(r)*sins(t))) = accumulator(r,t)
                    enddo 
                endif 
            enddo 
        enddo 
        call img_edge%kill()
    end subroutine hough_line

    ! otsu binarization for images, based on the implementation
    ! of otsu algo for 1D vectors
    subroutine otsu_img( img, thresh, mskrad, positive, tight, tighter )
        class(image),      intent(inout) :: img
        real,    optional, intent(inout) :: thresh
        real,    optional, intent(in)    :: mskrad
        logical, optional, intent(in)    :: positive ! consider just the positive range (specific for nanoparticles)
        logical, optional, intent(in)    :: tight    ! for generating a tight mask in automasking
        logical, optional, intent(in)    :: tighter  ! for generating a tighter mask
        real,    pointer     :: rmat(:,:,:)
        real,    allocatable :: x(:)
        logical, allocatable :: lmsk(:,:,:)
        type(image) :: tmpimg
        integer     :: ldim(3)
        logical     :: ppositive, ttight, ttighter
        real        :: selected_t, outlier_t, t1, t2
        ldim = img%get_ldim()
        if( present(mskrad) )then
            call tmpimg%disc(ldim, img%get_smpd(), mskrad, lmsk)
            call tmpimg%kill
        endif
        ppositive = .false.
        if( present(positive) ) ppositive = positive
        ttight = .false.
        if( present(tight) ) ttight = tight
        ttighter = .false.
        if( present(tighter) ) ttighter = tighter
        call img%get_rmat_ptr(rmat)
        if( ppositive ) then
            if( allocated(lmsk) )then
                x = pack(rmat(1:ldim(1),1:ldim(2),1:ldim(3)), rmat(1:ldim(1),1:ldim(2),1:ldim(3)) > 0. .and. lmsk)
            else
                x = pack(rmat(1:ldim(1),1:ldim(2),1:ldim(3)), rmat(1:ldim(1),1:ldim(2),1:ldim(3)) > 0. )
            endif
        else
            if( allocated(lmsk) )then
                x = pack(rmat(1:ldim(1),1:ldim(2),1:ldim(3)), lmsk)
            else
                x = pack(rmat(1:ldim(1),1:ldim(2),1:ldim(3)), .true.)
            endif
        endif
        call otsu(size(x), x, selected_t)
        if(present(thresh)) thresh= selected_t
        t1 = selected_t
        if( ttight )then
            x = pack(x, x > selected_t)
            call otsu(size(x), x, selected_t)
            if(present(thresh)) thresh = selected_t
            t2 = selected_t
            if( count(rmat(1:ldim(1),1:ldim(2),1:ldim(3)) >= selected_t) < ldim(1) )then
                selected_t = t1
                if(present(thresh)) thresh = selected_t
            endif
        endif
        if( ttighter )then
            x = pack(x, x > selected_t)
            call otsu(size(x), x, selected_t)
            if(present(thresh)) thresh = selected_t
            t2 = selected_t
            x = pack(x, x > selected_t)
            call otsu(size(x), x, selected_t)
            if(present(thresh)) thresh = selected_t
            if( count(rmat(1:ldim(1),1:ldim(2),1:ldim(3)) >= selected_t) < ldim(1) )then
                selected_t = t2
                if(present(thresh)) thresh = selected_t
            endif
        endif
        deallocate(x)
        where(rmat(1:ldim(1),1:ldim(2),1:ldim(3))>=selected_t)
            rmat(1:ldim(1),1:ldim(2),1:ldim(3)) = 1.
        elsewhere
            rmat(1:ldim(1),1:ldim(2),1:ldim(3)) = 0.
        endwhere
        if( allocated(lmsk) ) deallocate(lmsk)
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
    subroutine otsu_robust_fast( img, is2D, noneg, thresh )
        class(image), intent(inout) :: img
        logical,      intent(in)    :: is2D    ! is it a 2D image
        logical,      intent(in)    :: noneg   ! is it a nanoparticle
        real,         intent(inout) :: thresh(3)
        type(image)       :: img_copy
        type(image)       :: img_avg, img_med
        real, allocatable :: neigh_8_pixs(:)
        integer :: ldim(3)
        integer :: i,j,k,nsz
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
        ! Fetch pointers
        call img%get_rmat_ptr(rmat_t)
        call img_copy%get_rmat_ptr(rmat)
        call img_avg%get_rmat_ptr(rmat_avg)
        call img_med%get_rmat_ptr(rmat_med)
        ! Generate avg and median volumes
        ! 1) img_avg: image where in each pixel is contained
        ! the average value of the gray levels in the neighbours of the pixel.
        ! 2) img_med: image where in each pixel is contained
        ! the median value of the gray levels in the neighbours of the pixel.
        if(ldim(3) == 1) then
            allocate(neigh_8_pixs(9))
            !$omp parallel do collapse(2) default(shared) private(i,j,neigh_8_pixs,nsz)&
            !$omp schedule(static) proc_bind(close)
            do i = 1, ldim(1)
                do j = 1, ldim(2)
                    call neigh_8(ldim, rmat(:ldim(1),:ldim(2),:1), [i,j,1], neigh_8_pixs, nsz)
                    rmat_avg(i,j,1) = sum(neigh_8_pixs(1:nsz)) / real(nsz)
                    rmat_med(i,j,1) = median(neigh_8_pixs(1:nsz))
                enddo
            enddo
            !$omp end parallel do
        else
            allocate(neigh_8_pixs(27))
            !$omp parallel do collapse(3) default(shared) private(i,j,k,neigh_8_pixs,nsz)&
            !$omp schedule(static) proc_bind(close)
            do i = 1, ldim(1)
                do j = 1, ldim(2)
                    do k = 1, ldim(3)
                        call neigh_8_3D(ldim, rmat(:ldim(1),:ldim(2),:ldim(3)), [i,j,k], neigh_8_pixs, nsz)
                        rmat_avg(i,j,k) = sum(neigh_8_pixs(1:nsz)) / real(nsz)
                        rmat_med(i,j,k) = median(neigh_8_pixs(1:nsz))
                    enddo
                enddo
            enddo
            !$omp end parallel do
        endif
        deallocate(neigh_8_pixs)
        ! Apply otsu1D to each img
        if( noneg .eqv. .true. )then
            call otsu_img(img_copy, thresh(1), positive = .true.)
            call otsu_img(img_avg,  thresh(2), positive = .true.)
            call otsu_img(img_med,  thresh(3), positive = .true.)
        else
            call otsu_img(img_copy, thresh(1), positive = .false.)
            call otsu_img(img_avg,  thresh(2), positive = .false.)
            call otsu_img(img_med,  thresh(3), positive = .false.)
        endif
        !Find most frequent value in the binary version of each voxel
        !$omp parallel do collapse(3) default(shared) private(i,j,k,count_back)&
        !$omp schedule(static) proc_bind(close)
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
        !$omp end parallel do
        ! kill
        call img_copy%kill
        call img_avg%kill
        call img_med%kill
    end subroutine otsu_robust_fast

    subroutine sauvola( img, winsz, img_sdevs, bias )
        class(image),   intent(inout) :: img
        integer,        intent(in)    :: winsz
        class(image),   intent(inout) :: img_sdevs ! local standard deviations
        real, optional, intent(in)    :: bias
        real, allocatable :: rmat(:,:,:), sdevs(:,:,:), avgs(:,:,:)
        integer :: ldim(3), ir(2), jr(2), isz, jsz, npix, i, j
        real    :: bbias, smpd, sdev_max, t
        bbias = 0.34 ! [0.2,0.5]
        if( present(bias) ) bbias = bias
        ldim = img%get_ldim()
        smpd = img%get_smpd()
        rmat = img%get_rmat()
        ! calculate the local standard deviations
        allocate(avgs(ldim(1),ldim(2),ldim(3)), sdevs(ldim(1),ldim(2),ldim(3)), source=0.)
        if( ldim(3) == 1 )then ! 2d
            !$omp parallel do default(shared) private(i,ir,isz,j,jr,jsz,npix) schedule(static) proc_bind(close)
            do i = 1, ldim(1)
                ir(1)            = max(1,       i - winsz)
                ir(2)            = min(ldim(1), i + winsz)
                isz              = ir(2) - ir(1) + 1
                do j = 1, ldim(2)
                    jr(1)        = max(1,       j - winsz)
                    jr(2)        = min(ldim(2), j + winsz)
                    jsz          = jr(2) - jr(1) + 1
                    npix         = isz * jsz
                    avgs(i,j,1)  = sum(rmat(ir(1):ir(2),jr(1):jr(2),1)) / real(npix)
                    sdevs(i,j,1) = sqrt(sum((rmat(ir(1):ir(2),jr(1):jr(2),1) - avgs(i,j,1))**2.0) / real(npix - 1))
                end do
            end do
            !$omp end parallel do
        else ! 3d
            THROW_HARD('not yet implemented for 3d')
        endif
        call img_sdevs%new(ldim, smpd)
        call img_sdevs%set_rmat(sdevs, .false.)
        sdev_max = maxval(sdevs)
        ! do the thresholding
        if( ldim(3) == 1 )then ! 2d
            !$omp parallel do default(shared) private(i,j,t) schedule(static) proc_bind(close)
            do i = 1, ldim(1)
                do j = 1, ldim(2)
                    t = avgs(i,j,1) * (1.0 + bbias * (sdevs(i,j,1) / sdev_max - 1.0))
                    if( rmat(i,j,1) >= t )then
                        rmat(i,j,1) = 1.0
                    else
                        rmat(i,j,1) = 0.
                    endif
                end do
            end do
            !$omp end parallel do
        endif
        call img%set_rmat(rmat, .false.)
    end subroutine sauvola

end module simple_segmentation
