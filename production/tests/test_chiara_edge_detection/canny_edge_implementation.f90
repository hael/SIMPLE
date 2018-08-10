module canny_no_thresh_mod
  include 'simple_lib.f08'
  use simple_image, only : image
  use simple_ctf,  only : ctf
  implicit none

  public :: sobel, canny_edge, canny, sobel_nothresh, build_ctf
  private

  contains
   !To build a ctf and visualize it in the image img_ctf
    function build_ctf(smpd, kv, cs, fraca, dfx, dfy, angast) result(img_ctf)
        real,        intent(in) :: smpd, kv, cs, fraca, dfx, dfy, angast
        type(image)             :: ctf_image, img4viz, img_ctf
        type(ctf)               :: tfun
        integer, parameter      :: box=256
        call ctf_image%new([box,box,1], smpd)
        call img4viz%new  ([box,box,1], smpd)
        tfun = ctf(smpd, kv, cs, fraca)
        call tfun%ctf2img(ctf_image, dfx, dfy, angast)
        call ctf_image%ft2img('real', img4viz)
        img_ctf = img4viz
    end function build_ctf

 !To sort a vector (low to high values)
  recursive subroutine QsortC(A)
        real, intent(inout) :: A(:)
        integer              :: iq
        if(size(A) > 1) then
          call Partition(A, iq)
          call QsortC(A(:iq-1))
          call QsortC(A(iq:))
        endif
  end subroutine QsortC

!To help sorting a vector
  subroutine Partition(A, marker)
      real,    intent(inout) :: A(:)
      integer, intent(out)   :: marker
      integer :: i, j
      real    :: temp, x
      x = A(1)
      i = 0
      j = size(A) + 1
      do
        j = j-1
        do
          if (A(j) <= x) exit
          j = j-1
        end do
        i = i+1
        do
          if (A(i) >= x) exit
          i = i+1
        end do
        if (i < j) then
            temp = A(i)  ! exchange A(i) and A(j)
            A(i) = A(j)
            A(j) = temp
        elseif (i == j) then
            marker = i+1
            return
        else
           marker = i
           return
       endif
     end do
   end subroutine Partition

!To calculate the median value of the intensities of entries of the matrix
    function median_img(mat_in) result(m)
        real, intent(in)  :: mat_in(:,:,:)
        real, allocatable :: vect2(:), vect(:)
        real              :: m
        integer           :: ldim(3), dim
        ldim = shape(mat_in)
        allocate(vect(ldim(1)*ldim(2)), source = reshape(mat_in, [ldim(1)*ldim(2)]))  !serialize
        call elim_dup(vect,vect2)
        deallocate(vect)
        call QsortC(vect2)
        dim = size(vect2)
        if (mod(dim, 2) == 0) then !even
            m = (vect2(dim/2)+vect2(dim/2+1))/2.
        else
            m = vect2((dim+1)/2)   !odd
        endif
        deallocate(vect2)
    end function median_img

    !Returns 8-neighborhoods of the pixel px in the matrix mat. This
    !function is specific for Canny edge detection. If Canny will
    !be inserted in simple_image, it would be useful to merge this
    !function with calc_neigh_8 present in simple_image.
    function calc_neigh_8_canny(mat, px) result(neigh_8)
      real,    intent(in)  :: mat(:,:,:)
      integer, intent(in)  :: px(3)
      integer, allocatable :: neigh_8(:,:,:)
      integer              :: i, j, ldim(3)
      ldim = shape(mat)
      if(px(3) /= 1) then
          print *, "The matrix has to be 2D!"
          stop
      endif
      i = px(1)
      j = px(2)            !Assumes to have a 2-dim matrix
      if ( i-1 < 1 .and. j-1 < 1 ) then
        allocate(neigh_8(3,3,1), source = 0)
        neigh_8(1:3,1,1) = [i+1,j,1]
        neigh_8(1:3,2,1) = [i+1,j+1,1]
        neigh_8(1:3,3,1) = [i,j+1,1]
      else if (j+1 > ldim(2) .and. i+1 > ldim(1)) then
        allocate(neigh_8(3,3,1), source = 0)
        neigh_8(1:3,1,1) = [i-1,j,1]
        neigh_8(1:3,2,1) = [i-1,j-1,1]
        neigh_8(1:3,3,1) = [i,j-1,1]
      else if (j-1 < 1  .and. i+1 >ldim(1)) then
        allocate(neigh_8(3,3,1), source = 0)
        neigh_8(1:3,1,1) = [i-1,j,1]
        neigh_8(1:3,2,1) = [i-1,j+1,1]
        neigh_8(1:3,3,1) = [i,j+1,1]
      else if (j+1 > ldim(2) .and. i-1 < 1) then
        allocate(neigh_8(3,3,1), source = 0)
        neigh_8(1:3,1,1) = [i,j-1,1]
        neigh_8(1:3,2,1) = [i+1,j-1,1]
        neigh_8(1:3,3,1) = [i+1,j,1]
      else if( j-1 < 1 ) then
        allocate(neigh_8(3,5,1), source = 0)
        neigh_8(1:3,1,1) = [i-1,j,1]
        neigh_8(1:3,2,1) = [i-1,j+1,1]
        neigh_8(1:3,3,1) = [i,j+1,1]
        neigh_8(1:3,4,1) = [i+1,j+1,1]
        neigh_8(1:3,5,1) = [i+1,j,1]
      else if ( j+1 > ldim(2) ) then
        allocate(neigh_8(3,5,1), source = 0)
        neigh_8(1:3,1,1) = [i-1,j,1]
        neigh_8(1:3,2,1) = [i-1,j-1,1]
        neigh_8(1:3,3,1) = [i,j-1,1]
        neigh_8(1:3,4,1) = [i+1,j-1,1]
        neigh_8(1:3,5,1) = [i+1,j,1]
      else if ( i-1 < 1 ) then
        allocate(neigh_8(3,5,1), source = 0)
        neigh_8(1:3,1,1) = [i,j-1,1]
        neigh_8(1:3,2,1) = [i+1,j-1,1]
        neigh_8(1:3,3,1) = [i+1,j,1]
        neigh_8(1:3,4,1) = [i+1,j+1,1]
        neigh_8(1:3,5,1) = [i,j+1,1]
      else if ( i+1 > ldim(1) ) then
        allocate(neigh_8(3,5,1), source = 0)
        neigh_8(1:3,1,1) = [i,j+1,1]
        neigh_8(1:3,2,1) = [i-1,j+1,1]
        neigh_8(1:3,3,1) = [i-1,j,1]
        neigh_8(1:3,4,1) = [i-1,j-1,1]
        neigh_8(1:3,5,1) = [i,j-1,1]
      else
        allocate(neigh_8(3,8,1), source = 0)
        neigh_8(1:3,1,1) = [i-1,j-1,1]
        neigh_8(1:3,2,1) = [i,j-1,1]
        neigh_8(1:3,3,1) = [i+1,j-1,1]
        neigh_8(1:3,4,1) = [i+1,j,1]
        neigh_8(1:3,5,1) = [i+1,j+1,1]
        neigh_8(1:3,6,1) = [i,j+1,1]
        neigh_8(1:3,7,1) = [i-1,j+1,1]
        neigh_8(1:3,8,1) = [i-1,j,1]
      endif
    end function calc_neigh_8_canny

!This function returns a matrix in which the first columns contain colum derivates,
!the last columns contain row derivates of the input image
    subroutine calc_gradient(self, grad, Dc, Dr)
      type(image),                 intent(inout) :: self
      real, allocatable,           intent(out)   :: grad(:,:,:)
      real, allocatable, optional, intent(out)   :: Dc(:,:,:), Dr(:,:,:)
      type(image)        :: img_p                         !padded image
      real, allocatable  :: wc(:,:,:), wr(:,:,:)          !row and column Sobel masks
      integer, parameter :: L = 3                         !dimension of the masks
      integer            :: ldim(3)                       !dimension of the image
      integer            :: i,j,m,n                       !loop indeces
      real, allocatable  :: Ddc(:,:,:),Ddr(:,:,:)         !column and row derivates
      real, allocatable  :: mat_in(:,:,:), mat_p(:,:,:)   !images, just the matrix
      ldim = self%get_ldim()
      if(ldim(3) /= 1) then
        print *, "The image has to be 2D!"
        stop
      endif
      allocate( Ddc(ldim(1),ldim(2),1), Ddr(ldim(1),ldim(2),1), grad(ldim(1),ldim(2),1), &
               & wc(-(L-1)/2:(L-1)/2,-(L-1)/2:(L-1)/2,1),wr(-(L-1)/2:(L-1)/2,-(L-1)/2:(L-1)/2,1), source = 0.)
      wc = (1./8.)*reshape([-1,0,1,-2,0,2,-1,0,1],[3,3,1])      !Sobel masks
      wr = (1./8.)*reshape([-1,-2,-1,0,0,0,1,2,1],[3,3,1])
      mat_in = self%get_rmat()
      call img_p%new([ldim(1)+L-1,ldim(2)+L-1,1],1.)
      call self%pad(img_p)                                      !padding
      mat_p = img_p%get_rmat()
      do i = 1, ldim(1)
          do j = 1, ldim(2)
              do m = -(L-1)/2,(L-1)/2
                  do n = -(L-1)/2,(L-1)/2
                      Ddc(i,j,1) = Ddc(i,j,1)+mat_p(i+m+1,j+n+1,1)*wc(m,n,1)
                      Ddr(i,j,1) = Ddr(i,j,1)+mat_p(i+m+1,j+n+1,1)*wr(m,n,1)
                  end do
              end do
          end do
      end do
      deallocate(wc,wr)
      grad = sqrt(Ddc**2+Ddr**2)
      if(present(Dc)) allocate(Dc(ldim(1),ldim(2),1), source = Ddc)
      if(present(Dr)) allocate(Dr(ldim(1),ldim(2),1), source = Ddr)
      deallocate(Ddc,Ddr)
    end subroutine calc_gradient

    !NON_MAX_SUPPRESSION
    !using the estimates of the Gx and Gy image gradients and the edge direction angle
    !determines whether the magnitude of the gradient assumes a local  maximum in the gradient direction
    subroutine non_max_supp(mat_in,dir_mat)
      real, intent(inout) :: mat_in(:,:,:)
      real,  intent(in)   :: dir_mat(:,:,:)
      real, allocatable   :: rmat(:,:,:), temp_mat(:,:,:)
      integer             :: ldim(3), i, j
      ldim = shape(mat_in)
      if(ldim(3) /= 1) then
        print *, "The matrix has to be 2D!"
        stop
      endif
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

    !Edges stronger than thresh(2) are mantained, edges weaker than thresh(1) are discarded, in between edges are marked as "weak edges"
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

      !Edge detection, Sobel algorithm
      subroutine sobel(img_in,img_out,thresh)
          type(image), intent(inout) :: img_in,img_out         !image input and output
          real,        intent(in)    :: thresh                 !threshold for Sobel algorithm
          integer            :: ldim(3)                        !dimension of the image
          integer            :: i,j                            !loop indeces
          real,  allocatable :: grad(:,:,:)                    !matrices
          ldim = img_in%get_ldim()
          call img_out%new(ldim,1.)     !reset if not empty
          call calc_gradient(img_in, grad)
          do i=1, ldim(1)               !where( grad > thresh ) img_out = 1.
              do j=1, ldim(2)
                if(grad(i,j,1) > thresh) call img_out%set([i,j,1], 1.)
              enddo
          enddo
          deallocate(grad)
      end subroutine sobel

     !Sobel no thresh
      subroutine sobel_nothresh(img_in,img_out, sigma)
          type(image),    intent(inout) :: img_in,img_out         !image input and output
          real, optional, intent(in)    :: sigma
          real               :: threshold, ssigma, m               !threshold for Sobel algorithm
          real,  allocatable :: grad(:,:,:)                    !matrices
          type(image)        :: Gr                             !Gradient image
          real, allocatable  :: rmat(:,:,:)
          ssigma = 0.33              !empirical experiments (papers advices)
          if(present(sigma)) ssigma = sigma
          call img_in%scale_img([0.,255.])
          call img_out%new(img_in%get_ldim(),1.)     !reset if not empty
          img_out = 0.
          rmat = img_out%get_rmat()
          call calc_gradient(img_in, grad)
          call Gr%new(img_in%get_ldim(),img_in%get_smpd())
          call Gr%set_rmat(grad)
          m = Gr%median_value()
          threshold = ( max(0.,(1-ssigma)*m) + min(255.,(1+ssigma)*m) )/2
          where( grad > threshold ) rmat = 1.
          call img_out%set_rmat(rmat)
          deallocate(grad, rmat)
      end subroutine sobel_nothresh

      !Canny routine for edge detection, it requires a double threshold
        subroutine canny_edge(img_in,img_out, thresh)
            type(image), intent(inout) :: img_in, img_out                                             !input and output images
            real,        intent(in)    :: thresh(2)                                                   !low and high thresholds
            integer              :: ldim(3), s(3), r(3), k, i, j                                      !just for implementation
            real, allocatable    :: dir_mat(:,:,:), Dc(:,:,:), Dr(:,:,:), grad(:,:,:)                 !derivates, gradient
            integer, allocatable :: neigh_8(:,:,:)                                                    !8-neighborhoods of a pixel
            ldim = img_in%get_ldim()
            if(ldim(3) /= 1) then
              print *, "The image has to be 2D!"
              stop
            endif
            call img_out%new(ldim,1.)  !reset if not empty
            !STEP 1: SMOOTHING
            !Apply a Gaussian filter to reduce noise

            !STEP 2: FINDING GRADIENTS
            allocate(dir_mat(ldim(1),ldim(2),1), source = 0.)
            call calc_gradient(img_in, grad, Dc, Dr)
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
            do i = 1, ldim(1) !now grad is in {0,0.5,1}
                do j = 1, ldim(2)
                    if(grad(i,j,1) == 0.5) then
                        grad(i,j,1) = 0.                         !suppress this edge, later I might re-build it
                        neigh_8 = calc_neigh_8_canny(grad,[i,j,1])
                        s = shape(neigh_8)
                        do k = 1, s(2)
                            r = int(neigh_8(1:3,k,1))           !indexes of 8 neighbourhoods
                            if(grad(r(1),r(2),1) == 1.) then    !one of the 8 neigh is a strong edge
                               grad(i,j,1) = 1.                 !re-build edge
                               exit
                            endif
                        enddo
                    endif
                enddo
            enddo
            call img_out%set_rmat(grad)
            deallocate(grad)
        end subroutine canny_edge

   !Canny edge detection with no need of a threshold
    subroutine canny(img_in, img_out, sigma)
            type(image), intent(inout) :: img_in, img_out
            real                       :: thresh(2), ssigma, m
            real, allocatable          :: grad(:,:,:)
            real, optional             :: sigma     !It behaves like a parameter to decide "how much" has to be detected
            integer                    :: ldim(3)
            ssigma = 0.33              !empirical experiments (papers advices)
            if(present(sigma)) ssigma = sigma
            call img_in%scale_img([0.,255.])
            ldim = img_in%get_ldim()
            call calc_gradient(img_in, grad)
            m = median_img(grad) !I used the gradient, the source talks about the image itself
            !https://www.pyimagesearch.com/2015/04/06/zero-parameter-automatic-canny-edge-detection-with-python-and-opencv/
            thresh(1) = max(0.   ,(1-ssigma)*m) !lower
            thresh(2) = min(255., (1+ssigma)*m) !upper
            call canny_edge(img_in, img_out, thresh)
    end subroutine canny
end module canny_no_thresh_mod

program canny_no_thresh
use canny_no_thresh_mod
use simple_image, only : image
use simple_jpg
implicit none
type(image)  :: img_in, img_out, img_in2, img_in1, img_in05

call img_in%new([128,128,1],1.)
call img_in2%new([128,128,1],1.)
call img_in1%new([128,128,1],1.)
call img_in05%new([128,128,1],1.)
call img_in%read('/home/lenovoc30/Desktop/Edge Detection/one_projection.mrc')
call img_in2%read('/home/lenovoc30/Desktop/Edge Detection/one_projection2.mrc')
call img_in1%read('/home/lenovoc30/Desktop/Edge Detection/one_projection1.mrc')
call img_in05%read('/home/lenovoc30/Desktop/Edge Detection/one_projection05.mrc')
!Let's use this image and apply three different metods for edge detection on it:
!1) Sobel, threshold:
!2) Canny, threshold:
!3) Canny, threshold: NO
!NOISELESS IMAGE
    call sobel(img_in,img_out, 0.5)           !1)
    call img_out%write('Sobel.mrc')
    call canny_edge(img_in,img_out, [0.1,0.5]) !2)
    call img_out%write('Canny_thresh.mrc')
    call canny(img_in,img_out, 0.9)                 !3)
    call img_out%write('Canny_NO_thresh.mrc')
!IMAGE WITH SNR = 2
    call sobel(img_in2,img_out, 0.99)           !1)
    call img_out%write('Sobel2.mrc')
    call canny_edge(img_in2,img_out, [0.5,0.95]) !2)
    call img_out%write('Canny_thresh2.mrc')
    call canny(img_in2,img_out, 0.9)                 !3)
    call img_out%write('Canny_NO_thresh2.mrc')
!IMAGE WITH SNR = 1
    call sobel(img_in1,img_out, 0.99)           !1)
    call img_out%write('Sobel1.mrc')
    call canny_edge(img_in1,img_out, [0.5,0.95]) !2)
    call img_out%write('Canny_thresh1.mrc')
    call canny(img_in1,img_out, 0.9)                 !3)
    call img_out%write('Canny_NO_thresh1.mrc')
!IMAGE WITH SNR = 0.5
    call sobel(img_in05,img_out, 0.99)           !1)
    call img_out%write('Sobel05.mrc')
    call canny_edge(img_in05,img_out, [0.5,0.95]) !2)
    call img_out%write('Canny_thresh05.mrc')
    call canny(img_in05,img_out, 0.9)                 !3)
    call img_out%write('Canny_NO_thresh05.mrc')
end program canny_no_thresh
!If it doesn't work as you would like you can also check "Auto Canny" in GIT directory
