module try_mod
  include 'simple_lib.f08'
  use simple_image
  use simple_ctf,  only : ctf
  implicit none
  contains
    function calc_neigh_8(mat, px) result(neigh_8)
      real, intent(in) :: mat(:,:,:)
      integer, intent(in) :: px(3)
      integer, allocatable :: neigh_8(:,:,:)
      integer :: i, j, ldim(3)
      ldim = shape(mat)
      i = px(1)
      j = px(2)
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
    end function calc_neigh_8

    subroutine print_mat(matrix)
        real, intent(in) :: matrix(:,:,:)
        integer          :: j,s(3)
        s = shape(matrix)
        do j = 1, s(2)
            print *, matrix(:,j,1)
        enddo
    end subroutine print_mat

    function calc_gradient(self) result(Dcr)
      type(image), intent(inout)          :: self
      type(image)                         :: img_p , img_out!, Gr            !padded image
      real, dimension(:,:,:), allocatable :: wc, wr            !row and column masks
      integer, parameter                  :: L = 3             !dimension of the mask
      integer, dimension(3)               :: ldim              !dimension of the image
      integer                             :: i,j,m,n         !loop indexes
      real, dimension(:,:,:), allocatable :: Dc,Dr, Dcr        !Column and row derivates and both of them
      real, dimension(:,:,:), allocatable :: mat_in, mat_p, mat_out     !images, just the matrix
      !real, dimension(:,:,:), allocatable :: grad_mat          !gradient and padded image, just the matrix
      ldim = self%get_ldim()
      allocate( Dc(ldim(1),ldim(2),1), Dcr(2*ldim(1),ldim(2),1), mat_out(ldim(1),ldim(2),1), &
               & Dr(ldim(1),ldim(2),1),wc(-(L-1)/2:(L-1)/2,-(L-1)/2:(L-1)/2,1),wr(-(L-1)/2:(L-1)/2,-(L-1)/2:(L-1)/2,1), source = 0.)
      wc = (1./8.)*reshape([-1,0,1,-2,0,2,-1,0,1],[3,3,1])      !Sobel masks
      wr = (1./8.)*reshape([-1,-2,-1,0,0,0,1,2,1],[3,3,1])
      mat_in = self%get_rmat()
      call img_p%new([ldim(1)+L-1,ldim(2)+L-1,1],1.)
      call self%pad(img_p)                                     !padding
      mat_p = img_p%get_rmat()
      do i = 1, ldim(1)
          do j = 1, ldim(2)
              do m = -(L-1)/2,(L-1)/2
                  do n = -(L-1)/2,(L-1)/2
                      Dc(i,j,1) = Dc(i,j,1)+mat_p(i+m+1,j+n+1,1)*wc(m,n,1)
                      Dr(i,j,1) = Dr(i,j,1)+mat_p(i+m+1,j+n+1,1)*wr(m,n,1)
                  end do
              end do
          end do
      end do
      Dcr(1:ldim(1), 1:ldim(2), 1) = Dc(:,:,1)
      Dcr(ldim(1)+1:2*ldim(1),1:ldim(2),1) = Dr(:,:,1)
      !grad_mat= sqrt(Dc**2+Dr**2);  !gradient matrix
      !call Gr%new(ldim,1.)          !gradient image
      !call Gr%set_rmat(grad_mat)
    end function calc_gradient

    !NON_MAX_SUPPRESSION
    !using the estimates of the Gx and Gy image gradients and the edge direction angle
    !determines whether the magnitude of the gradient assumes a local  maximum in the gradient direction
    subroutine non_max_supp(img_in,dir_mat)
      type(image), intent(inout) :: img_in
      real,  intent(in)          :: dir_mat(:,:,:)
      real, allocatable          :: rmat(:,:,:), temp_mat(:,:,:)
      integer                    :: ldim(3), i, j
      ldim = img_in%get_ldim()
    !if the rounded edge direction angle is 0 degrees, checks the north and south directions
    !if the rounded edge direction angle is 45 degrees, checks the northwest and southeast directions
    !if the rounded edge direction angle is 90 degrees, checks the east and west directions
    allocate(rmat(ldim(1),ldim(2),1), temp_mat(ldim(1),ldim(2),1), source = 0.)
    rmat = img_in%get_rmat()
    do i = 1, ldim(1)
      do j = 1, ldim(2)
               if( pi/4. + pi/8. < dir_mat(i,j,1) .and. dir_mat(i,j,1) <= pi/2.) then         !case(90), north - south
                    if(j+1 > ldim(2)) then
                          if(rmat(i,j,1) > rmat(i,j-1,1)) temp_mat(i,j,1) = rmat(i,j,1)
                    else if(j-1 < 1) then
                          if(rmat(i,j,1) > rmat(i,j+1,1)) temp_mat(i,j,1) = rmat(i,j,1)
                    else
                          if(rmat(i,j,1) > rmat(i,j+1,1) .and. rmat(i,j,1) > rmat(i,j-1,1)) temp_mat(i,j,1) = rmat(i,j,1)
                    end if
            else if( 0. <= dir_mat(i,j,1) .and. dir_mat(i,j,1) <= pi/8.) then                !case(0), east - west
                  if(i+1 > ldim(1)) then
                        if(rmat(i,j,1) > rmat(i-1,j,1)) temp_mat(i,j,1) = rmat(i,j,1)
                  else if(i-1 < 1) then
                        if(rmat(i,j,1) > rmat(i+1,j,1)) temp_mat(i,j,1) = rmat(i,j,1)
                  else
                        if(rmat(i,j,1) > rmat(i+1,j,1) .and. rmat(i,j,1) > rmat(i-1,j,1)) temp_mat(i,j,1) = rmat(i,j,1)
                  endif
             else if(pi/8. < dir_mat(i,j,1) .and. dir_mat(i,j,1) <= pi/4. + pi/8.) then     !case(45), northwest - southeast
                 if( (i-1 < 1 .and. j-1 >0) .or.(j+1 > ldim(2) .and. i+1 <= ldim(1)) ) then
                      if(rmat(i,j,1) > rmat(i+1,j-1,1)) temp_mat(i,j,1) = rmat(i,j,1)   !just northwest
                 else if( (i-1 >0  .and. j-1 < 1) .or. ( i+1 > ldim(1) .and. j+1 <= ldim(2)) ) then
                      if(rmat(i,j,1) > rmat(i-1,j+1,1)) temp_mat(i,j,1) = rmat(i,j,1)  !just southeast
                 else
                      if(rmat(i,j,1) > rmat(i-1,j+1,1) .and. rmat(i,j,1) > rmat(i+1,j-1,1)) temp_mat(i,j,1) = rmat(i,j,1) !northwest - southeast
                !In the angles I decide not to to anything
                endif
            else !case default
             print *, "There is an error in the approximation of the gradient direction"
           end if
      enddo
    enddo
    call img_in%set_rmat(temp_mat)
    call img_in%vis
    end subroutine non_max_supp

    subroutine double_thresh(img_in, thresh)
      real, intent(in)           :: thresh(2)
      type(image), intent(inout) :: img_in
      real, allocatable          :: mat(:,:,:), mask(:,:,:)
      integer                    :: ldim(3), i, j
      ldim = img_in%get_ldim()
      allocate(mat(ldim(1),ldim(2),1),mask(ldim(1),ldim(2),1), source = 0.)
      mat = img_in%get_rmat()
      print *, "Maximum value of gradient bf thresh: ", maxval(mat)
      print *, "Minimum value of gradient bf thresh: ", minval(mat)
      do i = 1, ldim(1)
        do j = 1, ldim(2)
              if(mat(i,j,1) > thresh(2))                                 mask(i,j,1) = 1. !strong edges
              if(mat(i,j,1) < thresh(1))                                 mask(i,j,1) = 0. !not edges
              if(thresh(1) <= mat(i,j,1) .and. mat(i,j,1) <= thresh(2))  mask(i,j,1) = 0.5 !weak edges
        enddo
      enddo
      !where(mat > thresh(2))
      !  mask = 1.    !strong edges
      !elsewhere(thresh(1) <= mat .and. mat <= thresh(2))
      !  mask = 0.5   !weak edges
      !elsewhere(mat < thresh(1))
      !  mask = 0.    !not edges
      !endwhere
      !mat = mask !double thresholding
      call img_in%set_rmat(mask)     !img_in will be my gradient
    end subroutine double_thresh

      subroutine canny_edge(img_in,img_out, thresh)
          type(image), intent(inout) :: img_in, img_out
          type(image)                :: Gr
          real, intent(in)           :: thresh(2)
          integer                    :: ldim(3), s(3), k, i, j
          real, allocatable          :: d_mat(:,:,:), dir_mat(:,:,:), Dc(:,:,:), Dr(:,:,:), grad(:,:,:)
          integer, allocatable       :: neigh_8(:,:,:)
          integer                    :: r(3), cnt
          ldim = img_in%get_ldim()
          call img_out%new(ldim,1.)  !reset if not empty
          !STEP 1: SMOOTHING
          !call img_in%gauran(snr = 1.1, im_noise)   !Step 1, gaussian blurring
          !STEP 2: FINDING GRADIENTS
          allocate( d_mat(2*ldim(1),ldim(2),1), Dr(ldim(1),ldim(2),1), Dc(ldim(1),ldim(2),1), grad(ldim(1),ldim(2),1), &
                  & dir_mat(ldim(1),ldim(2),1), source = 0.)
          d_mat = calc_gradient(img_in)
          Dr(:,:,1) = d_mat(1:ldim(1), 1:ldim(2),1)
          Dc(:,:,1) = d_mat(ldim(1)+1:2*ldim(1), 1:ldim(2),1)
          do i = 1, ldim(1)
            do j = 1, ldim(2)
              if(Dc(i,j,1) /= 0.) then
                dir_mat(i,j,1) = atan(abs(Dr(i,j,1))/abs(Dc(i,j,1))) !Output of ATAN is in radians, not sure I have to put abs
              else
                dir_mat(i,j,1) = pi/2 !I chose this, let's see
              endif
            enddo
          enddo
          grad = sqrt(Dc**2+Dr**2)
          call Gr%new(ldim,1.)
          call Gr%set_rmat(grad)
          call Gr%vis
          !print *, "Maximum value: ", maxval(dir_mat)
          !print *, "Minimum value: ", minval(dir_mat)
          !print *, "pi/2: ", pi/2
          !STEP 3: NON-MAXIMUM SUPPRESSION
          call non_max_supp(Gr,dir_mat)
          !STEP 4: DOUBLE THRESHOLDING
          call double_thresh(Gr,thresh)
          !STEP 5: EDGE TRACKING BY HYSTERESIS
          cnt = 0
          do i = 1, ldim(1)
              do j = 1, ldim(2)
                  if(grad(i,j,1) == 0.5) then
                    cnt = cnt+1
                      grad(i,j,1) = 0.                        !suppress this edge, later I might re-build it
                      neigh_8 = calc_neigh_8(grad,[i,j,1])
                      s = shape(neigh_8)
                      do k = 1, s(2)
                          r = int(neigh_8(1:3,k,1))              !indexes of 8 neighbourhoods
                          if(grad(r(1),r(2),1) == 1.) then              !one of the 8 neigh is a strong edge
                             grad(i,j,1) = 1.                 !re-build edge
                             exit
                          endif
                      enddo
                  endif
              enddo
          enddo
          print *,"Worked?", cnt
          call img_out%set_rmat(grad)
      end subroutine canny_edge

      function build_ctf(smpd, kv, cs, fraca, dfx, dfy, angast) result(img_ctf)
          integer, parameter     :: box=256
          real                   :: smpd, kv, cs, fraca, dfx, dfy, angast
          type(image)            :: ctf_image, img4viz, img_ctf
          type(ctf)              :: tfun
          call ctf_image%new([box,box,1], smpd)
          call img4viz%new  ([box,box,1], smpd)
          tfun = ctf(smpd, kv, cs, fraca)
          call tfun%ctf2img(ctf_image, dfx, dfy, angast)
          call ctf_image%ft2img('real', img4viz)
          img_ctf = img4viz
      end function build_ctf

  end module try_mod

program try
use try_mod
implicit none
integer, allocatable :: neigh_8(:,:,:)
real, allocatable    :: rmat(:,:,:), matrix(:,:,:), dir_mat(:,:,:), Dc(:,:,:), Dr(:,:,:), d_mat(:,:,:), dir_mat_help(:,:,:)
integer              :: k, s(3), r(3), i, j
type(image)          :: img_in, img_out, img4viz
integer              :: box, ldim(3)
real                 :: thresh(2)
box    = 256
ldim   = [256,256,1]
thresh = [0.1,0.6]
!call img_in%new( ldim,1. )
!call img4viz%new( [ldim(1),ldim(2),1],1. )
img_in = build_ctf(smpd=1.0, kv=300., cs=2.7, fraca=0.1, dfx=0., dfy=0., angast=0.)
allocate(rmat(ldim(1),ldim(2),1), source = 0.)
!call canny_edge(img_in,img_out, thresh)
!rmat = img_out%get_rmat()
!print *, "Maximum value of result: ", maxval(rmat)
!print *, "Minimum value of result: ", minval(rmat)
!call img_out%vis
allocate(matrix(5,4,1), source = 0.)
matrix = reshape([2.,3.,5.,4.,6.,4.,5.,7.,6.,7.,5.,6.,4.,3.,2.,3.,4.,3.,1.,1.],[5,4,1])
print *, "Original matrix: "
call print_mat(matrix)
call img4viz%new( [5,4,1],1. )
call img4viz%set_rmat(matrix)
call img4viz%vis
allocate(Dc(5,4,1),Dr(5,4,1),dir_mat(5,4,1), d_mat(10,5,1),dir_mat_help(5,4,1), source = 0.)
d_mat = calc_gradient(img4viz)
Dr(:,:,1) = d_mat(1:5, 1:4,1)
Dc(:,:,1) = d_mat(6:10,1:4,1)
do i = 1, 5
  do j = 1, 4
    if(Dc(i,j,1) /= 0.) then
      dir_mat(i,j,1) = atan(abs(Dr(i,j,1))/abs(Dc(i,j,1))) !Output of ATAN is in radians, not sure I have to put abs
    else
      dir_mat(i,j,1) = pi/2. !I chose this, let's see
    endif
  enddo
enddo
print *, "Direction matrix radians: "
call print_mat(dir_mat)

!print *, "Direction matrix degrees: "
!do i = 1, 5
!  do j = 1, 4
!        if((pi/8. < dir_mat(i,j,1)) .and. (dir_mat(i,j,1) <= (pi/4. + pi/8.))) then     !case(45), northwest - southeast
!                dir_mat_help(i,j,1) = 45.
!
!        else if( ((pi/4. + pi/8.) < dir_mat(i,j,1)) .and. (dir_mat(i,j,1) <= pi/2.)) then         !case(90), north - south
!                dir_mat_help(i,j,1) = 90.
!
!        else if((0. <= dir_mat(i,j,1)) .and. (dir_mat(i,j,1) <= pi/8.)) then   ! .              !case(0), east - west
!                dir_mat_help(i,j,1) = 0.
!
!        else if ((dir_mat(i,j,1) > pi/2.) .or. (dir_mat(i,j,1) < 0.)) then                    !case default
!         print *, "There is an error in the approximation of the gradient direction"
!       end if
!  enddo
!enddo
!call print_mat(dir_mat_help)
!dir_mat = dir_mat_help
call non_max_supp(img4viz,dir_mat)
matrix = img4viz%get_rmat()
print *, "Non maximum suppression: "
call print_mat(matrix)

end program try
