module connected_components
  include 'simple_lib.f08'
  use simple_image, only : image
contains

  subroutine elimin_cc(img)
      type(image), intent(inout) :: img !image which contains connected components
      integer, allocatable :: sz(:,:), biggest_cc(:), biggest_val(:)
      real, allocatable    :: rmat(:,:,:)

      rmat = img%get_rmat()
      sz = size_connected_comps(img)
      biggest_cc  = maxloc(sz(2,:))!, dim = 2)
      biggest_val = real(sz(1, biggest_cc))
      ! print *, "Biggest cc",  biggest_cc
      ! print *, "Biggest val", biggest_val
      where(rmat /= biggest_val(1)) rmat = 0.
      call img%set_rmat(rmat)
      deallocate(rmat, sz,biggest_cc)
  end subroutine elimin_cc

  function size_connected_comps(img) result(sz)
      type(image), intent(in) :: img
      integer, allocatable :: sz(:,:)
      real,    allocatable :: rmat(:,:,:)
      integer :: n_cc

      rmat = img%get_rmat()
      allocate(sz(2,int(maxval(rmat))), source = 0) !In the first column it's stored the cc label, in the second one its size
      do n_cc = int(maxval(rmat)),1,-1
          sz(1, n_cc) = n_cc                  !I could avoid to save it
          sz(2, n_cc) = count(rmat == n_cc)
      enddo
  end function size_connected_comps

      subroutine elimin_isolated_points(bin_img)
        type(image), intent(inout) :: bin_img
        real, allocatable :: neigh_8(:), rmat(:,:,:)
        integer           :: ldim(3), i, j
        real              :: cn   !crossing number

        ldim = bin_img%get_ldim()
        rmat = bin_img%get_rmat()
        do i = 1, ldim(1)
           do j = 1, ldim(2)
              if(rmat(i,j,1) /= 0.) then         !Consider just white pixels
                  neigh_8 = calc_neigh_8(rmat,[i,j,1])
                  cn = crossing_number(neigh_8)
                  if(cn == 0.) rmat(i,j,1) = 0.  !Eliminate isolated points
                  deallocate(neigh_8)
              endif
           enddo
        enddo
        call bin_img%set_rmat(rmat)
        deallocate(rmat)
      end subroutine elimin_isolated_points

  !*****************************************************************************80
  !
  !! img_components assigns contiguous nonzero pixels to a common component.
  !
  !  Discussion:
  !
  !    On input, the A array contains values of 0 or 1.
  !
  !    The 0 pixels are to be ignored.  The 1 pixels are to be grouped
  !    into connected components.
  !
  !    The pixel A(I,J) is "connected" to the pixels A(I-1,J), A(I+1,J),
  !    A(I,J-1) and A(I,J+1), so most pixels have 4 neighbors.
  !
  !    (Another choice would be to assume that a pixel was connected
  !    to the other 8 pixels in the 3x3 block containing it.)
  !
  !    On output, COMPONENT_NUM reports the number of components of nonzero
  !    data, and the array C contains the component assignment for
  !    each nonzero pixel, and is 0 for zero pixels.
  !  Picture:
  !
  !    Input A:
  !
  !      0  2  0  0 17  0  3
  !      0  0  3  0  1  0  4
  !      1  0  4  8  8  0  7
  !      3  0  6 45  0  0  0
  !      3 17  0  5  9  2  5
  !
  !    Output:
  !
  !      COMPONENT_NUM = 4
  !
  !      C:
  !
  !      0  1  0  0  2  0  3
  !      0  0  2  0  2  0  3
  !      4  0  2  2  2  0  3
  !      4  0  2  2  0  0  0
  !      4  4  0  2  2  2  2
  subroutine img_components ( img_in, img_out, component_num )
    type(image), intent(in)  :: img_in  !it has to be binary! (or at least integer)
    type(image), intent(out) :: img_out
    integer,     intent(out) :: component_num
    integer :: m, n, ldim(3), b
    real, allocatable :: rmat(:,:,:), ccmat(:,:,:)
    integer, allocatable ::  p(:), q(:)
    integer :: component, i,j, north, west

    ldim = img_in%get_ldim()
    m = ldim(1)  !save them just for confort
    n = ldim(2)
    rmat = img_in%get_rmat()
    allocate(ccmat(ldim(1), ldim(2),1), source = 0.) !Initialization
    allocate(p(m*n), q(m*n),            source = 0 )
    component_num = 0
    call img_out%new(ldim,1.)
  !  p is used to store the component labels.  The dimension used is an absurd overestimate.
    p(1: m*n) = [(i, i = 1, m*n)]
  !  "Read" the array one pixel at a time.  If a (nonzero) pixel's north or
  !  west neighbor already has a label, the current pixel inherits it.
  !  In case the labels disagree, we need to adjust the P array so we can
  !  later deal with the fact that the two labels need to be merged.
    do i = 1, m
      do j = 1, n
        if ( i == 1 ) then
          north = 0
        else
          north = ccmat(i-1,j,1)
        end if
        if ( j == 1 ) then
          west = 0
        else
          west = ccmat(i,j-1,1)
        end if
        if ( rmat(i,j,1) /= 0 ) then
          if ( north == 0 ) then
            if ( west == 0 ) then
              component_num = component_num + 1
              ccmat(i,j,1) = component_num
            else
              ccmat(i,j,1) = west
            end if
          else if ( north /= 0 ) then
            if ( west == 0 .or. west == north ) then
              ccmat(i,j,1) = north
            else
              ccmat(i,j,1) = min ( north, west )
              if ( north < west ) then
                p(west) = north
              else
                p(north) = west
              end if
           end if
          end if
        end if
      end do
    end do
  !  When a component has multiple labels, have the higher labels point to the lowest one.
    do component = component_num, 1, -1
      b = component
      do while ( p(b) /= b )
        b = p(b)
      end do
      p(component) = b
    end do
  !  Locate the minimum label for each component. Assign these mininum labels new consecutive indices.
    q(0:component_num) = 0
    i = 0
    do component = 1, component_num
      if ( p(component) == component ) then
        i = i + 1
        q(component) = i
      end if
    end do
    component_num = i
  !  Replace the labels by consecutive labels.
    do i = 1, m
      do j = 1, n
        ccmat(i,j,1) = q ( p ( int(ccmat(i,j,1)) ) )
      end do
    end do
    call img_out%set_rmat(ccmat)
    return
  end subroutine img_components

  !The output of this function is the crossing number of a pixel as defined in
  !'Binary digital image processing' (book).
  !The crossing number indicates the number of 4-connected components in the 8-neighbourhoods.
  !It takes in input the 8-neighbourhoods of a px in a BINARY image.
  !It doesn't check if the image is binarised, but it should be for definition.
  function crossing_number(neigh_8) result(cn)
      real, intent(in) :: neigh_8(:)
      real    :: cn         !crossing number of the pixel px
      real    :: prod, sum  !just for comfort
      integer :: i          !loop

      prod = 1.
      sum  = 0.
      do i = 1, size(neigh_8)-1  !-1 because I don't consider the pixel itself
          prod  = prod*neigh_8(i)
          if(i<size(neigh_8)-1) then
              sum = sum + abs(neigh_8(i+1)-neigh_8(i))
          endif
      enddo
      sum = sum + abs(neigh_8(1) - neigh_8(size(neigh_8)-1))
      sum = 0.5*sum
      cn = prod + sum
    end function crossing_number

  subroutine print_mat(matrix)
    integer, intent(in) :: matrix(:,:)
    integer             :: j, s(2)
    s = shape(matrix)
    do j = 1, s(1)
      print *, matrix(j,:)
    enddo
  end subroutine print_mat

  subroutine enumerate_white_pixels(self, tmp_matrix)
    type(image), intent(in)        :: self
    real, allocatable, intent(out) :: tmp_matrix(:,:,:)
    real, allocatable :: rmat(:,:,:)
    integer           :: i, j, ldim(3), cnt

    rmat = self%get_rmat()
    ldim = self%get_ldim()
    allocate(tmp_matrix(ldim(1),ldim(2),1), source = 0.)
    cnt = 0
    do i = 1, ldim(1)
      do j = 1, ldim(2)
        if(rmat(i,j,1) /= 0) then
          cnt = cnt + 1
          tmp_matrix(i,j,1) = cnt
        endif
      enddo
    enddo
  end subroutine enumerate_white_pixels



  !Returns 8-neighborhoods of the pixel px in the matrix mat, in particular it returns
  !the intensity values of the 8-neigh ina CLOCKWISE order, starting from any 4-neigh
  !and the value of the pixel itself (the last one)
  function calc_neigh_8(mat, px) result(neigh_8)
    real,    intent(in)  :: mat(:,:,:)
    integer, intent(in)  :: px(3)
    integer, allocatable :: neigh_8(:)
    integer              :: i, j, ldim(3)
    ldim = shape(mat)
    if(px(3) /= 1) then
        print *, "The matrix has to be 2D!"
        stop
    endif
    i = px(1)
    j = px(2)            !Assumes to have a 2-dim matrix
    if ( i-1 < 1 .and. j-1 < 1 ) then
      allocate(neigh_8(4), source = 0)
      neigh_8(1) = mat(i+1,j,1)
      neigh_8(2) = mat(i+1,j+1,1)
      neigh_8(3) = mat(i,j+1,1)
      neigh_8(4) = mat(i,j,1)   !the pixel itself
    else if (j+1 > ldim(2) .and. i+1 > ldim(1)) then
      allocate(neigh_8(4), source = 0)
      neigh_8(1) = mat(i-1,j,1)
      neigh_8(2) = mat(i-1,j-1,1)
      neigh_8(3) = mat(i,j-1,1)
      neigh_8(4) = mat(i,j,1)   !the pixel itself
    else if (j-1 < 1  .and. i+1 >ldim(1)) then
      allocate(neigh_8(4), source = 0)
      neigh_8(3) = mat(i-1,j,1)
      neigh_8(2) = mat(i-1,j+1,1)
      neigh_8(1) = mat(i,j+1,1)
      neigh_8(4) = mat(i,j,1)   !the pixel itself
    else if (j+1 > ldim(2) .and. i-1 < 1) then
      allocate(neigh_8(4), source = 0)
      neigh_8(1) = mat(i,j-1,1)
      neigh_8(2) = mat(i+1,j-1,1)
      neigh_8(3) = mat(i+1,j,1)
      neigh_8(4) = mat(i,j,1)   !the pixel itself
    else if( j-1 < 1 ) then
      allocate(neigh_8(6), source = 0)
      neigh_8(5) = mat(i-1,j,1)
      neigh_8(4) = mat(i-1,j+1,1)
      neigh_8(3) = mat(i,j+1,1)
      neigh_8(2) = mat(i+1,j+1,1)
      neigh_8(1) = mat(i+1,j,1)
      neigh_8(6) = mat(i,j,1)   !the pixel itself
    else if ( j+1 > ldim(2) ) then
      allocate(neigh_8(6), source = 0)
      neigh_8(1) = mat(i-1,j,1)
      neigh_8(2) = mat(i-1,j-1,1)
      neigh_8(3) = mat(i,j-1,1)
      neigh_8(4) = mat(i+1,j-1,1)
      neigh_8(5) = mat(i+1,j,1)
      neigh_8(6) = mat(i,j,1)   !the pixel itself
    else if ( i-1 < 1 ) then
      allocate(neigh_8(6), source = 0)
      neigh_8(1) = mat(i,j-1,1)
      neigh_8(2) = mat(i+1,j-1,1)
      neigh_8(3) = mat(i+1,j,1)
      neigh_8(4) = mat(i+1,j+1,1)
      neigh_8(5) = mat(i,j+1,1)
      neigh_8(6) = mat(i,j,1)   !the pixel itself
    else if ( i+1 > ldim(1) ) then
      allocate(neigh_8(6), source = 0)
      neigh_8(1) = mat(i,j+1,1)
      neigh_8(2) = mat(i-1,j+1,1)
      neigh_8(3) = mat(i-1,j,1)
      neigh_8(4) = mat(i-1,j-1,1)
      neigh_8(5) = mat(i,j-1,1)
      neigh_8(6) = mat(i,j,1)   !the pixel itself
    else
      allocate(neigh_8(9), source = 0)
      neigh_8(1) = mat(i-1,j-1,1)
      neigh_8(2) = mat(i,j-1,1)
      neigh_8(3) = mat(i+1,j-1,1)
      neigh_8(4) = mat(i+1,j,1)
      neigh_8(5) = mat(i+1,j+1,1)
      neigh_8(6) = mat(i,j+1,1)
      neigh_8(7) = mat(i-1,j+1,1)
      neigh_8(8) = mat(i-1,j,1)
      neigh_8(9) = mat(i,j,1)   !the pixel itself
    endif
  end function calc_neigh_8
  !
  ! !Returns 8-neighborhoods of the pixel px in the matrix mat, it Returns
  ! !the intensity values and the value of the pixel itself
  ! function calc_neigh_8(mat, px) result(neigh_8)
  !   real,    intent(in)  :: mat(:,:,:)
  !   integer, intent(in)  :: px(3)
  !   integer, allocatable :: neigh_8(:)
  !   integer              :: i, j, ldim(3)
  !   ldim = shape(mat)
  !   if(px(3) /= 1) then
  !       print *, "The matrix has to be 2D!"
  !       stop
  !   endif
  !   i = px(1)
  !   j = px(2)            !Assumes to have a 2-dim matrix
  !   if ( i-1 < 1 .and. j-1 < 1 ) then
  !     allocate(neigh_8(4), source = 0)
  !     neigh_8(1) = mat(i+1,j,1)
  !     neigh_8(2) = mat(i+1,j+1,1)
  !     neigh_8(3) = mat(i,j+1,1)
  !     neigh_8(4) = mat(i,j,1)   !the pixel itself
  !   else if (j+1 > ldim(2) .and. i+1 > ldim(1)) then
  !     allocate(neigh_8(4), source = 0)
  !     neigh_8(1) = mat(i-1,j,1)
  !     neigh_8(2) = mat(i-1,j-1,1)
  !     neigh_8(3) = mat(i,j-1,1)
  !     neigh_8(4) = mat(i,j,1)   !the pixel itself
  !   else if (j-1 < 1  .and. i+1 >ldim(1)) then
  !     allocate(neigh_8(4), source = 0)
  !     neigh_8(1) = mat(i-1,j,1)
  !     neigh_8(2) = mat(i-1,j+1,1)
  !     neigh_8(3) = mat(i,j+1,1)
  !     neigh_8(4) = mat(i,j,1)   !the pixel itself
  !   else if (j+1 > ldim(2) .and. i-1 < 1) then
  !     allocate(neigh_8(4), source = 0)
  !     neigh_8(1) = mat(i,j-1,1)
  !     neigh_8(2) = mat(i+1,j-1,1)
  !     neigh_8(3) = mat(i+1,j,1)
  !     neigh_8(4) = mat(i,j,1)   !the pixel itself
  !   else if( j-1 < 1 ) then
  !     allocate(neigh_8(6), source = 0)
  !     neigh_8(1) = mat(i-1,j,1)
  !     neigh_8(2) = mat(i-1,j+1,1)
  !     neigh_8(3) = mat(i,j+1,1)
  !     neigh_8(4) = mat(i+1,j+1,1)
  !     neigh_8(5) = mat(i+1,j,1)
  !     neigh_8(6) = mat(i,j,1)   !the pixel itself
  !   else if ( j+1 > ldim(2) ) then
  !     allocate(neigh_8(6), source = 0)
  !     neigh_8(1) = mat(i-1,j,1)
  !     neigh_8(2) = mat(i-1,j-1,1)
  !     neigh_8(3) = mat(i,j-1,1)
  !     neigh_8(4) = mat(i+1,j-1,1)
  !     neigh_8(5) = mat(i+1,j,1)
  !     neigh_8(6) = mat(i,j,1)   !the pixel itself
  !   else if ( i-1 < 1 ) then
  !     allocate(neigh_8(6), source = 0)
  !     neigh_8(1) = mat(i,j-1,1)
  !     neigh_8(2) = mat(i+1,j-1,1)
  !     neigh_8(3) = mat(i+1,j,1)
  !     neigh_8(4) = mat(i+1,j+1,1)
  !     neigh_8(5) = mat(i,j+1,1)
  !     neigh_8(6) = mat(i,j,1)   !the pixel itself
  !   else if ( i+1 > ldim(1) ) then
  !     allocate(neigh_8(6), source = 0)
  !     neigh_8(1) = mat(i,j+1,1)
  !     neigh_8(2) = mat(i-1,j+1,1)
  !     neigh_8(3) = mat(i-1,j,1)
  !     neigh_8(4) = mat(i-1,j-1,1)
  !     neigh_8(5) = mat(i,j-1,1)
  !     neigh_8(6) = mat(i,j,1)   !the pixel itself
  !   else
  !     allocate(neigh_8(9), source = 0)
  !     neigh_8(1) = mat(i-1,j-1,1)
  !     neigh_8(2) = mat(i,j-1,1)
  !     neigh_8(3) = mat(i+1,j-1,1)
  !     neigh_8(4) = mat(i+1,j,1)
  !     neigh_8(5) = mat(i+1,j+1,1)
  !     neigh_8(6) = mat(i,j+1,1)
  !     neigh_8(7) = mat(i-1,j+1,1)
  !     neigh_8(8) = mat(i-1,j,1)
  !     neigh_8(9) = mat(i,j,1)   !the pixel itself
  !   endif
  ! end function calc_neigh_8

function cft_matrix(mat1,mat2) result(yes_no)
  real, intent(in) :: mat1(:,:,:), mat2(:,:,:)
  logical :: yes_no
  integer :: s(3), i, j, k
  s = shape(mat1)
  yes_no = .true.
  !if(s /= shape(mat2)) stop 'Input matrices have different dimensions'
  do i = 1, s(1)
     do j = 1, s(2)
        do k = 1, s(3)
            if(mat1(i,j,k) /= mat2(i,j,k)) then
                yes_no = .false.
                return
            endif
        enddo
      enddo
  enddo
end function cft_matrix

subroutine find_connected_comps(img_in, img_out)
type(image), intent(in)    :: img_in
type(image), intent(inout) :: img_out
real, allocatable :: rmat_in(:,:,:), tmp_mat(:,:,:),tmp_mat_cfr(:,:,:), neigh_8(:)
integer           :: i, j, n_it, ldim(3), n_maxit

rmat_in = img_in%get_rmat()
ldim    = img_in%get_ldim()
call img_out%new(ldim,1.)
allocate(tmp_mat_cfr(ldim(1),ldim(2),1), source = 0.)
call enumerate_white_pixels(img_in, tmp_mat)
n_maxit = maxval(tmp_mat) !maybe I can improve it
do n_it = 1, n_maxit
  if(cft_matrix(tmp_mat_cfr,tmp_mat)) then
    call img_out%set_rmat(tmp_mat)
    return
  endif
  tmp_mat_cfr = tmp_mat
  do i = 1, ldim(1)
    do j = 1, ldim(2)
      if(rmat_in(i,j,1) /= 0) then
        neigh_8 = calc_neigh_8(tmp_mat,[i,j,1])
        tmp_mat(i,j,1) = minval(neigh_8, neigh_8 /= 0)
      endif
    enddo
  enddo
enddo
call img_out%set_rmat(tmp_mat)
deallocate(rmat_in,tmp_mat, tmp_mat_cfr)
end subroutine find_connected_comps

function is_picked( new_coord,saved_coord, saved ) result(yes_no)
  integer, intent(in) :: new_coord(2)       !Coordinates of a new window to extract
  integer, intent(in) :: saved_coord(:,:)   !Coordinates of extracted windows
  integer, intent(in), optional :: saved              !How many particles have already been saved
  logical :: yes_no
  integer :: iwind, ssaved, s(2), thresh

  thresh = 1
  s = shape(saved_coord)
  if(s(2) /= 2) stop 'Dimension error'
  yes_no = .false.
  ssaved = s(1)
  if(present(saved)) ssaved = saved
  do iwind = 1, ssaved
      if(  (new_coord(1)-saved_coord(iwind,1))**2 + (new_coord(2)-saved_coord(iwind,2)) **2 < thresh) then
        yes_no = .true.
        return
      endif
  enddo
end function is_picked
end module connected_components

program simple_test_chiara_recentering
  include 'simple_lib.f08'
  use simple_micops
  use simple_image
  use simple_stackops
  use simple_math
  use connected_components
  logical              :: picked
  type(image)          :: img4try_in
  integer              :: matrix(25,2)

  matrix = reshape([ 1,1,1,0,0,6,6, &
                   & 1,1,0,0,6,6,6, &
                   & 1,0,0,2,0,6,0, &
                   & 0,0,2,2,0,0,4, &
                   & 0,5,0,0,0,4,4, &
                   & 0,5,5,5,0,0,0, &
                   & 0,5,5,0,0,3,3,3, 6],[25,2])

  call print_mat(matrix)
  picked = is_picked([3,4],matrix,3)
  print *, '3 4 is picked: ', picked

  picked = is_picked([1,0],matrix,3)
  print *, '1 0 is picked: ', picked

  picked = is_picked([1,0],matrix,1)
  print *, '1 0 is picked: ', picked

  picked = is_picked([0,0],matrix,1)
  print *, '0 0 is picked: ', picked

  picked = is_picked([0,0],matrix)
  print *, '0 0 is picked: ', picked
end program simple_test_chiara_recentering
