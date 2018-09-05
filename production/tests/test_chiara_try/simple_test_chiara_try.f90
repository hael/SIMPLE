module simple_test_chiara_try_mod
  include 'simple_lib.f08'
  use simple_image, only : image
  implicit none
  contains
      ! !!!!!!!!!!!!!!ADDED BY CHIARA!!!!!!!!!!!!!!!!!!
      ! ! To extract particles from micrograph self, basing on connected
      ! ! component (=: cc) image img_cc.
      ! ! The center of each particle is identified as the pixel which minimize
      ! ! the distance between itself and all the other pixels in the same cc.
      ! subroutine extract_windows_centered(self, img_cc, box)
      !     use simple_image
      !     type(image),  intent(inout) :: self        !original image
      !     type(image),  intent(inout) :: img_cc      !connected
      !     integer,      intent(in)    :: box         !size of the extracted windows
      !     type(image) :: imgwin, imgbin
      !     integer     :: n_cc               ! n_cc  = # cc in the input image
      !     integer     :: i, j, ldim(3)
      !     integer     :: n_window, n_px     !counters
      !     integer     :: idx(2)
      !     logical     :: outside
      !     real,    allocatable :: rmat_cc(:,:,:)   !corresponding matrix to img_cc
      !     real,    allocatable :: dist(:,:)        !to extract the window according to the px which minimize the dist
      !     logical, allocatable :: mask(:,:)        !to calc the min of an array along a specific dim
      !     integer, allocatable :: pos(:,:)         !position of the pixels of a fixed cc
      !     integer, allocatable :: rmat_mask(:,:,:) !indentify a specific cc
      !     ldim = self%get_ldim()
      !     rmat_cc = img_cc%get_rmat()
      !     allocate(rmat_mask(ldim(1), ldim(2), 1), source = 0)
      !     n_window = 0
      !     call imgwin%new([box,box,1], 1.)
      !     call imgbin%new([box,box,1], 1.)
      !     do n_cc = 1, int(maxval(rmat_cc))   !for every cc extract one window
      !         mask = .false.
      !         rmat_mask = 0
      !         where(abs(rmat_cc-real(n_cc))< TINY) rmat_mask = 1
      !         if(any(rmat_mask > 0.5)) then
      !             call get_pixel_pos(rmat_mask, pos)
      !             if(allocated(dist)) deallocate(dist)
      !             if(allocated(mask)) deallocate(mask)
      !             allocate(dist(3,size(pos, dim = 2)), source = 0.)
      !             allocate(mask(3,size(pos, dim = 2)), source = .false.)
      !             mask(3,:) = .true.
      !             n_px = 0
      !             do i = 1, ldim(1)
      !                 do j = 1, ldim(2)
      !                     if(rmat_mask(i,j,1) > 0.5) then
      !                         n_px = n_px + 1
      !                         dist( 3, n_px) = pixels_dist( [i,j,1], pos )
      !                         dist(:2, n_px) = [i,j]
      !                     endif
      !                 enddo
      !             enddo
      !             idx = minloc(dist, mask)
      !             call   self%window_slim(int(dist(:2, idx(2)))- box/2 , box, imgwin, outside )
      !           !  call img_cc%window_slim(int(dist(:2, idx(2)))- box/2 , box, imgbin, outside )
      !             if(.not. outside) then
      !                 n_window = n_window + 1
      !                 call imgwin%write('extract_windows_centered.mrc', n_window)
      !               !  call imgbin%write('extract_windows_centered_bin.mrc', n_window)
      !             endif
      !         endif
      !     enddo
      !     deallocate(pos, mask)
      ! end subroutine extract_windows_centered

      ! This function stores in pos the indeces corresponding to
      ! the pixels with value > 0 in the binary matrix rmat_masked.
      subroutine get_pixel_pos(rmat_masked, pos)
          integer,              intent(in)  :: rmat_masked(:,:,:)
          integer, allocatable, intent(out) :: pos(:,:)
          integer :: s(3), i, j, cnt
          ! if( any(rmat_masked > 1.0001) .or. any(rmat_masked < 0. ))&
          ! THROW_HARD('Input not binary; get_pixel_pos')
          s = shape(rmat_masked)
          allocate(pos(3, count(rmat_masked > 0.5)), source = 0)
          cnt = 0
          do i = 1, s(1)
                do j = 1, s(2)
                    if(rmat_masked(i,j,1) > 0.5) then !rmat_masked is binary
                        cnt = cnt + 1
                        pos(:,cnt) = [i,j,1]
                    endif
                enddo
          enddo
      end subroutine get_pixel_pos

      !>   calculates the euclidean distance between one pixel and a list of other pixels.
      ! if which == 'max' then distance is the maximum value of the distance between
      !              the selected pixel and all the others
      ! if which == 'sum' then distance is the sum of the distances between the
      !              selected pixel and all the others.
      function pixels_dist( px, vec, which) result( dist )
          integer, intent(in)           :: px(3)
          integer, intent(in)           :: vec(:,:)
          character(len=*),  intent(in) :: which
          real :: dist
          select case(which)
          case('max')
              dist =  maxval(sqrt((px(1)-vec(1,:))**2.+(px(2)-vec(2,:))**2+(px(3)-vec(3,:))**2))
          case('sum')
              dist =  maxval(sqrt((px(1)-vec(1,:))**2.+(px(2)-vec(2,:))**2+(px(3)-vec(3,:))**2))
          case DEFAULT
              write(*,*) 'Pixels_dist kind: ', trim(which)
              !THROW_HARD('Unsupported pixels_dist kind; pixels_dist')
          end select
      end function pixels_dist

      function part_diameter(img) result(diameter)
          type(image), intent(inout) :: img
          real, allocatable    :: rmat(:,:,:), dist(:)
          integer, allocatable :: rmat_masked(:,:,:), pos(:,:)
          real :: diameter
          integer :: ldim(3), i, j, cnt
          ldim = img%get_ldim()
          ! rmat = img%get_rmat()
          ! if(ldim(3) /= 1) THROW_HARD('This subroutine is for 2D images!; dilatation')
          ! if( any(rmat > 1.0001) .or. any(rmat < 0. ))&
          ! THROW_HARD('input for dilatation not binary; dilatation')
          ! deallocate(rmat)
          call img%prepare_connected_comps() !consider just one particle, not multiple
          allocate(rmat_masked(ldim(1),ldim(2),1), source = 0)
          rmat_masked = int(img%get_rmat())
          call get_pixel_pos(rmat_masked,pos)
          cnt = 0
          allocate(dist(count(rmat_masked > 0.5))) !one for each white px
          do i = 1, ldim(1)
              do j = 1, ldim(2)
                if(rmat_masked(i,j,1) > 0.5) then  !just ones
                  cnt = cnt + 1
                  dist(cnt) = pixels_dist([i,j,1],pos, 'max')
                endif
              enddo
          enddo
          diameter = maxval(dist)
      end function part_diameter
end module simple_test_chiara_try_mod

program simple_test_chiara_try
  include 'simple_lib.f08'
  use simple_oris
  use gnufor2
  use simple_picker_chiara
  use simple_micops
  use simple_image
  use simple_stackops
  use simple_math
  use simple_test_chiara_try_mod
  use simple_edge_detector, only : automatic_thresh_sobel
  use simple_parameters, only: parameters
  use simple_cmdline,    only: cmdline

real :: dat(10)
real, allocatable :: means(:), diameter(:)
real, allocatable :: data(:), ddiameter(:)
integer, allocatable :: labels(:), n_pixels(:)
real, allocatable ::  nn_pixels(:), rmat(:,:,:)
integer, parameter :: MAXITS = 100
integer :: i, j, n, cnt_noise, cnt_particle, cnt_aggregation
type(image) :: img, mic
character(len=15), allocatable:: my_label(:)

! dat = reshape(real([1,2,3,5,5,5,7,8,10,10]), [10])
! print *, 'dat = ', dat
! call hist(dat,10, pause = 1., color = 'red', persist ='yes')
! call elim_dup(dat, data)   !IN ORDER TO USE SORTMEANS I HAVE TO ELIM_DUP FROM THE VECTOR
! allocate(means(size(data)), source = 0.)
! call sortmeans(data,MAXITS,means,labels)
! call img%disc([64,64,1], 1.0 , 4.)
! diameter = part_diameter(img)
! print *, 'Diameter  4 = ', diameter
! call img%kill
! call img%disc([64,64,1], 1.0 , 3.)
! diameter = part_diameter(img)
! print *, 'Diameter 3 = ', diameter
! call img%kill
! call img%disc([64,64,1], 1.0 , 7.)
! diameter = part_diameter(img)
! print *, 'Diameter 7 = ', diameter

!!!!!!!!!!!!ATTENTION TO PREPARE_CC, YOU MIGHT NEED TO ENUMERATE THEM FIRST, OR USE ELIMIN_CC
    call img%new([70,70,1],1.)
    allocate(diameter(74), source = 0.)
    allocate(n_pixels(74), source = 0)
    allocate(my_label(74))
    do i = 1, 74
        call img%read('centered_particles_BIN.mrc', i)
        diameter(i) = part_diameter(img)
        n_pixels(i) = img%nforeground()
    enddo
    n = int(maxval(diameter)) + 1 - int(minval(diameter))!aint((maxval(diameter)+1-minval(diameter))/74)
    call hist(diameter,n, color = 'red', persist ='yes')
    call elim_dup(diameter, ddiameter)   !IN ORDER TO USE SORTMEANS I HAVE TO ELIM_DUP FROM THE VECTOR
    allocate(means(size(ddiameter)), source = 0.)
    call sortmeans(ddiameter,MAXITS,means,labels)
    deallocate(means, labels)

    n = int(maxval(n_pixels)) + 1 - int(minval(n_pixels))!aint((maxval(diameter)+1-minval(diameter))/74)
    call hist(real(n_pixels),n, color = 'blue', persist ='yes')
    call elim_dup(real(n_pixels), nn_pixels)   !IN ORDER TO USE SORTMEANS I HAVE TO ELIM_DUP FROM THE VECTOR
    allocate(means(size(nn_pixels)), source = 0.)
    call sortmeans(nn_pixels,MAXITS,means,labels)

    !Let s 'come back'
    where(diameter <= 22. .or. n_pixels < 200.)
      my_label='noise'
    elsewhere(diameter > 44. .or. n_pixels > 600.)
      my_label = 'aggregation'
    elsewhere
      my_label = 'particle'
    endwhere
    cnt_noise = 0
    cnt_particle = 0
    cnt_aggregation = 0
    do i = 1, 74
      if(my_label(i) .eq. 'particle') then
        cnt_particle = cnt_particle + 1
        call img%read('centered_particles.mrc', i)
        call img%write('Particle.mrc', cnt_particle)
      elseif(my_label(i) .eq. 'noise') then
        cnt_noise = cnt_noise + 1
        call img%read('centered_particles.mrc', i)
        call img%write('Noise.mrc', cnt_noise)
      elseif(my_label(i) .eq. 'aggregation') then
        cnt_aggregation = cnt_aggregation + 1
        call img%read('centered_particles.mrc', i)
        call img%write('Aggregation.mrc', cnt_aggregation)
      else
        print *, 'Error!'
      endif
    enddo
end program simple_test_chiara_try

! matrix = reshape(real([ 1,1,1,0,0,6,5, &
!                  & 1,1,0,0,6,6,6, &
!                  & 1,0,0,2,0,6,0, &
!                  & 0,0,2,2,0,0,4, &
!                  & 0,5,0,0,0,4,4, &
!                  & 0,5,5,5,0,0,0, &
!                  & 0,5,5,0,0,3,3]),[7,7,1])
! call img_in%disc( ldim, 1.0 , 10.)
! call build_ellipse(img_in,[10.,10.], [4.,4.], 0.)
