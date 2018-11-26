module simple_test_chiara_try_mod
  include 'simple_lib.f08'
  use simple_image, only : image
  implicit none
  public
  ! #include "simple_local_flags.inc"

  contains
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
    !This function plots the gray-level histogram of a 2D image
    ! COULD IMPLEMENT OTHER OPTIONS FOR THE HIST
    subroutine gray_level_hist(img)
        use gnufor2
        class(image), intent(inout) :: img
        ! real(kind=4),    optional   :: pause
        ! character(len=*),optional   :: color, terminal, filename, persist, input
        real, allocatable :: rmat(:,:,:)
        real, allocatable :: x(:) !vectorisation of the matrix
        integer :: n ! number of intervals
        rmat = img%get_rmat()
        x = pack(rmat(:,:,:), .true.)
        n = int(maxval(rmat))-int(minval(rmat)) + 1
        call hist(x, n)
        deallocate(rmat,x)
    end subroutine gray_level_hist

      ! subroutine hist(x,n,pause,color,terminal,filename,persist,input)
      !
      !     ! this subroutine plots the histogram of data contained in array x, using n bins
      !
      !     implicit none
      !     real(kind=4), intent(in)  :: x(:) !the data to plot
      !     integer, intent(in)       :: n !the number of intervals
      !     real(kind=4), optional    :: pause
      !     character(len=*),optional :: color, terminal, filename, persist, input
      !     integer                   :: i, j, ierror, ios, file_unit, nx
      !     character(len=100)        :: data_file_name, command_file_name, yrange, xrange1, xrange2, my_color, &
      !         & xtic_start, dxtic, xtic_end, my_pause, my_persist
      !     real(kind=4)              :: xmin, xmax, xhist(0:n), yhist(n+1), dx
end module simple_test_chiara_try_mod

program simple_test_chiara_try
  include 'simple_lib.f08'
  use simple_powerspec_analysis
  use simple_oris
  use gnufor2
  use simple_picker_chiara
  use simple_micops
  use simple_image
  use simple_stackops
  use simple_math
  use simple_test_chiara_try_mod
  use simple_edge_detector
  use simple_parameters, only: parameters
  use simple_cmdline,    only: cmdline
type(image) :: img, img_out
real, allocatable :: rmat(:,:,:)
integer :: i, winsz
logical ::yes_no
real :: matrix(7,7,1), thresh, lp

 ! img = build_ice_template(512, 1.41, winsz)
 ! print *, 'WINSZ = ', winsz
 ! call img_win%new([int(winsz/2),int(winsz/2),1], 1.41)
 ! call img%write('IceTemplate.mrc')
 ! call img%window_slim([int(winsz/4),int(winsz/4)], 10, img_win, outside )
 ! call img_win%write('IceTemplateWin.mrc')
 ! print *, 'OUTSIDE = ', outside

 ! call img%ellipse([256,256],[20.,20.], 'yes')
 ! call img%ellipse([352,312],[10.,5.], 'yes')
 ! call img%ellipse([23,25],[5.,10.], 'yes')
 ! call img%ellipse([112,53],[10.,10.], 'yes')
 !  call img%ellipse([220,153],[8.,8.], 'yes')
 ! call img%read('ToyImage.mrc')
 ! call img%find_connected_comps(img_cc)
 ! call img_cc%write('ToyImageCC.mrc')
 ! yes_no = is_symmetric(img_cc, 1)
 ! print *, 'CC ', 1, ' is symmetric ', yes_no
 ! yes_no = is_symmetric(img_cc, 2)
 ! print *, 'CC ', 2, ' is symmetric ', yes_no
 ! yes_no = is_symmetric(img_cc, 3)
 ! print *, 'CC ', 3, ' is symmetric ', yes_no
 ! yes_no = is_symmetric(img_cc, 4)
 ! print *, 'CC ', 4, ' is symmetric ', yes_no
 ! yes_no = is_symmetric(img_cc, 5)
 ! print *, 'CC ', 5, ' is symmetric ', yes_no

  ! !TO START AGAIN FROM HERE
  ! call img%new([512,512,1],1.)
  ! call img_cc%new([512,512,1],1.)
  ! call img%read('SAGAWhithICEBin.mrc')
  ! call img%find_connected_comps(img_cc)
  ! call img_cc%write('SAGAWithICEBinCC.mrc')
  ! rmat = img_cc%get_rmat()
  ! do i = 1, int(maxval(rmat))
  !   yes_no = is_symmetric(img_cc, i)
  !   print *, 'CC ', i, ' is symmetric ', yes_no
  ! enddo

 ! call process_ps_stack('pspecs_sphire_tstdat.mrc', 'analysed_pspecs_sphire.mrc', 1.41, 20.,1)
! call process_ps_stack('pspecs_saga_polii.mrc', 'analysed_pspecs_saga_polii.mrc', 1.14, 35.,2)


! call img%new([2048,2048,1],1.)
! !call img%read('/home/lenovoc30/Desktop/MassCenter/simple_test_chiara_mass_center/shrunken_hpassfiltered.mrc')
! call img%read('ImgNOnoise.mrc')
! call iterative_thresholding(img,img_out,thresh)
! call img%write('Img_inParticle.mrc')
! call img_out%write('Img_outParticle_iterative_thresholding.mrc')
! print *, 'selected threshold: ' , thresh
! call otsu(img,img_out,thresh)
! call img_out%write('Img_outParticle_otsu.mrc')
! print *, 'selected threshold: ' , thresh
! 
! call img%new([1024,1024,1],1.)
! call img%read('/home/lenovoc30/Desktop/PickingResults/SomeExamples/NegStainingWorking/try23Nov/bp_filtered.mrc')
! call img%NLmean()
! call img%write('NLmeanFiltered.mrc')
end program simple_test_chiara_try
!In here you have to implement the statistics.

! !!!!!!!!!!!!ATTENTION TO PREPARE_CC, YOU MIGHT NEED TO ENUMERATE THEM FIRST, OR USE ELIMIN_CC
!     call img%new([70,70,1],1.)
!     allocate(diameter(74), source = 0.)
!     allocate(n_pixels(74), source = 0)
!     allocate(my_label(74))
!     do i = 1, 74
!         call img%read('centered_particles_BIN.mrc', i)
!         diameter(i) = part_diameter(img)
!         n_pixels(i) = img%nforeground()
!     enddo
!     n = int(maxval(diameter)) + 1 - int(minval(diameter))!aint((maxval(diameter)+1-minval(diameter))/74)
!     call hist(diameter,n, color = 'red', persist ='yes')
!     call elim_dup(diameter, ddiameter)   !IN ORDER TO USE SORTMEANS I HAVE TO ELIM_DUP FROM THE VECTOR
!     allocate(means(size(ddiameter)), source = 0.)
!     call sortmeans(ddiameter,MAXITS,means,labels)
!     deallocate(means, labels)
!
!     n = int(maxval(n_pixels)) + 1 - int(minval(n_pixels))!aint((maxval(diameter)+1-minval(diameter))/74)
!     call hist(real(n_pixels),n, color = 'blue', persist ='yes')
!     call elim_dup(real(n_pixels), nn_pixels)   !IN ORDER TO USE SORTMEANS I HAVE TO ELIM_DUP FROM THE VECTOR
!     allocate(means(size(nn_pixels)), source = 0.)
!     call sortmeans(nn_pixels,MAXITS,means,labels)
!
!     !Let s 'come back'
!     where(diameter <= 22. .or. n_pixels < 200.)
!       my_label='noise'
!     elsewhere(diameter > 44. .or. n_pixels > 600.)
!       my_label = 'aggregation'
!     elsewhere
!       my_label = 'particle'
!     endwhere
!     cnt_noise = 0
!     cnt_particle = 0
!     cnt_aggregation = 0
!     do i = 1, 74
!       if(my_label(i) .eq. 'particle') then
!         cnt_particle = cnt_particle + 1
!         call img%read('centered_particles.mrc', i)
!         call img%write('Particle.mrc', cnt_particle)
!       elseif(my_label(i) .eq. 'noise') then
!         cnt_noise = cnt_noise + 1
!         call img%read('centered_particles.mrc', i)
!         call img%write('Noise.mrc', cnt_noise)
!       elseif(my_label(i) .eq. 'aggregation') then
!         cnt_aggregation = cnt_aggregation + 1
!         call img%read('centered_particles.mrc', i)
!         call img%write('Aggregation.mrc', cnt_aggregation)
!       else
!         print *, 'Error!'
!       endif
!     enddo


! matrix = reshape(real([ 1,1,1,0,0,6,5, &
!                  & 1,1,0,0,6,6,6, &
!                  & 1,0,0,2,0,6,0, &
!                  & 0,0,2,2,0,0,4, &
!                  & 0,5,0,0,0,4,4, &
!                  & 0,5,5,5,0,0,0, &
!                  & 0,5,5,0,0,3,3]),[7,7,1])
! call img_in%disc( ldim, 1.0 , 10.)
