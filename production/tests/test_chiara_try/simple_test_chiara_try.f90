module simple_test_chiara_try_mod
  include 'simple_lib.f08'
  use simple_image, only : image
  implicit none
  public
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

  call process_ps_stack('pspecs_saga_polii.mrc', 'analysed_pspecs_saga_polii.mrc', 1.14, 35.)

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
