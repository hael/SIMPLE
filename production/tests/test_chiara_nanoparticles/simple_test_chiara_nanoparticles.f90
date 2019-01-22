module simple_test_chiara_nanoparticles_mod
    use simple_image
    use simple_picker_chiara
    use simple_math
    implicit none
        real, allocatable :: centers(:,:)
        integer, allocatable ::  centers_chia(:,:)
        integer :: counter
    contains
  ! This subroutine takes in input a connected component (cc) image
  ! and the label of one of its ccs and calculates its aspect ratio, which
  ! is defined as the ratio of the width and the height.
  ! The idea behind this is that the center of the cc is calculated,
  ! than everything is deleted except the borders of the cc. Finally,
  ! in order to calculate the width and the height, the min/max
  ! distances between the center and the borders are calculated. The
  ! aspect ratio is the ratio of those 2 distances.
  subroutine calc_aspect_ratio(img, img_cc, label, ratio)
      type(image), intent(inout) :: img
      type(image), intent(inout) :: img_cc
      integer,     intent(in)    :: label
      real,        intent(out)   :: ratio
      integer, allocatable :: imat(:,:,:)
      integer, allocatable :: pos(:,:)
      real,    allocatable :: rmat(:,:,:), rmat_cc(:,:,:)
      logical, allocatable :: border(:,:,:)
      integer :: ldim(3)
      real :: shortest_dist, longest_dist
      type(image) :: img_aux
      logical, parameter :: DEBUG = .true.
      rmat = img%get_rmat()
      rmat_cc = img_cc%get_rmat()
      ldim = img%get_ldim()
      allocate(imat(ldim(1), ldim(2), ldim(3)), source = nint(rmat_cc))
      call img_aux%new(ldim, img%get_smpd())
      where(abs(rmat_cc-real(label)) > TINY)
          rmat = 0.
          rmat_cc = 0.
      endwhere
      where (imat .ne. label) imat = 0
      call img_aux%set_rmat(rmat)   !NEW
      call img_aux%masscen(centers(:3, counter))
      centers_chia(:3,counter) = center_cc(imat)
      !call img_cc%grow_bin()
      ! write(logfhandle,*)'GROWN BORDERS '
      call img_cc%border_mask(border, label)
      where(border .eqv. .true.)
          imat = 1
      elsewhere
          imat = 0
      endwhere
      call get_pixel_pos(imat,pos)
      shortest_dist = points_dist(centers(:,counter)+real(ldim/2.),real(pos),'min')
      longest_dist  = points_dist(centers(:,counter)+real(ldim/2.),real(pos),'max')
      ratio = shortest_dist/longest_dist
      if(DEBUG) then
          write(logfhandle,*) '>>>>>>>>'
          write(logfhandle,*) 'CC # ', label
          write(logfhandle,*) 'shortest dist = ', shortest_dist
          write(logfhandle,*) 'longest dist = ', longest_dist
          write(logfhandle,*) 'RATIO = ', ratio
      endif
      call img_aux%kill
      deallocate(rmat, rmat_cc, border, imat, pos)
  end subroutine calc_aspect_ratio

  ! This subrotuine takes in input a nanoparticle and
  ! binarizes it by thresholding. The gray level histogram is split
  ! in 20 parts, which corrispond to 20 possible threshold.
  ! Among those threshold, the selected one is the for which
  ! ths size of the connected component has the maximum median value.
  ! The idea is that if the threshold is wrong, than the binarization
  ! produces a few huge ccs (connected components) and a lot of very small
  ! ones (dust). So the threshold with the maximum median value would
  ! correspond to the most consistent one, meaning that the size of the ccs
  ! would have a gaussian distribution.
  subroutine bin_nanopart(img)
      type(image), intent(inout) :: img
      type(image) :: img_bin ! binary image
      ! type(image) :: img_pad ! padded image, for the calc_neigh_8 routine in 3D
      type(image) :: img_cc  ! ccs image
      real, allocatable :: rmat(:,:,:), rmat_t(:,:,:)
      integer :: i
      real    :: step   !histogram disretization step
      real    :: thresh !binarization threshold
      real    :: seleted_t(1)
      real    :: x_thresh(19), y_med(19)
      real    :: smpd
      integer, allocatable :: sz(:,:)
      integer :: ldim(3)
      ldim = img%get_ldim()
      smpd = img%get_smpd()
      allocate(rmat_t(ldim(1), ldim(2), ldim(3)), source = 0.)
      rmat = img%get_rmat()
      call img_bin%new(ldim, smpd)
      call img_cc%new (ldim, smpd)
      step = maxval(rmat)/20.
      do i = 1, 19
          thresh = real(i)*step
          where(rmat > thresh)
              rmat_t = 1.
          elsewhere
              rmat_t = 0.
          endwhere
          call img_bin%set_rmat(rmat_t)
          call img_bin%find_connected_comps(img_cc)
          sz          = img_cc%size_connected_comps()
          x_thresh(i) = thresh
          y_med(i)    = median(real(sz(2,:)))
      enddo
      seleted_t(:) = x_thresh(maxloc(y_med))
      if(DEBUG) write(logfhandle,*)  'SELECTED THRESHOLD = ', seleted_t(1)
      where(rmat > seleted_t(1))
          rmat_t = 1.
      elsewhere
          rmat_t = 0.
      endwhere
      call img%set_rmat(rmat_t)
      deallocate(rmat, rmat_t, sz)
  end subroutine bin_nanopart

  subroutine aspect_ratios_hist(img,img_cc)
    use gnufor2
    type(image), intent(inout):: img
    type(image), intent(inout):: img_cc
    type(image) :: img_bin
    integer :: ldim(3)
    integer :: i, j
    real, allocatable :: ratios(:)
    real, allocatable :: rmat_cc(:,:,:)
    real, allocatable :: paral(:,:)
    real, allocatable :: rmat_t(:,:,:)
    ldim = img_cc%get_ldim()
    call img_cc%elim_cc([1, 100000])
    call img_cc%bin(0.5) !binarize to re-calculate the cc, to have them in order again
    img_bin = img_cc
    call img_bin%find_connected_comps(img_cc) !cc re-calculation
    rmat_cc   = img_cc%get_rmat()
    rmat_t = rmat_cc
    allocate(ratios    (nint(maxval(rmat_cc))), source = 0.)
    allocate(centers (3,nint(maxval(rmat_cc))), source = 0.) !global variable
    allocate(centers_chia (3,nint(maxval(rmat_cc))), source = 0) !global variable
    allocate(paral (3,nint(maxval(rmat_cc))), source = 0.) !global variable
    counter = 0
    do i = 1, nint(maxval(rmat_cc))
        counter = counter + 1
        call calc_aspect_ratio(img, img_cc, i, ratios(i))
    enddo
    call hist(ratios, 20)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! centers visualization
rmat_t = 0.
do i = 1, nint(maxval(rmat_cc))
    rmat_t(ldim(1)/2+nint(centers(1,i)),ldim(2)/2+nint(centers(2,i)),ldim(3)/2+nint(centers(3,i))) = 1.
enddo
call img_cc%set_rmat(rmat_t)
call img_cc%write('CentersParticle1MassCenter.mrc')
rmat_t = 0.
do i = 1, nint(maxval(rmat_cc))
    rmat_t(centers_chia(1,i),centers_chia(2,i),centers_chia(3,i)) = 1.
enddo
call img_cc%set_rmat(rmat_t)
call img_cc%write('CentersParticle1Chia.mrc')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    ! print *, 'CENTERS MASS NEW = '
    ! call vis_mat(transpose(centers))
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ! Checking if there r parallel vectors
!     ! 3D PLOT OF IDENTIFIED CENTERS
!     ! rmat_t = 0.
!     ! do i = 1, nint(maxval(rmat))
!     !     rmat_t(ldim(1)/2+nint(centers(1,i)),ldim(2)/2+nint(centers(2,i)),ldim(3)/2+nint(centers(3,i))) = 1.
!     ! enddo
!     ! call img_cc%set_rmat(rmat_t)
!     ! call img_cc%write('CentersParticle1MassCenter.mrc')
!     rmat_t = 0.
!     i = 55
!     !do i = 1, size(centers, dim = 2)
!         do j = 1, size(centers, dim = 2)
!             !if (j > i)then
!                 if(sqrt((centers(1,i)-centers(1,j))**2.+(centers(2,i)-centers(2,j))**2.) < 0.8) then
!                     print *, 'i = ', i, '& j = ', j, 'PARALLEL'
! rmat_t(ldim(1)/2+nint(centers(1,i)),ldim(2)/2+nint(centers(2,i)),ldim(3)/2+nint(centers(3,i))) = 1.
! rmat_t(ldim(1)/2+nint(centers(1,j)),ldim(2)/2+nint(centers(2,j)),ldim(3)/2+nint(centers(3,j))) = 1.
!                 endif
!             !endif
!         enddo
!     !enddo
!     call img_cc%set_rmat(rmat_t)
!     call img_cc%write('ParallelCenters55.mrc')
    deallocate(ratios, rmat_cc, centers)
  end subroutine aspect_ratios_hist
end module simple_test_chiara_nanoparticles_mod

program simple_test_chiara_nanoparticles
  include 'simple_lib.f08'
  use gnufor2
  use simple_picker_chiara
  use simple_test_chiara_nanoparticles_mod
  type(image) :: img, img_cc
  integer :: n_vol
  real, allocatable :: rmat(:,:,:)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!INITIALIZATIONS AND ALLOCATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call img%new    ([160,160,160], 0.358)
  call img_cc%new ([160,160,160], 0.358)
  !!!!!!!!!!!!!!!!!!NANOPARTICLES BINARIZATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! do n_vol = 1, 7
  !     call img%read('particle'//int2str(n_vol)//'.mrc')
  !     call bin_nanopart(img)
  !     call img%write(int2str(n_vol)//'BINparticle.mrc')
  ! enddo
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!ASPECT RATIOS CALCULATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do n_vol = 1, 1!7
    call img%read('particle'//int2str(n_vol)//'.mrc')
    call img_cc%read(int2str(n_vol)//'BINparticle.mrc')
    call aspect_ratios_hist(img,img_cc)
enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D PLOT OF IDENTIFIED CENTERS
  ! rmat_t = 0.
  ! do i = 1, nint(maxval(rmat))
  !     rmat_t(ldim(1)/2+nint(centers(1,i)),ldim(2)/2+nint(centers(2,i)),ldim(3)/2+nint(centers(3,i))) = 1.
  ! enddo
  ! call img%set_rmat(rmat_t)
  ! call img%write('CentersParticle1MassCenter.mrc')
  !!!!!!!!!!!!!!!!!!!!!!!

 ! call img%read('1BINparticle.mrc')
 ! rmat = img%get_rmat()
 ! rmat(1:50, )

end program simple_test_chiara_nanoparticles
