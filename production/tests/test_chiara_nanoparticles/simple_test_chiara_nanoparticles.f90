module simple_test_chiara_nanoparticles_mod
    use simple_image
    use simple_picker_chiara
    use simple_math
    implicit none
        real, allocatable :: centers(:,:)
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
  subroutine calc_aspect_ratio(img_cc, label, ratio)
      type(image), intent(inout) :: img_cc
      integer,     intent(in)    :: label
      real,        intent(out)   :: ratio
      integer, allocatable :: imat(:,:,:)
      integer, allocatable :: pos(:,:)
      real,    allocatable :: rmat(:,:,:)
      logical, allocatable :: border(:,:,:)
      integer :: ldim(3)
      real :: shortest_dist, longest_dist
      type(image) :: img_aux
      logical, parameter :: DEBUG = .true.
      rmat = img_cc%get_rmat()
      ldim = img_cc%get_ldim()
      allocate(imat(ldim(1), ldim(2), ldim(3)), source = nint(rmat))
      call img_aux%new(ldim, img_cc%get_smpd())
      where(abs(rmat-real(label)) > TINY) rmat = 0.
      !where (imat .ne. label) imat = 0
      call img_aux%set_rmat(rmat)
      call img_aux%masscen(centers(:3, counter))
      !call img_cc%grow_bin()
      ! print *, 'GROWN BORDERS '
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
          print *, '>>>>>>>>'
          print *, 'CC # ', label
          print *, 'shortest dist = ', shortest_dist
          print *, 'longest dist = ', longest_dist
          print *, 'RATIO = ', ratio
      endif
      call img_aux%kill
      deallocate(rmat, border, imat, pos)
  end subroutine calc_aspect_ratio
end module simple_test_chiara_nanoparticles_mod

program simple_test_chiara_nanoparticles
  include 'simple_lib.f08'
  use gnufor2
  use simple_picker_chiara
  use simple_test_chiara_nanoparticles_mod
  type(image)       :: img, img_cc, img_bin, img_p, img_cc2
  real, allocatable :: rmat(:,:,:), rmat_t(:,:,:)
  integer :: i, j, ldim(3), n_vol
  real :: smpd, seleted_t(1)
  integer :: h, k
  real :: thresh
  integer, allocatable :: sz(:,:)
  real :: step
  integer ::ncc
  real :: x_thresh(20), y_med(20)
  integer, allocatable :: pos(:,:), imat(:,:,:)
  real :: dist
  real,    allocatable :: ratios(:)
  logical, allocatable :: border(:,:,:)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!INITIALIZATIONS AND ALLOCATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call img%new    ([160,160,160], 0.358)
  call img_bin%new([160,160,160], 0.358)
  call img_cc%new ([160,160,160], 0.358)
  ! call img_p%new  ([162,162,162], 0.358)
  ! call img_cc%new ([162,162,162], 0.358)
  ! call img_cc2%new ([162,162,162], 0.358)
  ldim = img%get_ldim()
  allocate(rmat_t(ldim(1),ldim(2),ldim(3)), source = 0.)
  !!!!!!!!!!!!!!!!!!NANOPARTICLES BINARIZATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! do n_vol = 1, 7
  !     call img%read('particle'//int2str(n_vol)//'.mrc')
  !     rmat = img%get_rmat()
  !     step = maxval(rmat)/20.
  !     do i = 1, 19
  !         thresh = real(i)*step
  !         where(rmat > thresh)
  !             rmat_t = 1.
  !         elsewhere
  !             rmat_t = 0.
  !         endwhere
  !         call img_bin%set_rmat(rmat_t)
  !         !call img_bin%pad(img_p)  !padding for the calculation of the neigh_8
  !         !call img_p%find_connected_comps(img_cc)
  !         call img_bin%find_connected_comps(img_cc)
  !         sz          = img_cc%size_connected_comps()
  !         x_thresh(i) = thresh
  !         y_med(i)    = median(real(sz(2,:)))
  !     enddo
  !     seleted_t(:) = x_thresh(maxloc(y_med))
  !     print *,'VOLUME ', n_vol,  'SELECTED THRESHOLD = ', seleted_t(1)
  !     where(rmat > seleted_t(1))
  !         rmat_t = 1.
  !     elsewhere
  !         rmat_t = 0.
  !     endwhere
  !     call img_bin%set_rmat(rmat_t)
  !     call img_bin%write(int2str(n_vol)//'BINparticle.mrc')
  ! enddo
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!ASPECT RATIOS CALCULATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do n_vol = 1, 7
    call img_cc%read(int2str(n_vol)//'BINparticle.mrc')
    call img_cc%elim_cc([1, 233540025])
    call img_cc%bin(0.5) !binarize to re-calculate the cc, to have them in order again
    img_bin = img_cc
    call img_bin%find_connected_comps(img_cc) !cc re-calculation
    call img_cc%write(int2str(n_vol)//'PolishedOrdered1801.mrc')
    rmat   = img_cc%get_rmat()
    allocate(ratios    (nint(maxval(rmat))), source = 0.)
    allocate(centers (3,nint(maxval(rmat))), source = 0.)
    counter = 0
    do i = 1, nint(maxval(rmat))
        counter = counter + 1
        call calc_aspect_ratio(img_cc, i, ratios(i))
    enddo
    call hist(ratios, 20)
    deallocate(ratios, centers)
enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! For 3d plot of the identified centers
  ! rmat_t = 0.
  ! do i = 1, nint(maxval(rmat))
  !     rmat_t(ldim(1)/2+nint(centers(1,i)),ldim(2)/2+nint(centers(2,i)),ldim(3)/2+nint(centers(3,i))) = 1.
  ! enddo
  ! call img%set_rmat(rmat_t)
  ! call img%write('CentersParticle1MassCenter.mrc')
  !!!!!!!!!!!!!!!!!!!!!!!
end program simple_test_chiara_nanoparticles
