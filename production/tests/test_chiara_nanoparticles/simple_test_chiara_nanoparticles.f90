module simple_test_chiara_nanoparticles_mod
    use simple_image
    use simple_picker_chiara
    use simple_math
    implicit none
        real, allocatable :: centers(:,:)
        logical, parameter :: DEBUGG = .false.
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
      real    :: shortest_dist, longest_dist
      rmat    = img%get_rmat()
      rmat_cc = img_cc%get_rmat()
      ldim    = img%get_ldim()
      centers(:3,label) = atom_masscen(img,img_cc,label)
      allocate(imat(ldim(1),ldim(2),ldim(3)),source = nint(rmat_cc)) !for function get_pixel_pos
      call img_cc%border_mask(border, label)
      where(border .eqv. .true.)
          imat = 1
      elsewhere
          imat = 0
      endwhere
      call get_pixel_pos(imat,pos)
      shortest_dist = pixels_dist(centers(:,label), real(pos),'min')
      longest_dist  = pixels_dist(centers(:,label), real(pos),'max')
      if(abs(longest_dist) > TINY) then
          ratio = shortest_dist/longest_dist
      else
          ratio = 0.
          if(DEBUGG) write(logfhandle,*) 'cc ', label, 'LONGEST DIST = 0'
      endif
      if(DEBUGG) then
          write(logfhandle,*) '>>>>>>>>'
          write(logfhandle,*) 'CC # ', label
          write(logfhandle,*) 'shortest dist = ', shortest_dist
          write(logfhandle,*) 'longest dist = ', longest_dist
          write(logfhandle,*) 'RATIO = ', ratio
      endif
      deallocate(rmat, rmat_cc, border, imat, pos)
  end subroutine calc_aspect_ratio

  ! This subrotuine takes in input a nanoparticle and
  ! binarizes it by thresholding. The gray level histogram is split
  ! in 20 parts, which corrispond to 20 possible threshold.
  ! Among those threshold, the selected one is the for which
  ! the size of the connected component of the correspondent
  ! thresholded nanoparticle has the maximum median value.
  ! The idea is that if the threshold is wrong, than the binarization
  ! produces a few huge ccs (connected components) and a lot of very small
  ! ones (dust). So the threshold with the maximum median value would
  ! correspond to the most consistent one, meaning that the size of the ccs
  ! would have a gaussian distribution.
  subroutine bin_nanopart(img)
      type(image), intent(inout) :: img
      type(image) :: img_bin   ! binary image
      type(image) :: img_cc    ! ccs image
      real, allocatable :: rmat(:,:,:), rmat_t(:,:,:)
      real    :: step          !histogram disretization step
      real    :: thresh        !binarization threshold
      real    :: seleted_t(1)  !selected threshold for nanoparticle binarization
      real    :: x_thresh(19), y_med(19)
      real    :: smpd
      integer, allocatable :: sz(:) !size of the ccs and correspondent label
      integer :: ldim(3), i
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
          y_med(i)    = median(real(sz))
      enddo
      seleted_t(:) = x_thresh(maxloc(y_med))
      if(DEBUGG) write(logfhandle,*)  'SELECTED THRESHOLD = ', seleted_t(1)
      where(rmat > seleted_t(1))
          rmat_t = 1.
      elsewhere
          rmat_t = 0.
      endwhere
      call img%set_rmat(rmat_t)
      deallocate(rmat, rmat_t, sz)
  end subroutine bin_nanopart

  ! subroutine full_atom_stripes(img_cc, img_stripes)
  !     type(image), intent(inout) :: img_cc
  !     type(image), intent(inout) :: img_stripes
  !     integer :: ldim(3), n_stripe
  !     real, allocatable :: rmat_cc(:,:,:), rmat_stripes(:,:,:)
  !     integer :: pos(3)
  !     integer :: label !integer for consistency with the other cases
  !
  !     ldim = img_cc%get_ldim()
  !     if(ldim .ne. img_stripes%get_ldim()) print *, "Non conforming dims for input images; full_atom_stripes"
  !     rmat_cc      = img_cc%get_rmat()
  !     rmat_stripes = img_stripes%get_rmat()
  !     do n_stripe = 1, int(maxval(rmat_stripes))
  !         pos(:3) = minloc(abs(rmat_stripes-real(n_stripe)))
  !         label = nint(rmat_cc(pos(1),pos(2),pos(3)))
  !     enddo
  !
  ! end subroutine full_atom_stripes

  !This subrotuine indentifies the 'stripes' of atoms in the z direction.
 ! The idea is to look at the projection of the nanoparticle on the
 ! xy plane and identify all the atoms whose center has (almost, see
 ! MAX_DIST_CENTERS) the same x and y coords.
 ! The inputs are:
 ! -) centers, coordinates of the centers of mass of the atoms;
 ! -) img_stripes, an empty (but ecisting) vol in which all the
 !    stripes will be saved
 ! For how it is built, centers(:3,i) contains the coords of the i-th
 ! cc.
  subroutine center_atom_stripes(img_cc,img_stripes,centers) !change the name
        type(image), intent(in)    :: img_cc
        type(image), intent(inout) :: img_stripes
        real,        intent(in)    :: centers(:,:)
        integer :: i, j, ldim(3)
        integer, allocatable :: imat(:,:,:)
        logical, allocatable :: flag(:) !not to check more than once the same center
        integer :: cnt, cnt_additional  !number of distinguished stripes
        real, parameter :: MAX_DIST_CENTERS = 2.5
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real, allocatable :: rmat_stripes(:,:,:), rmat_cc(:,:,:)
        integer :: label
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(.not. img_stripes%exists()) then
             print *, 'Input image has to be constructed before! center_atom_stripes'
             stop
         endif
        ldim = img_stripes%get_ldim()
        allocate(flag(size(centers, dim = 2)),  source = .true. )
        allocate(imat        (ldim(1),ldim(2),ldim(3)), source = 0 )
        allocate(rmat_stripes(ldim(1),ldim(2),ldim(3)), source = 0.)
        rmat_cc = img_cc%get_rmat()
        cnt = 0
        cnt_additional = 0
        do i = 1, size(centers, dim = 2)
            if(flag(i)) cnt = cnt + 1
            do j = 1, size(centers, dim = 2)
                if(j > i .and. flag(i) .and. flag(j)) then
                    if(sqrt((centers(1,i)-centers(1,j))**2.+(centers(2,i)-centers(2,j))**2.) < MAX_DIST_CENTERS) then
                       if(DEBUGG )print *, 'i = ', i, '& j = ', j, 'PARALLEL,&
                                  &  tot dist = ', &
                                  & sqrt((centers(1,i)-centers(1,j))**2.+(centers(2,i)-centers(2,j))**2.)
                       imat(nint(centers(1,i)),&
                            & nint(centers(2,i)),nint(centers(3,i))) = cnt
                       imat(nint(centers(1,j)),&
                            & nint(centers(2,j)),nint(centers(3,j))) = cnt
                       flag(j) = .false.
                       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                       ! Set rmat_stripes=cnt for the whole atom (not just the center)
                       where(abs(rmat_cc-rmat_cc(nint(centers(1,j)),nint(centers(2,j)),nint(centers(3,j))))<TINY)
                            rmat_stripes = cnt
                        endwhere
                    endif
                endif
            enddo
            label = nint(rmat_cc(nint(centers(1,i)),nint(centers(2,i)),nint(centers(3,i))))
            where(abs(rmat_cc-label)<TINY) rmat_stripes = cnt
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!SURF EACH STRIPE!!!!!!!!!!!!!
            ! where(imat .ne. cnt) imat = 0
            ! if(flag(i) .and. any(imat > 0)) then
            !     cnt_additional = cnt_additional + 1
            !     call img_stripes%set_rmat(real(imat))
            !     call img_stripes%write('Stripe'//int2str(cnt_additional)//'.mrc')
            ! endif
            where(abs(rmat_stripes - cnt) > TINY) rmat_stripes = 0.
            if(flag(i) .and. any(rmat_stripes > TINY)) then
                cnt_additional = cnt_additional + 1
                call img_stripes%set_rmat(rmat_stripes)
                call img_stripes%write('FULLStripe'//int2str(cnt_additional)//'.mrc')
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        enddo
        ! call img_stripes%set_rmat(real(imat))
        ! call img_stripes%write('StripeAll.mrc')
  end subroutine center_atom_stripes


  subroutine aspect_ratios_hist(img,img_cc)
    use gnufor2
    type(image), intent(inout):: img
    type(image), intent(inout):: img_cc
    type(image) :: img_bin, img_stripes
    integer :: ldim(3)
    integer :: i, j
    real    :: smpd
    real, allocatable :: ratios(:)
    real, allocatable :: rmat_cc(:,:,:)
    real, allocatable :: paral(:,:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real, allocatable :: rmat_t(:,:,:)
    integer :: cnt
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ldim = img_cc%get_ldim()
    smpd = img_cc%get_smpd()
    call img_stripes%new(ldim, smpd)
    if(DEBUGG) call img_cc%write('OrderedCC.mrc')
    call img_cc%elim_cc([5, 100000])
    call img_cc%bin(0.5) !binarize to re-calculate the cc, to have them in order again
    img_bin = img_cc
    call img_cc%kill
    call img_bin%find_connected_comps(img_cc) !cc re-calculation
    if(DEBUGG) call img_cc%write('PolishedOrderedCC.mrc')
    rmat_cc   = img_cc%get_rmat()
    allocate(ratios    (nint(maxval(rmat_cc))), source = 0.)
    allocate(centers (3,nint(maxval(rmat_cc))), source = 0.)     !global variable
    do i = 1, nint(maxval(rmat_cc))
        call calc_aspect_ratio(img, img_cc, i, ratios(i))
    enddo
    call hist(ratios, 20)
    call center_atom_stripes(img_cc,img_stripes,centers)
    deallocate(ratios, rmat_cc, centers)
  end subroutine aspect_ratios_hist

  !This function calculates the centers of mass of an
  !atom. It takes in input the image, its connected
  !component (cc) image and the label of the cc that
  !identifies the atom whose center of mass is to be
  !calculated. The result is stored in m.
  function atom_masscen(img, img_cc, label) result(m)
      type(image), intent(inout) :: img, img_cc
      integer,     intent(in)    :: label
      real                       :: m(3)  !mass center coords
      real,    allocatable :: rmat(:,:,:), rmat_cc(:,:,:)
      integer, allocatable :: imat(:,:,:)
      integer :: i, j, k, ldim(3)
      integer :: sz !sz of the cc in the label
      rmat    = img%get_rmat()
      rmat_cc = img_cc%get_rmat()
      ldim = img%get_ldim()
      where(     abs(rmat_cc-real(label)) > TINY) rmat = 0.
      sz = count(abs(rmat_cc-real(label)) < TINY)
      m = 0.
      do i = 1, ldim(1)
          do j = 1, ldim(2)
              do k = 1, ldim(3)
                  m = m + rmat(i,j,k)*[i,j,k]
              enddo
          enddo
      enddo
      m = m/real(sz)
      if(ldim(3) == 1) m(3) = 0.
  end function atom_masscen
end module simple_test_chiara_nanoparticles_mod

program simple_test_chiara_nanoparticles
  include 'simple_lib.f08'
  use gnufor2
  use simple_picker_chiara
  use simple_test_chiara_nanoparticles_mod
  type(image) :: img, img_cc
  integer :: n_vol, i
  real, allocatable :: rmat(:,:,:), rmat_cc(:,:,:)
  real :: m(3)
  !!!!!!!!!!!!!!!!!!!!!!!INITIALIZATIONS AND ALLOCATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call img%new    ([160,160,160], 0.358)
  call img_cc%new ([160,160,160], 0.358)
  !!!!!!!!!!!!!!!!!!NANOPARTICLES BINARIZATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! do n_vol = 1, 7
  !     call img%read('particle'//int2str(n_vol)//'.mrc')
  !     call bin_nanopart(img)
  !     call img%write(int2str(n_vol)//'BINparticle.mrc')
  ! enddo
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!ASPECT RATIOS CALCULATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do n_vol = 1, 1!7
    call img%read(int2str(n_vol)//'BINparticle.mrc')
    call img%find_connected_comps(img_cc)
    call aspect_ratios_hist(img,img_cc)
enddo
end program simple_test_chiara_nanoparticles
