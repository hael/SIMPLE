! Use simple_test_chiara_nanopart_new prg=dummy smpd=0.358
module simple_nanopart_new_mod
include 'simple_lib.f08'
use simple_image, only : image
use simple_picker_chiara, only: pixels_dist, get_pixel_pos
use simple_cmdline,       only: cmdline
use simple_parameters,    only: parameters


implicit none

public :: nanoparticle
private
#include "simple_local_flags.inc"

! module global constants
integer, parameter :: N_THRESH         = 20      !number of thresholds for binarization
real, parameter    :: MAX_DIST_CENTERS = 3.1     !to determin atoms belonging to the same row
logical, parameter :: DEBUGG           = .false. !for debuggin purposes
! module global variables
real,    allocatable :: centers(:,:)
integer, allocatable :: loc_longest_dist(:,:)!for indentific of the vxl that determins the longest dim of the atom

type :: nanoparticle
    private
    ! these image objects are part of the instance to avoid excessive memory re-allocations
    type(image) :: img, img_bin, img_cc, img_centers !INSERT ALSO CORRESPONDING MATRICES??
    type(image) :: img_rows !don't know if I need it.
    integer     :: ldim(3)
    real        :: smpd
    integer     :: nptcls
    type(parameters) :: p
    type(cmdline)    :: cline
    ! these strings are part of the instance for reporting purposes
    character(len=STDLEN) :: partname
  contains
    procedure ::  new => new_nanoparticle
    procedure ::  binarize => nanopart_binarization
    procedure, private :: calc_aspect_ratio
    procedure, private :: aspect_ratios_hist
    procedure, private :: atom_rows
    procedure, private :: run_nanoparticle_job
    procedure :: run
    procedure, private::  kill => kill_nanoparticle
end type nanoparticle

contains
    ! 3D line for regression for identification atom rows
    function line_4regression(p,n) result(r)
            real,    intent(in) :: p(:)
            integer, intent(in) :: n
            real :: r(n)
            real :: x, y, z
            x = p(1)
            y = p(2)
            z = p(3)
            r(1) = 1
            r(2) = x
            r(3) = y
            r(4) = z
    end function line_4regression

    subroutine new_nanoparticle(self)
        class(nanoparticle), intent(inout) :: self
        self%partname = 'particle1.mrc' !I set it for simplicity, it will have to be passed by cmdline
        call find_ldim_nptcls(self%partname, self%ldim, self%nptcls, self%smpd)
        call self%img%new         (self%ldim, self%p%smpd)
        call self%img_bin%new     (self%ldim, self%p%smpd)
        call self%img_cc%new      (self%ldim, self%p%smpd) !don't think I need to initializate it
        call self%img_centers%new (self%ldim, self%p%smpd)
        call self%img_rows%new    (self%ldim, self%p%smpd) !don't know if need it
    end subroutine new_nanoparticle

    !This subrotuine takes in input a nanoparticle and
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
    subroutine nanopart_binarization( self )
        class(nanoparticle), intent(inout) :: self
        real, allocatable :: rmat(:,:,:), rmat_t(:,:,:)
        real    :: step          !histogram disretization step
        real    :: thresh        !binarization threshold
        real    :: seleted_t(1)  !selected threshold for nanoparticle binarization
        real    :: x_thresh(N_THRESH-1), y_med(N_THRESH-1)
        integer, allocatable :: sz(:) !size of the ccs and correspondent label
        integer ::  i
        write(logfhandle,*) '****binarization, init'
        rmat = self%img%get_rmat()
        allocate(rmat_t(self%ldim(1),self%ldim(2),self%ldim(3)), source = 0.)
        step = maxval(rmat)/real(N_THRESH)
        do i = 1, N_THRESH-1
            thresh = real(i)*step
            where(rmat > thresh)
                rmat_t = 1.
            elsewhere
                rmat_t = 0.
            endwhere
            call self%img_bin%set_rmat(rmat_t)
            call self%img_bin%find_connected_comps(self%img_cc)
            sz          = self%img_cc%size_connected_comps()
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
        call self%img_bin%set_rmat(rmat_t)
        deallocate(rmat, rmat_t, sz)
        write(logfhandle,*) '****binarization, completed'
    end subroutine nanopart_binarization

    ! This subroutine takes in input a connected component (cc) image
    ! and the label of one of its ccs and calculates its aspect ratio, which
    ! is defined as the ratio of the width and the height.
    ! The idea behind this is that the center of the cc is calculated,
    ! than everything is deleted except the borders of the cc. Finally,
    ! in order to calculate the width and the height, the min/max
    ! distances between the center and the borders are calculated. The
    ! aspect ratio is the ratio of those 2 distances.
    subroutine calc_aspect_ratio(self, label, ratio)
        class(nanoparticle), intent(inout) :: self
        integer,             intent(in)    :: label
        real,                intent(out)   :: ratio
        integer, allocatable :: imat(:,:,:)
        integer, allocatable :: pos(:,:)
        real,    allocatable :: rmat(:,:,:), rmat_cc(:,:,:)
        logical, allocatable :: border(:,:,:)
        integer :: location(1) !location of the farest vxls of the atom from its center
        integer :: i
        real    :: shortest_dist, longest_dist
        rmat    = self%img%get_rmat()
        rmat_cc = self%img_cc%get_rmat()
        centers(:3,label) = atom_masscen(self,label)
        allocate(imat(self%ldim(1),self%ldim(2),self%ldim(3)),source = nint(rmat_cc)) !for function get_pixel_pos
        call self%img_cc%border_mask(border, label)
        where(border .eqv. .true.)
            imat = 1
        elsewhere
            imat = 0
        endwhere
        call get_pixel_pos(imat,pos)                       !pxls positions of the shell
        shortest_dist = pixels_dist(centers(:,label), real(pos),'min')
        longest_dist  = pixels_dist(centers(:,label), real(pos),'max',location)
        loc_longest_dist(:3, label) =  pos(:3,location(1))
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
            write(logfhandle,*) 'longest dist = ',   longest_dist
            write(logfhandle,*) 'RATIO = ', ratio
        endif
        deallocate(rmat, rmat_cc, border, imat, pos)
    contains

        !This function calculates the centers of mass of an
        !atom. It takes in input the image, its connected
        !component (cc) image and the label of the cc that
        !identifies the atom whose center of mass is to be
        !calculated. The result is stored in m.
        function atom_masscen(self, label) result(m)
            type(nanoparticle), intent(inout) :: self
            integer,            intent(in)    :: label
            real(sp)             :: m(3)  !mass center coords
            real,    allocatable :: rmat_in(:,:,:), rmat_cc_in(:,:,:)
            integer :: i, j, k
            integer :: sz !sz of the cc in the label
            rmat_in    = self%img_bin%get_rmat()
            rmat_cc_in = self%img_cc%get_rmat()
            where(     abs(rmat_cc_in-real(label)) > TINY) rmat_in = 0.
            sz = count(abs(rmat_cc_in-real(label)) < TINY)
            m = 0.
            !omp do collapse(3) private(i,j,k,m) schedule(static)
            do i = 1, self%ldim(1)
                do j = 1, self%ldim(2)
                    do k = 1, self%ldim(3)
                        if(abs(rmat_in(i,j,k))> TINY) m = m + rmat_in(i,j,k)*[i,j,k]
                    enddo
                enddo
            enddo
            !omp end do
            m = m/real(sz)
            if(self%ldim(3) == 1) m(3) = 0. !for 2D imgs
        end function atom_masscen
    end subroutine calc_aspect_ratio

    ! This subrotuine indentifies the 'rows' of atoms in the z direction.
    ! The idea is to look at the projection of the nanoparticle on the
    ! xy plane and identify all the atoms whose center has (almost, see
    ! MAX_DIST_CENTERS) the same x and y coords.
    ! The inputs are:
    ! -) centers, coordinates of the centers of mass of the atoms;
    ! For how it is built, centers(:3,i) contains the coords of the
    ! center of mass of the i-th cc.
    subroutine atom_rows(self,centers) !change the name
        use gnufor2
        class(nanoparticle), intent(inout) :: self
          real,             intent(in)    :: centers(:,:)
          integer :: i, j, h, k
          integer, allocatable :: imat(:,:,:)
          logical, allocatable :: flag(:)    !not to check more than once the same center
          integer :: cnt, cnt_additional
          real    :: label_i, label_j
          real, allocatable :: rmat_rows(:,:,:), rmat_cc(:,:,:)
          integer :: label_max, label1, label2
          integer :: m(3)
          real :: theta, mod_1, mod_2, dot_prod
          real, allocatable :: ang_var(:), cc_in_the_row(:)
          real, allocatable :: loc_ld_real(:,:), rmat_aux_cc(:,:,:)
          logical, allocatable :: do_it(:)   ! not to do it more than once, cnt sometimes is not increased
          type(image) :: img_aux, img_aux_cc !to store the cc of just one row of atoms
          integer :: n_row                   !number of atoms rows
          print *, 'atom rows analysis'
          call img_aux%new   (self%ldim, self%p%smpd)
          call img_aux_cc%new(self%ldim, self%p%smpd)
          allocate(flag(size(centers, dim = 2)),  source = .true. )
          allocate(imat     (self%ldim(1),self%ldim(2),self%ldim(3)), source = 0 )
          allocate(rmat_rows(self%ldim(1),self%ldim(2),self%ldim(3)), source = 0.)
          allocate(loc_ld_real(3,size(loc_longest_dist, dim=2)), source = 0. )
          rmat_cc = self%img_cc%get_rmat()
          cnt = 0
          cnt_additional = 0
          allocate(do_it(size(centers, dim = 2)), source = .true.)  !initialise
          do i = 1, size(centers, dim = 2)-1
              if(flag(i)) cnt = cnt + 1
              do j = 2, size(centers, dim = 2)
                  if(j > i .and. flag(i) .and. flag(j)) then
                      if(sqrt((centers(1,i)-centers(1,j))**2+(centers(2,i)-centers(2,j))**2) < MAX_DIST_CENTERS) then
                         !if(DEBUGG )
                         print *, 'i = ', i, '& j = ', j, 'PARALLEL,&
                                    &  tot dist = ', &
                                    & sqrt((centers(1,i)-centers(1,j))**2+(centers(2,i)-centers(2,j))**2+(centers(3,i)-centers(3,j))**2)
                         imat(nint(centers(1,i)),&
                              & nint(centers(2,i)),nint(centers(3,i))) = cnt
                         imat(nint(centers(1,j)),&
                              & nint(centers(2,j)),nint(centers(3,j))) = cnt
                         flag(j) = .false.
                         ! Set rmat_rows=cnt for the whole atom (not just the center)
                         label_j = rmat_cc(nint(centers(1,j)),nint(centers(2,j)),nint(centers(3,j)))
                         label_i = rmat_cc(nint(centers(1,i)),nint(centers(2,i)),nint(centers(3,i)))
                         if(abs(label_j) > TINY) then
                             where(abs(rmat_cc-label_j)<TINY)
                                  rmat_rows = cnt
                             endwhere
                         endif
                      endif
                  endif
              enddo
              if(abs(label_i)>TINY) then
                  where(abs(rmat_cc-label_i)<TINY) rmat_rows = real(cnt)
              endif
              !SURF EACH ROW!!!!!!!!!!!!!
              ! where(abs(rmat_rows - real(cnt)) > TINY) rmat_rows = 0.
              ! if(flag(i) .and. any(rmat_rows > TINY)) then
              !     cnt_additional = cnt_additional + 1
              !     call self%img_rows%set_rmat(rmat_rows)
              !     call self%img_rows%write('FULLRow'//int2str(cnt_additional)//'xy.mrc')
              ! endif
              !!!!!!!!!!!!!!!!!!!!!!!!!!!
          enddo
          print *, 'NEW PART'
          call self%img_rows%set_rmat(rmat_rows)
          call self%img_rows%write('ImgROWS.mrc')
          print *, 'Number of rows identified: ', cnt
          !take each atom row and analyse it for polarization
          do n_row = 1, cnt
              if(do_it(n_row)) then
                  do_it(n_row) = .false.
                  rmat_rows = self%img_rows%get_rmat()   ! restore rmat_rows
                  rmat_cc   = self%img_cc%get_rmat()     ! restore rmat_cc
                  where(abs(rmat_rows - real(n_row)) > TINY)
                      rmat_rows = 0.
                      rmat_cc = 0.
                  endwhere
                  call img_aux%set_rmat(rmat_cc)
                  call img_aux%find_connected_comps(img_aux_cc)
                  rmat_aux_cc = img_aux_cc%get_rmat()
                  cnt_additional = 0
                  label_max = nint(maxval(rmat_aux_cc)) ! cc which corresponds to the highest label in the row
                  !label_max has to be at least 4 to have significance (rows with enough number of atoms)
                  !print *, 'row= ', n_row, 'label_max = ', label_max
                  if(label_max > 3) then
                      allocate(ang_var      (label_max-1), source = 0.) !NOT SURE
                      allocate(cc_in_the_row(label_max-1), source = 0.) !NOT SURE
                      do k = 1, label_max-1
                          m = minloc(abs(rmat_aux_cc(:,:,:)-real(k))) !select one pxl labeled k
                          label1 = nint(rmat_cc(m(1),m(2),m(3)))      !identify the corresponding label in rmat_cc
                          ! calculate for k just the first time
                          if (k == 1) loc_longest_dist(:3,k) = loc_longest_dist(:3,k)-nint(centers(:3,label1))+1 !vxl identifying the longest dim of the atom translated into the origin
                          !MAYBE TRANSLATE IT WHEN U CALCULATE IT IN CALC_ASPECT RATIO
                          m = minloc(abs(rmat_aux_cc(:,:,:)-real(k+1))) !select one pxl labeled k+1
                          label2 = nint(rmat_cc(m(1),m(2),m(3)))
                          !translating initial point(centers) of the vector center-loc_longest_dist in the origin(1,1,1)
                          loc_longest_dist(:3,k+1) = loc_longest_dist(:3,k+1)-nint(centers(:3,label2))+1 !+1 because the origin is not (0,0,0)
                          ! https://www.intmath.com/vectors/7-vectors-in-3d-space.php FORMULAE HERE
                          loc_ld_real = real(loc_longest_dist) !translation in reals
                          dot_prod = loc_ld_real(1,k)*loc_ld_real(1,k+1)+ &
                          &          loc_ld_real(2,k)*loc_ld_real(2,k+1)+ &
                          &          loc_ld_real(3,k)*loc_ld_real(3,k+1)
                          mod_1 = sqrt(loc_ld_real(1,k  )**2+loc_ld_real(2,k  )**2+loc_ld_real(3,k  )**2)
                          mod_2 = sqrt(loc_ld_real(1,k+1)**2+loc_ld_real(2,k+1)**2+loc_ld_real(3,k+1)**2)
                          if(dot_prod/(mod_1*mod_2) > 1. .or. dot_prod/(mod_1*mod_2)< -1.) then
                              THROW_WARN('Out of the domain of definition of arccos; atom_rows')
                              theta = 0.
                          else
                              theta = acos(dot_prod/(mod_1*mod_2)) !the output of acos in RADIANS
                          endif
                          if(DEBUGG) then
                              print *, '>>>>>>>>>>>>>>>>>>>>>>'
                              print *, ' K = ', k
                              print *, 'dot_prod = ', dot_prod
                              print *, 'mod_1 = ', mod_1
                              print *, 'mod_2 = ', mod_2
                              print *, 'mod1*mod2=',(mod_1*mod_2), 'frac = ', dot_prod/(mod_1*mod_2)
                              print *, 'in radians theta =: ', theta, 'degrees= ', rad2deg(theta)
                          endif
                          ang_var(k) = rad2deg(theta) !radians
                          cc_in_the_row(k) = real(k)
                      enddo
                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      if(DEBUGG) then
                          print *, '>>>>>>>>>>>>>>>>>>>>>>'
                          print *, 'N_ROW = ', n_row
                          print *, 'cc_in_the_row=', cc_in_the_row
                          print *, 'ang_var=', ang_var
                      endif
                      open(129, file='Polarization_xy',position = 'append')
                      write (129,*) 'Ang_vars'//int2str(n_row)//'=[...'
                      do k = 1, label_max-1
                              write (129,'(A)', advance='no') trim(real2str(cc_in_the_row(k)))
                              write (129,'(A)', advance='no') ', '
                              write (129,'(A)', advance='no') trim(real2str(ang_var(k)))
                              write (129,'(A)')'; ...'
                      end do
                      write (129,*) '];'
                      close(129)
                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      ! call plot(cc_in_the_row, ang_var, xlabel='Connectec components of row'//int2str(n_row), ylabel='AngVar')
                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      deallocate(ang_var,cc_in_the_row)
                  endif
              endif
          enddo
          print *, 'atom rows analysis completed'
    end subroutine atom_rows

    subroutine aspect_ratios_hist(self)
      use gnufor2
      class(nanoparticle), intent(inout) :: self
      integer :: label
      real, allocatable :: ratios(:)
      real, allocatable :: rmat_cc(:,:,:)
      integer :: cnt
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !CENTERS PLOTTING
      ! type(image) :: img_centers
      ! real, allocatable :: rmat_centers(:,:,:)
      !DEBUGGING
      ! real :: dist
      ! real, allocatable :: other_centers(:,:)
      ! integer :: i
      !SVD FITTING
      !     real,  allocatable :: sig(:), y(:)
      !     real :: a(4), w(4)
      !     real :: v(4,4)
      !     real :: chisq
      ! !LINE PLOTTING
      ! type(image) :: img_line
      ! real, allocatable :: rmat_line(:,:,:), rmat(:,:,:)
      ! integer :: i, j, k
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call self%img_bin%find_connected_comps(self%img_cc)
      call self%img_cc%elim_cc([5, 1000]) !connected components clean up
      call self%img_cc%bin(0.5) !binarize to re-calculate the cc, to have them in order again
      self%img_bin = self%img_cc
      call self%img_cc%kill
      call self%img_bin%find_connected_comps(self%img_cc) !cc re-calculation
      if(DEBUGG) call self%img_cc%write('PolishedOrderedCC.mrc')
      rmat_cc   = self%img_cc%get_rmat()
      allocate(ratios          (  nint(maxval(rmat_cc))), source = 0.)
      allocate(centers         (3,nint(maxval(rmat_cc))), source = 0.)     !global variable
      allocate(loc_longest_dist(3,nint(maxval(rmat_cc))), source = 0 )     !global variable
      print *, 'aspect ratios calculations'
      do label = 1, nint(maxval(rmat_cc))
          call calc_aspect_ratio(self, label, ratios(label))
      enddo
      print *, 'aspect ratios calculations completed'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! allocate(other_centers(3,size(centers, dim = 2)-1), source = 0.)
      ! do label = 1, size(centers, dim = 2)
      !     cnt = 0
      !     do i = 1, size(centers, dim = 2)
      !         if(i .ne. label) then
      !             cnt = cnt + 1
      !             other_centers(:,cnt) = centers(:,i)
      !         endif
      !     enddo
      !     dist = pixels_dist(centers(:,label),other_centers(:,:),'min')
      !     print *, 'label = ', label, 'min_dist = ', dist
      ! enddo
      ! stop
      ! allocate(sig(size(centers, dim = 2)), source = 1.)
      ! allocate(y  (size(centers, dim = 2)), source = 0.) ! y is the grayvalue of the nanopart in the centers coords
      !                                                    !i suppose to do this on the binary version of the volume
      ! rmat = self%img%get_rmat()
      ! call img_line%new(self%ldim,self%p%smpd)
      ! call img_line%set_rmat(rmat)
      ! !!!!!!!!!!!!!!!!!!!!!! TO START AGAIN FROM HERE
      ! print *, 'centers = '
      ! call vis_mat(centers)
      ! do i = 1, size(centers, dim = 2)
      !         if( (centers(1,i)>0.) .and. (centers(1,i)<real(self%ldim(1))) .and. &
      !         &   (centers(2,i)>0.) .and. (centers(2,i)<real(self%ldim(2))) .and. &
      !         &   (centers(3,i)>0.) .and. (centers(3,i)<real(self%ldim(3))) ) then
      !         y(i) = rmat( nint(centers(1,i)),  nint(centers(2,i)), nint(centers(3,i)))
      !     else
      !         y(i) = 0.
      !     endif
      ! enddo
      ! The line I am looking for is called the 3D Orthogonal Distance Regression (ODR) line
      ! call svd_multifit(centers,y,sig,a,v,w,chisq,line_4regression)
      ! print *, 'a = ', a
      ! print *, '>>>>>>>>>>>> TESTING >>>>>>>>>>>>'
      ! do i = 1, size(centers, dim = 2)
      !     print *, 'for center', i, a(1)+a(2)*centers(1,i)+a(3)*centers(2,i)+a(4)*centers(3,i)-1. , 'IT SHOULD BE 0 if approximated'
      ! enddo
      ! rmat_line = img_line%get_rmat()
      ! rmat_line = 0.
      ! rmat = self%img%get_rmat()
      ! !omp do collapse(3) private(i,j,k) schedule(static)
      ! do i = 1, self%ldim(1)
      !     do j = 1, self%ldim(2)
      !         do k = 1, self%ldim(3)
      !             if(abs(a(1)+a(2)*i+a(3)*j+a(4)*k-rmat(i,j,k)) < TINY) rmat_line = 1.
      !         enddo
      !     enddo
      ! enddo
      ! !omp end do
      ! call img_line%set_rmat(rmat_line)
      ! call img_line%write('ImgLine.mrc')

      !stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!CENTERS PLOT 3D!!!!!!!!!!!!!!!!!!!!!!!
      ! call img_centers%new(self%ldim, self%smpd)
      ! rmat_centers = img_centers%get_rmat()
      ! do label = 1, size(centers, dim = 2)
      !     rmat_centers(nint(centers(1,label)), nint(centers(2,label)), nint(centers(3,label))) = 1.
      ! enddo
      ! call img_centers%set_rmat(rmat_centers)
      ! call img_centers%write('CENTERS.mrc')
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      cnt = 0
      open(125, file='AspectRatios', position='append')
      write (125,*) 'Ratios=[...'
      do cnt = 1, size(ratios)
          write (125,'(A)', advance='no') trim(real2str(ratios(cnt)))
          if(cnt < size(ratios)) write (125,'(A)', advance='no') ', '
      end do
      write (125,*) '];'
      close(125)
      !call hist(ratios, N_THRESH)
      call atom_rows(self,centers)
      deallocate(ratios, rmat_cc, centers)
    end subroutine aspect_ratios_hist

    subroutine run_nanoparticle_job(self, n_vol)
      class(nanoparticle), intent(inout) :: self
      integer,             intent(in)    :: n_vol !how many nanoparticle4s you want to analyse
      integer ::  i
      real, allocatable :: rmat(:,:,:), rmat_cc(:,:,:)
      real    :: m(3)
      !!!!!!!!!!!!!!!!!!NANOPARTICLES BINARIZATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i = 1, n_vol
        call self%img%read(self%partname)
        call self%binarize()
        call self%img_bin%write(int2str(i)//'BINparticle.mrc')
      enddo
      !!!!!!!!!!!!!!!!!!!!!!!ASPECT RATIOS CALCULATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 1) aspect ratio
      ! 2) polarization
       do i= 1, n_vol
        call aspect_ratios_hist(self)
       enddo
    end subroutine run_nanoparticle_job

    subroutine run(self)
        class(nanoparticle), intent(inout) :: self
        call self%cline%parse_oldschool()
        call self%cline%checkvar('smpd', 1)
        call self%cline%check()
        call self%p%new(self%cline)
        call self%new()
        call self%run_nanoparticle_job(1) !just one nanoparticle
        call self%kill()
    end subroutine run

    subroutine kill_nanoparticle(self)
        class(nanoparticle), intent(inout) :: self
        call self%img%kill()
        call self%img_bin%kill()
        call self%img_cc%kill()
        call self%img_centers%kill()
        call self%img_rows%kill()
    end subroutine kill_nanoparticle
end module simple_nanopart_new_mod

program simple_test_chiara_nanopart_new
    use simple_nanopart_new_mod
    implicit none
    type(nanoparticle) :: nanopart
    call nanopart%run()
end program simple_test_chiara_nanopart_new
