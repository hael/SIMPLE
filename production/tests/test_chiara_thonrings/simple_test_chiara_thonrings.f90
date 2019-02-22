! ! Use simple_test_chiara_thonrings prg=dummy smpd=1.14
module simple_test_chiara_thonrings_mod
include 'simple_lib.f08'
use simple_image,      only : image
use simple_stackops,   only : binarize_stack
use simple_cmdline,    only : cmdline
use simple_parameters, only : parameters

implicit none

public :: powerspectrum
private
#include "simple_local_flags.inc"

! module global variables
real, allocatable :: res_vec(:,:)
integer           :: nn_shells

type :: powerspectrum
    private
    ! these image objects are part of the instance to avoid excessive memory re-allocations
    type(image), allocatable :: img(:), img_bin(:)
    integer     :: ldim(3)
    real        :: smpd
    integer     :: nstack
    type(parameters) :: p
    type(cmdline)    :: cline
    ! these strings are part of the instance for reporting purposes
    character(len=STDLEN) :: powerspec_name
  contains
    procedure ::  new => new_powerspectrum
    procedure, private :: delete_shells
    procedure, private :: extract_3points
    procedure, private :: identify_center_radius
    procedure, private :: build_ice_template
    procedure, private :: is_symmetric
    procedure, private :: build_resolutions_vector
    procedure, private :: find_res
    procedure, private :: discard_ps
    procedure, private :: ice_presence
    procedure, private :: process_ps_stack
    procedure, private :: run_powerspectrum_job
    procedure :: run
    procedure, private::  kill => kill_powerspectrum
end type powerspectrum

contains

    subroutine new_powerspectrum(self)
        class(powerspectrum), intent(inout) :: self
        integer :: i
        ! self%powerspec_name = 'pspecs_sphire_tstdat.mrc' !I set it for simplicity, it will have to be passed by cmdline
        self%powerspec_name = 'pspecs_saga_polii.mrc'
        call find_ldim_nptcls(self%powerspec_name, self%ldim, self%nstack, self%smpd)
        self%ldim(3) = 1 !do not consider it as a stack
        allocate(self%img    (self%nstack))
        allocate(self%img_bin(self%nstack))
        do i = 1, self%nstack
            call self%img(i)%new     (self%ldim, self%smpd)
            call self%img_bin(i)%new (self%ldim, self%smpd)
        enddo
        do i = 1, self%nstack
            call self%img(i)%read    (self%powerspec_name,i)
        enddo
    end subroutine new_powerspectrum

    !This subroutine is meant to count the nb of visible and contiguous
    !rings are detected in the binarization of a power-spectra stack.
    !It saves the nb of vis rings in the i-th image of the stack
    !in the variable counts(i).
    subroutine count_visible_continuous_ring(self,n,lim,count)
        class(powerspectrum), intent(inout) :: self
        integer,              intent(in)    :: n
        integer,              intent(in)    :: lim  !pixel until the one the variations have to be checked
        integer,              intent(out)   :: count
        real, allocatable :: rmat_bin(:,:,:)
        rmat_bin = self%img_bin(n)%get_rmat()
        call count_grad_variations(rmat_bin,lim,count)
    contains
    !This subroutine takes in input the matrix associated
    !with a binary image and calculates the number of variations
    !white -> black or black -> white in the direction of the x-axis,
    !starting from te central pixel of the image. The result is
    !store in count.
    !It is meant to calculate how many visible and contiguous
    !rings are identified in the binarization of a power spectrum.
    subroutine count_grad_variations(rmat_bin,lim,count)
        real,    intent(in)  :: rmat_bin(:,:,:)
        integer, intent(in)  :: lim  !pixel until the one the variations have to be checked
        integer, intent(out) :: count
        integer :: s(3) !shape of the matrix
        integer :: j
        integer, allocatable :: imat(:,:,:) !input is meant to be binary, so I consider the integer matrix associated
        integer :: tmp_val
        integer :: from,to,fixed
        if(any(rmat_bin <-0.0001) .or. any(rmat_bin>1.0001)) THROW_HARD('Non implemented non binary images; count_grad_variations')
        s = shape(rmat_bin)
        allocate(imat(s(1),s(2),s(3)), source = nint(rmat_bin))
        if(s(3) .ne. 1) THROW_HARD('Non implemented for volumes; count_grad_variations')
        count   = 0               !counter for gradient variations
        fixed   = s(1)/2
        from    = s(2)/2
        to      = lim !from + lim
        tmp_val = imat(fixed,from-1,1) !initialize
        do j = from, to
            if(imat(fixed,j,1) .ne. tmp_val) then
                count = count + 1
                tmp_val = imat(fixed,j,1)
            endif
        enddo
        count = count/2
    end subroutine count_grad_variations
    end subroutine count_visible_continuous_ring

    !This function takes in input an image (suppsed to be
    !a connected components image) and deletes (sets to 0)
    !the internal shells (up tp 5A resolution) and the
    !external frame (up to 3. A).
    !It is meant to discard the connected components which
    !are not in the range to represent ice contaminations.
    subroutine delete_shells(self, n, box)
        class(powerspectrum), intent(inout) :: self
        integer,              intent(in)    :: n
        integer,     intent(in)    :: box
        integer :: x, i, j
        real, allocatable :: rmat(:,:,:)
        x = calc_fourier_index(5.,box,self%smpd)
        rmat = self%img_bin(n)%get_rmat()
        !Discard internal shell
        do i = 1, self%ldim(1)
            do j = 1, self%ldim(2)
                if((real(i)-real(box/2))**2/(real(x)**2) + (real(j)-real(box/2))**2/(real(x)**2) - 1 < TINY) rmat(i,j,1) = 0.
            enddo
        enddo
        !Discard external frame
        x = calc_fourier_index(3.,box,self%smpd)
        rmat(:box/2-x,:,1) = 0.
        rmat(:,:box/2-x,1) = 0.
        rmat(box/2+x:,:,1) = 0.
        rmat(:,box/2+x:,1) = 0.
        call self%img_bin(n)%set_rmat(rmat)
        call self%img_bin(n)%write('Deleted.mrc')
    end subroutine delete_shells

    !This subroutine identifies the location of 3 of the
    !white pixels belonging to the cc label and stores
    !them in the variable points.
    subroutine extract_3points(self, n, img_cc, label, points)
        class(powerspectrum), intent(inout) :: self
        type(image),          intent(inout) :: img_cc
        integer,              intent(in)    :: n
        integer,              intent(in)    :: label
        integer,              intent(out)   :: points(3,2)
        real, allocatable :: rmat(:,:,:), rmat_2d(:,:)
        integer :: i, j, cnt
        integer, allocatable :: pos_white_pxls(:,:)
        rmat = img_cc%get_rmat()
        where(abs(rmat-real(label)) > TINY) rmat = 0. !keep just the cc identified by label
        allocate(rmat_2d(self%ldim(1), self%ldim(2)), source = 0.)
        rmat_2d(:,:) = rmat(:,:,1)
        allocate(pos_white_pxls(count(rmat > 0.5), 2), source = 0)
        cnt = 0
        do j = 1, self%ldim(2)
            do i = 1, self%ldim(1)
                if(rmat_2d(i,j) > 0.5) then
                    cnt = cnt + 1
                    pos_white_pxls(cnt, :) = [i,j] !save the location of all the white pixels of the cc label
                endif
            enddo
        enddo
        !Extrat 3 of the identified point
        points(1,:) = pos_white_pxls(1,:)
        points(2,:) = pos_white_pxls(cnt/2,:)
        points(3,:) = pos_white_pxls(cnt,:)
        deallocate(rmat, rmat_2d, pos_white_pxls)
    end subroutine extract_3points

    ! subroutine extract_3points(self, n, label, points)
    !     class(powerspectrum), intent(inout) :: self
    !     integer,              intent(in)    :: n
    !     integer,              intent(in)    :: label
    !     integer,              intent(out)   :: points(3,2)
    !     real, allocatable :: rmat(:,:,:), rmat_2d(:,:)
    !     integer :: i
    !     rmat = self%img(n)%get_rmat()
    !     where(abs(rmat-real(label)) > TINY) rmat = 0. !keep just the cc identified by label
    !     allocate(rmat_2d(self%ldim(1), self%ldim(2)), source = 0.)
    !     rmat_2d(:,:) = rmat(:,:,1)
    !     do i = 1, 3
    !         points(i,:) = minloc(abs(rmat_2d(:,:)-real(label)))
    !         if(points(i,1)-1 > 0 .and. points(i,1)+1< self%ldim(1)) then
    !              rmat_2d(points(i,1)-1:points(i,1)+1,points(i,2)) = 0.
    !              rmat_2d(points(i,1)-1:points(i,1)+1,points(i,2)) = 0.
    !         else
    !              rmat_2d(points(i,1),points(i,2)) = 0.
    !         endif
    !         if(points(i,2)-1 > 0 .and. points(i,2)+1<  self%ldim(2)) then
    !              rmat_2d(points(i,1),points(i,2)-1:points(i,2)+1) = 0.
    !              rmat_2d(points(i,1),points(i,2)-1:points(i,2)+1) = 0.
    !         else
    !              rmat_2d(points(i,1),points(i,2)) = 0.
    !         endif
    !     enddo
    !     if(any(points(:,:) .eq. 1))then
    !         rmat_2d(:,:) = rmat(:,:,1) !restore
    !         do i =1, 3
    !             points(i,:) = minloc(abs(rmat_2d(:,:)-real(label)))
    !             rmat_2d(points(i,1),points(i,2)) = 0.
    !         enddo
    !     endif
    !     deallocate(rmat, rmat_2d)
    ! end subroutine extract_3points

    !This subroutine takes in input the coordinates of 3 points of
    !a cc and identifies the centers coords and the magnitude of
    !the radius of the circumference passing for the 3points.
    !Geometric forumula.
    subroutine identify_center_radius(self, points, center, radius, bad)
        class(powerspectrum), intent(inout) :: self
        integer,              intent(in)    :: points(3,2)
        integer,              intent(out)   :: center(2)
        real,                 intent(out)   :: radius
        logical,              intent(out)   :: bad
        real :: x, y
        bad = .false.
        if(   (2*(points(1,1)*(points(2,2)-points(3,2)) - &
            &     points(1,2)*(points(2,1)-points(3,1)) + &
            &     points(2,1)*points(3,2)-points(3,1)*points(2,2))) &
            &      .ne. 0) then
            center(1) = ((points(1,1)**2+points(1,2)**2)*(points(2,2)-points(3,2)) + &
            &            (points(2,1)**2+points(2,2)**2)*(points(3,2)-points(1,2)) + &
            &            (points(3,1)**2+points(3,2)**2)*(points(1,2)-points(2,2)))/ &
            & (2*(points(1,1)*(points(2,2)-points(3,2)) - &
            &     points(1,2)*(points(2,1)-points(3,1)) + &
            &     points(2,1)*points(3,2)-points(3,1)*points(2,2)))

            center(2) = ((points(1,1)**2+points(1,2)**2)*(points(3,1)-points(2,1)) + &
            &            (points(2,1)**2+points(2,2)**2)*(points(1,1)-points(3,1)) + &
            &            (points(3,1)**2+points(3,2)**2)*(points(2,1)-points(1,1)))/ &
            & (2*(points(1,1)*(points(2,2)-points(3,2)) - &
            &     points(1,2)*(points(2,1)-points(3,1)) + &
            &     points(2,1)*points(3,2)-points(3,1)*points(2,2)))
        else
            center(:) = points(1,:)
            bad = .true.
            return
        endif
        x = real(center(1))
        y = real(center(2))
        radius = sqrt((x-real(points(1,1)))**2+(y-real(points(1,2)))**2)
         write(logfhandle,*)'CENTER = ', center, 'RADIUS = ', radius
    end subroutine identify_center_radius

    !This function  takes in input the box sz and the smpd
    !of a power spectra image and consequentially builds an
    !ice template (circular)
    function build_ice_template(self, box, radius) result(img_templ)
        class(powerspectrum), intent(inout) :: self
        integer, intent(in) :: box
        real,    intent(in) :: radius
        type(image) :: img_templ
        call img_templ%new([box,box,1], self%smpd)
        call img_templ%ellipse([box/2,box/2], [radius,radius], 'yes')
    end function build_ice_template

    ! This function takes in input a connected component (cc) image
    ! and the label of a cc. If this cc is in the III or IV quadrant,
    ! then it discard is true and it doesn't do anything.
    ! Otherwise it checks whether the cc is symmetric wrt
    ! the origin and stores the answer in yes_no.
    function is_symmetric(self, n, img_cc, label, discard, label_mirror) result(yes_no)
        class(powerspectrum), intent(inout) :: self
        integer,              intent(in)    :: n       !number of image in the stack
        type(image),          intent(inout) :: img_cc
        integer,              intent(in)    :: label   !label of the cc has to be checked
        logical,              intent(out)   :: discard
        real, optional,       intent(out)   :: label_mirror !correspondent symmetric label (CHANGE TO INTEGER??)
        type(image) :: img1, img2
        logical     :: yes_no
        integer     :: sz_cc, center(2), i, j, h, k, sh
        integer, allocatable :: m(:)
        real,    allocatable :: rmat(:,:,:), rmat_t(:,:,:)
        logical, allocatable :: mask(:,:,:)
        real, parameter      :: THRESH = 0.05
        real                 :: r
        rmat = img_cc%get_rmat()
        center(1) = nint(real(self%ldim(1))/2.)
        center(2) = nint(real(self%ldim(2))/2.)
        yes_no  = .false. !initialize
        discard = .false.
        allocate(rmat_t(self%ldim(1),self%ldim(2),self%ldim(3)), source = 0.)
        m = minloc(abs(rmat - real(label)))
        if(m(1) < center(1) .and. m(2) > center(2)) then ! II quandrant
            do i = 1, center(1)
                do j = center(2), self%ldim(2)
                    if(abs(rmat(i,j,1) - real(label)) < TINY) then
                        rmat_t(i,j,1) = rmat(2*center(1)-i,2*center(2)-j,1)
                        if(rmat(self%ldim(1)-i,self%ldim(2)-j,1)> 0.5) label_mirror = rmat_t(i,j,1)
                        where(rmat_t > .5) rmat_t = real(label)
                    endif
                enddo
            enddo
         else if(m(1) > center(1) .and. m(2) > center(2)) then  ! I quadrant
            do i = center(1), self%ldim(1)
                do j = center(2), self%ldim(2)
                    if(abs(rmat(i,j,1) - real(label)) < TINY) then
                         rmat_t(i,j,1) = rmat(self%ldim(1)-i,self%ldim(2)-j,1)
                         if(rmat(self%ldim(1)-i,self%ldim(2)-j,1)> 0.5) label_mirror = rmat_t(i,j,1)
                         where(rmat_t > .5) rmat_t = real(label)
                    endif
                enddo
            enddo
        else
            discard = .true.
        endif
        if(.not. discard) then
            call img1%new(self%ldim, self%smpd)
            call img2%new(self%ldim, self%smpd)
            where(abs(rmat - real(label)) > TINY) rmat = 0.
            call img1%set_rmat(rmat)
            call img2%set_rmat(rmat_t)
            r =  img1%real_corr(img2)
            !write(logfhandle,*)'cc ', label, 'correaltion ', r
            if(abs(r-1.) < THRESH) yes_no = .true.
            call img1%kill
            call img2%kill
        endif
        deallocate(rmat,rmat_t,m)
    end function is_symmetric

    !This function builds a vector that splits up the images in n_shells shells
    !and stores the result in res_vec(:,2). In res_vec(:,1) will be saved the
    !number of white pixels in the shell corresponding to res_vec(:,2).
    subroutine build_resolutions_vector(self, box, res_vec, n_shells)
        class(powerspectrum), intent(in)  :: self
        integer,              intent(in)  :: box
        real, allocatable,    intent(out) :: res_vec(:,:)
        integer, optional :: n_shells
        integer :: step, nn_shells, i
        nn_shells = 10
        if(present(n_shells) .and. n_shells < 2) THROW_HARD('Too low number of shells; build_resolutions_res_vec')
        if(present(n_shells)) nn_shells = n_shells
        if(allocated(res_vec)) deallocate(res_vec)
        allocate(res_vec(nn_shells+1,2), source = 0.)
        step = box/(2*n_shells)
        do i = 2, nn_shells+1
            res_vec(i,2) = (i-1)*real(step)
        end do
    end subroutine build_resolutions_vector

    !This function uses the ratio (white pxls)/(black pxls) per shell
    !to estimate how far do the detected rings go.
    !The resolution is set to the shell for which the nb of white
    !pxls is less than half of the avg nb of white pxls per shell.
    function find_res(self,n_image, box) result(res)
        use gnufor2
      class(powerspectrum), intent(in) :: self
      integer,              intent(in) :: n_image
      integer,              intent(in) :: box
      integer :: i, cnt
      real, allocatable :: rmat_bin(:,:,:)
      real    :: res
      real, allocatable :: x_hist(:)
      if(.not. allocated(res_vec)) THROW_HARD ('You have to build the resolution vector first; find_res')
      rmat_bin = self%img_bin(n_image)%get_rmat()
      if(any(rmat_bin>1.001) .or. any(rmat_bin<-0.001)) THROW_HARD('Img has to be binarised first; find_res')!sanity check
      !White pixels histogram
      allocate(x_hist(size(res_vec, dim = 1)), source = 0.)
      do i = 1, size(res_vec, dim = 1)
          x_hist(i) = res_vec(i,1)
      enddo
      ! open(125, file='HistogramWhitePxls', position='append')
      ! write (125,*) 'ratios'//int2str(n_image)//'=[...'
      ! do i = 1, size(res_vec, dim = 1)
      !     write (125,'(A)', advance='no') trim(real2str(x_hist(i)))
      !     if(i < size(res_vec, dim = 1)) write (125,'(A)', advance='no') ', '
      ! end do
      ! write (125,*) '];'
      ! close(125, status = 'keep')
      do i = 1, size(res_vec,dim=1)
          if(res_vec(i,1) <  sum(x_hist(:))/(2*size(x_hist))) then
              res = calc_fourier_index(real(res_vec(i,2)), box, self%smpd)
              return
          endif
      enddo
    end function find_res

    !This subroutine is meant to discard empty power spectra images.
    !If after binarisation the # of white pixels detected in the central
    !zone of the image is less than 5% of the tot
    !# of central pixels, than it returns yes, otherwhise no.
    function discard_ps(self, n) result(yes_no)
        class(powerspectrum), intent(inout) :: self
        integer, intent(in) :: n
        logical :: yes_no
        real, allocatable :: rmat(:,:,:),rmat_central(:,:,:)
        yes_no = .false.
        rmat = self%img_bin(n)%get_rmat()
        if(any(rmat > 1.001) .or. any(rmat < 0.)) THROW_HARD('Expected binary image in input; discard_ps')
        allocate(rmat_central(self%ldim(1)/2+1,self%ldim(2)/2+1,1), source = 0.)
        rmat_central(1:self%ldim(1)/2+1,1:self%ldim(2)/2+1,1) = &
        &        rmat( self%ldim(1)/2-self%ldim(1)/4 : self%ldim(1)/2+self%ldim(1)/4 , &
                                    & self%ldim(2)/2-self%ldim(2)/4 : self%ldim(2)/2+self%ldim(2)/4 , 1)
        deallocate(rmat)
        if(count(rmat_central(:,:,:) > 0.5)< 2*self%ldim(1)*self%ldim(2)*self%ldim(3)/(2*2*100)) yes_no = .true.; return
        deallocate(rmat_central)
    end function discard_ps

   !This subroutine is meant to identify ice presence, given the index
   !of the power spectra image in the stack, it prints if it detects
   !ice.
    subroutine ice_presence(self, n)
        class(powerspectrum), intent(inout) :: self
        integer,              intent(in)    :: n
        type(image)       :: img_cc
        real, allocatable :: rmat(:,:,:)
        integer :: j
        logical :: discard, yes_no
        call self%img_bin(n)%find_connected_comps(img_cc)
        rmat = img_cc%get_rmat()
        do j = 1, int(maxval(rmat))
        ! do j = int(minval(rmat,rmat > 0.5)), int(maxval(rmat))
              yes_no = is_symmetric(self, n, img_cc, j, discard)
              if(.not. discard) then
                  write(logfhandle,*)'cc n ', j, 'is symm ', yes_no
                  if(yes_no) write(logfhandle,*)'DETECTED ICE'
              endif
        enddo
    end subroutine ice_presence

    !This function takes in input the name of a stack of power spectra images (fname2process),
    !the name of a stack in which store the results (fname), the smpd and a low-pass parameter.
    !It binarises all the images in the stack and estimates how far do the detected rings go.
    subroutine process_ps_stack(self, fname, lp, n_shells)
      use gnufor2
      use simple_stackops, only : prepare_stack
      class(powerspectrum), intent(inout) :: self
      character(len=*),     intent(in) :: fname
      real,                 intent(in) :: lp
      integer, optional,    intent(in) :: n_shells  !number of shells
      integer, allocatable :: counts(:)             !nb of visible and contiguous rings
      integer, allocatable :: counter(:)
      real, allocatable :: rmat(:,:,:)
      integer           :: box
      integer           :: sh, ind(2), nn_shells
      integer           :: h, k, i, j, n_image
      logical, allocatable :: mask(:,:)
      real                 :: res, ax
      type(image)          :: img
      character(len = 100) :: iom
      integer              :: status
      integer              :: limit    !to count nb of gradient variations until there
      logical              :: discard
      box = self%ldim(1)
      nn_shells = 10
      if(present(n_shells)) nn_shells = n_shells
      !wwinsz = 2 !default
      !if(present(winsz)) wwinsz = winsz
      call prepare_stack(self%powerspec_name, 'prepared_stack.mrc', self%smpd, lp)
      write(logfhandle,*) '>>>>>>>>>>>>>STACK PREPARED SUCCESSFULLY>>>>>>>>>>>>>'
      call binarize_stack('prepared_stack.mrc','binarised_stack.mrc', self%smpd, .true.) !true means I am gonna set to 0 external frame to manage border effects
      write(logfhandle,*) '>>>>>>>>>>>>>STACK BINARISED SUCCESSFULLY>>>>>>>>>>>>>'
      call build_resolutions_vector(self, box, res_vec, nn_shells)
      write(logfhandle,*) "Power spectra divided into ", nn_shells, ' shells'
      allocate(counter(nn_shells+1), mask(nn_shells+1,2))
      open(unit = 17, access = 'sequential', action = 'readwrite',file = "PowerSpectraAnalysis.txt", form = 'formatted', iomsg = iom, iostat = status, position = 'append', status = 'replace')
      write(unit = 17, fmt = '(a)') '>>>>>>>>>>>>>>>>>>>>POWER SPECTRA STATISTICS>>>>>>>>>>>>>>>>>>'
      write(unit = 17, fmt = '(a)') ''
      write(unit = 17, fmt = "(a,a)")  'Input stack  ', self%powerspec_name
      write(unit = 17, fmt = "(a,a)")  'Output stack ', fname
      write(unit = 17, fmt = "(a,i0,tr1,i0,tr1,i0)") 'Dimensions ', self%img(1)%get_ldim()
      write(unit = 17, fmt = "(a,f0.2)")  'Smpd ', self%smpd
      write(unit = 17, fmt = "(a,i0)")  'N images  ', self%nstack
      write(unit = 17, fmt = '(a)') ''
      write(unit = 17, fmt = "(a)")  '-----------SELECTED PARAMETERS------------- '
      write(unit = 17, fmt = '(a)') ''
      write(unit = 17, fmt = "(a, f0.0)")  'Low pass filter ', lp
      write(unit = 17, fmt = "(a,tr1,i0.0)") 'Number of shells: ', nn_shells
      write(unit = 17, fmt = "(a,i0.0,a,i0.0)") 'Pixel split up in shells in the interval  0 -  ', &
                              & calc_fourier_index(res_vec(2,2), box, self%smpd),  &
                              & ' with step ', calc_fourier_index(res_vec(size(res_vec, dim = 1),2), box, self%smpd)
      write(unit = 17, fmt = '(a)') ''
      write(unit = 17, fmt = "(a)")  '-----------IMAGE ANALYSIS------------- '
      mask(:,1) = .false.
      mask(:,2) = .true.
      allocate(counts(self%nstack), source = 0) !nb of visible contiguous
      do n_image = 1, self%nstack
          call self%img_bin(n_image)%read('binarised_stack.mrc', n_image)
          discard = discard_ps(self, n_image)
          if(discard) write(unit = 17, fmt = '(a)') 'Empty micrograph ', n_image ; continue
          if(.not. discard) then
              rmat = self%img_bin(n_image)%get_rmat()
              counter = 0
              do i = 1, self%ldim(1)
                  do j = 1, self%ldim(2)
                      h   = -int(self%ldim(1)/2) + i - 1
                      k   = -int(self%ldim(2)/2) + j - 1
                      sh  = nint(hyp(real(h),real(k)))        !shell to which px (i,j) belongs
                      ind = minloc(abs(res_vec-sh),mask)      !corresponding shell in res_vec
                      counter(ind(1)) = counter(ind(1)) + 1   !Number of pixels per shell, it is fixed, I could calculate aside
                      !to fixxxxxxxxxxxxxxxxxxxxxxxxxxx
                      if (rmat(i,j,1) > 0.5 .and. ind(1) <= nn_shells) then !binary image, discard edges
                          res_vec(ind(1),1) = res_vec(ind(1),1) + 1. !update # of white pixel in the shell
                      endif
                  enddo
              enddo
              where(counter > 0) res_vec(:,1) = res_vec(:,1) / counter(:)  ! normalise
              write(unit = 17, fmt = "(a,tr1,i0)") 'Image', n_image-1
              res = find_res(self,n_image,box)
          endif
          if(.not. discard) then
               ax = calc_fourier_index(res, box, self%smpd)
               rmat = self%img(n_image)%get_rmat() !restore
               rmat(self%ldim(1)/2+nint(ax):self%ldim(1)/2+nint(ax)+1,self%ldim(2)/2-4:self%ldim(2)/2+4,1) = maxval(rmat)               !
               call self%img(n_image)%set_rmat(rmat) !mark visible rings
               if(.not. discard) then
                   write(unit = 17, fmt = "(a,i0.0,a)") 'Visible rings until res ', int(res), 'A'
                   limit = box/2 + nint(ax)
                   !if(res_vec(limit,2) == res_vec(size(res_vec, dim =1),2)) print *, 'high defocus image'
                   call count_visible_continuous_ring(self, n_image,limit,counts(n_image))
                   write(unit = 17, fmt = "(a,i0.0,a)") 'Nb of contiguous visible rings ', counts(n_image)
                   print *, 'N_IMAGE = ', n_image, 'visib cont rings: ', counts(n_image)
               endif
           endif
          call self%img(n_image)%write(fname,n_image)
      enddo
      close(17, status = "keep")
      deallocate(mask, counter)
      if(allocated(rmat)) deallocate(rmat)
  end subroutine process_ps_stack

   subroutine run_powerspectrum_job(self)
    class(powerspectrum), intent(inout) :: self
    ! real, allocatable :: rmat(:,:,:)
    ! real, allocatable :: rmat_aux(:,:,:)
    ! type(image) :: img_win, img_templ, img_try, img_cc
    ! logical :: yes_no, discard
    ! integer :: i, box
    ! integer :: points(3,2), px(2), j
    ! integer :: center(2)
    ! integer :: h, k, m(3), sh, cnt, sz
    ! real    :: label_mirror
    ! real    :: radius
    ! logical :: outside, bad
    ! real :: r
    ! box = 512
    ! cnt = 0
    ! call img_try%new([512,512,1],1.)
    ! call self%img_bin(1)%read('SAGAIceBinary.mrc')
    ! call delete_shells(self, 1, 512)
    ! call self%img_bin(1)%find_connected_comps(img_cc)
    ! call img_cc%write('CC.mrc')
    ! call img_win%new([box/8, box/8, 1], self%smpd)
    ! rmat = img_cc%get_rmat()
    ! allocate(rmat_aux(self%ldim(1),self%ldim(2),1), source = 0.)
    ! do i = 1, int(maxval(rmat))
    !     sz = count(abs(rmat - real(i)) < TINY)
    !     if(sz < 7 ) cycle !discard too small cc
    !     yes_no = is_symmetric(self, 1, img_cc, i, discard, label_mirror)
    !     if(.not. discard) then !discard the cc that are not in the I or II quadrant
    !         if(yes_no) then    !if it is symmetric
    !             write(logfhandle,*)'CC number ', i, 'is symm with ', label_mirror
    !             call extract_3points(self,1,img_cc,i,points)
    !                 call identify_center_radius(self,points,center,radius, bad)
    !                 if(.not. bad) then
    !                     rmat_aux = 0.
    !                     where(abs(rmat(:,:,:)-real(i)) < TINY) rmat_aux = 1.
    !                     call img_try%set_rmat(rmat_aux)
    !                     call img_try%window_slim(center-box/16, box/8, img_win, outside)
    !                     if(.not. outside) then
    !                         cnt = cnt + 1
    !                         call img_win%write('img_win.mrc', cnt)
    !                         img_templ = build_ice_template(self, box/8, radius)
    !                         call img_templ%write('img_templ.mrc', cnt)
    !                         r = img_win%real_corr(img_templ)
    !                         write(logfhandle,*)'correlation = ', r
    !                         rmat = img_cc%get_rmat()
    !                         m(:) =  minloc(abs(rmat- real(i)))
    !                         if(r > 0.4) write(logfhandle,*)'DETECTED ICE! At px ', m
    !                         h   = -box/2 + m(1) - 1
    !                         k   = -box/2 + m(2) - 1
    !                         sh  = nint(hyp(real(h),real(k)))
    !                         !write(logfhandle,*)'shell = ', sh, ' Resolution = ', calc_fourier_index(real(sh),box,self%smpd)
    !                     endif
    !                 endif
    !            endif
    !     endif
    ! enddo
    call process_ps_stack(self, 'analised_stack.mrc', 35., 10)
  end subroutine run_powerspectrum_job

  subroutine run(self)
      class(powerspectrum), intent(inout) :: self
      call self%cline%parse_oldschool()
      call self%cline%checkvar('smpd', 1)
      call self%cline%check()
      call self%p%new(self%cline)
      call self%new()
      call self%run_powerspectrum_job()
      call self%kill()
  end subroutine run

  subroutine kill_powerspectrum(self)
      class(powerspectrum), intent(inout) :: self
      integer :: i
      do i = 1, self%nstack
          call self%img(i)%kill()
          call self%img_bin(i)%kill()
      enddo
      deallocate(self%img, self%img_bin)
  end subroutine kill_powerspectrum
end module simple_test_chiara_thonrings_mod

program simple_test_chiara_thonrings
    use simple_test_chiara_thonrings_mod
    implicit none
    type(powerspectrum) :: ps
    call ps%run()
end program simple_test_chiara_thonrings
