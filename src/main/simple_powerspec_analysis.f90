!USAGE : simple_exec prg=pspec_stats smpd=1.41 stk='pspecs_saga_polii.mrc' (optional lp=20)
module simple_powerspec_analysis
include 'simple_lib.f08'
use simple_image,      only : image
use simple_stackops,   only : binarize_stack
use simple_cmdline,    only : cmdline
use simple_parameters, only : parameters

implicit none

public :: powerspectrum
private
#include "simple_local_flags.inc"

type :: powerspectrum
    private

    character(len=STDLEN)    :: powerspec_name  = ''
    character(len=STDLEN)    :: fbody           = ''
    type(image), allocatable :: img(:), img_bin(:)
    real, allocatable :: res_vec(:,:)
    type(parameters)  :: p
    type(cmdline)     :: cline
    integer     :: ldim(3)    = 0
    real        :: smpd       = 0.
    integer     :: nstack     = 0
    integer     :: nn_shells  = 0
    real        :: avg_wb_ratio_per_shell = 0.

  contains
    procedure          :: new => new_powerspectrum
    procedure, private :: build_resolutions_vector
    procedure, private :: find_res
    procedure, private :: empty
    procedure, private :: is_close_to_focus
    procedure, private :: process_ps_stack
    procedure, private :: print_resolutions_vector
    procedure, private :: run_powerspectrum_job
    procedure          :: run
    procedure, private :: kill => kill_powerspectrum
end type powerspectrum

contains

    !constructor
    subroutine new_powerspectrum(self, name, n)
        class(powerspectrum), intent(inout) :: self
        character(len=*),     intent(in)    :: name
        integer,              intent(in)    :: n
        integer :: i
        !To set the number of shells in which Fourier space
        !has to be segmented in the analysis.
        self%nn_shells = n
        self%powerspec_name = name
        call find_ldim_nptcls(self%powerspec_name, self%ldim, self%nstack, self%smpd)
        self%fbody = get_fbody(trim(self%powerspec_name), trim(fname2ext(self%powerspec_name)))
        self%ldim(3) = 1 !it's the dim of every img of the stack
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
    !in the variable count.
    subroutine count_visible_continuous_ring(self,n,lim,count)
        class(powerspectrum), intent(inout) :: self
        integer,              intent(in)    :: n
        integer,              intent(in)    :: lim  !pixel until the one the variations have to be checked
        integer,              intent(out)   :: count
        real, allocatable :: rmat_bin(:,:,:)
        rmat_bin = self%img_bin(n)%get_rmat()
        !sanity check
        if(any(rmat_bin <-0.0001) .or. any(rmat_bin>1.0001)) THROW_HARD('Non implemented non binary images; count_visible_continuous_ring')
        call count_grad_variations(rmat_bin,lim,count)
        deallocate(rmat_bin)
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
            s = shape(rmat_bin)
            allocate(imat(s(1),s(2),s(3)), source = nint(rmat_bin))
            if(s(3) .ne. 1) THROW_HARD('Non implemented for volumes; count_grad_variations')
            count   = 0               !counter for gradient variations
            fixed   = s(1)/2
            from    = s(2)/2
            to      = lim
            tmp_val = imat(fixed,from-1,1) !initialize
            do j = from, to
                if(imat(fixed,j,1) .ne. tmp_val) then
                    count = count + 1
                    tmp_val = imat(fixed,j,1)
                endif
            enddo
            count = count/2 !not to count twice (black->white->black is just one ring)
            deallocate(imat)
        end subroutine count_grad_variations
    end subroutine count_visible_continuous_ring

    !This function builds a vector that splits up the images in nn_shells shells
    !and stores the result in res_vec(:,2). In res_vec(:,1) will be saved the
    !number of white pixels in the shell corresponding to res_vec(:,2).
    subroutine build_resolutions_vector(self, box)
        class(powerspectrum), intent(inout) :: self
        integer,              intent(in)    :: box
        real    :: step
        integer :: i
        if(self%nn_shells < 2) THROW_HARD('Too low number of shells; build_resolutions_res_vec')
        if(allocated(self%res_vec)) deallocate(self%res_vec)
        allocate(self%res_vec(self%nn_shells,2), source = 0.)
        step = real(box)/real(2*self%nn_shells)
        do i = 1, self%nn_shells
            self%res_vec(i,2) = real(i)*step
            !first coumn is left empty, it's going to contain the ratio between white and tot pixel in the shell
        end do
    end subroutine build_resolutions_vector

    subroutine print_resolutions_vector(self)
        class(powerspectrum), intent(inout) :: self
        if(.not. allocated(self%res_vec)) THROW_HARD('Resolution vector hasn t been allocated yet; print_resolutions_vector')
        write(logfhandle,*) 'RESOLUTION VECTOR: '
        call vis_mat(self%res_vec)
    end subroutine print_resolutions_vector

    !This function uses the ratio (white pxls)/(black pxls) per shell
    !to estimate how far do the detected rings go.
    !The resolution is set to the shell for which the nb of white
    !pxls is less than half of the avg nb of white pxls per shell.
    function find_res(self,n_image, box) result(res)
      use gnufor2
      class(powerspectrum), intent(in) :: self
      integer,              intent(in) :: n_image
      integer,              intent(in) :: box
      real, allocatable :: rmat_bin(:,:,:)
      real    :: res
      integer :: i, cnt
      if(.not. allocated(self%res_vec)) THROW_HARD ('You have to build the resolution vector first; find_res')
      rmat_bin = self%img_bin(n_image)%get_rmat()
      if(any(rmat_bin>1.001) .or. any(rmat_bin<-0.001)) THROW_HARD('Img has to be binarised first; find_res')!sanity check
      do i = 1, size(self%res_vec,dim=1)
          if(self%res_vec(i,1) >  .5*self%avg_wb_ratio_per_shell) then
              res = calc_fourier_index(real(self%res_vec(i,2)), box, self%smpd)
              return
          endif
      enddo
      deallocate(rmat_bin)
    end function find_res

    !This subroutine is meant to discard empty power spectra images.
    !If after binarisation the # of white pixels detected in the central
    !zone of the image is less than 2% of the tot
    !# of central pixels, than it returns yes, otherwhise no.
    function empty(self, n) result(yes_no)
        class(powerspectrum), intent(inout) :: self
        integer, intent(in) :: n
        logical             :: yes_no
        real, allocatable   :: rmat(:,:,:),rmat_central(:,:,:)
        yes_no = .false.
        rmat = self%img_bin(n)%get_rmat()
        if(any(rmat > 1.001) .or. any(rmat < 0.)) THROW_HARD('Expected binary image in input; discard_ps')
        allocate(rmat_central(self%ldim(1)/2+1,self%ldim(2)/2+1,1), source = 0.)
        rmat_central(1:self%ldim(1)/2+1,1:self%ldim(2)/2+1,1) = &
        &        rmat( self%ldim(1)/2-self%ldim(1)/4 : self%ldim(1)/2+self%ldim(1)/4 , &
                     & self%ldim(2)/2-self%ldim(2)/4 : self%ldim(2)/2+self%ldim(2)/4 , 1)
        if(count(rmat_central(:,:,:) > 0.5)< 2*self%ldim(1)*self%ldim(2)*self%ldim(3)/(2*2*100)) yes_no = .true.; return
        deallocate(rmat, rmat_central)
    end function empty

   !This subroutine analyses one of the images in the ps stack
   !and determines wether it is close to focus or not. If so,
   !this image needs to be re-processed because edge detection
   !in this conditions needs specific settings.
    subroutine is_close_to_focus(self, n, box, yes_no)
        class(powerspectrum), intent(inout) :: self
        integer,              intent(in)    :: n
        integer,              intent(in)    :: box
        logical,              intent(out)   :: yes_no
        integer, allocatable :: sz(:)
        real        :: ave                      !mean sz of the connected components
        type(image) :: img, img_cc
        call img%new([self%ldim(1),self%ldim(2),1],self%smpd)
        ! call img%read('binarised_stack.mrc',n) !fetch image
        ! call img%find_connected_comps(img_cc)
        ! sz   = img_cc%size_connected_comps()
        ! ave  = real(sum(sz(:)))/real(size(sz))
        ! yes_no=.false.
        ! if(ave < 20.) yes_no=.true.
        ! write(logfhandle, *)  'N_IMAGE = ', n-1, 'avg sz of its ccs = ', ave,'close to focus =', yes_no
        ! deallocate(sz)
        ! call img%kill
        ! call img_cc%kill
    end subroutine is_close_to_focus



    !This function takes in input the name of a stack of power spectra images (fname2process),
    !and a low-pass parameter.
    !It binarises all the images in the stack and estimates how far do the detected rings go.
    !It also gives an estimation of the number of countiguous and visible rings.
    subroutine process_ps_stack(self, lp)
      use gnufor2
      use simple_stackops, only : prepare_stack
      class(powerspectrum), intent(inout) :: self
      real,                 intent(in)    :: lp
      integer, allocatable :: counter(:) !total number of pixels per shell
      integer, allocatable :: counts(:)  !nb of visible and contiguous rings
      real, allocatable :: rmat(:,:,:)
      integer           :: box
      integer           :: sh, ind(2)
      integer           :: h, k, i, j, n_image
      logical, allocatable :: mask(:,:)     !for the minloc calculation
      real                 :: res, ax
      type(image)          :: img
      character(len = 100) :: iom
      integer              :: status
      integer              :: limit    !to count nb of gradient variations until there
      logical, allocatable :: close_to_focus(:)
      logical, allocatable :: empty(:) !to keep track of empty ps
      integer, allocatable :: sz(:)
      real,    allocatable :: avg_sz_ccs(:)
      type(image) :: img_cc
      real :: avg_sz_cc_stack
      box = self%ldim(1)
      ! allocations
      allocate(close_to_focus(self%nstack), source = .false.)
      allocate(empty         (self%nstack), source = .false.)
      allocate(avg_sz_ccs    (self%nstack), source = 0.)
      call prepare_stack(trim(self%powerspec_name), trim(self%fbody)//'_prepared_stack.mrc', self%smpd, lp)
      write(logfhandle,*) '>>>>>>>>>>>>>STACK PREPARED SUCCESSFULLY>>>>>>>>>>>>>'
      call binarize_stack(trim(self%fbody)//'_prepared_stack.mrc',trim(self%fbody)//'_binarised_stack.mrc', self%smpd, .true.)!true means set to 0 external frame to manage border effects
      write(logfhandle,*) '>>>>>>>>>>>>>STACK BINARISED SUCCESSFULLY>>>>>>>>>>>>>'
      write(logfhandle,*) '>>>>>>>>>>>>>CLOSE TO FOCUS REFINEMENT>>>>>>>>>>>>>>>>'
      ! calculation of avg sz of ccs in the stack
      do n_image = 1, self%nstack
          call self%img_bin(n_image)%read(trim(self%fbody)//'_binarised_stack.mrc', n_image)
          call self%img_bin(n_image)%find_connected_comps(img_cc)
          sz = img_cc%size_connected_comps()
          avg_sz_ccs(n_image)  = real(sum(sz(:)))/real(size(sz))
          deallocate(sz)
          empty(n_image) = self%empty(n_image)
      enddo
      avg_sz_cc_stack = sum(avg_sz_ccs)/real(self%nstack)
      print *, '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      print *, 'avg_sz_cc_stack = ', avg_sz_cc_stack
      do n_image = 1, self%nstack
          if(avg_sz_ccs(n_image) < 0.8*avg_sz_cc_stack) close_to_focus(n_image)=.true.
          write(logfhandle, *)  'N_IMAGE = ', n_image-1, 'avg sz of its ccs = ', avg_sz_ccs(n_image),'close to focus =', close_to_focus(n_image)
      enddo
      call prepare_stack (trim(self%powerspec_name), trim(self%fbody)//'_prepared_stack.mrc', self%smpd, lp,close_to_focus)
      call binarize_stack(trim(self%fbody)//'_prepared_stack.mrc',trim(self%fbody)//'_binarised_stack.mrc', self%smpd,close_to_focus,.true.) !true means I am gonna set to 0 external frame to manage border effects
      call build_resolutions_vector(self, box)
      allocate(counter(self%nn_shells), mask(self%nn_shells,2)) !nb of pixels per shell
      allocate(counts (self%nstack), source = 0)                !nb of visible contiguous
      open(unit = 17, access = 'sequential', action = 'readwrite',file = "PowerSpectraAnalysis.txt", form = 'formatted', iomsg = iom, iostat = status, position = 'append', status = 'replace')
      write(unit = 17, fmt = '(a)') '>>>>>>>>>>>>>>>>>>>>POWER SPECTRA STATISTICS>>>>>>>>>>>>>>>>>>'
      write(unit = 17, fmt = '(a)') ''
      write(unit = 17, fmt = "(a,a)")  'Input stack  ', trim(self%fbody)
      write(unit = 17, fmt = "(a,a)")  'Output stack ', trim(self%fbody)//'_analysed_stack.mrc'
      write(unit = 17, fmt = "(a,i0,tr1,i0,tr1,i0)") 'Logical dim ', self%img(1)%get_ldim()
      write(unit = 17, fmt = "(a,f0.2)")             'Smpd        ', self%smpd
      write(unit = 17, fmt = "(a,i0)")               'N images    ', self%nstack
      write(unit = 17, fmt = '(a)') ''
      write(unit = 17, fmt = "(a)")  '-----------SELECTED PARAMETERS------------- '
      write(unit = 17, fmt = '(a)') ''
      write(unit = 17, fmt = "(a, f0.0)")  'Low pass filter ', lp
      write(unit = 17, fmt = "(a,tr1,i0.0)") 'Number of shells: ', self%nn_shells
      write(unit = 17, fmt = "(a,i0.0,a,i0.0)") 'Pixel split up in shells in the interval  0 -  ', &
                              & calc_fourier_index(self%res_vec(2,2), box, self%smpd),  &
                              & ' with step ', calc_fourier_index(self%res_vec(size(self%res_vec, dim = 1),2), box, self%smpd)
      write(unit = 17, fmt = '(a)') ''
      write(unit = 17, fmt = "(a)")  '-----------IMAGE ANALYSIS------------- '
      mask(:,1) = .false.
      mask(:,2) = .true.
      do n_image = 1, self%nstack
          if(empty(n_image)) write(unit = 17, fmt = '(a)')  'Empty micrograph ', n_image-1
          call self%img_bin(n_image)%read(trim(self%fbody)//'_binarised_stack.mrc', n_image)
          if(.not. empty(n_image)) then
              rmat    = self%img_bin(n_image)%get_rmat()
              counter = 0
              do i = 1, self%ldim(1)
                  do j = 1, self%ldim(2)
                      h   = -int(self%ldim(1)/2) + i - 1
                      k   = -int(self%ldim(2)/2) + j - 1
                      sh  = nint(hyp(real(h),real(k)))        !shell to which px (i,j) belongs
                      ind = minloc(abs(self%res_vec-sh),mask) !corresponding shell in res_vec
                      counter(ind(1)) = counter(ind(1)) + 1
                      if (rmat(i,j,1) > 0.5) then ! .and. ind(1) <= self%nn_shells) then !binary image
                          self%res_vec(ind(1),1) = self%res_vec(ind(1),1) + 1.   !update # of white pixel in the shell
                      endif
                  enddo
              enddo
              where(counter > 0) self%res_vec(:,1) = counter(:)/self%res_vec(:,1) ! normalise
              self%avg_wb_ratio_per_shell = sum(self%res_vec(:,1))/real(size(self%res_vec, dim = 1))
              print *, 'AVG RATIO >>>>>>>>>>>>>>>>>>>>>>>>>>>> for',n_image-1, 'is', self%avg_wb_ratio_per_shell !overwrite for every img in the stack
              write(unit = 17, fmt = "(a,tr1,i0)") 'Image', n_image-1
              call self%img_bin(n_image)%write('BinForResCalculation.mrc',n_image)
              res  = find_res(self,n_image,box)
              ax   = calc_fourier_index(res, box, self%smpd)
              rmat = self%img(n_image)%get_rmat()   !restore
              rmat(self%ldim(1)/2+nint(ax):self%ldim(1)/2+nint(ax)+1,self%ldim(2)/2-4:self%ldim(2)/2+4,1) = maxval(rmat) !mark on the image
              call self%img(n_image)%set_rmat(rmat)
              write(unit = 17, fmt = "(a,i0.0,a)") 'Visible rings until res ', int(res), 'A'
              limit = box/2 + nint(ax)
              call count_visible_continuous_ring(self, n_image,limit,counts(n_image))
              write(unit = 17, fmt = "(a,i0.0,a)") 'Nb of contiguous visible rings ', counts(n_image)
              write(logfhandle, *) 'N_IMAGE = ', n_image-1, 'visib rings until res ', res, 'A, # visib cont rings: ', counts(n_image)
              call self%img(n_image)%write(trim(self%fbody)//'_analysed_stack.mrc',n_image)
          endif
      enddo
      call self%print_resolutions_vector()
      close(17, status = "keep")
      deallocate(mask, counter, counts, close_to_focus, empty)
      if(allocated(rmat)) deallocate(rmat)
  end subroutine process_ps_stack

   subroutine run_powerspectrum_job(self, lp)
    class(powerspectrum), intent(inout) :: self
    real,                 intent(in)    :: lp
    call process_ps_stack(self,lp)
  end subroutine run_powerspectrum_job

  subroutine run(self)
      class(powerspectrum), intent(inout) :: self
      real :: lp
      call self%p%new(self%cline)
      lp = 35. !default value
      if(self%cline%defined('lp')) lp = self%p%lp
      call self%run_powerspectrum_job(lp)
      call self%kill()
  end subroutine run

  subroutine kill_powerspectrum(self)
      class(powerspectrum), intent(inout) :: self
      integer :: i
      do i = 1, self%nstack
          call self%img(i)%kill()
          call self%img_bin(i)%kill()
      enddo
      if(allocated(self%img))     deallocate(self%img)
      if(allocated(self%img_bin)) deallocate(self%img_bin)
      if(allocated(self%res_vec)) deallocate(self%res_vec)
      self%powerspec_name = ''
      self%fbody          = ''
      self%ldim      = 0
      self%smpd      = 0.
      self%nstack    = 0
      self%nn_shells = 0
      self%avg_wb_ratio_per_shell = 0.
  end subroutine kill_powerspectrum
end module simple_powerspec_analysis
