!USAGE : simple_exec prg=pspec_stats smpd=1.41 stk='pspecs_saga_polii.mrc' (optional lp=20)
module simple_genpspec_and_statistics
include 'simple_lib.f08'
use simple_image,      only : image
use simple_cmdline,    only : cmdline
use simple_parameters, only : parameters

implicit none

public :: pspec_statistics
private
#include "simple_local_flags.inc"

! module global constants
logical, parameter :: DEBUG_HERE = .false.
integer, parameter :: BOX = 512

type :: pspec_statistics
    private
    character(len=STDLEN)    :: micname  = ''
    character(len=STDLEN)    :: fbody    = ''
    type(image)        :: mic          !micrograph
    type(image)        :: ps, ps_bin   !power spectrum and its binary version
    real, allocatable  :: res_vec(:,:)
    type(parameters)   :: p
    type(cmdline)      :: cline
    real        :: avg_wb_ratio  !avg ratio white/black pxls among shells
    integer     :: ldim_mic(3)                      = 0  !dim of the mic
    integer     :: ldim(3)                          = 0  !dim of the ps
    real        :: smpd                             = 0. !smpd
    integer     :: nn_shells                        = 0  !nb of shells in which to divide the ps
    real        :: lp                               = 0. !low pass filter
    logical     :: close_to_focus = .false.
    integer, allocatable :: npxls_per_shell(:)           !number of pixles per shell

  contains
    procedure          :: new => new_pspec_statistics
    procedure, private :: build_resolutions_vector
    procedure, private :: per_shell_entropy
    procedure, private :: find_res
    procedure, private :: empty
    procedure, private :: process_ps
    procedure, private :: print_resolutions_vector
    procedure, private :: print_info
    procedure          :: run
    procedure, private :: kill => kill_pspec_statistics
end type pspec_statistics

contains

    !constructor
    subroutine new_pspec_statistics(self, name, n)
        use simple_defs
        class(pspec_statistics), intent(inout) :: self
        character(len=*),        intent(in)    :: name
        integer,                 intent(in)    :: n
        integer :: nptcls
        !To set the number of shells in which Fourier space
        !has to be segmented in the analysis.
        self%nn_shells = n
        self%micname   = name
        call find_ldim_nptcls(self%micname, self%ldim_mic, nptcls, self%smpd)
        self%fbody   = get_fbody(trim(self%micname), trim(fname2ext(self%micname)))
        call self%mic%new(self%ldim_mic, self%smpd)
        self%ldim(1) = BOX
        self%ldim(2) = BOX
        self%ldim(3) = 1
        call self%ps_bin%new(self%ldim, self%smpd)
        call self%mic%read(self%micname)
        !power spectrum generation
        call self%mic%mic2spec(BOX, 'power', LP_PSPEC_BACKGR_SUBTR, self%ps) !LP_PSPEC_BACKGR_SUBTR is 20.
        call self%ps%scale_pixels([1.,real(BOX)*2.+1.]) ! to set n in entropy calc
        call self%ps%write(trim(self%fbody)//'_generated_ps.mrc')
    end subroutine new_pspec_statistics

    !This subroutine is meant to count the nb of visible and contiguous
    !rings are detected in the binarization of a power-spectra stack.
    !It saves the nb of vis rings in the i-th image of the stack
    !in the variable count.
    subroutine count_visible_continuous_ring(self,lim,count)
        class(pspec_statistics), intent(inout) :: self
        integer,                 intent(in)    :: lim  !pixel until the one the variations have to be checked
        integer,                 intent(out)   :: count
        real, allocatable :: rmat_bin(:,:,:)
        rmat_bin = self%ps_bin%get_rmat()
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
            integer, allocatable :: imat(:,:,:) !input is meant to be binary, so I consider the integer matrix associated
            integer :: s(3) !shape of the matrix
            integer :: j
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

    !This function builds a vector that splits up the ps in nn_shells shells
    !and stores the pixel limit in res_vec(:,2). In res_vec(:,1) will be saved the
    !number of white pixels in the shell corresponding to res_vec(:,2).
    !The interval considered is [Nyq, 35] A.
    subroutine build_resolutions_vector(self, box)
        class(pspec_statistics), intent(inout) :: self
        integer,                 intent(in)    :: box
        integer :: i
        real    :: step
        real :: min_val, max_val
        if(self%nn_shells < 2) THROW_HARD('Too low number of shells; build_resolutions_res_vec')
        if(allocated(self%res_vec)) deallocate(self%res_vec)
        allocate(self%res_vec(self%nn_shells-1,2), source = 0.)
        min_val = calc_fourier_index(35., box, self%smpd)
        max_val = calc_fourier_index(2.*self%smpd, box, self%smpd)
        step = (max_val-min_val)/real(self%nn_shells) !go from Nyq to 35 A
        do i = 1, self%nn_shells-1
            self%res_vec(i,2) = min_val + real(i)*step
            !first coumn is left empty, it's going to contain the ratio between white and tot pixel in the shell
        end do
    end subroutine build_resolutions_vector

   ! Reporting purposes
    subroutine print_resolutions_vector(self)
        class(pspec_statistics), intent(inout) :: self
        if(.not. allocated(self%res_vec)) THROW_HARD('Resolution vector hasn t been allocated yet; print_resolutions_vector')
        write(logfhandle,*) 'RESOLUTION VECTOR: '
        call vis_mat(self%res_vec)
    end subroutine print_resolutions_vector

    !This function uses the ratio (white pxls)/(black pxls) per shell
    !to estimate how far do the detected rings extend.
    !The resolution is set to the shell for which the nb of white
    !pxls is less than 0.3 of the avg nb of white pxls per shell.
    function find_res(self, box) result(res)
      use gnufor2
      class(pspec_statistics), intent(in) :: self
      integer,                 intent(in) :: box
      real, allocatable :: rmat_bin(:,:,:)
      real              :: res
      integer           :: i
      if(.not. allocated(self%res_vec)) THROW_HARD ('You have to build the resolution vector first; find_res')
      rmat_bin = self%ps_bin%get_rmat()
      if(any(rmat_bin>1.001) .or. any(rmat_bin<-0.001)) THROW_HARD('Img has to be binarised first; find_res')!sanity check
      do i = 1, size(self%res_vec,dim=1)
          if(.not. self%close_to_focus) then
              if(self%res_vec(i,1) < .3*self%avg_wb_ratio) then
                  res = calc_fourier_index(real(self%res_vec(i,2)), box, self%smpd)
                  return
              endif
          else  !high focus images have more noise in the background, need higher thresh
              if(self%res_vec(i,1) < .4*self%avg_wb_ratio) then
                  res = calc_fourier_index(real(self%res_vec(i,2)), box, self%smpd)
                  return
              endif
          endif
      enddo
      deallocate(rmat_bin)
    end function find_res

   ! This function calculates the entropy per shell.
    subroutine per_shell_entropy(self)
        use simple_math
        class(pspec_statistics), intent(inout) :: self
        real,    allocatable :: X(:)
        real,    allocatable :: rmat(:,:,:)
        logical, allocatable :: mask(:,:)
        real    :: e(self%nn_shells-1)
        real    :: avg_entropy_per_shell, stdev_entropy_per_shell
        real    :: res
        integer :: cnt(self%nn_shells-1)
        integer :: ind(2)
        integer :: i, j
        integer :: h, k, n_shell
        integer :: sh
        do i = 1, size(self%res_vec, dim = 1)
            res =  calc_fourier_index(real(self%res_vec(i,2)), BOX, self%smpd)
        enddo
        allocate(mask(self%nn_shells-1,2),  source = .false.) !nb of pixels per shell
        mask(:,2) = .true.
        rmat = self%ps%get_rmat()
        ! Calculate nb of pixels per shell
        if(allocated(self%npxls_per_shell)) deallocate(self%npxls_per_shell)
        allocate(self%npxls_per_shell(self%nn_shells-1), source  = 0)
        do i = 1, self%ldim(1)
            do j = 1, self%ldim(2)
                h   = -int(self%ldim(1)/2) + i - 1
                k   = -int(self%ldim(2)/2) + j - 1
                sh  =  nint(hyp(real(h),real(k)))       !shell to which px (i,j) belongs
                ind = minloc(abs(self%res_vec-sh),mask) !corresponding shell in res_vec
                self%npxls_per_shell(ind(1)) = self%npxls_per_shell(ind(1)) + 1
            enddo
        enddo
        cnt = 0 ! initialisation
        ! Calculate entropy per shell
        do n_shell = 1, self%nn_shells-1
            if(allocated(X)) deallocate(X)
            allocate(X(self%npxls_per_shell(n_shell)), source = 0.)
            do i = 1, self%ldim(1)
                do j = 1, self%ldim(2)
                    h   = -int(self%ldim(1)/2) + i - 1
                    k   = -int(self%ldim(2)/2) + j - 1
                    sh  =  nint(hyp(real(h),real(k)))       !shell to which px (i,j) belongs
                    ind = minloc(abs(self%res_vec-sh),mask) !corresponding shell in res_vec
                    if(ind(1) == n_shell) then
                        cnt(ind(1)) = cnt(ind(1)) + 1
                        X(cnt(ind(1))) = rmat(i,j,1)
                    endif
                 enddo
            enddo
            e(n_shell) = entropy(X,BOX/4)
            print *, 'shell ', n_shell, 'entropy ', e(n_shell)
        enddo
        ! Calculation of entropy statistics
        avg_entropy_per_shell   = 0.
        stdev_entropy_per_shell = 0.
        avg_entropy_per_shell = sum(e)/real(self%nn_shells-1)
        do n_shell = 1, self%nn_shells-1
            stdev_entropy_per_shell = stdev_entropy_per_shell+ (e(n_shell)-avg_entropy_per_shell)**2
        enddo
        stdev_entropy_per_shell = sqrt(stdev_entropy_per_shell/real(self%nn_shells-2))
        write(logfhandle, *) 'avg_entropy_per_shell   = ', avg_entropy_per_shell
        write(logfhandle, *) 'stdev_entropy_per_shell = ', stdev_entropy_per_shell
    end subroutine per_shell_entropy

    !This subroutine is meant to discard empty power spectra images.
    !If after binarisation the # of white pixels detected in the central
    !zone of the image is less than 2% of the tot
    !# of central pixels, than it returns yes, otherwhise no.
    !The central zone is a rectangle with dim self%ldim/2
    function empty(self) result(yes_no)
        class(pspec_statistics), intent(inout) :: self
        logical             :: yes_no
        real, allocatable   :: rmat(:,:,:),rmat_central(:,:,:)
        yes_no = .false.
        rmat = self%ps_bin%get_rmat()
        if(any(rmat > 1.001) .or. any(rmat < -0.001)) THROW_HARD('Expected binary image in input; discard_ps')
        allocate(rmat_central(self%ldim(1)/2+1,self%ldim(2)/2+1,1), source = 0.)
        rmat_central(1:self%ldim(1)/2+1,1:self%ldim(2)/2+1,1) = &
               & rmat( self%ldim(1)/2-self%ldim(1)/4 : self%ldim(1)/2+self%ldim(1)/4 , &
               &       self%ldim(2)/2-self%ldim(2)/4 : self%ldim(2)/2+self%ldim(2)/4 , 1)
        if(count(rmat_central(:,:,:) > 0.5)< 2*self%ldim(1)*self%ldim(2)*self%ldim(3)/(2*2*100)) yes_no = .true.; return
        deallocate(rmat, rmat_central)
    end function empty

    subroutine print_info(self)
        class(pspec_statistics), intent(inout) :: self
        character(len = 100) :: iom
        integer              :: status
        open(unit = 17, access = 'sequential', action = 'readwrite',file = trim(self%fbody)//"PowerSpectraAnalysis.txt", form = 'formatted', iomsg = iom, iostat = status, position = 'append', status = 'replace')
        write(unit = 17, fmt = '(a)') '>>>>>>>>>>>>>>>>>>>>POWER SPECTRA STATISTICS>>>>>>>>>>>>>>>>>>'
        write(unit = 17, fmt = '(a)') ''
        write(unit = 17, fmt = "(a,a)")  'Input stack  ', trim(self%fbody)
        write(unit = 17, fmt = "(a,a)")  'Output stack ', trim(self%fbody)//'_analysed_stack.mrc'
        write(unit = 17, fmt = "(a,i0,tr1,i0,tr1,i0)") 'Logical dim mic ', self%mic%get_ldim()
        write(unit = 17, fmt = "(a,i0,tr1,i0,tr1,i0)") 'Logical dim ps  ', self%ps%get_ldim()
        write(unit = 17, fmt = "(a,f0.2)")             'Smpd            ', self%smpd
        write(unit = 17, fmt = '(a)') ''
        write(unit = 17, fmt = "(a)")  '-----------SELECTED PARAMETERS------------- '
        write(unit = 17, fmt = '(a)') ''
        write(unit = 17, fmt = "(a, f0.0)")    'Low pass filter ', self%lp
        write(unit = 17, fmt = "(a,tr1,i0.0)") 'Nb of shells    ', self%nn_shells
        write(unit = 17, fmt = "(a,f5.2,a,i2)") 'Pixel split up in shells in the interval', 2.*self%smpd, ' - ', &
                                & calc_fourier_index(self%res_vec(1,2), BOX, self%smpd)
        write(unit = 17, fmt = '(a)') ''
        write(unit = 17, fmt = "(a)")  '-----------IMAGE ANALYSIS------------- '
    end subroutine print_info

    !This is the core subroutine of this module. Preprocessing and binarization of the
    !ps is performed here. Then estimation of how far do the detected rings extend.
    ! Finally estimation of the number of countiguous rings.
    subroutine process_ps(self, lp)
      class(pspec_statistics), intent(inout) :: self
      real,                    intent(in)    :: lp
      type(image)          :: ps_ccs
      real,    allocatable :: rmat(:,:,:)
      logical, allocatable :: mask(:,:)     !for the minloc calculation
      integer, allocatable :: sz(:)
      real                 :: res
      real                 :: avg_sz_ccs, stdev_sz_ccs
      real                 :: ax
      logical              :: empty         !to keep track of empty ps
      integer              :: ind(2)
      integer              :: h, k, i, j
      integer              :: sh
      integer              :: counts        !nb of visible and contiguous rings
      integer              :: limit         !to count nb of gradient variations until there
      self%lp = lp
      call build_resolutions_vector(self, BOX)
      call self%per_shell_entropy() !to discard at this stage fallacious data, TO WORK ON THIS,
      call prepare_ps(self)
      write(logfhandle,*) '>>>>>>>>>>>>>PREPARED SUCCESSFULLY>>>>>>>>>>>>>'
      call binarize_ps(self,.true.) !true means set to 0 external frame to manage border effects
      write(logfhandle,*) '>>>>>>>>>>>>>BINARISED SUCCESSFULLY>>>>>>>>>>>>>'
      ! calculation of avg sz of ccs in the stack
      call self%ps_bin%find_connected_comps(ps_ccs)
      sz = ps_ccs%size_connected_comps()
      call ps_ccs%kill
      avg_sz_ccs = 0. !initialise
      avg_sz_ccs  = real(sum(sz(:)))/real(size(sz))
      empty = self%empty()
      stdev_sz_ccs = 0.
      do i = 1, size(sz)
          stdev_sz_ccs = stdev_sz_ccs + (sz(i)-avg_sz_ccs)**2
      enddo
      stdev_sz_ccs = sqrt(stdev_sz_ccs/real(size(sz)-1))
      deallocate(sz)
      !If it's close to focus then the size of the ccs
      !is small everywhere (broken rings), so the deviation should be small. Viceversa,
      !if it's not close to focus, then there are some big ccs (the rings) and small ones
      !(noise); so the stdev should be bigger.
      if(stdev_sz_ccs< avg_sz_ccs) self%close_to_focus=.true.
      write(logfhandle, *)' low defocus =', self%close_to_focus
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! INSERT HERE DIFFERENT APPROACH FOR BINARIZATION AND
      ! PREPROCESSING OF CLOSE TO FOCUS PS
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(mask(self%nn_shells-1,2), source = .false.)
      counts = 0            !nb of visible contiguous
      mask(:,2) = .true.
      call self%print_info()
      if(empty) write(unit = 17, fmt = '(a)')  'Empty micrograph   '
      call self%ps_bin%read(trim(self%fbody)//'_binarized.mrc')
      if(.not. empty) then
          rmat    = self%ps_bin%get_rmat()
          do i = 1, self%ldim(1)
              do j = 1, self%ldim(2)
                  h   = -int(self%ldim(1)/2) + i - 1
                  k   = -int(self%ldim(2)/2) + j - 1
                  sh  =  nint(hyp(real(h),real(k)))       !shell to which px (i,j) belongs
                  ind = minloc(abs(self%res_vec-sh),mask) !corresponding shell in res_vec
                  if (rmat(i,j,1) > 0.5) self%res_vec(ind(1),1) = self%res_vec(ind(1),1) + 1.   !update # of white pixel in the shell
              enddo
          enddo
          self%res_vec(:,1) = self%res_vec(:,1)/self%npxls_per_shell(:)
          self%avg_wb_ratio = sum(self%res_vec(:,1))/real(size(self%res_vec, dim = 1))
          if(DEBUG_HERE) write(logfhandle,*) 'AVG RATIO', self%avg_wb_ratio
          res  = find_res(self,BOX)
          ax   = calc_fourier_index(res, BOX, self%smpd)
          rmat = self%ps%get_rmat() !overwrite
          ! Draw on the power spectrum the estimated resolution interval
          rmat(self%ldim(1)/2+nint(ax):self%ldim(1)/2+nint(ax)+1,self%ldim(2)/2-4:self%ldim(2)/2+4,1) = maxval(rmat)
          call self%ps%set_rmat(rmat)
          deallocate(rmat)
          write(unit = 17, fmt = "(a,f4.1,a)") 'Visible rings res   ', res,' A'
          limit = BOX/2 + nint(ax)
          call count_visible_continuous_ring(self,limit,counts)
          write(unit = 17, fmt = "(a,i2)")    'Nb contiguous rings  ', counts
          write(unit = 17, fmt = "(a)") ' '
          write(logfhandle, *) 'Visib rings res ', trim(real2str(res)), ' A, # visib cont rings: ', counts
          call self%ps%write(trim(self%fbody)//'_analysed.mrc')
      endif
      call self%print_resolutions_vector()
      close(17, status = "keep")
      deallocate(mask)
  contains

      !Preprocessing steps: -) elimination of central spot
      !                     -) tvfiltering
      !                     -) background subtraction
      subroutine prepare_ps(self)
          use simple_tvfilter, only: tvfilter
          class(pspec_statistics), intent(inout) :: self
          type(tvfilter) :: tvf
          real           :: lambda ! tvf parameter
          call manage_central_spot(self)
          lambda = 1.
          call tvf%new()
          call tvf%apply_filter(self%ps,lambda)
          call self%ps%write(trim(self%fbody)//'_tvfiltered_ps.mrcs')
          call tvf%kill
          call self%ps%subtr_backgr(self%lp)
          call self%ps%write(trim(self%fbody)//'_backsubtr_ps.mrcs')
      end subroutine prepare_ps

      ! This subroutine extract a circle from the central part of the
      ! power spectrum, until 30 A and replace these gray values with
      ! the average of the gray values in the ps outside that circle.
      ! It is meant to reduce the influence of the central pixel.
      subroutine manage_central_spot(self)
          class(pspec_statistics), intent(inout) :: self
          real, allocatable :: rmat(:,:,:)
          real, parameter   :: limit = 30.   !discard until 30 A
          real              :: avg
          integer           :: sh
          integer           :: i, j, k, h
          integer           :: cnt
          allocate(rmat(self%ldim(1),self%ldim(2),1), source = self%ps%get_rmat())
          cnt = 0
          do i = 1, self%ldim(1)
              do j = 1, self%ldim(2)
                  h   = -int(self%ldim(1)/2) + i - 1
                  k   = -int(self%ldim(2)/2) + j - 1
                  sh  = nint(hyp(real(h),real(k))) !shell to which px
                  if(sh < limit) then
                       rmat(i,j,1 ) = -10000.
                  else                       !calculate avg value pxls outside central shell
                      avg = avg + rmat(i,j,1)
                      cnt = cnt + 1
                  endif
              enddo
          enddo
          avg = avg/cnt
          where (abs(rmat+10000.) < TINY) rmat = avg
          call self%ps%set_rmat(rmat)
      end subroutine manage_central_spot

      ! This subroutine performs binarization of a ps image
      ! using Canny edge detection with automatic threshold
      ! selection. If input discard_borders is set to .true.
      ! then, after binarization, white pixels detected in the
      ! borders are set back to zero. It is done to get rid
      ! of border effects. The size of the borders is 1/10
      ! of the dimension of the image.
      subroutine binarize_ps(self, discard_borders)
          use simple_segmentation, only : canny
          class(pspec_statistics), intent(inout) :: self
          logical, optional,       intent(in) :: discard_borders
          real,    allocatable :: rmat_bin(:,:,:)
          logical     :: ddiscard_borders
          integer     :: ldim(3)
          type(image) :: tmp_img
          call tmp_img%new([self%ldim(1),self%ldim(2),1], self%smpd)
          call tmp_img%copy(self%ps)
          ddiscard_borders = .false.
          if(present(discard_borders)) ddiscard_borders = discard_borders
          call canny(self%ps)
          call self%ps_bin%copy(self%ps)
          !restore
          call self%ps%copy(tmp_img)
          call tmp_img%kill
          if(ddiscard_borders) then
              allocate(rmat_bin(self%ldim(1),self%ldim(2),1), source = self%ps_bin%get_rmat())
              rmat_bin(1:int(self%ldim(1)/10), :, 1) = 0.                        !bottom horizontal border
              rmat_bin(self%ldim(1)-int(self%ldim(1)/10):self%ldim(1),:, 1) = 0. !top horizontal border
              rmat_bin(:,1:int(self%ldim(2)/10), 1) = 0.                         !bottom vertical border
              rmat_bin(:,self%ldim(2)-int(self%ldim(2)/10):self%ldim(2), 1) = 0. !top vertical border
              call self%ps_bin%set_rmat(rmat_bin)
              deallocate(rmat_bin)
          endif
          call self%ps_bin%write(trim(self%fbody)//'_binarized.mrc')
      end subroutine binarize_ps
  end subroutine process_ps

  subroutine run(self)
      class(pspec_statistics), intent(inout) :: self
      real :: lp
      call self%p%new(self%cline)
      lp = 35. !default value
      if(self%cline%defined('lp')) lp = self%p%lp
      call process_ps(self,lp)
      call self%kill()
  end subroutine run

  subroutine kill_pspec_statistics(self)
      class(pspec_statistics), intent(inout) :: self
      call self%mic%kill()
      call self%ps%kill()
      call self%ps_bin%new (self%ldim, self%smpd)
      if(allocated(self%res_vec))         deallocate(self%res_vec)
      if(allocated(self%npxls_per_shell)) deallocate(self%npxls_per_shell)
      self%micname      = ''
      self%fbody        = ''
      self%ldim         = 0
      self%smpd         = 0.
      self%nn_shells    = 0
      self%lp           = 0.
      self%avg_wb_ratio = 0.
  end subroutine kill_pspec_statistics
end module simple_genpspec_and_statistics

! SUBROUTINES MIGHT BE USEFUL
! subroutine log_ps(self)
!     class(pspec_statistics), intent(inout) :: self
!     real, allocatable :: rmat(:,:,:)
!     call self%ps%scale_pixels([1.,real(self%ldim(1))*2.+1.]) !not to have problems in calculating the log
!     rmat = self%ps%get_rmat()
!     rmat = log(rmat)
!     call self%ps%set_rmat(rmat)
!     deallocate(rmat)
! end subroutine log_ps
