!USAGE : simple_exec prg=pspec_stats smpd=1.41 filetab='sjhskl.mrc'
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
logical, parameter :: DEBUG_HERE = .true.
integer, parameter :: BOX        = 512 !ps size
real,    parameter :: LOW_LIM    = 30. !30 A, lower limit for resolution. Before that, we discard
real,    parameter :: UP_LIM     = 4.  !upper limit for resolution in the entropy calculation.
integer, parameter :: N_BINS     = 64  !number of bins for hist
real,    parameter :: LAMBDA     = 1.  !for tv filtering

type :: pspec_statistics
    private
    character(len=STDLEN)    :: micname  = ''
    character(len=STDLEN)    :: fbody    = ''
    type(image)           :: mic                  !micrograph
    type(image)           :: ps, ps_bin, ps_ccs   !power spectrum, its binary version, its connected components
    integer, allocatable  :: res_vec(:)
    type(parameters)      :: p
    type(cmdline)         :: cline
    integer     :: ldim_mic(3) = 0       !dim of the mic
    integer     :: ldim(3)     = 0       !dim of the ps
    real        :: smpd                  !smpd
    integer     :: nn_shells   = 0       !nb of shells in which to divide the ps
    real        :: score       = 0.      !score, to keep/discard the mic
    logical     :: fallacious  = .false. !whether the mic has to be discarded or not

  contains
    procedure          :: new => new_pspec_statistics
    procedure, private :: build_resolutions_vector
    procedure, private :: per_shell_entropy
    procedure, private :: empty
    procedure, private :: calc_weighted_avg_sz_ccs
    procedure, private :: process_ps
    procedure, private :: print_resolutions_vector
    procedure, private :: print_info
    procedure          :: get_score
    procedure          :: run
    procedure          :: kill => kill_pspec_statistics
end type pspec_statistics

contains

    !constructor
    subroutine new_pspec_statistics(self, name, cline_smpd, n)
        use simple_defs
        class(pspec_statistics), intent(inout) :: self
        character(len=*),        intent(in)    :: name
        real,                    intent(in)    :: cline_smpd
        integer,                 intent(in)    :: n
        integer     :: nptcls
        real        :: smpd
        type(image) :: hist_stretch  !to perform histogram stretching (see Adiga s paper)
        !set parameters
        self%smpd      = cline_smpd
        self%nn_shells = n
        self%micname   = name
        call find_ldim_nptcls(self%micname, self%ldim_mic, nptcls, smpd)
        self%fbody   = get_fbody(trim(self%micname),trim(fname2ext(self%micname)))
        call self%mic%new(self%ldim_mic, self%smpd)
        call self%mic%read(self%micname)
        self%ldim(1) = BOX
        self%ldim(2) = BOX
        self%ldim(3) = 1
        call self%ps_bin%new(self%ldim, self%smpd)
        call self%ps_ccs%new(self%ldim, self%smpd)
        !power spectrum generation
        call self%mic%mic2spec(BOX, 'power',2.82*LOW_LIM, self%ps)
        if(DEBUG_HERE) call self%ps%write(trim(self%fbody)//'_generated_ps.mrc')
        ! call hist_stretch%new(self%ldim, self%smpd)  !histogram stretching, to try to reduce the influence of the central spot for entropy calculation
        ! call self%ps%hist_stretching(hist_stretch)
        ! if(DEBUG_HERE) call hist_stretch%write(trim(self%fbody)//'_hist_stretch_ps.mrc')
        ! call self%ps%copy(hist_stretch)
    end subroutine new_pspec_statistics

    !This function builds a vector that splits up the ps in nn_shells shells
    !and stores the pixel limit in res_vec. The interval considered is [Nyq, 30] A.
    subroutine build_resolutions_vector(self, box)
        class(pspec_statistics), intent(inout) :: self
        integer,                 intent(in)    :: box
        integer :: step, i
        integer :: min_val, max_val
        if(self%nn_shells < 2) THROW_HARD('Too low number of shells; build_resolutions_res_vec')
        if(allocated(self%res_vec)) deallocate(self%res_vec)
        allocate(self%res_vec(self%nn_shells), source = 0)
        min_val = calc_fourier_index(LOW_LIM,      BOX, self%smpd)
        max_val = calc_fourier_index(2.*self%smpd, BOX, self%smpd)
        step = (max_val-min_val)/self%nn_shells        !go from Nyq to 30 A
        do i = 1, self%nn_shells
            self%res_vec(i) = min_val-step/2 + i*step  !-step/2 to have the central spot, not the upper limit
        end do
    end subroutine build_resolutions_vector

   ! Reporting purposes
    subroutine print_resolutions_vector(self)
        class(pspec_statistics), intent(inout) :: self
        integer :: j
        if(.not. allocated(self%res_vec)) THROW_HARD('Resolution vector hasn t been allocated yet; print_resolutions_vector')
        write(logfhandle,*) 'RESOLUTION VECTOR: '
         do j = 1, size(self%res_vec)
            write(logfhandle,*) self%res_vec(j), calc_lowpass_lim(self%res_vec(j),BOX, self%smpd)
        enddo
    end subroutine print_resolutions_vector

   ! This function calculates the entropy per shell.
    subroutine per_shell_entropy(self)
        use simple_math
        class(pspec_statistics), intent(inout) :: self
        real,    allocatable :: X(:)
        real,    allocatable :: rmat(:,:,:)
        real    :: e(self%nn_shells)
        real    :: avg_entropy_per_shell, stdev_entropy_per_shell
        integer :: cnt(self%nn_shells)
        integer :: npxls_per_shell(self%nn_shells)
        integer :: ind(1)
        integer :: i, j
        integer :: h, k, sh, n_shell
        integer :: ULim_Findex, LLim_Findex
        ULim_Findex = calc_fourier_index(max(2.*self%smpd,UP_LIM), BOX,self%smpd)! Everything at a resolution higher than UP_LIM A is discarded
        LLim_Findex = calc_fourier_index(LOW_LIM,BOX,self%smpd)                  ! Everything at a resolution lower than LOW_LIM A is discarded
        rmat = self%ps%get_rmat()
        if(DEBUG_HERE) call self%ps%write(trim(self%fbody)//'_ps_for_entropy.mrc')
        ! Calculate nb of pixels per shell
        npxls_per_shell = 0
        do i = 1, BOX
            do j = 1, BOX
                h   = -int(BOX/2) + i - 1
                k   = -int(BOX/2) + j - 1
                sh  =  nint(hyp(real(h),real(k)))                !shell to which px (i,j) belongs
                if(sh < ULim_Findex .and. sh > LLim_Findex) then !discard outside LOW_LIM - UP_LIM A
                    ! It's inconsistent with the resolution vector, because the res vec goes from Nyq to LOW_LIM, not UP_LIM
                    ind = minloc(abs(self%res_vec-sh))           !corresponding shell in res_vec
                    npxls_per_shell(ind(1)) = npxls_per_shell(ind(1)) + 1
                endif
            enddo
        enddo
        cnt = 0 ! initialisation
        ! Calculate entropy per shell
        do n_shell = 1, self%nn_shells
            if(allocated(X)) deallocate(X)
            allocate(X(npxls_per_shell(n_shell)), source = 0.)
            do i = 1, BOX
                do j = 1, BOX
                    h   = -int(BOX/2) + i - 1
                    k   = -int(BOX/2) + j - 1
                    sh  =  nint(hyp(real(h),real(k)))
                    if(sh < ULim_Findex .and. sh > LLim_Findex) then
                        ind = minloc(abs(self%res_vec-sh))
                        if(ind(1) == n_shell) then  !The pixel [i,j,1] belongs to the considered shell
                            cnt(n_shell)    = cnt(n_shell) + 1
                            X(cnt(n_shell)) =  rmat(i,j,1)
                        endif
                    endif
                 enddo
            enddo
            e(n_shell) = entropy_shells(X,maxval(rmat),minval(rmat),4*N_BINS)
            write(logfhandle,*) 'shell ', n_shell, 'entropy per shells ', e(n_shell)
        enddo
        ! Calculation of entropy statistics
        avg_entropy_per_shell   = 0.
        stdev_entropy_per_shell = 0.
        avg_entropy_per_shell = sum(e)/real(self%nn_shells) !discard first shell
        do n_shell = 1, self%nn_shells
            stdev_entropy_per_shell = stdev_entropy_per_shell+ (e(self%nn_shells)-avg_entropy_per_shell)**2
        enddo
        stdev_entropy_per_shell = sqrt(stdev_entropy_per_shell/real(self%nn_shells-1))
        write(logfhandle, *) 'avg_entropy_per_shell   = ', avg_entropy_per_shell
        write(logfhandle, *) 'stdev_entropy_per_shell = ', stdev_entropy_per_shell
        write(logfhandle, *) 'avg + std deviation     = ', avg_entropy_per_shell + stdev_entropy_per_shell
    end subroutine per_shell_entropy

    !This subroutine is meant to discard fallacious power spectra images
    !characterized by producing an 'empty' binarization.
    !If after binarization the # of white pixels detected in the central
    !zone of the image is less than 2% of the tot
    !# of central pixels, than it returns yes, otherwhise no.
    !The central zone is a rectangle with dim BOX/2
    function empty(self) result(yes_no)
        class(pspec_statistics), intent(inout) :: self
        logical             :: yes_no
        real, allocatable   :: rmat(:,:,:),rmat_central(:,:,:)
        yes_no = .false.
        rmat = self%ps_bin%get_rmat()
        if(any(rmat > 1.001) .or. any(rmat < -0.001)) THROW_HARD('Expected binary image in input; discard_ps')
        allocate(rmat_central(BOX/2+1,BOX/2+1,1), source = 0.)
        rmat_central(1:BOX/2+1,1:BOX/2+1,1) = &
               & rmat( BOX/2-BOX/4 : BOX/2+BOX/4 , &
               &       BOX/2-BOX/4 : BOX/2+BOX/4 , 1)
        if(count(rmat_central(:,:,:) > 0.5)< 2*BOX*BOX/(2*2*100)) yes_no = .true.; return
        deallocate(rmat, rmat_central)
    end function empty

    subroutine print_info(self)
        class(pspec_statistics), intent(inout) :: self
        character(len = 100) :: iom
        integer              :: status
        open(unit = 17, access = 'sequential', action = 'read(write)',file = trim(self%fbody)//"PowerSpectraAnalysis.txt", form = 'formatted', iomsg = iom, iostat = status, position = 'append', status = 'replace')
        write(unit = 17, fmt = '(a)') '>>>>>>>>>>>>>>>>>>>>POWER SPECTRA STATISTICS>>>>>>>>>>>>>>>>>>'
        write(unit = 17, fmt = '(a)') ''
        write(unit = 17, fmt = "(a,a)")  'Input mic  ', trim(self%fbody)
        write(unit = 17, fmt = "(a,i0,tr1,i0,tr1,i0)") 'Logical dim mic ', self%mic%get_ldim()
        write(unit = 17, fmt = "(a,i0,tr1,i0,tr1,i0)") 'Logical dim ps  ', self%ps%get_ldim()
        write(unit = 17, fmt = "(a,f0.2)")             'Smpd            ', self%smpd
        write(unit = 17, fmt = '(a)') ''
        write(unit = 17, fmt = "(a)")  '-----------SELECTED PARAMETERS------------- '
        write(unit = 17, fmt = '(a)') ''
        write(unit = 17, fmt = "(a,tr1,i0.0)") 'Nb of shells    ', self%nn_shells
        write(unit = 17, fmt = "(a,f5.2,a,f5.2, a)") 'Pixel split up in shells in the interval', 2.*self%smpd, ' - ', LOW_LIM,' A'
    end subroutine print_info


   ! This subroutine calculates the weighted avg of the size
   ! of the connected components of the binary version of
   ! the power spectrum. The aim is to discard fallacious
   ! micrographs based on the analysis of their power
   ! spectra. The valid ones should have rings in the
   ! center, so big connected components toward the
   ! center of the image.
    subroutine calc_weighted_avg_sz_ccs(self, sz)
        class(pspec_statistics), intent(inout) :: self
        integer,                 intent(in)    :: sz(:)
        integer, allocatable :: imat(:,:,:)
        integer, allocatable :: imat_sz(:,:,:)
        integer :: ULim_Findex, LLim_Findex
        integer :: h,k,sh
        integer :: i, j
        real    :: denom  !denominator of the formula
        real    :: a
        ULim_Findex = calc_fourier_index(max(2.*self%smpd,UP_LIM), BOX,self%smpd)! Everything at a resolution higher than UP_LIM A is discarded
        LLim_Findex = calc_fourier_index(LOW_LIM,BOX,self%smpd)                  ! Everything at a resolution lower than LOW_LIM A is discarded
        denom = 0.
        call self%ps_ccs%elim_cc([1,BOX*4]) !eliminate connected components with size one
        call self%ps_ccs%write(trim(self%fbody)//'_polished_ccs.mrc')
        imat = nint(self%ps_ccs%get_rmat())
        allocate(imat_sz(BOX,BOX,1), source = 0)
        call generate_mat_sz_ccs(imat,imat_sz,sz)
        self%score = 0.  !initialise
        do i = 2, BOX-1  !discard border effects due to binarization
            do j = 2, BOX-1
                if(imat(i,j,1) > 0) then
                    h   = -int(BOX/2) + i - 1
                    k   = -int(BOX/2) + j - 1
                    sh  =  nint(hyp(real(h),real(k)))
                    if( sh > ULim_Findex .or. sh < LLim_Findex ) then
                         if(DEBUG_HERE) call self%ps_bin%set([i,j,1], 0.)
                         cycle ! Do not consider white pixels detected after outside frequency range [LOW_LIM,UP_LIM]
                    endif
                    a = (1.-real(sh)/real(ULim_Findex))
                    denom = denom + a
                    self%score = self%score + a*imat_sz(i,j,1)/(2.*PI*sh)
                endif
            enddo
        enddo
        if(DEBUG_HERE) call self%ps_bin%write(trim(self%fbody)//'_binarized_polished.mrc')
        !normalization
        if(abs(denom) > TINY) then
            self%score = self%score/denom
        else
            THROW_HARD('Denominator = 0! calc_weighted_avg_sz_ccs')
        endif
        print *, 'SCORE =       ', self%score
    contains

        ! This subroutine is meant to generate a matrix in which in each pixel is
        ! stored the size of the correspondent connected component.
        subroutine generate_mat_sz_ccs(imat_in,imat_out,sz)
            integer, intent(in)  :: imat_in (:,:,:)
            integer, intent(out) :: imat_out(:,:,:)
            integer, intent(in)  :: sz(:)
            integer :: i, j, k
            if(size(imat_in) .ne. size(imat_out)) THROW_HARD('Input and Output matrices have to have the same dim!')
            do i = 1,size(imat_in, dim = 1)
                do j = 1,size(imat_in, dim = 2)
                    do k = 1,size(imat_in, dim = 3)
                        if(imat_in(i,j,k) > 0) then
                            where(imat_in == imat_in(i,j,k)) imat_out = sz(imat_in(i,j,k))
                        endif
                    enddo
                enddo
            enddo
        end subroutine generate_mat_sz_ccs
    end subroutine calc_weighted_avg_sz_ccs

    !This is the core subroutine of this module. Preprocessing and binarization of the
    !ps is performed here. Then estimation of how far do the detected rings extend.
    ! Finally estimation of the number of countiguous rings.
    subroutine process_ps(self)
      use gnufor2
      class(pspec_statistics), intent(inout) :: self
      integer, allocatable :: sz(:)
      real                 :: res
      integer              :: ind(1)
      integer              :: h, k, sh
      integer              :: i, j
      call build_resolutions_vector(self, BOX)
      ! calculation of entropy to identify fallacious data
      call self%per_shell_entropy()
      call prepare_ps(self)
      call binarize_ps(self)
      ! calculation of the weighted average size of the ccs ->  score
      call self%ps_bin%find_connected_comps(self%ps_ccs)
      sz = self%ps_ccs%size_connected_comps()
      call self%calc_weighted_avg_sz_ccs(sz)
      ! connected components size histogram
      !call hist(real(sz),N_BINS)
      ! Matlab compatible file for histograms
      ! open(119, file=trim(self%fbody)//'CCsSizeHist')
      ! write (119,*) 'sz=[...'
      ! do i = 1, size(sz)
      !     write (119,'(A)', advance='no') trim(int2str(sz(i)))
      !     if(i < size(sz)) write (119,'(A)', advance='no') ', '
      ! end do
      ! write (119,*) '];'
      ! close(119)
      deallocate(sz)
      ! call self%print_info()
      !if(self%fallacious) write(unit = 17, fmt = '(a)')  'Micrograph is fallacious '
      !close(17, status = "keep")
      ! printing shell rings on ps image, for developing and debugging
      ! Representation of the division in shells on the ps image
      do i = 1, BOX
          do j = 1, BOX
              h   = -int(BOX/2) + i - 1
              k   = -int(BOX/2) + j - 1
              sh  =  nint(hyp(real(h),real(k)))
              do k =1,size(self%res_vec)
                  if(abs(real(sh)-self%res_vec(k))<1) call self%ps%set([i,j,1], real(BOX/3))
              enddo
        enddo
    enddo
    call self%ps%write(trim(self%fbody)//'_shells.mrc')
  contains

      !Preprocessing steps: -) elimination of central spot
      !??? to complete according to what we decide
      subroutine prepare_ps(self)
          use simple_tvfilter, only: tvfilter
          class(pspec_statistics), intent(inout) :: self
          type(tvfilter) :: tvf
          integer        :: winsz
          call manage_central_spot(self)
          call self%ps%scale_pixels([1.,real(N_BINS)]) ! to set facilitate entropy calculation
          ! call tvf%new()
          ! call tvf%apply_filter(self%ps,LAMBDA)
          ! if(DEBUG_HERE) call self%ps%write(trim(self%fbody)//'_tvfiltered_ps.mrcs')
          ! call tvf%kill
          ! median filtering
          ! winsz = 2
          ! call self%ps%real_space_filter(winsz,'median')
          ! if(DEBUG_HERE) call self%ps%write(trim(self%fbody)//'_medianfiltered_ps.mrcs')
      end subroutine prepare_ps

      ! This subroutine extract a circle from the central part of the
      ! power spectrum, until LOW_LIM A and replace these gray values with
      ! the weighted average of the gray values in the ps outside that circle.
      ! The weight is the dist from the center. It's done to make it
      ! smoother.
      ! It is meant to reduce the influence of the central pixel.
      subroutine manage_central_spot(self)
          class(pspec_statistics), intent(inout) :: self
          real, allocatable :: rmat(:,:,:), rmat_dist(:,:,:)
          real              :: avg
          integer           :: lims(3,2), mh, mk
          integer           :: LLim_Findex
          integer           :: i, j
          integer           :: k, h, sh
          integer           :: cnt
          logical           :: lmsk(BOX,BOX,1)
          lmsk = .true.
          allocate(rmat     (BOX,BOX,1), source = self%ps%get_rmat())
          allocate(rmat_dist(BOX,BOX,1), source = 0.)
          LLim_Findex = calc_fourier_index(LOW_LIM,BOX,self%smpd)
          lims = self%ps%loop_lims(3)
          mh   = abs(lims(1,1))
          mk   = abs(lims(2,1))
          avg  = 0. !initialisations
          cnt  = 0
          do h=lims(1,1),lims(1,2)
              do k=lims(2,1),lims(2,2)
                  sh  = nint(hyp(real(h),real(k))) !shell to which px
                  i = min(max(1,h+mh+1),BOX)
                  j = min(max(1,k+mk+1),BOX)
                  if(sh < LLim_Findex ) then
                      lmsk(i,j,1) = .false.
                  elseif(sh >= LLim_Findex .and. sh < LLim_Findex + 2) then
                      avg = avg + rmat(i,j,1)
                      cnt = cnt + 1
                  endif
              enddo
          enddo
          avg = avg / real(cnt)
          do h=lims(1,1),lims(1,2)
              do k=lims(2,1),lims(2,2)
                  sh  = nint(hyp(real(h),real(k))) !shell to which px
                  i = min(max(1,h+mh+1),BOX)
                  j = min(max(1,k+mk+1),BOX)
                  if(sh < LLim_Findex ) then
                      rmat_dist(i,j,1) = sqrt(real((i-BOX/2-1)**2+(j-BOX/2-1)**2)) !dist between [i,j,1] and [BOX/2,BOX/2,1] (distance from the center)
                  endif
              enddo
          enddo
          rmat_dist = rmat_dist/maxval(rmat_dist)
          where( .not.lmsk ) rmat = avg*rmat_dist
          call self%ps%set_rmat(rmat)
      end subroutine manage_central_spot

      ! This subroutine performs binarization of a ps image
      ! using Canny edge detection with automatic threshold
      ! selection.
      subroutine binarize_ps(self)
          use simple_segmentation, only : canny
          class(pspec_statistics), intent(inout) :: self
          real,    allocatable :: rmat_bin(:,:,:)
          integer     :: ldim(3),i
          real        :: scale_range(2)
          scale_range = [1.,real(BOX)]
          call canny(self%ps,self%ps_bin,scale_range)
          do i=1,BOX  !get rid of border effects
              call self%ps%set([i,1,1],0.)
              call self%ps%set([i,BOX,1],0.)
              call self%ps%set([1,i,1],0.)
              call self%ps%set([BOX,i,1],0.)
          enddo
          if(DEBUG_HERE) call self%ps_bin%write(trim(self%fbody)//'_binarized.mrc')
          self%fallacious =  self%empty() !check if the mic is fallacious
          if(self%fallacious) write(logfhandle,*) trim(self%fbody), ' to be discarded'
      end subroutine binarize_ps
  end subroutine process_ps

  subroutine get_score(self,score)
      class(pspec_statistics), intent(inout) :: self
      real,                    intent(out)   :: score
      score = self%score
  end subroutine get_score

  subroutine run(self)
      class(pspec_statistics), intent(inout) :: self
      call self%p%new(self%cline)
      write(logfhandle,*) '******> Processing PS ', trim(self%fbody)
      call process_ps(self)
  end subroutine run

  subroutine kill_pspec_statistics(self)
      class(pspec_statistics), intent(inout) :: self
      call self%mic%kill()
      call self%ps%kill()
      call self%ps_bin%kill()
      call self%ps_ccs%kill()
      if(allocated(self%res_vec)) deallocate(self%res_vec)
      self%micname      = ''
      self%fbody        = ''
      self%ldim         = 0
      self%smpd         = 0.
      self%nn_shells    = 0
      self%fallacious   = .false.
  end subroutine kill_pspec_statistics
end module simple_genpspec_and_statistics

! SUBROUTINES MIGHT BE USEFUL
! subroutine log_ps(self,img_out)
!     class(pspec_statistics), intent(inout) :: self
!     type(image), optional :: img_out
!     real, allocatable :: rmat(:,:,:)
!     ! call self%ps%scale_pixels([1.,real(BOX)*2.+1.]) !not to have problems in calculating the log
!     rmat = self%ps%get_rmat()
!     rmat = log(rmat)
!     if(present(img_out)) then
!         call img_out%new(self%ldim, self%smpd)
!         call img_out%set_rmat(rmat)
!     else
!         call self%ps%set_rmat(rmat)
!     endif
!     deallocate(rmat)
! end subroutine log_ps
