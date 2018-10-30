module simple_powerspec_analysis
include 'simple_lib.f08'
use simple_image,    only : image
use simple_stackops, only : calc_log_stack, binarize_stack
implicit none

public :: process_ps_stack
private
#include "simple_local_flags.inc"

real, allocatable :: res_vec(:,:)
integer           :: nn_shells

contains
    !This function builds a vector that splits up the images in n_shells shells
    !and stores the result in res_vec(:,2). In res_vec(:,1) will be saved the
    !number of white pixels in the shell corresponding to res_vec(:,2).
    subroutine build_resolutions_vector(box, smpd, res_vec, n_shells)
        integer, intent(in)  :: box
        real,    intent(in)  :: smpd
        real, allocatable, intent(out) :: res_vec(:,:)
        integer, optional :: n_shells
        integer :: step, nn_shells, i
        nn_shells = 10
        if(allocated(res_vec)) THROW_HARD ('You are trying to allocate again the resolution vector; build_resolutions_vector')
        if(present(n_shells) .and. n_shells < 2) THROW_HARD('Too low number of shells; build_resolutions_res_vec')
        if(present(n_shells)) nn_shells = n_shells
        allocate(res_vec(nn_shells,2), source = 0.)
        step = int(box/(2*(nn_shells-2)))
        do i = 2, nn_shells
            res_vec(i,2) = (i-1)*real(step)
        end do
    end subroutine build_resolutions_vector

    !This function calculates the ratio (white pxls)/(black pxls) per shell
    !and selects as a resolution the shell in which the ratio becomes lower
    !than thres.
    function find_res(box, smpd) result(res)
      real, intent(in)    :: smpd
      integer, intent(in) :: box
      real, PARAMETER     :: thres = 2*10.**(-2)
      integer :: i, cnt
      real    :: res
      if(.not. allocated(res_vec)) THROW_HARD ('You have to build the resolution vector first; find_res')
      cnt = 0
      res = 0.
      do i = 1, size(res_vec, dim = 1)
          cnt = cnt + 1
          if(res_vec(i, 1) < thres) then
              res = calc_fourier_index(real(res_vec(i,2)), box, smpd)
              return
          endif
      end do
      if(cnt == size(res_vec, dim = 1)) print *, 'Too low threshold for resolution estimation, find_res'
    end function find_res

    !This function takes in input the name of a stack of power spectra images (fname2process),
    !the name of a stack in which store the results (fname), the smpd and a low-pass parameter.
    !It binarises all the images in the stack and estimates how far do the rings go.
    subroutine process_ps_stack(fname2process, fname, smpd, lp)
      character(len=*), intent(in) :: fname2process, fname
      real,             intent(in) :: smpd, lp
      real, allocatable :: rmat(:,:,:)
      integer           :: n, ldim(3), box
      integer           :: sh, ind(2), nn_shells
      integer           :: h, k, i, j, n_image
      integer, allocatable :: counter(:)
      logical, allocatable :: mask(:,:)
      real                 :: res, ax
      type(image)          :: img
      character(len = 100) :: iom
      integer              :: status
      call find_ldim_nptcls(fname2process, ldim, n)
      ldim(3) = 1
      box = ldim(1)
      nn_shells = 8
      call img%new(ldim,smpd)
      call binarize_stack(fname2process, 'binarised_stack.mrc', smpd, lp)
      call build_resolutions_vector(box, smpd, res_vec,nn_shells)
      print *, "NN_SHELLS", nn_shells
      allocate(counter(nn_shells), mask(nn_shells,2))
      open(unit = 17, access = 'sequential', action = 'readwrite',file = "PowerSpectraAnalysis.txt", form = 'formatted', iomsg = iom, iostat = status, position = 'append', status = 'replace')
      write(unit = 17, fmt = '(a)') '>>>>>>>>>>>>>>>>>>>>POWER SPECTRA STATISTICS>>>>>>>>>>>>>>>>>>'
      write(unit = 17, fmt = '(a)') ''
      write(unit = 17, fmt = "(a,a)")  'Input stack ', fname2process
      write(unit = 17, fmt = "(a,a)")  'Output stack ', fname
      write(unit = 17, fmt = "(a,i0,tr1,i0,tr1,i0)") 'Dimensions ', ldim
      write(unit = 17, fmt = "(a,f0.2)")  'Smpd ', smpd
      write(unit = 17, fmt = "(a,i0)")  'N images  ', n
      write(unit = 17, fmt = '(a)') ''
      write(unit = 17, fmt = "(a)")  '-----------SELECTED PARAMETERS------------- '
      write(unit = 17, fmt = '(a)') ''
      write(unit = 17, fmt = "(a, f0.0)")  'Low pass filter ', lp
      !write(unit = 17, fmt = "(a, f0.0)")  'Threshold for resolution calculation ', thres
      write(unit = 17, fmt = "(a,tr1,i0.0)") 'Number of shells: ', nn_shells
      write(unit = 17, fmt = "(a,f0.0,a,f0.0)") 'Pixel split up in shells from 0 to : '&
                              &, res_vec(size(res_vec, dim = 1),2), 'with step ', res_vec(2,2)
      write(unit = 17, fmt = '(a)') ''
      write(unit = 17, fmt = "(a)")  '-----------IMAGE ANALYSIS------------- '
      mask(:,1) = .false.
      mask(:,2) = .true.
      do n_image = 1, n
          call img%read('binarised_stack.mrc', n_image)
          rmat = img%get_rmat()
          counter = 0
          do i = 1, ldim(1)
              do j = 1, ldim(2)
                  h   = -int(ldim(1)/2) + i - 1
                  k   = -int(ldim(2)/2) + j - 1
                  sh  = nint(hyp(real(h),real(k)))        !shell to which px (i,j) belongs
                  ind = minloc(abs(res_vec-sh),mask)      !corresponding shell in res_vec
                  counter(ind(1)) = counter(ind(1)) + 1   !Number of pixels per shell, it is fixed, I could calculate apart
                  if (rmat(i,j,1) > 0.5 .and. ind(1) < nn_shells) then !binary image, discard edges
                      res_vec(ind(1),1) = res_vec(ind(1),1) + 1. !update # of white pixel in the shell
                  endif
              enddo
          enddo
          where(counter > 0) res_vec(:,1) = res_vec(:,1) / counter(:)  ! normalise
          write(unit = 17, fmt = "(a,tr1,i0)") 'Image', n_image-1
          !call vis_mat(res_vec)
          res = find_res(box, smpd)
          call img%read(fname2process, n_image)
          ax = calc_fourier_index(res, box, smpd)
          call img%build_ellipse([real(ldim(1))/2., real(ldim(2))/2.],[ax,ax], 0.)
          call img%write(fname,n_image)
          !print *, 'Image ', n_image-1, 'res = ', res
          write(unit = 17, fmt = "(a,i0.0,a)") 'Visible rings until res ', int(res), 'A'
      enddo
      close(17, status = "keep")
      call img%kill
      deallocate(rmat, mask, counter)
    end subroutine process_ps_stack
end module simple_powerspec_analysis
