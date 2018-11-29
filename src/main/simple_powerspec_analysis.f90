module simple_powerspec_analysis
include 'simple_lib.f08'
use simple_image,    only : image
use simple_stackops, only : binarize_stack
implicit none

public
!public :: process_ps_stack
!private
#include "simple_local_flags.inc"

real, allocatable :: res_vec(:,:)
integer           :: nn_shells

contains

    !img has to be the connected component (cc) image
    function is_symmetric(img, label) result(yes_no)
        type(image), intent(inout) :: img
        integer,     intent(in)    :: label
        logical :: yes_no
        integer :: sz_cc, ldim(3)
        integer, allocatable :: vector(:), vector_sym(:)
        real,    allocatable :: rmat(:,:,:), rmat_t(:,:,:)
        logical, allocatable :: mask(:,:,:)
        real, parameter      :: THRESH = 0.2
        rmat = img%get_rmat()
        ldim = img%get_ldim()
        allocate(mask  (ldim(1),ldim(2),1), source = .false.)
        allocate(rmat_t(ldim(1),ldim(2),1), source = 0.)
        rmat_t(:,:,1) = transpose(rmat(:,:,1))
        where(abs(rmat - real(label)) < TINY) mask = .true.
        sz_cc = count(abs(rmat-real(label)) < TINY)
        allocate(vector(sz_cc),vector_sym(sz_cc), source = 0)
        vector     = pack(rmat,  mask)
        vector_sym = pack(rmat_t,mask)
        yes_no = .false.
        if((norm2(real(vector-vector_sym)))/(norm2(real(vector))+norm2(real(vector_sym)))< THRESH) yes_no = .true.
        deallocate(rmat,rmat_t,vector,vector_sym)
    end function is_symmetric


    !This function  takes in input the box sz and the smpd
    !of a power spectra image and consequentially builds an
    !ice template (circular)
    function build_ice_template(box, smpd, winsz) result(img_templ)
        integer, intent(in)  :: box
        real,    intent(in)  :: smpd
        integer, optional, intent(out) :: winsz
        type(image) :: img_templ
        integer :: wwinsz
        real :: axis
        wwinsz = int(box*4/100) ! 4% of the dim of the original image
        if(present(winsz)) winsz = wwinsz
        if(img_templ%exists()) call img_templ%kill
        call img_templ%new([winsz,winsz,1], smpd)
        axis = real(winsz*2/10) !20% of winsz
        call img_templ%ellipse([winsz/2,winsz/2], [axis,axis], 'yes')
    end function build_ice_template

    subroutine find_ice(bin_img,smpd)
        type(image), intent(inout) :: bin_img
        real,        intent(in)    :: smpd
        real,    allocatable :: rmat(:,:,:)
        logical, allocatable :: mask(:,:,:)
        integer :: ldim(3), i, j, winsz
        type(image) :: img_templ, bin_img_win
        real :: corr_coeff
        logical :: outside
        outside = .false.
        rmat = bin_img%get_rmat()
        if(any(rmat > 1.001) .or. any(rmat < 0.)) THROW_HARD('Expected binary image in input; find_ice')
        ldim = bin_img%get_ldim()
        allocate(mask(ldim(1),ldim(2),ldim(3)), source = .true.)
        img_templ = build_ice_template(ldim(1), smpd, winsz)
        print *, 'WINSZ = ', winsz
        call bin_img_win%new([winsz,winsz,1],smpd)
        do i = 1, ldim(1)
            do j = 1, ldim(2)
                if(rmat(i,j,1) > 0.5 .and. mask(i,j,1)) then !just around white pixels and just where mask == true
                   call bin_img%window_slim( [i+winsz/2,j+winsz/2], winsz, bin_img_win, outside )
                   !!!!!!!!!!!!!
                   ! if(i==221 .and. j==146) then
                   !   call bin_img%window_slim( [i+winsz,j+winsz], winsz, bin_img_win, outside )
                   !   call bin_img_win%write('SelectedWindow1.mrc')
                   !   call bin_img%window_slim( [i-winsz,j-winsz], winsz, bin_img_win, outside )
                   !   call bin_img_win%write('SelectedWindow2.mrc')
                   !   call bin_img%window_slim( [i,j], winsz, bin_img_win, outside )
                   !   call bin_img_win%write('SelectedWindow3.mrc')
                   !   call bin_img%window_slim( [i-winsz,j], winsz, bin_img_win, outside )
                   !   call bin_img_win%write('SelectedWindow4.mrc')
                   ! endif
                   !!!!!!!!!!!!
                   if(.not. outside) then
                        corr_coeff =  img_templ%real_corr(bin_img)
                        ! print *, 'pixel', i, j, 'corr ', corr_coeff
                   else
                        corr_coeff = 0.
                   endif
                   if(corr_coeff > 0.75) then
                      mask(i:i+winsz,j:j+winsz,1) = .false. !do not check again in this zone
                      print *, 'Detected ICE: pixel', i,j
                    endif
                   !build function is_simmetric(bin_img, coord) !coordinates of the detected ice
                  ! if(is_simmetric(bin_img, [i,j,1])) print *, "ICE detected at ", calc_fourier_index(i, ldim(1), smpd), ' A'  ! ?? not sure
                endif
            enddo
        enddo
        deallocate(rmat,mask)
    end subroutine find_ice

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

    !This function uses the ratio (white pxls)/(black pxls) per shell
    !to estimate how far do the detected rings go
    function find_res(box, smpd) result(res)
      real, intent(in)    :: smpd
      integer, intent(in) :: box
      integer :: i, cnt
      real    :: res
      real, parameter :: MIN_PERC = 0.02 !percentage of white pxls to understand how far do the rings go
      if(.not. allocated(res_vec)) THROW_HARD ('You have to build the resolution vector first; find_res')
      res = 0.
      do i = 1, size(res_vec,dim=1)
          if(res_vec(i,1) < MIN_PERC) then
              res = calc_fourier_index(real(res_vec(i,2)), box, smpd)
              return
          endif
      enddo
    end function find_res

    !This subroutine is meant to discard empty power spectra images.
    !If after binarisation the # of white pixels detected in the central
    !zone of the image is less than 2% of the tot
    !# of central pixels, than it returns yes, otherwhise no.
    function discard_ps(img,ldim) result(yes_no)
        type(image), intent(in) :: img
        integer,     intent(in) :: ldim(3)
        logical :: yes_no
        real, allocatable :: rmat(:,:,:),rmat_central(:,:,:)
        yes_no = .false.
        rmat = img%get_rmat()
        if(any(rmat > 1.001) .or. any(rmat < 0.)) THROW_HARD('Expected binary image in input; discard_ps')
        allocate(rmat_central(ldim(1)/2+1,ldim(2)/2+1,1), source = 0.)
        rmat_central(:,:,:) = rmat( ldim(1)/2-ldim(1)/4 : ldim(1)/2+ldim(1)/4 , &
                                  & ldim(2)/2-ldim(2)/4 : ldim(2)/2+ldim(2)/4 , :)
        deallocate(rmat)
        if(count(rmat_central(:,:,:) > 0.5)< 2*ldim(1)*ldim(2)*ldim(3)/(2*2*100)) yes_no = .true.; return
        deallocate(rmat_central)
    end function discard_ps


    !This function takes in input the name of a stack of power spectra images (fname2process),
    !the name of a stack in which store the results (fname), the smpd and a low-pass parameter.
    !It binarises all the images in the stack and estimates how far do the detected rings go.
    subroutine process_ps_stack(fname2process, fname, smpd, lp,winsz, n_shells)
      use gnufor2
      use simple_stackops, only : prepare_stack
      character(len=*),  intent(in) :: fname2process, fname
      real,              intent(in) :: smpd, lp
      integer, optional, intent(in) :: winsz     !for median filtering
      integer, optional, intent(in) :: n_shells  !number of shells
      integer              :: wwinsz
      integer, allocatable :: counter(:)
      real, allocatable :: rmat(:,:,:)
      integer           :: n, ldim(3), box
      integer           :: sh, ind(2), nn_shells
      integer           :: h, k, i, j, n_image
      logical, allocatable :: mask(:,:)
      real                 :: res, ax
      type(image)          :: img
      character(len = 100) :: iom
      integer              :: status
      logical              :: discard
      call find_ldim_nptcls(fname2process, ldim, n)
      ldim(3) = 1
      box = ldim(1)
      nn_shells = 10
      if(present(n_shells)) nn_shells = n_shells
      wwinsz = 2 !default
      if(present(winsz)) wwinsz = winsz
      call img%new(ldim,smpd)
      call prepare_stack(fname2process, 'prepared_stack.mrc', smpd, lp,wwinsz)
      print *, '>>>>>>>>>>>>>STACK PREPARED SUCCESSFULLY>>>>>>>>>>>>>'
      call binarize_stack('prepared_stack.mrc','binarised_stack.mrc', smpd)
      print *, '>>>>>>>>>>>>>STACK BINARISED SUCCESSFULLY>>>>>>>>>>>>>'
      call build_resolutions_vector(box, smpd, res_vec,nn_shells)
      print *, "Power spectra divided into ", nn_shells, ' shells'
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
      write(unit = 17, fmt = "(a,tr1,i0.0)") 'Number of shells: ', nn_shells
      write(unit = 17, fmt = "(a,i0.0,a,i0.0)") 'Pixel split up in shells in the interval  0 -  ', &
                              & calc_fourier_index(res_vec(2,2), box, smpd),  &
                              & ' with step ', calc_fourier_index(res_vec(size(res_vec, dim = 1),2), box, smpd)
      write(unit = 17, fmt = '(a)') ''
      write(unit = 17, fmt = "(a)")  '-----------IMAGE ANALYSIS------------- '
      mask(:,1) = .false.
      mask(:,2) = .true.
      do n_image = 1, n
          call img%read('binarised_stack.mrc', n_image)
          discard = discard_ps(img,ldim)
          if(discard) write(unit = 17, fmt = '(a)') 'Empty micrograph ', n_image ; continue
          if(.not. discard) then
              rmat = img%get_rmat()
              counter = 0
              do i = 1, ldim(1)
                  do j = 1, ldim(2)
                      h   = -int(ldim(1)/2) + i - 1
                      k   = -int(ldim(2)/2) + j - 1
                      sh  = nint(hyp(real(h),real(k)))        !shell to which px (i,j) belongs
                      ind = minloc(abs(res_vec-sh),mask)      !corresponding shell in res_vec
                      counter(ind(1)) = counter(ind(1)) + 1   !Number of pixels per shell, it is fixed, I could calculate aside
                      if (rmat(i,j,1) > 0.5 .and. ind(1) < nn_shells) then !binary image, discard edges
                          res_vec(ind(1),1) = res_vec(ind(1),1) + 1. !update # of white pixel in the shell
                      endif
                  enddo
              enddo
              where(counter > 0) res_vec(:,1) = res_vec(:,1) / counter(:)  ! normalise
              write(unit = 17, fmt = "(a,tr1,i0)") 'Image', n_image-1
              res = find_res(box, smpd)
          endif
          call img%read(fname2process, n_image)
          if(.not. discard) ax = calc_fourier_index(res, box, smpd)
          if(.not. discard) call img%build_ellipse([real(ldim(1))/2., real(ldim(2))/2.],[ax,ax], 0.)
          call img%write(fname,n_image)
          if(.not. discard) write(unit = 17, fmt = "(a,i0.0,a)") 'Visible rings until res ', int(res), 'A'
      enddo
      close(17, status = "keep")
      call img%kill
      deallocate(mask, counter)
      if(allocated(rmat)) deallocate(rmat)
    end subroutine process_ps_stack
end module simple_powerspec_analysis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TESTS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 ! img = build_ice_template(512, 1.41, winsz)
 ! print *, 'WINSZ = ', winsz
 ! call img_win%new([int(winsz/2),int(winsz/2),1], 1.41)
 ! call img%write('IceTemplate.mrc')
 ! call img%window_slim([int(winsz/4),int(winsz/4)], 10, img_win, outside )
 ! call img_win%write('IceTemplateWin.mrc')
 ! print *, 'OUTSIDE = ', outside

 ! call img%ellipse([256,256],[20.,20.], 'yes')
 ! call img%ellipse([352,312],[10.,5.], 'yes')
 ! call img%ellipse([23,25],[5.,10.], 'yes')
 ! call img%ellipse([112,53],[10.,10.], 'yes')
 !  call img%ellipse([220,153],[8.,8.], 'yes')
 ! call img%read('ToyImage.mrc')
 ! call img%find_connected_comps(img_cc)
 ! call img_cc%write('ToyImageCC.mrc')
 ! yes_no = is_symmetric(img_cc, 1)
 ! print *, 'CC ', 1, ' is symmetric ', yes_no
 ! yes_no = is_symmetric(img_cc, 2)
 ! print *, 'CC ', 2, ' is symmetric ', yes_no
 ! yes_no = is_symmetric(img_cc, 3)
 ! print *, 'CC ', 3, ' is symmetric ', yes_no
 ! yes_no = is_symmetric(img_cc, 4)
 ! print *, 'CC ', 4, ' is symmetric ', yes_no
 ! yes_no = is_symmetric(img_cc, 5)
 ! print *, 'CC ', 5, ' is symmetric ', yes_no

  ! !TO START AGAIN FROM HERE
  ! call img%new([512,512,1],1.)
  ! call img_cc%new([512,512,1],1.)
  ! call img%read('SAGAWhithICEBin.mrc')
  ! call img%find_connected_comps(img_cc)
  ! call img_cc%write('SAGAWithICEBinCC.mrc')
  ! rmat = img_cc%get_rmat()
  ! do i = 1, int(maxval(rmat))
  !   yes_no = is_symmetric(img_cc, i)
  !   print *, 'CC ', i, ' is symmetric ', yes_no
  ! enddo

 ! call process_ps_stack('pspecs_sphire_tstdat.mrc', 'analysed_pspecs_sphire.mrc', 1.41, 20.,1)
! call process_ps_stack('pspecs_saga_polii.mrc', 'analysed_pspecs_saga_polii.mrc', 1.14, 35.,2)
