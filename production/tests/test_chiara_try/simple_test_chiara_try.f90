module simple_test_chiara_try_mod
    use simple_powerspec_analysis
    use simple_image
    use simple_math
    use simple_segmentation
    use simple_picker_chiara
    implicit none
        real, allocatable :: centers(:,:)
        integer :: counter
    contains
! This subroutine performs laplacian filtering on the input image.
subroutine laplacian_filt(self)
    type(image), intent(inout) :: self
    real    :: k3(3,3,3), k2(3,3) !laplacian kernels (2D-3D)
    integer :: ldim(3)
    ldim = self%get_ldim()
    k2 = (1./8.)*reshape([0.,1.,0.,1.,-4., 1., 0., 1., 0.], [3,3])
    k3 = (1./12.)*reshape([0.,0.,0., 0.,1.,0., 0.,0.,0.,&
    &                     0.,1.,0., 1.,-6.,1., 0.,1.,0.,0.,0.,0., 0.,1.,0., 0.,0.,0.], [3,3,3])
    if(ldim(3) .ne. 1) then
        call self%imfilter(k3)
    else
        call self%imfilter(k2)
    endif
end subroutine laplacian_filt

  ! For Canny3D visit the following websites
  ! https://au.mathworks.com/matlabcentral/fileexchange/46260-3d-differential-canny-edge-detector
  ! https://en.wikipedia.org/wiki/Edge_detection#Second-order_approaches
end module simple_test_chiara_try_mod

program simple_test_chiara_try
    !$ use omp_lib
    !$ use omp_lib_kinds
  include 'simple_lib.f08'
  use simple_test_chiara_try_mod
  use simple_powerspec_analysis
  use gnufor2
  use simple_ctf
  use simple_micops
  use simple_image
  use simple_stackops
  use simple_math
  use simple_picker_chiara
  use simple_segmentation
  use simple_parameters, only: parameters
  use simple_cmdline,    only: cmdline
  use simple_tvfilter
  type(image)       :: img1, img2, pspec_img1, pspec_img2
  real, allocatable :: rmat(:,:,:), rmat_t(:,:,:)
  integer :: i,j, ldim(3), nptcls, box,n_vol
  type(ctf)       :: tfun
  type(ctfparams) :: ctfparms
  real :: smpd
  integer :: h, k, sh, cnt, px(3)
  real :: ave, sdev, maxv, minv, SumSQR
  real :: thresh(1)
  real, allocatable :: x(:)
  integer, allocatable :: sz(:,:)
  integer, allocatable :: pos(:,:), imat(:,:,:)
  real :: dist, ratio, corr_real, corr_ft
  type(tvfilter) :: tvf
  real,    allocatable :: xhist(:) !discretization of the values
  integer, allocatable :: yhist(:)
  real :: m(1)
  integer :: npxls_at_mode
  real    :: stretch_lim(2)
 !call process_ps_stack('pspecs_saga_polii.mrc', 'analisedSAGA.mrc', 1.14, 35., 1, 10) !winsz = 2
!call process_ps_stack('pspecs_saga_polii.mrc',     'saga_analysis_TVdenoising.mrc', 1.14, 50., 1, 10)
!call process_ps_stack('pspecs_sphire_tstdat.mrc', 'sphire_analysis_TVdenoising.mrc', 1.41, 20.,1, 10)

call img1%new([1854,1918,1], 1.)
call img1%read('AnisoResIteration1.mrc')
call img2%new([1854,1918,1], 1.)
call img2%read('AnisoResIteration2.mrc')
corr_real = img1%real_corr(img2)
print *, 'REAL SPACE CORRELATION ', corr_real
pspec_img1 = img1%mic2spec(512, 'sqrt', LP_PSPEC_BACKGR_SUBTR)
call pspec_img1%write('PowerSpectrum1.mrc')
pspec_img2 = img2%mic2spec(512, 'sqrt', LP_PSPEC_BACKGR_SUBTR)
call pspec_img2%write('PowerSpectrum2.mrc')
corr_ft =  pspec_img1%real_corr(pspec_img2)
print *, 'FT SPACE CORRELATION ', corr_ft

 end program simple_test_chiara_try
! !call find_ldim_nptcls('/home/chiara/Desktop/Chiara/ANTERGOS/forctf/0001_forctf.mrc', ldim, nptcls)

! ctfparms%smpd   = 1.32
! ctfparms%cs     = 2.7
! ctfparms%kv     = 300
! ctfparms%fraca  = 0.1
! ctfparms%dfx    = 2.62365627
! ctfparms%dfy    = 2.60851598
! ctfparms%angast = -29.8392296
! call find_ldim_nptcls('/home/chiara/Desktop/Chiara/ANTERGOS/forctf/0001_forctf.mrc', ldim, nptcls)
! call mic%new(ldim, ctfparms%smpd)
! tfun = ctf(ctfparms%smpd,ctfparms%kv,ctfparms%cs,ctfparms%fraca)


! matrix = reshape(real([ 1,1,1,0,0,6,5, &
!                  & 1,1,0,0,6,6,6, &
!                  & 1,0,0,2,0,6,0, &
!                  & 0,0,2,2,0,0,4, &
!                  & 0,5,0,0,0,4,4, &
!                  & 0,5,5,5,0,0,0, &
!                  & 0,5,5,0,0,3,3]),[7,7,1])
! call img%new([128,128,1], 1.)
! call img_cc%new([128,128,1], 1.)
! call img%ellipse([30,30], [11.,10.], 'no' )
! call img%ellipse([50,50], [11.,20.], 'no' )
