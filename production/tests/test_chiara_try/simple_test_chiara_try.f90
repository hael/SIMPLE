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
  type(image)       :: img, img_cc, img_bin
  real, allocatable :: rmat(:,:,:), rmat_t(:,:,:), lap(:,:,:)
  integer :: i,j, ldim(3), nptcls, box, nthr, n_vol
  type(ctf)       :: tfun
  type(ctfparams) :: ctfparms
  real :: smpd, seleted_t(1)
  integer :: h, k, sh, cnt
  real :: ave, sdev, maxv, minv, SumSQR
  real :: thresh
  real, allocatable :: x(:)
  integer, allocatable :: sz(:,:)
  real,    allocatable :: xhist(:) !discretization of the values
  integer, allocatable :: yhist(:) !number of occurences
  integer, allocatable :: pos(:,:), imat(:,:,:)
  real :: dist, diameter, ratio
  real, allocatable :: ratios(:), x_cc(:)
  logical, allocatable :: border(:,:,:)

  call img%new([4096,4096,1],1.)
  call img%read('/home/chiara/Desktop/Chiara/Denoising/TotalVariation/NegativeSTORIGINAL.mrc')

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

!FOR improved gradient calculation
! type(image) :: img, Gr_original, Gr_improved
! real, allocatable :: grad_original(:,:,:), grad_improved(:,:,:)
! real :: thresh(1)
! call img%new([128,128,1],1.)
! call Gr_original%new([128,128,1],1.)
! call Gr_improved%new([128,128,1],1.)
!
! call img%read('/home/chiara/Desktop/Chiara/Segmentation/EdgeDetection/one_projection.mrc')
! call img%calc_gradient(grad_original)
! call Gr_original%set_rmat(grad_original)
! call Gr_original%write('Grad_original_one_projection.mrc')
! call img%calc_gradient_improved(grad_improved)
! call Gr_improved%set_rmat(grad_improved)
! call Gr_improved%write('Grad_improved_one_projection.mrc')
! thresh(1) = 0.05
! call sobel(img, thresh)
! call img%write('Sobel_improved_thresh005.mrc')
