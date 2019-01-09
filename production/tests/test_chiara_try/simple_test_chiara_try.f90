module simple_test_chiara_try_mod
    use simple_powerspec_analysis
    use simple_image
    use simple_math
    use simple_segmentation
    implicit none

    contains
! This function returns a the laplacian matrix of the input image.
subroutine calc_laplacian(self, lap)
    type(image),                intent(inout) :: self
    real, allocatable,           intent(out)   :: lap(:,:,:)  !laplacian matrix
    type(image)        :: img_p                         !padded image
    real, allocatable  :: k(:,:,:)                      !Laplacian kernel
    real, allocatable  :: rmat_p(:,:,:)
    integer, parameter :: L = 3                         !dimension of the kernel
    integer            :: ldim(3)
    integer            :: i,j,m,n                       !loop indeces
    ldim = self%get_ldim()
    allocate(lap(ldim(1),ldim(2),1), k(-(L-1)/2:(L-1)/2,-(L-1)/2:(L-1)/2,1), source = 0.)
    k = (1./8.)*reshape([0,1,0,1,-4, 1, 0, 1, 0], [3,3,1]) ! Laplacian masks
    call img_p%new([ldim(1)+L-1,ldim(2)+L-1,1],1.)
    call self%pad(img_p)
    rmat_p = img_p%get_rmat()
    !$omp parallel do collapse(2) default(shared) private(i,j,m,n)&
    !$omp schedule(static) proc_bind(close)
    do i = 1, ldim(1)
      do j = 1, ldim(2)
          do m = -(L-1)/2,(L-1)/2
              do n = -(L-1)/2,(L-1)/2
                  lap(i,j,1) = lap(i,j,1)+rmat_p(i+m+1,j+n+1,1)*k(m,n,1)
              end do
          end do
      end do
    end do
    !omp end parallel do
    deallocate(k)
end subroutine calc_laplacian

subroutine laplacian_edge_detector(img)
    type(image), intent(inout) :: img
    type(image) :: img_lap
    real, allocatable :: lap(:,:,:), rmat(:,:,:), rmat_p(:,:,:), stdev_px(:,:,:)
    integer :: ldim(3), ldim_p(3)
    real    :: pixels(25)
    real :: thres, smpd
    integer :: i, j
    real ::  ave, sdev, maxv, minv
    call calc_laplacian(img,lap)
    ldim = img%get_ldim()
    call img_lap%new([ldim(1),ldim(2),ldim(3)],1.)
    call img_lap%set_rmat(lap)
    call img_lap%write('LaplacianImage.mrc')
    allocate(rmat(ldim(1), ldim(2), ldim(3)), source = 1.)
    call  img_lap%stats( ave, sdev, maxv, minv )
    print *, 'selected ', abs(minv+maxv)
    where(abs(lap) < abs(minv+maxv)) rmat = 0.
    call img%set_rmat(rmat)
    print *, 'minval ', minval(lap), 'maxval ', maxval(lap)
    ! smpd = img%get_smpd()
    ! ldim_p = ldim + 4
    ! ldim_p(3) = 1
    ! call img_p%new([ldim_p(1),ldim_p(2),ldim_p(3)], smpd)
    ! call img%pad(img_p)
    ! rmat_p = img_p%get_rmat()
    ! allocate(rmat(ldim(1), ldim(2), ldim(3)), source = 0.)
    ! allocate(stdev_px(ldim(1), ldim(2), ldim(3)), source = 0.)
    ! where(abs(lap)<TINY) rmat = 1.
    ! minv = 1. !initializaton
    ! maxv = 0.
    ! do i = 1, ldim(1)
    !     do j = 1, ldim(2)
    !         if(abs(lap(i,j,1))< TINY) then !just for white pxls
    !             pixels = reshape(rmat_p(i:i+4,j:j+4,1), [25])
    !             stdev_px(i,j,1) = stdev(pixels)
    !             rmat(i,j,1) = 1.
    !             if(stdev(pixels) > TINY) then
    !                 if(stdev_px(i,j,1) < minv) minv = stdev_px(i,j,1)
    !                 if(stdev_px(i,j,1) > maxv) maxv = stdev_px(i,j,1)
    !             endif
    !         endif
    !     enddo
    ! enddo
    ! thres = (minv + maxv)/2.
    ! !print *, 'thresh = ', thres
    ! where(stdev_px<thres) rmat = 0.
    ! call img%set_rmat(rmat)
    ! deallocate(lap, rmat, rmat_p, stdev_px)
end subroutine laplacian_edge_detector
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
  type(image)       :: img, img_try
  real, allocatable :: rmat(:,:,:), rmat_t(:,:,:), lap(:,:,:)
  integer, allocatable :: rmat_masked(:,:,:)
  integer :: i, ldim(3), nptcls, box, nthr
  real :: label_mirror
  type(ctf)       :: tfun
  type(ctfparams) :: ctfparms
  logical :: yes_no, discard
  real :: smpd
  integer :: points(3,2), px(2), j
  integer :: center(2)
  real :: radius
  logical :: outside
  type(image) :: img_win, img_templ
  real :: r
  integer :: h, k, m(3), sh, cnt, sz
  real :: ave, sdev, maxv, minv
  real :: thresh(1)
  real, allocatable :: grad3D(:,:,:), grad(:,:,:)

  call img%new([160,160,160], 1.)
  thresh(1) = 0.0034
  call img%read('/home/chiara/Desktop/Chiara/Segmentation/SObel3D/particle1.mrc')
  call sobel(img,thresh(1))
  call img%write('SobelPARTICLE1onow.mrc')

  ! call img%new([2048,2048,1],1.)
  ! call img%read('/home/chiara/Desktop/Chiara/Segmentation/ImgNOnoise.mrc')
  ! call laplacian_edge_detector(img)
  ! call img%write('LaplacianEdgeParticlesNOnoiseAvg.mrc')
  ! call img%new([1024,1024,1],1.)
  ! call img%read('/home/chiara/Desktop/Chiara/Segmentation/ImgNOISE.mrc')
  ! call laplacian_edge_detector(img)
  ! call img%write('LaplacianEdgeParticlesNOISEAvg.mrc')

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
