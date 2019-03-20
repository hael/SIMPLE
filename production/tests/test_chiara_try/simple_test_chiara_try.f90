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

  ! 3D line for regression for identification atom rows
  function fun_try(p,n) result(r)
          real,    intent(in) :: p(:)
          integer, intent(in) :: n
          real :: r(n)
          real :: x, y
          x = p(1)
          y = p(2)
          r(1) = 1
          r(2) = x*y
  end function fun_try

      SUBROUTINE PRINT_EIGENVALUES( DESC, N, WR, WI )
          CHARACTER(len= *)    :: DESC
          INTEGER ::  N
          REAL    ::  WR( : ), WI( : )
          REAL, parameter ::   ZERO = 0.0
          INTEGER ::  J
          WRITE(*,*)
          WRITE(*,*) DESC
          DO J = 1, N
             IF( WI( J ).EQ.ZERO ) THEN
                WRITE(*,9998,ADVANCE='NO') WR( J )
             ELSE
                WRITE(*,9999,ADVANCE='NO') WR( J ), WI( J )
             END IF
          END DO
          WRITE(*,*)
     9998 FORMAT( 11(:,1X,F6.2) )
     9999 FORMAT( 11(:,1X,'(',F6.2,',',F6.2,')') )
          RETURN
      END SUBROUTINE PRINT_EIGENVALUES

      SUBROUTINE PRINT_EIGENVECTORS( DESC, N, WI, V, LDV )
          CHARACTER(len = *) ::   DESC
          INTEGER            :: N, LDV
          REAL   :: WI( : ), V( :, : )
          REAL, parameter    :: ZERO = 0.0
          INTEGER   ::       I, J
          WRITE(*,*)
          WRITE(*,*) DESC
          DO I = 1, N
             J = 1
             DO WHILE( J.LE.N )
                IF( WI( J ).EQ.ZERO ) THEN
                   WRITE(*,9998,ADVANCE='NO') V( I, J )
                   J = J + 1
                ELSE
                   WRITE(*,9999,ADVANCE='NO') V( I, J ), V( I, J+1 )
                   WRITE(*,9999,ADVANCE='NO') V( I, J ), -V( I, J+1 )
                   J = J + 2
                END IF
             END DO
             WRITE(*,*)
          END DO
     9998 FORMAT( 11(:,1X,F6.2) )
     9999 FORMAT( 11(:,1X,'(',F6.2,',',F6.2,')') )
          RETURN
      END SUBROUTINE PRINT_EIGENVECTORS
end module simple_test_chiara_try_mod

program simple_test_chiara_try
    !$ use omp_lib
    !$ use omp_lib_kinds
  include 'simple_lib.f08'
  use simple_aff_prop
  use simple_commander_distr_wflows
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
  use simple_ctf
  use simple_ppca
  use simple_stat
  use simple_lapackblas, only : sgeev
  type(image)       :: img, img_cc
  real, allocatable :: rmat(:,:,:), rmat_t(:,:,:)
  integer :: i,j, ldim(3), nptcls,n_vol
  type(ctfparams) :: ctfparms
  real :: smpd
  integer :: h, k, sh, cnt, px(3)
  real :: ave, sdev, maxv, minv, SumSQR
  real :: thresh(1)
  real :: x(2,4), y(4), sig(4), v(2,2), w(2), chisq
  integer, allocatable :: sz(:)
  integer, allocatable :: pos(:,:), imat(:,:,:)
  real :: dist, ratio, corr_real, corr_ft
  type(tvfilter) :: tvf
  real,    allocatable :: xhist(:) !discretization of the values
  integer, allocatable :: yhist(:)
  real :: lambda
  integer :: npxls_at_mode
  real    :: stretch_lim(2)
  logical :: outside
  type(motion_correct_distr_commander) :: xmotion_correct_distr
  type(cmdline) :: cline
  type(ctf) :: tfun
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer              :: ldim_shrunken(3), n_images, D, recsz, ifeat, box_shrunken
  real,    allocatable :: avg(:), feat(:), dat(:), mat(:,:,:)
  logical, allocatable :: l_mask(:,:,:)
  integer, allocatable :: coords(:,:)
  integer, parameter   :: BOX = 256, OFFSET = BOX/4-1, BOFFSET = 20 !BOX/4-1
  real,    parameter   :: SHRINK = 1.
  type(ppca)           :: my_ppca
  type(image)          :: img_rev, img_msk, mic_denoised, mic_denoised_norm
  real                 :: mskrad, smpd_shrunken
  logical              :: err
  logical, allocatable :: border(:,:,:)
  !TRY BORDER MASK
  call img%new([50,50,50], 1.)
  rmat = img%get_rmat()
  rmat(20:30,20:30,20:30) = 1.
  call img%set_rmat(rmat)
  call img%write('Cube.mrc')
  print *, 'count = ', count(rmat > 0.)
  call img%border_mask(border,1)
  print *, 'count border = ', count(border .eqv. .true.)
  rmat = 0.
  where(border .eqv. .true.) rmat = 1.
  call img%set_rmat(rmat)
  call img%write('BorderCube.mrc')


  !eigenvectors and eigenvalues
  ! real    :: A(5,5)
  ! integer, parameter :: N = 5
  ! integer, parameter :: LDA = 5, LDVL = 5, LDVR = 5
  ! integer, parameter :: LWMAX = 1000
  ! integer :: INFO, LWORK
  ! real    :: matrix( 5, 5 ), VL(  5, 5  ), VR(  5, 5  ),    WR( 5 ), WI( 5 ), WORK(  1000 )
! ! EIGENVECTORS AND VALUES CALCULATION WITH LAPACK see
! ! http://www.netlib.org/lapack/explore-html/d3/dfb/group__real_g_eeigen_ga104525b749278774f7b7f57195aa6798.html https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/sgeev_ex.f.htm
! ! https://github.com/joe-of-all-trades/regionprops3/blob/master/regionprops3.m
!  A = reshape([-1.01, 0.86, -4.6, 3.31,-4.81, &
!                  &  3.98, 0.53,-7.04, 5.29, 3.55, &
!                  &  3.30, 8.26,-3.89, 8.20,-1.51, &
!                  &  4.43, 4.96,-7.66,-7.33, 6.18, &
!                  &  7.31,-6.43,-6.16, 2.47, 5.58],[5,5])
!  matrix = transpose(matrix)
!  LWORK = -1
!  CALL SGEEV( 'Vectors', 'Vectors', N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
!  LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
! !Solve eigenproblem.
!  CALL SGEEV( 'Vectors', 'Vectors', N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
! !Check for convergence.
!  IF( INFO.GT.0 ) THEN
!      WRITE(*,*)'The algorithm failed to compute eigenvalues.'
!      STOP
!  END IF
!  !Print eigenvalues.
! CALL PRINT_EIGENVALUES( 'Eigenvalues', N, WR, WI )
! !Print left eigenvectors.
!  CALL PRINT_EIGENVECTORS( 'Left eigenvectors', N, WI, VL, LDVL )
! !Print right eigenvectors.
!  CALL PRINT_EIGENVECTORS( 'Right eigenvectors', N, WI, VR, LDVR )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!WORKING ON POWER SPECTRA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! call process_ps_stack('pspecs_saga_polii.mrc', 'analisedSAGA.mrc', 1.14, 35., 1, 10) !winsz = 2
!call process_ps_stack('pspecs_saga_polii.mrc',     'saga_analysis_TVdenoising.mrc', 1.14, 50., 1, 10)
!call process_ps_stack('pspecs_sphire_tstdat.mrc', 'sphire_analysis_TVdenoising.mrc', 1.41, 20.,1, 10)

!!!!!!!!TEST CORRELATION IN REAL-F SPACE. TO CORRECT
! call img1%new([1854,1918,1], 1.)
! call img1%read('AnisoResIteration1.mrc')
! call img2%new([1854,1918,1], 1.)
! call img2%read('AnisoResIteration2.mrc')
! corr_real = img1%real_corr(img2)
! print *, 'REAL SPACE CORRELATION ', corr_real
! pspec_img1 = img1%mic2spec(512, 'sqrt', LP_PSPEC_BACKGR_SUBTR)
! call pspec_img1%write('PowerSpectrum1.mrc')
! pspec_img2 = img2%mic2spec(512, 'sqrt', LP_PSPEC_BACKGR_SUBTR)
! call pspec_img2%write('PowerSpectrum2.mrc')
! corr_ft =  pspec_img1%real_corr(pspec_img2)
! print *, 'FT SPACE CORRELATION ', corr_ft
! !!!!!!! ARTIFICIAL MOVIE CREATION!!!!!!!!!!
!  call img%new([512,512,1],1.)
!  call img1%new([256,256,1],1.)
! do i = 1, 6
!  call img1%read('Dd'//int2str(i)//'.mrc')
!  call img1%add_gauran(0.8)
!  call img1%fft()
!  call img1%bp(0.,10.)
!  call img1%ifft()
!  call img1%write('StackImagesNoiseLP2objects.mrc',i)
! enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !call xmotion_correct_distr%execute(cline)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !FRAME DENOISING FOR FEPT NANOPARTICLES
 !NLmean denoising
 ! call find_ldim_nptcls('/home/chiara/Desktop/Chiara/mrcsFePtNanoparticles/BiMetal6_dm4_Hour_00_Minute_00_Second_00_Frame_0328.mrc', ldim, nptcls, smpd)
 ! call img%new(ldim, smpd)
 ! call img%read('/home/chiara/Desktop/Chiara/mrcsFePtNanoparticles/BiMetal6_dm4_Hour_00_Minute_00_Second_00_Frame_0328.mrc')
 ! call img%NLmean()
 ! call img%write('NLmean0328.mrc')
 !Total variation denosing
 ! call img%read('/home/chiara/Desktop/Chiara/mrcsFePtNanoparticles/BiMetal6_dm4_Hour_00_Minute_00_Second_00_Frame_0000.mrc')
 ! call tvf%new()
 !call raise_exception( nptcls, ldim, 'apply_tvf_imgfile' )
 ! do i = 1, 100
 !     lambda = 1+(real(i)/100.)
 !     call tvf%apply_filter(img, lambda)
 !     call img%write('TV0000L'//real2str(real(i)/100.)//'.mrc')
 ! enddo
 !!!!!!!!!!!!!!!!!PCA UNDERSTANDING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LLOOOOOK HERE https://sebastianraschka.com/Articles/2014_pca_step_by_step.html
 ! call my_ppca%new(size(centers, dim = 2),3,1)
 ! call my_ppca%master('vecs4ppca.bin', recsz ,'feat_stk.bin', 10)
 !
 ! subroutine master( self, datastk, recsz, featstk, maxpcaits )
 !     class(ppca), intent(inout)             :: self
 !     character(len=*), intent(in)           :: datastk, featstk
 !     integer, intent(in)                    :: recsz, maxpcaits
 !     integer                                :: k, file_stat, funit2, recsz2, err
 !     real                                   :: p, p_prev
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!affinity propagation
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
