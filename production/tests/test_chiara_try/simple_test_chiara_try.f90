module simple_test_chiara_try_mod
    include 'simple_lib.f08'
    ! use simple_aff_prop
    ! use gnufor2
    ! use simple_ctf
    ! use simple_micops
    use simple_image, only : image
    ! use simple_stackops
    ! use simple_math
    ! use simple_segmentation
    ! use simple_parameters, only: parameters
    ! use simple_cmdline,    only: cmdline
    ! use simple_tvfilter
    ! use simple_ctf
    ! use simple_ppca
    ! use simple_stat
    ! use simple_lapackblas, only : sgeev
    implicit none
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

  subroutine circumference(img, center, rad, full)
      type(image),       intent(inout) :: img
      integer,           intent(in)    :: center(2)
      real,              intent(in)    :: rad
      logical, optional, intent(in)    :: full ! Default is false
      integer :: i, j, sh, h, k
      integer :: ldim(3)
      real    :: smpd
      logical :: ffull
      ffull = .false.
      if(present(full)) ffull = full
      ldim = img%get_ldim()
      smpd = img%get_smpd()
      do i = 1, ldim(1)
          do j = 1, ldim(2)
              h   =  i - 1 - center(1)
              k   =  j - 1 - center(2)
              sh  =  nint(hyp(real(h),real(k)))
              if(.not. ffull) then
                  if(abs(real(sh)-rad)<1) call img%set([i,j,1], 1.)
              else
                  if(abs(real(sh))<rad) call img%set([i,j,1], 1.)
              endif
            enddo
        enddo
    end subroutine circumference

    ! subroutine exec_symmetry_test_try(lp, msk)
    !     use simple_symanalyzer
    !     use simple_cmdline
    !     use simple_parameters
    !     use simple_builder
    !     real, intent(in) :: lp
    !     real, intent(in) :: msk ! radius in pixels
    !     type(parameters)      :: params
    !     type(image)           :: vol
    !     type(cmdline)         :: cline
    !     type(builder)         :: build
    !     character(len=STDLEN) :: fbody
    !     character(len=3)      :: pgrp
    !     real                  :: shvec(3), scale, smpd
    !     integer               :: ldim(3)
    !     integer, parameter    :: MAXBOX = 128
    !     !call cline%new()
    !     call cline%set('mkdir',  'yes')
    !     call cline%set('cenlp',     1.)
    !     call cline%set('center', 'yes')
    !     call cline%set('vol1', 'AutoSymmDetect/Convoluted.mrc')
    !     call cline%set('smpd', 0.358)
    !     call vol%new(self%ldim,self%smpd)
    !     call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
    !     call build%vol%read(params%vols(1))
    !     ! possible downscaling of input vol
    !     ldim = build%vol%get_ldim()
    !     scale = 1.
    !     params%msk = msk
    !     params%lp = lp
    !     if( ldim(1) > MAXBOX )then
    !         scale = real(MAXBOX) / real(ldim(1))
    !         call build%vol%fft
    !         call build%vol%clip_inplace([MAXBOX,MAXBOX,MAXBOX])
    !         call build%vol%ifft
    !         smpd         = build%vol%get_smpd()
    !         print *, 'smpd = ', smpd
    !         params%msk   = round2even(scale * params%msk)
    !         params%inner = round2even(scale * params%inner)
    !         params%width = scale * params%width
    !     endif
    !     ! low-pass limit safety
    !     params%lp = max(2. * smpd, params%lp)
    !     ! centering
    !     shvec = 0.
    !     if( params%center.eq.'yes' )then
    !         shvec = build%vol%calc_shiftcen(params%cenlp,params%msk)
    !         call build%vol%shift(shvec)
    !         ! fbody = get_fbody(params%vols(1),fname2ext(params%vols(1)))
    !         ! call build%vol%write(trim(fbody)//'_centered.mrc')
    !     endif
    !     ! mask volume
    !     if( params_glob%l_innermsk )then
    !         call build%vol%mask(params%msk, 'soft', inner=params%inner, width=params%width)
    !     else
    !         call build%vol%mask(params%msk, 'soft')
    !     endif
    !     ! run test
    !     print *, 'DEBUGGING'
    !     print *, 'params%msk = ', params%msk
    !     print *, 'params%lp  = ', params%lp
    !     print *, 'params%hp  = ', params%hp
    !     print *, 'params%cn_stop = ', params%cn_stop
    !     print *, 'params%platonic = ', params%platonic
    !     print *, 'pgrp ', pgrp
    !     call symmetry_tester(build%vol, params%msk, params%hp,&
    !     &params%lp, params%cn_stop, params%platonic .eq. 'yes', pgrp)
    !     ! end gracefully
    !     call simple_end('**** SIMPLE_SYMMETRY_TEST NORMAL STOP ****')
    ! end subroutine exec_symmetry_test_try

    subroutine generate_distribution(fname_coords_pdb, ldim, smpd)
        use simple_atoms, only: atoms
           character(len=*), intent(in) :: fname_coords_pdb
           integer,          intent(in) :: ldim(3)
           real,             intent(in) :: smpd
           type(atoms) :: atom
           type(image) :: vol !simulated distribution
           real        :: cutoff
           ! Generate distribution based on atomic position
           call vol%new(ldim,smpd)
           cutoff = 8.*smpd
           call atom%new(fname_coords_pdb)
           call atom%convolve(vol, cutoff)
           call vol%write('Convoluted.mrc')
   end subroutine generate_distribution


 ! My implementation of the Flood fill algorithm
 ! Soille, P., Morphological Image Analysis: Principles and Applications, Springer-Verlag, 1999, pp. 173–174.
 ! It's implemented for 8-neighbours, should I put the option for
 ! 4 neigh? To test.
   subroutine floodfill(img, px, tvalue, nb_modif)
     type(image), intent(inout) :: img
     integer,     intent(in)    :: px(3)
     integer,     intent(in)    :: tvalue ! target value
     integer,     intent(out)   :: nb_modif ! number of 8-neighbours which value has been modified
     integer :: pvalue ! pixel color
     integer :: ldim(3)
     integer :: i, nsz
     integer :: neigh_8(3,8,1) ! indeces of the 8 neigh
     integer :: neigh_px(3)
     nb_modif = 0
     ldim = img%get_ldim()
     if(ldim(3) > 1)  return!THROW_HARD('Not implemented for volumes! floodfill')
     pvalue = img%get([px(1),px(2),1]) ! For 2d images
     if (pvalue == tvalue) return
     ! set pixel value to tvalue
     call img%set(px(1:3),real(tvalue))
     call img%calc_neigh_8(px, neigh_8, nsz)
     do i = 1, nsz
         neigh_px(1:3) = neigh_8(1:3,i,1)
         if(nint(img%get(neigh_px(1:3))) /= tvalue ) then
             nb_modif = nb_modif + 1
             call img%set(neigh_px(1:3),real(tvalue))
         endif
     enddo
   end subroutine floodfill

   ! ! This subroutine calculates the diamenter of the
   ! ! connected component labelled n_cc in the connected
   ! ! component image img_cc
   ! subroutine diameter_cc(img_cc, n_cc, diam)
   !     type(image),  intent(inout) :: img_cc
   !     integer,      intent(in)    :: n_cc
   !     real,         intent(out)   :: diam
   !     integer, allocatable :: pos(:,:)         !position of the pixels of a fixed cc
   !     integer, allocatable :: imat_cc(:,:,:)
   !     logical, allocatable :: msk(:) ! For using function pixels_dist
   !     real  :: center_of_mass(3) ! geometrical center of mass
   !     real  :: radius
   !     imat_cc = int(img_cc%get_rmat())
   !     where(imat_cc .ne. n_cc) imat_cc = 0
   !     ! if(.not. any(imat_cc(:,:,:)) > 0) THROW_HARD('Inputted non-existent cc')
   !     ! Find center of mass of the cc
   !     call get_pixel_pos(imat_cc,pos)
   !     center_of_mass(1) = sum(pos(1,:))/real(size(pos,dim = 2))
   !     center_of_mass(2) = sum(pos(2,:))/real(size(pos,dim = 2))
   !     center_of_mass(3) = 1.
   !     if(allocated(imat_cc)) deallocate(imat_cc)
   !     allocate(msk(size(pos, dim =2)), source = .true.)
   !     ! Calculate maximim radius
   !     radius = pixels_dist(center_of_mass,real(pos),'max',msk)
   !     ! Return diameter
   !     diam = 2.*radius
   !     if(allocated(msk)) deallocate(msk)
   !     if(allocated(pos)) deallocate(pos)
   ! end subroutine diameter_cc
end module simple_test_chiara_try_mod

    program simple_test_chiara_try
       include 'simple_lib.f08'
       use simple_math
       use simple_segmentation
       use simple_test_chiara_try_mod
       use simple_image, only : image
       type(image) :: img, img_out
       integer     :: nb_modif
       integer     :: ldim(3)
       integer, allocatable :: imat_cc(:,:,:)
       real, allocatable :: rmat(:,:,:)
       real :: diam
       integer :: i
       real :: thresh(3)


       call img%new([160,160,1],0.358)
       do i = 1,887
           call progress(i,887)
           call img%read('chunk_avgs.mrcs',i)
           call otsu_robust_fast(img, img_out, stk=.true., nano=.false.,thresh=thresh)
           call img_out%write('otsuFastRobust.mrcs', i)
       enddo

 end program simple_test_chiara_try

     !BORDER EFFECTS IN PHASECORRELATION EXPLAINATION
     ! call img1%new([512,512,1],1.)
     ! call img2%new([512,512,1],1.)
     ! allocate(rmat(512,512,1), source=0.)
     ! rmat(256-50:256+50, 256-10:256+10,1) = 1.
     ! rmat(10:50,50:80,1) = 1.
     ! call img1%set_rmat(rmat)
     ! call img1%write('RectangleIMG.mrc')
     ! call img2%gauimg2D(50.,10.)
     ! call img2%write('GaussianIMG.mrc')
     ! pc = img1%phase_corr(img2,1)
     ! call pc%write('PhaseCorrelationIMG.mrc')
     ! VOLUM CLIUPPING FOR CHANGING ITS SMPD
     ! call find_ldim_nptcls('EMD-6287.mrc', vol_dim, nptcls, vol_smpd)
     ! call img%new(vol_dim,vol_smpd)
     ! call img%read('EMD-6287.mrc')
     ! print *, 'original dim = ', vol_dim, 'original smpd', vol_smpd
     ! shrink_factor = 2.64/vol_smpd
     ! print *, 'shrink_factor = ', shrink_factor
     ! ldim_clip(1) = round2even(real(vol_dim(1))/shrink_factor)
     ! ldim_clip(2) = round2even(real(vol_dim(2))/shrink_factor)
     ! ldim_clip(3) = round2even(real(vol_dim(3))/shrink_factor)
     ! print *, 'ldim_clip = ', ldim_clip
     ! call clip_vol('EMD-6287.mrc', 'EMD-6287_clip.mrc', ldim_clip, vol_smpd )
     ! call find_ldim_nptcls('EMD-6287_clip.mrc', ldim_clip, nptcls, clip_smpd)
     ! print *, 'Now smpd = ',clip_smpd

     ! CURVATURE ESTIMATION ON ARTOFICIAL DATA
     ! BOX = 256
     ! call img%new([BOX,BOX,1],1.)
     ! print * , '************************************************'
     ! print *, 'Curvature estimation on a circle'
     ! call img%ellipse([BOX/2,BOX/2],[22.,22.], 'yes')
     ! c = estimate_curvature(img,1)
     ! print * , '************************************************'
     ! print *, 'Curvature estimation on a arc of a circle'
     ! imat = nint(img%get_rmat())
     ! do i = 1,BOX
     !     do j = 1,BOX
     !         if(i < j) imat(i,j,1) = 0
     !     enddo
     ! enddo
     ! call img%set_rmat(real(imat))
     ! c = estimate_curvature(img,1)
     ! print * , '************************************************'
     ! imat = 0
     ! call img%set_rmat(real(imat))
     ! print *, 'Curvature estimation on a ellipse'
     ! call img%ellipse([BOX/2,BOX/2],[22.,17.], 'yes')
     ! imat = nint(img%get_rmat())
     ! c = estimate_curvature(img,1)
     ! print * , '************************************************'
     ! print *, 'Curvature estimation on a arc of an ellipse'
     ! do i = 1,BOX
     !     do j = 1,BOX
     !         if(i < j) imat(i,j,1) = 0
     !     enddo
     ! enddo
     ! call img%set_rmat(real(imat))
     ! c = estimate_curvature(img,1)
     ! print * , '************************************************'
     ! print *, 'Curvature estimation on a segment, length 20 pxls'
     ! imat = 0
     ! imat(BOX/2,140:160,1) = 1
     ! call img%set_rmat(real(imat))
     ! call img%ellipse([BOX/2,BOX/2],[46.,57.], 'yes')
     ! call img%find_connected_comps(img_cc)
     ! c = estimate_curvature(img_cc,1)
     ! c = estimate_curvature(img_cc,2)

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
    !Check for convergence.
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
