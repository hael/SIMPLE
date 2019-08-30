module simple_test_chiara_try_mod
    include 'simple_lib.f08'
    ! use simple_aff_prop
    ! use simple_commander_distr_wflows
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

    subroutine calc_percentage_phiber(imgfile_jpeg)
         use simple_jpg
         character(len=*), intent(inout) :: imgfile_jpeg
         type(image) :: img
         type(image) :: img_aux ! TO REMOVE
         integer     :: ldim(3)
         integer     :: n_images
         integer     :: i, j
         integer     :: cnt
         real        :: smpd
         real        :: means(3)
         logical, allocatable :: msk(:,:,:)
         integer, allocatable :: labels(:)
         integer, allocatable :: imat(:,:,:)
         real,    allocatable :: x(:)
         real,    allocatable :: rmat(:,:,:)
         integer, parameter   :: MAXIT = 100
         call img%new([256,256,1],1.)
         ! call img%read_jpg(imgfile_jpeg)
         !call find_ldim_nptcls(imgfile_jpeg, ldim, n_images, smpd)
         ! call img%new(ldim, smpd)
         call img%write_jpg('Fuck.jpg')
         ! call img%read_jpg(imgfile_jpeg)
         call img%write('InputImgMrcFormat.mrc')
         rmat = img%get_rmat()
         allocate(msk(ldim(1),ldim(2),1), source = .true.)
         x = pack(rmat, msk)
         call classify_3objects(x,MAXIT, means,labels)
         call img_aux%new(ldim, smpd)
         cnt = 0
         do i = 1, ldim(1)
             do j = 1, ldim(2)
                 cnt = cnt + 1
                 rmat(j,i,1) = labels(cnt)
             enddo
         enddo
         call img_aux%set_rmat(rmat)
         call img_aux%write('LabelsImg.mrc')
         imat = nint(img_aux%get_rmat())
         print *, 'Minval(imat) = ', minval(imat), 'maxval(imat) = ', maxval(imat)
         print *, 'Void   = ', count(imat == 1), 'percent ',real(count(imat == 1))/real(size(imat))
         print *, 'Phiber = ', count(imat == 2), 'percent ',real(count(imat == 2))/real(size(imat))
         print *, 'Brown  =',  count(imat == 3), 'percent ',real(count(imat == 3))/real(size(imat))
       contains
         !implements the sortmeans algorithm
         subroutine classify_3objects( dat, maxits, means, labels )
             real,                 intent(in)  :: dat(:)    !< array for input
             integer,              intent(in)  :: maxits    !< limit sort
             real,                 intent(out) :: means(:)  !< array for output
             integer, allocatable, intent(out) :: labels(:) !< labels for output
             logical, allocatable :: mask(:)
             real, parameter :: TOLL = 0.1
             integer :: ncls, ndat, clssz, i, j, cnt_means, loc(1), changes
             ncls = size(means)
             ndat = size(dat)
             if( allocated(labels) ) deallocate(labels)
             allocate(  mask(ndat), labels(ndat), stat=alloc_stat )
             if(alloc_stat /= 0) call allocchk("sortmeans; simple_math", alloc_stat)
             ! initialization by sorting
             !In this way I don't need to sort the data
             mask = .true.
             means(1) = minval(dat, mask) !void
             where(abs(dat-means(1))<TINY) mask = .false.
             means(2) = maxval(dat, mask) !phiber
             where(abs(dat-means(2))<TINY) mask = .false.
             means(3) = sum(dat, mask)/real(count(mask))
             print *, 'means = ', means
             ! the kmeans step
             labels = 1
             do j=1,maxits
                 changes = 0
                 ! find closest
                 do i=1,ndat
                     loc = minloc((means-dat(i))**2.0)
                     if( labels(i) /= loc(1) ) changes = changes + 1
                     labels(i) = loc(1)
                 end do
                 ! update means
                 do i=1,ncls
                     where( labels == i )
                         mask = .true.
                     else where
                         mask = .false.
                     end where
                     means(i) = sum(dat, mask=mask)/real(count(mask))
                 end do
                 if( changes == 0 ) exit
             end do
             deallocate(mask)
         end subroutine classify_3objects
       end subroutine calc_percentage_phiber

       subroutine generate_gaussian_refs(refs,ldim,smpd,min_rad,max_rad)
           type(image), allocatable, intent(inout) :: refs(:)
           integer,     intent(in) :: ldim(3)
           real,        intent(in) :: smpd
           real,        intent(in) :: min_rad, max_rad
           integer, parameter :: STEP = 5
           integer :: nrefs
           real    :: sigma_x, sigma_y
           integer :: i, j
           logical :: rotate_ref
           real, allocatable :: rmat_out(:,:,:) !to use rtsq_serial
           nrefs = 0
           do i = 0, STEP-1
               sigma_x = min_rad + i*STEP
               if(sigma_x < max_rad) then
                   do j = 0, STEP-1
                       sigma_y = min_rad + j*STEP
                       if(sigma_y < max_rad) then
                           nrefs = nrefs + 1
                           if(sigma_y/sigma_x > 1.5) rotate_ref = .true.
                           if(rotate_ref) nrefs = nrefs + 2 !total 3 references
                       endif
                       rotate_ref = .false. !restore
                   enddo
               endif
           enddo
           if(allocated(refs)) deallocate(refs)
           allocate( refs(nrefs) )
           allocate(rmat_out(ldim(1),ldim(2),1), source = 0.)
           nrefs = 0
           do i = 0, STEP-1
               sigma_x = min_rad + i*STEP
               if(sigma_x < max_rad) then
                   do j = 0, STEP-1
                       sigma_y = min_rad + j*STEP
                       if(sigma_y < max_rad) then
                            print *, 'generating ref with sigmax = ', sigma_x, 'sigmay = ', sigma_y
                            if(sigma_y/sigma_x > 1.5) rotate_ref = .true.
                            nrefs = nrefs + 1
                            call refs(nrefs)%new(ldim,smpd)
                            call refs(nrefs)%gauimg2D(sigma_x,sigma_y)
                            call refs(nrefs)%write('_GaussianReference.mrc',nrefs)
                            if(rotate_ref) then
                                call refs(nrefs)%rtsq_serial( 45., 0., 0., rmat_out )
                                nrefs = nrefs + 1
                                call refs(nrefs)%new(ldim,smpd)
                                call refs(nrefs)%set_rmat(rmat_out)
                                call refs(nrefs)%write('_GaussianReference.mrc',nrefs)
                                call refs(nrefs)%rtsq_serial( 90., 0., 0., rmat_out )
                                nrefs = nrefs + 1
                                call refs(nrefs)%new(ldim,smpd)
                                call refs(nrefs)%set_rmat(rmat_out)
                                call refs(nrefs)%write('_GaussianReference.mrc',nrefs)
                            endif
                       endif
                       rotate_ref = .false. !restore
                   enddo
               endif
           enddo
           deallocate(rmat_out)
           ! sigma_x = (self%min_rad)/2. ! half of the radius
           ! sigma_y = (self%max_rad)/2.
           ! rotate_ref = .false.
           ! if(sigma_y/sigma_x > 1.5) rotate_ref = .true.
           ! print *, 'sigma_x = ', sigma_x, 'sigma_y = ', sigma_y
           ! rot_step = 180./real(N_ROT) !it's gaussian, so rotate just 180 degrees
           ! call self%reference%gauimg2D(sigma_x,sigma_y)
           ! call self%reference%write(PATH_HERE//basename(trim(self%fbody))//'_GaussianReference.mrc')
           ! call field%fft()
           ! allocate(rmat_out(self%ldim_shrunken(1),self%ldim_shrunken(2),1), source = 0.)
           ! call self%phasecorr%zero_and_unflag_ft()
           ! print *, 'rotate_ref = ', rotate_ref
       end subroutine generate_gaussian_refs
       end module simple_test_chiara_try_mod

    program simple_test_chiara_try
       include 'simple_lib.f08'
       use simple_math
       use simple_jpg
       use simple_test_chiara_try_mod
       use simple_micops
       use simple_image, only : image
       use simple_tvfilter, only : tvfilter
       use simple_nanoparticles_mod, only : nanoparticle
       type(image), allocatable :: refs(:)
       integer :: i,j, cnt, ldim(3)
       character(len=100) :: fname
       real :: smpd, min_rad, max_rad
       smpd    = 1.
       min_rad = 10.
       max_rad = 40.
       ldim(1) = 256
       ldim(2) = 256
       ldim(3) = 1
       call generate_gaussian_refs(refs,ldim,smpd,min_rad,max_rad)

       !fname = 'ImgTest.jpg'
       !call calc_percentage_phiber(fname)

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
    !ARTIFICIAL DATA TEST FOR PHASECORRELATIONS
    ! sigma_x = 5.
    ! sigma_y = 10.
    ! call field%new([512,512,1],1.)
    ! call phasecorr%new([512,512,1],1.)
    ! call reference%new([512,512,1],1.)
    ! call aux%new([512,512,1],1.)
    !
    ! allocate(rmat(512,512,1), source=0.)
    ! rmat(256-50:256+50, 256-10:256+10,1) = 1.
    ! call field%set_rmat(rmat)
    ! call field%add_gauran(0.1)
    ! call field%write('RectangleNoise.mrc')
    !
    !
    ! call reference%gauimg2D(sigma_x,sigma_y)
    ! call reference%add_gauran(0.2)
    ! step = 180./real(N_ROT)
    ! call phasecorr%zero_and_unflag_ft()
    ! call field%fft()
    ! allocate(rmat_out(512,512,1), source = 0.)
    ! do n_ang = 1, N_ROT
    !     if(n_ang > 1) then !do not rotate the first time
    !         if(reference%is_ft()) call reference%ifft() ! in order to rotate it has to be real
    !         if(phasecorr%is_ft()) call phasecorr%ifft() ! in order to rotate it has to be real
    !         call reference%rtsq_serial( step, 0., 0., rmat_out )
    !         call reference%set_rmat(rmat_out)
    !         call reference%write(trim(int2str(n_ang))//'ReferenceNoise.mrc')
    !     else
    !         call reference%write('1ReferenceNoise.mrc')
    !     endif
    !     call reference%fft()
    !     call phasecorr%fft()    ! before sum they need to be in the same state
    !     aux = reference%conjg() ! necessary
    !     phasecorr = phasecorr + field*aux ! sum of them
    !     call phasecorr%ifft()
    !     ! FORMULA: phasecorr = ifft2(fft2(field).*conj(fft2(reference)));
    ! enddo
    ! call phasecorr%div(real(N_ROT-1))
    ! call phasecorr%write('SumPhaseCorrNoise.mrc')

    !REAL DATA TEST FOR PHASE CORRELATIONS
    ! ldim_shrunken(1) = 1854
    ! ldim_shrunken(2) = 1918
    ! ldim_shrunken(3) = 1
    ! smpd_shrunken = 0.66*SHRINK
    ! min_rad = 70.
    ! max_rad = 75.
    ! detector = 'bin'
    ! call img%new   (ldim_shrunken,smpd_shrunken)
    ! call img_cc%new(ldim_shrunken, smpd_shrunken)
    ! lambda = 5.
    ! hp_box =  4.*max_rad+2.*max_rad
    ! ! 0) Reading and saving original micrograph
    ! call read_micrograph('14sep05c_c_00003gr_00014sq_00002hl_00005es_c.mrc', smpd=0.66)
    ! ! 1) Shrink and high pass filtering
    ! call shrink_micrograph(SHRINK,ldim_shrunken ,smpd_shrunken)
    ! call set_box(int(SHRINK*(hp_box)), box_shrunken, micrograph_shrunken)
    ! call img%copy(micrograph_shrunken)
    ! ! To take care of shrinking
    ! min_rad = min_rad/SHRINK ! I am thingking I shouldn't multiply by the smpd cuz I am working in pxls
    ! max_rad = max_rad/SHRINK
    ! ! 2) Low pass filtering
    ! call micrograph_shrunken%bp(0.,20.)
    ! call micrograph_shrunken%ifft()
    ! !2.1) TV denoising
    ! call tvf%new()
    ! call tvf%apply_filter(micrograph_shrunken, lambda)
    ! call tvf%kill()
    ! call micrograph_shrunken%write(PATH_HERE//'_tvfiltered.mrc')
    ! ! 2.2) negative image, to obtain a binarization with white particles
    ! call micrograph_shrunken%neg() !TO REMOVE IN CASE OF NEGATIVE STAINING
    ! ! 3) Reference generation and Phase Correlation calculation
    ! sigma_x = 3.5
    ! sigma_y = 3.5
    ! call reference%new(ldim_shrunken,smpd_shrunken)
    ! call phasecorr%new(ldim_shrunken,smpd_shrunken)
    ! call reference%gauimg2D(sigma_x,sigma_y)
    ! call reference%add_gauran(0.2)
    ! call reference%write('_ReferenceGau.mrc')
    ! call micrograph_shrunken%fft()
    ! call reference%fft()
    ! aux = reference%conjg() ! Check if necessary
    ! phasecorr = micrograph_shrunken*aux
    ! call phasecorr%ifft()
    ! call phasecorr%write('_PhaseCorrelationsSNR02.mrc')
    ! ! 3) Binarization
    ! call phasecorr%stats( ave, sdev, maxv, minv )
    ! if(detector .eq. 'sobel') then
    !     thresh(1) = ave+.5*sdev !sobel needs lower thresh not to pick just edges
    !     call sobel(phasecorr,thresh)
    ! else if (detector .eq. 'bin') then
    !     call phasecorr%bin(ave+.8*sdev)
    ! else if (detector .eq. 'otsu') then
    !     rmat = phasecorr%get_rmat()
    !     x = pack(rmat, .true.)
    !     call otsu(x,x_out)
    !     rmat = reshape(x_out, [ldim_shrunken(1),ldim_shrunken(2),1])
    !     call phasecorr%set_rmat(rmat)
    !     call phasecorr%erosion() !morphological erosion
    !     deallocate(x,x_out,rmat)
    ! else
    !     print *, 'Invalid detector; preprocess_mic'
    ! endif
    ! call phasecorr%write(PATH_HERE//'_Bin.mrc')
    ! winsz = int(min_rad+max_rad)/4 !half of the avg between the dimensions of the particles
    ! call phasecorr%real_space_filter(winsz,'median') !median filtering allows easy calculation of cc
    ! call phasecorr%write(PATH_HERE//'_BinMedian.mrc')



    ! phasecorr = ifft2(fft2(field).*conj(fft2(reference)));
    ! figure;
    ! plot_im(phasecorr);
    ! end
    !PICKER TESTS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ldim_shrunken(1) = 1854
    ! ldim_shrunken(2) = 1918
    ! ldim_shrunken(3) = 1
    ! smpd_shrunken = 0.66*SHRINK
    ! min_rad = 70.
    ! max_rad = 75.
    ! detector = 'bin'
    ! call img%new   (ldim_shrunken,smpd_shrunken)
    ! call img_cc%new(ldim_shrunken, smpd_shrunken)
    ! lambda = 5.
    ! hp_box =  4.*max_rad+2.*max_rad
    ! ! 0) Reading and saving original micrograph
    ! call read_micrograph('14sep05c_c_00003gr_00014sq_00002hl_00005es_c.mrc', smpd=0.66)
    ! ! 1) Shrink and high pass filtering
    ! call shrink_micrograph(SHRINK,ldim_shrunken ,smpd_shrunken)
    ! call set_box(int(SHRINK*(hp_box)), box_shrunken, micrograph_shrunken)
    ! call img%copy(micrograph_shrunken)
    ! ! To take care of shrinking
    ! min_rad = min_rad/SHRINK ! I am thingking I shouldn't multiply by the smpd cuz I am working in pxls
    ! max_rad = max_rad/SHRINK
    ! ! 2) Low pass filtering
    ! call micrograph_shrunken%bp(0., 20.)
    ! call micrograph_shrunken%ifft()
    ! call micrograph_shrunken%write('_mic_shrunken.mrc')
    ! ! NEW APPROACH
    !
    ! call micrograph_shrunken%hist_stretching(micrograph_shrunken)
    ! call micrograph_shrunken%write(PATH_HERE//'_hist_stretch.mrc')
    ! call micrograph_shrunken%scale_pixels([1.,64.])
    ! ! rmat = micrograph_shrunken%get_rmat()
    ! ! rmat = sqrt(rmat)
    ! ! call micrograph_shrunken%set_rmat(rmat)
    ! ! call micrograph_shrunken%write(PATH_HERE//'_sqrt.mrc')
    !
    ! rmat = micrograph_shrunken%get_rmat()
    ! rmat = log(rmat)
    ! call micrograph_shrunken%set_rmat(rmat)
    ! call micrograph_shrunken%write(PATH_HERE//'log.mrc')
    !
    ! ! 2.1) TV denoising
    ! call tvf%new()
    ! call tvf%apply_filter(micrograph_shrunken, lambda)
    ! call tvf%kill()
    ! call micrograph_shrunken%write(PATH_HERE//'_log_tvfiltered.mrc')
    ! ! 2.2) negative image, to obtain a binarization with white particles
    ! call micrograph_shrunken%neg() !TO REMOVE IN CASE OF NEGATIVE STAINING
    ! ! 3) Binarization
    ! call micrograph_shrunken%stats( ave, sdev, maxv, minv )
    ! if(detector .eq. 'sobel') then
    !     thresh(1) = ave+.5*sdev !sobel needs lower thresh not to pick just edges
    !     call sobel(micrograph_shrunken,thresh)
    ! else if (detector .eq. 'bin') then
    !     call micrograph_shrunken%bin(ave+.8*sdev)
    ! else if (detector .eq. 'otsu') then
    !     rmat = micrograph_shrunken%get_rmat()
    !     x = pack(rmat, .true.)
    !     call otsu(x,x_out)
    !     rmat = reshape(x_out, [ldim_shrunken(1),ldim_shrunken(2),1])
    !     call micrograph_shrunken%set_rmat(rmat)
    !     call micrograph_shrunken%erosion() !morphological erosion
    !     deallocate(x,x_out,rmat)
    ! else
    !     print *, 'Invalid detector; preprocess_mic'
    ! endif
    ! call micrograph_shrunken%write(PATH_HERE//'_Bin.mrc')
    ! winsz = int(min_rad+max_rad)/4 !half of the avg between the dimensions of the particles
    ! call micrograph_shrunken%real_space_filter(winsz,'median') !median filtering allows easy calculation of cc
    ! call micrograph_shrunken%write(PATH_HERE//'_BinMedian.mrc')
    ! ! 5) Connected components (cc) identification
    ! ! call micrograph_shrunken%find_connected_comps(img_cc)
    ! ! 6) cc filtering
    ! ! print *, 'before polishing the ccs: ', self%get_n_ccs()
    ! ! call self%img_cc%polish_cc(self%min_rad,self%max_rad)
    ! ! print *, 'after polishing the ccs: ', self%get_n_ccs()
    ! ! if( DEBUG_HERE ) call self%img_cc%write_jpg(PATH_HERE//basename(trim(self%fbody))//'_ConnectedComponentsElimin.jpg')
    ! ! call self%img_cc%write(PATH_HERE//basename(trim(self%fbody))//'_ConnectedComponentsElimin.mrc')
    ! ! call micrograph_shrunken%kill
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! CURVATURE ESTIMATION ON REAL DATA
    ! call img%new([512,512,1],1.34)
    ! print *, 'Analysing 008_movieb_forctf_binarized_polished.mrc GOOD ONE'
    ! call img%read('008_movieb_forctf_binarized_polished.mrc')
    ! call img%find_connected_comps(img_cc)
    ! sz = img_cc%size_connected_comps()
    ! avg_curvature = 0.
    ! do i = 1, number_biggest_ccs
    !     cc(:) = maxloc(sz)
    !     c = estimate_curvature(img_cc,cc(1),i)
    !     avg_curvature = avg_curvature + c
    !     call img_aux%copy(img_cc)
    !     imat = nint(img_aux%get_rmat())
    !     where(imat /= cc(1)) imat = 0
    !     call img_aux%set_rmat(real(imat))
    !     call img_aux%write(trim(int2str(i))//'BiggestCC.mrc')
    !     sz(cc(1)) = 0 !discard
    ! enddo
    ! avg_curvature = avg_curvature/number_biggest_ccs
    ! print *, 'avg_curvature = ', avg_curvature
    !
    ! print *, 'Analysing May08_03.05.02.bin_forctf_binarized_polished.mrc APOFERRITIN'
    ! call img%read('May08_03.05.02.bin_forctf_binarized_polished.mrc')
    ! call img%find_connected_comps(img_cc)
    ! sz = img_cc%size_connected_comps()
    ! avg_curvature = 0.
    ! do i = 1, number_biggest_ccs
    !     cc(:) = maxloc(sz)
    !     c = estimate_curvature(img_cc,cc(1),i)
    !     avg_curvature = avg_curvature + c
    !     call img_aux%copy(img_cc)
    !     imat = nint(img_aux%get_rmat())
    !     where(imat /= cc(1)) imat = 0
    !     call img_aux%set_rmat(real(imat))
    !     call img_aux%write(trim(int2str(i))//'BiggestCC.mrc')
    !     sz(cc(1)) = 0 !discard
    ! enddo
    ! avg_curvature = avg_curvature/number_biggest_ccs
    ! print *, 'avg_curvature = ', avg_curvature
    !
    ! print *, 'Analysing 0315_intg_binarized_polished.mrc FALLACIOUS'
    ! call img%read('0315_intg_binarized_polished.mrc')
    ! call img%find_connected_comps(img_cc)
    ! sz = img_cc%size_connected_comps()
    ! avg_curvature = 0.
    ! do i = 1, number_biggest_ccs
    !     cc(:) = maxloc(sz)
    !     c = estimate_curvature(img_cc,cc(1),i)
    !     avg_curvature = avg_curvature + c
    !     call img_aux%copy(img_cc)
    !     imat = nint(img_aux%get_rmat())
    !     where(imat /= cc(1)) imat = 0
    !     call img_aux%set_rmat(real(imat))
    !     call img_aux%write(trim(int2str(i))//'BiggestCC.mrc')
    !     sz(cc(1)) = 0 !discard
    ! enddo
    ! avg_curvature = avg_curvature/number_biggest_ccs
    ! print *, 'avg_curvature = ', avg_curvature

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


        ! call img%new([512,512,1], 1.41)
    ! !call img%read('pspecs_saga_polii.mrc', 71)
    ! rmat = img%get_rmat()
    ! rmat = 5.
    ! call img%set_rmat(rmat)
    ! call img%write('RmatUniform.mrc')
    ! X = reshape(rmat, [512*512])

    ! Entropy calculation
    ! X  = [40.,50.,60.,17.,17.,17.]
    ! print *, 'entropy = ', entropy_shells(X,67.,4.)
    ! print *, ' before:', entropy(X,64)
    !
    !
    ! X  = [1./6.,1./6.,1./6.,1./6.,1./6.,1./6.]
    ! print *, 'entropy = ', entropy_shells(X,1.,0.)
    ! print *, ' before:', entropy(X,64)
    !
    ! X  = [0.,0.,60.,0.,0.,0.]
    ! print *, 'entropy = ', entropy_shells(X,67.,0.)
    ! print *, ' before:', entropy(X,64)
    !
    ! X  = [50.,40.,60.,30.,20.,10.]
    ! print *, 'entropy = ', entropy_shells(X,67.,0.)
    ! print *, ' before:', entropy(X,64)


    !print *, 'min X = ', minval(X), 'max(X) = ', maxval(X), 'shape(X)', shape(X)
    !e = entropy_try(X,6)
    !print *, 'e = ', e ! USE SCALE IMAGE???
    ! ! centers1 = reshape([1.,1.,1.,1.5,1.5,1.5,2.3,2.4,2.5,4.1,4.3,4.7],[3,4])
    ! print *, 'centers1(:3,1) = ',centers1(:3,1)
    ! call vis_mat(centers1)
    ! centers2 = reshape([2.1,2.3,2.5,1.,0.7,0.6,1.4,1.3,1.6,4.3,3.9,4.9],[3,4])
    ! print *, 'centers2 = '
    ! call vis_mat(centers2)
    ! call calc_rmsd(centers1,centers2)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! vec = reshape([1.,2.,0.,0.,5.,0.,0.,8.,9.,0., 0.,12.,13.,0.,0.,16.,17.,0.,0.,20.], [10,2])
    ! print *, 'Vec Before = '
   ! call vis_mat(vec)
    ! packed_vec = pack(vec,vec  > TINY)
    ! print *, 'Vec After = '
   ! call vis_mat(packed_vec)

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
