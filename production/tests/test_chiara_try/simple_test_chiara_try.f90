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
  ! subroutine aff_prop
  !     real                 :: centers(3,303)
  !     type(aff_prop)       :: apcls
  !     real                 :: simmat(900,900), simsum
  !     integer, allocatable :: centers(:), labels(:)
  !     integer              :: i, j, ncls, nerr
  !     write(logfhandle,'(a)') '**info(simple_aff_prop_unit_test): testing all functionality'
  !     ! make data
  !     do i=1,300
  !         datavecs(i,:) = 1.
  !     end do
  !     do i=301,600
  !         datavecs(i,:) = 5.
  !     end do
  !     do i=601,900
  !         datavecs(i,:) = 10.
  !     end do
  !     do i=1,900-1
  !         do j=i+1,900
  !             simmat(i,j) = -euclid(datavecs(i,:),datavecs(j,:))
  !             simmat(j,i) = simmat(i,j)
  !         end do
 !     end do
  !     call apcls%new(900, simmat)
  !     call apcls%propagate(centers, labels, simsum)
  !     ncls = size(centers)
  !     nerr = 0
  !     do i=1,299
  !         do j=i+1,300
  !             if( labels(i) /= labels(j) ) nerr = nerr+1
  !         end do
  !     end do
  !     do i=301,599
  !         do j=i+1,600
  !             if( labels(i) /= labels(j) ) nerr = nerr+1
  !         end do
  !     end do
  !     do i=601,899
  !         do j=i+1,900
  !             if( labels(i) /= labels(j) ) nerr = nerr+1
  !         end do
  !     end do
  !     write(logfhandle,*) 'NR OF CLUSTERS FOUND:', ncls
 !     write(logfhandle,*) 'NR OF ASSIGNMENT ERRORS:', nerr
  !     write(logfhandle,*) 'CENTERS'
  !     do i=1,size(centers)
  !         write(logfhandle,*) datavecs(centers(i),:)
  !     end do
  !     if( ncls == 3 .and. nerr == 0 )then
  !         write(logfhandle,'(a)') 'SIMPLE_AFF_PROP_UNIT_TEST COMPLETED ;-)'
  !     else
  !         write(logfhandle,'(a)') 'SIMPLE_AFF_PROP_UNIT_TEST FAILED!'
  !     endif
  ! end subroutine aff_prop

  ! This subroutine takes in input coords of the atomic positions
  ! in the nanoparticles to compare and calculates the rmsd
  ! between them.
  ! See formula
  ! https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions
  subroutine calc_rmsd(centers1,centers2,r)
      real, intent(in) :: centers1(:,:), centers2(:,:)
      real, optional, intent(out) :: r !rmsd calculated
      integer :: N, i, j
      real    :: sum, rmsd
      real, allocatable :: dist(:)
      logical, allocatable :: mask(:)
      integer :: location(1)
      ! If they don't have the same nb of atoms.
      if(size(centers1, dim = 2) <= size(centers2, dim = 2)) then
          N = size(centers1, dim = 2) !compare based on centers1
          print *, 'N = ', N
          allocate(dist(N), source = 0.)
          allocate(mask(size(centers2, dim = 2)), source = .true.)
          do i = 1, N
              dist(i)  = pixels_dist(centers1(:,i),centers2(:,:),'min',mask,location)
              print *, 'pixel = ', i, 'location = ', location
              dist(i) = dist(i)**2 !formula wants them square
              mask(location(1)) = .false. ! not to consider the same atom more than once
          enddo
          print *, 'dist = ', dist, 'sum_dist = ', sum(dist), 'sum(dist)/N = ', sum(dist)/real(N)
          rmsd = sqrt(sum(dist)/real(N))
      print *, 'RMSD = ', rmsd
      else
          N = size(centers2, dim = 2) !compare based on centers2
          allocate(dist(N), source = 0.)
          allocate(mask(size(centers1, dim = 2)), source = .true.)
          do i = 1, N
              dist(i) = pixels_dist(centers2(:,i),centers1(:,:),'min',mask,location)
              dist(i) = dist(i)**2        !formula wants them square
              mask(location(1)) = .false. ! not to consider the same atom more than once
      enddo
          rmsd = sqrt(sum(dist)/real(N))
          print *, 'RMSD = ', rmsd

      endif
      if(present(r)) r=rmsd
  end subroutine calc_rmsd

  subroutine calc_circ(img, i)
      type(image), intent(inout) :: img
      integer, optional, intent(in) :: i
      type(image) :: img_cc
      real,    allocatable :: rmat(:,:,:)
      logical, allocatable :: border(:,:,:)
      integer, allocatable :: imat(:,:,:)
  integer :: label
      integer :: v(1), s(1) !volumes and surfaces of each atom
      real, parameter :: pi = 3.14159265359
      real :: circularity
      integer, allocatable :: pos(:,:)
      logical, allocatable :: mask_dist(:) !for min and max dist calculation
      integer :: location(1) !location of the farest vxls of the atom from its center
      integer :: loc_min(1)
      real    :: longest_dist
      real :: volume, surface
      !call img%find_connected_comps(img_cc)
      ! call img_cc%write('ImgCc.mrc')
      ! rmat = img_cc%get_rmat()
      ! allocate(imat(50,50,50),source = nint(rmat)) !for function get_pixel_pos
      ! call img_cc%border_mask(border, 1, .true.) !use 4neigh instead of 8neigh
      ! where(border .eqv. .true.)
      !     imat = 1
      ! elsewhere
      !     imat = 0
      ! endwhere
      ! call get_pixel_pos(imat,pos)   !pxls positions of the shell
      ! if(allocated(mask_dist)) deallocate(mask_dist)
      ! allocate(mask_dist(size(pos, dim = 2)), source = .true. )
      ! longest_dist  = pixels_dist([24.5,24.5,24.5], real(pos),'max', mask_dist, location)
      ! print *, 'radius = ', longest_dist
      ! volume  = (4.*pi)/3.*(longest_dist)**3
      ! print *, 'volume =', volume
      ! surface = (4.*pi)*(longest_dist)**2
      ! circularity = (6.*sqrt(pi)*volume)/ sqrt((surface**3))
      ! print *, 'circularity = ',circularity
      ! stop
      !rmat = img_cc%get_rmat()
      rmat = img%get_rmat()
      allocate(imat(160,160,160),source = nint(rmat)) !for function get_pixel_pos
      v(1) = count(imat == 1) !just count the number of vxls in the cc
      call img%border_mask(border, 1)!, .true.)
      rmat = 0.
      where(border .eqv. .true.) rmat = 1.
      call img%set_rmat(rmat)
      if(present(i)) then
        call img%write(int2str(i)//'Border.mrc')
      else
          call img%write('Border.mrc')
      endif
      s(1) = count(border .eqv. .true.)
      circularity = (6.*sqrt(pi)*real(v(1)))/sqrt(real(s(1))**3)
      print*, 'vol = ', v(1), 'surf = ', s(1), 'circ = ', circularity
      deallocate(imat, border, rmat)
  end subroutine calc_circ


  ! function entropy(X,n) result(e)
  !     use simple_math
  !     real, intent(in)    :: X(:)
  !     integer, intent(in) :: N
  !     real                :: e !entropy value
  !     real,    allocatable :: xhist(:) !discretization of the values
  !     integer, allocatable :: yhist(:) !number of occurences
  !     real,    allocatable :: p(:), p_no_zero(:)
  !     integer :: i, cnt
  !     call create_hist_vector(X,n,xhist,yhist)
  !     allocate(p(size(yhist)), source = 0.)
  !     p =  real(yhist)/real(sum(yhist)) !probabilities
  !     print *, 'p = ', p
  !     cnt = count(p>TINY)
  !     print *, 'cnt = ', cnt
  !     allocate(p_no_zero(cnt), source = 0.)
  !     cnt = 0 !reset
  !     do i = 1, size(p)
  !         if(p(i) > TINY) then
  !             cnt = cnt + 1
  !             print *, 'i = ', i, 'p(i) = ', p(i)
  !             p_no_zero(cnt) = p(i)
  !         endif
  !     enddo
  !     e = -sum(p_no_zero*(log10(p_no_zero)/log10(2.))) !formula: sum(p*log2(p))
  !     deallocate(xhist,yhist)
  ! end function entropy

  ! function entropy_try(X,n) result(e)
  !     real, intent(in)    :: X(:)
  !     integer, intent(in) :: N
  !     real                :: e !entropy value
  !     real,    allocatable :: xhist(:) !discretization of the values
  !     integer, allocatable :: yhist(:) !number of occurences
  !     real,    allocatable :: p(:), p_no_zero(:)
  !     integer :: i, cnt
  !     call create_hist_vector(X,n,xhist,yhist)
  !     allocate(p(size(yhist)), source = 0.)
  !     p =  real(yhist)/real(sum(yhist)) !probabilities
  !     cnt = count(p>TINY)
  !     allocate(p_no_zero(cnt), source = 0.)
  !     cnt = 0 !reset
  !     do i = 1, size(p)
  !         if(p(i) > TINY) then
  !             cnt = cnt + 1
  !             p_no_zero(cnt) = p(i) !remove zeroes occurrencies
  !         endif
  !     enddo
  !     e = -sum(p_no_zero*(log10(p_no_zero)/log10(2.))) !formula: sum(p*log2(p))
  !     deallocate(xhist,yhist)
  ! end function entropy_try


  !From the paper 'Approximate Entropy for Testing Randomness'
  ! I am writing it to compare binary strings (vectors)
  ! function approx_entropy(X,m) result(e)
  !     integer,    intent(inout) :: X(:)
  !     integer, intent(in)    :: m
  !     real :: e ! approximate entropy
  !     integer, allocatable :: C(:)
  !     real, allocatable :: Phi(:)
  !     integer :: n, h, hh
  !     integer, allocatable :: cnt(:)
  !     n = size(X, dim = 1) ! number of elements in X
  !     allocate(C  (n-m+1), source = 0 )
  !     allocate(Phi(n-m+1), source = 0.)
  !     allocate(cnt(n-m+1), source = 0 )
  !     do h = 1, n-m+1
  !         do hh = 1, n-m+1
  !             if(vectors_are_equal(X(hh:hh+m-1), X(h:h+m-1)) .and. h /= hh) cnt(h) = cnt(h) + 1  !I don't know if h/=hh is necessary
  !         enddo
  !         C(h) = cnt(h)/(n-m+1)
  !     enddo
  !     print *, 'cnt = ', cnt
  !     print *, 'C = ', C
  !     do h = 1, n-m+1
  !         Phi(h) = Phi(h)+log(real(C(h)))
  !     enddo
  !     print *, 'Phi = ', Phi
  !     e = Phi(m)-Phi(m+1)
  ! end function approx_entropy


  ! This subroutine takes in input 2 2D vectors, centered in the origin
  ! and it gives as an output the angle between them, IN DEGREES.
  function ang2D_vecs(vec1, vec2) result(ang)
      real, intent(inout) :: vec1(2), vec2(2)
      real :: ang        !output angle
      real :: ang_rad    !angle in radians
      real :: mod1, mod2
      real :: dot_prod
      mod1 = sqrt(vec1(1)**2+vec1(2)**2)
      mod2 = sqrt(vec2(1)**2+vec2(2)**2)
      ! normalise
      vec1 = vec1/mod1
      vec2 = vec2/mod2
      ! dot product
      dot_prod = vec1(1)*vec2(1)+vec1(2)*vec2(2)
      ! sanity check
      if(dot_prod > 1. .or. dot_prod< -1.) then
         ! THROW_WARN('Out of the domain of definition of arccos; ang2D_vecs')
          ang_rad = 0.
      else
          ang_rad = acos(dot_prod)
      endif
      ! if(DEBUG_HERE) then
      !     write(logfhandle,*)'>>>>>>>>>>>>>>>>>>>>>>'
      !     write(logfhandle,*)'mod_1     = ', mod1
      !     write(logfhandle,*)'mod_2     = ', mod2
      !     write(logfhandle,*)'dot_prod  = ', dot_prod
      !     write(logfhandle,*)'ang in radians', acos(dot_prod)
      ! endif
      !output angle
      ang = rad2deg(ang_rad)
  end function ang2D_vecs

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

    function phase_corr_try(self1,self2,border) result(pc)
        type(image),      intent(inout) :: self1, self2
        integer, optional, intent(in)    :: border
        type(image) :: pc
        type(image) :: self1_crop, self2_crop, pc_crop
        integer     :: bborder
        integer     :: i, j
        integer :: ldim(3)
        real :: smpd
        real, pointer :: rmat1(:,:,:), rmat2(:,:,:)
        real, pointer :: rmat1_crop(:,:,:), rmat2_crop(:,:,:)
        ldim = self1%get_ldim()
        smpd = self1%get_smpd()
        bborder = 0
        if(present(border)) bborder = border
        call self1_crop%new([ldim(1)-2*bborder,ldim(2)-2*bborder,1], smpd)
        call self2_crop%new([ldim(1)-2*bborder,ldim(2)-2*bborder,1], smpd)
        call pc_crop%new   ([ldim(1)-2*bborder,ldim(2)-2*bborder,1], smpd)
        call pc%new(ldim, smpd)
        call self1%get_rmat_ptr(rmat1)
        call self2%get_rmat_ptr(rmat2)
        call self1_crop%get_rmat_ptr(rmat1_crop)
        call self2_crop%get_rmat_ptr(rmat2_crop)
        do i = bborder+1, ldim(1)-bborder
            do j = bborder+1, ldim(2)-bborder
                 rmat1_crop(i,j,1) =  rmat1(i,j,1) ! maybe instead of using the functions every time it's convinient to use matrices
                 rmat2_crop(i,j,1) =  rmat2(i,j,1) ! maybe instead of using the functions every time it's convinient to use matrices
            enddo
        enddo
        call self1_crop%fft()
        call self2_crop%fft()
        self2_crop = self2_crop%conjg()
        pc_crop = self1_crop*self2_crop
        call pc_crop%ifft()
        call pc_crop%pad(pc,0.)
    end function phase_corr_try
end module simple_test_chiara_try_mod

program simple_test_chiara_try
    include 'simple_lib.f08'
    use simple_math
    use simple_test_chiara_try_mod
    use simple_micops
    use simple_image, only : image
    use simple_tvfilter, only : tvfilter
    use simple_nanoparticles_mod, only : nanoparticle
    type(image) :: img1, img2, pc
    real, allocatable :: rmat(:,:,:)
    !BORDER EFFECTS IN PHASECORRELATION EXPLAINATION
    call img1%new([512,512,1],1.)
    call img2%new([512,512,1],1.)
    allocate(rmat(512,512,1), source=0.)
    rmat(256-50:256+50, 256-10:256+10,1) = 1.
    rmat(10:50,50:80,1) = 1.
    call img1%set_rmat(rmat)
    call img1%write('RectangleIMG.mrc')
    call img2%gauimg2D(50.,10.)
    call img2%write('GaussianIMG.mrc')
    pc = img1%phase_corr(img2,1)
    call pc%write('PhaseCorrelationIMG.mrc')
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
