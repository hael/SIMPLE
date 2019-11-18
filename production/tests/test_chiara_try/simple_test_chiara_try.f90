module simple_test_chiara_try_mod
    include 'simple_lib.f08'
    use simple_image, only : image

    public

    type emgaufit
        private
        integer           :: N    ! number of data
        integer           :: K    ! number of components
        real, allocatable :: avgs(:), vars(:), mix(:), gammas(:,:)
        real, pointer     :: dat(:)=>null()
        logical           :: exists=.false.
      contains
        ! CONSTRUCTOR
        procedure :: new
        ! GETTERS/SETTERS
        procedure :: set_data
        procedure :: set_avgs
        procedure :: get_avgs
        procedure :: get_vars
        ! MACHINERY
        procedure :: fit
        procedure, private :: init
        procedure, private :: init_chiara
        procedure, private :: logL
        procedure, private :: estep
        procedure, private :: mstep
        ! DESTRUCTOR
        procedure :: kill
    end type

    contains

        ! CONSTRUCTOR

        !>  \brief  is a constructor
        subroutine new( self, N, K )
            class(emgaufit), intent(inout) :: self
            integer, intent(in)            :: N, K
            integer :: alloc_stat
            if( self%exists ) call self%kill
            self%N = N
            self%K = K
            allocate( self%avgs(K), self%vars(K), self%mix(K), self%gammas(N,K), stat=alloc_stat )
        end subroutine

        ! GETTERS/SETTERS

        !>  \brief  is for setting pointer to data
        subroutine set_data( self, dat )
            class(emgaufit), intent(inout) :: self
            real, intent(in), target       :: dat(:)
            self%dat => dat
        end subroutine

        !>  \brief  is for setting the averages
        subroutine set_avgs( self, avgs)
            class(emgaufit), intent(inout) :: self
            real,            intent(in)    :: avgs(self%K)
            self%avgs = avgs
        end subroutine

        !>  \brief  is for getting the averages
        function get_avgs( self ) result( avgs )
            class(emgaufit), intent(in) :: self
            real :: avgs(self%K)
            avgs = self%avgs
        end function

        !>  \brief  is for getting the variances
        function get_vars( self ) result( vars )
            class(emgaufit), intent(in) :: self
            real :: vars(self%K)
            vars = self%vars
        end function

        ! MACHINERY

        !>  \brief  is the master fitting subroutine
        subroutine fit( self, mits, avgs)
            class(emgaufit), intent(inout) :: self
            integer,         intent(in)    :: mits
            real   ,         intent(in)    :: avgs(self%K)
            integer :: k
            real    :: L
            ! initialize
            ! call self%init
            call self%init_chiara(avgs)
            ! calculate initial log-likelihood
            L = self%logL()
            print *, 'initial log(L):', L
            ! iterate
            do k=1,mits
                call self%estep
                call self%mstep
                L = self%logL()
                if( k == 1 .or. mod(k,5) == 0 )then
                    write(*,"(1X,A,1X,I3,1X,A,1X,F10.0)") 'Iteration:', k, 'Log-likelihood:', L
                endif
            end do
        end subroutine

        !>  \brief  is for random initialization of the gammas
        subroutine init( self )
            use simple_rnd, only: ran3
            class(emgaufit), intent(inout) :: self
            integer :: i, j
            real    :: P(self%K), S
            do i=1,self%N
                S = 0.
                do j=1,self%K
                    P(j) = ran3()
                    S = S+P(j)
                end do
                do j=1,self%K
                    self%gammas(i,j) = P(j)/S
                end do
            end do
            call self%mstep
            print *, 'self%avgs        ',self%avgs
            print *, 'self%mix         ',self%mix
            print *, 'self%vars        ', self%vars
            ! print *, 'self%gammas      ', self%gammas
        end subroutine

        !>  \brief  is for random initialization of the gammas
        subroutine init_chiara( self, avgs )
            use simple_rnd, only: ran3
            class(emgaufit), intent(inout) :: self
            real,            intent(in)    :: avgs(self%K) ! input avgs are kmeans output
            integer :: i, j
            real    :: P(self%K), S
            ! set averages
            call self%set_avgs(avgs)
            print *, 'self%avgs   CHIARA',self%avgs
            ! set mixing coeffs
            self%mix = 1./real(self%K)
            print *, 'self%mix   CHIARA',self%mix
            ! set variances
            do i = 1, self%K
                self%vars(i) = 10.*ran3()
            enddo
            print *, 'self%vars  CHIARA', self%vars
            ! set gammas
            call self%estep
            ! print *, 'self%gammas CHIARA', self%gammas
        end subroutine

        !>  \brief  calculates the log likelihood
        function logL( self ) result ( L )
            use simple_math
            class(emgaufit), intent(inout) :: self
            integer :: i, j
            real    :: L, S
            L = 0.
            do i=1,self%N
                S = 0.
                do j=1,self%K
                    S = S+self%mix(j)*gaussian1D(self%dat(i),self%avgs(j),sqrt(self%vars(j)))
                end do
                L = L+log(S)
            end do
        end function

        !>  \brief  evaluates the responsibilities (gammas) using the current parameter values
        subroutine estep( self )
            use simple_math
            class(emgaufit), intent(inout) :: self
            real    :: P(self%K), S
            integer :: i, j
            do i=1,self%N
                S = 0.
                do j=1,self%K
                    P(j) = self%mix(j)*gaussian1D(self%dat(i),self%avgs(j),sqrt(self%vars(j)))
                    S = S+P(j)
                end do
                do j=1,self%K
                    self%gammas(i,j) = P(j)/S
                end do
            end do
        end subroutine

        !>  \brief  updates the statistics based on the new responsibilities
        subroutine mstep( self )
            class(emgaufit), intent(inout) :: self
            real    :: NK, dev
            integer :: i, j
            do j=1,self%K
                ! calculate average
                self%avgs(j) = 0.
                NK = 0.
                do i=1,self%N
                    self%avgs(j) = self%avgs(j)+self%dat(i)*self%gammas(i,j)
                    NK = NK+self%gammas(i,j)
                end do
                self%avgs(j) = self%avgs(j)/NK
                ! calculate variance
                self%vars(j) = 0.
                do i=1,self%N
                    dev = self%dat(i)-self%avgs(j)
                    self%vars(j) = self%vars(j)+self%gammas(i,j)*dev*dev
                end do
                self%vars(j) = self%vars(j)/NK
                ! calculate mixing coefficients
                self%mix(j) = NK/real(self%N)
            end do
        end subroutine

        ! DESTRUCTOR

        !>  \brief  is a destructor
        subroutine kill( self )
            class(emgaufit), intent(inout) :: self
            if( self%exists )then
                deallocate(self%avgs, self%vars, self%mix, self%gammas)
                self%exists = .false.
            endif
        end subroutine

        ! UNIT TEST

        !>  \brief  is the unit test for the class
        subroutine test_emgaufit
            use simple_rnd,  only: gasdev
            use simple_stat, only: moment
            use gnufor2
            real           :: mydata(1000), avgs(2), vars(2)
            real           :: ave, sdev, var
            integer        :: i
            logical        :: err
            type(emgaufit) :: emfit
            ! generate the data
            do i=1,1000
                if( i <= 500 )then
                    mydata(i) = gasdev(1.,1.)
                else
                    mydata(i) = gasdev(5.,2.)
                endif
            end do
            call moment(mydata(:500), ave, sdev, var, err)
            write(*,*) 'FACIT AVG/VAR 1:', ave, var
            call moment(mydata(501:), ave, sdev, var, err)
            write(*,*) 'FACIT AVG/VAR 2:', ave, var
            ! fit
            call emfit%new(1000,2)
            call emfit%set_data(mydata)
            avgs = 1.5
            avgs = 2.5
            call emfit%fit(50, avgs)
            avgs = emfit%get_avgs()
            vars = emfit%get_vars()
            write(*,*) 'AVG/VAR 1:', avgs(1), vars(1)
            write(*,*) 'AVG/VAR 2:', avgs(2), vars(2)
        end subroutine

        ! Source https://www.mathworks.com/help/stats/hierarchical-clustering.html#bq_679x-10
        subroutine hierarc_clust(vec,thresh,labels)
          real,    intent(in)  :: vec(:)    ! input vector
          real,    intent(in)  :: thresh    ! threshold for class merging
          integer, intent(out) :: labels(:) ! labels of the elements in vec
          real, allocatable :: mat(:,:)     ! matrix that contains the distances
          integer :: N                      ! number of data points
          real    :: d
          integer :: i, j, index(1), loc1(1), loc2(1), cnt, ncls
          logical, allocatable :: mask(:), outliers(:)
          real,    allocatable :: centroids(:)
          N = size(vec)
          ! 1) calc all the couples of distances, using euclid dist
          allocate(mat(N,N), source = 0.)
          do i = 1, N-1
            do j = i+1, N
              mat(i,j) = abs(vec(i)-vec(j))
              mat(j,i) = mat(i,j)
            enddo
          enddo
          ! print *, 'mat: '
          ! call vis_mat(mat)
          ! 2) Generate binary clusters
          allocate(mask(N),source = .true. )
          allocate(outliers(N), source = .false.)
          ncls = 0
          do i = 1, N ! fin the class for vec(i)
            if(mask(i)) then ! if it's not already been clustered
              mask(i) = .false.
              ! find the index of the couple
              d = minval(mat(i,:), mask)
              index(:) = minloc(mat(i,:), mask)
              ncls = ncls + 1
              ! assign lavels
              labels(i) = ncls
              if(d <= thresh) then ! if it's not an outlier (it has a couple)
                labels(index(1)) = labels(i)
                mask(index(1)) = .false. ! index(1) has already been clustered
              else
                outliers(i) = .true.
              endif
            endif
          enddo
          print *, 'vec:     ', nint(vec)
          print *, 'labels:  ', labels
          print *, 'outliers: ', outliers
          ! 3) Calculate centroids
          allocate(centroids(ncls), source = 0.)
          mask = .true. ! reset
          do i = 1, ncls
            ! find the first member of the class
            loc1(:) = minloc(abs(labels-i))
            if(.not. outliers(loc1(1))) then
              mask(loc1(1)) = .false.
              ! find the second member of the class
              loc2(:) = minloc(abs(labels-i), mask)
              mask(loc2(1)) = .false.
              centroids(i) = (vec(loc1(1)) + vec(loc2(1)))/2.
            else ! the class has just one element
              loc1(:) = minloc(abs(labels-i))
              centroids(i) = vec(loc1(1))
              mask(loc1(1)) = .false.
            endif
          enddo
          print *, 'centroids: ', centroids
          mask  = .true. ! reset
          ! 4) merge classes
          do i = 1, ncls-1
            do j = i+1, ncls
              if(abs(centroids(i)-centroids(j)) <= thresh) then ! merge classes
                ! update labels
                where(labels == j) labels = i
              endif
            enddo
          enddo
          print *, 'after merging labels:', labels
          ! 5) Reoder labels
          cnt = 0
          do i = 1, ncls
             if(any(labels== i)) then !there is a class labelled i
                cnt = cnt + 1                  !increasin order cc
                where(labels == i) labels = cnt
              endif
          enddo
          print *, 'after ordering labels:', labels
          ! 6) recalculate centroids
          deallocate(centroids)
          ncls = maxval(labels) ! the nb of classes is maxval(labels)
          allocate(centroids(ncls), source = 0.)
          mask = .true. ! reset
          do i = 1, ncls
            cnt = count(labels == i) ! counts the nb of elements in the class
            ! find the all the cnt member of the class and update the centroids
            do j = 1, cnt
              loc1(:) = minloc(abs(labels-i), mask)
              mask(loc1(1)) = .false. ! do not consider this member of the class anymore
              centroids(i) = centroids(i)+ vec(loc1(1))
            enddo
            centroids(i) = centroids(i)/real(cnt)
          enddo
          print *, 'final centroids: ', centroids
        end subroutine hierarc_clust

   !  ! This subroutine performs laplacian filtering on the input image.
   !  subroutine laplacian_filt(self)
   !      type(image), intent(inout) :: self
   !      real    :: k3(3,3,3), k2(3,3) !laplacian kernels (2D-3D)
   !      integer :: ldim(3)
   !      ldim = self%get_ldim()
   !      k2 = (1./8.)*reshape([0.,1.,0.,1.,-4., 1., 0., 1., 0.], [3,3])
   !      k3 = (1./12.)*reshape([0.,0.,0., 0.,1.,0., 0.,0.,0.,&
   !      &                     0.,1.,0., 1.,-6.,1., 0.,1.,0.,0.,0.,0., 0.,1.,0., 0.,0.,0.], [3,3,3])
   !      if(ldim(3) .ne. 1) then
   !          call self%imfilter(k3)
   !      else
   !          call self%imfilter(k2)
   !      endif
   !  end subroutine laplacian_filt
   !
   !  subroutine generate_distribution(fname_coords_pdb, ldim, smpd)
   !      use simple_atoms, only: atoms
   !         character(len=*), intent(in) :: fname_coords_pdb
   !         integer,          intent(in) :: ldim(3)
   !         real,             intent(in) :: smpd
   !         type(atoms) :: atom
   !         type(image) :: vol !simulated distribution
   !         real        :: cutoff
   !         ! Generate distribution based on atomic position
   !         call vol%new(ldim,smpd)
   !         cutoff = 8.*smpd
   !         call atom%new(fname_coords_pdb)
   !         call atom%convolve(vol, cutoff)
   !         call vol%write('Convoluted.mrc')
   ! end subroutine generate_distribution
end module simple_test_chiara_try_mod

    program simple_test_chiara_try
       include 'simple_lib.f08'
       use simple_math
       use simple_image, only : image
       use simple_atoms, only : atoms
       use simple_nanoparticles_mod, only: nanoparticle
       use simple_segmentation
       use simple_bin_cc_image, only : binimage
       use simple_test_chiara_try_mod
       real          :: thresh ! threshold for class definition, user inputted
       integer                  :: i, j, ncls
       integer, parameter       :: DIM = 10
       ! real     :: avg, stdev
       real, allocatable :: vector(:)
       integer  :: labels(DIM), centers(DIM)
       integer  :: cnt(DIM), not_found
       real     :: avg_within(DIM), stdev_within(DIM)



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
     ! call img2%write('gaussian1DIMG.mrc')
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
