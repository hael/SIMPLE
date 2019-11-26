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
       use simple_nanoparticle, only: nanoparticle
       use simple_segmentation
       use simple_binimage, only : binimage
       use simple_test_chiara_try_mod
       real          :: thresh ! threshold for class definition, user inputted
       integer                  :: i, j, ncls
       integer, parameter       :: DIM = 10
       real, allocatable :: vector(:), vector_read(:)
       integer  :: labels(DIM), centers(DIM)
       integer  :: io_stat, filnum
       real     :: avg_within(DIM), stdev_within(DIM)
       character(len=2) :: element
       character(len=6) :: pdbout
       type(atoms) :: atom
       type(nanoparticle) :: nano
       element = 'pt'
       pdbout = 'PDBout'
       call nano%new( '../'//'71-4_NP65A.mrc', 0.358, element)
       call nano%set_atomic_coords('71-4_NP65A_atom_centers.pdb')
       call nano%set_img('71-4_NP65ACC.mrc', 'img_cc')
       call atom%new('71-4_NP65A_atom_centers.pdb')
       call nano%center_on_atom('71-4_NP65A_atom_centers.pdb',pdbout )

 end program simple_test_chiara_try
