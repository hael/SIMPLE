module simple_test_chiara_try_mod
    include 'simple_lib.f08'
    use simple_image, only : image
    use simple_binimage, only : binimage

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

   ! TO GENERATE IMG CONTAINING JUST GRAPHENE
   ! ! This subroutine operates on img_spec on the basis of the
   ! ! provided logical mask. For each pxl in the msk, it subsitutes
   ! ! its value with the avg value in the neighbours of the pxl that
   ! ! do not belong to the msk. The neighbours are in a shpere of radius 3 pxls.
   ! subroutine filter_peaks(fname,nimage,lmsk)
   !   character(len=*), intent(in) :: fname
   !   integer,          intent(in) :: nimage ! img identifier in the stack
   !   integer,          intent(in) :: lmsk(:,:,:)
   !   real,    pointer   :: rmat(:,:,:)
   !   complex, pointer   :: cmat(:,:,:)
   !   real,    parameter :: RAD = 4.
   !   type(image) :: img_spec
   !   real        :: avg, smpd, rabs, ang
   !   integer     :: nptcls, i, j, ii, jj, ldim(3), cnt
   !   integer     :: hp,k,kp,phys(3),h
   !
   !
   !   call find_ldim_nptcls( fname, ldim, nptcls, smpd)
   !   ldim(3) = 1
   !   if(nimage < 1 .or. nimage > nptcls) stop !sanity check
   !   call img_spec%new(ldim, smpd)
   !   call img_spec%read(fname, nimage)
   !   call img_spec%norm()
   !   call img_spec%fft()
   !
   !   call img_spec%get_cmat_ptr(cmat)
   !   do i = ldim(1)/2+1, ldim(1)
   !      do j = 1, ldim(2)
   !        if(lmsk(i,j,1)==0) then
   !          h = i-(ldim(1)/2+1)
   !          k = j-(ldim(2)/2+1)
   !          phys = img_spec%comp_addr_phys(h,k,0)
   !          cmat(phys(1),phys(2),1) = cmplx(0.,0.)
   !        endif
   !      enddo
   !    enddo
   !    call img_spec%ifft()
   !    call img_spec%write('graphene.mrc')
   !    call img_spec%kill
   ! end subroutine filter_peaks
end module simple_test_chiara_try_mod

    program simple_test_chiara_try
       include 'simple_lib.f08'
       use simple_math
       use simple_image, only : image
       use simple_binimage, only : binimage
       use simple_test_chiara_try_mod
       use simple_nano_utils
       real,    allocatable :: points1(:,:), points2(:,:), ppoints1(:,:), dist(:), dist_sq(:),dist_close(:)
       logical, allocatable :: mask(:)
       integer :: i, j, N_min, N_max
       integer :: cnt, cnt2, cnt3, location(1)
       real    :: smpd, theoretical_radius, rmsd!,avg, m(3), tmp_max
       real, parameter :: CLOSE_THRESH = 1.
       real :: U(3,3), r(3), lrms
       theoretical_radius = 1.1
       smpd = 1.
       N_min = 5
       N_max = 7
       allocate(points1(3,N_max), points2(3,N_min), source  = 0.)
       do i = 1, N_max
         points1(:,i) = [0.8*real(i),2.*real(i)+0.1,real(i)+0.2]
         if(i <= N_min) points2(:3,i) = [0.77*real(i),2.*real(i)+0.02,real(i)+0.54]
       enddo
       allocate(dist(N_max), dist_sq(N_max), source = 0.)
       allocate(dist_close(N_max), source = 0.) ! there are going to be unused entry of the vector
       allocate(mask(N_min), source = .true.)
       cnt  = 0
       cnt2 = 0
       cnt3 = 0
       do i = 1, N_max !compare based on centers1
           if(cnt2+cnt3+1 <=N_min) then ! just N_min couples, starting from 0
               dist(i) = pixels_dist(points1(:,i),points2(:,:),'min',mask,location, keep_zero = .true.)
               if(dist(i)*smpd > 2.*theoretical_radius) then
                   dist(i) = 0.
                   cnt = cnt + 1
               elseif(dist(i)*smpd <= CLOSE_THRESH) then
                   cnt3 = cnt3 + 1
                   dist_close(i) = dist(i)**2
               elseif(dist(i)*smpd > CLOSE_THRESH .and. dist(i)*smpd<=2.*theoretical_radius ) then  !to save the atoms which correspond with a precision in the range [0,2*theoretical_radius] pm
                   cnt2 = cnt2 + 1
                   mask(location(1)) = .false. ! not to consider the same atom more than once
               endif
               dist_sq(i) = dist(i)**2 !formula wants them square
            endif
       enddo
       rmsd = sqrt(sum(dist_sq)/real(count(abs(dist_sq) > TINY)))
       write(*,*) 'RMSD CALCULATED CONSIDERING ALL ATOMS = ', rmsd*smpd, ' A'
       write(*,*) 'RMSD ATOMS THAT CORRESPOND WITHIN 1 A = ', (sqrt(sum(dist_close)/real(count(dist_close > TINY))))*smpd, ' A'
       ! Kabsch testing
       allocate(ppoints1(3,5), source = points1(:3,:5))
       call kabsch(ppoints1,points2,U,r,lrms)
       print *, 'Rotation matrix: '
       call vis_mat(U)
       print *, 'Translation vec: ', r
       print *, 'RMSD: ', lrms
       print *, 'coords after kabsch, ppoints1'
       call vis_mat(ppoints1)
       print *, 'points2'
       call vis_mat(points2)

       ! test function find_couples
       ! real :: P(3,6), Q(3,8), U(3,3), r(3), lrms
       ! real, allocatable :: PP(:,:), QQ(:,:)
       ! character(len=2) :: element
       ! element='PT'
       ! P = reshape([1.,6.,0.,2.,5.,9.,3.,4.,8.,4.,3.,7.,5.,2.,6.,6.,1.,5.],[3,6])
       ! Q = reshape([1.,6.,0.,2.,5.,9.,3.,4.,8.,4.,3.,7.,5.,2.,6.,6.,1.,5.,5.,6.,7.,8.,9.,0.],[3,8])
       ! print *, 'P: '
       ! call vis_mat(P)
       ! print *, 'Q: '
       ! call vis_mat(Q)
       ! call find_couples(P,Q,element,0.358,PP,QQ)
       ! print *, 'shape(PP): ', shape(PP)
       ! print *, 'shape(QQ): ', shape(QQ)
       ! print *, 'PP: '
       ! call vis_mat(PP)
       ! print *, 'QQ: '
       ! call vis_mat(QQ)


       ! real :: A(4,2), sig(4), copyA(4,2)
       ! type(image)    :: raw_img
       ! type(image) :: img_spec
       ! real :: ave, sdev, maxv, minv
       ! real, parameter :: UP_LIM=2.14, L_LIM = 1.42 !1.23
       ! real, allocatable :: x(:)
       ! integer :: UFlim, LFlim, lims(2), i, h, k, sh, cnt, loc(3)
       ! integer, parameter :: BOX = 160
       ! real :: smpd, thresh, radius, m1(3)
       ! integer :: center(3)
       ! real, pointer :: rmat(:,:,:)
       ! real, allocatable :: rmat_aux(:,:,:)
       ! integer, allocatable :: lmsk(:,:,:), lmat(:,:,:)
       !


       ! ! test for removing graphene peaks
       ! smpd=0.358
       ! call img_spec%new([BOX,BOX,1], smpd)
       ! call raw_img%new([BOX,BOX,1], 1.)
       ! do i = 1, 1000
       !     call progress(i,1000)
       !     call img_spec%read('NP_4_background_pspec.mrc',i)
       !     call raw_img%read('NP_4.mrc',i)
       !     call remove_graphene_peaks(raw_img,img_spec)
       !     call raw_img%write('RemovedPeaksNOroavgMASKED.mrc',i)
       ! enddo
       ! call img_spec%kill
       ! call raw_img%kill


       ! center(1:2) = BOX/2 + 1
       ! center(3)   = 1
       ! call img_spec%new([BOX,BOX,1],smpd)
       ! call img_spec%get_rmat_ptr(rmat)
       ! UFlim = calc_fourier_index(max(2.*smpd,UP_LIM), BOX,smpd)! Everything at a resolution higher than UP_LIM A is discarded
       ! lims(1) = UFlim-3
       ! lims(2) = UFlim+3
       ! cnt = 0
       ! do i = 1, BOX
       !     do j = 1, BOX
       !         h   = -int(BOX/2) + i - 1
       !         k   = -int(BOX/2) + j - 1
       !         sh  =  nint(hyp(real(h),real(k)))
       !         if(sh < lims(2) .and. sh > lims(1)) then
       !             cnt = cnt + 1                                !number of pixels in the selected zone
       !             rmat(i,j,1) = real(cnt)
       !         else
       !             rmat(i,j,1) = 0.
       !         endif
       !      enddo
       ! enddo
       ! print *, 'Nb of pxls in the selected zone: ', cnt
       ! call img_spec%write('OrigianalImg.mrc')
       ! ! Find maximum
       ! loc = maxloc(rmat(1:BOX,1:BOX,:))
       ! print *, 'location of the maximum :', loc, 'value: ', maxval(rmat(1:BOX,1:BOX,1))
       ! ! Calculate the radius
       ! radius = euclid(real(loc), real(center))
       ! print *, 'Maximum value identified at radius: ', radius
       ! !  Translate into the origin
       ! loc = loc - center
       ! m1(3) = 1
       ! m1(2) = sin(PI/3.)*loc(2)
       ! m1(1) = cos(PI/3.)*loc(1)
       ! print *, 'm1 : ', m1
       ! ! Come back
       ! m1(1:2)  = m1(1:2)  + real(center(1:2))
       ! loc(1:2) = loc(1:2) + center(1:2)
       ! call img_aux%new([BOX,BOX,1], smpd)
       ! rmat_aux = img_aux%get_rmat()
       ! rmat_aux(nint(m1(1)), nint(m1(2)),1) = 20.
       ! rmat_aux(loc(1), loc(2), 1) = 1.
       ! call img_aux%set_rmat(rmat_aux)
       ! call img_aux%write('RotatedPoint.mrc')
       ! call img_spec%kill
       ! call img_aux%kill
       ! ! Do they have the same dist from the origin?
       ! print *, 'euclid(loc, center):', euclid(real(loc), real(center))
       ! print *, 'euclid(m1, center):', euclid(m1, real(center))
       !
       ! stop
       ! x = pack(rmat(1:BOX,1:BOX,1), abs(rmat(1:BOX,1:BOX,1)) > 0.) ! pack where it's not 0
       ! call otsu(x, thresh)
       ! print *, 'selected threshold : ', thresh
       ! ! report results on the image
       ! where(rmat(1:BOX,1:BOX,1)>=thresh)
       !     rmat(1:BOX,1:BOX,1) = 1.
       ! elsewhere
       !     rmat(1:BOX,1:BOX,1) = 0.
       ! endwhere
       ! call img_spec%write('BinarizedZoneOtsu.mrc')


        ! TRASLATE  INTO SUBROUTINE
        ! call img_spec%new([BOX,BOX,1],smpd)
        ! call img_spec%read('RotationalAvgStk.mrc', 1)
        ! call img_spec%get_rmat_ptr(rmat)
        ! UFlim = calc_fourier_index(max(2.*smpd,UP_LIM), BOX,smpd)! Everything at a resolution higher than UP_LIM A is discarded
        ! lims(1) = UFlim-3
        ! lims(2) = UFlim+3
        ! cnt = 0
        ! do i = 1, BOX
        !     do j = 1, BOX
        !         h   = -int(BOX/2) + i - 1
        !         k   = -int(BOX/2) + j - 1
        !         sh  =  nint(hyp(real(h),real(k)))
        !         if(sh < lims(2) .and. sh > lims(1)) then
        !             cnt = cnt + 1                                !number of pixels in the selected zone
        !         else
        !             rmat(i,j,1) = 0.
        !         endif
        !      enddo
        ! enddo
        ! print *, 'Nb of pxls in the selected zone: ', cnt
        ! call img_spec%write('ROAVGSelectedZone.mrc')
        ! x = pack(rmat(1:BOX,1:BOX,1), abs(rmat(1:BOX,1:BOX,1)) > 0.) ! pack where it's not 0
        ! call otsu(x, thresh)
        ! print *, 'selected threshold : ', thresh
        ! ! report results on the image
        ! allocate(lmat(1:BOX,1:BOX), source = .false.)
        ! where(rmat(1:BOX,1:BOX,1)>=thresh)  lmat = .true.
        !
        ! call img_spec%write('ROAVGBinarizedZoneOtsu.mrc')
        !
        ! ! Do the same thing for the other band
        ! deallocate(x)
        ! call img_spec%read('RotationalAvgStk.mrc', 1) ! reset
        ! call img_spec%get_rmat_ptr(rmat)
        ! LFlim = calc_fourier_index(max(2.*smpd,L_LIM), BOX,smpd)! Everything at a resolution higher than UP_LIM A is discarded
        ! lims(1) = LFlim-3
        ! lims(2) = LFlim+3
        ! cnt = 0
        ! do i = 1, BOX
        !     do j = 1, BOX
        !         h   = -int(BOX/2) + i - 1
        !         k   = -int(BOX/2) + j - 1
        !         sh  =  nint(hyp(real(h),real(k)))
        !         if(sh < lims(2) .and. sh > lims(1)) then
        !             cnt = cnt + 1                                !number of pixels in the selected zone
        !         else
        !             rmat(i,j,1) = 0.
        !         endif
        !      enddo
        ! enddo
        ! print *, 'Nb of pxls in the selected zone: ', cnt
        ! call img_spec%write('ROAVGSelectedZone2.mrc')
        ! x = pack(rmat(1:BOX,1:BOX,1), abs(rmat(1:BOX,1:BOX,1)) > 0.) ! pack where it's not 0
        ! call otsu(x, thresh)
        ! print *, 'selected threshold : ', thresh
        ! ! report results on the image
        ! where(rmat(1:BOX,1:BOX,1)>=thresh)  lmat = .true.
        !
        ! call img_spec%write('ROAVGBinarizedZone2Otsu.mrc')

        ! mask merge generation



        ! call img_spec%stats( ave=ave, sdev=sdev, maxv=maxv, minv=minv)!,mskimg=mask_img)
        ! call img_spec%binarize(ave+2.5*sdev)
        ! call img_spec%write('BinarizedBIN.mrc')




       ! SVD tries
       ! A = reshape([1.,3.,5.,7.,2.,4.,6.,8.], [4,2])
       ! copyA = A
       ! print *, 'A'
       ! call vis_mat(A)
       ! call svdcmp(A,w,v)
       ! print *, 'u'
       ! call vis_mat(A)
       ! print *, 'w', w
       ! print *, 'v'
       ! call vis_mat(v)
       ! d = A(:,2)
       ! print *, 'd', d
       ! t = matmul(d,copyA)
       ! t1  = minval(t)
       ! t2  = maxval(t)
       ! print *, 't', t
       ! print *, 't1', t1
       ! print *, 't2', t2
 end program simple_test_chiara_try

 ! SVD tries
 ! A = reshape([1.,3.,5.,7.,2.,4.,6.,8.], [4,2])
 ! copyA = A
 ! print *, 'A'
 ! call vis_mat(A)
 ! call svdcmp(A,w,v)
 ! print *, 'u'
 ! call vis_mat(A)
 ! print *, 'w', w
 ! print *, 'v'
 ! call vis_mat(v)
 ! d = A(:,2)
 ! print *, 'd', d
 ! t = matmul(d,copyA)
 ! t1  = minval(t)
 ! t2  = maxval(t)
 ! print *, 't', t
 ! print *, 't1', t1
 ! print *, 't2', t2


! GRAPHENE PEAKS
! real, parameter :: UP_LIM=2.14, L_LIM = 1.23
! call img_spec%new([BOX,BOX,1], smpd)
! call img_spec%read('RotationalAvgStk.mrc',1)
! call raw_img%new([BOX,BOX,1], 1.)
! call raw_img%read('NP841_background_nn2.mrcs',1)
! call remove_graphene_peaks(raw_img,img_spec,'RemovedPeaksUTILS.mrc')
! call img_spec%kill
! call raw_img%kill
