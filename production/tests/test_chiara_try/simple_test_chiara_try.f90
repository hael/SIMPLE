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

   subroutine fcc_lattice(element,moldiam,box,smpd,m)
     use simple_atoms, only: atoms
     character(len=4), intent(inout) :: element
     integer,          intent(in)    :: box
     real,             intent(in)    :: moldiam,smpd
     real,             intent(inout) :: m(3) ! center of mass of the generated fcc
     type(image)      :: vol
     type(atoms)      :: atoms_obj
     character(len=2) :: el1, el2
     character(len=8) :: crystal_system
     real             :: center(3),ha,x,x1,x2,x3,y,y1,y2,y3,z,z1,z2,z3,msksq,cutoff
     real             :: a ! lattice parameter
     integer          :: i,j,k,n,ncubes
     if(.not.is_even(box))then
         print *, 'BOX must be even'
         stop
     endif
     el1 = element(1:2)
     el2 = element(3:4)
     msksq = (moldiam/2.)**2.
     m = 0.
     call vol%new([box,box,box], smpd)
     select case(uppercase(trim(element)))
         case('C')
             a = 3.567 ! diamond
         case('SI')
             a = 5.431020511
         case('GE')
             a = 5.658
         case('AL')
             a = 4.046
         case('NI')
             a = 3.499
         case('CU')
             a = 3.597
         case('PT')
             a = 3.912
         case('AU')
             a = 4.065
         case('AG')
             a = 4.079
         case('PD')
             a = 3.859
         case('PB')
             a = 4.920
         case('FE')
             a = 2.856; crystal_system = 'bcc     '
         case('MO')
             a = 3.142; crystal_system = 'bcc     '
         case('W')
             a = 3.155; crystal_system = 'bcc     '
         case('PBSE')
             a = 6.12;  crystal_system = 'rocksalt'
         case DEFAULT
             a = 3.76
     end select
     ha = a/2.
     ncubes = floor(real(box) * smpd / a)
     ! atoms at edges
     n = 0
     center = (real([box,box,box]/2 + 1)-1.)*smpd
      ! count
      do i=1,ncubes
          x  = real(i-1)*a
          x1 = x+ha; x2 = x; x3 = x+ha
          do j=1,ncubes
              y  = real(j-1)*a
              y1 = y+ha; y2 = y+ha; y3 = y
              do k=1,ncubes
                  z  = real(k-1)*a
                  z1 = z; z2 = z+ha; z3 = z+ha
                  ! edge
                  if( sum(([x ,y ,z ]-center)**2.) <= msksq ) n = n+1
                  ! faces
                  if( sum(([x1,y1,z1]-center)**2.) <= msksq ) n = n+1
                  if( sum(([x2,y2,z2]-center)**2.) <= msksq ) n = n+1
                  if( sum(([x3,y3,z3]-center)**2.) <= msksq ) n = n+1
              enddo
          enddo
      enddo
      ! generate atoms object
      call atoms_obj%new(n)
      do i = 1,n
          call atoms_obj%set_name(i,element//'  ')
          call atoms_obj%set_element(i,el1)
      enddo
      ! generate all coordinates
      n = 0
      do i=1,ncubes
          x  = real(i-1)*a ;x1 = x+ha; x2 = x; x3 = x+ha
          do j=1,ncubes
              y  = real(j-1)*a; y1 = y+ha; y2 = y+ha; y3 = y
              do k=1,ncubes
                  z  = real(k-1)*a; z1 = z; z2 = z+ha; z3 = z+ha
                  ! edge
                  if( sum(([x ,y ,z ]-center)**2.) <= msksq )then
                       n = n+1; call atoms_obj%set_coord(n, [x,y,z])
                       m = m + [x,y,z]
                  endif
                  ! faces
                  if( sum(([x1,y1,z1]-center)**2.) <= msksq )then
                      n = n+1; call atoms_obj%set_coord(n, [x1,y1,z1])
                      m = m + [x1,y1,z1]
                  endif
                  if( sum(([x2,y2,z2]-center)**2.) <= msksq )then
                      n = n+1; call atoms_obj%set_coord(n, [x2,y2,z2])
                      m = m + [x2,y2,z2]
                  endif
                  if( sum(([x3,y3,z3]-center)**2.) <= msksq )then
                      n = n+1; call atoms_obj%set_coord(n, [x3,y3,z3])
                      m = m + [x3,y3,z3]
                  endif
              enddo
          enddo
      enddo
      ! center of gravity
      m = m/real(n)
      ! write PDB
      do i=1,n
          call atoms_obj%set_num(i,i)
          call atoms_obj%set_resnum(i,i)
          call atoms_obj%set_chain(i,'A')
          call atoms_obj%set_occupancy(i,1.)
      enddo
      call atoms_obj%writePDB('FccLattice')
      ! for convolution
      cutoff = 3.*a
      call atoms_obj%convolve(vol,cutoff)
      call vol%write('FccLattice.mrc')
   end subroutine fcc_lattice

   subroutine read_table( fname, vals )
    character(len=*),  intent(in)  :: fname    !< input filename
    real, allocatable, intent(out) :: vals(:,:) !< array of values
    integer :: nl, funit, iline,io_stat
    nl = nlines(trim(fname))
    call fopen(funit,fname,'old','unknown',io_stat)
    call fileiochk("read_table failed to open file "//trim(fname),io_stat )
    allocate( vals(nl,2), stat=alloc_stat )
    if(alloc_stat /= 0) call allocchk ('In: read_filetable; simple_fileio  ', alloc_stat)
    do iline=1,nl
        read(funit,*) vals(iline,1), vals(iline,2)
    end do
    call fclose(funit,io_stat)
    call fileiochk("read_filetable failed to close",io_stat)
end subroutine read_table
end module simple_test_chiara_try_mod

program simple_test_chiara_try
   include 'simple_lib.f08'
   use simple_image, only : image
   use simple_nanoparticle_utils, only : dock_nanosPDB
   use simple_atoms, only : atoms
   use simple_dock_coords
   character(len=100) :: pdbfile1, pdbfile2, pdbfile_out1, pdbfile_out2, pdbfile
   real, allocatable :: model(:,:)
   real :: theta,t(3),a(3), trs_in, rot_trans(7), U(3,3), r(3), r_com(3)
   type(atoms) :: a_ref, a_targ
   integer :: i, n


   ! In /Users/chiara/Desktop/Nanoparticles/StrainAnalysis/TestLatticeFit
   pdbfile1= 'vol_Nov26_atom_centers.pdb'
   pdbfile2= 'vol_Nov28_atom_centers.pdb'
   pdbfile_out1 = 'vol_26Nov_atom_centers_docked'
   pdbfile_out2 = 'vol_28Nov_atom_centers_docked'
   !In /Users/chiara/Desktop/Nanoparticles/MapsToDOCK/TestDockPDBFiles
   ! pdbfile1= 'particle1_atom_centers.pdb'
   ! pdbfile2= 'particle1_c3sym_atom_centers.pdb'
   ! pdbfile_out1 = 'particle1_atom_centers_docked'
   ! pdbfile_out2 = 'particle1_c3sym_atom_centers_docked'
   call dock_nanosPDB(pdbfile1,pdbfile2,pdbfile_out1,pdbfile_out2)


end program simple_test_chiara_try
