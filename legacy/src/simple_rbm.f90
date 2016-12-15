!==Class simple_RBM
!
! RBM is a Restricted Bolzamnn Machine. Later we must implement support for having the weight matrix and hidden
! layers shared with the corresponding sigmodial layer of a MLP (Multi-Layer Perceptron) network. This is to be 
! able to use the RBM as a building block for the DBN (Deep Belief Network). The class defines the parameters of 
! the model along with basic operations for inferring hidden from visible (and vice-versa), as well as for 
! performing CD (Contrastive Divergence) updates.The code is distributed with the hope that it will be useful, 
! but _WITHOUT_ _ANY_ _WARRANTY_. Redistribution or modification is regulated by the GNU General Public License. 
! *Author:* Hans Elmlund, 2012-06-21.
!
!==Changes are documented below
!
module simple_RBM
use simple_rnd    ! singleton
use simple_jiffys ! singleton
use simple_stat   ! singleton
implicit none

public :: RBM
private

type :: RBM
    private
    integer               :: n_vis=0           !< nr of visible
    integer               :: n_hid=0           !< nr of hidden
    integer               :: N=0               !< nr of images
    real, allocatable     :: v0(:,:)           !< copy  of input data vector
    real, allocatable     :: W(:,:)            !< weight matrix
    real, allocatable     :: hbias(:,:)        !< hidden unit biases
    real, allocatable     :: vbias(:,:)        !< visible unit biases
    real, allocatable     :: h0_sigm_act(:,:)  !< sigmoid activation for hidden layer 0
    real, allocatable     :: h1_sigm_act(:,:)  !< sigmoid activation for hidden layer 1
    real, allocatable     :: v1_sigm_act(:,:)  !< sigmoid activation for visible layer 1
    real, allocatable     :: h0_sample(:,:)    !< Bernoulli vector for hidden layer 0
    real, allocatable     :: v1_sample(:,:)    !< Bernoulli vector for visible layer 1
    logical               :: existence=.false. !< to indicate existence
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! GIBBS SAMPLING & CONTRASTIVE DIVERGENCE
    procedure          :: generate
    procedure          :: master
    procedure          :: init
    procedure          :: RBMupdate
    procedure          :: sample
    procedure          :: sample_rev
    procedure          :: update
    ! METRICS
    procedure          :: energy
    procedure          :: free_energy
    ! GETTERS/SETTERS
    procedure          :: get_hlayer_sz
    procedure          :: get_v1_sigm_act
    procedure          :: set_v1_sigm_act
    procedure          :: get_h1_sigm_act
    procedure          :: set_h1_sigm_act
    procedure          :: get_vbias
    procedure          :: set_vbias
    procedure          :: get_hbias
    procedure          :: set_hbias
    procedure          :: get_vsample
    procedure          :: get_hsample
    ! I/O
    procedure          :: write
    procedure          :: read
    procedure, private :: serialize
    procedure, private :: unserialize
    ! DESTRUCTOR
    procedure          :: kill
end type

interface RBM
    module procedure constructor
end interface

!integer :: cnt=0

contains

    ! CONSTRUCTORS
    
    !>  \brief  is a constructor
    function constructor( n_vis, n_hid ) result( self )
        integer, intent(in) :: n_vis, n_hid
        type(RBM)           :: self
        call self%new( n_vis, n_hid )
    end function
    
    !>  \brief  is a constructor
    subroutine new( self, n_vis, n_hid )
        class(RBM), intent(inout) :: self
        integer, intent(in)       :: n_vis, n_hid
        integer                   :: alloc_stat
        self%n_vis = n_vis
        self%n_hid = n_hid
        allocate( self%v0(self%n_vis,1), self%W(self%n_hid,self%n_vis), self%hbias(self%n_hid,1),&
        self%vbias(self%n_vis,1), self%h0_sigm_act(self%n_hid,1), self%h1_sigm_act(self%n_hid,1),&
        self%v1_sigm_act(self%n_vis,1), self%h0_sample(self%N_hid,1), self%v1_sample(self%n_vis,1),&
        stat=alloc_stat )
        call alloc_err('new; simple_RBM', alloc_stat)
        self%v0          = 0.
        self%W           = 0.
        self%hbias       = 0.
        self%vbias       = 0.
        self%h0_sigm_act = 0.
        self%h1_sigm_act = 0.
        self%v1_sigm_act = 0.
        self%h0_sample   = 0.
        self%v1_sample   = 0.
        self%existence = .true.
    end subroutine
    
    ! GIBBS SAMPLING & CONTRASTIVE DIVERGENCE
    
    !>  \brief  is for sampling the generative model at a given image index
    function generate( self, datastk, recsz, i, AVG ) result( dat )
        class(RBM), intent(inout)    :: self
        character(len=*), intent(in) :: datastk
        integer, intent(in)          :: recsz, i
        real, intent(in)             :: AVG(self%n_vis)
        real, allocatable            :: dat(:)
        integer                      :: funit, file_stat
        open(unit=funit, status='old', action='read', file=datastk,&
        access='direct', form='unformatted', recl=recsz, iostat=file_stat)
        call fopen_err('generate; simple_RBM', file_stat)
        ! read data vec
        read(funit, rec=i) self%v0(:,1)
        call normalize_bin(self%v0(:,1))
        ! calculate activations
        call self%sample
        dat = self%get_v1_sigm_act()
        close(unit=funit)
    end function
    
    !>  \brief  doing it all
    subroutine master( self, datastk, recsz, N, maxits, eps_in )
        class(RBM), intent(inout)    :: self
        character(len=*), intent(in) :: datastk
        integer, intent(in)          :: recsz, N, maxits
        real, intent(in), optional   :: eps_in
        integer                      :: k, file_stat, funit, i
        real                         :: e, e_prev, x, eps
        eps = 0.1 ! default learning rate 
        if( present(eps_in) ) eps = eps_in
        write(*,'(A)') '>>> OPTIMIZING THE RBM'
        funit = get_fileunit()
        open(unit=funit, status='old', action='read', file=datastk,&
        access='direct', form='unformatted', recl=recsz, iostat=file_stat)
        call fopen_err('master; simple_RBM', file_stat)
        call self%init
        e = huge(x)
        k = 0
        do
            k = k+1
            e_prev = e
            e = 0.
            do i=1,N
                call progress(i,N)
                ! read data vec
                read(funit, rec=i) self%v0(:,1)
                call  normalize_bin(self%v0(:,1))
                ! update model
                call self%RBMupdate(eps)
                ! update energy fun
                e = e+self%energy()
            end do
            if( abs(e-e_prev) < 0.1 ) exit
!            if( k == 1 .or. mod(k,5) == 0 )then
                write(*,"(1XA,1XI3,1XA,1XF10.0)") 'Iteration:', k, 'Energy:', e
!            endif
            if( k == maxits ) exit
        end do
        close(unit=funit)
    end subroutine
    
    !>  \brief  is for initialization
    subroutine init( self )
        class(RBM), intent(inout) :: self
        real    :: l1, l2, delta
        integer :: i, j
        ! initialize W with uniform random samples from 
        ! -4.*sqrt(6./(n_visible+n_hidden)) and 4.*sqrt(6./(n_hidden+n_visible))
        l2 = 4.*sqrt(6./real(self%n_vis+self%n_hid))
        l1 = -l2
        delta = l2-l1
        do i=1,self%n_hid
            do j=1,self%n_vis
                self%W(i,j) = ran3()*delta+l1
            end do
        end do
        self%W = 0.
        ! initialize the biases with zero
        self%hbias = 0.
        self%vbias = 0.
    end subroutine
    
    !>  \brief  is for updating the RBM
    subroutine RBMupdate( self, eps )
        class(RBM), intent(inout)  :: self
        real, intent(in)           :: eps
        call self%sample
        call self%update(eps)
    end subroutine
    
    !>  \brief  is for sampling the nonlinear manifold provided by the RBM
    subroutine sample( self )
        class(RBM), intent(inout) :: self
        integer                   :: i, j
        real                      :: s, ran
        !$omp parallel default(shared) private(i,j,s)
        ! LOOP1, calculate h0 given v0 (input layer)
        ! Q=distribution of hidden layer
        !$omp workshare
        self%h0_sample = 0.
        !$omp end workshare
        !$omp do schedule(static)
        do i=1,self%n_hid
            s = 0.
            do j=1,self%n_vis
                s = s+self%W(i,j)*self%v0(j,1)
            end do
            self%h0_sigm_act(i,1) = sigm(self%hbias(i,1)+s)
            ! Bernoulli vector generation, here we interpret the sigmod activation 
            ! as the probability of success in a coin flipping experiment
            if( ran3() <= self%h0_sigm_act(i,1) ) self%h0_sample(i,1) = 1.
        end do
        !$omp end do
        ! LOOP2: calculate v1 given fresh h0
        ! P=distribution of visible layer
        !$omp workshare
        self%v1_sample = 0.
        !$omp end workshare
        !$omp do schedule(static) 
        do j=1,self%n_vis
            s = 0.
            do i=1,self%n_hid
                s = s+self%W(i,j)*self%h0_sigm_act(i,1)
            end do
            self%v1_sigm_act(j,1) = sigm(self%vbias(i,1)+s)
            ! Bernoulli vector generation, here we interpret the sigmod activation 
            ! as the probability of success in a coin flipping experiment
            if( ran3() <= self%v1_sigm_act(j,1) ) self%v1_sample(j,1) = 1.
        end do
        !$omp end do
        ! LOOP3: same as first, but now calculating h1 given fresh v1
        !$omp do schedule(static) 
        do i=1,self%n_hid
            s = 0.
            do j=1,self%n_vis   
                s = s+self%W(i,j)*self%v1_sigm_act(j,1)
            end do
            self%h1_sigm_act(i,1) = sigm(self%hbias(i,1)+s)
        end do
        !$omp end do nowait
        !$omp end parallel
    end subroutine
    
    !>  \brief  for sampling the nonlinear manifold provided by the RBM in reverse (required in the DBN)
    subroutine sample_rev( self, h1in )
        class(RBM), intent(inout) :: self
        real, intent(in)          :: h1in(self%n_hid)
        integer                   :: i, j
        real                      :: s, ran
        real                      :: Wt(self%n_vis,self%n_hid)
        !$omp parallel default(shared) private(i,j,s)
        !$omp workshare
        Wt = transpose(self%W)
        !$omp end workshare
        ! LOOP1, calculate v1 given h1
        !$omp do schedule(static) 
        do i=1,self%n_vis
            s = 0.
            do j=1,self%n_hid
                s = s+Wt(i,j)*h1in(j)
            end do
            self%v1_sigm_act(i,1) = sigm(self%vbias(i,1)+s)
        end do
        !$omp end do nowait
        !$omp end parallel
    end subroutine

    !>  \brief  This is the RBM update procedure for binomial units. It also works
    !!          for exponential and truncated exponential units, and for the linear parameters
    !!          of a Gaussian unit (using the appropriate procedure for Q and P). It can
    !!          readily be adapted for the variance parameter of Gaussian units.
    subroutine update( self, eps )
        class(RBM), intent(inout)  :: self
        real, intent(in)           :: eps
        integer                    :: i, j
        real                       :: s
        ! CONTRASTIVE DIVERGENCE UPDATES
        ! Weight matrix
        self%W = self%W+eps*matmul(self%h0_sample,transpose(self%v0))-&
        matmul(self%h1_sigm_act,transpose(self%v1_sigm_act))
        ! Hidden unit biases
        self%hbias = self%hbias+eps*(self%h0_sample-self%h1_sigm_act)
        ! Visible unit biases
        self%vbias = self%vbias+eps*(self%v0-self%v1_sample)
    end subroutine
    
    ! METRICS
    
    !>  \brief  for calculating the RBM free energy contribution for one data input
    !!          (take the negative logarithm of the summed contributions)
    function free_energy( self ) result( f )
        class(RBM), intent(inout) :: self
        real :: f
        f = exp(-energy(self))
    end function
    
    !>  \brief  for calculating the RBM energy for one data input
    !!          E(v,h) = - vbias'v - hbias'h - h'Wv
    function energy( self ) result( e )
        class(RBM), intent(inout) :: self
        integer :: i, j
        real    :: e !tmp(1,1)
        e = 0.
        do j=1,self%n_vis
            e = e-self%vbias(j,1)*self%v1_sigm_act(j,1)
        end do
        do i=1,self%n_hid
            e = e-self%hbias(i,1)*self%h1_sigm_act(i,1)
        end do
        do j=1,self%n_vis
            do i=1,self%n_hid
                e = e-self%W(i,j)*self%h1_sigm_act(i,1)*self%v1_sigm_act(j,1)
            end do
        end do
!        tmp = -matmul(transpose(self%vbias),self%v1_sigm_act)&
!        -matmul(transpose(self%hbias),self%h1_sigm_act)&
!        -matmul(transpose(self%h1_sigm_act),matmul(self%W,self%v1_sigm_act))
!        e = tmp(1,1)
    end function
    
    ! GETTERS/SETTERS
    
    !>  \brief  is a getter
    function get_hlayer_sz( self ) result( sz )
        class(RBM), intent(in) :: self
        integer :: sz
        inquire(iolength=sz) self%h1_sigm_act(:,1)
    end function
    
    !>  \brief  is a getter
    function get_v1_sigm_act( self ) result( v1_sigm_act )
        class(RBM), intent(in) :: self
        real, allocatable      :: v1_sigm_act(:)
        allocate( v1_sigm_act(self%n_vis), source=self%v1_sigm_act(:,1) )
    end function
    
    !>  \brief  is a setter
    subroutine set_v1_sigm_act( self, v1_sigm_act )
        class(RBM), intent(inout) :: self
        real, intent(in)          :: v1_sigm_act(self%n_vis)
        self%v1_sigm_act(:,1) = v1_sigm_act
    end subroutine
    
    !>  \brief  is a getter
    function get_h1_sigm_act( self ) result( h1_sigm_act )
        class(RBM), intent(in) :: self
        real, allocatable      :: h1_sigm_act(:)
        allocate( h1_sigm_act(self%n_hid), source=self%h1_sigm_act(:,1) )
    end function
    
    !>  \brief  is a setter
    subroutine set_h1_sigm_act( self, h1_sigm_act )
        class(RBM), intent(inout) :: self
        real, intent(in)          :: h1_sigm_act(self%n_hid)
        self%h1_sigm_act(:,1) = h1_sigm_act
    end subroutine
    
    !>  \brief  is a getter
    function get_hbias( self ) result( hbias )
        class(RBM), intent(in) :: self
        real, allocatable      :: hbias(:)
        allocate( hbias(self%n_hid), source=self%hbias(:,1) )
    end function
    
    !>  \brief  is a setter
    subroutine set_hbias( self, hbias )
        class(RBM), intent(inout) :: self
        real, intent(in)          :: hbias(self%n_hid)
        self%hbias(:,1) = hbias
    end subroutine
    
    !>  \brief  is a getter
    function get_vbias( self ) result( vbias )
        class(RBM), intent(in) :: self
        real, allocatable      :: vbias(:)
        allocate( vbias(self%n_vis), source=self%vbias(:,1) )
    end function
    
    !>  \brief  is a setter
    subroutine set_vbias( self, vbias )
        class(RBM), intent(inout) :: self
        real, intent(in)          :: vbias(self%n_vis)
        self%vbias(:,1) = vbias
    end subroutine
    
    !>  \brief  is a getter
    function get_vsample( self ) result( v1_sample )
        class(RBM), intent(in) :: self
        real, allocatable      :: v1_sample(:)
        real                   :: ran
        integer                :: j
        allocate( v1_sample(self%n_vis), source=self%v1_sample(:,1) )
    end function
    
    !>  \brief  is a getter
    function get_hsample( self ) result( h0_sample )
        class(RBM), intent(in) :: self
        real, allocatable      :: h0_sample(:)
        real                   :: ran
        integer                :: i
        allocate( h0_sample(self%n_hid), source=self%h0_sample(:,1) )
    end function
    
    ! I/O
    
    !>  \brief  for writing an RBM to disk
    subroutine write( self, fname )
        class(RBM), intent(in)       :: self
        character(len=*), intent(in) :: fname
        real, allocatable            :: series(:)
        integer                      :: recsz, handle, file_stat
        series = self%serialize()
        inquire( iolength=recsz ) series
        handle = get_fileunit()
        open( unit=handle, file=fname, status='replace', iostat=file_stat,&
        access='direct', action='write', form='unformatted', recl=recsz )
        call fopen_err( 'write; simple_RBM', file_stat )
        write(unit=handle,rec=1) series
        close(unit=handle) 
        deallocate( series )
    end subroutine
    
    !>  \brief  for reading an RBM from disk
    subroutine read( self, fname )
        class(RBM), intent(inout)    :: self
        character(len=*), intent(in) :: fname
        real, allocatable            :: series(:)
        integer                      :: recsz, handle, file_stat
        series = self%serialize()
        inquire( iolength=recsz ) series
        handle = get_fileunit()
        open( unit=handle, file=fname, status='old', iostat=file_stat,&
        access='direct', action='read', form='unformatted', recl=recsz )
        call fopen_err( 'read; simple_RBM', file_stat )
        read(unit=handle,rec=1) series
        close(unit=handle) 
        call self%unserialize( series )
        deallocate( series )
    end subroutine
    
    !>  \brief  for serializing an RBM
    function serialize( self ) result( series )
        class(RBM), intent(in) :: self
        real, allocatable :: series(:)
        integer :: cnt, alloc_stat, sz, i, j
        ! count the number of reals in (1) weight matrix (2) hidden biases 
        ! (3) visible biases (4) h1_sigm_act (5) v1_sigm_act
        ! allocate and store in series
        sz = self%n_hid*self%n_vis+2*self%n_hid+2*self%n_vis 
        allocate( series(sz), stat=alloc_stat )
        call alloc_err('serialize; simple_RBM', alloc_stat)
        cnt = 0
        do i=1,self%n_hid
            do j=1,self%n_vis
                cnt = cnt+1
                series(cnt) = self%W(i,j)
            end do
        end do
        do i=1,self%n_hid
            cnt = cnt+1
            series(cnt) = self%hbias(i,1)
        end do
        do i=1,self%n_vis
            cnt = cnt+1
            series(cnt) = self%vbias(i,1)
        end do
        do i=1,self%n_hid
            cnt = cnt+1
            series(cnt) = self%h1_sigm_act(i,1)
        end do
        do i=1,self%n_vis
            cnt = cnt+1
            series(cnt) = self%v1_sigm_act(i,1)
        end do
    end function
    
    !>  \brief  for unpacking a serialiized RBM
    subroutine unserialize( self, series )
        class(RBM), intent(inout) :: self
        real, intent(in)          :: series(:)
        integer :: sz, i, j, cnt
        sz = self%n_hid*self%n_vis+2*self%n_hid+2*self%n_vis
        if( sz /= size(series) ) stop 'Incompatible series; unserialize; simple_rbm'
        cnt = 0
        do i=1,self%n_hid
            do j=1,self%n_vis
                cnt = cnt+1
                self%W(i,j) = series(cnt)
            end do
        end do
        do i=1,self%n_hid
            cnt = cnt+1
            self%hbias(i,1) = series(cnt) 
        end do
        do i=1,self%n_vis
            cnt = cnt+1
            self%vbias(i,1) = series(cnt)
        end do
        do i=1,self%n_hid
            cnt = cnt+1
            self%h1_sigm_act(i,1) = series(cnt)
        end do
        do i=1,self%n_vis
            cnt = cnt+1
            self%v1_sigm_act(i,1) = series(cnt)
        end do
    end subroutine
    
    ! DESTRUCTOR
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(RBM), intent(inout) :: self
        if( self%existence )then
            deallocate(self%v0, self%W, self%hbias, self%vbias, self%h0_sigm_act,&
            self%h1_sigm_act, self%v1_sigm_act, self%h0_sample, self%v1_sample)
            self%existence = .false.
        endif
    end subroutine

end module simple_RBM