module simple_autoencoder
use simple_rnd    ! singleton
use simple_jiffys ! singleton
use simple_stat   ! singleton
implicit none

public :: autoencoder
private

type :: autoencoder
    private
    integer           :: nvis=0     !< nr of visible
    integer           :: nhid=0     !< nr of hidden
    real, allocatable :: W1(:,:)    !< weight matrix 1               (nhidxnvis)
    real, allocatable :: W2(:,:)    !< weight matrix 1               (nvisxnhid)
    real, allocatable :: hbias(:,:) !< hidden unit biases            (nhidx1)
    real, allocatable :: vbias(:,:) !< visible unit biases           (nvisx1)
    real, allocatable :: hact(:,:)  !< activations for hidden layer  (nhidx1) (y)
    real, allocatable :: vact(:,:)  !< activations for visible layer (nvisx1) (z)
    real, allocatable :: v0(:,:)    !< copy of input data vec        (nvisx1)
    logical           :: existence=.false. !< to indicate existence
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! CALCULATORS
    procedure          :: generate
    procedure          :: master
    procedure, private :: init
    procedure, private :: calc_hid_act
    procedure, private :: calc_vis_act
    procedure, private :: update
    procedure, private :: energy
end type

contains

    !>  \brief  is a constructor
!    function constructor( nvis, nhid ) result( self )
!        integer, intent(in) :: nvis, nhid
!        type(autoencoder)   :: self
!        call self%new( nvis, nhid )
!    end function

    !>  \brief  is a constructor
    subroutine new( self, nvis, nhid )
        class(autoencoder), intent(inout) :: self
        integer, intent(in)               :: nvis, nhid
        integer                           :: alloc_stat
        self%nvis = nvis
        self%nhid = nhid
        allocate( self%W1(self%nhid,self%nvis), self%W2(self%nvis,self%nhid),&
        self%hbias(self%nhid,1), self%vbias(self%nvis,1), self%hact(self%nhid,1),&
        self%vact(self%nvis,1), self%v0(self%nvis,1),&
        stat=alloc_stat )
        call alloc_err('new; simple_autoencoder', alloc_stat)
        self%W1    = 0.
        self%W2    = 0.
        self%hbias = 0.
        self%vbias = 0.
        self%hact  = 0.
        self%vact  = 0.
        self%existence = .true.
    end subroutine
    
    !>  \brief  is for sampling the generative model at a given image index
    function generate( self, datastk, recsz, i, AVG ) result( dat )
        class(autoencoder), intent(inout) :: self
        character(len=*), intent(in)      :: datastk
        integer, intent(in)               :: recsz, i
        real, intent(in)                  :: AVG(self%nvis)
        real, allocatable                 :: dat(:)
        integer                           :: funit, file_stat
        open(unit=funit, status='old', action='read', file=datastk,&
        access='direct', form='unformatted', recl=recsz, iostat=file_stat)
        call fopen_err('generate; simple_autoencoder', file_stat)
        ! read data vec
        read(funit, rec=i) self%v0(:,1)
        ! calculate activations
        call self%calc_hid_act
        call self%calc_vis_act
        allocate(dat(self%nvis), source=self%hact(:,1))
!        dat = dat+AVG
        close(unit=funit)
    end function
    
    !>  \brief  doing it all
    subroutine master( self, datastk, recsz, N, maxits, eps_in )
        class(autoencoder), intent(inout) :: self
        character(len=*), intent(in)      :: datastk
        integer, intent(in)               :: recsz, N, maxits
        real, intent(in), optional        :: eps_in
        integer                           :: k, file_stat, funit, i
        real                              :: e, e_prev, x, eps
        eps = 0.1 ! default learning rate 
        if( present(eps_in) ) eps = eps_in
        write(*,'(A)') '>>> OPTIMIZING THE AUTOASSOCIATIVE NET'
        funit = get_fileunit()
        open(unit=funit, status='old', action='read', file=datastk,&
        access='direct', form='unformatted', recl=recsz, iostat=file_stat)
        call fopen_err('master; simple_autoencoder', file_stat)
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
                ! calculate activations
                call self%calc_hid_act
                call self%calc_vis_act
                ! update model
                call self%update(eps)
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
        class(autoencoder), intent(inout) :: self
        real    :: l1, l2, delta
        integer :: j, k
        ! initialize W1 with uniform random samples from 
        ! -4.*sqrt(6./(nvisible+nhidden)) and 4.*sqrt(6./(nhidden+n_isible))
        l2 = 4.*sqrt(6./real(self%nvis+self%nhid))
        l1 = -l2
        delta = l2-l1
        do k=1,self%nhid
            do j=1,self%nvis
                self%W1(k,j) = ran3()*delta+l1
            end do
        end do
        self%W2 = transpose(self%W1) ! tied weights
        ! initialize the biases with zero
        self%hbias = 0.
        self%vbias = 0.
    end subroutine

    !>  \brief  is for calculating the hidden layer activations
    subroutine calc_hid_act( self )
        class(autoencoder), intent(inout) :: self
        integer :: k
        self%hact = matmul(self%W1,self%v0)+self%hbias
        do k=1,self%nhid
            self%hact(k,1) = tanh(self%hact(k,1))
        end do
    end subroutine
    
    !>  \brief  is for calculating the visible layer activations
    subroutine calc_vis_act( self )
        class(autoencoder), intent(inout) :: self
        integer :: j
        self%vact = matmul(self%W2,self%hact)+self%vbias
        do j=1,self%nvis
            self%vact(j,1) = tanh(self%vact(j,1))
        end do
    end subroutine
    
    !>  \brief  is the stochastic gradient descent update for the weight matrix 
    subroutine update( self, eps )
        class(autoencoder), intent(inout) :: self
        real, intent(in)                  :: eps
        integer :: j, k
        real :: deltaj(self%nvis), deltak(self%nhid)
        do j=1,self%nvis
            deltaj(j)       = self%vact(j,1)-self%v0(j,1)
            self%vbias(j,1) = self%vbias(j,1)-eps*deltaj(j)
        end do
        do k=1,self%nhid
            deltak(k)       = (1.-self%hact(k,1)**2.)*sum(self%W1(k,:)*deltaj(:))
            self%hbias(k,1) = self%hbias(k,1)-eps*deltak(k) 
        end do
        do k=1,self%nhid 
            do j=1,self%nvis
                self%W1(k,j) = self%W1(k,j)-eps*(deltak(k)*self%v0(j,1)+deltaj(j)*self%hact(k,1))
                self%W2(j,k) = self%W1(k,j)
!                self%W2(j,k) = self%W2(j,k)-eps*(deltaj(j)*self%hact(k,1))
            end do
        end do
    end subroutine
    
    !>  \brief  is the reconstruction error (energy)
    function energy( self ) result( e )
        class(autoencoder), intent(in) :: self
        real :: e
        integer :: j
        e = 0.
        do j=1,self%nvis
            e = e+(self%vact(j,1)-self%v0(j,1))**2.
        end do
    end function 

end module simple_autoencoder