!==Class simple_DBN
!
! DBN is a Deep Belief Network. The class defines the parameters of the model along with basic operations for 
! inferring hidden from visible (and vice-versa), as well as for performing CD (Contrastive Divergence) updates.
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. Redistribution 
! or modification is regulated by the GNU General Public License. 
! *Author:* Hans Elmlund, 2012-06-21.
!
!==Changes are documented below
!
module simple_DBN
use simple_RBM,      only: RBM
use simple_ran_tabu, only: ran_tabu 
use simple_jiffys    ! singleton
use simple_defs      ! singleton
implicit none

public :: DBN
private

type :: DBN
    private
    type(ran_tabu)                     :: rt                !< for random image sampling
    integer, allocatable               :: validimgs(:)      !< indices of the validation images
    integer                            :: nlayers=0         !< nr of layers
    integer                            :: N=0               !< nr of images
    integer                            :: nsample=0         !< nr of images to sample in each epoch
    integer                            :: nvalid=0          !< nr of images for validating the model (and check convergence)
    integer, allocatable               :: nunits(:)         !< nr of units per hidden layer
    type(RBM), allocatable             :: layers(:)         !< the layers are composed of RBMs
    character(len=STDLEN)              :: datastk=''        !< input, data stack
    integer                            :: D=0               !< dimension of data
    character(len=STDLEN), allocatable :: fnames(:)         !< hidden layer filenames (images)
    character(len=STDLEN), allocatable :: vnames(:)         !< hidden layer filenames (validation images)
    character(len=STDLEN), allocatable :: rbmnames(:)       !< file names for flushing the optimized RBMs
    integer, allocatable               :: fhandles(:)       !< automatic handles (images)
    integer, allocatable               :: vhandles(:)       !< automatic handles (validation images)
    logical                            :: existence=.false. !< to indicate existence
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! I/O
    procedure          :: open
    procedure          :: close
    procedure, private :: get_img_from_file
    procedure, private :: get_h1_img_from_file
    procedure, private :: get_h1_valid_from_file
    ! CD & GIBBS
    procedure          :: train
    procedure          :: sample
    ! GETTERS   
    procedure          :: get_vsample
!    procedure          :: get_hsample
    procedure          :: get_feature
    procedure          :: get_v1_sigm_act
    ! DESTRUCTOR
    procedure          :: kill
end type

! make structure for selecting data and validation data

interface DBN
    module procedure constructor
end interface

contains

    ! CONSTRUCTORS
    
    !>  \brief  is a constructor
    function constructor( datastk, N, D, nunits, nvalid, nsample ) result( self )
        character(len=*), intent(in)  :: datastk
        integer, intent(in)           :: N, D
        integer, intent(in)           :: nunits(:)
        integer, intent(in), optional :: nvalid, nsample
        type(DBN)                     :: self
        if( present(nvalid) )then
            if( present(nsample) )then
                call self%new(datastk, N, D, nunits, nvalid, nsample)
            else
                call self%new(datastk, N, D, nunits, nvalid)
            endif
        else
             if( present(nsample) )then
                call self%new(datastk, N, D, nunits, nsample)
            else
                call self%new(datastk, N, D, nunits)
            endif
        endif
    end function
    
    !>  \brief  is a constructor
    
    subroutine new( self, datastk, N, D, nunits, nvalid, nsample )
        class(DBN), intent(inout)     :: self
        character(len=*), intent(in)  :: datastk
        integer, intent(in)           :: N, D
        integer, intent(in)           :: nunits(:)
        integer, intent(in), optional :: nvalid, nsample
        integer                       :: alloc_stat, i, file_stat, ran
        character(len=32)             :: dig
        self%datastk = datastk
        self%N       = N
        self%D       = D 
        self%nlayers = size(nunits)
        self%nvalid = 0
        if( present(nvalid) ) self%nvalid  = nvalid
        self%nsample = nint(real(self%N-self%nvalid)*0.8)
        if( present(nsample) ) self%nsample = min(self%nsample,nsample)
        allocate( self%nunits(self%nlayers+1), self%validimgs(nvalid), self%layers(self%nlayers),&
        self%fnames(self%nlayers), self%vnames(self%nlayers), self%fhandles(self%nlayers),&
        self%vhandles(self%nlayers), self%rbmnames(self%nlayers), stat=alloc_stat )
        call alloc_err('new; simple_DBN', alloc_stat)
        self%nunits    = 0
        self%validimgs = 0
        self%fnames    = ''
        self%vnames    = ''
        self%fhandles  = 0
        self%vhandles  = 0
        self%rbmnames  = ''
        self%rt = ran_tabu(N)
        if( nvalid /= 0 )then
            do i=1,nvalid
                ran = self%rt%irnd()
                call self%rt%insert(ran)
                self%validimgs(i) = ran
            end do
        endif
        ! set nr of units
        self%nunits(1) = D
        do i=1,self%nlayers
            self%nunits(i+1) = nunits(i)
        end do
        ! make rbms and set filenames
        do i=1,self%nlayers
            self%layers(i) = RBM(self%nunits(i),self%nunits(i+1))
            write(dig,*) i
            self%fnames(i) = 'h_sigm_act_layer_'//trim(adjustl(dig))//'.bin'
            self%vnames(i) = 'h_sigm_act_layer_valid_'//trim(adjustl(dig))//'.bin'
            self%rbmnames(i) = 'rbm_'//trim(adjustl(dig))//'.bin'
        end do
        self%existence = .true.
    end subroutine
    
    !>  \brief  opens the files storing data and hidden layer activations for all images
    subroutine open( self )
        class(DBN), intent(inout) :: self
        integer :: i, file_stat
        do i=1,self%nlayers
            self%fhandles(i) = get_fileunit()
            open( unit=self%fhandles(i), file=self%fnames(i),  status='unknown', iostat=file_stat,&
            access='direct', action='readwrite', form='unformatted', recl=self%layers(i)%get_hlayer_sz() )
            call fopen_err( 'open_files; simple_DBN, 1', file_stat )
            self%vhandles(i) = get_fileunit()
            open( unit=self%vhandles(i), file=self%vnames(i),  status='unknown', iostat=file_stat,&
            access='direct', action='readwrite', form='unformatted', recl=self%layers(i)%get_hlayer_sz() )
            call fopen_err( 'open_files; simple_DBN, 2', file_stat )
        end do
    end subroutine
    
    !>  \brief  closes the files storing the hidden layer activations for all images
    subroutine close( self )
        class(DBN), intent(inout) :: self
        integer :: i
        close(self%datastk_unit)
        do i=1,self%nlayers
            close(self%fhandles(i))
            close(self%vhandles(i))
        end do
    end subroutine
    
    ! NOT NEEDED THE RBM READS ITS OWN IMAGE
     !>  \brief  is a getter, assumes files open
!    function get_img_from_file( self, j ) result( img )
!        class(DBN), intent(in) :: self
!        integer, intent(in)    :: j
!        real, allocatable      :: img(:)
!        integer                :: alloc_stat
!        allocate( img(self%nunits(1)), stat=alloc_stat )
!        call alloc_err('get_img_from_file; simple_DBN', alloc_stat)
!        read(unit=self%datastk_unit,rec=j) img
!    end function
    
    !>  \brief  is a getter, assumes files open
    function get_h1_img_from_file( self, i, j ) result( img )
        class(DBN), intent(in) :: self
        integer, intent(in)    :: i, j ! i is layer, j is image
        real, allocatable      :: img(:)
        integer                :: alloc_stat
        allocate( img(self%nunits(i+1)), stat=alloc_stat )
        call alloc_err('get_h1_img_from_file; simple_DBN', alloc_stat)
        read(unit=self%fhandle(i),rec=j) img
    end function
    
    !>  \brief  is a getter, assumes files open
    function get_h1_valid_from_file( self, i, j ) result( img )
        class(DBN), intent(in) :: self
        integer, intent(in)    :: i, j  ! i is layer, j is image
        real, allocatable      :: img(:)
        integer                :: alloc_stat
        allocate( img(self%nunits(i+1)), stat=alloc_stat )
        call alloc_err('get_h1_valid_from_file; simple_DBN', alloc_stat)
        read(unit=self%vhandle(i),rec=j) img
    end function
    
    !>  \brief  this is the master training subroutine. To avoid over-fitting the generative DBN, a fraction 
    !!          of the data is selected for validation only. In each iteration, we monitor a correlation and 
    !!          a free energy that tells us how well the generative model generalizes to data never used for 
    !!          training. The routine assumes that the files are open
    subroutine train( self, eps_in, MAXITS_in )
        class(DBN), intent(inout)     :: self
        real, intent(in), optional    :: eps_in
        integer, intent(in), optional :: MAXITS_in
        real                          :: eps, energy_img, energy_valid
        integer                       :: MAXITS, it, ran, i, j
        logical                       :: converged
        eps = 1.    ! default learning rate
        if( present(eps_in) ) eps = eps_in
        MAXITS = 15 ! default maximum iterations (per RBM) 
        if( present(MAXITS_in) ) MAXITS = MAXITS_in
        ! loop over layers
        do i=1,self%nlayers
            call self%layers(i)%init ! random initialization
            if( i > 1 )then
                ! take the hidden biases from below as the visible biases here
                call self%layers(i)%set_vbias(self%layers(i-1)%get_hbias())
            endif
            converged = .false.
            it = 0
            do while( .not. converged )
                it = it+1
                energy_img = 0.
                energy_valid = 0.
                ! Process nsample images
                do j=1,self%nsample
                    call progress(j, self%nsample)
                    ! Image sampling
                    ran = self%rt%irnd()
                    call self%rt%insert(ran)
                    ! RBM update
                    if( i == 1 )then
                        call self%layers(i)%RBMupdate(ind=ran, eps_in=eps)
                    else
                        call self%layers(i)%RBMupdate(v0in=self%get_h1_img_from_file(i-1,ran), eps_in=eps)
                    endif
                    ! Write hidden layer to file
                    write(unit=self%fhandle(i),rec=ran) self%layers(i)%get_h1_sigm_act()
                    ! Metric update
                    energy_img = energy_img+self%layers(i)%energy()
                end do
                if( self%nvalid > 0 )then
                    ! Process nvalid validation images
                    do j=1,self%nvalid
                        ! Gibbs sampling
                        if( i == 1 )then
                            call self%layers(i)%sample(ind=self%validimgs(j))
                        else
                            call self%layers(i)%sample(v0in=self%get_h1_img_from_file(i-1,self%validimgs(j)))    
                        endif
                        ! Write hidden layer to file
                        write(unit=self%vhandle(i),rec=self%validimgs(j)) self%layers(i)%get_h1_sigm_act()
                        ! Metric update
                        energy_valid = energy_valid+self%layers(i)%energy()
                    end do
                endif
                ! check if RBM converged
                energy_img = energy_img/real(self%nsample)
                energy_valid = energy_valid/real(self%nsample)
                write(*,*) '***************************************'
                write(*,*) 'ITERATION:', it
                write(*,*) 'ENERGY VALID/IMG:', energy_valid, energy_img
                if( it == MAXITS ) converged = .true. ! only MAXITS 4 now
                ! recreate rt
                call self%rt%reset
                if( self%nvalid > 0 )then
                    do j=1,self%nvalid
                        call self%rt%insert(self%validimgs(i))
                    end do
                endif
            end do
            ! write RBM to file
            call self%layers(i)%write( self%rbmnames(i) )
        end do
    end subroutine
    
    !>  \brief  is for sampling the stacked generative DBN model
    !!          assumes that files are open for reading
    subroutine sample( self, j )
        class(DBN), intent(inout) :: self
        integer, intent(in)       :: j ! image index
        integer                   :: i
        do i=self%nlayers,1,-1 ! reverse loop
            if( i == self%nlayers )then
                ! get the top level hidden layer activations from file
                ! these are the most abstract features
                call self%layers(i)%sample_rev(self%get_h1_img_from_file(i,j))
            else
                ! use the visible layer activations obtained by reversed sampling 
                ! from the above layer as hidden layer activations for this layer
                call self%layers(i)%sample_rev(self%layers(i+1)%get_v1_sigm_act())
            endif
        end do
        ! the bottom layer RBM visible activations now defines the generative model
        ! Use get_vsample( self, j ) (binary image) or  get_v1_sigm_act ([0,1] normalized
        ! continuous image) to create an image
    end subroutine
    
    ! GETTERS
    
    !>  \brief  is for sampling the input level RBM, assumes that the DBN is sampled
    function get_vsample( self, j ) result( v1_sample )
        class(DBN), intent(in) :: self
        integer, intent(in)    :: j ! image index
        real, allocatable      :: v1_sample(:)
        v1_sample = self%layers(1)%get_vsample()
    end function
    
!    RATHER USELESS
    !>  \brief  is for sampling the input level RBM, assumes that the DBN is sampled
!    function get_hsample( self, j ) result( h1_sample )
!        class(DBN), intent(in) :: self
!        integer, intent(in)    :: j ! image index
!        real, allocatable      :: h1_sample(:)
!        h1_sample = self%layers(1)%get_hsample()
!    end function
    
    !>  \brief  is for getting the the most abstract feature, assumes that the DBN is sampled
    function get_feature( self, j ) result( feat )
        class(DBN), intent(in) :: self
        integer, intent(in)    :: j
        real, allocatable      :: feat(:)
        feat = self%get_h1_img_from_file(self%nlayers,j)
    end function
    
    !>  \brief  is for getting the sigmoid activations of the input level RBM, assumes that the DBN is sampled
    function get_v1_sigm_act( self ) result( v1_sigm_act )
        class(DBN), intent(in) :: self
        real, allocatable      :: v1_sigm_act(:)
        v1_sigm_act = self%layers(1)%get_v1_sigm_act()
    end function
    
    ! DESTRUCTOR
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(DBN), intent(inout) :: self
        integer :: i
        if( self%existence )then
            do i=1,self%nlayers
                call self%layers(i)%kill
            end do
            deallocate( self%nunits, self%validimgs, self%layers,&
            self%fnames, self%vnames, self%fhandles, self%vhandles )
            call self%rt%kill
            self%existence = .false.
        endif
    end subroutine
   
end module simple_DBN