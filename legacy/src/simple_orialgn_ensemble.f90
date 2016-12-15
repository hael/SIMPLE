!>  \brief  SIMPLE orialgn_ensemble class
module simple_orialgn_ensemble
use simple_orialgn_pair, only: orialgn_pair
use simple_oris,         only: oris
use simple_defs          ! singleton
implicit none

public :: orialgn_ensemble
private

type :: orialgn_ensemble
    private
    type(oris), allocatable :: ensemble(:)       !< the alignment ensemble
    type(oris)              :: ref               !< reference 4 alignment
    type(orialgn_pair)      :: oap               !< orialgn_pair object
    real, allocatable       :: cost(:)           !< cost values
    integer                 :: n=0               !< nr of members
    logical                 :: existence=.false. !< to indicate existence
  contains
    ! CONSTRUCTORS
    procedure :: new
    procedure :: construct_ensemble
    ! GETTER/SETTER
    procedure :: get_ori
    procedure :: set_ori
    ! ALIGNMENT
    procedure          :: l1median
    procedure, private :: align_ensemble
    procedure, private :: align_member
    procedure, private :: gen_reference
    ! STATS
    procedure :: find_closest
    procedure :: dist2member
    ! I/O
    procedure :: read
    procedure :: write
    ! DESTRUCTOR
    procedure :: kill
end type

interface orialgn_ensemble
    module procedure constructor
end interface

real, parameter :: TOL=1e-6

contains

    ! CONSTRUCTOR

    !>  \brief  constructor
    function constructor( n ) result( self )
        integer, intent(in)    :: n
        type(orialgn_ensemble) :: self
        call self%new( n )
    end function
    
    !>  \brief  constructor
    subroutine new( self, n )
        use simple_jiffys, only: alloc_err
        class(orialgn_ensemble), intent(inout) :: self
        integer, intent(in) :: n
        integer :: alloc_stat
        self%n = n
        call self%oap%new
        call self%oap%set_reforis(self%ref)
        allocate( self%ensemble(self%n), self%cost(self%n), stat=alloc_stat )
        call alloc_err("In: new; simple_orialgn_ensemble", alloc_stat)
        self%existence = .true.
    end subroutine
    
    !>  \brief  ensemble constructor
    subroutine construct_ensemble( self, m )
        class(orialgn_ensemble), intent(inout) :: self
        integer, intent(in) :: m
        integer :: i 
        do i=1,self%n
            self%ensemble(i) = oris(m)
        end do
    end subroutine
    
    ! GETTER/SETTER
    
    !>  \brief  is a getter
    function get_ori( self, i, j ) result( o )
        use simple_ori, only: ori
        class(orialgn_ensemble), intent(inout) :: self
        integer, intent(in)                    :: i, j
        type(ori)                              :: o
        o = self%ensemble(i)%get_ori(j)
    end function
    
    !>  \brief  is a setter
    subroutine set_ori( self, i, j, o )
        use simple_ori, only: ori
        class(orialgn_ensemble), intent(inout) :: self
        integer, intent(in)                    :: i, j
        class(ori), intent(inout)              :: o
        call self%ensemble(i)%set_ori(j,o)
    end subroutine
    
    ! ALIGNMENT
    
    !>  \brief  generates the l1median alignment from the ensemble
    function l1median( self ) result( median )
        use simple_rnd, only: irnd_uni
        class(orialgn_ensemble), intent(inout) :: self
        type(oris)            :: median
        integer               :: i, iter
        real                  :: x, err, err_prev
        character(len=STDLEN) :: dig
        do i=1,self%n
            if( self%ensemble(i)%get_noris() > 1 )then
                ! alles gut
            else
                stop 'need to read oris into ensemble; l1median; simple_orialgn_ensemble'
            endif
        end do
        ! make initial reference
        self%ref = self%ensemble(1)
        call self%ref%spiral
        call self%ref%rnd_inpls 
        ! refine
        err  = huge(x)
        iter = 0
        do
            iter = iter+1
            write(dig,*) iter
            call self%align_ensemble
            err_prev = err
            err = sum(self%cost)
            write(*,'(a,1x,f7.2)') '>>> ERROR:', err
            call self%gen_reference
            call self%ref%write('consensus_iter_'//trim(adjustl(dig))//'.txt')
            if( abs(err-err_prev) < 1. ) exit
        end do
        median = self%ref
    end function
    
    !>  \brief  aligns all members of the ensemble to the reference
    subroutine align_ensemble( self )
        use simple_jiffys, only: progress
        integer :: i
        class(orialgn_ensemble), intent(inout) :: self
        write(*,'(a)') '>>> ALIGNING ENSEMBLE'
        do i=1,self%n
            call progress(i, self%n)
            call self%align_member( i )
        end do
    end subroutine
    
    !>  \brief  aligns one member of the ensemble to the reference
    subroutine align_member( self, i )
        class(orialgn_ensemble), intent(inout) :: self
        integer, intent(in) :: i
        call self%oap%set_targoris(self%ensemble(i))
        call self%oap%dock(self%cost(i))
        call self%oap%transform_oris(self%ensemble(i))
    end subroutine
    
    !>  \brief  generate a reference from the aligned ensemble
    subroutine gen_reference( self )
        use simple_ori, only: ori
        class(orialgn_ensemble), intent(inout) :: self
        type(oris) :: tmp
        integer    :: noris, i, j
        write(*,'(a)') '>>> CALCULATING THE L1-MEDIAN ALIGNMENT'
        ! make the necessary objects
        noris    = self%ensemble(1)%get_noris()
        self%ref = oris(noris)
        tmp      = oris(self%n)
        ! calculate the median orientations
        do i=1,noris
            do j=1,self%n
                call tmp%set_ori(j,self%ensemble(j)%get_ori(i))
            end do
            call self%ref%set_ori(i,tmp%l1median())
        end do
    end subroutine
    
    ! STATS
    
    !>  \brief  4 finding the closest ensemble member
    subroutine find_closest( self, o, angdist, which )
        class(orialgn_ensemble), intent(inout) :: self
        class(oris), intent(inout)             :: o
        real, intent(out)                      :: angdist
        integer, intent(out)                   :: which
        real                                   :: dist
        integer                                :: i
        angdist = self%ensemble(1).inpldist.o
        which   = 1
        do i=2,self%n
            dist = self%ensemble(i).inpldist.o
            if( dist < angdist )then
                angdist = dist
                which = i
            endif
        end do
    end subroutine
    
    !>  \brief  4 finding the closest ensemble member
    function dist2member( self, o, which ) result( angdist )
        class(orialgn_ensemble), intent(inout) :: self
        class(oris), intent(inout)             :: o
        integer, intent(in)                    :: which
        real                                   :: angdist
        angdist = self%ensemble(which).inpldist.o
    end function
    
    ! I/O
    
    !>  \brief  4 reading an ensemble member
    subroutine read( self, i, orifile )
        use simple_jiffys, only: nlines
        class(orialgn_ensemble), intent(inout) :: self
        integer, intent(in)                    :: i
        character(len=*), intent(in)           :: orifile
        self%ensemble(i) = oris(nlines(orifile))
        call self%ensemble(i)%read(orifile)
    end subroutine
    
    !>  \brief  4 writing an ensemble member
    subroutine write( self, i, orifile )
        use simple_jiffys, only: nlines
        class(orialgn_ensemble), intent(inout) :: self
        integer, intent(in)                    :: i
        character(len=*), intent(in)           :: orifile
        call self%ensemble(i)%write(orifile)
    end subroutine
    
    ! DESTRUCTOR
  
    !>  \brief  destructor
    subroutine kill( self )
        class(orialgn_ensemble), intent(inout) :: self
        integer :: i
        if( self%existence )then
            call self%oap%kill
            do i=1,self%n
                call self%ensemble(i)%kill
            end do 
            deallocate( self%ensemble, self%cost )
            self%existence = .false.
        endif
    end subroutine
  
end module simple_orialgn_ensemble
