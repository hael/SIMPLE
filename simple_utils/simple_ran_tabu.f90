!==Class simple_ran_tabu
!
! simple_ran_tabu contains routines for generation of directed random numbers. The code is distributed 
! with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. Redistribution
! or modification is regulated by the GNU General Public License. *Author:* Hans Elmlund, 2009-05-12.
! 
!==Changes are documented below
!* incorporated in the _SIMPLE_ library, HE 2009-06-25
!* OOD 2003, HE 2012-06-15
!
module simple_ran_tabu
implicit none

type :: ran_tabu
    private
    integer :: NP=0, N_tabus=0 ! integer ranges from 1 to NP
    logical, allocatable :: list(:)
  contains
    procedure :: reset
    procedure :: insert
    procedure :: is
    procedure :: remove
    procedure :: irnd
    procedure :: irnd_pair
    procedure :: irnd_gau
    procedure :: ne_ran_iarr
    procedure :: shuffle
    procedure :: balanced
    procedure :: kill
end type

interface ran_tabu
    module procedure constructor
end interface

contains

    !>  \brief  is a constructor
    function constructor( NP ) result( self )
        use simple_jiffys, only: alloc_err
        integer, intent(in) :: NP
        type(ran_tabu)      :: self
        integer             :: alloc_stat
        call self%kill
        self%NP = NP
        self%N_tabus = 0
        allocate( self%list(NP), stat=alloc_stat )
        call alloc_err('In: new_ran_tabu, module: simple_ran_tabu.f90', alloc_stat)
        self%list = .true. ! all integers from 1 to NP made available
    end function

    !>  \brief  is for clearing the tabu history
    subroutine reset( self )
        class(ran_tabu), intent(inout) :: self
        self%N_tabus = 0
        self%list    = .true. ! all integers from 1 to NP made available
    end subroutine

    !>  \brief  is for insertion of a tabu
    subroutine insert( self, i )
        class(ran_tabu), intent(inout) :: self
        integer, intent(in) :: i
        if( self%list(i) ) then
            self%N_tabus = self%N_tabus+1
            if( self%N_tabus > self%NP ) stop 'nr of tabus larger than NP; insert; simple_ran_tabu'
        endif
        self%list(i) = .false.
    end subroutine

    !>  \brief  is for checking tabu status
    function is( self, i ) result( yep )
        class(ran_tabu), intent(in) :: self
        integer, intent(in)        :: i
        logical :: yep
        if( allocated(self%list) )then
            yep = .not. self%list(i)
        else
            yep = .false.
        endif
    end function

    !>  \brief  is for removal of a tabu
    subroutine remove( self, i )
        class(ran_tabu), intent(inout) :: self
        integer, intent(in) :: i
        if( .not. self%list(i) )then
            self%N_tabus = self%N_tabus-1
        endif
        self%list(i) = .true.
    end subroutine

    !>  \brief  generates a uniform random integer [_1_,_NP_] not in the tabulist, 
    !!          used to direct Monte Carlo search out of forbidden regions.
    function irnd( self ) result( ir )
        use simple_rnd, only: irnd_uni
        class(ran_tabu), intent(in) :: self
        integer                     :: ir
        if( self%N_tabus == self%NP ) stop 'all numbers tabu; irnd; simple_ran_tabu'
        do
            ir = irnd_uni(self%NP)
            if( self%list(ir) ) exit
        end do
    end function

    !>  \brief  generates a random disjoint pair
    subroutine irnd_pair( self, irnd, jrnd )
        class(ran_tabu), intent(in) :: self
        integer, intent(out)        :: irnd, jrnd
        irnd = self%irnd()
        jrnd = irnd
        do while( irnd == jrnd ) 
            jrnd = self%irnd()
        end do
    end subroutine
    
    !>  \brief  generates a normal random integer [_1_,_NP_] not in the tabulist, 
    !!          used to direct Monte Carlo search out of forbidden regions.
    function irnd_gau( self, mean, stdev ) result( irnd )
        use simple_rnd, only: irnd_gasdev 
        class(ran_tabu), intent(in) :: self
        real, intent(in)            :: mean, stdev
        integer                     :: irnd
        if( self%N_tabus == self%NP ) stop 'all numbers tabu; irnd_gau; simple_ran_tabu'
        do
            irnd = irnd_gasdev( mean, stdev, self%NP )
            if( self%list(irnd) ) exit
        end do
    end function
    
    !>  \brief  generates sequence of uniform random numbers [_1_,_NP_] without repetition
    subroutine ne_ran_iarr( self, rndiarr )
        class(ran_tabu), intent(inout) :: self
        integer, intent(out)           :: rndiarr(:)
        integer                        :: i, szrndiarr
        szrndiarr = size(rndiarr)
        if( szrndiarr+self%N_tabus > self%NP ) then
            write( *,* ) 'Random numbers must be generated from a larger set'
            write( *,* ) 'In: ne_ran_iarr, module: simple_ran_tabu.f90'
            stop
        else
            do i=1,szrndiarr
                rndiarr(i) = self%irnd()
                call self%insert(rndiarr(i))  
            end do
        endif
    end subroutine
    
    !>  \brief  shuffles an integer array
    subroutine shuffle( self, shuffled ) 
        class(ran_tabu), intent(inout) :: self
        integer, intent(inout)         :: shuffled(self%NP)
        integer                        :: tmp(self%NP), irnd, i
        call self%reset
        do i=1,self%NP
            irnd = self%irnd()
            tmp(i) = shuffled(irnd)
            call self%insert(irnd)  
        end do
        shuffled = tmp
    end subroutine
    
    !>  \brief  creates a balanced randomised paritioning over nstates states
    subroutine balanced( self, nstates, iarr ) 
        class(ran_tabu), intent(inout) :: self
        integer, intent(in)            :: nstates
        integer, intent(inout)         :: iarr(self%NP)
        integer :: i, s
        i = 0
        do while( i < self%NP )
            do s=1,nstates
                i = i+1
                if( i > self%NP ) exit
                iarr(i) = s
            end do
        end do
        call self%shuffle(iarr)
    end subroutine
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(ran_tabu), intent(inout) :: self
        if( allocated(self%list) ) deallocate( self%list )
    end subroutine

end module simple_ran_tabu
