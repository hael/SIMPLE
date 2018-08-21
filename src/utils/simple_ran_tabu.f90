! routines for generation of directed random numbers
module simple_ran_tabu
use simple_defs
use simple_error, only: allocchk, simple_exception
use simple_rnd,   only: multinomal, irnd_gasdev, irnd_uni
implicit none

public :: ran_tabu
private
#include "simple_local_flags.inc"

type :: ran_tabu
    private
    integer :: NP=0       !< integer ranges from 1 to NP
    integer :: N_tabus=0
    logical, allocatable :: avail(:) !< flags for checking availability
  contains
    procedure          :: reset
    procedure          :: insert
    procedure          :: is
    procedure          :: remove
    procedure          :: irnd
    procedure          :: irnd_pair
    procedure          :: irnd_gau
    procedure          :: mnomal
    procedure          :: ne_ran_iarr
    procedure          :: ne_mnomal_iarr
    procedure          :: stoch_nnmat
    procedure, private :: shuffle_1
    procedure, private :: shuffle_2
    generic            :: shuffle => shuffle_1, shuffle_2
    procedure          :: balanced
    procedure          :: kill
end type ran_tabu

interface ran_tabu
    module procedure constructor
end interface ran_tabu

contains

    !>  \brief  is a constructor
    function constructor( NP ) result( self )
        integer, intent(in) :: NP    !< max number of tabus
        type(ran_tabu)      :: self
        call self%kill
        self%NP = NP
        self%N_tabus = 0
        allocate( self%avail(NP), stat=alloc_stat )
        if(alloc_stat /= 0) call allocchk('In: new_ran_tabu, module: simple_ran_tabu.f90', alloc_stat)
        self%avail = .true. ! all integers from 1 to NP made available
    end function constructor

    !>  \brief  is for clearing the tabu history
    subroutine reset( self )
        class(ran_tabu), intent(inout) :: self
        self%N_tabus = 0
        self%avail   = .true. ! all integers from 1 to NP made available
    end subroutine reset

    !>  \brief  is for insertion of a tabu
    subroutine insert( self, i )
        class(ran_tabu), intent(inout) :: self
        integer,         intent(in)    :: i !< input tabu
        if( self%avail(i) ) then
            self%N_tabus = self%N_tabus + 1
            if( self%N_tabus > self%NP ) THROW_HARD('nr of tabus larger than NP; insert')
        endif
        self%avail(i) = .false.
    end subroutine insert

    !>  \brief  is for checking tabu status
    function is( self, i ) result( yep )
        class(ran_tabu), intent(in) :: self
        integer,         intent(in) :: i !< query tabu
        logical :: yep
        if( allocated(self%avail) )then
            yep = .not. self%avail(i)
        else
            yep = .false.
        endif
    end function is

    !>  \brief  is for removal of a tabu
    subroutine remove( self, i )
        class(ran_tabu), intent(inout) :: self
        integer,         intent(in)    :: i  !< remove tabu
        if( .not. self%avail(i) )then
            self%N_tabus = self%N_tabus - 1
        endif
        self%avail(i) = .true.
    end subroutine remove

    !>  \brief  generates a uniform random integer [_1_,_NP_] not tabu,
    !!          used to direct Monte Carlo search out of forbidden regions.
    function irnd( self ) result( ir )
        class(ran_tabu), intent(in) :: self
        integer :: ir
        if( self%N_tabus == self%NP ) THROW_HARD('all numbers tabu; irnd')
        do
            ir = irnd_uni(self%NP)
            if( self%avail(ir) ) return
        end do
    end function irnd

    !>  \brief  generate random disjoint pair
    !! \param irnd,jrnd  random integers
    subroutine irnd_pair( self, irnd, jrnd )
        class(ran_tabu), intent(in)  :: self
        integer,         intent(out) :: irnd, jrnd
        irnd = self%irnd()
        jrnd = irnd
        do while( irnd == jrnd )
            jrnd = self%irnd()
        end do
    end subroutine irnd_pair

    !>  \brief  generates a normal random integer [_1_,_NP_] not tabu,
    !!          used to direct Monte Carlo search out of forbidden regions.
    !!  \param mean Gaussian mean
    !!  \param stdev Gaussian standard deviation
    function irnd_gau( self, mean, stdev ) result( irnd )
        class(ran_tabu), intent(in) :: self
        real,            intent(in) :: mean, stdev
        integer :: irnd
        if( self%N_tabus == self%NP ) THROW_HARD('all numbers tabu; irnd_gau')
        do
            irnd = irnd_gasdev( mean, stdev, self%NP )
            if( self%avail(irnd) ) return
        end do
    end function irnd_gau

    !>  \brief  generates a multinomal random integer not tabu
    function mnomal( self, pvec ) result( irnd )
        class(ran_tabu), intent(in) :: self
        real,            intent(in) :: pvec(self%NP) !< multinomal vector
        integer :: irnd, nrepeats
        if( self%N_tabus == self%NP ) THROW_HARD('all numbers tabu; mnomal')
        nrepeats = 0
        do
            irnd = multinomal(pvec)
            if( self%avail(irnd) )then
                return
            else
                nrepeats = nrepeats + 1
            endif
            if( nrepeats == self%NP )then
                irnd = 0
                return
            endif
        end do
    end function mnomal

    !>  \brief  generates sequence of uniform random numbers [_1_,_NP_] without repetition
    subroutine ne_ran_iarr( self, rndiarr )
        class(ran_tabu), intent(inout) :: self
        integer,         intent(out)   :: rndiarr(:) !< random integer array
        integer :: i, szrndiarr
        szrndiarr = size(rndiarr)
        if( szrndiarr + self%N_tabus > self%NP ) then
            write( *,* ) 'Random numbers must be generated from a larger set:',szrndiarr + self%N_tabus,self%NP
            write( *,* ) 'In: ne_ran_iarr, module: simple_ran_tabu.f90'
            stop
        else
            call self%reset
            do i=1,szrndiarr
                rndiarr(i) = self%irnd()
                call self%insert(rndiarr(i))
            end do
            call self%reset
        endif
    end subroutine ne_ran_iarr

    !>  \brief  generates sequence of uniform random numbers [_1_,_NP_] without repetition
    subroutine ne_mnomal_iarr( self, pvec, rndiarr )
        class(ran_tabu), intent(inout) :: self
        real,            intent(in)    :: pvec(self%NP) !< multinomal vector
        integer,         intent(out)   :: rndiarr(:) !< random integer array
        integer        :: i, szrndiarr, irnd, cnt
        type(ran_tabu) :: rt4shuffle
        szrndiarr = size(rndiarr)
        if( szrndiarr + self%N_tabus > self%NP ) then
            write( *,* ) 'Random numbers must be generated from a larger set'
            write( *,* ) 'In: ne_mnomal_iarr, module: simple_ran_tabu.f90'
            stop
        else
            call self%reset
            cnt = 0
            do i=1,self%NP
                irnd = self%mnomal(pvec)
                if( irnd > 0 )then
                    cnt = cnt + 1
                    rndiarr(cnt) = irnd
                    call self%insert(rndiarr(cnt))
                    if( cnt == szrndiarr ) return
                endif
            end do
            if( cnt <= szrndiarr )then
                do i=cnt,szrndiarr
                    rndiarr(i) = self%irnd()
                    call self%insert(rndiarr(i))
                end do
            endif
            call self%reset
        endif
        ! make sure the order is shuffled
        rt4shuffle = ran_tabu(szrndiarr)
        call rt4shuffle%shuffle(rndiarr)
        call rt4shuffle%kill
    end subroutine ne_mnomal_iarr

    !>  \brief  stochastic nearest neighbor generation
    function stoch_nnmat( self, pfromto, nnn, pmat ) result( nnmat )
        class(ran_tabu), intent(inout) :: self
        integer,         intent(in)    :: pfromto(2), nnn                     !< pmat range
        real,            intent(in)    :: pmat(pfromto(1):pfromto(2),self%NP) !< multinomal array
        integer, allocatable :: nnmat(:,:) !> output nearest neigh matrix
        integer :: iptcl
        allocate(nnmat(pfromto(1):pfromto(2),nnn), stat=alloc_stat)
        if(alloc_stat /= 0) call allocchk('In: simple_ran_tabu; stoch_nnmat', alloc_stat)
        do iptcl=pfromto(1),pfromto(2)
            call self%ne_mnomal_iarr( pmat(iptcl,:), nnmat(iptcl,:))
        end do
    end function stoch_nnmat

    !>  \brief  shuffles an integer array
    subroutine shuffle_1( self, shuffled )
        class(ran_tabu), intent(inout) :: self
        integer,         intent(inout) :: shuffled(self%NP)  !< integer array for shuffling
        integer :: tmp(self%NP), irnd, i
        call self%reset
        do i=1,self%NP
            irnd = self%irnd()
            tmp(i) = shuffled(irnd)
            call self%insert(irnd)
        end do
        shuffled = tmp
        call self%reset
    end subroutine shuffle_1

    !>  \brief  shuffles a real array
    subroutine shuffle_2( self, shuffled )
        class(ran_tabu), intent(inout) :: self
        real,            intent(inout) :: shuffled(self%NP)  !< integer array for shuffling
        real    :: tmp(self%NP)
        integer :: irnd, i
        call self%reset
        do i=1,self%NP
            irnd = self%irnd()
            tmp(i) = shuffled(irnd)
            call self%insert(irnd)
        end do
        shuffled = tmp
        call self%reset
    end subroutine shuffle_2

    !>  \brief  creates a balanced randomised paritioning over nstates states
    subroutine balanced( self, nstates, iarr )
        class(ran_tabu), intent(inout) :: self
        integer,         intent(in)    :: nstates !< num states
        integer,         intent(inout) :: iarr(self%NP) !< integer array for partitioning
        integer :: i, s
        i = 0
        ! random init
        do s = irnd_uni(nstates),nstates
            i = i+1
            if( i > self%NP ) exit
            iarr(i) = s
        end do
        ! elongation
        do while( i < self%NP )
            do s=1,nstates
                i = i+1
                if( i > self%NP ) exit
                iarr(i) = s
            end do
        end do
        call self%shuffle(iarr)
    end subroutine balanced

    !>  \brief  is a destructor
    subroutine kill( self )
        class(ran_tabu), intent(inout) :: self
        if( allocated(self%avail) )then
            deallocate( self%avail, stat=alloc_stat )
            if(alloc_stat /= 0) call allocchk('In: simple_ran_tabu; kill ', alloc_stat)
        end if
    end subroutine kill

end module simple_ran_tabu
