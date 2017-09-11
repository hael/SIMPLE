! clustering of a similarity matrix using stochastic hill-climbing
#include "simple_lib.f08"
module simple_shc_cluster
use simple_defs
use simple_oris, only: oris
use simple_rnd,  only: irnd_uni
use simple_syslib, only: alloc_errchk
implicit none

public :: shc_cluster, test_shc_cluster
private

type shc_cluster
    private
    integer              :: ncls              !< number of clusters
    integer              :: N                 !< number of elements to cluster
    class(oris), pointer :: o_ptr=>null()     !< pointer to orientation data struct
    integer, allocatable :: labels(:)         !< cluster labels
    real                 :: MINS=-1.          !< minimum similarity
    real, pointer        :: S(:,:)            !< pointer to similarity marix
    real, allocatable    :: SPS(:)            !< single-particle similarities
    logical              :: existence=.false. !< to indicate existence
  contains
    procedure          :: new
    procedure          :: shc
    procedure, private :: init
    procedure, private :: move
    procedure, private :: per_ptcl_sim
    procedure, private :: eval_change
    procedure          :: kill
end type shc_cluster

contains

    !>  \brief  is a constructor
    subroutine new( self, N_in, ncls_in, S_in, o, minsim )
        class(shc_cluster), intent(inout) :: self
        integer,            intent(in)    :: N_in, ncls_in
        real, target,       intent(in)    :: S_in(N_in,N_in)
        type(oris), target, intent(in)    :: o
        real, optional,     intent(in)    :: minsim
        call self%kill
        self%N     =  N_in
        self%ncls  =  ncls_in
        self%S     => S_in
        self%o_ptr => o
        if( present(minsim) )then
            self%MINS = minsim
        else
            self%MINS = -1.
        endif
        allocate(self%labels(self%N), self%SPS(self%N), stat=alloc_stat)
        if(alloc_stat /= 0) allocchk('In: nsimple_shc_cluster::new')
        self%labels = 0
        self%SPS    = 0.
        self%existence = .true.
    end subroutine new

    !>  \brief  is the stochastic hill climbing routine for clustering
    subroutine shc( self, doprint, label, sim )
        use simple_ran_tabu, only: ran_tabu
        class(shc_cluster), intent(inout) :: self
        logical,            intent(in)    :: doprint
        character(len=*),   intent(in)    :: label
        real,               intent(out)   :: sim
        integer             :: order(self%N)
        type(ran_tabu)      :: rt
        integer             :: iptcl, iter, convcnt
        real                :: sim_prev, sim_new
        integer, parameter  :: MAXITS=1000, CONVFAC=15
        rt = ran_tabu(self%N)
        call self%init
        convcnt = 0
        if( doprint ) write(*,*) 'SHC CLUSTERING ITERATION: ', 0, 'SIMILARITY: ', sum(self%SPS)/real(self%N)
        do iter=1,MAXITS
            ! make random order
            call rt%reset
            call rt%ne_ran_iarr(order)
            ! calculate previous similarity
            sim_prev = sum(self%SPS)/real(self%N)
            ! execute moves
            do iptcl=1,self%N
                call self%move(order(iptcl))
            end do
            ! calculate new similarity
            sim_new = sum(self%SPS)/real(self%N)
            if( doprint )then
                write(*,*) 'SHC CLUSTERING ITERATION: ', iter, 'SIMILARITY: ', sim_new
            endif
            if( sim_prev/sim_new > 0.995 )then
                convcnt = convcnt+1
            else
                convcnt = 0
            endif
            if( convcnt == CONVFAC ) exit
        end do
        ! set ouput classes
        do iptcl=1,self%N
            call self%o_ptr%set(iptcl, label, real(self%labels(iptcl)))
        end do
        call rt%kill
        sim = sim_new
    end subroutine shc

    !>  \brief  initialises randomly and calculates per-particle similarities
    subroutine init( self )
        class(shc_cluster), intent(inout) :: self
        integer :: iptcl
        ! random labels
        do iptcl=1,self%N
            self%labels(iptcl) = irnd_uni(self%ncls)
        end do
        ! per-particle similarities
        do iptcl=1,self%N
            call self%per_ptcl_sim(iptcl)
        end do
    end subroutine init

    !>  \brief  executes a move and updates the label and the per-particle similarities
    subroutine move( self, iptcl )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(shc_cluster), intent(inout) :: self
        integer,            intent(in)    :: iptcl
        real    :: new_sim
        integer :: cls, new_cls, jptcl
        cls = self%labels(iptcl)
        new_cls = cls
        do while( new_cls == cls )
            new_cls = irnd_uni(self%ncls)
        end do
        new_sim = self%eval_change(iptcl, new_cls)
        if( new_sim >= self%SPS(iptcl) )then
            self%labels(iptcl) = new_cls
            self%SPS(iptcl)    = new_sim
            !$omp parallel do default(shared) private(jptcl) schedule(static) proc_bind(close)
            do jptcl=1,self%N
                if( jptcl /= iptcl )then
                    if( self%labels(jptcl) == new_cls .or. self%labels(jptcl) == cls )then
                        call self%per_ptcl_sim(jptcl)
                    endif
                endif
            end do
            !$omp end parallel do
        endif
    end subroutine move

    !>  \brief  calculates the per-particle similarity
    subroutine per_ptcl_sim( self, iptcl )
        class(shc_cluster), intent(inout) :: self
        integer,            intent(in)    :: iptcl
        integer :: pop, jptcl
        pop = 0
        self%SPS(iptcl) = 0.
        do jptcl=1,self%N
            if( iptcl == jptcl ) cycle
            if( self%labels(iptcl) == self%labels(jptcl) )then
                self%SPS(iptcl) = self%SPS(iptcl)+self%S(iptcl,jptcl)
                pop = pop+1
            endif
        end do
        if( pop > 0 )then
            self%SPS(iptcl) = self%SPS(iptcl)/real(pop)
        else
            self%SPS(iptcl) = self%MINS
        endif
    end subroutine per_ptcl_sim

    !>  \brief  evaluates a solution element replacement
    function eval_change( self, iptcl, cls ) result( sim )
        class(shc_cluster), intent(inout) :: self
        integer,            intent(in)    :: iptcl, cls
        real    :: sim
        integer :: pop, jptcl
        pop = 0
        sim = 0.
        do jptcl=1,self%N
            if( iptcl == jptcl ) cycle
            if( cls == self%labels(jptcl) )then
                sim = sim+self%S(iptcl,jptcl)
                pop = pop+1
            endif
        end do
        if( pop > 0 )then
            sim = sim/real(pop)
        else
            self%SPS(iptcl) = self%MINS
        endif
    end function eval_change

    !>  \brief  is the unit test
    subroutine test_shc_cluster
        real              :: smat(10,10)
        integer           :: labels(10), iptcl, jptcl
        type(shc_cluster) :: shcc
        type(oris)        :: o
        real              :: sim
        o = oris(10)
        ! make known solution
        labels(1:5)  = 1
        labels(6:10) = 2
        do iptcl=1,9
            do jptcl=iptcl+1,10
                if( labels(iptcl) == labels(jptcl) )then
                    smat(iptcl,jptcl) = 1.
                    smat(jptcl,iptcl) = 1.
                else
                    smat(iptcl,jptcl) = 0.5
                    smat(jptcl,iptcl) = 0.5
                endif
            end do
        end do
        call shcc%new(10, 2, smat, o)
        call shcc%shc(.true., 'class', sim)
        call o%write('test_shc_cluster_doc.txt')
    end subroutine test_shc_cluster

    !>  \brief  is the destructor
    subroutine kill( self )
        class(shc_cluster), intent(inout) :: self
        if( self%existence )then
            self%o_ptr => null()
            self%S     => null()
            deallocate( self%labels, self%SPS , stat=alloc_stat)
            if(alloc_stat /= 0) allocchk('In: nsimple_shc_cluster::kill ')
            self%existence = .false.
        endif
    end subroutine kill

end module simple_shc_cluster
