!> Simple ori methods module: scatter orientation search 
module simple_scatter_orisrch
use simple_defs
use simple_oris, only: oris
use simple_ori,  only: ori
implicit none

public :: scatter_orisrch
private

type scatter_orisrch
    private
    integer    :: ndiverse = 50
    integer    :: nquality = 20
    integer    :: nopt     = 0
    type(oris) :: ref_set_diverse
    type(oris) :: ref_set_quality
    logical    :: found_better =.false.
    logical    :: quality_set  =.false.
    logical    :: exists       =.false.
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! GETTERS
    procedure          :: get_nopt
    procedure          :: get_diverse_ref_oris
    procedure          :: get_quality_ref_oris
    procedure          :: get_best_quality_ref_ori
    ! MODIFIERS
    procedure          :: reset
    procedure          :: diversify
    ! GENERATORS
    procedure          :: gen_initial_pool
    procedure          :: update_pool
    procedure          :: gen_ref_sets_from_pool
    procedure, private :: gen_subset
end type scatter_orisrch

real, parameter :: DIV_SCORE_LIMIT = 0.02

contains

    ! CONSTRUCTOR

    subroutine new( self, ndiverse, nquality, oprev, athres )
        use simple_jiffys, only: alloc_err
        class(scatter_orisrch), intent(inout) :: self
        integer,    optional,   intent(in)    :: ndiverse, nquality
        class(ori), optional,   intent(in)    :: oprev
        real,       optional,   intent(in)    :: athres
        if( present(ndiverse) ) self%ndiverse  = ndiverse
        if( present(nquality) ) self%nquality  = nquality
        call self%ref_set_diverse%new(self%ndiverse)
        call self%diversify(oprev, athres)
        call self%ref_set_quality%new(self%nquality)
        self%found_better = .false. 
        self%quality_set  = .false.
        self%exists       = .true.
    end subroutine new

    ! GETTERS

    integer function get_nopt( self )
        class(scatter_orisrch), intent(in) :: self
        get_nopt = self%nopt
    end function get_nopt

    function get_diverse_ref_oris( self ) result( os )
        class(scatter_orisrch), intent(in) :: self
        type(oris) :: os
        os = self%ref_set_diverse
    end function get_diverse_ref_oris

    function get_quality_ref_oris( self ) result( os )
        class(scatter_orisrch), intent(in) :: self
        type(oris) :: os
        os = self%ref_set_quality
    end function get_quality_ref_oris

    function get_best_quality_ref_ori( self ) result( o )
        class(scatter_orisrch), intent(inout) :: self
        type(ori) :: o
        real      :: curr_scores(self%nquality)
        integer   :: i, loc(1)
        ! extract current quality set scores
        do i=1,self%nquality
            curr_scores(i) = self%ref_set_quality%get(i,'corr')
        end do
        loc = maxloc(curr_scores)
        o   = self%ref_set_quality%get_ori(loc(1))
    end function get_best_quality_ref_ori

    ! MODIFIERS

    subroutine reset( self )
        class(scatter_orisrch), intent(inout) :: self
        self%nopt         = 0
        self%found_better = .false. 
        self%quality_set  = .false.
    end subroutine reset

    subroutine diversify( self, oprev, athres )
        class(scatter_orisrch), intent(inout) :: self
        class(ori), optional,   intent(in)    :: oprev
        real,       optional,   intent(in)    :: athres
        type(ori) :: o
        integer   :: nreplace, i
        integer, allocatable :: oriinds(:)
        call self%ref_set_diverse%gen_diverse
        if( present(oprev) )then
            if( present(athres))then
                ! replace 20% of the diverse oris with random samples
                ! around the previous best
                nreplace = nint(0.2*real(self%ndiverse))
                allocate( oriinds(nreplace) )
                call self%ref_set_diverse%find_closest_oris(oprev, oriinds)
                do i=1,nreplace
                    call o%rnd_euler(oprev, athres)
                    call self%ref_set_diverse%set_ori(oriinds(i), o)
                end do
                deallocate(oriinds)
            else
                stop 'need athres dummy arg in conjunction with&
                &oprev; simple_scatter_orisrch :: diversify'
            endif
        endif
    end subroutine diversify

    ! GENERATORS

    function gen_initial_pool( self ) result( pool )
        class(scatter_orisrch), intent(inout) :: self
        type(oris) :: pool, subset_div
        pool       = self%get_diverse_ref_oris()
        subset_div = self%gen_subset('diverse')
        call pool%merge(subset_div)
    end function gen_initial_pool

    function update_pool( self ) result( pool )
        class(scatter_orisrch), intent(inout) :: self
        type(oris) :: pool, subset_div
        pool       = self%gen_subset('quality')
        subset_div = self%gen_subset('diverse')
        call pool%merge(subset_div)
    end function update_pool

    subroutine gen_ref_sets_from_pool( self, pool )
        use simple_math, only: hpsort
        class(scatter_orisrch), intent(inout) :: self
        class(oris),            intent(inout) :: pool
        integer, allocatable :: quality_ranks(:), diversity_ranks(:)
        real,    allocatable :: diversity_scores(:)
        real      :: prev_scores(self%nquality), curr_scores(self%nquality)
        real      :: rswap, div_score
        type(ori) :: o
        integer   :: i, noris_in_pool, cnt, loc_worst_prev(1), loc_best_curr(1)
        quality_ranks = pool%order()
        if( self%quality_set )then
            ! extract previous and current quality set scores
            do i=1,self%nquality
                prev_scores(i) = self%ref_set_quality%get(i,'corr')
                curr_scores(i) = pool%get(quality_ranks(i),'corr')
            end do
            ! update quality set
            cnt = 0
            loc_worst_prev    = minloc(prev_scores)
            loc_best_curr     = maxloc(curr_scores)
            self%found_better = .false.
            do while(curr_scores(loc_best_curr(1)) > prev_scores(loc_worst_prev(1)))
                o = pool%get_ori(quality_ranks(loc_best_curr(1)))
                div_score = self%ref_set_quality%gen_diversity_score(o)
                if( div_score > DIV_SCORE_LIMIT )then
                    call self%ref_set_quality%set_ori(loc_worst_prev(1),o)
                    cnt = cnt + 1
                    self%found_better = .true.
                    self%nopt = self%nopt + 1
                endif
                rswap = prev_scores(loc_worst_prev(1))
                prev_scores(loc_worst_prev(1)) = curr_scores(loc_best_curr(1))
                curr_scores(loc_best_curr(1)) = rswap
                loc_worst_prev = minloc(prev_scores)
                loc_best_curr  = maxloc(curr_scores)
            end do
        else
            ! the nquality best solutions become part of the quality ref_set
            do i=1,self%nquality
                o = pool%get_ori(quality_ranks(i))
                call self%ref_set_quality%set_ori(i,o)
            end do
            self%nopt         = self%nquality
            self%found_better = .true.
            self%quality_set  = .true.
        endif
        ! the ndiverse most diverse solutions become part of the diverse ref_set
        diversity_scores = pool%gen_diversity_scores(self%ref_set_quality)
        noris_in_pool    = pool%get_noris()
        allocate(diversity_ranks(noris_in_pool))
        diversity_ranks = (/(i,i=1,noris_in_pool)/)
        call hpsort(noris_in_pool, diversity_scores, diversity_ranks)
        cnt = 0
        do i=noris_in_pool,1,-1
            o = pool%get_ori(diversity_ranks(i))
            div_score = self%ref_set_diverse%gen_diversity_score(o)
            if( div_score > DIV_SCORE_LIMIT )then
                cnt = cnt + 1
                call self%ref_set_diverse%set_ori(cnt,o)
            endif
            if( cnt == self%ndiverse ) exit
        end do
    end subroutine gen_ref_sets_from_pool

    function gen_subset( self, which ) result( os )
        class(scatter_orisrch), intent(inout) :: self
        character(len=*),       intent(in)    :: which
        type(oris) :: os
        select case(which)
            case('diverse')
                os = self%ref_set_diverse%gen_subset()
            case('quality')
                os = self%ref_set_quality%gen_subset()            
            case DEFAULT
                stop 'unsupported which flag; simple_scatter_orisrch :: gen_subset'
        end select
    end function gen_subset

end module simple_scatter_orisrch
